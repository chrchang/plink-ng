

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
#include <string>
#include <set>
#include <cmath>
#include <algorithm>

#include "options.h"
#include "helper.h"
#include "plink.h"
#include "phase.h"
#include "model.h"
#include "linear.h"
#include "stats.h"

//////////////////////////////////////////
// Helper classes to store and sort proxies

class ProxyResult {
public:

  string name;
  double f;
  double r2;
  double odds;
  double chisq;
  double pvalue;

  ProxyResult(string n, double frq, double r, double o, double c, double p)
    : name(n), f(frq), r2(r), odds(o), chisq(c), pvalue(p) { } 

  bool operator< (const ProxyResult & b) const
    {
      return ( pvalue == b.pvalue ? 
	       name < b.name : 
	       pvalue < b.pvalue );
    }
};


class LDPair {
public:
  int s1;
  int s2;
  double ld;
  bool operator< (const LDPair & b) const
  {
    return  ld > b.ld;
  }


}; 


///////////////////////////////////////////////
// Helper function to find n of m combinations

void combinations_recursive(const vector<int> &elems,
			    unsigned long req_len,
			    vector<unsigned long> &pos,
			    unsigned long depth,
			    unsigned long margin,
			    vector<vector<int> > & collection)
{
  
  if (depth >= req_len) {
    vector<int> t;
    for (unsigned long ii = 0; ii < pos.size(); ++ii)
      t.push_back( elems[pos[ii]] );
    collection.push_back(t);
    return;
  }

  if ((elems.size() - margin) < (req_len - depth))
    return;

  for (unsigned long ii = margin; ii < elems.size(); ++ii) {
    pos[depth] = ii;
    combinations_recursive(elems, req_len, pos, depth + 1, ii + 1, collection);
  }
  return;
}


///////////////////////////////////////////////
// For locus l, main proxy association function

void Plink::performProxyTests(int l)
{

  // Consider a particular SNP, l, and 
  // a) form haplotypes surrounding SNPs
  // b) calculate association with phenotype for these SNPs
  // c) calculate r^2 (haplotypic) with allele and haplotypes
  
  // If this is a rarer SNP, use a slightly broader search strategy
  
  if ( locus[l]->freq < par::proxy_planB_threshold )
    {
      par::proxy_kb = par::proxy_kb_planB;
      par::proxy_window = par::proxy_window_planB;
      par::proxy_snp_filter = par::proxy_snp_filter_planB;
      par::proxy_r2_filter_A = par::proxy_r2_filter_A_planB;
      par::proxy_r2_filter_B = par::proxy_r2_filter_B_planB;
      par::proxy_r2_filter_C = par::proxy_r2_filter_C_planB;
    }
  else
    {
      par::proxy_kb = par::proxy_kb_planA;
      par::proxy_window = par::proxy_window_planA;
      par::proxy_snp_filter = par::proxy_snp_filter_planA;
      par::proxy_r2_filter_A = par::proxy_r2_filter_A_planA;
      par::proxy_r2_filter_B = par::proxy_r2_filter_B_planA;
      par::proxy_r2_filter_C = par::proxy_r2_filter_C_planA;
    }



  // Form haplotypes based on the surrounding SNPs
  // Form phenotype based on the patterns of missingness for test SNP
  // Is there any association?

  bool old_silent = par::silent;

  
  ////////////////////////////////////////////////////////////////
  // If we are in 'impute' mode: this is to evaluate imputation, in a
  // 'leave-one-out' manner; here we assume the reference panel is
  // coded as 'missing phenotype' and the rest of the sample is coded
  // as non-missing phenotype; so we must first blank out (but later
  // replace) any genotype data for these individuals

  vector<bool> tmp1, tmp2;

  if ( par::proxy_impute || par::proxy_leave_out ) 
    {
      
      // Pretend that we do not have these genotypes, except for the
      // reference panel (i.e. individuals with a missing phenotype)

      for ( int i = 0 ; i < n ; i++ ) 
	if ( ! sample[i]->missing )
	  {
	    tmp1.push_back( SNP[l]->one[i] );
	    tmp2.push_back( SNP[l]->two[i] );
	    
	    SNP[l]->one[i] = true;
	    SNP[l]->two[i] = false;
	  }
    }
 

  ///////////////////////
  // Form test haplotypes
  
  CSNP * s = SNP[l];
  
  vector<int> proxyHaplotypePlusSNP;


  // Add reference SNP

  proxyHaplotypePlusSNP.push_back(l);


  // Either a fixed maximum number of SNPs left and right (allowing
  // for different filters, and chromosome ends) or read from a file
  // (these SNPs must be on same chromosome)
  
  if ( par::proxy_list )
    {
      checkFileExists( par::proxy_list_file );
      printLOG("Reading proxy list from [ " 
	       + par::proxy_list_file + " ]\n");
      ifstream PL( par::proxy_list_file.c_str(), ios::in);
      
      map<string,int> mlocus;
      for (int j=0;j<nl_all;j++)
	mlocus.insert(make_pair(locus[j]->name,j));

      while ( ! PL.eof() )
	{
	  string psnp;
	  PL >> psnp;

	  if ( psnp == "" )
	    continue;
	  
	  map<string,int>::iterator m = mlocus.find( psnp );
	  
	  if ( m == mlocus.end() )
	    continue;
	  
	  // Add if okay w/ MAF and genotyping thresholds
	  
	  int pn = m->second;
	  
	  if ( pn != l && 
	       locus[pn]->chr == locus[l]->chr && 
	       locus[pn]->freq >= par::proxy_maf && 
	       abs(double((locus[l]->bp - locus[pn]->bp)/1000.0)) 
	       <= par::proxy_kb &&
	       locus[pn]->pos <= par::proxy_geno )
	    {
	      proxyHaplotypePlusSNP.push_back( pn );	      
	    }
	  
	}

      PL.close();
    }
  
  else // ... use window approach
    {
      
      int i = l-1;
      int added = 0;
      
      while ( added < par::proxy_window )
	{
	  if ( i >= 0 && 
	       locus[i]->chr == locus[l]->chr )
	    {
	      
	      // Add MAF and genotyping thresholds here
	      if ( locus[i]->freq >= par::proxy_maf && 
		   abs(double((locus[l]->bp - locus[i]->bp)/1000.0)) 
		   <= par::proxy_kb &&
		   locus[i]->pos <= par::proxy_geno )
		{
		  proxyHaplotypePlusSNP.push_back(i);
		  added++;
		}
	      
	      // Shift left
	      --i;
	    }
	  else
	    {
	      // Cannot add any more
	      added = par::proxy_window;
	    }
	}
      
      // Now move right
      added = 0;
      i = l+1;
      while ( added < par::proxy_window )
	{
	  if ( i < nl_all && 
	       locus[i]->chr == locus[l]->chr )
	    {
	      
	      // Add MAF, kb and genotyping thresholds here
	      if ( locus[i]->freq >= par::proxy_maf && 
		   abs(double((locus[l]->bp - locus[i]->bp)/1000.0)) 
		   <= par::proxy_kb &&
		   locus[i]->pos <= par::proxy_geno )
		{
		  proxyHaplotypePlusSNP.push_back(i);
		  added++;
		}
	      
	      //Shift right
	      ++i;
	    }
	  else
	    {
	      // Cannot add any more
	      added = par::proxy_window;
	    }
	  
	}
    }
  
  
  ///////////////////////////////////////////////////////////////////////////
  //
  // Optionally, filter list based on LD with reference and with eachother
  //
  ///////////////////////////////////////////////////////////////////////////


  if ( par::proxy_r2_filter ) 
    {
      
      // Only use Reference Panel for these r-sq calculations at this stage; 
      // add flag to modify this

      if ( par::proxy_reference_only )
	haplo->reference_only = true;
      
      set<int> added;

      // Examine SNP number list: proxyHaplotypePlusSNP
      // First entry is always the reference SNP

      set<LDPair> proxies;
      
      // Skip first entry
      for (int i=1; i<proxyHaplotypePlusSNP.size(); i++)
	{
	  LDPair p;
	  p.s1 = l;
	  p.s2 = proxyHaplotypePlusSNP[i];
	  p.ld = haplo->rsq(p.s1,p.s2);
	  //p.ld = haplo->dprime(p.s1,p.s2);

	  proxies.insert(p);
	}
      
      set<LDPair>::iterator pi = proxies.begin();

      while ( pi != proxies.end() )
	{
	  
	  // Enough proxies already?
	  
	  if ( added.size() >= par::proxy_snp_filter )
	    break;
	  
	  if ( ( added.size() < 2 && pi->ld >= par::proxy_r2_filter_A ) // low filter
	       || pi->ld >= par::proxy_r2_filter_B )  // higher filter, once 2 proxies found
	    {
	      
	      // But does this already correlate too strongly with an
	      // existing proxy?
	      
	      set<int>::iterator si = added.begin();
	      
	      bool okay = true;
	      while ( si != added.end() )
		{
		  int2 snps;
		  
		  if ( snps.p1 > snps.p2 )
		    {
		      snps.p1 = *si;
		      snps.p2 = pi->s2;
		    }
		  else
		    {
		      snps.p1 = pi->s2;
		      snps.p2 = *si;
 		    }
		  
		  double ld;
		  
		  map<int2,double>::iterator f = proxyLD.find(snps); 
		  if ( f != proxyLD.end() )
		    {
		      ld = f->second;
		    }
		  else
		    {
		      ld = haplo->rsq(snps.p1,snps.p2);
		      //ld = haplo->dprime(snps.p1,snps.p2);
		      proxyLD.insert(make_pair(snps,ld));
		    }
		  
		  if ( ld > par::proxy_r2_filter_C )
		    {
		      okay = false;
		      break;
		    }
		  ++si;
		}
	      
	      if ( okay ) 
		{
		  added.insert( pi->s2 );
		}
	    }
	  
	  // Consider next proxy SNP
	  ++pi;
	}


      //////////////////////////////////
      // Update the fitlered proxy list
      
      proxyHaplotypePlusSNP.clear();
      proxyHaplotypePlusSNP.push_back(l);
      set<int>::iterator si = added.begin();
      while ( si != added.end() )
	{
	  proxyHaplotypePlusSNP.push_back( *si );
	  ++si;

	}
 

      // And remember to reset
      haplo->reference_only = false;
    }


  //////////////////////////////////////////////////////////////////
  //
  // Sort, and select reference SNP
  //
  //////////////////////////////////////////////////////////////////

  sort( proxyHaplotypePlusSNP.begin(), 
	proxyHaplotypePlusSNP.end());

  int cnt = proxyHaplotypePlusSNP.size();
  int ref;

  for (int i=0; i<cnt; i++)
    if ( proxyHaplotypePlusSNP[i] == l )
      {
	ref = i;
	break;
      }


  //////////////////////////////////////////////////////////////////////
  //
  // Do we actually have 1 or more proxy SNP selected? If not, leave now
  // if in verbose mode
  //
  //////////////////////////////////////////////////////////////////////
  
  if ( par::proxy_all && ( ! par::proxy_full_report ) )
    {
      if ( false && proxyHaplotypePlusSNP.size() == 1 )
 	{
	  
 	  if ( ! par::proxy_impute )
 	    {
 	      haplo->HTEST << setw( 4 ) << locus[l]->chr << " " 
 			   << setw(par::pp_maxsnp) << locus[l]->name << " "
 			   << setw(12) << locus[l]->bp << " "
 			   << setw(4) << locus[l]->allele1 << " "
 			   << setw(4) << locus[l]->allele2 << " "
 			   << setw(10) << 1 - locus[l]->pos << " "
 			   << setw(4) << proxyHaplotypePlusSNP.size() - 1 << " "
 			   << setw(8) << "NA" << " ";
	      
 	      // Display C/C or T/U for case/control and TDT, else F and BETA
	      
 	      if ( par::qt || par::proxy_glm ) 
 		haplo->HTEST << setw(8) << "NA" << " "
 			     << setw(8) << "NA" << " ";	      
 	      else
 		haplo->HTEST << setw(8) << "NA" << " "
 			     << setw(8) << "NA" << " "
 			     << setw(8) << "NA" << " ";
	      
 	      haplo->HTEST << setw(10) << "NA" << " ";
	      
 	      if ( par::proxy_list_proxies )
 		haplo->HTEST << "(NONE)";
	      
 	      haplo->HTEST << "\n";
	      haplo->HTEST.flush();

	      // Finally, replace any genotypes made temporarily missing

	      if ( par::proxy_leave_out )
		{
		  int cnt = 0;
		  for ( int i = 0 ; i < n ; i++ ) 
		    if ( ! sample[i]->missing )
		      {
			SNP[l]->one[i] = tmp1[cnt];
			SNP[l]->two[i] = tmp2[cnt++];
		      }
		  
		  tmp1.clear();
		  tmp2.clear();
		}

	      return;
 	    }

 	}
     }

  


  ///////////////////////////////////////////////////////////////////////////////
  //
  // Phase haplotypes
  //
  ///////////////////////////////////////////////////////////////////////////////  

  haplo->reset();

  haplo->new_pred_locus.resize(1);

  haplo->new_map.resize(1);

  haplo->new_pred_locus[0] = proxyHaplotypePlusSNP;

  haplo->new_map[0] = locus[l];

  par::silent = true;

  haplo->phaseAllHaplotypes(true,*pperm);

  haplo->hname = locus[l]->name;
  
  par::silent = old_silent;


  ///////////////////////////////////////////////////////////////////////////////
  //
  // Per-genotype error screen 
  //
  ///////////////////////////////////////////////////////////////////////////////

  if ( par::proxy_error ) 
    {

      // Identify test SNP with respect to 0..cnt phased region
      // (i.e. 'ref' and not 'l')
      
      haplo->queryGenotype( ref );
      return;
    }
   
 
  //////////////////////////////////////////////////////////////////////////////
  //
  // If we are in 'impute' mode: do not perform association tests, but
  // just compare imputed genotypes to actual (left-out) genotypes,
  // and report on this.
  //
  //////////////////////////////////////////////////////////////////////////////  

  if ( par::proxy_leave_out )
    {
      int cnt = 0;
      for ( int i = 0 ; i < n ; i++ ) 
	if ( ! sample[i]->missing )
	  {
	    SNP[l]->one[i] = tmp1[cnt];
	    SNP[l]->two[i] = tmp2[cnt++];
	  }

      tmp1.clear();
      tmp2.clear();
    }
  else if ( par::proxy_impute ) 
    {
      
      ////////////////////////////////////////
      // Replace genotypes, and record dosage
      
      // For imputation quality score
      boolvec_t m1(cnt,false);
      m1[ref] = true;
      boolvec_t a1(cnt,false);
      map<int,int> tests = haplo->makeTestSet(m1,a1);
      set<int> hs;
      map<int,int>::iterator i1 = tests.begin();
      while ( i1 != tests.end() )
	{
	  if ( i1->second == 0 )
	    hs.insert( i1->first);
	  ++i1;
	}

      haplo->calculateEmpiricalVariance(hs);
    

      int con[4][4];
      for (int j=0;j<4;j++)
	for (int k=0; k<4; k++)
	  con[j][k] = 0; 
      

      if ( par::proxy_record_dosage )
	OUTFILE << locus[l]->name << "\t"
		<< locus[l]->allele1 << "\t"
		<< locus[l]->allele2 << "\t"
		<< haplo->ratio << "\t";

      int cnt = 0;
      for ( int i = 0 ; i < n ; i++ ) 
	if ( ! sample[i]->missing )
	  {
	    
	    // Actual genotypes

	    bool a1 = tmp1[cnt];
	    bool a2 = tmp2[cnt++];

	    int og;

	    if ( a1 )
	      {
		if ( a2 ) 
		  og = 2;
		else
		  og = 3;
	      }
	    else
	      {
		if ( a2 ) 
		  og = 1;
		else
		  og = 0;
	      }
	    
	    // Imputed genotypes

	    vector_t g = haplo->imputeGenotype(i,ref);

	    // Call imputed

	    bool i1, i2;

	    if ( g[0] > par::proxy_impute_threshold ) 
	      {
		i1 = i2 = false;
		con[og][0]++;
	      }
	    else if ( g[1] > par::proxy_impute_threshold ) 
	      { 
		i1 = false; 
		i2 = true; 
		con[og][1]++;
	      }
	    else if ( g[2] > par::proxy_impute_threshold ) 
	      {
		i1 = i2 = true;
		con[og][2]++;
	      }
	    else 
	      { 
		i1 = true; i2 = false; 
		con[og][3]++;
	      }

	    if ( par::proxy_full_report ) 
	      haplo->HTEST << locus[l]->name << "\t" 
			   << sample[i]->fid << " " 
			   << sample[i]->iid << "\t" 
			   << a1<<a2 << " " 
			   << i1<<i2 << " " 
			   << g[0] << " " 
			   << g[1] << " " 
			   << g[2] << "\n";
	    
	    if ( par::proxy_record_dosage )
	      OUTFILE << g[0] + 0.5 * g[1] << "\t";
	    
	    
	    ////////////////////////////////////////////////
	    // Replace missing slots with imputed genotypes
	    
	    // Note: because genotyping rate is pre-calculated, this 
	    // will not be a problem in practice (as long as panel is 
	    // sufficiently large) -- i.e. in smaller cases, we need 
	    // to worry about whether this imputed SNP now might be 
	    // used as a proxy for another SNP -- wouldn't be terrible, 
	    // but probably better to try to use real SNPs as 
	    // proxies whenever possible.

//  	    if ( ! par::proxy_impute_preserve_genotyped ) 
// 	      {
		if ( par::proxy_impute_replace || ( a1 && ! a2 ) ) 
		  {
		    SNP[l]->one[i] = i1;
		    SNP[l]->two[i] = i2;
		  }
		else
		  {
		    SNP[l]->one[i] = a1;
		    SNP[l]->two[i] = a2;		
		  }
//	      }
	  }


      if ( par::proxy_record_dosage )
	OUTFILE << "\n";

      // Report matrix of concordance

      int total = 0;
      for (int j=0;j<4;j++)
	for (int k=0;k<4;k++)
	  total += con[j][k];

      int concordant = con[0][0] + con[1][1] + con[2][2];

      int both_geno = 0;
      for (int j=0;j<3;j++)
	for (int k=0;k<3;k++)
	  both_geno += con[j][k];
      
      int observed_geno = total;
      for (int j=0; j<4; j++)
	observed_geno -= con[3][j];
      
      int imputed_geno = total;
      for (int j=0; j<4; j++)
	imputed_geno -= con[j][3];
      
      double rate = both_geno == 0 ? 
	-1 : 
	(double)concordant/(double)both_geno ;
      
      double rate_obs = (double)observed_geno/(double)total;
      double rate_imp = (double)imputed_geno/(double)total;
      double rate_ovr = (double)both_geno/(double)total;
            
      if ( par::proxy_full_report ) 
	{
	  haplo->HTEST << "\nImputation matrix "
		       << "(rows observed, columns imputed)\n\n";
	  
	  for (int j=0;j<4;j++)
	    {
	      
	      if ( j == 0 ) 
		haplo->HTEST << locus[l]->allele1 << "/"
			     << locus[l]->allele1 << "\t";
	      else if ( j == 1 ) 
		haplo->HTEST << locus[l]->allele1 << "/"
			     << locus[l]->allele2 << "\t";
	      else if ( j == 2 ) 
		haplo->HTEST << locus[l]->allele2 << "/"
			     << locus[l]->allele2 << "\t";
	      else
		haplo->HTEST << par::missing_genotype << "/"
			     << par::missing_genotype << "\t";	  
	      
	      for (int k=0; k<4; k++)
		haplo->HTEST << con[j][k] << "\t";
	      
	      haplo->HTEST << "\n";
	      
	    }
	  
	  haplo->HTEST << setw(4) << "CHR" << " " 
		       << setw(par::pp_maxsnp) << "SNP" << " "
		       << setw(8) << "INFO" << " " 
		       << setw(4) << "NPRX" << " "
		       << setw(8) << "TOTAL_N" << " "
		       << setw(8) << "OBSERVD" << " "
		       << setw(8) << "IMPUTED" << " "
		       << setw(8) << "OVERLAP" << " "
		       << setw(8) << "CONCORD" << " ";
	  
	  if ( par::proxy_impute_genotypic_concordance )
	    haplo->HTEST << setw(8) << "F_AA" << " "
			 << setw(8) << "I_AA" << " "
			 << setw(8) << "C_AA" << " "
			 << setw(8) << "F_AB" << " "
			 << setw(8) << "I_AB" << " "
			 << setw(8) << "C_AB" << " "
			 << setw(8) << "F_BB" << " "
			 << setw(8) << "I_BB" << " "
			 << setw(8) << "C_BB" << " ";
	  
	  if ( par::proxy_list_proxies )
	    haplo->HTEST << "SNPS";

	  haplo->HTEST << "\n";
	  
	}

            
      haplo->HTEST << setw(4) << locus[l]->chr << " " 
		   << setw(par::pp_maxsnp) << locus[l]->name << " "
		   << setw(4) << proxyHaplotypePlusSNP.size() - 1 << " "
		   << setw(8) << haplo->ratio << " " 
		   << setw(8) << total << " "
		   << setw(8) << rate_obs << " "
		   << setw(8) << rate_imp << " "
		   << setw(8) << rate_ovr << " ";

      if ( rate >= 0 )
	haplo->HTEST << setw(8) <<  rate << " ";
      else
	haplo->HTEST << setw(8) << "NA" << " ";

      if ( par::proxy_impute_genotypic_concordance )
	{

	  int impAA = con[0][0]+con[0][1]+con[0][2];
	  int impAB = con[1][0]+con[1][1]+con[1][2];
	  int impBB = con[2][0]+con[2][1]+con[2][2];
	  
	  double fAA = (double)(con[0][0]+con[1][0]+con[2][0]) / (double)observed_geno;
	  double fAB = (double)(con[0][1]+con[1][1]+con[2][1]) / (double)observed_geno;
	  double fBB = (double)(con[0][2]+con[1][2]+con[2][2]) / (double)observed_geno;
	  
	  if ( observed_geno == 0 ) 
	    fAA = fAB = fBB = 0;
	  
	  int obsAA = impAA+con[0][3];
	  int obsAB = impAB+con[1][3];
	  int obsBB = impBB+con[2][3];

	  haplo->HTEST << setw(8) << fAA << " ";

	  if ( obsAA == 0 ) 
	    haplo->HTEST << setw(8) << "NA" << " ";
	  else
	    haplo->HTEST << setw(8) << (double)(impAA)/(double)(impAA+con[0][3]) << " ";	  

	  if ( impAA==0 )
	    haplo->HTEST << setw(8) << "NA" << " ";
	  else
	    haplo->HTEST << setw(8) << (double)(con[0][0])/(double)impAA << " ";
	  
	  haplo->HTEST << setw(8) << fAB << " ";
	  if ( obsAB == 0 ) 
	    haplo->HTEST << setw(8) << "NA" << " ";
	  else
	  haplo->HTEST << setw(8) << (double)(impAB)/(double)(impAB+con[1][3]) << " ";
	  if ( impAB==0 )
	    haplo->HTEST << setw(8) << "NA" << " ";
	  else
	    haplo->HTEST << setw(8) << (double)(con[1][1])/(double)impAB << " ";

	  haplo->HTEST << setw(8) << fBB << " ";
	  if ( obsBB == 0 ) 
	    haplo->HTEST << setw(8) << "NA" << " ";
	  else
	  haplo->HTEST << setw(8) << (double)(impBB)/(double)(impBB+con[2][3]) << " ";
	  if ( impBB==0 )
	    haplo->HTEST << setw(8) << "NA" << " ";
	  else
	    haplo->HTEST << setw(8) << (double)(con[2][2])/(double)impBB << " ";
	  
	}

      
      if ( par::proxy_list_proxies )
	{
	  bool printed = false;
	  for (int l0=0; l0< proxyHaplotypePlusSNP.size(); l0++)
	    {
	      if ( proxyHaplotypePlusSNP[ l0 ] != l ) 
		{
		  if ( printed ) 
		    haplo->HTEST << "|";		  
		  haplo->HTEST << locus[ proxyHaplotypePlusSNP[ l0 ] ]->name;
		  printed = true;
		  
		}
	    }
	}
      
      
      haplo->HTEST << endl;
      haplo->HTEST.flush();

      tmp1.clear();
      tmp2.clear();
      
      return;
    }




  ////////////////////////////////////////////////////////////////
  //                                                            //
  // Consider all subsets of subhaplotypes, if in verbose mode  //
  //                                                            //
  ////////////////////////////////////////////////////////////////

  
  if ( ! ( par::proxy_all || par::proxy_full_report ) )
    printLOG("Estimated haplotype frequencies: now considering combinations...\n");




  ////////////////////////////////////////////////////
  // Display haplotype frequencies and individual r^2

  if ( ( ! par::proxy_all ) || par::proxy_full_report )
    haplo->HTEST << "\n"
		 << "     *** Proxy haplotype association report for " 
		 << haplo->hname << " *** \n\n";
  

  ////////////////////////////////////////
  // Report SNPs and single-SNP r-squared
	    
  if ( ( ! par::proxy_all) || par::proxy_full_report )
    {
      haplo->HTEST << setw(par::pp_maxsnp) << "SNP" << " " 
		   << setw(8) << "MAF" << " " 
		   << setw(8) << "GENO" << " "
		   << setw(8) << "KB" << " "
		   << setw(8) << "RSQ" << " ";
      if ( par::qt ) 
	haplo->HTEST << setw(8) << "BETA" << " "
		     << setw(8) << "STAT" << " "	  
		     << setw(8) << "P" << "\n";
      else 
	haplo->HTEST << setw(8) << "OR" << " "
		     << setw(8) << "CHISQ" << " "	  
		     << setw(8) << "P" << "\n";


      for ( int s = 0 ; s < cnt ; s++)
	{	  
	  haplo->HTEST << setw(par::pp_maxsnp) 
		       << locus[ haplo->new_pred_locus[0][s] ]->name << " "
		       << setw(8) 
		       << locus[ haplo->new_pred_locus[0][s] ]->freq << " "
		       << setw(8) 
		       << 1- locus[ haplo->new_pred_locus[0][s] ]->pos  << " "
		       << setw(8) 
		       << (double)(locus[ haplo->new_pred_locus[0][s]]->bp - 
				   locus[ haplo->new_pred_locus[0][ref]]->bp)/1000<<" ";
	  
	  // R-squared

	  if ( s == ref ) 
	    haplo->HTEST << setw(8) << "*" << " ";
	  else
	    haplo->HTEST << setw(8) << haplo->rsq_internal(s,ref) << " ";

	  // Single SNP association

	  boolvec_t snpmask(cnt,false);
	  boolvec_t dummy_allele(cnt,true);
	  snpmask[s] = true;
	  
	  haplo->testSet = haplo->makeTestSet(snpmask,dummy_allele);

	  if (par::proxy_CC)
	    {
	      // GLM or standard test?
	      if ( par::proxy_glm )
		{
		  glmAssoc(false,*pperm);
		}
	      else
		{		   
		  if ( par::qt )
		    haplo->haplotypicQTL(haplo->testSet,2,false);
		  else
		    haplo->haplotypicCC(haplo->testSet,2,false);
		}
	    }
	  else if (par::proxy_TDT)
	    {
	      
	      haplo->subhaplotypes = true;
	      haplo->downcoding = haplo->testSet;
	      
	      haplo->trans.clear();
	      haplo->untrans.clear();
	      
	      haplo->trans.resize(2,0);
	      haplo->untrans.resize(2,0);
	      haplo->nt = 2;

	      // First rescore transmissions
	      
	      for (int i=0; i<n; i++)
		{
		  if ( (!sample[i]->founder) && 
		       haplo->include[i] )
		    {		      
		      haplo->transmissionCount(i,haplo->phasemap[i]);
		    }
		}
	      

	      // Then recount T:U based on revised transmission counts
	      
	      haplo->haplotypicTDT(haplo->testSet,2,false);
	      
	      haplo->subhaplotypes = false;
	      haplo->downcoding.clear();
	      
	    }
	  
	  // Recover main results from Model, if GLM used
	  if ( par::proxy_glm )
	    {
	      vector_t coef = model->getCoefs();
	      haplo->odds = par::bt ? exp(coef[1]) : coef[1];
	      haplo->result = model->isValid() ? model->getStatistic() : 0;
	      haplo->pvalue = par::bt ? chiprobP(haplo->result,1) : ((LinearModel*)model)->getPValue();
	      delete model;
	    }

	  haplo->HTEST << setw(8) << haplo->odds << " "
		       << setw(8) << haplo->result << " "
		       << setw(8) << haplo->pvalue << "\n";
	  
	}

      haplo->HTEST << "\n\n";
    }

  
  ///////////////////////////////////////////////////
  // Report haplotypes and frequencies, and also all
  // haplotype-specific tests and omnibus test result
  
  if ( ( ! par::proxy_all ) || par::proxy_full_report )
    {
      string str = "";
      for ( int i = 0; i < cnt ; i++)
	if ( i == ref ) 
	  str += "*";
	else
	  str += ".";
      
      haplo->HTEST << setw(14) << str << " "
		   << setw(10) << "FREQ" << " ";
      if ( par::qt ) 
	haplo->HTEST << setw(8) << "BETA" << " "
		     << setw(8) << "STAT" << " "	  
		     << setw(8) << "P" << "\n";
      else 
	haplo->HTEST << setw(10) << "OR" << " "
		     << setw(10) << "CHISQ" << " " 
		     << setw(10) << "P" << "\n";

      for (int h=0; h< haplo->nh; h++)      
	{
	  if (haplo->f[h] >= par::proxy_mhf )
	    {	  
	      haplo->HTEST << setw(14) << haplo->haplotypeName(h) << " " 
			   << setw(10) << haplo->f[h] << " ";
	      

	      haplo->testSet.clear();
	      
	      for (int h2=0; h2 < haplo->nh; h2++)
		{
		  if ( haplo->f[h2] >= par::proxy_mhf)
		    {
		      if (h==h2) 
			{
			  haplo->testSet.insert(make_pair(h2,0));
			}
		      else 
			{
			  haplo->testSet.insert(make_pair(h2,1));
			}
		    }
		}

	      if (par::proxy_CC)      
		{
		  // GLM or standard test?
		  if ( par::proxy_glm )
		    {
		      glmAssoc(false,*pperm);
		    }
		  else
		    {		   
		      
		      if ( par::qt)
			haplo->haplotypicQTL(haplo->testSet,2,false);		
		      else
			haplo->haplotypicCC(haplo->testSet,2,false);		
		    }
		}
	      else if (par::proxy_TDT)
		{

		  haplo->subhaplotypes = true;
		  haplo->downcoding = haplo->testSet;
		  
		  haplo->trans.clear();
		  haplo->untrans.clear();
		  
		  haplo->trans.resize(2,0);
		  haplo->untrans.resize(2,0);
		  haplo->nt = 2;

		  // First rescore transmissions
		  
		  for (int i=0; i<n; i++)
		    {
		      if ( (!sample[i]->founder) && 
			   haplo->include[i] )
			{		      
			  haplo->transmissionCount(i,haplo->phasemap[i]);
			}
		    }
		  
		  // Then recount T:U based on revised transmission counts
		  
		  haplo->haplotypicTDT(haplo->testSet,2,false);
		  
		  haplo->subhaplotypes = false;
		  haplo->downcoding.clear();
		  
		  
		}
	      
	      // Recover main results from Model, if GLM used
	      if ( par::proxy_glm )
		{
		  vector_t coef = model->getCoefs();
		  haplo->odds = par::bt ? exp(coef[1]) : coef[1];
		  haplo->result = model->isValid() ? model->getStatistic() : 0;
		  haplo->pvalue = par::bt ? chiprobP(haplo->result,1) : ((LinearModel*)model)->getPValue();
		  delete model;
		}	      
	      
	      haplo->HTEST << setw(10) << haplo->odds << " "
			   << setw(10) << haplo->result << " "
			   << setw(10) << haplo->pvalue << "\n";
	      haplo->HTEST.flush();

	    }
	}


      haplo->HTEST << "\nHaplotype frequency estimation based on " 
		   << haplo->validN; 

      if ( haplo->X ) 
	{
	  int found_chr = 0;
	  for (int i=0; i<n; i++)
	    if ( sample[i]->founder )
	      {
		if ( sample[i]->sex )
		  found_chr++;
		else
		  found_chr+=2;
	      }
	  haplo->HTEST << " of " << found_chr << " founder chromosomes\n";
	}
      else if ( haplo->haploid )
	haplo->HTEST << " of " << haplo->cnt_f << " founder chromosomes\n";
      else
	haplo->HTEST << " of " << haplo->cnt_f * 2 << " founder chromosomes\n";

		   

      ///////////////////////////////////
      // Omnibus test: C/C only

      if ( par::proxy_CC && ! par::qt )
	{
	   
	  map<int,int> tests;
	  int nch=0;
	  for (int h=0; h < haplo->nh; h++)
	    if ( haplo->f[h] >= par::proxy_mhf)
	      tests.insert(make_pair(h,nch++));
	  
	  if (nch>2)
	    {
	      haplo->haplotypicCC(tests,nch,false);
	      
	      haplo->HTEST << "Omnibus haplotype test statistic: " 
			   << haplo->result << ", df = " << nch-1 << ", "
			   << "p = " 
			   << chiprobP( haplo->result , nch-1 ) << "\n\n";
	    }
	  
	} 
    }



  ///////////////////////////////////////////
  // Create masks

  // Reference SNP (fix here)
  boolvec_t m1(cnt,false);
  m1[ref] = true;
  boolvec_t a1(cnt,true);

  // Proxy haplotype (populated below)
  boolvec_t m2(cnt,false);
  boolvec_t a2(cnt); 
  

  ///////////////////////////////////////////
  //         Report actual SNP             //
  ///////////////////////////////////////////

  haplo->testSet = haplo->makeTestSet(m1,a1);
  
  if (par::proxy_CC)
    {
      // GLM or standard test?
      if ( par::proxy_glm )
	{
	  glmAssoc(false,*pperm);
	}
      else
	{		   	  
	  if ( par::qt )
	    haplo->haplotypicQTL(haplo->testSet,2,false);
	  else
	    haplo->haplotypicCC(haplo->testSet,2,false);
	}
    }
  else if (par::proxy_TDT)
    {      

      haplo->subhaplotypes = true;
      haplo->downcoding = haplo->testSet;
      
      haplo->trans.clear();
      haplo->untrans.clear();
      
      haplo->trans.resize(2,0);
      haplo->untrans.resize(2,0);
      haplo->nt = 2;

      // First rescore transmissions
      
      for (int i=0; i<n; i++)
	{
	  if ( (!sample[i]->founder) && 
	       haplo->include[i] )
	    {		      
	      haplo->transmissionCount(i,haplo->phasemap[i]);
	    }
	}

    
      // Then recount T:U based on revised transmission counts
      
      haplo->haplotypicTDT(haplo->testSet,2,false);
      
      haplo->subhaplotypes = false;
      haplo->downcoding.clear();
     
      // Also, calculate the information score
      set<int> t1 = haplo->makeSetFromMap(haplo->testSet);
      haplo->calculateEmpiricalVariance(t1);

    }
  
  
  // Recover main results from Model, if GLM used
  if ( par::proxy_glm )
    {
      vector_t coef = model->getCoefs();
      haplo->odds = par::bt ? exp(coef[1]) : coef[1];
      haplo->result = model->isValid() ? model->getStatistic() : 0;
      haplo->pvalue = par::bt ? chiprobP(haplo->result,1) : ((LinearModel*)model)->getPValue();
      delete model;
      
      // Also, calculate the information score
      set<int> t1 = haplo->makeSetFromMap(haplo->testSet);
      haplo->calculateEmpiricalVariance(t1);
    }
  

  //////////////////////////////////////////////////////////////////////
  //
  // Just report the single SNP result (based on haplotype test)?
  //
  //////////////////////////////////////////////////////////////////////

  if ( par::proxy_all && ( ! par::proxy_full_report ) )
    {

      
      
      haplo->HTEST << setw( 4 ) << locus[l]->chr << " " 
		   << setw(par::pp_maxsnp) << locus[l]->name << " "
		   << setw(12) << locus[l]->bp << " "
		   << setw(4) << locus[l]->allele1 << " "
		   << setw(4) << locus[l]->allele2 << " "
		   << setw(10) << 1 - locus[l]->pos << " " 
		   << setw(4) << proxyHaplotypePlusSNP.size() - 1 << " "
		   << setw(8) << haplo->ratio << " ";

      // Display C/C or T/U for case/control and TDT, else F and BETA

      if ( haplo->pvalue < -1 ) 
	{
	  // Not a valid test, e.g. monomorphic, and so p value is returned as -9
	  if ( par::qt || par::proxy_glm ) 
	    haplo->HTEST << setw(8) << "NA" << " "
			 << setw(8) << "NA" << " "
			 << setw(8) << "NA" << " ";
	  else
	    haplo->HTEST << setw(8) << "NA" << " "
			 << setw(8) << "NA" << " "
			 << setw(8) << "NA" << " "
			 << setw(8) << "NA" << " ";	    	  
	}
      else
	{
      
	  if ( par::qt || par::proxy_glm ) 
	    {
	      // Get population haplotype frequency for imputed SNP
	      double f = haplo->freq(m1,a1);
	      
	      haplo->HTEST << setw(8) << f << " "
			   << setw(8) << haplo->odds << " ";
	    }
	  else
	    haplo->HTEST << setw(8) << haplo->case_freq << " "
			 << setw(8) << haplo->control_freq << " "
			 << setw(8) << haplo->odds << " ";
	  
	  haplo->HTEST << setw(10) << haplo->pvalue << " ";
	}
      
      
      if ( par::proxy_list_proxies )
	{
	  bool printed = false;
	  for (int l0=0; l0< proxyHaplotypePlusSNP.size(); l0++)
	    {
	      if ( proxyHaplotypePlusSNP[ l0 ] != l ) 
		{
		  if ( printed )
		    haplo->HTEST << "|";				  
		  haplo->HTEST << locus[ proxyHaplotypePlusSNP[ l0 ] ]->name;
		  printed = true;
 		}
	    }
	}


      haplo->HTEST << endl;
      haplo->HTEST.flush();

      return;
    }


  



  ////////////////////////////////////////////////////////////////////////
  //
  // Rest of this function is for the extended report mode  
  //
  ////////////////////////////////////////////////////////////////////////        



  ///////////////////////////////////////////
  // Consider from cnt-1 to 1 SNP haplotypes

  int search_cnt = par::proxy_include_reference ? cnt : cnt - 1;
  int maxsnp = search_cnt < par::proxy_maxhap ? search_cnt : par::proxy_maxhap;
  
  int num_subhaps_total = 0;
  int num_subhaps_valid = 0;

  set<ProxyResult> presults;
  
  for (int i=1; i<=maxsnp; i++)
    {

      // For an i-SNP of cnt-SNP, consider all the permutations

      vector<unsigned long> pos1(i);
      vector<int> d(search_cnt);
      int i2=0;
      for (int z=0; z<cnt; z++) 
	if ( z != ref ) 
	  {
	    d[i2]=i2;
	    i2++;
	  }

      vector<vector<int> > collection;
      combinations_recursive(d,i,pos1,0,0,collection);

      vector<int> mapback;
      for (int s=0;s<cnt;s++)
	if ( par::proxy_include_reference || s != ref )
	  mapback.push_back(s);
      
      // Consider each combination of i of cnt-1 SNP haplotypes
      // then consider each possible haplotype within that
      
      for (int c1 = 0; c1 < collection.size() ; c1++ )
	{
	  
	  // For biallelic haplotypes (SNPs)
	  // skip the second allele
	  
	  // Which SNPs do we want to look up?
	  
	  vector<int> posit = collection[c1];
	  
	  // Now consider all search_cnt SNP haplotypes 
	  // i.e. excluding reference SNP
	  
	  int hapcnt = (int)pow((double)2,i);
	  int h = 0;
	  while ( h < hapcnt )
	    {
	      
	      // Skip redundant second allele of 
	      // SNPs

	      if ( collection[c1].size() == 1 && 
		   h == 1 )
		{
		  h++;
		  continue;
		}
	      
	      
	      vector<bool> tmp;
	      unsigned int p=1;
	      for (int s=0;s<i;s++)
		{
		  if ( h & p ) tmp.push_back(false);
		  else tmp.push_back(true);
		  p <<= 1;
		}
	      
	      
	      ///////////////////////
	      // Set proxy haplotype
	      
	      // Insert i-SNP haplotype back into 
	      // full cnt-SNP space
	      
	      fill(m2.begin(), m2.end(), false);
	      
	      int t1=0;
	      for (int t2=0; t2<i; t2++)
		{
		  m2[ mapback[ posit[t2] ] ] = true;
		  a2[ mapback[ posit[t2] ] ] = tmp[ t2 ];
		}
	      
	      
	      ////////////////////////////////
	      // Check frequency of haplotype
	      
	      double f = haplo->freq(m2,a2);
	      
	      
	      ////////////////////////////////////////////////////////////
	      // Calculate r^2 between haplotype and reference SNP allele
	      
	      double r2 = haplo->rsq_internal(m1,a1,m2,a2);
	     

	      	      
	      ////////////////////////////////////////////
	      // Is this haplotype not worth considering?
	     

	      ++num_subhaps_total;
 
 	      if ( r2 < par::proxy_r2 ||
		   f < par::proxy_mhf )
 		{		  
 		  h++;
 		  continue;
		}


	      ++num_subhaps_valid;
 

	      ////////////////////////////////////////////////////////////
	      // Calculte association between proxy haplotype and disease
	      
	      // If only two haplotypes, report only 1
	      // (note may only be two *common* haplotypes, but
	      // in that case, we should report both
	      
	      haplo->testSet = haplo->makeTestSet(m2,a2);
		      
	      if (par::proxy_CC)
		{
		  // GLM or standard test?
		  if ( par::proxy_glm )
		    {
		      glmAssoc(false,*pperm);
		    }
		  else
		    {		   
		      if ( par::qt )
			haplo->haplotypicQTL(haplo->testSet,2,false);
		      else
			haplo->haplotypicCC(haplo->testSet,2,false);
		    }
		}
	      else if (par::proxy_TDT)
 		{

		  haplo->subhaplotypes = true;
		  haplo->downcoding = haplo->testSet;
		  
		  haplo->trans.clear();
		  haplo->untrans.clear();
		  
		  haplo->trans.resize(2,0);
		  haplo->untrans.resize(2,0);
		  haplo->nt = 2;
		  
		  // First rescore transmissions
		  
		  for (int i=0; i<n; i++)
		    {
		      if ( (!sample[i]->founder) && 
			   haplo->include[i] )
			{		      
			  haplo->transmissionCount(i,haplo->phasemap[i]);
			}
		    }
		  
		  // Then recount T:U based on revised transmission counts
		  
		  haplo->haplotypicTDT(haplo->testSet,2,false);
		  
		  haplo->subhaplotypes = false;
		  haplo->downcoding.clear();
		}
	      
	      string str = par::proxy_include_reference ? 
		haplo->getSubHaplotypeName(m2,a2,-1) :
		haplo->getSubHaplotypeName(m2,a2,ref);
	      

	      // Recover main results from Model, if GLM used
	      if ( par::proxy_glm )
		{
		  vector_t coef = model->getCoefs();
		  haplo->odds = par::bt ? exp(coef[1]) : coef[1];
		  haplo->result = model->isValid() ? model->getStatistic() : 0;
		  haplo->pvalue = par::bt ? chiprobP(haplo->result,1) : ((LinearModel*)model)->getPValue();
		  delete model;
		}
	      
	      
	      //////////////////////
	      // Store this result

	      ProxyResult r(str,f,r2,
			    haplo->odds,
			    haplo->result,
			    haplo->pvalue);
	      presults.insert(r);
	      
	      // Consider next haplotype
	      
	      h++;
	      
	    }
	}      
    }
  

  haplo->HTEST << "Of " << num_subhaps_total << " subhaplotypes considered, " 
	       << num_subhaps_valid << " met proxy criteria\n\n";
  
  // Report results

  if ( presults.size() == 0 )
    haplo->HTEST << "No proxies found above r-sq " << par::proxy_r2 << "\n";
  else
    {
      
      haplo->HTEST << setw(14) << "HAP" << " "
		   << setw(10) << "FREQ" << " " 
		   << setw(10) << "RSQ" << " ";
      if ( par::qt ) 
	haplo->HTEST << setw(8) << "BETA" << " "
		     << setw(8) << "STAT" << " "	  
		     << setw(8) << "P" << "\n";
      else 
	haplo->HTEST << setw(10) << "OR" << " "
		     << setw(10) << "CHISQ" << " "
		     << setw(10) << "P" << "\n";
      
      
      set<ProxyResult>::iterator i = presults.begin();
      
      while ( i != presults.end() )
	{
	  haplo->HTEST << setw(14) << i->name << " "
		       << setw(10) << i->f << " " 
		       << setw(10) << i->r2 << " " 
		       << setw(10) << i->odds << " " 
		       << setw(10) << i->chisq << " " 
		       << setw(10) << i->pvalue << "\n";
	  i++;
	}
    }

  if ( par::proxy_full_report ) 
    haplo->HTEST << "\n+--------------------------------------------------------------------+\n\n";

  return;
    
}


