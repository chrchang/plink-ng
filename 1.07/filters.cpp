

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
#include "crandom.h"


#define MISSING(i,l) ( SNP[l]->one[i] && ( ! SNP[l]->two[i] ) ) 

void Plink::filterSNPs()
{


  //////////////////////////////////////////////////////
  // This functions applies the following filters and 
  // functions:
  
  // Counts of number of founders, nonfounders 
  // Per-individual genotyping rate
  // Read in, or calculate, allele frequencies (and save these)
  // Optionally write allele frequencies, then close
  // Exclude SNPs with too many missing genotypes
  // Identify/correct heterozygote haploid
  // Identify SNPs with no founder genotypes
  // Calculate/report/filter on HWE tests
  // Calculate/report genotyping rate per SNP/per individual
  // Filter on MAF 
  // Remove filtered-out SNPs

  
  bool original_SNP_major = par::SNP_major;
  
  if ( ! par::SNP_major )
    Ind2SNP();



  // Which SNPs to delete
  vector<bool> del(locus.size(),false);

  // Which individuals to delete
  vector<bool> indel(sample.size(),false);
  

  //////////////////////////////////////////
  // Display number of founders/nonfounders 

  cnt_f=0;
  vector<Individual*>::iterator person = sample.begin();
  while ( person != sample.end() ) 
    {  
      if ( (*person)->founder ) cnt_f++;	    
      person++;
    }
  printLOG(int2str(cnt_f)+" founders and "+int2str(n-cnt_f)+
	   " non-founders found\n");

  if (cnt_f<n) par::has_nonfounders = true;


  ////////////////////////////////////////////
  // If we have obligatory missing genotypes:
  // ensure they really are missng
  
  if ( par::oblig_missing ) 
    {
      set<int2>::iterator p = oblig_missing.begin();
      while ( p != oblig_missing.end() )
	{
	  int l = p->p1;
	  int k = p->p2;

	  for (int i=0; i<sample.size(); i++)
	    {
	      Individual * person = sample[i];
	      if ( person->sol == k )
		{
		  SNP[l]->one[i] = true;
		  SNP[l]->two[i] = false;
		}
	    }
	  ++p;
	}
    }


  /////////////////////////////////////////////////
  // Remove individuals with too many missing calls
  
  double total_genotyping = 0;

  if ( par::MAX_IND_MISSING < 1 )
    {
  
      int n_removed = 0;
      int n_orig = n;

      // Consider each individual
      if ( ! par::oblig_missing )
	{
	  for (int i=0;i<sample.size();i++)
	    {

	      bool female = !sample[i]->sex;

	      // Sum missingness over all SNPs	      
	      int m=0;       // Missing SNPs
	      int nsnps=0;   // All non-obligatory missing SNPs
	      
	      for (int l=0; l<locus.size();l++)
		{
		  
		  // Skip female Y chromosomes
		  if ( female && par::chr_Y[locus[l]->chr] )
		    {
		      continue;
		    }
		  
		  ++nsnps;
		  if ( MISSING(i,l) ) m++;
		}
	      
	      // Too much missingness?
	      if ( (double)m/(double)nsnps > par::MAX_IND_MISSING )
		{
		  indel[i] = true;
		  n_removed++;
		}  	      
	    } // next individual	  
	}
      else // ... allow oblig missing values
	{
	  for (int i=0;i<sample.size();i++)
	    {	      
	      bool female = ! sample[i]->sex;

	      // Sum missingness over all SNPs
	      int m=0;       // Missing SNPs
	      int nsnps=0;   // All non-obligatory missing SNPs
	      
	      for (int l=0; l<locus.size();l++)
		{

		  // Skip female Y chromosomes
		  if ( female && par::chr_Y[locus[l]->chr] )
		    continue;

		  if ( ! obligMissing(i,l) )
		    {
		      if ( MISSING(i,l) )
			++m;
		      ++nsnps;
		    }
		}
	      
	      // Too much missingness?
	      if ( (double)m/(double)nsnps > par::MAX_IND_MISSING )
		{
		  indel[i] = true;
		  n_removed++;
		}
	      
	    } // next individual
	} // end if oblig-missing section
      

      ////////////////////////////////////////
      // Save list of any removed individuals
      
      if (n_removed>0)
	{
	  string f = par::output_file_name + ".irem";
	  printLOG("Writing list of removed individuals to [ " + f + " ]\n");	
	  ofstream REM;
	  REM.open(f.c_str(), ifstream::out);
	  for (int i=0;i<sample.size();i++)
	    if (indel[i])
	      REM << sample[i]->fid << "\t" << sample[i]->iid << "\n";
	  REM.close();      
	  
         
	  // And now remove these individuals, so that 
	  // SNP-based statistics are calculated with 
	  // these samples already excluded
	  
	  n_removed = deleteIndividuals(indel);
	  
	}
      
      printLOG(int2str(n_removed)+" of "+int2str(n_orig));
      printLOG(" individuals removed for low genotyping ( MIND > ");
      printLOG(dbl2str(par::MAX_IND_MISSING)+" )\n");
            
      
    } // end of remove people conditional


  /////////////////////////////////
  // Calculate or read from file? 

  if (par::af_read)
    {
      checkFileExists(par::af_file);
      printLOG( "Reading allele frequencies from [ " + par::af_file + " ] \n");
      
      // Make hash of original SNP names
      map<string,int> mlocus;
      map<string,int>::iterator ilocus;
      
      vector<Locus*>::iterator loc = locus.begin();
      int l=0;
      while ( loc != locus.end() ) 
	{ 
	  mlocus.insert(make_pair( (*loc)->name,l));
	  loc++;
	  l++;
	}


      // Read allele frequencies
      ifstream FRQ;
      FRQ.open(par::af_file.c_str());
      FRQ.clear();
      
      string dum1, dum2, dum3, dum4, dum5, dum6;
      string snpname;
      double freq;
      int nm;
      
      loc = locus.begin();
      while ( loc != locus.end() ) 
	{
	  (*loc)->freq = -1;
	  (*loc)->nm = 0;
	  loc++;
	}
      
      // Skip header line
      FRQ >> dum1 >> dum2 >> dum3 >> dum4 >> dum5 >> dum6;

      while(!FRQ.eof())
	{

	  vector<string> tokens = tokenizeLine(FRQ);
	  
	  if (tokens.size() == 0)
	    continue;
	  else if (tokens.size() != 6)
	    {
	      string sline="";
	      for (int i=0; i<tokens.size(); i++)
		sline += tokens[i] + " ";
	      error("Problem with allele frequency line: 6 fields required:\n"
		    +sline+"\n");
	    }

	  ilocus = mlocus.find(tokens[1]);
	  if (ilocus != mlocus.end())
	    {
	      string rareAllele = tokens[2];
	      string commonAllele = tokens[3];
	      Locus * loc = locus[ilocus->second];
	      
	      if( ! from_string<double>( loc->freq, tokens[4],std::dec))
		{
		  loc->freq = 0;
		  loc->nm = 0;
		}
	      else if( ! from_string<int>(loc->nm,tokens[5],std::dec))
		{
		  loc->freq = 0;
		  loc->nm = 0;
		}

	      // But was that pointing to the correct allele?	      
	      
	      if ( rareAllele == loc->allele2 && 
		   rareAllele != par::missing_genotype &&
		   loc->allele2 != par::missing_genotype )
		loc->freq = 1 - loc->freq;
	      else if ( commonAllele == loc->allele1 && 
			commonAllele != par::missing_genotype &&
			loc->allele1 != par::missing_genotype )
		loc->freq = 1 - loc->freq;

	    }
	}
      FRQ.clear();
      FRQ.close();
      
    }
  

  /////////////////////////////////
  // Calculate allele frequencies
  
  vector<string> hetlist(0);
  
  vector<bool>::iterator d = del.begin();
  vector<Locus*>::iterator loc = locus.begin();
  vector<CSNP*>::iterator s = SNP.begin();
  int l = 0; // Main locus counter

  int exc_maf = 0;
  int exc_miss = 0;
  
  vector<Locus*> no_founders_found_list; 
  
  while ( loc != locus.end() ) 
    {
      
      if (!par::af_read)      
	{
	  (*loc)->freq = 0;
	  // count 1 per allele, for frequency
	  (*loc)->nm = 0; 
	}
      
      // count 1 per genotype, for missingness
      int geno_nm = 0; 
      
      // count 1 per non-obligatory missing genotype 
      // (or set to N individuals)
      int geno_real = 0;
      
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
      int i = 0;

      while ( person != sample.end() ) 
	{
	  
	  bool s1 = *i1;
	  bool s2 = *i2;

	  // Check female Y genotypes
	  if ( par::chr_Y[(*loc)->chr] && ! (*person)->sex )
	    {
	      // Set to missing, unless in a RECODE mode
	      if ( ! par::preserve_all_genotypes )
		{
		  s1 = *i1 = true;
		  s2 = *i2 = false;
		}		  

	      // But in any case, do not include this marker in 
	      // any genotype counts: skip to next person
	      
	      ++person;
	      ++i;
	      ++i1;
	      ++i2;
	      
	      continue;
	    }
	  
	  // For haploid heterozygosity check, also consider all individuals
	  if ( haploid || ( X && (*person)->sex ) )
	    {
	      if ( (!s1) && s2 )
		{
		  hetlist.push_back( (*person)->fid + "\t" + 
				     (*person)->iid + "\t" + 
				     (*loc)->name );
		  
		  // Set to missing, unless in a RECODE mode
		  if ( ! par::preserve_all_genotypes )
		    {
		      s1 = *i1 = true;
		      s2 = *i2 = false;
		    }
		}
	    } 

	  	  
	  
	  // For missing genotypes
	  if ( ! ( s1 && (!s2) ) ) 
	    geno_nm++;
	  
	  // But is this a real genotype in any case?
	  if ( par::oblig_missing )
	    {
	      if ( ! obligMissing(i,l) )
		++geno_real;
	    }
	  else
	    ++geno_real;
	  

	  // Do not recount alleles if we have read in allele frequencies
	  if (!par::af_read)
	    {	      
	      // For allele frequencies
	      // only consider founders	
	      if ( par::summ_nonfounders || (*person)->founder ) 
		{
		  
		  if ( haploid || ( X && (*person)->sex ) )
		    {
		      
		      //////////////////
		      // Haploid counts
		      
		      // "1" allele count
		      
		      if ( (!s1) && (!s2) )   //  FF = hom(11)
			{
			  (*loc)->freq++;
			  (*loc)->nm++;
			}	
		      else if ( s1 && s2 )   //  TT = hom(22)
			{
			  (*loc)->nm++;
			}
		      
		    }
		  else
		    {
		      
		      //////////////////
		      // Autosomal count
		      
		      // "1" allele count
		      
		      if (!s1)
			{ 
			  if (!s2) //   00 = hom(11)
			    {
			      (*loc)->freq+=2;
			      (*loc)->nm+=2;
			    }	
			  else //   01 = het(12)
			    {
			      (*loc)->freq+=1;
			      (*loc)->nm+=2;
			    }
			}
		      else if ( s2 ) // 11 = hom(22)
			{
			  (*loc)->nm+=2;
			}
		    }
		  
		}
	    }
	  
	  // Next individual
	  ++person;
	  ++i;
	  ++i1;
	  ++i2;
	}
      
  
      ////////////////////////////////
      // Calculate allele frequencies
      
      if (!par::af_read)
	{
	  if ( par::af_count) // Allele counts...
	    {
	      // Use freq to store count (keep as is)	

	      // Use "bp" to store number of allele 2
	      (*loc)->bp = (long int)((*loc)->nm - (*loc)->freq);

	      // Use "pos" to store number missing genotypes
	      (*loc)->pos = geno_real - geno_nm;
	    }
	  else // ... or frequencies
	    {		
	      if ((*loc)->nm>0)
		(*loc)->freq /= (double)(*loc)->nm;
	      else
		{
		  (*loc)->freq = 1; 
		  // If we aren't getting rid of it anyway
		  if ( (double)geno_nm/(double)geno_real >= (1-par::MAX_GENO_MISSING))
		    no_founders_found_list.push_back(*loc);
		}
	    }
	}


      //////////////////////////////////////////
      // Record total proportion of missingness
      
      double snp_genotyping = n>0 ? (double)geno_nm/(double)geno_real : 0;
      total_genotyping += snp_genotyping;


      /////////////////////////////////////////////////
      // Exclude if SNP has too many missing genotypes

      if ( snp_genotyping < (1-par::MAX_GENO_MISSING) )
	{	  
	  *d = true;
	  exc_miss++;
	}
      

      ////////////////////////////////////////////////
      // Make allele1 always the least common allele
            
      if ( par::make_minor_allele && (!par::af_count) && (*loc)->freq > 0.5 ) 
	{
	  
	  // then we need to swap alleles

	  (*loc)->freq = 1 - (*loc)->freq;
	  
	  string tmp = (*loc)->allele2;
	  (*loc)->allele2 = (*loc)->allele1;
	  (*loc)->allele1 = tmp;
	  
	  vector<bool>::iterator i1 = (*s)->one.begin();
          vector<bool>::iterator i2 = (*s)->two.begin();

          while ( i1 != (*s)->one.end() ) 
	    {	      

	      if ( (*i1) == (*i2) ) 
		{
		  *i1 = ! (*i1);
		  *i2 = ! (*i2);
		}
	      
	      i1++;
	      i2++;
	    }
	  
	}
      
      // Next SNP
      ++d;
      ++loc;
      ++l;
      ++s;
    }
  



  /////////////////////////////////////////////////
  // Save list of any heterozygous haploid alleles
  
  if (hetlist.size()>0)
    {
      printLOG(int2str( hetlist.size()) 
	       + " heterozygous haploid genotypes; set to missing\n");
      string f = par::output_file_name + ".hh";
      printLOG("Writing list of heterozygous haploid genotypes to [ " 
	       + f + " ]\n");
      ofstream REM;
      REM.open(f.c_str(), ifstream::out);
      for (int i=0; i<hetlist.size(); i++)
	REM << hetlist[i] << "\n";
      REM.close();      
    }
  hetlist.clear();



  /////////////////////////////////////////////////
  // Save list of SNPs with no founders observed
  
  if (no_founders_found_list.size()>0)
    {
      printLOG(int2str( no_founders_found_list.size()) 
	       + " SNPs with no founder genotypes observed\n");
      printLOG("Warning, MAF set to 0 for these SNPs (see --nonfounders)\n"); 
      string f = par::output_file_name + ".nof";
      printLOG( "Writing list of these SNPs to [ " + f + " ]\n");
      ofstream NOF;
      NOF.open(f.c_str(), ifstream::out);
      for (int i=0; i<no_founders_found_list.size(); i++)
	NOF << no_founders_found_list[i]->name << "\n";
      NOF.close();      
    }
  no_founders_found_list.clear();


  
  //////////////////////////
  // Write allele freq file
  
  if (par::af_write)
    {
      if (par::include_cluster_from_file)
	calcStratifiedAlleleFreqs();
      else
	{
	  ofstream FRQ;
	  string f = par::output_file_name + ".frq";
	  if (par::af_count) f += ".count";

	  if (par::summ_nonfounders)
	    printLOG("Writing allele frequencies (all individuals) to [ " 
		     + f + " ] \n");
	  else
	    printLOG("Writing allele frequencies (founders-only) to [ " 
		     + f + " ] \n");
	  
	  if (par::af_count)
	    printLOG("Display counts rather than frequencies\n");

	  FRQ.open(f.c_str(), ifstream::out);
	  FRQ.precision(4);

	  FRQ << setw(4) << "CHR" << " "
	      << setw(par::pp_maxsnp) << "SNP" << " "
	      << setw(4) << "A1" << " "
	      << setw(4) << "A2" << " ";
	  if (par::af_count)
	    FRQ << setw(6) << "C1" << " "
		<< setw(6) << "C2" << " "
		<< setw(6) << "G0" << "\n";
	  else
	    FRQ << setw(12) << "MAF" << " " 
		<< setw(8) << "NCHROBS" 
		<< "\n";	  
	  
	  vector<Locus*>::iterator loc = locus.begin();
	  while (loc != locus.end() )
	    {
	      string a1 = (*loc)->allele1;
	      string a2 = (*loc)->allele2;
	      if (a1=="") a1="0";
	      if (a2=="") a2="0";
	      FRQ << setw(4)  << (*loc)->chr  << " "
		  << setw(par::pp_maxsnp) << (*loc)->name  << " "
		  << setw(4)  << a1  << " "
		  << setw(4)  << a2  << " ";

	      if (par::af_count)
		{
		  FRQ << setw(6) << int( (*loc)->freq ) << " "
		      << setw(6) << int( (*loc)->bp   ) << " " 
		      << setw(6) << int( (*loc)->pos  ) << "\n";
		}
	      else
		{
		  if ( (*loc)->nm > 0 )
		    FRQ << setw(12) << (*loc)->freq << " ";
		  else
		    FRQ << setw(12) << "NA" << " ";
		  FRQ << setw(8) << (*loc)->nm << "\n";
		}
	      loc++;
	    }
	  FRQ.close();
	  
	}
      
      // Close after we've done alle freqs,
      shutdown();
    }



  /////////////////////////
  // Write HWE statistics
  
  if (par::HWD_test || par::HWD_report)
    {
      
      ofstream HWD;
      if (par::HWD_report)
	{
	  
	  if (par::summ_nonfounders)
	    printLOG("Writing Hardy-Weinberg tests (all individuals) to [ " +
		     par::output_file_name + ".hwe ] \n");
	  else
	    printLOG("Writing Hardy-Weinberg tests (founders-only) to [ " +
		     par::output_file_name + ".hwe ] \n");
	  
	  string f = par::output_file_name + ".hwe";
	  HWD.open(f.c_str(), ifstream::out);
	  
	  HWD.precision(4);
	  
	  HWD << setw(4) << "CHR" << " "
	      << setw(par::pp_maxsnp) << "SNP" << " "
	      << setw(8) << "TEST" << " " 
	      << setw(4) << "A1" << " "
	      << setw(4) << "A2" << " "
	      << setw(20) << "GENO" << " "
	      << setw(8) << "O(HET)" << " "
	      << setw(8) << "E(HET)" << " "
	      << setw(12) << "P" << " "	      
	      << "\n";
	}
      
      int cnt=0, cnt_a=0, cnt_u=0;
      

      ////////////////////////
      // Consider each locus
      
      vector<bool>::iterator d = del.begin();
      vector<Locus*>::iterator loc = locus.begin();

      for ( int l = 0 ; l < locus.size() ; l++ )
	{
	  
	  // Compute p-values for HWE test in cases, controls & all
	  
	  // Only consider founders
	  
	  int a11, a12, a22;
	  int u11, u12, u22;
	  int b11, b12, b22;
	  
	  a11=a12=a22=0;
	  u11=u12=u22=0;
	  b11=b12=b22=0;
	  
	  bool X = false, haploid = false;
	  if (par::chr_sex[(*loc)->chr]) X=true;
	  else if (par::chr_haploid[(*loc)->chr]) haploid=true;
	  
	  
	  ///////////////////////////////////////////////
	  // Iterate over each individual, founders only
	 
	  for ( int i = 0 ; i < sample.size() ; i++ ) 
	    {
	      
	      Individual * person = sample[i];

	      ///////////////////////////////////////////////
	      // Only consider founders, & diploid genotypes
	      
	      if ( par::summ_nonfounders || person->founder )
		if ( ! ( haploid || ( X && person->sex ) ) )
		  {	
		    
		    bool s1 = SNP[l]->one[i];
		    bool s2 = SNP[l]->two[i];

		    // Consider everybody, irrespective of phenotype
		    // (QT, C/C or missing)
		    
		    if (!s1)
		      { 
			if (!s2) b11++;   //   00 = hom(11)
			else b12++;       //   01 = het(12)
		      }
		    else if ( s2 ) b22++; // 11 = hom(22)
		    
		    
		    if (par::bt)  // for binary trait, separately for cases/controls
		      {
			if (person->phenotype == 1)
			  {
			    
			    if (!s1)
			      { 
				if (!s2) u11++;   //   00 = hom(11)
				else u12++;       //   01 = het(12)
			      }
			    else if ( s2 ) u22++; //   11 = hom(22)
			    
			  }
			else if (person->phenotype == 2)
			  {
			    if (!s1)
			      { 
				if (!s2) a11++;   //   00 = hom(11)
			      else a12++;         //   01 = het(12)
			      }
			    else if ( s2 ) a22++; //   11 = hom(22)
			  }
			
		      }
		  }
	    
	      
	      // Next individual
	    }


	  // Allele frequencies
	  double afreq = 0, ufreq = 0, freq = 0;
	  
	  bool include_cases = true;
	  bool include_controls = true;
	  if (par::qt)
	    freq = ( b11 + (double)b12/2.0 ) / (double)( b11+b12+b22 );
	  else
	    {
	      afreq = ( a11 + (double)a12/2.0 ) / (double)( a11+a12+a22 );
	      ufreq = ( u11 + (double)u12/2.0 ) / (double)( u11+u12+u22 );
	      freq =  ( b11 + (double)b12/2.0 ) / (double)( b11+b12+b22 );

	      if ( a11+a12+a22 == 0 ) include_cases = false;
	      if ( u11+u12+u22 == 0 ) include_controls = false;		
	    }
	  
	  
	  if (par::qt)
	    {
	      
	      double p;
	      
	      if (par::HWD_standard)
		{
		  double tot = b11 + b12 + b22;
		  double exp_11 = freq * freq * tot;
		  double exp_12 = 2 * freq * (1-freq) * tot;
		  double exp_22 = (1-freq) * (1-freq) * tot;
		  
		  double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11 
		    + ( (b12-exp_12)*(b12-exp_12) ) / exp_12 
		    + ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;
	      
		  p = chiprobP(chisq,1);
		}
	      else
		p = SNPHWE( b12, b11, b22 );
	      

	      if (par::HWD_report)
		{
		  HWD << setw(4) << (*loc)->chr << " " 
		      << setw(par::pp_maxsnp) << (*loc)->name << " "
		      << setw(8) << "ALL(QT)" << " "
		      << setw(4) << (*loc)->allele1 << " "
		      << setw(4) << (*loc)->allele2 << " "
		      << setw(20) << (int2str(b11)+
				      "/"+int2str(b12)+
				      "/"+int2str(b22)) << " "
		      << setw(8) << (double)b12/(double)(b11+b12+b22) << " "
		      << setw(8) << 2 * freq * (1-freq)  << " ";
		  if ( realnum(p) )
		    HWD << setw(12) << p << "\n";
		  else
		    HWD << setw(12) << "NA" << "\n";
		}
	      
	      if ( p <= par::HWD_limit && p > -1 ) 
		{
		  cnt++;
		  *d = true;
		}
	    }
	  else
	    {
	      // For case/control data

	      double p, p_a, p_u;
	      
	      if (par::HWD_standard)
		{
		  double exp_a11 = afreq * afreq * (a11+a12+a22);
		  double exp_a12 = 2 * afreq * (1-afreq) * (a11+a12+a22);
		  double exp_a22 = (1-afreq) * (1-afreq) * (a11+a12+a22);

		  double exp_u11 = ufreq * ufreq * (u11+u12+u22);
		  double exp_u12 = 2 * ufreq * (1-ufreq) * (u11+u12+u22);
		  double exp_u22 = (1-ufreq) * (1-ufreq) * (u11+u12+u22);
		  
		  double exp_11 = freq * freq * (b11+b12+b22);
		  double exp_12 = 2 * freq * (1-freq) * (b11+b12+b22);
		  double exp_22 = (1-freq) * (1-freq) * (b11+b12+b22);
		  
		  double chisq_a = ( (a11-exp_a11)*(a11-exp_a11) ) / exp_a11 
		    + ( (a12-exp_a12)*(a12-exp_a12) ) / exp_a12 
		    + ( (a22-exp_a22)*(a22-exp_a22) ) / exp_a22 ;
		  
		  double chisq_u = ( (u11-exp_u11)*(u11-exp_u11) ) / exp_u11 
		    + ( (u12-exp_u12)*(u12-exp_u12) ) / exp_u12 
		    + ( (u22-exp_u22)*(u22-exp_u22) ) / exp_u22 ;
		  
		  double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11 
		    + ( (b12-exp_12)*(b12-exp_12) ) / exp_12 
		    + ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;
		  
		  p = chiprobP(chisq,1);
		  p_a = chiprobP(chisq_a,1);
		  p_u = chiprobP(chisq_u,1);
		}
	      else
		{
		  p = SNPHWE( b12, b11, b22 );
		  p_a = SNPHWE( a12, a11, a22 );
		  p_u = SNPHWE( u12, u11, u22 );
		}

	      if (par::HWD_report)
		{
		  
		  HWD << setw(4) << (*loc)->chr << " "
		      << setw(par::pp_maxsnp) << (*loc)->name  << " "
		      << setw(8) << "ALL" << " "		    
		      << setw(4) << (*loc)->allele1 << " "
		      << setw(4) << (*loc)->allele2 << " "
		      << setw(20) 
		      << int2str(b11)+"/"+int2str(b12)+"/"+int2str(b22) << " "
		      << setw(8) << (double)b12/(double)(b11+b12+b22) << " "
		      << setw(8) << 2 * freq * (1-freq)  << " ";
		  if ( p > -1 ) 
		    HWD << setw(12) << p  << "\n";
		  else 
		    HWD << setw(12) << "NA"  << "\n";

		  HWD << setw(4) << (*loc)->chr << " "
		      << setw(par::pp_maxsnp) << (*loc)->name  << " "
		      << setw(8) << "AFF" << " "
		      << setw(4) << (*loc)->allele1 << " "
		      << setw(4) << (*loc)->allele2 << " "
		      << setw(20) 
		      << int2str(a11)+"/"+int2str(a12)+"/"+int2str(a22) << " "
		      << setw(8) << (double)a12/(double)(a11+a12+a22) << " "
		      << setw(8) << 2 * afreq * (1-afreq)  << " ";
	
		  
		  if (include_cases && p_a > -1 )
		    HWD << setw(12) << p_a  << "\n";
		  else
		    HWD << setw(12) << "NA" << "\n";


		  HWD << setw(4) << (*loc)->chr << " "
		      << setw(par::pp_maxsnp) << (*loc)->name  << " "
		      << setw(8) << "UNAFF" << " "
		      << setw(4) << (*loc)->allele1 << " "
		      << setw(4) << (*loc)->allele2 << " "
		      << setw(20) 
		      << int2str(u11)+"/"+int2str(u12)+"/"+int2str(u22) << " "
		      << setw(8) << (double)u12/(double)(u11+u12+u22) << " "
		      << setw(8) << 2 * ufreq * (1-ufreq)  << " ";

		    
		  if (include_controls && p_u > -1 )
		    HWD << setw(12) << p_u  << "\n";
		  else
		    HWD << setw(12) << "NA" << "\n";		  
		  		  
		}
	      
	      // Increase counts: in cases
	      if ( include_cases && p_a < par::HWD_limit && p_a > -1 ) cnt_a++;
	      
	      // Controls (and, if possible, exclude on this value)
	      if ( include_controls &&
		   p_u < par::HWD_limit && p_u > -1 ) 
		{
		  cnt_u++;	      

		  if ( ! par::HWD_filter_on_all ) 
		    { 
		      *d = true;
		      cnt++;
		    }
		}
	      
	      // In total sample, and if needed, exclude here
	      if ( p < par::HWD_limit && p>-1 ) 
		{ 
		  if ( par::HWD_filter_on_all || ! include_controls )
		    {
		      *d = true;
		      cnt++;
		    }
		}
	      



	    }


	  // next locus
	  ++loc;
	  ++d;
	}
      

      // Finish the report...
      if (par::HWD_report)
	HWD.close();
      
      // ...or finish pruning
      
      printLOG( int2str(cnt) + 
		" markers to be excluded based on HWE test ( p <= " +
		dbl2str(par::HWD_limit) + " )\n");
      
      if (par::bt)
	{
	  printLOG("\t" + int2str(cnt_a) + " markers failed HWE test in cases\n");
	  printLOG("\t" + int2str(cnt_u) + " markers failed HWE test in controls\n");
	}
    }
  
  

  
  ///////////////////////////////////////////////////
  // Summary statistics for genotyping/missing rates

  if (par::report_missing)
    {
      
      ///////////////////////////////////////////
      // Report by genotyping rate by individual
      // possibly allowing for obligatory missingness
      
      printLOG( "Writing individual missingness information to [ " +
		par::output_file_name + ".imiss ] \n");
      
      ofstream MIS;
      string f = par::output_file_name + ".imiss";
      MIS.open(f.c_str(), ifstream::out);
      MIS.precision(4);

      MIS << setw(par::pp_maxfid) << "FID" << " "
	  << setw(par::pp_maxiid) << "IID" << " "
	  << setw(10) << "MISS_PHENO" << " "
	  << setw(8) << "N_MISS" << " ";
      MIS << setw(8) << "N_GENO" << " ";
      MIS << setw(8) << "F_MISS" << "\n";
      
      for (int i=0; i<n; i++)
	{
	  
 	  MIS << setw(par::pp_maxfid) << sample[i]->fid << " " 
	      << setw(par::pp_maxiid) << sample[i]->iid << " ";
	  if (sample[i]->missing) MIS << setw(10) << "Y" << " ";
	  else MIS << setw(10) << "N" << " " ;
	  
	  
	  // Sum missingness over all SNPs
	  int m=0;       // Missing SNPs
	  int nsnps=0;   // All non-obligatory missing SNPs

	  bool female = ! sample[i]->sex;

	  if ( ! par::oblig_missing )
	    {
	      for (int l=0; l<locus.size();l++)
		{
		  // Skip female Y chromosomes
		  if ( female && par::chr_Y[locus[l]->chr] )
		    continue;
		  
		  if ( MISSING(i,l) ) ++m;
		  ++nsnps;
		}
	    }
	  else // ... allow oblig missing values
	    {
	      
	      for (int l=0; l<locus.size();l++)
		{
		  // Skip female Y chromosomes
		  if ( female && par::chr_Y[locus[l]->chr] )
		    continue;
		  
		  if ( ! obligMissing(i,l) )
		    {
		      if ( MISSING(i,l) )
			++m;
		      ++nsnps;
		    }
		}
	    }
	  
	  MIS << setw(8) << m << " ";
	  MIS << setw(8) << nsnps << " ";
	  MIS << setw(8) << (double)m/(double)nsnps << "\n";
	  
	}
      MIS.close();
      


      ///////////////////////////////////////////
      // Report by genotyping rate by locus
      // possibly allowing for sample strata
      // possibly allowing for obligatory missingness

      printLOG("Writing locus missingness information to [ " +
	       par::output_file_name +".lmiss ] \n");
      f = par::output_file_name + ".lmiss";
      MIS.open(f.c_str(), ifstream::out);
      MIS.clear();
      MIS.precision(4);

      MIS << setw(4) << "CHR" << " " 
	  << setw(par::pp_maxsnp) << "SNP" << " ";
	  
      if (par::include_cluster_from_file)
	MIS << setw(10) << "CLST" << " ";

      MIS << setw(8) << "N_MISS" << " ";
      
      MIS << setw(8) << "N_GENO" << " ";
      
      if (par::include_cluster_from_file)
	MIS << setw(8) << "N_CLST" << " ";
      
      MIS << setw(8) << "F_MISS" << "\n";
      
      
      for (int l=0; l<locus.size(); l++)
	{
	  
	  Locus * loc = locus[l];
	  bool chrY = par::chr_Y[locus[l]->chr];
	  
	  // nk==1 for basic missingness (i.e. not stratified by
	  // cluster)
	  
	  for (int k=0; k<nk; k++)
	    {
	      
	      MIS << setw(4) << loc->chr << " "
		  << setw(par::pp_maxsnp) << loc->name << " ";
	      
	      if (par::include_cluster_from_file)
		MIS << setw(10) << kname[k] << " ";
	      
	      int m=0;     // Number of missing genotypes
	      int c=0;     // Number of people in cluster
	      int nsnps=0; // Number of actual genotypes in cluster
	      
	      for ( int i=0; i<sample.size(); i++)
		{
		  
		  // Skip female Y chromosome calls
		  if ( chrY && ! sample[i]->sex )
		    continue;

		  if (par::include_cluster_from_file)
		    {
		      if ( sample[i]->sol == k ) 
			{
			  if ( ( ! par::oblig_missing ) ||
			       ( ! obligMissing(i,l) ) )
			    {
			      if ( MISSING(i,l) ) ++m;
			      ++nsnps;
			    }
			  ++c;
			}
		    }
		  else // ... ignore cluster strata
		    {
		      if ( ( ! par::oblig_missing ) ||
			   ( ! obligMissing(i,l) ) )
			{
			  if ( MISSING(i,l) ) ++m;
			  ++nsnps;
			}		      
		    }
		  
		  // Next individual
		}
	      	      
	      MIS << setw(8) << m << " ";

	      if (par::include_cluster_from_file)
		MIS << setw(8) << c << " ";	    

	      MIS << setw(8) << nsnps << " ";	    	      
	      MIS << setw(8) << (double)m / (double)nsnps << "\n";
	    }
	  
	  // Next SNP
	}
      
      MIS.close();

    }
  



  /////////////////////////////////
  // Remove rare SNPs 

   
  loc = locus.begin();
  d = del.begin();

  while ( loc != locus.end() )
    {

      // Note epsilon correction for MAF, due to floating point 
      // issues: only apply to the lower MAF range

      if ( (*loc)->freq < 0 || 
	   (*loc)->freq + par::epsilon < par::min_af ||
	   (*loc)->freq > par::max_af )        	 
	{
	  *d = true;
	  exc_maf++;
	}

      d++;
      loc++;
    }
  
  

  /////////////////////////////////////////
  // Remove SNPs based on thresholds

  if ( locus.size() > 0 ) 
    printLOG("Total genotyping rate in remaining individuals is " 
	     + dbl2str(total_genotyping/(double)locus.size())+"\n");
  printLOG(int2str(exc_miss)+" SNPs failed missingness test ( GENO > "
	   +dbl2str(par::MAX_GENO_MISSING)+" )\n");
  printLOG(int2str(exc_maf)+" SNPs failed frequency test ( MAF < "+dbl2str(par::min_af));
  if (par::max_af < 0.5 ) printLOG(" or MAF > " + dbl2str(par::max_af));
  printLOG(" )\n");


  int tmp = deleteSNPs(del);


  //////////////////////////////////////////
  // Need to make back to individual major?

  if ( ! original_SNP_major )
    SNP2Ind();


  return;

}


void Plink::thinSNPs()
{
  if ( par::thin_param <= 0 || par::thin_param >= 1 )
    error("Parameter for --thin must be 0<x<1");
  printLOG("Thinning SNP list, per-SNP probability of keeping is " + dbl2str( par::thin_param ) + "\n"); 
  int x = 0;
  
  set<Locus*> toKeep;
  for (int l=0;l<nl_all;l++)
    if ( CRandom::rand() < par::thin_param ) 
      toKeep.insert(locus[l]);

  int rem = keepSNPs( toKeep );
  printLOG("Removed " + int2str(rem) + " SNPs from thinning\n");

}

