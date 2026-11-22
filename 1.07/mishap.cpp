

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
#include <string>
#include "options.h"
#include "helper.h"
#include "plink.h"
#include "phase.h"
#include "stats.h"

extern Plink * PP;

void Plink::performMisHapTests()
{

  printLOG("\nPerforming haplotype-based tests for non-random missingness\n");
    
  if (!par::SNP_major) 
    Ind2SNP();
  
  // Consider 1 SNP at a time
  // Form haplotypes based on the surrounding SNPs
  // Form phenotype based on the patterns of missingness for test SNP
  // Is there any association?

  bool old_silent = par::silent;

  string f = par::output_file_name + ".missing.hap";
  haplo->HTEST.open(f.c_str(),ios::out);
  haplo->HTEST.precision(3);
  printLOG("Writing haplotype-based missingness results to [ " + f + " ] \n");
 
  haplo->HTEST << setw(par::pp_maxsnp) << "SNP" << " " 
	       << setw(10) << "HAPLOTYPE" << " "
	       << setw(8) << "F_0" << " "
	       << setw(8) << "F_1" << " " 
	       << setw(20) << "M_H1" << " "
	       << setw(20) << "M_H2" << " " 
	       << setw(8) << "CHISQ" << " "
	       << setw(8) << "P" << " "
	       << "FLANKING" << "\n";
 
  //////////////////////////////////
  // Set up a single haplotype entry

  haplo->new_pred_locus.resize(1);
    
  // Make new entry in MAP file
  Locus * loc = new Locus;
  loc->name = "";
  loc->chr = 0;
  loc->pos = 1;
  loc->bp = 1;
  loc->allele1 = "1";
  loc->allele2 = "2";
  
  // Add this single dummy locus to the list
  haplo->new_map.clear();
  haplo->new_map.push_back(loc);


  vector<CSNP*>::iterator s = SNP.begin();
  int l=0;
  
  while ( s != SNP.end() )
    {

      if (!par::silent)
	cout << l+1 << " of " << nl_all << " SNPs tested        \r";
	  

      /////////////////////////////////
      // Currently, skipped haplod SNPs
      
      if (par::chr_sex[locus[l]->chr] || 
	  par::chr_haploid[locus[l]->chr] )
	{
	  s++;
	  l++;
	  continue;
	}


      ////////////////////////////////////////////
      // Form phenotype (missingness for this SNP)

      vector<Individual*>::iterator person = sample.begin();
      vector<bool>::iterator i1 = (*s)->one.begin();
      vector<bool>::iterator i2 = (*s)->two.begin();
      int nmiss=0;

      while ( person != sample.end() )
	{
	  
	  // Missing at test SNP?
	  if ( *i1 && ! *i2 )
	    {
	      nmiss++;
	      (*person)->aff = true;
	    }
	  else
	    (*person)->aff = false;
	  
	  // Include all individuals in this analysis
	  (*person)->missing = false;
	  
	  i1++;
	  i2++;
	  person++;
	}      
      
      
      ////////////////////////////////////////
      // Enough missingness to warrant a test?
 
      if (nmiss<5) 
	{
	  s++;
	  l++;
	  continue;
	}
      

      
      ///////////////////////
      // Form test haplotypes
      
      haplo->reset();
      
      // Add which allele to look for (corresponding to allele1)

      vector<int> tmp;
      
      for ( int i = l - par::mishap_window ; 
	    i <= l + par::mishap_window ; i++ )
	{
	  if ( i >= 0 && i < nl_all && i != l ) 
	    if ( locus[i]->chr == locus[l]->chr )
	      tmp.push_back(i);
	}
      
      haplo->new_pred_locus[0] = tmp;
      haplo->new_map[0] = locus[l];
      
      ///////////////////
      // Phase haplotypes
      
      par::silent = true;
      haplo->new_pred_allele = listPossibleHaplotypes(*this,
						      haplo->new_pred_locus[0]);
      haplo->phaseAllHaplotypes(false,*PP->pperm);
      haplo->hname = locus[l]->name;
      par::silent = old_silent;
      
      ///////////////////////////////////////
      // Test association with each haplotype

      map<int,int> tests;
      
      for (int h=0; h < haplo->nh; h++)
	{

	  if (haplo->f[h] >= par::min_af)
	    {
	      tests.clear();

	      for (int h2=0; h2 < haplo->nh; h2++)
		{
		  if (haplo->f[h2] >= par::min_af)
		    {
		      if (h==h2) 
			{
			  tests.insert(make_pair(h2,0));
			}
		      else tests.insert(make_pair(h2,1));
		    }
		}
	      

	      //////////////////////
	      // Test each haplotype
	      
	      int nt = 2;

	      vector<double> caseN(nt,0);
	      vector<double> controlN(nt,0);
  
	      // Consider each individual
	      for (int i=0; i<n; i++)
		{
		  
		  Individual * person = sample[i];
		  
		  if ( person->aff )
		    {
		      
		      for (int z = 0 ; z < haplo->hap1[i].size(); z++)
			{
			  
			  map<int,int>::iterator i1 = tests.find(haplo->hap1[i][z]);
			  map<int,int>::iterator i2 = tests.find(haplo->hap2[i][z]);
			  
			  if ( i1 != tests.end() )
			    {
			      if (!haplo->ambig[i]) 
				caseN[i1->second]++;
			      else
				caseN[i1->second] += haplo->pp[i][z];
			    }
			  
			  if ( i2 != tests.end() )
			    {
			      if (!haplo->ambig[i]) 
				caseN[i2->second]++;
			      else
				caseN[i2->second] += haplo->pp[i][z];
			    }
			  
			}
		    }
		  // Or control?
		  else 
		    {
		      for (int z = 0 ; z < haplo->hap1[i].size(); z++)
			{
			  
			  map<int,int>::iterator i1 = tests.find(haplo->hap1[i][z]);
			  map<int,int>::iterator i2 = tests.find(haplo->hap2[i][z]);
			  
			  if ( i1 != tests.end() )
			    {
			      if (!haplo->ambig[i]) 
				controlN[i1->second]++;
			      else
				controlN[i1->second] += haplo->pp[i][z];
			    }
			  
			  if ( i2 != tests.end() )
			    {
			      if (!haplo->ambig[i]) 
				controlN[i2->second]++;
			      else
				controlN[i2->second] += haplo->pp[i][z];
			    }
			  
			}  	  
		    }
		  
		  
		} // next individual
	      
	      
	      haplo->HTEST << setw(par::pp_maxsnp) << haplo->hname << " ";
	      
		  
	      // Find test haplotype
	      int hh=0;
	      map<int,int>::iterator i1 = tests.begin();
	      while ( i1 != tests.end() )
		{
		  if ( i1->second == 0 )
		    hh = i1->first;
		  i1++;
		}

	      haplo->HTEST << setw(10) << haplo->haplotypeName(hh) << " ";

	      if ( caseN[0] + caseN[1] == 0 )
		haplo->HTEST << setw(8) << "NA" << " ";
	      else
		haplo->HTEST << setw(8) << caseN[0] / ( caseN[0] + caseN[1] ) << " ";
	      
	      if ( controlN[0] + controlN[1] == 0 )
		haplo->HTEST << setw(8) << "NA" << " ";
	      else
		haplo->HTEST << setw(8) << controlN[0] / ( controlN[0] + controlN[1] ) << " ";
	      
	      haplo->HTEST << setw(20) << (dbl2str(caseN[0])+"/"+dbl2str(controlN[0])) << " "
			   << setw(20) << (dbl2str(caseN[1])+"/"+dbl2str(controlN[1])) << " ";
	      
	      
	      vector<double> rowT(nt);
	      double caseT = 0;
	      double controlT = 0;
	      
	      for (int h=0; h<nt; h++)
		{
		  rowT[h] = caseN[h] + controlN[h];
		  caseT += caseN[h];
		  controlT += controlN[h];       
		}
	      
	      double chi2 = 0;
	      for (int h=0; h<nt; h++)
		{
		  double exp = ( rowT[h] * caseT ) / (caseT + controlT);
		  chi2 += ( ( caseN[h] - exp ) * ( caseN[h] - exp ) ) / exp ; 
		  
		  exp = ( rowT[h] * controlT ) / (caseT + controlT);
		  chi2 += ( ( controlN[h] - exp ) * ( controlN[h] - exp ) ) / exp ;  
		}
	      
	      if ( realnum(chi2) )
		haplo->HTEST << setw(8) << chi2 << " "
			     << setw(8) << chiprobP(chi2,nt-1) << " ";
	      else
		haplo->HTEST << setw(8) << "NA" << " "
			     << setw(8) << "NA" << " ";
	      
	      for (int snps=0; snps<haplo->ns-1; snps++)
		haplo->HTEST << locus[haplo->S[snps]]->name << "|";
	      
	      haplo->HTEST << locus[haplo->S[haplo->ns-1]]->name << "\n";
	      
	      	    
	    } // next haplotype
	}
    
      //////////////////////
      // Heterozygsote test

      double caseHET=0;
      double controlHET=0;
      double caseHOM=0;
      double controlHOM=0;
      
      
      // Consider each individual
      for (int i=0; i<n; i++)
	{
	  
	  Individual * person = sample[i];
	  
	  if (person->aff )
	    {
	      
	      for (int z = 0 ; z < haplo->hap1[i].size(); z++)
		{
		 		  
		  if ( haplo->hap1[i][z] != haplo->hap2[i][z] )
		    {
		      if (!haplo->ambig[i]) 
			caseHET++;
		      else
			caseHET += haplo->pp[i][z];			      
		    }
		  else
		    {
		      if (!haplo->ambig[i]) 
			caseHOM++;
		      else
			caseHOM += haplo->pp[i][z];			      
		      
		    }
		  
		}
	    }
	  // Or control?
	  else 
	    {
	      for (int z = 0 ; z < haplo->hap1[i].size(); z++)
		{
		  
		  if ( haplo->hap1[i][z] != haplo->hap2[i][z] )
		    {
		      if (!haplo->ambig[i]) 
			controlHET++;
		      else
			controlHET += haplo->pp[i][z];			      
		    }
		  else
		    {
		      if (!haplo->ambig[i]) 
			controlHOM++;
		      else
			controlHOM += haplo->pp[i][z];			      
		      
		    }
		  
		}  	  
	    }
	  
	  
	} // next individual
      
      double chi2 = 0;

      double total = caseHET + caseHOM + controlHET + controlHOM;

      double total_case = caseHET + caseHOM;
      double total_control = controlHET + controlHOM;
      double total_het = caseHET + controlHET;
      double total_hom = caseHOM + controlHOM;

      double exp_caseHET    = ( total_case * total_het ) / total;
      double exp_controlHET = ( total_control * total_het ) / total;
      double exp_caseHOM    = ( total_case * total_hom ) / total;
      double exp_controlHOM = ( total_control * total_hom ) / total;
      
      chi2 += ( ( caseHET - exp_caseHET ) * ( caseHET - exp_caseHET ) ) / exp_caseHET
	+ ( ( caseHOM - exp_caseHOM ) * ( caseHOM - exp_caseHOM ) ) / exp_caseHOM
	+ ( ( controlHET - exp_controlHET ) * ( controlHET - exp_controlHET ) ) / exp_controlHET
	+ ( ( controlHOM - exp_controlHOM ) * ( controlHOM - exp_controlHOM ) ) / exp_controlHOM;
	
    

      haplo->HTEST << setw(par::pp_maxsnp) << haplo->hname << " ";
	      
      haplo->HTEST << setw(10) << "HETERO" << " ";
      
      if (  caseHET + caseHOM == 0 )
	haplo->HTEST << setw(8) << "NA" << " ";
      else
	haplo->HTEST << setw(8) << caseHET / ( caseHET + caseHOM ) << " ";
      
      
      if ( controlHET + controlHOM == 0 )
	haplo->HTEST << setw(8) << "NA" << " ";
      else
	haplo->HTEST << setw(8) << controlHET / ( controlHET + controlHOM ) << " ";
      
      haplo->HTEST << setw(20) << (dbl2str(caseHET)+"/"+dbl2str(controlHET)) << " "
		   << setw(20) << (dbl2str(caseHOM)+"/"+dbl2str(controlHOM)) << " ";

      if ( realnum(chi2) )
	haplo->HTEST << setw(8) << chi2 << " "
		     << setw(8) << chiprobP(chi2,1) << " ";
      else
	haplo->HTEST << setw(8) << "NA" << " "
		     << setw(8) << "NA" << " ";
      
      for (int snps=0; snps<haplo->ns-1; snps++)
	haplo->HTEST << locus[haplo->S[snps]]->name << "|";
      
      haplo->HTEST << locus[haplo->S[haplo->ns-1]]->name << "\n";

      
      // Next test SNP

      s++;
      l++;

    } 
	
  if (!par::silent)
    cout << "\n";
      
  
  // Shut down results file
  haplo->HTEST.close();
  
  
}

