

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2008 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////

#include <iomanip>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "phase.h"

void Plink::proxyWrapper()
{

  if (!par::SNP_major) 
    Ind2SNP();
  
  if ( par::proxy_glm )
    par::assoc_glm_without_main_snp = true;

  /////////////////////////////////////////////////////
  // Use 'pos' slot in Locus to store genotyping
  // rate information, as we need access to this often

  for (int l=0 ; l<nl_all; l++)
    locus[l]->pos = genotypingRate(*this,l);
  
  string f = par::output_file_name;
  string f2 = par::output_file_name;

  if ( par::proxy_impute )
    {
      f += ".proxy.impute";
      f2 += ".proxy.impute.dosage";
    }
  else if ( par::proxy_error )
    {
      f += ".proxy.genocheck";
    }
  else if ( par::proxy_all && ! par::proxy_full_report )
    {
      if ( par::proxy_CC )
	{
	  if ( par::qt )
	    f += ".qassoc.proxy";
	  else
	    f += ".assoc.proxy";
	}
      else
	f += ".tdt.proxy";
    }
  else 
    f += ".proxy.report";
  
  haplo->HTEST.open(f.c_str(),ios::out);
  haplo->HTEST.precision(3);
  
  if ( par::proxy_record_dosage )
    {
      OUTFILE.open(f2.c_str(),ios::out);
      OUTFILE.precision(3);
    }

  printLOG("\n");

/*  printLOG("Criteria for selecting proxy SNPs:\n");
  printLOG("   Selecting at most " + int2str( par::proxy_snp_filter ) + " proxy SNPs (--proxy-maxsnp)\n");
  printLOG("   Searching up to " + int2str( par::proxy_window ) + " SNPs around reference (--proxy-window)\n");
  printLOG("   Searching within " + dbl2str(par::proxy_kb) +" kb around reference (--proxy-kb)\n");
  printLOG("   Proxy genotype missingness threshold is "+dbl2str(par::proxy_geno)+" (--proxy-geno)\n");
  printLOG("   Proxy MAF threshold is " + dbl2str(par::proxy_maf) +" (--proxy-maf)\n");


  if ( par::proxy_r2_filter )
  {
      printLOG("   Proxy r-sq filters of "
               + dbl2str( par::proxy_r2_filter_A ) + ", "
               + dbl2str( par::proxy_r2_filter_B ) + ", "
               + dbl2str( par::proxy_r2_filter_C ) + " (--proxy-r2)\n");
  }
  else
  {
      printLOG("   No proxy r-sq filter selected (--proxy-no-r2-filter)\n");
  }
*/
  
  printLOG("Criteria for selecting proxy SNPs using frequency based metrics:\n");
  printLOG("For SNPs with MAF above, then below " + dbl2str(par::proxy_planB_threshold) + ", \n");
  printLOG("   Selecting at most " + int2str( par::proxy_snp_filter_planA ) + "," 
	   + int2str( par::proxy_snp_filter_planB )
	   + " proxy SNPs (--proxy-maxsnp)\n");
  printLOG("   Searching up to " + int2str( par::proxy_window_planA ) + ","
	   + int2str( par::proxy_window_planB ) + ","
	   + " SNPs around reference (--proxy-window)\n");
  printLOG("   Searching within " + dbl2str(par::proxy_kb_planA) +","
	   + dbl2str(par::proxy_kb_planB) +","
	   +" kb around reference (--proxy-kb)\n");
  printLOG("   Proxy genotype missingness threshold is "+dbl2str(par::proxy_geno)+" (--proxy-geno)\n");
  printLOG("   Proxy MAF threshold is " + dbl2str(par::proxy_maf) +" (--proxy-maf)\n");
  printLOG("   Proxy r-sq filters of "
	   + dbl2str( par::proxy_r2_filter_A_planA ) + ", "
	   + dbl2str( par::proxy_r2_filter_B_planA ) + ", "
	   + dbl2str( par::proxy_r2_filter_C_planA ) + " and " 
	   + dbl2str( par::proxy_r2_filter_A_planB ) + ", "
	   + dbl2str( par::proxy_r2_filter_B_planB ) + ", "
	   + dbl2str( par::proxy_r2_filter_C_planB )  	 
	   + " (--proxy-r2)\n");
    
  if ( par::proxy_impute ) 
    {
      int reference_panel = 0;
      for ( int i = 0 ; i < n ; i++ )
	if ( sample[i]->missing )
	  reference_panel++;
      if ( reference_panel == 0 )
	error("No reference panel for imputation (i.e. individuals with missing phenotypes)\n");
	  printLOG("Imputation reference panel of " 
		     + int2str(reference_panel) + " individuals\n");
    }
  
  if ( ( ! par::proxy_impute ) && ( ( !par::proxy_all ) || par::proxy_full_report ) ) 
    {
      printLOG("Criteria for selecting proxy subhaplotypes:\n");
      printLOG("   Haplotype frequency threshold " +dbl2str(par::proxy_mhf)+" (--proxy-mhf)\n");
      printLOG("   r-squared threshold to reference " + dbl2str(par::proxy_r2) +" (--proxy-sub-r2)\n");
      printLOG("   Maximum of " + int2str(par::proxy_maxhap) + " SNPs per subhaplotype (--proxy-sub-maxsnp)\n");
    }
  
  printLOG("\n");
  
  printLOG("Writing haplotype-based proxy tests to [ " + f + " ] \n");
  if ( par::proxy_record_dosage )
    printLOG("Writing dosage imputation information to [ " + f2 + " ]\n");
  
  // Either condsider all markers, with reduced output (i.e. 
  // just the haplotype based test for that SNP, or a single 
  // marker with extended outout
  
  if ( par::proxy_all )
    {
      
      set<int> plocus;
      
      if ( par::proxy_all_list ) 
	{
	  // read list of SNPs to use as reference SNP
	  // for proxy
	  
	  map<string,int> mlocus;
	  for (int l=0; l<nl_all; l++)
	    mlocus.insert(make_pair( locus[l]->name,l));
	  
	  ifstream PLIST;
	  checkFileExists(par::proxy_all_list_file);
	  PLIST.open(par::proxy_all_list_file.c_str(), ios::in);
	  while ( ! PLIST.eof() )
	    {
	      string marker ;
	      PLIST >> marker;
	      map<string,int>::iterator mi = mlocus.find( marker );
	      if ( mi != mlocus.end() )
		plocus.insert( mi->second );
	    }
	  PLIST.close();
	  
	}
      
      // And we might want a full report for each?
      // Turn off 'all' mode in that case to get 
      // verbose output
      
      if ( par::proxy_full_report ) 
	par::proxy_all = false;
      else
	{
	  if ( par::proxy_impute )
	    {
	      haplo->HTEST << setw(4) << "CHR" << " " 
			   << setw(par::pp_maxsnp) << "SNP" << " "
			   << setw(4) << "NPRX" << " "
			   << setw(8) << "INFO" << " "
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
		haplo->HTEST << "PROXIES";
	      
	      haplo->HTEST << "\n";
	      
	      
	    }
	  else
	    {
	      //haplo->HTEST.setf(ios::scientific);
	      haplo->HTEST << setw(4) << "CHR" << " " 
			     << setw(par::pp_maxsnp) << "SNP" << " "
			     << setw(12) << "BP" << " "
			     << setw(4) << "A1" << " "
			     << setw(4) << "A2" << " "
			     << setw(10) << "GENO" << " " 
			     << setw(4) << "NPRX" << " "
			     << setw(8) << "INFO" << " ";
	      
	      if ( par::proxy_CC )
		{
		  if ( par::qt || par::proxy_glm ) 
		    {
		      haplo->HTEST << setw(8) << "F" << " ";
		      if ( par::qt ) 
			haplo->HTEST << setw(8) << "BETA" << " ";			
		      else
			haplo->HTEST << setw(8) << "ODDS" << " ";			
		    }
		  else
		    haplo->HTEST << setw(8) << "F_A" << " "
				 << setw(8) << "F_U" << " "
				 << setw(8) << "OR" << " ";		      
		}
	      else if ( par::proxy_TDT )
		haplo->HTEST << setw(8) << "T" << " "
			     << setw(8) << "U" << " "
			     << setw(8) << "OR" << " ";
	      
	      haplo->HTEST << setw(10) << "P" << " ";
	      
	      if ( par::proxy_list_proxies )
		haplo->HTEST << "PROXIES";
	      
	      haplo->HTEST << "\n";	      
	    }
	}
      
      
      ////////////////////////////////////////
      // Set-up cache for LD values
      
      proxyLD.clear();
      
      
	  
      ////////////////////////////////////////
      // Iterate over all SNPs specified 

      int similar = 0;
      
      for ( int l = 0 ; l < nl_all ; l++ ) 
	{
	  

	  //	  cout << "testing " << locus[l]->name << "\n";

	  /////////////////////////
	  // Skip if not on list
	  
	  if ( par::proxy_all_list ) 
	    if ( plocus.find(l) == plocus.end() )
	      continue;
	  

	  /////////////////////////////////////////
	  // Have we already seen an identical SNP?
	  // (SKIP THIS PART FOR NOW -- DOESN'T 
	  // PRACTICALLY ADD MUCH AT ALL)

// 	  if ( l != 0 && identicalSNPs( this, l-1, l ) )
// 	    {
// 	      ++similar;
// 	      continue;
// 	    }
	  

	  if ( ! par::silent )
	    {
	      cout << l+1 << " of " << nl_all << " performed         \r";
	      cout.flush();
	    }
	  
	  

	    ///////////////////////////////////
	    // Perform actual proxy-based test

	    performProxyTests(l);
	    
	    //	    cout << "done!\n";
	  
	    ///////////////////////////////////
            // Clear the LD cache occassionally
	      
	    if ( proxyLD.size() > 50000 )
	      proxyLD.clear();
	    

	}

      if ( ! par::silent )
	cout << "\n";	  

    }
  else
    {
      int l = getMarkerNumber(*this,par::proxy_assoc_snp);
      if ( l < 0 )
	error("Cannot find proxy SNP [ " + par::proxy_assoc_snp + " ]\n");
      performProxyTests(l);
      
    }
  
  haplo->HTEST.close();

  if ( par::proxy_record_dosage )
    OUTFILE.close();

}
