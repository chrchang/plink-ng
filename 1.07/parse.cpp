


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
#include <string>
#include <cmath>
#include <cstdlib>

#include "options.h"
#include "helper.h"
#include "stats.h"
#include "nlist.h"

void getOutputFilename(CArgs & a)
{

  ////////////////////////////////////////
  // Insert commands from a file
  
  if (a.find("--script"))
    {      
      string f = a.value("--script");
      a.fromScript(f);
    }

  if (a.find("--rerun"))
    {      
      string f = a.value("--rerun");
      a.fromPriorLog(f);
    }
  

  ////////////////////////////////////////////////////
  // Start processing commands





  if (a.find("--out")) 
    {
      par::output_file_name = a.value("--out");
    }
  
  if (a.find("--gplink")) 
    {
      if ( par::cli )
	error("Cannot specify --gplink and --interactive");
      par::silent = true;
      par::gplink = true;
    }
  
  if (a.find("--silent")) 
    {
      if ( par::cli )
	error("Cannot specify --silent and --interactive");

      par::silent = true;
    }

  if (a.find("--interactive")) 
    {
      if ( a.find("--script") )
	error("Cannot specify --script and --interactive");      
      par::cli = true;
    }
}

void setOptions(CArgs & a)
{


  /////////////////////////////////////
  // Web-based functions
  
  if (a.find("--noweb")) par::web_check = false;
  
  if (a.find("--lookup"))
    {
      par::lookup = true;
      par::lookup_single_snp = true;
      par::lookup_snp = a.value("--lookup");
     }
  
   if (a.find("--lookup-list"))
     {
       par::lookup = true;
       par::lookup_to_file = true;
       par::lookup_single_snp = false;
       par::lookup_snp = a.value("--lookup-list");
     }

   if (a.find("--lookup-save"))
     par::lookup_to_file = true;

   if (a.find("--lookup-gene"))
     {
       par::lookup = true;
       par::lookup_to_file = true;
       par::lookup_gene = true;
       par::lookup_gene_name = a.value("--lookup-gene");
     }

   if (a.find("--lookup-gene-list"))
     {
       par::lookup = true;
       par::lookup_to_file = true;
       par::lookup_multiple_genes = true;
       par::lookup_gene = true;
       par::lookup_gene_name = a.value("--lookup-gene-list");
     }

   if (a.find("--lookup-gene-kb"))
     {
       par::lookup_gene_kb_window = a.value_int("--lookup-gene-kb");       
     }
   
   if (a.find("--lookup-kb"))
     {
       par::lookup_snp_kb_window = a.value_int("--lookup-kb");       
     }

   if (a.find("--lookup2"))
     {
       par::lookup2 = true;       
       par::lookup2_cmd = a.value("--lookup2");
     }

   if (a.find("--lookup-gene2"))
     {
       par::lookup2 = true;
       par::lookup_gene = true;
       par::lookup2_cmd = a.value("--lookup-gene2");
     }


   
   if (a.find("--id-dict"))
     {
       par::idhelp = true;
       par::idhelp_dictionary = a.value("--id-dict");

       if ( a.find("--id-auto-alias") )
	 par::idhelp_auto_alias = true;

       // Just print each exact line matching field=value
       if ( a.find("--id-dump") )
	 {
	   par::idhelp_dump_from_dict = true;
	   par::idhelp_dump_from_dict_cmd = a.value("--id-dump");
	 }

       // Default is to dump whole table
       
       if ( a.find("--id-table"))
	 {
	   par::idhelp_subset = "dump_subtable";
	   par::idhelp_subset_string = a.value("--id-table");
	 }
       
       if ( a.find("--id-lookup"))
	 {
	   par::idhelp_lookup = true;
	   par::idhelp_lookup_string = a.value("--id-lookup");
	 }

       if ( a.find("--id-replace"))
	 {
	   par::idhelp_replace = true;
	   vector<string> s = a.value("--id-replace",3);
	   par::idhelp_replace_string = s[0] + "," + "," + s[1] + "," + s[2];
	 }

       if ( a.find("--id-match"))
	 {
	   if ( a.find("--id-replace"))
	     error("You cannot ID match and ID replace together");
	   
	   par::idhelp_match = true;	   
	   par::idhelp_match_string = a.varValue("--id-match");
	   
	 }
       
       if (a.find("--id-delimit"))
	 {
	   par::idhelp_output_delimit = a.value("--id-delimit");

	   if ( par::idhelp_output_delimit == "tab" )
	     par::idhelp_output_delimit = "\t";

	   if ( par::idhelp_output_delimit == "bar" )
	     par::idhelp_output_delimit = "|";

	   if ( par::idhelp_output_delimit == "semi-colon" )
	     par::idhelp_output_delimit = ";";

	   if ( par::idhelp_output_delimit == "comma" )
	     par::idhelp_output_delimit = ",";

	   if ( par::idhelp_output_delimit == "space" )
	     par::idhelp_output_delimit = " ";
	 }

       if ( a.find("--id-alias"))
	 par::idhelp_list_aliases = true;
       
       if ( par::idhelp_match && par::idhelp_replace )
	 error("Cannot specify --id-replace and --id-match together");
       if ( par::idhelp_replace && par::idhelp_lookup )
	 error("Cannot specify --id-replace and --id-lookup together");
       if ( par::idhelp_replace && par::idhelp_subset )
	 error("Cannot specify --id-replace and --id-table together"); 
       if ( par::idhelp_match && par::idhelp_lookup )
	 error("Cannot specify --id-match and --id-lookup together"); 
     }
   
   if ( a.find("--id-match") && !a.find("--id-dict") )
     {
       par::idhelp = true;
       par::idhelp_match = true;	   
       par::idhelp_no_dict = true;
       par::idhelp_match_string = a.varValue("--id-match");       
     }


   if (a.find("--R"))
     {
       par::run_R_script = true;
       par::R_script = a.value("--R");
     }
   
   if (a.find("--R-port"))
     par::R_port = a.value_int("--R-port");
   
   if (a.find("--R-debug"))
     par::run_R_write_script = true;

   if (a.find("--R-chisq"))
     par::run_R_chisq = true;

   if (a.find("--R-z"))
     par::run_R_z = true;

   if (a.find("--R-nsnps"))
     par::run_R_nsnps = a.value_int("--R-nsnps");

   if (a.find("--seed"))
     par::random_seed = a.value_lui("--seed");      

   ////////////////////////////////////////
   // A data-generation option? 

   bool makedata = false;

   if (a.find("--recode")) 
     {
       makedata = true;
       par::recode = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recode12")) 
     {
       makedata = true;
       par::recode_12 = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recodeHV")) 
     {
       makedata = true;
       par::recode_HV = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recode-lgen")) 
     {
       makedata = true;
       par::recode_long = true;

       if (a.find("--with-reference"))
	 par::recode_long_ref = true;

       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recode-rlist")) 
     {
       makedata = true;
       par::recode_mutlist = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recode-whap")) 
     {
       makedata = true;
       par::recode_whap = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recode-fastphase")) 
     {
       makedata = true;
       par::recode_fastphase = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recode-structure")) 
     {
       makedata = true;
       par::recode_structure = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recode-bimbam")) 
     {
       makedata = true;
       par::recode_bimbam = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recodeAD")) 
     {
       makedata = true;
       par::recode_AD = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }


   if (a.find("--recodeA")) 
     {
       makedata = true;
       par::recode_AD = true;
       par::recode_AD_Aonly = true;
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     }

   if (a.find("--recode-allele"))
     {
       if ( ! par::recode_AD ) 
	 error("You need to specify --recodeA or --recodeAD also");
       par::recode_allele_coding = true;
       par::recode_allele_coding_file = a.value("--recode-allele");
     }

   if (a.find("--make-bed")) 
     { 
       makedata = true;
       par::write_bitfile = true; 
    
       // and unless otherwise specified, set GENO = 1 and MAF = 0
       if (!a.find("--maf")) par::min_af = 0.0;
       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;
       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
     } 
  

   //////////////////////////////////////////////
   // Modifiers of main data generation commands

   if (a.find("--alleleACGT") || a.find("--allele-ACGT")) par::recode_ACGT = true;

   if (a.find("--allele1234") || a.find("--allele-1234")) par::recode_1234 = true;

   if (a.find("--reference-allele")) 
     {
       par::set_reference_allele = true;       
       par::set_reference_allele_file = a.value("--reference-allele");
       par::make_minor_allele = false;
     }

   if (a.find("--transpose"))
     {
       if ( ! ( a.find("--recode") || a.find("--recode12")  ) )
 	error("--transpose requires --recode or --recode12");
       par::recode_transpose = true;
     }

   if (a.find("--tab"))
     {
       if ( ! makedata ) 
 	error("You can only specify --tab with a --recode* option");
       par::recode_delimit = "\t";
     }
  

   ////////////////////////////////////////////////////
   // Convenience functions and obligatory missingness

   // Zero out genotypes for a particular SNP / cluster

   if ( a.find("--zero-cluster") )
     {
       if ( ! a.find("--within") )
 	error("You must specify --within with --zero-cluster");
       par::zero_cluster = true;
       par::zero_cluster_filename = a.value("--zero-cluster");
       checkFileExists(par::zero_cluster_filename);
     }

   // Perform genotyping rate calculations, but allow for 
   // obligatory missingness
  
   if ( a.find("--oblig-missing") )
     {
       if ( ! a.find("--oblig-clusters") )
 	error("You must specify --oblig-clusters with --oblig-missing");
       par::oblig_missing = true;
       par::oblig_missing_filename = a.value("--oblig-missing");
       checkFileExists(par::oblig_missing_filename);
     }

   if ( a.find("--oblig-clusters") )
     {
       if ( ! a.find("--oblig-missing") )
 	error("You must specify --oblig-missing with --oblig-clusters");
       par::oblig_clusters_filename = a.value("--oblig-clusters");
       checkFileExists(par::oblig_clusters_filename );	      
     }

   // Multiple category phenotype

   if ( a.find("--loop-assoc") )
     {
       if ( a.find("--within") )
 	error("You cannot specify --within with --loop-assoc");
    
       par::all_pheno = true;
       par::assoc_test = true;
       par::loop_over = true;
       par::include_cluster_from_file = true;
       par::include_cluster_filename = a.value("--loop-assoc");
       checkFileExists(par::include_cluster_filename);
     }

  

   ////////////////////////////////////////
   // TODO list

   if (a.find("--todo"))
     {
       cout << "TODO list\n"
	    << " * Make --min/--max thresholds apply w/ read-segment\n"
	    << " * Odds ratio for DFAM\n"
	    << " * Check Hotelling's T2 imputation of missing genotypes\n"
	    << " * Add improved familial phasing to haplotype tests\n"
	    << " * Add OR and CI to haplotype tests\n"
	    << " * Add CMH and/or permutation to haplotype tests\n"
	    << " \n";
       exit(0);
     }


   ///////////////////////////////////////
   // Output options

   if (a.find("--flag")) par::flag = true;

   if (a.find("--verbose")) par::verbose = true;

   if (a.find("--pedigree")) par::dumpped = true;

   if (a.find("--tucc")) par::tucc = true;

   if (a.find("--debug")) par::debug = true;


   ///////////////////////////////////////
     // IBD analyses, output options
   
     if (a.find("--multi")) par::multi_output = true;
  
     if (a.find("--gmulti")) par::gmulti_output = true;
  
     if (a.find("--summarise-ibd")) par::summary_ibd_output = true;
     
     if (a.find("--genome") || 
	 a.find("--Z-genome") ||
	 a.find("--genome-minimal") )
       {
	 // By default, include everybody
	 if (!a.find("--min"))
	   {
	     par::MIN_PIHAT = 0;
	     par::pihat_filter = false;
	   }
	 
	 par::genome_output = true;
	 
	 if (a.find("--Z-genome"))
	   par::compress_genome = true;
	 
	 if (a.find("--genome-minimal"))
	   par::genome_output_minimal = true;
	 else if (a.find("--genome-full"))
	   par::genome_output_full = true;	 
       }

     
     if (a.find("--rel-check"))
       par::genome_only_check_rels = true;
     
     
   if (a.find("--impossible")) 
     {
       if (!a.find("--genome"))
 	error("Can only specify --impossible with --genome");
       par::show_impossible_IBD = false;
     }

   if (a.find("--nudge")) 
     {
       if (!a.find("--genome"))
 	error("Can only specify --nudge with --genome");
       if (a.find("--impossible"))
 	error("Cannot specify --impossible and --nudge together");
       par::nudge = true;
     }

   if (a.find("--unbounded")) 
     { 
 	par::bound = false; 
 	par::MIN_PIHAT = -1;
 	par::MAX_PIHAT = 1;
     }

   if (a.find("--genome-lists")) 
     {
       if ( ! ( a.find("--genome") || a.find("--segment") ) )
 	error("Must specify --genome or --segment with --genome-lists");      
       par::genome_2sets = true;
       vector<string> s = a.value("--genome-lists",2);
       par::genome_setlist1 = s[0];
       par::genome_setlist2 = s[1];
     }

   if (a.find("--segment-within"))
     {
       if ( !a.find("--within"))
 	error("You need to specify a --within {clusterfile}");
       par::IBD_within = true;
     }

   if (a.find("--genome-test")) 
     {
       par::genome_test = true;
       par::genome_test_threshold = a.value_double("--genome-test");
       par::plink = true;
       par::permute = true;
       par::adaptive_perm = false;
       par::replicates = 100000;
     }

   if (a.find("--ibs-test")) 
     {
       par::ibs_test = true;
       par::permute = true;
       par::adaptive_perm = false;
       par::replicates = 100000;
     }

   if (a.find("--ibs-test2")) 
     {
       par::ibs_test = true;
       par::ibs_test_method2 = true;
       par::permute = true;
       par::adaptive_perm = false;
       par::replicates = 100000;
     }


   if (a.find("--segment-match-snp")) 
     par::genome_test_min_snp = a.value_int("--segment-match-snp");
  



   ///////////////////////////////////////////
   // WGAS main options: summary stats and QC


   if (a.find("--missing")) { par::report_missing = true; } 

   if (a.find("--mendel")) { par::MENDEL_report = true; } 

   if (a.find("--hardy")) par::HWD_report = true;

   if (a.find("--hardy2")) par::HWD_report = par::HWD_standard = true;

   if (a.find("--het")) par::inbreeding = true;

   if (a.find("--Fst") || a.find("--fst"))
     {
       if ( ! a.find("--within") )
	 error("Need to specify --within with --Fst");
       par::calcFst = true;
     }

   if (a.find("--check-sex"))
     par::check_sex = true;

   if (a.find("--impute-sex"))
     par::check_sex = par::impute_sex = true;

   if (a.find("--test-missing")) { par::test_missing = true; } 

   if (a.find("--test-mishap")) { par::mishap_test = true; } 

   if (a.find("--mishap-window")) 
     { 
       par::mishap_test = true; 
       par::mishap_window = a.value_int("--mishap-window");
     } 

   if (a.find("--score"))
     {
       par::score_risk = true;
       par::score_risk_file = a.value("--score");
       
       if (a.find("--score-no-mean-imputation"))
	 par::score_impute_expected = false;

       if (a.find("--score-ranges"))
	 {
	   par::score_risk_ranges = true;
	   par::score_risk_ranges_file = a.value("--score-ranges");
	 }

       if (a.find("--score-ranges-min"))
	 par::score_risk_ranges_min = a.value_int("--score-ranges-min");

       if (a.find("--score-ranges-border"))
	 par::make_set_border = a.value_int("--score-ranges-border") * 1000;      

       if (a.find("--q-score-file"))
	 {
	   if ( ! a.find("--q-score-range") )
	     error("Must specify --q-score-range with --q-score-file");
	   
	   par::score_risk_on_qrange = true;
	   par::score_qrange_file = a.value("--q-score-range");
	   par::score_qfile = a.value("--q-score-file");
	 }
       
       if ( a.find("--score-test") )
	 {
	   par::score_test = true; 
	 }

       if ( a.find("--set") || a.find("--make-set") ) 
	 par::profile_sets = true;

     }


  ////////////////////////////////////////////////////////////////
  // Proxy association methods

  if (a.find("--proxy-assoc"))
    {
      par::proxy_assoc = true;
      par::proxy_CC = true;
      par::proxy_assoc_snp = a.value("--proxy-assoc");
      if ( par::proxy_assoc_snp == "all" || par::proxy_assoc_snp == "ALL" )
 	par::proxy_all = true;
      
      if ( a.find("--proxy-glm") )
 	par::proxy_glm = true;
    }
  
  
  if (a.find("--proxy-drop"))
    {
      if ( a.find("--proxy-impute"))
 	error("Cannot have --proxy-drop and --proxy-impute (already implied)");
      par::proxy_leave_out = true;
    }
  
  if (a.find("--proxy-exclude"))
    {
      par::proxy_exclude = true;
      par::proxy_exclude_from_file = false;
      par::proxy_exclude_list = a.value("--proxy-exclude");
    }
  else if (a.find("--proxy-exclude-list"))
     {
       par::proxy_exclude = true;
       par::proxy_exclude_from_file = true;
       par::proxy_exclude_list = a.value("--proxy-exclude-list");
     }
  
  if (a.find("--proxy-include-reference"))
    par::proxy_include_reference = true;
  
  if (a.find("--proxy-impute"))
    {
      par::proxy_assoc = true;
      par::proxy_impute = true;
      par::proxy_assoc_snp = a.value("--proxy-impute");
      if ( par::proxy_assoc_snp == "all" || par::proxy_assoc_snp == "ALL" )
 	par::proxy_all = true;
      
      if ( a.find("--proxy-dosage") )
 	par::proxy_record_dosage = true;
      
      if ( a.find("--proxy-replace") )
 	par::proxy_impute_replace = true;
      
      //       if ( a.find("--proxy-preserve") )
      // 	par::proxy_impute_preserve_genotyped = true;
      
      if ( a.find("--proxy-genotypic-concordance") )
 	par::proxy_impute_genotypic_concordance = true;
      
    }
  
  if (a.find("--proxy-error")) par::proxy_error = true;
  
  
  if (a.find("--proxy-impute-threshold"))
    par::proxy_impute_threshold = a.value_double("--proxy-impute-threshold");
  
  if (a.find("--proxy-tdt"))
    {
      if (a.find("--hap-tdt"))
	error("Cannot specify --proxy-tdt and --hap together");
      
      par::proxy_assoc = true;
      par::proxy_TDT = true;
      par::proxy_assoc_snp = a.value("--proxy-tdt");
      if ( par::proxy_assoc_snp == "all" || par::proxy_assoc_snp == "ALL" )
 	par::proxy_all = true;
    }

  if ( par::proxy_assoc ) 
    {
      // Reset default missing rate per individual
      // per haplotype (more relevant if the window
      // size is small (i.e. to goal is not to discard
      // too many / any individuals, but rely on E-M 
      // phasing to reconstruct missing genotypes. Normally
      // the default value here is 0.5;  it can be modified as
      // --hap-miss option is below this one
      
      par::hap_missing_geno = 0.9;
      
    }
  
  if (a.find("--proxy-verbose"))
    par::proxy_full_report = true;
  

   if (a.find("--proxy-flanking"))
     {
       if ( par::proxy_all ) 
 	error("Cannt specify --proxy-flanking with >1 reference SNP");

       par::proxy_list = true;
       par::proxy_list_file = a.value("--proxy-flanking");
     }

   if (a.find("--proxy-list"))
     {
       if ( a.find("--proxy-flanking"))
 	error("Cannt specify --proxy-flanking with >1 reference SNP");

       par::proxy_all_list = true;
       par::proxy_all_list_file = a.value("--proxy-list");
     }


   if (a.find("--proxy-sub-r2"))
     par::proxy_r2 = a.value_double("--proxy-sub-r2");

   if (a.find("--proxy-maf"))
     par::proxy_maf = a.value_double("--proxy-maf");

   if (a.find("--proxy-geno"))
     par::proxy_geno = a.value_double("--proxy-geno");

   if (a.find("--proxy-mhf"))
     par::proxy_mhf = a.value_double("--proxy-mhf");

   if (a.find("--proxy-sub-maxsnp"))
     par::proxy_maxhap = a.value_int("--proxy-sub-maxsnp");

   if (a.find("--proxy-no-r2-filter"))
     par::proxy_r2_filter = false;

   if (a.find("--proxy-show-proxies"))
     par::proxy_list_proxies = true;

   if (a.find("--proxy-r2-reference-only"))
     par::proxy_reference_only = true;

   //////////////////////////////////////
   // Frequency dependent proxy filters
   
   // Plans A and B

   if (a.find("--proxy-r2"))
     {
       par::proxy_r2_filter = true;
       vector<string> s = a.value("--proxy-r2",3);
       par::proxy_r2_filter_A_planA = getDouble(s[0], "--proxy-r2");
       par::proxy_r2_filter_B_planA = getDouble(s[1], "--proxy-r2");
       par::proxy_r2_filter_C_planA = getDouble(s[2], "--proxy-r2");           
     }

   if (a.find("--proxy-maxsnp"))
     par::proxy_snp_filter_planA = a.value_int("--proxy-maxsnp");

   if (a.find("--proxy-window"))
     par::proxy_window_planA = a.value_int("--proxy-window");

   if (a.find("--proxy-kb"))
     par::proxy_kb_planA = a.value_double("--proxy-kb");
   

   // And plan B

   if (a.find("--proxy-b-threshold"))
     par::proxy_planB_threshold = a.value_double("--proxy-b-threshold");
   
   if (a.find("--proxy-b-r2"))
     {
       par::proxy_r2_filter = true;
       vector<string> s = a.value("--proxy-b-r2",3);
       par::proxy_r2_filter_A_planB = getDouble(s[0], "--proxy-b-r2");
       par::proxy_r2_filter_B_planB = getDouble(s[1], "--proxy-b-r2");
       par::proxy_r2_filter_C_planB = getDouble(s[2], "--proxy-b-r2");           
     }

   if (a.find("--proxy-b-maxsnp"))
     par::proxy_snp_filter_planB = a.value_int("--proxy-b-maxsnp");

   if (a.find("--proxy-b-window"))
     par::proxy_window_planB = a.value_int("--proxy-b-window");

   if (a.find("--proxy-b-kb"))
     par::proxy_kb_planB = a.value_double("--proxy-b-kb");




  ////////////////////////////////////////
  // Segmental options
  
  //  --segment-match-snp
  //  --homozyg
  //  --homozyg-window-kb
  //  --homozyg-window-snp
  //  --homozyg-window-het
  //  --homozyg-window-missing
  //  --homozyg-snp
  //  --homozyg-kb
  //  --homozyg-density
  //  --homozyg-gap
  //  --homozyg-group
  //  --homozyg-match
  //  --homozyg-het
  //  --homozyg-verbose
   
  //  --segment
  //  --segment-gap
  //  --segment-length
  //  --segment-thresholds
  //  --segment-minimal

  //  --segment-group
  //  --segment-spanning
  //  --segment-from
  //  --segment-to
  //  --segment-force
  //  --segment-match
  //  --segment-verbose

  //  --pool-size 

  // Read segments back in
  //  --read-segment
  //  --read-segment-minimal


  //////////////////////////
  // CNV/other segment types
 
  if (a.find("--cnv-list"))
    {
      par::cnv_list = true;
      par::cnv_listname = a.value("--cnv-list");      
    }
  
  if (a.find("--cfile"))
    {
      if ( a.find("--cnv-list"))
 	error("Cannot specify --cfile and --cnv-list together");
      if ( a.find("--fam"))
 	error("Cannot specify --cfile and --fam together");
      if ( a.find("--map"))
 	error("Cannot specify --cfile and --map together");
      
      par::cnv_list = true;
      par::fileroot = a.value("--cfile");
      par::cnv_listname = par::fileroot + ".cnv";
      par::famfile = par::fileroot + ".fam";
      par::mapfile = par::fileroot + ".cnv.map";
    }
  
  if ( par::cnv_list )
    {
      // Remove any missing individuals
      if (!a.find("--cnv-missing-phenotypes"))
	par::ignore_phenotypes = false;
      
      if (a.find("--cnv-make-map"))
	par::cnv_makemap = true;

      if (a.find("--cnv-write"))
	{
	  par::cnv_writelist = true;
	  if (a.find("--with-phenotype"))
	    par::dump_covar_with_phenotype = true;
	}
      
      if ( a.find("--cnv-disrupt" ) )
	{
	  if (a.find("--cnv-overlap") || a.find("--cnv-union-overlap") || a.find("--cnv-region-overlap"))
	    error("Cannot specify --cnv-overlap and --cnv-disrupt together");
	  par::cnv_disrupt = true;
	} 

      if (a.find("--cnv-intersect"))
	{
	  if ( a.find("--cnv-exclude") || 
	       a.find("--cnv-count") )
	    error("Cannot specify --cnv-count/intersect/exclude/disrupt together");
	  par::cnv_intersect = true;
	  par::cnv_intersect_file = a.value("--cnv-intersect");
	}

      if (a.find("--cnv-subset"))
	{
	  if ( ! ( a.find("--cnv-exclude") || 
		   a.find("--cnv-count") ||
		   a.find("--cnv-intersect")) )
	    error("Must use --cnv-intersect/exclude/count with --cnv-subset");
	  par::cnv_intersect_subset = true;
	  par::cnv_intersect_subset_file = a.value("--cnv-subset");
	}
      
      
      if (a.find("--cnv-check-no-overlap"))
	par::cnv_check_overlap = true;
      
      if (a.find("--cnv-border"))
	{
	  par::cnv_region_border = 1000 * a.value_int("--cnv-border");
	}

      if (a.find("--cnv-exclude"))
	{
	  if ( a.find("--cnv-count") || 
	       a.find("--cnv-intersect"))
	    error("Cannot specify --cnv-count/intersect/exclude/disrupt together");

	  par::cnv_intersect = par::cnv_exclude = true;
	  par::cnv_intersect_file = a.value("--cnv-exclude");
	}

      if (a.find("--cnv-count"))
	{
	  if ( a.find("--cnv-exclude") || 
	       a.find("--cnv-intersect"))
	    error("Cannot specify --cnv-count/intersect/exclude together");
	  par::cnv_intersect = par::cnv_count = true;
	  par::cnv_intersect_file = a.value("--cnv-count");

	  if ( a.find("--cnv-count-baseline"))
	    {
	      par::cnv_count_baseline = true;
	      par::cnv_count_baseline_file = a.value("--cnv-count-baseline");
	    }
	  
	  if ( a.find("--cnv-weighted-count"))
	    par::cnv_weighted_gene_test = true;
	  
	}
      
      if (a.find("--cnv-freq-method2"))
	{
	  par::cnv_freq_method2 = true;
	  par::cnv_freq_method2_threshold = a.value_double("--cnv-freq-method2");

	  if (a.find("--cnv-unique"))
	    error("Cannot specify --cnv-unique and --cnv-method2 together");
	  
	  if (a.find("--cnv-write-freq"))
	    par::cnv_write_freq = true;
	}


      if (a.find("--cnv-freq-exclude-above"))
	{
	  par::cnv_freq_include = true;
	  par::cnv_freq_include_below = true;
	  par::cnv_freq_include_cnt = a.value_int("--cnv-freq-exclude-above");
	}
      
      if (a.find("--cnv-freq-exclude-below"))
	{
	  par::cnv_freq_include = true;
	  par::cnv_freq_include_below = false;
	  par::cnv_freq_include_cnt = a.value_int("--cnv-freq-exclude-below");
	}      

      if (a.find("--cnv-freq-include-exact"))
	{
	  par::cnv_freq_include = true;
	  par::cnv_freq_include_exact = true;
	  par::cnv_freq_include_cnt = a.value_int("--cnv-freq-include-exact");
	}

      if (a.find("--cnv-freq-exclude-exact"))
	{
	  par::cnv_freq_include = true;
	  par::cnv_freq_include_exact = true;
	  par::cnv_freq_include_exact_exclude = true;
	  par::cnv_freq_include_cnt = a.value_int("--cnv-freq-exclude-exact");
	}


      if (a.find("--cnv-unique"))
	{
	  par::cnv_unique = true;
	}
      
      if (a.find("--cnv-report-regions") )
	{
	  if ( ! ( a.find("--cnv-intersect") || 
		   a.find("--cnv-exclude") || 
		   a.find("--cnv-disrupt") ) )
	    error("Must specify --cnv-intersect/exclude/disrupt with --cnv-report-regions");
	  par::cnv_intersect_writeback = true;
	}
      else if (a.find("--cnv-verbose-report-regions") )
	{
	  if ( ! ( a.find("--cnv-intersect") || 
		   a.find("--cnv-exclude") || 
		   a.find("--cnv-disrupt") ) )
	    error("Must specify --cnv-intersect/exclude/disrupt with --cnv-report-regions");
	  par::cnv_intersect_writeback = true;
	  par::cnv_intersect_writeback_verbose = true;
	}

      
      

      if (a.find("--cnv-overlap"))
	{
	  par::cnv_defined_overlap = true;
	  par::cnv_overlap = a.value_double("--cnv-overlap");
	}
      else if (a.find("--cnv-union-overlap"))
	{
	  par::cnv_defined_overlap = true;
	  par::cnv_overlap = a.value_double("--cnv-union-overlap");
	  par::cnv_union_overlap = true;
	}
      else if (a.find("--cnv-region-overlap"))
	{
	  par::cnv_defined_overlap = true;
	  par::cnv_overlap = a.value_double("--cnv-region-overlap");
	  par::cnv_region_overlap = true;
	}

      if (a.find("--cnv-kb"))
	par::cnv_min_kb = a.value_int("--cnv-kb");
      if (a.find("--cnv-score"))
	par::cnv_min_score = a.value_double("--cnv-score");
      if (a.find("--cnv-sites"))
	par::cnv_min_sites = a.value_int("--cnv-sites");

      if (a.find("--cnv-max-kb"))
	par::cnv_max_kb = a.value_int("--cnv-max-kb");
      if (a.find("--cnv-max-score"))
	par::cnv_max_score = a.value_double("--cnv-max-score");
      if (a.find("--cnv-max-sites"))
	par::cnv_max_sites = a.value_int("--cnv-max-sites");
      
      if (a.find("--cnv-del"))
	par::cnv_del_only = true;
      if (a.find("--cnv-dup"))
	par::cnv_dup_only = true;
     
      if (a.find("--cnv-test-window"))
	{
	  par::seg_test_window = true;
	  par::seg_test_window_bp = a.value_double("--cnv-test-window") * 1000;
	}
      
      if (a.find("--cnv-test-region"))
	{
	  if ( par::seg_test_window )
	    error("Cannot specify both --cnv-test-window and --cnv-test-region");
	  par::seg_test_region = true;
	}

      if (a.find("--cnv-test-2sided"))
	par::segment_test_1sided = false;
      
      if (a.find("--cnv-test-1sided"))
	par::segment_test_force_1sided = true;
      
//       if (a.find("--cnv-glm"))
// 	par::cnv_glm = true;
      
      if (a.find("--cnv-indiv-perm"))
	par::cnv_indiv_perm = true;

      if (a.find("--cnv-enrichment-test"))
	{
	  if ( ! a.find("--cnv-count") )
	    error("The --cnv-enrichment-test option requires --cnv-count");
	  par::cnv_enrichment_test = true;
	  if ( a.find("--cnv-model") )
	    par::cnv_en_model = a.value_int("--cnv-model");
	}

      if (a.find("--cnv-position-perm"))
	par::cnv_pos_perm = true;
      
      if (a.find("--cnv-seglist"))
	par::display_segment_long = true;

      if (a.find("--cnv-track"))
	par::display_cnv_track = true;

      if (a.find("--cnv-blue"))
	par::cnv_col = 1;

      if (a.find("--cnv-green"))
	par::cnv_col = 2;

      if (a.find("--cnv-red"))
	par::cnv_col = 3;

      if (a.find("--cnv-brown"))
	par::cnv_col = 4;


      if (a.find("--cnv-drop-no-segment"))
	par::cnv_drop_no_segment = true;

      if (a.find("--merge") || a.find("--bmerge"))
	error("Cannot specify --cnv-list and merge SNP data together");
 
      if (a.find("--exclude") || a.find("--extract"))
	error("Cannot specify --cnv-list and exclude/extract markers SNP data together");

      if (a.find("--maf") || a.find("--geno") || a.find("--mind"))
	error("Cannot specify --cnv-list and filter on SNP data together");

    }

  ///////////////////////
  // Runs of homozygosity

  if (a.find("--homozyg"))
    {
      par::homo_run = true;
    }

  if (a.find("--read-homozyg"))
    {
      par::homo_run = true;
      par::read_segment_filename = a.value("--read-homozyg");
      par::read_segment_file = true;
    }

  if (a.find("--homozyg-snp"))
    {
      par::homo_run = true;
      par::homo_run_snps = true;
      par::homo_run_length_snps = a.value_int("--homozyg-snp");
    }

  if (a.find("--homozyg-kb")) 
  { 
      par::homo_run = true;
      par::homo_run_kb = true;
      par::homo_run_length_kb = a.value_int("--homozyg-kb");
    }

  if (a.find("--homozyg-density")) 
    { 
      par::homo_run = true;
      par::homo_run_density = a.value_double("--homozyg-density");
    }

  if (a.find("--homozyg-gap")) 
    { 
      par::homo_run = true;
      par::homo_run_gap = a.value_int("--homozyg-gap");
    }

  if (a.find("--homozyg-window-snp")) 
    { 
      par::homo_run = true;
      par::homo_windowSize = a.value_int("--homozyg-window-snp");
    }
  
  if (a.find("--homozyg-window-kb")) 
    { 
      par::homo_run = true;
      par::homo_windowKB = a.value_int("--homozyg-window-kb");
    }
    
  if (a.find("--homozyg-window-het")) 
    { 
      par::homo_run = true;
      par::homo_windowAllowedHet = a.value_int("--homozyg-window-het");
    }
    
  if (a.find("--homozyg-window-missing")) 
    { 
      par::homo_run = true;
      par::homo_windowAllowedMissing = a.value_int("--homozyg-window-missing");
    }

  if (a.find("--homozyg-window-threshold")) 
    { 
      par::homo_run = true;
      par::homo_threshold = a.value_double("--homozyg-window-threshold");
    }

  if (a.find("--homozyg-group"))
    {
      par::homo_summary_allelic_match = true;
      par::fuzzy_homo = 0.95;
    }
    
  if (a.find("--homozyg-match"))
    {
      par::homo_summary_allelic_match = true;
      par::fuzzy_homo = a.value_double("--homozyg-match");
    }

  if (a.find("--homozyg-het")) 
    { 
      if (! ( a.find("--homozyg-snp") ||
	      a.find("--homozyg-kb") ) )
	error("Must specify --homozyg-snp or --homozyg-kb with --homozyg-het");
      par::homo_run_het = a.value_int("--homozyg-het");
    }

  if (a.find("--homozyg-verbose"))
    par::homozyg_verbose = true;

  if (a.find("--consensus-match"))
    {
      par::homo_run_consensus_match = true;
    }

  if (a.find("--homozyg-include-missing"))
    {
      par::homo_miss_as_hom = true;
    }

  if (a.find("--pool-size"))
    {
      par::pool_size_min = a.value_int("--pool-size");
    }
  
  if (a.find("--ibs")) par::ibs_run = true;

  if (a.find("--ibs2")) par::ibs_run = par::ibs_2only = true;

  if (a.find("--ibs-density")) 
    {
      par::ibs_run = true;
      par::ibs_run_density = a.value_double("--ibs-density");
    }

  if (a.find("--ibs-kb")) 
    {
      par::ibs_run = true;
      par::ibs_run_length_kb = a.value_int("--ibs-kb");
    }

  if (a.find("--ibs-snp")) 
    {
      par::ibs_run = true;
      par::ibs_run_length_snps = a.value_int("--ibs-snp");
    }

//   static int ibs_inner_run_length_kb;

//   if (a.find("--ibs-join-snp")) 
//     {
//       par::ibs_run = true;
//       par::ibs_join_snp = a.value_int("--ibs-join-snp");
//     }

//   if (a.find("--ibs-join-snp")) 
//     {
//       par::ibs_run = true;
//       par::ibs_join_snp = a.value_int("--ibs-join-snp");
//     }

//   if (a.find("--ibs-join-kb")) 
//     {
//       par::ibs_run = true;
//       par::ibs_join_kb = a.value_int("--ibs-join-kb") * 1000;
//     }


  if (a.find("--ibs-gap")) 
    {
      par::ibs_run = true;
      par::ibs_inter_snp_distance = a.value_int("--ibs-gap") * 1000;
   }

  if (a.find("--ibs-miss")) 
    {
      par::ibs_run = true;
      par::ibs_run_missing = a.value_int("--ibs-miss");
    }

  if (a.find("--ibs-err")) 
    {
      par::ibs_run = true;
      par::ibs_run_0 = a.value_int("--ibs-err");
    }
  
  if (a.find("--miss-run-snps"))
    {
      par::miss_run = true;
      par::miss_run_length = a.value_int("--miss-run-snps");
    }
  
  if (a.find("--miss-run")) 
    { 
      cerr << "\n*** WARNING -- use --miss-run-snps N option instead\n\n";
      par::miss_run = true;
      par::miss_run_length_kb = true;
      par::miss_run_length = a.value_int("--miss-run");
    }

  if (a.find("--miss-run-level"))
    {
      par::miss_run = true;
      par::miss_run_level = a.value_double("--miss-run-level");
    }


  ///////////////////////////
  // Segmental sharing tests
  
//   if (a.find("--plink")) 
//     {
//       par::plink = true;
//       par::nudge = true;
//     }

  if (a.find("--segment")) 
    {
      par::plink = true;
      par::nudge = true;
      par::segment_output = true;
    }

  if (a.find("--segment-ibs"))
    {
      par::segment_validate = true; 
    }

  if (a.find("--segment-minimal")) 
    {
      par::plink = true;
      par::nudge = true;
      par::segment_output = true;
      par::segment_minimal = true;
    }
  
  if (a.find("--segment-test-individual"))
    {
      // Instead of pairwise Case-Case versus non-Case-Case, 
      // just count the number of segs that any one individual
      // has, and compare count in cases to count in controls
      par::segment_test_individual = true;
      if (a.find("--specific"))
	par::segment_test_specific_segs = true;
    }

  if (a.find("--segment-test-ignore-discordant"))
    {
      // Pairwise Case-Case versus Control-Control
      par::segment_test_ignore_discordant = true;
    }

  if (a.find("--segment-test-fisher"))
    {
      par::segment_test_fisher = true;
    }

  if (a.find("--segment-test-2sided"))
    par::segment_test_1sided = false;

  if (a.find("--segment-group")) 
    {
      // Use HBD segment match routine, but without
      // allelic identity option
      par::segment_overlap = true;
      par::homo_summary_allelic_match = false;
      par::fuzzy_homo = 0.95; 
    }

   if (a.find("--segment-spanning")) 
     {
       
//        if (!a.find("--segment-group"))
// 	 error("Must specify --segment-group for --segment-spanning\n");
       
       par::segment_overlap = true;
       par::homo_summary_allelic_match = false;
       par::fuzzy_homo = 0.95;        
       par::segment_m1 = par::segment_m2 = a.value("--segment-spanning");
     }

   if (a.find("--segment-from")) 
     {
       if (!a.find("--segment-group"))
	 error("Must specify --segment-group for --segment-from\n");
       par::segment_m1 = a.value("--segment-from");
     }

   if (a.find("--segment-to")) 
     {
       if (!a.find("--segment-group"))
	 error("Must specify --segment-group for --segment-to\n");
       par::segment_m2 = a.value("--segment-to");
     }

   if (a.find("--segment-force"))
     {
       if (! ( a.find("--segment-from") && a.find("--segment-to") ) )
	 error("Can only use --segment-force with --segment-from/to \n");
       par::force_span = true;
     }

  if (a.find("--segment-match")) 
    {
      // Use HBD segment match routine, but without
      // allelic identity option
      par::segment_overlap = true;
      par::homo_summary_allelic_match = true;
      par::fuzzy_homo = a.value_double("--segment-match"); 
    }

  if (a.find("--segment-verbose"))
    par::segment_verbose = true;

  if (a.find("--read-segment"))
    {
      if (par::segment_output)
	error("Cannot specify both --segment and --read-segment\n");
      par::read_segment_filename = a.value("--read-segment");
      par::read_segment_file = true;
      par::plink = true;
      par::segment_output = true;
    }

  if (a.find("--read-segment-minimal"))
    {
      if (par::segment_output) 
	error("Cannot specify both --segment and --read-segment\n");
      par::read_segment_filename = a.value("--read-segment-minimal");
      par::read_segment_file = true;
      par::plink = true;
      par::segment_output = true;
      par::segment_minimal = true;
    }

  if (a.find("--segment-gap")) 
    {
      par::segment_output = true;
      par::segment_inter_snp_distance = a.value_int("--segment-gap");
    }

  if (a.find("--segment-length")) 
    {
      // Min length in kb, convert to bp
      par::segment_length = a.value_int("--segment-length") * 1000;
    }
  
  if (a.find("--segment-snp")) 
    {
      // Min length in SNPs
      par::segment_snp = a.value_int("--segment-snp");
    }

  if (a.find("--segment-thresholds"))
    {
      par::segment_output = true;
      vector<string> s = a.value("--segment-thresholds",2);
      par::segment_threshold_start = getDouble(s[0], "--segment-thresholds");
      par::segment_threshold_finish = getDouble(s[1], "--segment-thresholds");
    }


  //////////////////////////
  // Misc, external functions

  if (a.find("--elf-test"))
    {
      vector<string> s = a.value("--elf-test",2);
      par::rarer_maf_threshold = getDouble(s[0], "--elf-test");

      if ( par::rarer_maf_threshold < 0  || 
	   par::rarer_maf_threshold > 1 ) 
	error("Frequency thresholds not valid for --elf-test");

      par::rarer_dist_threshold = getDouble(s[1], "--elf-test") * 1000;
      par::rare_test = true;

      if (a.find("--elf-weight"))
	par::rare_test_weight1 = true;

      if (a.find("--elf-interval"))
	par::rarer_interval = a.value_int("--elf-interval");

      if (a.find("--elf-eigen"))
	par::elf_pcmode = true;

      if (a.find("--elf-2sided"))
	par::elf_pcmode_2sided = true;

    }

  if (a.find("--elf-baseline"))
    {
      par::rarer_maf_threshold = a.value_double("--elf-baseline");
      par::rare_test = true;
      par::elf_baseline = true;
    }

  if (a.find("--elf-summary"))
    {
      par::rare_test_score_range = true;

      vector<string> s = a.value("--elf-summary",2);
      par::rare_test_score_range_threshold = getDouble(s[0], "--elf-summary");
      par::rare_test_score_results_file = s[1];
      if (a.find("--elf-range"))
	{
	  par::rare_test_score_range_file = a.value("--elf-range");
	}
      else
	error("Need to specify a --elf-range file also");
      if ( a.find("--elf-controls") )
	par::rare_test_summary_controls = true;
    }
  
  if (a.find("--elf-details"))
    {
      par::rare_test_print_details = true;
      par::rare_test_print_details_snp = a.value("--elf-details");
    }

  //////////////////////////
  // Association testing
  
  if (a.find("--assoc")) 
    {
      par::assoc_test = true; 
      if ( a.find("--counts"))
	par::assoc_counts = true;
    }

  if (a.find("--qt-means")) 
    {
      if ( ! a.find("--assoc"))
	error("Can only specify --qt-means with --assoc for a QTL test");
      par::qt_means = true;
    }

  if (a.find("--logistic")) par::assoc_test = par::assoc_glm = true;
  
  if (a.find("--beta")) par::return_beta = true;

  if (a.find("--linear")) par::assoc_test = par::assoc_glm = true;

  if (a.find("--standard-beta"))
    {
      if ( ! a.find("--linear") )
	error("Must specify --linear with --standard-beta");
      par::standard_beta = true;
    }

  if (a.find("--no-snp")) par::assoc_glm_without_main_snp = true;

  if (a.find("--vif")) par::vif_threshold = a.value_double("--vif");
  
  if (a.find("--genotypic")) 
    {
      if ( a.find("--hethom"))
	par::twoDFmodel_hethom = true;

      par::twoDFmodel = true;
    }

  if (a.find("--interaction")) par::simple_interaction = true;

  if (a.find("--parameters"))
    {
      string ilist = a.value("--parameters");
            
      // NList nl(0);
      // par::parameter_list = nl.deparseNumberList(ilist);
      
      par::parameter_list = parse2int(ilist);
      par::glm_user_parameters = true;
    }

  if (a.find("--tests"))
    {
      string ilist = a.value("--tests");
      // NList nl(0);
      // par::test_list = nl.deparseNumberList(ilist);

      par::test_list = parse2int(ilist);
      par::glm_user_test = true;
    }

  if (a.find("--condition"))
    {
      if (! ( a.find("--linear") || 
	      a.find("--logistic") || 
	      a.find("--proxy-glm") || 
	      a.find("--hap-linear") ||
	      a.find("--hap-logistic") ||
	      a.find("--elf-test") ||
	      a.find("--chap")) )
	error("Can only use --condition with --linear, --logistic, --proxy-glm, --chap, --elf-test, --hap-linear or --hap-logistic");
      
      par::conditioning_snp_single = true;
      par::conditioning_snps = true;
      par::conditioning_snp_name = a.value("--condition");
    }

  if (a.find("--condition-list"))
    {
      if ( ! ( a.find("--linear") || a.find("--logistic") || a.find("--hap-linear") || a.find("--hap-logistic") || a.find("--chap")) )
	error("Can only use --condition-list with --linear, --logistic, --hap-linear, --hap-logistic or --chap");
      
      par::conditioning_snps = true;
      par::conditioning_snps_file = a.value("--condition-list");
    }

  if (a.find("--test-all")) par::test_full_model = true;

  if (a.find("--sex"))
    par::glm_sex_effect = true;

  if (a.find("--no-x-sex"))
    par::glm_no_auto_sex_effect = true;
  
  if (a.find("--xchr-model"))
    {
      par::xchr_model = a.value_int("--xchr-model");
      if (par::xchr_model < 0 || par::xchr_model >4)
	error("--xchr-model must have a value between 1 and 4");
    }

  if (a.find("--dominant"))
    {
      if ( ! ( a.find("--linear") || a.find("--logistic")) )
	error("Can only use --condition-list with --linear or --logistic");
      if (a.find("--genotypic"))
	error("Cannot specify --dominant and --genotypic together");
      par::glm_dominant = true;
      
      // Only consider autosomes
      par::xchr_model = 0; 
    }

  if (a.find("--recessive"))
    {
      if ( ! ( a.find("--linear") || a.find("--logistic")) )
	error("Can only use --condition-list with --linear or --logistic");
      if (a.find("--genotypic"))
	error("Cannot specify --recessive and --genotypic together");
      if (a.find("--dominant"))
	error("Cannot specify --recessive and --dominant together");
      par::glm_recessive = true;
      
      // Only consider autosomes
      par::xchr_model = 0; 
    }

  if (a.find("--qfam")) 
    {

      if (a.find("--within"))
	error("Cannot specify --within and --qfam together");
      if (a.find("--family"))
	error("Cannot specify --family and --qfam together");	
      if (a.find("--genedrop"))
	error("Cannot specify --genedrop and --qfam together");

      par::QTDT_test = true;
      par::assoc_test = true;
      par::QFAM_within1 = true;
    } 

  if (a.find("--qfam-total")) 
    {
   
      if (a.find("--within"))
	error("Cannot specify --within and --qfam-total together");
      if (a.find("--family"))
	error("Cannot specify --family and --qfam-total together");	
      if (a.find("--genedrop"))
	error("Cannot specify --genedrop and --qfam-total together");

      par::QTDT_test = true;
      par::assoc_test = true; 
      par::QFAM_total = true;
    }

  if (a.find("--qfam-between")) 
    {
   
      if (a.find("--within"))
	error("Cannot specify --within and --qfam-between together");
      if (a.find("--family"))
	error("Cannot specify --family and --qfam-between together");	
      if (a.find("--genedrop"))
	error("Cannot specify --genedrop and --qfam-between together");

      par::QTDT_test = true;
      par::assoc_test = true; 
      par::QFAM_between = true;
    }

  if (a.find("--qfam-parents")) 
    {
 
      if (a.find("--within"))
	error("Cannot specify --within and --qfam-parents together");
      if (a.find("--family"))
	error("Cannot specify --family and --qfam-parents together");	
      if (a.find("--genedrop"))
	error("Cannot specify --genedrop and --qfam-parents together");
      
      par::QTDT_test = true;
      par::assoc_test = true; 
      par::QFAM_within2 = true;
    }

      
  if (a.find("--mh") || a.find("--cmh") || a.find("--mh1")) par::CMH_test_1 = par::assoc_test = true;
  if (a.find("--mh-ord")) par::CMH_test_ORD = par::assoc_test = true;
  if (a.find("--mh2")) par::CMH_test_2 = par::assoc_test = true;

  if (a.find("--bd")) 
    {
      par::CMH_test_1 = par::assoc_test = true;
      par::breslowday = true;
    }

  if (a.find("--homog")) par::assoc_test = par::OR_homog_test = true;

  if (a.find("--model") || 
      a.find("--trend") || 
      a.find("--model-dom") ||
      a.find("--model-rec") ||
      a.find("--model-trend") ||
      a.find("--model-gen") )    
    { 
      par::assoc_test = true; 
      par::full_model_assoc = true; 
      
      if ( a.find("--trend") )
	{
	  par::trend_only = true;
	  par::model_perm_trend = true;
	}
      else
	{
	  if (a.find("--model-gen")) par::model_perm_gen = true;
	  else if (a.find("--model-dom")) par::model_perm_dom = true;
	  else if (a.find("--model-rec")) par::model_perm_rec = true;
	  else if (a.find("--model-trend")) par::model_perm_trend = true;
	  else par::model_perm_best = true;
	  
	  if ( (par::model_perm_gen || par::model_perm_best ) 
	       && a.find("--adjust") )
	    error("Cannot specific --model-gen or --model-best with --adjust\n"
		  "Add --model-dom, --model-rec or --model-trend instead");
	}

    }



  if (a.find("--fisher"))
    {
      par::assoc_test = true; 
      par::fisher_test = true;
      par::min_geno_cell = 0;
    }

  if (a.find("--cell"))
    par::min_geno_cell = a.value_int("--cell");

  if (a.find("--gxe")) 
    { 
      if (!a.find("--covar")) 
	error ("--covar {filename} must be specified with --gxe");
      par::assoc_gxe = true; 
    } 
  
  if (a.find("--tdt")) 
    {
      par::TDT_test = true;
      par::perm_TDT_basic = true;
      par::perm_TDT_parent = false;
    }

  if (a.find("--mating"))
    {
      par::TDT_test = true;
      par::mating_tests = true;
    }

  if (a.find("--dfam"))
    {
      if ( a.find("--family") )
	error("Cannot specify --family with --dfam\n");

      par::TDT_test = true;
      par::sibTDT_test = true;
    }

  if (a.find("--dfam-no-tdt"))
    par::dfam_tdt = false;

  if (a.find("--dfam-no-sibs"))
    par::dfam_sibs = false;
  
  if (a.find("--dfam-no-unrelateds"))
    par::dfam_unrelateds = false;
  
  if (a.find("--parentdt1")) 
    {
      par::TDT_test = true;
      par::perm_TDT_basic = false;
      par::perm_TDT_parent = true;
    }

  if (a.find("--parentdt2")) 
    {
      par::TDT_test = true;
      par::perm_TDT_basic = false;
      par::perm_TDT_parent = false;
    }

  if (a.find("--poo"))
    {
      if (!a.find("--tdt"))
	error("Parent-of-origin analysis requires --tdt option");
      
      if (a.find("--hap-tdt"))
	error("Parent-of-origin analysis not yet implemented for haplotypic TDT");
      
      par::TDT_test = true;
      par::parent_of_origin = true;

      // flavour of permutation?
      if (a.find("--pat"))
	{
	  par::perm_POO_poo = false;
	  par::perm_POO_pat = true;
	}
      else if (a.find("--mat"))
	{
	  par::perm_POO_poo = false;
	  par::perm_POO_mat = true;
	}
      else if (a.find("--best"))
	{
	  par::perm_POO_poo = false;
	  par::perm_POO_best = true;
	}
    }

  

      
  if (a.find("--sharing")) par::ibs_sharing_test = true;
  
  if (a.find("--boot")) par::boot = true;
  
  
  if (a.find("--blocks"))
    {
      par::make_blocks = true;
      // This default can be over-ridden below (200kb)
      par::disp_r_window_kb = 200 * 1000;      
    }

  if (a.find("--ld"))
    {
      par::calc_SNPSNP_LD = true;
      vector<string> s = a.value("--ld",2);
      par::ld_SNP1 = s[0];
      par::ld_SNP2 = s[1];
    }
  
  
  if (a.find("--ld-snp"))
    {
      par::disp_r2 = true;
      par::ld_anchor = true;
      par::ld_SNP1 = a.value("--ld-snp");
    }
  else if (a.find("--ld-snp-list"))
    {
      par::disp_r2 = true;
      par::ld_anchor = true;
      par::ld_anchor_list = true;
      par::ld_SNP1_file = a.value("--ld-snp-list");
      checkFileExists(par::ld_SNP1_file);
    }

  if (a.find("--ld-window"))
    {
      par::disp_r2 = true;
      if (a.find("--matrix"))
	error("Cannot specify --matrix and --ld-window together\n");
      par::disp_r_window = true;
      par::disp_r_window_snp = a.value_int("--ld-window");      
    }

  if (a.find("--ld-window-kb"))
    {
      par::disp_r2 = true;
      if (a.find("--matrix"))
	error("Cannot specify --matrix and --ld-window together\n");
      par::disp_r_window = true;
      // Store in base-pair units
      par::disp_r_window_kb = a.value_int("--ld-window-kb") * 1000;      
    }


  if (a.find("--ld-window-r2"))
    {
      par::disp_r2 = true;
      if (a.find("--matrix"))
	error("Cannot specify --matrix and --ld-window together\n");
      par::disp_r_window_r2 = a.value_double("--ld-window-r2");      
    }

  if (a.find("--r2")) 
    {
      par::disp_r2 = true;
      if ( ! ( a.find("--matrix") || a.find("--inter-chr") ) )
	par::disp_r_window = true;
    }

  // By default, we assume r^2 is often interest, but here at the end, if 
  // specified we can swap to make the basic correlation
  if (a.find("--r")) 
    {
      par::disp_r1 = true;
      par::disp_r2 = false;
      if ( ! ( a.find("--matrix") || a.find("--inter-chr") ) )
	par::disp_r_window = true;
    }


  if (a.find("--flip-scan"))
    {
      par::flip_scan = true;
      if ( a.find("--flip-scan-threshold") )
	par::flip_scan_threshold = a.value_double("--flip-scan-threshold");
      if ( a.find("--flip-scan-verbose") )
	par::flip_scan_verbose = true;
      
    }
  
  if (a.find("--indep"))
    {
      if (makedata)
	{
	  string msg = "Cannot specify --indep with --make-bed or --recode\n";
	  msg += "    use --extract/--exclude with *.prune.in/*.prune.out files";
 	  error(msg);
	}
      
      par::prune_ld = true;
      vector<string> s = a.value("--indep",3);
      par::prune_ld_win = getInt(s[0].c_str(),"--indep");
      par::prune_ld_step = getInt(s[1].c_str(),"--indep");
      par::prune_ld_vif = getDouble(s[2].c_str(),"--indep");
      
      if (par::prune_ld_win<2)
	error("Cannot have a window size < 2 for --indep {window} {step} {VIF}");

      if (par::prune_ld_step<1)
	error("Cannot have a window step < 1 for --indep {window} {step} {VIF}");

      if (par::prune_ld_vif<1)
	error("Cannot have a VIF threshold < 1 for --indep {window} {step} {VIF}");

    }

  if (a.find("--indep-pairwise"))
    {
      if (makedata) 
        {
          string msg = "Cannot specify --indep with --make-bed or --recode\n";
          msg += "    use --extract/--exclude with *.prune.in/*.prune.out files";
          error(msg);
        }

      par::prune_ld = true;
      par::prune_ld_pairwise = true;
      vector<string> s = a.value("--indep-pairwise",3);
      par::prune_ld_win = getInt(s[0].c_str(),"--indep-pairwise");
      par::prune_ld_step = getInt(s[1].c_str(),"--indep-pairwise");
      par::prune_ld_r2 = getDouble(s[2].c_str(),"--indep-pairwise");
      
      if (par::prune_ld_win<2)
        error("Cannot have a window size < 2 for --indep-pairwise {window} {step} {R^2}");
      if (par::prune_ld_step<1)
        error("Cannot have a window step < 1 for --indep-pairwise {window} {step} {R^2}");
      if (par::prune_ld_r2<0)
        error("Cannot have an R^2 threshold < 0 for --indep-pairwise {window} {step} {R^2}");
      if (par::prune_ld_r2>1)
        error("Cannot have an R^2 threshold > 1 for --indep-pairwise {window} {step} {R^2}");

      // Compare to correlation
      par::prune_ld_r2 = sqrt(par::prune_ld_r2);

      if (a.find("--indep-prefer"))
	{
	  par::prune_r2_prefer = true;
	  par::prune_r2_prefer_list = a.value("--indep-prefer");
	}

      if (a.find("--indep-fixed"))
	{
	  par::prune_r2_fixed = true;
	  par::prune_r2_fixed_list = a.value("--indep-fixed");
	}


    }
  
  if (a.find("--T2") || a.find("--t2") ) 
    {

      error("This command is disabled in v1.04");

      if ( a.find("--set-test") || a.find("--assoc") )
	error("Cannot specify --T2 and other association commands");
      par::set_test = true;
      par::hotel = true;
    }

  if (a.find("--set")) 
    {
      par::read_set = true;
      par::setfile = a.value("--set");
    }
  
  if (a.find("--set-test"))
    {
      par::set_test = true;
      
      // Force use of LD-aware test
      par::set_r2 = true;
      
      if ( (!a.find("--gene")) && ( a.find("--assoc") || a.find("--tdt") ) && 
	   (  ! ( a.find("--mperm") || a.find("--perm") ) ) )
	error("Must use --mperm N or --perm with set association tests");

    }
  
  if (a.find("--set-screen"))
    {
      if ( ! ( a.find("--set") || a.find("--make-set") ) )
	error("Must specify a --set or --make-set\n");
      par::set_screen = true;
      par::set_screen_resultfile = a.value("--set-screen");
    }


  if (a.find("--set-step"))
    {
      par::set_test = true;
      par::set_step = true;
      par::set_step_in = a.value_double("--set-step");
    }

  if (a.find("--set-score"))
    {
      par::set_test = true;
      par::set_score = true;
      par::set_score_p = a.value_double("--set-score");
      par::set_max = par::set_min = 1;
      if ( a.find("--set-min") || a.find("--set-max") )
	error("Cannot specify --set-min or --set-max with --set-score");
    }
  
  if (a.find("--set-p2"))
    {
      par::set_p2 = true;
    }

  if (a.find("--set-min"))
    {
      if (!a.find("--set"))
	error ("You need to specify --set also");
      par::set_min = a.value_int("--set-min");
    }
  
  if (a.find("--set-max"))
    {
      if (! ( a.find("--set") || a.find("--make-set")))
	error ("You need to specify --set or --make-set also");
      if (!a.find("--set-test"))
	error ("You need to specify --set-test also");
      par::set_max = a.value_int("--set-max");
    }

  if (a.find("--set-r2"))
    {
      if (! ( a.find("--set") || a.find("--make-set")))
	error ("You need to specify --set or --make-set also");
      if (!a.find("--set-test"))
	error ("You need to specify --set-test also");
      if (a.find("--set-r2-phase"))
	par::set_r2_phase = true;
      par::set_r2_val = a.value_double("--set-r2");
      par::set_r2 = true;
    }

  if (a.find("--write-set-r2"))
    {
      par::set_r2_write = true;
    }

  if (a.find("--read-set-r2"))
    {
      par::set_r2_read = true;
      par::set_r2_read_file = a.value("--read-set-r2");
    }
  
  if (a.find("--set-p"))
    {
      if (!a.find("--set-test"))
	error ("You need to specify --set-test also");

      double p = a.value_double("--set-p");
      if ( p <= 0 || p > 1 ) 
	error("P-value for --set-p must be 0<p<=1");
      par::set_chisq_threshold = inverse_chiprob(p,1);
    }

  if (a.find("--subset"))
    {
      if ( ! ( a.find("--set") || a.find("--make-set")))
	error("Must also specify --set or --make-set");
      par::subsetfile = a.value("--subset");
      par::use_subset = true;
    }

  if (a.find("--set-table"))
    {
      if ( ! ( a.find("--set") || a.find("--make-set") ))
	error("Must also specify a --set {file} or --make-set {file}");
      par::set_table = true;
      par::drop_sets = false;
    }
  

  //////////////////////////
  // Impute tagged SNPs
  
  if (a.find("--hap-impute")) 
    { 
      par::impute_tags = true;
      par::phase_snps = true;
      //      par::hap_missing_geno = 1;
    }

  if (a.find("--hap-impute-verbose"))
    {
      par::impute_tags = true;    
      par::phase_snps = true;

      //      par::hap_missing_geno = 1;
      par::impute_verbose = true;
    }

  if (a.find("--hap")) 
    { 
      par::phase_snps = true; 
      par::tagfile = a.value("--hap"); 
    }

  if (a.find("--hap-window"))
    {
      par::phase_snps = true;
      par::sliding_window = true;
      par::phase_hap_all = true; 
      if (a.find("--hap"))
	error("Cannot specify both --hap-window and --hap {file}\n");
      par::sliding_window_size = a.value("--hap-window");            
    }


  // Using all SNPs, perform a basic sliding window analysis and look
  // for haplotypically shared regions, i.e. given that we now have
  // much more highly heterozygous markers

  if (a.find("--homozyg-haplo-track"))
    {
      if ( ! a.find("--hap-window"))
	error("The 'haplo-track' option requires --hap-window {s} to be specified");
      
      par::segment_haplotrack = true;
      vector<string> s = a.value("--homozyg-haplo-track",2);
      par::segment_haplotrack_fid1 = s[0];
      par::segment_haplotrack_iid1 = s[1];
      par::segment_haplotrack_fid2 = par::segment_haplotrack_fid1;
      par::segment_haplotrack_iid2 = par::segment_haplotrack_iid1;      
    }

  if (a.find("--segment-haplo-track"))
    {
      if ( ! a.find("--hap-window"))
	error("The 'haplo-track' option requires --hap-window {s} to be specified");

      par::segment_haplotrack = true;
      vector<string> s = a.value("--segment-haplo-track",4);
      par::segment_haplotrack_fid1 = s[0];
      par::segment_haplotrack_iid1 = s[1];
      par::segment_haplotrack_fid2 = s[2];
      par::segment_haplotrack_iid2 = s[3];
    }



  /////////////////////////////////
  // EM Phasing options 

  if (a.find("--em-verbose"))
    par::haplo_plem_verbose = true;

  if (a.find("--em-follow"))
    {
      par::haplo_plem_follow = true;
      vector<string> s = a.value("--em-follow",2);
      par::haplo_plem_follow_fid = s[0];
      par::haplo_plem_follow_iid = s[1];
    }

  if (a.find("--em-window"))
    par::haplo_plem_window = a.value_int("--em-window");

  if (a.find("--em-overlap"))
    par::haplo_plem_original_overlap = par::haplo_plem_overlap = a.value_int("--em-overlap");

  if (a.find("--em-window-iter"))
    par::haplo_plem_iter = a.value_int("--em-window-iter");  

  if (a.find("--em-window-prune-phase"))
    par::haplo_plem_window_prune_phase = 
      a.value_double("--em-window-prune-phase");

  if (a.find("--em-window-likelihood"))
    par::haplo_plem_likelihood_iter = a.value_int("--em-window-likilood");
    
  if (a.find("--em-window-tol"))
    par::haplo_plem_window_tol = a.value_double("--em-window-tol");
  
  if (a.find("--em-window-prune-haplotype"))
    {
      par::haplo_plem_zero_threshold = a.value_double("--em-window-prune-haplotype");
      if ( par::haplo_plem_zero_threshold == 0 )
	par::haplo_plem_nonzero_threshold = false;
      else
	par::haplo_plem_nonzero_threshold = true;
    }
  


  //////////////////////////
  // Meta EM parameters

  if (a.find("--em-meta-window"))
    {
      par::haplo_plem_meta_window = a.value_int("--em-meta-window");
      if ( par::haplo_plem_meta_window < 2 )
	error("--em-meta-window must be >1");
    }
  
  if (a.find("--em-meta-prune-haplotype"))
    par::haplo_plem_meta_prune_haplotype = 
      a.value_double("--em-meta-prune-haplotype");

  if (a.find("--em-meta-prune-phase"))
    par::haplo_plem_meta_prune_phase = 
      a.value_double("--em-meta-prune-phase");

  if (a.find("--em-meta-iter"))
    par::haplo_plem_meta_iter = a.value_int("--em-meta-iter");  

  if (a.find("--em-meta-likilood"))
    par::haplo_plem_meta_likelihood_iter = a.value_int("--em-meta-likilood");

  if (a.find("--em-meta-tol"))
    par::haplo_plem_meta_tol = a.value_double("--em-meta-tol");
  

  /////////////////////////////
  // Other haplotype options

  if (a.find("--hap-all")) 
    { 
      if (!a.find("--hap"))
	error("--hap-all modifies --hap, but you have not specified --hap");
      par::phase_hap_all = true; 
    }
  
  if (a.find("--whap")) 
    { 
      par::phase_snps = true; 
      par::tagfile = a.value("--whap"); 
      par::weighted_mm = true;
    }

  if (a.find("--hap-pp")) { par::hap_post_prob = a.value_double("--hap-pp"); }
  if (a.find("--hap-miss")) { par::hap_missing_geno = a.value_double("--hap-miss"); }
  
  if (a.find("--hap-freq")) 
    { 
      if ( par::impute_tags ) 
	error("Cannot specify --hap-impute and --hap-freq\n");
      par::display_hap_freqs = true;     
    }

  if (a.find("--hap-assoc")) 
    { 
      if ( par::impute_tags ) 
	error("Cannot specify --hap-impute and --hap-assoc\n");
      par::test_hap_CC = true; 
    }

  if ( a.find("--hap-logistic") || a.find("--hap-linear") )
    {
      par::test_hap_GLM = true;
      if ( a.find("--hap-omnibus"))
	par::test_hap_GLM_omnibus = true;
      if ( a.find("--perm"))
	error("Cannot currently use --perm and --hap-logistic or --hap-linear. Use --mperm N instead");
    }

  if (a.find("--chap"))
    {
      if ( ! a.find("--hap-snps"))
	error("--chap requires --hap-snps");

      par::chap_test = true;

      if ( a.find("--null-group") )
	{
	  par::chap_specified_groups = true;
	  par::chap_model0 = a.value("--null-group");
	}
      
      if ( a.find("--alt-group") )
	{
	  par::chap_specified_groups = true;
	  par::chap_model1 = a.value("--alt-group");
	}

      if ( a.find("--null-snp") )
	{
	  if ( par::chap_specified_groups )
	    error("Cannot specify SNPs and groups for --chap tests");
	  par::chap_specified_snps = true;
	  par::chap_model0 = a.value("--null-snp");
	}

      if ( a.find("--alt-snp") )
	{
	  if ( par::chap_specified_groups )
	    error("Cannot specify SNPs and groups for --chap tests");
	  par::chap_specified_snps = true;
	  par::chap_model1 = a.value("--alt-snp");
	}

      if ( a.find("--control") )
	{
	  if ( par::chap_specified_groups ||
	       par::chap_specified_snps )
	    error("Cannot use --control and other --chap options");
	  par::chap_sole_variant = true;
	  par::chap_entity = a.value("--control");
	  
	  if ( a.find("--control-alleles"))
	    {
	      par::chap_sole_variant_specific_alleles = true;
	      par::chap_sole_variant_specific_allele_list = a.value("--control-alleles");
	    }
	}
      else if ( a.find("--independent-effect") )
	{
	  if ( par::chap_specified_groups ||
	       par::chap_specified_snps )
	    error("Cannot use --independent-effect and other --chap options");
	  par::chap_independent_effect = true;
	  par::chap_entity = a.value("--independent-effect");
	}      
      else if ( a.find("--specific-haplotype") )
	{
	  if ( par::chap_specified_groups ||
	       par::chap_specified_snps )
	    error("Cannot use --specific-haplotype and other --chap options");
	  par::chap_haplotype_specific = true;
	  par::chap_entity = a.value("--specific-haplotype");	  
	}

      if ( a.find("--test-snp") )
	{
	  if ( ! ( a.find("--condition") || a.find("--condition-list") ) )
	    error("You must first specify conditioning SNPs");
	  par::chap_drop_snps = true;
	  par::chap_drop_snps_list = a.value("--test-snp");
	}

      if ( a.find("--each-vs-others") )
	par::chap_add_grp_specifics = true;

      if ( a.find("--each-versus-others") )
	par::chap_add_grp_specifics = true;

    }

  
  if (a.find("--hap-snps"))
    {
      if ( a.find("--hap") || a.find("--hap-window"))
	error("Cannot specify --hap-snps with --hap or --hap-window");
      par::phase_snps = true;
      par::hap_specific_snps = true;
      par::hap_specific_snps_list = a.value("--hap-snps");
      par::phase_hap_all = true; 
    }
  
  if (a.find("--hap-tdt")) 
    { 
      if ( par::impute_tags ) 
	error("Cannot specify --hap-impute and --hap-tdt\n");
      par::test_hap_TDT = true; 
    }
  
  if (a.find("--hap-phase")) 
    { 
      if ( par::impute_tags ) 
	error("Cannot specify --hap-impute and --hap-phase\n");
      par::display_phase_probs = true;
    }
  
  if (a.find("--hap-phase-wide")) 
    { 
      if ( par::impute_tags ) 
	error("Cannot specify --hap-impute and --hap-phase-wide\n");
      par::display_phase_probs = par::display_phase_probs_wide = true; 
    }
  
  if (a.find("--hap-only")) { par::test_hap_only = true; } 
  
  if (a.find("--hap-max-phase")) 
    { 
      par::hap_max_nf_phases = a.value_int("--hap-max-phase"); 
      if ( par::hap_max_nf_phases < 1 ) error("Invalid --hap-max-phase value (should be >0)");
    }
  
  if (a.find("--hap-min-phase-prob")) 
    { 
      par::hap_min_phase_prob = a.value_double("--hap-min-phase-prob"); 
      if ( par::hap_min_phase_prob < 0 || par::hap_min_phase_prob > 1 ) 
	error("Invalid --hap-min-phase-prob value, should be between 0 and 1");
    }
  



  //////////////////////////
  // Epistasis

  if (a.find("--epistasis")) 
    {
      if (a.find("--set") || a.find("--make-set"))
	par::set_test = true;
      par::epistasis = true;
    }

  // Use odds-ratio test as default fast-epistasis method
  if (a.find("--fast-epistasis")) 
    { 
      par::fast_epistasis = par::epistasis = true; 
      if ( a.find("--set") || a.find("--make-set"))
	par::set_test = true;
    }

  if (a.find("--case-only")) 
    { 
      par::epistasis = true; 
      par::fast_epistasis = true;
      par::epi_caseonly = true; 
      if (a.find("--epistasis"))
	error("--case-only requires --fast-epistasis");
    } 

  if (a.find("--gap")) 
    { 
      if (!a.find("--case-only")) error("--gap option only valid when --caseonly is in effect");
      par::epi_caseonly_kb_gap = a.value_double("--gap"); 
    }
  if (a.find("--nop")) par::epi_quickscan = true;
  if (a.find("--set-by-all")) par::set_by_set = par::drop_sets = false;
  if (a.find("--epi1")) par::epi_alpha1 = a.value_double("--epi1");
  if (par::epi_alpha1 > 1) par::epi_filter = false;
  if (a.find("--epi2")) par::epi_alpha2 = a.value_double("--epi2");
  if (a.find("--twolocus"))
    {
      par::list_twolocus = true;
      vector<string> s = a.value("--twolocus",2);
      par::twolocus_snp1 = s[0];
      par::twolocus_snp2 = s[1];      
    }
  
  if (a.find("--genepi")) 
    {
      par::set_test = true;
      par::epi_genebased = true;
    }


  
  ///////////////////////////////////////
  // Gene X environment / heterogeneity
  



  //////////////////////////////////////
  // File output options

  if (a.find("--freq")) 
    {
      par::af_write = true;
      // and unless otherwise specified, set GENO = 1 and MAF = 0
//       if (!a.find("--maf")) par::min_af = 0.0;
//       if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;
//       if (!a.find("--mind")) par::MAX_IND_MISSING = 1;

      // display MAF counts instead of freqs?
      if (a.find("--counts"))
	{
	  if (a.find("--within"))
	    error("Cannot specify --counts and --within together\n");
	  par::af_count = true;
	}
   }  


  
  if (a.find("--nonfounders")) par::summ_nonfounders = true;

  if (a.find("--make-founders")) par::make_founders = true;

  if (a.find("--allow-no-sex")) par::ignore_missing_sex = true;
  
  if (a.find("--must-have-sex"))
    {
      if ( ! makedata )
	error("Can only specify --must-have-sex with a data generation command");
      par::ignore_missing_sex = false;
    }
  else if ( makedata ) 
    {
      par::ignore_missing_sex = true;
    }


  if (a.find("--read-freq")) 
    {
      par::af_read = true;
      par::af_file = a.value("--read-freq");
    }

  if (a.find("--read-genome")) 
    {
      par::ibd_read = true;
      par::ibd_file = a.value("--read-genome");
      checkFileExists( par::ibd_file );
    }
  else if (a.find("--read-genome-list"))
    {
      par::ibd_read = par::ibd_read_list = true;
      par::ibd_file_list = a.value("--read-genome-list");
      checkFileExists( par::ibd_file_list ); 
    }

  if (a.find("--Z-genome"))
    par::compress_genome = true;

  if ( a.find("--genome-groups") )
    {
      if ( ! a.find("--genome-groups"))
	error("You must specify a .genome file with --read-genome");
      if ( ! a.find("--within"))
	error("You must specify a cluster file with --within");
      par::genome_groups = true;
    }
  
  if (a.find("--read-genome-minimal"))
    {
      if (a.find("--read-genome")) 
	error("Cannot specify both --read-genome and --read-genome-minimal");
      par::ibd_read = true;
      par::ibd_read_minimal = true;
      par::ibd_file = a.value("--read-genome-minimal");
      checkFileExists( par::ibd_file );
    }
  

  if (a.find("--list")) 
    {
      par::list_by_allele = true;
      // and unless otherwise specified, set GENO = 1 and MAF = 0
      if (!a.find("--maf")) par::min_af = 0.0;
      if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
      if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
    }
  
  if (a.find("--report"))
    {
      par::indiv_report = true;
      vector<string> s = a.value("--report",2);
      par::indiv_report_fid = s[0];
      par::indiv_report_iid = s[1];
      if (!a.find("--maf")) par::min_af = 0.0;
      if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
      if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
    }

  if (a.find("--plist")) 
    {
      par::plist = true;
      // and unless otherwise specified, set GENO = 1 and MAF = 0
      if (!a.find("--maf")) par::min_af = 0.0;
      if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;     
      if (!a.find("--mind")) par::MAX_IND_MISSING = 1;

      vector<string> s = a.value("--plist",4);
      par::plist_fid1 = s[0];
      par::plist_iid1 = s[1];
      par::plist_fid2 = s[2];
      par::plist_iid2 = s[3];
      
    }


  if (a.find("--fix-allele"))
    {
      if (! (a.find("--recodeAD") || a.find("--recodeA")) )
	error("--fix-allele option only works with --recodeA/--recodeAD options");
      else par::recode_AD_fixed = true;
    }

  if (a.find("--merge"))
    {
      if (a.find("--merge-list") || a.find("--bmerge") )
	error("Can only specify --merge or --bmerge or --merge-list");

      par::merge_data = true;
      vector<string> s = a.value("--merge",2);
      par::merge_pedfile = s[0];
      par::merge_mapfile = s[1];
      checkFileExists( par::merge_pedfile );
      checkFileExists( par::merge_mapfile );
    }

  if (a.find("--bmerge"))
    {
      if (a.find("--merge-list") || a.find("--merge") )
	error("Can only specify --bmerge or --merge or --merge-list");

      par::merge_data = true;
      par::merge_binary = true;

      vector<string> s = a.value("--bmerge",3);
      par::merge_bedfile = s[0];
      par::merge_bimfile = s[1];
      par::merge_famfile = s[2];

      checkFileExists( par::merge_bedfile );
      checkFileExists( par::merge_bimfile );
      checkFileExists( par::merge_famfile );

    }

  if (a.find("--merge-flip"))
    par::merge_force_strand = true;

  if (a.find("--merge-list"))
    {
      if (a.find("--merge") || a.find("--bmerge") )
	error("Can only specify --merge or --bmerge or --merge-list");
      
      par::merge_data = true;
      par::merge_list = true;
      par::merge_list_filename = a.value("--merge-list");
    }

  if (a.find("--merge-mode"))
    {
      if (! (a.find("--merge") || a.find("--merge-list") || a.find("--bmerge") ) )
	error("Can only specify --merge-mode when --merge or --bmerge or --merge-list is used");
      par::merge_mode = a.value_int("--merge-mode");
      if (par::merge_mode < 1 || par::merge_mode > 7)
	error("--merge-mode N, where N must be between 1 and 7"); 
      
      if (par::merge_list && par::merge_mode >= 6)
	error("Can not specify --merge-mode 6/7 (diff) and --merge-list");
    }


  if (a.find("--flip")) 
    {
      par::flip_strand = true;     
      par::flip_file = a.value("--flip");
 
      if (a.find("--flip-subset"))
	{
	  par::flip_subset = true;
	  par::flip_subset_file = a.value("--flip-subset");
	  checkFileExists( par::flip_subset_file );
	}
    }

  // Remove these individuals from a file
  if (a.find("--remove")) 
    { 
      par::remove_indiv = true; 
      par::remove_indiv_list = a.value("--remove");
      checkFileExists( par::remove_indiv_list );
    }

  // Keep only these individuals from a file 
  if (a.find("--keep")) 
    { 
      par::keep_indiv = true;
      par::keep_indiv_list = a.value("--keep");
      checkFileExists( par::keep_indiv_list );
    }
  
  // By default, remove then keep
  if (a.find("--keep-before-remove"))
    par::remove_before_keep = false;

  // Extract only these SNPs from a file
  if (a.find("--extract")) 
    { 
      par::extract_set = true; 
      par::extract_file = a.value("--extract");
      checkFileExists( par::extract_file );
    }

  // Exclude these SNPs from a file
  if (a.find("--exclude")) 
    { 
      par::exclude_set = true; 
      par::exclude_file = a.value("--exclude");
      checkFileExists( par::exclude_file );  
    }
  
  if (a.find("--thin"))
    {
      par::thin_snps = true;
      par::thin_param = a.value_double("--thin");
    }

  // Select a set of SNPs based on different physical 
  // positions: this modifies the behavior of --extract 
  // and --exclude

  if (a.find("--range"))
    {
      if ( ! ( a.find("--extract") || a.find("--exclude") ) )
	error("Must specify --extract or --exclude with --range");
      par::snp_range_list = true;
    }

  if (a.find("--border"))
    {
      par::make_set_border = a.value_int("--border") * 1000;
    }

  if (a.find("--make-set"))
    {
      par::make_set = true;
      par::make_set_file = a.value("--make-set");

      if ( a.find("--make-set-border"))
	par::make_set_border = a.value_int("--make-set-border") * 1000;      

      if ( a.find("--make-set-collapse-group"))
	{
	  par::make_set_collapse = true;
	}
      else if ( a.find("--make-set-collapse-all"))
	{
	  par::make_set_collapse = true;
	  par::make_set_collapse_label = a.value("--make-set-collapse-all");	  
	  par::make_set_ignore_group = true;
	}
      else if ( a.find("--make-set-complement-group") )
	{
	  par::make_set_complement = true;	  
	}
      else if ( a. find("--make-set-complement-all"))
	{
	  par::make_set_complement = true; 
	  par::make_set_ignore_group = true;
	  par::make_set_collapse_label = a.value("--make-set-complement-all");
	}
      
    }

  if (a.find("--write-set"))
    {
      if ( ! ( a.find("--set") || a.find("--make-set")))
	error("Must specify --set or --make-set with --write-set");
      par::write_set = true;
    }

  // By default, extract before exclude
  if (a.find("--exclude-before-extract"))
    par::extract_before_exclude = false;
  

  // Extract SNPs in a certain GENE (specified in the SET file)
  if (a.find("--gene"))
    {
      if (! (a.find("--set") || a.find("--make-set")))
	error("You must also specify --set or --make-set with --gene\n");
      par::extract_set = true;
      par::dump_gene = true;
      par::dump_genename = a.value("--gene");
    }
  

  if (a.find("--ind-major"))
    {
      if (!a.find("--make-bed"))
	error("You can only specify --ind-major when --make-bed is in effect");
      par::out_SNP_major = false;
    }

  
  if (a.find("--include")) par::inc_write = true;
  
  if (a.find("--read-include")) 
    {
      par::inc_read = true;
      par::inc_file = a.value("--read-include");
    }



  ////////////////////////////////
  // Genotype quality score files
  
  if (a.find("--qual-scores"))
    {
      par::read_snp_qual = true;
      par::snp_qual_file = a.value("--qual-scores");

      if (a.find("--qual-threshold"))
	par::snp_qual_min = a.value_double("--qual-threshold");
      
      if (a.find("--qual-max-threshold"))
	par::snp_qual_max = a.value_double("--qual-max-threshold");
      
    }
  
  if (a.find("--qual-geno-scores"))
    {
      par::read_geno_qual = true;
      par::geno_qual_file = a.value("--qual-geno-scores");

      if (a.find("--qual-geno-threshold"))
	par::geno_qual_min = a.value_double("--qual-geno-threshold");
      
      if (a.find("--qual-geno-max-threshold"))
	par::geno_qual_max = a.value_double("--qual-geno-max-threshold");
      
    }

  

  ///////////////////////////
  // Plink phenotype definition
//   // Squared differences (quantitative trait)
//   if (a.find("--SD")) { par::SD = true; par::CP = false; }
//   // Cross-product (quantitative trait)
//   if (a.find("--CP")) { par::CP = true; par::SD = false; } 
//   // Fix the prevalence of a binary trait
//   if (a.find("--prev")) {
//     par::fix_prev=true; 
//     par::fixed_prev = a.value_double("--prev"); 
//   }
//   // Remove unaffected pairs from analysis 
//   // i.e. so only discordant versus affected concordant
//   // i.e. need at least 1 affected "1A"
//   if (a.find("--1aff")) { par::remove_unaffected_pairs = true; }



  ///////////////////////////
  // Some basic filters

  if (a.find("--prune")) 
    par::ignore_phenotypes = false;  
  if (a.find("--filter-cases"))
    par::filter_cases = true;
  if (a.find("--filter-controls"))
    par::filter_controls = true;
  if (a.find("--filter-females"))
    par::filter_females = true;
  if (a.find("--filter-males"))
    par::filter_males = true;
  if (a.find("--filter-founders"))
    par::filter_founders = true;
  if (a.find("--filter-nonfounders"))
    par::filter_nonfounders = true;

  if (a.find("--filter-cases") && a.find("--filter-controls"))
    error("Cannot filter on both cases and controls");
  if (a.find("--filter-males") && a.find("--filter-females"))
    error("Cannot filter on both males and females");
  if (a.find("--filter-founders") && a.find("--filter-nonfounders"))
    error("Cannot filter on both founders and nonfounders");

  if (a.find("--attrib")) 
    {
      par::snp_attrib_filter = true;
      vector<string> s = a.value("--attrib",2);
      par::snp_attrib_file = s[0];
      par::snp_attrib_value = s[1];
    }
  
  if (a.find("--attrib-indiv")) 
    {
      par::ind_attrib_filter = true;
      vector<string> s = a.value("--attrib-indiv",2);
      par::ind_attrib_file = s[0];
      par::ind_attrib_value = s[1];
    }


  /////////////////////////////////
  // Basic input file processing 
  
  if (a.find("--dummy"))
    {
      vector<string> s = a.value("--dummy",2);      
      par::dummy = true;
      par::dummy_nind = getInt(s[0].c_str(),"--dummy");
      par::dummy_nsnp = getInt(s[1].c_str(),"--dummy");
    }  
  
  if (a.find("--simulate"))
    {
      par::simul = true;
      par::simul_file = a.value("--simulate");
    }
  else if (a.find("--simulate-qt"))
    {
      par::simul = par::simul_qt = true;
      par::simul_file = a.value("--simulate-qt");
      if (a.find("--simulate-n"))
	par::simul_ncases = a.value_int("--simulate-n");
    }

  if (a.find("--simulate-label"))
    par::simul_label = a.value("--simulate-label");
  
  if (a.find("--simulate-tags"))
    {
      if (!a.find("--simulate"))
	error("Requires --simulate for --simulate-tags");
      par::simul_tags = true;
    }
  
  if (a.find("--simulate-haps"))
    {
      if (!a.find("--simulate"))
	error("Requires --simulate for --simulate-haps");
      par::simul_tags = true;
      par::simul_haps = true;
    }


  if (a.find("--simulate-ncases"))
    par::simul_ncases = a.value_int("--simulate-ncases");

  if (a.find("--simulate-ncontrols"))
    par::simul_ncontrols = a.value_int("--simulate-ncontrols");

  if (a.find("--simulate-prevalence"))
    par::simul_prevalence = a.value_double("--simulate-prevalence");
            

  ////////////////////////////////////////////////////
  // Main file input options

  if (a.find("--compress"))
    {
      par::compress_file = true;
      par::compress_filename = a.value("--compress");
    }

  if (a.find("--decompress"))
    {
      par::uncompress_file = true;
      par::compress_filename = a.value("--decompress");
    }

  if (a.find("--file")) 
    { 
      if (a.find("--map") || a.find("--ped") ) 
	error("Use either --file {root} OR --ped {name} --map {name}");
      par::read_ped = true;
      par::fileroot = a.value("--file");
      par::pedfile = par::fileroot + ".ped";
      par::mapfile = par::fileroot + ".map";
    }


  if (a.find("--tfile")) 
    { 
      if (a.find("--tfam") || a.find("--tped") ) 
	error("Use either --tfile {root} OR --tped {name} --tfam {name}");

      par::tfile_input = true;
      par::fileroot = a.value("--tfile");
      par::tpedfile = par::fileroot + ".tped";
      par::tfamfile = par::fileroot + ".tfam";
    }

  if (a.find("--tped")) 
    {
      par::tpedfile = a.value("--tped");
      if (par::tpedfile == "-")
	  par::ped_from_stdin = true;
      par::tfile_input = true;
    }
  
  if (a.find("--tfam")) 
    {
      par::tfamfile = a.value("--tfam");
      par::tfile_input = true;
    }
      
  // Long file format

  if (a.find("--lfile")) 
    { 
      if (a.find("--fam") || a.find("--lgen") || a.find("--map") )  
	error("Use either --lfile {root} OR --lgen {name} --map (name) --fam {name}");

      par::lfile_input = true;
      par::fileroot = a.value("--lfile");
      par::lpedfile = par::fileroot + ".lgen";
      par::famfile = par::fileroot + ".fam";
      par::mapfile = par::fileroot + ".map";
    }

  if (a.find("--lgen")) 
    {
      par::lpedfile = a.value("--lgen");
      if (par::lpedfile == "-")
	  par::ped_from_stdin = true;
      par::lfile_input = true;
    }
  
  ////////////////////////
  // Reference allele file

  if (a.find("--reference"))
    {
      par::ref_file = true;
      par::ref_file_name = a.value("--reference");
    }

  /////////////////////
  // Generic variant 

  if ( a.find("--gfile") || a.find("--gvar") )
    {
      // Analyse generic variants
      par::gvar = true;
      
      // Primarily load these; this will be changed
      // if another load is not specified
      par::load_gvar = false;
      
      if ( a.find("--gfile") )
	{
	  par::fileroot = a.value("--gfile");
	  par::gvarfile = par::fileroot + ".gvar";
	  par::gmapfile = par::fileroot + ".map";
	  par::gfamfile = par::fileroot + ".fam";
	}
      else
	{
	  par::gvarfile = a.value("--gvar");
	}

      if (a.find("--gvar-verbose"))
	par::gvar_verbose_association = true;

      if (a.find("--gvar-all"))
	par::gvar_include_all_variants = true;
      
      if (a.find("--gvar-convert"))
	par::gvar_to_standard = true;

      if (a.find("--gvar-verbose"))
	par::gvar_full_report = true;      
    }
  
  if (a.find("--gvar-write"))
    par::gvar_write = true;
      
  // Text-file modifiers

  if (a.find("--map3"))
    par::map3 = true;

  if (a.find("--no-sex"))
    par::ped_skip_sex = true;

  if (a.find("--no-parents"))
    par::ped_skip_parents = true;

  if (a.find("--no-fid"))
    par::ped_skip_fid = true;
  
  if (a.find("--no-pheno"))
    par::ped_skip_pheno = true;

  if (a.find("--liability"))
    par::liability = true;
      
  if (a.find("--bfile")) 
    { 
      if (a.find("--bim") || a.find("--bed") || a.find("--fam") ) 
	error("Use either --bfile {root} OR --bed {name} --bim {name} --fam {name}");
      
      par::read_bitfile = true;
      par::fileroot = a.value("--bfile");
      par::bitfilename = par::pedfile = par::fileroot + ".bed";
      par::bitfilename_map = par::mapfile = par::fileroot + ".bim";
      par::famfile = par::fileroot + ".fam";
    }

  if (a.find("--no-snps"))
    par::do_not_load_snps = true;

  if (a.find("--bfile-faster"))
    par::fast_binary = true;
  
  if (a.find("--ped")) 
      {
	par::read_ped = true;
	par::pedfile = a.value("--ped");
	if (par::pedfile == "-")
	  par::ped_from_stdin = true;
      }
  
  if (a.find("--map")) par::mapfile = a.value("--map");
    

  if (a.find("--bed"))
    {
      par::read_bitfile = true;
      par::bitfilename = par::pedfile = a.value("--bed");
    }
  
  if (a.find("--fam"))
  {
    par::famfile = a.value("--fam");
  }
  
  if (a.find("--bim"))
  {
    par::read_bitfile = true;
    par::bitfilename_map = par::mapfile = a.value("--bim");
  }
  

  // Single phenotype in file specified
  if (a.find("--pheno")) 
    {
      par::pheno_file = true;
      par::pheno_filename = a.value("--pheno");
    }

  // Multiple phenotypes 
  if (a.find("--mult-pheno"))
    {
      par::multiple_phenotypes = true;
      par::multiple_phenotype_file = a.value("--mult-pheno");      
    }

  if (a.find("--mult-pheno-number"))
    {
      par::plist_selection = par::plist_selection_number = true;
      par::plist_selection_string = a.value("--mult-pheno-number");
    }
  
  if ( a.find("--mult-pheno-name") )
    {
      par::plist_selection = par::plist_selection_name = true;
      par::plist_selection_string = a.value("--mult-pheno-name");
    }

  

  // Get phenotype from a cluster file
  if (a.find("--make-pheno"))
    {
      par::make_pheno = true;
      if ( a.find("--pheno") || a.find("--all-pheno"))
	error("Incompatible phenotype selection commands specified");

      vector<string> s = a.value("--make-pheno",2);      
      par::make_pheno_filename  = s[0];
      par::make_pheno_value = s[1];

      if ( par::make_pheno_value == "*" )
	par::make_pheno_present = true;
      
    }

 
  // Binary 0/1 coding instead of 1/2
  if (a.find("--1"))
    par::coding01 = true;


  // Multiple phenotypes in a file specified
  if (a.find("--mpheno")) 
    {
      if (!a.find("--pheno"))
	error("You need to specify --pheno {file} with --mpheno {N}");
      par::mult_pheno = a.value_int("--mpheno");
    }


  // Select phenotype ny name
  if (a.find("--pheno-name")) 
    {
      if (!a.find("--pheno"))
	error("You need to specify --pheno {file} with --pheno-name {name}");
      if (a.find("--mpheno"))
	error("You cannot specify --mpheno and --pheno-name together");

      par::name_pheno = a.value("--pheno-name");
      
    }

  if (a.find("--values"))
    {
      par::number_list_string = a.value("--values");
    }

  if (a.find("--valueless"))
    {
      par::number_list_string = a.value("--valueless");
      par::number_list_positive = false;
    }

  // Loop over all phenotypes
  if (a.find("--all-pheno")) 
    {
      if (!a.find("--pheno"))
	error("You need to specify --pheno {file} with --all-pheno");
      
      if (a.find("--mpheno"))
	error("You cannot specify --mpheno {N} with --all-pheno");
      par::mult_pheno = 1;
      par::all_pheno = true;

    }

  
  // Single covariate in file specified
  if (a.find("--covar")) 
    { 

      par::covar_file = true;

      if ( a.find("--gxe") )
	{
	  par::covar_filename = a.value("--covar");	  
	  checkFileExists(par::covar_filename);
	  par::clist = false;
	}
      else
	{
	  // Aside from old "--gxe" method, all other
	  // options now use the new clist format

	  // Multiple covariates in file specified, read all of them

	  par::clist = true;
	  
	  par::clist_filename = a.value("--covar");
	  checkFileExists(par::clist_filename);

	  // Request to dump back out all covariates
	  if (makedata)
	    {
	      par::clist = true;
	      if (a.find("--dummy-coding"))
		par::dump_covar_dummy_coding = true;
	    }

	  if (a.find("--with-phenotype"))
	    par::dump_covar_with_phenotype = true;
	}
    }
  

  // Multiple covariates in file specified, select one
  if (a.find("--mcovar")) 
    {
      if (!a.find("--covar"))
	error("You need to specify --covar {file} with --mcovar {N}");
      par::mult_covar = a.value_int("--mcovar");
    }

  // Selection fields for covariates
  if (a.find("--covar-number"))
    {
      par::clist_selection = par::clist_selection_number = true;
      par::clist_selection_string = a.value("--covar-number");
    }
  
  if ( a.find("--covar-name") )
    {
      par::clist_selection = par::clist_selection_name = true;
      par::clist_selection_string = a.value("--covar-name");
    }
    
  // Request to dump back out all covariates
  if (a.find("--write-covar"))
  {
    if (makedata)
      error("No need to specify --write-covar separately");
    if (!a.find("--covar")) 
      error("You must specify a --covar {file} with --write-covar");

    if (a.find("--with-phenotype"))
      par::dump_covar_with_phenotype = true;

    if (a.find("--dummy-coding"))
      par::dump_covar_dummy_coding = true;
    
    par::dump_covar = true;
    par::clist = true;
  }
  
  // Request to dump back out all covariates
   if (a.find("--write-cluster"))
   {
     if (!a.find("--within")) 
       error("You must specify a --within {file} with --write-covar");   
     par::dump_clst = true;
   }

  if (a.find("--write-snplist"))
    par::write_snplist = true;

  if (a.find("--update-map"))
    {
      if (a.find("--update-cm"))
	par::update_cm = true;
      else if (a.find("--update-chr"))
	par::update_chr = true;
      else if ( a.find("--update-name"))
	par::update_name = true;

      par::update_map = true;
      par::update_mapfile = a.value("--update-map");
    }

  if (a.find("--update-ids"))
    {
      par::update_ids = true;
      par::update_ids_file = a.value("--update-ids");
    }
  

  if (a.find("--update-sex"))
    {
      if ( a.find("--update-ids"))
	error("Cannot --update-ids at same time as --update-sex");
      par::update_sex = true;
      par::update_sex_file = a.value("--update-sex");
    }

  if (a.find("--update-parents"))
    {
      if ( a.find("--update-ids"))
	error("Cannot --update-ids at same time as --update-parents");
      par::update_parents = true;
      par::update_parents_file = a.value("--update-parents");
    }

  if (a.find("--update-pheno"))
    {
      if ( a.find("--update-ids"))
	error("Cannot --update-ids at same time as --update-pheno");
      par::update_pheno = true;
      par::update_pheno_file = a.value("--update-pheno");
    }

  if (a.find("--update-alleles"))
    {
      if ( a.find("--update-name"))
	error("Cannot --update-alleles at same time as --update-name");

      par::update_alleles = true;
      par::update_allele_file = a.value("--update-alleles");
    }


  // Examine only a subset of the data?
  if (a.find("--filter"))
    {
      par::filter_on_covar = true;

      vector<string> s = a.value("--filter",2);
      
      par::filter_filename  = s[0];
      par::filter_value = s[1];
      checkFileExists( par::filter_filename );
    }
  
  if (a.find("--mfilter"))
    {
      if (!a.find("--filter"))
	error("You can only specify --mfilter with --filter\n");
      par::mult_filter = a.value_int("--mfilter");
    }
  
  // Different species other than human?
  // i.e. alters chromosome definitions
  
  if (a.find("--dog")) par::species_dog = true;
  if (a.find("--cow")) par::species_cow = true;
  if (a.find("--mouse")) par::species_mouse = true;
  if (a.find("--sheep")) par::species_sheep = true;
  if (a.find("--horse")) par::species_horse = true;
  if (a.find("--rice")) par::species_rice = true;


  //////////////////////////////////
  // Multipoint and singlepoint 

  //  if (a.find("--singlepoint")) par::singlepoint = true; 

  if (a.find("--fringe")) 
    {
      par::singlepoint = false;
      par::fringe = a.value_double("--fringe");
    }  

  if (a.find("--grid")) 
    {
      par::singlepoint = false;
      par::grid = a.value_double("--grid");
      par::inter_grid = 0;
    }

  if (a.find("--step")) 
    {
      par::singlepoint = false;
      par::inter_grid = a.value_int("--step");
    }

  if (a.find("--cm")) par::cm_map = true;
    
  if (a.find("--ci")) 
    {
      par::display_ci = true;
      par::ci_level = a.value_double("--ci");
      if ( par::ci_level < 0.01 || par::ci_level >= 1 )
	error("CI level (--ci) must be between 0 and 1\n");
      par::ci_zt = ltqnorm( 1 - (1 - par::ci_level) / 2  );
    }

  if (a.find("--pfilter")) 
    {
      par::pfilter = true;
      par::pfvalue = a.value_double("--pfilter");
    }

  if (a.find("--hide-covar"))
    par::no_show_covar = true;


  if (a.find("--meta-analysis"))
    {
      par::meta_analysis = true;
      par::meta_files = a.varValue("--meta-analysis");      
    }

  if (a.find("--annotate"))
    {
      par::annot_file = true;
      par::annot_filename = a.value("--annotate");
    }
    
  if (a.find("--gene-report"))
    {
      par::greport = true;

      par::greport_results = a.value("--gene-report");
      
      if (!a.find("--gene-list") )
	error("You must specify a --gene-list");
      else
	par::greport_gene_list = a.value("--gene-list");
      
      if ( a.find("--gene-list-border") )
	par::make_set_border = a.value_int("--gene-list-border") * 1000;      

      if ( a.find("--gene-subset") )
	{
	  par::greport_subset = true;
	  par::greport_subset_file = a.value("--gene-subset");
	}	

      if ( a.find("--gene-report-empty") )
	par::greport_display_empty = true;

    }
  

  if (a.find("--show-tags"))
    {
      par::gettag_mode = true;
      par::gettag_file = a.value("--show-tags");

      if (a.find("--list-all") || par::gettag_file == "all" )
	par::gettag_listall = true;

      if ( a.find("--tag-mode2") )
	{
	  par::gettag_mode1 = false;
	  par::gettag_mode2 = true;
	}

      if ( a.find("--tag-r2"))
	par::gettag_r2 = a.value_double("--tag-r2");
      if ( a.find("--tag-kb"))
	par::gettag_kb = a.value_int("--tag-kb") * 1000;
    }


  if (a.find("--clump"))
    {
      par::clumpld = true;
      par::clumpld_results = a.value("--clump");
      
      if (a.find("--clump-best"))
	par::clumpld_best = true;

      if (a.find("--clump-field"))
	par::clumpld_column = a.value("--clump-field");      

      if (a.find("--clump-verbose"))
	par::clumpld_verbose = true;
      
      if (a.find("--clump-p1"))
	par::clumpld_p1 = a.value_double("--clump-p1");

      if (a.find("--clump-p2"))
	par::clumpld_p2 = a.value_double("--clump-p2");

      if (a.find("--clump-r2"))
	par::clumpld_r2 = a.value_double("--clump-r2");

      if (a.find("--clump-kb"))
	par::clumpld_kb = a.value_int("--clump-kb") * 1000;

      if (a.find("--clump-index-first"))
	  par::clumpld_index1 = true;
      
      if (a.find("--clump-replicate"))
	  par::clumpld_only_show_replications = true;
      
      if (a.find("--clump-range"))
      {
	  
	  if ( a.find("--make-set"))
	      error("Cannot specify --make-set and --clump-range together");
	  
	  par::clumpld_range_annotate = true;
	  par::clumpld_range_file = a.value("--clump-range");
	  if ( a.find("--clump-range-border"))
	      par::make_set_border = a.value_int("--clump-range-border") * 1000;      
	  
      }
      
      if (a.find("--clump-only-non-index"))
	{
	  par::clumpld_only_show_replications = true;
	  par::clumpld_only_show_replications_list = true;
	}

      if (a.find("--clump-annotate"))
	{
	  par::clumpld_annot = true;
	  par::clumpld_annot_fields = a.value("--clump-annotate");
	}
      
      if (a.find("--clump-allow-overlap"))
        par::clumpld_indep = false;

    }


  if (a.find("--adjust")) 
    { par::multtest = true; }

  if (a.find("--log10")) 
    { par::logscale = true; }

  if (a.find("--qq-plot")) 
    { par::qq_plot = true; }

  if (a.find("--lambda"))
    {
      par::fix_lambda = true;
      par::lambda = a.value_double("--lambda");
      if ( par::lambda < 1 ) 
	par::lambda = 1;
    }

  if (a.find("--gc"))
    {
      if (!a.find("--adjust"))
	error("Must specify --adjust to use --gc");
      par::use_GC = true;
    }

  ///////////////////////
  // Permutation options


  // Use permutations (default is adaptive)
  if (a.find("--perm"))
    {
      par::permute = true;
    }

  // Return counts not p-values (i.e. number of times exceeded)
  if ( a.find("--perm-count") )
    {
      par::perm_count = true;
    }

  // Specify parameters for adaptive permutation
  if (a.find("--aperm"))
    {
      if (a.find("--segment"))
	error("--segment options requires --pperm option");
      
      if (a.find("--set")||a.find("--make-set"))
	error("Cannot use --aperm with SET options (use --mperm N instead)");
      
      par::permute = true;
      par::adaptive_perm = true;
      vector<string> s = a.value("--aperm",6);

      par::adaptive_min       = getInt(s[0].c_str(),"--aperm");
      par::adaptive_max       = getInt(s[1].c_str(),"--aperm");
      par::adaptive_alpha     = getDouble(s[2].c_str(),"--aperm");
      par::adaptive_ci        = getDouble(s[3].c_str(),"--aperm");
      par::adaptive_interval  = getInt(s[4].c_str(),"--aperm");
      par::adaptive_interval2 = getDouble(s[5].c_str(),"--aperm");
    }


  // Non-adaptive (maxT) permutations
  if (a.find("--mperm")) 
    {
      par::permute = true;
      par::adaptive_perm = false;
      par::replicates = a.value_int("--mperm");

      if (a.find("--mperm-save"))
	par::mperm_save_best = true;
      else if (a.find("--mperm-save-all"))
	par::mperm_save_all = true;


      // But make special fix for QFAM tests
      if ( par::QTDT_test )
	{
	  par::adaptive_perm = true;
	  par::adaptive_min       = par::replicates;
	  par::adaptive_max       = par::replicates;
	  par::adaptive_alpha     = 0;
	  par::adaptive_ci        = 0;
	  par::adaptive_interval  = par::replicates+1;
	  par::adaptive_interval2 = 0;
	  par::QFAM_adaptive = true;
	}
    }

  if (a.find("--make-perm-pheno"))
    {
      if ( a.find("--mperm") || a.find("--perm") )
	error("Cannot specify --make-perm-pheno with other permutation options");
      par::output_pheno_perm = true;
      par::permute = true;
      par::adaptive_perm = false;
      par::replicates = a.value_int("--make-perm-pheno");
    }


  if (a.find("--rank"))
    {
      if (! a.find("--mperm") )
	error("--rank requires --mperm to be specified");
      par::mperm_rank = true;
    }

  // PLINK permutations

  if (a.find("--pperm")) 
    {
      if (! ( a.find("--segment") || a.find("--read-segment") ) )
	error("--pperm options requires --segment or --read-segment option");
      par::permute = true;
      par::adaptive_perm = false;
      par::replicates = a.value_int("--pperm");
    }


  if (a.find("--p2"))
    {
      if ((!a.find("--perm")) 
	  && (!a.find("--mperm"))
	  && (!a.find("--aperm")))
	error("--p2 option also requires--perm, --aperm or --mperm"); 
      else if (!a.find("--assoc"))
	error("The --p2 option can only be specified with --assoc");
      else par::assoc_test_alt_perm = true;
    }


  /////////////////////
  // Gene-dropping

  if (a.find("--genedrop")) 
    par::perm_genedrop = par::permute = true;

  if (a.find("--swap-sibs")) 
    {
      if (!a.find("--genedrop"))
	error("--swap-sibs only makes sense when --genedrop specified");
      par::perm_genedrop_sibships = true;      
      par::perm_genedrop_and_swap = true;
    }

  if (a.find("--swap-parents")) 
    {
      if (!a.find("--genedrop"))
	error("--swap-parents only makes sense when --genedrop specified");
      par::perm_genedrop_parents = true;      
      par::perm_genedrop_and_swap = true;
    }

  if (a.find("--swap-unrel")) 
    {
      if (!a.find("--genedrop"))
	error("--swap-unrel only makes sense when --genedrop specified");

      par::perm_genedrop_unrel = true;
      par::perm_genedrop_and_swap = true;
    }


  ///////////////////////
  // Misc. options

  if (a.find("--compound-genotypes"))
    {
      if ( ! ( a.find("--file") || a.find("--ped") || a.find("--lfile") || a.find("--lgen") ) )
	error("--compound-genotype only works with PED/MAP or LGEN filesets currently");
      par::compound_genotype_code = true;
    }

  if (a.find("--allele-count"))
    {
      if ( ! ( ( a.find("--lfile") || a.find("--lgen")  ) && a.find("--reference") ) )
	error("Can only use --allele-count with --lgen and --reference\n");
      par::lfile_allele_count = true;
      // expect either 1 or 2, to indicate # of non-reference alleles (i.e. # of mutations)	   
    }

  if (a.find("--missing-genotype")) 
    {
      par::missing_genotype = a.value("--missing-genotype");
      par::out_missing_genotype = par::missing_genotype;
      par::missing_genotype_explicit = true;
    }

  
  if (a.find("--missing-phenotype")) 
    {
      par::missing_phenotype = a.value("--missing-phenotype");
      par::out_missing_phenotype = par::missing_phenotype;
      par::missing_phenotype_explicit = true;
    }
  
  if (a.find("--output-missing-genotype")) 
    {
      par::out_missing_genotype = a.value("--output-missing-genotype");
    }

  if (a.find("--output-missing-phenotype")) 
    {
      par::out_missing_phenotype = a.value("--output-missing-phenotype");
    }


  if (a.find("--FIX")) {
    par::FIXED = par::FIXED_p = true;
    vector<string> p = a.value("--FIX",4);
    
    par::FIX_IBD.z0 = getDouble(p[0].c_str(),"--FIX");
    par::FIX_IBD.z1 = getDouble(p[1].c_str(),"--FIX");
    par::FIX_IBD.z2 = getDouble(p[2].c_str(),"--FIX");
    par::FIX_p = getDouble(p[3].c_str(),"--FIX");
    cout << "Fixing Z0, Z1, Z2 and p to " 
	 << par::FIX_IBD.z0 << " "
	 << par::FIX_IBD.z1 << " "
	 << par::FIX_IBD.z2 << " "
	 << par::FIX_p << "\n(p must refer to '1' allele in '1/2' genotype)\n";
  }


  if (a.find("--fix-ibd")) {
    par::FIXED = true;
    vector<string> p = a.value("--fix-ibd",3);
    
    par::FIX_IBD.z0 = getDouble(p[0].c_str(),"--fix-ibd");
    par::FIX_IBD.z1 = getDouble(p[1].c_str(),"--fix-ibd");
    par::FIX_IBD.z2 = getDouble(p[2].c_str(),"--fib-ibd");

  }


  if (a.find("--batch")) 
    par::BATCH_SIZE = a.value_int("--batch");

  if (a.find("--min")) 
    par::MIN_PIHAT = a.value_double("--min"); 
  
  if (a.find("--max")) 
    par::MAX_PIHAT = a.value_double("--max"); 
  
  if (a.find("--all-pairs"))
    {
      if (a.find("--min"))
	error("Cannot specify --min and --all-pairs\n");      
      par::include_all_pairs = true;      
    }


//   if (a.find("--lock"))
//     { par::locked = true; }
  
//   if (a.find("--unlock"))
//     { par::locked = false; }


  ///////////////////////////////
  // Basic filters: make this the 
  // default now...

  if ( a.find("--all") || true )
    {
      par::min_af = 0.0;
      par::MAX_GENO_MISSING = 1;
      par::MAX_IND_MISSING = 1;
    }
 
  if (a.find("--geno")) 
    par::MAX_GENO_MISSING = a.value_double("--geno"); 

  if (a.find("--mind")) 
    par::MAX_IND_MISSING = a.value_double("--mind"); 


  if (a.find("--maf")) 
    par::min_af = a.value_double("--maf");

  if (a.find("--max-maf")) 
    {
      par::max_af = a.value_double("--max-maf");
      if (par::max_af < par::min_af)
	error("Cannot set --max-maf less than --maf\n");
    }

  if (a.find("--keep-allele-order"))
    {
      par::make_minor_allele = false;
    }
  
  if (a.find("--mhf"))
    par::min_hf = a.value_double("--mhf");
  
  if (a.find("--max-mhf"))
    {
      par::max_hf = a.value_double("--max-mhf");
      if (par::max_hf < par::min_hf)
	error("Cannot set --max-mhf less than --mhf\n");
    }

  if (a.find("--hwe")) 
    {
      par::HWD_test = true;
      par::HWD_limit = a.value_double("--hwe");
    }

  if (a.find("--hwe2")) 
    {
      par::HWD_test = true;
      par::HWD_standard = true;
      par::HWD_limit = a.value_double("--hwe2");
    }

  if (a.find("--hwe-all")) 
    {
      par::HWD_filter_on_all = true;
    }


  if (a.find("--me")) 
    {
      par::MENDEL_test = true;
      vector<string> s = a.value("--me",2);
      par::MENDEL_ind = getDouble(s[0],"--me");
      par::MENDEL_snp = getDouble(s[1],"--me");
    }


  //////////////////////////////////
  // Reading a dosage file

  if (a.find("--dosage"))
    {
      par::dosage_assoc = true;      
      par::dosage_file = a.value("--dosage");
      if ( ! a.find("--fam") )
	error("You need to also specify a FAM (--fam) file");
      if ( a.find("--map"))
	par::dosage_hasMap = true;
      if ( a.find("--hard-call"))
	{
	  if ( ! a.find("--map") )
	    error("Need to specify --map with --hard-call");
	  
	  par::dosage_hard_call = true;

	  vector<string> s = a.value("--hard-call",2);
	  par::dosage_hard_call_thresh = getDouble(s[0].c_str(),"--hard-call");
	  par::dosage_hard_call_thresh2 = getInt(s[1].c_str(),"--hard-call");
	}
      else if ( a.find("--write-dosage") )
	par::write_dosage = true;
    }


  //////////////////////////////////
  // IBS clustering

  if (a.find("--cluster")) 
    {
      par::cluster = true;
      if (a.find("--within"))
	par::force_initial_cluster = true;

      if (a.find("--group-avg") || a.find("--group-average"))
	{
	  par::cluster_group_avg = true;
	}      
    }



  if (a.find("--euclidean")) 
    {
      if (!a.find("--cluster"))
	error("Cannot specify --euclidean without --cluster");
      par::cluster_euclidean = true;
    }



  if (a.find("--pick1")) 
    {
      par::cluster_selcon = true;
      par::cluster_selcon_file = a.value("--pick1");
    }

  if (a.find("--cluster-missing")) 
    {
      par::cluster = true;
      par::cluster_missing = true;
      par::matrix = true;
      
      if (!a.find("--maf")) par::min_af = 0.0;
      if (!a.find("--geno")) par::MAX_GENO_MISSING = 1;
      if (!a.find("--mind")) par::MAX_IND_MISSING = 1;
 
      par::merge_p = 0;
      if (a.find("--ppc"))
	error("Cannot specify --ppc with --cluster-missing");
    }

  if (a.find("--K")) 
    {
      if (!a.find("--cluster"))
	error("Must specify --cluster also if --K used");
      par::max_cluster_N = a.value_int("--K");     
    }
  
  if (a.find("--neighbour")) 
    {
      vector<string> s = a.value("--neighbour",2);
      par::min_neighbour = getInt(s[0].c_str(),"--neighbour");
      par::max_neighbour = getInt(s[1].c_str(),"--neighbour");
      par::outlier_detection = true;
    }

  if (a.find("--matrix")) par::matrix = true;
  if (a.find("--distance-matrix")) { par::distance_matrix = par::matrix = true; }

  if (a.find("--mds-plot")) 
    {
      par::cluster_plot = true;
      par::cluster_mds_dim = a.value_int("--mds-plot");
    }
  
  if (a.find("--mds-cluster"))
    par::mds_by_individual = false;  
      
  if (a.find("--pmerge"))
    error("--pmerge is depreciated: use --ppc instead\n");
  
  if (a.find("--ppc")) 
    {
      if (!par::cluster) error("--ppc options requires --cluster");
      par::merge_p = a.value_double("--ppc");
    }

  if (a.find("--pibs-gap"))
    error("--pibs-gap is depreciated: please use --ppc-gap\n");
  
  if (a.find("--ppc-gap"))
    {
      par::ibstest_gap = 1000 * a.value_int("--ppc-gap");
    }
  
  if (a.find("--ibm"))
    {
      if (! a.find("--cluster"))
	error("Can only use --ibm with --cluster\n");
      par::cluster_ibm_constraint = true;
      par::cluster_ibm_constraint_value = a.value_double("--ibm");      
    }

  if (a.find("--mc")) 
    {
      if (!par::cluster) error("--mc options requires --cluster");
      par::max_cluster_size = a.value_int("--mc");
    }

  if (a.find("--cc")) 
    {
      if (!par::cluster) error("--cc options requires --cluster");
      par::cluster_on_phenotype = true;
    }
  
  if (a.find("--mcc")) 
    {
      if (!par::cluster) error("--mcc options requires --cluster");
      par::cluster_on_mcc = true;
      vector<string> s = a.value("--mcc",2);
      par::max_cluster_case = getInt(s[0].c_str(),"--mcc");
      par::max_cluster_control = getInt(s[1].c_str(),"--mcc");
      if (a.find("--mc") || a.find("-cc"))
	error("Cannot specify --mc N and/or --cc as well as --mcc N1 N2\n");
    }



  /////////////////////////////////
  // External criteria to match on

  // Categorical binary traits,
  // by default
  //  e.g. { A, A } is a match and so are pairable
  //       { A, B } is not 
  // 
  //  if match-type file is also specifed, then matches 
  //  can potentially be otherwise, e.g. 
  //       { A, B } are pairable
  //       { A, A } are not

  if (a.find("--match")) 
  {
     par::bmatch = true;     
     par::bmatch_filename = a.value("--match");
  }

  if (a.find("--match-type"))
    {
      if (!a.find("--match"))
	error("Must specify a --match {file} with the --match-type {file} option");
      par::bmatch_usertype = true;
      par::bmatch_direction_filename = a.value("--match-type");
    }

  // Quantitative trait match
  // Based on difference exceeding a certain threshold
  // e.g. (X-Y)>T    => no match
  //      (X-Y)<=T   => match

  // T is specified by including an extra individual in the qmatch file
  // with the Family ID and Individual ID "_T_" 

  if (a.find("--qmatch"))
  {
     if (!a.find("--qt"))
      error("You need to specify a --qt file when using --qmatch");
     par::qmatch_threshold_filename = a.value("--qt");    
     par::qmatch = true;     
     par::qmatch_filename = a.value("--qmatch");
  }



  //////////////////////////
  // Permutation clustering

  if (a.find("--family")) 
    {
      par::sol_family = true;
      par::permute_within_sol = true;
      par::include_cluster = true;
    }
  
  if (a.find("--within"))
    {
      par::permute_within_sol = par::include_cluster = true;            	
      par::include_cluster_from_file = true;
      par::include_cluster_filename = a.value("--within");
      checkFileExists(par::include_cluster_filename);
    }

  if (a.find("--mwithin"))
    {
      if (!a.find("--within"))
	error("You can only specify --mwithin with --within");
      par::mult_clst = a.value_int("--mwithin");
    }

  
  //////////////////////////////////
  // Specific scan region selected
  
  if (!a.find("--from"))
    {
      
      // Specify a specific chromosome
      string c="0";  // Default all chromosomes
      
      if (a.find("--chr")) 
	c = a.value("--chr");
      
      if (c=="X" || c=="x") par::run_chr = 23;
      else if (c=="Y" || c=="y") par::run_chr = 24;
      else par::run_chr = getInt(c,"--chr");
    }
      
  if (a.find("--from")) 
    {
      
      if (!a.find("--to"))
	error("Must also specify --to {marker} when using --from {marker}");
	
      par::m1 = a.value("--from");
      par::m2 = a.value("--to");
      par::run_chr = -1;
    }
  
  if (a.find("--snp"))
    {
      par::m1 = a.value("--snp");
      par::m2 = a.value("--snp");
      par::run_chr = -1;      
    }

  if (a.find("--snps"))
    {
      par::extract_set = true; 
      par::snp_include_from_cl = true;
      par::snp_include_range = a.value("--snps");

      if ( a.find("--snp") || a.find("--window") || a.find("--extract") || a.find("--exclude") )
	error("Cannot specify multiple SNP-selection options with --snps");
      
    }

  if ( a.find("--d") )
    {
      par::range_delimiter = a.value("--d");
      if ( par::range_delimiter.length() > 1 )
	error("Range delimiter can only be 1 character");
      if ( par::range_delimiter == "," )
	error("Cannot set range delimiter to comma");
    }

  if (a.find("--window"))
    {
      if (!a.find("--snp")) error("Must specify --snp with --window");
      par::window = a.value_double("--window");
    }

  if (a.find("--from-bp"))
    {
      if (!a.find("--to-bp")) error("Must specify --to-bp with --from-bp");
      par::from_window = a.value_int("--from-bp");
      par::position_window = true;
    }

  if (a.find("--from-kb"))
    {
      if (!a.find("--to-kb")) error("Must specify --to-kb with --from-kb");
      par::from_window = int(a.value_double("--from-kb") * 1000);
      par::position_window = true;
    }
  
  if (a.find("--from-mb"))
    {
      if (!a.find("--to-mb")) error("Must specify --to-mb with --from-mb");
      double v = a.value_double("--from-mb");
      if (v>1000) error("Too large a value for --from-mb");
      par::from_window = int(v * 1000 * 1000);
      par::position_window = true;
    }
  
  
  if (a.find("--to-bp"))
    {
      if (!a.find("--from-bp")) error("Must specify --from-bp with --to-bp");
      par::to_window = a.value_int("--to-bp");
      par::position_window = true;
    }
  
  if (a.find("--to-kb"))
    {
      if (!a.find("--from-kb")) error("Must specify --from-kb with --to-kb");
      par::to_window = int(a.value_double("--to-kb") * 1000);
      par::position_window = true;
    }

  if (a.find("--to-mb"))
    {
      if (!a.find("--from-mb")) error("Must specify --from-mb with --to-mb");
      double v = a.value_double("--to-mb");
      if (v>1000) error("Too large a value for --to-mb");
      par::to_window = int(v * 1000 * 1000);
      par::position_window = true;
    }

  if (par::position_window)
    if (!a.find("--chr"))
      error("You must specify which chromosome (--chr N) also");
  


  ////////////////////////////////////////////
  // General warnings

  if ( a.find("--assoc") && a.find("--covar") )
    error("Cannot specify --covar with --assoc");

  if ( a.find("--model") && a.find("--covar") )
    error("Cannot specify --covar with --model");

  if ( a.find("--assoc") && a.find("--linear") )
    error("Cannot specify --assoc with --linear");
  
  if ( a.find("--assoc") && a.find("--logistic") )
    error("Cannot specify --assoc with --logistic");

  if ( a.find("--model") && a.find("--linear") )
    error("Cannot specify --model with --linear");
  
  if ( a.find("--model") && a.find("--logistic") )
    error("Cannot specify --model with --logistic");

  if ( a.find("--model") && a.find("--assoc") )
    error("Cannot specify --model with --assoc");


  /////////////////////////////////////////////
  //  Help -- display all options

  
  if (a.find("--help") || a.find("-h"))
    {

      cout << "\n"
	   << "Please visit the PLINK website for a complete list of options\n"
	   << "\n"
	   << "A few common options are listed here:\n"
	   << "\n";

      cout << "plink --file {fileroot}         Specify .ped and .map files \n"
	   << "      --bfile {fileroot}        Specify .bed, .fam and .map \n"
	   << "\n"
	   << "      --out {fileroot}          Specify output root filename\n"
	   << "\n"
	   << "      --missing-genotype {0}    Missing genotype code       \n"
	   << "      --missing-phenotype {-9}  Missing phenotype code      \n"
	   << "\n"
	   << "      --pheno {phenofile}       Specify .phe file           \n"
	   << "      --within {file}           Specify cluster file        \n"
	   << "      --cov {covarfile}         Specify .cov file           \n"
	   << "\n"
	   << "      --extract {snplist}       Extract list of SNPs        \n"
	   << "      --exclude {snplist}       Exclude list of SNPs        \n"
	   << "      --remove {indlist}        Remove these individuals    \n"
	   << "      --keep {indlist}          Keep these individuals      \n"
	   << "\n"
	   << "      --make-bed                Make .bed, .fam and .bim    \n"
	   << "      --recode                  Output new PED and MAP files\n"
	   << "      --recode12                As above, with 1/2 alleles  \n"
	   << "      --recodeAD                As above, but: 1/0/-1, 0/1/0\n"
	   << "      --recodeA                 As above, but: 1/0/-1 only  \n"
	   << "\n"
	   << "      --snp {marker}            Specify this single SNP     \n"
	   << "      --snps {marker list}      Specify list,range of SNPs  \n"
	   << "      --window {kb}             Select +/- kb around --snp  \n"
	   << "      --chr {N}                 Analyse chromosome          \n"
	   << "      --from-kb {KB}            Start scan here (kilobase)  \n"
	   << "      --to-kb {KB}              End scan here               \n"
	   << "\n"
           << "      --all                     Set filters to include all  \n"
	   << "      --maf {0.01}              Minor allele frequency      \n"
	   << "      --geno {0.1}              Maximum per-SNP missing     \n"
	   << "      --mind {0.1}              Maximum per-person missing  \n"
	   << "\n"
	   << "      --freq                    Output allele frequencies   \n"
	   << "      --hardy                   Hardy-Weinberg tests        \n"
	   << "      --missing                 Genotyping rate information \n"
	   << "      --het                     Individual inbreeding       \n"
	   << "      --genome                  Genome-wide IBS/IBD         \n"
	   << "      --cluster                 Perform IBS clustering      \n"
	   << "\n"
	   << "      --assoc                   Case/control, QT association\n"
	   << "      --model                   Full-model C/C association  \n"
	   << "      --tdt                     Family-based TDT association\n"
	   << "      --linear                  Linear regression model     \n"
	   << "      --logistic                Logistic regression model   \n"
	   << "\n"
 	   << "      --perm                    Apaptive permutations       \n"
 	   << "      --mperm {1000}            max(T) permutations         \n"
	   << "\n"
	   << "      --hap {tagfilename}       Multimarker predictor list  \n"
	   << "      --hap-window {N}          Phase sliding window        \n"
	   << "      --hap-snps {snp list}     Phase this set of SNPs      \n" 
	   << "      --hap-assoc               Haplotype-based association \n"
	   << "      --hap-tdt                 Haplotype-based TDT         \n"
	   << "      --chap                    Conditional haplotype tests \n"
	   << "      --hap-phase               Report haplotype phases     \n"
	   << "      --hap-freq                Report haplotype frequencies\n"
	   << "\n";      	
      
      cout << "\nPlease visit the PLINK website for a complete list of options\n\n";

      shutdown();
    }
  
  
  // By default, most tests are SNP major
  
  par::SNP_major = true;
  
  // Exceptions are:
  //  TDT  ( family structure confuses things)
  //  Whole genome / IBS clustering 
  //  PLINK
  
  if (par::TDT_test || 
      par::MENDEL_test || 
      par::MENDEL_report || 
      par::genome_output || 
      par::cluster || 
      par::plink )
    par::SNP_major = false;
  
  if (a.find("--ind"))
    par::SNP_major = false;
  

  // If recoding data, the default will be not to set heterozygous
  // haploid genotypes to missing. Likewise for Mendel errors. Merge
  // operations will also specify a recode/make-bed, so they are also
  // captured here. The one special case where we want to allow to
  // preserve males hets on the X is the --check-sex

  
  if ( par::check_sex )
    par::preserve_all_genotypes = true;
	  
  if ( par::write_bitfile || 
       par::recode || 
       par::recode_HV ||
       par::recode_whap ||
       par::recode_12 ||
       par::recode_AD )
    {
      
      // Unless flag given, these options will not replace haploid
      // heterozygotes with a missing genotype
      
      if ( a.find("--set-hh-missing") )
	par::preserve_all_genotypes = false;
      else
	par::preserve_all_genotypes = true;

      if ( a.find("--set-me-missing") )
	{
	  par::preserve_mendel_errors = false;
	}
      else
	par::preserve_mendel_errors = true;
      
    }


    }

