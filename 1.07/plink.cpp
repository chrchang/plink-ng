

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
#include <vector>
#include <map>
#include <list>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>

#include "plink.h"
#include "options.h"
#include "helper.h"
#include "crandom.h"
#include "perm.h"
#include "sets.h"
#include "linear.h"
#include "logistic.h"
#include "phase.h"
#include "clumpld.h"
#include "nlist.h"
#include "sets.h"
#include "stats.h"
#include "idhelp.h"
#include "zed.h"

using namespace std;

ofstream LOG;
string PVERSION;
string PDATE;
string PREL;
Plink * PP;
map<string,int> Range::groupNames;

int main(int argc, char* argv[]) 
{


  
  /////////////////////////
  // Setup, display title

  cout.setf(ios::fixed);
  cout.precision(8);
  
  set_new_handler(NoMem);

  PVERSION = "1.07";        // 4 chars
  PREL = " ";               // space or p (full, or prelease) 
  PDATE    = "10/Aug/2009"; // 11 chars


  //////////////////
  // The major class

  Plink P;
  PP = &P;



  /////////////////////////////////////////////////////
  // General class for all haplotype-related functions

  P.haplo = new HaploPhase(P);


  //////////////////////////
  // Command line arguments

  CArgs a(argc,argv);
  getOutputFilename(a);


  //////////////////////////
  // Start logging, title

  LOG.open(string(par::output_file_name + ".log").c_str());
  P.printLOG("\n"
	     "@----------------------------------------------------------@\n"
	     "|        PLINK!       |     v"+PVERSION+PREL+"     |   "+PDATE+"     |\n"
	     "|----------------------------------------------------------|\n"
	     "|  (C) 2009 Shaun Purcell, GNU General Public License, v2  |\n"
	     "|----------------------------------------------------------|\n"
	     "|  For documentation, citation & bug-report instructions:  |\n"
             "|        http://pngu.mgh.harvard.edu/purcell/plink/        |\n"
             "@----------------------------------------------------------@\n"
	     "\n");



  //////////////////////////
  // Fully parse command line

  setOptions(a);


  /////////////////////
  // Permutation class


  if ( par::random_seed == 0 )
    CRandom::srand(time(0));
  else
    CRandom::srand( par::random_seed );

  Perm perm(P);
  P.pperm = & perm;


  ////////////////
  // Check version 
  
  if (par::web_check)
    P.webcheck(a);
  else
    P.printLOG("Skipping web check... [ --noweb ] \n");
  

  /////////////
  // Time stamp

  P.printLOG("Writing this text to log file [ "+
	     par::output_file_name + ".log ]\n");
    
  time_t curr=time(0);
  string tdstamp = (string)ctime(&curr);
  P.printLOG("Analysis started: " + tdstamp +"\n");



  /////////////////////////////////////
  // Validate and record all arguments

  a.check_unused_options(P);

  if ( par::output_file_name.find(".",0) != string::npos )
    P.printLOG("** For gPLINK compatibility, do not use '.' in --out **\n");

  

  //////////////////////////
  // Some basic definitions

  if (par::species_dog) defineDogChromosomes();
  else if (par::species_sheep) defineSheepChromosomes();
  else if (par::species_cow) defineCowChromosomes();
  else if (par::species_horse) defineHorseChromosomes();
  else if (par::species_rice) defineRiceChromosomes();
  else if (par::species_mouse) defineMouseChromosomes();
  else defineHumanChromosomes();


  ///////////////////////////////
  // Web-based SNPServer lookup?
  
  if ( par::lookup )
    {      
      P.lookup();
      shutdown();
    }

  if ( par::lookup2 )
    {      
      P.lookup2();
      shutdown();
    }


  /////////////////////////
  // ID helper?

  if ( par::idhelp )
    {
      IDHelper ID;
      ID.idHelp();
      shutdown();
    }


  /////////////////////////
  // File compression utility

  if ( par::compress_file )
    {
      fileCompress();
      shutdown();
    }

  if ( par::uncompress_file )
    {
      fileUncompress();
      shutdown();
    }




  //////////////////////////////////////////////////
  // Main Input files

  // Simulate or read in data:
  // Binary or ASCII format; transposed/long/generic

  if (par::dummy) P.dummyLoader();
  else if (par::greport) P.displayGeneReport();
  else if (par::annot_file) P.annotateFile();
  else if (par::meta_analysis) P.metaAnalysis();
  else if (par::rare_test_score_range) P.displayRareRange();
  else if (par::simul) 
    {
      if ( par::simul_qt )
	P.simulateSNPs_QT();
      else
	P.simulateSNPs();      
    }
  else if (par::cnv_list) P.setUpForCNVList();
  else if (par::read_bitfile) P.readBinData();
  else if (par::lfile_input) P.readDataLongFormat();
  else if (par::tfile_input) P.readTransposedData();
  else if (par::read_ped) P.readData();  
  else if (par::gvar) {  
    par::load_gvar=true; 
    P.readGenericVariantData();
  }
  else if ( par::dosage_assoc ) 
    {
      P.readFamFile(par::famfile);
      if ( par::dosage_hasMap )
	{
	  checkFileExists( par::mapfile );
	  vector<bool> include;
	  vector<int> include_pos(0);
	  int nvar = 0;
	  P.readMapFile(par::mapfile,
			include,
			include_pos,
			nvar);
	}
    }

  // Set number of individuals
  P.n = P.sample.size();
  
  // Set number of pairs
  P.np = (int)((double)(P.n*(P.n-1))/(double)2);
  
  // Total number of all (test+background) loci
  P.nl_all = P.locus.size();

  // Number of permutations to store
  P.maxr2.resize(par::replicates);


  // Check for duplicate individual or SNP names
  checkDupes(P);


  /////////////////////////////////////
  // Merge with a secondary data file 
  // Standard (non-list) mode

  if (par::merge_data && !par::merge_list)
    {
      
      if (par::merge_binary)
	P.mergeBinaryData();
      else
	P.mergeData();
      
      // Reset number of individuals
      P.n = P.sample.size();
      
      // Set number of pairs
      P.np = (int)((double)(P.n*(P.n-1))/(double)2);
      
      // Total number of all (test+background) loci
      P.nl_all = P.locus.size();
    }
  

  
  /////////////////////////////////////
  // Merge with a secondary data file 
  // List mode
  
  if (par::merge_list)
    P.mergeList();
      
  
  //////////////////////////////////////////  
  // A different phenotype file specified?
  
  if (par::pheno_file) P.readPhenoFile();
  else if (par::make_pheno) P.makePhenotype();
  else if (par::multiple_phenotypes) P.readMultiplePhenoFile();

  
  ////////////////////////////////
  // Remove any individuals with 
  // missing phenotypes
  
  if (!par::ignore_phenotypes)
    removeMissingPhenotypes(P);
  

  //////////////////////////////////
  // Binary affection status coding

  if (par::bt)
    affCoding(P);
        


  /////////////////////////////////
  // Update MAP file information?

  if (par::update_map)
    P.updateMapFile();

  /////////////////////////////////
  // Update FAM information?

  if (par::update_ids || par::update_parents || par::update_sex || par::update_pheno )
    P.updateFamFile();

  /////////////////////////////////
  // Update allele file information?

  if (par::update_alleles)
    P.updateAlleles();


  /////////////////////////////////
  // Flip DNA strand for any SNPs? 

  if (par::flip_strand)
    P.flipStrand();

      
  /////////////////////////////////
  // Recode any alleles? 

  if (par::recode_ACGT || par::recode_1234)
    P.alleleRecoding(); 
 



  //////////////////////////////////////////////////////////
  // Output a specific set of SNPs (--extract or --exclude)
  
  if ( par::extract_before_exclude )
    {
      if (par::extract_set)
	P.extractExcludeSet(false);      
      if (par::exclude_set)
	P.extractExcludeSet(true);
    }
  else
    {
      if (par::exclude_set)
	P.extractExcludeSet(true);
      if (par::extract_set)
	P.extractExcludeSet(false);      
    }
  

  /////////////////////////////////////////////////////////////
  // Output a specific set of individuals --remove or --keep
    
  if ( par::remove_before_keep )
    {
      if (par::remove_indiv)
	P.removeIndividuals(false);
      if (par::keep_indiv)
	P.removeIndividuals(true);
    }
  else
    {
      if (par::keep_indiv)
	P.removeIndividuals(true);
      if (par::remove_indiv)
	P.removeIndividuals(false);
    }

  ///////////////////////////////////////////////
  // Filter based on attribute files

  if ( par::snp_attrib_filter )
    P.attribFilterSNP();
  
  if ( par::ind_attrib_filter )
    P.attribFilterInd();


  ///////////////////////////////////////////////
  // Filter based on qualiy scores

  if ( par::read_snp_qual )
    P.filterQualSNPs();
  
  if ( par::read_geno_qual )
    P.filterQualGenotypes();
  



  //////////////////////////////////////////////////
  // Pull a random subset of SNPs?
  
  if ( par::thin_snps )
    P.thinSNPs();
    


  /////////////////////////////////////////////////////////////
  // If in --genome list mode, keep the two lists of individuals
  
  if (par::genome_2sets)
    {
      P.keep2SetsForGenome();
    }
  

  ///////////////////////////////////////////////
  // Read a list of obligatory missing genotypes?

  if (par::oblig_missing) 
    P.setObligMissing();
    

  //////////////////////////////////////////////////
  // Filter individuals based on external covariate? 
  
  if (par::filter_on_covar)
    {
      P.filterOnCovariate();
      // Reset number of individuals
      P.n = P.sample.size();
      P.np = (int)((double)(P.n*(P.n-1))/(double)2);
    }


  ////////////////////////////
  // Any simple preset filters
  
  if (par::filter_males)
    P.filterOnMale();
  else if (par::filter_females)
    P.filterOnFemale();

  if (par::filter_cases)
    P.filterOnCase();
  else if (par::filter_controls)
    P.filterOnControl();
  
  if (par::filter_founders)
    P.filterOnFounder();
  else if (par::filter_nonfounders)
    P.filterOnNonFounder();


  //////////////////////////////// 
  // A covariate file specified?

  if (par::covar_file) 
    {
      // Multiple covariates?
      if (par::clist)
	{
	  if (!P.readCovListFile()) 
	    error("Problem reading the covariates");
	}
      else // a single covariate
	{
	  if (!P.readCovariateFile()) 
	    error("Problem reading the specified covariate from the covariate file");
	}
    }
  


  //////////////////////////////////////
  // Assign cluster solution from file 
  
  if (par::include_cluster_from_file)
    {
      P.printLOG("Reading clusters from [ " + 
		 par::include_cluster_filename+" ]\n");
      if (!P.readClusterFile())
	error("Problem reading from [ "+par::include_cluster_filename+" ]");
    }
  else if ( par::sol_family ) 
    {
      P.printLOG("Setting clusters based on family IDs\n");
      vector<string> famlist;
      P.kname.resize(0);
      for (int i=0; i<P.n; i++)
	{
	  Individual * person = P.sample[i];
	  bool match = false;
	  for (unsigned int j=0; j<famlist.size(); j++)
	    if (person->fid == famlist[j]) 
	      {
  		match=true;
		person->sol=j;
	      }	  
	  if (!match)
	    {
	      famlist.push_back(person->fid);
	      person->sol=famlist.size()-1;
	      P.kname.push_back(person->fid);
	    }
	}

      // Set number of clusters/families
      P.nk = famlist.size();     

      // Set klist variable
      
      P.klist.clear();
      for (int j=0;j<P.nk;j++)
	P.klist.push_back( new Cluster );
      
      for (int i=0; i<P.n; i++)
	if ( P.sample[i]->sol > -1 )
	  P.klist[P.sample[i]->sol]->person.push_back(P.sample[i]);
            
    }
  else
    {
      P.klist.clear();
      P.klist.push_back( new Cluster );
      for (int i=0; i<P.n; i++)
	P.klist[0]->person.push_back(P.sample[i]);      
    }



  /////////////////////////////////////////
  // Zero-out specific sets of genotypes?

  if ( par::zero_cluster )
    P.zeroOnCluster();
  
    
  /////////////////////////////////
  // Fix reference allele? 

  if ( par::set_reference_allele ) 
    P.setReferenceAllele();



  //////////////////////////////////
  // Determine formats for printing

  P.prettyPrintLengths();
  


 
  //////////////////////////////////////////////////
  //                                              //
  // Process a dosage file                        //
  //                                              //
  //////////////////////////////////////////////////

  if ( par::dosage_assoc )
    {
      // Normal behavior is to load data, and perform 
      // analysis; if the hard-call option is specified, 
      // then this will generate a dataset, that we can
      // subsequent filter and save, etc, as usual, i.e.
      // in that case, do not halt

      P.processDosageFile();
      
      if ( ! par::dosage_hard_call )
	shutdown();

    }




  //////////////////////////////////////////////////
  //                                              //
  // Handle CNV segments separately               //
  //                                              //
  //////////////////////////////////////////////////

  
  if ( par::cnv_list )
    {
      P.readCNVList();
      P.processCNVList();
      shutdown();
    }
  


  //////////////////////////////////////////////////
  //                                              //
  // Handle non-SNP data separately               //
  //                                              //
  //////////////////////////////////////////////////
  
  if ( par::gvar || par::gvar_write )
    {
      
      // We might want to load generic variants on top
      // of existing SNP data; or afresh if none of the 
      // above have been specified
      
      if ( ! par::load_gvar )
	P.readGenericVariantData();
      
      if ( par::gvar_write )
	{
	  P.outputGenericVariantFile();
	  shutdown();
	}
      
      P.processGVAR();
      shutdown();
    }
 


  //////////////////////////////////////////////////
  //                                              //
  // Misc. .genome grouper utility                //
  //                                              //
  //////////////////////////////////////////////////
  
  if ( par::genome_groups )    
    {
      P.groupGenome();
      shutdown();
    }


  //////////////////////////////////
  // Missing code

//   if ( par::bt && ! par::missing_phenotype_explicit ) 
//     par::missing_phenotype = "0";
  


  //////////////////////////////////////////////////
  //                                              //
  // Basic MAF, genotyping filters & HWE/ME       //
  //                                              //
  //////////////////////////////////////////////////
  
  
  P.printLOG("Before frequency and genotyping pruning, there are "
	     +int2str(P.nl_all)+" SNPs\n");
  
  if (!par::FIXED_p)
    {
      P.filterSNPs();
    }
  else
    for (int i=0; i<P.nl_all;i++)
      {
	if (P.locus[i]->allele1=="1")
	  P.locus[i]->freq = par::FIX_p;
	else 
	  P.locus[i]->freq = 1-par::FIX_p;
      }
  
  
  P.printLOG("After frequency and genotyping pruning, there are "
	     +int2str(P.nl_all)+" SNPs\n");
  
  if ( P.nl_all == 0 )
    error("Stopping as there are no SNPs left for analysis\n");
  
  if ( P.n == 0 )
    error("Stopping as there are no individuals left for analysis\n");



  //////////////////////////////////////////////////
  // Re-report basic case/control counts, etc

  summaryBasics(P);
  


  //////////////////////////////////////////////////
  // Any null allele codes (monomorhpics)?
  
  for (int l=0; l<P.nl_all; l++)
    {
      if ( P.locus[l]->allele1 == "" )
	P.locus[l]->allele1 = "0";
    }


  /////////////////////////////////////////
  // SET statistics?
  
  if ( par::read_set ) 
    P.readSet();
  else if (par::make_set) 
    P.outputSetFile();
  
  Set S( P.snpset );  
  P.pS = & S;
 
  // Remove any SNPs not in a set
  // unless using particular commands
  // (set-by-all epistasis, set-table)
  
  if ( par::read_set || par::make_set )
    {
      if ( par::drop_sets ) 
	P.pS->dropNotSet(P);
    }
  

  //////////////////////////////////////////////////
  // Build final marker scaffold 

  makeScaffold(P);


  //////////////////////////////////////////////////
  //                                              //
  // Create family units?                         //
  //                                              //
  //////////////////////////////////////////////////


  ///////////////////////
  // Create family units?
  
  if (par::MENDEL_test || 
      par::MENDEL_report || 
      par::TDT_test ||
      par::QTDT_test ||      
      par::make_founders &&
      !par::built_families)
    {
      map<string,Individual*> fnd;
      map<Individual*,int> idmap;
      P.linkRelateds(idmap, fnd);
      P.parseTrios();      
      par::built_families = true;
      
      // Perform now, so that the user has an option to 
      // save a new fileset with mendel errors removed
      
      if (par::MENDEL_report || par::MENDEL_test)
	P.checkMendel();     
     
    }


  ////////////////////////////////////////////////
  // Reset PAT/MAT codes of any non- nonfounders?
  // i.e. if parents not actually present in sample?
  
  if ( par::make_founders )
    {
      P.makeFounders();
    }
  

  //////////////////////////////////////
  // Sex check
 
  if (par::check_sex)
    {
      P.sexCheck();      
    }



  //////////////////////////////////////
  // Split TDT units to case/controls

  if ( par::tucc )
    {
      if ( !par::built_families )
	{
	  map<string,Individual*> fnd;
	  map<Individual*,int> idmap;
	  P.linkRelateds(idmap, fnd);
	  P.parseTrios();      
	  par::built_families = true;
	}

      P.checkMendel();     
      P.pseudoCaseControl();

    }




  //////////////////////////////////////////////////
  //                                              //
  // Haplotype imputation methods                 //
  //                                              //
  //////////////////////////////////////////////////

  // Do not use this old IMPUTATION method
  // Restrict to --proxy-impute, or original
  // --hap-impute (i.e. based on multi-marker list)

  if ( par::meta_large_phase )
    {
      
      // Automatically try to impute all one window per chromosome
      // We can put in some other restraints here if need be
      
      if (par::has_nonfounders && !par::built_families)
	{
	  map<string,Individual*> fnd;
	  map<Individual*,int> idmap;
	  P.linkRelateds(idmap, fnd);
	  P.parseTrios();
	  P.checkMendel();      
	  par::built_families = true;
	}

      P.printLOG("Estimating haplotype frequencies/phases ( MHF >= " 
		 + dbl2str(par::min_hf)+" )\n");
      P.printLOG("Considering phases P(H|G) >= " 
		 +dbl2str(par::hap_min_phase_prob)+"\n");
      P.printLOG("Requiring per individual per haplotype missingness < " 
		 +dbl2str(par::hap_missing_geno)+" \n");

      P.printLOG("Initial EM window size " + int2str(par::haplo_plem_window) 
		 + " SNPs with " + int2str(par::haplo_plem_overlap) + " SNP overlap\n");
      

      // Count number of founders
      
      P.haplo->cnt_f = 0;
      vector<Individual*>::iterator person = P.sample.begin();
      while ( person != P.sample.end() ) 
	{  
	  if ( (*person)->founder ) P.haplo->cnt_f++;	    
	  person++;
	}
      
      if (P.haplo->cnt_f<P.n)
	{
	  P.haplo->nonfounders = true;
	  P.printLOG("Initial phasing based on "+
		     int2str(P.haplo->cnt_f)+" founders ("+
		     int2str(P.n-P.haplo->cnt_f)+
		     " non-founders)\n");
	}
      

      // Start off just with the autosomes
      // We assume that "--chr" has been specified on the command line, 
      // and so we are only dealing with a single chromosome here
            
      if ( par::impute_verbose )
	{
	  P.printLOG("Writing verbose imputation output to [ " 
		     +par::output_file_name + ".phased.out ]\n");     
	  P.haplo->HIMPUTE.open((par::output_file_name+".phased.out").c_str(), ios::out);
	  P.haplo->HIMPUTE.setf(ios::fixed);
	  P.haplo->HIMPUTE.precision(2);
	}
      

      // Run imputation in blocks of up to 1000 SNPs
      
      P.haplo->makeSlidingWindow( "20+20" );
      
      P.haplo->phaseAllHaplotypes(true,perm);

      if ( par::impute_verbose )
	P.haplo->HIMPUTE.close();

    }


  
  ////////////////////////////////////
  // Proxy-based haplotype imputation 
  
  if (par::proxy_impute)
    {
      P.proxyWrapper();
      // Do not shut down: we assume a --make-bed will
      // be called below
    }



  //////////////////////////////////////////////////
  //                                              //
  // Generate dummy permuted phenotype file       //
  //                                              //
  //////////////////////////////////////////////////

  if ( par::output_pheno_perm )
    {
      P.outputPermedPhenotypes(perm);
      shutdown();
    } 


  //////////////////////////////////////////////////
  //                                              //
  // Output formats and transformations           //
  //                                              //
  //////////////////////////////////////////////////


  // Covariate files can also be output (--covar) 
  // for the major options: --make-bed, --recode* 
  // and also just --write-covar option

  if (par::set_table)
    {
      P.setTable();
      shutdown();
    }

  if (par::write_set)
    {
      P.writeSetFile();
      shutdown();
    }

  if (par::dump_covar)
    {
      P.write_covariates();
      shutdown();
    }

  if (par::dump_clst)
    {
      P.write_clusters();
      shutdown();
    }

  if (par::write_snplist)
    {
      P.write_snplist();
      shutdown();
    }
  
  if (par::write_bitfile) 
    {
      P.write_BITFILE();
      if (par::clist)
	P.write_covariates();
      shutdown();
    }

  if ( par::recode_fastphase )
  {
    P.output_fastphase_format();
    shutdown();
  }

  if ( par::recode_bimbam )
  {
    P.output_bimbam_format();
    shutdown();
  }

  if ( par::recode_structure )
  {
    P.output_structure_format();
    shutdown();
  }


  if (par::recode || par::recode_HV || par::recode_12 || par::recode_whap ) 
    {
      if ( ! par::recode_transpose )
	P.display_recoded_PEDFILE();
      else
	P.display_recoded_PEDFILE_transpose();
      if (par::clist)
	P.write_covariates();
      shutdown();
    }

  if (par::recode_AD) 
    {
      P.display_recoded_PEDFILE_AD();
      if (par::clist)
	P.write_covariates();
      shutdown();
    }

  if (par::recode_long) 
    {
      P.display_recoded_LONG();
      if (par::clist)
	P.write_covariates();
      shutdown();
    }

  if (par::recode_mutlist) 
    {
      P.display_recoded_MUTLIST();
      if (par::clist)
	P.write_covariates();
      shutdown();
    }
  

  if (par::list_by_allele)
    {
      P.display_listByAllele();
      shutdown();
    }

  if (par::plist)
    {
      P.display_pairList();
    }

  if (par::indiv_report)
    {
      P.display_indivReport();
      shutdown();
    }

  if (par::list_twolocus)
    {
      P.display_twolocus();
      shutdown();
    }
  

  //////////////////////////////////////////////////
  //                                              //
  // LD-based lookups                             //
  //                                              //
  //////////////////////////////////////////////////
  
  

  ////////////////////////////////////////////
  // Set summary statistics

  if (par::set_screen) 
    {
      P.setAssocSummary();
      shutdown();
    }


  ////////////////////////////////////////////
  // LD-based clumping

  if ( par::clumpld )
    {
      clump_LD cld(&P,P.haplo,
		   par::clumpld_p1, 
		   par::clumpld_kb, 
		   par::clumpld_p2, 
		   par::clumpld_r2);
      cld.clump();
      shutdown();
    }


  ////////////////////////////////////////////
  // Show tags

  if ( par::gettag_mode )
    {
      P.tagMode();
      shutdown();
    }
  

  ////////////////////////////////////////////
  // Haplotype block action

  if ( par::make_blocks )
    {
      P.mkBlks(0,P.nl_all-1);
      shutdown();
    }


  //////////////////////////////////////////////////
  //                                              //
  // Main set of whole-genome tests               //
  //                                              //
  //////////////////////////////////////////////////
  


  /////////////////////////////////
  // Some initial set-up work here


  /////////////////////
  // Conditioning SNPs

  if (par::conditioning_snps)
    {
      if (par::conditioning_snp_single)
	{

	  // ** todo ** change this to allow a NList

	  int x = getMarkerNumber( P, par::conditioning_snp_name );

	  if (x<0) error("Marker "
			 +par::conditioning_snp_name
			 +" does not exist in filtered data\n");

	  P.conditioner.push_back( x );
	  P.conditioner_mask.push_back( false );
	}
      else
	P.readConditioningList();
    }


  //////////////////////////////////////////
  // Warn if not enough markers in analysis

  if ( par::plink 
       || par::cluster 
       || par::cluster_plot 
       || par::outlier_detection
       || par::genome_output 
       || par::inbreeding )
    {
      if (P.nl_all < 10000) 
	P.printLOG("\n **Warning** this analysis typically requires whole-genome level data\n"
		   "             to give accurate results \n\n");
    }



  //////////////////////////////////////////
  // Arbitrary external functions

  if (par::myfunction)
    {
      if (1)
	{
	  if (par::has_nonfounders && !par::built_families)
	    {
	      map<string,Individual*> fnd;
	      map<Individual*,int> idmap;
	      P.linkRelateds(idmap, fnd);
	      P.parseTrios();
	      P.checkMendel();      
	      par::built_families = true;
	    }        	  
	}

      // P.callMe();
      shutdown();
    }



  //////////////////////////////////////////////////
  //                                              //
  // IBS and IBD genome-wide analyses             //
  //                                              //
  //////////////////////////////////////////////////


  //////////////////////////////////////////////
  // Perform a cluster analysis and/or MDS plot 

  if (par::cluster || par::cluster_plot || par::outlier_detection ) 
    { 
      P.buildCluster();
      shutdown();
    }
  

  ///////////////////////////////////
  // Permutation test between groups
  // based on IBS diffeences
  
  if (par::ibs_test) 
  { 
      P.permutationIBSTest(perm);
      shutdown();
  }

  
  ////////////////////////////////////////////////////
  // Precalculate frequency-averaged P(IBD|IBS) table

  if (par::plink || par::genome_output) 
    {
      if (par::has_nonfounders && !par::built_families)
	{
	  map<string,Individual*> fnd;
	  map<Individual*,int> idmap;
	  P.linkRelateds(idmap, fnd);
	  P.parseTrios();
	  // P.checkMendel(); // skip this when in --rel-check mode
	  par::built_families = true;
	}      

      
      // So that correct IBD expectation is calculated, 
      // we need to fill in empty slots for missing parents
      
      P.makeMissingParents();        

      P.preCalcGenomeIBD();
    }
  
  //////////////////////////////
  // Genome-wide output only
 
  if (par::genome_output)
    {
      P.displayGenomeWideInfo();
      
      if ( par::genome_test ) 
	P.testGenomeIBDByCovariate(perm);

      shutdown();
    }
  
  
  //////////////////////////////////////
  // Genome-wide inbreeding output only
 
  if (par::inbreeding)
    {
      
      if (par::SNP_major) 
	P.SNP2Ind();

      ofstream HET;
      string f = par::output_file_name + ".het";
      HET.open(f.c_str(),ios::out);
      HET.precision(4);
        
      P.printLOG("Writing individual heterozygosity information to [ "+f+" ] \n");
      HET << setw(par::pp_maxfid) << "FID" << " "
	  << setw(par::pp_maxiid) << "IID" << " "
	  << setw(12) << "O(HOM)" << " " 
	  << setw(12) << "E(HOM)" << " " 
	  << setw(12) << "N(NM)" << " " 
	  << setw(12) << "F" << "\n";
      
      for (int i1=0; i1<P.n; i1++)
	P.calcInbreeding(P.sample[i1],0,P.nl_all-1,HET);
      
      HET.close();
      
  }



  ///////////////////////
  // Runs of homozygosity
  
  if (par::homo_run)
    {
      
      P.findAllHomozygousRuns(perm);
      
      if ( par::segment_test_individual )
	P.segmentIndividualTest(perm);
      
      shutdown();
    }


  ///////////////////////////////////
  // Runs of missingness (deletions)

   if (par::miss_run)
     {

       if (par::SNP_major) P.SNP2Ind();
       
       ofstream RUN;
       string f = par::output_file_name + ".rum";
       RUN.open(f.c_str(),ios::out);
    
       P.printLOG("Writing run-of-missings information to [ "+f+" ] \n");
       
       string msg = "Run defined as " + int2str(par::miss_run_length);
       if (par::miss_run_length_kb) msg += " kb\n";
       else msg += " SNPs\n";
       P.printLOG(msg);

       stringstream s2;
       s2 << "With at least " << par::miss_run_level << " missingness\n";
       P.printLOG(s2.str());

       for (int i1=0; i1<P.n; i1++)
 	{
 	  if (!par::silent)
 	    cout << i1+1 << " of " << P.n << " individuals      \r";
 	  P.findMissRuns(P.sample[i1],RUN);
 	}

       if (!par::silent)
 	cout << "\n\n";
    
       RUN.close();

       shutdown();
   }




  //////////////////////////////////////////////////
  //                                              //
  // LD and haplotype-based analyses              //
  //                                              //
  //////////////////////////////////////////////////


  ///////////////////
  // LD-based pruning
  
  if (par::prune_ld)
    {  
      P.pruneLD();
      shutdown();
    }


  //////////////////////////////
  // Flip-scan

  if ( par::flip_scan )
    {
      P.calcFlipScan();
      shutdown();
    }

  //////////////////////////////
  // LD statistics

  if (par::calc_SNPSNP_LD)
    {
      P.calcPairwiseLD();
      shutdown();
    }

  if (par::disp_r1 || par::disp_r2)
    {
      P.calcLDStatistics();
      shutdown();
    }
 

  ///////////////////////////////////////////////////
  // General class for haplotype phasing and testing
  
  if ( par::test_hap_CC && par::qt )
    {
      par::test_hap_CC = false;
      par::test_hap_QTL = true;
    }

  // In case families are included, build family structure if not
  // already done
  
  if ( par::phase_snps || par::mishap_test || par::proxy_assoc) 
    {

      // Read in list of tests, or make sliding window?
      
      if (par::phase_snps) 
	{
	  if (par::sliding_window)
	    P.haplo->makeSlidingWindow(par::sliding_window_size);
	  else if (par::hap_specific_snps)
	    P.haplo->setSpecificSNPs(par::hap_specific_snps_list);
	  else
	    P.haplo->readTagFile();
	}
      
      if (par::has_nonfounders && !par::built_families)
	{
	  map<string,Individual*> fnd;
	  map<Individual*,int> idmap;
	  P.linkRelateds(idmap, fnd);
	  P.parseTrios();
	  P.checkMendel();      
	  par::built_families = true;
	}

      P.printLOG("Estimating haplotype frequencies/phases ( MHF >= " 
		 + dbl2str(par::min_hf)+" )\n");
      P.printLOG("Considering phases P(H|G) >= " 
		 +dbl2str(par::hap_min_phase_prob)+"\n");
      P.printLOG("Requiring per individual per haplotype missingness < " 
		 +dbl2str(par::hap_missing_geno)+" \n");

      // Count number of founders
      
      P.haplo->cnt_f = 0;
      vector<Individual*>::iterator person = P.sample.begin();
      while ( person != P.sample.end() ) 
	{  
	  if ( (*person)->founder ) P.haplo->cnt_f++;	    
	  person++;
	}

      if (P.haplo->cnt_f<P.n)
	{
	  
	  P.haplo->nonfounders = true;
	  P.printLOG("Initial phasing based on "+
		     int2str(P.haplo->cnt_f)+" founders ("+
		     int2str(P.n-P.haplo->cnt_f)+
		     " non-founders)\n");
	}
      
      if (P.n == P.haplo->cnt_f && 
	  ( par::test_hap_TDT || par::proxy_TDT ) )
	error("Can not perform TDT in sample with no non-founders");
      
    }
  
  
  /////////////////////////
  // Haplotype frequencies
  
  if (par::phase_snps && par::display_hap_freqs)
    {
      P.haplo->calculateHaplotypeFrequencies();
      shutdown();
    }

  
  ////////////////////////////////
  // Haplotype phase probabilities
  
  if (par::phase_snps && par::display_phase_probs)
    {
      P.haplo->calculateHaplotypeFrequencies();
      shutdown();
    }
    



  /////////////////////////////////////////////
  // Haplotypic test of non-random missing data
  
  if (par::mishap_test)
    {
      P.performMisHapTests();
      shutdown();
    }
  


  ////////////////////////////////////////////////////
  // Haplotype tracking of an extended region, for an 
  // individual or pair
  
  if (par::phase_snps && par::segment_haplotrack)
    {
      P.haplo->trackSharedHaplotypes();
      shutdown();
    }



  ////////////////////////////////////////////////////
  // Haplotype tracking of an extended region, for an 
  // individual or pair
  
  if (par::phase_snps && par::impute_tags )
    {
      P.haplo->imputeAllHaplotypes();
      shutdown();
    }


  ///////////////////////////////////////////////////////
  // Haplotypic test of SNP proxy (convenience function)
  // (we've already done the imputation step above)

  if (par::proxy_assoc && ! par::proxy_impute)
    {
      P.proxyWrapper();
      shutdown();
    }



  //////////////////////////////////////////////////
  //                                              //
  // Misc tests that do not fall within the       //
  // main phenotype loop                          //
  //                                              //
  //////////////////////////////////////////////////


  //////////////////////////////
  // Genome-wide IBS sharing test

  if (par::ibs_sharing_test)
    {
      P.perm_sharingIBSTest(perm);
      shutdown();
    }


  //////////////////////////////
  // Gene-based test of epistasis

  if (par::epi_genebased)
    {
      P.driverSCREEPI();
      shutdown();
    }


  //////////////////////////////
  // Genome-wide epistasis tests

  if (par::epistasis)
    {
      P.calcEpistasis();      
      shutdown();
    }



  /////////////////////////////////////////////
  // Determine per-individual risk profiles
  
  if (par::score_risk)
    {
      P.scoreIndividuals();
      shutdown();
    }


  //////////////////////////////////
  // Apply an R-script to the data?

  if ( par::run_R_script )
    {
     
#ifdef WITH_R_PLUGINS
      P.Rfunc();
      shutdown();
#else
      error("R plugin support has not been compiled in");
#endif  
    }



  //////////////////////////////////////////////////
  //                                              //
  // Genome-wide association tests                //
  //                                              //
  //////////////////////////////////////////////////
  
  // Allow for the fact that we might be iterating 
  // over multiple phenotypes

  string original_file_root = par::output_file_name;

  
  if ( ! par::plink ) 
    while ( 1 )
      {
	
	if ( par::all_pheno )
	  {
	    if ( par::loop_over )
	      {
		P.phenoLabel = P.kname[ par::loop_counter ];
		par::output_file_name = original_file_root + "." + P.phenoLabel;
		par::bt = true; par::qt = false;

		for (int i=0; i<P.n; i++)
		  {
		    // Include all samples
		    P.sample[i]->missing = false;
		    P.sample[i]->aff = P.sample[i]->sol == par::loop_counter ? true : false ;		
		  }
	      }
	    else
	      {
		if ( P.phenotype_name == "" )
		  P.phenoLabel = "P"+int2str(par::mult_pheno);
		else
		  P.phenoLabel = P.phenotype_name;

		par::output_file_name = original_file_root + "." + P.phenoLabel;
	      }
	  }
	
	if (par::assoc_test)
	  {
	    
	    if (par::CMH_test_2)
	      P.calcMH();
	    else if (par::OR_homog_test)
	      P.calcHomog();
	    else if (par::QTDT_test)
	      {
		// Force a Mendel error check
		if (! (par::MENDEL_report || par::MENDEL_test) )
		  P.checkMendel();
		P.perm_testQTDT(perm);
	      }
	    else if (par::boot)
	      {
		// Redundant
		error("Bootstrap option is no longer supported\n");
		
		P.calcAssociationWithBootstrap();
	      }
	    else
	      {
		
		// Includes 
		//        basic allelic test
		//        model-based tests
		//        linear & logistic models
		//        2x2xK Cochran-Mantel-Haenszel
		
		P.calcAssociationWithPermutation(perm);
	      }
	    
	    if ( ! par::all_pheno )
	      shutdown();
	    
	  }
	
	
	/////////////////////////////////
	// Haplotype association analysis
	
	if (par::phase_snps && ( par::test_hap_CC ||
				 par::test_hap_GLM ||
				 par::test_hap_QTL ||
				 par::test_hap_TDT ) )
	  {
	    
	    // This is done separaytely, via the main
	    // assoc. loop

	    if (par::test_hap_GLM)
	      P.calcAssociationWithPermutation(perm);
	    else
	      {

		////////////////////////////////////////////////
		// Perform omnibus and haplotype-specific tests
		
		string f;
		
		if ( par::test_hap_CC ) 
		  f =  par::output_file_name + ".assoc.hap";
		else if ( par::test_hap_QTL ) 
		  f =  par::output_file_name + ".qassoc.hap";
		else if ( par::test_hap_TDT ) 
		  f = par::output_file_name + ".tdt.hap";
		
		if (par::test_hap_CC)
		  {
		    
		    P.printLOG("Writing haplotype association statistics to [ " + f + " ]\n");
		    P.haplo->HTEST.open(f.c_str(), ios::out);
		    P.haplo->HTEST.precision(4);
		    
		    P.haplo->HTEST << setw(10) << "LOCUS" << " " 
				   << setw(12) << "HAPLOTYPE" << " "
				   << setw(10) << "F_A" << " "
				   << setw(10) << "F_U" << " " 
				   << setw(10) << "CHISQ" << " "
				   << setw(4) << "DF" << " "
				   << setw(10) << "P" << " "
				   << "SNPS" << "\n";
		  }
		
		
		if ( par::test_hap_QTL )
		  {
		    P.printLOG("Writing haplotype association statistics to [ " + f + " ]\n");
		    P.haplo->HTEST.open(f.c_str(), ios::out);
		    P.haplo->HTEST.precision(4);
		    
		    P.haplo->HTEST << setw(10) << "LOCUS" << " " 
				   << setw(12) << "HAPLOTYPE" << " "
				   << setw(8) << "NANAL" << " "
				   << setw(10) << "BETA" << " "
				   << setw(10) << "R2" << " " 
				   << setw(8) << "STAT" << " "
				   << setw(10) << "P" << " "
				   << "SNPS" << "\n";		  
		    
		  }
		
		if (par::test_hap_TDT)
		  {
		    P.printLOG("Writing haplotype TDT statistics to [ " + f + " ]\n");
		    P.haplo->HTEST.open(f.c_str(), ios::out);
		    P.haplo->HTEST.precision(4);
		    
		    P.haplo->HTEST << setw(10) << "LOCUS" << " " 
				   << setw(12) << "HAPLOTYPE" << " "
				   << setw(10) << "T" << " "
				   << setw(10) << "U" << " " 
				   << setw(10) << "CHISQ" << " "
				   << setw(10) << "P" << " "
				   << "SNPS" << "\n";
		  }
	    
		
		P.haplo->phaseAllHaplotypes(true,perm);
	    
		P.haplo->HTEST.close();
	      }

	    if ( ! par::all_pheno )
	      shutdown();
	    
	  }
	

	//////////////////////////////////////////////////////
	// Haplotypic conditional tests (WHAP implementation, 
	// now called CHAP, for conditional haplotype
	
	if (par::chap_test)
	  {
	    
	    P.conditionalHaplotypeTest(true,perm);
	    
	    if ( ! par::all_pheno )
	      shutdown();
	    
	  }


	
	//////////////////////////////
	// QTL interaction test
      
	if (par::assoc_gxe)
	  {
	    P.perm_testGXE2(perm);
	    if ( ! par::all_pheno )
	      shutdown();
	  }
	
	
	
	/////////////////////////////
	// Rare allele test 
	
	if ( par::elf_baseline )
	  {
	    P.elfBaseline();
	    shutdown();
	  }
	
	if (par::rare_test)
	  {
	    P.permTestRareDistribution(perm);
	    if ( ! par::all_pheno )
	      shutdown();
	  }
	
	
	/////////////////////////
	// Hotelling's T^2 test
	
	if (par::hotel)
	  {
	    P.perm_testHotel(perm);
	    if ( ! par::all_pheno )
	      shutdown();
	  }
	
	
	///////////////////////////////////
	// Test difference in missing rates
	
	if (par::test_missing)
	  {	
	    P.calcAssociationWithPermutation(perm);
	    if ( ! par::all_pheno )
	      shutdown();
	  }
	
	
	//////////////////////////////////
	// Genome-wide family-based (TDT)
	// and Parent-of-origin analysis
	
	if (par::TDT_test)
	  {
	  
	    // Force a Mendel error check, if we have not 
	    // already
	    
	    if (! (par::MENDEL_report || par::MENDEL_test) )
	      P.checkMendel();
	    
	    // Either basic TDT or Parent-Of-Origin analysis
	    
	    if (par::parent_of_origin)
	      P.perm_testTDT_POO(perm);
	    else if (par::sibTDT_test)	  
	      P.perm_testTDT(perm);
	    else
	      P.perm_testTDT(perm);
	    
	    if ( ! par::all_pheno )
	      shutdown();
	    
	  }
	
	
	// Read next phenotype: repeat, or shutdown
	
	if ( par::all_pheno )
	  {
	    
	    if ( par::loop_over )
	      {
		// Construct next phenotype from cluster file
		
		par::loop_counter++;
		if ( par::loop_counter == P.nk ) 
		  shutdown();
		
	      }
	    else
	      {
		// Read next phenotype from file
		
		par::mult_pheno++;
		if ( ! P.readPhenoFile() )
		  shutdown();

		// and recode, if a binary affection status coding

		if (par::bt)
		  affCoding(P);

	      }
	  }
	
	if ( ! par::all_pheno )
	  shutdown();
      
      } // Next potential phenotype
  
  
  

  //////////////////////////////////////////////////
  //                                              //
  // PLINK segmental sharing analyses             //
  //                                              //
  //////////////////////////////////////////////////

  // Stop now, unless a plink analysis is specified
  if (!par::plink)
    shutdown();

  if (par::SNP_major) 
    P.SNP2Ind();



  //////////////////////////////////////////////
  // Read pre-computed segment list and perform 
  // segmental tests?

  if (par::read_segment_file)
    {

      ifstream SEG;
      SEG.open(par::read_segment_filename.c_str(),ios::in);
      P.printLOG("Reading IBD-segment information from [ "
		 +par::read_segment_filename+" ]\n");
      checkFileExists(par::read_segment_filename);

      if (par::segment_minimal)
	P.readSegmentFileMinimal(SEG);
      else
	P.readSegmentFile(SEG);
      SEG.close();      
      
      // IBS validation of segments (i.e. possibly in a larger
      // datafile? but one that must be a superset of all SNPs in
      // segment file)
      
      if (false)
	{
	  P.validateSegments();
	  shutdown();
	}
      
      // Find overlap in segments?
      if (par::segment_overlap)
	P.summariseHomoRuns();  
      
      // Per-individual summary/test?
      if ( par::segment_test_individual )
	{
	  P.segmentIndividualTest(perm);
	  shutdown();
	}
      
      // Perform pairwise summary/analysis of segments?

      
      P.summaryIBSsegments(perm);
      

      P.printLOG("Writing segment summary to [ " + par::output_file_name 
	       + ".segment.indiv ]\n\n");

      P.indivSegmentSummary();

      shutdown();
    }


  
  //////////////////////////////
  // Pair inclusion/exclusion
  
  // Number of informative pairs
  int c=0; 
  
  // Read or calculate informative pairs?
  if (par::ibd_read) 
    c = P.readInformative();
  else
    c = P.calcInformative();


  ////////////////////////////////////////////////
  // Test of genome-wide relatedness by covariate

  if ( par::genome_test ) 
    {
      P.testGenomeIBDByCovariate(perm);
      shutdown();
    }
  
  
  ///////////////////////////////////
  // Save pairs to be included? i.e. 
  // after removing all pairs for 

  // a) low IBD
  // b) not being an affected pair
  // c) being a concordant unaffected pair

  if (par::inc_write) 
    P.writeInformative();


  /////////////////////////////////////////////////////////////////
  // Get and display information on chromosomal range to be tested
  
  //  else if (par::singlepoint)
  //    P.printLOG("Using singlepoint analysis mode\n");
  else if (par::inter_grid>0)
    {
      stringstream s2;
      s2 << "Using multipoint analysis: step = " 
	 << par::inter_grid 
	 << " and fringe = " 
	 << par::fringe << " cM\n";
      P.printLOG(s2.str());
    }
  else
    {
      stringstream s2;
      s2 << "Using multipoint analysis: grid = " 
	 << par::grid 
	 << " and fringe = " 
	 << par::fringe << " cM\n";
      P.printLOG(s2.str());
    }

  vector<int> chrs;
  if (par::run_chr==0) {

    vector<int> r = getChromosomeRange(P);
    
    P.printLOG("\nScanning from autosomes from chromosome "+
	       chromosomeName( r[0] ) + " to "+
	       chromosomeName( r[1] ) + "\n\n");
    for (int i=r[0];i<=r[1];i++) 
      if ( ( !par::chr_haploid[i] ) && 
	   ( !par::chr_sex[i] ) )
	chrs.push_back(i);
  } else chrs.push_back(par::run_chr);
  

  ofstream SEG;
  if (par::segment_output)
    {
      string f = par::output_file_name + ".segment";
      SEG.open(f.c_str(),ios::out);
      P.printLOG("Writing IBD-segment information to [ "+f+" ]\n");
      if (par::segment_minimal) P.printLOG("Minimal segment file format\n");

      // Header row for non-minimal format
      if (! par::segment_minimal)
	{
	  SEG << setw(par::pp_maxfid) << "FID1" << " "
	      << setw(par::pp_maxiid) << "IID1" << " "
	      << setw(par::pp_maxfid) << "FID2" << " "
	      << setw(par::pp_maxiid) << "IID2" << " ";
	  
	  if (par::bt) SEG << setw(4) << "PHE" << " ";
	  
	  SEG << setw(4) << "CHR" << " "
	      << setw(10) << "BP1" << " "
	      << setw(10) << "BP2" << " " 
	      << setw(par::pp_maxsnp) << "SNP1" << " "
	      << setw(par::pp_maxsnp) << "SNP2" << " "
	      << setw(6) << "NSNP" << " "
	      << setw(10) << "KB" << "\n";
	}
      
      f = par::output_file_name + ".segment.summary";
      P.printLOG("Writing IBD-segment summary to [ "+f+" ]\n\n");
      
      P.printLOG("Minimum segment length is "
		 +dbl2str((double)par::segment_length/(double)1000)
		 +" kb and "+int2str(par::segment_snp)+" SNPs\n");
      P.printLOG("Segment thresholds are "+dbl2str(par::segment_threshold_start)
		 +" and "+dbl2str(par::segment_threshold_finish)+"\n");
      P.printLOG("Maximum intra-segment inter-SNP distance is "
		 +int2str(par::segment_inter_snp_distance)
		 +"\n");


    }


  ofstream MP;
  if (par::multi_output)
    {
      string f = par::output_file_name + ".multi";
      MP.open(f.c_str(), ios::out);
      MP.setf(ios::fixed);
      MP.precision(5);
      P.printLOG("Writing multipoint IBD estimates to [ "+ f+" ]\n");
    }

  ofstream GMULTI;
  if (par::gmulti_output)
    {
      string f = par::output_file_name + ".gmulti";
      GMULTI.open(f.c_str(), ios::out);
      GMULTI.precision(4);
      P.printLOG("Writing genotype/multipoint IBD estimates to [ "+ f +" ]\n");
    }



  //////////////////////////////
  // Consider each chromosome

  for (int i=0;i<chrs.size();i++)
    {

      //////////////////////////
      // Reset main variables

      P.phenotype.resize(0);
      P.pihat.resize(0);
      P.Zlocus.resize(0);
      P.m_pihat.resize(0);
      P.v_pihat.resize(0);
      P.pair1.resize(0);
      P.pair2.resize(0);
      
      // Set chromosome
      par::run_chr = chrs[i];
      
      // Find scan range
      P.setMarkerRange();
      
      // Total number of all (test+background) loci
      P.nl_all = P.locus.size();
      
      // Number of loci for test loci
      P.nl = par::run_end - par::run_start + 1; 
      P.printLOG(int2str(P.nl)+" markers in this scan\n");
      
      // Set up (single) multipoint marker map
      //  if (par::singlepoint)
      //   P.preCalcSinglePoint();

      P.preCalcMultiPoint();

           
      /////////////////////////////////
      //  For each pair of individuals
      //  calculate mutlipoint pihats
      
      
      // reset counter
      int c1=0;
      int c2=0;
      
      for (int i1=0; i1<P.n-1; i1++)
	for (int i2=i1+1; i2<P.n; i2++)
	  {
	    
	    if (!P.skip_pair[c1++])
	      {
		Individual * p1 = P.sample[i1];
		Individual * p2 = P.sample[i2];
		

		/////////////////////////////////
		// Skip if genome sets specified

		if ( par::genome_2sets )
		  {	    
		    if ( ! ( (P.gset1.find(p1) != P.gset1.end() && 
			      P.gset2.find(p2) != P.gset2.end() ) ||
			     (P.gset1.find(p2) != P.gset1.end() && 
			      P.gset2.find(p1) != P.gset2.end() ) ) )
		      continue;
		  }


		/////////////////////////////////
		// Skip if within-cluster analysis
		// specified

		if ( par::IBD_within )
		  {
		    if ( p1->sol != p2->sol )
		      continue;
		  }
	      
		
		/////////////////////////////////
		//  1. Calculate IBD(g) | IBS(g) 
		  
	        Z IBDg = P.saved_IBDg[c2];
		
		if (!par::silent)
		  {
		    cout << "IBD calculation: " 
			 << ++c2
			 << " of " 
			 << c 
			 << "                  \r";
		    cout.flush();
		  }
		
		/////////////////////////////////
		//  2. Calculate IBD(l) - IBD(g)
		
		vector<Z> IBDl = P.calcLocusIBD(p1,p2,IBDg);
		
		
		/////////////////////////////////
		//  3. Multipoint calculation
				
		P.pairid = itoa((int)p1->phenotype,10)
		  +" "+itoa((int)p2->phenotype,10)+" ";
		P.pairid += itoa((int)c2,10) + " ";
		P.pairid += p1->fid+"_"+p1->iid+"_ ";
		P.pairid += p2->fid+"_"+p2->iid+"_";
		
		vector_t p;
		

		//////////////////////////////////
		// Perform either using
		//     Singlepoint analysis
		//     Multipoint analysis (default)
		
		// if (par::singlepoint) 
		//  p = P.calcSinglePoint(IBDl,IBDg);
		
		p = P.calcMultiPoint(IBDl,IBDg,MP);

		

		////////////////////////////////
		//   3b. Verbose output:
		//       genotypes for each pair
		
		if (par::gmulti_output)
		  {
		    for (int l=par::run_start; l<=par::run_end; l++)
		      P.displayGMULTI(p1,p2,l,GMULTI);			
		  }
		
				

		///////////////////////////////
		//  4. Scan for segments of IBD
		  
		P.findSegments(i1,i2,p,SEG);


		/////////////////////////////
		//  5. Add to list 

		// Do not bother saving for now...
		// only save segments...
		if (false) P.pihat.push_back(p);
		
		// And (A,B) pair to list
		P.pair1.push_back(i1);
		P.pair2.push_back(i2);
		
	      }
  	    	    
	  }
      
      if (!par::silent)
	cout << "\n";
      
      
      /////////////////////////////////
      //  Make list of unique individuals
      
      // copy first set of individuals
      P.in_anal = P.pair1;
      for (unsigned int ind=0; ind<P.pair2.size(); ind++) 
	P.in_anal.push_back(P.pair2[ind]);
      sort(P.in_anal.begin(),P.in_anal.end());
      vector<int>::iterator new_end=
	unique(P.in_anal.begin(), P.in_anal.end());
      // delete all elements past new_end 
      P.in_anal.erase(new_end, P.in_anal.end());
      
      P.printLOG(int2str(P.in_anal.size())+
		 " unique, informative individuals in analysis\n");
      
      if ( P.in_anal.size() == 0 )
	{
	  error("No individuals left in analysis: halting");
	}
      
      /////////////////////////////////
      //  Verbose output: summarise IBD
	
      if (par::segment_output)
	{
	  // P.summaryIBDsegments(perm);
	}
      else
	if (par::summary_ibd_output)
	  P.summaryIBD();
      
	
      /////////////////////////////////
      //  Next chromosome
      
      par::done_global_pihat = true;
      
      if (!par::silent)
	cout << "\n";
      
    }
  

  // Now do IBD segment (as IBS...)

  if (par::segment_output)
    {
      P.summaryIBSsegments(perm);
      P.indivSegmentSummary();
    }
  

  //////////////////////////////
  // Find overlap in segments?

  if (par::segment_overlap)
    P.summariseHomoRuns();  


  //////////////////////////////////////
  // Shut segment and multipoint files

  if (par::segment_output)
    SEG.close();
  if (par::multi_output)
    MP.close();
  if (par::gmulti_output)
    GMULTI.close();
  
  
  ////////////////////////////////
  //  Output genome-wide p-values

  if (par::permute)
    if (chrs.size()>=1 && (!par::ignore_phenotypes))
      P.displayGenomePV();
  
    

  ////////////////////////////////
  //  We're definitely done now
  
  shutdown();

}


////////////////////////////////
// Clean-up

void Plink::cleanUp()
{

  for (int i=0; i<sample.size(); i++)
    delete sample[i];
  
  for (int l=0; l<locus.size(); l++)
    delete locus[l];      
  
  if (par::SNP_major)
    for (int l=0; l<SNP.size(); l++)
      delete SNP[l];      

}

