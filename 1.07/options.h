

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distibuted under the GNU General Public         //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


#ifndef __OPTIONS_H__
#define __OPTIONS_H__

#include <string>
#include <vector>
#include "plink.h"

using namespace std;

/*   static Options lookup2_options; */
//   static Options idhelp_replace_options;
//  static Options idhelp_match_options;

//   static Options annot_options;
//      static Options dosage_options;
//


class OptionSet {
  
 public:
  
  map<string,vector<string> > val;

  bool isSet(string s)
    {
      return val.find(s) != val.end();
    }
  
  vector<string> getValues(string s)
    {
      vector<string> sv;
      map<string,vector<string> >::iterator i = val.find(s);
      if ( i == val.end() )
	return sv;
      return i->second;
    }
  
  string getValue(string s)
    {
      
      map<string,vector<string> >::iterator i = val.find(s);
      if ( i == val.end() )
	return "";
      if ( i->second.size() > 0 )
	return i->second[0];
      else 
	return "";
    }

  void display()
    {
      map<string,vector<string> >::iterator i = val.begin();
      while ( i != val.end() )
	{
	  cout << i->first;
	  if ( i->second.size() > 0 )
	    {
	      cout << " : ";
	      for ( int k = 0 ; k < i->second.size(); k++)
		{
		  cout << " " << i->second[k]  ;
		}
	    }
	  cout << "\n";
	  ++i;
	}
    }
};


class Options {
  
  map<string,OptionSet*> opt;
  
 public:
  
  OptionSet * addOption(string s)
    {
      map<string,OptionSet*>::iterator i = opt.find(s);
      if ( i != opt.end() )
	return i->second;
      OptionSet * o = new OptionSet;
      opt.insert(make_pair(s,o));
      return o;
    }
  
  ~Options()
    {
      map<string,OptionSet*>::iterator i = opt.begin();
      while (i != opt.end() )
	{
	  delete i->second;
	  ++i;
	}
    }
  
  OptionSet * getOptions(string s)
    {
      map<string,OptionSet*>::iterator i = opt.find(s);
      if ( i != opt.end() )
	return i->second;
      else
	{
	  OptionSet * o = new OptionSet;
	  return o;	  
	}
    }
};

  
class par {
  
 public:

  static bool myfunction;
  
  static Options opt;

  static bool verbose;
  static bool flag;
  static bool dumpped;
  static bool debug;
  static bool dummy;
  static int dummy_nind;
  static int dummy_nsnp;
  static bool web_check;
  static bool tucc;
  static bool do_not_load_snps;
  
  static const double epsilon;
  static long unsigned int random_seed;

  static int simul_ncases;
  static int simul_ncontrols;
  static string simul_label;
  static double simul_prevalence;
  static bool simul;
  static string simul_file;
  static bool simul_tags;
  static bool simul_haps;
  static bool simul_qt;
  static double simul_qt_var;

  static bool lookup;
  static bool lookup_single_snp;
  static bool lookup_to_file;
  static string lookup_snp;
  static string lookup_gene_name;
  static bool lookup_gene;
  static bool lookup_multiple_genes;
  static int lookup_gene_kb_window;
  static int lookup_snp_kb_window;

  static bool lookup2;
  static string lookup2_cmd;


  static bool idhelp;
  static string idhelp_output_delimit;
  static string idhelp_dictionary;
  static bool idhelp_dump_from_dict;
  static string idhelp_dump_from_dict_cmd;
  static bool idhelp_auto_alias;
  static bool idhelp_lookup;
  static string idhelp_lookup_string;
  static bool idhelp_subset;
  static string idhelp_subset_string;
  static bool idhelp_replace;
  static string idhelp_replace_string;
  static bool idhelp_match;
  static vector<string> idhelp_match_string;

  static bool idhelp_no_dict;
  static bool idhelp_list_aliases;
  static bool idhelp_alias_update;

  static string idhelp_command;
  static string idhelp_input;

  static bool run_R_script;
  static bool run_R_write_script;
  static string R_script;
  static bool run_R_chisq;
  static bool run_R_z;
  static int run_R_nsnps;
  static int R_port;

  static bool recode;
  static bool recode_transpose;
  static bool recode_long;
  static bool recode_long_ref;
  static bool recode_mutlist;
  static bool recode_12;
  static bool recode_AD;
  static bool recode_AD_Aonly;
  static bool recode_AD_fixed;
  static bool recode_allele_coding;
  static string recode_allele_coding_file;

  static bool recode_1234;
  static bool recode_ACGT;
  
  static bool set_reference_allele;
  static string set_reference_allele_file;
  static bool lfile_allele_count;

  static string recode_delimit;
  static string recode_indelimit;
  static bool recode_HV;
  static bool recode_whap;
  static bool recode_fastphase;
  static bool recode_structure;
  static bool recode_bimbam;
  static bool preserve_all_genotypes;
  static bool preserve_mendel_errors;
  static bool zero_cluster;
  static string zero_cluster_filename;
  static bool oblig_missing;
  static string oblig_missing_filename;
  static string oblig_clusters_filename;
  static bool loop_over;
  static string loop_over_label;
  static int loop_counter;
  static string loop_over_filename;
  static bool list_by_allele; 
  static bool list_twolocus;
  static string twolocus_snp1;
  static string twolocus_snp2;
  static bool indiv_report;
  static string indiv_report_fid;
  static string indiv_report_iid;
  static bool plist;
  static string plist_fid1;
  static string plist_iid1;
  static string plist_fid2;
  static string plist_iid2;
  static bool merge_data;
  static bool merge_force_strand;
  static int  merge_mode;
  static bool merge_binary;
  static bool merge_list;
  static string merge_list_filename;
  static string merge_pedfile;
  static string merge_mapfile;
  static string merge_bedfile;
  static string merge_bimfile;
  static string merge_famfile;
  static bool write_snplist;
  static bool update_map;
  static bool update_cm;
  static bool update_chr;
  static bool update_name;
  static bool update_ids;
  static string update_ids_file;
  static bool update_sex;
  static string update_sex_file;
  static bool update_parents;
  static string update_parents_file;
  static bool update_pheno;
  static string update_pheno_file;

  static string update_mapfile;
  static string range_delimiter;

  static bool update_alleles;
  static string update_allele_file;

  static bool compound_genotype_code;

  static string tpedfile;
  static string tfamfile;
  static bool tfile_input;

  static string lpedfile;
  static string lfamfile;  
  static bool lfile_input;
  
  static bool ref_file;
  static string ref_file_name;
  
  static bool gvar;
  static bool gvar_write;
  static bool gvar_to_standard;
  static bool load_gvar;
  static bool gvar_verbose_association;
  static string gmapfile;
  static string gfamfile;
  static string gvarfile;
  static bool gvar_include_all_variants;
  static bool gvar_full_report;

  static bool flip_strand;
  static string flip_file;
  static bool flip_subset;
  static string flip_subset_file;
  


  static bool read_bitfile;
  static bool write_bitfile;
  static bool fast_binary;
  static string bitfilename;
  static string famfile;
  static string bitfilename_map;

  static bool SNP_major;
  static bool out_SNP_major;

  static bool compress_file;
  static bool uncompress_file;
  static string compress_filename;

  static bool read_ped;
  static string pedfile;
  static string mapfile;
  static bool ped_from_stdin;
  static string fileroot;
  static bool map3;
  static bool liability;

  static bool ped_skip_sex;
  static bool ped_skip_parents;
  static bool ped_skip_fid;
  static bool ped_skip_pheno;

  static string output_file_name;
  static bool silent;
  static bool gplink;
  static bool cli;

  static string missing_genotype;
  static string out_missing_genotype;
  static string missing_phenotype;
  static string out_missing_phenotype;
  static bool missing_genotype_explicit;
  static bool missing_phenotype_explicit;
  
  static bool ignore_missing_sex;

  static bool pheno_file;
  static bool covar_file;
  static bool clist;
  static bool no_show_covar;
  static bool dump_covar;
  static bool dump_covar_with_phenotype;
  static bool dump_covar_dummy_coding;
  static bool filter_on_covar;
  static int clist_number;
  static int plist_number;

  static bool snp_attrib_filter;
  static string snp_attrib_value;
  static string snp_attrib_file;
  static bool ind_attrib_filter;
  static string ind_attrib_value;  
  static string ind_attrib_file;

  static bool multiple_phenotypes;
  static string multiple_phenotype_file;

  static string make_pheno_filename;
  static string make_pheno_value;
  static bool make_pheno;
  static bool make_pheno_present;

  static bool dump_clst;

  static bool clist_selection; 
  static bool clist_selection_name; 
  static bool clist_selection_number;
  static string clist_selection_string;

  static bool plist_selection; 
  static bool plist_selection_name; 
  static bool plist_selection_number;
  static string plist_selection_string;

  static int mult_pheno;
  static string name_pheno;
  static bool all_pheno;
  static int mult_covar;
  static int mult_clst;
  static int mult_filter;
  static string filter_value;

  static string number_list_string;
  static bool   number_list_positive;

  static string pheno_filename;
  static string covar_filename;
  static string clist_filename;
  static string filter_filename;

  static bool cm_map;
  static double grid;
  static double fringe;
  static bool singlepoint;
  static int inter_grid;

  static bool done_global_pihat;

  static bool sol_family;

  static bool summ_nonfounders;
  static bool make_founders;
  static bool has_nonfounders;
  static bool make_missing_parents;

  static bool score_risk;
  static string score_risk_file;    
  static bool score_risk_ranges;
  static string score_risk_ranges_file;
  static int score_risk_ranges_min;
  static bool score_impute_expected;
  static bool score_risk_on_qrange;
  static string score_qrange_file;
  static string score_qfile;
  static bool score_test;
  static bool profile_sets;

  static bool report_missing;
  static bool test_missing;
  static bool mishap_test;
  static int  mishap_window;

  static bool calcFst;

  static bool proxy_assoc;
  static bool proxy_glm;
  static bool proxy_all;
  static bool proxy_full_report;
  static bool proxy_error;
  static bool proxy_impute;
  static bool proxy_impute_replace;
  static bool proxy_impute_preserve_genotyped;
  static bool proxy_record_dosage;
  static bool proxy_impute_genotypic_concordance;
  static double proxy_impute_threshold;
  static double proxy_info_threshold;
  static bool impute_verbose;
  static bool proxy_exclude;
  static string proxy_exclude_list;
  static bool proxy_exclude_from_file;
  static bool proxy_reference_only;

  static bool proxy_leave_out;
  static bool proxy_include_reference;
  static bool proxy_CC;
  static bool proxy_TDT;
  static string proxy_assoc_snp;
  static int proxy_window;
  static bool proxy_list;
  static string proxy_list_file;
  static bool proxy_all_list;
  static string proxy_all_list_file;
  static double proxy_kb;
  static double proxy_r2;
  static double proxy_maf;
  static double proxy_mhf;
  static double proxy_geno;
  static bool proxy_list_proxies;
  static int proxy_maxhap;
  static bool proxy_r2_filter;
  static double proxy_r2_filter_A;
  static double proxy_r2_filter_B;
  static double proxy_r2_filter_C;
  static int proxy_snp_filter;

  static double proxy_kb_planA;
  static int proxy_window_planA;
  static int proxy_snp_filter_planA;
  static double proxy_r2_filter_A_planA;
  static double proxy_r2_filter_B_planA;
  static double proxy_r2_filter_C_planA;

  static double proxy_planB_threshold;
  static double proxy_kb_planB;
  static int proxy_window_planB;
  static int proxy_snp_filter_planB;
  static double proxy_r2_filter_A_planB;
  static double proxy_r2_filter_B_planB;
  static double proxy_r2_filter_C_planB;

  static bool   greport;
  static string greport_results;
  static string greport_gene_list;
  static bool   greport_subset;
  static string greport_subset_file;
  static bool   greport_display_empty;

  static bool   annot_file;
  static string annot_filename;
  
  static bool   meta_analysis;
  static vector<string> meta_files;

  static bool   set_screen;
  static string set_screen_resultfile;

  static bool gettag_mode;
  static bool gettag_mode1;
  static bool gettag_mode2;
  static string gettag_file;
  static double gettag_r2;
  static int gettag_kb;
  static bool gettag_listall;

  static bool clumpld;
  static bool clumpld_best;
  static string clumpld_results;
  static string clumpld_column;
  static bool    clumpld_verbose;
  static bool    clumpld_indep;
  static int     clumpld_kb;
  static double  clumpld_r2;
  static double  clumpld_p1;
  static double  clumpld_p2;
  static bool    clumpld_index1;
  static bool    clumpld_only_show_replications;
  static bool    clumpld_only_show_replications_list;
  static bool	 clumpld_annot;
  static string	 clumpld_annot_fields;
  static string  clumpld_range_file;
  static bool    clumpld_range_annotate;
  static int     clumpld_min;

  static double min_af;
  static double max_af;
  static bool make_minor_allele;

  static double min_hf;
  static double max_hf;

  static int min_geno_cell;

  static double rarer_maf_threshold;
  static double rarer_dist_threshold;
  static int rarer_interval;
  static bool rare_test;
  static bool rare_test_weight1;
  static bool rare_test_print_details;
  static string rare_test_print_details_snp;
  static bool elf_pcmode;
  static bool elf_pcmode_2sided;
  static bool elf_baseline;

  static bool rare_test_score_range;
  static double rare_test_score_range_threshold;
  static string rare_test_score_results_file;
  static string rare_test_score_range_file;
  static bool rare_test_summary_controls;

  static vector<bool> chr_haploid;
  static vector<bool> chr_sex;
  static vector<bool> chr_Y;
  static vector<string> chr_code;
  static map<string,int> chr_map;

  static bool species_dog;
  static bool species_cow;
  static bool species_sheep;
  static bool species_horse;
  static bool species_rice;
  static bool species_mouse;

  static int run_start;
  static int run_end;
  static int run_chr;
  static string m1;
  static string m2;
  static double window;
  static bool position_window;
  static int from_window;
  static int to_window;

  static bool qt;
  static bool bt;
  static bool coding01;

  static bool ignore_phenotypes;
  static bool filter_cases;
  static bool filter_controls;
  static bool filter_males;
  static bool filter_females;
  static bool filter_founders;
  static bool filter_nonfounders;

  static bool SD;
  static bool CP;
  static bool affpair;  
  static bool remove_unaffected_pairs;
  static bool fix_prev;
  static double fixed_prev;

  static string tagfile;
  static string mapfile_impute;
  static bool make_tags;
  static bool impute_tags;
  static bool sliding_window;
  static string sliding_window_size;

  static bool make_blocks;

  static bool meta_large_phase;

  static bool phase_snps;
  static bool phase_hap_all;
  static double hap_post_prob;
  static double hap_missing_geno;
  static double hap_min_phase_prob;
  static int hap_max_nf_phases;
  static bool display_hap_freqs;

  static int haplo_plem_window;
  static int haplo_plem_overlap;
  static int haplo_plem_original_overlap;
  static int haplo_plem_iter;
  static bool haplo_plem_verbose;

  static bool haplo_plem_follow;
  static int haplo_plem_follow_ind;
  static string haplo_plem_follow_fid;
  static string haplo_plem_follow_iid;
  
  static int haplo_plem_likelihood_iter;
  static double haplo_plem_window_prune_phase;
 
  static double haplo_plem_window_tol;
  static double haplo_plem_zero_threshold;
  static bool haplo_plem_nonzero_threshold;
  
  static int haplo_plem_meta_window;
  static double haplo_plem_meta_prune_haplotype;
  static double haplo_plem_meta_prune_phase;
  static int haplo_plem_meta_iter;
  static int haplo_plem_meta_likelihood_iter;
  static double haplo_plem_meta_tol;
 
  static bool test_hap_CC;
  static bool test_hap_TDT;
  static bool test_hap_QTL;
  static bool test_hap_only;
  static bool test_hap_GLM;
  static bool test_hap_GLM_omnibus;
  static bool display_phase_probs;
  static bool display_phase_probs_wide;
  static bool weighted_mm;
  
  static bool chap_test;
  static bool chap_sole_variant;
  static bool chap_sole_variant_specific_alleles;
  static string chap_sole_variant_specific_allele_list;
  static bool chap_independent_effect;
  static bool chap_haplotype_specific;
  static string chap_entity;
  static bool chap_specified_groups;
  static bool chap_specified_snps;
  static string chap_model1;
  static string chap_model0;
  static bool chap_drop_snps;
  static string chap_drop_snps_list;
  static bool chap_add_grp_specifics;

  static bool assoc_test;
  static bool assoc_counts;
  static bool assoc_glm;
  static bool standard_beta;
  static bool assoc_glm_without_main_snp;
  static bool assoc_test_alt_perm;
  static bool full_model_assoc;
  static bool trend_only;
  static bool fisher_test;
  static bool return_beta;
  static bool hap_specific_snps;
  static string hap_specific_snps_list;

  static bool output_pheno_perm;

  static bool qt_means;
  static bool conditioning_snp_single;
  static string conditioning_snp_name;
  static bool conditioning_snps;
  static string conditioning_snps_file;

  static int xchr_model;
  static bool glm_sex_effect;
  static bool glm_no_auto_sex_effect;
  static bool glm_dominant;
  static bool glm_recessive;

  static double vif_threshold;

  static bool twoDFmodel;
  static bool twoDFmodel_hethom;
  static bool test_full_model;
  static bool simple_interaction;
  static vector<int> parameter_list;
  static vector<int> test_list;
  static bool glm_user_test;
  static bool glm_user_parameters;

  static bool qt_with_covariates;

  static bool model_perm_best;
  static bool model_perm_gen;
  static bool model_perm_dom;
  static bool model_perm_rec;
  static bool model_perm_trend;

  static bool assoc_gxe;

  static bool QTDT_test;
  static bool QFAM_total;
  static bool QFAM_between;
  static bool QFAM_within1;
  static bool QFAM_within2;
  static bool QFAM_adaptive;

  static bool TDT_test;
  static bool sibTDT_test;
  static bool mating_tests;
  static bool dfam_tdt;
  static bool dfam_sibs;
  static bool dfam_unrelateds;

  static bool perm_TDT_basic;
  static bool perm_TDT_parent;
  static bool discordant_parents;
  static bool parent_of_origin;
  static bool perm_POO_poo;
  static bool perm_POO_pat;
  static bool perm_POO_mat;
  static bool perm_POO_best;
  static bool built_families;

  static bool MENDEL_test;
  static bool MENDEL_report;
  static double MENDEL_snp;
  static double MENDEL_ind;

  static bool HWD_test;
  static bool HWD_report;
  static double HWD_limit;
  static bool HWD_standard;
  static bool HWD_filter_on_all;

  static bool CMH_test_1;
  static bool CMH_test_2;
  static bool CMH_test_ORD;
  static bool breslowday;

  static bool OR_homog_test;

  static double ci_level;
  static double ci_zt;
  static bool display_ci;

  static bool pfilter;
  static double pfvalue;

  static bool multtest;
  static bool use_GC;
  static bool fix_lambda;
  static double lambda;
  static bool qq_plot;
  static bool logscale;

  static bool ibs_sharing_test;

  static bool extract_set;
  static bool exclude_set;
  static bool snp_range_list;

  static bool thin_snps;
  static double thin_param;

  static bool make_set;
  static string make_set_file;
  static int make_set_border;
  static bool make_set_collapse;
  static bool make_set_ignore_group;
  static string make_set_collapse_label;
  static bool make_set_complement;
  static bool write_set;
  static bool read_set;

  static string exclude_file;
  static string extract_file;
  static string keep_file;
  static string remove_file;

  static bool read_snp_qual;
  static string snp_qual_file;
  static double snp_qual_min;
  static double snp_qual_max;
  static bool read_geno_qual;
  static string geno_qual_file;
  static double geno_qual_min;
  static double geno_qual_max;

  static bool snp_include_from_cl;
  static string snp_include_range;

  static bool dump_gene;
  static string dump_genename;
  
  static bool hotel;

  static bool set_test;
  static bool set_p2;
  static int set_min;
  static int set_max;
  static bool set_r2;
  static double set_r2_val;
  static bool set_r2_phase;
  static double set_chisq_threshold;
  static bool set_r2_write;
  static bool set_r2_read;
  static string set_r2_read_file;

  static string subsetfile;
  static bool use_subset;

  static string setfile;
  static bool set_score;
  static double set_score_p;
  static double set_step_in;
  static bool set_step;
  static bool set_table;

  static bool permute_within_sol;
  static bool boot;
  static bool disp_r1;
  static bool disp_r2;
  static bool disp_r_window;
  static int disp_r_window_snp;
  static int disp_r_window_kb;
  static double disp_r_window_r2;
  static bool ld_anchor;
  static bool ld_anchor_list;

  static bool flip_scan;
  static double flip_scan_threshold;
  static bool flip_scan_verbose;

  static bool prune_ld;
  static bool prune_ld_pairwise;
  static bool prune_ld_pairwise_maf;
  static double prune_ld_vif;
  static double prune_ld_r2;
  static int prune_ld_win;
  static int prune_ld_step;
  static bool prune_r2_prefer;
  static string prune_r2_prefer_list;
  static bool prune_r2_fixed;
  static string prune_r2_fixed_list;

  static bool calc_SNPSNP_LD;
  static string ld_SNP1;
  static string ld_SNP1_file;
  static string ld_SNP2;

  static bool epistasis;
  static bool fast_epistasis;

  static bool epi_caseonly;
  static double epi_caseonly_kb_gap;
  static bool epi_filter;
  static double epi_alpha1;  
  static double epi_alpha2;
  static bool set_by_set;
  static bool epi_genebased;
  static bool epi_quickscan;

  static bool drop_sets;

  static bool inbreeding;
  static bool check_sex;
  static bool impute_sex;
  static double sex_threshold_male;
  static double sex_threshold_female;

  static bool homo_run;
  static bool homo_run_consensus_match;
  static bool homo_run_kb;
  static bool homo_run_snps;
  static double homo_run_density;
  static int homo_run_gap;
  static bool homo_miss_as_hom;

  static int homo_windowSize;
  static int homo_windowKB; 
  static int homo_windowAllowedHet;
  static int homo_windowAllowedMissing;
  static double homo_threshold;

  static int homo_run_length_kb;
  static int homo_run_length_snps;
  static int homo_run_het;
  static bool homo_summary_allelic_match;
  static double fuzzy_homo;
  static bool homozyg_verbose;
  static int pool_size_min;

  static bool ibs_run;
  static int ibs_run_length_snps;
  static int ibs_run_length_kb;
  static double ibs_run_density;
  static int ibs_inner_run_length_kb;
  static int ibs_inner_run_length_snp;
  static int ibs_join_kb;
  static int ibs_join_snp;
  static int ibs_run_missing;
  static int ibs_run_0;
  static int ibs_inter_snp_distance;
  static bool ibs_2only;

  static bool miss_run;
  static int miss_run_length;
  static bool miss_run_length_kb;
  static double miss_run_level;

  static bool segment_haplotrack;
  static string segment_haplotrack_fid1;
  static string segment_haplotrack_iid1;
  static string segment_haplotrack_fid2;
  static string segment_haplotrack_iid2;

  static bool mk_datfile;
  static bool segment_output;
  static bool segment_minimal;
  static bool segment_silently_return_groups;
  static int segment_current_focal_snp;
  static bool segment_overlap;
  static bool segment_verbose;
  static bool segment_validate;
  static bool segment_test_individual;
  static bool segment_test_specific_segs;
  static bool segment_test_fisher;
  static bool segment_test_1sided;
  static bool segment_test_force_1sided;
  static bool segment_test_ignore_discordant;
  static int segment_snp1;
  static int segment_snp2;
  static string segment_m1;
  static string segment_m2;
  static bool force_span;
  static int segment_length;
  static int segment_snp;
  static bool segment_output_started;
  static bool read_segment_file;
  static string read_segment_filename;

  static int segment_inter_snp_distance;
  static bool multi_output;
  static bool gmulti_output;
  static bool pihat_filter;
  static bool genome_output;
  static bool compress_genome;
  static bool genome_only_check_rels;
  static bool genome_output_minimal;
  static bool genome_output_full;
  static bool genome_2sets;
  static string genome_setlist1;
  static string genome_setlist2;
  static bool genome_test;
  static double genome_test_threshold;
  static int genome_test_min_snp;
  static bool ibs_test;
  static int ibs_test_min_snp;
  static bool ibs_test_method2; 
  static bool summary_ibd_output;
  static double IBD_threshold;
  static double segment_threshold_start;
  static double segment_threshold_finish;
  static bool nudge;
  static bool bound;
  static bool show_impossible_IBD;
  static bool IBD_within;
  
  
  static bool permute;
  static int replicates;
  static bool perm_count;
  static bool mperm_save_best;
  static bool mperm_save_all;
  static bool mperm_rank;
  static bool adaptive_perm;
  static int adaptive_min;
  static int adaptive_max;
  static int adaptive_interval;
  static double adaptive_interval2;
  static double adaptive_alpha;
  static double adaptive_ci;

  static bool perm_genedrop;
  static bool perm_genedrop_and_swap;
  static bool perm_genedrop_unrel;
  static bool perm_genedrop_parents;
  static bool perm_genedrop_sibships;


  static bool FIXED;
  static bool FIXED_p;
  static Z FIX_IBD;
  static double FIX_p;

  static bool matrix;
  static bool distance_matrix;
  static bool cluster;
  static bool cluster_euclidean;
  static bool cluster_group_avg;
  static bool cluster_plot;
  static bool force_initial_cluster;
  static int cluster_mds_dim;
  static bool mds_by_individual;
  static bool genome_groups;
  static bool cluster_ibm_constraint;
  static double cluster_ibm_constraint_value;
  static bool cluster_missing;
  static bool cluster_selcon;
  static string cluster_selcon_file;
  static int max_cluster_N;
  static double merge_p;
  static int ibstest_gap;
  static int max_cluster_size;
  static int max_cluster_case;
  static int max_cluster_control;
  static bool include_cluster;
  static bool include_cluster_from_file;
  static string include_cluster_filename;
  static int analyse_cluster;
  static bool cluster_on_phenotype;
  static bool cluster_on_mcc;
  static int min_neighbour;
  static int max_neighbour;
  static bool outlier_detection;
  static bool bmatch;
  static bool bmatch_usertype;
  static bool qmatch;
  static string bmatch_filename;
  static string bmatch_direction_filename;  
  static string qmatch_filename;
  static string qmatch_threshold_filename;  
  
  static bool include_all_pairs;
  static double include_all_z1;    

  static double MIN_PIHAT;
  static double MAX_PIHAT;
  static double MAX_CORR_PIHAT_PIHAT_G;
  static double MAX_GENO_MISSING;
  static double MAX_IND_MISSING;
  static int    MAX_LINE_LENGTH;

  static bool remove_indiv;
  static string remove_indiv_list;
  static string keep_indiv_list;
  static bool keep_indiv;
  static bool extract_before_exclude;
  static bool remove_before_keep;

  static bool locked;
  
  static bool af_read;
  static bool af_write;

  static bool ibd_read;
  static string ibd_file;
  static bool ibd_read_minimal;
  static bool ibd_read_list;
  static string ibd_file_list;

  static string af_file;
  static bool af_count;
  
  static bool inc_write;
  static bool inc_read;
  static string inc_file;

  static int pp_maxsnp;
  static int pp_maxfid;
  static int pp_maxiid;
  
  static int BATCH_SIZE;
  
  static bool plink;

  static bool display_segment_long;
  static bool display_cnv_track;
  static int cnv_col;
  static bool cnv_makemap;
  static bool cnv_writelist;
  static bool cnv_list;
  static string cnv_listname;
  static int cnv_min_kb;
  static double cnv_min_score;
  static int cnv_min_sites;
  static int cnv_max_kb;
  static double cnv_max_score;
  static int cnv_max_sites;
  static bool cnv_del_only;
  static bool cnv_dup_only;
  static int cnv_type;
  static bool cnv_intersect;
  static bool cnv_exclude;
  static string cnv_intersect_file;
  static bool cnv_intersect_subset;
  static string cnv_intersect_subset_file;
  static bool cnv_count;
  static double cnv_overlap;
  static bool cnv_defined_overlap;
  static bool cnv_indiv_perm;
  static bool cnv_pos_perm;
  static bool cnv_drop_no_segment;
  static bool cnv_freq_method2;
  static double cnv_freq_method2_threshold;
  static bool cnv_write_freq;
  static bool cnv_freq_include;
  static bool cnv_freq_include_below;
  static bool cnv_freq_include_exact;
  static bool cnv_freq_include_exact_exclude;
  static int cnv_freq_include_cnt;
  static bool cnv_unique;
  static bool cnv_intersect_writeback;
  static bool cnv_intersect_writeback_verbose;
  static bool cnv_disrupt;
  static int cnv_region_border;
  static bool cnv_union_overlap;
  static bool cnv_region_overlap;
  static bool cnv_check_overlap;
  static bool cnv_count_baseline;
  static string cnv_count_baseline_file;
  static bool cnv_weighted_gene_test;
  static bool cnv_enrichment_test;
  static int cnv_en_model;
  static bool cnv_glm;

  static bool seg_test_window;
  static double seg_test_window_bp;
  static bool seg_test_region;

  static bool dosage_assoc;
  static string dosage_file;

  static bool dosage_hard_call;
  static double dosage_hard_call_thresh;
  static int dosage_hard_call_thresh2;
  static bool dosage_hasMap;
  static bool write_dosage;
  
};


void setOptions(CArgs &);
void getOutputFilename(CArgs &);



#endif
