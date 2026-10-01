

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
#include "options.h"

using namespace std;

// Temporary dummy function
bool par::myfunction = false;

Options par::opt;

bool   par::verbose = false;
bool   par::flag = false;
bool   par::dumpped = false;
bool   par::debug = false;
bool   par::dummy = false;
int    par::dummy_nind = 0;
int    par::dummy_nsnp = 0;
bool   par::web_check = true;
bool   par::tucc = false;
bool   par::do_not_load_snps = false;

double const par::epsilon = 1e-12;

long unsigned int par::random_seed = 0;

int par::simul_ncases = 1000;
int par::simul_ncontrols = 1000;
string par::simul_label = "";
double par::simul_prevalence = 0.01;
bool par::simul = false;
string par::simul_file = "";
bool par::simul_tags = false;
bool par::simul_haps = false;
bool par::simul_qt = false;
double par::simul_qt_var = 0.05;

bool   par::lookup = false;
bool   par::lookup_single_snp = false;
bool   par::lookup_to_file = false;
bool   par::lookup_gene = false;
bool   par::lookup_multiple_genes = false;
string par::lookup_gene_name = "GENE1";
string par::lookup_snp = "rs1234";
int    par::lookup_gene_kb_window = 20;
int    par::lookup_snp_kb_window = 100;

bool    par::lookup2 = false;
string  par::lookup2_cmd = "";

bool par::idhelp = false;
string par::idhelp_output_delimit = " ";
bool   par::idhelp_dump_from_dict = false;
string par::idhelp_dump_from_dict_cmd = "";
string par::idhelp_dictionary = "";
bool   par::idhelp_auto_alias = false;
bool   par::idhelp_lookup = false;
string par::idhelp_lookup_string = "";
bool   par::idhelp_subset = false;
string par::idhelp_subset_string ="";
bool   par::idhelp_replace = false;
string par::idhelp_replace_string = "";

bool   par::idhelp_match = false;
vector<string> par::idhelp_match_string;

bool   par::idhelp_no_dict = false;
bool   par::idhelp_list_aliases = false;
bool   par::idhelp_alias_update = true;

bool   par::run_R_script = false;
bool   par::run_R_write_script = false;
string par::R_script = "script.R";
bool   par::run_R_chisq = false;
bool   par::run_R_z = false;
int    par::run_R_nsnps = 100;
int    par::R_port = 6311;

bool   par::recode = false;
bool   par::recode_transpose = false;
bool   par::recode_long = false;
bool   par::recode_long_ref = false;
bool   par::recode_mutlist = false;
bool   par::recode_12 = false;
bool   par::recode_AD = false;
bool   par::recode_AD_fixed = false;
bool   par::recode_AD_Aonly = false;
bool   par::recode_allele_coding = false;
string par::recode_allele_coding_file = "file.lst";
string par::recode_delimit = " ";
string par::recode_indelimit = " ";
bool   par::recode_HV = false;
bool   par::recode_whap = false;
bool   par::recode_fastphase = false;
bool   par::recode_structure = false;
bool   par::recode_bimbam = false;
bool   par::recode_1234 = false;
bool   par::recode_ACGT = false;

bool   par::set_reference_allele = false;
string par::set_reference_allele_file = "dummy.file";
bool   par::lfile_allele_count = false;

bool   par::preserve_all_genotypes = false;
bool   par::preserve_mendel_errors = false;
bool   par::zero_cluster = false;
string par::zero_cluster_filename = "plink.zero";
bool   par::oblig_missing = false;
string par::oblig_missing_filename = "plink.zero";
string par::oblig_clusters_filename = "plink.clst";
bool   par::loop_over = false;
string par::loop_over_label = "";
int    par::loop_counter = 0;
string par::loop_over_filename = "plink.clst";
bool   par::list_by_allele = false;
bool   par::list_twolocus = false;
string par::twolocus_snp1 = "";
string par::twolocus_snp2 = "";
bool   par::indiv_report = false;
string par::indiv_report_fid = "fid1";
string par::indiv_report_iid = "iid1";
bool   par::plist = false;
string par::plist_fid1 = "";
string par::plist_iid1 = "";
string par::plist_fid2 = "";
string par::plist_iid2 = "";
bool   par::merge_data = false;
bool   par::merge_force_strand = false;
int    par::merge_mode = 1;
bool   par::merge_binary = false;
bool   par::merge_list = false;
string par::merge_list_filename = "merge.list";
string par::merge_pedfile = "merge.ped";
string par::merge_mapfile = "merge.map";
string par::merge_bedfile = "merge.bed";
string par::merge_bimfile = "merge.bim";
string par::merge_famfile = "merge.fam";
bool par::write_snplist = false;
bool par::update_map = false;
bool par::update_cm = false;
bool par::update_chr = false;
bool par::update_name = false;
bool par::update_ids = false;
string par::update_ids_file = "";
bool par::update_sex = false;
string par::update_sex_file = "";
bool par::update_parents = false;
string par::update_parents_file = "";
bool par::update_pheno = false;
string par::update_pheno_file = "";

string par::update_mapfile = "new.map";
string par::range_delimiter = "-";

bool par::update_alleles = false;
string par::update_allele_file = "dummy";

bool par::compound_genotype_code = false;

string par::tpedfile = "plink.tped";
string par::tfamfile = "plink.tfam";
bool par::tfile_input = false;

string par::lpedfile = "plink.lgen";
bool par::lfile_input = false;

bool par::ref_file = false;
string par::ref_file_name = "";

bool par::gvar = false;
bool par::gvar_write = false;
bool par::gvar_to_standard = false;
bool par::load_gvar = false;
bool par::gvar_include_all_variants = false;
bool par::gvar_verbose_association = false;
string par::gmapfile = "plink.map";
string par::gfamfile = "plink.fam";
string par::gvarfile = "plink.gvar";

bool par::gvar_full_report = false;

bool   par::flip_strand = false;
string par::flip_file = "plink.flip";
bool par::flip_subset = false;
string par::flip_subset_file = "plink.file";

bool par::compress_file = false;
bool par::uncompress_file = false;
string par::compress_filename = "";

bool par::read_ped = false;
string par::pedfile  = "plink.ped";
string par::mapfile  = "plink.map";
string par::fileroot = "plink";
bool par::ped_from_stdin = false;
bool par::map3 = false;
bool par::liability = false;

bool par::ped_skip_sex = false;
bool par::ped_skip_parents = false;
bool par::ped_skip_fid = false;
bool par::ped_skip_pheno = false;

bool par::SNP_major = true;
bool par::out_SNP_major = true;

string par::output_file_name = "plink";
bool par::silent = false;
bool par::gplink = false;
bool par::cli = false;

bool par::fast_binary = false;
string par::bitfilename = "plink.bed";
string par::famfile = "plink.fam";
string par::bitfilename_map = "plink.bim";

bool par::write_bitfile = false;
bool par::read_bitfile = false;

bool par::pheno_file = false;
bool par::covar_file = false;
bool par::clist = false;
bool par::no_show_covar = false;
bool par::dump_covar = false;
bool par::dump_covar_with_phenotype = false;
bool par::dump_covar_dummy_coding = false;
bool par::filter_on_covar = false;
int par::clist_number = 0;
int par::plist_number = 0;

bool par::snp_attrib_filter = false;
string par::snp_attrib_value = "";
string par::snp_attrib_file = "";

bool par::ind_attrib_filter = false;
string par::ind_attrib_value = "";
string par::ind_attrib_file = "";

bool par::multiple_phenotypes = false;
string par::multiple_phenotype_file = ""; 

string par::make_pheno_filename  = "";
string par::make_pheno_value = "";
bool par::make_pheno = false;
bool par::make_pheno_present = false;

bool par::dump_clst = false;

bool par::clist_selection = false;
bool par::clist_selection_name = false;
bool par::clist_selection_number = false;
string par::clist_selection_string = "";

bool par::plist_selection = false;
bool par::plist_selection_name = false;
bool par::plist_selection_number = false;
string par::plist_selection_string = "";

int par::mult_pheno = 1;
string par::name_pheno = "";
bool par::all_pheno = false;
int par::mult_covar = 1;
int par::mult_clst = 1;
int par::mult_filter = 1;
string par::filter_value = "1";

string par::number_list_string = "";
bool   par::number_list_positive = true;

string par::pheno_filename = "plink.phe";
string par::covar_filename = "plink.cov";
string par::clist_filename = "plink.cov";
string par::filter_filename = "plink.cov";

string par::missing_genotype = "0";
string par::missing_phenotype = "-9";
string par::out_missing_genotype = "0";
string par::out_missing_phenotype = "-9";
bool par::missing_phenotype_explicit = false;
bool par::missing_genotype_explicit = false;
bool par::ignore_missing_sex = false;

bool par::cm_map = false;
double par::grid = 0.005; // 0.5cM = 500kb grid
double par::fringe = .01; // 1 cM fringe 
bool par::singlepoint = false;
int par::inter_grid = 2;

bool par::done_global_pihat = false;

bool par::summ_nonfounders = false;
bool par::make_founders = false;
bool par::has_nonfounders = false;
bool par::make_missing_parents = false;

bool par::report_missing = false;
bool par::test_missing = false;
bool par::mishap_test = false;
int par::mishap_window = 1;

bool par::calcFst = false;

bool par::score_risk = false;
string par::score_risk_file = "plink.risk";
bool par::score_risk_ranges = false;
string par::score_risk_ranges_file = "plink.ranges";
int par::score_risk_ranges_min = 0;
bool par::score_impute_expected = true;
bool par::score_risk_on_qrange = false;
string par::score_qrange_file = "";
string par::score_qfile = "";
bool par::score_test = false; 
bool par::profile_sets = false;

bool par::proxy_assoc = false;
bool par::proxy_glm = false;
bool par::proxy_all = false;
bool par::proxy_full_report = false;
bool par::proxy_impute = false;
bool par::proxy_impute_replace = false;
bool par::proxy_impute_preserve_genotyped = false;
bool par::proxy_impute_genotypic_concordance = false;
bool par::proxy_record_dosage = false;
bool par::proxy_error = false;
bool par::proxy_leave_out = false;
bool par::proxy_include_reference = false;
bool par::proxy_CC = false;
bool par::proxy_TDT = false;
string par::proxy_assoc_snp = "rs1234";
bool par::proxy_list = false;
string par::proxy_list_file = "proxy.hap";
bool par::proxy_all_list = false;
string par::proxy_all_list_file = "proxy.list";
bool par::proxy_list_proxies = false;
bool par::proxy_exclude = false;
string par::proxy_exclude_list = "pexclude.list";
bool par::proxy_exclude_from_file = false;
bool par::proxy_reference_only = false;

int par::proxy_maxhap = 3;
double par::proxy_r2 = 0.5;
double par::proxy_info_threshold = 0.5;
bool   par::impute_verbose = false;

double par::proxy_maf = 0.005;
double par::proxy_mhf = 0.01;
double par::proxy_geno = 0.2;
double par::proxy_impute_threshold = 0.9;
bool   par::make_minor_allele = true;

double par::proxy_planB_threshold = 0.1;
double par::proxy_kb_planB = 500;
int par::proxy_window_planB = 30;
int par::proxy_snp_filter_planB = 10;
double par::proxy_r2_filter_A_planB = 0.00;
double par::proxy_r2_filter_B_planB = 0.01;
double par::proxy_r2_filter_C_planB = 0.50;

double par::proxy_kb_planA = 250;
int par::proxy_window_planA = 15;
int par::proxy_snp_filter_planA = 5;
double par::proxy_r2_filter_A_planA = 0.00;
double par::proxy_r2_filter_B_planA = 0.25;
double par::proxy_r2_filter_C_planA = 0.50;

double par::proxy_kb = 250;
int par::proxy_window = 15;
int par::proxy_snp_filter = 5;
bool par::proxy_r2_filter = true;
double par::proxy_r2_filter_A = 0.00;
double par::proxy_r2_filter_B = 0.05;
double par::proxy_r2_filter_C = 0.50;

bool   par::greport = false;
string par::greport_results = "file1";
string par::greport_gene_list = "file2";
bool   par::greport_subset = false;
string par::greport_subset_file = "file3";
bool   par::greport_display_empty = false;

bool   par::annot_file = false;
string par::annot_filename = "";

bool   par::meta_analysis = false;
vector<string> par::meta_files;

bool   par::set_screen = false;
string par::set_screen_resultfile = "";

bool   par::gettag_mode = false;
bool   par::gettag_mode1 = true;
bool   par::gettag_mode2 = false;
string par::gettag_file = "";
double par::gettag_r2 = 0.8;
int    par::gettag_kb = 250000; // 250kb default, in BP
bool   par::gettag_listall = false;

bool    par::clumpld = false;
bool    par::clumpld_best = false;
string  par::clumpld_results = "plink.assoc";
string  par::clumpld_column = "P";
bool    par::clumpld_verbose = false;
bool    par::clumpld_indep = true;
int     par::clumpld_kb = 250000;
double  par::clumpld_r2 = 0.5;
double  par::clumpld_p1 = 1e-4;
double  par::clumpld_p2 = 1e-2;
bool    par::clumpld_index1 = false;
bool    par::clumpld_only_show_replications = false;
bool    par::clumpld_only_show_replications_list = false;
bool	par::clumpld_annot = false;
string	par::clumpld_annot_fields = "";
string  par::clumpld_range_file = "range.list";
bool    par::clumpld_range_annotate = false;
int     par::clumpld_min = 0; // NOT USED

double par::min_af = 0.01;
double par::max_af = 1;  // max minor allele freq

double par::min_hf = 0.01;
double par::max_hf = 1;

int par::min_geno_cell = 5;

double par::rarer_maf_threshold = 0.1;
double par::rarer_dist_threshold = 100000; // 100 kb
int par::rarer_interval = 100; // in bp
bool par::rare_test = false;
bool par::rare_test_weight1 = false;
bool par::rare_test_print_details = false;
string par::rare_test_print_details_snp = "";
bool par::elf_pcmode = false;
bool par::elf_pcmode_2sided = false;
bool par::elf_baseline = false;

bool par::rare_test_score_range = false;
double par::rare_test_score_range_threshold = 0.01;
string par::rare_test_score_results_file = "";
string par::rare_test_score_range_file = "";
bool par::rare_test_summary_controls = false;

vector<bool> par::chr_haploid(0);
vector<bool> par::chr_sex(0);
vector<bool> par::chr_Y(0);
vector<string> par::chr_code(0);
map<string,int> par::chr_map;

bool par::species_dog = false;
bool par::species_cow = false;
bool par::species_horse = false;
bool par::species_sheep = false;
bool par::species_rice = false;
bool par::species_mouse = false;

int par::run_chr = 0;
int par::run_start = 0;
int par::run_end = 0;
string par::m1 = "";
string par::m2 = "";
double par::window = 0; // kb 
bool par::position_window = false;
int par::from_window = 0; // bp 
int par::to_window = 0;   // bp 

bool par::mk_datfile = false;
bool par::qt = false;
bool par::bt = true;
bool par::coding01 = false;
bool par::ignore_phenotypes = true;
bool par::filter_cases = false;
bool par::filter_controls = false;
bool par::filter_males = false;
bool par::filter_females = false;
bool par::filter_founders = false;
bool par::filter_nonfounders = false;

bool par::segment_haplotrack = false;
string par::segment_haplotrack_fid1 = "1";
string par::segment_haplotrack_iid1 = "1";
string par::segment_haplotrack_fid2 = "2";
string par::segment_haplotrack_iid2 = "2";

bool par::segment_output = false;
bool par::segment_minimal = false;
bool par::segment_silently_return_groups = false;
int  par::segment_current_focal_snp = -1;
bool par::segment_overlap = false;
bool par::segment_verbose = false;
bool par::segment_validate = false;
bool par::segment_test_individual = false;
bool par::segment_test_specific_segs = false;
bool par::segment_test_fisher = false;
bool par::segment_test_1sided = true;
bool par::segment_test_force_1sided = false;
bool par::segment_test_ignore_discordant = false;
int par::segment_snp1 = -1;
int par::segment_snp2 = -1;
string par::segment_m1 = "";
string par::segment_m2 = "";
bool par::force_span = false;
int  par::segment_length = 1000000; // 1000kb default min length
int  par::segment_snp = 100; // 100 SNPs 
bool par::segment_output_started = false;
bool par::read_segment_file = false;
string par::read_segment_filename = "";
int par::segment_inter_snp_distance = 1000; //unit = kb
bool par::multi_output = false;
bool par::gmulti_output = false;
bool par::pihat_filter = true;
bool par::genome_output = false;
bool par::compress_genome = false;
bool par::genome_only_check_rels = false;
bool par::genome_output_minimal = false;
bool par::genome_output_full = false;
bool par::genome_2sets = false;
string par::genome_setlist1 = "plink.set1";
string par::genome_setlist2 = "plink.set2";
bool par::genome_test = false;
double par::genome_test_threshold = 0.01;
int par::genome_test_min_snp = 20;
bool par::ibs_test = false;
int par::ibs_test_min_snp = 20;
bool par::ibs_test_method2 = false;
bool par::summary_ibd_output = false;
double par::IBD_threshold = 0.2;
double par::segment_threshold_start = 0.25;
double par::segment_threshold_finish = 0.25;
bool par::nudge = false;
bool par::bound = true;
bool par::show_impossible_IBD = true;
bool par::IBD_within = false;


bool par::SD = true;
bool par::CP = false;
bool par::affpair = false;
bool par::remove_unaffected_pairs = false;
bool par::fix_prev = false;
double par::fixed_prev = 0;

bool par::sol_family = false;

string par::tagfile = "plink.tag";
string par::mapfile_impute = "plink.impute.map";
bool par::impute_tags = false;
bool par::sliding_window = false;
string par::sliding_window_size = "2";

bool par::make_blocks = false;

bool par::meta_large_phase = false;
bool par::phase_snps = false;
bool par::phase_hap_all = false;
double par::hap_post_prob = 0.8;
double par::hap_missing_geno = 0.5;
int par::hap_max_nf_phases = 1024;
double par::hap_min_phase_prob = 1e-2;
bool par::display_hap_freqs = false;

bool par::haplo_plem_verbose = false;
bool par::haplo_plem_follow = false;
int par::haplo_plem_follow_ind = -1;
string par::haplo_plem_follow_fid = "FID1";
string par::haplo_plem_follow_iid = "IID1";

int par::haplo_plem_window = 6;
int par::haplo_plem_overlap = 2;
int par::haplo_plem_original_overlap = 2;
int par::haplo_plem_iter = 20;
int par::haplo_plem_likelihood_iter = 5;
double par::haplo_plem_window_prune_phase = 1e-10;
double par::haplo_plem_window_tol = 1e-4;
double par::haplo_plem_zero_threshold = -1;
bool par::haplo_plem_nonzero_threshold = true;

int par::haplo_plem_meta_window = 2;
double par::haplo_plem_meta_prune_haplotype = 1e-6;
double par::haplo_plem_meta_prune_phase = 0.01;
int par::haplo_plem_meta_iter = 200;
int par::haplo_plem_meta_likelihood_iter = 5;
double par::haplo_plem_meta_tol = 1e-4;


bool par::test_hap_CC = false;
bool par::test_hap_TDT = false;
bool par::test_hap_QTL = false;
bool par::test_hap_GLM = false;
bool par::test_hap_GLM_omnibus = false;
bool par::test_hap_only = false;
bool par::display_phase_probs = false;
bool par::display_phase_probs_wide = false;
bool par::weighted_mm = false;

bool par::chap_test = false;
bool par::chap_sole_variant = false;
bool par::chap_independent_effect = false;
bool par::chap_sole_variant_specific_alleles = false;
string par::chap_sole_variant_specific_allele_list = "";
bool par::chap_haplotype_specific = false;
string par::chap_entity = "";
bool par::chap_specified_groups = false;
bool par::chap_specified_snps = false;
string par::chap_model1 = "";
string par::chap_model0 = "";
bool par::chap_drop_snps = false;
string par::chap_drop_snps_list = "";
bool par::chap_add_grp_specifics = false;

bool par::assoc_test = false;
bool par::assoc_counts = false;
bool par::assoc_glm = false;
bool par::standard_beta = false;
bool par::assoc_glm_without_main_snp = false;
bool par::assoc_test_alt_perm = false;
bool par::full_model_assoc = false;
bool par::trend_only = false;
bool par::fisher_test = false;
bool par::return_beta = false;

bool par::hap_specific_snps = false;
string par::hap_specific_snps_list = "";

bool par::qt_means = false;

bool par::conditioning_snp_single = false;
string par::conditioning_snp_name = "rs1234";
bool par::conditioning_snps = false;
string par::conditioning_snps_file = "plink.list";

int par::xchr_model = 1;
bool par::glm_sex_effect = false;
bool par::glm_no_auto_sex_effect = false;
bool par::glm_dominant = false;
bool par::glm_recessive = false;

double par::vif_threshold = 50; 

bool par::twoDFmodel = false;
bool par::twoDFmodel_hethom = false;
bool par::test_full_model = false;
bool par::simple_interaction = false;
vector<int> par::parameter_list(0);
vector<int> par::test_list(0);
bool par::glm_user_test = false;
bool par::glm_user_parameters = false;

bool par::qt_with_covariates = false;

bool par::model_perm_best = false;
bool par::model_perm_gen = false;
bool par::model_perm_dom = false;
bool par::model_perm_rec = false;
bool par::model_perm_trend = false;

bool par::output_pheno_perm = false;

bool par::assoc_gxe = false;

bool par::QTDT_test = false;
bool par::QFAM_total = false;
bool par::QFAM_between = false;
bool par::QFAM_within1 = false;
bool par::QFAM_within2 = false;
bool par::QFAM_adaptive = false;

bool par::TDT_test = false;
bool par::sibTDT_test = false;
bool par::mating_tests = false;
bool par::dfam_tdt = true;
bool par::dfam_sibs = true;
bool par::dfam_unrelateds = true;

bool par::perm_TDT_basic = true;
bool par::perm_TDT_parent = false;
bool par::discordant_parents = false;
bool par::parent_of_origin = false;
bool par::perm_POO_poo = true;
bool par::perm_POO_pat = false;
bool par::perm_POO_mat = false;
bool par::perm_POO_best = false;
bool par::built_families = false;

bool par::HWD_test = false;
bool par::HWD_report = false;
double par::HWD_limit = 0.001;
bool par::HWD_standard = false;
bool par::HWD_filter_on_all = false;

bool par::MENDEL_test = false;
bool par::MENDEL_report = false;
double par::MENDEL_ind = 0.1;
double par::MENDEL_snp = 0.1;

bool par::CMH_test_1 = false;
bool par::CMH_test_2 = false;
bool par::CMH_test_ORD = false;
bool par::breslowday = false;

bool par::OR_homog_test = false;

double par::ci_level = 0.95;
double par::ci_zt = 0;
bool par::display_ci = false;

bool par::pfilter = false;
double par::pfvalue = 1e-5;

bool par::multtest = false;
bool par::use_GC = false;
bool par::fix_lambda = false;
double par::lambda = 1;

bool par::qq_plot = false;
bool par::logscale = false;

bool par::ibs_sharing_test = false;

string par::keep_file = "plink.list";
string par::remove_file = "plink.list";

bool par::extract_set = false;
bool par::exclude_set = false;
string par::exclude_file = "plink.list";
string par::extract_file = "plink.list";
bool par::snp_range_list = false;

bool par::thin_snps = false;
double par::thin_param = 0;

bool par::read_snp_qual = false;
string par::snp_qual_file = "dummy";
double par::snp_qual_min = 0;
double par::snp_qual_max = 1;

bool par::read_geno_qual = false;
string par::geno_qual_file = "dummy";
double par::geno_qual_min = 0;
double par::geno_qual_max = 1;


bool par::make_set = false;
string par::make_set_file = "plink.set";
int par::make_set_border = 0;
bool par::make_set_collapse = false;
bool par::make_set_ignore_group = false;
string par::make_set_collapse_label = "SET";
bool par::make_set_complement = false;
bool par::write_set = false;
bool par::read_set = false;

bool par::drop_sets = true;

bool par::snp_include_from_cl = false;
string par::snp_include_range = "";

bool par::dump_gene = false;
string par::dump_genename = "";

bool   par::permute           = false;
bool   par::perm_count        = false;
bool   par::mperm_save_all    = false;
bool   par::mperm_save_best   = false;
bool   par::mperm_rank        = false;
int    par::replicates        = 1000;
bool   par::adaptive_perm     = true;
int    par::adaptive_min      = 5;
int    par::adaptive_max      = 1000000;
int    par::adaptive_interval = 1;
double par::adaptive_interval2= 0.001;
double par::adaptive_alpha    = 0.00;
double par::adaptive_ci       = 0.0001;

bool   par::perm_genedrop          = false;
bool   par::perm_genedrop_and_swap = false;
bool   par::perm_genedrop_unrel    = false;
bool   par::perm_genedrop_parents  = false;
bool   par::perm_genedrop_sibships = false;

bool par::hotel = false;

bool   par::set_test              = false;
bool   par::set_p2                = false;
int    par::set_min               = -1;   
int    par::set_max               = 5;
double par::set_r2_val            = 0.5;
bool   par::set_r2                = false; 
bool   par::set_r2_phase          = false;
double par::set_chisq_threshold   = 3.84146;
bool   par::set_r2_write          = false;
bool   par::set_r2_read           = false;
string par::set_r2_read_file      = "plink.ldset";

string par::subsetfile            = "dummy.file";
bool   par::use_subset            = false;

string par::setfile               = "plink.set";
bool par::set_score               = false;
double par::set_score_p           = 1;
bool par::set_table               = false;

double par::set_step_in           = 0.05;
bool   par::set_step              = false;

bool par::permute_within_sol      = false;

bool par::boot                    = false;

bool par::disp_r1                 = false;
bool par::disp_r2                 = false;
bool par::disp_r_window           = false;
int par::disp_r_window_snp        = 10;
int par::disp_r_window_kb         = 1000000;
double par::disp_r_window_r2      = 0.2;
bool par::ld_anchor               = false;
bool par::ld_anchor_list          = false;

bool par::flip_scan               = false;
double par::flip_scan_threshold   = 0.5;
bool   par::flip_scan_verbose     = false;

bool par::prune_ld = false;
bool par::prune_ld_pairwise = false;
bool par::prune_ld_pairwise_maf = true;
double par::prune_ld_vif = 2;
double par::prune_ld_r2 = 1 - 1e-6;
int par::prune_ld_win = 100;
int par::prune_ld_step = 50;
bool par::prune_r2_prefer = false;
string par::prune_r2_prefer_list = "dummy";
bool par::prune_r2_fixed = false;
string par::prune_r2_fixed_list = "dummy";

bool par::calc_SNPSNP_LD = false;
string par::ld_SNP1 = "";
string par::ld_SNP1_file ="";
string par::ld_SNP2 = "";

double par::epi_alpha1 = 0.0001;
double par::epi_alpha2 = 0.01;
bool par::epi_filter = true;
bool par::set_by_set = true;
bool par::epistasis = false;
bool par::fast_epistasis = false;
bool par::epi_caseonly = false;
double par::epi_caseonly_kb_gap = 1000; // 1Mb default gap in SNPxSNP tests
bool par::epi_genebased = false;
bool par::epi_quickscan = false;

bool par::inbreeding = false;
bool par::check_sex = false;
bool par::impute_sex = false;
double par::sex_threshold_male = 0.8;
double par::sex_threshold_female = 0.2;

bool par::homo_run = false;
bool par::homo_run_snps = false;
bool par::homo_run_kb = false;
double par::homo_run_density = 50;
int par::homo_run_gap = 1000;  // 1Mb
bool par::homo_miss_as_hom = false;

int par::homo_run_length_snps= 100;
int par::homo_run_length_kb = 1000; // 1Mb in kb
int par::homo_run_het = 1;

int par::homo_windowSize = 50;  // SNPs
int par::homo_windowKB = 5000;  // 
int par::homo_windowAllowedHet = 1; // 1 SNP per 20
int par::homo_windowAllowedMissing = 5;
double par::homo_threshold = 0.05;

bool par::homo_summary_allelic_match = false;
bool par::homo_run_consensus_match = false;
double par::fuzzy_homo = 0.99;
bool par::homozyg_verbose = false;
int par::pool_size_min = 2;

bool par::ibs_run = false;
int par::ibs_run_length_snps = 100;
int par::ibs_run_length_kb = 100;
double par::ibs_run_density = 0.01;  // 1 SNP per 100kb average
int par::ibs_inner_run_length_kb = 100;
int par::ibs_inner_run_length_snp = 20;
int par::ibs_join_kb  = 100;
int par::ibs_join_snp = 1;
int par::ibs_run_missing = 2;
int par::ibs_run_0 = 1;
int par::ibs_inter_snp_distance = 1000000; // units=bp, 1=Mb
bool par::ibs_2only = false;

bool par::miss_run = false;
int par::miss_run_length = 100;
bool par::miss_run_length_kb = false;
double par::miss_run_level = 0.80;

bool par::FIXED = false;
bool par::FIXED_p = false;
Z par::FIX_IBD;
double par::FIX_p = 0.5;

bool par::matrix = false;
bool par::distance_matrix = false;
bool par::cluster = false;
bool par::cluster_euclidean = false;
bool par::cluster_group_avg = false;
bool par::force_initial_cluster = false;
bool par::cluster_plot = false;
int  par::cluster_mds_dim = 2;
bool par::mds_by_individual = true;
bool par::genome_groups = false;
bool par::cluster_ibm_constraint = false;
double par::cluster_ibm_constraint_value = 0;
bool par::cluster_missing = false;
bool par::cluster_selcon = false;
string par::cluster_selcon_file = "plink.clst";
int par::max_cluster_N = -1;
double par:: merge_p = 0;
int par::ibstest_gap = 500000; // 500 kb
int par::max_cluster_size = 0;
int par::max_cluster_case = 0;
int par::max_cluster_control = 0;
bool par::cluster_on_phenotype = false;
bool par::cluster_on_mcc = false;
int par::min_neighbour = 1;
int par::max_neighbour = 10;
bool par::outlier_detection = false;
bool par::bmatch = false;
bool par::bmatch_usertype = false;
bool par::qmatch = false;
string par::bmatch_filename = "plink.bmatch";
string par::qmatch_filename = "plink.qmatch";
string par::bmatch_direction_filename = "plink.bm";
string par::qmatch_threshold_filename = "plink.qt";

bool par::include_cluster = false;
bool par::include_cluster_from_file = false;
string par::include_cluster_filename = "plink.clst";
int par::analyse_cluster = 0;

bool par::af_write = false;
bool par::af_count = false;
bool par::af_read = false;

bool par::ibd_read = false;
string par::ibd_file = "plink.genome";
bool par::ibd_read_minimal = false;
bool par::ibd_read_list = false;
string par::ibd_file_list = "plink.genome.list";

bool par::inc_write = false;
bool par::inc_read = false;
string par::inc_file = "plink.inc";
string par::af_file = "plink.frq";

bool par::locked = false;

bool par::include_all_pairs = false;
double par::include_all_z1 = 0.001;

double par::MIN_PIHAT = 0.0025;
double par::MAX_PIHAT = 1.0000;
double par::MAX_CORR_PIHAT_PIHAT_G = 0.9;
double par::MAX_GENO_MISSING = 0.1;
double par::MAX_IND_MISSING = 0.1;
int par::MAX_LINE_LENGTH = 1000000;

bool par::remove_indiv = false;
bool par::keep_indiv = false;

bool par::extract_before_exclude = true;
bool par::remove_before_keep = true;
string par::remove_indiv_list = "plink.list";
string par::keep_indiv_list = "plink.list";

int par::pp_maxsnp = 6;
int par::pp_maxfid = 6;
int par::pp_maxiid = 6;

int par::BATCH_SIZE = 500000;

bool par::plink = false;

bool par::display_segment_long = false;
bool par::cnv_makemap = false;
bool par::cnv_writelist = false;
bool par::cnv_list = false;
bool par::display_cnv_track = false;
int par::cnv_col = 0;
string par::cnv_listname = "plink.cnv";
int par::cnv_min_kb = -1;
double par::cnv_min_score = -1;
int par::cnv_min_sites = -1;
int par::cnv_max_kb = -1;
double par::cnv_max_score = -1;
int par::cnv_max_sites = -1;
bool par::cnv_del_only = false;
bool par::cnv_dup_only = false;
int par::cnv_type = -1;
bool par::cnv_intersect = false;
bool par::cnv_exclude = false;
string par::cnv_intersect_file = "plink.file";
bool par::cnv_intersect_subset = false;
string par::cnv_intersect_subset_file = "plink.file";
double par::cnv_overlap = -1;
bool par::cnv_count = false;
bool par::cnv_defined_overlap = false;
bool par::cnv_indiv_perm = false;
bool par::cnv_pos_perm = false;
bool par::cnv_drop_no_segment = false;
bool par::cnv_freq_method2 = false;
double par::cnv_freq_method2_threshold = 0.8;
bool par::cnv_write_freq = false;
bool par::cnv_freq_include = false;
bool par::cnv_freq_include_below = true;
bool par::cnv_freq_include_exact = false;
bool par::cnv_freq_include_exact_exclude = false;
int par::cnv_freq_include_cnt = -1;
bool par::cnv_unique = false;
bool par::cnv_intersect_writeback = false;
bool par::cnv_intersect_writeback_verbose = false;
bool par::cnv_disrupt = false;
int par::cnv_region_border = 0; // kb
bool par::cnv_union_overlap = false;
bool par::cnv_region_overlap = false;
bool par::cnv_check_overlap = false;
bool par::cnv_count_baseline = false;
string par::cnv_count_baseline_file = "";
bool par::cnv_weighted_gene_test = false;
bool par::cnv_enrichment_test = false;
int par::cnv_en_model = 4;
bool par::cnv_glm = false;

bool par::seg_test_window = false;
double par::seg_test_window_bp = 100000;
bool par::seg_test_region = false;

bool par::dosage_assoc = false;
string par::dosage_file = "";
bool par::dosage_hard_call = false;
double par::dosage_hard_call_thresh = 0.99;
int par::dosage_hard_call_thresh2 = 0;
bool par::dosage_hasMap = false;
bool par::write_dosage = false;
