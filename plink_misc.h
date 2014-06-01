#ifndef __PLINK_MISC_H__
#define __PLINK_MISC_H__

#define SCORE_HEADER 1
#define SCORE_SUM 2
#define SCORE_NO_MEAN_IMPUTATION 4
#define SCORE_CENTER 8
#define SCORE_DATA_HEADER 0x10

typedef struct {
  char* fname;
  char* range_fname;
  char* data_fname;
  uint32_t modifier;
  uint32_t varid_col;
  uint32_t allele_col;
  uint32_t effect_col;
  uint32_t data_varid_col;
  uint32_t data_col;
} Score_info;

void misc_init(Score_info* sc_ip);

void misc_cleanup(Score_info* sc_ip);

int32_t make_founders(uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t require_two, uintptr_t* indiv_exclude, uintptr_t* founder_info);

int32_t write_nosex(char* outname, char* outname_end, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* sex_nm, uintptr_t gender_unk_ct, char* person_ids, uintptr_t max_person_id_len);

int32_t makepheno_load(FILE* phenofile, char* makepheno_str, uintptr_t unfiltered_indiv_ct, char* sorted_person_ids, uintptr_t max_person_id_len, uint32_t* id_map, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr);

int32_t load_pheno(FILE* phenofile, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_exclude_ct, char* sorted_person_ids, uintptr_t max_person_id_len, uint32_t* id_map, int32_t missing_pheno, uint32_t affection_01, uint32_t mpheno_col, char* phenoname_str, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, char* phenoname_load, uintptr_t max_pheno_name_len);

int32_t convert_tail_pheno(uint32_t unfiltered_indiv_ct, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, double tail_bottom, double tail_top, double missing_phenod);

int32_t apply_cm_map(char* cm_map_fname, char* cm_map_chrname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, double* marker_cms, Chrom_info* chrom_info_ptr);

int32_t update_marker_cms(Two_col_params* update_cm, char* sorted_marker_ids, uintptr_t sorted_ids_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, double* marker_cms);

int32_t update_marker_pos(Two_col_params* update_map, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t* marker_pos, uint32_t* map_is_unsorted_ptr, Chrom_info* chrom_info_ptr);

int32_t update_marker_names(Two_col_params* update_name, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char* true_marker_ids);

int32_t update_marker_alleles(char* update_alleles_fname, char* sorted_marker_ids, uintptr_t sorted_ids_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, char* outname, char* outname_end);

int32_t flip_strand(char* flip_fname, char* sorted_marker_ids, uintptr_t sorted_ids_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char** marker_allele_ptrs);

int32_t update_indiv_ids(char* update_ids_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, char* person_ids);

int32_t update_indiv_parents(char* update_parents_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info);

int32_t update_indiv_sexes(char* update_sex_fname, uint32_t update_sex_col, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, uintptr_t* sex_nm, uintptr_t* sex_male);

uint32_t calc_plink_maxsnp(uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len);

int32_t read_external_freqs(char* freqname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, char** marker_allele_ptrs, double* set_allele_freqs, uint32_t maf_succ, double exponent, uint32_t wt_needed, double* marker_weights);

int32_t load_ax_alleles(Two_col_params* axalleles, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, double* set_allele_freqs, uint32_t is_a2);

int32_t write_stratified_freqs(FILE* bedfile, uintptr_t bed_offset, char* outname, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uint32_t indiv_f_ct, uintptr_t* founder_info, uint32_t nonfounders, uintptr_t* sex_male, uint32_t indiv_f_male_ct, uintptr_t* marker_reverse, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len);

int32_t write_freqs(char* outname, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, double* set_allele_freqs, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, int32_t* hapl_cts, int32_t* haph_cts, uint32_t indiv_f_ct, uint32_t indiv_f_male_ct, uint64_t misc_flags, uintptr_t* marker_reverse);

void calc_marker_weights(double exponent, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, double* marker_weights);

int32_t sexcheck(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uint64_t misc_flags, double check_sex_fthresh, double check_sex_mthresh, uint32_t max_f_yobs, uint32_t min_m_yobs, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* gender_unk_ct_ptr);

int32_t write_snplist(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t list_23_indels);

int32_t het_report(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs);

int32_t score_report(Score_info* sc_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, double* set_allele_freqs, uintptr_t indiv_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t hh_exists, Chrom_info* chrom_info_ptr, char* outname, char* outname_end);

#endif // __PLINK_MISC_H__
