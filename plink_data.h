#ifndef __PLINK_DATA_H__
#define __PLINK_DATA_H__

int32_t load_pheno(FILE* phenofile, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_exclude_ct, char* sorted_person_ids, uintptr_t max_person_id_len, uint32_t* id_map, int32_t missing_pheno, uint32_t missing_pheno_len, uint32_t affection_01, uint32_t mpheno_col, char* phenoname_str, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, char* phenoname_load, uintptr_t max_pheno_name_len);

int32_t sort_item_ids_nx(char** sorted_ids_ptr, uint32_t** id_map_ptr, uintptr_t item_ct, char* item_ids, uintptr_t max_id_len);

int32_t indiv_major_to_snp_major(char* indiv_major_fname, char* outname, uintptr_t unfiltered_marker_ct);

int32_t load_bim(char* bimname, uint32_t* map_cols_ptr, uintptr_t* unfiltered_marker_ct_ptr, uintptr_t* marker_exclude_ct_ptr, uintptr_t* max_marker_id_len_ptr, uintptr_t** marker_exclude_ptr, double** set_allele_freqs_ptr, char*** marker_allele_pp, uintptr_t* max_marker_allele_len_ptr, char** marker_ids_ptr, Chrom_info* chrom_info_ptr, double** marker_cms_ptr, uint32_t** marker_pos_ptr, char* freqname, uint64_t calculation_type, uint64_t misc_flags, uint32_t recode_modifier, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, Range_list* sf_range_list, uint32_t* map_is_unsorted_ptr, uint32_t marker_pos_needed, uint32_t marker_cms_needed, uint32_t marker_alleles_needed, const char* extension, const char* split_chrom_cmd);

int32_t apply_cm_map(char* cm_map_fname, char* cm_map_chrname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, double* marker_cms, Chrom_info* chrom_info_ptr);

int32_t update_marker_cms(Two_col_params* update_cm, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, double* marker_cms);

int32_t update_marker_pos(Two_col_params* update_map, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t* marker_pos, uint32_t* map_is_unsorted_ptr, Chrom_info* chrom_info_ptr);

int32_t update_marker_names(Two_col_params* update_name, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char* true_marker_ids);

int32_t update_marker_alleles(char* update_alleles_fname, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, char* outname, char* outname_end);

int32_t update_indiv_ids(char* update_ids_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, char* person_ids);

int32_t update_indiv_parents(char* update_parents_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info);

int32_t update_indiv_sexes(char* update_sex_fname, char* sorted_person_ids, uintptr_t indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, uintptr_t* sex_nm, uintptr_t* sex_male);

int32_t flip_strand(char* flip_fname, char* sorted_marker_ids, uintptr_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char** marker_allele_ptrs);

int32_t include_or_exclude(char* fname, char* sorted_ids, uintptr_t sorted_ids_len, uintptr_t max_id_len, uint32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr, uint32_t flags);

int32_t filter_attrib(char* fname, char* condition_str, char* sorted_ids, uintptr_t sorted_ids_ct, uintptr_t max_id_len, uint32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr, uint32_t is_indiv);

int32_t read_dists(char* dist_fname, char* id_fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* indiv_to_cluster, uint32_t for_cluster_flag, uint32_t is_max_dist, double* dists, uint32_t neighbor_n2, double* neighbor_quantiles, uint32_t* neighbor_qindices);

int32_t load_covars(char* covar_fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, double missing_phenod, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t gxe_mcovar, uintptr_t* covar_ct_ptr, char** covar_names_ptr, uintptr_t* max_covar_name_len_ptr, uintptr_t* pheno_nm, uintptr_t** covar_nm_ptr, double** covar_d_ptr, uintptr_t** gxe_covar_nm_ptr, uintptr_t** gxe_covar_c_ptr);

int32_t write_covars(char* outname, char* outname_end, uint32_t write_covar_modifier, uint32_t write_covar_dummy_max_categories, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, double missing_phenod, char* output_missing_pheno, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d);

int32_t make_bed(FILE* bedfile, uintptr_t bed_offset, char* bimname, uint32_t map_cols, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t map_is_unsorted, uint32_t* indiv_sort_map, uint64_t misc_flags, Two_col_params* update_chr, char* flip_subset_fname, uint32_t hh_exists, Chrom_info* chrom_info_ptr);

int32_t load_fam(FILE* famfile, uint32_t buflen, uint32_t fam_cols, uint32_t tmp_fam_col_6, int32_t missing_pheno, uint32_t missing_pheno_len, uint32_t affection_01, uintptr_t* unfiltered_indiv_ct_ptr, char** person_ids_ptr, uintptr_t* max_person_id_len_ptr, char** paternal_ids_ptr, uintptr_t* max_paternal_id_len_ptr, char** maternal_ids_ptr, uintptr_t* max_maternal_id_len_ptr, uintptr_t** sex_nm_ptr, uintptr_t** sex_male_ptr, uint32_t* affection_ptr, uintptr_t** pheno_nm_ptr, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, uintptr_t** founder_info_ptr, uintptr_t** indiv_exclude_ptr);

int32_t oxford_to_bed(char* genname, char* samplename, char* outname, char* outname_end, double hard_call_threshold, char* missing_code, int32_t missing_pheno, uint64_t misc_flags, Chrom_info* chrom_info_ptr);

int32_t ped_to_bed(char* pedname, char* mapname, char* outname, char* outname_end, uint32_t fam_cols, uint32_t affection_01, int32_t missing_pheno, Chrom_info* chrom_info_ptr);

int32_t lgen_to_bed(char* lgen_namebuf, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, uint32_t lgen_modifier, char* lgen_reference_fname, Chrom_info* chrom_info_ptr);

int32_t transposed_to_bed(char* tpedname, char* tfamname, char* outname, char* outname_end, uint64_t misc_flags, Chrom_info* chrom_info_ptr);

int32_t vcf_to_bed(char* vcfname, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, char* const_fid, char id_delim, double vcf_min_qual, char* vcf_filter_exceptions_flattened, Chrom_info* chrom_info_ptr);

int32_t bcf_to_bed(char* bcfname, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, char* const_fid, char id_delim, double vcf_min_qual, char* vcf_filter_exceptions_flattened, Chrom_info* chrom_info_ptr);

int32_t bed_from_23(char* fname, char* outname, char* outname_end, uint32_t modifier_23, char* fid_23, char* iid_23, double pheno_23, char* paternal_id_23, char* maternal_id_23, char* convert_xy_fname, Chrom_info* chrom_info_ptr);

int32_t generate_dummy(char* outname, char* outname_end, uint32_t flags, uintptr_t marker_ct, uintptr_t indiv_ct, double geno_mrate, double pheno_mrate);

int32_t simulate_dataset(char* outname, char* outname_end, uint32_t flags, char* simulate_fname, uint32_t case_ct, uint32_t ctrl_ct, double prevalence, uint32_t qt_indiv_ct, double missing_freq, char* name_prefix);

int32_t recode(uint32_t recode_modifier, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, char* recode_allele_name, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uint32_t* marker_pos, uintptr_t* marker_reverse, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t map_is_unsorted, uint64_t misc_flags, uint32_t hh_exists, Chrom_info* chrom_info_ptr);

int32_t merge_datasets(char* bedname, char* bimname, char* famname, char* outname, char* outname_end, char* mergename1, char* mergename2, char* mergename3, uint64_t calculation_type, uint32_t merge_type, uint32_t indiv_sort, uint64_t misc_flags, Chrom_info* chrom_info_ptr);

#endif // __PLINK_DATA_H__
