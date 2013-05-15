int32_t sort_item_ids_nx(char** sorted_ids_ptr, int32_t** id_map_ptr, uintptr_t item_ct, char* item_ids, uintptr_t max_id_len);

int32_t sort_item_ids(char** sorted_ids_ptr, int32_t** id_map_ptr, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t exclude_ct, char* item_ids, uintptr_t max_id_len, int(* comparator_deref)(const void*, const void*), int32_t* duplicate_fail_ptr);

int32_t indiv_major_to_snp_major(char* indiv_major_fname, char* outname, FILE** outfile_ptr, uintptr_t unfiltered_marker_ct);

int32_t load_bim(FILE** bimfile_ptr, char* mapname, int32_t* map_cols_ptr, uintptr_t* unfiltered_marker_ct_ptr, uintptr_t* marker_exclude_ct_ptr, uintptr_t* max_marker_id_len_ptr, uint32_t* plink_maxsnp_ptr, uintptr_t** marker_exclude_ptr, double** set_allele_freqs_ptr, char** marker_alleles_ptr, uintptr_t* max_marker_allele_len_ptr, char** marker_ids_ptr, Chrom_info* chrom_info_ptr, double** marker_cms_ptr, uint32_t** marker_pos_ptr, char* freqname, uint64_t calculation_type, uint32_t recode_modifier, int32_t marker_pos_start, int32_t marker_pos_end, uint32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, char* sf_markers, unsigned char* sf_starts_range, uint32_t sf_ct, uint32_t sf_max_len, Two_col_params* update_map, int32_t* map_is_unsorted_ptr, uint32_t marker_cms_needed, uint32_t marker_alleles_needed, const char* extension);

int32_t update_marker_cms(Two_col_params* update_cm, char* sorted_marker_ids, uint32_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, double* marker_cms);

int32_t update_marker_names(Two_col_params* update_name, char* sorted_marker_ids, uint32_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char* true_marker_ids);

int32_t update_marker_alleles(char* update_alleles_fname, char* sorted_marker_ids, uint32_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char* marker_alelles, uintptr_t max_marker_allele_len, char* outname, char* outname_end);

int32_t flip_strand(char* flip_fname, char* sorted_marker_ids, uint32_t marker_ct, uintptr_t max_marker_id_len, uint32_t* marker_id_map, char* marker_alleles, uintptr_t max_marker_allele_len);

int32_t include_or_exclude(char* fname, char* sorted_ids, int32_t sorted_ids_len, uintptr_t max_id_len, int32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr, int32_t indivs, int32_t do_exclude, int32_t duplicates_present);

int32_t make_bed(FILE* bedfile, int32_t bed_offset, FILE* bimfile, int32_t map_cols, FILE** bedoutfile_ptr, FILE** famoutfile_ptr, FILE** bimoutfile_ptr, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, double* marker_cms, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, double missing_phenod, char* output_missing_pheno, uintptr_t max_marker_id_len, int32_t map_is_unsorted, uint32_t* indiv_sort_map, uint32_t set_hh_missing, Two_col_params* update_chr, Two_col_params* update_map, char* flip_subset_fname, Chrom_info* chrom_info_ptr);

int32_t load_fam(FILE* famfile, uint32_t buflen, uint32_t fam_cols, uint32_t tmp_fam_col_6, int32_t missing_pheno, int32_t missing_pheno_len, int32_t affection_01, uintptr_t* unfiltered_indiv_ct_ptr, char** person_ids_ptr, uintptr_t* max_person_id_len_ptr, char** paternal_ids_ptr, uintptr_t* max_paternal_id_len_ptr, char** maternal_ids_ptr, uintptr_t* max_maternal_id_len_ptr, uintptr_t** sex_nm_ptr, uintptr_t** sex_male_ptr, int32_t* affection_ptr, uintptr_t** pheno_nm_ptr, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, uintptr_t** founder_info_ptr, uintptr_t** indiv_exclude_ptr);

int32_t ped_to_bed(char* pedname, char* mapname, char* outname, char* outname_end, uint32_t fam_cols, int32_t affection_01, int32_t missing_pheno, Chrom_info* chrom_info_ptr);

int32_t lgen_to_bed(char* lgen_namebuf, char* outname, char* outname_end, int32_t missing_pheno, int32_t affection_01, uint32_t lgen_modifier, char* lgen_reference_fname, Chrom_info* chrom_info_ptr);

int32_t transposed_to_bed(char* tpedname, char* tfamname, char* outname, char* outname_end, char missing_geno, Chrom_info* chrom_info_ptr);

int32_t generate_dummy(char* outname, char* outname_end, uint32_t flags, uintptr_t marker_ct, uintptr_t indiv_ct, double geno_mrate, double pheno_mrate);

int32_t simulate_dataset(char* outname, char* outname_end, uint32_t flags, char* simulate_fname, uint32_t case_ct, uint32_t ctrl_ct, double prevalence, uint32_t qt_indiv_ct, double missing_freq, char* name_prefix);

int32_t recode(uint32_t recode_modifier, FILE* bedfile, int32_t bed_offset, FILE* famfile, FILE* bimfile, FILE** outfile_ptr, char* outname, char* outname_end, char* recode_allele_name, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, uint32_t* marker_pos, uintptr_t* marker_reverse, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, double missing_phenod, char output_missing_geno, char* output_missing_pheno, uint32_t set_hh_missing, uint32_t xmhh_exists, uint32_t nxmhh_exists, Chrom_info* chrom_info_ptr);

int32_t merge_datasets(char* bedname, char* bimname, char* famname, char* outname, char* outname_end, char* mergename1, char* mergename2, char* mergename3, uint64_t calculation_type, int32_t merge_type, int32_t indiv_sort, int32_t keep_allele_order, Chrom_info* chrom_info_ptr);
