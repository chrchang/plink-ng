int sort_item_ids_nx(char** sorted_ids_ptr, int** id_map_ptr, int item_ct, char* item_ids, unsigned long max_id_len);

int sort_item_ids(char** sorted_ids_ptr, int** id_map_ptr, int unfiltered_ct, unsigned long* exclude_arr, unsigned int exclude_ct, char* item_ids, unsigned long max_id_len, int(* comparator_deref)(const void*, const void*), int* duplicate_fail_ptr);

int indiv_major_to_snp_major(char* indiv_major_fname, char* outname, FILE** outfile_ptr, int unfiltered_marker_ct);

int load_bim(FILE** bimfile_ptr, char* mapname, int* map_cols_ptr, unsigned int* unfiltered_marker_ct_ptr, unsigned int* marker_exclude_ct_ptr, unsigned long* max_marker_id_len_ptr, unsigned int* plink_maxsnp_ptr, unsigned long** marker_exclude_ptr, double** set_allele_freqs_ptr, char** marker_alleles_ptr, char** marker_ids_ptr, Chrom_info* chrom_info_ptr, unsigned int** marker_pos_ptr, char* extractname, char* excludename, char* freqname, char* refalleles, int calculation_type, int recode_modifier, int allelexxxx, int marker_pos_start, int marker_pos_end, unsigned int snp_window_size, char* markername_from, char* markername_to, char* markername_snp, char* sf_markers, unsigned char* sf_starts_range, unsigned int sf_ct, unsigned int sf_max_len, int* map_is_unsorted_ptr);

int sort_and_write_bim(int** map_reverse_ptr, FILE* mapfile, int map_cols, FILE** bimfile_ptr, char* outname, char* outname_end, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_ct, unsigned long max_marker_id_len, int species, char* marker_alleles);

int make_bed(FILE* bedfile, int bed_offset, FILE* bimfile, int map_cols, FILE** bedoutfile_ptr, FILE** famoutfile_ptr, FILE** bimoutfile_ptr, char* outname, char* outname_end, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_ct, char* marker_alleles, unsigned long* marker_reverse, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int indiv_ct, char* person_ids, unsigned int max_person_id_len, char* paternal_ids, unsigned int max_paternal_id_len, char* maternal_ids, unsigned int max_maternal_id_len, unsigned char* sex_info, char* pheno_c, double* pheno_d, double missing_phenod, char* output_missing_pheno, unsigned long max_marker_id_len, int map_is_unsorted, unsigned int* indiv_sort_map, int species);

int load_fam(FILE* famfile, unsigned long buflen, int fam_col_1, int fam_col_34, int fam_col_5, int fam_col_6, int true_fam_col_6, int missing_pheno, int missing_pheno_len, int affection_01, unsigned int* unfiltered_indiv_ct_ptr, char** person_ids_ptr, unsigned int* max_person_id_len_ptr, char** paternal_ids_ptr, unsigned int* max_paternal_id_len_ptr, char** maternal_ids_ptr, unsigned int* max_maternal_id_len_ptr, unsigned char** sex_info_ptr, int* affection_ptr, char** pheno_c_ptr, double** pheno_d_ptr, unsigned long** founder_info_ptr, unsigned long** indiv_exclude_ptr, int binary_files, unsigned long long** line_locs_ptr, unsigned long long** line_mids_ptr, int* pedbuflen_ptr);

int ped_to_bed(char* pedname, char* mapname, char* outname, int fam_col_1, int fam_col_34, int fam_col_5, int fam_col_6, int affection_01, int missing_pheno, Chrom_info* chrom_info_ptr);

int lgen_to_bed(char* lgen_namebuf, char* outname, int missing_pheno, int affection_01, Chrom_info* chrom_info_ptr);

int transposed_to_bed(char* tpedname, char* tfamname, char* outname, char missing_geno, Chrom_info* chrom_info_ptr);

int generate_dummy(char* outname, unsigned int dummy_flags, unsigned int dummy_marker_ct, unsigned int dummy_indiv_ct, double dummy_missing_geno, double dummy_missing_pheno);

int recode(int recode_modifier, FILE* bedfile, int bed_offset, FILE* famfile, FILE* bimfile, FILE** outfile_ptr, char* outname, char* outname_end, unsigned int unfiltered_marker_ct, unsigned long* marker_exclude, unsigned int marker_ct, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int indiv_ct, char* marker_ids, unsigned long max_marker_id_len, char* marker_alleles, char* person_ids, unsigned int max_person_id_len, char* paternal_ids, unsigned int max_paternal_id_len, char* maternal_ids, unsigned int max_maternal_id_len, unsigned char* sex_info, char* pheno_c, double* pheno_d, double missing_phenod, char output_missing_geno, char* output_missing_pheno, double* set_allele_freqs, int binary_files);

int merge_datasets(char* bedname, char* bimname, char* famname, char* outname, char* outname_end, char* mergename1, char* mergename2, char* mergename3, int calculation_type, int merge_type, int indiv_sort, int keep_allele_order, Chrom_info* chrom_info_ptr);
