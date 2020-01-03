#ifndef __PLINK_DATA_H__
#define __PLINK_DATA_H__

// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


int32_t sample_major_to_snp_major(char* sample_major_fname, char* outname, uintptr_t unfiltered_marker_ct, uintptr_t unfiltered_sample_ct, uint64_t fsize);

int32_t load_bim(char* bimname, uintptr_t* unfiltered_marker_ct_ptr, uintptr_t* marker_exclude_ct_ptr, uintptr_t* max_marker_id_blen_ptr, uintptr_t** marker_exclude_ptr, double** set_allele_freqs_ptr, uint32_t** nchrobs_ptr, char*** marker_allele_pp, uintptr_t* max_marker_allele_blen_ptr, char** marker_ids_ptr, char* missing_mid_template, uint32_t new_id_max_allele_slen, const char* missing_marker_id_match, Chrom_info* chrom_info_ptr, double** marker_cms_ptr, uint32_t** marker_pos_ptr, uint64_t misc_flags, uint64_t filter_flags, int32_t marker_pos_start, int32_t marker_pos_end, int32_t snp_window_size, char* markername_from, char* markername_to, char* markername_snp, uint32_t* map_is_unsorted_ptr, uint32_t marker_pos_needed, uint32_t marker_cms_needed, uint32_t marker_alleles_needed, const char* split_chrom_cmd, const char* ftype_str, uint32_t* max_bim_linelen_ptr);

int32_t load_covars(char* covar_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, double missing_phenod, uint32_t covar_modifier, Range_list* covar_range_list_ptr, uint32_t gxe_mcovar, uintptr_t* covar_ctx_ptr, char** covar_names_ptr, uintptr_t* max_covar_name_len_ptr, uintptr_t* pheno_nm, uintptr_t** covar_nm_ptr, double** covar_d_ptr, uintptr_t** gxe_covar_nm_ptr, uintptr_t** gxe_covar_c_ptr);

int32_t write_covars(char* outname, char* outname_end, uint32_t write_covar_modifier, uint32_t write_covar_dummy_max_categories, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, double missing_phenod, char* output_missing_pheno, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d);

int32_t make_bed(FILE* bedfile, uintptr_t bed_offset, char* bimname, char* outname, char* outname_end, uint64_t calculation_type, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t map_is_unsorted, uint32_t* sample_sort_map, uint64_t misc_flags, uint32_t splitx_bound1, uint32_t splitx_bound2, Two_col_params* update_chr, char* flip_fname, char* flip_subset_fname, char* zerofname, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t hh_exists, Chrom_info* chrom_info_ptr, uint32_t mendel_modifier, uint32_t max_bim_linelen);

int32_t load_fam(char* famname, uint32_t fam_cols, uint32_t tmp_fam_col_6, int32_t missing_pheno, uint32_t affection_01, uintptr_t* unfiltered_sample_ct_ptr, char** sample_ids_ptr, uintptr_t* max_sample_id_len_ptr, char** paternal_ids_ptr, uintptr_t* max_paternal_id_len_ptr, char** maternal_ids_ptr, uintptr_t* max_maternal_id_len_ptr, uintptr_t** sex_nm_ptr, uintptr_t** sex_male_ptr, uint32_t* affection_ptr, uintptr_t** pheno_nm_ptr, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, uintptr_t** founder_info_ptr, uintptr_t** sample_exclude_ptr, uint32_t allow_no_samples);

int32_t oxford_to_bed(char* genname, char* samplename, char* outname, char* outname_end, char* single_chr, char* pheno_name, double hard_call_threshold, char* missing_code, int32_t missing_pheno, uint64_t misc_flags, uint32_t is_bgen, Chrom_info* chrom_info_ptr);

int32_t ped_to_bed(char* pedname, char* mapname, char* outname, char* outname_end, uint32_t fam_cols, uint64_t misc_flags, int32_t missing_pheno, Chrom_info* chrom_info_ptr);

uint32_t realpath_identical(const char* outname, const char* read_realpath, char* write_realpath_buf);

int32_t lgen_to_bed(char* lgenname, char* mapname, char* famname, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, uint32_t lgen_modifier, char* lgen_reference_fname, Chrom_info* chrom_info_ptr);

int32_t transposed_to_bed(char* tpedname, char* tfamname, char* outname, char* outname_end, uint64_t misc_flags, Chrom_info* chrom_info_ptr);

int32_t vcf_to_bed(char* vcfname, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, char* const_fid, char id_delim, char vcf_idspace_to, double vcf_min_qual, char* vcf_filter_exceptions_flattened, double vcf_min_gq, double vcf_min_gp, uint32_t vcf_half_call, Chrom_info* chrom_info_ptr);

int32_t bcf_to_bed(char* bcfname, char* outname, char* outname_end, int32_t missing_pheno, uint64_t misc_flags, char* const_fid, char id_delim, char vcf_idspace_to, double vcf_min_qual, char* vcf_filter_exceptions_flattened, Chrom_info* chrom_info_ptr);

int32_t bed_from_23(char* fname, char* outname, char* outname_end, uint32_t modifier_23, char* fid_23, char* iid_23, double pheno_23, uint64_t misc_flags, char* paternal_id_23, char* maternal_id_23, Chrom_info* chrom_info_ptr);

int32_t generate_dummy(char* outname, char* outname_end, uint32_t flags, uintptr_t marker_ct, uintptr_t sample_ct, double geno_mrate, double pheno_mrate, int32_t missing_pheno);

int32_t simulate_dataset(char* outname, char* outname_end, uint32_t flags, char* simulate_fname, uint32_t case_ct, uint32_t ctrl_ct, double prevalence, uint32_t qt_sample_ct, double missing_freq, char* name_prefix);

int32_t recode(uint32_t recode_modifier, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, char* recode_allele_name, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* marker_ids, uintptr_t max_marker_id_len, double* marker_cms, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uint32_t* marker_pos, uintptr_t* marker_reverse, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t map_is_unsorted, uint64_t misc_flags, uint32_t hh_exists, Chrom_info* chrom_info_ptr);

int32_t sample_sort_file_map(char* sample_sort_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uint32_t** sample_sort_map_ptr);

int32_t merge_datasets(char* bedname, char* bimname, char* famname, char* outname, char* outname_end, char* mergename1, char* mergename2, char* mergename3, char* sample_sort_fname, uint64_t calculation_type, uint32_t merge_type, uint32_t sample_sort, uint64_t misc_flags, Chrom_info* chrom_info_ptr);

#endif // __PLINK_DATA_H__
