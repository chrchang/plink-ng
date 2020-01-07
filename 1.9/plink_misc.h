#ifndef __PLINK_MISC_H__
#define __PLINK_MISC_H__

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


#define SCORE_HEADER 1
#define SCORE_SUM 2
#define SCORE_NO_MEAN_IMPUTATION 4
#define SCORE_CENTER 8
#define SCORE_DATA_HEADER 0x10

#define METAANAL_STUDY 1
#define METAANAL_NO_MAP 2
#define METAANAL_NO_ALLELE 4
#define METAANAL_REPORT_ALL 8
#define METAANAL_LOGSCALE 0x10
#define METAANAL_QT 0x20
#define METAANAL_WEIGHTED_Z 0x40
#define METAANAL_REPORT_DUPS 0x80

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

int32_t make_founders(uintptr_t unfiltered_sample_ct, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t require_two, uintptr_t* sample_exclude, uintptr_t* founder_info);

int32_t write_nosex(char* outname, char* outname_end, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sex_nm, uintptr_t gender_unk_ct, char* sample_ids, uintptr_t max_sample_id_len);

int32_t makepheno_load(FILE* phenofile, char* makepheno_str, uintptr_t unfiltered_sample_ct, char* sorted_sample_ids, uintptr_t max_sample_id_len, uint32_t* id_map, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr);

int32_t load_pheno(FILE* phenofile, uintptr_t unfiltered_sample_ct, uintptr_t sample_exclude_ct, char* sorted_sample_ids, uintptr_t max_sample_id_len, uint32_t* id_map, int32_t missing_pheno, uint32_t affection_01, uint32_t mpheno_col, char* phenoname_str, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, char* phenoname_load, uintptr_t max_pheno_name_len);

int32_t convert_tail_pheno(uint32_t unfiltered_sample_ct, uintptr_t* pheno_nm, uintptr_t** pheno_c_ptr, double** pheno_d_ptr, double tail_bottom, double tail_top, double missing_phenod);

int32_t apply_cm_map(char* cm_map_fname, char* cm_map_chrname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, double* marker_cms, Chrom_info* chrom_info_ptr);

int32_t update_marker_cms(Two_col_params* update_cm, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, double* marker_cms);

int32_t update_marker_pos(Two_col_params* update_map, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t* marker_pos, uint32_t* map_is_unsorted_ptr, Chrom_info* chrom_info_ptr);

int32_t update_marker_names(Two_col_params* update_name, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct);

int32_t update_marker_alleles(char* update_alleles_fname, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, char* outname, char* outname_end);

int32_t flip_strand(char* flip_fname, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char** marker_allele_ptrs);

int32_t update_sample_ids(char* update_ids_fname, char* sorted_sample_ids, uintptr_t sample_ct, uintptr_t max_sample_id_len, uint32_t* sample_id_map, char* sample_ids);

int32_t update_sample_parents(char* update_parents_fname, char* sorted_sample_ids, uintptr_t sample_ct, uintptr_t max_sample_id_len, uint32_t* sample_id_map, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info);

int32_t update_sample_sexes(char* update_sex_fname, uint32_t update_sex_col, char* sorted_sample_ids, uintptr_t sample_ct, uintptr_t max_sample_id_len, uint32_t* sample_id_map, uintptr_t* sex_nm, uintptr_t* sex_male);

void calc_plink_maxfid(uint32_t unfiltered_sample_ct, uintptr_t* sample_exclude, uint32_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uint32_t* plink_maxfid_ptr, uint32_t* plink_maxiid_ptr);

uint32_t calc_plink_maxsnp(uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len);

int32_t read_external_freqs(char* freqname, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, char** marker_allele_ptrs, double* set_allele_freqs, uint32_t* nchrobs, uint32_t maf_succ);

int32_t load_ax_alleles(Two_col_params* axalleles, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char** marker_allele_ptrs, uintptr_t* max_marker_allele_len_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, double* set_allele_freqs, uint32_t is_a2);

int32_t write_stratified_freqs(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t output_gz, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_sample_ct, uintptr_t sample_ct, uint32_t sample_f_ct, uintptr_t* founder_info, uint32_t nonfounders, uintptr_t* sex_male, uint32_t sample_f_male_ct, uintptr_t* marker_reverse, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len);

int32_t write_cc_freqs(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t output_gz, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_sample_ct, uintptr_t* founder_info, uint32_t nonfounders, uintptr_t* sex_male, uintptr_t* marker_reverse, uintptr_t* pheno_nm, uintptr_t* pheno_c);

int32_t write_freqs(char* outname, char* outname_end, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, double* set_allele_freqs, Chrom_info* chrom_info_ptr, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, int32_t* ll_cts, int32_t* lh_cts, int32_t* hh_cts, int32_t* hapl_cts, int32_t* haph_cts, uint32_t sample_f_ct, uint32_t sample_f_male_ct, uint32_t nonfounders, uint64_t misc_flags, uintptr_t* marker_reverse);

int32_t sexcheck(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, uint64_t misc_flags, double check_sex_fthresh, double check_sex_mthresh, uint32_t max_f_yobs, uint32_t min_m_yobs, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* gender_unk_ct_ptr);

int32_t write_snplist(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t list_23_indels);

int32_t write_var_ranges(char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t write_var_range_ct);

int32_t list_duplicate_vars(char* outname, char* outname_end, uint32_t dupvar_modifier, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, Chrom_info* chrom_info_ptr, char** marker_allele_ptrs);

int32_t het_report(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t output_gz, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uintptr_t* founder_info, Chrom_info* chrom_info_ptr, double* set_allele_freqs);

int32_t fst_report(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts);

int32_t score_report(Score_info* sc_ip, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, double* set_allele_freqs, uintptr_t sample_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, char* output_missing_pheno, uint32_t hh_exists, Chrom_info* chrom_info_ptr, char* outname, char* outname_end);

int32_t meta_analysis(char* input_fnames, char* chrfield_search_order, char* snpfield_search_order, char* bpfield_search_order, char* a1field_search_order, char* a2field_search_order, char* pfield_search_order, char* sefield_search_order, char* essfield_search_order, uint32_t flags, char* extractname, char* outname, char* outname_end, double output_min_p, Chrom_info* chrom_info_ptr);

#endif // __PLINK_MISC_H__
