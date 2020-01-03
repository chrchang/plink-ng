#ifndef __PLINK_FILTER_H__
#define __PLINK_FILTER_H__

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


#include "plink_common.h"

typedef struct {
  uintptr_t cluster_ct;
  uintptr_t entry_ct;
  uint64_t* entries;
  uint32_t* cluster_ref_cts;
  uint32_t* sample_lookup;
  char* marker_fname;
  char* sample_fname;
} Oblig_missing_info;

void oblig_missing_init(Oblig_missing_info* om_ip);

void oblig_missing_cleanup(Oblig_missing_info* om_ip);

int32_t keep_or_remove(char* fname, char* sorted_ids, uintptr_t sorted_ids_len, uintptr_t max_id_len, uint32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr, uint32_t flags, uint32_t allow_no_samples);

int32_t snps_flag(const char* variant_ids, const uint32_t* variant_id_htable, const Range_list* snps_range_list_ptr, uint32_t raw_variant_ct, uintptr_t max_variant_id_blen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_exclude, uintptr_t* exclude_ct_ptr);

int32_t extract_exclude_flag_norange(char* fname, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, uint32_t do_exclude, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t allow_no_variants);

int32_t filter_attrib(char* fname, char* condition_str, uint32_t* id_htable, uint32_t id_htable_size, uint32_t allow_no_variants, char* item_ids, uintptr_t max_id_len, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr);

int32_t filter_attrib_sample(char* fname, char* condition_str, char* sorted_ids, uintptr_t sorted_ids_ct, uintptr_t max_id_len, uint32_t* id_map, uintptr_t unfiltered_ct, uint32_t allow_no_samples, uintptr_t* exclude_arr, uintptr_t* exclude_ct_ptr);

int32_t filter_qual_scores(Two_col_params* qual_filter, double qual_min_thresh, double qual_max_thresh, uint32_t* marker_id_htable, uint32_t marker_id_htable_size, uint32_t allow_no_variants, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr);

uint32_t random_thin_markers(double thin_keep_prob, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t allow_no_variants);

int32_t random_thin_markers_ct(uint32_t thin_keep_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr);

uint32_t random_thin_samples(double thin_keep_prob, uintptr_t unfiltered_sample_ct, uint32_t allow_no_samples, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr);

int32_t random_thin_samples_ct(uint32_t thin_keep_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr);

int32_t load_oblig_missing(FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, char* sorted_sample_ids, uintptr_t sorted_sample_ct, uintptr_t max_sample_id_len, uint32_t* sample_id_map, uintptr_t unfiltered_sample_ct, uintptr_t* sex_male, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip);

int32_t filter_samples_file(char* filtername, char* sorted_sample_ids, uintptr_t sorted_ids_len, uintptr_t max_sample_id_len, uint32_t* id_map, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, char* filtervals_flattened, uint32_t mfilter_col, uint32_t allow_no_samples);

void filter_samples_bitfields(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, uintptr_t* orfield, int32_t orfield_flip, uintptr_t* ornot);

int32_t mind_filter(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double mind_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* sex_male, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, uint32_t allow_no_samples);

int32_t calc_freqs_and_hwe(FILE* bedfile, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_exclude_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* founder_info, int32_t nonfounders, int32_t maf_succ, double* set_allele_freqs, uintptr_t bed_offset, uint32_t hwe_needed, uint32_t hwe_all, uint32_t hardy_needed, uint32_t min_ac, uint32_t max_ac, double geno_thresh, uintptr_t* pheno_nm, uintptr_t* pheno_c, int32_t** hwe_lls_ptr, int32_t** hwe_lhs_ptr, int32_t** hwe_hhs_ptr, int32_t** hwe_ll_cases_ptr, int32_t** hwe_lh_cases_ptr, int32_t** hwe_hh_cases_ptr, int32_t** hwe_ll_allfs_ptr, int32_t** hwe_lh_allfs_ptr, int32_t** hwe_hh_allfs_ptr, int32_t** hwe_hapl_allfs_ptr, int32_t** hwe_haph_allfs_ptr, uintptr_t** geno_excl_bitfield_ptr, uintptr_t** ac_excl_bitfield_ptr, uint32_t* sample_male_ct_ptr, uint32_t* sample_f_ct_ptr, uint32_t* sample_f_male_ct_ptr, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t is_split_chrom, uint32_t* hh_exists_ptr);

int32_t write_missingness_reports(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t output_gz, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_sample_ct, uintptr_t sample_ct, uintptr_t* sample_exclude, uintptr_t* pheno_nm, uintptr_t* sex_male, uint32_t sample_male_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t hh_exists);

int32_t hardy_report(char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, uint32_t hwe_modifier, uint32_t nonfounders, int32_t* hwe_ll_cases, int32_t* hwe_lh_cases, int32_t* hwe_hh_cases, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs, uint32_t pheno_nm_ct, uintptr_t* pheno_c, Chrom_info* chrom_info_ptr);

uint32_t enforce_hwe_threshold(double hwe_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, uint32_t hwe_modifier, uint32_t allow_no_variants, int32_t* hwe_ll_allfs, int32_t* hwe_lh_allfs, int32_t* hwe_hh_allfs, Chrom_info* chrom_info_ptr);

uint32_t enforce_minor_allele_thresholds(double min_maf, double max_maf, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* ac_excl_bitfield, uintptr_t* marker_exclude_ct_ptr, double* set_allele_freqs, uint32_t allow_no_variants);

void enforce_min_bp_space(int32_t min_bp_space, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, uintptr_t* marker_exclude_ct_ptr, Chrom_info* chrom_info_ptr);

#endif // __PLINK_FILTER_H__
