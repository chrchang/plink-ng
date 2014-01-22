#ifndef __PLINK_FILTER_H__
#define __PLINK_FILTER_H__

#include "plink_common.h"

typedef struct {
  uintptr_t cluster_ct;
  uintptr_t entry_ct;
  uint64_t* entries;
  uintptr_t* cluster_zmasks;
  uint32_t* cluster_ref_cts;
  uint32_t* indiv_lookup;
  char* marker_fname;
  char* indiv_fname;
} Oblig_missing_info;

void oblig_missing_init(Oblig_missing_info* om_ip);

void oblig_missing_cleanup(Oblig_missing_info* om_ip);

int32_t load_oblig_missing(FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, char* marker_ids, uintptr_t max_marker_id_len, char* sorted_person_ids, uintptr_t sorted_indiv_ct, uintptr_t max_person_id_len, uint32_t* indiv_id_map, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* sex_male, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip);

int32_t filter_indivs_file(char* filtername, char* sorted_person_ids, uintptr_t sorted_ids_len, uintptr_t max_person_id_len, uint32_t* id_map, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* filtervals_flattened, uint32_t mfilter_col);

void filter_indivs_bitfields(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, uintptr_t* orfield, int32_t orfield_flip, uintptr_t* ornot);

int32_t mind_filter(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double mind_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_exclude_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* person_ids, uintptr_t max_person_id_len, uintptr_t* sex_male, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip);

int32_t calc_freqs_and_hwe(FILE* bedfile, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_exclude_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t* founder_info, int32_t nonfounders, int32_t maf_succ, double* set_allele_freqs, uintptr_t** marker_reverse_ptr, uint32_t** marker_allele_cts_ptr, uintptr_t bed_offset, uint32_t hwe_or_geno_needed, uint32_t hwe_all, uint32_t hardy_needed, uintptr_t* pheno_nm, uintptr_t* pheno_c, int32_t** hwe_lls_ptr, int32_t** hwe_lhs_ptr, int32_t** hwe_hhs_ptr, int32_t** hwe_ll_cases_ptr, int32_t** hwe_lh_cases_ptr, int32_t** hwe_hh_cases_ptr, int32_t** hwe_ll_allfs_ptr, int32_t** hwe_lh_allfs_ptr, int32_t** hwe_hh_allfs_ptr, int32_t** hwe_hapl_allfs_ptr, int32_t** hwe_haph_allfs_ptr, uint32_t* indiv_male_ct_ptr, uint32_t* indiv_f_ct_ptr, uint32_t* indiv_f_male_ct_ptr, uint32_t wt_needed, uintptr_t* topsize_ptr, double** marker_weights_ptr, double exponent, Chrom_info* chrom_info_ptr, uintptr_t* sex_nm, uintptr_t* sex_male, uint32_t is_split_chrom, uint32_t* hh_exists_ptr);

int32_t write_missingness_reports(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, uintptr_t* indiv_exclude, uintptr_t* pheno_nm, uintptr_t* sex_male, uint32_t indiv_male_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t hh_exists);

uintptr_t binary_geno_filter(double geno_thresh, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t indiv_ct, uintptr_t male_ct, int32_t* hwe_lls, int32_t* hwe_lhs, int32_t* hwe_hhs, Chrom_info* chrom_info_ptr, Oblig_missing_info* om_ip);

#endif // __PLINK_ASSOC_H__
