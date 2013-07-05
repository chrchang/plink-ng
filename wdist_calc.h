#ifndef __WDIST_CALC_H__
#define __WDIST_CALC_H__

#include "wdist_cluster.h"

extern double* g_rel_dists;
extern uint32_t* g_indiv_missing_unwt;
extern uint32_t* g_missing_dbl_excluded;
extern double* g_dists;
extern uint32_t* g_genome_main;

// int32_t invert_matrix(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf);

// int32_t invert_matrix_trunc_singular(__CLPK_integer dim, double* matrix, MATRIX_INVERT_BUF1_TYPE* int_1d_buf, double* dbl_2d_buf, __CLPK_integer min_dim);

int32_t regress_rel_main(uintptr_t* indiv_exclude, uintptr_t indiv_ct, uintptr_t regress_rel_iters, uint32_t regress_rel_d, pthread_t* threads, double* pheno_d);

#ifndef NOLAPACK
int32_t calc_unrelated_herit(uint64_t calculation_type, int32_t ibc_type, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, double* pheno_d, double* rel_ibc, uint32_t is_strict, double unrelated_herit_covg, double unrelated_herit_covr, double unrelated_herit_tol);

int32_t unrelated_herit_batch(uint32_t load_grm_bin, char* grmname, char* phenoname, uint32_t mpheno_col, char* phenoname_str, int32_t missing_pheno, uint32_t is_strict, double unrelated_herit_tol, double unrelated_herit_covg, double unrelated_herit_covr);
#endif

int32_t ibs_test_calc(pthread_t* threads, char* read_dists_fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t perm_ct, uintptr_t pheno_nm_ct, uintptr_t pheno_ctrl_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c);

int32_t groupdist_calc(pthread_t* threads, uint32_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t groupdist_iters, uint32_t groupdist_d, uint32_t pheno_nm_ct, uint32_t pheno_ctrl_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c);

int32_t calc_regress_pcs(char* evecname, uint32_t regress_pcs_modifier, uint32_t max_pcs, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t indiv_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uintptr_t max_person_id_len, uintptr_t* sex_nm, uintptr_t* sex_male, double* pheno_d, double missing_phenod, char* outname, char* outname_end);

int32_t calc_genome(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, double* set_allele_freqs, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, char* outname, char* outname_end, int32_t nonfounders, uint64_t calculation_type, uint32_t genome_modifier, uint32_t ppc_gap, uintptr_t* pheno_nm, uintptr_t* pheno_c, Pedigree_rel_info pri, uint32_t skip_write);

int32_t ld_prune(FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uint32_t* marker_pos, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* sex_male, uint32_t ld_window_size, uint32_t ld_window_kb, uint32_t ld_window_incr, double ld_last_param, char* outname, char* outname_end, uint64_t calculation_type);

int32_t rel_cutoff_batch(uint32_t load_grm_bin, char* grmname, char* outname, char* outname_end, double rel_cutoff, uint32_t rel_calc_type);

int32_t calc_rel(pthread_t* threads, uint32_t parallel_idx, uint32_t parallel_tot, uint64_t calculation_type, uint32_t rel_calc_type, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* person_ids, uintptr_t max_person_id_len, int32_t ibc_type, double rel_cutoff, double* set_allele_freqs, double** rel_ibc_ptr, Chrom_info* chrom_info_ptr);

int32_t calc_rel_f(pthread_t* threads, uint32_t parallel_idx, uint32_t parallel_tot, uint64_t calculation_type, uint32_t rel_calc_type, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, char* person_ids, uintptr_t max_person_id_len, int32_t ibc_type, float rel_cutoff, double* set_allele_freqs, Chrom_info* chrom_info_ptr);

int32_t calc_distance(pthread_t* threads, uint32_t parallel_idx, uint32_t parallel_tot, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, uint32_t dist_calc_type, uintptr_t* marker_exclude, uint32_t marker_ct, double* set_allele_freqs, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uintptr_t max_person_id_len, Chrom_info* chrom_info_ptr, uint32_t wt_needed, uint32_t* marker_weights_i, double exponent);

int32_t calc_cluster_neighbor(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, char* read_dists_fname, char* read_dists_id_fname, char* read_genome_fname, char* outname, char* outname_end, uint64_t calculation_type, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Cluster_info* cp, int32_t missing_pheno, uint32_t neighbor_n1, uint32_t neighbor_n2, uint32_t ppc_gap, uintptr_t* pheno_c, double* mds_plot_dmatrix_copy, uintptr_t* cluster_merge_prevented, double* cluster_sorted_ibs, unsigned char* wkspace_mark_precluster, unsigned char* wkspace_mark_postcluster);

#endif // __WDIST_CALC_H__
