#ifndef __PLINK_CALC_H__
#define __PLINK_CALC_H__

#include "plink_cluster.h"
#include "plink_family.h"

#define REL_CALC_COV 1
#define REL_CALC_SQ 2
#define REL_CALC_SQ0 4
#define REL_CALC_TRI 6
#define REL_CALC_SHAPEMASK 6
#define REL_CALC_GZ 8
#define REL_CALC_BIN 0x10
#define REL_CALC_BIN4 0x20
#define REL_CALC_GRM 0x40
#define REL_CALC_GRM_BIN 0x80
#define REL_CALC_MASK 0xff
#define REL_UNRELATED_HERITABILITY_STRICT 0x100
#define REL_PCA_HEADER 0x200
#define REL_PCA_TABS 0x400
#define REL_PCA_VAR_WTS 0x800

typedef struct {
  uint32_t modifier;
  uint32_t regress_rel_d;
  uintptr_t regress_rel_iters;
  double unrelated_herit_tol;
  double unrelated_herit_covg;
  double unrelated_herit_covr;
  double cutoff;
  char* pca_cluster_names_flattened;
  char* pca_clusters_fname;
  int32_t ibc_type;
  uint32_t pc_ct;
} Rel_info;

extern double* g_rel_dists;
extern uint32_t* g_sample_missing_unwt;
extern uint32_t* g_missing_dbl_excluded;
extern double* g_dists;

void rel_init(Rel_info* relip);

void rel_cleanup(Rel_info* relip);

int32_t regress_rel_main(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, Rel_info* relip, pthread_t* threads, double* pheno_d);

#ifndef NOLAPACK
int32_t calc_unrelated_herit(uint64_t calculation_type, Rel_info* relip, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, double* pheno_d, double* rel_ibc);

int32_t unrelated_herit_batch(uint32_t load_grm_bin, char* grmname, char* phenoname, uint32_t mpheno_col, char* phenoname_str, int32_t missing_pheno, Rel_info* relip);
#endif

int32_t ibs_test_calc(pthread_t* threads, char* read_dists_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t perm_ct, uintptr_t pheno_nm_ct, uintptr_t pheno_ctrl_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c);

int32_t groupdist_calc(pthread_t* threads, uint32_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t groupdist_iters, uint32_t groupdist_d, uint32_t pheno_nm_ct, uint32_t pheno_ctrl_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c);

int32_t calc_genome(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, double* set_allele_freqs, uint32_t* nchrobs, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info, uint32_t parallel_idx, uint32_t parallel_tot, char* outname, char* outname_end, int32_t nonfounders, uint64_t calculation_type, uint32_t genome_modifier, uint32_t ppc_gap, double min_pi_hat, double max_pi_hat, uintptr_t* pheno_nm, uintptr_t* pheno_c, Pedigree_rel_info pri, uint32_t skip_write);

int32_t rel_cutoff_batch(uint32_t load_grm_bin, char* grmname, char* outname, char* outname_end, Rel_info* relip);

int32_t calc_rel(pthread_t* threads, uint32_t parallel_idx, uint32_t parallel_tot, uint64_t calculation_type, Rel_info* relip, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, char* distance_wts_fname, uint32_t distance_wts_noheader, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t* marker_reverse, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, char* sample_ids, uintptr_t max_sample_id_len, double* set_allele_freqs, double** rel_ibc_ptr, Chrom_info* chrom_info_ptr);

#ifndef NOLAPACK
int32_t calc_pca(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, Rel_info* relip, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* pca_sample_exclude, uintptr_t pca_sample_ct, char* sample_ids, uintptr_t max_sample_id_len, double* set_allele_freqs, Chrom_info* chrom_info_ptr, double* rel_ibc);
#endif

int32_t calc_distance(pthread_t* threads, uint32_t parallel_idx, uint32_t parallel_tot, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, char* read_dists_fname, char* distance_wts_fname, double distance_exp, uint64_t calculation_type, uint32_t dist_calc_type, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, Chrom_info* chrom_info_ptr);

int32_t calc_cluster_neighbor(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, char* read_dists_fname, char* read_dists_id_fname, char* read_genome_fname, char* outname, char* outname_end, uint64_t calculation_type, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Cluster_info* cp, int32_t missing_pheno, uint32_t neighbor_n1, uint32_t neighbor_n2, uint32_t ppc_gap, uintptr_t* pheno_c, double* mds_plot_dmatrix_copy, uintptr_t* cluster_merge_prevented, double* cluster_sorted_ibs, unsigned char* bigstack_mark_precluster, unsigned char* bigstack_mark_postcluster);

int32_t regress_distance(pthread_t* threads, uint64_t calculation_type, double* pheno_d, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t thread_ct, uintptr_t regress_iters, uint32_t regress_d);

#endif // __PLINK_CALC_H__
