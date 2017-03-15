#ifndef __PLINK_CLUSTER_H__
#define __PLINK_CLUSTER_H__

#define CLUSTER_CC 1
#define CLUSTER_GROUP_AVG 2
#define CLUSTER_MISSING 4
#define CLUSTER_ONLY2 8
#define CLUSTER_OLD_TIEBREAKS 0x10
#define CLUSTER_MDS 0x20
#define CLUSTER_MDS_EIGENDECOMP 0x40
#define CLUSTER_MDS_EIGVALS 0x80
#define CLUSTER_CMH_BD 0x100
#define CLUSTER_CMH_PERM 0x200
#define CLUSTER_CMH_PERM_BD 0x400
#define CLUSTER_CMH_MPERM 0x800
#define CLUSTER_CMH_PERM_COUNT 0x1000
#define CLUSTER_CMH_SET_TEST 0x2000
#define CLUSTER_CMH2 0x4000

typedef struct {
  char* fname;
  char* match_fname;
  char* match_missing_str;
  char* match_type_fname;
  char* qmatch_fname;
  char* qmatch_missing_str;
  char* qt_fname;
  char* keep_fname;
  char* remove_fname;
  char* keep_flattened;
  char* remove_flattened;
  char* zerofname;
  uint32_t modifier;
  double ppc;
  uint32_t max_size;
  uint32_t max_cases;
  uint32_t max_ctrls;
  uint32_t min_ct;
  uint32_t mds_dim_ct;
  uint32_t cmh_mperm_val;
  double min_ibm;
} Cluster_info;

void cluster_init(Cluster_info* cluster_ptr);

void cluster_cleanup(Cluster_info* cluster_ptr);

int32_t load_clusters(char* fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, char* sample_ids, uintptr_t max_sample_id_len, uint32_t mwithin_col, uint32_t keep_na, uintptr_t* cluster_ct_ptr, uint32_t** cluster_map_ptr, uint32_t** cluster_starts_ptr, char** cluster_ids_ptr, uintptr_t* max_cluster_id_len_ptr, char* keep_fname, char* keep_flattened, char* remove_fname, char* remove_flattened, uint32_t allow_no_samples);

void fill_unfiltered_sample_to_cluster(uintptr_t unfiltered_sample_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* sample_to_cluster);

int32_t fill_sample_to_cluster(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* sample_to_cluster, uint32_t* late_clidx_to_sample_uidx);

int32_t write_clusters(char* outname, char* outname_end, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uint32_t omit_unassigned, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len);

int32_t extract_clusters(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, char* cluster_names_flattened, char* clusters_fname, uintptr_t** new_sample_exclude_ptr, uintptr_t* new_sample_ct_ptr);

int32_t cluster_include_and_reindex(uintptr_t unfiltered_sample_ct, uintptr_t* sample_include, uint32_t remove_size1, uintptr_t* pheno_c, uintptr_t sample_ct, uint32_t is_perm1, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* new_cluster_ct_ptr, uint32_t** new_cluster_map_ptr, uint32_t** new_cluster_starts_ptr, uint32_t** cluster_case_cts_ptr, uintptr_t** cluster_cc_perm_preimage_ptr);

int32_t cluster_alloc_and_populate_magic_nums(uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t** tot_quotients_ptr, uint64_t** totq_magics_ptr, uint32_t** totq_preshifts_ptr, uint32_t** totq_postshifts_ptr, uint32_t** totq_incrs_ptr);

int32_t read_dists(char* dist_fname, char* id_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* sample_to_cluster, uint32_t for_cluster_flag, uint32_t is_max_dist, double* dists, uint32_t neighbor_n2, double* neighbor_quantiles, uint32_t* neighbor_qindices);

void update_neighbor(uintptr_t sample_ct, uint32_t neighbor_n2, uintptr_t sample_idx1, uintptr_t sample_idx2, double cur_ibs, double* neighbor_quantiles, uint32_t* neighbor_qindices);

int32_t read_genome(char* read_genome_fname, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* cluster_merge_prevented, double* cluster_sorted_ibs, uint32_t neighbor_n2, double* neighbor_quantiles, uint32_t* neighbor_qindices, uint32_t* ppc_fail_counts, double min_ppc, uint32_t is_max_dist, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* sample_to_cluster);

int32_t cluster_enforce_match(Cluster_info* cp, int32_t missing_pheno, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* sample_to_cluster, uintptr_t* merge_prevented);

uint32_t cluster_main(uintptr_t cluster_ct, uintptr_t* merge_prevented, uintptr_t list_size, uint32_t* sorted_ibs_indices, uint32_t* cluster_index, uint32_t* cur_cluster_sizes, uint32_t sample_ct, uint32_t* cur_cluster_case_cts, uint32_t case_ct, uint32_t ctrl_ct, uint32_t* cur_cluster_remap, Cluster_info* cp, uintptr_t* ibs_ties, uint32_t* merge_sequence);

uint32_t cluster_group_avg_main(uint32_t cluster_ct, uintptr_t* merge_prevented, uint32_t heap_size, double* heap_vals, uint32_t* val_to_cindices, uint32_t* cluster_index, uint32_t* cur_cluster_sizes, uint32_t sample_ct, uint32_t* cur_cluster_case_cts, uint32_t case_ct, uint32_t ctrl_ct, uint32_t* cur_cluster_remap, Cluster_info* cp, uint32_t* merge_sequence);

int32_t write_cluster_solution(char* outname, char* outname_end, uint32_t* orig_sample_to_cluster, uintptr_t sample_ct, uint32_t* orig_cluster_map, uint32_t* orig_cluster_starts, uint32_t* late_clidx_to_sample_uidx, uint32_t orig_within_ct, uint32_t orig_cluster_ct, char* sample_ids, uintptr_t max_sample_id_len, uintptr_t* pheno_c, uint32_t* sample_idx_to_uidx, Cluster_info* cp, uint32_t* cluster_remap, uint32_t* clidx_table_space, uint32_t merge_ct, uint32_t* merge_sequence);

#ifndef NOLAPACK
int32_t mds_plot(char* outname, char* outname_end, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t* sample_idx_to_uidx, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uint32_t cur_cluster_ct, uint32_t merge_ct, uint32_t* orig_sample_to_cluster, uint32_t* cur_cluster_remap, uint32_t dim_ct, uint32_t is_mds_cluster, uint32_t dump_eigvals, double* dists);

int32_t mds_plot_eigendecomp(char* outname, char* outname_end, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t* sample_idx_to_uidx, char* sample_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_sample_id_len, uint32_t cur_cluster_ct, uint32_t merge_ct, uint32_t* orig_sample_to_cluster, uint32_t* cur_cluster_remap, uint32_t dim_ct, uint32_t is_mds_cluster, uint32_t dump_eigvals, double* dists);
#endif

#endif // __PLINK_CLUSTER_H__
