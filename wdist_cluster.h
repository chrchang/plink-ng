#ifndef __WDIST_CLUSTER_H__
#define __WDIST_CLUSTER_H__

#include "wdist_common.h"

#define CLUSTER_CC 1
#define CLUSTER_GROUP_AVG 2
#define CLUSTER_MISSING 4
#define CLUSTER_ONLY2 8
#define CLUSTER_MDS 0x10

typedef struct {
  char* fname;
  char* match_fname;
  char* match_missing_str;
  char* match_type_fname;
  char* qmatch_fname;
  char* qmatch_missing_str;
  char* qt_fname;
  uint32_t modifier;
  double ppc;
  uint32_t max_size;
  uint32_t max_cases;
  uint32_t max_controls;
  uint32_t min_ct;
  uint32_t mds_dim_ct;
  double min_ibm;
} Cluster_info;

void cluster_init(Cluster_info* cluster_ptr);

int32_t load_clusters(char* fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uint32_t mwithin_col, uint32_t keep_na, uintptr_t* cluster_ct_ptr, uint32_t** cluster_map_ptr, uint32_t** cluster_starts_ptr, char** cluster_ids_ptr, uintptr_t* max_cluster_id_len_ptr);

void fill_unfiltered_indiv_to_cluster(uintptr_t unfiltered_indiv_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* indiv_to_cluster);

int32_t fill_indiv_to_cluster(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* indiv_to_cluster);

int32_t write_clusters(char* outname, char* outname_end, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uint32_t omit_unassigned, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len);

int32_t cluster_include_and_reindex(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_include, uint32_t remove_size1, uintptr_t* pheno_c, uintptr_t indiv_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t* new_cluster_ct_ptr, uint32_t** new_cluster_map_ptr, uint32_t** new_cluster_starts_ptr, uint32_t** cluster_case_cts_ptr, uintptr_t** cluster_cc_perm_preimage_ptr);

int32_t cluster_alloc_and_populate_magic_nums(uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t** tot_quotients_ptr, uint64_t** totq_magics_ptr, uint32_t** totq_preshifts_ptr, uint32_t** totq_postshifts_ptr, uint32_t** totq_incrs_ptr);

int32_t read_genome(char* read_genome_fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t* cluster_merge_prevented, double* cluster_sdistances, double* mds_plot_dmatrix_copy, uint32_t neighbor_n2, double* neighbor_quantiles, uint32_t* neighbor_qindices, uint32_t* ppc_fail_counts, double min_ppc, uintptr_t cluster_ct, uint32_t* cluster_starts, uint32_t* indiv_to_cluster);

#endif // __WDIST_CLUSTER_H__
