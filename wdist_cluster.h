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
  char* match_type_fname;
  char* qmatch_fname;
  char* qt_fname;
  uint32_t modifier;
  double ppc;
  uint32_t max_size;
  uint32_t max_cases;
  uint32_t max_controls;
  uint32_t min_ct;
  uint32_t mds_dim_ct;
  uint32_t neighbor_n1; // yeah, American spelling, deal with it
  uint32_t neighbor_n2;
  double max_missing_discordance;
} Cluster_info;

void cluster_init(Cluster_info* cluster_ptr);

int32_t load_clusters(char* fname, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uint32_t mwithin_col, uint32_t keep_na, uintptr_t* cluster_ct_ptr, uint32_t** cluster_map_ptr, uint32_t** cluster_starts_ptr, char** cluster_ids_ptr, uintptr_t* max_cluster_id_len_ptr);

int32_t write_clusters(char* outname, char* outname_end, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uint32_t omit_unassigned, uintptr_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len);

#endif // __WDIST_CLUSTER_H__
