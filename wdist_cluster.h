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

#endif // __WDIST_CLUSTER_H__
