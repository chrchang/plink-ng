#ifndef __WDIST_HOMOZYG_H__
#define __WDIST_HOMOZYG_H__

#include "wdist_common.h"

typedef struct {
  uint32_t modifier;
  uint32_t min_snp;
  uint32_t min_bases;
  double max_bases_per_snp;
  uint32_t max_gap;
  uint32_t max_hets;
  uint32_t window_size;
  uint32_t window_max_hets;
  uint32_t window_max_missing;
  double hit_threshold;
  double overlap_min;
  uint32_t pool_size_min;
} Homozyg_info;

#define HOMOZYG_GROUP 1
#define HOMOZYG_GROUP_VERBOSE 2
#define HOMOZYG_CONSENSUS_MATCH 4
#define HOMOZYG_INCLUDE_MISSING 8

void homozyg_init(Homozyg_info* homozyg_ptr);

#endif // __WDIST_HOMOZYG_H__
