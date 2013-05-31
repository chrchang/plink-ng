#include "wdist_homozyg.h"

void homozyg_init(Homozyg_info* homozyg_ptr) {
  homozyg_ptr->modifier = 0;
  homozyg_ptr->min_snp = 0;
  homozyg_ptr->min_bases = 0;
  homozyg_ptr->max_bases_per_snp = INFINITY;
  homozyg_ptr->max_gap = 0x7fffffff;
  homozyg_ptr->max_hets = 0x7fffffff;
  homozyg_ptr->window_size = 50;
  homozyg_ptr->window_max_hets = 1;
  homozyg_ptr->window_max_missing = 5;
  homozyg_ptr->hit_threshold = 0.05;
  homozyg_ptr->overlap_min = 0.95;
  homozyg_ptr->pool_size_min = 2;
}
