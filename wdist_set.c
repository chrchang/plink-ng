#include "wdist_set.h"

void set_init(Set_info* sip) {
  sip->fname = NULL;
  sip->subset_fname = NULL;
  sip->merged_set_name = NULL;
  sip->genekeep_flattened = NULL;
}

void set_cleanup(Set_info* sip) {
  free_cond(sip->fname);
  free_cond(sip->subset_fname);
  free_cond(sip->merged_set_name);
  free_cond(sip->genekeep_flattened);
}

int32_t define_sets(Set_info* sip, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, char* marker_ids, uintptr_t max_marker_id_len) {
  logprint("Error: --set and --make-set are currently under development.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
}
