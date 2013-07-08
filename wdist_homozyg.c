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

int32_t calc_homozyg(Homozyg_info* hp, FILE* bedfile, uintptr_t bed_offset, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, uint32_t* marker_pos, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, char* person_ids, uint32_t plink_maxfid, uint32_t plink_maxiid, uintptr_t max_person_id_len, char* outname, char* outname_end, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d) {
  int32_t retval = 0;
  logprint("Error: --homozyg is currently under development.\n");
  retval = RET_CALC_NOT_YET_SUPPORTED;
  goto calc_homozyg_ret_1;
  while (0) {
  }
 calc_homozyg_ret_1:
  return retval;
}
