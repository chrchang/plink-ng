#include "wdist_common.h"

int32_t model_assoc(pthread_t* threads, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, int32_t calculation_type, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double pfilter, uint32_t mtest_adjust, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude) {
  unsigned char* wkspace_mark = wkspace_base;
  int32_t retval = 0;
  FILE* outfile = NULL;

  logprint("Error: --assoc and --model are currently under development.\n");
  retval = RET_CALC_NOT_YET_SUPPORTED;
  goto model_assoc_ret_1;

  while (0) {
  }
 model_assoc_ret_1:
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}
