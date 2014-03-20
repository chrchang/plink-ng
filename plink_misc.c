#include "plink_misc.h"

void misc_init(Score_info* sc_ip) {
  sc_ip->fname = NULL;
  sc_ip->range_fname = NULL;
  sc_ip->data_fname = NULL;
  sc_ip->modifier = 0;
  sc_ip->varid_col = 1;
  sc_ip->allele_col = 2;
  sc_ip->effect_col = 3;
  sc_ip->data_varid_col = 1;
  sc_ip->data_col = 2;
}

void misc_cleanup(Score_info* sc_ip) {
  free_cond(sc_ip->fname);
  free_cond(sc_ip->range_fname);
  free_cond(sc_ip->data_fname);
}

int32_t score_report(FILE* bedfile, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, char* person_ids, uintptr_t max_person_id_len, uintptr_t* sex_male, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uint32_t hh_exists, Chrom_info* chrom_info_ptr, Score_info* sc_ip) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* outfile = NULL;
  int32_t retval = 0;
  if (sc_ip->range_fname) {
    logprint("Error: --q-score-range is currently under development.\n");
    return RET_CALC_NOT_YET_SUPPORTED;
  }
  logprint("Error: --score is currently under development.\n");
  retval = RET_CALC_NOT_YET_SUPPORTED;
  while (0) {
  score_report_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  }
  wkspace_reset(wkspace_mark);
  fclose_cond(outfile);
  return retval;
}
