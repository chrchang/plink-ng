#include "wdist_common.h"

int32_t wdist_cnv(char* outname, char* outname_end, char* cnvname, char* mapname, char* famname, char* phenoname, uint32_t cnv_calc_type, uint32_t cnv_kb, uint32_t cnv_max_kb, uint32_t cnv_score, uint32_t cnv_max_score, uint32_t cnv_sites, uint32_t cnv_max_sites, uint32_t cnv_intex_type, char* cnv_intex_fname, char* cnv_subset_fname, uint32_t cnv_overlap_type, double cnv_overlap_val, uint32_t cnv_freq_type, uint32_t cnv_freq_val, uint32_t mperm_val, uint32_t cnv_test_window, uint32_t segment_modifier, char* segment_spanning_fname, uint32_t cnv_indiv_mperms, uint32_t cnv_test_mperms, uint32_t cnv_test_region_mperms, uint32_t cnv_enrichment_test_mperms, Chrom_info* chrom_info_ptr) {
  unsigned char* wkspace_mark = wkspace_base;
  FILE* cnvfile = NULL;
  FILE* famfile = NULL;
  FILE* mapfile = NULL;
  FILE* outfile = NULL;
  int32_t retval = 0;
  while (0) {
  wdist_cnv_ret_NOMEM:
    retval = RET_NOMEM;
    break;
  wdist_cnv_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  wdist_cnv_ret_READ_FAIL:
    retval = RET_READ_FAIL;
    break;
  wdist_cnv_ret_WRITE_FAIL:
    retval = RET_WRITE_FAIL;
    break;
  }
  fclose_cond(cnvfile);
  fclose_cond(famfile);
  fclose_cond(mapfile);
  fclose_cond(outfile);
  wkspace_reset(wkspace_mark);
  return 0;
}

int32_t wdist_gvar(char* outname, char* outname_end, char* gvarname, char* mapname, char* famname) {
  logprint("Error: Common CNP analysis not yet supported.\n");
  return RET_CALC_NOT_YET_SUPPORTED;
}
