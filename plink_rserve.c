// Rconnection.cc is part of a continuously updated package and uses
// C++-specific language constructs, so we omit this module when compiling with
// gcc instead of g++.
#if defined __cplusplus && !defined _WIN32

#include "plink_common.h"

#define MAIN
#define SOCK_ERRORS
#include "sisocks.h"
#include "Rconnection.h"

int32_t rserve_call(char* rplugin_fname, uint32_t rplugin_port, FILE* bedfile, uintptr_t bed_offset, uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uint32_t* marker_pos, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t* pheno_nm, uint32_t pheno_nm_ct, uintptr_t* pheno_c, double* pheno_d, char* outname, char* outname_end) {
  FILE* infile = NULL;
  int32_t retval = 0;
  if (fopen_checked(&infile, rplugin_fname, "r")) {
    goto rserve_call_ret_OPEN_FAIL;
  }
  while (0) {
  rserve_call_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
    break;
  }
  fclose_cond(infile);
  return retval;
}

#endif // __cplusplus, !_WIN32
