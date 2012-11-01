#include "wdist_common.h"
#include "wdist_dosage.h"

int wdist_dosage(int calculation_type, char* genname, char* samplename, int distance_3d) {
  FILE* genfile = NULL;
  FILE* samplefile = NULL;
  int retval = 0;
  if (fopen_checked(&genfile, genname, "r")) {
    return RET_OPEN_FAIL;
  }
  if (fopen_checked(&samplefile, samplename, "r")) {
    goto wdist_dosage_ret_OPEN_FAIL;
  }
  printf("wdist_dosage() incomplete.\n");
  while (0) {
  wdist_dosage_ret_OPEN_FAIL:
    retval = RET_OPEN_FAIL;
  }
  fclose_cond(genfile);
  fclose_cond(samplefile);
  return retval;
}
