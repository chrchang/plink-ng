#ifndef __WDIST_DOSAGE_H__
#define __WDIST_DOSAGE_H__

#include "wdist_common.h"

#ifndef STABLE_BUILD
int32_t wdist_dosage(uint64_t calculation_type, uint32_t dist_calc_type, char* genname, char* samplename, char* outname, char* outname_end, char* missing_code, double exponent, uint32_t maf_succ, uintptr_t regress_iters, uint32_t regress_d, uint32_t thread_ct, uint32_t parallel_idx, uint32_t parallel_tot);
#endif

#endif
