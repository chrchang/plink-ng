#ifndef __WDIST_DOSAGE_H__
#define __WDIST_DOSAGE_H__

int32_t wdist_dosage(uint64_t calculation_type, int32_t dist_calc_type, char* genname, char* samplename, char* outname, char* outname_end, char* missing_code, int32_t distance_3d, int32_t distance_flat_missing, double exponent, int32_t maf_succ, uintptr_t regress_iters, uint32_t regress_d, uint32_t thread_ct, int32_t parallel_idx, uint32_t parallel_tot);

#endif
