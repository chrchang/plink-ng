#ifndef __WDIST_DOSAGE_H__
#define __WDIST_DOSAGE_H__

int wdist_dosage(int calculation_type, char* genname, char* samplename, char* missing_code, int distance_3d, int distance_flat_missing, double exponent, int maf_succ, unsigned int thread_ct, int parallel_idx, unsigned int parallel_tot);

#endif
