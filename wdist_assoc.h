#ifndef __WDIST_ASSOC_H__
#define __WDIST_ASSOC_H__

#include "wdist_common.h"

int32_t model_assoc(pthread_t* threads, FILE* bedfile, int32_t bed_offset, char* outname, char* outname_end, uint64_t calculation_type, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double ci_zt, double pfilter, uint32_t mtest_adjust, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uint32_t aperm_min, uint32_t aperm_max, double aperm_alpha, double aperm_beta, double aperm_init_interval, double aperm_interval_slope, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t* sex_nm, uintptr_t* sex_male);

#endif // __WDIST_ASSOC_H__
