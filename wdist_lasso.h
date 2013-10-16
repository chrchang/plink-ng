#ifndef __WDIST_LASSO_H__
#define __WDIST_LASSO_H__

#include "wdist_common.h"

int32_t lasso(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char* marker_alleles, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_indiv_ct, uintptr_t indiv_ct, double lasso_h2, uint32_t report_zeroes, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* sex_male, uint32_t hh_exists);

#endif // __WDIST_LASSO_H__
