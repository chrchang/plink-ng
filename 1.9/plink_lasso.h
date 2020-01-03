#ifndef __PLINK_LASSO_H__
#define __PLINK_LASSO_H__

// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


int32_t lasso_lambda(const uintptr_t* marker_exclude, const uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, const uintptr_t* sample_exclude, uintptr_t* sex_male, uintptr_t* pheno_nm, const uintptr_t* covar_nm, uintptr_t bed_offset, uintptr_t unfiltered_marker_ct, uintptr_t marker_ct, uintptr_t unfiltered_sample_ct, uint32_t sample_ct, uintptr_t pheno_nm_ct, uint32_t hh_or_mt_exists, uint32_t lasso_lambda_iters, double lasso_h2, FILE* bedfile, char* outname, char* outname_end, double* lasso_minlambda_ptr);

int32_t lasso(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uint32_t sample_ct, uintptr_t pheno_nm_ct, double lasso_h2, double lasso_minlambda, Range_list* select_covars_range_list_ptr, uint64_t misc_flags, const uintptr_t* sample_exclude, uintptr_t* pheno_nm, uintptr_t* pheno_c, double* pheno_d, uintptr_t covar_ct, char* covar_names, uintptr_t max_covar_name_len, uintptr_t* covar_nm, double* covar_d, uintptr_t* sex_male, uint32_t hh_or_mt_exists);

#endif // __PLINK_LASSO_H__
