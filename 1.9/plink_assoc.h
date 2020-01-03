#ifndef __PLINK_ASSOC_H__
#define __PLINK_ASSOC_H__

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


#define TESTMISS_PERM 1
#define TESTMISS_MPERM 2
#define TESTMISS_PERM_COUNT 4
#define TESTMISS_MIDP 8

#include "plink_set.h"

void aperm_init(Aperm_info* apip);

int32_t multcomp(char* outname, char* outname_end, uint32_t* marker_uidxs, uintptr_t chi_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, double* chi, double pfilter, double output_min_p, uint32_t mtest_adjust, uint32_t skip_gc, double adjust_lambda, uint32_t* tcnt, double* pvals);

int32_t model_assoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t model_modifier, uint32_t model_cell_ct, uint32_t model_mperm_val, double ci_size, double ci_zt, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* founder_info, uintptr_t* sex_male, uint32_t hh_exists, uint32_t ld_ignore_x, uint32_t perm_batch_size, Set_info* sip);

int32_t qassoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t model_modifier, uint32_t model_mperm_val, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* founder_info, uintptr_t* sex_male, uint32_t hh_exists, uint32_t ld_ignore_x, uint32_t perm_batch_size, Set_info* sip);

int32_t gxe_assoc(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uintptr_t sample_ct, uintptr_t* sample_exclude, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* gxe_covar_nm, uintptr_t* gxe_covar_c, uintptr_t* sex_male, uint32_t hh_or_mt_exists);

int32_t testmiss(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t testmiss_mperm_val, uint32_t testmiss_modifier, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male, uint32_t hh_exists);

int32_t cmh_assoc(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t cmh_mperm_val, uint32_t cmh_modifier, double ci_size, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male, uint32_t hh_or_mt_exists, Set_info* sip);

int32_t cmh2_assoc(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male, uint32_t hh_exists);

int32_t homog_assoc(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double output_min_p, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, Chrom_info* chrom_info_ptr, double* set_allele_freqs, uintptr_t unfiltered_sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, char* cluster_ids, uintptr_t max_cluster_id_len, uint32_t pheno_nm_ct, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* sex_male, uint32_t hh_or_mt_exists);

#endif // __PLINK_ASSOC_H__
