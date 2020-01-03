#ifndef __PLINK_FAMILY_H__
#define __PLINK_FAMILY_H__

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


#include "plink_set.h"

#define TDT_EXACT 1
#define TDT_MIDP 2
#define TDT_POO 4
#define TDT_PERM 8
#define TDT_MPERM 0x10
#define TDT_PERM_COUNT 0x20
#define TDT_PARENPERM1 0x40
#define TDT_PARENPERM2 0x80
#define TDT_POOPERM_PAT 0x100
#define TDT_POOPERM_MAT 0x200
#define TDT_SET_TEST 0x400

#define DFAM_NO_UNRELATEDS 1
#define DFAM_PERM 2
#define DFAM_MPERM 4
#define DFAM_PERM_COUNT 8
#define DFAM_SET_TEST 0x10

#define QFAM_WITHIN1 1
#define QFAM_WITHIN2 2
#define QFAM_BETWEEN 4
#define QFAM_TOTAL 8
#define QFAM_TEST 15
#define QFAM_PERM 0x10
#define QFAM_MPERM 0x20
#define QFAM_PERM_COUNT 0x40
#define QFAM_EMP_SE 0x80

extern const uint32_t mendel_error_table[];
extern const uint32_t mendel_error_table_x[];

typedef struct {
  double mendel_max_trio_error;
  double mendel_max_var_error;
  double mendel_exclude_one_ratio;
  uint32_t mendel_modifier;
  uint32_t tdt_modifier;
  uint32_t tdt_mperm_val;
  uint32_t dfam_modifier;
  uint32_t dfam_mperm_val;
  uint32_t qfam_modifier;
  uint32_t qfam_mperm_val;
  uint32_t tucc_bed;
} Family_info;

void family_init(Family_info* fam_ip);

int32_t get_trios_and_families(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, char** fids_ptr, uintptr_t* max_fid_len_ptr, char** iids_ptr, uintptr_t* max_iid_len_ptr, uint64_t** family_list_ptr, uint32_t* family_ct_ptr, uint64_t** trio_list_ptr, uintptr_t* trio_ct_ptr, uint32_t** trio_lookup_ptr, uint32_t include_duos, uint32_t toposort);

uint32_t erase_mendel_errors(uintptr_t unfiltered_sample_ct, uintptr_t* loadbuf, uintptr_t* workbuf, uintptr_t* sex_male, uint32_t* trio_lookup, uint32_t trio_ct, uint32_t is_x, uint32_t multigen);

int32_t mendel_error_scan(Family_info* fam_ip, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uint32_t allow_no_variants, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t* sample_exclude_ct_ptr, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t hh_exists, Chrom_info* chrom_info_ptr, uint32_t calc_mendel);

typedef struct {
  char* family_ids;
  uintptr_t max_family_id_len; // includes trailing null
  uint32_t* family_sizes;

  uintptr_t* family_rel_space_offsets; // offset for rel_space lookup
  uint32_t* family_founder_cts;
  // direct sample uidx -> family idx lookup, to reduce number of bsearches
  uint32_t* family_idxs;

  // truncated triangular arrays of pedigree coefficient of relationship
  double* rel_space;

  // direct sample idx -> rel_space idx lookup
  uint32_t* family_rel_nf_idxs;

  // following three variables are technically unnecessary for --genome, but we
  // get them for "free" in the process of calculating everything else, and
  // they'll be nice to use if we ever need to iterate by family in the future.
  uint32_t family_id_ct;
  // list of idxs of all samples in first family, then second family, etc.
  uint32_t* family_info_space;
  uint32_t* family_info_offsets; // offset in family_info_space
} Pedigree_rel_info;

int32_t populate_pedigree_rel_info(Pedigree_rel_info* pri_ptr, uintptr_t unfiltered_sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info);

int32_t tdt(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double ci_size, double ci_zt, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, Aperm_info* apip, uint32_t mperm_save, uintptr_t* pheno_nm, uintptr_t* pheno_c, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, Family_info* fam_ip);

int32_t dfam(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, double pfilter, double output_min_p, uint32_t mtest_adjust, double adjust_lambda, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uint32_t cluster_ct, uint32_t* cluster_map, uint32_t* cluster_starts, Aperm_info* apip, uint32_t mperm_save, uintptr_t* pheno_c, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uint32_t within_cmdflag, uint32_t perm_batch_size, Family_info* fam_ip, Set_info* sip);

int32_t qfam(pthread_t* threads, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_len, uint32_t plink_maxsnp, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, Aperm_info* apip, uintptr_t* pheno_nm, double* pheno_d, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uint32_t perm_batch_size, Family_info* fam_ip);

int32_t make_pseudocontrols(FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t marker_ct, char* marker_ids, uintptr_t max_marker_id_blen, double* marker_cms, uint32_t* marker_pos, char** marker_allele_ptrs, uintptr_t max_marker_allele_blen, uintptr_t* marker_reverse, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uintptr_t sample_ct, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* sample_ids, uintptr_t max_sample_id_blen, char* paternal_ids, uintptr_t max_paternal_id_blen, char* maternal_ids, uintptr_t max_maternal_id_blen, Chrom_info* chrom_info_ptr, Family_info* fam_ip);

#endif // __PLINK_FAMILY_H__
