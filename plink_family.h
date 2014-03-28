#ifndef __PLINK_FAMILY_H__
#define __PLINK_FAMILY_H__

#include "plink_common.h"

#define TDT_EXACT 1
#define TDT_MIDP 2
#define TDT_POO 4
#define TDT_MPERM 8
#define TDT_PARENPERM1 0x10
#define TDT_PARENPERM2 0x20
#define TDT_POOPERM_PAT 0x40
#define TDT_POOPERM_MAT 0x80
#define TDT_SET_TEST 0x100

extern const uint32_t mendel_error_table[];
extern const uint32_t mendel_error_table_x[];

typedef struct {
  uint32_t modifier;
  double max_trio_error;
  double max_var_error;
  double exclude_one_ratio;
} Mendel_info;

void family_init(Mendel_info* me_ip);

int32_t get_trios_and_families(uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t indiv_ct, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, char** fids_ptr, uintptr_t* max_fid_len_ptr, char** iids_ptr, uintptr_t* max_iid_len_ptr, uint64_t** family_list_ptr, uint32_t* family_ct_ptr, uint64_t** trio_list_ptr, uintptr_t* trio_ct_ptr, uint32_t** trio_lookup_ptr, uint32_t include_duos, uint32_t toposort);

int32_t mendel_error_scan(Mendel_info* me_ip, FILE* bedfile, uintptr_t bed_offset, char* outname, char* outname_end, uint32_t plink_maxfid, uint32_t plink_maxiid, uint32_t plink_maxsnp, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uintptr_t* marker_reverse, char* marker_ids, uintptr_t max_marker_id_len, char** marker_allele_ptrs, uintptr_t max_marker_allele_len, uintptr_t unfiltered_indiv_ct, uintptr_t* indiv_exclude, uintptr_t* indiv_exclude_ct_ptr, uintptr_t* founder_info, uintptr_t* sex_nm, uintptr_t* sex_male, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uint32_t hh_exists, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr, uint32_t calc_mendel);

typedef struct {
  char* family_ids;
  uintptr_t max_family_id_len; // includes trailing null
  uint32_t* family_sizes;

  uintptr_t* family_rel_space_offsets; // offset for rel_space lookup
  uint32_t* family_founder_cts;
  // direct indiv uidx -> family idx lookup, to reduce number of bsearches
  uint32_t* family_idxs;

  // truncated triangular arrays of pedigree coefficient of relationship
  double* rel_space;

  // direct indiv idx -> rel_space idx lookup
  uint32_t* family_rel_nf_idxs;

  // following three variables are technically unnecessary for --genome, but we
  // get them for "free" in the process of calculating everything else, and
  // they'll be nice to use if we ever need to iterate by family in the future.
  uint32_t family_id_ct;
  // list of idxs of all individuals in first family, then second family, etc.
  uint32_t* family_info_space;
  uint32_t* family_info_offsets; // offset in family_info_space
} Pedigree_rel_info;

int32_t populate_pedigree_rel_info(Pedigree_rel_info* pri_ptr, uintptr_t unfiltered_indiv_ct, char* person_ids, uintptr_t max_person_id_len, char* paternal_ids, uintptr_t max_paternal_id_len, char* maternal_ids, uintptr_t max_maternal_id_len, uintptr_t* founder_info);

#endif // __PLINK_FAMILY_H__
