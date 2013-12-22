#ifndef __PLINK_SET_H__
#define __PLINK_SET_H__

#include "plink_common.h"

#define SET_MAKE_FROM_RANGES 1
#define SET_COMPLEMENTS 2
#define SET_MAKE_COLLAPSE_GROUP 4
#define SET_C_PREFIX 8
#define SET_GENE_ALL 0x10
#define SET_WRITE_LIST 0x20
#define SET_WRITE_TABLE 0x40
#define SET_R2_WRITE 0x80

typedef struct {
  // command-line, allocated on heap and freed by main()
  char* fname;
  char* setnames_flattened;
  char* subset_fname;
  char* merged_set_name;
  char* genekeep_flattened;
  uint32_t modifier;
  uint32_t make_set_border;
  double set_r2;
  double set_p;
  uint32_t set_max;

  // main data structure, allocated on stack
  uintptr_t ct;
  char* names;
  uintptr_t max_name_len;

  // The simplest set representation is a bitfield of length
  // unfiltered_marker_ct.  But with, say, 30k genes and 35m variants, pretty
  // soon you're talking about some real memory consumption.
  // However, sets are normally used to represent genes, and genes are
  // localized.  So, if we track the location of the first and last set
  // entries, we can usually get away with a much smaller bitfield.  In fact,
  // often we won't need a bitfield at all; we can check if marker_idx is in a
  // single range or small list of them (assuming we don't query marker_idx
  // values which have been excluded elsewhere).

  // setdefs[r] oints to a 16-byte-aligned array of uint32_ts.  If the first
  // element is not 0xffffffffU, its contents are
  //   [0]: number of ranges in set #r
  //   [2k+1], [2k+2]: start and end of range k (half-open); sorted
  // Otherwise, set #r is represented as a bitfield starting from [4], and
  //   [1]: offset of first bit (always divisible by 128)
  //   [2]: number of bits (divisible by 128 unless very last bit included)
  //   [3]: 1 if all out-of-bounds bits are set, 0 otherwise (other flags may
  //        be added later)
  uint32_t** setdefs;
} Set_info;

void set_init(Set_info* sip);

void set_cleanup(Set_info* sip);

int32_t define_sets(Set_info* sip, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, uintptr_t* marker_exclude_ct_ptr, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr);

int32_t write_set(Set_info* sip, char* outname, char* outname_end, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, uint32_t zero_extra_chroms, Chrom_info* chrom_info_ptr);

void unpack_set_unfiltered(uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* setdef, uintptr_t* new_exclude);

uint32_t extract_set_union_unfiltered(Set_info* sip, uintptr_t* set_incl, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t** union_marker_exclude_ptr, uintptr_t* union_marker_ct_ptr);

uint32_t setdefs_compress(Set_info* sip, uintptr_t* set_incl, uintptr_t set_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, uintptr_t* marker_exclude, uintptr_t marker_ct, uint32_t*** new_setdefs_ptr);

#endif // __PLINK_SET_H__
