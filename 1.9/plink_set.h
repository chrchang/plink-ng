#ifndef __PLINK_SET_H__
#define __PLINK_SET_H__

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


#define SET_MAKE_FROM_RANGES 1
#define SET_COMPLEMENTS 2
#define SET_MAKE_COLLAPSE_GROUP 4
#define SET_C_PREFIX 8
#define SET_GENE_ALL 0x10
#define SET_WRITE_LIST 0x20
#define SET_WRITE_TABLE 0x40
#define SET_R2_WRITE 0x80

typedef struct make_set_range_struct {
  struct make_set_range_struct* next;
  uint32_t uidx_start;
  uint32_t uidx_end;
} Make_set_range;

typedef struct {
  // command-line, allocated on heap and freed by main()
  char* fname;
  char* setnames_flattened;
  char* subset_fname;
  char* merged_set_name;
  char* genekeep_flattened;
  uint32_t modifier;

  // main data structure, allocated on stack
  uint32_t make_set_border;
  double set_r2;
  double set_p;
  double set_test_lambda;
  uint32_t set_max;

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

  // setdefs[r] points to a 16-byte-aligned array of uint32_ts.  If the first
  // element is not 0xffffffffU, its contents are
  //   [0]: number of ranges in set #r
  //   [2k+1], [2k+2]: start and end of range k (half-open); sorted
  // Otherwise, set #r is represented as a bitfield starting from [4], and
  //   [1]: offset of first bit (always divisible by 128)
  //   [2]: number of bits (divisible by 128 unless a variant in the last block
  //        is included)
  //   [3]: 1 if all out-of-bounds bits are set, 0 otherwise (other flags may
  //        be added later)
  // Empty sets are always stored in the first format, with [0] == 0.
  uint32_t** setdefs;
} Set_info;

#define ANNOT_NA 1
#define ANNOT_PRUNE 2
#define ANNOT_BLOCK 4
#define ANNOT_MINIMAL 8
#define ANNOT_DISTANCE 0x10

typedef struct {
  char* fname;
  char* attrib_fname;
  char* ranges_fname;
  char* filter_fname;
  char* snps_fname;
  char* subset_fname;
  char* snpfield;
  uint32_t modifier;
  uint32_t border;
} Annot_info;

void set_init(Set_info* sip, Annot_info* aip);

void set_cleanup(Set_info* sip, Annot_info* aip);

uint32_t in_setdef(uint32_t* setdef, uint32_t marker_idx);

uint32_t interval_in_setdef(uint32_t* setdef, uint32_t marker_idx_start, uint32_t marker_idx_end);

uint32_t setdef_size(uint32_t* setdef, uint32_t marker_ct);

void setdef_iter_init(uint32_t* setdef, uint32_t marker_ct, uint32_t start_idx, uint32_t* cur_idx_ptr, uint32_t* aux_ptr);

uint32_t setdef_iter(uint32_t* setdef, uint32_t* cur_idx_ptr, uint32_t* aux_ptr);

uint32_t alloc_and_populate_nonempty_set_incl(Set_info* sip, uint32_t* nonempty_set_ct_ptr, uintptr_t** nonempty_set_incl_ptr);

int32_t extract_exclude_range(char* fname, uint32_t* marker_pos, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t* marker_exclude_ct_ptr, uint32_t is_exclude, uint32_t allow_no_variants, Chrom_info* chrom_info_ptr);

uint32_t save_set_bitfield(uintptr_t* marker_bitfield_tmp, uint32_t marker_ct, uint32_t range_start, uint32_t range_end, uint32_t complement_sets, uint32_t** set_range_pp);

uint32_t save_set_range(uint64_t* range_sort_buf, uint32_t marker_ct, uint32_t rsb_last_idx, uint32_t complement_sets, uint32_t** set_range_pp);

int32_t define_sets(Set_info* sip, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* marker_pos, uintptr_t* marker_exclude_ct_ptr, char* marker_ids, uintptr_t max_marker_id_len, Chrom_info* chrom_info_ptr, uint32_t allow_no_variants);

int32_t write_set(Set_info* sip, char* outname, char* outname_end, uint32_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, char* marker_ids, uintptr_t max_marker_id_len, uint32_t* marker_pos, Chrom_info* chrom_info_ptr);

void unpack_set(uintptr_t marker_ct, uint32_t* setdef, uintptr_t* include_bitfield);

void unpack_set_unfiltered(uintptr_t marker_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uint32_t* setdef, uintptr_t* new_exclude);

uint32_t extract_set_union(uint32_t** setdefs, uintptr_t set_ct, uintptr_t* set_incl, uintptr_t* filtered_union, uintptr_t marker_ct);

uint32_t extract_set_union_unfiltered(Set_info* sip, uintptr_t* set_incl, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude, uintptr_t** union_marker_exclude_ptr, uintptr_t* union_marker_ct_ptr);

uint32_t setdefs_compress(Set_info* sip, uintptr_t* set_incl, uintptr_t set_ct, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uintptr_t marker_ct_orig, uintptr_t* marker_exclude, uintptr_t marker_ct, uint32_t*** new_setdefs_ptr);

int32_t load_range_list_sortpos(char* fname, uint32_t border_extend, uintptr_t subset_ct, char* sorted_subset_ids, uintptr_t max_subset_id_len, Chrom_info* chrom_info_ptr, uintptr_t* gene_ct_ptr, char** gene_names_ptr, uintptr_t* max_gene_id_len_ptr, uintptr_t** chrom_bounds_ptr, uint32_t*** genedefs_ptr, uintptr_t* chrom_max_gene_ct_ptr, const char* file_descrip);

int32_t annotate(const Annot_info* aip, uint32_t allow_extra_chroms, char* outname, char* outname_end, double pfilter, Chrom_info* chrom_info_ptr);

int32_t gene_report(char* fname, char* glist, char* subset_fname, uint32_t border, uint32_t allow_extra_chroms, char* extractname, const char* snp_field, char* outname, char* outname_end, double pfilter, Chrom_info* chrom_info_ptr);

#endif // __PLINK_SET_H__
