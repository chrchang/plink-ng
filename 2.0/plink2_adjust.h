#ifndef __PLINK2_ADJUST_H__
#define __PLINK2_ADJUST_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfAdjust0,
  kfAdjustGc = (1 << 0),
  kfAdjustLog10 = (1 << 1),
  kfAdjustZs = (1 << 2),
  kfAdjustInputLog10 = (1 << 3),

  kfAdjustColChrom = (1 << 4),
  kfAdjustColPos = (1 << 5),
  kfAdjustColRef = (1 << 6),
  kfAdjustColAlt1 = (1 << 7),
  kfAdjustColAlt = (1 << 8),
  kfAdjustColUnadj = (1 << 9),
  kfAdjustColGc = (1 << 10),
  kfAdjustColQq = (1 << 11),
  kfAdjustColBonf = (1 << 12),
  kfAdjustColHolm = (1 << 13),
  kfAdjustColSidakss = (1 << 14),
  kfAdjustColSidaksd = (1 << 15),
  kfAdjustColFdrbh = (1 << 16),
  kfAdjustColFdrby = (1 << 17),
  kfAdjustColDefault = (kfAdjustColChrom | kfAdjustColUnadj | kfAdjustColGc | kfAdjustColBonf | kfAdjustColHolm | kfAdjustColSidakss | kfAdjustColSidaksd | kfAdjustColFdrbh | kfAdjustColFdrby),
  kfAdjustColAll = ((kfAdjustColFdrby * 2) - kfAdjustColChrom)
FLAGSET_DEF_END(adjust_flags_t);

typedef struct adjust_info_struct {
  adjust_flags_t flags;
  double lambda;
} adjust_info_t;

typedef struct adjust_file_info_struct {
  adjust_info_t base;
  char* fname;
  char* test_name;
  char* chr_field;
  char* pos_field;
  char* id_field;
  char* ref_field;
  char* alt_field;
  char* test_field;
  char* p_field;
  char* chisq_field;
} adjust_file_info_t;

void init_adjust(adjust_info_t* adjust_info_ptr, adjust_file_info_t* adjust_file_info_ptr);

void cleanup_adjust(adjust_file_info_t* adjust_file_info_ptr);

pglerr_t multcomp(const uintptr_t* variant_include, const chr_info_t* cip, const char* const* chr_ids, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* variant_allele_idxs, const char* const* allele_storage, const adjust_info_t* adjust_info_ptr, const double* pvals, const double* chisqs, uint32_t orig_variant_ct, uint32_t max_allele_slen, double pfilter, double output_min_p, uint32_t skip_gc, uint32_t max_thread_ct, char* outname, char* outname_end);

pglerr_t adjust_file(const adjust_file_info_t* afip, double pfilter, double output_min_p, uint32_t max_thread_ct, char* outname, char* outname_end);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_ADJUST_H__
