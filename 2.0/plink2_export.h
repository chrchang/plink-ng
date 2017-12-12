#ifndef __PLINK2_EXPORT_H__
#define __PLINK2_EXPORT_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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


#include "plink2_data.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfIdpaste0,
  kfIdpasteFid = (1 << 0),
  kfIdpasteIid = (1 << 1),
  kfIdpasteMaybesid = (1 << 2),
  kfIdpasteSid = (1 << 3),
  kfIdpasteDefault = (kfIdpasteFid | kfIdpasteIid | kfIdpasteMaybesid)
FLAGSET_DEF_END(idpaste_t);

pglerr_t exportf(char* xheader, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen,  uint32_t max_thread_ct, make_plink2_t make_plink2_modifier, exportf_flags_t exportf_modifier, idpaste_t exportf_id_paste, char exportf_id_delim, uint32_t exportf_bits, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_EXPORT_H__
