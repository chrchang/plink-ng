#ifndef __PLINK2_EXPORT_H__
#define __PLINK2_EXPORT_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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

ENUM_U31_DEF_START()
  kVcfExport0,
  kVcfExportGp,
  kVcfExportDs,
  kVcfExportDsForce,
  kVcfExportHds,
  kVcfExportHdsForce
ENUM_U31_DEF_END(VcfExportMode);

FLAGSET_DEF_START()
  kfIdpaste0,
  kfIdpasteMaybefid = (1 << 0),
  kfIdpasteFid = (1 << 1),
  kfIdpasteIid = (1 << 2),
  kfIdpasteMaybesid = (1 << 3),
  kfIdpasteSid = (1 << 4),
  kfIdpasteDefault = (kfIdpasteMaybefid | kfIdpasteIid | kfIdpasteMaybesid)
FLAGSET_DEF_END(IdpasteFlags);

typedef struct ExportfStruct {
  ExportfFlags flags;
  IdpasteFlags idpaste_flags;
  char id_delim;
  uint32_t bgen_bits;
  VcfExportMode vcf_mode;
  char* export_allele_fname;
} ExportfInfo;

void InitExportf(ExportfInfo* exportf_info_ptr);

void CleanupExportf(ExportfInfo* exportf_info_ptr);

PglErr Exportf(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const ExportfInfo* eip, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, UnsortedVar vpos_sortstatus, uint32_t max_thread_ct, MakePlink2Flags make_plink2_flags, uintptr_t pgr_alloc_cacheline_ct, char* xheader, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_EXPORT_H__
