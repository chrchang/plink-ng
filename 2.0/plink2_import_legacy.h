#ifndef __PLINK2_IMPORT_LEGACY_H__
#define __PLINK2_IMPORT_LEGACY_H__

// This file is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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

#include "include/plink2_base.h"
#include "plink2_common.h"
#include "plink2_data.h"

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr LoadMap(const char* mapname, MiscFlags misc_flags, LoadFilterLogFlags load_filter_log_import_flags, ChrInfo* cip, uint32_t* max_variant_id_slen_ptr, ChrIdx** variant_chr_codes_ptr, uint32_t** variant_bps_ptr, char*** variant_ids_ptr, double** variant_cms_ptr, uint32_t* variant_ct_ptr);

PglErr TpedToPgen(const char* tpedname, const char* tfamname, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, FamCol fam_cols, int32_t missing_pheno, char input_missing_geno_char, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* psam_generated_ptr);

PglErr Plink1SampleMajorToPgen(const char* pgenname, const uintptr_t* allele_flips, uintptr_t variant_ct, uintptr_t sample_ct, uint32_t real_ref_alleles, uint32_t max_thread_ct, FILE* infile);

PglErr PedmapToPgen(const char* pedname, const char* mapname, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, uint32_t psam_01, FamCol fam_cols, int32_t missing_pheno, char input_missing_geno_char, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_IMPORT_LEGACY_H__
