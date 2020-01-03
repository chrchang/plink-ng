#ifndef __PLINK2_PSAM_H__
#define __PLINK2_PSAM_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.


#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Key .psam properties: (i) .fam files are valid .psam files; (ii) .fam-like
// files with an additional 'SID' column ("source ID", disambiguating multiple
// samples from e.g. the same cancer patient) are valid; (iii) zero, or many,
// phenotypes are now ok; (iv) Oxford .sample files can be converted to .psam
// without loss of information; and (v) loader should work on .ped files as
// well.

// File starts with an arbitrary (possibly zero) number of header lines marked
// by a leading '#'.  All lines which don't start with '#FID' or '#IID' are
// currently ignored.  The #FID/IID line specifies the columns in the .psam
// file; the following column headers are recognized:
//   IID (required)
//   SID
//   PAT
//   MAT
//   SEX
// FID must either be the first column, or absent.  If it's absent, all FID
// values are now set to '0' by default, and FID is omitted by default from
// output files when there's a choice.  (This is a change from both v1.9 and
// v2.0 alpha 1.)
// Any other value is treated as a phenotype/covariate name.
//
// The loader will error out of multiple #FID/IID lines are in the header for
// some bizarre reason.  If no #FID/IID line is present, fam_cols controls the
// default (e.g. fam_cols == FAM_COL_13456 means "#FID IID PAT MAT SEX PHENO").


// Memory for all the return arrays is allocated off the bottom of g_bigstack.

// todo: the new indiv_represent[] bitvector tracks which sample should be
// considered the primary one for an individual, when there are multiple
// samples with identical FID/IID but differing SID.  (Default is the first
// sample for each individual, but if a later sample has SID == '0' that takes
// precedence.)
// plink2 order of operations should be arranged so that FID/IID-insensitive
// sample filters happen first, then indiv_represent is computed, then
// FID/IID-sensitive sample filters are applied.

// chosen to be likely to fit in L3 cache
CONSTI32(kCatHtableSize, 524287);
static_assert(kCatHtableSize >= kMaxPhenoCt, "kCatHtableSize cannot be smaller than kMaxPhenoCt.");

PglErr LoadPsam(const char* psamname, const RangeList* pheno_range_list_ptr, FamCol fam_cols, uint32_t pheno_ct_max, int32_t missing_pheno, uint32_t affection_01, uint32_t max_thread_ct, PedigreeIdInfo* piip, uintptr_t** sample_include_ptr, uintptr_t** founder_info_ptr, uintptr_t** sex_nm_ptr, uintptr_t** sex_male_ptr, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* raw_sample_ct_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr);

HEADER_INLINE BoolErr IsReservedPhenoName(const char* pheno_name, uint32_t pheno_name_slen) {
  if (pheno_name_slen != 3) {
    return 0;
  }
  // tolerate "SEX" column in phenotype/covariate files; just impose some
  // restrictions on it when writing .psam files.
  return memequal_k(pheno_name, "FID", 3) || memequal_k(pheno_name, "IID", 3) || memequal_k(pheno_name, "SID", 3) || memequal_k(pheno_name, "PAT", 3) || memequal_k(pheno_name, "MAT", 3);
}

// also for loading covariates.  set affection_01 to 2 to prohibit case/control
PglErr LoadPhenos(const char* pheno_fname, const RangeList* pheno_range_list_ptr, const uintptr_t* sample_include, const char* sample_ids, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, int32_t missing_pheno, uint32_t affection_01, uint32_t iid_only, uint32_t numeric_ranges, uint32_t max_thread_ct, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_PSAM_H__
