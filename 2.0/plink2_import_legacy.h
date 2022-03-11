#ifndef __PLINK2_IMPORT_LEGACY_H__
#define __PLINK2_IMPORT_LEGACY_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2022 Shaun Purcell,
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

PglErr PedmapToPgen(const char* pedname, const char* mapname, MiscFlags misc_flags, ImportFlags import_flags, uint32_t psam_01, FamCol fam_cols, int32_t missing_pheno, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip);

PglErr TpedToPgen(const char* tpedname, const char* tfamname, MiscFlags misc_flags, ImportFlags import_flags, FamCol fam_cols, int32_t missing_pheno, char input_missing_geno_char, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* psam_generated_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_IMPORT_LEGACY_H__
