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


#include "plink2_import.h"

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr PedmapToPgen(__attribute__((unused)) const char* pedname, __attribute__((unused)) const char* mapname, __attribute__((unused)) MiscFlags misc_flags, __attribute__((unused)) FamCol fam_cols, __attribute__((unused)) int32_t missing_pheno, __attribute__((unused)) uint32_t max_thread_ct, __attribute__((unused)) char* outname, __attribute__((unused)) char* outname_end, __attribute__((unused)) ChrInfo* cip) {
  logerrputs("Error: .ped import is under development.\n");
  return kPglRetNotYetSupported;
}

PglErr TpedToPgen(__attribute__((unused)) const char* tpedname, __attribute__((unused)) const char* tfamname, __attribute__((unused)) MiscFlags misc_flags, __attribute__((unused)) FamCol fam_cols, __attribute__((unused)) int32_t missing_pheno, __attribute__((unused)) uint32_t max_thread_ct, __attribute__((unused)) char* outname, __attribute__((unused)) char* outname_end, __attribute__((unused)) ChrInfo* cip) {
  logerrputs("Error: .tped import is under development.\n");
  return kPglRetNotYetSupported;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
