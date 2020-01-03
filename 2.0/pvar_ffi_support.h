#ifndef __PVAR_FFI_SUPPORT_H__
#define __PVAR_FFI_SUPPORT_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "include/pgenlib_misc.h"

#ifdef __cplusplus
namespace plink2 {
#endif

struct RefcountedWptrStruct {
  uintptr_t ref_ct;
  uintptr_t p[];
};

typedef struct RefcountedWptrStruct RefcountedWptr;

RefcountedWptr* CreateRefcountedWptr(uintptr_t size);

void CondReleaseRefcountedWptr(RefcountedWptr** rwpp);

// Minimal .pvar loader, using malloc/free instead of bigstack.  Necessary for
// clean multiallelic-variant support.
// Doesn't use plink2_decompress for now, since that has too many dependencies.
// (todo: remove zlibWrapper and plink2_cmdline dependencies from
// plink2_decompress.)
struct MinimalPvarStruct {
  const char** variant_ids;
  const char** allele_storage;
  RefcountedWptr* allele_idx_offsetsp;
  uint32_t variant_ct;
  uint32_t max_allele_ct;
};

typedef struct MinimalPvarStruct MinimalPvar;

void PreinitMinimalPvar(MinimalPvar* mpp);

PglErr LoadMinimalPvar(const char* fname, MinimalPvar* mpp, char* errstr_buf);

void CleanupMinimalPvar(MinimalPvar* mpp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PVAR_FFI_SUPPORT_H__
