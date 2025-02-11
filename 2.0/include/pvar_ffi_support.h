#ifndef __PVAR_FFI_SUPPORT_H__
#define __PVAR_FFI_SUPPORT_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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

#include "plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Would prefer to make p a flexible array member, but that doesn't adhere to
// current CRAN coding standards, so we have an extra dereference here.
struct RefcountedWptrStruct {
  uintptr_t ref_ct;
  uintptr_t* p;
};

typedef struct RefcountedWptrStruct RefcountedWptr;

RefcountedWptr* CreateRefcountedWptr(uintptr_t size);

void CondReleaseRefcountedWptr(RefcountedWptr** rwpp);

// Only enforced for chr_names for now.
CONSTI32(kMaxIdSlen, 16000);

// This is based on plink2 defaults: 65376 = 1 unplaced + 95 autosomes + 6
// special cases (X/Y/XY/MT/PAR1/PAR2) + 65274 contigs
CONSTI32(kMaxChromosomes, 65376);
// We use a simple fixed-size hash table for the chromosome names.  Ensure load
// factor doesn't exceed 0.5 (of course it's almost always tiny in practice).
CONSTI32(kChrHtableSize, 130752);
typedef uint16_t ChrIdx;

// Minimal .pvar loader (CHROM/POS optional; no QUAL/FILTER/INFO), using
// malloc/free instead of bigstack.  Necessary for clean multiallelic-variant
// support.
struct MinimalPvarStruct {
  const char** chr_names;
  const ChrIdx* chr_idxs;
  const int32_t* variant_bps;
  const char** variant_ids;
  const char** allele_storage;  // part of varisnt_ids allocation

  // reference-counted since ownership can be shared with .pgen object
  RefcountedWptr* allele_idx_offsetsp;

  uint32_t chr_ct;
  uint32_t variant_ct;
  uint32_t max_allele_ct;
};

typedef struct MinimalPvarStruct MinimalPvar;

void PreinitMinimalPvar(MinimalPvar* mpp);

FLAGSET_DEF_START()
  kfLoadMinimalPvar0,
  kfLoadMinimalPvarOmitChrom = (1 << 0),
  kfLoadMinimalPvarOmitPos = (1 << 1)
FLAGSET_DEF_END(LoadMinimalPvarFlags);

PglErr LoadMinimalPvarEx(const char* fname, LoadMinimalPvarFlags flags, MinimalPvar* mpp, char* errstr_buf);

HEADER_INLINE PglErr LoadMinimalPvar(const char* fname, MinimalPvar* mpp, char* errstr_buf) {
  return LoadMinimalPvarEx(fname, kfLoadMinimalPvarOmitChrom | kfLoadMinimalPvarOmitPos, mpp, errstr_buf);
}

void CleanupMinimalPvar(MinimalPvar* mpp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PVAR_FFI_SUPPORT_H__
