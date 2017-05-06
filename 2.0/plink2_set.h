#ifndef __PLINK2_SET_H__
#define __PLINK2_SET_H__

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


#include "plink2_common.h"

#ifdef __cplusplus
namespace plink2 {
#endif

typedef struct make_set_range_struct {
  struct make_set_range_struct* next;
  uint32_t uidx_start;
  uint32_t uidx_end;
} make_set_range_t;

pglerr_t extract_exclude_range(const char* fname, const chr_info_t* cip, const uint32_t* variant_bps, uint32_t raw_variant_ct, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr);

#ifdef __cplusplus
}
#endif

#endif // __PLINK2_SET_H__
