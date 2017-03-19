#ifndef __PLINK2_ADJUST_H__
#define __PLINK2_ADJUST_H__

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

FLAGSET_DEF_START()
  kfAdjust0,
  kfAdjustGc = (1 << 0),
  kfAdjustLog10 = (1 << 1),

  kfAdjustColChrom = (1 << 2),
  kfAdjustColPos = (1 << 3),
  kfAdjustColRef = (1 << 4),
  kfAdjustColAlt1 = (1 << 5),
  kfAdjustColAlt = (1 << 6),
  kfAdjustColUnadj = (1 << 7),
  kfAdjustColGc = (1 << 8),
  kfAdjustColQq = (1 << 9),
  kfAdjustColBonf = (1 << 10),
  kfAdjustColHolm = (1 << 11),
  kfAdjustColSidakss = (1 << 12),
  kfAdjustColSidaksd = (1 << 13),
  kfAdjustColFdrbh = (1 << 14),
  kfAdjustColFdrby = (1 << 15),
  kfAdjustColDefault = (kfAdjustColChrom | kfAdjustColUnadj | kfAdjustColGc | kfAdjustColBonf | kfAdjustColHolm | kfAdjustColSidakss | kfAdjustColSidaksd | kfAdjustColFdrbh | kfAdjustColFdrby),
  kfAdjustColAll = ((kfAdjustColFdrby * 2) - kfAdjustColChrom)
FLAGSET_DEF_END(adjust_flags_t);

typedef struct adjust_info_struct {
  adjust_flags_t flags;
  double lambda;
} adjust_info_t;

void init_adjust(adjust_info_t* adjust_info_ptr);

#ifdef __cplusplus
} // namespace plink2
#endif
 
#endif // __PLINK2_ADJUST_H__
