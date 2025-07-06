#ifndef __PLINK2_FAMILY_H__
#define __PLINK2_FAMILY_H__

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

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfMendel0,
  kfMendelDuos = (1 << 0),
  kfMendelMultigen = (1 << 1),
  kfMendelFilterVarFirst = (1 << 2),
  kfMendelRptZs = (1 << 3),
  kfMendelRptSummariesOnly = (1 << 4),
  kfMendelRptColMaybefid = (1 << 5),
  kfMendelRptColFid = (1 << 6),
  kfMendelRptColMaybesid = (1 << 7),
  kfMendelRptColSid = (1 << 8),
  kfMendelRptColChrom = (1 << 9),
  kfMendelRptColPos = (1 << 10),
  kfMendelRptColId = (1 << 11),
  kfMendelRptColCode = (1 << 12),
  kfMendelRptColError = (1 << 13),
  kfMendelRptColDefault = (kfMendelRptColMaybefid | kfMendelRptColMaybesid | kfMendelRptColChrom | kfMendelRptColId | kfMendelRptColCode | kfMendelRptColError),
  kfMendelRptColAll = ((kfMendelRptColError * 2) - kfMendelRptColMaybefid)
FLAGSET_DEF_END(MendelFlags);

typedef struct MendelStruct {
  NONCOPYABLE(MendelStruct);
  MendelFlags flags;
  double max_trio_error;
  double max_var_error;
  double exclude_one_ratio;  // 0 = off, DBL_MAX = always child
} MendelInfo;

void InitMendel(MendelInfo* mendel_info_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_FAMILY_H__
