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

#include "plink2_family.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitMendel(MendelInfo* mendel_info_ptr) {
  mendel_info_ptr->flags = kfMendel0;
  mendel_info_ptr->max_trio_error = 1.0;
  mendel_info_ptr->max_var_error = 1.0;
  mendel_info_ptr->exclude_one_ratio = 0.0;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
