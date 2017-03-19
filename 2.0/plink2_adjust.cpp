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

#include "plink2_glm.h"
#include "plink2_matrix.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void init_adjust(adjust_info_t* adjust_info_ptr) {
  adjust_info_ptr->flags = kfAdjust0;
  adjust_info_ptr->lambda = 0.0;
}

typedef struct adjustable_assoc_result_struct {
  double pval;
  uint32_t variant_uidx;
#ifdef __cplusplus
  bool operator<(const struct adjustable_assoc_result_struct& rhs) const {
    return pval < rhs.pval;
  }
#endif
} adjustable_assoc_result_t;

#ifdef __cplusplus
} // namespace plink2
#endif
