#ifndef __PLINK2_HELP_H__
#define __PLINK2_HELP_H__

// This file is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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


#include "plink2_cmdline.h"

#ifdef __cplusplus
namespace plink2 {
#endif

extern const char kCmdlineFormatStr[];

PglErr DispHelp(const char* const* argvk, uint32_t param_ct);

#ifdef __cplusplus
}
#endif

#endif  // __PLINK2_HELP_H__
