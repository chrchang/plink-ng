// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
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

#include "plink2_float.h"

#ifdef __APPLE__
#  include <fenv.h>
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

void flush_denormals() {
#if defined(__x86_64__) || defined(__i386__)
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#else
#  ifdef __APPLE__
  fesetenv(FE_DFL_DISABLE_DENORMS_ENV);
#  else
  // generic arm64
  uint64_t fpcr;
  asm volatile("mrs %0, fpcr" : "=r"(fpcr));
  fpcr |= (1LLU << 24);
  asm volatile("msr fpcr, %0" : : "r"(fpcr));
#  endif
#endif
}

#ifdef __cplusplus
}
#endif
