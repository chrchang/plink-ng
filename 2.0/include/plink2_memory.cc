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

#include "plink2_memory.h"

#include <unistd.h>  // sysconf()
#ifdef __APPLE__
#  include <sys/sysctl.h>  // sysctl()
#endif

#include "plink2_string.h"

#ifdef __cplusplus
namespace plink2 {
#endif

uint64_t DetectMib() {
  int64_t llxx;
  // return zero if detection failed
  // see e.g. http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system .
#ifdef __APPLE__
  int32_t mib[2];
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  llxx = 0;
  size_t sztmp = sizeof(int64_t);
  sysctl(&mib[0], 2, &llxx, &sztmp, nullptr, 0);
  llxx /= 1048576;
#else
#  ifdef _WIN32
  MEMORYSTATUSEX memstatus;
  memstatus.dwLength = sizeof(memstatus);
  GlobalMemoryStatusEx(&memstatus);
  llxx = memstatus.ullTotalPhys / 1048576;
#  else
  // TODO: check whether getrlimit(RLIMIT_AS, ...) returns a per-user limit on
  // shared systems.
  llxx = S_CAST(uint64_t, sysconf(_SC_PHYS_PAGES)) * S_CAST(size_t, sysconf(_SC_PAGESIZE)) / 1048576;
#  endif
#endif
  return llxx;
}

uint64_t GetMemAvailableKib(__attribute__((unused)) uint32_t textbuf_blen, __attribute__((unused)) char* textbuf) {
  // On more recent Linux builds, the MemAvailable value reported by
  // "cat /proc/meminfo" is a pretty accurate estimate of how much we can
  // afford to allocate before introducing substantial risk of OOM-killer
  // action.
  // Returns (~0LLU) if this estimate is unavailable.
#if defined(__APPLE__) || defined(_WIN32)
  return (~0LLU);
#else
  FILE* meminfo = fopen("/proc/meminfo", FOPEN_RB);
  if (!meminfo) {
    return (~0LLU);
  }
  {
    textbuf[textbuf_blen - 1] = ' ';
    do {
      if ((!fgets(textbuf, textbuf_blen, meminfo)) || (!textbuf[textbuf_blen - 1])) {
        goto GetMemAvailableKib_ret_FAIL;
      }
    } while (!StrStartsWithUnsafe(textbuf, "MemAvailable:"));
    const char* textbuf_iter = &(textbuf[strlen("MemAvailable:")]);
    textbuf_iter = FirstNonTspace(textbuf_iter);
    uint64_t kib_free;
    if (ScanmovU64Capped(1LLU << 54, &textbuf_iter, &kib_free)) {
      goto GetMemAvailableKib_ret_FAIL;
    }
    textbuf_iter = FirstNonTspace(textbuf_iter);
    if (!memequal_sk(textbuf_iter, "kB\n")) {
      goto GetMemAvailableKib_ret_FAIL;
    }
    return kib_free;
  }
 GetMemAvailableKib_ret_FAIL:
  fclose(meminfo);
  return (~0LLU);
#endif
}

BoolErr arena_calloc_w(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, uintptr_t** w_arr_ptr) {
  *w_arr_ptr = S_CAST(uintptr_t*, arena_alloc(arena_top, ct * sizeof(intptr_t), arena_bottom_ptr));
  if (unlikely(!(*w_arr_ptr))) {
    return 1;
  }
  ZeroWArr(ct, *w_arr_ptr);
  return 0;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
