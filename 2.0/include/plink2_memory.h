#ifndef __PLINK2_MEMORY_H__
#define __PLINK2_MEMORY_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.


// Double-ended stack allocator, and a few support functions (e.g.
// GetMemAvailableKib to estimate how much memory we can afford to use).

#include <assert.h>

#include "plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

uint64_t DetectMib();

// On more recent Linux builds, the MemAvailable value reported by
// "cat /proc/meminfo" is a pretty accurate estimate of how much we can afford
// to allocate before introducing substantial risk of OOM-killer action.
// Returns (~0LLU) if this estimate is unavailable.
// Requires a text buffer large enough to fit longest line in /proc/meminfo.
uint64_t GetMemAvailableKib(uint32_t textbuf_blen, char* textbuf);

// ensure vector-alignment
CONSTI32(kEndAllocAlign, MAXV(kBytesPerVec, 16));

HEADER_INLINE void* arena_alloc(unsigned char* arena_top, uintptr_t size, unsigned char** arena_bottom_ptr) {
  size = RoundUpPow2(size, kCacheline);
  if (S_CAST(uintptr_t, arena_top - (*arena_bottom_ptr)) < size) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return arena_alloc_raw(size, arena_bottom_ptr);
}

HEADER_INLINE BoolErr arena_alloc_c(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, char** c_arr_ptr) {
  *c_arr_ptr = S_CAST(char*, arena_alloc(arena_top, ct, arena_bottom_ptr));
  return !(*c_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_d(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, double** d_arr_ptr) {
  *d_arr_ptr = S_CAST(double*, arena_alloc(arena_top, ct * sizeof(double), arena_bottom_ptr));
  return !(*d_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_f(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, float** f_arr_ptr) {
  *f_arr_ptr = S_CAST(float*, arena_alloc(arena_top, ct * sizeof(float), arena_bottom_ptr));
  return !(*f_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_i32(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, int32_t** i32_arr_ptr) {
  *i32_arr_ptr = S_CAST(int32_t*, arena_alloc(arena_top, ct * sizeof(int32_t), arena_bottom_ptr));
  return !(*i32_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_uc(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = S_CAST(unsigned char*, arena_alloc(arena_top, ct, arena_bottom_ptr));
  return !(*uc_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_u32(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, uint32_t** u32_arr_ptr) {
  *u32_arr_ptr = S_CAST(uint32_t*, arena_alloc(arena_top, ct * sizeof(int32_t), arena_bottom_ptr));
  return !(*u32_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_w(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, uintptr_t** w_arr_ptr) {
  *w_arr_ptr = S_CAST(uintptr_t*, arena_alloc(arena_top, ct * sizeof(intptr_t), arena_bottom_ptr));
  return !(*w_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_i64(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, int64_t** i64_arr_ptr) {
  *i64_arr_ptr = S_CAST(int64_t*, arena_alloc(arena_top, ct * sizeof(int64_t), arena_bottom_ptr));
  return !(*i64_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_u64(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, uint64_t** u64_arr_ptr) {
  *u64_arr_ptr = S_CAST(uint64_t*, arena_alloc(arena_top, ct * sizeof(int64_t), arena_bottom_ptr));
  return !(*u64_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_cp(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, char*** cp_arr_ptr) {
  *cp_arr_ptr = S_CAST(char**, arena_alloc(arena_top, ct * sizeof(intptr_t), arena_bottom_ptr));
  return !(*cp_arr_ptr);
}

HEADER_INLINE BoolErr arena_alloc_kcp(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, const char*** kcp_arr_ptr) {
  *kcp_arr_ptr = S_CAST(const char**, arena_alloc(arena_top, ct * sizeof(intptr_t), arena_bottom_ptr));
  return !(*kcp_arr_ptr);
}

BoolErr arena_calloc_w(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, uintptr_t** w_arr_ptr);

HEADER_INLINE void ArenaBaseSet(const void* unaligned_base, unsigned char** arena_bottom_ptr) {
  *arena_bottom_ptr = R_CAST(unsigned char*, RoundUpPow2(R_CAST(uintptr_t, unaligned_base), kCacheline));
}

HEADER_INLINE void ArenaEndSet(const void* unaligned_end, unsigned char** arena_top_ptr) {
  *arena_top_ptr = R_CAST(unsigned char*, RoundDownPow2(R_CAST(uintptr_t, unaligned_end), kEndAllocAlign));
}

HEADER_INLINE void* arena_end_alloc_raw(uintptr_t size, unsigned char** arena_top_ptr) {
  assert(!(size % kEndAllocAlign));
  unsigned char* alloc_ptr = *arena_top_ptr;
  alloc_ptr -= size;
  *arena_top_ptr = alloc_ptr;
  return alloc_ptr;
}

HEADER_INLINE void* arena_end_alloc(unsigned char* arena_bottom, uintptr_t size, unsigned char** arena_top_ptr) {
  size = RoundUpPow2(size, kEndAllocAlign);
  if (unlikely(S_CAST(uintptr_t, (*arena_top_ptr) - arena_bottom) < size)) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return arena_end_alloc_raw(size, arena_top_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_c(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, char** c_arr_ptr) {
  *c_arr_ptr = S_CAST(char*, arena_end_alloc(arena_bottom, ct, arena_top_ptr));
  return !(*c_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_d(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, double** d_arr_ptr) {
  *d_arr_ptr = S_CAST(double*, arena_end_alloc(arena_bottom, ct * sizeof(double), arena_top_ptr));
  return !(*d_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_f(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, float** f_arr_ptr) {
  *f_arr_ptr = S_CAST(float*, arena_end_alloc(arena_bottom, ct * sizeof(float), arena_top_ptr));
  return !(*f_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_i32(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, int32_t** i32_arr_ptr) {
  *i32_arr_ptr = S_CAST(int32_t*, arena_end_alloc(arena_bottom, ct * sizeof(int32_t), arena_top_ptr));
  return !(*i32_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_uc(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = S_CAST(unsigned char*, arena_end_alloc(arena_bottom, ct, arena_top_ptr));
  return !(*uc_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_u32(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, uint32_t** u32_arr_ptr) {
  *u32_arr_ptr = S_CAST(uint32_t*, arena_end_alloc(arena_bottom, ct * sizeof(int32_t), arena_top_ptr));
  return !(*u32_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_w(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, uintptr_t** w_arr_ptr) {
  *w_arr_ptr = S_CAST(uintptr_t*, arena_end_alloc(arena_bottom, ct * sizeof(intptr_t), arena_top_ptr));
  return !(*w_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_i64(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, int64_t** i64_arr_ptr) {
  *i64_arr_ptr = S_CAST(int64_t*, arena_end_alloc(arena_bottom, ct * sizeof(int64_t), arena_top_ptr));
  return !(*i64_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_u64(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, uint64_t** u64_arr_ptr) {
  *u64_arr_ptr = S_CAST(uint64_t*, arena_end_alloc(arena_bottom, ct * sizeof(int64_t), arena_top_ptr));
  return !(*u64_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_cp(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, char*** cp_arr_ptr) {
  *cp_arr_ptr = S_CAST(char**, arena_end_alloc(arena_bottom, ct * sizeof(intptr_t), arena_top_ptr));
  return !(*cp_arr_ptr);
}

HEADER_INLINE BoolErr arena_end_alloc_kcp(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, const char*** kcp_arr_ptr) {
  *kcp_arr_ptr = S_CAST(const char**, arena_end_alloc(arena_bottom, ct * sizeof(intptr_t), arena_top_ptr));
  return !(*kcp_arr_ptr);
}

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_MEMORY_H__
