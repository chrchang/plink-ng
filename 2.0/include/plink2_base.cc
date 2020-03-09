// This library is part of PLINK 2, copyright (C) 2005-2020 Shaun Purcell,
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

uintptr_t g_failed_alloc_attempt_size = 0;

#if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 7) && !defined(__APPLE__)
BoolErr pgl_malloc(uintptr_t size, void* pp) {
  *S_CAST(unsigned char**, pp) = S_CAST(unsigned char*, malloc(size));
  if (likely(*S_CAST(unsigned char**, pp))) {
    return 0;
  }
  g_failed_alloc_attempt_size = size;
  return 1;
}
#endif

BoolErr fwrite_checked(const void* buf, uintptr_t len, FILE* outfile) {
  while (len > kMaxBytesPerIO) {
    // OS X fwrite() doesn't support 2GiB+ writes
    // typical disk block size is 4kb, so 0x7ffff000 is the largest sensible
    // write size
    // bugfix (9 Mar 2018): forgot a 'not' here...
    if (unlikely(!fwrite_unlocked(buf, kMaxBytesPerIO, 1, outfile))) {
      return 1;
    }
    buf = &(S_CAST(const unsigned char*, buf)[kMaxBytesPerIO]);
    len -= kMaxBytesPerIO;
  }
  uintptr_t written_byte_ct = fwrite_unlocked(buf, 1, len, outfile);
  // must do the right thing when len == 0
  return (written_byte_ct != len);
}

/*
IntErr fread_checked2(void* buf, uintptr_t len, FILE* infile, uintptr_t* bytes_read_ptr) {
  uintptr_t bytes_read = 0;
  while (len > kMaxBytesPerIO) {
    const uintptr_t cur_bytes_read = fread_unlocked(buf, 1, kMaxBytesPerIO, infile);
    bytes_read += cur_bytes_read;
    if (cur_bytes_read != kMaxBytesPerIO) {
      *bytes_read_ptr = bytes_read;
      return ferror_unlocked(infile);
    }
    buf = &(((char*)buf)[kMaxBytesPerIO]);
    len -= kMaxBytesPerIO;
  }
  bytes_read += fread_unlocked(buf, 1, len, infile);
  *bytes_read_ptr = bytes_read;
  // could skip ferror_unlocked call if bytes_read == original len
  return ferror_unlocked(infile);
}
*/

BoolErr fread_checked(void* buf, uintptr_t len, FILE* infile) {
  while (len > kMaxBytesPerIO) {
    const uintptr_t cur_bytes_read = fread_unlocked(buf, 1, kMaxBytesPerIO, infile);
    if (unlikely(cur_bytes_read != kMaxBytesPerIO)) {
      return 1;
    }
    buf = &(S_CAST(unsigned char*, buf)[kMaxBytesPerIO]);
    len -= kMaxBytesPerIO;
  }
  const uintptr_t cur_bytes_read = fread_unlocked(buf, 1, len, infile);
  return (cur_bytes_read != len);
}

#ifdef __LP64__
static inline BoolErr ScanUintCappedFinish(const char* str_iter, uint64_t cap, uint32_t* valp) {
  uint64_t val = *valp;
  while (1) {
    // a little bit of unrolling seems to help
    const uint64_t cur_digit = ctou64(*str_iter++) - 48;
    if (cur_digit >= 10) {
      break;
    }
    // val = val * 10 + cur_digit;
    const uint64_t cur_digit2 = ctou64(*str_iter++) - 48;
    if (cur_digit2 >= 10) {
      val = val * 10 + cur_digit;
      if (unlikely(val > cap)) {
        return 1;
      }
      break;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
    if (unlikely(val > cap)) {
      return 1;
    }
  }
  *valp = val;
  return 0;
}

BoolErr ScanPosintCapped(const char* str_iter, uint64_t cap, uint32_t* valp) {
  // '0' has ascii code 48
  assert(ctou32(str_iter[0]) > 32);
  *valp = ctou32(*str_iter++) - 48;
  if (*valp >= 10) {
    // permit leading '+' (ascii 43), but not '++' or '+-'
    // reasonable to use unlikely() here since these functions aren't used for
    // numeric vs. non-numeric classification anyway due to erroring out on
    // overflow
    if (unlikely(*valp != 0xfffffffbU)) {
      return 1;
    }
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely(*valp >= 10)) {
      return 1;
    }
  }
  while (!(*valp)) {
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely((*valp) >= 10)) {
      return 1;
    }
  }
  return ScanUintCappedFinish(str_iter, cap, valp);
}

// Note that NumericRangeListToBitarr() can call this in an ignore-overflow
// mode.  If similar logic ever goes into an inner loop, remove all unlikely()
// annotations in this function and its children.
BoolErr ScanUintCapped(const char* str_iter, uint64_t cap, uint32_t* valp) {
  // Reads an integer in [0, cap].  Assumes first character is nonspace.
  assert(ctou32(str_iter[0]) > 32);
  uint32_t val = ctou32(*str_iter++) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      // '-' has ascii code 45, so unsigned 45 - 48 = 0xfffffffdU
      if (unlikely((val != 0xfffffffdU) || (*str_iter != '0'))) {
        return 1;
      }
      // accept "-0", "-00", etc.
      while (*(++str_iter) == '0');
      *valp = 0;
      return (ctou32(*str_iter) - 48) < 10;
    }
    // accept leading '+'
    val = ctou32(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  *valp = val;
  return ScanUintCappedFinish(str_iter, cap, valp);
}

BoolErr ScanIntAbsBounded(const char* str_iter, uint64_t bound, int32_t* valp) {
  // Reads an integer in [-bound, bound].  Assumes first character is nonspace.
  assert(ctou32(str_iter[0]) > 32);
  *valp = ctou32(*str_iter++) - 48;
  int32_t sign = 1;
  if (ctou32(*valp) >= 10) {
    if (*valp == -3) {
      sign = -1;
    } else if (unlikely(*valp != -5)) {
      return 1;
    }
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely(*valp >= 10)) {
      return 1;
    }
  }
  if (unlikely(ScanUintCappedFinish(str_iter, bound, R_CAST(uint32_t*, valp)))) {
    return 1;
  }
  *valp *= sign;
  return 0;
}
#else  // not __LP64__
BoolErr ScanPosintCapped32(const char* str_iter, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp) {
  // '0' has ascii code 48
  assert(ctou32(str_iter[0]) > 32);
  uint32_t val = ctou32(*str_iter++) - 48;
  if (val >= 10) {
    if (unlikely(val != 0xfffffffbU)) {
      return 1;
    }
    val = ctou32(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  while (!val) {
    val = ctou32(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  for (; ; ++str_iter) {
    const uint32_t cur_digit = ctou32(*str_iter) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      return 0;
    }
    // avoid integer overflow in middle of computation
    if (unlikely((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10)))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}

BoolErr ScanUintCapped32(const char* str_iter, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp) {
  // Reads an integer in [0, cap].  Assumes first character is nonspace.
  assert(ctou32(str_iter[0]) > 32);
  uint32_t val = ctou32(*str_iter++) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      if (unlikely((val != 0xfffffffdU) || (*str_iter != '0'))) {
        return 1;
      }
      while (*(++str_iter) == '0');
      *valp = 0;
      return (ctou32(*str_iter) - 48) < 10;
    }
    val = ctou32(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  for (; ; ++str_iter) {
    const uint32_t cur_digit = ctou32(*str_iter) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      return 0;
    }
    if (unlikely((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10)))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}

BoolErr ScanIntAbsBounded32(const char* str_iter, uint32_t bound_div_10, uint32_t bound_mod_10, int32_t* valp) {
  // Reads an integer in [-bound, bound].  Assumes first character is nonspace.
  assert(ctou32(str_iter[0]) > 32);
  uint32_t val = ctou32(*str_iter++) - 48;
  int32_t sign = 1;
  if (val >= 10) {
    if (val == 0xfffffffdU) {
      sign = -1;
    } else if (unlikely(val != 0xfffffffbU)) {
      return 1;
    }
    val = ctou32(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  for (; ; ++str_iter) {
    const uint32_t cur_digit = ctou32(*str_iter) - 48;
    if (cur_digit >= 10) {
      *valp = sign * S_CAST(int32_t, val);
      return 0;
    }
    if (unlikely((val >= bound_div_10) && ((val > bound_div_10) || (cur_digit > bound_mod_10)))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}
#endif

BoolErr aligned_malloc(uintptr_t size, uintptr_t alignment, void* aligned_pp) {
  // Assumes malloc returns word-aligned addresses.
  assert(alignment);
  assert(!(alignment % kBytesPerWord));
  uintptr_t malloc_addr;
  if (unlikely(pgl_malloc(size + alignment, &malloc_addr))) {
    return 1;
  }
  assert(!(malloc_addr % kBytesPerWord));
  uintptr_t** casted_aligned_pp = S_CAST(uintptr_t**, aligned_pp);
  *casted_aligned_pp = R_CAST(uintptr_t*, RoundDownPow2(malloc_addr + alignment, alignment));
  (*casted_aligned_pp)[-1] = malloc_addr;
  return 0;
}

#ifdef __LP64__
int32_t memequal(const void* m1, const void* m2, uintptr_t byte_ct) {
  const unsigned char* m1_uc = S_CAST(const unsigned char*, m1);
  const unsigned char* m2_uc = S_CAST(const unsigned char*, m2);
  if (byte_ct < 16 + (kBytesPerVec / 2)) {
    if (byte_ct < kBytesPerWord) {
      if (byte_ct < 4) {
        if (byte_ct < 2) {
          return (!byte_ct) || (m1_uc[0] == m2_uc[0]);
        }
        if ((*S_CAST(const uint16_t*, m1)) != (*S_CAST(const uint16_t*, m2))) {
          return 0;
        }
        if ((byte_ct == 3) && (m1_uc[2] != m2_uc[2])) {
          return 0;
        }
        return 1;
      }
      if ((*R_CAST(const uint32_t*, m1_uc)) != (*R_CAST(const uint32_t*, m2_uc))) {
        return 0;
      }
      if (byte_ct > 4) {
        const uintptr_t final_offset = byte_ct - 4;
        if ((*R_CAST(const uint32_t*, &(m1_uc[final_offset]))) != (*R_CAST(const uint32_t*, &(m2_uc[final_offset])))) {
          return 0;
        }
      }
      return 1;
    }
    const uintptr_t* m1_alias = R_CAST(const uintptr_t*, m1_uc);
    const uintptr_t* m2_alias = R_CAST(const uintptr_t*, m2_uc);
    if (m1_alias[0] != m2_alias[0]) {
      return 0;
    }
    if (byte_ct >= 16) {
      if (m1_alias[1] != m2_alias[1]) {
        return 0;
      }
#  ifdef USE_AVX2
      if (byte_ct >= 24) {
        if (m1_alias[2] != m2_alias[2]) {
          return 0;
        }
      }
#  endif
    }
    if (byte_ct % kBytesPerWord) {
      const uintptr_t final_offset = byte_ct - kBytesPerWord;
      if ((*R_CAST(const uintptr_t*, &(m1_uc[final_offset]))) != (*R_CAST(const uintptr_t*, &(m2_uc[final_offset])))) {
        return 0;
      }
    }
    return 1;
  }
  // Don't use VecW since _mm_cmpeq_epi64() not defined until SSE4.1.
  const VecUc* m1_alias = S_CAST(const VecUc*, m1);
  const VecUc* m2_alias = S_CAST(const VecUc*, m2);
  const uintptr_t vec_ct = byte_ct / kBytesPerVec;
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    // tried unrolling this, doesn't make a difference
    const VecUc v1 = vecuc_loadu(&(m1_alias[vidx]));
    const VecUc v2 = vecuc_loadu(&(m2_alias[vidx]));
    if (vecuc_movemask(v1 == v2) != kVec8thUintMax) {
      return 0;
    }
  }
  if (byte_ct % kBytesPerVec) {
    // put this last instead of first, for better behavior when inputs are
    // aligned
    const uintptr_t final_offset = byte_ct - kBytesPerVec;
    const VecUc v1 = vecuc_loadu(&(m1_uc[final_offset]));
    const VecUc v2 = vecuc_loadu(&(m2_uc[final_offset]));
    if (vecuc_movemask(v1 == v2) != kVec8thUintMax) {
      return 0;
    }
  }
  return 1;
}

// clang/gcc memcmp is not that well-optimized for the short strings we usually
// compare.
int32_t Memcmp(const void* m1, const void* m2, uintptr_t byte_ct) {
  const unsigned char* m1_uc = S_CAST(const unsigned char*, m1);
  const unsigned char* m2_uc = S_CAST(const unsigned char*, m2);
  // tried larger crossover threshold, doesn't help
  if (byte_ct < kBytesPerVec) {
    if (byte_ct < kBytesPerWord) {
      if (byte_ct < 4) {
        for (uintptr_t pos = 0; pos != byte_ct; ++pos) {
          const unsigned char ucc1 = m1_uc[pos];
          const unsigned char ucc2 = m2_uc[pos];
          if (ucc1 != ucc2) {
            return (ucc1 < ucc2)? -1 : 1;
          }
        }
        return 0;
      }
      uint32_t m1_u32 = *S_CAST(const uint32_t*, m1);
      uint32_t m2_u32 = *S_CAST(const uint32_t*, m2);
      if (m1_u32 != m2_u32) {
        return (__builtin_bswap32(m1_u32) < __builtin_bswap32(m2_u32))? -1 : 1;
      }
      if (byte_ct > 4) {
        const uintptr_t final_offset = byte_ct - 4;
        m1_u32 = *R_CAST(const uint32_t*, &(m1_uc[final_offset]));
        m2_u32 = *R_CAST(const uint32_t*, &(m2_uc[final_offset]));
        if (m1_u32 != m2_u32) {
          return (__builtin_bswap32(m1_u32) < __builtin_bswap32(m2_u32))? -1 : 1;
        }
      }
      return 0;
    }
    const uintptr_t* m1_alias = R_CAST(const uintptr_t*, m1_uc);
    const uintptr_t* m2_alias = R_CAST(const uintptr_t*, m2_uc);
    uintptr_t m1_word = m1_alias[0];
    uintptr_t m2_word = m2_alias[0];
    if (m1_word != m2_word) {
      return (__builtin_bswap64(m1_word) < __builtin_bswap64(m2_word))? -1 : 1;
    }
#  ifdef USE_AVX2
    if (byte_ct >= 16) {
      m1_word = m1_alias[1];
      m2_word = m2_alias[1];
      if (m1_word != m2_word) {
        return (__builtin_bswap64(m1_word) < __builtin_bswap64(m2_word))? -1 : 1;
      }
      if (byte_ct >= 24) {
        m1_word = m1_alias[2];
        m2_word = m2_alias[2];
        if (m1_word != m2_word) {
          return (__builtin_bswap64(m1_word) < __builtin_bswap64(m2_word))? -1 : 1;
        }
      }
    }
#  endif
    if (byte_ct % kBytesPerWord) {
      const uintptr_t final_offset = byte_ct - kBytesPerWord;
      m1_word = *R_CAST(const uintptr_t*, &(m1_uc[final_offset]));
      m2_word = *R_CAST(const uintptr_t*, &(m2_uc[final_offset]));
      if (m1_word != m2_word) {
        return (__builtin_bswap64(m1_word) < __builtin_bswap64(m2_word))? -1 : 1;
      }
    }
    return 0;
  }
  const VecUc* m1_alias = S_CAST(const VecUc*, m1);
  const VecUc* m2_alias = S_CAST(const VecUc*, m2);
  const uintptr_t fullvec_ct = byte_ct / kBytesPerVec;
  // uh, clang/LLVM -O2 optimizes this better when comparison is != instead of
  // <?  ugh, time to change all of the for loops...
  // (and yes, both -O3 configurations generate worse code here)
  // at least for loop is better than do-while loop even when 1 iteration is
  // guaranteed...
  for (uintptr_t vidx = 0; vidx != fullvec_ct; ++vidx) {
    const VecUc v1 = vecuc_loadu(&(m1_alias[vidx]));
    const VecUc v2 = vecuc_loadu(&(m2_alias[vidx]));
    // is this even worthwhile now in non-AVX2 case?
    const uint32_t movemask_result = vecuc_movemask(v1 == v2);
    if (movemask_result != kVec8thUintMax) {
      const uintptr_t diff_pos = vidx * kBytesPerVec + ctzu32(~movemask_result);
      return (m1_uc[diff_pos] < m2_uc[diff_pos])? -1 : 1;
    }
  }
  if (byte_ct % kBytesPerVec) {
    const uintptr_t final_offset = byte_ct - kBytesPerVec;
    const VecUc v1 = vecuc_loadu(&(m1_uc[final_offset]));
    const VecUc v2 = vecuc_loadu(&(m2_uc[final_offset]));
    const uint32_t movemask_result = vecuc_movemask(v1 == v2);
    if (movemask_result != kVec8thUintMax) {
      const uintptr_t diff_pos = final_offset + ctzu32(~movemask_result);
      return (m1_uc[diff_pos] < m2_uc[diff_pos])? -1 : 1;
    }
  }
  return 0;
}
#endif

#ifdef __cplusplus
}  // namespace plink2
#endif
