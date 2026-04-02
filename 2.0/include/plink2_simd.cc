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


#include "plink2_simd.h"

#ifdef __cplusplus
namespace plink2 {
#endif

#if defined(USE_SSE2) && !defined(NO_UNALIGNED)
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
    if (!vecucs_are_equal(v1, v2)) {
      return 0;
    }
  }
  if (byte_ct % kBytesPerVec) {
    // put this last instead of first, for better behavior when inputs are
    // aligned
    const uintptr_t final_offset = byte_ct - kBytesPerVec;
    const VecUc v1 = vecuc_loadu(&(m1_uc[final_offset]));
    const VecUc v2 = vecuc_loadu(&(m2_uc[final_offset]));
    if (!vecucs_are_equal(v1, v2)) {
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
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    // is this even worthwhile now in non-AVX2 case?
    const uint32_t movemask_result = vecuc_movemask(v1 == v2);
    if (movemask_result != kVec8thUintMax) {
      // todo: check if this is faster to do outside of vector-space
      const uintptr_t diff_pos = vidx * kBytesPerVec + ctzu32(~movemask_result);
      return (m1_uc[diff_pos] < m2_uc[diff_pos])? -1 : 1;
    }
#  else
    const uint64_t eq_nybbles = arm_shrn4_uc(v1 == v2);
    if (eq_nybbles != UINT64_MAX) {
      const uintptr_t diff_pos = vidx * kBytesPerVec + ctzw(~eq_nybbles) / 4;
      return (m1_uc[diff_pos] < m2_uc[diff_pos])? -1 : 1;
    }
#  endif
  }
  if (byte_ct % kBytesPerVec) {
    const uintptr_t final_offset = byte_ct - kBytesPerVec;
    const VecUc v1 = vecuc_loadu(&(m1_uc[final_offset]));
    const VecUc v2 = vecuc_loadu(&(m2_uc[final_offset]));
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    const uint32_t movemask_result = vecuc_movemask(v1 == v2);
    if (movemask_result != kVec8thUintMax) {
      const uintptr_t diff_pos = final_offset + ctzu32(~movemask_result);
      return (m1_uc[diff_pos] < m2_uc[diff_pos])? -1 : 1;
    }
#  else
    const uint64_t eq_nybbles = arm_shrn4_uc(v1 == v2);
    if (eq_nybbles != UINT64_MAX) {
      const uintptr_t diff_pos = final_offset + ctzw(~eq_nybbles) / 4;
      return (m1_uc[diff_pos] < m2_uc[diff_pos])? -1 : 1;
    }
#  endif
  }
  return 0;
}

uintptr_t FirstUnequal4(const void* arr1, const void* arr2, uintptr_t nbytes) {
  // Similar to memequal().
  if (nbytes < kBytesPerVec) {
    if (nbytes < kBytesPerWord) {
      uint32_t xor_result = (*S_CAST(const uint32_t*, arr1)) ^ (*S_CAST(const uint32_t*, arr2));
      if (xor_result) {
        return ctzu32(xor_result) / CHAR_BIT;
      }
      if (nbytes > 4) {
        const uintptr_t final_offset = nbytes - 4;
        const char* s1 = S_CAST(const char*, arr1);
        const char* s2 = S_CAST(const char*, arr2);
        xor_result = (*R_CAST(const uint32_t*, &(s1[final_offset]))) ^ (*R_CAST(const uint32_t*, &(s2[final_offset])));
        if (xor_result) {
          return final_offset + ctzu32(xor_result) / CHAR_BIT;
        }
      }
      return nbytes;
    }
    const uintptr_t* arr1_alias = S_CAST(const uintptr_t*, arr1);
    const uintptr_t* arr2_alias = S_CAST(const uintptr_t*, arr2);
    const uintptr_t word_ct = nbytes / kBytesPerWord;
    for (uint32_t widx = 0; widx != word_ct; ++widx) {
      const uintptr_t xor_result = arr1_alias[widx] ^ arr2_alias[widx];
      if (xor_result) {
        return widx * kBytesPerWord + ctzw(xor_result) / CHAR_BIT;
      }
    }
    if (nbytes % kBytesPerWord) {
      const uintptr_t final_offset = nbytes - kBytesPerWord;
      const char* s1 = S_CAST(const char*, arr1);
      const char* s2 = S_CAST(const char*, arr2);
      const uintptr_t xor_result = (*R_CAST(const uintptr_t*, &(s1[final_offset]))) ^ (*R_CAST(const uintptr_t*, &(s2[final_offset])));
      if (xor_result) {
        return final_offset + ctzw(xor_result) / CHAR_BIT;
      }
    }
    return nbytes;
  }
  const VecUc* arr1_alias = S_CAST(const VecUc*, arr1);
  const VecUc* arr2_alias = S_CAST(const VecUc*, arr2);
  const uintptr_t vec_ct = nbytes / kBytesPerVec;
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    const VecUc v1 = vecuc_loadu(&(arr1_alias[vidx]));
    const VecUc v2 = vecuc_loadu(&(arr2_alias[vidx]));
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    const uint32_t eq_result = vecuc_movemask(v1 == v2);
    if (eq_result != kVec8thUintMax) {
      return vidx * kBytesPerVec + ctzu32(~eq_result);
    }
#  else
    const uint64_t eq_nybbles = arm_shrn4_uc(v1 == v2);
    if (eq_nybbles != UINT64_MAX) {
      return vidx * kBytesPerVec + ctzw(~eq_nybbles) / 4;
    }
#  endif
  }
  if (nbytes % kBytesPerVec) {
    const uintptr_t final_offset = nbytes - kBytesPerVec;
    const char* s1 = S_CAST(const char*, arr1);
    const char* s2 = S_CAST(const char*, arr2);
    // bugfix (5 Jul 2025): this must be VecUc on ARM, not VecW
    // (disturbing that arm_shrn4_uc() call compiled at all...)
    const VecUc v1 = vecuc_loadu(&(s1[final_offset]));
    const VecUc v2 = vecuc_loadu(&(s2[final_offset]));
#  ifndef SIMDE_ARM_NEON_A32V8_NATIVE
    const uint32_t eq_result = vecuc_movemask(v1 == v2);
    if (eq_result != kVec8thUintMax) {
      return final_offset + ctzu32(~eq_result);
    }
#  else
    const uint64_t eq_nybbles = arm_shrn4_uc(v1 == v2);
    if (eq_nybbles != UINT64_MAX) {
      return final_offset + ctzw(~eq_nybbles) / 4;
    }
#  endif
  }
  return nbytes;
}
#else // !(defined(USE_SSE2) && !defined(NO_UNALIGNED))
uintptr_t FirstUnequalW(const void* arr1, const void* arr2, uintptr_t nbytes) {
  const unsigned char* arr1b = S_CAST(const unsigned char*, arr1);
  const unsigned char* arr2b = S_CAST(const unsigned char*, arr2);
  const uintptr_t word_ct = nbytes / kBytesPerWord;
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    uintptr_t arr1_word;
    uintptr_t arr2_word;
    CopyFromUnalignedOffsetW(&arr1_word, arr1b, widx);
    CopyFromUnalignedOffsetW(&arr2_word, arr2b, widx);
    const uintptr_t xor_result = arr1_word ^ arr2_word;
    if (xor_result) {
      return widx * kBytesPerWord + ctzw(xor_result) / CHAR_BIT;
    }
  }
  if (nbytes % kBytesPerWord) {
    const uintptr_t final_offset = nbytes - kBytesPerWord;
    uintptr_t arr1_word;
    uintptr_t arr2_word;
    CopyFromUnalignedW(&arr1_word, &(arr1b[final_offset]));
    CopyFromUnalignedW(&arr2_word, &(arr2b[final_offset]));
    const uintptr_t xor_result = arr1_word ^ arr2_word;
    if (xor_result) {
      return final_offset + ctzw(xor_result) / CHAR_BIT;
    }
  }
  return nbytes;
}
#endif

#ifdef __LP64__
uintptr_t CountVintsNonempty(const unsigned char* buf, const unsigned char* buf_end) {
  const uintptr_t starting_addr = R_CAST(uintptr_t, buf);
  const VecUc* buf_viter = R_CAST(const VecUc*, RoundDownPow2(starting_addr, kBytesPerVec));
  const uintptr_t ending_addr = R_CAST(uintptr_t, buf_end);
  const VecUc* buf_vlast = R_CAST(const VecUc*, RoundDownPow2(ending_addr - 1, kBytesPerVec));
  const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, buf_viter);
  // todo: better ARM implementation
  Vec8thUint vint_ends = (UINT32_MAX << leading_byte_ct) & (~vecuc_movemask(*buf_viter));
  uintptr_t total = 0;
  while (buf_viter != buf_vlast) {
    total += PopcountVec8thUint(vint_ends);
    ++buf_viter;
    vint_ends = ~vecuc_movemask(*buf_viter);
  }
  const uint32_t trailing_byte_ct = ending_addr - R_CAST(uintptr_t, buf_vlast);
  vint_ends &= (k1LU << trailing_byte_ct) - 1;
  total += PopcountVec8thUint(vint_ends);
  return total;
}
#else
uintptr_t CountVints(const unsigned char* buf, const unsigned char* buf_end) {
  // Could check one word at a time.
  const uintptr_t len = buf_end - buf;
  uintptr_t inv_result = 0;
  for (uintptr_t ulii = 0; ulii != len; ++ulii) {
    inv_result += buf[ulii] >> 7;
  }
  return len - inv_result;
}
#endif

#ifdef __cplusplus
}  // namespace plink2
#endif
