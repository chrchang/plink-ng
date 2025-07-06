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


#include "plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

const char kErrprintfFopen[] = "Error: Failed to open %s : %s.\n";
const char kErrprintfFread[] = "Error: %s read failure: %s.\n";
const char kErrprintfRewind[] = "Error: %s could not be scanned twice. (Process-substitution/named-pipe input is not permitted in this use case.)\n";
const char kErrstrNomem[] = "Error: Out of memory.  The --memory flag may be helpful.\n";
const char kErrstrWrite[] = "Error: File write failure: %s.\n";
const char kErrprintfDecompress[] = "Error: %s decompression failure: %s.\n";

uint64_t g_failed_alloc_attempt_size = 0;

#if (((__GNUC__ == 4) && (__GNUC_MINOR__ < 7)) || (__GNUC__ >= 11)) && !defined(__APPLE__)
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
  if (unlikely(ScanUintCappedFinish(str_iter, bound, I32ToU32(valp)))) {
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
#endif // defined(USE_SSE2) && !defined(NO_UNALIGNED)

const uint16_t kDigitPair[] = {
  0x3030, 0x3130, 0x3230, 0x3330, 0x3430, 0x3530, 0x3630, 0x3730, 0x3830, 0x3930,
  0x3031, 0x3131, 0x3231, 0x3331, 0x3431, 0x3531, 0x3631, 0x3731, 0x3831, 0x3931,
  0x3032, 0x3132, 0x3232, 0x3332, 0x3432, 0x3532, 0x3632, 0x3732, 0x3832, 0x3932,
  0x3033, 0x3133, 0x3233, 0x3333, 0x3433, 0x3533, 0x3633, 0x3733, 0x3833, 0x3933,
  0x3034, 0x3134, 0x3234, 0x3334, 0x3434, 0x3534, 0x3634, 0x3734, 0x3834, 0x3934,
  0x3035, 0x3135, 0x3235, 0x3335, 0x3435, 0x3535, 0x3635, 0x3735, 0x3835, 0x3935,
  0x3036, 0x3136, 0x3236, 0x3336, 0x3436, 0x3536, 0x3636, 0x3736, 0x3836, 0x3936,
  0x3037, 0x3137, 0x3237, 0x3337, 0x3437, 0x3537, 0x3637, 0x3737, 0x3837, 0x3937,
  0x3038, 0x3138, 0x3238, 0x3338, 0x3438, 0x3538, 0x3638, 0x3738, 0x3838, 0x3938,
  0x3039, 0x3139, 0x3239, 0x3339, 0x3439, 0x3539, 0x3639, 0x3739, 0x3839, 0x3939};

char* u32toa(uint32_t uii, char* start) {
  // Memory-efficient fast integer writer.  (You can do a bit better sometimes
  // by using a larger lookup table, but on average I doubt that pays off.)
  // Returns a pointer to the end of the integer (not null-terminated).
  //
  // Nearly identical to 'branchlut' from
  // https://github.com/miloyip/itoa-benchmark , except that the hardcoded
  // binary search is more balanced (start by comparing 6+ digits vs. <6,
  // instead of 9+ digits vs. <8).  This tends to be slightly better unless the
  // integers are almost uniformly distributed over [0, 2^32).
  //
  // Todo: compare against an_itoa in https://github.com/appnexus/acf/ .
  //
  // (Making the first comparison 7+ digits vs. <7 would seem to make sense,
  // but it seems to benchmark slightly worse on my Mac?)
  //
  // (Since we want to execute different code depending on the number of
  // digits, the UintSlen() approach doesn't pay off.)
  uint32_t quotient;
  if (uii < 100000) {
    if (uii < 100) {
      if (uii >= 10) {
        goto u32toa_just2;
      }
      *start++ = '0' + uii;
      return start;
    }
    if (uii < 10000) {
      if (uii >= 1000) {
        goto u32toa_just4;
      }
      quotient = uii / 100;
      *start++ = '0' + quotient;
      goto u32toa_2left;
    }
    quotient = uii / 10000;
    *start++ = '0' + quotient;
    goto u32toa_4left;
  }
  if (uii < 100000000) {
    if (uii < 1000000) {
      goto u32toa_just6;
    }
    if (uii >= 10000000) {
      goto u32toa_just8;
    }
    quotient = uii / 1000000;
    *start++ = '0' + quotient;
    goto u32toa_6left;
  }
  quotient = uii / 100000000;
  if (uii < 1000000000) {
    *start++ = '0' + quotient;
  } else {
    start = memcpya_k(start, &(kDigitPair[quotient]), 2);
  }
  uii -= quotient * 100000000;
 u32toa_just8:
  quotient = uii / 1000000;
  start = memcpya_k(start, &(kDigitPair[quotient]), 2);
 u32toa_6left:
  uii -= quotient * 1000000;
 u32toa_just6:
  quotient = uii / 10000;
  start = memcpya_k(start, &(kDigitPair[quotient]), 2);
 u32toa_4left:
  uii -= quotient * 10000;
 u32toa_just4:
  quotient = uii / 100;
  start = memcpya_k(start, &(kDigitPair[quotient]), 2);
 u32toa_2left:
  uii -= quotient * 100;
 u32toa_just2:
  return memcpya_k(start, &(kDigitPair[uii]), 2);
}

char* i64toa(int64_t llii, char* start) {
  uint64_t ullii = llii;
  uint64_t top_digits;
  uint32_t bottom_eight;
  uint32_t middle_eight;
  if (llii < 0) {
    *start++ = '-';
    ullii = -ullii;
  }
  if (ullii <= 0xffffffffLLU) {
    return u32toa(S_CAST(uint32_t, ullii), start);
  }
  top_digits = ullii / 100000000;
  bottom_eight = S_CAST(uint32_t, ullii - (top_digits * 100000000));
  if (top_digits <= 0xffffffffLLU) {
    start = u32toa(S_CAST(uint32_t, top_digits), start);
    return uitoa_z8(bottom_eight, start);
  }
  ullii = top_digits / 100000000;
  middle_eight = S_CAST(uint32_t, top_digits - (ullii * 100000000));
  start = u32toa(S_CAST(uint32_t, ullii), start);
  start = uitoa_z8(middle_eight, start);
  return uitoa_z8(bottom_eight, start);
}

#if defined(USE_SSE2) && !defined(NO_UNALIGNED)
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
    const uint32_t eq_result = vecw_movemask(v1 == v2);
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
