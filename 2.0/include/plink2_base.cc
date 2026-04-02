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
    start = memcpya_k2(start, &(kDigitPair[quotient]));
  }
  uii -= quotient * 100000000;
 u32toa_just8:
  quotient = uii / 1000000;
  start = memcpya_k2(start, &(kDigitPair[quotient]));
 u32toa_6left:
  uii -= quotient * 1000000;
 u32toa_just6:
  quotient = uii / 10000;
  start = memcpya_k2(start, &(kDigitPair[quotient]));
 u32toa_4left:
  uii -= quotient * 10000;
 u32toa_just4:
  quotient = uii / 100;
  start = memcpya_k2(start, &(kDigitPair[quotient]));
 u32toa_2left:
  uii -= quotient * 100;
 u32toa_just2:
  return memcpya_k2(start, &(kDigitPair[uii]));
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

uint32_t MaxElementU32(const uint32_t* u32arr, uintptr_t entry_ct) {
  // confirmed that *std::max_element() has horrible performance on at least
  // macOS
  // also confirmed that macOS compiler autovectorizes this; main loop
  // processes 4 vectors at a time
  uint32_t result = u32arr[0];
  for (uintptr_t entry_idx = 1; entry_idx != entry_ct; ++entry_idx) {
    const uint32_t cur_element = u32arr[entry_idx];
    if (cur_element > result) {
      result = cur_element;
    }
  }
  return result;
}

double MaxElementD(const double* darr, uintptr_t entry_ct) {
  double result = darr[0];
  for (uintptr_t entry_idx = 1; entry_idx != entry_ct; ++entry_idx) {
    const double cur_element = darr[entry_idx];
    if (cur_element > result) {
      result = cur_element;
    }
  }
  return result;
}

double MinElementD(const double* darr, uintptr_t entry_ct) {
  double result = darr[0];
  for (uintptr_t entry_idx = 1; entry_idx != entry_ct; ++entry_idx) {
    const double cur_element = darr[entry_idx];
    if (cur_element < result) {
      result = cur_element;
    }
  }
  return result;
}


int32_t u32cmp(const void* aa, const void* bb) {
  const uint32_t uaa = *S_CAST(const uint32_t*, aa);
  const uint32_t ubb = *S_CAST(const uint32_t*, bb);
  if (uaa < ubb) {
    return -1;
  }
  return (uaa > ubb);
}

#ifdef __cplusplus
}  // namespace plink2
#endif
