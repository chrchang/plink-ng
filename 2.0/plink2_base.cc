// This library is part of PLINK 2, copyright (C) 2005-2018 Shaun Purcell,
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
    } else if (*valp != -5) {
      return 1;
    }
    *valp = ctou32(*str_iter++) - 48;
    if (unlikely(ctou32(*valp) >= 10)) {
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
  while (1) {
    const uint32_t cur_digit = ctou32(*str_iter++) - 48;
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
  while (1) {
    const uint32_t cur_digit = ctou32(*str_iter++) - 48;
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
    } else if (val != 0xfffffffbU) {
      return 1;
    }
    val = ctou32(*str_iter++) - 48;
    if (unlikely(val >= 10)) {
      return 1;
    }
  }
  while (1) {
    const uint32_t cur_digit = ctou32(*str_iter++) - 48;
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

void SetAllBits(uintptr_t ct, uintptr_t* bitarr) {
  // leaves bits beyond the end unset
  // ok for ct == 0
  uintptr_t quotient = ct / kBitsPerWord;
  uintptr_t remainder = ct % kBitsPerWord;
  SetAllWArr(quotient, bitarr);
  if (remainder) {
    bitarr[quotient] = (k1LU << remainder) - k1LU;
  }
}

void BitvecAnd(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec AND arg_bitvec
#ifdef __LP64__
  VecW* main_bitvvec_iter = R_CAST(VecW*, main_bitvec);
  const VecW* arg_bitvvec_iter = R_CAST(const VecW*, arg_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    *main_bitvvec_iter++ &= *arg_bitvvec_iter++;
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter++ &= *arg_bitvvec_iter++;
    *main_bitvvec_iter++ &= *arg_bitvvec_iter++;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter++ &= *arg_bitvvec_iter++;
    *main_bitvvec_iter++ &= *arg_bitvvec_iter++;
    *main_bitvvec_iter++ &= *arg_bitvvec_iter++;
    *main_bitvvec_iter++ &= *arg_bitvvec_iter++;
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] &= arg_bitvec[base_idx];
    main_bitvec[base_idx + 1] &= arg_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] &= arg_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    main_bitvec[widx] &= arg_bitvec[widx];
  }
#endif
}

void BitvecAndNot(const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec ANDNOT exclude_bitvec
  // note that this is the reverse of the _mm_andnot() operand order
#ifdef __LP64__
  VecW* main_bitvvec_iter = R_CAST(VecW*, main_bitvec);
  const VecW* exclude_bitvvec_iter = R_CAST(const VecW*, exclude_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    *main_bitvvec_iter++ &= ~(*exclude_bitvvec_iter++);
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter++ &= ~(*exclude_bitvvec_iter++);
    *main_bitvvec_iter++ &= ~(*exclude_bitvvec_iter++);
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter++ &= ~(*exclude_bitvvec_iter++);
    *main_bitvvec_iter++ &= ~(*exclude_bitvvec_iter++);
    *main_bitvvec_iter++ &= ~(*exclude_bitvvec_iter++);
    *main_bitvvec_iter++ &= ~(*exclude_bitvvec_iter++);
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] &= ~exclude_bitvec[base_idx];
    main_bitvec[base_idx + 1] &= ~exclude_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] &= ~exclude_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    main_bitvec[widx] &= ~exclude_bitvec[widx];
  }
#endif
}

void BitvecInvert(uintptr_t word_ct, uintptr_t* main_bitvec) {
#ifdef __LP64__
  VecW* main_bitvvec_iter = R_CAST(VecW*, main_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  const VecW all1 = VCONST_W(~k0LU);
  if (full_vec_ct & 1) {
    *main_bitvvec_iter++ ^= all1;
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter++ ^= all1;
    *main_bitvvec_iter++ ^= all1;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter++ ^= all1;
    *main_bitvvec_iter++ ^= all1;
    *main_bitvvec_iter++ ^= all1;
    *main_bitvvec_iter++ ^= all1;
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] ^= ~k0LU;
    main_bitvec[base_idx + 1] ^= ~k0LU;
  }
#  endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] ^= ~k0LU;
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    main_bitvec[widx] ^= ~k0LU;
  }
#endif
}

uintptr_t AdvTo1Bit(const uintptr_t* bitarr, uintptr_t loc) {
  const uintptr_t* bitarr_iter = &(bitarr[loc / kBitsPerWord]);
  uintptr_t ulii = (*bitarr_iter) >> (loc % kBitsPerWord);
  if (ulii) {
    return loc + ctzw(ulii);
  }
  do {
    ulii = *(++bitarr_iter);
  } while (!ulii);
  return S_CAST(uintptr_t, bitarr_iter - bitarr) * kBitsPerWord + ctzw(ulii);
}

uintptr_t AdvTo0Bit(const uintptr_t* bitarr, uintptr_t loc) {
  const uintptr_t* bitarr_iter = &(bitarr[loc / kBitsPerWord]);
  uintptr_t ulii = (~(*bitarr_iter)) >> (loc % kBitsPerWord);
  if (ulii) {
    return loc + ctzw(ulii);
  }
  do {
    ulii = *(++bitarr_iter);
  } while (ulii == ~k0LU);
  return S_CAST(uintptr_t, bitarr_iter - bitarr) * kBitsPerWord + ctzw(~ulii);
}

/*
uintptr_t NextNonmissingUnsafe(const uintptr_t* genoarr, uintptr_t loc) {
  const uintptr_t* genoarr_iter = &(genoarr[loc / kBitsPerWordD2]);
  uintptr_t ulii = (~(*genoarr_iter)) >> (2 * (loc % kBitsPerWordD2));
  if (ulii) {
    return loc + (ctzw(ulii) / 2);
  }
  do {
    ulii = *(++genoarr_iter);
  } while (ulii == ~k0LU);
  return S_CAST(uintptr_t, genoarr_iter - genoarr) * kBitsPerWordD2 + (ctzw(~ulii) / 2);
}
*/

uint32_t AdvBoundedTo1Bit(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil) {
  // safe version.
  const uintptr_t* bitarr_iter = &(bitarr[loc / kBitsPerWord]);
  uintptr_t ulii = (*bitarr_iter) >> (loc % kBitsPerWord);
  uint32_t rval;
  if (ulii) {
    rval = loc + ctzw(ulii);
    return MINV(rval, ceil);
  }
  const uintptr_t* bitarr_last = &(bitarr[(ceil - 1) / kBitsPerWord]);
  do {
    if (bitarr_iter >= bitarr_last) {
      return ceil;
    }
    ulii = *(++bitarr_iter);
  } while (!ulii);
  rval = S_CAST(uintptr_t, bitarr_iter - bitarr) * kBitsPerWord + ctzw(ulii);
  return MINV(rval, ceil);
}

uint32_t FindLast1BitBefore(const uintptr_t* bitarr, uint32_t loc) {
  // unlike the next_{un}set family, this always returns a STRICTLY earlier
  // position
  const uintptr_t* bitarr_iter = &(bitarr[loc / kBitsPerWord]);
  const uint32_t remainder = loc % kBitsPerWord;
  uintptr_t ulii;
  if (remainder) {
    ulii = bzhi(*bitarr_iter, remainder);
    if (ulii) {
      return loc - remainder + bsrw(ulii);
    }
  }
  do {
    ulii = *(--bitarr_iter);
  } while (!ulii);
  return S_CAST(uintptr_t, bitarr_iter - bitarr) * kBitsPerWord + bsrw(ulii);
}

#ifdef USE_AVX2
// void CopyBitarrSubsetEx(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t bit_idx_start, uint32_t bit_idx_end, uintptr_t* __restrict output_bitarr) {
void CopyBitarrSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t bit_idx_end, uintptr_t* __restrict output_bitarr) {
  const uint32_t bit_idx_end_lowbits = bit_idx_end % kBitsPerWord;
  uintptr_t* output_bitarr_iter = output_bitarr;
  uintptr_t* output_bitarr_last = &(output_bitarr[bit_idx_end / kBitsPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = 0;
  while ((output_bitarr_iter != output_bitarr_last) || (write_idx_lowbits != bit_idx_end_lowbits)) {
    uintptr_t cur_mask_word;
    // sparse subset_mask optimization
    // guaranteed to terminate since there's at least one more set bit
    do {
      cur_mask_word = subset_mask[++read_widx];
    } while (!cur_mask_word);
    uintptr_t extracted_bits = raw_bitarr[read_widx];
    uint32_t set_bit_ct = kBitsPerWord;
    if (cur_mask_word != ~k0LU) {
      extracted_bits = _pext_u64(extracted_bits, cur_mask_word);
      set_bit_ct = PopcountWord(cur_mask_word);
    }
    cur_output_word |= extracted_bits << write_idx_lowbits;
    const uint32_t new_write_idx_lowbits = write_idx_lowbits + set_bit_ct;
    if (new_write_idx_lowbits >= kBitsPerWord) {
      *output_bitarr_iter++ = cur_output_word;
      // ...and these are the bits that fell off
      // bugfix: unsafe to right-shift 64
      if (write_idx_lowbits) {
        cur_output_word = extracted_bits >> (kBitsPerWord - write_idx_lowbits);
      } else {
        cur_output_word = 0;
      }
    }
    write_idx_lowbits = new_write_idx_lowbits % kBitsPerWord;
  }
  if (write_idx_lowbits) {
    *output_bitarr_iter = cur_output_word;
  }
}

uintptr_t PopcountVecsAvx2(const VecW* bit_vvec, uintptr_t vec_ct) {
  // See popcnt_avx2() in libpopcnt.
  VecW cnt = vecw_setzero();
  VecW ones = vecw_setzero();
  VecW twos = vecw_setzero();
  VecW fours = vecw_setzero();
  VecW eights = vecw_setzero();
  for (uintptr_t vec_idx = 0; vec_idx < vec_ct; vec_idx += 16) {
    VecW twos_a = Csa256(bit_vvec[vec_idx + 0], bit_vvec[vec_idx + 1], &ones);
    VecW twos_b = Csa256(bit_vvec[vec_idx + 2], bit_vvec[vec_idx + 3], &ones);
    VecW fours_a = Csa256(twos_a, twos_b, &twos);

    twos_a = Csa256(bit_vvec[vec_idx + 4], bit_vvec[vec_idx + 5], &ones);
    twos_b = Csa256(bit_vvec[vec_idx + 6], bit_vvec[vec_idx + 7], &ones);
    VecW fours_b = Csa256(twos_a, twos_b, &twos);
    const VecW eights_a = Csa256(fours_a, fours_b, &fours);

    twos_a = Csa256(bit_vvec[vec_idx + 8], bit_vvec[vec_idx + 9], &ones);
    twos_b = Csa256(bit_vvec[vec_idx + 10], bit_vvec[vec_idx + 11], &ones);
    fours_a = Csa256(twos_a, twos_b, &twos);

    twos_a = Csa256(bit_vvec[vec_idx + 12], bit_vvec[vec_idx + 13], &ones);
    twos_b = Csa256(bit_vvec[vec_idx + 14], bit_vvec[vec_idx + 15], &ones);
    fours_b = Csa256(twos_a, twos_b, &twos);
    const VecW eights_b = Csa256(fours_a, fours_b, &fours);
    const VecW sixteens = Csa256(eights_a, eights_b, &eights);
    cnt = cnt + PopcountVecAvx2(sixteens);
  }
  cnt = vecw_slli(cnt, 4);
  cnt = cnt + vecw_slli(PopcountVecAvx2(eights), 3);
  cnt = cnt + vecw_slli(PopcountVecAvx2(fours), 2);
  cnt = cnt + vecw_slli(PopcountVecAvx2(twos), 1);
  cnt = cnt + PopcountVecAvx2(ones);
  return Hsum64(cnt);
}

void ExpandBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, uint32_t word_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t leading_byte_ct = 1 + (expand_sizex_m1 % kBitsPerWord) / CHAR_BIT;
  uintptr_t compact_word = SubwordLoad(compact_bitarr, leading_byte_ct) >> read_start_bit;
  const uintptr_t* compact_bitarr_iter = R_CAST(const uintptr_t*, &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]));
  uint32_t compact_idx_lowbits = read_start_bit + CHAR_BIT * (sizeof(intptr_t) - leading_byte_ct);
  for (uint32_t widx = 0; widx < word_ct; ++widx) {
    const uintptr_t mask_word = expand_mask[widx];
    uintptr_t write_word = 0;
    if (mask_word) {
      const uint32_t mask_set_ct = PopcountWord(mask_word);
      uint32_t next_compact_idx_lowbits = compact_idx_lowbits + mask_set_ct;
      if (next_compact_idx_lowbits <= kBitsPerWord) {
        write_word = _pdep_u64(compact_word, mask_word);
        if (mask_set_ct != kBitsPerWord) {
          compact_word = compact_word >> mask_set_ct;
        } else {
          // avoid nasal demons
          compact_word = 0;
        }
      } else {
#ifdef __arm__
#  error "Unaligned accesses in ExpandBytearr()."
#endif
        uintptr_t next_compact_word = *compact_bitarr_iter++;
        next_compact_idx_lowbits -= kBitsPerWord;
        compact_word |= next_compact_word << (kBitsPerWord - compact_idx_lowbits);
        write_word = _pdep_u64(compact_word, mask_word);
        if (next_compact_idx_lowbits != kBitsPerWord) {
          compact_word = next_compact_word >> next_compact_idx_lowbits;
        } else {
          compact_word = 0;
        }
      }
      compact_idx_lowbits = next_compact_idx_lowbits;
    }
    target[widx] = write_word;
  }
}

void ExpandThenSubsetBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, const uintptr_t* __restrict subset_mask, uint32_t expand_size, uint32_t subset_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t leading_byte_ct = 1 + (expand_sizex_m1 % kBitsPerWord) / CHAR_BIT;
  uintptr_t compact_word = SubwordLoad(compact_bitarr, leading_byte_ct) >> read_start_bit;
  const uintptr_t* compact_bitarr_alias = R_CAST(const uintptr_t*, &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]));
  uint32_t compact_widx = UINT32_MAX;  // deliberate overflow
  uint32_t compact_idx_lowbits = read_start_bit + CHAR_BIT * (sizeof(uintptr_t) - leading_byte_ct);
  const uint32_t subset_size_lowbits = subset_size % kBitsPerWord;
  uintptr_t* target_iter = target;
  uintptr_t* target_last = &(target[subset_size / kBitsPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = 0;

  // bugfix (5 Feb 2018): missed a case in sparse subset_mask optimization
  uint32_t expand_bit_ct_skip = 0;
  while ((target_iter != target_last) || (write_idx_lowbits != subset_size_lowbits)) {
    uintptr_t expand_word;
    uintptr_t subset_word;
    uint32_t expand_bit_ct;
    while (1) {
      ++read_widx;
      expand_word = expand_mask[read_widx];
      subset_word = subset_mask[read_widx];
      expand_bit_ct = PopcountWord(expand_word);
      if (subset_word) {
        break;
      }
      expand_bit_ct_skip += expand_bit_ct;
    }
    uintptr_t extracted_bits = 0;
    const uint32_t set_bit_ct = PopcountWord(subset_word);
    if (expand_word & subset_word) {
      // lazy load
      compact_idx_lowbits += expand_bit_ct_skip;
      if (compact_idx_lowbits >= kBitsPerWord) {
        compact_widx += compact_idx_lowbits / kBitsPerWord;
        compact_idx_lowbits = compact_idx_lowbits % kBitsPerWord;
#ifdef __arm__
#  error "Unaligned accesses in ExpandThenSubsetBytearr()."
#endif
        compact_word = compact_bitarr_alias[compact_widx] >> compact_idx_lowbits;
      } else {
        compact_word = compact_word >> expand_bit_ct_skip;
      }
      uint32_t next_compact_idx_lowbits = compact_idx_lowbits + expand_bit_ct;
      uintptr_t expanded_bits;
      if (next_compact_idx_lowbits <= kBitsPerWord) {
        expanded_bits = _pdep_u64(compact_word, expand_word);
        if (expand_bit_ct != kBitsPerWord) {
          compact_word = compact_word >> expand_bit_ct;
        }
      } else {
        uintptr_t next_compact_word = compact_bitarr_alias[++compact_widx];
        next_compact_idx_lowbits -= kBitsPerWord;
        compact_word |= next_compact_word << (kBitsPerWord - compact_idx_lowbits);
        expanded_bits = _pdep_u64(compact_word, expand_word);
        if (next_compact_idx_lowbits != kBitsPerWord) {
          compact_word = next_compact_word >> next_compact_idx_lowbits;
        }
      }
      extracted_bits = _pext_u64(expanded_bits, subset_word);
      compact_idx_lowbits = next_compact_idx_lowbits;
      cur_output_word |= extracted_bits << write_idx_lowbits;
      expand_bit_ct_skip = 0;
    } else {
      expand_bit_ct_skip += expand_bit_ct;
    }
    const uint32_t new_write_idx_lowbits = write_idx_lowbits + set_bit_ct;
    if (new_write_idx_lowbits >= kBitsPerWord) {
      *target_iter++ = cur_output_word;
      // ...and these are the bits that fell off
      if (write_idx_lowbits) {
        cur_output_word = extracted_bits >> (kBitsPerWord - write_idx_lowbits);
      } else {
        cur_output_word = 0;
      }
    }
    write_idx_lowbits = new_write_idx_lowbits % kBitsPerWord;
  }
  if (write_idx_lowbits) {
    *target_iter = cur_output_word;
  }
}

void ExpandBytearrNested(const void* __restrict compact_bitarr, const uintptr_t* __restrict mid_bitarr, const uintptr_t* __restrict top_expand_mask, uint32_t word_ct, uint32_t mid_popcount, uint32_t mid_start_bit, uintptr_t* __restrict mid_target, uintptr_t* __restrict compact_target) {
  assert(mid_popcount);
  const uint32_t leading_byte_ct = 1 + ((mid_popcount - 1) % kBitsPerWord) / CHAR_BIT;
  uintptr_t compact_read_word = SubwordLoad(compact_bitarr, leading_byte_ct);
  uint32_t compact_idx_lowbits = CHAR_BIT * (sizeof(intptr_t) - leading_byte_ct);
  const uintptr_t* compact_bitarr_iter = R_CAST(const uintptr_t*, &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]));
  const uintptr_t* mid_bitarr_iter = mid_bitarr;
  uint32_t mid_idx_lowbits = mid_start_bit;
  uintptr_t mid_read_word = (*mid_bitarr_iter) >> mid_start_bit;
  for (uint32_t widx = 0; widx < word_ct; ++widx) {
    const uintptr_t top_word = top_expand_mask[widx];
    uintptr_t mid_write_word = 0;
    uintptr_t compact_write_word = 0;
    if (top_word) {
      const uint32_t top_set_ct = PopcountWord(top_word);
      uint32_t next_mid_idx_lowbits = mid_idx_lowbits + top_set_ct;
      if (next_mid_idx_lowbits <= kBitsPerWord) {
        mid_write_word = _pdep_u64(mid_read_word, top_word);
        if (top_set_ct != kBitsPerWord) {
          mid_read_word = mid_read_word >> top_set_ct;
        } else {
          // avoid nasal demons
          mid_read_word = 0;
        }
      } else {
        uintptr_t next_mid_read_word = *(++mid_bitarr_iter);
        next_mid_idx_lowbits -= kBitsPerWord;
        mid_read_word |= next_mid_read_word << (kBitsPerWord - mid_idx_lowbits);
        mid_write_word = _pdep_u64(mid_read_word, top_word);
        if (next_mid_idx_lowbits != kBitsPerWord) {
          mid_read_word = next_mid_read_word >> next_mid_idx_lowbits;
        } else {
          mid_read_word = 0;
        }
      }
      mid_idx_lowbits = next_mid_idx_lowbits;
      if (mid_write_word) {
        const uint32_t mid_set_ct = PopcountWord(mid_write_word);
        uint32_t next_compact_idx_lowbits = compact_idx_lowbits + mid_set_ct;
        if (next_compact_idx_lowbits <= kBitsPerWord) {
          compact_write_word = _pdep_u64(compact_read_word, mid_write_word);
          if (mid_set_ct != kBitsPerWord) {
            compact_read_word = compact_read_word >> mid_set_ct;
          } else {
            compact_read_word = 0;
          }
        } else {
#ifdef __arm__
#  error "Unaligned accesses in ExpandBytearrNested()."
#endif
          uintptr_t next_compact_word = *compact_bitarr_iter++;
          next_compact_idx_lowbits -= kBitsPerWord;
          compact_read_word |= next_compact_word << (kBitsPerWord - compact_idx_lowbits);
          compact_write_word = _pdep_u64(compact_read_word, mid_write_word);
          if (next_compact_idx_lowbits != kBitsPerWord) {
            compact_read_word = next_compact_word >> next_compact_idx_lowbits;
          } else {
            compact_read_word = 0;
          }
        }
        compact_idx_lowbits = next_compact_idx_lowbits;
      }
    }
    mid_target[widx] = mid_write_word;
    compact_target[widx] = compact_write_word;
  }
}

void ExpandThenSubsetBytearrNested(const void* __restrict compact_bitarr, const uintptr_t* __restrict mid_bitarr, const uintptr_t* __restrict top_expand_mask, const uintptr_t* __restrict subset_mask, uint32_t subset_size, uint32_t mid_popcount, uint32_t mid_start_bit, uintptr_t* __restrict mid_target, uintptr_t* __restrict compact_target) {
  assert(mid_popcount);
  const uint32_t leading_byte_ct = 1 + ((mid_popcount - 1) % kBitsPerWord) / CHAR_BIT;
  uintptr_t compact_read_word = SubwordLoad(compact_bitarr, leading_byte_ct);
  uint32_t compact_idx_lowbits = CHAR_BIT * (sizeof(intptr_t) - leading_byte_ct);
  const uintptr_t* compact_bitarr_alias = R_CAST(const uintptr_t*, &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]));
  const uintptr_t* mid_bitarr_iter = mid_bitarr;
  const uint32_t subset_size_lowbits = subset_size % kBitsPerWord;
  const uint32_t write_widx_last = subset_size / kBitsPerWord;
  uintptr_t mid_read_word = (*mid_bitarr_iter) >> mid_start_bit;
  uintptr_t mid_output_word = 0;
  uintptr_t compact_output_word = 0;
  uint32_t mid_idx_lowbits = mid_start_bit;
  uint32_t compact_widx = UINT32_MAX;  // deliberate overflow
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = 0;
  uint32_t write_widx = 0;

  // bugfix (5 Feb 2018): missed a case in sparse subset_mask optimization
  uint32_t mid_set_skip = 0;
  while ((write_widx != write_widx_last) || (write_idx_lowbits != subset_size_lowbits)) {
    uintptr_t subset_word;
    uintptr_t mid_expanded_bits;
    uint32_t mid_set_ct;
    while (1) {
      ++read_widx;
      uintptr_t top_word = top_expand_mask[read_widx];
      subset_word = subset_mask[read_widx];
      mid_expanded_bits = 0;
      if (top_word) {
        uint32_t top_set_ct = PopcountWord(top_word);
        uint32_t next_mid_idx_lowbits = mid_idx_lowbits + top_set_ct;
        if (next_mid_idx_lowbits <= kBitsPerWord) {
          mid_expanded_bits = _pdep_u64(mid_read_word, top_word);
          if (top_set_ct != kBitsPerWord) {
            mid_read_word = mid_read_word >> top_set_ct;
          } else {
            // avoid nasal demons
            mid_read_word = 0;
          }
        } else {
          uintptr_t next_mid_read_word = *(++mid_bitarr_iter);
          next_mid_idx_lowbits -= kBitsPerWord;
          mid_read_word |= next_mid_read_word << (kBitsPerWord - mid_idx_lowbits);
          mid_expanded_bits = _pdep_u64(mid_read_word, top_word);
          if (next_mid_idx_lowbits != kBitsPerWord) {
            mid_read_word = next_mid_read_word >> next_mid_idx_lowbits;
          } else {
            mid_read_word = 0;
          }
        }
        mid_idx_lowbits = next_mid_idx_lowbits;
      }
      mid_set_ct = PopcountWord(mid_expanded_bits);
      if (subset_word) {
        break;
      }
      mid_set_skip += mid_set_ct;
    }

    uintptr_t mid_extracted_bits = 0;
    uintptr_t compact_extracted_bits = 0;
    uint32_t set_bit_ct = PopcountWord(subset_word);
    if (mid_expanded_bits & subset_word) {
      // lazy load
      compact_idx_lowbits += mid_set_skip;
      if (compact_idx_lowbits >= kBitsPerWord) {
        compact_widx += compact_idx_lowbits / kBitsPerWord;
        compact_idx_lowbits = compact_idx_lowbits % kBitsPerWord;
#ifdef __arm__
#  error "Unaligned accesses in ExpandThenSubsetBytearrNested()."
#endif
        compact_read_word = compact_bitarr_alias[compact_widx] >> compact_idx_lowbits;
      } else {
        compact_read_word = compact_read_word >> mid_set_skip;
      }
      uint32_t next_compact_idx_lowbits = compact_idx_lowbits + mid_set_ct;
      uintptr_t compact_expanded_bits;
      if (next_compact_idx_lowbits <= kBitsPerWord) {
        compact_expanded_bits = _pdep_u64(compact_read_word, mid_expanded_bits);
        if (mid_set_ct != kBitsPerWord) {
          compact_read_word = compact_read_word >> mid_set_ct;
        }
      } else {
        uintptr_t next_compact_word = compact_bitarr_alias[++compact_widx];
        next_compact_idx_lowbits -= kBitsPerWord;
        compact_read_word |= next_compact_word << (kBitsPerWord - compact_idx_lowbits);
        compact_expanded_bits = _pdep_u64(compact_read_word, mid_expanded_bits);
        if (next_compact_idx_lowbits != kBitsPerWord) {
          compact_read_word = next_compact_word >> next_compact_idx_lowbits;
        }
      }
      compact_extracted_bits = _pext_u64(compact_expanded_bits, subset_word);
      mid_extracted_bits = _pext_u64(mid_expanded_bits, subset_word);
      compact_idx_lowbits = next_compact_idx_lowbits;
      compact_output_word |= compact_extracted_bits << write_idx_lowbits;
      mid_output_word |= mid_extracted_bits << write_idx_lowbits;
      mid_set_skip = 0;
    } else {
      mid_set_skip += mid_set_ct;
    }
    const uint32_t new_write_idx_lowbits = write_idx_lowbits + set_bit_ct;
    if (new_write_idx_lowbits >= kBitsPerWord) {
      mid_target[write_widx] = mid_output_word;
      compact_target[write_widx] = compact_output_word;
      ++write_widx;
      if (write_idx_lowbits) {
        mid_output_word = mid_extracted_bits >> (kBitsPerWord - write_idx_lowbits);
        compact_output_word = compact_extracted_bits >> (kBitsPerWord - write_idx_lowbits);
      } else {
        mid_output_word = 0;
        compact_output_word = 0;
      }
    }
    write_idx_lowbits = new_write_idx_lowbits % kBitsPerWord;
  }
  if (write_idx_lowbits) {
    mid_target[write_widx] = mid_output_word;
    compact_target[write_widx] = compact_output_word;
  }
}
#else  // !USE_AVX2
void CopyBitarrSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t bit_idx_end, uintptr_t* __restrict output_bitarr) {
  const uint32_t bit_idx_end_lowbits = bit_idx_end % kBitsPerWord;
  uintptr_t* output_bitarr_iter = output_bitarr;
  uintptr_t* output_bitarr_last = &(output_bitarr[bit_idx_end / kBitsPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = 0;
  while ((output_bitarr_iter != output_bitarr_last) || (write_idx_lowbits != bit_idx_end_lowbits)) {
    uintptr_t cur_mask_word;
    // sparse subset_mask optimization
    // guaranteed to terminate since there's at least one more set bit
    do {
      cur_mask_word = subset_mask[++read_widx];
    } while (!cur_mask_word);
    uintptr_t cur_masked_input_word = raw_bitarr[read_widx] & cur_mask_word;
    const uint32_t cur_mask_popcount = PopcountWord(cur_mask_word);
    uintptr_t subsetted_input_word = 0;
    if (cur_masked_input_word) {
      const uintptr_t cur_inv_mask = ~cur_mask_word;
      do {
        const uint32_t read_uidx_nz_start_lowbits = ctzw(cur_masked_input_word);
        const uintptr_t cur_inv_mask_shifted = cur_inv_mask >> read_uidx_nz_start_lowbits;
        if (!cur_inv_mask_shifted) {
          subsetted_input_word |= cur_masked_input_word >> (kBitsPerWord - cur_mask_popcount);
          break;
        }
        const uint32_t cur_read_end = ctzw(cur_inv_mask_shifted) + read_uidx_nz_start_lowbits;
        // this seems to optimize better than (k1LU << cur_read_end) - k1LU
        // todo: check if/when that's true elsewhere
        const uintptr_t lowmask = (~k0LU) >> (kBitsPerWord - cur_read_end);
        const uintptr_t bits_to_copy = cur_masked_input_word & lowmask;
        cur_masked_input_word -= bits_to_copy;
        // todo: check if a less-popcounty implementation should be used in
        // non-SSE4.2 case
        const uint32_t cur_write_end = PopcountWord(cur_mask_word & lowmask);
        subsetted_input_word |= bits_to_copy >> (cur_read_end - cur_write_end);
      } while (cur_masked_input_word);
    }
    cur_output_word |= subsetted_input_word << write_idx_lowbits;
    const uint32_t new_write_idx_lowbits = write_idx_lowbits + cur_mask_popcount;
    if (new_write_idx_lowbits >= kBitsPerWord) {
      *output_bitarr_iter++ = cur_output_word;
      // ...and these are the bits that fell off
      // bugfix: unsafe to right-shift 64
      if (write_idx_lowbits) {
        cur_output_word = subsetted_input_word >> (kBitsPerWord - write_idx_lowbits);
      } else {
        cur_output_word = 0;
      }
    }
    write_idx_lowbits = new_write_idx_lowbits % kBitsPerWord;
  }
  if (write_idx_lowbits) {
    *output_bitarr_iter = cur_output_word;
  }
}

// Basic SSE2 implementation of Lauradoux/Walisch popcount.
uintptr_t PopcountVecsNoSse42(const VecW* bit_vvec, uintptr_t vec_ct) {
  // popcounts vptr[0..(vec_ct-1)].  Assumes vec_ct is a multiple of 3 (0 ok).
  assert(!(vec_ct % 3));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  const VecW* bit_vvec_iter = bit_vvec;
  uintptr_t tot = 0;
  while (1) {
    UniVec acc;
    acc.vw = vecw_setzero();
    const VecW* bit_vvec_stop;
    if (vec_ct < 30) {
      if (!vec_ct) {
        return tot;
      }
      bit_vvec_stop = &(bit_vvec_iter[vec_ct]);
      vec_ct = 0;
    } else {
      bit_vvec_stop = &(bit_vvec_iter[30]);
      vec_ct -= 30;
    }
    do {
      VecW count1 = *bit_vvec_iter++;
      VecW count2 = *bit_vvec_iter++;
      VecW half1 = *bit_vvec_iter++;
      VecW half2 = vecw_srli(half1, 1) & m1;
      half1 = half1 & m1;
      // Two bits can represent values from 0-3, so make each pair in count1
      // count2 store a partial bitcount covering themselves AND another bit
      // from elsewhere.
      count1 = count1 - (vecw_srli(count1, 1) & m1);
      count2 = count2 - (vecw_srli(count2, 1) & m1);
      count1 = count1 + half1;
      count2 = count2 + half2;
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = (count1 & m2) + (vecw_srli(count1, 2) & m2);
      count1 = count1 + (count2 & m2) + (vecw_srli(count2, 2) & m2);
      // Accumulator stores sixteen 0-255 counts in parallel.
      // (32 in AVX2 case, 4 in 32-bit case)
      acc.vw = acc.vw + (count1 & m4) + (vecw_srli(count1, 4) & m4);
    } while (bit_vvec_iter < bit_vvec_stop);
    acc.vw = (acc.vw & m8) + (vecw_srli(acc.vw, 8) & m8);
    tot += UniVecHsum16(acc);
  }
}

void ExpandBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, uint32_t word_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  ZeroWArr(word_ct, target);
  const uintptr_t* compact_bitarr_alias = S_CAST(const uintptr_t*, compact_bitarr);
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t compact_widx_last = expand_sizex_m1 / kBitsPerWord;
  uint32_t compact_widx = 0;
  uint32_t compact_idx_lowbits = read_start_bit;
  uint32_t loop_len = kBitsPerWord;
  uint32_t write_uidx = 0;
  while (1) {
    uintptr_t compact_word;
    if (compact_widx >= compact_widx_last) {
      if (compact_widx > compact_widx_last) {
        return;
      }
      loop_len = 1 + (expand_sizex_m1 % kBitsPerWord);
      // avoid possible segfault
      compact_word = SubwordLoad(&(compact_bitarr_alias[compact_widx]), DivUp(loop_len, CHAR_BIT));
    } else {
#ifdef __arm__
#  error "Unaligned accesses in ExpandBytearr()."
#endif
      compact_word = compact_bitarr_alias[compact_widx];
    }
    for (; compact_idx_lowbits < loop_len; ++compact_idx_lowbits, ++write_uidx) {
      write_uidx = AdvTo1Bit(expand_mask, write_uidx);
      // bugfix: can't just use (compact_word & 1) and compact_word >>= 1,
      // since we may skip the first bit on the first loop iteration
      if ((compact_word >> compact_idx_lowbits) & 1) {
        SetBit(write_uidx, target);
      }
    }
    compact_idx_lowbits = 0;
    ++compact_widx;
  }
}

void ExpandThenSubsetBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, const uintptr_t* __restrict subset_mask, uint32_t expand_size, uint32_t subset_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t leading_byte_ct = 1 + (expand_sizex_m1 % kBitsPerWord) / CHAR_BIT;
  uint32_t read_idx_lowbits = CHAR_BIT * (sizeof(intptr_t) - leading_byte_ct);
  uintptr_t compact_read_word = SubwordLoad(compact_bitarr, leading_byte_ct) << read_idx_lowbits;
  read_idx_lowbits += read_start_bit;
  const uintptr_t* compact_bitarr_iter = R_CAST(const uintptr_t*, &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]));
  const uint32_t subset_size_lowbits = subset_size % kBitsPerWord;
  uintptr_t* target_iter = target;
  uintptr_t* target_last = &(target[subset_size / kBitsPerWord]);
  uintptr_t compact_write_word = 0;
  uint32_t read_widx = 0;
  // further improvement is probably possible (e.g. use AVX2 lazy-load), but
  // I'll postpone for now
  uint32_t write_idx_lowbits = 0;
  while ((target_iter != target_last) || (write_idx_lowbits != subset_size_lowbits)) {
    const uintptr_t subset_word = subset_mask[read_widx];
    const uintptr_t expand_word = expand_mask[read_widx];
    ++read_widx;
    uintptr_t tmp_compact_write_word = 0;
    if (expand_word) {
      const uint32_t expand_bit_ct = PopcountWord(expand_word);
      uint32_t read_idx_lowbits_end = read_idx_lowbits + expand_bit_ct;
      uintptr_t tmp_compact_read_word = 0;
      if (read_idx_lowbits != kBitsPerWord) {
        tmp_compact_read_word = compact_read_word >> read_idx_lowbits;
      }
      if (read_idx_lowbits_end > kBitsPerWord) {
#ifdef __arm__
#  error "Unaligned accesses in ExpandThenSubsetBytearr()."
#endif
        compact_read_word = *compact_bitarr_iter++;
        tmp_compact_read_word |= compact_read_word << (kBitsPerWord - read_idx_lowbits);
        read_idx_lowbits_end -= kBitsPerWord;
      }
      tmp_compact_read_word = bzhi_max(tmp_compact_read_word, expand_bit_ct);
      read_idx_lowbits = read_idx_lowbits_end;
      if (tmp_compact_read_word) {
        uintptr_t cur_intersect = subset_word & expand_word;
        while (cur_intersect) {
          const uintptr_t cur_intersect_and_arg = cur_intersect - k1LU;
          const uintptr_t lowmask = (cur_intersect ^ cur_intersect_and_arg) >> 1;
          const uint32_t read_idx_offset = PopcountWord(expand_word & lowmask);
          uintptr_t shifted_compact_read_word = tmp_compact_read_word >> read_idx_offset;
          if (shifted_compact_read_word & 1) {
            tmp_compact_write_word |= (k1LU << PopcountWord(subset_word & lowmask));
            if (shifted_compact_read_word == 1) {
              break;
            }
          }
          cur_intersect &= cur_intersect_and_arg;
        }
      }
      compact_write_word |= tmp_compact_write_word << write_idx_lowbits;
    }
    const uint32_t write_idx_lowbits_end = write_idx_lowbits + PopcountWord(subset_word);
    if (write_idx_lowbits_end >= kBitsPerWord) {
      *target_iter++ = compact_write_word;
      if (write_idx_lowbits) {
        compact_write_word = tmp_compact_write_word >> (kBitsPerWord - write_idx_lowbits);
      } else {
        compact_write_word = 0;
      }
    }
    write_idx_lowbits = write_idx_lowbits_end % kBitsPerWord;
  }
  if (write_idx_lowbits) {
    *target_iter = compact_write_word;
  }
}

// compact_bitarr := phaseinfo
// mid_bitarr := phasepresent, [1 + het_ct]
// top_expand_mask := all_hets, [raw_sample_ct]
void ExpandBytearrNested(const void* __restrict compact_bitarr, const uintptr_t* __restrict mid_bitarr, const uintptr_t* __restrict top_expand_mask, uint32_t word_ct, uint32_t mid_popcount, uint32_t mid_start_bit, uintptr_t* __restrict mid_target, uintptr_t* __restrict compact_target) {
  ZeroWArr(word_ct, mid_target);
  ZeroWArr(word_ct, compact_target);
  const uintptr_t* compact_bitarr_alias = S_CAST(const uintptr_t*, compact_bitarr);
  const uint32_t mid_popcount_m1 = mid_popcount - 1;
  const uint32_t compact_widx_last = mid_popcount_m1 / kBitsPerWord;
  uint32_t mid_idx = mid_start_bit;
  uint32_t compact_widx = 0;
  // can allow compact_idx_lowbits to be initialized to nonzero
  uint32_t loop_len = kBitsPerWord;
  uint32_t write_uidx = 0;
  while (1) {
    uintptr_t compact_word;
    if (compact_widx >= compact_widx_last) {
      if (compact_widx > compact_widx_last) {
        return;
      }
      loop_len = 1 + (mid_popcount_m1 % kBitsPerWord);
      // avoid possible segfault
      compact_word = SubwordLoad(&(compact_bitarr_alias[compact_widx]), DivUp(loop_len, CHAR_BIT));
    } else {
#ifdef __arm__
#  error "Unaligned accesses in ExpandBytearrNested()."
#endif
      compact_word = compact_bitarr_alias[compact_widx];
    }
    for (uint32_t compact_idx_lowbits = 0; compact_idx_lowbits < loop_len; ++mid_idx, ++write_uidx) {
      write_uidx = AdvTo1Bit(top_expand_mask, write_uidx);
      if (IsSet(mid_bitarr, mid_idx)) {
        const uintptr_t new_bit = k1LU << (write_uidx % kBitsPerWord);
        const uint32_t sample_widx = write_uidx / kBitsPerWord;
        mid_target[sample_widx] |= new_bit;
        compact_target[sample_widx] |= new_bit * (compact_word & 1);
        compact_word >>= 1;
        ++compact_idx_lowbits;
      }
    }
    ++compact_widx;
  }
}

void ExpandThenSubsetBytearrNested(const void* __restrict compact_bitarr, const uintptr_t* __restrict mid_bitarr, const uintptr_t* __restrict top_expand_mask, const uintptr_t* __restrict subset_mask, uint32_t subset_size, uint32_t mid_popcount, uint32_t mid_start_bit, uintptr_t* __restrict mid_target, uintptr_t* __restrict compact_target) {
  assert(mid_popcount);
  const uint32_t leading_byte_ct = 1 + ((mid_popcount - 1) % kBitsPerWord) / CHAR_BIT;
  uint32_t compact_idx_lowbits = CHAR_BIT * (sizeof(intptr_t) - leading_byte_ct);
  uintptr_t compact_read_word = SubwordLoad(compact_bitarr, leading_byte_ct) << compact_idx_lowbits;
  const uintptr_t* compact_bitarr_iter = R_CAST(const uintptr_t*, &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]));
  // bugfix (12 Apr 2018): need to round down here
  const uint32_t subset_size_dl = subset_size / kBitsPerWord;
  const uint32_t subset_size_lowbits = subset_size % kBitsPerWord;
  const uintptr_t* mid_read_iter = mid_bitarr;
  uintptr_t mid_read_word = *mid_read_iter++;
  uintptr_t mid_write_word = 0;
  uintptr_t compact_write_word = 0;
  uint32_t mid_idx_lowbits = mid_start_bit;
  uint32_t write_idx_lowbits = 0;
  uint32_t write_widx = 0;
  uint32_t read_widx = 0;
  while ((write_widx != subset_size_dl) || (write_idx_lowbits != subset_size_lowbits)) {
    const uintptr_t subset_word = subset_mask[read_widx];
    const uintptr_t top_word = top_expand_mask[read_widx];
    ++read_widx;
    uintptr_t tmp_mid_write_word = 0;
    uintptr_t tmp_compact_write_word = 0;
    if (top_word) {
      const uint32_t top_set_ct = PopcountWord(top_word);
      uint32_t mid_idx_lowbits_end = mid_idx_lowbits + top_set_ct;
      uintptr_t tmp_mid_read_word = 0;
      if (mid_idx_lowbits != kBitsPerWord) {
        tmp_mid_read_word = mid_read_word >> mid_idx_lowbits;
      }
      if (mid_idx_lowbits_end > kBitsPerWord) {
        // be paranoid for now re: reading an extra word off the end of
        // mid_bitarr
        mid_read_word = *mid_read_iter++;
        tmp_mid_read_word |= mid_read_word << (kBitsPerWord - mid_idx_lowbits);
        mid_idx_lowbits_end -= kBitsPerWord;
      }
      tmp_mid_read_word = bzhi_max(tmp_mid_read_word, top_set_ct);
      mid_idx_lowbits = mid_idx_lowbits_end;
      if (tmp_mid_read_word) {
        const uint32_t mid_set_ct = PopcountWord(tmp_mid_read_word);
        uintptr_t tmp_compact_read_word;
        if (compact_idx_lowbits != kBitsPerWord) {
          const uint32_t compact_idx_lowbits_end = compact_idx_lowbits + mid_set_ct;
          tmp_compact_read_word = compact_read_word >> compact_idx_lowbits;
          // avoid reading off end of compact_bitarr here
          if (compact_idx_lowbits_end <= kBitsPerWord) {
            compact_idx_lowbits = compact_idx_lowbits_end;
          } else {
#ifdef __arm__
#  error "Unaligned accesses in ExpandThenSubsetBytearrNested()."
#endif
            compact_read_word = *compact_bitarr_iter++;
            tmp_compact_read_word |= compact_read_word << (kBitsPerWord - compact_idx_lowbits);
            compact_idx_lowbits = compact_idx_lowbits_end - kBitsPerWord;
          }
        } else {
          // special case, can't right-shift 64
          compact_read_word = *compact_bitarr_iter++;
          compact_idx_lowbits = mid_set_ct;
          tmp_compact_read_word = compact_read_word;
        }
        tmp_compact_read_word = bzhi_max(tmp_compact_read_word, mid_set_ct);

        uintptr_t cur_masked_top = subset_word & top_word;
        while (cur_masked_top) {
          const uintptr_t cur_masked_top_and_arg = cur_masked_top - k1LU;
          const uintptr_t lowmask = (cur_masked_top ^ cur_masked_top_and_arg) >> 1;
          const uint32_t read_idx_offset = PopcountWord(top_word & lowmask);
          uintptr_t shifted_mid_read_word = tmp_mid_read_word >> read_idx_offset;
          if (shifted_mid_read_word & 1) {
            // bugfix (7 Sep 2017): forgot the "k1LU << " part of this
            const uintptr_t cur_bit = k1LU << PopcountWord(subset_word & lowmask);
            tmp_mid_write_word |= cur_bit;
            tmp_compact_write_word += cur_bit * ((tmp_compact_read_word >> (mid_set_ct - PopcountWord(shifted_mid_read_word))) & 1);
            if (shifted_mid_read_word == 1) {
              break;
            }
          }
          cur_masked_top &= cur_masked_top_and_arg;
        }
      }
      mid_write_word |= tmp_mid_write_word << write_idx_lowbits;
      compact_write_word |= tmp_compact_write_word << write_idx_lowbits;
    }
    const uint32_t write_idx_lowbits_end = write_idx_lowbits + PopcountWord(subset_word);
    if (write_idx_lowbits_end >= kBitsPerWord) {
      mid_target[write_widx] = mid_write_word;
      compact_target[write_widx] = compact_write_word;
      ++write_widx;
      if (write_idx_lowbits) {
        const uint32_t rshift = kBitsPerWord - write_idx_lowbits;
        mid_write_word = tmp_mid_write_word >> rshift;
        compact_write_word = tmp_compact_write_word >> rshift;
      } else {
        mid_write_word = 0;
        compact_write_word = 0;
      }
    }
    write_idx_lowbits = write_idx_lowbits_end % kBitsPerWord;
  }
  if (write_idx_lowbits) {
    mid_target[write_widx] = mid_write_word;
    compact_target[write_widx] = compact_write_word;
  }
}
#endif
uintptr_t PopcountBytes(const unsigned char* bitarr, uintptr_t byte_ct) {
  const uint32_t lead_byte_ct = (-R_CAST(uintptr_t, bitarr)) % kBytesPerVec;
  uintptr_t tot = 0;
  const uintptr_t* bitarr_iter;
  uint32_t trail_byte_ct;
  // bugfix: had wrong condition here
  if (byte_ct >= lead_byte_ct) {
#ifdef __LP64__
    const uint32_t word_rem = lead_byte_ct % kBytesPerWord;
    if (word_rem) {
      tot = PopcountWord(ProperSubwordLoad(bitarr, word_rem));
    }
    bitarr_iter = R_CAST(const uintptr_t*, &(bitarr[word_rem]));
    if (lead_byte_ct >= kBytesPerWord) {
      tot += PopcountWord(*bitarr_iter++);
#  ifdef USE_AVX2
      if (lead_byte_ct >= 2 * kBytesPerWord) {
        tot += PopcountWord(*bitarr_iter++);
        if (lead_byte_ct >= 3 * kBytesPerWord) {
          tot += PopcountWord(*bitarr_iter++);
        }
      }
#  endif
    }
#else
    if (lead_byte_ct) {
      tot = PopcountWord(ProperSubwordLoad(bitarr, lead_byte_ct));
    }
    bitarr_iter = R_CAST(const uintptr_t*, &(bitarr[lead_byte_ct]));
#endif
    byte_ct -= lead_byte_ct;
    const uintptr_t word_ct = byte_ct / kBytesPerWord;
    // vec-alignment required here
    tot += PopcountWords(bitarr_iter, word_ct);
    bitarr_iter = &(bitarr_iter[word_ct]);
    trail_byte_ct = byte_ct % kBytesPerWord;
  } else {
    bitarr_iter = R_CAST(const uintptr_t*, bitarr);
    // this may still be >= kBytesPerWord, so can't remove loop
    trail_byte_ct = byte_ct;
  }
  while (1) {
    uintptr_t cur_word;
    if (trail_byte_ct < kBytesPerWord) {
      if (!trail_byte_ct) {
        return tot;
      }
      cur_word = ProperSubwordLoad(bitarr_iter, trail_byte_ct);
      trail_byte_ct = 0;
    } else {
      cur_word = *bitarr_iter++;
      trail_byte_ct -= kBytesPerWord;
    }
    tot += PopcountWord(cur_word);
  }
}

uintptr_t PopcountBytesMasked(const unsigned char* bitarr, const uintptr_t* mask_arr, uintptr_t byte_ct) {
  // could detect aligned case, but that shouldn't happen often enough?
  const uintptr_t word_ct = byte_ct / kBytesPerWord;
#ifdef USE_SSE42
  const uintptr_t* bitarr_alias = R_CAST(const uintptr_t*, bitarr);
  uintptr_t tot = 0;
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    tot += PopcountWord(bitarr_alias[widx] & mask_arr[widx]);
  }
  const uint32_t trail_byte_ct = byte_ct % kBytesPerWord;
  if (trail_byte_ct) {
    uintptr_t cur_word = ProperSubwordLoad(&(bitarr_alias[word_ct]), trail_byte_ct);
    tot += PopcountWord(cur_word & mask_arr[word_ct]);
  }
  return tot;
#else
  const uintptr_t* bitarr_iter = R_CAST(const uintptr_t*, bitarr);
  const uintptr_t mainblock_word_ct = word_ct - (word_ct % (24 / kBytesPerWord));
  const uintptr_t* bitarr_24b_end = &(bitarr_iter[mainblock_word_ct]);
  const uintptr_t* mask_arr_iter = mask_arr;
  uintptr_t tot = 0;
  while (bitarr_iter < bitarr_24b_end) {
    uintptr_t loader = (*bitarr_iter++) & (*mask_arr_iter++);
    uintptr_t ulii = loader - ((loader >> 1) & kMask5555);
    loader = (*bitarr_iter++) & (*mask_arr_iter++);
    uintptr_t uljj = loader - ((loader >> 1) & kMask5555);
    loader = (*bitarr_iter++) & (*mask_arr_iter++);
    ulii += (loader >> 1) & kMask5555;
    uljj += loader & kMask5555;
    ulii = (ulii & kMask3333) + ((ulii >> 2) & kMask3333);
    ulii += (uljj & kMask3333) + ((uljj >> 2) & kMask3333);
    uintptr_t tmp_stor = (ulii & kMask0F0F) + ((ulii >> 4) & kMask0F0F);

#  ifndef __LP64__
    loader = (*bitarr_iter++) & (*mask_arr_iter++);
    ulii = loader - ((loader >> 1) & kMask5555);
    loader = (*bitarr_iter++) & (*mask_arr_iter++);
    uljj = loader - ((loader >> 1) & kMask5555);
    loader = (*bitarr_iter++) & (*mask_arr_iter++);
    ulii += (loader >> 1) & kMask5555;
    uljj += loader & kMask5555;
    ulii = (ulii & kMask3333) + ((ulii >> 2) & kMask3333);
    ulii += (uljj & kMask3333) + ((uljj >> 2) & kMask3333);
    tmp_stor += (ulii & kMask0F0F) + ((ulii >> 4) & kMask0F0F);
#  endif

    // 32-bit case: each 8-bit slot stores a number in 0..48.  Multiplying by
    // 0x01010101 is equivalent to the left-shifts and adds we need to sum
    // those four 8-bit numbers in the high-order slot.
    // 64-bit case: each 8-bit slot stores a number in 0..24.
    tot += (tmp_stor * kMask0101) >> (kBitsPerWord - 8);
  }
  uint32_t trail_byte_ct = byte_ct - (mainblock_word_ct * kBytesPerWord);
  while (1) {
    uintptr_t cur_word;
    if (trail_byte_ct < kBytesPerWord) {
      if (!trail_byte_ct) {
        return tot;
      }
      cur_word = ProperSubwordLoad(bitarr_iter, trail_byte_ct);
      trail_byte_ct = 0;
    } else {
      cur_word = *bitarr_iter++;
      trail_byte_ct -= kBytesPerWord;
    }
    tot += PopcountWord(cur_word & (*mask_arr_iter++));
  }
#endif
}

void FillCumulativePopcounts(const uintptr_t* subset_mask, uint32_t word_ct, uint32_t* cumulative_popcounts) {
  assert(word_ct);
  const uint32_t word_ct_m1 = word_ct - 1;
  uint32_t cur_sum = 0;
  for (uint32_t widx = 0; widx < word_ct_m1; ++widx) {
    cumulative_popcounts[widx] = cur_sum;
    cur_sum += PopcountWord(subset_mask[widx]);
  }
  cumulative_popcounts[word_ct_m1] = cur_sum;
}

void UidxsToIdxs(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, const uint32_t idx_list_len, uint32_t* idx_list) {
  uint32_t* idx_list_end = &(idx_list[idx_list_len]);
  for (uint32_t* idx_list_iter = idx_list; idx_list_iter != idx_list_end; ++idx_list_iter) {
    *idx_list_iter = RawToSubsettedPos(subset_mask, subset_cumulative_popcounts, *idx_list_iter);
  }
}


static_assert(kPglBitTransposeBatch == S_CAST(uint32_t, kBitsPerCacheline), "TransposeBitblockInternal() needs to be updated.");
#ifdef __LP64__
void TransposeBitblockInternal(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, void* buf0) {
  // buf must be vector-aligned and have space for 512 bytes
  const uint32_t block_ct = DivUp(write_batch_size, CHAR_BIT);
  // fold the first 6 shuffles into the initial ingestion loop
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t write_mmui_stride = kMovemaskUintPerWord * write_ul_stride;
  const uint32_t read_batch_rem = kBitsPerCacheline - read_batch_size;
  MovemaskUint* target_iter0 = R_CAST(MovemaskUint*, write_iter);
  const uint32_t full_block_ct = write_batch_size / 8;
  const uint32_t loop_vec_ct = 4 * DivUp(read_batch_size, kBytesPerVec * 4);
  for (uint32_t block_idx = 0; block_idx < block_ct; ++block_idx) {
    const unsigned char* read_iter_tmp = R_CAST(const unsigned char*, read_iter);
    read_iter_tmp = &(read_iter_tmp[block_idx]);
    unsigned char* initial_target_iter = S_CAST(unsigned char*, buf0);
    for (uint32_t ujj = 0; ujj < read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_byte_stride]);
    }
    memset(initial_target_iter, 0, read_batch_rem);
    const VecW* source_iter = S_CAST(VecW*, buf0);
    if (block_idx == full_block_ct) {
      const uint32_t remainder = write_batch_size % 8;
      const uint32_t remainder_from8 = 8 - remainder;
      const uint32_t remainder_m1 = remainder - 1;
      for (uint32_t vec_idx = 0; vec_idx < loop_vec_ct; ++vec_idx) {
        VecW loader = *source_iter++;
        loader = vecw_slli(loader, remainder_from8);
        uint32_t write_idx_lowbits = remainder_m1;
        while (1) {
          target_iter0[write_mmui_stride * write_idx_lowbits] = vecw_movemask(loader);
          if (!write_idx_lowbits) {
            break;
          }
          loader = vecw_slli(loader, 1);
          --write_idx_lowbits;
        }
        ++target_iter0;
      }
      break;
    }

    MovemaskUint* target_iter1 = &(target_iter0[write_mmui_stride]);
    MovemaskUint* target_iter2 = &(target_iter1[write_mmui_stride]);
    MovemaskUint* target_iter3 = &(target_iter2[write_mmui_stride]);
    MovemaskUint* target_iter4 = &(target_iter3[write_mmui_stride]);
    MovemaskUint* target_iter5 = &(target_iter4[write_mmui_stride]);
    MovemaskUint* target_iter6 = &(target_iter5[write_mmui_stride]);
    MovemaskUint* target_iter7 = &(target_iter6[write_mmui_stride]);
    for (uint32_t vec_idx = 0; vec_idx < loop_vec_ct; ++vec_idx) {
      VecW loader = source_iter[vec_idx];
      target_iter7[vec_idx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      target_iter6[vec_idx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      target_iter5[vec_idx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      target_iter4[vec_idx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      target_iter3[vec_idx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      target_iter2[vec_idx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      target_iter1[vec_idx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      target_iter0[vec_idx] = vecw_movemask(loader);
    }
    target_iter0 = &(target_iter7[write_mmui_stride]);
  }
}
#else  // !__LP64__
static_assert(kWordsPerVec == 1, "TransposeBitblockInternal() needs to be updated.");
void TransposeBitblockInternal(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecW* __restrict buf0, VecW* __restrict buf1) {
  // buf must be vector-aligned and have size 64k
  const uint32_t initial_read_byte_ct = DivUp(write_batch_size, CHAR_BIT);
  // fold the first 6 shuffles into the initial ingestion loop
  const unsigned char* initial_read_iter = R_CAST(const unsigned char*, read_iter);
  const unsigned char* initial_read_end = &(initial_read_iter[initial_read_byte_ct]);
  unsigned char* initial_target_iter = R_CAST(unsigned char*, buf0);
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t read_batch_rem = kBitsPerCacheline - read_batch_size;
  for (; initial_read_iter < initial_read_end; ++initial_read_iter) {
    const unsigned char* read_iter_tmp = initial_read_iter;
    for (uint32_t ujj = 0; ujj < read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_byte_stride]);
    }
    initial_target_iter = memsetua(initial_target_iter, 0, read_batch_rem);
  }

  // third-to-last shuffle, 8 bit spacing -> 4
  const VecW* source_iter = buf0;
  uintptr_t* target_iter0 = buf1;
  const uint32_t write_word_ct = BitCtToWordCt(read_batch_size);
  const uint32_t first_inner_loop_iter_ct = 4 * write_word_ct;
  uint32_t cur_write_skip = 4 * kWordsPerCacheline - first_inner_loop_iter_ct;
  // coincidentally, this also needs to run DivUp(write_batch_size, CHAR_BIT)
  // times
  for (uint32_t uii = 0; uii < initial_read_byte_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 4]);
    for (uint32_t ujj = 0; ujj < first_inner_loop_iter_ct; ++ujj) {
      const uintptr_t source_word_lo = *source_iter++;
      const uintptr_t source_word_hi = *source_iter++;
      uintptr_t target_word0_lo = source_word_lo & kMask0F0F;
      uintptr_t target_word1_lo = (source_word_lo >> 4) & kMask0F0F;
      uintptr_t target_word0_hi = source_word_hi & kMask0F0F;
      uintptr_t target_word1_hi = (source_word_hi >> 4) & kMask0F0F;
      target_word0_lo = (target_word0_lo | (target_word0_lo >> 4)) & kMask00FF;
      target_word1_lo = (target_word1_lo | (target_word1_lo >> 4)) & kMask00FF;
      target_word0_hi = (target_word0_hi | (target_word0_hi >> 4)) & kMask00FF;
      target_word1_hi = (target_word1_hi | (target_word1_hi >> 4)) & kMask00FF;
      target_word0_lo = target_word0_lo | (target_word0_lo >> kBitsPerWordD4);
      target_word1_lo = target_word1_lo | (target_word1_lo >> kBitsPerWordD4);
      target_word0_hi = target_word0_hi | (target_word0_hi >> kBitsPerWordD4);
      target_word1_hi = target_word1_hi | (target_word1_hi >> kBitsPerWordD4);
      *target_iter0++ = S_CAST(Halfword, target_word0_lo) | (target_word0_hi << kBitsPerWordD2);
      *target_iter1++ = S_CAST(Halfword, target_word1_lo) | (target_word1_hi << kBitsPerWordD2);
    }
    source_iter = &(source_iter[2 * cur_write_skip]);
    target_iter0 = &(target_iter1[cur_write_skip]);
  }

  // second-to-last shuffle, 4 bit spacing -> 2
  source_iter = buf1;
  target_iter0 = buf0;
  const uint32_t second_outer_loop_iter_ct = DivUp(write_batch_size, 4);
  const uint32_t second_inner_loop_iter_ct = 2 * write_word_ct;
  cur_write_skip = 2 * kWordsPerCacheline - second_inner_loop_iter_ct;
  for (uint32_t uii = 0; uii < second_outer_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 2]);
    for (uint32_t ujj = 0; ujj < second_inner_loop_iter_ct; ++ujj) {
      const uintptr_t source_word_lo = *source_iter++;
      const uintptr_t source_word_hi = *source_iter++;
      uintptr_t target_word0_lo = source_word_lo & kMask3333;
      uintptr_t target_word1_lo = (source_word_lo >> 2) & kMask3333;
      uintptr_t target_word0_hi = source_word_hi & kMask3333;
      uintptr_t target_word1_hi = (source_word_hi >> 2) & kMask3333;
      target_word0_lo = (target_word0_lo | (target_word0_lo >> 2)) & kMask0F0F;
      target_word1_lo = (target_word1_lo | (target_word1_lo >> 2)) & kMask0F0F;
      target_word0_hi = (target_word0_hi | (target_word0_hi >> 2)) & kMask0F0F;
      target_word1_hi = (target_word1_hi | (target_word1_hi >> 2)) & kMask0F0F;
      target_word0_lo = (target_word0_lo | (target_word0_lo >> 4)) & kMask00FF;
      target_word1_lo = (target_word1_lo | (target_word1_lo >> 4)) & kMask00FF;
      target_word0_hi = (target_word0_hi | (target_word0_hi >> 4)) & kMask00FF;
      target_word1_hi = (target_word1_hi | (target_word1_hi >> 4)) & kMask00FF;
      target_word0_lo = target_word0_lo | (target_word0_lo >> kBitsPerWordD4);
      target_word1_lo = target_word1_lo | (target_word1_lo >> kBitsPerWordD4);
      target_word0_hi = target_word0_hi | (target_word0_hi >> kBitsPerWordD4);
      target_word1_hi = target_word1_hi | (target_word1_hi >> kBitsPerWordD4);
      *target_iter0++ = S_CAST(Halfword, target_word0_lo) | (target_word0_hi << kBitsPerWordD2);
      *target_iter1++ = S_CAST(Halfword, target_word1_lo) | (target_word1_hi << kBitsPerWordD2);
    }
    source_iter = &(source_iter[2 * cur_write_skip]);
    target_iter0 = &(target_iter1[cur_write_skip]);
  }
  // last shuffle, 2 bit spacing -> 1
  source_iter = buf0;
  target_iter0 = write_iter;
  const uint32_t last_loop_iter_ct = DivUp(write_batch_size, 2);
  for (uint32_t uii = 0; uii < last_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    for (uint32_t ujj = 0; ujj < write_word_ct; ++ujj) {
      const uintptr_t source_word_lo = S_CAST(uintptr_t, *source_iter++);
      const uintptr_t source_word_hi = S_CAST(uintptr_t, *source_iter++);
      uintptr_t target_word0_lo = source_word_lo & kMask5555;
      uintptr_t target_word1_lo = (source_word_lo >> 1) & kMask5555;
      uintptr_t target_word0_hi = source_word_hi & kMask5555;
      uintptr_t target_word1_hi = (source_word_hi >> 1) & kMask5555;
      target_word0_lo = (target_word0_lo | (target_word0_lo >> 1)) & kMask3333;
      target_word1_lo = (target_word1_lo | (target_word1_lo >> 1)) & kMask3333;
      target_word0_hi = (target_word0_hi | (target_word0_hi >> 1)) & kMask3333;
      target_word1_hi = (target_word1_hi | (target_word1_hi >> 1)) & kMask3333;
      target_word0_lo = (target_word0_lo | (target_word0_lo >> 2)) & kMask0F0F;
      target_word1_lo = (target_word1_lo | (target_word1_lo >> 2)) & kMask0F0F;
      target_word0_hi = (target_word0_hi | (target_word0_hi >> 2)) & kMask0F0F;
      target_word1_hi = (target_word1_hi | (target_word1_hi >> 2)) & kMask0F0F;
      target_word0_lo = (target_word0_lo | (target_word0_lo >> 4)) & kMask00FF;
      target_word1_lo = (target_word1_lo | (target_word1_lo >> 4)) & kMask00FF;
      target_word0_hi = (target_word0_hi | (target_word0_hi >> 4)) & kMask00FF;
      target_word1_hi = (target_word1_hi | (target_word1_hi >> 4)) & kMask00FF;
      target_word0_lo = target_word0_lo | (target_word0_lo >> kBitsPerWordD4);
      target_word1_lo = target_word1_lo | (target_word1_lo >> kBitsPerWordD4);
      target_word0_hi = target_word0_hi | (target_word0_hi >> kBitsPerWordD4);
      target_word1_hi = target_word1_hi | (target_word1_hi >> kBitsPerWordD4);
      target_iter0[ujj] = S_CAST(Halfword, target_word0_lo) | (target_word0_hi << kBitsPerWordD2);
      target_iter1[ujj] = S_CAST(Halfword, target_word1_lo) | (target_word1_hi << kBitsPerWordD2);
    }
    source_iter = &(source_iter[2 * (kWordsPerCacheline - write_word_ct)]);
    target_iter0 = &(target_iter1[write_ul_stride]);
  }
}
#endif  // !__LP64__

#ifdef __cplusplus
}  // namespace plink2
#endif
