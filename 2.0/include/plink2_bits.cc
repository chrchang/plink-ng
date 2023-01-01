// This library is part of PLINK 2, copyright (C) 2005-2023 Shaun Purcell,
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


#include "plink2_bits.h"

#ifdef __cplusplus
namespace plink2 {
#endif

#if defined(__LP64__) && !defined(USE_AVX2)
// No alignment assumptions.
void Pack32bTo16bMask(const void* words, uintptr_t ct_32b, void* dest) {
  // This is also competitive in the AVX2 case, but never quite beats the
  // simple loop.  (We'd want to enable a similar function for Ryzen,
  // processing one 32-byte vector instead of two 16-byte vectors at a time in
  // the main loop since _mm256_packus_epi16() doesn't do what we want.)
  const VecW m1 = VCONST_W(kMask5555);
#  ifdef USE_SSE42
  const VecW swap12 = vecw_setr8(
      0, 1, 4, 5, 2, 3, 6, 7,
      8, 9, 12, 13, 10, 11, 14, 15);
#  else
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  const VecW* words_valias = R_CAST(const VecW*, words);
  __m128i* dest_alias = R_CAST(__m128i*, dest);
  for (uintptr_t vidx = 0; vidx != ct_32b; ++vidx) {
    VecW vec_lo = vecw_loadu(&(words_valias[2 * vidx])) & m1;
    VecW vec_hi = vecw_loadu(&(words_valias[2 * vidx + 1])) & m1;
#  ifdef USE_SSE42
    // this right-shift-3 + shuffle shortcut saves two operations.
    vec_lo = (vec_lo | vecw_srli(vec_lo, 3)) & m4;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 3)) & m4;
    vec_lo = vecw_shuffle8(swap12, vec_lo);
    vec_hi = vecw_shuffle8(swap12, vec_hi);
#  else
    vec_lo = (vec_lo | vecw_srli(vec_lo, 1)) & m2;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 1)) & m2;
    vec_lo = (vec_lo | vecw_srli(vec_lo, 2)) & m4;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 2)) & m4;
#  endif
    vec_lo = (vec_lo | vecw_srli(vec_lo, 4)) & m8;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 4)) & m8;
    const __m128i vec_packed = _mm_packus_epi16(R_CAST(__m128i, vec_lo), R_CAST(__m128i, vec_hi));
    _mm_storeu_si128(&(dest_alias[vidx]), vec_packed);
  }
}
#endif

#ifdef __x86_64__
VecW vecw_slli_variable_ct(VecW vv, uint32_t ct) {
  return vecw_slli(vv, ct);
}
#else
// Using a lookup table because NEON bit shift functions can only be called with compile-time constants
// https://eigen.tuxfamily.org/bz/show_bug.cgi?id=1631
// https://github.com/VectorCamp/vectorscan/issues/21
VecW vecw_slli_variable_ct(VecW vv, uint32_t ct) {
  switch(ct) {
    default: return vv;
    case 1: return vecw_slli(vv, 1);
    case 2: return vecw_slli(vv, 2);
    case 3: return vecw_slli(vv, 3);
    case 4: return vecw_slli(vv, 4);
    case 5: return vecw_slli(vv, 5);
    case 6: return vecw_slli(vv, 6);
    case 7: return vecw_slli(vv, 7);
  }
}
#endif

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

void FillBitsNz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr) {
  assert(end_idx > start_idx);
  uintptr_t maj_start = start_idx / kBitsPerWord;
  uintptr_t maj_end = end_idx / kBitsPerWord;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bitarr[maj_start] |= (k1LU << (end_idx % kBitsPerWord)) - (k1LU << (start_idx % kBitsPerWord));
  } else {
    bitarr[maj_start] |= ~((k1LU << (start_idx % kBitsPerWord)) - k1LU);
    SetAllWArr(maj_end - maj_start - 1, &(bitarr[maj_start + 1]));
    minor = end_idx % kBitsPerWord;
    if (minor) {
      bitarr[maj_end] |= (k1LU << minor) - k1LU;
    }
  }
}

void ClearBitsNz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr) {
  assert(end_idx > start_idx);
  uintptr_t maj_start = start_idx / kBitsPerWord;
  uintptr_t maj_end = end_idx / kBitsPerWord;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bitarr[maj_start] &= ~((k1LU << (end_idx % kBitsPerWord)) - (k1LU << (start_idx % kBitsPerWord)));
  } else {
    bitarr[maj_start] = bzhi(bitarr[maj_start], start_idx % kBitsPerWord);
    ZeroWArr(maj_end - maj_start - 1, &(bitarr[maj_start + 1]));
    minor = end_idx % kBitsPerWord;
    if (minor) {
      bitarr[maj_end] &= ~((k1LU << minor) - k1LU);
    }
  }
}

void BitvecAnd(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec AND arg_bitvec
#ifdef __LP64__
  VecW* main_bitvvec_iter = R_CAST(VecW*, main_bitvec);
  const VecW* arg_bitvvec_iter = R_CAST(const VecW*, arg_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  // ok, retested this explicit unroll (Jun 2018) and it's still noticeably
  // faster for small cases than the simple loop.  sigh.
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
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    main_bitvec[widx] &= arg_bitvec[widx];
  }
#endif
}

void BitvecInvmask(const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec ANDNOT exclude_bitvec
  // note that this is the reverse of the _mm_andnot() operand order
#ifdef __LP64__
  VecW* main_bitvvec_iter = R_CAST(VecW*, main_bitvec);
  const VecW* exclude_bitvvec_iter = R_CAST(const VecW*, exclude_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    *main_bitvvec_iter = vecw_and_notfirst(*exclude_bitvvec_iter++, *main_bitvvec_iter);
    ++main_bitvvec_iter;
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter = vecw_and_notfirst(*exclude_bitvvec_iter++, *main_bitvvec_iter);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = vecw_and_notfirst(*exclude_bitvvec_iter++, *main_bitvvec_iter);
    ++main_bitvvec_iter;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter = vecw_and_notfirst(*exclude_bitvvec_iter++, *main_bitvvec_iter);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = vecw_and_notfirst(*exclude_bitvvec_iter++, *main_bitvvec_iter);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = vecw_and_notfirst(*exclude_bitvvec_iter++, *main_bitvvec_iter);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = vecw_and_notfirst(*exclude_bitvvec_iter++, *main_bitvvec_iter);
    ++main_bitvvec_iter;
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
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    main_bitvec[widx] &= ~exclude_bitvec[widx];
  }
#endif
}

void BitvecOr(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec) {
  // main_bitvec := main_bitvec OR arg_bitvec
#ifdef __LP64__
  VecW* main_bitvvec_iter = R_CAST(VecW*, main_bitvec);
  const VecW* arg_bitvvec_iter = R_CAST(const VecW*, arg_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= (*arg_bitvvec_iter++);
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] |= arg_bitvec[base_idx];
    main_bitvec[base_idx + 1] |= arg_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] |= arg_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    main_bitvec[widx] |= arg_bitvec[widx];
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
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    main_bitvec[widx] ^= ~k0LU;
  }
#endif
}

void BitvecXorCopy(const uintptr_t* __restrict source1_bitvec, const uintptr_t* __restrict source2_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec) {
#ifdef __LP64__
  VecW* target_bitvvec = R_CAST(VecW*, target_bitvec);
  const VecW* source1_bitvvec = R_CAST(const VecW*, source1_bitvec);
  const VecW* source2_bitvvec = R_CAST(const VecW*, source2_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  for (uintptr_t ulii = 0; ulii != full_vec_ct; ++ulii) {
    target_bitvvec[ulii] = source1_bitvvec[ulii] ^ source2_bitvvec[ulii];
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    target_bitvec[base_idx] = source1_bitvec[base_idx] ^ source2_bitvec[base_idx];
    target_bitvec[base_idx + 1] = source1_bitvec[base_idx + 1] ^ source2_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    target_bitvec[word_ct - 1] = source1_bitvec[word_ct - 1] ^ source2_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    target_bitvec[widx] = source1_bitvec[widx] ^ source2_bitvec[widx];
  }
#endif
}

void BitvecInvertCopy(const uintptr_t* __restrict source_bitvec, uintptr_t word_ct, uintptr_t* __restrict target_bitvec) {
#ifdef __LP64__
  const VecW* source_bitvvec_iter = R_CAST(const VecW*, source_bitvec);
  VecW* target_bitvvec_iter = R_CAST(VecW*, target_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  const VecW all1 = VCONST_W(~k0LU);
  // As of Apple clang 11, this manual unroll is no longer relevant.  todo:
  // check Linux performance, and remove all of these unrolls if perf is good
  // enough without them.
  if (full_vec_ct & 1) {
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
  }
  if (full_vec_ct & 2) {
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    target_bitvec[base_idx] = ~source_bitvec[base_idx];
    target_bitvec[base_idx + 1] = ~source_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    target_bitvec[word_ct - 1] = ~source_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    target_bitvec[widx] = ~source_bitvec[widx];
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
  if (ulii) {
    const uint32_t rval = loc + ctzw(ulii);
    return MINV(rval, ceil);
  }
  const uintptr_t* bitarr_last = &(bitarr[(ceil - 1) / kBitsPerWord]);
  do {
    if (bitarr_iter >= bitarr_last) {
      return ceil;
    }
    ulii = *(++bitarr_iter);
  } while (!ulii);
  const uint32_t rval = S_CAST(uintptr_t, bitarr_iter - bitarr) * kBitsPerWord + ctzw(ulii);
  return MINV(rval, ceil);
}

uintptr_t AdvBoundedTo0Bit(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil) {
  assert(ceil >= 1);
  const uintptr_t* bitarr_ptr = &(bitarr[loc / kBitsPerWord]);
  uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % kBitsPerWord);
  if (ulii) {
    loc += ctzw(ulii);
    return MINV(loc, ceil);
  }
  const uintptr_t* bitarr_last = &(bitarr[(ceil - 1) / kBitsPerWord]);
  do {
    if (bitarr_ptr >= bitarr_last) {
      return ceil;
    }
    ulii = *(++bitarr_ptr);
  } while (ulii == ~k0LU);
  loc = S_CAST(uintptr_t, bitarr_ptr - bitarr) * kBitsPerWord + ctzw(~ulii);
  return MINV(loc, ceil);
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

uint32_t AllBytesAreX(const unsigned char* bytes, unsigned char match, uintptr_t byte_ct) {
  if (byte_ct < kBytesPerWord) {
    for (uint32_t uii = 0; uii != byte_ct; ++uii) {
      if (bytes[uii] != match) {
        return 0;
      }
    }
    return 1;
  }
  const uintptr_t* bytes_alias = R_CAST(const uintptr_t*, bytes);
  const uintptr_t word_match = S_CAST(uintptr_t, match) * kMask0101;
  uintptr_t word_ct_m1 = (byte_ct - 1) / kBytesPerWord;
  // todo: try movemask in AVX2 case
  for (uintptr_t widx = 0; widx != word_ct_m1; ++widx) {
    if (bytes_alias[widx] != word_match) {
      return 0;
    }
  }
  const uintptr_t last_word = *R_CAST(const uintptr_t*, &(bytes[byte_ct - kBytesPerWord]));
  if (last_word != word_match) {
    return 0;
  }
  return 1;
}

#ifdef USE_AVX2
// void CopyBitarrSubsetEx(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t bit_idx_start, uint32_t output_bit_idx_end, uintptr_t* __restrict output_bitarr) {
void CopyBitarrSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t output_bit_idx_end, uintptr_t* __restrict output_bitarr) {
  const uint32_t output_bit_idx_end_lowbits = output_bit_idx_end % kBitsPerWord;
  uintptr_t* output_bitarr_iter = output_bitarr;
  uintptr_t* output_bitarr_last = &(output_bitarr[output_bit_idx_end / kBitsPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = 0;
  while ((output_bitarr_iter != output_bitarr_last) || (write_idx_lowbits != output_bit_idx_end_lowbits)) {
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
  VecW prev_sad_result = vecw_setzero();
  const uintptr_t vec_ct_a16 = RoundDownPow2(vec_ct, 16);
  for (uintptr_t vec_idx = 0; vec_idx != vec_ct_a16; vec_idx += 16) {
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
    cnt = cnt + prev_sad_result;
    // work around high SAD latency
    prev_sad_result = PopcountVecAvx2(sixteens);
  }
  bit_vvec = &(bit_vvec[vec_ct_a16]);
  const uintptr_t remainder = vec_ct % 16;
  cnt = cnt + prev_sad_result;
  if (remainder < 12) {
    cnt = vecw_slli(cnt, 4);
    if (remainder) {
      VecW popcnt1_acc = vecw_setzero();
      VecW popcnt2_acc = vecw_setzero();
      const VecW lookup1 = vecw_setr8(4, 5, 5, 6, 5, 6, 6, 7,
                                      5, 6, 6, 7, 6, 7, 7, 8);
      const VecW lookup2 = vecw_setr8(4, 3, 3, 2, 3, 2, 2, 1,
                                      3, 2, 2, 1, 2, 1, 1, 0);

      const VecW m4 = VCONST_W(kMask0F0F);
      for (uintptr_t vec_idx = 0; vec_idx != remainder; ++vec_idx) {
        const VecW vv = bit_vvec[vec_idx];
        const VecW lo = vv & m4;
        const VecW hi = vecw_srli(vv, 4) & m4;
        popcnt1_acc = popcnt1_acc + vecw_shuffle8(lookup1, lo);
        popcnt2_acc = popcnt2_acc + vecw_shuffle8(lookup2, hi);
      }
      cnt = cnt + vecw_sad(popcnt1_acc, popcnt2_acc);
    }
  } else {
    VecW twos_a = Csa256(bit_vvec[0], bit_vvec[1], &ones);
    VecW twos_b = Csa256(bit_vvec[2], bit_vvec[3], &ones);
    VecW fours_a = Csa256(twos_a, twos_b, &twos);
    twos_a = Csa256(bit_vvec[4], bit_vvec[5], &ones);
    twos_b = Csa256(bit_vvec[6], bit_vvec[7], &ones);
    VecW fours_b = Csa256(twos_a, twos_b, &twos);
    const VecW eights_a = Csa256(fours_a, fours_b, &fours);
    twos_a = Csa256(bit_vvec[8], bit_vvec[9], &ones);
    twos_b = Csa256(bit_vvec[10], bit_vvec[11], &ones);
    fours_a = Csa256(twos_a, twos_b, &twos);
    twos_a = vecw_setzero();
    if (remainder & 2) {
      twos_a = Csa256(bit_vvec[12], bit_vvec[13], &ones);
    }
    twos_b = vecw_setzero();
    if (remainder & 1) {
      twos_b = CsaOne256(bit_vvec[remainder - 1], &ones);
    }
    fours_b = Csa256(twos_a, twos_b, &twos);
    const VecW eights_b = Csa256(fours_a, fours_b, &fours);
    const VecW sixteens = Csa256(eights_a, eights_b, &eights);
    cnt = cnt + PopcountVecAvx2(sixteens);
    cnt = vecw_slli(cnt, 4);
  }
  // Appears to be counterproductive to put multiple SAD instructions in
  // flight.
  // Compiler is smart enough that it's pointless to manually inline
  // PopcountVecAvx2.  (Tried combining the 4 SAD calls into one, didn't help.)
  cnt = cnt + vecw_slli(PopcountVecAvx2(eights), 3);
  cnt = cnt + vecw_slli(PopcountVecAvx2(fours), 2);
  cnt = cnt + vecw_slli(PopcountVecAvx2(twos), 1);
  cnt = cnt + PopcountVecAvx2(ones);
  return HsumW(cnt);
}

uintptr_t PopcountVecsAvx2Intersect(const VecW* __restrict vvec1_iter, const VecW* __restrict vvec2_iter, uintptr_t vec_ct) {
  // See popcnt_avx2() in libpopcnt.  vec_ct must be a multiple of 16.
  VecW cnt = vecw_setzero();
  VecW ones = vecw_setzero();
  VecW twos = vecw_setzero();
  VecW fours = vecw_setzero();
  VecW eights = vecw_setzero();
  for (uintptr_t vec_idx = 0; vec_idx < vec_ct; vec_idx += 16) {
    VecW twos_a = Csa256(vvec1_iter[vec_idx + 0] & vvec2_iter[vec_idx + 0], vvec1_iter[vec_idx + 1] & vvec2_iter[vec_idx + 1], &ones);
    VecW twos_b = Csa256(vvec1_iter[vec_idx + 2] & vvec2_iter[vec_idx + 2], vvec1_iter[vec_idx + 3] & vvec2_iter[vec_idx + 3], &ones);
    VecW fours_a = Csa256(twos_a, twos_b, &twos);

    twos_a = Csa256(vvec1_iter[vec_idx + 4] & vvec2_iter[vec_idx + 4], vvec1_iter[vec_idx + 5] & vvec2_iter[vec_idx + 5], &ones);
    twos_b = Csa256(vvec1_iter[vec_idx + 6] & vvec2_iter[vec_idx + 6], vvec1_iter[vec_idx + 7] & vvec2_iter[vec_idx + 7], &ones);
    VecW fours_b = Csa256(twos_a, twos_b, &twos);
    const VecW eights_a = Csa256(fours_a, fours_b, &fours);

    twos_a = Csa256(vvec1_iter[vec_idx + 8] & vvec2_iter[vec_idx + 8], vvec1_iter[vec_idx + 9] & vvec2_iter[vec_idx + 9], &ones);
    twos_b = Csa256(vvec1_iter[vec_idx + 10] & vvec2_iter[vec_idx + 10], vvec1_iter[vec_idx + 11] & vvec2_iter[vec_idx + 11], &ones);
    fours_a = Csa256(twos_a, twos_b, &twos);

    twos_a = Csa256(vvec1_iter[vec_idx + 12] & vvec2_iter[vec_idx + 12], vvec1_iter[vec_idx + 13] & vvec2_iter[vec_idx + 13], &ones);
    twos_b = Csa256(vvec1_iter[vec_idx + 14] & vvec2_iter[vec_idx + 14], vvec1_iter[vec_idx + 15] & vvec2_iter[vec_idx + 15], &ones);
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
  return HsumW(cnt);
}

uintptr_t PopcountWordsIntersect(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2_iter, uintptr_t word_ct) {
  const uintptr_t* bitvec1_end = &(bitvec1_iter[word_ct]);
  const uintptr_t block_ct = word_ct / (16 * kWordsPerVec);
  uintptr_t tot = 0;
  if (block_ct) {
    tot = PopcountVecsAvx2Intersect(R_CAST(const VecW*, bitvec1_iter), R_CAST(const VecW*, bitvec2_iter), block_ct * 16);
    bitvec1_iter = &(bitvec1_iter[block_ct * (16 * kWordsPerVec)]);
    bitvec2_iter = &(bitvec2_iter[block_ct * (16 * kWordsPerVec)]);
  }
  while (bitvec1_iter < bitvec1_end) {
    tot += PopcountWord((*bitvec1_iter++) & (*bitvec2_iter++));
  }
  return tot;
}

uintptr_t PopcountVecsAvx2Xor(const VecW* __restrict vvec1_iter, const VecW* __restrict vvec2_iter, uintptr_t vec_ct) {
  // vec_ct must be a multiple of 16.
  VecW cnt = vecw_setzero();
  VecW ones = vecw_setzero();
  VecW twos = vecw_setzero();
  VecW fours = vecw_setzero();
  VecW eights = vecw_setzero();
  for (uintptr_t vec_idx = 0; vec_idx < vec_ct; vec_idx += 16) {
    VecW twos_a = Csa256(vvec1_iter[vec_idx + 0] ^ vvec2_iter[vec_idx + 0], vvec1_iter[vec_idx + 1] ^ vvec2_iter[vec_idx + 1], &ones);
    VecW twos_b = Csa256(vvec1_iter[vec_idx + 2] ^ vvec2_iter[vec_idx + 2], vvec1_iter[vec_idx + 3] ^ vvec2_iter[vec_idx + 3], &ones);
    VecW fours_a = Csa256(twos_a, twos_b, &twos);

    twos_a = Csa256(vvec1_iter[vec_idx + 4] ^ vvec2_iter[vec_idx + 4], vvec1_iter[vec_idx + 5] ^ vvec2_iter[vec_idx + 5], &ones);
    twos_b = Csa256(vvec1_iter[vec_idx + 6] ^ vvec2_iter[vec_idx + 6], vvec1_iter[vec_idx + 7] ^ vvec2_iter[vec_idx + 7], &ones);
    VecW fours_b = Csa256(twos_a, twos_b, &twos);
    const VecW eights_a = Csa256(fours_a, fours_b, &fours);

    twos_a = Csa256(vvec1_iter[vec_idx + 8] ^ vvec2_iter[vec_idx + 8], vvec1_iter[vec_idx + 9] ^ vvec2_iter[vec_idx + 9], &ones);
    twos_b = Csa256(vvec1_iter[vec_idx + 10] ^ vvec2_iter[vec_idx + 10], vvec1_iter[vec_idx + 11] ^ vvec2_iter[vec_idx + 11], &ones);
    fours_a = Csa256(twos_a, twos_b, &twos);

    twos_a = Csa256(vvec1_iter[vec_idx + 12] ^ vvec2_iter[vec_idx + 12], vvec1_iter[vec_idx + 13] ^ vvec2_iter[vec_idx + 13], &ones);
    twos_b = Csa256(vvec1_iter[vec_idx + 14] ^ vvec2_iter[vec_idx + 14], vvec1_iter[vec_idx + 15] ^ vvec2_iter[vec_idx + 15], &ones);
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
  return HsumW(cnt);
}

uintptr_t PopcountWordsXor(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2_iter, uintptr_t word_ct) {
  const uintptr_t* bitvec1_end = &(bitvec1_iter[word_ct]);
  const uintptr_t block_ct = word_ct / (16 * kWordsPerVec);
  uintptr_t tot = 0;
  if (block_ct) {
    tot = PopcountVecsAvx2Xor(R_CAST(const VecW*, bitvec1_iter), R_CAST(const VecW*, bitvec2_iter), block_ct * 16);
    bitvec1_iter = &(bitvec1_iter[block_ct * (16 * kWordsPerVec)]);
    bitvec2_iter = &(bitvec2_iter[block_ct * (16 * kWordsPerVec)]);
  }
  while (bitvec1_iter < bitvec1_end) {
    tot += PopcountWord((*bitvec1_iter++) ^ (*bitvec2_iter++));
  }
  return tot;
}

void ExpandBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, uint32_t word_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t leading_byte_ct = 1 + (expand_sizex_m1 % kBitsPerWord) / CHAR_BIT;
  uintptr_t compact_word = SubwordLoad(compact_bitarr, leading_byte_ct) >> read_start_bit;
  const uintptr_t* compact_bitarr_iter = R_CAST(const uintptr_t*, &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]));
  uint32_t compact_idx_lowbits = read_start_bit + CHAR_BIT * (sizeof(intptr_t) - leading_byte_ct);
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
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
#  ifdef NO_UNALIGNED
#    error "Unaligned accesses in ExpandBytearr()."
#  endif
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
#  ifdef NO_UNALIGNED
#    error "Unaligned accesses in ExpandThenSubsetBytearr()."
#  endif
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
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
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
#  ifdef NO_UNALIGNED
#    error "Unaligned accesses in ExpandBytearrNested()."
#  endif
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
#  ifdef NO_UNALIGNED
#    error "Unaligned accesses in ExpandThenSubsetBytearrNested()."
#  endif
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
void CopyBitarrSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t output_bit_idx_end, uintptr_t* __restrict output_bitarr) {
  const uint32_t output_bit_idx_end_lowbits = output_bit_idx_end % kBitsPerWord;
  uintptr_t* output_bitarr_iter = output_bitarr;
  uintptr_t* output_bitarr_last = &(output_bitarr[output_bit_idx_end / kBitsPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = 0;
  while ((output_bitarr_iter != output_bitarr_last) || (write_idx_lowbits != output_bit_idx_end_lowbits)) {
    uintptr_t cur_mask_word;
    // sparse subset_mask optimization
    // guaranteed to terminate since there's at least one more set bit
    do {
      cur_mask_word = subset_mask[++read_widx];
    } while (!cur_mask_word);
    uintptr_t cur_masked_input_word = raw_bitarr[read_widx] & cur_mask_word;
    const uint32_t cur_mask_popcount = PopcountWord(cur_mask_word);
    uintptr_t subsetted_input_word = 0;
    while (cur_masked_input_word) {
      const uintptr_t mask_word_high = (cur_mask_word | (cur_masked_input_word ^ (cur_masked_input_word - 1))) + 1;
      if (!mask_word_high) {
        subsetted_input_word |= cur_masked_input_word >> (kBitsPerWord - cur_mask_popcount);
        break;
      }
      const uint32_t cur_read_end = ctzw(mask_word_high);
      const uintptr_t bits_to_copy = cur_masked_input_word & (~mask_word_high);
      cur_masked_input_word ^= bits_to_copy;
      const uint32_t cur_write_end = PopcountWord(cur_mask_word & (~mask_word_high));
      subsetted_input_word |= bits_to_copy >> (cur_read_end - cur_write_end);
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
uintptr_t PopcountVecsNoAvx2(const VecW* bit_vvec, uintptr_t vec_ct) {
  // popcounts vptr[0..(vec_ct-1)].  Assumes vec_ct is a multiple of 3 (0 ok).
  assert(!(vec_ct % 3));
  const VecW m0 = vecw_setzero();
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* bit_vvec_iter = bit_vvec;
  VecW prev_sad_result = vecw_setzero();
  VecW acc = vecw_setzero();
  uintptr_t cur_incr = 30;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 30) {
      if (!vec_ct) {
        acc = acc + prev_sad_result;
        return HsumW(acc);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc = vecw_setzero();
    const VecW* bit_vvec_stop = &(bit_vvec_iter[cur_incr]);
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
      inner_acc = inner_acc + (count1 & m4) + (vecw_srli(count1, 4) & m4);
    } while (bit_vvec_iter < bit_vvec_stop);
    // _mm_sad_epu8() has better throughput than the previous method of
    // horizontal-summing the bytes in inner_acc, by enough to compensate for
    // the loop length being reduced from 30 to 15 vectors, but it has high
    // latency.  We work around that by waiting till the end of the next full
    // loop iteration to actually use the SAD result.
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

static inline uintptr_t PopcountVecsNoAvx2Intersect(const VecW* __restrict vvec1_iter, const VecW* __restrict vvec2_iter, uintptr_t vec_ct) {
  // popcounts vvec1 AND vvec2[0..(ct-1)].  ct is a multiple of 3.
  assert(!(vec_ct % 3));
  const VecW m0 = vecw_setzero();
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  VecW prev_sad_result = vecw_setzero();
  VecW acc = vecw_setzero();
  uintptr_t cur_incr = 30;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 30) {
      if (!vec_ct) {
        acc = acc + prev_sad_result;
        return HsumW(acc);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc = vecw_setzero();
    const VecW* vvec1_stop = &(vvec1_iter[cur_incr]);
    do {
      VecW count1 = (*vvec1_iter++) & (*vvec2_iter++);
      VecW count2 = (*vvec1_iter++) & (*vvec2_iter++);
      VecW half1 = (*vvec1_iter++) & (*vvec2_iter++);
      const VecW half2 = vecw_srli(half1, 1) & m1;
      half1 = half1 & m1;
      count1 = count1 - (vecw_srli(count1, 1) & m1);
      count2 = count2 - (vecw_srli(count2, 1) & m1);
      count1 = count1 + half1;
      count2 = count2 + half2;
      count1 = (count1 & m2) + (vecw_srli(count1, 2) & m2);
      count1 = count1 + (count2 & m2) + (vecw_srli(count2, 2) & m2);
      inner_acc = inner_acc + (count1 & m4) + (vecw_srli(count1, 4) & m4);
    } while (vvec1_iter < vvec1_stop);
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

uintptr_t PopcountWordsIntersect(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2_iter, uintptr_t word_ct) {
  uintptr_t tot = 0;
  const uintptr_t* bitvec1_end = &(bitvec1_iter[word_ct]);
  const uintptr_t trivec_ct = word_ct / (3 * kWordsPerVec);
  tot += PopcountVecsNoAvx2Intersect(R_CAST(const VecW*, bitvec1_iter), R_CAST(const VecW*, bitvec2_iter), trivec_ct * 3);
  bitvec1_iter = &(bitvec1_iter[trivec_ct * (3 * kWordsPerVec)]);
  bitvec2_iter = &(bitvec2_iter[trivec_ct * (3 * kWordsPerVec)]);
  while (bitvec1_iter < bitvec1_end) {
    tot += PopcountWord((*bitvec1_iter++) & (*bitvec2_iter++));
  }
  return tot;
}

static inline uintptr_t PopcountVecsNoAvx2Xor(const VecW* __restrict vvec1_iter, const VecW* __restrict vvec2_iter, uintptr_t vec_ct) {
  // popcounts vvec1 XOR vvec2[0..(ct-1)].  ct is a multiple of 3.
  assert(!(vec_ct % 3));
  const VecW m0 = vecw_setzero();
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  VecW prev_sad_result = vecw_setzero();
  VecW acc = vecw_setzero();
  uintptr_t cur_incr = 30;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 30) {
      if (!vec_ct) {
        acc = acc + prev_sad_result;
        return HsumW(acc);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc = vecw_setzero();
    const VecW* vvec1_stop = &(vvec1_iter[cur_incr]);
    do {
      VecW count1 = (*vvec1_iter++) ^ (*vvec2_iter++);
      VecW count2 = (*vvec1_iter++) ^ (*vvec2_iter++);
      VecW half1 = (*vvec1_iter++) ^ (*vvec2_iter++);
      const VecW half2 = vecw_srli(half1, 1) & m1;
      half1 = half1 & m1;
      count1 = count1 - (vecw_srli(count1, 1) & m1);
      count2 = count2 - (vecw_srli(count2, 1) & m1);
      count1 = count1 + half1;
      count2 = count2 + half2;
      count1 = (count1 & m2) + (vecw_srli(count1, 2) & m2);
      count1 = count1 + (count2 & m2) + (vecw_srli(count2, 2) & m2);
      inner_acc = inner_acc + (count1 & m4) + (vecw_srli(count1, 4) & m4);
    } while (vvec1_iter < vvec1_stop);
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

uintptr_t PopcountWordsXor(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2_iter, uintptr_t word_ct) {
  uintptr_t tot = 0;
  const uintptr_t* bitvec1_end = &(bitvec1_iter[word_ct]);
  const uintptr_t trivec_ct = word_ct / (3 * kWordsPerVec);
  tot += PopcountVecsNoAvx2Xor(R_CAST(const VecW*, bitvec1_iter), R_CAST(const VecW*, bitvec2_iter), trivec_ct * 3);
  bitvec1_iter = &(bitvec1_iter[trivec_ct * (3 * kWordsPerVec)]);
  bitvec2_iter = &(bitvec2_iter[trivec_ct * (3 * kWordsPerVec)]);
  while (bitvec1_iter < bitvec1_end) {
    tot += PopcountWord((*bitvec1_iter++) ^ (*bitvec2_iter++));
  }
  return tot;
}

void ExpandBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, uint32_t word_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  ZeroWArr(word_ct, target);
  const uintptr_t* compact_bitarr_alias = S_CAST(const uintptr_t*, compact_bitarr);
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t compact_widx_last = expand_sizex_m1 / kBitsPerWord;
  uint32_t compact_idx_lowbits = read_start_bit;
  uint32_t loop_len = kBitsPerWord;
  uintptr_t write_widx = 0;
  uintptr_t expand_mask_bits = expand_mask[0];
  for (uint32_t compact_widx = 0; ; ++compact_widx) {
    uintptr_t compact_word;
    if (compact_widx >= compact_widx_last) {
      if (compact_widx > compact_widx_last) {
        return;
      }
      loop_len = 1 + (expand_sizex_m1 % kBitsPerWord);
      // avoid possible segfault
      compact_word = SubwordLoad(&(compact_bitarr_alias[compact_widx]), DivUp(loop_len, CHAR_BIT));
    } else {
#  ifdef NO_UNALIGNED
#    error "Unaligned accesses in ExpandBytearr()."
#  endif
      compact_word = compact_bitarr_alias[compact_widx];
    }
    for (; compact_idx_lowbits != loop_len; ++compact_idx_lowbits) {
      const uintptr_t lowbit = BitIter1y(expand_mask, &write_widx, &expand_mask_bits);
      // bugfix: can't just use (compact_word & 1) and compact_word >>= 1,
      // since we may skip the first bit on the first loop iteration
      if ((compact_word >> compact_idx_lowbits) & 1) {
        target[write_widx] |= lowbit;
      }
    }
    compact_idx_lowbits = 0;
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
#  ifdef NO_UNALIGNED
#    error "Unaligned accesses in ExpandThenSubsetBytearr()."
#  endif
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
  // can allow compact_idx_lowbits to be initialized to nonzero
  uint32_t loop_len = kBitsPerWord;
  uintptr_t write_widx = 0;
  uintptr_t top_expand_mask_bits = top_expand_mask[0];
  for (uint32_t compact_widx = 0; ; ++compact_widx) {
    uintptr_t compact_word;
    if (compact_widx >= compact_widx_last) {
      if (compact_widx > compact_widx_last) {
        return;
      }
      loop_len = 1 + (mid_popcount_m1 % kBitsPerWord);
      // avoid possible segfault
      compact_word = SubwordLoad(&(compact_bitarr_alias[compact_widx]), DivUp(loop_len, CHAR_BIT));
    } else {
#ifdef NO_UNALIGNED
#  error "Unaligned accesses in ExpandBytearrNested()."
#endif
      compact_word = compact_bitarr_alias[compact_widx];
    }
    for (uint32_t compact_idx_lowbits = 0; compact_idx_lowbits != loop_len; ++mid_idx) {
      const uintptr_t lowbit = BitIter1y(top_expand_mask, &write_widx, &top_expand_mask_bits);
      if (IsSet(mid_bitarr, mid_idx)) {
        mid_target[write_widx] |= lowbit;
        compact_target[write_widx] |= lowbit * (compact_word & 1);
        compact_word >>= 1;
        ++compact_idx_lowbits;
      }
    }
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
#ifdef NO_UNALIGNED
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
uintptr_t PopcountBytes(const void* bitarr, uintptr_t byte_ct) {
  const unsigned char* bitarr_uc = S_CAST(const unsigned char*, bitarr);
  const uint32_t lead_byte_ct = (-R_CAST(uintptr_t, bitarr_uc)) % kBytesPerVec;
  uintptr_t tot = 0;
  const uintptr_t* bitarr_iter;
  uint32_t trail_byte_ct;
  // bugfix: had wrong condition here
  if (byte_ct >= lead_byte_ct) {
#ifdef __LP64__
    const uint32_t word_rem = lead_byte_ct % kBytesPerWord;
    if (word_rem) {
      tot = PopcountWord(ProperSubwordLoad(bitarr_uc, word_rem));
    }
    bitarr_iter = R_CAST(const uintptr_t*, &(bitarr_uc[word_rem]));
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
      tot = PopcountWord(ProperSubwordLoad(bitarr_uc, lead_byte_ct));
    }
    bitarr_iter = R_CAST(const uintptr_t*, &(bitarr_uc[lead_byte_ct]));
#endif
    byte_ct -= lead_byte_ct;
    const uintptr_t word_ct = byte_ct / kBytesPerWord;
    // vec-alignment required here
    tot += PopcountWords(bitarr_iter, word_ct);
    bitarr_iter = &(bitarr_iter[word_ct]);
    trail_byte_ct = byte_ct % kBytesPerWord;
  } else {
    bitarr_iter = R_CAST(const uintptr_t*, bitarr_uc);
    // this may still be >= kBytesPerWord, so can't remove loop
    trail_byte_ct = byte_ct;
  }
  for (uint32_t bytes_remaining = trail_byte_ct; ; ) {
    uintptr_t cur_word;
    if (bytes_remaining < kBytesPerWord) {
      if (!bytes_remaining) {
        return tot;
      }
      cur_word = ProperSubwordLoad(bitarr_iter, bytes_remaining);
      bytes_remaining = 0;
    } else {
      cur_word = *bitarr_iter++;
      bytes_remaining -= kBytesPerWord;
    }
    tot += PopcountWord(cur_word);
  }
}

uintptr_t PopcountBytesMasked(const void* bitarr, const uintptr_t* mask_arr, uintptr_t byte_ct) {
  // todo: try modifying PopcountWordsIntersect() to use unaligned load
  // instructions; then, if there is no performance penalty, try modifying this
  // main loop to call it.
  const uintptr_t word_ct = byte_ct / kBytesPerWord;
#ifdef USE_SSE42
  const uintptr_t* bitarr_w = S_CAST(const uintptr_t*, bitarr);
  uintptr_t tot = 0;
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    tot += PopcountWord(bitarr_w[widx] & mask_arr[widx]);
  }
  const uint32_t trail_byte_ct = byte_ct % kBytesPerWord;
  if (trail_byte_ct) {
    uintptr_t cur_word = ProperSubwordLoad(&(bitarr_w[word_ct]), trail_byte_ct);
    tot += PopcountWord(cur_word & mask_arr[word_ct]);
  }
  return tot;
#else
  const uintptr_t* bitarr_iter = S_CAST(const uintptr_t*, bitarr);
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
  for (uint32_t trail_byte_ct = byte_ct - (mainblock_word_ct * kBytesPerWord); ; ) {
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
  for (uint32_t widx = 0; widx != word_ct_m1; ++widx) {
    cumulative_popcounts[widx] = cur_sum;
    cur_sum += PopcountWord(subset_mask[widx]);
  }
  cumulative_popcounts[word_ct_m1] = cur_sum;
}

void UidxsToIdxs(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, const uintptr_t idx_list_len, uint32_t* idx_list) {
  uint32_t* idx_list_end = &(idx_list[idx_list_len]);
  for (uint32_t* idx_list_iter = idx_list; idx_list_iter != idx_list_end; ++idx_list_iter) {
    *idx_list_iter = RawToSubsettedPos(subset_mask, subset_cumulative_popcounts, *idx_list_iter);
  }
}

void Expand1bitTo8(const void* __restrict bytearr, uint32_t input_bit_ct, uint32_t incr, uintptr_t* __restrict dst) {
  const unsigned char* bytearr_uc = S_CAST(const unsigned char*, bytearr);
  const uint32_t input_bit_ct_plus = input_bit_ct + kBytesPerWord - 1;
#ifdef USE_SSE42
  const uint32_t input_byte_ct = input_bit_ct_plus / 8;
  const uint32_t fullvec_ct = input_byte_ct / (kBytesPerVec / 8);
  uint32_t byte_idx = 0;
  if (fullvec_ct) {
    const Vec8thUint* bytearr_alias = R_CAST(const Vec8thUint*, bytearr);
#  ifdef USE_AVX2
    const VecUc byte_gather = R_CAST(VecUc, _mm256_setr_epi64x(0, kMask0101, 2 * kMask0101, 3 * kMask0101));
    const VecUc bit_mask = R_CAST(VecUc, _mm256_set1_epi64x(0x7fbfdfeff7fbfdfeLL));
#  else
    const VecUc byte_gather = R_CAST(VecUc, _mm_setr_epi32(0, 0, 0x01010101, 0x01010101));
    const VecUc bit_mask = R_CAST(VecUc, _mm_set1_epi64x(0x7fbfdfeff7fbfdfeLL));
#  endif
    const VecUc all1 = vecuc_set1(255);
    const VecUc subfrom = vecuc_set1(incr);
    VecUc* dst_alias = R_CAST(VecUc*, dst);
    for (uint32_t vec_idx = 0; vec_idx != fullvec_ct; ++vec_idx) {
#  ifdef USE_AVX2
      VecUc vmask = R_CAST(VecUc, _mm256_set1_epi32(bytearr_alias[vec_idx]));
#  else
      VecUc vmask = R_CAST(VecUc, _mm_set1_epi16(bytearr_alias[vec_idx]));
#  endif
      vmask = vecuc_shuffle8(vmask, byte_gather);
      vmask = vmask | bit_mask;
      vmask = (vmask == all1);
      const VecUc result = subfrom - vmask;
      vecuc_storeu(&(dst_alias[vec_idx]), result);
    }
    byte_idx = fullvec_ct * (kBytesPerVec / 8);
  }
  const uintptr_t incr_word = incr * kMask0101;
  for (; byte_idx != input_byte_ct; ++byte_idx) {
    const uintptr_t input_byte = bytearr_uc[byte_idx];
#  ifdef USE_AVX2
    const uintptr_t input_byte_scatter = _pdep_u64(input_byte, kMask0101);
#  else
    const uintptr_t input_byte_scatter = (((input_byte & 0xfe) * 0x2040810204080LLU) & kMask0101) | (input_byte & 1);
#  endif
    dst[byte_idx] = incr_word + input_byte_scatter;
  }
#else
  const uintptr_t incr_word = incr * kMask0101;
#  ifdef __LP64__
  const uint32_t input_byte_ct = input_bit_ct_plus / 8;
  for (uint32_t uii = 0; uii != input_byte_ct; ++uii) {
    // this operation maps binary hgfedcba to h0000000g0000000f...
    //                                        ^       ^       ^
    //                                        |       |       |
    //                                       56      48      40
    // 1. (cur_variant_include_word & 0xfe) gives us hgfedcb0; necessary to
    //    avoid carryover.
    // 2. multiply by the number with bits 7, 14, 21, ..., 49 set, to get
    //    hgfedcbhgfedcbhgf...
    //    ^       ^       ^
    //    |       |       |
    //   56      48      40
    // 3. mask out all but bits 8, 16, 24, ..., 56
    // todo: test if this actually beats the per-character loop...
    const uintptr_t input_byte = bytearr_uc[uii];
    const uintptr_t input_byte_scatter = (((input_byte & 0xfe) * 0x2040810204080LLU) & kMask0101) | (input_byte & 1);
    dst[uii] = incr_word + input_byte_scatter;
  }
#  else
  const uint32_t fullbyte_ct = input_bit_ct_plus / 8;
  for (uint32_t uii = 0; uii != fullbyte_ct; ++uii) {
    // dcba -> d0000000c0000000b0000000a
    const uintptr_t input_byte = bytearr_uc[uii];
    uintptr_t input_byte_scatter = ((input_byte & 0xf) * 0x204081) & kMask0101;
    dst[2 * uii] = incr_word + input_byte_scatter;
    input_byte_scatter = ((input_byte >> 4) * 0x204081) & kMask0101;
    dst[2 * uii + 1] = incr_word + input_byte_scatter;
  }
  if (input_bit_ct_plus & 4) {
    uintptr_t input_byte = bytearr_uc[fullbyte_ct];
    // input_bit_ct mod 8 in 1..4, so high bits zeroed out
    uintptr_t input_byte_scatter = (input_byte * 0x204081) & kMask0101;
    dst[2 * fullbyte_ct] = incr_word + input_byte_scatter;
  }
#  endif
#endif
}

void Expand1bitTo16(const void* __restrict bytearr, uint32_t input_bit_ct, uint32_t incr, uintptr_t* __restrict dst) {
  const unsigned char* bytearr_uc = S_CAST(const unsigned char*, bytearr);
#ifdef USE_SSE42
  const uint32_t input_nybble_ct = DivUp(input_bit_ct, 4);
  const uint32_t fullvec_ct = input_nybble_ct / (kBytesPerVec / 8);
  uint32_t byte_idx = 0;
  if (fullvec_ct) {
    const Vec16thUint* bytearr_alias = R_CAST(const Vec16thUint*, bytearr);
#  ifdef USE_AVX2
    const VecU16 byte_gather = R_CAST(VecU16, _mm256_setr_epi64x(0, 0, kMask0101, kMask0101));
    const VecU16 bit_mask = R_CAST(VecU16, _mm256_set_epi32(0xff7fffbfU, 0xffdfffefU, 0xfff7fffbU, 0xfffdfffeU, 0xff7fffbfU, 0xffdfffefU, 0xfff7fffbU, 0xfffdfffeU));
#  else
    const VecU16 byte_gather = VCONST_S(0);
    const VecU16 bit_mask = R_CAST(VecU16, _mm_set_epi32(0xff7fffbfU, 0xffdfffefU, 0xfff7fffbU, 0xfffdfffeU));
#  endif
    const VecU16 all1 = VCONST_S(0xffff);
    const VecU16 subfrom = vecu16_set1(incr);
    VecU16* dst_alias = R_CAST(VecU16*, dst);
    // todo: check whether this is actually any better than the non-vectorized
    // loop
    for (uint32_t vec_idx = 0; vec_idx != fullvec_ct; ++vec_idx) {
#  ifdef USE_AVX2
      VecU16 vmask = R_CAST(VecU16, _mm256_set1_epi16(bytearr_alias[vec_idx]));
#  else
      VecU16 vmask = R_CAST(VecU16, _mm_set1_epi8(bytearr_alias[vec_idx]));
#  endif
      vmask = vecu16_shuffle8(vmask, byte_gather);
      vmask = vmask | bit_mask;
      vmask = (vmask == all1);
      const VecU16 result = subfrom - vmask;
      vecu16_storeu(&(dst_alias[vec_idx]), result);
    }
    byte_idx = fullvec_ct * (kBytesPerVec / 16);
  }
  const uintptr_t incr_word = incr * kMask0001;
  const uint32_t fullbyte_ct = input_nybble_ct / 2;
  for (; byte_idx != fullbyte_ct; ++byte_idx) {
    const uintptr_t input_byte = bytearr_uc[byte_idx];
    const uintptr_t input_byte_scatter = input_byte * 0x200040008001LLU;
    const uintptr_t write0 = input_byte_scatter & kMask0001;
    const uintptr_t write1 = (input_byte_scatter >> 4) & kMask0001;
    dst[2 * byte_idx] = incr_word + write0;
    dst[2 * byte_idx + 1] = incr_word + write1;
  }
  if (input_nybble_ct % 2) {
    const uintptr_t input_byte = bytearr_uc[byte_idx];
    const uintptr_t write0 = (input_byte * 0x200040008001LLU) & kMask0001;
    dst[input_nybble_ct - 1] = incr_word + write0;
  }
#else
  const uintptr_t incr_word = incr * kMask0001;
#  ifdef __LP64__
  const uint32_t input_nybble_ct = DivUp(input_bit_ct, 4);
  const uint32_t fullbyte_ct = input_nybble_ct / 2;
  for (uint32_t uii = 0; uii != fullbyte_ct; ++uii) {
    const uintptr_t input_byte = bytearr_uc[uii];
    const uintptr_t input_byte_scatter = input_byte * 0x200040008001LLU;
    const uintptr_t write0 = input_byte_scatter & kMask0001;
    const uintptr_t write1 = (input_byte_scatter >> 4) & kMask0001;
    dst[2 * uii] = incr_word + write0;
    dst[2 * uii + 1] = incr_word + write1;
  }
  if (input_nybble_ct % 2) {
    const uintptr_t input_byte = bytearr_uc[fullbyte_ct];
    const uintptr_t write0 = (input_byte * 0x200040008001LLU) & kMask0001;
    dst[input_nybble_ct - 1] = incr_word + write0;
  }
#  else
  const uint32_t fullbyte_ct = input_bit_ct / 8;
  for (uint32_t uii = 0; uii != fullbyte_ct; ++uii) {
    uintptr_t input_byte = bytearr_uc[uii];
    const uintptr_t input_byte_scatter = input_byte * 0x8001;
    dst[4 * uii] = (input_byte_scatter & kMask0001) + incr_word;
    dst[4 * uii + 1] = ((input_byte_scatter >> 2) & kMask0001) + incr_word;
    dst[4 * uii + 2] = ((input_byte_scatter >> 4) & kMask0001) + incr_word;
    dst[4 * uii + 3] = ((input_byte_scatter >> 6) & kMask0001) + incr_word;
  }
  const uint32_t remainder = input_bit_ct % 8;
  if (remainder) {
    uintptr_t input_byte = bytearr_uc[fullbyte_ct];
    uint16_t* dst_alias = R_CAST(uint16_t*, &(dst[4 * fullbyte_ct]));
    for (uint32_t uii = 0; uii < remainder; ++uii) {
      dst_alias[uii] = (input_byte & 1) + incr;
      input_byte = input_byte >> 1;
    }
  }
#  endif
#endif
}

static_assert(kPglBitTransposeBatch == S_CAST(uint32_t, kBitsPerCacheline), "TransposeBitblock64() needs to be updated.");
#ifdef __LP64__
void TransposeBitblock64(const uintptr_t* read_iter, uintptr_t read_ul_stride, uintptr_t write_ul_stride, uint32_t read_row_ct, uint32_t write_row_ct, uintptr_t* write_iter, VecW* __restrict buf0, VecW* __restrict buf1) {
  // We need to perform the equivalent of 9 shuffles (assuming a full-size
  // 512x512 bitblock).
  // The first shuffles are performed by the ingestion loop: we write the first
  // word from every row to buf0, then the second word from every row, etc.,
  // yielding
  //   (0,0) ...   (0,63)  (1,0) ...   (1,63)  (2,0) ...   (511,63)
  //   (0,64) ...  (0,127) (1,64) ...  (1,127) (2,64) ...  (511,127)
  //   ...
  //   (0,448) ... (0,511) (1,448) ... (1,511) (2,448) ... (511,511)
  // in terms of the original bit positions.
  // Since each input row has 8 words, this amounts to 3 shuffles.
  //
  // The second step writes
  //   (0,0) (0,1) ... (0,7)   (1,0) (1,1) ... (1,7) ...   (511,7)
  //   (0,8) (0,9) ... (0,15)  (1,8) (1,9) ... (1,15) ...  (511,15)
  //   ...
  //   (0,504) ...     (0,511) (1,504) ...     (1,511) ... (511,511)
  // to buf1, performing the equivalent of 3 shuffles, and the third step
  // finishes the transpose using movemask.
  //
  // buf0 and buf1 must both be 32KiB vector-aligned buffers.

  const uint32_t buf0_row_ct = DivUp(write_row_ct, 64);
  {
    uintptr_t* buf0_ul = R_CAST(uintptr_t*, buf0);
    const uint32_t zfill_ct = (-read_row_ct) & 63;
    for (uint32_t bidx = 0; bidx != buf0_row_ct; ++bidx) {
      const uintptr_t* read_iter_tmp = &(read_iter[bidx]);
      uintptr_t* buf0_row_start = &(buf0_ul[512 * bidx]);
      for (uint32_t uii = 0; uii != read_row_ct; ++uii) {
        buf0_row_start[uii] = *read_iter_tmp;
        read_iter_tmp = &(read_iter_tmp[read_ul_stride]);
      }
      // This is a simple way of fulfilling the trailing-zero part of the
      // function contract.
      // (   buf0 rows zeroed out to 512 bytes
      //  -> buf1 rows zeroed out to 64 bytes
      //  -> output rows zeroed out to 8 bytes)
      ZeroWArr(zfill_ct, &(buf0_row_start[read_row_ct]));
    }
  }
  // Each width-unit corresponds to 64 input rows.
  const uint32_t buf_row_xwidth = DivUp(read_row_ct, 64);
  {
    const VecW* buf0_read_iter = buf0;
    uintptr_t* write_iter0 = R_CAST(uintptr_t*, buf1);
#  ifdef USE_SSE42
    const VecW gather_u16s = vecw_setr8(0, 8, 1, 9, 2, 10, 3, 11,
                                        4, 12, 5, 13, 6, 14, 7, 15);
#    ifdef USE_AVX2
    const VecW gather_u32s = vecw_setr8(0, 1, 8, 9, 2, 3, 10, 11,
                                        4, 5, 12, 13, 6, 7, 14, 15);
#    endif
#  else
    const VecW m8 = VCONST_W(kMask00FF);
#  endif
    const uint32_t buf0_row_clwidth = buf_row_xwidth * 8;
    for (uint32_t bidx = 0; bidx != buf0_row_ct; ++bidx) {
      uintptr_t* write_iter1 = &(write_iter0[64]);
      uintptr_t* write_iter2 = &(write_iter1[64]);
      uintptr_t* write_iter3 = &(write_iter2[64]);
      uintptr_t* write_iter4 = &(write_iter3[64]);
      uintptr_t* write_iter5 = &(write_iter4[64]);
      uintptr_t* write_iter6 = &(write_iter5[64]);
      uintptr_t* write_iter7 = &(write_iter6[64]);
      for (uint32_t clidx = 0; clidx != buf0_row_clwidth; ++clidx) {
#  ifdef USE_AVX2
        VecW loader0 = buf0_read_iter[clidx * 2];
        VecW loader1 = buf0_read_iter[clidx * 2 + 1];
        //    (0,0) (0,1) ... (0,7) (1,0) (1,1) ... (1,7) (2,0) ... (3,7)
        // -> (0,0) (1,0) (0,1) (1,1) (0,2) .... (1,7) (2,0) (3,0) (2,1) ...
        loader0 = vecw_shuffle8(loader0, gather_u16s);
        loader1 = vecw_shuffle8(loader1, gather_u16s);
        // -> (0,0) (1,0) (0,1) (1,1) (0,2) (1,2) (0,3) (1,3) (2,0) (3,0) ...
        VecW vec_lo = vecw_permute0xd8_if_avx2(loader0);
        VecW vec_hi = vecw_permute0xd8_if_avx2(loader1);
        // -> (0,0) (1,0) (2,0) (3,0) (0,1) (1,1) (2,1) (3,1) (0,2) ...
        vec_lo = vecw_shuffle8(vec_lo, gather_u32s);
        // -> (4,0) (5,0) (6,0) (7,0) (4,1) (5,1) (6,1) (7,1) (4,2) ...
        vec_hi = vecw_shuffle8(vec_hi, gather_u32s);
        const VecW final0145 = vecw_unpacklo32(vec_lo, vec_hi);
        const VecW final2367 = vecw_unpackhi32(vec_lo, vec_hi);
        write_iter0[clidx] = vecw_extract64_0(final0145);
        write_iter1[clidx] = vecw_extract64_1(final0145);
        write_iter2[clidx] = vecw_extract64_0(final2367);
        write_iter3[clidx] = vecw_extract64_1(final2367);
        write_iter4[clidx] = vecw_extract64_2(final0145);
        write_iter5[clidx] = vecw_extract64_3(final0145);
        write_iter6[clidx] = vecw_extract64_2(final2367);
        write_iter7[clidx] = vecw_extract64_3(final2367);
#  else  // !USE_AVX2
        VecW loader0 = buf0_read_iter[clidx * 4];
        VecW loader1 = buf0_read_iter[clidx * 4 + 1];
        VecW loader2 = buf0_read_iter[clidx * 4 + 2];
        VecW loader3 = buf0_read_iter[clidx * 4 + 3];
        //    (0,0) (0,1) ... (0,7) (1,0) (1,1) ... (1,7)
        // -> (0,0) (1,0) (0,1) (1,1) (0,2) ... (1,7)
#    ifdef USE_SSE42
        loader0 = vecw_shuffle8(loader0, gather_u16s);
        loader1 = vecw_shuffle8(loader1, gather_u16s);
        loader2 = vecw_shuffle8(loader2, gather_u16s);
        loader3 = vecw_shuffle8(loader3, gather_u16s);
#    else
        VecW tmp_lo = vecw_unpacklo8(loader0, loader1);
        VecW tmp_hi = vecw_unpackhi8(loader0, loader1);
        loader0 = vecw_blendv(vecw_slli(tmp_hi, 8), tmp_lo, m8);
        loader1 = vecw_blendv(tmp_hi, vecw_srli(tmp_lo, 8), m8);
        tmp_lo = vecw_unpacklo8(loader2, loader3);
        tmp_hi = vecw_unpackhi8(loader2, loader3);
        loader2 = vecw_blendv(vecw_slli(tmp_hi, 8), tmp_lo, m8);
        loader3 = vecw_blendv(tmp_hi, vecw_srli(tmp_lo, 8), m8);
#    endif
        // -> (0,0) (1,0) (2,0) (3,0) (0,1) (1,1) (2,1) (3,1) (0,2) ...
        const VecW lo_0123 = vecw_unpacklo16(loader0, loader1);
        // -> (0,4) (1,4) (2,4) (3,4) (0,5) (1,5) (2,5) (3,5) (0,6) ...
        const VecW lo_4567 = vecw_unpackhi16(loader0, loader1);
        const VecW hi_0123 = vecw_unpacklo16(loader2, loader3);
        const VecW hi_4567 = vecw_unpackhi16(loader2, loader3);

        VecW final01 = vecw_unpacklo32(lo_0123, hi_0123);
        VecW final23 = vecw_unpackhi32(lo_0123, hi_0123);
        VecW final45 = vecw_unpacklo32(lo_4567, hi_4567);
        VecW final67 = vecw_unpackhi32(lo_4567, hi_4567);
        write_iter0[clidx] = vecw_extract64_0(final01);
        write_iter1[clidx] = vecw_extract64_1(final01);
        write_iter2[clidx] = vecw_extract64_0(final23);
        write_iter3[clidx] = vecw_extract64_1(final23);
        write_iter4[clidx] = vecw_extract64_0(final45);
        write_iter5[clidx] = vecw_extract64_1(final45);
        write_iter6[clidx] = vecw_extract64_0(final67);
        write_iter7[clidx] = vecw_extract64_1(final67);
#  endif  // !USE_AVX2
      }
      buf0_read_iter = &(buf0_read_iter[512 / kWordsPerVec]);
      write_iter0 = &(write_iter7[64]);
    }
  }
  const VecW* buf1_read_iter = buf1;
  const uint32_t write_v8ui_stride = kVec8thUintPerWord * write_ul_stride;
  const uint32_t buf1_fullrow_ct = write_row_ct / 8;
  const uint32_t buf1_row_vecwidth = buf_row_xwidth * (8 / kWordsPerVec);
  Vec8thUint* write_iter0 = R_CAST(Vec8thUint*, write_iter);
  for (uint32_t bidx = 0; bidx != buf1_fullrow_ct; ++bidx) {
    Vec8thUint* write_iter1 = &(write_iter0[write_v8ui_stride]);
    Vec8thUint* write_iter2 = &(write_iter1[write_v8ui_stride]);
    Vec8thUint* write_iter3 = &(write_iter2[write_v8ui_stride]);
    Vec8thUint* write_iter4 = &(write_iter3[write_v8ui_stride]);
    Vec8thUint* write_iter5 = &(write_iter4[write_v8ui_stride]);
    Vec8thUint* write_iter6 = &(write_iter5[write_v8ui_stride]);
    Vec8thUint* write_iter7 = &(write_iter6[write_v8ui_stride]);
    for (uint32_t vidx = 0; vidx != buf1_row_vecwidth; ++vidx) {
      VecW loader = buf1_read_iter[vidx];
      write_iter7[vidx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      write_iter6[vidx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      write_iter5[vidx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      write_iter4[vidx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      write_iter3[vidx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      write_iter2[vidx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      write_iter1[vidx] = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      write_iter0[vidx] = vecw_movemask(loader);
    }
    buf1_read_iter = &(buf1_read_iter[64 / kWordsPerVec]);
    write_iter0 = &(write_iter7[write_v8ui_stride]);
  }
  const uint32_t row_ct_rem = write_row_ct % 8;
  if (!row_ct_rem) {
    return;
  }
  const uint32_t lshift = 8 - row_ct_rem;
  Vec8thUint* write_iter_last = &(write_iter0[write_v8ui_stride * (row_ct_rem - 1)]);
  for (uint32_t vidx = 0; vidx != buf1_row_vecwidth; ++vidx) {
    VecW loader = buf1_read_iter[vidx];
    loader = vecw_slli_variable_ct(loader, lshift);
    Vec8thUint* inner_write_iter = &(write_iter_last[vidx]);
    for (uint32_t uii = 0; uii != row_ct_rem; ++uii) {
      *inner_write_iter = vecw_movemask(loader);
      loader = vecw_slli(loader, 1);
      inner_write_iter -= write_v8ui_stride;
    }
  }
}
#else  // !__LP64__
static_assert(kWordsPerVec == 1, "TransposeBitblock32() needs to be updated.");
void TransposeBitblock32(const uintptr_t* read_iter, uintptr_t read_ul_stride, uintptr_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecW* __restrict buf0, VecW* __restrict buf1) {
  // buf must be vector-aligned and have size 64k
  const uint32_t initial_read_byte_ct = DivUp(write_batch_size, CHAR_BIT);
  // fold the first 6 shuffles into the initial ingestion loop
  const unsigned char* initial_read_iter = R_CAST(const unsigned char*, read_iter);
  const unsigned char* initial_read_end = &(initial_read_iter[initial_read_byte_ct]);
  unsigned char* initial_target_iter = R_CAST(unsigned char*, buf0);
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t read_batch_rem = kBitsPerCacheline - read_batch_size;
  for (; initial_read_iter != initial_read_end; ++initial_read_iter) {
    const unsigned char* read_iter_tmp = initial_read_iter;
    for (uint32_t ujj = 0; ujj != read_batch_size; ++ujj) {
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
  for (uint32_t uii = 0; uii != initial_read_byte_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 4]);
    for (uint32_t ujj = 0; ujj != first_inner_loop_iter_ct; ++ujj) {
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
  for (uint32_t uii = 0; uii != second_outer_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 2]);
    for (uint32_t ujj = 0; ujj != second_inner_loop_iter_ct; ++ujj) {
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
  for (uint32_t uii = 0; uii != last_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    for (uint32_t ujj = 0; ujj != write_word_ct; ++ujj) {
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

#ifdef __LP64__
void TransposeNybbleblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, VecW* vecaligned_buf) {
  // Very similar to TransposeNypblock64() in pgenlib_internal.
  // vecaligned_buf must be vector-aligned and have size 8k
  const uint32_t buf_row_ct = DivUp(write_batch_size, 8);
  // fold the first 4 shuffles into the initial ingestion loop
  const uint32_t* initial_read_iter = R_CAST(const uint32_t*, read_iter);
  const uint32_t* initial_read_end = &(initial_read_iter[buf_row_ct]);
  uint32_t* initial_target_iter = R_CAST(uint32_t*, vecaligned_buf);
  const uint32_t read_u32_stride = read_ul_stride * (kBytesPerWord / 4);
  const uint32_t read_batch_rem = kNybblesPerCacheline - read_batch_size;
  for (; initial_read_iter != initial_read_end; ++initial_read_iter) {
    const uint32_t* read_iter_tmp = initial_read_iter;
    for (uint32_t ujj = 0; ujj != read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_u32_stride]);
    }
    if (!read_batch_rem) {
      continue;
    }
    memset(initial_target_iter, 0, read_batch_rem * 4);
    initial_target_iter = &(initial_target_iter[read_batch_rem]);
  }

  // 32 bit spacing -> 4
  const VecW* source_iter = vecaligned_buf;
  const VecW m4 = VCONST_W(kMask0F0F);
  const uint32_t buf_fullrow_ct = write_batch_size / 8;
  const uint32_t eightword_ct = DivUp(read_batch_size, 16);
  uintptr_t* target_iter0 = write_iter;
  uint32_t cur_dst_row_ct = 8;
#  ifdef USE_SSE42
  const VecW gather_u16s = vecw_setr8(0, 8, 1, 9, 2, 10, 3, 11,
                                      4, 12, 5, 13, 6, 14, 7, 15);
#  else
  const VecW m8 = VCONST_W(kMask00FF);
#  endif
#  ifdef USE_AVX2
  // movemask is slower even in AVX2 case
  const VecW gather_u32s = vecw_setr8(0, 1, 8, 9, 2, 3, 10, 11,
                                      4, 5, 12, 13, 6, 7, 14, 15);
  for (uint32_t buf_row_idx = 0; ; ++buf_row_idx) {
    if (buf_row_idx >= buf_fullrow_ct) {
      if (buf_row_idx == buf_row_ct) {
        return;
      }
      cur_dst_row_ct = write_batch_size % 8;
    }
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    uintptr_t* target_iter2 = &(target_iter1[write_ul_stride]);
    uintptr_t* target_iter3 = &(target_iter2[write_ul_stride]);
    uintptr_t* target_iter4 = &(target_iter3[write_ul_stride]);
    uintptr_t* target_iter5 = &(target_iter4[write_ul_stride]);
    uintptr_t* target_iter6 = &(target_iter5[write_ul_stride]);
    uintptr_t* target_iter7 = &(target_iter6[write_ul_stride]);
    for (uint32_t dvidx = 0; dvidx != eightword_ct; ++dvidx) {
      const VecW loader0 = source_iter[dvidx * 2];
      const VecW loader1 = source_iter[dvidx * 2 + 1];
      VecW even_nybbles0 = loader0 & m4;
      VecW odd_nybbles0 = vecw_and_notfirst(m4, loader0);
      VecW even_nybbles1 = loader1 & m4;
      VecW odd_nybbles1 = vecw_and_notfirst(m4, loader1);
      even_nybbles0 = even_nybbles0 | vecw_srli(even_nybbles0, 28);
      odd_nybbles0 = vecw_slli(odd_nybbles0, 28) | odd_nybbles0;
      even_nybbles1 = even_nybbles1 | vecw_srli(even_nybbles1, 28);
      odd_nybbles1 = vecw_slli(odd_nybbles1, 28) | odd_nybbles1;
      // Label the bytes in even_nybbles0 (0, 1, 2, ..., 31), and the bytes in
      // even_nybbles1 (32, 33, ..., 63).  We wish to generate the following
      // lane-and-vector-crossing permutation:
      //   (0, 8, 16, 24, 32, 40, 48, 56, 1, 9, 17, 25, 33, 41, 49, 57)
      //   (2, 10, 18, 26, 34, 42, 50, 58, 3, 11, 19, 27, 35, 43, 51, 59)

      // first shuffle:
      //   (0, 8, 1, 9, 2, 10, 3, 11, _, _, _, _, _, _, _, _,
      //    16, 24, 17, 25, 18, 26, 19, 27, _, _, _, _, _, _, _, _)
      //
      //   (32, 40, 33, 41, 34, 42, 35, 43, _, _, _, _, _, _, _, _,
      //    48, 56, 49, 57, 50, 58, 51, 59, _, _, _, _, _, _, _, _)
      //
      // _mm256_unpacklo_epi16:
      //   (0, 8, 32, 40, 1, 9, 33, 41, 2, 10, 34, 42, 3, 11, 35, 43,
      //    16, 24, 48, 56, 17, 25, 49, 57, 18, 26, 50, 58, 19, 27, 51, 59)
      //
      // {0, 2, 1, 3} permute:
      //   (0, 8, 32, 40, 1, 9, 33, 41, 16, 24, 48, 56, 17, 25, 49, 57,
      //    2, 10, 34, 42, 3, 11, 35, 43, 18, 26, 50, 58, 19, 27, 51, 59)
      //
      // final shuffle gives us what we want.
      even_nybbles0 = vecw_shuffle8(even_nybbles0, gather_u16s);
      odd_nybbles0 = vecw_shuffle8(odd_nybbles0, gather_u16s);
      even_nybbles1 = vecw_shuffle8(even_nybbles1, gather_u16s);
      odd_nybbles1 = vecw_shuffle8(odd_nybbles1, gather_u16s);

      VecW target_even = vecw_unpacklo16(even_nybbles0, even_nybbles1);
      VecW target_odd = vecw_unpackhi16(odd_nybbles0, odd_nybbles1);

      target_even = vecw_permute0xd8_if_avx2(target_even);
      target_odd = vecw_permute0xd8_if_avx2(target_odd);

      target_even = vecw_shuffle8(target_even, gather_u32s);
      target_odd = vecw_shuffle8(target_odd, gather_u32s);

      // tried using _mm_stream_si64 here, that totally sucked
      switch (cur_dst_row_ct) {
        case 8:
          target_iter7[dvidx] = vecw_extract64_3(target_odd);
          // fall through
        case 7:
          target_iter6[dvidx] = vecw_extract64_3(target_even);
          // fall through
        case 6:
          target_iter5[dvidx] = vecw_extract64_2(target_odd);
          // fall through
        case 5:
          target_iter4[dvidx] = vecw_extract64_2(target_even);
          // fall through
        case 4:
          target_iter3[dvidx] = vecw_extract64_1(target_odd);
          // fall through
        case 3:
          target_iter2[dvidx] = vecw_extract64_1(target_even);
          // fall through
        case 2:
          target_iter1[dvidx] = vecw_extract64_0(target_odd);
          // fall through
        default:
          target_iter0[dvidx] = vecw_extract64_0(target_even);
      }
    }
    source_iter = &(source_iter[(4 * kPglNybbleTransposeBatch) / kBytesPerVec]);
    target_iter0 = &(target_iter7[write_ul_stride]);
  }
#  else  // !USE_AVX2
  for (uint32_t buf_row_idx = 0; ; ++buf_row_idx) {
    if (buf_row_idx >= buf_fullrow_ct) {
      if (buf_row_idx == buf_row_ct) {
        return;
      }
      cur_dst_row_ct = write_batch_size % 8;
    }
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    uintptr_t* target_iter2 = &(target_iter1[write_ul_stride]);
    uintptr_t* target_iter3 = &(target_iter2[write_ul_stride]);
    uintptr_t* target_iter4 = &(target_iter3[write_ul_stride]);
    uintptr_t* target_iter5 = &(target_iter4[write_ul_stride]);
    uintptr_t* target_iter6 = &(target_iter5[write_ul_stride]);
    uintptr_t* target_iter7 = &(target_iter6[write_ul_stride]);
    for (uint32_t qvidx = 0; qvidx != eightword_ct; ++qvidx) {
      const VecW loader0 = source_iter[qvidx * 4];
      const VecW loader1 = source_iter[qvidx * 4 + 1];
      const VecW loader2 = source_iter[qvidx * 4 + 2];
      const VecW loader3 = source_iter[qvidx * 4 + 3];
      VecW even_nybbles0 = loader0 & m4;
      VecW odd_nybbles0 = vecw_and_notfirst(m4, loader0);
      VecW even_nybbles1 = loader1 & m4;
      VecW odd_nybbles1 = vecw_and_notfirst(m4, loader1);
      VecW even_nybbles2 = loader2 & m4;
      VecW odd_nybbles2 = vecw_and_notfirst(m4, loader2);
      VecW even_nybbles3 = loader3 & m4;
      VecW odd_nybbles3 = vecw_and_notfirst(m4, loader3);
      even_nybbles0 = even_nybbles0 | vecw_srli(even_nybbles0, 28);
      odd_nybbles0 = vecw_slli(odd_nybbles0, 28) | odd_nybbles0;
      even_nybbles1 = even_nybbles1 | vecw_srli(even_nybbles1, 28);
      odd_nybbles1 = vecw_slli(odd_nybbles1, 28) | odd_nybbles1;
      even_nybbles2 = even_nybbles2 | vecw_srli(even_nybbles2, 28);
      odd_nybbles2 = vecw_slli(odd_nybbles2, 28) | odd_nybbles2;
      even_nybbles3 = even_nybbles3 | vecw_srli(even_nybbles3, 28);
      odd_nybbles3 = vecw_slli(odd_nybbles3, 28) | odd_nybbles3;
      // Label the bytes in even_nybbles0 (0, 1, 2, ..., 15), the bytes in
      // even_nybbles1 (16, 17, ..., 31), ..., up to even_nybbles3 being (48,
      // 49, ..., 63).  We wish to generate the following vector-crossing
      // permutation:
      //   (0, 8, 16, 24, 32, 40, 48, 56, 1, 9, 17, 25, 33, 41, 49, 57)
      //   (2, 10, 18, 26, 34, 42, 50, 58, 3, 11, 19, 27, 35, 43, 51, 59)

      // first shuffle:
      //   (0, 8, 1, 9, 2, 10, 3, 11, _, _, _, _, _, _, _, _)
      //   (16, 24, 17, 25, 18, 26, 19, 27, _, _, _, _, _, _, _, _)
      //   (32, 40, 33, 41, 34, 42, 35, 43, _, _, _, _, _, _, _, _)
      //   (48, 56, 49, 57, 50, 58, 51, 59, _, _, _, _, _, _, _, _)

      // _mm_unpacklo_epi16:
      //   (0, 8, 16, 24, 1, 9, 17, 25, 2, 10, 18, 26, 3, 11, 19, 27)
      //   (32, 40, 48, 56, 33, 41, 49, 57, 34, 42, 50, 58, 35, 43, 51, 59)
      //
      // finish with _mm_unpack{lo,hi}_epi32
#    ifdef USE_SSE42
      even_nybbles0 = vecw_shuffle8(even_nybbles0, gather_u16s);
      odd_nybbles0 = vecw_shuffle8(odd_nybbles0, gather_u16s);
      even_nybbles1 = vecw_shuffle8(even_nybbles1, gather_u16s);
      odd_nybbles1 = vecw_shuffle8(odd_nybbles1, gather_u16s);
      even_nybbles2 = vecw_shuffle8(even_nybbles2, gather_u16s);
      odd_nybbles2 = vecw_shuffle8(odd_nybbles2, gather_u16s);
      even_nybbles3 = vecw_shuffle8(even_nybbles3, gather_u16s);
      odd_nybbles3 = vecw_shuffle8(odd_nybbles3, gather_u16s);
#    else
      VecW tmp_lo = vecw_unpacklo8(even_nybbles0, odd_nybbles0);
      VecW tmp_hi = vecw_unpackhi8(even_nybbles0, odd_nybbles0);
      even_nybbles0 = vecw_blendv(vecw_slli(tmp_hi, 8), tmp_lo, m8);
      odd_nybbles0 = vecw_blendv(tmp_hi, vecw_srli(tmp_lo, 8), m8);
      tmp_lo = vecw_unpacklo8(even_nybbles1, odd_nybbles1);
      tmp_hi = vecw_unpackhi8(even_nybbles1, odd_nybbles1);
      even_nybbles1 = vecw_blendv(vecw_slli(tmp_hi, 8), tmp_lo, m8);
      odd_nybbles1 = vecw_blendv(tmp_hi, vecw_srli(tmp_lo, 8), m8);
      tmp_lo = vecw_unpacklo8(even_nybbles2, odd_nybbles2);
      tmp_hi = vecw_unpackhi8(even_nybbles2, odd_nybbles2);
      even_nybbles2 = vecw_blendv(vecw_slli(tmp_hi, 8), tmp_lo, m8);
      odd_nybbles2 = vecw_blendv(tmp_hi, vecw_srli(tmp_lo, 8), m8);
      tmp_lo = vecw_unpacklo8(even_nybbles3, odd_nybbles3);
      tmp_hi = vecw_unpackhi8(even_nybbles3, odd_nybbles3);
      even_nybbles3 = vecw_blendv(vecw_slli(tmp_hi, 8), tmp_lo, m8);
      odd_nybbles3 = vecw_blendv(tmp_hi, vecw_srli(tmp_lo, 8), m8);
#    endif

      const VecW even_lo = vecw_unpacklo16(even_nybbles0, even_nybbles1);
      const VecW odd_lo = vecw_unpackhi16(odd_nybbles0, odd_nybbles1);
      const VecW even_hi = vecw_unpacklo16(even_nybbles2, even_nybbles3);
      const VecW odd_hi = vecw_unpackhi16(odd_nybbles2, odd_nybbles3);

      const VecW final02 = vecw_unpacklo32(even_lo, even_hi);
      const VecW final13 = vecw_unpacklo32(odd_lo, odd_hi);
      const VecW final46 = vecw_unpackhi32(even_lo, even_hi);
      const VecW final57 = vecw_unpackhi32(odd_lo, odd_hi);
      switch (cur_dst_row_ct) {
        case 8:
          target_iter7[qvidx] = vecw_extract64_1(final57);
          // fall through
        case 7:
          target_iter6[qvidx] = vecw_extract64_1(final46);
          // fall through
        case 6:
          target_iter5[qvidx] = vecw_extract64_0(final57);
          // fall through
        case 5:
          target_iter4[qvidx] = vecw_extract64_0(final46);
          // fall through
        case 4:
          target_iter3[qvidx] = vecw_extract64_1(final13);
          // fall through
        case 3:
          target_iter2[qvidx] = vecw_extract64_1(final02);
          // fall through
        case 2:
          target_iter1[qvidx] = vecw_extract64_0(final13);
          // fall through
        default:
          target_iter0[qvidx] = vecw_extract64_0(final02);
      }
    }
    source_iter = &(source_iter[(4 * kPglNybbleTransposeBatch) / kBytesPerVec]);
    target_iter0 = &(target_iter7[write_ul_stride]);
  }
#  endif  // !USE_AVX2
}
#else  // !__LP64__
static_assert(kWordsPerVec == 1, "TransposeNybbleblock() needs to be updated.");
void TransposeNybbleblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, VecW* vecaligned_buf) {
  // Very similar to TransposeNypblock32() in pgenlib_internal.
  // vecaligned_buf must be vector-aligned and have size 8k
  const uint32_t buf_row_ct = NybbleCtToByteCt(write_batch_size);
  // fold the first 6 shuffles into the initial ingestion loop
  const unsigned char* initial_read_iter = R_CAST(const unsigned char*, read_iter);
  const unsigned char* initial_read_end = &(initial_read_iter[buf_row_ct]);
  unsigned char* initial_target_iter = R_CAST(unsigned char*, vecaligned_buf);
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t read_batch_rem = kNybblesPerCacheline - read_batch_size;
  for (; initial_read_iter != initial_read_end; ++initial_read_iter) {
    const unsigned char* read_iter_tmp = initial_read_iter;
    for (uint32_t ujj = 0; ujj != read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_byte_stride]);
    }
    initial_target_iter = memsetua(initial_target_iter, 0, read_batch_rem);
  }

  // 8 bit spacing -> 4
  const VecW* source_iter = vecaligned_buf;
  uintptr_t* target_iter0 = write_iter;
  const uint32_t buf_fullrow_ct = write_batch_size / 2;
  const uint32_t write_word_ct = NybbleCtToWordCt(read_batch_size);
  for (uint32_t uii = 0; uii != buf_fullrow_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    for (uint32_t ujj = 0; ujj != write_word_ct; ++ujj) {
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
      target_iter0[ujj] = S_CAST(Halfword, target_word0_lo) | (target_word0_hi << kBitsPerWordD2);
      target_iter1[ujj] = S_CAST(Halfword, target_word1_lo) | (target_word1_hi << kBitsPerWordD2);
    }
    source_iter = &(source_iter[2 * (kWordsPerCacheline - write_word_ct)]);
    target_iter0 = &(target_iter1[write_ul_stride]);
  }
  const uint32_t remainder = write_batch_size % 2;
  if (!remainder) {
    return;
  }
  for (uint32_t ujj = 0; ujj != write_word_ct; ++ujj) {
    const uintptr_t source_word_lo = *source_iter++;
    const uintptr_t source_word_hi = *source_iter++;
    uintptr_t target_word0_lo = source_word_lo & kMask0F0F;
    uintptr_t target_word0_hi = source_word_hi & kMask0F0F;
    target_word0_lo = (target_word0_lo | (target_word0_lo >> 4)) & kMask00FF;
    target_word0_hi = (target_word0_hi | (target_word0_hi >> 4)) & kMask00FF;
    target_word0_lo = target_word0_lo | (target_word0_lo >> kBitsPerWordD4);
    target_word0_hi = target_word0_hi | (target_word0_hi >> kBitsPerWordD4);
    target_iter0[ujj] = S_CAST(Halfword, target_word0_lo) | (target_word0_hi << kBitsPerWordD2);
  }
}
#endif  // !__LP64__

#ifdef __LP64__
#  ifdef USE_AVX2
const unsigned char kLeadMask[2 * kBytesPerVec] __attribute__ ((aligned (64))) =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
   255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};
#  else
const unsigned char kLeadMask[2 * kBytesPerVec] __attribute__ ((aligned (32))) =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};
#  endif

uintptr_t BytesumArr(const void* bytearr, uintptr_t byte_ct) {
  uintptr_t tot = 0;
  if (byte_ct < kBytesPerVec) {
    const unsigned char* bytearr_uc = S_CAST(const unsigned char*, bytearr);
    for (uintptr_t ulii = 0; ulii != byte_ct; ++ulii) {
      tot += bytearr_uc[ulii];
    }
    return tot;
  }
  const unsigned char* bytearr_uc_iter = S_CAST(const unsigned char*, bytearr);
  const unsigned char* bytearr_uc_final = &(bytearr_uc_iter[byte_ct - kBytesPerVec]);
  const VecW m0 = vecw_setzero();
  VecW acc = vecw_setzero();
  while (bytearr_uc_iter < bytearr_uc_final) {
    const VecW cur_vec = vecw_loadu(bytearr_uc_iter);
    acc = acc + vecw_sad(cur_vec, m0);
    bytearr_uc_iter = &(bytearr_uc_iter[kBytesPerVec]);
  }
  VecW cur_vec = vecw_loadu(bytearr_uc_final);
  const uintptr_t overlap_byte_ct = bytearr_uc_iter - bytearr_uc_final;
  const VecW mask_vec = vecw_loadu(&(kLeadMask[kBytesPerVec - overlap_byte_ct]));
  cur_vec = cur_vec & mask_vec;
  acc = acc + vecw_sad(cur_vec, m0);
  return HsumW(acc);
}

#else  // !__LP64__
uintptr_t BytesumArr(const void* bytearr, uintptr_t byte_ct) {
  // Assumes sum < 2^32.
#  ifdef NO_UNALIGNED
#    error "Unaligned accesses in BytesumArr()."
#  endif
  const uint32_t word_ct = byte_ct / kBytesPerWord;
  const uintptr_t* bytearr_alias_iter = S_CAST(const uintptr_t*, bytearr);
  const uint32_t wordblock_idx_trail = word_ct / 256;
  const uint32_t wordblock_idx_end = DivUp(word_ct, 256);
  uint32_t wordblock_len = 256;
  uintptr_t tot = 0;
  for (uint32_t wordblock_idx = 0; ; ++wordblock_idx) {
    if (wordblock_idx >= wordblock_idx_trail) {
      if (wordblock_idx == wordblock_idx_end) {
        byte_ct = byte_ct % kBytesPerWord;
        const unsigned char* bytearr_alias_iter2 = R_CAST(const unsigned char*, bytearr_alias_iter);
        for (uint32_t uii = 0; uii != byte_ct; ++uii) {
          tot += bytearr_alias_iter2[uii];
        }
        return tot;
      }
      wordblock_len = word_ct % 256;
    }
    const uintptr_t* bytearr_alias_stop = &(bytearr_alias_iter[wordblock_len]);
    uintptr_t acc_even = 0;
    uintptr_t acc_odd = 0;
    do {
      uintptr_t cur_word = *bytearr_alias_iter++;
      acc_even += cur_word & kMask00FF;
      acc_odd += (cur_word >> 8) & kMask00FF;
    } while (bytearr_alias_iter < bytearr_alias_stop);
    acc_even = S_CAST(Halfword, acc_even) + (acc_even >> kBitsPerWordD2);
    acc_odd = S_CAST(Halfword, acc_odd) + (acc_odd >> kBitsPerWordD2);
    tot += acc_even + acc_odd;
  }
}
#endif  // !__LP64__

uintptr_t CountByte(const void* bytearr, unsigned char ucc, uintptr_t byte_ct) {
#ifdef __LP64__
  if (byte_ct < kBytesPerVec) {
#endif
    const unsigned char* bytearr_uc = S_CAST(const unsigned char*, bytearr);
    uintptr_t tot = 0;
    for (uintptr_t ulii = 0; ulii != byte_ct; ++ulii) {
      tot += (bytearr_uc[ulii] == ucc);
    }
    return tot;
#ifdef __LP64__
  }
  const unsigned char* bytearr_uc_iter = S_CAST(const unsigned char*, bytearr);
  const VecW m0 = vecw_setzero();
  const VecUc match_vvec = vecuc_set1(ucc);
  VecW acc = vecw_setzero();
  while (byte_ct > 255 * kBytesPerVec) {
    VecUc inner_acc = vecuc_setzero();
    for (uint32_t uii = 0; uii != 255; ++uii) {
      const VecUc cur_vvec = vecuc_loadu(bytearr_uc_iter);
      bytearr_uc_iter = &(bytearr_uc_iter[kBytesPerVec]);
      inner_acc = inner_acc - (cur_vvec == match_vvec);
    }
    acc = acc + vecw_sad(R_CAST(VecW, inner_acc), m0);
    byte_ct -= 255 * kBytesPerVec;
  }
  const unsigned char* bytearr_uc_final = &(bytearr_uc_iter[byte_ct - kBytesPerVec]);
  VecUc inner_acc = vecuc_setzero();
  while (bytearr_uc_iter < bytearr_uc_final) {
    const VecUc cur_vvec = vecuc_loadu(bytearr_uc_iter);
    bytearr_uc_iter = &(bytearr_uc_iter[kBytesPerVec]);
    inner_acc = inner_acc - (cur_vvec == match_vvec);
  }
  VecUc cur_vvec = vecuc_loadu(bytearr_uc_final);
  const uintptr_t overlap_byte_ct = bytearr_uc_iter - bytearr_uc_final;
  const VecUc mask_vvec = vecuc_loadu(&(kLeadMask[kBytesPerVec - overlap_byte_ct]));
  cur_vvec = (cur_vvec == match_vvec) & mask_vvec;
  inner_acc = inner_acc - cur_vvec;
  acc = acc + vecw_sad(R_CAST(VecW, inner_acc), m0);
  return HsumW(acc);
#endif  // __LP64__
}

uintptr_t CountU16(const void* u16arr, uint16_t usii, uintptr_t u16_ct) {
#ifdef __LP64__
  if (u16_ct < (kBytesPerVec / 2)) {
#endif
    const uint16_t* u16arr_alias = S_CAST(const uint16_t*, u16arr);
    uintptr_t tot = 0;
    for (uintptr_t ulii = 0; ulii != u16_ct; ++ulii) {
      tot += (u16arr_alias[ulii] == usii);
    }
    return tot;
#ifdef __LP64__
  }
  const uint16_t* u16arr_iter = S_CAST(const uint16_t*, u16arr);
  const VecW m0 = vecw_setzero();
  const VecU16 match_vvec = vecu16_set1(usii);
  VecW acc = vecw_setzero();
  // can also use larger loop and a slightly different accumulation algorithm,
  // but it should make practically no difference; lets keep these loops as
  // similar as possible for now.
  while (u16_ct > 255 * (kBytesPerVec / 2)) {
    VecU16 inner_acc = vecu16_setzero();
    for (uint32_t uii = 0; uii != 255; ++uii) {
      const VecU16 cur_vvec = vecu16_loadu(u16arr_iter);
      u16arr_iter = &(u16arr_iter[kBytesPerVec / 2]);
      inner_acc = inner_acc - (cur_vvec == match_vvec);
    }
    acc = acc + vecw_sad(R_CAST(VecW, inner_acc), m0);
    u16_ct -= 255 * (kBytesPerVec / 2);
  }
  const uint16_t* u16arr_final = &(u16arr_iter[u16_ct - (kBytesPerVec / 2)]);
  VecU16 inner_acc = vecu16_setzero();
  while (u16arr_iter < u16arr_final) {
    const VecU16 cur_vvec = vecu16_loadu(u16arr_iter);
    u16arr_iter = &(u16arr_iter[kBytesPerVec / 2]);
    inner_acc = inner_acc - (cur_vvec == match_vvec);
  }
  VecU16 cur_vvec = vecu16_loadu(u16arr_final);
  const uintptr_t overlap_u16_ct = u16arr_iter - u16arr_final;
  const VecU16 mask_vvec = vecu16_loadu(&(kLeadMask[kBytesPerVec - 2 * overlap_u16_ct]));
  cur_vvec = (cur_vvec == match_vvec) & mask_vvec;
  inner_acc = inner_acc - cur_vvec;
  acc = acc + vecw_sad(R_CAST(VecW, inner_acc), m0);
  return HsumW(acc);
#endif  // __LP64__
}

uint32_t Copy1bit8Subset(const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uintptr_t* __restrict sample_include, uint32_t src_subset_size, uint32_t sample_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals) {
  if (!src_subset_size) {
    return 0;
  }
  CopyBitarrSubset(src_subset, sample_include, sample_ct, dst_subset);
  const unsigned char* src_vals_uc = S_CAST(const unsigned char*, src_vals);
  unsigned char* dst_vals_uc = S_CAST(unsigned char*, dst_vals);
  unsigned char* dst_vals_iter = dst_vals_uc;
  uintptr_t sample_widx = 0;
  uintptr_t src_subset_bits = src_subset[0];
  for (uint32_t src_idx = 0; src_idx != src_subset_size; ++src_idx) {
    const uintptr_t lowbit = BitIter1y(src_subset, &sample_widx, &src_subset_bits);
    if (sample_include[sample_widx] & lowbit) {
      *dst_vals_iter++ = src_vals_uc[src_idx];
    }
  }
  return dst_vals_iter - dst_vals_uc;
}

uint32_t Copy1bit16Subset(const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uintptr_t* __restrict sample_include, uint32_t src_subset_size, uint32_t sample_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals) {
  if (!src_subset_size) {
    return 0;
  }
  CopyBitarrSubset(src_subset, sample_include, sample_ct, dst_subset);
  const uint16_t* src_vals_u16 = S_CAST(const uint16_t*, src_vals);
  uint16_t* dst_vals_u16 = S_CAST(uint16_t*, dst_vals);
  uint16_t* dst_vals_iter = dst_vals_u16;
  uintptr_t sample_widx = 0;
  uintptr_t src_subset_bits = src_subset[0];
  for (uint32_t src_idx = 0; src_idx != src_subset_size; ++src_idx) {
    const uintptr_t lowbit = BitIter1y(src_subset, &sample_widx, &src_subset_bits);
    if (sample_include[sample_widx] & lowbit) {
      *dst_vals_iter++ = src_vals_u16[src_idx];
    }
  }
  return dst_vals_iter - dst_vals_u16;
}

// 'Unsafe' because it assumes high bits of every byte are 0.
void Reduce8to4bitInplaceUnsafe(uintptr_t entry_ct, uintptr_t* arr) {
#ifdef __LP64__
  const uintptr_t fullvec_ct = entry_ct / (kBytesPerVec * 2);
  const VecW m8 = VCONST_W(kMask00FF);
  VecW* varr = R_CAST(VecW*, arr);
  for (uintptr_t write_vidx = 0; write_vidx != fullvec_ct; ++write_vidx) {
    VecW v0 = varr[write_vidx * 2];
    VecW v1 = varr[write_vidx * 2 + 1];
    v0 = v0 | vecw_srli(v0, 4);
    v1 = v1 | vecw_srli(v1, 4);
    varr[write_vidx] = vecw_gather_even(v0, v1, m8);
  }
  uintptr_t write_idx = fullvec_ct * kWordsPerVec;
  if (write_idx == entry_ct * 2) {
    return;
  }
#else
  uintptr_t write_idx = 0;
#endif
  // Read two words at a time and write one.
  // We could instead read one word and write a Halfword at a time, but I'd
  // rather not worry about the strict-aliasing issues involved there.
  const uintptr_t write_idx_last = (entry_ct - 1) / (kBytesPerWord * 2);
  uintptr_t write_word;
  for (; ; ++write_idx) {
    uintptr_t inword0 = arr[2 * write_idx];
    uintptr_t inword1 = arr[2 * write_idx + 1];
#ifdef USE_AVX2
    inword0 = _pext_u64(inword0, kMask0F0F);
    inword1 = _pext_u64(inword1, kMask0F0F);
#else
    // 0 . 1 . 2 . 3 . 4 . 5 . 6 . 7 .
    // (or just 0 . 1 . 2 . 3 . in 32-bit case)
    inword0 = (inword0 | (inword0 >> 4)) & kMask00FF;
    inword1 = (inword1 | (inword1 >> 4)) & kMask00FF;
    // 0 1 . . 2 3 . . 4 5 . . 6 7 . .
#  ifdef __LP64__
    inword0 = (inword0 | (inword0 >> 8)) & kMask0000FFFF;
    inword1 = (inword1 | (inword1 >> 8)) & kMask0000FFFF;
    // 0 1 2 3 . . . . 4 5 6 7 . . . .
#  endif
    inword0 = S_CAST(Halfword, inword0 | (inword0 >> kBitsPerWordD4));
    inword1 = S_CAST(Halfword, inword1 | (inword1 >> kBitsPerWordD4));
#endif
    write_word = inword0 | (inword1 << kBitsPerWordD2);
    if (write_idx == write_idx_last) {
      break;
    }
    arr[write_idx] = write_word;
  }
  const uint32_t remaining_entry_ct = ModNz(entry_ct, kBytesPerWord * 2);
  arr[write_idx] = bzhi_max(write_word, remaining_entry_ct * 4);
}

#ifdef __cplusplus
}  // namespace plink2
#endif
