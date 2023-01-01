#ifndef __PLINK2_BITS_H__
#define __PLINK2_BITS_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2023 Shaun Purcell,
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


// Bitarray support.  (Inline single-word operations are in plink2_base.h.)

#include "plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

#if defined(__LP64__) && !defined(USE_AVX2)
// may also want a version which doesn't always apply kMask5555
void Pack32bTo16bMask(const void* words, uintptr_t ct_32b, void* dest);

HEADER_INLINE void PackWordsToHalfwordsMask(const uintptr_t* words, uintptr_t word_ct, Halfword* dest) {
  uintptr_t widx = 0;
  if (word_ct >= (32 / kBytesPerWord)) {
    const uintptr_t ct_32b = word_ct / (32 / kBytesPerWord);
    Pack32bTo16bMask(words, ct_32b, dest);
    widx = ct_32b * (32 / kBytesPerWord);
  }
  for (; widx != word_ct; ++widx) {
    dest[widx] = PackWordToHalfwordMask5555(words[widx]);
  }
}
#else
HEADER_INLINE void PackWordsToHalfwordsMask(const uintptr_t* words, uintptr_t word_ct, Halfword* dest) {
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    dest[widx] = PackWordToHalfwordMask5555(words[widx]);
  }
}
#endif

// ok for ct == 0
void SetAllBits(uintptr_t ct, uintptr_t* bitarr);

// "Nz" added to names to make it obvious these require positive len
void FillBitsNz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr);
void ClearBitsNz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr);

void BitvecAnd(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

void BitvecInvmask(const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

void BitvecOr(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec);

void BitvecInvert(uintptr_t word_ct, uintptr_t* main_bitvec);

void BitvecXorCopy(const uintptr_t* __restrict source1_bitvec, const uintptr_t* __restrict source2_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec);

void BitvecInvertCopy(const uintptr_t* __restrict source_bitvec, uintptr_t word_ct, uintptr_t* __restrict target_bitvec);

// These ensure the trailing bits are zeroed out.
// 'AlignedBitarr' instead of Bitvec since this takes bit_ct instead of word_ct
// as the size argument, and zeroes trailing bits.
HEADER_INLINE void AlignedBitarrInvert(uintptr_t bit_ct, uintptr_t* main_bitvec) {
  const uintptr_t fullword_ct = bit_ct / kBitsPerWord;
  BitvecInvert(fullword_ct, main_bitvec);
  const uint32_t trail_ct = bit_ct % kBitsPerWord;
  if (trail_ct) {
    main_bitvec[fullword_ct] = bzhi(~main_bitvec[fullword_ct], trail_ct);
  }
}

HEADER_INLINE void AlignedBitarrInvertCopy(const uintptr_t* __restrict source_bitvec, uintptr_t bit_ct, uintptr_t* __restrict target_bitvec) {
  const uintptr_t fullword_ct = bit_ct / kBitsPerWord;
  BitvecInvertCopy(source_bitvec, fullword_ct, target_bitvec);
  const uint32_t trail_ct = bit_ct % kBitsPerWord;
  if (trail_ct) {
    target_bitvec[fullword_ct] = bzhi(~source_bitvec[fullword_ct], trail_ct);
  }
}

// Functions with "adv" in the name generally take an index or char-pointer as
// an argument and return its new value, while "mov" functions take a
// pointer-to-index or pointer-to-char-pointer and move it.

// These return the current index if the corresponding bit satisfies the
// condition.
uintptr_t AdvTo1Bit(const uintptr_t* bitarr, uintptr_t loc);

uintptr_t AdvTo0Bit(const uintptr_t* bitarr, uintptr_t loc);

// uintptr_t NextNonmissingUnsafe(const uintptr_t* genoarr, uintptr_t loc);

uint32_t AdvBoundedTo1Bit(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil);

uintptr_t AdvBoundedTo0Bit(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil);

uint32_t FindLast1BitBefore(const uintptr_t* bitarr, uint32_t loc);

// possible todo: check if movemask-based solution is better in AVX2 case
HEADER_INLINE uint32_t AllWordsAreZero(const uintptr_t* word_arr, uintptr_t word_ct) {
  while (word_ct--) {
    if (*word_arr++) {
      return 0;
    }
  }
  return 1;
}

HEADER_INLINE uint32_t AllBitsAreOne(const uintptr_t* bitarr, uintptr_t bit_ct) {
  const uintptr_t fullword_ct = bit_ct / kBitsPerWord;
  for (uintptr_t widx = 0; widx != fullword_ct; ++widx) {
    if (~(bitarr[widx])) {
      return 0;
    }
  }
  const uint32_t trailing_bit_ct = bit_ct % kBitsPerWord;
  return (!trailing_bit_ct) || ((~(bitarr[fullword_ct])) << (kBitsPerWord - trailing_bit_ct));
}

uint32_t AllBytesAreX(const unsigned char* bytes, unsigned char match, uintptr_t byte_ct);

// Updated PopcountWords() code is based on
// https://github.com/kimwalisch/libpopcnt .  libpopcnt license text follows.

/*
 * libpopcnt.h - C/C++ library for counting the number of 1 bits (bit
 * population count) in an array as quickly as possible using
 * specialized CPU instructions i.e. POPCNT, AVX2, AVX512, NEON.
 *
 * Copyright (c) 2016 - 2017, Kim Walisch
 * Copyright (c) 2016 - 2017, Wojciech Mula
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifdef USE_AVX2
// 'Csa' = carry, save, add
// If bb, cc, and *lp are bitvectors, this returns the carry bitvector and sets
// *lp to contain the low-order bits of the sums.  I.e. for each position:
//   if none of bb, cc, and *lp are set, *lp bit is zero and carry bit is zero
//   if exactly 1 is set, *lp bit becomes one and carry bit is zero
//   if exactly 2 are set, *lp bit becomes zero and carry bit is one
//   if all 3 are set, *lp bit becomes one and carry bit is one
HEADER_INLINE VecW Csa256(VecW bb, VecW cc, VecW* lp) {
  const VecW aa = *lp;
  const VecW uu = aa ^ bb;
  *lp = uu ^ cc;
  return (aa & bb) | (uu & cc);
}

HEADER_INLINE VecW CsaOne256(VecW bb, VecW* lp) {
  const VecW aa = *lp;
  *lp = aa ^ bb;
  return aa & bb;
}

HEADER_INLINE VecW PopcountVecAvx2(VecW vv) {
  const VecW lookup1 = vecw_setr8(4, 5, 5, 6, 5, 6, 6, 7,
                                  5, 6, 6, 7, 6, 7, 7, 8);
  const VecW lookup2 = vecw_setr8(4, 3, 3, 2, 3, 2, 2, 1,
                                  3, 2, 2, 1, 2, 1, 1, 0);

  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW lo = vv & m4;
  const VecW hi = vecw_srli(vv, 4) & m4;
  const VecW popcnt1 = vecw_shuffle8(lookup1, lo);
  const VecW popcnt2 = vecw_shuffle8(lookup2, hi);
  return vecw_sad(popcnt1, popcnt2);
}

HEADER_INLINE uintptr_t HsumW(VecW vv) {
  UniVec vu;
  vu.vw = vv;
  return vu.w[0] + vu.w[1] + vu.w[2] + vu.w[3];
  // _mm256_extract_epi64() only worth it if we don't need to extract all the
  // values.
  // (also, I wouldn't be surprised if the compiler recognized the pattern
  // above)
}

// This no longer has any restrictions on vec_ct, though it isn't worth the
// overhead for vec_ct < 16.
uintptr_t PopcountVecsAvx2(const VecW* bit_vvec, uintptr_t vec_ct);

HEADER_INLINE uintptr_t PopcountWords(const uintptr_t* bitvec, uintptr_t word_ct) {
  // Efficiently popcounts bitvec[0..(word_ct - 1)].  In the 64-bit case,
  // bitvec[] must be 16-byte aligned.
  // The PopcountWordsNzbase() wrapper takes care of starting from a later
  // index.
  uintptr_t tot = 0;
  if (word_ct >= 76) {
    assert(VecIsAligned(bitvec));
    const uintptr_t remainder = word_ct % kWordsPerVec;
    const uintptr_t main_block_word_ct = word_ct - remainder;
    tot = PopcountVecsAvx2(R_CAST(const VecW*, bitvec), main_block_word_ct / kWordsPerVec);
    bitvec = &(bitvec[main_block_word_ct]);
    word_ct = remainder;
  }
  // note that recent clang versions automatically expand this to a
  // full-service routine; takes ~50% longer than PopcountVecsAvx2 on >1kb
  // arrays, but way better than the naive loop
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    tot += PopcountWord(bitvec[widx]);
  }
  return tot;
}
#else  // !USE_AVX2
#  ifdef __LP64__
HEADER_INLINE uintptr_t HsumW(VecW vv) {
  UniVec vu;
  vu.vw = vv;
  return vu.w[0] + vu.w[1];
}
#  else
HEADER_INLINE uintptr_t HsumW(VecW vv) {
  return vv;
}
#  endif

// assumes vec_ct is a multiple of 3
uintptr_t PopcountVecsNoAvx2(const VecW* bit_vvec, uintptr_t vec_ct);

HEADER_INLINE uintptr_t PopcountWords(const uintptr_t* bitvec, uintptr_t word_ct) {
  uintptr_t tot = 0;
#ifndef USE_SSE42
  if (word_ct >= (3 * kWordsPerVec)) {
    // This has an asymptotic ~10% advantage in the SSE4.2 case, but word_ct
    // needs to be in the hundreds before the initial comparison even starts to
    // pay for itself.
    assert(VecIsAligned(bitvec));
    const uintptr_t remainder = word_ct % (3 * kWordsPerVec);
    const uintptr_t main_block_word_ct = word_ct - remainder;
    tot = PopcountVecsNoAvx2(R_CAST(const VecW*, bitvec), main_block_word_ct / kWordsPerVec);
    word_ct = remainder;
    bitvec = &(bitvec[main_block_word_ct]);
  }
#endif
  for (uintptr_t trailing_word_idx = 0; trailing_word_idx != word_ct; ++trailing_word_idx) {
    tot += PopcountWord(bitvec[trailing_word_idx]);
  }
  return tot;
}
#endif  // !USE_AVX2

uintptr_t PopcountWordsIntersect(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2_iter, uintptr_t word_ct);

uintptr_t PopcountWordsXor(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2_iter, uintptr_t word_ct);

// requires positive word_ct
// stay agnostic a bit longer re: word_ct := DIV_UP(entry_ct, kBitsPerWord)
// vs. word_ct := 1 + (entry_ct / kBitsPerWord)
// (this is a source of bugs, though; interface should probably be changed to
// use entry_ct once multiallelic/dosage implementation is done)
void FillCumulativePopcounts(const uintptr_t* subset_mask, uint32_t word_ct, uint32_t* cumulative_popcounts);

// If idx_list is a list of valid unfiltered indexes, this converts them
// in-place to corresponding filtered indexes.
void UidxsToIdxs(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, const uintptr_t idx_list_len, uint32_t* idx_list);

// These functions do not overread, but may write extra bytes up to the word
// boundary.
void Expand1bitTo8(const void* __restrict bytearr, uint32_t input_bit_ct, uint32_t incr, uintptr_t* __restrict dst);

void Expand1bitTo16(const void* __restrict bytearr, uint32_t input_bit_ct, uint32_t incr, uintptr_t* __restrict dst);


// might rename this to IsSet01 (guaranteeing 0/1 return value), and change
// IsSet() to bitarr[idx / kBitsPerWord] & (k1LU << (idx % kBitsPerWord)) since
// I'd expect that to play better with out-of-order execution.  but need to
// benchmark first.
HEADER_INLINE uintptr_t IsSet(const uintptr_t* bitarr, uintptr_t idx) {
  return (bitarr[idx / kBitsPerWord] >> (idx % kBitsPerWord)) & 1;
}

HEADER_INLINE void SetBit(uintptr_t idx, uintptr_t* bitarr) {
  bitarr[idx / kBitsPerWord] |= k1LU << (idx % kBitsPerWord);
}

HEADER_INLINE void ClearBit(uintptr_t idx, uintptr_t* bitarr) {
  bitarr[idx / kBitsPerWord] &= ~(k1LU << (idx % kBitsPerWord));
}

HEADER_INLINE void AssignBit(uintptr_t idx, uintptr_t newbit, uintptr_t* bitarr) {
  const uintptr_t inv_mask = k1LU << (idx % kBitsPerWord);
  uintptr_t* cur_word_ptr = &(bitarr[idx / kBitsPerWord]);
  *cur_word_ptr = ((*cur_word_ptr) & (~inv_mask)) | (inv_mask * newbit);
}

/*
HEADER_INLINE uintptr_t BitInnerIter1(uintptr_t uidx_base, uintptr_t* cur_bitsp, uintptr_t* cur_uidx_stopp) {
  const uintptr_t cur_bits = *cur_bitsp;
  const uint32_t uidx_start_lowbits = ctzw(*cur_bitsp);
  // Key idea is to iterate over sub-blocks of set bits in a single word, in
  // essentially the same manner as non-AVX2 CopyBitarrSubset() was doing.
  // This particular expression 'finds' the end of the current sub-block.
  const uintptr_t cur_bits_lfill_p1 = (cur_bits | (cur_bits - 1)) + 1;
  *cur_bitsp = cur_bits & cur_bits_lfill_p1;
  uint32_t uidx_stop_lowbits = kBitsPerWord;
  if (cur_bits_lfill_p1) {
    uidx_stop_lowbits = ctzw(cur_bits_lfill_p1);
  }
  *cur_uidx_stopp = uidx_base + uidx_stop_lowbits;
  return uidx_base + uidx_start_lowbits;
}
*/

HEADER_INLINE uintptr_t BitIter1(const uintptr_t* __restrict bitarr, uintptr_t* __restrict uidx_basep, uintptr_t* __restrict cur_bitsp) {
  uintptr_t cur_bits = *cur_bitsp;
  if (!cur_bits) {
    uintptr_t widx = (*uidx_basep) / kBitsPerWord;
    do {
      cur_bits = bitarr[++widx];
    } while (!cur_bits);
    *uidx_basep = widx * kBitsPerWord;
  }
  *cur_bitsp = cur_bits & (cur_bits - 1);
  return (*uidx_basep) + ctzw(cur_bits);
}

// Returns lowbit index instead of the full index.
HEADER_INLINE uint32_t BitIter1x(const uintptr_t* __restrict bitarr, uintptr_t* __restrict widxp, uintptr_t* __restrict cur_bitsp) {
  uintptr_t cur_bits = *cur_bitsp;
  while (!cur_bits) {
    cur_bits = bitarr[++(*widxp)];
  }
  *cur_bitsp = cur_bits & (cur_bits - 1);
  return ctzw(cur_bits);
}

// Returns isolated lowbit.
HEADER_INLINE uintptr_t BitIter1y(const uintptr_t* __restrict bitarr, uintptr_t* __restrict widxp, uintptr_t* __restrict cur_bitsp) {
  uintptr_t cur_bits = *cur_bitsp;
  while (!cur_bits) {
    cur_bits = bitarr[++(*widxp)];
  }
  const uintptr_t shifted_bit = cur_bits & (-cur_bits);
  *cur_bitsp = cur_bits ^ shifted_bit;
  return shifted_bit;
}

HEADER_INLINE void BitIter1Start(const uintptr_t* __restrict bitarr, uintptr_t restart_uidx, uintptr_t* __restrict uidx_basep, uintptr_t* __restrict cur_bitsp) {
  const uintptr_t widx = restart_uidx / kBitsPerWord;
  *cur_bitsp = bitarr[widx] & ((~k0LU) << (restart_uidx % kBitsPerWord));
  *uidx_basep = widx * kBitsPerWord;
}

HEADER_INLINE uintptr_t BitIter1NoAdv(const uintptr_t* __restrict bitarr, uintptr_t* __restrict uidx_basep, uintptr_t* __restrict cur_bitsp) {
  uintptr_t cur_bits = *cur_bitsp;
  if (!cur_bits) {
    uintptr_t widx = (*uidx_basep) / kBitsPerWord;
    do {
      cur_bits = bitarr[++widx];
    } while (!cur_bits);
    *uidx_basep = widx * kBitsPerWord;
    *cur_bitsp = cur_bits;
  }
  return (*uidx_basep) + ctzw(cur_bits);
}

HEADER_INLINE uintptr_t BitIter0(const uintptr_t* __restrict bitarr, uintptr_t* __restrict uidx_basep, uintptr_t* __restrict cur_inv_bitsp) {
  uintptr_t cur_inv_bits = *cur_inv_bitsp;
  if (!cur_inv_bits) {
    uintptr_t widx = (*uidx_basep) / kBitsPerWord;
    do {
      cur_inv_bits = ~bitarr[++widx];
    } while (!cur_inv_bits);
    *uidx_basep = widx * kBitsPerWord;
  }
  *cur_inv_bitsp = cur_inv_bits & (cur_inv_bits - 1);
  return (*uidx_basep) + ctzw(cur_inv_bits);
}

HEADER_INLINE void BitIter0Start(const uintptr_t* __restrict bitarr, uintptr_t restart_uidx, uintptr_t* __restrict uidx_basep, uintptr_t* __restrict cur_inv_bitsp) {
  const uintptr_t widx = restart_uidx / kBitsPerWord;
  *cur_inv_bitsp = (~bitarr[widx]) & ((~k0LU) << (restart_uidx % kBitsPerWord));
  *uidx_basep = widx * kBitsPerWord;
}

HEADER_INLINE uintptr_t BitIter0NoAdv(const uintptr_t* __restrict bitarr, uintptr_t* __restrict uidx_basep, uintptr_t* __restrict cur_inv_bitsp) {
  uintptr_t cur_inv_bits = *cur_inv_bitsp;
  if (!cur_inv_bits) {
    uintptr_t widx = (*uidx_basep) / kBitsPerWord;
    do {
      cur_inv_bits = ~bitarr[++widx];
    } while (!cur_inv_bits);
    *uidx_basep = widx * kBitsPerWord;
    *cur_inv_bitsp = cur_inv_bits;
  }
  return (*uidx_basep) + ctzw(cur_inv_bits);
}

// todo: test this against extracting a nonmissing bitarr first
/*
HEADER_INLINE void NextNonmissingUnsafeCk32(const uintptr_t* __restrict genoarr, uint32_t* __restrict loc_ptr) {
  if (GetNyparrEntry(genoarr, *loc_ptr) == 3) {
    *loc_ptr = NextNonmissingUnsafe(genoarr, *loc_ptr);
  }
}
*/

// Equivalent to popcount_bit_idx(subset_mask, 0, raw_idx).
HEADER_INLINE uint32_t RawToSubsettedPos(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, uint32_t raw_idx) {
  // this should be much better than keeping a uidx_to_idx array!
  // (update: there are more compact indexes, but postpone for now, this is
  // is nice and simple and gets us most of what we need.)
  const uint32_t raw_widx = raw_idx / kBitsPerWord;
  return subset_cumulative_popcounts[raw_widx] + PopcountWord(bzhi(subset_mask[raw_widx], raw_idx % kBitsPerWord));
}

HEADER_INLINE void ZeroTrailingBits(uintptr_t bit_ct, uintptr_t* bitarr) {
  const uint32_t trail_ct = bit_ct % kBitsPerWord;
  if (trail_ct) {
    bitarr[bit_ct / kBitsPerWord] = bzhi(bitarr[bit_ct / kBitsPerWord], trail_ct);
  }
}

#ifdef __LP64__
HEADER_INLINE void ZeroTrailingWords(uint32_t word_ct, uintptr_t* bitvec) {
  const uint32_t remainder = word_ct % kWordsPerVec;
  if (remainder) {
    ZeroWArr(kWordsPerVec - remainder, &(bitvec[word_ct]));
  }
}
#else
HEADER_INLINE void ZeroTrailingWords(__maybe_unused uint32_t word_ct, __maybe_unused uintptr_t* bitvec) {
}
#endif

HEADER_INLINE void CopyBitarr(const uintptr_t* __restrict src, uintptr_t bit_ct, uintptr_t* __restrict dst) {
  memcpy(dst, src, BitCtToWordCt(bit_ct) * kBytesPerWord);
}

// output_bit_idx_end is practically always subset_size
void CopyBitarrSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t output_bit_idx_end, uintptr_t* __restrict output_bitarr);

// expand_size + read_start_bit must be positive.
void ExpandBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, uint32_t word_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target);

// equivalent to calling ExpandBytearr() followed by CopyBitarrSubset()
void ExpandThenSubsetBytearr(const void* __restrict compact_bitarr, const uintptr_t* __restrict expand_mask, const uintptr_t* __restrict subset_mask, uint32_t expand_size, uint32_t subset_size, uint32_t read_start_bit, uintptr_t* __restrict target);

// mid_popcount must be positive
void ExpandBytearrNested(const void* __restrict compact_bitarr, const uintptr_t* __restrict mid_bitarr, const uintptr_t* __restrict top_expand_mask, uint32_t word_ct, uint32_t mid_popcount, uint32_t mid_start_bit, uintptr_t* __restrict mid_target, uintptr_t* __restrict compact_target);

// mid_popcount must be positive
// if mid_start_bit == 1, mid_popcount should not include that bit
void ExpandThenSubsetBytearrNested(const void* __restrict compact_bitarr, const uintptr_t* __restrict mid_bitarr, const uintptr_t* __restrict top_expand_mask, const uintptr_t* __restrict subset_mask, uint32_t subset_size, uint32_t mid_popcount, uint32_t mid_start_bit, uintptr_t* __restrict mid_target, uintptr_t* __restrict compact_target);

// these don't read past the end of bitarr
uintptr_t PopcountBytes(const void* bitarr, uintptr_t byte_ct);
uintptr_t PopcountBytesMasked(const void* bitarr, const uintptr_t* mask_arr, uintptr_t byte_ct);


// TransposeNypblock(), which is more plink-specific, is in pgenlib_misc
CONSTI32(kPglBitTransposeBatch, kBitsPerCacheline);
CONSTI32(kPglBitTransposeWords, kWordsPerCacheline);
// * Up to 512x512; vecaligned_buf must have size 64k
// * write_iter must be allocated up to at least
//   RoundUpPow2(write_batch_size, 2) rows
// * We use pointers with different types to read from and write to buf0/buf1,
//   so defining the base type as unsigned char* is theoretically necessary to
//   avoid breaking strict-aliasing rules, while the restrict qualifiers should
//   tell the compiler it doesn't need to be paranoid about writes to one of
//   the buffers screwing with reads from the other.
#ifdef __LP64__
CONSTI32(kPglBitTransposeBufbytes, (kPglBitTransposeBatch * kPglBitTransposeBatch) / (CHAR_BIT / 2));
void TransposeBitblock64(const uintptr_t* read_iter, uintptr_t read_ul_stride, uintptr_t write_ul_stride, uint32_t read_row_ct, uint32_t write_row_ct, uintptr_t* write_iter, VecW* __restrict buf0, VecW* __restrict buf1);

HEADER_INLINE void TransposeBitblock(const uintptr_t* read_iter, uintptr_t read_ul_stride, uintptr_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecW* vecaligned_buf) {
  TransposeBitblock64(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, vecaligned_buf, &(vecaligned_buf[kPglBitTransposeBufbytes / (2 * kBytesPerVec)]));
}

#else  // !__LP64__
CONSTI32(kPglBitTransposeBufbytes, (kPglBitTransposeBatch * kPglBitTransposeBatch) / (CHAR_BIT / 2));
void TransposeBitblock32(const uintptr_t* read_iter, uintptr_t read_ul_stride, uintptr_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecW* __restrict buf0, VecW* __restrict buf1);

// If this ever needs to be called on an input byte array, read_iter could be
// changed to const void*; in that case, read_ul_stride should be changed to a
// byte count.
HEADER_INLINE void TransposeBitblock(const uintptr_t* read_iter, uintptr_t read_ul_stride, uintptr_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecW* vecaligned_buf) {
  TransposeBitblock32(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, vecaligned_buf, &(vecaligned_buf[kPglBitTransposeBufbytes / (2 * kBytesPerVec)]));
}
#endif

CONSTI32(kPglBitTransposeBufwords, kPglBitTransposeBufbytes / kBytesPerWord);
CONSTI32(kPglBitTransposeBufvecs, kPglBitTransposeBufbytes / kBytesPerVec);

CONSTI32(kPglNybbleTransposeBatch, kNybblesPerCacheline);
CONSTI32(kPglNybbleTransposeWords, kWordsPerCacheline);

CONSTI32(kPglNybbleTransposeBufbytes, (kPglNybbleTransposeBatch * kPglNybbleTransposeBatch) / 2);

// up to 128x128; vecaligned_buf must have size 8k
// now ok for write_iter to not be padded when write_batch_size odd
void TransposeNybbleblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, VecW* vecaligned_buf);

#ifdef __LP64__
#  ifdef USE_AVX2
extern const unsigned char kLeadMask[2 * kBytesPerVec] __attribute__ ((aligned (64)));
#  else
extern const unsigned char kLeadMask[2 * kBytesPerVec] __attribute__ ((aligned (32)));
#  endif
#endif

uintptr_t BytesumArr(const void* bytearr, uintptr_t byte_ct);

uintptr_t CountByte(const void* bytearr, unsigned char ucc, uintptr_t byte_ct);

uintptr_t CountU16(const void* u16arr, uint16_t usii, uintptr_t u16_ct);


// Applies sample_include to {src_subset, src_vals}.
uint32_t Copy1bit8Subset(const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uintptr_t* __restrict sample_include, uint32_t src_subset_size, uint32_t sample_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals);

uint32_t Copy1bit16Subset(const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uintptr_t* __restrict sample_include, uint32_t src_subset_size, uint32_t sample_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals);

// more verbose than (val + 3) / 4, but may as well make semantic meaning
// obvious; any explicit DivUp(val, 4) expressions should have a different
// meaning
// (not needed for bitct -> bytect, DivUp(val, CHAR_BIT) is clear enough)
HEADER_INLINE uintptr_t NypCtToByteCt(uintptr_t val) {
  return DivUp(val, 4);
}

HEADER_INLINE uintptr_t NypCtToVecCt(uintptr_t val) {
  return DivUp(val, kNypsPerVec);
}

HEADER_INLINE uintptr_t NypCtToWordCt(uintptr_t val) {
  return DivUp(val, kBitsPerWordD2);
}

HEADER_INLINE uintptr_t NypCtToAlignedWordCt(uintptr_t val) {
  return kWordsPerVec * NypCtToVecCt(val);
}

HEADER_INLINE uintptr_t NypCtToCachelineCt(uintptr_t val) {
  return DivUp(val, kNypsPerCacheline);
}

HEADER_INLINE uintptr_t GetNyparrEntry(const uintptr_t* nyparr, uint32_t idx) {
  return (nyparr[idx / kBitsPerWordD2] >> (2 * (idx % kBitsPerWordD2))) & 3;
}

// todo: check if this optimizes newval=0 out
HEADER_INLINE void AssignNyparrEntry(uint32_t idx, uintptr_t newval, uintptr_t* nyparr) {
  const uint32_t bit_shift_ct = 2 * (idx % kBitsPerWordD2);
  uintptr_t* wordp = &(nyparr[idx / kBitsPerWordD2]);
  *wordp = ((*wordp) & (~((3 * k1LU) << bit_shift_ct))) | (newval << bit_shift_ct);
}

HEADER_INLINE void ClearNyparrEntry(uint32_t idx, uintptr_t* nyparr) {
  nyparr[idx / kBitsPerWordD2] &= ~((3 * k1LU) << (idx % kBitsPerWordD2));
}

// Assumes arr is vector-aligned.
// 'Unsafe' because it assumes high bits of every byte are 0 and entry_ct is
// positive.
void Reduce8to4bitInplaceUnsafe(uintptr_t entry_ct, uintptr_t* arr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_BITS_H__
