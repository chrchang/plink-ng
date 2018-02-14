// This library is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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


#include "pgenlib_internal.h"

#ifndef NO_MMAP
#  include <sys/types.h>  // fstat()
#  include <sys/stat.h>  // open(), fstat()
#  include <sys/mman.h>  // mmap()
#  include <fcntl.h>  // open()
#  include <unistd.h>  // fstat()
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

#ifdef USE_AVX2
void CopyQuaterarrNonemptySubset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_quaterarr) {
  if (subset_entry_ct == raw_quaterarr_entry_ct) {
    memcpy(output_quaterarr, raw_quaterarr, DivUp(subset_entry_ct, kBitsPerWordD2) * sizeof(intptr_t));
    ZeroTrailingQuaters(subset_entry_ct, output_quaterarr);
    return;
  }
  assert(subset_entry_ct);
  uintptr_t cur_output_word = 0;

  uintptr_t* output_quaterarr_iter = output_quaterarr;

  uintptr_t* output_quaterarr_last = &(output_quaterarr[subset_entry_ct / kBitsPerWordD2]);
  const uint32_t word_write_shift_end = 2 * (subset_entry_ct % kBitsPerWordD2);
  uint32_t word_write_shift = 0;
  uint32_t subset_mask_widx = 0;
  while (1) {
    const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
    if (cur_include_word) {
      uint32_t wordhalf_idx = 0;
      uint32_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
      while (1) {
        if (cur_include_halfword) {
          uintptr_t extracted_bits = raw_quaterarr[subset_mask_widx * 2 + wordhalf_idx];
          uint32_t set_bit_ct = kBitsPerWord;
          if (cur_include_halfword != UINT32_MAX) {
            const uintptr_t pext_mask = 3 * UnpackHalfwordToWord(cur_include_halfword);
            extracted_bits = _pext_u64(extracted_bits, pext_mask);
            set_bit_ct = PopcountWord(pext_mask);
          }
          cur_output_word |= extracted_bits << word_write_shift;
          word_write_shift += set_bit_ct;
          if (word_write_shift >= kBitsPerWord) {
            *output_quaterarr_iter++ = cur_output_word;
            word_write_shift -= kBitsPerWord;
            cur_output_word = 0;
            if (word_write_shift) {
              cur_output_word = extracted_bits >> (set_bit_ct - word_write_shift);
            }
          }
        }
        if (wordhalf_idx) {
          break;
        }
        ++wordhalf_idx;
        cur_include_halfword = cur_include_word >> kBitsPerWordD2;
      }
      if (output_quaterarr_iter == output_quaterarr_last) {
        if (word_write_shift == word_write_shift_end) {
          if (word_write_shift_end) {
            *output_quaterarr_last = cur_output_word;
          }
          return;
        }
      }
    }
    ++subset_mask_widx;
  }
}
#else  // !USE_AVX2
void CopyQuaterarrNonemptySubset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_quaterarr) {
  // in plink 2.0, we probably want (0-based) bit raw_quaterarr_entry_ct of
  // subset_mask to be always allocated and unset.  This removes a few special
  // cases re: iterating past the end of arrays.
  assert(subset_entry_ct);
  assert(raw_quaterarr_entry_ct >= subset_entry_ct);
  uintptr_t cur_output_word = 0;

  uintptr_t* output_quaterarr_iter = output_quaterarr;

  uintptr_t* output_quaterarr_last = &(output_quaterarr[subset_entry_ct / kBitsPerWordD2]);
  const uint32_t word_write_halfshift_end = subset_entry_ct % kBitsPerWordD2;
  uint32_t word_write_halfshift = 0;
  // if <= 2/3-filled, use sparse copy algorithm
  // (tried CopyBitarrSubset() approach, that actually worsened things)
  if (subset_entry_ct * (3 * k1LU) <= raw_quaterarr_entry_ct * (2 * k1LU)) {
    uint32_t subset_mask_widx = 0;
    while (1) {
      const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
      if (cur_include_word) {
        uint32_t wordhalf_idx = 0;
        uint32_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
        while (1) {
          if (cur_include_halfword) {
            uintptr_t raw_quaterarr_word = raw_quaterarr[subset_mask_widx * 2 + wordhalf_idx];
            do {
              uint32_t rqa_idx_lowbits = __builtin_ctz(cur_include_halfword);
              cur_output_word |= ((raw_quaterarr_word >> (rqa_idx_lowbits * 2)) & 3) << (word_write_halfshift * 2);
              if (++word_write_halfshift == kBitsPerWordD2) {
                *output_quaterarr_iter++ = cur_output_word;
                word_write_halfshift = 0;
                cur_output_word = 0;
              }
              cur_include_halfword &= cur_include_halfword - 1;
            } while (cur_include_halfword);
          }
          if (wordhalf_idx) {
            break;
          }
          ++wordhalf_idx;
          cur_include_halfword = cur_include_word >> kBitsPerWordD2;
        }
        if (output_quaterarr_iter == output_quaterarr_last) {
          if (word_write_halfshift == word_write_halfshift_end) {
            if (word_write_halfshift_end) {
              *output_quaterarr_last = cur_output_word;
            }
            return;
          }
        }
      }
      ++subset_mask_widx;
    }
  }
  // blocked copy
  const uintptr_t* raw_quaterarr_iter = raw_quaterarr;
  while (1) {
    const uintptr_t cur_include_word = *subset_mask++;
    uint32_t wordhalf_idx = 0;
    uintptr_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
    while (1) {
      uintptr_t raw_quaterarr_word = *raw_quaterarr_iter++;
      while (cur_include_halfword) {
        uint32_t rqa_idx_lowbits = CTZLU(cur_include_halfword);
        uintptr_t halfword_invshifted = (~cur_include_halfword) >> rqa_idx_lowbits;
        uintptr_t raw_quaterarr_curblock_unmasked = raw_quaterarr_word >> (rqa_idx_lowbits * 2);
        uint32_t rqa_block_len = CTZLU(halfword_invshifted);
        uint32_t block_len_limit = kBitsPerWordD2 - word_write_halfshift;
        cur_output_word |= raw_quaterarr_curblock_unmasked << (2 * word_write_halfshift);
        if (rqa_block_len < block_len_limit) {
          word_write_halfshift += rqa_block_len;
          cur_output_word &= (k1LU << (word_write_halfshift * 2)) - k1LU;
        } else {
          // no need to mask, extra bits vanish off the high end
          *output_quaterarr_iter++ = cur_output_word;
          word_write_halfshift = rqa_block_len - block_len_limit;
          if (word_write_halfshift) {
            cur_output_word = (raw_quaterarr_curblock_unmasked >> (2 * block_len_limit)) & ((k1LU << (2 * word_write_halfshift)) - k1LU);
          } else {
            // avoid potential right-shift-[word length]
            cur_output_word = 0;
          }
        }
        cur_include_halfword &= (~(k1LU << (rqa_block_len + rqa_idx_lowbits))) + k1LU;
      }
      if (wordhalf_idx) {
        break;
      }
      ++wordhalf_idx;
      cur_include_halfword = cur_include_word >> kBitsPerWordD2;
    }
    if (output_quaterarr_iter == output_quaterarr_last) {
      if (word_write_halfshift == word_write_halfshift_end) {
        if (word_write_halfshift_end) {
          *output_quaterarr_last = cur_output_word;
        }
        return;
      }
    }
  }
}
#endif

// Harley-Seal algorithm only works for bitarrays, not quaterarrays, so don't
// add an AVX2 specialization here.
// ...unless something like the interleaved_vec strategy is used?  hmm.  should
// test this on basic frequency counter.
void Count2FreqVec3(const VecUL* geno_vvec, uint32_t vec_ct, uint32_t* __restrict alt1_plus_bothset_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 3));
  // Increments bothset_ct by the number of 0b11 in the current block, and
  // alt1_ct by twice the number of 0b10 plus the number of 0b01.
  const VecUL m1 = VCONST_UL(kMask5555);
  const VecUL m2 = VCONST_UL(kMask3333);
  const VecUL m4 = VCONST_UL(kMask0F0F);
  const VecUL* geno_vvec_iter = geno_vvec;
  uint32_t alt1_plus_bothset_ct = 0;
  uint32_t bothset_ct = 0;

  while (1) {
    UniVec acc_alt1_plus_bothset;
    UniVec acc_bothset;
    acc_alt1_plus_bothset.vi = vecul_setzero();
    acc_bothset.vi = vecul_setzero();
    const VecUL* geno_vvec_stop;
    if (vec_ct < 30) {
      if (!vec_ct) {
        *alt1_plus_bothset_ctp = alt1_plus_bothset_ct;
        *bothset_ctp = bothset_ct;
        return;
      }
      geno_vvec_stop = &(geno_vvec_iter[vec_ct]);
      vec_ct = 0;
    } else {
      geno_vvec_stop = &(geno_vvec_iter[30]);
      vec_ct -= 30;
    }
    do {
      VecUL cur_geno_vword1 = *geno_vvec_iter++;
      // process first two vwords simultaneously to minimize linear dependence
      VecUL cur_geno_vword2 = *geno_vvec_iter++;
      VecUL cur_geno_vword_low_lshifted1 = vecul_slli(cur_geno_vword1 & m1, 1);
      VecUL cur_geno_vword_low_lshifted2 = vecul_slli(cur_geno_vword2 & m1, 1);

      // 00 -> 00; 01 -> 01; 10 -> 10; 11 -> 01
      // note that _mm_andnot_si128 flips the *first* argument before the AND
      // operation.
      VecUL alt1_plus_bothset1 = (~cur_geno_vword_low_lshifted1) & cur_geno_vword1;
      VecUL alt1_plus_bothset2 = (~cur_geno_vword_low_lshifted2) & cur_geno_vword2;

      VecUL bothset1 = vecul_srli(cur_geno_vword_low_lshifted1 & cur_geno_vword1, 1);
      VecUL bothset2 = vecul_srli(cur_geno_vword_low_lshifted2 & cur_geno_vword2, 1);

      cur_geno_vword1 = *geno_vvec_iter++;
      alt1_plus_bothset1 = (alt1_plus_bothset1 & m2) + (vecul_srli(alt1_plus_bothset1, 2) & m2);
      bothset2 = bothset1 + bothset2;
      alt1_plus_bothset2 = (alt1_plus_bothset2 & m2) + (vecul_srli(alt1_plus_bothset2, 2) & m2);
      cur_geno_vword_low_lshifted1 = vecul_slli(cur_geno_vword1 & m1, 1);

      alt1_plus_bothset2 = alt1_plus_bothset1 + alt1_plus_bothset2;
      // alt1_plus_bothset2 now contains 4-bit values from 0-8, while bothset2
      // contains 2-bit values from 0-2
      // (todo: check whether this is faster if we use double_bothsetx
      // variables instead of bothset1/bothset2)
      bothset1 = vecul_srli(cur_geno_vword_low_lshifted1 & cur_geno_vword1, 1);
      alt1_plus_bothset1 = (~cur_geno_vword_low_lshifted1) & cur_geno_vword1;
      bothset2 = bothset1 + bothset2;
      alt1_plus_bothset1 = (alt1_plus_bothset1 & m2) + (vecul_srli(alt1_plus_bothset1, 2) & m2);

      bothset2 = (bothset2 & m2) + (vecul_srli(bothset2, 2) & m2);
      alt1_plus_bothset2 = alt1_plus_bothset1 + alt1_plus_bothset2;
      // alt1_plus_bothset2 now contains 4-bit values from 0-12, while bothset2
      // contains 4-bit values from 0-6.  aggregate both into 8-bit values.
      bothset2 = (bothset2 & m4) + (vecul_srli(bothset2, 4) & m4);
      alt1_plus_bothset2 = (alt1_plus_bothset2 & m4) + (vecul_srli(alt1_plus_bothset2, 4) & m4);

      acc_bothset.vi = acc_bothset.vi + bothset2;
      acc_alt1_plus_bothset.vi = acc_alt1_plus_bothset.vi + alt1_plus_bothset2;
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecUL m8 = VCONST_UL(kMask00FF);
    acc_bothset.vi = (acc_bothset.vi + vecul_srli(acc_bothset.vi, 8)) & m8;
    acc_alt1_plus_bothset.vi = (acc_alt1_plus_bothset.vi & m8) + (vecul_srli(acc_alt1_plus_bothset.vi, 8) & m8);
    bothset_ct += UniVecHsum16(acc_bothset);
    alt1_plus_bothset_ct += UniVecHsum16(acc_alt1_plus_bothset);
  }
}

void Count3FreqVec6(const VecUL* geno_vvec, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 6));
  // Sets even_ct to the number of set low bits in the current block, odd_ct to
  // the number of set high bits, and bothset_ct by the number of 0b11s.
  // Easy to adapt this to take a subset quatervec parameter.
  const VecUL m1 = VCONST_UL(kMask5555);
  const VecUL m2 = VCONST_UL(kMask3333);
  const VecUL m4 = VCONST_UL(kMask0F0F);
  const VecUL* geno_vvec_iter = geno_vvec;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  while (1) {
    UniVec acc_even;
    UniVec acc_odd;
    UniVec acc_bothset;
    acc_even.vi = vecul_setzero();
    acc_odd.vi = vecul_setzero();
    acc_bothset.vi = vecul_setzero();
    const VecUL* geno_vvec_stop;
    if (vec_ct < 60) {
      if (!vec_ct) {
        *even_ctp = even_ct;
        *odd_ctp = odd_ct;
        *bothset_ctp = bothset_ct;
        return;
      }
      geno_vvec_stop = &(geno_vvec_iter[vec_ct]);
      vec_ct = 0;
    } else {
      geno_vvec_stop = &(geno_vvec_iter[60]);
      vec_ct -= 60;
    }
    do {
      // hmm, this seems to have more linear dependence than I'd want, but the
      // reorderings I tried just made the code harder to read without helping,
      // so I'll leave this alone
      VecUL cur_geno_vword = *geno_vvec_iter++;
      VecUL odd1 = m1 & vecul_srli(cur_geno_vword, 1);
      VecUL even1 = m1 & cur_geno_vword;
      VecUL bothset1 = odd1 & cur_geno_vword;

      cur_geno_vword = *geno_vvec_iter++;
      VecUL cur_geno_vword_high = m1 & vecul_srli(cur_geno_vword, 1);
      even1 = even1 + (m1 & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high;
      bothset1 = bothset1 + (cur_geno_vword_high & cur_geno_vword);

      cur_geno_vword = *geno_vvec_iter++;
      cur_geno_vword_high = m1 & vecul_srli(cur_geno_vword, 1);
      even1 = even1 + (m1 & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high;
      bothset1 = bothset1 + (cur_geno_vword_high & cur_geno_vword);

      even1 = (even1 & m2) + (vecul_srli(even1, 2) & m2);
      odd1 = (odd1 & m2) + (vecul_srli(odd1, 2) & m2);
      bothset1 = (bothset1 & m2) + (vecul_srli(bothset1, 2) & m2);

      cur_geno_vword = *geno_vvec_iter++;
      VecUL odd2 = m1 & vecul_srli(cur_geno_vword, 1);
      VecUL even2 = m1 & cur_geno_vword;
      VecUL bothset2 = odd2 & cur_geno_vword;

      cur_geno_vword = *geno_vvec_iter++;
      cur_geno_vword_high = m1 & vecul_srli(cur_geno_vword, 1);
      even2 = even2 + (m1 & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high;
      bothset2 = bothset2 + (cur_geno_vword_high & cur_geno_vword);

      cur_geno_vword = *geno_vvec_iter++;
      cur_geno_vword_high = m1 & vecul_srli(cur_geno_vword, 1);
      even2 = even2 + (m1 & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high;
      bothset2 = bothset2 + (cur_geno_vword_high & cur_geno_vword);

      even1 = even1 + (even2 & m2) + (vecul_srli(even2, 2) & m2);
      odd1 = odd1 + (odd2 & m2) + (vecul_srli(odd2, 2) & m2);
      bothset1 = bothset1 + (bothset2 & m2) + (vecul_srli(bothset2, 2) & m2);
      // these now contain 4-bit values from 0-12

      acc_even.vi = acc_even.vi + (even1 & m4) + (vecul_srli(even1, 4) & m4);
      acc_odd.vi = acc_odd.vi + (odd1 & m4) + (vecul_srli(odd1, 4) & m4);
      acc_bothset.vi = acc_bothset.vi + (bothset1 & m4) + (vecul_srli(bothset1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecUL m8 = VCONST_UL(kMask00FF);
    acc_even.vi = (acc_even.vi & m8) + (vecul_srli(acc_even.vi, 8) & m8);
    acc_odd.vi = (acc_odd.vi & m8) + (vecul_srli(acc_odd.vi, 8) & m8);
    acc_bothset.vi = (acc_bothset.vi & m8) + (vecul_srli(acc_bothset.vi, 8) & m8);
    even_ct += UniVecHsum16(acc_even);
    odd_ct += UniVecHsum16(acc_odd);
    bothset_ct += UniVecHsum16(acc_bothset);
  }
}

void CountSubset3FreqVec6(const VecUL* __restrict geno_vvec, const VecUL* __restrict interleaved_mask_vvec, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 6));
  // Sets even_ct to the number of set low bits in the current block, odd_ct to
  // the number of set high bits, and bothset_ct by the number of 0b11s.
  // Easy to adapt this to take a subset quatervec parameter.
  const VecUL m1 = VCONST_UL(kMask5555);
  const VecUL m2 = VCONST_UL(kMask3333);
  const VecUL m4 = VCONST_UL(kMask0F0F);
  const VecUL* geno_vvec_iter = geno_vvec;
  const VecUL* interleaved_mask_vvec_iter = interleaved_mask_vvec;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  while (1) {
    UniVec acc_even;
    UniVec acc_odd;
    UniVec acc_bothset;
    acc_even.vi = vecul_setzero();
    acc_odd.vi = vecul_setzero();
    acc_bothset.vi = vecul_setzero();
    const VecUL* geno_vvec_stop;
    if (vec_ct < 60) {
      if (!vec_ct) {
        *even_ctp = even_ct;
        *odd_ctp = odd_ct;
        *bothset_ctp = bothset_ct;
        return;
      }
      geno_vvec_stop = &(geno_vvec_iter[vec_ct]);
      vec_ct = 0;
    } else {
      geno_vvec_stop = &(geno_vvec_iter[60]);
      vec_ct -= 60;
    }
    do {
      VecUL interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      VecUL cur_geno_vword = *geno_vvec_iter++;
      VecUL cur_mask = interleaved_mask_vword & m1;
      VecUL odd1 = cur_mask & vecul_srli(cur_geno_vword, 1);
      VecUL even1 = cur_mask & cur_geno_vword;
      VecUL bothset1 = odd1 & cur_geno_vword;

      cur_mask = vecul_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = *geno_vvec_iter++;
      VecUL cur_geno_vword_high_masked = cur_mask & vecul_srli(cur_geno_vword, 1);
      even1 = even1 + (cur_mask & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high_masked;
      bothset1 = bothset1 + (cur_geno_vword_high_masked & cur_geno_vword);

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = *geno_vvec_iter++;
      cur_geno_vword_high_masked = cur_mask & vecul_srli(cur_geno_vword, 1);
      even1 = even1 + (cur_mask & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high_masked;
      bothset1 = bothset1 + (cur_geno_vword_high_masked & cur_geno_vword);

      even1 = (even1 & m2) + (vecul_srli(even1, 2) & m2);
      odd1 = (odd1 & m2) + (vecul_srli(odd1, 2) & m2);
      bothset1 = (bothset1 & m2) + (vecul_srli(bothset1, 2) & m2);

      cur_mask = vecul_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = *geno_vvec_iter++;
      VecUL odd2 = cur_mask & vecul_srli(cur_geno_vword, 1);
      VecUL even2 = cur_mask & cur_geno_vword;
      VecUL bothset2 = odd2 & cur_geno_vword;

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = *geno_vvec_iter++;
      cur_geno_vword_high_masked = cur_mask & vecul_srli(cur_geno_vword, 1);
      even2 = even2 + (cur_mask & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high_masked;
      bothset2 = bothset2 + (cur_geno_vword_high_masked & cur_geno_vword);

      cur_mask = vecul_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = *geno_vvec_iter++;
      cur_geno_vword_high_masked = cur_mask & vecul_srli(cur_geno_vword, 1);
      even2 = even2 + (cur_mask & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high_masked;
      bothset2 = bothset2 + (cur_geno_vword_high_masked & cur_geno_vword);

      even1 = even1 + (even2 & m2) + (vecul_srli(even2, 2) & m2);
      odd1 = odd1 + (odd2 & m2) + (vecul_srli(odd2, 2) & m2);
      bothset1 = bothset1 + (bothset2 & m2) + (vecul_srli(bothset2, 2) & m2);
      // these now contain 4-bit values from 0-12

      acc_even.vi = acc_even.vi + (even1 & m4) + (vecul_srli(even1, 4) & m4);
      acc_odd.vi = acc_odd.vi + (odd1 & m4) + (vecul_srli(odd1, 4) & m4);
      acc_bothset.vi = acc_bothset.vi + (bothset1 & m4) + (vecul_srli(bothset1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecUL m8 = VCONST_UL(kMask00FF);
    acc_even.vi = (acc_even.vi & m8) + (vecul_srli(acc_even.vi, 8) & m8);
    acc_odd.vi = (acc_odd.vi & m8) + (vecul_srli(acc_odd.vi, 8) & m8);
    acc_bothset.vi = (acc_bothset.vi & m8) + (vecul_srli(acc_bothset.vi, 8) & m8);
    even_ct += UniVecHsum16(acc_even);
    odd_ct += UniVecHsum16(acc_odd);
    bothset_ct += UniVecHsum16(acc_bothset);
  }
}

// possible todo: AVX2
uint32_t Count01Vec6(const VecUL* geno_vvec, uint32_t vec_ct) {
  assert(!(vec_ct % 6));
  const VecUL m1 = VCONST_UL(kMask5555);
  const VecUL m2 = VCONST_UL(kMask3333);
  const VecUL m4 = VCONST_UL(kMask0F0F);
  const VecUL m8 = VCONST_UL(kMask00FF);
  const VecUL* geno_vvec_iter = geno_vvec;
  uint32_t tot = 0;
  while (1) {
    UniVec acc;
    acc.vi = vecul_setzero();
    const VecUL* geno_vvec_stop;
    if (vec_ct < 60) {
      if (!vec_ct) {
        return tot;
      }
      geno_vvec_stop = &(geno_vvec_iter[vec_ct]);
      vec_ct = 0;
    } else {
      geno_vvec_stop = &(geno_vvec_iter[60]);
      vec_ct -= 60;
    }
    do {
      VecUL loader1 = *geno_vvec_iter++;
      VecUL loader2 = *geno_vvec_iter++;
      VecUL count1 = ((~vecul_srli(loader1, 1)) & loader1) & m1;
      VecUL count2 = ((~vecul_srli(loader2, 1)) & loader2) & m1;

      loader1 = *geno_vvec_iter++;
      loader2 = *geno_vvec_iter++;
      count1 = count1 + (((~vecul_srli(loader1, 1)) & loader1) & m1);
      count2 = count2 + (((~vecul_srli(loader2, 1)) & loader2) & m1);

      loader1 = *geno_vvec_iter++;
      loader2 = *geno_vvec_iter++;
      count1 = count1 + (((~vecul_srli(loader1, 1)) & loader1) & m1);
      count2 = count2 + (((~vecul_srli(loader2, 1)) & loader2) & m1);

      count1 = (count1 & m2) + (vecul_srli(count1, 2) & m2);
      count1 = count1 + (count2 & m2) + (vecul_srli(count2, 2) & m2);
      acc.vi = acc.vi + (count1 & m4) + (vecul_srli(count1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    acc.vi = (acc.vi & m8) + (vecul_srli(acc.vi, 8) & m8);
    tot += UniVecHsum16(acc);
  }
}

// This type of code can probably be improved; see
//   https://lemire.me/blog/2018/01/09/how-fast-can-you-bit-interleave-32-bit-integers-simd-edition/ .
// (Also, it may be worthwhile to add a compilation option which generates AVX2
// instructions but mostly avoids _pdep_u64() and _pext_u64(); apparently even
// the latest AMD processors have poor implementations of those two
// operations?)
void FillInterleavedMaskVec(const uintptr_t* __restrict subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec) {
#ifdef __LP64__
  const uintptr_t* subset_mask_iter = subset_mask;
  uintptr_t* interleaved_mask_vec_iter = interleaved_mask_vec;
#  ifdef USE_AVX2
  uintptr_t orig_word1 = 0;
  uintptr_t orig_word3 = 0;
#  endif
  for (uint32_t vec_idx = 0; vec_idx < base_vec_ct; ++vec_idx) {
#  ifdef USE_AVX2
    // 0 128 1 129 2 130 ...
    for (uint32_t widx = 0; widx < 4; ++widx) {
      uintptr_t ww_even;
      uintptr_t ww_odd;
      if (!(widx % 2)) {
        orig_word1 = subset_mask_iter[0];
        orig_word3 = subset_mask_iter[2];
        ++subset_mask_iter;
        ww_even = S_CAST(uint32_t, orig_word1);
        ww_odd = S_CAST(uint32_t, orig_word3);
      } else {
        ww_even = orig_word1 >> 32;
        ww_odd = orig_word3 >> 32;
      }
      ww_even = UnpackHalfwordToWord(ww_even);
      ww_odd = UnpackHalfwordToWord(ww_odd);
      *interleaved_mask_vec_iter++ = ww_even | (ww_odd << 1);
    }
    subset_mask_iter = &(subset_mask_iter[2]);
#  else  // not USE_AVX2
    // 0 64 1 65 2 66 ...
    const uintptr_t orig_word1 = *subset_mask_iter++;
    const uintptr_t orig_word2 = *subset_mask_iter++;
    for (uint32_t widx = 0; widx < 2; ++widx) {
      uintptr_t ww_even;
      uintptr_t ww_odd;
      // todo: check if there's a better way to organize this loop
      if (!widx) {
        ww_even = S_CAST(uint32_t, orig_word1);
        ww_odd = S_CAST(uint32_t, orig_word2);
      } else {
        ww_even = orig_word1 >> 32;
        ww_odd = orig_word2 >> 32;
      }
      ww_even = UnpackHalfwordToWord(ww_even);
      ww_odd = UnpackHalfwordToWord(ww_odd);
      *interleaved_mask_vec_iter++ = ww_even | (ww_odd << 1);
    }
#  endif  // not USE_AVX2
  }
#else
  for (uint32_t widx = 0; widx < base_vec_ct; ++widx) {
    const uintptr_t orig_word = subset_mask[widx];
    uintptr_t ww_even = S_CAST(uint16_t, orig_word);
    uintptr_t ww_odd = orig_word >> 16;
    ww_even = UnpackHalfwordToWord(ww_even);
    ww_odd = UnpackHalfwordToWord(ww_odd);
    interleaved_mask_vec[widx] = ww_even | (ww_odd << 1);
  }
#endif
}

void GenovecAlleleCtsUnsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* __restrict allele_cts, uint32_t* __restrict bothset_ctp) {
  // assumes trailing bits of last genovec word are zeroed out.
  // sets allele_cts[0] to the number of observed ref alleles, and
  // allele_cts[1] to the number of observed alt1s.
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (3 * kWordsPerVec));
  uint32_t alt1_plus_bothset_ct;
  uint32_t bothset_ct;
  assert(IsVecAligned(genovec));
  Count2FreqVec3(R_CAST(const VecUL*, genovec), word_idx / kWordsPerVec, &alt1_plus_bothset_ct, &bothset_ct);
  for (; word_idx < sample_ctl2; ++word_idx) {
    const uintptr_t cur_geno_word = genovec[word_idx];
    const uintptr_t cur_geno_word_low_lshifted = (cur_geno_word & kMask5555) << 1;
    alt1_plus_bothset_ct += QuatersumWord((~cur_geno_word_low_lshifted) & cur_geno_word);
    bothset_ct += QuatersumWord(cur_geno_word_low_lshifted & cur_geno_word);
  }
  const uint32_t alt1_ct = alt1_plus_bothset_ct - bothset_ct;
  allele_cts[0] = (sample_ct - bothset_ct) * 2 - alt1_ct;
  allele_cts[1] = alt1_ct;
  *bothset_ctp = bothset_ct;
}

void GenovecCountFreqsUnsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* genocounts) {
  // fills genocounts[0] with the number of 00s, genocounts[1] with the number
  // of 01s, etc.
  // assumes trailing bits of last genovec word are zeroed out.
  // sample_ct == 0 ok.
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  uint32_t even_ct;
  uint32_t odd_ct;
  uint32_t bothset_ct;
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (6 * kWordsPerVec));
  assert(IsVecAligned(genovec));
  Count3FreqVec6(R_CAST(const VecUL*, genovec), word_idx / kWordsPerVec, &even_ct, &odd_ct, &bothset_ct);
  for (; word_idx < sample_ctl2; ++word_idx) {
    const uintptr_t cur_geno_word = genovec[word_idx];
    const uintptr_t cur_geno_word_high = kMask5555 & (cur_geno_word >> 1);
    even_ct += Popcount01Word(cur_geno_word & kMask5555);
    odd_ct += Popcount01Word(cur_geno_word_high);
    bothset_ct += Popcount01Word(cur_geno_word & cur_geno_word_high);
  }
  genocounts[0] = sample_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}

void GenovecCountSubsetFreqs(const uintptr_t* __restrict genovec, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t* genocounts) {
  // fills genocounts[0] with the number of 00s, genocounts[1] with the number
  // of 01s, etc.
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctv2 = QuaterCtToVecCt(raw_sample_ct);
  uint32_t even_ct;
  uint32_t odd_ct;
  uint32_t bothset_ct;
#ifdef __LP64__
  uint32_t vec_idx = raw_sample_ctv2 - (raw_sample_ctv2 % 6);
  assert(IsVecAligned(genovec));
  CountSubset3FreqVec6(R_CAST(const VecUL*, genovec), R_CAST(const VecUL*, sample_include_interleaved_vec), vec_idx, &even_ct, &odd_ct, &bothset_ct);
  const uintptr_t* genovec_iter = &(genovec[kWordsPerVec * vec_idx]);
  const uintptr_t* interleaved_mask_iter = &(sample_include_interleaved_vec[(kWordsPerVec / 2) * vec_idx]);
#  ifdef USE_AVX2
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  uintptr_t mask_base3 = 0;
  uintptr_t mask_base4 = 0;
  for (; vec_idx < raw_sample_ctv2; ++vec_idx) {
    uintptr_t mask_word1;
    uintptr_t mask_word2;
    uintptr_t mask_word3;
    uintptr_t mask_word4;
    if (!(vec_idx % 2)) {
      mask_base1 = *interleaved_mask_iter++;
      mask_base2 = *interleaved_mask_iter++;
      mask_base3 = *interleaved_mask_iter++;
      mask_base4 = *interleaved_mask_iter++;
      mask_word1 = mask_base1 & kMask5555;
      mask_word2 = mask_base2 & kMask5555;
      mask_word3 = mask_base3 & kMask5555;
      mask_word4 = mask_base4 & kMask5555;
    } else {
      mask_word1 = (mask_base1 >> 1) & kMask5555;
      mask_word2 = (mask_base2 >> 1) & kMask5555;
      mask_word3 = (mask_base3 >> 1) & kMask5555;
      mask_word4 = (mask_base4 >> 1) & kMask5555;
    }
    uint32_t uii = 0;
    while (1) {
      const uintptr_t cur_geno_word1 = *genovec_iter++;
      const uintptr_t cur_geno_word2 = *genovec_iter++;
      const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
      const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
      even_ct += PopcountWord(((cur_geno_word1 & mask_word1) << 1) | (cur_geno_word2 & mask_word2));
      odd_ct += PopcountWord((cur_geno_word1_high_masked << 1) | cur_geno_word2_high_masked);
      bothset_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
      if (uii) {
        break;
      }
      ++uii;
      mask_word1 = mask_word3;
      mask_word2 = mask_word4;
    }
  }
#  else  // not USE_AVX2
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  for (; vec_idx < raw_sample_ctv2; ++vec_idx) {
    uintptr_t mask_word1;
    uintptr_t mask_word2;
    if (!(vec_idx % 2)) {
      mask_base1 = *interleaved_mask_iter++;
      mask_base2 = *interleaved_mask_iter++;
      mask_word1 = mask_base1 & kMask5555;
      mask_word2 = mask_base2 & kMask5555;
    } else {
      mask_word1 = (mask_base1 >> 1) & kMask5555;
      mask_word2 = (mask_base2 >> 1) & kMask5555;
    }
    const uintptr_t cur_geno_word1 = *genovec_iter++;
    const uintptr_t cur_geno_word2 = *genovec_iter++;
    const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
    const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
#    ifdef USE_SSE42
    even_ct += PopcountWord(((cur_geno_word1 & mask_word1) << 1) | (cur_geno_word2 & mask_word2));
    odd_ct += PopcountWord((cur_geno_word1_high_masked << 1) | cur_geno_word2_high_masked);
    bothset_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
#    else
    even_ct += QuatersumWord((cur_geno_word1 & mask_word1) + (cur_geno_word2 & mask_word2));
    odd_ct += QuatersumWord(cur_geno_word1_high_masked + cur_geno_word2_high_masked);
    bothset_ct += QuatersumWord((cur_geno_word1 & cur_geno_word1_high_masked) + (cur_geno_word2 & cur_geno_word2_high_masked));
#    endif
  }
#  endif  // not USE_AVX2
#else  // not __LP64__
  uint32_t word_idx = raw_sample_ctv2 - (raw_sample_ctv2 % 6);
  CountSubset3FreqVec6(R_CAST(const VecUL*, genovec), R_CAST(const VecUL*, sample_include_interleaved_vec), word_idx, &even_ct, &odd_ct, &bothset_ct);
  const uintptr_t* interleaved_mask_iter = &(sample_include_interleaved_vec[word_idx / 2]);
  uintptr_t mask_base = 0;
  for (; word_idx < raw_sample_ctv2; ++word_idx) {
    uintptr_t mask_word;
    if (!(word_idx % 2)) {
      mask_base = *interleaved_mask_iter++;
      mask_word = mask_base & kMask5555;
    } else {
      mask_word = (mask_base >> 1) & kMask5555;
    }
    const uintptr_t cur_geno_word = genovec[word_idx];
    const uintptr_t cur_geno_word_high_masked = mask_word & (cur_geno_word >> 1);
    even_ct += Popcount01Word(cur_geno_word & mask_word);
    odd_ct += Popcount01Word(cur_geno_word_high_masked);
    bothset_ct += Popcount01Word(cur_geno_word & cur_geno_word_high_masked);
  }
#endif
  genocounts[0] = sample_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}

uint32_t GenovecCount01Unsafe(const uintptr_t* genovec, uint32_t sample_ct) {
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (6 * kWordsPerVec));
  assert(IsVecAligned(genovec));
  uint32_t tot = Count01Vec6(R_CAST(const VecUL*, genovec), word_idx / kWordsPerVec);
  for (; word_idx < sample_ctl2; ++word_idx) {
    const uintptr_t cur_geno_word = genovec[word_idx];
    tot += Popcount01Word(cur_geno_word & (~(cur_geno_word >> 1)) & kMask5555);
  }
  return tot;
}

#ifdef __arm__
#  error "Unaligned accesses in SmallGenoarrCount3FreqIncr()."
#endif
void SmallGenoarrCount3FreqIncr(const uintptr_t* genoarr_iter, uint32_t byte_ct, uint32_t* even_ctp, uint32_t* odd_ctp, uint32_t* bothset_ctp) {
  while (1) {
    uintptr_t cur_geno_word;
    if (byte_ct < kBytesPerWord) {
      if (!byte_ct) {
        return;
      }
      cur_geno_word = ProperSubwordLoad(genoarr_iter, byte_ct);
      byte_ct = 0;
    } else {
      cur_geno_word = *genoarr_iter++;
      byte_ct -= kBytesPerWord;
    }
    const uintptr_t cur_geno_word_high = kMask5555 & (cur_geno_word >> 1);
    *even_ctp += Popcount01Word(cur_geno_word & kMask5555);
    *odd_ctp += Popcount01Word(cur_geno_word_high);
    *bothset_ctp += Popcount01Word(cur_geno_word & cur_geno_word_high);
  }
}

void GenoarrCountFreqs(const unsigned char* genoarr, uint32_t sample_ct, uint32_t* genocounts) {
  // does not read past the end of genoarr
  uint32_t lead_byte_ct = (-R_CAST(uintptr_t, genoarr)) % kBytesPerVec;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  const uintptr_t* genoarr_iter;
  uint32_t trail_ct;
  if (sample_ct > lead_byte_ct * 4 + (6 * kQuatersPerVec)) {
    const uint32_t remaining_sample_ct = sample_ct - 4 * lead_byte_ct;
    // strictly speaking, this relies on undefined behavior: see e.g.
    // http://pzemtsov.github.io/2016/11/06/bug-story-alignment-on-x86.html
    // probably want to search out all instances of __arm__ and make the code
    // standard-compliant, if that can be done without a speed penalty
    SmallGenoarrCount3FreqIncr(R_CAST(const uintptr_t*, genoarr), lead_byte_ct, &even_ct, &odd_ct, &bothset_ct);
    genoarr_iter = R_CAST(const uintptr_t*, &(genoarr[lead_byte_ct]));
    const uint32_t remaining_full_vec_ct = remaining_sample_ct / kQuatersPerVec;
    uint32_t even_ct_incr;
    uint32_t odd_ct_incr;
    uint32_t bothset_ct_incr;
    const uint32_t vec_ct = remaining_full_vec_ct - (remaining_full_vec_ct % 6);
    Count3FreqVec6(R_CAST(const VecUL*, genoarr_iter), vec_ct, &even_ct_incr, &odd_ct_incr, &bothset_ct_incr);
    even_ct += even_ct_incr;
    odd_ct += odd_ct_incr;
    bothset_ct += bothset_ct_incr;
    genoarr_iter = &(genoarr_iter[kWordsPerVec * vec_ct]);
    trail_ct = remaining_sample_ct - (vec_ct * kQuatersPerVec);
  } else {
    genoarr_iter = R_CAST(const uintptr_t*, genoarr);
    trail_ct = sample_ct;
  }
  const uint32_t trail_byte_ct = QuaterCtToByteCt(trail_ct);
  SmallGenoarrCount3FreqIncr(genoarr_iter, trail_byte_ct, &even_ct, &odd_ct, &bothset_ct);
  genocounts[0] = sample_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}

#ifdef __arm__
#  error "Unaligned accesses in GenoarrCountSubsetFreqs()."
#endif
void GenoarrCountSubsetFreqs(const unsigned char* genoarr, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t* genocounts) {
  // does not read past the end of genoarr
  const uintptr_t* genoarr_iter = R_CAST(const uintptr_t*, genoarr);
  const uintptr_t* interleaved_mask_iter = sample_include_interleaved_vec;
  const uint32_t raw_sample_ctv2 = QuaterCtToVecCt(raw_sample_ct);
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
#ifdef USE_AVX2
  const uint32_t halfvec_idx_trail = (raw_sample_ct + 3) / (kBitsPerVec / 4);
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  uintptr_t mask_base3 = 0;
  uintptr_t mask_base4 = 0;
  for (uint32_t vec_idx = 0; vec_idx < raw_sample_ctv2; ++vec_idx) {
    uintptr_t mask_word1;
    uintptr_t mask_word2;
    uintptr_t mask_word3;
    uintptr_t mask_word4;
    if (!(vec_idx % 2)) {
      mask_base1 = *interleaved_mask_iter++;
      mask_base2 = *interleaved_mask_iter++;
      mask_base3 = *interleaved_mask_iter++;
      mask_base4 = *interleaved_mask_iter++;
      mask_word1 = mask_base1 & kMask5555;
      mask_word2 = mask_base2 & kMask5555;
      mask_word3 = mask_base3 & kMask5555;
      mask_word4 = mask_base4 & kMask5555;
    } else {
      mask_word1 = (mask_base1 >> 1) & kMask5555;
      mask_word2 = (mask_base2 >> 1) & kMask5555;
      mask_word3 = (mask_base3 >> 1) & kMask5555;
      mask_word4 = (mask_base4 >> 1) & kMask5555;
    }
    uint32_t uii = 0;
    while (1) {
      uintptr_t cur_geno_word1;
      uintptr_t cur_geno_word2;
      if (2 * vec_idx + uii < halfvec_idx_trail) {
        cur_geno_word1 = *genoarr_iter++;
        cur_geno_word2 = *genoarr_iter++;
      } else {
        const uint32_t remaining_byte_ct = QuaterCtToByteCt(raw_sample_ct) % kBytesPerVec;
        uii = 1;  // todo: check if this harms usual-case loop efficiency
        if (remaining_byte_ct < kBytesPerWord) {
          cur_geno_word1 = ProperSubwordLoad(genoarr_iter, remaining_byte_ct);
          cur_geno_word2 = 0;
        } else {
          cur_geno_word1 = *genoarr_iter++;
          cur_geno_word2 = ProperSubwordLoad(genoarr_iter, remaining_byte_ct - kBytesPerWord);
        }
      }
      const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
      const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
      even_ct += PopcountWord(((cur_geno_word1 & mask_word1) << 1) | (cur_geno_word2 & mask_word2));
      odd_ct += PopcountWord((cur_geno_word1_high_masked << 1) | cur_geno_word2_high_masked);
      bothset_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
      if (uii) {
        break;
      }
      ++uii;
      mask_word1 = mask_word3;
      mask_word2 = mask_word4;
    }
  }
#else  // not USE_AVX2
  const uint32_t vec_idx_trail = (raw_sample_ct + 3) / kQuatersPerVec;
#  ifdef __LP64__
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  for (uint32_t vec_idx = 0; vec_idx < raw_sample_ctv2; ++vec_idx) {
    uintptr_t mask_word1;
    uintptr_t mask_word2;
    if (!(vec_idx % 2)) {
      mask_base1 = *interleaved_mask_iter++;
      mask_base2 = *interleaved_mask_iter++;
      mask_word1 = mask_base1 & kMask5555;
      mask_word2 = mask_base2 & kMask5555;
    } else {
      mask_word1 = (mask_base1 >> 1) & kMask5555;
      mask_word2 = (mask_base2 >> 1) & kMask5555;
    }
    uintptr_t cur_geno_word1;
    uintptr_t cur_geno_word2;
    if (vec_idx < vec_idx_trail) {
      cur_geno_word1 = *genoarr_iter++;
      cur_geno_word2 = *genoarr_iter++;
    } else {
      const uint32_t remaining_byte_ct = QuaterCtToByteCt(raw_sample_ct) % kBytesPerVec;
      if (remaining_byte_ct < kBytesPerWord) {
        cur_geno_word1 = ProperSubwordLoad(genoarr_iter, remaining_byte_ct);
        cur_geno_word2 = 0;
      } else {
        cur_geno_word1 = *genoarr_iter++;
        cur_geno_word2 = ProperSubwordLoad(genoarr_iter, remaining_byte_ct - kBytesPerWord);
      }
    }
    const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
    const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
#    ifdef USE_SSE42
    even_ct += PopcountWord(((cur_geno_word1 & mask_word1) << 1) | (cur_geno_word2 & mask_word2));
    odd_ct += PopcountWord((cur_geno_word1_high_masked << 1) | cur_geno_word2_high_masked);
    bothset_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
#    else
    even_ct += QuatersumWord((cur_geno_word1 & mask_word1) + (cur_geno_word2 & mask_word2));
    odd_ct += QuatersumWord(cur_geno_word1_high_masked + cur_geno_word2_high_masked);
    bothset_ct += QuatersumWord((cur_geno_word1 & cur_geno_word1_high_masked) + (cur_geno_word2 & cur_geno_word2_high_masked));
#    endif
  }
#  else  // not __LP64__
  uintptr_t mask_base = 0;
  for (uint32_t word_idx = 0; word_idx < raw_sample_ctv2; ++word_idx) {
    uintptr_t mask_word;
    if (!(word_idx % 2)) {
      mask_base = *interleaved_mask_iter++;
      mask_word = mask_base & kMask5555;
    } else {
      mask_word = (mask_base >> 1) & kMask5555;
    }
    uintptr_t cur_geno_word;
    if (word_idx < vec_idx_trail) {
      cur_geno_word = *genoarr_iter++;
    } else {
      const uint32_t remaining_byte_ct = QuaterCtToByteCt(raw_sample_ct) % kBytesPerVec;
      cur_geno_word = ProperSubwordLoad(genoarr_iter, remaining_byte_ct);
    }
    const uintptr_t cur_geno_word_high_masked = mask_word & (cur_geno_word >> 1);
    even_ct += Popcount01Word(cur_geno_word & mask_word);
    odd_ct += Popcount01Word(cur_geno_word_high_masked);
    bothset_ct += Popcount01Word(cur_geno_word & cur_geno_word_high_masked);
  }
#  endif  // not __LP64__
#endif  // not USE_AVX2
  genocounts[0] = sample_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}

void GenoarrCountSubsetFreqs2(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t* genocounts) {
  // slower GenoarrCountSubsetFreqs() which does not require
  // sample_include_interleaved_vec to be precomputed.
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctl2 = QuaterCtToWordCt(raw_sample_ct);
  const uint32_t fullword_ct = raw_sample_ctl2 / 2;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  for (uint32_t widx = 0; widx < fullword_ct; ++widx) {
    const uintptr_t mask_word = sample_include[widx];
    if (mask_word) {
      uintptr_t geno_word = genoarr[2 * widx];
      uintptr_t geno_even = PackWordToHalfword(geno_word & kMask5555);
      uintptr_t geno_odd = PackWordToHalfword((geno_word >> 1) & kMask5555);
      geno_word = genoarr[2 * widx + 1];
      geno_even |= S_CAST(uintptr_t, PackWordToHalfword(geno_word & kMask5555)) << kBitsPerWordD2;
      geno_odd |= S_CAST(uintptr_t, PackWordToHalfword((geno_word >> 1) & kMask5555)) << kBitsPerWordD2;
      const uintptr_t geno_even_masked = geno_even & mask_word;
      even_ct += PopcountWord(geno_even_masked);
      odd_ct += PopcountWord(geno_odd & mask_word);
      bothset_ct += PopcountWord(geno_odd & geno_even_masked);
    }
  }
  if (raw_sample_ctl2 % 2) {
    const uintptr_t mask_hw = sample_include[fullword_ct];
    if (mask_hw) {
      const uintptr_t geno_word = genoarr[2 * fullword_ct];
      // todo: benchmark main loop unpack vs. pack
      const uintptr_t mask_word = UnpackHalfwordToWord(mask_hw);
      const uintptr_t geno_word_shifted = geno_word >> 1;
      const uintptr_t geno_word_masked = geno_word & mask_word;
      even_ct += Popcount01Word(geno_word_masked);
      odd_ct += Popcount01Word(geno_word_shifted & mask_word);
      bothset_ct += Popcount01Word(geno_word_masked & geno_word_shifted);
    }
  }
  genocounts[0] = sample_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}

void GenoarrCountSubsetIntersectFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict subset1, const uintptr_t* __restrict subset2, uint32_t raw_sample_ct, uint32_t* genocounts) {
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctl2 = QuaterCtToWordCt(raw_sample_ct);
  const uint32_t fullword_ct = raw_sample_ctl2 / 2;
  uint32_t subset_intersect_ct = 0;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  for (uint32_t widx = 0; widx < fullword_ct; ++widx) {
    const uintptr_t mask_word = subset1[widx] & subset2[widx];
    if (mask_word) {
      uintptr_t geno_word = genoarr[2 * widx];
      uintptr_t geno_even = PackWordToHalfword(geno_word & kMask5555);
      uintptr_t geno_odd = PackWordToHalfword((geno_word >> 1) & kMask5555);
      geno_word = genoarr[2 * widx + 1];
      geno_even |= S_CAST(uintptr_t, PackWordToHalfword(geno_word & kMask5555)) << kBitsPerWordD2;
      geno_odd |= S_CAST(uintptr_t, PackWordToHalfword((geno_word >> 1) & kMask5555)) << kBitsPerWordD2;
      const uintptr_t geno_even_masked = geno_even & mask_word;
      subset_intersect_ct += PopcountWord(mask_word);
      even_ct += PopcountWord(geno_even_masked);
      odd_ct += PopcountWord(geno_odd & mask_word);
      bothset_ct += PopcountWord(geno_odd & geno_even_masked);
    }
  }
  if (raw_sample_ctl2 % 2) {
    const uintptr_t mask_hw = subset1[fullword_ct] & subset2[fullword_ct];
    if (mask_hw) {
      const uintptr_t geno_word = genoarr[fullword_ct * 2];
      const uintptr_t mask_word = UnpackHalfwordToWord(mask_hw);
      const uintptr_t geno_word_shifted = geno_word >> 1;
      const uintptr_t geno_word_masked = geno_word & mask_word;
      subset_intersect_ct += Popcount01Word(mask_word);
      even_ct += Popcount01Word(geno_word_masked);
      odd_ct += Popcount01Word(geno_word_shifted & mask_word);
      bothset_ct += Popcount01Word(geno_word_masked & geno_word_shifted);
    }
  }
  genocounts[0] = subset_intersect_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}

void GenovecCountFreqs(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* genocounts) {
  // ok to read trailing genovec bytes, but must mask them out
  const uint32_t sample_ct_remainder = sample_ct % kBitsPerWordD2;
  GenovecCountFreqsUnsafe(genovec, sample_ct - sample_ct_remainder, genocounts);
  if (sample_ct_remainder) {
    uintptr_t cur_geno_word = bzhi(genovec[sample_ct / kBitsPerWordD2], 2 * sample_ct_remainder);
    const uintptr_t cur_geno_word_high = kMask5555 & (cur_geno_word >> 1);
    const uint32_t even_ct = Popcount01Word(cur_geno_word & kMask5555);
    const uint32_t odd_ct = Popcount01Word(cur_geno_word_high);
    const uint32_t bothset_ct = Popcount01Word(cur_geno_word & cur_geno_word_high);
    genocounts[0] += sample_ct_remainder + bothset_ct - even_ct - odd_ct;
    genocounts[1] += even_ct - bothset_ct;
    genocounts[2] += odd_ct - bothset_ct;
    genocounts[3] += bothset_ct;
  }
}

void GenovecInvertUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // flips 0 to 2 and vice versa.
  // "unsafe" because trailing bits are not zeroed out.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(IsVecAligned(genovec));
  const VecUL not_m1 = VCONST_UL(kMaskAAAA);
  VecUL* vptr = R_CAST(VecUL*, genovec);
  for (uint32_t vidx = 0; vidx < vec_ct; ++vidx) {
    VecUL cur_vec = vptr[vidx];
    // flip high bit iff low bit is unset
    vptr[vidx] = cur_vec ^ ((~vecul_slli(cur_vec, 1)) & not_m1);
  }
}

void GenovecInvertCopyUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict genovec_inverted_copy) {
  // flips 0 to 2 and vice versa.
  // "unsafe" because trailing bits are not zeroed out.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(IsVecAligned(genovec));
  const VecUL not_m1 = VCONST_UL(kMaskAAAA);
  const VecUL* vin_ptr = R_CAST(const VecUL*, genovec);
  VecUL* vout_ptr = R_CAST(VecUL*, genovec_inverted_copy);
  for (uint32_t vidx = 0; vidx < vec_ct; ++vidx) {
    VecUL cur_vec = vin_ptr[vidx];
    // flip high bit iff low bit is unset
    vout_ptr[vidx] = cur_vec ^ ((~vecul_slli(cur_vec, 1)) & not_m1);
  }
}

void GenovecNonmissingToZeroUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // sets 1 and 2 to zero; leaves 3s untouched.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(IsVecAligned(genovec));
  const VecUL m1 = VCONST_UL(kMask5555);
  VecUL* vptr = R_CAST(VecUL*, genovec);
  for (uint32_t vidx = 0; vidx < vec_ct; ++vidx) {
    VecUL cur_vec = vptr[vidx];
    VecUL cur_vec_rshifted = vecul_srli(cur_vec, 1);
    cur_vec = cur_vec & m1;
    cur_vec = cur_vec & cur_vec_rshifted;
    vptr[vidx] = cur_vec | vecul_slli(cur_vec, 1);
  }
}

void GenovecNonzeroToMissingUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // converts 1s and 2s to 3s, leaves zeroes untouched.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(IsVecAligned(genovec));
  const VecUL m1 = VCONST_UL(kMask5555);
  VecUL* vptr = R_CAST(VecUL*, genovec);
  for (uint32_t vidx = 0; vidx < vec_ct; ++vidx) {
    VecUL cur_vec = vptr[vidx];
    VecUL cur_vec_rshifted = vecul_srli(cur_vec, 1);
    cur_vec = cur_vec | cur_vec_rshifted;
    cur_vec = cur_vec & m1;
    vptr[vidx] = cur_vec | vecul_slli(cur_vec, 1);
  }
}

void DifflistCountSubsetFreqs(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t common_geno, uint32_t difflist_len, uint32_t sample_ct, uint32_t* genocounts) {
  ZeroUiArr(4, genocounts);
  uint32_t common_geno_ct = sample_ct;
  for (uint32_t difflist_idx = 0; difflist_idx < difflist_len; ++difflist_idx) {
    const uint32_t raw_sample_idx = difflist_sample_ids[difflist_idx];
    if (IsSet(sample_include, raw_sample_idx)) {
      genocounts[GetQuaterarrEntry(raregeno, difflist_idx)] += 1;
      --common_geno_ct;
    }
  }
  genocounts[common_geno] = common_geno_ct;
}


static_assert(kPglQuaterTransposeBatch == S_CAST(uint32_t, kQuatersPerCacheline), "TransposeQuaterblock() needs to be updated.");
#ifdef USE_AVX2
void TransposeQuaterblockAvx2(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, void* __restrict buf0) {
  // buf must be vector-aligned and have size 16k
  const uint32_t block_ct = QuaterCtToByteCt(write_batch_size);
  // fold the first 6 shuffles into the initial ingestion loop
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t write_ul_ct = QuaterCtToWordCt(read_batch_size);
  const uint32_t read_batch_rem = write_ul_ct * 32 - read_batch_size;
  uintptr_t* target_iter0 = write_iter;
  const uint32_t full_block_ct = write_batch_size / 4;
  // Tried _pext_u64() too, that was slightly worse.
  // Tried _mm256_i32gather_epi32(), didn't help.
  // But seems to be fine to just process 4 write-rows at a time.
  for (uint32_t block_idx = 0; block_idx < block_ct; ++block_idx) {
    const unsigned char* read_iter_tmp = R_CAST(const unsigned char*, read_iter);
    read_iter_tmp = &(read_iter_tmp[block_idx]);
    unsigned char* initial_target_iter = S_CAST(unsigned char*, buf0);
    for (uint32_t ujj = 0; ujj < read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_byte_stride]);
    }
    memset(initial_target_iter, 0, read_batch_rem);
    const __m256i* source_iter = S_CAST(__m256i*, buf0);
    if (block_idx == full_block_ct) {
      const uint32_t remainder = write_batch_size % 4;
      const uint32_t initial_lshift = 2 * (4 - remainder);
      for (uint32_t widx = 0; widx < write_ul_ct; ++widx) {
        __m256i loader = *source_iter++;
        loader = _mm256_slli_epi64(loader, initial_lshift);
        uint32_t write_idx_lowbits = remainder - 1;
        while (1) {
          const uintptr_t uint_hi = _mm256_movemask_epi8(loader);
          loader = _mm256_slli_epi64(loader, 1);
          const uintptr_t uint_lo = _mm256_movemask_epi8(loader);
          target_iter0[write_ul_stride * write_idx_lowbits] = _pdep_u64(uint_hi, kMaskAAAA) | _pdep_u64(uint_lo, kMask5555);
          if (!write_idx_lowbits) {
            break;
          }
          loader = _mm256_slli_epi64(loader, 1);
          --write_idx_lowbits;
        }
        ++target_iter0;
      }
      break;
    }

    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    uintptr_t* target_iter2 = &(target_iter1[write_ul_stride]);
    uintptr_t* target_iter3 = &(target_iter2[write_ul_stride]);
    for (uint32_t widx = 0; widx < write_ul_ct; ++widx) {
      __m256i loader = source_iter[widx];
      uintptr_t uint_hi = _mm256_movemask_epi8(loader);
      loader = _mm256_slli_epi64(loader, 1);
      uintptr_t uint_lo = _mm256_movemask_epi8(loader);
      target_iter3[widx] = _pdep_u64(uint_hi, kMaskAAAA) | _pdep_u64(uint_lo, kMask5555);
      loader = _mm256_slli_epi64(loader, 1);
      uint_hi = _mm256_movemask_epi8(loader);
      loader = _mm256_slli_epi64(loader, 1);
      uint_lo = _mm256_movemask_epi8(loader);
      target_iter2[widx] = _pdep_u64(uint_hi, kMaskAAAA) | _pdep_u64(uint_lo, kMask5555);
      loader = _mm256_slli_epi64(loader, 1);
      uint_hi = _mm256_movemask_epi8(loader);
      loader = _mm256_slli_epi64(loader, 1);
      uint_lo = _mm256_movemask_epi8(loader);
      target_iter1[widx] = _pdep_u64(uint_hi, kMaskAAAA) | _pdep_u64(uint_lo, kMask5555);
      loader = _mm256_slli_epi64(loader, 1);
      uint_hi = _mm256_movemask_epi8(loader);
      loader = _mm256_slli_epi64(loader, 1);
      uint_lo = _mm256_movemask_epi8(loader);
      target_iter0[widx] = _pdep_u64(uint_hi, kMaskAAAA) | _pdep_u64(uint_lo, kMask5555);
    }
    target_iter0 = &(target_iter3[write_ul_stride]);
  }
}
#else  // !USE_AVX2
#  ifdef __LP64__
static_assert(kWordsPerVec == 2, "TransposeQuaterblock() needs to be updated.");
#  else
static_assert(kWordsPerVec == 1, "TransposeQuaterblock() needs to be updated.");
#  endif
void TransposeQuaterblockNoAvx2(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1) {
  // buf must be vector-aligned and have size 32k
  // defining them as unsigned char* might prevent a strict-aliasing issue?
  // (might need to go through greater contortions to actually be safe?)
  const uint32_t initial_read_byte_ct = QuaterCtToByteCt(write_batch_size);
  // fold the first 6 shuffles into the initial ingestion loop
  const unsigned char* initial_read_iter = R_CAST(const unsigned char*, read_iter);
  const unsigned char* initial_read_end = &(initial_read_iter[initial_read_byte_ct]);
  unsigned char* initial_target_iter = buf0;
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t read_batch_rem = kQuatersPerCacheline - read_batch_size;
  for (; initial_read_iter < initial_read_end; ++initial_read_iter) {
    const unsigned char* read_iter_tmp = initial_read_iter;
    for (uint32_t ujj = 0; ujj < read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_byte_stride]);
    }
    initial_target_iter = R_CAST(unsigned char*, memseta(initial_target_iter, 0, read_batch_rem));
  }

  // second-to-last shuffle, 8 bit spacing -> 4
  const VecUL* source_iter = R_CAST(VecUL*, buf0);
  uintptr_t* target_iter0 = R_CAST(uintptr_t*, buf1);
#  ifdef __LP64__
  const VecUL m4 = VCONST_UL(kMask0F0F);
  const VecUL m8 = VCONST_UL(kMask00FF);
  const VecUL m16 = VCONST_UL(kMask0000FFFF);
#  endif
  const uint32_t write_word_ct = QuaterCtToWordCt(read_batch_size);
  const uint32_t penult_inner_loop_iter_ct = 2 * write_word_ct;
  const uint32_t cur_write_skip = 2 * kWordsPerCacheline - penult_inner_loop_iter_ct;
  // coincidentally, this also needs to run DivUp(write_batch_size, 4) times
  for (uint32_t uii = 0; uii < initial_read_byte_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 2]);
    for (uint32_t ujj = 0; ujj < penult_inner_loop_iter_ct; ++ujj) {
#  ifdef __LP64__
      const VecUL loader = *source_iter++;
      VecUL target0 = loader & m4;
      VecUL target1 = (vecul_srli(loader, 4)) & m4;
      target0 = (target0 | (vecul_srli(target0, 4))) & m8;
      target1 = (target1 | (vecul_srli(target1, 4))) & m8;
      target0 = (target0 | (vecul_srli(target0, 8))) & m16;
      target1 = (target1 | (vecul_srli(target1, 8))) & m16;
      UniVec target0u;
      UniVec target1u;
      target0u.vi = target0 | (vecul_srli(target0, 16));
      target1u.vi = target1 | (vecul_srli(target1, 16));
      *target_iter0++ = S_CAST(uint32_t, target0u.u8[0]) | (target0u.u8[1] << 32);
      *target_iter1++ = S_CAST(uint32_t, target1u.u8[0]) | (target1u.u8[1] << 32);
#  else
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
#  endif
    }
#  ifdef __LP64__
    source_iter = &(source_iter[cur_write_skip]);
#  else
    source_iter = &(source_iter[2 * cur_write_skip]);
#  endif
    target_iter0 = &(target_iter1[cur_write_skip]);
  }

  // last shuffle, 4 bit spacing -> 2
  source_iter = R_CAST(VecUL*, buf1);
  target_iter0 = write_iter;
#  ifdef __LP64__
  const VecUL m2 = VCONST_UL(kMask3333);
#  endif
  const uint32_t last_loop_iter_ct = DivUp(write_batch_size, 2);
  for (uint32_t uii = 0; uii < last_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    for (uint32_t ujj = 0; ujj < write_word_ct; ++ujj) {
#  ifdef __LP64__
      // in AVX2 case, use write_dword_ct instead of write_word_ct, etc.
      const VecUL loader = *source_iter++;
      VecUL target0 = loader & m2;
      VecUL target1 = (vecul_srli(loader, 2)) & m2;
      target0 = (target0 | (vecul_srli(target0, 2))) & m4;
      target1 = (target1 | (vecul_srli(target1, 2))) & m4;
      target0 = (target0 | (vecul_srli(target0, 4))) & m8;
      target1 = (target1 | (vecul_srli(target1, 4))) & m8;
      target0 = (target0 | (vecul_srli(target0, 8))) & m16;
      target1 = (target1 | (vecul_srli(target1, 8))) & m16;
      UniVec target0u;
      UniVec target1u;
      target0u.vi = target0 | (vecul_srli(target0, 16));
      target1u.vi = target1 | (vecul_srli(target1, 16));
      target_iter0[ujj] = S_CAST(uint32_t, target0u.u8[0]) | (target0u.u8[1] << 32);
      target_iter1[ujj] = S_CAST(uint32_t, target1u.u8[0]) | (target1u.u8[1] << 32);
#  else
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
      target_iter0[ujj] = S_CAST(Halfword, target_word0_lo) | (target_word0_hi << kBitsPerWordD2);
      target_iter1[ujj] = S_CAST(Halfword, target_word1_lo) | (target_word1_hi << kBitsPerWordD2);
#  endif
    }
#  ifdef __LP64__
    source_iter = &(source_iter[kWordsPerCacheline - write_word_ct]);
#  else
    source_iter = &(source_iter[2 * (kWordsPerCacheline - write_word_ct)]);
#  endif
    target_iter0 = &(target_iter1[write_ul_stride]);
  }
}
#endif  // !USE_AVX2

void BiallelicDosage16Invert(uint32_t dosage_ct, uint16_t* dosage_vals) {
  // replace each x with (32768 - x).
  // uses vector operations, but does not require dosage_vals to be
  // vec-aligned.
  const VecUL subvec = VCONST_UL(32768 * kMask0001);
  const uint32_t lead_usi_ct = ((-R_CAST(uintptr_t, dosage_vals)) % kBytesPerVec) / sizeof(int16_t);
  if (dosage_ct >= lead_usi_ct) {
    for (uint32_t uii = 0; uii < lead_usi_ct; ++uii) {
      *dosage_vals = 32768 - (*dosage_vals);
      ++dosage_vals;
    }
    dosage_ct -= lead_usi_ct;
    const uint32_t vec_ct = dosage_ct / (kBytesPerVec / sizeof(int16_t));
    VecUL* dosage_vals_iter = R_CAST(VecUL*, dosage_vals);
    for (uint32_t vec_idx = 0; vec_idx < vec_ct; ++vec_idx) {
      const VecUL cur_vec = *dosage_vals_iter;
      *dosage_vals_iter++ = subvec - cur_vec;
    }
    dosage_vals = &(dosage_vals[vec_ct * (kBytesPerVec / sizeof(int16_t))]);
    dosage_ct -= vec_ct * (kBytesPerVec / sizeof(int16_t));
  }
  for (uint32_t uii = 0; uii < dosage_ct; ++uii) {
    dosage_vals[uii] = 32768 - dosage_vals[uii];
  }
}

void GenovecToMissingnessUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict missingness) {
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  Halfword* missingness_alias = R_CAST(Halfword*, missingness);
  for (uint32_t widx = 0; widx < sample_ctl2; ++widx) {
    const uintptr_t cur_geno_word = genovec[widx];
    missingness_alias[widx] = PackWordToHalfword(cur_geno_word & (cur_geno_word >> 1) & kMask5555);
  }
  if (sample_ctl2 % 2) {
    missingness_alias[sample_ctl2] = 0;
  }
}

void GenovecToNonmissingnessUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict nonmissingness) {
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  Halfword* nonmissingness_alias = R_CAST(Halfword*, nonmissingness);
  for (uint32_t widx = 0; widx < sample_ctl2; ++widx) {
    const uintptr_t cur_geno_word = genovec[widx];
    nonmissingness_alias[widx] = PackWordToHalfword((~(cur_geno_word & (cur_geno_word >> 1))) & kMask5555);
  }
}

void SplitHomRef2hetUnsafeW(const uintptr_t* genoarr, uint32_t inword_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf) {
  Halfword* hom_alias = R_CAST(Halfword*, hom_buf);
  Halfword* ref2het_alias = R_CAST(Halfword*, ref2het_buf);
  for (uint32_t widx = 0; widx < inword_ct; ++widx) {
    const uintptr_t geno_word = genoarr[widx];
    hom_alias[widx] = PackWordToHalfword(kMask5555 & (~geno_word));
    ref2het_alias[widx] = PackWordToHalfword(kMask5555 & (~(geno_word >> 1)));
  }
}

void SplitHomRef2het(const uintptr_t* genoarr, uint32_t sample_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf) {
  const uint32_t full_outword_ct = sample_ct / kBitsPerWord;
  SplitHomRef2hetUnsafeW(genoarr, full_outword_ct * 2, hom_buf, ref2het_buf);
  const uint32_t remainder = sample_ct % kBitsPerWord;
  if (remainder) {
    uintptr_t geno_word = genoarr[full_outword_ct * 2];
    uintptr_t hom_word = PackWordToHalfword(kMask5555 & (~geno_word));
    uintptr_t ref2het_word = PackWordToHalfword(kMask5555 & (~(geno_word >> 1)));
    if (remainder > kBitsPerWordD2) {
      geno_word = genoarr[full_outword_ct * 2 + 1];
      hom_word |= S_CAST(uintptr_t, PackWordToHalfword(kMask5555 & (~geno_word))) << kBitsPerWordD2;
      ref2het_word |= S_CAST(uintptr_t, PackWordToHalfword(kMask5555 & (~(geno_word >> 1)))) << kBitsPerWordD2;
    }
    const uintptr_t cur_mask = (k1LU << remainder) - 1;
    hom_buf[full_outword_ct] = hom_word & cur_mask;
    ref2het_buf[full_outword_ct] = ref2het_word & cur_mask;
  }
}


void PreinitPgfi(PgenFileInfo* pgfip) {
  pgfip->shared_ff = nullptr;
  pgfip->block_base = nullptr;
  // we want this for proper handling of e.g. sites-only VCFs
  pgfip->nonref_flags = nullptr;
}

uint32_t CountPgfiAllocCachelinesRequired(uint32_t raw_variant_ct) {
  // assumes variable-width variant records, otherwise pgfi.vrtypes and
  // pgfi.vr_fpos can just be nullptr.

  // vrtypes: 1 byte per entry, (raw_variant_ct + 1) entries
  uint32_t cachelines_required = 1 + (raw_variant_ct / kCacheline);

  // var_fpos: 8 bytes per entry, (raw_variant_ct + 1) entries
  cachelines_required += 1 + (raw_variant_ct / kInt64PerCacheline);
  return cachelines_required;
}

uint32_t CountPgrAllocCachelinesRequired(uint32_t raw_sample_ct, PgenGlobalFlags gflags, uint32_t max_alt_allele_ct, uint32_t fread_buf_byte_ct) {
  // workspace_vec: always needed, 2 bits per entry, up to raw_sample_ct
  // entries
  const uint32_t genovec_cacheline_req = QuaterCtToCachelineCt(raw_sample_ct);
  const uint32_t bitvec_cacheline_req = BitCtToCachelineCt(raw_sample_ct);
  uint32_t cachelines_required = genovec_cacheline_req;
  // fread_buf.  fread_buf_byte_ct should be zero if mmap() is being used.
  // DivUp() won't overflow since fread_buf_byte_ct requirement can't exceed
  // kPglMaxBytesPerVariant, which is sufficiently far from 2^32.
  cachelines_required += DivUp(fread_buf_byte_ct, kCacheline);

  const uint32_t ld_compression_present = (gflags / kfPgenGlobalLdCompressionPresent) & 1;
  if (gflags & kfPgenGlobalDifflistOrLdPresent) {
    const uint32_t max_difflist_entry_ct_base = (raw_sample_ct / kPglMaxDifflistLenDivisor);
    // const uint32_t max_difflist_entry_ct = max_difflist_entry_ct_base * (1 + ld_compression_present);
    // workspace_raregeno_vec
    cachelines_required += QuaterCtToCachelineCt(max_difflist_entry_ct_base);

    // workspace_difflist_sample_ids
    // bugfix: must add 1 since several routines add a terminator element
    cachelines_required += 1 + (max_difflist_entry_ct_base / kInt32PerCacheline);

    // workspace_raregeno_tmp_loadbuf
    cachelines_required += QuaterCtToCachelineCt(max_difflist_entry_ct_base);

    // workspace_difflist_sample_ids_tmp
    cachelines_required += 1 + (max_difflist_entry_ct_base / kInt32PerCacheline);

    if (ld_compression_present) {
      // ldbase_genovec
      cachelines_required += genovec_cacheline_req;

      // ldbase_raregeno
      cachelines_required += QuaterCtToCachelineCt(max_difflist_entry_ct_base);

      // ldbase_difflist_sample_ids
      cachelines_required += 1 + (max_difflist_entry_ct_base / kInt32PerCacheline);
    }
  }
  if (max_alt_allele_ct > 1) {
    // workspace_aux1_nonmissing_vec
    cachelines_required += bitvec_cacheline_req;

    // workspace_aux1_code_vec
    // prepare for worst-case scenario for now.  todo: use loaded variant
    // record lengths to bound this when appropriate
    uintptr_t aux1_allele_bytect = GetAux1AlleleByteCt(max_alt_allele_ct, raw_sample_ct);
    if (aux1_allele_bytect > kPglMaxBytesPerVariant) {
      // but assume the file isn't actually invalid
      aux1_allele_bytect = kPglMaxBytesPerVariant;
    }

    cachelines_required += DivUp(aux1_allele_bytect, kCacheline);

    // workspace_ambig_sample_ids
    cachelines_required += Int32CtToCachelineCt(raw_sample_ct);
  }
  if (gflags & kfPgenGlobalHardcallPhasePresent) {
    // workspace_all_hets, possibly ldbase_all_hets
    cachelines_required += bitvec_cacheline_req * (1 + ld_compression_present);
  }
  if (gflags & kfPgenGlobalDosagePresent) {
    // aux track #3: usually bitarray tracking which samples have dosage info
    // (may be stored on disk as a dosage list)
    cachelines_required += bitvec_cacheline_req;
    if (gflags & kfPgenGlobalDosagePhasePresent) {
      // aux track #5: bitarray tracking which dosage entries are phased
      cachelines_required += bitvec_cacheline_req;

      // phased aux tracks #4,6: max_alt_allele_ct * 4 bytes per sample
      // There may be overflow risk here in the future.
      // (commented out since caller always provides this buffer for now)
      // cachelines_required += DivUp(max_alt_allele_ct * 4 * k1LU * raw_sample_ct, kCacheline);
    }
    // unphased aux track #4: max_alt_allele_ct * 2 bytes per sample
    // cachelines_required += DivUp(max_alt_allele_ct * 2 * k1LU * raw_sample_ct, kCacheline);
  }
  return cachelines_required;
}

static_assert(kPglMaxAltAlleleCt == 254, "Need to update PgfiInitPhase1().");
PglErr PgfiInitPhase1(const char* fname, uint32_t raw_variant_ct, uint32_t raw_sample_ct, uint32_t use_mmap, PgenHeaderCtrl* header_ctrl_ptr, PgenFileInfo* pgfip, uintptr_t* pgfi_alloc_cacheline_ct_ptr, char* errstr_buf) {
  pgfip->var_fpos = nullptr;
  pgfip->vrtypes = nullptr;
  pgfip->allele_idx_offsets = nullptr;
  pgfip->nonref_flags = nullptr;

  pgfip->max_alt_allele_ct = 1;
  pgfip->max_dosage_alt_allele_ct = 0;

  pgfip->block_base = nullptr;
  // this should force overflow when value is uninitialized.
  pgfip->block_offset = 1LLU << 63;

  uint64_t fsize;
  const unsigned char* fread_ptr;
  FILE* shared_ff = nullptr;
  unsigned char small_readbuf[3];
  if (use_mmap) {
    pgfip->shared_ff = nullptr;  // this must be initialized before block_base
#ifdef NO_MMAP
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase1() use_mmap parameter is nonzero, but pgenlib was not compiled with mmap support.\n");
    return kPglRetImproperFunctionCall;
#else
    int32_t file_handle = open(fname, O_RDONLY);
    if (file_handle < 0) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Failed to open %s.\n", fname);
      return kPglRetOpenFail;
    }
    struct stat statbuf;
    if (fstat(file_handle, &statbuf) < 0) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Failed to open %s.\n", fname);
      return kPglRetOpenFail;
    }
    fsize = statbuf.st_size;
    pgfip->block_offset = 0;
    pgfip->file_size = fsize;
    pgfip->block_base = S_CAST(const unsigned char*, mmap(0, pgfip->file_size, PROT_READ, MAP_SHARED, file_handle, 0));
    if (R_CAST(uintptr_t, pgfip->block_base) == (~k0LU)) {
      pgfip->block_base = nullptr;
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
    // this provided less than a ~5% boost on OS X; mmap still took >80% longer
    // than fread on an 85GB file there
    // try MAP_POPULATE on Linux?
    // madvise((unsigned char*)(pgfip->block_base), fsize, MADV_SEQUENTIAL);
    close(file_handle);
    // update (7 Jan 2018): drop support for zero-sample and zero-variant
    // files, not worth the development cost
    if (fsize < 4) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s is too small to be a valid .pgen file.\n", fname);
      return kPglRetMalformedInput;
    }
    fread_ptr = pgfip->block_base;
#endif
  } else {
    shared_ff = fopen(fname, FOPEN_RB);
    pgfip->shared_ff = shared_ff;
    if (!shared_ff) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Failed to open %s.\n", fname);
      return kPglRetOpenFail;
    }
    if (fseeko(shared_ff, 0, SEEK_END)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
    fsize = ftello(shared_ff);
    if (fsize < 4) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s is too small to be a valid .pgen file.\n", fname);
      return kPglRetMalformedInput;
    }
    rewind(shared_ff);
    if (!fread_unlocked(small_readbuf, 3, 1, shared_ff)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
    fread_ptr = small_readbuf;
  }
  // deliberate underflow
  if (((raw_variant_ct - 1) > 0x7ffffffc) && (raw_variant_ct != UINT32_MAX)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid raw_variant_ct function parameter.\n");
    return kPglRetImproperFunctionCall;
  }
  if (((raw_sample_ct - 1) > 0x7ffffffd) && (raw_sample_ct != UINT32_MAX)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid raw_sample_ct function parameter.\n");
    return kPglRetImproperFunctionCall;
  }
  if (memcmp(fread_ptr, "l\x1b", 2)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s is not a .pgen file (first two bytes don't match the magic number).\n", fname);
    return kPglRetMalformedInput;
  }
  const uint32_t file_type_code = fread_ptr[2];
  *header_ctrl_ptr = 0;
  if (file_type_code < 2) {
    // plink 1 binary
    if (!file_type_code) {
      // sample-major.  validate file size here so we don't have to recheck it
      if ((raw_sample_ct != UINT32_MAX) && (raw_variant_ct != UINT32_MAX)) {
        const uint64_t fsize_expected = 3 + S_CAST(uint64_t, raw_sample_ct) * QuaterCtToByteCt(raw_variant_ct);
        if (fsize != fsize_expected) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Unexpected PLINK 1 sample-major .bed file size (%" PRIu64 " bytes expected).\n", fsize_expected);
          return kPglRetMalformedInput;
        }
      }
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: pgenlib does not directly support sample-major PLINK 1 .bed files.\n(However, PLINK 2 automatically transposes and compresses them for you.)\n");
      return kPglRetSampleMajorBed;
    }
    if (raw_sample_ct == UINT32_MAX) {
      // either .fam must be loaded first, or user must provide sample count
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase1() must be called with an accurate raw_sample_ct value, since %s is a PLINK 1 .bed file.\n", fname);
      return kPglRetImproperFunctionCall;
    }
    const uint32_t const_vrec_width = QuaterCtToByteCt(raw_sample_ct);
    if (raw_variant_ct == UINT32_MAX) {
      // allow raw_variant_ct to be inferred
      uint64_t quotient = (fsize - 3) / const_vrec_width;
      if ((quotient > 0x7fffffffU) || (quotient * const_vrec_width + 3 != fsize)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Unexpected PLINK 1 .bed file size (since raw_sample_ct was %u, [file size - 3] should be divisible by %u and the quotient should be smaller than 2^31).\n", raw_sample_ct, const_vrec_width);
        return kPglRetMalformedInput;
      }
      raw_variant_ct = quotient;
    } else {
      if (S_CAST(uint64_t, raw_variant_ct) * const_vrec_width + 3 != fsize) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Unexpected PLINK 1 .bed file size (expected %" PRIu64 " bytes).\n", S_CAST(uint64_t, raw_variant_ct) * const_vrec_width + 3);
        return kPglRetMalformedInput;
      }
    }
    pgfip->raw_variant_ct = raw_variant_ct;
    pgfip->raw_sample_ct = raw_sample_ct;
    pgfip->const_fpos_offset = 3;

    pgfip->const_vrtype = kPglVrtypePlink1;
    pgfip->const_vrec_width = const_vrec_width;
    pgfip->gflags = kfPgenGlobalAllNonref;
    *pgfi_alloc_cacheline_ct_ptr = 0;
    return kPglRetSuccess;
  }

  if (fsize < 12) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s is too small to be a valid .pgen file.\n", fname);
    return kPglRetMalformedInput;
  }
#ifndef NO_MMAP
  if (use_mmap) {
    memcpy(&(pgfip->raw_variant_ct), &(fread_ptr[3]), sizeof(int32_t));
    memcpy(&(pgfip->raw_sample_ct), &(fread_ptr[7]), sizeof(int32_t));
    memcpy(header_ctrl_ptr, &(fread_ptr[11]), 1);
  } else {
#endif
    if ((!fread_unlocked(&(pgfip->raw_variant_ct), sizeof(int32_t), 1, shared_ff)) ||
        (!fread_unlocked(&(pgfip->raw_sample_ct), sizeof(int32_t), 1, shared_ff)) ||
        (!fread_unlocked(header_ctrl_ptr, 1, 1, shared_ff))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
#ifndef NO_MMAP
  }
#endif
  PgenHeaderCtrl header_ctrl = *header_ctrl_ptr;
  if (raw_variant_ct == UINT32_MAX) {
    raw_variant_ct = pgfip->raw_variant_ct;
    // deliberate underflow
    if ((raw_variant_ct - 1) > 0x7ffffffc) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid variant count in .pgen file.\n");
      return kPglRetMalformedInput;
    }
  } else if (raw_variant_ct != pgfip->raw_variant_ct) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase1() was called with raw_variant_ct == %u, but %s contains %u variant%s.\n", raw_variant_ct, fname, pgfip->raw_variant_ct, (pgfip->raw_variant_ct == 1)? "" : "s");
    return kPglRetInconsistentInput;
  }
  if (raw_sample_ct == UINT32_MAX) {
    raw_sample_ct = pgfip->raw_sample_ct;
    // deliberate underflow
    if ((raw_sample_ct - 1) > 0x7ffffffd) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid sample count in .pgen file.\n");
      return kPglRetMalformedInput;
    }
  } else if (raw_sample_ct != pgfip->raw_sample_ct) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase1() was called with raw_sample_ct == %u, but %s contains %u sample%s.\n", raw_sample_ct, fname, pgfip->raw_sample_ct, (pgfip->raw_sample_ct == 1)? "" : "s");
    return kPglRetInconsistentInput;
  }
  pgfip->gflags = kfPgenGlobal0;
  pgfip->const_fpos_offset = 12;

  // explicit storage of "is this reference allele untrusted?"
  // need caller to allocate this
  uint32_t nonref_flags_storage = header_ctrl >> 6;
  if (nonref_flags_storage == 3) {
    pgfip->const_fpos_offset += DivUp(raw_variant_ct, CHAR_BIT);
  } else if (nonref_flags_storage == 2) {
    pgfip->gflags |= kfPgenGlobalAllNonref;
  }

  if (file_type_code < 16) {
    // plink 2 binary, single constant-width vrtype
    if (file_type_code > 4) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Third byte of %s does not correspond to a storage mode supported by this version of pgenlib.\n", fname);
      return kPglRetNotYetSupported;
    }
    if (header_ctrl & 63) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Third byte of %s corresponds to a fixed-width storage mode, but twelfth byte is only consistent with a variable-width mode.\n", fname);
      return kPglRetMalformedInput;
    }
    uint32_t vrtype = 0;
    uintptr_t const_vrec_width = QuaterCtToByteCt(raw_sample_ct);
    if (file_type_code == 3) {
      vrtype = 0x40;
      const_vrec_width += raw_sample_ct * 2;
      pgfip->gflags |= kfPgenGlobalDosagePresent;
    } else if (file_type_code == 4) {
      vrtype = 0xc0;
      const_vrec_width += raw_sample_ct * 4;
      pgfip->gflags |= kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent;
    }
    if (S_CAST(uint64_t, raw_variant_ct) * const_vrec_width + pgfip->const_fpos_offset != fsize) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Unexpected .pgen file size (expected %" PRIu64 " bytes).\n", S_CAST(uint64_t, raw_variant_ct) * const_vrec_width + pgfip->const_fpos_offset);
      return kPglRetMalformedInput;
    }
    pgfip->const_vrtype = vrtype;
    pgfip->const_vrec_width = const_vrec_width;
    *pgfi_alloc_cacheline_ct_ptr = 0;
    return kPglRetSuccess;
  }
  if (file_type_code >= 0x11) {
    // todo: 0x11 phase sets
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Third byte of %s does not correspond to a storage mode supported by this version of pgenlib.\n", fname);
    return kPglRetNotYetSupported;
  }
  // plink 2 binary, general-purpose
  pgfip->const_vrtype = UINT32_MAX;
  pgfip->const_vrec_width = 0;
  const uintptr_t alt_allele_ct_byte_ct = (header_ctrl >> 4) & 3;
  if (alt_allele_ct_byte_ct > 1) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: This version of pgenlib does not support >254 alternate alleles for a single variant.\n");
    return kPglRetNotYetSupported;
  }

  // 8 extra bytes per vblock, to support fast random access
  const uintptr_t vblock_ct = DivUp(raw_variant_ct, kPglVblockSize);

  uint64_t vrtype_and_vrec_len_quarterbyte_cost;
  if (header_ctrl & 8) {
    const uint32_t header_ctrl_lowbits = header_ctrl & 15;
    if (header_ctrl_lowbits > 9) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Twelfth byte of %s does not correspond to a format supported by this version of pgenlib.\n", fname);
      return kPglRetNotYetSupported;
    }
    vrtype_and_vrec_len_quarterbyte_cost = header_ctrl_lowbits - 7;
  } else {
    // set this to *2* if true, 0 if false
    const uint32_t phase_or_dosage_present = (header_ctrl >> 1) & 2;
    // vrtype entries = 2 quarterbytes if no phase/dosage, 4 otherwise
    // var_fpos entries = 4 + (4 * (header_ctrl & 3)) quarterbytes
    vrtype_and_vrec_len_quarterbyte_cost = 6 + phase_or_dosage_present + 4 * (header_ctrl & 3);
  }
  pgfip->const_fpos_offset += (raw_sample_ct * vrtype_and_vrec_len_quarterbyte_cost + 3) / 4 + (raw_sample_ct * alt_allele_ct_byte_ct) + (8 * vblock_ct);
  *pgfi_alloc_cacheline_ct_ptr = CountPgfiAllocCachelinesRequired(raw_variant_ct);
  return kPglRetSuccess;
}

static_assert(kPglMaxAltAlleleCt == 254, "Need to update PgfiInitPhase2().");
PglErr PgfiInitPhase2(PgenHeaderCtrl header_ctrl, uint32_t allele_cts_already_loaded, uint32_t nonref_flags_already_loaded, uint32_t use_blockload, uint32_t vblock_idx_start, uint32_t vidx_end, uint32_t* max_vrec_width_ptr, PgenFileInfo* pgfip, unsigned char* pgfi_alloc, uintptr_t* pgr_alloc_cacheline_ct_ptr, char* errstr_buf) {
  // *max_vrec_width_ptr technically only needs to be set in single-variant
  // fread() mode, but its computation is not currently optimized out in the
  // other two modes.

  // possible todo: add option to skip validation when
  // allele_cts/nonref_flags are already loaded.  but let's play it
  // safe for now.
  const uint32_t raw_variant_ct = pgfip->raw_variant_ct;
  const uint32_t const_vrec_width = pgfip->const_vrec_width;
  *pgr_alloc_cacheline_ct_ptr = 0;
  unsigned char loadbuf[kPglVblockSize * 4];
  uintptr_t* allele_idx_offsets_iter = pgfip->allele_idx_offsets;
  uintptr_t prev_allele_idx_offset = 0;
  if (allele_idx_offsets_iter) {
    if (!allele_cts_already_loaded) {
      *allele_idx_offsets_iter = 0;
    } else {
      prev_allele_idx_offset = *allele_idx_offsets_iter;
    }
    ++allele_idx_offsets_iter;
  }
  if (!raw_variant_ct) {
    return kPglRetSuccess;
  }
  const uint32_t nonref_flags_stored = ((header_ctrl >> 6) == 3);
  unsigned char* nonref_flags_iter = R_CAST(unsigned char*, pgfip->nonref_flags);
  const unsigned char* fread_ptr = nullptr;  // maybe-uninitialized warning
  FILE* shared_ff = pgfip->shared_ff;
  if (const_vrec_width) {
    // no allele counts to verify if fixed-width
    // always need workspace_vec
    *pgr_alloc_cacheline_ct_ptr = QuaterCtToCachelineCt(pgfip->raw_sample_ct);
    *max_vrec_width_ptr = const_vrec_width;
#ifdef NO_MMAP
    assert(shared_ff);
#else
    if (!shared_ff) {
      if (use_blockload) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase2() cannot be called with use_blockload set when PgfiInitPhase1() had use_mmap set.\n");
        return kPglRetImproperFunctionCall;
      }
      if ((!(header_ctrl & 192)) || (pgfip->const_vrtype == kPglVrtypePlink1)) {
        return kPglRetSuccess;
      }
      fread_ptr = &(pgfip->block_base[12]);
      const uint32_t nonref_flags_byte_ct = DivUp(raw_variant_ct, CHAR_BIT);
      if (!nonref_flags_already_loaded) {
        if (nonref_flags_stored) {
          memcpy(nonref_flags_iter, fread_ptr, nonref_flags_byte_ct);
        }
        return kPglRetSuccess;
      }
      if (nonref_flags_stored) {
        if (memcmp(nonref_flags_iter, fread_ptr, nonref_flags_byte_ct)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
          return kPglRetInconsistentInput;
        }
        return kPglRetSuccess;
      }
      if (header_ctrl & 64) {
        // all ref
        if (!AllWordsAreZero(pgfip->nonref_flags, BitCtToWordCt(raw_variant_ct))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
          return kPglRetInconsistentInput;
        }
        return kPglRetSuccess;
      }
      // all nonref
      if (!AllBitsAreOne(pgfip->nonref_flags, raw_variant_ct)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
        return kPglRetInconsistentInput;
      }
      return kPglRetSuccess;
    }
#endif
    if (!use_blockload) {
      // using fread() single-variant-at-a-time, need pgr.fread_buf
      *pgr_alloc_cacheline_ct_ptr += DivUp(const_vrec_width, kCacheline);
    }
    if ((!(header_ctrl & 192)) || (pgfip->const_vrtype == kPglVrtypePlink1)) {
      return kPglRetSuccess;
    }
    if ((header_ctrl >> 6) == 1) {
      // all ref
      if (nonref_flags_already_loaded) {
        if (!AllWordsAreZero(pgfip->nonref_flags, BitCtToWordCt(raw_variant_ct))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
          return kPglRetInconsistentInput;
        }
      }
      return kPglRetSuccess;
    }
    if ((header_ctrl >> 6) == 2) {
      // all nonref
      if (nonref_flags_already_loaded) {
        if (!AllBitsAreOne(pgfip->nonref_flags, raw_variant_ct)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
          return kPglRetInconsistentInput;
        }
      }
      return kPglRetSuccess;
    }
    // _last more useful than _end iff we just refer to the number of elements
    // in the block and have no use for a _stop pointer
    unsigned char* nonref_flags_last = &(nonref_flags_iter[((raw_variant_ct - 1) / (kPglVblockSize * 32)) * (kPglVblockSize * 4)]);
    uint32_t cur_byte_ct = kPglVblockSize * 4;
    while (1) {
      if (nonref_flags_iter >= nonref_flags_last) {
        if (nonref_flags_iter > nonref_flags_last) {
          return kPglRetSuccess;
        }
        cur_byte_ct = 1 + ((raw_variant_ct - 1) % (kPglVblockSize * 32)) / CHAR_BIT;
      }
      unsigned char* loadptr = nonref_flags_already_loaded? loadbuf : nonref_flags_iter;
      if (!fread_unlocked(loadptr, cur_byte_ct, 1, shared_ff)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
        return kPglRetReadFail;
      }
      if (nonref_flags_already_loaded) {
        if (memcmp(nonref_flags_iter, loadbuf, cur_byte_ct)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
          return kPglRetInconsistentInput;
        }
      }
      nonref_flags_iter = &(nonref_flags_iter[cur_byte_ct]);
    }
  }

  const uint32_t raw_sample_ct = pgfip->raw_sample_ct;
  unsigned char* vrtypes_iter = pgfi_alloc;
  pgfip->vrtypes = vrtypes_iter;
  uint64_t* var_fpos_iter = R_CAST(uint64_t*, &(vrtypes_iter[RoundUpPow2(raw_variant_ct + 1, kCacheline)]));
  pgfip->var_fpos = var_fpos_iter;
  uint32_t vblock_ct_m1 = (raw_variant_ct - 1) / kPglVblockSize;
  uint32_t max_vrec_width = 0;
  uint64_t cur_fpos;
#ifdef NO_MMAP
  assert(shared_ff);
#else
  if (!shared_ff) {
    if (use_blockload) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase2() cannot be called with use_blockload set when PgfiInitPhase1() had use_mmap set.\n");
      return kPglRetImproperFunctionCall;
    }
    fread_ptr = &(pgfip->block_base[12 + 8 * vblock_idx_start]);
    memcpy(&cur_fpos, fread_ptr, sizeof(int64_t));
    fread_ptr = &(fread_ptr[(vblock_ct_m1 + 1 - vblock_idx_start) * sizeof(int64_t)]);
  } else {
#endif
    if (vblock_idx_start) {
      if (fseeko(shared_ff, vblock_idx_start * sizeof(int64_t), SEEK_CUR)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
        return kPglRetReadFail;
      }
    }
    if (!fread_unlocked(&cur_fpos, sizeof(int64_t), 1, shared_ff)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
    // May also need to load the rest of these values in the future, if we want
    // to support dynamic insertion into a memory-mapped file.  But skip them
    // for now.
    if (fseeko(shared_ff, (vblock_ct_m1 - vblock_idx_start) * sizeof(int64_t), SEEK_CUR)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
#ifndef NO_MMAP
  }
#endif
  const uint32_t vrtype_and_fpos_storage = header_ctrl & 15;
  const uint32_t alt_allele_ct_byte_ct = (header_ctrl >> 4) & 3;
  if (alt_allele_ct_byte_ct) {
    assert(alt_allele_ct_byte_ct == 1);
    if (!allele_idx_offsets_iter) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: pgfip->allele_idx_offsets must be allocated before PgfiInitPhase2() is called.\n");
      return kPglRetImproperFunctionCall;
    }
  }
  uint32_t vblock_idx = vblock_idx_start;
  vblock_ct_m1 = (vidx_end - 1) / kPglVblockSize;
  if (vblock_idx) {
    uintptr_t header_vblock_byte_ct = kPglVblockSize * alt_allele_ct_byte_ct;
    if (nonref_flags_stored) {
      header_vblock_byte_ct += kPglVblockSize / CHAR_BIT;
    }
    if (vrtype_and_fpos_storage & 8) {
      header_vblock_byte_ct += kPglVblockSize >> (10 - vrtype_and_fpos_storage);
    } else {
      if (!(vrtype_and_fpos_storage & 4)) {
        header_vblock_byte_ct += kPglVblockSize / 2;
      } else {
        header_vblock_byte_ct += kPglVblockSize;
      }
      header_vblock_byte_ct += kPglVblockSize * (1 + (vrtype_and_fpos_storage & 3));
    }
#ifndef NO_MMAP
    if (!shared_ff) {
      fread_ptr = &(fread_ptr[header_vblock_byte_ct * S_CAST(uint64_t, vblock_idx)]);
    } else {
#endif
      if (fseeko(shared_ff, header_vblock_byte_ct * S_CAST(uint64_t, vblock_idx), SEEK_CUR)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
        return kPglRetReadFail;
      }
#ifndef NO_MMAP
    }
#endif
  }
  uint32_t cur_vblock_variant_ct = kPglVblockSize;
  uint32_t max_alt_allele_ct = 1;
  while (1) {
    if (vblock_idx >= vblock_ct_m1) {
      if (vblock_idx > vblock_ct_m1) {
        // finish up
#ifndef NO_MMAP
        // now > instead of != to allow additional information to be stored
        // between header and first variant record
        if (!shared_ff) {
          if (S_CAST(uintptr_t, fread_ptr - pgfip->block_base) > pgfip->var_fpos[0]) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid .pgen header.\n");
            return kPglRetMalformedInput;
          }
        } else {
#endif
          if (S_CAST(uint64_t, ftello(shared_ff)) > pgfip->var_fpos[0]) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid .pgen header.\n");
            return kPglRetMalformedInput;
          }
#ifndef NO_MMAP
        }
#endif
        pgfip->var_fpos[vidx_end] = cur_fpos;
        pgfip->max_alt_allele_ct = max_alt_allele_ct;
        // if difflist/LD might be present, scan for them in a way that's
        // likely to terminate quickly
        PgenGlobalFlags new_gflags = kfPgenGlobal0;
        if (vrtype_and_fpos_storage < 8) {
          uintptr_t* vrtypes_alias_start = R_CAST(uintptr_t*, pgfip->vrtypes);
          uintptr_t* vrtypes_alias_end = &(vrtypes_alias_start[DivUp(vidx_end, kBytesPerWord)]);
          if (vblock_idx_start) {
            vrtypes_alias_start = &(vrtypes_alias_start[vblock_idx_start * (kPglVblockSize / kBytesPerWord)]);
          }
          uintptr_t* vrtypes_alias_iter = vrtypes_alias_start;
          const uint32_t trailing_byte_ct = vidx_end & (kBytesPerWord - 1);
          if (trailing_byte_ct) {
            vrtypes_alias_end[-1] = bzhi(vrtypes_alias_end[-1], trailing_byte_ct * CHAR_BIT);
          }
          for (; vrtypes_alias_iter < vrtypes_alias_end; ++vrtypes_alias_iter) {
            const uintptr_t cur_word = *vrtypes_alias_iter;
            const uintptr_t cur_word_shifted = cur_word >> 1;
            // check if any vrtype has bit 1 set and bit 2 clear
            if (cur_word & (~cur_word_shifted) & (2 * kMask0101)) {
              new_gflags |= kfPgenGlobalLdCompressionPresent | kfPgenGlobalDifflistOrLdPresent;
              break;
            }
            if (cur_word & (5 * kMask0101)) {
              // this catches onebit
              new_gflags |= kfPgenGlobalDifflistOrLdPresent;
            }
          }
          if (!(vrtype_and_fpos_storage & 3)) {
            // 1 byte per vrec_len entry, don't bother to determine true
            // maximum
            max_vrec_width = 255;
          }
          if (vrtype_and_fpos_storage & 4) {
            // likely for one of {hphase, dosage} to be present without the
            // other; make this scan faster in that case, at the cost of
            // failing to early-exit when both are present
            uintptr_t or_word = 0;  // just bitwise-or everything together
            for (vrtypes_alias_iter = vrtypes_alias_start; vrtypes_alias_iter < vrtypes_alias_end; ++vrtypes_alias_iter) {
              or_word |= *vrtypes_alias_iter;
            }
            if (or_word & (0x10 * kMask0101)) {
              new_gflags |= kfPgenGlobalHardcallPhasePresent;
            }
            if (or_word & (0x60 * kMask0101)) {
              new_gflags |= kfPgenGlobalDosagePresent;
              if (or_word & (0x80 * kMask0101)) {
                new_gflags |= kfPgenGlobalDosagePhasePresent;
              }
            }
          }
          pgfip->gflags |= new_gflags;
        } else {
          // just assume worst case here.  the funny-looking
          // (vrtype_and_fpos_storage * 12) - 93 expression evaluates to 3 for
          // the 2-bit case and 15 for the 4-bit case.
          assert(vrtype_and_fpos_storage < 10);
          max_vrec_width = QuaterCtToByteCt(raw_sample_ct) + (vrtype_and_fpos_storage * 12) - 93;
        }
        *pgr_alloc_cacheline_ct_ptr = CountPgrAllocCachelinesRequired(raw_sample_ct, new_gflags, max_alt_allele_ct, (shared_ff && (!use_blockload))? max_vrec_width : 0);
        *max_vrec_width_ptr = max_vrec_width;
        return kPglRetSuccess;
      }
      cur_vblock_variant_ct = ModNz(vidx_end, kPglVblockSize);
    }
    ++vblock_idx;
    // 1. handle vrtypes and var_fpos.
    if (vrtype_and_fpos_storage & 8) {
      // vrtype_and_fpos_storage == 8 -> 2-bit storage -> right-shift 2
      //                         == 9 -> 4-bit storage -> right-shift 1
      const uint32_t log2_entry_bit_width = vrtype_and_fpos_storage - 7;
      const uint32_t entry_bit_width = log2_entry_bit_width * 2;
      const uintptr_t entry_mask = (1 << entry_bit_width) - 1;
      const uint32_t log2_entries_per_word = kBitsPerWordLog2 - log2_entry_bit_width;
      const uintptr_t base_vrec_len = QuaterCtToByteCt(raw_sample_ct);
      const uint32_t cur_byte_ct = 1 + ((cur_vblock_variant_ct - 1) >> (10 - vrtype_and_fpos_storage));
      uint32_t block_len = 1 << log2_entries_per_word;
      const uintptr_t* loadbuf_iter;
#ifdef __arm__
#  error "Unaligned accesses in PgfiInitPhase2()."
#endif
#ifndef NO_MMAP
      if (!shared_ff) {
        loadbuf_iter = R_CAST(const uintptr_t*, fread_ptr);
        fread_ptr = &(fread_ptr[cur_byte_ct]);
      } else {
#endif
        if (!fread_unlocked(loadbuf, cur_byte_ct, 1, shared_ff)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
          return kPglRetReadFail;
        }
        loadbuf_iter = R_CAST(const uintptr_t*, loadbuf);
#ifndef NO_MMAP
      }
#endif
      uint32_t cur_vblock_idx = 0;
      uint32_t cur_vblock_idx_stop = 0;
      while (1) {
        cur_vblock_idx_stop += block_len;
        if (cur_vblock_idx_stop > cur_vblock_variant_ct) {
          if (cur_vblock_idx == cur_vblock_variant_ct) {
            break;
          }
          cur_vblock_idx_stop = cur_vblock_variant_ct;
        }
        uintptr_t input_word = *loadbuf_iter++;
        for (; cur_vblock_idx < cur_vblock_idx_stop; ++cur_vblock_idx) {
          const uintptr_t input_word_masked = input_word & entry_mask;
          *vrtypes_iter++ = input_word_masked? 8 : 0;
          *var_fpos_iter++ = cur_fpos;
          cur_fpos += base_vrec_len + input_word_masked;
          input_word >>= entry_bit_width;
        }
      }
    } else {
      if (!(vrtype_and_fpos_storage & 4)) {
        // no phase or dosage present, 4-bit vrtypes
        const uint32_t cur_byte_ct = DivUp(cur_vblock_variant_ct, 2);
#ifndef NO_MMAP
        if (shared_ff) {
#endif
          if (!fread_unlocked(loadbuf, cur_byte_ct, 1, shared_ff)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
            return kPglRetReadFail;
          }
          fread_ptr = loadbuf;
#ifndef NO_MMAP
        }
#endif
        const uint32_t word_write_ct = DivUp(cur_vblock_variant_ct, kBytesPerWord);
        uintptr_t* vrtypes_alias_fullword = R_CAST(uintptr_t*, vrtypes_iter);
        const Halfword* loadbuf_alias_halfword = R_CAST(const Halfword*, fread_ptr);
        for (uint32_t widx = 0; widx < word_write_ct; ++widx) {
          uintptr_t ww = loadbuf_alias_halfword[widx];
#ifdef __LP64__
          ww = (ww | (ww << 16)) & kMask0000FFFF;
#endif
          ww = (ww | (ww << 8)) & kMask00FF;
          vrtypes_alias_fullword[widx] = (ww | (ww << 4)) & kMask0F0F;
        }
        const uint32_t last_word_byte_ct = cur_vblock_variant_ct % kBytesPerWord;
        vrtypes_iter = &(vrtypes_iter[cur_vblock_variant_ct]);
        if (last_word_byte_ct) {
          memset(vrtypes_iter, 0, kBytesPerWord - last_word_byte_ct);
        } else {
          // must guarantee a trailing zero for is_ldbase check to work
          vrtypes_iter[0] = 0;
        }
#ifndef NO_MMAP
        if (!shared_ff) {
          fread_ptr = &(fread_ptr[cur_byte_ct]);
        }
#endif
      } else {
        // phase and dosage
#ifndef NO_MMAP
        if (shared_ff) {
#endif
          if (!fread_unlocked(vrtypes_iter, cur_vblock_variant_ct, 1, shared_ff)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
            return kPglRetReadFail;
          }
#ifndef NO_MMAP
        } else {
          memcpy(vrtypes_iter, fread_ptr, cur_vblock_variant_ct);
        }
#endif
        const uint32_t last_word_byte_ct = cur_vblock_variant_ct % kBytesPerWord;
        vrtypes_iter = &(vrtypes_iter[cur_vblock_variant_ct]);
        if (last_word_byte_ct) {
          memset(vrtypes_iter, 0, kBytesPerWord - last_word_byte_ct);
        } else {
          // must guarantee a trailing zero for is_ldbase check to work
          vrtypes_iter[0] = 0;
        }
#ifndef NO_MMAP
        if (!shared_ff) {
          fread_ptr = &(fread_ptr[cur_vblock_variant_ct]);
        }
#endif
      }
      const uint32_t bytes_per_entry = 1 + (vrtype_and_fpos_storage & 3);
      const uint32_t cur_byte_ct = cur_vblock_variant_ct * bytes_per_entry;
#ifndef NO_MMAP
      if (shared_ff) {
#endif
        if (!fread_unlocked(loadbuf, cur_byte_ct, 1, shared_ff)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
          return kPglRetReadFail;
        }
        fread_ptr = loadbuf;
#ifndef NO_MMAP
      }
#endif
      if (bytes_per_entry == 1) {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx < cur_vblock_variant_ct; ++cur_vblock_vidx) {
          var_fpos_iter[cur_vblock_vidx] = cur_fpos;
          uint32_t cur_vrec_len = fread_ptr[cur_vblock_vidx];
          cur_fpos += cur_vrec_len;
          // no need for correct max_vrec_width
        }
      } else if (bytes_per_entry == 2) {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx < cur_vblock_variant_ct; ++cur_vblock_vidx) {
          var_fpos_iter[cur_vblock_vidx] = cur_fpos;
          uint32_t cur_vrec_len = R_CAST(const uint16_t*, fread_ptr)[cur_vblock_vidx];
          cur_fpos += cur_vrec_len;
          if (cur_vrec_len > max_vrec_width) {
            // todo: check whether we're better off just assuming 2^16 - 1
            max_vrec_width = cur_vrec_len;
          }
        }
      } else if (bytes_per_entry == 3) {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx < cur_vblock_variant_ct; ++cur_vblock_vidx) {
          var_fpos_iter[cur_vblock_vidx] = cur_fpos;
          uint32_t cur_vrec_len = (*R_CAST(const uint32_t*, &(fread_ptr[cur_vblock_vidx * bytes_per_entry]))) & 0xffffff;
          cur_fpos += cur_vrec_len;
          if (cur_vrec_len > max_vrec_width) {
            max_vrec_width = cur_vrec_len;
          }
        }
      } else {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx < cur_vblock_variant_ct; ++cur_vblock_vidx) {
          var_fpos_iter[cur_vblock_vidx] = cur_fpos;
          uint32_t cur_vrec_len = R_CAST(const uint32_t*, fread_ptr)[cur_vblock_vidx];
          cur_fpos += cur_vrec_len;
          if (cur_vrec_len > max_vrec_width) {
            max_vrec_width = cur_vrec_len;
          }
        }
#ifdef __LP64__
        if (max_vrec_width > kPglMaxBytesPerVariant) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid .pgen header.\n");
          return kPglRetMalformedInput;
        }
#else
        if (max_vrec_width > kMaxBytesPerIO) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Variant records too large for 32-bit pgenlib.\n");
          return kPglRetNomem;
        }
#endif
      }
      var_fpos_iter = &(var_fpos_iter[cur_vblock_variant_ct]);
#ifndef NO_MMAP
      if (!shared_ff) {
        fread_ptr = &(fread_ptr[cur_byte_ct]);
      }
#endif
    }
    // 2. allele counts?
    if (alt_allele_ct_byte_ct) {
      assert(alt_allele_ct_byte_ct == 1);
#ifndef NO_MMAP
      if (shared_ff) {
#endif
        if (!fread_unlocked(loadbuf, cur_vblock_variant_ct * alt_allele_ct_byte_ct, 1, shared_ff)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
          return kPglRetReadFail;
        }
        fread_ptr = loadbuf;
#ifndef NO_MMAP
      }
#endif
      if (allele_cts_already_loaded) {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx < cur_vblock_variant_ct; ++cur_vblock_vidx) {
          uintptr_t cur_allele_idx_offset = allele_idx_offsets_iter[cur_vblock_vidx];
          uint32_t cur_alt_allele_ct = fread_ptr[cur_vblock_vidx];
          if ((cur_allele_idx_offset - prev_allele_idx_offset) != (cur_alt_allele_ct + 1)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded allele_idx_offsets do not match values in .pgen file.\n");
            return kPglRetInconsistentInput;
          }
          prev_allele_idx_offset = cur_allele_idx_offset;
          if (cur_alt_allele_ct > max_alt_allele_ct) {
            max_alt_allele_ct = cur_alt_allele_ct;
          }
        }
      } else {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx < cur_vblock_variant_ct; ++cur_vblock_vidx) {
          uint32_t cur_alt_allele_ct = fread_ptr[cur_vblock_vidx];
          prev_allele_idx_offset += cur_alt_allele_ct + 1;
          allele_idx_offsets_iter[cur_vblock_vidx] = prev_allele_idx_offset;
          if (cur_alt_allele_ct > max_alt_allele_ct) {
            max_alt_allele_ct = cur_alt_allele_ct;
          }
        }
      }
      allele_idx_offsets_iter = &(allele_idx_offsets_iter[cur_vblock_variant_ct]);
#ifndef NO_MMAP
      if (!shared_ff) {
        fread_ptr = &(fread_ptr[cur_vblock_variant_ct * alt_allele_ct_byte_ct]);
      }
#endif
    }
    // 3. nonref flags?
    if (nonref_flags_stored) {
      const uint32_t cur_byte_ct = DivUp(cur_vblock_variant_ct, CHAR_BIT);
#ifndef NO_MMAP
      if (!shared_ff) {
        if (nonref_flags_already_loaded) {
          if (memcmp(nonref_flags_iter, fread_ptr, cur_byte_ct)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
            return kPglRetInconsistentInput;
          }
        } else {
          memcpy(nonref_flags_iter, fread_ptr, cur_byte_ct);
        }
        fread_ptr = &(fread_ptr[cur_byte_ct]);
      } else {
#endif
        unsigned char* loadptr = nonref_flags_already_loaded? loadbuf : nonref_flags_iter;
        if (!fread_unlocked(loadptr, cur_byte_ct, 1, shared_ff)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
          return kPglRetReadFail;
        }
        if (nonref_flags_already_loaded) {
          if (memcmp(nonref_flags_iter, loadbuf, cur_byte_ct)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
            return kPglRetInconsistentInput;
          }
        }
#ifndef NO_MMAP
      }
#endif
      nonref_flags_iter = &(nonref_flags_iter[cur_byte_ct]);
    }
  }
}

uint32_t GetLdbaseVidx(const unsigned char* vrtypes, uint32_t cur_vidx) {
  const uintptr_t* vrtypes_walias = R_CAST(const uintptr_t*, vrtypes);
  const uint32_t cur_vidx_orig_remainder = cur_vidx % kBytesPerWord;
  uint32_t vidx_word_idx = (cur_vidx - 1) / kBytesPerWord;
  uintptr_t cur_vrtypes_word = vrtypes_walias[vidx_word_idx];
  if (cur_vidx_orig_remainder) {
    // make sure we don't detect a byte after the current position.
    cur_vrtypes_word = bzhi(cur_vrtypes_word, CHAR_BIT * cur_vidx_orig_remainder);
    cur_vrtypes_word |= (kMask0101 * 2) << (CHAR_BIT * cur_vidx_orig_remainder);
  }
  while (1) {
    // ((bit 2) OR (NOT bit 1)) for each byte.  (possible experiment: see if
    // the same assembly is generated if this expression is rewritten to use
    // ands/nots.)
    const uintptr_t detect_non_ld_word = ((cur_vrtypes_word >> 1) | (~cur_vrtypes_word)) & (kMask0101 * 2);
    if (detect_non_ld_word) {
      // find the highest-order set bit in detect_non_ld_word; this corresponds
      // to the last non-LD-compressed byte (assuming little-endian).
      const uint32_t new_ldbase_vidx_loworder = kBytesPerWord - 1 - (CLZLU(detect_non_ld_word) / CHAR_BIT);
      return (vidx_word_idx * kBytesPerWord) + new_ldbase_vidx_loworder;
    }
    // everything LD-compressed in the current block.  move back 8 bytes in the
    // array (or 4-bytes for 32-bit build).
    cur_vrtypes_word = vrtypes_walias[--vidx_word_idx];
  }
}

uint64_t PgfiMultireadGetCachelineReq(const uintptr_t* variant_include, const PgenFileInfo* pgfip, uint32_t variant_ct, uint32_t block_size) {
  // if block_size < kPglVblockSize, it should be a power of 2 (to avoid
  // unnecessary vblock crossing), but that's not required.
  const uint32_t raw_variant_ct = pgfip->raw_variant_ct;
  if (variant_ct == raw_variant_ct) {
    variant_include = nullptr;
  }
  uint32_t block_ct_m1 = 0;
  if (raw_variant_ct < block_size) {
    block_size = raw_variant_ct;
  } else {
    block_ct_m1 = (raw_variant_ct - 1) / block_size;
  }
  const uint64_t* var_fpos = pgfip->var_fpos;
  if ((!variant_include) && (!var_fpos)) {
    return DivUpU64(S_CAST(uint64_t, pgfip->const_vrec_width) * block_size, kCacheline);
  }
  uint64_t max_block_byte_ct = 0;
  uint32_t max_block_variant_ct = 0;
  uint32_t block_idx = 0;
  while (1) {
    uint32_t variant_uidx_start = block_idx * block_size;
    uint32_t variant_uidx_end = variant_uidx_start + block_size;
    if (block_idx >= block_ct_m1) {
      if (block_idx > block_ct_m1) {
        break;
      }
      variant_uidx_end = raw_variant_ct;
    }
    if (variant_include) {
      variant_uidx_start = FindFirst1BitFromBounded(variant_include, variant_uidx_start, variant_uidx_end);
      if (variant_uidx_start == variant_uidx_end) {
        ++block_idx;
        continue;
      }
      variant_uidx_end = 1 + FindLast1BitBefore(variant_include, variant_uidx_end);
    }
    if (var_fpos) {
      if (pgfip->vrtypes && ((pgfip->vrtypes[variant_uidx_start] & 6) == 2)) {
        // need to start loading from LD-buddy
        variant_uidx_start = GetLdbaseVidx(pgfip->vrtypes, variant_uidx_start);
      }
      uint64_t cur_block_byte_ct = var_fpos[variant_uidx_end] - var_fpos[variant_uidx_start];
      if (cur_block_byte_ct > max_block_byte_ct) {
        max_block_byte_ct = cur_block_byte_ct;
      }
    } else {
      // no LD compression here
      const uint32_t cur_block_variant_ct = variant_uidx_end - variant_uidx_start;
      if (cur_block_variant_ct > max_block_variant_ct) {
        max_block_variant_ct = cur_block_variant_ct;
        if (cur_block_variant_ct == block_size) {
          // no larger value possible, terminate search
          break;
        }
      }
    }
    ++block_idx;
  }
  if (!var_fpos) {
    max_block_byte_ct = max_block_variant_ct * S_CAST(uint64_t, pgfip->const_vrec_width);
  }
  return DivUpU64(max_block_byte_ct, kCacheline);
}

PglErr PgfiMultiread(const uintptr_t* variant_include, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t load_variant_ct, PgenFileInfo* pgfip) {
  // we could permit 0, but that encourages lots of unnecessary thread wakeups
  assert(load_variant_ct);
  if (variant_include) {
    FindFirst1BitFromU32(variant_include, &variant_uidx_start);
  }
  assert(variant_uidx_start < pgfip->raw_variant_ct);
  uint64_t block_offset;
  if (pgfip->vrtypes && ((pgfip->vrtypes[variant_uidx_start] & 6) == 2)) {
    // need to start loading from LD-buddy
    // assume for now that we can't skip any variants between the LD-buddy and
    // the actual first variant; should remove this assumption later
    block_offset = pgfip->var_fpos[GetLdbaseVidx(pgfip->vrtypes, variant_uidx_start)];
  } else {
    block_offset = GetPgfiFpos(pgfip, variant_uidx_start);
  }
  pgfip->block_offset = block_offset;
  uint64_t next_read_start_fpos = block_offset;
  // break this up into multiple freads whenever this lets us skip an entire
  // disk block
  // (possible todo: make the disk block size a parameter of this function)
  do {
    const uint64_t cur_read_start_fpos = next_read_start_fpos;
    uint32_t cur_read_uidx_end;
    uint64_t cur_read_end_fpos;
    while (1) {
      cur_read_uidx_end = variant_uidx_end;
      if (cur_read_uidx_end - variant_uidx_start == load_variant_ct) {
        cur_read_end_fpos = GetPgfiFpos(pgfip, cur_read_uidx_end);
        load_variant_ct = 0;
        break;
      }
      cur_read_uidx_end = FindFirst0BitFrom(variant_include, variant_uidx_start);
      cur_read_end_fpos = GetPgfiFpos(pgfip, cur_read_uidx_end);
      load_variant_ct -= cur_read_uidx_end - variant_uidx_start;
      if (!load_variant_ct) {
        break;
      }
      variant_uidx_start = FindFirst1BitFrom(variant_include, cur_read_uidx_end);
      next_read_start_fpos = GetPgfiFpos(pgfip, variant_uidx_start);
      if (pgfip->vrtypes && ((pgfip->vrtypes[variant_uidx_start] & 6) == 2)) {
        const uint32_t variant_read_uidx_start = GetLdbaseVidx(pgfip->vrtypes, variant_uidx_start);
        if (variant_read_uidx_start <= cur_read_uidx_end) {
          continue;
        }
        next_read_start_fpos = pgfip->var_fpos[variant_read_uidx_start];
      }
      // bugfix: can't use do..while, since previous "continue" needs to skip
      // this check
      if (RoundDownPow2U64(cur_read_end_fpos + kDiskBlockSize + 1LLU, kDiskBlockSize) < RoundDownPow2U64(next_read_start_fpos, kDiskBlockSize)) {
        // minor bugfix (7 Jul 2017): break, not continue
        break;
      }
    }
    if (fseeko(pgfip->shared_ff, cur_read_start_fpos, SEEK_SET)) {
      return kPglRetReadFail;
    }
    uintptr_t len = cur_read_end_fpos - cur_read_start_fpos;
    if (fread_checked(K_CAST(unsigned char*, &(pgfip->block_base[cur_read_start_fpos - block_offset])), len, pgfip->shared_ff)) {
      return kPglRetReadFail;
    }
  } while (load_variant_ct);
  return kPglRetSuccess;
}


void PreinitPgr(PgenReader* pgrp) {
  pgrp->ff = nullptr;
}

PglErr PgrInit(const char* fname, uint32_t max_vrec_width, PgenFileInfo* pgfip, PgenReader* pgrp, unsigned char* pgr_alloc) {
  // See CountPgrAllocCachelinesRequired().
  // Could add a debug mode.

  // Mode 1 (mmap): block_base initialized, shared_ff == nullptr.  fname must
  //   be nullptr.
  // Mode 2 (block-fread): block_base initialized, shared_ff != nullptr.  fname
  //   must be nullptr.
  // Mode 3 (per-variant fread): block_base == nullptr.  fname must be
  //   non-null, though it isn't actually referenced during the first
  //   PgenReader initialization (instead shared_ff is moved).
  unsigned char* pgr_alloc_iter = pgr_alloc;
  if (pgfip->block_base != nullptr) {
    if (fname != nullptr) {
      return kPglRetImproperFunctionCall;
    }
    pgrp->ff = nullptr;  // make sure PgrCleanup() doesn't break
  } else {
    if (pgfip->shared_ff != nullptr) {
      if (fname == nullptr) {
        return kPglRetImproperFunctionCall;
      }
      // move instead of close/reopen.
      pgrp->ff = pgfip->shared_ff;
      pgfip->shared_ff = nullptr;
    } else {
      pgrp->ff = fopen(fname, FOPEN_RB);
      if (!pgrp->ff) {
        return kPglRetOpenFail;
      }
    }
    // now that arbitrary info can be stored between header and first variant
    // record, always seek.
    uint64_t seek_pos;
    if (pgfip->var_fpos) {
      seek_pos = pgfip->var_fpos[0];
    } else {
      seek_pos = pgfip->const_fpos_offset;
    }
    if (fseeko(pgrp->ff, seek_pos, SEEK_SET)) {
      return kPglRetReadFail;
    }
  }
  pgrp->fi = *pgfip;  // struct copy
  if (fname) {
    // Mode 3 per-reader load buffer
    pgrp->fread_buf = pgr_alloc_iter;
    pgr_alloc_iter = &(pgr_alloc_iter[RoundUpPow2(max_vrec_width, kCacheline)]);
  }
  pgrp->fp_vidx = 0;
  pgrp->ldbase_vidx = UINT32_MAX;
  pgrp->ldbase_stypes = kfPgrLdcache0;
  pgrp->ldbase_genovec = nullptr;
  pgrp->ldbase_raregeno = nullptr;
  pgrp->ldbase_difflist_sample_ids = nullptr;

  const PgenGlobalFlags gflags = pgrp->fi.gflags;
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  pgrp->workspace_vec = R_CAST(uintptr_t*, pgr_alloc_iter);
  const uint32_t genovec_bytes_req = QuaterCtToCachelineCt(raw_sample_ct) * kCacheline;
  pgr_alloc_iter = &(pgr_alloc_iter[genovec_bytes_req]);
  const uint32_t bitvec_bytes_req = BitCtToCachelineCt(raw_sample_ct) * kCacheline;
  const uint32_t ld_compression_present = (gflags / kfPgenGlobalLdCompressionPresent) & 1;
  if (gflags & kfPgenGlobalDifflistOrLdPresent) {
    const uint32_t max_difflist_entry_ct_base = (raw_sample_ct / kPglMaxDifflistLenDivisor);
    // const uint32_t max_difflist_entry_ct = max_difflist_entry_ct_base * (1 + ld_compression_present);

    pgrp->workspace_raregeno_vec = R_CAST(uintptr_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[QuaterCtToCachelineCt(max_difflist_entry_ct_base) * kCacheline]);

    pgrp->workspace_difflist_sample_ids = R_CAST(uint32_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[(1 + (max_difflist_entry_ct_base / kInt32PerCacheline)) * (kCacheline * k1LU)]);

    pgrp->workspace_raregeno_tmp_loadbuf = R_CAST(uintptr_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[QuaterCtToCachelineCt(max_difflist_entry_ct_base) * kCacheline]);

    pgrp->workspace_difflist_sample_ids_tmp = R_CAST(uint32_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[(1 + (max_difflist_entry_ct_base / kInt32PerCacheline)) * (kCacheline * k1LU)]);

    if (ld_compression_present) {
      pgrp->ldbase_genovec = R_CAST(uintptr_t*, pgr_alloc_iter);
      pgr_alloc_iter = &(pgr_alloc_iter[genovec_bytes_req]);

      pgrp->ldbase_raregeno = R_CAST(uintptr_t*, pgr_alloc_iter);
      pgr_alloc_iter = &(pgr_alloc_iter[QuaterCtToCachelineCt(max_difflist_entry_ct_base) * kCacheline]);

      pgrp->ldbase_difflist_sample_ids = R_CAST(uint32_t*, pgr_alloc_iter);
      pgr_alloc_iter = &(pgr_alloc_iter[(1 + (max_difflist_entry_ct_base / kInt32PerCacheline)) * (kCacheline * k1LU)]);
    }
  } else {
    pgrp->workspace_raregeno_vec = nullptr;
    pgrp->workspace_difflist_sample_ids = nullptr;
    pgrp->workspace_raregeno_tmp_loadbuf = nullptr;
    pgrp->workspace_difflist_sample_ids_tmp = nullptr;
  }
  const uint32_t max_alt_allele_ct = pgrp->fi.max_alt_allele_ct;
  if (max_alt_allele_ct > 1) {
    pgrp->workspace_aux1_nonmissing_vec = R_CAST(uintptr_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);

    uintptr_t aux1_allele_bytect = GetAux1AlleleByteCt(max_alt_allele_ct, raw_sample_ct);
    if (aux1_allele_bytect > kPglMaxBytesPerVariant) {
      aux1_allele_bytect = kPglMaxBytesPerVariant;
    }
    pgrp->workspace_aux1_code_vec = R_CAST(uintptr_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[RoundUpPow2(aux1_allele_bytect, kCacheline)]);

    pgrp->workspace_ambig_sample_ids = R_CAST(uint32_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[Int32CtToCachelineCt(raw_sample_ct) * (kCacheline * k1LU)]);
  } else {
    pgrp->workspace_aux1_nonmissing_vec = nullptr;
    pgrp->workspace_aux1_code_vec = nullptr;
    pgrp->workspace_ambig_sample_ids = nullptr;
  }
  pgrp->workspace_all_hets = nullptr;
  pgrp->ldbase_all_hets = nullptr;
  if (gflags & kfPgenGlobalHardcallPhasePresent) {
    pgrp->workspace_all_hets = R_CAST(uintptr_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);
    if (ld_compression_present) {
      pgrp->ldbase_all_hets = R_CAST(uintptr_t*, pgr_alloc_iter);
      pgrp->ldbase_all_hets[(raw_sample_ct - 1) / kBitsPerWord] = 0;
      pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);
    }
  }
  pgrp->workspace_dosage_present = nullptr;
  pgrp->workspace_dosage_phased = nullptr;
  if (gflags & kfPgenGlobalDosagePresent) {
    pgrp->workspace_dosage_present = R_CAST(uintptr_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);
    if (gflags & kfPgenGlobalDosagePhasePresent) {
      pgrp->workspace_dosage_phased = R_CAST(uintptr_t*, pgr_alloc_iter);
    }
    // pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);
  }
  return kPglRetSuccess;
}

void PgrPlink1ToPlink2InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // 00 -> 10, 01 -> 11, 10 -> 01, 11 -> 00
  // new low bit  = [old low] ^ [old high]
  // new high bit = ~[old high]
  // "unsafe" because trailing bits are not zeroed out.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  const VecUL m1 = VCONST_UL(kMask5555);
  const VecUL not_m1 = VCONST_UL(kMaskAAAA);
  VecUL* vptr = R_CAST(VecUL*, genovec);
  for (uint32_t vidx = 0; vidx < vec_ct; vidx++) {
    const VecUL not_cur_vec_high = (~vptr[vidx]) & not_m1;
    vptr[vidx] = (((~vptr[vidx]) & m1) ^ vecul_srli(not_cur_vec_high, 1)) | not_cur_vec_high;
  }
}

void PgrPlink2ToPlink1InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // 00 -> 11, 01 -> 10, 10 -> 00, 11 -> 01
  // new low bit  = [old low] ^ (~[old high])
  // new high bit = ~[old high]
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  const VecUL not_m1 = VCONST_UL(kMaskAAAA);
  VecUL* vptr = R_CAST(VecUL*, genovec);
  for (uint32_t vidx = 0; vidx < vec_ct; vidx++) {
    VecUL cur_vec = vptr[vidx];
    VecUL not_cur_vec_high = (~cur_vec) & not_m1;
    vptr[vidx] = (((~not_m1) & cur_vec) ^ vecul_srli(not_cur_vec_high, 1)) | not_cur_vec_high;
  }
}

PglErr ParseDifflistHeader(const unsigned char* fread_end, uint32_t raw_sample_ct, const unsigned char** fread_pp, uintptr_t* raregeno_buf, const unsigned char** difflist_group_info_ptr, uint32_t* difflist_len_ptr) {
  // Can be used for deltalists as well: pass raregeno_buf == nullptr.
  // Trailing bits of raregeno may not be zeroed out.
  const uint32_t difflist_len = GetVint31(fread_end, fread_pp);
  // moved here to address maybe-uninitialized warnings
  *difflist_group_info_ptr = *fread_pp;
  *difflist_len_ptr = difflist_len;
  if (!difflist_len) {
    return kPglRetSuccess;
  }
  if (difflist_len > raw_sample_ct / kPglMaxDifflistLenDivisor) {
    // automatically catches GetVint31() failure
    return kPglRetMalformedInput;
  }
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  const uint32_t difflist_index_byte_ct = group_ct * (sample_id_byte_ct + 1) - 1;
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < difflist_index_byte_ct) {
    return kPglRetMalformedInput;
  }
  *fread_pp += difflist_index_byte_ct;
  if (!raregeno_buf) {
    // for sample ID lists without 2-bit genotype info, used for sparse dosage
    return kPglRetSuccess;
  }
  const uint32_t raregeno_byte_ct = QuaterCtToByteCt(difflist_len);
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < raregeno_byte_ct) {
    return kPglRetMalformedInput;
  }
  const unsigned char* raregeno_end = &((*fread_pp)[raregeno_byte_ct]);
  // possible todo: just return a pointer to the beginning of the raregeno
  // segment, and let the caller perform this copy
  memcpy(raregeno_buf, *fread_pp, raregeno_byte_ct);
  *fread_pp = raregeno_end;
  return kPglRetSuccess;
}

PglErr ParseAndSaveDifflist(const unsigned char* fread_end, uint32_t raw_sample_ct, const unsigned char** fread_pp, uintptr_t* __restrict raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr) {
  // Appropriate when we need to iterate through the difflist multiple times.
  // Other functions are more efficient if we only need to process the list
  // once.
  // Trailing bits of raregeno may not be zeroed out.
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno, &group_info_iter, difflist_len_ptr);
  uint32_t difflist_remaining = *difflist_len_ptr;
  // todo: check if difflist_len == 0 early exit is a net positive or negative
  // on a few test datasets
  if (reterr || (!difflist_remaining)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  uint32_t* difflist_sample_ids_iter = difflist_sample_ids;
  while (1) {
    const uint32_t* difflist_sample_ids_stop;
    if (difflist_remaining < kPglDifflistGroupSize) {
      if (!difflist_remaining) {
        return kPglRetSuccess;
      }
      difflist_sample_ids_stop = &(difflist_sample_ids_iter[difflist_remaining]);
      difflist_remaining = 0;
    } else {
      difflist_sample_ids_stop = &(difflist_sample_ids_iter[kPglDifflistGroupSize]);
      difflist_remaining -= kPglDifflistGroupSize;
    }
    // can't use uint32_t assignment trick for now since there's a corner case
    // where that would read past the end of the mapped address range
    uintptr_t raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
    group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    while (1) {
#ifndef __LP64__
      // perform more frequent checks in 32-bit build since raw_sample_idx may
      // overflow
      // misses "small negative" malformed input, but it'll catch data
      // corruption with very high probability
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      *difflist_sample_ids_iter++ = raw_sample_idx;
      if (difflist_sample_ids_iter == difflist_sample_ids_stop) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
#ifdef __LP64__
    if (raw_sample_idx >= raw_sample_ct) {
      return kPglRetMalformedInput;
    }
#endif
  }
  return kPglRetSuccess;
}

void GetDifflistAmbigIdsUnsafe(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_len, uint32_t* __restrict ambig_sample_ids, uint32_t* ambig_id_ct_ptr) {
  // assumes trailing bits of raregeno are zeroed out
  const uint32_t difflist_wct = QuaterCtToWordCt(difflist_len);
  uint32_t ambig_id_ct = 0;
  for (uint32_t widx = 0; widx < difflist_wct; ++widx) {
    uintptr_t detect_11 = raregeno[widx] & (raregeno[widx] >> 1) & kMask5555;
    // now detect_11 has a set bit iff the raregeno entry is 0b11
    if (detect_11) {
      const uint32_t* difflist_sample_ids_base = &(difflist_sample_ids[widx * kBitsPerWordD2]);
      do {
        uint32_t difflist_idx_lowbits = CTZLU(detect_11) / 2;
        ambig_sample_ids[ambig_id_ct++] = difflist_sample_ids_base[difflist_idx_lowbits];
        detect_11 &= detect_11 - 1;
      } while (detect_11);
    }
  }
  *ambig_id_ct_ptr = ambig_id_ct;
}

PglErr ParseAndSaveDifflistProperSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t raw_sample_ct, const unsigned char** fread_pp, uintptr_t* __restrict raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr, uintptr_t* __restrict raregeno_workspace) {
  // Requires a PROPER subset, and does not save ambig_sample_ids.  Use the
  // more generic parse_and_save_difflist_subset() if the latter might be
  // needed.
  // Might want to just merge this with ParseAndSaveDifflist() and rename
  // appropriately.
  // Trailing bits of raregeno are zeroed out.
  uint32_t raw_difflist_len;
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno_workspace, &group_info_iter, &raw_difflist_len);
  if (reterr || (!raw_difflist_len)) {
    *difflist_len_ptr = 0;
    return reterr;
  }
  const uint32_t subgroup_idx_last = (raw_difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  uintptr_t* raregeno_workspace_iter = raregeno_workspace;
  uintptr_t* raregeno_iter = raregeno;
  uint32_t* difflist_sample_ids_iter = difflist_sample_ids;

  // technically doesn't need to be initialized, but I have principles
  uintptr_t raw_sample_idx = 0;

  uintptr_t raregeno_word = 0;
  uint32_t subgroup_idx = 0;
  uint32_t subgroup_len_m1 = kBitsPerWordD2 - 1;
  uint32_t difflist_len_lowbits = 0;
  while (1) {
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        if (difflist_len_lowbits) {
          *raregeno_iter = raregeno_word;
        }
        *difflist_len_ptr = S_CAST(uintptr_t, difflist_sample_ids_iter - difflist_sample_ids) + difflist_len_lowbits;
        return kPglRetSuccess;
      }
      subgroup_len_m1 &= raw_difflist_len - 1;
    }
    // We need to consume a new rare genotype word every 32 entries, and pull a
    // raw sample index from the difflist header every 64 entries.  So it's
    // best to make the inner loop have a period of 32 (call this a "subgroup",
    // where "group" refers to a set of 64 entries).
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
#ifdef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    ++subgroup_idx;
    uintptr_t raregeno_workspace_word = *raregeno_workspace_iter++;
    uint32_t raw_difflist_idx_lowbits = 0;
    while (1) {
#ifndef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      if (IsSet(sample_include, raw_sample_idx)) {
        raregeno_word |= ((raregeno_workspace_word >> (2 * raw_difflist_idx_lowbits)) & 3) << (difflist_len_lowbits * 2);
        difflist_sample_ids_iter[difflist_len_lowbits] = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, raw_sample_idx);
        if (difflist_len_lowbits++ == (kBitsPerWordD2 - 1)) {
          *raregeno_iter++ = raregeno_word;
          raregeno_word = 0;
          difflist_len_lowbits = 0;
          difflist_sample_ids_iter = &(difflist_sample_ids_iter[kBitsPerWordD2]);
        }
      }
      if (raw_difflist_idx_lowbits == subgroup_len_m1) {
        break;
      }
      ++raw_difflist_idx_lowbits;
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
  }
}

PglErr ParseAndSaveDifflistSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t raw_sample_ct, const unsigned char** fread_pp, uint32_t* __restrict ambig_sample_ids, uintptr_t* __restrict raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr, uint32_t* __restrict ambig_id_ct_ptr, uintptr_t* __restrict raregeno_workspace) {
  // Generic interface.  sample_include should be nullptr if
  // sample_ct == raw_sample_ct.
  // Trailing bits of raregeno are zeroed out.
  if (!ambig_sample_ids) {
    if (!sample_include) {
      return ParseAndSaveDifflist(fread_end, raw_sample_ct, fread_pp, raregeno, difflist_sample_ids, difflist_len_ptr);
    }
    return ParseAndSaveDifflistProperSubset(fread_end, sample_include, sample_include_cumulative_popcounts, raw_sample_ct, fread_pp, raregeno, difflist_sample_ids, difflist_len_ptr, raregeno_workspace);
  }
  uint32_t raw_difflist_len;
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno_workspace, &group_info_iter, &raw_difflist_len);
  if (reterr || (!raw_difflist_len)) {
    *difflist_len_ptr = 0;
    *ambig_id_ct_ptr = 0;
    return reterr;
  }
  const uint32_t subgroup_idx_last = (raw_difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  uintptr_t* raregeno_workspace_iter = raregeno_workspace;
  uintptr_t* raregeno_iter = raregeno;
  uint32_t* difflist_sample_ids_iter = difflist_sample_ids;

  // technically doesn't need to be initialized, but I have principles
  uintptr_t raw_sample_idx = 0;

  uintptr_t raregeno_word = 0;
  uint32_t ambig_id_ct = 0;
  uint32_t subgroup_idx = 0;
  uint32_t subgroup_len_m1 = kBitsPerWordD2 - 1;
  uint32_t difflist_len_lowbits = 0;
  while (1) {
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        if (difflist_len_lowbits) {
          *raregeno_iter = raregeno_word;
        }
        *difflist_len_ptr = S_CAST(uintptr_t, difflist_sample_ids_iter - difflist_sample_ids) + difflist_len_lowbits;
        *ambig_id_ct_ptr = ambig_id_ct;
        return kPglRetSuccess;
      }
      subgroup_len_m1 &= raw_difflist_len - 1;
    }
    // We need to consume a new rare genotype word every 32 entries, and pull a
    // raw sample index from the difflist header every 64 entries.  So it's
    // best to make the inner loop have a period of 32 (call this a "subgroup",
    // where "group" refers to a set of 64 entries).
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
#ifdef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    ++subgroup_idx;
    uintptr_t raregeno_workspace_word = *raregeno_workspace_iter++;
    uint32_t raw_difflist_idx_lowbits = 0;
    while (1) {
#ifndef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      const uintptr_t cur_geno = raregeno_workspace_word & 3;
      if ((!sample_include) || IsSet(sample_include, raw_sample_idx)) {
        uint32_t sample_idx = raw_sample_idx;
        if (sample_include) {
          sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, raw_sample_idx);
        }
        raregeno_word |= cur_geno << (difflist_len_lowbits * 2);
        difflist_sample_ids_iter[difflist_len_lowbits] = sample_idx;
        if (difflist_len_lowbits++ == (kBitsPerWordD2 - 1)) {
          *raregeno_iter++ = raregeno_word;
          raregeno_word = 0;
          difflist_len_lowbits = 0;
          difflist_sample_ids_iter = &(difflist_sample_ids_iter[kBitsPerWordD2]);
        }
      }
      if (cur_geno == 3) {
        ambig_sample_ids[ambig_id_ct++] = raw_sample_idx;
      }
      if (raw_difflist_idx_lowbits == subgroup_len_m1) {
        break;
      }
      ++raw_difflist_idx_lowbits;
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      raregeno_workspace_word >>= 2;
    }
  }
}

PglErr ParseLdAndMergeDifflistSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict ldbase_raregeno, const uint32_t* __restrict ldbase_difflist_sample_ids, uint32_t ldbase_difflist_len, uintptr_t ldbase_common_geno, uint32_t raw_sample_ct, uint32_t sample_ct, const unsigned char** fread_pp, uint32_t* __restrict ambig_sample_ids, uintptr_t* __restrict merged_raregeno, uint32_t* __restrict merged_difflist_sample_ids, uint32_t* __restrict merged_difflist_len_ptr, uint32_t* __restrict ambig_id_ct_ptr, uintptr_t* __restrict diff_from_ldbase_raregeno_iter) {
  // Used when the ldbase variant was saved as a difflist, and it's useful to
  // process the current variant as a difflist.
  // * If the current variant is multiallelic, 0b11 entries in ldbase_raregeno
  //   do NOT have associated aux1 entries; only freshly loaded 0b11 values do.
  //   ambig_sample_ids keeps track of this.
  // * ambig_sample_ids should be nullptr if the current variant is not
  //   multiallelic.  (Hence its positioning in the argument list: it's in/out,
  //   everything after it is a pure outparemeter.)
  // * ambig_sample_ids is NOT subsetted; otherwise it wouldn't support
  //   subsequent loading of the aux1 data track.
  // * Assumes ldbase_difflist_sample_ids[ldbase_difflist_len]==sample_ct.
  // * Assumes sample_include == nullptr if no subsetting needed.  (Otherwise,
  //   it'll still work, but performance will be slower.)
  // Trailing bits of merged_raregeno may not be zeroed out.
  // Caller is responsible for inverting ldbase_common_geno and merged_raregeno
  // afterward if necessary.
  assert(ldbase_difflist_sample_ids[ldbase_difflist_len] == sample_ct);
  uint32_t diff_from_ldbase_len;
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, diff_from_ldbase_raregeno_iter, &group_info_iter, &diff_from_ldbase_len);
  if (reterr) {
    return reterr;
  }
  if (!diff_from_ldbase_len) {
    memcpy(merged_difflist_sample_ids, ldbase_difflist_sample_ids, ldbase_difflist_len * sizeof(int32_t));
    *ambig_id_ct_ptr = 0;
    *merged_difflist_len_ptr = ldbase_difflist_len;
    CopyQuaterarr(ldbase_raregeno, ldbase_difflist_len, merged_raregeno);
    return kPglRetSuccess;
  }
  if (ldbase_common_geno == 3) {
    ldbase_common_geno = 4;  // force these to be saved
  }
  const uint32_t subgroup_idx_last = (diff_from_ldbase_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  uintptr_t* merged_raregeno_iter = merged_raregeno;
  uint32_t* merged_difflist_sample_ids_iter = merged_difflist_sample_ids;
  uintptr_t merged_raregeno_word = 0;
  uintptr_t ldbase_raregeno_word = 0;
  uintptr_t diff_from_ldbase_raregeno_word = 0;
  uint32_t ldbase_sample_idx = ldbase_difflist_sample_ids[0];
  uintptr_t raw_sample_idx = 0;
  uintptr_t cur_geno = 0;
  uint32_t sample_idx = 0;
  uint32_t ldbase_difflist_idx = 0;
  uint32_t ambig_id_ct = 0;
  uint32_t done = 0;
  uint32_t subgroup_idx = 0;
  uint32_t subgroup_len_m1 = kBitsPerWordD2 - 1;
  uint32_t merge_idx_lowbits = 0;
  while (1) {
    uint32_t diff_from_ldbase_idx_lowbits = 0;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        done = 1;
        sample_idx = sample_ct;
        goto ParseLdAndMergeDifflistSubset_finish;
      }
      subgroup_len_m1 &= diff_from_ldbase_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    diff_from_ldbase_raregeno_word = *diff_from_ldbase_raregeno_iter++;
    ++subgroup_idx;
    while (1) {
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
      cur_geno = diff_from_ldbase_raregeno_word & 3;
      if ((!sample_include) || IsSet(sample_include, raw_sample_idx)) {
        sample_idx = sample_include? RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, raw_sample_idx) : raw_sample_idx;
      ParseLdAndMergeDifflistSubset_finish:
        while (ldbase_sample_idx < sample_idx) {
          // replace with blocked copy?
          if (!(ldbase_difflist_idx % kBitsPerWordD2)) {
            ldbase_raregeno_word = ldbase_raregeno[ldbase_difflist_idx / kBitsPerWordD2];
          }
          *merged_difflist_sample_ids_iter++ = ldbase_sample_idx;
          merged_raregeno_word |= (ldbase_raregeno_word & 3) << (2 * merge_idx_lowbits);
          if (merge_idx_lowbits++ == (kBitsPerWordD2 - 1)) {
            *merged_raregeno_iter++ = merged_raregeno_word;
            merged_raregeno_word = 0;
            merge_idx_lowbits = 0;
          }
          ++ldbase_difflist_idx;
          ldbase_raregeno_word >>= 2;
          ldbase_sample_idx = ldbase_difflist_sample_ids[ldbase_difflist_idx];
        }
        if (ldbase_sample_idx == sample_idx) {
          if (done) {
            if (merge_idx_lowbits) {
              *merged_raregeno_iter = merged_raregeno_word;
            }
            *ambig_id_ct_ptr = ambig_id_ct;
            *merged_difflist_len_ptr = merged_difflist_sample_ids_iter - merged_difflist_sample_ids;
            return kPglRetSuccess;
          }
          if (!(ldbase_difflist_idx % kBitsPerWordD2)) {
            ldbase_raregeno_word = ldbase_raregeno[ldbase_difflist_idx / kBitsPerWordD2];
          }
          ++ldbase_difflist_idx;
          ldbase_raregeno_word >>= 2;
          ldbase_sample_idx = ldbase_difflist_sample_ids[ldbase_difflist_idx];
        }
        if (cur_geno != ldbase_common_geno) {
          *merged_difflist_sample_ids_iter++ = sample_idx;
          merged_raregeno_word |= cur_geno << (2 * merge_idx_lowbits);
          if (merge_idx_lowbits++ == (kBitsPerWordD2 - 1)) {
            *merged_raregeno_iter++ = merged_raregeno_word;
            merged_raregeno_word = 0;
            merge_idx_lowbits = 0;
          }
        }
      }
      if (ambig_sample_ids && (cur_geno == 3)) {
        ambig_sample_ids[ambig_id_ct++] = raw_sample_idx;
      }
      if (diff_from_ldbase_idx_lowbits == subgroup_len_m1) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      ++diff_from_ldbase_idx_lowbits;
      diff_from_ldbase_raregeno_word >>= 2;
    }
  }
}

void PgrDifflistToGenovecUnsafe(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uintptr_t difflist_common_geno, uint32_t sample_ct, uint32_t difflist_len, uintptr_t* __restrict genovec) {
  // Ok for trailing bits of raregeno to be nonzero.  Does not zero out
  // trailing bits of genovec.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  // could just memset up to word boundary; this should be a bit more
  // vector-instruction-friendly, though?
  memset(genovec, S_CAST(unsigned int, difflist_common_geno * 0x55), vec_ct * kBytesPerVec);
  const uintptr_t* raregeno_incr = raregeno;
  uint32_t difflist_idx = 0;
  uint32_t difflist_idx_stop = 0;
  if (!difflist_common_geno) {
    // faster inner loop since there's no existing value to mask out
    // todo: verify that memset with "unknown" parameter, set to zero, is only
    // a tiny bit slower than hardcoded memset zero
    // todo: check if this should just be deleted since the code bloat causes
    // too many more cache misses
    while (1) {
      difflist_idx_stop += kBitsPerWordD2;
      if (difflist_idx_stop > difflist_len) {
        if (difflist_idx == difflist_len) {
          return;
        }
        difflist_idx_stop = difflist_len;
      }
      uintptr_t raregeno_word = *raregeno_incr++;
      for (; difflist_idx < difflist_idx_stop; ++difflist_idx) {
        const uint32_t cur_sample_idx = difflist_sample_ids[difflist_idx];
        genovec[cur_sample_idx / kBitsPerWordD2] |= (raregeno_word & 3) << (2 * (cur_sample_idx % kBitsPerWordD2));
        raregeno_word >>= 2;
      }
    }
  }
  while (1) {
    difflist_idx_stop += kBitsPerWordD2;
    if (difflist_idx_stop > difflist_len) {
      if (difflist_idx == difflist_len) {
        return;
      }
      difflist_idx_stop = difflist_len;
    }
    uintptr_t raregeno_word = *raregeno_incr++;
    for (; difflist_idx < difflist_idx_stop; ++difflist_idx) {
      const uint32_t cur_sample_idx = difflist_sample_ids[difflist_idx];
      AssignQuaterarrEntry(cur_sample_idx, raregeno_word & 3, genovec);
      raregeno_word >>= 2;
    }
  }
}

/*
void PrunedDifflistToGenovecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t sample_ct, uint32_t difflist_common_geno, uint32_t difflist_len, uintptr_t* __restrict genovec) {
  // Designed to be used after genovec subsetting.  Assumes all difflist
  // entries are valid.  Ok for trailing bits of raregeno to be nonzero.  Does
  // not zero out trailing bits of genovec.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  memset(genovec, difflist_common_geno * 0x55, vec_ct * kBytesPerVec);
  if (!difflist_len) {
    return;
  }
  const uintptr_t* raregeno_incr = raregeno;
  const uint32_t* difflist_sample_ids_iter = difflist_sample_ids;
  const uint32_t* difflist_sample_ids_end = &(difflist_sample_ids[difflist_len]);
  // don't think there's a point to separating out the
  // difflist_common_geno == 0 case here, since the RawToSubsettedPos
  // operation is a bit expensive
  while (1) {
    // er, get rid of this undefined behavior if we uncomment this function
    const uint32_t* difflist_sample_ids_stop = &(difflist_sample_ids_iter[kBitsPerWordD2]);
    uintptr_t raregeno_word = *raregeno_incr++;
    if (difflist_sample_ids_stop > difflist_sample_ids_end) {
      if (difflist_sample_ids_iter == difflist_sample_ids_end) {
        return;
      }
      difflist_sample_ids_stop = difflist_sample_ids_end;
    }
    while (1) {
      const uint32_t cur_sample_idx = *difflist_sample_ids_iter;
      const uint32_t cur_subsetted_pos = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, cur_sample_idx);
      AssignQuaterarrEntry(cur_subsetted_pos, raregeno_word & 3, genovec);
      if (difflist_sample_ids_iter++ == difflist_sample_ids_stop) {
        break;
      }
      raregeno_word >>= 2;
    }
  }
}
*/

PglErr ParseAndApplyDifflist(const unsigned char* fread_end, uint32_t multiallelic_relevant, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effects: uses pgr.workspace_raregeno_tmp_loadbuf.
  // Cannot occur after genovec subsetting since the difflist sample indexes
  // will be incorrect.
  // If multiallelic_relevant is true, a list of sample indices with freshly
  // loaded raregeno value 0b11 is saved to pgr.workspace_ambig_sample_ids, and
  // pgr.workspace_ambig_id_ct is set to the length of the list.
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  uintptr_t* cur_raregeno_iter = pgrp->workspace_raregeno_tmp_loadbuf;
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, cur_raregeno_iter, &group_info_iter, &difflist_len);
  if (reterr || (!difflist_len)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  uint32_t* ambig_sample_ids = multiallelic_relevant? pgrp->workspace_ambig_sample_ids : nullptr;
  uintptr_t raw_sample_idx = 0;
  uint32_t ambig_id_ct = 0;
  uint32_t subgroup_idx = 0;
  while (1) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        pgrp->workspace_ambig_id_ct = ambig_id_ct;
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    ++subgroup_idx;
    uintptr_t cur_raregeno_word = *cur_raregeno_iter++;
    // This loop tends to be the decompression bottleneck.  Tried to modify it
    // to process 4 entries at a time, but that didn't end up helping.
    while (1) {
      // always check, since otherwise AssignQuaterarrEntry() can scribble
      // over arbitrary memory
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
      const uintptr_t cur_geno = cur_raregeno_word & 3;
      AssignQuaterarrEntry(raw_sample_idx, cur_geno, genovec);
      if (multiallelic_relevant && (cur_geno == 3)) {
        ambig_sample_ids[ambig_id_ct++] = raw_sample_idx;
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      --remaining_deltas_in_subgroup;
      cur_raregeno_word >>= 2;
    }
  }
}

// could merge ParseAndApplyDifflist() with this?
PglErr ParseAndApplyDifflistSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t multiallelic_relevant, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effects: uses pgr.workspace_raregeno_tmp_loadbuf.
  // If multiallelic_relevant is true, a list of sample indices with freshly
  // loaded raregeno value 0b11 is saved to pgr.workspace_ambig_sample_ids, and
  // pgr.workspace_ambig_id_ct is set to the length of the list.
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  if (sample_ct == raw_sample_ct) {
    return ParseAndApplyDifflist(fread_end, multiallelic_relevant, fread_pp, pgrp, genovec);
  }
  uintptr_t* cur_raregeno_iter = pgrp->workspace_raregeno_tmp_loadbuf;
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, cur_raregeno_iter, &group_info_iter, &difflist_len);
  if (reterr || (!difflist_len)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  uint32_t* ambig_sample_ids = multiallelic_relevant? pgrp->workspace_ambig_sample_ids : nullptr;
  uintptr_t raw_sample_idx = 0;
  uint32_t ambig_id_ct = 0;
  uint32_t subgroup_idx = 0;
  while (1) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        pgrp->workspace_ambig_id_ct = ambig_id_ct;
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    ++subgroup_idx;
    uintptr_t cur_raregeno_word = *cur_raregeno_iter++;
    // This loop tends to be the decompression bottleneck.  Tried to modify it
    // to process 4 entries at a time, but that didn't end up helping.
    while (1) {
      // always check, since otherwise AssignQuaterarrEntry() can scribble
      // over arbitrary memory
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
      const uintptr_t cur_geno = cur_raregeno_word & 3;
      if (IsSet(sample_include, raw_sample_idx)) {
        AssignQuaterarrEntry(RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, raw_sample_idx), cur_geno, genovec);
      }
      if (multiallelic_relevant && (cur_geno == 3)) {
        ambig_sample_ids[ambig_id_ct++] = raw_sample_idx;
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      --remaining_deltas_in_subgroup;
      cur_raregeno_word >>= 2;
    }
  }
}

PglErr ParseOnebitUnsafe(const unsigned char* fread_end, uint32_t difflist_ambig_ids_needed, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // doesn't zero out trailing genovec bits
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t common2_and_bitarray_byte_ct = (raw_sample_ct + 15) / CHAR_BIT;
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < common2_and_bitarray_byte_ct) {
    return kPglRetMalformedInput;
  }
  const unsigned char* fread_difflist_start = &((*fread_pp)[common2_and_bitarray_byte_ct]);
  const uintptr_t common2_code = *((*fread_pp)++);
  const uintptr_t word_base = (common2_code / 4) * kMask5555;
  const uintptr_t common_code_delta = common2_code & 3;
  const uint32_t genovec_widx_trail = (raw_sample_ct + 7) / kBitsPerWordD2;
  const uint32_t genovec_widx_end = QuaterCtToWordCt(raw_sample_ct);
  uint32_t genovec_widx = 0;
#ifdef __arm__
#  error "Unaligned accesses in ParseOnebitUnsafe()."
#endif
  const Halfword* fread_alias = R_CAST(const Halfword*, *fread_pp);
  while (1) {
    uintptr_t ww;
    if (genovec_widx >= genovec_widx_trail) {
      // might want to modify to not go here if last read is an entire halfword
      if (genovec_widx == genovec_widx_end) {
        break;
      }
      ww = ProperSubwordLoad(&(fread_alias[genovec_widx_trail]), 1 + (((raw_sample_ct - 1) % kBitsPerWordD2) / CHAR_BIT));
    } else {
      ww = fread_alias[genovec_widx];
    }
    // apply middle-out operation
    // 64-bit:
    //   const uintptr_t middle_out_result = (ww | (ww << 31)) & kMask5555;
    // 32-bit:
    //   *genovec_iter++ = word_base + (ww & kMask5555) * common_code_delta;
    //   *genovec_iter++ = word_base + ((ww >> 1) & kMask5555) * common_code_delta;
    // (scrapped for now since the time savings don't seem to be worth the
    // extra end-of-vector corner cases, apparently the extra operations here
    // are sufficiently cheap)
    ww = UnpackHalfwordToWord(ww);
    genovec[genovec_widx++] = word_base + ww * common_code_delta;
  }
  *fread_pp = fread_difflist_start;
  return ParseAndApplyDifflist(fread_end, difflist_ambig_ids_needed, fread_pp, pgrp, genovec);
}

void ExtractGenoarrAmbigIds(const uintptr_t* genoarr, uint32_t raw_sample_ct, uint32_t* __restrict ambig_sample_ids, uint32_t* ambig_id_ct_ptr) {
#ifdef __arm__
#  error "Unaligned accesses in ExtractGenoarrAmbigIds() (genoarr may not be aligned)."
#endif
  // does not read trailing bytes of genoarr
  const uint32_t word_ct_trail = (raw_sample_ct + 3) / kBitsPerWordD2;
  const uint32_t word_ct_end = QuaterCtToWordCt(raw_sample_ct);
  uint32_t ambig_id_ct = 0;
  uint32_t widx = 0;
  while (1) {
    uintptr_t detect_11;
    if (widx >= word_ct_trail) {
      if (widx == word_ct_end) {
        *ambig_id_ct_ptr = ambig_id_ct;
        return;
      }
      detect_11 = SubwordLoad(&(genoarr[widx]), QuaterCtToByteCt(raw_sample_ct % kBitsPerWordD2));
    } else {
      detect_11 = genoarr[widx];
    }
    detect_11 = detect_11 & (detect_11 >> 1) & kMask5555;
    // now detect_11 has a set bit iff the genoarr entry is 0b11
    if (detect_11) {
      const uint32_t sample_idx_base = widx * kBitsPerWordD2;
      do {
        const uint32_t sample_idx_lowbits = CTZLU(detect_11) / 2;
        ambig_sample_ids[ambig_id_ct++] = sample_idx_base + sample_idx_lowbits;
        detect_11 &= detect_11 - 1;
      } while (detect_11);
    }
  }
}

PglErr Parse1or2bitGenovecUnsafe(const unsigned char* fread_end, uint32_t vrtype, uint32_t difflist_ambig_ids_needed, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effect: may use pgrp->workspace_raregeno_tmp_loadbuf.
  // Does not update fp_vidx, does not rotate plink1-formatted data (since it's
  // better to do that post-subsetting)
  // If difflist_ambig_ids_needed is set, pgr.workspace_ambig_sample_ids and
  // pgr.workspace_ambig_id_ct are filled with with the sample IDs
  // corresponding to aux track 1.
  assert((!difflist_ambig_ids_needed) || VrtypeMultiallelic(vrtype));
  if (vrtype & 3) {
    return ParseOnebitUnsafe(fread_end, difflist_ambig_ids_needed, fread_pp, pgrp, genovec);
  }
  // uncompressed storage
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t genovec_byte_ct = QuaterCtToByteCt(raw_sample_ct);
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < genovec_byte_ct) {
    return kPglRetMalformedInput;
  }
  const unsigned char* new_fread_ptr = &((*fread_pp)[genovec_byte_ct]);
  memcpy(genovec, *fread_pp, genovec_byte_ct);
  *fread_pp = new_fread_ptr;
  if (difflist_ambig_ids_needed) {
    ExtractGenoarrAmbigIds(genovec, raw_sample_ct, pgrp->workspace_ambig_sample_ids, &(pgrp->workspace_ambig_id_ct));
  }
  return kPglRetSuccess;
}

PglErr ParseNonLdGenovecSubsetUnsafe(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vrtype, uint32_t difflist_ambig_ids_needed, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effects:
  //   may use pgrp->workspace_raregeno_tmp_loadbuf
  //   may use pgrp->workspace_vec (subsetting)
  // See comments on Parse1or2bitGenovecUnsafe().
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  if (!VrtypeDifflist(vrtype)) {
    const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
    uintptr_t* raw_genovec = subsetting_required? pgrp->workspace_vec : genovec;
    PglErr reterr = Parse1or2bitGenovecUnsafe(fread_end, vrtype, difflist_ambig_ids_needed, fread_pp, pgrp, raw_genovec);
    if ((!subsetting_required) || reterr) {
      return reterr;
    }
    CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, genovec);
    return kPglRetSuccess;
  }
  assert((!difflist_ambig_ids_needed) || VrtypeMultiallelic(vrtype));
  const uint32_t vrtype_low2 = vrtype & 3;
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);

  // This memset is frequently the limiting operation.  Probably not worth
  // worrying about before 2019 or so, but it suggests that the
  // _difflist_or_genovec interface shouldn't be abandoned.
  memset(genovec, vrtype_low2 * 0x55, vec_ct * kBytesPerVec);
  return ParseAndApplyDifflistSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, difflist_ambig_ids_needed, fread_pp, pgrp, genovec);
}

PglErr InitReadPtrs(uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp) {
  const unsigned char* block_base = pgrp->fi.block_base;
  if (block_base != nullptr) {
    // possible todo: special handling of end of vblock
    const uint64_t block_offset = pgrp->fi.block_offset;
    *fread_pp = &(block_base[GetPgfiFpos(&(pgrp->fi), vidx) - block_offset]);
    *fread_endp = &(block_base[GetPgfiFpos(&(pgrp->fi), vidx + 1) - block_offset]);
    pgrp->fp_vidx = vidx + 1;
    return kPglRetSuccess;
  }
  if (pgrp->fp_vidx != vidx) {
    if (fseeko(pgrp->ff, GetPgfiFpos(&(pgrp->fi), vidx), SEEK_SET)) {
      return kPglRetReadFail;
    }
  }
  const uintptr_t cur_vrec_width = GetPgfiVrecWidth(&(pgrp->fi), vidx);
#ifdef __LP64__
  if (fread_checked(pgrp->fread_buf, cur_vrec_width, pgrp->ff)) {
    return kPglRetReadFail;
  }
#else
  // cur_vrec_width < 2^31 since otherwise we error out on initialization
  if (!fread_unlocked(pgrp->fread_buf, cur_vrec_width, 1, pgrp->ff)) {
    return kPglRetReadFail;
  }
#endif
  *fread_pp = pgrp->fread_buf;
  *fread_endp = &(pgrp->fread_buf[cur_vrec_width]);
  pgrp->fp_vidx = vidx + 1;
  return kPglRetSuccess;
}

uint32_t LdLoadNecessary(uint32_t cur_vidx, PgenReader* pgrp) {
  // Determines whether LD base variant needs to be loaded (in addition to the
  // current variant).  Updates pgrp->ldbase_vidx when necessary.
  // possible todo: add the vidx_word_stop optimization to GetLdbaseVidx() so
  // this can call it instead of duplicating code (this forces
  // pgfi_block_load() to pass a fp_vidx parameter of zero).
  const uint32_t fp_vidx = pgrp->ldbase_stypes? pgrp->fp_vidx : 0;
  if (cur_vidx == fp_vidx) {
    // ldbase variant guaranteed to be up-to-date if we didn't skip the last
    // variant, and cache wasn't cleared
    return 0;
  }
  // Find the last vrtypes[] value before vrtypes[cur_vidx] with bit 1 unset or
  // bit 2 set.
  const uintptr_t* vrtypes_walias = R_CAST(const uintptr_t*, pgrp->fi.vrtypes);
  const uint32_t cur_vidx_orig_remainder = cur_vidx % kBytesPerWord;
  uint32_t vidx_word_idx = (cur_vidx - 1) / kBytesPerWord;
  uintptr_t cur_vrtypes_word = vrtypes_walias[vidx_word_idx];
  if (cur_vidx_orig_remainder) {
    // make sure we don't detect a byte after the current position.
    cur_vrtypes_word = bzhi(cur_vrtypes_word, CHAR_BIT * cur_vidx_orig_remainder);
    cur_vrtypes_word |= (kMask0101 * 2) << (CHAR_BIT * cur_vidx_orig_remainder);
  }
  const uint32_t vidx_word_stop = (fp_vidx < cur_vidx)? (fp_vidx / kBytesPerWord) : 0;
  while (1) {
    // ((bit 2) OR (NOT bit 1)) for each byte.  (possible experiment: see if
    // the same assembly is generated if this expression is rewritten to use
    // ands/nots.)
    uintptr_t detect_non_ld_word = ((cur_vrtypes_word >> 1) | (~cur_vrtypes_word)) & (kMask0101 * 2);

    if (detect_non_ld_word) {
      // find the highest-order set bit in detect_non_ld_word; this corresponds
      // to the last non-LD-compressed byte (assuming little-endian).
      const uint32_t old_ldbase_vidx = pgrp->ldbase_vidx;
      const uint32_t new_ldbase_vidx_loworder = kBytesPerWord - 1 - (CLZLU(detect_non_ld_word) / CHAR_BIT);
      const uint32_t new_ldbase_vidx = (vidx_word_idx * kBytesPerWord) + new_ldbase_vidx_loworder;
      pgrp->ldbase_vidx = new_ldbase_vidx;
      return (old_ldbase_vidx != new_ldbase_vidx);
    }
    // everything LD-compressed in the current block.  move back 8 bytes in the
    // array (or 4-bytes for 32-bit build).
    if (vidx_word_idx == vidx_word_stop) {
      return 0;
    }
    --vidx_word_idx;
    cur_vrtypes_word = vrtypes_walias[vidx_word_idx];
  }
}

// loads ldbase variant if necessary, unpacks to genovec if necessary (the
// latter happens even if the variant was already loaded)
// may use workspace_vec
PglErr LdLoadGenovecSubsetIfNecessary(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp) {
  if (LdLoadNecessary(vidx, pgrp)) {
    const uint32_t ldbase_vidx = pgrp->ldbase_vidx;
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (InitReadPtrs(ldbase_vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    pgrp->ldbase_stypes = kfPgrLdcacheQuater;
    return ParseNonLdGenovecSubsetUnsafe(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, pgrp->fi.vrtypes[ldbase_vidx], 0, &fread_ptr, pgrp, pgrp->ldbase_genovec);
  }
  if (!(pgrp->ldbase_stypes & kfPgrLdcacheQuater)) {
    assert(pgrp->ldbase_stypes & kfPgrLdcacheDifflist);
    PgrDifflistToGenovecUnsafe(pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, pgrp->fi.vrtypes[pgrp->ldbase_vidx] & 3, sample_ct, pgrp->ldbase_difflist_len, pgrp->ldbase_genovec);
    pgrp->ldbase_stypes |= kfPgrLdcacheQuater;
  }
  return kPglRetSuccess;
}

// fread_pp should be non-null iff this is being called by an internal function
// as part of a multiallelic variant read
PglErr ReadRefalt1GenovecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict genovec) {
  // Side effects:
  //   may use pgr.workspace_vec iff subsetting required
  //   may use pgr.workspace_raregeno_tmp_loadbuf (any difflist)
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t maintrack_vrtype = vrtype & 7;
  const uint32_t multiallelic_relevant = fread_pp && VrtypeMultiallelic(vrtype);
  if (VrtypeLdCompressed(maintrack_vrtype)) {
    // LD compression
    PglErr reterr = LdLoadGenovecSubsetIfNecessary(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp);
    if (reterr) {
      return reterr;
    }
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    CopyQuaterarr(pgrp->ldbase_genovec, sample_ct, genovec);
    reterr = ParseAndApplyDifflistSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, multiallelic_relevant, &fread_ptr, pgrp, genovec);
    if (reterr) {
      return reterr;
    }
    if (maintrack_vrtype == 3) {
      GenovecInvertUnsafe(sample_ct, genovec);
    }
    if (fread_pp) {
      *fread_pp = fread_ptr;
      *fread_endp = fread_end;
    }
    return kPglRetSuccess;
  }
  const unsigned char* fread_ptr;
  const unsigned char* fread_end = nullptr;  // maybe-uninitialized warning
  // tried inserting special-case code for the plink1 case to avoid a copy, and
  // it was actually slower
  if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
    return kPglRetReadFail;
  }
  // tried to add more sophisticated caching, but turns out it isn't worth it
  PglErr reterr = ParseNonLdGenovecSubsetUnsafe(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, maintrack_vrtype, multiallelic_relevant, &fread_ptr, pgrp, genovec);
  if (reterr) {
    return reterr;
  }
  const uint32_t is_ldbase = pgrp->fi.vrtypes && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
  if (is_ldbase) {
    CopyQuaterarr(genovec, sample_ct, pgrp->ldbase_genovec);
    pgrp->ldbase_vidx = vidx;
    pgrp->ldbase_stypes = kfPgrLdcacheQuater;
  }
  if (vrtype == kPglVrtypePlink1) {
    PgrPlink1ToPlink2InplaceUnsafe(sample_ct, genovec);
  } else if (fread_pp) {
    *fread_pp = fread_ptr;
    *fread_endp = fread_end;
  }
  return kPglRetSuccess;
}

PglErr PgrReadRefalt1GenovecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  return ReadRefalt1GenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genovec);
}

/*
void CopyAndSubsetDifflist(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict raw_raregeno, const uint32_t* __restrict raw_difflist_sample_ids, uint32_t raw_difflist_len, uintptr_t* __restrict new_raregeno, uint32_t* __restrict new_difflist_sample_ids, uint32_t* __restrict new_difflist_len_ptr) {
  // Trailing bits of new_raregeno are zeroed out.
  if (!raw_difflist_len) {
    *new_difflist_len_ptr = 0;
    return;
  }
  const uintptr_t* raw_raregeno_incr = raw_raregeno;
  const uint32_t* raw_difflist_sample_ids_iter = raw_difflist_sample_ids;
  const uint32_t* raw_difflist_sample_ids_last = &(raw_difflist_sample_ids[RoundDownPow2(raw_difflist_len - 1, kBitsPerWordD2)]);
  uintptr_t* new_raregeno_incr = new_raregeno;
  uintptr_t new_raregeno_word = 0;
  uint32_t new_difflist_len = 0;
  uint32_t block_len_m1 = kBitsPerWordD2 - 1;
  while (1) {
    if (raw_difflist_sample_ids_iter >= raw_difflist_sample_ids_last) {
      if (raw_difflist_sample_ids_iter > raw_difflist_sample_ids_last) {
        if (new_difflist_len % kBitsPerWordD2) {
          *new_raregeno_incr = new_raregeno_word;
        }
        *new_difflist_len_ptr = new_difflist_len;
        return;
      }
      block_len_m1 &= raw_difflist_len - 1;
    }
    uintptr_t raw_raregeno_word = *raw_raregeno_incr++;
    uint32_t raw_difflist_idx_lowbits = 0;
    while (1) {
      const uint32_t raw_sample_idx = raw_difflist_sample_ids_iter[raw_difflist_idx_lowbits];
      if (IsSet(sample_include, raw_sample_idx)) {
        new_difflist_sample_ids[new_difflist_len] = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, raw_sample_idx);
        new_raregeno_word |= ((raw_raregeno_word >> (2 * raw_difflist_idx_lowbits)) & 3) << (2 * (new_difflist_len % kBitsPerWordD2));
        ++new_difflist_len;
        if (!(new_difflist_len % kBitsPerWordD2)) {
          *new_raregeno_incr++ = new_raregeno_word;
          new_raregeno_word = 0;
        }
      }
      if (raw_difflist_idx_lowbits == block_len_m1) {
        break;
      }
      ++raw_difflist_idx_lowbits;
    }
    raw_difflist_sample_ids_iter = &(raw_difflist_sample_ids_iter[kBitsPerWordD2]);
  }
}
*/

// populates pgrp->ldbase_genovec or
// pgrp->ldbase_{raregeno,difflist_sample_ids,difflist_len}, depending on
// storage type
// requires workspace_vec
PglErr LdLoadMinimalSubsetIfNecessary(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp) {
  if (!LdLoadNecessary(vidx, pgrp)) {
    return kPglRetSuccess;
  }
  const uint32_t ldbase_vidx = pgrp->ldbase_vidx;
  const uint64_t cur_vidx_fpos = pgrp->fi.var_fpos[ldbase_vidx];
  const uint32_t ldbase_vrtype = pgrp->fi.vrtypes[ldbase_vidx];
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* raw_genovec = subsetting_required? pgrp->workspace_vec : pgrp->ldbase_genovec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  const unsigned char* block_base = pgrp->fi.block_base;
  PglErr reterr = kPglRetSuccess;
  if (block_base != nullptr) {
    {
      const uint64_t block_offset = pgrp->fi.block_offset;
      fread_ptr = &(block_base[cur_vidx_fpos - block_offset]);
      fread_end = &(block_base[pgrp->fi.var_fpos[ldbase_vidx + 1] - block_offset]);
    }
    if (!VrtypeDifflist(ldbase_vrtype)) {
      pgrp->ldbase_stypes = kfPgrLdcacheQuater;
      reterr = Parse1or2bitGenovecUnsafe(fread_end, ldbase_vrtype, 0, &fread_ptr, pgrp, raw_genovec);
    LdLoadMinimalSubsetIfNecessary_genovec_finish:
      if ((!subsetting_required) || reterr) {
        return reterr;
      }
      CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, pgrp->ldbase_genovec);
      return kPglRetSuccess;
    }
  } else {
    if (fseeko(pgrp->ff, pgrp->fi.var_fpos[ldbase_vidx], SEEK_SET)) {
      return kPglRetReadFail;
    }
    pgrp->ldbase_stypes = kfPgrLdcacheQuater;
    if (!(ldbase_vrtype & 7)) {
      // don't actually need to fread the whole record in this case
      if (!fread_unlocked(raw_genovec, QuaterCtToByteCt(raw_sample_ct), 1, pgrp->ff)) {
        return kPglRetReadFail;
      }
      goto LdLoadMinimalSubsetIfNecessary_genovec_finish;
    }
    const uintptr_t cur_vrec_width = pgrp->fi.var_fpos[ldbase_vidx + 1] - cur_vidx_fpos;
    if (!fread_unlocked(pgrp->fread_buf, cur_vrec_width, 1, pgrp->ff)) {
      return kPglRetReadFail;
    }
    fread_ptr = pgrp->fread_buf;
    fread_end = &(pgrp->fread_buf[cur_vrec_width]);
    if (!VrtypeDifflist(ldbase_vrtype)) {
      reterr = ParseOnebitUnsafe(fread_end, 0, &fread_ptr, pgrp, raw_genovec);
      goto LdLoadMinimalSubsetIfNecessary_genovec_finish;
    }
  }
  uint32_t ldbase_difflist_len;
  if (!subsetting_required) {
    reterr = ParseAndSaveDifflist(fread_end, raw_sample_ct, &fread_ptr, pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, &ldbase_difflist_len);
  } else {
    reterr = ParseAndSaveDifflistProperSubset(fread_end, sample_include, sample_include_cumulative_popcounts, raw_sample_ct, &fread_ptr, pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, &ldbase_difflist_len, pgrp->workspace_raregeno_tmp_loadbuf);
  }
  if (reterr) {
    return reterr;
  }
  pgrp->ldbase_difflist_len = ldbase_difflist_len;
  pgrp->ldbase_difflist_sample_ids[ldbase_difflist_len] = sample_ct;
  pgrp->ldbase_stypes = kfPgrLdcacheDifflist;
  pgrp->fp_vidx = ldbase_vidx + 1;
  return kPglRetSuccess;
}

PglErr ReadRefalt1DifflistOrGenovecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t max_simple_difflist_len, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict genovec, uint32_t* difflist_common_geno_ptr, uintptr_t* __restrict main_raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  assert(sample_ct);
  // Side effects:
  //   may use pgr.workspace_vec or workspace_difflist_sample_ids_tmp iff
  //     subsetting required.
  //   may use pgr.workspace_raregeno_tmp_loadbuf
  // Trailing bits of genovec/main_raregeno may not be zeroed out.
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t maintrack_vrtype = vrtype & 7;
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  const uint32_t multiallelic_relevant = fread_pp && VrtypeMultiallelic(vrtype);
  if (VrtypeLdCompressed(maintrack_vrtype)) {
    // LD compression

    // note that this can currently load a difflist longer than
    // max_simple_difflist_len
    PglErr reterr = LdLoadMinimalSubsetIfNecessary(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp);
    if (reterr) {
      return reterr;
    }
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    const uint32_t ld_invert = (maintrack_vrtype == 3);
    if (pgrp->ldbase_stypes & kfPgrLdcacheDifflist) {
      const uint32_t ldbase_common_geno = pgrp->fi.vrtypes[pgrp->ldbase_vidx] & 3;
      // unnecessary for this to branch on LD difflist length, since that's
      // limited to 3/4 of the ldbase difflist length.
      *difflist_common_geno_ptr = ldbase_common_geno;
      reterr = ParseLdAndMergeDifflistSubset(fread_end, subsetting_required? sample_include : nullptr, sample_include_cumulative_popcounts, pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, pgrp->ldbase_difflist_len, ldbase_common_geno, raw_sample_ct, sample_ct, &fread_ptr, multiallelic_relevant? pgrp->workspace_ambig_sample_ids : nullptr, main_raregeno, difflist_sample_ids, difflist_len_ptr, &(pgrp->workspace_ambig_id_ct), pgrp->workspace_raregeno_tmp_loadbuf);
      if (reterr) {
        return reterr;
      }
      if (ld_invert) {
        *difflist_common_geno_ptr = (6 - ldbase_common_geno) & 3;
        GenovecInvertUnsafe(*difflist_len_ptr, main_raregeno);
      }
      return kPglRetSuccess;
    }
    assert(pgrp->ldbase_stypes & kfPgrLdcacheQuater);
    *difflist_common_geno_ptr = UINT32_MAX;
    CopyQuaterarr(pgrp->ldbase_genovec, sample_ct, genovec);
    reterr = ParseAndApplyDifflistSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, multiallelic_relevant, &fread_ptr, pgrp, genovec);
    if (reterr) {
      return reterr;
    }
    if (ld_invert) {
      GenovecInvertUnsafe(sample_ct, genovec);
    }
    if (fread_pp) {
      *fread_pp = fread_ptr;
      *fread_endp = fread_end;
    }
    return kPglRetSuccess;
  }
  const unsigned char* fread_ptr;
  const unsigned char* fread_end = nullptr;  // maybe-uninitialized warning
  if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
    return kPglRetReadFail;
  }
  const uint32_t is_ldbase = pgrp->fi.vrtypes && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
  const uint32_t saved_difflist_len = VrtypeDifflist(vrtype)? PeekVint31(fread_ptr, fread_end) : raw_sample_ct;
  pgrp->ldbase_vidx = vidx;
  // no limit is slightly better than /16 but substantially worse than /32 on
  // the large test dataset (/64 is slightly worse than /32)
  // no limit is best on the small test dataset
  if (saved_difflist_len > max_simple_difflist_len) {
    *difflist_common_geno_ptr = UINT32_MAX;
    PglErr reterr = ParseNonLdGenovecSubsetUnsafe(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, vrtype, multiallelic_relevant, &fread_ptr, pgrp, genovec);
    if (reterr) {
      return reterr;
    }
    if (is_ldbase) {
      CopyQuaterarr(genovec, sample_ct, pgrp->ldbase_genovec);
      pgrp->ldbase_stypes = kfPgrLdcacheQuater;
    }
    if (vrtype == kPglVrtypePlink1) {
      PgrPlink1ToPlink2InplaceUnsafe(sample_ct, genovec);
    }
    if (fread_pp) {
      *fread_pp = fread_ptr;
      *fread_endp = fread_end;
    }
    return kPglRetSuccess;
  }
  *difflist_common_geno_ptr = vrtype & 3;
  if (ParseAndSaveDifflistSubset(fread_end, subsetting_required? sample_include : nullptr, sample_include_cumulative_popcounts, raw_sample_ct, &fread_ptr, multiallelic_relevant? pgrp->workspace_ambig_sample_ids : nullptr, main_raregeno, difflist_sample_ids, difflist_len_ptr, &(pgrp->workspace_ambig_id_ct), pgrp->workspace_raregeno_tmp_loadbuf)) {
    return kPglRetMalformedInput;
  }
  if (is_ldbase) {
    const uint32_t difflist_len = *difflist_len_ptr;
    pgrp->ldbase_stypes = kfPgrLdcacheDifflist;
    pgrp->ldbase_difflist_len = difflist_len;
    CopyQuaterarr(main_raregeno, difflist_len, pgrp->ldbase_raregeno);
    memcpy(pgrp->ldbase_difflist_sample_ids, difflist_sample_ids, difflist_len * sizeof(int32_t));
    pgrp->ldbase_difflist_sample_ids[difflist_len] = sample_ct;
  }
  if (fread_pp) {
    *fread_pp = fread_ptr;
    *fread_endp = fread_end;
  }
  return kPglRetSuccess;
}

PglErr PgrReadRefalt1DifflistOrGenovecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t max_simple_difflist_len, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uint32_t* difflist_common_geno_ptr, uintptr_t* __restrict main_raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    *difflist_common_geno_ptr = UINT32_MAX;
    return kPglRetSuccess;
  }
  return ReadRefalt1DifflistOrGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, max_simple_difflist_len, vidx, pgrp, nullptr, nullptr, genovec, difflist_common_geno_ptr, main_raregeno, difflist_sample_ids, difflist_len_ptr);
}

PglErr LdSubsetAdjustGenocounts(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict ldbase_genovec, uint32_t raw_sample_ct, const unsigned char** fread_pp, uint32_t* __restrict genocounts, uint32_t* __restrict ambig_sample_ids, uint32_t* __restrict ambig_id_ct_ptr, uint32_t* __restrict ambig_id_ct_filtered_ptr, uintptr_t* __restrict raregeno_workspace) {
  // * Assumes genocounts[] is initialized to the proper values for the LD
  //   reference variant (including subsetting).
  // * Tried a hybrid implementation which allowed the base variant to be saved
  //   as a difflist; turns out it's practically always better to unpack to a
  //   genovec first.
  // * ambig_sample_ids should be nullptr if it doesn't need to be filled.
  //   Note that, for a multiallelic variant, we don't need ambig_sample_ids to
  //   be filled unless we're also looking at a subset of the samples.  If we
  //   skip filling ambig_sample_ids, ambig_id_ct will be zero, but
  //   ambig_id_ct_filtered will contain the correct value.  (Strangely, the
  //   function slows down substantially if I try to conditionally assign it to
  //   ambig_id_ct instead of returning it separately.)
  // * This is the main frequency-counting bottleneck.
  uint32_t raw_difflist_len;
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno_workspace, &group_info_iter, &raw_difflist_len);
  if (reterr || (!raw_difflist_len)) {
    *ambig_id_ct_ptr = 0;
    // assumes ambig_id_ct_filtered is initialized to zero
    return reterr;
  }
  const uint32_t subgroup_idx_last = (raw_difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  uintptr_t* raregeno_workspace_iter = raregeno_workspace;
  uintptr_t raw_sample_idx = 0;
  uint32_t subgroup_idx = 0;
  uint32_t ambig_id_ct = 0;
  uint32_t delta_counts[16];
  ZeroUiArr(16, delta_counts);
  while (1) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        *ambig_id_ct_ptr = ambig_id_ct;
        *ambig_id_ct_filtered_ptr = delta_counts[12] + delta_counts[13] + delta_counts[14] + delta_counts[15];
        const int32_t incr0 = delta_counts[1] + delta_counts[2] + delta_counts[3] - delta_counts[4] - delta_counts[8] - delta_counts[12];
        const int32_t incr1 = delta_counts[4] + delta_counts[6] + delta_counts[7] - delta_counts[1] - delta_counts[9] - delta_counts[13];
        const int32_t incr2 = delta_counts[8] + delta_counts[9] + delta_counts[11] - delta_counts[2] - delta_counts[6] - delta_counts[14];
        genocounts[0] += incr0;
        genocounts[1] += incr1;
        genocounts[2] += incr2;
        genocounts[3] -= incr0 + incr1 + incr2;
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= raw_difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
#ifdef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    ++subgroup_idx;
    uintptr_t cur_raregeno_word = *raregeno_workspace_iter++;
    while (1) {
#ifndef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      const uintptr_t cur_geno = cur_raregeno_word & 3;
      if (!sample_include) {
        delta_counts[cur_geno * 4 + GetQuaterarrEntry(ldbase_genovec, raw_sample_idx)] += 1;
      } else if (IsSet(sample_include, raw_sample_idx)) {
        const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, raw_sample_idx);
        delta_counts[cur_geno * 4 + GetQuaterarrEntry(ldbase_genovec, sample_idx)] += 1;
      }
      if (ambig_sample_ids && (cur_geno == 3)) {
        ambig_sample_ids[ambig_id_ct++] = raw_sample_idx;
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      --remaining_deltas_in_subgroup;
      cur_raregeno_word >>= 2;
    }
  }
}

uint32_t BytesumArr(const void* bytearr, uint32_t byte_ct) {
  // Assumes sum < 2^32.
  // This is only slightly slower than SSE2 code, while tolerating an unaligned
  // starting address.
#ifdef __arm__
#  error "Unaligned accesses in BytesumArr()."
#endif
  const uint32_t word_ct = byte_ct / kBytesPerWord;
  const uintptr_t* bytearr_alias_iter = S_CAST(const uintptr_t*, bytearr);
  const uint32_t wordblock_idx_trail = word_ct / 256;
  const uint32_t wordblock_idx_end = DivUp(word_ct, 256);
  uint32_t wordblock_idx = 0;
  uint32_t wordblock_len = 256;
  uint32_t tot = 0;
  while (1) {
    if (wordblock_idx >= wordblock_idx_trail) {
      if (wordblock_idx == wordblock_idx_end) {
        byte_ct = byte_ct % kBytesPerWord;
        const unsigned char* bytearr_alias_iter2 = R_CAST(const unsigned char*, bytearr_alias_iter);
        for (uint32_t uii = 0; uii < byte_ct; ++uii) {
          tot += bytearr_alias_iter2[uii];
        }
        return tot;
      }
      wordblock_len = word_ct % 256;
    }
    ++wordblock_idx;
    const uintptr_t* bytearr_alias_stop = &(bytearr_alias_iter[wordblock_len]);
    uintptr_t acc_even = 0;
    uintptr_t acc_odd = 0;
    do {
      uintptr_t cur_word = *bytearr_alias_iter++;
      acc_even += cur_word & kMask00FF;
      acc_odd += (cur_word >> 8) & kMask00FF;
    } while (bytearr_alias_iter < bytearr_alias_stop);
    acc_even += acc_odd;
#ifdef __LP64__
    acc_even = (acc_even & kMask0000FFFF) + ((acc_even >> 16) & kMask0000FFFF);
#endif
    tot += S_CAST(Halfword, acc_even) + (acc_even >> kBitsPerWordD2);
  }
}

PglErr SkipDifflistIds(const unsigned char* fread_end, const unsigned char* group_info, uint32_t difflist_len, uint32_t raw_sample_ct, const unsigned char** fread_pp) {
  assert(difflist_len);
  // fread_pp is a pure output parameter here
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  const unsigned char* extra_byte_cts = &(group_info[group_ct * sample_id_byte_ct]);
  const uint32_t extra_byte_tot = BytesumArr(extra_byte_cts, group_ct - 1);

  // (group_ct - 1) for extra_byte_cts
  // (difflist_len + 3) / 4 for raregeno
  // (group_ct - 1) * (kPglDifflistGroupSize - 1) + extra_byte_tot for
  //   all but last ID block
  // total = (group_ct - 1) * kPglDifflistGroupSize + extra_byte_tot +
  //         (difflist_len + 3) / 4
#ifdef __arm__
#  error "Unaligned accesses in SkipDifflistIds()."
#endif
  const uintptr_t* fread_alias = R_CAST(const uintptr_t*, &(extra_byte_cts[(group_ct - 1) * kPglDifflistGroupSize + extra_byte_tot + QuaterCtToByteCt(difflist_len)]));
  const uintptr_t* fread_alias_stop = R_CAST(const uintptr_t*, &(fread_end[-S_CAST(int32_t, kBytesPerWord)]));
  uint32_t remaining_id_ct = (difflist_len - 1) % kPglDifflistGroupSize;
  while (remaining_id_ct >= kBytesPerWord) {
    // scan a word at a time, count number of high bits set
    if (fread_alias > fread_alias_stop) {
      return kPglRetMalformedInput;
    }
#ifdef USE_SSE42
    const uintptr_t ww = (*fread_alias++) & (0x80 * kMask0101);
    remaining_id_ct -= kBytesPerWord - PopcountWord(ww);
#else
    const uintptr_t ww = ((*fread_alias++) >> 7) & kMask0101;
    remaining_id_ct -= kBytesPerWord - ((ww * kMask0101) >> (kBitsPerWord - 8));
#endif
  }
  const unsigned char* fread_ptr = R_CAST(const unsigned char*, fread_alias);
  if (!remaining_id_ct) {
    *fread_pp = fread_ptr;
    return kPglRetSuccess;
  }
  --remaining_id_ct;
  while (fread_ptr < fread_end) {
    if ((*fread_ptr++) <= 127) {
      if (!remaining_id_ct) {
        *fread_pp = fread_ptr;
        return kPglRetSuccess;
      }
      --remaining_id_ct;
    }
  }
  return kPglRetMalformedInput;
}

PglErr CountparseDifflistSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, uint32_t common_geno, uint32_t raw_sample_ct, uint32_t sample_ct, const unsigned char** fread_pp, uint32_t* __restrict ambig_sample_ids, uint32_t* __restrict ambig_id_ct_ptr, uint32_t* __restrict genocounts, uintptr_t* __restrict raregeno_workspace) {
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno_workspace, &group_info_iter, &difflist_len);
  *ambig_id_ct_ptr = 0;
  ZeroUiArr(4, genocounts);
  if (reterr || (!difflist_len)) {
    genocounts[common_geno] = sample_ct;
    return reterr;
  }
  if (raw_sample_ct == sample_ct) {
    ZeroTrailingQuaters(difflist_len, raregeno_workspace);
    GenovecCountFreqsUnsafe(raregeno_workspace, difflist_len, genocounts);
    if (ambig_sample_ids && genocounts[3]) {
      // no need to update ambig_sample_ids[], but necessary to set ambig_id_ct
      // and fread_pp to enable rarealt counting.
      reterr = SkipDifflistIds(fread_end, group_info_iter, difflist_len, raw_sample_ct, fread_pp);
      *ambig_id_ct_ptr = genocounts[3];
    }
    genocounts[common_geno] += sample_ct - difflist_len;
    return kPglRetSuccess;
  }
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  uintptr_t* raregeno_workspace_iter = raregeno_workspace;
  uintptr_t raw_sample_idx = 0;
  uint32_t subgroup_idx = 0;
  uint32_t ambig_id_ct = 0;
  uint32_t common_decr = 0;
  while (1) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        *ambig_id_ct_ptr = ambig_id_ct;
        genocounts[common_geno] = sample_ct - common_decr;
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
#ifdef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    ++subgroup_idx;
    uintptr_t cur_raregeno_word = *raregeno_workspace_iter++;
    while (1) {
#ifndef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      const uintptr_t cur_geno = cur_raregeno_word & 3;
      if (IsSet(sample_include, raw_sample_idx)) {
        genocounts[cur_geno] += 1;
        ++common_decr;
      }
      if (ambig_sample_ids && (cur_geno == 3)) {
        ambig_sample_ids[ambig_id_ct++] = raw_sample_idx;
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      --remaining_deltas_in_subgroup;
      cur_raregeno_word >>= 2;
    }
  }
}

// 1-bit, unsubsetted: count 1-bit array, then count raregeno
// 1-bit, subsetted: count [1-bit array AND sample_include], iterate through
//   difflist
PglErr CountparseOnebitSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, const unsigned char** fread_pp, uint32_t* __restrict ambig_sample_ids, uint32_t* __restrict ambig_id_ct_ptr, uint32_t* __restrict genocounts, uintptr_t* __restrict raregeno_workspace) {
  const uint32_t initial_bitarray_byte_ct = DivUp(raw_sample_ct, CHAR_BIT);
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) <= initial_bitarray_byte_ct) {
    return kPglRetMalformedInput;
  }
  const unsigned char* fread_difflist_start = &((*fread_pp)[1 + initial_bitarray_byte_ct]);
  const uint32_t common2_code = *((*fread_pp)++);
  const uint32_t geno_code_low = common2_code / 4;
  const uint32_t geno_code_high = (common2_code & 3) + geno_code_low;
#ifdef __arm__
#  error "Unaligned accesses in CountparseOnebitSubset()."
#endif
  const uintptr_t* onebitarr = R_CAST(const uintptr_t*, *fread_pp);
  uint32_t high_geno_ct;
  if (raw_sample_ct == sample_ct) {
    high_geno_ct = PopcountBytes(*fread_pp, initial_bitarray_byte_ct);
  } else {
    high_geno_ct = PopcountBytesMasked(*fread_pp, sample_include, initial_bitarray_byte_ct);
  }
  *fread_pp = fread_difflist_start;
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno_workspace, &group_info_iter, &difflist_len);
  *ambig_id_ct_ptr = 0;
  ZeroUiArr(4, genocounts);
  if (reterr || (!difflist_len)) {
    genocounts[geno_code_low] = sample_ct - high_geno_ct;
    genocounts[geno_code_high] = high_geno_ct;
    return reterr;
  }
  if (raw_sample_ct == sample_ct) {
    ZeroTrailingQuaters(difflist_len, raregeno_workspace);
    GenovecCountFreqsUnsafe(raregeno_workspace, difflist_len, genocounts);
    if (ambig_sample_ids && genocounts[3]) {
      // no need to update ambig_sample_ids[], but necessary to set ambig_id_ct
      // and fread_pp to enable rarealt counting.
      reterr = SkipDifflistIds(fread_end, group_info_iter, difflist_len, raw_sample_ct, fread_pp);
      *ambig_id_ct_ptr = genocounts[3];
    }
    genocounts[geno_code_low] += sample_ct - difflist_len - high_geno_ct;
    genocounts[geno_code_high] += high_geno_ct;
    return kPglRetSuccess;
  }
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  uintptr_t* raregeno_workspace_iter = raregeno_workspace;
  uintptr_t raw_sample_idx = 0;
  uint32_t subgroup_idx = 0;
  uint32_t ambig_id_ct = 0;
  uint32_t common_decr = 0;
  while (1) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        *ambig_id_ct_ptr = ambig_id_ct;
        genocounts[geno_code_low] += sample_ct - common_decr - high_geno_ct;
        genocounts[geno_code_high] += high_geno_ct;
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
#ifdef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    ++subgroup_idx;
    uintptr_t cur_raregeno_word = *raregeno_workspace_iter++;
    while (1) {
#ifndef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      const uintptr_t cur_geno = cur_raregeno_word & 3;
      if (IsSet(sample_include, raw_sample_idx)) {
        genocounts[cur_geno] += 1;
        ++common_decr;
        high_geno_ct -= IsSet(onebitarr, raw_sample_idx);
      }
      if (ambig_sample_ids && (cur_geno == 3)) {
        ambig_sample_ids[ambig_id_ct++] = raw_sample_idx;
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      --remaining_deltas_in_subgroup;
      cur_raregeno_word >>= 2;
    }
  }
}

// fread_pp should be non-null iff this is being called by an internal function
// gathering rarealt counts, or dosages, as well
PglErr GetRefalt1GenotypeCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uint32_t* genocounts) {
  // genocounts[0] := ref/ref, genocounts[1] := ref/alt1,
  // genocounts[2] := alt1/alt1, genocounts[3] := missing/other
  assert(vidx < pgrp->fi.raw_variant_ct);
  assert(sample_ct);
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  const uint32_t multiallelic_relevant = fread_pp && VrtypeMultiallelic(vrtype);
  if (VrtypeLdCompressed(vrtype)) {
    // LD compression
    PglErr reterr = LdLoadGenovecSubsetIfNecessary(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp);
    if (reterr) {
      return reterr;
    }
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    if (!(pgrp->ldbase_stypes & kfPgrLdcacheRefalt1Genocounts)) {
      ZeroTrailingQuaters(sample_ct, pgrp->ldbase_genovec);
      GenovecCountFreqsUnsafe(pgrp->ldbase_genovec, sample_ct, pgrp->ldbase_refalt1_genocounts);
      pgrp->ldbase_stypes |= kfPgrLdcacheRefalt1Genocounts;
    }
    memcpy(genocounts, pgrp->ldbase_refalt1_genocounts, 4 * sizeof(int32_t));
    uint32_t ambig_id_ct_filtered = 0;
    reterr = LdSubsetAdjustGenocounts(fread_end, subsetting_required? sample_include : nullptr, sample_include_cumulative_popcounts, pgrp->ldbase_genovec, raw_sample_ct, &fread_ptr, genocounts, (subsetting_required && multiallelic_relevant)? pgrp->workspace_ambig_sample_ids : nullptr, &(pgrp->workspace_ambig_id_ct), &ambig_id_ct_filtered, pgrp->workspace_raregeno_tmp_loadbuf);
    if (!subsetting_required) {
      pgrp->workspace_ambig_id_ct = ambig_id_ct_filtered;
    }
    if (vrtype & 1) {
      // inverted
      const uint32_t tmpval = genocounts[0];
      genocounts[0] = genocounts[2];
      genocounts[2] = tmpval;
    }
    if (fread_pp) {
      *fread_pp = fread_ptr;
      *fread_endp = fread_end;
    }
    return reterr;
  }
  const unsigned char* fread_ptr;
  const unsigned char* fread_end = nullptr;  // maybe-uninitialized warning
  if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
    return kPglRetReadFail;
  }
  const uint32_t is_ldbase = pgrp->fi.vrtypes && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
  if (is_ldbase) {
    // difflists are very efficient to count directly when not subsetting
    // (since we can entirely ignore the sample IDs), but it's often better to
    // unpack them first when subsetting.

    // ...er, the statement above is a lie, unpack-first almost always seems to
    // be better.
    pgrp->ldbase_vidx = vidx;
    // this may be slowed down by the LD caching change.
    PglErr reterr = ParseNonLdGenovecSubsetUnsafe(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, vrtype, multiallelic_relevant, &fread_ptr, pgrp, pgrp->ldbase_genovec);
    if (reterr) {
      return reterr;
    }
    ZeroTrailingQuaters(sample_ct, pgrp->ldbase_genovec);
    GenovecCountFreqsUnsafe(pgrp->ldbase_genovec, sample_ct, genocounts);
    memcpy(pgrp->ldbase_refalt1_genocounts, genocounts, 4 * sizeof(int32_t));
    pgrp->ldbase_stypes = kfPgrLdcacheQuater | kfPgrLdcacheRefalt1Genocounts;
    if (fread_pp) {
      *fread_pp = fread_ptr;
      *fread_endp = fread_end;
    }
    return kPglRetSuccess;
  }
  if (VrtypeDifflist(vrtype)) {
    PglErr reterr = CountparseDifflistSubset(fread_end, sample_include, vrtype & 3, raw_sample_ct, sample_ct, &fread_ptr, multiallelic_relevant? pgrp->workspace_ambig_sample_ids : nullptr, &(pgrp->workspace_ambig_id_ct), genocounts, pgrp->workspace_raregeno_tmp_loadbuf);
    if (fread_pp) {
      *fread_pp = fread_ptr;
      *fread_endp = fread_end;
    }
    return reterr;
  }
  if (vrtype & 1) {
    PglErr reterr = CountparseOnebitSubset(fread_end, sample_include, raw_sample_ct, sample_ct, &fread_ptr, multiallelic_relevant? pgrp->workspace_ambig_sample_ids : nullptr, &(pgrp->workspace_ambig_id_ct), genocounts, pgrp->workspace_raregeno_tmp_loadbuf);
    if (fread_pp) {
      *fread_pp = fread_ptr;
      *fread_endp = fread_end;
    }
    return reterr;
  }
  const uint32_t genovec_byte_ct = QuaterCtToByteCt(raw_sample_ct);
  if (S_CAST(uintptr_t, fread_end - fread_ptr) < genovec_byte_ct) {
    return kPglRetMalformedInput;
  }
  const unsigned char* fread_2bit_end = &(fread_ptr[genovec_byte_ct]);
  const uint32_t fread_ptr_unaligned = R_CAST(uintptr_t, fread_ptr) % kBytesPerVec;
  if (!subsetting_required) {
    if (fread_ptr_unaligned) {
      GenoarrCountFreqs(fread_ptr, raw_sample_ct, genocounts);
    } else {
      GenovecCountFreqs(R_CAST(const uintptr_t*, fread_ptr), raw_sample_ct, genocounts);
    }
  } else {
    if (fread_ptr_unaligned) {
      GenoarrCountSubsetFreqs(fread_ptr, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
    } else {
      GenovecCountSubsetFreqs(R_CAST(const uintptr_t*, fread_ptr), sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
    }
  }
  if (multiallelic_relevant) {
    ExtractGenoarrAmbigIds(R_CAST(const uintptr_t*, fread_ptr), raw_sample_ct, pgrp->workspace_ambig_sample_ids, &(pgrp->workspace_ambig_id_ct));
    *fread_pp = fread_2bit_end;
    *fread_endp = fread_end;
  } else if (vrtype == kPglVrtypePlink1) {
    // [3] -> [0]
    // [2] -> [1]
    // [1] -> [3]
    // [0] -> [2]
    const uint32_t save2 = genocounts[0];
    const uint32_t save3 = genocounts[1];
    genocounts[0] = genocounts[3];
    genocounts[1] = genocounts[2];
    genocounts[2] = save2;
    genocounts[3] = save3;
  }
  return kPglRetSuccess;
}

PglErr PgrGetRefalt1GenotypeCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uint32_t* genocounts) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    ZeroUiArr(4, genocounts);
    return kPglRetSuccess;
  }
  return GetRefalt1GenotypeCounts(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genocounts);
}


PglErr ParseAux1(const unsigned char* fread_end, uint32_t alt_allele_ct, const unsigned char** fread_pp, PgenReader* pgrp, uint32_t* aux1_nonmissing_ct_ptr) {
  // Assumes pgr.workspace_ambig_id_ct has the correct entry count.
  // Fills pgr.workspace_aux1_nonmissing_vec and (usually)
  //   pgr.workspace_aux1_code_vec, zeroes trailing bits of the former.
  // aux1_nonmissing_ct_ptr can be set to nullptr to skip past this track
  //   instead of copying it into aux1_code_vec.
  const uint32_t ambig_id_ct = pgrp->workspace_ambig_id_ct;
  uintptr_t* aux1_nonmissing_vec = pgrp->workspace_aux1_nonmissing_vec;
  const uint32_t aux1_nonmissing_byte_ct = DivUp(ambig_id_ct, CHAR_BIT);
  const unsigned char* fread_ptr = *fread_pp;
  if (S_CAST(uintptr_t, fread_end - fread_ptr) < aux1_nonmissing_byte_ct) {
    return kPglRetMalformedInput;
  }
  memcpy(aux1_nonmissing_vec, fread_ptr, aux1_nonmissing_byte_ct);
  ZeroTrailingBits(ambig_id_ct, aux1_nonmissing_vec);
  fread_ptr = &(fread_ptr[aux1_nonmissing_byte_ct]);
  const uint32_t aux1_nonmissing_ct = PopcountWords(aux1_nonmissing_vec, BitCtToWordCt(ambig_id_ct));
  const uint32_t aux1_allele_entry_bytect = GetAux1AlleleByteCt(alt_allele_ct, aux1_nonmissing_ct);
  if (S_CAST(uintptr_t, fread_end - fread_ptr) < aux1_allele_entry_bytect) {
    return kPglRetMalformedInput;
  }
  if (aux1_nonmissing_ct_ptr) {
    *aux1_nonmissing_ct_ptr = aux1_nonmissing_ct;
    memcpy(pgrp->workspace_aux1_code_vec, fread_ptr, aux1_allele_entry_bytect);
  }
  *fread_pp = &(fread_ptr[aux1_allele_entry_bytect]);
  return kPglRetSuccess;
}

static_assert(kPglMaxAltAlleleCt == 254, "Need to update Aux1UpdateAlleleCounts().");
void Aux1UpdateAlleleCounts(uint32_t alt_allele_ct, uint32_t aux1_nonmissing_ct, uintptr_t* aux1_code_vec, uint32_t* allele_ct_buf) {
  // aux1_code_vec not const since we might zero the trailing bits
  // todo: validate?
  ZeroUiArr(alt_allele_ct - 1, &(allele_ct_buf[2]));
  if (!aux1_nonmissing_ct) {
    return;
  }
  if (alt_allele_ct < 4) {
    assert(alt_allele_ct >= 2);
    uint32_t aux1_counts[4];
    uint32_t code_vec_entry_ct = aux1_nonmissing_ct * (alt_allele_ct - 1);
    ZeroTrailingQuaters(code_vec_entry_ct, aux1_code_vec);
    GenovecCountFreqsUnsafe(aux1_code_vec, code_vec_entry_ct, aux1_counts);
    allele_ct_buf[0] += aux1_counts[0];
    allele_ct_buf[1] += aux1_counts[1];
    if (alt_allele_ct == 2) {
      allele_ct_buf[2] += aux1_counts[0] + aux1_counts[1] + 2 * aux1_counts[2];
    } else {
      allele_ct_buf[2] += aux1_counts[2];
      allele_ct_buf[3] += aux1_counts[3];
    }
    return;
  }
  const uintptr_t* aux1_code_vec_iter = aux1_code_vec;
  const uint32_t aux1_nonmissing_allele_ct = 2 * aux1_nonmissing_ct;
  // Slightly different code must be used for 256 <= alt_allele_ct < 4096, and
  // 65536 <= alt_allele_ct < 2^24.
  assert(alt_allele_ct <= kPglMaxAltAlleleCt);
  uint32_t halfcode_bit_width;
  uint32_t log2_halfcodes_per_word;
  if (alt_allele_ct < 16) {
    halfcode_bit_width = 4;
    log2_halfcodes_per_word = kBitsPerWordLog2 - 2;
  } else {
    halfcode_bit_width = 8;
    log2_halfcodes_per_word = kBitsPerWordLog2 - 3;
  }
  const uintptr_t* aux1_code_vec_last = &(aux1_code_vec[aux1_nonmissing_allele_ct >> log2_halfcodes_per_word]);
  const uint32_t halfcode_mask = (1 << halfcode_bit_width) - 1;
  uint32_t block_len_m1 = (1 << log2_halfcodes_per_word) - 1;
  while (1) {
    if (aux1_code_vec_iter >= aux1_code_vec_last) {
      if (aux1_code_vec_iter > aux1_code_vec_last) {
        return;
      }
      block_len_m1 &= aux1_nonmissing_allele_ct - 1;
    }
    uintptr_t cur_aux_word = *aux1_code_vec_iter++;
    uint32_t aux_idx_lowbits = 0;
    while (1) {
      allele_ct_buf[cur_aux_word & halfcode_mask] += 1;
      if (aux_idx_lowbits == block_len_m1) {
        break;
      }
      ++aux_idx_lowbits;
      cur_aux_word >>= halfcode_bit_width;
    }
  }
}

static_assert(kPglMaxAltAlleleCt == 254, "Need to update Aux1SubsetUpdateAlleleCounts().");
void Aux1SubsetUpdateAlleleCounts(const uint32_t* __restrict ambig_sample_ids, const uintptr_t* __restrict aux1_nonmissing_vec, const uintptr_t* __restrict aux1_code_vec, const uintptr_t* __restrict sample_include, uint32_t alt_allele_ct, uint32_t aux1_nonmissing_ct, uint32_t* allele_ct_buf) {
  // todo: validate?
  ZeroUiArr(alt_allele_ct - 1, &(allele_ct_buf[2]));
  uint32_t ambig_idx = 0;
  if (alt_allele_ct == 2) {
    for (uint32_t aux_idx = 0; aux_idx < aux1_nonmissing_ct; ++aux_idx, ++ambig_idx) {
      FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
      uint32_t sample_idx = ambig_sample_ids[ambig_idx];
      if (IsSet(sample_include, sample_idx)) {
        allele_ct_buf[GetQuaterarrEntry(aux1_code_vec, aux_idx)] += 1;
      }
    }
    allele_ct_buf[2] += aux1_nonmissing_ct;
    return;
  }
  assert(alt_allele_ct <= kPglMaxAltAlleleCt);
  uint32_t log2_codes_per_word;
  uint32_t halfcode_bit_width;
  if (alt_allele_ct == 3) {
    log2_codes_per_word = kBitsPerWordLog2 - 2;
    halfcode_bit_width = 2;
  } else if (alt_allele_ct < 16) {
    log2_codes_per_word = kBitsPerWordLog2 - 3;
    halfcode_bit_width = 4;
  } else {
    log2_codes_per_word = kBitsPerWordLog2 - 4;
    halfcode_bit_width = 8;
  }
  const uint32_t idx_mask = (1 << log2_codes_per_word) - 1;
  const uint32_t entry_bit_ct = halfcode_bit_width * 2;
  const uintptr_t halfcode_mask = (k1LU << halfcode_bit_width) - 1;
  for (uint32_t aux_idx = 0; aux_idx < aux1_nonmissing_ct; ++aux_idx, ++ambig_idx) {
    FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
    const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
    if (IsSet(sample_include, sample_idx)) {
      const uintptr_t cur_code_unmasked = aux1_code_vec[aux_idx >> log2_codes_per_word] >> (entry_bit_ct * (aux_idx & idx_mask));
      allele_ct_buf[cur_code_unmasked & halfcode_mask] += 1;
      allele_ct_buf[(cur_code_unmasked >> halfcode_bit_width) & halfcode_mask] += 1;
    }
  }
}

void Aux1UpdateRefOrAlt1Countvec(const uint32_t* __restrict ambig_sample_ids, const uintptr_t* aux1_nonmissing_vec, const uintptr_t* __restrict aux1_code_vec, uint32_t alt_allele_ct, uint32_t aux1_nonmissing_ct, uint32_t allele_idx, uintptr_t* __restrict allele_countvec) {
  if (!aux1_nonmissing_ct) {
    return;
  }
  uint32_t log2_codes_per_word;
  uint32_t code_bit_width;
  if (alt_allele_ct == 2) {
    log2_codes_per_word = kBitsPerWordLog2 - 1;
    code_bit_width = 2;
  } else if (alt_allele_ct == 3) {
    log2_codes_per_word = kBitsPerWordLog2 - 2;
    code_bit_width = 4;
  } else if (alt_allele_ct < 16) {
    log2_codes_per_word = kBitsPerWordLog2 - 3;
    code_bit_width = 8;
  } else {
    log2_codes_per_word = kBitsPerWordLog2 - 4;
    code_bit_width = 16;
  }

  // The "+ (code_bit_width == 2)" is needed to handle alt_allele_ct == 2
  // correctly.
  // This code may need to be changed when we increase the alt allele count
  // limit, since a natural 256-4095 alt allele idx representation uses 12 bits
  // per halfcode, which is not a nice power of two.  We might decide at that
  // point that code simplicity is worth bloating this part of the file by 33%
  // (i.e. use 16 bits per halfcode), but the necessary code isn't very
  // complicated...
  const uint32_t halfcode_mask = (1 << ((code_bit_width / 2) + (code_bit_width == 2))) - 1;
  const uintptr_t* aux1_code_vec_iter = aux1_code_vec;
  const uintptr_t* aux1_code_vec_last = &(aux1_code_vec[(aux1_nonmissing_ct - 1) >> log2_codes_per_word]);
  uint32_t ambig_idx = 0;
  uint32_t block_len_m1 = (1 << log2_codes_per_word) - 1;
  while (1) {
    if (aux1_code_vec_iter >= aux1_code_vec_last) {
      if (aux1_code_vec_iter > aux1_code_vec_last) {
        return;
      }
      block_len_m1 &= aux1_nonmissing_ct - 1;
    }
    uintptr_t aux1_code_word = *aux1_code_vec_iter++;
    uint32_t aux_idx_lowbits = 0;
    while (1) {
      FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
      const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
      const uintptr_t cur_allele_ct = ((aux1_code_word & halfcode_mask) == allele_idx);
      AssignQuaterarrEntry(sample_idx, cur_allele_ct, allele_countvec);
      if (aux_idx_lowbits == block_len_m1) {
        break;
      }
      ++aux_idx_lowbits;
      aux1_code_word >>= code_bit_width;
    }
  }
}

void Aux1UpdateRefOrAlt1CountvecSubset(const uint32_t* __restrict ambig_sample_ids, const uintptr_t* aux1_nonmissing_vec, const uintptr_t* __restrict aux1_code_vec, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t alt_allele_ct, uint32_t aux1_nonmissing_ct, uint32_t allele_idx, uintptr_t* __restrict allele_countvec) {
  uint32_t log2_codes_per_word;
  uint32_t code_bit_width;
  if (alt_allele_ct == 2) {
    log2_codes_per_word = kBitsPerWordLog2 - 1;
    code_bit_width = 2;
  } else if (alt_allele_ct == 3) {
    log2_codes_per_word = kBitsPerWordLog2 - 2;
    code_bit_width = 4;
  } else if (alt_allele_ct < 16) {
    log2_codes_per_word = kBitsPerWordLog2 - 3;
    code_bit_width = 8;
  } else {
    log2_codes_per_word = kBitsPerWordLog2 - 4;
    code_bit_width = 16;
  }
  const uint32_t halfcode_mask = (1 << ((code_bit_width / 2) + (code_bit_width == 2))) - 1;
  uint32_t ambig_idx = 0;
  for (uint32_t aux_idx = 0; aux_idx < aux1_nonmissing_ct; ++aux_idx, ++ambig_idx) {
    FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
    const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
    if (IsSet(sample_include, sample_idx)) {
      const uint32_t cur_code = (aux1_code_vec[aux_idx >> log2_codes_per_word] >> ((code_bit_width * aux_idx) & (kBitsPerWord - 1))) & halfcode_mask;
      const uintptr_t cur_allele_ct = (cur_code == allele_idx);
      AssignQuaterarrEntry(RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_idx), cur_allele_ct, allele_countvec);
    }
  }
}

// "rarealt" = alt2/alt3/etc.
void Aux1UpdateRarealtCountvec(const uint32_t* __restrict ambig_sample_ids, const uintptr_t* aux1_nonmissing_vec, const uintptr_t* __restrict aux1_code_vec, uint32_t alt_allele_ct, uint32_t aux1_nonmissing_ct, uintptr_t allele_idx, uintptr_t* __restrict allele_countvec) {
  // todo: check whether promoting allele_idx to uintptr_t actually helps
  if (!aux1_nonmissing_ct) {
    return;
  }
  const uint32_t aux1_nonmissing_ct_m1 = aux1_nonmissing_ct - 1;
  const uintptr_t* aux1_code_vec_iter = aux1_code_vec;
  uint32_t ambig_idx = 0;
  if (alt_allele_ct == 2) {
    assert(allele_idx == 2);
    const uintptr_t* aux1_code_vec_last = &(aux1_code_vec_iter[aux1_nonmissing_ct_m1 / kBitsPerWordD2]);
    uint32_t block_len_m1 = kBitsPerWordD2 - 1;
    while (1) {
      if (aux1_code_vec_iter >= aux1_code_vec_last) {
        if (aux1_code_vec_iter > aux1_code_vec_last) {
          return;
        }
        block_len_m1 &= aux1_nonmissing_ct_m1;
      }
      uintptr_t aux1_code_word = *aux1_code_vec_iter++;
      uint32_t aux_idx_lowbits = 0;
      while (1) {
        FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
        const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
        uintptr_t cur_allele_ct = 1 + ((aux1_code_word & 3) == allele_idx);
        AssignQuaterarrEntry(sample_idx, cur_allele_ct, allele_countvec);
        if (aux_idx_lowbits == block_len_m1) {
          break;
        }
        ++aux_idx_lowbits;
        aux1_code_word >>= 2;
      }
    }
  }
  uint32_t log2_codes_per_word;
  uint32_t halfcode_bit_width;
  if (alt_allele_ct == 3) {
    log2_codes_per_word = kBitsPerWordLog2 - 2;
    halfcode_bit_width = 4;
  } else if (alt_allele_ct < 16) {
    log2_codes_per_word = kBitsPerWordLog2 - 3;
    halfcode_bit_width = 8;
  } else {
    log2_codes_per_word = kBitsPerWordLog2 - 4;
    halfcode_bit_width = 16;
  }
  const uint32_t code_bit_width = halfcode_bit_width * 2;
  const uintptr_t halfcode_mask = (1 << halfcode_bit_width) - 1;
  const uintptr_t* aux1_code_vec_last = &(aux1_code_vec_iter[aux1_nonmissing_ct_m1 >> log2_codes_per_word]);
  uint32_t block_len_m1 = (1 << log2_codes_per_word) - 1;
  while (1) {
    if (aux1_code_vec_iter >= aux1_code_vec_last) {
      if (aux1_code_vec_iter > aux1_code_vec_last) {
        return;
      }
      block_len_m1 &= aux1_nonmissing_ct_m1;
    }
    uintptr_t aux1_code_word = *aux1_code_vec_iter++;
    uint32_t aux_idx_lowbits = 0;
    while (1) {
      FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
      const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
      const uintptr_t cur_allele_ct = ((aux1_code_word & halfcode_mask) == allele_idx) + (((aux1_code_word >> halfcode_bit_width) & halfcode_mask) == allele_idx);
      AssignQuaterarrEntry(sample_idx, cur_allele_ct, allele_countvec);
      if (aux_idx_lowbits == block_len_m1) {
        break;
      }
      ++aux_idx_lowbits;
      aux1_code_word >>= code_bit_width;
    }
  }
}

void Aux1UpdateRarealtCountvecSubset(const uint32_t* __restrict ambig_sample_ids, const uintptr_t* aux1_nonmissing_vec, const uintptr_t* __restrict aux1_code_vec, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t alt_allele_ct, uint32_t aux1_nonmissing_ct, uint32_t allele_idx, uintptr_t* __restrict allele_countvec) {
  uint32_t ambig_idx = 0;
  if (alt_allele_ct == 2) {
    assert(allele_idx == 2);
    for (uint32_t aux_idx = 0; aux_idx < aux1_nonmissing_ct; ++aux_idx, ++ambig_idx) {
      FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
      const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
      if (IsSet(sample_include, sample_idx)) {
        const uint32_t cur_code = GetQuaterarrEntry(aux1_code_vec, aux_idx);
        const uintptr_t cur_allele_ct = 1 + (cur_code == allele_idx);
        AssignQuaterarrEntry(RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_idx), cur_allele_ct, allele_countvec);
      }
    }
  }
  uint32_t log2_codes_per_word;
  uint32_t halfcode_bit_width;
  if (alt_allele_ct == 3) {
    log2_codes_per_word = kBitsPerWordLog2 - 2;
    halfcode_bit_width = 2;
  } else if (alt_allele_ct < 16) {
    log2_codes_per_word = kBitsPerWordLog2 - 3;
    halfcode_bit_width = 4;
  } else {
    log2_codes_per_word = kBitsPerWordLog2 - 4;
    halfcode_bit_width = 8;
  }
  const uint32_t code_bit_width = 2 * halfcode_bit_width;
  const uintptr_t halfcode_mask = (k1LU << halfcode_bit_width) - 1;
  for (uint32_t aux_idx = 0; aux_idx < aux1_nonmissing_ct; ++aux_idx, ++ambig_idx) {
    FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
    const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
    if (IsSet(sample_include, sample_idx)) {
      const uintptr_t cur_code_unmasked = aux1_code_vec[aux_idx >> log2_codes_per_word] >> ((code_bit_width * aux_idx) % kBitsPerWord);
      const uintptr_t cur_allele_ct = ((cur_code_unmasked & halfcode_mask) == allele_idx) + (((cur_code_unmasked >> halfcode_bit_width) & halfcode_mask) == allele_idx);
      AssignQuaterarrEntry(RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_idx), cur_allele_ct, allele_countvec);
    }
  }
}

// See comments toward end of PgrReadGenovecSubsetThenCommon2().
void Aux1UpdateGenovecMatch2Unsafe(const uint32_t* __restrict ambig_sample_ids, const uintptr_t* aux1_nonmissing_vec, const uintptr_t* __restrict aux1_code_vec, uint32_t alt_allele_ct, uint32_t most_common_idx, uint32_t second_most_common_idx, uint32_t aux1_nonmissing_ct, uintptr_t* __restrict genovec) {
  assert(aux1_nonmissing_ct);
  // One rarealt, one ref/alt1.  Only two codes to match, since
  // homozygous-ref/alt1 calls are not stored in this data track.
  // This code is separate from the ...match3() functions below since, if we
  // iterate over a variable-length array, we probably don't get to take
  // advantage of registers.  But maybe the compiler is actually smart enough;
  // need to test the simple implementation...
  const uint32_t rarealt_is_minor = (most_common_idx < 2);
  uint32_t nonrare_idx;
  uint32_t rarealt_idx;
  if (rarealt_is_minor) {
    nonrare_idx = most_common_idx;
    rarealt_idx = second_most_common_idx;
  } else {
    nonrare_idx = second_most_common_idx;
    rarealt_idx = most_common_idx;
  }
  uint32_t log2_codes_per_word;
  uint32_t code_bit_width;
  if (alt_allele_ct == 2) {
    log2_codes_per_word = kBitsPerWordLog2 - 1;
    code_bit_width = 2;
  } else if (alt_allele_ct == 3) {
    log2_codes_per_word = kBitsPerWordLog2 - 2;
    code_bit_width = 4;
  } else if (alt_allele_ct < 16) {
    log2_codes_per_word = kBitsPerWordLog2 - 3;
    code_bit_width = 8;
  } else {
    log2_codes_per_word = kBitsPerWordLog2 - 4;
    code_bit_width = 16;
  }
  const uintptr_t* aux1_code_vec_iter = aux1_code_vec;
  const uintptr_t* aux1_code_vec_last = &(aux1_code_vec[(aux1_nonmissing_ct - 1) >> log2_codes_per_word]);
  const uint32_t halfcode_bit_width = (code_bit_width / 2) + (code_bit_width == 2);
  const uint32_t code1 = (rarealt_idx << halfcode_bit_width) + nonrare_idx;
  const uint32_t code02 = rarealt_idx * ((1 << halfcode_bit_width) + 1);
  const uintptr_t store02 = rarealt_is_minor * 2;
  assert(code_bit_width < 32);
  const uint32_t code_mask = (1 << code_bit_width) - 1;
  uint32_t ambig_idx = 0;
  uint32_t block_len_m1 = (1 << log2_codes_per_word) - 1;
  while (1) {
    if (aux1_code_vec_iter >= aux1_code_vec_last) {
      if (aux1_code_vec_iter > aux1_code_vec_last) {
        return;
      }
      block_len_m1 &= aux1_nonmissing_ct - 1;
    }
    uintptr_t aux1_code_word = *aux1_code_vec_iter++;
    uint32_t aux_idx_lowbits = 0;
    while (1) {
      FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
      const uint32_t cur_code = aux1_code_word & code_mask;
      const uint32_t match1 = (cur_code == code1);
      if (match1 || (cur_code == code02)) {
        const uint32_t sample_idx = ambig_sample_ids[ambig_idx];

        // todo: check if there's a better way to perform this assignment
        // e.g. (rarealt_is_minor ^ match1) + rarealt_is_minor
        const uintptr_t new_geno = match1? 1 : store02;

        AssignQuaterarrEntry(sample_idx, new_geno, genovec);
      }
      if (aux_idx_lowbits == block_len_m1) {
        break;
      }
      ++aux_idx_lowbits;
      aux1_code_word >>= code_bit_width;
    }
  }
}

void Aux1UpdateGenovecSubsetMatch2(const uint32_t* __restrict ambig_sample_ids, const uintptr_t* aux1_nonmissing_vec, const uintptr_t* __restrict aux1_code_vec, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t alt_allele_ct, uint32_t most_common_idx, uint32_t second_most_common_idx, uint32_t aux1_nonmissing_ct, uintptr_t* __restrict genovec) {
  const uint32_t rarealt_is_minor = (most_common_idx < 2);
  uint32_t nonrare_idx;
  uint32_t rarealt_idx;
  if (rarealt_is_minor) {
    nonrare_idx = most_common_idx;
    rarealt_idx = second_most_common_idx;
  } else {
    nonrare_idx = second_most_common_idx;
    rarealt_idx = most_common_idx;
  }
  uint32_t log2_codes_per_word;
  uint32_t code_bit_width;
  if (alt_allele_ct == 2) {
    log2_codes_per_word = kBitsPerWordLog2 - 1;
    code_bit_width = 2;
  } else if (alt_allele_ct == 3) {
    log2_codes_per_word = kBitsPerWordLog2 - 2;
    code_bit_width = 4;
  } else if (alt_allele_ct < 16) {
    log2_codes_per_word = kBitsPerWordLog2 - 3;
    code_bit_width = 8;
  } else {
    log2_codes_per_word = kBitsPerWordLog2 - 4;
    code_bit_width = 16;
  }
  const uint32_t halfcode_bit_width = (code_bit_width / 2) + (code_bit_width == 2);
  const uint32_t code1 = (rarealt_idx << halfcode_bit_width) + nonrare_idx;
  const uint32_t code02 = rarealt_idx * ((1 << halfcode_bit_width) + 1);
  const uintptr_t store02 = rarealt_is_minor * 2;
  assert(code_bit_width < 32);
  const uint32_t code_mask = (1 << code_bit_width) - 1;
  uint32_t ambig_idx = 0;
  for (uint32_t aux_idx = 0; aux_idx < aux1_nonmissing_ct; ++aux_idx, ++ambig_idx) {
    FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
    const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
    if (IsSet(sample_include, sample_idx)) {
      const uint32_t cur_code = (aux1_code_vec[aux_idx >> log2_codes_per_word] >> ((code_bit_width * aux_idx) & (kBitsPerWord - 1))) & code_mask;
      const uint32_t match1 = (cur_code == code1);
      if (match1 || (cur_code == code02)) {
        const uintptr_t new_geno = match1? 1 : store02;
        AssignQuaterarrEntry(RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_idx), new_geno, genovec);
      }
    }
  }
}

void Aux1UpdateGenovecMatch3Unsafe(const uint32_t* __restrict ambig_sample_ids, const uintptr_t* aux1_nonmissing_vec, const uintptr_t* __restrict aux1_code_vec, uint32_t alt_allele_ct, uint32_t most_common_idx, uint32_t second_most_common_idx, uint32_t aux1_nonmissing_ct, uintptr_t* __restrict genovec) {
  assert(aux1_nonmissing_ct);
  // can't have two rarealts be the most common if there is only one rarealt...
  assert(alt_allele_ct > 2);
  uint32_t log2_codes_per_word;
  uint32_t code_bit_width;
  if (alt_allele_ct == 3) {
    log2_codes_per_word = kBitsPerWordLog2 - 2;
    code_bit_width = 4;
  } else if (alt_allele_ct < 16) {
    log2_codes_per_word = kBitsPerWordLog2 - 3;
    code_bit_width = 8;
  } else {
    log2_codes_per_word = kBitsPerWordLog2 - 4;
    code_bit_width = 16;
  }
  const uintptr_t* aux1_code_vec_iter = aux1_code_vec;
  const uintptr_t* aux1_code_vec_last = &(aux1_code_vec[(aux1_nonmissing_ct - 1) >> log2_codes_per_word]);
  const uint32_t halfcode_bit_width = code_bit_width / 2;
  assert(halfcode_bit_width <= 16);
  const uint32_t code0 = most_common_idx * ((1 << halfcode_bit_width) + 1);
  const uint32_t code1 = (most_common_idx << halfcode_bit_width) + second_most_common_idx;
  const uint32_t code2 = second_most_common_idx * ((1 << halfcode_bit_width) + 1);
  assert(code_bit_width < 32);
  const uint32_t code_mask = (1 << code_bit_width) - 1;
  uint32_t ambig_idx = 0;
  uint32_t block_len_m1 = (1 << log2_codes_per_word) - 1;
  while (1) {
    if (aux1_code_vec_iter >= aux1_code_vec_last) {
      if (aux1_code_vec_iter > aux1_code_vec_last) {
        return;
      }
      block_len_m1 &= aux1_nonmissing_ct - 1;
    }
    uintptr_t aux1_code_word = *aux1_code_vec_iter++;
    uint32_t aux_idx_lowbits = 0;
    while (1) {
      FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
      const uint32_t cur_code = aux1_code_word & code_mask;
      const uint32_t match0 = (cur_code == code0);
      const uint32_t match1 = (cur_code == code1);
      if (match0 || match1 || (cur_code == code2)) {
        const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
        const uintptr_t new_geno = match0? 0 : (2 - match1);
        AssignQuaterarrEntry(sample_idx, new_geno, genovec);
      }
      if (aux_idx_lowbits == block_len_m1) {
        break;
      }
      ++aux_idx_lowbits;
      aux1_code_word >>= code_bit_width;
    }
  }
}

void Aux1UpdateGenovecSubsetMatch3(const uint32_t* __restrict ambig_sample_ids, const uintptr_t* aux1_nonmissing_vec, const uintptr_t* __restrict aux1_code_vec, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t alt_allele_ct, uint32_t most_common_idx, uint32_t second_most_common_idx, uint32_t aux1_nonmissing_ct, uintptr_t* __restrict genovec) {
  assert(aux1_nonmissing_ct);
  // can't have two rarealts be the most common if there is only one rarealt...
  assert(alt_allele_ct > 2);
  uint32_t log2_codes_per_word;
  uint32_t code_bit_width;
  if (alt_allele_ct == 3) {
    log2_codes_per_word = kBitsPerWordLog2 - 2;
    code_bit_width = 4;
  } else if (alt_allele_ct < 16) {
    log2_codes_per_word = kBitsPerWordLog2 - 3;
    code_bit_width = 8;
  } else {
    log2_codes_per_word = kBitsPerWordLog2 - 4;
    code_bit_width = 16;
  }
  const uint32_t halfcode_bit_width = code_bit_width / 2;
  assert(halfcode_bit_width <= 16);
  const uint32_t code0 = most_common_idx * ((1 << halfcode_bit_width) + 1);
  const uint32_t code1 = (most_common_idx << halfcode_bit_width) + second_most_common_idx;
  const uint32_t code2 = second_most_common_idx * ((1 << halfcode_bit_width) + 1);
  assert(code_bit_width < 32);
  const uint32_t code_mask = (1 << code_bit_width) - 1;
  uint32_t ambig_idx = 0;
  for (uint32_t aux_idx = 0; aux_idx < aux1_nonmissing_ct; ++aux_idx, ++ambig_idx) {
    FindFirst1BitFromU32(aux1_nonmissing_vec, &ambig_idx);
    const uint32_t sample_idx = ambig_sample_ids[ambig_idx];
    if (IsSet(sample_include, sample_idx)) {
      const uint32_t cur_code = (aux1_code_vec[aux_idx >> log2_codes_per_word] >> ((code_bit_width * aux_idx) & (kBitsPerWord - 1))) & code_mask;
      const uint32_t match0 = (cur_code == code0);
      const uint32_t match1 = (cur_code == code1);
      if (match0 || match1 || (cur_code == code2)) {
        const uintptr_t new_geno = match0? 0 : (2 - match1);
        AssignQuaterarrEntry(RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_idx), new_geno, genovec);
      }
    }
  }
}

PglErr PgrReadGenovecSubsetThenCommon2(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uint32_t* __restrict maj_allele_idx_ptr, uint32_t* __restrict second_allele_idx_ptr, uint32_t* __restrict allele_ct_buf) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  const uint32_t allele_ct = pgrp->fi.allele_idx_offsets[vidx + 1] - pgrp->fi.allele_idx_offsets[vidx];
  if (!sample_ct) {
    *maj_allele_idx_ptr = 0;
    *second_allele_idx_ptr = 0;
    ZeroUiArr(allele_ct, allele_ct_buf);
    return kPglRetSuccess;
  }
  // major allele corresponds to 0 bits, second-most-common allele corresponds
  // to 1.
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRefalt1GenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, &fread_ptr, &fread_end, genovec);
  if (reterr) {
    return reterr;
  }
  ZeroTrailingQuaters(sample_ct, genovec);
  uint32_t bothset_ct;
  GenovecAlleleCtsUnsafe(genovec, sample_ct, allele_ct_buf, &bothset_ct);
  uint32_t second_most_common_idx;
  if (allele_ct == 2) {
  PgrReadGenovecSubsetThenCommon2_refalt1_finish:
    second_most_common_idx = (allele_ct_buf[0] >= allele_ct_buf[1]);
    *maj_allele_idx_ptr = 1 - second_most_common_idx;
    *second_allele_idx_ptr = second_most_common_idx;
    if (!second_most_common_idx) {
      GenovecInvertUnsafe(sample_ct, genovec);
      ZeroTrailingQuaters(sample_ct, genovec);
    }
    return kPglRetSuccess;
  }
  const uint32_t ambig_id_ct = pgrp->workspace_ambig_id_ct;
  second_most_common_idx = (allele_ct_buf[0] >= allele_ct_buf[1]);
  uint32_t included_ambig_id_ct_threshold = DivUp(allele_ct_buf[second_most_common_idx], 2);
  // avoid processing the aux1 data track if possible.
  const uint32_t subsetting_required = (sample_ct != pgrp->fi.raw_sample_ct);
  if (ambig_id_ct < included_ambig_id_ct_threshold) {
    goto PgrReadGenovecSubsetThenCommon2_refalt1_finish;
  }
  if (subsetting_required) {
    const uint32_t* __restrict ambig_sample_ids = pgrp->workspace_ambig_sample_ids;
    uint32_t included_ambig_id_ct = 0;
    uint32_t ambig_idx = 0;
    // minor optimization: check included_ambig_id_ct halfway through, we might
    // be able to skip the second half of the list.
    const uint32_t half_ambig_id_ct = ambig_id_ct / 2;
    for (; ambig_idx < half_ambig_id_ct; ++ambig_idx) {
      included_ambig_id_ct += IsSet(sample_include, ambig_sample_ids[ambig_idx]);
    }
    if (included_ambig_id_ct < included_ambig_id_ct_threshold) {
      for (; ambig_idx < ambig_id_ct; ++ambig_idx) {
        included_ambig_id_ct += IsSet(sample_include, ambig_sample_ids[ambig_idx]);
      }
      if (included_ambig_id_ct < included_ambig_id_ct_threshold) {
        goto PgrReadGenovecSubsetThenCommon2_refalt1_finish;
      }
    }
  }
  const uint32_t alt_allele_ct = allele_ct - 1;
  uint32_t aux1_nonmissing_ct;
  if (ParseAux1(fread_end, alt_allele_ct, &fread_ptr, pgrp, &aux1_nonmissing_ct)) {
    return kPglRetMalformedInput;
  }
  if (subsetting_required) {
    Aux1SubsetUpdateAlleleCounts(pgrp->workspace_ambig_sample_ids, pgrp->workspace_aux1_nonmissing_vec, pgrp->workspace_aux1_code_vec, sample_include, alt_allele_ct, aux1_nonmissing_ct, allele_ct_buf);
  } else {
    Aux1UpdateAlleleCounts(alt_allele_ct, aux1_nonmissing_ct, pgrp->workspace_aux1_code_vec, allele_ct_buf);
  }
  uint32_t most_common_idx = 0;
  uint32_t most_common_allele_ct = allele_ct_buf[0];  // ref count
  second_most_common_idx = 1;
  uint32_t second_most_common_allele_ct = 0;
  for (uint32_t allele_idx = 1; allele_idx <= alt_allele_ct; ++allele_idx) {
    uint32_t cur_allele_ct = allele_ct_buf[allele_idx];
    if (cur_allele_ct > second_most_common_allele_ct) {
      if (cur_allele_ct > most_common_allele_ct) {
        second_most_common_allele_ct = most_common_allele_ct;
        second_most_common_idx = most_common_idx;
        most_common_allele_ct = cur_allele_ct;
        most_common_idx = allele_idx;
      } else {
        second_most_common_allele_ct = cur_allele_ct;
        second_most_common_idx = allele_idx;
      }
    }
  }
  if (most_common_idx + second_most_common_idx == 1) {
    goto PgrReadGenovecSubsetThenCommon2_refalt1_finish;
  }
  uint32_t rarealt_is_minor = (most_common_idx < 2);
  if (rarealt_is_minor || (second_most_common_idx < 2)) {
    // One of the most common alleles is ref or alt1, and the other is a
    // rarealt.  Suppose for clarity that ref and alt2 are the two most common.
    // 1. Keep just the hom ref calls in the base genotype vector.  het
    //    ref/alt1 and hom alt1 are converted to missing.
    // 2. Search aux1_code_vec for het ref/alt2 and hom alt2 genotypes.  Update
    //    the corresponding positions in genovec.
    uint32_t ref_or_alt1_idx = rarealt_is_minor? most_common_idx : second_most_common_idx;
    if (ref_or_alt1_idx == 1) {
      GenovecInvertUnsafe(sample_ct, genovec);
    }
    GenovecNonzeroToMissingUnsafe(sample_ct, genovec);
    if (!rarealt_is_minor) {
      GenovecInvertUnsafe(sample_ct, genovec);
    }
    if (subsetting_required) {
      Aux1UpdateGenovecSubsetMatch2(pgrp->workspace_ambig_sample_ids, pgrp->workspace_aux1_nonmissing_vec, pgrp->workspace_aux1_code_vec, sample_include, sample_include_cumulative_popcounts, alt_allele_ct, most_common_idx, second_most_common_idx, aux1_nonmissing_ct, genovec);
    } else {
      Aux1UpdateGenovecMatch2Unsafe(pgrp->workspace_ambig_sample_ids, pgrp->workspace_aux1_nonmissing_vec, pgrp->workspace_aux1_code_vec, alt_allele_ct, most_common_idx, second_most_common_idx, aux1_nonmissing_ct, genovec);
    }
  }
  // Both of the most common alleles are rarealts.  Supposing for clarity that
  // they are alt2 and alt3,
  // 1. Initialize genovec to all-missing.
  // 2. Search aux1_code_vec for hom alt2, het alt2/alt3, and hom alt3
  //    genotypes.  Update the corresponding positions in genovec.
  SetAllBits(sample_ct * 2, genovec);
  if (subsetting_required) {
    Aux1UpdateGenovecSubsetMatch3(pgrp->workspace_ambig_sample_ids, pgrp->workspace_aux1_nonmissing_vec, pgrp->workspace_aux1_code_vec, sample_include, sample_include_cumulative_popcounts, alt_allele_ct, most_common_idx, second_most_common_idx, aux1_nonmissing_ct, genovec);
  } else {
    Aux1UpdateGenovecMatch3Unsafe(pgrp->workspace_ambig_sample_ids, pgrp->workspace_aux1_nonmissing_vec, pgrp->workspace_aux1_code_vec, alt_allele_ct, most_common_idx, second_most_common_idx, aux1_nonmissing_ct, genovec);
  }
  return kPglRetSuccess;
}

PglErr ParseDifflistJustAmbigIds(const unsigned char* fread_end, PgenReader* pgrp, const unsigned char** fread_pp) {
  // Side effects: uses pgr.workspace_raregeno_tmp_loadbuf.
  uintptr_t* __restrict raregeno_iter = pgrp->workspace_raregeno_tmp_loadbuf;
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, pgrp->fi.raw_sample_ct, fread_pp, raregeno_iter, &group_info_iter, &difflist_len);
  if (reterr || (!difflist_len)) {
    pgrp->workspace_ambig_id_ct = 0;
    return reterr;
  }
  // variant is guaranteed to be multiallelic, so little point in optimizing
  // for sparsity.
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  uint32_t* __restrict ambig_sample_ids = pgrp->workspace_ambig_sample_ids;
  uint32_t subgroup_idx = 0;
  uint32_t subgroup_len_m1 = kBitsPerWordD2 - 1;
  uint32_t ambig_id_ct = 0;
  uintptr_t raw_sample_idx = 0;
  while (1) {
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        pgrp->workspace_ambig_id_ct = ambig_id_ct;
        return kPglRetSuccess;
      }
      subgroup_len_m1 &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
#ifdef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    ++subgroup_idx;
    // invert so we can compare vs. zero
    uintptr_t cur_raregeno_word_inv = ~(*raregeno_iter++);
    uint32_t difflist_idx_lowbits = 0;
    while (1) {
#ifndef __LP64__
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
#endif
      if (!(cur_raregeno_word_inv & 3)) {
        ambig_sample_ids[ambig_id_ct++] = raw_sample_idx;
      }
      if (difflist_idx_lowbits == subgroup_len_m1) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      ++difflist_idx_lowbits;
      cur_raregeno_word_inv >>= 2;
    }
  }
}

PglErr ParseJustAmbigIds(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp) {
  assert(pgrp->fi.vrtypes);
  // Just initializes pgr.workspace_ambig_sample_ids and
  // pgr.workspace_ambig_id_ct.  Avoids some unnecessary work when we're just
  // interested in rare alt(s).  May use pgrp->workspace_vec.
  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  assert(VrtypeMultiallelic(vrtype));
  const uint32_t maintrack_vrtype = vrtype & 7;
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  const uint32_t is_ld_compressed = VrtypeLdCompressed(maintrack_vrtype);
  const uint32_t is_ldbase = pgrp->fi.vrtypes && (!is_ld_compressed) && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
  if (!maintrack_vrtype) {
    const uint32_t genovec_byte_ct = QuaterCtToByteCt(raw_sample_ct);
    if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < genovec_byte_ct) {
      return kPglRetMalformedInput;
    }
    const unsigned char* fread_ptr_new = &((*fread_pp)[genovec_byte_ct]);
#ifdef __arm__
#  error "Unaligned accesses in ParseJustAmbigIds()."
#endif
    const uintptr_t* cur_genoarr = R_CAST(const uintptr_t*, *fread_pp);
    if (is_ldbase) {
      if (!subsetting_required) {
        memcpy(pgrp->ldbase_genovec, cur_genoarr, genovec_byte_ct);
        cur_genoarr = pgrp->ldbase_genovec;  // may as well guarantee alignment
      } else {
        CopyQuaterarrNonemptySubset(cur_genoarr, sample_include, raw_sample_ct, sample_ct, pgrp->ldbase_genovec);
      }
    }
    ExtractGenoarrAmbigIds(cur_genoarr, raw_sample_ct, pgrp->workspace_ambig_sample_ids, &(pgrp->workspace_ambig_id_ct));
    *fread_pp = fread_ptr_new;
    return kPglRetSuccess;
  }
  if (is_ld_compressed) {
    // LD compression
    PglErr reterr = LdLoadMinimalSubsetIfNecessary(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp);
    if (reterr) {
      return reterr;
    }
  }
  if (!is_ldbase) {
    // difflist storage, unnecessary to save ldbase info
    if (maintrack_vrtype == 1) {
      *fread_pp += (raw_sample_ct + 15) / CHAR_BIT;
    }
    return ParseDifflistJustAmbigIds(fread_end, pgrp, fread_pp);
  }
  if (maintrack_vrtype == 1) {
    uintptr_t* cur_genoarr = subsetting_required? pgrp->workspace_vec : pgrp->ldbase_genovec;
    if (ParseOnebitUnsafe(fread_end, 1, fread_pp, pgrp, cur_genoarr)) {
      return kPglRetMalformedInput;
    }
    if (subsetting_required) {
      CopyQuaterarrNonemptySubset(cur_genoarr, sample_include, raw_sample_ct, sample_ct, pgrp->ldbase_genovec);
    }
    return kPglRetSuccess;
  }
  PglErr reterr = ParseAndSaveDifflistSubset(fread_end, sample_include, sample_include_cumulative_popcounts, raw_sample_ct, fread_pp, pgrp->workspace_ambig_sample_ids, pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, &(pgrp->ldbase_difflist_len), &(pgrp->workspace_ambig_id_ct), pgrp->workspace_vec);
  if (reterr) {
    return reterr;
  }
  pgrp->ldbase_difflist_sample_ids[pgrp->ldbase_difflist_len] = sample_ct;
  return kPglRetSuccess;
}

PglErr PgrReadAlleleCountvecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_countvec) {
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  if (allele_idx < 2) {
    const uint32_t is_multiallelic = VrtypeMultiallelic(vrtype);
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    PglErr reterr = ReadRefalt1GenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, is_multiallelic? (&fread_ptr) : nullptr, &fread_end, allele_countvec);
    if (reterr) {
      return reterr;
    }
    if (!allele_idx) {
      GenovecInvertUnsafe(sample_ct, allele_countvec);
    }
    if (!is_multiallelic) {
      return kPglRetSuccess;
    }
    const uint32_t alt_allele_ct = pgrp->fi.allele_idx_offsets[vidx + 1] - pgrp->fi.allele_idx_offsets[vidx];
    uint32_t aux1_nonmissing_ct;
    if (ParseAux1(fread_end, alt_allele_ct, &fread_ptr, pgrp, &aux1_nonmissing_ct)) {
      return kPglRetReadFail;
    }
    if (subsetting_required) {
      Aux1UpdateRefOrAlt1CountvecSubset(pgrp->workspace_ambig_sample_ids, pgrp->workspace_aux1_nonmissing_vec, pgrp->workspace_aux1_code_vec, sample_include, sample_include_cumulative_popcounts, alt_allele_ct, aux1_nonmissing_ct, allele_idx, allele_countvec);
    } else {
      Aux1UpdateRefOrAlt1Countvec(pgrp->workspace_ambig_sample_ids, pgrp->workspace_aux1_nonmissing_vec, pgrp->workspace_aux1_code_vec, alt_allele_ct, aux1_nonmissing_ct, allele_idx, allele_countvec);
    }
    return kPglRetSuccess;
  }
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
    return kPglRetReadFail;
  }
  assert(VrtypeMultiallelic(vrtype));
  const uint32_t alt_allele_ct = pgrp->fi.allele_idx_offsets[vidx + 1] - pgrp->fi.allele_idx_offsets[vidx];
  assert(allele_idx <= alt_allele_ct);
  uint32_t* ambig_sample_ids = pgrp->workspace_ambig_sample_ids;
  if ((vrtype & 7) == 7) {
    // most values missing, 0b11 entries
    SetAllUlArr(QuaterCtToWordCt(sample_ct), allele_countvec);
    uintptr_t* __restrict raregeno_vec = pgrp->workspace_raregeno_vec;
    uint32_t* __restrict difflist_sample_ids = pgrp->workspace_difflist_sample_ids;
    uint32_t difflist_len;
    PglErr reterr = ParseAndSaveDifflist(fread_end, raw_sample_ct, &fread_ptr, raregeno_vec, difflist_sample_ids, &difflist_len);
    if (reterr || (!difflist_len)) {
      return reterr;
    }
    uintptr_t* raregeno_iter = raregeno_vec;
    uintptr_t raregeno_word = *raregeno_iter++;
    uint32_t difflist_idx = 0;
    uint32_t ambig_id_ct = 0;
    if (subsetting_required) {
      while (1) {
        const uint32_t cur_raregeno = raregeno_word & 3;
        const uint32_t sample_idx = difflist_sample_ids[difflist_idx];
        ++difflist_idx;
        if (cur_raregeno == 3) {
          ambig_sample_ids[ambig_id_ct++] = sample_idx;
        } else if (IsSet(sample_include, sample_idx)) {
          const uint32_t subsetted_pos = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_idx);
          ClearQuaterarrEntry(subsetted_pos, allele_countvec);
        }
        if (difflist_idx == difflist_len) {
          break;
        }
        raregeno_word >>= 2;
        if (!(difflist_idx % kBitsPerWordD2)) {
          raregeno_word = *raregeno_iter++;
        }
      }
    } else {
      while (1) {
        const uint32_t cur_raregeno = raregeno_word & 3;
        const uint32_t sample_idx = difflist_sample_ids[difflist_idx];
        ++difflist_idx;
        if (cur_raregeno == 3) {
          ambig_sample_ids[ambig_id_ct++] = sample_idx;
        } else {
          ClearQuaterarrEntry(sample_idx, allele_countvec);
        }
        if (difflist_idx == difflist_len) {
          break;
        }
        raregeno_word >>= 2;
        if (!(difflist_idx % kBitsPerWordD2)) {
          raregeno_word = *raregeno_iter++;
        }
      }
    }
    pgrp->workspace_ambig_id_ct = ambig_id_ct;
  } else {
    PglErr reterr = ParseJustAmbigIds(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, &fread_ptr);
    if (reterr) {
      return reterr;
    }
    if (VrtypeLdCompressed(vrtype)) {
      // todo: optimize difflist case
      if (!(pgrp->ldbase_stypes & kfPgrLdcacheQuater)) {
        assert(pgrp->ldbase_stypes & kfPgrLdcacheDifflist);
        PgrDifflistToGenovecUnsafe(pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, pgrp->fi.vrtypes[pgrp->ldbase_vidx] & 3, sample_ct, pgrp->ldbase_difflist_len, pgrp->ldbase_genovec);
        pgrp->ldbase_stypes |= kfPgrLdcacheQuater;
      }
      CopyQuaterarr(pgrp->ldbase_genovec, sample_ct, allele_countvec);
      GenovecNonmissingToZeroUnsafe(sample_ct, allele_countvec);
    } else {
      // most values nonmissing, 0b00 entries
      ZeroUlArr(QuaterCtToWordCt(sample_ct), allele_countvec);
    }
  }
  uint32_t aux1_nonmissing_ct;
  if (ParseAux1(fread_end, alt_allele_ct, &fread_ptr, pgrp, &aux1_nonmissing_ct)) {
    return kPglRetReadFail;
  }
  uint32_t ambig_id_ct = pgrp->workspace_ambig_id_ct;
  uint32_t ambig_missing_ct = ambig_id_ct - aux1_nonmissing_ct;
  const uintptr_t* __restrict aux1_nonmissing_vec = pgrp->workspace_aux1_nonmissing_vec;
  if (subsetting_required) {
    uint32_t ambig_idx = 0;
    for (uint32_t ambig_missing_idx = 0; ambig_missing_idx < ambig_missing_ct; ++ambig_missing_idx) {
      FindFirst0BitFromU32(aux1_nonmissing_vec, &ambig_idx);
      AssignQuaterarrEntry(RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, ambig_sample_ids[ambig_idx]), 3, allele_countvec);
    }
    Aux1UpdateRarealtCountvecSubset(ambig_sample_ids, aux1_nonmissing_vec, pgrp->workspace_aux1_code_vec, sample_include, sample_include_cumulative_popcounts, alt_allele_ct, aux1_nonmissing_ct, allele_idx, allele_countvec);
  } else {
    uint32_t ambig_idx = 0;
    for (uint32_t ambig_missing_idx = 0; ambig_missing_idx < ambig_missing_ct; ++ambig_missing_idx) {
      FindFirst0BitFromU32(aux1_nonmissing_vec, &ambig_idx);
      AssignQuaterarrEntry(ambig_sample_ids[ambig_idx], 3, allele_countvec);
    }
    Aux1UpdateRarealtCountvec(ambig_sample_ids, aux1_nonmissing_vec, pgrp->workspace_aux1_code_vec, alt_allele_ct, aux1_nonmissing_ct, allele_idx, allele_countvec);
  }
  return kPglRetSuccess;
}

void DetectGenovecHetsHw(const uintptr_t*__restrict genovec, uint32_t raw_sample_ctl2, Halfword* all_hets_hw) {
  // requires trailing bits of genovec to be zeroed out.  does not update last
  // all_hets[] halfword if raw_sample_ctl2 is odd.
  for (uint32_t widx = 0; widx < raw_sample_ctl2; ++widx) {
    const uintptr_t cur_word = genovec[widx];
    uintptr_t ww = (~(cur_word >> 1)) & cur_word & kMask5555;  // low 1, high 0
    all_hets_hw[widx] = PackWordToHalfword(ww);
  }
}

PglErr ParseAndApplyDifflistHphaseSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t multiallelic_relevant, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict all_hets) {
  // Side effects: uses pgr.workspace_raregeno_tmp_loadbuf.
  // If multiallelic_relevant is true, a list of sample indices with freshly
  // loaded raregeno value 0b11 is saved to pgr.workspace_ambig_sample_ids, and
  // pgr.workspace_ambig_id_ct is set to the length of the list.
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  if (sample_ct == raw_sample_ct) {
    PglErr reterr = ParseAndApplyDifflist(fread_end, multiallelic_relevant, fread_pp, pgrp, genovec);
    if (reterr) {
      return reterr;
    }
    PgrDetectGenovecHets(genovec, raw_sample_ct, all_hets);
    return kPglRetSuccess;
  }
  uintptr_t* cur_raregeno_iter = pgrp->workspace_raregeno_tmp_loadbuf;
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, cur_raregeno_iter, &group_info_iter, &difflist_len);
  if (reterr || (!difflist_len)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  uint32_t* ambig_sample_ids = multiallelic_relevant? pgrp->workspace_ambig_sample_ids : nullptr;
  uintptr_t raw_sample_idx = 0;
  uint32_t ambig_id_ct = 0;
  uint32_t subgroup_idx = 0;
  while (1) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        pgrp->workspace_ambig_id_ct = ambig_id_ct;
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
      raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    ++subgroup_idx;
    uintptr_t cur_raregeno_word = *cur_raregeno_iter++;
    // This loop tends to be the decompression bottleneck.  Tried to modify it
    // to process 4 entries at a time, but that didn't end up helping.
    while (1) {
      // always check, since otherwise AssignQuaterarrEntry() can scribble
      // over arbitrary memory
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
      const uintptr_t cur_geno = cur_raregeno_word & 3;
      if (IsSet(sample_include, raw_sample_idx)) {
        AssignQuaterarrEntry(RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, raw_sample_idx), cur_geno, genovec);
      }
      if (cur_geno == 1) {
        SetBit(raw_sample_idx, all_hets);
      } else {
        ClearBit(raw_sample_idx, all_hets);  // needed for LD decompression
        if (multiallelic_relevant && (cur_geno == 3)) {
          ambig_sample_ids[ambig_id_ct++] = raw_sample_idx;
        }
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      --remaining_deltas_in_subgroup;
      cur_raregeno_word >>= 2;
    }
  }
}

PglErr ParseNonLdGenovecHphaseSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vrtype, uint32_t difflist_ambig_ids_needed, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict all_hets) {
  // If all_hets is nullptr, this is essentially identical to
  // ParseNonLdGenovecSubsetUnsafe().
  // Side effects:
  //   may use pgrp->workspace_raregeno_tmp_loadbuf
  //   may use pgrp->workspace_vec (subsetting)
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  if (!VrtypeDifflist(vrtype)) {
    const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
    uintptr_t* raw_genovec = subsetting_required? pgrp->workspace_vec : genovec;
    PglErr reterr = Parse1or2bitGenovecUnsafe(fread_end, vrtype, difflist_ambig_ids_needed, fread_pp, pgrp, raw_genovec);
    // can technically avoid this return and just plow forward in error case,
    // but that tiny optimization is not worth the associated maintenance
    // problems
    if (reterr) {
      return reterr;
    }
    ZeroTrailingQuaters(raw_sample_ct, raw_genovec);
    if (all_hets) {
      PgrDetectGenovecHetsUnsafe(raw_genovec, QuaterCtToWordCt(raw_sample_ct), all_hets);
    }
    if (subsetting_required) {
      CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, genovec);
    }
    return kPglRetSuccess;
  }
  const uint32_t vrtype_low2 = vrtype & 3;
  const uint32_t word_ct = QuaterCtToWordCt(sample_ct);
  memset(genovec, vrtype_low2 * 0x55, word_ct * kBytesPerWord);
  ZeroTrailingQuaters(sample_ct, genovec);
  // common genotype can't be het
  if (all_hets) {
    ZeroUlArr(BitCtToWordCt(raw_sample_ct), all_hets);
    return ParseAndApplyDifflistHphaseSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, difflist_ambig_ids_needed, fread_pp, pgrp, genovec, all_hets);
  } else {
    return ParseAndApplyDifflistSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, difflist_ambig_ids_needed, fread_pp, pgrp, genovec);
  }
}

// may use workspace_vec
PglErr ld_load_genovec_hphase_subset_if_necessary(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp) {
  // bugfix: this conditional was in the other order, which was wrong since
  // we depend on LdLoadNecessary() to set pgrp->ldbase_vidx as a side effect
  // todo: make that explicit instead of a side effect...
  if (LdLoadNecessary(vidx, pgrp) || (!(pgrp->ldbase_stypes & kfPgrLdcacheAllHets))) {
    const uint32_t ldbase_vidx = pgrp->ldbase_vidx;
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (InitReadPtrs(ldbase_vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    pgrp->ldbase_stypes = kfPgrLdcacheQuater;
    return ParseNonLdGenovecHphaseSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, pgrp->fi.vrtypes[ldbase_vidx], 0, &fread_ptr, pgrp, pgrp->ldbase_genovec, pgrp->ldbase_all_hets);
  }
  if (!(pgrp->ldbase_stypes & kfPgrLdcacheQuater)) {
    assert(pgrp->ldbase_stypes & kfPgrLdcacheDifflist);
    PgrDifflistToGenovecUnsafe(pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, pgrp->fi.vrtypes[pgrp->ldbase_vidx] & 3, sample_ct, pgrp->ldbase_difflist_len, pgrp->ldbase_genovec);
    pgrp->ldbase_stypes |= kfPgrLdcacheQuater;
  }
  return kPglRetSuccess;
}

// "h" in hphase = "hardcall"
// No need for fread_pp/fread_endp, since any function which needed them would
// use a multiallelic variant-supporting phase loader.
// Iff *phasepresent_ct_ptr is nonzero, returned phasepresent is guaranteed to
// only have set bits for het calls where phase information is present.
/*
PglErr PgrReadRefalt1GenovecHphaseRawUnsafe(uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phaseraw, uint32_t* phasepresent_ct_ptr) {
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  if (!VrtypeHphase(vrtype)) {
    // don't bother updating ldbase_all_hets, too much of a performance
    // penalty, and too likely that we won't need it
    *phasepresent_ct_ptr = 0;
    return ReadRefalt1GenovecSubsetUnsafe(nullptr, nullptr, raw_sample_ct, vidx, pgrp, nullptr, nullptr, genovec);
  }
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t multiallelic_relevant = VrtypeMultiallelic(vrtype);
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  uintptr_t* all_hets = pgrp->workspace_all_hets;
  if (VrtypeLdCompressed(vrtype)) {
    // ldbase_all_hets not needed in this case
    PglErr reterr = LdLoadGenovecSubsetIfNecessary(nullptr, nullptr, raw_sample_ct, vidx, pgrp);
    if (reterr) {
      return reterr;
    }
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    CopyQuaterarr(pgrp->ldbase_genovec, raw_sample_ct, genovec);

    reterr = ParseAndApplyDifflist(fread_end, multiallelic_relevant, &fread_ptr, pgrp, genovec);
    if (reterr) {
      return reterr;
    }
    PgrDetectGenovecHets(genovec, raw_sample_ct, all_hets);
    if ((vrtype & 7) == 3) {
      GenovecInvertUnsafe(raw_sample_ct, genovec);
    }
  } else {
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    PglErr reterr = ParseNonLdGenovecHphaseSubset(fread_end, nullptr, nullptr, raw_sample_ct, vrtype, multiallelic_relevant, &fread_ptr, pgrp, genovec, all_hets);
    if (reterr) {
      return reterr;
    }
    const uint32_t is_ldbase = pgrp->fi.vrtypes && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
    if (is_ldbase) {
      CopyQuaterarr(genovec, raw_sample_ct, pgrp->ldbase_genovec);
      pgrp->ldbase_vidx = vidx;
      pgrp->ldbase_stypes = kfPgrLdcacheQuater;
    }
  }
  if (multiallelic_relevant) {
    // todo
    // can't ignore multiallelic track, since it may contain additional het
    // calls.  (these additional calls must *not* be saved to ldbase_all_hets.)
    return kPglRetNotYetSupported;
  }
  const uint32_t het_ct = PopcountWords(all_hets, raw_sample_ctl);
  if (!het_ct) {
    // there shouldn't be a hphase track at all in this case
    return kPglRetMalformedInput;
  }
  const uint32_t het_ctdl = het_ct / kBitsPerWord;
  phaseraw[het_ctdl] = 0;
  const uint32_t first_half_byte_ct = 1 + (het_ct / CHAR_BIT);
  memcpy(phaseraw, fread_ptr, first_half_byte_ct);
  if (!(fread_ptr[0] & 1)) {
    // phase always present, phasepresent not stored
    *phasepresent_ct_ptr = het_ct;
    return kPglRetSuccess;
  }
  const uint32_t raw_phasepresent_ct = PopcountWords(phaseraw, het_ctdl + 1) - 1;
  if (!raw_phasepresent_ct) {
    // there shouldn't be a hphase track at all in this case, either
    return kPglRetMalformedInput;
  }
  // put this in a phasepresent-independent location, to make things more
  // convenient for the caller
  memcpy(&(phaseraw[1 + (raw_sample_ct / kBitsPerWord)]), &(fread_ptr[first_half_byte_ct]), DivUp(raw_phasepresent_ct, CHAR_BIT));
  *phasepresent_ct_ptr = raw_phasepresent_ct;
  return kPglRetSuccess;
}
*/

// If fread_pp/fread_endp are non-null, this always moves fread_ptr to the end
// of aux2.  Set phasepresent/phaseinfo to nullptr when you don't actually care
// about the contents of aux2.
PglErr ReadRefalt1GenovecHphaseSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* phasepresent_ct_ptr) {
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t multiallelic_relevant = VrtypeMultiallelic(vrtype);
  if (!VrtypeHphase(vrtype)) {
    // don't bother updating ldbase_all_hets, too much of a performance
    // penalty, and too likely that we won't need it
    *phasepresent_ct_ptr = 0;
    PglErr reterr = ReadRefalt1GenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, fread_pp, fread_endp, genovec);
    if ((!multiallelic_relevant) || (!fread_pp) || reterr) {
      return reterr;
    }
    return ParseAux1(*fread_endp, pgrp->fi.allele_idx_offsets[vidx + 1] - pgrp->fi.allele_idx_offsets[vidx] - 1, fread_pp, pgrp, nullptr);
  }
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  uintptr_t* all_hets = pgrp->workspace_all_hets;
  if (VrtypeLdCompressed(vrtype)) {
    PglErr reterr;
    if (!subsetting_required) {
      // ldbase_all_hets not needed in this case
      reterr = LdLoadGenovecSubsetIfNecessary(nullptr, nullptr, sample_ct, vidx, pgrp);
    } else {
      reterr = ld_load_genovec_hphase_subset_if_necessary(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp);
      memcpy(all_hets, pgrp->ldbase_all_hets, raw_sample_ctl * sizeof(intptr_t));
    }
    if (reterr) {
      return reterr;
    }
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    CopyQuaterarr(pgrp->ldbase_genovec, sample_ct, genovec);
    reterr = ParseAndApplyDifflistHphaseSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, multiallelic_relevant, &fread_ptr, pgrp, genovec, all_hets);
    if (reterr) {
      return reterr;
    }
    if ((vrtype & 7) == 3) {
      GenovecInvertUnsafe(sample_ct, genovec);
    }
  } else {
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    PglErr reterr = ParseNonLdGenovecHphaseSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, vrtype, multiallelic_relevant, &fread_ptr, pgrp, genovec, all_hets);
    if (reterr) {
      return reterr;
    }
    const uint32_t is_ldbase = pgrp->fi.vrtypes && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
    if (is_ldbase) {
      CopyQuaterarr(genovec, sample_ct, pgrp->ldbase_genovec);
      pgrp->ldbase_vidx = vidx;
      pgrp->ldbase_stypes = kfPgrLdcacheQuater;
      if (subsetting_required) {
        pgrp->ldbase_stypes |= kfPgrLdcacheAllHets;
        memcpy(pgrp->ldbase_all_hets, all_hets, raw_sample_ctl * sizeof(intptr_t));
      }
    }
  }
  if (multiallelic_relevant) {
    // todo
    // can't ignore multiallelic track, since it may contain additional het
    // calls.  (these additional calls must *not* be saved to ldbase_all_hets.)
    return kPglRetNotYetSupported;
  }
  const uint32_t het_ct = PopcountWords(all_hets, raw_sample_ctl);
  if (!het_ct) {
    // there shouldn't be a hphase track at all in this case, het_ct is not
    // computed off a subset
    return kPglRetMalformedInput;
  }
  if (fread_pp) {
    *fread_endp = fread_end;
  }
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  if (!(fread_ptr[0] & 1)) {
    // phase always present
    if (fread_pp) {
      *fread_pp = &(fread_ptr[1 + (het_ct / CHAR_BIT)]);
      if (!phaseinfo) {
        // for internal callers which just want to skip aux2
        return kPglRetSuccess;
      }
    }
    if (!subsetting_required) {
      memcpy(phasepresent, all_hets, raw_sample_ctl * kBytesPerWord);
      *phasepresent_ct_ptr = het_ct;
      ExpandBytearr(fread_ptr, all_hets, raw_sample_ctl, het_ct, 1, phaseinfo);
    } else {
      // we could drop the "phasepresent bit can only be set at hets" guarantee
      // and speed up this case, but I doubt it's worth it
      CopyBitarrSubset(all_hets, sample_include, sample_ct, phasepresent);
      *phasepresent_ct_ptr = PopcountWords(phasepresent, sample_ctl);
      if (!(*phasepresent_ct_ptr)) {
        // bugfix (7 Dec 2017): clear sample_ctl words here, not raw_sample_ctl
        ZeroUlArr(sample_ctl, phaseinfo);
      } else {
        ExpandThenSubsetBytearr(fread_ptr, all_hets, sample_include, het_ct, sample_ct, 1, phaseinfo);
      }
    }
    return kPglRetSuccess;
  }
  const uint32_t het_ctdl = het_ct / kBitsPerWord;

  // explicit phasepresent
  const uintptr_t* aux2_first_half = R_CAST(const uintptr_t*, fread_ptr);
  uintptr_t* aux2_first_half_copy = pgrp->workspace_vec;
  aux2_first_half_copy[het_ctdl] = 0;
  memcpy(aux2_first_half_copy, aux2_first_half, 1 + (het_ct / CHAR_BIT));
  const uint32_t raw_phasepresent_ct = PopcountWords(aux2_first_half_copy, het_ctdl + 1) - 1;
  if (!raw_phasepresent_ct) {
    // there shouldn't be a hphase track at all in this case
    return kPglRetMalformedInput;
  }
  if (fread_pp) {
    *fread_pp = &(fread_ptr[1 + (het_ct / CHAR_BIT) + DivUp(raw_phasepresent_ct, CHAR_BIT)]);
    if (!phaseinfo) {
      return kPglRetSuccess;
    }
  }
  const unsigned char* aux2_second_half = &(fread_ptr[1 + (het_ct / CHAR_BIT)]);
  if (!subsetting_required) {
    ExpandBytearrNested(aux2_second_half, aux2_first_half_copy, all_hets, sample_ctl, raw_phasepresent_ct, 1, phasepresent, phaseinfo);
    *phasepresent_ct_ptr = raw_phasepresent_ct;
  } else {
    ExpandThenSubsetBytearrNested(aux2_second_half, aux2_first_half_copy, all_hets, sample_include, sample_ct, raw_phasepresent_ct, 1, phasepresent, phaseinfo);
    *phasepresent_ct_ptr = PopcountWords(phasepresent, sample_ctl);
  }
  return kPglRetSuccess;
}

PglErr PgrReadRefalt1GenovecHphaseSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* phasepresent_ct_ptr) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    *phasepresent_ct_ptr = 0;
    return kPglRetSuccess;
  }
  return ReadRefalt1GenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genovec, phasepresent, phaseinfo, phasepresent_ct_ptr);
}


// similar to ParseAndSaveDifflist()
PglErr ParseAndSaveDeltalistAsBitarr(const unsigned char* fread_end, uint32_t raw_sample_ct, const unsigned char** fread_pp, uintptr_t* deltalist_include, uint32_t* deltalist_len_ptr) {
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, deltalist_len_ptr);
  const uint32_t deltalist_len = *deltalist_len_ptr;
  if (reterr || (!deltalist_len)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(raw_sample_ct);
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t group_idx_last = (deltalist_len - 1) / kPglDifflistGroupSize;
  ZeroUlArr(raw_sample_ctl, deltalist_include);
  uint32_t group_len_m1 = kPglDifflistGroupSize - 1;
  uint32_t group_idx = 0;
  while (1) {
    if (group_idx >= group_idx_last) {
      if (group_idx > group_idx_last) {
        return kPglRetSuccess;
      }
      group_len_m1 &= deltalist_len - 1;
    }
    uintptr_t raw_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
    group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    ++group_idx;
    uint32_t raw_deltalist_idx_lowbits = 0;
    while (1) {
      // always check, otherwise we may scribble over arbitrary memory
      if (raw_sample_idx >= raw_sample_ct) {
        return kPglRetMalformedInput;
      }
      SetBit(raw_sample_idx, deltalist_include);
      if (raw_deltalist_idx_lowbits == group_len_m1) {
        break;
      }
      ++raw_deltalist_idx_lowbits;
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
  }
}

PglErr ParseDosage16(const unsigned char* fread_ptr, const unsigned char* fread_end, const uintptr_t* __restrict sample_include, uint32_t sample_ct, uint32_t vidx, uint32_t alt_allele_ct, PgenReader* pgrp, uint32_t* dosage_ct_ptr, uintptr_t* __restrict dosage_present, uint16_t* dosage_vals) {
  // Side effect: may use pgrp->workspace_dosage_present
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* raw_dosage_present = subsetting_required? pgrp->workspace_dosage_present : dosage_present;
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t is_unconditional_dosage = ((vrtype & 0x60) == 0x40);
  uint32_t raw_dosage_ct;
  if ((vrtype & 0x60) == 0x20) {
    // case 1: dosage list
    if (ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, &fread_ptr, raw_dosage_present, &raw_dosage_ct)) {
      return kPglRetMalformedInput;
    }
  } else if (is_unconditional_dosage) {
    // case 2: unconditional dosage.  handle separately from other two cases
    // since missing values may be present.
    SetAllBits(raw_sample_ct, raw_dosage_present);
    raw_dosage_ct = raw_sample_ct;
  } else {
    // case 3: dosage bitarray
    raw_dosage_present[raw_sample_ctl - 1] = 0;
    const uint32_t raw_sample_ctb = DivUp(raw_sample_ct, CHAR_BIT);
    memcpy(raw_dosage_present, fread_ptr, raw_sample_ctb);
    fread_ptr = &(fread_ptr[raw_sample_ctb]);
    raw_dosage_ct = PopcountWords(raw_dosage_present, raw_sample_ctl);
  }
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  uint32_t dosage_ct;
  if (subsetting_required) {
    CopyBitarrSubset(raw_dosage_present, sample_include, sample_ct, dosage_present);
    dosage_ct = PopcountWords(dosage_present, sample_ctl);
  } else {
    dosage_ct = raw_dosage_ct;
  }
  if (dosage_ct_ptr) {
    *dosage_ct_ptr = dosage_ct;
  }
  if (!dosage_ct) {
    return kPglRetSuccess;
  }
#ifdef __arm__
#  error "Unaligned accesses in ParseDosage16()."
#endif
  const uint16_t* dosage_vals_read_iter = R_CAST(const uint16_t*, fread_ptr);
  uint16_t* dosage_vals_write_iter = dosage_vals;
  if (!(vrtype & 0x80)) {
    if (alt_allele_ct == 1) {
      if (dosage_ct == raw_dosage_ct) {
        if (!is_unconditional_dosage) {
          memcpy(dosage_vals_write_iter, dosage_vals_read_iter, dosage_ct * sizeof(int16_t));
        } else {
          for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
            const uint16_t cur_dosage = *dosage_vals_read_iter++;
            if (cur_dosage != 65535) {
              *dosage_vals_write_iter++ = cur_dosage;
            } else {
              ClearBit(sample_idx, raw_dosage_present);
            }
          }
          if (dosage_ct_ptr) {
            *dosage_ct_ptr = dosage_vals_write_iter - dosage_vals;
          }
        }
      } else {
        uint32_t sample_uidx = 0;
        // bugfix (22 May 2017): dosage_entry_idx needs to iterate up to
        // raw_dosage_ct, not dosage_ct
        for (uint32_t dosage_entry_idx = 0; dosage_entry_idx < raw_dosage_ct; ++dosage_entry_idx, ++sample_uidx, ++dosage_vals_read_iter) {
          FindFirst1BitFromU32(raw_dosage_present, &sample_uidx);
          if (!IsSet(sample_include, sample_uidx)) {
            continue;
          }
          *dosage_vals_write_iter++ = *dosage_vals_read_iter;
        }
      }
    } else {
      // todo: multiallelic dosage
      // need to support downcode to ref/alt1 as well as raw load
      // (dosage_ct_ptr should be nullptr iff we're doing a raw load)
      return kPglRetNotYetSupported;
    }
    return kPglRetSuccess;
  } else {
    // todo: phased dosage
    return kPglRetNotYetSupported;
  }
}

PglErr PgrReadRefalt1GenovecDosage16SubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr, uint32_t* is_explicit_alt1_ptr) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    *dosage_ct_ptr = 0;
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  if ((!VrtypeDosage(vrtype)) || (!dosage_present)) {
    PglErr reterr = ReadRefalt1GenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genovec);
    *dosage_ct_ptr = 0;
    return reterr;
  }
  const unsigned char* fread_ptr = nullptr;
  const unsigned char* fread_end = nullptr;
  uint32_t phasepresent_ct;
  PglErr reterr = ReadRefalt1GenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, &fread_ptr, &fread_end, genovec, nullptr, nullptr, &phasepresent_ct);
  if (reterr) {
    return reterr;
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t alt_allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx] - 1) : 1;
  *is_explicit_alt1_ptr = (alt_allele_ct > 1);
  return ParseDosage16(fread_ptr, fread_end, sample_include, sample_ct, vidx, alt_allele_ct, pgrp, dosage_ct_ptr, dosage_present, dosage_vals);
}

uint64_t U16VecSum(const uint16_t* __restrict uint16_vec, uint32_t entry_ct) {
#ifdef __LP64__
  // UniVecHsum32() could overflow once we exceed this
  const uint32_t max_loop_len = (131072 / kInt32PerVec) - 1;

  const VecUL m16 = VCONST_UL(kMask0000FFFF);
  const VecUL* uint16_vvec_iter = R_CAST(const VecUL*, uint16_vec);
  uint32_t full_vecs_remaining = entry_ct / (kBytesPerVec / sizeof(int16_t));
  uint64_t sum = 0;
  while (1) {
    UniVec acc_even;
    UniVec acc_odd;
    acc_even.vi = vecul_setzero();
    acc_odd.vi = vecul_setzero();
    const VecUL* uint16_vvec_stop;
    if (full_vecs_remaining < max_loop_len) {
      if (!full_vecs_remaining) {
        const uint32_t trail_ct = entry_ct % (kBytesPerVec / sizeof(int16_t));
        uint16_vec = R_CAST(const uint16_t*, uint16_vvec_iter);
        for (uint32_t uii = 0; uii < trail_ct; ++uii) {
          sum += uint16_vec[uii];
        }
        return sum;
      }
      uint16_vvec_stop = &(uint16_vvec_iter[full_vecs_remaining]);
      full_vecs_remaining = 0;
    } else {
      uint16_vvec_stop = &(uint16_vvec_iter[max_loop_len]);
      full_vecs_remaining -= max_loop_len;
    }
    do {
      const VecUL cur_vec = *uint16_vvec_iter++;
      acc_even.vi = acc_even.vi + (cur_vec & m16);
      acc_odd.vi = acc_odd.vi + (vecul_srli(cur_vec, 16) & m16);
    } while (uint16_vvec_iter < uint16_vvec_stop);
    sum += UniVecHsum32(acc_even);
    sum += UniVecHsum32(acc_odd);
  }
#else
  uint64_t sum = 0;
  for (uint32_t uii = 0; uii < entry_ct; ++uii) {
    sum += uint16_vec[uii];
  }
  return sum;
#endif
}

PglErr GetRefNonrefGenotypeCountsAndDosage16s(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, double* mach_r2_ptr, uint32_t* genocounts, uint64_t* all_dosages) {
  // genocounts[0] := ref/ref, genocounts[1] := ref/altx,
  // genocounts[2] := altx/alty, genocounts[3] := missing
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  // to avoid LD cache thrashing, we either always keep a subsetted cache, or
  // never do so.
  // todo: can't take the shortcut in the multiallelic variant case
  if ((!(pgrp->fi.gflags & (kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent))) || ((!(vrtype & 0x68)) && (!subsetting_required))) {
    {
      PglErr reterr = GetRefalt1GenotypeCounts(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genocounts);
      if (reterr) {
        return reterr;
      }
    }
  GetRefNonrefGenotypeCountsAndDosage16s_basic_finish:
    all_dosages[0] = (genocounts[0] * 2 + genocounts[1]) * 16384LLU;
    all_dosages[1] = (genocounts[2] * 2 + genocounts[1]) * 16384LLU;
    if (!mach_r2_ptr) {
      return kPglRetSuccess;
    }
    // yeah, it's sinful to implement mach-r2 here...
    const uint32_t nm_sample_ct = sample_ct - genocounts[3];
    double mach_r2 = 1.0;
    if (nm_sample_ct) {
      const uintptr_t dosage_sum = genocounts[2] * 2 + genocounts[1];
      const int64_t dosage_ssq = dosage_sum + genocounts[2] * 2LLU;
      const double dosage_sumd = dosage_sum;
      const double dosage_avg = dosage_sumd / u31tod(nm_sample_ct);
      const double dosage_variance = dosage_ssq - dosage_sumd * dosage_avg;
      mach_r2 = 2 * dosage_variance / (dosage_sumd * (2 - dosage_avg));
    }
    *mach_r2_ptr = mach_r2;
    return kPglRetSuccess;
  }
  uintptr_t* tmp_genovec = pgrp->workspace_vec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRefalt1GenovecSubsetUnsafe(nullptr, nullptr, raw_sample_ct, vidx, pgrp, &fread_ptr, &fread_end, tmp_genovec);
  if (reterr) {
    return reterr;
  }
  if (!subsetting_required) {
    ZeroTrailingQuaters(raw_sample_ct, tmp_genovec);
    GenovecCountFreqsUnsafe(tmp_genovec, raw_sample_ct, genocounts);
  } else {
    GenovecCountSubsetFreqs(tmp_genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
  }
  if (!(vrtype & 0x68)) {
    goto GetRefNonrefGenotypeCountsAndDosage16s_basic_finish;
  }
  if (vrtype & 8) {
    // todo: multiallelic case
    assert(0);
    if (!(vrtype & 0x60)) {
      return kPglRetSuccess;
    }
    // update raw_het_ct if hphase present
  }
  if (vrtype & 0x10) {
    uint32_t raw_het_ct;
    if (!subsetting_required) {
      raw_het_ct = genocounts[1];
    } else {
      ZeroTrailingQuaters(raw_sample_ct, tmp_genovec);
      raw_het_ct = GenovecCount01Unsafe(tmp_genovec, raw_sample_ct);
    }
    // skip phase info
    // probably make this its own function...
    // bugfix: need to use raw het ct, not subsetted
    const uint32_t first_half_byte_ct = 1 + (raw_het_ct / CHAR_BIT);
    const uint32_t explicit_phasepresent = fread_ptr[0] & 1;
    if (explicit_phasepresent) {
      const uint32_t raw_phasepresent_ct = PopcountBytes(fread_ptr, first_half_byte_ct) - 1;
      const uint32_t second_half_byte_ct = DivUp(raw_phasepresent_ct, CHAR_BIT);
      fread_ptr = &(fread_ptr[first_half_byte_ct + second_half_byte_ct]);
    } else {
      fread_ptr = &(fread_ptr[first_half_byte_ct]);
    }
  }

  // todo: phased dosage
#ifndef NDEBUG
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t alt_allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx] - 1) : 1;
  assert(alt_allele_ct == 1);
#endif
  uint64_t alt1_dosage = 0;
  uint64_t alt1_dosage_sq_sum = 0;
  uint32_t dosage_ct = 0;
  uint32_t replaced_genocounts[4];
  if ((vrtype & 0x60) == 0x40) {
    // unconditional dosage.  needs to be handled separately from the other
    // cases due to possible presence of missing values.
    // note that this code will also need to be adjusted when multiallelic
    // support is added.
    uint32_t sample_uidx = 0;
#ifdef __arm__
#  error "Unaligned accesses in GetRefNonrefGenotypeCountsAndDosage16s()."
#endif
    ZeroUiArr(4, replaced_genocounts);
    const uint16_t* dosage_vals = R_CAST(const uint16_t*, fread_ptr);
    if (subsetting_required) {
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
        FindFirst1BitFromU32(sample_include, &sample_uidx);
        const uintptr_t cur_dosage_val = dosage_vals[sample_uidx];
        if (cur_dosage_val != 65535) {
          alt1_dosage += cur_dosage_val;

          // todo: check if this is slow enough to justify removing it from the
          // main loop
          alt1_dosage_sq_sum += cur_dosage_val * cur_dosage_val;

          ++dosage_ct;
          const uint32_t hardcall_code = GetQuaterarrEntry(tmp_genovec, sample_uidx);
          replaced_genocounts[hardcall_code] += 1;
        }
      }
    } else {
      for (; sample_uidx < sample_ct; ++sample_uidx) {
        const uintptr_t cur_dosage_val = dosage_vals[sample_uidx];
        if (cur_dosage_val != 65535) {
          alt1_dosage += cur_dosage_val;
          alt1_dosage_sq_sum += cur_dosage_val * cur_dosage_val;
          ++dosage_ct;
          const uint32_t hardcall_code = GetQuaterarrEntry(tmp_genovec, sample_uidx);
          replaced_genocounts[hardcall_code] += 1;
        }
      }
    }
  } else {
    uintptr_t* raw_dosage_present = pgrp->workspace_dosage_present;
    uint32_t raw_dosage_ct;
    if (!(vrtype & 0x40)) {
      // dosage list
      if (ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, &fread_ptr, raw_dosage_present, &raw_dosage_ct)) {
        return kPglRetMalformedInput;
      }
    } else {
      // dosage bitarray
      const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
      raw_dosage_present[raw_sample_ctl - 1] = 0;
      const uint32_t raw_sample_ctb = DivUp(raw_sample_ct, CHAR_BIT);
      memcpy(raw_dosage_present, fread_ptr, raw_sample_ctb);
      fread_ptr = &(fread_ptr[raw_sample_ctb]);
      raw_dosage_ct = PopcountWords(raw_dosage_present, raw_sample_ctl);
    }
    const uint16_t* dosage_vals_iter = R_CAST(const uint16_t*, fread_ptr);
    uint32_t sample_uidx = 0;
    if (subsetting_required) {
      for (uint32_t dosage_idx = 0; dosage_idx < raw_dosage_ct; ++dosage_idx, ++sample_uidx) {
        FindFirst1BitFromU32(raw_dosage_present, &sample_uidx);
        if (IsSet(sample_include, sample_uidx)) {
          const uintptr_t cur_dosage_val = dosage_vals_iter[dosage_idx];
          alt1_dosage += cur_dosage_val;
          alt1_dosage_sq_sum += cur_dosage_val * cur_dosage_val;
          ++dosage_ct;
        }
      }
      GenoarrCountSubsetIntersectFreqs(tmp_genovec, raw_dosage_present, sample_include, raw_sample_ct, replaced_genocounts);
    } else {
      if (!mach_r2_ptr) {
        for (uint32_t dosage_idx = 0; dosage_idx < raw_dosage_ct; ++dosage_idx) {
          alt1_dosage += dosage_vals_iter[dosage_idx];
        }
      } else {
        for (uint32_t dosage_idx = 0; dosage_idx < raw_dosage_ct; ++dosage_idx) {
          const uintptr_t cur_dosage_val = dosage_vals_iter[dosage_idx];
          alt1_dosage += cur_dosage_val;
          alt1_dosage_sq_sum += cur_dosage_val * cur_dosage_val;
        }
      }
      dosage_ct = raw_dosage_ct;
      GenoarrCountSubsetFreqs2(tmp_genovec, raw_dosage_present, raw_sample_ct, raw_dosage_ct, replaced_genocounts);
    }
  }
  const uint32_t replaced_ct = replaced_genocounts[0] + replaced_genocounts[1] + replaced_genocounts[2];
  const uint32_t remaining_het_ct = genocounts[1] - replaced_genocounts[1];
  const uint32_t remaining_hom_alt_ct = genocounts[2] - replaced_genocounts[2];
  const uint32_t alt1_ct = 2 * remaining_hom_alt_ct + remaining_het_ct;
  alt1_dosage += alt1_ct * 16384LLU;
  all_dosages[1] = alt1_dosage;
  const uint32_t nondosage_nm_ct = sample_ct - genocounts[3] - replaced_ct;
  const uint32_t new_sample_nm_ct = dosage_ct + nondosage_nm_ct;
  all_dosages[0] = new_sample_nm_ct * 32768LLU - alt1_dosage;
  if (!mach_r2_ptr) {
    return kPglRetSuccess;
  }
  double mach_r2 = 1.0;
  if (new_sample_nm_ct) {
    // 16384^2, 32768^2
    alt1_dosage_sq_sum += remaining_het_ct * 0x10000000LLU + remaining_hom_alt_ct * 0x40000000LLU;
    const double dosage_sumd = u63tod(alt1_dosage);
    const double dosage_avg = dosage_sumd / u31tod(new_sample_nm_ct);
    const double dosage_variance = u63tod(alt1_dosage_sq_sum) - dosage_sumd * dosage_avg;
    mach_r2 = 2 * dosage_variance / (dosage_sumd * (32768 - dosage_avg));
  }
  *mach_r2_ptr = mach_r2;
  return kPglRetSuccess;
}

PglErr PgrGetRefNonrefGenotypeCountsAndDosage16s(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, double* mach_r2_ptr, uint32_t* genocounts, uint64_t* all_dosages) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    ZeroUiArr(4, genocounts);
    const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
    const uint32_t cur_allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
    ZeroU64Arr(cur_allele_ct, all_dosages);
    if (mach_r2_ptr) {
      *mach_r2_ptr = 1.0;
    }
    return kPglRetSuccess;
  }
  return GetRefNonrefGenotypeCountsAndDosage16s(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, mach_r2_ptr, genocounts, all_dosages);
}

PglErr PgrReadRefalt1GenovecHphaseDosage16SubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* phasepresent_ct_ptr, uintptr_t* __restrict dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr, uint32_t* is_explicit_alt1_ptr) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    *phasepresent_ct_ptr = 0;
    *dosage_ct_ptr = 0;
    return kPglRetSuccess;
  }
  const unsigned char* fread_ptr = nullptr;
  const unsigned char* fread_end = nullptr;
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t dosage_is_present = VrtypeDosage(vrtype);
  PglErr reterr = ReadRefalt1GenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, dosage_is_present? (&fread_ptr) : nullptr, dosage_is_present? (&fread_end) : nullptr, genovec, phasepresent, phaseinfo, phasepresent_ct_ptr);
  if (reterr || (!dosage_is_present)) {
    *dosage_ct_ptr = 0;
    return reterr;
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t alt_allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx] - 1) : 1;
  *is_explicit_alt1_ptr = (alt_allele_ct > 1);
  return ParseDosage16(fread_ptr, fread_end, sample_include, sample_ct, vidx, alt_allele_ct, pgrp, dosage_ct_ptr, dosage_present, dosage_vals);
}

PglErr PgrReadRaw(uint32_t vidx, PgenGlobalFlags read_gflags, PgenReader* pgrp, uintptr_t** loadbuf_iter_ptr, unsigned char* loaded_vrtype_ptr) {
  // currently handles hardcall phase and unphased dosage
  // todo: multiallelic variants, phased dosage
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  uintptr_t* genovec = (*loadbuf_iter_ptr);
  uintptr_t* loadbuf_iter = &(genovec[QuaterCtToAlignedWordCt(raw_sample_ct)]);
  const uint32_t hphase_is_present = (vrtype / 0x10) & 1;
  const uint32_t save_hphase = hphase_is_present && (read_gflags & kfPgenGlobalHardcallPhasePresent);
  const uint32_t dosage_is_present = (vrtype & 0x60)? 1 : 0;
  const uint32_t save_dosage = dosage_is_present && (read_gflags & kfPgenGlobalDosagePresent);
  if (loaded_vrtype_ptr) {
    *loaded_vrtype_ptr = save_hphase * 0x10 + save_dosage * 0x60;
  }
  if (!(save_hphase || save_dosage)) {
    // don't bother updating ldbase_all_hets, too much of a performance
    // penalty, and too likely that we won't need it
    *loadbuf_iter_ptr = loadbuf_iter;
    return ReadRefalt1GenovecSubsetUnsafe(nullptr, nullptr, raw_sample_ct, vidx, pgrp, nullptr, nullptr, genovec);
  }

  // todo: main multiallelic track goes here
  // (i) nonmissing bitarray
  // (ii) appropriate-width allele codes

  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t multiallelic_relevant = VrtypeMultiallelic(vrtype);
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  uintptr_t* all_hets = hphase_is_present? pgrp->workspace_all_hets : nullptr;
  if (VrtypeLdCompressed(vrtype)) {
    // ldbase_all_hets not needed in this case
    PglErr reterr = LdLoadGenovecSubsetIfNecessary(nullptr, nullptr, raw_sample_ct, vidx, pgrp);
    if (reterr) {
      return reterr;
    }
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    CopyQuaterarr(pgrp->ldbase_genovec, raw_sample_ct, genovec);

    reterr = ParseAndApplyDifflist(fread_end, multiallelic_relevant, &fread_ptr, pgrp, genovec);
    if (reterr) {
      return reterr;
    }
    if (all_hets) {
      PgrDetectGenovecHets(genovec, raw_sample_ct, all_hets);
    }
    if ((vrtype & 7) == 3) {
      GenovecInvertUnsafe(raw_sample_ct, genovec);
    }
  } else {
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      return kPglRetReadFail;
    }
    PglErr reterr = ParseNonLdGenovecHphaseSubset(fread_end, nullptr, nullptr, raw_sample_ct, vrtype, multiallelic_relevant, &fread_ptr, pgrp, genovec, all_hets);
    if (reterr) {
      return reterr;
    }
    const uint32_t is_ldbase = pgrp->fi.vrtypes && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
    if (is_ldbase) {
      CopyQuaterarr(genovec, raw_sample_ct, pgrp->ldbase_genovec);
      pgrp->ldbase_vidx = vidx;
      pgrp->ldbase_stypes = kfPgrLdcacheQuater;
    }
  }
  if (multiallelic_relevant) {
    // todo
    return kPglRetNotYetSupported;
  }
  if (all_hets) {
    const uint32_t het_ct = PopcountWords(all_hets, raw_sample_ctl);
    if (!het_ct) {
      // there shouldn't be a hphase track at all in this case
      return kPglRetMalformedInput;
    }
    const uint32_t het_ctdl = het_ct / kBitsPerWord;
    uintptr_t* phaseraw = loadbuf_iter;
    const uint32_t first_half_byte_ct = 1 + (het_ct / CHAR_BIT);
    if (save_hphase) {
      // this needs to be synced with phaseraw_word_ct in make_pgen_thread()
      loadbuf_iter = &(loadbuf_iter[kWordsPerVec + RoundDownPow2(raw_sample_ct / kBitsPerWordD2, kWordsPerVec)]);
      phaseraw[het_ctdl] = 0;

      memcpy(phaseraw, fread_ptr, first_half_byte_ct);
    }
    const uint32_t explicit_phasepresent = fread_ptr[0] & 1;
    fread_ptr = &(fread_ptr[first_half_byte_ct]);
    if (explicit_phasepresent) {
      const uint32_t raw_phasepresent_ct = PopcountWords(phaseraw, het_ctdl + 1) - 1;
      if (!raw_phasepresent_ct) {
        // there shouldn't be a hphase track at all in this case, either
        return kPglRetMalformedInput;
      }
      const uint32_t second_half_byte_ct = DivUp(raw_phasepresent_ct, CHAR_BIT);
      if (save_hphase) {
        // put this in a phasepresent-independent location, to make things more
        // convenient for the caller
        memcpy(&(phaseraw[1 + (raw_sample_ct / kBitsPerWord)]), fread_ptr, second_half_byte_ct);
      }
      fread_ptr = &(fread_ptr[second_half_byte_ct]);
    }
  }
  if (!save_dosage) {
    *loadbuf_iter_ptr = loadbuf_iter;
    return kPglRetSuccess;
  }
  uintptr_t* dosage_present = loadbuf_iter;
  loadbuf_iter = &(loadbuf_iter[BitCtToAlignedWordCt(raw_sample_ct)]);
  uint16_t* dosage_vals = R_CAST(uint16_t*, loadbuf_iter);
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t alt_allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx] - 1) : 1;
  *loadbuf_iter_ptr = &(loadbuf_iter[kWordsPerVec * DivUp(raw_sample_ct, (kBytesPerVec / sizeof(int16_t)))]);
  return ParseDosage16(fread_ptr, fread_end, nullptr, raw_sample_ct, vidx, alt_allele_ct, pgrp, nullptr, dosage_present, dosage_vals);
}


// tried to have more custom code, turned out to not be worth it
PglErr ReadMissingness(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict missingness, uintptr_t* __restrict hets, uintptr_t* __restrict genovec_buf) {
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRefalt1GenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, &fread_ptr, &fread_end, genovec_buf);
  ZeroTrailingQuaters(sample_ct, genovec_buf);
  GenovecToMissingnessUnsafe(genovec_buf, sample_ct, missingness);
  if (hets) {
    PgrDetectGenovecHetsUnsafe(genovec_buf, QuaterCtToWordCt(sample_ct), hets);
  }
  if (fread_pp) {
    *fread_pp = fread_ptr;
    *fread_endp = fread_end;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t is_multiallelic = VrtypeMultiallelic(vrtype);
  if (reterr || (!is_multiallelic)) {
    return reterr;
  }
  // todo: multiallelic case
  assert(0);
  return kPglRetSuccess;
}

PglErr PgrReadMissingness(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict missingness, uintptr_t* __restrict genovec_buf) {
  // may as well add a hets parameter?
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  return ReadMissingness(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, missingness, nullptr, genovec_buf);
}

/*
PglErr PgrReadMissingnessDosage(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict missingness, uintptr_t* __restrict genovec_buf) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t dosage_is_present = VrtypeDosage(vrtype);
  const uint32_t need_to_skip_hphase = dosage_is_present && VrtypeHphase(vrtype);
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  const unsigned char* fread_ptr = nullptr;
  const unsigned char* fread_end = nullptr;
  if (!need_to_skip_hphase) {
    PglErr reterr = ReadMissingness(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, dosage_is_present? (&fread_ptr) : nullptr, dosage_is_present? (&fread_end) : nullptr, missingness, genovec_buf);
    if (reterr || (!dosage_is_present)) {
      return reterr;
    }
  } else {
    uint32_t dummy;
    PglErr reterr = ReadRefalt1GenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, &fread_ptr, &fread_end, genovec_buf, nullptr, nullptr, &dummy);
    if (reterr) {
      return reterr;
    }
    ZeroTrailingQuaters(sample_ct, genovec_buf);
    GenovecToMissingnessUnsafe(genovec_buf, sample_ct, missingness);
  }
  // now perform bitwise andnot with dosage_present
  if ((vrtype & 0x60) == 0x40) {
    // unconditional dosage.  spot-check the appropriate entries for equality
    // to 65535.
#ifdef __arm__
#  error "Unaligned accesses in PgrReadMissingnessDosage()."
#endif
    const uint16_t* dosage_vals = (const uint16_t*)fread_ptr;
    uint32_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      FindFirst1BitFromU32(sample_include, &sample_uidx);
      if (!IsSet(missingness, sample_idx)) {
        continue;
      }
      if (dosage_vals[sample_uidx] != 65535) {
        ClearBit(sample_idx, missingness);
      }
    }
    return kPglRetSuccess;
  }
  uintptr_t* dosage_present = pgrp->workspace_dosage_present;
  if ((vrtype & 0x60) == 0x20) {
    // dosage list
    uint32_t dummy;
    if (ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, &fread_ptr, dosage_present, &dummy)) {
      return kPglRetMalformedInput;
    }
  } else {
    // dosage bitarray
    dosage_present[raw_sample_ctl - 1] = 0;
    const uint32_t raw_sample_ctb = DivUp(raw_sample_ct, CHAR_BIT);
    memcpy(dosage_present, fread_ptr, raw_sample_ctb);
  }
  if (subsetting_required) {
    CopyBitarrSubset(dosage_present, sample_include, sample_ct, pgrp->workspace_vec);
    dosage_present = pgrp->workspace_vec;
  }
  BitvecAndNot(dosage_present, BitCtToWordCt(sample_ct), missingness);
  return kPglRetSuccess;
}
*/

PglErr PgrReadMissingnessMulti(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict missingness_hc, uintptr_t* __restrict missingness_dosage, uintptr_t* __restrict hets, uintptr_t* __restrict genovec_buf) {
  // either missingness_hc or missingness_dosage must be non-null
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t dosage_is_relevant = missingness_dosage && VrtypeDosage(vrtype);
  const uint32_t need_to_skip_hphase = dosage_is_relevant && VrtypeHphase(vrtype);
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  const unsigned char* fread_ptr = nullptr;
  const unsigned char* fread_end = nullptr;
  uintptr_t* missingness_base = missingness_hc? missingness_hc : missingness_dosage;
  if (!need_to_skip_hphase) {
    PglErr reterr = ReadMissingness(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, dosage_is_relevant? (&fread_ptr) : nullptr, dosage_is_relevant? (&fread_end) : nullptr, missingness_base, hets, genovec_buf);
    if (missingness_dosage && missingness_hc) {
      memcpy(missingness_dosage, missingness_hc, BitCtToWordCt(sample_ct) * sizeof(intptr_t));
    }
    if (reterr || (!dosage_is_relevant)) {
      return reterr;
    }
  } else {
    uint32_t dummy;
    // will need to switch to a different function when multiallelic variants
    // are implemented.
    PglErr reterr = ReadRefalt1GenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, &fread_ptr, &fread_end, genovec_buf, nullptr, nullptr, &dummy);
    if (reterr) {
      return reterr;
    }
    ZeroTrailingQuaters(sample_ct, genovec_buf);
    GenovecToMissingnessUnsafe(genovec_buf, sample_ct, missingness_base);
    if (hets) {
      PgrDetectGenovecHetsUnsafe(genovec_buf, QuaterCtToWordCt(sample_ct), hets);
    }
    if (missingness_hc) {
      memcpy(missingness_dosage, missingness_hc, BitCtToWordCt(sample_ct) * sizeof(intptr_t));
    }
  }
  // now perform bitwise andnot with dosage_present
  if ((vrtype & 0x60) == 0x40) {
    // unconditional dosage.  spot-check the appropriate entries for equality
    // to 65535.
#ifdef __arm__
#  error "Unaligned accesses in PgrReadMissingnessMulti()."
#endif
    const uint16_t* dosage_vals = R_CAST(const uint16_t*, fread_ptr);
    uint32_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      FindFirst1BitFromU32(sample_include, &sample_uidx);
      if (!IsSet(missingness_dosage, sample_idx)) {
        continue;
      }
      if (dosage_vals[sample_uidx] != 65535) {
        ClearBit(sample_idx, missingness_dosage);
      }
    }
    return kPglRetSuccess;
  }
  uintptr_t* dosage_present = pgrp->workspace_dosage_present;
  if ((vrtype & 0x60) == 0x20) {
    // dosage list
    uint32_t dummy;
    if (ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, &fread_ptr, dosage_present, &dummy)) {
      return kPglRetMalformedInput;
    }
  } else {
    // dosage bitarray
    dosage_present[raw_sample_ctl - 1] = 0;
    const uint32_t raw_sample_ctb = DivUp(raw_sample_ct, CHAR_BIT);
    memcpy(dosage_present, fread_ptr, raw_sample_ctb);
  }
  if (subsetting_required) {
    CopyBitarrSubset(dosage_present, sample_include, sample_ct, pgrp->workspace_vec);
    dosage_present = pgrp->workspace_vec;
  }
  BitvecAndNot(dosage_present, BitCtToWordCt(sample_ct), missingness_dosage);
  return kPglRetSuccess;
}

static inline BoolErr ValidateVint31(const unsigned char* buf_end, const unsigned char** bufpp, uint32_t* val_ptr) {
  if (buf_end <= (*bufpp)) {
    return 1;
  }
  uint32_t vint32 = *((*bufpp)++);
  if (vint32 <= 127) {
    *val_ptr = vint32;
    return 0;
  }
  vint32 &= 127;
  for (uint32_t shift = 7; shift < 28; shift += 7) {
    if (buf_end == (*bufpp)) {
      return 1;
    }
    uint32_t uii = *((*bufpp)++);
    vint32 |= (uii & 127) << shift;
    if (uii <= 127) {
      *val_ptr = vint32;
      return 0;
    }
  }
  if (buf_end == (*bufpp)) {
    return 1;
  }
  uint32_t uii = *((*bufpp)++);
  if (uii > 7) {
    return 1;
  }
  vint32 |= uii << 28;
  *val_ptr = vint32;
  return 0;
}

BoolErr ValidateDifflistHeader(const unsigned char* fread_end, uint32_t sample_ct, const unsigned char** fread_pp, uintptr_t* raregeno_buf, const unsigned char** difflist_group_info_ptr, uint32_t* difflist_len_ptr) {
  // can be used for deltalists: pass raregeno_buf == nullptr.
  if (ValidateVint31(fread_end, fread_pp, difflist_len_ptr)) {
    // todo: ensure fread_pp points to a problematic byte whenever a validate_
    // function returns an error, so the error message can provide an accurate
    // byte offset.
    return 1;
  }
  const uint32_t difflist_len = *difflist_len_ptr;
  *difflist_group_info_ptr = *fread_pp;
  if (!difflist_len) {
    return 0;
  }
  if (difflist_len > sample_ct / kPglMaxDifflistLenDivisor) {
    return 1;
  }
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(sample_ct);
  const uint32_t difflist_index_byte_ct = group_ct * (sample_id_byte_ct + 1) - 1;
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < difflist_index_byte_ct) {
    return 1;
  }
  *fread_pp += difflist_index_byte_ct;
  if (!raregeno_buf) {
    return 0;
  }
  const uint32_t raregeno_byte_ct = QuaterCtToByteCt(difflist_len);
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < raregeno_byte_ct) {
    return 1;
  }
  const unsigned char* raregeno_end = &((*fread_pp)[raregeno_byte_ct]);
  memcpy(raregeno_buf, *fread_pp, raregeno_byte_ct);
  *fread_pp = raregeno_end;
  const uint32_t difflist_len_mod4 = difflist_len % 4;
  if (difflist_len_mod4) {
    const uint32_t last_raregeno_byte = (*fread_pp)[-1];
    if (last_raregeno_byte >> (2 * difflist_len_mod4)) {
      return 1;
    }
  }
  return 0;
}

BoolErr ValidateAndApplyDifflist(const unsigned char* fread_end, uint32_t common2_code, uint32_t multiallelic_relevant, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effects: uses pgr.workspace_raregeno_tmp_loadbuf.
  // Similar to ParseAndApplyDifflist(), but with exhaustive input
  // validation.
  // If multiallelic_relevant is true, a list of sample indices with freshly
  // loaded raregeno value 0b11 is saved to pgr.workspace_ambig_sample_ids, and
  // pgr.workspace_ambig_id_ct is set to the length of the list.
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  uintptr_t* cur_raregeno_iter = pgrp->workspace_raregeno_tmp_loadbuf;
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  if (ValidateDifflistHeader(fread_end, sample_ct, fread_pp, cur_raregeno_iter, &group_info_iter, &difflist_len)) {
    return 1;
  }
  if (!difflist_len) {
    return 0;
  }
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  if (common2_code) {
    // 1-bit format + list of exceptions.  In this case,
    //   (i) the length of the exception list must be < (sample_ct / 16)
    //   (ii) every raregeno entry must either be one of the two rare genotype
    //        values, or involve a rare alt allele.
    if (difflist_len >= (sample_ct / (2 * kPglMaxDifflistLenDivisor))) {
      return 1;
    }
    const uintptr_t common_code_delta = common2_code & 3;
    const uintptr_t inv_common_word1 = (3 - common2_code / 4) * kMask5555;
    const uintptr_t inv_common_word2 = inv_common_word1 - (common_code_delta * kMask5555);
    uint32_t subgroup_idx = 0;
    while (1) {
      uintptr_t cur_raregeno_word = cur_raregeno_iter[subgroup_idx];
      uintptr_t match1 = cur_raregeno_word ^ inv_common_word1;
      match1 = match1 & (match1 >> 1) & kMask5555;
      uintptr_t match2 = cur_raregeno_word ^ inv_common_word2;
      match2 = match2 & (match2 >> 1) & kMask5555;
      if (subgroup_idx == subgroup_idx_last) {
        // ignore trailing bits
        const uint32_t lshift = ((-difflist_len) % kBitsPerWordD2) * 2;
        if ((match1 << lshift) || (match2 << lshift)) {
          return 1;
        }
        break;
      }
      if (match1 || match2) {
        // todo: if (multiallelic_relevant && (!inv_common_word2)), record
        // might be fine; but we need to verify these are actually rare alt
        // alleles.
        return 1;
      }
      ++subgroup_idx;
    }
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(sample_ct);
  uint32_t* ambig_sample_ids = multiallelic_relevant? pgrp->workspace_ambig_sample_ids : nullptr;
  const unsigned char* group_byte_cts_iter = &(group_info_iter[DivUp(difflist_len, kPglDifflistGroupSize) * sample_id_byte_ct]);
  const unsigned char* prev_group_start = *fread_pp;

  uintptr_t sample_idx = 0;
  uint32_t ambig_id_ct = 0;
  uint32_t subgroup_idx = 0;
  while (1) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        pgrp->workspace_ambig_id_ct = ambig_id_ct;
        return 0;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
      uintptr_t new_sample_idx_start = SubUintLoad(group_info_iter, sample_id_byte_ct);
      if (subgroup_idx) {
        if (sample_idx >= new_sample_idx_start) {
          return 1;
        }
        const uint32_t group_byte_ct = S_CAST(uint32_t, *group_byte_cts_iter++) + 63;
        if (S_CAST(uintptr_t, (*fread_pp) - prev_group_start) != group_byte_ct) {
          return 1;
        }
        prev_group_start = *fread_pp;
      }
      sample_idx = new_sample_idx_start;
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      uint32_t sample_idx_incr;
      if (ValidateVint31(fread_end, fread_pp, &sample_idx_incr) || (!sample_idx_incr)) {
        return 1;
      }
      sample_idx += sample_idx_incr;
    }
    ++subgroup_idx;
    uintptr_t cur_raregeno_word = *cur_raregeno_iter++;
    while (1) {
      if (sample_idx >= sample_ct) {
        return 1;
      }
      const uintptr_t cur_geno = cur_raregeno_word & 3;
      AssignQuaterarrEntry(sample_idx, cur_geno, genovec);
      if (multiallelic_relevant && (cur_geno == 3)) {
        ambig_sample_ids[ambig_id_ct++] = sample_idx;
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      uint32_t sample_idx_incr;
      if (ValidateVint31(fread_end, fread_pp, &sample_idx_incr) || (!sample_idx_incr)) {
        return 1;
      }
      sample_idx += sample_idx_incr;
      --remaining_deltas_in_subgroup;
      cur_raregeno_word >>= 2;
    }
  }
}

BoolErr ValidateOnebit(const unsigned char* fread_end, uint32_t difflist_ambig_ids_needed, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // ParseOnebitUnsafe() with exhaustive input validation.
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t common2_and_bitarray_byte_ct = (sample_ct + 15) / CHAR_BIT;
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < common2_and_bitarray_byte_ct) {
    return 1;
  }
  const unsigned char* fread_difflist_start = &((*fread_pp)[common2_and_bitarray_byte_ct]);
  const uintptr_t common2_code = *((*fread_pp)++);
  const uintptr_t common_code_delta = common2_code & 3;
  uintptr_t word_base = common2_code / 4;
  if ((!common_code_delta) || (word_base + common_code_delta > 3)) {
    return 1;
  }
  word_base *= kMask5555;
  const uint32_t genovec_widx_trail = (sample_ct + 7) / kBitsPerWordD2;
  const uint32_t genovec_widx_end = QuaterCtToWordCt(sample_ct);
  uint32_t genovec_widx = 0;
#ifdef __arm__
#  error "Unaligned accesses in ValidateOnebit()."
#endif
  const Halfword* fread_alias = R_CAST(const Halfword*, *fread_pp);
  while (1) {
    uintptr_t ww;
    if (genovec_widx >= genovec_widx_trail) {
      if (genovec_widx == genovec_widx_end) {
        break;
      }
      const uint32_t nontrail_byte_ct = ((sample_ct - 1) % kBitsPerWordD2) / CHAR_BIT;
      ww = ProperSubwordLoad(&(fread_alias[genovec_widx_trail]), 1 + nontrail_byte_ct);
      const uint32_t sample_ct_mod8 = sample_ct % 8;
      if (sample_ct_mod8) {
        if (ww >> (nontrail_byte_ct * 8 + sample_ct_mod8)) {
          return 1;
        }
      }
    } else {
      ww = fread_alias[genovec_widx];
    }
    ww = UnpackHalfwordToWord(ww);
    genovec[genovec_widx++] = word_base + ww * common_code_delta;
  }
  *fread_pp = fread_difflist_start;
  return ValidateAndApplyDifflist(fread_end, common2_code, difflist_ambig_ids_needed, fread_pp, pgrp, genovec);
}

// assumes that we aren't dealing with the trivial fixed-width case.
// saves main genotype array to workspace_vec.  does not zero out trailing
// bits.
BoolErr ValidateGeno(const unsigned char* fread_end, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, char* errstr_buf) {
  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t multiallelic_relevant = VrtypeMultiallelic(vrtype);
  uintptr_t* genovec = pgrp->workspace_vec;
  if (VrtypeLdCompressed(vrtype)) {
    CopyQuaterarr(pgrp->ldbase_genovec, sample_ct, genovec);
    if (ValidateAndApplyDifflist(fread_end, 0, multiallelic_relevant, fread_pp, pgrp, genovec)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid LD difflist for (0-based) variant #%u.\n", vidx);
      return 1;
    }
    if (vrtype & 1) {
      // do we actually need this?
      GenovecInvertUnsafe(sample_ct, genovec);
    }
    return 0;
  }
  const uint32_t is_ldbase = VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
  if (!VrtypeDifflist(vrtype)) {
    if (vrtype & 1) {
      if (ValidateOnebit(fread_end, multiallelic_relevant, fread_pp, pgrp, genovec)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid 1-bit genotype record for (0-based) variant #%u.\n", vidx);
        return 1;
      }
    } else {
      const uint32_t genovec_byte_ct = DivUp(sample_ct, 4);
      if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < genovec_byte_ct) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid 2-bit genotype record for (0-based) variant #%u\n", vidx);
        return 1;
      }
      memcpy(genovec, *fread_pp, genovec_byte_ct);
      *fread_pp += genovec_byte_ct;
      const uint32_t sample_ct_mod4 = sample_ct % 4;
      if (sample_ct_mod4) {
        const uint32_t last_geno_byte = (*fread_pp)[-1];
        if (last_geno_byte >> (2 * sample_ct_mod4)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Last genotype byte for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
          return 1;
        }
      }
      if (VrtypeMultiallelic(vrtype)) {
        ExtractGenoarrAmbigIds(genovec, sample_ct, pgrp->workspace_ambig_sample_ids, &(pgrp->workspace_ambig_id_ct));
      }
    }
  } else {
    const uint32_t vrtype_low2 = vrtype & 3;
    const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
    memset(genovec, vrtype_low2 * 0x55, vec_ct * kBytesPerVec);
    if (ValidateAndApplyDifflist(fread_end, 0, multiallelic_relevant, fread_pp, pgrp, genovec)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid genotype difflist for (0-based) variant #%u.\n", vidx);
      return 1;
    }
  }
  if (is_ldbase) {
    CopyQuaterarr(genovec, sample_ct, pgrp->ldbase_genovec);
  }
  return 0;
}

BoolErr ValidateHphase(const unsigned char* fread_end, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, char* errstr_buf) {
  const uintptr_t* all_hets = pgrp->workspace_all_hets;
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t het_ct = PopcountWords(all_hets, sample_ctl);
  if (!het_ct) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Hardcall phase track present for (0-based) variant #%u, but there were no heterozygous calls.\n", vidx);
    return 1;
  }
  const uint32_t aux2_first_part_byte_ct = 1 + (het_ct / CHAR_BIT);
  const unsigned char* aux2_first_part = *fread_pp;
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < aux2_first_part_byte_ct) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid hardcall phase track present for (0-based) variant #%u.\n", vidx);
    return 1;
  }
  *fread_pp += aux2_first_part_byte_ct;
  const uint32_t het_ct_p1_mod8 = (het_ct + 1) % CHAR_BIT;
  if (het_ct_p1_mod8) {
    // verify trailing bits are zero
    if ((*fread_pp)[-1] >> het_ct_p1_mod8) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Hardcall phase track for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
      return 1;
    }
  }
  if (!((*aux2_first_part) & 1)) {
    // phase always present, "first part" is only part
    return 0;
  }
  const uint32_t phasepresent_ct = PopcountBytes(aux2_first_part, aux2_first_part_byte_ct) - 1;
  if (!phasepresent_ct) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Hardcall phase track for (0-based) variant #%u does not have any actual phase information.\n", vidx);
    return 1;
  }
  const uint32_t phaseinfo_byte_ct = DivUp(phasepresent_ct, CHAR_BIT);
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < phaseinfo_byte_ct) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid hardcall phase track present for (0-based) variant #%u.\n", vidx);
    return 1;
  }
  *fread_pp += phaseinfo_byte_ct;
  const uint32_t phasepresent_ct_mod8 = phasepresent_ct % 8;
  if (phasepresent_ct_mod8) {
    if ((*fread_pp)[-1] >> phasepresent_ct_mod8) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Hardcall phase track for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
      return 1;
    }
  }
  return 0;
}

BoolErr ValidateAndCountDeltalist(const unsigned char* fread_end, uint32_t sample_ct, const unsigned char** fread_pp, uint32_t* deltalist_len_ptr) {
  // we only need to know the number of entries in the list, not the actual bit
  // positions for now
  // (if we do need the bit positions, copy
  // ParseAndSaveDeltalistAsBitarr().)
  const unsigned char* group_info_iter;
  if (ValidateDifflistHeader(fread_end, sample_ct, fread_pp, nullptr, &group_info_iter, deltalist_len_ptr)) {
    return 1;
  }
  const uint32_t deltalist_len = *deltalist_len_ptr;
  if (!deltalist_len) {
    return 0;
  }
  // not an appropriate error, since this is just a tuning parameter for the
  // compressor; readers are expected to handle at least
  // (sample_ct / kPglMaxDifflistLenDivisor) entries.
  // if (deltalist_len > sample_ct / kPglMaxDeltalistLenDivisor) {
  //   return 1;
  // }
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(sample_ct);
  const uint32_t group_idx_last = (deltalist_len - 1) / kPglDifflistGroupSize;
  const unsigned char* group_byte_cts_iter = &(group_info_iter[DivUp(deltalist_len, kPglDifflistGroupSize) * sample_id_byte_ct]);
  const unsigned char* prev_group_start = *fread_pp;
  uint32_t group_len_m1 = kPglDifflistGroupSize - 1;
  uint32_t group_idx = 0;
  uintptr_t sample_idx = 0;
  while (1) {
    if (group_idx >= group_idx_last) {
      if (group_idx > group_idx_last) {
        return 0;
      }
      group_len_m1 &= deltalist_len - 1;
    }
    uintptr_t new_sample_idx = SubUintLoad(group_info_iter, sample_id_byte_ct);
    if (group_idx) {
      if (sample_idx >= new_sample_idx) {
        return 1;
      }
      const uint32_t group_byte_ct = S_CAST(uint32_t, *group_byte_cts_iter++) + 63;
      if (S_CAST(uintptr_t, (*fread_pp) - prev_group_start) != group_byte_ct) {
        return 1;
      }
      prev_group_start = *fread_pp;
    }
    sample_idx = new_sample_idx;
    group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    ++group_idx;
    uint32_t deltalist_idx_lowbits = 0;
    while (1) {
      if (sample_idx >= sample_ct) {
        return 1;
      }
      if (deltalist_idx_lowbits == group_len_m1) {
        break;
      }
      ++deltalist_idx_lowbits;
      uint32_t sample_idx_incr;
      if (ValidateVint31(fread_end, fread_pp, &sample_idx_incr) || (!sample_idx_incr)) {
        return 1;
      }
      sample_idx += sample_idx_incr;
    }
  }
}

PglErr ValidateDosage16(const unsigned char* fread_end, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, char* errstr_buf) {
  // similar to ParseDosage16().  doesn't support multiallelic data yet.
  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  if (vrtype & 0x80) {
    // this should be trivial, just multiply array lengths by 2...
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Phased dosage validation is not implemented yet.\n");
    return kPglRetNotYetSupported;
  }
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  if ((vrtype & 0x60) == 0x40) {
    // unconditional dosage.  handle separately from other two cases since
    // 65535 is valid.
#ifdef __arm__
#  error "Unaligned accesses in ValidateDosage16()."
#endif
    if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < sample_ct * sizeof(int16_t)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid unconditional dosage track for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
    const uint16_t* dosage_vals_read_iter = R_CAST(const uint16_t*, *fread_pp);
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      uint16_t cur_dosage_val_p1 = *dosage_vals_read_iter++;
      cur_dosage_val_p1 += 1;  // intentional overflow on 65535
      if (cur_dosage_val_p1 > 32769) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid unconditional dosage track for (0-based) variant #%u (dosage is greater than 2).\n", vidx);
        return kPglRetMalformedInput;
      }
    }
    *fread_pp += sample_ct * sizeof(int16_t);
    return kPglRetSuccess;
  }
  uint32_t dosage_ct;
  if ((vrtype & 0x60) == 0x20) {
    // dosage list
    if (ValidateAndCountDeltalist(fread_end, sample_ct, fread_pp, &dosage_ct)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid dosage list for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
  } else {
    const uint32_t sample_ctb = DivUp(sample_ct, CHAR_BIT);
    if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < sample_ctb) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid dosage subset for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
    dosage_ct = PopcountBytes(*fread_pp, sample_ctb);
    *fread_pp += sample_ctb;
    const uint32_t sample_ct_mod8 = sample_ct % 8;
    if (sample_ct_mod8) {
      if ((*fread_pp)[-1] >> sample_ct_mod8) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Dosage subset bitarray for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
        return kPglRetMalformedInput;
      }
    }
  }
  if (S_CAST(uintptr_t, fread_end - (*fread_pp)) < dosage_ct * sizeof(int16_t)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid dosage track for (0-based) variant #%u.\n", vidx);
    return kPglRetMalformedInput;
  }
  const uint16_t* dosage_vals_read_iter = R_CAST(const uint16_t*, *fread_pp);
  for (uint32_t dosage_idx = 0; dosage_idx < dosage_ct; ++dosage_idx) {
    if ((*dosage_vals_read_iter++) > 32768) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid dosage track for (0-based) variant #%u (dosage is greater than 2).\n", vidx);
      return kPglRetMalformedInput;
    }
  }
  *fread_pp += dosage_ct * sizeof(int16_t);
  return kPglRetSuccess;
}

static_assert(kPglVblockSize == 65536, "PgrValidate() needs to have an error message updated.");
PglErr PgrValidate(PgenReader* pgrp, char* errstr_buf) {
  // Performs all validation which isn't done by pgfi_init_phase{1,2}() and
  // PgrInit().
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t variant_ct = pgrp->fi.raw_variant_ct;
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t const_vrtype = pgrp->fi.const_vrtype;
  if (const_vrtype != UINT32_MAX) {
    if (allele_idx_offsets && (allele_idx_offsets[variant_ct] != 2 * variant_ct)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pvar file contains multiallelic variant(s), but .%s file does not.\n", (const_vrtype == kPglVrtypePlink1)? "bed" : "pgen");
      return kPglRetInconsistentInput;
    }
    // const uintptr_t const_vrec_width = pgrp->fi.const_vrec_width;
    if ((!const_vrtype) || (const_vrtype == kPglVrtypePlink1)) {
      // only thing that can go wrong is nonzero trailing bits
      const uint32_t dbl_sample_ct_mod4 = 2 * (sample_ct % 4);
      if (!dbl_sample_ct_mod4) {
        return kPglRetSuccess;
      }
      for (uint32_t vidx = 0; vidx < variant_ct; ++vidx) {
        const unsigned char* fread_ptr;
        const unsigned char* fread_end = nullptr;
        if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
          return kPglRetReadFail;
        }
        const uint32_t last_byte_in_record = fread_end[-1];
        if (last_byte_in_record >> dbl_sample_ct_mod4) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Last byte of (0-based) variant #%u has nonzero trailing bits.\n", vidx);
          return kPglRetMalformedInput;
        }
      }
      return kPglRetSuccess;
    }
    // todo: 16-bit dosage entries can't be in [32769,65534]
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Validation of fixed-width dosage formats is not implemented yet.\n");
    return kPglRetNotYetSupported;
  }
  const unsigned char* vrtypes = pgrp->fi.vrtypes;
  for (uint32_t vidx = 0; vidx < variant_ct; vidx += kPglVblockSize) {
    if (VrtypeLdCompressed(vrtypes[vidx])) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: (0-based) variant #%u is LD-compressed; this is prohibited when the variant index is a multiple of 65536.\n", vidx);
      return kPglRetMalformedInput;
    }
  }
  // file size may not be validated yet.
  uint64_t fsize;
  FILE* ff = pgrp->ff;
#ifndef NO_MMAP
  if (ff == nullptr) {
    // mmap case
    fsize = pgrp->fi.file_size;
  } else {
#endif
    if (fseeko(ff, 0, SEEK_END)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
    fsize = ftello(ff);
    pgrp->fp_vidx = 1;  // force fseek when loading first variant
#ifndef NO_MMAP
  }
#endif
  // todo: modify this check when phase sets are implemented
  const uint64_t expected_fsize = pgrp->fi.var_fpos[variant_ct];
  if (expected_fsize != fsize) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen header indicates that file size should be %" PRIu64 " bytes, but actual file size is %" PRIu64 " bytes.\n", expected_fsize, fsize);
    return kPglRetMalformedInput;
  }
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  uint32_t header_ctrl = 0;
#ifndef NO_MMAP
  if (ff == nullptr) {
#  ifdef __arm__
#    error "Unaligned accesses in PgrValidate()."
#  endif
    memcpy(&header_ctrl, &(pgrp->fi.block_base[11]), 1);
    // validate the random-access index.
    const uint64_t* fpos_index = R_CAST(const uint64_t*, &(pgrp->fi.block_base[12]));
    for (uint32_t vblock_idx = 0; vblock_idx < vblock_ct; ++vblock_idx) {
      if (fpos_index[vblock_idx] != pgrp->fi.var_fpos[vblock_idx * kPglVblockSize]) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen header vblock-start index is inconsistent with variant record length index.\n");
        return kPglRetMalformedInput;
      }
    }
  } else {
#endif
    if (fseeko(ff, 11, SEEK_SET)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
    header_ctrl = getc_unlocked(ff);
    if (header_ctrl > 255) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
    for (uint32_t vblock_idx = 0; vblock_idx < vblock_ct; ++vblock_idx) {
      uint64_t vblock_start_fpos;
      if (!fread_unlocked(&vblock_start_fpos, sizeof(int64_t), 1, ff)) {
        return kPglRetReadFail;
      }
      if (vblock_start_fpos != pgrp->fi.var_fpos[vblock_idx * kPglVblockSize]) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen header vblock-start index is inconsistent with variant record length index.\n");
        return kPglRetMalformedInput;
      }
    }
#ifndef NO_MMAP
  }
#endif
  const uint32_t vrtype_and_fpos_storage = header_ctrl & 15;
  const uint32_t alt_allele_ct_byte_ct = (header_ctrl >> 4) & 3;
  const uint32_t nonref_flags_stored = ((header_ctrl >> 6) == 3);

  // does not include vrtypes yet
  uint64_t vblock_index_byte_ct = kPglVblockSize * (1 + (vrtype_and_fpos_storage & 3) + alt_allele_ct_byte_ct);
  if (nonref_flags_stored) {
    vblock_index_byte_ct += kPglVblockSize / CHAR_BIT;
  }
  uint64_t last_vrtype_byte_offset = 0;
  uint32_t trailing_shift = 4;
  if (vrtype_and_fpos_storage & 8) {
    vblock_index_byte_ct += kPglVblockSize >> (10 - vrtype_and_fpos_storage);
    if (vrtype_and_fpos_storage == 8) {
      const uint32_t variant_ct_mod4 = variant_ct % 4;
      if (variant_ct_mod4) {
        last_vrtype_byte_offset = 20 + (vblock_ct - 1) * (vblock_index_byte_ct + sizeof(int64_t)) + ((variant_ct % kPglVblockSize) / 4);
        trailing_shift = variant_ct_mod4 * 2;
      }
    } else {
      assert(vrtype_and_fpos_storage == 9);
      if (variant_ct % 2) {
        last_vrtype_byte_offset = 20 + (vblock_ct - 1) * (vblock_index_byte_ct + sizeof(int64_t)) + ((variant_ct % kPglVblockSize) / 2);
      }
    }
  } else if (!(vrtype_and_fpos_storage & 4)) {
    vblock_index_byte_ct += kPglVblockSize / 2;
    if (variant_ct % 2) {
      // bugfix (22 Nov 2017): forgot to add offset in last block
      last_vrtype_byte_offset = 20 + (vblock_ct - 1) * (vblock_index_byte_ct + sizeof(int64_t)) + ((variant_ct % kPglVblockSize) / 2);
    }
    /*
  } else {
    vblock_index_byte_ct += kPglVblockSize;
    */
  }
  if (last_vrtype_byte_offset) {
    uint32_t last_vrtype_byte = 0;
#ifndef NO_MMAP
    if (ff == nullptr) {
      memcpy(&last_vrtype_byte, &(pgrp->fi.block_base[last_vrtype_byte_offset]), 1);
    } else {
#endif
      if (fseeko(ff, last_vrtype_byte_offset, SEEK_SET)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
        return kPglRetReadFail;
      }
      last_vrtype_byte = getc_unlocked(ff);
      if (last_vrtype_byte > 255) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
        return kPglRetReadFail;
      }
#ifndef NO_MMAP
    }
#endif
    if (last_vrtype_byte >> trailing_shift) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Nonzero trailing bits in last vrtype index byte.\n");
      return kPglRetMalformedInput;
    }
  }
  const uintptr_t* nonref_flags = pgrp->fi.nonref_flags;
  if (nonref_flags) {
    const uint32_t variant_ct_modl = variant_ct % kBitsPerWord;
    if (variant_ct % CHAR_BIT) {
      if (nonref_flags[variant_ct / kBitsPerWord] >> variant_ct_modl) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Nonzero trailing bits in last nonref_flags byte.\n");
        return kPglRetMalformedInput;
      }
    }
  }

  // could move most of this into plink2_common and make it multithreaded, if
  // speed is ever an issue.
  for (uint32_t vidx = 0; vidx < variant_ct; ++vidx) {
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: File read failure.\n");
      return kPglRetReadFail;
    }
    if (ValidateGeno(fread_end, vidx, pgrp, &fread_ptr, errstr_buf)) {
      return kPglRetMalformedInput;
    }
    const uint32_t vrtype = vrtypes[vidx];
    if (VrtypeHphase(vrtype)) {
      PgrDetectGenovecHets(pgrp->workspace_vec, sample_ct, pgrp->workspace_all_hets);
    }
    if (VrtypeMultiallelic(vrtype)) {
      // todo
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Validation of multiallelic data track is not implemented yet.\n");
      return kPglRetNotYetSupported;
    }
    // don't need pgrp->workspace_vec to store main genotypes past this point.
    if (VrtypeHphase(vrtype)) {
      if (ValidateHphase(fread_end, vidx, pgrp, &fread_ptr, errstr_buf)) {
        return kPglRetMalformedInput;
      }
    }
    if (vrtype & 0xe0) {
      if ((vrtype & 0xe0) == 0x80) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid record type for (0-based) variant #%u (phased dosage bit set, but main dosage bits unset).\n", vidx);
        return kPglRetMalformedInput;
      }
      PglErr reterr = ValidateDosage16(fread_end, vidx, pgrp, &fread_ptr, errstr_buf);
      if (reterr) {
        return reterr;
      }
    }
    if (fread_ptr != fread_end) {
      // possible todo: tolerate this at the end of a vblock.
      assert(fread_ptr < fread_end);
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Extra byte(s) in (0-based) variant record #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
  }
  return kPglRetSuccess;
}


BoolErr PgfiCleanup(PgenFileInfo* pgfip) {
  // memory is the responsibility of the caller
  if (pgfip->shared_ff) {
    if (fclose_null(&pgfip->shared_ff)) {
      return 1;
    }
#ifndef NO_MMAP
  } else if (pgfip->block_base != nullptr) {
    munmap(K_CAST(unsigned char*, pgfip->block_base), pgfip->file_size);
#endif
  }
  return 0;
}

BoolErr PgrCleanup(PgenReader* pgrp) {
  // assume file is open if pgr.ff is not null
  // memory is the responsibility of the caller for now
  if (!pgrp->ff) {
    return 0;
  }
  return fclose_null(&(pgrp->ff));
}


// ***** end PgenReader, begin {st,mt}_pgen_writer_t *****


void PreinitSpgw(STPgenWriter* spgwp) {
  spgwp->pgen_outfile = nullptr;
}

PglErr PwcInitPhase1(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, PgenWriterCommon* pwcp, FILE** pgen_outfile_ptr) {
  pwcp->allele_idx_offsets = allele_idx_offsets;
  pwcp->explicit_nonref_flags = nullptr;
  if (nonref_flags_storage == 3) {
    if (!explicit_nonref_flags) {
      return kPglRetImproperFunctionCall;
    }
    pwcp->explicit_nonref_flags = explicit_nonref_flags;
  }
  pwcp->variant_ct = variant_ct;
  pwcp->sample_ct = sample_ct;
  pwcp->phase_dosage_gflags = phase_dosage_gflags;
#ifndef NDEBUG
  pwcp->vblock_fpos = nullptr;
  pwcp->vrec_len_buf = nullptr;
  pwcp->vrtype_buf = nullptr;
  pwcp->fwrite_buf = nullptr;
  pwcp->fwrite_bufp = nullptr;
  pwcp->genovec_invert_buf = nullptr;
  pwcp->ldbase_genovec = nullptr;
  pwcp->ldbase_raregeno = nullptr;
  pwcp->ldbase_difflist_sample_ids = nullptr;
#endif
  pwcp->vidx = 0;

  FILE* pgen_outfile = fopen(fname, FOPEN_WB);
  *pgen_outfile_ptr = pgen_outfile;
  if (!pgen_outfile) {
    return kPglRetOpenFail;
  }
  fwrite_unlocked("l\x1b\x10", 3, 1, pgen_outfile);
  fwrite_unlocked(&(pwcp->variant_ct), sizeof(int32_t), 1, pgen_outfile);
  fwrite_unlocked(&(pwcp->sample_ct), sizeof(int32_t), 1, pgen_outfile);

  const unsigned char control_byte = (vrec_len_byte_ct - 1) + (4 * (phase_dosage_gflags != 0)) + (nonref_flags_storage << 6);
  pwcp->vrec_len_byte_ct = vrec_len_byte_ct;
  fwrite_unlocked(&control_byte, 1, 1, pgen_outfile);
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  uintptr_t header_bytes_left = vblock_ct * sizeof(int64_t) + variant_ct * vrec_len_byte_ct;
  if (phase_dosage_gflags) {
    // 8-bit vrtypes
    header_bytes_left += variant_ct;
  } else {
    // 4-bit vrtypes
    header_bytes_left += DivUp(variant_ct, 2);
  }
  if (nonref_flags_storage == 3) {
    header_bytes_left += DivUp(variant_ct, CHAR_BIT);
  }

  // this should be the position of the first variant
  pwcp->vblock_fpos_offset = 12 + header_bytes_left;

  uintptr_t zeroed_cachelines_needed = DivUp(header_bytes_left, kCacheline);
  if (zeroed_cachelines_needed > (kPglFwriteBlockSize / kCacheline)) {
    zeroed_cachelines_needed = kPglFwriteBlockSize / kCacheline;
  }
  // could wait until fwrite_buf is allocated, and make sure it's aligned?
  unsigned char zerobuf[kPglFwriteBlockSize];
  memset(zerobuf, 0, zeroed_cachelines_needed * kCacheline);
  while (header_bytes_left > kPglFwriteBlockSize) {
    fwrite_unlocked(zerobuf, kPglFwriteBlockSize, 1, pgen_outfile);
    header_bytes_left -= kPglFwriteBlockSize;
  }
  if (fwrite_checked(zerobuf, header_bytes_left, pgen_outfile)) {
    return kPglRetWriteFail;
  }
  return kPglRetSuccess;
}

uint32_t CountSpgwAllocCachelinesRequired(uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uint32_t max_vrec_len) {
  // vblock_fpos
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  uint32_t cachelines_required = Int64CtToCachelineCt(vblock_ct);

  // vrec_len_buf
  // overlapping uint32_t writes used, so (variant_ct * vrec_len_byte_ct) might
  // not be enough
  const uintptr_t vrec_len_byte_ct = BytesToRepresentUi(max_vrec_len);
  cachelines_required += DivUp((variant_ct - 1) * vrec_len_byte_ct + sizeof(int32_t), kCacheline);

  // vrtype_buf
  if (phase_dosage_gflags) {
    cachelines_required += DivUp(variant_ct, kCacheline);
  } else {
    cachelines_required += DivUp(variant_ct, kCacheline * 2);
  }

  // genovec_invert_buf, ldbase_genovec
  cachelines_required += 2 * QuaterCtToCachelineCt(sample_ct);

  const uint32_t max_difflist_len = 2 * (sample_ct / kPglMaxDifflistLenDivisor);
  // ldbase_raregeno
  cachelines_required += QuaterCtToCachelineCt(max_difflist_len);

  // ldbase_difflist_sample_ids
  cachelines_required += 1 + (max_difflist_len / kInt32PerCacheline);

  // fwrite_buf
  cachelines_required += DivUp(max_vrec_len + kPglFwriteBlockSize - 1, kCacheline);
  if (phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) {
    // phasepresent, phaseinfo
    cachelines_required += 2 * BitCtToCachelineCt(sample_ct);
  }
  // possible todo: dosage (doesn't currently need an allocation, but that's
  // unlikely to remain true--e.g. get_ref_nonref_genotype_counts_and_dosages
  // tends to use workspace_vec when a function it calls doesn't use it...)
  return cachelines_required;
}

static_assert(kPglMaxAltAlleleCt == 254, "Need to update SpgwInitPhase1().");
PglErr SpgwInitPhase1(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, STPgenWriter* spgwp, uintptr_t* alloc_cacheline_ct_ptr, uint32_t* max_vrec_len_ptr) {
  assert(variant_ct);
  assert(sample_ct);

  // separate from MpgwInitPhase1's version of this computation since the
  // latter wants a better bound on the compressed size of an entire vblock
  // than max_vrec_len * kPglVblockSize...
  uint64_t max_vrec_len = QuaterCtToByteCt(sample_ct);
  uintptr_t max_alt_ct_p1 = 2;
  if (allele_idx_offsets && (allele_idx_offsets[variant_ct] != 2 * variant_ct)) {
    assert(allele_idx_offsets[0] == 0);
    assert(allele_idx_offsets[variant_ct] > 2 * variant_ct);
    // could add this as a parameter, since caller should know...
    max_alt_ct_p1 = 3;
    uintptr_t prev_offset = 0;
    for (uint32_t vidx = 1; vidx <= variant_ct; ++vidx) {
      const uintptr_t cur_offset = allele_idx_offsets[vidx];
      if (cur_offset - prev_offset > max_alt_ct_p1) {
        max_alt_ct_p1 = cur_offset - prev_offset;
      }
      prev_offset = cur_offset;
    }
    // nonmissingness array
    max_vrec_len += DivUp(sample_ct, CHAR_BIT) + GetAux1AlleleByteCt(max_alt_ct_p1 - 1, sample_ct);
    // try to permit uncompressed records to be larger than this, only error
    // out when trying to write a larger compressed record?  (might not be
    // worth it.)
  }
  if (phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) {
    max_vrec_len += 2 * DivUp(sample_ct, CHAR_BIT);
  }
  if (phase_dosage_gflags & kfPgenGlobalDosagePresent) {
    const uint32_t dosage_phase_gflag = (phase_dosage_gflags / kfPgenGlobalDosagePhasePresent) & 1;
    // aux3, aux4
    max_vrec_len += (1 + dosage_phase_gflag) * DivUp(sample_ct, 8);
    // aux5
    max_vrec_len += (2 + 2 * dosage_phase_gflag) * S_CAST(uint64_t, sample_ct) * (max_alt_ct_p1 - 1);

  }
  if (max_vrec_len >= kPglMaxBytesPerVariant) {
#ifdef __LP64__
    max_vrec_len = kPglMaxBytesPerVariant;
#else
    return kPglRetNomem;
#endif
  }
  *max_vrec_len_ptr = max_vrec_len;
  const uintptr_t vrec_len_byte_ct = BytesToRepresentUi(max_vrec_len);

  PglErr reterr = PwcInitPhase1(fname, allele_idx_offsets, explicit_nonref_flags, variant_ct, sample_ct, phase_dosage_gflags, nonref_flags_storage, vrec_len_byte_ct, &(spgwp->pwc), &(spgwp->pgen_outfile));
  if (!reterr) {
    *alloc_cacheline_ct_ptr = CountSpgwAllocCachelinesRequired(variant_ct, sample_ct, phase_dosage_gflags, max_vrec_len);
  }
  return reterr;
}

static_assert(kPglMaxAltAlleleCt == 254, "Need to update MpgwInitPhase1().");
void MpgwInitPhase1(const uintptr_t* __restrict allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uintptr_t* alloc_base_cacheline_ct_ptr, uint64_t* alloc_per_thread_cacheline_ct_ptr, uint32_t* vrec_len_byte_ct_ptr, uint64_t* vblock_cacheline_ct_ptr) {
  assert(variant_ct);
  assert(sample_ct);
  // vblock_fpos
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  uint32_t alloc_base_cacheline_ct = Int64CtToCachelineCt(vblock_ct);

  // vrtype_buf
  if (phase_dosage_gflags) {
    alloc_base_cacheline_ct += DivUp(variant_ct, kCacheline);
  } else {
    alloc_base_cacheline_ct += DivUp(variant_ct, kCacheline * 2);
  }

  // pwcs
  uint64_t alloc_per_thread_cacheline_ct = DivUp(sizeof(PgenWriterCommon), kCacheline);

  // genovec_invert_buf, ldbase_genovec
  alloc_per_thread_cacheline_ct += 2 * QuaterCtToCachelineCt(sample_ct);

  const uint32_t max_difflist_len = 2 * (sample_ct / kPglMaxDifflistLenDivisor);
  // ldbase_raregeno
  alloc_per_thread_cacheline_ct += QuaterCtToCachelineCt(max_difflist_len);

  // ldbase_difflist_sample_ids
  alloc_per_thread_cacheline_ct += 1 + (max_difflist_len / kInt32PerCacheline);

  uint64_t max_vrec_len = QuaterCtToByteCt(sample_ct);
  if (phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) {
    max_vrec_len += 2 * DivUp(sample_ct, CHAR_BIT);
  }
  const uint32_t dosage_gflag = (phase_dosage_gflags / kfPgenGlobalDosagePresent) & 1;
  const uint32_t dosage_phase_gflag = (phase_dosage_gflags / kfPgenGlobalDosagePhasePresent) & 1;
  if (dosage_gflag) {
    max_vrec_len += ((1 + dosage_phase_gflag) * DivUp(sample_ct, CHAR_BIT)) + (2 + 2 * dosage_phase_gflag) * S_CAST(uint64_t, sample_ct);
  }
  const uint32_t max_vblock_size = MINV(variant_ct, kPglVblockSize);
  // max_vrec_len is a uint64_t
  uint64_t max_vblock_byte_ct = max_vrec_len * max_vblock_size;
  if (max_vrec_len >= kPglMaxBytesPerVariant) {
    max_vrec_len = kPglMaxBytesPerVariant;
    max_vblock_byte_ct = kPglMaxBytesPerVariant * S_CAST(uint64_t, max_vblock_size);
  } else if (allele_idx_offsets && (allele_idx_offsets[variant_ct] != 2 * variant_ct)) {
    assert(allele_idx_offsets[0] == 0);
    assert(allele_idx_offsets[variant_ct] > 2 * variant_ct);
    // When multiallelic variants are present, larger write buffers are needed.
    // we compute the largest possible size here.
    //
    // For aux1, a nonmissingness array with (sample_ct + 7) / 8 bytes is
    // always needed.  on top of that,
    //   alt ct  additional bytes required
    //   ------  -------------------------
    //        2        (sample_ct + 3) / 4
    //        3        (sample_ct + 1) / 2
    //     4-15                  sample_ct
    //   16-255              2 * sample_ct
    //
    // For aux5, (2 + 2 * dosage_phase_gflag) additional bytes are needed per
    // sample x additional alt allele (yes, it isn't hard to exceed the ~4GB
    // variant record size limit here).
    //
    // Between the two, we have a piecewise linear function with up to 5
    // segments (last segment could correspond to the record size limit).
    // Okay, the last segment means "out of memory" unless we have something
    // like 256TB RAM, but still.
    uintptr_t prev_offset = 0;
    uint32_t vidx = 0;
    const uint32_t extra_bytes_base = DivUp(sample_ct, CHAR_BIT);
    const uint64_t extra_bytes_max = kPglMaxBytesPerVariant - max_vrec_len;
    const uint64_t extra_dosage_bytes_per_alt = dosage_phase_gflag * (2 + 2 * dosage_phase_gflag) * S_CAST(uint64_t, sample_ct);
    uint64_t extra_byte_cts[4];
    uint32_t extra_alt_ceil = kPglMaxAltAlleleCt + 1;

    // alt_ct == 2
    uint64_t cur_extra_byte_ct = extra_bytes_base + DivUp(sample_ct, 4) + extra_dosage_bytes_per_alt;
    extra_byte_cts[0] = cur_extra_byte_ct;
    extra_byte_cts[1] = 0;  // force initialization
    extra_byte_cts[2] = 0;
    extra_byte_cts[3] = 0;
    if (cur_extra_byte_ct >= extra_bytes_max) {
      extra_alt_ceil = 2;
    } else {
      // alt_ct == 3
      cur_extra_byte_ct = extra_bytes_base + DivUp(sample_ct, 4) + 2 * extra_dosage_bytes_per_alt;
      extra_byte_cts[1] = cur_extra_byte_ct;
      if (cur_extra_byte_ct >= extra_bytes_max) {
        extra_alt_ceil = 3;
      } else {
        // alt_ct in [4, 15]
        cur_extra_byte_ct = extra_bytes_base + sample_ct + 3 * extra_dosage_bytes_per_alt;
        extra_byte_cts[2] = cur_extra_byte_ct;
        if (cur_extra_byte_ct >= extra_bytes_max) {
          extra_alt_ceil = 4;
        } else if (cur_extra_byte_ct + 11 * extra_dosage_bytes_per_alt >= extra_bytes_max) {
          extra_alt_ceil = 5 + (extra_bytes_max - cur_extra_byte_ct - 1) / extra_dosage_bytes_per_alt;
        } else {
          // alt_ct in [16, 254]
          cur_extra_byte_ct = extra_bytes_base + 2 * sample_ct + 15 * extra_dosage_bytes_per_alt;
          extra_byte_cts[3] = cur_extra_byte_ct;
          if (cur_extra_byte_ct >= extra_bytes_max) {
            extra_alt_ceil = 16;
          } else if (cur_extra_byte_ct + 238 * extra_dosage_bytes_per_alt >= extra_bytes_max) {
            extra_alt_ceil = 17 + (extra_bytes_max - cur_extra_byte_ct - 1) / extra_dosage_bytes_per_alt;
          }
        }
      }
    }
    uint64_t extra_nonceil_altp1_total = 0;
    uint32_t extra_alt_ceil_ct = 0;
    const uint64_t uncompressed_biallelic_vrec_len = max_vrec_len;
    uint32_t altx_seen_mask = 0;
    uint32_t max_alt_ct_p1 = 3;
    while (1) {
      uint32_t vblock_end = vidx + kPglVblockSize;
      if (vblock_end > variant_ct) {
        if (vidx == variant_ct) {
          break;
        }
        vblock_end = variant_ct;
      }
      uint32_t altx_seen[4];
      ZeroUiArr(4, altx_seen);
      for (; vidx < vblock_end;) {
        const uintptr_t cur_offset = allele_idx_offsets[++vidx];
        const uint32_t alt_ct_p1 = cur_offset - prev_offset;
        if (alt_ct_p1 > 2) {
          if (alt_ct_p1 >= extra_alt_ceil) {
            ++extra_alt_ceil_ct;
          } else {
            // don't need to track this when we hit the ceiling
            if (alt_ct_p1 > max_alt_ct_p1) {
              max_alt_ct_p1 = alt_ct_p1;
            }

            extra_nonceil_altp1_total += alt_ct_p1;
            if (alt_ct_p1 < 5) {
              altx_seen[alt_ct_p1 - 3] += 1;
            } else {
              altx_seen[2 + (alt_ct_p1 >= 16)] += 1;
            }
          }
        }
        prev_offset = cur_offset;
      }
      uint64_t cur_vblock_byte_ct = uncompressed_biallelic_vrec_len * (vblock_end - vidx);
      cur_vblock_byte_ct += extra_alt_ceil_ct * extra_bytes_max;
      for (uint32_t uii = 0; uii < 4; ++uii) {
        if (altx_seen[uii]) {
          const uint32_t cur_seen_ct = altx_seen[uii];
          altx_seen_mask |= 1 << uii;
          cur_vblock_byte_ct += cur_seen_ct * extra_byte_cts[uii];
        }
      }
      if (dosage_gflag) {
        cur_vblock_byte_ct += (extra_nonceil_altp1_total - altx_seen[0] * 3 - altx_seen[1] * 4 - altx_seen[2] * 5 - altx_seen[3] * 17) * extra_dosage_bytes_per_alt;
      }
      if (cur_vblock_byte_ct > max_vblock_byte_ct) {
        max_vblock_byte_ct = cur_vblock_byte_ct;
      }
    }
    if (extra_alt_ceil_ct) {
      max_vrec_len = kPglMaxBytesPerVariant;
    } else {
      max_vrec_len = uncompressed_biallelic_vrec_len + extra_byte_cts[31 - __builtin_clz(altx_seen_mask)];
      if (dosage_gflag && (max_alt_ct_p1 >= 6)) {
        if (max_alt_ct_p1 >= 17) {
          max_vrec_len += (max_alt_ct_p1 - 17) * extra_dosage_bytes_per_alt;
        } else {
          max_vrec_len += (max_alt_ct_p1 - 5) * extra_dosage_bytes_per_alt;
        }
      }
    }
  }
  // vrec_len_buf
  // previously used overlapping uint32_t writes-to-memory, but that was
  // incompatible with multithreaded compression
  *vrec_len_byte_ct_ptr = BytesToRepresentUi(max_vrec_len);
  *alloc_base_cacheline_ct_ptr = alloc_base_cacheline_ct + DivUp(S_CAST(uintptr_t, variant_ct) * (*vrec_len_byte_ct_ptr), kCacheline);

  // main write buffer
  *vblock_cacheline_ct_ptr = DivUpU64(max_vblock_byte_ct, kCacheline);
  *alloc_per_thread_cacheline_ct_ptr = alloc_per_thread_cacheline_ct + (*vblock_cacheline_ct_ptr);
}


void PwcInitPhase2(uintptr_t fwrite_cacheline_ct, uint32_t thread_ct, PgenWriterCommon** pwcs, unsigned char* pwc_alloc) {
  const uint32_t variant_ct = pwcs[0]->variant_ct;
  unsigned char* alloc_iter = pwc_alloc;
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  const PgenGlobalFlags phase_dosage_gflags = pwcs[0]->phase_dosage_gflags;
  uint32_t vrtype_buf_bytes;
  if (phase_dosage_gflags) {
    vrtype_buf_bytes = RoundUpPow2(variant_ct, kCacheline);
  } else {
    vrtype_buf_bytes = DivUp(variant_ct, kCacheline * 2) * kCacheline;
  }
  pwcs[0]->vblock_fpos = R_CAST(uint64_t*, alloc_iter);
  alloc_iter = &(alloc_iter[Int64CtToCachelineCt(vblock_ct) * kCacheline]);

  pwcs[0]->vrec_len_buf = alloc_iter;
  alloc_iter = &(alloc_iter[RoundUpPow2(variant_ct * pwcs[0]->vrec_len_byte_ct, kCacheline)]);

  pwcs[0]->vrtype_buf = R_CAST(uintptr_t*, alloc_iter);
  // spgw_append() assumes these bytes are zeroed out
  memset(pwcs[0]->vrtype_buf, 0, vrtype_buf_bytes);
  alloc_iter = &(alloc_iter[vrtype_buf_bytes]);

  const uint32_t sample_ct = pwcs[0]->sample_ct;
  const uint32_t genovec_byte_alloc = QuaterCtToCachelineCt(sample_ct) * kCacheline;
  const uint32_t max_difflist_len = 2 * (sample_ct / kPglMaxDifflistLenDivisor);
  for (uint32_t tidx = 0; tidx < thread_ct; ++tidx) {
    if (tidx) {
      pwcs[tidx]->vblock_fpos = pwcs[0]->vblock_fpos;
      pwcs[tidx]->vrec_len_buf = pwcs[0]->vrec_len_buf;
      pwcs[tidx]->vrtype_buf = pwcs[0]->vrtype_buf;
    }
    pwcs[tidx]->genovec_invert_buf = R_CAST(uintptr_t*, alloc_iter);
    alloc_iter = &(alloc_iter[genovec_byte_alloc]);
    pwcs[tidx]->ldbase_genovec = R_CAST(uintptr_t*, alloc_iter);
    alloc_iter = &(alloc_iter[genovec_byte_alloc]);

    pwcs[tidx]->ldbase_raregeno = R_CAST(uintptr_t*, alloc_iter);
    alloc_iter = &(alloc_iter[QuaterCtToCachelineCt(max_difflist_len) * kCacheline]);
    pwcs[tidx]->ldbase_difflist_sample_ids = R_CAST(uint32_t*, alloc_iter);
    alloc_iter = &(alloc_iter[(1 + (max_difflist_len / kInt32PerCacheline)) * kCacheline]);

    pwcs[tidx]->fwrite_buf = alloc_iter;
    pwcs[tidx]->fwrite_bufp = alloc_iter;
    alloc_iter = &(alloc_iter[fwrite_cacheline_ct * kCacheline]);
  }
}

void SpgwInitPhase2(uint32_t max_vrec_len, STPgenWriter* spgwp, unsigned char* spgw_alloc) {
  uintptr_t fwrite_cacheline_ct = DivUp(max_vrec_len + kPglFwriteBlockSize - 1, kCacheline);
  PgenWriterCommon* pwcp = &(spgwp->pwc);
  if (pwcp->phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) {
    // er, why do we need this extra space?  Document this.
    fwrite_cacheline_ct += 2 * BitCtToCachelineCt(pwcp->sample_ct);
  }
  PwcInitPhase2(fwrite_cacheline_ct, 1, &pwcp, spgw_alloc);
}

PglErr MpgwInitPhase2(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, uintptr_t vblock_cacheline_ct, uint32_t thread_ct, unsigned char* mpgw_alloc, MTPgenWriter* mpgwp) {
  assert(thread_ct);
  const uintptr_t pwc_byte_ct = RoundUpPow2(sizeof(PgenWriterCommon), kCacheline);
  for (uint32_t tidx = 0; tidx < thread_ct; ++tidx) {
    mpgwp->pwcs[tidx] = R_CAST(PgenWriterCommon*, &(mpgw_alloc[tidx * pwc_byte_ct]));
  }
  PglErr reterr = PwcInitPhase1(fname, allele_idx_offsets, explicit_nonref_flags, variant_ct, sample_ct, phase_dosage_gflags, nonref_flags_storage, vrec_len_byte_ct, mpgwp->pwcs[0], &(mpgwp->pgen_outfile));
  if (!reterr) {
    mpgwp->thread_ct = thread_ct;
    if (thread_ct > 1) {
      for (uint32_t tidx = 1; tidx < thread_ct; ++tidx) {
        memcpy(mpgwp->pwcs[tidx], mpgwp->pwcs[0], sizeof(PgenWriterCommon));
        mpgwp->pwcs[tidx]->vidx = tidx * kPglVblockSize;
      }
    }
    PwcInitPhase2(vblock_cacheline_ct, thread_ct, mpgwp->pwcs, &(mpgw_alloc[thread_ct * pwc_byte_ct]));
  }
  return reterr;
}


void CountLdAndInvertedLdDiffs(const uintptr_t* __restrict ldbase_genovec, const uintptr_t* __restrict genovec, uint32_t sample_ct, uint32_t* ld_diff_ctp, uint32_t* ld_inv_diff_ctp) {
  // Requires trailing bits to be zeroed out.
  const uint32_t word_ct = QuaterCtToWordCt(sample_ct);
  const uintptr_t* genovec_end = &(genovec[word_ct]);
  uint32_t ld_diff_ct = 0;
  uint32_t ld_inv_diff_ct = 0;
  // construct the words we want to popcount_quatervec_01 on the fly
  const VecUL m1 = VCONST_UL(kMask5555);
  const VecUL m2 = VCONST_UL(kMask3333);
  const VecUL m4 = VCONST_UL(kMask0F0F);
  const VecUL* ldbase_vvec_iter = R_CAST(const VecUL*, ldbase_genovec);
  const VecUL* geno_vvec_iter = R_CAST(const VecUL*, genovec);
  uint32_t full_vecs_left = 3 * (word_ct / (3 * kWordsPerVec));
  UniVec acc_ld;
  UniVec acc_ld_inv;
  while (1) {
    acc_ld.vi = vecul_setzero();
    acc_ld_inv.vi = vecul_setzero();
    const VecUL* geno_vvec_stop;
    if (full_vecs_left < 60) {
      if (!full_vecs_left) {
        break;
      }
      geno_vvec_stop = &(geno_vvec_iter[full_vecs_left]);
      full_vecs_left = 0;
    } else {
      geno_vvec_stop = &(geno_vvec_iter[60]);
      full_vecs_left -= 60;
    }
    do {
      VecUL loader_ldbase1 = *ldbase_vvec_iter++;
      VecUL loader_geno1 = *geno_vvec_iter++;
      VecUL loader_ldbase2 = *ldbase_vvec_iter++;
      VecUL loader_geno2 = *geno_vvec_iter++;
      VecUL xor1 = loader_ldbase1 ^ loader_geno1;
      VecUL xor2 = loader_ldbase2 ^ loader_geno2;
      VecUL xor_shifted1 = vecul_srli(xor1, 1);
      VecUL xor_shifted2 = vecul_srli(xor2, 1);
      // xor(_low)  xor_shifted  loader_geno   result
      //         1                                  1
      //         0            0            0        1
      //         0            0            1        0
      //         0            1            0        0
      //         0            1            1        1
      // gah, don't see a way to avoid throwing in an extra xor for
      // loader_geno...
      VecUL count_ld_inv = (xor1 | (xor_shifted1 ^ loader_geno1 ^ m1)) & m1;
      loader_ldbase1 = *ldbase_vvec_iter++;
      VecUL count_ld = (xor1 | xor_shifted1) & m1;
      loader_geno1 = *geno_vvec_iter++;
      count_ld_inv = count_ld_inv + ((xor2 | (xor_shifted2 ^ loader_geno2 ^ m1)) & m1);
      xor1 = loader_ldbase1 ^ loader_geno1;
      count_ld = count_ld + ((xor2 | xor_shifted2) & m1);
      xor_shifted1 = vecul_srli(xor1, 1);
      count_ld_inv = count_ld_inv + ((xor1 | (xor_shifted1 ^ loader_geno1 ^ m1)) & m1);
      count_ld = count_ld + ((xor1 | xor_shifted1) & m1);
      // now count_ld and count_ld_inv each have 64 2-bit values from 0-3

      count_ld_inv = (count_ld_inv & m2) + (vecul_srli(count_ld_inv, 2) & m2);
      count_ld = (count_ld & m2) + (vecul_srli(count_ld, 2) & m2);
      // now they have 32 4-bit values from 0-6

      acc_ld_inv.vi = acc_ld_inv.vi + ((count_ld_inv + vecul_srli(count_ld_inv, 4)) & m4);
      acc_ld.vi = acc_ld.vi + ((count_ld + vecul_srli(count_ld, 4)) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecUL m8 = VCONST_UL(kMask00FF);
    acc_ld_inv.vi = (acc_ld_inv.vi & m8) + (vecul_srli(acc_ld_inv.vi, 8) & m8);
    acc_ld.vi = (acc_ld.vi & m8) + (vecul_srli(acc_ld.vi, 8) & m8);
    ld_inv_diff_ct += UniVecHsum16(acc_ld_inv);
    ld_diff_ct += UniVecHsum16(acc_ld);
  }
  const uintptr_t* ldbase_iter = R_CAST(const uintptr_t*, ldbase_vvec_iter);
  const uintptr_t* genovec_iter = R_CAST(const uintptr_t*, geno_vvec_iter);
  while (genovec_iter < genovec_end) {
    uintptr_t ldbase_word = *ldbase_iter++;
    uintptr_t geno_word = *genovec_iter++;
    uintptr_t xor_result = ldbase_word ^ geno_word;
    uintptr_t xor_result_shifted = xor_result >> 1;
    ld_diff_ct += Popcount01Word((xor_result | xor_result_shifted) & kMask5555);
    ld_inv_diff_ct += Popcount01Word((xor_result | (xor_result_shifted ^ (~geno_word))) & kMask5555);
  }
  *ld_diff_ctp = ld_diff_ct;
  // trailing entries in last word are always "different"
  *ld_inv_diff_ctp = ld_inv_diff_ct - ((-sample_ct) & (kBitsPerWordD2 - 1));
}

uint32_t CountLdAndInvertedLdDiffsList(const uintptr_t* __restrict ldbase_raregeno, const uint32_t* __restrict ldbase_difflist_sample_ids, const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t ldbase_difflist_len, uint32_t difflist_len, uint32_t* ld_diff_ctp, uint32_t* ld_inv_diff_ctp) {
  // assumes ldbase_difflist_sample_ids[ldbase_difflist_len] == sample_ct
  // assumes variant isn't multiallelic

  // only the count(s) with aligned common_geno values are valid.  e.g. if
  // ldbase_common_geno and difflist_common_geno are both zero, the ld_inv_diff
  // return value can be anything, while if they're both three, ld_diff and
  // ld_inv_diff are both accurate.

  // some similarities to ParseLdAndMergeDifflistSubset(), but much
  // simpler.
  // noticeably slower than CountLdAndInvertedLdDiffs() when the lists
  // aren't tiny.
  // possible todo: take threshold into account?

  uint32_t collision_ct = 0;
  uint32_t ld_diff_ct = 0;
  uint32_t ld_inv_diff_ct = 0;
  uint32_t ldbase_sample_idx = ldbase_difflist_sample_ids[0];
  uint32_t ldbase_difflist_idx = 1;
  // this loop is a bit slow.  attempt to bail halfway through?
  for (uint32_t difflist_idx = 0; difflist_idx < difflist_len; ++difflist_idx) {
    const uint32_t raw_sample_idx = difflist_sample_ids[difflist_idx];
    while (ldbase_sample_idx < raw_sample_idx) {
      ldbase_sample_idx = ldbase_difflist_sample_ids[ldbase_difflist_idx++];
    }
    if (ldbase_sample_idx > raw_sample_idx) {
      continue;
    }
    const uint32_t cur_raregeno = GetQuaterarrEntry(raregeno, difflist_idx);
    const uint32_t cur_ldbase_raregeno = GetQuaterarrEntry(ldbase_raregeno, ldbase_difflist_idx - 1);
    const uint32_t cur_inv_raregeno = (6 - cur_raregeno) & 3;
    ld_diff_ct += (cur_ldbase_raregeno != cur_raregeno);
    ldbase_sample_idx = ldbase_difflist_sample_ids[ldbase_difflist_idx++];
    ++collision_ct;
    ld_inv_diff_ct += (cur_ldbase_raregeno != cur_inv_raregeno);
  }
  // no more collisions, don't actually need to look at rest of
  // ldbase_difflist
  const uint32_t base_diff_ct = ldbase_difflist_len + difflist_len - 2 * collision_ct;
  *ld_diff_ctp = base_diff_ct + ld_diff_ct;
  *ld_inv_diff_ctp = base_diff_ct + ld_inv_diff_ct;
  return 1;
}

uint32_t SaveLdDifflist(const uintptr_t* __restrict genovec, const uintptr_t* __restrict ldbase_genovec, uintptr_t common_geno, uint32_t difflist_len, PgenWriterCommon* pwcp) {
  unsigned char* fwrite_bufp = pwcp->fwrite_bufp;
  if (!difflist_len) {
    *fwrite_bufp = 0;
    pwcp->fwrite_bufp = &(fwrite_bufp[1]);
    return 1;
  }
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(difflist_len, fwrite_bufp);
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(pwcp->sample_ct);
  const uintptr_t common_geno_word = common_geno * kMask5555;
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
#ifdef __arm__
#  error "Unaligned accesses in SaveLdDifflist()."
#endif
  uintptr_t* raregeno_iter = R_CAST(uintptr_t*, &(extra_byte_cts_iter[group_ct - 1]));
  fwrite_bufp = &(extra_byte_cts_iter[group_ct + (difflist_len - 1) / 4]);
  unsigned char* last_group_vint_start = fwrite_bufp;
  uintptr_t raregeno_word = 0;
  uint32_t last_sample_idx = 0;
  uint32_t difflist_idx = 0;
  uint32_t widx = 0;
  while (1) {
    const uintptr_t cur_geno_word = genovec[widx];
    uintptr_t xor_word = ldbase_genovec? ldbase_genovec[widx] : common_geno_word;
    xor_word ^= cur_geno_word;
    if (xor_word) {
      const uint32_t sample_idx_base = widx * kBitsPerWordD2;
      do {
        const uint32_t sample_idx_lowbits = CTZLU(xor_word) / 2;
        const uint32_t new_sample_idx = sample_idx_base + sample_idx_lowbits;
        raregeno_word |= ((cur_geno_word >> (2 * sample_idx_lowbits)) & 3) << (2 * (difflist_idx % kBitsPerWordD2));
        if (!(difflist_idx % kPglDifflistGroupSize)) {
          group_first_sample_ids_iter = R_CAST(unsigned char*, memcpya(group_first_sample_ids_iter, &new_sample_idx, sample_id_byte_ct));
          if (difflist_idx) {
            *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
          }
          last_group_vint_start = fwrite_bufp;
        } else {
          assert(new_sample_idx >= last_sample_idx + 1);
          fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
        }
        ++difflist_idx;
        last_sample_idx = new_sample_idx;
        if (difflist_idx == difflist_len) {
          SubwordStore(raregeno_word, 1 + (((difflist_len - 1) / 4) % sizeof(intptr_t)), raregeno_iter);
          pwcp->fwrite_bufp = fwrite_bufp;
          return fwrite_bufp - fwrite_bufp_start;
        }
        if (!(difflist_idx % kBitsPerWordD2)) {
          *raregeno_iter++ = raregeno_word;
          raregeno_word = 0;
        }
        xor_word &= (~(3 * k1LU)) << (2 * sample_idx_lowbits);
      } while (xor_word);
    }
    ++widx;
  }
}

void OnebitPreprocessBuf(const uintptr_t* __restrict genovec, uint32_t sample_ct, uint32_t common2_code, uintptr_t* __restrict genovec_buf) {
  assert(sample_ct);
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  // todo: look for better ways to perform some of these operations
  const VecUL* geno_vvec_iter = R_CAST(const VecUL*, genovec);
  const VecUL* geno_vvec_end = &(geno_vvec_iter[vec_ct]);
  VecUL* write_iter = R_CAST(VecUL*, genovec_buf);
  const VecUL m1 = VCONST_UL(kMask5555);
  if (common2_code < 5) {
    if (common2_code == 1) {
      // 11 -> 10, everything else unchanged
      // todo: check if these loops are actually faster as simple while loops
      // todo: check if it's better to unroll these loops to process 2 __m128is
      //       at a time
      do {
        const VecUL cur_geno = *geno_vvec_iter++;
        *write_iter++ = (~(m1 & vecul_srli(cur_geno, 1))) & cur_geno;
      } while (geno_vvec_iter < geno_vvec_end);
    } else if (common2_code == 3) {
      // 00 -> 00, 01 -> 10, 10 -> 10, 11 -> 01
      do {
        const VecUL cur_geno = *geno_vvec_iter++;
        const VecUL cur_geno_rshift = vecul_srli(cur_geno, 1);
        const VecUL cur_geno_xor_masked = (cur_geno ^ cur_geno_rshift) & m1;
        const VecUL cur_geno_or_masked = (cur_geno | cur_geno_rshift) & m1;
        *write_iter++ = cur_geno_xor_masked + cur_geno_or_masked;
      } while (geno_vvec_iter < geno_vvec_end);
    } else {
      assert(common2_code == 2);
      // 00 -> 00, 01 -> 10, 10 -> 01, 11 -> 10
      do {
        const VecUL cur_geno = *geno_vvec_iter++;
        const VecUL cur_geno_or_masked = (cur_geno | vecul_srli(cur_geno, 1)) & m1;
        const VecUL cur_geno_lowbits = cur_geno & m1;
        *write_iter++ = cur_geno_lowbits + cur_geno_or_masked;
      } while (geno_vvec_iter < geno_vvec_end);
    }
  } else {
    if (common2_code == 5) {
      // 00 -> 10, 01 -> 00, 10 -> 01, 11 -> 10
      do {
        const VecUL cur_geno = *geno_vvec_iter++;
        const VecUL cur_geno_rshift = vecul_srli(cur_geno, 1);
        const VecUL cur_geno_not_xor_masked = (~(cur_geno ^ cur_geno_rshift)) & m1;
        const VecUL cur_geno_rshift_masked = cur_geno_rshift & m1;
        *write_iter++ = cur_geno_not_xor_masked + (cur_geno_not_xor_masked | cur_geno_rshift_masked);
      } while (geno_vvec_iter < geno_vvec_end);
    } else if (common2_code == 9) {
      // 00 -> 10, 01 -> 10, 10 -> 00, 11 -> 01
      const VecUL not_m1 = VCONST_UL(kMaskAAAA);
      do {
        const VecUL cur_geno = *geno_vvec_iter++;
        *write_iter++ = (cur_geno ^ not_m1) - ((~not_m1) & ((~vecul_srli(cur_geno, 1)) & cur_geno));
      } while (geno_vvec_iter < geno_vvec_end);
    } else {
      assert(common2_code == 6);
      // 00 -> 10, 01 -> 00, 10 -> 10, 11 -> 01
      do {
        const VecUL cur_geno = *geno_vvec_iter++;
        const VecUL cur_geno_not_lowbits = (~cur_geno) & m1;
        const VecUL cur_geno_rshift_masked = vecul_srli(cur_geno, 1) & m1;
        *write_iter++ = cur_geno_not_lowbits + (cur_geno_not_lowbits | cur_geno_rshift_masked);
      } while (geno_vvec_iter < geno_vvec_end);
    }
  }
}

uint32_t SaveOnebit(const uintptr_t* __restrict genovec, uint32_t common2_code, uint32_t onebit_difflist_len, PgenWriterCommon* pwcp) {
  // Uses ldbase_genovec as a temporary buffer.

  // common2_code is expected to have the difference between the common
  // genotype values in bits 0-1, and the smaller common genotype value in bits
  // 2-3.
  unsigned char* fwrite_bufp_start = pwcp->fwrite_bufp;
  *fwrite_bufp_start = common2_code;
  const uint32_t sample_ct = pwcp->sample_ct;
  uintptr_t* __restrict genovec_buf = pwcp->ldbase_genovec;
  // There's a 4-byte-interleaved format which is slightly more efficient for
  // unsubsetted handling (~10 fewer compression/decompression operations per
  // 32 genotypes), but that's only a 1-2% speedup, which probably isn't worth
  // the more annoying subsetting.
  //
  // Any 10s and 11s are saved as 00 in this part.
  // Similar logic is used to handle the other five possibilities (00/10,
  // 00/11, 01/10, 01/11, 10/11); all of them should be expected to actually
  // happen.  (E.g. 01/11 can happen at a high MAF variant when there's lots of
  // missing data.)  To reduce branching, we preprocess genovec_buf to have
  // even bit set iff the corresponding genotype is equal to the high common
  // genotype value, and odd bit set iff the corresponding genotype is one of
  // the two uncommon values.  (There may be a better way to do this, analogous
  // to the simpler decompression algorithm.)
  OnebitPreprocessBuf(genovec, sample_ct, common2_code, genovec_buf);
  ZeroTrailingQuaters(sample_ct, genovec_buf);
  const uint32_t word_read_ct = QuaterCtToWordCt(sample_ct);
#ifdef __arm__
#  error "Unaligned accesses in SaveOnebit()."
#endif
  Halfword* fwrite_bufp_alias_halfword = R_CAST(Halfword*, &(fwrite_bufp_start[1]));
  for (uint32_t widx = 0; widx < word_read_ct; ++widx) {
    const uintptr_t cur_buf_word = genovec_buf[widx] & kMask5555;
    fwrite_bufp_alias_halfword[widx] = PackWordToHalfword(cur_buf_word);
  }
  const uint32_t onebit_block_len = (sample_ct + 15) / CHAR_BIT;
  unsigned char* fwrite_bufp = Vint32Append(onebit_difflist_len, &(fwrite_bufp_start[onebit_block_len]));
  // the rest is almost identical to SaveLdDifflist()
  if (!onebit_difflist_len) {
    pwcp->fwrite_bufp = fwrite_bufp;
    return (onebit_block_len + 1);
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(onebit_difflist_len, kPglDifflistGroupSize);
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
  uintptr_t* raregeno_iter = R_CAST(uintptr_t*, &(extra_byte_cts_iter[group_ct - 1]));
  fwrite_bufp = &(extra_byte_cts_iter[group_ct + (onebit_difflist_len - 1) / 4]);
  unsigned char* last_group_vint_start = fwrite_bufp;
  uintptr_t raregeno_word = 0;
  uint32_t last_sample_idx = 0;
  uint32_t difflist_idx = 0;
  uint32_t widx = 0;
  while (1) {
    uintptr_t xor_word = genovec_buf[widx] & kMaskAAAA;
    if (xor_word) {
      const uintptr_t cur_geno_word = genovec[widx];
      const uint32_t sample_idx_base = widx * kBitsPerWordD2;

      // enable stronger loop optimizations
      const uint32_t difflist_idx_end = difflist_idx + PopcountWord(xor_word);
      while (1) {
        const uint32_t sample_idx_lowbits = CTZLU(xor_word) / 2;
        const uint32_t new_sample_idx = sample_idx_base + sample_idx_lowbits;
        if (!(difflist_idx % kBitsPerWordD2)) {
          if (!(difflist_idx % kPglDifflistGroupSize)) {
            SubUintStoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
            if (difflist_idx) {
              *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
              *raregeno_iter++ = raregeno_word;
              raregeno_word = 0;
            }
            last_group_vint_start = fwrite_bufp;
            goto SaveOnebit_skip_delta_write;
          }
          *raregeno_iter++ = raregeno_word;
          raregeno_word = 0;
        }
        assert(new_sample_idx >= last_sample_idx + 1);
        fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
      SaveOnebit_skip_delta_write:
        raregeno_word |= ((cur_geno_word >> (2 * sample_idx_lowbits)) & 3) << (2 * (difflist_idx % kBitsPerWordD2));
        ++difflist_idx;
        last_sample_idx = new_sample_idx;
        if (difflist_idx == difflist_idx_end) {
          break;
        }
        xor_word &= xor_word - 1;
      }
      // trailing bits of genovec_buf guaranteed to be zeroed out
      if (difflist_idx == onebit_difflist_len) {
        SubwordStore(raregeno_word, 1 + (((onebit_difflist_len - 1) / 4) % sizeof(intptr_t)), raregeno_iter);
        pwcp->fwrite_bufp = fwrite_bufp;
        return fwrite_bufp - fwrite_bufp_start;
      }
    }
    ++widx;
  }
}

uint32_t PwcAppendBiallelicGenovecMain(const uintptr_t* __restrict genovec, uint32_t vidx, PgenWriterCommon* pwcp, uint32_t* het_ct_ptr, unsigned char* vrtype_ptr) {
#ifndef NDEBUG
  if (pwcp->allele_idx_offsets) {
    assert(pwcp->allele_idx_offsets[vidx + 1] == pwcp->allele_idx_offsets[vidx] + 2);
  }
#endif
  const uint32_t sample_ct = pwcp->sample_ct;
  assert((!(sample_ct % kBitsPerWordD2)) || (!(genovec[sample_ct / kBitsPerWordD2] >> (2 * (sample_ct % kBitsPerWordD2)))));
  uint32_t genocounts[4];
  GenovecCountFreqsUnsafe(genovec, sample_ct, genocounts);
  *het_ct_ptr = genocounts[1];
  uint32_t most_common_geno = (genocounts[1] > genocounts[0]);
  uint32_t second_most_common_geno = 1 - most_common_geno;
  uint32_t largest_geno_ct = genocounts[most_common_geno];
  uint32_t second_largest_geno_ct = genocounts[second_most_common_geno];
  for (uint32_t cur_geno = 2; cur_geno < 4; ++cur_geno) {
    const uint32_t cur_geno_ct = genocounts[cur_geno];
    if (cur_geno_ct > second_largest_geno_ct) {
      if (cur_geno_ct > largest_geno_ct) {
        second_largest_geno_ct = largest_geno_ct;
        second_most_common_geno = most_common_geno;
        largest_geno_ct = cur_geno_ct;
        most_common_geno = cur_geno;
      } else {
        second_largest_geno_ct = cur_geno_ct;
        second_most_common_geno = cur_geno;
      }
    }
  }
  const uint32_t difflist_len = sample_ct - largest_geno_ct;
  const uint32_t rare_2_geno_ct_sum = difflist_len - second_largest_geno_ct;
  // average of 10-11 bits per difflist entry
  const uint32_t sample_ctd8 = sample_ct / 8;
  const uint32_t sample_ctd64 = sample_ct / 64;
  uint32_t max_difflist_len = sample_ctd8 - 2 * sample_ctd64 + rare_2_geno_ct_sum;
  if (max_difflist_len > sample_ctd8) {
    max_difflist_len = sample_ctd8;
  }
  const uint32_t difflist_viable = (most_common_geno != 1) && (difflist_len <= max_difflist_len);

  uintptr_t* ldbase_genovec = pwcp->ldbase_genovec;
  uint32_t* ldbase_genocounts = pwcp->ldbase_genocounts;
  if (!(vidx % kPglVblockSize)) {
    // beginning of a variant block.  save raw fpos in header; LD compression
    // prohibited.

    // er, need to use a relative offset in the multithreaded case, absolute
    // position isn't known
    pwcp->vblock_fpos[vidx / kPglVblockSize] = pwcp->vblock_fpos_offset + S_CAST(uintptr_t, pwcp->fwrite_bufp - pwcp->fwrite_buf);
  } else if (difflist_len > sample_ctd64) {
    // do not use LD compression if there are at least this many differences.
    // tune this threshold in the future.
    const uint32_t ld_diff_threshold = difflist_viable? (difflist_len - sample_ctd64) : max_difflist_len;
    // number of changes between current genovec and LD reference is bounded
    // below by sum(genocounts[x] - ldbase_genocounts[x]) / 2
    const int32_t count02_limit = 2 * ld_diff_threshold - abs_i32(genocounts[1] - ldbase_genocounts[1]) + abs_i32(genocounts[3] - ldbase_genocounts[3]);
    if ((S_CAST(int32_t, abs_i32(genocounts[0] - ldbase_genocounts[0]) + abs_i32(genocounts[2] - ldbase_genocounts[2])) < count02_limit) || (S_CAST(int32_t, abs_i32(genocounts[0] - ldbase_genocounts[2]) + abs_i32(genocounts[2] - ldbase_genocounts[0])) < count02_limit)) {
      uint32_t ld_diff_ct;
      uint32_t ld_inv_diff_ct;
      // okay, perform a brute-force diff
      // (could check LD vs. inverted LD separately?)
      if (pwcp->ldbase_common_geno < 4) {
        // unpack to ldbase_genovec
        PgrDifflistToGenovecUnsafe(pwcp->ldbase_raregeno, pwcp->ldbase_difflist_sample_ids, pwcp->ldbase_common_geno, sample_ct, pwcp->ldbase_difflist_len, ldbase_genovec);
        ZeroTrailingQuaters(sample_ct, ldbase_genovec);
        pwcp->ldbase_common_geno = UINT32_MAX;
      }
      CountLdAndInvertedLdDiffs(ldbase_genovec, genovec, sample_ct, &ld_diff_ct, &ld_inv_diff_ct);
      if ((ld_diff_ct < ld_diff_threshold) || (ld_inv_diff_ct < ld_diff_threshold)) {
        const uintptr_t invert_before_compressing = (ld_inv_diff_ct < ld_diff_ct);
        *vrtype_ptr = 2 + invert_before_compressing;
        if (invert_before_compressing) {
          GenovecInvertCopyUnsafe(genovec, sample_ct, pwcp->genovec_invert_buf);
          ld_diff_ct = ld_inv_diff_ct;
        }
        return SaveLdDifflist(invert_before_compressing? pwcp->genovec_invert_buf : genovec, ldbase_genovec, 0, ld_diff_ct, pwcp);
      }
    }
  }
  const uint32_t genovec_word_ct = QuaterCtToWordCt(sample_ct);
  memcpy(ldbase_genocounts, genocounts, 4 * sizeof(int32_t));
  pwcp->ldbase_common_geno = UINT32_MAX;
  if ((!difflist_viable) && (rare_2_geno_ct_sum < sample_ct / (2 * kPglMaxDifflistLenDivisor))) {
    *vrtype_ptr = 1;
    uint32_t larger_common_geno = second_most_common_geno;
    uint32_t smaller_common_geno = most_common_geno;
    if (most_common_geno > second_most_common_geno) {
      larger_common_geno = most_common_geno;
      smaller_common_geno = second_most_common_geno;
    }
    const uint32_t vrec_len = SaveOnebit(genovec, larger_common_geno + (smaller_common_geno * 3), rare_2_geno_ct_sum, pwcp);
    memcpy(ldbase_genovec, genovec, genovec_word_ct * sizeof(intptr_t));
    return vrec_len;
  }
  memcpy(ldbase_genovec, genovec, genovec_word_ct * sizeof(intptr_t));
  if (difflist_viable) {
    *vrtype_ptr = 4 + most_common_geno;
    return SaveLdDifflist(genovec, nullptr, most_common_geno, difflist_len, pwcp);
  }
  *vrtype_ptr = 0;
  const uint32_t vrec_len = QuaterCtToByteCt(sample_ct);
  pwcp->fwrite_bufp = R_CAST(unsigned char*, memcpya(pwcp->fwrite_bufp, genovec, vrec_len));
  return vrec_len;
}

void PwcAppendBiallelicGenovec(const uintptr_t* __restrict genovec, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  uint32_t het_ct;  // dummy
  unsigned char vrtype;
  const uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, &vrtype);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  pwcp->vidx += 1;
  SubUintStore(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * pwcp->vrec_len_byte_ct]));
  // could have a single expression which branchlessly handles both cases, but
  // doubt that's worthwhile
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
}

PglErr SpgwAppendBiallelicGenovec(const uintptr_t* __restrict genovec, STPgenWriter* spgwp) {
  // flush write buffer if necessary
  if (spgwp->pwc.fwrite_bufp >= &(spgwp->pwc.fwrite_buf[kPglFwriteBlockSize])) {
    const uintptr_t cur_byte_ct = spgwp->pwc.fwrite_bufp - spgwp->pwc.fwrite_buf;
    if (fwrite_checked(spgwp->pwc.fwrite_buf, cur_byte_ct, spgwp->pgen_outfile)) {
      return kPglRetWriteFail;
    }
    spgwp->pwc.vblock_fpos_offset += cur_byte_ct;
    spgwp->pwc.fwrite_bufp = spgwp->pwc.fwrite_buf;
  }
  PwcAppendBiallelicGenovec(genovec, &(spgwp->pwc));
  return kPglRetSuccess;
}

uint32_t SaveLdTwoListDelta(const uintptr_t* __restrict difflist_raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t ld_diff_ct, PgenWriterCommon* pwcp) {
  // assumes ldbase_difflist_sample_ids[ldbase_difflist_len] == sample_ct, and
  // difflist_sample_ids[ldbase_difflist_len] == sample_ct.
  // assumes biallelic data.

  // similar to SaveLdDifflist() and, to a lesser degree,
  // ParseLdAndMergeDifflistSubset()
  unsigned char* fwrite_bufp = pwcp->fwrite_bufp;
  if (!ld_diff_ct) {
    *fwrite_bufp = 0;
    pwcp->fwrite_bufp = &(fwrite_bufp[1]);
    return 1;
  }
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(ld_diff_ct, fwrite_bufp);
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(ld_diff_ct, kPglDifflistGroupSize);
  const uint32_t ldbase_common_geno = pwcp->ldbase_common_geno;
  assert(ldbase_common_geno < 4);
  const uintptr_t* __restrict ldbase_raregeno = pwcp->ldbase_raregeno;
  const uint32_t* __restrict ldbase_sample_ids = pwcp->ldbase_difflist_sample_ids;
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
#ifdef __arm__
#  error "Unaligned accesses in SaveLdTwoListDelta()."
#endif
  uintptr_t* raregeno_write_iter = R_CAST(uintptr_t*, &(extra_byte_cts_iter[group_ct - 1]));
  fwrite_bufp = &(extra_byte_cts_iter[group_ct + (ld_diff_ct - 1) / 4]);
  unsigned char* last_group_vint_start = fwrite_bufp;
  uintptr_t ldbase_raregeno_word = 0;
  uintptr_t difflist_raregeno_word = 0;
  uintptr_t raregeno_write_word = 0;
  uint32_t last_sample_idx = 0;

  uint32_t next_ldbase_sample_idx = ldbase_sample_ids[0];
  uint32_t next_difflist_sample_idx = difflist_sample_ids[0];
  uint32_t ldbase_idx = 0;
  uint32_t difflist_idx = 0;
  uint32_t diff_written_ct = 0;
  while (diff_written_ct < ld_diff_ct) {
    uintptr_t cur_geno;
    uint32_t new_sample_idx;
    if (next_ldbase_sample_idx <= next_difflist_sample_idx) {
      ldbase_raregeno_word >>= 2;
      if (!(ldbase_idx % kBitsPerWordD2)) {
        ldbase_raregeno_word = ldbase_raregeno[ldbase_idx / kBitsPerWordD2];
      }
      ++ldbase_idx;
    }
    if (next_difflist_sample_idx <= next_ldbase_sample_idx) {
      difflist_raregeno_word >>= 2;
      if (!(difflist_idx % kBitsPerWordD2)) {
        difflist_raregeno_word = difflist_raregeno[difflist_idx / kBitsPerWordD2];
      }
      new_sample_idx = next_difflist_sample_idx;
      ++difflist_idx;
      cur_geno = difflist_raregeno_word & 3;
      next_difflist_sample_idx = difflist_sample_ids[difflist_idx];
      if (next_ldbase_sample_idx == new_sample_idx) {
        next_ldbase_sample_idx = ldbase_sample_ids[ldbase_idx];
        if (cur_geno == (ldbase_raregeno_word & 3)) {
          continue;
        }
      }
    } else {
      cur_geno = ldbase_common_geno;
      new_sample_idx = next_ldbase_sample_idx;
      next_ldbase_sample_idx = ldbase_sample_ids[ldbase_idx];
    }
    raregeno_write_word |= cur_geno << (2 * (diff_written_ct % kBitsPerWordD2));
    if (!(diff_written_ct % kPglDifflistGroupSize)) {
      SubUintStoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
      if (diff_written_ct) {
        *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
      }
      last_group_vint_start = fwrite_bufp;
    } else {
      fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
    }
    last_sample_idx = new_sample_idx;
    ++diff_written_ct;
    if (!(diff_written_ct % kBitsPerWordD2)) {
      *raregeno_write_iter++ = raregeno_write_word;
      raregeno_write_word = 0;
    }
  }
  if (diff_written_ct % kBitsPerWordD2) {
    SubwordStore(raregeno_write_word, 1 + (((ld_diff_ct - 1) / 4) % kBytesPerWord), raregeno_write_iter);
  }
  pwcp->fwrite_bufp = fwrite_bufp;
  return fwrite_bufp - fwrite_bufp_start;
}

uint32_t SaveLdInputList(PgenWriterCommon* pwcp) {
  // simply "copies" ldbase_{raregeno,difflist_sample_ids,difflist_len} to the
  // write buffer.
  unsigned char* fwrite_bufp = pwcp->fwrite_bufp;
  const uint32_t difflist_len = pwcp->ldbase_difflist_len;
  if (!difflist_len) {
    *fwrite_bufp = 0;
    pwcp->fwrite_bufp = &(fwrite_bufp[1]);
    return 1;
  }
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(difflist_len, fwrite_bufp);
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t* __restrict difflist_sample_ids = pwcp->ldbase_difflist_sample_ids;
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
  fwrite_bufp = R_CAST(unsigned char*, memcpya(&(extra_byte_cts_iter[group_ct - 1]), pwcp->ldbase_raregeno, QuaterCtToByteCt(difflist_len)));
  unsigned char* last_group_vint_start = nullptr;
  uint32_t last_sample_idx = 0;
  for (uint32_t difflist_idx = 0; difflist_idx < difflist_len; ++difflist_idx) {
    const uint32_t new_sample_idx = difflist_sample_ids[difflist_idx];
    if (!(difflist_idx % kPglDifflistGroupSize)) {
      SubUintStoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
      if (difflist_idx) {
        *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
      }
      last_group_vint_start = fwrite_bufp;
    } else {
      // assert(new_sample_idx >= last_sample_idx + 1);
      fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
    }
    last_sample_idx = new_sample_idx;
  }
  pwcp->fwrite_bufp = fwrite_bufp;
  return fwrite_bufp - fwrite_bufp_start;
}

uint32_t PwcAppendBiallelicDifflistLimitedMain(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t vidx, uint32_t difflist_common_geno, uint32_t difflist_len, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr) {
  const uint32_t sample_ct = pwcp->sample_ct;
  // caller's responsibility not to exceed this limit
  assert(difflist_len <= 2 * (sample_ct / kPglMaxDifflistLenDivisor));

  // trailing bits of raregeno must be zeroed out

  assert(difflist_common_geno < 4);
#ifndef NDEBUG
  if (pwcp->allele_idx_offsets) {
    assert(pwcp->allele_idx_offsets[vidx + 1] == pwcp->allele_idx_offsets[vidx] + 2);
  }
#endif
  assert((!(difflist_len % kBitsPerWordD2)) || (!(raregeno[difflist_len / kBitsPerWordD2] >> (2 * (difflist_len % kBitsPerWordD2)))));
  assert(difflist_sample_ids[difflist_len] == sample_ct);
  uint32_t genocounts[4];
  GenovecCountFreqsUnsafe(raregeno, difflist_len, genocounts);
  assert(!genocounts[difflist_common_geno]);
  genocounts[difflist_common_geno] = sample_ct - difflist_len;
  uint32_t second_most_common_geno = difflist_common_geno? 0 : 1;
  uint32_t second_largest_geno_ct = genocounts[second_most_common_geno];
  for (uint32_t cur_geno = second_most_common_geno + 1; cur_geno < 4; ++cur_geno) {
    if (cur_geno == difflist_common_geno) {
      continue;
    }
    const uint32_t cur_geno_ct = genocounts[cur_geno];
    if (cur_geno_ct > second_largest_geno_ct) {
      second_most_common_geno = cur_geno;
      second_largest_geno_ct = cur_geno_ct;
    }
  }
  const uint32_t rare_2_geno_ct_sum = difflist_len - second_largest_geno_ct;
  const uint32_t sample_ctd8 = sample_ct / 8;
  const uint32_t sample_ctd64 = sample_ct / 64;
  uint32_t max_difflist_len = sample_ctd8 - 2 * sample_ctd64 + rare_2_geno_ct_sum;
  if (max_difflist_len > sample_ctd8) {
    max_difflist_len = sample_ctd8;
  }
  const uint32_t difflist_viable = (difflist_common_geno != 1) && (difflist_len <= max_difflist_len);
  uint32_t* ldbase_genocounts = pwcp->ldbase_genocounts;
  if (!(vidx % kPglVblockSize)) {
    pwcp->vblock_fpos[vidx / kPglVblockSize] = pwcp->vblock_fpos_offset + S_CAST(uintptr_t, pwcp->fwrite_bufp - pwcp->fwrite_buf);
  } else if (difflist_len > sample_ctd64) {
    const uint32_t ld_diff_threshold = difflist_viable? (difflist_len - sample_ctd64) : max_difflist_len;
    // number of changes between current genovec and LD reference is bounded
    // below by sum(genocounts[x] - ldbase_genocounts[x]) / 2
    const int32_t count02_limit = 2 * ld_diff_threshold - abs_i32(genocounts[1] - ldbase_genocounts[1]) + abs_i32(genocounts[3] - ldbase_genocounts[3]);
    if ((S_CAST(int32_t, abs_i32(genocounts[0] - ldbase_genocounts[0]) + abs_i32(genocounts[2] - ldbase_genocounts[2])) < count02_limit) || (S_CAST(int32_t, abs_i32(genocounts[0] - ldbase_genocounts[2]) + abs_i32(genocounts[2] - ldbase_genocounts[0])) < count02_limit)) {
      uint32_t ld_diff_ct;
      uint32_t ld_inv_diff_ct;
      if (pwcp->ldbase_common_geno < 4) {
        pwcp->ldbase_difflist_sample_ids[pwcp->ldbase_difflist_len] = sample_ct;
        if (CountLdAndInvertedLdDiffsList(pwcp->ldbase_raregeno, pwcp->ldbase_difflist_sample_ids, raregeno, difflist_sample_ids, pwcp->ldbase_difflist_len, difflist_len, &ld_diff_ct, &ld_inv_diff_ct)) {
          const uint32_t difflist_common_geno_inv = (6 - difflist_common_geno) & 3;
          if (pwcp->ldbase_common_geno != difflist_common_geno) {
            ld_diff_ct = ld_diff_threshold;
          }
          if (pwcp->ldbase_common_geno != difflist_common_geno_inv) {
            ld_inv_diff_ct = ld_diff_threshold;
          }
          if ((ld_diff_ct < ld_diff_threshold) || (ld_inv_diff_ct < ld_diff_threshold)) {
            const uintptr_t invert_before_compressing = (ld_inv_diff_ct < ld_diff_ct);
            *vrtype_ptr = 2 + invert_before_compressing;
            if (invert_before_compressing) {
              GenovecInvertCopyUnsafe(raregeno, difflist_len, pwcp->genovec_invert_buf);
              // difflist_common_geno = difflist_common_geno_inv;
              ld_diff_ct = ld_inv_diff_ct;
            }
            return SaveLdTwoListDelta(invert_before_compressing? pwcp->genovec_invert_buf : raregeno, difflist_sample_ids, ld_diff_ct, pwcp);
          }
        }
      } else {
        uintptr_t* __restrict genobuf = pwcp->genovec_invert_buf;
        PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, sample_ct, difflist_len, genobuf);
        ZeroTrailingQuaters(sample_ct, genobuf);
        CountLdAndInvertedLdDiffs(pwcp->ldbase_genovec, genobuf, sample_ct, &ld_diff_ct, &ld_inv_diff_ct);
        if ((ld_diff_ct < ld_diff_threshold) || (ld_inv_diff_ct < ld_diff_threshold)) {
          const uintptr_t invert_before_compressing = (ld_inv_diff_ct < ld_diff_ct);
          *vrtype_ptr = 2 + invert_before_compressing;
          if (invert_before_compressing) {
            GenovecInvertUnsafe(sample_ct, genobuf);
            ld_diff_ct = ld_inv_diff_ct;
          }
          return SaveLdDifflist(genobuf, pwcp->ldbase_genovec, 0, ld_diff_ct, pwcp);
        }
      }
    }
  }
  memcpy(ldbase_genocounts, genocounts, 4 * sizeof(int32_t));
  if (difflist_viable) {
    *vrtype_ptr = 4 + difflist_common_geno;
    memcpy(pwcp->ldbase_raregeno, raregeno, QuaterCtToByteCt(difflist_len));
    memcpy(pwcp->ldbase_difflist_sample_ids, difflist_sample_ids, difflist_len * sizeof(int32_t));
    // memcpy(pwcp->ldbase_difflist_sample_ids, difflist_sample_ids, (difflist_len + 1) * sizeof(int32_t));
    pwcp->ldbase_common_geno = difflist_common_geno;
    pwcp->ldbase_difflist_len = difflist_len;
    return SaveLdInputList(pwcp);
  }
  pwcp->ldbase_common_geno = UINT32_MAX;
  const uint32_t use_onebit = (rare_2_geno_ct_sum < sample_ct / (2 * kPglMaxDifflistLenDivisor));
  uintptr_t* genobuf = use_onebit? pwcp->genovec_invert_buf : pwcp->ldbase_genovec;
  PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, sample_ct, difflist_len, genobuf);
  ZeroTrailingQuaters(sample_ct, genobuf);
  *vrtype_ptr = use_onebit;
  if (use_onebit) {
    uint32_t larger_common_geno = second_most_common_geno;
    uint32_t smaller_common_geno = difflist_common_geno;
    if (difflist_common_geno > second_most_common_geno) {
      larger_common_geno = difflist_common_geno;
      smaller_common_geno = second_most_common_geno;
    }
    const uint32_t vrec_len = SaveOnebit(genobuf, larger_common_geno + (smaller_common_geno * 3), rare_2_geno_ct_sum, pwcp);
    memcpy(pwcp->ldbase_genovec, genobuf, QuaterCtToWordCt(sample_ct) * sizeof(uintptr_t));
    return vrec_len;
  }
  const uint32_t vrec_len = QuaterCtToByteCt(sample_ct);
  pwcp->fwrite_bufp = R_CAST(unsigned char*, memcpya(pwcp->fwrite_bufp, genobuf, vrec_len));
  return vrec_len;
}

void PwcAppendBiallelicDifflistLimited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  const uint32_t vrec_len = PwcAppendBiallelicDifflistLimitedMain(raregeno, difflist_sample_ids, vidx, difflist_common_geno, difflist_len, pwcp, &vrtype);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  pwcp->vidx += 1;
  SubUintStore(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * pwcp->vrec_len_byte_ct]));
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
}

PglErr SpgwAppendBiallelicDifflistLimited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, STPgenWriter* spgwp) {
  // trailing bits of raregeno must be zeroed out

  // flush write buffer if necessary
  if (spgwp->pwc.fwrite_bufp >= &(spgwp->pwc.fwrite_buf[kPglFwriteBlockSize])) {
    const uintptr_t cur_byte_ct = spgwp->pwc.fwrite_bufp - spgwp->pwc.fwrite_buf;
    if (fwrite_checked(spgwp->pwc.fwrite_buf, cur_byte_ct, spgwp->pgen_outfile)) {
      return kPglRetWriteFail;
    }
    spgwp->pwc.vblock_fpos_offset += cur_byte_ct;
    spgwp->pwc.fwrite_bufp = spgwp->pwc.fwrite_buf;
  }
  PwcAppendBiallelicDifflistLimited(raregeno, difflist_sample_ids, difflist_common_geno, difflist_len, &(spgwp->pwc));
  return kPglRetSuccess;
}


PglErr SpgwAppendMultiallelicCounts(__attribute__((unused)) const uintptr_t** __restrict alt_countvecs) {
  // todo
  return kPglRetNotYetSupported;
}


void AppendHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, uint32_t het_ct, uint32_t phasepresent_ct, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  assert(phasepresent_ct);
  const uint32_t sample_ct = pwcp->sample_ct;
  *vrtype_ptr += 16;
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
#ifdef __arm__
#  error "Unaligned accesses in AppendHphase()."
#endif
  uintptr_t* fwrite_bufp_alias = R_CAST(uintptr_t*, pwcp->fwrite_bufp);
  uintptr_t phaseinfo_write_word = 0;
  uint32_t phaseinfo_write_idx_lowbits;
  unsigned char* fwrite_bufp_final;
  if (het_ct == phasepresent_ct) {
    // no need to write phasepresent; just write phaseinfo directly to output
    // buffer
    phaseinfo_write_idx_lowbits = 1;
    for (uint32_t widx = 0; widx < sample_ctl2; ++widx) {
      const uintptr_t geno_word = genovec[widx];
      uintptr_t geno_hets = (~(geno_word >> 1)) & geno_word & kMask5555;
      if (geno_hets) {
        const uint32_t phaseinfo_halfword = R_CAST(const Halfword*, phaseinfo)[widx];
        do {
          const uint32_t sample_idx_lowbits = CTZLU(geno_hets) / 2;
          phaseinfo_write_word |= S_CAST(uintptr_t, (phaseinfo_halfword >> sample_idx_lowbits) & k1LU) << phaseinfo_write_idx_lowbits;
          if (++phaseinfo_write_idx_lowbits == kBitsPerWord) {
            *fwrite_bufp_alias++ = phaseinfo_write_word;
            phaseinfo_write_word = 0;
            phaseinfo_write_idx_lowbits = 0;
          }
          geno_hets &= geno_hets - 1;
        } while (geno_hets);
      }
    }
    fwrite_bufp_final = R_CAST(unsigned char*, fwrite_bufp_alias);
  } else {
    uintptr_t* phaseinfo_tmp = pwcp->genovec_invert_buf;
    uintptr_t* phaseinfo_tmp_iter = phaseinfo_tmp;
    uint32_t phasepresent_write_idx_lowbits = 1;
    phaseinfo_write_idx_lowbits = 0;
    uintptr_t phasepresent_write_word = 1;
    for (uint32_t widx = 0; widx < sample_ctl2; ++widx) {
      const uintptr_t geno_word = genovec[widx];
      uintptr_t geno_hets = (~(geno_word >> 1)) & geno_word & kMask5555;
      if (geno_hets) {
        const uint32_t phasepresent_halfword = R_CAST(const Halfword*, phasepresent)[widx];
        if (phasepresent_halfword) {
          const uint32_t phaseinfo_halfword = R_CAST(const Halfword*, phaseinfo)[widx];
          do {
            const uint32_t sample_idx_lowbits = CTZLU(geno_hets) / 2;
            if ((phasepresent_halfword >> sample_idx_lowbits) & 1) {
              phasepresent_write_word |= k1LU << phasepresent_write_idx_lowbits;
              phaseinfo_write_word |= S_CAST(uintptr_t, (phaseinfo_halfword >> sample_idx_lowbits) & k1LU) << phaseinfo_write_idx_lowbits;
              if (++phaseinfo_write_idx_lowbits == kBitsPerWord) {
                *phaseinfo_tmp_iter++ = phaseinfo_write_word;
                phaseinfo_write_word = 0;
                phaseinfo_write_idx_lowbits = 0;
              }
            }
            if (++phasepresent_write_idx_lowbits == kBitsPerWord) {
              *fwrite_bufp_alias++ = phasepresent_write_word;
              phasepresent_write_word = 0;
              phasepresent_write_idx_lowbits = 0;
            }
            geno_hets &= geno_hets - 1;
          } while (geno_hets);
        } else {
          phasepresent_write_idx_lowbits += PopcountWord(geno_hets);
          if (phasepresent_write_idx_lowbits >= kBitsPerWord) {
            *fwrite_bufp_alias++ = phasepresent_write_word;
            phasepresent_write_word = 0;
            phasepresent_write_idx_lowbits -= kBitsPerWord;
          }
        }
      }
    }
    fwrite_bufp_final = R_CAST(unsigned char*, fwrite_bufp_alias);
    if (phasepresent_write_idx_lowbits) {
      const uint32_t cur_byte_ct = DivUp(phasepresent_write_idx_lowbits, CHAR_BIT);
      SubwordStoreMov(phasepresent_write_word, cur_byte_ct, &fwrite_bufp_final);
    }
    fwrite_bufp_final = R_CAST(unsigned char*, memcpya(fwrite_bufp_final, phaseinfo_tmp, sizeof(intptr_t) * (phaseinfo_tmp_iter - phaseinfo_tmp)));
  }
  if (phaseinfo_write_idx_lowbits) {
    const uint32_t cur_byte_ct = DivUp(phaseinfo_write_idx_lowbits, CHAR_BIT);
    SubwordStoreMov(phaseinfo_write_word, cur_byte_ct, &fwrite_bufp_final);
  }
#ifdef __LP64__
  assert(((*vrec_len_ptr) + S_CAST(uintptr_t, fwrite_bufp_final - pwcp->fwrite_bufp)) <= kPglMaxBytesPerVariant);
#endif
  *vrec_len_ptr += fwrite_bufp_final - pwcp->fwrite_bufp;
  pwcp->fwrite_bufp = fwrite_bufp_final;
}

void PwcAppendBiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, PgenWriterCommon* pwcp) {
  // assumes phase_dosage_gflags is nonzero
  const uint32_t vidx = pwcp->vidx;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  uint32_t het_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, vrtype_dest);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  const uint32_t phasepresent_ct = phasepresent? PopcountWords(phasepresent, sample_ctl) : het_ct;
  if (phasepresent_ct) {
    AppendHphase(genovec, phasepresent, phaseinfo, het_ct, phasepresent_ct, pwcp, vrtype_dest, &vrec_len);
  }
  SubUintStore(vrec_len, vrec_len_byte_ct, vrec_len_dest);
}

PglErr SpgwAppendBiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, STPgenWriter* spgwp) {
  // flush write buffer if necessary
  if (spgwp->pwc.fwrite_bufp >= &(spgwp->pwc.fwrite_buf[kPglFwriteBlockSize])) {
    const uintptr_t cur_byte_ct = spgwp->pwc.fwrite_bufp - spgwp->pwc.fwrite_buf;
    if (fwrite_checked(spgwp->pwc.fwrite_buf, cur_byte_ct, spgwp->pgen_outfile)) {
      return kPglRetWriteFail;
    }
    spgwp->pwc.vblock_fpos_offset += cur_byte_ct;
    spgwp->pwc.fwrite_bufp = spgwp->pwc.fwrite_buf;
  }
  PwcAppendBiallelicGenovecHphase(genovec, phasepresent, phaseinfo, &(spgwp->pwc));
  return kPglRetSuccess;
}


uint32_t PwcAppendDeltalist(const uintptr_t* delta_bitarr, uint32_t deltalist_len, PgenWriterCommon* pwcp) {
  assert(deltalist_len);
  unsigned char* fwrite_bufp = pwcp->fwrite_bufp;
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(deltalist_len, fwrite_bufp);
  const uint32_t sample_id_byte_ct = BytesToRepresentUi(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(deltalist_len, kPglDifflistGroupSize);
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
  fwrite_bufp = &(extra_byte_cts_iter[group_ct - 1]);
  unsigned char* last_group_vint_start = nullptr;
  uint32_t last_sample_idx = 0;
  uint32_t new_sample_idx = 0;
  for (uint32_t deltalist_idx = 0; deltalist_idx < deltalist_len; ++deltalist_idx, ++new_sample_idx) {
    FindFirst1BitFromU32(delta_bitarr, &new_sample_idx);
    if (!(deltalist_idx % kPglDifflistGroupSize)) {
      SubUintStoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
      if (deltalist_idx) {
        *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
      }
      last_group_vint_start = fwrite_bufp;
    } else {
      assert(new_sample_idx >= last_sample_idx + 1);
      fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
    }
    last_sample_idx = new_sample_idx;
  }
  pwcp->fwrite_bufp = fwrite_bufp;
  return fwrite_bufp - fwrite_bufp_start;
}

void AppendDosage16(const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t max_deltalist_entry_ct = sample_ct / kPglMaxDeltalistLenDivisor;
  if (dosage_ct <= max_deltalist_entry_ct) {
    // case 1: store dosage IDs as deltalist.
    *vrec_len_ptr += PwcAppendDeltalist(dosage_present, dosage_ct, pwcp);
    *vrtype_ptr += 0x20;
  } else if (dosage_ct == sample_ct) {
    // case 2: fixed-width, no need to store dosage IDs at all.
    // dosage_vals permitted to have 65535 = missing
    *vrtype_ptr += 0x40;
  } else {
    // case 3: save dosage_present bitarray directly.
    const uint32_t sample_ctb = DivUp(sample_ct, CHAR_BIT);
    *vrec_len_ptr += sample_ctb;
    pwcp->fwrite_bufp = R_CAST(unsigned char*, memcpya(pwcp->fwrite_bufp, dosage_present, sample_ctb));
    *vrtype_ptr += 0x60;
  }
  pwcp->fwrite_bufp = R_CAST(unsigned char*, memcpya(pwcp->fwrite_bufp, dosage_vals, dosage_ct * sizeof(int16_t)));
  *vrec_len_ptr += dosage_ct * sizeof(int16_t);
}

void PwcAppendBiallelicGenovecDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, PgenWriterCommon* pwcp) {
  // safe to call this even when entire file has no phase/dosage info
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  uint32_t het_ct;  // dummy
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, &vrtype);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  if (dosage_ct) {
    AppendDosage16(dosage_present, dosage_vals, dosage_ct, pwcp, &vrtype, &vrec_len);
  }
  SubUintStore(vrec_len, vrec_len_byte_ct, vrec_len_dest);
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
}

PglErr SpgwAppendBiallelicGenovecDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, STPgenWriter* spgwp) {
  // flush write buffer if necessary
  if (spgwp->pwc.fwrite_bufp >= &(spgwp->pwc.fwrite_buf[kPglFwriteBlockSize])) {
    const uintptr_t cur_byte_ct = spgwp->pwc.fwrite_bufp - spgwp->pwc.fwrite_buf;
    if (fwrite_checked(spgwp->pwc.fwrite_buf, cur_byte_ct, spgwp->pgen_outfile)) {
      return kPglRetWriteFail;
    }
    spgwp->pwc.vblock_fpos_offset += cur_byte_ct;
    spgwp->pwc.fwrite_bufp = spgwp->pwc.fwrite_buf;
  }
  PwcAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_vals, dosage_ct, &(spgwp->pwc));
  return kPglRetSuccess;
}

void PwcAppendBiallelicGenovecHphaseDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, PgenWriterCommon* pwcp) {
  // assumes there is phase and/or dosage data in output file, otherwise
  // vrtype_dest needs to be replaced

  // this mostly overlaps with PwcAppendBiallelicGenovecHphase(); probably
  // get rid of the latter
  const uint32_t vidx = pwcp->vidx;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  uint32_t het_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, vrtype_dest);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  const uint32_t phasepresent_ct = phasepresent? PopcountWords(phasepresent, sample_ctl) : het_ct;
  if (phasepresent_ct) {
    AppendHphase(genovec, phasepresent, phaseinfo, het_ct, phasepresent_ct, pwcp, vrtype_dest, &vrec_len);
  }
  if (dosage_ct) {
    AppendDosage16(dosage_present, dosage_vals, dosage_ct, pwcp, vrtype_dest, &vrec_len);
  }
  SubUintStore(vrec_len, vrec_len_byte_ct, vrec_len_dest);
}

PglErr SpgwAppendBiallelicGenovecHphaseDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, STPgenWriter* spgwp) {
  // flush write buffer if necessary
  if (spgwp->pwc.fwrite_bufp >= &(spgwp->pwc.fwrite_buf[kPglFwriteBlockSize])) {
    const uintptr_t cur_byte_ct = spgwp->pwc.fwrite_bufp - spgwp->pwc.fwrite_buf;
    if (fwrite_checked(spgwp->pwc.fwrite_buf, cur_byte_ct, spgwp->pgen_outfile)) {
      return kPglRetWriteFail;
    }
    spgwp->pwc.vblock_fpos_offset += cur_byte_ct;
    spgwp->pwc.fwrite_bufp = spgwp->pwc.fwrite_buf;
  }
  PwcAppendBiallelicGenovecHphaseDosage16(genovec, phasepresent, phaseinfo, dosage_present, dosage_vals, dosage_ct, &(spgwp->pwc));
  return kPglRetSuccess;
}

/*
// need to rewrite these to reflect reordering of tracks #4-5, and addition of
// track #6.
void AppendDphase16(const uintptr_t* __restrict dosage_present, const uintptr_t* __restrict dphase_present, const uint16_t* __restrict dosage_vals, uint32_t dosage_ct, uint32_t dphase_ct, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  if (!dphase_ct) {
    AppendDosage16(dosage_present, dosage_vals, dosage_ct, pwcp, vrtype_ptr, vrec_len_ptr);
  }
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t max_deltalist_entry_ct = sample_ct / kPglMaxDeltalistLenDivisor;
  if (dosage_ct <= max_deltalist_entry_ct) {
    // case 1: store dosage IDs as deltalist.
    *vrec_len_ptr += PwcAppendDeltalist(dosage_present, dosage_ct, pwcp);
    *vrtype_ptr += 0x20;
  } else if (dosage_ct == sample_ct) {
    // case 2: fixed-width, no need to store dosage IDs at all.
    // dosage_vals permitted to have 65535 = missing
    *vrtype_ptr += 0x40;
  } else {
    // case 3: save dosage_present bitarray directly.
    const uint32_t sample_ctb = DivUp(sample_ct, CHAR_BIT);
    *vrec_len_ptr += sample_ctb;
    pwcp->fwrite_bufp = (unsigned char*)memcpya(pwcp->fwrite_bufp, dosage_present, sample_ctb);
    *vrtype_ptr += 0x60;
  }
  *vrtype_ptr += 0x80;
  if (dosage_ct == dphase_ct) {
    *(pwcp->fwrite_bufp)++ = 0;
    pwcp->fwrite_bufp = (unsigned char*)memcpya(pwcp->fwrite_bufp, dosage_vals, dphase_ct * 2 * sizeof(int16_t));
    *vrec_len_ptr += 1 + (dphase_ct * 2 * sizeof(int16_t));
  } else {
    uintptr_t* dphase_present_tmp_write_iter = pwcp->genovec_invert_buf;
    const uint32_t dosage_ctp1b = 1 + (dosage_ct / CHAR_BIT);
    const uint32_t widx_last = dosage_ct / kBitsPerWord;
    uintptr_t dphase_present_write_word = 1;
    uint32_t sample_idx = 0;
    uint32_t dosage_idx_lowbits = 1;
    uint32_t widx = 0;
    uint32_t loop_end = kBitsPerWord;
    while (1) {
      if (widx >= widx_last) {
        if (widx > widx_last) {
          break;
        }
        loop_end = 1 + (dosage_ct % kBitsPerWord);
      }
      for (; dosage_idx_lowbits < loop_end; ++dosage_idx_lowbits, ++sample_idx) {
        FindFirst1BitFromU32(dosage_present, &sample_idx);
        if (IsSet(dphase_present, sample_idx)) {
          dphase_present_write_word |= k1LU << dosage_idx_lowbits;
        }
      }
      *dphase_present_tmp_write_iter++ = dphase_present_write_word;
      dphase_present_write_word = 0;
      dosage_idx_lowbits = 0;
      ++widx;
    }
    char* cur_write_iter = memcpya(pwcp->fwrite_bufp, pwcp->genovec_invert_buf, dosage_ctp1b);
    cur_write_iter = memcpya(cur_write_iter, dosage_vals, (dosage_ct + dphase_ct) * sizeof(int16_t));
    *vrec_len_ptr += (uintptr_t)(cur_write_iter - pwcp->fwrite_bufp);
    pwcp->fwrite_bufp = (unsigned char*)cur_write_iter;
  }
}

void PwcAppendBiallelicGenovecDphase16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uintptr_t* __restrict dphase_present, const uint16_t* __restrict dosage_vals, uint32_t dosage_ct, uint32_t dphase_ct, PgenWriterCommon* pwcp) {
  // assumes there is phase and/or dosage data in output file, otherwise
  // vrtype_dest needs to be replaced
  const uint32_t vidx = pwcp->vidx;
  unsigned char* vrtype_dest = &(((unsigned char*)pwcp->vrtype_buf)[vidx]);
  uint32_t het_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, vrtype_dest);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  const uint32_t phasepresent_ct = phasepresent? PopcountWords(phasepresent, sample_ctl) : het_ct;
  if (phasepresent_ct) {
    AppendHphase(genovec, phasepresent, phaseinfo, het_ct, phasepresent_ct, pwcp, vrtype_dest, &vrec_len);
  }
  if (dosage_ct) {
    AppendDphase16(dosage_present, dphase_present, dosage_vals, dosage_ct, dphase_ct, pwcp, vrtype_dest, &vrec_len);
  }
  memcpy(vrec_len_dest, &vrec_len, vrec_len_byte_ct);
}

PglErr SpgwAppendBiallelicGenovecDphase16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uintptr_t* dphase_present, const uint16_t* dosage_vals, uint32_t dosage_ct, uint32_t dphase_ct, STPgenWriter* spgwp) {
  // flush write buffer if necessary
  if (spgwp->pwc.fwrite_bufp >= &(spgwp->pwc.fwrite_buf[kPglFwriteBlockSize])) {
    const uintptr_t cur_byte_ct = (uintptr_t)(spgwp->pwc.fwrite_bufp - spgwp->pwc.fwrite_buf);
    if (fwrite_checked(spgwp->pwc.fwrite_buf, cur_byte_ct, spgwp->pgen_outfile)) {
      return kPglRetWriteFail;
    }
    spgwp->pwc.vblock_fpos_offset += cur_byte_ct;
    spgwp->pwc.fwrite_bufp = spgwp->pwc.fwrite_buf;
  }
  PwcAppendBiallelicGenovecDphase16(genovec, phasepresent, phaseinfo, dosage_present, dphase_present, dosage_vals, dosage_ct, dphase_ct, &(spgwp->pwc));
  return kPglRetSuccess;
}
*/

PglErr PwcFinish(PgenWriterCommon* pwcp, FILE** pgen_outfile_ptr) {
  const uint32_t variant_ct = pwcp->variant_ct;
  assert(pwcp->vidx == variant_ct);
  FILE* pgen_outfile = *pgen_outfile_ptr;
  if (fseeko(pgen_outfile, 12, SEEK_SET)) {
    return kPglRetWriteFail;
  }
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  fwrite_unlocked(pwcp->vblock_fpos, vblock_ct * sizeof(int64_t), 1, pgen_outfile);
  const unsigned char* vrtype_buf_iter = R_CAST(unsigned char*, pwcp->vrtype_buf);
  const uint32_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  const unsigned char* vrec_len_buf_iter = pwcp->vrec_len_buf;
  const PgenGlobalFlags phase_dosage_gflags = pwcp->phase_dosage_gflags;
  uint32_t vrec_iter_incr = kPglVblockSize * vrec_len_byte_ct;
  uint32_t vrtype_buf_iter_incr = phase_dosage_gflags? kPglVblockSize : (kPglVblockSize / 2);
  uint32_t nonref_flags_write_byte_ct = kPglVblockSize / CHAR_BIT;
  const unsigned char* vrec_len_buf_last = &(vrec_len_buf_iter[S_CAST(uintptr_t, vblock_ct - 1) * vrec_iter_incr]);
  uintptr_t* explicit_nonref_flags = pwcp->explicit_nonref_flags;
  uintptr_t* explicit_nonref_flags_iter = explicit_nonref_flags;
  while (1) {
    if (vrec_len_buf_iter >= vrec_len_buf_last) {
      if (vrec_len_buf_iter > vrec_len_buf_last) {
        return fclose_null(pgen_outfile_ptr)? kPglRetWriteFail : kPglRetSuccess;
      }
      const uint32_t vblock_size = ModNz(variant_ct, kPglVblockSize);
      vrtype_buf_iter_incr = phase_dosage_gflags? vblock_size : DivUp(vblock_size, 2);
      vrec_iter_incr = vblock_size * vrec_len_byte_ct;
      nonref_flags_write_byte_ct = DivUp(vblock_size, CHAR_BIT);
    }
    // 4b(i): array of 4-bit or 1-byte vrtypes
    fwrite_unlocked(vrtype_buf_iter, vrtype_buf_iter_incr, 1, pgen_outfile);
    vrtype_buf_iter = &(vrtype_buf_iter[vrtype_buf_iter_incr]);

    // 4b(ii): array of variant record lengths
    if (fwrite_checked(vrec_len_buf_iter, vrec_iter_incr, pgen_outfile)) {
      return kPglRetWriteFail;
    }
    vrec_len_buf_iter = &(vrec_len_buf_iter[vrec_iter_incr]);

    // 4b(iii): alt allele counts
    // not yet supported

    // 4b(iv): explicit nonref flags
    if (explicit_nonref_flags) {
      if (fwrite_checked(explicit_nonref_flags_iter, nonref_flags_write_byte_ct, pgen_outfile)) {
        return kPglRetWriteFail;
      }
      explicit_nonref_flags_iter = &(explicit_nonref_flags_iter[kPglVblockSize / kBitsPerWord]);
    }
  }
}

PglErr SpgwFinish(STPgenWriter* spgwp) {
  if (fwrite_checked(spgwp->pwc.fwrite_buf, spgwp->pwc.fwrite_bufp - spgwp->pwc.fwrite_buf, spgwp->pgen_outfile)) {
    return kPglRetWriteFail;
  }
  return PwcFinish(&(spgwp->pwc), &(spgwp->pgen_outfile));
}

PglErr MpgwFlush(MTPgenWriter* mpgwp) {
  PgenWriterCommon* pwcp = mpgwp->pwcs[0];
  uint32_t vidx = RoundDownPow2(pwcp->vidx - 1, kPglVblockSize);
  uint32_t thread_ct = mpgwp->thread_ct;
  const uint32_t variant_ct = pwcp->variant_ct;
  const uint32_t is_last_flush = ((vidx + thread_ct * kPglVblockSize) >= variant_ct);
  if (is_last_flush) {
    thread_ct = DivUp(variant_ct - vidx, kPglVblockSize);
  }
  uint64_t* vblock_fpos = pwcp->vblock_fpos;
  FILE* pgen_outfile = mpgwp->pgen_outfile;
  const uint32_t vidx_incr = (thread_ct - 1) * kPglVblockSize;
  uint64_t cur_vblock_fpos = ftello(pgen_outfile);
  for (uint32_t tidx = 0; tidx < thread_ct; ++tidx) {
    vblock_fpos[(vidx / kPglVblockSize) + tidx] = cur_vblock_fpos;
    PgenWriterCommon* cur_pwcp = mpgwp->pwcs[tidx];
    uintptr_t cur_vblock_byte_ct = cur_pwcp->fwrite_bufp - cur_pwcp->fwrite_buf;
    if (fwrite_checked(cur_pwcp->fwrite_buf, cur_vblock_byte_ct, pgen_outfile)) {
      return kPglRetWriteFail;
    }
    cur_pwcp->vidx += vidx_incr;
    cur_pwcp->fwrite_bufp = cur_pwcp->fwrite_buf;
    cur_vblock_fpos += cur_vblock_byte_ct;
  }
  if (!is_last_flush) {
    return kPglRetSuccess;
  }
  pwcp->vidx = variant_ct;
  return PwcFinish(pwcp, &(mpgwp->pgen_outfile));
}

BoolErr SpgwCleanup(STPgenWriter* spgwp) {
  // assume file is open if spgw.pgen_outfile is not null
  // memory is the responsibility of the caller for now
  if (!spgwp->pgen_outfile) {
    return 0;
  }
  return fclose_null(&(spgwp->pgen_outfile));
}

BoolErr MpgwCleanup(MTPgenWriter* mpgwp) {
  if ((!mpgwp) || (!mpgwp->pgen_outfile)) {
    return 0;
  }
  return fclose_null(&(mpgwp->pgen_outfile));
}

#ifdef __cplusplus
}  // namespace plink2
#endif
