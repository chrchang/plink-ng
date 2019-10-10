// This library is part of PLINK 2.00, copyright (C) 2005-2019 Shaun Purcell,
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

#include <errno.h>

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
  for (uint32_t subset_mask_widx = 0; ; ++subset_mask_widx) {
    const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
    if (cur_include_word) {
      uint32_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
      for (uint32_t wordhalf_idx = 0; ; ++wordhalf_idx) {
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
  }
}

// bit_idx_start assumed to be < kBitsPerWord
void CopyGenomatchSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict genovec, uintptr_t match_word, uint32_t write_bit_idx_start, uint32_t bit_ct, uintptr_t* __restrict output_bitarr) {
  const uint32_t bit_idx_end = bit_ct + write_bit_idx_start;
  const uint32_t bit_idx_end_lowbits = bit_idx_end % kBitsPerWord;
  const Halfword* raw_bitarr_alias = R_CAST(const Halfword*, raw_bitarr);
  uintptr_t* output_bitarr_iter = output_bitarr;
  uintptr_t* output_bitarr_last = &(output_bitarr[bit_idx_end / kBitsPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = write_bit_idx_start;
  while ((output_bitarr_iter != output_bitarr_last) || (write_idx_lowbits != bit_idx_end_lowbits)) {
    uintptr_t cur_mask_word;
    // sparse genovec optimization
    // guaranteed to terminate since there's at least one more set bit
    do {
      // todo: try reading two genovec words at a time.  would need to be very
      // careful with the possible trailing word, though.
      // more important to optimize this function now that regular phased-call
      // handling code is using it.
      cur_mask_word = genovec[++read_widx] ^ match_word;
      cur_mask_word = (~(cur_mask_word | (cur_mask_word >> 1))) & kMask5555;
    } while (!cur_mask_word);
    uintptr_t extracted_bits = raw_bitarr_alias[read_widx];
    uint32_t set_bit_ct = kBitsPerWordD2;
    if (cur_mask_word != kMask5555) {
      const uintptr_t cur_mask_hw = PackWordToHalfword(cur_mask_word);
      set_bit_ct = PopcountWord(cur_mask_word);
      extracted_bits = _pext_u64(extracted_bits, cur_mask_hw);
    }
    cur_output_word |= extracted_bits << write_idx_lowbits;
    const uint32_t new_write_idx_lowbits = write_idx_lowbits + set_bit_ct;
    if (new_write_idx_lowbits >= kBitsPerWord) {
      *output_bitarr_iter++ = cur_output_word;
      // ...and these are the bits that fell off
      // impossible for write_idx_lowbits to be zero here
      cur_output_word = extracted_bits >> (kBitsPerWord - write_idx_lowbits);
    }
    write_idx_lowbits = new_write_idx_lowbits % kBitsPerWord;
  }
  if (write_idx_lowbits) {
    *output_bitarr_iter = cur_output_word;
  }
}

// Variant of ExpandBytearr() which is based off a target 2-bit value instead
// of single expand_mask bits.  expand_size must be the number of instance of
// the target value in genovec.
void ExpandBytearrFromGenovec(const void* __restrict compact_bitarr, const uintptr_t* __restrict genovec, uintptr_t match_word, uint32_t genoword_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t leading_byte_ct = 1 + (expand_sizex_m1 % kBitsPerWord) / CHAR_BIT;
  const uint32_t genoword_ct_m1 = genoword_ct - 1;
  uintptr_t compact_word = SubwordLoad(compact_bitarr, leading_byte_ct) >> read_start_bit;
  const uintptr_t* compact_bitarr_iter = R_CAST(const uintptr_t*, &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]));
  uint32_t compact_idx_lowbits = read_start_bit + CHAR_BIT * (sizeof(intptr_t) - leading_byte_ct);
  for (uint32_t widx = 0; ; widx += 2) {
    uintptr_t mask_word;
    if (widx >= genoword_ct_m1) {
      if (widx > genoword_ct_m1) {
        return;
      }
      mask_word = 0;
    } else {
      const uintptr_t geno_word1 = genovec[widx + 1] ^ match_word;
      mask_word = PackWordToHalfwordMask5555(~(geno_word1 | (geno_word1 >> 1)));
      mask_word = mask_word << 32;
    }
    const uintptr_t geno_word0 = genovec[widx] ^ match_word;
    mask_word |= PackWordToHalfwordMask5555(~(geno_word0 | (geno_word0 >> 1)));
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
#  ifdef __arm__
#    error "Unaligned accesses in ExpandBytearrFromGenovec()."
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
    target[widx / 2] = write_word;
  }
}
#else  // !USE_AVX2
void CopyQuaterarrNonemptySubset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_quaterarr) {
  if (subset_entry_ct == raw_quaterarr_entry_ct) {
    // subset_mask may be nullptr in this case
    memcpy(output_quaterarr, raw_quaterarr, DivUp(subset_entry_ct, kBitsPerWordD2) * sizeof(intptr_t));
    ZeroTrailingQuaters(subset_entry_ct, output_quaterarr);
    return;
  }
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
    for (uint32_t subset_mask_widx = 0; ; ++subset_mask_widx) {
      const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
      if (cur_include_word) {
        uint32_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
        for (uint32_t wordhalf_idx = 0; ; ++wordhalf_idx) {
          if (cur_include_halfword) {
            uintptr_t raw_quaterarr_word = raw_quaterarr[subset_mask_widx * 2 + wordhalf_idx];
            do {
              uint32_t rqa_idx_lowbits = ctzu32(cur_include_halfword);
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
  }
  // blocked copy
  const uintptr_t* raw_quaterarr_iter = raw_quaterarr;
  for (; ; ++subset_mask) {
    const uintptr_t cur_include_word = *subset_mask;
    uintptr_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
    for (uint32_t wordhalf_idx = 0; ; ++wordhalf_idx) {
      uintptr_t raw_quaterarr_word = *raw_quaterarr_iter++;
      while (cur_include_halfword) {
        const uint32_t rqa_idx_lowbits = ctzw(cur_include_halfword);

        // TAOCP, 7.1.3, (43).
        const uintptr_t bottom_block_remover = (cur_include_halfword | (cur_include_halfword - 1)) + 1;

        const uintptr_t raw_quaterarr_curblock_unmasked = raw_quaterarr_word >> (rqa_idx_lowbits * 2);
        const uint32_t rqa_block_len = ctzw(bottom_block_remover) - rqa_idx_lowbits;
        const uint32_t block_len_limit = kBitsPerWordD2 - word_write_halfshift;
        cur_output_word |= raw_quaterarr_curblock_unmasked << (2 * word_write_halfshift);
        if (rqa_block_len < block_len_limit) {
          word_write_halfshift += rqa_block_len;
          cur_output_word = bzhi(cur_output_word, word_write_halfshift * 2);
        } else {
          // no need to mask, extra bits vanish off the high end
          *output_quaterarr_iter++ = cur_output_word;
          word_write_halfshift = rqa_block_len - block_len_limit;
          if (word_write_halfshift) {
            cur_output_word = bzhi(raw_quaterarr_curblock_unmasked >> (2 * block_len_limit), 2 * word_write_halfshift);
          } else {
            // avoid potential right-shift-[word length]
            cur_output_word = 0;
          }
        }
        cur_include_halfword &= bottom_block_remover;
      }
      if (wordhalf_idx) {
        break;
      }
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

void CopyGenomatchSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict genovec, uintptr_t match_word, uint32_t write_bit_idx_start, uint32_t bit_ct, uintptr_t* __restrict output_bitarr) {
  const uint32_t bit_idx_end = write_bit_idx_start + bit_ct;
  const uint32_t bit_idx_end_lowbits = bit_idx_end % kBitsPerWord;
  const Halfword* raw_bitarr_alias = R_CAST(const Halfword*, raw_bitarr);
  uintptr_t* output_bitarr_iter = output_bitarr;
  uintptr_t* output_bitarr_last = &(output_bitarr[bit_idx_end / kBitsPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = write_bit_idx_start;
  while ((output_bitarr_iter != output_bitarr_last) || (write_idx_lowbits != bit_idx_end_lowbits)) {
    uintptr_t geno_word;
    // sparse genovec optimization
    // guaranteed to terminate since there's at least one more set bit
    do {
      geno_word = genovec[++read_widx] ^ match_word;
      geno_word = (~(geno_word | (geno_word >> 1))) & kMask5555;
    } while (!geno_word);
    // screw it, just iterate over set bits
    const uint32_t cur_halfword = raw_bitarr_alias[read_widx];
    do {
      const uint32_t sample_idx_lowbits = ctzw(geno_word) / 2;
      cur_output_word |= S_CAST(uintptr_t, (cur_halfword >> sample_idx_lowbits) & k1LU) << write_idx_lowbits;
      if (++write_idx_lowbits == kBitsPerWord) {
        *output_bitarr_iter++ = cur_output_word;
        cur_output_word = 0;
        write_idx_lowbits = 0;
      }
      geno_word &= geno_word - 1;
    } while (geno_word);
  }
  if (write_idx_lowbits) {
    *output_bitarr_iter = cur_output_word;
  }
}

void ExpandBytearrFromGenovec(const void* __restrict compact_bitarr, const uintptr_t* __restrict genovec, uintptr_t match_word, uint32_t genoword_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  Halfword* target_alias = R_CAST(Halfword*, target);
  ZeroHwArr(RoundUpPow2(genoword_ct, 2), target_alias);
  const uintptr_t* compact_bitarr_alias = S_CAST(const uintptr_t*, compact_bitarr);
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t compact_widx_last = expand_sizex_m1 / kBitsPerWord;
  uint32_t compact_idx_lowbits = read_start_bit;
  uint32_t loop_len = kBitsPerWord;
  uintptr_t write_hwidx = 0;
  uintptr_t genomatch_bits = genovec[0] ^ match_word;
  genomatch_bits = (~(genomatch_bits | (genomatch_bits >> 1))) & kMask5555;
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
#ifdef __arm__
#  error "Unaligned accesses in ExpandBytearrFromGenovec()."
#endif
      compact_word = compact_bitarr_alias[compact_widx];
    }
    for (; compact_idx_lowbits != loop_len; ++compact_idx_lowbits) {
      while (!genomatch_bits) {
        genomatch_bits = genovec[++write_hwidx] ^ match_word;
        genomatch_bits = (~(genomatch_bits | (genomatch_bits >> 1))) & kMask5555;
      }
      if (compact_word & (k1LU << compact_idx_lowbits)) {
        const uint32_t lowbit_idx = ctzw(genomatch_bits);
        target_alias[write_hwidx] |= 1U << (lowbit_idx / 2);
      }
      genomatch_bits &= genomatch_bits - 1;
    }
    compact_idx_lowbits = 0;
  }
}
#endif

// Harley-Seal algorithm only works for bitarrays, not quaterarrays, so don't
// add an AVX2 specialization here.
// ...unless something like the interleaved_vec strategy is used?  hmm.  should
// test this on basic frequency counter.
/*
void Count2FreqVec3(const VecW* geno_vvec, uint32_t vec_ct, uint32_t* __restrict alt1_plus_bothset_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 3));
  // Increments bothset_ct by the number of 0b11 in the current block, and
  // alt1_ct by twice the number of 0b10 plus the number of 0b01.
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* geno_vvec_iter = geno_vvec;
  uint32_t alt1_plus_bothset_ct = 0;
  uint32_t bothset_ct = 0;

  while (1) {
    UniVec acc_alt1_plus_bothset;
    UniVec acc_bothset;
    acc_alt1_plus_bothset.vw = vecw_setzero();
    acc_bothset.vw = vecw_setzero();
    const VecW* geno_vvec_stop;
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
      VecW cur_geno_vword1 = *geno_vvec_iter++;
      // process first two vwords simultaneously to minimize linear dependence
      VecW cur_geno_vword2 = *geno_vvec_iter++;
      VecW cur_geno_vword_low_lshifted1 = vecw_slli(cur_geno_vword1 & m1, 1);
      VecW cur_geno_vword_low_lshifted2 = vecw_slli(cur_geno_vword2 & m1, 1);

      // 00 -> 00; 01 -> 01; 10 -> 10; 11 -> 01
      // note that _mm_andnot_si128 flips the *first* argument before the AND
      // operation.
      VecW alt1_plus_bothset1 = (~cur_geno_vword_low_lshifted1) & cur_geno_vword1;
      VecW alt1_plus_bothset2 = (~cur_geno_vword_low_lshifted2) & cur_geno_vword2;

      VecW bothset1 = vecw_srli(cur_geno_vword_low_lshifted1 & cur_geno_vword1, 1);
      VecW bothset2 = vecw_srli(cur_geno_vword_low_lshifted2 & cur_geno_vword2, 1);

      cur_geno_vword1 = *geno_vvec_iter++;
      alt1_plus_bothset1 = (alt1_plus_bothset1 & m2) + (vecw_srli(alt1_plus_bothset1, 2) & m2);
      bothset2 = bothset1 + bothset2;
      alt1_plus_bothset2 = (alt1_plus_bothset2 & m2) + (vecw_srli(alt1_plus_bothset2, 2) & m2);
      cur_geno_vword_low_lshifted1 = vecw_slli(cur_geno_vword1 & m1, 1);

      alt1_plus_bothset2 = alt1_plus_bothset1 + alt1_plus_bothset2;
      // alt1_plus_bothset2 now contains 4-bit values from 0-8, while bothset2
      // contains 2-bit values from 0-2
      // (todo: check whether this is faster if we use double_bothsetx
      // variables instead of bothset1/bothset2)
      bothset1 = vecw_srli(cur_geno_vword_low_lshifted1 & cur_geno_vword1, 1);
      alt1_plus_bothset1 = (~cur_geno_vword_low_lshifted1) & cur_geno_vword1;
      bothset2 = bothset1 + bothset2;
      alt1_plus_bothset1 = (alt1_plus_bothset1 & m2) + (vecw_srli(alt1_plus_bothset1, 2) & m2);

      bothset2 = (bothset2 & m2) + (vecw_srli(bothset2, 2) & m2);
      alt1_plus_bothset2 = alt1_plus_bothset1 + alt1_plus_bothset2;
      // alt1_plus_bothset2 now contains 4-bit values from 0-12, while bothset2
      // contains 4-bit values from 0-6.  aggregate both into 8-bit values.
      bothset2 = (bothset2 & m4) + (vecw_srli(bothset2, 4) & m4);
      alt1_plus_bothset2 = (alt1_plus_bothset2 & m4) + (vecw_srli(alt1_plus_bothset2, 4) & m4);

      acc_bothset.vw = acc_bothset.vw + bothset2;
      acc_alt1_plus_bothset.vw = acc_alt1_plus_bothset.vw + alt1_plus_bothset2;
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m8 = VCONST_W(kMask00FF);
    acc_bothset.vw = (acc_bothset.vw + vecw_srli(acc_bothset.vw, 8)) & m8;
    acc_alt1_plus_bothset.vw = (acc_alt1_plus_bothset.vw & m8) + (vecw_srli(acc_alt1_plus_bothset.vw, 8) & m8);
    bothset_ct += UniVecHsum16(acc_bothset);
    alt1_plus_bothset_ct += UniVecHsum16(acc_alt1_plus_bothset);
  }
}
*/

// todo: benchmark against general-purpose counter
void Count12Vec6(const VecW* geno_vvec, uint32_t vec_ct, uint32_t* __restrict raw_01_ctp, uint32_t* __restrict raw_both_ctp) {
  assert(!(vec_ct % 6));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* geno_vvec_iter = geno_vvec;
  VecW acc_01 = vecw_setzero();
  VecW acc_both = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        *raw_01_ctp = HsumW(acc_01);
        *raw_both_ctp = HsumW(acc_both);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_01 = vecw_setzero();
    VecW inner_acc_both = vecw_setzero();
    const VecW* geno_vvec_stop = &(geno_vvec_iter[cur_incr]);
    do {
      VecW cur_vvec = *geno_vvec_iter++;
      VecW vvec_rshift = vecw_srli(cur_vvec, 1);
      VecW quater_both = m1 & (cur_vvec ^ vvec_rshift);
      VecW quater_01 = quater_both & cur_vvec;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      VecW vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      quater_01 = quater_01 + (vvec_both & cur_vvec);
      quater_both = quater_both + vvec_both;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      quater_01 = quater_01 + (vvec_both & cur_vvec);
      quater_both = quater_both + vvec_both;

      VecW nibble_01 = (quater_01 & m2) + (vecw_srli(quater_01, 2) & m2);
      VecW nibble_both = (quater_both & m2) + (vecw_srli(quater_both, 2) & m2);

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      quater_both = m1 & (cur_vvec ^ vvec_rshift);
      quater_01 = quater_both & cur_vvec;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      quater_01 = quater_01 + (vvec_both & cur_vvec);
      quater_both = quater_both + vvec_both;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      quater_01 = quater_01 + (vvec_both & cur_vvec);
      quater_both = quater_both + vvec_both;

      nibble_01 = nibble_01 + (quater_01 & m2) + (vecw_srli(quater_01, 2) & m2);
      nibble_both = nibble_both + (quater_both & m2) + (vecw_srli(quater_both, 2) & m2);

      inner_acc_01 = inner_acc_01 + (nibble_01 & m4) + (vecw_srli(nibble_01, 4) & m4);
      inner_acc_both = inner_acc_both + (nibble_both & m4) + (vecw_srli(nibble_both, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_01 = acc_01 + vecw_bytesum(inner_acc_01, m0);
    acc_both = acc_both + vecw_bytesum(inner_acc_both, m0);
  }
}

void GenovecCount12Unsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* __restrict raw_01_ctp, uint32_t* __restrict raw_10_ctp) {
  // assumes trailing bits of last genovec word are zeroed out.
  // sample_ct == 0 ok.
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  uint32_t raw_01_ct;
  uint32_t raw_both_ct;
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (6 * kWordsPerVec));
  assert(VecIsAligned(genovec));
  Count12Vec6(R_CAST(const VecW*, genovec), word_idx / kWordsPerVec, &raw_01_ct, &raw_both_ct);
  for (; word_idx != sample_ctl2; ++word_idx) {
    const uintptr_t cur_geno_word = genovec[word_idx];
    const uintptr_t cur_rshift = cur_geno_word >> 1;
    const uintptr_t cur_both = kMask5555 & (cur_geno_word ^ cur_rshift);
    raw_01_ct += Popcount01Word(cur_both & cur_geno_word);
    raw_both_ct += Popcount01Word(cur_both);
  }
  *raw_01_ctp = raw_01_ct;
  *raw_10_ctp = raw_both_ct - raw_01_ct;
}

void Count3FreqVec6(const VecW* geno_vvec, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 6));
  // Sets even_ct to the number of set low bits in the current block, odd_ct to
  // the number of set high bits, and bothset_ct by the number of 0b11s.
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* geno_vvec_iter = geno_vvec;
  VecW acc_even = vecw_setzero();
  VecW acc_odd = vecw_setzero();
  VecW acc_bothset = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        *even_ctp = HsumW(acc_even);
        *odd_ctp = HsumW(acc_odd);
        *bothset_ctp = HsumW(acc_bothset);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_even = vecw_setzero();
    VecW inner_acc_odd = vecw_setzero();
    VecW inner_acc_bothset = vecw_setzero();
    const VecW* geno_vvec_stop = &(geno_vvec_iter[cur_incr]);
    do {
      // hmm, this seems to have more linear dependence than I'd want, but the
      // reorderings I tried just made the code harder to read without helping,
      // so I'll leave this alone
      VecW cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW odd1 = m1 & vecw_srli(cur_geno_vword, 1);
      VecW even1 = m1 & cur_geno_vword;
      VecW bothset1 = odd1 & cur_geno_vword;

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (m1 & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high;
      bothset1 = bothset1 + (cur_geno_vword_high & cur_geno_vword);

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (m1 & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high;
      bothset1 = bothset1 + (cur_geno_vword_high & cur_geno_vword);

      even1 = (even1 & m2) + (vecw_srli(even1, 2) & m2);
      odd1 = (odd1 & m2) + (vecw_srli(odd1, 2) & m2);
      bothset1 = (bothset1 & m2) + (vecw_srli(bothset1, 2) & m2);

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW odd2 = m1 & vecw_srli(cur_geno_vword, 1);
      VecW even2 = m1 & cur_geno_vword;
      VecW bothset2 = odd2 & cur_geno_vword;

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (m1 & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high;
      bothset2 = bothset2 + (cur_geno_vword_high & cur_geno_vword);

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (m1 & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high;
      bothset2 = bothset2 + (cur_geno_vword_high & cur_geno_vword);

      even1 = even1 + (even2 & m2) + (vecw_srli(even2, 2) & m2);
      odd1 = odd1 + (odd2 & m2) + (vecw_srli(odd2, 2) & m2);
      bothset1 = bothset1 + (bothset2 & m2) + (vecw_srli(bothset2, 2) & m2);
      // these now contain 4-bit values from 0-12

      inner_acc_even = inner_acc_even + (even1 & m4) + (vecw_srli(even1, 4) & m4);
      inner_acc_odd = inner_acc_odd + (odd1 & m4) + (vecw_srli(odd1, 4) & m4);
      inner_acc_bothset = inner_acc_bothset + (bothset1 & m4) + (vecw_srli(bothset1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_even = acc_even + vecw_bytesum(inner_acc_even, m0);
    acc_odd = acc_odd + vecw_bytesum(inner_acc_odd, m0);
    acc_bothset = acc_bothset + vecw_bytesum(inner_acc_bothset, m0);
  }
}

void FillInterleavedMaskVec(const uintptr_t* __restrict subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec) {
#ifdef __LP64__
  // This is a cousin of github.com/KWillets/simd_interleave , which was
  // written in response to
  //   https://lemire.me/blog/2018/01/09/how-fast-can-you-bit-interleave-32-bit-integers-simd-edition/
  // AVX2 implementation takes ~40% less time than before, and SSE4.2 takes
  // ~65% less.  This also avoids the Ryzen-screwing _pdep_u64()/_pext_u64()
  // operations.
  const VecW m4 = VCONST_W(kMask0F0F);
#  ifdef USE_SSE42
  const VecW lookup0 = vecw_setr8(
      0, 1, 4, 5, 16, 17, 20, 21,
      64, 65, 68, 69, 80, 81, 84, 85);
  const VecW lookup1 = vecw_slli(lookup0, 1);
#  else
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW* subset_mask_valias = R_CAST(const VecW*, subset_mask);
  VecW* interleaved_mask_valias = R_CAST(VecW*, interleaved_mask_vec);

  for (uint32_t vidx = 0; vidx != base_vec_ct; ++vidx) {
    // I'll assume the compiler can handle this register allocation job.
    VecW cur_vec = subset_mask_valias[vidx];

    // This is the critical lane-crossing operation.  May wrap this (in a way
    // that reduces to a null-op in the SSE4.2 case?), but let's first see how
    // often we use it, and in what ways.
    //   0xd8: {0, 2, 1, 3}
    cur_vec = vecw_permute0xd8_if_avx2(cur_vec);
    const VecW vec_even = cur_vec & m4;
    const VecW vec_odd = vecw_srli(cur_vec, 4) & m4;
    VecW vec_lo = vecw_unpacklo8(vec_even, vec_odd);
    VecW vec_hi = vecw_unpackhi8(vec_even, vec_odd);
#  ifdef USE_SSE42
    vec_lo = vecw_shuffle8(lookup0, vec_lo);
    vec_hi = vecw_shuffle8(lookup1, vec_hi);
#  else
    vec_lo = (vec_lo | vecw_slli(vec_lo, 2)) & m2;
    vec_hi = (vec_hi | vecw_slli(vec_hi, 2)) & m2;
    vec_lo = (vec_lo | vecw_slli(vec_lo, 1)) & m1;
    vec_hi = (vec_hi | vecw_slli(vec_hi, 1)) & m1;
    vec_hi = vecw_slli(vec_hi, 1);
#  endif
    interleaved_mask_valias[vidx] = vec_lo | vec_hi;
  }
#else  // !__LP64__
  for (uint32_t widx = 0; widx != base_vec_ct; ++widx) {
    const uintptr_t orig_word = subset_mask[widx];
    uintptr_t ww_even = S_CAST(uint16_t, orig_word);
    uintptr_t ww_odd = orig_word >> 16;
    ww_even = UnpackHalfwordToWord(ww_even);
    ww_odd = UnpackHalfwordToWord(ww_odd);
    interleaved_mask_vec[widx] = ww_even | (ww_odd << 1);
  }
#endif
}

/*
void GenovecAlleleCtsUnsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* __restrict allele_cts, uint32_t* __restrict bothset_ctp) {
  // assumes trailing bits of last genovec word are zeroed out.
  // sets allele_cts[0] to the number of observed ref alleles, and
  // allele_cts[1] to the number of observed alt1s.
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (3 * kWordsPerVec));
  uint32_t alt1_plus_bothset_ct;
  uint32_t bothset_ct;
  assert(VecIsAligned(genovec));
  Count2FreqVec3(R_CAST(const VecW*, genovec), word_idx / kWordsPerVec, &alt1_plus_bothset_ct, &bothset_ct);
  for (; word_idx != sample_ctl2; ++word_idx) {
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
*/

void GenovecCountFreqsUnsafe(const uintptr_t* genovec, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // fills genocounts[0] with the number of 00s, genocounts[1] with the number
  // of 01s, etc.
  // assumes trailing bits of last genovec word are zeroed out.
  // sample_ct == 0 ok.
  // no longer requires vector-alignment.
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  uint32_t even_ct;
  uint32_t odd_ct;
  uint32_t bothset_ct;
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (6 * kWordsPerVec));
  Count3FreqVec6(R_CAST(const VecW*, genovec), word_idx / kWordsPerVec, &even_ct, &odd_ct, &bothset_ct);
  for (; word_idx != sample_ctl2; ++word_idx) {
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

// geno_vvec now allowed to be unaligned.
void CountSubset3FreqVec6(const VecW* __restrict geno_vvec, const VecW* __restrict interleaved_mask_vvec, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 6));
  // Sets even_ct to the number of set low bits in the current block after
  // subsetting, odd_ct to the number of set high bits, and bothset_ct by the
  // number of 0b11s.
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* geno_vvec_iter = geno_vvec;
  const VecW* interleaved_mask_vvec_iter = interleaved_mask_vvec;
  VecW acc_even = vecw_setzero();
  VecW acc_odd = vecw_setzero();
  VecW acc_bothset = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        *even_ctp = HsumW(acc_even);
        *odd_ctp = HsumW(acc_odd);
        *bothset_ctp = HsumW(acc_bothset);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_even = vecw_setzero();
    VecW inner_acc_odd = vecw_setzero();
    VecW inner_acc_bothset = vecw_setzero();
    const VecW* geno_vvec_stop = &(geno_vvec_iter[cur_incr]);
    do {
      VecW interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      VecW cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW cur_mask = interleaved_mask_vword & m1;
      VecW odd1 = cur_mask & vecw_srli(cur_geno_vword, 1);
      VecW even1 = cur_mask & cur_geno_vword;
      VecW bothset1 = odd1 & cur_geno_vword;

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (cur_mask & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high_masked;
      bothset1 = bothset1 + (cur_geno_vword_high_masked & cur_geno_vword);

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (cur_mask & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high_masked;
      bothset1 = bothset1 + (cur_geno_vword_high_masked & cur_geno_vword);

      even1 = (even1 & m2) + (vecw_srli(even1, 2) & m2);
      odd1 = (odd1 & m2) + (vecw_srli(odd1, 2) & m2);
      bothset1 = (bothset1 & m2) + (vecw_srli(bothset1, 2) & m2);

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW odd2 = cur_mask & vecw_srli(cur_geno_vword, 1);
      VecW even2 = cur_mask & cur_geno_vword;
      VecW bothset2 = odd2 & cur_geno_vword;

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (cur_mask & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high_masked;
      bothset2 = bothset2 + (cur_geno_vword_high_masked & cur_geno_vword);

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (cur_mask & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high_masked;
      bothset2 = bothset2 + (cur_geno_vword_high_masked & cur_geno_vword);

      even1 = even1 + (even2 & m2) + (vecw_srli(even2, 2) & m2);
      odd1 = odd1 + (odd2 & m2) + (vecw_srli(odd2, 2) & m2);
      bothset1 = bothset1 + (bothset2 & m2) + (vecw_srli(bothset2, 2) & m2);
      // these now contain 4-bit values from 0-12

      inner_acc_even = inner_acc_even + (even1 & m4) + (vecw_srli(even1, 4) & m4);
      inner_acc_odd = inner_acc_odd + (odd1 & m4) + (vecw_srli(odd1, 4) & m4);
      inner_acc_bothset = inner_acc_bothset + (bothset1 & m4) + (vecw_srli(bothset1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_even = acc_even + vecw_bytesum(inner_acc_even, m0);
    acc_odd = acc_odd + vecw_bytesum(inner_acc_odd, m0);
    acc_bothset = acc_bothset + vecw_bytesum(inner_acc_bothset, m0);
  }
}

void GenovecCountSubsetFreqs(const uintptr_t* __restrict genovec, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // fills genocounts[0] with the number of 00s, genocounts[1] with the number
  // of 01s, etc.
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctv2 = QuaterCtToVecCt(raw_sample_ct);
  uint32_t even_ct;
  uint32_t odd_ct;
  uint32_t bothset_ct;
#ifdef __LP64__
  uint32_t vec_idx = raw_sample_ctv2 - (raw_sample_ctv2 % 6);
  CountSubset3FreqVec6(R_CAST(const VecW*, genovec), R_CAST(const VecW*, sample_include_interleaved_vec), vec_idx, &even_ct, &odd_ct, &bothset_ct);
  const uintptr_t* genovec_iter = &(genovec[kWordsPerVec * vec_idx]);
  const uintptr_t* interleaved_mask_iter = &(sample_include_interleaved_vec[(kWordsPerVec / 2) * vec_idx]);
#  ifdef USE_AVX2
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  uintptr_t mask_base3 = 0;
  uintptr_t mask_base4 = 0;
  for (; vec_idx != raw_sample_ctv2; ++vec_idx) {
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
    for (uint32_t vechalf_idx = 0; ; ++vechalf_idx) {
      const uintptr_t cur_geno_word1 = *genovec_iter++;
      const uintptr_t cur_geno_word2 = *genovec_iter++;
      const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
      const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
      even_ct += PopcountWord(((cur_geno_word1 & mask_word1) << 1) | (cur_geno_word2 & mask_word2));
      odd_ct += PopcountWord((cur_geno_word1_high_masked << 1) | cur_geno_word2_high_masked);
      bothset_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
      if (vechalf_idx) {
        break;
      }
      mask_word1 = mask_word3;
      mask_word2 = mask_word4;
    }
  }
#  else  // not USE_AVX2
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  for (; vec_idx != raw_sample_ctv2; ++vec_idx) {
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
  CountSubset3FreqVec6(R_CAST(const VecW*, genovec), R_CAST(const VecW*, sample_include_interleaved_vec), word_idx, &even_ct, &odd_ct, &bothset_ct);
  const uintptr_t* interleaved_mask_iter = &(sample_include_interleaved_vec[word_idx / 2]);
  uintptr_t mask_base = 0;
  for (; word_idx != raw_sample_ctv2; ++word_idx) {
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

#ifdef __arm__
#  error "Unaligned accesses in SmallGenoarrCount3FreqIncr()."
#endif
void SmallGenoarrCount3FreqIncr(const uintptr_t* genoarr_iter, uint32_t byte_ct, uint32_t* even_ctp, uint32_t* odd_ctp, uint32_t* bothset_ctp) {
  for (uint32_t bytes_left = byte_ct; ; ) {
    uintptr_t cur_geno_word;
    if (bytes_left < kBytesPerWord) {
      if (!bytes_left) {
        return;
      }
      cur_geno_word = ProperSubwordLoad(genoarr_iter, bytes_left);
      bytes_left = 0;
    } else {
      cur_geno_word = *genoarr_iter++;
      bytes_left -= kBytesPerWord;
    }
    const uintptr_t cur_geno_word_high = kMask5555 & (cur_geno_word >> 1);
    *even_ctp += Popcount01Word(cur_geno_word & kMask5555);
    *odd_ctp += Popcount01Word(cur_geno_word_high);
    *bothset_ctp += Popcount01Word(cur_geno_word & cur_geno_word_high);
  }
}

void GenoarrCountFreqs(const unsigned char* genoarr, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
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
    // Probably want to search out all instances of __arm__ and make the code
    // standard-compliant, if that can be done without a speed penalty.  Though
    // it makes sense to wait until more is known about Apple's MacBook
    // processor plans...
    SmallGenoarrCount3FreqIncr(R_CAST(const uintptr_t*, genoarr), lead_byte_ct, &even_ct, &odd_ct, &bothset_ct);
    genoarr_iter = R_CAST(const uintptr_t*, &(genoarr[lead_byte_ct]));
    const uint32_t remaining_full_vec_ct = remaining_sample_ct / kQuatersPerVec;
    uint32_t even_ct_incr;
    uint32_t odd_ct_incr;
    uint32_t bothset_ct_incr;
    const uint32_t vec_ct = remaining_full_vec_ct - (remaining_full_vec_ct % 6);
    Count3FreqVec6(R_CAST(const VecW*, genoarr_iter), vec_ct, &even_ct_incr, &odd_ct_incr, &bothset_ct_incr);
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
void GenoarrCountSubsetFreqs(const unsigned char* genoarr, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // does not read past the end of genoarr
  const uint32_t raw_sample_ctv2 = QuaterCtToVecCt(raw_sample_ct);
  uint32_t even_ct;
  uint32_t odd_ct;
  uint32_t bothset_ct;
  uint32_t vec_idx = raw_sample_ctv2 - (raw_sample_ctv2 % 6);
  CountSubset3FreqVec6(R_CAST(const VecW*, genoarr), R_CAST(const VecW*, sample_include_interleaved_vec), vec_idx, &even_ct, &odd_ct, &bothset_ct);
  const uintptr_t* genoarr_iter = &(R_CAST(const uintptr_t*, genoarr)[kWordsPerVec * vec_idx]);
#ifdef __LP64__
  const uintptr_t* interleaved_mask_iter = &(sample_include_interleaved_vec[vec_idx * (kWordsPerVec / 2)]);
#else
  // bugfix (19 Jul 2018): (kWordsPerVec / 2) doesn't work in 32-bit case
  const uintptr_t* interleaved_mask_iter = &(sample_include_interleaved_vec[(vec_idx * kWordsPerVec) / 2]);
#endif
#ifdef USE_AVX2
  const uint32_t halfvec_idx_trail = (raw_sample_ct + 3) / (kBitsPerVec / 4);
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  uintptr_t mask_base3 = 0;
  uintptr_t mask_base4 = 0;
  for (; vec_idx != raw_sample_ctv2; ++vec_idx) {
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
    uint32_t vechalf_idx = 0;
    while (1) {
      uintptr_t cur_geno_word1;
      uintptr_t cur_geno_word2;
      if (2 * vec_idx + vechalf_idx < halfvec_idx_trail) {
        cur_geno_word1 = *genoarr_iter++;
        cur_geno_word2 = *genoarr_iter++;
      } else {
        const uint32_t remaining_byte_ct = QuaterCtToByteCt(raw_sample_ct) % kBytesPerVec;
        // todo: check if this harms usual-case loop efficiency
        vechalf_idx = 1;
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
      if (vechalf_idx) {
        break;
      }
      ++vechalf_idx;
      mask_word1 = mask_word3;
      mask_word2 = mask_word4;
    }
  }
#else  // not USE_AVX2
  const uint32_t vec_idx_trail = (raw_sample_ct + 3) / kQuatersPerVec;
#  ifdef __LP64__
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  for (; vec_idx != raw_sample_ctv2; ++vec_idx) {
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
  for (; vec_idx != raw_sample_ctv2; ++vec_idx) {
    uintptr_t mask_word;
    if (!(vec_idx % 2)) {
      mask_base = *interleaved_mask_iter++;
      mask_word = mask_base & kMask5555;
    } else {
      mask_word = (mask_base >> 1) & kMask5555;
    }
    uintptr_t cur_geno_word;
    if (vec_idx < vec_idx_trail) {
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

void GenovecCountSubsetFreqs2(const uintptr_t* __restrict genovec, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // slower GenovecCountSubsetFreqs() which does not require
  // sample_include_interleaved_vec to be precomputed.
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctl2 = QuaterCtToWordCt(raw_sample_ct);
  const uint32_t fullword_ct = raw_sample_ctl2 / 2;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  for (uint32_t widx = 0; widx != fullword_ct; ++widx) {
    // possible todo: try vectorizing this loop in SSE4.2+ high-sample-ct case
    // with shuffle-based dynamic unpacking of sample_include?
    const uintptr_t mask_word = sample_include[widx];
    if (mask_word) {
      uintptr_t geno_word = genovec[2 * widx];
      uintptr_t geno_even = PackWordToHalfwordMask5555(geno_word);
      uintptr_t geno_odd = PackWordToHalfwordMaskAAAA(geno_word);
      geno_word = genovec[2 * widx + 1];
      geno_even |= S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_word)) << kBitsPerWordD2;
      geno_odd |= S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(geno_word)) << kBitsPerWordD2;
      const uintptr_t geno_even_masked = geno_even & mask_word;
      even_ct += PopcountWord(geno_even_masked);
      odd_ct += PopcountWord(geno_odd & mask_word);
      bothset_ct += PopcountWord(geno_odd & geno_even_masked);
    }
  }
  if (raw_sample_ctl2 % 2) {
    const uintptr_t mask_hw = sample_include[fullword_ct];
    if (mask_hw) {
      const uintptr_t geno_word = genovec[2 * fullword_ct];
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

void GenoarrCountSubsetIntersectFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict subset1, const uintptr_t* __restrict subset2, uint32_t raw_sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctl2 = QuaterCtToWordCt(raw_sample_ct);
  const uint32_t fullword_ct = raw_sample_ctl2 / 2;
  uint32_t subset_intersect_ct = 0;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  for (uint32_t widx = 0; widx != fullword_ct; ++widx) {
    // hmm, there may be little point to vectorizing this
    const uintptr_t mask_word = subset1[widx] & subset2[widx];
    if (mask_word) {
      uintptr_t geno_word = genoarr[2 * widx];
      uintptr_t geno_even = PackWordToHalfwordMask5555(geno_word);
      uintptr_t geno_odd = PackWordToHalfwordMaskAAAA(geno_word);
      geno_word = genoarr[2 * widx + 1];
      geno_even |= S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_word)) << kBitsPerWordD2;
      geno_odd |= S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(geno_word)) << kBitsPerWordD2;
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

void GenovecCountFreqs(const uintptr_t* genovec, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // this masks out trailing genovec bits
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
  assert(VecIsAligned(genovec));
  const VecW not_m1 = VCONST_W(kMaskAAAA);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    VecW cur_vec = vptr[vidx];
    // flip high bit iff low bit is unset
    vptr[vidx] = cur_vec ^ ((~vecw_slli(cur_vec, 1)) & not_m1);
  }
}

void GenovecInvertCopyUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict genovec_inverted_copy) {
  // flips 0 to 2 and vice versa.
  // "unsafe" because trailing bits are not zeroed out.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(VecIsAligned(genovec));
  const VecW not_m1 = VCONST_W(kMaskAAAA);
  const VecW* vin_ptr = R_CAST(const VecW*, genovec);
  VecW* vout_ptr = R_CAST(VecW*, genovec_inverted_copy);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    const VecW cur_vec = vin_ptr[vidx];
    // flip high bit iff low bit is unset
    vout_ptr[vidx] = cur_vec ^ ((~vecw_slli(cur_vec, 1)) & not_m1);
  }
}

void GenovecNonmissingToZeroUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // sets 1 and 2 to zero; leaves 3s untouched.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(VecIsAligned(genovec));
  const VecW m1 = VCONST_W(kMask5555);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    VecW cur_vec = vptr[vidx];
    const VecW cur_vec_rshifted = vecw_srli(cur_vec, 1);
    cur_vec = cur_vec & m1;
    cur_vec = cur_vec & cur_vec_rshifted;
    vptr[vidx] = cur_vec | vecw_slli(cur_vec, 1);
  }
}

void GenovecNonzeroToMissingUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // converts 1s and 2s to 3s, leaves zeroes untouched.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(VecIsAligned(genovec));
  const VecW m1 = VCONST_W(kMask5555);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    VecW cur_vec = vptr[vidx];
    const VecW cur_vec_rshifted = vecw_srli(cur_vec, 1);
    cur_vec = cur_vec | cur_vec_rshifted;
    cur_vec = cur_vec & m1;
    vptr[vidx] = cur_vec | vecw_slli(cur_vec, 1);
  }
}

void GenovecNontwoToMissingUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // 0 -> 3, 1 -> 3.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(VecIsAligned(genovec));
  const VecW not_m1 = VCONST_W(kMaskAAAA);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    const VecW cur_vec = vptr[vidx];
    const VecW cur_vec_hi = not_m1 & (~cur_vec);
    const VecW cur_or = cur_vec_hi | vecw_srli(cur_vec_hi, 1);
    vptr[vidx] = cur_vec | cur_or;
  }
}

void GenovecNonzeroToMissingThenInvertUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // 0 -> 2, 1 -> 3, 2 -> 3
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(VecIsAligned(genovec));
  const VecW not_m1 = VCONST_W(kMaskAAAA);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    const VecW cur_vec = vptr[vidx];
    vptr[vidx] = cur_vec | vecw_srli(cur_vec, 1) | not_m1;
  }
}

void GenovecInvertThenNonzeroToMissingUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // 0 -> 3, 1 -> 3, 2 -> 0
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  assert(VecIsAligned(genovec));
  const VecW m1 = VCONST_W(kMask5555);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    const VecW cur_vec = vptr[vidx];
    const VecW cur_vec_rshifted = vecw_srli(cur_vec, 1);
    const VecW not2 = m1 & (~(cur_vec_rshifted & (~cur_vec)));
    vptr[vidx] = not2 | vecw_slli(not2, 1);
  }
}

void DifflistCountSubsetFreqs(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t common_geno, uint32_t difflist_len, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  STD_ARRAY_REF_FILL0(4, genocounts);
  for (uint32_t difflist_idx = 0; difflist_idx != difflist_len; ++difflist_idx) {
    const uint32_t raw_sample_idx = difflist_sample_ids[difflist_idx];
    if (IsSet(sample_include, raw_sample_idx)) {
      genocounts[GetQuaterarrEntry(raregeno, difflist_idx)] += 1;
    }
  }
  genocounts[common_geno] = sample_ct - genocounts[0] - genocounts[1] - genocounts[2] - genocounts[3];
}


static_assert(kPglQuaterTransposeBatch == S_CAST(uint32_t, kQuatersPerCacheline), "TransposeQuaterblock64() needs to be updated.");
#ifdef __LP64__
void TransposeQuaterblock64(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1) {
  // buf0 and buf1 must each be vector-aligned and have size 16k
  // Tried using previous AVX2 small-buffer approach, but that was a bit
  // worse... though maybe it should be revisited?
  const uint32_t buf0_row_ct = DivUp(write_batch_size, 32);
  {
    // Fold the first 3 shuffles into the ingestion loop.
    // Can fold 2 additional shuffles here by ingesting uint16s instead.  That
    // removes the need for the second loop, halving the workspace requirement.
    // Benchmark results of that approach are slightly but consistently worse,
    // though.
    const uintptr_t* initial_read_iter = R_CAST(const uintptr_t*, read_iter);
    const uintptr_t* initial_read_end = &(initial_read_iter[buf0_row_ct]);
    uintptr_t* initial_target_iter = R_CAST(uintptr_t*, buf0);
    const uint32_t read_batch_rem = kQuatersPerCacheline - read_batch_size;
    // Tried fancier vector-based ingestion, didn't help.
    for (; initial_read_iter != initial_read_end; ++initial_read_iter) {
      const uintptr_t* read_iter_tmp = initial_read_iter;
      for (uint32_t ujj = 0; ujj != read_batch_size; ++ujj) {
        *initial_target_iter++ = *read_iter_tmp;
        read_iter_tmp = &(read_iter_tmp[read_ul_stride]);
      }
      ZeroWArr(read_batch_rem, initial_target_iter);
      initial_target_iter = &(initial_target_iter[read_batch_rem]);
    }
  }

  // shuffle from 64 -> 16
  const uint32_t buf1_row_ct = DivUp(write_batch_size, 8);
  {
    // full buf0 row is 256 * 8 bytes
    // full buf1 row is 512 bytes
    const VecW* buf0_read_iter = R_CAST(const VecW*, buf0);
    __m128i* write_iter0 = R_CAST(__m128i*, buf1);
    const uint32_t buf0_row_clwidth = DivUp(read_batch_size, 8);
#  ifdef USE_SSE42
    const VecW gather_u32s = vecw_setr8(0, 1, 8, 9, 2, 3, 10, 11,
                                        4, 5, 12, 13, 6, 7, 14, 15);
#  else
    const VecW m16 = VCONST_W(kMask0000FFFF);
#  endif
    for (uint32_t bidx = 0; bidx != buf0_row_ct; ++bidx) {
      __m128i* write_iter1 = &(write_iter0[32]);
      __m128i* write_iter2 = &(write_iter1[32]);
      __m128i* write_iter3 = &(write_iter2[32]);
      for (uint32_t clidx = 0; clidx != buf0_row_clwidth; ++clidx) {
#  ifdef USE_AVX2
        VecW loader0 = buf0_read_iter[clidx * 2];
        VecW loader1 = buf0_read_iter[clidx * 2 + 1];
        loader0 = vecw_shuffle8(loader0, gather_u32s);
        loader1 = vecw_shuffle8(loader1, gather_u32s);
        //    (0,0) (0,1) (1,0) ... (7,1) (0,2) ... (7,3) (8,0) ... (31,3)
        //      (0,4) ... (31,7)
        // -> (0,0) ... (7,1) (0,2) ... (7,3) (8,0) ... (15,3) (0,4) ... (15,7)
        //      (16,0) ... (31,3) (16,4) ...
        loader0 = vecw_permute0xd8_if_avx2(loader0);
        loader1 = vecw_permute0xd8_if_avx2(loader1);
        // -> (0,0) ... (7,1) (0,2) ... (7,3) (0,4) ... (7,7) (8,0) ... (15,7)
        //      (16,0) ... (23,7) (24,0) ... (31,7)
        const __m256i vec_lo = _mm256_shuffle_epi32(R_CAST(__m256i, loader0), 0xd8);
        const __m256i vec_hi = _mm256_shuffle_epi32(R_CAST(__m256i, loader1), 0xd8);
        const __m256i final0145 = _mm256_unpacklo_epi64(vec_lo, vec_hi);
        const __m256i final2367 = _mm256_unpackhi_epi64(vec_lo, vec_hi);
        // GCC doesn't support _mm256_storeu_si128i as of this writing.
        write_iter0[clidx] = _mm256_castsi256_si128(final0145);
        write_iter1[clidx] = _mm256_castsi256_si128(final2367);
        write_iter2[clidx] = _mm256_extracti128_si256(final0145, 1);
        write_iter3[clidx] = _mm256_extracti128_si256(final2367, 1);
#  else
        VecW loader0 = buf0_read_iter[clidx * 4];
        VecW loader1 = buf0_read_iter[clidx * 4 + 1];
        VecW loader2 = buf0_read_iter[clidx * 4 + 2];
        VecW loader3 = buf0_read_iter[clidx * 4 + 3];
#    ifdef USE_SSE42
        loader0 = vecw_shuffle8(loader0, gather_u32s);
        loader1 = vecw_shuffle8(loader1, gather_u32s);
        loader2 = vecw_shuffle8(loader2, gather_u32s);
        loader3 = vecw_shuffle8(loader3, gather_u32s);
#    else
        VecW tmp_lo = vecw_unpacklo16(loader0, loader1);
        VecW tmp_hi = vecw_unpackhi16(loader0, loader1);
        loader0 = (tmp_lo & m16) | vecw_and_notfirst(m16, vecw_slli(tmp_hi, 16));
        loader1 = (vecw_srli(tmp_lo, 16) & m16) | (vecw_and_notfirst(m16, tmp_hi));
        tmp_lo = vecw_unpacklo16(loader2, loader3);
        tmp_hi = vecw_unpackhi16(loader2, loader3);
        loader2 = (tmp_lo & m16) | vecw_and_notfirst(m16, vecw_slli(tmp_hi, 16));
        loader3 = (vecw_srli(tmp_lo, 16) & m16) | (vecw_and_notfirst(m16, tmp_hi));
#    endif
        //    (0,0) (0,1) (1,0) ... (7,1) (0,2) ... (7,3) (8,0) ... (31,3)
        //  + (0,4) ... (31,7)
        // -> (0,0) ... (7,3) (0,4) ... (7,7) (8,0) ... (15,7)
        const VecW lo_0_15 = vecw_unpacklo32(loader0, loader1);
        const VecW lo_16_31 = vecw_unpackhi32(loader0, loader1);
        const VecW hi_0_15 = vecw_unpacklo32(loader2, loader3);
        const VecW hi_16_31 = vecw_unpackhi32(loader2, loader3);
        write_iter0[clidx] = R_CAST(__m128i, vecw_unpacklo64(lo_0_15, hi_0_15));
        write_iter1[clidx] = R_CAST(__m128i, vecw_unpackhi64(lo_0_15, hi_0_15));
        write_iter2[clidx] = R_CAST(__m128i, vecw_unpacklo64(lo_16_31, hi_16_31));
        write_iter3[clidx] = R_CAST(__m128i, vecw_unpackhi64(lo_16_31, hi_16_31));
#  endif
      }
      buf0_read_iter = &(buf0_read_iter[2048 / kBytesPerVec]);
      write_iter0 = &(write_iter3[32]);
    }
  }

  // movemask from 16 -> 2
  const VecW* source_iter = R_CAST(VecW*, buf1);
  const VecW m8 = VCONST_W(kMask00FF);

  // Take advantage of current function contract.
  const uint32_t buf1_fullrow_ct = (write_batch_size + 3) / 8;

  const uint32_t write_v8ui_stride = kVec8thUintPerWord * write_ul_stride;
  const uint32_t vec_ct = DivUp(read_batch_size, (kBytesPerVec / 2));
  Vec8thUint* target_iter0 = R_CAST(Vec8thUint*, write_iter);
  for (uint32_t uii = 0; uii != buf1_fullrow_ct; ++uii) {
    Vec8thUint* target_iter1 = &(target_iter0[write_v8ui_stride]);
    Vec8thUint* target_iter2 = &(target_iter1[write_v8ui_stride]);
    Vec8thUint* target_iter3 = &(target_iter2[write_v8ui_stride]);
    Vec8thUint* target_iter4 = &(target_iter3[write_v8ui_stride]);
    Vec8thUint* target_iter5 = &(target_iter4[write_v8ui_stride]);
    Vec8thUint* target_iter6 = &(target_iter5[write_v8ui_stride]);
    Vec8thUint* target_iter7 = &(target_iter6[write_v8ui_stride]);
    for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
      const VecW loader = source_iter[vidx];
      // Using goal bit-coordinates, where '_' indicates irrelevant content, we
      // set target_0123 to
      //   _ (0, 0) _ (1, 0) _ (2, 0) _ (3, 0) _ (0, 1) _ (1, 1) _ (2, 1) ...
      // and target_4567 to
      //   _ (4, 0) _ (5, 0) _ (6, 0) _ (7, 0) _ (4, 1) _ (5, 1) _ (6, 1) ...
      // This is perfectly arranged for movemask.
      VecW target_4567 = (m8 & vecw_srli(loader, 7)) | vecw_and_notfirst(m8, loader);
      target_iter7[vidx] = vecw_movemask(target_4567);
      target_4567 = vecw_slli(target_4567, 2);
      target_iter6[vidx] = vecw_movemask(target_4567);
      target_4567 = vecw_slli(target_4567, 2);
      target_iter5[vidx] = vecw_movemask(target_4567);
      target_4567 = vecw_slli(target_4567, 2);
      target_iter4[vidx] = vecw_movemask(target_4567);
      VecW target_0123 = (vecw_slli(loader, 1) & m8) | vecw_and_notfirst(m8, vecw_slli(loader, 8));
      target_iter3[vidx] = vecw_movemask(target_0123);
      target_0123 = vecw_slli(target_0123, 2);
      target_iter2[vidx] = vecw_movemask(target_0123);
      target_0123 = vecw_slli(target_0123, 2);
      target_iter1[vidx] = vecw_movemask(target_0123);
      target_0123 = vecw_slli(target_0123, 2);
      target_iter0[vidx] = vecw_movemask(target_0123);
    }
    source_iter = &(source_iter[(2 * kPglQuaterTransposeBatch) / kBytesPerVec]);
    target_iter0 = &(target_iter7[write_v8ui_stride]);
  }
  if (buf1_fullrow_ct == buf1_row_ct) {
    return;
  }
  Vec8thUint* target_iter1 = &(target_iter0[write_v8ui_stride]);
  Vec8thUint* target_iter2 = &(target_iter1[write_v8ui_stride]);
  Vec8thUint* target_iter3 = &(target_iter2[write_v8ui_stride]);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    const VecW loader = source_iter[vidx];
    VecW target_0123 = (vecw_slli(loader, 1) & m8) | vecw_and_notfirst(m8, vecw_slli(loader, 8));
    target_iter3[vidx] = vecw_movemask(target_0123);
    target_0123 = vecw_slli(target_0123, 2);
    target_iter2[vidx] = vecw_movemask(target_0123);
    target_0123 = vecw_slli(target_0123, 2);
    target_iter1[vidx] = vecw_movemask(target_0123);
    target_0123 = vecw_slli(target_0123, 2);
    target_iter0[vidx] = vecw_movemask(target_0123);
  }
}
#else  // !__LP64__
static_assert(kWordsPerVec == 1, "TransposeQuaterblock32() needs to be updated.");
void TransposeQuaterblock32(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1) {
  // buf0 and buf1 must each be vector-aligned and have size 16k
  // defining them as unsigned char* might prevent a strict-aliasing issue?
  // (might need to go through greater contortions to actually be safe?)
  const uint32_t buf_row_ct = QuaterCtToByteCt(write_batch_size);
  // fold the first 6 shuffles into the initial ingestion loop
  const unsigned char* initial_read_iter = R_CAST(const unsigned char*, read_iter);
  const unsigned char* initial_read_end = &(initial_read_iter[buf_row_ct]);
  unsigned char* initial_target_iter = buf0;
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t read_batch_rem = kQuatersPerCacheline - read_batch_size;
  for (; initial_read_iter != initial_read_end; ++initial_read_iter) {
    const unsigned char* read_iter_tmp = initial_read_iter;
    for (uint32_t ujj = 0; ujj != read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_byte_stride]);
    }
    initial_target_iter = memsetua(initial_target_iter, 0, read_batch_rem);
  }

  // second-to-last shuffle, 8 bit spacing -> 4
  const VecW* source_iter = R_CAST(VecW*, buf0);
  uintptr_t* target_iter0 = R_CAST(uintptr_t*, buf1);
  const uint32_t write_word_ct = QuaterCtToWordCt(read_batch_size);
  const uint32_t penult_inner_loop_iter_ct = 2 * write_word_ct;
  const uint32_t cur_write_skip = 2 * kWordsPerCacheline - penult_inner_loop_iter_ct;
  for (uint32_t uii = 0; uii != buf_row_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 2]);
    for (uint32_t ujj = 0; ujj != penult_inner_loop_iter_ct; ++ujj) {
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

  // last shuffle, 4 bit spacing -> 2
  source_iter = R_CAST(VecW*, buf1);
  target_iter0 = write_iter;
  const uint32_t last_loop_iter_ct = DivUp(write_batch_size, 2);
  for (uint32_t uii = 0; uii != last_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    for (uint32_t ujj = 0; ujj != write_word_ct; ++ujj) {
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
    }
    source_iter = &(source_iter[2 * (kWordsPerCacheline - write_word_ct)]);
    target_iter0 = &(target_iter1[write_ul_stride]);
  }
}
#endif  // !__LP64__

void BiallelicDosage16Invert(uint32_t dosage_ct, uint16_t* dosage_main) {
  // replace each x with (32768 - x).
  // compiler is smart enough to vectorize this.
  for (uint32_t uii = 0; uii != dosage_ct; ++uii) {
    dosage_main[uii] = 32768 - dosage_main[uii];
  }
}

void BiallelicDphase16Invert(uint32_t dphase_ct, int16_t* dphase_delta) {
  for (uint32_t uii = 0; uii != dphase_ct; ++uii) {
    dphase_delta[uii] = -dphase_delta[uii];
  }
}

void GenovecToMissingnessUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict missingness) {
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  Halfword* missingness_alias = R_CAST(Halfword*, missingness);
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uintptr_t cur_geno_word = genovec[widx];
    missingness_alias[widx] = PackWordToHalfwordMask5555(cur_geno_word & (cur_geno_word >> 1));
  }
  if (sample_ctl2 % 2) {
    missingness_alias[sample_ctl2] = 0;
  }
}

void GenovecToNonmissingnessUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict nonmissingness) {
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  Halfword* nonmissingness_alias = R_CAST(Halfword*, nonmissingness);
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uintptr_t cur_geno_word = genovec[widx];
    nonmissingness_alias[widx] = PackWordToHalfwordMask5555(~(cur_geno_word & (cur_geno_word >> 1)));
  }
}

void SplitHomRef2hetUnsafeW(const uintptr_t* genoarr, uint32_t inword_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf) {
  Halfword* hom_alias = R_CAST(Halfword*, hom_buf);
  Halfword* ref2het_alias = R_CAST(Halfword*, ref2het_buf);
  for (uint32_t widx = 0; widx != inword_ct; ++widx) {
    const uintptr_t geno_word = genoarr[widx];
    hom_alias[widx] = PackWordToHalfwordMask5555(~geno_word);
    ref2het_alias[widx] = PackWordToHalfwordMaskAAAA(~geno_word);
  }
}

void SplitHomRef2het(const uintptr_t* genoarr, uint32_t sample_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf) {
  const uint32_t full_outword_ct = sample_ct / kBitsPerWord;
  SplitHomRef2hetUnsafeW(genoarr, full_outword_ct * 2, hom_buf, ref2het_buf);
  const uint32_t remainder = sample_ct % kBitsPerWord;
  if (remainder) {
    uintptr_t geno_word = genoarr[full_outword_ct * 2];
    uintptr_t hom_word = PackWordToHalfwordMask5555(~geno_word);
    uintptr_t ref2het_word = PackWordToHalfwordMaskAAAA(~geno_word);
    if (remainder > kBitsPerWordD2) {
      geno_word = genoarr[full_outword_ct * 2 + 1];
      hom_word |= S_CAST(uintptr_t, PackWordToHalfwordMask5555(~geno_word)) << kBitsPerWordD2;
      ref2het_word |= S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(~geno_word)) << kBitsPerWordD2;
    }
    const uintptr_t cur_mask = (k1LU << remainder) - 1;
    hom_buf[full_outword_ct] = hom_word & cur_mask;
    ref2het_buf[full_outword_ct] = ref2het_word & cur_mask;
  }
}


/*
void GenoarrLookup256x1bx4(const uintptr_t* genoarr, const void* table256x1bx4, uint32_t sample_ct, void* __restrict result) {
  const uint32_t* table_alias = S_CAST(const uint32_t*, table256x1bx4);
  const unsigned char* genoarr_alias = R_CAST(const unsigned char*, genoarr);
  uint32_t* result_alias = S_CAST(uint32_t*, result);
  const uint32_t full_byte_ct = sample_ct / 4;
  for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
    result_alias[byte_idx] = table_alias[genoarr_alias[byte_idx]];
  }
  const uint32_t remainder = sample_ct % 4;
  if (remainder) {
    unsigned char* result_last = R_CAST(unsigned char*, &(result_alias[full_byte_ct]));
    uintptr_t geno_byte = genoarr_alias[full_byte_ct];
    for (uint32_t uii = 0; uii != remainder; ++uii) {
      result_last[uii] = table_alias[geno_byte & 3];
      geno_byte >>= 2;
    }
  }
}
*/

void GenoarrLookup16x4bx2(const uintptr_t* genoarr, const void* table16x4bx2, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table16x4bx2);
  uint64_t* result_iter = S_CAST(uint64_t*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
          memcpy(result_iter, &(table_alias[geno_word & 3]), 4);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
      *result_iter++ = table_alias[cur_2geno];
      geno_word >>= 4;
    }
  }
}

// this might be important for genovec -> AlleleCode expansion
void GenoarrLookup256x2bx4(const uintptr_t* genoarr, const void* table256x2bx4, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table256x2bx4);
  const unsigned char* genoarr_alias = R_CAST(const unsigned char*, genoarr);
  uint64_t* result_alias = S_CAST(uint64_t*, result);
  const uint32_t full_byte_ct = sample_ct / 4;
  for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
    result_alias[byte_idx] = table_alias[genoarr_alias[byte_idx]];
  }
  const uint32_t remainder = sample_ct % 4;
  if (remainder) {
    uint16_t* result_last = R_CAST(uint16_t*, &(result_alias[full_byte_ct]));
    uintptr_t geno_byte = genoarr_alias[full_byte_ct];
    for (uint32_t uii = 0; uii != remainder; ++uii) {
      memcpy_k(&(result_last[uii]), &(table_alias[geno_byte & 3]), 2);
      geno_byte >>= 2;
    }
  }
}

#ifdef __LP64__
void GenoarrLookup4x16b(const uintptr_t* genoarr, const void* table4x16b, uint32_t sample_ct, void* result) {
  const __m128i* table_alias = S_CAST(const __m128i*, table4x16b);
  __m128i* result_iter = S_CAST(__m128i*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD2;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      _mm_storeu_si128(result_iter, table_alias[geno_word & 3]);
      ++result_iter;
      geno_word >>= 2;
    }
  }
}

void GenoarrLookup16x8bx2(const uintptr_t* genoarr, const void* table16x8bx2, uint32_t sample_ct, void* __restrict result) {
  const __m128i* table_alias = S_CAST(const __m128i*, table16x8bx2);
  __m128i* result_iter = S_CAST(__m128i*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
          memcpy(result_iter, &(table_alias[geno_word & 3]), 8);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
      _mm_storeu_si128(result_iter, table_alias[cur_2geno]);
      ++result_iter;
      geno_word >>= 4;
    }
  }
}

void GenoarrLookup256x4bx4(const uintptr_t* genoarr, const void* table256x4bx4, uint32_t sample_ct, void* __restrict result) {
  const __m128i* table_alias = S_CAST(const __m128i*, table256x4bx4);
  const unsigned char* genoarr_alias = R_CAST(const unsigned char*, genoarr);
  __m128i* result_alias = S_CAST(__m128i*, result);
  const uint32_t full_byte_ct = sample_ct / 4;
  for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
    _mm_storeu_si128(&(result_alias[byte_idx]), table_alias[genoarr_alias[byte_idx]]);
  }
  const uint32_t remainder = sample_ct % 4;
  if (remainder) {
    uint32_t* result_last = R_CAST(uint32_t*, &(result_alias[full_byte_ct]));
    uintptr_t geno_byte = genoarr_alias[full_byte_ct];
    for (uint32_t uii = 0; uii != remainder; ++uii) {
      memcpy(&(result_last[uii]), &(table_alias[geno_byte & 3]), 4);
      geno_byte >>= 2;
    }
  }
}
#else
void GenoarrLookup4x16b(const uintptr_t* genoarr, const void* table4x16b, uint32_t sample_ct, void* result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table4x16b);
  uint64_t* result_iter = S_CAST(uint64_t*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD2;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      memcpy(result_iter, &(table_alias[(geno_word & 3) * 2]), 16);
      result_iter = &(result_iter[2]);
      geno_word >>= 2;
    }
  }
}

void GenoarrLookup16x8bx2(const uintptr_t* genoarr, const void* table16x8bx2, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table16x8bx2);
  uint64_t* result_iter = S_CAST(uint64_t*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
          memcpy(result_iter, &(table_alias[(geno_word & 3) * 2]), 8);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
      memcpy(result_iter, &(table_alias[cur_2geno * 2]), 16);
      result_iter = &(result_iter[2]);
      geno_word >>= 4;
    }
  }
}

void GenoarrLookup256x4bx4(const uintptr_t* genoarr, const void* table256x4bx4, uint32_t sample_ct, void* __restrict result) {
  const uint32_t* table_alias = S_CAST(const uint32_t*, table256x4bx4);
  const unsigned char* genoarr_alias = R_CAST(const unsigned char*, genoarr);
  uint32_t* result_alias = S_CAST(uint32_t*, result);
  const uint32_t full_byte_ct = sample_ct / 4;
  for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
    memcpy(&(result_alias[byte_idx * 4]), &(table_alias[genoarr_alias[byte_idx] * 4]), 16);
  }
  const uint32_t remainder = sample_ct % 4;
  if (remainder) {
    uint32_t* result_last = &(result_alias[full_byte_ct * 4]);
    uintptr_t geno_byte = genoarr_alias[full_byte_ct];
    for (uint32_t uii = 0; uii != remainder; ++uii) {
      memcpy(&(result_last[uii]), &(table_alias[(geno_byte & 3) * 4]), 4);
      geno_byte >>= 2;
    }
  }
}
#endif

void InitLookup16x4bx2(void* table16x4bx2) {
  uint32_t* table_iter = S_CAST(uint32_t*, table16x4bx2);
  uint32_t vals[4];
  vals[0] = table_iter[0];
  table_iter[1] = vals[0];
  vals[1] = table_iter[2];
  table_iter[3] = vals[0];
  vals[2] = table_iter[4];
  table_iter[5] = vals[0];
  vals[3] = table_iter[6];
  table_iter[7] = vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
}

void InitLookup16x8bx2(void* table16x8bx2) {
  uint64_t* table_iter = S_CAST(uint64_t*, table16x8bx2);
  uint64_t vals[4];
  vals[0] = table_iter[0];
  table_iter[1] = vals[0];
  vals[1] = table_iter[2];
  table_iter[3] = vals[0];
  vals[2] = table_iter[4];
  table_iter[5] = vals[0];
  vals[3] = table_iter[6];
  table_iter[7] = vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    // bugfix (20 Jun 2018): cur_high needs to be a uint64_t, not a uint32_t
    const uint64_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
}

void InitLookup256x4bx4(void* table256x4bx4) {
  uint32_t* table_iter = S_CAST(uint32_t*, table256x4bx4);
  uint32_t vals[4];
  vals[0] = table_iter[0];
  vals[1] = table_iter[4];
  vals[2] = table_iter[8];
  vals[3] = table_iter[12];
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = vals[high_idx];
    for (uint32_t second_idx = 0; second_idx != 4; ++second_idx) {
      const uint32_t cur_second = vals[second_idx];
      for (uint32_t third_idx = 0; third_idx != 4; ++third_idx) {
        const uint32_t cur_third = vals[third_idx];
        for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
          *table_iter++ = vals[low_idx];
          *table_iter++ = cur_third;
          *table_iter++ = cur_second;
          *table_iter++ = cur_high;
        }
      }
    }
  }
}

void PhaseLookup4b(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const void* table56x4bx2, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table56x4bx2);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const Halfword* phasepresent_alias = R_CAST(const Halfword*, phasepresent);
  const Halfword* phaseinfo_alias = R_CAST(const Halfword*, phaseinfo);
  uint64_t* result_iter = S_CAST(uint64_t*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  uintptr_t phasepresent_hw_shifted = 0;
  uintptr_t phaseinfo_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word & 3);
          // assume trailing bits of phasepresent/phaseinfo clear
          // phaseinfo_hw_shifted not necessarily updated, so need if-statement
          // bugfix (25 Jun 2018): must only consider active bit of
          // phasepresent_hw_shifted, not the already-processed ones
          if (phasepresent_hw_shifted & 16) {
            cur_idx ^= 16 | (phaseinfo_hw_shifted & 2);
          }
          memcpy(result_iter, &(table_alias[cur_idx]), 4);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    phasepresent_hw_shifted = phasepresent_alias[widx];
    if (!phasepresent_hw_shifted) {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        *result_iter++ = table_alias[geno_word & 15];
        geno_word >>= 4;
      }
    } else {
      phasepresent_hw_shifted = phasepresent_hw_shifted << 4;
      phaseinfo_hw_shifted = phaseinfo_alias[widx];

      // note that this must be on a separate line (or we have to static_cast)
      phaseinfo_hw_shifted = phaseinfo_hw_shifted << 1;

      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uintptr_t cur_idx = ((geno_word & 15) | (phasepresent_hw_shifted & 48)) ^ (phaseinfo_hw_shifted & 6);
        *result_iter++ = table_alias[cur_idx];
        geno_word >>= 4;
        phasepresent_hw_shifted >>= 2;
        phaseinfo_hw_shifted >>= 2;
      }
    }
  }
}

void InitPhaseLookup4b(void* table56x4bx2) {
  uint32_t* table_iter = S_CAST(uint32_t*, table56x4bx2);
  uint32_t vals[4];
  vals[0] = table_iter[0];
  table_iter[1] = vals[0];
  vals[1] = table_iter[2];
  table_iter[3] = vals[0];
  vals[2] = table_iter[4];
  table_iter[5] = vals[0];
  vals[3] = table_iter[6];
  table_iter[7] = vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [16][0]..[31][1]: bit 4 is set
  // low bits must be 01 or 11
  const uint32_t val_phaseinfo0 = table_iter[2];
  table_iter[3] = vals[0];
  const uint32_t val_phaseinfo1 = table_iter[6];
  table_iter[7] = vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = vals[high_idx];
    table_iter[2] = val_phaseinfo0;
    table_iter[3] = cur_high;
    table_iter[6] = val_phaseinfo1;
    table_iter[7] = cur_high;
    table_iter = &(table_iter[8]);
  }
  // [32][0]..[39][1]: bit 5 set, bit 4 unset
  // high bits must be 00 or 01
  for (uint32_t high_idx = 0; high_idx != 2; ++high_idx) {
    const uint32_t cur_high = high_idx? val_phaseinfo0 : val_phaseinfo1;
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  table_iter = &(table_iter[16]);
  // [48][0]..[55][1]: bits 4 and 5 set
  for (uint32_t high_idx = 0; high_idx != 2; ++high_idx) {
    const uint32_t cur_high = high_idx? val_phaseinfo0 : val_phaseinfo1;
    table_iter[2] = val_phaseinfo0;
    table_iter[3] = cur_high;
    table_iter[6] = val_phaseinfo1;
    table_iter[7] = cur_high;
    table_iter = &(table_iter[8]);
  }
}

// bits 0..3: two genotypes
// bits 4..5: two (phasepresent | sex_male) bits
// bits 1,3: unpacked phaseinfo xor
void PhaseXNohhLookup4b(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const uintptr_t* sex_male, const void* table64x4bx2, uint32_t sample_ct, void* result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table64x4bx2);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const Halfword* phasepresent_alias = R_CAST(const Halfword*, phasepresent);
  const Halfword* phaseinfo_alias = R_CAST(const Halfword*, phaseinfo);
  const Halfword* sex_male_alias = R_CAST(const Halfword*, sex_male);
  uint64_t* result_iter = S_CAST(uint64_t*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word_xored = 0;
  uintptr_t male_or_phasepresent_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word_xored & 3) | (male_or_phasepresent_hw_shifted & 16);
          memcpy(result_iter, &(table_alias[cur_idx]), 4);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word_xored = genoarr[widx];
    male_or_phasepresent_hw_shifted = sex_male_alias[widx];
    const uintptr_t phasepresent_hw = phasepresent_alias[widx];
    male_or_phasepresent_hw_shifted |= phasepresent_hw;
    male_or_phasepresent_hw_shifted <<= 4;
    if (!phasepresent_hw) {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        *result_iter++ = table_alias[(geno_word_xored & 15) | (male_or_phasepresent_hw_shifted & 48)];
        geno_word_xored >>= 4;
        male_or_phasepresent_hw_shifted >>= 2;
      }
    } else {
      geno_word_xored ^= UnpackHalfwordToWordShift1(phaseinfo_alias[widx]);
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uintptr_t cur_idx = (geno_word_xored & 15) | (male_or_phasepresent_hw_shifted & 48);
        *result_iter++ = table_alias[cur_idx];
        geno_word_xored >>= 4;
        male_or_phasepresent_hw_shifted >>= 2;
      }
    }
  }
}

void InitPhaseXNohhLookup4b(void* table64x4bx2) {
  uint32_t* table_iter = S_CAST(uint32_t*, table64x4bx2);
  uint32_t vals[4];
  vals[0] = table_iter[0];
  table_iter[1] = vals[0];
  vals[1] = table_iter[2];
  table_iter[3] = vals[0];
  vals[2] = table_iter[4];
  table_iter[5] = vals[0];
  vals[3] = table_iter[6];
  table_iter[7] = vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [16][0]..[31][1]: bit 4 is set
  uint32_t male_or_phasepresent_vals[4];
  for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
    male_or_phasepresent_vals[low_idx] = *table_iter++;
    *table_iter++ = vals[0];
  }
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = male_or_phasepresent_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [32][0]..[47][1]: bit 5 set, bit 4 unset
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = male_or_phasepresent_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [48][0]..[63][1]: bits 4 and 5 set
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = male_or_phasepresent_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = male_or_phasepresent_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
}

void GenoarrSexLookup4b(const uintptr_t* genoarr, const uintptr_t* sex_male, const void* table64x4bx2, uint32_t sample_ct, void* result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table64x4bx2);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const Halfword* sex_male_alias = R_CAST(const Halfword*, sex_male);
  uint64_t* result_iter = S_CAST(uint64_t*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  uintptr_t male_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word & 3) | (male_hw_shifted & 16);
          memcpy(result_iter, &(table_alias[cur_idx]), 4);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    male_hw_shifted = sex_male_alias[widx];
    male_hw_shifted <<= 4;
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      *result_iter++ = table_alias[(geno_word & 15) | (male_hw_shifted & 48)];
      geno_word >>= 4;
      male_hw_shifted >>= 2;
    }
  }
}

void ClearGenovecMissing1bit8Unsafe(const uintptr_t* __restrict genovec, uint32_t* subset_sizep, uintptr_t* __restrict subset, void* __restrict sparse_vals) {
  const uint32_t orig_subset_size = *subset_sizep;
  Halfword* subset_alias = R_CAST(Halfword*, subset);
  uint32_t read_idx = 0;
  // deliberate overflow
  for (uint32_t read_widx = UINT32_MAX; ; ) {
    uint32_t subset_bits;
    do {
      subset_bits = subset_alias[++read_widx];
    } while (!subset_bits);
    uintptr_t detect_11 = genovec[read_widx];
    detect_11 = detect_11 & (detect_11 >> 1) & kMask5555;
    if (detect_11) {
      uint32_t detect_11_hw = PackWordToHalfword(detect_11);
      const uint32_t joint_u32 = subset_bits & detect_11_hw;
      if (joint_u32) {
        uintptr_t lowbit = joint_u32 & (-joint_u32);
        uint32_t write_idx = read_idx + PopcountWord(subset_bits & (lowbit - 1));
        read_idx = write_idx + 1;
        uint32_t subset_bits_write = subset_bits ^ lowbit;
        unsigned char* sparse_vals_uc = S_CAST(unsigned char*, sparse_vals);
        subset_bits &= -(2 * lowbit);
        for (; read_idx != orig_subset_size; ++read_idx) {
#ifdef USE_AVX2
          if (!subset_bits) {
            subset_alias[read_widx] = subset_bits_write;
            do {
              subset_bits = subset_alias[++read_widx];
            } while (!subset_bits);
            subset_bits_write = subset_bits;
            detect_11 = genovec[read_widx];
            detect_11 = detect_11 & (detect_11 >> 1);
            detect_11_hw = PackWordToHalfwordMask5555(detect_11);
          }
          lowbit = subset_bits & (-subset_bits);
          subset_bits ^= lowbit;
          if (lowbit & detect_11_hw) {
            subset_bits_write ^= lowbit;
            continue;
          }
#else
          if (!subset_bits) {
            subset_alias[read_widx] = subset_bits_write;
            do {
              subset_bits = subset_alias[++read_widx];
            } while (!subset_bits);
            subset_bits_write = subset_bits;
            detect_11 = genovec[read_widx];
            detect_11 = detect_11 & (detect_11 >> 1);
          }
          lowbit = subset_bits & (-subset_bits);
          subset_bits ^= lowbit;
          if ((lowbit * lowbit) & detect_11) {
            subset_bits_write ^= lowbit;
            continue;
          }
#endif
          sparse_vals_uc[write_idx++] = sparse_vals_uc[read_idx];
        }
        subset_alias[read_widx] = subset_bits_write;
        *subset_sizep = write_idx;
        return;
      }
    }
    read_idx += PopcountWord(subset_bits);
    if (read_idx == orig_subset_size) {
      return;
    }
  }
}

void ClearGenovecMissing1bit16Unsafe(const uintptr_t* __restrict genovec, uint32_t* subset_sizep, uintptr_t* __restrict subset, void* __restrict sparse_vals) {
  const uint32_t orig_subset_size = *subset_sizep;
  Halfword* subset_alias = R_CAST(Halfword*, subset);
  uint32_t read_idx = 0;
  // deliberate overflow
  for (uint32_t read_widx = UINT32_MAX; ; ) {
    uint32_t subset_bits;
    do {
      subset_bits = subset_alias[++read_widx];
    } while (!subset_bits);
    uintptr_t detect_11 = genovec[read_widx];
    detect_11 = detect_11 & (detect_11 >> 1) & kMask5555;
    if (detect_11) {
      uint32_t detect_11_hw = PackWordToHalfword(detect_11);
      const uint32_t joint_u32 = subset_bits & detect_11_hw;
      if (joint_u32) {
        uintptr_t lowbit = joint_u32 & (-joint_u32);
        uint32_t write_idx = read_idx + PopcountWord(subset_bits & (lowbit - 1));
        read_idx = write_idx + 1;
        uint32_t subset_bits_write = subset_bits ^ lowbit;
        uint16_t* sparse_vals_u16 = S_CAST(uint16_t*, sparse_vals);
        subset_bits &= -(2 * lowbit);
        for (; read_idx != orig_subset_size; ++read_idx) {
#ifdef USE_AVX2
          if (!subset_bits) {
            subset_alias[read_widx] = subset_bits_write;
            do {
              subset_bits = subset_alias[++read_widx];
            } while (!subset_bits);
            subset_bits_write = subset_bits;
            detect_11 = genovec[read_widx];
            detect_11 = detect_11 & (detect_11 >> 1);
            detect_11_hw = PackWordToHalfwordMask5555(detect_11);
          }
          lowbit = subset_bits & (-subset_bits);
          subset_bits ^= lowbit;
          if (lowbit & detect_11_hw) {
            subset_bits_write ^= lowbit;
            continue;
          }
#else
          if (!subset_bits) {
            subset_alias[read_widx] = subset_bits_write;
            do {
              subset_bits = subset_alias[++read_widx];
            } while (!subset_bits);
            subset_bits_write = subset_bits;
            detect_11 = genovec[read_widx];
            detect_11 = detect_11 & (detect_11 >> 1);
          }
          lowbit = subset_bits & (-subset_bits);
          subset_bits ^= lowbit;
          if ((lowbit * lowbit) & detect_11) {
            subset_bits_write ^= lowbit;
            continue;
          }
#endif
          sparse_vals_u16[write_idx++] = sparse_vals_u16[read_idx];
        }
        subset_alias[read_widx] = subset_bits_write;
        *subset_sizep = write_idx;
        return;
      }
    }
    read_idx += PopcountWord(subset_bits);
    if (read_idx == orig_subset_size) {
      return;
    }
  }
}

double MultiallelicDiploidMinimac3R2(const uint64_t* __restrict sums, const uint64_t* __restrict hap_ssqs_x2, uint32_t nm_sample_ct, uint32_t allele_ct, uint32_t extra_phased_het_ct) {
  // sums[k] == sum_i [left_dosage_{ik} + right_dosage_{ik}]
  // hap_ssqs_x2[k] ==
  //   2 * sum_i [(left_dosage_{ik})^2 + (right_dosage_{ik})^2]
  //   This may be odd, since it's computed as
  //     (left + right)^2 + (left - right)^2
  //   and the parities of the two integers can be different.
  // For phased hardcalls, it is fine for the hap_ssqs_x2[k] values to
  // correspond to unphased hardcalls iff extra_phased_het_ct is the number of
  // phased-hets that weren't accounted for in hap_ssqs_x2[]; this makes it
  // straightforward for GetMultiallelicCountsAndDosage16s to stick to the
  // custom internal multiallelic-count functions.
  if (!nm_sample_ct) {
    return (0.0 / 0.0);
  }
  // Allelic phased-dosages are on a (k-1)-dimensional simplex; embed this in
  // R^k as the (1, 0, ..., 0), (0, 1, ..., 0), ..., (0, 0, ..., 1) polytope.
  // Define
  //   m_k := (1/2n) * sum_i [left_dosage_{ik} + right_dosage_{ik}]
  // Minimac3-r2 is defined as empirical phased-dosage variance divided by
  // expected-under-allele-frequencies variance.
  // Expected sum-of-squared-Euclidean-distance with perfect imputation is
  //   2n("1"^2 - sum_k ((m_k)^2))
  // and observed sum-of-squared-distance is
  //   sum_k (sum_i [(left_dosage_{ik})^2 + (right_dosage_{ik})^2] -
  //          2n((m_k)^2))

  // ssq_sum_x2 can be as large as 2^31 * nm_sample_ct; meansq_sum can cancel
  // as little as (1 / allele_ct) of it
  if (nm_sample_ct < 92682) {
    uint64_t ssq_sum_x2 = extra_phased_het_ct * 0x20000000LLU;
    uint64_t meansq_sum = 0;
    for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
      const uint64_t cur_allele_dosage = sums[allele_idx];
      ssq_sum_x2 += hap_ssqs_x2[allele_idx];
      // cur_allele_dosage == 2n * m_k
      // -> meansq_sum becomes 2n * sum_k [2n((m_k)^2)]
      meansq_sum += cur_allele_dosage * cur_allele_dosage;
    }
    const uint64_t observed_variance_times_2n = ssq_sum_x2 * nm_sample_ct - meansq_sum;
    // "1"^2 -> 16384^2 in our units.  So 2n * 2n * "1"^2 is equal to
    //   n * n * 16384^2 * 4.
    const uint64_t expected_variance_times_2n = nm_sample_ct * 0x40000000LLU * nm_sample_ct - meansq_sum;
    // mach_r2 == 1 test cases:
    // - AA, AB, BB: 1, 4, 4
    //   sums[0] = 6 * 2^14
    //   sums[1] = 12 * 2^14
    //   ssqs[0] = 8 * 2^28
    //   ssqs[1] = 20 * 2^28
    //   ssq_sum = (8 + 20) * 2^28
    //   meansq_sum = (6 * 6 + 12 * 12) * 2^28
    //   observed_variance = 28 * 9 * 2^28 - 180 * 2^28
    //   expected_variance = (9 * 9 * 4 * 2^28 - 180 * 2^28) / 2
    // - AA, AB, BB, AC, BC, CC: 1, 4, 4, 6, 12, 9
    //   sums[0] = 12 * 2^14
    //   sums[1] = 24 * 2^14
    //   sums[2] = 36 * 2^14
    //   ssqs[0] = 14 * 2^28
    //   ssqs[1] = 32 * 2^28
    //   ssqs[2] = 54 * 2^28
    //   ssq_sum = (14 + 32 + 54) * 2^28
    //   meansq_sum = (12 * 12 + 24 * 24 + 36 * 36) * 2^28
    //   observed_variance = 100 * 36 * 2^28 - 56 * 36 * 2^28
    //   expected_variance = (36 * 36 * 4 * 2^28 - 56 * 36 * 2^28) / 2
    return S_CAST(double, observed_variance_times_2n) / S_CAST(double, expected_variance_times_2n);
  }
  // Need to avoid catastrophic cancellation here.
  const uint64_t nm_sample_ct_x32768 = nm_sample_ct * 0x8000LLU;
  double expected_variance_times_2n = 0.0;
  double observed_variance_times_2n_hi = 0.0;
  int64_t observed_variance_times_2n_lo = 0;
  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
    const uint64_t cur_allele_dosage = sums[allele_idx];
    const uint64_t cur_ssq_x2 = hap_ssqs_x2[allele_idx];
    const double cur_allele_dosaged = u63tod(cur_allele_dosage);
    // Restructure expected_variance as sum/product of
    // guaranteed-nonnegative numbers.
    //
    // cur_allele_dosage == 2n * m_k
    // cur_allele_dosage[k] * ("1" * 2n - cur_allele_dosage[k])
    // = "1" * 4n^2 * m_k - 4n^2 * m_k^2
    //
    // sum_k (cur_allele_dosage[k] * ("1" * 2n - cur_allele_dosage[k]))
    // = "1" * 4n^2 * sum_k m_k - 4n^2 * sum_k m_k^2
    // = "1" * 4n^2 * "1" - 4n^2 * sum_k m_k^2
    // = 4n^2 * ("1"^2 - sum_k m_k^2)
    // This is 2n times the expected variance.
    expected_variance_times_2n += cur_allele_dosaged * u63tod(nm_sample_ct_x32768 - cur_allele_dosage);
    // Emulate 128-bit integers to make observed_variance computation work.
    // "left" corresponds to the ssq_sum_x2 * nm_sample_ct term, "right"
    // corresponds to meansq_sum.
    const uint64_t cur_ssq_x2_hi = cur_ssq_x2 >> 32;
    uint64_t left_lo = (cur_ssq_x2 & 0xffffffffLLU) * nm_sample_ct;
    const uint64_t left_hi = (left_lo >> 32) + cur_ssq_x2_hi * nm_sample_ct;
    left_lo &= 0xffffffffU;
    const uint64_t cur_allele_dosage_lo = cur_allele_dosage & 0xffffffffLLU;
    const uint64_t cur_allele_dosage_hi = cur_allele_dosage >> 32;
    uint64_t right_lo = cur_allele_dosage_lo * cur_allele_dosage_lo;
    const uint64_t right_hi = (right_lo >> 32) + (cur_allele_dosage_lo + cur_allele_dosage) * cur_allele_dosage_hi;
    right_lo &= 0xffffffffU;
    observed_variance_times_2n_hi += u63tod(left_hi - right_hi);
    observed_variance_times_2n_lo += S_CAST(int64_t, left_lo) - S_CAST(int64_t, right_lo);
  }
  const double observed_variance_times_2n = (observed_variance_times_2n_hi * 4294967296.0) + (u31tod(extra_phased_het_ct) * 536870912.0) + observed_variance_times_2n_lo;
  return observed_variance_times_2n / expected_variance_times_2n;
}

double BiallelicDiploidMinimac3R2(uint64_t alt1_dosage, uint64_t hap_alt1_ssq_x2, uint32_t nm_sample_ct) {
  if (!nm_sample_ct) {
    return (0.0 / 0.0);
  }

  const uint64_t nm_sample_ct_x32768 = nm_sample_ct * 0x8000LLU;
  if (nm_sample_ct < 131072) {
    const uint64_t alt1_dosage_sq = alt1_dosage * alt1_dosage;
    const uint64_t observed_variance_times_2n = hap_alt1_ssq_x2 * nm_sample_ct - alt1_dosage * alt1_dosage;
    const uint64_t expected_variance_times_2n = nm_sample_ct_x32768 * alt1_dosage - alt1_dosage_sq;
    return S_CAST(double, observed_variance_times_2n) / S_CAST(double, expected_variance_times_2n);
  }
  // Need to avoid catastrophic cancellation here.
  const double alt1_dosaged = u63tod(alt1_dosage);
  const double expected_variance_times_2n = alt1_dosaged * u63tod(nm_sample_ct_x32768 - alt1_dosage);
  const uint64_t hap_alt1_ssq_x2_hi = hap_alt1_ssq_x2 >> 32;
  uint64_t left_lo = (hap_alt1_ssq_x2 & 0xffffffffLLU) * nm_sample_ct;
  const uint64_t left_hi = (left_lo >> 32) + hap_alt1_ssq_x2_hi * nm_sample_ct;
  left_lo &= 0xffffffffU;
  const uint64_t alt1_dosage_lo = alt1_dosage & 0xffffffffLLU;
  const uint64_t alt1_dosage_hi = alt1_dosage >> 32;
  uint64_t right_lo = alt1_dosage_lo * alt1_dosage_lo;
  const uint64_t right_hi = (right_lo >> 32) + (alt1_dosage_lo + alt1_dosage) * alt1_dosage_hi;
  right_lo &= 0xffffffffU;
  const double observed_variance_times_2n_hi = u63tod(left_hi - right_hi);
  const int64_t observed_variance_times_2n_lo = S_CAST(int64_t, left_lo) - S_CAST(int64_t, right_lo);
  const double observed_variance_times_2n = (observed_variance_times_2n_hi * 4294967296.0) + observed_variance_times_2n_lo;
  return observed_variance_times_2n / expected_variance_times_2n;
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

uint32_t CountPgrAllocCachelinesRequired(uint32_t raw_sample_ct, PgenGlobalFlags gflags, uint32_t max_allele_ct, uint32_t fread_buf_byte_ct) {
  // ldbase_raw_genovec: always needed, 2 bits per entry, up to raw_sample_ct
  // entries
  const uint32_t genovec_cacheline_req = QuaterCtToCachelineCt(raw_sample_ct);
  const uint32_t bitvec_cacheline_req = BitCtToCachelineCt(raw_sample_ct);
  uint32_t cachelines_required = genovec_cacheline_req;
  // fread_buf.  fread_buf_byte_ct should be zero if mmap() is being used.
  // DivUp() won't overflow since fread_buf_byte_ct requirement can't exceed
  // kPglMaxBytesPerVariant, which is sufficiently far from 2^32.
  cachelines_required += DivUp(fread_buf_byte_ct, kCacheline);

  const uint32_t ld_compression_present = (gflags / kfPgenGlobalLdCompressionPresent) & 1;
  const uint32_t max_difflist_entry_ct_base = (raw_sample_ct / kPglMaxDifflistLenDivisor);
  if ((gflags & kfPgenGlobalDifflistOrLdPresent) || (max_allele_ct > 2)) {
    // workspace_difflist_sample_ids
    // bugfix: must add 1 since several routines add a terminator element
    cachelines_required += 1 + (max_difflist_entry_ct_base / kInt32PerCacheline);
  }
  if (gflags & kfPgenGlobalDifflistOrLdPresent) {
    // const uint32_t max_difflist_entry_ct = max_difflist_entry_ct_base * (1 + ld_compression_present);
    // workspace_raregeno_vec
    cachelines_required += QuaterCtToCachelineCt(max_difflist_entry_ct_base);

    // workspace_raregeno_tmp_loadbuf
    cachelines_required += QuaterCtToCachelineCt(max_difflist_entry_ct_base);

    if (ld_compression_present) {
      // ldbase_genovec
      cachelines_required += genovec_cacheline_req;

      // ldbase_raregeno
      cachelines_required += QuaterCtToCachelineCt(max_difflist_entry_ct_base);

      // ldbase_difflist_sample_ids
      cachelines_required += 1 + (max_difflist_entry_ct_base / kInt32PerCacheline);
    }
  }
  const PgenGlobalFlags gflags_hphase_dosage = gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent);
  if ((max_allele_ct > 2) || gflags_hphase_dosage) {
    cachelines_required += genovec_cacheline_req;  // workspace_vec
    if (max_allele_ct > 2) {
      // workspace_aux1x_present
      cachelines_required += bitvec_cacheline_req;
      // workspace_imp_r2
      cachelines_required += Int64CtToCachelineCt(2 * max_allele_ct);
    }
    if (gflags & kfPgenGlobalHardcallPhasePresent) {
      // workspace_all_hets, workspace_subset
      cachelines_required += bitvec_cacheline_req * 2;
    }
    if (gflags & kfPgenGlobalDosagePresent) {
      // aux track #3: usually bitarray tracking which samples have dosage info
      // (may be stored on disk as a dosage list)
      cachelines_required += bitvec_cacheline_req;
      if (gflags & kfPgenGlobalDosagePhasePresent) {
        // aux track #7: bitarray tracking which dosage entries are phased
        cachelines_required += bitvec_cacheline_req;

        // phased aux tracks #4,8: 2 bytes per sample
        // There may be overflow risk here in the future.
        // (commented out since caller always provides this buffer for now)
        // cachelines_required += DivUp(2 * k1LU * raw_sample_ct, kCacheline);
      }
      // unphased aux track #4: 2 bytes per sample
      // cachelines_required += DivUp(2 * k1LU * raw_sample_ct, kCacheline);

      // may need deltalist64 workspace in multiallelic dosage case
    }
  }
  return cachelines_required;
}

static_assert(kPglMaxAltAlleleCt == 254, "Need to update PgfiInitPhase1().");
PglErr PgfiInitPhase1(const char* fname, uint32_t raw_variant_ct, uint32_t raw_sample_ct, uint32_t use_mmap, PgenHeaderCtrl* header_ctrl_ptr, PgenFileInfo* pgfip, uintptr_t* pgfi_alloc_cacheline_ct_ptr, char* errstr_buf) {
  pgfip->var_fpos = nullptr;
  pgfip->vrtypes = nullptr;
  pgfip->allele_idx_offsets = nullptr;
  pgfip->nonref_flags = nullptr;

  // Caller is currently expected to reset max_allele_ct if allele_idx_offsets
  // is preloaded... need to fix this interface.
  pgfip->max_allele_ct = 2;
  // pgfip->max_dosage_allele_ct = 0;

  pgfip->block_base = nullptr;
  // this should force overflow when value is uninitialized.
  pgfip->block_offset = 1LLU << 63;

  uint64_t fsize;
  const unsigned char* fread_ptr;
  FILE* shared_ff = nullptr;
  unsigned char small_readbuf[3];
#ifdef NO_MMAP
  if (unlikely(use_mmap)) {
    pgfip->shared_ff = nullptr;  // this must be initialized before block_base
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase1() use_mmap parameter is nonzero, but pgenlib was not compiled with mmap support.\n");
    return kPglRetImproperFunctionCall;
  }
#else
  if (use_mmap) {
    pgfip->shared_ff = nullptr;  // this must be initialized before block_base
    int32_t file_handle = open(fname, O_RDONLY);
    if (unlikely(file_handle < 0)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Failed to open %s : %s.\n", fname, strerror(errno));
      return kPglRetOpenFail;
    }
    struct stat statbuf;
    if (unlikely(fstat(file_handle, &statbuf) < 0)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Failed to open %s : %s.\n", fname, strerror(errno));
      return kPglRetOpenFail;
    }
    fsize = statbuf.st_size;
    pgfip->block_offset = 0;
    pgfip->file_size = fsize;
    pgfip->block_base = S_CAST(const unsigned char*, mmap(0, pgfip->file_size, PROT_READ, MAP_SHARED, file_handle, 0));
    if (unlikely(R_CAST(uintptr_t, pgfip->block_base) == (~k0LU))) {
      pgfip->block_base = nullptr;
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s read failure: %s.\n", fname, strerror(errno));
      return kPglRetReadFail;
    }
    // this provided less than a ~5% boost on OS X; mmap still took >80% longer
    // than fread on an 85GB file there
    // try MAP_POPULATE on Linux?
    // madvise((unsigned char*)(pgfip->block_base), fsize, MADV_SEQUENTIAL);
    close(file_handle);
    // update (7 Jan 2018): drop support for zero-sample and zero-variant
    // files, not worth the development cost
    if (unlikely(fsize < 4)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s is too small to be a valid .pgen file.\n", fname);
      return kPglRetMalformedInput;
    }
    fread_ptr = pgfip->block_base;
  }
#endif
  else {
    shared_ff = fopen(fname, FOPEN_RB);
    pgfip->shared_ff = shared_ff;
    if (unlikely(!shared_ff)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Failed to open %s : %s.\n", fname, strerror(errno));
      return kPglRetOpenFail;
    }
    if (unlikely(fseeko(shared_ff, 0, SEEK_END))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s read failure: %s.\n", fname, strerror(errno));
      return kPglRetReadFail;
    }
    fsize = ftello(shared_ff);
    if (unlikely(fsize < 4)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s is too small to be a valid .pgen file.\n", fname);
      return kPglRetMalformedInput;
    }
    rewind(shared_ff);
    if (unlikely(!fread_unlocked(small_readbuf, 3, 1, shared_ff))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s read failure: %s.\n", fname, strerror(errno));
      return kPglRetReadFail;
    }
    fread_ptr = small_readbuf;
  }
  // deliberate underflow
  if (unlikely(((raw_variant_ct - 1) > 0x7ffffffc) && (raw_variant_ct != UINT32_MAX))) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid raw_variant_ct function parameter.\n");
    return kPglRetImproperFunctionCall;
  }
  if (unlikely(((raw_sample_ct - 1) > 0x7ffffffd) && (raw_sample_ct != UINT32_MAX))) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid raw_sample_ct function parameter.\n");
    return kPglRetImproperFunctionCall;
  }
  if (unlikely(!memequal_k(fread_ptr, "l\x1b", 2))) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s is not a .pgen file (first two bytes don't match the magic number).\n", fname);
    return kPglRetMalformedInput;
  }
  const uint32_t file_type_code = fread_ptr[2];
  *header_ctrl_ptr = 0;
  if (file_type_code < 2) {
    // plink 1 binary
    if (unlikely(!file_type_code)) {
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
    if (unlikely(raw_sample_ct == UINT32_MAX)) {
      // either .fam must be loaded first, or user must provide sample count
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase1() must be called with an accurate raw_sample_ct value, since %s is a PLINK 1 .bed file.\n", fname);
      return kPglRetImproperFunctionCall;
    }
    const uint32_t const_vrec_width = QuaterCtToByteCt(raw_sample_ct);
    if (raw_variant_ct == UINT32_MAX) {
      // allow raw_variant_ct to be inferred
      uint64_t quotient = (fsize - 3) / const_vrec_width;
      if (unlikely((quotient > 0x7fffffffU) || (quotient * const_vrec_width + 3 != fsize))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Unexpected PLINK 1 .bed file size (since raw_sample_ct was %u, [file size - 3] should be divisible by %u and the quotient should be smaller than 2^31).\n", raw_sample_ct, const_vrec_width);
        return kPglRetMalformedInput;
      }
      raw_variant_ct = quotient;
    } else {
      if (unlikely(S_CAST(uint64_t, raw_variant_ct) * const_vrec_width + 3 != fsize)) {
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

  if (unlikely(fsize < 12)) {
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
    if (unlikely(
            (!fread_unlocked(&(pgfip->raw_variant_ct), sizeof(int32_t), 1, shared_ff)) ||
            (!fread_unlocked(&(pgfip->raw_sample_ct), sizeof(int32_t), 1, shared_ff)) ||
            (!fread_unlocked(header_ctrl_ptr, 1, 1, shared_ff)))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s read failure: %s.\n", fname, strerror(errno));
      return kPglRetReadFail;
    }
#ifndef NO_MMAP
  }
#endif
  PgenHeaderCtrl header_ctrl = *header_ctrl_ptr;
  if (raw_variant_ct == UINT32_MAX) {
    raw_variant_ct = pgfip->raw_variant_ct;
    // deliberate underflow
    if (unlikely((raw_variant_ct - 1) > 0x7ffffffc)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid variant count in .pgen file.\n");
      return kPglRetMalformedInput;
    }
  } else if (unlikely(raw_variant_ct != pgfip->raw_variant_ct)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase1() was called with raw_variant_ct == %u, but %s contains %u variant%s.\n", raw_variant_ct, fname, pgfip->raw_variant_ct, (pgfip->raw_variant_ct == 1)? "" : "s");
    return kPglRetInconsistentInput;
  }
  if (raw_sample_ct == UINT32_MAX) {
    raw_sample_ct = pgfip->raw_sample_ct;
    // deliberate underflow
    if (unlikely((raw_sample_ct - 1) > 0x7ffffffd)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid sample count in .pgen file.\n");
      return kPglRetMalformedInput;
    }
  } else if (unlikely(raw_sample_ct != pgfip->raw_sample_ct)) {
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
    if (unlikely(file_type_code > 4)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Third byte of %s does not correspond to a storage mode supported by this version of pgenlib.\n", fname);
      return kPglRetNotYetSupported;
    }
    if (unlikely(header_ctrl & 63)) {
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
    if (unlikely(S_CAST(uint64_t, raw_variant_ct) * const_vrec_width + pgfip->const_fpos_offset != fsize)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Unexpected .pgen file size (expected %" PRIu64 " bytes).\n", S_CAST(uint64_t, raw_variant_ct) * const_vrec_width + pgfip->const_fpos_offset);
      return kPglRetMalformedInput;
    }
    pgfip->const_vrtype = vrtype;
    pgfip->const_vrec_width = const_vrec_width;
    *pgfi_alloc_cacheline_ct_ptr = 0;
    return kPglRetSuccess;
  }
  if (unlikely(file_type_code >= 0x11)) {
    // todo: 0x11 phase sets (maybe not before 2021, though)
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Third byte of %s does not correspond to a storage mode supported by this version of pgenlib.\n", fname);
    return kPglRetNotYetSupported;
  }
  // plink 2 binary, general-purpose
  pgfip->const_vrtype = UINT32_MAX;
  pgfip->const_vrec_width = 0;
  const uintptr_t alt_allele_ct_byte_ct = (header_ctrl >> 4) & 3;
  if (unlikely(alt_allele_ct_byte_ct > 1)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: This version of pgenlib does not support >254 alternate alleles for a single variant.\n");
    return kPglRetNotYetSupported;
  }

  // 8 extra bytes per vblock, to support fast random access
  const uintptr_t vblock_ct = DivUp(raw_variant_ct, kPglVblockSize);

  uint64_t vrtype_and_vrec_len_bit_cost;
  if (header_ctrl & 8) {
    // Special header_ctrl modes:
    //   8: 1 bit per fused vrtype-length.  Unset = vrtype 5, set = vrtype 0.
    //   9: 2 bits, multiallelic.  0 = vrtype 5, 1 = vrtype 0, 2-3 = vrtype
    //      8 with that many more bytes than vrtype 0.  Note that this is
    //      limited to 16 ALT alleles.
    //   10: 2 bits, phased.  0 = vrtype 5, 1 = vrtype 0, 2-3 = vrtype 16
    //       with that many minus 1 bytes beyond vrtype 0.  While this is also
    //       aimed at the single-sample use case, it technically supports up to
    //       15 always-phased or 7 partially-phased samples.
    //   11: 4 bits, multiallelic + phased.  0 = vrtype 5, 1 = vrtype 0,
    //       2-7 = vrtype 8 with that many bytes beyond vrtype 0, 9 = vrtype 16
    //       phase info requiring just 1 byte, 10-15 = vrtype 24 with (x-7)
    //       extra bytes required between multiallelic and phased tracks.
    //   12: 2 bits, dosage, must be single-sample.  0 = vrtype 5,
    //       1 = vrtype 0, 2 = vrtype 0x45 with 2 bytes, 3 = vrtype 0x40 with 3
    //       total bytes.
    //   13: reserved for single-sample multiallelic + dosage.
    //   14: 4 bits, phased + dosage, must be single-sample.  0 and 1 as usual,
    //       3 = vrtype 16 with 1 phaseinfo byte, 4 = vrtype 0x45 with 2 bytes,
    //       5 = vrtype 0x40 with 3 total bytes, 12 = vrtype 0xc5 with 4 total
    //       bytes, 13 = vrtype 0xc0 with 5 total bytes, 15 = vrtype 0xe0 with
    //       6 total bytes
    //   15: reserved for single-sample multiallelic + phased dosage.
    const uint32_t header_ctrl_low3 = header_ctrl & 7;
    // this can be a table lookup once 13/15 are implemented
    if (!header_ctrl_low3) {
      vrtype_and_vrec_len_bit_cost = 1;
    } else if ((header_ctrl_low3 == 3) || (header_ctrl_low3 == 6)) {
      vrtype_and_vrec_len_bit_cost = 4;
    } else if (likely(header_ctrl_low3 <= 4)) {
      vrtype_and_vrec_len_bit_cost = 2;
    } else {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Twelfth byte of %s does not correspond to a format supported by this version of pgenlib.\n", fname);
      return kPglRetNotYetSupported;
    }
  } else {
    // set this to *4* if true, 0 if false
    const uint32_t phase_or_dosage_present_x4 = header_ctrl & 4;
    // vrtype entries = 4 bits if no phase/dosage, 8 otherwise
    // var_fpos entries = 8 + (8 * (header_ctrl & 3)) bits
    vrtype_and_vrec_len_bit_cost = 12 + phase_or_dosage_present_x4 + 8 * (header_ctrl & 3);
  }
  pgfip->const_fpos_offset += (raw_sample_ct * vrtype_and_vrec_len_bit_cost + 7) / 8 + (raw_sample_ct * alt_allele_ct_byte_ct) + (8 * vblock_ct);
  *pgfi_alloc_cacheline_ct_ptr = CountPgfiAllocCachelinesRequired(raw_variant_ct);
  return kPglRetSuccess;
}

static_assert(kPglMaxAltAlleleCt == 254, "Need to update PgfiInitPhase2().");
PglErr PgfiInitPhase2(PgenHeaderCtrl header_ctrl, uint32_t allele_cts_already_loaded, uint32_t nonref_flags_already_loaded, uint32_t use_blockload, uint32_t vblock_idx_start, uint32_t vidx_end, uint32_t* max_vrec_width_ptr, PgenFileInfo* pgfip, unsigned char* pgfi_alloc, uintptr_t* pgr_alloc_cacheline_ct_ptr, char* errstr_buf) {
  // *max_vrec_width_ptr technically only needs to be set in single-variant
  // fread() mode, but its computation is not currently optimized out in the
  // other two modes.

  // possible todo: add option to skip validation when allele_cts/nonref_flags
  // are already loaded.  but let's play it safe for now.
  const uint32_t raw_variant_ct = pgfip->raw_variant_ct;
  const uint32_t const_vrec_width = pgfip->const_vrec_width;
  *pgr_alloc_cacheline_ct_ptr = 0;

  // Note that this is a rather hefty stack allocation.
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
    // always need ldbase_raw_genovec
    *pgr_alloc_cacheline_ct_ptr = QuaterCtToCachelineCt(pgfip->raw_sample_ct);
    *max_vrec_width_ptr = const_vrec_width;
#ifdef NO_MMAP
    assert(shared_ff);
#else
    if (!shared_ff) {
      if (unlikely(use_blockload)) {
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
        if (unlikely(!memequal(nonref_flags_iter, fread_ptr, nonref_flags_byte_ct))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
          return kPglRetInconsistentInput;
        }
        return kPglRetSuccess;
      }
      if (header_ctrl & 64) {
        // all ref
        if (unlikely(!AllWordsAreZero(pgfip->nonref_flags, BitCtToWordCt(raw_variant_ct)))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
          return kPglRetInconsistentInput;
        }
        return kPglRetSuccess;
      }
      // all nonref
      if (unlikely(!AllBitsAreOne(pgfip->nonref_flags, raw_variant_ct))) {
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
        if (unlikely(!AllWordsAreZero(pgfip->nonref_flags, BitCtToWordCt(raw_variant_ct)))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
          return kPglRetInconsistentInput;
        }
      }
      return kPglRetSuccess;
    }
    if ((header_ctrl >> 6) == 2) {
      // all nonref
      if (nonref_flags_already_loaded) {
        if (unlikely(!AllBitsAreOne(pgfip->nonref_flags, raw_variant_ct))) {
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
    for (; ; nonref_flags_iter = &(nonref_flags_iter[cur_byte_ct])) {
      if (nonref_flags_iter >= nonref_flags_last) {
        if (nonref_flags_iter > nonref_flags_last) {
          return kPglRetSuccess;
        }
        cur_byte_ct = 1 + ((raw_variant_ct - 1) % (kPglVblockSize * 32)) / CHAR_BIT;
      }
      unsigned char* loadptr = nonref_flags_already_loaded? loadbuf : nonref_flags_iter;
      if (unlikely(!fread_unlocked(loadptr, cur_byte_ct, 1, shared_ff))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
        return kPglRetReadFail;
      }
      if (nonref_flags_already_loaded) {
        if (unlikely(!memequal(nonref_flags_iter, loadbuf, cur_byte_ct))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded nonref_flags do not match values in .pgen file.\n");
          return kPglRetInconsistentInput;
        }
      }
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
    if (unlikely(use_blockload)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: PgfiInitPhase2() cannot be called with use_blockload set when PgfiInitPhase1() had use_mmap set.\n");
      return kPglRetImproperFunctionCall;
    }
    fread_ptr = &(pgfip->block_base[12 + 8 * vblock_idx_start]);
    memcpy(&cur_fpos, fread_ptr, sizeof(int64_t));
    fread_ptr = &(fread_ptr[(vblock_ct_m1 + 1 - vblock_idx_start) * sizeof(int64_t)]);
  } else {
#endif
    if (vblock_idx_start) {
      if (unlikely(fseeko(shared_ff, vblock_idx_start * sizeof(int64_t), SEEK_CUR))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
        return kPglRetReadFail;
      }
    }
    if (unlikely(!fread_unlocked(&cur_fpos, sizeof(int64_t), 1, shared_ff))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
      return kPglRetReadFail;
    }
    // May also need to load the rest of these values in the future, if we want
    // to support dynamic insertion into a memory-mapped file.  But skip them
    // for now.
    if (unlikely(fseeko(shared_ff, (vblock_ct_m1 - vblock_idx_start) * sizeof(int64_t), SEEK_CUR))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
      return kPglRetReadFail;
    }
#ifndef NO_MMAP
  }
#endif
  const uint32_t vrtype_and_fpos_storage = header_ctrl & 15;
  const uint32_t alt_allele_ct_byte_ct = (header_ctrl >> 4) & 3;
  if (alt_allele_ct_byte_ct) {
    assert(alt_allele_ct_byte_ct == 1);
    if (unlikely(!allele_idx_offsets_iter)) {
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
      if (unlikely(fseeko(shared_ff, header_vblock_byte_ct * S_CAST(uint64_t, vblock_idx), SEEK_CUR))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
        return kPglRetReadFail;
      }
#ifndef NO_MMAP
    }
#endif
  }
  uint32_t cur_vblock_variant_ct = kPglVblockSize;
  uint32_t max_allele_ct = pgfip->max_allele_ct;
  for (; ; ++vblock_idx) {
    if (vblock_idx >= vblock_ct_m1) {
      if (vblock_idx > vblock_ct_m1) {
        // finish up
#ifndef NO_MMAP
        // now > instead of != to allow additional information to be stored
        // between header and first variant record
        if (!shared_ff) {
          if (unlikely(S_CAST(uintptr_t, fread_ptr - pgfip->block_base) > pgfip->var_fpos[0])) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid .pgen header.\n");
            return kPglRetMalformedInput;
          }
        } else {
#endif
          if (unlikely(S_CAST(uint64_t, ftello(shared_ff)) > pgfip->var_fpos[0])) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid .pgen header.\n");
            return kPglRetMalformedInput;
          }
#ifndef NO_MMAP
        }
#endif
        pgfip->var_fpos[vidx_end] = cur_fpos;
        pgfip->max_allele_ct = max_allele_ct;
        // if difflist/LD might be present, scan for them in a way that's
        // likely to terminate quickly
        PgenGlobalFlags new_gflags = kfPgenGlobal0;
        if (vrtype_and_fpos_storage != 8) {
          const uint32_t trailing_byte_ct = vidx_end & (kBytesPerVec - 1);
          if (trailing_byte_ct) {
            memset(&(pgfip->vrtypes[vidx_end]), 0, kBytesPerVec - trailing_byte_ct);
          }
          const VecW* vrtypes_alias_start = R_CAST(VecW*, pgfip->vrtypes);
          const VecW* vrtypes_alias_end = &(vrtypes_alias_start[DivUp(vidx_end, kBytesPerVec)]);
          if (vblock_idx_start) {
            vrtypes_alias_start = &(vrtypes_alias_start[vblock_idx_start * (kPglVblockSize / kBytesPerVec)]);
          }
          const VecW* vrtypes_alias_iter = vrtypes_alias_start;
          if (vrtype_and_fpos_storage < 8) {
            for (; vrtypes_alias_iter != vrtypes_alias_end; ++vrtypes_alias_iter) {
              const VecW cur_vvec = *vrtypes_alias_iter;
#ifdef __LP64__
              const VecW cur_vvec_bit2 = vecw_slli(cur_vvec, 5);
              const VecW cur_vvec_bit1 = vecw_slli(cur_vvec, 6);
              // check if any vrtype has bit 1 set and bit 2 clear
              if (vecw_movemask(cur_vvec_bit1 & (~cur_vvec_bit2))) {
                new_gflags |= kfPgenGlobalLdCompressionPresent | kfPgenGlobalDifflistOrLdPresent;
                break;
              }
              const VecW cur_vvec_bit0 = vecw_slli(cur_vvec, 7);
              if (vecw_movemask(cur_vvec_bit0 | cur_vvec_bit2)) {
                // this catches onebit
                new_gflags |= kfPgenGlobalDifflistOrLdPresent;
              }
#else
              const uintptr_t cur_vvec_shifted = cur_vvec >> 1;
              // check if any vrtype has bit 1 set and bit 2 clear
              if (cur_vvec & (~cur_vvec_shifted) & (2 * kMask0101)) {
                new_gflags |= kfPgenGlobalLdCompressionPresent | kfPgenGlobalDifflistOrLdPresent;
                break;
              }
              if (cur_vvec & (5 * kMask0101)) {
                // this catches onebit
                new_gflags |= kfPgenGlobalDifflistOrLdPresent;
              }
#endif
            }
          }
          if (vrtype_and_fpos_storage >= 4) {
            // Likely for one of {hphase, dosage} to be present without the
            // other; make this scan faster in that case, at the cost of
            // failing to early-exit when both are present.
            // This is also suboptimal for the vrtype_and_fpos_storage > 8
            // special encodings.
            VecW or_vvec = vecw_setzero();
            for (vrtypes_alias_iter = vrtypes_alias_start; vrtypes_alias_iter != vrtypes_alias_end; ++vrtypes_alias_iter) {
              or_vvec |= *vrtypes_alias_iter;
            }
#ifdef __LP64__
            const VecW or_vvec_bit3 = vecw_slli(or_vvec, 4);
            if (vecw_movemask(or_vvec_bit3)) {
              // note that, if no phase or dosage data is present, we don't
              // look for multiallelic hardcalls.
              new_gflags |= kfPgenGlobalMultiallelicHardcallFound;
            }
            const VecW or_vvec_bit4 = vecw_slli(or_vvec, 3);
            if (vecw_movemask(or_vvec_bit4)) {
              new_gflags |= kfPgenGlobalHardcallPhasePresent;
            }
            const VecW or_vvec_bit5 = vecw_slli(or_vvec, 2);
            const VecW or_vvec_bit6 = vecw_slli(or_vvec, 1);
            if (vecw_movemask(or_vvec_bit5 | or_vvec_bit6)) {
              new_gflags |= kfPgenGlobalDosagePresent;
              if (vecw_movemask(or_vvec)) {
                new_gflags |= kfPgenGlobalDosagePhasePresent;
              }
            }
#else
            if (or_vvec & (8 * kMask0101)) {
              new_gflags |= kfPgenGlobalMultiallelicHardcallFound;
            }
            if (or_vvec & (0x10 * kMask0101)) {
              new_gflags |= kfPgenGlobalHardcallPhasePresent;
            }
            if (or_vvec & (0x60 * kMask0101)) {
              new_gflags |= kfPgenGlobalDosagePresent;
              if (or_vvec & (0x80 * kMask0101)) {
                new_gflags |= kfPgenGlobalDosagePhasePresent;
              }
            }
#endif
          }
          if (vrtype_and_fpos_storage > 8) {
            if (vrtype_and_fpos_storage == 12) {
              max_vrec_width = 3;
            } else if (vrtype_and_fpos_storage == 14) {
              max_vrec_width = 6;
            } else {
              max_vrec_width = QuaterCtToByteCt(raw_sample_ct);
              if (vrtype_and_fpos_storage == 9) {
                max_vrec_width += 3;
              } else if (vrtype_and_fpos_storage == 10) {
                max_vrec_width += 2;
              } else {
                // 11
                max_vrec_width += 8;
              }
              // 13 and 15 not specified yet
            }
          } else if (!(vrtype_and_fpos_storage & 3)) {
            // 1 byte per vrec_len entry, don't bother to determine true
            // maximum
            max_vrec_width = 255;
          }
          pgfip->gflags |= new_gflags;
        } else {
          // vrtype_and_fpos_storage == 8.
          max_vrec_width = QuaterCtToByteCt(raw_sample_ct);
        }
        *pgr_alloc_cacheline_ct_ptr = CountPgrAllocCachelinesRequired(raw_sample_ct, new_gflags, max_allele_ct, (shared_ff && (!use_blockload))? max_vrec_width : 0);
        *max_vrec_width_ptr = max_vrec_width;
        return kPglRetSuccess;
      }
      cur_vblock_variant_ct = ModNz(vidx_end, kPglVblockSize);
    }
    // 1. handle vrtypes and var_fpos.
    if (vrtype_and_fpos_storage >= 8) {
      // Special encodings.
      uint32_t log2_entry_bit_width = 1;
      unsigned char vrtype_table[16];
      uint32_t vrec_len_table[16];
      vrtype_table[0] = 5;
      vrtype_table[1] = 0;
      vrec_len_table[0] = 0;
      const uint32_t raw_sample_ct4 = QuaterCtToByteCt(raw_sample_ct);
      vrec_len_table[1] = raw_sample_ct4;
      if (vrtype_and_fpos_storage == 8) {
        log2_entry_bit_width = 0;
      } else if (vrtype_and_fpos_storage == 9) {
        vrtype_table[2] = 8;
        vrtype_table[3] = 8;
        vrec_len_table[2] = raw_sample_ct4 + 2;
        vrec_len_table[3] = raw_sample_ct4 + 3;
      } else if (vrtype_and_fpos_storage == 10) {
        vrtype_table[2] = 16;
        vrtype_table[3] = 16;
        vrec_len_table[2] = raw_sample_ct4 + 1;
        vrec_len_table[3] = raw_sample_ct4 + 2;
      } else if (vrtype_and_fpos_storage == 11) {
        log2_entry_bit_width = 2;
        vrtype_table[2] = 8;
        vrtype_table[3] = 8;
        vrtype_table[4] = 8;
        vrtype_table[5] = 8;
        vrtype_table[6] = 8;
        vrtype_table[7] = 8;
        // 8 invalid
        vrtype_table[9] = 16;
        vrtype_table[10] = 24;
        vrtype_table[11] = 24;
        vrtype_table[12] = 24;
        vrtype_table[13] = 24;
        vrtype_table[14] = 24;
        vrtype_table[15] = 24;
        vrec_len_table[9] = raw_sample_ct4 + 1;
        for (uint32_t uii = 2; uii < 8; ++uii) {
          vrec_len_table[uii] = raw_sample_ct4 + uii;
          vrec_len_table[uii + 8] = raw_sample_ct4 + 1 + uii;
        }
      } else if (vrtype_and_fpos_storage == 12) {
        assert(raw_sample_ct == 1);
        vrtype_table[2] = 0x45;
        vrtype_table[3] = 0x40;
        vrec_len_table[2] = 2;
        vrec_len_table[3] = 3;
      } else {
        // 14 is only remaining possibility for now
        assert(raw_sample_ct == 1);
        log2_entry_bit_width = 2;
        vrtype_table[3] = 0x10;
        vrtype_table[4] = 0x45;
        vrtype_table[5] = 0x40;
        vrtype_table[12] = 0xc5;
        vrtype_table[13] = 0xc0;
        vrtype_table[15] = 0xe0;
        vrec_len_table[3] = 2;
        vrec_len_table[4] = 2;
        vrec_len_table[5] = 3;
        vrec_len_table[12] = 4;
        vrec_len_table[13] = 5;
        vrec_len_table[15] = 6;
      }
      const uint32_t entry_bit_width = 1 << log2_entry_bit_width;
      const uint32_t entry_mask = (1 << entry_bit_width) - 1;
      const uint32_t cur_byte_ct = 1 + ((cur_vblock_variant_ct - 1) >> (3 - log2_entry_bit_width));
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
        if (unlikely(!fread_unlocked(loadbuf, cur_byte_ct, 1, shared_ff))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
          return kPglRetReadFail;
        }
        loadbuf_iter = R_CAST(const uintptr_t*, loadbuf);
#ifndef NO_MMAP
      }
#endif
      const uint32_t log2_entries_per_word = kBitsPerWordLog2 - log2_entry_bit_width;
      const uint32_t block_len = 1 << log2_entries_per_word;
      uint32_t cur_vblock_idx = 0;
      uint32_t cur_vblock_idx_stop = block_len;
      for (; ; cur_vblock_idx_stop += block_len) {
        if (cur_vblock_idx_stop > cur_vblock_variant_ct) {
          if (cur_vblock_idx == cur_vblock_variant_ct) {
            break;
          }
          cur_vblock_idx_stop = cur_vblock_variant_ct;
        }
        uintptr_t input_word = *loadbuf_iter++;
        for (; cur_vblock_idx != cur_vblock_idx_stop; ++cur_vblock_idx) {
          const uint32_t input_word_masked = input_word & entry_mask;
          *vrtypes_iter++ = vrtype_table[input_word_masked];
          *var_fpos_iter++ = cur_fpos;
          cur_fpos += vrec_len_table[input_word_masked];
          input_word >>= entry_bit_width;
        }
      }
    } else {
      if (vrtype_and_fpos_storage < 4) {
        // no phase or dosage present, 4-bit vrtypes
        const uint32_t cur_byte_ct = DivUp(cur_vblock_variant_ct, 2);
#ifndef NO_MMAP
        if (shared_ff) {
#endif
          if (unlikely(!fread_unlocked(loadbuf, cur_byte_ct, 1, shared_ff))) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
            return kPglRetReadFail;
          }
          fread_ptr = loadbuf;
#ifndef NO_MMAP
        }
#endif
        const uint32_t word_write_ct = DivUp(cur_vblock_variant_ct, kBytesPerWord);
        uintptr_t* vrtypes_alias_fullword = R_CAST(uintptr_t*, vrtypes_iter);
        const Halfword* loadbuf_alias_halfword = R_CAST(const Halfword*, fread_ptr);
        for (uint32_t widx = 0; widx != word_write_ct; ++widx) {
          uintptr_t ww = loadbuf_alias_halfword[widx];
#ifdef USE_AVX2
          // speed advantage is small on my Mac since compiler auto-vectorizes
          // the code below?
          vrtypes_alias_fullword[widx] = _pdep_u64(ww, kMask0F0F);
#else
#  ifdef __LP64__
          ww = (ww | (ww << 16)) & kMask0000FFFF;
#  endif
          ww = (ww | (ww << 8)) & kMask00FF;
          vrtypes_alias_fullword[widx] = (ww | (ww << 4)) & kMask0F0F;
#endif  // !USE_AVX2
        }
        const uint32_t last_word_byte_ct = cur_vblock_variant_ct % kBytesPerWord;
        vrtypes_iter = &(vrtypes_iter[cur_vblock_variant_ct]);
        if (last_word_byte_ct) {
          ProperSubwordStore(0, kBytesPerWord - last_word_byte_ct, vrtypes_iter);
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
          if (unlikely(!fread_unlocked(vrtypes_iter, cur_vblock_variant_ct, 1, shared_ff))) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
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
          ProperSubwordStore(0, kBytesPerWord - last_word_byte_ct, vrtypes_iter);
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
        if (unlikely(!fread_unlocked(loadbuf, cur_byte_ct, 1, shared_ff))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
          return kPglRetReadFail;
        }
        fread_ptr = loadbuf;
#ifndef NO_MMAP
      }
#endif
      if (bytes_per_entry == 1) {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx != cur_vblock_variant_ct; ++cur_vblock_vidx) {
          var_fpos_iter[cur_vblock_vidx] = cur_fpos;
          uint32_t cur_vrec_len = fread_ptr[cur_vblock_vidx];
          cur_fpos += cur_vrec_len;
          // no need for correct max_vrec_width
        }
      } else if (bytes_per_entry == 2) {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx != cur_vblock_variant_ct; ++cur_vblock_vidx) {
          var_fpos_iter[cur_vblock_vidx] = cur_fpos;
          uint16_t cur_vrec_len;
          memcpy_k(&cur_vrec_len, &(fread_ptr[cur_vblock_vidx * 2]), 2);
          cur_fpos += cur_vrec_len;
          if (cur_vrec_len > max_vrec_width) {
            // todo: check whether we're better off just assuming 2^16 - 1
            max_vrec_width = cur_vrec_len;
          }
        }
      } else if (bytes_per_entry == 3) {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx != cur_vblock_variant_ct; ++cur_vblock_vidx) {
          var_fpos_iter[cur_vblock_vidx] = cur_fpos;
          uint32_t cur_vrec_len;
          // safe to read a byte past the end, since that's either in loadbuf
          // or, in mmap case, we can't be at the end of a valid file
          memcpy(&cur_vrec_len, &(fread_ptr[cur_vblock_vidx * 3]), sizeof(int32_t));
          cur_vrec_len &= 0xffffff;
          cur_fpos += cur_vrec_len;
          if (cur_vrec_len > max_vrec_width) {
            max_vrec_width = cur_vrec_len;
          }
        }
      } else {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx != cur_vblock_variant_ct; ++cur_vblock_vidx) {
          var_fpos_iter[cur_vblock_vidx] = cur_fpos;
          uint32_t cur_vrec_len;
          memcpy(&cur_vrec_len, &(fread_ptr[cur_vblock_vidx * 4]), 4);
          cur_fpos += cur_vrec_len;
          if (cur_vrec_len > max_vrec_width) {
            max_vrec_width = cur_vrec_len;
          }
        }
#ifdef __LP64__
        if (unlikely(max_vrec_width > kPglMaxBytesPerVariant)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid .pgen header.\n");
          return kPglRetMalformedInput;
        }
#else
        if (unlikely(max_vrec_width > kMaxBytesPerIO)) {
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
        if (unlikely(!fread_unlocked(loadbuf, cur_vblock_variant_ct * alt_allele_ct_byte_ct, 1, shared_ff))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
          return kPglRetReadFail;
        }
        fread_ptr = loadbuf;
#ifndef NO_MMAP
      }
#endif
      // max_allele_ct scan can probably be sped up with _mm{256}_max_epu8()?
      // probably can't do much for main loop (at least in sizeof(AlleleCode)
      // == 1 case)
      if (allele_cts_already_loaded) {
        // todo: update this for multibyte AlleleCode
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx != cur_vblock_variant_ct; ++cur_vblock_vidx) {
          const uintptr_t cur_allele_idx_offset = allele_idx_offsets_iter[cur_vblock_vidx];
          const uint32_t cur_allele_ct = fread_ptr[cur_vblock_vidx];
          if (unlikely((cur_allele_idx_offset - prev_allele_idx_offset) != cur_allele_ct)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Loaded allele_idx_offsets do not match values in .pgen file.\n");
            return kPglRetInconsistentInput;
          }
          prev_allele_idx_offset = cur_allele_idx_offset;
          if (cur_allele_ct > max_allele_ct) {
            max_allele_ct = cur_allele_ct;
          }
        }
      } else {
        for (uint32_t cur_vblock_vidx = 0; cur_vblock_vidx != cur_vblock_variant_ct; ++cur_vblock_vidx) {
          const uint32_t cur_allele_ct = fread_ptr[cur_vblock_vidx];
          allele_idx_offsets_iter[cur_vblock_vidx] = prev_allele_idx_offset;
          prev_allele_idx_offset += cur_allele_ct;
          if (cur_allele_ct > max_allele_ct) {
            max_allele_ct = cur_allele_ct;
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
          if (unlikely(!memequal(nonref_flags_iter, fread_ptr, cur_byte_ct))) {
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
        if (unlikely(!fread_unlocked(loadptr, cur_byte_ct, 1, shared_ff))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
          return kPglRetReadFail;
        }
        if (nonref_flags_already_loaded) {
          if (unlikely(!memequal(nonref_flags_iter, loadbuf, cur_byte_ct))) {
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
#ifdef __LP64__
  const VecW* vrtypes_valias = R_CAST(const VecW*, vrtypes);
  const uint32_t cur_vidx_orig_remainder = cur_vidx % kBytesPerVec;
  uint32_t vidx_vec_idx = cur_vidx / kBytesPerVec;
  Vec8thUint v8ui = 0;
  if (cur_vidx_orig_remainder) {
    const VecW cur_vvec = vrtypes_valias[vidx_vec_idx];
    // non-ld: ((bit 2) OR (NOT bit 1))
    const VecW cur_vvec_bit2 = vecw_slli(cur_vvec, 5);
    const VecW inv_cur_vvec_bit1 = ~vecw_slli(cur_vvec, 6);
    v8ui = vecw_movemask(cur_vvec_bit2 | inv_cur_vvec_bit1);
    v8ui = bzhi(v8ui, cur_vidx_orig_remainder);
  }
  while (!v8ui) {
    const VecW cur_vvec = vrtypes_valias[--vidx_vec_idx];
    const VecW cur_vvec_bit2 = vecw_slli(cur_vvec, 5);
    const VecW inv_cur_vvec_bit1 = ~vecw_slli(cur_vvec, 6);
    v8ui = vecw_movemask(cur_vvec_bit2 | inv_cur_vvec_bit1);
  }
  return (vidx_vec_idx * kBytesPerVec) + bsru32(v8ui);
#else
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
      const uint32_t new_ldbase_vidx_loworder = bsrw(detect_non_ld_word) / CHAR_BIT;
      return (vidx_word_idx * kBytesPerWord) + new_ldbase_vidx_loworder;
    }
    // everything LD-compressed in the current block.  move back 8 bytes in the
    // array (or 4-bytes for 32-bit build).
    cur_vrtypes_word = vrtypes_walias[--vidx_word_idx];
  }
#endif
}

uint64_t PgfiMultireadGetCachelineReq(const uintptr_t* variant_include, const PgenFileInfo* pgfip, uint32_t variant_ct, uint32_t block_size) {
  // if block_size < kPglVblockSize, it's ideal for it to be a power of 2 (to
  // avoid unnecessary vblock crossing), but that's not required.
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
  for (uint32_t block_idx = 0; ; ++block_idx) {
    uint32_t variant_uidx_start = block_idx * block_size;
    uint32_t variant_uidx_end = variant_uidx_start + block_size;
    if (block_idx >= block_ct_m1) {
      if (block_idx > block_ct_m1) {
        break;
      }
      variant_uidx_end = raw_variant_ct;
    }
    if (variant_include) {
      variant_uidx_start = AdvBoundedTo1Bit(variant_include, variant_uidx_start, variant_uidx_end);
      if (variant_uidx_start == variant_uidx_end) {
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
    variant_uidx_start = AdvTo1Bit(variant_include, variant_uidx_start);
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
      cur_read_uidx_end = AdvTo0Bit(variant_include, variant_uidx_start);
      cur_read_end_fpos = GetPgfiFpos(pgfip, cur_read_uidx_end);
      load_variant_ct -= cur_read_uidx_end - variant_uidx_start;
      if (!load_variant_ct) {
        break;
      }
      variant_uidx_start = AdvTo1Bit(variant_include, cur_read_uidx_end);
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
    if (unlikely(fseeko(pgfip->shared_ff, cur_read_start_fpos, SEEK_SET))) {
      return kPglRetReadFail;
    }
    uintptr_t len = cur_read_end_fpos - cur_read_start_fpos;
    if (unlikely(fread_checked(K_CAST(unsigned char*, &(pgfip->block_base[cur_read_start_fpos - block_offset])), len, pgfip->shared_ff))) {
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
    if (unlikely(fname != nullptr)) {
      return kPglRetImproperFunctionCall;
    }
    pgrp->ff = nullptr;  // make sure CleanupPgr() doesn't break
  } else {
    if (pgfip->shared_ff != nullptr) {
      if (unlikely(fname == nullptr)) {
        return kPglRetImproperFunctionCall;
      }
      // move instead of close/reopen.
      pgrp->ff = pgfip->shared_ff;
      pgfip->shared_ff = nullptr;
    } else {
      pgrp->ff = fopen(fname, FOPEN_RB);
      if (unlikely(!pgrp->ff)) {
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
    if (unlikely(fseeko(pgrp->ff, seek_pos, SEEK_SET))) {
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
  const uint32_t genovec_bytes_req = QuaterCtToCachelineCt(raw_sample_ct) * kCacheline;
  pgrp->ldbase_raw_genovec = R_CAST(uintptr_t*, pgr_alloc_iter);
  pgr_alloc_iter = &(pgr_alloc_iter[genovec_bytes_req]);
  const uint32_t bitvec_bytes_req = BitCtToCachelineCt(raw_sample_ct) * kCacheline;
  const uint32_t ld_compression_present = (gflags / kfPgenGlobalLdCompressionPresent) & 1;
  const uint32_t max_difflist_entry_ct_base = (raw_sample_ct / kPglMaxDifflistLenDivisor);
  const uint32_t max_allele_ct = pgrp->fi.max_allele_ct;
  pgrp->workspace_difflist_sample_ids = nullptr;
  if ((gflags & kfPgenGlobalDifflistOrLdPresent) || (max_allele_ct > 2)) {
    pgrp->workspace_difflist_sample_ids = R_CAST(uint32_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[(1 + (max_difflist_entry_ct_base / kInt32PerCacheline)) * (kCacheline * k1LU)]);
  }
  if (gflags & kfPgenGlobalDifflistOrLdPresent) {
    // const uint32_t max_difflist_entry_ct = max_difflist_entry_ct_base * (1 + ld_compression_present);

    pgrp->workspace_raregeno_vec = R_CAST(uintptr_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[QuaterCtToCachelineCt(max_difflist_entry_ct_base) * kCacheline]);

    pgrp->workspace_raregeno_tmp_loadbuf = R_CAST(uintptr_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[QuaterCtToCachelineCt(max_difflist_entry_ct_base) * kCacheline]);

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
    pgrp->workspace_raregeno_tmp_loadbuf = nullptr;
  }
  pgrp->workspace_vec = nullptr;
  pgrp->workspace_aux1x_present = nullptr;
  pgrp->workspace_imp_r2 = nullptr;
  pgrp->workspace_all_hets = nullptr;
  pgrp->workspace_subset = nullptr;
  const PgenGlobalFlags gflags_hphase_dosage = gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent);
  if ((max_allele_ct > 2) || gflags_hphase_dosage) {
    pgrp->workspace_vec = R_CAST(uintptr_t*, pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[genovec_bytes_req]);
    if (max_allele_ct > 2) {
      pgrp->workspace_aux1x_present = R_CAST(uintptr_t*, pgr_alloc_iter);
      pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);
      pgrp->workspace_imp_r2 = R_CAST(uint64_t*, pgr_alloc_iter);
      pgr_alloc_iter = &(pgr_alloc_iter[Int64CtToCachelineCt(2 * max_allele_ct) * (kCacheline * k1LU)]);
    }
    if (gflags & kfPgenGlobalHardcallPhasePresent) {
      pgrp->workspace_all_hets = R_CAST(uintptr_t*, pgr_alloc_iter);
      pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);
      pgrp->workspace_subset = R_CAST(uintptr_t*, pgr_alloc_iter);
      pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);
    }
    pgrp->workspace_dosage_present = nullptr;
    pgrp->workspace_dphase_present = nullptr;
    if (gflags & kfPgenGlobalDosagePresent) {
      pgrp->workspace_dosage_present = R_CAST(uintptr_t*, pgr_alloc_iter);
      pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);
      if (gflags & kfPgenGlobalDosagePhasePresent) {
        pgrp->workspace_dphase_present = R_CAST(uintptr_t*, pgr_alloc_iter);
      }
      // pgr_alloc_iter = &(pgr_alloc_iter[bitvec_bytes_req]);
    }
  }
  return kPglRetSuccess;
}

void PgrPlink1ToPlink2InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // 00 -> 10, 01 -> 11, 10 -> 01, 11 -> 00
  // new low bit  = [old low] ^ [old high]
  // new high bit = ~[old high]
  // "unsafe" because trailing bits are not zeroed out.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  const VecW m1 = VCONST_W(kMask5555);
  const VecW not_m1 = VCONST_W(kMaskAAAA);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; vidx++) {
    const VecW not_cur_vec_high = (~vptr[vidx]) & not_m1;
    vptr[vidx] = (((~vptr[vidx]) & m1) ^ vecw_srli(not_cur_vec_high, 1)) | not_cur_vec_high;
  }
}

void PgrPlink2ToPlink1InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // 00 -> 11, 01 -> 10, 10 -> 00, 11 -> 01
  // new low bit  = [old low] ^ (~[old high])
  // new high bit = ~[old high]
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  const VecW not_m1 = VCONST_W(kMaskAAAA);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; vidx++) {
    VecW cur_vec = vptr[vidx];
    VecW not_cur_vec_high = (~cur_vec) & not_m1;
    vptr[vidx] = (((~not_m1) & cur_vec) ^ vecw_srli(not_cur_vec_high, 1)) | not_cur_vec_high;
  }
}

PglErr ParseDifflistHeader(const unsigned char* fread_end, uint32_t raw_sample_ct, const unsigned char** fread_pp, uintptr_t* raregeno_buf, const unsigned char** difflist_group_info_ptr, uint32_t* difflist_len_ptr) {
  // Can be used for deltalists as well: pass raregeno_buf == nullptr.
  // Trailing bits of raregeno may not be zeroed out.
  // Will need a separate 64-bit version of this for multiallelic dosages.
  const uint32_t difflist_len = GetVint31(fread_end, fread_pp);
  // moved here to address maybe-uninitialized warnings
  *difflist_group_info_ptr = *fread_pp;
  *difflist_len_ptr = difflist_len;
  if (!difflist_len) {
    return kPglRetSuccess;
  }
  if (unlikely(difflist_len > raw_sample_ct / kPglMaxDifflistLenDivisor)) {
    // automatically catches GetVint31() failure
    return kPglRetMalformedInput;
  }
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  const uint32_t difflist_index_byte_ct = group_ct * (sample_id_byte_ct + 1) - 1;
  if (PtrAddCk(fread_end, difflist_index_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  if (!raregeno_buf) {
    // for sample ID lists without 2-bit genotype info, used for sparse dosage
    return kPglRetSuccess;
  }
  const uint32_t raregeno_byte_ct = QuaterCtToByteCt(difflist_len);
  const unsigned char* raregeno_start = *fread_pp;
  if (PtrAddCk(fread_end, raregeno_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  // possible todo: just return raregeno_start, and let the caller perform this
  // copy
  memcpy(raregeno_buf, raregeno_start, raregeno_byte_ct);
  return kPglRetSuccess;
}

PglErr ParseAndSaveDifflist(const unsigned char* fread_end, uint32_t raw_sample_ct, const unsigned char** fread_pp, uintptr_t* __restrict raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr) {
  // Appropriate when we need to iterate through the difflist multiple times.
  // Other functions are more efficient if we only need to process the list
  // once.
  // Trailing bits of raregeno may not be zeroed out.
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno, &group_info_iter, difflist_len_ptr);
  uint32_t difflist_len = *difflist_len_ptr;
  // todo: check if difflist_len == 0 early exit is a net positive or negative
  // on a few test datasets
  if (reterr || (!difflist_len)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  uint32_t* difflist_sample_ids_iter = difflist_sample_ids;
  for (uint32_t difflist_remaining = difflist_len; ; ) {
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
    uintptr_t raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
    group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    while (1) {
#ifndef __LP64__
      // perform more frequent checks in 32-bit build since raw_sample_idx may
      // overflow
      // misses "small negative" malformed input, but it'll catch data
      // corruption with very high probability
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
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
    if (unlikely(raw_sample_idx >= raw_sample_ct)) {
      return kPglRetMalformedInput;
    }
#endif
  }
  return kPglRetSuccess;
}

PglErr ParseAndSaveDifflistProperSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t raw_sample_ct, const unsigned char** fread_pp, uintptr_t* __restrict raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr, uintptr_t* __restrict raregeno_workspace) {
  // Requires a PROPER subset.  Might want to just merge this with
  // ParseAndSaveDifflist() and rename appropriately.
  // Trailing bits of raregeno are zeroed out.
  uint32_t raw_difflist_len;
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno_workspace, &group_info_iter, &raw_difflist_len);
  if (reterr || (!raw_difflist_len)) {
    *difflist_len_ptr = 0;
    return reterr;
  }
  const uint32_t subgroup_idx_last = (raw_difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  uintptr_t* raregeno_workspace_iter = raregeno_workspace;
  uintptr_t* raregeno_iter = raregeno;
  uint32_t* difflist_sample_ids_iter = difflist_sample_ids;

  // technically doesn't need to be initialized, but I have principles
  uintptr_t raw_sample_idx = 0;

  uintptr_t raregeno_word = 0;
  uint32_t subgroup_len_m1 = kBitsPerWordD2 - 1;
  uint32_t difflist_len_lowbits = 0;
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
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
    // best to make the inner loop have a period of 32 (call this a 'subgroup',
    // where 'group' refers to a set of 64 entries).
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
#ifdef __LP64__
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    uintptr_t raregeno_workspace_word = *raregeno_workspace_iter++;
    for (uint32_t raw_difflist_idx_lowbits = 0; ; ++raw_difflist_idx_lowbits) {
#ifndef __LP64__
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
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
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
  }
}

PglErr ParseLdAndMergeDifflistSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict ldbase_raregeno, const uint32_t* __restrict ldbase_difflist_sample_ids, uint32_t ldbase_difflist_len, uintptr_t ldbase_common_geno, uint32_t raw_sample_ct, uint32_t sample_ct, const unsigned char** fread_pp, uintptr_t* __restrict merged_raregeno, uint32_t* __restrict merged_difflist_sample_ids, uint32_t* __restrict merged_difflist_len_ptr, uintptr_t* __restrict diff_from_ldbase_raregeno_iter) {
  // Used when the ldbase variant was saved as a difflist, and it's useful to
  // process the current variant as a difflist.
  // * Assumes ldbase_difflist_sample_ids[ldbase_difflist_len]==sample_ct.
  // * Assumes sample_include == nullptr if no subsetting needed.  (Otherwise,
  //   it'll still work, but performance will be worseslowerw.)
  // Trailing bits of merged_raregeno may not be zeroed out.
  // Caller is responsible for inverting ldbase_common_geno and merged_raregeno
  // afterward if necessary.
  assert(ldbase_difflist_sample_ids[ldbase_difflist_len] == sample_ct);
  uint32_t diff_from_ldbase_len;
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, diff_from_ldbase_raregeno_iter, &group_info_iter, &diff_from_ldbase_len);
  if (unlikely(reterr)) {
    return reterr;
  }
  if (!diff_from_ldbase_len) {
    memcpy(merged_difflist_sample_ids, ldbase_difflist_sample_ids, ldbase_difflist_len * sizeof(int32_t));
    *merged_difflist_len_ptr = ldbase_difflist_len;
    CopyQuaterarr(ldbase_raregeno, ldbase_difflist_len, merged_raregeno);
    return kPglRetSuccess;
  }
  const uint32_t subgroup_idx_last = (diff_from_ldbase_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
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
  uint32_t done = 0;
  uint32_t subgroup_len_m1 = kBitsPerWordD2 - 1;
  uint32_t merge_idx_lowbits = 0;
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
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
      raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    diff_from_ldbase_raregeno_word = *diff_from_ldbase_raregeno_iter++;
    for (; ; ++diff_from_ldbase_idx_lowbits) {
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
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
      if (diff_from_ldbase_idx_lowbits == subgroup_len_m1) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      diff_from_ldbase_raregeno_word >>= 2;
    }
  }
}

void PgrDifflistToGenovecUnsafe(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uintptr_t difflist_common_geno, uint32_t sample_ct, uint32_t difflist_len, uintptr_t* __restrict genovec) {
  // Ok for trailing bits of raregeno to be nonzero.  Does not zero out
  // trailing bits of genovec.
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  // could just memset up to word boundary; this should be a bit more
  // vector-instruction-friendly, though
  vecset(genovec, difflist_common_geno * kMask5555, vec_ct);
  const uintptr_t* raregeno_incr = raregeno;
  uint32_t difflist_idx = 0;
  uint32_t difflist_idx_stop = kBitsPerWordD2;
  if (!difflist_common_geno) {
    // faster inner loop since there's no existing value to mask out
    // todo: check if this should just be deleted since the code bloat causes
    // too many more cache misses
    for (; ; difflist_idx_stop += kBitsPerWordD2) {
      if (difflist_idx_stop > difflist_len) {
        if (difflist_idx == difflist_len) {
          return;
        }
        difflist_idx_stop = difflist_len;
      }
      uintptr_t raregeno_word = *raregeno_incr++;
      for (; difflist_idx != difflist_idx_stop; ++difflist_idx) {
        const uint32_t cur_sample_idx = difflist_sample_ids[difflist_idx];
        genovec[cur_sample_idx / kBitsPerWordD2] |= (raregeno_word & 3) << (2 * (cur_sample_idx % kBitsPerWordD2));
        raregeno_word >>= 2;
      }
    }
  }
  for (; ; difflist_idx_stop += kBitsPerWordD2) {
    if (difflist_idx_stop > difflist_len) {
      if (difflist_idx == difflist_len) {
        return;
      }
      difflist_idx_stop = difflist_len;
    }
    uintptr_t raregeno_word = *raregeno_incr++;
    for (; difflist_idx != difflist_idx_stop; ++difflist_idx) {
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
  vecset(genovec, difflist_common_geno * kMask5555, vec_ct);
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

PglErr ParseAndApplyDifflist(const unsigned char* fread_end, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effects: uses pgr.workspace_raregeno_tmp_loadbuf.
  // Cannot occur after genovec subsetting since the difflist sample indexes
  // will be incorrect.
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  uintptr_t* cur_raregeno_iter = pgrp->workspace_raregeno_tmp_loadbuf;
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, cur_raregeno_iter, &group_info_iter, &difflist_len);
  if (reterr || (!difflist_len)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  uintptr_t raw_sample_idx = 0;
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
      raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    uintptr_t cur_raregeno_word = *cur_raregeno_iter++;
    // This loop tends to be the decompression bottleneck.  Tried to modify it
    // to process 4 entries at a time, but that didn't end up helping.
    for (; ; --remaining_deltas_in_subgroup) {
      // always check, since otherwise AssignQuaterarrEntry() can scribble
      // over arbitrary memory
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
      const uintptr_t cur_geno = cur_raregeno_word & 3;
      AssignQuaterarrEntry(raw_sample_idx, cur_geno, genovec);
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      cur_raregeno_word >>= 2;
    }
  }
}

// could merge ParseAndApplyDifflist() with this?
PglErr ParseAndApplyDifflistSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effects: uses pgr.workspace_raregeno_tmp_loadbuf.
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  if (sample_ct == raw_sample_ct) {
    return ParseAndApplyDifflist(fread_end, fread_pp, pgrp, genovec);
  }
  uintptr_t* cur_raregeno_iter = pgrp->workspace_raregeno_tmp_loadbuf;
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, cur_raregeno_iter, &group_info_iter, &difflist_len);
  if (reterr || (!difflist_len)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  uintptr_t raw_sample_idx = 0;
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
      raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    uintptr_t cur_raregeno_word = *cur_raregeno_iter++;
    // This loop tends to be the decompression bottleneck.  Tried to modify it
    // to process 4 entries at a time, but that didn't end up helping.
    for (; ; --remaining_deltas_in_subgroup) {
      // always check, since otherwise AssignQuaterarrEntry() can scribble
      // over arbitrary memory
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
      if (IsSet(sample_include, raw_sample_idx)) {
        const uintptr_t cur_geno = cur_raregeno_word & 3;
        AssignQuaterarrEntry(RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, raw_sample_idx), cur_geno, genovec);
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      cur_raregeno_word >>= 2;
    }
  }
}

PglErr ParseOnebitUnsafe(const unsigned char* fread_end, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // doesn't zero out trailing genovec bits
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t common2_and_bitarray_byte_ct = (raw_sample_ct + 15) / CHAR_BIT;
  const unsigned char* onebit_main_iter = *fread_pp;
  if (PtrAddCk(fread_end, common2_and_bitarray_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uintptr_t common2_code = *onebit_main_iter++;
  const uintptr_t word_base = (common2_code / 4) * kMask5555;
  const uintptr_t common_code_delta = common2_code & 3;
  uint32_t genovec_widx = 0;
#if defined(__LP64__) && !defined(USE_AVX2)
  // this is slower in AVX2 case
  const uint32_t read_hw_ct = raw_sample_ct / kBitsPerWordD2;
  if (read_hw_ct >= 2 * kWordsPerVec) {
    const uint32_t read_vec_ct = raw_sample_ct / kBitsPerVec;
    const VecW* onebit_main_valias = R_CAST(const VecW*, onebit_main_iter);
    const VecW m4 = VCONST_W(kMask0F0F);
#  ifdef USE_SSE42
    // 0, 1, 4, 5, 16, 17, 20, 21, 64, 65, 68, 69, 80, 81, 84, 85 if the codes
    // are 0 and 1
    const VecW lookup = {word_base + common_code_delta * 0x1514111005040100LLU,
                         word_base + common_code_delta * 0x5554515045444140LLU};
#  else
    const VecW m1 = VCONST_W(kMask5555);
    const VecW m2 = VCONST_W(kMask3333);
    const VecW vec_base = VCONST_W(word_base);
    const VecW vec_delta = VCONST_W(common_code_delta * kMask5555);
#  endif
    VecW* genovec_valias = R_CAST(VecW*, genovec);
    for (uint32_t vidx = 0; vidx != read_vec_ct; ++vidx) {
      const VecW cur_vec = vecw_loadu(&(onebit_main_valias[vidx]));
      const VecW vec_even = cur_vec & m4;
      const VecW vec_odd = vecw_srli(cur_vec, 4) & m4;
      VecW vec_lo = vecw_unpacklo8(vec_even, vec_odd);
      VecW vec_hi = vecw_unpackhi8(vec_even, vec_odd);
#  ifdef USE_SSE42
      vec_lo = vecw_shuffle8(lookup, vec_lo);
      vec_hi = vecw_shuffle8(lookup, vec_hi);
#  else
      // unpack bytes, then use as mask for vec_add.
      vec_lo = (vec_lo | vecw_slli(vec_lo, 2)) & m2;
      vec_hi = (vec_hi | vecw_slli(vec_hi, 2)) & m2;
      vec_lo = (vec_lo | vecw_slli(vec_lo, 1)) & m1;
      vec_hi = (vec_hi | vecw_slli(vec_hi, 1)) & m1;
      vec_lo = vec_lo | vecw_slli(vec_lo, 1);
      vec_hi = vec_hi | vecw_slli(vec_hi, 1);
      vec_lo = vec_base + (vec_delta & vec_lo);
      vec_hi = vec_base + (vec_delta & vec_hi);
#  endif
      genovec_valias[2 * vidx] = vec_lo;
      genovec_valias[2 * vidx + 1] = vec_hi;
    }
    genovec_widx = read_vec_ct * (2 * kWordsPerVec);
  }
#endif
  const uint32_t genovec_widx_trail = (raw_sample_ct + 7) / kBitsPerWordD2;
  const uint32_t genovec_widx_end = QuaterCtToWordCt(raw_sample_ct);
#  ifdef __arm__
#    error "Unaligned accesses in ParseOnebitUnsafe()."
#  endif
  const Halfword* onebit_main_alias = R_CAST(const Halfword*, onebit_main_iter);
  for (; ; ++genovec_widx) {
    uintptr_t ww;
    if (genovec_widx >= genovec_widx_trail) {
      // might want to modify to not go here if last read is an entire halfword
      if (genovec_widx == genovec_widx_end) {
        break;
      }
      ww = ProperSubwordLoad(&(onebit_main_alias[genovec_widx_trail]), 1 + (((raw_sample_ct - 1) % kBitsPerWordD2) / CHAR_BIT));
    } else {
      ww = onebit_main_alias[genovec_widx];
    }
    // apply middle-out operation
    // 64-bit:
    //   const uintptr_t middle_out_result = (ww | (ww << 31)) & kMask5555;
    // 32-bit:
    //   *genovec_iter++ = word_base + (ww & kMask5555) * common_code_delta;
    //   *genovec_iter++ = word_base + ((ww >> 1) & kMask5555) * common_code_delta;
    // (scrapped since the time savings don't seem to be worth the extra
    // end-of-vector corner cases, apparently the extra operations here are
    // sufficiently cheap, or even negative-cost in AVX2 case)

    ww = UnpackHalfwordToWord(ww);
    genovec[genovec_widx] = word_base + ww * common_code_delta;
  }
  return ParseAndApplyDifflist(fread_end, fread_pp, pgrp, genovec);
}

PglErr Parse1or2bitGenovecUnsafe(const unsigned char* fread_end, uint32_t vrtype, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effect: may use pgrp->workspace_raregeno_tmp_loadbuf.
  // Does not update fp_vidx, does not rotate plink1-formatted data (since it's
  // better to do that post-subsetting)
  if (vrtype & 3) {
    return ParseOnebitUnsafe(fread_end, fread_pp, pgrp, genovec);
  }
  // uncompressed storage
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t genovec_byte_ct = QuaterCtToByteCt(raw_sample_ct);
  const unsigned char* src_genodata = *fread_pp;
  if (PtrAddCk(fread_end, genovec_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  memcpy(genovec, src_genodata, genovec_byte_ct);
  return kPglRetSuccess;
}

PglErr ParseNonLdGenovecSubsetUnsafe(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vrtype, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effects:
  //   may use pgrp->workspace_raregeno_tmp_loadbuf
  //   fills pgrp->ldbase_raw_genovec iff (!(vrtype & 4)) and
  //     subsetting_required (does not update ldbase_stypes, caller's
  //     responsibility to care)
  // See comments on Parse1or2bitGenovecUnsafe().
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  if (!(vrtype & 4)) {
    const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
    uintptr_t* raw_genovec = subsetting_required? pgrp->ldbase_raw_genovec : genovec;
    PglErr reterr = Parse1or2bitGenovecUnsafe(fread_end, vrtype, fread_pp, pgrp, raw_genovec);
    if ((!subsetting_required) || reterr) {
      return reterr;
    }
    CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, genovec);
    return kPglRetSuccess;
  }
  const uint32_t vrtype_low2 = vrtype & 3;
  if (vrtype_low2 != 1) {
    const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);

    // This memset is frequently the limiting operation.  This suggests that we
    // should eventually make more use of the DifflistOrGenovec interface.
    vecset(genovec, vrtype_low2 * kMask5555, vec_ct);
    return ParseAndApplyDifflistSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, fread_pp, pgrp, genovec);
  }
  // all homozygous-ref special case
  ZeroWArr(QuaterCtToWordCt(sample_ct), genovec);
  return kPglRetSuccess;
}

BoolErr InitReadPtrs(uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp) {
  const unsigned char* block_base = pgrp->fi.block_base;
  if (block_base != nullptr) {
    // possible todo: special handling of end of vblock
    const uint64_t block_offset = pgrp->fi.block_offset;
    *fread_pp = &(block_base[GetPgfiFpos(&(pgrp->fi), vidx) - block_offset]);
    *fread_endp = &(block_base[GetPgfiFpos(&(pgrp->fi), vidx + 1) - block_offset]);

    // still a useful hint to LdLoadNecessary()
    pgrp->fp_vidx = vidx + 1;

    return 0;
  }
  if (pgrp->fp_vidx != vidx) {
    if (unlikely(fseeko(pgrp->ff, GetPgfiFpos(&(pgrp->fi), vidx), SEEK_SET))) {
      return 1;
    }
  }
  const uintptr_t cur_vrec_width = GetPgfiVrecWidth(&(pgrp->fi), vidx);
#ifdef __LP64__
  if (unlikely(fread_checked(pgrp->fread_buf, cur_vrec_width, pgrp->ff))) {
    return 1;
  }
#else
  // cur_vrec_width < 2^31 since otherwise we error out on initialization
  if (unlikely(!fread_unlocked(pgrp->fread_buf, cur_vrec_width, 1, pgrp->ff))) {
    return 1;
  }
#endif
  *fread_pp = pgrp->fread_buf;
  *fread_endp = &(pgrp->fread_buf[cur_vrec_width]);
  pgrp->fp_vidx = vidx + 1;
  return 0;
}

uint32_t LdLoadNecessary(uint32_t cur_vidx, PgenReader* pgrp) {
  // Determines whether LD base variant needs to be loaded (in addition to the
  // current variant), assuming we need (possibly subsetted) hardcalls.
  // Important: this updates pgrp->ldbase_vidx when necessary, as a side
  // effect.
  // bugfix (22 May 2018): this only checked whether ldbase_stypes was nonzero;
  // there was an AllHets + cache-clear edge case where that's not good enough.
  // now that AllHets has been removed, though, it should be safe again.
  if (pgrp->ldbase_stypes && (cur_vidx == pgrp->fp_vidx)) {
    assert(pgrp->ldbase_stypes & (kfPgrLdcacheQuater | kfPgrLdcacheDifflist | kfPgrLdcacheRawQuater));
    // ldbase variant guaranteed to be up-to-date if we didn't skip the last
    // variant, and cache wasn't cleared
    return 0;
  }
  // Find the last vrtypes[] value before vrtypes[cur_vidx] with bit 1 unset or
  // bit 2 set.
  const uint32_t old_ldbase_vidx = pgrp->ldbase_vidx;
  const uint32_t new_ldbase_vidx = GetLdbaseVidx(pgrp->fi.vrtypes, cur_vidx);
  if (old_ldbase_vidx == new_ldbase_vidx) {
    return 0;
  }
  pgrp->ldbase_vidx = new_ldbase_vidx;
  return 1;
}

// Fills dest with subsetted ldbase contents, and ensures ldcache is filled so
// no explicit reload of ldbase is needed for next variant if we're extracting
// the same sample subset.  (Reload is occasionally needed if next variant is
// multiallelic or phased, we only prevent that when convenient.)
PglErr LdLoadAndCopyGenovecSubset(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* dest) {
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  if (LdLoadNecessary(vidx, pgrp)) {
    const uint32_t ldbase_vidx = pgrp->ldbase_vidx;
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (unlikely(InitReadPtrs(ldbase_vidx, pgrp, &fread_ptr, &fread_end))) {
      return kPglRetReadFail;
    }
    const uint32_t vrtype = pgrp->fi.vrtypes[ldbase_vidx];
    PglErr reterr = ParseNonLdGenovecSubsetUnsafe(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, vrtype, &fread_ptr, pgrp, dest);
    pgrp->ldbase_stypes = ((sample_ct != raw_sample_ct) && (!(vrtype & 4)))? (kfPgrLdcacheQuater | kfPgrLdcacheRawQuater) : kfPgrLdcacheQuater;
    CopyQuaterarr(dest, sample_ct, pgrp->ldbase_genovec);
    return reterr;
  }
  if (pgrp->ldbase_stypes & kfPgrLdcacheQuater) {
    CopyQuaterarr(pgrp->ldbase_genovec, sample_ct, dest);
  } else {
    if ((pgrp->ldbase_stypes & kfPgrLdcacheRawQuater) && (sample_ct == raw_sample_ct)) {
      CopyQuaterarr(pgrp->ldbase_raw_genovec, sample_ct, dest);
    } else if (pgrp->ldbase_stypes & kfPgrLdcacheDifflist) {
      // rematerialize-from-difflist is cheap.
      PgrDifflistToGenovecUnsafe(pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, pgrp->fi.vrtypes[pgrp->ldbase_vidx] & 3, sample_ct, pgrp->ldbase_difflist_len, dest);
    } else {
      CopyQuaterarrNonemptySubset(pgrp->ldbase_raw_genovec, sample_include, pgrp->fi.raw_sample_ct, sample_ct, dest);
      CopyQuaterarr(dest, sample_ct, pgrp->ldbase_genovec);
      pgrp->ldbase_stypes |= kfPgrLdcacheQuater;
    }
  }
  return kPglRetSuccess;
}

// fread_pp should be non-null iff this is being called by an internal function
// as part of a more complex read.
// in multiallelic case:
//   hom-ref = 0
//   het-ref = 1
//   two nonref = 2
//   missing = 3
PglErr ReadGenovecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict genovec) {
  // Side effects:
  //   may use pgr.workspace_raregeno_tmp_loadbuf (any difflist)
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t maintrack_vrtype = vrtype & 7;
  if (VrtypeLdCompressed(maintrack_vrtype)) {
    // LD compression
    PglErr reterr = LdLoadAndCopyGenovecSubset(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, genovec);
    if (unlikely(reterr)) {
      return reterr;
    }
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (unlikely(InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end))) {
      return kPglRetReadFail;
    }
    reterr = ParseAndApplyDifflistSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, &fread_ptr, pgrp, genovec);
    if (unlikely(reterr)) {
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
  if (unlikely(InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end))) {
    return kPglRetReadFail;
  }
  // tried to add more sophisticated caching, but turns out it isn't worth it
  PglErr reterr = ParseNonLdGenovecSubsetUnsafe(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, maintrack_vrtype, &fread_ptr, pgrp, genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  if (vrtype == kPglVrtypePlink1) {
    PgrPlink1ToPlink2InplaceUnsafe(sample_ct, genovec);
  } else {
    const uint32_t is_ldbase = pgrp->fi.vrtypes && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
    const uint32_t ldbase_raw_genovec_saved = (sample_ct != pgrp->fi.raw_sample_ct) && (!(maintrack_vrtype & 4));
    if (is_ldbase) {
      CopyQuaterarr(genovec, sample_ct, pgrp->ldbase_genovec);
      pgrp->ldbase_vidx = vidx;
      // may be better to just always set to kfPgrLdcacheQuater?  this depends
      // on multiallelic code
      pgrp->ldbase_stypes = ldbase_raw_genovec_saved? (kfPgrLdcacheQuater | kfPgrLdcacheRawQuater) : kfPgrLdcacheQuater;
    } else if (ldbase_raw_genovec_saved) {
      // bugfix (22 Sep 2018): when accessing variants out of order, need to
      // note that we just clobbered the cache
      pgrp->ldbase_stypes &= ~kfPgrLdcacheRawQuater;
    }
  }
  if (fread_pp) {
    *fread_pp = fread_ptr;
    *fread_endp = fread_end;
  }
  return kPglRetSuccess;
}

PglErr PgrGet(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  return ReadGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genovec);
}

// Fills dest with ldbase contents, and ensures ldcache is filled so no
// explicit reload of ldbase is needed for next variant.
PglErr LdLoadAndCopyRawGenovec(uint32_t subsetting_required, uint32_t vidx, PgenReader* pgrp, uintptr_t* dest) {
  const uint32_t genovec_byte_ct = QuaterCtToVecCt(pgrp->fi.raw_sample_ct) * kBytesPerVec;
  if (LdLoadNecessary(vidx, pgrp) || (subsetting_required && (!(pgrp->ldbase_stypes & kfPgrLdcacheRawQuater)))) {
    const uint32_t ldbase_vidx = pgrp->ldbase_vidx;
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (unlikely(InitReadPtrs(ldbase_vidx, pgrp, &fread_ptr, &fread_end))) {
      return kPglRetReadFail;
    }
    const uint32_t vrtype = pgrp->fi.vrtypes[ldbase_vidx];
    pgrp->ldbase_stypes = kfPgrLdcacheRawQuater;
    assert((vrtype & 7) != 5); // all-hom-ref can't be ldbase
    uintptr_t* raw_genovec = pgrp->ldbase_raw_genovec;
    PglErr reterr;
    if (!(vrtype & 4)) {
      reterr = Parse1or2bitGenovecUnsafe(fread_end, vrtype, &fread_ptr, pgrp, raw_genovec);
    } else {
      const uint32_t vrtype_low2 = vrtype & 3;
      vecset(raw_genovec, vrtype_low2 * kMask5555, DivUp(genovec_byte_ct, kBytesPerVec));
      reterr = ParseAndApplyDifflist(fread_end, &fread_ptr, pgrp, raw_genovec);
    }
    memcpy(dest, raw_genovec, genovec_byte_ct);
    return reterr;
  }
  if (pgrp->ldbase_stypes & kfPgrLdcacheRawQuater) {
    memcpy(dest, pgrp->ldbase_raw_genovec, genovec_byte_ct);
  } else {
    // no subsetting, can use regular Ldcache entries
    const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
    if (pgrp->ldbase_stypes & kfPgrLdcacheQuater) {
      memcpy(dest, pgrp->ldbase_genovec, genovec_byte_ct);
    } else {
      PgrDifflistToGenovecUnsafe(pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, pgrp->fi.vrtypes[pgrp->ldbase_vidx] & 3, raw_sample_ct, pgrp->ldbase_difflist_len, dest);
    }
  }
  return kPglRetSuccess;
}

// Does not zero out trailing bits.
// Requires fread_pp and fread_endp to be non-null for now.
PglErr ReadRawGenovec(uint32_t subsetting_required, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* raw_genovec) {
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t maintrack_vrtype = vrtype & 7;
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  if (VrtypeLdCompressed(maintrack_vrtype)) {
    // LD compression
    PglErr reterr = LdLoadAndCopyRawGenovec(subsetting_required, vidx, pgrp, raw_genovec);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (unlikely(InitReadPtrs(vidx, pgrp, fread_pp, fread_endp))) {
      return kPglRetReadFail;
    }
    reterr = ParseAndApplyDifflist(*fread_endp, fread_pp, pgrp, raw_genovec);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (maintrack_vrtype == 3) {
      GenovecInvertUnsafe(raw_sample_ct, raw_genovec);
    }
    return kPglRetSuccess;
  }
  if (unlikely(InitReadPtrs(vidx, pgrp, fread_pp, fread_endp))) {
    return kPglRetReadFail;
  }
  const unsigned char* fread_end = *fread_endp;
  PglErr reterr;
  if (!(vrtype & 4)) {
    reterr = Parse1or2bitGenovecUnsafe(fread_end, vrtype, fread_pp, pgrp, raw_genovec);
  } else {
    const uint32_t vrtype_low2 = vrtype & 3;
    if (vrtype_low2 == 1) {
      ZeroWArr(QuaterCtToWordCt(raw_sample_ct), raw_genovec);
      // all-hom-ref can't be ldbase
      return kPglRetSuccess;
    }
    const uint32_t vec_ct = QuaterCtToVecCt(raw_sample_ct);
    vecset(raw_genovec, vrtype_low2 * kMask5555, vec_ct);
    reterr = ParseAndApplyDifflist(fread_end, fread_pp, pgrp, raw_genovec);
  }
  if (vrtype == kPglVrtypePlink1) {
    PgrPlink1ToPlink2InplaceUnsafe(raw_sample_ct, raw_genovec);
  } else {
    const uint32_t is_ldbase = pgrp->fi.vrtypes && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
    if (is_ldbase) {
      CopyQuaterarr(raw_genovec, raw_sample_ct, pgrp->ldbase_raw_genovec);
      pgrp->ldbase_vidx = vidx;
      pgrp->ldbase_stypes = kfPgrLdcacheRawQuater;
    }
  }
  return reterr;
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

// Populates pgrp->ldbase_genovec or
// pgrp->ldbase_{raregeno,difflist_sample_ids,difflist_len}, depending on
// storage type.
// Currently just called by ReadDifflistOrGenovecSubsetUnsafe(), which isn't
// exploited by plink2 yet.
PglErr LdLoadMinimalSubsetIfNecessary(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp) {
  if (!LdLoadNecessary(vidx, pgrp)) {
    return kPglRetSuccess;
  }
  const uint32_t ldbase_vidx = pgrp->ldbase_vidx;
  const uint64_t cur_vidx_fpos = pgrp->fi.var_fpos[ldbase_vidx];
  const uint32_t ldbase_vrtype = pgrp->fi.vrtypes[ldbase_vidx];
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* raw_genovec = subsetting_required? pgrp->ldbase_raw_genovec : pgrp->ldbase_genovec;
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
    if (!(ldbase_vrtype & 4)) {
      reterr = Parse1or2bitGenovecUnsafe(fread_end, ldbase_vrtype, &fread_ptr, pgrp, raw_genovec);
    LdLoadMinimalSubsetIfNecessary_genovec_finish:
      pgrp->ldbase_stypes = subsetting_required? (kfPgrLdcacheQuater | kfPgrLdcacheRawQuater) : kfPgrLdcacheQuater;
      if ((!subsetting_required) || reterr) {
        return reterr;
      }
      CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, pgrp->ldbase_genovec);
      return kPglRetSuccess;
    }
    pgrp->fp_vidx = ldbase_vidx + 1;
  } else {
    if (unlikely(fseeko(pgrp->ff, pgrp->fi.var_fpos[ldbase_vidx], SEEK_SET))) {
      return kPglRetReadFail;
    }
    const uintptr_t cur_vrec_width = pgrp->fi.var_fpos[ldbase_vidx + 1] - cur_vidx_fpos;
    pgrp->fp_vidx = ldbase_vidx + 1;
    if (!(ldbase_vrtype & 7)) {
      // don't actually need to fread the whole record in this case
      const uint32_t raw_sample_ct4 = QuaterCtToByteCt(raw_sample_ct);
      if (unlikely(!fread_unlocked(raw_genovec, raw_sample_ct4, 1, pgrp->ff))) {
        return kPglRetReadFail;
      }
      if (raw_sample_ct4 != cur_vrec_width) {
        // ensure this doesn't match
        pgrp->fp_vidx = 0;
      }
      goto LdLoadMinimalSubsetIfNecessary_genovec_finish;
    }
    if (unlikely(!fread_unlocked(pgrp->fread_buf, cur_vrec_width, 1, pgrp->ff))) {
      return kPglRetReadFail;
    }
    fread_ptr = pgrp->fread_buf;
    fread_end = &(pgrp->fread_buf[cur_vrec_width]);
    if (!(ldbase_vrtype & 4)) {
      reterr = ParseOnebitUnsafe(fread_end, &fread_ptr, pgrp, raw_genovec);
      goto LdLoadMinimalSubsetIfNecessary_genovec_finish;
    }
  }
  uint32_t ldbase_difflist_len;
  if (!subsetting_required) {
    reterr = ParseAndSaveDifflist(fread_end, raw_sample_ct, &fread_ptr, pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, &ldbase_difflist_len);
  } else {
    reterr = ParseAndSaveDifflistProperSubset(fread_end, sample_include, sample_include_cumulative_popcounts, raw_sample_ct, &fread_ptr, pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, &ldbase_difflist_len, pgrp->workspace_raregeno_tmp_loadbuf);
  }
  if (unlikely(reterr)) {
    return reterr;
  }
  pgrp->ldbase_difflist_len = ldbase_difflist_len;
  pgrp->ldbase_difflist_sample_ids[ldbase_difflist_len] = sample_ct;
  pgrp->ldbase_stypes = kfPgrLdcacheDifflist;
  return kPglRetSuccess;
}

PglErr ReadDifflistOrGenovecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t max_simple_difflist_len, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict genovec, uint32_t* difflist_common_geno_ptr, uintptr_t* __restrict main_raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  assert(sample_ct);
  assert(max_simple_difflist_len < sample_ct);
  // Side effects:
  //   may use pgr.workspace_raregeno_tmp_loadbuf
  // Trailing bits of genovec/main_raregeno may not be zeroed out.
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t maintrack_vrtype = vrtype & 7;
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  // const uint32_t multiallelic_hc_present = fread_pp && VrtypeMultiallelic(vrtype);
  if (VrtypeLdCompressed(maintrack_vrtype)) {
    // LD compression

    // note that this can currently load a difflist longer than
    // max_simple_difflist_len
    PglErr reterr = LdLoadMinimalSubsetIfNecessary(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp);
    if (unlikely(reterr)) {
      return reterr;
    }
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (unlikely(InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end))) {
      return kPglRetReadFail;
    }
    const uint32_t ld_invert = (maintrack_vrtype == 3);
    if (pgrp->ldbase_stypes & kfPgrLdcacheDifflist) {
      const uint32_t ldbase_common_geno = pgrp->fi.vrtypes[pgrp->ldbase_vidx] & 3;
      // unnecessary for this to branch on LD difflist length, since that's
      // limited to 3/4 of the ldbase difflist length.
      *difflist_common_geno_ptr = ldbase_common_geno;
      reterr = ParseLdAndMergeDifflistSubset(fread_end, subsetting_required? sample_include : nullptr, sample_include_cumulative_popcounts, pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, pgrp->ldbase_difflist_len, ldbase_common_geno, raw_sample_ct, sample_ct, &fread_ptr, main_raregeno, difflist_sample_ids, difflist_len_ptr, pgrp->workspace_raregeno_tmp_loadbuf);
      if (unlikely(reterr)) {
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
    reterr = ParseAndApplyDifflistSubset(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, &fread_ptr, pgrp, genovec);
    if (unlikely(reterr)) {
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
  if (unlikely(InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end))) {
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
    PglErr reterr = ParseNonLdGenovecSubsetUnsafe(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, vrtype, &fread_ptr, pgrp, genovec);
    if (unlikely(reterr)) {
      return reterr;
    }
    const uint32_t ldbase_raw_genovec_saved = (subsetting_required && (!(vrtype & 4)));
    if (is_ldbase) {
      CopyQuaterarr(genovec, sample_ct, pgrp->ldbase_genovec);
      pgrp->ldbase_stypes = ldbase_raw_genovec_saved? (kfPgrLdcacheQuater | kfPgrLdcacheRawQuater) : kfPgrLdcacheQuater;
    } else if (ldbase_raw_genovec_saved) {
      // bugfix (22 Sep 2018)
      pgrp->ldbase_stypes &= ~kfPgrLdcacheRawQuater;
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
  PglErr reterr;
  if (!subsetting_required) {
    reterr = ParseAndSaveDifflist(fread_end, raw_sample_ct, &fread_ptr, main_raregeno, difflist_sample_ids, difflist_len_ptr);
  } else {
    reterr = ParseAndSaveDifflistProperSubset(fread_end, sample_include, sample_include_cumulative_popcounts, raw_sample_ct, &fread_ptr, main_raregeno, difflist_sample_ids, difflist_len_ptr, pgrp->workspace_raregeno_tmp_loadbuf);
  }
  if (unlikely(reterr)) {
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

PglErr PgrGetDifflistOrGenovec(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t max_simple_difflist_len, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uint32_t* difflist_common_geno_ptr, uintptr_t* __restrict main_raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    *difflist_common_geno_ptr = UINT32_MAX;
    return kPglRetSuccess;
  }
  return ReadDifflistOrGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, max_simple_difflist_len, vidx, pgrp, nullptr, nullptr, genovec, difflist_common_geno_ptr, main_raregeno, difflist_sample_ids, difflist_len_ptr);
}

PglErr LdSubsetAdjustGenocounts(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict ldbase_genovec, uint32_t raw_sample_ct, const unsigned char** fread_pp, STD_ARRAY_REF(uint32_t, 4) genocounts, uintptr_t* __restrict raregeno_workspace) {
  // * sample_include assumed to be nullptr if no subsetting required
  // * Assumes genocounts[] is initialized to the proper values for the LD
  //   reference variant (including subsetting).
  // * Tried a hybrid implementation which allowed the base variant to be saved
  //   as a difflist; turns out it's practically always better to unpack to a
  //   genovec first.
  // * There are two modes:
  //   1. If sample_include is nullptr, we're not selecting a sample subset.
  //   2. If sample_include and sample_include_cumulative_popcounts are both
  //      non-null, we're computing counts over a sample subset, and
  //      ldbase_genovec is assumed to be subsetted.
  //   Experimented with a third mode where ldbase_genovec was replaced with
  //   ldbase_raw_genovec in the subsetted case, but that didn't seem to pay
  //   off.
  // * This is the main frequency-counting bottleneck.
  uint32_t raw_difflist_len;
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno_workspace, &group_info_iter, &raw_difflist_len);
  if (reterr || (!raw_difflist_len)) {
    return reterr;
  }
  const uint32_t subgroup_idx_last = (raw_difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  uintptr_t* raregeno_workspace_iter = raregeno_workspace;
  uintptr_t raw_sample_idx = 0;
  STD_ARRAY_DECL(uint32_t, 16, delta_counts);
  STD_ARRAY_FILL0(delta_counts);
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
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
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    uintptr_t cur_raregeno_word = *raregeno_workspace_iter++;
    if (!sample_include) {
      for (; ; --remaining_deltas_in_subgroup) {
#ifndef __LP64__
        if (unlikely(raw_sample_idx >= raw_sample_ct)) {
          return kPglRetMalformedInput;
        }
#endif
        const uintptr_t cur_geno = cur_raregeno_word & 3;
        delta_counts[cur_geno * 4 + GetQuaterarrEntry(ldbase_genovec, raw_sample_idx)] += 1;
        if (!remaining_deltas_in_subgroup) {
          break;
        }
        raw_sample_idx += GetVint31(fread_end, fread_pp);
        cur_raregeno_word >>= 2;
      }
    } else {
      for (; ; --remaining_deltas_in_subgroup) {
#ifndef __LP64__
        if (unlikely(raw_sample_idx >= raw_sample_ct)) {
          return kPglRetMalformedInput;
        }
#endif
        if (IsSet(sample_include, raw_sample_idx)) {
          const uintptr_t cur_geno = cur_raregeno_word & 3;
          const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, raw_sample_idx);
          delta_counts[cur_geno * 4 + GetQuaterarrEntry(ldbase_genovec, sample_idx)] += 1;
        }
        if (!remaining_deltas_in_subgroup) {
          break;
        }
        raw_sample_idx += GetVint31(fread_end, fread_pp);
        cur_raregeno_word >>= 2;
      }
    }
  }
}

PglErr SkipDeltalistIds(const unsigned char* fread_end, const unsigned char* group_info, uint32_t difflist_len, uint32_t raw_sample_ct, uint32_t has_genotypes, const unsigned char** fread_pp) {
  assert(difflist_len);
  // fread_pp is a pure output parameter here
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  const unsigned char* extra_byte_cts = &(group_info[group_ct * sample_id_byte_ct]);
  const uint32_t extra_byte_tot = BytesumArr(extra_byte_cts, group_ct - 1);

  // (group_ct - 1) for extra_byte_cts
  // (difflist_len + 3) / 4 for raregeno
  // (group_ct - 1) * (kPglDifflistGroupSize - 1) + extra_byte_tot for
  //   all but last ID block
  // total = (group_ct - 1) * kPglDifflistGroupSize + extra_byte_tot +
  //         (difflist_len + 3) / 4
#ifdef __arm__
#  error "Unaligned accesses in SkipDeltalistIds()."
#endif
  const unsigned char* iddiff_start = &(extra_byte_cts[(group_ct - 1) * kPglDifflistGroupSize + extra_byte_tot]);
  if (has_genotypes) {
    iddiff_start = &(iddiff_start[QuaterCtToByteCt(difflist_len)]);
  }
  const uintptr_t* fread_alias = R_CAST(const uintptr_t*, iddiff_start);
  const uintptr_t* fread_alias_stop = R_CAST(const uintptr_t*, &(fread_end[-S_CAST(int32_t, kBytesPerWord)]));
  uint32_t remaining_id_ct = (difflist_len - 1) % kPglDifflistGroupSize;
#ifdef __LP64__
  while (remaining_id_ct >= kBytesPerVec) {
    if (unlikely(fread_alias > fread_alias_stop)) {
      return kPglRetMalformedInput;
    }
    const VecW vv = vecw_loadu(R_CAST(const VecW*, fread_alias));
    fread_alias = &(fread_alias[kWordsPerVec]);
    const uint32_t highbits = vecw_movemask(vv);
    remaining_id_ct -= kBytesPerVec - PopcountVec8thUint(highbits);
  }
#endif
  while (remaining_id_ct >= kBytesPerWord) {
    // scan a word at a time, count number of high bits set
    if (unlikely(fread_alias > fread_alias_stop)) {
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
  while (likely(fread_ptr < fread_end)) {
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

PglErr CountparseDifflistSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, uint32_t common_geno, uint32_t raw_sample_ct, uint32_t sample_ct, const unsigned char** fread_pp, STD_ARRAY_REF(uint32_t, 4) genocounts, uintptr_t* __restrict raregeno_workspace) {
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno_workspace, &group_info_iter, &difflist_len);
  STD_ARRAY_REF_FILL0(4, genocounts);
  if (reterr || (!difflist_len)) {
    genocounts[common_geno] = sample_ct;
    return reterr;
  }
  if (raw_sample_ct == sample_ct) {
    ZeroTrailingQuaters(difflist_len, raregeno_workspace);
    GenovecCountFreqsUnsafe(raregeno_workspace, difflist_len, genocounts);
    genocounts[common_geno] = sample_ct - difflist_len;
    // bugfix (26 Mar 2019): forgot to advance fread_pp
    return SkipDeltalistIds(fread_end, group_info_iter, difflist_len, raw_sample_ct, 1, fread_pp);
  }
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  uintptr_t* raregeno_workspace_iter = raregeno_workspace;
  uintptr_t raw_sample_idx = 0;
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        genocounts[common_geno] = sample_ct - genocounts[0] - genocounts[1] - genocounts[2] - genocounts[3];
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
#ifdef __LP64__
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    uintptr_t cur_raregeno_word = *raregeno_workspace_iter++;
    for (; ; --remaining_deltas_in_subgroup) {
#ifndef __LP64__
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
#endif
      if (IsSet(sample_include, raw_sample_idx)) {
        const uintptr_t cur_geno = cur_raregeno_word & 3;
        genocounts[cur_geno] += 1;
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      cur_raregeno_word >>= 2;
    }
  }
}

// 1-bit, unsubsetted: count 1-bit array, then count raregeno
// 1-bit, subsetted: count [1-bit array AND sample_include], iterate through
//   difflist
PglErr CountparseOnebitSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, const unsigned char** fread_pp, STD_ARRAY_REF(uint32_t, 4) genocounts, uintptr_t* __restrict raregeno_workspace) {
  const uint32_t initial_bitarray_byte_ct = DivUp(raw_sample_ct, CHAR_BIT);
  const unsigned char* onebit_main_iter = *fread_pp;
  if (PtrAddCk(fread_end, initial_bitarray_byte_ct + 1, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uint32_t common2_code = *onebit_main_iter++;
  const uint32_t geno_code_low = common2_code / 4;
  const uint32_t geno_code_high = (common2_code & 3) + geno_code_low;
#ifdef __arm__
#  error "Unaligned accesses in CountparseOnebitSubset()."
#endif
  uint32_t high_geno_ct;
  if (raw_sample_ct == sample_ct) {
    high_geno_ct = PopcountBytes(onebit_main_iter, initial_bitarray_byte_ct);
  } else {
    high_geno_ct = PopcountBytesMasked(onebit_main_iter, sample_include, initial_bitarray_byte_ct);
  }
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, raregeno_workspace, &group_info_iter, &difflist_len);
  STD_ARRAY_REF_FILL0(4, genocounts);
  if (reterr || (!difflist_len)) {
    genocounts[geno_code_low] = sample_ct - high_geno_ct;
    genocounts[geno_code_high] = high_geno_ct;
    return reterr;
  }
  if (raw_sample_ct == sample_ct) {
    ZeroTrailingQuaters(difflist_len, raregeno_workspace);
    GenovecCountFreqsUnsafe(raregeno_workspace, difflist_len, genocounts);
    genocounts[geno_code_low] = sample_ct - difflist_len - high_geno_ct;
    genocounts[geno_code_high] = high_geno_ct;
    // bugfix (26 Mar 2019): forgot to advance fread_pp
    return SkipDeltalistIds(fread_end, group_info_iter, difflist_len, raw_sample_ct, 1, fread_pp);
  }
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kBitsPerWordD2;
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  const uintptr_t* onebitarr = R_CAST(const uintptr_t*, onebit_main_iter);
  uintptr_t* raregeno_workspace_iter = raregeno_workspace;
  uintptr_t raw_sample_idx = 0;
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        // avoid read-after-write dependency?
        genocounts[geno_code_low] = sample_ct - high_geno_ct - genocounts[0] - genocounts[1] - genocounts[2] - genocounts[3];
        genocounts[geno_code_high] = high_geno_ct;
        return kPglRetSuccess;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
#ifdef __LP64__
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
#endif
      raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    uintptr_t cur_raregeno_word = *raregeno_workspace_iter++;
    for (; ; --remaining_deltas_in_subgroup) {
#ifndef __LP64__
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
#endif
      if (IsSet(sample_include, raw_sample_idx)) {
        const uintptr_t cur_geno = cur_raregeno_word & 3;
        genocounts[cur_geno] += 1;
        high_geno_ct -= IsSet(onebitarr, raw_sample_idx);
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
      cur_raregeno_word >>= 2;
    }
  }
}

// loads ldbase variant if necessary, guarantees pgrp->ldbase_genovec is filled
// on return
// only called by GetBasicGenotypeCounts(), usually LdLoadAndCopy... is better
PglErr LdLoadGenovecSubsetIfNecessary(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp) {
  if (LdLoadNecessary(vidx, pgrp)) {
    const uint32_t ldbase_vidx = pgrp->ldbase_vidx;
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (unlikely(InitReadPtrs(ldbase_vidx, pgrp, &fread_ptr, &fread_end))) {
      return kPglRetReadFail;
    }
    const uint32_t vrtype = pgrp->fi.vrtypes[ldbase_vidx];
    // bugfix (6 Mar 2019): ldbase_raw_genovec is only filled in (!difflist) &&
    //   subsetting_required case; (!difflist) isn't enough
    pgrp->ldbase_stypes = ((vrtype & 4) || (sample_ct == pgrp->fi.raw_sample_ct))? kfPgrLdcacheQuater : (kfPgrLdcacheQuater | kfPgrLdcacheRawQuater);
    return ParseNonLdGenovecSubsetUnsafe(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, vrtype, &fread_ptr, pgrp, pgrp->ldbase_genovec);
  }
  if (!(pgrp->ldbase_stypes & kfPgrLdcacheQuater)) {
    if (pgrp->ldbase_stypes & kfPgrLdcacheDifflist) {
      PgrDifflistToGenovecUnsafe(pgrp->ldbase_raregeno, pgrp->ldbase_difflist_sample_ids, pgrp->fi.vrtypes[pgrp->ldbase_vidx] & 3, sample_ct, pgrp->ldbase_difflist_len, pgrp->ldbase_genovec);
    } else {
      assert(pgrp->ldbase_stypes & kfPgrLdcacheRawQuater);
      CopyQuaterarrNonemptySubset(pgrp->ldbase_raw_genovec, sample_include, pgrp->fi.raw_sample_ct, sample_ct, pgrp->ldbase_genovec);
    }
    pgrp->ldbase_stypes |= kfPgrLdcacheQuater;
  }
  return kPglRetSuccess;
}

PglErr GetBasicGenotypeCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uint32_t* unphased_het_ctp, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // genocounts[0] := ref/ref, genocounts[1] := ref/altx,
  // genocounts[2] := altx/alty, genocounts[3] := missing
  // If unphased_het_ctp is non-null, this assumes multiallelic hardcalls are
  // not present, phased hardcalls are present, we aren't subsetting, and
  // unphased_het_ct is initialized to zero.
  assert(vidx < pgrp->fi.raw_variant_ct);
  assert(sample_ct);
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  const unsigned char* fread_ptr;
  const unsigned char* fread_end = nullptr;  // maybe-uninitialized warning
  PglErr reterr;
  if (VrtypeLdCompressed(vrtype)) {
    // LD compression
    reterr = LdLoadGenovecSubsetIfNecessary(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (unlikely(InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end))) {
      return kPglRetReadFail;
    }
    if (!(pgrp->ldbase_stypes & kfPgrLdcacheBasicGenocounts)) {
      ZeroTrailingQuaters(sample_ct, pgrp->ldbase_genovec);
      GenovecCountFreqsUnsafe(pgrp->ldbase_genovec, sample_ct, pgrp->ldbase_basic_genocounts);
      pgrp->ldbase_stypes |= kfPgrLdcacheBasicGenocounts;
    }
    STD_ARRAY_COPY(pgrp->ldbase_basic_genocounts, 4, genocounts);
    reterr = LdSubsetAdjustGenocounts(fread_end, subsetting_required? sample_include : nullptr, sample_include_cumulative_popcounts, pgrp->ldbase_genovec, raw_sample_ct, &fread_ptr, genocounts, pgrp->workspace_raregeno_tmp_loadbuf);
    if (vrtype & 1) {
      // inverted
      const uint32_t tmpval = genocounts[0];
      genocounts[0] = genocounts[2];
      genocounts[2] = tmpval;
    }
  } else {
    if (unlikely(InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end))) {
      return kPglRetReadFail;
    }
    const uint32_t is_ldbase = pgrp->fi.vrtypes && VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
    if (is_ldbase) {
      // difflists are very efficient to count directly when not subsetting
      // (since we can entirely ignore the sample IDs), but it's often better
      // to unpack them first when subsetting.

      // ...er, the statement above is a lie, unpack-first almost always seems
      // to be better.
      pgrp->ldbase_vidx = vidx;
      // this may be slowed down by the LD caching change.
      reterr = ParseNonLdGenovecSubsetUnsafe(fread_end, sample_include, sample_include_cumulative_popcounts, sample_ct, vrtype, &fread_ptr, pgrp, pgrp->ldbase_genovec);
      ZeroTrailingQuaters(sample_ct, pgrp->ldbase_genovec);
      GenovecCountFreqsUnsafe(pgrp->ldbase_genovec, sample_ct, genocounts);
      STD_ARRAY_COPY(genocounts, 4, pgrp->ldbase_basic_genocounts);
      pgrp->ldbase_stypes = (subsetting_required && (!(vrtype & 4)))? (kfPgrLdcacheQuater | kfPgrLdcacheRawQuater | kfPgrLdcacheBasicGenocounts) : (kfPgrLdcacheQuater | kfPgrLdcacheBasicGenocounts);
    } else if (vrtype & 4) {
      const uint32_t vrtype_low2 = vrtype & 3;
      if (vrtype_low2 != 1) {
        reterr = CountparseDifflistSubset(fread_end, sample_include, vrtype & 3, raw_sample_ct, sample_ct, &fread_ptr, genocounts, pgrp->workspace_raregeno_tmp_loadbuf);
      } else {
        genocounts[0] = sample_ct;
        genocounts[1] = 0;
        genocounts[2] = 0;
        genocounts[3] = 0;
        reterr = kPglRetSuccess;
      }
    } else if (vrtype & 1) {
      reterr = CountparseOnebitSubset(fread_end, sample_include, raw_sample_ct, sample_ct, &fread_ptr, genocounts, pgrp->workspace_raregeno_tmp_loadbuf);
    } else {
      const uint32_t genovec_byte_ct = QuaterCtToByteCt(raw_sample_ct);
      const unsigned char* raw_genovec = fread_ptr;
      if (PtrAddCk(fread_end, genovec_byte_ct, &fread_ptr)) {
        return kPglRetMalformedInput;
      }
      const uint32_t raw_genovec_unaligned = R_CAST(uintptr_t, raw_genovec) % kBytesPerVec;
      if (!subsetting_required) {
        if (raw_genovec_unaligned) {
          GenoarrCountFreqs(raw_genovec, raw_sample_ct, genocounts);
        } else {
          GenovecCountFreqs(R_CAST(const uintptr_t*, raw_genovec), raw_sample_ct, genocounts);
        }
      } else {
        GenoarrCountSubsetFreqs(raw_genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
      }
      if (vrtype == kPglVrtypePlink1) {
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
      reterr = kPglRetSuccess;
    }
  }
  if ((!unphased_het_ctp) || reterr) {
    return reterr;
  }
  assert((!subsetting_required) && ((vrtype & 0x18) == 0x10));
  const uint32_t het_ct = genocounts[1];
  const uint32_t aux2_first_part_byte_ct = 1 + (het_ct / CHAR_BIT);
  if (PtrCheck(fread_end, fread_ptr, aux2_first_part_byte_ct)) {
    return kPglRetMalformedInput;
  }
  const uint32_t explicit_phasepresent = fread_ptr[0] & 1;
  if (explicit_phasepresent) {
    // otherwise initial value if 0 is correct
    *unphased_het_ctp = het_ct + 1 - PopcountBytes(fread_ptr, aux2_first_part_byte_ct);
  }
  return kPglRetSuccess;
}

PglErr PgrGetCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    STD_ARRAY_REF_FILL0(4, genocounts);
    return kPglRetSuccess;
  }
  return GetBasicGenotypeCounts(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, genocounts);
}

// Ok for quater_vvec to be unaligned.
uint32_t CountQuaterVec6(const VecW* quater_vvec, uintptr_t quater_word, uint32_t vec_ct) {
  assert(!(vec_ct % 6));
  const VecW m0 = vecw_setzero();
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW xor_vvec = vecw_set1(quater_word);
  const VecW* quater_vvec_iter = quater_vvec;
  VecW prev_sad_result = vecw_setzero();
  VecW acc = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        acc = acc + prev_sad_result;
        return HsumW(acc);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc = vecw_setzero();
    const VecW* quater_vvec_stop = &(quater_vvec_iter[cur_incr]);
    do {
      VecW loader1 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      VecW loader2 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      VecW count1 = (~(vecw_srli(loader1, 1) | loader1)) & m1;
      VecW count2 = (~(vecw_srli(loader2, 1) | loader2)) & m1;

      loader1 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      loader2 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      count1 = count1 + ((~(vecw_srli(loader1, 1) | loader1)) & m1);
      count2 = count2 + ((~(vecw_srli(loader2, 1) | loader2)) & m1);

      loader1 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      loader2 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      count1 = count1 + ((~(vecw_srli(loader1, 1) | loader1)) & m1);
      count2 = count2 + ((~(vecw_srli(loader2, 1) | loader2)) & m1);

      count1 = (count1 & m2) + (vecw_srli(count1, 2) & m2);
      count1 = count1 + (count2 & m2) + (vecw_srli(count2, 2) & m2);
      inner_acc = inner_acc + (count1 & m4) + (vecw_srli(count1, 4) & m4);
    } while (quater_vvec_iter < quater_vvec_stop);
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

// Ok for quatervec to be unaligned.  Ok if unsafe to read trailing bytes of
// quatervec.
uint32_t CountQuater(const void* quaterarr, uintptr_t quater_word, uint32_t quater_ct) {
  const uint32_t fullword_ct = quater_ct / kBitsPerWordD2;
  uint32_t word_idx = fullword_ct - (fullword_ct % (6 * kWordsPerVec));
  uint32_t tot = CountQuaterVec6(S_CAST(const VecW*, quaterarr), quater_word, word_idx / kWordsPerVec);
  const uintptr_t* quatervec = S_CAST(const uintptr_t*, quaterarr);
  for (; word_idx != fullword_ct; ++word_idx) {
    const uintptr_t cur_word = quatervec[word_idx] ^ quater_word;
    tot += Popcount01Word(Word00(cur_word));
  }
  const uint32_t trailing_quater_ct = quater_ct % kBitsPerWordD2;
  if (trailing_quater_ct) {
    const uint32_t trailing_byte_ct = DivUp(trailing_quater_ct, (CHAR_BIT / 2));
    uintptr_t cur_word = SubwordLoad(&(quatervec[fullword_ct]), trailing_byte_ct) ^ quater_word;
    cur_word = bzhi(Word00(cur_word), trailing_quater_ct * 2);
    tot += Popcount01Word(cur_word);
  }
  return tot;
}

/*
uint32_t CountQuaterSubsetVec6(const VecW* __restrict quater_vvec, const VecW* __restrict interleaved_mask_vvec, uintptr_t quater_word, uint32_t vec_ct) {
  assert(!(vec_ct % 6));
  const VecW m0 = vecw_setzero();
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW xor_vvec = vecw_set1(quater_word);
  const VecW* quater_vvec_iter = quater_vvec;
  const VecW* interleaved_mask_vvec_iter = interleaved_mask_vvec;
  VecW prev_sad_result = vecw_setzero();
  VecW acc = vecw_setzero();
  uintptr_t cur_incr = 60;
  while (1) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        acc = acc + prev_sad_result;
        return HsumW(acc);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc = vecw_setzero();
    const VecW* quater_vvec_stop = &(quater_vvec_iter[cur_incr]);
    vec_ct -= cur_incr;
    do {
      VecW mask1 = *interleaved_mask_vvec_iter++;
      VecW loader1 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      VecW mask2 = vecw_srli(mask1, 1) & m1;
      VecW loader2 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      mask1 = mask1 & m1;
      VecW count1 = (~(vecw_srli(loader1, 1) | loader1)) & mask1;
      VecW count2 = (~(vecw_srli(loader2, 1) | loader2)) & mask2;

      mask1 = *interleaved_mask_vvec_iter++;
      loader1 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      mask2 = vecw_srli(mask1, 1) & m1;
      loader2 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      mask1 = mask1 & m1;
      count1 = count1 + ((~(vecw_srli(loader1, 1) | loader1)) & mask1);
      count2 = count2 + ((~(vecw_srli(loader2, 1) | loader2)) & mask2);

      mask1 = *interleaved_mask_vvec_iter++;
      loader1 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      mask2 = vecw_srli(mask2, 1) & m1;
      loader2 = vecw_loadu(quater_vvec_iter++) ^ xor_vvec;
      mask1 = mask1 & m1;
      count1 = count1 + ((~(vecw_srli(loader1, 1) | loader1)) & mask1);
      count2 = count2 + ((~(vecw_srli(loader2, 1) | loader2)) & mask2);

      count1 = (count1 & m2) + (vecw_srli(count1, 2) & m2);
      count1 = count1 + (count2 & m2) + (vecw_srli(count2, 2) & m2);
      inner_acc = inner_acc + (count1 & m4) + (vecw_srli(count1, 4) & m4);
    } while (quater_vvec_iter < quater_vvec_stop);
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

uint32_t CountQuaterSubset(const uintptr_t* __restrict quatervec, const uintptr_t* __restrict interleaved_vec, uintptr_t quater_word, uint32_t raw_quater_ct) {
  // simplified GenovecCountSubsetFreqs()
  const uint32_t raw_quater_ctv2 = QuaterCtToVecCt(raw_quater_ct);
#ifdef __LP64__
  uint32_t vec_idx = raw_quater_ctv2 - (raw_quater_ctv2 % 6);
  uint32_t tot = CountQuaterSubsetVec6(R_CAST(const VecW*, quatervec), R_CAST(const VecW*, interleaved_vec), quater_word, vec_idx);
  const uintptr_t* quatervec_iter = &(quatervec[kWordsPerVec * vec_idx]);
  const uintptr_t* interleaved_mask_iter = &(interleaved_vec[(kWordsPerVec / 2) * vec_idx]);
#  ifdef USE_AVX2
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  uintptr_t mask_base3 = 0;
  uintptr_t mask_base4 = 0;
  for (; vec_idx != raw_quater_ctv2; ++vec_idx) {
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
      const uintptr_t cur_geno_word1 = (*quatervec_iter++) ^ quater_word;
      const uintptr_t cur_geno_word2 = (*quatervec_iter++) ^ quater_word;
      const uintptr_t masked1 = mask_word1 & (~(cur_geno_word1 | (cur_geno_word1 >> 1)));
      const uintptr_t masked2 = mask_word2 & (~(cur_geno_word2 | (cur_geno_word2 >> 1)));
      tot += PopcountWord((masked1 << 1) | masked2);
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
  for (; vec_idx != raw_quater_ctv2; ++vec_idx) {
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
    const uintptr_t cur_geno_word1 = (*quatervec_iter++) ^ quater_word;
    const uintptr_t cur_geno_word2 = (*quatervec_iter++) ^ quater_word;
    const uintptr_t masked1 = mask_word1 & (~(cur_geno_word1 | (cur_geno_word1 >> 1)));
    const uintptr_t masked2 = mask_word2 & (~(cur_geno_word2 | (cur_geno_word2 >> 1)));
#    ifdef USE_SSE42
    tot += PopcountWord((masked1 << 1) | masked2);
#    else
    tot += QuatersumWord(masked1 + masked2);
#    endif
  }
#  endif  // not USE_AVX2
#else  // not __LP64__
  uint32_t word_idx = raw_quater_ctv2 - (raw_quater_ctv2 % 6);
  uint32_t tot = CountQuaterSubsetVec6(R_CAST(const VecW*, quatervec), R_CAST(const VecW*, interleaved_vec), quater_word, word_idx);
  const uintptr_t* interleaved_mask_iter = &(interleaved_vec[word_idx / 2]);
  uintptr_t mask_base = 0;
  for (; word_idx != raw_quater_ctv2; ++word_idx) {
    uintptr_t mask_word;
    if (!(word_idx % 2)) {
      mask_base = *interleaved_mask_iter++;
      mask_word = mask_base & kMask5555;
    } else {
      mask_word = (mask_base >> 1) & kMask5555;
    }
    const uintptr_t cur_geno_word = quatervec[word_idx] ^ quater_word;
    const uintptr_t masked = mask_word & (~(cur_geno_word | (cur_geno_word >> 1)));
    tot += Popcount01Word(masked);
  }
#endif
  return tot;
}
*/

// Ok for nibble_vvec to be unaligned.
uint32_t CountNibbleVec(const VecW* nibble_vvec, uintptr_t nibble_word, uint32_t vec_ct) {
  const VecW m0 = vecw_setzero();
  const VecW alld15 = VCONST_W(kMask1111);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW xor_vvec = vecw_set1(nibble_word);
  const VecW* nibble_vvec_iter = nibble_vvec;
  VecW prev_sad_result = vecw_setzero();
  VecW acc = vecw_setzero();
  uintptr_t cur_incr = 15;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 15) {
      if (!vec_ct) {
        acc = acc + prev_sad_result;
        return HsumW(acc);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc = vecw_setzero();
    const VecW* nibble_vvec_stop = &(nibble_vvec_iter[cur_incr]);
    do {
      VecW loader = vecw_loadu(nibble_vvec_iter++) ^ xor_vvec;
      // DetectAllZeroNibbles() followed by right-shift-3 is the same number of
      // operations, can see if that's any faster in practice
      loader = vecw_srli(loader, 1) | loader;
      loader = vecw_srli(loader, 2) | loader;
      inner_acc = inner_acc + ((~loader) & alld15);
    } while (nibble_vvec_iter < nibble_vvec_stop);
    inner_acc = (inner_acc & m4) + (vecw_srli(inner_acc, 4) & m4);
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

uint32_t CountNibble(const void* nibblearr, uintptr_t nibble_word, uintptr_t nibble_ct) {
  const uint32_t fullword_ct = nibble_ct / kBitsPerWordD4;
  uint32_t tot = CountNibbleVec(S_CAST(const VecW*, nibblearr), nibble_word, fullword_ct / kWordsPerVec);
  const uintptr_t* nibblevec = S_CAST(const uintptr_t*, nibblearr);
#ifdef __LP64__
  for (uint32_t word_idx = RoundDownPow2(fullword_ct, kWordsPerVec); word_idx != fullword_ct; ++word_idx) {
    uintptr_t cur_word = nibblevec[word_idx] ^ nibble_word;
    cur_word = cur_word | (cur_word >> 1);
    cur_word = cur_word | (cur_word >> 2);
    tot += Popcount0001Word((~cur_word) & kMask1111);
  }
#endif
  const uint32_t trailing_nibble_ct = nibble_ct % kBitsPerWordD4;
  if (trailing_nibble_ct) {
    const uint32_t trailing_byte_ct = DivUp(trailing_nibble_ct, (CHAR_BIT / 4));
    uintptr_t cur_word = SubwordLoad(&(nibblevec[fullword_ct]), trailing_byte_ct) ^ nibble_word;
    cur_word = cur_word | (cur_word >> 1);
    cur_word = cur_word | (cur_word >> 2);
    cur_word = bzhi((~cur_word) & kMask1111, trailing_nibble_ct * 4);
#if defined(USE_SSE42) || !defined(__LP64__)
    tot += Popcount0001Word(cur_word);
#else
    // minor optimization, can't overflow
    tot += (cur_word * kMask1111) >> 60;
#endif
  }
  return tot;
}

// similar to ParseAndSaveDifflist()
PglErr ParseAndSaveDeltalist(const unsigned char* fread_end, uint32_t raw_sample_ct, const unsigned char** fread_pp, uint32_t* __restrict deltalist, uint32_t* __restrict deltalist_len_ptr) {
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, deltalist_len_ptr);
  const uint32_t deltalist_len = *deltalist_len_ptr;
  if (reterr || (!deltalist_len)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  const uint32_t group_idx_last = (deltalist_len - 1) / kPglDifflistGroupSize;
  uint32_t* deltalist_iter = deltalist;
  uint32_t group_len_m1 = kPglDifflistGroupSize - 1;
  for (uint32_t group_idx = 0; ; ++group_idx) {
    if (group_idx >= group_idx_last) {
      if (group_idx > group_idx_last) {
        return kPglRetSuccess;
      }
      group_len_m1 &= deltalist_len - 1;
    }
    uintptr_t raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
    group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    for (uint32_t raw_deltalist_idx_lowbits = 0; ; ++raw_deltalist_idx_lowbits) {
      // always check, otherwise we may scribble over arbitrary memory
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
      deltalist_iter[raw_deltalist_idx_lowbits] = raw_sample_idx;
      if (raw_deltalist_idx_lowbits == group_len_m1) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
    deltalist_iter = &(deltalist_iter[group_len_m1 + 1]);
  }
}

PglErr CountDeltalistIntersect(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, const unsigned char** fread_pp, uint32_t* __restrict intersect_ctp, uint32_t* __restrict raw_deltalist_len_ptr) {
  // Requires a PROPER subset.
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, raw_deltalist_len_ptr);
  const uint32_t raw_deltalist_len = *raw_deltalist_len_ptr;
  if (reterr || (!raw_deltalist_len)) {
    *intersect_ctp = 0;
    return reterr;
  }
  const uint32_t group_idx_last = (raw_deltalist_len - 1) / kPglDifflistGroupSize;
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  uintptr_t intersect_ct = 0;

  // technically doesn't need to be initialized, but I have principles
  uintptr_t raw_sample_idx = 0;

  uint32_t group_len_m1 = kPglDifflistGroupSize - 1;
  for (uint32_t group_idx = 0; ; ++group_idx) {
    if (group_idx >= group_idx_last) {
      if (group_idx > group_idx_last) {
        *intersect_ctp = intersect_ct;
        return kPglRetSuccess;
      }
      group_len_m1 &= raw_deltalist_len - 1;
    }
    // We need to pull a raw sample index from the deltalist header every 64
    // entries.
#ifdef __LP64__
    if (unlikely(raw_sample_idx >= raw_sample_ct)) {
      return kPglRetMalformedInput;
    }
#endif
    raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
    group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    for (uint32_t raw_deltalist_idx_lowbits = 0; ; ++raw_deltalist_idx_lowbits) {
#ifndef __LP64__
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
#endif
      intersect_ct += IsSet(sample_include, raw_sample_idx);
      if (raw_deltalist_idx_lowbits == group_len_m1) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
  }
}

uint32_t CountAux1aDense(const void* patch_01_fvals, uint32_t allele_ct, uint32_t allele_idx, uint32_t raw_01_ct, uint32_t rare01_ct) {
  // The 'f' in patch_01_fset/patch_01_fvals is to distinguish the in-file
  // representation from the returned AlleleCode*-based representation.
  if (allele_idx == 1) {
    // safe to ignore allele codes
    return raw_01_ct - rare01_ct;
  }
  if (allele_ct < 5) {
    if (allele_ct == 3) {
      return rare01_ct;
    }
    // need to count matches
    const uint32_t allele_code_byte_ct = DivUp(rare01_ct, 8);
    const uint32_t alt3_ct = PopcountBytes(patch_01_fvals, allele_code_byte_ct);
    if (allele_idx == 3) {
      return alt3_ct;
    }
    return rare01_ct - alt3_ct;
  }
  if (allele_ct < 19) {
    if (allele_ct < 7) {
      return CountQuater(patch_01_fvals, (allele_idx - 2) * kMask5555, rare01_ct);
    }
    return CountNibble(patch_01_fvals, (allele_idx - 2) * kMask1111, rare01_ct);
  }
  return CountByte(patch_01_fvals, allele_idx - 2, rare01_ct);
}

uint32_t GetAux1aWidth(uint32_t allele_ct) {
  if (allele_ct < 7) {
    if (allele_ct < 5) {
      return allele_ct - 3;
    }
    return 2;
  }
  if (allele_ct < 19) {
    return 4;
  }
  return 8;
}

// Returns allele_code_width.  Other return values are inaccurate for allele_ct
// == 3, since it's assumed that they're unused in that case.
uint32_t GetAux1aConsts(uint32_t allele_ct, uintptr_t* detect_mask_hi_ptr, uintptr_t* detect_mask_lo_ptr, uint32_t* allele_code_logwidth_ptr) {
  if (allele_ct < 7) {
    if (allele_ct < 5) {
      *detect_mask_hi_ptr = ~k0LU;
      *detect_mask_lo_ptr = ~k0LU;
      *allele_code_logwidth_ptr = 0;
      return allele_ct - 3;
    }
    *detect_mask_hi_ptr = kMaskAAAA;
    *detect_mask_lo_ptr = kMask5555;
    *allele_code_logwidth_ptr = 1;
    return 2;
  }
  if (allele_ct < 19) {
    *detect_mask_hi_ptr = kMask1111 * 8;
    *detect_mask_lo_ptr = kMask1111;
    *allele_code_logwidth_ptr = 2;
    return 4;
  }
  *detect_mask_hi_ptr = kMask0101 * 0x80;
  *detect_mask_lo_ptr = kMask0101;
  *allele_code_logwidth_ptr = 3;
  return 8;
}

// Advances *fread_pp past aux1a, and sets *het_ctp to the number of ref-altx
// hets where x == allele_idx in sample_include.  (If allele_idx == 1, *het_ctp
// is raw_01_ct - [# of aux1a entries] when there's no subsetting.)
// Note that raw_01_ct must be an un-subsetted count.
// Ok for subsetted_01_ct to be uninitialized if not subsetting, or allele_idx
// != 1.
// sample_include assumed to be nullptr if no subsetting required
PglErr CountAux1a(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uintptr_t* __restrict raw_genovec, uint32_t aux1a_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t allele_idx, uint32_t raw_01_ct, uint32_t subsetted_01_ct, const unsigned char** fread_pp, uint32_t* __restrict het_ctp, uint32_t* __restrict deltalist_workspace) {
  if (aux1a_mode == 15) {
    if (allele_idx == 1) {
      if (sample_include) {
        *het_ctp = subsetted_01_ct;
      } else {
        *het_ctp = raw_01_ct;
      }
    } else {
      *het_ctp = 0;
    }
    return kPglRetSuccess;
  }
  const uint32_t ignore_01_fvals = (allele_idx == 1) || (allele_ct == 3);
  uintptr_t detect_mask_hi;
  uintptr_t detect_mask_lo;
  uint32_t allele_code_logwidth;
  const uint32_t allele_code_width = GetAux1aConsts(allele_ct, &detect_mask_hi, &detect_mask_lo, &allele_code_logwidth);
  const uintptr_t xor_word = (allele_idx - 2) * detect_mask_lo;
  if (!aux1a_mode) {
    // 01-collapsed bitarray
    const uint32_t fset_byte_ct = DivUp(raw_01_ct, CHAR_BIT);
    const uint32_t rare01_ct = PopcountBytes(*fread_pp, fset_byte_ct);
#ifdef __arm__
#  error "Unaligned accesses in CountAux1a()."
#endif
    const uintptr_t* patch_01_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    *fread_pp += fset_byte_ct;
    const uintptr_t* patch_01_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    if (!sample_include) {
      *het_ctp = CountAux1aDense(patch_01_fvalsw, allele_ct, allele_idx, raw_01_ct, rare01_ct);
      return kPglRetSuccess;
    }
    const Halfword* sample_include_hw = R_CAST(const Halfword*, sample_include);
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_hets = Word01(raw_genovec[0]);
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t subsetted_hetx_ct = 0;
    uint32_t loop_len = kBitsPerWord;
    uint32_t rare01_lowbits = kBitsPerWord;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          break;
        }
        fset_bits = SubwordLoad(&(patch_01_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_01_ct, kBitsPerWord);
      } else {
        fset_bits = patch_01_fsetw[fset_widx];
      }
      // format 0, sample_include non-null
      if (ignore_01_fvals) {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_hets) {
            cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            // Considered replacing cur_raw_genovec_hets with the result of
            // two PackWordToHalfword() operations, since that keeps all
            // the sample word-indexes aligned.  Couldn't justify it given
            // the expected sparsity of this case, though.
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_hets) / 2;
            subsetted_hetx_ct += (sample_include_hw[sample_hwidx] >> sample_uidx_lowbits) & 1;
          }
          cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
          fset_bits = fset_bits >> 1;
        }
      } else {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_hets) {
            cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare01_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_01_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_01_fvalsw[fvals_widx];
              }
              fvals_bits = fvals_bits ^ xor_word;
              fvals_bits = (detect_mask_hi & (~(fvals_bits | ((fvals_bits | detect_mask_hi) - detect_mask_lo)))) >> (allele_code_width - 1);
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare01_lowbits = 0;
            }
            if (fvals_bits & (k1LU << rare01_lowbits)) {
              const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_hets) / 2;
              subsetted_hetx_ct += (sample_include_hw[sample_hwidx] >> sample_uidx_lowbits) & 1;
            }
            rare01_lowbits += allele_code_width;
          }
          cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
          fset_bits = fset_bits >> 1;
        }
      }
    }
    if (allele_idx == 1) {
      *het_ctp = subsetted_01_ct - subsetted_hetx_ct;
    } else {
      *het_ctp = subsetted_hetx_ct;
    }
    return kPglRetSuccess;
  }
  // mode 1: difflist.
  if (!sample_include) {
    const unsigned char* group_info_iter;
    uint32_t rare01_ct;
    PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, &rare01_ct);
    // rare01_ct == 0 should be impossible
    if (unlikely(reterr)) {
      return reterr;
    }
    reterr = SkipDeltalistIds(fread_end, group_info_iter, rare01_ct, raw_sample_ct, 1, fread_pp);
    if (unlikely(reterr)) {
      return reterr;
    }
    const unsigned char* patch_01_fvals = *fread_pp;
    const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }

    *het_ctp = CountAux1aDense(patch_01_fvals, allele_ct, allele_idx, raw_01_ct, rare01_ct);
    return kPglRetSuccess;
  }
  if (ignore_01_fvals) {
    // Don't need to save deltalist contents in this case.
    uint32_t subsetted_hetx_ct;
    uint32_t rare01_ct;
    PglErr reterr = CountDeltalistIntersect(fread_end, sample_include, raw_sample_ct, fread_pp, &subsetted_hetx_ct, &rare01_ct);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (allele_idx == 1) {
      *het_ctp = subsetted_01_ct - subsetted_hetx_ct;
      const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
      if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
        return kPglRetMalformedInput;
      }
    } else {
      *het_ctp = subsetted_hetx_ct;
    }
    return kPglRetSuccess;
  }
  // Save deltalist elements, iterate.
  uint32_t rare01_ct;
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare01_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_01_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  uint32_t subsetted_hetx_ct = 0;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        break;
      }
      fvals_bits = SubwordLoad(&(patch_01_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
    } else {
      fvals_bits = patch_01_fvalsw[fvals_widx];
    }
    fvals_bits = fvals_bits ^ xor_word;
    fvals_bits = detect_mask_hi & (~(fvals_bits | ((fvals_bits | detect_mask_hi) - detect_mask_lo)));
    if (fvals_widx == fvals_word_ct_m1) {
      fvals_bits = bzhi_max(fvals_bits, ModNz(rare01_ct << allele_code_logwidth, kBitsPerWord));
    }
    if (!fvals_bits) {
      continue;
    }
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - allele_code_logwidth)]);
    do {
      const uint32_t rare01_idx_lowbits = ctzw(fvals_bits) >> allele_code_logwidth;
      const uint32_t sample_uidx = cur_deltalist_base[rare01_idx_lowbits];
      subsetted_hetx_ct += IsSet(sample_include, sample_uidx);
      fvals_bits &= fvals_bits - 1;
    } while (fvals_bits);
  }
  *het_ctp = subsetted_hetx_ct;
  return kPglRetSuccess;
}

void CountAux1bDense(const void* patch_10_fvals, uint32_t allele_ct, uint32_t allele_idx_m1, uint32_t raw_10_ct, uint32_t rare10_ct, uint32_t* __restrict het_ctp, uint32_t* __restrict hom_ctp) {
  uint32_t matching_hom_ct = 0;
  uint32_t het_incr;
  if (allele_ct < 6) {
    if (allele_ct == 3) {
      const uint32_t allele_code_byte_ct = DivUp(rare10_ct, 8);
      matching_hom_ct = PopcountBytes(patch_10_fvals, allele_code_byte_ct);
      het_incr = rare10_ct - matching_hom_ct;
    } else {
      // 2+2 bits
      het_incr = CountQuater(patch_10_fvals, allele_idx_m1 * kMask5555, rare10_ct * 2);
      if (allele_idx_m1) {
        matching_hom_ct = CountNibble(patch_10_fvals, allele_idx_m1 * kMask5555, rare10_ct);
      }
    }
  } else {
    if (allele_ct < 18) {
      // 4+4 bits
      het_incr = CountNibble(patch_10_fvals, allele_idx_m1 * kMask1111, rare10_ct * 2);
      if (allele_idx_m1) {
        matching_hom_ct = CountByte(patch_10_fvals, allele_idx_m1 * 0x11, rare10_ct);
      }
    } else {
      // 8+8 bits
      het_incr = CountByte(patch_10_fvals, allele_idx_m1 * 0x11, rare10_ct * 2);
      if (allele_idx_m1) {
        matching_hom_ct = CountU16(patch_10_fvals, allele_idx_m1 * 0x1111, rare10_ct);
      }
    }
  }
  if (!allele_idx_m1) {
    *hom_ctp = raw_10_ct - rare10_ct;
  } else {
    het_incr -= 2 * matching_hom_ct;
    *hom_ctp = matching_hom_ct;
  }
  *het_ctp += het_incr;
}

// Returns allele_code_logwidth.
uint32_t GetAux1bConsts(uint32_t allele_ct, uintptr_t* detect_hom_mask_lo_ptr) {
  if (allele_ct < 6) {
    if (allele_ct == 3) {
      *detect_hom_mask_lo_ptr = ~k0LU;
      return 0;
    }
    *detect_hom_mask_lo_ptr = kMask1111;
    return 1;
  }
  if (allele_ct < 18) {
    *detect_hom_mask_lo_ptr = kMask0101;
    return 2;
  }
  *detect_hom_mask_lo_ptr = kMask0001;
  return 3;
}

// Advances *fread_pp past aux1b; increments *het_ctp by the number of
// altx-alty genotypes in aux1b and sample_include with one allele ==
// allele_idx; and sets *hom_ctp to the number of such hom-allele_idx genotypes
// present.  (For allele_idx == 1, *hom_ctp is equal to raw_10_ct -
// [# of aux1b entries] when there's no subsetting.)
// Trailing bits of raw_genovec must be cleared.
// Ok for subsetted_10_ct to be uninitialized if not subsetting, or allele_idx
// != 1.
// sample_include assumed to be nullptr if no subsetting required
PglErr CountAux1b(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uintptr_t* __restrict raw_genovec, uint32_t aux1b_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t allele_idx, uint32_t raw_10_ct, uint32_t subsetted_10_ct, const unsigned char** fread_pp, uint32_t* __restrict het_ctp, uint32_t* __restrict hom_ctp, uint32_t* __restrict deltalist_workspace) {
  if (aux1b_mode == 15) {
    if (allele_idx == 1) {
      if (sample_include) {
        *hom_ctp = subsetted_10_ct;
      } else {
        *hom_ctp = raw_10_ct;
      }
    } else {
      *hom_ctp = 0;
    }
    return kPglRetSuccess;
  }
  uintptr_t detect_hom_mask_lo;
  const uint32_t allele_code_logwidth = GetAux1bConsts(allele_ct, &detect_hom_mask_lo);
  const uint32_t code10_logwidth = allele_code_logwidth + (allele_code_logwidth != 0);
  const uint32_t code10_width = 1U << code10_logwidth;
  const uint32_t allele_idx_m1 = allele_idx - 1;
  uint32_t rare10_lowbits = kBitsPerWord;
  if (!aux1b_mode) {
    // 10-collapsed bitarray
    const uint32_t fset_byte_ct = DivUp(raw_10_ct, CHAR_BIT);
    const uint32_t rare10_ct = PopcountBytes(*fread_pp, fset_byte_ct);
#ifdef __arm__
#  error "Unaligned accesses in CountAux1b()."
#endif
    const uintptr_t* patch_10_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    *fread_pp += fset_byte_ct;
    const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) * code10_width, 8);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    if (!sample_include) {
      CountAux1bDense(patch_10_fvalsw, allele_ct, allele_idx_m1, raw_10_ct, rare10_ct, het_ctp, hom_ctp);
      return kPglRetSuccess;
    }
    const Halfword* sample_include_hw = R_CAST(const Halfword*, sample_include);
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_xys = Word10(raw_genovec[0]);
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t loop_len = kBitsPerWord;
    if ((!allele_idx_m1) || (allele_ct == 3)) {
      const uintptr_t detect_hom_mask_hi = detect_hom_mask_lo << (code10_width - 1);
      uint32_t subsetted_rare10_ct = 0;
      uint32_t het_incr = 0;
      for (uint32_t fset_widx = 0; ; ++fset_widx) {
        uintptr_t fset_bits;
        if (fset_widx >= fset_word_ct_m1) {
          if (fset_widx > fset_word_ct_m1) {
            break;
          }
          fset_bits = SubwordLoad(&(patch_10_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
          loop_len = ModNz(raw_10_ct, kBitsPerWord);
        } else {
          fset_bits = patch_10_fsetw[fset_widx];
        }
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_xys) {
            cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare10_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_10_fvalsw[fvals_widx];
              }
              if (!allele_idx_m1) {
                fvals_bits = (detect_hom_mask_hi & (~(fvals_bits | ((fvals_bits | detect_hom_mask_hi) - detect_hom_mask_lo)))) >> (code10_width - 1);
              }
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare10_lowbits = 0;
            }
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
            if (sample_include_hw[sample_hwidx] & (1U << sample_uidx_lowbits)) {
              ++subsetted_rare10_ct;
              het_incr += (fvals_bits >> rare10_lowbits) & 1;
            }
            rare10_lowbits += code10_width;
          }
          cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
          fset_bits = fset_bits >> 1;
        }
      }
      if (allele_ct == 3) {
        if (allele_idx_m1) {
          *hom_ctp = het_incr;
          *het_ctp += subsetted_rare10_ct - het_incr;
          return kPglRetSuccess;
        }
        het_incr = subsetted_rare10_ct - het_incr;
      }
      *hom_ctp = subsetted_10_ct - subsetted_rare10_ct;
      *het_ctp += het_incr;
      return kPglRetSuccess;
    }
    // allele_idx > 1, allele_ct > 3
    const uint32_t allele_code_width = 1U << allele_code_logwidth;
    const uintptr_t detect_all_mask_lo = detect_hom_mask_lo | (detect_hom_mask_lo << allele_code_width);
    const uintptr_t detect_all_mask_hi = detect_all_mask_lo << (allele_code_width - 1);
    const uintptr_t xor_word = allele_idx_m1 * detect_all_mask_lo;
    uint32_t matching_allele_ct = 0;  // 2x hom + 1x het
    uint32_t matching_het_or_hom_ct = 0;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          break;
        }
        fset_bits = SubwordLoad(&(patch_10_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_10_ct, kBitsPerWord);
      } else {
        fset_bits = patch_10_fsetw[fset_widx];
      }
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        while (!cur_raw_genovec_xys) {
          cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
        }
        if (fset_bits & 1) {
          if (rare10_lowbits == kBitsPerWord) {
            if (fvals_widx == fvals_word_ct_m1) {
              fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
            } else {
              fvals_bits = patch_10_fvalsw[fvals_widx];
            }
            fvals_bits ^= xor_word;
            fvals_bits = (detect_all_mask_hi & (~(fvals_bits | ((fvals_bits | detect_all_mask_hi) - detect_all_mask_lo)))) >> (allele_code_width - 1);
            // unnecessary to apply bzhi or detect_hom_mask_lo here
            fvals_bits = fvals_bits + (fvals_bits >> allele_code_width);
            ++fvals_widx;
            rare10_lowbits = 0;
          }
          const uintptr_t cur_hit_ct = (fvals_bits >> rare10_lowbits) & 3;
          rare10_lowbits += code10_width;
          if (cur_hit_ct) {
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
            if (sample_include_hw[sample_hwidx] & (1U << sample_uidx_lowbits)) {
              ++matching_het_or_hom_ct;
              matching_allele_ct += cur_hit_ct;
            }
          }
        }
        cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
        fset_bits = fset_bits >> 1;
      }
    }
    const uint32_t matching_hom_ct = matching_allele_ct - matching_het_or_hom_ct;
    *hom_ctp = matching_hom_ct;
    *het_ctp += matching_het_or_hom_ct - matching_hom_ct;
    return kPglRetSuccess;
  }
  // mode 1: difflist.
  if (!sample_include) {
    const unsigned char* group_info_iter;
    uint32_t rare10_ct;
    PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, &rare10_ct);
    // rare10_ct == 0 should be impossible
    if (unlikely(reterr)) {
      return reterr;
    }
    reterr = SkipDeltalistIds(fread_end, group_info_iter, rare10_ct, raw_sample_ct, 0, fread_pp);
    if (unlikely(reterr)) {
      return reterr;
    }
    const unsigned char* patch_10_fvals = *fread_pp;
    const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) << code10_logwidth, CHAR_BIT);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    CountAux1bDense(patch_10_fvals, allele_ct, allele_idx_m1, raw_10_ct, rare10_ct, het_ctp, hom_ctp);
    return kPglRetSuccess;
  }
  // Save deltalist elements, iterate.
  uint32_t rare10_ct;
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare10_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) << code10_logwidth, CHAR_BIT);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  if ((!allele_idx_m1) || (allele_ct == 3)) {
    const uintptr_t detect_hom_mask_hi = detect_hom_mask_lo << (code10_width - 1);
    uint32_t subsetted_rare10_ct = 0;
    uint32_t het_incr = 0;
    uint32_t loop_len = kBitsPerWord >> code10_logwidth;
    for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
      uintptr_t fvals_bits;
      if (fvals_widx >= fvals_word_ct_m1) {
        if (fvals_widx > fvals_word_ct_m1) {
          break;
        }
        fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
        loop_len = 1 + ((rare10_ct - 1) & ((kBitsPerWord >> code10_logwidth) - 1));
      } else {
        fvals_bits = patch_10_fvalsw[fvals_widx];
      }
      if (!allele_idx_m1) {
        fvals_bits = (detect_hom_mask_hi & (~(fvals_bits | ((fvals_bits | detect_hom_mask_hi) - detect_hom_mask_lo)))) >> (code10_width - 1);
      }
      const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - code10_logwidth)]);
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uint32_t sample_uidx = cur_deltalist_base[uii];
        if (IsSet(sample_include, sample_uidx)) {
          ++subsetted_rare10_ct;
          het_incr += (fvals_bits >> (uii << code10_logwidth)) & 1;
        }
      }
    }
    if (allele_ct == 3) {
      if (allele_idx_m1) {
        *hom_ctp = het_incr;
        *het_ctp += subsetted_rare10_ct - het_incr;
        return kPglRetSuccess;
      }
      het_incr = subsetted_rare10_ct - het_incr;
    }
    *hom_ctp = subsetted_10_ct - subsetted_rare10_ct;
    *het_ctp += het_incr;
    return kPglRetSuccess;
  }
  // allele_idx > 1, allele_ct > 3
  const uint32_t allele_code_width = 1U << allele_code_logwidth;
  const uintptr_t detect_all_mask_lo = detect_hom_mask_lo | (detect_hom_mask_lo << allele_code_width);
  const uintptr_t detect_all_mask_hi = detect_all_mask_lo << (allele_code_width - 1);
  detect_hom_mask_lo = detect_hom_mask_lo * 3;
  const uintptr_t xor_word = allele_idx_m1 * detect_all_mask_lo;
  uint32_t matching_het_or_hom_ct = 0;
  uint32_t matching_hom_ct = 0;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        break;
      }
      fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
    } else {
      fvals_bits = patch_10_fvalsw[fvals_widx];
    }
    fvals_bits = fvals_bits ^ xor_word;
    fvals_bits = detect_all_mask_hi & (~(fvals_bits | ((fvals_bits | detect_all_mask_hi) - detect_all_mask_lo)));
    if (fvals_widx == fvals_word_ct_m1) {
      fvals_bits = bzhi_max(fvals_bits, ModNz(rare10_ct << code10_logwidth, kBitsPerWord));
    }
    if (!fvals_bits) {
      continue;
    }
    fvals_bits = fvals_bits >> (allele_code_width - 1);
    fvals_bits = (fvals_bits + (fvals_bits >> allele_code_width)) & detect_hom_mask_lo;
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - code10_logwidth)]);
    do {
      const uint32_t bit_idx = ctzw(fvals_bits);
      const uint32_t sample_uidx = cur_deltalist_base[bit_idx >> code10_logwidth];
      if (IsSet(sample_include, sample_uidx)) {
        ++matching_het_or_hom_ct;
        matching_hom_ct += bit_idx & 1;
      }
      fvals_bits &= fvals_bits - 1;
    } while (fvals_bits);
  }
  *hom_ctp = matching_hom_ct;
  *het_ctp += matching_het_or_hom_ct - matching_hom_ct;
  return kPglRetSuccess;
}

PglErr PgrGetInv1Counts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // May use workspace_vec and workspace_difflist_sample_ids.
  if (!sample_ct) {
    STD_ARRAY_REF_FILL0(4, genocounts);
    return kPglRetSuccess;
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  PglErr reterr;
  if ((!allele_idx) || (!allele_idx_offsets)) {
  PgrGetInv1Counts_biallelic:
    reterr = GetBasicGenotypeCounts(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, genocounts);
    if (allele_idx) {
      const uint32_t homref_ct = genocounts[0];
      genocounts[0] = genocounts[2];
      genocounts[2] = homref_ct;
    }
    return reterr;
  }
  const uint32_t allele_ct = allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx];
  if (allele_ct == 2) {
    goto PgrGetInv1Counts_biallelic;
  }
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* tmp_genovec = pgrp->workspace_vec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  reterr = ReadRawGenovec(subsetting_required, vidx, pgrp, &fread_ptr, &fread_end, tmp_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  ZeroTrailingQuaters(raw_sample_ct, tmp_genovec);
  const uint32_t aux1_first_byte = *fread_ptr++;
  const uint32_t aux1a_mode = aux1_first_byte & 15;
  const uint32_t aux1b_mode = aux1_first_byte >> 4;
  // raw_01_ct not needed when aux1a uses difflist form and subsetting is
  // occurring; same applies to raw_10_ct.
  uint32_t raw_01_ct = 0;
  uint32_t raw_10_ct = 0;
  if ((!subsetting_required) || (!aux1a_mode) || (!aux1b_mode)) {
    GenovecCountFreqsUnsafe(tmp_genovec, raw_sample_ct, genocounts);
    raw_01_ct = genocounts[1];
    raw_10_ct = genocounts[2];
  }
  uint32_t subsetted_01_ct = 0;
  uint32_t subsetted_10_ct = 0;
  if (subsetting_required) {
    // need accurate subsetted missing count for allele_idx > 1 case
    GenovecCountSubsetFreqs(tmp_genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
    subsetted_01_ct = genocounts[1];
    subsetted_10_ct = genocounts[2];
  } else {
    sample_include = nullptr;
  }
  uint32_t het_ct;
  reterr = CountAux1a(fread_end, sample_include, tmp_genovec, aux1a_mode, raw_sample_ct, allele_ct, allele_idx, raw_01_ct, subsetted_01_ct, &fread_ptr, &het_ct, pgrp->workspace_difflist_sample_ids);
  if (unlikely(reterr)) {
    return reterr;
  }
  uint32_t hom_ct;
  reterr = CountAux1b(fread_end, sample_include, tmp_genovec, aux1b_mode, raw_sample_ct, allele_ct, allele_idx, raw_10_ct, subsetted_10_ct, &fread_ptr, &het_ct, &hom_ct, pgrp->workspace_difflist_sample_ids);
  genocounts[0] = hom_ct;
  genocounts[1] = het_ct;
  genocounts[2] = sample_ct - genocounts[3] - hom_ct - het_ct;
  return reterr;
}

// sample_include assumed to be nullptr if no subsetting required
PglErr GenovecAux1aUpdate(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict raw_genovec, uint32_t aux1a_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t allele_idx, uintptr_t lshifted_bit, uint32_t raw_01_ct, const unsigned char** fread_pp, uintptr_t* __restrict target_genovec, uint32_t* __restrict deltalist_workspace) {
  if (aux1a_mode == 15) {
    return kPglRetSuccess;
  }
  const uint32_t ignore_01_fvals = (allele_idx == 1) || (allele_ct == 3);
  uintptr_t detect_mask_hi;
  uintptr_t detect_mask_lo;
  uint32_t allele_code_logwidth;
  const uint32_t allele_code_width = GetAux1aConsts(allele_ct, &detect_mask_hi, &detect_mask_lo, &allele_code_logwidth);
  const uintptr_t xor_word = (allele_idx - 2) * detect_mask_lo;
  if (!aux1a_mode) {
#ifdef __arm__
#  error "Unaligned accesses in GenovecAux1aUpdate()."
#endif
    const uintptr_t* patch_01_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fset_byte_ct = DivUp(raw_01_ct, 8);
    uint32_t rare01_ct = 0;
    if (allele_ct > 3) {
      rare01_ct = PopcountBytes(*fread_pp, fset_byte_ct);
    }
    *fread_pp += fset_byte_ct;
    const uintptr_t* patch_01_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_hets = Word01(raw_genovec[0]);
    uint32_t loop_len = kBitsPerWord;
    const uintptr_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    const uint32_t lshift = lshifted_bit - 1;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t rare01_lowbits = kBitsPerWord;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          return kPglRetSuccess;
        }
        fset_bits = SubwordLoad(&(patch_01_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_01_ct, kBitsPerWord);
      } else {
        fset_bits = patch_01_fsetw[fset_widx];
      }
      if (!sample_include) {
        if (ignore_01_fvals) {
          for (uint32_t uii = 0; uii != loop_len; ++uii) {
            while (!cur_raw_genovec_hets) {
              cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
            }
            if (fset_bits & 1) {
              // ref/altx present for x>1.  Change genovec entry from 01 to 11
              // (or 11 -> 01 in allele_idx == 2, allele_ct == 3 case; same xor
              // operation works for that)
              const uintptr_t lowbit = cur_raw_genovec_hets & (-cur_raw_genovec_hets);
              target_genovec[sample_hwidx] ^= lowbit << lshift;
            }
            cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
            fset_bits = fset_bits >> 1;
          }
        } else {
          for (uint32_t uii = 0; uii != loop_len; ++uii) {
            while (!cur_raw_genovec_hets) {
              cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
            }
            if (fset_bits & 1) {
              if (rare01_lowbits == kBitsPerWord) {
                if (fvals_widx == fvals_word_ct_m1) {
                  fvals_bits = SubwordLoad(&(patch_01_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
                } else {
                  fvals_bits = patch_01_fvalsw[fvals_widx];
                }
                fvals_bits = fvals_bits ^ xor_word;
                fvals_bits = (detect_mask_hi & (~(fvals_bits | ((fvals_bits | detect_mask_hi) - detect_mask_lo)))) >> (allele_code_width - 1);
                // unnecessary to apply bzhi here
                ++fvals_widx;
                rare01_lowbits = 0;
              }
              if (fvals_bits & (k1LU << rare01_lowbits)) {
                const uintptr_t lowbit = cur_raw_genovec_hets & (-cur_raw_genovec_hets);
                target_genovec[sample_hwidx] ^= lowbit << lshift;
              }
              rare01_lowbits += allele_code_width;
            }
            cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
            fset_bits = fset_bits >> 1;
          }
        }
      } else {
        // format 0, sample_include non-null
        if (ignore_01_fvals) {
          for (uint32_t uii = 0; uii != loop_len; ++uii) {
            while (!cur_raw_genovec_hets) {
              cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
            }
            if (fset_bits & 1) {
              // Considered replacing cur_raw_genovec_hets with the result of
              // two PackWordToHalfword() operations, since that keeps all
              // the sample word-indexes aligned.  Couldn't justify it given
              // the expected sparsity of this case, though.
              const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_hets) / 2;
              if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
                const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
                target_genovec[sample_idx / kBitsPerWordD2] ^= lshifted_bit << (2 * (sample_idx % kBitsPerWordD2));
              }
            }
            cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
            fset_bits = fset_bits >> 1;
          }
        } else {
          for (uint32_t uii = 0; uii != loop_len; ++uii) {
            while (!cur_raw_genovec_hets) {
              cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
            }
            if (fset_bits & 1) {
              if (rare01_lowbits == kBitsPerWord) {
                if (fvals_widx == fvals_word_ct_m1) {
                  fvals_bits = SubwordLoad(&(patch_01_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
                } else {
                  fvals_bits = patch_01_fvalsw[fvals_widx];
                }
                fvals_bits = fvals_bits ^ xor_word;
                fvals_bits = (detect_mask_hi & (~(fvals_bits | ((fvals_bits | detect_mask_hi) - detect_mask_lo)))) >> (allele_code_width - 1);
                // unnecessary to apply bzhi here
                ++fvals_widx;
                rare01_lowbits = 0;
              }
              if (fvals_bits & (k1LU << rare01_lowbits)) {
                const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_hets) / 2;
                if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
                  const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
                  target_genovec[sample_idx / kBitsPerWordD2] ^= lshifted_bit << (2 * (sample_idx % kBitsPerWordD2));
                }
              }
              rare01_lowbits += allele_code_width;
            }
            cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
            fset_bits = fset_bits >> 1;
          }
        }
      }
    }
  }
  // aux1a_mode == 1
  uint32_t rare01_ct;
  // Might hardcode the ParseAndSaveDeltalist logic later, but lets get
  // this working first.
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare01_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_01_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uintptr_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  if (ignore_01_fvals) {
    if (!sample_include) {
      for (uint32_t rare01_idx = 0; rare01_idx != rare01_ct; ++rare01_idx) {
        const uint32_t sample_uidx = deltalist_workspace[rare01_idx];
        // todo: benchmark against k1LU << (lshift + ...)
        target_genovec[sample_uidx / kBitsPerWordD2] ^= lshifted_bit << (2 * (sample_uidx % kBitsPerWordD2));
      }
      return kPglRetSuccess;
    }
    for (uint32_t rare01_idx = 0; rare01_idx != rare01_ct; ++rare01_idx) {
      const uint32_t sample_uidx = deltalist_workspace[rare01_idx];
      // could wrap this boilerplate
      const uint32_t sample_widx = sample_uidx / kBitsPerWord;
      const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
      const uintptr_t sample_include_word = sample_include[sample_widx];
      if (sample_include_word & lowbit) {
        const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
        target_genovec[sample_idx / kBitsPerWordD2] ^= lshifted_bit << (2 * (sample_idx % kBitsPerWordD2));
      }
    }
    return kPglRetSuccess;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        return kPglRetSuccess;
      }
      fvals_bits = SubwordLoad(&(patch_01_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
    } else {
      fvals_bits = patch_01_fvalsw[fvals_widx];
    }
    fvals_bits = fvals_bits ^ xor_word;
    fvals_bits = detect_mask_hi & (~(fvals_bits | ((fvals_bits | detect_mask_hi) - detect_mask_lo)));
    if (fvals_widx == fvals_word_ct_m1) {
      fvals_bits = bzhi_max(fvals_bits, ModNz(rare01_ct << allele_code_logwidth, kBitsPerWord));
    }
    if (!fvals_bits) {
      continue;
    }
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - allele_code_logwidth)]);
    if (!sample_include) {
      do {
        const uint32_t rare01_idx_lowbits = ctzw(fvals_bits) >> allele_code_logwidth;
        const uint32_t sample_uidx = cur_deltalist_base[rare01_idx_lowbits];
        target_genovec[sample_uidx / kBitsPerWordD2] ^= lshifted_bit << (2 * (sample_uidx % kBitsPerWordD2));
        fvals_bits &= fvals_bits - 1;
      } while (fvals_bits);
    } else {
      do {
        const uint32_t rare01_idx_lowbits = ctzw(fvals_bits) >> allele_code_logwidth;
        const uint32_t sample_uidx = cur_deltalist_base[rare01_idx_lowbits];
        const uint32_t sample_widx = sample_uidx / kBitsPerWord;
        const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
        const uintptr_t sample_include_word = sample_include[sample_widx];
        if (sample_include_word & lowbit) {
          const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
          target_genovec[sample_idx / kBitsPerWordD2] ^= lshifted_bit << (2 * (sample_idx % kBitsPerWordD2));
        }
        fvals_bits &= fvals_bits - 1;
      } while (fvals_bits);
    }
  }
}

// sample_include assumed to be nullptr if no subsetting required
PglErr GenovecAux1bStandardUpdate(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict raw_genovec, uint32_t aux1b_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t allele_idx, uint32_t raw_10_ct, const unsigned char** fread_pp, uintptr_t* __restrict target_genovec, uint32_t* __restrict deltalist_workspace) {
  if (aux1b_mode == 15) {
    return kPglRetSuccess;
  }
  const uint32_t allele_idx_m1 = allele_idx - 1;
  uintptr_t detect_hom_mask_lo;
  const uint32_t allele_code_logwidth = GetAux1bConsts(allele_ct, &detect_hom_mask_lo);
  const uint32_t code10_logwidth = allele_code_logwidth + (allele_code_logwidth != 0);
  const uint32_t code10_width = 1U << code10_logwidth;
  uint32_t rare10_lowbits = kBitsPerWord;
  if (!aux1b_mode) {
#ifdef __arm__
#  error "Unaligned accesses in GenovecAux1bStandardUpdate()."
#endif
    const uintptr_t* patch_10_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fset_byte_ct = DivUp(raw_10_ct, 8);
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t rare10_ct = PopcountBytes(*fread_pp, fset_byte_ct);
    *fread_pp += fset_byte_ct;
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_xys = Word10(raw_genovec[0]);
    const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) * code10_width, CHAR_BIT);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t loop_len = kBitsPerWord;
    if ((!allele_idx_m1) || (allele_ct == 3)) {
      const uintptr_t detect_hom_mask_hi = detect_hom_mask_lo << (code10_width - 1);
      const uint32_t lowcode_add = 2 - allele_idx_m1;
      for (uint32_t fset_widx = 0; ; ++fset_widx) {
        uintptr_t fset_bits;
        if (fset_widx >= fset_word_ct_m1) {
          if (fset_widx > fset_word_ct_m1) {
            break;
          }
          fset_bits = SubwordLoad(&(patch_10_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
          loop_len = ModNz(raw_10_ct, kBitsPerWord);
        } else {
          fset_bits = patch_10_fsetw[fset_widx];
        }
        if (!sample_include) {
          for (uint32_t uii = 0; uii != loop_len; ++uii) {
            while (!cur_raw_genovec_xys) {
              cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
            }
            if (fset_bits & 1) {
              if (rare10_lowbits == kBitsPerWord) {
                if (fvals_widx == fvals_word_ct_m1) {
                  fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
                } else {
                  fvals_bits = patch_10_fvalsw[fvals_widx];
                }
                // modify to het 1/x = 1, otherwise 0, except in allele_idx ==
                // 2 special case.
                if (!allele_idx_m1) {
                  fvals_bits = (detect_hom_mask_hi & (~(fvals_bits | ((fvals_bits | detect_hom_mask_hi) - detect_hom_mask_lo)))) >> (code10_width - 1);
                }
                // unnecessary to apply bzhi here
                ++fvals_widx;
                rare10_lowbits = 0;
              }
              const uint32_t cur_lowcode0 = (fvals_bits >> rare10_lowbits) & 1;
              rare10_lowbits += code10_width;
              const uintptr_t lowbit = cur_raw_genovec_xys & (-cur_raw_genovec_xys);
              target_genovec[sample_hwidx] ^= lowbit * (lowcode_add + cur_lowcode0);
            }
            cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
            fset_bits = fset_bits >> 1;
          }
        } else {
          // sample_include non-null
          for (uint32_t uii = 0; uii != loop_len; ++uii) {
            while (!cur_raw_genovec_xys) {
              cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
            }
            if (fset_bits & 1) {
              if (rare10_lowbits == kBitsPerWord) {
                if (fvals_widx == fvals_word_ct_m1) {
                  fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
                } else {
                  fvals_bits = patch_10_fvalsw[fvals_widx];
                }
                // modify to het 1/x = 1, otherwise 0, except in allele_idx ==
                // 2 special case
                if (!allele_idx_m1) {
                  fvals_bits = (detect_hom_mask_hi & (~(fvals_bits | ((fvals_bits | detect_hom_mask_hi) - detect_hom_mask_lo)))) >> (code10_width - 1);
                }
                // unnecessary to apply bzhi here
                ++fvals_widx;
                rare10_lowbits = 0;
              }
              const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
              if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
                const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
                const uintptr_t cur_lowcode0 = (fvals_bits >> rare10_lowbits) & 1;
                const uintptr_t shifted_xor_mult = (lowcode_add + cur_lowcode0) << (2 * (sample_idx % kBitsPerWordD2));
                target_genovec[sample_idx / kBitsPerWordD2] ^= shifted_xor_mult;
              }
              rare10_lowbits += code10_width;
            }
            cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
            fset_bits = fset_bits >> 1;
          }
        }
      }
      return kPglRetSuccess;
    }
    // allele_idx > 1, allele_ct > 3
    const uint32_t allele_code_width = 1U << allele_code_logwidth;
    const uintptr_t detect_all_mask_lo = detect_hom_mask_lo | (detect_hom_mask_lo << allele_code_width);
    const uintptr_t detect_all_mask_hi = detect_all_mask_lo << (allele_code_width - 1);
    const uintptr_t xor_word = allele_idx_m1 * detect_all_mask_lo;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          break;
        }
        fset_bits = SubwordLoad(&(patch_10_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_10_ct, kBitsPerWord);
      } else {
        fset_bits = patch_10_fsetw[fset_widx];
      }
      if (!sample_include) {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_xys) {
            cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare10_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_10_fvalsw[fvals_widx];
              }
              // modify to hom = 2, het = 1, neither = 0
              fvals_bits = fvals_bits ^ xor_word;
              fvals_bits = (detect_all_mask_hi & (~(fvals_bits | ((fvals_bits | detect_all_mask_hi) - detect_all_mask_lo)))) >> (allele_code_width - 1);
              // unnecessary to apply bzhi or detect_hom_mask_lo here
              fvals_bits = fvals_bits + (fvals_bits >> allele_code_width);
              ++fvals_widx;
              rare10_lowbits = 0;
            }
            const uintptr_t cur_hit_ct = (fvals_bits >> rare10_lowbits) & 3;
            rare10_lowbits += code10_width;
            if (cur_hit_ct) {
              const uintptr_t lowbit = cur_raw_genovec_xys & (-cur_raw_genovec_xys);
              target_genovec[sample_hwidx] ^= lowbit * cur_hit_ct;
            }
          }
          cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
          fset_bits = fset_bits >> 1;
        }
      } else {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_xys) {
            cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare10_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_10_fvalsw[fvals_widx];
              }
              // modify to hom = 2, het = 1, neither = 0
              fvals_bits = fvals_bits ^ xor_word;
              fvals_bits = (detect_all_mask_hi & (~(fvals_bits | ((fvals_bits | detect_all_mask_hi) - detect_all_mask_lo)))) >> (allele_code_width - 1);
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = bzhi_max(fvals_bits, ModNz(rare10_ct * code10_width, kBitsPerWord));
              }
              fvals_bits = fvals_bits + (fvals_bits >> allele_code_width);
              ++fvals_widx;
              rare10_lowbits = 0;
            }
            const uintptr_t cur_hit_ct = (fvals_bits >> rare10_lowbits) & 3;
            rare10_lowbits += code10_width;
            if (cur_hit_ct) {
              const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
              if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
                const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
                target_genovec[sample_idx / kBitsPerWordD2] ^= cur_hit_ct << (2 * (sample_idx % kBitsPerWordD2));
              }
            }
          }
          cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
          fset_bits = fset_bits >> 1;
        }
      }
    }
    return kPglRetSuccess;
  }
  // aux1b_mode == 1
  uint32_t rare10_ct;
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare10_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) << code10_logwidth, CHAR_BIT);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  if ((!allele_idx_m1) || (allele_ct == 3)) {
    const uintptr_t detect_hom_mask_hi = detect_hom_mask_lo << (code10_width - 1);
    const uintptr_t lowcode_add = 2 - allele_idx_m1;
    uint32_t loop_len = kBitsPerWord >> code10_logwidth;
    for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
      uintptr_t fvals_bits;
      if (fvals_widx >= fvals_word_ct_m1) {
        if (fvals_widx > fvals_word_ct_m1) {
          break;
        }
        fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
        loop_len = 1 + ((rare10_ct - 1) & ((kBitsPerWord >> code10_logwidth) - 1));
      } else {
        fvals_bits = patch_10_fvalsw[fvals_widx];
      }
      if (!allele_idx_m1) {
        fvals_bits = (detect_hom_mask_hi & (~(fvals_bits | ((fvals_bits | detect_hom_mask_hi) - detect_hom_mask_lo)))) >> (code10_width - 1);
      }
      const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - code10_logwidth)]);
      if (!sample_include) {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          const uint32_t sample_uidx = cur_deltalist_base[uii];
          const uintptr_t cur_lowcode0 = fvals_bits & 1;
          const uintptr_t shifted_xor_mult = (lowcode_add + cur_lowcode0) << (2 * (sample_uidx % kBitsPerWordD2));
          target_genovec[sample_uidx / kBitsPerWordD2] ^= shifted_xor_mult;
          fvals_bits = fvals_bits >> code10_width;
        }
      } else {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          const uint32_t sample_uidx = cur_deltalist_base[uii];
          const uint32_t sample_widx = sample_uidx / kBitsPerWord;
          const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
          const uintptr_t sample_include_word = sample_include[sample_widx];
          if (sample_include_word & lowbit) {
            const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
            const uintptr_t cur_lowcode0 = fvals_bits & 1;
            const uintptr_t shifted_xor_mult = (lowcode_add + cur_lowcode0) << (2 * (sample_idx % kBitsPerWordD2));
            target_genovec[sample_idx / kBitsPerWordD2] ^= shifted_xor_mult;
          }
          fvals_bits = fvals_bits >> code10_width;
        }
      }
    }
    return kPglRetSuccess;
  }
  // allele_idx > 1, allele_ct > 3
  const uint32_t allele_code_width = 1U << allele_code_logwidth;
  const uintptr_t detect_all_mask_lo = detect_hom_mask_lo | (detect_hom_mask_lo << allele_code_width);
  const uintptr_t detect_all_mask_hi = detect_all_mask_lo << (allele_code_width - 1);
  detect_hom_mask_lo = detect_hom_mask_lo * 3;
  const uintptr_t xor_word = allele_idx_m1 * detect_all_mask_lo;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        break;
      }
      fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
    } else {
      fvals_bits = patch_10_fvalsw[fvals_widx];
    }
    fvals_bits = fvals_bits ^ xor_word;
    fvals_bits = detect_all_mask_hi & (~(fvals_bits | ((fvals_bits | detect_all_mask_hi) - detect_all_mask_lo)));
    if (fvals_widx == fvals_word_ct_m1) {
      fvals_bits = bzhi_max(fvals_bits, ModNz(rare10_ct << code10_logwidth, kBitsPerWord));
    }
    if (!fvals_bits) {
      continue;
    }
    fvals_bits = fvals_bits >> (allele_code_width - 1);
    fvals_bits = (fvals_bits + (fvals_bits >> allele_code_width)) & detect_hom_mask_lo;
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - code10_logwidth)]);
    if (!sample_include) {
      do {
        const uint32_t bit_idx = ctzw(fvals_bits);
        const uint32_t sample_uidx = cur_deltalist_base[bit_idx >> code10_logwidth];
        target_genovec[sample_uidx / kBitsPerWordD2] ^= k1LU << ((bit_idx % 2) + 2 * (sample_uidx % kBitsPerWordD2));
        fvals_bits &= fvals_bits - 1;
      } while (fvals_bits);
    } else {
      do {
        const uint32_t bit_idx = ctzw(fvals_bits);
        const uint32_t sample_uidx = cur_deltalist_base[bit_idx >> code10_logwidth];
        const uint32_t sample_widx = sample_uidx / kBitsPerWord;
        const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
        const uintptr_t sample_include_word = sample_include[sample_widx];
        if (sample_include_word & lowbit) {
          const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
          target_genovec[sample_idx / kBitsPerWordD2] ^= k1LU << ((bit_idx % 2) + 2 * (sample_idx % kBitsPerWordD2));
        }
        fvals_bits &= fvals_bits - 1;
      } while (fvals_bits);
    }
  }
  return kPglRetSuccess;
}

// if aux1b_het_present is true, aux1b_hets becomes a 1-bit-per-sample bitarray
// with the positions of altx/alty hets in aux1b.
PglErr GetAux1bHets(const unsigned char* fread_end, const uintptr_t* __restrict raw_genovec, uint32_t aux1b_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t raw_10_ct, const unsigned char** fread_pp, uintptr_t* __restrict aux1b_hets, uint32_t* __restrict aux1b_het_presentp, uint32_t* __restrict deltalist_workspace) {
  if (aux1b_mode == 15) {
    *aux1b_het_presentp = 0;
    return kPglRetSuccess;
  }
  uintptr_t detect_hom_mask_lo;
  const uint32_t allele_code_logwidth = GetAux1bConsts(allele_ct, &detect_hom_mask_lo);
  const uint32_t code10_logwidth = allele_code_logwidth + (allele_code_logwidth != 0);
  const uint32_t code10_width = 1U << code10_logwidth;
  const uint32_t allele_code_width = 1U << allele_code_logwidth;
  const uintptr_t detect_all_mask_lo = detect_hom_mask_lo | (detect_hom_mask_lo << allele_code_width);
  const uintptr_t detect_all_mask_hi = detect_all_mask_lo << (allele_code_width - 1);
  Halfword* aux1b_hets_alias = R_CAST(Halfword*, aux1b_hets);
  uint32_t rare10_lowbits = kBitsPerWord;
  uint32_t aux1b_het_present = 0;
  if (!aux1b_mode) {
#ifdef __arm__
#  error "Unaligned accesses in GetAux1bHets()."
#endif
    const uintptr_t* patch_10_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fset_byte_ct = DivUp(raw_10_ct, 8);
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t rare10_ct = PopcountBytes(*fread_pp, fset_byte_ct);
    *fread_pp += fset_byte_ct;
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_xys = Word10(raw_genovec[0]);
    const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) * code10_width, CHAR_BIT);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t loop_len = kBitsPerWord;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          break;
        }
        fset_bits = SubwordLoad(&(patch_10_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_10_ct, kBitsPerWord);
      } else {
        fset_bits = patch_10_fsetw[fset_widx];
      }
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        while (!cur_raw_genovec_xys) {
          cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
        }
        if (fset_bits & 1) {
          if (rare10_lowbits == kBitsPerWord) {
            if (fvals_widx == fvals_word_ct_m1) {
              fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
            } else {
              fvals_bits = patch_10_fvalsw[fvals_widx];
            }
            // allele_ct == 3: just invert raw fvals_bits
            // allele_ct > 3: shift by allele_code_width, xor with self so that
            // 0 == hom, detect nonzero by inverting the usual check
            if (allele_ct == 3) {
              fvals_bits = ~fvals_bits;
            } else {
              fvals_bits = fvals_bits ^ (fvals_bits << allele_code_width);
              // conveniently, removing a ~ here is equivalent to inverting the
              // relevant bits of the final result
              fvals_bits = detect_hom_mask_lo & ((fvals_bits | ((fvals_bits | detect_all_mask_hi) - detect_all_mask_lo)) >> (code10_width - 1));
            }
            // bzhi only relevant for detecting if there are any hets at all
            if (!aux1b_het_present) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = bzhi_max(fvals_bits, ModNz(rare10_ct * code10_width, kBitsPerWord));
              }
              if (fvals_bits) {
                // lazy-initialize
                aux1b_het_present = 1;
                ZeroHwArr(2 * BitCtToWordCt(raw_sample_ct), aux1b_hets_alias);
              }
            }
            ++fvals_widx;
            rare10_lowbits = 0;
          }
          if (fvals_bits & (k1LU << rare10_lowbits)) {
            const uint32_t bit_idx = ctzw(cur_raw_genovec_xys) / 2;
            aux1b_hets_alias[sample_hwidx] |= 1U << bit_idx;
          }
          rare10_lowbits += code10_width;
        }
        cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
        fset_bits = fset_bits >> 1;
      }
    }
    *aux1b_het_presentp = aux1b_het_present;
    return kPglRetSuccess;
  }
  // aux1b_mode == 1
  uint32_t rare10_ct;
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare10_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) << code10_logwidth, CHAR_BIT);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        break;
      }
      fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
    } else {
      fvals_bits = patch_10_fvalsw[fvals_widx];
    }
    if (allele_ct == 3) {
      fvals_bits = ~fvals_bits;
    } else {
      fvals_bits = fvals_bits ^ (fvals_bits << allele_code_width);
      fvals_bits = detect_hom_mask_lo & ((fvals_bits | ((fvals_bits | detect_all_mask_hi) - detect_all_mask_lo)) >> (code10_width - 1));
    }
    if (fvals_widx == fvals_word_ct_m1) {
      fvals_bits = bzhi_max(fvals_bits, ModNz(rare10_ct << code10_logwidth, kBitsPerWord));
    }
    if (!fvals_bits) {
      continue;
    }
    if (!aux1b_het_present) {
      aux1b_het_present = 1;
      ZeroHwArr(2 * BitCtToWordCt(raw_sample_ct), aux1b_hets_alias);
    }
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - code10_logwidth)]);
    do {
      const uint32_t bit_idx = ctzw(fvals_bits);
      const uint32_t sample_uidx = cur_deltalist_base[bit_idx >> code10_logwidth];
      aux1b_hets_alias[sample_uidx / kBitsPerWordD2] |= 1U << (sample_uidx % kBitsPerWordD2);
      fvals_bits &= fvals_bits - 1;
    } while (fvals_bits);
  }
  *aux1b_het_presentp = aux1b_het_present;
  return kPglRetSuccess;
}

PglErr Get1Multiallelic(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict all_hets, uintptr_t* __restrict allele_countvec, uintptr_t** subsetted_10hetp) {
  // sample_ct > 0; either allele_idx > 1 or ((allele_idx == 1) &&
  // multiallelic_hc_present)
  // subsetted_10het assumed to be initialized to nullptr, if present at all
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* raw_genovec = pgrp->workspace_vec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRawGenovec(subsetting_required, vidx, pgrp, &fread_ptr, &fread_end, raw_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }

  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  const uint32_t multiallelic_hc_present = VrtypeMultiallelicHc(vrtype);
  CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, allele_countvec);
  ZeroTrailingQuaters(raw_sample_ct, raw_genovec);
  if (fread_pp) {
    *fread_endp = fread_end;
    if (all_hets) {
      PgrDetectGenovecHets(raw_genovec, raw_sample_ct, all_hets);
    }
  }
  if (allele_idx != 1) {
    GenovecNonmissingToZeroUnsafe(sample_ct, allele_countvec);
    if (!multiallelic_hc_present) {
      if (fread_pp) {
        *fread_pp = fread_ptr;
      }
      return kPglRetSuccess;
    }
  }
  const uint32_t aux1_first_byte = *fread_ptr++;
  const uint32_t aux1a_mode = aux1_first_byte & 15;
  const uint32_t aux1b_mode = aux1_first_byte >> 4;
  // only need to initialize these in dense modes
  uint32_t raw_01_ct = 0;
  uint32_t raw_10_ct = 0;
  if ((!aux1a_mode) || (!aux1b_mode)) {
    GenovecCount12Unsafe(raw_genovec, raw_sample_ct, &raw_01_ct, &raw_10_ct);
  }

  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx];
  if (!subsetting_required) {
    sample_include = nullptr;
  }
  // allele_idx == 1 case:
  //   allele_countvec currently contains ALT counts; we want to reduce them to
  //   ALT1 counts.  This can be done with the following steps:
  //   1. For every element of patch_01_fset, reduce the value from 1 to 0.  We
  //      don't actually need to look at patch_01_fvals.
  //   2. For every element of patch_10_fset, reduce the value from 2 depending
  //      on the low bit(s) of the patch_01_fvals entry (reduce to 0 unless low
  //      bit(s) are all zero).
  // allele_idx > 1 case:
  //   1. For every element of patch_01_fset, set a 1 for each matching value
  //      of patch_01_fvals.
  //   2. For every element of patch_10_fset, set a 1 for each het-matching
  //      value of patch_10_fvals, and a 2 for each hom-match.
  uint32_t* deltalist_workspace = pgrp->workspace_difflist_sample_ids;
  // Two cases:
  // - If allele_idx == 1, convert all aux1a entries from 01 to 00.
  // - Otherwise, for each matching aux1a entry, convert from 00 to 01.
  reterr = GenovecAux1aUpdate(fread_end, sample_include, sample_include_cumulative_popcounts, raw_genovec, aux1a_mode, raw_sample_ct, allele_ct, allele_idx, 1, raw_01_ct, &fread_ptr, allele_countvec, deltalist_workspace);
  if (unlikely(reterr)) {
    return reterr;
  }
  const unsigned char* aux1b_start = fread_ptr;
  reterr = GenovecAux1bStandardUpdate(fread_end, sample_include, sample_include_cumulative_popcounts, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, allele_idx, raw_10_ct, &fread_ptr, allele_countvec, deltalist_workspace);
  if ((!fread_pp) || reterr) {
    return reterr;
  }
  *fread_pp = fread_ptr;
  if (all_hets) {
    // can merge this with GenovecAux1bStandardUpdate if this is ever a
    // significant bottleneck
    uintptr_t* aux1b_hets = pgrp->workspace_aux1x_present;
    uint32_t aux1b_het_present;
    reterr = GetAux1bHets(fread_end, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &aux1b_start, aux1b_hets, &aux1b_het_present, deltalist_workspace);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (aux1b_het_present) {
      BitvecOr(aux1b_hets, BitCtToWordCt(raw_sample_ct), all_hets);
      if (!sample_include) {
        *subsetted_10hetp = aux1b_hets;
      } else {
        // Don't need raw_genovec any more.
        CopyBitarrSubset(aux1b_hets, sample_include, sample_ct, raw_genovec);
        *subsetted_10hetp = raw_genovec;
      }
    }
  }
  return kPglRetSuccess;
}

PglErr PgrGet1(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_countvec) {
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t multiallelic_hc_present = VrtypeMultiallelicHc(vrtype);
  if ((!allele_idx) || ((allele_idx == 1) && (!multiallelic_hc_present))) {
    PglErr reterr = ReadGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, allele_countvec);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (!allele_idx) {
      GenovecInvertUnsafe(sample_ct, allele_countvec);
    }
    return kPglRetSuccess;
  }
  return Get1Multiallelic(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, nullptr, nullptr, nullptr, allele_countvec, nullptr);
}

PglErr PgrGetInv1(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_invcountvec) {
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t multiallelic_hc_present = VrtypeMultiallelicHc(vrtype);
  if ((!allele_idx) || ((allele_idx == 1) && (!multiallelic_hc_present))) {
    PglErr reterr = ReadGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, allele_invcountvec);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (allele_idx) {
      GenovecInvertUnsafe(sample_ct, allele_invcountvec);
    }
    return kPglRetSuccess;
  }
  PglErr reterr = Get1Multiallelic(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, nullptr, nullptr, nullptr, allele_invcountvec, nullptr);
  GenovecInvertUnsafe(sample_ct, allele_invcountvec);
  return reterr;
}

// Assumes allele_idx0 < allele_idx1, and allele_idx0 < 2.  Rotates hardcalls
// such that, if no multiallelic hardcalls are present, 0 = 0/0, 1 = 0/1,
// 2 = 1/1, and 3 = anything else.
void Rotate2(uint32_t allele_idx0, uint32_t allele_idx1, uint32_t sample_ct, uintptr_t* genovec) {
  if (!allele_idx0) {
    if (allele_idx1 > 1) {
      GenovecNonzeroToMissingUnsafe(sample_ct, genovec);
    }
  } else {
    GenovecInvertThenNonzeroToMissingUnsafe(sample_ct, genovec);
  }
}

PglErr SkipAux1a(const unsigned char* fread_end, uint32_t aux1a_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t raw_01_ct, const unsigned char** fread_pp) {
  if (aux1a_mode == 15) {
    return kPglRetSuccess;
  }
  uint32_t rare01_ct;
  if (!aux1a_mode) {
    const uint32_t fset_byte_ct = DivUp(raw_01_ct, CHAR_BIT);
    rare01_ct = PopcountBytes(*fread_pp, fset_byte_ct);
    *fread_pp += fset_byte_ct;
  } else {
    const unsigned char* group_info_iter;
    PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, &rare01_ct);
    if (unlikely(reterr)) {
      return reterr;
    }
    reterr = SkipDeltalistIds(fread_end, group_info_iter, rare01_ct, raw_sample_ct, 0, fread_pp);
    if (unlikely(reterr)) {
      return reterr;
    }
  }
  const uint32_t fvals_byte_ct = GetAux1aAlleleEntryByteCt(allele_ct, rare01_ct);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  return kPglRetSuccess;
}

// sample_include assumed to be nullptr if no subsetting required
PglErr GenovecAux1bUpdate2(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict raw_genovec, uint32_t aux1b_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t allele_idx0, uint32_t allele_idx1, uint32_t raw_10_ct, const unsigned char** fread_pp, uintptr_t* __restrict target_genovec, uint32_t* __restrict deltalist_workspace) {
  // Possible aux1b updates:
  // - allele_idx0 == 0:
  //     allele_idx1 == 1: all altx/alty including a rarealt from 10 to 11
  //     allele_idx1 > 1: set one rarealtx/rarealtx from 11 to 10
  //
  // - allele_idx0 == 1: change all alt1/rarealtx from 00 to 01,
  //   rarealtx/rarealtx from 00 to 10, and all other aux1b entries to missing.
  //   This can use the same driver as Get1Multiallelic.
  //
  // - allele_idx0 > 1: change all rarealtx/rarealtx from missing to 00,
  //   rarealtx/rarealty to 01, and rarealty/rarealty to 10.
  if (aux1b_mode == 15) {
    return kPglRetSuccess;
  }
  if (allele_idx0 == 1) {
    return GenovecAux1bStandardUpdate(fread_end, sample_include, sample_include_cumulative_popcounts, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, allele_idx1, raw_10_ct, fread_pp, target_genovec, deltalist_workspace);
  }
  uintptr_t detect_hom_mask_lo;
  const uint32_t allele_code_logwidth = GetAux1bConsts(allele_ct, &detect_hom_mask_lo);
  const uint32_t code10_logwidth = allele_code_logwidth + (allele_code_logwidth != 0);
  const uint32_t code10_width = 1U << code10_logwidth;
  const uintptr_t detect_hom_mask_hi = detect_hom_mask_lo << (code10_width - 1);
  uintptr_t xor_word2 = allele_idx1 - 1;
  // fortunately, this sequence of operations happens to work for allele_ct ==
  // 3
  xor_word2 = xor_word2 | (xor_word2 << (code10_width / 2));
  xor_word2 = xor_word2 * detect_hom_mask_lo;
  uint32_t rare10_lowbits = kBitsPerWord;
  if (!aux1b_mode) {
#ifdef __arm__
#  error "Unaligned accesses in GenovecAux1bUpdate2()."
#endif
    const uintptr_t* patch_10_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fset_byte_ct = DivUp(raw_10_ct, 8);
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t rare10_ct = PopcountBytes(*fread_pp, fset_byte_ct);
    *fread_pp += fset_byte_ct;
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_xys = Word10(raw_genovec[0]);
    const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) * code10_width, CHAR_BIT);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t loop_len = kBitsPerWord;
    if (!allele_idx0) {
      for (uint32_t fset_widx = 0; ; ++fset_widx) {
        uintptr_t fset_bits;
        if (fset_widx >= fset_word_ct_m1) {
          if (fset_widx > fset_word_ct_m1) {
            return kPglRetSuccess;
          }
          fset_bits = SubwordLoad(&(patch_10_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
          loop_len = ModNz(raw_10_ct, kBitsPerWord);
        } else {
          fset_bits = patch_10_fsetw[fset_widx];
        }
        if (!sample_include) {
          if (allele_idx1 == 1) {
            // All aux1b 10 -> 11.  Ignore aux1b_fvals.
            for (uint32_t uii = 0; uii != loop_len; ++uii) {
              while (!cur_raw_genovec_xys) {
                cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
              }
              if (fset_bits & 1) {
                const uintptr_t lowbit = cur_raw_genovec_xys & (-cur_raw_genovec_xys);
                target_genovec[sample_hwidx] ^= lowbit;
              }
              cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
              fset_bits = fset_bits >> 1;
            }
          } else {
            // hom-altx 11 -> 10.
            for (uint32_t uii = 0; uii != loop_len; ++uii) {
              while (!cur_raw_genovec_xys) {
                cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
              }
              if (fset_bits & 1) {
                if (rare10_lowbits == kBitsPerWord) {
                  if (fvals_widx == fvals_word_ct_m1) {
                    fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
                  } else {
                    fvals_bits = patch_10_fvalsw[fvals_widx];
                  }
                  fvals_bits = fvals_bits ^ xor_word2;
                  fvals_bits = (detect_hom_mask_hi & (~(fvals_bits | ((fvals_bits | detect_hom_mask_hi) - detect_hom_mask_lo)))) >> (code10_width - 1);
                  // unnecessary to apply bzhi here
                  ++fvals_widx;
                  rare10_lowbits = 0;
                }
                if (fvals_bits & (k1LU << rare10_lowbits)) {
                  const uintptr_t lowbit = cur_raw_genovec_xys & (-cur_raw_genovec_xys);
                  target_genovec[sample_hwidx] ^= lowbit;
                }
                rare10_lowbits += code10_width;
              }
              cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
              fset_bits = fset_bits >> 1;
            }
          }
        } else {
          // sample_include non-null
          if (allele_idx1 == 1) {
            for (uint32_t uii = 0; uii != loop_len; ++uii) {
              while (!cur_raw_genovec_xys) {
                cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
              }
              if (fset_bits & 1) {
                const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
                if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
                  const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
                  target_genovec[sample_idx / kBitsPerWordD2] ^= k1LU << (2 * (sample_idx % kBitsPerWordD2));
                }
                rare10_lowbits += code10_width;
              }
              cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
              fset_bits = fset_bits >> 1;
            }
          } else {
            for (uint32_t uii = 0; uii != loop_len; ++uii) {
              while (!cur_raw_genovec_xys) {
                cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
              }
              if (fset_bits & 1) {
                if (rare10_lowbits == kBitsPerWord) {
                  if (fvals_widx == fvals_word_ct_m1) {
                    fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
                  } else {
                    fvals_bits = patch_10_fvalsw[fvals_widx];
                  }
                  fvals_bits = fvals_bits ^ xor_word2;
                  fvals_bits = (detect_hom_mask_hi & (~(fvals_bits | ((fvals_bits | detect_hom_mask_hi) - detect_hom_mask_lo)))) >> (code10_width - 1);
                  // unnecessary to apply bzhi here
                  ++fvals_widx;
                  rare10_lowbits = 0;
                }
                if (fvals_bits & (k1LU << rare10_lowbits)) {
                  const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
                  if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
                    const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
                    target_genovec[sample_idx / kBitsPerWordD2] ^= k1LU << (2 * (sample_idx % kBitsPerWordD2));
                  }
                }
                rare10_lowbits += code10_width;
              }
              cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
              fset_bits = fset_bits >> 1;
            }
          }
        }
      }
    }
    // 2 <= allele_idx0 < allele_idx1 (so allele_ct > 3 guaranteed)
    uintptr_t xor_word1 = allele_idx1 - 1;
    uintptr_t xor_word0 = allele_idx0 - 1;
    xor_word1 = xor_word0 | (xor_word1 << (code10_width / 2));
    xor_word0 = xor_word0 | (xor_word0 << (code10_width / 2));
    xor_word1 *= detect_hom_mask_lo;
    xor_word0 *= detect_hom_mask_lo;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          return kPglRetSuccess;
        }
        fset_bits = SubwordLoad(&(patch_10_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_10_ct, kBitsPerWord);
      } else {
        fset_bits = patch_10_fsetw[fset_widx];
      }
      if (!sample_include) {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_xys) {
            cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare10_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_10_fvalsw[fvals_widx];
              }
              uintptr_t match0 = fvals_bits ^ xor_word0;
              uintptr_t match1 = fvals_bits ^ xor_word1;
              uintptr_t match2 = fvals_bits ^ xor_word2;
              match0 = detect_hom_mask_hi & (~(match0 | ((match0 | detect_hom_mask_hi) - detect_hom_mask_lo)));
              match1 = detect_hom_mask_hi & (~(match1 | ((match1 | detect_hom_mask_hi) - detect_hom_mask_lo)));
              match2 = detect_hom_mask_hi & (~(match2 | ((match2 | detect_hom_mask_hi) - detect_hom_mask_lo)));
              // Now want match0 -> 11, match1 -> 10, and match2 -> 01.
              fvals_bits = ((match0 | match1) >> (code10_width - 2)) | ((match0 | match2) >> (code10_width - 1));
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare10_lowbits = 0;
            }
            const uintptr_t xor_val = (fvals_bits >> rare10_lowbits) & 3;
            if (xor_val) {
              const uintptr_t lowbit = cur_raw_genovec_xys & (-cur_raw_genovec_xys);
              target_genovec[sample_hwidx] ^= lowbit * xor_val;
            }
            rare10_lowbits += code10_width;
          }
          cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
          fset_bits = fset_bits >> 1;
        }
      } else {
        // sample_include non-null
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_xys) {
            cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare10_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_10_fvalsw[fvals_widx];
              }
              uintptr_t match0 = fvals_bits ^ xor_word0;
              uintptr_t match1 = fvals_bits ^ xor_word1;
              uintptr_t match2 = fvals_bits ^ xor_word2;
              match0 = detect_hom_mask_hi & (~(match0 | ((match0 | detect_hom_mask_hi) - detect_hom_mask_lo)));
              match1 = detect_hom_mask_hi & (~(match1 | ((match1 | detect_hom_mask_hi) - detect_hom_mask_lo)));
              match2 = detect_hom_mask_hi & (~(match2 | ((match2 | detect_hom_mask_hi) - detect_hom_mask_lo)));
              fvals_bits = ((match0 | match1) >> (code10_width - 2)) | ((match0 | match2) >> (code10_width - 1));
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare10_lowbits = 0;
            }
            const uintptr_t xor_val = (fvals_bits >> rare10_lowbits) & 3;
            if (xor_val) {
              const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
              if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
                const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
                target_genovec[sample_idx / kBitsPerWordD2] ^= xor_val << (2 * (sample_idx % kBitsPerWordD2));
              }
            }
            rare10_lowbits += code10_width;
          }
          cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
          fset_bits = fset_bits >> 1;
        }
      }
    }
  }
  // aux1b_mode == 1
  uint32_t rare10_ct;
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare10_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) << code10_logwidth, CHAR_BIT);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  if (allele_idx1 == 1) {
    if (!sample_include) {
      for (uint32_t rare10_idx = 0; rare10_idx != rare10_ct; ++rare10_idx) {
        const uint32_t sample_uidx = deltalist_workspace[rare10_idx];
        target_genovec[sample_uidx / kBitsPerWordD2] ^= k1LU << (2 * (sample_uidx % kBitsPerWordD2));
      }
      return kPglRetSuccess;
    }
    for (uint32_t rare10_idx = 0; rare10_idx != rare10_ct; ++rare10_idx) {
      const uint32_t sample_uidx = deltalist_workspace[rare10_idx];
      const uint32_t sample_widx = sample_uidx / kBitsPerWord;
      const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
      const uintptr_t sample_include_word = sample_include[sample_widx];
      if (sample_include_word & lowbit) {
        const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
        target_genovec[sample_idx / kBitsPerWordD2] ^= k1LU << (2 * (sample_idx % kBitsPerWordD2));
      }
    }
    return kPglRetSuccess;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  if (!allele_idx0) {
    for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
      uintptr_t fvals_bits;
      if (fvals_widx >= fvals_word_ct_m1) {
        if (fvals_widx > fvals_word_ct_m1) {
          return kPglRetSuccess;
        }
        fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
      } else {
        fvals_bits = patch_10_fvalsw[fvals_widx];
      }
      fvals_bits = fvals_bits ^ xor_word2;
      fvals_bits = detect_hom_mask_hi & (~(fvals_bits | ((fvals_bits | detect_hom_mask_hi) - detect_hom_mask_lo)));
      if (fvals_widx == fvals_word_ct_m1) {
        fvals_bits = bzhi_max(fvals_bits, ModNz(rare10_ct << code10_logwidth, kBitsPerWord));
      }
      if (!fvals_bits) {
        continue;
      }
      const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - code10_logwidth)]);
      if (!sample_include) {
        do {
          const uint32_t bit_idx = ctzw(fvals_bits);
          const uint32_t sample_uidx = cur_deltalist_base[bit_idx >> code10_logwidth];
          target_genovec[sample_uidx / kBitsPerWordD2] ^= k1LU << (2 * (sample_uidx % kBitsPerWordD2));
          fvals_bits &= fvals_bits - 1;
        } while (fvals_bits);
      } else {
        do {
          const uint32_t bit_idx = ctzw(fvals_bits);
          const uint32_t sample_uidx = cur_deltalist_base[bit_idx >> code10_logwidth];
          const uint32_t sample_widx = sample_uidx / kBitsPerWord;
          const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
          const uintptr_t sample_include_word = sample_include[sample_widx];
          if (sample_include_word & lowbit) {
            const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
            target_genovec[sample_idx / kBitsPerWordD2] ^= k1LU << (2 * (sample_idx % kBitsPerWordD2));
          }
          fvals_bits &= fvals_bits - 1;
        } while (fvals_bits);
      }
    }
  }
  // 2 <= allele_idx0 < allele_idx1
  uintptr_t xor_word1 = allele_idx1 - 1;
  uintptr_t xor_word0 = allele_idx0 - 1;
  xor_word1 = xor_word0 | (xor_word1 << (code10_width / 2));
  xor_word0 = xor_word0 | (xor_word0 << (code10_width / 2));
  xor_word1 *= detect_hom_mask_lo;
  xor_word0 *= detect_hom_mask_lo;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        return kPglRetSuccess;
      }
      fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
    } else {
      fvals_bits = patch_10_fvalsw[fvals_widx];
    }
    uintptr_t match0 = fvals_bits ^ xor_word0;
    uintptr_t match1 = fvals_bits ^ xor_word1;
    uintptr_t match2 = fvals_bits ^ xor_word2;
    match0 = detect_hom_mask_hi & (~(match0 | ((match0 | detect_hom_mask_hi) - detect_hom_mask_lo)));
    match1 = detect_hom_mask_hi & (~(match1 | ((match1 | detect_hom_mask_hi) - detect_hom_mask_lo)));
    match2 = detect_hom_mask_hi & (~(match2 | ((match2 | detect_hom_mask_hi) - detect_hom_mask_lo)));
    // since code10_width >= 4, we can use match0 == 3 (mod 4), match1 == 2
    // (mod 4), match2 == 1 (mod 4) representation.
    fvals_bits = (match0 >> (code10_width - 4)) | (match1 >> (code10_width - 3)) | (match2 >> (code10_width - 2));
    if (fvals_widx == fvals_word_ct_m1) {
      fvals_bits = bzhi_max(fvals_bits, ModNz(rare10_ct << code10_logwidth, kBitsPerWord));
    }
    if (!fvals_bits) {
      continue;
    }
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - code10_logwidth)]);
    if (!sample_include) {
      do {
        const uintptr_t bit_idx = ctzw(fvals_bits);
        const uint32_t sample_uidx = cur_deltalist_base[bit_idx >> code10_logwidth];
        target_genovec[sample_uidx / kBitsPerWordD2] ^= (bit_idx & 3) << (2 * (sample_uidx % kBitsPerWordD2));
        fvals_bits &= fvals_bits - 1;
      } while (fvals_bits);
    } else {
      do {
        const uintptr_t bit_idx = ctzw(fvals_bits);
        const uint32_t sample_uidx = cur_deltalist_base[bit_idx >> code10_logwidth];
        const uint32_t sample_widx = sample_uidx / kBitsPerWord;
        const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
        const uintptr_t sample_include_word = sample_include[sample_widx];
        if (sample_include_word & lowbit) {
          const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
          target_genovec[sample_idx / kBitsPerWordD2] ^= (bit_idx & 3) << (2 * (sample_idx % kBitsPerWordD2));
        }
        fvals_bits &= fvals_bits - 1;
      } while (fvals_bits);
    }
  }
}

PglErr PgrGet2(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx0, uint32_t allele_idx1, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  assert(allele_idx0 != allele_idx1);
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t multiallelic_hc_present = VrtypeMultiallelicHc(vrtype);
  if (!multiallelic_hc_present) {
    if ((allele_idx0 > 1) && (allele_idx1 > 1)) {
      // Trivial all-missing case.
      SetAllBits(2 * sample_ct, genovec);
      return kPglRetSuccess;
    }
    PglErr reterr = ReadGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genovec);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (allele_idx0 < allele_idx1) {
      Rotate2(allele_idx0, allele_idx1, sample_ct, genovec);
      return kPglRetSuccess;
    }
    if (allele_idx0 == 1) {
      GenovecInvertUnsafe(sample_ct, genovec);
      return kPglRetSuccess;
    }
    if (!allele_idx1) {
      GenovecNonzeroToMissingThenInvertUnsafe(sample_ct, genovec);
      return kPglRetSuccess;
    }
    GenovecNontwoToMissingUnsafe(sample_ct, genovec);
    return kPglRetSuccess;
  }
  uintptr_t* raw_genovec = pgrp->workspace_vec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRawGenovec(subsetting_required, vidx, pgrp, &fread_ptr, &fread_end, raw_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  ZeroTrailingQuaters(raw_sample_ct, raw_genovec);

  uint32_t invert = 0;
  if (allele_idx0 > allele_idx1) {
    const uint32_t swap = allele_idx0;
    allele_idx0 = allele_idx1;
    allele_idx1 = swap;
    invert = 1;
  }
  if (allele_idx0 > 1) {
    SetAllBits(2 * sample_ct, genovec);
  } else {
    CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, genovec);
    Rotate2(allele_idx0, allele_idx1, sample_ct, genovec);
  }
  const uint32_t aux1_first_byte = *fread_ptr++;
  const uint32_t aux1a_mode = aux1_first_byte & 15;
  const uint32_t aux1b_mode = aux1_first_byte >> 4;
  uint32_t raw_01_ct = 0;
  uint32_t raw_10_ct = 0;
  if ((!aux1a_mode) || (!aux1b_mode)) {
    GenovecCount12Unsafe(raw_genovec, raw_sample_ct, &raw_01_ct, &raw_10_ct);
  }
  if (!subsetting_required) {
    sample_include = nullptr;
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx];
  uint32_t* deltalist_workspace = pgrp->workspace_difflist_sample_ids;
  if (!allele_idx0) {
    // Two cases:
    // - If allele_idx == 1, convert all aux1a entries from 01 to 11.
    // - Otherwise, for each matching aux1a entry, convert from 11 to 01.
    reterr = GenovecAux1aUpdate(fread_end, sample_include, sample_include_cumulative_popcounts, raw_genovec, aux1a_mode, raw_sample_ct, allele_ct, allele_idx1, 2, raw_01_ct, &fread_ptr, genovec, deltalist_workspace);
  } else {
    reterr = SkipAux1a(fread_end, aux1a_mode, raw_sample_ct, allele_ct, raw_01_ct, &fread_ptr);
  }
  if (unlikely(reterr)) {
    return reterr;
  }
  reterr = GenovecAux1bUpdate2(fread_end, sample_include, sample_include_cumulative_popcounts, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, allele_idx0, allele_idx1, raw_10_ct, &fread_ptr, genovec, deltalist_workspace);
  if (unlikely(reterr)) {
    return reterr;
  }
  if (invert) {
    GenovecInvertUnsafe(sample_ct, genovec);
  }
  return kPglRetSuccess;
}

void PreinitPgv(PgenVariant* pgvp) {
  pgvp->genovec = nullptr;
  pgvp->patch_01_set = nullptr;
  pgvp->patch_01_vals = nullptr;
  pgvp->patch_10_set = nullptr;
  pgvp->patch_10_vals = nullptr;
  pgvp->phasepresent = nullptr;
  pgvp->phaseinfo = nullptr;
  pgvp->dosage_present = nullptr;
  pgvp->dosage_main = nullptr;
  pgvp->multidosage_present = nullptr;
  pgvp->multidosage_cts = nullptr;
  pgvp->multidosage_codes = nullptr;
  pgvp->multidosage_vals = nullptr;
  pgvp->dphase_present = nullptr;
  pgvp->dphase_delta = nullptr;
  pgvp->multidphase_present = nullptr;
  pgvp->multidphase_cts = nullptr;
  pgvp->multidphase_codes = nullptr;
  pgvp->multidphase_delta = nullptr;

  pgvp->patch_01_ct = 0;
  pgvp->patch_10_ct = 0;
  pgvp->phasepresent_ct = 0;
  pgvp->dosage_ct = 0;
  pgvp->multidosage_sample_ct = 0;
  pgvp->dphase_ct = 0;
  pgvp->multidphase_sample_ct = 0;
}

// similar to ParseAndSaveDifflist()
PglErr ParseAndSaveDeltalistAsBitarr(const unsigned char* fread_end, uint32_t raw_sample_ct, const unsigned char** fread_pp, uintptr_t* deltalist_include, uint32_t* deltalist_len_ptr) {
  const unsigned char* group_info_iter;
  PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, deltalist_len_ptr);
  const uint32_t deltalist_len = *deltalist_len_ptr;
  if (reterr || (!deltalist_len)) {
    return reterr;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(raw_sample_ct);
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t group_idx_last = (deltalist_len - 1) / kPglDifflistGroupSize;
  ZeroWArr(raw_sample_ctl, deltalist_include);
  uint32_t group_len_m1 = kPglDifflistGroupSize - 1;
  for (uint32_t group_idx = 0; ; ++group_idx) {
    if (group_idx >= group_idx_last) {
      if (group_idx > group_idx_last) {
        return kPglRetSuccess;
      }
      group_len_m1 &= deltalist_len - 1;
    }
    uintptr_t raw_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
    group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    for (uint32_t raw_deltalist_idx_lowbits = 0; ; ++raw_deltalist_idx_lowbits) {
      // always check, otherwise we may scribble over arbitrary memory
      if (unlikely(raw_sample_idx >= raw_sample_ct)) {
        return kPglRetMalformedInput;
      }
      SetBit(raw_sample_idx, deltalist_include);
      if (raw_deltalist_idx_lowbits == group_len_m1) {
        break;
      }
      raw_sample_idx += GetVint31(fread_end, fread_pp);
    }
  }
}

// These functions do not overread, but may write extra bytes up to the word
// boundary.
void Expand2bitTo8(const void* __restrict bytearr, uint32_t input_quater_ct, uint32_t incr, uintptr_t* __restrict dst) {
  const unsigned char* src_iter = S_CAST(const unsigned char*, bytearr);
  const uint32_t input_byte_ct = DivUp(input_quater_ct, 4);
#ifdef __arm__
#  error "Unaligned accesses in Expand2bitTo8()."
#endif
#ifdef __LP64__
  const uint32_t input_vec_ct = input_byte_ct / kBytesPerVec;
  unsigned char* dst_iter = R_CAST(unsigned char*, dst);
  if (input_vec_ct) {
    const VecW mincr = R_CAST(VecW, vecuc_set1(incr));
    const VecW m03 = VCONST_W(kMask0303);
    for (uint32_t vec_idx = 0; vec_idx != input_vec_ct; ++vec_idx) {
      VecW cur_vec = vecw_loadu(src_iter);
      src_iter = &(src_iter[kBytesPerVec]);
#  ifdef USE_AVX2
      // (todo: benchmark against just reading 8 bytes at a time and
      // broadcasting.)
      // midswapped_vec contains {0-1-2-3, 4-5-6-7, ..., 12-13-14-15,
      //                          32-33-34-35, ..., 44-45-46-47,
      //                          16-17-18-19, ..., 28-29-30-31,
      //                          48-49-50-51, ..., 60-61-62-63,
      //                          64-65-66-67, ..., 76-77-78-79,
      //                          96-97-98-99, ..., 108-109-110-111,
      //                          80-81-82-83, ..., 92-93-94-95,
      //                          112-113-114-115, ..., 124-125-126-127}
      // 0xd8: {0, 2, 1, 3}
      const __m256i midswapped_vec = _mm256_shuffle_epi32(R_CAST(__m256i, cur_vec), 0xd8);
      // This operation is also used in FillInterleavedMaskVec().
      // cur_vec now contains {0-1-2-3, 4-5-6-7, 8-9-10-11, 12-13-14-15,
      //                       32-33-34-35, ..., 44-45-46-47,
      //                       64-65-66-67, ..., 76-77-78-79,
      //                       96-97-98-99, ..., 108-109-110-111,
      //                       16-17-18-19, ..., 28-29-30-31,
      //                       48-49-50-51, ..., 60-61-62-63,
      //                       80-81-82-83, ..., 92-93-94-95,
      //                       112-113-114-115, ..., 124-125-126-127}
      cur_vec = vecw_permute0xd8_if_avx2(R_CAST(VecW, midswapped_vec));
#  endif
      // AVX2:
      //   vec_even contains {0-1, 4-5, 8-9, 12-13, 32-33, ..., 44-45,
      //                      64-65, ..., 76-77, 96-97, ..., 108-109,
      //                      16-17, ..., 28-29, 48-49, ..., 60-61,
      //                      80-81, ..., 92-93, 112-113, ..., 124-125}
      //   vec_odd contains {2-3, 6-7, 10-11, 14-15, 34-35, ..., 46-47,
      //                     66-67, ..., 78-79, 98-99, ..., 110-111,
      //                     18-19, ..., 30-31, 50-51, ..., 62-63,
      //                     82-83, ..., 94-95, 114-115, ..., 126-127}
      // SSE2:
      //   vec_even contains {0-1, 4-5, 8-9, ..., 60-61}
      //   vec_odd contains {2-3, 6-7, 10-11, ..., 62-63}
      const VecW vec_even = cur_vec;
      const VecW vec_odd = vecw_srli(cur_vec, 4);

      // AVX2:
      //   vec01 contains {0-1, 2-3, 4-5, ..., 14-15, 32-33, ..., 46-47,
      //                   16-17, ..., 30-31, 48-49, ..., 62-63}
      //   vec23 contains {64-65, 66-67, ..., 78-79, 96-97, ..., 110-111,
      //                   80-81, ..., 94-95, 112-113, ..., 126-127}
      // SSE2:
      //   vec01 contains {0-1, 2-3, 4-5, 6-7, ..., 30-31}
      //   vec23 contains {32-33, 34-35, 36-37, 38-39, ..., 62-63}
      const VecW vec01 = vecw_unpacklo8(vec_even, vec_odd);
      const VecW vec23 = vecw_unpackhi8(vec_even, vec_odd);

      // AVX2:
      //   vec01_even contains {0, 2, 4, ..., 14, 32, 34, ..., 46,
      //                        16, 18, ..., 30, 48, 50, ..., 62}
      //   vec01_odd contains {1, 3, 5, ..., 15, 33, 35, ..., 47,
      //                       17, 19, ..., 31, 49, 51, ..., 63}
      // SSE2:
      //   vec01_even contains {0, 2, 4, 6, ..., 30}
      //   vec01_odd contains {1, 3, 5, 7, ..., 31}
      const VecW vec01_even = vec01 & m03;
      const VecW vec01_odd = vecw_srli(vec01, 2) & m03;

      // AVX2:
      //   vecw_unpacklo8() contains {0, 1, ..., 15, 16, ..., 31}
      //   vecw_unpachhi8() contains {32, 33, ..., 47, 48, ..., 63}
      // SSE2:
      //   vecw_unpacklo8() contains {0, 1, ..., 15}
      //   vecw_unpachhi8() contains {16, 17, ..., 31}
      vecw_storeu(dst_iter, mincr + vecw_unpacklo8(vec01_even, vec01_odd));
      dst_iter = &(dst_iter[kBytesPerVec]);
      vecw_storeu(dst_iter, mincr + vecw_unpackhi8(vec01_even, vec01_odd));
      dst_iter = &(dst_iter[kBytesPerVec]);
      const VecW vec23_odd = vecw_srli(vec23, 2) & m03;
      const VecW vec23_even = vec23 & m03;
      vecw_storeu(dst_iter, mincr + vecw_unpacklo8(vec23_even, vec23_odd));
      dst_iter = &(dst_iter[kBytesPerVec]);
      vecw_storeu(dst_iter, mincr + vecw_unpackhi8(vec23_even, vec23_odd));
      dst_iter = &(dst_iter[kBytesPerVec]);
    }
  }
  const uint32_t remainder = input_byte_ct % kBytesPerVec;
  if (remainder) {
    const uint32_t full_qw_ct = remainder / sizeof(Quarterword);
    const Quarterword* src_alias = R_CAST(const Quarterword*, src_iter);
    const uintptr_t incr_word = kMask0101 * incr;
    uintptr_t* dstw = R_CAST(uintptr_t*, dst_iter);
    for (uint32_t uii = 0; uii != full_qw_ct; ++uii) {
      const uintptr_t cur_2byte = src_alias[uii];
      dstw[uii] = incr_word + Unpack0303(cur_2byte);
    }
    if (input_byte_ct % 2) {
      uintptr_t cur_byte = src_iter[remainder - 1];
#  ifdef USE_AVX2
      cur_byte = _pdep_u64(cur_byte, kMask0303);
#  else
      cur_byte = cur_byte | (cur_byte << 12);
      cur_byte = (cur_byte | (cur_byte << 6)) & kMask0303;
#  endif
      dstw[full_qw_ct] = incr_word + cur_byte;
    }
  }
#else  // !__LP64__
  const Quarterword* src_alias = R_CAST(const Quarterword*, src_iter);
  const uintptr_t incr_word = kMask0101 * incr;
  uintptr_t* dstw = R_CAST(uintptr_t*, dst);
  for (uint32_t uii = 0; uii != input_byte_ct; ++uii) {
    const uintptr_t cur_2byte = src_alias[uii];
    dstw[uii] = incr_word + Unpack0303(cur_2byte);
  }
#endif
}

void Expand4bitTo8(const void* __restrict bytearr, uint32_t input_nibble_ct, uint32_t incr, uintptr_t* __restrict dst) {
  const unsigned char* src_iter = R_CAST(const unsigned char*, bytearr);
  const uint32_t input_byte_ct = DivUp(input_nibble_ct, 2);
#ifdef __LP64__
  const uint32_t input_vec_ct = input_byte_ct / kBytesPerVec;
  unsigned char* dst_iter = R_CAST(unsigned char*, dst);
  if (input_vec_ct) {
    const VecW mincr = R_CAST(VecW, vecuc_set1(incr));
    const VecW m4 = VCONST_W(kMask0F0F);
    for (uint32_t vec_idx = 0; vec_idx != input_vec_ct; ++vec_idx) {
      VecW cur_vec = vecw_loadu(src_iter);
      src_iter = &(src_iter[kBytesPerVec]);
      cur_vec = vecw_permute0xd8_if_avx2(cur_vec);
      // AVX2:
      //   vec_even contains {0, 2, 4, ..., 14, 32, 34, ..., 46,
      //                      16, 18, ..., 30, 48, ... 62}
      //   vec_odd contains {1, 3, 5, ..., 15, 33, 35, ..., 47,
      //                     17, 19, ..., 31, 49, ..., 63}
      // SSE2:
      //   vec_even contains {0, 2, 4, ..., 30}
      //   vec_odd contains {1, 3, 5, ..., 31}
      const VecW vec_even = cur_vec & m4;
      const VecW vec_odd = vecw_srli(cur_vec, 4) & m4;

      // AVX2:
      //   vec_lo contains {0, 1, ..., 31}
      //   vec_hi contains {32, 33, ..., 63}
      // SSE2:
      //   vec_lo contains {0, 1, 2, ..., 15}
      //   vec_hi contains {16, 17, 18, ..., 31}
      const VecW vec_lo = vecw_unpacklo8(vec_even, vec_odd);
      const VecW vec_hi = vecw_unpackhi8(vec_even, vec_odd);
      vecw_storeu(dst_iter, mincr + vec_lo);
      dst_iter = &(dst_iter[kBytesPerVec]);
      vecw_storeu(dst_iter, mincr + vec_hi);
      dst_iter = &(dst_iter[kBytesPerVec]);
    }
  }
  const uint32_t remainder = input_byte_ct % kBytesPerVec;
  if (remainder) {
    const Halfword* src_alias = R_CAST(const Halfword*, src_iter);
    uintptr_t incr_word = kMask0101 * incr;
    const uint32_t hw_ct_m1 = (remainder - 1) / sizeof(Halfword);
    uintptr_t* dstw = R_CAST(uintptr_t*, dst_iter);
    for (uint32_t hwidx = 0; ; ++hwidx) {
      uint32_t cur_4byte;
      if (hwidx >= hw_ct_m1) {
        if (hwidx > hw_ct_m1) {
          break;
        }
        cur_4byte = SubU32Load(&(src_alias[hwidx]), ModNz(remainder, 4));
      } else {
        cur_4byte = src_alias[hwidx];
      }
      dstw[hwidx] = incr_word + Unpack0F0F(cur_4byte);
    }
  }
#else
  unsigned char* dst_iter = R_CAST(unsigned char*, dst);
  for (uint32_t uii = 0; uii < input_byte_ct; ++uii) {
    uint32_t cur_byte = src_iter[uii];
    *dst_iter++ = (cur_byte & 15) + incr;
    *dst_iter++ = (cur_byte >> 4) + incr;
  }
#endif
}

static_assert(sizeof(AlleleCode) == 1, "GetAux1aCodes() must be updated.");
PglErr GetAux1aCodes(const unsigned char* fread_end, uint32_t rare01_ct, uint32_t allele_ct, const unsigned char** fread_pp, AlleleCode* __restrict patch_01_vals) {
  if (allele_ct == 3) {
    memset(patch_01_vals, 2, rare01_ct);
    return kPglRetSuccess;
  }
  const unsigned char* patch_01_fvals = *fread_pp;
  if (allele_ct == 4) {
    const uint32_t patch_01_fvals_byte_ct = DivUp(rare01_ct, CHAR_BIT);
    if (PtrAddCk(fread_end, patch_01_fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    Expand1bitTo8(patch_01_fvals, rare01_ct, 2, R_CAST(uintptr_t*, patch_01_vals));
    return kPglRetSuccess;
  }
  if (allele_ct < 7) {
    const uint32_t patch_01_fvals_byte_ct = DivUp(rare01_ct, 4);
    if (PtrAddCk(fread_end, patch_01_fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    Expand2bitTo8(patch_01_fvals, rare01_ct, 2, R_CAST(uintptr_t*, patch_01_vals));
    return kPglRetSuccess;
  }
  if (allele_ct < 19) {
    const uint32_t patch_01_fvals_byte_ct = DivUp(rare01_ct, 2);
    if (PtrAddCk(fread_end, patch_01_fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    Expand4bitTo8(patch_01_fvals, rare01_ct, 2, R_CAST(uintptr_t*, patch_01_vals));
    return kPglRetSuccess;
  }
  if (PtrAddCk(fread_end, rare01_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  // todo: verify the compiler recognizes this
  for (uint32_t uii = 0; uii < rare01_ct; ++uii) {
    patch_01_vals[uii] = patch_01_fvals[uii] + 2;
  }
  return kPglRetSuccess;
}

// Assumes aux1a_mode != 15.
PglErr ExportAux1a(const unsigned char* fread_end, const uintptr_t* __restrict raw_genovec, uint32_t aux1a_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t raw_01_ct, const unsigned char** fread_pp, uintptr_t* __restrict patch_01_set, AlleleCode* __restrict patch_01_vals, uint32_t* __restrict rare01_ctp) {
  uint32_t rare01_ct;
  if (!aux1a_mode) {
    const unsigned char* patch_01_fset = *fread_pp;
    const uint32_t fset_byte_ct = DivUp(raw_01_ct, CHAR_BIT);
    if (PtrAddCk(fread_end, fset_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    rare01_ct = PopcountBytes(patch_01_fset, fset_byte_ct);
    ExpandBytearrFromGenovec(patch_01_fset, raw_genovec, kMask5555, QuaterCtToWordCt(raw_sample_ct), raw_01_ct, 0, patch_01_set);
  } else {
    if (unlikely(ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, fread_pp, patch_01_set, &rare01_ct))) {
      return kPglRetMalformedInput;
    }
  }
  *rare01_ctp = rare01_ct;
  return GetAux1aCodes(fread_end, rare01_ct, allele_ct, fread_pp, patch_01_vals);
}

PglErr ExportAux1aProperSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict raw_genovec, uint32_t aux1a_mode, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t allele_ct, uint32_t raw_01_ct, const unsigned char** fread_pp, uintptr_t* __restrict dst_01_set, AlleleCode* __restrict dst_01_vals, uint32_t* __restrict rare01_ctp, uint32_t* __restrict deltalist_workspace) {
  const uint32_t allele_code_width = GetAux1aWidth(allele_ct);
  const uintptr_t allele_code_mask = (1U << allele_code_width) - 1;
  memset(dst_01_set, 0, BitCtToWordCt(sample_ct) * sizeof(intptr_t));
  AlleleCode* dst_01_vals_iter = dst_01_vals;
  if (!aux1a_mode) {
#ifdef __arm__
#  error "Unaligned accesses in ExportAux1aProperSubset()."
#endif
    // similar to GenovecAux1aUpdate()
    const uintptr_t* patch_01_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fset_byte_ct = DivUp(raw_01_ct, CHAR_BIT);
    const uint32_t rare01_ct = PopcountBytes(patch_01_fsetw, fset_byte_ct);
    if (PtrAddCk(fread_end, fset_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const uintptr_t* patch_01_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_hets = Word01(raw_genovec[0]);
    uint32_t loop_len = kBitsPerWord;
    const uintptr_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t rare01_lowbits = kBitsPerWord;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          break;
        }
        fset_bits = SubwordLoad(&(patch_01_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_01_ct, kBitsPerWord);
      } else {
        fset_bits = patch_01_fsetw[fset_widx];
      }
      if (allele_ct == 3) {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_hets) {
            cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_hets) / 2;
            if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
              const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
              SetBit(sample_idx, dst_01_set);
              *dst_01_vals_iter++ = 2;
            }
          }
          cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
          fset_bits = fset_bits >> 1;
        }
      } else {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_hets) {
            cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare01_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_01_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_01_fvalsw[fvals_widx];
              }
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare01_lowbits = 0;
            }
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_hets) / 2;
            if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
              const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
              SetBit(sample_idx, dst_01_set);
              *dst_01_vals_iter++ = 2 + ((fvals_bits >> rare01_lowbits) & allele_code_mask);
            }
            rare01_lowbits += allele_code_width;
          }
          cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
          fset_bits = fset_bits >> 1;
        }
      }
    }
    *rare01_ctp = dst_01_vals_iter - dst_01_vals;
    return kPglRetSuccess;
  }
  // aux1a_mode == 1
  uint32_t rare01_ct;
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare01_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_01_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uintptr_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  if (allele_ct == 3) {
    for (uint32_t rare01_idx = 0; rare01_idx != rare01_ct; ++rare01_idx) {
      const uint32_t sample_uidx = deltalist_workspace[rare01_idx];
      // could wrap this boilerplate
      const uint32_t sample_widx = sample_uidx / kBitsPerWord;
      const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
      const uintptr_t sample_include_word = sample_include[sample_widx];
      if (sample_include_word & lowbit) {
        const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
        SetBit(sample_idx, dst_01_set);
        *dst_01_vals_iter++ = 2;
      }
    }
    *rare01_ctp = dst_01_vals_iter - dst_01_vals;
    return kPglRetSuccess;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  const uint32_t allele_code_logwidth = ctzu32(allele_code_width);
  uint32_t loop_len = kBitsPerWord >> allele_code_logwidth;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        break;
      }
      fvals_bits = SubwordLoad(&(patch_01_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
      loop_len = 1 + ((rare01_ct - 1) & (loop_len - 1));
    } else {
      fvals_bits = patch_01_fvalsw[fvals_widx];
    }
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - allele_code_logwidth)]);
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uint32_t sample_uidx = cur_deltalist_base[uii];
      const uint32_t sample_widx = sample_uidx / kBitsPerWord;
      const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
      const uintptr_t sample_include_word = sample_include[sample_widx];
      if (sample_include_word & lowbit) {
        const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
        SetBit(sample_idx, dst_01_set);
        *dst_01_vals_iter++ = 2 + ((fvals_bits >> (uii << allele_code_logwidth)) & allele_code_mask);
      }
    }
  }
  *rare01_ctp = dst_01_vals_iter - dst_01_vals;
  return kPglRetSuccess;
}

static_assert(sizeof(AlleleCode) == 1, "GetAux1bCodes() must be updated.");
PglErr GetAux1bCodes(const unsigned char* fread_end, uint32_t rare10_ct, uint32_t allele_ct, const unsigned char** fread_pp, AlleleCode* __restrict patch_10_vals) {
  const unsigned char* patch_10_fvals = *fread_pp;
  if (allele_ct == 3) {
    // 1 bit, distinguishes between 0x0201 and 0x0202
    const uint32_t patch_10_fvals_byte_ct = DivUp(rare10_ct, CHAR_BIT);
    if (PtrAddCk(fread_end, patch_10_fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    Expand1bitTo16(patch_10_fvals, rare10_ct, 0x0201, R_CAST(uintptr_t*, patch_10_vals));
    return kPglRetSuccess;
  }
  const uint32_t rare10_ct_x2 = rare10_ct * 2;
  if (allele_ct < 6) {
    // 2+2 bits, add 1
    const uint32_t patch_10_fvals_byte_ct = DivUp(rare10_ct, 2);
    if (PtrAddCk(fread_end, patch_10_fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    Expand2bitTo8(patch_10_fvals, rare10_ct_x2, 1, R_CAST(uintptr_t*, patch_10_vals));
    return kPglRetSuccess;
  }
  if (allele_ct < 18) {
    // 4+4 bits
    if (PtrAddCk(fread_end, rare10_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    Expand4bitTo8(patch_10_fvals, rare10_ct_x2, 1, R_CAST(uintptr_t*, patch_10_vals));
    return kPglRetSuccess;
  }
  if (PtrAddCk(fread_end, rare10_ct_x2, fread_pp)) {
    return kPglRetMalformedInput;
  }
  // todo: verify the compiler recognizes this
  for (uint32_t uii = 0; uii < rare10_ct_x2; ++uii) {
    patch_10_vals[uii] = patch_10_fvals[uii] + 1;
  }
  return kPglRetSuccess;
}

// Assumes aux1b_mode != 15.
PglErr ExportAux1b(const unsigned char* fread_end, const uintptr_t* __restrict raw_genovec, uint32_t aux1b_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t raw_10_ct, const unsigned char** fread_pp, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals, uint32_t* __restrict rare10_ctp) {
  uint32_t rare10_ct;
  if (!aux1b_mode) {
    const unsigned char* patch_10_fset = *fread_pp;
    const uint32_t fset_byte_ct = DivUp(raw_10_ct, CHAR_BIT);
    if (PtrAddCk(fread_end, fset_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    rare10_ct = PopcountBytes(patch_10_fset, fset_byte_ct);
    ExpandBytearrFromGenovec(patch_10_fset, raw_genovec, kMaskAAAA, QuaterCtToWordCt(raw_sample_ct), raw_10_ct, 0, patch_10_set);
  } else {
    if (unlikely(ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, fread_pp, patch_10_set, &rare10_ct))) {
      return kPglRetMalformedInput;
    }
  }
  *rare10_ctp = rare10_ct;
  return GetAux1bCodes(fread_end, rare10_ct, allele_ct, fread_pp, patch_10_vals);
}

PglErr ExportAux1bProperSubset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uintptr_t* __restrict raw_genovec, uint32_t aux1b_mode, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t allele_ct, uint32_t raw_10_ct, const unsigned char** fread_pp, uintptr_t* __restrict dst_10_set, AlleleCode* __restrict dst_10_vals, uint32_t* __restrict rare10_ctp, uint32_t* __restrict deltalist_workspace) {
  uintptr_t detect_hom_mask_lo;  // unused
  const uint32_t allele_code_logwidth = GetAux1bConsts(allele_ct, &detect_hom_mask_lo);
  const uint32_t allele_code_width = 1U << allele_code_logwidth;
  const uintptr_t allele_code_mask = (1U << allele_code_width) - 1;
  const uint32_t code10_logwidth = allele_code_logwidth + (allele_code_logwidth != 0);
  const uint32_t code10_width = 1U << code10_logwidth;
  memset(dst_10_set, 0, BitCtToWordCt(sample_ct) * sizeof(intptr_t));
  AlleleCode* dst_10_vals_iter = dst_10_vals;
  if (!aux1b_mode) {
#ifdef __arm__
#  error "Unaligned accesses in ExportAux1bProperSubset()."
#endif
    const uintptr_t* patch_10_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fset_byte_ct = DivUp(raw_10_ct, CHAR_BIT);
    const uint32_t rare10_ct = PopcountBytes(patch_10_fsetw, fset_byte_ct);
    if (PtrAddCk(fread_end, fset_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_xys = Word10(raw_genovec[0]);
    uint32_t loop_len = kBitsPerWord;
    const uintptr_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) * code10_width, 8);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t rare10_lowbits = kBitsPerWord;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          break;
        }
        fset_bits = SubwordLoad(&(patch_10_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_10_ct, kBitsPerWord);
      } else {
        fset_bits = patch_10_fsetw[fset_widx];
      }
      if (allele_ct == 3) {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_xys) {
            cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare10_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_10_fvalsw[fvals_widx];
              }
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare10_lowbits = 0;
            }
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
            if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
              const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
              SetBit(sample_idx, dst_10_set);
              *dst_10_vals_iter++ = 1 + ((fvals_bits >> rare10_lowbits) & 1);
              *dst_10_vals_iter++ = 2;
            }
            ++rare10_lowbits;
          }
          cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
          fset_bits = fset_bits >> 1;
        }
      } else {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_xys) {
            cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare10_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_10_fvalsw[fvals_widx];
              }
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare10_lowbits = 0;
            }
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
            if ((R_CAST(const Halfword*, sample_include)[sample_hwidx]) & (1U << sample_uidx_lowbits)) {
              const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_hwidx * kBitsPerWordD2 + sample_uidx_lowbits);
              SetBit(sample_idx, dst_10_set);
              const uintptr_t cur_code_pair = fvals_bits >> rare10_lowbits;
              const uint32_t cur_code_hi = (cur_code_pair >> allele_code_width) & allele_code_mask;
              const uint32_t cur_code_lo = cur_code_pair & allele_code_mask;
              *dst_10_vals_iter++ = 1 + cur_code_lo;
              *dst_10_vals_iter++ = 1 + cur_code_hi;
            }
            rare10_lowbits += code10_width;
          }
          cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
          fset_bits = fset_bits >> 1;
        }
      }
    }
    *rare10_ctp = S_CAST(uintptr_t, dst_10_vals_iter - dst_10_vals) / 2;
    return kPglRetSuccess;
  }
  // aux1b_mode == 1
  uint32_t rare10_ct;
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare10_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uintptr_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) << code10_logwidth, 8);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  uint32_t loop_len = kBitsPerWord >> code10_logwidth;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        break;
      }
      fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
      loop_len = 1 + ((rare10_ct - 1) & (loop_len - 1));
    } else {
      fvals_bits = patch_10_fvalsw[fvals_widx];
    }
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - code10_logwidth)]);
    if (allele_ct == 3) {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uint32_t sample_uidx = cur_deltalist_base[uii];
        const uint32_t sample_widx = sample_uidx / kBitsPerWord;
        const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
        const uintptr_t sample_include_word = sample_include[sample_widx];
        if (sample_include_word & lowbit) {
          const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
          SetBit(sample_idx, dst_10_set);
          *dst_10_vals_iter++ = 1 + ((fvals_bits >> uii) & 1);
          *dst_10_vals_iter++ = 2;
        }
      }
    } else {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uint32_t sample_uidx = cur_deltalist_base[uii];
        const uint32_t sample_widx = sample_uidx / kBitsPerWord;
        const uintptr_t lowbit = k1LU << (sample_uidx % kBitsPerWord);
        const uintptr_t sample_include_word = sample_include[sample_widx];
        if (sample_include_word & lowbit) {
          const uint32_t sample_idx = sample_include_cumulative_popcounts[sample_widx] + PopcountWord(sample_include_word & (lowbit - 1));
          SetBit(sample_idx, dst_10_set);
          const uintptr_t cur_code_pair = fvals_bits >> (uii << code10_logwidth);
          const uint32_t cur_code_hi = (cur_code_pair >> allele_code_width) & allele_code_mask;
          const uint32_t cur_code_lo = cur_code_pair & allele_code_mask;
          *dst_10_vals_iter++ = 1 + cur_code_lo;
          *dst_10_vals_iter++ = 1 + cur_code_hi;
        }
      }
    }
  }
  *rare10_ctp = S_CAST(uintptr_t, dst_10_vals_iter - dst_10_vals) / 2;
  return kPglRetSuccess;
}

// Assumes sample_ct > 0, multiallelic-hc track is present, and patch_01_ct and
// patch_10_ct are zero-initialized.
PglErr GetMultiallelicCodes(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict all_hets, PgenVariant* pgvp) {
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* raw_genovec = pgrp->workspace_vec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRawGenovec(subsetting_required, vidx, pgrp, &fread_ptr, &fread_end, raw_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, pgvp->genovec);
  ZeroTrailingQuaters(raw_sample_ct, raw_genovec);
  const uint32_t aux1_first_byte = *fread_ptr++;
  const uint32_t aux1a_mode = aux1_first_byte & 15;
  const uint32_t aux1b_mode = aux1_first_byte >> 4;
  uint32_t raw_01_ct = 0;
  uint32_t raw_10_ct = 0;
  if ((!aux1a_mode) || (!aux1b_mode)) {
    GenovecCount12Unsafe(raw_genovec, raw_sample_ct, &raw_01_ct, &raw_10_ct);
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx];
  uint32_t* deltalist_workspace = pgrp->workspace_difflist_sample_ids;
  if (aux1a_mode != 15) {
    if (!subsetting_required) {
      reterr = ExportAux1a(fread_end, raw_genovec, aux1a_mode, raw_sample_ct, allele_ct, raw_01_ct, &fread_ptr, pgvp->patch_01_set, pgvp->patch_01_vals, &(pgvp->patch_01_ct));
    } else {
      reterr = ExportAux1aProperSubset(fread_end, sample_include, sample_include_cumulative_popcounts, raw_genovec, aux1a_mode, raw_sample_ct, sample_ct, allele_ct, raw_01_ct, &fread_ptr, pgvp->patch_01_set, pgvp->patch_01_vals, &(pgvp->patch_01_ct), deltalist_workspace);
    }
    if (unlikely(reterr)) {
      return reterr;
    }
  }
  const unsigned char* aux1b_start = fread_ptr;
  if (aux1b_mode != 15) {
    if (!subsetting_required) {
      reterr = ExportAux1b(fread_end, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &fread_ptr, pgvp->patch_10_set, pgvp->patch_10_vals, &(pgvp->patch_10_ct));
    } else {
      reterr = ExportAux1bProperSubset(fread_end, sample_include, sample_include_cumulative_popcounts, raw_genovec, aux1b_mode, raw_sample_ct, sample_ct, allele_ct, raw_10_ct, &fread_ptr, pgvp->patch_10_set, pgvp->patch_10_vals, &(pgvp->patch_10_ct), deltalist_workspace);
    }
    if (unlikely(reterr)) {
      return reterr;
    }
  }
  if (fread_pp) {
    *fread_pp = fread_ptr;
    *fread_endp = fread_end;
    if (all_hets) {
      PgrDetectGenovecHets(raw_genovec, raw_sample_ct, all_hets);
      if (aux1b_mode != 15) {
        // can merge this with ExportAux1b functions later
        uintptr_t* aux1b_hets = pgrp->workspace_aux1x_present;
        uint32_t aux1b_het_present;
        reterr = GetAux1bHets(fread_end, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &aux1b_start, aux1b_hets, &aux1b_het_present, deltalist_workspace);
        if (unlikely(reterr)) {
          return reterr;
        }
        if (aux1b_het_present) {
          BitvecOr(aux1b_hets, BitCtToWordCt(raw_sample_ct), all_hets);
        }
      }
    }
  }
  return kPglRetSuccess;
}

PglErr PgrGetM(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp) {
  pgvp->patch_01_ct = 0;
  pgvp->patch_10_ct = 0;
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t multiallelic_hc_present = VrtypeMultiallelicHc(vrtype);
  if (!multiallelic_hc_present) {
    return ReadGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, pgvp->genovec);
  }
  return GetMultiallelicCodes(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, nullptr, pgvp);
}

void DetectGenovecHetsHw(const uintptr_t*__restrict genovec, uint32_t raw_sample_ctl2, Halfword* all_hets_hw) {
  // requires trailing bits of genovec to be zeroed out.  does not update last
  // all_hets[] halfword if raw_sample_ctl2 is odd.
  for (uint32_t widx = 0; widx != raw_sample_ctl2; ++widx) {
    const uintptr_t cur_word = genovec[widx];
    uintptr_t ww = (~(cur_word >> 1)) & cur_word;  // low 1, high 0
    all_hets_hw[widx] = PackWordToHalfwordMask5555(ww);
  }
}

void PgrDetectGenovecHetsMultiallelic(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t raw_sample_ct, uintptr_t* __restrict all_hets) {
  const Halfword* patch_10_set_alias = R_CAST(const Halfword*, patch_10_set);
  const AlleleCode* patch_10_vals_iter = patch_10_vals;
  const uint32_t word_ct_m1 = (raw_sample_ct - 1) / kBitsPerWordD2;
  Halfword* all_hets_hw = R_CAST(Halfword*, all_hets);
  for (uint32_t widx = 0; ; ++widx) {
    uintptr_t cur_geno_word;
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        if (widx % 2) {
          all_hets_hw[widx] = 0;
        }
        return;
      }
      const uint32_t final_ct = ModNz(raw_sample_ct, kBitsPerWordD2);
      cur_geno_word = bzhi_max(genovec[widx], 2 * final_ct);
    } else {
      cur_geno_word = genovec[widx];
    }
    uint32_t patch_10_hw = patch_10_set_alias[widx];
    uint32_t cur_hets = Pack01ToHalfword(cur_geno_word);
    while (patch_10_hw) {
      const AlleleCode code1 = *patch_10_vals_iter++;
      const AlleleCode code2 = *patch_10_vals_iter++;
      const uint32_t lowbit = patch_10_hw & (-patch_10_hw);
      if (code1 != code2) {
        cur_hets |= lowbit;
      }
      patch_10_hw ^= lowbit;
    }
    all_hets_hw[widx] = cur_hets;
  }
}

PglErr SkipAux1b(const unsigned char* fread_end, uint32_t aux1b_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t raw_10_ct, const unsigned char** fread_pp) {
  if (aux1b_mode == 15) {
    return kPglRetSuccess;
  }
  uint32_t rare10_ct;
  if (!aux1b_mode) {
    const uint32_t fset_byte_ct = DivUp(raw_10_ct, CHAR_BIT);
    rare10_ct = PopcountBytes(*fread_pp, fset_byte_ct);
    *fread_pp += fset_byte_ct;
  } else {
    const unsigned char* group_info_iter;
    PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, &rare10_ct);
    if (unlikely(reterr)) {
      return reterr;
    }
    reterr = SkipDeltalistIds(fread_end, group_info_iter, rare10_ct, raw_sample_ct, 0, fread_pp);
    if (unlikely(reterr)) {
      return reterr;
    }
  }
  const uint32_t fvals_byte_ct = GetAux1bAlleleEntryByteCt(allele_ct, rare10_ct);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  return kPglRetSuccess;
}

PglErr SkipAux1(const unsigned char* fread_end, const uintptr_t* __restrict raw_genovec, uint32_t raw_sample_ct, uint32_t allele_ct, const unsigned char** fread_pp) {
  const uint32_t aux1_first_byte = **fread_pp;
  (*fread_pp) += 1;
  const uint32_t aux1a_mode = aux1_first_byte & 15;
  const uint32_t aux1b_mode = aux1_first_byte >> 4;
  uint32_t raw_01_ct = 0;
  uint32_t raw_10_ct = 0;
  if ((!aux1a_mode) || (!aux1b_mode)) {
    GenovecCount12Unsafe(raw_genovec, raw_sample_ct, &raw_01_ct, &raw_10_ct);
  }
  PglErr reterr = SkipAux1a(fread_end, aux1a_mode, raw_sample_ct, allele_ct, raw_01_ct, fread_pp);
  if (unlikely(reterr)) {
    return reterr;
  }
  return SkipAux1b(fread_end, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, fread_pp);
}

// sample_include assumed to be nullptr if no subsetting required
// subsetted_10het should only be provided when you explicitly want to exclude
// those phase entries
// set phasepresent == phaseinfo == nullptr if you want to skip the entire
// track; ok for phasepresent_ct_ptr to be nullptr too in that case
// (also see SkipAux2() and GetPhasepresentAndSkipPhaseinfo() below)
PglErr ParseAux2Subset(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uintptr_t* __restrict all_hets, const uintptr_t* __restrict subsetted_10het, uint32_t raw_sample_ct, uint32_t sample_ct, const unsigned char** fread_pp, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr, uintptr_t* __restrict workspace_subset) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t het_ct = PopcountWords(all_hets, raw_sample_ctl);
  if (unlikely(!het_ct)) {
    // there shouldn't be a hphase track at all in this case, het_ct is not
    // computed off a subset
    return kPglRetMalformedInput;
  }
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const unsigned char* aux2_start = *fread_pp;
  if (!(aux2_start[0] & 1)) {
    // phase always present
    if (PtrAddCk(fread_end, 1 + (het_ct / CHAR_BIT), fread_pp)) {
      return kPglRetMalformedInput;
    }
    if (!phaseinfo) {
      // for internal callers which just want to skip aux2
      return kPglRetSuccess;
    }
    if (!sample_include) {
      memcpy(phasepresent, all_hets, raw_sample_ctl * kBytesPerWord);
      ExpandBytearr(aux2_start, all_hets, raw_sample_ctl, het_ct, 1, phaseinfo);
      if (!subsetted_10het) {
        *phasepresent_ct_ptr = het_ct;
        return kPglRetSuccess;
      }
    } else {
      CopyBitarrSubset(all_hets, sample_include, sample_ct, phasepresent);
      if (AllWordsAreZero(phasepresent, sample_ctl)) {
        *phasepresent_ct_ptr = 0;
        // bugfix (7 Dec 2017): clear sample_ctl words here, not raw_sample_ctl
        ZeroWArr(sample_ctl, phaseinfo);
        return kPglRetSuccess;
      }
      ExpandThenSubsetBytearr(aux2_start, all_hets, sample_include, het_ct, sample_ct, 1, phaseinfo);
    }
    *phasepresent_ct_ptr = PopcountWords(phasepresent, sample_ctl);
    return kPglRetSuccess;
  }
  const uint32_t het_ctdl = het_ct / kBitsPerWord;

  // explicit phasepresent
  const uintptr_t* aux2_first_part = R_CAST(const uintptr_t*, aux2_start);
  uintptr_t* aux2_first_part_copy = workspace_subset;
  aux2_first_part_copy[het_ctdl] = 0;
  memcpy(aux2_first_part_copy, aux2_first_part, 1 + (het_ct / CHAR_BIT));
  const uint32_t raw_phasepresent_ct = PopcountWords(aux2_first_part_copy, het_ctdl + 1) - 1;
  if (unlikely(!raw_phasepresent_ct)) {
    // there shouldn't be a hphase track at all in this case
    return kPglRetMalformedInput;
  }
  const unsigned char* aux2_second_part = &(aux2_start[1 + (het_ct / CHAR_BIT)]);
  *fread_pp = aux2_second_part;
  if (PtrAddCk(fread_end, DivUp(raw_phasepresent_ct, CHAR_BIT), fread_pp)) {
    return kPglRetMalformedInput;
  }
  if (!phaseinfo) {
    return kPglRetSuccess;
  }
  if (!sample_include) {
    ExpandBytearrNested(aux2_second_part, aux2_first_part_copy, all_hets, sample_ctl, raw_phasepresent_ct, 1, phasepresent, phaseinfo);
    if (!subsetted_10het) {
      *phasepresent_ct_ptr = raw_phasepresent_ct;
      return kPglRetSuccess;
    }
  } else {
    // could skip if intersection of phasepresent with sample_include is empty,
    // but this function call should be fast enough there anyway?
    ExpandThenSubsetBytearrNested(aux2_second_part, aux2_first_part_copy, all_hets, sample_include, sample_ct, raw_phasepresent_ct, 1, phasepresent, phaseinfo);
  }
  if (subsetted_10het) {
    BitvecInvmask(subsetted_10het, sample_ctl, phasepresent);
  }
  *phasepresent_ct_ptr = PopcountWords(phasepresent, sample_ctl);
  return kPglRetSuccess;
}

PglErr SkipAux2(const unsigned char* fread_end, uint32_t het_ct, const unsigned char** fread_pp, uint32_t* __restrict phasepresent_ctp) {
  const unsigned char* aux2_start = *fread_pp;
  const uint32_t aux2_first_part_byte_ct = 1 + (het_ct / CHAR_BIT);
  if (PtrAddCk(fread_end, aux2_first_part_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  if (!(aux2_start[0] & 1)) {
    if (phasepresent_ctp) {
      *phasepresent_ctp = het_ct;
    }
    return kPglRetSuccess;
  }
  const uint32_t phasepresent_ct = PopcountBytes(aux2_start, aux2_first_part_byte_ct) - 1;
  if (phasepresent_ctp) {
    *phasepresent_ctp = phasepresent_ct;
  }
  if (PtrAddCk(fread_end, DivUp(phasepresent_ct, CHAR_BIT), fread_pp)) {
    return kPglRetMalformedInput;
  }
  return kPglRetSuccess;
}

// If fread_pp/fread_endp are non-null, this always moves fread_ptr to the end
// of aux2.  Set phasepresent/phaseinfo to nullptr when you don't actually care
// about the contents of aux2.
// In multiallelic case, this guarantees phasepresent bits are only set at
// ref/altx hets, not at altx/alty hets.  (We don't currently guarantee this
// for phaseinfo, since popcounts on that array are meaningless.)  Yes, this is
// mildly annoying, but the code would be messier if the ordering of
// multiallelic-hardcall and hardcall-phase info were swapped.
PglErr ReadGenovecHphaseSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* phasepresent_ct_ptr) {
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  if ((!(vrtype & 0x18)) || ((!fread_pp) && (!VrtypeHphase(vrtype)))) {
    *phasepresent_ct_ptr = 0;
    return ReadGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, fread_pp, fread_endp, genovec);
  }
  // Either hphase track is present; or if it's absent, multiallelic track is
  // present and we were asked to advance fread_ptr to the end of aux2.
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* raw_genovec = (subsetting_required || VrtypeMultiallelicHc(vrtype))? pgrp->workspace_vec : genovec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRawGenovec(subsetting_required, vidx, pgrp, &fread_ptr, &fread_end, raw_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  ZeroTrailingQuaters(raw_sample_ct, raw_genovec);
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  if (raw_genovec != genovec) {
    CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, genovec);
    if (!VrtypeHphase(vrtype)) {
      // only possible if multiallelic track present and fread_ptr must be
      // advanced to end of aux2
      *fread_pp = fread_ptr;
      *fread_endp = fread_end;
      return SkipAux1(fread_end, raw_genovec, raw_sample_ct, allele_ct, fread_pp);
    }
  }
  uintptr_t* all_hets = pgrp->workspace_all_hets;
  PgrDetectGenovecHets(raw_genovec, raw_sample_ct, all_hets);
  uintptr_t* subsetted_10het = nullptr;
  if (VrtypeMultiallelicHc(vrtype)) {
    const uint32_t aux1_first_byte = *fread_ptr++;
    const uint32_t aux1a_mode = aux1_first_byte & 15;
    const uint32_t aux1b_mode = aux1_first_byte >> 4;
    uint32_t raw_01_ct = 0;
    uint32_t raw_10_ct = 0;
    if ((!aux1a_mode) || (!aux1b_mode)) {
      GenovecCount12Unsafe(raw_genovec, raw_sample_ct, &raw_01_ct, &raw_10_ct);
    }
    reterr = SkipAux1a(fread_end, aux1a_mode, raw_sample_ct, allele_ct, raw_01_ct, &fread_ptr);
    if (unlikely(reterr)) {
      return reterr;
    }
    // 1. fill workspace_aux1x_present with aux1b
    // 2. clear bit for each hom-altx call in aux1b
    // 3. bitvec-or to set new workspace_all_hets bits
    // 4. if not subsetting, set subsetted_10het := workspace_all_hets
    //    if subsetting, copy-subset to pgrp->workspace_vec and set to that
    //    if AllWordsAreZero, keep as nullptr
    uintptr_t* aux1b_hets = pgrp->workspace_aux1x_present;
    uint32_t* deltalist_workspace = pgrp->workspace_difflist_sample_ids;
    uint32_t aux1b_het_present;
    reterr = GetAux1bHets(fread_end, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &fread_ptr, aux1b_hets, &aux1b_het_present, deltalist_workspace);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (aux1b_het_present) {
      BitvecOr(aux1b_hets, BitCtToWordCt(raw_sample_ct), all_hets);
      if (!subsetting_required) {
        subsetted_10het = aux1b_hets;
      } else {
        // Don't need raw_genovec any more.
        CopyBitarrSubset(aux1b_hets, sample_include, sample_ct, raw_genovec);
        subsetted_10het = raw_genovec;
      }
    }
  }
  reterr = ParseAux2Subset(fread_end, subsetting_required? sample_include : nullptr, all_hets, subsetted_10het, raw_sample_ct, sample_ct, &fread_ptr, phasepresent, phaseinfo, phasepresent_ct_ptr, pgrp->workspace_subset);
  if (fread_pp) {
    *fread_pp = fread_ptr;
    *fread_endp = fread_end;
  }
  return reterr;
}

PglErr PgrGetP(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    *phasepresent_ct_ptr = 0;
    return kPglRetSuccess;
  }
  return ReadGenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genovec, phasepresent, phaseinfo, phasepresent_ct_ptr);
}

// eventually want to return fread_ptr/fread_end, but not relevant until
// multiallelic dosage working
PglErr Get1MP(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_countvec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr) {
  // sample_ct > 0; either allele_idx > 1 or ((allele_idx == 1) &&
  // multiallelic_hc_present)
  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  if (!VrtypeHphase(vrtype)) {
    *phasepresent_ct_ptr = 0;
    return PgrGet1(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, allele_countvec);
  }
  uintptr_t* all_hets = pgrp->workspace_all_hets;
  uintptr_t* subsetted_10het = nullptr;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = Get1Multiallelic(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, &fread_ptr, &fread_end, all_hets, allele_countvec, &subsetted_10het);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  reterr = ParseAux2Subset(fread_end, (sample_ct != raw_sample_ct)? sample_include : nullptr, all_hets, subsetted_10het, raw_sample_ct, sample_ct, &fread_ptr, phasepresent, phaseinfo, phasepresent_ct_ptr, pgrp->workspace_subset);
  // bugfix (7 Sep 2018): Need to postprocess phasepresent when collapsing
  // multiple alleles.
  if (reterr || (!(*phasepresent_ct_ptr))) {
    return reterr;
  }

  // Might want to make this its own function.
  const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
  Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
  for (uint32_t hwidx = 0; hwidx != sample_ctl2; ++hwidx) {
    phasepresent_alias[hwidx] &= Pack01ToHalfword(allele_countvec[hwidx]);
  }
  *phasepresent_ct_ptr = PopcountWords(phasepresent, BitCtToWordCt(sample_ct));

  return kPglRetSuccess;
}

PglErr PgrGet1P(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_countvec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr) {
  if (!sample_ct) {
    *phasepresent_ct_ptr = 0;
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t multiallelic_hc_present = VrtypeMultiallelicHc(vrtype);
  if ((!allele_idx) || ((allele_idx == 1) && (!multiallelic_hc_present))) {
    PglErr reterr = ReadGenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, allele_countvec, phasepresent, phaseinfo, phasepresent_ct_ptr);
    if (allele_idx) {
      GenovecInvertUnsafe(sample_ct, allele_countvec);
      if (*phasepresent_ct_ptr) {
        BitvecInvert(BitCtToWordCt(sample_ct), phaseinfo);
      }
    }
    return reterr;
  }
  return Get1MP(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, allele_countvec, phasepresent, phaseinfo, phasepresent_ct_ptr);
}

PglErr PgrGetInv1P(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_invcountvec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr) {
  if (!sample_ct) {
    *phasepresent_ct_ptr = 0;
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t multiallelic_hc_present = VrtypeMultiallelicHc(vrtype);
  if ((!allele_idx) || ((allele_idx == 1) && (!multiallelic_hc_present))) {
    PglErr reterr = ReadGenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, allele_invcountvec, phasepresent, phaseinfo, phasepresent_ct_ptr);
    if (!allele_idx) {
      GenovecInvertUnsafe(sample_ct, allele_invcountvec);
      if (*phasepresent_ct_ptr) {
        BitvecInvert(BitCtToWordCt(sample_ct), phaseinfo);
      }
    }
    return reterr;
  }
  PglErr reterr = Get1MP(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, allele_invcountvec, phasepresent, phaseinfo, phasepresent_ct_ptr);
  if (unlikely(reterr)) {
    return reterr;
  }
  GenovecInvertUnsafe(sample_ct, allele_invcountvec);
  if (*phasepresent_ct_ptr) {
    BitvecInvert(BitCtToWordCt(sample_ct), phaseinfo);
  }
  return kPglRetSuccess;
}

PglErr PgrGet2P(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx0, uint32_t allele_idx1, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr) {
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  if (!VrtypeHphase(vrtype)) {
    *phasepresent_ct_ptr = 0;
    return PgrGet2(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx0, allele_idx1, pgrp, genovec);
  }
  if (!sample_ct) {
    *phasepresent_ct_ptr = 0;
    return kPglRetSuccess;
  }
  if (allele_idx0 + allele_idx1 == 1) {
    PglErr reterr = ReadGenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genovec, phasepresent, phaseinfo, phasepresent_ct_ptr);
    if (allele_idx0) {
      GenovecInvertUnsafe(sample_ct, genovec);
      if (*phasepresent_ct_ptr) {
        BitvecInvert(BitCtToWordCt(sample_ct), phaseinfo);
      }
    }
    return reterr;
  }
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* raw_genovec = pgrp->workspace_vec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRawGenovec(subsetting_required, vidx, pgrp, &fread_ptr, &fread_end, raw_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  ZeroTrailingQuaters(raw_sample_ct, raw_genovec);

  uint32_t invert = 0;
  if (allele_idx0 > allele_idx1) {
    const uint32_t swap = allele_idx0;
    allele_idx0 = allele_idx1;
    allele_idx1 = swap;
    invert = 1;
  }
  if (allele_idx0 > 1) {
    SetAllBits(2 * sample_ct, genovec);
  } else {
    CopyQuaterarrNonemptySubset(raw_genovec, sample_include, raw_sample_ct, sample_ct, genovec);
    // allele_idx1 > 1 guaranteed
    if (!allele_idx0) {
      GenovecNonzeroToMissingUnsafe(sample_ct, genovec);
    } else {
      GenovecInvertThenNonzeroToMissingUnsafe(sample_ct, genovec);
    }
  }
  uintptr_t* all_hets = pgrp->workspace_all_hets;
  PgrDetectGenovecHets(raw_genovec, raw_sample_ct, all_hets);
  uintptr_t* subsetted_10het = nullptr;
  if (!subsetting_required) {
    sample_include = nullptr;
  }

  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx];
  if (VrtypeMultiallelicHc(vrtype)) {
    // This combines ReadGenovecHphaseSubsetUnsafe() and PgrGet2()'s logic.
    const uint32_t aux1_first_byte = *fread_ptr++;
    const uint32_t aux1a_mode = aux1_first_byte & 15;
    const uint32_t aux1b_mode = aux1_first_byte >> 4;
    uint32_t raw_01_ct = 0;
    uint32_t raw_10_ct = 0;
    if ((!aux1a_mode) || (!aux1b_mode)) {
      GenovecCount12Unsafe(raw_genovec, raw_sample_ct, &raw_01_ct, &raw_10_ct);
    }
    uint32_t* deltalist_workspace = pgrp->workspace_difflist_sample_ids;
    if (!allele_idx0) {
      // Two cases:
      // - If allele_idx == 1, convert all aux1a entries from 01 to 11.
      // - Otherwise, for each matching aux1a entry, convert from 11 to 01.
      reterr = GenovecAux1aUpdate(fread_end, sample_include, sample_include_cumulative_popcounts, raw_genovec, aux1a_mode, raw_sample_ct, allele_ct, allele_idx1, 2, raw_01_ct, &fread_ptr, genovec, deltalist_workspace);
    } else {
      reterr = SkipAux1a(fread_end, aux1a_mode, raw_sample_ct, allele_ct, raw_01_ct, &fread_ptr);
    }
    if (unlikely(reterr)) {
      return reterr;
    }
    const unsigned char* aux1b_start = fread_ptr;
    reterr = GenovecAux1bUpdate2(fread_end, sample_include, sample_include_cumulative_popcounts, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, allele_idx0, allele_idx1, raw_10_ct, &fread_ptr, genovec, deltalist_workspace);
    if (unlikely(reterr)) {
      return reterr;
    }
    // Can have a modified version of GenovecAux1bUpdate2() which only requires
    // one pass, but let's keep the logic simpler for now since I don't expect
    // this function to be used frequently.
    uintptr_t* aux1b_hets = pgrp->workspace_aux1x_present;
    uint32_t aux1b_het_present;
    reterr = GetAux1bHets(fread_end, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &aux1b_start, aux1b_hets, &aux1b_het_present, deltalist_workspace);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (aux1b_het_present) {
      BitvecOr(aux1b_hets, BitCtToWordCt(raw_sample_ct), all_hets);
      if (!subsetting_required) {
        subsetted_10het = aux1b_hets;
      } else {
        // Don't need raw_genovec any more.
        CopyBitarrSubset(aux1b_hets, sample_include, sample_ct, raw_genovec);
        subsetted_10het = raw_genovec;
      }
    }
  }
  reterr = ParseAux2Subset(fread_end, sample_include, all_hets, subsetted_10het, raw_sample_ct, sample_ct, &fread_ptr, phasepresent, phaseinfo, phasepresent_ct_ptr, pgrp->workspace_subset);
  if (unlikely(reterr)) {
    return reterr;
  }
  if (VrtypeMultiallelicHc(vrtype) && (*phasepresent_ct_ptr)) {
    const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
    Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
    for (uint32_t hwidx = 0; hwidx != sample_ctl2; ++hwidx) {
      phasepresent_alias[hwidx] &= Pack01ToHalfword(genovec[hwidx]);
    }
    *phasepresent_ct_ptr = PopcountWords(phasepresent, BitCtToWordCt(sample_ct));
  }
  if (invert) {
    GenovecInvertUnsafe(sample_ct, genovec);
    if (*phasepresent_ct_ptr) {
      BitvecInvert(BitCtToWordCt(sample_ct), phaseinfo);
    }
  }
  return kPglRetSuccess;
}

PglErr PgrGetMP(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp) {
  pgvp->patch_01_ct = 0;
  pgvp->patch_10_ct = 0;
  if (!sample_ct) {
    pgvp->phasepresent_ct = 0;
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t multiallelic_hc_present = VrtypeMultiallelicHc(vrtype);
  if (!multiallelic_hc_present) {
    return ReadGenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, pgvp->genovec, pgvp->phasepresent, pgvp->phaseinfo, &(pgvp->phasepresent_ct));
  }
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  uintptr_t* all_hets = VrtypeHphase(vrtype)? pgrp->workspace_all_hets : nullptr;
  PglErr reterr = GetMultiallelicCodes(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, all_hets? (&fread_ptr) : nullptr, all_hets? (&fread_end) : nullptr, all_hets, pgvp);
  if (reterr || (!all_hets)) {
    return reterr;
  }
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  return ParseAux2Subset(fread_end, (sample_ct != raw_sample_ct)? sample_include : nullptr, all_hets, nullptr, raw_sample_ct, sample_ct, &fread_ptr, pgvp->phasepresent, pgvp->phaseinfo, &(pgvp->phasepresent_ct), pgrp->workspace_subset);
}

// ok for sample_include to be nullptr if not subsetting, though this is not
// required
PglErr ParseDosage16(const unsigned char* fread_ptr, const unsigned char* fread_end, const uintptr_t* __restrict sample_include, uint32_t sample_ct, uint32_t vidx, uint32_t allele_ct, PgenReader* pgrp, uint32_t* __restrict dosage_ct_ptr, uintptr_t* __restrict dphase_present, int16_t* dphase_delta, uint32_t* __restrict dphase_ct_ptr, uintptr_t* __restrict dosage_present, uint16_t* dosage_main) {
  // Side effect: may use pgrp->workspace_dosage_present and
  // pgrp->workspace_dphase_present
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* raw_dosage_present = subsetting_required? pgrp->workspace_dosage_present : dosage_present;
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t is_unconditional_dosage = ((vrtype & 0x60) == 0x40);
  uint32_t raw_dosage_ct;
  if ((vrtype & 0x60) == 0x20) {
    // case 1: dosage list
    if (unlikely(ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, &fread_ptr, raw_dosage_present, &raw_dosage_ct))) {
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
    if (dphase_ct_ptr) {
      *dphase_ct_ptr = 0;
    }
    return kPglRetSuccess;
  }
#ifdef __arm__
#  error "Unaligned accesses in ParseDosage16()."
#endif
  const uint16_t* dosage_main_read_iter = R_CAST(const uint16_t*, fread_ptr);
  uint16_t* dosage_main_write_iter = dosage_main;
  uint32_t raw_dphase_ct = 0;
  uint32_t dphase_ct = 0;
  uintptr_t* raw_dphase_present = nullptr;
  if (dphase_present && (vrtype & 0x80)) {
    fread_ptr = &(fread_ptr[raw_dosage_ct * 2]);
    if (!is_unconditional_dosage) {
      const uintptr_t* file_dphase_present = R_CAST(const uintptr_t*, fread_ptr);
      fread_ptr = &(fread_ptr[DivUp(raw_dosage_ct, CHAR_BIT)]);
      raw_dphase_present = subsetting_required? pgrp->workspace_dphase_present : dphase_present;
      ExpandBytearr(file_dphase_present, raw_dosage_present, raw_sample_ctl, raw_dosage_ct, 0, raw_dphase_present);
      raw_dphase_ct = PopcountWords(raw_dphase_present, raw_sample_ctl);
      dphase_ct = raw_dphase_ct;
      if (subsetting_required) {
        CopyBitarrSubset(raw_dphase_present, sample_include, sample_ct, dphase_present);
        dphase_ct = PopcountWords(dphase_present, sample_ctl);
      }
    } else {
      // raw_dphase_present = raw_dosage_present;
      dphase_ct = dosage_ct;
      SetAllBits(sample_ct, dphase_present);
    }
  }
  if (!dphase_ct) {
    if (allele_ct == 2) {
      if (!is_unconditional_dosage) {
        if (dosage_ct == raw_dosage_ct) {
          memcpy(dosage_main_write_iter, dosage_main_read_iter, dosage_ct * sizeof(int16_t));
        } else {
          // bugfix (22 May 2017): dosage_entry_idx needs to iterate up to
          // raw_dosage_ct, not dosage_ct
          uintptr_t widx = ~k0LU;
          uint32_t dosage_entry_idx = 0;
          do {
            uintptr_t cur_bits;
            do {
              cur_bits = raw_dosage_present[++widx];
            } while (!cur_bits);
            const uintptr_t sample_include_word = sample_include[widx];
            do {
              const uintptr_t low_bit = cur_bits & (-cur_bits);
              if (sample_include_word & low_bit) {
                *dosage_main_write_iter++ = dosage_main_read_iter[dosage_entry_idx];
              }
              ++dosage_entry_idx;
              cur_bits ^= low_bit;
            } while (cur_bits);
          } while (dosage_entry_idx != raw_dosage_ct);
        }
      } else {
        if (!subsetting_required) {
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            const uint16_t cur_dosage = *dosage_main_read_iter++;
            if (cur_dosage != 65535) {
              *dosage_main_write_iter++ = cur_dosage;
            } else {
              ClearBit(sample_idx, dosage_present);
            }
          }
        } else {
          uintptr_t widx = ~k0LU;
          uint32_t sample_idx = 0;
          do {
            uintptr_t cur_bits;
            do {
              cur_bits = sample_include[++widx];
            } while (!cur_bits);
            const uintptr_t sample_uidx_base = widx * kBitsPerWord;
            const uint16_t* dosage_main_readp = &(dosage_main_read_iter[sample_uidx_base]);
            do {
              const uint32_t sample_uidx_lowbits = ctzw(cur_bits);
              const uint16_t cur_dosage = dosage_main_readp[sample_uidx_lowbits];
              if (cur_dosage != 65535) {
                *dosage_main_write_iter++ = cur_dosage;
              } else {
                ClearBit(sample_idx, dosage_present);
              }
              ++sample_idx;
              cur_bits &= cur_bits - 1;
            } while (cur_bits);
          } while (sample_idx != sample_ct);
        }
        if (dosage_ct_ptr) {
          *dosage_ct_ptr = dosage_main_write_iter - dosage_main;
        }
      }
    } else {
      // todo: multiallelic dosage
      // need to support downcode to ref/nonref as well as raw load
      // (dosage_ct_ptr should be nullptr iff we're doing a raw load)
      fputs("multiallelic variants not yet supported by ParseDosage16()\n", stderr);
      exit(S_CAST(int32_t, kPglRetNotYetSupported));
      return kPglRetSuccess;
    }
    if (dphase_ct_ptr) {
      *dphase_ct_ptr = 0;
    }
  } else {
    // phased dosage
    if (allele_ct == 2) {
      if (!is_unconditional_dosage) {
        if (dphase_ct == raw_dphase_ct) {
          memcpy(dosage_main_write_iter, dosage_main_read_iter, dosage_ct * sizeof(int16_t));
          memcpy(dphase_delta, fread_ptr, dphase_ct * sizeof(int16_t));
          if (dphase_ct_ptr) {
            *dphase_ct_ptr = dphase_ct;
          }
        } else {
          uintptr_t widx = ~k0LU;
          uint32_t dosage_entry_idx = 0;
          do {
            uintptr_t cur_bits;
            do {
              cur_bits = raw_dosage_present[++widx];
            } while (!cur_bits);
            const uintptr_t sample_include_word = sample_include[widx];
            do {
              const uintptr_t low_bit = cur_bits & (-cur_bits);
              if (sample_include_word & low_bit) {
                *dosage_main_write_iter++ = dosage_main_read_iter[dosage_entry_idx];
              }
              ++dosage_entry_idx;
              cur_bits ^= low_bit;
            } while (cur_bits);
          } while (dosage_entry_idx != raw_dosage_ct);
          widx = ~k0LU;
          uint32_t dphase_entry_idx = 0;
          const int16_t* dphase_delta_read_alias = R_CAST(const int16_t*, fread_ptr);
          int16_t* dphase_delta_write_iter = dphase_delta;
          do {
            uintptr_t cur_bits;
            do {
              cur_bits = raw_dphase_present[++widx];
            } while (!cur_bits);
            const uintptr_t sample_include_word = sample_include[widx];
            do {
              const uintptr_t low_bit = cur_bits & (-cur_bits);
              if (sample_include_word & low_bit) {
                *dphase_delta_write_iter++ = dphase_delta_read_alias[dphase_entry_idx];
              }
              ++dphase_entry_idx;
              cur_bits ^= low_bit;
            } while (cur_bits);
          } while (dphase_entry_idx != raw_dphase_ct);
          if (dphase_ct_ptr) {
            *dphase_ct_ptr = dphase_delta_write_iter - dphase_delta;
          }
        }
      } else {
        const int16_t* dphase_delta_read = R_CAST(const int16_t*, fread_ptr);
        int16_t* dphase_delta_write_iter = dphase_delta;
        if (!subsetting_required) {
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            const uint16_t cur_dosage = *dosage_main_read_iter++;
            if (cur_dosage != 65535) {
              *dosage_main_write_iter++ = cur_dosage;
              const int16_t dphase_delta_val = dphase_delta_read[sample_idx];
              if (dphase_delta_val) {
                *dphase_delta_write_iter++ = dphase_delta_val;
              } else {
                ClearBit(sample_idx, dphase_present);
              }
            } else {
              // assert(dphase_delta_read[sample_idx] == -32768);
              ClearBit(sample_idx, dosage_present);
            }
          }
        } else {
          uintptr_t sample_uidx_base = 0;
          uintptr_t sample_include_bits = sample_include[0];
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
            const uint16_t cur_dosage = dosage_main_read_iter[sample_uidx];
            if (cur_dosage != 65535) {
              *dosage_main_write_iter++ = cur_dosage;
              const int16_t dphase_delta_val = dphase_delta_read[sample_uidx];
              if (dphase_delta_val) {
                *dphase_delta_write_iter++ = dphase_delta_val;
              } else {
                ClearBit(sample_idx, dphase_present);
              }
            } else {
              // assert(dphase_delta_read[sample_uidx] == -32768);
              ClearBit(sample_idx, dosage_present);
            }
          }
        }
        dosage_ct = dosage_main_write_iter - dosage_main;
        if (dosage_ct != sample_ct) {
          BitvecAnd(dosage_present, sample_ctl, dphase_present);
        }
        if (dosage_ct_ptr) {
          *dosage_ct_ptr = dosage_ct;
        }
        if (dphase_ct_ptr) {
          *dphase_ct_ptr = dphase_delta_write_iter - dphase_delta;
        }
      }
    } else {
      // multiallelic subcase
      fputs("multiallelic variants not yet supported by ParseDosage16()\n", stderr);
      exit(S_CAST(int32_t, kPglRetNotYetSupported));
      return kPglRetSuccess;
    }
  }
  return kPglRetSuccess;
}

PglErr PgrGetD(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    *dosage_ct_ptr = 0;
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  if ((!VrtypeDosage(vrtype)) || (!dosage_present)) {
    *dosage_ct_ptr = 0;
    return ReadGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, genovec);
  }
  const unsigned char* fread_ptr = nullptr;
  const unsigned char* fread_end = nullptr;
  uint32_t phasepresent_ct;
  PglErr reterr = ReadGenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, &fread_ptr, &fread_end, genovec, nullptr, nullptr, &phasepresent_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  return ParseDosage16(fread_ptr, fread_end, sample_include, sample_ct, vidx, allele_ct, pgrp, dosage_ct_ptr, nullptr, nullptr, nullptr, dosage_present, dosage_main);
}

PglErr PgrGet1D(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_countvec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr) {
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  if ((allele_ct == 2) || (!allele_idx)) {
    uint32_t dosage_ct;
    PglErr reterr = PgrGetD(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, allele_countvec, dosage_present, dosage_main, &dosage_ct);
    if (!allele_idx) {
      GenovecInvertUnsafe(sample_ct, allele_countvec);
      if (dosage_ct) {
        BiallelicDosage16Invert(dosage_ct, dosage_main);
      }
    }
    *dosage_ct_ptr = dosage_ct;
    return reterr;
  }
  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  if (!VrtypeDosage(vrtype)) {
    *dosage_ct_ptr = 0;
    return PgrGet1(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, allele_countvec);
  }
  fputs("multiallelic variants not yet supported by PgrGet1D()\n", stderr);
  exit(S_CAST(int32_t, kPglRetNotYetSupported));
  return kPglRetSuccess;
}

PglErr PgrGetInv1D(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_invcountvec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr) {
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  if ((allele_ct == 2) || (!allele_idx)) {
    uint32_t dosage_ct;
    PglErr reterr = PgrGetD(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, allele_invcountvec, dosage_present, dosage_main, &dosage_ct);
    if (allele_idx) {
      GenovecInvertUnsafe(sample_ct, allele_invcountvec);
      if (dosage_ct) {
        BiallelicDosage16Invert(dosage_ct, dosage_main);
      }
    }
    *dosage_ct_ptr = dosage_ct;
    return reterr;
  }
  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  if (!VrtypeDosage(vrtype)) {
    *dosage_ct_ptr = 0;
    return PgrGetInv1(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, allele_invcountvec);
  }
  fputs("multiallelic variants not yet supported by PgrGetInv1D()\n", stderr);
  exit(S_CAST(int32_t, kPglRetNotYetSupported));
  return kPglRetSuccess;
}

PglErr GetAux1bHetIncr(const unsigned char* fread_end, uint32_t aux1b_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t raw_10_ct, const unsigned char** fread_pp, uint32_t* __restrict raw_het_ctp) {
  if (aux1b_mode == 15) {
    return kPglRetSuccess;
  }
  uint32_t rare10_ct;
  if (!aux1b_mode) {
    const uint32_t fset_byte_ct = DivUp(raw_10_ct, 8);
    rare10_ct = PopcountBytes(*fread_pp, fset_byte_ct);
    *fread_pp += fset_byte_ct;
  } else {
    // aux1b_mode == 1
    const unsigned char* group_info_iter;
    PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, &rare10_ct);
    if (unlikely(reterr)) {
      return reterr;
    }
    reterr = SkipDeltalistIds(fread_end, group_info_iter, rare10_ct, raw_sample_ct, 0, fread_pp);
    if (unlikely(reterr)) {
      return reterr;
    }
  }
  uintptr_t detect_hom_mask_lo;
  const uint32_t allele_code_logwidth = GetAux1bConsts(allele_ct, &detect_hom_mask_lo);
  const uint32_t code10_logwidth = allele_code_logwidth + (allele_code_logwidth != 0);
#ifdef __arm__
#  error "Unaligned accesses in GetAux1bHetIncr()."
#endif
  const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) << code10_logwidth, CHAR_BIT);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  if (allele_ct == 3) {
    const uint32_t hom22_ct = PopcountBytes(patch_10_fvalsw, fvals_byte_ct);
    *raw_het_ctp += rare10_ct - hom22_ct;
    return kPglRetSuccess;
  }
  // possible todo: vectorized het-counter, analogous to CountAux1bDense()
  const uint32_t code10_width = 1U << code10_logwidth;
  const uint32_t allele_code_width = 1U << allele_code_logwidth;
  const uintptr_t detect_all_mask_lo = detect_hom_mask_lo | (detect_hom_mask_lo << allele_code_width);
  const uintptr_t detect_all_mask_hi = detect_all_mask_lo << (allele_code_width - 1);
  const uintptr_t detect_hom_mask_hi = detect_hom_mask_lo << (code10_width - 1);
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  uint32_t het_incr = 0;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        break;
      }
      fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
    } else {
      fvals_bits = patch_10_fvalsw[fvals_widx];
    }
    // allele_ct > 3 guaranteed
    fvals_bits = fvals_bits ^ (fvals_bits << allele_code_width);
    fvals_bits = detect_hom_mask_hi & (fvals_bits | ((fvals_bits | detect_all_mask_hi) - detect_all_mask_lo));
    if (fvals_widx == fvals_word_ct_m1) {
      fvals_bits = bzhi_max(fvals_bits, ModNz(rare10_ct << code10_logwidth, kBitsPerWord));
    }
    het_incr += PopcountWord(fvals_bits);
  }
  *raw_het_ctp += het_incr;
  return kPglRetSuccess;
}

uint64_t U16VecSum(const uint16_t* __restrict uint16_vec, uint32_t entry_ct) {
#ifdef __LP64__
  // UniVecHsum32() could overflow once we exceed this
  const uint32_t max_loop_len = (131072 / kInt32PerVec) - 1;

  const VecW m16 = VCONST_W(kMask0000FFFF);
  const VecW* uint16_vvec_iter = R_CAST(const VecW*, uint16_vec);
  uint64_t sum = 0;
  for (uint32_t full_vecs_remaining = entry_ct / (kBytesPerVec / sizeof(int16_t)); ; ) {
    UniVec acc_even;
    UniVec acc_odd;
    acc_even.vw = vecw_setzero();
    acc_odd.vw = vecw_setzero();
    const VecW* uint16_vvec_stop;
    if (full_vecs_remaining < max_loop_len) {
      if (!full_vecs_remaining) {
        const uint32_t trail_ct = entry_ct % (kBytesPerVec / sizeof(int16_t));
        uint16_vec = R_CAST(const uint16_t*, uint16_vvec_iter);
        for (uint32_t uii = 0; uii != trail_ct; ++uii) {
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
      const VecW cur_vec = *uint16_vvec_iter++;
      acc_even.vw = acc_even.vw + (cur_vec & m16);
      acc_odd.vw = acc_odd.vw + (vecw_srli(cur_vec, 16) & m16);
    } while (uint16_vvec_iter < uint16_vvec_stop);
    sum += UniVecHsum32(acc_even);
    sum += UniVecHsum32(acc_odd);
  }
#else
  uint64_t sum = 0;
  for (uint32_t uii = 0; uii != entry_ct; ++uii) {
    sum += uint16_vec[uii];
  }
  return sum;
#endif
}

PglErr GetPhasepresentAndSkipPhaseinfo(const unsigned char* fread_end, const uintptr_t* __restrict all_hets, uint32_t raw_sample_ct, uint32_t het_ct, const unsigned char** fread_pp, uintptr_t* __restrict phasepresent, uint32_t* __restrict phasepresent_ctp) {
  const unsigned char* aux2_start = *fread_pp;
  const uint32_t aux2_first_part_byte_ct = 1 + (het_ct / CHAR_BIT);
  if (PtrAddCk(fread_end, aux2_first_part_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  if (!(aux2_start[0] & 1)) {
    memcpy(phasepresent, all_hets, raw_sample_ctl * kBytesPerWord);
    *phasepresent_ctp = het_ct;
    return kPglRetSuccess;
  }
  const uint32_t phasepresent_ct = PopcountBytes(aux2_start, aux2_first_part_byte_ct) - 1;
  if (PtrAddCk(fread_end, DivUp(phasepresent_ct, CHAR_BIT), fread_pp)) {
    return kPglRetMalformedInput;
  }
  *phasepresent_ctp = phasepresent_ct;
  ExpandBytearr(aux2_start, all_hets, raw_sample_ctl, het_ct, 1, phasepresent);
  return kPglRetSuccess;
}

PglErr GetUnphasedBiallelicHetCt(const uintptr_t* __restrict sample_include, const uintptr_t* raw_genovec, const unsigned char* fread_ptr, const unsigned char* fread_end, uint32_t subsetted_het_ct, PgenReader* pgrp, uint32_t* unphased_het_ctp) {
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  uint32_t raw_het_ct;
  if (!sample_include) {
    raw_het_ct = subsetted_het_ct;
  } else {
    raw_het_ct = CountQuater(raw_genovec, kMask5555, raw_sample_ct);
  }
  const uint32_t aux2_first_part_byte_ct = 1 + (raw_het_ct / CHAR_BIT);
  if (PtrCheck(fread_end, fread_ptr, aux2_first_part_byte_ct)) {
    return kPglRetMalformedInput;
  }
  const uint32_t explicit_phasepresent = fread_ptr[0] & 1;
  if (!explicit_phasepresent) {
    // initial value of 0 is correct
    return kPglRetSuccess;
  }
  if (raw_het_ct == subsetted_het_ct) {
    *unphased_het_ctp = raw_het_ct + 1 - PopcountBytes(fread_ptr, aux2_first_part_byte_ct);
    return kPglRetSuccess;
  }
  // A dedicated counting function would be faster, but this case
  // should rarely come up.
  uintptr_t* all_hets = pgrp->workspace_all_hets;
  PgrDetectGenovecHets(raw_genovec, raw_sample_ct, all_hets);
  uintptr_t* raw_phasepresent = pgrp->workspace_subset;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  // todo: compare against ExpandThenSubsetBytearr followed by simple popcount
  ExpandBytearr(fread_ptr, all_hets, raw_sample_ctl, raw_het_ct, 1, raw_phasepresent);
  *unphased_het_ctp = subsetted_het_ct - PopcountWordsIntersect(raw_phasepresent, sample_include, raw_sample_ctl);
  return kPglRetSuccess;
}

PglErr GetPhasedBiallelicGenotypeSubsetCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uint32_t* unphased_het_ctp, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // Currently much less optimized than the other count functions.  (This case
  // shouldn't come up much, the user has to be computing minimac-r2 on a file
  // with no dosages...)
  uintptr_t* raw_genovec = pgrp->workspace_vec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRawGenovec(1, vidx, pgrp, &fread_ptr, &fread_end, raw_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  ZeroTrailingQuaters(raw_sample_ct, raw_genovec);
  GenovecCountSubsetFreqs(raw_genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
  return GetUnphasedBiallelicHetCt(sample_include, raw_genovec, fread_ptr, fread_end, genocounts[1], pgrp, unphased_het_ctp);
}

// Imputation r^2 computation:
// * This function assumes the biallelic diploid case.  Divide by two to get
//   the biallelic haploid value, for whatever that's worth.
// * chrX requires sex information, so that's handled directly in
//   LoadAlleleAndGenoCountsThread()... er, actually, we just give up on that
//   for now.
// * See PgrGetMDCounts() support functions below for multiallelic-diploid
//   notes.
PglErr GetBasicGenotypeCountsAndDosage16s(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t is_minimac3_r2, PgenReader* pgrp, double* imp_r2_ptr, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* __restrict all_dosages) {
  // genocounts[0] := ref/ref, genocounts[1] := ref/altx,
  // genocounts[2] := altx/alty, genocounts[3] := missing
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uint32_t unphased_het_ct = 0;
  // To avoid LD cache thrashing, we try to either always keep a subsetted
  // cache, or never do so.  (Always, when only hardcalls are present;
  // otherwise never.)
  if ((!(pgrp->fi.gflags & kfPgenGlobalDosagePresent)) ||
      ((!(vrtype & 0x60)) && (!subsetting_required))) {
    {
      const uint32_t need_unphased_het_ct = is_minimac3_r2 && VrtypeHphase(vrtype);
      PglErr reterr;
      if (!(subsetting_required && need_unphased_het_ct)) {
        reterr = GetBasicGenotypeCounts(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, need_unphased_het_ct? (&unphased_het_ct) : nullptr, genocounts);
      } else {
        reterr = GetPhasedBiallelicGenotypeSubsetCounts(sample_include, sample_include_interleaved_vec, sample_ct, vidx, pgrp, &unphased_het_ct, genocounts);
      }
      if (unlikely(reterr)) {
        return reterr;
      }
    }
  GetBasicGenotypeCountsAndDosage16s_basic_finish:
    all_dosages[0] = (genocounts[0] * 2 + genocounts[1]) * 16384LLU;
    all_dosages[1] = (genocounts[2] * 2 + genocounts[1]) * 16384LLU;
    if (!imp_r2_ptr) {
      return kPglRetSuccess;
    }
    // yeah, it's sinful to implement imputation r2 here...
    const uint32_t nm_sample_ct = sample_ct - genocounts[3];
    const uint64_t alt1_dosage = genocounts[2] * 0x8000LLU + genocounts[1] * 0x4000LLU;
    uint64_t hap_alt1_ssq_x2 = genocounts[2] * 0x40000000LLU + genocounts[1] * 0x10000000LLU;
    if (is_minimac3_r2) {
      if (!VrtypeHphase(vrtype)) {
        unphased_het_ct = genocounts[1];
      }
      hap_alt1_ssq_x2 += (genocounts[1] - unphased_het_ct) * 0x10000000LLU;
    }
    *imp_r2_ptr = BiallelicDiploidMinimac3R2(alt1_dosage, hap_alt1_ssq_x2, nm_sample_ct);
    if (!is_minimac3_r2) {
      *imp_r2_ptr *= 2;
    }
    return kPglRetSuccess;
  }
  uintptr_t* raw_genovec = pgrp->workspace_vec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRawGenovec(subsetting_required, vidx, pgrp, &fread_ptr, &fread_end, raw_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  ZeroTrailingQuaters(raw_sample_ct, raw_genovec);
  if (!subsetting_required) {
    GenovecCountFreqsUnsafe(raw_genovec, raw_sample_ct, genocounts);
  } else {
    GenovecCountSubsetFreqs(raw_genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
  }
  if (!(vrtype & 0x60)) {
    if (is_minimac3_r2 && VrtypeHphase(vrtype)) {
      assert(!VrtypeMultiallelicHc(vrtype));
      reterr = GetUnphasedBiallelicHetCt(subsetting_required? sample_include : nullptr, raw_genovec, fread_ptr, fread_end, genocounts[1], pgrp, &unphased_het_ct);
      if (unlikely(reterr)) {
        return reterr;
      }
    }
    goto GetBasicGenotypeCountsAndDosage16s_basic_finish;
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  uintptr_t* raw_phasepresent = pgrp->workspace_subset;
  uint32_t raw_phasepresent_ct = 0;
  if (VrtypeHphase(vrtype)) {
    uint32_t raw_het_ct = genocounts[1];  // inaccurate if subsetting_required
    if (!is_minimac3_r2) {
      if (VrtypeMultiallelicHc(vrtype)) {
        const uint32_t aux1_first_byte = *fread_ptr++;
        const uint32_t aux1a_mode = aux1_first_byte & 15;
        const uint32_t aux1b_mode = aux1_first_byte >> 4;
        uint32_t raw_10_ct = 0;
        if ((!aux1a_mode) || (!aux1b_mode) || subsetting_required) {
          GenovecCount12Unsafe(raw_genovec, raw_sample_ct, &raw_het_ct, &raw_10_ct);
        }
        reterr = SkipAux1a(fread_end, aux1a_mode, raw_sample_ct, allele_ct, raw_het_ct, &fread_ptr);
        if (unlikely(reterr)) {
          return reterr;
        }
        reterr = GetAux1bHetIncr(fread_end, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &fread_ptr, &raw_het_ct);
        if (unlikely(reterr)) {
          return reterr;
        }
      } else if (subsetting_required) {
        raw_het_ct = CountQuater(raw_genovec, kMask5555, raw_sample_ct);
      }
      reterr = SkipAux2(fread_end, raw_het_ct, &fread_ptr, nullptr);
      if (unlikely(reterr)) {
        return reterr;
      }
    } else {
      assert(!VrtypeMultiallelicHc(vrtype));
      uintptr_t* all_hets = pgrp->workspace_all_hets;
      PgrDetectGenovecHets(raw_genovec, raw_sample_ct, all_hets);
      if (subsetting_required) {
        raw_het_ct = PopcountWords(all_hets, raw_sample_ctl);
      }
      const uint32_t first_half_byte_ct = 1 + (raw_het_ct / CHAR_BIT);
      const uint32_t explicit_phasepresent = fread_ptr[0] & 1;
      if (explicit_phasepresent) {
        ExpandBytearr(fread_ptr, all_hets, raw_sample_ctl, raw_het_ct, 1, raw_phasepresent);
        raw_phasepresent_ct = PopcountBytes(fread_ptr, first_half_byte_ct) - 1;
        const uint32_t second_half_byte_ct = DivUp(raw_phasepresent_ct, CHAR_BIT);
        fread_ptr = &(fread_ptr[first_half_byte_ct + second_half_byte_ct]);
      } else {
        raw_phasepresent_ct = raw_het_ct;
        memcpy(raw_phasepresent, all_hets, raw_sample_ctl * sizeof(intptr_t));
        fread_ptr = &(fread_ptr[first_half_byte_ct]);
      }
    }
  } else if (VrtypeMultiallelicHc(vrtype)) {
    reterr = SkipAux1(fread_end, raw_genovec, raw_sample_ct, allele_ct, &fread_ptr);
    if (unlikely(reterr)) {
      return reterr;
    }
  }
  if (allele_ct != 2) {
    // Maybe make this an invalid function call?  If that happens, the
    // VrtypeMultiallelicHc() branch above can be removed.
    fputs("multiallelic dosages not yet supported by GetBasicGenotypeCountsAndDosage16s()\n", stderr);
    exit(S_CAST(int32_t, kPglRetNotYetSupported));
    return kPglRetSuccess;
  }

  const uint32_t is_unconditional_dosage = ((vrtype & 0x60) == 0x40);
  uint64_t alt1_dosage = 0;
  uint32_t dosage_ct = 0;
  STD_ARRAY_DECL(uint32_t, 4, replaced_genocounts);
  if ((!is_minimac3_r2) || (!(vrtype & 0x90))) {
    uint64_t alt1_dosage_sq_sum = 0;
    if (is_unconditional_dosage) {
      // needs to be handled separately from the other cases due to possible
      // presence of missing values.
      // note that this code will also need to be adjusted when multiallelic
      // support is added.
#ifdef __arm__
#  error "Unaligned accesses in GetBasicGenotypeCountsAndDosage16s()."
#endif
      STD_ARRAY_FILL0(replaced_genocounts);
      const uint16_t* dosage_main = R_CAST(const uint16_t*, fread_ptr);
      if (PtrAddCk(fread_end, raw_sample_ct * sizeof(int16_t), &fread_ptr)) {
        return kPglRetMalformedInput;
      }
      if (subsetting_required) {
        uintptr_t sample_uidx_base = 0;
        uintptr_t sample_include_bits = sample_include[0];
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
          const uintptr_t cur_dosage_val = dosage_main[sample_uidx];
          if (cur_dosage_val != 65535) {
            alt1_dosage += cur_dosage_val;

            // todo: check if this is slow enough to justify removing it from
            // the main loop
            alt1_dosage_sq_sum += cur_dosage_val * cur_dosage_val;
            ++dosage_ct;
          }
        }
      } else {
        for (uint32_t sample_uidx = 0; sample_uidx != sample_ct; ++sample_uidx) {
          const uintptr_t cur_dosage_val = dosage_main[sample_uidx];
          if (cur_dosage_val != 65535) {
            alt1_dosage += cur_dosage_val;
            alt1_dosage_sq_sum += cur_dosage_val * cur_dosage_val;
            ++dosage_ct;
          }
        }
      }
      // update (20 Mar 2019): .pgen specification tightened to remove the need
      // to update replaced_genocounts in the main loops above.
      STD_ARRAY_COPY(genocounts, 4, replaced_genocounts);
      replaced_genocounts[3] = replaced_genocounts[3] + dosage_ct - sample_ct;
    } else {
      uintptr_t* raw_dosage_present = pgrp->workspace_dosage_present;
      uint32_t raw_dosage_ct;
      if (!(vrtype & 0x40)) {
        // dosage list
        if (unlikely(ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, &fread_ptr, raw_dosage_present, &raw_dosage_ct))) {
          return kPglRetMalformedInput;
        }
      } else {
        // dosage bitarray
        raw_dosage_present[raw_sample_ctl - 1] = 0;
        const uint32_t raw_sample_ctb = DivUp(raw_sample_ct, CHAR_BIT);
        memcpy(raw_dosage_present, fread_ptr, raw_sample_ctb);
        fread_ptr = &(fread_ptr[raw_sample_ctb]);
        raw_dosage_ct = PopcountWords(raw_dosage_present, raw_sample_ctl);
      }
      const uint16_t* dosage_main_iter = R_CAST(const uint16_t*, fread_ptr);
      if (PtrAddCk(fread_end, raw_dosage_ct * sizeof(int16_t), &fread_ptr)) {
        return kPglRetMalformedInput;
      }
      if (subsetting_required) {
        uintptr_t sample_widx = 0;
        uintptr_t dosage_present_bits = raw_dosage_present[0];
        for (uint32_t dosage_idx = 0; dosage_idx != raw_dosage_ct; ++dosage_idx) {
          const uintptr_t lowbit = BitIter1y(raw_dosage_present, &sample_widx, &dosage_present_bits);
          if (sample_include[sample_widx] & lowbit) {
            const uintptr_t cur_dosage_val = dosage_main_iter[dosage_idx];
            alt1_dosage += cur_dosage_val;
            alt1_dosage_sq_sum += cur_dosage_val * cur_dosage_val;
            ++dosage_ct;
          }
        }
        GenoarrCountSubsetIntersectFreqs(raw_genovec, raw_dosage_present, sample_include, raw_sample_ct, replaced_genocounts);
      } else {
        if (!imp_r2_ptr) {
          for (uint32_t dosage_idx = 0; dosage_idx != raw_dosage_ct; ++dosage_idx) {
            alt1_dosage += dosage_main_iter[dosage_idx];
          }
        } else {
          for (uint32_t dosage_idx = 0; dosage_idx != raw_dosage_ct; ++dosage_idx) {
            const uintptr_t cur_dosage_val = dosage_main_iter[dosage_idx];
            alt1_dosage += cur_dosage_val;
            alt1_dosage_sq_sum += cur_dosage_val * cur_dosage_val;
          }
        }
        dosage_ct = raw_dosage_ct;
        GenovecCountSubsetFreqs2(raw_genovec, raw_dosage_present, raw_sample_ct, raw_dosage_ct, replaced_genocounts);
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
    if (!imp_r2_ptr) {
      return kPglRetSuccess;
    }
    // possible todo: also move all-hardcall-phase-present, no-dosage
    // is_minimac3_r2 case under this branch, since we can just set imp_r2 to
    // NaN or 1.
    // 16384^2, 32768^2
    alt1_dosage_sq_sum += remaining_het_ct * 0x10000000LLU + remaining_hom_alt_ct * 0x40000000LLU;
    *imp_r2_ptr = BiallelicDiploidMinimac3R2(alt1_dosage, alt1_dosage_sq_sum, new_sample_nm_ct);
    if (!is_minimac3_r2) {
      *imp_r2_ptr *= 2;
    }
    return kPglRetSuccess;
  }
  // Need to deal with implicitly phased dosages.  Best to have raw_genovec,
  // raw_phasepresent, dosage_present, and dosage_main all available, then loop
  // over everything at once.
  // (phaseinfo is irrelevant since only absolute value of (left - right)
  // matters.)

  // We have the following 2x2x3 cases to deal with:
  // - Subsetted vs. un-subsetted.  Un-subsetted comes up a lot, so we have an
  //   optimized code path for it.
  // - Unconditional vs. conditional dosage.  Unconditional should not come up
  //   much, so we just mock up raw_dosage_present... er, actually, that
  //   doesn't work because dosage_main would also need to be collapsed.  Sigh.
  //   Ok, it's still handled separately.
  // - Only hardcall-phase, vs. only dosage-phase, vs. both.  At least we can
  //   merge the "only dosage-phase" and "both" cases.
  // So we end up with 8 primary code paths.
  // This is kind of a nightmare; it would obviously be nicer to move this
  // out of pgenlib_internal, and that may eventually happen.  But we don't
  // want users to be discouraged from running --minimac3-r2-filter when it's
  // appropriate just because it's a lot slower than other standard filters;
  // and this also serves as a testing ground for efficient phased-dosage
  // handling strategies.
  if (!VrtypeHphase(vrtype)) {
    ZeroWArr(raw_sample_ctl, raw_phasepresent);
  }
  uintptr_t* raw_dosage_present = nullptr;
  const uint16_t* dosage_main;
  uint32_t raw_dosage_ct = 0;
  if (is_unconditional_dosage) {
    dosage_main = R_CAST(const uint16_t*, fread_ptr);
    if (PtrAddCk(fread_end, raw_sample_ct * sizeof(int16_t), &fread_ptr)) {
      return kPglRetMalformedInput;
    }
    // raw_dosage_ct unused in this case.
  } else {
    // could move some duplicate code before the big branch
    raw_dosage_present = pgrp->workspace_dosage_present;
    if (!(vrtype & 0x40)) {
      // dosage list
      if (unlikely(ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, &fread_ptr, raw_dosage_present, &raw_dosage_ct))) {
        return kPglRetMalformedInput;
      }
    } else {
      // dosage bitarray
      raw_dosage_present[raw_sample_ctl - 1] = 0;
      const uint32_t raw_sample_ctb = DivUp(raw_sample_ct, CHAR_BIT);
      memcpy(raw_dosage_present, fread_ptr, raw_sample_ctb);
      fread_ptr = &(fread_ptr[raw_sample_ctb]);
      raw_dosage_ct = PopcountWords(raw_dosage_present, raw_sample_ctl);
    }
    dosage_main = R_CAST(const uint16_t*, fread_ptr);
    if (PtrAddCk(fread_end, raw_dosage_ct * sizeof(int16_t), &fread_ptr)) {
      return kPglRetMalformedInput;
    }
  }
  const uint16_t* dosage_main_iter = dosage_main;
  uint64_t hap_ssq_x2 = 0;
  uint32_t phased_hc_het_ct = 0;
  if (!(vrtype & 0x80)) {
    if (is_unconditional_dosage) {
      if (!subsetting_required) {
        const uint32_t raw_sample_ctl_m1 = raw_sample_ctl - 1;
        uint32_t loop_len = kBitsPerWord;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= raw_sample_ctl_m1) {
            if (widx > raw_sample_ctl_m1) {
              break;
            }
            loop_len = ModNz(raw_sample_ct, kBitsPerWord);
          }
          uintptr_t phasepresent_word = raw_phasepresent[widx];
          for (uint32_t uii = 0; uii != loop_len; ++uii) {
            const uintptr_t cur_dosage_val = *dosage_main_iter++;
            if (cur_dosage_val != 65535) {
              alt1_dosage += cur_dosage_val;
              hap_ssq_x2 += cur_dosage_val * cur_dosage_val;
              ++dosage_ct;
              if (phasepresent_word & 1) {
                // For each dosage, when phasepresent bit is set, implicit
                // dphase_delta value is 16384 - |16384 - x|.
                const uintptr_t homdist = 16384 - abs_i32(16384 - cur_dosage_val);
                hap_ssq_x2 += homdist * homdist;
              }
            }
            phasepresent_word = phasepresent_word >> 1;
          }
        }
      } else {
        for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
          uintptr_t sample_include_word = sample_include[widx];
          if (!sample_include_word) {
            continue;
          }
          const uintptr_t phasepresent_word = raw_phasepresent[widx];
          const uint16_t* cur_dosage_main = &(dosage_main[widx * kBitsPerWord]);
          do {
            const uint32_t sample_idx_lowbits = ctzw(sample_include_word);
            const uintptr_t cur_dosage_val = cur_dosage_main[sample_idx_lowbits];
            const uintptr_t lowbit = sample_include_word & (-sample_include_word);
            if (cur_dosage_val != 65535) {
              alt1_dosage += cur_dosage_val;
              hap_ssq_x2 += cur_dosage_val * cur_dosage_val;
              ++dosage_ct;
              if (lowbit & phasepresent_word) {
                const uintptr_t homdist = 16384 - abs_i32(16384 - cur_dosage_val);
                hap_ssq_x2 += homdist * homdist;
              }
            }
            sample_include_word ^= lowbit;
          } while (sample_include_word);
        }
      }
      STD_ARRAY_COPY(genocounts, 4, replaced_genocounts);
      replaced_genocounts[3] = replaced_genocounts[3] + dosage_ct - sample_ct;
    } else {  // !is_unconditional_dosage
      if (!subsetting_required) {
        // phased_hc_het_ct := popcount(phasepresent & (~dosage_present))
        phased_hc_het_ct = raw_phasepresent_ct - PopcountWordsIntersect(raw_phasepresent, raw_dosage_present, raw_sample_ctl);

        for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
          uintptr_t dosage_present_word = raw_dosage_present[widx];
          if (dosage_present_word) {
            const uintptr_t phasepresent_word = raw_phasepresent[widx];
            do {
              const uintptr_t cur_dosage_val = *dosage_main_iter++;
              alt1_dosage += cur_dosage_val;
              const uintptr_t lowbit = dosage_present_word & (-dosage_present_word);
              hap_ssq_x2 += cur_dosage_val * cur_dosage_val;
              if (lowbit & phasepresent_word) {
                const uintptr_t homdist = 16384 - abs_i32(16384 - cur_dosage_val);
                hap_ssq_x2 += homdist * homdist;
              }
              dosage_present_word ^= lowbit;
            } while (dosage_present_word);
          }
        }
        dosage_ct = raw_dosage_ct;
        GenovecCountSubsetFreqs2(raw_genovec, raw_dosage_present, raw_sample_ct, raw_dosage_ct, replaced_genocounts);
      } else {
        for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
          const uintptr_t sample_include_word = sample_include[widx];
          uintptr_t dosage_present_word = raw_dosage_present[widx];
          if (!sample_include_word) {
            dosage_main_iter = &(dosage_main_iter[PopcountWord(dosage_present_word)]);
            continue;
          }
          const uintptr_t phasepresent_word = raw_phasepresent[widx];
          phased_hc_het_ct += PopcountWord(sample_include_word & phasepresent_word & (~dosage_present_word));
          while (dosage_present_word) {
            const uintptr_t lowbit = dosage_present_word & (-dosage_present_word);
            if (lowbit & sample_include_word) {
              const uintptr_t cur_dosage_val = *dosage_main_iter;
              alt1_dosage += cur_dosage_val;
              hap_ssq_x2 += cur_dosage_val * cur_dosage_val;
              ++dosage_ct;
              if (lowbit & phasepresent_word) {
                const uintptr_t homdist = 16384 - abs_i32(16384 - cur_dosage_val);
                hap_ssq_x2 += homdist * homdist;
              }
            }
            dosage_present_word ^= lowbit;
            ++dosage_main_iter;
          }
        }
        GenoarrCountSubsetIntersectFreqs(raw_genovec, raw_dosage_present, sample_include, raw_sample_ct, replaced_genocounts);
      }
    }
  } else {
    if (is_unconditional_dosage) {
      if (PtrCheck(fread_end, fread_ptr, raw_sample_ct * sizeof(int16_t))) {
        return kPglRetMalformedInput;
      }
      const int16_t* dphase_delta = R_CAST(const int16_t*, fread_ptr);
      if (!subsetting_required) {
        for (uint32_t sample_uidx = 0; sample_uidx != raw_sample_ct; ++sample_uidx) {
          const uintptr_t cur_dosage_val = dosage_main[sample_uidx];
          if (cur_dosage_val != 65535) {
            alt1_dosage += cur_dosage_val;
            hap_ssq_x2 += cur_dosage_val * cur_dosage_val;
            ++dosage_ct;
            // .pgen specification now requires this value to never be missing.
            const intptr_t dphase_delta_val = dphase_delta[sample_uidx];
            hap_ssq_x2 += dphase_delta_val * dphase_delta_val;
          }
        }
      } else {
        for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
          uintptr_t sample_include_word = sample_include[widx];
          if (!sample_include_word) {
            continue;
          }
          const uint16_t* cur_dosage_main = &(dosage_main[widx * kBitsPerWord]);
          const int16_t* cur_dphase_delta = &(dphase_delta[widx * kBitsPerWord]);
          do {
            const uint32_t sample_idx_lowbits = ctzw(sample_include_word);
            const uintptr_t cur_dosage_val = cur_dosage_main[sample_idx_lowbits];
            if (cur_dosage_val != 65535) {
              alt1_dosage += cur_dosage_val;
              hap_ssq_x2 += cur_dosage_val * cur_dosage_val;
              ++dosage_ct;
              const intptr_t dphase_delta_val = cur_dphase_delta[sample_idx_lowbits];
              hap_ssq_x2 += dphase_delta_val * dphase_delta_val;
            }
            sample_include_word &= sample_include_word - 1;
          } while (sample_include_word);
        }
      }
      STD_ARRAY_COPY(genocounts, 4, replaced_genocounts);
      replaced_genocounts[3] = replaced_genocounts[3] + dosage_ct - sample_ct;
    } else {
      const uintptr_t* file_dphase_present = R_CAST(const uintptr_t*, fread_ptr);
      const uint32_t raw_dosage_ctb = DivUp(raw_dosage_ct, CHAR_BIT);
      if (PtrAddCk(fread_end, raw_dosage_ctb, &fread_ptr)) {
        return kPglRetMalformedInput;
      }
      const uint32_t raw_dphase_ct = PopcountBytes(file_dphase_present, raw_dosage_ctb);
      if (PtrCheck(fread_end, fread_ptr, raw_dphase_ct * sizeof(int16_t))) {
        return kPglRetMalformedInput;
      }
      uintptr_t* raw_dphase_present = pgrp->workspace_dphase_present;
      ExpandBytearr(file_dphase_present, raw_dosage_present, raw_sample_ctl, raw_dosage_ct, 0, raw_dphase_present);
      const int16_t* dphase_delta_iter = R_CAST(const int16_t*, fread_ptr);
      if (!subsetting_required) {
        phased_hc_het_ct = raw_phasepresent_ct - PopcountWordsIntersect(raw_phasepresent, raw_dosage_present, raw_sample_ctl);

        for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
          uintptr_t dosage_present_word = raw_dosage_present[widx];
          if (dosage_present_word) {
            const uintptr_t phasepresent_word = raw_phasepresent[widx];
            const uintptr_t dphase_present_word = raw_dphase_present[widx];
            do {
              const uintptr_t cur_dosage_val = *dosage_main_iter++;
              alt1_dosage += cur_dosage_val;
              const uintptr_t lowbit = dosage_present_word & (-dosage_present_word);
              hap_ssq_x2 += cur_dosage_val * cur_dosage_val;
              if (lowbit & dphase_present_word) {
                const intptr_t dphase_delta_val = *dphase_delta_iter++;
                hap_ssq_x2 += dphase_delta_val * dphase_delta_val;
              } else if (lowbit & phasepresent_word) {
                const uintptr_t homdist = 16384 - abs_i32(16384 - cur_dosage_val);
                hap_ssq_x2 += homdist * homdist;
              }
              dosage_present_word ^= lowbit;
            } while (dosage_present_word);
          }
        }
        dosage_ct = raw_dosage_ct;
        GenovecCountSubsetFreqs2(raw_genovec, raw_dosage_present, raw_sample_ct, raw_dosage_ct, replaced_genocounts);
      } else {
        for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
          const uintptr_t sample_include_word = sample_include[widx];
          const uintptr_t dphase_present_word = raw_dphase_present[widx];
          uintptr_t dosage_present_word = raw_dosage_present[widx];
          if (!sample_include_word) {
            dosage_main_iter = &(dosage_main_iter[PopcountWord(dosage_present_word)]);
            dphase_delta_iter = &(dphase_delta_iter[PopcountWord(dphase_present_word)]);
            continue;
          }
          const uintptr_t phasepresent_word = raw_phasepresent[widx];
          phased_hc_het_ct += PopcountWord(sample_include_word & phasepresent_word & (~dosage_present_word));
          while (dosage_present_word) {
            const uintptr_t lowbit = dosage_present_word & (-dosage_present_word);
            const uintptr_t dphase_here = lowbit & dphase_present_word;
            if (lowbit & sample_include_word) {
              const uintptr_t cur_dosage_val = *dosage_main_iter;
              alt1_dosage += cur_dosage_val;
              hap_ssq_x2 += cur_dosage_val * cur_dosage_val;
              ++dosage_ct;
              if (dphase_here) {
                const intptr_t dphase_delta_val = *dphase_delta_iter;
                hap_ssq_x2 += dphase_delta_val * dphase_delta_val;
              } else if (lowbit & phasepresent_word) {
                const uintptr_t homdist = 16384 - abs_i32(16384 - cur_dosage_val);
                hap_ssq_x2 += homdist * homdist;
              }
            }
            dphase_delta_iter += (dphase_here != 0);
            dosage_present_word ^= lowbit;
            ++dosage_main_iter;
          }
        }
        GenoarrCountSubsetIntersectFreqs(raw_genovec, raw_dosage_present, sample_include, raw_sample_ct, replaced_genocounts);
      }
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
  hap_ssq_x2 += (remaining_het_ct + phased_hc_het_ct) * 0x10000000LLU + remaining_hom_alt_ct * 0x40000000LLU;
  *imp_r2_ptr = BiallelicDiploidMinimac3R2(alt1_dosage, hap_ssq_x2, new_sample_nm_ct);
  return kPglRetSuccess;
}

PglErr PgrGetDCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t is_minimac3_r2, PgenReader* pgrp, double* imp_r2_ptr, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* __restrict all_dosages) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    STD_ARRAY_REF_FILL0(4, genocounts);
    all_dosages[0] = 0;
    all_dosages[1] = 0;
    if (imp_r2_ptr) {
      *imp_r2_ptr = 0.0 / 0.0;
    }
    return kPglRetSuccess;
  }
  return GetBasicGenotypeCountsAndDosage16s(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, vidx, is_minimac3_r2, pgrp, imp_r2_ptr, genocounts, all_dosages);
}

// Does not zero-initialize results[].
void CountAllBytes64(const void* bytearr, uintptr_t byte_ct, uint64_t* __restrict results) {
  const unsigned char* bytearr_uc = S_CAST(const unsigned char*, bytearr);
  for (uintptr_t ulii = 0; ulii != byte_ct; ++ulii) {
    results[bytearr_uc[ulii]] += 1;
  }
}

// Does not zero-initialize results[].
void CountAllNibbles64(const void* nibblearr, uintptr_t nibble_ct, uint64_t* __restrict results) {
  // possible todo: for sufficiently large nibble_ct, use CountAllBytes and
  // then postprocess
  const uintptr_t fullbyte_ct = nibble_ct / 2;
  const unsigned char* nibblearr_uc = S_CAST(const unsigned char*, nibblearr);
  for (uintptr_t ulii = 0; ulii != fullbyte_ct; ++ulii) {
    const uint32_t uii = nibblearr_uc[ulii];
    results[uii & 15] += 1;
    results[uii >> 4] += 1;
  }
  if (nibble_ct % 2) {
    results[nibblearr_uc[fullbyte_ct] & 15] += 1;
  }
}

void CountAllAux1aDense(const void* patch_01_fvals, uint32_t allele_ct, uint32_t rare01_ct, uint64_t* __restrict one_cts) {
  one_cts[1] -= rare01_ct;
  if (allele_ct < 5) {
    if (allele_ct == 3) {
      // all entries are 0/1 -> 0/2
      one_cts[2] = rare01_ct;
      return;
    }
    const uint32_t allele_code_byte_ct = DivUp(rare01_ct, 8);
    const uint32_t alt3_ct = PopcountBytes(patch_01_fvals, allele_code_byte_ct);
    one_cts[2] = rare01_ct - alt3_ct;
    one_cts[3] = alt3_ct;
    return;
  }
  if (allele_ct < 19) {
    if (allele_ct < 7) {
      STD_ARRAY_DECL(uint32_t, 4, rare0het_counts);
      GenovecCountFreqs(R_CAST(const uintptr_t*, patch_01_fvals), rare01_ct, rare0het_counts);
      for (uint32_t allele_idx_p2 = 2; allele_idx_p2 != allele_ct; ++allele_idx_p2) {
        one_cts[allele_idx_p2] = rare0het_counts[allele_idx_p2 - 2];
      }
      return;
    }
    CountAllNibbles64(patch_01_fvals, rare01_ct, &(one_cts[2]));
    return;
  }
  CountAllBytes64(patch_01_fvals, rare01_ct, &(one_cts[2]));
}

// assumes one_cts[1] initialized to genocounts[1]
// sample_include should be nullptr if we aren't subsetting
PglErr CountAllAux1a(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uintptr_t* __restrict raw_genovec, uint32_t aux1a_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t raw_01_ct, const unsigned char** fread_pp, uint64_t* __restrict one_cts, uint32_t* __restrict deltalist_workspace) {
  if (aux1a_mode == 15) {
    return kPglRetSuccess;
  }
  if (!sample_include) {
    uint32_t rare01_ct;
    if (!aux1a_mode) {
      const uint32_t fset_byte_ct = DivUp(raw_01_ct, CHAR_BIT);
      rare01_ct = PopcountBytes(*fread_pp, fset_byte_ct);
      *fread_pp += fset_byte_ct;
    } else {
      const unsigned char* group_info_iter;
      PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, &rare01_ct);
      if (unlikely(reterr)) {
        return reterr;
      }
      reterr = SkipDeltalistIds(fread_end, group_info_iter, rare01_ct, raw_sample_ct, 0, fread_pp);
      if (unlikely(reterr)) {
        return reterr;
      }
    }
    const unsigned char* patch_01_fvals = *fread_pp;
    const uint32_t fvals_byte_ct = GetAux1aAlleleEntryByteCt(allele_ct, rare01_ct);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    CountAllAux1aDense(patch_01_fvals, allele_ct, rare01_ct, one_cts);
    return kPglRetSuccess;
  }
  const uint32_t allele_code_width = GetAux1aWidth(allele_ct);
  const uintptr_t allele_code_mask = (1U << allele_code_width) - 1;
  uint64_t* one_cts_offset2 = &(one_cts[2]);
  if (!aux1a_mode) {
    const uint32_t fset_byte_ct = DivUp(raw_01_ct, CHAR_BIT);
    const uint32_t rare01_ct = PopcountBytes(*fread_pp, fset_byte_ct);
#ifdef __arm__
#  error "Unaligned accesses in CountAllAux1a()."
#endif
    const uintptr_t* patch_01_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    *fread_pp += fset_byte_ct;
    const uintptr_t* patch_01_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const Halfword* sample_include_hw = R_CAST(const Halfword*, sample_include);
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_hets = Word01(raw_genovec[0]);
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t subsetted_rare01_ct = 0;
    uint32_t loop_len = kBitsPerWord;
    uint32_t rare01_lowbits = kBitsPerWord;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          break;
        }
        fset_bits = SubwordLoad(&(patch_01_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_01_ct, kBitsPerWord);
      } else {
        fset_bits = patch_01_fsetw[fset_widx];
      }
      if (allele_ct == 3) {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_hets) {
            cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_hets) / 2;
            subsetted_rare01_ct += (sample_include_hw[sample_hwidx] >> sample_uidx_lowbits) & 1;
          }
          cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
          fset_bits = fset_bits >> 1;
        }
      } else {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_hets) {
            cur_raw_genovec_hets = Word01(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare01_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_01_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_01_fvalsw[fvals_widx];
              }
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare01_lowbits = 0;
            }
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_hets) / 2;
            if (sample_include_hw[sample_hwidx] & (1U << sample_uidx_lowbits)) {
              ++subsetted_rare01_ct;
              one_cts_offset2[(fvals_bits >> rare01_lowbits) & allele_code_mask] += 1;
            }
            rare01_lowbits += allele_code_width;
          }
          cur_raw_genovec_hets &= cur_raw_genovec_hets - 1;
          fset_bits = fset_bits >> 1;
        }
      }
    }
    one_cts_offset2[-1] -= subsetted_rare01_ct;
    if (allele_ct == 3) {
      one_cts_offset2[0] = subsetted_rare01_ct;
    }
    return kPglRetSuccess;
  }
  // mode 1: difflist.
  if (allele_ct == 3) {
    // Use CountDeltalistIntersect shortcut here.
    uint32_t subsetted_02_ct;
    uint32_t rare01_ct;
    PglErr reterr = CountDeltalistIntersect(fread_end, sample_include, raw_sample_ct, fread_pp, &subsetted_02_ct, &rare01_ct);
    if (unlikely(reterr)) {
      return reterr;
    }
    one_cts_offset2[-1] -= subsetted_02_ct;
    one_cts_offset2[0] = subsetted_02_ct;
    return kPglRetSuccess;
  }
  // Save deltalist elements, iterate.
  uint32_t rare01_ct;
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare01_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_01_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare01_ct) * allele_code_width, 8);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  const uint32_t allele_code_logwidth = ctzu32(allele_code_width);
  uint32_t subsetted_rare01_ct = 0;
  uint32_t loop_len = kBitsPerWord >> allele_code_logwidth;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        break;
      }
      fvals_bits = SubwordLoad(&(patch_01_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
      loop_len = 1 + ((rare01_ct - 1) & (loop_len - 1));
    } else {
      fvals_bits = patch_01_fvalsw[fvals_widx];
    }
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - allele_code_logwidth)]);
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uint32_t sample_uidx = cur_deltalist_base[uii];
      if (IsSet(sample_include, sample_uidx)) {
        ++subsetted_rare01_ct;
        one_cts_offset2[(fvals_bits >> (uii << allele_code_logwidth)) & allele_code_mask] += 1;
      }
    }
  }
  one_cts_offset2[-1] -= subsetted_rare01_ct;
  return kPglRetSuccess;
}

void CountAllAux1bDense(const void* __restrict patch_10_fvals, uint32_t allele_ct, uint32_t rare10_ct, uint64_t* __restrict one_cts_offset1, uint64_t* __restrict two_cts_offset1) {
  // probable todo: faster path if two_cts_offset1 == nullptr
  const uint32_t allele_ct_m1 = allele_ct - 1;
  two_cts_offset1[0] -= rare10_ct;
  if (allele_ct_m1 < 5) {
    if (allele_ct_m1 == 2) {
      const uint32_t allele_code_byte_ct = DivUp(rare10_ct, 8);
      const uint32_t hom22_ct = PopcountBytes(patch_10_fvals, allele_code_byte_ct);
      const uint32_t het12_ct = rare10_ct - hom22_ct;
      one_cts_offset1[0] += het12_ct;
      one_cts_offset1[1] += het12_ct;
      two_cts_offset1[1] = hom22_ct;
      return;
    }
    STD_ARRAY_DECL(uint32_t, 4, alt_counts);
    GenovecCountFreqs(R_CAST(const uintptr_t*, patch_10_fvals), rare10_ct * 2, alt_counts);
    one_cts_offset1[0] += alt_counts[0];
    for (uint32_t allele_idx_m1 = 1; allele_idx_m1 != allele_ct_m1; ++allele_idx_m1) {
      const uint32_t homxx_ct = CountNibble(patch_10_fvals, allele_idx_m1 * kMask5555, rare10_ct);
      one_cts_offset1[allele_idx_m1] += alt_counts[allele_idx_m1] - 2 * homxx_ct;
      two_cts_offset1[allele_idx_m1] = homxx_ct;
    }
    return;
  }
  const unsigned char* patch_10_fvals_uc = S_CAST(const unsigned char*, patch_10_fvals);
  if (allele_ct_m1 < 17) {
    // for larger rare10_ct, this should use a byte counter
    for (uint32_t uii = 0; uii != rare10_ct; ++uii) {
      const uint32_t cur_byte = patch_10_fvals_uc[uii];
      const uint32_t cur_byte_hi = cur_byte >> 4;
      const uint32_t cur_byte_lo = cur_byte & 15;
      if (cur_byte_hi == cur_byte_lo) {
        two_cts_offset1[cur_byte_lo] += 1;
      } else {
        one_cts_offset1[cur_byte_lo] += 1;
        one_cts_offset1[cur_byte_hi] += 1;
      }
    }
    return;
  }
  for (uint32_t uii = 0; uii != rare10_ct; ++uii) {
    const uint32_t cur_byte_lo = patch_10_fvals_uc[2 * uii];
    const uint32_t cur_byte_hi = patch_10_fvals_uc[2 * uii + 1];
    if (cur_byte_hi == cur_byte_lo) {
      two_cts_offset1[cur_byte_lo] += 1;
    } else {
      one_cts_offset1[cur_byte_lo] += 1;
      one_cts_offset1[cur_byte_hi] += 1;
    }
  }
}

PglErr CountAllAux1b(const unsigned char* fread_end, const uintptr_t* __restrict sample_include, const uintptr_t* __restrict raw_genovec, uint32_t aux1b_mode, uint32_t raw_sample_ct, uint32_t allele_ct, uint32_t raw_10_ct, const unsigned char** fread_pp, uint64_t* __restrict one_cts, uint64_t* __restrict two_cts, uint32_t* __restrict deltalist_workspace) {
  if (aux1b_mode == 15) {
    return kPglRetSuccess;
  }
  uint64_t* one_cts_offset1 = &(one_cts[1]);
  uint64_t* two_cts_offset1 = &(two_cts[1]);
  if (!sample_include) {
    uint32_t rare10_ct;
    if (!aux1b_mode) {
      const uint32_t fset_byte_ct = DivUp(raw_10_ct, CHAR_BIT);
      rare10_ct = PopcountBytes(*fread_pp, fset_byte_ct);
      *fread_pp += fset_byte_ct;
    } else {
      const unsigned char* group_info_iter;
      PglErr reterr = ParseDifflistHeader(fread_end, raw_sample_ct, fread_pp, nullptr, &group_info_iter, &rare10_ct);
      if (unlikely(reterr)) {
        return reterr;
      }
      reterr = SkipDeltalistIds(fread_end, group_info_iter, rare10_ct, raw_sample_ct, 0, fread_pp);
      if (unlikely(reterr)) {
        return reterr;
      }
    }
    const unsigned char* patch_10_fvals = *fread_pp;
    const uint32_t fvals_byte_ct = GetAux1bAlleleEntryByteCt(allele_ct, rare10_ct);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    CountAllAux1bDense(patch_10_fvals, allele_ct, rare10_ct, one_cts_offset1, two_cts_offset1);
    return kPglRetSuccess;
  }
  uintptr_t detect_hom_mask_lo;  // unused
  const uint32_t allele_code_logwidth = GetAux1bConsts(allele_ct, &detect_hom_mask_lo);
  const uint32_t code10_logwidth = allele_code_logwidth + (allele_code_logwidth != 0);
  const uint32_t allele_code_width = 1U << allele_code_logwidth;
  const uint32_t allele_code_mask = (1U << allele_code_width) - 1;
  const uint32_t allele_ct_m1 = allele_ct - 1;
  uint32_t rare10_lowbits = kBitsPerWord;
  // probable todo: faster paths when two_cts_offset1 == nullptr
  if (!aux1b_mode) {
    const uint32_t fset_byte_ct = DivUp(raw_10_ct, CHAR_BIT);
    const uint32_t rare10_ct = PopcountBytes(*fread_pp, fset_byte_ct);
#ifdef __arm__
#  error "Unaligned accesses in CountAllAux1b()."
#endif
    const uintptr_t* patch_10_fsetw = R_CAST(const uintptr_t*, *fread_pp);
    *fread_pp += fset_byte_ct;
    const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) << code10_logwidth, 8);
    if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
      return kPglRetMalformedInput;
    }
    const Halfword* sample_include_hw = R_CAST(const Halfword*, sample_include);
    uintptr_t sample_hwidx = 0;
    uintptr_t cur_raw_genovec_xys = Word10(raw_genovec[0]);
    const uint32_t fset_word_ct_m1 = (fset_byte_ct - 1) / kBytesPerWord;
    const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
    const uint32_t code10_width = 1U << code10_logwidth;
    uintptr_t fvals_bits = 0;
    uint32_t fvals_widx = 0;
    uint32_t subsetted_rare10_ct = 0;
    uint32_t loop_len = kBitsPerWord;
    for (uint32_t fset_widx = 0; ; ++fset_widx) {
      uintptr_t fset_bits;
      if (fset_widx >= fset_word_ct_m1) {
        if (fset_widx > fset_word_ct_m1) {
          break;
        }
        fset_bits = SubwordLoad(&(patch_10_fsetw[fset_word_ct_m1]), ModNz(fset_byte_ct, kBytesPerWord));
        loop_len = ModNz(raw_10_ct, kBitsPerWord);
      } else {
        fset_bits = patch_10_fsetw[fset_widx];
      }
      if (allele_ct_m1 == 2) {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_xys) {
            cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare10_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_10_fvalsw[fvals_widx];
              }
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare10_lowbits = 0;
            }
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
            if (sample_include_hw[sample_hwidx] & (1U << sample_uidx_lowbits)) {
              ++subsetted_rare10_ct;
              two_cts_offset1[1] += (fvals_bits >> rare10_lowbits) & 1;
            }
            ++rare10_lowbits;
          }
          cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
          fset_bits = fset_bits >> 1;
        }
      } else {
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          while (!cur_raw_genovec_xys) {
            cur_raw_genovec_xys = Word10(raw_genovec[++sample_hwidx]);
          }
          if (fset_bits & 1) {
            if (rare10_lowbits == kBitsPerWord) {
              if (fvals_widx == fvals_word_ct_m1) {
                fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
              } else {
                fvals_bits = patch_10_fvalsw[fvals_widx];
              }
              // unnecessary to apply bzhi here
              ++fvals_widx;
              rare10_lowbits = 0;
            }
            const uint32_t sample_uidx_lowbits = ctzw(cur_raw_genovec_xys) / 2;
            if (sample_include_hw[sample_hwidx] & (1U << sample_uidx_lowbits)) {
              ++subsetted_rare10_ct;
              const uintptr_t cur_code_pair = fvals_bits >> rare10_lowbits;
              const uint32_t cur_code_hi = (cur_code_pair >> allele_code_width) & allele_code_mask;
              const uint32_t cur_code_lo = cur_code_pair & allele_code_mask;
              if (cur_code_hi == cur_code_lo) {
                two_cts_offset1[cur_code_lo] += 1;
              } else {
                one_cts_offset1[cur_code_lo] += 1;
                one_cts_offset1[cur_code_hi] += 1;
              }
            }
            rare10_lowbits += code10_width;
          }
          cur_raw_genovec_xys &= cur_raw_genovec_xys - 1;
          fset_bits = fset_bits >> 1;
        }
      }
    }
    two_cts_offset1[0] -= subsetted_rare10_ct;
    if (allele_ct == 3) {
      const uint32_t subsetted_het12_ct = subsetted_rare10_ct - two_cts_offset1[1];
      one_cts_offset1[0] += subsetted_het12_ct;
      one_cts_offset1[1] += subsetted_het12_ct;
    }
    return kPglRetSuccess;
  }
  // Save deltalist elements, iterate.
  uint32_t rare10_ct;
  PglErr reterr = ParseAndSaveDeltalist(fread_end, raw_sample_ct, fread_pp, deltalist_workspace, &rare10_ct);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uintptr_t* patch_10_fvalsw = R_CAST(const uintptr_t*, *fread_pp);
  const uint32_t fvals_byte_ct = DivUpU64(S_CAST(uint64_t, rare10_ct) << code10_logwidth, 8);
  if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  const uint32_t fvals_word_ct_m1 = (fvals_byte_ct - 1) / kBytesPerWord;
  uint32_t subsetted_rare10_ct = 0;
  uint32_t loop_len = kBitsPerWord >> code10_logwidth;
  for (uint32_t fvals_widx = 0; ; ++fvals_widx) {
    uintptr_t fvals_bits;
    if (fvals_widx >= fvals_word_ct_m1) {
      if (fvals_widx > fvals_word_ct_m1) {
        break;
      }
      fvals_bits = SubwordLoad(&(patch_10_fvalsw[fvals_widx]), ModNz(fvals_byte_ct, kBytesPerWord));
      loop_len = 1 + ((rare10_ct - 1) & (loop_len - 1));
    } else {
      fvals_bits = patch_10_fvalsw[fvals_widx];
    }
    const uint32_t* cur_deltalist_base = &(deltalist_workspace[fvals_widx << (kBitsPerWordLog2 - code10_logwidth)]);
    if (allele_ct == 3) {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uint32_t sample_uidx = cur_deltalist_base[uii];
        if (IsSet(sample_include, sample_uidx)) {
          ++subsetted_rare10_ct;
          two_cts_offset1[1] += (fvals_bits >> uii) & 1;
        }
      }
    } else {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uint32_t sample_uidx = cur_deltalist_base[uii];
        if (IsSet(sample_include, sample_uidx)) {
          ++subsetted_rare10_ct;
          const uintptr_t cur_code_pair = fvals_bits >> (uii << code10_logwidth);
          const uint32_t cur_code_hi = (cur_code_pair >> allele_code_width) & allele_code_mask;
          const uint32_t cur_code_lo = cur_code_pair & allele_code_mask;
          if (cur_code_hi == cur_code_lo) {
            two_cts_offset1[cur_code_lo] += 1;
          } else {
            one_cts_offset1[cur_code_lo] += 1;
            one_cts_offset1[cur_code_hi] += 1;
          }
        }
      }
    }
  }
  two_cts_offset1[0] -= subsetted_rare10_ct;
  if (allele_ct == 3) {
    const uint32_t subsetted_het12_ct = subsetted_rare10_ct - two_cts_offset1[1];
    one_cts_offset1[0] += subsetted_het12_ct;
    one_cts_offset1[1] += subsetted_het12_ct;
  }
  return kPglRetSuccess;
}

PglErr GetMultiallelicCountsAndDosage16s(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t sample_ct, uint32_t vidx, uint32_t allele_ct, __maybe_unused uint32_t is_minimac3_r2, PgenReader* pgrp, double* __restrict imp_r2_ptr, uint32_t* __restrict het_ctp, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* all_dosages) {
  // only called on multiallelic variants
  // no dosages for now
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  uintptr_t* raw_genovec = pgrp->workspace_vec;
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRawGenovec(subsetting_required, vidx, pgrp, &fread_ptr, &fread_end, raw_genovec);
  if (unlikely(reterr)) {
    return reterr;
  }
  ZeroTrailingQuaters(raw_sample_ct, raw_genovec);
  if (!subsetting_required) {
    GenovecCountFreqsUnsafe(raw_genovec, raw_sample_ct, genocounts);
    sample_include = nullptr;
  } else {
    GenovecCountSubsetFreqs(raw_genovec, sample_include_interleaved_vec, raw_sample_ct, sample_ct, genocounts);
  }
  uint64_t* one_cts = pgrp->workspace_imp_r2;
  uint64_t* two_cts = &(one_cts[allele_ct]);
  one_cts[0] = genocounts[1];
  one_cts[1] = genocounts[1];
  ZeroU64Arr(allele_ct - 2, &(one_cts[2]));
  two_cts[0] = genocounts[0];
  two_cts[1] = genocounts[2];
  ZeroU64Arr(allele_ct - 2, &(two_cts[2]));
  // Cases:
  // - No hardcall-phase present.  Then we don't need to know raw_het_ct.
  // - No multiallelic dosages present, not computing minimac3-r2.  Then we
  //   still don't need to know raw_het_ct.
  // - Otherwise, we need to know raw_het_ct, either for the minimac3-r2
  //   computation or to locate the beginning of aux3/aux4.
  //   If we're computing minimac3-r2, AND
  //     (i) we're subsetting, or
  //     (ii) multiallelic dosages are present,
  //   it's also necessary to compute all_hets, either to compute correct
  //   subsetted minimac3-r2 or to know how many phased-hardcalls are
  //   overridden by phased dosages.
  const uint32_t raw_het_ct_needed = VrtypeHphase(vrtype) && (is_minimac3_r2 || (vrtype & 0x60));
  uintptr_t* all_hets = nullptr;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  uint32_t raw_het_ct = genocounts[1]; // inaccurate, corrected later if needed
  if (VrtypeMultiallelicHc(vrtype)) {
    const uint32_t aux1_first_byte = *fread_ptr++;
    const uint32_t aux1a_mode = aux1_first_byte & 15;
    const uint32_t aux1b_mode = aux1_first_byte >> 4;
    uint32_t raw_10_ct = 0;
    if ((!aux1a_mode) || (!aux1b_mode) || sample_include) {
      GenovecCount12Unsafe(raw_genovec, raw_sample_ct, &raw_het_ct, &raw_10_ct);
    }
    uint32_t* deltalist_workspace = pgrp->workspace_difflist_sample_ids;
    reterr = CountAllAux1a(fread_end, sample_include, raw_genovec, aux1a_mode, raw_sample_ct, allele_ct, raw_het_ct, &fread_ptr, one_cts, deltalist_workspace);
    if (unlikely(reterr)) {
      return reterr;
    }
    const unsigned char* aux1b_start = fread_ptr;
    reterr = CountAllAux1b(fread_end, sample_include, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &fread_ptr, one_cts, two_cts, deltalist_workspace);
    if (unlikely(reterr)) {
      return reterr;
    }
    if (raw_het_ct_needed) {
      if (!sample_include) {
        raw_het_ct += genocounts[2];
        for (uint32_t aidx = 1; aidx != allele_ct; ++aidx) {
          raw_het_ct -= two_cts[aidx];
        }
      }
      if (sample_include || (is_minimac3_r2 && (vrtype & 0x60))) {
        all_hets = pgrp->workspace_all_hets;
        PgrDetectGenovecHets(raw_genovec, raw_sample_ct, all_hets);
        if (aux1b_mode != 15) {
          uintptr_t* aux1b_hets = pgrp->workspace_aux1x_present;
          uint32_t aux1b_het_present;
          reterr = GetAux1bHets(fread_end, raw_genovec, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &aux1b_start, aux1b_hets, &aux1b_het_present, deltalist_workspace);
          if (unlikely(reterr)) {
            return reterr;
          }
          if (aux1b_het_present) {
            BitvecOr(aux1b_hets, raw_sample_ctl, all_hets);
          }
        }
        if (sample_include) {
          raw_het_ct = PopcountWords(all_hets, raw_sample_ctl);
        }
      }
    }
  }
  uintptr_t* raw_phasepresent = nullptr;
  uint32_t extra_phased_het_ct = 0;
  if (raw_het_ct_needed) {
    if (!all_hets) {
      reterr = SkipAux2(fread_end, raw_het_ct, &fread_ptr, is_minimac3_r2? (&extra_phased_het_ct) : nullptr);
      if (unlikely(reterr)) {
        return reterr;
      }
    } else {
      raw_phasepresent = pgrp->workspace_subset;
      reterr = GetPhasepresentAndSkipPhaseinfo(fread_end, all_hets, raw_sample_ct, raw_het_ct, &fread_ptr, raw_phasepresent, &extra_phased_het_ct);
      if (unlikely(reterr)) {
        return reterr;
      }
      if (sample_include) {
        extra_phased_het_ct = PopcountWordsIntersect(raw_phasepresent, sample_include, raw_sample_ctl);
      }
    }
  }
  if (!(vrtype & 0x60)) {
    uint32_t hom_hc_ct = 0;
    for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
      const uint64_t cur_hom_ct = two_cts[allele_idx];
      hom_hc_ct += cur_hom_ct;
      const uint64_t two_dosage = cur_hom_ct * 0x8000LLU;
      const uint64_t dosage_sum = one_cts[allele_idx] * 0x4000LLU + two_dosage;
      all_dosages[allele_idx] = dosage_sum;
      // Repurpose two_cts[] to store ssqs.
      two_cts[allele_idx] = (dosage_sum + two_dosage) * 0x4000LLU;
    }
    const uint32_t nm_sample_ct = sample_ct - genocounts[3];
    *het_ctp = nm_sample_ct - hom_hc_ct;
    if (!imp_r2_ptr) {
      return kPglRetSuccess;
    }
    *imp_r2_ptr = MultiallelicDiploidMinimac3R2(all_dosages, two_cts, nm_sample_ct, allele_ct, extra_phased_het_ct);
    if (!is_minimac3_r2) {
      *imp_r2_ptr *= 2;
    }
    return kPglRetSuccess;
  }
  fputs("dosages not yet supported by GetMultiallelicCountsAndDosage16s()\n", stderr);
  exit(S_CAST(int32_t, kPglRetNotYetSupported));
  return kPglRetNotYetSupported;
}

PglErr PgrGetMDCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t is_minimac3_r2, PgenReader* pgrp, double* __restrict imp_r2_ptr, uint32_t* __restrict het_ctp, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* __restrict all_dosages) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  if (!sample_ct) {
    STD_ARRAY_REF_FILL0(4, genocounts);
    ZeroU64Arr(allele_ct, all_dosages);
    if (imp_r2_ptr) {
      *imp_r2_ptr = 0.0 / 0.0;
    }
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  if ((allele_ct == 2) || (!(vrtype & 0x68))) {
    PglErr reterr = GetBasicGenotypeCountsAndDosage16s(sample_include, sample_include_interleaved_vec, sample_include_cumulative_popcounts, sample_ct, vidx, is_minimac3_r2, pgrp, imp_r2_ptr, genocounts, all_dosages);
    *het_ctp = genocounts[1];
    ZeroU64Arr(allele_ct - 2, &(all_dosages[2]));
    return reterr;
  }
  return GetMultiallelicCountsAndDosage16s(sample_include, sample_include_interleaved_vec, sample_ct, vidx, allele_ct, is_minimac3_r2, pgrp, imp_r2_ptr, het_ctp, genocounts, all_dosages);
}

PglErr PgrGetMD(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp) {
  pgvp->patch_01_ct = 0;
  pgvp->patch_10_ct = 0;
  pgvp->dosage_ct = 0;
  pgvp->multidosage_sample_ct = 0;
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  if ((allele_ct == 2) || (!(vrtype & 0x68))) {
    return PgrGetD(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, pgvp->genovec, pgvp->dosage_present, pgvp->dosage_main, &(pgvp->dosage_ct));
  }
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  uintptr_t* all_hets = VrtypeHphase(vrtype)? pgrp->workspace_all_hets : nullptr;
  if (VrtypeMultiallelicHc(vrtype)) {
    PglErr reterr = GetMultiallelicCodes(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, all_hets? (&fread_ptr) : nullptr, all_hets? (&fread_end) : nullptr, all_hets, pgvp);
    if (!(vrtype & 0x60)) {
      return reterr;
    }
  } else {
    // todo: ReadRawGenovec, etc.
  }
  fputs("true multiallelic dosages not yet supported by PgrGetMD()\n", stderr);
  exit(S_CAST(int32_t, kPglRetNotYetSupported));
  return kPglRetSuccess;
}

PglErr PgrGetDp(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp) {
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    pgvp->phasepresent_ct = 0;
    pgvp->dosage_ct = 0;
    pgvp->dphase_ct = 0;
    return kPglRetSuccess;
  }
  const unsigned char* fread_ptr = nullptr;
  const unsigned char* fread_end = nullptr;
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t dosage_is_present = VrtypeDosage(vrtype);
  PglErr reterr = ReadGenovecHphaseSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, dosage_is_present? (&fread_ptr) : nullptr, dosage_is_present? (&fread_end) : nullptr, pgvp->genovec, pgvp->phasepresent, pgvp->phaseinfo, &(pgvp->phasepresent_ct));
  if (reterr || (!dosage_is_present)) {
    pgvp->dosage_ct = 0;
    pgvp->dphase_ct = 0;
    return reterr;
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  return ParseDosage16(fread_ptr, fread_end, sample_include, sample_ct, vidx, allele_ct, pgrp, &(pgvp->dosage_ct), pgvp->dphase_present, pgvp->dphase_delta, &(pgvp->dphase_ct), pgvp->dosage_present, pgvp->dosage_main);
}

PglErr PgrGetInv1Dp(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReader* pgrp, PgenVariant* pgvp) {
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  if ((allele_ct == 2) || (!allele_idx)) {
    PglErr reterr = PgrGetDp(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, pgvp);
    if (allele_idx) {
      GenovecInvertUnsafe(sample_ct, pgvp->genovec);
      if (pgvp->phasepresent_ct) {
        BitvecInvert(BitCtToWordCt(sample_ct), pgvp->phaseinfo);
      }
      if (pgvp->dosage_ct) {
        BiallelicDosage16Invert(pgvp->dosage_ct, pgvp->dosage_main);
        if (pgvp->dphase_ct) {
          BiallelicDphase16Invert(pgvp->dphase_ct, pgvp->dphase_delta);
        }
      }
    }
    return reterr;
  }
  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  if (!VrtypeDosage(vrtype)) {
    pgvp->dosage_ct = 0;
    pgvp->dphase_ct = 0;
    return PgrGetInv1P(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, pgvp->genovec, pgvp->phasepresent, pgvp->phaseinfo, &(pgvp->phasepresent_ct));
  }
  fputs("multiallelic dosage not yet supported by PgrGetInv1Dp()\n", stderr);
  exit(S_CAST(int32_t, kPglRetNotYetSupported));
  return kPglRetSuccess;
}

PglErr PgrGetMDp(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp) {
  pgvp->patch_01_ct = 0;
  pgvp->patch_10_ct = 0;
  pgvp->phasepresent_ct = 0;
  pgvp->dosage_ct = 0;
  pgvp->multidosage_sample_ct = 0;
  pgvp->dphase_ct = 0;
  pgvp->multidphase_sample_ct = 0;
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  if ((allele_ct == 2) || (!(vrtype & 0x68))) {
    return PgrGetDp(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, pgvp);
  }
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  uintptr_t* all_hets = VrtypeHphase(vrtype)? pgrp->workspace_all_hets : nullptr;
  if (VrtypeMultiallelicHc(vrtype)) {
    PglErr reterr = GetMultiallelicCodes(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, all_hets? (&fread_ptr) : nullptr, all_hets? (&fread_end) : nullptr, all_hets, pgvp);
    if (reterr || (!all_hets)) {
      return reterr;
    }
    if (!(vrtype & 0x60)) {
      const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
      return ParseAux2Subset(fread_end, (sample_ct != raw_sample_ct)? sample_include : nullptr, all_hets, nullptr, raw_sample_ct, sample_ct, &fread_ptr, pgvp->phasepresent, pgvp->phaseinfo, &(pgvp->phasepresent_ct), pgrp->workspace_subset);
    }
  } else {
    // todo: ReadRawGenovec, etc.
  }
  fputs("true multiallelic dosages not yet supported by PgrGetMDp()\n", stderr);
  fprintf(stderr, "%u\n", vidx);
  exit(S_CAST(int32_t, kPglRetNotYetSupported));
  return kPglRetSuccess;

}

static_assert(sizeof(AlleleCode) == 1, "CountAux1bHets() must be updated.");
uintptr_t CountAux1bHets(const AlleleCode* patch_10_vals, uintptr_t rare10_ct) {
  // Similar to CountByte().
  uintptr_t byte_ct = rare10_ct * 2;
#ifdef __LP64__
  if (byte_ct < kBytesPerVec) {
#endif
    uintptr_t tot = 0;
    for (uintptr_t offset = 0; offset < byte_ct; offset += 2) {
      tot += (patch_10_vals[offset] != patch_10_vals[offset + 1]);
    }
    return tot;
#ifdef __LP64__
  }
  const unsigned char* bytearr_uc_iter = R_CAST(const unsigned char*, patch_10_vals);
  const VecW m0 = vecw_setzero();
  const VecW m8 = VCONST_W(kMask00FF);
  VecW acc = vecw_setzero();
  while (byte_ct > 255 * kBytesPerVec) {
    VecUc inner_acc = vecuc_setzero();
    for (uint32_t uii = 0; uii != 255; ++uii) {
      const VecUc cur_vvec = vecuc_loadu(bytearr_uc_iter);
      bytearr_uc_iter = &(bytearr_uc_iter[kBytesPerVec]);
      const VecUc shifted_vvec = R_CAST(const VecUc, vecw_srli(R_CAST(const VecW, cur_vvec), 8));
      inner_acc = inner_acc - (cur_vvec == shifted_vvec);
    }
    const VecW partial_sums = R_CAST(VecW, inner_acc) & m8;
    acc = acc + vecw_sad(partial_sums, m0);
    byte_ct -= 255 * kBytesPerVec;
  }
  const unsigned char* bytearr_uc_final = &(bytearr_uc_iter[byte_ct - kBytesPerVec]);
  VecUc inner_acc = vecuc_setzero();
  while (bytearr_uc_iter < bytearr_uc_final) {
    const VecUc cur_vvec = vecuc_loadu(bytearr_uc_iter);
    bytearr_uc_iter = &(bytearr_uc_iter[kBytesPerVec]);
    const VecUc shifted_vvec = R_CAST(const VecUc, vecw_srli(R_CAST(const VecW, cur_vvec), 8));
    inner_acc = inner_acc - (cur_vvec == shifted_vvec);
  }
  VecUc cur_vvec = vecuc_loadu(bytearr_uc_final);
  const uintptr_t overlap_byte_ct = bytearr_uc_iter - bytearr_uc_final;
  const VecUc shifted_vvec = R_CAST(const VecUc, vecw_srli(R_CAST(const VecW, cur_vvec), 8));
  const VecUc mask_vvec = vecuc_loadu(&(kLeadMask[kBytesPerVec - overlap_byte_ct]));
  cur_vvec = (cur_vvec == shifted_vvec) & mask_vvec;
  inner_acc = inner_acc - cur_vvec;
  const VecW partial_sums = R_CAST(VecW, inner_acc) & m8;
  acc = acc + vecw_sad(partial_sums, m0);
  const uintptr_t tot = HsumW(acc);
  return rare10_ct - tot;
#endif
}

PglErr PgrGetRaw(uint32_t vidx, PgenGlobalFlags read_gflags, PgenReader* pgrp, uintptr_t** loadbuf_iter_ptr, unsigned char* loaded_vrtype_ptr) {
  // currently handles multiallelic hardcalls, hardcall phase, and biallelic
  // dosage (both unphased and phased)
  // todo: multiallelic dosage
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  uintptr_t* genovec = (*loadbuf_iter_ptr);
  uintptr_t* loadbuf_iter = &(genovec[QuaterCtToAlignedWordCt(raw_sample_ct)]);
  const uint32_t multiallelic_hc_present = (vrtype / 8) & 1;
  const uint32_t save_multiallelic_hc = multiallelic_hc_present && (read_gflags & kfPgenGlobalMultiallelicHardcallFound);
  const uint32_t hphase_is_present = (vrtype / 0x10) & 1;
  const uint32_t save_hphase = hphase_is_present && (read_gflags & kfPgenGlobalHardcallPhasePresent);
  const uint32_t dosage_is_present = (vrtype & 0x60)? 1 : 0;
  const uint32_t save_dosage = dosage_is_present && (read_gflags & kfPgenGlobalDosagePresent);

  const uint32_t save_dphase = (vrtype & 0x80) && (read_gflags & kfPgenGlobalDosagePhasePresent);
  assert(save_dosage || (!save_dphase));

  if (loaded_vrtype_ptr) {
    *loaded_vrtype_ptr = save_multiallelic_hc * 8 + save_hphase * 0x10 + save_dosage * 0x60 + save_dphase * 0x80;
  }
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadRawGenovec(0, vidx, pgrp, &fread_ptr, &fread_end, genovec);
  if ((!(multiallelic_hc_present || save_hphase || save_dosage)) || reterr) {
    *loadbuf_iter_ptr = loadbuf_iter;
    return reterr;
  }

  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  ZeroTrailingQuaters(raw_sample_ct, genovec);
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
  uint32_t het_ct = 0;
  if (multiallelic_hc_present) {
    if (!save_multiallelic_hc) {
      // todo: erase-alt2+ fast path
      // mostly mirror PgrGet2P(0, 1), but a bit of extra logic is needed to
      // throw out phased-10het entries
      return kPglRetNotYetSupported;
    }
    // assume we always save multiallelic info
    // raw format:
    //   rare01_ct, padded out to a word
    //   rare10_ct, padded out to a word
    //   [round up to vector boundary, for patch_01_set]
    //   aux1a, if not mode 15:
    //     patch_01_set as bitarray, raw_sample_ctl words
    //     patch_01_vals, round up to word boundary
    //     [round up to vector boundary, for patch_10_set]
    //   aux1b, if not mode 15:
    //     patch_10_set as bitarray, raw_sample_ctl words
    //     patch_10_vals, round up to word boundary
    // round up to vector boundary at end
    const uint32_t aux1_first_byte = *fread_ptr++;
    const uint32_t aux1a_mode = aux1_first_byte & 15;
    const uint32_t aux1b_mode = aux1_first_byte >> 4;
    uint32_t raw_10_ct = 0;
    if ((!aux1a_mode) || hphase_is_present) {
      if (!aux1b_mode) {
        GenovecCount12Unsafe(genovec, raw_sample_ct, &het_ct, &raw_10_ct);
      } else {
        het_ct = CountQuater(genovec, kMask5555, raw_sample_ct);
      }
    } else if (!aux1b_mode) {
      raw_10_ct = CountQuater(genovec, kMaskAAAA, raw_sample_ct);
    }
    uintptr_t* multihc_raw = loadbuf_iter;
    loadbuf_iter = &(loadbuf_iter[RoundUpPow2(2, kWordsPerVec)]);
    uint32_t rare01_ct = 0;
    if (aux1a_mode != 15) {
      uintptr_t* patch_01_set = loadbuf_iter;
      loadbuf_iter = &(loadbuf_iter[raw_sample_ctl]);
      // (could decide to vector-align patch_01_vals later)
      AlleleCode* patch_01_vals = R_CAST(AlleleCode*, loadbuf_iter);
      reterr = ExportAux1a(fread_end, genovec, aux1a_mode, raw_sample_ct, allele_ct, het_ct, &fread_ptr, patch_01_set, patch_01_vals, &rare01_ct);
      if (unlikely(reterr)) {
        return reterr;
      }
      loadbuf_iter = &(loadbuf_iter[DivUp(rare01_ct, kBytesPerWord / sizeof(AlleleCode))]);
      VecAlignUp64(&loadbuf_iter);
    }
    uint32_t rare10_ct = 0;
    if (aux1b_mode != 15) {
      uintptr_t* patch_10_set = loadbuf_iter;
      loadbuf_iter = &(loadbuf_iter[raw_sample_ctl]);
      AlleleCode* patch_10_vals = R_CAST(AlleleCode*, loadbuf_iter);
      reterr = ExportAux1b(fread_end, genovec, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &fread_ptr, patch_10_set, patch_10_vals, &rare10_ct);
      if (unlikely(reterr)) {
        return reterr;
      }
      loadbuf_iter = &(loadbuf_iter[DivUp(rare10_ct, kBytesPerWord / (2 * sizeof(AlleleCode)))]);
      VecAlignUp64(&loadbuf_iter);
      if (hphase_is_present) {
        het_ct += CountAux1bHets(patch_10_vals, rare10_ct);
      }
    }
    multihc_raw[0] = rare01_ct;
    multihc_raw[1] = rare10_ct;
  } else if (hphase_is_present) {
    het_ct = CountQuater(genovec, kMask5555, raw_sample_ct);
  }

  if (hphase_is_present) {
    if (unlikely(!het_ct)) {
      // there shouldn't be a hphase track at all in this case
      return kPglRetMalformedInput;
    }
    const uint32_t het_ctdl = het_ct / kBitsPerWord;
    uintptr_t* phaseraw = loadbuf_iter;
    const uint32_t first_half_byte_ct = 1 + (het_ct / CHAR_BIT);
    if (save_hphase) {
      // this needs to be synced with MakePgenThread()
#ifdef __LP64__
      // save het_ct later so we can use PopcountWords() below
      phaseraw[0] = 0;
#else
      phaseraw[0] = het_ct;
      phaseraw[1] = 0;
#endif
      loadbuf_iter = &(loadbuf_iter[8 / kBytesPerWord]);
      loadbuf_iter[het_ctdl] = 0;
      memcpy(loadbuf_iter, fread_ptr, first_half_byte_ct);
      loadbuf_iter = &(loadbuf_iter[1 + het_ctdl]);
    }
    const uint32_t explicit_phasepresent = fread_ptr[0] & 1;
    const unsigned char* aux2_start = fread_ptr;
    fread_ptr = &(fread_ptr[first_half_byte_ct]);
    if (explicit_phasepresent) {
      uint32_t raw_phasepresent_ct;
      if (save_hphase) {
#ifdef __LP64__
        raw_phasepresent_ct = PopcountWords(phaseraw, het_ctdl + 2);
#else
        raw_phasepresent_ct = PopcountWords(&(phaseraw[2]), het_ctdl + 1);
#endif
      } else {
        // bugfix (11 Apr 2018): not copied to phaseraw in this case
        raw_phasepresent_ct = PopcountBytes(aux2_start, first_half_byte_ct);
      }
      --raw_phasepresent_ct;
      if (unlikely(!raw_phasepresent_ct)) {
        // there shouldn't be a hphase track at all in this case, either
        return kPglRetMalformedInput;
      }
      const uint32_t second_half_byte_ct = DivUp(raw_phasepresent_ct, CHAR_BIT);
      if (save_hphase) {
#ifdef __LP64__
        phaseraw[0] = het_ct | (S_CAST(uint64_t, raw_phasepresent_ct) << 32);
#else
        phaseraw[1] = raw_phasepresent_ct;
#endif
        memcpy(loadbuf_iter, fread_ptr, second_half_byte_ct);
        loadbuf_iter = &(loadbuf_iter[BitCtToWordCt(raw_phasepresent_ct)]);
      }
      fread_ptr = &(fread_ptr[second_half_byte_ct]);
    }
#ifdef __LP64__
    if (save_hphase) {
      if (!explicit_phasepresent) {
        phaseraw[0] = het_ct;
      }
      VecAlignUp(&loadbuf_iter);
    }
#endif
  }
  if (!save_dosage) {
    *loadbuf_iter_ptr = loadbuf_iter;
    return kPglRetSuccess;
  }
  uintptr_t* dosage_present = loadbuf_iter;
  const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
  loadbuf_iter = &(loadbuf_iter[raw_sample_ctaw]);
  uint16_t* dosage_main = R_CAST(uint16_t*, loadbuf_iter);
  // probable todo: pack this more tightly in the future
  const uintptr_t dosage_main_aligned_wordct = kWordsPerVec * DivUp(raw_sample_ct, (kBytesPerVec / sizeof(int16_t)));
  loadbuf_iter = &(loadbuf_iter[dosage_main_aligned_wordct]);
  uintptr_t* dphase_present = nullptr;
  int16_t* dphase_delta = nullptr;
  if (save_dphase) {
    dphase_present = loadbuf_iter;
    loadbuf_iter = &(loadbuf_iter[raw_sample_ctaw]);
    dphase_delta = R_CAST(int16_t*, loadbuf_iter);
    loadbuf_iter = &(loadbuf_iter[dosage_main_aligned_wordct]);
  }
  *loadbuf_iter_ptr = loadbuf_iter;
  return ParseDosage16(fread_ptr, fread_end, nullptr, raw_sample_ct, vidx, allele_ct, pgrp, nullptr, dphase_present, dphase_delta, nullptr, dosage_present, dosage_main);
}


// Currently assumes no phase or multiallelic hardcalls.
// tried to have more custom code, turned out to not be worth it
PglErr ReadMissingness(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, const unsigned char** fread_endp, uintptr_t* __restrict missingness, uintptr_t* __restrict hets, uintptr_t* __restrict genovec_buf) {
  const unsigned char* fread_ptr;
  const unsigned char* fread_end;
  PglErr reterr = ReadGenovecSubsetUnsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, &fread_ptr, &fread_end, genovec_buf);
  ZeroTrailingQuaters(sample_ct, genovec_buf);
  GenovecToMissingnessUnsafe(genovec_buf, sample_ct, missingness);
  if (hets) {
    PgrDetectGenovecHetsUnsafe(genovec_buf, QuaterCtToWordCt(sample_ct), hets);
  }
  if (fread_pp) {
    *fread_pp = fread_ptr;
    *fread_endp = fread_end;
  }
  return reterr;
}

PglErr PgrGetMissingness(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict missingness, uintptr_t* __restrict genovec_buf) {
  // may as well add a hets parameter?
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  return ReadMissingness(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, nullptr, nullptr, missingness, nullptr, genovec_buf);
}

PglErr PgrGetMissingnessD(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict missingness_hc, uintptr_t* __restrict missingness_dosage, uintptr_t* __restrict hets, uintptr_t* __restrict genovec_buf) {
  // sample_include can't be null
  // either missingness_hc or missingness_dosage must be non-null
  assert(vidx < pgrp->fi.raw_variant_ct);
  if (!sample_ct) {
    return kPglRetSuccess;
  }
  const uint32_t vrtype = GetPgfiVrtype(&(pgrp->fi), vidx);
  const uint32_t dosage_is_relevant = missingness_dosage && VrtypeDosage(vrtype);
  const uint32_t need_to_skip_aux1or2 = dosage_is_relevant && (vrtype & 0x18);
  const uint32_t raw_sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uint32_t subsetting_required = (sample_ct != raw_sample_ct);
  const unsigned char* fread_ptr = nullptr;
  const unsigned char* fread_end = nullptr;
  uintptr_t* missingness_base = missingness_hc? missingness_hc : missingness_dosage;
  if (!need_to_skip_aux1or2) {
    PglErr reterr = ReadMissingness(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, dosage_is_relevant? (&fread_ptr) : nullptr, dosage_is_relevant? (&fread_end) : nullptr, missingness_base, hets, genovec_buf);
    if (missingness_dosage && missingness_hc) {
      memcpy(missingness_dosage, missingness_hc, BitCtToWordCt(sample_ct) * sizeof(intptr_t));
    }
    if (reterr || (!dosage_is_relevant)) {
      return reterr;
    }
  } else {
    PglErr reterr = ReadRawGenovec(subsetting_required, vidx, pgrp, &fread_ptr, &fread_end, genovec_buf);
    if (unlikely(reterr)) {
      return reterr;
    }
    ZeroTrailingQuaters(raw_sample_ct, genovec_buf);
    uintptr_t* subsetted_genovec = pgrp->workspace_vec;
    CopyQuaterarrNonemptySubset(genovec_buf, sample_include, raw_sample_ct, sample_ct, subsetted_genovec);
    GenovecToMissingnessUnsafe(subsetted_genovec, sample_ct, missingness_base);
    if (missingness_hc) {
      memcpy(missingness_dosage, missingness_hc, BitCtToWordCt(sample_ct) * sizeof(intptr_t));
    }

    const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
    const uint32_t allele_ct = allele_idx_offsets? (allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx]) : 2;
    if (VrtypeHphase(vrtype) || hets) {
      uintptr_t* all_hets = pgrp->workspace_all_hets;
      PgrDetectGenovecHets(genovec_buf, raw_sample_ct, all_hets);
      if (VrtypeMultiallelicHc(vrtype)) {
        // see analogous branch in ReadGenovecHphaseSubsetUnsafe()
        // probable todo: make this a separate function
        const uint32_t aux1_first_byte = *fread_ptr++;
        const uint32_t aux1a_mode = aux1_first_byte & 15;
        const uint32_t aux1b_mode = aux1_first_byte >> 4;
        uint32_t raw_01_ct = 0;
        uint32_t raw_10_ct = 0;
        if ((!aux1a_mode) || (!aux1b_mode)) {
          GenovecCount12Unsafe(genovec_buf, raw_sample_ct, &raw_01_ct, &raw_10_ct);
        }
        reterr = SkipAux1a(fread_end, aux1a_mode, raw_sample_ct, allele_ct, raw_01_ct, &fread_ptr);
        if (unlikely(reterr)) {
          return reterr;
        }
        uintptr_t* aux1b_hets = pgrp->workspace_aux1x_present;
        uint32_t* deltalist_workspace = pgrp->workspace_difflist_sample_ids;
        uint32_t aux1b_het_present;
        reterr = GetAux1bHets(fread_end, genovec_buf, aux1b_mode, raw_sample_ct, allele_ct, raw_10_ct, &fread_ptr, aux1b_hets, &aux1b_het_present, deltalist_workspace);
        if (unlikely(reterr)) {
          return reterr;
        }
        if (aux1b_het_present) {
          BitvecOr(aux1b_hets, raw_sample_ctl, all_hets);
        }
      }
      if (hets) {
        CopyBitarrSubset(all_hets, sample_include, sample_ct, hets);
      }
      if (VrtypeHphase(vrtype)) {
        reterr = SkipAux2(fread_end, PopcountWords(all_hets, raw_sample_ctl), &fread_ptr, nullptr);
        if (unlikely(reterr)) {
          return reterr;
        }
      }
    } else {
      SkipAux1(fread_end, genovec_buf, raw_sample_ct, allele_ct, &fread_ptr);
    }
  }
  // now perform bitwise andnot with dosage_present
  if ((vrtype & 0x60) == 0x40) {
    // unconditional dosage.  spot-check the appropriate entries for equality
    // to 65535.
#ifdef __arm__
#  error "Unaligned accesses in PgrGetMissingnessPD()."
#endif
    const uint16_t* dosage_main = R_CAST(const uint16_t*, fread_ptr);
    // bugfix (18 Feb 2019): sample_include is permitted to be nullptr here
    if (!subsetting_required) {
      // probable todo: faster iteration over set bits
      for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
        uintptr_t missing_dosage_bits = missingness_dosage[widx];
        if (missing_dosage_bits) {
          const uint16_t* cur_dosage_main = &(dosage_main[widx * kBitsPerWord]);
          do {
            uint32_t sample_idx_lowbits = ctzw(missing_dosage_bits);
            if (cur_dosage_main[sample_idx_lowbits] != 65535) {
              missingness_dosage[widx] ^= missing_dosage_bits & (-missing_dosage_bits);
            }
            missing_dosage_bits &= missing_dosage_bits - 1;
          } while (missing_dosage_bits);
        }
      }
    } else {
      uintptr_t sample_uidx_base = 0;
      uintptr_t sample_include_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
        if (!IsSet(missingness_dosage, sample_idx)) {
          continue;
        }
        if (dosage_main[sample_uidx] != 65535) {
          ClearBit(sample_idx, missingness_dosage);
        }
      }
    }
    return kPglRetSuccess;
  }
  uintptr_t* dosage_present = pgrp->workspace_dosage_present;
  if ((vrtype & 0x60) == 0x20) {
    // dosage list
    uint32_t dummy;
    if (unlikely(ParseAndSaveDeltalistAsBitarr(fread_end, raw_sample_ct, &fread_ptr, dosage_present, &dummy))) {
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
  BitvecInvmask(dosage_present, BitCtToWordCt(sample_ct), missingness_dosage);
  return kPglRetSuccess;
}

static inline BoolErr ValidateVint31(const unsigned char* buf_end, const unsigned char** bufpp, uint32_t* val_ptr) {
  if (unlikely(buf_end <= (*bufpp))) {
    return 1;
  }
  uint32_t vint32 = *((*bufpp)++);
  if (vint32 <= 127) {
    *val_ptr = vint32;
    return 0;
  }
  vint32 &= 127;
  for (uint32_t shift = 7; shift != 28; shift += 7) {
    if (unlikely(buf_end == (*bufpp))) {
      return 1;
    }
    uint32_t uii = *((*bufpp)++);
    vint32 |= (uii & 127) << shift;
    if (uii <= 127) {
      *val_ptr = vint32;
      return 0;
    }
  }
  if (unlikely(buf_end == (*bufpp))) {
    return 1;
  }
  uint32_t uii = *((*bufpp)++);
  if (unlikely(uii > 7)) {
    return 1;
  }
  vint32 |= uii << 28;
  *val_ptr = vint32;
  return 0;
}

BoolErr ValidateDifflistHeader(const unsigned char* fread_end, uint32_t sample_ct, const unsigned char** fread_pp, uintptr_t* raregeno_buf, const unsigned char** difflist_group_info_ptr, uint32_t* difflist_len_ptr) {
  // can be used for deltalists: pass raregeno_buf == nullptr.
  if (unlikely(ValidateVint31(fread_end, fread_pp, difflist_len_ptr))) {
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
  if (unlikely(difflist_len > sample_ct / kPglMaxDifflistLenDivisor)) {
    return 1;
  }
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(sample_ct);
  const uint32_t difflist_index_byte_ct = group_ct * (sample_id_byte_ct + 1) - 1;
  if (PtrAddCk(fread_end, difflist_index_byte_ct, fread_pp)) {
    return 1;
  }
  if (!raregeno_buf) {
    return 0;
  }
  const uint32_t raregeno_byte_ct = QuaterCtToByteCt(difflist_len);
  const unsigned char* raregeno_start = *fread_pp;
  if (PtrAddCk(fread_end, raregeno_byte_ct, fread_pp)) {
    return 1;
  }
  memcpy(raregeno_buf, raregeno_start, raregeno_byte_ct);
  const uint32_t difflist_len_mod4 = difflist_len % 4;
  if (difflist_len_mod4) {
    const uint32_t last_raregeno_byte = (*fread_pp)[-1];
    if (unlikely(last_raregeno_byte >> (2 * difflist_len_mod4))) {
      return 1;
    }
  }
  return 0;
}

BoolErr ValidateAndApplyDifflist(const unsigned char* fread_end, uint32_t common2_code, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // Side effects: uses pgr.workspace_raregeno_tmp_loadbuf.
  // Similar to ParseAndApplyDifflist(), but with exhaustive input
  // validation.
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  uintptr_t* cur_raregeno_iter = pgrp->workspace_raregeno_tmp_loadbuf;
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  if (unlikely(ValidateDifflistHeader(fread_end, sample_ct, fread_pp, cur_raregeno_iter, &group_info_iter, &difflist_len))) {
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
    if (unlikely(difflist_len >= (sample_ct / (2 * kPglMaxDifflistLenDivisor)))) {
      return 1;
    }
    const uintptr_t common_code_delta = common2_code & 3;
    const uintptr_t inv_common_word1 = (3 - common2_code / 4) * kMask5555;
    const uintptr_t inv_common_word2 = inv_common_word1 - (common_code_delta * kMask5555);
    for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
      uintptr_t cur_raregeno_word = cur_raregeno_iter[subgroup_idx];
      const uintptr_t match1 = Word11(cur_raregeno_word ^ inv_common_word1);
      const uintptr_t match2 = Word11(cur_raregeno_word ^ inv_common_word2);
      if (subgroup_idx == subgroup_idx_last) {
        // ignore trailing bits
        const uint32_t lshift = ((-difflist_len) % kBitsPerWordD2) * 2;
        if (unlikely((match1 << lshift) || (match2 << lshift))) {
          return 1;
        }
        break;
      }
      if (unlikely(match1 || match2)) {
        // todo: if (multiallelic_hc_present && (!inv_common_word2)), record
        // might be fine; but we need to verify these are actually rare alt
        // alleles.
        // (er, above comment is obsolete)
        return 1;
      }
    }
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(sample_ct);
  const unsigned char* group_byte_cts_iter = &(group_info_iter[DivUp(difflist_len, kPglDifflistGroupSize) * sample_id_byte_ct]);
  const unsigned char* prev_group_start = *fread_pp;

  uintptr_t sample_idx = 0;
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
    uint32_t remaining_deltas_in_subgroup = kBitsPerWordD2 - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        return 0;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    if (!(subgroup_idx % (kPglDifflistGroupSize / kBitsPerWordD2))) {
      uintptr_t new_sample_idx_start = SubU32Load(group_info_iter, sample_id_byte_ct);
      if (subgroup_idx) {
        if (unlikely(sample_idx >= new_sample_idx_start)) {
          return 1;
        }
        const uint32_t group_byte_ct = S_CAST(uint32_t, *group_byte_cts_iter++) + 63;
        if (unlikely(S_CAST(uintptr_t, (*fread_pp) - prev_group_start) != group_byte_ct)) {
          return 1;
        }
        prev_group_start = *fread_pp;
      }
      sample_idx = new_sample_idx_start;
      group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    } else {
      uint32_t sample_idx_incr;
      if (unlikely(ValidateVint31(fread_end, fread_pp, &sample_idx_incr) || (!sample_idx_incr))) {
        return 1;
      }
      sample_idx += sample_idx_incr;
    }
    uintptr_t cur_raregeno_word = *cur_raregeno_iter++;
    for (; ; --remaining_deltas_in_subgroup) {
      if (unlikely(sample_idx >= sample_ct)) {
        return 1;
      }
      const uintptr_t cur_geno = cur_raregeno_word & 3;
      AssignQuaterarrEntry(sample_idx, cur_geno, genovec);
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      uint32_t sample_idx_incr;
      if (unlikely(ValidateVint31(fread_end, fread_pp, &sample_idx_incr) || (!sample_idx_incr))) {
        return 1;
      }
      sample_idx += sample_idx_incr;
      cur_raregeno_word >>= 2;
    }
  }
}

BoolErr ValidateOnebit(const unsigned char* fread_end, const unsigned char** fread_pp, PgenReader* pgrp, uintptr_t* __restrict genovec) {
  // ParseOnebitUnsafe() with exhaustive input validation.
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t common2_and_bitarray_byte_ct = (sample_ct + 15) / CHAR_BIT;
  const unsigned char* onebit_main_iter = *fread_pp;
  if (PtrAddCk(fread_end, common2_and_bitarray_byte_ct, fread_pp)) {
    return 1;
  }
  const uintptr_t common2_code = *onebit_main_iter++;
  const uintptr_t common_code_delta = common2_code & 3;
  uintptr_t word_base = common2_code / 4;
  if (unlikely((!common_code_delta) || (word_base + common_code_delta > 3))) {
    return 1;
  }
  word_base *= kMask5555;
  const uint32_t genovec_widx_trail = (sample_ct + 7) / kBitsPerWordD2;
  const uint32_t genovec_widx_end = QuaterCtToWordCt(sample_ct);
#ifdef __arm__
#  error "Unaligned accesses in ValidateOnebit()."
#endif
  const Halfword* onebit_main = R_CAST(const Halfword*, onebit_main_iter);
  for (uint32_t genovec_widx = 0; ; ++genovec_widx) {
    uintptr_t ww;
    if (genovec_widx >= genovec_widx_trail) {
      if (genovec_widx == genovec_widx_end) {
        break;
      }
      const uint32_t nontrail_byte_ct = ((sample_ct - 1) % kBitsPerWordD2) / CHAR_BIT;
      ww = ProperSubwordLoad(&(onebit_main[genovec_widx_trail]), 1 + nontrail_byte_ct);
      const uint32_t sample_ct_mod8 = sample_ct % 8;
      if (sample_ct_mod8) {
        if (unlikely(ww >> (nontrail_byte_ct * 8 + sample_ct_mod8))) {
          return 1;
        }
      }
    } else {
      ww = onebit_main[genovec_widx];
    }
    ww = UnpackHalfwordToWord(ww);
    genovec[genovec_widx] = word_base + ww * common_code_delta;
  }
  return ValidateAndApplyDifflist(fread_end, common2_code, fread_pp, pgrp, genovec);
}

// assumes that we aren't dealing with the trivial fixed-width case.
// saves main genotype array to genovec.  does not zero out trailing bits.
BoolErr ValidateGeno(const unsigned char* fread_end, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, uintptr_t* genovec, char* errstr_buf) {
  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  if (VrtypeLdCompressed(vrtype)) {
    CopyQuaterarr(pgrp->ldbase_genovec, sample_ct, genovec);
    if (unlikely(ValidateAndApplyDifflist(fread_end, 0, fread_pp, pgrp, genovec))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid LD difflist for (0-based) variant #%u.\n", vidx);
      return 1;
    }
    if (vrtype & 1) {
      GenovecInvertUnsafe(sample_ct, genovec);
    }
    return 0;
  }
  const uint32_t is_ldbase = VrtypeLdCompressed(pgrp->fi.vrtypes[vidx + 1]);
  if (!(vrtype & 4)) {
    if (vrtype & 1) {
      if (unlikely(ValidateOnebit(fread_end, fread_pp, pgrp, genovec))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid 1-bit genotype record for (0-based) variant #%u.\n", vidx);
        return 1;
      }
    } else {
      const uint32_t genovec_byte_ct = DivUp(sample_ct, 4);
      const unsigned char* src_genodata = *fread_pp;
      if (PtrAddCk(fread_end, genovec_byte_ct, fread_pp)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid 2-bit genotype record for (0-based) variant #%u\n", vidx);
        return 1;
      }
      memcpy(genovec, src_genodata, genovec_byte_ct);
      const uint32_t sample_ct_mod4 = sample_ct % 4;
      if (sample_ct_mod4) {
        const uint32_t last_geno_byte = (*fread_pp)[-1];
        if (unlikely(last_geno_byte >> (2 * sample_ct_mod4))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Last genotype byte for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
          return 1;
        }
      }
    }
  } else {
    const uint32_t vrtype_low2 = vrtype & 3;
    if (vrtype_low2 != 1) {
      const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
      vecset(genovec, vrtype_low2 * kMask5555, vec_ct);
      if (unlikely(ValidateAndApplyDifflist(fread_end, 0, fread_pp, pgrp, genovec))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid genotype difflist for (0-based) variant #%u.\n", vidx);
        return 1;
      }
    } else {
      if (is_ldbase) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid LD back-reference from variant #%u to all-hom-ref variant #%u.\n", vidx + 1, vidx);
        return 1;
      }
      ZeroWArr(QuaterCtToWordCt(sample_ct), genovec);
    }
  }
  if (is_ldbase) {
    CopyQuaterarr(genovec, sample_ct, pgrp->ldbase_genovec);
  }
  return 0;
}

BoolErr ValidateAndCountDeltalist(const unsigned char* fread_end, uint32_t sample_ct, const unsigned char** fread_pp, uint32_t* __restrict deltalist, uint32_t* deltalist_len_ptr) {
  // pass deltalist == nullptr when actual bit positions aren't needed
  const unsigned char* group_info_iter;
  if (unlikely(ValidateDifflistHeader(fread_end, sample_ct, fread_pp, nullptr, &group_info_iter, deltalist_len_ptr))) {
    return 1;
  }
  const uint32_t deltalist_len = *deltalist_len_ptr;
  if (!deltalist_len) {
    return 0;
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(sample_ct);
  const uint32_t group_idx_last = (deltalist_len - 1) / kPglDifflistGroupSize;
  const unsigned char* group_byte_cts_iter = &(group_info_iter[DivUp(deltalist_len, kPglDifflistGroupSize) * sample_id_byte_ct]);
  const unsigned char* prev_group_start = *fread_pp;
  uint32_t* deltalist_iter = deltalist;
  uint32_t group_len_m1 = kPglDifflistGroupSize - 1;
  uintptr_t sample_idx = 0;
  for (uint32_t group_idx = 0; ; ++group_idx) {
    if (group_idx >= group_idx_last) {
      if (group_idx > group_idx_last) {
        return 0;
      }
      group_len_m1 &= deltalist_len - 1;
    }
    uintptr_t new_sample_idx = SubU32Load(group_info_iter, sample_id_byte_ct);
    if (group_idx) {
      if (unlikely(sample_idx >= new_sample_idx)) {
        return 1;
      }
      const uint32_t group_byte_ct = S_CAST(uint32_t, *group_byte_cts_iter++) + 63;
      if (unlikely(S_CAST(uintptr_t, (*fread_pp) - prev_group_start) != group_byte_ct)) {
        return 1;
      }
      prev_group_start = *fread_pp;
    }
    sample_idx = new_sample_idx;
    group_info_iter = &(group_info_iter[sample_id_byte_ct]);
    for (uint32_t deltalist_idx_lowbits = 0; ; ++deltalist_idx_lowbits) {
      if (unlikely(sample_idx >= sample_ct)) {
        return 1;
      }
      if (deltalist_iter) {
        *deltalist_iter++ = sample_idx;
      }
      if (deltalist_idx_lowbits == group_len_m1) {
        break;
      }
      uint32_t sample_idx_incr;
      if (unlikely(ValidateVint31(fread_end, fread_pp, &sample_idx_incr) || (!sample_idx_incr))) {
        return 1;
      }
      sample_idx += sample_idx_incr;
    }
  }
}

BoolErr ValidateMultiallelicHc(const unsigned char* fread_end, const uintptr_t* __restrict raw_genovec, uint32_t vidx, uint32_t allele_ct, PgenReader* pgrp, const unsigned char** fread_pp, uint32_t* __restrict het_ctp, char* __restrict errstr_buf) {
  if (unlikely(allele_ct <= 2)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Multiallelic hardcall track present for (0-based) variant #%u, but it apparently has only %u allele%s.\n", vidx, allele_ct, (allele_ct == 1)? "" : "s");
    return 1;
  }
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t aux1_first_byte = **fread_pp;
  *fread_pp += 1;
  if (unlikely(
          aux1_first_byte &&
          (aux1_first_byte != 1) &&
          (aux1_first_byte != 15) &&
          (aux1_first_byte != 16) &&
          (aux1_first_byte != 17) &&
          (aux1_first_byte != 31) &&
          (aux1_first_byte != 240) &&
          (aux1_first_byte != 241))) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic hardcall track mode byte (%u; must be in {0, 1, 15, 16, 17, 31, 240, 241}) in (0-based) variant #%u.\n", aux1_first_byte, vidx);
    return 1;
  }
  const uint32_t aux1a_mode = aux1_first_byte & 15;
  const uint32_t aux1b_mode = aux1_first_byte >> 4;
  uint32_t raw_01_ct;
  uint32_t raw_10_ct;
  GenovecCount12Unsafe(raw_genovec, sample_ct, &raw_01_ct, &raw_10_ct);
  uint32_t* deltalist_workspace = pgrp->workspace_difflist_sample_ids;
  if (aux1a_mode != 15) {
    if (unlikely(!raw_01_ct)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Multiallelic het-ref hardcall track present for (0-based) variant #%u, but no het-ref calls exist.\n", vidx);
      return 1;
    }
    uint32_t rare01_ct;
    if (!aux1a_mode) {
      const uint32_t subset_byte_ct = DivUp(raw_01_ct, CHAR_BIT);
      if (PtrCheck(fread_end, *fread_pp, subset_byte_ct)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall bitarray-subset for (0-based) variant #%u.\n", vidx);
        return 1;
      }
      rare01_ct = PopcountBytes(*fread_pp, subset_byte_ct);
      if (unlikely(!rare01_ct)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Empty multiallelic het-ref hardcall bitarray-subset for (0-based) variant #%u.\n", vidx);
        return 1;
      }
      *fread_pp += subset_byte_ct;
      const uint32_t raw_01_ct_mod8 = raw_01_ct % 8;
      if (raw_01_ct_mod8) {
        if (unlikely((*fread_pp)[-1] >> raw_01_ct_mod8)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Multiallelic het-ref hardcall bitarray-subset for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
          return 1;
        }
      }
    } else {
      if (unlikely(ValidateAndCountDeltalist(fread_end, sample_ct, fread_pp, deltalist_workspace, &rare01_ct) || (!rare01_ct))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall deltalist-subset for (0-based) variant #%u.\n", vidx);
        return 1;
      }
      for (uint32_t uii = 0; uii != rare01_ct; ++uii) {
        if (unlikely(GetQuaterarrEntry(raw_genovec, deltalist_workspace[uii]) != 1)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall deltalist-subset for (0-based) variant #%u (an index doesn't correspond to a het-ref call).\n", vidx);
          return 1;
        }
      }
    }
    if (allele_ct < 5) {
      // Nothing to do for allele_ct == 3.
      if (allele_ct == 4) {
        // 1-bit entries.  Contents must be in range, so just validate trailing
        // bits.
        const uint32_t fvals_byte_ct = DivUp(rare01_ct, 8);
        if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (shorter than expected).\n", vidx);
          return 1;
        }
        const uint32_t rare01_ct_mod8 = rare01_ct % 8;
        if (rare01_ct_mod8) {
          if (unlikely((*fread_pp)[-1] >> rare01_ct_mod8)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (nonzero trailing bits).\n", vidx);
            return 1;
          }
        }
      }
    } else {
      const unsigned char* fvals = *fread_pp;
      if (allele_ct < 19) {
        if (allele_ct < 7) {
          // 2-bit entries.
          const uint32_t fvals_byte_ct = DivUp(rare01_ct, 4);
          if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (shorter than expected).\n", vidx);
            return 1;
          }
          if (allele_ct == 5) {
            // Contents may be out-of-range.
            const uint32_t fullword_ct = fvals_byte_ct / kBytesPerWord;
            uint32_t widx = 0;
            if (fullword_ct) {
              const uintptr_t* fvals_alias = R_CAST(const uintptr_t*, fvals);
              for (; widx != fullword_ct; ++widx) {
                const uintptr_t cur_word = fvals_alias[widx];
                if (unlikely(cur_word & (cur_word >> 1) & kMask5555)) {
                  snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (out-of-range allele code).\n", vidx);
                  return 1;
                }
              }
            }
            for (uint32_t uii = widx * kBytesPerWord; uii != fvals_byte_ct; ++uii) {
              const uint32_t cur_byte = fvals[uii];
              if (unlikely(cur_byte & (cur_byte >> 1) & 0x55)) {
                snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (out-of-range allele code).\n", vidx);
                return 1;
              }
            }
          }
          // Validate trailing bits.
          const uint32_t rare01_ct_mod4 = rare01_ct % 4;
          if (rare01_ct_mod4) {
            if (unlikely((*fread_pp)[-1] >> (2 * rare01_ct_mod4))) {
              snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (nonzero trailing bits).\n", vidx);
              return 1;
            }
          }
        } else {
          // 4-bit entries.
          const uint32_t fvals_byte_ct = DivUp(rare01_ct, 2);
          if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (shorter than expected).\n", vidx);
            return 1;
          }
          if (allele_ct != 18) {
            // Contents may be out-of-range.
            // (Can optimize this loop later.)
            const uint32_t max_code = allele_ct - 3;
            for (uint32_t uii = 0; uii != fvals_byte_ct; ++uii) {
              const uint32_t cur_byte = fvals[uii];
              if (unlikely(((cur_byte & 15) > max_code) || ((cur_byte >> 4) > max_code))) {
                snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (out-of-range allele code).\n", vidx);
                return 1;
              }
            }
          }
          // Validate trailing bits.
          const uint32_t rare01_ct_mod2 = rare01_ct % 2;
          if (rare01_ct_mod2) {
            if (unlikely((*fread_pp)[-1] >> 4)) {
              snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (nonzero trailing bits).\n", vidx);
              return 1;
            }
          }
        }
      } else {
        // 8-bit entries.
        if (PtrAddCk(fread_end, rare01_ct, fread_pp)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (shorter than expected).\n", vidx);
          return 1;
        }
        // Can optimize this loop later.
        const uint32_t max_code = allele_ct - 3;
        for (uint32_t uii = 0; uii != rare01_ct; ++uii) {
          const uint32_t cur_byte = fvals[uii];
          if (unlikely(cur_byte > max_code)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic het-ref hardcall track for (0-based) variant #%u (out-of-range allele code).\n", vidx);
            return 1;
          }
        }
      }
    }
  }
  if (aux1b_mode != 15) {
    if (unlikely(!raw_10_ct)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Multiallelic altxy hardcall track present for (0-based) variant #%u, but no altxy calls exist.\n", vidx);
      return 1;
    }
    uint32_t rare10_ct;
    if (!aux1b_mode) {
      const uint32_t subset_byte_ct = DivUp(raw_10_ct, CHAR_BIT);
      if (PtrCheck(fread_end, *fread_pp, subset_byte_ct)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall bitarray-subset for (0-based) variant #%u.\n", vidx);
        return 1;
      }
      rare10_ct = PopcountBytes(*fread_pp, subset_byte_ct);
      if (unlikely(!rare10_ct)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Empty multiallelic altxy hardcall bitarray-subset for (0-based) variant #%u.\n", vidx);
        return 1;
      }
      *fread_pp += subset_byte_ct;
      const uint32_t raw_10_ct_mod8 = raw_10_ct % 8;
      if (raw_10_ct_mod8) {
        if (unlikely((*fread_pp)[-1] >> raw_10_ct_mod8)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Multiallelic altxy hardcall bitarray-subset for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
          return 1;
        }
      }
    } else {
      if (unlikely(ValidateAndCountDeltalist(fread_end, sample_ct, fread_pp, deltalist_workspace, &rare10_ct) || (!rare10_ct))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall deltalist-subset for (0-based) variant #%u.\n", vidx);
        return 1;
      }
      for (uint32_t uii = 0; uii != rare10_ct; ++uii) {
        if (unlikely(GetQuaterarrEntry(raw_genovec, deltalist_workspace[uii]) != 2)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall deltalist-subset for (0-based) variant #%u (an index doesn't correspond to an altxy call).\n", vidx);
          return 1;
        }
      }
    }
    const unsigned char* fvals = *fread_pp;
    uint32_t het_incr;
    if (allele_ct < 6) {
      if (allele_ct == 3) {
        // 1-bit entries.
        const uint32_t fvals_byte_ct = DivUp(rare10_ct, 8);
        if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (shorter than expected).\n", vidx);
          return 1;
        }
        const uint32_t rare10_ct_mod8 = rare10_ct % 8;
        if (rare10_ct_mod8) {
          if (unlikely((*fread_pp)[-1] >> rare10_ct_mod8)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (nonzero trailing bits).\n", vidx);
            return 1;
          }
        }
        het_incr = rare10_ct - PopcountBytes(fvals, fvals_byte_ct);
      } else {
        // 2+2 bit entries.
        const uint32_t fvals_byte_ct = DivUp(rare10_ct, 2);
        if (PtrAddCk(fread_end, fvals_byte_ct, fread_pp)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (shorter than expected).\n", vidx);
          return 1;
        }
        // Can optimize this later.
        uint64_t nibble_cts[16];
        ZeroU64Arr(16, nibble_cts);
        CountAllNibbles64(fvals, rare10_ct, nibble_cts);
        // 1/1 is invalid here
        if (nibble_cts[0]) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (out-of-range allele code pair).\n", vidx);
          return 1;
        }
        const uint32_t max_code = allele_ct - 2;
        for (uint32_t hi_code = 0; hi_code != 4; ++hi_code) {
          uint32_t lo_code = hi_code + 1;
          if (hi_code > max_code) {
            lo_code = 0;
          }
          const uint64_t* nibble_cts_offset = &(nibble_cts[hi_code * 4]);
          for (; lo_code != 4; ++lo_code) {
            if (nibble_cts_offset[lo_code]) {
              snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (out-of-range allele code pair).\n", vidx);
              return 1;
            }
          }
        }
        const uintptr_t rarehom_ct = nibble_cts[5] + nibble_cts[10] + nibble_cts[15];
        het_incr = rare10_ct - rarehom_ct;
        const uint32_t rare10_ct_mod2 = rare10_ct % 2;
        if (rare10_ct_mod2) {
          if (unlikely((*fread_pp)[-1] >> 4)) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (nonzero trailing bits).\n", vidx);
            return 1;
          }
        }
      }
    } else {
      if (allele_ct < 18) {
        // 4+4 bit entries.
        if (PtrAddCk(fread_end, rare10_ct, fread_pp)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (shorter than expected).\n", vidx);
          return 1;
        }
        const uint32_t max_code = allele_ct - 2;
        het_incr = 0;
        for (uint32_t uii = 0; uii != rare10_ct; ++uii) {
          const uint32_t cur_byte = fvals[uii];
          const uint32_t lo_code = cur_byte & 15;
          const uint32_t hi_code = cur_byte >> 4;
          if (unlikely((!hi_code) || (hi_code > max_code) || (lo_code > hi_code))) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (out-of-range or misordered allele code pair).\n", vidx);
            return 1;
          }
          het_incr += (lo_code != hi_code);
        }
      } else {
        // 8+8 bit entries
        if (PtrAddCk(fread_end, 2 * rare10_ct, fread_pp)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (shorter than expected).\n", vidx);
          return 1;
        }
        const uint32_t max_code = allele_ct - 2;
        het_incr = 0;
        for (uint32_t uii = 0; uii != rare10_ct; ++uii) {
          const AlleleCode lo_code = fvals[2 * uii];
          const AlleleCode hi_code = fvals[2 * uii + 1];
          if (unlikely((!hi_code) || (hi_code > max_code) || (lo_code > hi_code))) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid multiallelic altxy hardcall track for (0-based) variant #%u (out-of-range or misordered allele code pair).\n", vidx);
            return 1;
          }
          het_incr += (lo_code != hi_code);
        }
      }
    }
    *het_ctp += het_incr;
  }
  return 0;
}

BoolErr ValidateHphase(const unsigned char* fread_end, uint32_t vidx, uint32_t het_ct, const unsigned char** fread_pp, char* errstr_buf) {
  if (unlikely(!het_ct)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Hardcall phase track present for (0-based) variant #%u, but there were no heterozygous calls.\n", vidx);
    return 1;
  }
  const uint32_t aux2_first_part_byte_ct = 1 + (het_ct / CHAR_BIT);
  const unsigned char* aux2_first_part = *fread_pp;
  if (PtrAddCk(fread_end, aux2_first_part_byte_ct, fread_pp)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid hardcall phase track present for (0-based) variant #%u.\n", vidx);
    return 1;
  }
  const uint32_t het_ct_p1_mod8 = (het_ct + 1) % CHAR_BIT;
  if (het_ct_p1_mod8) {
    // verify trailing bits are zero
    if (unlikely((*fread_pp)[-1] >> het_ct_p1_mod8)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Hardcall phase track for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
      return 1;
    }
  }
  if (!((*aux2_first_part) & 1)) {
    // phase always present, "first part" is only part
    return 0;
  }
  const uint32_t phasepresent_ct = PopcountBytes(aux2_first_part, aux2_first_part_byte_ct) - 1;
  if (unlikely(!phasepresent_ct)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Hardcall phase track for (0-based) variant #%u does not have any actual phase information.\n", vidx);
    return 1;
  }
  const uint32_t phaseinfo_byte_ct = DivUp(phasepresent_ct, CHAR_BIT);
  if (PtrAddCk(fread_end, phaseinfo_byte_ct, fread_pp)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid hardcall phase track present for (0-based) variant #%u.\n", vidx);
    return 1;
  }
  const uint32_t phasepresent_ct_mod8 = phasepresent_ct % 8;
  if (phasepresent_ct_mod8) {
    if (unlikely((*fread_pp)[-1] >> phasepresent_ct_mod8)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Hardcall phase track for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
      return 1;
    }
  }
  return 0;
}

PglErr ValidateDosage16(const unsigned char* fread_end, uint32_t vidx, PgenReader* pgrp, const unsigned char** fread_pp, char* errstr_buf) {
  // similar to ParseDosage16().  doesn't support multiallelic data yet.
  const uint32_t vrtype = pgrp->fi.vrtypes[vidx];
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  if ((vrtype & 0x60) == 0x40) {
    // unconditional dosage.  handle separately from other two cases since
    // 65535 is valid.
#ifdef __arm__
#  error "Unaligned accesses in ValidateDosage16()."
#endif
    const uint16_t* dosage_main = R_CAST(const uint16_t*, *fread_pp);
    if (PtrAddCk(fread_end, sample_ct * sizeof(int16_t), fread_pp)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid unconditional dosage track for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
    // todo: verify genotype and dosage are consistent
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      uint16_t cur_dosage_val_p1 = dosage_main[sample_idx];
      cur_dosage_val_p1 += 1;  // intentional overflow on 65535
      if (unlikely(cur_dosage_val_p1 > 32769)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid unconditional dosage track for (0-based) variant #%u (dosage is greater than 2).\n", vidx);
        return kPglRetMalformedInput;
      }
    }
    if (vrtype & 0x80) {
      const int16_t* dphase_delta = R_CAST(const int16_t*, *fread_pp);
      if (PtrAddCk(fread_end, sample_ct * sizeof(int16_t), fread_pp)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid unconditional phased-dosages for (0-based) variant #%u.\n", vidx);
        return kPglRetMalformedInput;
      }
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uint16_t dosage_val = dosage_main[sample_idx];
        const int16_t dphase_delta_val = dphase_delta[sample_idx];
        const uint16_t dpiece0_x2 = dosage_val + dphase_delta_val;
        const uint16_t dpiece1_x2 = dosage_val - dphase_delta_val;
        // Update (11 May 2018): parity condition removed.
        if ((dpiece0_x2 > 32768) || (dpiece1_x2 > 32768)) {
          if (unlikely((dphase_delta_val != -32768) || (dosage_val != 65535))) {
            snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid unconditional phased-dosages for (0-based) variant #%u.\n", vidx);
            return kPglRetMalformedInput;
          }
        }
      }
    }
    return kPglRetSuccess;
  }
  uint32_t dosage_ct;
  if ((vrtype & 0x60) == 0x20) {
    // dosage list
    if (unlikely(ValidateAndCountDeltalist(fread_end, sample_ct, fread_pp, nullptr, &dosage_ct))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid dosage list for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
  } else {
    const uint32_t sample_ctb = DivUp(sample_ct, CHAR_BIT);
    if (PtrCheck(fread_end, *fread_pp, sample_ctb)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid dosage subset for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
    dosage_ct = PopcountBytes(*fread_pp, sample_ctb);
    *fread_pp += sample_ctb;
    const uint32_t sample_ct_mod8 = sample_ct % 8;
    if (sample_ct_mod8) {
      if (unlikely((*fread_pp)[-1] >> sample_ct_mod8)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Dosage subset bitarray for (0-based) variant #%u has nonzero trailing bits.\n", vidx);
        return kPglRetMalformedInput;
      }
    }
  }
  const uint16_t* dosage_main = R_CAST(const uint16_t*, *fread_pp);
  if (PtrAddCk(fread_end, dosage_ct * sizeof(int16_t), fread_pp)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid dosage track for (0-based) variant #%u.\n", vidx);
    return kPglRetMalformedInput;
  }
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    if (unlikely(dosage_main[dosage_idx] > 32768)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid dosage track for (0-based) variant #%u (dosage is greater than 2).\n", vidx);
      return kPglRetMalformedInput;
    }
  }
  if (vrtype & 0x80) {
    const uintptr_t* file_dphase_present = R_CAST(const uintptr_t*, *fread_pp);
    const uint32_t dphase_present_byte_ct = DivUp(dosage_ct, CHAR_BIT);
    if (PtrAddCk(fread_end, dphase_present_byte_ct, fread_pp)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid phased-dosage track for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
    const uint32_t trailing_bit_ct = dosage_ct % CHAR_BIT;
    if (unlikely(trailing_bit_ct && ((*fread_pp)[-1] & (255 << trailing_bit_ct)))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid phased-dosage track for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
    const uint16_t* dosage_main_read_iter = dosage_main;
    const int16_t* dphase_delta_read_iter = R_CAST(const int16_t*, *fread_pp);
    const uint32_t dphase_widx_last = (dphase_present_byte_ct - 1) / kBytesPerWord;
    uint32_t loop_end = kBitsPerWord;
    for (uint32_t dphase_widx = 0; ; ++dphase_widx) {
      uintptr_t ww;
      if (dphase_widx >= dphase_widx_last) {
        if (dphase_widx > dphase_widx_last) {
          break;
        }
        loop_end = 1 + ((dosage_ct - 1) % kBitsPerWord);
        const uint32_t final_byte_ct = DivUp(loop_end, CHAR_BIT);
        ww = SubwordLoad(&(file_dphase_present[dphase_widx]), final_byte_ct);
      } else {
        ww = file_dphase_present[dphase_widx];
      }
      for (uint32_t dphase_lowbits = 0; dphase_lowbits != loop_end; ++dphase_lowbits, ++dosage_main_read_iter) {
        if (!((ww >> dphase_lowbits) & 1)) {
          continue;
        }
        const uint16_t dosage_val = *dosage_main_read_iter;
        const int16_t dphase_delta_val = *dphase_delta_read_iter++;
        const uint16_t dpiece0_x2 = dosage_val + dphase_delta_val;
        const uint16_t dpiece1_x2 = dosage_val - dphase_delta_val;
        if (unlikely((dpiece0_x2 > 32768) || (dpiece1_x2 > 32768))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid phased-dosage track for (0-based) variant #%u.\n", vidx);
          return kPglRetMalformedInput;
        }
      }
    }
    if (unlikely(dphase_delta_read_iter == R_CAST(const int16_t*, *fread_pp))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid phased-dosage track for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
    *fread_pp = R_CAST(const unsigned char*, dphase_delta_read_iter);
    if (PtrCheck(fread_end, *fread_pp, 0)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid phased-dosage track for (0-based) variant #%u.\n", vidx);
      return kPglRetMalformedInput;
    }
  }
  return kPglRetSuccess;
}

static_assert(kPglVblockSize == 65536, "PgrValidate() needs to have an error message updated.");
PglErr PgrValidate(PgenReader* pgrp, uintptr_t* genovec_buf, char* errstr_buf) {
  // Performs all validation which isn't done by pgfi_init_phase{1,2}() and
  // PgrInit().
  const uintptr_t* allele_idx_offsets = pgrp->fi.allele_idx_offsets;
  const uint32_t variant_ct = pgrp->fi.raw_variant_ct;
  const uint32_t sample_ct = pgrp->fi.raw_sample_ct;
  const uint32_t const_vrtype = pgrp->fi.const_vrtype;
  if (const_vrtype != UINT32_MAX) {
    if (unlikely(allele_idx_offsets && (allele_idx_offsets[variant_ct] != 2 * variant_ct))) {
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
      for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
        const unsigned char* fread_ptr;
        const unsigned char* fread_end = nullptr;
        if (unlikely(InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end))) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
          return kPglRetReadFail;
        }
        const uint32_t last_byte_in_record = fread_end[-1];
        if (unlikely(last_byte_in_record >> dbl_sample_ct_mod4)) {
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
    if (unlikely(VrtypeLdCompressed(vrtypes[vidx]))) {
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
    if (unlikely(fseeko(ff, 0, SEEK_END))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
      return kPglRetReadFail;
    }
    fsize = ftello(ff);
    pgrp->fp_vidx = 1;  // force fseek when loading first variant
#ifndef NO_MMAP
  }
#endif
  // todo: modify this check when phase sets are implemented
  const uint64_t expected_fsize = pgrp->fi.var_fpos[variant_ct];
  if (unlikely(expected_fsize != fsize)) {
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
    for (uint32_t vblock_idx = 0; vblock_idx != vblock_ct; ++vblock_idx) {
      if (unlikely(fpos_index[vblock_idx] != pgrp->fi.var_fpos[vblock_idx * kPglVblockSize])) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen header vblock-start index is inconsistent with variant record length index.\n");
        return kPglRetMalformedInput;
      }
    }
  } else {
#endif
    if (unlikely(fseeko(ff, 11, SEEK_SET))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
      return kPglRetReadFail;
    }
    header_ctrl = getc_unlocked(ff);
    if (unlikely(header_ctrl > 255)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
      return kPglRetReadFail;
    }
    for (uint32_t vblock_idx = 0; vblock_idx != vblock_ct; ++vblock_idx) {
      uint64_t vblock_start_fpos;
      if (unlikely(!fread_unlocked(&vblock_start_fpos, sizeof(int64_t), 1, ff))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
        return kPglRetReadFail;
      }
      if (unlikely(vblock_start_fpos != pgrp->fi.var_fpos[vblock_idx * kPglVblockSize])) {
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
      if (unlikely(fseeko(ff, last_vrtype_byte_offset, SEEK_SET))) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
        return kPglRetReadFail;
      }
      last_vrtype_byte = getc_unlocked(ff);
      if (unlikely(last_vrtype_byte > 255)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
        return kPglRetReadFail;
      }
#ifndef NO_MMAP
    }
#endif
    if (unlikely(last_vrtype_byte >> trailing_shift)) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Nonzero trailing bits in last vrtype index byte.\n");
      return kPglRetMalformedInput;
    }
  }
  const uintptr_t* nonref_flags = pgrp->fi.nonref_flags;
  if (nonref_flags) {
    const uint32_t variant_ct_modl = variant_ct % kBitsPerWord;
    if (variant_ct % CHAR_BIT) {
      if (unlikely(nonref_flags[variant_ct / kBitsPerWord] >> variant_ct_modl)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Nonzero trailing bits in last nonref_flags byte.\n");
        return kPglRetMalformedInput;
      }
    }
  }

  // could move most of this into plink2_common and make it multithreaded, if
  // speed is ever an issue.
  uint32_t allele_ct = 2;
  for (uint32_t vidx = 0; vidx != variant_ct; ++vidx) {
    const unsigned char* fread_ptr;
    const unsigned char* fread_end;
    if (unlikely(InitReadPtrs(vidx, pgrp, &fread_ptr, &fread_end))) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: .pgen read failure: %s.\n", strerror(errno));
      return kPglRetReadFail;
    }
    const unsigned char* fread_ptr_start = fread_ptr;
    if (unlikely(ValidateGeno(fread_end, vidx, pgrp, &fread_ptr, genovec_buf, errstr_buf))) {
      return kPglRetMalformedInput;
    }
    ZeroTrailingQuaters(sample_ct, genovec_buf);
    const uint32_t vrtype = vrtypes[vidx];
    uint32_t het_ct = CountQuater(genovec_buf, kMask5555, sample_ct);
    if (allele_idx_offsets) {
      allele_ct = allele_idx_offsets[vidx + 1] - allele_idx_offsets[vidx];
    }
    if (VrtypeMultiallelicHc(vrtype)) {
      if (unlikely(ValidateMultiallelicHc(fread_end, genovec_buf, vidx, allele_ct, pgrp, &fread_ptr, &het_ct, errstr_buf))) {
        return kPglRetMalformedInput;
      }
    }
    // don't need genovec_buf to store main genotypes past this point.
    if (VrtypeHphase(vrtype)) {
      if (unlikely(ValidateHphase(fread_end, vidx, het_ct, &fread_ptr, errstr_buf))) {
        return kPglRetMalformedInput;
      }
    }
    if (vrtype & 0xe0) {
      if (unlikely((vrtype & 0xe0) == 0x80)) {
        snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid record type for (0-based) variant #%u (phased dosage bit set, but main dosage bits unset).\n", vidx);
        return kPglRetMalformedInput;
      }
      PglErr reterr = ValidateDosage16(fread_end, vidx, pgrp, &fread_ptr, errstr_buf);
      if (unlikely(reterr)) {
        return reterr;
      }
    }
    if (unlikely(fread_ptr != fread_end)) {
      // possible todo: tolerate this at the end of a vblock.
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Extra byte(s) in (0-based) variant record #%u. (record type = %u; expected length = %" PRIuPTR ", actual = %" PRIuPTR ")\n", vidx, vrtype, S_CAST(uintptr_t, fread_ptr - fread_ptr_start), S_CAST(uintptr_t, fread_end - fread_ptr_start));
      return kPglRetMalformedInput;
    }
  }
  return kPglRetSuccess;
}


BoolErr CleanupPgfi(PgenFileInfo* pgfip, PglErr* reterrp) {
  // memory is the responsibility of the caller
  if (pgfip->shared_ff) {
    if (unlikely(fclose_null(&pgfip->shared_ff))) {
      if (*reterrp == kPglRetSuccess) {
        *reterrp = kPglRetReadFail;
        return 1;
      }
    }
#ifndef NO_MMAP
  } else if (pgfip->block_base != nullptr) {
    munmap(K_CAST(unsigned char*, pgfip->block_base), pgfip->file_size);
#endif
  }
  return 0;
}

BoolErr CleanupPgr(PgenReader* pgrp, PglErr* reterrp) {
  // assume file is open if pgr.ff is not null
  // memory is the responsibility of the caller for now
  if (!pgrp->ff) {
    return 0;
  }
  if (fclose_null(&(pgrp->ff))) {
    if (*reterrp == kPglRetSuccess) {
      *reterrp = kPglRetReadFail;
      return 1;
    }
  }
  return 0;
}


// ***** end PgenReader, begin {st,mt}_pgen_writer_t *****


void PreinitSpgw(STPgenWriter* spgwp) {
  spgwp->pgen_outfile = nullptr;
}

PglErr PwcInitPhase1(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, PgenWriterCommon* pwcp, FILE** pgen_outfile_ptr) {
  pwcp->allele_idx_offsets = allele_idx_offsets;
  pwcp->explicit_nonref_flags = nullptr;
  if (nonref_flags_storage == 3) {
    if (unlikely(!explicit_nonref_flags)) {
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
  pwcp->genovec_hets_buf = nullptr;
  pwcp->genovec_invert_buf = nullptr;
  pwcp->ldbase_genovec = nullptr;
  pwcp->ldbase_raregeno = nullptr;
  pwcp->ldbase_difflist_sample_ids = nullptr;
#endif
  pwcp->vidx = 0;

  FILE* pgen_outfile = fopen(fname, FOPEN_WB);
  *pgen_outfile_ptr = pgen_outfile;
  if (unlikely(!pgen_outfile)) {
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
  if (unlikely(fwrite_checked(zerobuf, header_bytes_left, pgen_outfile))) {
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
  const uintptr_t vrec_len_byte_ct = BytesToRepresentNzU32(max_vrec_len);
  cachelines_required += DivUp((variant_ct - 1) * vrec_len_byte_ct + sizeof(int32_t), kCacheline);

  // vrtype_buf
  if (phase_dosage_gflags) {
    cachelines_required += DivUp(variant_ct, kCacheline);
  } else {
    cachelines_required += DivUp(variant_ct, kCacheline * 2);
  }

  // genovec_hets_buf, genovec_invert_buf, ldbase_genovec
  cachelines_required += 3 * QuaterCtToCachelineCt(sample_ct);

  const uint32_t max_difflist_len = 2 * (sample_ct / kPglMaxDifflistLenDivisor);
  // ldbase_raregeno
  cachelines_required += QuaterCtToCachelineCt(max_difflist_len);

  // ldbase_difflist_sample_ids
  cachelines_required += 1 + (max_difflist_len / kInt32PerCacheline);

  // fwrite_buf
  // + (5 + sizeof(AlleleCode)) * kPglDifflistGroupSize to avoid buffer
  // overflow in middle of difflist writing
  cachelines_required += DivUp(max_vrec_len + kPglFwriteBlockSize + (5 + sizeof(AlleleCode)) * kPglDifflistGroupSize, kCacheline);
  // possible todo: dosage (doesn't currently need an allocation, but that's
  // unlikely to remain true--e.g. get_ref_nonref_genotype_counts_and_dosages
  // tends to use workspace_vec when a function it calls doesn't use it...)

  // probable todo: multiallelic (phased) dosage
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
    // see comments in middle of MpgwInitPhase1()
    max_vrec_len += 2 + sizeof(AlleleCode) + GetAux1bAlleleEntryByteCt(max_alt_ct_p1, sample_ct - 1);
    // try to permit uncompressed records to be larger than this, only error
    // out when trying to write a larger compressed record.
  }
  if (phase_dosage_gflags & kfPgenGlobalHardcallPhasePresent) {
    // phasepresent, phaseinfo
    max_vrec_len += 2 * DivUp(sample_ct, CHAR_BIT);
  }
  if (phase_dosage_gflags & kfPgenGlobalDosagePresent) {
    const uint32_t dphase_gflag = (phase_dosage_gflags / kfPgenGlobalDosagePhasePresent) & 1;
    // aux3, aux7
    max_vrec_len += (1 + dphase_gflag) * DivUp(sample_ct, 8);
    // aux4, aux8
    max_vrec_len += (2 + 2 * dphase_gflag) * S_CAST(uint64_t, sample_ct);
    // todo: multiallelic dosage
  }
#ifdef __LP64__
  if (max_vrec_len >= kPglMaxBytesPerVariant) {
    max_vrec_len = kPglMaxBytesPerVariant;
  }
#else
  if (unlikely(max_vrec_len > kMaxBytesPerIO + 1 - kPglFwriteBlockSize)) {
    return kPglRetNomem;
  }
#endif
  *max_vrec_len_ptr = max_vrec_len;
  const uintptr_t vrec_len_byte_ct = BytesToRepresentNzU32(max_vrec_len);

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

  // genovec_hets_buf, genovec_invert_buf, ldbase_genovec
  // (could make genovec_hets_buf allocation conditional, but not worth the
  // additional load this puts on import functions)
  alloc_per_thread_cacheline_ct += 3 * QuaterCtToCachelineCt(sample_ct);

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
  if (allele_idx_offsets && (allele_idx_offsets[variant_ct] != 2 * variant_ct)) {
    assert(allele_idx_offsets[0] == 0);
    assert(allele_idx_offsets[variant_ct] > 2 * variant_ct);
    // When multiallelic variants are present, larger write buffers are needed.
    // we compute the largest possible size here.
    //
    // One way to hit the maximum size is to have exactly 1 ref/altx genotype,
    // all other genotypes are altx/alty, and all genotypes include a rarealt.
    // Then, we need 2 + sizeof(AlleleCode) bytes for the format byte and
    // aux1a, (sample_ct + 6) / 8 bytes for the bitarray at the front of aux1b,
    // and
    //   alt ct  aux1_last_quarter bytes required
    //   ------  --------------------------------
    //        2           (sample_ct - 1 + 7) / 8
    //      3-4           (sample_ct - 1 + 1) / 2
    //     5-16                   (sample_ct - 1)
    //   17-256               2 * (sample_ct - 1)
    //
    // (todo: also take multiallelic dosage into account)
    assert(!dosage_gflag);
    uintptr_t prev_offset = 0;
    uint32_t vidx = 0;
    const uint32_t extra_bytes_base = 2 + sizeof(AlleleCode) + (sample_ct + 6) / 8;
    const uint64_t extra_bytes_max = kPglMaxBytesPerVariant - max_vrec_len;
    uint64_t extra_byte_cts[4];
    uint32_t extra_alt_ceil = kPglMaxAltAlleleCt + 1;

    // alt_ct == 2
    uint64_t cur_extra_byte_ct = extra_bytes_base + (sample_ct + 6) / 8;
    extra_byte_cts[0] = cur_extra_byte_ct;
    extra_byte_cts[1] = 0;  // force initialization
    extra_byte_cts[2] = 0;
    extra_byte_cts[3] = 0;
    if (cur_extra_byte_ct >= extra_bytes_max) {
      extra_alt_ceil = 2;
    } else {
      // alt_ct in [3, 4]
      cur_extra_byte_ct = extra_bytes_base + (sample_ct / 2);
      extra_byte_cts[1] = cur_extra_byte_ct;
      if (cur_extra_byte_ct >= extra_bytes_max) {
        extra_alt_ceil = 3;
      } else {
        // alt_ct in [5, 16]
        cur_extra_byte_ct = extra_bytes_base + sample_ct - 1;
        extra_byte_cts[2] = cur_extra_byte_ct;
        if (cur_extra_byte_ct >= extra_bytes_max) {
          extra_alt_ceil = 5;
        } else {
          // alt_ct in [17, 256]
          cur_extra_byte_ct = extra_bytes_base + 2 * (sample_ct - 1);
          extra_byte_cts[3] = cur_extra_byte_ct;
          if (cur_extra_byte_ct >= extra_bytes_max) {
            extra_alt_ceil = 16;
          }
        }
      }
    }
    uint32_t extra_alt_ceil_ct = 0;
    const uint64_t uncompressed_biallelic_vrec_len = max_vrec_len;
    uint32_t altx_seen_mask = 0;
    uint32_t max_alt_ct_p1 = 3;
    for (uint32_t vblock_end = vidx + kPglVblockSize; ; vblock_end += kPglVblockSize) {
      if (vblock_end > variant_ct) {
        if (vidx == variant_ct) {
          break;
        }
        vblock_end = variant_ct;
      }
      uint64_t cur_vblock_byte_ct = uncompressed_biallelic_vrec_len * (vblock_end - vidx);
      uint32_t altx_seen[4];
      ZeroU32Arr(4, altx_seen);
      for (; vidx < vblock_end; ) {
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

            if (alt_ct_p1 < 6) {
              altx_seen[(alt_ct_p1 != 3)] += 1;
            } else {
              altx_seen[2 + (alt_ct_p1 >= 17)] += 1;
            }
          }
        }
        prev_offset = cur_offset;
      }
      cur_vblock_byte_ct += extra_alt_ceil_ct * extra_bytes_max;
      for (uint32_t uii = 0; uii != 4; ++uii) {
        if (altx_seen[uii]) {
          const uint32_t cur_seen_ct = altx_seen[uii];
          altx_seen_mask |= 1 << uii;
          cur_vblock_byte_ct += cur_seen_ct * extra_byte_cts[uii];
        }
      }
      if (cur_vblock_byte_ct > max_vblock_byte_ct) {
        max_vblock_byte_ct = cur_vblock_byte_ct;
      }
    }
    if (extra_alt_ceil_ct) {
      max_vrec_len = kPglMaxBytesPerVariant;
    } else {
      max_vrec_len = uncompressed_biallelic_vrec_len + extra_byte_cts[bsru32(altx_seen_mask)];
    }
  }
  if (max_vrec_len >= kPglMaxBytesPerVariant) {
    max_vrec_len = kPglMaxBytesPerVariant;
    max_vblock_byte_ct = kPglMaxBytesPerVariant * S_CAST(uint64_t, max_vblock_size);
  }

  // vrec_len_buf
  // previously used overlapping uint32_t writes-to-memory, but that was
  // incompatible with multithreaded compression
  *vrec_len_byte_ct_ptr = BytesToRepresentNzU32(max_vrec_len);
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
  for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
    if (tidx) {
      pwcs[tidx]->vblock_fpos = pwcs[0]->vblock_fpos;
      pwcs[tidx]->vrec_len_buf = pwcs[0]->vrec_len_buf;
      pwcs[tidx]->vrtype_buf = pwcs[0]->vrtype_buf;
    }
    pwcs[tidx]->genovec_hets_buf = R_CAST(uintptr_t*, alloc_iter);
    alloc_iter = &(alloc_iter[genovec_byte_alloc]);
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
  const uintptr_t fwrite_cacheline_ct = DivUp(max_vrec_len + kPglFwriteBlockSize + (5 + sizeof(AlleleCode)) * kPglDifflistGroupSize, kCacheline);
  PgenWriterCommon* pwcp = &(spgwp->pwc);
  PwcInitPhase2(fwrite_cacheline_ct, 1, &pwcp, spgw_alloc);
}

PglErr MpgwInitPhase2(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, uintptr_t vblock_cacheline_ct, uint32_t thread_ct, unsigned char* mpgw_alloc, MTPgenWriter* mpgwp) {
  assert(thread_ct);
  const uintptr_t pwc_byte_ct = RoundUpPow2(sizeof(PgenWriterCommon), kCacheline);
  for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
    mpgwp->pwcs[tidx] = R_CAST(PgenWriterCommon*, &(mpgw_alloc[tidx * pwc_byte_ct]));
  }
  PglErr reterr = PwcInitPhase1(fname, allele_idx_offsets, explicit_nonref_flags, variant_ct, sample_ct, phase_dosage_gflags, nonref_flags_storage, vrec_len_byte_ct, mpgwp->pwcs[0], &(mpgwp->pgen_outfile));
  if (unlikely(reterr)) {
    return reterr;
  }
  mpgwp->thread_ct = thread_ct;
  for (uint32_t tidx = 1; tidx != thread_ct; ++tidx) {
    memcpy(mpgwp->pwcs[tidx], mpgwp->pwcs[0], sizeof(PgenWriterCommon));
    mpgwp->pwcs[tidx]->vidx = tidx * kPglVblockSize;
  }
  PwcInitPhase2(vblock_cacheline_ct, thread_ct, mpgwp->pwcs, &(mpgw_alloc[thread_ct * pwc_byte_ct]));
  return kPglRetSuccess;
}


void CountLdAndInvertedLdDiffs(const uintptr_t* __restrict ldbase_genovec, const uintptr_t* __restrict genovec, uint32_t sample_ct, uint32_t* ld_diff_ctp, uint32_t* ld_inv_diff_ctp) {
  // Requires trailing bits to be zeroed out.
  const uint32_t word_ct = QuaterCtToWordCt(sample_ct);
  const uintptr_t* genovec_end = &(genovec[word_ct]);
  uint32_t ld_diff_ct = 0;
  uint32_t ld_inv_diff_ct = 0;
  // construct the words we want to popcount_quatervec_01 on the fly
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* ldbase_vvec_iter = R_CAST(const VecW*, ldbase_genovec);
  const VecW* geno_vvec_iter = R_CAST(const VecW*, genovec);
  VecW acc_ld = vecw_setzero();
  VecW acc_ld_inv = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (uintptr_t full_vecs_left = 3 * (word_ct / (3 * kWordsPerVec)); ; full_vecs_left -= cur_incr) {
    if (full_vecs_left < 60) {
      if (!full_vecs_left) {
        ld_diff_ct = HsumW(acc_ld);
        ld_inv_diff_ct = HsumW(acc_ld_inv);
        break;
      }
      cur_incr = full_vecs_left;
    }
    VecW inner_acc_ld = vecw_setzero();
    VecW inner_acc_ld_inv = vecw_setzero();
    const VecW* geno_vvec_stop = &(geno_vvec_iter[cur_incr]);
    do {
      VecW loader_ldbase1 = *ldbase_vvec_iter++;
      VecW loader_geno1 = *geno_vvec_iter++;
      VecW loader_ldbase2 = *ldbase_vvec_iter++;
      VecW loader_geno2 = *geno_vvec_iter++;
      VecW xor1 = loader_ldbase1 ^ loader_geno1;
      VecW xor2 = loader_ldbase2 ^ loader_geno2;
      VecW xor_shifted1 = vecw_srli(xor1, 1);
      VecW xor_shifted2 = vecw_srli(xor2, 1);
      // xor(_low)  xor_shifted  loader_geno   result
      //         1                                  1
      //         0            0            0        1
      //         0            0            1        0
      //         0            1            0        0
      //         0            1            1        1
      // gah, don't see a way to avoid throwing in an extra xor for
      // loader_geno...
      VecW count_ld_inv = (xor1 | (xor_shifted1 ^ loader_geno1 ^ m1)) & m1;
      loader_ldbase1 = *ldbase_vvec_iter++;
      VecW count_ld = (xor1 | xor_shifted1) & m1;
      loader_geno1 = *geno_vvec_iter++;
      count_ld_inv = count_ld_inv + ((xor2 | (xor_shifted2 ^ loader_geno2 ^ m1)) & m1);
      xor1 = loader_ldbase1 ^ loader_geno1;
      count_ld = count_ld + ((xor2 | xor_shifted2) & m1);
      xor_shifted1 = vecw_srli(xor1, 1);
      count_ld_inv = count_ld_inv + ((xor1 | (xor_shifted1 ^ loader_geno1 ^ m1)) & m1);
      count_ld = count_ld + ((xor1 | xor_shifted1) & m1);
      // now count_ld and count_ld_inv each have 64 2-bit values from 0-3

      count_ld_inv = (count_ld_inv & m2) + (vecw_srli(count_ld_inv, 2) & m2);
      count_ld = (count_ld & m2) + (vecw_srli(count_ld, 2) & m2);
      // now they have 32 4-bit values from 0-6

      inner_acc_ld_inv = inner_acc_ld_inv + ((count_ld_inv + vecw_srli(count_ld_inv, 4)) & m4);
      inner_acc_ld = inner_acc_ld + ((count_ld + vecw_srli(count_ld, 4)) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_ld_inv = acc_ld_inv + vecw_bytesum(inner_acc_ld_inv, m0);
    acc_ld = acc_ld + vecw_bytesum(inner_acc_ld, m0);
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

void CountLdAndInvertedLdDiffsList(const uintptr_t* __restrict ldbase_raregeno, const uint32_t* __restrict ldbase_difflist_sample_ids, const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t ldbase_difflist_len, uint32_t difflist_len, uint32_t* ld_diff_ctp, uint32_t* ld_inv_diff_ctp) {
  // assumes ldbase_difflist_sample_ids[ldbase_difflist_len] == sample_ct

  // only the count(s) with aligned common_geno values are valid.  e.g. if
  // ldbase_common_geno and difflist_common_geno are both zero, the ld_inv_diff
  // return value can be anything, while if they're both three, ld_diff and
  // ld_inv_diff are both accurate.

  // some similarities to ParseLdAndMergeDifflistSubset(), but much simpler.
  // noticeably slower than CountLdAndInvertedLdDiffs() when the lists aren't
  // tiny.
  // possible todo: take threshold into account?

  uint32_t collision_ct = 0;
  uint32_t ld_diff_ct = 0;
  uint32_t ld_inv_diff_ct = 0;
  uint32_t ldbase_sample_idx = ldbase_difflist_sample_ids[0];
  uint32_t ldbase_difflist_idx = 1;
  // this loop is a bit slow.  attempt to bail halfway through?
  for (uint32_t difflist_idx = 0; difflist_idx != difflist_len; ++difflist_idx) {
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
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
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
  for (uint32_t widx = 0; ; ++widx) {
    const uintptr_t cur_geno_word = genovec[widx];
    uintptr_t xor_word = ldbase_genovec? ldbase_genovec[widx] : common_geno_word;
    xor_word ^= cur_geno_word;
    if (xor_word) {
      const uint32_t sample_idx_base = widx * kBitsPerWordD2;
      do {
        const uint32_t sample_idx_lowbits = ctzw(xor_word) / 2;
        const uint32_t new_sample_idx = sample_idx_base + sample_idx_lowbits;
        raregeno_word |= ((cur_geno_word >> (2 * sample_idx_lowbits)) & 3) << (2 * (difflist_idx % kBitsPerWordD2));
        if (!(difflist_idx % kPglDifflistGroupSize)) {
          group_first_sample_ids_iter = memcpyua(group_first_sample_ids_iter, &new_sample_idx, sample_id_byte_ct);
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
  }
}

void OnebitPreprocessBuf(const uintptr_t* __restrict genovec, uint32_t sample_ct, uint32_t common2_code, uintptr_t* __restrict genovec_buf) {
  assert(sample_ct);
  const uint32_t vec_ct = QuaterCtToVecCt(sample_ct);
  // todo: look for better ways to perform some of these operations
  const VecW* geno_vvec_iter = R_CAST(const VecW*, genovec);
  const VecW* geno_vvec_end = &(geno_vvec_iter[vec_ct]);
  VecW* write_iter = R_CAST(VecW*, genovec_buf);
  const VecW m1 = VCONST_W(kMask5555);
  if (common2_code < 5) {
    if (common2_code == 1) {
      // 11 -> 10, everything else unchanged
      // todo: check if these loops are actually faster as simple while loops
      // todo: check if it's better to unroll these loops to process 2 __m128is
      //       at a time
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        *write_iter++ = (~(m1 & vecw_srli(cur_geno, 1))) & cur_geno;
      } while (geno_vvec_iter < geno_vvec_end);
    } else if (common2_code == 3) {
      // 00 -> 00, 01 -> 10, 10 -> 10, 11 -> 01
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        const VecW cur_geno_rshift = vecw_srli(cur_geno, 1);
        const VecW cur_geno_xor_masked = (cur_geno ^ cur_geno_rshift) & m1;
        const VecW cur_geno_or_masked = (cur_geno | cur_geno_rshift) & m1;
        *write_iter++ = cur_geno_xor_masked + cur_geno_or_masked;
      } while (geno_vvec_iter < geno_vvec_end);
    } else {
      assert(common2_code == 2);
      // 00 -> 00, 01 -> 10, 10 -> 01, 11 -> 10
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        const VecW cur_geno_or_masked = (cur_geno | vecw_srli(cur_geno, 1)) & m1;
        const VecW cur_geno_lowbits = cur_geno & m1;
        *write_iter++ = cur_geno_lowbits + cur_geno_or_masked;
      } while (geno_vvec_iter < geno_vvec_end);
    }
  } else {
    if (common2_code == 5) {
      // 00 -> 10, 01 -> 00, 10 -> 01, 11 -> 10
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        const VecW cur_geno_rshift = vecw_srli(cur_geno, 1);
        const VecW cur_geno_not_xor_masked = (~(cur_geno ^ cur_geno_rshift)) & m1;
        const VecW cur_geno_rshift_masked = cur_geno_rshift & m1;
        *write_iter++ = cur_geno_not_xor_masked + (cur_geno_not_xor_masked | cur_geno_rshift_masked);
      } while (geno_vvec_iter < geno_vvec_end);
    } else if (common2_code == 9) {
      // 00 -> 10, 01 -> 10, 10 -> 00, 11 -> 01
      const VecW not_m1 = VCONST_W(kMaskAAAA);
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        *write_iter++ = (cur_geno ^ not_m1) - ((~not_m1) & ((~vecw_srli(cur_geno, 1)) & cur_geno));
      } while (geno_vvec_iter < geno_vvec_end);
    } else {
      assert(common2_code == 6);
      // 00 -> 10, 01 -> 00, 10 -> 10, 11 -> 01
      do {
        const VecW cur_geno = *geno_vvec_iter++;
        const VecW cur_geno_not_lowbits = (~cur_geno) & m1;
        const VecW cur_geno_rshift_masked = vecw_srli(cur_geno, 1) & m1;
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
  PackWordsToHalfwordsMask(genovec_buf, word_read_ct, fwrite_bufp_alias_halfword);
  const uint32_t onebit_block_len = (sample_ct + 15) / CHAR_BIT;
  unsigned char* fwrite_bufp = Vint32Append(onebit_difflist_len, &(fwrite_bufp_start[onebit_block_len]));
  // the rest is almost identical to SaveLdDifflist()
  if (!onebit_difflist_len) {
    pwcp->fwrite_bufp = fwrite_bufp;
    return (onebit_block_len + 1);
  }
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(onebit_difflist_len, kPglDifflistGroupSize);
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
  uintptr_t* raregeno_iter = R_CAST(uintptr_t*, &(extra_byte_cts_iter[group_ct - 1]));
  fwrite_bufp = &(extra_byte_cts_iter[group_ct + (onebit_difflist_len - 1) / 4]);
  unsigned char* last_group_vint_start = fwrite_bufp;
  uintptr_t raregeno_word = 0;
  uint32_t last_sample_idx = 0;
  uint32_t difflist_idx = 0;
  for (uint32_t widx = 0; ; ++widx) {
    uintptr_t xor_word = genovec_buf[widx] & kMaskAAAA;
    if (xor_word) {
      const uintptr_t cur_geno_word = genovec[widx];
      const uint32_t sample_idx_base = widx * kBitsPerWordD2;

      do {
        const uint32_t sample_idx_lowbits = ctzw(xor_word) / 2;
        const uint32_t new_sample_idx = sample_idx_base + sample_idx_lowbits;
        if (!(difflist_idx % kBitsPerWordD2)) {
          if (!(difflist_idx % kPglDifflistGroupSize)) {
            SubU32StoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
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
        xor_word &= xor_word - 1;
      } while (xor_word);
      // trailing bits of genovec_buf guaranteed to be zeroed out
      if (difflist_idx == onebit_difflist_len) {
        SubwordStore(raregeno_word, 1 + (((onebit_difflist_len - 1) / 4) % sizeof(intptr_t)), raregeno_iter);
        pwcp->fwrite_bufp = fwrite_bufp;
        return fwrite_bufp - fwrite_bufp_start;
      }
    }
  }
}

// returns vrec_len
uint32_t PwcAppendBiallelicGenovecMain(const uintptr_t* __restrict genovec, uint32_t vidx, PgenWriterCommon* pwcp, uint32_t* het_ct_ptr, uint32_t* altxy_ct_ptr, unsigned char* vrtype_ptr) {
  const uint32_t sample_ct = pwcp->sample_ct;
  assert((!(sample_ct % kBitsPerWordD2)) || (!(genovec[sample_ct / kBitsPerWordD2] >> (2 * (sample_ct % kBitsPerWordD2)))));
  STD_ARRAY_DECL(uint32_t, 4, genocounts);
  GenovecCountFreqsUnsafe(genovec, sample_ct, genocounts);
  if (het_ct_ptr) {
    *het_ct_ptr = genocounts[1];
    if (altxy_ct_ptr) {
      *altxy_ct_ptr = genocounts[2];
    }
  }
  uint32_t most_common_geno = (genocounts[1] > genocounts[0]);
  uint32_t second_most_common_geno = 1 - most_common_geno;
  uint32_t largest_geno_ct = genocounts[most_common_geno];
  uint32_t second_largest_geno_ct = genocounts[second_most_common_geno];
  for (uint32_t cur_geno = 2; cur_geno != 4; ++cur_geno) {
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
  STD_ARRAY_REF(uint32_t, 4) ldbase_genocounts = pwcp->ldbase_genocounts;
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
  STD_ARRAY_COPY(genocounts, 4, ldbase_genocounts);
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
#ifdef FUTURE_ENCODER
    if ((!difflist_len) && (!most_common_geno)) {
      *vrtype_ptr = 5;
      return 0;
    }
#endif
    *vrtype_ptr = 4 + most_common_geno;
    return SaveLdDifflist(genovec, nullptr, most_common_geno, difflist_len, pwcp);
  }
  *vrtype_ptr = 0;
  const uint32_t vrec_len = QuaterCtToByteCt(sample_ct);
  pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, genovec, vrec_len);
  return vrec_len;
}

void PwcAppendBiallelicGenovec(const uintptr_t* __restrict genovec, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  const uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, nullptr, nullptr, &vrtype);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  pwcp->vidx += 1;
  SubU32Store(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]));
  // could have a single expression which branchlessly handles both cases, but
  // doubt that's worthwhile
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
}

BoolErr SpgwFlush(STPgenWriter* spgwp) {
  if (spgwp->pwc.fwrite_bufp >= &(spgwp->pwc.fwrite_buf[kPglFwriteBlockSize])) {
    const uintptr_t cur_byte_ct = spgwp->pwc.fwrite_bufp - spgwp->pwc.fwrite_buf;
    if (unlikely(fwrite_checked(spgwp->pwc.fwrite_buf, cur_byte_ct, spgwp->pgen_outfile))) {
      return 1;
    }
    spgwp->pwc.vblock_fpos_offset += cur_byte_ct;
    spgwp->pwc.fwrite_bufp = spgwp->pwc.fwrite_buf;
  }
  return 0;
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
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
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
      SubU32StoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
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
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t* __restrict difflist_sample_ids = pwcp->ldbase_difflist_sample_ids;
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
  fwrite_bufp = memcpyua(&(extra_byte_cts_iter[group_ct - 1]), pwcp->ldbase_raregeno, QuaterCtToByteCt(difflist_len));
  unsigned char* last_group_vint_start = nullptr;
  uint32_t last_sample_idx = 0;
  for (uint32_t difflist_idx = 0; difflist_idx != difflist_len; ++difflist_idx) {
    const uint32_t new_sample_idx = difflist_sample_ids[difflist_idx];
    if (!(difflist_idx % kPglDifflistGroupSize)) {
      SubU32StoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
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
  STD_ARRAY_DECL(uint32_t, 4, genocounts);
  GenovecCountFreqsUnsafe(raregeno, difflist_len, genocounts);
  assert(!genocounts[difflist_common_geno]);
  genocounts[difflist_common_geno] = sample_ct - difflist_len;
  uint32_t second_most_common_geno = difflist_common_geno? 0 : 1;
  uint32_t second_largest_geno_ct = genocounts[second_most_common_geno];
  for (uint32_t cur_geno = second_most_common_geno + 1; cur_geno != 4; ++cur_geno) {
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
  STD_ARRAY_REF(uint32_t, 4) ldbase_genocounts = pwcp->ldbase_genocounts;
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
        CountLdAndInvertedLdDiffsList(pwcp->ldbase_raregeno, pwcp->ldbase_difflist_sample_ids, raregeno, difflist_sample_ids, pwcp->ldbase_difflist_len, difflist_len, &ld_diff_ct, &ld_inv_diff_ct);
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
  STD_ARRAY_COPY(genocounts, 4, ldbase_genocounts);
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
  pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, genobuf, vrec_len);
  return vrec_len;
}

void PwcAppendBiallelicDifflistLimited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  const uint32_t vrec_len = PwcAppendBiallelicDifflistLimitedMain(raregeno, difflist_sample_ids, vidx, difflist_common_geno, difflist_len, pwcp, &vrtype);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  pwcp->vidx += 1;
  SubU32Store(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]));
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
}


static inline BoolErr CheckedVrecLenIncr(uintptr_t incr, uint32_t* vrec_len_ptr) {
  // maybe track vrec_left instead of vrec_len...
#ifdef __LP64__
  const uintptr_t new_vrec_len = *vrec_len_ptr + incr;
  if (unlikely(new_vrec_len > kPglMaxBytesPerVariant)) {
    return 1;
  }
  *vrec_len_ptr = new_vrec_len;
  return 0;
#else
  *vrec_len_ptr += incr;
  return 0;
#endif
}

// ok if delta_bitarr trailing bits set
BoolErr PwcAppendDeltalist(const uintptr_t* delta_bitarr, uint32_t deltalist_len, PgenWriterCommon* pwcp, uint32_t* vrec_len_ptr) {
  assert(deltalist_len);
  unsigned char* fwrite_bufp = pwcp->fwrite_bufp;
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(deltalist_len, fwrite_bufp);
  const uint32_t sample_id_byte_ct = BytesToRepresentNzU32(pwcp->sample_ct);
  const uint32_t group_ct = DivUp(deltalist_len, kPglDifflistGroupSize);
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  if (unlikely(CheckedVrecLenIncr(S_CAST(uintptr_t, fwrite_bufp - fwrite_bufp_start) + group_ct * sample_id_byte_ct + (group_ct - 1), vrec_len_ptr))) {
    return 1;
  }
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * sample_id_byte_ct]);
  fwrite_bufp_start = &(extra_byte_cts_iter[group_ct - 1]);
  fwrite_bufp = fwrite_bufp_start;
#ifdef __LP64__
  const intptr_t vrec_left = kPglMaxBytesPerVariant - (*vrec_len_ptr);
#endif
  unsigned char* last_group_vint_start = nullptr;
  uintptr_t last_sample_idx = 0;
  uintptr_t new_sample_idx_base = 0;
  uintptr_t delta_bits = delta_bitarr[0];
  for (uint32_t deltalist_idx = 0; deltalist_idx != deltalist_len; ++deltalist_idx) {
    const uintptr_t new_sample_idx = BitIter1(delta_bitarr, &new_sample_idx_base, &delta_bits);
    if (!(deltalist_idx % kPglDifflistGroupSize)) {
      SubU32StoreMov(new_sample_idx, sample_id_byte_ct, &group_first_sample_ids_iter);
      if (deltalist_idx) {
        *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
      }
#ifdef __LP64__
      // Note that this inequality goes in the opposite direction from the
      // usual PtrCheck.
      if (unlikely(S_CAST(intptr_t, fwrite_bufp - fwrite_bufp_start) > vrec_left)) {
        return 1;
      }
#endif
      last_group_vint_start = fwrite_bufp;
    } else {
      assert(new_sample_idx >= last_sample_idx + 1);
      fwrite_bufp = Vint32Append(new_sample_idx - last_sample_idx, fwrite_bufp);
    }
    last_sample_idx = new_sample_idx;
  }
  pwcp->fwrite_bufp = fwrite_bufp;
  return CheckedVrecLenIncr(fwrite_bufp - fwrite_bufp_start, vrec_len_ptr);
}

static_assert(sizeof(AlleleCode) == 1, "PwcAppendMultiallelicMain() needs to be updated.");
BoolErr PwcAppendMultiallelicMain(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t patch_01_ct, uint32_t patch_10_ct, uint32_t vidx, PgenWriterCommon* pwcp, const uintptr_t** genovec_hets_ptr, uint32_t* het_ct_ptr, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  uint32_t genovec_het_ct;
  uint32_t genovec_altxy_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &genovec_het_ct, &genovec_altxy_ct, vrtype_ptr);
  if ((!patch_01_ct) && (!patch_10_ct)) {
    if (genovec_hets_ptr) {
      *genovec_hets_ptr = genovec;
      *het_ct_ptr = genovec_het_ct;
    }
    *vrec_len_ptr = vrec_len;
    return 0;
  }
  const uint32_t allele_ct = pwcp->allele_idx_offsets[vidx + 1] - pwcp->allele_idx_offsets[vidx];
  assert(allele_ct > 2);
  unsigned char* format_byte_ptr = pwcp->fwrite_bufp;
  ++vrec_len;
  pwcp->fwrite_bufp += 1;
  unsigned char format_byte = 0;
  if (!patch_01_ct) {
    format_byte = 15;
  } else {
    // We could use a different bitarray vs. deltalist threshold, since the
    // usual threshold was chosen under the assumption that
    //   (maximum allowed deltalist index + 1) == bitarray length.
    // However, the subsetting-speed advantage of the deltalist representation
    // can be expected to make up for its larger size, so we may actually want
    // to increase instead of decrease the current 1/9 fraction.  For now, I'll
    // punt on this decision.
    const uint32_t max_deltalist_entry_ct = genovec_het_ct / kPglMaxDeltalistLenDivisor;
    if (patch_01_ct < max_deltalist_entry_ct) {
      // impossible for this to run out of space
      PwcAppendDeltalist(patch_01_set, patch_01_ct, pwcp, &vrec_len);
      format_byte = 1;
    } else {
      const uint32_t genovec_het_ctb = DivUp(genovec_het_ct, CHAR_BIT);
#ifdef __arm__
#  error "Unaligned accesses in PwcAppendMultiallelicMain()."
#endif
      CopyGenomatchSubset(patch_01_set, genovec, kMask5555, 0, genovec_het_ct, R_CAST(uintptr_t*, pwcp->fwrite_bufp));
      pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[genovec_het_ctb]);
      vrec_len += genovec_het_ctb;
    }
    if (allele_ct > 3) {
      if (allele_ct <= 18) {
        const AlleleCode* patch_01_vals_iter = patch_01_vals;

        // er, no, we want this to be a different pointer type depending on the
        // allele count and build type.  but make those optimizations later,
        // after a simple unified version of this code has been tested.
        uintptr_t* write_alias = R_CAST(uintptr_t*, pwcp->fwrite_bufp);
        uint32_t written_byte_ct;
        if (allele_ct == 4) {
          // 1 bit per entry, clear = ref/alt2, set = ref/alt3
#ifdef __LP64__
          const uint32_t fullvec_ct = patch_01_ct / kBytesPerVec;
          Vec8thUint* write_aliasv = R_CAST(Vec8thUint*, write_alias);
          if (fullvec_ct) {
            // parallel-equality-to-3 check, movemask
            // (+125 or <<7 followed by movemask also works)
            for (uint32_t vec_idx = 0; vec_idx != fullvec_ct; ++vec_idx) {
              VecUc cur_vec = vecuc_loadu(patch_01_vals_iter);
              patch_01_vals_iter = &(patch_01_vals_iter[kBytesPerVec]);
              write_aliasv[vec_idx] = vecuc_movemask(R_CAST(VecUc, vecw_slli(R_CAST(VecW, cur_vec), 7)));
            }
          }
          const uint32_t remainder = patch_01_ct % kBytesPerVec;
          if (remainder) {
            const uintptr_t* cur_read_iter = R_CAST(const uintptr_t*, patch_01_vals_iter);
            unsigned char* write_uc = R_CAST(unsigned char*, &(write_aliasv[fullvec_ct]));
            const uint32_t word_ct_m1 = (remainder - 1) / kBytesPerWord;
            for (uint32_t widx = 0; ; ++widx) {
              uintptr_t cur_word;
              if (widx >= word_ct_m1) {
                if (widx > word_ct_m1) {
                  break;
                }
                cur_word = SubwordLoad(cur_read_iter, ModNz(remainder, kBytesPerWord));
              } else {
                cur_word = *cur_read_iter++;
              }
#  ifdef USE_AVX2
              cur_word = _pext_u64(cur_word, kMask0101);
#  else
              // extract the low bit from each byte
              cur_word = cur_word & kMask0101;
              cur_word = (cur_word * 0x102040810204080LLU) >> 56;
#  endif
              write_uc[widx] = cur_word;
            }
          }
#else  // !__LP64__
          const uint32_t word_ct_m1 = (patch_01_ct - 1) / kBitsPerWord;
          uint32_t loop_len = kBitsPerWord;
          for (uint32_t widx = 0; ; ++widx) {
            if (widx >= word_ct_m1) {
              if (widx > word_ct_m1) {
                break;
              }
              loop_len = ModNz(patch_01_ct, kBitsPerWord);
            }
            uintptr_t cur_word = 0;
            for (uint32_t uii = loop_len - 1; ; --uii) {
              cur_word |= patch_01_vals_iter[uii] & 1;
              if (!uii) {
                break;
              }
              cur_word <<= 1;
            }
            write_alias[widx] = cur_word;
            patch_01_vals_iter = &(patch_01_vals_iter[loop_len]);
          }
#endif  // !__LP64__
          written_byte_ct = DivUp(patch_01_ct, 8);
        } else if (allele_ct <= 6) {
          // 2 bits per entry
          const uint32_t word_ct_m1 = (patch_01_ct - 1) / kBitsPerWordD2;
          uint32_t loop_len = kBitsPerWordD2;
          for (uint32_t widx = 0; ; ++widx) {
            if (widx >= word_ct_m1) {
              if (widx > word_ct_m1) {
                break;
              }
              loop_len = ModNz(patch_01_ct, kBitsPerWordD2);
            }
            uintptr_t cur_word = 0;
            for (uint32_t uii = loop_len - 1; ; --uii) {
              // todo:
              // 1. parallel-subtract (now all bytes are in 0..3)
              // 2. bitwise or with right-shift-6
              // 3. bitwise or with right-shift-12
              // 4. _mm{256}_shuffle_epi8() gather from bytes 0, 4, 8, ...
              // SSE4.2:
              // 5. store result of _mm_extract_epi32(, 0)
              // AVX2:
              // 5. store results of _mm256_extract_epi32(, 0) and (, 4)
              cur_word |= patch_01_vals_iter[uii] - 2;
              if (!uii) {
                break;
              }
              cur_word <<= 2;
            }
            write_alias[widx] = cur_word;
            patch_01_vals_iter = &(patch_01_vals_iter[loop_len]);
          }
          written_byte_ct = DivUp(patch_01_ct, 4);
        } else {
          // 4 bits per entry
          const uint32_t word_ct_m1 = (patch_01_ct - 1) / kBitsPerWordD4;
          uint32_t loop_len = kBitsPerWordD4;
          for (uint32_t widx = 0; ; ++widx) {
            if (widx >= word_ct_m1) {
              if (widx > word_ct_m1) {
                break;
              }
              loop_len = ModNz(patch_01_ct, kBitsPerWordD4);
            }
            uintptr_t cur_word = 0;
            for (uint32_t uii = loop_len - 1; ; --uii) {
              // todo: replace this with parallel-subtract, followed by bitwise
              // or with right-shift-4, followed by _mm{256}_shuffle_epi8()
              // gather from even bytes, followed by _mm{256}_extract_epi64()
              // (or _mm256_permute4x64_epi64() followed by
              // _mm256_extracti128_si256()?  test this.)
              // Or in non-SSE4.2 case, use end of PackWordToHalfword logic.
              cur_word |= patch_01_vals_iter[uii] - 2;
              if (!uii) {
                break;
              }
              cur_word <<= 4;
            }
            write_alias[widx] = cur_word;
            patch_01_vals_iter = &(patch_01_vals_iter[loop_len]);
          }
          written_byte_ct = DivUp(patch_01_ct, 2);
        }
        pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[written_byte_ct]);
        vrec_len += written_byte_ct;
      } else {
        // 1 byte per entry
        // need 2-byte and 3-byte cases here for larger AlleleCode, those also
        // need to check for vrec_len overflow
        unsigned char* payload = pwcp->fwrite_bufp;
        for (uint32_t patch_idx = 0; patch_idx != patch_01_ct; ++patch_idx) {
          // todo: check if the compiler vectorizes this
          payload[patch_idx] = patch_01_vals[patch_idx] - 2;
        }
        pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[patch_01_ct]);
        vrec_len += patch_01_ct;
      }
    }
  }
  if (!patch_10_ct) {
    format_byte |= 240;
    if (genovec_hets_ptr) {
      *genovec_hets_ptr = genovec;
      *het_ct_ptr = genovec_het_ct;
    }
  } else {
    const uint32_t max_deltalist_entry_ct = genovec_altxy_ct / kPglMaxDeltalistLenDivisor;
    if (patch_10_ct < max_deltalist_entry_ct) {
      // can't actually fail for now, but cannot ignore overflow check if
      // sizeof(AlleleCode) > 1
      if (unlikely(PwcAppendDeltalist(patch_10_set, patch_10_ct, pwcp, &vrec_len))) {
        return 1;
      }
      format_byte |= 16;
    } else {
      const uintptr_t genovec_altxy_ctb = DivUp(genovec_altxy_ct, CHAR_BIT);
      if (unlikely(CheckedVrecLenIncr(genovec_altxy_ctb, &vrec_len))) {
        return 1;
      }
      CopyGenomatchSubset(patch_10_set, genovec, kMaskAAAA, 0, genovec_altxy_ct, R_CAST(uintptr_t*, pwcp->fwrite_bufp));
      pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[genovec_altxy_ctb]);
    }
    if (allele_ct <= 17) {
      const AlleleCode* patch_10_vals_iter = patch_10_vals;
      uintptr_t* write_alias = R_CAST(uintptr_t*, pwcp->fwrite_bufp);
      uint32_t bytes_to_write;
      if (allele_ct == 3) {
        // 1 bit per entry, clear = alt1/alt2, set = alt2/alt2
        bytes_to_write = DivUp(patch_10_ct, 8);
        if (unlikely(CheckedVrecLenIncr(bytes_to_write, &vrec_len))) {
          return 1;
        }
#ifdef __LP64__
        const uint32_t fullvec_ct = patch_10_ct / (kBytesPerVec / 2);
        Vec16thUint* write_aliasv = R_CAST(Vec16thUint*, write_alias);
        if (fullvec_ct) {
#  if defined(USE_SSE42) && !defined(USE_AVX2)
          // SSE4.2: _mm_shuffle_epi8() to gather even bytes, parallel
          //         equality-to-2 check, movemask
          // (+126 or <<6 followed by movemask also works)
          const VecUc gather_even = vecuc_setr8(0, 2, 4, 6, 8, 10, 12, 14,
                                                -1, -1, -1, -1, -1, -1, -1, -1);
          for (uint32_t vec_idx = 0; vec_idx != fullvec_ct; ++vec_idx) {
            const VecUc cur_vec = vecuc_loadu(patch_10_vals_iter);
            patch_10_vals_iter = &(patch_10_vals_iter[kBytesPerVec]);
            const VecUc even_vec = vecuc_shuffle8(cur_vec, gather_even);
            const uint32_t cur_u32 = vecuc_movemask(R_CAST(VecUc, vecw_slli(R_CAST(VecW, even_vec), 6)));
            write_aliasv[vec_idx] = cur_u32;
          }
#  else
          // SSE2/AVX2: parallel equality-to-{2, 2} check, movemask,
          //            PackWordToHalfword
          // (uint16_t +126 followed by movemask also works, and avoids an
          // extra masking step in SSE2 case)
          // (<<6 doesn't work here)
          const VecU16 all126 = vecu16_set1(126);
          for (uint32_t vec_idx = 0; vec_idx != fullvec_ct; ++vec_idx) {
            VecU16 cur_vec = vecu16_loadu(patch_10_vals_iter);
            patch_10_vals_iter = &(patch_10_vals_iter[kBytesPerVec]);
            const uint32_t cur_u32 = vecu16_movemask(cur_vec + all126);
            write_aliasv[vec_idx] = PackVec8thUintTo16th(cur_u32);
          }
#  endif
        }
        const uint32_t remainder = patch_10_ct % (kBytesPerVec / 2);
        if (remainder) {
          uint32_t cur_u32 = 0;
          for (uint32_t uii = remainder - 1; ; --uii) {
            cur_u32 += patch_10_vals_iter[2 * uii];
            if (!uii) {
              break;
            }
            cur_u32 <<= 1;
          }
          // this effectively subtracts 1 from each entry
          cur_u32 += 1 - (1U << remainder);
          write_aliasv[fullvec_ct] = cur_u32;
          patch_10_vals_iter = &(patch_10_vals_iter[remainder * 2]);
        }
#else  // !__LP64__
        const uint32_t word_ct_m1 = (patch_10_ct - 1) / kBitsPerWord;
        uint32_t loop_len = kBitsPerWord;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= word_ct_m1) {
            if (widx > word_ct_m1) {
              break;
            }
            loop_len = ModNz(patch_10_ct, kBitsPerWord);
          }
          uintptr_t cur_word = 0;
          for (uint32_t uii = loop_len - 1; ; --uii) {
            cur_word |= patch_10_vals_iter[2 * uii] - 1;
            if (!uii) {
              break;
            }
            cur_word <<= 1;
          }
          write_alias[widx] = cur_word;
          patch_10_vals_iter = &(patch_10_vals_iter[loop_len * 2]);
        }
#endif  // !__LP64__
      } else if (allele_ct <= 5) {
        // 2 bits per half-entry
        bytes_to_write = DivUp(patch_10_ct, 2);
        if (unlikely(CheckedVrecLenIncr(bytes_to_write, &vrec_len))) {
          return 1;
        }
        // this may be worse than simple loop, check this later
        const uint32_t word_ct_m1 = (patch_10_ct - 1) / kBitsPerWordD4;
        uint32_t loop_len = kBitsPerWordD2;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= word_ct_m1) {
            if (widx > word_ct_m1) {
              break;
            }
            loop_len = 2 * ModNz(patch_10_ct, kBitsPerWordD4);
          }
          uintptr_t cur_word = 0;
          for (uint32_t uii = loop_len - 1; ; --uii) {
            // todo: same as 2-bit patch_01 case, except starting with
            // parallel-subtract 1 instead of 2
            cur_word |= patch_10_vals_iter[uii] - 1;
            if (!uii) {
              break;
            }
            cur_word <<= 2;
          }
          write_alias[widx] = cur_word;
          patch_10_vals_iter = &(patch_10_vals_iter[loop_len]);
        }
      } else {
        // 4 bits per half-entry
        bytes_to_write = patch_10_ct;
        if (unlikely(CheckedVrecLenIncr(bytes_to_write, &vrec_len))) {
          return 1;
        }
        const uint32_t word_ct_m1 = (patch_10_ct - 1) / kBytesPerWord;
        uint32_t loop_len = kBitsPerWordD4;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= word_ct_m1) {
            if (widx > word_ct_m1) {
              break;
            }
            loop_len = 2 * ModNz(patch_10_ct, kBytesPerWord);
          }
          uintptr_t cur_word = 0;
          for (uint32_t uii = loop_len - 1; ; --uii) {
            cur_word |= patch_10_vals_iter[uii] - 1;
            if (!uii) {
              break;
            }
            cur_word <<= 4;
          }
          write_alias[widx] = cur_word;
          patch_10_vals_iter = &(patch_10_vals_iter[loop_len]);
        }
      }
      pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[bytes_to_write]);
    } else {
      // 1 byte per half-entry
      // need 2- and 3-byte cases for larger AlleleCode
      unsigned char* payload = pwcp->fwrite_bufp;
      if (unlikely(CheckedVrecLenIncr(patch_10_ct * (2 * k1LU), &vrec_len))) {
        return 1;
      }
      const uint32_t patch_10_ct_x2 = patch_10_ct * 2;
      // hopefully the compiler automatically vectorizes this when
      // sizeof(AlleleCode) == 1
      for (uint32_t uii = 0; uii != patch_10_ct_x2; ++uii) {
        payload[uii] = patch_10_vals[uii] - 1;
      }
      pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[patch_10_ct * 2]);
    }
    if (genovec_hets_ptr) {
      // This is somewhat redundant with the earlier steps, but let's get
      // everything working before doing more optimization.
      const uint32_t sample_ct = pwcp->sample_ct;
      const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
      const Halfword* patch_10_set_alias = R_CAST(const Halfword*, patch_10_set);
      const AlleleCode* patch_10_vals_iter = patch_10_vals;
      const AlleleCode* patch_10_vals_end = &(patch_10_vals[patch_10_ct * 2]);
      uintptr_t* genovec_hets = nullptr;
      uint32_t het_ct = genovec_het_ct;
      // todo: try inverting this loop, at least in 64-bit case.  perform
      // movemask-based scan for heterozygous patch_10_vals entries, then
      // advance widx/detect_10 as needed.
      for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
        uint32_t patch_10_hw = patch_10_set_alias[widx];
        if (patch_10_hw) {
          do {
            const AlleleCode val1 = *patch_10_vals_iter++;
            const AlleleCode val2 = *patch_10_vals_iter++;
            const uintptr_t lowbit = patch_10_hw & (-patch_10_hw);
            if (val1 != val2) {
              if (!genovec_hets) {
                genovec_hets = pwcp->genovec_hets_buf;
                memcpy(genovec_hets, genovec, sample_ctl2 * sizeof(intptr_t));
              }
              // clear high bit, set low bit
              genovec_hets[widx] ^= lowbit * lowbit * 3;
              ++het_ct;
            }
            patch_10_hw ^= lowbit;
          } while (patch_10_hw);
          if (patch_10_vals_iter == patch_10_vals_end) {
            break;
          }
        }
      }
      *het_ct_ptr = het_ct;
      if (genovec_hets) {
        *genovec_hets_ptr = genovec_hets;
      } else {
        *genovec_hets_ptr = genovec;
      }
    }
  }
  *format_byte_ptr = format_byte;
  *vrtype_ptr |= 8;
  *vrec_len_ptr = vrec_len;
  return 0;
}

BoolErr PwcAppendMultiallelicSparse(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t patch_01_ct, uint32_t patch_10_ct, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  uint32_t vrec_len;
  if (unlikely(PwcAppendMultiallelicMain(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, patch_01_ct, patch_10_ct, vidx, pwcp, nullptr, nullptr, &vrtype, &vrec_len))) {
    return 1;
  }
  pwcp->vidx += 1;
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  SubU32Store(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]));
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
  return 0;
}

void PglMultiallelicDenseToSparse(const AlleleCode* __restrict wide_codes, uint32_t sample_ct, uintptr_t* __restrict genovec, uintptr_t* __restrict patch_01_set, AlleleCode* __restrict patch_01_vals, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals, uint32_t* __restrict patch_01_ct_ptr, uint32_t* __restrict patch_10_ct_ptr) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const AlleleCode* wide_codes_iter = wide_codes;
  Halfword* patch_01_set_alias = R_CAST(Halfword*, patch_01_set);
  Halfword* patch_10_set_alias = R_CAST(Halfword*, patch_10_set);
  AlleleCode* patch_01_vals_iter = patch_01_vals;
  AlleleCode* patch_10_vals_iter = patch_10_vals;
  uint32_t loop_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        *patch_01_ct_ptr = patch_01_vals_iter - patch_01_vals;
        *patch_10_ct_ptr = S_CAST(uintptr_t, patch_10_vals_iter - patch_10_vals) / 2;
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = 0;
    uint32_t patch_01_hw = 0;
    uint32_t patch_10_hw = 0;
    // ascending loop may be worse than descending here
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
      // todo: try making this movemask-based in 64-bit case
      // parallel checks for equality-to-255, equality-to-0, inequality-to-1
      // rarealts are inequality-to-1 and (not 0 or 255); can use ctzu32() and
      // &= x - 1 to iterate
      // hom-ref iff double-equality to 0
      // het-ref iff low equality-to-0 bit set, high clear
      // missing iff equality-to-255 (better be double)
      // altx-alty otherwise
      const AlleleCode first_code = *wide_codes_iter++;
      const AlleleCode second_code = *wide_codes_iter++;
      uintptr_t cur_geno;
      if (first_code) {
        if (first_code == kMissingAlleleCode) {
          cur_geno = 3;
        } else {
          cur_geno = 2;
          if (second_code > 1) {
            patch_10_hw |= 1U << sample_idx_lowbits;
            *patch_10_vals_iter++ = first_code;
            *patch_10_vals_iter++ = second_code;
          }
        }
      } else if (!second_code) {
        continue;
      } else {
        cur_geno = 1;
        if (second_code > 1) {
          patch_01_hw |= 1U << sample_idx_lowbits;
          *patch_01_vals_iter++ = second_code;
        }
      }
      geno_word |= cur_geno << (2 * sample_idx_lowbits);
    }
    genovec[widx] = geno_word;
    patch_01_set_alias[widx] = patch_01_hw;
    patch_10_set_alias[widx] = patch_10_hw;
  }
}

static const uint16_t kHcToAlleleCodes[1024] = QUAD_TABLE256(0, 0x100, 0x101, 0xffff);

static_assert(sizeof(AlleleCode) == 1, "PglMultiallelicSparseToDense() needs to be updated.");
void PglMultiallelicSparseToDense(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, const AlleleCode* __restrict remap, uint32_t sample_ct, uint32_t patch_01_ct, uint32_t patch_10_ct, uintptr_t* __restrict flipped, AlleleCode* __restrict wide_codes) {
  if (flipped) {
    ZeroWArr(BitCtToWordCt(sample_ct), flipped);
  }
  if ((!remap) || ((!remap[0]) && (remap[1] == 1))) {
    GenoarrLookup256x2bx4(genovec, kHcToAlleleCodes, sample_ct, wide_codes);
    if (!remap) {
      if (patch_01_ct) {
        uintptr_t sample_idx_base = 0;
        uintptr_t cur_bits = patch_01_set[0];
        AlleleCode* wide_codes1 = &(wide_codes[1]);
        for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
          const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
          wide_codes1[2 * sample_idx] = patch_01_vals[uii];
        }
      }
      if (patch_10_ct) {
        uintptr_t sample_idx_base = 0;
        uintptr_t cur_bits = patch_10_set[0];
        const DoubleAlleleCode* patch_10_vals_alias = R_CAST(const DoubleAlleleCode*, patch_10_vals);
        DoubleAlleleCode* wide_codes_alias = R_CAST(DoubleAlleleCode*, wide_codes);
        for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
          const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
          wide_codes_alias[sample_idx] = patch_10_vals_alias[uii];
        }
      }
      return;
    }
  } else {
    // 0 -> (remap[0], remap[0])
    // 1 -> (remap[0], remap[1])
    // 2 -> (remap[1], remap[1])
    // 3 -> (255, 255)
    // 256x4 lookup table usually too expensive to reinitialize for every
    // single variant.  16x2 should be a good compromise?
    // todo: try Expand1bitTo8 strategy in SSE4.2+ case, see if that's fast
    // enough to remove 0/1 special case
    const uint32_t remap0 = remap[0];
    const uint32_t remap1 = remap[1];
    uint32_t table4[4];
    table4[0] = remap0 * 257;
    table4[1] = remap0 + remap1 * 256;
    table4[2] = remap1 * 257;
    table4[3] = 0xffff;
    if (remap1 < remap0) {
      table4[1] = remap1 + remap0 * 256;
      if (flipped) {
        // See GenovecToMissingnessUnsafe().
        const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
        Halfword* flipped_alias = R_CAST(Halfword*, flipped);
        for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
          const uintptr_t cur_geno_word = genovec[widx];
          flipped_alias[widx] = PackWordToHalfwordMask5555(cur_geno_word & (~(cur_geno_word >> 1)));
        }
      }
    }
    // could have an InitLookup16x2bx2 function for this
    uint32_t table16[16];
    for (uint32_t uii = 0; uii != 4; ++uii) {
      const uint32_t high_bits = table4[uii] << 16;
      uint32_t* table_row = &(table16[uii * 4]);
      for (uint32_t ujj = 0; ujj != 4; ++ujj) {
        table_row[ujj] = high_bits | table4[ujj];
      }
    }
    const uint32_t sample_ctd4 = sample_ct / 4;
    const unsigned char* genovec_alias = R_CAST(const unsigned char*, genovec);
    uint32_t* wide_codes_alias_iter = R_CAST(uint32_t*, wide_codes);
    for (uint32_t byte_idx = 0; byte_idx != sample_ctd4; ++byte_idx) {
      const uint32_t geno_byte = genovec_alias[byte_idx];
      *wide_codes_alias_iter++ = table16[geno_byte & 15];
      *wide_codes_alias_iter++ = table16[geno_byte >> 4];
    }
    const uint32_t remainder = sample_ct % 4;
    if (remainder) {
      uint32_t geno_byte = genovec_alias[sample_ctd4];
      if (remainder & 2) {
        *wide_codes_alias_iter++ = table16[geno_byte & 15];
        geno_byte = geno_byte >> 4;
      }
      if (remainder & 1) {
        uint16_t* last_code_pair_ptr = R_CAST(uint16_t*, wide_codes_alias_iter);
        // guess it's easy enough to defend against dirty trailing bits for now
        *last_code_pair_ptr = table4[geno_byte & 3];
      }
    }
  }
  if (patch_01_vals) {
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = patch_01_set[0];
    const AlleleCode remap0 = remap[0];
    const AlleleCode remap1 = remap[1];
    // TODO: need special flipped handling for remap1 < remap0
    if ((!remap0) || ((!flipped) && (remap0 == 1) && (!remap1))) {
      // no flips possible
      AlleleCode* wide_codes1 = &(wide_codes[1]);
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
        wide_codes1[2 * sample_idx] = remap[patch_01_vals[uii]];
      }
    } else if (!flipped) {
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac = remap[patch_01_vals[uii]];
        if (ac > remap0) {
          wide_codes[2 * sample_idx + 1] = ac;
        } else {
          wide_codes[2 * sample_idx] = ac;
          wide_codes[2 * sample_idx + 1] = remap0;
        }
      }
    } else if (remap1 >= remap0) {
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac = remap[patch_01_vals[uii]];
        if (ac > remap0) {
          wide_codes[2 * sample_idx + 1] = ac;
        } else {
          wide_codes[2 * sample_idx] = ac;
          wide_codes[2 * sample_idx + 1] = remap0;
          SetBit(sample_idx, flipped);
        }
      }
    } else {
      for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac = remap[patch_01_vals[uii]];
        if (ac > remap0) {
          wide_codes[2 * sample_idx + 1] = ac;
          ClearBit(sample_idx, flipped);
        } else {
          wide_codes[2 * sample_idx] = ac;
          wide_codes[2 * sample_idx + 1] = remap0;
        }
      }
    }
  }
  if (patch_10_ct) {
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = patch_10_set[0];
    if (!flipped) {
      for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac0 = remap[patch_10_vals[2 * uii]];
        const AlleleCode ac1 = remap[patch_10_vals[2 * uii + 1]];
        if (ac0 <= ac1) {
          wide_codes[2 * sample_idx] = ac0;
          wide_codes[2 * sample_idx + 1] = ac1;
        } else {
          wide_codes[2 * sample_idx] = ac1;
          wide_codes[2 * sample_idx + 1] = ac0;
        }
      }
    } else {
      for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
        const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
        const AlleleCode ac0 = remap[patch_10_vals[2 * uii]];
        const AlleleCode ac1 = remap[patch_10_vals[2 * uii + 1]];
        if (ac0 <= ac1) {
          wide_codes[2 * sample_idx] = ac0;
          wide_codes[2 * sample_idx + 1] = ac1;
        } else {
          wide_codes[2 * sample_idx] = ac1;
          wide_codes[2 * sample_idx + 1] = ac0;
          SetBit(sample_idx, flipped);
        }
      }
    }
  }
}

static_assert(sizeof(AlleleCode) == 1, "PglMultiallelicSparseToDenseMiss() needs to be updated.");
void PglMultiallelicSparseToDenseMiss(const PgenVariant* pgvp, uint32_t sample_ct, AlleleCode* __restrict wide_codes) {
  GenoarrLookup256x2bx4(pgvp->genovec, kHcToAlleleCodes, sample_ct, wide_codes);
  const uint32_t patch_01_ct = pgvp->patch_01_ct;
  if (patch_01_ct) {
    const uintptr_t* patch_01_set = pgvp->patch_01_set;
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = patch_01_set[0];
    const AlleleCode* patch_01_vals = pgvp->patch_01_vals;
    AlleleCode* wide_codes1 = &(wide_codes[1]);
    for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
      wide_codes1[2 * sample_idx] = patch_01_vals[uii];
    }
  }
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  if (patch_10_ct) {
    const uintptr_t* patch_10_set = pgvp->patch_10_set;
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = patch_10_set[0];
    const DoubleAlleleCode* patch_10_vals_alias = R_CAST(const DoubleAlleleCode*, pgvp->patch_10_vals);
    DoubleAlleleCode* wide_codes_alias = R_CAST(DoubleAlleleCode*, wide_codes);
    for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
      wide_codes_alias[sample_idx] = patch_10_vals_alias[uii];
    }
  }
}

// tolerates extraneous phaseinfo bits
BoolErr AppendHphase(const uintptr_t* __restrict genovec_hets, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, uint32_t het_ct, uint32_t phasepresent_ct, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  assert(phasepresent_ct);
  const uint32_t sample_ct = pwcp->sample_ct;
  *vrtype_ptr += 16;
  const uint32_t het_ctp1_8 = 1 + (het_ct / CHAR_BIT);
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
    if (unlikely(CheckedVrecLenIncr(het_ctp1_8, vrec_len_ptr))) {
      return 1;
    }
    CopyGenomatchSubset(phaseinfo, genovec_hets, kMask5555, 1, het_ct, fwrite_bufp_alias);
    fwrite_bufp_final = &(pwcp->fwrite_bufp[het_ctp1_8]);
  } else {
    // this is a minor variant of ExpandThenSubsetBytearr()
    if (unlikely(CheckedVrecLenIncr(het_ctp1_8 + DivUp(phasepresent_ct, 8), vrec_len_ptr))) {
      return 1;
    }
    const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
    uintptr_t* phaseinfo_tmp = pwcp->genovec_invert_buf;
    uintptr_t* phaseinfo_tmp_iter = phaseinfo_tmp;
    uint32_t phasepresent_write_idx_lowbits = 1;
    phaseinfo_write_idx_lowbits = 0;
    uintptr_t phasepresent_write_word = 1;  // first bit set
    for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
      const uintptr_t geno_word = genovec_hets[widx];
      uintptr_t geno_hets = Word01(geno_word);
      if (geno_hets) {
        const uint32_t phasepresent_halfword = R_CAST(const Halfword*, phasepresent)[widx];
        if (phasepresent_halfword) {
          const uint32_t phaseinfo_halfword = R_CAST(const Halfword*, phaseinfo)[widx];
          do {
            const uint32_t sample_idx_lowbits = ctzw(geno_hets) / 2;
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
      // er, safe to write the entire word...
      SubwordStoreMov(phasepresent_write_word, cur_byte_ct, &fwrite_bufp_final);
    }
    fwrite_bufp_final = memcpyua(fwrite_bufp_final, phaseinfo_tmp, sizeof(intptr_t) * (phaseinfo_tmp_iter - phaseinfo_tmp));
    if (phaseinfo_write_idx_lowbits) {
      const uint32_t cur_byte_ct = DivUp(phaseinfo_write_idx_lowbits, CHAR_BIT);
      SubwordStoreMov(phaseinfo_write_word, cur_byte_ct, &fwrite_bufp_final);
    }
    assert(S_CAST(uintptr_t, fwrite_bufp_final - pwcp->fwrite_bufp) == het_ctp1_8 + DivUp(phasepresent_ct, 8));
  }
  pwcp->fwrite_bufp = fwrite_bufp_final;
  return 0;
}

void PwcAppendBiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, PgenWriterCommon* pwcp) {
  // assumes phase_dosage_gflags is nonzero
  const uint32_t vidx = pwcp->vidx;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  uint32_t het_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, nullptr, vrtype_dest);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  const uint32_t phasepresent_ct = phasepresent? PopcountWords(phasepresent, sample_ctl) : het_ct;
  if (phasepresent_ct) {
    AppendHphase(genovec, phasepresent, phaseinfo, het_ct, phasepresent_ct, pwcp, vrtype_dest, &vrec_len);
  }
  SubU32Store(vrec_len, vrec_len_byte_ct, vrec_len_dest);
}

BoolErr PwcAppendMultiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, uint32_t patch_01_ct, uint32_t patch_10_ct, PgenWriterCommon* pwcp) {
  const uint32_t vidx = pwcp->vidx;
  const uintptr_t* genovec_hets;
  unsigned char vrtype;
  uint32_t het_ct;
  uint32_t vrec_len;
  if (unlikely(PwcAppendMultiallelicMain(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, patch_01_ct, patch_10_ct, vidx, pwcp, &genovec_hets, &het_ct, &vrtype, &vrec_len))) {
    return 1;
  }
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t phasepresent_ct = phasepresent? PopcountWords(phasepresent, sample_ctl) : het_ct;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  *vrtype_dest = vrtype;
  if (phasepresent_ct) {
    if (unlikely(AppendHphase(genovec_hets, phasepresent, phaseinfo, het_ct, phasepresent_ct, pwcp, vrtype_dest, &vrec_len))) {
      return 1;
    }
  }
  pwcp->vidx += 1;
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  SubU32Store(vrec_len, vrec_len_byte_ct, &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]));
  return 0;
}

// ok if dosage_present trailing bits set
BoolErr AppendDosage16(const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, uint32_t dphase_ct, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  const uint32_t sample_ct = pwcp->sample_ct;
  const uint32_t max_deltalist_entry_ct = sample_ct / kPglMaxDeltalistLenDivisor;
  if (dosage_ct < max_deltalist_entry_ct) {
    // case 1: store dosage IDs as deltalist.
    if (unlikely(PwcAppendDeltalist(dosage_present, dosage_ct, pwcp, vrec_len_ptr))) {
      return 1;
    }
    *vrtype_ptr += 0x20;
  } else if ((dosage_ct == sample_ct) && ((!dphase_ct) || (dphase_ct == sample_ct))) {
    // case 2: fixed-width, no need to store dosage IDs at all.
    // dosage_main permitted to have 65535 = missing
    // bugfix (9 Dec 2018): since this forces fixed-width dosage-phase storage,
    // also require dphase_ct == 0 or dphase_ct == sample_ct.
    *vrtype_ptr += 0x40;
  } else {
    // case 3: save dosage_present bitarray directly.
    const uint32_t sample_ctb = DivUp(sample_ct, CHAR_BIT);
    if (unlikely(CheckedVrecLenIncr(sample_ctb, vrec_len_ptr))) {
      return 1;
    }
    pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, dosage_present, sample_ctb);
    *vrtype_ptr += 0x60;
  }
  if (unlikely(CheckedVrecLenIncr(dosage_ct * sizeof(int16_t), vrec_len_ptr))) {
    return 1;
  }
  pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, dosage_main, dosage_ct * sizeof(int16_t));
  return 0;
}

// ok if dosage_present trailing bits set
BoolErr PwcAppendBiallelicGenovecDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, PgenWriterCommon* pwcp) {
  // safe to call this even when entire file has no phase/dosage info
  const uint32_t vidx = pwcp->vidx;
  unsigned char vrtype;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, nullptr, nullptr, &vrtype);
  const uintptr_t vrec_len_byte_ct = pwcp->vrec_len_byte_ct;
  pwcp->vidx += 1;
  unsigned char* vrec_len_dest = &(pwcp->vrec_len_buf[vidx * vrec_len_byte_ct]);
  if (dosage_ct) {
    if (unlikely(AppendDosage16(dosage_present, dosage_main, dosage_ct, 0, pwcp, &vrtype, &vrec_len))) {
      return 1;
    }
  }
  SubU32Store(vrec_len, vrec_len_byte_ct, vrec_len_dest);
  if (!pwcp->phase_dosage_gflags) {
    pwcp->vrtype_buf[vidx / kBitsPerWordD4] |= S_CAST(uintptr_t, vrtype) << (4 * (vidx % kBitsPerWordD4));
  } else {
    R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx] = vrtype;
  }
  return 0;
}

// NOT ok if phasepresent trailing bits set, but ok for dosage_present
BoolErr PwcAppendBiallelicGenovecHphaseDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, PgenWriterCommon* pwcp) {
  // assumes there is phase and/or dosage data in output file, otherwise
  // vrtype_dest needs to be replaced

  // this mostly overlaps with PwcAppendBiallelicGenovecHphase(); probably
  // get rid of the latter
  const uint32_t vidx = pwcp->vidx;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  uint32_t het_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, nullptr, vrtype_dest);
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
    if (unlikely(AppendDosage16(dosage_present, dosage_main, dosage_ct, 0, pwcp, vrtype_dest, &vrec_len))) {
      return 1;
    }
  }
  SubU32Store(vrec_len, vrec_len_byte_ct, vrec_len_dest);
  return 0;
}

BoolErr AppendDphase16(const uintptr_t* __restrict dosage_present, const uintptr_t* __restrict dphase_present, const int16_t* dphase_delta, uint32_t dosage_ct, uint32_t dphase_ct, PgenWriterCommon* pwcp, unsigned char* vrtype_ptr, uint32_t* vrec_len_ptr) {
  assert(dphase_ct);
  *vrtype_ptr += 0x80;
  if (dphase_ct != pwcp->sample_ct) {
    // bugfix (9 Dec 2018): forgot to separate fixed-width phased-dosage case
    // bugfix (12 Feb 2019): double-added dphase_ct * sizeof(int16_t)
    const uint32_t dphase_present_byte_ct = DivUp(dosage_ct, CHAR_BIT);
    if (unlikely(CheckedVrecLenIncr(dphase_present_byte_ct, vrec_len_ptr))) {
      return 1;
    }
    CopyBitarrSubset(dphase_present, dosage_present, dosage_ct, R_CAST(uintptr_t*, pwcp->fwrite_bufp));
    pwcp->fwrite_bufp = &(pwcp->fwrite_bufp[dphase_present_byte_ct]);
  }
  // bugfix (2 Jan 2019): forgot to update vrec_len here
  if (unlikely(CheckedVrecLenIncr(dphase_ct * sizeof(int16_t), vrec_len_ptr))) {
    return 1;
  }
  pwcp->fwrite_bufp = memcpyua(pwcp->fwrite_bufp, dphase_delta, dphase_ct * sizeof(int16_t));
  return 0;
}

BoolErr PwcAppendBiallelicGenovecDphase16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uintptr_t* __restrict dphase_present, const uint16_t* dosage_main, const int16_t* dphase_delta, uint32_t dosage_ct, uint32_t dphase_ct, PgenWriterCommon* pwcp) {
  // assumes there is phase and/or dosage data in output file, otherwise
  // vrtype_dest needs to be replaced
  const uint32_t vidx = pwcp->vidx;
  unsigned char* vrtype_dest = &(R_CAST(unsigned char*, pwcp->vrtype_buf)[vidx]);
  uint32_t het_ct;
  uint32_t vrec_len = PwcAppendBiallelicGenovecMain(genovec, vidx, pwcp, &het_ct, nullptr, vrtype_dest);
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
    if (unlikely(AppendDosage16(dosage_present, dosage_main, dosage_ct, dphase_ct, pwcp, vrtype_dest, &vrec_len))) {
      return 1;
    }
    if (dphase_ct) {
      if (unlikely(AppendDphase16(dosage_present, dphase_present, dphase_delta, dosage_ct, dphase_ct, pwcp, vrtype_dest, &vrec_len))) {
        return 1;
      }
    }
  }
  SubU32Store(vrec_len, vrec_len_byte_ct, vrec_len_dest);
  return 0;
}

PglErr PwcFinish(PgenWriterCommon* pwcp, FILE** pgen_outfile_ptr) {
  const uint32_t variant_ct = pwcp->variant_ct;
  assert(pwcp->vidx == variant_ct);
  FILE* pgen_outfile = *pgen_outfile_ptr;
  if (unlikely(fseeko(pgen_outfile, 12, SEEK_SET))) {
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
  for (; ; vrec_len_buf_iter = &(vrec_len_buf_iter[vrec_iter_incr])) {
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
    if (unlikely(fwrite_checked(vrec_len_buf_iter, vrec_iter_incr, pgen_outfile))) {
      return kPglRetWriteFail;
    }

    // 4b(iii): alt allele counts
    // not yet supported

    // 4b(iv): explicit nonref flags
    if (explicit_nonref_flags) {
      if (unlikely(fwrite_checked(explicit_nonref_flags_iter, nonref_flags_write_byte_ct, pgen_outfile))) {
        return kPglRetWriteFail;
      }
      explicit_nonref_flags_iter = &(explicit_nonref_flags_iter[kPglVblockSize / kBitsPerWord]);
    }
  }
}

PglErr SpgwFinish(STPgenWriter* spgwp) {
  if (unlikely(fwrite_checked(spgwp->pwc.fwrite_buf, spgwp->pwc.fwrite_bufp - spgwp->pwc.fwrite_buf, spgwp->pgen_outfile))) {
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
  for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
    vblock_fpos[(vidx / kPglVblockSize) + tidx] = cur_vblock_fpos;
    PgenWriterCommon* cur_pwcp = mpgwp->pwcs[tidx];
    uintptr_t cur_vblock_byte_ct = cur_pwcp->fwrite_bufp - cur_pwcp->fwrite_buf;
    if (unlikely(fwrite_checked(cur_pwcp->fwrite_buf, cur_vblock_byte_ct, pgen_outfile))) {
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
