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

#include "pgenlib_misc.h"

#include <limits.h>

#ifdef __cplusplus
namespace plink2 {
#endif

#ifndef NDEBUG
char* g_pgl_errbuf = nullptr;
char* g_pgl_errbuf_write_iter = nullptr;
char* g_pgl_errbuf_end = nullptr;
#endif

#ifdef USE_AVX2
void CopyNyparrNonemptySubset(const uintptr_t* __restrict raw_nyparr, const uintptr_t* __restrict subset_mask, uint32_t raw_nyparr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_nyparr) {
  if (subset_entry_ct == raw_nyparr_entry_ct) {
    memcpy(output_nyparr, raw_nyparr, DivUp(subset_entry_ct, kBitsPerWordD2) * sizeof(intptr_t));
    ZeroTrailingNyps(subset_entry_ct, output_nyparr);
    return;
  }
  assert(subset_entry_ct);
  uintptr_t cur_output_word = 0;

  uintptr_t* output_nyparr_iter = output_nyparr;

  uintptr_t* output_nyparr_last = &(output_nyparr[subset_entry_ct / kBitsPerWordD2]);
  const uint32_t word_write_shift_end = 2 * (subset_entry_ct % kBitsPerWordD2);
  uint32_t word_write_shift = 0;
  for (uint32_t subset_mask_widx = 0; ; ++subset_mask_widx) {
    const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
    if (cur_include_word) {
      uint32_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
      for (uint32_t wordhalf_idx = 0; ; ++wordhalf_idx) {
        if (cur_include_halfword) {
          uintptr_t extracted_bits = raw_nyparr[subset_mask_widx * 2 + wordhalf_idx];
          uint32_t set_bit_ct = kBitsPerWord;
          if (cur_include_halfword != UINT32_MAX) {
            const uintptr_t pext_mask = 3 * UnpackHalfwordToWord(cur_include_halfword);
            extracted_bits = _pext_u64(extracted_bits, pext_mask);
            set_bit_ct = PopcountWord(pext_mask);
          }
          cur_output_word |= extracted_bits << word_write_shift;
          word_write_shift += set_bit_ct;
          if (word_write_shift >= kBitsPerWord) {
            *output_nyparr_iter++ = cur_output_word;
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
      if (output_nyparr_iter == output_nyparr_last) {
        if (word_write_shift == word_write_shift_end) {
          if (word_write_shift_end) {
            *output_nyparr_last = cur_output_word;
          }
          return;
        }
      }
    }
  }
}

// bit_idx_start assumed to be < kBitsPerWord
void CopyGenomatchSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict genoarr, uintptr_t match_word, uint32_t write_bit_idx_start, uint32_t bit_ct, void* __restrict output) {
  const uint32_t bit_idx_end = bit_ct + write_bit_idx_start;
  const uint32_t bit_idx_end_lowbits = bit_idx_end % kBitsPerWord;
  const Halfword* raw_bitarr_alias = DowncastKWToHW(raw_bitarr);
  unsigned char* output_biter = S_CAST(unsigned char*, output);
  unsigned char* output_last = &(output_biter[(bit_idx_end / kBitsPerWord) * kBytesPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = write_bit_idx_start;
  while ((output_biter != output_last) || (write_idx_lowbits != bit_idx_end_lowbits)) {
    uintptr_t cur_mask_word;
    // sparse genoarr optimization
    // guaranteed to terminate since there's at least one more set bit
    do {
      // todo: try reading two genoarr words at a time.  would need to be very
      // careful with the possible trailing word, though.
      // more important to optimize this function now that regular phased-call
      // handling code is using it.
      cur_mask_word = genoarr[++read_widx] ^ match_word;
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
      AppendW(cur_output_word, &output_biter);
      // ...and these are the bits that fell off
      // impossible for write_idx_lowbits to be zero here
      cur_output_word = extracted_bits >> (kBitsPerWord - write_idx_lowbits);
    }
    write_idx_lowbits = new_write_idx_lowbits % kBitsPerWord;
  }
  if (write_idx_lowbits) {
    CopyToUnalignedW(output_biter, &cur_output_word);
  }
}

// Variant of ExpandBytearr() which is based off a target 2-bit value instead
// of single expand_mask bits.  expand_size must be the number of instances of
// the target value in genovec.
void ExpandBytearrFromGenoarr(const void* __restrict compact_bitarr, const uintptr_t* __restrict genoarr, uintptr_t match_word, uint32_t genoword_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t leading_byte_ct = 1 + (expand_sizex_m1 % kBitsPerWord) / CHAR_BIT;
  const uint32_t genoword_ct_m1 = genoword_ct - 1;
  uintptr_t compact_word = SubwordLoad(compact_bitarr, leading_byte_ct) >> read_start_bit;
  const unsigned char* compact_bitarr_biter = &(S_CAST(const unsigned char*, compact_bitarr)[leading_byte_ct]);
  uint32_t compact_idx_lowbits = read_start_bit + CHAR_BIT * (sizeof(intptr_t) - leading_byte_ct);
  for (uint32_t widx = 0; ; widx += 2) {
    uintptr_t mask_word;
    if (widx >= genoword_ct_m1) {
      if (widx > genoword_ct_m1) {
        return;
      }
      mask_word = 0;
    } else {
      const uintptr_t geno_word1 = genoarr[widx + 1] ^ match_word;
      mask_word = PackWordToHalfwordMask5555(~(geno_word1 | (geno_word1 >> 1)));
      mask_word = mask_word << 32;
    }
    const uintptr_t geno_word0 = genoarr[widx] ^ match_word;
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
        uintptr_t next_compact_word;
        CopyFromUnalignedIncrW(&next_compact_word, &compact_bitarr_biter);
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
void CopyNyparrNonemptySubset(const uintptr_t* __restrict raw_nyparr, const uintptr_t* __restrict subset_mask, uint32_t raw_nyparr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_nyparr) {
  if (subset_entry_ct == raw_nyparr_entry_ct) {
    // subset_mask may be nullptr in this case
    memcpy(output_nyparr, raw_nyparr, DivUp(subset_entry_ct, kBitsPerWordD2) * sizeof(intptr_t));
    ZeroTrailingNyps(subset_entry_ct, output_nyparr);
    return;
  }
  assert(subset_entry_ct);
  assert(raw_nyparr_entry_ct >= subset_entry_ct);
  uintptr_t cur_output_word = 0;

  uintptr_t* output_nyparr_iter = output_nyparr;

  uintptr_t* output_nyparr_last = &(output_nyparr[subset_entry_ct / kBitsPerWordD2]);
  const uint32_t word_write_halfshift_end = subset_entry_ct % kBitsPerWordD2;
  uint32_t word_write_halfshift = 0;
  // if <= 2/3-filled, use sparse copy algorithm
  // (tried CopyBitarrSubset() approach, that actually worsened things)
  if (subset_entry_ct * (3 * k1LU) <= raw_nyparr_entry_ct * (2 * k1LU)) {
    for (uint32_t subset_mask_widx = 0; ; ++subset_mask_widx) {
      const uintptr_t cur_include_word = subset_mask[subset_mask_widx];
      if (cur_include_word) {
        uint32_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
        for (uint32_t wordhalf_idx = 0; ; ++wordhalf_idx) {
          if (cur_include_halfword) {
            uintptr_t raw_nyparr_word = raw_nyparr[subset_mask_widx * 2 + wordhalf_idx];
            do {
              uint32_t rqa_idx_lowbits = ctzu32(cur_include_halfword);
              cur_output_word |= ((raw_nyparr_word >> (rqa_idx_lowbits * 2)) & 3) << (word_write_halfshift * 2);
              if (++word_write_halfshift == kBitsPerWordD2) {
                *output_nyparr_iter++ = cur_output_word;
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
        if (output_nyparr_iter == output_nyparr_last) {
          if (word_write_halfshift == word_write_halfshift_end) {
            if (word_write_halfshift_end) {
              *output_nyparr_last = cur_output_word;
            }
            return;
          }
        }
      }
    }
  }
  // blocked copy
  const uintptr_t* raw_nyparr_iter = raw_nyparr;
  for (; ; ++subset_mask) {
    const uintptr_t cur_include_word = *subset_mask;
    uintptr_t cur_include_halfword = S_CAST(Halfword, cur_include_word);
    for (uint32_t wordhalf_idx = 0; ; ++wordhalf_idx) {
      uintptr_t raw_nyparr_word = *raw_nyparr_iter++;
      while (cur_include_halfword) {
        const uint32_t rqa_idx_lowbits = ctzw(cur_include_halfword);

        // TAOCP, 7.1.3, (43).
        const uintptr_t bottom_block_remover = (cur_include_halfword | (cur_include_halfword - 1)) + 1;

        const uintptr_t raw_nyparr_curblock_unmasked = raw_nyparr_word >> (rqa_idx_lowbits * 2);
        const uint32_t rqa_block_len = ctzw(bottom_block_remover) - rqa_idx_lowbits;
        const uint32_t block_len_limit = kBitsPerWordD2 - word_write_halfshift;
        cur_output_word |= raw_nyparr_curblock_unmasked << (2 * word_write_halfshift);
        if (rqa_block_len < block_len_limit) {
          word_write_halfshift += rqa_block_len;
          cur_output_word = bzhi(cur_output_word, word_write_halfshift * 2);
        } else {
          // no need to mask, extra bits vanish off the high end
          *output_nyparr_iter++ = cur_output_word;
          word_write_halfshift = rqa_block_len - block_len_limit;
          if (word_write_halfshift) {
            cur_output_word = bzhi(raw_nyparr_curblock_unmasked >> (2 * block_len_limit), 2 * word_write_halfshift);
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
    if (output_nyparr_iter == output_nyparr_last) {
      if (word_write_halfshift == word_write_halfshift_end) {
        if (word_write_halfshift_end) {
          *output_nyparr_last = cur_output_word;
        }
        return;
      }
    }
  }
}

void CopyGenomatchSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict genovec, uintptr_t match_word, uint32_t write_bit_idx_start, uint32_t bit_ct, void* __restrict output) {
  const uint32_t bit_idx_end = write_bit_idx_start + bit_ct;
  const uint32_t bit_idx_end_lowbits = bit_idx_end % kBitsPerWord;
  const Halfword* raw_bitarr_alias = DowncastKWToHW(raw_bitarr);
  unsigned char* output_biter = S_CAST(unsigned char*, output);
  unsigned char* output_last = &(output_biter[(bit_idx_end / kBitsPerWord) * kBytesPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = UINT32_MAX;  // deliberate overflow
  uint32_t write_idx_lowbits = write_bit_idx_start;
  while ((output_biter != output_last) || (write_idx_lowbits != bit_idx_end_lowbits)) {
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
        AppendW(cur_output_word, &output_biter);
        cur_output_word = 0;
        write_idx_lowbits = 0;
      }
      geno_word &= geno_word - 1;
    } while (geno_word);
  }
  if (write_idx_lowbits) {
    CopyToUnalignedW(output_biter, &cur_output_word);
  }
}

void ExpandBytearrFromGenoarr(const void* __restrict compact_bitarr, const uintptr_t* __restrict genoarr, uintptr_t match_word, uint32_t genoword_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target) {
  Halfword* target_alias = DowncastWToHW(target);
  ZeroHwArr(RoundUpPow2(genoword_ct, 2), target_alias);
  const unsigned char* compact_bitarr_alias = S_CAST(const unsigned char*, compact_bitarr);
  const uint32_t expand_sizex_m1 = expand_size + read_start_bit - 1;
  const uint32_t compact_widx_last = expand_sizex_m1 / kBitsPerWord;
  uint32_t compact_idx_lowbits = read_start_bit;
  uint32_t loop_len = kBitsPerWord;
  uintptr_t write_hwidx = 0;
  uintptr_t genomatch_bits = genoarr[0] ^ match_word;
  genomatch_bits = (~(genomatch_bits | (genomatch_bits >> 1))) & kMask5555;
  for (uint32_t compact_widx = 0; ; ++compact_widx) {
    uintptr_t compact_word;
    if (compact_widx >= compact_widx_last) {
      if (compact_widx > compact_widx_last) {
        return;
      }
      loop_len = 1 + (expand_sizex_m1 % kBitsPerWord);
      // avoid possible segfault
      compact_word = SubwordLoad(&(compact_bitarr_alias[compact_widx * kBytesPerWord]), DivUp(loop_len, CHAR_BIT));
    } else {
      CopyFromUnalignedOffsetW(&compact_word, compact_bitarr_alias, compact_widx);
    }
    for (; compact_idx_lowbits != loop_len; ++compact_idx_lowbits) {
      while (!genomatch_bits) {
        genomatch_bits = genoarr[++write_hwidx] ^ match_word;
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

// Harley-Seal algorithm only works for bitarrays, not nyparrays, so don't
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
      VecW alt1_plus_bothset1 = vecw_and_notfirst(cur_geno_vword_low_lshifted1, cur_geno_vword1);
      VecW alt1_plus_bothset2 = vecw_and_notfirst(cur_geno_vword_low_lshifted2, cur_geno_vword2);

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
      alt1_plus_bothset1 = vecw_and_notfirst(cur_geno_vword_low_lshifted1, cur_geno_vword1);
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
      VecW nyp_both = m1 & (cur_vvec ^ vvec_rshift);
      VecW nyp_01 = nyp_both & cur_vvec;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      VecW vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_01 + (vvec_both & cur_vvec);
      nyp_both = nyp_both + vvec_both;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_01 + (vvec_both & cur_vvec);
      nyp_both = nyp_both + vvec_both;

      VecW nybble_01 = (nyp_01 & m2) + (vecw_srli(nyp_01, 2) & m2);
      VecW nybble_both = (nyp_both & m2) + (vecw_srli(nyp_both, 2) & m2);

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      nyp_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_both & cur_vvec;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_01 + (vvec_both & cur_vvec);
      nyp_both = nyp_both + vvec_both;

      cur_vvec = *geno_vvec_iter++;
      vvec_rshift = vecw_srli(cur_vvec, 1);
      vvec_both = m1 & (cur_vvec ^ vvec_rshift);
      nyp_01 = nyp_01 + (vvec_both & cur_vvec);
      nyp_both = nyp_both + vvec_both;

      nybble_01 = nybble_01 + (nyp_01 & m2) + (vecw_srli(nyp_01, 2) & m2);
      nybble_both = nybble_both + (nyp_both & m2) + (vecw_srli(nyp_both, 2) & m2);

      inner_acc_01 = inner_acc_01 + (nybble_01 & m4) + (vecw_srli(nybble_01, 4) & m4);
      inner_acc_both = inner_acc_both + (nybble_both & m4) + (vecw_srli(nybble_both, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_01 = acc_01 + vecw_bytesum(inner_acc_01, m0);
    acc_both = acc_both + vecw_bytesum(inner_acc_both, m0);
  }
}

void GenovecCount12Unsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* __restrict raw_01_ctp, uint32_t* __restrict raw_10_ctp) {
  // assumes trailing bits of last genovec word are zeroed out.
  // sample_ct == 0 ok.
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uint32_t raw_01_ct;
  uint32_t raw_both_ct;
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (6 * kWordsPerVec));
  assert(IsVecAligned(genovec));
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

void Count3FreqVec6(const void* genoarr, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 6));
  // Sets even_ct to the number of set low bits in the current block, odd_ct to
  // the number of set high bits, and bothset_ct by the number of 0b11s.
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const unsigned char* geno_vvec_biter = S_CAST(const unsigned char*, genoarr);
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
    const unsigned char* geno_vvec_stop = &(geno_vvec_biter[cur_incr * kBytesPerVec]);
    do {
      // hmm, this seems to have more linear dependence than I'd want, but the
      // reorderings I tried just made the code harder to read without helping,
      // so I'll leave this alone
      VecW cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW odd1 = m1 & vecw_srli(cur_geno_vword, 1);
      VecW even1 = m1 & cur_geno_vword;
      VecW bothset1 = odd1 & cur_geno_vword;

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (m1 & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high;
      bothset1 = bothset1 + (cur_geno_vword_high & cur_geno_vword);

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (m1 & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high;
      bothset1 = bothset1 + (cur_geno_vword_high & cur_geno_vword);

      even1 = (even1 & m2) + (vecw_srli(even1, 2) & m2);
      odd1 = (odd1 & m2) + (vecw_srli(odd1, 2) & m2);
      bothset1 = (bothset1 & m2) + (vecw_srli(bothset1, 2) & m2);

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW odd2 = m1 & vecw_srli(cur_geno_vword, 1);
      VecW even2 = m1 & cur_geno_vword;
      VecW bothset2 = odd2 & cur_geno_vword;

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high = m1 & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (m1 & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high;
      bothset2 = bothset2 + (cur_geno_vword_high & cur_geno_vword);

      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
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
    } while (geno_vvec_biter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_even = acc_even + vecw_bytesum(inner_acc_even, m0);
    acc_odd = acc_odd + vecw_bytesum(inner_acc_odd, m0);
    acc_bothset = acc_bothset + vecw_bytesum(inner_acc_bothset, m0);
  }
}

void FillInterleavedMaskVec(const uintptr_t* __restrict subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec) {
#ifdef USE_SSE2
  // This is a cousin of github.com/KWillets/simd_interleave , which was
  // written in response to
  //   https://lemire.me/blog/2018/01/09/how-fast-can-you-bit-interleave-32-bit-integers-simd-edition/
  // AVX2 implementation takes ~40% less time than before, and SSE4.2 takes
  // ~65% less.  This also avoids the Ryzen-screwing _pdep_u64()/_pext_u64()
  // operations.
  const VecW m4 = VCONST_W(kMask0F0F);
#  ifdef USE_SHUFFLE8
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
    VecW vec_lo;
    VecW vec_hi;
    vecw_lo_and_hi_nybbles(subset_mask_valias[vidx], m4, &vec_lo, &vec_hi);
#  ifdef USE_SHUFFLE8
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
#else  // !USE_SSE2
  for (uint32_t widx = 0; widx != base_vec_ct; ++widx) {
    const uintptr_t orig_word = subset_mask[widx];
    uintptr_t ww_even = S_CAST(Halfword, orig_word);
    uintptr_t ww_odd = orig_word >> kBitsPerWordD2;
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
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (3 * kWordsPerVec));
  uint32_t alt1_plus_bothset_ct;
  uint32_t bothset_ct;
  assert(IsVecAligned(genovec));
  Count2FreqVec3(R_CAST(const VecW*, genovec), word_idx / kWordsPerVec, &alt1_plus_bothset_ct, &bothset_ct);
  for (; word_idx != sample_ctl2; ++word_idx) {
    const uintptr_t cur_geno_word = genovec[word_idx];
    const uintptr_t cur_geno_word_low_lshifted = (cur_geno_word & kMask5555) << 1;
    alt1_plus_bothset_ct += NypsumWord((~cur_geno_word_low_lshifted) & cur_geno_word);
    bothset_ct += NypsumWord(cur_geno_word_low_lshifted & cur_geno_word);
  }
  const uint32_t alt1_ct = alt1_plus_bothset_ct - bothset_ct;
  allele_cts[0] = (sample_ct - bothset_ct) * 2 - alt1_ct;
  allele_cts[1] = alt1_ct;
  *bothset_ctp = bothset_ct;
}
*/

void GenoarrCountFreqsUnsafe(const uintptr_t* genoarr, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // fills genocounts[0] with the number of 00s, genocounts[1] with the number
  // of 01s, etc.
  // assumes trailing bits of last genoarr word are zeroed out.
  // sample_ct == 0 ok.
  // no longer requires vector-alignment.
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uint32_t even_ct;
  uint32_t odd_ct;
  uint32_t bothset_ct;
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (6 * kWordsPerVec));
  Count3FreqVec6(genoarr, word_idx / kWordsPerVec, &even_ct, &odd_ct, &bothset_ct);
  for (; word_idx != sample_ctl2; ++word_idx) {
    const uintptr_t cur_geno_word = genoarr[word_idx];
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

uintptr_t MostCommonGenoUnsafe(const uintptr_t* genoarr, uint32_t sample_ct) {
  STD_ARRAY_DECL(uint32_t, 4, genocounts);
  GenoarrCountFreqsUnsafe(genoarr, sample_ct, genocounts);
  uint32_t most_common_geno_ct = genocounts[0];
  if (most_common_geno_ct * 2 >= sample_ct) {
    return 0;
  }
  uintptr_t most_common_geno = 0;
  for (uintptr_t cur_geno = 1; cur_geno != 4; ++cur_geno) {
    if (most_common_geno_ct < genocounts[cur_geno]) {
      most_common_geno = cur_geno;
      most_common_geno_ct = genocounts[cur_geno];
    }
  }
  return most_common_geno;
}

// genoarr now allowed to be unaligned.
void CountSubset3FreqVec6(const void* __restrict genoarr, const VecW* __restrict interleaved_mask_vvec, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp) {
  assert(!(vec_ct % 6));
  // Sets even_ct to the number of set low bits in the current block after
  // subsetting, odd_ct to the number of set high bits, and bothset_ct by the
  // number of 0b11s.
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const unsigned char* geno_vvec_biter = S_CAST(const unsigned char*, genoarr);
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
    const unsigned char* geno_vvec_stop = &(geno_vvec_biter[cur_incr * kBytesPerVec]);
    do {
      VecW interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      VecW cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW cur_mask = interleaved_mask_vword & m1;
      VecW odd1 = cur_mask & vecw_srli(cur_geno_vword, 1);
      VecW even1 = cur_mask & cur_geno_vword;
      VecW bothset1 = odd1 & cur_geno_vword;

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (cur_mask & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high_masked;
      bothset1 = bothset1 + (cur_geno_vword_high_masked & cur_geno_vword);

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even1 = even1 + (cur_mask & cur_geno_vword);
      odd1 = odd1 + cur_geno_vword_high_masked;
      bothset1 = bothset1 + (cur_geno_vword_high_masked & cur_geno_vword);

      even1 = (even1 & m2) + (vecw_srli(even1, 2) & m2);
      odd1 = (odd1 & m2) + (vecw_srli(odd1, 2) & m2);
      bothset1 = (bothset1 & m2) + (vecw_srli(bothset1, 2) & m2);

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      VecW odd2 = cur_mask & vecw_srli(cur_geno_vword, 1);
      VecW even2 = cur_mask & cur_geno_vword;
      VecW bothset2 = odd2 & cur_geno_vword;

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
      cur_geno_vword_high_masked = cur_mask & vecw_srli(cur_geno_vword, 1);
      even2 = even2 + (cur_mask & cur_geno_vword);
      odd2 = odd2 + cur_geno_vword_high_masked;
      bothset2 = bothset2 + (cur_geno_vword_high_masked & cur_geno_vword);

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_biter);
      geno_vvec_biter += kBytesPerVec;
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
    } while (geno_vvec_biter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_even = acc_even + vecw_bytesum(inner_acc_even, m0);
    acc_odd = acc_odd + vecw_bytesum(inner_acc_odd, m0);
    acc_bothset = acc_bothset + vecw_bytesum(inner_acc_bothset, m0);
  }
}

void GenoarrCountSubsetFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // fills genocounts[0] with the number of 00s, genocounts[1] with the number
  // of 01s, etc.
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctv2 = NypCtToVecCt(raw_sample_ct);
  uint32_t even_ct;
  uint32_t odd_ct;
  uint32_t bothset_ct;
#ifdef USE_SSE2
  uint32_t vec_idx = raw_sample_ctv2 - (raw_sample_ctv2 % 6);
  CountSubset3FreqVec6(genoarr, R_CAST(const VecW*, sample_include_interleaved_vec), vec_idx, &even_ct, &odd_ct, &bothset_ct);
  const uintptr_t* genoarr_iter = &(genoarr[kWordsPerVec * vec_idx]);
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
      const uintptr_t cur_geno_word1 = *genoarr_iter++;
      const uintptr_t cur_geno_word2 = *genoarr_iter++;
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
    const uintptr_t cur_geno_word1 = *genoarr_iter++;
    const uintptr_t cur_geno_word2 = *genoarr_iter++;
    const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
    const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
#    ifdef USE_SSE42
    even_ct += PopcountWord(((cur_geno_word1 & mask_word1) << 1) | (cur_geno_word2 & mask_word2));
    odd_ct += PopcountWord((cur_geno_word1_high_masked << 1) | cur_geno_word2_high_masked);
    bothset_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
#    else
    even_ct += NypsumWord((cur_geno_word1 & mask_word1) + (cur_geno_word2 & mask_word2));
    odd_ct += NypsumWord(cur_geno_word1_high_masked + cur_geno_word2_high_masked);
    bothset_ct += NypsumWord((cur_geno_word1 & cur_geno_word1_high_masked) + (cur_geno_word2 & cur_geno_word2_high_masked));
#    endif
  }
#  endif  // not USE_AVX2
#else  // not USE_SSE2
  uint32_t word_idx = raw_sample_ctv2 - (raw_sample_ctv2 % 6);
  CountSubset3FreqVec6(genoarr, R_CAST(const VecW*, sample_include_interleaved_vec), word_idx, &even_ct, &odd_ct, &bothset_ct);
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
    const uintptr_t cur_geno_word = genoarr[word_idx];
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

void GenoarrCountSubsetFreqs2(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // slower GenoarrCountSubsetFreqs() which does not require
  // sample_include_interleaved_vec to be precomputed.
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctl2 = NypCtToWordCt(raw_sample_ct);
  const uint32_t fullword_ct = raw_sample_ctl2 / 2;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  for (uint32_t widx = 0; widx != fullword_ct; ++widx) {
    // possible todo: try vectorizing this loop in SSE4.2+ high-sample-ct case
    // with shuffle-based dynamic unpacking of sample_include?
    const uintptr_t mask_word = sample_include[widx];
    if (mask_word) {
      uintptr_t geno_word = genoarr[2 * widx];
      uintptr_t geno_even = PackWordToHalfwordMask5555(geno_word);
      uintptr_t geno_odd = PackWordToHalfwordMaskAAAA(geno_word);
      geno_word = genoarr[2 * widx + 1];
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

void GenoarrCountInvsubsetFreqs2(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_exclude, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // {raw_}sample_ct == 0 ok.
  // ugh, 'fullword' is overloaded.  probable todo: keep halfword/fullword,
  // switch more common use case to bodyword/trailword.
  const uint32_t bodyword_ct = raw_sample_ct / kBitsPerWord;
  uint32_t even_ct = 0;
  uint32_t odd_ct = 0;
  uint32_t bothset_ct = 0;
  for (uint32_t widx = 0; widx != bodyword_ct; ++widx) {
    // possible todo: try vectorizing this loop in SSE4.2+ high-sample-ct case
    // with shuffle-based dynamic unpacking of sample_exclude?
    const uintptr_t mask_word = ~sample_exclude[widx];
    if (mask_word) {
      uintptr_t geno_word = genoarr[2 * widx];
      uintptr_t geno_even = PackWordToHalfwordMask5555(geno_word);
      uintptr_t geno_odd = PackWordToHalfwordMaskAAAA(geno_word);
      geno_word = genoarr[2 * widx + 1];
      geno_even |= S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_word)) << kBitsPerWordD2;
      geno_odd |= S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(geno_word)) << kBitsPerWordD2;
      const uintptr_t geno_even_masked = geno_even & mask_word;
      even_ct += PopcountWord(geno_even_masked);
      odd_ct += PopcountWord(geno_odd & mask_word);
      bothset_ct += PopcountWord(geno_odd & geno_even_masked);
    }
  }
  const uint32_t remainder = raw_sample_ct % kBitsPerWord;
  if (remainder) {
    const uintptr_t mask_word = bzhi(~sample_exclude[bodyword_ct], remainder);
    if (mask_word) {
      uintptr_t geno_word = genoarr[2 * bodyword_ct];
      uintptr_t geno_even = PackWordToHalfwordMask5555(geno_word);
      uintptr_t geno_odd = PackWordToHalfwordMaskAAAA(geno_word);
      if (remainder > kBitsPerWordD2) {
        geno_word = genoarr[2 * bodyword_ct + 1];
        geno_even |= S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_word)) << kBitsPerWordD2;
        geno_odd |= S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(geno_word)) << kBitsPerWordD2;
      }
      const uintptr_t geno_even_masked = geno_even & mask_word;
      even_ct += PopcountWord(geno_even_masked);
      odd_ct += PopcountWord(geno_odd & mask_word);
      bothset_ct += PopcountWord(geno_odd & geno_even_masked);
    }
  }
  genocounts[0] = sample_ct + bothset_ct - even_ct - odd_ct;
  genocounts[1] = even_ct - bothset_ct;
  genocounts[2] = odd_ct - bothset_ct;
  genocounts[3] = bothset_ct;
}

void GenoarrCountSubsetIntersectFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict subset1, const uintptr_t* __restrict subset2, uint32_t raw_sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  // {raw_}sample_ct == 0 ok.
  const uint32_t raw_sample_ctl2 = NypCtToWordCt(raw_sample_ct);
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

void GenovecInvertUnsafe(uint32_t sample_ct, uintptr_t* genovec) {
  // flips 0 to 2 and vice versa.
  // "unsafe" because trailing bits are not zeroed out.
  const uint32_t vec_ct = NypCtToVecCt(sample_ct);
  assert(IsVecAligned(genovec));
  const VecW not_m1 = VCONST_W(kMaskAAAA);
  VecW* vptr = R_CAST(VecW*, genovec);
  for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
    VecW cur_vec = vptr[vidx];
    // flip high bit iff low bit is unset
    vptr[vidx] = cur_vec ^ vecw_and_notfirst(vecw_slli(cur_vec, 1), not_m1);
  }
}

void DifflistCountSubsetFreqs(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t common_geno, uint32_t difflist_len, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts) {
  STD_ARRAY_REF_FILL0(4, genocounts);
  for (uint32_t difflist_idx = 0; difflist_idx != difflist_len; ++difflist_idx) {
    const uint32_t raw_sample_idx = difflist_sample_ids[difflist_idx];
    if (IsSet(sample_include, raw_sample_idx)) {
      genocounts[GetNyparrEntry(raregeno, difflist_idx)] += 1;
    }
  }
  genocounts[common_geno] = sample_ct - genocounts[0] - genocounts[1] - genocounts[2] - genocounts[3];
}


#ifdef USE_SSE2
static_assert(kPglNypTransposeBatch == S_CAST(uint32_t, kNypsPerCacheline), "TransposeNypblock64() needs to be updated.");
#  ifdef CACHELINE64
void TransposeNypblock64(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1) {
  // buf0 and buf1 must each be vector-aligned and have size 16k
  // Tried using previous AVX2 small-buffer approach, but that was a bit
  // worse... though maybe it should be revisited?

  // Each input row has 256 nyps, across 8 words.
  // First word of each row goes into first buf0 row, etc.
  const uint32_t buf0_row_ct = DivUp(write_batch_size, 32);
  {
    // Fold the first 3 shuffles into the ingestion loop.
    // Can fold 2 additional shuffles here by ingesting uint16s instead.  That
    // removes the need for the second loop, halving the workspace requirement.
    // Benchmark results of that approach are slightly but consistently worse,
    // though.
    const uintptr_t* initial_read_iter = read_iter;
    const uintptr_t* initial_read_end = &(initial_read_iter[buf0_row_ct]);
    uintptr_t* initial_target_iter = R_CAST(uintptr_t*, buf0);
    const uint32_t read_batch_rem = kNypsPerCacheline - read_batch_size;
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
#    ifdef USE_SHUFFLE8
    const VecW gather_u32s = vecw_setr8(0, 1, 8, 9, 2, 3, 10, 11,
                                        4, 5, 12, 13, 6, 7, 14, 15);
#    else
    const VecW m16 = VCONST_W(kMask0000FFFF);
#    endif
    for (uint32_t bidx = 0; bidx != buf0_row_ct; ++bidx) {
      __m128i* write_iter1 = &(write_iter0[32]);
      __m128i* write_iter2 = &(write_iter1[32]);
      __m128i* write_iter3 = &(write_iter2[32]);
      for (uint32_t clidx = 0; clidx != buf0_row_clwidth; ++clidx) {
#    ifdef USE_AVX2
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
        const __m256i vec_lo = _mm256_shuffle_epi32(WToVec(loader0), 0xd8);
        const __m256i vec_hi = _mm256_shuffle_epi32(WToVec(loader1), 0xd8);
        const __m256i final0145 = _mm256_unpacklo_epi64(vec_lo, vec_hi);
        const __m256i final2367 = _mm256_unpackhi_epi64(vec_lo, vec_hi);
        // GCC doesn't support _mm256_storeu_si128i as of this writing.
        write_iter0[clidx] = _mm256_castsi256_si128(final0145);
        write_iter1[clidx] = _mm256_castsi256_si128(final2367);
        write_iter2[clidx] = _mm256_extracti128_si256(final0145, 1);
        write_iter3[clidx] = _mm256_extracti128_si256(final2367, 1);
#    else
        VecW loader0 = buf0_read_iter[clidx * 4];
        VecW loader1 = buf0_read_iter[clidx * 4 + 1];
        VecW loader2 = buf0_read_iter[clidx * 4 + 2];
        VecW loader3 = buf0_read_iter[clidx * 4 + 3];
#      ifdef USE_SHUFFLE8
        loader0 = vecw_shuffle8(loader0, gather_u32s);
        loader1 = vecw_shuffle8(loader1, gather_u32s);
        loader2 = vecw_shuffle8(loader2, gather_u32s);
        loader3 = vecw_shuffle8(loader3, gather_u32s);
#      else
        VecW tmp_lo = vecw_unpacklo16(loader0, loader1);
        VecW tmp_hi = vecw_unpackhi16(loader0, loader1);
        loader0 = vecw_blendv(vecw_slli(tmp_hi, 16), tmp_lo, m16);
        loader1 = vecw_blendv(tmp_hi, vecw_srli(tmp_lo, 16), m16);
        tmp_lo = vecw_unpacklo16(loader2, loader3);
        tmp_hi = vecw_unpackhi16(loader2, loader3);
        loader2 = vecw_blendv(vecw_slli(tmp_hi, 16), tmp_lo, m16);
        loader3 = vecw_blendv(tmp_hi, vecw_srli(tmp_lo, 16), m16);
#      endif
        //    (0,0) (0,1) (1,0) ... (7,1) (0,2) ... (7,3) (8,0) ... (31,3)
        //  + (0,4) ... (31,7)
        // -> (0,0) ... (7,3) (0,4) ... (7,7) (8,0) ... (15,7)
        const VecW lo_0_15 = vecw_unpacklo32(loader0, loader1);
        const VecW lo_16_31 = vecw_unpackhi32(loader0, loader1);
        const VecW hi_0_15 = vecw_unpacklo32(loader2, loader3);
        const VecW hi_16_31 = vecw_unpackhi32(loader2, loader3);
        write_iter0[clidx] = WToVec(vecw_unpacklo64(lo_0_15, hi_0_15));
        write_iter1[clidx] = WToVec(vecw_unpackhi64(lo_0_15, hi_0_15));
        write_iter2[clidx] = WToVec(vecw_unpacklo64(lo_16_31, hi_16_31));
        write_iter3[clidx] = WToVec(vecw_unpackhi64(lo_16_31, hi_16_31));
#    endif
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
  Vec8thUint* target_iter0 = DowncastWToV8(write_iter);
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
      // todo: better ARM implementation
      VecW target_4567 = vecw_blendv(loader, vecw_srli(loader, 7), m8);
      target_iter7[vidx] = vecw_movemask(target_4567);
      target_4567 = vecw_slli(target_4567, 2);
      target_iter6[vidx] = vecw_movemask(target_4567);
      target_4567 = vecw_slli(target_4567, 2);
      target_iter5[vidx] = vecw_movemask(target_4567);
      target_4567 = vecw_slli(target_4567, 2);
      target_iter4[vidx] = vecw_movemask(target_4567);
      VecW target_0123 = vecw_blendv(vecw_slli(loader, 8), vecw_slli(loader, 1), m8);
      target_iter3[vidx] = vecw_movemask(target_0123);
      target_0123 = vecw_slli(target_0123, 2);
      target_iter2[vidx] = vecw_movemask(target_0123);
      target_0123 = vecw_slli(target_0123, 2);
      target_iter1[vidx] = vecw_movemask(target_0123);
      target_0123 = vecw_slli(target_0123, 2);
      target_iter0[vidx] = vecw_movemask(target_0123);
    }
    source_iter = &(source_iter[(2 * kPglNypTransposeBatch) / kBytesPerVec]);
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
    VecW target_0123 = vecw_blendv(vecw_slli(loader, 8), vecw_slli(loader, 1), m8);
    target_iter3[vidx] = vecw_movemask(target_0123);
    target_0123 = vecw_slli(target_0123, 2);
    target_iter2[vidx] = vecw_movemask(target_0123);
    target_0123 = vecw_slli(target_0123, 2);
    target_iter1[vidx] = vecw_movemask(target_0123);
    target_0123 = vecw_slli(target_0123, 2);
    target_iter0[vidx] = vecw_movemask(target_0123);
  }
}
#  else
#    ifndef CACHELINE128
#      error "CACHELINE64 or CACHELINE128 expected."
#    endif
// assumes USE_SHUFFLE8, !USE_AVX2
void TransposeNypblock64(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1) {
  // buf0 and buf1 must each be vector-aligned and have size 64k
  // Each input row has 512 nyps, across 16 words.
  // First word of each row goes into first buf0 row, etc.
  const uint32_t buf0_row_ct = DivUp(write_batch_size, 32);
  {
    // Fold the first 4 shuffles into the ingestion loop.
    const uintptr_t* initial_read_iter = read_iter;
    const uintptr_t* initial_read_end = &(initial_read_iter[buf0_row_ct]);
    uintptr_t* initial_target_iter = R_CAST(uintptr_t*, buf0);
    const uint32_t read_batch_rem = kNypsPerCacheline - read_batch_size;
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

  // First buf0 row now corresponds to a 512x32 nyp matrix (512 * 8 bytes) that
  // we wish to transpose.  We split this into eight 512x4 matrices.
  // (ARMv8 doesn't have efficient movemask, so this should be better than four
  // 512x8 matrices.)
  // This is nearly identical to the middle step in TransposeBitblock64().
  {
    const VecW* buf0_read_iter = R_CAST(const VecW*, buf0);
    uintptr_t* write_iter0 = R_CAST(uintptr_t*, buf1);
    const VecW gather_u16s = vecw_setr8(0, 8, 1, 9, 2, 10, 3, 11,
                                        4, 12, 5, 13, 6, 14, 7, 15);
    const uint32_t buf0_row_b64width = DivUp(read_batch_size, 8);
    for (uint32_t ridx = 0; ridx != buf0_row_ct; ++ridx) {
      uintptr_t* write_iter1 = &(write_iter0[64]);
      uintptr_t* write_iter2 = &(write_iter1[64]);
      uintptr_t* write_iter3 = &(write_iter2[64]);
      uintptr_t* write_iter4 = &(write_iter3[64]);
      uintptr_t* write_iter5 = &(write_iter4[64]);
      uintptr_t* write_iter6 = &(write_iter5[64]);
      uintptr_t* write_iter7 = &(write_iter6[64]);
      for (uint32_t b64idx = 0; b64idx != buf0_row_b64width; ++b64idx) {
        VecW loader0 = buf0_read_iter[b64idx * 4];
        VecW loader1 = buf0_read_iter[b64idx * 4 + 1];
        VecW loader2 = buf0_read_iter[b64idx * 4 + 2];
        VecW loader3 = buf0_read_iter[b64idx * 4 + 3];
        loader0 = vecw_shuffle8(loader0, gather_u16s);
        loader1 = vecw_shuffle8(loader1, gather_u16s);
        loader2 = vecw_shuffle8(loader2, gather_u16s);
        loader3 = vecw_shuffle8(loader3, gather_u16s);
        const VecW lo_0123 = vecw_unpacklo16(loader0, loader1);
        const VecW lo_4567 = vecw_unpackhi16(loader0, loader1);
        const VecW hi_0123 = vecw_unpacklo16(loader2, loader3);
        const VecW hi_4567 = vecw_unpackhi16(loader2, loader3);

        const VecW final01 = vecw_unpacklo32(lo_0123, hi_0123);
        const VecW final23 = vecw_unpackhi32(lo_0123, hi_0123);
        const VecW final45 = vecw_unpacklo32(lo_4567, hi_4567);
        const VecW final67 = vecw_unpackhi32(lo_4567, hi_4567);
        write_iter0[b64idx] = vecw_extract64_0(final01);
        write_iter1[b64idx] = vecw_extract64_1(final01);
        write_iter2[b64idx] = vecw_extract64_0(final23);
        write_iter3[b64idx] = vecw_extract64_1(final23);
        write_iter4[b64idx] = vecw_extract64_0(final45);
        write_iter5[b64idx] = vecw_extract64_1(final45);
        write_iter6[b64idx] = vecw_extract64_0(final67);
        write_iter7[b64idx] = vecw_extract64_1(final67);
      }
      buf0_read_iter = &(buf0_read_iter[256]);
      write_iter0 = &(write_iter7[64]);
    }
  }

  // 8 -> 2
  // This is similar to the main TransposeNybbleblock() loop.
  const VecW* source_iter = R_CAST(VecW*, buf1);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  const VecW gather_even = vecw_setr8(0, 2, 4, 6, 8, 10, 12, 14,
                                      -1, -1, -1, -1, -1, -1, -1, -1);
  // Take advantage of current function contract.
  const uint32_t buf1_row_ct = (write_batch_size + 3) / 4;

  const uint32_t fourword_ct = DivUp(read_batch_size, 32);
  uintptr_t* target_iter0 = write_iter;
  for (uint32_t ridx = 0; ridx != buf1_row_ct; ++ridx) {
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    uintptr_t* target_iter2 = &(target_iter1[write_ul_stride]);
    uintptr_t* target_iter3 = &(target_iter2[write_ul_stride]);
    for (uint32_t dvidx = 0; dvidx != fourword_ct; ++dvidx) {
      const VecW loader0 = source_iter[dvidx * 2];
      const VecW loader1 = source_iter[dvidx * 2 + 1];

      VecW even_nyps0 = loader0 & m2;
      VecW even_nyps1 = loader1 & m2;
      VecW odd_nyps0 = vecw_srli(loader0, 2) & m2;
      VecW odd_nyps1 = vecw_srli(loader1, 2) & m2;
      even_nyps0 = even_nyps0 | vecw_srli(even_nyps0, 6);
      even_nyps1 = even_nyps1 | vecw_srli(even_nyps1, 6);
      odd_nyps0 = odd_nyps0 | vecw_srli(odd_nyps0, 6);
      odd_nyps1 = odd_nyps1 | vecw_srli(odd_nyps1, 6);
      // Low four bits of even_nyps{0,1}[0], [2], ..., [14] are destined for
      // target_iter0; high four bits of those bytes are destined for
      // target_iter2.
      const VecW even_nyps = vecw_gather_even(even_nyps0, even_nyps1, m8);
      const VecW odd_nyps = vecw_gather_even(odd_nyps0, odd_nyps1, m8);

      VecW mod0_nyps = even_nyps & m4;
      VecW mod1_nyps = odd_nyps & m4;
      VecW mod2_nyps = vecw_srli(even_nyps, 4) & m4;
      VecW mod3_nyps = vecw_srli(odd_nyps, 4) & m4;
      mod0_nyps = mod0_nyps | vecw_srli(mod0_nyps, 4);
      mod1_nyps = mod1_nyps | vecw_srli(mod1_nyps, 4);
      mod2_nyps = mod2_nyps | vecw_srli(mod2_nyps, 4);
      mod3_nyps = mod3_nyps | vecw_srli(mod3_nyps, 4);
      mod0_nyps = vecw_shuffle8(mod0_nyps, gather_even);
      mod1_nyps = vecw_shuffle8(mod1_nyps, gather_even);
      mod2_nyps = vecw_shuffle8(mod2_nyps, gather_even);
      mod3_nyps = vecw_shuffle8(mod3_nyps, gather_even);
      target_iter0[dvidx] = vecw_extract64_0(mod0_nyps);
      target_iter1[dvidx] = vecw_extract64_0(mod1_nyps);
      target_iter2[dvidx] = vecw_extract64_0(mod2_nyps);
      target_iter3[dvidx] = vecw_extract64_0(mod3_nyps);
    }
    source_iter = &(source_iter[32]);
    target_iter0 = &(target_iter3[write_ul_stride]);
  }
}
#  endif
#else  // !USE_SSE2
#  ifdef __LP64__
static_assert(kWordsPerVec == 1, "TransposeNypblock64() needs to be updated.");
void TransposeNypblock64(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1) {
  // buf0 and buf1 must each be vector-aligned and have size 16k
  // defining them as unsigned char* might prevent a strict-aliasing issue?
  // (might need to go through greater contortions to actually be safe?)
  const uint32_t buf_row_ct = NypCtToByteCt(write_batch_size);
  // fold the first 6 shuffles into the initial ingestion loop
  const unsigned char* initial_read_iter = DowncastKToUc(read_iter);
  const unsigned char* initial_read_end = &(initial_read_iter[buf_row_ct]);
  unsigned char* initial_target_iter = buf0;
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t read_batch_rem = kNypsPerCacheline - read_batch_size;
  for (; initial_read_iter != initial_read_end; ++initial_read_iter) {
    const unsigned char* read_iter_tmp = initial_read_iter;
    for (uint32_t ujj = 0; ujj != read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_byte_stride]);
    }
    initial_target_iter = memsetua(initial_target_iter, 0, read_batch_rem);
  }

  // second-to-last shuffle, 8 bit spacing -> 4
  const uintptr_t* source_iter = R_CAST(uintptr_t*, buf0);
  uintptr_t* target_iter0 = R_CAST(uintptr_t*, buf1);
  const uint32_t write_word_ct = NypCtToWordCt(read_batch_size);
  const uint32_t penult_inner_loop_iter_ct = 2 * write_word_ct;
  const uint32_t cur_write_skip = 2 * kWordsPerCacheline - penult_inner_loop_iter_ct;
  for (uint32_t uii = 0; uii != buf_row_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 2]);
    for (uint32_t ujj = 0; ujj != penult_inner_loop_iter_ct; ++ujj) {
      const uintptr_t source_word_lo = *source_iter++;
      const uintptr_t source_word_hi = *source_iter++;
      *target_iter0++ = PackTwo0F0FMask(source_word_lo, source_word_hi);
      *target_iter1++ = PackTwo0F0FMask(source_word_lo >> 4, source_word_hi >> 4);
    }
    source_iter = &(source_iter[2 * cur_write_skip]);
    target_iter0 = &(target_iter1[cur_write_skip]);
  }

  // last shuffle, 4 bit spacing -> 2
  source_iter = R_CAST(uintptr_t*, buf1);
  target_iter0 = write_iter;
  const uint32_t last_loop_iter_ct = DivUp(write_batch_size, 2);
  for (uint32_t uii = 0; uii != last_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    for (uint32_t ujj = 0; ujj != write_word_ct; ++ujj) {
      const uintptr_t source_word_lo = *source_iter++;
      const uintptr_t source_word_hi = *source_iter++;
      target_iter0[ujj] = PackTwo3333Mask(source_word_lo, source_word_hi);
      target_iter1[ujj] = PackTwo3333Mask(source_word_lo >> 2, source_word_hi >> 2);
    }
    source_iter = &(source_iter[2 * (kWordsPerCacheline - write_word_ct)]);
    target_iter0 = &(target_iter1[write_ul_stride]);
  }
}
#  else
static_assert(kWordsPerVec == 1, "TransposeNypblock32() needs to be updated.");
void TransposeNypblock32(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1) {
  // buf0 and buf1 must each be vector-aligned and have size 16k
  // defining them as unsigned char* might prevent a strict-aliasing issue?
  // (might need to go through greater contortions to actually be safe?)
  const uint32_t buf_row_ct = NypCtToByteCt(write_batch_size);
  // fold the first 6 shuffles into the initial ingestion loop
  const unsigned char* initial_read_iter = DowncastKToUc(read_iter);
  const unsigned char* initial_read_end = &(initial_read_iter[buf_row_ct]);
  unsigned char* initial_target_iter = buf0;
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t read_batch_rem = kNypsPerCacheline - read_batch_size;
  for (; initial_read_iter != initial_read_end; ++initial_read_iter) {
    const unsigned char* read_iter_tmp = initial_read_iter;
    for (uint32_t ujj = 0; ujj != read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_byte_stride]);
    }
    initial_target_iter = memsetua(initial_target_iter, 0, read_batch_rem);
  }

  // second-to-last shuffle, 8 bit spacing -> 4
  const uintptr_t* source_iter = R_CAST(uintptr_t*, buf0);
  uintptr_t* target_iter0 = R_CAST(uintptr_t*, buf1);
  const uint32_t write_word_ct = NypCtToWordCt(read_batch_size);
  const uint32_t penult_inner_loop_iter_ct = 2 * write_word_ct;
  const uint32_t cur_write_skip = 2 * kWordsPerCacheline - penult_inner_loop_iter_ct;
  for (uint32_t uii = 0; uii != buf_row_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 2]);
    for (uint32_t ujj = 0; ujj != penult_inner_loop_iter_ct; ++ujj) {
      const uintptr_t source_word_lo = *source_iter++;
      const uintptr_t source_word_hi = *source_iter++;
      *target_iter0++ = PackTwo0F0FMask(source_word_lo, source_word_hi);
      *target_iter1++ = PackTwo0F0FMask(source_word_lo >> 4, source_word_hi >> 4);
    }
    source_iter = &(source_iter[2 * cur_write_skip]);
    target_iter0 = &(target_iter1[cur_write_skip]);
  }

  // last shuffle, 4 bit spacing -> 2
  source_iter = R_CAST(uintptr_t*, buf1);
  target_iter0 = write_iter;
  const uint32_t last_loop_iter_ct = DivUp(write_batch_size, 2);
  for (uint32_t uii = 0; uii != last_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    for (uint32_t ujj = 0; ujj != write_word_ct; ++ujj) {
      const uintptr_t source_word_lo = *source_iter++;
      const uintptr_t source_word_hi = *source_iter++;
      target_iter0[ujj] = PackTwo3333Mask(source_word_lo, source_word_hi);
      target_iter1[ujj] = PackTwo3333Mask(source_word_lo >> 2, source_word_hi >> 2);
    }
    source_iter = &(source_iter[2 * (kWordsPerCacheline - write_word_ct)]);
    target_iter0 = &(target_iter1[write_ul_stride]);
  }
}
#  endif  // !__LP64__
#endif  // !USE_SSE2

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

#if defined(USE_SSE2) && !defined(USE_AVX2)
void PackWordsToHalfwordsInvmatch(const uintptr_t* __restrict genoarr, uintptr_t inv_match_word, uint32_t inword_ct, uintptr_t* __restrict dst) {
  // In shuffle8 case, this takes ~30% less time than a
  // PackWordToHalfwordMask5555 loop.
  const uint32_t out_fullvec_ct = inword_ct / (kWordsPerVec * 2);
  const VecW xor_vec = vecw_set1(inv_match_word);
  const VecW m1 = VCONST_W(kMask5555);
#  ifdef USE_SHUFFLE8
  const VecW swap12 = vecw_setr8(
      0, 1, 4, 5, 2, 3, 6, 7,
      8, 9, 12, 13, 10, 11, 14, 15);
#  else
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  for (uintptr_t vidx = 0; vidx != out_fullvec_ct; ++vidx) {
    VecW vec_lo = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx])) ^ xor_vec;
    VecW vec_hi = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx + kWordsPerVec])) ^ xor_vec;
    vec_lo = vec_lo & vecw_srli(vec_lo, 1) & m1;
    vec_hi = vec_hi & vecw_srli(vec_hi, 1) & m1;
#  ifdef USE_SHUFFLE8
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
    vec_lo = vec_lo | vecw_srli(vec_lo, 4);
    vec_hi = vec_hi | vecw_srli(vec_hi, 4);
    const VecW vec_packed = vecw_gather_even(vec_lo, vec_hi, m8);
    vecw_storeu(&(dst[kWordsPerVec * vidx]), vec_packed);
  }
  Halfword* dst_alias = DowncastWToHW(dst);
  uint32_t widx = RoundDownPow2(inword_ct, kWordsPerVec * 2);
  for (; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ inv_match_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word & (cur_word >> 1));
    dst_alias[widx] = hw;
  }
}

void PackWordsToHalfwordsMismatch(const uintptr_t* __restrict genoarr, uintptr_t mismatch_word, uint32_t inword_ct, uintptr_t* __restrict dst) {
  const uint32_t out_fullvec_ct = inword_ct / (kWordsPerVec * 2);
  const VecW xor_vec = vecw_set1(mismatch_word);
  const VecW m1 = VCONST_W(kMask5555);
#  ifdef USE_SHUFFLE8
  const VecW swap12 = vecw_setr8(
      0, 1, 4, 5, 2, 3, 6, 7,
      8, 9, 12, 13, 10, 11, 14, 15);
#  else
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  for (uintptr_t vidx = 0; vidx != out_fullvec_ct; ++vidx) {
    VecW vec_lo = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx])) ^ xor_vec;
    VecW vec_hi = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx + kWordsPerVec])) ^ xor_vec;
    vec_lo = (vec_lo | vecw_srli(vec_lo, 1)) & m1;
    vec_hi = (vec_hi | vecw_srli(vec_hi, 1)) & m1;
#  ifdef USE_SHUFFLE8
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
    vec_lo = vec_lo | vecw_srli(vec_lo, 4);
    vec_hi = vec_hi | vecw_srli(vec_hi, 4);
    const VecW vec_packed = vecw_gather_even(vec_lo, vec_hi, m8);
    vecw_storeu(&(dst[kWordsPerVec * vidx]), vec_packed);
  }
  Halfword* dst_alias = DowncastWToHW(dst);
  uint32_t widx = RoundDownPow2(inword_ct, kWordsPerVec * 2);
  for (; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ mismatch_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word | (cur_word >> 1));
    dst_alias[widx] = hw;
  }
}

void MaskWordsToHalfwordsInvmatch(const uintptr_t* __restrict genoarr, uintptr_t inv_match_word, uint32_t inword_ct, uintptr_t* src, uintptr_t* dst) {
  const uint32_t out_fullvec_ct = inword_ct / (kWordsPerVec * 2);
  const VecW xor_vec = vecw_set1(inv_match_word);
  const VecW m1 = VCONST_W(kMask5555);
#  ifdef USE_SHUFFLE8
  const VecW swap12 = vecw_setr8(
      0, 1, 4, 5, 2, 3, 6, 7,
      8, 9, 12, 13, 10, 11, 14, 15);
#  else
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  for (uintptr_t vidx = 0; vidx != out_fullvec_ct; ++vidx) {
    VecW vec_lo = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx])) ^ xor_vec;
    VecW vec_hi = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx + kWordsPerVec])) ^ xor_vec;
    VecW src_vec = vecw_loadu(&(src[kWordsPerVec * vidx]));
    vec_lo = vec_lo & vecw_srli(vec_lo, 1) & m1;
    vec_hi = vec_hi & vecw_srli(vec_hi, 1) & m1;
#  ifdef USE_SHUFFLE8
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
    vec_lo = vec_lo | vecw_srli(vec_lo, 4);
    vec_hi = vec_hi | vecw_srli(vec_hi, 4);
    const VecW vec_packed = vecw_gather_even(vec_lo, vec_hi, m8);
    vecw_storeu(&(dst[kWordsPerVec * vidx]), src_vec & vec_packed);
  }
  Halfword* src_alias = DowncastWToHW(src);
  Halfword* dst_alias = DowncastWToHW(dst);
  uint32_t widx = RoundDownPow2(inword_ct, kWordsPerVec * 2);
  for (; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ inv_match_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word & (cur_word >> 1));
    dst_alias[widx] = src_alias[widx] & hw;
  }
}
#else
void PackWordsToHalfwordsInvmatch(const uintptr_t* __restrict genoarr, uintptr_t inv_match_word, uint32_t inword_ct, uintptr_t* __restrict dst) {
  Halfword* dst_alias = DowncastWToHW(dst);
  for (uint32_t widx = 0; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ inv_match_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word & (cur_word >> 1));
    dst_alias[widx] = hw;
  }
}

void PackWordsToHalfwordsMismatch(const uintptr_t* __restrict genoarr, uintptr_t mismatch_word, uint32_t inword_ct, uintptr_t* __restrict dst) {
  Halfword* dst_alias = DowncastWToHW(dst);
  for (uint32_t widx = 0; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ mismatch_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word | (cur_word >> 1));
    dst_alias[widx] = hw;
  }
}

void MaskWordsToHalfwordsInvmatch(const uintptr_t* __restrict genoarr, uintptr_t inv_match_word, uint32_t inword_ct, uintptr_t* src, uintptr_t* dst) {
  Halfword* src_alias = DowncastWToHW(src);
  Halfword* dst_alias = DowncastWToHW(dst);
  for (uint32_t widx = 0; widx != inword_ct; ++widx) {
    const uintptr_t cur_word = genoarr[widx] ^ inv_match_word;
    const Halfword hw = PackWordToHalfwordMask5555(cur_word & (cur_word >> 1));
    dst_alias[widx] = src_alias[widx] & hw;
  }
}
#endif

void SparseToMissingness(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uint32_t sample_ct, uint32_t difflist_common_geno, uint32_t difflist_len, uintptr_t* __restrict missingness) {
  if (difflist_common_geno != 3) {
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    ZeroWArr(sample_ctl, missingness);
    if (!difflist_len) {
      return;
    }
    const uint32_t raregeno_word_ct = NypCtToWordCt(difflist_len);
    for (uint32_t widx = 0; widx != raregeno_word_ct; ++widx) {
      const uintptr_t raregeno_word = raregeno[widx];
      uintptr_t raregeno_11 = raregeno_word & (raregeno_word >> 1) & kMask5555;
      if (raregeno_11) {
        const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
        do {
          const uint32_t sample_idx_lowbits = ctzw(raregeno_11) / 2;
          const uint32_t cur_sample_id = cur_difflist_sample_ids[sample_idx_lowbits];
          SetBit(cur_sample_id, missingness);
          raregeno_11 &= raregeno_11 - 1;
        } while (raregeno_11);
      }
    }
  } else {
    SetAllBits(sample_ct, missingness);
    // Don't need to look at raregeno, all cases are nonmissing.
    for (uint32_t uii = 0; uii != difflist_len; ++uii) {
      const uint32_t cur_sample_id = difflist_sample_ids[uii];
      ClearBit(cur_sample_id, missingness);
    }
  }
}

#if defined(USE_SSE2) && !defined(USE_AVX2)
void SplitHomRef2hetUnsafeW(const uintptr_t* __restrict genoarr, uint32_t inword_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf) {
  const uint32_t out_fullvec_ct = inword_ct / (kWordsPerVec * 2);
  // In shuffle8 case, this takes ~55% less time than unvectorized loop.
  // homozygous: geno = 0 or 2
  // ref2het: geno = 0 or 1
  const VecW m1 = VCONST_W(kMask5555);
#  ifdef USE_SHUFFLE8
  const VecW swap12 = vecw_setr8(
      0, 1, 4, 5, 2, 3, 6, 7,
      8, 9, 12, 13, 10, 11, 14, 15);
#  else
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  for (uintptr_t vidx = 0; vidx != out_fullvec_ct; ++vidx) {
    const VecW vec_lo = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx]));
    const VecW vec_hi = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx + kWordsPerVec]));
    VecW hom_lo = vecw_and_notfirst(vec_lo, m1);
    VecW r2h_lo = vecw_and_notfirst(vecw_srli(vec_lo, 1), m1);
    VecW hom_hi = vecw_and_notfirst(vec_hi, m1);
    VecW r2h_hi = vecw_and_notfirst(vecw_srli(vec_hi, 1), m1);
#  ifdef USE_SHUFFLE8
    hom_lo = (hom_lo | vecw_srli(hom_lo, 3)) & m4;
    r2h_lo = (r2h_lo | vecw_srli(r2h_lo, 3)) & m4;
    hom_hi = (hom_hi | vecw_srli(hom_hi, 3)) & m4;
    r2h_hi = (r2h_hi | vecw_srli(r2h_hi, 3)) & m4;
    hom_lo = vecw_shuffle8(swap12, hom_lo);
    r2h_lo = vecw_shuffle8(swap12, r2h_lo);
    hom_hi = vecw_shuffle8(swap12, hom_hi);
    r2h_hi = vecw_shuffle8(swap12, r2h_hi);
#  else
    hom_lo = (hom_lo | vecw_srli(hom_lo, 1)) & m2;
    r2h_lo = (r2h_lo | vecw_srli(r2h_lo, 1)) & m2;
    hom_hi = (hom_hi | vecw_srli(hom_hi, 1)) & m2;
    r2h_hi = (r2h_hi | vecw_srli(r2h_hi, 1)) & m2;
    hom_lo = (hom_lo | vecw_srli(hom_lo, 2)) & m4;
    r2h_lo = (r2h_lo | vecw_srli(r2h_lo, 2)) & m4;
    hom_hi = (hom_hi | vecw_srli(hom_hi, 2)) & m4;
    r2h_hi = (r2h_hi | vecw_srli(r2h_hi, 2)) & m4;
#  endif
    hom_lo = hom_lo | vecw_srli(hom_lo, 4);
    r2h_lo = r2h_lo | vecw_srli(r2h_lo, 4);
    hom_hi = hom_hi | vecw_srli(hom_hi, 4);
    r2h_hi = r2h_hi | vecw_srli(r2h_hi, 4);
    const VecW hom_packed = vecw_gather_even(hom_lo, hom_hi, m8);
    const VecW r2h_packed = vecw_gather_even(r2h_lo, r2h_hi, m8);
    vecw_storeu(&(hom_buf[kWordsPerVec * vidx]), hom_packed);
    vecw_storeu(&(ref2het_buf[kWordsPerVec * vidx]), r2h_packed);
  }
  Halfword* hom_alias = DowncastWToHW(hom_buf);
  Halfword* r2h_alias = DowncastWToHW(ref2het_buf);
  uint32_t widx = RoundDownPow2(inword_ct, kWordsPerVec * 2);
  for (; widx != inword_ct; ++widx) {
    const uintptr_t inv_geno_word = ~genoarr[widx];
    hom_alias[widx] = PackWordToHalfwordMask5555(inv_geno_word);
    r2h_alias[widx] = PackWordToHalfwordMaskAAAA(inv_geno_word);
  }
}
#else
void SplitHomRef2hetUnsafeW(const uintptr_t* genoarr, uint32_t inword_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf) {
  Halfword* hom_alias = DowncastWToHW(hom_buf);
  Halfword* ref2het_alias = DowncastWToHW(ref2het_buf);
  for (uint32_t widx = 0; widx != inword_ct; ++widx) {
    const uintptr_t inv_geno_word = ~genoarr[widx];
    hom_alias[widx] = PackWordToHalfwordMask5555(inv_geno_word);
    ref2het_alias[widx] = PackWordToHalfwordMaskAAAA(inv_geno_word);
  }
}
#endif

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

BoolErr HapsplitMustPhased(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, uint32_t sample_ct, uint32_t phase_exists, uintptr_t* hap_arr, uintptr_t* nm_arr) {
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uintptr_t detect_unphased = 0;
  if (!phase_exists) {
    // error out if het encountered
    for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
      const uintptr_t geno_word = genoarr[widx];
      const uintptr_t nm_word = 3 * (kMask5555 & (~(geno_word & (geno_word >> 1))));
      const uintptr_t geno_nm = geno_word & nm_word;
      // geno_nm is now {00, 01, 10, 00}.  In particular, if we ever see 01, we
      // error out.
      const uintptr_t geno_nm_hi = (geno_nm >> 1) & kMask5555;
      nm_arr[widx] = nm_word;
      hap_arr[widx] = geno_nm | geno_nm_hi;
      detect_unphased |= geno_nm;
    }
    detect_unphased &= kMask5555;
  } else {
    // error out if het encountered, and not covered by phasepresent
    const Halfword* phasepresent_alias = DowncastKWToHW(phasepresent);
    const Halfword* phaseinfo_alias = DowncastKWToHW(phaseinfo);
    for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
      const uintptr_t geno_word = genoarr[widx];
      const uintptr_t geno_rshift = geno_word >> 1;
      const uintptr_t nm_word = 3 * (kMask5555 & (~(geno_word & geno_rshift)));
      const uintptr_t geno_nm = geno_word & nm_word;
      const uintptr_t geno_nm_hi = (geno_nm >> 1) & kMask5555;
      const uintptr_t het_word = geno_nm & kMask5555;
      const uintptr_t phasepresent_word = UnpackHalfwordToWord(phasepresent_alias[widx]);
      const uintptr_t phaseinfo_word = phasepresent_word & UnpackHalfwordToWord(phaseinfo_alias[widx]);
      nm_arr[widx] = nm_word;
      hap_arr[widx] = geno_nm + geno_nm_hi + phaseinfo_word;
      detect_unphased |= het_word & (~phasepresent_word);
    }
  }
  const uint32_t trailing_nyp_ct = sample_ct % kBitsPerWordD2;
  if (trailing_nyp_ct) {
    const uint32_t trailing_bit_ct = trailing_nyp_ct * 2;
    const uint32_t last_word_idx = sample_ctl2 - 1;
    nm_arr[last_word_idx] = bzhi(nm_arr[last_word_idx], trailing_bit_ct);
    hap_arr[last_word_idx] = bzhi(hap_arr[last_word_idx], trailing_bit_ct);
  }
  return (detect_unphased != 0);
}

#if defined(USE_SSE2) && !defined(USE_AVX2)
void HapsplitHaploid(const uintptr_t* __restrict genoarr, uint32_t sample_ct, uintptr_t* __restrict hap_arr, uintptr_t* __restrict nm_arr) {
  // In shuffle8 case, this takes ~40% less time than unvectorized loop.
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t sample_ctl2_is_odd = NypCtToWordCt(sample_ct) & 1;
  const uint32_t wordpair_ct = sample_ctl - sample_ctl2_is_odd;
  const uint32_t out_fullvec_ct = wordpair_ct / 2;
  const VecW m1 = VCONST_W(kMask5555);
#  ifdef USE_SHUFFLE8
  const VecW swap12 = vecw_setr8(
      0, 1, 4, 5, 2, 3, 6, 7,
      8, 9, 12, 13, 10, 11, 14, 15);
#  else
  const VecW m2 = VCONST_W(kMask3333);
#  endif
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  for (uintptr_t vidx = 0; vidx != out_fullvec_ct; ++vidx) {
    const VecW vec_lo = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx]));
    const VecW vec_hi = vecw_loadu(&(genoarr[2 * kWordsPerVec * vidx + kWordsPerVec]));
    VecW nm_lo = vecw_and_notfirst(vec_lo, m1);
    VecW nm_hi = vecw_and_notfirst(vec_hi, m1);
    VecW hap_lo = nm_lo & vecw_srli(vec_lo, 1);
    VecW hap_hi = nm_hi & vecw_srli(vec_hi, 1);
#  ifdef USE_SHUFFLE8
    nm_lo = (nm_lo | vecw_srli(nm_lo, 3)) & m4;
    nm_hi = (nm_hi | vecw_srli(nm_hi, 3)) & m4;
    hap_lo = (hap_lo | vecw_srli(hap_lo, 3)) & m4;
    hap_hi = (hap_hi | vecw_srli(hap_hi, 3)) & m4;
    nm_lo = vecw_shuffle8(swap12, nm_lo);
    nm_hi = vecw_shuffle8(swap12, nm_hi);
    hap_lo = vecw_shuffle8(swap12, hap_lo);
    hap_hi = vecw_shuffle8(swap12, hap_hi);
#  else
    nm_lo = (nm_lo | vecw_srli(nm_lo, 1)) & m2;
    nm_hi = (nm_hi | vecw_srli(nm_hi, 1)) & m2;
    hap_lo = (hap_lo | vecw_srli(hap_lo, 1)) & m2;
    hap_hi = (hap_hi | vecw_srli(hap_hi, 1)) & m2;
    nm_lo = (nm_lo | vecw_srli(nm_lo, 2)) & m4;
    nm_hi = (nm_hi | vecw_srli(nm_hi, 2)) & m4;
    hap_lo = (hap_lo | vecw_srli(hap_lo, 2)) & m4;
    hap_hi = (hap_hi | vecw_srli(hap_hi, 2)) & m4;
#  endif
    nm_lo = nm_lo | vecw_srli(nm_lo, 4);
    nm_hi = nm_hi | vecw_srli(nm_hi, 4);
    hap_lo = hap_lo | vecw_srli(hap_lo, 4);
    hap_hi = hap_hi | vecw_srli(hap_hi, 4);
    const VecW nm_packed = vecw_gather_even(nm_lo, nm_hi, m8);
    const VecW hap_packed = vecw_gather_even(hap_lo, hap_hi, m8);
    vecw_storeu(&(nm_arr[kWordsPerVec * vidx]), nm_packed);
    vecw_storeu(&(hap_arr[kWordsPerVec * vidx]), hap_packed);
  }
  if (wordpair_ct % 2) {
    const uint32_t widx = wordpair_ct - 1;
    const uintptr_t geno_word0 = genoarr[widx * 2];
    const uintptr_t geno_word1 = genoarr[widx * 2 + 1];
    const uintptr_t nm_word0 = ~geno_word0;
    const uintptr_t nm_word1 = ~geno_word1;
    const uintptr_t hap_word0 = nm_word0 & (geno_word0 >> 1);
    const uintptr_t hap_word1 = nm_word1 & (geno_word1 >> 1);
    nm_arr[widx] = PackTwo5555Mask(nm_word0, nm_word1);
    hap_arr[widx] = PackTwo5555Mask(hap_word0, hap_word1);
  }
  if (sample_ctl2_is_odd) {
    const uintptr_t geno_word0 = genoarr[wordpair_ct * 2];
    const uintptr_t nm_word0 = ~geno_word0;
    const uintptr_t hap_word0 = nm_word0 & (geno_word0 >> 1);
    nm_arr[wordpair_ct] = PackWordToHalfwordMask5555(nm_word0);
    hap_arr[wordpair_ct] = PackWordToHalfwordMask5555(hap_word0);
  }
  const uint32_t trailing_bit_ct = sample_ct % kBitsPerWord;
  if (trailing_bit_ct) {
    const uint32_t last_word_idx = sample_ctl - 1;
    nm_arr[last_word_idx] = bzhi(nm_arr[last_word_idx], trailing_bit_ct);
    hap_arr[last_word_idx] = bzhi(hap_arr[last_word_idx], trailing_bit_ct);
  }
}
#else
void HapsplitHaploid(const uintptr_t* __restrict genoarr, uint32_t sample_ct, uintptr_t* __restrict hap_arr, uintptr_t* __restrict nm_arr) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t sample_ctl2_is_odd = NypCtToWordCt(sample_ct) & 1;
  const uint32_t wordpair_ct = sample_ctl - sample_ctl2_is_odd;
  for (uint32_t widx = 0; widx != wordpair_ct; ++widx) {
    const uintptr_t geno_word0 = genoarr[widx * 2];
    const uintptr_t geno_word1 = genoarr[widx * 2 + 1];
    const uintptr_t nm_word0 = ~geno_word0;
    const uintptr_t nm_word1 = ~geno_word1;
    const uintptr_t hap_word0 = nm_word0 & (geno_word0 >> 1);
    const uintptr_t hap_word1 = nm_word1 & (geno_word1 >> 1);
    nm_arr[widx] = PackTwo5555Mask(nm_word0, nm_word1);
    hap_arr[widx] = PackTwo5555Mask(hap_word0, hap_word1);
  }
  if (sample_ctl2_is_odd) {
    const uintptr_t geno_word0 = genoarr[wordpair_ct * 2];
    const uintptr_t nm_word0 = ~geno_word0;
    const uintptr_t hap_word0 = nm_word0 & (geno_word0 >> 1);
    nm_arr[wordpair_ct] = PackWordToHalfwordMask5555(nm_word0);
    hap_arr[wordpair_ct] = PackWordToHalfwordMask5555(hap_word0);
  }
  const uint32_t trailing_bit_ct = sample_ct % kBitsPerWord;
  if (trailing_bit_ct) {
    const uint32_t last_word_idx = sample_ctl - 1;
    nm_arr[last_word_idx] = bzhi(nm_arr[last_word_idx], trailing_bit_ct);
    hap_arr[last_word_idx] = bzhi(hap_arr[last_word_idx], trailing_bit_ct);
  }
}
#endif

#ifdef USE_SHUFFLE8
// ~70% less time than per-byte lookup.
void GenoarrLookup256x1bx4(const uintptr_t* genoarr, const void* table256x1bx4, uint32_t sample_ct, void* __restrict result) {
  const uint32_t* table_alias = S_CAST(const uint32_t*, table256x1bx4);
  const unsigned char* genoarr_alias = DowncastKToUc(genoarr);
  unsigned char* resultb = S_CAST(unsigned char*, result);
  const uint32_t full_byte_ct = sample_ct / 4;
  if (full_byte_ct >= kBytesPerVec) {
    const uint32_t last_genoarr_offset = full_byte_ct - kBytesPerVec;
    const VecW lookup_even = vecw_loadu(&(table_alias[256]));
    const VecW lookup_odd = vecw_loadu(&(table_alias[256 + kInt32PerVec]));
    const VecW m4 = VCONST_W(kMask0F0F);
    for (uint32_t genoarr_offset = 0; ; genoarr_offset += kBytesPerVec) {
      if (genoarr_offset >= last_genoarr_offset) {
        if (genoarr_offset == full_byte_ct) {
          break;
        }
        genoarr_offset = last_genoarr_offset;
      }
      VecW cur_vec = vecw_loadu(&(genoarr_alias[genoarr_offset]));
      VecW vec_lo;
      VecW vec_hi;
      vecw_lo_and_hi_nybbles(cur_vec, m4, &vec_lo, &vec_hi);
      VecW result_lo_even = vecw_shuffle8(lookup_even, vec_lo);
      VecW result_hi_even = vecw_shuffle8(lookup_even, vec_hi);
      VecW result_lo_odd = vecw_shuffle8(lookup_odd, vec_lo);
      VecW result_hi_odd = vecw_shuffle8(lookup_odd, vec_hi);
      // In AVX2 case, result_lo_even has (0, 2, 4, ..., 62), and result_lo_odd
      // has (1, 3, 5, ..., 63).
      //   even -> (0, 2, 4, ..., 14, 32, ..., 46, 16, ..., 30, 48, ..., 62)
      //   odd ->  (1, 3, 5, ..., 15, 33, ..., 47, 17, ..., 31, 49, ..., 63)
      // Then unpacklo yields (0, 1, 2, ..., 31).
      result_lo_even = vecw_permute0xd8_if_avx2(result_lo_even);
      result_lo_odd = vecw_permute0xd8_if_avx2(result_lo_odd);
      result_hi_even = vecw_permute0xd8_if_avx2(result_hi_even);
      result_hi_odd = vecw_permute0xd8_if_avx2(result_hi_odd);
      const VecW result0 = vecw_unpacklo8(result_lo_even, result_lo_odd);
      const VecW result1 = vecw_unpackhi8(result_lo_even, result_lo_odd);
      const VecW result2 = vecw_unpacklo8(result_hi_even, result_hi_odd);
      const VecW result3 = vecw_unpackhi8(result_hi_even, result_hi_odd);
      unsigned char* cur_resultb = &(resultb[genoarr_offset * 4]);
      vecw_storeu(cur_resultb, result0);
      vecw_storeu(&(cur_resultb[kBytesPerVec]), result1);
      vecw_storeu(&(cur_resultb[2 * kBytesPerVec]), result2);
      vecw_storeu(&(cur_resultb[3 * kBytesPerVec]), result3);
    }
  } else {
    for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
      CopyToUnalignedOffsetU32(resultb, &(table_alias[genoarr_alias[byte_idx]]), byte_idx);
    }
  }
  const uint32_t remainder = sample_ct % 4;
  if (remainder) {
    unsigned char* result_last = &(resultb[full_byte_ct * 4]);
    uintptr_t geno_byte = genoarr_alias[full_byte_ct];
    for (uint32_t uii = 0; uii != remainder; ++uii) {
      result_last[uii] = table_alias[geno_byte & 3];
      geno_byte >>= 2;
    }
  }
}
#else
void GenoarrLookup256x1bx4(const uintptr_t* genoarr, const void* table256x1bx4, uint32_t sample_ct, void* __restrict result) {
  const uint32_t* table_alias = S_CAST(const uint32_t*, table256x1bx4);
  const unsigned char* genoarr_alias = DowncastKToUc(genoarr);
  unsigned char* resultb = S_CAST(unsigned char*, result);
  const uint32_t full_byte_ct = sample_ct / 4;
  for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
    CopyToUnalignedOffsetU32(resultb, &(table_alias[genoarr_alias[byte_idx]]), byte_idx);
  }
  const uint32_t remainder = sample_ct % 4;
  if (remainder) {
    unsigned char* result_last = &(resultb[full_byte_ct * 4]);
    uintptr_t geno_byte = genoarr_alias[full_byte_ct];
    for (uint32_t uii = 0; uii != remainder; ++uii) {
      result_last[uii] = table_alias[geno_byte & 3];
      geno_byte >>= 2;
    }
  }
}
#endif

#ifndef NO_UNALIGNED
void GenoarrLookup16x4bx2(const uintptr_t* genoarr, const void* table16x4bx2, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table16x4bx2);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
          // not a regular CopyToUnaligned, we're only copying the bottom half
          // of the element
          memcpy(result_biter, &(table_alias[geno_word & 3]), 4);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
      AppendU64(table_alias[cur_2geno], &result_biter);
      geno_word >>= 4;
    }
  }
}
#else
// plink2_glm_logistic does not guarantee 8 byte table alignment.
void GenoarrLookup16x4bx2(const uintptr_t* genoarr, const void* table16x4bx2, uint32_t sample_ct, void* __restrict result) {
  const uint32_t* table_u32 = S_CAST(const uint32_t*, table16x4bx2);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
          // not a regular CopyToUnaligned, we're only copying the bottom half
          // of the element
          memcpy(result_biter, &(table_u32[(geno_word & 3) * 2]), 4);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
      memcpy(result_biter, &(table_u32[cur_2geno * 2]), 8);
      result_biter += 8;
      geno_word >>= 4;
    }
  }
}
#endif

// this might be important for genovec -> AlleleCode expansion
void GenoarrLookup256x2bx4(const uintptr_t* genoarr, const void* table256x2bx4, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table256x2bx4);
  const unsigned char* genoarr_alias = DowncastKToUc(genoarr);
  unsigned char* resultb = S_CAST(unsigned char*, result);
  const uint32_t full_byte_ct = sample_ct / 4;
  for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
    CopyToUnalignedOffsetU64(resultb, &(table_alias[genoarr_alias[byte_idx]]), byte_idx);
  }
  const uint32_t remainder = sample_ct % 4;
  if (remainder) {
    unsigned char* result_last = &(resultb[full_byte_ct * sizeof(int64_t)]);
    uintptr_t geno_byte = genoarr_alias[full_byte_ct];
    for (uint32_t uii = 0; uii != remainder; ++uii) {
      CopyToUnalignedOffsetU16(result_last, DowncastKU64ToU16(&(table_alias[geno_byte & 3])), uii);
      geno_byte >>= 2;
    }
  }
}

#ifdef USE_SSE2
void GenoarrLookup4x16b(const uintptr_t* genoarr, const void* table4x16b, uint32_t sample_ct, void* result) {
  const __m128i* table_alias = S_CAST(const __m128i*, table4x16b);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
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
      _mm_storeu_si128(R_CAST(__m128i*, result_biter), table_alias[geno_word & 3]);
      result_biter += 16;
      geno_word >>= 2;
    }
  }
}

void GenoarrLookup16x8bx2(const uintptr_t* genoarr, const void* table16x8bx2, uint32_t sample_ct, void* __restrict result) {
  const __m128i* table_alias = S_CAST(const __m128i*, table16x8bx2);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
          memcpy(result_biter, &(table_alias[geno_word & 3]), 8);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
      _mm_storeu_si128(R_CAST(__m128i*, result_biter), table_alias[cur_2geno]);
      result_biter += 16;
      geno_word >>= 4;
    }
  }
}

void GenoarrLookup256x4bx4(const uintptr_t* genoarr, const void* table256x4bx4, uint32_t sample_ct, void* __restrict result) {
  const __m128i* table_alias = S_CAST(const __m128i*, table256x4bx4);
  const unsigned char* genoarr_alias = DowncastKToUc(genoarr);
  unsigned char* resultb = S_CAST(unsigned char*, result);
  const uint32_t full_byte_ct = sample_ct / 4;
  for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
    _mm_storeu_si128(R_CAST(__m128i*, &(resultb[byte_idx * 16])), table_alias[genoarr_alias[byte_idx]]);
  }
  const uint32_t remainder = sample_ct % 4;
  if (remainder) {
    unsigned char* result_last = &(resultb[full_byte_ct * 16]);
    uintptr_t geno_byte = genoarr_alias[full_byte_ct];
    for (uint32_t uii = 0; uii != remainder; ++uii) {
      CopyToUnalignedOffsetU32(result_last, R_CAST(const uint32_t*, &(table_alias[geno_byte & 3])), uii);
      geno_byte >>= 2;
    }
  }
}
#else // !USE_SSE2
void GenoarrLookup4x16b(const uintptr_t* genoarr, const void* table4x16b, uint32_t sample_ct, void* result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table4x16b);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
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
      memcpy(result_biter, &(table_alias[(geno_word & 3) * 2]), 16);
      result_biter += 16;
      geno_word >>= 2;
    }
  }
}

void GenoarrLookup16x8bx2(const uintptr_t* genoarr, const void* table16x8bx2, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table16x8bx2);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        if (sample_ct % 2) {
          CopyToUnalignedU64(result_biter, &(table_alias[(geno_word & 3) * 2]));
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      const uintptr_t cur_2geno = geno_word & 15;
      memcpy(result_biter, &(table_alias[cur_2geno * 2]), 16);
      result_biter += 16;
      geno_word >>= 4;
    }
  }
}

void GenoarrLookup256x4bx4(const uintptr_t* genoarr, const void* table256x4bx4, uint32_t sample_ct, void* __restrict result) {
  const uint32_t* table_alias = S_CAST(const uint32_t*, table256x4bx4);
  const unsigned char* genoarr_alias = DowncastKToUc(genoarr);
  unsigned char* resultb = S_CAST(unsigned char*, result);
  const uint32_t full_byte_ct = sample_ct / 4;
  for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
    memcpy(&(resultb[byte_idx * 16]), &(table_alias[genoarr_alias[byte_idx] * 4]), 16);
  }
  const uint32_t remainder = sample_ct % 4;
  if (remainder) {
    unsigned char* result_last = &(resultb[full_byte_ct * 16]);
    uintptr_t geno_byte = genoarr_alias[full_byte_ct];
    for (uint32_t uii = 0; uii != remainder; ++uii) {
      CopyToUnalignedOffsetU32(result_last, &(table_alias[(geno_byte & 3) * 4]), uii);
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

void InitLookup256x1bx4(void* table256x1bx4) {
  unsigned char* table_iter = S_CAST(unsigned char*, table256x1bx4);
  unsigned char vals[4];
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
#ifdef USE_SHUFFLE8
  // first vector: 0, 1, 2, 3
  // second vector: 0, 4, 8, 12
  uint32_t* final_u32s = R_CAST(uint32_t*, table_iter);
  uint32_t val0123;
  memcpy(&val0123, vals, 4);
  for (uint32_t uii = 0; uii != kInt32PerVec; ++uii) {
    final_u32s[uii] = val0123;
  }
  final_u32s[kInt32PerVec] = S_CAST(uint32_t, vals[0]) * 0x1010101;
  final_u32s[kInt32PerVec + 1] = S_CAST(uint32_t, vals[1]) * 0x1010101;
  final_u32s[kInt32PerVec + 2] = S_CAST(uint32_t, vals[2]) * 0x1010101;
  final_u32s[kInt32PerVec + 3] = S_CAST(uint32_t, vals[3]) * 0x1010101;
#  ifdef USE_AVX2
  memcpy(&(final_u32s[12]), &(final_u32s[8]), 16);
#  endif
#endif
}

void InitLookup256x2bx4(void* table256x2bx4) {
  uint16_t* table_iter = S_CAST(uint16_t*, table256x2bx4);
  uint16_t vals[4];
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
  const Halfword* phasepresent_alias = DowncastKWToHW(phasepresent);
  const Halfword* phaseinfo_alias = DowncastKWToHW(phaseinfo);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
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
          memcpy(result_biter, &(table_alias[cur_idx]), 4);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    phasepresent_hw_shifted = phasepresent_alias[widx];
    if (!phasepresent_hw_shifted) {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        AppendU64(table_alias[geno_word & 15], &result_biter);
        geno_word >>= 4;
      }
    } else {
      phasepresent_hw_shifted = phasepresent_hw_shifted << 4;
      phaseinfo_hw_shifted = phaseinfo_alias[widx];

      // note that this must be on a separate line (or we have to static_cast)
      phaseinfo_hw_shifted = phaseinfo_hw_shifted << 1;

      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uintptr_t cur_idx = ((geno_word & 15) | (phasepresent_hw_shifted & 48)) ^ (phaseinfo_hw_shifted & 6);
        AppendU64(table_alias[cur_idx], &result_biter);
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

#ifdef USE_SSE2
void PhaseLookup8b(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const void* table56x8bx2, uint32_t sample_ct, void* __restrict result) {
  const __m128i* table_alias = S_CAST(const __m128i*, table56x8bx2);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const Halfword* phasepresent_alias = DowncastKWToHW(phasepresent);
  const Halfword* phaseinfo_alias = DowncastKWToHW(phaseinfo);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  uintptr_t phasepresent_hw_shifted = 0;
  uintptr_t phaseinfo_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word & 3);
          if (phasepresent_hw_shifted & 16) {
            cur_idx ^= 16 | (phaseinfo_hw_shifted & 2);
          }
          memcpy(result_biter, &(table_alias[cur_idx]), 8);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    phasepresent_hw_shifted = phasepresent_alias[widx];
    if (!phasepresent_hw_shifted) {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        _mm_storeu_si128(R_CAST(__m128i*, result_biter), table_alias[geno_word & 15]);
        result_biter += 16;
        geno_word >>= 4;
      }
    } else {
      phasepresent_hw_shifted = phasepresent_hw_shifted << 4;
      phaseinfo_hw_shifted = phaseinfo_alias[widx];
      phaseinfo_hw_shifted = phaseinfo_hw_shifted << 1;
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uintptr_t cur_idx = ((geno_word & 15) | (phasepresent_hw_shifted & 48)) ^ (phaseinfo_hw_shifted & 6);
        _mm_storeu_si128(R_CAST(__m128i*, result_biter), table_alias[cur_idx]);
        result_biter += 16;
        geno_word >>= 4;
        phasepresent_hw_shifted >>= 2;
        phaseinfo_hw_shifted >>= 2;
      }
    }
  }
}
#else // !USE_SSE2
void PhaseLookup8b(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const void* table56x8bx2, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table56x8bx2);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const Halfword* phasepresent_alias = DowncastKWToHW(phasepresent);
  const Halfword* phaseinfo_alias = DowncastKWToHW(phaseinfo);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  uintptr_t phasepresent_hw_shifted = 0;
  uintptr_t phaseinfo_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word & 3);
          if (phasepresent_hw_shifted & 16) {
            cur_idx ^= 16 | (phaseinfo_hw_shifted & 2);
          }
          CopyToUnalignedU64(result_biter, &(table_alias[cur_idx * 2]));
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    phasepresent_hw_shifted = phasepresent_alias[widx];
    if (!phasepresent_hw_shifted) {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        memcpy(result_biter, &(table_alias[(geno_word & 15) * 2]), 16);
        result_biter = &(result_biter[16]);
        geno_word >>= 4;
      }
    } else {
      phasepresent_hw_shifted = phasepresent_hw_shifted << 4;
      phaseinfo_hw_shifted = phaseinfo_alias[widx];
      phaseinfo_hw_shifted = phaseinfo_hw_shifted << 1;
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uintptr_t cur_idx = ((geno_word & 15) | (phasepresent_hw_shifted & 48)) ^ (phaseinfo_hw_shifted & 6);
        memcpy(result_biter, &(table_alias[cur_idx * 2]), 16);
        // bugfix (12 Jun 2023): forgot this
        result_biter = &(result_biter[16]);

        geno_word >>= 4;
        phasepresent_hw_shifted >>= 2;
        phaseinfo_hw_shifted >>= 2;
      }
    }
  }
}
#endif

void InitPhaseLookup8b(void* table56x8bx2) {
  uint64_t* table_iter = S_CAST(uint64_t*, table56x8bx2);
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
    const uint64_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [16][0]..[31][1]: bit 4 is set
  // low bits must be 01 or 11
  const uint64_t val_phaseinfo0 = table_iter[2];
  table_iter[3] = vals[0];
  const uint64_t val_phaseinfo1 = table_iter[6];
  table_iter[7] = vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint64_t cur_high = vals[high_idx];
    table_iter[2] = val_phaseinfo0;
    table_iter[3] = cur_high;
    table_iter[6] = val_phaseinfo1;
    table_iter[7] = cur_high;
    table_iter = &(table_iter[8]);
  }
  // [32][0]..[39][1]: bit 5 set, bit 4 unset
  // high bits must be 00 or 01
  for (uint32_t high_idx = 0; high_idx != 2; ++high_idx) {
    const uint64_t cur_high = high_idx? val_phaseinfo0 : val_phaseinfo1;
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  table_iter = &(table_iter[16]);
  // [48][0]..[55][1]: bits 4 and 5 set
  for (uint32_t high_idx = 0; high_idx != 2; ++high_idx) {
    const uint64_t cur_high = high_idx? val_phaseinfo0 : val_phaseinfo1;
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
  const Halfword* phasepresent_alias = DowncastKWToHW(phasepresent);
  const Halfword* phaseinfo_alias = DowncastKWToHW(phaseinfo);
  const Halfword* sex_male_alias = DowncastKWToHW(sex_male);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word_xored = 0;
  uintptr_t male_or_phasepresent_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word_xored & 3) | (male_or_phasepresent_hw_shifted & 16);
          memcpy(result_biter, &(table_alias[cur_idx]), 4);
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
        AppendU64(table_alias[(geno_word_xored & 15) | (male_or_phasepresent_hw_shifted & 48)], &result_biter);
        geno_word_xored >>= 4;
        male_or_phasepresent_hw_shifted >>= 2;
      }
    } else {
      geno_word_xored ^= UnpackHalfwordToWordShift1(phaseinfo_alias[widx]);
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uintptr_t cur_idx = (geno_word_xored & 15) | (male_or_phasepresent_hw_shifted & 48);
        AppendU64(table_alias[cur_idx], &result_biter);
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
  const Halfword* sex_male_alias = DowncastKWToHW(sex_male);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  uintptr_t male_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word & 3) | (male_hw_shifted & 16);
          memcpy(result_biter, &(table_alias[cur_idx]), 4);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    male_hw_shifted = sex_male_alias[widx];
    male_hw_shifted <<= 4;
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      AppendU64(table_alias[(geno_word & 15) | (male_hw_shifted & 48)], &result_biter);
      geno_word >>= 4;
      male_hw_shifted >>= 2;
    }
  }
}

void InitPhaseXNohhLookup8b(void* table64x8bx2) {
  uint64_t* table_iter = S_CAST(uint64_t*, table64x8bx2);
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
    const uint64_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [16][0]..[31][1]: bit 4 is set
  uint64_t male_or_phasepresent_vals[4];
  for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
    male_or_phasepresent_vals[low_idx] = *table_iter++;
    *table_iter++ = vals[0];
  }
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint64_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = male_or_phasepresent_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [32][0]..[47][1]: bit 5 set, bit 4 unset
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    const uint64_t cur_high = male_or_phasepresent_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [48][0]..[63][1]: bits 4 and 5 set
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    const uint64_t cur_high = male_or_phasepresent_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = male_or_phasepresent_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
}

#ifdef USE_SSE2
void GenoarrSexLookup8b(const uintptr_t* genoarr, const uintptr_t* sex_male, const void* table64x8bx2, uint32_t sample_ct, void* result) {
  const __m128i* table_alias = S_CAST(const __m128i*, table64x8bx2);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const Halfword* sex_male_alias = DowncastKWToHW(sex_male);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  uintptr_t male_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word & 3) | (male_hw_shifted & 16);
          memcpy(result_biter, &(table_alias[cur_idx]), 8);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    male_hw_shifted = sex_male_alias[widx];
    male_hw_shifted <<= 4;
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      _mm_storeu_si128(R_CAST(__m128i*, result_biter), table_alias[(geno_word & 15) | (male_hw_shifted & 48)]);
      result_biter += 16;
      geno_word >>= 4;
      male_hw_shifted >>= 2;
    }
  }
}
#else // !USE_SSE2
void GenoarrSexLookup8b(const uintptr_t* genoarr, const uintptr_t* sex_male, const void* table64x8bx2, uint32_t sample_ct, void* result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table64x8bx2);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const Halfword* sex_male_alias = DowncastKWToHW(sex_male);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  uintptr_t male_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          const uintptr_t cur_idx = (geno_word & 3) | (male_hw_shifted & 16);
          memcpy(result_biter, &(table_alias[cur_idx * 2]), 8);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    male_hw_shifted = sex_male_alias[widx];
    male_hw_shifted <<= 4;
    for (uint32_t uii = 0; uii != loop_len; ++uii) {
      memcpy(result_biter, &(table_alias[((geno_word & 15) | (male_hw_shifted & 48)) * 2]), 16);
      result_biter = &(result_biter[16]);
      geno_word >>= 4;
      male_hw_shifted >>= 2;
    }
  }
}
#endif

void VcfPhaseLookup4b(const uintptr_t* genoarr, const uintptr_t* cur_phased, const uintptr_t* phaseinfo, const void* table246x4bx2, uint32_t sample_ct, void* __restrict result) {
  const uint64_t* table_alias = S_CAST(const uint64_t*, table246x4bx2);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const Halfword* cur_phased_alias = DowncastKWToHW(cur_phased);
  const Halfword* phaseinfo_alias = DowncastKWToHW(phaseinfo);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  uintptr_t cur_phased_hw_shifted = 0;
  uintptr_t phaseinfo_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word & 3) | (cur_phased_hw_shifted & 16) | (phaseinfo_hw_shifted & 64);
          memcpy(result_biter, &(table_alias[cur_idx]), 4);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    cur_phased_hw_shifted = cur_phased_alias[widx];
    if (!cur_phased_hw_shifted) {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        AppendU64(table_alias[geno_word & 15], &result_biter);
        geno_word >>= 4;
      }
    } else {
      cur_phased_hw_shifted = cur_phased_hw_shifted << 4;
      phaseinfo_hw_shifted = phaseinfo_alias[widx];

      // note that this must be on a separate line (or we have to static_cast)
      phaseinfo_hw_shifted = phaseinfo_hw_shifted << 6;

      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uintptr_t cur_idx = (geno_word & 15) | (cur_phased_hw_shifted & 48) | (phaseinfo_hw_shifted & 192);
        AppendU64(table_alias[cur_idx], &result_biter);
        geno_word >>= 4;
        cur_phased_hw_shifted >>= 2;
        phaseinfo_hw_shifted >>= 2;
      }
    }
  }
}

void InitVcfPhaseLookup4b(void* table246x4bx2) {
  uint32_t* table_iter = S_CAST(uint32_t*, table246x4bx2);
  uint32_t unphased_vals[4];
  unphased_vals[0] = table_iter[0];
  table_iter[1] = unphased_vals[0];
  unphased_vals[1] = table_iter[2];
  table_iter[3] = unphased_vals[0];
  unphased_vals[2] = table_iter[4];
  table_iter[5] = unphased_vals[0];
  unphased_vals[3] = table_iter[6];
  table_iter[7] = unphased_vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = unphased_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = unphased_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [16][0]..[31][1]: first entry is phased and unflipped, second is unphased
  uint32_t phased_unflipped_vals[4];
  phased_unflipped_vals[0] = table_iter[0];
  table_iter[1] = unphased_vals[0];
  phased_unflipped_vals[1] = table_iter[2];
  table_iter[3] = unphased_vals[0];
  phased_unflipped_vals[2] = table_iter[4];
  table_iter[5] = unphased_vals[0];
  phased_unflipped_vals[3] = table_iter[6];
  table_iter[7] = unphased_vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = unphased_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = phased_unflipped_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [32][0]..[63][1]: second entry is phased and unflipped
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = phased_unflipped_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = unphased_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = phased_unflipped_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = phased_unflipped_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [64][0]..[79][1] should be impossible
  table_iter = &(table_iter[32]);
  // [80][0]..[95][1]: first entry is phased and flipped, second is unphased
  // genotype must be 01
  const uint32_t phased_flipped_01 = table_iter[2];
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    table_iter[2] = phased_flipped_01;
    table_iter[3] = unphased_vals[high_idx];
    table_iter = &(table_iter[8]);
  }
  // [96][0]..[111][1] should be impossible
  table_iter = &(table_iter[32]);
  // [112][0]..[127][1]: first entry phased-flipped, second phased-unflipped
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    table_iter[2] = phased_flipped_01;
    table_iter[3] = phased_unflipped_vals[high_idx];
    table_iter = &(table_iter[8]);
  }
  // [128][0]..[163][1] should be impossible
  table_iter = &(table_iter[72]);
  // [164][0]..[167][1]: second entry phased-flipped, first entry unphased
  for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
    *table_iter++ = unphased_vals[low_idx];
    *table_iter++ = phased_flipped_01;
  }
  // [168][0]..[179][1] should be impossible
  table_iter = &(table_iter[24]);
  // [180][0]..[183][1]: second entry phased-flipped, first phased-unflipped
  for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
    *table_iter++ = phased_unflipped_vals[low_idx];
    *table_iter++ = phased_flipped_01;
  }
  // [184][0]..[244][1] should be impossible
  // [245][0]..[245][1]: both phased-flipped
  table_iter[122] = phased_flipped_01;
  table_iter[123] = phased_flipped_01;
}

void VcfPhaseLookup2b(const uintptr_t* genoarr, const uintptr_t* cur_phased, const uintptr_t* phaseinfo, const void* table246x2bx2, uint32_t sample_ct, void* __restrict result) {
  const uint32_t* table_alias = S_CAST(const uint32_t*, table246x2bx2);
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const Halfword* cur_phased_alias = DowncastKWToHW(cur_phased);
  const Halfword* phaseinfo_alias = DowncastKWToHW(phaseinfo);
  unsigned char* result_biter = S_CAST(unsigned char*, result);
  uint32_t loop_len = kBitsPerWordD4;
  uintptr_t geno_word = 0;
  uintptr_t cur_phased_hw_shifted = 0;
  uintptr_t phaseinfo_hw_shifted = 0;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        if (sample_ct % 2) {
          uintptr_t cur_idx = (geno_word & 3) | (cur_phased_hw_shifted & 16) | (phaseinfo_hw_shifted & 64);
          memcpy(result_biter, &(table_alias[cur_idx]), 2);
        }
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2) / 2;
    }
    geno_word = genoarr[widx];
    cur_phased_hw_shifted = cur_phased_alias[widx];
    if (!cur_phased_hw_shifted) {
      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        AppendU32(table_alias[geno_word & 15], &result_biter);
        geno_word >>= 4;
      }
    } else {
      cur_phased_hw_shifted = cur_phased_hw_shifted << 4;
      phaseinfo_hw_shifted = phaseinfo_alias[widx];

      // note that this must be on a separate line (or we have to static_cast)
      phaseinfo_hw_shifted = phaseinfo_hw_shifted << 6;

      for (uint32_t uii = 0; uii != loop_len; ++uii) {
        const uintptr_t cur_idx = (geno_word & 15) | (cur_phased_hw_shifted & 48) | (phaseinfo_hw_shifted & 192);
        AppendU32(table_alias[cur_idx], &result_biter);
        geno_word >>= 4;
        cur_phased_hw_shifted >>= 2;
        phaseinfo_hw_shifted >>= 2;
      }
    }
  }
}

void InitVcfPhaseLookup2b(void* table246x2bx2) {
  uint16_t* table_iter = S_CAST(uint16_t*, table246x2bx2);
  uint16_t unphased_vals[4];
  unphased_vals[0] = table_iter[0];
  table_iter[1] = unphased_vals[0];
  unphased_vals[1] = table_iter[2];
  table_iter[3] = unphased_vals[0];
  unphased_vals[2] = table_iter[4];
  table_iter[5] = unphased_vals[0];
  unphased_vals[3] = table_iter[6];
  table_iter[7] = unphased_vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = unphased_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = unphased_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [16][0]..[31][1]: first entry is phased and unflipped, second is unphased
  uint16_t phased_unflipped_vals[4];
  phased_unflipped_vals[0] = table_iter[0];
  table_iter[1] = unphased_vals[0];
  phased_unflipped_vals[1] = table_iter[2];
  table_iter[3] = unphased_vals[0];
  phased_unflipped_vals[2] = table_iter[4];
  table_iter[5] = unphased_vals[0];
  phased_unflipped_vals[3] = table_iter[6];
  table_iter[7] = unphased_vals[0];
  table_iter = &(table_iter[8]);
  for (uint32_t high_idx = 1; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = unphased_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = phased_unflipped_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [32][0]..[63][1]: second entry is phased and unflipped
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = phased_unflipped_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = unphased_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    const uint32_t cur_high = phased_unflipped_vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
      *table_iter++ = phased_unflipped_vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
  // [64][0]..[79][1] should be impossible
  table_iter = &(table_iter[32]);
  // [80][0]..[95][1]: first entry is phased and flipped, second is unphased
  // genotype must be 01
  const uint32_t phased_flipped_01 = table_iter[2];
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    table_iter[2] = phased_flipped_01;
    table_iter[3] = unphased_vals[high_idx];
    table_iter = &(table_iter[8]);
  }
  // [96][0]..[111][1] should be impossible
  table_iter = &(table_iter[32]);
  // [112][0]..[127][1]: first entry phased-flipped, second phased-unflipped
  for (uint32_t high_idx = 0; high_idx != 4; ++high_idx) {
    table_iter[2] = phased_flipped_01;
    table_iter[3] = phased_unflipped_vals[high_idx];
    table_iter = &(table_iter[8]);
  }
  // [128][0]..[163][1] should be impossible
  table_iter = &(table_iter[72]);
  // [164][0]..[167][1]: second entry phased-flipped, first entry unphased
  for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
    *table_iter++ = unphased_vals[low_idx];
    *table_iter++ = phased_flipped_01;
  }
  // [168][0]..[179][1] should be impossible
  table_iter = &(table_iter[24]);
  // [180][0]..[183][1]: second entry phased-flipped, first phased-unflipped
  for (uint32_t low_idx = 0; low_idx != 4; ++low_idx) {
    *table_iter++ = phased_unflipped_vals[low_idx];
    *table_iter++ = phased_flipped_01;
  }
  // [184][0]..[244][1] should be impossible
  // [245][0]..[245][1]: both phased-flipped
  table_iter[122] = phased_flipped_01;
  table_iter[123] = phased_flipped_01;
}


void ClearGenoarrMissing1bit8Unsafe(const uintptr_t* __restrict genoarr, uint32_t* subset_sizep, uintptr_t* __restrict subset, void* __restrict sparse_vals) {
  const uint32_t orig_subset_size = *subset_sizep;
  Halfword* subset_alias = DowncastWToHW(subset);
  uint32_t read_idx = 0;
  // deliberate overflow
  for (uint32_t read_widx = UINT32_MAX; ; ) {
    uint32_t subset_bits;
    do {
      subset_bits = subset_alias[++read_widx];
    } while (!subset_bits);
    uintptr_t detect_11 = genoarr[read_widx];
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
            detect_11 = genoarr[read_widx];
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
            detect_11 = genoarr[read_widx];
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

double u127prod_diff_d(uint64_t plus_term0, uint64_t plus_term1, uint64_t minus_term0, uint64_t minus_term1) {
  uint64_t plus_hi;
  const uint64_t plus_lo = multiply64to128(plus_term0, plus_term1, &plus_hi);
  uint64_t minus_hi;
  const uint64_t minus_lo = multiply64to128(minus_term0, minus_term1, &minus_hi);
  const uint64_t result_lo = plus_lo - minus_lo;
  const uint64_t result_hi = plus_hi - minus_hi - (plus_lo < minus_lo);
  return u127tod(result_hi, result_lo);
}

double i127prod_diff_d(uint64_t plus_term0, uint64_t plus_term1, uint64_t minus_term0, uint64_t minus_term1) {
  uint64_t plus_hi;
  const uint64_t plus_lo = multiply64to128(plus_term0, plus_term1, &plus_hi);
  uint64_t minus_hi;
  const uint64_t minus_lo = multiply64to128(minus_term0, minus_term1, &minus_hi);
  const uint64_t result_lo = plus_lo - minus_lo;
  const uint64_t result_hi = plus_hi - minus_hi - (plus_lo < minus_lo);
  if (!(result_hi & (1LLU << 63))) {
    return u127tod(result_hi, result_lo);
  } else {
    const uint64_t minus_result_lo = -result_lo;
    const uint64_t minus_result_hi = (-result_hi) - (result_lo != 0);
    return -u127tod(minus_result_hi, minus_result_lo);
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
  uint64_t ssq_sum_x2 = extra_phased_het_ct * 0x20000000LLU;
  uint64_t meansq_sum_lo = 0;
  uint64_t meansq_sum_hi = 0;
  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
    const uint64_t cur_allele_dosage = sums[allele_idx];
    const uint64_t cur_ssq_x2 = hap_ssqs_x2[allele_idx];
    ssq_sum_x2 += cur_ssq_x2;
    uint64_t incr_hi;
    uint64_t incr_lo = multiply64to128(cur_allele_dosage, cur_allele_dosage, &incr_hi);
    meansq_sum_lo += incr_lo;
    meansq_sum_hi += incr_hi + (meansq_sum_lo < incr_lo);
  }
  uint64_t prod_hi;
  uint64_t prod_lo = multiply64to128(ssq_sum_x2, nm_sample_ct, &prod_hi);
  const uint64_t observed_variance_times_2n_lo = prod_lo - meansq_sum_lo;
  const uint64_t observed_variance_times_2n_hi = prod_hi - meansq_sum_hi - (prod_lo < meansq_sum_lo);
  const double observed_variance_times_2n = u127tod(observed_variance_times_2n_hi, observed_variance_times_2n_lo);

  prod_lo = multiply64to128(nm_sample_ct * 0x40000000LLU, nm_sample_ct, &prod_hi);
  const uint64_t expected_variance_times_2n_lo = prod_lo - meansq_sum_lo;
  const uint64_t expected_variance_times_2n_hi = prod_hi - meansq_sum_hi - (prod_lo < meansq_sum_lo);
  const double expected_variance_times_2n = u127tod(expected_variance_times_2n_hi, expected_variance_times_2n_lo);
  return observed_variance_times_2n / expected_variance_times_2n;
}

void PgrDifflistToGenovecUnsafe(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uintptr_t difflist_common_geno, uint32_t sample_ct, uint32_t difflist_len, uintptr_t* __restrict genovec) {
  // Ok for trailing bits of raregeno to be nonzero.  Does not zero out
  // trailing bits of genovec.
  const uint32_t vec_ct = NypCtToVecCt(sample_ct);
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
      AssignNyparrEntry(cur_sample_idx, raregeno_word & 3, genovec);
      raregeno_word >>= 2;
    }
  }
}

const uint16_t kHcToAlleleCodes[1024] = QUAD_TABLE256(0, 0x100, 0x101, 0xffff);

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

uintptr_t PglComputeMaxAlleleCt(const uintptr_t* allele_idx_offsets, uint32_t variant_ct) {
  if ((!allele_idx_offsets) || (allele_idx_offsets[variant_ct] == 2 * variant_ct)) {
    return 2;
  }
  // todo: try vectorizing this
  uintptr_t max_allele_ct = 2;
  uintptr_t prev_offset = allele_idx_offsets[0];
  const uintptr_t* shifted_offsets = &(allele_idx_offsets[1]);
  for (uintptr_t uii = 0; uii != variant_ct; ++uii) {
    const uintptr_t cur_offset = shifted_offsets[uii];
    const uintptr_t cur_allele_ct = cur_offset - prev_offset;
    if (cur_allele_ct > max_allele_ct) {
      max_allele_ct = cur_allele_ct;
    }
    prev_offset = cur_offset;
  }
  return max_allele_ct;
}

// Ok for nybble_vvec to be unaligned.
uint32_t CountNybbleVec(const unsigned char* nybble_vvec_biter, uintptr_t nybble_word, uint32_t vec_ct) {
  const VecW m0 = vecw_setzero();
  const VecW alld15 = VCONST_W(kMask1111);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW xor_vvec = vecw_set1(nybble_word);
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
    const unsigned char* nybble_vvec_stop = &(nybble_vvec_biter[cur_incr * kBytesPerVec]);
    do {
      VecW loader = vecw_loadu(nybble_vvec_biter) ^ xor_vvec;
      nybble_vvec_biter += kBytesPerVec;
      // DetectAllZeroNybbles() followed by right-shift-3 is the same number of
      // operations, can see if that's any faster in practice
      loader = vecw_srli(loader, 1) | loader;
      loader = vecw_srli(loader, 2) | loader;
      inner_acc = inner_acc + vecw_and_notfirst(loader, alld15);
    } while (nybble_vvec_biter < nybble_vvec_stop);
    inner_acc = (inner_acc & m4) + (vecw_srli(inner_acc, 4) & m4);
    acc = acc + prev_sad_result;
    prev_sad_result = vecw_bytesum(inner_acc, m0);
  }
}

uint32_t CountNybble(const void* nybblearr, uintptr_t nybble_word, uintptr_t nybble_ct) {
  const unsigned char* nybblearr_uc = S_CAST(const unsigned char*, nybblearr);
  const uint32_t fullword_ct = nybble_ct / kBitsPerWordD4;
  uint32_t tot = CountNybbleVec(nybblearr_uc, nybble_word, fullword_ct / kWordsPerVec);
#ifdef __LP64__
  for (uint32_t word_idx = RoundDownPow2(fullword_ct, kWordsPerVec); word_idx != fullword_ct; ++word_idx) {
    uintptr_t cur_word;
    CopyFromUnalignedOffsetW(&cur_word, nybblearr_uc, word_idx);
    cur_word ^= nybble_word;
    cur_word = cur_word | (cur_word >> 1);
    cur_word = cur_word | (cur_word >> 2);
    tot += Popcount0001Word((~cur_word) & kMask1111);
  }
#endif
  const uint32_t trailing_nybble_ct = nybble_ct % kBitsPerWordD4;
  if (trailing_nybble_ct) {
    const uint32_t trailing_byte_ct = DivUp(trailing_nybble_ct, (CHAR_BIT / 4));
    uintptr_t cur_word = SubwordLoad(&(nybblearr_uc[fullword_ct * kBytesPerWord]), trailing_byte_ct) ^ nybble_word;
    cur_word = cur_word | (cur_word >> 1);
    cur_word = cur_word | (cur_word >> 2);
    cur_word = bzhi((~cur_word) & kMask1111, trailing_nybble_ct * 4);
#if defined(USE_SSE42) || !defined(__LP64__)
    tot += Popcount0001Word(cur_word);
#else
    // minor optimization, can't overflow
    tot += (cur_word * kMask1111) >> 60;
#endif
  }
  return tot;
}

uint64_t PglHeaderBaseEndOffset(uint32_t variant_ct, uintptr_t vrec_len_byte_ct, uint32_t phase_or_dosage_present, uint32_t explicit_nonref_flags) {
  const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
  uint64_t offset = 12 + vblock_ct * sizeof(int64_t) + variant_ct * vrec_len_byte_ct;
  if (phase_or_dosage_present) {
    // 8-bit vrtypes
    offset += variant_ct;
  } else {
    // 4-bit vrtypes
    offset += DivUp(variant_ct, 2);
  }
  if (explicit_nonref_flags) {
    offset += DivUp(variant_ct, CHAR_BIT);
  }
  return offset;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
