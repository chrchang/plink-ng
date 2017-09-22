// This library is part of PLINK 2, copyright (C) 2005-2017 Shaun Purcell,
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
boolerr_t pgl_malloc(uintptr_t size, void* pp) {
  *((unsigned char**)pp) = (unsigned char*)malloc(size);
  if (*((unsigned char**)pp)) {
    return 0;
  }
  g_failed_alloc_attempt_size = size;
  return 1;
}
#endif

interr_t fwrite_checked(const void* buf, uintptr_t len, FILE* outfile) {
  while (len > kMaxBytesPerIO) {
    // OS X can't perform 2GB+ writes
    // typical disk block size is 4kb, so 0x7ffff000 is the largest sensible
    // write size
    fwrite(buf, kMaxBytesPerIO, 1, outfile);
    buf = &(((const unsigned char*)buf)[kMaxBytesPerIO]);
    len -= kMaxBytesPerIO;
  }
  fwrite(buf, len, 1, outfile);
  return ferror(outfile);
}

interr_t fread_checked2(void* buf, uintptr_t len, FILE* infile, uintptr_t* bytes_read_ptr) {
  uintptr_t bytes_read = 0;
  while (len > kMaxBytesPerIO) {
    const uintptr_t cur_bytes_read = fread(buf, 1, kMaxBytesPerIO, infile);
    bytes_read += cur_bytes_read;
    if (cur_bytes_read != kMaxBytesPerIO) {
      *bytes_read_ptr = bytes_read;
      return ferror(infile);
    }
    buf = &(((char*)buf)[kMaxBytesPerIO]);
    len -= kMaxBytesPerIO;
  }
  bytes_read += fread(buf, 1, len, infile);
  *bytes_read_ptr = bytes_read;
  return ferror(infile);
}

#ifdef __LP64__
static inline boolerr_t scan_uint_capped_finish(const char* ss, uint64_t cap, uint32_t* valp) {
  uint64_t val = *valp;
  while (1) {
    // a little bit of unrolling seems to help
    const uint64_t cur_digit = (uint64_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit >= 10) {
      break;
    }
    // val = val * 10 + cur_digit;
    const uint64_t cur_digit2 = (uint64_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit2 >= 10) {
      val = val * 10 + cur_digit;
      if (val > cap) {
	return 1;
      }
      break;
    }
    val = val * 100 + cur_digit * 10 + cur_digit2;
    if (val > cap) {
      return 1;
    }
  }
  *valp = (uint32_t)val;
  return 0;
}

boolerr_t scan_posint_capped(const char* ss, uint64_t cap, uint32_t* valp) {
  // '0' has ascii code 48
  assert(((unsigned char)ss[0]) > 32);
  *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
  if (*valp >= 10) {
    // permit leading '+' (ascii 43), but not '++' or '+-'
    if (*valp != 0xfffffffbU) {
      return 1;
    }
    *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (*valp >= 10) {
      return 1;
    }
  }
  while (!(*valp)) {
    *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
    if ((*valp) >= 10) {
      return 1;
    }
  }
  return scan_uint_capped_finish(ss, cap, valp);
}

boolerr_t scan_uint_capped(const char* ss, uint64_t cap, uint32_t* valp) {
  // Reads an integer in [0, cap].  Assumes first character is nonspace. 
  assert(((unsigned char)ss[0]) > 32);
  uint32_t val = (uint32_t)((unsigned char)(*ss++)) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      // '-' has ascii code 45, so unsigned 45 - 48 = 0xfffffffdU
      if ((val != 0xfffffffdU) || (*ss != '0')) {
	return 1;
      }
      // accept "-0", "-00", etc.
      while (*(++ss) == '0');
      *valp = 0;
      return ((uint32_t)((unsigned char)(*ss)) - 48) < 10;      
    }
    // accept leading '+'
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  *valp = val;
  return scan_uint_capped_finish(ss, cap, valp);
}

boolerr_t scan_int_abs_bounded(const char* ss, uint64_t bound, int32_t* valp) {
  // Reads an integer in [-bound, bound].  Assumes first character is nonspace.
  assert(((unsigned char)ss[0]) > 32);
  *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
  int32_t sign = 1;
  if (((uint32_t)*valp) >= 10) {
    if (*valp == -3) {
      sign = -1;
    } else if (*valp != -5) {
      return 1;
    }
    *valp = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (((uint32_t)*valp) >= 10) {
      return 1;
    }
  }
  if (scan_uint_capped_finish(ss, bound, (uint32_t*)valp)) {
    return 1;
  }
  *valp *= sign;
  return 0;
}
#else // not __LP64__
boolerr_t scan_posint_capped32(const char* ss, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp) {
  // '0' has ascii code 48
  assert(((unsigned char)ss[0]) > 32);
  uint32_t val = (uint32_t)((unsigned char)(*ss++)) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      return 1;
    }
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (!val) {
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (1) {
    const uint32_t cur_digit = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      return 0;
    }
    // avoid integer overflow in middle of computation
    if ((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}

boolerr_t scan_uint_capped32(const char* ss, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp) {
  // Reads an integer in [0, cap].  Assumes first character is nonspace. 
  assert(((unsigned char)ss[0]) > 32);
  uint32_t val = (uint32_t)((unsigned char)(*ss++)) - 48;
  if (val >= 10) {
    if (val != 0xfffffffbU) {
      if ((val != 0xfffffffdU) || (*ss != '0')) {
	return 1;
      }
      while (*(++ss) == '0');
      *valp = 0;
      return ((uint32_t)((unsigned char)(*ss)) - 48) < 10;
    }
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (1) {
    const uint32_t cur_digit = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit >= 10) {
      *valp = val;
      return 0;
    }
    if ((val >= cap_div_10) && ((val > cap_div_10) || (cur_digit > cap_mod_10))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}

boolerr_t scan_int_abs_bounded32(const char* ss, uint32_t bound_div_10, uint32_t bound_mod_10, int32_t* valp) {
  // Reads an integer in [-bound, bound].  Assumes first character is nonspace.
  assert(((unsigned char)ss[0]) > 32);
  uint32_t val = (uint32_t)((unsigned char)(*ss++)) - 48;
  int32_t sign = 1;
  if (val >= 10) {
    if (val == 0xfffffffdU) {
      sign = -1;
    } else if (val != 0xfffffffbU) {
      return 1;
    }
    val = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (val >= 10) {
      return 1;
    }
  }
  while (1) {
    const uint32_t cur_digit = (uint32_t)((unsigned char)(*ss++)) - 48;
    if (cur_digit >= 10) {
      *valp = sign * ((int32_t)val);
      return 0;
    }
    if ((val >= bound_div_10) && ((val > bound_div_10) || (cur_digit > bound_mod_10))) {
      return 1;
    }
    val = val * 10 + cur_digit;
  }
}
#endif

boolerr_t aligned_malloc(uintptr_t size, uintptr_t alignment, void* aligned_pp) {
  // Assumes malloc returns word-aligned addresses.
  assert(alignment);
  assert(!(alignment % kBytesPerWord));
  uintptr_t malloc_addr;
  if (pgl_malloc(size + alignment, &malloc_addr)) {
    return 1;
  }
  assert(!(malloc_addr % kBytesPerWord));
  uintptr_t** casted_aligned_pp = (uintptr_t**)aligned_pp;
  *casted_aligned_pp = (uintptr_t*)round_down_pow2(malloc_addr + alignment, alignment);
  (*casted_aligned_pp)[-1] = malloc_addr;
  return 0;
}

void fill_all_bits(uintptr_t ct, uintptr_t* bitarr) {
  // leaves bits beyond the end unset
  // ok for ct == 0
  uintptr_t quotient = ct / kBitsPerWord;
  uintptr_t remainder = ct % kBitsPerWord;
  fill_ulong_one(quotient, bitarr);
  if (remainder) {
    bitarr[quotient] = (k1LU << remainder) - k1LU;
  }
}

void bitvec_and(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec AND arg_bitvec
#ifdef __LP64__
  vul_t* main_bitvvec_iter = (vul_t*)main_bitvec;
  const vul_t* arg_bitvvec_iter = (const vul_t*)arg_bitvec;
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
  #ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] &= arg_bitvec[base_idx];
    main_bitvec[base_idx + 1] &= arg_bitvec[base_idx + 1];
  }
  #endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] &= arg_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    main_bitvec[widx] &= arg_bitvec[widx];
  }
#endif
}

void bitvec_andnot(const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := main_bitvec ANDNOT exclude_bitvec
  // note that this is the reverse of the _mm_andnot() operand order
#ifdef __LP64__
  vul_t* main_bitvvec_iter = (vul_t*)main_bitvec;
  const vul_t* exclude_bitvvec_iter = (const vul_t*)exclude_bitvec;
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
  #ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] &= ~exclude_bitvec[base_idx];
    main_bitvec[base_idx + 1] &= ~exclude_bitvec[base_idx + 1];
  }
  #endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] &= ~exclude_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    main_bitvec[widx] &= ~exclude_bitvec[widx];
  }
#endif
}

uint32_t next_set_unsafe(const uintptr_t* bitarr, uint32_t loc) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / kBitsPerWord]);
  uintptr_t ulii = (*bitarr_ptr) >> (loc % kBitsPerWord);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bitarr_ptr);
  } while (!ulii);
  return (uint32_t)(((uintptr_t)(bitarr_ptr - bitarr)) * kBitsPerWord + CTZLU(ulii));
}

uint32_t next_unset_unsafe(const uintptr_t* bitarr, uint32_t loc) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / kBitsPerWord]);
  uintptr_t ulii = (~(*bitarr_ptr)) >> (loc % kBitsPerWord);
  if (ulii) {
    return loc + CTZLU(ulii);
  }
  do {
    ulii = *(++bitarr_ptr);
  } while (ulii == ~k0LU);
  return (uint32_t)(((uintptr_t)(bitarr_ptr - bitarr)) * kBitsPerWord + CTZLU(~ulii));
}

/*
uint32_t next_nonmissing_unsafe(const uintptr_t* genoarr, uint32_t loc) {
  const uintptr_t* genoarr_ptr = &(genoarr[loc / kBitsPerWordD2]);
  uintptr_t ulii = (~(*genoarr_ptr)) >> (2 * (loc % kBitsPerWordD2));
  if (ulii) {
    return loc + (CTZLU(ulii) / 2);
  }
  do {
    ulii = *(++genoarr_ptr);
  } while (ulii == ~k0LU);
  return ((uintptr_t)(genoarr_ptr - genoarr)) * kBitsPerWordD2 + (CTZLU(~ulii) / 2);
}
*/

uint32_t next_set(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil) {
  // safe version.
  const uintptr_t* bitarr_ptr = &(bitarr[loc / kBitsPerWord]);
  uintptr_t ulii = (*bitarr_ptr) >> (loc % kBitsPerWord);
  uint32_t rval;
  if (ulii) {
    rval = loc + CTZLU(ulii);
    return MINV(rval, ceil);
  }
  const uintptr_t* bitarr_last = &(bitarr[(ceil - 1) / kBitsPerWord]);
  do {
    if (bitarr_ptr >= bitarr_last) {
      return ceil;
    }
    ulii = *(++bitarr_ptr);
  } while (!ulii);
  rval = (uint32_t)(((uintptr_t)(bitarr_ptr - bitarr)) * kBitsPerWord + CTZLU(ulii));
  return MINV(rval, ceil);
}

uint32_t prev_set_unsafe(const uintptr_t* bitarr, uint32_t loc) {
  // unlike the next_{un}set family, this always returns a STRICTLY earlier
  // position
  const uintptr_t* bitarr_ptr = &(bitarr[loc / kBitsPerWord]);
  uint32_t remainder = loc % kBitsPerWord;
  uintptr_t ulii;
  if (remainder) {
    ulii = (*bitarr_ptr) & ((k1LU << remainder) - k1LU);
    if (ulii) {
      return (loc | (kBitsPerWord - 1)) - CLZLU(ulii);
    }
  }
  do {
    ulii = *(--bitarr_ptr);
  } while (!ulii);
  return (uint32_t)(((uintptr_t)(bitarr_ptr - bitarr)) * kBitsPerWord + kBitsPerWord - 1 - CLZLU(ulii));
}

void copy_bitarr_subset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t subset_size, uintptr_t* __restrict output_bitarr) {
  // could try exploiting _pext_u64() intrinsic, but probably not worthwhile
  // until 2020ish
  const uint32_t subset_size_lowbits = subset_size % kBitsPerWord;
  uintptr_t* output_bitarr_iter = output_bitarr;
  uintptr_t* output_bitarr_last = &(output_bitarr[subset_size / kBitsPerWord]);
  uintptr_t cur_output_word = 0;
  uint32_t read_widx = 0xffffffffU; // deliberate overflow
  uint32_t write_idx_lowbits = 0;
  while ((output_bitarr_iter != output_bitarr_last) || (write_idx_lowbits != subset_size_lowbits)) {
    uintptr_t cur_mask_word;
    // sparse subset_mask optimization
    // guaranteed to terminate since there's at least one more set bit
    do {
      cur_mask_word = subset_mask[++read_widx];
    } while (!cur_mask_word);
    uintptr_t cur_masked_input_word = raw_bitarr[read_widx] & cur_mask_word;
    const uint32_t cur_mask_popcount = popcount_long(cur_mask_word);
    uintptr_t subsetted_input_word = 0;
    if (cur_masked_input_word) {
      const uintptr_t cur_inv_mask = ~cur_mask_word;
      do {
	const uint32_t read_uidx_nz_start_lowbits = CTZLU(cur_masked_input_word);
	const uintptr_t cur_inv_mask_shifted = cur_inv_mask >> read_uidx_nz_start_lowbits;
	if (!cur_inv_mask_shifted) {
	  subsetted_input_word |= cur_masked_input_word >> (kBitsPerWord - cur_mask_popcount);
	  break;
	}
	const uint32_t cur_read_end = CTZLU(cur_inv_mask_shifted) + read_uidx_nz_start_lowbits;
	// this seems to optimize better than (k1LU << cur_read_end) - k1LU
	// todo: check if/when that's true elsewhere
        const uintptr_t lowmask = (~k0LU) >> (kBitsPerWord - cur_read_end);
	const uintptr_t bits_to_copy = cur_masked_input_word & lowmask;
	cur_masked_input_word -= bits_to_copy;
	// todo: check if a less-popcounty implementation should be used in
	// non-SSE4.2 case
	const uint32_t cur_write_end = popcount_long(cur_mask_word & lowmask);
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
uintptr_t popcount_vecs(const vul_t* bit_vvec, uintptr_t vec_ct) {
  // popcounts vptr[0..(vec_ct-1)].  Assumes vec_ct is a multiple of 3 (0 ok).
  assert(!(vec_ct % 3));
  const vul_t m1 = VCONST_UL(kMask5555);
  const vul_t m2 = VCONST_UL(kMask3333);
  const vul_t m4 = VCONST_UL(kMask0F0F);
  const vul_t m8 = VCONST_UL(kMask00FF);
  const vul_t* bit_vvec_iter = bit_vvec;
  uintptr_t tot = 0;
  while (1) {
    univec_t acc;
    acc.vi = vul_setzero();
    const vul_t* bit_vvec_stop;
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
      vul_t count1 = *bit_vvec_iter++;
      vul_t count2 = *bit_vvec_iter++;
      vul_t half1 = *bit_vvec_iter++;
      vul_t half2 = vul_rshift(half1, 1) & m1;
      half1 = half1 & m1;
      // Two bits can represent values from 0-3, so make each pair in count1
      // count2 store a partial bitcount covering themselves AND another bit
      // from elsewhere.
      count1 = count1 - (vul_rshift(count1, 1) & m1);
      count2 = count2 - (vul_rshift(count2, 1) & m1);
      count1 = count1 + half1;
      count2 = count2 + half2;
      // Four bits represent 0-15, so we can safely add four 0-3 partial
      // bitcounts together.
      count1 = (count1 & m2) + (vul_rshift(count1, 2) & m2);
      count1 = count1 + (count2 & m2) + (vul_rshift(count2, 2) & m2);
      // Accumulator stores sixteen 0-255 counts in parallel.
      // (32 in AVX2 case, 4 in 32-bit case)
      acc.vi = acc.vi + (count1 & m4) + (vul_rshift(count1, 4) & m4);
    } while (bit_vvec_iter < bit_vvec_stop);
    acc.vi = (acc.vi & m8) + (vul_rshift(acc.vi, 8) & m8);
    tot += univec_hsum_16bit(acc);
  }
}

uintptr_t popcount_bytes(const unsigned char* bitarr, uintptr_t byte_ct) {
  const uint32_t lead_byte_ct = ((uintptr_t)(-((uintptr_t)bitarr))) % kBytesPerVec;
  uintptr_t tot = 0;
  const uintptr_t* bitarr_iter;
  uint32_t trail_byte_ct;
  // bugfix: had wrong condition here
  if (byte_ct >= lead_byte_ct) {
#ifdef __LP64__
    const uint32_t word_rem = lead_byte_ct % kBytesPerWord;
    if (word_rem) {
      uintptr_t cur_word = 0;
      memcpy(&cur_word, bitarr, word_rem);
      tot = popcount_long(cur_word);
    }
    bitarr_iter = (const uintptr_t*)(&(bitarr[word_rem]));
    if (lead_byte_ct / kBytesPerWord) {
      tot += popcount_long(*bitarr_iter++);
    }
#else
    if (lead_byte_ct) {
      uintptr_t cur_word = 0;
      memcpy(&cur_word, bitarr, lead_byte_ct);
      tot = popcount_long(cur_word);
    }
    bitarr_iter = (const uintptr_t*)(&(bitarr[lead_byte_ct]));
#endif
    byte_ct -= lead_byte_ct;
    const uintptr_t word_ct = byte_ct / kBytesPerWord;
    // vec-alignment required here
    tot += popcount_longs(bitarr_iter, word_ct);
    bitarr_iter = &(bitarr_iter[word_ct]);
    trail_byte_ct = byte_ct % kBytesPerWord;
  } else {
    bitarr_iter = (const uintptr_t*)bitarr;
    // this may still be >= kBytesPerWord, so can't remove loop
    trail_byte_ct = (uint32_t)byte_ct;
  }
  while (1) {
    uintptr_t cur_word;
    if (trail_byte_ct < kBytesPerWord) {
      if (!trail_byte_ct) {
	return tot;
      }
      cur_word = 0;
      memcpy(&cur_word, bitarr_iter, trail_byte_ct);
      trail_byte_ct = 0;
    } else {
      cur_word = *bitarr_iter++;
      trail_byte_ct -= kBytesPerWord;
    }
    tot += popcount_long(cur_word);
  }
}

uintptr_t popcount_bytes_masked(const unsigned char* bitarr, const uintptr_t* mask_arr, uintptr_t byte_ct) {
  // could detect aligned case, but that shouldn't happen often enough?
  const uintptr_t word_ct = byte_ct / kBytesPerWord;
#ifdef USE_SSE42
  const uintptr_t* bitarr_alias = (const uintptr_t*)bitarr;
  uintptr_t tot = 0;
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    tot += popcount_long(bitarr_alias[widx] & mask_arr[widx]);
  }
  const uint32_t trail_byte_ct = byte_ct % kBytesPerWord;
  if (trail_byte_ct) {
    uintptr_t cur_word = 0;
    memcpy(&cur_word, &(bitarr_alias[word_ct]), trail_byte_ct);
    tot += popcount_long(cur_word & mask_arr[word_ct]);
  }
  return tot;
#else
  const uintptr_t* bitarr_iter = (const uintptr_t*)bitarr;
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

  #ifndef __LP64__
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
  #endif

    // 32-bit case: each 8-bit slot stores a number in 0..48.  Multiplying by
    // 0x01010101 is equivalent to the left-shifts and adds we need to sum
    // those four 8-bit numbers in the high-order slot.
    // 64-bit case: each 8-bit slot stores a number in 0..24.
    tot += (tmp_stor * kMask0101) >> (kBitsPerWord - 8);
  }
  uint32_t trail_byte_ct = (uint32_t)(byte_ct - (mainblock_word_ct * kBytesPerWord));
  while (1) {
    uintptr_t cur_word;
    if (trail_byte_ct < kBytesPerWord) {
      if (!trail_byte_ct) {
	return tot;
      }
      cur_word = 0;
      memcpy(&cur_word, bitarr_iter, trail_byte_ct);
      trail_byte_ct = 0;
    } else {
      cur_word = *bitarr_iter++;
      trail_byte_ct -= kBytesPerWord;
    }
    tot += popcount_long(cur_word & (*mask_arr_iter++));
  }
#endif
}

void fill_cumulative_popcounts(const uintptr_t* subset_mask, uint32_t word_ct, uint32_t* cumulative_popcounts) {
  assert(word_ct);
  const uint32_t word_ct_m1 = word_ct - 1;
  uint32_t cur_sum = 0;
  for (uint32_t widx = 0; widx < word_ct_m1; ++widx) {
    cumulative_popcounts[widx] = cur_sum;
    cur_sum += popcount_long(subset_mask[widx]);
  }
  cumulative_popcounts[word_ct_m1] = cur_sum;
}

void uidxs_to_idxs(const uintptr_t* subset_mask, const uint32_t* subset_cumulative_popcounts, const uint32_t idx_list_len, uint32_t* idx_list) {
  uint32_t* idx_list_end = &(idx_list[idx_list_len]);
  for (uint32_t* idx_list_iter = idx_list; idx_list_iter != idx_list_end; ++idx_list_iter) {
    *idx_list_iter = raw_to_subsetted_pos(subset_mask, subset_cumulative_popcounts, *idx_list_iter);
  }
}


static_assert(kPglBitTransposeBatch == ((uint32_t)kBitsPerCacheline), "transpose_bitblock() needs to be updated.");
#ifdef __LP64__
static_assert(kWordsPerVec == 2, "transpose_bitblock() needs to be updated.");
#else
static_assert(kWordsPerVec == 1, "transpose_bitblock() needs to be updated.");
#endif
void transpose_bitblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, vul_t* vecaligned_buf) {
  // buf must be vector-aligned and have size 64k
  const uint32_t initial_read_byte_ct = DIV_UP(write_batch_size, CHAR_BIT);
  // fold the first 6 shuffles into the initial ingestion loop
  const unsigned char* initial_read_iter = (const unsigned char*)read_iter;
  const unsigned char* initial_read_end = &(initial_read_iter[initial_read_byte_ct]);
  unsigned char* initial_target_iter = (unsigned char*)vecaligned_buf;
  const uint32_t read_byte_stride = read_ul_stride * kBytesPerWord;
  const uint32_t read_batch_rem = kBitsPerCacheline - read_batch_size;
  for (; initial_read_iter < initial_read_end; ++initial_read_iter) {
    const unsigned char* read_iter_tmp = initial_read_iter;
    for (uint32_t ujj = 0; ujj < read_batch_size; ++ujj) {
      *initial_target_iter++ = *read_iter_tmp;
      read_iter_tmp = &(read_iter_tmp[read_byte_stride]);
    }
    initial_target_iter = memseta(initial_target_iter, 0, read_batch_rem);
  }

  // third-to-last shuffle, 8 bit spacing -> 4
  const vul_t* source_iter = vecaligned_buf;
  uintptr_t* target_iter0 = (uintptr_t*)(&(vecaligned_buf[kPglBitTransposeBufwords / (2 * kWordsPerVec)]));
#ifdef __LP64__
  const vul_t m4 = VCONST_UL(kMask0F0F);
  const vul_t m8 = VCONST_UL(kMask00FF);
  const vul_t m16 = VCONST_UL(kMask0000FFFF);
#endif
  const uint32_t write_word_ct = BITCT_TO_WORDCT(read_batch_size);
  const uint32_t first_inner_loop_iter_ct = 4 * write_word_ct;
  uint32_t cur_write_skip = 4 * kWordsPerCacheline - first_inner_loop_iter_ct;
  // coincidentally, this also needs to run DIV_UP(write_batch_size, CHAR_BIT)
  // times
  for (uint32_t uii = 0; uii < initial_read_byte_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 4]);
    for (uint32_t ujj = 0; ujj < first_inner_loop_iter_ct; ++ujj) {
#ifdef __LP64__
      const vul_t loader = *source_iter++;
      vul_t target0 = loader & m4;
      vul_t target1 = (vul_rshift(loader, 4)) & m4;
      target0 = (target0 | (vul_rshift(target0, 4))) & m8;
      target1 = (target1 | (vul_rshift(target1, 4))) & m8;
      target0 = (target0 | (vul_rshift(target0, 8))) & m16;
      target1 = (target1 | (vul_rshift(target1, 8))) & m16;
      univec_t target0u;
      univec_t target1u;
      target0u.vi = target0 | (vul_rshift(target0, 16));
      target1u.vi = target1 | (vul_rshift(target1, 16));
      *target_iter0++ = ((uint32_t)target0u.u8[0]) | (target0u.u8[1] << 32);
      *target_iter1++ = ((uint32_t)target1u.u8[0]) | (target1u.u8[1] << 32);
#else
      const uintptr_t source_word_lo = (uintptr_t)(*source_iter++);
      const uintptr_t source_word_hi = (uintptr_t)(*source_iter++);
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
      *target_iter0++ = ((halfword_t)target_word0_lo) | (target_word0_hi << kBitsPerWordD2);
      *target_iter1++ = ((halfword_t)target_word1_lo) | (target_word1_hi << kBitsPerWordD2);
#endif
    }
#ifdef __LP64__
    source_iter = &(source_iter[cur_write_skip]);
#else
    source_iter = &(source_iter[2 * cur_write_skip]);
#endif
    target_iter0 = &(target_iter1[cur_write_skip]);
  }

  // second-to-last shuffle, 4 bit spacing -> 2
  source_iter = (&(vecaligned_buf[kPglBitTransposeBufwords / (2 * kWordsPerVec)]));
  target_iter0 = (uintptr_t*)vecaligned_buf;
#ifdef __LP64__
  const vul_t m2 = VCONST_UL(kMask3333);
#endif
  const uint32_t second_outer_loop_iter_ct = DIV_UP(write_batch_size, 4);
  const uint32_t second_inner_loop_iter_ct = 2 * write_word_ct;
  cur_write_skip = 2 * kWordsPerCacheline - second_inner_loop_iter_ct;
  for (uint32_t uii = 0; uii < second_outer_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[kWordsPerCacheline * 2]);
    for (uint32_t ujj = 0; ujj < second_inner_loop_iter_ct; ++ujj) {
#ifdef __LP64__
      // in AVX2 case, use write_dword_ct instead of write_word_ct, etc.
      const vul_t loader = *source_iter++;
      vul_t target0 = loader & m2;
      vul_t target1 = (vul_rshift(loader, 2)) & m2;
      target0 = (target0 | (vul_rshift(target0, 2))) & m4;
      target1 = (target1 | (vul_rshift(target1, 2))) & m4;
      target0 = (target0 | (vul_rshift(target0, 4))) & m8;
      target1 = (target1 | (vul_rshift(target1, 4))) & m8;
      target0 = (target0 | (vul_rshift(target0, 8))) & m16;
      target1 = (target1 | (vul_rshift(target1, 8))) & m16;
      univec_t target0u;
      univec_t target1u;
      target0u.vi = target0 | (vul_rshift(target0, 16));
      target1u.vi = target1 | (vul_rshift(target1, 16));
      *target_iter0++ = ((uint32_t)target0u.u8[0]) | (target0u.u8[1] << 32);
      *target_iter1++ = ((uint32_t)target1u.u8[0]) | (target1u.u8[1] << 32);
#else
      const uintptr_t source_word_lo = (uintptr_t)(*source_iter++);
      const uintptr_t source_word_hi = (uintptr_t)(*source_iter++);
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
      *target_iter0++ = ((halfword_t)target_word0_lo) | (target_word0_hi << kBitsPerWordD2);
      *target_iter1++ = ((halfword_t)target_word1_lo) | (target_word1_hi << kBitsPerWordD2);
#endif
    }
#ifdef __LP64__
    source_iter = &(source_iter[cur_write_skip]);
#else
    source_iter = &(source_iter[2 * cur_write_skip]);
#endif
    target_iter0 = &(target_iter1[cur_write_skip]);
  }
  // last shuffle, 2 bit spacing -> 1
  source_iter = vecaligned_buf;
  target_iter0 = write_iter;
#ifdef __LP64__
  const vul_t m1 = VCONST_UL(kMask5555);
#endif
  const uint32_t last_loop_iter_ct = DIV_UP(write_batch_size, 2);
  for (uint32_t uii = 0; uii < last_loop_iter_ct; ++uii) {
    uintptr_t* target_iter1 = &(target_iter0[write_ul_stride]);
    for (uint32_t ujj = 0; ujj < write_word_ct; ++ujj) {
#ifdef __LP64__
      // in AVX2 case, use write_dword_ct instead of write_word_ct, etc.
      const vul_t loader = *source_iter++;
      vul_t target0 = loader & m1;
      vul_t target1 = (vul_rshift(loader, 1)) & m1;
      target0 = (target0 | (vul_rshift(target0, 1))) & m2;
      target1 = (target1 | (vul_rshift(target1, 1))) & m2;
      target0 = (target0 | (vul_rshift(target0, 2))) & m4;
      target1 = (target1 | (vul_rshift(target1, 2))) & m4;
      target0 = (target0 | (vul_rshift(target0, 4))) & m8;
      target1 = (target1 | (vul_rshift(target1, 4))) & m8;
      target0 = (target0 | (vul_rshift(target0, 8))) & m16;
      target1 = (target1 | (vul_rshift(target1, 8))) & m16;
      univec_t target0u;
      univec_t target1u;
      target0u.vi = target0 | (vul_rshift(target0, 16));
      target1u.vi = target1 | (vul_rshift(target1, 16));
      target_iter0[ujj] = ((uint32_t)target0u.u8[0]) | (target0u.u8[1] << 32);
      target_iter1[ujj] = ((uint32_t)target1u.u8[0]) | (target1u.u8[1] << 32);
#else
      const uintptr_t source_word_lo = (uintptr_t)(*source_iter++);
      const uintptr_t source_word_hi = (uintptr_t)(*source_iter++);
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
      target_iter0[ujj] = ((halfword_t)target_word0_lo) | (target_word0_hi << kBitsPerWordD2);
      target_iter1[ujj] = ((halfword_t)target_word1_lo) | (target_word1_hi << kBitsPerWordD2);
#endif
    }
#ifdef __LP64__
    source_iter = &(source_iter[kWordsPerCacheline - write_word_ct]);
#else
    source_iter = &(source_iter[2 * (kWordsPerCacheline - write_word_ct)]);
#endif
    target_iter0 = &(target_iter1[write_ul_stride]);
  }
}

#ifdef __cplusplus
} // namespace plink2
#endif
