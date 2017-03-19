// This library is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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


#include "pgenlib_python_support.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void genoarr_to_bytes_minus9(const uintptr_t* genoarr, uint32_t sample_ct, int8_t* genobytes) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBytesPerWord;
  const quarterword_t* read_alias = (const quarterword_t*)genoarr;
  uintptr_t* write_walias = (uintptr_t*)genobytes;
  uint32_t widx = 0;
  while (1) {
    uintptr_t qw = read_alias[widx];
#ifdef __LP64__
    qw = (qw | (qw << 24)) & kMask000000FF;
#endif
    qw = (qw | (qw << 12)) & kMask000F;
    qw = (qw | (qw << 6)) & kMask0303;
    // now each byte is in {0, 1, 2, 3}.  Convert the 3s to -9s in a branchless
    // manner.
    // (-9) - 3 = -12, which is represented as 244 in a uint8_t
    const uintptr_t geno_missing = qw & (qw >> 1) & kMask0101;
    qw += geno_missing * 244;
    if (widx == word_ct_m1) {
      memcpy(&(write_walias[widx]), &qw, MOD_NZ(sample_ct, kBytesPerWord));
      return;
    }
    write_walias[widx++] = qw;
  }
}

// could have a size-16 lookup table in 64-bit builds, etc.
static const int32_t geno_to_int32[4] = {0, 1, 2, -9};

void genoarr_to_int32_minus9(const uintptr_t* genoarr, uint32_t sample_ct, int32_t* geno_int32) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  int32_t* write_iter = geno_int32;
  uint32_t subgroup_len = kBitsPerWordD2;
  uint32_t widx = 0;
  while (1) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
	return;
      }
      subgroup_len = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
      *write_iter++ = geno_to_int32[geno_word & 3];
      geno_word >>= 2;
    }
    ++widx;
  }
}

static const int64_t geno_to_int64[4] = {0, 1, 2, -9};

void genoarr_to_int64_minus9(const uintptr_t* genoarr, uint32_t sample_ct, int64_t* geno_int64) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  int64_t* write_iter = geno_int64;
  uint32_t subgroup_len = kBitsPerWordD2;
  uint32_t widx = 0;
  while (1) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
	return;
      }
      subgroup_len = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
      *write_iter++ = geno_to_int64[geno_word & 3];
      geno_word >>= 2;
    }
    ++widx;
  }
}

// missing = -9
static const uint64_t geno_to_intcode_pair[4] = {0, 0x100000000LLU, 0x100000001LLU, 0xfffffff7fffffff7LLU};

void genoarr_to_allele_codes(const uintptr_t* genoarr, uint32_t sample_ct, int32_t* allele_codes) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint64_t* write_iter = (uint64_t*)allele_codes;
  uint32_t subgroup_len = kBitsPerWordD2;
  uint32_t widx = 0;
  while (1) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
	return;
      }
      subgroup_len = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
      *write_iter++ = geno_to_intcode_pair[geno_word & 3];
      geno_word >>= 2;
    }
    ++widx;
  }
}

void genoarr_phased_to_allele_codes(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, uint32_t sample_ct, uint32_t phasepresent_ct, unsigned char* phasebytes, int32_t* allele_codes) {
  // phasebytes can be nullptr
  genoarr_to_allele_codes(genoarr, sample_ct, allele_codes);
  uint64_t* allele_codes_alias64 = (uint64_t*)allele_codes;
  uint32_t sample_uidx = 0;
  if (!phasebytes) {
    for (uint32_t phased_idx = 0; phased_idx < phasepresent_ct; ++phased_idx, ++sample_uidx) {
      next_set_unsafe_ck(phasepresent, &sample_uidx);
      if (IS_SET(phaseinfo, sample_uidx)) {
	// 1|0
	allele_codes_alias64[sample_uidx] = 1;
      }
    }
    return;
  }
  // 0 and 2 = homozygous, automatically phased; otherwise patch in from
  // phaseinfo if phasepresent_ct is nonzero
  // so, start off by extracting low bit from each pair and flipping it
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBytesPerWord;
  const quarterword_t* read_alias = (const quarterword_t*)genoarr;
  uintptr_t* write_walias = (uintptr_t*)phasebytes;
  uint32_t widx = 0;
  while (1) {
    uintptr_t qw = read_alias[widx];
#ifdef __LP64__
    qw = (qw | (qw << 24)) & kMask000000FF;
#endif
    qw = (qw | (qw << 12)) & kMask000F;
    qw = (~(qw | (qw << 6))) & kMask0101;
    if (widx == word_ct_m1) {
      memcpy(&(write_walias[widx]), &qw, MOD_NZ(sample_ct, kBytesPerWord));
      break;
    }
    write_walias[widx++] = qw;
  }
  for (uint32_t phased_idx = 0; phased_idx < phasepresent_ct; ++phased_idx, ++sample_uidx) {
    next_set_unsafe_ck(phasepresent, &sample_uidx);
    phasebytes[sample_uidx] = 1;
    if (IS_SET(phaseinfo, sample_uidx)) {
      allele_codes_alias64[sample_uidx] = 1;
    }
  }
}

// missing = -9
static const int32_t geno_to_hap0_code[6] = {0, 0, 1, -9, 0, 1};
static const int32_t geno_to_hap1_code[6] = {0, 1, 1, -9, 0, 0};

// todo: write version of this which fills phasebytes
void genoarr_phased_to_hap_codes(const uintptr_t* genoarr, const uintptr_t* phaseinfo, uint32_t variant_batch_size, int32_t* hap0_codes_iter, int32_t* hap1_codes_iter) {
  // assumes genoarr and phaseinfo have already been transposed
  const uint32_t word_ct_m1 = (variant_batch_size - 1) / kBitsPerWordD2;
  const halfword_t* phaseinfo_alias = (const halfword_t*)phaseinfo;
  uint32_t subgroup_len = kBitsPerWordD2;
  uint32_t widx = 0;
  while (1) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
	return;
      }
      subgroup_len = MOD_NZ(variant_batch_size, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    uintptr_t phaseinfo_hw = phaseinfo_alias[widx];
    for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
      const uintptr_t cur_pgeno_code = (geno_word & 3) + 4 * (phaseinfo_hw & 1);
      *hap0_codes_iter++ = geno_to_hap0_code[cur_pgeno_code];
      *hap1_codes_iter++ = geno_to_hap1_code[cur_pgeno_code];
      geno_word >>= 2;
      phaseinfo_hw >>= 1;
    }
    ++widx;
  }
}

void bytes_to_bits_unsafe(const uint8_t* boolbytes, uint32_t sample_ct, uintptr_t* bitarr) {
  const uint32_t ull_ct_m1 = (sample_ct - 1) / 8;
  const uint64_t* read_alias = (const uint64_t*)boolbytes;
  unsigned char* write_alias = (unsigned char*)bitarr;
  uint32_t ullidx = 0;
  while (1) {
    uint64_t cur_ull;
    if (ullidx >= ull_ct_m1) {
      if (ullidx > ull_ct_m1) {
	return;
      }
      cur_ull = 0;
      memcpy(&cur_ull, &(read_alias[ullidx]), MOD_NZ(sample_ct, 8));
    } else {
      cur_ull = read_alias[ullidx];
    }
    // assuming boolbytes is 0/1-valued, this multiply-and-shift maps binary
    //  h0000000g0000000f... to binary hgfedcba.
    //  ^       ^       ^
    //  |       |       |
    // 56      48      40
    // (the constant has bits 0, 7, 14, 21, 28, 35, 42, and 49 set)
    write_alias[ullidx++] = (unsigned char)((cur_ull * 0x2040810204081LLU) >> 49);
  }
}

void bytes_to_genoarr_unsafe(const int8_t* genobytes, uint32_t sample_ct, uintptr_t* genoarr) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBytesPerWord;
  const uintptr_t* read_walias = (const uintptr_t*)genobytes;
  quarterword_t* write_alias = (quarterword_t*)genoarr;
  uint32_t widx = 0;
  while (1) {
    uintptr_t ww;
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
	return;
      }
      ww = 0;
      memcpy(&ww, &(read_walias[widx]), MOD_NZ(sample_ct, kBytesPerWord));
    } else {
      ww = read_walias[widx];
    }
    ww &= kMask0303;
    ww = (ww | (ww >> 6)) & kMask000F;
#ifdef __LP64__
    ww = (ww | (ww >> 12)) & kMask000000FF;
    write_alias[widx] = (quarterword_t)(ww | (ww >> 24));
#else
    write_alias[widx] = (quarterword_t)(ww | (ww >> 12));
#endif
    ++widx;
  }
}

void allele_codes_to_genoarr_unsafe(const int32_t* allele_codes, const unsigned char* phasepresent_bytes, uint32_t sample_ct, uintptr_t* genoarr, uintptr_t* phasepresent, uintptr_t* phaseinfo) {
  // - If phasepresent_bytes is nullptr, phasepresent is not updated.  In this
  //   case, phaseinfo is updated iff it's not nullptr.  It's okay for both
  //   phasepresent and phaseinfo to be nullptr here.
  // - Otherwise, phasepresent and phaseinfo are always updated; neither can be
  //   nullptr.
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBytesPerWord;
  uint32_t subgroup_len = kBitsPerWordD2;
  uint32_t widx = 0;
  const int32_t* read_alias = allele_codes;
  halfword_t* phaseinfo_alias = (halfword_t*)phaseinfo;
  if (!phasepresent_bytes) {
    while (1) {
      if (widx >= word_ct_m1) {
	if (widx > word_ct_m1) {
	  return;
	}
	subgroup_len = MOD_NZ(sample_ct, kBitsPerWordD2);
      }
      uintptr_t geno_write_word = 0;
      if (!phaseinfo) {
	for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
	  // 0,0 -> 0
	  // 0,1 or 1,0 -> 1
	  // 1,1 -> 2
	  // -9,-9 -> 3
	  // undefined behavior on e.g. 0,2
	  const uint32_t first_code = (uint32_t)(*read_alias++);
	  const uint32_t second_code = (uint32_t)(*read_alias++);
	  uintptr_t cur_geno;
	  if (first_code <= 1) {
	    cur_geno = first_code + second_code;
	  } else {
	    // todo: test whether branchless is better
	    // (in practice, this will usually be predictable?)
	    cur_geno = 3;
	  }
	  geno_write_word |= (cur_geno << (uii * 2));
	}
      } else {
	halfword_t phaseinfo_write_hw = 0;
	for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
	  // set phaseinfo_write_hw bit iff 1,0
	  const uint32_t first_code = (uint32_t)(*read_alias++);
	  const uint32_t second_code = (uint32_t)(*read_alias++);
	  uintptr_t cur_geno;
	  if (first_code <= 1) {
	    cur_geno = first_code + second_code;
	    phaseinfo_write_hw |= (cur_geno & first_code) << uii;
	  } else {
	    // todo: test whether branchless is better
	    // (in practice, this will usually be predictable?)
	    cur_geno = 3;
	  }
	  geno_write_word |= (cur_geno << (uii * 2));
	}
	phaseinfo_alias[widx] = phaseinfo_write_hw;
      }
      genoarr[widx] = geno_write_word;
      ++widx;
    }
  }
  const unsigned char* phasepresent_bytes_iter = phasepresent_bytes;
  halfword_t* phasepresent_alias = (halfword_t*)phasepresent;
  while (1) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
	return;
      }
      subgroup_len = MOD_NZ(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_write_word = 0;
    halfword_t phasepresent_write_hw = 0;
    halfword_t phaseinfo_write_hw = 0;
    for (uint32_t uii = 0; uii < subgroup_len; ++uii) {
      const uint32_t first_code = (uint32_t)(*read_alias++);
      const uint32_t second_code = (uint32_t)(*read_alias++);
      uintptr_t cur_geno;
      if (first_code <= 1) {
	cur_geno = first_code + second_code;
	const uint32_t cur_phasepresent = cur_geno & phasepresent_bytes_iter[uii];
	phasepresent_write_hw |= cur_phasepresent << uii;
	phaseinfo_write_hw |= (cur_phasepresent & first_code) << uii;
      } else {
	cur_geno = 3;
      }
      geno_write_word |= (cur_geno << (uii * 2));
    }
    phasepresent_bytes_iter = &(phasepresent_bytes_iter[subgroup_len]);
    phasepresent_alias[widx] = phasepresent_write_hw;
    phaseinfo_alias[widx] = phaseinfo_write_hw;
    genoarr[widx] = geno_write_word;
    ++widx;
  }
}

#ifdef __cplusplus
} // namespace plink2
#endif
