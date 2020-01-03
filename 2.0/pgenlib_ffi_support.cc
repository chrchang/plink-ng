// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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

#include "pgenlib_ffi_support.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void GenoarrToBytesMinus9(const uintptr_t* genoarr, uint32_t sample_ct, int8_t* genobytes) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBytesPerWord;
  const Quarterword* read_alias = R_CAST(const Quarterword*, genoarr);
  uintptr_t* write_walias = R_CAST(uintptr_t*, genobytes);
  for (uint32_t widx = 0; ; ++widx) {
    uintptr_t qw = Unpack0303(read_alias[widx]);
    // now each byte is in {0, 1, 2, 3}.  Convert the 3s to -9s in a branchless
    // manner.
    // (-9) - 3 = -12, which is represented as 244 in a uint8_t
    const uintptr_t geno_missing = qw & (qw >> 1) & kMask0101;
    qw += geno_missing * 244;
    if (widx == word_ct_m1) {
      SubwordStore(qw, ModNz(sample_ct, kBytesPerWord), &(write_walias[widx]));
      return;
    }
    write_walias[widx] = qw;
  }
}

static const int32_t kGenoInt32Quads[1024] ALIGNV16 = QUAD_TABLE256(0, 1, 2, -9);

void GenoarrToInt32sMinus9(const uintptr_t* genoarr, uint32_t sample_ct, int32_t* geno_int32) {
  GenoarrLookup256x4bx4(genoarr, kGenoInt32Quads, sample_ct, geno_int32);
}

// todo: use GenoarrLookup16x8bx2()
static const int64_t kGenoToInt64[4] = {0, 1, 2, -9};

void GenoarrToInt64sMinus9(const uintptr_t* genoarr, uint32_t sample_ct, int64_t* geno_int64) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  int64_t* write_iter = geno_int64;
  uint32_t subgroup_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        return;
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      *write_iter++ = kGenoToInt64[geno_word & 3];
      geno_word >>= 2;
    }
  }
}

// missing = -9
const double kGenoDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, -9.0);

const uint64_t kGenoToIntcodeDPairs[32] ALIGNV16 = PAIR_TABLE16(0, 0x100000000LLU, 0x100000001LLU, 0xfffffff7fffffff7LLU);

void GenoarrPhasedToAlleleCodes(const uint64_t* genoarr_to_intcode_pair_table, const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, uint32_t sample_ct, uint32_t phasepresent_ct, unsigned char* phasebytes, int32_t* allele_codes) {
  // phasebytes can be nullptr, phasepresent cannot
  GenoarrToAlleleCodes(genoarr_to_intcode_pair_table, genoarr, sample_ct, allele_codes);
  uint64_t* allele_codes_alias64 = R_CAST(uint64_t*, allele_codes);
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = phasepresent[0];
  if (!phasebytes) {
    for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
      const uintptr_t sample_uidx = BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
      if (IsSet(phaseinfo, sample_uidx)) {
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
  const Quarterword* read_alias = R_CAST(const Quarterword*, genoarr);
  uintptr_t* write_walias = R_CAST(uintptr_t*, phasebytes);
  for (uint32_t widx = 0; ; ++widx) {
    uintptr_t qw = Unpack0303(read_alias[widx]);
    qw = (~qw) & kMask0101;
    if (widx == word_ct_m1) {
      SubwordStore(qw, ModNz(sample_ct, kBytesPerWord), &(write_walias[widx]));
      break;
    }
    write_walias[widx] = qw;
  }
  for (uint32_t phased_idx = 0; phased_idx != phasepresent_ct; ++phased_idx) {
    const uintptr_t sample_uidx = BitIter1(phasepresent, &sample_uidx_base, &cur_bits);
    phasebytes[sample_uidx] = 1;
    if (IsSet(phaseinfo, sample_uidx)) {
      allele_codes_alias64[sample_uidx] = 1;
    }
  }
}

// missing = -9
// may want a double-lookup function for this
static const int32_t kGenoToHap0Code[6] = {0, 0, 1, -9, 0, 1};
static const int32_t kGenoToHap1Code[6] = {0, 1, 1, -9, 0, 0};

// todo: write version of this which fills phasebytes
void GenoarrPhasedToHapCodes(const uintptr_t* genoarr, const uintptr_t* phaseinfo, uint32_t variant_batch_size, int32_t* hap0_codes_iter, int32_t* hap1_codes_iter) {
  // assumes genoarr and phaseinfo have already been transposed
  const uint32_t word_ct_m1 = (variant_batch_size - 1) / kBitsPerWordD2;
  const Halfword* phaseinfo_alias = R_CAST(const Halfword*, phaseinfo);
  uint32_t subgroup_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        return;
      }
      subgroup_len = ModNz(variant_batch_size, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    uintptr_t phaseinfo_hw = phaseinfo_alias[widx];
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      const uintptr_t cur_pgeno_code = (geno_word & 3) + 4 * (phaseinfo_hw & 1);
      *hap0_codes_iter++ = kGenoToHap0Code[cur_pgeno_code];
      *hap1_codes_iter++ = kGenoToHap1Code[cur_pgeno_code];
      geno_word >>= 2;
      phaseinfo_hw >>= 1;
    }
  }
}

// todo: use GenoarrLookup256x4bx4()
static const float kGenoToFloat[4] = {0.0, 1.0, 2.0, -9.0};

void Dosage16ToFloatsMinus9(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, float* geno_float) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  float* write_iter = geno_float;
  uint32_t subgroup_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        break;
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      *write_iter++ = kGenoToFloat[geno_word & 3];
      geno_word >>= 2;
    }
  }
  if (dosage_ct) {
    const uint16_t* dosage_main_iter = dosage_main;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = dosage_present[0];
    for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
      const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
      // multiply by 2^{-14}
      geno_float[sample_uidx] = S_CAST(float, *dosage_main_iter++) * S_CAST(float, 0.00006103515625);
    }
  }
}

void Dosage16ToDoubles(const double* geno_double_pair_table, const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double) {
  GenoarrLookup16x8bx2(genoarr, geno_double_pair_table, sample_ct, geno_double);
  if (dosage_ct) {
    const uint16_t* dosage_main_iter = dosage_main;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = dosage_present[0];
    for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
      const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
      geno_double[sample_uidx] = S_CAST(double, *dosage_main_iter++) * 0.00006103515625;
    }
  }
}

BoolErr Dosage16ToDoublesMeanimpute(const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, double* geno_double) {
  STD_ARRAY_DECL(uint32_t, 4, genocounts);
  double geno_double_pair_buf[32];
  if (!dosage_ct) {
    GenoarrCountFreqsUnsafe(genoarr, sample_ct, genocounts);
    const double* geno_double_pair_table = kGenoDoublePairs;
    if (genocounts[3]) {
      const uint32_t denom = sample_ct - genocounts[3];
      if (!denom) {
        return 1;
      }
      const uint32_t numer = genocounts[1] + 2 * genocounts[2];
      const double missing_val = u63tod(numer) / u31tod(denom);
      geno_double_pair_buf[0] = 0.0;
      geno_double_pair_buf[2] = 1.0;
      geno_double_pair_buf[4] = 2.0;
      geno_double_pair_buf[6] = missing_val;
      InitLookup16x8bx2(geno_double_pair_buf);
      geno_double_pair_table = geno_double_pair_buf;
    }
    GenoarrLookup16x8bx2(genoarr, geno_double_pair_table, sample_ct, geno_double);
    return 0;
  }
  // In the generic case, it may be faster to check for the existence of a
  // missing value before calling GenoarrCountInvsubsetFreqs2() (since if there
  // are no missing values, we don't need to count at all).  However, we assume
  // the caller is using this function over Dosage16ToDoubles() for a reason.
  GenoarrCountInvsubsetFreqs2(genoarr, dosage_present, sample_ct, sample_ct - dosage_ct, genocounts);
  const double* geno_double_pair_table = kGenoDoublePairs;
  if (genocounts[3]) {
    uint64_t denom = sample_ct - genocounts[3];
    if (!denom) {
      return 1;
    }
    denom *= 16384LLU;
    uint64_t numer = 0;
    for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
      numer += dosage_main[dosage_idx];
    }
    numer += 16384LLU * (genocounts[1] + 2 * genocounts[2]);
    const double missing_val = u63tod(numer) / u63tod(denom);
    geno_double_pair_buf[0] = 0.0;
    geno_double_pair_buf[2] = 1.0;
    geno_double_pair_buf[4] = 2.0;
    geno_double_pair_buf[6] = missing_val;
    InitLookup16x8bx2(geno_double_pair_buf);
    geno_double_pair_table = geno_double_pair_buf;
  }
  GenoarrLookup16x8bx2(genoarr, geno_double_pair_table, sample_ct, geno_double);
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
    geno_double[sample_uidx] = S_CAST(double, dosage_main[dosage_idx]) * 0.00006103515625;
  }
  return 0;
}

double LinearCombinationMeanimpute(const double* weights, const uintptr_t* genoarr, const uintptr_t* dosage_present, const uint16_t* dosage_main, uint32_t sample_ct, uint32_t dosage_ct) {
  const uint32_t word_ct = DivUp(sample_ct, kBitsPerWordD2);
  double result = 0.0;
  double result2 = 0.0;
  double miss_weight = 0.0;
  if (!dosage_ct) {
    for (uint32_t widx = 0; widx != word_ct; ++widx) {
      const uintptr_t geno_word = genoarr[widx];
      if (!geno_word) {
        continue;
      }
      const double* cur_weights = &(weights[widx * kBitsPerWordD2]);
      uintptr_t geno_word1 = geno_word & kMask5555;
      uintptr_t geno_word2 = (geno_word >> 1) & kMask5555;
      uintptr_t geno_missing_word = geno_word1 & geno_word2;
      geno_word1 ^= geno_missing_word;
      while (geno_word1) {
        const uint32_t sample_idx_lowbits = ctzw(geno_word1) / 2;
        result += cur_weights[sample_idx_lowbits];
        geno_word1 &= geno_word1 - 1;
      }
      geno_word2 ^= geno_missing_word;
      while (geno_word2) {
        const uint32_t sample_idx_lowbits = ctzw(geno_word2) / 2;
        result2 += cur_weights[sample_idx_lowbits];
        geno_word2 &= geno_word2 - 1;
      }
      while (geno_missing_word) {
        const uint32_t sample_idx_lowbits = ctzw(geno_missing_word) / 2;
        miss_weight += cur_weights[sample_idx_lowbits];
        geno_missing_word &= geno_missing_word - 1;
      }
    }
    result += 2 * result2;
    if (miss_weight != 0.0) {
      // bugfix (29 Oct 2019): previous mean-imputation formula was based on
      // *weighted* MAF, which was obviously nonsense when negative weights
      // were present.
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      GenoarrCountFreqsUnsafe(genoarr, sample_ct, genocounts);
      const double numer = u63tod(genocounts[1] + 2 * genocounts[2]);
      const double denom = u31tod(sample_ct - genocounts[3]);
      result += miss_weight * (numer / denom);
    }
    return result;
  }
  const Halfword* dosage_present_hws = R_CAST(const Halfword*, dosage_present);
  uint32_t onealt_ct = 0;
  uint32_t twoalt_ct = 0;
  uint32_t missing_ct = 0;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t geno_word = genoarr[widx];
    if (geno_word) {
      const double* cur_weights = &(weights[widx * kBitsPerWordD2]);
      uintptr_t geno_word1 = geno_word & kMask5555;
      uintptr_t geno_word2 = (geno_word >> 1) & kMask5555;
      uintptr_t geno_missing_word = geno_word1 & geno_word2;
      const uintptr_t mask_word = ~(geno_missing_word | UnpackHalfwordToWord(dosage_present_hws[widx]));
      geno_word1 &= mask_word;
      while (geno_word1) {
        const uint32_t sample_idx_lowbits = ctzw(geno_word1) / 2;
        result += cur_weights[sample_idx_lowbits];
        // probably sparse enough that this is faster than popcount?
        ++onealt_ct;
        geno_word1 &= geno_word1 - 1;
      }
      geno_word2 &= mask_word;
      while (geno_word2) {
        const uint32_t sample_idx_lowbits = ctzw(geno_word2) / 2;
        result2 += cur_weights[sample_idx_lowbits];
        ++twoalt_ct;
        geno_word2 &= geno_word2 - 1;
      }
      while (geno_missing_word) {
        const uint32_t sample_idx_lowbits = ctzw(geno_missing_word) / 2;
        miss_weight += cur_weights[sample_idx_lowbits];
        ++missing_ct;
        geno_missing_word &= geno_missing_word - 1;
      }
    }
  }
  result += result2 * 2;
  const uint16_t* dosage_main_iter = dosage_main;
  double resultx = 0.0;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  if (miss_weight == 0.0) {
    for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
      const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
      resultx += S_CAST(double, *dosage_main_iter++) * weights[sample_uidx];
    }
    result += 0.00006103515625 * resultx;
    return result;
  }
  // Also need to track dosage-sum for mean-imputation.
  uint64_t dosage_sum = 0;
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
    const uint32_t cur_dosage = *dosage_main_iter++;
    dosage_sum += cur_dosage;
    resultx += u31tod(cur_dosage) * weights[sample_uidx];
  }
  result += 0.00006103515625 * resultx;
  const double numer = u63tod(16384 * S_CAST(uint64_t, onealt_ct + 2 * twoalt_ct) + dosage_sum);
  const double denom = 16384 * u31tod(sample_ct - missing_ct);
  result += miss_weight * (numer / denom);
  return result;
}

void BytesToBitsUnsafe(const uint8_t* boolbytes, uint32_t sample_ct, uintptr_t* bitarr) {
  const uint32_t ull_ct_m1 = (sample_ct - 1) / 8;
  const uint64_t* read_alias = R_CAST(const uint64_t*, boolbytes);
  unsigned char* write_alias = R_CAST(unsigned char*, bitarr);
  for (uint32_t ullidx = 0; ; ++ullidx) {
    uint64_t cur_ull;
    if (ullidx >= ull_ct_m1) {
      if (ullidx > ull_ct_m1) {
        return;
      }
      cur_ull = SubwordLoad(&(read_alias[ullidx]), ModNz(sample_ct, 8));
    } else {
      cur_ull = read_alias[ullidx];
    }
    // assuming boolbytes is 0/1-valued, this multiply-and-shift maps binary
    //  h0000000g0000000f... to binary hgfedcba.
    //  ^       ^       ^
    //  |       |       |
    // 56      48      40
    // (the constant has bits 0, 7, 14, 21, 28, 35, 42, and 49 set)
    // (can also use _pext_u64() in AVX2 case)
    write_alias[ullidx] = S_CAST(unsigned char, (cur_ull * 0x2040810204081LLU) >> 49);
  }
}

void BytesToGenoarrUnsafe(const int8_t* genobytes, uint32_t sample_ct, uintptr_t* genoarr) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBytesPerWord;
  const uintptr_t* read_walias = R_CAST(const uintptr_t*, genobytes);
  Quarterword* write_alias = R_CAST(Quarterword*, genoarr);
  for (uint32_t widx = 0; ; ++widx) {
    uintptr_t ww;
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        return;
      }
      ww = SubwordLoad(&(read_walias[widx]), ModNz(sample_ct, kBytesPerWord));
    } else {
      ww = read_walias[widx];
    }
    write_alias[widx] = Pack0303Mask(ww);
  }
}

void AlleleCodesToGenoarrUnsafe(const int32_t* allele_codes, const unsigned char* phasepresent_bytes, uint32_t sample_ct, uintptr_t* genoarr, uintptr_t* phasepresent, uintptr_t* phaseinfo) {
  // - If phasepresent_bytes is nullptr, phasepresent is not updated.  In this
  //   case, phaseinfo is updated iff it's not nullptr.  It's okay for both
  //   phasepresent and phaseinfo to be nullptr here.
  // - Otherwise, phasepresent and phaseinfo are always updated; neither can be
  //   nullptr.
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t subgroup_len = kBitsPerWordD2;
  const uint32_t* read_alias = R_CAST(const uint32_t*, allele_codes);
  Halfword* phaseinfo_alias = R_CAST(Halfword*, phaseinfo);
  if (!phasepresent_bytes) {
    for (uint32_t widx = 0; ; ++widx) {
      if (widx >= word_ct_m1) {
        if (widx > word_ct_m1) {
          return;
        }
        subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
      }
      uintptr_t geno_write_word = 0;
      if (!phaseinfo) {
        for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
          // 0,0 -> 0
          // 0,1 or 1,0 -> 1
          // 1,1 -> 2
          // -9,-9 -> 3
          // undefined behavior on e.g. 0,2
          const uint32_t first_code = *read_alias++;
          const uint32_t second_code = *read_alias++;
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
        Halfword phaseinfo_write_hw = 0;
        for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
          // set phaseinfo_write_hw bit iff 1,0
          const uint32_t first_code = *read_alias++;
          const uint32_t second_code = *read_alias++;
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
    }
  }
  const unsigned char* phasepresent_bytes_iter = phasepresent_bytes;
  Halfword* phasepresent_alias = R_CAST(Halfword*, phasepresent);
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        return;
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_write_word = 0;
    Halfword phasepresent_write_hw = 0;
    Halfword phaseinfo_write_hw = 0;
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      const uint32_t first_code = *read_alias++;
      const uint32_t second_code = *read_alias++;
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
  }
}

static inline uint32_t BiallelicDosage16Halfdist(uint32_t dosage_int) {
  const uint32_t dosage_int_rem = dosage_int & 16383;
  return abs_i32(S_CAST(int32_t, dosage_int_rem) - 8192);
}

void FloatsToDosage16(const float* floatarr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const float* read_iter = floatarr;
  Halfword* dosage_present_alias = R_CAST(Halfword*, dosage_present);
  uint16_t* dosage_main_iter = dosage_main;
  uint32_t subgroup_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        if (widx % 2) {
          dosage_present_alias[widx] = 0;
        }
        break;
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = 0;
    uint32_t dosage_present_hw = 0;
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != subgroup_len; ++sample_idx_lowbits) {
      // 0..2 -> 0..32768
      const float fxx = (*read_iter++) * 16384 + 0.5;
      uintptr_t cur_geno = 3;
      if ((fxx >= 0.0) && (fxx < 32769)) {
        uint32_t dosage_int = S_CAST(int32_t, fxx);
        const uint32_t cur_halfdist = BiallelicDosage16Halfdist(dosage_int);
        if (cur_halfdist >= hard_call_halfdist) {
          cur_geno = (dosage_int + (8192 * k1LU)) / 16384;
        }
        if (cur_halfdist != 8192) {
          dosage_present_hw |= 1U << sample_idx_lowbits;
          *dosage_main_iter++ = dosage_int;
        }
      }
      geno_word |= cur_geno << (2 * sample_idx_lowbits);
    }
    genoarr[widx] = geno_word;
    dosage_present_alias[widx] = dosage_present_hw;
  }
  *dosage_ct_ptr = dosage_main_iter - dosage_main;
}

void DoublesToDosage16(const double* doublearr, uint32_t sample_ct, uint32_t hard_call_halfdist, uintptr_t* genoarr, uintptr_t* dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const double* read_iter = doublearr;
  Halfword* dosage_present_alias = R_CAST(Halfword*, dosage_present);
  uint16_t* dosage_main_iter = dosage_main;
  uint32_t subgroup_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= word_ct_m1) {
      if (widx > word_ct_m1) {
        if (widx % 2) {
          dosage_present_alias[widx] = 0;
        }
        break;
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = 0;
    uint32_t dosage_present_hw = 0;
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != subgroup_len; ++sample_idx_lowbits) {
      // 0..2 -> 0..32768
      const double dxx = (*read_iter++) * 16384 + 0.5;
      uintptr_t cur_geno = 3;
      if ((dxx >= 0.0) && (dxx < 32769)) {
        uint32_t dosage_int = S_CAST(int32_t, dxx);
        const uint32_t cur_halfdist = BiallelicDosage16Halfdist(dosage_int);
        if (cur_halfdist >= hard_call_halfdist) {
          cur_geno = (dosage_int + (8192 * k1LU)) / 16384;
        }
        if (cur_halfdist != 8192) {
          dosage_present_hw |= 1U << sample_idx_lowbits;
          *dosage_main_iter++ = dosage_int;
        }
      }
      geno_word |= cur_geno << (2 * sample_idx_lowbits);
    }
    genoarr[widx] = geno_word;
    dosage_present_alias[widx] = dosage_present_hw;
  }
  *dosage_ct_ptr = dosage_main_iter - dosage_main;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
