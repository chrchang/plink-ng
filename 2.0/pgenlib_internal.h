#ifndef __PGENLIB_INTERNAL_H__
#define __PGENLIB_INTERNAL_H__

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


// Low-level C99/C++03/C++11 library for reading .pgen (PLINK 2.0 binary) files
// (designed to produce good lowest-common-denominator binaries across
// Windows/OS X/Linux).
//
// File format design:
// - With the header loaded, it is possible to efficiently access a variant by
//   its index.  Since records can now be variable-length, this sometimes
//   requires storage of record lengths.
// - Due to the power of LD-based compression, we permit a variant record to
//   just store a list of differences from an earlier, fully stored variant.
//   However, only short-range dependence is permitted; sequential processing
//   of the file only requires caching of the most recent explicitly stored
//   variant.
// - Like the plink1 format, this is balanced for relatively easy reading and
//   writing; in particular, the mode-0x10/0x11 header is not read-optimized,
//   it passes up some obvious compression opportunities which would make it
//   more difficult to write e.g. an efficient file merger.  This isn't a big
//   deal if we don't have a huge number of one-sample .pgen files sharing a
//   single .bim file (or equivalent).  (If they don't share the same .bim
//   file, .bim overhead > .pgen overhead.)  If we ever do, we can define an
//   additional mode to handle that case more efficiently.
// - Building blocks are arrays of 1-bit, 2-bit, 4-bit, 1-byte, 2-byte, 3-byte,
//   and 4-byte values.  3/5/6/7(/9...)-bit values don't play well with
//   bitwise operations, and when it's important, there's usually a natural way
//   to split them into power-of-2-bit components.
//   (unsigned integers known to be smaller than 2^24, but not known to be
//   smaller than 2^16, are stored as 3-byte values on disk and "decompressed"
//   to uint32_t during loading.)
// - Missing value is usually all-1s.  (Only exceptions right now: plink1
//   backward compatibility mode; presence/absence of rare alts for variants
//   with >2 alt alleles is an array of 1-bit values, where absence = 0; and
//   presence/absence of phasing info is similar.)  Time to move away from 01
//   nonsense.
// - Individual variant records are prohibited from being >= 4GiB, to reduce
//   integer overflow issues.  (This may be reduced to 2GiB later, but I'll
//   attempt to handle the 2-4GiB range properly for now since it's conceivable
//   for multiallelic records in very large datasets to reach that size.)
// - (later todo: include stuff like file creation command in .pvar header;
//   that doesn't really belong in a binary file.)

// Additional parameter conventions:
// - "quaterarr" indicates a word-aligned, packed array of 2-bit values, while
//   "quatervec" is the vector-aligned equivalent.  "nibblearr" marks the much
//   rarer case of a packed array of 4-bit values.
// - "quatervec_01" indicates a packed, vector-aligned array of 2-bit values
//   where each value is zero or one.  This data structure was used quite a bit
//   bit by plink 1.9 for operating on a subset of a 2-bit-genotype array.
// - "genovec" indicates a quatervec containing genotype information.
// - "interleaved_vec" is the plink 2.0 replacement for quatervec_01: we
//   basically stack pairs of adjacent vectors on top of each other and unpack
//   on the fly, since that tends to be faster than having to access twice as
//   much memory.

#include "plink2_base.h"

// 10000 * major + 100 * minor + patch
// Exception to CONSTI32, since we want the preprocessor to have access to this
// value.  Named with all caps as a consequence.
#define PGENLIB_INTERNAL_VERNUM 1200

#ifdef __cplusplus
namespace plink2 {
#endif

// other configuration-ish values needed by plink2_common subset
typedef unsigned char AlleleCode;
typedef uint16_t DoubleAlleleCode;
static_assert(sizeof(DoubleAlleleCode) == 2 * sizeof(AlleleCode), "Inconsistent AlleleCode and DoubleAlleleCode definitions.");
// Set this to 65534 if AlleleCode is uint16_t, 2^24 - 1 if uint32_t.
CONSTI32(kPglMaxAltAlleleCt, 254);
#ifdef __cplusplus
#  define kMissingAlleleCode S_CAST(plink2::AlleleCode, -1)
#else
#  define kMissingAlleleCode S_CAST(AlleleCode, -1)
#endif
CONSTI32(kAlleleCodesPerVec, kBytesPerVec / sizeof(AlleleCode));

// more verbose than (val + 3) / 4, but may as well make semantic meaning
// obvious; any explicit DivUp(val, 4) expressions should have a different
// meaning
// (not needed for bitct -> bytect, DivUp(val, CHAR_BIT) is clear enough)
HEADER_INLINE uintptr_t QuaterCtToByteCt(uintptr_t val) {
  return DivUp(val, 4);
}

HEADER_INLINE uintptr_t QuaterCtToVecCt(uintptr_t val) {
  return DivUp(val, kQuatersPerVec);
}

HEADER_INLINE uintptr_t QuaterCtToWordCt(uintptr_t val) {
  return DivUp(val, kBitsPerWordD2);
}

HEADER_INLINE uintptr_t QuaterCtToAlignedWordCt(uintptr_t val) {
  return kWordsPerVec * QuaterCtToVecCt(val);
}

HEADER_INLINE uintptr_t QuaterCtToCachelineCt(uintptr_t val) {
  return DivUp(val, kQuatersPerCacheline);
}

HEADER_INLINE uintptr_t AlleleCodeCtToVecCt(uintptr_t val) {
  return DivUp(val, kAlleleCodesPerVec);
}

HEADER_INLINE uintptr_t AlleleCodeCtToAlignedWordCt(uintptr_t val) {
  return kWordsPerVec * AlleleCodeCtToVecCt(val);
}

HEADER_INLINE uintptr_t GetQuaterarrEntry(const uintptr_t* quaterarr, uint32_t idx) {
  return (quaterarr[idx / kBitsPerWordD2] >> (2 * (idx % kBitsPerWordD2))) & 3;
}

// todo: check if this optimizes newval=0 out
HEADER_INLINE void AssignQuaterarrEntry(uint32_t idx, uintptr_t newval, uintptr_t* quaterarr) {
  const uint32_t bit_shift_ct = 2 * (idx % kBitsPerWordD2);
  uintptr_t* wordp = &(quaterarr[idx / kBitsPerWordD2]);
  *wordp = ((*wordp) & (~((3 * k1LU) << bit_shift_ct))) | (newval << bit_shift_ct);
}

HEADER_INLINE void ClearQuaterarrEntry(uint32_t idx, uintptr_t* quaterarr) {
  quaterarr[idx / kBitsPerWordD2] &= ~((3 * k1LU) << (idx % kBitsPerWordD2));
}

// returns a word with low bit in each pair set at each 00.
HEADER_INLINE uintptr_t Word00(uintptr_t ww) {
  return (~(ww | (ww >> 1))) & kMask5555;
}

HEADER_INLINE uintptr_t Word01(uintptr_t ww) {
  return ww & (~(ww >> 1)) & kMask5555;
}

// returns a word with *low* bit in each pair set at each 10.
HEADER_INLINE uintptr_t Word10(uintptr_t ww) {
  return (~ww) & (ww >> 1) & kMask5555;
}

HEADER_INLINE uintptr_t Word11(uintptr_t ww) {
  return ww & (ww >> 1) & kMask5555;
}

HEADER_INLINE Halfword Pack01ToHalfword(uintptr_t ww) {
  return PackWordToHalfwordMask5555(ww & (~(ww >> 1)));
}

HEADER_INLINE Halfword Pack11ToHalfword(uintptr_t ww) {
  return PackWordToHalfwordMask5555(ww & (ww >> 1));
}

#ifdef USE_SSE42
HEADER_INLINE uint32_t Popcount01Word(uintptr_t val) {
  return PopcountWord(val);
}

HEADER_INLINE uint32_t Popcount0001Word(uintptr_t val) {
  return PopcountWord(val);
}
#else
HEADER_INLINE uint32_t Popcount01Word(uintptr_t val) {
  return QuatersumWord(val);
}

HEADER_INLINE uint32_t Popcount0001Word(uintptr_t val) {
#  ifdef __LP64__
  // (val * kMask1111) >> 60 can barely overflow, sigh
  const uintptr_t val0 = val & 1;
  return (((val - val0) * kMask1111) >> 60) + val0;
#  else
  return (val * kMask1111) >> 28;
#  endif
}
#endif

// safe errstr_buf size for pgen_init_phase{1,2}()
CONSTI32(kPglErrstrBufBlen, kPglFnamesize + 256);

// assumes subset_mask has trailing zeroes up to the next vector boundary
void FillInterleavedMaskVec(const uintptr_t* __restrict subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec);

HEADER_INLINE void CopyQuaterarr(const uintptr_t* __restrict source_quaterarr, uint32_t quaterarr_entry_ct, uintptr_t* __restrict target_quaterarr) {
  memcpy(target_quaterarr, source_quaterarr, QuaterCtToWordCt(quaterarr_entry_ct) * kBytesPerWord);
}

// may want bit past the end of subset_mask (i.e. position
// raw_quaterarr_entry_ct) to always be allocated and unset.  This removes the
// need for some explicit end-of-bitarray checks.
void CopyQuaterarrNonemptySubset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_quaterarr);

// Copies a bit from raw_bitarr for each genovec entry matching match_word.
// (match_word must be a multiple of kMask5555.)
void CopyGenomatchSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict genovec, uintptr_t match_word, uint32_t write_bit_idx_start, uint32_t bit_ct, uintptr_t* __restrict output_bitarr);

void ExpandBytearrFromGenovec(const void* __restrict compact_bitarr, const uintptr_t* __restrict genovec, uintptr_t match_word, uint32_t genoword_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target);


// These two functions are "unsafe" since they assume trailing bits of last
// genovec word are zeroed out.
void GenovecCount12Unsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* __restrict raw_01_ctp, uint32_t* __restrict raw_10_ctp);

void GenovecCountFreqsUnsafe(const uintptr_t* genovec, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts);


void GenovecCountSubsetFreqs(const uintptr_t* __restrict genovec, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts);

// slower GenovecCountSubsetFreqs() which does not require
// sample_include_interleaved_vec to be precomputed
void GenovecCountSubsetFreqs2(const uintptr_t* __restrict genovec, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts);

void GenoarrCountSubsetIntersectFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict subset1, const uintptr_t* __restrict subset2, uint32_t raw_sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts);

void GenovecInvertUnsafe(uint32_t sample_ct, uintptr_t* genovec);

HEADER_INLINE uintptr_t InvertGenoWordUnsafe(uintptr_t geno_word) {
  return (geno_word ^ ((~(geno_word << 1)) & kMaskAAAA));
}

// too easy to forget to multiply by 2
HEADER_INLINE void ZeroTrailingQuaters(uintptr_t quater_ct, uintptr_t* bitarr) {
  ZeroTrailingBits(quater_ct * 2, bitarr);
}

HEADER_INLINE void SetTrailingQuaters(uintptr_t quater_ct, uintptr_t* bitarr) {
  const uintptr_t trail_ct = quater_ct % kBitsPerWordD2;
  if (trail_ct) {
    bitarr[quater_ct / kBitsPerWordD2] |= (~k0LU) << (quater_ct * 2);
  }
}

// A VINT is a sequence of bytes where each byte stores just 7 bits of an
// an integer, and the high bit is set when the integer has more nonzero bits.
// See e.g.
//   https://developers.google.com/protocol-buffers/docs/encoding#varints
// (Note that protocol buffers used "group varints" at one point, but then
// abandoned them.  I suspect they'd be simultaneously slower and less
// compact here.)

HEADER_INLINE unsigned char* Vint32Append(uint32_t uii, unsigned char* buf) {
  while (uii > 127) {
    *buf++ = (uii & 127) + 128;
    uii >>= 7;
  }
  *buf++ = uii;
  return buf;
}

// Returns 0x80000000U on read-past-end instead of UINT32_MAX so overflow check
// works properly in 32-bit build.  Named "GetVint31" to make it more obvious
// that a 2^31 return value can't be legitimate.
HEADER_INLINE uint32_t GetVint31(const unsigned char* buf_end, const unsigned char** buf_iterp) {
  if (likely(buf_end > (*buf_iterp))) {
    uint32_t vint32 = *((*buf_iterp)++);
    if (vint32 <= 127) {
      return vint32;
    }
    vint32 &= 127;
    uint32_t shift = 7;
    while (likely(buf_end > (*buf_iterp))) {
      uint32_t uii = *((*buf_iterp)++);
      vint32 |= (uii & 127) << shift;
      if (uii <= 127) {
        return vint32;
      }
      shift += 7;
      // currently don't check for shift >= 32 (that's what ValidateVint31()
      // is for).
    }
  }
  return 0x80000000U;
}

// Input must be validated, or bufp must be >= 5 characters before the end of
// the read buffer.  Currently unused.
// todo: check if this has enough of a speed advantage over GetVint31() to
// justify using this in the main loops and catching SIGSEGV.  (seems to be no
// more than 3%?)
/*
HEADER_INLINE uint32_t GetVint31Unsafe(const unsigned char** buf_iterp) {
  uint32_t vint32 = *(*buf_iterp)++;
  if (vint32 <= 127) {
    return vint32;
  }
  vint32 &= 127;
  for (uint32_t shift = 7; shift != 35; shift += 7) {
    uint32_t uii = *(*buf_iterp)++;
    vint32 |= (uii & 127) << shift;
    if (uii <= 127) {
      return vint32;
    }
  }
  return 0x80000000U;
}
*/

// Does not update buf_iter.
HEADER_INLINE uint32_t PeekVint31(const unsigned char* buf_iter, const unsigned char* buf_end) {
  if (likely(buf_end > buf_iter)) {
    uint32_t vint32 = *buf_iter++;
    if (vint32 <= 127) {
      return vint32;
    }
    vint32 &= 127;
    uint32_t shift = 7;
    while (likely(buf_end > buf_iter)) {
      uint32_t uii = *buf_iter++;
      vint32 |= (uii & 127) << shift;
      if (uii <= 127) {
        return vint32;
      }
      shift += 7;
    }
  }
  return 0x80000000U;
}

/*
HEADER_INLINE uint32_t FGetVint31(FILE* ff) {
  // Can't be used when multiple threads are reading from ff.
  uint32_t vint32 = getc_unlocked(ff);
  if (vint32 <= 127) {
    return vint32;
  }
  vint32 &= 127;
  for (uint32_t shift = 7; shift != 35; shift += 7) {
    uint32_t uii = getc_unlocked(ff);
    vint32 |= (uii & 127) << shift;
    if (uii <= 127) {
      return vint32;
    }
  }
  return 0x80000000U;
}

HEADER_INLINE void FPutVint31(uint32_t uii, FILE* ff) {
  // caller's responsibility to periodically check ferror
  while (uii > 127) {
    putc_unlocked((uii & 127) + 128, ff);
    uii >>= 7;
  }
  putc_unlocked(uii, ff);
}
*/

// Need this for sparse multiallelic dosage.
HEADER_INLINE unsigned char* Vint64Append(uint64_t ullii, unsigned char* buf) {
  while (ullii > 127) {
    *buf++ = (ullii & 127) + 128;
    ullii >>= 7;
  }
  *buf++ = ullii;
  return buf;
}

// Returns 2^63 on read-past-end, and named GetVint63 to make it more obvious
// that a 2^63 return value can't be legitimate.
HEADER_INLINE uint64_t GetVint63(const unsigned char* buf_end, const unsigned char** buf_iterp) {
  if (likely(buf_end > (*buf_iterp))) {
    uint64_t vint64 = *((*buf_iterp)++);
    if (vint64 <= 127) {
      return vint64;
    }
    vint64 &= 127;
    uint32_t shift = 7;
    while (likely(buf_end > (*buf_iterp))) {
      uint64_t ullii = *((*buf_iterp)++);
      vint64 |= (ullii & 127) << shift;
      if (ullii <= 127) {
        return vint64;
      }
      shift += 7;
      // currently don't check for shift >= 64 (that's what ValidateVint63()
      // will be for).
    }
  }
  return (1LLU << 63);
}

// main batch size
CONSTI32(kPglQuaterTransposeBatch, kQuatersPerCacheline);

// word width of each matrix row
CONSTI32(kPglQuaterTransposeWords, kWordsPerCacheline);

#ifdef USE_SSE42
CONSTI32(kPglQuaterTransposeBufbytes, (kPglQuaterTransposeBatch * kPglQuaterTransposeBatch) / 4);

void TransposeQuaterblockShuffle(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, void* __restrict buf0);
#else  // !USE_SSE42
CONSTI32(kPglQuaterTransposeBufbytes, (kPglQuaterTransposeBatch * kPglQuaterTransposeBatch) / 2);

void TransposeQuaterblockNoshuffle(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1);
#endif
CONSTI32(kPglQuaterTransposeBufwords, kPglQuaterTransposeBufbytes / kBytesPerWord);

// up to 256x256; vecaligned_buf must have size 16k (shuffle) or 32k
// (noshuffle)
// important: write_iter must be allocated up to at least
//   RoundUpPow2(write_batch_size, 4) rows

HEADER_INLINE void TransposeQuaterblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecW* vecaligned_buf) {
#ifdef USE_SSE42
  // assert(!(write_ul_stride % 2));
  TransposeQuaterblockShuffle(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, vecaligned_buf);
#else
  TransposeQuaterblockNoshuffle(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, R_CAST(unsigned char*, vecaligned_buf), &(R_CAST(unsigned char*, vecaligned_buf)[kPglQuaterTransposeBufbytes / 2]));
#endif
}


// replaces each x with (32768 - x)
// okay for dosage_main to be nullptr if dosage_ct == 0
void BiallelicDosage16Invert(uint32_t dosage_ct, uint16_t* dosage_main);

// replaces each x with -x
void BiallelicDphase16Invert(uint32_t dphase_ct, int16_t* dphase_delta);

// currently does zero trailing halfword
void GenovecToMissingnessUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict missingness);

// currently does not zero trailing halfword
void GenovecToNonmissingnessUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict nonmissingness);

// hom_buf gets set bits when genoarr value is 0 or 2.
// ref2het_buf gets set bits when genoarr value is 0 or 1.
// N.B. assumes trailing bits of loadbuf have been filled with 1s, not 0s
// Also takes genoarr word count instead of sample count.
void SplitHomRef2hetUnsafeW(const uintptr_t* genoarr, uint32_t inword_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf);


void SplitHomRef2het(const uintptr_t* genoarr, uint32_t sample_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf);


// These functions use 16- or 256-element lookup tables to apply functions of
// the form
//   f: {0,1,2,3} -> x
// to genoarr, saving the output to result[].
// 256-element tables result in a substantially faster inner loop, but they are
// more expensive to set up and consume a non-negligible fractino of L1 cache,
// they aren't always the right choice.
// When lookup table rows are 16 bytes, they are assumed to be 16-byte aligned
// in 64-bit builds.  result[] is not assumed to be aligned.
void GenoarrLookup16x4bx2(const uintptr_t* genoarr, const void* table16x4bx2, uint32_t sample_ct, void* result);

void GenoarrLookup256x2bx4(const uintptr_t* genoarr, const void* table256x2bx4, uint32_t sample_ct, void* result);

void GenoarrLookup16x8bx2(const uintptr_t* genoarr, const void* table16x8bx2, uint32_t sample_ct, void* result);

void GenoarrLookup256x4bx4(const uintptr_t* genoarr, const void* table256x4bx4, uint32_t sample_ct, void* result);

// Lookup table initialization functions.  table[0][0], [1][0], [2][0], and
// [3][0] must be initialized to f(0), f(1), f(2), and f(3) respectively.
void InitLookup16x4bx2(void* table16x4bx2);

void InitLookup16x8bx2(void* table16x8bx2);

void InitLookup256x4bx4(void* table256x4bx4);

void PhaseLookup4b(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const void* table56x4bx2, uint32_t sample_ct, void* result);

// [0][0]..[3][0], [17][0], and [19][0] should contain the relevant values
void InitPhaseLookup4b(void* table56x4bx2);

// het-haploid prohibited.  64-entry table suffices: we use the same bits for
// phasepresent and sex_male since they can't be true simultaneously.
// phaseinfo is xor'd with bits 1 and 3 instead of 1 and 2.
void PhaseXNohhLookup4b(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const uintptr_t* sex_male, const void* table64x4bx2, uint32_t sample_ct, void* result);

// [0][0]..[3][0], [16][0]..[19][0]
void InitPhaseXNohhLookup4b(void* table64x4bx2);

// uses same table as PhaseXNohhLookup
void GenoarrSexLookup4b(const uintptr_t* genoarr, const uintptr_t* sex_male, const void* table64x4bx2, uint32_t sample_ct, void* result);

// Analogue of BitIter1x.
HEADER_INLINE uint32_t GenoIter1x(const uintptr_t* __restrict genoarr, uintptr_t match_word, uintptr_t* __restrict widxp, uintptr_t* __restrict cur_bitsp) {
  uintptr_t cur_bits = *cur_bitsp;
  while (!cur_bits) {
    cur_bits = genoarr[++(*widxp)] ^ match_word;
    cur_bits = (~(cur_bits | (cur_bits >> 1))) & kMask5555;
  }
  *cur_bitsp = cur_bits & (cur_bits - 1);
  return ctzw(cur_bits);
}

// For every missing entry in genovec, clear the corresponding subset and
// sparse_vals entries.
void ClearGenovecMissing1bit8Unsafe(const uintptr_t* __restrict genovec, uint32_t* subset_sizep, uintptr_t* __restrict subset, void* __restrict sparse_vals);

void ClearGenovecMissing1bit16Unsafe(const uintptr_t* __restrict genovec, uint32_t* subset_sizep, uintptr_t* __restrict subset, void* __restrict sparse_vals);

double MultiallelicDiploidMachR2(const uint64_t* __restrict sums, const uint64_t* __restrict ssqs, uint32_t nm_sample_ct, uint32_t allele_ct);

// ----- end plink2_common subset -----

// other configuration-ish values
// this part of the specification is set in stone.

CONSTI32(kPglVblockSize, 65536);

// Currently chosen so that it plus kPglFwriteBlockSize + kCacheline - 2 is
// < 2^32, so DivUp(kPglMaxBytesPerVariant + kPglFwriteBlockSize - 1,
// kCacheline) doesn't overflow.
static const uint32_t kPglMaxBytesPerVariant = 0xfffdffc0U;
// CONSTI32(kPglMaxBytesPerDataTrack, 0x7ffff000);
// static_assert(kMaxBytesPerIO >= (int32_t)kPglMaxBytesPerDataTrack, "pgenlib_internal assumes a single variant data track always fits in one fread/fwrite operation.");

// mmap is a horrible idea for 32-bit builds, and as long as we have non-mmap
// code we may as well not worry about Win64 CreateFileMapping.

// also, OS X mmap implementation seems to be crappy for large sequentially
// accessed files, compared to Linux.

// possible todo: SIGBUS handling?  do we ever want to try to recover from an
// I/O error?
#if defined(_WIN32) || !defined(__LP64__)
#  define NO_MMAP
#endif

// currently must be power of 2, and multiple of (kBitsPerWord / 2)
CONSTI32(kPglDifflistGroupSize, 64);

FLAGSET_DEF_START()
  kfPgenGlobal0,
  kfPgenGlobalLdCompressionPresent = (1 << 0),
  kfPgenGlobalDifflistOrLdPresent = (1 << 1),

  // set opportunistically for now; may still be present if unset.
  kfPgenGlobalMultiallelicHardcallFound = (1 << 2),

  kfPgenGlobalHardcallPhasePresent = (1 << 3),
  kfPgenGlobalDosagePresent = (1 << 4),
  kfPgenGlobalDosagePhasePresent = (1 << 5),
  kfPgenGlobalAllNonref = (1 << 6)
FLAGSET_DEF_END(PgenGlobalFlags);

FLAGSET_DEF_START()
  kfPgrLdcache0,
  kfPgrLdcacheQuater = (1 << 0),
  kfPgrLdcacheDifflist = (1 << 1),
  kfPgrLdcacheRawQuater = (1 << 2),
  // may also want RawDifflist
  kfPgrLdcacheBasicGenocounts = (1 << 3)
FLAGSET_DEF_END(PgrLdcacheFlags);

// difflist/LD compression should never involve more than
//   raw_sample_ct / kPglMaxDifflistLenDivisor
// entries.  (however, returned difflists can have up to twice as many entries,
// when a variant is LD-compressed and the reference variant is
// difflist-compressed.)
// This should be considered set in stone.
CONSTI32(kPglMaxDifflistLenDivisor, 8);

// Threshold for using a deltalist to represent a bitarray on disk (currently
// relevant for dosage data).  This is a tunable parameter, but must be >=
// kPglMaxDifflistLenDivisor.
CONSTI32(kPglMaxDeltalistLenDivisor, 9);

// The actual format:
// 1. 2 magic bytes 0x6c 0x1b.
//
// 2. Mode byte.
//      0x01 = plink1 variant-major.
//      0x02 = plink2 basic variant-major.  variant/sample counts in header,
//             00 = hom ref, 01 = het, 10 = hom alt, 11 = missing.  (vrtype 0)
//      0x03 = plink2 basic unphased dosage (vrtype 0x40)
//      0x04 = plink2 basic phased dosage (vrtype 0xc0)
//      These are designed to be easy to write.  Note that the dosage formats
//      require hardcalls to be stored as well; however, you can just set them
//      to all-missing and then use
//        plink2 --hard-call-threshold [...] --make-pgen
//      to populate them.
//
//      0x10 = variable-type and/or variable-length records present.
//      0x11 = mode 0x10, but with phase set information at the end of the
//             file.
//      0x05..0x0f and 0x12..0x7f are reserved for possible use by future
//      versions of the PGEN specification, and 0 is off-limits (PLINK 1
//      sample-major .bed).
//      0x80..0xff can be safely used by developers for their own purposes.
//
// 3. If not plink1-format,
//    a. 4-byte # of variants; call this M.
//    b. 4-byte # of samples, call this N.
//    c. Additional 1-byte header 'control' value (PgenHeaderCtrl).  May be
//       extended in the future.
//       bits 0-3: Indicates vrtype and variant record length storage widths.
//         If bit 3 is unset, bits 0-1 store (vrec_len_byte_ct - 1), while bit
//         2 is set iff phase or dosage info is present (requiring 8 bits
//         instead of 4 bits for vrtypes).
//         If bit 3 is set, a specialized encoding is used which combines the
//         two pieces of information (reducing the overhead for files with few
//         samples).  The following encodings are currently defined:
//         1000: No difflist/LD/onebit compression, 2 bit
//               (vrec_len - ceil(sample_ct / 4))) value.  vrtype is zero if
//               the entry is zero, and 8 (multiallelic-hardcall) if the record
//               has 1-3 extra bytes.  Designed for single-sample files sharing
//               a single .bim-like file (note that if they don't share a .bim,
//               .bim size will dominate).
//         1001: No difflist/LD/onebit compression, 4 bit
//               (vrec_len - ceil(sample_ct / 4)) value.  vrtype is zero if the
//               entry is zero, and 8 if the record has 1-15 extra bytes.
//       bits 4-5: allele count storage (00 = unstored, 01-11 = bytes per ct)
//       bits 6-7: nonref flags info (00 = unstored, 01 = all ref/alt, 10 =
//                 never ref/alt, 11 = explicitly stored)
//       Bits 0-5 do not apply to the fixed-length modes (currently 0x02-0x04)
//       and should be zeroed out in that case.
//
// 4. If mode 0x10/0x11,
//    a. Array of 8-byte fpos values for the first variant in each vblock.
//       (Note that this suggests a way to support in-place insertions: some
//       unused space can be left between the vblocks.)
//    b. Sequence of header blocks, each containing information about
//       kPglVblockSize variants (except the last may be shorter).  All values
//       are known-width, to allow e.g. plink2 --make-pgen/--pmerge to compress
//       all variant records first, then fseek to the beginning of the output
//       file and write the header.
//         i. array of 4-bit or 1-byte vrtypes.
//        ii. array of variant record lengths (each occupying vrec_len_byte_ct
//            bytes, or 2-4 bits).
//       iii. if bits 4-5 of {3c} aren't 00, array of alt allele counts.
//        iv. nonref flags info, if explicitly stored
//      (this representation allows more efficient random access)
//    If mode 0x02-0x04, and nonref flags info explicitly stored, just that
//    bitarray.
//
// 5. The variant records.  See below for details.

// Difflist format:
//   a. [difflist_len VINT]
//   If difflist_len is zero, that's it.  Otherwise, the difflist is organized
//   into 64-element groups (the last group will usually be smaller), to make
//   extraction of e.g. a single sample less painful.  Note that with 20k
//   samples, a difflist is space-saving even with MAF 5%:
//     ~1/400 hom alt + ~38/400 het = (~39/400) * 20k
//                                  = ~1950 sample IDs.
//     that's 31 groups, requiring about 2 + 62 + 30 + 488 + 1919 = 2501 bytes
//     (can be slightly higher since a few ID deltas may be larger than 127);
//     uncompressed storage requires 5000 bytes.
//   b. [array of group start sample IDs, each of sample_id_byte_ct]
//   c. [array of 1-byte [delta segment lengths minus 63], with last entry
//      omitted]
//   d. [optional payload of fixed-width genotype values]
//      (in retrospect, this should have been positioned after (e) and omitted
//      from a lower-level packed-bitarray definition, oh well...)
//   e. one "delta segment"/group: [array of [group size - 1] VINT values,
//      each indicating the difference between the current and previous sample
//      IDs; i.e. value is 1 for two consecutive samples]

// PgenFileInfo and PgenReader are the main exported "classes".
// Exported functions involving these data structure should all have
// "pgfi"/"pgr" in their names.

// Note that this can be default-copied.
struct PgenFileInfoStruct {
  // ----- Header information, constant after initialization -----
  uint32_t raw_variant_ct;
  uint32_t raw_sample_ct;

  // 0 if variant records aren't all the same length.
  // If they are (e.g. PLINK 1 encoding; or vrtype bits 0-5 unset), we just
  // fseek to
  //   const_fpos_offset + const_vrec_width * ((uint64_t)variant_idx).
  uint64_t const_fpos_offset;

  uint32_t const_vrec_width;

  // see below.  positioned here instead of slightly later due to struct
  // packing behavior.
  uint32_t const_vrtype;  // 256 for plink 1 encoding, UINT32_MAX for nonconst

  // size (raw_variant_ct + 1), so that the number of bytes of (zero-based)
  // variant n is var_fpos[n+1] - var_fpos[n].  nullptr if
  // const_vrec_width is nonzero.
  // It's not difficult to save some memory here (e.g. unless we're dealing
  // with >256 TB files, it's trivial to go from 8 bytes down to 6 bytes per
  // entry), but I doubt that's worth the trouble; let's worry about
  // O(mn)-or-worse stuff, and on-disk stuff, first.
  uint64_t* var_fpos;

  // representation type codes.
  //
  // bits 0-2:
  //   000 = Simple 2-bit encoding.
  //   100, 110, 111 = Simple difflist.  Low two bits store the base value.
  //         (101 is currently reserved for future use since
  //         almost-all-ref/altx variants shouldn't really exist, thanks to
  //         Hardy-Weinberg equilibrium.)
  //   010 = Differences-from-earlier-variant encoding ("LD compression").  The
  //         last variant without this type of encoding is the base.
  //         To simplify random access logic, the first variant in each vblock
  //         is prohibited from using this encoding.
  //   011 = Inverted differences-from-earlier-variant encoding.  (This covers
  //         the case where a reference allele is 'wrong'.)  When decoding, the
  //         difflist should be processed first, then the entire genovec should
  //         be flipped.
  //   001 = 1-bit + difflist representation.  Suppose most calls are
  //         hom ref or het (e.g. a 20% MAF variant with ~4% hom alt1, ~36%
  //         het ref/alt1, ~64% hom ref), then the main datatrack has just the
  //         low bits of the usual 2-bit codes.  This is followed by a difflist
  //         containing the hom alt1 and missing genotypes.
  //         The main datatrack is preceded by a single byte indicating what
  //         the two common values are: 2 low bits = [set value - unset value],
  //         next 2 bits = unset value (6 possibilities).  Top 4 bits are
  //         reserved.
  // bit 3: multiallelic hardcalls present with alt2/alt3/... present?  If yes,
  //        auxiliary data track #1 disambiguates the 0b01 (ref/altx) and 0b10
  //        (altx/alty, x may equal y) hardcalls.  This contains a format byte,
  //        followed by a list of ref/altx patches, then a list of altx/alty
  //        patches.  All unpatched genotypes are ref/alt1 or alt1/alt1.
  //        The bottom 4 bits of the format byte describe how the ref/altx
  //        patch set is stored.
  //   0 = Starts with a bitarray with [total ref/altx count] bits (divide by 8
  //       and round up to get byte count of this component; any trailing bits
  //       in the last byte must be 0), where each set bit corresponds to
  //       presence of a rarealt (i.e. alt2/alt3/...).
  //       ExpandBytearr(aux1_first_quarter, all_01, raw_sample_ctl, ct_01, 0,
  //                     patch_01);
  //       can be used to convert this into a set of sample IDs, though we may
  //       want to avoid an intermediate unpacking step in practice.  Note that
  //       when we're passing in sample_ct < raw_sample_ct and the main
  //       datatrack is LD-compressed, we'd like ldbase_raw_genovec to be
  //       cached.
  //       This is followed by a packed array of fixed-width [allele idx - 2]
  //       values, where the width depends on the total number of alt alleles.
  //       2 alts: width ZERO.  All set bits in the first bitarray correspond
  //               to ref/alt2.
  //       3 alts: width 1 bit.  Set bits correspond to ref/alt3, clear bits
  //               correspond to ref/alt2.
  //       4-5 alts: width 2 bits.  0b00 corresponds to ref/alt2, 0b01 =
  //                 ref/alt3, 0b10 = ref/alt4, etc.
  //       6-17 alts: width 4 bits.
  //       18-257 alts: width 8 bits.
  //       258-65537 alts: width 16 bits.
  //       65538-16777215 alts: width 24 bits.  Reasonable to prohibit more
  //                            than 2^24 - 1 = 16777215, since variant records
  //                            are limited to 4 GiB.  I can imagine some
  //                            applications of >65534 in highly variable
  //                            regions, though, and it doesn't actually cost
  //                            us anything to define a way to represent it.
  //                            (A plink2 binary compiled with AlleleCode
  //                            typedef'd as uint32_t will run more slowly, of
  //                            course, but most binaries will not be compiled
  //                            that way.)
  //   1 = Same as mode 0, except the initial bitarray is replaced by a
  //       difflist with sample IDs.  (We could make that piece somewhat
  //       smaller by storing 0-based ref/altx indexes instead, but I'm pretty
  //       sure that isn't worth the performance penalty of requiring all_01
  //       and more complicated unpacking.  Though we'll need to peek at
  //       aux1[0] before decompressing the main datatrack to exploit this.)
  //   15 = Empty patch set.  Might remove this (storing this as mode 1 just
  //        takes 1 more byte), but it may enable some relevant performance
  //        optimizations.
  //   2-14 are reserved for future use.  We don't define an efficient way to
  //   represent a variant that e.g. has more alt2s than alt1s for now, since
  //   alt alleles will usually be sorted in order of decreasing frequency, but
  //   maybe this will matter in the future.
  //
  //   The top 4 bits describe how the altx/alty patch set is stored.  0/1/15
  //   have the same meaning as they do for the ref/altx patch set; the only
  //   thing that changes is the format of the packed array of values at the
  //   end.
  //   2 alts: width 1.  This is treated as a special case.  Set bits
  //           correspond to alt2/alt2, clear = alt1/alt2.
  //   3-4 alts: width 2+2 bits.  Each stores [allele idx - 1], with the
  //             smaller number in the lower bits.  E.g. alt1/alt2 is stored as
  //             0b0100; alt3/alt3 is stored as 0b1010.
  //   5-16 alts: width 4+4 bits.
  //   17-256 alts: width 8+8 bits.
  //   257-65536 alts: width 16+16 bits.
  //   65537-16777215 alts: width 24+24 bits.
  //
  // bit 4: hardcall phased?  If yes, auxiliary data track #2 contains
  //        phasing information for heterozygous calls.
  //        The first *bit* of the track indicates whether an explicit
  //        "phasepresent" bitarray is stored.  If it's set, the next het_ct
  //        bits are 1-bit values, where 0 = no phasing info known, and 1 =
  //        phasing info present.  If it's unset, phasing info is present for
  //        every het call.
  //        This is followed by a "phaseinfo" bitarray, where 0 = unswapped,
  //        1 = swapped (e.g. "1|0" in VCF).
  //        This track is normally unpacked into fixed-size bitarrays when
  //        loaded, but a raw mode is also provided (which doesn't support
  //        subsetting).
  //        By default, entire chromosomes/contigs are assumed to be phased
  //        together.  (Todo: support contiguous phase sets.)
  //
  // bits 5-6:
  //   00 = no dosage data.
  //   01 = dosage list.  Auxiliary data track #3 contains a delta-encoded list
  //        of sample IDs (like a difflist, but with no genotypes).  Track #4
  //        contains a 16-bit (0..2^15; 65535 missing value is only permitted
  //        in unconditional-dosage case) value expressing the sum of all alt
  //        allele dosages.  (Yes, making this the ref allele dosage would have
  //        been a bit cleaner, but it's too late now.)
  //        If the variant is multiallelic, nonzero alt2/alt3/... dosages are
  //        likely to be sparse.  So,
  //        - track #5 contains a delta-encoded list describing which
  //          [sample_uidx x rarealt dosage] entries are nonzero, where rarealt
  //          index is in the lowest bits and sample_uidx can be computed via
  //          right-shift (to avoid integer-division headaches, especially
  //          important since indexes in this list can be larger than 2^32).
  //          We use sample_uidx here to make subsetting less painful.
  //          Since each nonzero dosage value requires 16 bits, and
  //          delta-encoding of a dense list adds less than 9 bits per entry,
  //          there isn't much point in supporting a dense bitarray mode here.
  //          To keep memory requirements sane for biobank-scale datasets when
  //          the number of alt alleles is very large, each sample is
  //          prohibited from having more than 255 nonzero allele dosages (this
  //          makes no difference when sizeof(AlleleCode) == 1, but it may
  //          matter later).
  //        - track #6 contains the rarealt nonzero dosage values.
  //        Note that this and the other dosage modes are in ADDITION to
  //        hardcalls.  This increases filesize by up to 12.5%, but makes the
  //        reader substantially simpler; --hard-call-threshold logic is nicely
  //        compartmentalized.
  //   10 = unconditional dosage (just track #4).
  //   11 = dosage bitarray.  In this case, auxiliary data track #3 contains an
  //        array of 1-bit values indicating which samples have dosages.  If
  //        the variant is multiallelic, tracks #5 and 6 are as described
  //        above.
  //   bgen 1.2 format no longer permits fractional missingness, so no good
  //   reason for us to support it.
  //   Considered putting *all* dosage data at the end of the file (like I will
  //   do for phase set info); this could actually be worthwhile for
  //   unconditional dosages, but it doesn't work well when only some samples
  //   have dosage data.
  // bit 7: some dosages explicitly phased?  If yes, and dosages are not
  //        unconditionally present, auxiliary data track #7 is a bitarray of
  //        length dosage_ct indicating whether dphase_delta exists for that
  //        sample.  Note that this is technically independent of bit 4; either
  //        can be set without the other.  (However, in practice, bit 4 is
  //        almost always set when bit 7 is, since that enables more efficient
  //        storage of 0|0.99, 1|0.02, and similar pairs.)
  //        When phased dosages are present, track #8 contains values
  //        representing [(hap1 alt prob) - (hap2 alt prob)], etc., where the
  //        underlying values are represented in [0..16384] (so the signed
  //        difference is in [-16384..16384]).  Track #4 contains the
  //        corresponding sums; parity does NOT need to match (necessary to
  //        allow this since we treat omitted values as zero; and since we are
  //        allowing it, no point in making parity match in other situations
  //        either).  In fixed-width case, -32768 should be stored in track #8
  //        when the entire call is missing, while 0 and missing-phase are
  //        considered synonymous.
  //        In the biallelic case, if a hardcall is phased, a dosage is
  //        present, and no explicit dosage-phase is, we define it to mean the
  //        unique dphase_delta sequence with maximal absolute value, and
  //        --make-pgen takes advantage of it.  This definition technically
  //        works for triallelic variants as well, but it breaks down with 4
  //        alleles, so we prohibit hardcall-phase + dosage + no-dosage-phase
  //        with more than 2 alleles.
  //        In the multiallelic case, tracks #9 and #10 are analogous to #5 and
  //        #6.
  //
  // Representation of variable ploidy (MT) was considered, but rejected since
  // dosages should be at least as appropriate for MT.
  // Oxford/VCF-style storage of separate probabilities for every possible
  // genotype (e.g. P(AA), P(AB), P(BB) instead of just 2P(AA) + P(AB) and
  // 2P(BB) + P(AB)) is tentatively rejected due to (i) lack of relevance to
  // PLINK's analysis functions and (ii) high storage cost where we can afford
  // it least.  In principle, this is subject to reevaluation if (i) changes,
  // but given the poor interaction with phased dosages, it's probably better
  // to just think of them as permanently outside PLINK's scope.
  //
  //
  // base pointer is null if mode is 0x01-0x04 (const_vrtype != UINT32_MAX).
  // if not nullptr, required to be length >=
  //   max(raw_variant_ct + 1, RoundUpPow2(raw_variant_ct, kBytesPerWord))
  unsigned char* vrtypes;

  // alt allele counts.

  // This can be nullptr if all alt allele counts are 1.
  // (actually, we store the allele index offsets, so
  // (allele_idx_offsets[n+1] - allele_idx_offsets[n]) is the number of alleles
  // for variant n.  Otherwise, we'd need another data structure to support
  // fast allele name lookup.)
  uintptr_t* allele_idx_offsets;

  uintptr_t* nonref_flags;

  // If pgr.nonref_flags is nullptr and kfPgenGlobalAllNonref is unset, all
  // reference alleles are assumed to be correct.
  PgenGlobalFlags gflags;

  uint32_t max_allele_ct;
  // uint32_t max_dosage_allele_ct;  // might need this later

  // * nullptr if using mmap
  // * if using per-variant fread(), this is non-null during Pgen_file_info
  //   initialization, but it's then "moved" to the first Pgen_reader and set
  //   to nullptr.
  FILE* shared_ff;

  const unsigned char* block_base;  // nullptr if using per-variant fread()
  uint64_t block_offset;  // 0 for mmap
#ifndef NO_MMAP
  uint64_t file_size;
#endif
};

typedef struct PgenFileInfoStruct PgenFileInfo;

struct PgenReaderStruct {
  NONCOPYABLE(PgenReaderStruct);
  // would like to make this const, but that makes initialization really
  // annoying in C99
  struct PgenFileInfoStruct fi;

  // ----- Mutable state -----
  // If we don't fseek, what's the next variant we'd read?  (Still relevant
  // with mmap due to how LD decompression is implemented.)
  uint32_t fp_vidx;

  // ** per-variant fread()-only **
  FILE* ff;
  unsigned char* fread_buf;
  // ** end per-variant fread()-only **

  // if LD compression is present, cache the last non-LD-compressed variant
  uint32_t ldbase_vidx;

  // flags indicating which base_variant buffers are populated
  PgrLdcacheFlags ldbase_stypes;

  uint32_t ldbase_difflist_len;

  // these should be treated as private after initial allocation.
  // not currently guaranteed to have trailing zeroes.
  uintptr_t* ldbase_raw_genovec;  // now allocated even with no LD compression
  uintptr_t* ldbase_genovec;
  uintptr_t* ldbase_raregeno;

  // when ldbase_difflist_ids[] is initialized, element [ldbase_difflist_len]
  // must be set to sample_ct.
  uint32_t* ldbase_difflist_sample_ids;

  // common genotype can be looked up from vrtypes[]

  STD_ARRAY_DECL(uint32_t, 4, ldbase_basic_genocounts);

  // now only allocated if multiallelic variants, phase, and/or dosage present
  // most commonly used for unsubsetted genovec; all_hets can be computed from
  // this and patch_10_{set,vals}, and then aux2 can be interpreted.
  // can also be used for other purposes after we're done processing aux2.
  uintptr_t* workspace_vec;

  // currently must hold (raw_sample_ct / kPglMaxDifflistLenDivisor)
  // entries; may need to double the sizes later
  // some top-level interface functions use these, so several lower-level
  // functions cannot
  uintptr_t* workspace_raregeno_vec;
  uint32_t* workspace_difflist_sample_ids;

  // must hold (raw_sample_ct / kPglMaxDifflistLenDivisor) entries
  uintptr_t* workspace_raregeno_tmp_loadbuf;

  uintptr_t* workspace_aux1x_present;
  uint64_t* workspace_mach_r2;  // needed in multiallelic case
  uintptr_t* workspace_all_hets;
  uintptr_t* workspace_subset;  // currently used for hphase decoding

  uintptr_t* workspace_dosage_present;
  uintptr_t* workspace_dphase_present;

  // phase set loading (mode 0x11) unimplemented for now; should be a sequence
  // of (sample ID, [uint32_t phase set begin, set end), [set begin, set end),
  // ...).
};

typedef struct PgenReaderStruct PgenReader;

// might want this value to be typed...
CONSTI32(kPglVrtypePlink1, 256);

HEADER_INLINE uint32_t GetPgfiVrtype(const PgenFileInfo* pgfip, uint32_t vidx) {
  if (pgfip->vrtypes) {
    return pgfip->vrtypes[vidx];
  }
  return pgfip->const_vrtype;
}

HEADER_INLINE uint64_t GetPgfiFpos(const PgenFileInfo* pgfip, uintptr_t vidx) {
  if (pgfip->var_fpos) {
    return pgfip->var_fpos[vidx];
  }
  return pgfip->const_fpos_offset + pgfip->const_vrec_width * S_CAST(uint64_t, vidx);
}

HEADER_INLINE uint32_t GetPgfiVrecWidth(const PgenFileInfo* pgfip, uint32_t vidx) {
  if (pgfip->var_fpos) {
    return pgfip->var_fpos[vidx + 1] - pgfip->var_fpos[vidx];
  }
  return pgfip->const_vrec_width;
}

HEADER_INLINE uint32_t PgfiIsSimpleFormat(const PgenFileInfo* pgfip) {
  return (pgfip->const_vrtype != UINT32_MAX);
}

HEADER_INLINE uint32_t VrtypeDifflist(uint32_t vrtype) {
  return (vrtype & 4);
}

HEADER_INLINE uint32_t VrtypeLdCompressed(uint32_t vrtype) {
  return (vrtype & 6) == 2;
}

// Only checks for rarealt-containing hardcall.  Multiallelic dosage may still
// be present when this returns zero.
HEADER_INLINE uint32_t VrtypeMultiallelicHc(uint32_t vrtype) {
  return (vrtype & 8);
}

HEADER_INLINE uint32_t VrtypeHphase(uint32_t vrtype) {
  return (vrtype & 0x10);
}

HEADER_INLINE uint32_t VrtypeAuxTracksPresent(uint32_t vrtype) {
  return (vrtype & 0x78);
}

HEADER_INLINE uint32_t VrtypeVariableWidth(uint32_t vrtype) {
  return (vrtype & 0x3e);
}

HEADER_INLINE uint32_t VrtypeDosage(uint32_t vrtype) {
  return (vrtype & 0x60);
}

static_assert(kPglMaxAltAlleleCt <= 254, "GetAux1xAlleleEntryByteCt() needs to be updated.");
HEADER_INLINE uintptr_t GetAux1aAlleleEntryByteCt(uint32_t allele_ct, uint32_t rare01_ct) {
  assert(allele_ct >= 3);
  if (allele_ct == 3) {
    return 0;
  }
  if (allele_ct == 4) {
    return DivUp(rare01_ct, 8);
  }
  if (allele_ct <= 6) {
    return DivUp(rare01_ct, 4);
  }
  if (allele_ct <= 18) {
    return DivUp(rare01_ct, 2);
  }
  return rare01_ct;
}

HEADER_INLINE uintptr_t GetAux1bAlleleEntryByteCt(uint32_t allele_ct, uint32_t rare10_ct) {
  assert(allele_ct >= 3);
  if (allele_ct == 3) {
    return DivUp(rare10_ct, 8);
  }
  if (allele_ct < 6) {
    return DivUp(rare10_ct, 2);
  }
  // one byte per entry for allele_ct <= 17, two bytes for 18..256
  return ((allele_ct >= 18) + 1) * rare10_ct;
  // todo: allele_ct > 257
}

// PgenFileInfo initialization is split into two phases, to decouple
// plink2's arena allocator from this library.
//
// Phase 1: Open the .pgen; verify that the initial bytes are consistent with
//   the file format; load/verify sample and variant counts, initialize
//   pgfi.const_vrtype, pgfi.const_vrec_width, and pgfi.const_fpos_offset;
//   determine initial memory allocation requirement.  first_alloc_cacheline_ct
//   does not include allele counts and nonref flags, since it may be more
//   appropriate to allocate those arrays earlier (during loading of a
//   .bim-like file).
//
//   pgfi.var_fpos is set to nullptr if pgfi.const_vrec_width is nonzero.
//   pgfi.vrtypes/var_allele_cts are set to nullptr in the plink1-format case.
//
//   raw_sample_ct and raw_variant_ct should be UINT32_MAX if not previously
//   known.
//
// Intermission: Caller obtains a block of pgfi_alloc_cacheline_ct * 64 bytes,
//   64-byte aligned.  The cachealigned_malloc() function can be used for this
//   purpose.  If necessary, pgfi.allele_idx_offsets and pgfi.nonref_flags
//   should be pointed at already-loaded data, or allocated so they can be
//   loaded during phase 2.
//
// Phase 2: Initialize most pointers in the PgenReader struct to appropriate
//   positions in first_alloc.  For modes 0x10-0x11, load pgfi.var_fpos and
//   pgfi.vrtypes, load/validate pgfi.allele_idx_offsets and pgfi.nonref_flags
//   if appropriate, and initialize pgfi.gflags, pgfi.max_allele_ct, and
//   pgfi.max_dosage_allele_ct.
//
// Finally, if block-fread mode is being used, pgfi.block_base must be
//   initialized to point to a memory large enough to handle the largest
//   pgfi_block_read() operation that will be attempted.
//   pgfi_blockload_get_cacheline_req() can be used to determine the necessary
//   buffer size.

// This type may change if we introduce a more read-optimized format in the
// future.  Right now it just tracks the presence/absence of two optional
// pieces of information: allele counts and nonref flags.
typedef uint32_t PgenHeaderCtrl;

void PreinitPgfi(PgenFileInfo* pgfip);

// There are three modes of operation:
// 1. mmaped file.  Appropriate for handling multiple queries across different
//    parts of the genome in parallel.  Suboptimal for whole-genome queries.
//    Doesn't currently run on Windows.
// 2. fread block-load.  Block-load operations are single-threaded, while
//    decompression/counting is multithreaded.  Appropriate for whole-genome
//    queries, since even with a SSD, reading from multiple parts of a file
//    simultaneously doesn't work well.
// 3. fread single-variant-at-a-time.  Simpler interface than block-load, and
//    doesn't share its inability to handle multiple queries at a time, but
//    less performant for CPU-heavy operations on the whole genome.
//
// To specify mode 1, pass in use_mmap == 1 here.
// To specify mode 2, pass in use_mmap == 0 here, and use_blockload == 1 during
//   phase2.
// To specify mode 3, pass in use_mmap == 0 here, and use_blockload == 0 during
//   phase2.
//
// Update (7 Jan 2018): raw_variant_ct must be in [1, 2^31 - 3], and
//   raw_sample_ct must be in [1, 2^31 - 2].
PglErr PgfiInitPhase1(const char* fname, uint32_t raw_variant_ct, uint32_t raw_sample_ct, uint32_t use_mmap, PgenHeaderCtrl* header_ctrl_ptr, PgenFileInfo* pgfip, uintptr_t* pgfi_alloc_cacheline_ct_ptr, char* errstr_buf);

// If allele_cts_already_loaded is set, but they're present in the file,
// they'll be validated; similarly for nonref_flags_already_loaded.
PglErr PgfiInitPhase2(PgenHeaderCtrl header_ctrl, uint32_t allele_cts_already_loaded, uint32_t nonref_flags_already_loaded, uint32_t use_blockload, uint32_t vblock_idx_start, uint32_t vidx_end, uint32_t* max_vrec_width_ptr, PgenFileInfo* pgfip, unsigned char* pgfi_alloc, uintptr_t* pgr_alloc_cacheline_ct_ptr, char* errstr_buf);


uint64_t PgfiMultireadGetCachelineReq(const uintptr_t* variant_include, const PgenFileInfo* pgfip, uint32_t variant_ct, uint32_t block_size);

// variant_include can be nullptr; in that case, we simply load all the
// variants (load_variant_ct must be variant_uidx_end - variant_uidx_start).)
// IMPORTANT: pgfi.block_offset must be manually copied to each reader for now.
//   (todo: probably replace pgr.fi with a pointer.  when doing that, need to
//   ensure multiple per-variant readers still works.)
PglErr PgfiMultiread(const uintptr_t* variant_include, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t load_variant_ct, PgenFileInfo* pgfip);


void PreinitPgr(PgenReader* pgrp);

// Before PgrInit() is called, the caller must obtain a block of
// pgr_alloc_cacheline_ct * 64 bytes (this value is returned by
// pgfi_init_phase2), 64-byte aligned; this is the pgr_alloc parameter.
//
// There's also a modal usage difference:
//
// * Modes 1-2 (mmap, block-fread): There is one PgenFileInfo per file
//   which doesn't belong to any reader.  After it's initialized, multiple
//   PgenReaders can be based off of it.  When the PgenFileInfo is
//   destroyed, those PgenReaders are invalidated and should be destroyed if
//   that hasn't already happened.
//
//   fname parameter must be nullptr.
//
// * Mode 3 (per-variant fread): Destruction of the original PgenFileInfo
//   struct does not invalidate any extant PgenReader instances (at least
//   from pgenlib_internal's perspective).  Instead, destruction of the
//   corresponding memory block or allele_idx_offsets/nonref_flags invalidates
//   the associated PgenReaders.
//
//   The only difference between the first reader and later readers of the same
//   file is that the first reader steals the shared_ff used to read the
//   header.
//
//   fname parameter must be non-null.

// max_vrec_width ignored when using mode 1 or 2.
PglErr PgrInit(const char* fname, uint32_t max_vrec_width, PgenFileInfo* pgfip, PgenReader* pgrp, unsigned char* pgr_alloc);

// practically all these functions require genovec to be allocated up to
// vector, not word, boundary
void PgrPlink1ToPlink2InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec);

void PgrPlink2ToPlink1InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec);

void PgrDifflistToGenovecUnsafe(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uintptr_t difflist_common_geno, uint32_t sample_ct, uint32_t difflist_len, uintptr_t* __restrict genovec);

// Function names for the main reader functions were getting ridiculous.
// New naming scheme:
// * PgrGet() is the basic two-bit genovec loader.  All ALT alleles are treated
//   as equivalent.  (00 = hom ref, 01 = het ref, 10 = two alt alleles, 11 =
//   missing.)
// * PgrGetInv1() is similar, except that the allele index to treat as REF can
//   be changed.
// * PgrGet1() only counts the specified allele.  To minimize inversion costs,
//   GetInv1() should be called on major alleles and Get1() should be called on
//   minor ones.
// * PgrGetM() is the multiallelic loader which doesn't collapse multiple
//   alleles into one.  This retrieves a sparse form identical to what
//   PwcAppendMultiallelicSparse takes.
//   Multiallelic-dosage read functions (PgrReadRaw() included) will probably
//   fill a 3-part data structure of the following form:
//   1. Bitarray indicating which samples have at least one rarealt dosage.
//   2. unsigned char array where, if bits a, b, and c are the only set ones in
//      the first array, the first three elements of the second array are
//      rarealt dosage counts (1..255) for those three samples.  (Could also
//      put those in positions [a], [b], and [c], but that produces worse
//      memory access locality, and it makes sense to treat multiallelic
//      dosages as fundamentally sparse.)
//   3. Let R := MINV(255, allele_ct - 2).
//      a. Length-(sample_ct x R) array of AlleleCodes.
//      b. Length-(sample_ct x R) array of uint16_t dosage (or int16_t dphase)
//         values.
//      Again we use the sparse representation, with payload values packed at
//      the beginning.
//   (--indiv-sort algorithm: initialize an array of uintptr_ts of length
//   sample_ct where [k] has that sample's start index in the payload arrays.)
// * PgrGetDifflistOrGenovec() opportunistically returns the sparse genotype
//   representation ('difflist'), for functions capable of taking advantage of
//   it.  I don't plan to use this in plink2 before at least 2019, but the
//   pgen_compress demo program illustrates its usage.
// * PgrGetCounts() is equivalent to calling PgrGet() and then counting the
//   number of 00s, 01s, 10s, and 11s, without the overhead of fully expanding
//   the compressed data, etc.
// * P suffix = also returns hardcall-phase information.
// * D suffix = also returns dosage information.
// * Dp suffix = also returns hardcall-phase, dosage and phased-dosage
//   information.
// * PgrGet2() and PgrGet2P() loads biallelic (possibly phased) hardcalls from
//   a possibly-multiallelic variant.  Any hardcall where either allele is not
//   one of the specified two alleles is set to missing.
//   There is no dosage-supporting version of this because rescaling sucks.

// This will normally extract only the genotype indexes corresponding to set
// bits in sample_include.  Set sample_ct == raw_sample_ct if you don't want
// any subsetting to occur (in this case sample_include is ignored, can be
// nullptr).
// sample_ct cannot be zero.  Trailing bits of genovec are not zeroed out.
// Ok if genovec only has space for sample_ct values.
PglErr PgrGet(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec);

// Loads the specified variant as a difflist if that's more efficient, setting
// difflist_common_geno to the common genotype value in that case.  Otherwise,
// genovec is populated and difflist_common_geno is set to UINT32_MAX.
//
// Note that the returned difflist_len can be much larger than
// max_simple_difflist_len when the variant is LD-encoded; it's bounded by
//   2 * (raw_sample_ct / kPglMaxDifflistLenDivisor).
PglErr PgrGetDifflistOrGenovec(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t max_simple_difflist_len, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uint32_t* difflist_common_geno_ptr, uintptr_t* __restrict main_raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr);

// This is necessary when changing sample_include, unless the new query is
// iterating from the first variant.  (Which can almost never be assumed in
// plink2 since variant_include[] may not include the first variant.)
HEADER_INLINE void PgrClearLdCache(PgenReader* pgrp) {
  pgrp->ldbase_stypes &= kfPgrLdcacheRawQuater;

  // bugfix, LdLoadNecessary() was otherwise claiming that reload wasn't
  // necessary in certain cases
  pgrp->ldbase_vidx = 0x80000000U;
}

// genocounts[0] = # hom ref, [1] = # het ref, [2] = two alts, [3] = missing
PglErr PgrGetCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, STD_ARRAY_REF(uint32_t, 4) genocounts);

// genocounts[0] = # of hardcalls with two copies of specified allele
// genocounts[1] = # of hardcalls with exactly one copy of specified allele
// genocounts[2] = # of hardcalls with no copies
// genocounts[3] = missing
PglErr PgrGetInv1Counts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, STD_ARRAY_REF(uint32_t, 4) genocounts);

// Loads a quatervec with counts of a single allele (allele_idx 0 corresponds
// to the reference allele, allele_idx 1 corresponds to alt1, etc.).  0b11 ==
// missing call.
// Note that calling this with allele_idx == 0 is similar to a plink1 load
// (except with missing == 0b11, of course).
// todo: provide a difflist interface once anyone wants it.
PglErr PgrGet1(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_countvec);

PglErr PgrGetInv1(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_invcountvec);

PglErr PgrGet2(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx0, uint32_t allele_idx1, PgenReader* pgrp, uintptr_t* __restrict genovec);

// This covers all the possibilities.  Todo: switch all functions exporting
// multiallelic codes and/or phased dosage to use this struct.  (Biallelic
// phased hardcalls and unphased dosages are simple enough for this to be
// overkill, though.)
struct PgenVariantStruct {
  NONCOPYABLE(PgenVariantStruct);
  uintptr_t* genovec;
  uintptr_t* patch_01_set;
  AlleleCode* patch_01_vals;
  uintptr_t* patch_10_set;
  AlleleCode* patch_10_vals;
  uintptr_t* phasepresent;
  uintptr_t* phaseinfo;
  uintptr_t* dosage_present;
  uint16_t* dosage_main;
  uintptr_t* multidosage_present;
  unsigned char* multidosage_cts;
  AlleleCode* multidosage_codes;
  uint16_t* multidosage_vals;
  uintptr_t* dphase_present;
  int16_t* dphase_delta;
  uintptr_t* multidphase_present;
  unsigned char* multidphase_cts;
  AlleleCode* multidphase_codes;
  int16_t* multidphase_delta;

  uint32_t patch_01_ct;
  uint32_t patch_10_ct;
  uint32_t phasepresent_ct;
  uint32_t dosage_ct;
  uint32_t multidosage_sample_ct;
  uint32_t dphase_ct;
  uint32_t multidphase_sample_ct;
};

typedef struct PgenVariantStruct PgenVariant;

void PreinitPgv(PgenVariant* pgvp);

PglErr PgrGetM(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp);

// possible todo: add functions which directly support MAF-based queries.  Note
// that when the difflist representation is used, we can disqualify some
// low-MAF variants without actually loading the genotype data, since the size
// of the record puts an upper bound on the alt allele frequency.

// requires trailing bits of genovec to be zeroed out, AND does not update high
// bits of last word if raw_sample_ctl2 is odd.
void DetectGenovecHetsHw(const uintptr_t*__restrict genovec, uint32_t raw_sample_ctl2, Halfword* __restrict all_hets_hw);

// requires trailing bits of genovec to be zeroed out.
HEADER_INLINE void PgrDetectGenovecHetsUnsafe(const uintptr_t*__restrict genovec, uint32_t raw_sample_ctl2, uintptr_t* __restrict all_hets) {
  Halfword* all_hets_alias = R_CAST(Halfword*, all_hets);
  DetectGenovecHetsHw(genovec, raw_sample_ctl2, all_hets_alias);
  if (raw_sample_ctl2 % 2) {
    all_hets_alias[raw_sample_ctl2] = 0;
  }
}

HEADER_INLINE void PgrDetectGenovecHets(const uintptr_t* __restrict genovec, uint32_t raw_sample_ct, uintptr_t* __restrict all_hets) {
  DetectGenovecHetsHw(genovec, QuaterCtToWordCt(raw_sample_ct), R_CAST(Halfword*, all_hets));
  ZeroTrailingBits(raw_sample_ct, all_hets);
}

// sample_ct > 0.  ok for trailing bits of genovec to not be zeroed out.
void PgrDetectGenovecHetsMultiallelic(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t raw_sample_ct, uintptr_t* __restrict all_hets);

// cannot assume phaseinfo bit is clear when phasepresent is clear.
PglErr PgrGetP(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr);

PglErr PgrGet1P(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_countvec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr);

PglErr PgrGetInv1P(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_invcountvec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr);

PglErr PgrGet2P(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx0, uint32_t allele_idx1, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr);

PglErr PgrGetMP(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp);

// if dosage_present and dosage_main are nullptr, dosage data is ignored
PglErr PgrGetD(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr);

PglErr PgrGet1D(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_countvec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr);

PglErr PgrGetInv1D(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_invcountvec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr);

PglErr PgrGetDCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t is_minimac3_r2, PgenReader* pgrp, double* imp_r2_ptr, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* __restrict all_dosages);

PglErr PgrGetMDCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t is_minimac3_r2, PgenReader* pgrp, double* __restrict imp_r2_ptr, uint32_t* __restrict het_ctp, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* __restrict all_dosages);

PglErr PgrGetMD(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp);

// ok for both dosage_present and dosage_main to be nullptr when no dosage data
// is present
// ok for dphase_present/dphase_delta to be nullptr; dphase_ct always set to 0
// in that case
PglErr PgrGetDp(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp);

// pgvp->genovec filled with inverse-counts for specified allele
PglErr PgrGetInv1Dp(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReader* pgrp, PgenVariant* pgvp);

PglErr PgrGetMDp(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, PgenVariant* pgvp);

// interface used by --make-pgen, just performs basic LD/difflist decompression
// to maximize parallelism
PglErr PgrGetRaw(uint32_t vidx, PgenGlobalFlags read_gflags, PgenReader* pgrp, uintptr_t** loadbuf_iter_ptr, unsigned char* loaded_vrtype_ptr);

PglErr PgrValidate(PgenReader* pgrp, uintptr_t* genovec_buf, char* errstr_buf);

// missingness bit is set iff hardcall is not present (even if dosage info *is*
// present)
PglErr PgrGetMissingness(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict missingness, uintptr_t* __restrict genovec_buf);

// either missingness_hc (hardcall) or missingness_dosage must be non-null for
// now
// missingness_dosage must be vector-aligned
PglErr PgrGetMissingnessD(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict missingness_hc, uintptr_t* __restrict missingness_dosage, uintptr_t* __restrict hets, uintptr_t* __restrict genovec_buf);


// failure = kPglRetReadFail
BoolErr CleanupPgfi(PgenFileInfo* pgfip);

BoolErr CleanupPgr(PgenReader* pgrp);


struct PgenWriterCommonStruct {
  NONCOPYABLE(PgenWriterCommonStruct);
  uint32_t variant_ct;
  uint32_t sample_ct;
  PgenGlobalFlags phase_dosage_gflags;  // subset of gflags

  // there should be a single copy of these arrays shared by all threads.
  // allele_idx_offsets is read-only.
  uint64_t* vblock_fpos;
  unsigned char* vrec_len_buf;
  uintptr_t* vrtype_buf;
  const uintptr_t* allele_idx_offsets;
  uintptr_t* explicit_nonref_flags;  // usually nullptr

  // needed for multiallelic-phased case
  uintptr_t* genovec_hets_buf;

  STD_ARRAY_DECL(uint32_t, 4, ldbase_genocounts);

  // should match ftello() return value in singlethreaded case, but be set to
  // zero in multithreaded case
  uint64_t vblock_fpos_offset;

  // these must hold sample_ct entries
  // genovec_invert_buf also used as phaseinfo and dphase_present temporary
  // storage
  uintptr_t* genovec_invert_buf;
  uintptr_t* ldbase_genovec;

  // these must hold 2 * (sample_ct / kPglMaxDifflistLenDivisor) entries
  uintptr_t* ldbase_raregeno;
  uint32_t* ldbase_difflist_sample_ids;  // 1 extra entry, == sample_ct

  // this must fit 64k variants in multithreaded case
  unsigned char* fwrite_buf;
  unsigned char* fwrite_bufp;

  uint32_t ldbase_common_geno;  // UINT32_MAX if ldbase_genovec present
  uint32_t ldbase_difflist_len;

  // I'll cache this for now
  uintptr_t vrec_len_byte_ct;

  uint32_t vidx;
};

typedef struct PgenWriterCommonStruct PgenWriterCommon;

CONSTI32(kPglFwriteBlockSize, 131072);

// Given packed arrays of unphased biallelic genotypes in uncompressed plink2
// binary format (00 = hom ref, 01 = het ref/alt1, 10 = hom alt1, 11 =
// missing), {Single,Multi}threaded_pgen_writer performs difflist (sparse
// variant), one bit (mostly-two-value), and LD compression before writing to
// disk, and backfills the header at the end.  CPRA -> CPR merging is under
// development.
// The major difference between the two interfaces is that
// Multithreaded_pgen_writer forces you to process large blocks of variants at
// a time (64k per thread).  So Singlethreaded_pgen_writer is still worth using
// in some cases (memory is very limited, I/O is slow, no programmer time to
// spare for the additional complexity).

struct STPgenWriterStruct {
  NONCOPYABLE(STPgenWriterStruct);
  struct PgenWriterCommonStruct pwc;
  FILE* pgen_outfile;
};

struct MTPgenWriterStruct {
  NONCOPYABLE(MTPgenWriterStruct);
  FILE* pgen_outfile;
  uint32_t thread_ct;
  struct PgenWriterCommonStruct* pwcs[];
};

typedef struct STPgenWriterStruct STPgenWriter;
typedef struct MTPgenWriterStruct MTPgenWriter;

void PreinitSpgw(STPgenWriter* spgwp);

// phase_dosage_gflags zero vs. nonzero is most important: this determines size
// of header.  Otherwise, setting more flags than necessary just increases
// memory requirements.
//
// nonref_flags_storage values:
//   0 = no info stored
//   1 = always trusted
//   2 = always untrusted
//   3 = use explicit_nonref_flags
PglErr SpgwInitPhase1(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, STPgenWriter* spgwp, uintptr_t* alloc_cacheline_ct_ptr, uint32_t* max_vrec_len_ptr);

void SpgwInitPhase2(uint32_t max_vrec_len, STPgenWriter* spgwp, unsigned char* spgw_alloc);

// moderately likely that there isn't enough memory to use the maximum number
// of threads, so this returns per-thread memory requirements before forcing
// the caller to specify thread count
// (eventually should write code which falls back on STPgenWriter
// when there isn't enough memory for even a single 64k variant block, at least
// for the most commonly used plink 2.0 functions)
void MpgwInitPhase1(const uintptr_t* __restrict allele_idx_offsets, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uintptr_t* alloc_base_cacheline_ct_ptr, uint64_t* alloc_per_thread_cacheline_ct_ptr, uint32_t* vrec_len_byte_ct_ptr, uint64_t* vblock_cacheline_ct_ptr);

PglErr MpgwInitPhase2(const char* __restrict fname, const uintptr_t* __restrict allele_idx_offsets, uintptr_t* __restrict explicit_nonref_flags, uint32_t variant_ct, uint32_t sample_ct, PgenGlobalFlags phase_dosage_gflags, uint32_t nonref_flags_storage, uint32_t vrec_len_byte_ct, uintptr_t vblock_cacheline_ct, uint32_t thread_ct, unsigned char* mpgw_alloc, MTPgenWriter* mpgwp);


// trailing bits of genovec must be zeroed out
void PwcAppendBiallelicGenovec(const uintptr_t* __restrict genovec, PgenWriterCommon* pwcp);

BoolErr SpgwFlush(STPgenWriter* spgwp);

HEADER_INLINE PglErr SpgwAppendBiallelicGenovec(const uintptr_t* __restrict genovec, STPgenWriter* spgwp) {
  if (unlikely(SpgwFlush(spgwp))) {
    return kPglRetWriteFail;
  }
  PwcAppendBiallelicGenovec(genovec, &(spgwp->pwc));
  return kPglRetSuccess;
}

// trailing bits of raregeno must be zeroed out
// all raregeno entries assumed to be unequal to difflist_common_geno; the
// difflist should be compacted first if this isn't true
// difflist_len must be <= 2 * (sample_ct / kPglMaxDifflistLenDivisor);
// there's an assert checking this
void PwcAppendBiallelicDifflistLimited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, PgenWriterCommon* pwcp);

HEADER_INLINE PglErr SpgwAppendBiallelicDifflistLimited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, STPgenWriter* spgwp) {
  if (unlikely(SpgwFlush(spgwp))) {
    return kPglRetWriteFail;
  }
  PwcAppendBiallelicDifflistLimited(raregeno, difflist_sample_ids, difflist_common_geno, difflist_len, &(spgwp->pwc));
  return kPglRetSuccess;
}

// Two interfaces for appending multiallelic hardcalls:
// 1. sparse: genovec, bitarray+values describing ref/altx hardcalls which
//    aren't ref/alt1, bitarray+values describing altx/alty hardcalls which
//    aren't alt1/alt1.
//    patch_01_vals[] contains one entry per relevant sample (only need altx
//    index), patch_10_vals[] contains two.
//    Ok if patch_01_ct == patch_10_ct == 0; in this case no aux1 track is
//    saved and bit 3 of vrtype is not set.  (Note that multiallelic dosage may
//    still be present when vrtype bit 3 is unset.)
// 2. generic dense: takes a length-2n array of AlleleCode allele codes.
//    Assumes [2k] <= [2k+1] for each k.  Instead of providing direct API
//    functions for this, we just provide a dense -> sparse helper function.
//
// All arrays should be vector-aligned.
BoolErr PwcAppendMultiallelicSparse(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t patch_01_ct, uint32_t patch_10_ct, PgenWriterCommon* pwcp);

HEADER_INLINE PglErr SpgwAppendMultiallelicSparse(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t patch_01_ct, uint32_t patch_10_ct, STPgenWriter* spgwp) {
  if (unlikely(SpgwFlush(spgwp))) {
    return kPglRetWriteFail;
  }
  if (unlikely(PwcAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, patch_01_ct, patch_10_ct, &(spgwp->pwc)))) {
    return kPglRetVarRecordTooLarge;
  }
  return kPglRetSuccess;
}

// This may not zero out trailing halfword of patch_{01,10}_set.
void PglMultiallelicDenseToSparse(const AlleleCode* __restrict wide_codes, uint32_t sample_ct, uintptr_t* __restrict genovec, uintptr_t* __restrict patch_01_set, AlleleCode* __restrict patch_01_vals, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals, uint32_t* __restrict patch_01_ct_ptr, uint32_t* __restrict patch_10_ct_ptr);

// If remap is not nullptr, this simultaneously performs a rotation operation:
// wide_codes[2n] and [2n+1] are set to remap[geno[n]] rather than geno[n], and
// bit n of flipped (if flipped non-null) is set iff phase orientation is
// flipped (i.e. wide_codes[2n] was larger than wide_codes[2n+1] before the
// final reordering pass; caller needs to know this to properly update
// phaseinfo, dphase_delta, multidphase_delta).
// It currently assumes no alleles are being mapped to 'missing'.
void PglMultiallelicSparseToDense(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, const AlleleCode* __restrict remap, uint32_t sample_ct, uint32_t patch_01_ct, uint32_t patch_10_ct, uintptr_t* __restrict flipped, AlleleCode* __restrict wide_codes);

// Permits missing codes, does not remap.
void PglMultiallelicSparseToDenseMiss(const PgenVariant* pgvp, uint32_t sample_ct, AlleleCode* __restrict wide_codes);

// phasepresent == nullptr ok, that indicates that ALL heterozygous calls are
// phased.  Caller should use e.g. PwcAppendBiallelicGenovec() if it's known
// in advance that no calls are phased.
// Ok for phaseinfo to have bits set at non-het calls, NOT currently okay for
//   phasepresent
void PwcAppendBiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, PgenWriterCommon* pwcp);

// phasepresent == nullptr ok
// ok for trailing bits of phaseinfo to not be zeroed out, NOT currently ok for
//   phasepresent
HEADER_INLINE PglErr SpgwAppendBiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, STPgenWriter* spgwp) {
  if (unlikely(SpgwFlush(spgwp))) {
    return kPglRetWriteFail;
  }
  PwcAppendBiallelicGenovecHphase(genovec, phasepresent, phaseinfo, &(spgwp->pwc));
  return kPglRetSuccess;
}

BoolErr PwcAppendMultiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, uint32_t patch_01_ct, uint32_t patch_10_ct, PgenWriterCommon* pwcp);

HEADER_INLINE PglErr SpgwAppendMultiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict patch_01_set, const AlleleCode* __restrict patch_01_vals, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, uint32_t patch_01_ct, uint32_t patch_10_ct, STPgenWriter* spgwp) {
  if (unlikely(SpgwFlush(spgwp))) {
    return kPglRetWriteFail;
  }
  if (unlikely(PwcAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, phasepresent, phaseinfo, patch_01_ct, patch_10_ct, &(spgwp->pwc)))) {
    return kPglRetVarRecordTooLarge;
  }
  return kPglRetSuccess;
}


// dosage_main[] has length dosage_ct, not sample_ct
// ok for traling bits of dosage_present to not be zeroed out
BoolErr PwcAppendBiallelicGenovecDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, PgenWriterCommon* pwcp);

HEADER_INLINE PglErr SpgwAppendBiallelicGenovecDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, STPgenWriter* spgwp) {
  if (unlikely(SpgwFlush(spgwp))) {
    return kPglRetWriteFail;
  }
  if (unlikely(PwcAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_main, dosage_ct, &(spgwp->pwc)))) {
    return kPglRetVarRecordTooLarge;
  }
  return kPglRetSuccess;
}

BoolErr PwcAppendBiallelicGenovecHphaseDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, PgenWriterCommon* pwcp);

HEADER_INLINE PglErr SpgwAppendBiallelicGenovecHphaseDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_main, uint32_t dosage_ct, STPgenWriter* spgwp) {
  if (unlikely(SpgwFlush(spgwp))) {
    return kPglRetWriteFail;
  }
  if (unlikely(PwcAppendBiallelicGenovecHphaseDosage16(genovec, phasepresent, phaseinfo, dosage_present, dosage_main, dosage_ct, &(spgwp->pwc)))) {
    return kPglRetVarRecordTooLarge;
  }
  return kPglRetSuccess;
}

// dosage_present cannot be null for nonzero dosage_ct
// could make dosage_main[] has length dosage_ct + dphase_ct instead of having
// separate dphase_delta[]?
BoolErr PwcAppendBiallelicGenovecDphase16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uintptr_t* __restrict dphase_present, const uint16_t* dosage_main, const int16_t* dphase_delta, uint32_t dosage_ct, uint32_t dphase_ct, PgenWriterCommon* pwcp);

HEADER_INLINE PglErr SpgwAppendBiallelicGenovecDphase16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uintptr_t* dphase_present, const uint16_t* dosage_main, const int16_t* dphase_delta, uint32_t dosage_ct, uint32_t dphase_ct, STPgenWriter* spgwp) {
  if (unlikely(SpgwFlush(spgwp))) {
    return kPglRetWriteFail;
  }
  if (unlikely(PwcAppendBiallelicGenovecDphase16(genovec, phasepresent, phaseinfo, dosage_present, dphase_present, dosage_main, dphase_delta, dosage_ct, dphase_ct, &(spgwp->pwc)))) {
    return kPglRetVarRecordTooLarge;
  }
  return kPglRetSuccess;
}

// Backfills header info, then closes the file.
PglErr SpgwFinish(STPgenWriter* spgwp);

// Last flush automatically backfills header info and closes the file.
// (caller should set mpgwp = nullptr after that)
PglErr MpgwFlush(MTPgenWriter* mpgwp);


// these close the file if open, but do not free any memory
// MpgwCleanup() handles mpgwp == nullptr, since it shouldn't be allocated on
// the stack
BoolErr SpgwCleanup(STPgenWriter* spgwp);
BoolErr MpgwCleanup(MTPgenWriter* mpgwp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PGENLIB_INTERNAL_H__
