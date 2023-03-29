#ifndef __PGENLIB_MISC_H__
#define __PGENLIB_MISC_H__

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
// See the bottom of this header file, and the pgen_spec/ subdirectory of
// plink-ng on GitHub, for details.

// Additional parameter conventions:
// - "nyparr" indicates a word-aligned, packed array of 2-bit values, while
//   "nypvec" is the vector-aligned equivalent.  "nybblearr" marks the much
//   rarer case of a packed array of 4-bit values.
// - "nypvec_01" indicates a packed, vector-aligned array of 2-bit values where
//   each value is zero or one.  This data structure was used quite a bit by
//   plink 1.9 for operating on a subset of a 2-bit-genotype array.
// - "genoarr"/"genovec" indicates a nyparr/nypvec containing genotype
//   information.
// - "interleaved_vec" is plink 2.0's preferred alternative to nypvec_01: we
//   basically stack pairs of adjacent vectors on top of each other and unpack
//   on the fly, since that tends to be faster than having to access twice as
//   much memory.

#include "plink2_bits.h"

// 10000 * major + 100 * minor + patch
// Exception to CONSTI32, since we want the preprocessor to have access to this
// value.  Named with all caps as a consequence.
#define PGENLIB_INTERNAL_VERNUM 1906

#ifdef __cplusplus
namespace plink2 {
#endif

// other configuration-ish values needed by plink2_common subset
typedef unsigned char AlleleCode;
typedef uint16_t DoubleAlleleCode;
static_assert(sizeof(DoubleAlleleCode) == 2 * sizeof(AlleleCode), "Inconsistent AlleleCode and DoubleAlleleCode definitions.");
// Set this to 65534 if AlleleCode is uint16_t, 2^24 - 1 if uint32_t.
CONSTI32(kPglMaxAltAlleleCt, 254);

CONSTI32(kPglMaxAlleleCt, kPglMaxAltAlleleCt + 1);
#define PGL_MAX_ALT_ALLELE_CT_STR "254"
#define PGL_MAX_ALLELE_CT_STR "255"
#ifdef __cplusplus
#  define kMissingAlleleCode S_CAST(plink2::AlleleCode, -1)
#  define kMissingDoubleAlleleCode S_CAST(plink2::DoubleAlleleCode, -1)
#else
#  define kMissingAlleleCode S_CAST(AlleleCode, -1)
#  define kMissingDoubleAlleleCode S_CAST(DoubleAlleleCode, -1)
#endif
CONSTI32(kAlleleCodesPerVec, kBytesPerVec / sizeof(AlleleCode));

HEADER_INLINE uintptr_t AlleleCodeCtToVecCt(uintptr_t val) {
  return DivUp(val, kAlleleCodesPerVec);
}

HEADER_INLINE uintptr_t AlleleCodeCtToAlignedWordCt(uintptr_t val) {
  return kWordsPerVec * AlleleCodeCtToVecCt(val);
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
  return NypsumWord(val);
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

// assumes subset_mask has trailing zeroes up to the next vector boundary
void FillInterleavedMaskVec(const uintptr_t* __restrict subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec);

HEADER_INLINE void CopyNyparr(const uintptr_t* __restrict source_nyparr, uint32_t nyparr_entry_ct, uintptr_t* __restrict target_nyparr) {
  memcpy(target_nyparr, source_nyparr, NypCtToWordCt(nyparr_entry_ct) * kBytesPerWord);
}

// may want bit past the end of subset_mask (i.e. position
// raw_nyparr_entry_ct) to always be allocated and unset.  This removes the
// need for some explicit end-of-bitarray checks.
void CopyNyparrNonemptySubset(const uintptr_t* __restrict raw_nyparr, const uintptr_t* __restrict subset_mask, uint32_t raw_nyparr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_nyparr);

// Copies a bit from raw_bitarr for each genoarr entry matching match_word.
// (match_word must be a multiple of kMask5555.)
void CopyGenomatchSubset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict genoarr, uintptr_t match_word, uint32_t write_bit_idx_start, uint32_t bit_ct, uintptr_t* __restrict output_bitarr);

void ExpandBytearrFromGenoarr(const void* __restrict compact_bitarr, const uintptr_t* __restrict genoarr, uintptr_t match_word, uint32_t genoword_ct, uint32_t expand_size, uint32_t read_start_bit, uintptr_t* __restrict target);


// These functions are "unsafe" since they assume trailing bits of last
// genovec/genoarr word are zeroed out.
void GenovecCount12Unsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* __restrict raw_01_ctp, uint32_t* __restrict raw_10_ctp);

void Count3FreqVec6(const VecW* geno_vvec, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp);

// vector-alignment preferred.
void GenoarrCountFreqsUnsafe(const uintptr_t* genoarr, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts);

// GenoarrCountFreqsUnsafe() wrapper that returns most common genotype,
// breaking ties in favor of the lower value.
uintptr_t MostCommonGenoUnsafe(const uintptr_t* genoarr, uint32_t sample_ct);

// geno_vvec now allowed to be unaligned.
void CountSubset3FreqVec6(const VecW* __restrict geno_vvec, const VecW* __restrict interleaved_mask_vvec, uint32_t vec_ct, uint32_t* __restrict even_ctp, uint32_t* __restrict odd_ctp, uint32_t* __restrict bothset_ctp);

// genoarr vector-alignment preferred.
void GenoarrCountSubsetFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts);

// slower GenoarrCountSubsetFreqs() which does not require
// sample_include_interleaved_vec to be precomputed
void GenoarrCountSubsetFreqs2(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts);

void GenoarrCountInvsubsetFreqs2(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_exclude, uint32_t raw_sample_ct, uint32_t sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts);

void GenoarrCountSubsetIntersectFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict subset1, const uintptr_t* __restrict subset2, uint32_t raw_sample_ct, STD_ARRAY_REF(uint32_t, 4) genocounts);

void GenovecInvertUnsafe(uint32_t sample_ct, uintptr_t* genovec);

HEADER_INLINE uintptr_t InvertGenoWordUnsafe(uintptr_t geno_word) {
  return (geno_word ^ ((~(geno_word << 1)) & kMaskAAAA));
}

// too easy to forget to multiply by 2
HEADER_INLINE void ZeroTrailingNyps(uintptr_t nyp_ct, uintptr_t* bitarr) {
  ZeroTrailingBits(nyp_ct * 2, bitarr);
}

HEADER_INLINE void SetTrailingNyps(uintptr_t nyp_ct, uintptr_t* bitarr) {
  const uintptr_t trail_ct = nyp_ct % kBitsPerWordD2;
  if (trail_ct) {
    bitarr[nyp_ct / kBitsPerWordD2] |= (~k0LU) << (nyp_ct * 2);
  }
}

// GetVint31 and Vint32Append moved to plink2_base.

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
CONSTI32(kPglNypTransposeBatch, kNypsPerCacheline);

// word width of each matrix row
CONSTI32(kPglNypTransposeWords, kWordsPerCacheline);

#ifdef __LP64__
CONSTI32(kPglNypTransposeBufbytes, (kPglNypTransposeBatch * kPglNypTransposeBatch) / 2);

void TransposeNypblock64(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1);
#else  // !__LP64__
CONSTI32(kPglNypTransposeBufbytes, (kPglNypTransposeBatch * kPglNypTransposeBatch) / 2);

void TransposeNypblock32(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1);
#endif
CONSTI32(kPglNypTransposeBufwords, kPglNypTransposeBufbytes / kBytesPerWord);

// - up to 256x256; vecaligned_buf must have size 16k (64-bit) or 32k (32-bit)
// - does NOT zero out trailing bits, because main application is ind-major-bed
//   <-> plink2 format conversion, where the zeroing would be undone...
// - important: write_iter must be allocated up to at least
//   RoundUpPow2(write_batch_size, 4) rows (may want to remove this
//   requirement)

HEADER_INLINE void TransposeNypblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecW* vecaligned_buf) {
#ifdef __LP64__
  // assert(!(write_ul_stride % 2));
  TransposeNypblock64(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, R_CAST(unsigned char*, vecaligned_buf), &(R_CAST(unsigned char*, vecaligned_buf)[kPglNypTransposeBufbytes / 2]));
#else
  TransposeNypblock32(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, R_CAST(unsigned char*, vecaligned_buf), &(R_CAST(unsigned char*, vecaligned_buf)[kPglNypTransposeBufbytes / 2]));
#endif
}


// replaces each x with (32768 - x)
// okay for dosage_main to be nullptr if dosage_ct == 0
void BiallelicDosage16Invert(uint32_t dosage_ct, uint16_t* dosage_main);

// replaces each x with -x
void BiallelicDphase16Invert(uint32_t dphase_ct, int16_t* dphase_delta);

// currently does zero trailing halfword
void GenoarrToMissingnessUnsafe(const uintptr_t* __restrict genoarr, uint32_t sample_ct, uintptr_t* __restrict missingness);

// currently does not zero trailing halfword
void GenoarrToNonmissingnessUnsafe(const uintptr_t* __restrict genoarr, uint32_t sample_ct, uintptr_t* __restrict nonmissingness);

void SparseToMissingness(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uint32_t sample_ct, uint32_t difflist_common_geno, uint32_t difflist_len, uintptr_t* __restrict missingness);

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
// more expensive to set up and consume a non-negligible fraction of L1 cache,
// so they aren't always the right choice.
// When lookup table rows are 16 bytes, they are assumed to be 16-byte aligned
// in 64-bit builds.  result[] is not assumed to be aligned.
void GenoarrLookup256x1bx4(const uintptr_t* genoarr, const void* table256x1bx4, uint32_t sample_ct, void* __restrict result);

void GenoarrLookup16x4bx2(const uintptr_t* genoarr, const void* table16x4bx2, uint32_t sample_ct, void* result);

void GenoarrLookup256x2bx4(const uintptr_t* genoarr, const void* table256x2bx4, uint32_t sample_ct, void* result);

void GenoarrLookup4x16b(const uintptr_t* genoarr, const void* table4x16b, uint32_t sample_ct, void* result);

#define PAIR_TABLE16(a, b, c, d) \
  {(a), (a), (b), (a), (c), (a), (d), (a), \
  (a), (b), (b), (b), (c), (b), (d), (b), \
  (a), (c), (b), (c), (c), (c), (d), (c), \
  (a), (d), (b), (d), (c), (d), (d), (d)}

void GenoarrLookup16x8bx2(const uintptr_t* genoarr, const void* table16x8bx2, uint32_t sample_ct, void* result);

#define QUAD_TABLE256_INTERNAL2(a, b, c, d, f2, f3, f4) \
  (a), (f2), (f3), (f4), \
  (b), (f2), (f3), (f4), \
  (c), (f2), (f3), (f4), \
  (d), (f2), (f3), (f4)
#define QUAD_TABLE256_INTERNAL3(a, b, c, d, f3, f4) \
  QUAD_TABLE256_INTERNAL2((a), (b), (c), (d), (a), (f3), (f4)), \
  QUAD_TABLE256_INTERNAL2((a), (b), (c), (d), (b), (f3), (f4)), \
  QUAD_TABLE256_INTERNAL2((a), (b), (c), (d), (c), (f3), (f4)), \
  QUAD_TABLE256_INTERNAL2((a), (b), (c), (d), (d), (f3), (f4))
#define QUAD_TABLE256_INTERNAL4(a, b, c, d, f4) \
  QUAD_TABLE256_INTERNAL3((a), (b), (c), (d), (a), (f4)), \
  QUAD_TABLE256_INTERNAL3((a), (b), (c), (d), (b), (f4)), \
  QUAD_TABLE256_INTERNAL3((a), (b), (c), (d), (c), (f4)), \
  QUAD_TABLE256_INTERNAL3((a), (b), (c), (d), (d), (f4))
#define QUAD_TABLE256(a, b, c, d) \
  {QUAD_TABLE256_INTERNAL4((a), (b), (c), (d), (a)), \
   QUAD_TABLE256_INTERNAL4((a), (b), (c), (d), (b)), \
   QUAD_TABLE256_INTERNAL4((a), (b), (c), (d), (c)), \
   QUAD_TABLE256_INTERNAL4((a), (b), (c), (d), (d))}

void GenoarrLookup256x4bx4(const uintptr_t* genoarr, const void* table256x4bx4, uint32_t sample_ct, void* result);

// Lookup table initialization functions.  table[0][0], [1][0], [2][0], and
// [3][0] must be initialized to f(0), f(1), f(2), and f(3) respectively.
void InitLookup16x4bx2(void* table16x4bx2);

void InitLookup16x8bx2(void* table16x8bx2);

void InitLookup256x1bx4(void* table256x1bx4);

void InitLookup256x2bx4(void* table256x2bx4);

void InitLookup256x4bx4(void* table256x4bx4);

void PhaseLookup4b(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const void* table56x4bx2, uint32_t sample_ct, void* result);

// [0][0]..[3][0], [17][0], and [19][0] should contain the relevant values
void InitPhaseLookup4b(void* table56x4bx2);

void PhaseLookup8b(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const void* table56x8bx2, uint32_t sample_ct, void* result);

void InitPhaseLookup8b(void* table56x8bx2);

// het-haploid prohibited.  64-entry table suffices: we use the same bits for
// phasepresent and sex_male since they can't be true simultaneously.
// phaseinfo is xor'd with bits 1 and 3 instead of 1 and 2.
void PhaseXNohhLookup4b(const uintptr_t* genoarr, const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const uintptr_t* sex_male, const void* table64x4bx2, uint32_t sample_ct, void* result);

// [0][0]..[3][0], [16][0]..[19][0]
void InitPhaseXNohhLookup4b(void* table64x4bx2);

// uses same table as PhaseXNohhLookup
void GenoarrSexLookup4b(const uintptr_t* genoarr, const uintptr_t* sex_male, const void* table64x4bx2, uint32_t sample_ct, void* result);

void InitPhaseXNohhLookup8b(void* table64x8bx2);

void GenoarrSexLookup8b(const uintptr_t* genoarr, const uintptr_t* sex_male, const void* table64x8bx2, uint32_t sample_ct, void* result);

// Unlike PhaseLookup4b(), this allows the cur_phased bit to be set when the
// genoarr entry is not 01 (het).
void VcfPhaseLookup4b(const uintptr_t* genoarr, const uintptr_t* cur_phased, const uintptr_t* phaseinfo, const void* table246x4bx2, uint32_t sample_ct, void* __restrict result);

// Precondition:
//   [0], [2], [4], [6] initialized with unphased entries
//   [32], [34], [36], [38] initialized with phased-unflipped entries
//   [162] initialized with phased-flipped case
void InitVcfPhaseLookup4b(void* table246x4bx2);

void VcfPhaseLookup2b(const uintptr_t* genoarr, const uintptr_t* cur_phased, const uintptr_t* phaseinfo, const void* table246x2bx2, uint32_t sample_ct, void* __restrict result);

void InitVcfPhaseLookup2b(void* table246x2bx2);


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

// For every missing entry in genoarr, clear the corresponding subset and
// sparse_vals entries.
void ClearGenoarrMissing1bit8Unsafe(const uintptr_t* __restrict genoarr, uint32_t* subset_sizep, uintptr_t* __restrict subset, void* __restrict sparse_vals);

void ClearGenoarrMissing1bit16Unsafe(const uintptr_t* __restrict genoarr, uint32_t* subset_sizep, uintptr_t* __restrict subset, void* __restrict sparse_vals);

double MultiallelicDiploidMinimac3R2(const uint64_t* __restrict sums, const uint64_t* __restrict hap_ssqs_x2, uint32_t nm_sample_ct, uint32_t allele_ct, uint32_t extra_phased_het_ct);

HEADER_INLINE double MultiallelicDiploidMachR2(const uint64_t* __restrict sums, const uint64_t* __restrict ssqs, uint32_t nm_sample_ct, uint32_t allele_ct) {
  return 2 * MultiallelicDiploidMinimac3R2(sums, ssqs, nm_sample_ct, allele_ct, 0);
}

// ----- end plink2_common subset -----

// other configuration-ish values
// this part of the specification is set in stone.

CONSTI32(kPglVblockSize, 65536);

// kPglDifflistGroupSize defined in plink2_base

// Currently chosen so that it plus kPglFwriteBlockSize + kCacheline - 2 is
// < 2^32, so DivUp(kPglMaxBytesPerVariant + kPglFwriteBlockSize - 1,
// kCacheline) doesn't overflow.
static const uint32_t kPglMaxBytesPerVariant = 0xfffdffc0U;
// CONSTI32(kPglMaxBytesPerDataTrack, 0x7ffff000);
// static_assert(kMaxBytesPerIO >= (int32_t)kPglMaxBytesPerDataTrack, "pgenlib assumes a single variant data track always fits in one fread/fwrite operation.");

FLAGSET_DEF_START()
  kfPgenGlobal0,
  kfPgenGlobalLdCompressionPresent = (1 << 0),
  kfPgenGlobalDifflistOrLdPresent = (1 << 1),

  // Only guaranteed to be set when present if phase or dosage also present.
  kfPgenGlobalMultiallelicHardcallFound = (1 << 2),

  kfPgenGlobalHardcallPhasePresent = (1 << 3),
  kfPgenGlobalDosagePresent = (1 << 4),
  kfPgenGlobalDosagePhasePresent = (1 << 5),
  kfPgenGlobalAllNonref = (1 << 6)
FLAGSET_DEF_END(PgenGlobalFlags);

// difflist/LD compression must not involve more than
//   raw_sample_ct / kPglMaxDifflistLenDivisor
// entries.  (however, returned difflists can have up to twice as many entries,
// when a variant is LD-compressed and the reference variant is
// difflist-compressed.)
// This value can be considered set in stone.
CONSTI32(kPglMaxDifflistLenDivisor, 8);

// Threshold for using a deltalist to represent a bitarray on disk (currently
// relevant for dosage data).  This is a tunable parameter, but must be >=
// kPglMaxDifflistLenDivisor.
CONSTI32(kPglMaxDeltalistLenDivisor, 9);

void PgrDifflistToGenovecUnsafe(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uintptr_t difflist_common_geno, uint32_t sample_ct, uint32_t difflist_len, uintptr_t* __restrict genovec);

// This covers all the possibilities.  Todo: switch all functions exporting
// multiallelic codes and/or phased dosage to use this struct.  (Biallelic
// phased hardcalls and unphased dosages are simple enough for this to be
// overkill, though.)
typedef struct PgenVariantStruct {
  MOVABLE_BUT_NONCOPYABLE(PgenVariantStruct);
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
} PgenVariant;

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

extern const uint16_t kHcToAlleleCodes[1024];

// Permits missing codes, does not remap.
void PglMultiallelicSparseToDenseMiss(const PgenVariant* pgvp, uint32_t sample_ct, AlleleCode* __restrict wide_codes);

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
//        plink2 --hard-call-threshold <...> --make-pgen
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
//         samples).  The following encodings are now defined (note that there
//         was a change of plans in Mar 2019):
//           8: 1 bit per fused vrtype-length.  Unset = vrtype 5, set = vrtype
//              0.
//           9: 2 bits, multiallelic.  0 = vrtype 5, 1 = vrtype 0, 2-3 = vrtype
//              8 with that many more bytes than vrtype 0.  Note that this is
//              limited to 16 ALT alleles.
//          10: 2 bits, phased.  0 = vrtype 5, 1 = vrtype 0, 2-3 = vrtype 16
//              with that many minus 1 bytes beyond vrtype 0.  While this is
//              also aimed at the single-sample use case, it technically
//              supports up to 15 always-phased or 7 partially-phased samples.
//          11: 4 bits, multiallelic + phased.  0 = vrtype 5, 1 = vrtype 0, 2-7
//              = vrtype 8 with that many bytes beyond vrtype 0, 9 = vrtype 16
//              phase info requiring just 1 byte, 10-15 = vrtype 24 with (x-7)
//              extra bytes required between multiallelic and phased tracks.
//          12: 2 bits, dosage, must be single-sample.  0 = vrtype 5, 1 =
//              vrtype 0, 2 = vrtype 0x45 with 2 bytes, 3 = vrtype 0x40 with 3
//              total bytes.
//          13: reserved for single-sample multiallelic + dosage.
//          14: 4 bits, phased + dosage, must be single-sample.  0 and 1 as
//              usual, 3 = vrtype 16 with 1 phaseinfo byte, 4 = vrtype 0x45
//              with 2 bytes, 5 = vrtype 0x40 with 3 total bytes, 12 = vrtype
//              0xc5 with 4 total bytes, 13 = vrtype 0xc0 with 5 total bytes,
//              15 = vrtype 0xe0 with 6 total bytes
//          15: reserved for single-sample multiallelic + phased dosage.
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
//   a. <difflist_len VINT>
//   If difflist_len is zero, that's it.  Otherwise, the difflist is organized
//   into 64-element groups (the last group will usually be smaller), to make
//   extraction of e.g. a single sample less painful.  Note that with 20k
//   samples, a difflist is space-saving even with MAF 5%:
//     ~1/400 hom alt + ~38/400 het = (~39/400) * 20k
//                                  = ~1950 sample IDs.
//     that's 31 groups, requiring about 2 + 62 + 30 + 488 + 1919 = 2501 bytes
//     (can be slightly higher since a few ID deltas may be larger than 127);
//     uncompressed storage requires 5000 bytes.
//   b. <array of group start sample IDs, each of sample_id_byte_ct>
//   c. <array of 1-byte <delta segment lengths minus 63>, with last entry
//      omitted>
//   d. <optional payload of fixed-width genotype values>
//      (in retrospect, it might have been better to position this after (e)
//      to avoid entanglement with the packed-bitarray definition, oh well...)
//   e. one "delta segment"/group: <array of <group size - 1> VINT values,
//      each indicating the difference between the current and previous sample
//      IDs; i.e. value is 1 for two consecutive samples>

// Variant record type ('vrtype') coding:
// bits 0-2:
//   000 = Simple 2-bit encoding.
//   100, 110, 111 = Simple difflist.  Low two bits store the base value.  (101
//                   isn't here since Hardy-Weinberg equilibrium prevents
//                   all het ref/alt from arising much in practice, outside of
//                   alignment/variant-calling technical artifacts that should
//                   be removed.)
//   010 = Differences-from-earlier-variant encoding ("LD compression").  The
//         last variant without this type of encoding is the base.
//         To simplify random access logic, the first variant in each vblock is
//         prohibited from using this encoding.
//   011 = Inverted differences-from-earlier-variant encoding.  (This covers
//         the case where a reference allele is 'wrong'.)  When decoding, the
//         difflist should be processed first, then the entire genovec should
//         be flipped.
//   001 = 1-bit + difflist representation.  Suppose most calls are hom ref or
//         het (e.g. a 20% MAF variant with ~4% hom alt1, ~36% het ref/alt1,
//         ~64% hom ref), then the main datatrack has just the low bits of the
//         usual 2-bit codes.  This is followed by a difflist containing the
//         hom alt1 and missing genotypes.
//         The main datatrack is preceded by a single byte indicating what
//         the two common values are: 2 low bits = <set value - unset value>,
//         next 2 bits = unset value (6 possibilities).  Top 4 bits are
//         reserved.
//   101 = All alleles are reference, no missing data.  The main datatrack is
//         empty in this case.  Although this saves only 1 byte per variant
//         over vrtype 100, this actually makes a huge difference for
//         single-sample files.
//         Since this was not defined until mid-2019, the standard plink2
//         alpha-test binaries will not use this encoding.  However,
//         alpha-2-final and later binaries interpret this encoding correctly.
//         If your workflow makes heavy use of single-sample .pgen files, you
//         can add -DFUTURE_ENCODER during compilation to unlock this feature.
//
// bit 3: multiallelic hardcalls present with alt2/alt3/... present?  If yes,
//        auxiliary data track #1 disambiguates the 0b01 (ref/altx) and 0b10
//        (altx/alty, x may equal y) hardcalls.  This contains a format byte,
//        followed by a list of ref/altx patches, then a list of altx/alty
//        patches.  All unpatched genotypes are ref/alt1 or alt1/alt1.
//        The bottom 4 bits of the format byte describe how the ref/altx patch
//        set is stored.
//   0 = Starts with a bitarray with <total ref/altx count> bits (divide by 8
//       and round up to get byte count of this component; any trailing bits in
//       the last byte must be 0), where each set bit corresponds to presence
//       of a rarealt (i.e. alt2/alt3/...).
//       ExpandBytearr(aux1_first_quarter, all_01, raw_sample_ctl, ct_01, 0,
//                     patch_01);
//       can be used to convert this into a set of sample IDs, though we may
//       want to avoid an intermediate unpacking step in practice.  Note that
//       when we're passing in sample_ct < raw_sample_ct and the main datatrack
//       is LD-compressed, we'd like ldbase_raw_genovec to be cached.
//       This is followed by a packed array of fixed-width <allele idx - 2>
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
//       65538-16777215 alts: width 24 bits.  Reasonable to prohibit more than
//                            2^24 - 1 = 16777215, since variant records are
//                            limited to 4 GiB.  I can imagine some
//                            applications of >65534 in highly variable
//                            regions, though, and it doesn't actually cost us
//                            anything to define a way to represent it.  (A
//                            plink2 binary compiled with AlleleCode typedef'd
//                            as uint32_t will run more slowly, of course, but
//                            most binaries will not be compiled that way.)
//   1 = Same as mode 0, except the initial bitarray is replaced by a difflist
//       with sample IDs.  (We could make that piece somewhat smaller by
//       storing 0-based ref/altx indexes instead, but I'm pretty sure that
//       isn't worth the performance penalty of requiring all_01 and more
//       complicated unpacking.  Though we'll need to peek at aux1[0] before
//       decompressing the main datatrack to exploit this.)
//   15 = Empty patch set.  Might remove this (storing this as mode 1 just
//        takes 1 more byte), but it may enable some relevant performance
//        optimizations.
//   2-14 are reserved for future use.  We don't define an efficient way to
//   represent a variant that e.g. has more alt2s than alt1s for now, since alt
//   alleles will usually be sorted in order of decreasing frequency, but maybe
//   this will matter in the future.
//
//   The top 4 bits describe how the altx/alty patch set is stored.  0/1/15
//   have the same meaning as they do for the ref/altx patch set; the only
//   thing that changes is the format of the packed array of values at the end.
//   2 alts: width 1.  This is treated as a special case.  Set bits correspond
//           to alt2/alt2, clear = alt1/alt2.
//   3-4 alts: width 2+2 bits.  Each stores <allele idx - 1>, with the smaller
//             number in the lower bits.  E.g. alt1/alt2 is stored as 0b0100;
//             alt3/alt3 is stored as 0b1010.
//   5-16 alts: width 4+4 bits.
//   17-256 alts: width 8+8 bits.
//   257-65536 alts: width 16+16 bits.
//   65537-16777215 alts: width 24+24 bits.
//
// bit 4: hardcall phased?  If yes, auxiliary data track #2 contains phasing
//        information for heterozygous calls.
//        The first *bit* of the track indicates whether an explicit
//        "phasepresent" bitarray is stored.  If it's set, the next het_ct bits
//        are 1-bit values, where 0 = no phasing info known, and 1 = phasing
//        info present.  If it's unset, phasing info is present for every het
//        call.
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
//        contains a 16-bit (0..2^15; 65535 missing value is only permitted in
//        unconditional-dosage case) value expressing the sum of all alt allele
//        dosages.  (Yes, making this the ref allele dosage would have been a
//        bit cleaner, but it's too late now.)
//        If the variant is multiallelic, nonzero alt2/alt3/... dosages are
//        likely to be sparse.  So,
//        - track #5 contains a delta-encoded list describing which
//          <sample_uidx x rarealt dosage> entries are nonzero, where rarealt
//          index is in the lowest bits and sample_uidx can be computed via
//          right-shift (to avoid integer-division headaches, especially
//          important since indexes in this list can be larger than 2^32).
//          We use sample_uidx here to make subsetting less painful.
//          Since each nonzero dosage value requires 16 bits, and
//          delta-encoding of a dense list adds less than 9 bits per entry,
//          there isn't much point in supporting a dense bitarray mode here.
//          To keep memory requirements sane for biobank-scale datasets when
//          the number of alt alleles is very large, each sample is prohibited
//          from having more than 255 nonzero allele dosages (this makes no
//          difference when sizeof(AlleleCode) == 1, but it may matter later).
//        - track #6 contains the rarealt nonzero dosage values.
//        Note that this and the other dosage modes are in ADDITION to
//        hardcalls.  This increases filesize by up to 12.5%, but makes the
//        reader substantially simpler; --hard-call-threshold logic is nicely
//        compartmentalized.
//   10 = unconditional dosage (just track #4).
//   11 = dosage bitarray.  In this case, auxiliary data track #3 contains an
//        array of 1-bit values indicating which samples have dosages.  If the
//        variant is multiallelic, tracks #5 and 6 are as described above.
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
//        can be set without the other.  (However, in practice, bit 4 is almost
//        always set when bit 7 is, since that enables more efficient storage
//        of 0|0.99, 1|0.02, and similar pairs.)
//        When phased dosages are present, track #8 contains values
//        representing <(hap1 alt prob) - (hap2 alt prob)>, etc., where the
//        underlying values are represented in [0..16384] (so the signed
//        difference is in [-16384..16384]).  Track #4 contains the
//        corresponding sums; parity does NOT need to match (necessary to allow
//        this since we treat omitted values as zero; and since we are allowing
//        it, no point in making parity match in other situations either).  In
//        fixed-width case, -32768 should be stored in track #8 when the entire
//        call is missing, while 0 and missing-phase are considered synonymous.
//        In the biallelic case, if a hardcall is phased, a dosage is present,
//        and no explicit dosage-phase is, we define it to mean the unique
//        dphase_delta sequence with maximal absolute value, and --make-pgen
//        takes advantage of it.  This definition technically works for
//        triallelic variants as well, but it breaks down with 4 alleles, so we
//        prohibit hardcall-phase + dosage + no-dosage-phase with more than 2
//        alleles.
//        In the multiallelic case, tracks #9 and #10 are analogous to #5 and
//        #6.
//
// Representation of variable ploidy (MT) was considered, but rejected since
// dosages should be at least as appropriate for MT.
// Oxford/VCF-style storage of separate probabilities for every possible
// genotype (e.g. P(AA), P(AB), P(BB) instead of just 2P(AA) + P(AB) and
// 2P(BB) + P(AB)) is tentatively rejected due to (i) lack of relevance to
// PLINK's analysis functions and (ii) high storage cost where we can afford it
// least.  In principle, this is subject to reevaluation if (i) changes, but
// given the poor interaction with phased dosages, it's probably better to just
// think of them as permanently outside PLINK's scope.

// maximum prime < 2^32 is 4294967291; quadratic hashing guarantee breaks down
// past that divided by 2.
CONSTI32(kPglMaxVariantCt, 0x7ffffffd);

CONSTI32(kPglMaxSampleCt, 0x7ffffffe);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PGENLIB_MISC_H__
