#ifndef __PGENLIB_INTERNAL_H__
#define __PGENLIB_INTERNAL_H__

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
// - Individual variant records are prohibited from being >= 4GB, to reduce
//   integer overflow issues.  (This may be reduced to 2GB later, but I'll
//   attempt to handle the 2-4GB range properly for now since it's conceivable
//   for multiallelic records in very large datasets to reach that size.)
// - (later todo: include stuff like file creation command in .bim successor
//   header; that doesn't really belong in a binary file.)

// Additional parameter conventions:
// - "quaterarr" indicates a word-aligned, packed array of 2-bit values, while
//   "quatervec" is the vector-aligned equivalent.  Similarly, "hexadecarr"
//   marks the much rarer case of a packed array of 4-bit values, etc.
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
// Exception to CONSTU31, since we want the preprocessor to have access to this
// value.  Named with all caps as a consequence.
#define PGENLIB_INTERNAL_VERNUM 800

// other configuration-ish values needed by plink2_common subset
typedef unsigned char AltAlleleCt;
// don't use CONSTU31 for this since it may need the 32nd bit in the future
#define kPglMaxAltAlleleCt S_CAST(uint32_t, S_CAST(AltAlleleCt, -2))

#ifdef __cplusplus
namespace plink2 {
#endif

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

HEADER_INLINE uintptr_t GetHexadecarrEntry(const uintptr_t* hexadecarr, uint32_t idx) {
  return (hexadecarr[idx / kBitsPerWordD4] >> (4 * (idx % kBitsPerWordD4))) & 15;
}

#ifdef USE_SSE42
HEADER_INLINE uint32_t Popcount01Word(uintptr_t val) {
  return PopcountWord(val);
}
#else
HEADER_INLINE uint32_t Popcount01Word(uintptr_t val) {
  return QuatersumWord(val);
}
#endif

// safe errstr_buf size for pgen_init_phase{1,2}()
CONSTU31(kPglErrstrBufBlen, kPglFnamesize + 256);

// assumes subset_mask has trailing zeroes up to the next vector boundary
void FillInterleavedMaskVec(const uintptr_t* __restrict subset_mask, uint32_t base_vec_ct, uintptr_t* interleaved_mask_vec);

HEADER_INLINE void CopyQuaterarr(const uintptr_t* __restrict source_quaterarr, uint32_t quaterarr_entry_ct, uintptr_t* __restrict target_quaterarr) {
  memcpy(target_quaterarr, source_quaterarr, QuaterCtToWordCt(quaterarr_entry_ct) * kBytesPerWord);
}

// may want bit past the end of subset_mask (i.e. position
// raw_quaterarr_entry_ct) to always be allocated and unset.  This removes the
// need for some explicit end-of-bitarray checks.
void CopyQuaterarrNonemptySubset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_entry_ct, uint32_t subset_entry_ct, uintptr_t* __restrict output_quaterarr);

void GenovecCountFreqsUnsafe(const uintptr_t* genovec, uint32_t sample_ct, uint32_t* counts);

void GenovecCountSubsetFreqs(const uintptr_t* __restrict genovec, const uintptr_t* __restrict sample_include_interleaved_vec, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t* genocounts);

// slower GenovecCountSubsetFreqs() which does not require
// sample_include_interleaved_vec to be precomputed (and incidentally doesn't
// require vector alignment)
void GenoarrCountSubsetFreqs2(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t* genocounts);

void GenoarrCountSubsetIntersectFreqs(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict subset1, const uintptr_t* __restrict subset2, uint32_t raw_sample_ct, uint32_t* genocounts);

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

// Returns 0x80000000U instead of UINT32_MAX so overflow check works properly
// in 32-bit build.  Named "GetVint31" to make it more obvious that a 2^31
// return value can't be legitimate.
HEADER_INLINE uint32_t GetVint31(const unsigned char* buf_end, const unsigned char** buf_iterp) {
  if (buf_end > (*buf_iterp)) {
    uint32_t vint32 = *((*buf_iterp)++);
    if (vint32 <= 127) {
      return vint32;
    }
    vint32 &= 127;
    uint32_t shift = 7;
    while (buf_end > (*buf_iterp)) {
      uint32_t uii = *((*buf_iterp)++);
      vint32 |= (uii & 127) << shift;
      if (uii <= 127) {
        return vint32;
      }
      shift += 7;
      // currently don't check for shift >= 32 (that's what validate_vint31()
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
HEADER_INLINE uint32_t GetVint31Unsafe(const unsigned char** buf_iterp) {
  uint32_t vint32 = *(*buf_iterp)++;
  if (vint32 <= 127) {
    return vint32;
  }
  vint32 &= 127;
  for (uint32_t shift = 7; shift < 32; shift += 7) {
    uint32_t uii = *(*buf_iterp)++;
    vint32 |= (uii & 127) << shift;
    if (uii <= 127) {
      return vint32;
    }
  }
  return 0x80000000U;
}

// Does not update buf_iter.
HEADER_INLINE uint32_t PeekVint31(const unsigned char* buf_iter, const unsigned char* buf_end) {
  if (buf_end > buf_iter) {
    uint32_t vint32 = *buf_iter++;
    if (vint32 <= 127) {
      return vint32;
    }
    vint32 &= 127;
    uint32_t shift = 7;
    while (buf_end > buf_iter) {
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
  for (uint32_t shift = 7; shift < 32; shift += 7) {
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

// main batch size
CONSTU31(kPglQuaterTransposeBatch, kQuatersPerCacheline);

// word width of each matrix row
CONSTU31(kPglQuaterTransposeWords, kWordsPerCacheline);

#ifdef USE_AVX2
CONSTU31(kPglQuaterTransposeBufbytes, kPglQuaterTransposeBatch);

void TransposeQuaterblockAvx2(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, void* __restrict buf0);
#else
CONSTU31(kPglQuaterTransposeBufbytes, (kPglQuaterTransposeBatch * kPglQuaterTransposeBatch) / 2);

void TransposeQuaterblockNoAvx2(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* __restrict write_iter, unsigned char* __restrict buf0, unsigned char* __restrict buf1);
#endif
CONSTU31(kPglQuaterTransposeBufwords, kPglQuaterTransposeBufbytes / kBytesPerWord);
// up to 256x256; vecaligned_buf must have size 32k (256 bytes in AVX2 case)
// write_iter must be allocated up to at least
//   RoundUpPow2(write_batch_size, 2) rows

HEADER_INLINE void TransposeQuaterblock(const uintptr_t* read_iter, uint32_t read_ul_stride, uint32_t write_ul_stride, uint32_t read_batch_size, uint32_t write_batch_size, uintptr_t* write_iter, VecUL* vecaligned_buf) {
#ifdef USE_AVX2
  // assert(!(write_ul_stride % 2));
  TransposeQuaterblockAvx2(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, vecaligned_buf);
#else
  TransposeQuaterblockNoAvx2(read_iter, read_ul_stride, write_ul_stride, read_batch_size, write_batch_size, write_iter, R_CAST(unsigned char*, vecaligned_buf), &(R_CAST(unsigned char*, vecaligned_buf)[kPglQuaterTransposeBufbytes / 2]));
#endif
}


// replaces each x with (32768 - x)
// okay for dosage_vals to be nullptr if dosage_ct == 0
void BiallelicDosage16Invert(uint32_t dosage_ct, uint16_t* dosage_vals);

// currently does zero trailing halfword
void GenovecToMissingnessUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict missingness);

// currently does not zero trailing halfword
void GenovecToNonmissingnessUnsafe(const uintptr_t* __restrict genovec, uint32_t sample_ct, uintptr_t* __restrict nonmissingness);

// N.B. assumes trailing bits of loadbuf have been filled with 1s, not 0s
// Also takes genoarr word count instead of sample count.
void SplitHomRef2hetUnsafeW(const uintptr_t* genoarr, uint32_t inword_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf);


void SplitHomRef2het(const uintptr_t* genoarr, uint32_t sample_ct, uintptr_t* __restrict hom_buf, uintptr_t* __restrict ref2het_buf);

// ----- end plink2_common subset -----

// other configuration-ish values
// this part of the specification is set in stone.

CONSTU31(kPglVblockSize, 65536);

// Currently chosen so that it plus kPglFwriteBlockSize + kCacheline - 2 is
// < 2^32, so DivUp(kPglMaxBytesPerVariant + kPglFwriteBlockSize - 1,
// kCacheline) doesn't overflow.
static const uint32_t kPglMaxBytesPerVariant = 0xfffdffc0U;
// CONSTU31(kPglMaxBytesPerDataTrack, 0x7ffff000);
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
CONSTU31(kPglDifflistGroupSize, 64);

FLAGSET_DEF_START()
  kfPgenGlobal0,
  kfPgenGlobalLdCompressionPresent = (1 << 0),
  kfPgenGlobalDifflistOrLdPresent = (1 << 1),
  kfPgenGlobalHardcallPhasePresent = (1 << 2),
  kfPgenGlobalDosagePresent = (1 << 3),
  kfPgenGlobalDosagePhasePresent = (1 << 4),
  kfPgenGlobalAllNonref = (1 << 5)
FLAGSET_DEF_END(PgenGlobalFlags);

FLAGSET_DEF_START()
  kfPgrLdcache0,
  kfPgrLdcacheQuater = (1 << 0),
  kfPgrLdcacheDifflist = (1 << 1),
  kfPgrLdcacheAllHets = (1 << 2),
  kfPgrLdcacheRefalt1Genocounts = (1 << 3)
FLAGSET_DEF_END(PgrLdcacheFlags);

// difflist/LD compression should never involve more than
//   raw_sample_ct / kPglMaxDifflistLenDivisor
// entries.  (however, returned difflists can have up to twice as many entries,
// when a variant is LD-compressed and the reference variant is
// difflist-compressed.)
CONSTU31(kPglMaxDifflistLenDivisor, 8);

// threshold for using a deltalist to represent a bitarray on disk (currently
// relevant for dosage data)
CONSTU31(kPglMaxDeltalistLenDivisor, 9);

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
//      larger values, and 0x05..0x0f, reserved for now.
//
// 3. If not plink1-format,
//    a. 4-byte # of variants; call this M.
//    b. 4-byte # of samples, call this N.
//    c. Additional 1-byte header "control" value (PgenHeaderCtrl).  May be
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
//               the entry is zero, and 8 (multiallelic) if the record has 1-3
//               extra bytes.  Designed for single-sample files sharing a
//               single .bim-like file (note that if they don't share a .bim,
//               .bim size will dominate), but it's usable whenever there's no
//               variant where >2 samples have a rare alternate allele
//               (assuming <16 alt alleles).
//         1001: No difflist/LD/onebit compression, 4 bit
//               (vrec_len - ceil(sample_ct / 4)) value.  vrtype is zero if the
//               entry is zero, and 8 if the record has 1-15 extra bytes.
//       bits 4-5: alt allele count storage (00 = unstored, 01-11 = bytes per
//                 count)
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
//            bytes, or 4 bits).
//       iii. if bits 4-5 of {3c} aren't 00, array of alt allele counts.
//        iv. nonref flags info, if explicitly stored
//      (this representation allows more efficient random access)
//    If mode 0x02-0x04, and nonref flags info explicitly stored, just that
//    bitarray.
//
// 5. The variant records.  See below for details.

// Difflist format (used for onebit, sparse variant, and LD compression):
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
//   d. [array of 2-bit replacement genotype values]
//   e. one "delta segment"/group: [array of [group size - 1] VINT values,
//      each indicating the difference between the current and previous sample
//      IDs; i.e. value is 1 for two consecutive samples]
//   f. [if multiallelic, array of appropriate-length disambiguation values]


// PgenFileInfo and PgenReader are the main exported "classes".
// Exported functions involving these data structure should all have
// "pgfi"/"pgr" in their names.

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
  //         for multiallelic variants, if the base value is 0b11 (missing),
  //         auxiliary data track #1 only contains entries for explicitly
  //         listed 0b11 values, the rest are assumed to be actual missing
  //         data.  (101 should practically never happen--gross violation of
  //         Hardy-Weinberg equilibrium--so it's reserved for future use.)
  //   010 = Differences-from-earlier-variant encoding ("LD compression").  The
  //         last variant without this type of encoding is the base.
  //         To simplify random access logic, the first variant in each vblock
  //         is prohibited from using this encoding.
  //   011 = Inverted differences-from-earlier-variant encoding.  (This covers
  //         the case where a reference allele is "wrong".)  When decoding, the
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
  //         reserved.  When the set value is 3, it does NOT represent
  //         rarealts; those must be explicitly spelled out in the difflist.
  // bit 3: more than 1 alt allele?
  // bit 4: hardcall phased?  if yes, auxiliary data track #2 contains
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
  //   01 = dosage list.  auxiliary data track #3 contains a delta-encoded list
  //        of sample IDs (like a difflist, but with no genotypes).  track #4
  //        contains a 16-bit (0..2^15; 65535 missing value is only permitted
  //        in unconditional-dosage case) value for each allele except the ref
  //        (intended to omit the last alt instead, but forgot, and it's too
  //        late now, oops).  When the variant is multiallelic, order is
  //        [sample0 alt1 prob] [sample0 alt2 prob] ... [sample1 alt1 prob] ...
  //        Note that this and the other dosage modes are in ADDITION to
  //        hardcalls.  This increases filesize by up to 12.5%, but makes the
  //        reader substantially simpler; --hard-call-threshold logic is nicely
  //        compartmentalized.
  //   10 = unconditional dosage (just track #4).
  //   11 = dosage bitarray.  in this case, auxiliary data track #3 contains an
  //        array of 1-bit values indicating which samples have dosages.
  //   bgen 1.2 format no longer permits fractional missingness, so no good
  //   reason for us to support it.
  //   considered putting *all* dosage data at the end of the file (like I will
  //   do for phase set info); this could actually be worthwhile for
  //   unconditional dosages, but it doesn't work well when only some samples
  //   have dosage data.
  // bit 7: some dosages phased?  if yes, and dosages are not unconditionally
  //        present, auxiliary data track #5 is either a single zero byte
  //        (indicating that all dosages are phased), or a bitarray of length
  //        (dosage_ct + 1) where the first bit is set, and the other bits
  //        indicate whether phase info is present for that sample (unset = no
  //        phasing info)
  //        note that this is independent of bit 4; either can be set without
  //        the other.
  //        when phased dosages are present, track #6 contains values
  //        representing [(hap1 alt1 prob) - (hap2 alt1 prob)],
  //        [(hap1 alt2 prob) - (hap2 alt2 prob)], etc., where the underlying
  //        values are represented in [0..16384] (so the signed difference is
  //        in [-16384..16384]).  track #4 contains the corresponding sums;
  //        parity should always match.  In fixed-width case, -32768 should be
  //        stored in track #6 when the entire call is missing; may use this to
  //        also represent missing-phase later, but postpone that decision.
  //
  // Representation of variable ploidy (MT) was considered, but rejected since
  // dosages should be at least as appropriate for MT.
  // Oxford/VCF-style storage of separate probabilities for every possible
  // genotype (e.g. P(AA), P(AB), P(BB) instead of just 2P(AA) + P(AB) and
  // 2P(BB) + P(AB)) is tentatively rejected due to (i) lack of relevance to
  // PLINK's analysis functions and (ii) high storage cost where we can afford
  // it least.  However, this is subject to reevaluation if (i) changes.
  //
  //
  // base pointer is null if mode is 0x01-0x04 (const_vrtype != UINT32_MAX).
  // if not nullptr, required to be length >=
  //   max(raw_variant_ct + 1, RoundUpPow2(raw_variant_ct, kBytesPerWord))
  unsigned char* vrtypes;

  // alt allele counts.  if >1, auxiliary data track #1 disambiguates all the
  // "missing or rare alt" explicit hardcalls.  genotype representation is:
  //   low bits: smaller [1-based alt allele idx], 0 = ref
  //   high bits: larger [1-based alt allele idx]
  //   ...
  //   2 alts: 1-bit array with 0 = missing, 1 = nonmissing.  Then, for the
  //     nonmissing subset, 2 low bits.  The high bits entry is omitted because
  //     the value has to be alt2; optimize the common case!
  //   3 alts: 1-bit nonmissingness array,  2 low bits, 2 high bits
  //   4-15 alts: 1-bit nonmissingness array, then 4 low bits, 4 high bits
  //   16-255 alts: 1-bit nonmissingness array; then 8 low bits, 8 high bits
  //   (the following is also defined, but not implemented for now:
  //   256-4095 alts: 1-bit nonmissingness array; 12 low bits, 12 high bits
  //   4096-65535 alts: 1-bit nonmissingness array; 16 low bits, 16 high bits
  //   65536-16777215 alts: 1-bit nonmissingness array; 24 low bits, 24 high
  //   bits; the latter might be necessary in the most variable regions if we
  //   use tiles...)
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

  uint32_t max_alt_allele_ct;
  uint32_t max_dosage_alt_allele_ct;

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
  uintptr_t* ldbase_genovec;
  uintptr_t* ldbase_raregeno;

  // when ldbase_difflist_ids[] is initialized, element [ldbase_difflist_len]
  // must be set to sample_ct.
  uint32_t* ldbase_difflist_sample_ids;

  uintptr_t* ldbase_all_hets;

  // common genotype can be looked up from vrtypes[]

  uint32_t ldbase_refalt1_genocounts[4];

  uintptr_t* workspace_vec;  // must hold raw_sample_ct entries

  // currently must hold (raw_sample_ct / kPglMaxDifflistLenDivisor)
  // entries; may need to double the sizes later
  // some top-level interface functions use these, so several lower-level
  // functions cannot
  uintptr_t* workspace_raregeno_vec;
  uint32_t* workspace_difflist_sample_ids;

  // must hold (raw_sample_ct / kPglMaxDifflistLenDivisor) entries
  uintptr_t* workspace_raregeno_tmp_loadbuf;
  uint32_t* workspace_difflist_sample_ids_tmp;

  uintptr_t* workspace_aux1_nonmissing_vec;
  uintptr_t* workspace_aux1_code_vec;

  uintptr_t* workspace_all_hets;

  uint32_t* workspace_ambig_sample_ids;
  uint32_t workspace_ambig_id_ct;

  uintptr_t* workspace_dosage_present;
  uintptr_t* workspace_dosage_phased;

  // phase set loading (mode 0x11) unimplemented for now; should be a sequence
  // of (sample ID, [uint32_t phase set begin, set end), [set begin, set end),
  // ...).
};

typedef struct PgenReaderStruct PgenReader;

// might want this value to be typed...
CONSTU31(kPglVrtypePlink1, 256);

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

HEADER_INLINE uint32_t VrtypeMultiallelic(uint32_t vrtype) {
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

HEADER_INLINE uintptr_t GetAux1AlleleByteCt(uint32_t alt_allele_ct, uint32_t aux1_nonmissing_ct) {
  assert(alt_allele_ct >= 2);
  if (alt_allele_ct == 2) {
    return DivUp(aux1_nonmissing_ct, 4);
  }
  if (alt_allele_ct == 3) {
    return DivUp(aux1_nonmissing_ct, 2);
  }
  // one byte per entry for alt_allele_ct < 16, two bytes for 16..255
  return ((alt_allele_ct >= 16) + 1) * aux1_nonmissing_ct;
  // todo: alt_allele_ct > 255
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
//   if appropriate, and initialize pgfi.gflags, pgfi.max_alt_allele_ct, and
//   pgfi.max_dosage_alt_allele_ct.
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

PglErr PgrInit(const char* fname, uint32_t max_vrec_width, PgenFileInfo* pgfip, PgenReader* pgrp, unsigned char* pgr_alloc);

// practically all these functions require genovec to be allocated up to
// vector, not word, boundary
void PgrPlink1ToPlink2InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec);

void PgrPlink2ToPlink1InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec);

void PgrDifflistToGenovecUnsafe(const uintptr_t* __restrict raregeno, const uint32_t* difflist_sample_ids, uintptr_t difflist_common_geno, uint32_t sample_ct, uint32_t difflist_len, uintptr_t* __restrict genovec);

// This will normally extract only the genotype indexes corresponding to set
// bits in sample_include.  Set sample_ct == raw_sample_ct if you don't want
// any subsetting to occur (in this case sample_include is ignored, can be
// nullptr).
// Only the maintrack is loaded.  00 = hom ref, 01 = het ref/alt1,
// 10 = hom alt1, 11 = missing or anything else.
// If multiallelic_relevant is set, and the current variant is multiallelic,
// pgr.workspace_ambig_sample_ids and pgr.workspace_ambig_id_ct are updated.
// "unsafe": sample_ct cannot be zero.  Trailing bits of genovec are not zeroed
// out.
// Ok if genovec only has space for sample_ct values.
PglErr PgrReadRefalt1GenovecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec);

// Loads the specified variant as a difflist if that's more efficient, setting
// difflist_common_geno to the common genotype value in that case.  Otherwise,
// genovec is populated and difflist_common_geno is set to UINT32_MAX.
//
// Note that the returned difflist_len can be much larger than
// max_simple_difflist_len when the variant is LD-encoded; it's bounded by
//   2 * (raw_sample_ct / kPglMaxDifflistLenDivisor).
PglErr PgrReadRefalt1DifflistOrGenovecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t max_simple_difflist_len, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uint32_t* difflist_common_geno_ptr, uintptr_t* __restrict main_raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr);

// This is necessary when changing sample_include, unless the new query is
// iterating from the first variant.  (Which can almost never be assumed in
// plink2 since variant_include[] may not include the first variant.)
HEADER_INLINE void PgrClearLdCache(PgenReader* pgrp) {
  pgrp->ldbase_stypes &= kfPgrLdcacheAllHets;

  // bugfix, ld_load_necessary() was otherwise claiming that reload wasn't
  // necessary in certain cases
  pgrp->ldbase_vidx = 0x80000000U;
}

PglErr PgrGetRefalt1GenotypeCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uint32_t* genocounts);

// allele_idx is set to 0 for ref, 1 for alt1, 2 for alt2, etc.
// frequencies are computed on the fly.  ties are broken in favor of the
// lower-indexed allele.
// possible todo: also provide ..._common2_then_subset() function.
// better default than the functions above for machine learning/GWAS, since the
// reference allele is "wrong" sometimes.
PglErr PgrReadGenovecSubsetThenCommon2(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uint32_t* __restrict maj_allele_idx_ptr, uint32_t* __restrict second_allele_idx_ptr, uint32_t* __restrict allele_ct_buf);

// Loads a quatervec with counts of a single allele (allele_idx 0 corresponds
// to the reference allele, allele_idx 1 corresponds to alt1, etc.).  0b11 ==
// missing call.
// Note that calling this with allele_idx == 0 is similar to a plink1 load
// (except with missing == 0b11, of course).
// todo: provide a difflist interface once anyone wants it.
PglErr PgrReadAlleleCountvecSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgrp, uintptr_t* __restrict allele_countvec);

// todo: add functions which directly support MAF-based queries.  Note that
// when the difflist representation is used, we can disqualify some low-MAF
// variants without actually loading the genotype data, since the size of the
// record puts an upper bound on the alt allele frequency.

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

// PglErr PgrReadRefalt1GenovecHphaseRawUnsafe(uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phaseraw, uint32_t* phasepresent_ct_ptr);

PglErr PgrReadRefalt1GenovecHphaseSubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* phasepresent_ct_ptr);

// if dosage_present and dosage_vals are nullptr, dosage data is ignored
PglErr PgrReadRefalt1GenovecDosage16SubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr, uint32_t* is_explicit_alt1_ptr);

PglErr PgrGetRefNonrefGenotypeCountsAndDosage16s(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, double* mach_r2_ptr, uint32_t* genocounts, uint64_t* all_dosages);

// ok for both dosage_present and dosage_vals to be nullptr when no dosage data
// is present
PglErr PgrReadRefalt1GenovecHphaseDosage16SubsetUnsafe(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* phasepresent_ct_ptr, uintptr_t* __restrict dosage_present, uint16_t* dosage_vals, uint32_t* dosage_ct_ptr, uint32_t* is_explicit_alt1_ptr);

// interface used by --make-pgen, just performs basic LD/difflist decompression
// (still needs multiallelic and dosage-phase extensions)
PglErr PgrReadRaw(uint32_t vidx, PgenGlobalFlags read_gflags, PgenReader* pgrp, uintptr_t** loadbuf_iter_ptr, unsigned char* loaded_vrtype_ptr);

PglErr PgrValidate(PgenReader* pgrp, char* errstr_buf);

// missingness bit is set iff hardcall is not present (even if dosage info *is*
// present)
PglErr PgrReadMissingness(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict missingness, uintptr_t* __restrict genovec_buf);

// either missingness_hc (hardcall) or missingness_dosage must be non-null
// missingness_dosage must be vector-aligned
PglErr PgrReadMissingnessMulti(const uintptr_t* __restrict sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReader* pgrp, uintptr_t* __restrict missingness_hc, uintptr_t* __restrict missingness_dosage, uintptr_t* __restrict hets, uintptr_t* __restrict genovec_buf);


// failure = kPglRetReadFail
BoolErr PgfiCleanup(PgenFileInfo* pgfip);

BoolErr PgrCleanup(PgenReader* pgrp);


struct PgenWriterCommonStruct {
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

  // you can get a ~1-2% speedup by writing directly to genovec and swapping
  // it with ldbase_genovec when appropriate; don't think that's worth
  // supporting, given the messier API.
  // uintptr_t* genovec;

  uint32_t ldbase_genocounts[4];

  // should match ftello() return value in singlethreaded case, but be set to
  // zero in multithreaded case
  uint64_t vblock_fpos_offset;

  // these must hold sample_ct entries (could be fewer if not subsetting, but
  // let's play it safe)
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

CONSTU31(kPglFwriteBlockSize, 131072);

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
  struct PgenWriterCommonStruct pwc;
  FILE* pgen_outfile;
};

struct MTPgenWriterStruct {
  FILE* pgen_outfile;
  uint32_t thread_ct;
  struct PgenWriterCommonStruct* pwcs[];
};

typedef struct STPgenWriterStruct STPgenWriter;
typedef struct MTPgenWriterStruct MTPgenWriter;

void PreinitSpgw(STPgenWriter* spgwp);

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

PglErr SpgwAppendBiallelicGenovec(const uintptr_t* __restrict genovec, STPgenWriter* spgwp);

// trailing bits of raregeno must be zeroed out
// all raregeno entries assumed to be unequal to difflist_common_geno; the
// difflist should be compacted first if this isn't true (might be possible
// with multiallelic projections?)
// difflist_len must be <= 2 * (sample_ct / kPglMaxDifflistLenDivisor);
// there's an assert checking this
void PwcAppendBiallelicDifflistLimited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, PgenWriterCommon* pwcp);

PglErr SpgwAppendBiallelicDifflistLimited(const uintptr_t* __restrict raregeno, const uint32_t* __restrict difflist_sample_ids, uint32_t difflist_common_geno, uint32_t difflist_len, STPgenWriter* spgwp);

// trailing bits of refalt1_genovec must be zeroed out
// not implemented yet
PglErr SpgwAppendMultiallelicCounts(const uintptr_t** __restrict alt_countvecs);

// phasepresent == nullptr ok, that indicates that ALL heterozygous calls are
// phased.  Caller should use e.g. PwcAppendBiallelicGenovec() if it's known
// in advance that no calls are phased.
// Ok for phaseinfo to have bits set at non-het calls, NOT currently okay for
//   phasepresent
// void pwc_append_biallelic_genovec_hphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, PgenWriterCommon* pwcp);

// phasepresent == nullptr ok
// ok for trailing bits of phaseinfo to not be zeroed out, NOT currently ok for
//   phasepresent
PglErr SpgwAppendBiallelicGenovecHphase(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, STPgenWriter* spgwp);

// dosage_vals[] has length dosage_ct, not sample_ct
void PwcAppendBiallelicGenovecDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, PgenWriterCommon* pwcp);

PglErr SpgwAppendBiallelicGenovecDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, STPgenWriter* spgwp);

void PwcAppendBiallelicGenovecHphaseDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, PgenWriterCommon* pwcp);

PglErr SpgwAppendBiallelicGenovecHphaseDosage16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uint16_t* dosage_vals, uint32_t dosage_ct, STPgenWriter* spgwp);

// dphase_present can be nullptr if dosage_ct == dphase_ct
// dosage_present cannot be null for nonzero dosage_ct
// dosage_vals[] has length dosage_ct + dphase_ct
// PglErr spgw_append_biallelic_genovec_dphase16(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uintptr_t* __restrict dosage_present, const uintptr_t* dphase_present, const uint16_t* dosage_vals, uint32_t dosage_ct, uint32_t dphase_ct, STPgenWriter* spgwp);


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
