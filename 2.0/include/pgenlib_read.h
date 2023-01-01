#ifndef __PGENLIB_READ_H__
#define __PGENLIB_READ_H__

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


// pgenlib_read contains reader-specific code.

#include "pgenlib_misc.h"

#ifdef __cplusplus
namespace plink2 {
#endif

FLAGSET_DEF_START()
  kfPgrLdcache0,
  kfPgrLdcacheNyp = (1 << 0),
  kfPgrLdcacheDifflist = (1 << 1),
  kfPgrLdcacheRawNyp = (1 << 2),
  // may also want RawDifflist
  kfPgrLdcacheBasicGenocounts = (1 << 3)
FLAGSET_DEF_END(PgrLdcacheFlags);

// PgenFileInfo and PgenReader are the main exported "classes".
// Exported functions involving these data structure should all have
// "pgfi"/"pgr" in their names.

// Note that this can be default-copied.
typedef struct PgenFileInfoStruct {
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

  // Variant record type codes.
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

  // if using per-variant fread(), this is non-null during PgenFileInfo
  // initialization, but it's then "moved" to the first Pgen_reader and set to
  // nullptr.
  FILE* shared_ff;

  // can only be non-null after PgfiInitPhase1 and before PgfiInitPhase2, and
  // only if the external-index-file representation is used.
  FILE* pgi_ff;

  const unsigned char* block_base;  // nullptr if using per-variant fread()
  uint64_t block_offset;
} PgenFileInfo;

typedef struct PgenReaderMainStruct {
  MOVABLE_BUT_NONCOPYABLE(PgenReaderMainStruct);
  // would like to make this const, but that makes initialization really
  // annoying in C99
  struct PgenFileInfoStruct fi;

  // ----- Mutable state -----
  // If we don't fseek, what's the next variant we'd read?
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

  // when ldbase_difflist_sample_ids[] is initialized, element
  // [ldbase_difflist_len] must be set to sample_ct.
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
  uint64_t* workspace_imp_r2;  // needed in multiallelic case

  uintptr_t* workspace_all_hets;
  uintptr_t* workspace_subset;  // currently used for hphase decoding

  uintptr_t* workspace_dosage_present;
  uintptr_t* workspace_dphase_present;

  // phase set loading (mode 0x11) unimplemented for now; should be a sequence
  // of (sample ID, [uint32_t phase set begin, set end), [set begin, set end),
  // ...).
} PgenReaderMain;

typedef struct PgenReaderStruct {
#ifdef __cplusplus
  PgenReaderMain& GET_PRIVATE_m() { return m; }
  PgenReaderMain const& GET_PRIVATE_m() const { return m; }
 private:
#endif
  PgenReaderMain m;
} PgenReader;

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
  return (vrtype & 4) && ((vrtype & 3) != 1);
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

HEADER_INLINE unsigned char* PgrGetFreadBuf(PgenReader* pgr_ptr) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  return pgrp->fread_buf;
}

HEADER_INLINE unsigned char* PgrGetVrtypes(PgenReader* pgr_ptr) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  return pgrp->fi.vrtypes;
}

HEADER_INLINE uint32_t PgrGetVrtype(const PgenReader* pgr_ptr, uint32_t vidx) {
  const PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  if (pgrp->fi.vrtypes) {
    return pgrp->fi.vrtypes[vidx];
  }
  return pgrp->fi.const_vrtype;
}

HEADER_INLINE uintptr_t* PgrGetNonrefFlags(PgenReader* pgr_ptr) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  return pgrp->fi.nonref_flags;
}

HEADER_INLINE PgenGlobalFlags PgrGetGflags(const PgenReader* pgr_ptr) {
  const PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  return pgrp->fi.gflags;
}

HEADER_INLINE uint32_t PgrGetMaxAlleleCt(const PgenReader* pgr_ptr) {
  const PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  return pgrp->fi.max_allele_ct;
}

HEADER_INLINE void PgrSetFreadBuf(unsigned char* fread_buf, PgenReader* pgr_ptr) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  pgrp->fread_buf = fread_buf;
}

HEADER_INLINE void PgrCopyBaseAndOffset(const PgenFileInfo* pgfip, uint32_t thread_ct, PgenReader** pgr_ptr_arr) {
  for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
    PgenReaderMain* pgrp = &GET_PRIVATE(*(pgr_ptr_arr[tidx]), m);
    pgrp->fi.block_base = pgfip->block_base;
    pgrp->fi.block_offset = pgfip->block_offset;
  }
}

// This is necessary when changing sample_include, unless the new query is
// iterating from the first variant.  (Which can almost never be assumed in
// plink2 since variant_include[] may not include the first variant.)
HEADER_INLINE void PgrClearLdCache(PgenReader* pgr_ptr) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  pgrp->ldbase_stypes &= kfPgrLdcacheRawNyp;

  // bugfix, LdLoadNecessary() was otherwise claiming that reload wasn't
  // necessary in certain cases
  pgrp->ldbase_vidx = 0x80000000U;
}

// Design change (30 Nov 2019): It is easy to forget to call PgrClearLdCache
// when changing sample_include.  However, each sample_include change must be
// accompanied by a sample_include_cumulative_popcounts update.  So, if we
// define a sample_include_cumulative_popcounts wrapper-type which can only be
// initialized by a function that also clears a PgenReader LD cache, and modify
// all PgrGet... functions to require this wrapper-type, the frequency of
// foot-shooting should go down.
//
// The key usage rule is: only use this as a local variable type, and define
// only one of these per function (unless you're using multiple PgenReaders
// simultaneously, anyway).  If you're changing the sample-subset when entering
// and exiting chrY, call PgrSetSampleSubsetIndex on your single
// PgrSampleSubsetIndex at the time you're crossing a chrY boundary.  Don't
// define two preinitialized PgrSetSampleSubsetIndexes...
// (possible todo: if compiling as C++ and NDEBUG isn't defined, add a counter
// field to PgenReader which is initialized to zero, asserted to be zero and
// then incremented by PgrSetSampleSubsetIndex, and decremented by the
// PgrSampleSubsetIndex destructor.)
typedef struct PgrSampleSubsetIndexStruct {
#ifdef __cplusplus
  const uint32_t*& GET_PRIVATE_cumulative_popcounts() { return cumulative_popcounts; }
  const uint32_t* const& GET_PRIVATE_cumulative_popcounts() const { return cumulative_popcounts; }
 private:
#endif
  const uint32_t* cumulative_popcounts;
} PgrSampleSubsetIndex;

HEADER_INLINE void PgrSetSampleSubsetIndex(const uint32_t* sample_include_cumulative_popcounts, PgenReader* pgr_ptr, PgrSampleSubsetIndex* pssi_ptr) {
  GET_PRIVATE(*pssi_ptr, cumulative_popcounts) = sample_include_cumulative_popcounts;
  PgrClearLdCache(pgr_ptr);
}

HEADER_INLINE void PgrClearSampleSubsetIndex(PgenReader* pgr_ptr, PgrSampleSubsetIndex* pssi_ptr) {
  GET_PRIVATE(*pssi_ptr, cumulative_popcounts) = nullptr;
  if (pgr_ptr) {
    PgrClearLdCache(pgr_ptr);
  }
}

HEADER_INLINE void PgrSetBaseAndOffset0(unsigned char* block_base, uint32_t thread_ct, PgenReader** pgr_ptr_arr) {
  for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
    PgenReader* pgr_ptr = pgr_ptr_arr[tidx];
    PgrClearLdCache(pgr_ptr);
    PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
    pgrp->fi.block_base = block_base;
    pgrp->fi.block_offset = 0;
  }
}

// PgenFileInfo initialization is split into two phases, to decouple
// plink2's arena allocator from this library.  (obvious todo: provide a simple
// malloc-using PgenReader constructor for anyone who doesn't want to worry
// about these details.)
//
// Phase 1: Open the .pgen (and .pgen.pgi, if relevant); verify that the
//   initial bytes are consistent with the file format; load/verify sample and
//   variant counts, initialize pgfi.const_vrtype, pgfi.const_vrec_width, and
//   pgfi.const_fpos_offset; determine initial memory allocation requirement.
//   pgfi_alloc_cacheline_ct does not include allele counts and nonref flags,
//   since it may be more appropriate to allocate those arrays earlier (during
//   loading of a .bim-like file).
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

// There are two modes of operation:
// 1. fread block-load.  Block-load operations are single-threaded, while
//    decompression/counting is multithreaded.  Appropriate for whole-genome
//    queries, since even with a SSD, reading from multiple parts of a file
//    simultaneously doesn't work well.
// 2. fread single-variant-at-a-time.  Simpler interface than block-load, and
//    doesn't share its inability to handle multiple queries at a time, but
//    less performant for CPU-heavy operations on the whole genome.
// First mode corresponds to use_blockload == 1 in phase2, and second mode
// corresponds to use_blockload == 0.
//
// There was originally a third mmap-based mode, which was removed on 14 Mar
// 2022.  If you are interested in building e.g. a webserver backend that can
// address multiple queries in parallel, refer to plink-ng commit c470317,
// which captures the state of the codebase immediately preceding removal of
// the mmap mode.
//
// Other notes:
// - If pgi_fname is nullptr but the .pgen has an external index file, the
//   index file name is assumed to be the .pgen filename with .pgi appended.
// - pgi_fname is ignored if the .pgen does not have an external index file.
// - raw_variant_ct must be in [1, 2^31 - 3], and raw_sample_ct must be in [1,
//   2^31 - 2].
PglErr PgfiInitPhase1(const char* fname, const char* pgi_fname, uint32_t raw_variant_ct, uint32_t raw_sample_ct, PgenHeaderCtrl* header_ctrl_ptr, PgenFileInfo* pgfip, uintptr_t* pgfi_alloc_cacheline_ct_ptr, char* errstr_buf);

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


void PreinitPgr(PgenReader* pgr_ptr);

// Before PgrInit() is called, the caller must obtain a block of
// pgr_alloc_cacheline_ct * 64 bytes (this value is returned by
// pgfi_init_phase2), 64-byte aligned; this is the pgr_alloc parameter.
//
// There's also a modal usage difference:
//
// * Mode 1 (block-fread): There is one PgenFileInfo per file which doesn't
//   belong to any reader.  After it's initialized, multiple PgenReaders can be
//   based off of it.  When the PgenFileInfo is destroyed, those PgenReaders
//   are invalidated and should be destroyed if that hasn't already happened.
//
//   fname parameter must be nullptr.
//
// * Mode 2 (per-variant fread): Destruction of the original PgenFileInfo
//   struct does not invalidate any extant PgenReader instances (at least
//   from pgenlib_read's perspective).  Instead, destruction of the
//   corresponding memory block or allele_idx_offsets/nonref_flags invalidates
//   the associated PgenReaders.
//
//   The only difference between the first reader and later readers of the same
//   file is that the first reader steals the shared_ff used to read the
//   header.
//
//   fname parameter must be non-null.

// max_vrec_width ignored when using mode 1.
PglErr PgrInit(const char* fname, uint32_t max_vrec_width, PgenFileInfo* pgfip, PgenReader* pgr_ptr, unsigned char* pgr_alloc);

// practically all these functions require genovec to be allocated up to
// vector, not word, boundary
void PgrPlink1ToPlink2InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec);

void PgrPlink2ToPlink1InplaceUnsafe(uint32_t sample_ct, uintptr_t* genovec);

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
//   it.  See SampleCountsThread() in plink2_misc for a usage example.
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
PglErr PgrGet(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, uintptr_t* __restrict genovec);

// Loads the specified variant as a difflist if that's more efficient, setting
// difflist_common_geno to the common genotype value in that case.  Otherwise,
// genovec is populated and difflist_common_geno is set to UINT32_MAX.
//
// max_simple_difflist_len must be smaller than sample_ct.
//
// Note that the returned difflist_len can be much larger than
// max_simple_difflist_len when the variant is LD-encoded; it's bounded by
//   2 * (raw_sample_ct / kPglMaxDifflistLenDivisor).
PglErr PgrGetDifflistOrGenovec(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t max_simple_difflist_len, uint32_t vidx, PgenReader* pgr_ptr, uintptr_t* __restrict genovec, uint32_t* difflist_common_geno_ptr, uintptr_t* __restrict main_raregeno, uint32_t* __restrict difflist_sample_ids, uint32_t* __restrict difflist_len_ptr);

// genocounts[0] = # hom ref, [1] = # het ref, [2] = two alts, [3] = missing
PglErr PgrGetCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, STD_ARRAY_REF(uint32_t, 4) genocounts);

// genocounts[0] = # of hardcalls with two copies of specified allele
// genocounts[1] = # of hardcalls with exactly one copy of specified allele
// genocounts[2] = # of hardcalls with no copies
// genocounts[3] = missing
PglErr PgrGetInv1Counts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgr_ptr, STD_ARRAY_REF(uint32_t, 4) genocounts);

PglErr IMPLPgrGet1(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReaderMain* pgrp, uintptr_t* __restrict allele_countvec);

// Loads a nypvec with counts of a single allele (allele_idx 0 corresponds to
// the reference allele, allele_idx 1 corresponds to alt1, etc.).  0b11 ==
// missing call.
// Note that calling this with allele_idx == 0 is similar to a plink1 load
// (except with missing == 0b11, of course).
// todo: provide a difflist interface once anyone wants it.
HEADER_INLINE PglErr PgrGet1(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgr_ptr, uintptr_t* __restrict allele_countvec) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  const uint32_t* sample_include_cumulative_popcounts = GET_PRIVATE(pssi, cumulative_popcounts);
  return IMPLPgrGet1(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, allele_countvec);
}

PglErr IMPLPgrGetInv1(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReaderMain* pgrp, uintptr_t* __restrict allele_invcountvec);

HEADER_INLINE PglErr PgrGetInv1(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgr_ptr, uintptr_t* __restrict allele_invcountvec) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  const uint32_t* sample_include_cumulative_popcounts = GET_PRIVATE(pssi, cumulative_popcounts);
  return IMPLPgrGetInv1(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, allele_invcountvec);
}

PglErr IMPLPgrGet2(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx0, uint32_t allele_idx1, PgenReaderMain* pgrp, uintptr_t* __restrict genovec);

HEADER_INLINE PglErr PgrGet2(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx0, uint32_t allele_idx1, PgenReader* pgr_ptr, uintptr_t* __restrict genovec) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  const uint32_t* sample_include_cumulative_popcounts = GET_PRIVATE(pssi, cumulative_popcounts);
  return IMPLPgrGet2(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx0, allele_idx1, pgrp, genovec);
}

void PreinitPgv(PgenVariant* pgvp);

PglErr PgrGetM(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, PgenVariant* pgvp);

// possible todo: add functions which directly support MAF-based queries.  Note
// that when the difflist representation is used, we can disqualify some
// low-MAF variants without actually loading the genotype data, since the size
// of the record puts an upper bound on the alt allele frequency.

// requires trailing bits of genoarr to be zeroed out, AND does not update high
// bits of last word if raw_sample_ctl2 is odd.
void DetectGenoarrHetsHw(const uintptr_t*__restrict genoarr, uint32_t raw_sample_ctl2, Halfword* __restrict all_hets_hw);

// requires trailing bits of genoarr to be zeroed out.
HEADER_INLINE void PgrDetectGenoarrHetsUnsafe(const uintptr_t*__restrict genoarr, uint32_t raw_sample_ctl2, uintptr_t* __restrict all_hets) {
  Halfword* all_hets_alias = R_CAST(Halfword*, all_hets);
  DetectGenoarrHetsHw(genoarr, raw_sample_ctl2, all_hets_alias);
  if (raw_sample_ctl2 % 2) {
    all_hets_alias[raw_sample_ctl2] = 0;
  }
}

HEADER_INLINE void PgrDetectGenoarrHets(const uintptr_t* __restrict genoarr, uint32_t raw_sample_ct, uintptr_t* __restrict all_hets) {
  DetectGenoarrHetsHw(genoarr, NypCtToWordCt(raw_sample_ct), R_CAST(Halfword*, all_hets));
  ZeroTrailingBits(raw_sample_ct, all_hets);
}

// sample_ct > 0.  ok for trailing bits of genoarr to not be zeroed out.
void PgrDetectGenoarrHetsMultiallelic(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t raw_sample_ct, uintptr_t* __restrict all_hets);

// cannot assume phaseinfo bit is clear when phasepresent is clear.
PglErr PgrGetP(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr);

PglErr PgrGet1P(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgr_ptr, uintptr_t* __restrict allele_countvec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr);

PglErr IMPLPgrGetInv1P(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReaderMain* pgrp, uintptr_t* __restrict allele_invcountvec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr);

HEADER_INLINE PglErr PgrGetInv1P(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx, PgenReader* pgr_ptr, uintptr_t* __restrict allele_invcountvec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  const uint32_t* sample_include_cumulative_popcounts = GET_PRIVATE(pssi, cumulative_popcounts);
  return IMPLPgrGetInv1P(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, allele_idx, pgrp, allele_invcountvec, phasepresent, phaseinfo, phasepresent_ct_ptr);
}

PglErr PgrGet2P(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, uint32_t allele_idx0, uint32_t allele_idx1, PgenReader* pgr_ptr, uintptr_t* __restrict genovec, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo, uint32_t* __restrict phasepresent_ct_ptr);

PglErr PgrGetMP(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, PgenVariant* pgvp);

PglErr IMPLPgrGetD(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReaderMain* pgrp, uintptr_t* __restrict genovec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr);

// if dosage_present and dosage_main are nullptr, dosage data is ignored
HEADER_INLINE PglErr PgrGetD(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, uintptr_t* __restrict genovec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  const uint32_t* sample_include_cumulative_popcounts = GET_PRIVATE(pssi, cumulative_popcounts);
  return IMPLPgrGetD(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, genovec, dosage_present, dosage_main, dosage_ct_ptr);
}

PglErr PgrGet1D(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReader* pgr_ptr, uintptr_t* __restrict allele_countvec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr);

PglErr PgrGetInv1D(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReader* pgr_ptr, uintptr_t* __restrict allele_invcountvec, uintptr_t* __restrict dosage_present, uint16_t* dosage_main, uint32_t* dosage_ct_ptr);

// When computing either form of imputation-r2, this function requires the
// variant to be biallelic; PgrGetMDCounts must be called in that multiallelic
// case.
// imp_r2_ptr must be non-null when is_minimac3_r2 is set.
PglErr PgrGetDCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, uint32_t is_minimac3_r2, PgenReader* pgr_ptr, double* imp_r2_ptr, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* __restrict all_dosages);

PglErr PgrGetMDCounts(const uintptr_t* __restrict sample_include, const uintptr_t* __restrict sample_include_interleaved_vec, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, uint32_t is_minimac3_r2, PgenReader* pgr_ptr, double* __restrict imp_r2_ptr, uint32_t* __restrict het_ctp, STD_ARRAY_REF(uint32_t, 4) genocounts, uint64_t* __restrict all_dosages);

PglErr PgrGetMD(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, PgenVariant* pgvp);

PglErr IMPLPgrGetDp(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, uint32_t sample_ct, uint32_t vidx, PgenReaderMain* pgrp, PgenVariant* pgvp);

// ok for both dosage_present and dosage_main to be nullptr when no dosage data
// is present
// ok for dphase_present/dphase_delta to be nullptr; dphase_ct always set to 0
// in that case
HEADER_INLINE PglErr PgrGetDp(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, PgenVariant* pgvp) {
  PgenReaderMain* pgrp = &GET_PRIVATE(*pgr_ptr, m);
  const uint32_t* sample_include_cumulative_popcounts = GET_PRIVATE(pssi, cumulative_popcounts);
  return IMPLPgrGetDp(sample_include, sample_include_cumulative_popcounts, sample_ct, vidx, pgrp, pgvp);
}

// pgvp->genovec filled with inverse-counts for specified allele
PglErr PgrGetInv1Dp(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, AlleleCode allele_idx, PgenReader* pgr_ptr, PgenVariant* pgvp);

PglErr PgrGetMDp(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, PgenVariant* pgvp);

// interface used by --make-pgen, just performs basic LD/difflist decompression
// to maximize parallelism
PglErr PgrGetRaw(uint32_t vidx, PgenGlobalFlags read_gflags, PgenReader* pgr_ptr, uintptr_t** loadbuf_iter_ptr, unsigned char* loaded_vrtype_ptr);

PglErr PgrValidate(PgenReader* pgr_ptr, uintptr_t* genovec_buf, char* errstr_buf);

// missingness bit is set iff hardcall is not present (even if dosage info *is*
// present)
PglErr PgrGetMissingness(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, uintptr_t* __restrict missingness, uintptr_t* __restrict genovec_buf);

// either missingness_hc (hardcall) or missingness_dosage must be non-null for
// now
// missingness_dosage must be vector-aligned
PglErr PgrGetMissingnessD(const uintptr_t* __restrict sample_include, PgrSampleSubsetIndex pssi, uint32_t sample_ct, uint32_t vidx, PgenReader* pgr_ptr, uintptr_t* __restrict missingness_hc, uintptr_t* __restrict missingness_dosage, uintptr_t* __restrict hets, uintptr_t* __restrict genovec_buf);


// error-return iff reterr was success and was changed to kPglRetReadFail (i.e.
// an error message should be printed).
BoolErr CleanupPgfi(PgenFileInfo* pgfip, PglErr* reterrp);

BoolErr CleanupPgr(PgenReader* pgr_ptr, PglErr* reterrp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PGENLIB_READ_H__
