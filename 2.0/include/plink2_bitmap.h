#ifndef __PLINK2_BITMAP_H__
#define __PLINK2_BITMAP_H__

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

// Compressed-bitmap-file support.  Similar to, but much simpler than, .pgen.
//
// Header format:
// - Two initial magic bytes, 0x6c 0x1c.
// - Number of rows, then number of columns, both as little-endian uint32s.
// - Number of bytes (in 1..4) used for each row-record-byte-length.
// - .pgen-style row block offsets, all but the last block have 2^16 rows.
// - .pgen-style alternation of 2^16 row record types (2 bits each) and row
//   record byte-lengths.
// The row record types are as follows:
//   0 = Direct storage.
//   1 = Differences-from-earlier-row encoding.  The last row without this type
//       of encoding is the base.  (There is of course no need for an analogue
//       of the "raregeno" field here, or in the last two record types.)
//   2 = Simple difflist, base value is 0.
//   3 = Simple difflist, base value is 1.

#include "plink2_bits.h"

#ifdef __cplusplus
namespace plink2 {
#endif

CONSTI32(kPglRblockSize, 65536);
CONSTI32(kPglMaxBitmapDifflistLenDivisor, 16);
CONSTI32(kPglMaxBitmapBytesPerRow, 0x20000000);

FLAGSET_DEF_START()
  kfPgrPrevdiffCache0,
  kfPgrPrevdiffCacheBit = (1 << 0),
  kfPgrPrevdiffCacheRawBit = (1 << 1),
  kfPgrPrevdiffCachePopcount = (1 << 2)
FLAGSET_DEF_END(PgrPrevdiffCacheFlags);

typedef struct BitmapReaderMainStruct {
  NONCOPYABLE(BitmapReaderMainStruct);
  // ----- Header information, constant after initialization -----
  uint32_t row_ct;
  uint32_t col_ct;

  // size (row_ct + 1) so that the number of bytes of (zero-based) row n is
  // row_fpos[n+1] - row_fpos[n].
  // --pmerge can *almost* use a sequential reader, which would let us get away
  // with a much smaller memory allocation here.  Unfortunately, same-position
  // variants may need to be reordered by ID, and that makes it simpler to
  // stick to a random-access API with this entire array preloaded.
  uint64_t* row_fpos;

  uintptr_t* rrtype_nyparr;

  // ----- Mutable state -----
  // If we don't fseek, what's the next row we'd read?
  uint32_t fp_ridx;

  FILE* ff;
  // allocated to word boundary, any trailing bytes always zeroed
  unsigned char* fread_buf;

  uint32_t prevdiff_ridx;
  PgrPrevdiffCacheFlags prevdiff_stypes;

  uintptr_t* prevdiff_base_raw_bitvec; // trailing bytes always zeroed
  uintptr_t* prevdiff_base_bitvec; // after col_include applied
  // should add prevdiff_base_col_ids (and uint32_t list length) if we add a
  // DifflistOrBitvec getter

  // might want workspace and/or prevdiff_base_raw_bitvec popcount
} BitmapReaderMain;

typedef struct BitmapReaderStruct {
#ifdef __cplusplus
  BitmapReaderMain& GET_PRIVATE_m() { return m; }
  BitmapReaderMain const& GET_PRIVATE_m() const { return m; }
 private:
#endif
  BitmapReaderMain m;
} BitmapReader;

HEADER_INLINE uintptr_t GetRrtype(const BitmapReaderMain* brp, uint32_t ridx) {
  return GetNyparrEntry(brp->rrtype_nyparr, ridx);
}

HEADER_INLINE uint64_t GetBrpFpos(const BitmapReaderMain* brp, uint32_t ridx) {
  return brp->row_fpos[ridx];
}

HEADER_INLINE uint64_t GetBrpRrecWidth(const BitmapReaderMain* brp, uint32_t ridx) {
  return brp->row_fpos[ridx + 1] - brp->row_fpos[ridx];
}

HEADER_INLINE uint32_t RrtypeIsPrevdiff(uint32_t rrtype) {
  return (rrtype == 1);
}

HEADER_INLINE void PgrClearPrevdiffCache(BitmapReader* br_ptr) {
  BitmapReaderMain* brp = &GET_PRIVATE(*br_ptr, m);
  brp->prevdiff_stypes &= kfPgrPrevdiffCacheRawBit;
  brp->prevdiff_ridx = 0x80000000U;
}

typedef struct PgrColSubsetIndexStruct {
#ifdef __cplusplus
  const uint32_t*& GET_PRIVATE_cumulative_popcounts() { return cumulative_popcounts; }
  const uint32_t* const& GET_PRIVATE_cumulative_popcounts() const { return cumulative_popcounts; }
 private:
#endif
  const uint32_t* cumulative_popcounts;
} PgrColSubsetIndex;

HEADER_INLINE void PgrSetColSubsetIndex(const uint32_t* col_include_cumulative_popcounts, BitmapReader* br_ptr, PgrColSubsetIndex* pcsi_ptr) {
  GET_PRIVATE(*pcsi_ptr, cumulative_popcounts) = col_include_cumulative_popcounts;
  PgrClearPrevdiffCache(br_ptr);
}

HEADER_INLINE void PgrClearColSubsetIndex(BitmapReader* br_ptr, PgrColSubsetIndex* pcsi_ptr) {
  GET_PRIVATE(*pcsi_ptr, cumulative_popcounts) = nullptr;
  if (br_ptr) {
    PgrClearPrevdiffCache(br_ptr);
  }
}

void PreinitBitmapReader(BitmapReader* br_ptr);

PglErr BitmapReaderInitPhase1(const char* fname, BitmapReader* br_ptr, uintptr_t* br_alloc_cacheline_ct_ptr, char* errstr_buf);

PglErr BitmapReaderInitPhase2(BitmapReader* br_ptr, unsigned char* br_alloc, char* errstr_buf);

HEADER_INLINE uint32_t BitmapReaderRowCt(const BitmapReader* br_ptr) {
  const BitmapReaderMain* brp = &GET_PRIVATE(*br_ptr, m);
  return brp->row_ct;
}

HEADER_INLINE uint32_t BitmapReaderColCt(const BitmapReader* br_ptr) {
  const BitmapReaderMain* brp = &GET_PRIVATE(*br_ptr, m);
  return brp->col_ct;
}

PglErr BitmapGet(const uintptr_t* __restrict col_include, PgrColSubsetIndex pcsi, uint32_t col_ct, uint32_t ridx, BitmapReader* br_ptr, uintptr_t* __restrict dst);

typedef struct BitmapWriterMainStruct {
  NONCOPYABLE(BitmapWriterMainStruct);
  uint32_t row_ct;
  uint32_t col_ct;

  uint64_t* rblock_fpos;
  unsigned char* rrec_len_buf;
  uintptr_t* rrtype_buf;

  uint32_t prevdiff_base_set_bit_ct;

  // should match ftello() return value
  uint64_t rblock_fpos_offset;

  // these must hold col_ct entries
  uintptr_t* difflist_bitvec_buf;
  uintptr_t* prevdiff_base_bitvec;

  unsigned char* fwrite_buf;
  unsigned char* fwrite_bufp;

  uint32_t prevdiff_list_len;

  uintptr_t rrec_len_byte_ct;

  uint32_t ridx;

  FILE* ff;
} BitmapWriterMain;

typedef struct BitmapWriterStruct {
#ifdef __cplusplus
  BitmapWriterMain& GET_PRIVATE_m() { return m; }
  BitmapWriterMain const& GET_PRIVATE_m() const { return m; }
 private:
#endif
  BitmapWriterMain m;
} BitmapWriter;

void PreinitBitmapWriter(BitmapWriter* bw_ptr);

PglErr BitmapWriterInitPhase1(const char* fname, uint32_t row_ct, uint32_t col_ct, BitmapWriter* bw_ptr, uintptr_t* bw_alloc_cacheline_ct_ptr);

PglErr BitmapWriterInitPhase2(BitmapWriter* bw_ptr, unsigned char* bw_alloc);

// trailing bits of bitvec must be zeroed out
PglErr BitmapAppend(const uintptr_t* bitvec, BitmapWriter* bw_ptr);

PglErr BitmapWriterFinish(BitmapWriter* bw_ptr);

BoolErr CleanupBitmapReader(BitmapReader* br_ptr, PglErr* reterrp);

BoolErr CleanupBitmapWriter(BitmapWriter* bw_ptr, PglErr* reterrp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_BITMAP_H__
