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


#include "plink2_bitmap.h"

#include <errno.h>

#ifdef __cplusplus
namespace plink2 {
#endif

static inline BitmapReaderMain* GetBrp(BitmapReader* br_ptr) {
  return &GET_PRIVATE(*br_ptr, m);
}

static inline const uint32_t* GetCicp(PgrColSubsetIndex pcsi) {
  return GET_PRIVATE(pcsi, cumulative_popcounts);
}

static inline BitmapWriterMain* GetBwp(BitmapWriter* bw_ptr) {
  return &GET_PRIVATE(*bw_ptr, m);
}

void PreinitBitmapReader(BitmapReader* br_ptr) {
  BitmapReaderMain* brp = GetBrp(br_ptr);
  brp->ff = nullptr;
}

void PreinitBitmapWriter(BitmapWriter* bw_ptr) {
  BitmapWriterMain* bwp = GetBwp(bw_ptr);
  bwp->ff = nullptr;
}

uintptr_t CountBrAllocCachelinesRequired(uint32_t row_ct, uint32_t col_ct) {
  // row_fpos: sizeof(int64_t) * (row_ct + 1) bytes
  // rrtype_nyparr: DivUp(row_ct + 1, 4) bytes
  // fread_buf, prevdiff_base_raw_bitvec, prevdiff_base_bitvec:
  //   DivUp(col_ct, CHAR_BIT) bytes each

  uintptr_t cachelines_required = DivUp(row_ct + 1, kInt64PerCacheline);
  cachelines_required += DivUp(row_ct + 1, kNypsPerCacheline);
  cachelines_required += 3 * DivUp(col_ct, kBitsPerCacheline);
  return cachelines_required;
}

PglErr BitmapReaderInitPhase1(const char* fname, BitmapReader* br_ptr, uintptr_t* br_alloc_cacheline_ct_ptr, char* errstr_buf) {
  // Open the bitmap; verify that the initial bytes are consistent with the
  // file format; load row/column counts; determine initial memory allocation
  // requirement.
  BitmapReaderMain* brp = GetBrp(br_ptr);
  FILE* ff = fopen(fname, FOPEN_RB);
  brp->ff = ff;
  if (unlikely(!ff)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Failed to open %s : %s.\n", fname, strerror(errno));
    return kPglRetOpenFail;
  }
  unsigned char small_readbuf[10];
  if (unlikely(!fread_unlocked(small_readbuf, 10, 1, ff))) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s read failure: %s.\n", fname, strerror(errno));
    return kPglRetReadFail;
  }
  if (unlikely(!memequal_k(small_readbuf, "l\x1c", 2))) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: %s is not a plink-bitmap (first two bytes don't match the magic number).\n", fname);
    return kPglRetMalformedInput;
  }
  memcpy(&(brp->row_ct), &(small_readbuf[2]), sizeof(int32_t));
  memcpy(&(brp->col_ct), &(small_readbuf[6]), sizeof(int32_t));

  *br_alloc_cacheline_ct_ptr = CountBrAllocCachelinesRequired(brp->row_ct, brp->col_ct);
  return kPglRetSuccess;
}

void FillBrReadErrstrFromNzErrno(char* errstr_buf) {
  snprintf(errstr_buf, kPglErrstrBufBlen, "Error: plink-bitmap read failure: %s.\n", strerror(errno));
}

void FillBrReadErrstrFromErrno(char* errstr_buf) {
  if (errno) {
    FillBrReadErrstrFromNzErrno(errstr_buf);
  } else {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: plink-bitmap read failure: File appears to be corrupted.\n");
  }
}

void FillBrReadErrstr(FILE* ff, char* errstr_buf) {
  if (feof_unlocked(ff)) {
    errno = 0;
  }
  FillBrReadErrstrFromErrno(errstr_buf);
}

PglErr BitmapReaderInitPhase2(BitmapReader* br_ptr, unsigned char* br_alloc, char* errstr_buf) {
  BitmapReaderMain* brp = GetBrp(br_ptr);
  const uintptr_t row_ct = brp->row_ct;
  const uintptr_t col_ct = brp->col_ct;
  // Allocate data structures from the br_alloc workspace.
  unsigned char* br_alloc_iter = br_alloc;
  brp->row_fpos = S_CAST(uint64_t*, arena_alloc_raw_rd((row_ct + 1) * sizeof(int64_t), &br_alloc_iter));
  brp->rrtype_nyparr = S_CAST(uintptr_t*, arena_alloc_raw((1 + (row_ct / kNypsPerCacheline)) * kCacheline, &br_alloc_iter));
  const uintptr_t bitvec_byte_ct = DivUp(col_ct, kBitsPerCacheline) * kCacheline;
  brp->fread_buf = S_CAST(unsigned char*, arena_alloc_raw(bitvec_byte_ct, &br_alloc_iter));
  const uint32_t raw_byte_ct = DivUp(col_ct, CHAR_BIT);
  memset(&(brp->fread_buf[raw_byte_ct]), 0, (-raw_byte_ct) & (kBytesPerWord - 1));
  brp->prevdiff_base_raw_bitvec = S_CAST(uintptr_t*, arena_alloc_raw(bitvec_byte_ct, &br_alloc_iter));
  if (col_ct % kBitsPerWord) {
    brp->prevdiff_base_raw_bitvec[(col_ct - 1) / kBitsPerWord] = 0;
  }
  brp->prevdiff_base_bitvec = S_CAST(uintptr_t*, arena_alloc_raw(bitvec_byte_ct, &br_alloc_iter));

  brp->fp_ridx = 0;
  brp->prevdiff_ridx = UINT32_MAX;
  brp->prevdiff_stypes = kfPgrPrevdiffCache0;

  FILE* ff = brp->ff;
  if (!row_ct) {
    // Don't need to make this an error for now.
    // (Advance file position to end of header so that fp_ridx is still
    // accurate.)
    if (unlikely(fseeko(ff, 1, SEEK_CUR))) {
      FillBrReadErrstrFromNzErrno(errstr_buf);
      return kPglRetReadFail;
    }
    return kPglRetSuccess;
  }

  // Note that this is a rather hefty stack allocation.
  unsigned char loadbuf[kPglRblockSize * 4];

  if (unlikely(!fread_unlocked(&(loadbuf[0]), 9, 1, ff))) {
    FillBrReadErrstr(ff, errstr_buf);
    return kPglRetReadFail;
  }
  const uint32_t row_record_byte_length = loadbuf[0];
  if (unlikely((row_record_byte_length - 1) >= 4)) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid plink-bitmap header.\n");
    return kPglRetMalformedInput;
  }

  uint64_t cur_fpos;
  memcpy(&cur_fpos, &(loadbuf[1]), sizeof(int64_t));
  const uint32_t rblock_ct_m1 = (row_ct - 1) / kPglRblockSize;
  // Don't need to load the other row block offsets, they're just present to
  // support cheaper random access.
  if (unlikely(fseeko(ff, rblock_ct_m1 * sizeof(int64_t), SEEK_CUR))) {
    FillBrReadErrstrFromNzErrno(errstr_buf);
    return kPglRetReadFail;
  }
  uintptr_t* rrtype_nyparr_iter = brp->rrtype_nyparr;
  // Must guarantee a trailing zero for is_prevdiff check to work.
  rrtype_nyparr_iter[row_ct / kBitsPerWordD2] = 0;
  uint64_t* row_fpos_iter = brp->row_fpos;
  uint32_t cur_rblock_row_ct = kPglRblockSize;
  for (uint32_t rblock_idx = 0; ; ++rblock_idx) {
    if (rblock_idx >= rblock_ct_m1) {
      if (rblock_idx > rblock_ct_m1) {
        break;
      }
      cur_rblock_row_ct = ModNz(row_ct, kPglRblockSize);
    }
    // Load row record types (2 bits each, 4 entries packed per byte)
    uint32_t cur_byte_ct = DivUp(cur_rblock_row_ct, 4);
    if (unlikely(!fread_unlocked(rrtype_nyparr_iter, cur_byte_ct, 1, ff))) {
      FillBrReadErrstr(ff, errstr_buf);
      return kPglRetReadFail;
    }
    rrtype_nyparr_iter = &(rrtype_nyparr_iter[cur_rblock_row_ct / kBitsPerWordD2]);

    // Load row record lengths.
    cur_byte_ct = cur_rblock_row_ct * row_record_byte_length;
    if (unlikely(!fread_unlocked(loadbuf, cur_byte_ct, 1, ff))) {
      FillBrReadErrstr(ff, errstr_buf);
      return kPglRetReadFail;
    }

    // Identical to fpos-array-filling code in PgfiInitPhase2().
    if (row_record_byte_length == 1) {
      for (uint32_t cur_rblock_ridx = 0; cur_rblock_ridx != cur_rblock_row_ct; ++cur_rblock_ridx) {
        row_fpos_iter[cur_rblock_ridx] = cur_fpos;
        uint32_t cur_rrec_len = loadbuf[cur_rblock_ridx];
        cur_fpos += cur_rrec_len;
      }
    } else if (row_record_byte_length == 2) {
      for (uint32_t cur_rblock_ridx = 0; cur_rblock_ridx != cur_rblock_row_ct; ++cur_rblock_ridx) {
        row_fpos_iter[cur_rblock_ridx] = cur_fpos;
        uint16_t cur_rrec_len;
        memcpy_k(&cur_rrec_len, &(loadbuf[cur_rblock_ridx * 2]), 2);
        cur_fpos += cur_rrec_len;
      }
    } else if (row_record_byte_length == 3) {
      for (uint32_t cur_rblock_ridx = 0; cur_rblock_ridx != cur_rblock_row_ct; ++cur_rblock_ridx) {
        row_fpos_iter[cur_rblock_ridx] = cur_fpos;
        uint32_t cur_rrec_len;
        // safe to read a byte past the end, since that's in loadbuf
        memcpy(&cur_rrec_len, &(loadbuf[cur_rblock_ridx * 3]), sizeof(int32_t));
        cur_rrec_len &= 0xffffff;
        cur_fpos += cur_rrec_len;
      }
    } else {
      for (uint32_t cur_rblock_ridx = 0; cur_rblock_ridx != cur_rblock_row_ct; ++cur_rblock_ridx) {
        row_fpos_iter[cur_rblock_ridx] = cur_fpos;
        uint32_t cur_rrec_len;
        memcpy(&cur_rrec_len, &(loadbuf[cur_rblock_ridx * 4]), 4);
        cur_fpos += cur_rrec_len;
        if (unlikely(cur_rrec_len > kPglMaxBitmapBytesPerRow)) {
          snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid plink-bitmap header.\n");
          return kPglRetMalformedInput;
        }
      }
    }
    row_fpos_iter = &(row_fpos_iter[cur_rblock_row_ct]);
  }

  if (unlikely(S_CAST(uint64_t, ftello(ff)) > brp->row_fpos[0])) {
    snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid plink-bitmap header.\n");
    return kPglRetMalformedInput;
  }

  const uint64_t actual_fpos = ftello(ff);
  if (actual_fpos != brp->row_fpos[0]) {
    if (unlikely(actual_fpos > brp->row_fpos[0])) {
      snprintf(errstr_buf, kPglErrstrBufBlen, "Error: Invalid plink-bitmap header.\n");
      return kPglRetMalformedInput;
    }
    if (unlikely(fseeko(ff, brp->row_fpos[0], SEEK_SET))) {
      FillBrReadErrstrFromNzErrno(errstr_buf);
      return kPglRetReadFail;
    }
  }
  brp->row_fpos[row_ct] = cur_fpos;
  return kPglRetSuccess;
}

uint32_t GetPrevdiffRidx(const uintptr_t* rrtype_nyparr, uint32_t cur_ridx) {
  // Assuming the current row is prevdiff-encoded, this function determines
  // the index of the base row by looking backward in rrtype_nyparr.
  // More precisely, we find the first preceding nyp where either the low bit
  // is unset, or the high bit is set.
  uint32_t widx = (cur_ridx - 1) / kBitsPerWordD2;
  uintptr_t cur_word = rrtype_nyparr[widx];
  uintptr_t non_01 = ((cur_word >> 1) | (~cur_word)) & kMask5555;
  // On the first iteration, we need to mask out high bits corresponding to
  // indexes >= cur_ridx.
  const uint32_t neg_remainder = (-cur_ridx) & (kBitsPerWordD2 - 1);
  non_01 = non_01 & ((~k0LU) >> (2 * neg_remainder));
  while (!non_01) {
    cur_word = rrtype_nyparr[--widx];
    non_01 = ((cur_word >> 1) | (~cur_word)) & kMask5555;
  }
  return (widx * kBitsPerWordD2) + (bsrw(non_01) / 2);
}

uint32_t PrevdiffLoadNecessary(uint32_t cur_ridx, BitmapReaderMain* brp) {
  // Assuming the current row is prevdiff-encoded, this function determines
  // whether the base row also needs to be loaded.
  // Important: this updates brp->prevdiff_ridx when necessary, as a side
  // effect.

  // If the last row we read was (cur_ridx - 1), and the cache contains a row,
  // it must be the correct row.
  if (brp->prevdiff_stypes && (cur_ridx == brp->fp_ridx)) {
    assert(brp->prevdiff_stypes & (kfPgrPrevdiffCacheBit | kfPgrPrevdiffCacheRawBit));
    return 0;
  }

  // Ok, we need to compute the index of the base row, and compare it to
  // brp->prevdiff_ridx.
  const uint32_t old_prevdiff_ridx = brp->prevdiff_ridx;
  brp->prevdiff_ridx = GetPrevdiffRidx(brp->rrtype_nyparr, cur_ridx);
  return (brp->prevdiff_ridx != old_prevdiff_ridx);
}

BoolErr ReadRawBitmapRow(uint32_t ridx, BitmapReaderMain* brp, uint32_t* rrec_width_ptr) {
  // Fills brp->fread_buf.
  if (brp->fp_ridx != ridx) {
    if (unlikely(fseeko(brp->ff, GetBrpFpos(brp, ridx), SEEK_SET))) {
      return 1;
    }
  }
  const uint32_t rrec_width = GetBrpRrecWidth(brp, ridx);
  if (unlikely(!fread_unlocked(brp->fread_buf, rrec_width, 1, brp->ff))) {
    if (feof_unlocked(brp->ff)) {
      errno = 0;
    }
    return 1;
  }
  *rrec_width_ptr = rrec_width;
  brp->fp_ridx = ridx + 1;
  return 0;
}

PglErr ParseBitmapDifflistHeader(const unsigned char* fread_end, uint32_t raw_col_ct, const unsigned char** fread_pp, const unsigned char** difflist_group_info_ptr, uint32_t* difflist_len_ptr) {
  // Identical to ParseDifflistHeader with raregeno_buf == nullptr, except for
  // divisor constant.
  const uint32_t difflist_len = GetVint31(fread_end, fread_pp);
  // moved here to address maybe-uninitialized warnings
  *difflist_group_info_ptr = *fread_pp;
  *difflist_len_ptr = difflist_len;
  if (!difflist_len) {
    return kPglRetSuccess;
  }
  if (unlikely(difflist_len > raw_col_ct / kPglMaxBitmapDifflistLenDivisor)) {
    // automatically catches GetVint31() failure
    return kPglRetMalformedInput;
  }
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  const uint32_t col_id_byte_ct = BytesToRepresentNzU32(raw_col_ct);
  const uint32_t difflist_index_byte_ct = group_ct * (col_id_byte_ct + 1) - 1;
  if (PtrAddCk(fread_end, difflist_index_byte_ct, fread_pp)) {
    return kPglRetMalformedInput;
  }
  return kPglRetSuccess;
}

// dst is in/out
PglErr ParseAndApplyBitmapDifflist(uint32_t rrec_width, BitmapReaderMain* brp, uintptr_t* dst) {
  // Cannot occur after bitvec subsetting since the difflist sample indexes
  // will be incorrect.
  const uint32_t raw_col_ct = brp->col_ct;
  const unsigned char* fread_ptr = brp->fread_buf;
  const unsigned char* fread_end = &(fread_ptr[rrec_width]);
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseBitmapDifflistHeader(fread_end, raw_col_ct, &fread_ptr, &group_info_iter, &difflist_len);
  if (reterr || (!difflist_len)) {
    return reterr;
  }
  const uint32_t col_id_byte_ct = BytesToRepresentNzU32(raw_col_ct);
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kPglDifflistGroupSize;
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
    uint32_t remaining_deltas_in_subgroup = kPglDifflistGroupSize - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        break;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    uint32_t raw_col_idx = SubU32Load(group_info_iter, col_id_byte_ct);
    group_info_iter = &(group_info_iter[col_id_byte_ct]);
    for (; ; --remaining_deltas_in_subgroup) {
      // always check, since otherwise subsequent assignment can mutate
      // arbitrary memory
      if (unlikely(raw_col_idx >= raw_col_ct)) {
        return kPglRetMalformedInput;
      }
      dst[raw_col_idx / kBitsPerWord] ^= k1LU << (raw_col_idx % kBitsPerWord);
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_col_idx += GetVint31(fread_end, &fread_ptr);
    }
  }
  return kPglRetSuccess;
}

// dst is in/out
PglErr ParseAndApplyBitmapDifflistSubset(const uintptr_t* __restrict col_include, const uint32_t* __restrict col_include_cumulative_popcounts, uint32_t rrec_width, uint32_t col_ct, BitmapReaderMain* brp, uintptr_t* __restrict dst) {
  const uint32_t raw_col_ct = brp->col_ct;
  if (col_ct == raw_col_ct) {
    return ParseAndApplyBitmapDifflist(rrec_width, brp, dst);
  }
  const unsigned char* fread_ptr = brp->fread_buf;
  const unsigned char* fread_end = &(fread_ptr[rrec_width]);
  const unsigned char* group_info_iter;
  uint32_t difflist_len;
  PglErr reterr = ParseBitmapDifflistHeader(fread_end, raw_col_ct, &fread_ptr, &group_info_iter, &difflist_len);
  if (reterr || (!difflist_len)) {
    return reterr;
  }
  const uint32_t col_id_byte_ct = BytesToRepresentNzU32(raw_col_ct);
  const uint32_t subgroup_idx_last = (difflist_len - 1) / kPglDifflistGroupSize;
  for (uint32_t subgroup_idx = 0; ; ++subgroup_idx) {
    uint32_t remaining_deltas_in_subgroup = kPglDifflistGroupSize - 1;
    if (subgroup_idx >= subgroup_idx_last) {
      if (subgroup_idx > subgroup_idx_last) {
        break;
      }
      remaining_deltas_in_subgroup &= difflist_len - 1;
    }
    uint32_t raw_col_idx = SubU32Load(group_info_iter, col_id_byte_ct);
    group_info_iter = &(group_info_iter[col_id_byte_ct]);
    for (; ; --remaining_deltas_in_subgroup) {
      // always check, since otherwise subsequent assignment can mutate
      // arbitrary memory
      if (unlikely(raw_col_idx >= raw_col_ct)) {
        return kPglRetMalformedInput;
      }
      if (IsSet(col_include, raw_col_idx)) {
        const uint32_t col_idx = RawToSubsettedPos(col_include, col_include_cumulative_popcounts, raw_col_idx);
        dst[col_idx / kBitsPerWord] ^= k1LU << (col_idx % kBitsPerWord);
      }
      if (!remaining_deltas_in_subgroup) {
        break;
      }
      raw_col_idx += GetVint31(fread_end, &fread_ptr);
    }
  }
  return kPglRetSuccess;
}

PglErr ParseNonPrevdiffBitvecSubset(const uintptr_t* __restrict col_include, const uint32_t* __restrict col_include_cumulative_popcounts, uint32_t rrec_width, uint32_t col_ct, uint32_t rrtype, BitmapReaderMain* brp, uintptr_t* __restrict dst) {
  // Fills brp->prevdiff_base_raw_bitvec iff rrtype == 0 and
  // subsetting_required.
  // (Does not update prevdiff_stypes, caller's responsibility to care)
  const uint32_t raw_col_ct = brp->col_ct;
  if (!rrtype) {
    // Uncompressed storage.
    // assert(rrec_width == DivUp(raw_col_ct, CHAR_BIT));
    const uintptr_t* raw_bitvec = R_CAST(const uintptr_t*, brp->fread_buf);
    if (col_ct == raw_col_ct) {
      // No subsetting required.
      CopyBitarr(raw_bitvec, col_ct, dst);
    } else {
      CopyBitarr(raw_bitvec, raw_col_ct, brp->prevdiff_base_raw_bitvec);
      CopyBitarrSubset(raw_bitvec, col_include, col_ct, dst);
    }
    return kPglRetSuccess;
  }
  if (rrtype == 2) {
    ZeroWArr(BitCtToWordCt(col_ct), dst);
  } else {
    SetAllBits(col_ct, dst);
  }
  return ParseAndApplyBitmapDifflistSubset(col_include, col_include_cumulative_popcounts, rrec_width, col_ct, brp, dst);
}

PglErr PrevdiffLoadAndCopyBitvecSubset(const uintptr_t* __restrict col_include, const uint32_t* __restrict col_include_cumulative_popcounts, uint32_t col_ct, uint32_t ridx, BitmapReaderMain* brp, uintptr_t* __restrict dst) {
  const uint32_t raw_col_ct = brp->col_ct;
  if (PrevdiffLoadNecessary(ridx, brp)) {
    const uint32_t prevdiff_ridx = brp->prevdiff_ridx;
    uint32_t rrec_width;
    if (unlikely(ReadRawBitmapRow(prevdiff_ridx, brp, &rrec_width))) {
      return kPglRetReadFail;
    }
    const uint32_t rrtype = GetRrtype(brp, prevdiff_ridx);
    PglErr reterr = ParseNonPrevdiffBitvecSubset(col_include, col_include_cumulative_popcounts, rrec_width, col_ct, rrtype, brp, dst);
    brp->prevdiff_stypes = ((col_ct != raw_col_ct) && (rrtype == 0))? (kfPgrPrevdiffCacheBit | kfPgrPrevdiffCacheRawBit) : kfPgrPrevdiffCacheBit;
    CopyBitarr(dst, col_ct, brp->prevdiff_base_bitvec);
    return reterr;
  }
  if (brp->prevdiff_stypes & kfPgrPrevdiffCacheBit) {
    CopyBitarr(brp->prevdiff_base_bitvec, col_ct, dst);
  } else {
    assert(brp->prevdiff_stypes & kfPgrPrevdiffCacheRawBit);
    if (col_ct == raw_col_ct) {
      CopyBitarr(brp->prevdiff_base_bitvec, col_ct, dst);
    } else {
      CopyBitarrSubset(brp->prevdiff_base_raw_bitvec, col_include, col_ct, dst);
      CopyBitarr(dst, col_ct, brp->prevdiff_base_bitvec);
      brp->prevdiff_stypes |= kfPgrPrevdiffCacheBit;
    }
  }
  return kPglRetSuccess;
}

PglErr BitmapGet(const uintptr_t* __restrict col_include, PgrColSubsetIndex pcsi, uint32_t col_ct, uint32_t ridx, BitmapReader* br_ptr, uintptr_t* __restrict dst) {
  if (!col_ct) {
    return kPglRetSuccess;
  }
  BitmapReaderMain* brp = GetBrp(br_ptr);
  assert(ridx < brp->row_ct);
  const uint32_t* col_include_cumulative_popcounts = GetCicp(pcsi);

  const uint32_t rrtype = GetRrtype(brp, ridx);
  if (RrtypeIsPrevdiff(rrtype)) {
    PglErr reterr = PrevdiffLoadAndCopyBitvecSubset(col_include, col_include_cumulative_popcounts, col_ct, ridx, brp, dst);
    if (unlikely(reterr)) {
      return reterr;
    }
    uint32_t rrec_width;
    if (unlikely(ReadRawBitmapRow(ridx, brp, &rrec_width))) {
      return kPglRetReadFail;
    }
    return ParseAndApplyBitmapDifflistSubset(col_include, col_include_cumulative_popcounts, rrec_width, col_ct, brp, dst);
  }
  uint32_t rrec_width;
  if (unlikely(ReadRawBitmapRow(ridx, brp, &rrec_width))) {
    return kPglRetReadFail;
  }
  PglErr reterr = ParseNonPrevdiffBitvecSubset(col_include, col_include_cumulative_popcounts, rrec_width, col_ct, rrtype, brp, dst);
  if (unlikely(reterr)) {
    return reterr;
  }
  const uint32_t is_prevdiff_base = RrtypeIsPrevdiff(GetRrtype(brp, ridx + 1));
  const uint32_t prevdiff_raw_bitvec_saved = (col_ct != brp->col_ct) && (rrtype == 0);
  if (is_prevdiff_base) {
    CopyBitarr(dst, col_ct, brp->prevdiff_base_bitvec);
    brp->prevdiff_ridx = ridx;
    brp->prevdiff_stypes = prevdiff_raw_bitvec_saved? (kfPgrPrevdiffCacheBit | kfPgrPrevdiffCacheRawBit) : kfPgrPrevdiffCacheBit;
  } else if (prevdiff_raw_bitvec_saved) {
    // we just clobbered the cache
    brp->prevdiff_stypes &= ~kfPgrPrevdiffCacheRawBit;
  }
  return kPglRetSuccess;
}

uintptr_t CountBwAllocCachelinesRequired(uint32_t row_ct, uint32_t col_ct) {
  // rblock_fpos
  const uint32_t rblock_ct = DivUp(row_ct, kPglRblockSize);
  uintptr_t cachelines_required = Int64CtToCachelineCt(rblock_ct);

  // rrec_len_buf
  const uint32_t max_rrec_len = DivUp(col_ct, CHAR_BIT);
  const uintptr_t rrec_len_byte_ct = BytesToRepresentNzU32(max_rrec_len);
  cachelines_required += DivUp(row_ct * rrec_len_byte_ct, kCacheline);

  // rrtype_buf
  cachelines_required += NypCtToCachelineCt(row_ct);

  // difflist_bitvec_buf, prevdiff_base_bitvec
  cachelines_required += 2 * BitCtToCachelineCt(col_ct);

  // fwrite_buf
  cachelines_required += DivUp(max_rrec_len + kPglFwriteBlockSize, kCacheline);
  return cachelines_required;
}

PglErr BitmapWriterInitPhase1(const char* fname, uint32_t row_ct, uint32_t col_ct, BitmapWriter* bw_ptr, uintptr_t* bw_alloc_cacheline_ct_ptr) {
  assert(row_ct);
  assert(col_ct);
  BitmapWriterMain* bwp = GetBwp(bw_ptr);

  bwp->row_ct = row_ct;
  bwp->col_ct = col_ct;
#ifndef NDEBUG
  bwp->rblock_fpos = nullptr;
  bwp->rrec_len_buf = nullptr;
  bwp->rrtype_buf = nullptr;
  bwp->difflist_bitvec_buf = nullptr;
  bwp->prevdiff_base_bitvec = nullptr;
  bwp->fwrite_buf = nullptr;
  bwp->fwrite_bufp = nullptr;
#endif
  bwp->ridx = 0;

  FILE* ff = fopen(fname, FOPEN_WB);
  bwp->ff = ff;
  if (unlikely(!ff)) {
    return kPglRetOpenFail;
  }
  fwrite_unlocked("l\x1c", 2, 1, ff);
  fwrite_unlocked(&row_ct, sizeof(int32_t), 1, ff);
  fwrite_unlocked(&col_ct, sizeof(int32_t), 1, ff);
  const uint32_t max_rrec_len = DivUp(col_ct, CHAR_BIT);
  const unsigned char rrec_len_byte_ct = BytesToRepresentNzU32(max_rrec_len);
  bwp->rrec_len_byte_ct = rrec_len_byte_ct;
  if (unlikely(putc_checked(rrec_len_byte_ct, ff))) {
    return kPglRetWriteFail;
  }
  const uint32_t rblock_ct = DivUp(row_ct, kPglRblockSize);
  uintptr_t header_bytes_left = rblock_ct * sizeof(int64_t) + row_ct * S_CAST(uintptr_t, rrec_len_byte_ct) + DivUp(row_ct, 4);
  // this should be the position of the first row
  bwp->rblock_fpos_offset = 11 + header_bytes_left;

  uintptr_t zeroed_cachelines_needed = DivUp(header_bytes_left, kCacheline);
  if (zeroed_cachelines_needed > (kPglFwriteBlockSize / kCacheline)) {
    zeroed_cachelines_needed = kPglFwriteBlockSize / kCacheline;
  }
  // could wait until fwrite_buf is allocated, and make sure it's aligned?
  unsigned char zerobuf[kPglFwriteBlockSize];
  memset(zerobuf, 0, zeroed_cachelines_needed * kCacheline);
  while (header_bytes_left > kPglFwriteBlockSize) {
    fwrite_unlocked(zerobuf, kPglFwriteBlockSize, 1, ff);
    header_bytes_left -= kPglFwriteBlockSize;
  }
  if (unlikely(fwrite_checked(zerobuf, header_bytes_left, ff))) {
    return kPglRetWriteFail;
  }

  *bw_alloc_cacheline_ct_ptr = CountBwAllocCachelinesRequired(row_ct, col_ct);
  return kPglRetSuccess;
}

PglErr BitmapWriterInitPhase2(BitmapWriter* bw_ptr, unsigned char* bw_alloc) {
  BitmapWriterMain* bwp = GetBwp(bw_ptr);
  const uint32_t row_ct = bwp->row_ct;
  unsigned char* alloc_iter = bw_alloc;
  const uint32_t rblock_ct = DivUp(row_ct, kPglRblockSize);
  const uint32_t rrtype_buf_byte_ct = NypCtToCachelineCt(row_ct) * kCacheline;
  bwp->rblock_fpos = S_CAST(uint64_t*, arena_alloc_raw_rd(rblock_ct * sizeof(int64_t), &alloc_iter));
  bwp->rrec_len_buf = S_CAST(unsigned char*, arena_alloc_raw_rd(row_ct * bwp->rrec_len_byte_ct, &alloc_iter));
  bwp->rrtype_buf = S_CAST(uintptr_t*, arena_alloc_raw(rrtype_buf_byte_ct, &alloc_iter));
  // the BwAppend... functions assume these bytes are zeroed out
  memset(bwp->rrtype_buf, 0, rrtype_buf_byte_ct);

  const uint32_t col_ct = bwp->col_ct;
  const uint32_t bitvec_byte_alloc = BitCtToCachelineCt(col_ct) * kCacheline;
  bwp->difflist_bitvec_buf = S_CAST(uintptr_t*, arena_alloc_raw(bitvec_byte_alloc, &alloc_iter));
  bwp->prevdiff_base_bitvec = S_CAST(uintptr_t*, arena_alloc_raw(bitvec_byte_alloc, &alloc_iter));
  bwp->fwrite_buf = alloc_iter;
  bwp->fwrite_bufp = alloc_iter;
  // alloc_iter = &(alloc_iter[fwrite_cacheline_ct * kCacheline]);
  return kPglRetSuccess;
}

uint32_t SaveBitvecDifflist(const uintptr_t* bitvec, uint32_t difflist_len, BitmapWriterMain* bwp) {
  unsigned char* fwrite_bufp = bwp->fwrite_bufp;
  if (!difflist_len) {
    *fwrite_bufp = 0;
    bwp->fwrite_bufp = &(fwrite_bufp[1]);
    return 1;
  }
  unsigned char* fwrite_bufp_start = fwrite_bufp;
  fwrite_bufp = Vint32Append(difflist_len, fwrite_bufp);
  const uint32_t col_id_byte_ct = BytesToRepresentNzU32(bwp->col_ct);
  const uint32_t group_ct = DivUp(difflist_len, kPglDifflistGroupSize);
  unsigned char* group_first_sample_ids_iter = fwrite_bufp;
  unsigned char* extra_byte_cts_iter = &(fwrite_bufp[group_ct * col_id_byte_ct]);
  fwrite_bufp = &(extra_byte_cts_iter[group_ct - 1]);
  unsigned char* last_group_vint_start = fwrite_bufp;
  uint32_t last_col_idx = 0;
  uint32_t difflist_idx = 0;
  for (uint32_t widx = 0; ; ++widx) {
    // probable todo: benchmark against SaveLdDifflist()'s strategy on typical
    // workloads.
    uintptr_t xor_word = bitvec[widx];
    if (!xor_word) {
      continue;
    }
    const uint32_t col_idx_base = widx * kBitsPerWord;
    do {
      const uint32_t col_idx_lowbits = ctzw(xor_word);
      const uint32_t new_col_idx = col_idx_base + col_idx_lowbits;
      if (!(difflist_idx % kPglDifflistGroupSize)) {
        group_first_sample_ids_iter = memcpyua(group_first_sample_ids_iter, &new_col_idx, col_id_byte_ct);
        if (difflist_idx) {
          *extra_byte_cts_iter++ = S_CAST(uintptr_t, fwrite_bufp - last_group_vint_start) - (kPglDifflistGroupSize - 1);
        }
        last_group_vint_start = fwrite_bufp;
      } else {
        assert(new_col_idx >= last_col_idx + 1);
        fwrite_bufp = Vint32Append(new_col_idx - last_col_idx, fwrite_bufp);
      }
      ++difflist_idx;
      last_col_idx = new_col_idx;
      if (difflist_idx == difflist_len) {
        bwp->fwrite_bufp = fwrite_bufp;
        return fwrite_bufp - fwrite_bufp_start;
      }
      xor_word &= xor_word - 1;
    } while (xor_word);
  }
}

PglErr BitmapAppend(const uintptr_t* bitvec, BitmapWriter* bw_ptr) {
  BitmapWriterMain* bwp = GetBwp(bw_ptr);

  if (bwp->fwrite_bufp >= &(bwp->fwrite_buf[kPglFwriteBlockSize])) {
    const uintptr_t cur_byte_ct = bwp->fwrite_bufp - bwp->fwrite_buf;
    if (unlikely(fwrite_checked(bwp->fwrite_buf, cur_byte_ct, bwp->ff))) {
      return kPglRetWriteFail;
    }
    bwp->rblock_fpos_offset += cur_byte_ct;
    bwp->fwrite_bufp = bwp->fwrite_buf;
  }

  const uint32_t ridx = bwp->ridx;
  const uint32_t col_ct = bwp->col_ct;
  const uint32_t word_ct = DivUp(col_ct, kBitsPerWord);
  assert((!(col_ct % kBitsPerWord)) || (!(bitvec[word_ct - 1] >> (col_ct % kBitsPerWord))));
  const uint32_t set_bit_ct = PopcountWords(bitvec, word_ct);
  const uint32_t unset_bit_ct = col_ct - set_bit_ct;
  uint32_t difflist_len = set_bit_ct;
  uint32_t common_bit = 0;
  if (set_bit_ct > unset_bit_ct) {
    difflist_len = unset_bit_ct;
    common_bit = 1;
  }
  const uint32_t max_difflist_len = col_ct / kPglMaxBitmapDifflistLenDivisor;
  const uint32_t difflist_viable = (difflist_len <= max_difflist_len);

  uintptr_t* prevdiff_base_bitvec = bwp->prevdiff_base_bitvec;
  const uint32_t prevdiff_base_set_bit_ct = bwp->prevdiff_base_set_bit_ct;
  uintptr_t rrtype;
  uint32_t rrec_len;
  if (!(ridx % kPglRblockSize)) {
    // beginning of a row block.  save raw fpos in header; prevdiff compression
    // prohibited.
    bwp->rblock_fpos[ridx / kPglRblockSize] = bwp->rblock_fpos_offset + S_CAST(uintptr_t, bwp->fwrite_bufp - bwp->fwrite_buf);
  } else {
    const uint32_t col_ctd64 = col_ct / 64;
    if (difflist_len > col_ctd64) {
      // do not use prevdiff compression if there are at least this many
      // differences.  could tune this threshold in the future.
      const uint32_t prevdiff_threshold = difflist_viable? (difflist_len - col_ctd64) : max_difflist_len;
      if (abs_i32(prevdiff_base_set_bit_ct - set_bit_ct) < prevdiff_threshold) {
        // possible todo: benchmark against a single BitvecXorCopyAndCount
        // function on typical workloads, if we don't switch to
        // SaveLdDifflist's strategy of avoiding this array-write
        const uint32_t diff_bit_ct = PopcountWordsXor(prevdiff_base_bitvec, bitvec, word_ct);
        if (diff_bit_ct < prevdiff_threshold) {
          BitvecXorCopy(prevdiff_base_bitvec, bitvec, word_ct, bwp->difflist_bitvec_buf);
          rrec_len = SaveBitvecDifflist(bwp->difflist_bitvec_buf, diff_bit_ct, bwp);
          rrtype = 1;
          goto BitmapAppend_finish;
        }
      }
    }
  }
  {
    memcpy(prevdiff_base_bitvec, bitvec, word_ct * sizeof(intptr_t));
    if (difflist_viable) {
      const uintptr_t* difflist_bitvec = bitvec;
      if (common_bit) {
        BitvecInvertCopy(bitvec, word_ct, bwp->difflist_bitvec_buf);
        difflist_bitvec = bwp->difflist_bitvec_buf;
      }
      rrec_len = SaveBitvecDifflist(difflist_bitvec, difflist_len, bwp);
      rrtype = 2 + common_bit;
    } else {
      rrec_len = DivUp(col_ct, CHAR_BIT);
      bwp->fwrite_bufp = memcpyua(bwp->fwrite_bufp, bitvec, rrec_len);
      rrtype = 0;
    }
  }
 BitmapAppend_finish:
  ;  // needed for pure-C compiler
  const uintptr_t rrec_len_byte_ct = bwp->rrec_len_byte_ct;
  SubU32Store(rrec_len, rrec_len_byte_ct, &(bwp->rrec_len_buf[ridx * rrec_len_byte_ct]));
  bwp->prevdiff_base_set_bit_ct = set_bit_ct;
  bwp->rrtype_buf[ridx / kBitsPerWordD2] |= rrtype << (2 * (ridx % kBitsPerWordD2));
  bwp->ridx = ridx + 1;
  return kPglRetSuccess;
}

PglErr BitmapWriterFinish(BitmapWriter* bw_ptr) {
  BitmapWriterMain* bwp = GetBwp(bw_ptr);
  const uint32_t row_ct = bwp->row_ct;
  assert(bwp->ridx == row_ct);
  FILE* ff = bwp->ff;
  if (unlikely(fwrite_checked(bwp->fwrite_buf, bwp->fwrite_bufp - bwp->fwrite_buf, ff) ||
               fseeko(ff, 11, SEEK_SET))) {
    return kPglRetWriteFail;
  }
  const uint32_t rblock_ct = DivUp(row_ct, kPglRblockSize);
  fwrite_unlocked(bwp->rblock_fpos, rblock_ct * sizeof(int64_t), 1, ff);
  const unsigned char* rrtype_buf_iter = R_CAST(unsigned char*, bwp->rrtype_buf);
  const uint32_t rrec_len_byte_ct = bwp->rrec_len_byte_ct;
  const unsigned char* rrec_len_buf_iter = bwp->rrec_len_buf;
  uint32_t rrec_iter_incr = kPglRblockSize * rrec_len_byte_ct;
  uint32_t rrtype_buf_iter_incr = (kPglRblockSize / 4);
  const unsigned char* rrec_len_buf_last = &(rrec_len_buf_iter[S_CAST(uintptr_t, rblock_ct - 1) * rrec_iter_incr]);
  for (; ; rrec_len_buf_iter = &(rrec_len_buf_iter[rrec_iter_incr])) {
    if (rrec_len_buf_iter >= rrec_len_buf_last) {
      if (rrec_len_buf_iter > rrec_len_buf_last) {
        break;
      }
      const uint32_t rblock_size = ModNz(row_ct, kPglRblockSize);
      rrtype_buf_iter_incr = DivUp(rblock_size, 4);
      rrec_iter_incr = rblock_size * rrec_len_byte_ct;
    }
    fwrite_unlocked(rrtype_buf_iter, rrtype_buf_iter_incr, 1, ff);
    rrtype_buf_iter = &(rrtype_buf_iter[rrtype_buf_iter_incr]);

    if (unlikely(fwrite_checked(rrec_len_buf_iter, rrec_iter_incr, ff))) {
      return kPglRetWriteFail;
    }
  }
  return fclose_null(&(bwp->ff))? kPglRetWriteFail : kPglRetSuccess;
}


BoolErr CleanupBitmapReader(BitmapReader* br_ptr, PglErr* reterrp) {
  BitmapReaderMain* brp = GetBrp(br_ptr);
  if (!brp->ff) {
    return 0;
  }
  if (!fclose_null(&brp->ff)) {
    return 0;
  }
  if (*reterrp != kPglRetSuccess) {
    return 0;
  }
  *reterrp = kPglRetReadFail;
  return 1;
}

BoolErr CleanupBitmapWriter(BitmapWriter* bw_ptr, PglErr* reterrp) {
  BitmapWriterMain* bwp = GetBwp(bw_ptr);
  if (!bwp->ff) {
    return 0;
  }
  if (!fclose_null(&bwp->ff)) {
    return 0;
  }
  if (*reterrp != kPglRetSuccess) {
    return 0;
  }
  *reterrp = kPglRetReadFail;
  return 1;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
