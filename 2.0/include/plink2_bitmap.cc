// This library is part of PLINK 2.00, copyright (C) 2005-2021 Shaun Purcell,
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

/*
static inline const uint32_t* GetCicp(PgrColSubsetIndex pcsi) {
  return GET_PRIVATE(pcsi, cumulative_popcounts);
}
*/

void PreinitBitmapReader(BitmapReader* brp) {
  brp->ff = nullptr;
}

void PreinitBitmapWriter(BitmapWriter* bwp) {
  bwp->ff = nullptr;
}

uintptr_t CountBrAllocCachelinesRequired(uint32_t row_ct, uint32_t col_ct) {
  // row_fpos: sizeof(int64_t) * (row_ct + 1) bytes
  // rrtype_nyparr: DivUp(row_ct + 1, 4) bytes
  // fread_buf, prevdiff_base_raw_bitvec, prevdiff_base_bitvec:
  //   DivUp(col_ct, CHAR_BIT) bytes each
  // prevdiff_list_col_ids: DivUp(col_ct, kPglMaxBitmapDifflistLenDivisor) *
  //                        sizeof(int32_t)

  uintptr_t cachelines_required = DivUp(row_ct + 1, kInt64PerCacheline);
  cachelines_required += DivUp(row_ct + 1, kNypsPerCacheline);
  cachelines_required += 3 * DivUp(col_ct, kBitsPerCacheline);
  cachelines_required += DivUp(col_ct, kPglMaxBitmapDifflistLenDivisor * kInt32PerCacheline);
  return cachelines_required;
}

PglErr BitmapReaderInitPhase1(const char* fname, BitmapReader* brp, uintptr_t* br_alloc_cacheline_ct_ptr, char* errstr_buf) {
  // Open the bitmap; verify that the initial bytes are consistent with the
  // file format; load row/column counts; determine initial memory allocation
  // requirement.
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

PglErr BitmapReaderInitPhase2(BitmapReader* brp, unsigned char* br_alloc, char* errstr_buf) {
  const uintptr_t row_ct = brp->row_ct;
  const uintptr_t col_ct = brp->col_ct;
  // Allocate data structures from the br_alloc workspace.
  unsigned char* br_alloc_iter = br_alloc;
  brp->row_fpos = S_CAST(uint64_t*, arena_alloc_raw_rd((row_ct + 1) * sizeof(int64_t), &br_alloc_iter));
  brp->rrtype_nyparr = S_CAST(uintptr_t*, arena_alloc_raw((1 + (row_ct / kNypsPerCacheline)) * kCacheline, &br_alloc_iter));
  const uintptr_t bitvec_byte_ct = DivUp(col_ct, kBitsPerCacheline) * kCacheline;
  brp->fread_buf = S_CAST(unsigned char*, arena_alloc_raw(bitvec_byte_ct, &br_alloc_iter));
  brp->prevdiff_base_raw_bitvec = S_CAST(uintptr_t*, arena_alloc_raw(bitvec_byte_ct, &br_alloc_iter));
  brp->prevdiff_base_bitvec = S_CAST(uintptr_t*, arena_alloc_raw(bitvec_byte_ct, &br_alloc_iter));
  brp->prevdiff_base_col_ids = S_CAST(uint32_t*, arena_alloc_raw_rd(DivUp(col_ct, kPglMaxBitmapDifflistLenDivisor) * sizeof(int32_t), &br_alloc_iter));

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
    if (unlikely(fread_unlocked(rrtype_nyparr_iter, cur_byte_ct, 1, ff))) {
      FillBrReadErrstr(ff, errstr_buf);
      return kPglRetReadFail;
    }
    rrtype_nyparr_iter = &(rrtype_nyparr_iter[cur_rblock_row_ct / kBytesPerWord]);

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

/*
PglErr BitmapGet(const uintptr_t* __restrict col_include, PgrColSubsetIndex pcsi, uint32_t col_ct, uint32_t ridx, BitmapReader* brp, uintptr_t* __restrict bitvec) {
  if (!ridx) {
    return kPglRetSuccess;
  }
  assert(ridx < brp->row_ct);
  const uint32_t* col_include_cumulative_popcounts = GetCicp(pcsi);

  const uint32_t rrtype = GetRrtype(brp, ridx);
  if (RrtypeIsPrevdiff(rrtype)) {
    // PglErr reterr = PrevdiffLoadAndCopyBitvecSubset(col_include, col_include_cumulative_popcounts, );

  }
  return kPglRetSuccess;
}
*/

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

  // prevdiff_base_bitvec
  cachelines_required += 3 * BitCtToCachelineCt(col_ct);

  const uint32_t max_difflist_len = 2 * (col_ct / kPglMaxBitmapDifflistLenDivisor);

  // prevdiff_base_col_ids
  cachelines_required += 1 + (max_difflist_len / kInt32PerCacheline);

  // fwrite_buf
  // TODO: determine appropriate overallocation to simplify difflist writing
  cachelines_required += DivUp(max_rrec_len + kPglFwriteBlockSize, kCacheline);
  return cachelines_required;
}

PglErr BitmapWriterInitPhase1(const char* fname, uint32_t row_ct, uint32_t col_ct, BitmapWriter* bwp, uintptr_t* bw_alloc_cacheline_ct_ptr) {
  assert(row_ct);
  assert(col_ct);

  bwp->row_ct = row_ct;
  bwp->col_ct = col_ct;
#ifndef NDEBUG
  bwp->rblock_fpos = nullptr;
  bwp->rrec_len_buf = nullptr;
  bwp->rrtype_buf = nullptr;
  bwp->prevdiff_base_bitvec = nullptr;
  bwp->prevdiff_base_col_ids = nullptr;
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

PglErr BitmapWriterInitPhase2(BitmapWriter* bwp, unsigned char* bw_alloc) {
  const uint32_t row_ct = bwp->row_ct;
  unsigned char* alloc_iter = bw_alloc;
  const uint32_t rblock_ct = DivUp(row_ct, kPglRblockSize);
  const uint32_t rrtype_buf_byte_ct = BitCtToCachelineCt(row_ct) * kCacheline;
  bwp->rblock_fpos = S_CAST(uint64_t*, arena_alloc_raw_rd(rblock_ct * sizeof(int64_t), &alloc_iter));
  bwp->rrec_len_buf = S_CAST(unsigned char*, arena_alloc_raw_rd(row_ct * bwp->rrec_len_byte_ct, &alloc_iter));
  bwp->rrtype_buf = S_CAST(uintptr_t*, arena_alloc_raw(rrtype_buf_byte_ct, &alloc_iter));
  // the BwAppend... functions assume these bytes are zeroed out
  memset(bwp->rrtype_buf, 0, rrtype_buf_byte_ct);

  const uint32_t col_ct = bwp->col_ct;
  const uint32_t bitvec_byte_alloc = BitCtToCachelineCt(col_ct) * kCacheline;
  const uint32_t max_difflist_len = 2 * (col_ct / kPglMaxBitmapDifflistLenDivisor);
  bwp->prevdiff_base_bitvec = S_CAST(uintptr_t*, arena_alloc_raw(bitvec_byte_alloc, &alloc_iter));
  bwp->prevdiff_base_col_ids = S_CAST(uint32_t*, arena_alloc_raw_rd(max_difflist_len * sizeof(int32_t), &alloc_iter));
  bwp->fwrite_buf = alloc_iter;
  bwp->fwrite_bufp = alloc_iter;
  // alloc_iter = &(alloc_iter[fwrite_cacheline_ct * kCacheline]);
  return kPglRetSuccess;
}

BoolErr CleanupBitmapReader(BitmapReader* brp, PglErr* reterrp) {
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

BoolErr CleanupBitmapWriter(BitmapWriter* bwp, PglErr* reterrp) {
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
