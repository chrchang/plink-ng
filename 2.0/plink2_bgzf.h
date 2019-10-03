#ifndef __PLINK2_BGZF_H__
#define __PLINK2_BGZF_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2019 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.


// Binary BGZF readers and writers, based on libdeflate.
//
// htslib dependency has been dropped due to impedance mismatches with plink2
// in several areas:
// - input streams (in particular, being forced to open the file once to detect
//   file type, close it, and then reopen it as a BGZF* breaks pipe file
//   descriptors; and while PLINK 2.0 avoids guaranteeing that pipe file
//   descriptors work in any given scenario, there's still a big practical
//   difference between "usually works, with a few exceptions where two-pass
//   processing is simpler" and "always fails")
// - Windows multithreading
// - error handling
// (Note that htslib has great BGZF compression/decompression performance, as
// long as v1.8+ is used with libdeflate.)

#include "plink2_string.h"
#include "plink2_thread.h"
#include "libdeflate/libdeflate.h"

#ifdef __cplusplus
namespace plink2 {
#endif

HEADER_INLINE int32_t IsBgzfHeader(const void* buf) {
  const uint32_t magic4 = *S_CAST(const uint32_t*, buf);
  return ((magic4 & 0x4ffffff) == 0x4088b1f) && memequal_k(&(S_CAST(const unsigned char*, buf)[10]), "\6\0BC\2", 6);
}

typedef struct BgzfRawDecompressStreamStruct {
  unsigned char* in;
  struct libdeflate_decompressor* ldc;
  uint32_t in_size;
  uint32_t in_pos;
} BgzfRawDecompressStream;

// (tested a few different values for this, 1 MiB appears to work well on the
// systems we care most about)
CONSTI32(kDecompressChunkSizeX, 1048576);
static_assert(!(kDecompressChunkSizeX % kCacheline), "kDecompressChunkSize must be a multiple of kCacheline.");
static_assert(kDecompressChunkSizeX >= kMaxMediumLine, "kDecompressChunkSize too small.");

CONSTI32(kMaxBgzfDecompressThreads, 5);
CONSTI32(kMaxBgzfCompressedBlockSize, 65536);
static_assert(kMaxBgzfDecompressThreads * 3 * kMaxBgzfCompressedBlockSize < kDecompressChunkSizeX, "kMaxBgzfDecompressThreads too large relative to kDecompressChunkSize.");
CONSTI32(kMaxBgzfDecompressedBlockSize, 65536);
static_assert((kMaxBgzfDecompressedBlockSize % kCacheline) == 0, "kMaxBgzfDecompressedBlockSize is assumed to be a cacheline multiple.");

// Communication with reader thread (i.e. worker thread 0).
typedef struct BgzfMtReadCommWithRStruct {
  // Reader -> consumer.  Loaded-but-not-decompressed byte interval.
  uint32_t remaining_start;
  uint32_t remaining_end;
  uint32_t remaining_end_is_eof;
#if __cplusplus >= 201103L
  // Can't just use PglErr, since BgzfRawMtDecompressStream is part of a union,
  // hence BgzfMtReadBody must have a trivial default constructor.
  PglErr::ec reterr;
#else
  int32_t reterr;
#endif
  const char* errmsg;

  // Consumer -> reader.  Currently-being-decompressed byte interval.
  uint32_t locked_start;
  uint32_t locked_end;
} BgzfMtReadCommWithR;

// Communication with decompressor threads.
typedef struct BgzfMtReadCommWithDStruct {
  // Decompressors -> consumer.
  unsigned char* overflow;
  // Consumer is responsible for tracking
  // overflow_start/overflow_end/target_written.

  // This can be set to 1 by any decompressor thread, but no other changes are
  // possible, so we don't need an extra write-lock (or expand this into a
  // kMaxBgzfDecompressThreads-length array).  It takes precedence over
  // reader-thread errors, since it refers to earlier data.
  uint32_t invalid_bgzf;

  // Consumer -> decompressors.
  uint32_t target_capacity;
  unsigned char* target;
  uint32_t in_offsets[kMaxBgzfDecompressThreads + 1];
  uint32_t out_offsets[kMaxBgzfDecompressThreads];
} BgzfMtReadCommWithD;

typedef struct BgzfMtReadBodyStruct {
  // Thread 0 performs read-ahead into trs.in.
  // Threads 1..(thread_ct-1) perform decompression from trs.in directly to
  // target when there's sufficient space, and to overflow[] when there isn't.
  // overflow actually points 64KiB after the start of the associated buffer;
  // this simplifies handling of half-overflowing blocks.
  struct libdeflate_decompressor* ldcs[kMaxBgzfDecompressThreads];

  // Borrowed from consumer, not closed by CleanupBgzfRawMtStream().
  FILE* ff;

  unsigned char* in;

  BgzfMtReadCommWithR* cwr[2];
  BgzfMtReadCommWithD* cwd[2];

  uint32_t initial_compressed_byte_ct;
} BgzfMtReadBody;

typedef struct BgzfRawMtDecompressStreamStruct {
  BgzfMtReadBody body;
  ThreadGroup tg;  // stores thread_ct

  // Between thread-group-joins, the consumer can safely copy from
  // body.cwd[consumer_parity]->overflow, while the decompressor threads
  // continue writing to body.cwd[1 - consumer_parity].
  // overflow_start[consumer_parity]..overflow_end[consumer_parity] indicate
  // which body.cwd[consumer_parity]->overflow bytes remain to be copied, while
  // overflow_start[1 - consumer_parity]..overflow_end[1 - consumer_parity]
  // indicate the initial body.cwd[1 - consumer_parity]->overflow bounds after
  // the next thread-group-join.
  uint32_t overflow_start[2];
  uint32_t overflow_end[2];
  uint32_t consumer_parity;
  uint32_t eof;
} BgzfRawMtDecompressStream;

extern const char kShortErrInvalidBgzf[];

// Two modes:
// - Regular: ff must point 16 bytes into the file, and header[] must contain
//   the first 16 bytes.  bgzf_st_ptr must be nullptr.
// - Move-construction: ff is assumed to point to &(bgzf_st_ptr->in[in_size]).
//   header must be nullptr.
// decompress_thread_ct must be positive.  It is automatically reduced to
// kMaxBgzfDecompressThreads if necessary.
PglErr BgzfRawMtStreamInit(const char* header, uint32_t decompress_thread_ct, FILE* ff, BgzfRawDecompressStream* bgzf_st_ptr, BgzfRawMtDecompressStream* bgzfp, const char** errmsgp);

PglErr BgzfRawMtStreamRead(unsigned char* dst_end, BgzfRawMtDecompressStream* bgzfp, unsigned char** dst_iterp, const char** errmsgp);

PglErr BgzfRawMtStreamRetarget(BgzfRawMtDecompressStream* bgzfp, FILE* next_ff, const char** errmsgp);

HEADER_INLINE PglErr BgzfRawMtStreamRewind(BgzfRawMtDecompressStream* bgzfp, const char** errmsgp) {
  return BgzfRawMtStreamRetarget(bgzfp, nullptr, errmsgp);
}

void CleanupBgzfRawMtStream(BgzfRawMtDecompressStream* bgzfp);


#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_BGZF_H__
