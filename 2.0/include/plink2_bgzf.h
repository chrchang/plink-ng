#ifndef __PLINK2_BGZF_H__
#define __PLINK2_BGZF_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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
#include "../libdeflate/libdeflate.h"

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
CONSTI32(kDecompressChunkSize, 1048576);
static_assert(!(kDecompressChunkSize % kCacheline), "kDecompressChunkSize must be a multiple of kCacheline.");
static_assert(kDecompressChunkSize >= kMaxMediumLine, "kDecompressChunkSize too small.");

CONSTI32(kMaxBgzfDecompressThreads, 5);
CONSTI32(kMaxBgzfCompressedBlockSize, 65536);
static_assert(kMaxBgzfDecompressThreads * 3 * kMaxBgzfCompressedBlockSize < kDecompressChunkSize, "kMaxBgzfDecompressThreads too large relative to kDecompressChunkSize.");
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

void PreinitBgzfRawMtStream(BgzfRawMtDecompressStream* bgzfp);

// Two modes:
// - Regular: ff must point 16 bytes into the file, and header[] must contain
//   the first 16 bytes.  bgzf_st_ptr must be nullptr.
// - Move-construction: ff is assumed to point to &(bgzf_st_ptr->in[in_size]).
//   header must be nullptr.
// decompress_thread_ct must be positive.  It is automatically reduced to
// kMaxBgzfDecompressThreads if necessary.
PglErr BgzfRawMtStreamInit(const char* header, uint32_t decompress_thread_ct, FILE* ff, BgzfRawDecompressStream* bgzf_st_ptr, BgzfRawMtDecompressStream* bgzfp, const char** errmsgp);

PglErr BgzfRawMtStreamRead(unsigned char* dst_end, BgzfRawMtDecompressStream* bgzfp, unsigned char** dst_iterp, const char** errmsgp);

PglErr BgzfRawMtStreamRetarget(const char* header, BgzfRawMtDecompressStream* bgzfp, FILE* next_ff, const char** errmsgp);

HEADER_INLINE PglErr BgzfRawMtStreamRewind(BgzfRawMtDecompressStream* bgzfp, const char** errmsgp) {
  return BgzfRawMtStreamRetarget(nullptr, bgzfp, nullptr, errmsgp);
}

void CleanupBgzfRawMtStream(BgzfRawMtDecompressStream* bgzfp);


// Compression strategy:
// - We have N compression-job memory slots, where N is the smallest power of 2
//   >= 4 * compressor_thread_ct.  (This could be adjusted and/or separately
//   configurable, any value >= 2 * compressor_thread_ct is reasonable.)
// - For the most common plink2 use case (VCF export), compression is over 90%
//   of the compute cost.  Thus, a good implementation must be optimized around
//   minimizing compressor-thread waits, in a setting where one core's workload
//   is split between production and compression and the other cores are
//   compressing all the time.
//   I first implemented the trivial load-balancing strategy (compressor thread
//   k processes blocks k, n+k, 2n+k, ..., where n is the number of compressor
//   threads), but that proved to be slightly slower than htslib's bgzf writer
//   in my testing; it didn't do a good enough job of managing the
//   split-between-production-and-compression core.
//   So I've switched to a more flexible thread pool using __atomic_fetch_add()
//   to distribute compression jobs.
// - There's a dedicated writer thread which flushes the compression results.
// I will look into adding direct support for this parallelization pattern to
// plink2_thread; the benefit is small (~5%), but applies to a lot of
// workloads.  It may also be worthwhile to implement some form of
// work-stealing.
//
// I tried making this less granular (each job contains two blocks instead of
// one) to see if there was still meaningful room for improvement re: reducing
// thread-synchronization overhead; that didn't seem to make any better use of
// the extra memory than keeping the basic granularity level and doubling the
// number of slots.

CONSTI32(kMaxBgzfCompressThreads, 15);
CONSTI32(kBgzfInputBlockSize, 0xff00);  // htslib BGZF_BLOCK_SIZE

typedef struct BgzfCompressStreamMainStruct BgzfCompressStreamMain;

typedef struct BgzfCompressCommWithWStruct {
  // Compressor -> writer.  One per block slot.
  unsigned char cbuf[kMaxBgzfCompressedBlockSize];
  uint32_t nbytes;  // UINT32_MAX = open, otherwise filled
  uint32_t eof;
#ifdef _WIN32
  HANDLE cbuf_filled_event;
  HANDLE cbuf_open_event;
#else
  pthread_mutex_t cbuf_mutex;
  pthread_cond_t cbuf_condvar;
#endif
} BgzfCompressCommWithW;

typedef struct BgzfCompressCommWithPStruct {
  // Producer -> compressor.  One per block slot.
  char ucbuf[kBgzfInputBlockSize];
#ifdef _WIN32
  HANDLE ucbuf_filled_event;
  HANDLE ucbuf_open_event;
#else
  pthread_mutex_t ucbuf_mutex;
  pthread_cond_t ucbuf_condvar;
#endif
  uint32_t nbytes;  // UINT32_MAX = open, otherwise filled
} BgzfCompressCommWithP;

typedef struct BgzfCompressorContextStruct {
  BgzfCompressStreamMain* parent;
  struct libdeflate_compressor* lc;
  // don't need tidx any more...
} BgzfCompressorContext;

struct BgzfCompressStreamMainStruct {
  NONCOPYABLE(BgzfCompressStreamMainStruct);
  FILE* ff;
  pthread_t* threads;  // n+1 elements, last element is writer
  BgzfCompressCommWithP** cwps;  // N elements
  BgzfCompressCommWithW** cwws;  // N elements

  BgzfCompressorContext* compressor_args;  // n elements

  // Atomically read/updated, guaranteed to be on its own cacheline.
  uintptr_t* next_job_idxp;

  int32_t write_errno;

  uint16_t slot_ct;
  uint16_t compressor_thread_ct;  // 0 if no compression
  uint16_t partial_slot_idx;  // this ucbuf must be open
  uint16_t partial_nbytes;  // always < kBgzfInputBlockSize

  // If nonzero, initialization didn't complete, and this value contains enough
  // (platform-specific) information to clean up properly.
  // Not an enum since we perform a bit of arithmetic on it.
  uint16_t unfinished_init_state;
};

typedef struct BgzfCompressStreamStruct {
  NONCOPYABLE(BgzfCompressStreamStruct);
#ifdef __cplusplus
  BgzfCompressStreamMain& GET_PRIVATE_m() { return m; }
  BgzfCompressStreamMain const& GET_PRIVATE_m() const { return m; }
 private:
#endif
  BgzfCompressStreamMain m;
} BgzfCompressStream;

void PreinitBgzfCompressStream(BgzfCompressStream* cstream_ptr);

CONSTI32(kBgzfDefaultClvl, 6);

// - Compression level 0 = plaintext (not a BGZF file at all), so application
//   code no longer has to explicitly branch on bgzf/non-bgzf output.
// - thread_ct only applies when clvl != 0.  When that's true, it's internally
//   clipped to [2, kMaxBgzfCompressThreads].
// - errno is set on open-fail.
PglErr InitBgzfCompressStreamEx(const char* out_fname, uint32_t do_append, uint32_t clvl, uint32_t thread_ct, BgzfCompressStream* cstream_ptr);

HEADER_INLINE PglErr InitBgzfCompressStream(const char* out_fname, uint32_t thread_ct, BgzfCompressStream* cstream_ptr) {
  return InitBgzfCompressStreamEx(out_fname, 0, kBgzfDefaultClvl, thread_ct, cstream_ptr);
}

// errno is set on write-fail.
BoolErr BgzfWrite(const char* buf, uintptr_t len, BgzfCompressStream* cstream_ptr);

BoolErr CleanupBgzfCompressStream(BgzfCompressStream* cstream_ptr, PglErr* reterrp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_BGZF_H__
