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

#include <errno.h>
#include "plink2_bgzf.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static inline BgzfCompressStreamMain* GetBgzfp(BgzfCompressStream* cstream_ptr) {
  return &GET_PRIVATE(*cstream_ptr, m);
}

void PreinitBgzfRawMtStream(BgzfRawMtDecompressStream* bgzfp) {
  PreinitThreads(&bgzfp->tg);
  bgzfp->body.in = nullptr;
  // Main unconditional allocation starts from cfr[0]; cfr[1], twc[0], twc[1],
  // and the twc[] overflow buffers are fixed offsets within that allocation,
  // so we only have to cache-align once.
  bgzfp->body.cwr[0] = nullptr;
}

const char kShortErrInvalidBgzf[] = "Malformed BGZF block";

CONSTI32(kBgzfRawMtStreamRetargetCode, 0x7fffffff);
static_assert(kBgzfRawMtStreamRetargetCode > kMaxBgzfCompressedBlockSize, "kBgzfRawMtStreamRetargetCode must be outside the valid locked_start range.");

CONSTI32(kBgzfRawMtStreamMaxCapacity, 0x7fffffc0);

// This function usually repeatedly joins and respawns the reader and
// decompressor worker threads, until *dst_iterp reaches dst_end or we've hit
// EOF.  In the EOF case, *dst_iterp is set to 1 past the last read byte;
// otherwise it's set to dst_end (assuming no error).
// There's one special case: if dst_end and dst_iterp are nullptr, this
// performs exactly one join-and-respawn (relevant during initialization and
// rewind).
//
// Preconditions:
// - bgzfp->overflow_start[bgzfp->consumer_parity] ==
//   bgzfp->overflow_end[bgzfp->consumer_parity], i.e. all previously-dumped
//   overflow bytes were previously appended to dst_iter.
// - In the non-null case, *dst_iterp points to where
//   bgzfp->body.cwr[1 - bgzfp->consumer_parity]->overflow[0] should be
//   appended.
// - Thread-group is unjoined, and bgzfp->eof must be false.
PglErr BgzfReadJoinAndRespawn(unsigned char* dst_end, BgzfRawMtDecompressStream* bgzfp, unsigned char** dst_iterp, const char** errmsgp) {
  PglErr reterr = kPglRetSuccess;
  ThreadGroup* tgp = &bgzfp->tg;
  unsigned char* dst_iter = nullptr;
  if (dst_iterp) {
    dst_iter = *dst_iterp;
  }
  unsigned char* next_target;
  do {
    JoinThreads(tgp);

    // 1. Check for decompression and read errors.
    const uint32_t next_producer_parity = bgzfp->consumer_parity;
    const uint32_t prev_producer_parity = 1 - next_producer_parity;
    BgzfMtReadBody* bodyp = &bgzfp->body;
    BgzfMtReadCommWithD* cwd = bodyp->cwd[prev_producer_parity];
    if (unlikely(cwd->invalid_bgzf)) {
      goto BgzfReadJoinAndRespawn_ret_INVALID_BGZF;
    }
    BgzfMtReadCommWithR* cwr = bodyp->cwr[prev_producer_parity];
    if (unlikely(cwr->reterr != kPglRetSuccess)) {
      *errmsgp = cwr->errmsg;
      reterr = S_CAST(PglErr, cwr->reterr);
      goto BgzfReadJoinAndRespawn_ret_1;
    }

    // 2. Determine amount of remaining dst space after existing-overflow-copy.
    const uint32_t remaining_start = cwr->remaining_start;
    const uint32_t remaining_end = cwr->remaining_end;
    const uint32_t remaining_end_is_eof = cwr->remaining_end_is_eof;
    const uint32_t new_overflow_ct = bgzfp->overflow_end[prev_producer_parity];
    unsigned char* overflow_dst_start = dst_iter;
    next_target = nullptr;
    uint32_t overflow_copy_ct = 0;
    uint32_t target_capacity = 0;
    if (dst_iter) {
      uintptr_t dst_capacity = dst_end - dst_iter;
      if (dst_capacity > new_overflow_ct) {
        overflow_copy_ct = new_overflow_ct;
        next_target = &(dst_iter[new_overflow_ct]);
        dst_capacity -= new_overflow_ct;
        if (dst_capacity <= kBgzfRawMtStreamMaxCapacity) {
          target_capacity = dst_capacity;
        } else {
          // May as well clip to int32.
          target_capacity = kBgzfRawMtStreamMaxCapacity;
        }
      } else {
        overflow_copy_ct = dst_capacity;
      }
    }
    // 3. Skip thread relaunch on eof.
    if (remaining_start == remaining_end) {
      assert(remaining_end_is_eof);
      bgzfp->eof = 1;
      // bugfix (2 Oct 2019): new_overflow_ct -> overflow_copy_ct
      dst_iter = &(dst_iter[overflow_copy_ct]);
      next_target = nullptr;
    } else {
      // 4. Determine block boundaries in in[remaining_start, remaining_end),
      //    and mark a multiple-of-decompress_thread_ct blocks for concurrent
      //    decompression, taking advantage of any remaining dst space.  (If
      //    there's no dst space left, the overflow buffer still lets us
      //    decompress the next decompress_thread_ct blocks in the background.)
      const uint32_t decompress_thread_ct = GetThreadCtTg(tgp) - 1;
      // We actually iterate through the blocks twice.  The first iteration
      // counts the number of blocks that target_capacity allows for, and the
      // second iteration fills next_cwd->{in_offsets, out_offsets}.
      unsigned char* in = bodyp->in;
      uint32_t n_blocks_per_thread = 0;
      unsigned char* in_iter = &(in[remaining_start]);
      unsigned char* in_end = &(in[remaining_end]);
      uint32_t write_offset = 0;
      uint32_t n_eof_blocks = 0;
      while (write_offset <= target_capacity) {
        uint32_t uii = 0;
        for (; uii != decompress_thread_ct; ++uii) {
          const uint32_t n_inbytes = in_end - in_iter;
          if (n_inbytes <= 25) {
            if (unlikely(remaining_end_is_eof && n_inbytes)) {
              goto BgzfReadJoinAndRespawn_ret_INVALID_BGZF;
            }
            break;
          }
          if (unlikely(!IsBgzfHeader(in_iter))) {
            goto BgzfReadJoinAndRespawn_ret_INVALID_BGZF;
          }
#  ifdef __arm__
#    error "Unaligned accesses in BgzfReadJoinAndRespawn()."
#  endif
          const uint32_t bsize_minus1 = *R_CAST(uint16_t*, &(in_iter[16]));
          if (unlikely(bsize_minus1 < 25)) {
            goto BgzfReadJoinAndRespawn_ret_INVALID_BGZF;
          }
          if (bsize_minus1 >= n_inbytes) {
            if (unlikely(remaining_end_is_eof)) {
              goto BgzfReadJoinAndRespawn_ret_INVALID_BGZF;
            }
            break;
          }
          const uint32_t in_size = bsize_minus1 - 25;
          const uint32_t out_size = *R_CAST(uint32_t*, &(in_iter[in_size + 22]));
          if (unlikely(out_size > 65536)) {
            goto BgzfReadJoinAndRespawn_ret_INVALID_BGZF;
          }
          in_iter = &(in_iter[bsize_minus1 + 1]);
          write_offset += out_size;
        }
        if (uii != decompress_thread_ct) {
          if (remaining_end_is_eof && (in_iter == in_end)) {
            n_eof_blocks = uii;
          }
          break;
        }
        ++n_blocks_per_thread;
      }

      // Second iteration.
      in_iter = &(in[remaining_start]);
      write_offset = 0;
      BgzfMtReadCommWithD* next_cwd = bodyp->cwd[next_producer_parity];
      next_cwd->target = next_target;
      next_cwd->target_capacity = target_capacity;
      uint32_t* in_offsets = next_cwd->in_offsets;
      uint32_t* out_offsets = next_cwd->out_offsets;
      for (uint32_t out_tidx = 0; out_tidx != decompress_thread_ct; ++out_tidx) {
        in_offsets[out_tidx] = in_iter - in;
        out_offsets[out_tidx] = write_offset;
        const uint32_t nblocks = n_blocks_per_thread + (n_eof_blocks > out_tidx);
        for (uint32_t uii = 0; uii != nblocks; ++uii) {
          const uint32_t bsize_minus1 = *R_CAST(uint16_t*, &(in_iter[16]));
          const uint32_t in_size = bsize_minus1 - 25;
          const uint32_t out_size = *R_CAST(uint32_t*, &(in_iter[in_size + 22]));
          in_iter = &(in_iter[bsize_minus1 + 1]);
          write_offset += out_size;
        }
      }
      const uint32_t locked_end = in_iter - in;
      in_offsets[decompress_thread_ct] = locked_end;

      BgzfMtReadCommWithR* next_cwr = bodyp->cwr[next_producer_parity];
      next_cwr->locked_start = remaining_start;
      next_cwr->locked_end = locked_end;
      SpawnThreads(tgp);

      bgzfp->overflow_start[next_producer_parity] = 0;
      uint32_t next_overflow_end = 0;
      if (write_offset < target_capacity) {
        dst_iter = &(next_target[write_offset]);
      } else {
        next_overflow_end = write_offset - target_capacity;
        dst_iter = dst_end;
      }
      bgzfp->overflow_end[next_producer_parity] = next_overflow_end;
    }
    bgzfp->consumer_parity = prev_producer_parity;
    if (overflow_copy_ct) {
      // Critical for this to happen after SpawnThreads() call.
      memcpy(overflow_dst_start, cwd->overflow, overflow_copy_ct);
      bgzfp->overflow_start[prev_producer_parity] = overflow_copy_ct;
    }
  } while (next_target != nullptr);
  if (dst_iterp) {
    *dst_iterp = dst_iter;
  }
  while (0) {
  BgzfReadJoinAndRespawn_ret_INVALID_BGZF:
    *errmsgp = kShortErrInvalidBgzf;
    reterr = kPglRetDecompressFail;
    break;
  }
 BgzfReadJoinAndRespawn_ret_1:
  return reterr;
}

THREAD_FUNC_DECL BgzfRawMtStreamThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  BgzfRawMtDecompressStream* context = S_CAST(BgzfRawMtDecompressStream*, arg->sharedp->context);
  BgzfMtReadBody* bodyp = &context->body;
  unsigned char* in = bodyp->in;
  uint32_t tidx = arg->tidx;
  uint32_t parity = 0;
  if (!tidx) {
    // Thread 0 reads raw compressed bytes into bodyp->in.  This only uses a
    // small fraction of a processor core, but fread has enough latency
    // (especially when the input file isn't in cache) that we don't want the
    // consumer thread to block on it.
    FILE* ff = bodyp->ff;
    // in[] has space for kDecompressChunkSize bytes.
    // The current buffer-usage logic assumes that thresh1 <= 1/3 of the buffer
    // size, and thresh2 >= 2/3.
    const uint32_t thresh1 = (GetThreadCtTg(&context->tg) - 1) * (26 + kMaxBgzfCompressedBlockSize);
    const uint32_t thresh2 = kDecompressChunkSize - thresh1;
    uint32_t remaining_read_start = bodyp->initial_compressed_byte_ct;
    uint32_t is_eof = 0;
    do {
      BgzfMtReadCommWithR* cwr = bodyp->cwr[parity];
      uint32_t locked_start = cwr->locked_start;
      uint32_t locked_end = cwr->locked_end;
      // in[locked_start, locked_end) is being decompressed by the other
      // threads, and must not be touched on this iteration.
      //
      // We have the following three regular states:
      // 1. locked_start <= locked_end < thresh1 (always true on function
      //    entry)  Try to load up to &(in[thresh2]).
      // 2. locked_start < thresh1 <= locked_end <= thresh2.  Try to load up to
      //    &(in[kDecompressChunkSize]) (the end of the buffer).
      // 3. thresh1 <= locked_start <= thresh2 < locked_end.  Copy
      //    [locked_end, kDecompressChunkSize) back to the beginning of the
      //    buffer, and try to load up to &(in[prev_start_offset]).
      // We also have the following special case:
      // 4. locked_start == kBgzfRawMtStreamRetargetCode indicates that the
      //    consumer has rewound or retargeted ff.
      if (locked_start == kBgzfRawMtStreamRetargetCode) {
        ff = bodyp->ff;
        locked_start = 0;
        // Consumer is expected to set cwr->locked_end = 0 here.
        remaining_read_start = 16;  // bugfix
        is_eof = 0;
      }
      uint32_t remaining_end;
      if (locked_end < thresh1) {
        // state 1
        remaining_end = thresh2;
      } else if (locked_end <= thresh2) {
        // state 2
        remaining_end = kDecompressChunkSize;
      } else {
        // state 3
        remaining_read_start -= locked_end;
        memcpy(in, &(in[locked_end]), remaining_read_start);
        locked_end = 0;
        remaining_end = locked_start;
      }
      uint32_t nbytes = 0;
      if (remaining_end > remaining_read_start) {
        if (!is_eof) {
          nbytes = fread_unlocked(&(in[remaining_read_start]), 1, remaining_end - remaining_read_start, ff);
          if (unlikely(ferror_unlocked(ff))) {
            cwr->errmsg = strerror(errno);
            cwr->reterr = kPglRetReadFail;
            continue;
          }
          is_eof = feof_unlocked(ff);
        }
        remaining_end = remaining_read_start + nbytes;
      }
      cwr->remaining_start = locked_end;
      cwr->remaining_end = remaining_end;
      cwr->remaining_end_is_eof = is_eof;
      remaining_read_start = remaining_end;
      parity = 1 - parity;
    } while (!THREAD_BLOCK_FINISH(arg));
  } else {
    // Threads 1..decompress_thread_ct decompress from bodyp->in to
    // bodyp->target, storing overflow bytes in the cwd->overflow buffers.
    --tidx;  // 0-based decompressor-thread indexes are more convenient here.
    struct libdeflate_decompressor* ldc = bodyp->ldcs[tidx];

    // Note that we do nothing on the first iteration: in_offsets[] is
    // initialized to all-zero, while we wait for enough raw bytes to load.
    do {
      BgzfMtReadCommWithD* cwd = bodyp->cwd[parity];
      unsigned char* overflow = cwd->overflow;
      unsigned char* target = cwd->target;
      uint32_t in_offset = cwd->in_offsets[tidx];
      const uint32_t in_offset_stop = cwd->in_offsets[tidx + 1];
      uint32_t out_offset = cwd->out_offsets[tidx];
      const uint32_t out_capacity = cwd->target_capacity;
      while (in_offset != in_offset_stop) {
        const uint32_t in_size = (*R_CAST(uint16_t*, &(in[in_offset + 16]))) - 25;
        const uint32_t out_size = *R_CAST(uint32_t*, &(in[in_offset + in_size + 22]));
        const uint32_t out_offset_end = out_offset + out_size;
        unsigned char* dst;
        if (out_offset_end > out_capacity) {
          dst = &(overflow[S_CAST(int32_t, out_offset - out_capacity)]);
        } else {
          dst = &(target[out_offset]);
        }
        if (unlikely(libdeflate_deflate_decompress(ldc, &(in[in_offset + 18]), in_size, dst, out_size, nullptr))) {
          cwd->invalid_bgzf = 1;
          break;
        }
        if ((out_offset_end > out_capacity) && (out_offset < out_capacity)) {
          memcpy(&(target[out_offset]), dst, out_capacity - out_offset);
        }
        in_offset += in_size + 26;
        out_offset = out_offset_end;
      }
      parity = 1 - parity;
    } while (!THREAD_BLOCK_FINISH(arg));
  }
  THREAD_RETURN;
}

PglErr BgzfRawMtStreamInit(const char* header, uint32_t decompress_thread_ct, FILE* ff, BgzfRawDecompressStream* bgzf_st_ptr, BgzfRawMtDecompressStream* bgzfp, const char** errmsgp) {
  PglErr reterr = kPglRetSuccess;
  {
    PreinitBgzfRawMtStream(bgzfp);
    assert(decompress_thread_ct);
    if (decompress_thread_ct > kMaxBgzfDecompressThreads) {
      decompress_thread_ct = kMaxBgzfDecompressThreads;
    }
    BgzfMtReadBody* bodyp = &bgzfp->body;
    ZeroPtrArr(decompress_thread_ct, bodyp->ldcs);
    ThreadGroup* tgp = &bgzfp->tg;
    if (unlikely(SetThreadCt(decompress_thread_ct + 1, tgp))) {
      // May as well avoid these leaks.
      if (bgzf_st_ptr) {
        free(bgzf_st_ptr->in);
        libdeflate_free_decompressor(bgzf_st_ptr->ldc);
      }
      goto BgzfRawMtStreamInit_ret_NOMEM;
    }
    if (bgzf_st_ptr) {
      bodyp->in = bgzf_st_ptr->in;
      const uint32_t in_pos = bgzf_st_ptr->in_pos;
      const uint32_t initial_compressed_byte_ct = bgzf_st_ptr->in_size - in_pos;
      memmove(bodyp->in, &(bodyp->in[in_pos]), initial_compressed_byte_ct);
      bodyp->initial_compressed_byte_ct = initial_compressed_byte_ct;
    } else {
      bodyp->in = S_CAST(unsigned char*, malloc(kDecompressChunkSize));
      if (unlikely(!bodyp->in)) {
        goto BgzfRawMtStreamInit_ret_NOMEM;
      }
      memcpy(bodyp->in, header, 16);
      bodyp->initial_compressed_byte_ct = 16;
    }
    for (uint32_t tidx = 0; tidx < decompress_thread_ct; ++tidx) {
      if ((!tidx) && bgzf_st_ptr) {
        bodyp->ldcs[0] = bgzf_st_ptr->ldc;
      } else {
        bodyp->ldcs[tidx] = libdeflate_alloc_decompressor();
        if (!bodyp->ldcs[tidx]) {
          goto BgzfRawMtStreamInit_ret_NOMEM;
        }
      }
    }
    assert(!bodyp->cwr[0]);
    const uint32_t cwr_aligned_size = RoundUpPow2(sizeof(BgzfMtReadCommWithR), kCacheline);
    const uint32_t cwd_aligned_size = RoundUpPow2(sizeof(BgzfMtReadCommWithD), kCacheline);
    const uint32_t overflow_buf_size = kMaxBgzfDecompressedBlockSize * (decompress_thread_ct + 1);
    unsigned char* raw_alloc;
    if (unlikely(cachealigned_malloc(2 * (cwr_aligned_size + cwd_aligned_size + overflow_buf_size), &raw_alloc))) {
      goto BgzfRawMtStreamInit_ret_NOMEM;
    }
    bodyp->cwr[0] = R_CAST(BgzfMtReadCommWithR*, raw_alloc);
    raw_alloc = &(raw_alloc[cwr_aligned_size]);
    bodyp->cwr[1] = R_CAST(BgzfMtReadCommWithR*, raw_alloc);
    raw_alloc = &(raw_alloc[cwr_aligned_size]);
    bodyp->cwd[0] = R_CAST(BgzfMtReadCommWithD*, raw_alloc);
    raw_alloc = &(raw_alloc[cwd_aligned_size]);
    bodyp->cwd[1] = R_CAST(BgzfMtReadCommWithD*, raw_alloc);
    raw_alloc = &(raw_alloc[cwd_aligned_size]);
    bodyp->cwd[0]->overflow = &(raw_alloc[kMaxBgzfDecompressedBlockSize]);
    raw_alloc = &(raw_alloc[overflow_buf_size]);
    bodyp->cwd[1]->overflow = &(raw_alloc[kMaxBgzfDecompressedBlockSize]);
    // raw_alloc = &(raw_alloc[overflow_buf_size]);

    bodyp->ff = ff;
    for (uint32_t parity = 0; parity != 2; ++parity) {
      bodyp->cwr[parity]->errmsg = nullptr;
      bodyp->cwr[parity]->reterr = kPglRetSuccess;
      bodyp->cwr[parity]->locked_start = 0;
      bodyp->cwr[parity]->locked_end = 0;
      bodyp->cwd[parity]->invalid_bgzf = 0;
      bodyp->cwd[parity]->target_capacity = 0;
      bodyp->cwd[parity]->target = nullptr;
      ZeroU32Arr(kMaxBgzfDecompressThreads + 1, bodyp->cwd[parity]->in_offsets);
      // out_offsets doesn't matter when in_offsets is zeroed
    }

    SetThreadFuncAndData(BgzfRawMtStreamThread, bgzfp, tgp);
    if (unlikely(SpawnThreads(tgp))) {
      goto BgzfRawMtStreamInit_ret_THREAD_CREATE_FAIL;
    }
    bgzfp->overflow_start[0] = 0;
    bgzfp->overflow_start[1] = 0;
    bgzfp->overflow_end[0] = 0;
    bgzfp->overflow_end[1] = 0;
    bgzfp->consumer_parity = 1;
    bgzfp->eof = 0;

    // Bottleneck is usually decompression, so we want it to be happening in
    // the background when BgzfRawMtStreamInit() returns.  This
    // currently requires a join-and-respawn.
    reterr = BgzfReadJoinAndRespawn(nullptr, bgzfp, nullptr, errmsgp);
  }
  while (0) {
  BgzfRawMtStreamInit_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  BgzfRawMtStreamInit_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  return reterr;
}

PglErr BgzfRawMtStreamRead(unsigned char* dst_end, BgzfRawMtDecompressStream* bgzfp, unsigned char** dst_iterp, const char** errmsgp) {
  const uint32_t consumer_parity = bgzfp->consumer_parity;
  const uint32_t overflow_start = bgzfp->overflow_start[consumer_parity];
  const uint32_t overflow_end = bgzfp->overflow_end[consumer_parity];
  const uint32_t overflow_remaining = overflow_end - overflow_start;
  const uintptr_t dst_capacity = dst_end - (*dst_iterp);
  unsigned char* overflow_src_start = &(bgzfp->body.cwd[consumer_parity]->overflow[overflow_start]);
  if (overflow_remaining >= dst_capacity) {
    memcpy(*dst_iterp, overflow_src_start, dst_capacity);
    bgzfp->overflow_start[consumer_parity] += dst_capacity;
    *dst_iterp = dst_end;
    return kPglRetSuccess;
  }
  bgzfp->overflow_start[consumer_parity] += overflow_remaining;
  *dst_iterp = memcpyua(*dst_iterp, overflow_src_start, overflow_remaining);
  if (bgzfp->eof) {
    return kPglRetSuccess;
  }
  return BgzfReadJoinAndRespawn(dst_end, bgzfp, dst_iterp, errmsgp);
}

PglErr BgzfRawMtStreamRetarget(const char* header, BgzfRawMtDecompressStream* bgzfp, FILE* next_ff, const char** errmsgp) {
  BgzfMtReadBody* bodyp = &bgzfp->body;
  ThreadGroup* tgp = &bgzfp->tg;
  if (!bgzfp->eof) {
    JoinThreads(tgp);
    const uint32_t prev_producer_parity = 1 - bgzfp->consumer_parity;
    BgzfMtReadCommWithD* prev_cwd = bodyp->cwd[prev_producer_parity];
    if (unlikely(prev_cwd->invalid_bgzf)) {
      *errmsgp = kShortErrInvalidBgzf;
      return kPglRetDecompressFail;
    }
    BgzfMtReadCommWithR* prev_cwr = bodyp->cwr[prev_producer_parity];
    if (unlikely(prev_cwr->reterr != kPglRetSuccess)) {
      *errmsgp = prev_cwr->errmsg;
      return S_CAST(PglErr, prev_cwr->reterr);
    }
    bgzfp->consumer_parity = prev_producer_parity;
  }
  const uint32_t next_producer_parity = 1 - bgzfp->consumer_parity;
  for (uint32_t parity = 0; parity != 2; ++parity) {
    bodyp->cwr[parity]->locked_start = 0;
    bodyp->cwr[parity]->locked_end = 0;
    bodyp->cwd[parity]->target_capacity = 0;
    bodyp->cwd[parity]->target = nullptr;
    ZeroU32Arr(kMaxBgzfDecompressThreads + 1, bodyp->cwd[parity]->in_offsets);
    // bugfix (3 Oct 2019): forgot this
    bgzfp->overflow_start[parity] = 0;
    bgzfp->overflow_end[parity] = 0;
  }
  BgzfMtReadCommWithR* next_cwr = bodyp->cwr[next_producer_parity];
  next_cwr->locked_start = kBgzfRawMtStreamRetargetCode;
  if (next_ff == nullptr) {
    rewind(bodyp->ff);
    // bugfix (8 Feb 2020): need to explicitly read the first 16 bytes.
    if (unlikely(!fread_unlocked(bodyp->in, 16, 1, bodyp->ff))) {
      return kPglRetRewindFail;
    }
  } else {
    // Caller is responsible for closing previous bodyp->ff, etc.
    bodyp->ff = next_ff;
    // bugfix (5 Oct 2019): forgot this
    memcpy(bodyp->in, header, 16);
  }
  SpawnThreads(tgp);
  bgzfp->eof = 0;
  // Turn the crank once, for the same reason we do so during stream creation.
  return BgzfReadJoinAndRespawn(nullptr, bgzfp, nullptr, errmsgp);
}

void CleanupBgzfRawMtStream(BgzfRawMtDecompressStream* bgzfp) {
  uint32_t decompress_thread_ct = GetThreadCtTg(&bgzfp->tg);
  if (decompress_thread_ct) {
    --decompress_thread_ct;
  }
  CleanupThreads(&bgzfp->tg);
  BgzfMtReadBody* bodyp = &bgzfp->body;
  for (uint32_t tidx = 0; tidx < decompress_thread_ct; ++tidx) {
    if (bodyp->ldcs[tidx]) {
      libdeflate_free_decompressor(bodyp->ldcs[tidx]);
    }
  }
  if (bodyp->in) {
    free(bodyp->in);
    bodyp->in = nullptr;
  }
  if (bodyp->cwr[0]) {
    aligned_free(bodyp->cwr[0]);
    bodyp->cwr[0] = nullptr;
  }
}


void PreinitBgzfCompressStream(BgzfCompressStream* cstream_ptr) {
  BgzfCompressStreamMain* bgzfp = GetBgzfp(cstream_ptr);
  bgzfp->ff = nullptr;
  bgzfp->threads = nullptr;
}

static const unsigned char kBgzfEofBlock[] = "\37\213\10\4\0\0\0\0\0\377\6\0\102\103\2\0\33\0\3\0\0\0\0\0\0\0\0";

THREAD_FUNC_DECL BgzfCompressorThread(void* raw_arg) {
  BgzfCompressorContext* context = S_CAST(BgzfCompressorContext*, raw_arg);
  struct libdeflate_compressor* compressor = context->lc;
  BgzfCompressStreamMain* bgzfp = context->parent;
  const uint32_t slot_mask = bgzfp->slot_ct - 1;
  BgzfCompressCommWithP** cwps = bgzfp->cwps;
  BgzfCompressCommWithW** cwws = bgzfp->cwws;
  uintptr_t* next_job_idxp = bgzfp->next_job_idxp;
  while (1) {
    const uint32_t slot_idx = __atomic_fetch_add(next_job_idxp, 1, __ATOMIC_ACQ_REL) & slot_mask;
    BgzfCompressCommWithP* cwp = cwps[slot_idx];
    BgzfCompressCommWithW* cww = cwws[slot_idx];
#ifdef _WIN32
    HANDLE wait_pair[2];
    wait_pair[0] = cwp->ucbuf_filled_event;
    wait_pair[1] = cww->cbuf_open_event;
    WaitForMultipleObjects(2, wait_pair, 1, INFINITE);
#else
    pthread_mutex_lock(&(cww->cbuf_mutex));
    while (cww->nbytes != UINT32_MAX) {
      pthread_cond_wait(&(cww->cbuf_condvar), &(cww->cbuf_mutex));
    }
    pthread_mutex_lock(&(cwp->ucbuf_mutex));
    while (cwp->nbytes == UINT32_MAX) {
      pthread_cond_wait(&(cwp->ucbuf_condvar), &(cwp->ucbuf_mutex));
    }
    // hold both mutexes during compress operation
#endif
    uint32_t src_byte_ct = cwp->nbytes;
    uint32_t nwrite = 0;
    if (src_byte_ct) {
      char* src = cwp->ucbuf;
      // 0..15 = prefilled header
      // 16..17 = compressed size + 25
      // 18..(compressed size + 17) = compressed data
      // (compressed size + 18..21) = crc32
      // (compressed size + 22..25) = src_byte_ct
      unsigned char* dst = cww->cbuf;
      uint16_t bsize = libdeflate_deflate_compress(compressor, src, src_byte_ct, &(dst[18]), kMaxBgzfCompressedBlockSize - 26) + 25;
      const uint32_t crc = libdeflate_crc32(0, src, src_byte_ct);
      memcpy(&(dst[16]), &bsize, 2);
      memcpy(&(dst[bsize - 7]), &crc, 4);
      memcpy(&(dst[bsize - 3]), &src_byte_ct, 4);
      nwrite = S_CAST(uint32_t, bsize) + 1;
    }
    cww->nbytes = nwrite;
    const uint32_t eof = (src_byte_ct != kBgzfInputBlockSize);
    cww->eof = eof;
    cwp->nbytes = UINT32_MAX;
#ifdef _WIN32
    SetEvent(cwp->ucbuf_open_event);
    SetEvent(cww->cbuf_filled_event);
#else
    pthread_cond_signal(&(cwp->ucbuf_condvar));
    pthread_mutex_unlock(&(cwp->ucbuf_mutex));
    pthread_cond_signal(&(cww->cbuf_condvar));
    pthread_mutex_unlock(&(cww->cbuf_mutex));
#endif
    if (eof) {
      THREAD_RETURN;
    }
  }
}

THREAD_FUNC_DECL BgzfCompressWriterThread(void* raw_arg) {
  BgzfCompressStreamMain* context = S_CAST(BgzfCompressStreamMain*, raw_arg);
  const uint32_t slot_ct = context->slot_ct;
  FILE* ff = context->ff;
  BgzfCompressCommWithW** cwws = context->cwws;
  for (uint32_t slot_idx = 0; ; ++slot_idx) {
    if (slot_idx == slot_ct) {
      slot_idx = 0;
    }
    BgzfCompressCommWithW* cww = cwws[slot_idx];
#ifdef _WIN32
    WaitForSingleObject(cww->cbuf_filled_event, INFINITE);
#else
    pthread_mutex_lock(&(cww->cbuf_mutex));
    while (cww->nbytes == UINT32_MAX) {
      pthread_cond_wait(&(cww->cbuf_condvar), &(cww->cbuf_mutex));
    }
    // hold mutex during write operation
#endif
    if (ff) {
      const uint32_t nbytes = cww->nbytes;
      if (nbytes) {
        if (unlikely(!fwrite_unlocked(cww->cbuf, nbytes, 1, ff))) {
          context->write_errno = errno;
          fclose(ff);
          ff = nullptr;
          context->ff = nullptr;
        }
      }
    }
    const uint32_t eof = cww->eof;
    cww->nbytes = UINT32_MAX;
#ifdef _WIN32
    SetEvent(cww->cbuf_open_event);
#else
    pthread_cond_signal(&(cww->cbuf_condvar));
    pthread_mutex_unlock(&(cww->cbuf_mutex));
#endif
    if (eof) {
      if (ff) {
        if (unlikely(!fwrite_unlocked(kBgzfEofBlock, 28, 1, ff))) {
          context->write_errno = errno;
        }
        fclose(ff);
        ff = nullptr;
        context->ff = nullptr;
      }
      THREAD_RETURN;
    }
  }
}

// smallest power of 2 larger than kMaxBgzfCompressThreads * 4
CONSTI32(kMaxBgzfSlotCt, 64);
static_assert(kMaxBgzfCompressThreads < 16, "InitBgzfCompressStreamEx and CleanupBgzfCompressStream must be updated.");

PglErr InitBgzfCompressStreamEx(const char* out_fname, uint32_t do_append, uint32_t clvl, uint32_t compressor_thread_ct, BgzfCompressStream* cstream_ptr) {
  BgzfCompressStreamMain* bgzfp = GetBgzfp(cstream_ptr);
  if (bgzfp->ff || bgzfp->threads || (clvl > 12)) {
    // Unlike BgzfRawMtDecompressStream (which should be wrapped by
    // TextStream), this is designed for direct use in application code, so we
    // include defensive checks.
    return kPglRetImproperFunctionCall;
  }
  bgzfp->slot_ct = 0;
  bgzfp->ff = fopen(out_fname, do_append? FOPEN_AB : FOPEN_WB);
  if (unlikely(!bgzfp->ff)) {
    // slot_ct set to 0, so we'll immediately segfault on write attempt
    return kPglRetOpenFail;
  }
  if (!clvl) {
    bgzfp->write_errno = 0;
    return kPglRetSuccess;
  }
  if (compressor_thread_ct > kMaxBgzfCompressThreads) {
    compressor_thread_ct = kMaxBgzfCompressThreads;
  } else if (!compressor_thread_ct) {
    compressor_thread_ct = 1;
  }
  // smallest power of 2 >= (compressor_thread_ct * 4)
  const uint32_t slot_ct = 1 << bsru32(compressor_thread_ct * 8 - 1);
  bgzfp->slot_ct = slot_ct;
  bgzfp->compressor_thread_ct = compressor_thread_ct;
  bgzfp->write_errno = -1;  // in case of initialization failure early-exit

  const uintptr_t ptrs_aligned_size = RoundUpPow2((compressor_thread_ct + 1 + 2 * slot_ct) * sizeof(intptr_t), kCacheline);
  const uintptr_t compressor_args_aligned_size = RoundUpPow2(compressor_thread_ct * sizeof(BgzfCompressorContext), kCacheline);
  const uintptr_t cwp_aligned_size = RoundUpPow2(sizeof(BgzfCompressCommWithP), kCacheline);
  const uintptr_t cww_aligned_size = RoundUpPow2(sizeof(BgzfCompressCommWithW), kCacheline);
  unsigned char* raw_alloc;
  if (cachealigned_malloc(ptrs_aligned_size + compressor_args_aligned_size + kCacheline + slot_ct * (cwp_aligned_size + cww_aligned_size), &raw_alloc)) {
    return kPglRetNomem;
  }
  bgzfp->threads = R_CAST(pthread_t*, raw_alloc);
  ZeroPtrArr(compressor_thread_ct + 1 + 2 * slot_ct, bgzfp->threads);
  bgzfp->cwps = R_CAST(BgzfCompressCommWithP**, &(bgzfp->threads[compressor_thread_ct + 1]));
  bgzfp->cwws = R_CAST(BgzfCompressCommWithW**, &(bgzfp->cwps[slot_ct]));
  raw_alloc = &(raw_alloc[ptrs_aligned_size]);

  bgzfp->compressor_args = R_CAST(BgzfCompressorContext*, raw_alloc);
  raw_alloc = &(raw_alloc[compressor_args_aligned_size]);

  bgzfp->next_job_idxp = R_CAST(uintptr_t*, raw_alloc);
  raw_alloc = &(raw_alloc[kCacheline]);

  for (uint32_t slot_idx = 0; slot_idx != slot_ct; ++slot_idx) {
    BgzfCompressCommWithP* cwp = R_CAST(BgzfCompressCommWithP*, raw_alloc);
    bgzfp->cwps[slot_idx] = cwp;
    raw_alloc = &(raw_alloc[cwp_aligned_size]);

#ifdef _WIN32
    HANDLE new_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (unlikely(!new_event)) {
      bgzfp->unfinished_init_state = (slot_idx << 3) | 1;
      return kPglRetThreadCreateFail;
    }
    cwp->ucbuf_filled_event = new_event;

    new_event = CreateEvent(nullptr, FALSE, slot_idx? TRUE : FALSE, nullptr);
    if (unlikely(!new_event)) {
      bgzfp->unfinished_init_state = (slot_idx << 3) | 2;
      return kPglRetThreadCreateFail;
    }
    cwp->ucbuf_open_event = new_event;
#else
    if (unlikely(pthread_mutex_init(&cwp->ucbuf_mutex, nullptr))) {
      bgzfp->unfinished_init_state = (slot_idx << 3) | 1;
      return kPglRetThreadCreateFail;
    }
    if (unlikely(pthread_cond_init(&cwp->ucbuf_condvar, nullptr))) {
      bgzfp->unfinished_init_state = (slot_idx << 3) | 2;
      return kPglRetThreadCreateFail;
    }
#endif
    cwp->nbytes = UINT32_MAX;

    BgzfCompressCommWithW* cww = R_CAST(BgzfCompressCommWithW*, raw_alloc);
    bgzfp->cwws[slot_idx] = cww;
    raw_alloc = &(raw_alloc[cww_aligned_size]);
#ifdef _WIN32
    new_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (unlikely(!new_event)) {
      bgzfp->unfinished_init_state = (slot_idx << 3) | 3;
      return kPglRetThreadCreateFail;
    }
    cww->cbuf_filled_event = new_event;
    new_event = CreateEvent(nullptr, FALSE, TRUE, nullptr);
    if (unlikely(!new_event)) {
      bgzfp->unfinished_init_state = (slot_idx << 3) | 4;
      return kPglRetThreadCreateFail;
    }
    cww->cbuf_open_event = new_event;
#else
    if (unlikely(pthread_mutex_init(&cww->cbuf_mutex, nullptr))) {
      bgzfp->unfinished_init_state = (slot_idx << 3) | 3;
      return kPglRetThreadCreateFail;
    }
    if (unlikely(pthread_cond_init(&cww->cbuf_condvar, nullptr))) {
      bgzfp->unfinished_init_state = (slot_idx << 3) | 4;
      return kPglRetThreadCreateFail;
    }
#endif
    memcpy(cww->cbuf, kBgzfEofBlock, 16);
    cww->nbytes = UINT32_MAX;
    cww->eof = 0;
  }

  for (uint32_t tidx = 0; tidx != compressor_thread_ct; ++tidx) {
    struct libdeflate_compressor* lc = libdeflate_alloc_compressor(clvl);
    if (unlikely(!lc)) {
      bgzfp->unfinished_init_state = (kMaxBgzfSlotCt << 3) | tidx;
      return kPglRetNomem;
    }
    bgzfp->compressor_args[tidx].parent = bgzfp;
    bgzfp->compressor_args[tidx].lc = lc;
  }
  *(bgzfp->next_job_idxp) = 0;
  bgzfp->partial_slot_idx = 0;
  bgzfp->partial_nbytes = 0;
#ifndef _WIN32
#  ifndef __cplusplus
  pthread_attr_t smallstack_thread_attr;
  if (unlikely(pthread_attr_init(&smallstack_thread_attr))) {
    bgzfp->unfinished_init_state = (kMaxBgzfSlotCt << 3);
    return kPglRetThreadCreateFail;
  }
#  endif
#endif
  for (uint32_t tidx = 0; tidx != compressor_thread_ct; ++tidx) {
#ifdef _WIN32
    bgzfp->threads[tidx] = R_CAST(HANDLE, _beginthreadex(nullptr, kDefaultThreadStack, BgzfCompressorThread, &(bgzfp->compressor_args[tidx]), 0, nullptr));
    if (unlikely(!bgzfp->threads[tidx])) {
      bgzfp->unfinished_init_state = (kMaxBgzfSlotCt << 4) | tidx;
      return kPglRetThreadCreateFail;
    }
#else
    if (unlikely(pthread_create(&(bgzfp->threads[tidx]),
#  ifdef __cplusplus
                                &g_thread_startup.smallstack_thread_attr,
#  else
                                &smallstack_thread_attr,
#  endif
                                BgzfCompressorThread, &(bgzfp->compressor_args[tidx])))) {
#  ifndef __cplusplus
      pthread_attr_destroy(&smallstack_thread_attr);
#  endif
      bgzfp->unfinished_init_state = (kMaxBgzfSlotCt << 4) | tidx;
      return kPglRetThreadCreateFail;
    }
#endif
  }
#ifdef _WIN32
  bgzfp->threads[compressor_thread_ct] = R_CAST(HANDLE, _beginthreadex(nullptr, kDefaultThreadStack, BgzfCompressWriterThread, bgzfp, 0, nullptr));
  if (unlikely(!bgzfp->threads[compressor_thread_ct])) {
    bgzfp->unfinished_init_state = (kMaxBgzfSlotCt << 4) | compressor_thread_ct;
    return kPglRetThreadCreateFail;
  }
#else
  if (unlikely(pthread_create(&(bgzfp->threads[compressor_thread_ct]),
#  ifdef __cplusplus
                              &g_thread_startup.smallstack_thread_attr,
#  else
                              &smallstack_thread_attr,
#  endif
                              BgzfCompressWriterThread, bgzfp))) {
#  ifndef __cplusplus
    pthread_attr_destroy(&smallstack_thread_attr);
#  endif
    bgzfp->unfinished_init_state = (kMaxBgzfSlotCt << 4) | compressor_thread_ct;
    return kPglRetThreadCreateFail;
  }
#  if !defined(__cplusplus)
  pthread_attr_destroy(&smallstack_thread_attr);
#  endif
#endif
  bgzfp->write_errno = 0;
  bgzfp->unfinished_init_state = 0;
  return kPglRetSuccess;
}

BoolErr BgzfWrite(const char* buf, uintptr_t len, BgzfCompressStream* cstream_ptr) {
  BgzfCompressStreamMain* bgzfp = GetBgzfp(cstream_ptr);
  const uint32_t slot_ct = bgzfp->slot_ct;
  if (!slot_ct) {
    // No compression.
    if (unlikely(fwrite_checked(buf, len, bgzfp->ff))) {
      // If caller ignores error-return and keeps trying to write, force
      // segfault instead of silently proceeding.
      fclose(bgzfp->ff);
      bgzfp->ff = nullptr;
      return 1;
    }
    return 0;
  }
  if (bgzfp->write_errno) {
    errno = bgzfp->write_errno;
    return 1;
  }
  uint32_t slot_idx = bgzfp->partial_slot_idx;
  BgzfCompressCommWithP* cwp = bgzfp->cwps[slot_idx];
  uint32_t dst_offset = bgzfp->partial_nbytes;
  while (dst_offset + len >= kBgzfInputBlockSize) {
    const uint32_t ncopy = kBgzfInputBlockSize - dst_offset;
    memcpy(&(cwp->ucbuf[dst_offset]), buf, ncopy);
#ifdef _WIN32
    cwp->nbytes = kBgzfInputBlockSize;
    SetEvent(cwp->ucbuf_filled_event);
#else
    pthread_mutex_lock(&(cwp->ucbuf_mutex));
    cwp->nbytes = kBgzfInputBlockSize;
    pthread_cond_signal(&(cwp->ucbuf_condvar));
    pthread_mutex_unlock(&(cwp->ucbuf_mutex));
#endif
    if (++slot_idx == slot_ct) {
      slot_idx = 0;
    }
    cwp = bgzfp->cwps[slot_idx];
#ifdef _WIN32
    WaitForSingleObject(cwp->ucbuf_open_event, INFINITE);
#else
    pthread_mutex_lock(&(cwp->ucbuf_mutex));
    while (cwp->nbytes != UINT32_MAX) {
      pthread_cond_wait(&(cwp->ucbuf_condvar), &(cwp->ucbuf_mutex));
    }
    pthread_mutex_unlock(&(cwp->ucbuf_mutex));
#endif
    if (bgzfp->write_errno) {
      errno = bgzfp->write_errno;
      bgzfp->partial_slot_idx = slot_idx;
      bgzfp->partial_nbytes = 0;
      return 1;
    }
    len -= ncopy;
    buf = &(buf[ncopy]);
    dst_offset = 0;
  }
  memcpy(&(cwp->ucbuf[dst_offset]), buf, len);
  bgzfp->partial_slot_idx = slot_idx;
  bgzfp->partial_nbytes = dst_offset + len;
  return 0;
}

BoolErr CleanupBgzfCompressStream(BgzfCompressStream* cstream_ptr, PglErr* reterrp) {
  BgzfCompressStreamMain* bgzfp = GetBgzfp(cstream_ptr);
  pthread_t* threads = bgzfp->threads;
  if (threads) {
    uint32_t slot_ct = bgzfp->slot_ct;
    const uint32_t unfinished_init = bgzfp->unfinished_init_state;
    if ((!unfinished_init) || (unfinished_init >= (kMaxBgzfSlotCt << 3))) {
      const uint32_t cthread_ct = bgzfp->compressor_thread_ct;
      uint32_t launched_cthread_ct = 0;
      uint32_t launched_thread_ct = 0;
      if (!unfinished_init) {
        launched_cthread_ct = cthread_ct;
        launched_thread_ct = cthread_ct + 1;
      } else if (unfinished_init >= (kMaxBgzfSlotCt << 4)) {
        launched_cthread_ct = unfinished_init - (kMaxBgzfSlotCt << 4);
        launched_thread_ct = launched_cthread_ct;
      }
      if (launched_cthread_ct) {
        uint32_t slot_idx = bgzfp->partial_slot_idx;
        uint32_t slot_idx_stop = slot_idx + launched_cthread_ct;
        if (slot_idx_stop >= slot_ct) {
          slot_idx_stop -= slot_ct;
        }
        uint32_t nbytes = bgzfp->partial_nbytes;
        BgzfCompressCommWithP* cwp = bgzfp->cwps[slot_idx];
        do {
          // ok to send nbytes > 0 even in error-return case: the precise value
          // will be ignored, the worker will just see that it's less than
          // kBgzfInputBlockSize and return
#ifdef _WIN32
          cwp->nbytes = nbytes;
          SetEvent(cwp->ucbuf_filled_event);
#else
          pthread_mutex_lock(&(cwp->ucbuf_mutex));
          cwp->nbytes = nbytes;
          pthread_cond_signal(&(cwp->ucbuf_condvar));
          pthread_mutex_unlock(&(cwp->ucbuf_mutex));
#endif
          if (++slot_idx == slot_ct) {
            slot_idx = 0;
          }
          cwp = bgzfp->cwps[slot_idx];
#ifdef _WIN32
          WaitForSingleObject(cwp->ucbuf_open_event, INFINITE);
#else
          pthread_mutex_lock(&(cwp->ucbuf_mutex));
          while (cwp->nbytes != UINT32_MAX) {
            pthread_cond_wait(&(cwp->ucbuf_condvar), &(cwp->ucbuf_mutex));
          }
          pthread_mutex_unlock(&(cwp->ucbuf_mutex));
#endif
          nbytes = 0;
        } while (slot_idx != slot_idx_stop);
#ifdef _WIN32
        WaitForAllObjects(launched_thread_ct, threads);
        for (uint32_t tidx = 0; tidx != launched_thread_ct; ++tidx) {
          CloseHandle(threads[tidx]);
        }
#else
        for (uint32_t tidx = 0; tidx != launched_thread_ct; ++tidx) {
          pthread_join(threads[tidx], nullptr);
        }
#endif
      }
      uint32_t compressor_ct = cthread_ct;
      if (unfinished_init & (kMaxBgzfSlotCt << 3)) {
        compressor_ct = unfinished_init - (kMaxBgzfSlotCt << 3);
      }
      for (uint32_t tidx = 0; tidx != compressor_ct; ++tidx) {
        libdeflate_free_compressor(bgzfp->compressor_args[tidx].lc);
      }
    } else {
      slot_ct = unfinished_init >> 3;
      const uint32_t slot_partial_init = unfinished_init & 7;
      do {
        if (slot_partial_init == 1) {
          break;
        }
        BgzfCompressCommWithP* cwp = bgzfp->cwps[slot_ct];
#ifdef _WIN32
        CloseHandle(cwp->ucbuf_filled_event);
        if (slot_partial_init == 2) {
          break;
        }
        CloseHandle(cwp->ucbuf_open_event);
        if (slot_partial_init == 3) {
          break;
        }
        BgzfCompressCommWithW* cww = bgzfp->cwws[slot_ct];
        CloseHandle(cww->cbuf_filled_event);
#else
        pthread_mutex_destroy(&cwp->ucbuf_mutex);
        if (slot_partial_init == 2) {
          break;
        }
        pthread_cond_destroy(&cwp->ucbuf_condvar);
        if (slot_partial_init == 3) {
          break;
        }
        BgzfCompressCommWithW* cww = bgzfp->cwws[slot_ct];
        pthread_mutex_destroy(&cww->cbuf_mutex);
#endif
        assert(slot_partial_init == 4);
      } while (0);
    }
    for (uint32_t slot_idx = 0; slot_idx != slot_ct; ++slot_idx) {
      BgzfCompressCommWithP* cwp = bgzfp->cwps[slot_idx];
      BgzfCompressCommWithW* cww = bgzfp->cwws[slot_idx];
#ifdef _WIN32
      CloseHandle(cwp->ucbuf_filled_event);
      CloseHandle(cwp->ucbuf_open_event);
      CloseHandle(cww->cbuf_filled_event);
      CloseHandle(cww->cbuf_open_event);
#else
      pthread_mutex_destroy(&cwp->ucbuf_mutex);
      pthread_cond_destroy(&cwp->ucbuf_condvar);
      pthread_mutex_destroy(&cww->cbuf_mutex);
      pthread_cond_destroy(&cww->cbuf_condvar);
#endif
    }
    aligned_free(threads);
    bgzfp->threads = nullptr;
    if (unfinished_init) {
      assert(*reterrp != kPglRetSuccess);
      fclose(bgzfp->ff);
      bgzfp->ff = nullptr;
      // don't need to clear bgzfp->unfinished_init_state
    } else {
      errno = bgzfp->write_errno;
    }
  } else if (bgzfp->ff) {  // uncompressed case
    if (unlikely(fclose(bgzfp->ff))) {
      if (*reterrp == kPglRetSuccess) {
        bgzfp->ff = nullptr;
        *reterrp = kPglRetWriteFail;
        return 1;
      }
    }
    bgzfp->ff = nullptr;
  }
  return 0;
}

#ifdef __cplusplus
}
#endif
