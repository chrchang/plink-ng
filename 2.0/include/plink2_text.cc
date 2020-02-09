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
#include "plink2_text.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static inline textFILEMain* GetTxfp(textFILE* txf_ptr) {
  return &GET_PRIVATE(*txf_ptr, m);
}

static inline TextStreamMain* GetTxsp(TextStream* txs_ptr) {
  return &GET_PRIVATE(*txs_ptr, m);
}

static inline const TextStreamMain* GetTxspK(const TextStream* txs_ptr) {
  return &GET_PRIVATE(*txs_ptr, m);
}

PglErr GetFileType(const char* fname, FileCompressionType* ftype_ptr) {
  FILE* infile = fopen(fname, FOPEN_RB);
  if (unlikely(!infile)) {
    // Note that this does not print an error message (since it may be called
    // by a worker thread).
    return kPglRetOpenFail;
  }
  unsigned char buf[16];
  const uint32_t nbytes = fread_unlocked(buf, 1, 16, infile);
  if (unlikely(ferror_unlocked(infile) || fclose(infile))) {
    return kPglRetReadFail;
  }
  if (nbytes < 4) {
    *ftype_ptr = kFileUncompressed;
    return kPglRetSuccess;
  }
  uint32_t magic4;
  memcpy(&magic4, buf, 4);

  if (IsZstdFrame(magic4)) {
    *ftype_ptr = kFileZstd;
    return kPglRetSuccess;
  }
  if (S_CAST(uint16_t, magic4) != 0x8b1f) { // gzip ID1/ID2 bytes
    *ftype_ptr = kFileUncompressed;
    return kPglRetSuccess;
  }
  if ((nbytes == 16) && IsBgzfHeader(buf)) {
    *ftype_ptr = kFileBgzf;
  } else {
    *ftype_ptr = kFileGzip;
  }
  return kPglRetSuccess;
}

void EraseTextFileBase(TextFileBase* trbp) {
  trbp->consume_iter = nullptr;
  trbp->consume_stop = nullptr;
  trbp->errmsg = nullptr;
  trbp->reterr = kPglRetEof;
  trbp->ff = nullptr;
  trbp->dst = nullptr;
}

void PreinitTextFile(textFILE* txf_ptr) {
  EraseTextFileBase(&GetTxfp(txf_ptr)->base);
}

BoolErr GzRawInit(const void* buf, uint32_t nbytes, GzRawDecompressStream* gzp) {
  gzp->ds_initialized = 0;
  gzp->in = S_CAST(unsigned char*, malloc(kDecompressChunkSize));
  if (!gzp->in) {
    return 1;
  }
  z_stream* dsp = &gzp->ds;
  memcpy(gzp->in, buf, nbytes);
  dsp->next_in = gzp->in;
  dsp->avail_in = nbytes;
  dsp->zalloc = nullptr;
  dsp->zfree = nullptr;
  dsp->opaque = nullptr;
  if (unlikely(inflateInit2(dsp, MAX_WBITS | 16) != Z_OK)) {
    return 1;
  }
  gzp->ds_initialized = 1;
  return 0;
}

BoolErr ZstRawInit(const void* buf, uint32_t nbytes, ZstRawDecompressStream* zstp) {
  zstp->ib.src = malloc(kDecompressChunkSize);
  if (unlikely(!zstp->ib.src)) {
    zstp->ds = nullptr;
    return 1;
  }
  zstp->ds = ZSTD_createDStream();
  if (unlikely(!zstp->ds)) {
    return 1;
  }
  memcpy(K_CAST(void*, zstp->ib.src), buf, nbytes);
  zstp->ib.size = nbytes;
  zstp->ib.pos = 0;
  return 0;
}

const char kShortErrRfileAlreadyOpen[] = "TextFileOpenInternal can't be called on an already-open file";
const char kShortErrRfileEnforcedMaxBlenTooSmall[] = "TextFileOpenInternal: enforced_max_line_blen too small (must be at least max(1 MiB, dst_capacity - 1 MiB))";
const char kShortErrRfileDstCapacityTooSmall[] = "TextFileOpenInternal: dst_capacity too small (2 MiB minimum)";

PglErr TextFileOpenInternal(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, char* dst, textFILEMain* txfp, TextStreamMain* txsp) {
  PglErr reterr = kPglRetSuccess;
  TextFileBase* trbp;
  if (txfp) {
    trbp = &txfp->base;
  } else {
    trbp = &txsp->base;
  }
  {
    // 1. Open file, get type.
    if (unlikely(trbp->ff)) {
      reterr = kPglRetImproperFunctionCall;
      trbp->errmsg = kShortErrRfileAlreadyOpen;
      goto TextFileOpenInternal_ret_1;
    }
    if (enforced_max_line_blen || txfp) {
      if (unlikely(enforced_max_line_blen < kDecompressMinBlen)) {
        reterr = kPglRetImproperFunctionCall;
        trbp->errmsg = kShortErrRfileEnforcedMaxBlenTooSmall;
        goto TextFileOpenInternal_ret_1;
      }
      if (dst) {
        if (unlikely(dst_capacity < kDecompressMinCapacity)) {
          reterr = kPglRetImproperFunctionCall;
          trbp->errmsg = kShortErrRfileDstCapacityTooSmall;
          goto TextFileOpenInternal_ret_1;
        }
        if (unlikely(enforced_max_line_blen + kDecompressChunkSize < dst_capacity)) {
          reterr = kPglRetImproperFunctionCall;
          trbp->errmsg = kShortErrRfileEnforcedMaxBlenTooSmall;
          goto TextFileOpenInternal_ret_1;
        }
      }
    } else {
      // token-reading mode.  dst == nullptr not currently supported.
      assert(dst && (dst_capacity == kTokenStreamBlen));
    }
    trbp->ff = fopen(fname, FOPEN_RB);
    if (unlikely(!trbp->ff)) {
      goto TextFileOpenInternal_ret_OPEN_FAIL;
    }
    trbp->file_type = kFileUncompressed;
    if (dst) {
      trbp->dst_owned_by_consumer = 1;
      trbp->dst_capacity = dst_capacity;
    } else {
      dst = S_CAST(char*, malloc(kDecompressMinCapacity));
      if (unlikely(dst == nullptr)) {
        goto TextFileOpenInternal_ret_NOMEM;
      }
      trbp->dst_owned_by_consumer = 0;
      trbp->dst_capacity = kDecompressMinCapacity;
    }
    trbp->dst = dst;
    uint32_t nbytes = fread_unlocked(dst, 1, 16, trbp->ff);
    trbp->dst_len = nbytes;
    trbp->enforced_max_line_blen = enforced_max_line_blen;
    trbp->consume_iter = dst;
    trbp->consume_stop = dst;
    if (nbytes >= 4) {
      const uint32_t magic4 = *R_CAST(uint32_t*, dst);
      if (IsZstdFrame(magic4)) {
        trbp->dst_len = 0;
        trbp->file_type = kFileZstd;
        ZstRawDecompressStream* zstp;
        if (txfp) {
          zstp = &txfp->rds.zst;
        } else {
          zstp = &txsp->rds.zst;
        }
        if (unlikely(ZstRawInit(dst, nbytes, zstp))) {
          goto TextFileOpenInternal_ret_NOMEM;
        }
      } else if ((magic4 << 8) == 0x088b1f00) {
        // gzip ID1/ID2 bytes, deflate compression method
        trbp->dst_len = 0;
        if ((nbytes == 16) && IsBgzfHeader(dst)) {
          trbp->file_type = kFileBgzf;
          if (txfp) {
            BgzfRawDecompressStream* bgzfp = &txfp->rds.bgzf;
            bgzfp->in = S_CAST(unsigned char*, malloc(kDecompressChunkSize));
            if (unlikely(!bgzfp->in)) {
              bgzfp->ldc = nullptr;
              goto TextFileOpenInternal_ret_NOMEM;
            }
            bgzfp->ldc = libdeflate_alloc_decompressor();
            if (!bgzfp->ldc) {
              goto TextFileOpenInternal_ret_NOMEM;
            }
            memcpy(bgzfp->in, dst, nbytes);
            bgzfp->in_size = nbytes;
            bgzfp->in_pos = 0;
          } else {
            reterr = BgzfRawMtStreamInit(dst, txsp->decompress_thread_ct, trbp->ff, nullptr, &txsp->rds.bgzf, &trbp->errmsg);
            if (unlikely(reterr)) {
              goto TextFileOpenInternal_ret_1;
            }
          }
        } else {
          trbp->file_type = kFileGzip;
          GzRawDecompressStream* gzp;
          if (txfp) {
            gzp = &txfp->rds.gz;
          } else {
            gzp = &txsp->rds.gz;
          }
          if (unlikely(GzRawInit(dst, nbytes, gzp))) {
            goto TextFileOpenInternal_ret_NOMEM;
          }
        }
      }
    } else if (!nbytes) {
      if (unlikely(!feof_unlocked(trbp->ff))) {
        goto TextFileOpenInternal_ret_READ_FAIL;
      }
      // May as well accept this.
      // Don't jump to ret_1 since we're setting txfp->reterr to a different
      // value than we're returning.
      trbp->reterr = kPglRetEof;
      return kPglRetSuccess;
    }
  }
  while (0) {
  TextFileOpenInternal_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TextFileOpenInternal_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    trbp->errmsg = strerror(errno);
    break;
  TextFileOpenInternal_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    trbp->errmsg = strerror(errno);
    break;
  }
 TextFileOpenInternal_ret_1:
  trbp->reterr = reterr;
  return reterr;
}

PglErr TextFileOpenEx(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, char* dst, textFILE* txf_ptr) {
  return TextFileOpenInternal(fname, enforced_max_line_blen, dst_capacity, dst, GetTxfp(txf_ptr), nullptr);
}

// Set enforced_max_line_blen == 0 in the token-reading case.
BoolErr IsPathologicallyLongLineOrToken(const char* line_start, const char* load_start, const char* known_line_end, uint32_t enforced_max_line_blen) {
  if (enforced_max_line_blen) {
    // Preconditions:
    // * No \n in [line_start, load_start).
    // * (known_line_end - load_start) is usually <= enforced_max_line_blen,
    //   and never much larger.  Not a hard requirement, but it's better to
    //   enforce the line-length limit during line iteration outside this
    //   regime to avoid duplicating work.
    if (S_CAST(uintptr_t, known_line_end - line_start) <= enforced_max_line_blen) {
      return 0;
    }
    const uint32_t already_scanned_byte_ct = load_start - line_start;
    if (unlikely(already_scanned_byte_ct >= enforced_max_line_blen)) {
      return 1;
    }
    const char* memchr_result = S_CAST(const char*, memchr(load_start, '\n', enforced_max_line_blen - already_scanned_byte_ct));
    if (unlikely(!memchr_result)) {
      return 1;
    }
    // If we've found a line with terminal \n at or after this address, there
    // are <= enforced_max_line_blen bytes left, so no remaining line can be
    // longer.
    const char* memchr_result_thresh = known_line_end - (enforced_max_line_blen + 1);
    while (1) {
      if (memchr_result >= memchr_result_thresh) {
        return 0;
      }
      memchr_result = S_CAST(const char*, memchr(&(memchr_result[1]), '\n', enforced_max_line_blen));
      if (unlikely(!memchr_result)) {
        return 1;
      }
    }
  }
  if (S_CAST(uintptr_t, known_line_end - line_start) <= kMaxTokenBlen) {
    return 0;
  }
  const uint32_t already_scanned_byte_ct = load_start - line_start;
  if (unlikely(already_scanned_byte_ct >= kMaxTokenBlen)) {
    return 1;
  }
  // No loop needed for now, since token-scanning buffer sizes are hardcoded.
  //
  // Replace with a forward-scanning version of this functionality when
  // available ("FirstPostspaceBoundedFar"?)
  return (LastSpaceOrEoln(load_start, kMaxTokenBlen - already_scanned_byte_ct) == nullptr);
}

const char kShortErrRfileTruncatedGz[] = "GzRawStreamRead: gzipped file appears to be truncated";

PglErr GzRawStreamRead(char* dst_end, FILE* ff, GzRawDecompressStream* gzp, char** dst_iterp, const char** errmsgp) {
  z_stream* dsp = &gzp->ds;
  if ((!dsp->avail_in) && feof_unlocked(ff)) {
    return kPglRetSuccess;
  }
  char* dst_iter = *dst_iterp;
  do {
    int zerr = Z_OK;
    if (dsp->avail_in) {  // can be zero after TextRewind()
      dsp->next_out = R_CAST(unsigned char*, dst_iter);
      dsp->avail_out = dst_end - dst_iter;
      zerr = inflate(dsp, Z_SYNC_FLUSH);
      if (unlikely((zerr < 0) || (zerr == Z_NEED_DICT))) {
        if (dsp->msg) {
          *errmsgp = dsp->msg;
        } else {
          *errmsgp = zError(zerr);
        }
        return kPglRetDecompressFail;
      }
      dst_iter = R_CAST(char*, dsp->next_out);
      if (dsp->avail_in) {
        assert(dst_iter == dst_end);
        break;
      }
    }
    const uint32_t nbytes = fread_unlocked(gzp->in, 1, kDecompressChunkSize, ff);
    dsp->next_in = gzp->in;
    dsp->avail_in = nbytes;
    if (!nbytes) {
      if (unlikely(!feof_unlocked(ff))) {
        *errmsgp = strerror(errno);
        return kPglRetReadFail;
      }
      if (unlikely(zerr == Z_OK)) {
        *errmsgp = kShortErrRfileTruncatedGz;
        return kPglRetDecompressFail;
      }
      // Normal EOF.
      break;
    }
  } while (dst_iter != dst_end);
  *dst_iterp = dst_iter;
  return kPglRetSuccess;
}

PglErr ZstRawStreamRead(char* dst_end, FILE* ff, ZstRawDecompressStream* zstp, char** dst_iterp, const char** errmsgp) {
  if ((!zstp->ib.size) && feof_unlocked(ff)) {
    return kPglRetSuccess;
  }
  // Sequentially dependent blocks limited to ~128 KiB.
  char* dst_iter = *dst_iterp;
  while (1) {
    ZSTD_outBuffer zob = {R_CAST(unsigned char*, dst_iter), S_CAST(size_t, dst_end - dst_iter), 0};
    // ib.size == 0 ok, no need to special-case rewind.
    const uintptr_t read_size_hint = ZSTD_decompressStream(zstp->ds, &zob, &zstp->ib);
    if (unlikely(ZSTD_isError(read_size_hint))) {
      *errmsgp = ZSTD_getErrorName(read_size_hint);
      return kPglRetDecompressFail;
    }
    dst_iter = &(dst_iter[zob.pos]);
    if (dst_iter == dst_end) {
      break;
    }
    // Decoder has flushed everything it could.  Either we're at EOF, or we
    // must load more.
    unsigned char* in = S_CAST(unsigned char*, K_CAST(void*, zstp->ib.src));
    const uint32_t n_inbytes = zstp->ib.size - zstp->ib.pos;
    memmove(in, &(in[zstp->ib.pos]), n_inbytes);
    unsigned char* load_start = &(in[n_inbytes]);
    const uint32_t nbytes = fread_unlocked(load_start, 1, kDecompressChunkSize - n_inbytes, ff);
    if (unlikely(ferror_unlocked(ff))) {
      *errmsgp = strerror(errno);
      return kPglRetReadFail;
    }
    zstp->ib.pos = 0;
    zstp->ib.size = nbytes + n_inbytes;
    if (!nbytes) {
      if (unlikely(n_inbytes)) {
        *errmsgp = kShortErrZstdPrefixUnknown;
        return kPglRetDecompressFail;
      }
      break;
    }
  }
  *dst_iterp = dst_iter;
  return kPglRetSuccess;
}

const char kShortErrLongLine[] = "Pathologically long line";
const char kShortErrInteriorEmptyLine[] = "Unexpected interior empty line";

PglErr TextFileAdvance(textFILE* txf_ptr) {
  textFILEMain* txfp = GetTxfp(txf_ptr);
  TextFileBase* basep = &txfp->base;
  if (basep->reterr) {
    return basep->reterr;
  }
  PglErr reterr = kPglRetSuccess;
  {
    char* line_start = basep->consume_stop;
    assert(basep->consume_iter == line_start);
    char* dst = basep->dst;
    char* dst_load_start;
    while (1) {
      const uint32_t dst_offset = line_start - dst;
      const uint32_t dst_rem = basep->dst_len - dst_offset;
      // (dst_rem guaranteed to be < basep->enforced_max_line_blen here, since
      // otherwise we error out earlier.)
      // Two cases:
      // 1. Move (possibly empty) unfinished line to the beginning of the
      //    buffer.
      // 2. Resize the buffer/report out-of-memory.
      if (dst_rem < basep->dst_capacity - kDecompressChunkSize) {
        if (dst_offset) {
          memmove(dst, line_start, dst_rem);
        }
      } else {
        if (unlikely(basep->dst_owned_by_consumer)) {
          goto TextFileAdvance_ret_NOMEM;
        }
        uint32_t next_dst_capacity = basep->enforced_max_line_blen + kDecompressChunkSize;
        if ((next_dst_capacity / 2) > basep->dst_capacity) {
          next_dst_capacity = basep->dst_capacity * 2;
        }
#ifndef __LP64__
        if (next_dst_capacity >= 0x80000000U) {
          goto TextFileAdvance_ret_NOMEM;
        }
#endif
        char* dst_next;
        if (!dst_offset) {
          dst_next = S_CAST(char*, realloc(dst, next_dst_capacity));
          if (unlikely(!dst_next)) {
            goto TextFileAdvance_ret_NOMEM;
          }
        } else {
          dst_next = S_CAST(char*, malloc(next_dst_capacity));
          if (unlikely(!dst_next)) {
            goto TextFileAdvance_ret_NOMEM;
          }
          memcpy(dst_next, line_start, dst_rem);
        }
        basep->dst = dst_next;
        dst = dst_next;
      }
      line_start = dst;
      dst_load_start = &(dst[dst_rem]);
      FILE* ff = basep->ff;
      char* dst_iter = dst_load_start;
      // We don't want to always fill the entire buffer here.  The main plink2
      // use case of textFILE is to just peek at an unknown-length header line
      // with a maximal-length line-load buffer, compute a
      // legitimate-line-length bound, and then (move-)construct a TextStream
      // with the shorter buffer size.
      // Instead, we load up to the smallest power of 2 >= (dst_rem + 1 MiB).
      uintptr_t stop_offset = (2 * k1LU) << bsru32(dst_rem + kDecompressChunkSize - 1);
      if (stop_offset > basep->dst_capacity) {
        stop_offset = basep->dst_capacity;
      }
      char* dst_stop = &(dst[stop_offset]);
      basep->consume_iter = dst;
      switch (basep->file_type) {
      case kFileUncompressed:
        {
          uint32_t rlen = dst_stop - dst_iter;
          if (rlen > kMaxBytesPerIO) {
            // We need to know how many bytes were read, so fread_checked()
            // doesn't work.
            // This is an if-statement instead of a while loop since rlen can
            // never be larger than 2 * kMaxBytesPerIO.
            const uint32_t nbytes = fread_unlocked(dst_iter, 1, kMaxBytesPerIO, ff);
            if (nbytes < kMaxBytesPerIO) {
              if (unlikely(ferror_unlocked(ff))) {
                goto TextFileAdvance_ret_READ_FAIL;
              }
              basep->dst_len = nbytes + dst_rem;
              break;
            }
            rlen -= kMaxBytesPerIO;
            dst_iter = &(dst_iter[kMaxBytesPerIO]);
          }
          const uint32_t nbytes = fread_unlocked(dst_iter, 1, rlen, ff);
          if (unlikely(ferror_unlocked(ff))) {
            goto TextFileAdvance_ret_READ_FAIL;
          }
          dst_iter = &(dst_iter[nbytes]);
          break;
        }
      case kFileGzip:
        {
          reterr = GzRawStreamRead(dst_stop, ff, &txfp->rds.gz, &dst_iter, &basep->errmsg);
          if (unlikely(reterr)) {
            goto TextFileAdvance_ret_1;
          }
          break;
        }
      case kFileBgzf:
        {
          // Fully independent blocks limited to 64 KiB.
          // probable todo: move this to a BgzfRawStreamRead() function in
          // plink2_bgzf (and move ZstRawStreamRead() to plink2_zstfile).
          BgzfRawDecompressStream* bgzfp = &txfp->rds.bgzf;
          if ((!bgzfp->in_size) && feof_unlocked(ff)) {
            break;
          }
          struct libdeflate_decompressor* ldc = bgzfp->ldc;
          unsigned char* in = bgzfp->in;
          unsigned char* in_iter = &(in[bgzfp->in_pos]);
          unsigned char* in_end = &(in[bgzfp->in_size]);
          while (1) {
            uint32_t n_inbytes = in_end - in_iter;
            if (n_inbytes > 25) {
              if (unlikely(!IsBgzfHeader(in_iter))) {
                goto TextFileAdvance_ret_INVALID_BGZF;
              }
#  ifdef __arm__
#    error "Unaligned accesses in TextFileAdvance()."
#  endif
              const uint32_t bsize_minus1 = *R_CAST(uint16_t*, &(in_iter[16]));
              if (unlikely(bsize_minus1 < 25)) {
                goto TextFileAdvance_ret_INVALID_BGZF;
              }
              if (bsize_minus1 < n_inbytes) {
                // We have at least one fully-loaded compressed block.
                // Decompress it if we have enough space.
                const uint32_t in_size = bsize_minus1 - 25;
                const uint32_t out_size = *R_CAST(uint32_t*, &(in_iter[in_size + 22]));
                if (unlikely(out_size > 65536)) {
                  goto TextFileAdvance_ret_INVALID_BGZF;
                }
                if (out_size > S_CAST(uintptr_t, dst_stop - dst_iter)) {
                  break;
                }
                if (unlikely(libdeflate_deflate_decompress(ldc, &(in_iter[18]), in_size, dst_iter, out_size, nullptr))) {
                  goto TextFileAdvance_ret_INVALID_BGZF;
                }
                in_iter = &(in_iter[bsize_minus1 + 1]);
                dst_iter = &(dst_iter[out_size]);
                continue;
              }
            }
            // Either we're at EOF, or we must load more.
            memmove(in, in_iter, n_inbytes);
            unsigned char* load_start = &(in[n_inbytes]);
            const uint32_t nbytes = fread_unlocked(load_start, 1, kDecompressChunkSize - n_inbytes, ff);
            if (unlikely(ferror_unlocked(ff))) {
              goto TextFileAdvance_ret_READ_FAIL;
            }
            in_iter = in;
            in_end = &(load_start[nbytes]);
            bgzfp->in_size = in_end - in;
            if (!nbytes) {
              if (unlikely(n_inbytes)) {
                goto TextFileAdvance_ret_INVALID_BGZF;
              }
              break;
            }
          }
          bgzfp->in_pos = in_iter - in;
          dst_stop = dst_iter;
          break;
        }
      case kFileZstd:
        {
          reterr = ZstRawStreamRead(dst_stop, ff, &txfp->rds.zst, &dst_iter, &basep->errmsg);
          if (unlikely(reterr)) {
            goto TextFileAdvance_ret_1;
          }
          break;
        }
      }
      basep->dst_len = dst_iter - dst;
      if (!basep->dst_len) {
        goto TextFileAdvance_ret_EOF;
      }
      if (dst_iter != dst_stop) {
        // If last character of file isn't a newline, append one to simplify
        // downstream code.
        if (dst_iter[-1] != '\n') {
          *dst_iter++ = '\n';
          basep->dst_len += 1;
        }
        basep->consume_stop = dst_iter;
        break;
      }
      char* last_byte_ptr = Memrchr(dst_load_start, '\n', dst_iter - dst_load_start);
      if (last_byte_ptr) {
        basep->consume_stop = &(last_byte_ptr[1]);
        break;
      }
      // Buffer is full, and no '\n' is present.  Restart the loop and try to
      // load more data (extending the buffer if necessary), if we aren't
      // already at/past the line-length limit.
      if (basep->dst_len >= basep->enforced_max_line_blen) {
        goto TextFileAdvance_ret_LONG_LINE;
      }
    }
    if (unlikely(IsPathologicallyLongLineOrToken(dst, dst_load_start, basep->consume_stop, basep->enforced_max_line_blen))) {
      goto TextFileAdvance_ret_LONG_LINE;
    }
  }
  while (0) {
  TextFileAdvance_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TextFileAdvance_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    basep->errmsg = strerror(errno);
    break;
  TextFileAdvance_ret_LONG_LINE:
    basep->errmsg = kShortErrLongLine;
    reterr = kPglRetMalformedInput;
    break;
  TextFileAdvance_ret_INVALID_BGZF:
    basep->errmsg = kShortErrInvalidBgzf;
    reterr = kPglRetDecompressFail;
    break;
  TextFileAdvance_ret_EOF:
    reterr = kPglRetEof;
    break;
  }
 TextFileAdvance_ret_1:
  basep->reterr = reterr;
  return reterr;
}

PglErr TextFileOnlyEmptyLinesLeft(textFILE* txf_ptr) {
  TextFileBase* basep = &GetTxfp(txf_ptr)->base;
  char* line_start = basep->consume_iter;
  while (1) {
    if (line_start == basep->consume_stop) {
      basep->consume_iter = line_start;
      PglErr reterr = TextFileAdvance(txf_ptr);
      if (reterr) {
        return reterr;
      }
      line_start = basep->consume_iter;
    }
    line_start = FirstNonTspace(line_start);
    if (unlikely(!IsEolnKns(*line_start))) {
      basep->reterr = kPglRetMalformedInput;
      basep->errmsg = kShortErrInteriorEmptyLine;
      return kPglRetMalformedInput;
    }
    line_start = AdvPastDelim(line_start, '\n');
  }
}

void TextFileRewind(textFILE* txf_ptr) {
  textFILEMain* txfp = GetTxfp(txf_ptr);
  TextFileBase* basep = &txfp->base;
  if ((!basep->ff) || ((basep->reterr) && (basep->reterr != kPglRetEof))) {
    return;
  }
  rewind(basep->ff);
  basep->reterr = kPglRetSuccess;
  basep->dst_len = 0;
  basep->consume_iter = basep->dst;
  basep->consume_stop = basep->dst;
  if (basep->file_type != kFileUncompressed) {
    if (basep->file_type == kFileGzip) {
      txfp->rds.gz.ds.avail_in = 0;
#ifdef NDEBUG
      inflateReset(&txfp->rds.gz.ds);
#else
      const int errcode = inflateReset(&txfp->rds.gz.ds);
      assert(errcode == Z_OK);
#endif
    } else if (basep->file_type == kFileBgzf) {
      txfp->rds.bgzf.in_size = 0;
      txfp->rds.bgzf.in_pos = 0;
    } else {
      // kFileZstd
      txfp->rds.zst.ib.size = 0;
      txfp->rds.zst.ib.pos = 0;
      ZSTD_DCtx_reset(txfp->rds.zst.ds, ZSTD_reset_session_only);
    }
  }
}

BoolErr CleanupTextFile(textFILE* txf_ptr, PglErr* reterrp) {
  textFILEMain* txfp = GetTxfp(txf_ptr);
  TextFileBase* basep = &txfp->base;
  basep->consume_iter = nullptr;
  basep->consume_stop = nullptr;
  basep->reterr = kPglRetEof;
  basep->errmsg = nullptr;
  if (basep->dst && (!basep->dst_owned_by_consumer)) {
    free(basep->dst);
    basep->dst = nullptr;
  }
  if (basep->ff) {
    if (basep->file_type != kFileUncompressed) {
      if (basep->file_type == kFileZstd) {
        if (txfp->rds.zst.ib.src) {
          free_const(txfp->rds.zst.ib.src);
          txfp->rds.zst.ib.src = nullptr;
        }
        if (txfp->rds.zst.ds) {
          ZSTD_freeDStream(txfp->rds.zst.ds);
          txfp->rds.zst.ds = nullptr;
        }
      } else if (basep->file_type == kFileBgzf) {
        if (txfp->rds.bgzf.in) {
          free(txfp->rds.bgzf.in);
          txfp->rds.bgzf.in = nullptr;
        }
        if (txfp->rds.bgzf.ldc) {
          libdeflate_free_decompressor(txfp->rds.bgzf.ldc);
          txfp->rds.bgzf.ldc = nullptr;
        }
      } else {
        // plain gzip
        if (txfp->rds.gz.in) {
          free(txfp->rds.gz.in);
          txfp->rds.gz.in = nullptr;
        }
        if (txfp->rds.gz.ds_initialized) {
          inflateEnd(&txfp->rds.gz.ds);
        }
      }
    }
    if (unlikely(fclose_null(&basep->ff))) {
      if (!reterrp) {
        return 1;
      }
      if (*reterrp == kPglRetSuccess) {
        *reterrp = kPglRetReadFail;
        return 1;
      }
    }
  }
  return 0;
}


void PreinitTextStream(TextStream* txs_ptr) {
  TextStreamMain* txsp = GetTxsp(txs_ptr);
  EraseTextFileBase(&txsp->base);
  txsp->syncp = nullptr;
}

// This type of code is especially bug-prone (ESR would call it a "defect
// attractor").  Goal is to get it right, and fast enough to be a major win
// over gzgets()... and then not worry about it again for years.
THREAD_FUNC_DECL TextStreamThread(void* raw_arg) {
  TextStreamMain* context = S_CAST(TextStreamMain*, raw_arg);
  TextFileBase* basep = &context->base;
  TextStreamSync* syncp = context->syncp;
  FileCompressionType file_type = basep->file_type;
  RawMtDecompressStream* rdsp = &context->rds;
  FILE* ff = basep->ff;
  char* buf = basep->dst;
  char* buf_end = &(buf[basep->dst_capacity]);
  char* cur_block_start = basep->consume_stop;
  char* read_head = &(buf[basep->dst_len]);

  // We can either be reading/decompressing into memory past the bytes passed
  // to the consumer, or we can be doing it before those bytes.
  // In the first case, read_stop is buf_end, but it gets changed to the
  // latest value of consume_tail when we return to the front of the buffer.
  // In the second case, read_stop is the position of the first passed byte.
  char* read_stop = buf_end;
#ifdef _WIN32
  CRITICAL_SECTION* critical_sectionp = &syncp->critical_section;
  HANDLE reader_progress_event = syncp->reader_progress_event;
  HANDLE consumer_progress_event = syncp->consumer_progress_event;
#else
  pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
  pthread_cond_t* reader_progress_condvarp = &syncp->reader_progress_condvar;
  pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
#endif
  const uint32_t enforced_max_line_blen = basep->enforced_max_line_blen;
  const char* new_fname = nullptr;
  const uint32_t is_token_stream = (enforced_max_line_blen == 0);
  while (1) {
    TxsInterrupt interrupt = kTxsInterruptNone;
    PglErr reterr;
    TxsInterrupt min_interrupt;
    while (1) {
      uintptr_t read_attempt_size = read_stop - read_head;
      if (!read_attempt_size) {
        const uint32_t memmove_required = (read_stop == buf_end);
        if (unlikely((cur_block_start == buf) && memmove_required)) {
          // May need to modify this predicate if we ever allow is_token_stream
          // && !dst_owned_by_consumer.
          const uint32_t prev_capacity = buf_end - buf;
          if (basep->dst_owned_by_consumer || (prev_capacity >= enforced_max_line_blen)) {
            goto TextStreamThread_LONG_LINE;
          }
          // Try to expand buffer.
          uint32_t next_dst_capacity = enforced_max_line_blen + kDecompressChunkSize;
          if ((next_dst_capacity / 2) > basep->dst_capacity) {
            next_dst_capacity = basep->dst_capacity * 2;
          }
#ifndef __LP64__
          if (next_dst_capacity >= 0x80000000U) {
            goto TextStreamThread_NOMEM;
          }
#endif
          char* dst_next = S_CAST(char*, realloc(buf, next_dst_capacity));
          if (unlikely(!dst_next)) {
            goto TextStreamThread_NOMEM;
          }
#ifdef _WIN32
          EnterCriticalSection(critical_sectionp);
#else
          pthread_mutex_lock(sync_mutexp);
#endif
          basep->dst = dst_next;
          basep->dst_capacity = next_dst_capacity;
          syncp->consume_tail = dst_next;
          syncp->available_end = dst_next;
          syncp->dst_reallocated = 1;
#ifdef _WIN32
          LeaveCriticalSection(critical_sectionp);
#else
          pthread_mutex_unlock(sync_mutexp);
#endif
          buf = dst_next;
          buf_end = &(buf[next_dst_capacity]);
          cur_block_start = buf;
          read_head = &(buf[prev_capacity]);
          read_stop = &(buf[next_dst_capacity]);
          continue;
        }
        // We cannot continue reading forward.  Cases:
        // 1. read_stop == buf_end, cur_block_start != buf.  This means we're
        //    in the middle of reading/decompressing a long line, and want to
        //    wait for consume_tail == cur_block_start, so we can memmove all
        //    the bytes back and continue reading forward.  (Tried
        //    relaxing this to
        //      consume_tail >= (buf_end - cur_block_start) + margin
        //    for various values of margin, but that didn't make a meaningful
        //    difference.)
        // 2. read_stop == buf_end, cur_block_start == buf.  We failed with a
        //    long-line error here.
        // 3. read_stop < buf_end (usual case).  This means the consumer may
        //    not be done handling some bytes-in-front we handed off earlier.
        //    We are waiting for consume_tail <= cur_block_start, which means
        //    all bytes in front have been consumed and we're free to continue
        //    reading forward.
        char* latest_consume_tail;
#ifdef _WIN32
        // bugfix (7 May 2018): when consumer thread is waiting with
        // syncp->consume_tail == cur_block_start, read_stop is near but not at
        // buf_end, and there's no '\n' in the subsequent read, we can reach
        // here a second time without releasing the consumer, so we'd enter
        // deadlock if we unconditionally wait on consumer_progress_event (and
        // in the Linux/OS X case, we'd be waiting for a spurious wakeup to
        // save us).
        // However, if memmove_required isn't true, we have to wait first; see
        // the 21 Mar bugfix.
        if (!memmove_required) {
          goto TextStreamThread_wait_first;
        }
        while (1) {
          EnterCriticalSection(critical_sectionp);
          interrupt = syncp->interrupt;
          if (interrupt != kTxsInterruptNone) {
            goto TextStreamThread_INTERRUPT;
          }
          latest_consume_tail = syncp->consume_tail;
          if (memmove_required) {
            if (latest_consume_tail == cur_block_start) {
              syncp->consume_tail = buf;
              syncp->available_end = buf;
              break;
            }
          } else if (latest_consume_tail <= cur_block_start) {
            break;
          }
          LeaveCriticalSection(critical_sectionp);
        TextStreamThread_wait_first:
          WaitForSingleObject(consumer_progress_event, INFINITE);
        }
        // bugfix (23 Mar 2018): didn't always leave the critical section
        LeaveCriticalSection(critical_sectionp);
#else
        pthread_mutex_lock(sync_mutexp);
        if (!memmove_required) {
          // Wait for all bytes in front of read_stop to be consumed.
          goto TextStreamThread_wait_first;
        }
        while (1) {
          interrupt = syncp->interrupt;
          if (interrupt != kTxsInterruptNone) {
            goto TextStreamThread_INTERRUPT;
          }
          latest_consume_tail = syncp->consume_tail;
          if (memmove_required) {
            if (latest_consume_tail == cur_block_start) {
              // All bytes have been consumed; memmove is now safe.
              // bugfix (2 Oct 2018): Previously, this just set
              // syncp->cur_circular_end = cur_block_start, but that created
              // TWO consume_iter == available_end == cur_circular_end cases,
              // one of which was handled incorrectly.
              syncp->consume_tail = buf;
              syncp->available_end = buf;
              break;
            }
            // There are bytes behind cur_block_start that haven't been
            // consumed yet.  This is possible on the first iteration through
            // the loop, since consumer_progress_state may have been set for a
            // reason we aren't interested in.

          } else if (latest_consume_tail <= cur_block_start) {
            // All bytes in front of read_stop have been consumed.
            break;
          }
        TextStreamThread_wait_first:
          while (!syncp->consumer_progress_state) {
            pthread_cond_wait(consumer_progress_condvarp, sync_mutexp);
          }
          syncp->consumer_progress_state = 0;
        }
        pthread_mutex_unlock(sync_mutexp);
#endif
        if (read_stop == buf_end) {
          const uint32_t cur_memmove_len = buf_end - cur_block_start;
          memmove(buf, cur_block_start, cur_memmove_len);
          cur_block_start = buf;
          read_head = &(buf[cur_memmove_len]);
        } else {
          read_stop = buf_end;
        }
        continue;
      }
      if (read_attempt_size > kDecompressChunkSize) {
        read_attempt_size = kDecompressChunkSize;
      }
      char* cur_read_end = read_head;
      char* cur_read_stop = &(read_head[read_attempt_size]);
      switch (file_type) {
      case kFileUncompressed:
        {
          cur_read_end += fread_unlocked(read_head, 1, read_attempt_size, ff);
          if (unlikely(ferror_unlocked(ff))) {
            goto TextStreamThread_READ_FAIL;
          }
          break;
        }
      case kFileGzip:
        {
          reterr = GzRawStreamRead(cur_read_stop, ff, &rdsp->gz, &cur_read_end, &syncp->errmsg);
          if (unlikely(reterr)) {
            goto TextStreamThread_MISC_FAIL;
          }
          break;
        }
      case kFileBgzf:
        {
          reterr = BgzfRawMtStreamRead(R_CAST(unsigned char*, cur_read_stop), &rdsp->bgzf, R_CAST(unsigned char**, &cur_read_end), &syncp->errmsg);
          if (unlikely(reterr)) {
            goto TextStreamThread_MISC_FAIL;
          }
          break;
        }
      case kFileZstd:
        {
          reterr = ZstRawStreamRead(cur_read_stop, ff, &rdsp->zst, &cur_read_end, &syncp->errmsg);
          if (unlikely(reterr)) {
            goto TextStreamThread_MISC_FAIL;
          }
          break;
        }
      }
      if (cur_read_end < cur_read_stop) {
        char* final_read_head = cur_read_end;
        if (cur_block_start != final_read_head) {
          if (final_read_head[-1] != '\n') {
            // Append '\n' so consumer can always use rawmemchr(., '\n') to
            // find the end of the current line.
            *final_read_head++ = '\n';
          }
        }
        // Still want to consistently enforce max line/token length.
        if (unlikely(IsPathologicallyLongLineOrToken(cur_block_start, read_head, final_read_head, enforced_max_line_blen))) {
          goto TextStreamThread_LONG_LINE;
        }
        read_head = final_read_head;
        goto TextStreamThread_EOF;
      }
      char* last_byte_ptr;
      if (!is_token_stream) {
        last_byte_ptr = Memrchr(read_head, '\n', read_attempt_size);
      } else {
        last_byte_ptr = LastSpaceOrEoln(read_head, read_attempt_size);
      }
      if (last_byte_ptr) {
        char* next_available_end = &(last_byte_ptr[1]);
        if (unlikely(IsPathologicallyLongLineOrToken(cur_block_start, read_head, next_available_end, enforced_max_line_blen))) {
          goto TextStreamThread_LONG_LINE;
        }
#ifdef _WIN32
        EnterCriticalSection(critical_sectionp);
#else
        pthread_mutex_lock(sync_mutexp);
#endif
        interrupt = syncp->interrupt;
        if (interrupt != kTxsInterruptNone) {
          goto TextStreamThread_INTERRUPT;
        }
        char* latest_consume_tail = syncp->consume_tail;
        const uint32_t all_later_bytes_consumed = (latest_consume_tail <= cur_block_start);
        const uint32_t return_to_start = all_later_bytes_consumed && (latest_consume_tail >= &(buf[kDecompressChunkSize]));
        if (return_to_start) {
          // bugfix (2 Oct 2018): This was previously setting
          // syncp->available_end = next_available_end too, and that was being
          // handled as a special case which conflicted with a rare legitimate
          // case.
          syncp->cur_circular_end = next_available_end;
          syncp->available_end = buf;
        } else {
          syncp->available_end = next_available_end;
        }
#ifdef _WIN32
        // bugfix (23 Mar 2018): this needs to be in the critical section,
        // otherwise there's a danger of this resetting legitimate progress
        ResetEvent(consumer_progress_event);
        SetEvent(reader_progress_event);
        LeaveCriticalSection(critical_sectionp);
#else
        // bugfix (21 Mar 2018): must force consumer_progress_state to 0 (or
        // ResetEvent(consumer_progress_event); otherwise the other wait loop's
        // read_stop = buf_end assignment may occur before all later bytes are
        // actually consumed, in the next_available_end == latest_consume_tail
        // edge case.
        syncp->consumer_progress_state = 0;
        pthread_cond_signal(reader_progress_condvarp);
        pthread_mutex_unlock(sync_mutexp);
#endif
        if (return_to_start) {
          // Best to return to the beginning of the buffer.
          // (Note that read_attempt_size is guaranteed to be
          // <= kDecompressChunkSize.)
          const uintptr_t trailing_byte_ct = cur_read_end - next_available_end;
          memcpy(buf, next_available_end, trailing_byte_ct);
          cur_block_start = buf;
          read_head = &(buf[trailing_byte_ct]);
          // May as well reduce false sharing risk.
          read_stop = R_CAST(char*, RoundDownPow2(R_CAST(uintptr_t, latest_consume_tail), kCacheline));
          continue;
        }
        if (all_later_bytes_consumed) {
          read_stop = buf_end;
        } else {
          read_stop = R_CAST(char*, RoundDownPow2(R_CAST(uintptr_t, latest_consume_tail), kCacheline));
        }
        cur_block_start = next_available_end;
      }
      read_head = cur_read_end;
    }
    while (0) {
    TextStreamThread_NOMEM:
      min_interrupt = kTxsInterruptShutdown;
      reterr = kPglRetNomem;
      break;
    TextStreamThread_OPEN_FAIL:
      min_interrupt = kTxsInterruptShutdown;
      syncp->errmsg = strerror(errno);
      reterr = kPglRetOpenFail;
      break;
    TextStreamThread_READ_FAIL:
      min_interrupt = kTxsInterruptShutdown;
      syncp->errmsg = strerror(errno);
      reterr = kPglRetReadFail;
      break;
    TextStreamThread_LONG_LINE:
      min_interrupt = kTxsInterruptShutdown;
      syncp->errmsg = kShortErrLongLine;
      reterr = kPglRetMalformedInput;
      break;
    TextStreamThread_EOF:
      min_interrupt = kTxsInterruptRetarget;
      reterr = kPglRetEof;
      break;
    TextStreamThread_MISC_FAIL:
      min_interrupt = kTxsInterruptShutdown;
      break;
    }
    // We need to wait for a message from the consumer before we can usefully
    // proceed.
    // More precisely:
    // * In the eof subcase, we're waiting for either a rewind or shutdown
    //   request.
    // * In the error subcase, we're just waiting for the shutdown request.

    // Pass the error code back.
#ifdef _WIN32
    EnterCriticalSection(critical_sectionp);
#else
    pthread_mutex_lock(sync_mutexp);
#endif
    syncp->reterr = reterr;
    interrupt = syncp->interrupt;
    if (interrupt >= min_interrupt) {
      // It's our lucky day: we don't need to wait again.
      goto TextStreamThread_INTERRUPT;
    }
    if (reterr == kPglRetEof) {
      syncp->available_end = read_head;
    }
#ifdef _WIN32
    SetEvent(reader_progress_event);
    LeaveCriticalSection(critical_sectionp);
    while (1) {
      WaitForSingleObject(consumer_progress_event, INFINITE);
      EnterCriticalSection(critical_sectionp);
      interrupt = syncp->interrupt;
      if (interrupt >= min_interrupt) {
        break;
      }
      LeaveCriticalSection(critical_sectionp);
    }
#else
    pthread_cond_signal(reader_progress_condvarp);
    do {
      while (!syncp->consumer_progress_state) {
        pthread_cond_wait(consumer_progress_condvarp, sync_mutexp);
      }
      syncp->consumer_progress_state = 0;
      interrupt = syncp->interrupt;
    } while (interrupt < min_interrupt);
#endif
  TextStreamThread_INTERRUPT:
    // must be in critical section here, or be holding the mutex.
    if (interrupt == kTxsInterruptRetarget) {
      new_fname = syncp->new_fname;
      syncp->interrupt = kTxsInterruptNone;
      syncp->reterr = kPglRetSuccess;
    }
#ifdef _WIN32
    LeaveCriticalSection(critical_sectionp);
#else
    pthread_mutex_unlock(sync_mutexp);
#endif
    if (interrupt == kTxsInterruptShutdown) {
      // possible todo: close the file here
      THREAD_RETURN;
    }
    assert(interrupt == kTxsInterruptRetarget);
    read_head = buf;
    if (!new_fname) {
      if (file_type == kFileBgzf) {
        reterr = BgzfRawMtStreamRewind(&rdsp->bgzf, &syncp->errmsg);
        if (unlikely(reterr)) {
          goto TextStreamThread_MISC_FAIL;
        }
      } else {
        // See TextFileRewind().
        rewind(ff);
        if (file_type != kFileUncompressed) {
          if (file_type == kFileGzip) {
            rdsp->gz.ds.avail_in = 0;
#ifdef NDEBUG
            inflateReset(&rdsp->gz.ds);
#else
            const int errcode = inflateReset(&rdsp->gz.ds);
            assert(errcode == Z_OK);
#endif
          } else {
            // kFileZstd
            rdsp->zst.ib.size = 0;
            rdsp->zst.ib.pos = 0;
            ZSTD_DCtx_reset(rdsp->zst.ds, ZSTD_reset_session_only);
          }
        }
      }
    } else {
      // Switch to another file, with less creation/destruction of resources.
      FILE* next_ff = fopen(new_fname, FOPEN_RB);
      if (unlikely(!next_ff)) {
        goto TextStreamThread_OPEN_FAIL;
      }
      // See TextFileOpenInternal().
      uint32_t nbytes = fread_unlocked(buf, 1, 16, next_ff);
      FileCompressionType next_file_type = kFileUncompressed;
      if (nbytes >= 4) {
        const uint32_t magic4 = *R_CAST(uint32_t*, buf);
        if (IsZstdFrame(magic4)) {
          next_file_type = kFileZstd;
        } else if ((magic4 << 8) == 0x088b1f00) {
          if ((nbytes == 16) && IsBgzfHeader(buf)) {
            next_file_type = kFileBgzf;
          } else {
            next_file_type = kFileGzip;
          }
        }
      }
      if (file_type != next_file_type) {
        // Destroy old type-specific resources, and allocate new ones.
        if (file_type == kFileGzip) {
          free(rdsp->gz.in);
          inflateEnd(&rdsp->gz.ds);
        } else if (file_type == kFileBgzf) {
          CleanupBgzfRawMtStream(&rdsp->bgzf);
        } else if (file_type == kFileZstd) {
          free_const(rdsp->zst.ib.src);
          ZSTD_freeDStream(rdsp->zst.ds);
        }

        if (unlikely(fclose(ff))) {
          fclose(next_ff);
          goto TextStreamThread_READ_FAIL;
        }
        ff = next_ff;
        basep->ff = ff;
        file_type = next_file_type;
        basep->file_type = file_type;
        switch (file_type) {
        case kFileUncompressed:
          read_head = &(read_head[nbytes]);
          break;
        case kFileGzip:
          if (unlikely(GzRawInit(buf, nbytes, &rdsp->gz))) {
            goto TextStreamThread_NOMEM;
          }
          break;
        case kFileBgzf:
          reterr = BgzfRawMtStreamInit(buf, context->decompress_thread_ct, ff, nullptr, &rdsp->bgzf, &syncp->errmsg);
          if (unlikely(reterr)) {
            goto TextStreamThread_MISC_FAIL;
          }
          // bugfix (5 Oct 2019): forgot this break
          break;
        case kFileZstd:
          if (unlikely(ZstRawInit(buf, nbytes, &rdsp->zst))) {
            goto TextStreamThread_NOMEM;
          }
          break;
        }
      } else {
        switch (file_type) {
        case kFileUncompressed:
          read_head = &(read_head[nbytes]);
          break;
          // Rest of this is similar to rewind.
        case kFileGzip:
          {
            GzRawDecompressStream* gzp = &rdsp->gz;
            z_stream* dsp = &gzp->ds;
#ifdef NDEBUG
            inflateReset(dsp);
#else
            const int errcode = inflateReset(dsp);
            assert(errcode == Z_OK);
#endif
            memcpy(gzp->in, buf, nbytes);
            dsp->next_in = gzp->in;
            dsp->avail_in = nbytes;
            break;
          }
        case kFileBgzf:
          {
            reterr = BgzfRawMtStreamRetarget(buf, &rdsp->bgzf, next_ff, &syncp->errmsg);
            if (unlikely(reterr)) {
              fclose(next_ff);
              goto TextStreamThread_MISC_FAIL;
            }
            break;
          }
        case kFileZstd:
          {
            ZstRawDecompressStream* zstp = &rdsp->zst;
            ZSTD_DCtx_reset(zstp->ds, ZSTD_reset_session_only);
            memcpy(K_CAST(void*, zstp->ib.src), buf, nbytes);
            zstp->ib.size = nbytes;
            zstp->ib.pos = 0;
            break;
          }
        }
        if (unlikely(fclose(ff))) {
          fclose(next_ff);
          goto TextStreamThread_READ_FAIL;
        }
        ff = next_ff;
        basep->ff = ff;
      }
    }
    cur_block_start = buf;
    read_stop = buf_end;
  }
}

const char kShortErrRfileInvalid[] = "TextStreamOpenEx can't be called with a closed or error-state textFILE";

PglErr TextStreamOpenEx(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, uint32_t decompress_thread_ct, textFILE* txf_ptr, char* dst, TextStream* txs_ptr) {
  TextStreamMain* txsp = GetTxsp(txs_ptr);
  TextFileBase* txs_basep = &txsp->base;
  PglErr reterr = kPglRetSuccess;
  {
    txsp->decompress_thread_ct = decompress_thread_ct;
    if (txf_ptr) {
      // Move-construct (unless there was an error, or file is not opened)
      if (unlikely((!TextFileIsOpen(txf_ptr)) || TextFileErrcode(txf_ptr))) {
        reterr = kPglRetImproperFunctionCall;
        txs_basep->errmsg = kShortErrRfileInvalid;
        goto TextStreamOpenEx_ret_1;
      }
      if (unlikely(TextIsOpen(txs_ptr))) {
        reterr = kPglRetImproperFunctionCall;
        txs_basep->errmsg = kShortErrRfileAlreadyOpen;
        goto TextStreamOpenEx_ret_1;
      }
      textFILEMain* txfp = GetTxfp(txf_ptr);
      *txs_basep = txfp->base;  // struct copy
      // Simplify TextStreamThread() initialization.
      const uint32_t backfill_ct = txs_basep->consume_iter - txs_basep->dst;
      if (backfill_ct) {
        txs_basep->dst_len -= backfill_ct;
        memmove(txs_basep->dst, txs_basep->consume_iter, txs_basep->dst_len);
        txs_basep->consume_iter = txs_basep->dst;
        txs_basep->consume_stop -= backfill_ct;
      }
      txs_basep->enforced_max_line_blen = enforced_max_line_blen;
      assert(txs_basep->dst_len <= dst_capacity);
      txs_basep->dst_capacity = dst_capacity;
      reterr = txfp->base.reterr;
      const FileCompressionType file_type = txfp->base.file_type;
      if (file_type != kFileUncompressed) {
        if (file_type == kFileGzip) {
          txsp->rds.gz = txfp->rds.gz;
        } else if (file_type == kFileZstd) {
          txsp->rds.zst = txfp->rds.zst;
        } else {
          reterr = BgzfRawMtStreamInit(nullptr, decompress_thread_ct, txs_basep->ff, &txfp->rds.bgzf, &txsp->rds.bgzf, &txs_basep->errmsg);
          if (unlikely(reterr)) {
            EraseTextFileBase(&txfp->base);
            goto TextStreamOpenEx_ret_1;
          }
        }
      }
      EraseTextFileBase(&txfp->base);
    } else {
      reterr = TextFileOpenInternal(fname, enforced_max_line_blen, dst_capacity, dst, nullptr, txsp);
    }
    if (reterr) {
      if (reterr == kPglRetEof) {
        txs_basep->reterr = kPglRetEof;
        return kPglRetSuccess;
      }
      goto TextStreamOpenEx_ret_1;
    }
    assert(!txsp->syncp);
    TextStreamSync* syncp;
    if (unlikely(cachealigned_malloc(RoundUpPow2(sizeof(TextStreamSync), kCacheline), &syncp))) {
      goto TextStreamOpenEx_ret_NOMEM;
    }
    txsp->syncp = syncp;
    dst = txs_basep->dst;
    syncp->consume_tail = dst;
    syncp->cur_circular_end = nullptr;
    syncp->available_end = txs_basep->consume_stop;
    syncp->errmsg = nullptr;
    syncp->reterr = kPglRetSuccess;
    syncp->dst_reallocated = 0;
    syncp->interrupt = kTxsInterruptNone;
    syncp->new_fname = nullptr;
#ifdef _WIN32
    syncp->read_thread = nullptr;
    // apparently this can raise a low-memory exception in older Windows
    // versions, but that's not really our problem.
    InitializeCriticalSection(&syncp->critical_section);

    syncp->reader_progress_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (unlikely(!syncp->reader_progress_event)) {
      DeleteCriticalSection(&syncp->critical_section);
      goto TextStreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    syncp->consumer_progress_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (unlikely(!syncp->consumer_progress_event)) {
      DeleteCriticalSection(&syncp->critical_section);
      CloseHandle(syncp->reader_progress_event);
      goto TextStreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    syncp->read_thread = R_CAST(HANDLE, _beginthreadex(nullptr, kDefaultThreadStack, TextStreamThread, txsp, 0, nullptr));
    if (unlikely(!syncp->read_thread)) {
      DeleteCriticalSection(&syncp->critical_section);
      CloseHandle(syncp->consumer_progress_event);
      CloseHandle(syncp->reader_progress_event);
      goto TextStreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
#else
    syncp->sync_init_state = 0;
    if (unlikely(pthread_mutex_init(&syncp->sync_mutex, nullptr))) {
      goto TextStreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    syncp->sync_init_state = 1;
    if (unlikely(pthread_cond_init(&syncp->reader_progress_condvar, nullptr))) {
      goto TextStreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    syncp->sync_init_state = 2;
    syncp->consumer_progress_state = 0;
    if (unlikely(pthread_cond_init(&syncp->consumer_progress_condvar, nullptr))) {
      goto TextStreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    syncp->sync_init_state = 3;
#  ifndef __cplusplus
    pthread_attr_t smallstack_thread_attr;
    if (unlikely(pthread_attr_init(&smallstack_thread_attr))) {
      goto TextStreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    pthread_attr_setstacksize(&smallstack_thread_attr, kDefaultThreadStack);
#  endif
    if (unlikely(pthread_create(&syncp->read_thread,
#  ifdef __cplusplus
                                &g_thread_startup.smallstack_thread_attr,
#  else
                                &smallstack_thread_attr,
#  endif
                                TextStreamThread, txsp))) {
#  ifndef __cplusplus
      pthread_attr_destroy(&smallstack_thread_attr);
#  endif
      goto TextStreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
#  ifndef __cplusplus
    pthread_attr_destroy(&smallstack_thread_attr);
#  endif
    syncp->sync_init_state = 4;
#endif
  }
  while (0) {
  TextStreamOpenEx_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TextStreamOpenEx_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 TextStreamOpenEx_ret_1:
  txs_basep->reterr = reterr;
  return reterr;
}

uint32_t TextDecompressThreadCt(const TextStream* txs_ptr) {
  const TextStreamMain* txsp = GetTxspK(txs_ptr);
  FileCompressionType file_type = txsp->base.file_type;
  if (file_type == kFileUncompressed) {
    return 0;
  }
  if (file_type != kFileBgzf) {
    return 1;
  }
  return GetThreadCtTg(&txsp->rds.bgzf.tg);
}

PglErr TextAdvance(TextStream* txs_ptr) {
  TextStreamMain* txsp = GetTxsp(txs_ptr);
  TextFileBase* basep = &txsp->base;
  char* consume_iter = basep->consume_iter;
  TextStreamSync* syncp = txsp->syncp;
#ifdef _WIN32
  CRITICAL_SECTION* critical_sectionp = &syncp->critical_section;
  HANDLE consumer_progress_event = syncp->consumer_progress_event;
  while (1) {
    EnterCriticalSection(critical_sectionp);
    const PglErr reterr = syncp->reterr;
    if (unlikely((reterr != kPglRetSuccess) && (reterr != kPglRetEof))) {
      basep->errmsg = syncp->errmsg;
      LeaveCriticalSection(critical_sectionp);
      basep->reterr = reterr;
      // No need to set consumer_progress event here, just let the cleanup
      // routine take care of that.
      return reterr;
    }
    char* available_end = syncp->available_end;
    char* cur_circular_end = syncp->cur_circular_end;
    if (consume_iter == cur_circular_end) {
      char* buf = basep->dst;
      consume_iter = buf;
      basep->consume_iter = buf;
      cur_circular_end = nullptr;
      syncp->cur_circular_end = nullptr;
      if (consume_iter != available_end) {
        SetEvent(consumer_progress_event);
      }
    }
    if (syncp->dst_reallocated) {
      consume_iter = basep->dst;
      syncp->dst_reallocated = 0;
    }
    syncp->consume_tail = consume_iter;
    if ((consume_iter != available_end) || cur_circular_end) {
      if (cur_circular_end) {
        basep->consume_stop = cur_circular_end;
      } else {
        basep->consume_stop = available_end;
      }
      LeaveCriticalSection(critical_sectionp);
      // We could set the consumer_progress event here, but it's not really
      // necessary?
      // SetEvent(consumer_progress_event);
      return kPglRetSuccess;
    }
    SetEvent(consumer_progress_event);
    LeaveCriticalSection(critical_sectionp);
    // We've processed all the consume-ready bytes...
    if (reterr != kPglRetSuccess) {
      // ...and we're at eof.  Don't set consumer_progress event here; let that
      // wait until cleanup or rewind/retarget.
      basep->reterr = kPglRetEof;
      return kPglRetEof;
    }
    // ...and there's probably more.
    WaitForSingleObject(syncp->reader_progress_event, INFINITE);
    // bugfix (2 Oct 2018)
    consume_iter = syncp->consume_tail;
    basep->consume_iter = consume_iter;
  }
#else
  pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
  pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
  pthread_cond_t* reader_progress_condvarp = &syncp->reader_progress_condvar;
  pthread_mutex_lock(sync_mutexp);
  while (1) {
    const PglErr reterr = syncp->reterr;
    if (unlikely((reterr != kPglRetSuccess) && (reterr != kPglRetEof))) {
      basep->errmsg = syncp->errmsg;
      pthread_mutex_unlock(sync_mutexp);
      basep->reterr = reterr;
      return reterr;
    }
    char* available_end = syncp->available_end;
    // bugfix (2 Oct 2018): There were TWO consume_iter == available_end ==
    // cur_circular_end cases.
    // printf("checking for more to consume: %lx %lx %lx\n", (uintptr_t)consume_iter, (uintptr_t)syncp->cur_circular_end, (uintptr_t)available_end);
    if (consume_iter == syncp->cur_circular_end) {
      char* buf = basep->dst;
      consume_iter = buf;
      basep->consume_iter = buf;
      syncp->cur_circular_end = nullptr;
      // File-reader could be waiting on either "all bytes in front have been
      // consumed, some bytes behind may remain" or "all bytes have been
      // consumed".  Signal in case it's the first.
      if (consume_iter != available_end) {
        syncp->consumer_progress_state = 1;
        pthread_cond_signal(consumer_progress_condvarp);
      }
    }
    if (syncp->dst_reallocated) {
      consume_iter = basep->dst;
      syncp->dst_reallocated = 0;
    }
    syncp->consume_tail = consume_iter;
    // If cur_circular_end is still non-null here, there must be bytes
    // available even when consume_iter == available_end.  (Is the latter
    // still possible?  Check this.)
    if ((consume_iter != available_end) || syncp->cur_circular_end) {
      if (syncp->cur_circular_end) {
        basep->consume_stop = syncp->cur_circular_end;
      } else {
        basep->consume_stop = available_end;
      }
      // pthread_cond_signal(consumer_progress_condvarp);
      pthread_mutex_unlock(sync_mutexp);
      // printf("consuming %lx..%lx\n", (uintptr_t)(*consume_iterp), (uintptr_t)rlsp->consume_stop);
      return kPglRetSuccess;
    }
    // We've processed all the consume-ready bytes...
    if (reterr != kPglRetSuccess) {
      // ...and we're at eof.
      pthread_mutex_unlock(sync_mutexp);
      basep->reterr = kPglRetEof;
      return kPglRetEof;
    }
    // ...and there's probably more.
    syncp->consumer_progress_state = 1;
    pthread_cond_signal(consumer_progress_condvarp);
    // no need for an explicit spurious-wakeup check, we'll check the progress
    // condition (available_end has advanced, or we have a read error) anyway
    // and get back here if it isn't satisfied
    pthread_cond_wait(reader_progress_condvarp, sync_mutexp);
    // bugfix (2 Oct 2018)
    consume_iter = syncp->consume_tail;
    basep->consume_iter = syncp->consume_tail;
  }
#endif
}

PglErr TextOnlyEmptyLinesLeft(TextStream* txs_ptr) {
  TextFileBase* basep = &GetTxsp(txs_ptr)->base;
  char* line_start = basep->consume_iter;
  while (1) {
    if (line_start == basep->consume_stop) {
      basep->consume_iter = line_start;
      PglErr reterr = TextAdvance(txs_ptr);
      if (reterr) {
        return reterr;
      }
      line_start = basep->consume_iter;
    }
    line_start = FirstNonTspace(line_start);
    if (unlikely(!IsEolnKns(*line_start))) {
      basep->reterr = kPglRetMalformedInput;
      basep->errmsg = kShortErrInteriorEmptyLine;
      return kPglRetMalformedInput;
    }
    line_start = AdvPastDelim(line_start, '\n');
  }
}

PglErr TextSkipNz(uintptr_t skip_ct, TextStream* txs_ptr) {
  TextFileBase* basep = &GetTxsp(txs_ptr)->base;
#ifdef __LP64__
  char* consume_iter = basep->consume_iter;
  // Minor extension of AdvToNthDelimChecked().
  const VecUc vvec_all_lf = vecuc_set1('\n');
  while (1) {
    uintptr_t starting_addr = R_CAST(uintptr_t, consume_iter);
    VecUc* consume_viter = R_CAST(VecUc*, RoundDownPow2(starting_addr, kBytesPerVec));
    uintptr_t ending_addr = R_CAST(uintptr_t, basep->consume_stop);
    VecUc* consume_vstop = R_CAST(VecUc*, RoundDownPow2(ending_addr, kBytesPerVec));
    VecUc cur_vvec = *consume_viter;
    VecUc lf_vvec = (cur_vvec == vvec_all_lf);
    uint32_t lf_bytes = vecuc_movemask(lf_vvec);
    const uint32_t leading_byte_ct = starting_addr - R_CAST(uintptr_t, consume_viter);
    const uint32_t leading_mask = UINT32_MAX << leading_byte_ct;
    lf_bytes &= leading_mask;
    uint32_t cur_lf_ct;
    for (; consume_viter != consume_vstop; ) {
      cur_lf_ct = PopcountVec8thUint(lf_bytes);
      if (cur_lf_ct >= skip_ct) {
        goto TextSkipNz_finish;
      }
      skip_ct -= cur_lf_ct;
      // bugfix (28 Sep 2019): forgot to update cur_vvec/lf_vvec/lf_bytes?!
      ++consume_viter;
      cur_vvec = *consume_viter;
      lf_vvec = (cur_vvec == vvec_all_lf);
      lf_bytes = vecuc_movemask(lf_vvec);
    }
    lf_bytes &= (1U << (ending_addr % kBytesPerVec)) - 1;
    cur_lf_ct = PopcountVec8thUint(lf_bytes);
    if (cur_lf_ct >= skip_ct) {
    TextSkipNz_finish:
      lf_bytes = ClearBottomSetBits(skip_ct - 1, lf_bytes);
      const uint32_t byte_offset_in_vec = ctzu32(lf_bytes) + 1;
      const uintptr_t result_addr = R_CAST(uintptr_t, consume_viter) + byte_offset_in_vec;
      basep->consume_iter = R_CAST(char*, result_addr);
      return kPglRetSuccess;
    }
    skip_ct -= cur_lf_ct;
    // bugfix (30 Oct 2019)
    basep->consume_iter = basep->consume_stop;
    PglErr reterr = TextAdvance(txs_ptr);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
    consume_iter = basep->consume_iter;
  }
#else
  char* consume_iter = basep->consume_iter;
  char* consume_stop = basep->consume_stop;
  for (uintptr_t ulii = 0; ulii != skip_ct; ++ulii) {
    if (consume_iter == consume_stop) {
      basep->consume_iter = consume_iter;
      PglErr reterr = TextAdvance(txs_ptr);
      if (reterr) {
        return reterr;
      }
      consume_iter = basep->consume_iter;
      consume_stop = basep->consume_stop;
    }
    consume_iter = AdvPastDelim(consume_iter, '\n');
  }
  basep->consume_iter = consume_iter;
  return kPglRetSuccess;
#endif
}

PglErr TextRetarget(const char* new_fname, TextStream* txs_ptr) {
  TextStreamMain* txsp = GetTxsp(txs_ptr);
  TextFileBase* basep = &txsp->base;
  TextStreamSync* syncp = txsp->syncp;
#ifdef _WIN32
  CRITICAL_SECTION* critical_sectionp = &syncp->critical_section;
  EnterCriticalSection(critical_sectionp);
  const PglErr reterr = syncp->reterr;
  if (reterr != kPglRetSuccess) {
    if (unlikely(reterr != kPglRetEof)) {
      basep->errmsg = syncp->errmsg;
      LeaveCriticalSection(critical_sectionp);
      basep->reterr = reterr;
      return reterr;
    }
    // clear eof
    syncp->reterr = kPglRetSuccess;
  }
  basep->reterr = kPglRetSuccess;
  // bugfix (5 Mar 2018): need to reset these here, can't wait for reader
  // thread to receive signal
  char* buf = basep->dst;
  syncp->consume_tail = buf;
  syncp->cur_circular_end = nullptr;
  syncp->available_end = buf;
  syncp->dst_reallocated = 0;
  syncp->interrupt = kTxsInterruptRetarget;
  // Could also just open the file in this function (before acquiring the
  // mutex) and pass a gzFile.  Advantages: nothing bad happens if new_fname
  // is overwritten before it's read, RLstreamErrPrint() no longer has to deal
  // with OpenFail error.  Disadvantage: peak resource usage is a bit higher if
  // we open the second file before closing the first one.  Advantages probably
  // outweigh disadvantages, but I'll wait till --pmerge development to make a
  // decision since that's the main function that actually cares.
  syncp->new_fname = new_fname;
  SetEvent(syncp->consumer_progress_event);
  LeaveCriticalSection(critical_sectionp);
#else
  pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
  pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
  pthread_mutex_lock(sync_mutexp);
  const PglErr reterr = syncp->reterr;
  if (reterr != kPglRetSuccess) {
    if (unlikely(reterr != kPglRetEof)) {
      basep->errmsg = syncp->errmsg;
      pthread_mutex_unlock(sync_mutexp);
      basep->reterr = reterr;
      return reterr;
    }
    // clear eof
    syncp->reterr = kPglRetSuccess;
  }
  // bugfix (4 Oct 2019): also need to clear eof here
  basep->reterr = kPglRetSuccess;
  char* buf = basep->dst;
  syncp->consume_tail = buf;
  syncp->cur_circular_end = nullptr;
  syncp->available_end = buf;
  syncp->dst_reallocated = 0;
  syncp->interrupt = kTxsInterruptRetarget;
  syncp->new_fname = new_fname;
  syncp->consumer_progress_state = 1;
  pthread_cond_signal(consumer_progress_condvarp);
  pthread_mutex_unlock(sync_mutexp);
#endif
  basep->consume_iter = buf;
  basep->consume_stop = buf;
  return kPglRetSuccess;
}

BoolErr CleanupTextStream(TextStream* txs_ptr, PglErr* reterrp) {
  TextStreamMain* txsp = GetTxsp(txs_ptr);
  TextFileBase* basep = &txsp->base;
  TextStreamSync* syncp = txsp->syncp;
  if (syncp) {
#ifdef _WIN32
    if (syncp->read_thread) {
      CRITICAL_SECTION* critical_sectionp = &syncp->critical_section;
      EnterCriticalSection(critical_sectionp);
      syncp->interrupt = kTxsInterruptShutdown;
      SetEvent(syncp->consumer_progress_event);
      LeaveCriticalSection(critical_sectionp);
      WaitForSingleObject(syncp->read_thread, INFINITE);
      DeleteCriticalSection(critical_sectionp);
      CloseHandle(syncp->consumer_progress_event);
      CloseHandle(syncp->reader_progress_event);
    }
#else
    const uint32_t sync_init_state = syncp->sync_init_state;
    if (sync_init_state) {
      pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
      pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
      if (sync_init_state == 4) {
        pthread_mutex_lock(sync_mutexp);
        syncp->interrupt = kTxsInterruptShutdown;
        syncp->consumer_progress_state = 1;
        pthread_cond_signal(consumer_progress_condvarp);
        pthread_mutex_unlock(sync_mutexp);
        pthread_join(syncp->read_thread, nullptr);
      }
      pthread_mutex_destroy(sync_mutexp);
      if (sync_init_state > 1) {
        pthread_cond_destroy(&syncp->reader_progress_condvar);
        if (sync_init_state > 2) {
          pthread_cond_destroy(consumer_progress_condvarp);
        }
      }
    }
#endif
    aligned_free(txsp->syncp);
    txsp->syncp = nullptr;
  }
  basep->consume_iter = nullptr;
  basep->consume_stop = nullptr;
  basep->reterr = kPglRetEof;
  basep->errmsg = nullptr;
  if (basep->dst && (!basep->dst_owned_by_consumer)) {
    free(basep->dst);
    basep->dst = nullptr;
  }
  if (basep->ff) {
    if (basep->file_type != kFileUncompressed) {
      if (basep->file_type == kFileZstd) {
        if (txsp->rds.zst.ib.src) {
          free_const(txsp->rds.zst.ib.src);
          txsp->rds.zst.ib.src = nullptr;
        }
        if (txsp->rds.zst.ds) {
          ZSTD_freeDStream(txsp->rds.zst.ds);
          txsp->rds.zst.ds = nullptr;
        }
      } else if (basep->file_type == kFileBgzf) {
        CleanupBgzfRawMtStream(&txsp->rds.bgzf);
      } else {
        // plain gzip
        if (txsp->rds.gz.in) {
          free(txsp->rds.gz.in);
          txsp->rds.gz.in = nullptr;
        }
        if (txsp->rds.gz.ds_initialized) {
          inflateEnd(&txsp->rds.gz.ds);
        }
      }
      basep->file_type = kFileUncompressed;
    }
    if (unlikely(fclose_null(&basep->ff))) {
      if (!reterrp) {
        return 1;
      }
      if (*reterrp == kPglRetSuccess) {
        // Note that we don't set basep->reterr or ->errmsg here.
        *reterrp = kPglRetReadFail;
        return 1;
      }
    }
  }
  return 0;
}


PglErr TksNext(TokenStream* tksp, uint32_t shard_ct, char** shard_boundaries) {
  TextStreamMain* txsp = GetTxsp(&tksp->txs);
  txsp->base.consume_iter = txsp->base.consume_stop;
  PglErr reterr = TextAdvance(&(tksp->txs));
  if (reterr) { // not unlikely due to eof
    return reterr;
  }
  char* consume_iter = txsp->base.consume_iter;
  char* consume_stop = txsp->base.consume_stop;
  shard_boundaries[0] = consume_iter;
  shard_boundaries[shard_ct] = consume_stop;
  if (shard_ct > 1) {
    const uintptr_t shard_size_target = S_CAST(uintptr_t, consume_stop - consume_iter) / shard_ct;
    char* boundary_min = consume_iter;
    char* cur_boundary = consume_iter;
    for (uint32_t boundary_idx = 1; boundary_idx < shard_ct; ++boundary_idx) {
      boundary_min = &(boundary_min[shard_size_target]);
      if (boundary_min > cur_boundary) {
        // last character must be token separator
        cur_boundary = FirstSpaceOrEoln(boundary_min);
        ++cur_boundary;
      }
      shard_boundaries[boundary_idx] = cur_boundary;
    }
  }
  return kPglRetSuccess;
}

#ifdef __cplusplus
}
#endif
