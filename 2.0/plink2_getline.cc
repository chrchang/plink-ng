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

#include <errno.h>
#include "plink2_getline.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void EraseTextRfileBase(TextRfileBase* trbp) {
  trbp->consume_iter = nullptr;
  trbp->consume_stop = nullptr;
  trbp->errmsg = nullptr;
  trbp->reterr = kPglRetEof;
  trbp->ff = nullptr;
  trbp->dst = nullptr;
}

void PreinitTextRfile(textRFILE* trfp) {
  EraseTextRfileBase(&trfp->base);
}

const char kShortErrRfileAlreadyOpen[] = "TextRfileOpenInternal can't be called on an already-open file";
const char kShortErrRfileEnforcedMaxBlenTooSmall[] = "TextRfileOpenInternal: enforced_max_line_blen too small (must be at least max(1 MiB, dst_capacity - 1 MiB))";
const char kShortErrRfileDstCapacityTooSmall[] = "TextRfileOpenInternal: dst_capacity too small (2 MiB minimum)";

PglErr TextRfileOpenInternal(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, uint32_t decompress_thread_ct, char* dst, textRFILE* trfp, TextRstream* trsp) {
  PglErr reterr = kPglRetSuccess;
  TextRfileBase* trbp;
  if (trfp) {
    trbp = &trfp->base;
  } else {
    trbp = &trsp->base;
  }
  {
    // 1. Open file, get type.
    if (unlikely(trbp->ff)) {
      reterr = kPglRetImproperFunctionCall;
      trbp->errmsg = kShortErrRfileAlreadyOpen;
      goto TextRfileOpenInternal_ret_1;
    }
    if (enforced_max_line_blen || trfp) {
      if (unlikely(enforced_max_line_blen < kDecompressChunkSize)) {
        reterr = kPglRetImproperFunctionCall;
        trbp->errmsg = kShortErrRfileEnforcedMaxBlenTooSmall;
        goto TextRfileOpenInternal_ret_1;
      }
      if (dst) {
        if (unlikely(dst_capacity < 2 * kDecompressChunkSize)) {
          reterr = kPglRetImproperFunctionCall;
          trbp->errmsg = kShortErrRfileDstCapacityTooSmall;
          goto TextRfileOpenInternal_ret_1;
        }
        if (unlikely(enforced_max_line_blen + kDecompressChunkSize < dst_capacity)) {
          reterr = kPglRetImproperFunctionCall;
          trbp->errmsg = kShortErrRfileEnforcedMaxBlenTooSmall;
          goto TextRfileOpenInternal_ret_1;
        }
      }
    } else {
      // token-reading mode.  dst == nullptr not currently supported.
      assert(dst && (dst_capacity == kTokenRstreamBlen));
    }
    trbp->ff = fopen(fname, FOPEN_RB);
    if (unlikely(!trbp->ff)) {
      goto TextRfileOpenInternal_ret_OPEN_FAIL;
    }
    trbp->file_type = kFileUncompressed;
    if (dst) {
      trbp->dst_owned_by_caller = 1;
      trbp->dst_capacity = dst_capacity;
    } else {
      dst = S_CAST(char*, malloc(2 * kDecompressChunkSize));
      if (unlikely(dst == nullptr)) {
        goto TextRfileOpenInternal_ret_NOMEM;
      }
      trbp->dst_owned_by_caller = 0;
      trbp->dst_capacity = 2 * kDecompressChunkSize;
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
        if (trfp) {
          zstp = &trfp->rds.zst;
        } else {
          zstp = &trsp->rds.zst;
        }
        zstp->ib.src = malloc(kDecompressChunkSize);
        if (unlikely(!zstp->ib.src)) {
          zstp->ds = nullptr;
          goto TextRfileOpenInternal_ret_NOMEM;
        }
        zstp->ds = ZSTD_createDStream();
        if (unlikely(!zstp->ds)) {
          goto TextRfileOpenInternal_ret_NOMEM;
        }
        memcpy(K_CAST(void*, zstp->ib.src), dst, nbytes);
        zstp->ib.size = nbytes;
        zstp->ib.pos = 0;
      } else if ((magic4 << 8) == 0x088b1f00) {
        // gzip ID1/ID2 bytes, deflate compression method
        trbp->dst_len = 0;
        if ((nbytes == 16) && IsBgzfHeader(dst)) {
          trbp->file_type = kFileBgzf;
          if (trfp) {
            BgzfRawDecompressStream* bgzfp = &trfp->rds.bgzf;
            bgzfp->in = S_CAST(unsigned char*, malloc(kDecompressChunkSize));
            if (unlikely(!bgzfp->in)) {
              bgzfp->ldc = nullptr;
              goto TextRfileOpenInternal_ret_NOMEM;
            }
            bgzfp->ldc = libdeflate_alloc_decompressor();
            if (!bgzfp->ldc) {
              goto TextRfileOpenInternal_ret_NOMEM;
            }
            memcpy(bgzfp->in, dst, nbytes);
            bgzfp->in_size = nbytes;
            bgzfp->in_pos = 0;
          } else {
            // TODO
            assert(decompress_thread_ct != 0);
            exit(1);
          }
        } else {
          trbp->file_type = kFileGzip;
          GzRawDecompressStream* gzp;
          if (trfp) {
            gzp = &trfp->rds.gz;
          } else {
            gzp = &trsp->rds.gz;
          }
          gzp->ds_initialized = 0;
          gzp->in = S_CAST(unsigned char*, malloc(kDecompressChunkSize));
          if (!gzp->in) {
            goto TextRfileOpenInternal_ret_NOMEM;
          }
          z_stream* dsp = &gzp->ds;
          memcpy(gzp->in, dst, nbytes);
          dsp->next_in = gzp->in;
          dsp->avail_in = nbytes;
          dsp->zalloc = nullptr;
          dsp->zfree = nullptr;
          dsp->opaque = nullptr;
          if (unlikely(inflateInit2(dsp, MAX_WBITS | 16) != Z_OK)) {
            // todo: verify no other errors possible
            goto TextRfileOpenInternal_ret_NOMEM;
          }
          gzp->ds_initialized = 1;
        }
      }
    } else if (!nbytes) {
      if (unlikely(!feof_unlocked(trbp->ff))) {
        goto TextRfileOpenInternal_ret_READ_FAIL;
      }
      // May as well accept this.
      // Don't jump to ret_1 since we're setting trfp->reterr to a different
      // value than we're returning.
      trbp->reterr = kPglRetEof;
      return kPglRetSuccess;
    }
  }
  while (0) {
  TextRfileOpenInternal_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TextRfileOpenInternal_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    trbp->errmsg = strerror(errno);
    break;
  TextRfileOpenInternal_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    trbp->errmsg = strerror(errno);
    break;
  }
 TextRfileOpenInternal_ret_1:
  trbp->reterr = reterr;
  return reterr;
}

PglErr TextRfileOpenEx(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, char* dst, textRFILE* trfp) {
  return TextRfileOpenInternal(fname, enforced_max_line_blen, dst_capacity, 0, dst, trfp, nullptr);
}

// Set enforced_max_line_blen == 0 in the token-reading case.
// trailing 'X' is temporary, to avoid duplicate-symbol error
BoolErr IsPathologicallyLongLineOrTokenX(const char* line_start, const char* load_start, const char* known_line_end, uint32_t enforced_max_line_blen) {
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

const char kShortErrInvalidBgzf[] = "Malformed BGZF block";
const char kShortErrRfileTruncatedGz[] = "TextRfileAdvance: gzipped file appears to be truncated";
const char kShortErrLongLine[] = "Pathologically long line";

PglErr TextRfileAdvance(textRFILE* trfp) {
  if (trfp->base.reterr) {
    return trfp->base.reterr;
  }
  PglErr reterr = kPglRetSuccess;
  {
    char* orig_line_start = trfp->base.consume_stop;
    assert(trfp->base.consume_iter == orig_line_start);
    char* dst = trfp->base.dst;
    char* dst_load_start;
    while (1) {
      const uint32_t dst_offset = orig_line_start - dst;
      const uint32_t dst_rem = trfp->base.dst_len - dst_offset;
      // (dst_rem guaranteed to be < trfp->base.enforced_max_line_blen here,
      // since otherwise we error out earlier.)
      // Two cases:
      // 1. Move (possibly empty) unfinished line to the beginning of the
      //    buffer.
      // 2. Resize the buffer/report out-of-memory.
      if (dst_rem < trfp->base.dst_capacity - kDecompressChunkSize) {
        memmove(dst, orig_line_start, dst_rem);
      } else {
        if (unlikely(trfp->base.dst_owned_by_caller)) {
          goto TextRfileAdvance_ret_NOMEM;
        }
        uint32_t next_dst_capacity = trfp->base.enforced_max_line_blen + kDecompressChunkSize;
        if ((next_dst_capacity / 2) > trfp->base.dst_capacity) {
          next_dst_capacity = trfp->base.dst_capacity * 2;
        }
#ifndef __LP64__
        if (next_dst_capacity >= 0x80000000U) {
          goto TextRfileAdvance_ret_NOMEM;
        }
#endif
        char* dst_next;
        if (!dst_offset) {
          dst_next = S_CAST(char*, realloc(dst, next_dst_capacity));
          if (unlikely(!dst_next)) {
            goto TextRfileAdvance_ret_NOMEM;
          }
        } else {
          dst_next = S_CAST(char*, malloc(next_dst_capacity));
          if (unlikely(!dst_next)) {
            goto TextRfileAdvance_ret_NOMEM;
          }
          memcpy(dst_next, orig_line_start, dst_rem);
        }
        trfp->base.dst = dst_next;
        dst = dst_next;
      }
      dst_load_start = &(dst[dst_rem]);
      char* dst_iter = dst_load_start;
      char* dst_end = &(dst[trfp->base.dst_capacity]);
      trfp->base.consume_iter = dst;
      FILE* ff = trfp->base.ff;
      switch (trfp->base.file_type) {
      case kFileUncompressed:
        {
          uint32_t rlen = dst_end - dst_iter;
          if (rlen > kMaxBytesPerIO) {
            // We need to know how many bytes were read, so fread_checked()
            // doesn't work.
            // This is an if-statement instead of a while loop since rlen can
            // never be larger than 2 * kMaxBytesPerIO.
            const uint32_t nbytes = fread_unlocked(dst_iter, 1, kMaxBytesPerIO, ff);
            if (nbytes < kMaxBytesPerIO) {
              if (unlikely(ferror_unlocked(ff))) {
                goto TextRfileAdvance_ret_READ_FAIL;
              }
              trfp->base.dst_len = nbytes + dst_rem;
              break;
            }
            rlen -= kMaxBytesPerIO;
            dst_iter = &(dst_iter[kMaxBytesPerIO]);
          }
          const uint32_t nbytes = fread_unlocked(dst_iter, 1, rlen, ff);
          if (unlikely(ferror_unlocked(ff))) {
            goto TextRfileAdvance_ret_READ_FAIL;
          }
          dst_iter = &(dst_iter[nbytes]);
          break;
        }
      case kFileGzip:
        {
          z_stream* dsp = &trfp->rds.gz.ds;
          if ((!dsp->avail_in) && feof_unlocked(ff)) {
            break;
          }
          do {
            int zerr = Z_OK;
            if (dsp->avail_in) {  // can be zero after TextRewind()
              dsp->next_out = R_CAST(unsigned char*, dst_iter);
              dsp->avail_out = dst_end - dst_iter;
              zerr = inflate(dsp, Z_SYNC_FLUSH);
              if (unlikely((zerr < 0) || (zerr == Z_NEED_DICT))) {
                if (dsp->msg) {
                  trfp->base.errmsg = dsp->msg;
                } else {
                  trfp->base.errmsg = zError(zerr);
                }
                goto TextRfileAdvance_ret_DECOMPRESS_FAIL;
              }
              dst_iter = R_CAST(char*, dsp->next_out);
              if (dsp->avail_in) {
                assert(dst_iter == dst_end);
                break;
              }
            }
            const uint32_t nbytes = fread_unlocked(trfp->rds.gz.in, 1, kDecompressChunkSize, ff);
            dsp->next_in = trfp->rds.gz.in;
            dsp->avail_in = nbytes;
            if (!nbytes) {
              if (unlikely(!feof_unlocked(ff))) {
                goto TextRfileAdvance_ret_READ_FAIL;
              }
              if (unlikely(zerr == Z_OK)) {
                trfp->base.errmsg = kShortErrRfileTruncatedGz;
                goto TextRfileAdvance_ret_DECOMPRESS_FAIL;
              }
              // Normal EOF.
              break;
            }
          } while (dst_iter != dst_end);
          break;
        }
      case kFileBgzf:
        {
          // Fully independent blocks limited to 64 KiB.
          if ((!trfp->rds.bgzf.in_size) && feof_unlocked(ff)) {
            break;
          }
          struct libdeflate_decompressor* ldc = trfp->rds.bgzf.ldc;
          unsigned char* in = trfp->rds.bgzf.in;
          unsigned char* in_iter = &(in[trfp->rds.bgzf.in_pos]);
          unsigned char* in_end = &(in[trfp->rds.bgzf.in_size]);
          while (1) {
            uint32_t n_inbytes = in_end - in_iter;
            if (n_inbytes > 25) {
              if (unlikely(!IsBgzfHeader(in_iter))) {
                goto TextRfileAdvance_ret_INVALID_BGZF;
              }
              const uint32_t bsize_minus1 = *R_CAST(uint16_t*, &(in_iter[16]));
              if (unlikely(bsize_minus1 < 25)) {
                goto TextRfileAdvance_ret_INVALID_BGZF;
              }
              if (bsize_minus1 < n_inbytes) {
                // We have at least one fully-loaded compressed block.
                // Decompress it if we have enough space.
                const uint32_t in_size = bsize_minus1 - 25;
#  ifdef __arm__
#    error "Unaligned accesses in TextRfileAdvance()."
#  endif
                const uint32_t out_size = *R_CAST(uint32_t*, &(in_iter[in_size + 22]));
                if (unlikely(out_size > 65536)) {
                  goto TextRfileAdvance_ret_INVALID_BGZF;
                }
                if (out_size > S_CAST(uintptr_t, dst_end - dst_iter)) {
                  break;
                }
                if (unlikely(libdeflate_deflate_decompress(ldc, &(in_iter[18]), in_size, dst_iter, out_size, nullptr))) {
                  goto TextRfileAdvance_ret_INVALID_BGZF;
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
              goto TextRfileAdvance_ret_READ_FAIL;
            }
            in_iter = in;
            in_end = &(load_start[nbytes]);
            trfp->rds.bgzf.in_size = in_end - in;
            if (!nbytes) {
              if (unlikely(n_inbytes)) {
                goto TextRfileAdvance_ret_INVALID_BGZF;
              }
              break;
            }
          }
          trfp->rds.bgzf.in_pos = in_iter - in;
          dst_end = dst_iter;
          break;
        }
      case kFileZstd:
        {
          if ((!trfp->rds.zst.ib.size) && feof_unlocked(ff)) {
            break;
          }
          // Sequentially dependent blocks limited to ~128 KiB.
          while (1) {
            ZSTD_outBuffer zob = {R_CAST(unsigned char*, dst_iter), S_CAST(size_t, dst_end - dst_iter), 0};
            // ib.size == 0 ok, no need to special-case rewind.
            const uintptr_t read_size_hint = ZSTD_decompressStream(trfp->rds.zst.ds, &zob, &trfp->rds.zst.ib);
            if (unlikely(ZSTD_isError(read_size_hint))) {
              trfp->base.errmsg = ZSTD_getErrorName(read_size_hint);
              goto TextRfileAdvance_ret_DECOMPRESS_FAIL;
            }
            dst_iter = &(dst_iter[zob.pos]);
            if (dst_iter == dst_end) {
              break;
            }
            // Decoder has flushed everything it could.  Either we're at EOF,
            // or we must load more.
            unsigned char* in = S_CAST(unsigned char*, K_CAST(void*, trfp->rds.zst.ib.src));
            const uint32_t n_inbytes = trfp->rds.zst.ib.size - trfp->rds.zst.ib.pos;
            memmove(in, &(in[trfp->rds.zst.ib.pos]), n_inbytes);
            unsigned char* load_start = &(in[n_inbytes]);
            const uint32_t nbytes = fread_unlocked(load_start, 1, kDecompressChunkSize - n_inbytes, ff);
            if (unlikely(ferror_unlocked(ff))) {
              goto TextRfileAdvance_ret_READ_FAIL;
            }
            trfp->rds.zst.ib.pos = 0;
            trfp->rds.zst.ib.size = nbytes + n_inbytes;
            if (!nbytes) {
              if (unlikely(n_inbytes)) {
                trfp->base.errmsg = kShortErrZstdPrefixUnknown;
                goto TextRfileAdvance_ret_DECOMPRESS_FAIL;
              }
              break;
            }
          }
          break;
        }
      }
      trfp->base.dst_len = dst_iter - dst;
      if (!trfp->base.dst_len) {
        goto TextRfileAdvance_ret_EOF;
      }
      if (dst_iter != dst_end) {
        // If last character of file isn't a newline, append one to simplify
        // downstream code.
        if (dst_iter[-1] != '\n') {
          *dst_iter++ = '\n';
          trfp->base.dst_len += 1;
        }
        trfp->base.consume_stop = dst_iter;
        break;
      }
      char* last_byte_ptr = Memrchr(dst_load_start, '\n', dst_iter - dst_load_start);
      if (last_byte_ptr) {
        trfp->base.consume_stop = &(last_byte_ptr[1]);
        break;
      }
      // Buffer is full, and no '\n' is present.  Restart the loop and try to
      // extend the buffer, if we aren't already at/past the line-length limit.
      if (trfp->base.dst_len >= trfp->base.enforced_max_line_blen) {
        goto TextRfileAdvance_ret_LONG_LINE;
      }
    }
    if (unlikely(IsPathologicallyLongLineOrTokenX(dst, dst_load_start, trfp->base.consume_stop, trfp->base.enforced_max_line_blen))) {
      goto TextRfileAdvance_ret_LONG_LINE;
    }
  }
  while (0) {
  TextRfileAdvance_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TextRfileAdvance_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    trfp->base.errmsg = strerror(errno);
    break;
  TextRfileAdvance_ret_LONG_LINE:
    trfp->base.errmsg = kShortErrLongLine;
    reterr = kPglRetMalformedInput;
    break;
  TextRfileAdvance_ret_INVALID_BGZF:
    trfp->base.errmsg = kShortErrInvalidBgzf;
  TextRfileAdvance_ret_DECOMPRESS_FAIL:
    reterr = kPglRetDecompressFail;
    break;
  TextRfileAdvance_ret_EOF:
    reterr = kPglRetEof;
    break;
  }
  trfp->base.reterr = reterr;
  return reterr;
}

void TextRfileRewind(textRFILE* trfp) {
  if ((!trfp->base.ff) || ((trfp->base.reterr) && (trfp->base.reterr != kPglRetEof))) {
    return;
  }
  rewind(trfp->base.ff);
  trfp->base.reterr = kPglRetSuccess;
  trfp->base.dst_len = 0;
  trfp->base.consume_iter = trfp->base.dst;
  trfp->base.consume_stop = trfp->base.dst;
  if (trfp->base.file_type != kFileUncompressed) {
    if (trfp->base.file_type == kFileGzip) {
      trfp->rds.gz.ds.avail_in = 0;
#ifdef NDEBUG
      inflateReset(&trfp->rds.gz.ds);
#else
      const int errcode = inflateReset(&trfp->rds.gz.ds);
      assert(errcode == Z_OK);
#endif
    } else if (trfp->base.file_type == kFileBgzf) {
      trfp->rds.bgzf.in_size = 0;
      trfp->rds.bgzf.in_pos = 0;
    } else {
      // kFileZstd
      trfp->rds.zst.ib.size = 0;
      trfp->rds.zst.ib.pos = 0;
      ZSTD_DCtx_reset(trfp->rds.zst.ds, ZSTD_reset_session_only);
    }
  }
}

BoolErr CleanupTextRfile(textRFILE* trfp, PglErr* reterrp) {
  trfp->base.consume_iter = nullptr;
  trfp->base.consume_stop = nullptr;
  trfp->base.reterr = kPglRetEof;
  trfp->base.errmsg = nullptr;
  if (trfp->base.dst && (!trfp->base.dst_owned_by_caller)) {
    free(trfp->base.dst);
    trfp->base.dst = nullptr;
  }
  if (trfp->base.ff) {
    if (trfp->base.file_type != kFileUncompressed) {
      if (trfp->base.file_type == kFileZstd) {
        if (trfp->rds.zst.ib.src) {
          free_const(trfp->rds.zst.ib.src);
          trfp->rds.zst.ib.src = nullptr;
        }
        if (trfp->rds.zst.ds) {
          ZSTD_freeDStream(trfp->rds.zst.ds);
          trfp->rds.zst.ds = nullptr;
        }
      } else if (trfp->base.file_type == kFileBgzf) {
        if (trfp->rds.bgzf.in) {
          free(trfp->rds.bgzf.in);
          trfp->rds.bgzf.in = nullptr;
        }
        if (trfp->rds.bgzf.ldc) {
          libdeflate_free_decompressor(trfp->rds.bgzf.ldc);
          trfp->rds.bgzf.ldc = nullptr;
        }
      } else {
        // plain gzip
        if (trfp->rds.gz.in) {
          free(trfp->rds.gz.in);
          trfp->rds.gz.in = nullptr;
        }
        if (trfp->rds.gz.ds_initialized) {
          inflateEnd(&trfp->rds.gz.ds);
        }
      }
    }
    if (unlikely(fclose_null(&trfp->base.ff))) {
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


void PreinitTextRstream(TextRstream* trsp) {
  EraseTextRfileBase(&trsp->base);
  trsp->syncp = nullptr;
#ifdef _WIN32
  trsp->read_thread = nullptr;
#else
  trsp->sync_init_state = 0;
#endif
}

// This type of code is especially bug-prone (ESR would call it a "defect
// attractor").  Goal is to get it right, and fast enough to be a major win
// over gzgets()... and then not worry about it again for years.
THREAD_FUNC_DECL TextRstreamThread(void* arg) {
  THREAD_RETURN;
  /*
  TextRstream* context = S_CAST(TextRstream*, arg);
  TextRstreamSync* syncp = context->syncp;
  FileCompressionType file_type = context->base.file_type;
  FILE* ff = context->base.ff;
  RawMtDecompressStream* rdsp = &context->rds;
  char* buf = context->base.dst;
  char* buf_end = &(buf[context->base.dst_capacity]);
  char* cur_block_start = buf;
  char* read_head = buf;

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
  const uint32_t enforced_max_line_blen = context->base.enforced_max_line_blen;
  const char* new_fname = nullptr;
  const uint32_t is_token_stream = (enforced_max_line_blen == 0);
  while (1) {
    TrsInterrupt interrupt = kTrsInterruptNone;
    PglErr reterr;
    TrsInterrupt min_interrupt;
    while (1) {
      uintptr_t read_attempt_size = read_stop - read_head;
      if (!read_attempt_size) {
        const uint32_t memmove_required = (read_stop == buf_end);
        if (unlikely((cur_block_start == buf) && memmove_required)) {
          // TODO: if !dst_owned_by_caller, try to expand buffer
          if (!context->base.dst_owned_by_caller) {
            exit(1);
          }
          goto TextRstreamThread_LONG_LINE;
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
          goto TextRstreamThread_wait_first;
        }
        while (1) {
          EnterCriticalSection(critical_sectionp);
          interrupt = syncp->interrupt;
          if (interrupt != kTrsInterruptNone) {
            goto TextRstreamThread_INTERRUPT;
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
        TextRstreamThread_wait_first:
          WaitForSingleObject(consumer_progress_event, INFINITE);
        }
        // bugfix (23 Mar 2018): didn't always leave the critical section
        LeaveCriticalSection(critical_sectionp);
#else
        pthread_mutex_lock(sync_mutexp);
        if (!memmove_required) {
          // Wait for all bytes in front of read_stop to be consumed.
          goto TextRstreamThread_wait_first;
        }
        while (1) {
          interrupt = syncp->interrupt;
          if (interrupt != kTrsInterruptNone) {
            goto TextRstreamThread_INTERRUPT;
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
        TextRstreamThread_wait_first:
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
      int32_t bytes_read = 0;
      // printf("reading %lx..%lx\n", (uintptr_t)read_head, (uintptr_t)(read_head + read_attempt_size));
      switch (file_type) {
      case kFileUncompressed:
        {
          bytes_read = fread_unlocked(read_head, 1, read_attempt_size, ff);
          if (unlikely(ferror_unlocked(ff))) {
            goto TextRstreamThread_READ_FAIL;
          }
          break;
        }
      case kFileGzip:
        {
          // todo
          break;
        }
      case kFileBgzf:
        {
          // TODO
          exit(1);
        }
      case kFileZstd:
        {
          if ((!trsp->rds.zst.ib.size) && feof_unlocked(ff)) {
            break;
          }
          while (1) {
          }
          break;
        }
      }
      char* cur_read_end = &(read_head[S_CAST(uint32_t, bytes_read)]);
      if (bytes_read < S_CAST(int32_t, S_CAST(uint32_t, read_attempt_size))) {
        char* final_read_head = cur_read_end;
        if (cur_block_start != final_read_head) {
          if (final_read_head[-1] != '\n') {
            // Append '\n' so consumer can always use rawmemchr(., '\n') to
            // find the end of the current line.
            *final_read_head++ = '\n';
          }
        }
        // Still want to consistently enforce max line/token length.
        if (unlikely(IsPathologicallyLongLineOrTokenX(cur_block_start, read_head, final_read_head, enforced_max_line_blen))) {
          goto TextRstreamThread_LONG_LINE;
        }
        read_head = final_read_head;
        goto TextRstreamThread_EOF;
      }
      char* last_byte_ptr;
      if (!is_token_stream) {
        last_byte_ptr = Memrchr(read_head, '\n', read_attempt_size);
      } else {
        last_byte_ptr = LastSpaceOrEoln(read_head, read_attempt_size);
      }
      if (last_byte_ptr) {
        char* next_available_end = &(last_byte_ptr[1]);
        if (unlikely(IsPathologicallyLongLineOrTokenX(cur_block_start, read_head, next_available_end, enforced_max_line_blen))) {
          goto TextRstreamThread_LONG_LINE;
        }
#ifdef _WIN32
        EnterCriticalSection(critical_sectionp);
#else
        pthread_mutex_lock(sync_mutexp);
#endif
        interrupt = syncp->interrupt;
        if (interrupt != kTrsInterruptNone) {
          goto TextRstreamThread_INTERRUPT;
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
    TextRstreamThread_OPEN_FAIL:
      min_interrupt = kTrsInterruptShutdown;
      reterr = kPglRetOpenFail;
      break;
    TextRstreamThread_READ_FAIL:
      min_interrupt = kTrsInterruptShutdown;
      reterr = kPglRetReadFail;
      break;
    TextRstreamThread_NOMEM:
      min_interrupt = kTrsInterruptShutdown;
      reterr = kPglRetNomem;
      break;
    TextRstreamThread_LONG_LINE:
      min_interrupt = kTrsInterruptShutdown;
      reterr = kPglRetLongLine;
      break;
    TextRstreamThread_EOF:
      min_interrupt = kTrsInterruptRetarget;
      reterr = kPglRetEof;
      break;
    TextRstreamThread_OPEN_OR_READ_FAIL:
      min_interrupt = kTrsInterruptShutdown;
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
      goto TextRstreamThread_INTERRUPT;
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
  TextRstreamThread_INTERRUPT:
    // must be in critical section here, or be holding the mutex.
    if (interrupt == kTrsInterruptRetarget) {
      new_fname = syncp->new_fname;
      syncp->interrupt = kTrsInterruptNone;
      syncp->reterr = kPglRetSuccess;
    }
#ifdef _WIN32
    LeaveCriticalSection(critical_sectionp);
#else
    pthread_mutex_unlock(sync_mutexp);
#endif
    if (interrupt == kTrsInterruptShutdown) {
      // possible todo: close the file here
      THREAD_RETURN;
    }
    assert(interrupt == kTrsInterruptRetarget);
    if (!new_fname) {
    } else {
    }
    cur_block_start = buf;
    read_head = buf;
    read_stop = buf_end;
  }
  */
}

const char kShortErrRfileInvalid[] = "TextRstreamOpenEx can't be called with a closed or error-state textRFILE";

PglErr TextRstreamOpenEx(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, uint32_t decompress_thread_ct, textRFILE* trfp, char* dst, TextRstream* trsp) {
  PglErr reterr = kPglRetSuccess;
  {
    if (trfp) {
      // Move-construct (unless there was an error, or file is not opened)
      if (unlikely((!TextRfileIsOpen(trfp)) || TextRfileErrcode(trfp))) {
        reterr = kPglRetImproperFunctionCall;
        trsp->base.errmsg = kShortErrRfileInvalid;
        goto TextRstreamOpenEx_ret_1;
      }
      if (unlikely(TextRstreamIsOpen(trsp))) {
        reterr = kPglRetImproperFunctionCall;
        trsp->base.errmsg = kShortErrRfileAlreadyOpen;
        goto TextRstreamOpenEx_ret_1;
      }
      trsp->base = trfp->base;
      reterr = trfp->base.reterr;
      const FileCompressionType file_type = trfp->base.file_type;
      if (file_type != kFileUncompressed) {
        if (file_type == kFileGzip) {
          trsp->rds.gz = trfp->rds.gz;
        } else if (file_type == kFileZstd) {
          trsp->rds.zst = trfp->rds.zst;
        } else {
          // bgzf.  todo: BgzfRawMtDecompressStream move-construction from
          // BgzfRawDecompressStream.
          exit(1);
        }
      }
      EraseTextRfileBase(&trfp->base);
    } else {
      reterr = TextRfileOpenInternal(fname, enforced_max_line_blen, dst_capacity, decompress_thread_ct, dst, nullptr, trsp);
    }
    if (reterr) {
      if (reterr == kPglRetEof) {
        trsp->base.reterr = kPglRetEof;
        return kPglRetSuccess;
      }
      goto TextRstreamOpenEx_ret_1;
    }
    assert(!trsp->syncp);
    TextRstreamSync* syncp;
    if (unlikely(cachealigned_malloc(RoundUpPow2(sizeof(TextRstreamSync), kCacheline), &syncp))) {
      goto TextRstreamOpenEx_ret_NOMEM;
    }
    trsp->syncp = syncp;
    dst = trsp->base.dst;
    syncp->consume_tail = dst;
    syncp->cur_circular_end = nullptr;
    syncp->available_end = dst;
    syncp->reterr = kPglRetSuccess;
    syncp->interrupt = kTrsInterruptNone;
    syncp->new_fname = nullptr;
#ifdef _WIN32
    // apparently this can raise a low-memory exception in older Windows
    // versions, but that's not really our problem.
    InitializeCriticalSection(&syncp->critical_section);

    syncp->reader_progress_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (unlikely(!trsp->reader_progress_event)) {
      DeleteCriticalSection(&syncp->critical_section);
      goto TextRstreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    syncp->consumer_progress_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (unlikely(!syncp->consumer_progress_event)) {
      DeleteCriticalSection(&syncp->critical_section);
      CloseHandle(syncp->reader_progress_event);
      goto TextRstreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    trsp->read_thread = R_CAST(HANDLE, _beginthreadex(nullptr, kDefaultThreadStack, TextRstreamThread, trsp, 0, nullptr));
    if (unlikely(!trsp->read_thread)) {
      DeleteCriticalSection(&syncp->critical_section);
      CloseHandle(syncp->consumer_progress_event);
      CloseHandle(syncp->reader_progress_event);
      goto TextRstreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
#else
    if (unlikely(pthread_mutex_init(&syncp->sync_mutex, nullptr))) {
      goto TextRstreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    trsp->sync_init_state = 1;
    if (unlikely(pthread_cond_init(&syncp->reader_progress_condvar, nullptr))) {
      goto TextRstreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    trsp->sync_init_state = 2;
    syncp->consumer_progress_state = 0;
    if (unlikely(pthread_cond_init(&syncp->consumer_progress_condvar, nullptr))) {
      goto TextRstreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    trsp->sync_init_state = 3;
    if (unlikely(pthread_create(&trsp->read_thread, &g_smallstack_thread_attrx, TextRstreamThread, trsp))) {
      goto TextRstreamOpenEx_ret_THREAD_CREATE_FAIL;
    }
    trsp->sync_init_state = 4;
#endif
    // *consume_iterp = &(buf[-1]);
  }
  while (0) {
  TextRstreamOpenEx_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TextRstreamOpenEx_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 TextRstreamOpenEx_ret_1:
  trsp->base.reterr = reterr;
  return reterr;
}

uint32_t TextDecompressThreadCt(const TextRstream* trsp) {
  FileCompressionType file_type = trsp->base.file_type;
  if (file_type == kFileUncompressed) {
    return 0;
  }
  if (file_type != kFileBgzf) {
    return 1;
  }
  // TODO
  return 1;
}

BoolErr CleanupTextRstream(TextRstream* trsp, PglErr* reterrp) {
#ifdef _WIN32
  if (trsp->read_thread) {
    TextRstreamSync* syncp = trsp->syncp;
    CRITICAL_SECTION* critical_sectionp = &syncp->critical_section;
    EnterCriticalSection(critical_sectionp);
    syncp->interrupt = kTrsInterruptShutdown;
    SetEvent(trsp->consumer_progress_event);
    LeaveCriticalSection(critical_sectionp);
    WaitForSingleObject(trsp->read_thread, INFINITE);
    DeleteCriticalSection(critical_sectionp);
    trsp->read_thread = nullptr;  // make it safe to call this multiple times
    CloseHandle(syncp->consumer_progress_event);
    CloseHandle(syncp->reader_progress_event);
  }
#else
  const uint32_t sync_init_state = trsp->sync_init_state;
  if (sync_init_state) {
    TextRstreamSync* syncp = trsp->syncp;
    pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
    pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
    if (sync_init_state == 4) {
      pthread_mutex_lock(sync_mutexp);
      syncp->interrupt = kTrsInterruptShutdown;
      syncp->consumer_progress_state = 1;
      pthread_cond_signal(consumer_progress_condvarp);
      pthread_mutex_unlock(sync_mutexp);
      pthread_join(trsp->read_thread, nullptr);
    }
    pthread_mutex_destroy(sync_mutexp);
    if (sync_init_state > 1) {
      pthread_cond_destroy(&syncp->reader_progress_condvar);
      if (sync_init_state > 2) {
        pthread_cond_destroy(consumer_progress_condvarp);
      }
    }
    trsp->sync_init_state = 0;  // make it safe to call this multiple times
  }
#endif
  if (trsp->syncp) {
    aligned_free(trsp->syncp);
    trsp->syncp = nullptr;
  }
  trsp->base.consume_iter = nullptr;
  trsp->base.consume_stop = nullptr;
  trsp->base.reterr = kPglRetEof;
  trsp->base.errmsg = nullptr;
  if (trsp->base.dst && (!trsp->base.dst_owned_by_caller)) {
    free(trsp->base.dst);
    trsp->base.dst = nullptr;
  }
  if (trsp->base.ff) {
    if (trsp->base.file_type != kFileUncompressed) {
      if (trsp->base.file_type == kFileZstd) {
        if (trsp->rds.zst.ib.src) {
          free_const(trsp->rds.zst.ib.src);
          trsp->rds.zst.ib.src = nullptr;
        }
        if (trsp->rds.zst.ds) {
          ZSTD_freeDStream(trsp->rds.zst.ds);
          trsp->rds.zst.ds = nullptr;
        }
      } else if (trsp->base.file_type == kFileBgzf) {
        // TODO
        exit(1);
      } else {
        // plain gzip
        if (trsp->rds.gz.in) {
          free(trsp->rds.gz.in);
          trsp->rds.gz.in = nullptr;
        }
        if (trsp->rds.gz.ds_initialized) {
          inflateEnd(&trsp->rds.gz.ds);
        }
      }
      trsp->base.file_type = kFileUncompressed;
    }
    if (unlikely(fclose_null(&trsp->base.ff))) {
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

#ifdef __cplusplus
}
#endif
