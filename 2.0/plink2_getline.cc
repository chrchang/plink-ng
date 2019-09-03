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

void PreinitTextRfile(textRFILE* trfp) {
  trfp->consume_iter = nullptr;
  trfp->consume_stop = nullptr;
  trfp->dst = nullptr;
  trfp->ff = nullptr;
  trfp->errmsg = nullptr;
  trfp->reterr = kPglRetEof;
  trfp->in = nullptr;
}

const char kShortErrRfileAlreadyOpen[] = "TextRfileOpenEx can't be called on an already-open file";
const char kShortErrRfileEnforcedMaxBlenTooSmall[] = "TextRfileOpenEx: enforced_max_line_blen too small (must be at least max(1 MiB, dst_capacity - 1 MiB))";
const char kShortErrRfileDstCapacityTooSmall[] = "TextRfileOpenEx: dst_capacity too small (2 MiB minimum)";

PglErr TextRfileOpenEx(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, char* dst, textRFILE* trfp) {
  PglErr reterr = kPglRetSuccess;
  {
    // 1. Open file, get type.
    if (unlikely(trfp->ff)) {
      reterr = kPglRetImproperFunctionCall;
      trfp->errmsg = kShortErrRfileAlreadyOpen;
      goto TextRfileOpenEx_ret_1;
    }
    if (unlikely(enforced_max_line_blen < kDecompressChunkSize)) {
      reterr = kPglRetImproperFunctionCall;
      trfp->errmsg = kShortErrRfileEnforcedMaxBlenTooSmall;
      goto TextRfileOpenEx_ret_1;
    }
    if (dst) {
      if (unlikely(dst_capacity < 2 * kDecompressChunkSize)) {
        reterr = kPglRetImproperFunctionCall;
        trfp->errmsg = kShortErrRfileDstCapacityTooSmall;
        goto TextRfileOpenEx_ret_1;
      }
      if (unlikely(enforced_max_line_blen + kDecompressChunkSize < dst_capacity)) {
        reterr = kPglRetImproperFunctionCall;
        trfp->errmsg = kShortErrRfileEnforcedMaxBlenTooSmall;
        goto TextRfileOpenEx_ret_1;
      }
    }
    trfp->ff = fopen(fname, FOPEN_RB);
    if (unlikely(!trfp->ff)) {
      goto TextRfileOpenEx_ret_OPEN_FAIL;
    }
    if (dst) {
      trfp->dst_owned_by_caller = 1;
      trfp->dst_capacity = dst_capacity;
    } else {
      dst = S_CAST(char*, malloc(2 * kDecompressChunkSize));
      if (unlikely(dst == nullptr)) {
        goto TextRfileOpenEx_ret_NOMEM;
      }
      trfp->dst_owned_by_caller = 0;
      trfp->dst_capacity = 2 * kDecompressChunkSize;
    }
    trfp->dst = dst;
    uint32_t nbytes = fread_unlocked(dst, 1, 16, trfp->ff);
    trfp->file_type = kFileUncompressed;
    trfp->dst_len = nbytes;
    trfp->consume_iter = dst;
    trfp->consume_stop = dst;
    if (nbytes >= 4) {
      const uint32_t magic4 = *R_CAST(uint32_t*, dst);
      if (IsZstdFrame(magic4)) {
        trfp->dst_len = 0;
        trfp->file_type = kFileZstd;
        trfp->raw.zst.ds = ZSTD_createDStream();
        if (unlikely(!trfp->raw.zst.ds)) {
          trfp->raw.zst.ib.src = nullptr;
          goto TextRfileOpenEx_ret_NOMEM;
        }
        trfp->in = S_CAST(unsigned char*, malloc(kDecompressChunkSize));
        if (unlikely(!trfp->in)) {
          goto TextRfileOpenEx_ret_NOMEM;
        }
        memcpy(trfp->in, dst, nbytes);
        trfp->raw.zst.ib.src = trfp->in;
        trfp->raw.zst.ib.size = nbytes;
        trfp->raw.zst.ib.pos = 0;
      } else if ((magic4 << 8) == 0x088b1f00) {
        // gzip ID1/ID2 bytes, deflate compression method
        trfp->in = S_CAST(unsigned char*, malloc(kDecompressChunkSize));
        if (unlikely(!trfp->in)) {
          goto TextRfileOpenEx_ret_NOMEM;
        }
        memcpy(trfp->in, dst, nbytes);
        trfp->dst_len = 0;
        if ((nbytes == 16) && IsBgzfHeader(trfp->in)) {
          trfp->file_type = kFileBgzf;
          trfp->raw.bgzf.ldc = libdeflate_alloc_decompressor();
          if (!trfp->raw.bgzf.ldc) {
            goto TextRfileOpenEx_ret_NOMEM;
          }
          trfp->raw.bgzf.in_size = nbytes;
          trfp->raw.bgzf.in_pos = 0;
        } else {
          trfp->file_type = kFileGzip;
          z_stream* dsp = &trfp->raw.gz.ds;
          dsp->next_in = trfp->in;
          dsp->avail_in = nbytes;
          dsp->zalloc = nullptr;
          dsp->zfree = nullptr;
          dsp->opaque = nullptr;
          if (unlikely(inflateInit2(dsp, MAX_WBITS | 16) != Z_OK)) {
            // todo: verify no other errors possible
            trfp->raw.gz.ds_initialized = 0;
            goto TextRfileOpenEx_ret_NOMEM;
          }
          trfp->raw.gz.ds_initialized = 1;
        }
      }
    } else if (!nbytes) {
      if (unlikely(!feof_unlocked(trfp->ff))) {
        goto TextRfileOpenEx_ret_READ_FAIL;
      }
      // May as well accept this.
      // Don't jump to ret_1 since we're setting trfp->reterr to a different
      // value than we're returning.
      trfp->reterr = kPglRetEof;
      return kPglRetSuccess;
    }
  }
  while (0) {
  TextRfileOpenEx_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TextRfileOpenEx_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    trfp->errmsg = strerror(errno);
    break;
  TextRfileOpenEx_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    trfp->errmsg = strerror(errno);
    break;
  }
 TextRfileOpenEx_ret_1:
  trfp->reterr = reterr;
  return reterr;
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
  if (trfp->reterr) {
    return trfp->reterr;
  }
  PglErr reterr = kPglRetSuccess;
  {
    char* orig_line_start = trfp->consume_stop;
    assert(trfp->consume_iter == orig_line_start);
    char* dst = trfp->dst;
    char* dst_load_start;
    while (1) {
      const uint32_t dst_offset = orig_line_start - dst;
      const uint32_t dst_rem = trfp->dst_len - dst_offset;
      // (dst_rem guaranteed to be < trfp->enforced_max_line_blen here, since
      // otherwise we error out earlier.)
      // Two cases:
      // 1. Move (possibly empty) unfinished line to the beginning of the
      //    buffer.
      // 2. Resize the buffer/report out-of-memory.
      if (dst_rem < trfp->dst_capacity - kDecompressChunkSize) {
        memmove(dst, orig_line_start, dst_rem);
      } else {
        if (unlikely(trfp->dst_owned_by_caller)) {
          goto TextRfileAdvance_ret_NOMEM;
        }
        uint32_t next_dst_capacity = trfp->enforced_max_line_blen + kDecompressChunkSize;
        if ((next_dst_capacity / 2) > trfp->dst_capacity) {
          next_dst_capacity = trfp->dst_capacity * 2;
        }
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
        trfp->dst = dst_next;
        dst = dst_next;
      }
      dst_load_start = &(dst[dst_rem]);
      char* dst_iter = dst_load_start;
      char* dst_end = &(dst[trfp->dst_capacity]);
      trfp->consume_iter = dst;
      FILE* ff = trfp->ff;
      switch (trfp->file_type) {
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
              if (ferror_unlocked(ff)) {
                goto TextRfileAdvance_ret_READ_FAIL;
              }
              trfp->dst_len = nbytes + dst_rem;
              break;
            }
            rlen -= kMaxBytesPerIO;
            dst_iter = &(dst_iter[kMaxBytesPerIO]);
          }
          const uint32_t nbytes = fread_unlocked(dst_iter, 1, rlen, ff);
          if (ferror_unlocked(ff)) {
            goto TextRfileAdvance_ret_READ_FAIL;
          }
          dst_iter = &(dst_iter[nbytes]);
          break;
        }
      case kFileGzip:
        {
          z_stream* dsp = &trfp->raw.gz.ds;
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
                  trfp->errmsg = dsp->msg;
                } else {
                  trfp->errmsg = zError(zerr);
                }
                goto TextRfileAdvance_ret_DECOMPRESS_FAIL;
              }
              dst_iter = R_CAST(char*, dsp->next_out);
              if (dsp->avail_in) {
                assert(dst_iter == dst_end);
                break;
              }
            }
            const uint32_t nbytes = fread_unlocked(trfp->in, 1, kDecompressChunkSize, ff);
            dsp->next_in = trfp->in;
            dsp->avail_in = nbytes;
            if (!nbytes) {
              if (unlikely(!feof_unlocked(ff))) {
                goto TextRfileAdvance_ret_READ_FAIL;
              }
              if (unlikely(zerr == Z_OK)) {
                trfp->errmsg = kShortErrRfileTruncatedGz;
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
          if ((!trfp->raw.bgzf.in_size) && feof_unlocked(ff)) {
            break;
          }
          struct libdeflate_decompressor* ldc = trfp->raw.bgzf.ldc;
          unsigned char* in = trfp->in;
          unsigned char* in_iter = &(in[trfp->raw.bgzf.in_pos]);
          unsigned char* in_end = &(in[trfp->raw.bgzf.in_size]);
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
            trfp->raw.bgzf.in_size = in_end - in;
            if (!nbytes) {
              if (unlikely(n_inbytes)) {
                goto TextRfileAdvance_ret_INVALID_BGZF;
              }
              break;
            }
          }
          trfp->raw.bgzf.in_pos = in_iter - in;
          dst_end = dst_iter;
          break;
        }
      case kFileZstd:
        {
          if ((!trfp->raw.zst.ib.size) && feof_unlocked(ff)) {
            break;
          }
          // Sequentially dependent blocks limited to ~128 KiB.
          while (1) {
            ZSTD_outBuffer zob = {R_CAST(unsigned char*, dst_iter), S_CAST(size_t, dst_end - dst_iter), 0};
            // ib.size == 0 ok, no need to special-case rewind.
            const uintptr_t read_size_hint = ZSTD_decompressStream(trfp->raw.zst.ds, &zob, &trfp->raw.zst.ib);
            if (unlikely(ZSTD_isError(read_size_hint))) {
              trfp->errmsg = ZSTD_getErrorName(read_size_hint);
              goto TextRfileAdvance_ret_DECOMPRESS_FAIL;
            }
            dst_iter = &(dst_iter[zob.pos]);
            if (dst_iter == dst_end) {
              break;
            }
            // Decoder has flushed everything it could.  Either we're at EOF,
            // or we must load more.
            unsigned char* in = trfp->in;
            const uint32_t n_inbytes = trfp->raw.zst.ib.size - trfp->raw.zst.ib.pos;
            memmove(in, &(in[trfp->raw.zst.ib.pos]), n_inbytes);
            unsigned char* load_start = &(in[n_inbytes]);
            const uint32_t nbytes = fread_unlocked(load_start, 1, kDecompressChunkSize - n_inbytes, ff);
            if (unlikely(ferror_unlocked(ff))) {
              goto TextRfileAdvance_ret_READ_FAIL;
            }
            trfp->raw.zst.ib.pos = 0;
            trfp->raw.zst.ib.size = nbytes + n_inbytes;
            if (!nbytes) {
              if (unlikely(n_inbytes)) {
                trfp->errmsg = kShortErrZstdPrefixUnknown;
                goto TextRfileAdvance_ret_DECOMPRESS_FAIL;
              }
              break;
            }
          }
          break;
        }
      }
      trfp->dst_len = dst_iter - dst;
      if (!trfp->dst_len) {
        goto TextRfileAdvance_ret_EOF;
      }
      if (dst_iter != dst_end) {
        // If last character of file isn't a newline, append one to simplify
        // downstream code.
        if (dst_iter[-1] != '\n') {
          *dst_iter++ = '\n';
          trfp->dst_len += 1;
        }
        trfp->consume_stop = dst_iter;
        break;
      }
      char* last_byte_ptr = Memrchr(dst_load_start, '\n', dst_iter - dst_load_start);
      if (last_byte_ptr) {
        trfp->consume_stop = &(last_byte_ptr[1]);
        break;
      }
      // Buffer is full, and no '\n' is present.  Restart the loop and try to
      // extend the buffer, if we aren't already at/past the line-length limit.
      if (trfp->dst_len >= trfp->enforced_max_line_blen) {
        goto TextRfileAdvance_ret_LONG_LINE;
      }
    }
    if (unlikely(IsPathologicallyLongLineOrTokenX(dst, dst_load_start, trfp->consume_stop, trfp->enforced_max_line_blen))) {
      goto TextRfileAdvance_ret_LONG_LINE;
    }
  }
  while (0) {
  TextRfileAdvance_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TextRfileAdvance_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    trfp->errmsg = strerror(errno);
    break;
  TextRfileAdvance_ret_LONG_LINE:
    trfp->errmsg = kShortErrLongLine;
    reterr = kPglRetMalformedInput;
    break;
  TextRfileAdvance_ret_INVALID_BGZF:
    trfp->errmsg = kShortErrInvalidBgzf;
  TextRfileAdvance_ret_DECOMPRESS_FAIL:
    reterr = kPglRetDecompressFail;
    break;
  TextRfileAdvance_ret_EOF:
    reterr = kPglRetEof;
    break;
  }
  trfp->reterr = reterr;
  return reterr;
}

void TextRfileRewind(textRFILE* trfp) {
  if ((!trfp->ff) || ((trfp->reterr) && (trfp->reterr != kPglRetEof))) {
    return;
  }
  rewind(trfp->ff);
  trfp->reterr = kPglRetSuccess;
  trfp->dst_len = 0;
  trfp->consume_iter = trfp->dst;
  trfp->consume_stop = trfp->dst;
  if (trfp->file_type != kFileUncompressed) {
    if (trfp->file_type == kFileGzip) {
      trfp->raw.gz.ds.avail_in = 0;
#ifdef NDEBUG
      inflateReset(&trfp->raw.gz.ds);
#else
      const int errcode = inflateReset(&trfp->raw.gz.ds);
      assert(errcode == Z_OK);
#endif
    } else if (trfp->file_type == kFileBgzf) {
      trfp->raw.bgzf.in_size = 0;
      trfp->raw.bgzf.in_pos = 0;
    } else {
      // kFileZstd
      trfp->raw.zst.ib.size = 0;
      trfp->raw.zst.ib.pos = 0;
      ZSTD_DCtx_reset(trfp->raw.zst.ds, ZSTD_reset_session_only);
    }
  }
}

BoolErr CleanupTextRfile(textRFILE* trfp, PglErr* reterrp) {
  trfp->reterr = kPglRetEof;
  trfp->errmsg = nullptr;
  if (trfp->in) {
    free(trfp->in);
    trfp->in = nullptr;
  }
  if (trfp->dst && (!trfp->dst_owned_by_caller)) {
    free(trfp->dst);
    trfp->dst = nullptr;
  }
  if (trfp->ff) {
    if (trfp->file_type != kFileUncompressed) {
      if (trfp->file_type == kFileZstd) {
        if (trfp->raw.zst.ds) {
          ZSTD_freeDStream(trfp->raw.zst.ds);
          trfp->raw.zst.ds = nullptr;
        }
      } else if (trfp->file_type == kFileBgzf) {
        if (trfp->raw.bgzf.ldc) {
          libdeflate_free_decompressor(trfp->raw.bgzf.ldc);
          trfp->raw.bgzf.ldc = nullptr;
        }
      } else {
        // plain gzip
        if (trfp->raw.gz.ds_initialized) {
          inflateEnd(&trfp->raw.gz.ds);
        }
      }
      trfp->file_type = kFileUncompressed;
    }
    if (unlikely(fclose_null(&trfp->ff))) {
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


/*
void PreinitText(TextStream* tsp) {
  tsp->consume_iter = nullptr;
  tsp->consume_stop = nullptr;
  tsp->dst = nullptr;
  tsp->ff = nullptr;
  tsp->errmsg = nullptr;
  tsp->reterr = kPglRetEof;
  tsp->in = nullptr;
}

BoolErr CleanupText(TextStream* tsp, PglErr* reterrp) {
  tsp->reterr = kPglRetEof;
  tsp->errmsg = nullptr;
  if (tsp->in) {
    free(tsp->in);
    tsp->in = nullptr;
  }
  if (tsp->dst && (!tsp->dst_owned_by_caller)) {
    free(tsp->dst);
    tsp->dst = nullptr;
  }
  if (tsp->ff) {
    if (tsp->file_type != kFileUncompressed) {
      if (tsp->file_type == kFileZstd) {
        if (tsp->raw.zst.ds) {
          ZSTD_freeDStream(tsp->raw.zst.ds);
          tsp->raw.zst.ds = nullptr;
        }
      } else if (tsp->file_type == kFileBgzf) {
        // TODO
      } else {
        // plain gzip
        if (tsp->raw.gz.ds_initialized) {
          inflateEnd(&tsp->raw.gz.ds);
        }
      }
      tsp->file_type = kFileUncompressed;
    }
    if (unlikely(fclose_null(&tsp->ff))) {
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
*/

#ifdef __cplusplus
}
#endif
