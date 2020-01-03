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
#include "plink2_compress_stream.h"

#ifdef __cplusplus
namespace plink2 {
#endif

uint32_t g_zst_level = 0;

PglErr InitCstreamNoop(const char* out_fname, uint32_t do_append, char* overflow_buf, CompressStreamState* css_ptr) {
  // css_ptr->z_outfile = nullptr;
  css_ptr->cctx = nullptr;
  // can't use fopen_checked since we need to be able to append
  css_ptr->outfile = fopen(out_fname, do_append? FOPEN_AB : FOPEN_WB);
  if (unlikely(!css_ptr->outfile)) {
    logputs("\n");
    logerrprintfww(kErrprintfFopen, out_fname, strerror(errno));
    return kPglRetOpenFail;
  }
  css_ptr->overflow_buf = overflow_buf;
  return kPglRetSuccess;
}

PglErr InitCstreamZstd(const char* out_fname, uint32_t do_append, __maybe_unused uint32_t thread_ct, uintptr_t overflow_buf_size, char* overflow_buf, unsigned char* compress_wkspace, CompressStreamState* css_ptr) {
  css_ptr->outfile = nullptr;
  css_ptr->cctx = ZSTD_createCCtx();
  if (unlikely(!css_ptr->cctx)) {
    return kPglRetNomem;
  }
  __maybe_unused size_t retval = ZSTD_CCtx_setParameter(css_ptr->cctx, ZSTD_c_compressionLevel, g_zst_level);
  assert(!ZSTD_isError(retval));
#ifdef ZSTD_MULTITHREAD
  // ignore failure; if zstd is dynamically linked and was built without MT
  // support, so be it
  ZSTD_CCtx_setParameter(css_ptr->cctx, ZSTD_c_nbWorkers, thread_ct);
#endif
  css_ptr->outfile = fopen(out_fname, do_append? FOPEN_AB : FOPEN_WB);
  if (unlikely(!css_ptr->outfile)) {
    logputs("\n");
    logerrprintfww(kErrprintfFopen, out_fname, strerror(errno));
    ZSTD_freeCCtx(css_ptr->cctx);  // might return an error later?
    return kPglRetOpenFail;
  }
  css_ptr->output.dst = compress_wkspace;
  css_ptr->output.size = CstreamWkspaceReq(overflow_buf_size);
  css_ptr->output.pos = 0;
  css_ptr->overflow_buf = overflow_buf;
  return kPglRetSuccess;
}

// possible todo: replace output_zst with an enum which permits gzipping
PglErr InitCstream(const char* out_fname, uint32_t do_append, uint32_t output_zst, uint32_t thread_ct, uintptr_t overflow_buf_size, char* overflow_buf, unsigned char* compress_wkspace, CompressStreamState* css_ptr) {
  if (!output_zst) {
    return InitCstreamNoop(out_fname, do_append, overflow_buf, css_ptr);
  }
  return InitCstreamZstd(out_fname, do_append, thread_ct, overflow_buf_size, overflow_buf, compress_wkspace, css_ptr);
}

PglErr InitCstreamAlloc(const char* out_fname, uint32_t do_append, uint32_t output_zst, uint32_t thread_ct, uintptr_t overflow_buf_size, CompressStreamState* css_ptr, char** cswritepp) {
  char* overflow_buf;
  if (unlikely(bigstack_alloc_c(overflow_buf_size, &overflow_buf))) {
    return kPglRetNomem;
  }
  unsigned char* compress_wkspace = nullptr;
  if (output_zst) {
    if (unlikely(bigstack_alloc_uc(CstreamWkspaceReq(overflow_buf_size), &compress_wkspace))) {
      return kPglRetNomem;
    }
  }
  PglErr reterr = InitCstream(out_fname, do_append, output_zst, thread_ct, overflow_buf_size, overflow_buf, compress_wkspace, css_ptr);
  *cswritepp = overflow_buf;
  return reterr;
}

BoolErr ForceUncompressedCswrite(CompressStreamState* css_ptr, char** writep_ptr) {
  char* writep = *writep_ptr;
  if (css_ptr->overflow_buf != writep) {
    if (unlikely(!fwrite_unlocked(css_ptr->overflow_buf, writep - css_ptr->overflow_buf, 1, css_ptr->outfile))) {
      return 1;
    }
    *writep_ptr = css_ptr->overflow_buf;
  }
  return 0;
}

BoolErr ForceCompressedCswrite(CompressStreamState* css_ptr, char** writep_ptr) {
  char* overflow_buf = css_ptr->overflow_buf;
  char* writep = *writep_ptr;
  if (overflow_buf != writep) {
    const uintptr_t in_size = writep - overflow_buf;
    ZSTD_inBuffer input = {overflow_buf, in_size, 0};
    while (1) {
      // todo: conditionally support seekable files
      __maybe_unused size_t retval = ZSTD_compressStream2(css_ptr->cctx, &css_ptr->output, &input, ZSTD_e_continue);
      assert(!ZSTD_isError(retval));
      if (css_ptr->output.pos >= kCompressStreamBlock) {
        if (unlikely(!fwrite_unlocked(css_ptr->output.dst, css_ptr->output.pos, 1, css_ptr->outfile))) {
          return 1;
        }
        css_ptr->output.pos = 0;
      }
      const uintptr_t bytes_left = input.size - input.pos;
      if (bytes_left < kCompressStreamBlock) {
        memmove(overflow_buf, &(overflow_buf[in_size - bytes_left]), bytes_left);
        *writep_ptr = &(overflow_buf[bytes_left]);
        break;
      }
    }
  }
  return 0;
}

BoolErr CsputsStd(const char* readp, uint32_t byte_ct, CompressStreamState* css_ptr, char** writep_ptr) {
  char* writep = *writep_ptr;
  char* overflow_buf = css_ptr->overflow_buf;
  uint32_t cur_write_space = 2 * kCompressStreamBlock - S_CAST(uintptr_t, writep - overflow_buf);
  if (IsUncompressedCstream(css_ptr)) {
    while (byte_ct > cur_write_space) {
      memcpy(writep, readp, cur_write_space);
      if (unlikely(!fwrite_unlocked(overflow_buf, 2 * kCompressStreamBlock, 1, css_ptr->outfile))) {
        return 1;
      }
      writep = overflow_buf;
      readp = &(readp[cur_write_space]);
      byte_ct -= cur_write_space;
      cur_write_space = 2 * kCompressStreamBlock;
    }
  } else {
    while (byte_ct > cur_write_space) {
      memcpy(writep, readp, cur_write_space);
      ZSTD_inBuffer input = {overflow_buf, 2 * kCompressStreamBlock, 0};
      while (1) {
        // todo: conditionally support seekable files
        __maybe_unused size_t retval = ZSTD_compressStream2(css_ptr->cctx, &css_ptr->output, &input, ZSTD_e_continue);
        assert(!ZSTD_isError(retval));
        if (css_ptr->output.pos >= kCompressStreamBlock) {
          if (unlikely(!fwrite_unlocked(css_ptr->output.dst, css_ptr->output.pos, 1, css_ptr->outfile))) {
            return 1;
          }
          css_ptr->output.pos = 0;
        }
        const uintptr_t bytes_left = input.size - input.pos;
        if (bytes_left < kCompressStreamBlock) {
          memmove(overflow_buf, &(overflow_buf[2 * kCompressStreamBlock - bytes_left]), bytes_left);
          writep = &(overflow_buf[retval]);
          break;
        }
      }
      readp = &(readp[cur_write_space]);
      byte_ct -= cur_write_space;
      cur_write_space = 2 * kCompressStreamBlock - S_CAST(uintptr_t, writep - overflow_buf);
    }
  }
  memcpy(writep, readp, byte_ct);
  *writep_ptr = &(writep[byte_ct]);
  return Cswrite(css_ptr, writep_ptr);
}

BoolErr UncompressedCswriteCloseNull(CompressStreamState* css_ptr, char* writep) {
  ForceUncompressedCswrite(css_ptr, &writep);
  css_ptr->overflow_buf = nullptr;
  int32_t ii = ferror_unlocked(css_ptr->outfile);
  int32_t jj = fclose(css_ptr->outfile);
  return ii || jj;
}

BoolErr CompressedCswriteCloseNull(CompressStreamState* css_ptr, char* writep) {
  char* overflow_buf = css_ptr->overflow_buf;
  const uintptr_t in_size = writep - overflow_buf;
  ZSTD_inBuffer input = {overflow_buf, in_size, 0};
  BoolErr reterr = 0;
  while (1) {
    __maybe_unused size_t retval = ZSTD_compressStream2(css_ptr->cctx, &css_ptr->output, &input, ZSTD_e_end);
    assert(!ZSTD_isError(retval));
    if (css_ptr->output.pos) {
      if (unlikely(!fwrite_unlocked(css_ptr->output.dst, css_ptr->output.pos, 1, css_ptr->outfile))) {
        reterr = 1;
      }
      css_ptr->output.pos = 0;
    }
    if (!retval) {
      break;
    }
  }
  ZSTD_freeCCtx(css_ptr->cctx);  // might return an error later?
  css_ptr->overflow_buf = nullptr;
  int32_t ii = ferror_unlocked(css_ptr->outfile);
  int32_t jj = fclose(css_ptr->outfile);
  return reterr || ii || jj;
}

BoolErr CswriteCloseNull(CompressStreamState* css_ptr, char* writep) {
  if (IsUncompressedCstream(css_ptr)) {
    return UncompressedCswriteCloseNull(css_ptr, writep);
  }
  return CompressedCswriteCloseNull(css_ptr, writep);
}

#ifdef __cplusplus
}
#endif
