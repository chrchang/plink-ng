// This library is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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


#include "plink2_compress_stream.h"

#ifdef __cplusplus
namespace plink2 {
#endif

uint32_t g_zst_level = 0;

pglerr_t uncompressed_cswrite_init(const char* out_fname, uint32_t do_append, char* overflow_buf, compress_stream_state_t* css_ptr) {
  // css_ptr->z_outfile = nullptr;
  css_ptr->cctx = nullptr;
  // can't use fopen_checked since we need to be able to append
  css_ptr->outfile = fopen(out_fname, do_append? FOPEN_AB : FOPEN_WB);
  if (!css_ptr->outfile) {
    logprint("\n");
    LOGERRPRINTFWW(g_errstr_fopen, out_fname);
    return kPglRetOpenFail;
  }
  css_ptr->overflow_buf = overflow_buf;
  return kPglRetSuccess;
}

pglerr_t zstd_cswrite_init(const char* out_fname, uint32_t do_append, __maybe_unused uint32_t thread_ct, uintptr_t overflow_buf_size, char* overflow_buf, unsigned char* compress_wkspace, compress_stream_state_t* css_ptr) {
  css_ptr->outfile = nullptr;
  css_ptr->cctx = ZSTD_createCCtx();
  if (!css_ptr->cctx) {
    return kPglRetNomem;
  }
  __maybe_unused size_t retval = ZSTD_CCtx_setParameter(css_ptr->cctx, ZSTD_p_compressionLevel, g_zst_level);
  assert(!ZSTD_isError(retval));
#ifdef ZSTD_MULTITHREAD
  retval = ZSTD_CCtx_setParameter(css_ptr->cctx, ZSTD_p_nbThreads, thread_ct);
  if (ZSTD_isError(retval)) {
    ZSTD_freeCCtx(css_ptr->cctx);
    return kPglRetNomem;
  }
#endif
  css_ptr->outfile = fopen(out_fname, do_append? FOPEN_AB : FOPEN_WB);
  if (!css_ptr->outfile) {
    logprint("\n");
    LOGERRPRINTFWW(g_errstr_fopen, out_fname);
    ZSTD_freeCCtx(css_ptr->cctx);  // might return an error later?
    return kPglRetOpenFail;
  }
  css_ptr->output.dst = compress_wkspace;
  css_ptr->output.size = css_wkspace_req(overflow_buf_size);
  css_ptr->output.pos = 0;
  css_ptr->overflow_buf = overflow_buf;
  return kPglRetSuccess;
}

// possible todo: replace output_zst with an enum which permits gzipping
pglerr_t cswrite_init(const char* out_fname, uint32_t do_append, uint32_t output_zst, uint32_t thread_ct, uintptr_t overflow_buf_size, char* overflow_buf, unsigned char* compress_wkspace, compress_stream_state_t* css_ptr) {
  if (!output_zst) {
    return uncompressed_cswrite_init(out_fname, do_append, overflow_buf, css_ptr);
  }
  return zstd_cswrite_init(out_fname, do_append, thread_ct, overflow_buf_size, overflow_buf, compress_wkspace, css_ptr);
}

pglerr_t cswrite_init2(const char* out_fname, uint32_t do_append, uint32_t output_zst, uint32_t thread_ct, uintptr_t overflow_buf_size, compress_stream_state_t* css_ptr, char** cswritepp) {
  char* overflow_buf;
  if (bigstack_alloc_c(overflow_buf_size, &overflow_buf)) {
    return kPglRetNomem;
  }
  unsigned char* compress_wkspace = nullptr;
  if (output_zst) {
    if (bigstack_alloc_uc(css_wkspace_req(overflow_buf_size), &compress_wkspace)) {
      return kPglRetNomem;
    }
  }
  pglerr_t reterr = cswrite_init(out_fname, do_append, output_zst, thread_ct, overflow_buf_size, overflow_buf, compress_wkspace, css_ptr);
  *cswritepp = overflow_buf;
  return reterr;
}

boolerr_t force_uncompressed_cswrite(compress_stream_state_t* css_ptr, char** writep_ptr) {
  char* writep = *writep_ptr;
  if (css_ptr->overflow_buf != writep) {
    if (!fwrite_unlocked(css_ptr->overflow_buf, writep - css_ptr->overflow_buf, 1, css_ptr->outfile)) {
      return 1;
    }
    *writep_ptr = css_ptr->overflow_buf;
  }
  return 0;
}

boolerr_t force_compressed_cswrite(compress_stream_state_t* css_ptr, char** writep_ptr) {
  char* overflow_buf = css_ptr->overflow_buf;
  char* writep = *writep_ptr;
  if (overflow_buf != writep) {
    const uintptr_t in_size = writep - overflow_buf;
    ZSTD_inBuffer input = {overflow_buf, in_size, 0};
    while (1) {
      // todo: conditionally support seekable files
      size_t retval = ZSTD_compress_generic(css_ptr->cctx, &css_ptr->output, &input, ZSTD_e_continue);
      if (ZSTD_isError(retval)) {
        // is this actually possible?  well, play it safe for now
        return 1;
      }
      if (css_ptr->output.pos >= kCompressStreamBlock) {
        if (!fwrite_unlocked(css_ptr->output.dst, css_ptr->output.pos, 1, css_ptr->outfile)) {
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

boolerr_t csputs_std(const char* readp, uint32_t byte_ct, compress_stream_state_t* css_ptr, char** writep_ptr) {
  char* writep = *writep_ptr;
  char* overflow_buf = css_ptr->overflow_buf;
  uint32_t cur_write_space = 2 * kCompressStreamBlock - S_CAST(uintptr_t, writep - overflow_buf);
  if (is_uncompressed_cswrite(css_ptr)) {
    while (byte_ct > cur_write_space) {
      memcpy(writep, readp, cur_write_space);
      if (!fwrite_unlocked(overflow_buf, 2 * kCompressStreamBlock, 1, css_ptr->outfile)) {
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
        size_t retval = ZSTD_compress_generic(css_ptr->cctx, &css_ptr->output, &input, ZSTD_e_continue);
        if (ZSTD_isError(retval)) {
          return 1;
        }
        if (css_ptr->output.pos >= kCompressStreamBlock) {
          if (!fwrite_unlocked(css_ptr->output.dst, css_ptr->output.pos, 1, css_ptr->outfile)) {
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
  return cswrite(css_ptr, writep_ptr);
}

boolerr_t uncompressed_cswrite_close_null(compress_stream_state_t* css_ptr, char* writep) {
  force_uncompressed_cswrite(css_ptr, &writep);
  css_ptr->overflow_buf = nullptr;
  int32_t ii = ferror_unlocked(css_ptr->outfile);
  int32_t jj = fclose(css_ptr->outfile);
  return ii || jj;
}

boolerr_t compressed_cswrite_close_null(compress_stream_state_t* css_ptr, char* writep) {
  char* overflow_buf = css_ptr->overflow_buf;
  const uintptr_t in_size = writep - overflow_buf;
  ZSTD_inBuffer input = {overflow_buf, in_size, 0};
  boolerr_t reterr = 0;
  while (1) {
    size_t retval = ZSTD_compress_generic(css_ptr->cctx, &css_ptr->output, &input, ZSTD_e_end);
    if (ZSTD_isError(retval)) {
      reterr = 1;
      break;
    }
    if (css_ptr->output.pos) {
      if (!fwrite_unlocked(css_ptr->output.dst, css_ptr->output.pos, 1, css_ptr->outfile)) {
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

boolerr_t cswrite_close_null(compress_stream_state_t* css_ptr, char* writep) {
  if (is_uncompressed_cswrite(css_ptr)) {
    return uncompressed_cswrite_close_null(css_ptr, writep);
  }
  return compressed_cswrite_close_null(css_ptr, writep);
}

#ifdef __cplusplus
}
#endif
