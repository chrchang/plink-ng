#ifndef __PLINK2_COMPRESS_STREAM_H__
#define __PLINK2_COMPRESS_STREAM_H__

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


// Successor to plink 1.9 pigz.h.  Provides a basic manually-buffered output
// stream interface for zstd compression.  Not multithreaded yet, but the
// interface is identical to the old multithreaded gzipper so we'll be able to
// upgrade the backend later without making significant changes to other code.

#include "include/plink2_zstfile.h"
#include "plink2_cmdline.h"

#ifdef __cplusplus
namespace plink2 {
#endif

extern uint32_t g_zst_level;

// This should be at least as large as zstd's internal block size.
// todo: test different values, may want to increase on at least OS X...
CONSTI32(kCompressStreamBlock, 131072);

typedef struct CompressStreamStateStruct {
  NONCOPYABLE(CompressStreamStateStruct);
  // Usually compress text, so appropriate to define this as char*.
  char* overflow_buf;

  FILE* outfile;
  ZSTD_CCtx* cctx;
  ZSTD_outBuffer output;
} CompressStreamState;

HEADER_INLINE uint32_t IsUncompressedCstream(const CompressStreamState* css_ptr) {
  return (css_ptr->cctx == nullptr);
}

HEADER_INLINE void PreinitCstream(CompressStreamState* css_ptr) {
  css_ptr->overflow_buf = nullptr;
}

PglErr InitCstreamNoop(const char* out_fname, uint32_t do_append, char* overflow_buf, CompressStreamState* css_ptr);

HEADER_INLINE uintptr_t CstreamWkspaceReq(uintptr_t overflow_buf_size) {
  return kCompressStreamBlock + MAXV(ZSTD_compressBound(overflow_buf_size), ZSTD_CStreamOutSize());
}

// overflow_buf must have space for at least kCompressStreamBlock + [max bytes
// added between cswrite() calls] bytes.
// compress_wkspace can be nullptr in no-compression case; otherwise it must
// have space for cswrite_wkspace_req(overflow_buf size) bytes.
PglErr InitCstream(const char* out_fname, uint32_t do_append, uint32_t output_zst, uint32_t thread_ct, uintptr_t overflow_buf_size, char* overflow_buf, unsigned char* compress_wkspace, CompressStreamState* css_ptr);

// Convenience interface which allocates from the bottom of g_bigstack.
PglErr InitCstreamAlloc(const char* out_fname, uint32_t do_append, uint32_t output_zst, uint32_t thread_ct, uintptr_t overflow_buf_size, CompressStreamState* css_ptr, char** cswritepp);

BoolErr ForceUncompressedCswrite(CompressStreamState* css_ptr, char** writep_ptr);

// No longer guaranteed to consume entire input buffer, only reduces it to
// <128k.
BoolErr ForceCompressedCswrite(CompressStreamState* css_ptr, char** writep_ptr);

HEADER_INLINE BoolErr Cswrite(CompressStreamState* css_ptr, char** writep_ptr) {
  if (S_CAST(uintptr_t, (*writep_ptr) - css_ptr->overflow_buf) >= kCompressStreamBlock + 1) {
    if (IsUncompressedCstream(css_ptr)) {
      return ForceUncompressedCswrite(css_ptr, writep_ptr);
    } else {
      return ForceCompressedCswrite(css_ptr, writep_ptr);
    }
  }
  return 0;
}

// assumes overflow_buf has size >= 2 * kCompressStreamBlock.
BoolErr CsputsStd(const char* readp, uint32_t byte_ct, CompressStreamState* css_ptr, char** writep_ptr);

BoolErr UncompressedCswriteCloseNull(CompressStreamState* css_ptr, char* writep);

BoolErr CompressedCswriteCloseNull(CompressStreamState* css_ptr, char* writep);

BoolErr CswriteCloseNull(CompressStreamState* css_ptr, char* writep);

HEADER_INLINE void UncompressedCswriteCloseCond(CompressStreamState* css_ptr, char* writep) {
  if (css_ptr->overflow_buf) {
    UncompressedCswriteCloseNull(css_ptr, writep);
  }
}

HEADER_INLINE void CompressedCswriteCloseCond(CompressStreamState* css_ptr, char* writep) {
  if (css_ptr->overflow_buf) {
    CompressedCswriteCloseNull(css_ptr, writep);
  }
}

HEADER_INLINE void CswriteCloseCond(CompressStreamState* css_ptr, char* writep) {
  if (css_ptr->overflow_buf) {
    CswriteCloseNull(css_ptr, writep);
  }
}


#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_COMPRESS_STREAM_H__
