#ifndef __PLINK2_COMPRESS_STREAM_H__
#define __PLINK2_COMPRESS_STREAM_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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
#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// todo: test different values, may want to increase on at least OS X...
CONSTU31(kCompressStreamBlock, 131072);

typedef struct {
  unsigned char* overflow_buf;
  FILE* outfile;
  gzFile z_outfile;
} compress_stream_state_t;

HEADER_INLINE uint32_t is_uncompressed_cswrite(const compress_stream_state_t* css_ptr) {
  return (css_ptr->outfile != nullptr);
}

HEADER_INLINE void cswrite_init_null(compress_stream_state_t* css_ptr) {
  css_ptr->overflow_buf = nullptr;
}

boolerr_t uncompressed_cswrite_init(const char* out_fname, uint32_t do_append, unsigned char* overflow_buf, compress_stream_state_t* css_ptr);

// overflow_buf must have space for at least kCompressStreamBlock + [max bytes
// added between cswrite() calls] bytes.
boolerr_t cswrite_init(const char* out_fname, uint32_t do_append, uint32_t output_zst, unsigned char* overflow_buf, compress_stream_state_t* css_ptr);

boolerr_t force_uncompressed_cswrite(compress_stream_state_t* css_ptr, char** writep_ptr);

// may or may not need the write_min parameter; let's leave it until we
// actually try to implement multithreaded zstd
boolerr_t force_compressed_cswrite(uint32_t write_min, compress_stream_state_t* css_ptr, char** writep_ptr);

HEADER_INLINE boolerr_t cswrite(compress_stream_state_t* css_ptr, char** writep_ptr) {
  if ((uintptr_t)(((unsigned char*)(*writep_ptr)) - css_ptr->overflow_buf) >= kCompressStreamBlock + 1) {
    if (is_uncompressed_cswrite(css_ptr)) {
      return force_uncompressed_cswrite(css_ptr, writep_ptr);
    } else {
      return force_compressed_cswrite(kCompressStreamBlock + 1, css_ptr, writep_ptr);
    }
  }
  return 0;
}

// assumes overflow_buf has size >= 2 * kCompressStreamBlock.
boolerr_t csputs_std(const char* ss, uint32_t sslen, compress_stream_state_t* css_ptr, char** writep_ptr);

boolerr_t uncompressed_cswrite_close_null(compress_stream_state_t* css_ptr, char* writep);

boolerr_t compressed_cswrite_close_null(compress_stream_state_t* css_ptr, char* writep);

boolerr_t cswrite_close_null(compress_stream_state_t* css_ptr, char* writep);

HEADER_INLINE void uncompressed_cswrite_close_cond(compress_stream_state_t* css_ptr, char* writep) {
  if (css_ptr->overflow_buf) {
    uncompressed_cswrite_close_null(css_ptr, writep);
  }
}

HEADER_INLINE void compressed_cswrite_close_cond(compress_stream_state_t* css_ptr, char* writep) {
  if (css_ptr->overflow_buf) {
    compressed_cswrite_close_null(css_ptr, writep);
  }
}

HEADER_INLINE void cswrite_close_cond(compress_stream_state_t* css_ptr, char* writep) {
  if (css_ptr->overflow_buf) {
    cswrite_close_null(css_ptr, writep);
  }
}


#ifdef __cplusplus
} // namespace plink2
#endif
 
#endif // __PLINK2_COMPRESS_STREAM_H__
