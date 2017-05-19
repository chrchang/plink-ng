#ifndef __PLINK2_DECOMPRESS_H__
#define __PLINK2_DECOMPRESS_H__

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


// This has been separated from plink2_common due to the relatively heavyweight
// dependence on zstd/zlibWrapper.
#include "plink2_common.h"

// documentation on ZWRAP_USE_ZSTD is incorrect as of 11 Jan 2017, necessary to
// edit zstd_zlibwrapper.c or use compile flag.
#include "zstd/zlibWrapper/zstd_zlibwrapper.h"
#ifndef STATIC_ZLIB
  #if !defined(ZLIB_VERNUM) || ZLIB_VERNUM < 0x1240
    #error "plink2_decompress requires zlib 1.2.4 or later."
  #endif
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

// Also sets 128k read buffer.
pglerr_t gzopen_read_checked(const char* fname, gzFile* gzf_ptr);

// plink2_compress_stream interface should be used for writing .gz files.

HEADER_INLINE boolerr_t gzclose_null(gzFile* gzf_ptr) {
  const int32_t ii = gzclose(*gzf_ptr);
  *gzf_ptr = nullptr;
  return (ii != Z_OK);
}

HEADER_INLINE void gzclose_cond(gzFile gz_infile) {
  if (gz_infile) {
    gzclose(gz_infile);
  }
}


// may return kPglRetLongLine or kPglRetEmptyFile
// loadbuf_iter_ptr can be nullptr
// line_idx must be zero unless initial lines were skipped
pglerr_t load_xid_header(const char* flag_name, sid_detect_mode_t sid_detect_mode, uintptr_t loadbuf_size, char* loadbuf, char** loadbuf_iter_ptr, uintptr_t* line_idx_ptr, char** loadbuf_first_token_ptr, gzFile* gz_infile_ptr, xid_mode_t* xid_mode_ptr);

// sets last character of loadbuf to ' '
pglerr_t open_and_load_xid_header(const char* fname, const char* flag_name, sid_detect_mode_t sid_detect_mode, uintptr_t loadbuf_size, char* loadbuf, char** loadbuf_iter_ptr, uintptr_t* line_idx_ptr, char** loadbuf_first_token_ptr, gzFile* gz_infile_ptr, xid_mode_t* xid_mode_ptr);


// currently hardcoded to have maximum token length = kMaxMediumLine, buffer
// size = 2 * kMaxMediumLine * 2.
typedef struct gz_token_stream_struct {
  gzFile gz_infile;
  char* buf_start;
  char* read_iter;
  char* buf_end;
} gz_token_stream_t;

void gz_token_stream_preinit(gz_token_stream_t* gtsp);

pglerr_t gz_token_stream_init(const char* fname, gz_token_stream_t* gtsp, char* buf_start);

// sets token_slen to 0xfffffffeU on read fail, 0xffffffffU on too-long token
// safe to null-terminate token between calls
char* gz_token_stream_advance(gz_token_stream_t* gtsp, uint32_t* token_slen_ptr);

// ok if already closed
boolerr_t gz_token_stream_close(gz_token_stream_t* gtsp);

#ifdef __cplusplus
} // namespace plink2
#endif
 
#endif // __PLINK2_DECOMPRESS_H__
