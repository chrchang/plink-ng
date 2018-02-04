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


#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

pglerr_t gzopen_read_checked(const char* fname, gzFile* gzf_ptr) {
  *gzf_ptr = gzopen(fname, FOPEN_RB);
  if (!(*gzf_ptr)) {
    logprint("\n");
    LOGERRPRINTFWW(g_errstr_fopen, fname);
    return kPglRetOpenFail;
  }
  if (gzbuffer(*gzf_ptr, 131072)) {
    return kPglRetNomem;
  }
  return kPglRetSuccess;
}

pglerr_t gzopen_and_skip_first_lines(const char* fname, uint32_t lines_to_skip, uintptr_t loadbuf_size, char* loadbuf, gzFile* gzf_ptr) {
  pglerr_t reterr = gzopen_read_checked(fname, gzf_ptr);
  if (reterr) {
    return reterr;
  }
  loadbuf[loadbuf_size - 1] = ' ';
  for (uint32_t line_idx = 1; line_idx <= lines_to_skip; ++line_idx) {
    if (!gzgets(*gzf_ptr, loadbuf, loadbuf_size)) {
      if (gzeof(*gzf_ptr)) {
        LOGERRPRINTFWW("Error: Fewer lines than expected in %s.\n", fname);
        return kPglRetInconsistentInput;
      }
      return kPglRetReadFail;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if ((loadbuf_size == kMaxMediumLine) || (loadbuf_size == kMaxLongLine)) {
        LOGERRPRINTFWW("Error: Line %u of %s is pathologically long.\n", line_idx, fname);
        return kPglRetMalformedInput;
      }
      return kPglRetNomem;
    }
  }
  return kPglRetSuccess;
}


void gz_token_stream_preinit(gz_token_stream_t* gtsp) {
  gtsp->gz_infile = nullptr;
}

pglerr_t gz_token_stream_init(const char* fname, gz_token_stream_t* gtsp, char* buf_start) {
  pglerr_t reterr = gzopen_read_checked(fname, &(gtsp->gz_infile));
  if (reterr) {
    return reterr;
  }
  gtsp->buf_start = buf_start;
  gtsp->read_iter = &(buf_start[kMaxMediumLine]);
  gtsp->buf_end = gtsp->read_iter;
  gtsp->buf_end[0] = '0'; // force initial load
  return kPglRetSuccess;
}

char* gz_token_stream_advance(gz_token_stream_t* gtsp, uint32_t* token_slen_ptr) {
  char* token_start = gtsp->read_iter;
  char* buf_end = gtsp->buf_end;
 gz_token_stream_advance_restart:
  while (ctou32(*token_start) <= ' ') {
    ++token_start;
  }
  while (1) {
    if (token_start < buf_end) {
      char* token_end = &(token_start[1]);
      while (ctou32(*token_end) > ' ') {
        ++token_end;
      }
      const uint32_t token_slen = token_end - token_start;
      if (token_end < buf_end) {
        *token_slen_ptr = token_slen;
        gtsp->read_iter = &(token_end[1]);
        return token_start;
      }
      if (token_slen > kMaxMediumLine) {
        *token_slen_ptr = UINT32_MAX;
        return nullptr;
      }
      char* new_token_start = &(gtsp->buf_start[kMaxMediumLine - token_slen]);
      memcpy(new_token_start, token_start, token_slen);
      token_start = new_token_start;
    } else {
      token_start = &(gtsp->buf_start[kMaxMediumLine]);
    }
    char* load_start = &(gtsp->buf_start[kMaxMediumLine]);
    const int32_t bufsize = gzread(gtsp->gz_infile, load_start, kMaxMediumLine);
    if (bufsize < 0) {
      *token_slen_ptr = 0xfffffffeU;
      return nullptr;
    }
    buf_end = &(load_start[S_CAST(uint32_t, bufsize)]);
    buf_end[0] = ' ';
    buf_end[1] = '0';
    gtsp->buf_end = buf_end;
    if (!bufsize) {
      if (!gzeof(gtsp->gz_infile)) {
        *token_slen_ptr = 0xfffffffeU;
        return nullptr;
      }
      // bufsize == 0, eof
      if (token_start == load_start) {
        *token_slen_ptr = 0;
        return nullptr;
      }
      gtsp->read_iter = load_start;
      *token_slen_ptr = load_start - token_start;
      return token_start;
    }
    if (token_start == load_start) {
      goto gz_token_stream_advance_restart;
    }
  }
}

boolerr_t gz_token_stream_close(gz_token_stream_t* gtsp) {
  if (!gtsp->gz_infile) {
    return 0;
  }
  return gzclose_null(&(gtsp->gz_infile));
}


#ifdef __cplusplus
}
#endif
