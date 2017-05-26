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


pglerr_t load_xid_header(const char* flag_name, sid_detect_mode_t sid_detect_mode, uintptr_t loadbuf_size, char* loadbuf, char** loadbuf_iter_ptr, uintptr_t* line_idx_ptr, char** loadbuf_first_token_ptr, gzFile* gz_infile_ptr, xid_mode_t* xid_mode_ptr) {
  // possible todo: support comma delimiter
  uintptr_t line_idx = *line_idx_ptr;
  uint32_t is_header_line;
  char* loadbuf_first_token;
  do {
    ++line_idx;
    if (!gzgets(*gz_infile_ptr, loadbuf, loadbuf_size)) {
      if (!gzeof(*gz_infile_ptr)) {
	return kPglRetReadFail;
      }
      return kPglRetEmptyFile;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      return kPglRetLongLine;
    }
    loadbuf_first_token = skip_initial_spaces(loadbuf);
    is_header_line = (loadbuf_first_token[0] == '#');
  } while (is_header_line && strcmp_se(&(loadbuf_first_token[1]), "FID", 3) && strcmp_se(&(loadbuf_first_token[1]), "IID", 3));
  xid_mode_t xid_mode = kfXidMode0;
  char* loadbuf_iter;
  if (is_header_line) {
    // the following header leading columns are supported:
    // #FID IID (sid_detect_mode can't be FORCE)
    // #FID IID SID (SID ignored on sid_detect_mode NOT_LOADED)
    // #IID
    // #IID SID
    loadbuf_iter = &(loadbuf_first_token[4]);
    if (loadbuf_first_token[1] == 'I') {
      xid_mode = kfXidModeFlagNeverFid;
    } else {
      loadbuf_iter = skip_initial_spaces(loadbuf_iter);
      if (strcmp_se(loadbuf_iter, "IID", 3)) {
	LOGERRPRINTF("Error: No IID column on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
	return kPglRetMalformedInput;
      }
      loadbuf_iter = &(loadbuf_iter[3]);
    }
    loadbuf_iter = skip_initial_spaces(loadbuf_iter);
    if (!strcmp_se(loadbuf_iter, "SID", 3)) {
      if ((uint32_t)sid_detect_mode >= kSidDetectModeLoaded) {
	xid_mode |= kfXidModeFlagSid;
      }
      loadbuf_iter = skip_initial_spaces(&(loadbuf_iter[3]));
    } else if (sid_detect_mode == kSidDetectModeForce) {
      LOGERRPRINTFWW("Error: No SID column on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
      return kPglRetMalformedInput;
    }
  } else {
    xid_mode = (sid_detect_mode == kSidDetectModeForce)? kfXidModeFidiidSid : kfXidModeFidiidOrIid;
    loadbuf_iter = loadbuf_first_token;
  }
  if (loadbuf_iter_ptr) {
    *loadbuf_iter_ptr = loadbuf_iter;
  }
  *loadbuf_first_token_ptr = loadbuf_first_token;
  *line_idx_ptr = line_idx;
  *xid_mode_ptr = xid_mode;
  return kPglRetSuccess;
}

pglerr_t open_and_load_xid_header(const char* fname, const char* flag_name, sid_detect_mode_t sid_detect_mode, uintptr_t loadbuf_size, char* loadbuf, char** loadbuf_iter_ptr, uintptr_t* line_idx_ptr, char** loadbuf_first_token_ptr, gzFile* gz_infile_ptr, xid_mode_t* xid_mode_ptr) {
  pglerr_t reterr = gzopen_read_checked(fname, gz_infile_ptr);
  if (reterr) {
    return reterr;
  }
  loadbuf[loadbuf_size - 1] = ' ';
  return load_xid_header(flag_name, sid_detect_mode, loadbuf_size, loadbuf, loadbuf_iter_ptr, line_idx_ptr, loadbuf_first_token_ptr, gz_infile_ptr, xid_mode_ptr);
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
  while ((unsigned char)(*token_start) <= ' ') {
    ++token_start;
  }
  while (1) {
    if (token_start < buf_end) {
      char* token_end = &(token_start[1]);
      while ((unsigned char)(*token_end) > ' ') {
	++token_end;
      }
      const uint32_t token_slen = (uintptr_t)(token_end - token_start);
      if (token_end < buf_end) {
	*token_slen_ptr = token_slen;
	gtsp->read_iter = &(token_end[1]);
	return token_start;
      }
      if (token_slen > kMaxMediumLine) {
	*token_slen_ptr = 0xffffffffU;
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
    buf_end = &(load_start[(uint32_t)bufsize]);
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
      *token_slen_ptr = (uintptr_t)(load_start - token_start);
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
