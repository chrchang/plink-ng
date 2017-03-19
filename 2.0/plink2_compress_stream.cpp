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


#include "plink2_compress_stream.h"

#ifdef __cplusplus
namespace plink2 {
#endif

boolerr_t uncompressed_cswrite_init(const char* out_fname, uint32_t do_append, unsigned char* overflow_buf, compress_stream_state_t* css_ptr) {
  css_ptr->z_outfile = nullptr;
  // can't use fopen_checked since we need to be able to append
  css_ptr->outfile = fopen(out_fname, do_append? FOPEN_AB : FOPEN_WB);
  if (!css_ptr->outfile) {
    logprint("\n");
    LOGERRPRINTFWW(g_errstr_fopen, out_fname);
    return 1;
  }
  css_ptr->overflow_buf = overflow_buf;
  return 0;
}

boolerr_t zstd_cswrite_init(const char* out_fname, uint32_t do_append, unsigned char* overflow_buf, compress_stream_state_t* css_ptr) {
  css_ptr->outfile = nullptr;
  css_ptr->z_outfile = gzopen(out_fname, do_append? FOPEN_AB : FOPEN_WB);
  if (!css_ptr->z_outfile) {
    logprint("\n");
    LOGERRPRINTFWW(g_errstr_fopen, out_fname);
    return 1;
  }
  css_ptr->overflow_buf = overflow_buf;
  return 0;
}

// possible todo: replace output_zst with an enum which permits gzipping
boolerr_t cswrite_init(const char* out_fname, uint32_t do_append, uint32_t output_zst, unsigned char* overflow_buf, compress_stream_state_t* css_ptr) {
  if (!output_zst) {
    return uncompressed_cswrite_init(out_fname, do_append, overflow_buf, css_ptr);
  } else {
    return zstd_cswrite_init(out_fname, do_append, overflow_buf, css_ptr);
  }
}

boolerr_t force_uncompressed_cswrite(compress_stream_state_t* css_ptr, char** writep_ptr) {
  unsigned char* writep = (unsigned char*)(*writep_ptr);
  if (css_ptr->overflow_buf != writep) {
    if (!fwrite(css_ptr->overflow_buf, writep - css_ptr->overflow_buf, 1, css_ptr->outfile)) {
      return 1;
    }
    *writep_ptr = (char*)(css_ptr->overflow_buf);
  }
  return 0;
}

boolerr_t force_compressed_cswrite(__attribute__((unused)) uint32_t write_min, compress_stream_state_t* css_ptr, char** writep_ptr) {
  unsigned char* writep = (unsigned char*)(*writep_ptr);
  if (css_ptr->overflow_buf != writep) {
    if (!gzwrite(css_ptr->z_outfile, css_ptr->overflow_buf, writep - css_ptr->overflow_buf)) {
      return 1;
    }
    *writep_ptr = (char*)(css_ptr->overflow_buf);
  }
  return 0;
}

boolerr_t csputs_std(const char* ss, uint32_t sslen, compress_stream_state_t* css_ptr, char** writep_ptr) {
  unsigned char* writep = (unsigned char*)(*writep_ptr);
  const unsigned char* readp = (const unsigned char*)ss;
  uint32_t cur_write_space = 2 * kCompressStreamBlock - ((uintptr_t)(writep - css_ptr->overflow_buf));
  while (sslen > cur_write_space) {
    memcpy(writep, readp, cur_write_space);
    if (is_uncompressed_cswrite(css_ptr)) {
      if (!fwrite(css_ptr->overflow_buf, 2 * kCompressStreamBlock, 1, css_ptr->outfile)) {
	return 1;
      }
    } else {
      if (!gzwrite(css_ptr->z_outfile, css_ptr->overflow_buf, 2 * kCompressStreamBlock)) {
	return 1;
      }
    }
    writep = css_ptr->overflow_buf;
    readp = &(readp[cur_write_space]);
    sslen -= cur_write_space;
    cur_write_space = 2 * kCompressStreamBlock;
  }
  memcpy(writep, readp, sslen);
  *writep_ptr = (char*)(&(writep[sslen]));
  return cswrite(css_ptr, writep_ptr);
}

boolerr_t uncompressed_cswrite_close_null(compress_stream_state_t* css_ptr, char* writep) {
  force_uncompressed_cswrite(css_ptr, &writep);
  css_ptr->overflow_buf = nullptr;
  int32_t ii = ferror(css_ptr->outfile);
  int32_t jj = fclose(css_ptr->outfile);
  return ii || jj;
}

boolerr_t compressed_cswrite_close_null(compress_stream_state_t* css_ptr, char* writep) {
  force_compressed_cswrite(0, css_ptr, &writep);
  css_ptr->overflow_buf = nullptr;
  return (gzclose(css_ptr->z_outfile) != Z_OK);
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
