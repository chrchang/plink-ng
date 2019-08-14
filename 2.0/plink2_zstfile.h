#ifndef __PLINK2_READLINE_H__
#define __PLINK2_READLINE_H__

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


// This file defines a simple zstRFILE type, which intentionally does not try
// to offer everything zlib does.  (It doesn't even pass uncompressed data
// through, since we can already use gzFile for that.)  The typename has
// all-caps FILE at the end to indicate that, like FILE and unlike gzFile, it's
// a struct instead of a pointer-to-struct.

#include "plink2_base.h"

#ifdef STATIC_ZSTD
#  include "zstd/lib/zstd.h"
#else
#  include <zstd.h>
#  if !defined(ZSTD_VERSION_NUMBER) || (ZSTD_VERSION_NUMBER < 10401)
#    error "plink2_zstfile requires zstd 1.4.1 or later."
#  endif
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

// may add kFileZstdSeekable later
ENUM_U31_DEF_START()
  kFileUncompressed,
  kFileGzip,
  kFileBgzf,
  kFileZstd
ENUM_U31_DEF_END(FileCompressionType);

PglErr GetFileType(const char* fname, FileCompressionType* ftype_ptr);

// (yes, it might be better to make this opaque at some point.)
typedef struct zstRFILEStruct {
  FILE* ff;
  ZSTD_DStream* zds;
  ZSTD_inBuffer zib;
  uint32_t zstd_eof;
} zstRFILE;

void PreinitZstRfile(zstRFILE* zrfp);

// Assumes file type is already known to be kFileZstd.
// Can return nomem, open-fail, or read-fail.
PglErr ZstRfileOpen(const char* fname, zstRFILE* zrfp);

// This has essentially the same behavior as gzread.  (to determine: can we use
// zlib's negative error codes, or is it necessary to make those
// zstd-specific?)
int32_t zstread(zstRFILE* zrfp, void* dst, uint32_t len);

void zstrewind(zstRFILE* zrfp);

HEADER_INLINE int32_t ZstRfileIsOpen(zstRFILE* zrfp) {
  return (zrfp->ff != nullptr);
}

HEADER_INLINE int32_t zsteof(zstRFILE* zrfp) {
  return zrfp->zstd_eof;
}

// Error should be interpreted as read-fail; errno always set in this case.
BoolErr CleanupZstRfile(zstRFILE* zrfp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_ZSTFILE_H__
