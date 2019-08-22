#ifndef __PLINK2_ZSTFILE_H__
#define __PLINK2_ZSTFILE_H__

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

HEADER_INLINE int32_t IsBgzfHeader(const void* buf) {
  const uint32_t magic4 = *S_CAST(const uint32_t*, buf);
  return ((magic4 & 0x4ffffff) == 0x4088b1f) && memequal_k(&(S_CAST(const unsigned char*, buf)[10]), "\6\0BC\2", 6);
}

HEADER_INLINE uint32_t IsZstdFrame(uint32_t magic4) {
  return (magic4 == ZSTD_MAGICNUMBER) || ((magic4 & ZSTD_MAGIC_SKIPPABLE_MASK) == ZSTD_MAGIC_SKIPPABLE_START);
}

PglErr GetFileType(const char* fname, FileCompressionType* ftype_ptr);

// (yes, it might be better to make this opaque at some point.)
typedef struct zstRFILEStruct {
  NONCOPYABLE(zstRFILEStruct);
  FILE* ff;
  ZSTD_DStream* zds;
  ZSTD_inBuffer zib;
  const char* errmsg;
  PglErr reterr;  // kPglRetSkipped == ordinary eof
} zstRFILE;

void PreinitZstRfile(zstRFILE* zrfp);

extern const char kShortErrZstdPrefixUnknown[];
extern const char kShortErrZstdCorruptionDetected[];

// Assumes file type is already known to be kFileZstd.
// Can return nomem, open-fail, or read-fail.
PglErr ZstRfileOpen(const char* fname, zstRFILE* zrfp);

// This has essentially the same behavior as gzread.  (to determine: can we use
// zlib's negative error codes, or is it necessary to make those
// zstd-specific?)
int32_t zstread(zstRFILE* zrfp, void* dst, uint32_t len);

void zstrewind(zstRFILE* zrfp);

HEADER_INLINE int32_t ZstRfileIsOpen(const zstRFILE* zrfp) {
  return (zrfp->ff != nullptr);
}

HEADER_INLINE int32_t zsteof(const zstRFILE* zrfp) {
  return (zrfp->reterr == kPglRetEof);
}

HEADER_INLINE const char* zsterror(const zstRFILE* zrfp) {
  return zrfp->errmsg;
}

HEADER_INLINE PglErr ZstRfileErrcode(const zstRFILE* zrfp) {
  if (zrfp->reterr == kPglRetEof) {
    return kPglRetSuccess;
  }
  return zrfp->reterr;
}

// If you need clearerr functionality, call CleanupZstRfile().  (Could provide
// a ZstRewindAndClearerr function too.  Don't think it's worth the trouble to
// support any other clearerr use patterns, though.)

// Ok to pass reterrp == nullptr.
// Returns nonzero iff file-close fails, and either reterrp == nullptr or
// *reterrp == kPglRetSuccess.  In the latter case, *reterrp is set to
// kPglRetReadFail.
BoolErr CleanupZstRfile(zstRFILE* zrfp, PglErr* reterrp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_ZSTFILE_H__
