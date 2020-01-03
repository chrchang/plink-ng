#ifndef __PLINK2_ZSTFILE_H__
#define __PLINK2_ZSTFILE_H__

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


// This file defines a simple zstRFILE type, which intentionally does not try
// to offer everything zlib does.  (It doesn't even pass uncompressed data
// through, since we can already use gzFile for that.)  The typename has
// all-caps FILE at the end to indicate that, like FILE and unlike gzFile, it's
// a struct instead of a pointer-to-struct.

#include "plink2_base.h"

#ifdef STATIC_ZSTD
#  include "../zstd/lib/zstd.h"
#else
#  include <zstd.h>
#  if !defined(ZSTD_VERSION_NUMBER) || (ZSTD_VERSION_NUMBER < 10404)
#    error "plink2_zstfile requires zstd 1.4.4 or later."
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

HEADER_INLINE uint32_t IsZstdFrame(uint32_t magic4) {
  return (magic4 == ZSTD_MAGICNUMBER) || ((magic4 & ZSTD_MAGIC_SKIPPABLE_MASK) == ZSTD_MAGIC_SKIPPABLE_START);
}

typedef struct zstRFILEMainStruct {
  FILE* ff;
  ZSTD_DStream* zds;
  ZSTD_inBuffer zib;
  const char* errmsg;
  PglErr reterr;  // kPglRetSkipped == ordinary eof
} zstRFILEMain;

typedef struct zstRFILEStruct {
#ifdef __cplusplus
  zstRFILEMain& GET_PRIVATE_m() { return m; }
  zstRFILEMain const& GET_PRIVATE_m() const { return m; }
 private:
#endif
  zstRFILEMain m;
} zstRFILE;

void PreinitZstRfile(zstRFILE* zrf_ptr);

extern const char kShortErrZstdPrefixUnknown[];
extern const char kShortErrZstdCorruptionDetected[];

// Assumes file type is already known to be kFileZstd.
// Can return nomem, open-fail, or read-fail.
PglErr ZstRfileOpen(const char* fname, zstRFILE* zrf_ptr);

// This has essentially the same behavior as gzread.  (to determine: can we use
// zlib's negative error codes, or is it necessary to make those
// zstd-specific?)
int32_t zstread(zstRFILE* zrf_ptr, void* dst, uint32_t len);

void zstrewind(zstRFILE* zrf_ptr);

HEADER_INLINE int32_t ZstRfileIsOpen(const zstRFILE* zrf_ptr) {
  return (GET_PRIVATE(*zrf_ptr, m).ff != nullptr);
}

HEADER_INLINE int32_t zsteof(const zstRFILE* zrf_ptr) {
  return (GET_PRIVATE(*zrf_ptr, m).reterr == kPglRetEof);
}

HEADER_INLINE const char* zsterror(const zstRFILE* zrf_ptr) {
  return GET_PRIVATE(*zrf_ptr, m).errmsg;
}

HEADER_INLINE PglErr ZstRfileErrcode(const zstRFILE* zrf_ptr) {
  const PglErr reterr = GET_PRIVATE(*zrf_ptr, m).reterr;
  if (reterr == kPglRetEof) {
    return kPglRetSuccess;
  }
  return reterr;
}

// If you need clearerr functionality, call CleanupZstRfile().  (Could provide
// a ZstRewindAndClearerr function too.  Don't think it's worth the trouble to
// support any other clearerr use patterns, though.)

// Ok to pass reterrp == nullptr.
// Returns nonzero iff file-close fails, and either reterrp == nullptr or
// *reterrp == kPglRetSuccess.  In the latter case, *reterrp is set to
// kPglRetReadFail.  (Note that this does *not* retrieve the existing
// zrfp->reterr value; caller is responsible for checking ZstRfileErrcode()
// first when they care.)
BoolErr CleanupZstRfile(zstRFILE* zrf_ptr, PglErr* reterrp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_ZSTFILE_H__
