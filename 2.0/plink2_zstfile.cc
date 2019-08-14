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

#include "plink2_zstfile.h"

#ifdef __cplusplus
namespace plink2 {
#endif

uint32_t IsZstdFrame(uint32_t magic4) {
  return (magic4 == ZSTD_MAGICNUMBER) || ((magic4 & ZSTD_MAGIC_SKIPPABLE_MASK) == ZSTD_MAGIC_SKIPPABLE_START);
}

PglErr GetFileType(const char* fname, FileCompressionType* ftype_ptr) {
  FILE* infile = fopen(fname, FOPEN_RB);
  if (unlikely(!infile)) {
    // Note that this does not print an error message (since it may be called
    // by a worker thread).
    return kPglRetOpenFail;
  }
  unsigned char buf[16];
  const uint32_t nbytes = fread_unlocked(buf, 1, 16, infile);
  if (unlikely(ferror_unlocked(infile) || fclose(infile))) {
    return kPglRetReadFail;
  }
  if (nbytes < 4) {
    *ftype_ptr = kFileUncompressed;
    return kPglRetSuccess;
  }
  const uint32_t magic4 = *R_CAST(uint32_t*, buf);
  if (IsZstdFrame(magic4)) {
    *ftype_ptr = kFileZstd;
    return kPglRetSuccess;
  }
  if (S_CAST(uint16_t, magic4) != 0x8b1f) { // gzip ID1/ID2 bytes
    *ftype_ptr = kFileUncompressed;
    return kPglRetSuccess;
  }
  if ((nbytes == 16) && (buf[2] == '\10') && (buf[3] & 4) && memequal_k(&(buf[10]), "\6\0BC\2", 6)) {
    *ftype_ptr = kFileBgzf;
  } else {
    *ftype_ptr = kFileGzip;
  }
  return kPglRetSuccess;
}

void PreinitZstRfile(zstRFILE* zrfp) {
  zrfp->ff = nullptr;
  zrfp->zds = nullptr;
  zrfp->zib.src = nullptr;
  zrfp->zstd_eof = 1;
}

PglErr ZstRfileOpen(const char* fname, zstRFILE* zrfp) {
  if (unlikely(zrfp->ff)) {
    return kPglRetImproperFunctionCall;
  }
  zrfp->ff = fopen(fname, FOPEN_RB);
  if (unlikely(!zrfp->ff)) {
    return kPglRetOpenFail;
  }
  zrfp->zds = ZSTD_createDStream();
  if (unlikely(!zrfp->zds)) {
    return kPglRetNomem;
  }
  zrfp->zib.src = malloc(ZSTD_DStreamInSize());
  if (!zrfp->zib.src) {
    return kPglRetNomem;
  }
  zrfp->zstd_eof = 0;
  uint32_t nbytes = fread_unlocked(K_CAST(void*, zrfp->zib.src), 1, 4, zrfp->ff);
  if (nbytes < 4) {
    return kPglRetReadFail;
  }
  zrfp->zib.size = 4;
  zrfp->zib.pos = 0;
  return kPglRetSuccess;
}

int32_t zstread(zstRFILE* zrfp, void* dst, uint32_t len) {
  if (zrfp->zstd_eof) {
    return 0;
  }
  assert(len < 0x80000000U);
  uint32_t wpos = 0;
  while (wpos != len) {
    ZSTD_outBuffer zob = {&(S_CAST(unsigned char*, dst)[wpos]), len - wpos, 0};
    const uintptr_t read_size_hint = ZSTD_decompressStream(zrfp->zds, &zob, &zrfp->zib);
    if (ZSTD_isError(read_size_hint)) {
      // zlib Z_STREAM_ERROR.  may want to return Z_DATA_ERROR when appropriate
      return -2;
    }
    wpos += zob.pos;
    if (!read_size_hint) {
      // We've finished decompressing the most recent input frame (or are at
      // the beginning of the file).  Finish loading the next frame if
      // necessary, or exit on eof.
      assert(zrfp->zib.size == zrfp->zib.pos);
      const uint32_t nbytes = fread_unlocked(K_CAST(void*, zrfp->zib.src), 1, 4, zrfp->ff);
      zrfp->zib.size = nbytes;
      zrfp->zib.pos = 0;
      if (nbytes < 4) {
        if (unlikely(nbytes)) {
          return -2;
        }
        if (unlikely(!feof_unlocked(zrfp->ff))) {
          // zlib Z_ERRNO
          return -1;
        }
        zrfp->zstd_eof = 1;
        break;
      }
      if (unlikely(!IsZstdFrame(*R_CAST(const uint32_t*, zrfp->zib.src)))) {
        return -2;
      }
      // impossible for this to fail
      ZSTD_DCtx_reset(zrfp->zds, ZSTD_reset_session_only);
      continue;
    }
    if (zrfp->zib.pos != zrfp->zib.size) {
      assert(wpos == len);
      break;
    }
    const uint32_t to_decode = MINV(read_size_hint, ZSTD_DStreamInSize());
    unsigned char* buf = S_CAST(unsigned char*, K_CAST(void*, zrfp->zib.src));
    const uint32_t nbytes = fread_unlocked(buf, 1, to_decode, zrfp->ff);
    if (unlikely(!nbytes)) {
      if (feof_unlocked(zrfp->ff)) {
        return -2;
      }
      return -1;
    }
    zrfp->zib.size = nbytes;
    zrfp->zib.pos = 0;
  }
  return wpos;
}

void zstrewind(zstRFILE* zrfp) {
  rewind(zrfp->ff);
  ZSTD_DCtx_reset(zrfp->zds, ZSTD_reset_session_only);
  zrfp->zib.size = 0;
  zrfp->zib.pos = 0;
  zrfp->zstd_eof = 0;
}

BoolErr CleanupZstRfile(zstRFILE* zrfp) {
  if (zrfp->zib.src) {
    free_const(zrfp->zib.src);
    zrfp->zib.src = nullptr;
  }
  if (zrfp->zds) {
    ZSTD_freeDStream(zrfp->zds); // this should never fail in practice
    zrfp->zds = nullptr;
  }
  if (zrfp->ff) {
    if (unlikely(fclose_null(&zrfp->ff))) {
      return 1;
    }
  }
  return 0;
}

#ifdef __cplusplus
}
#endif
