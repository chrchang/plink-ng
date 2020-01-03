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

#include <errno.h>
#include "plink2_zstfile.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static inline zstRFILEMain* GetZrfp(zstRFILE* zrf_ptr) {
  return &GET_PRIVATE(*zrf_ptr, m);
}

void PreinitZstRfile(zstRFILE* zrf_ptr) {
  zstRFILEMain* zrfp = GetZrfp(zrf_ptr);
  zrfp->ff = nullptr;
  zrfp->zds = nullptr;
  zrfp->zib.src = nullptr;
  zrfp->errmsg = nullptr;
  zrfp->reterr = kPglRetEof;
}

// zstd prefix_unknown error string
const char kShortErrZstdPrefixUnknown[] = "Unknown frame descriptor";
const char kShortErrZstdCorruptionDetected[] = "Corrupted block detected";

const char kShortErrZstAlreadyOpen[] = "ZstRfileOpen can't be called on an already-open file";

PglErr ZstRfileOpen(const char* fname, zstRFILE* zrf_ptr) {
  zstRFILEMain* zrfp = GetZrfp(zrf_ptr);
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely(zrfp->ff)) {
      reterr = kPglRetImproperFunctionCall;
      zrfp->errmsg = kShortErrZstAlreadyOpen;
      goto ZstRfileOpen_ret_1;
    }
    zrfp->ff = fopen(fname, FOPEN_RB);
    if (unlikely(!zrfp->ff)) {
      goto ZstRfileOpen_ret_OPEN_FAIL;
    }
    zrfp->zds = ZSTD_createDStream();
    if (unlikely(!zrfp->zds)) {
      goto ZstRfileOpen_ret_NOMEM;
    }
    zrfp->zib.src = malloc(ZSTD_DStreamInSize());
    if (!zrfp->zib.src) {
      goto ZstRfileOpen_ret_NOMEM;
    }
    uint32_t nbytes = fread_unlocked(K_CAST(void*, zrfp->zib.src), 1, 4, zrfp->ff);
    if (nbytes < 4) {
      if (!feof_unlocked(zrfp->ff)) {
        goto ZstRfileOpen_ret_READ_FAIL;
      }
      if (!nbytes) {
        // May as well accept this.
        // Don't jump to ret_1 since we're setting zrfp->reterr to a
        // different value than what we're returning.
        zrfp->reterr = kPglRetEof;
        return kPglRetSuccess;
      }
      reterr = kPglRetDecompressFail;
      zrfp->errmsg = kShortErrZstdPrefixUnknown;
      goto ZstRfileOpen_ret_1;
    }
    zrfp->zib.size = 4;
    zrfp->zib.pos = 0;
  }
  while (0) {
  ZstRfileOpen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ZstRfileOpen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    zrfp->errmsg = strerror(errno);
    break;
  ZstRfileOpen_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    zrfp->errmsg = strerror(errno);
    break;
  }
 ZstRfileOpen_ret_1:
  zrfp->reterr = reterr;
  return reterr;
}

int32_t zstread(zstRFILE* zrf_ptr, void* dst, uint32_t len) {
  zstRFILEMain* zrfp = GetZrfp(zrf_ptr);
  if (zrfp->reterr) {
    return 0;
  }
  assert(len < 0x80000000U);
  uint32_t wpos = 0;
  while (wpos != len) {
    ZSTD_outBuffer zob = {&(S_CAST(unsigned char*, dst)[wpos]), len - wpos, 0};
    const uintptr_t read_size_hint = ZSTD_decompressStream(zrfp->zds, &zob, &zrfp->zib);
    if (ZSTD_isError(read_size_hint)) {
      zrfp->reterr = kPglRetDecompressFail;
      zrfp->errmsg = ZSTD_getErrorName(read_size_hint);
      // zlib Z_STREAM_ERROR.  may want to return Z_DATA_ERROR when appropriate
      return -2;
    }
    wpos += zob.pos;
    if (!read_size_hint) {
      // We've finished decompressing the most recent input frame.  Finish
      // loading the next frame if necessary, or exit on eof.
      assert(zrfp->zib.size == zrfp->zib.pos);
      const uint32_t nbytes = fread_unlocked(K_CAST(void*, zrfp->zib.src), 1, 4, zrfp->ff);
      zrfp->zib.size = nbytes;
      zrfp->zib.pos = 0;
      if (nbytes < 4) {
        if (unlikely(!feof_unlocked(zrfp->ff))) {
          // zlib Z_ERRNO
          zrfp->reterr = kPglRetReadFail;
          zrfp->errmsg = strerror(errno);
          return -1;
        }
        if (unlikely(nbytes)) {
          zrfp->reterr = kPglRetDecompressFail;
          zrfp->errmsg = kShortErrZstdPrefixUnknown;
          return -2;
        }
        zrfp->reterr = kPglRetEof;
        break;
      }
      if (unlikely(!IsZstdFrame(*R_CAST(const uint32_t*, zrfp->zib.src)))) {
        zrfp->reterr = kPglRetDecompressFail;
        zrfp->errmsg = kShortErrZstdPrefixUnknown;
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
    if (unlikely(!fread_unlocked(buf, to_decode, 1, zrfp->ff))) {
      if (feof_unlocked(zrfp->ff)) {
        zrfp->reterr = kPglRetDecompressFail;
        zrfp->errmsg = kShortErrZstdCorruptionDetected;
        return -2;
      }
      zrfp->reterr = kPglRetReadFail;
      zrfp->errmsg = strerror(errno);
      return -1;
    }
    zrfp->zib.size = to_decode;
    zrfp->zib.pos = 0;
  }
  return wpos;
}

void zstrewind(zstRFILE* zrf_ptr) {
  zstRFILEMain* zrfp = GetZrfp(zrf_ptr);
  if ((zrfp->reterr == kPglRetSuccess) || (zrfp->reterr == kPglRetEof)) {
    rewind(zrfp->ff);
    ZSTD_DCtx_reset(zrfp->zds, ZSTD_reset_session_only);
    zrfp->zib.size = 0;
    zrfp->zib.pos = 0;
    zrfp->reterr = kPglRetSuccess;
  }
}

BoolErr CleanupZstRfile(zstRFILE* zrf_ptr, PglErr* reterrp) {
  zstRFILEMain* zrfp = GetZrfp(zrf_ptr);
  zrfp->reterr = kPglRetEof;
  zrfp->errmsg = nullptr;
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
      if (!reterrp) {
        return 1;
      }
      if (*reterrp == kPglRetSuccess) {
        *reterrp = kPglRetReadFail;
        return 1;
      }
    }
  }
  return 0;
}

#ifdef __cplusplus
}
#endif
