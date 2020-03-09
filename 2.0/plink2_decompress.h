#ifndef __PLINK2_DECOMPRESS_H__
#define __PLINK2_DECOMPRESS_H__

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


// This has been separated from plink2_cmdline due to the relatively
// heavyweight dependence on zstd.
#include "include/plink2_text.h"
#include "plink2_cmdline.h"

#ifdef __cplusplus
namespace plink2 {
#endif

extern const char kErrprintfDecompress[];

HEADER_INLINE BoolErr CleanupTextFile2(const char* file_descrip, textFILE* txfp, PglErr* reterrp) {
  if (unlikely(CleanupTextFile(txfp, reterrp))) {
    logerrprintfww(kErrprintfFread, file_descrip, strerror(errno));
    return 1;
  }
  return 0;
}

PglErr InitTextStreamEx(const char* fname, uint32_t alloc_at_end, uint32_t enforced_max_line_blen, uint32_t max_line_blen, uint32_t decompress_thread_ct, TextStream* txsp);

HEADER_INLINE PglErr InitTextStream(const char* fname, uint32_t max_line_blen, uint32_t decompress_thread_ct, TextStream* txsp) {
  return InitTextStreamEx(fname, 0, kMaxLongLine, max_line_blen, decompress_thread_ct, txsp);
}

// required_byte_ct can't be greater than kMaxLongLine.
// Now ok for unstandardized_byte_ct to be bigstack_left(), since other
// allocations are made on the heap instead of the arena (to be more usable in
// non-plink2 software).
// Note that the actual buffer size is max_line_blen + kDecompressChunkSize,
// not max_line_blen.
HEADER_INLINE BoolErr StandardizeMaxLineBlenEx(uintptr_t unstandardized_byte_ct, uint32_t required_byte_ct, uint32_t* max_line_blenp) {
#ifdef __LP64__
  if (unstandardized_byte_ct >= S_CAST(uintptr_t, kMaxLongLine) + S_CAST(uintptr_t, kDecompressChunkSize)) {
    *max_line_blenp = kMaxLongLine;
    return 0;
  }
#endif
  if (unlikely(unstandardized_byte_ct < kDecompressChunkSize + RoundUpPow2(MAXV(kDecompressChunkSize, required_byte_ct), kCacheline))) {
    return 1;
  }
  *max_line_blenp = RoundDownPow2(unstandardized_byte_ct, kCacheline) - kDecompressChunkSize;
  return 0;
}

HEADER_INLINE BoolErr StandardizeMaxLineBlen(uintptr_t unstandardized_byte_ct, uint32_t* max_line_blenp) {
  return StandardizeMaxLineBlenEx(unstandardized_byte_ct, kMaxMediumLine + 1, max_line_blenp);
}

HEADER_INLINE PglErr SizeAndInitTextStream(const char* fname, uintptr_t unstandardized_byte_ct, uint32_t decompress_thread_ct, TextStream* txsp) {
  // plink 1.9 immediately failed with an out-of-memory error if a "long line"
  // buffer would be smaller than kMaxMediumLine + 1 bytes, so may as well make
  // that the default lower bound.  (The precise value is currently irrelevant
  // since kTextStreamBlenLowerBound is larger and we take the maximum of the
  // two at compile time, but it's useful to distinguish "minimum acceptable
  // potentially-long-line buffer size" from "load/decompression block size
  // which generally has good performance".)
  uint32_t max_line_blen;
  if (unlikely(StandardizeMaxLineBlen(unstandardized_byte_ct, &max_line_blen))) {
    return kPglRetNomem;
  }
  return InitTextStream(fname, max_line_blen, decompress_thread_ct, txsp);
}

HEADER_INLINE unsigned char* TextStreamMemStart(TextStream* txs_ptr) {
  // placed here instead of plink2_text.h since it's pretty specific to arena
  // memory-management
  return R_CAST(unsigned char*, GET_PRIVATE(*txs_ptr, m).base.dst);
}

// TODO: logputs("\n") first when necessary
void TextErrPrint(const char* file_descrip, const char* errmsg, PglErr reterr);

HEADER_INLINE void TextFileErrPrint(const char* file_descrip, const textFILE* txfp) {
  TextErrPrint(file_descrip, TextFileError(txfp), TextFileErrcode(txfp));
}

HEADER_INLINE void TextStreamErrPrint(const char* file_descrip, const TextStream* txsp) {
  TextErrPrint(file_descrip, TextStreamError(txsp), TextStreamErrcode(txsp));
}

HEADER_INLINE void TextStreamErrPrintRewind(const char* file_descrip, const TextStream* txsp, PglErr* reterrp) {
  if ((*reterrp == kPglRetOpenFail) || (*reterrp == kPglRetEof)) {
    // attempting to rewind/reopen a pipe file descriptor should manifest as
    // one of these two errors.  (todo: verify that open-fail is possible.)
    *reterrp = kPglRetRewindFail;
  }
  if (*reterrp == kPglRetRewindFail) {
    logerrprintfww(kErrprintfRewind, file_descrip);
  } else {
    TextStreamErrPrint(file_descrip, txsp);
  }
}

HEADER_INLINE BoolErr CleanupTextStream2(const char* file_descrip, TextStream* txsp, PglErr* reterrp) {
  if (unlikely(CleanupTextStream(txsp, reterrp))) {
    logerrprintfww(kErrprintfFread, file_descrip, strerror(errno));
    return 1;
  }
  return 0;
}


HEADER_INLINE PglErr InitTokenStreamEx(const char* fname, uint32_t alloc_at_end, uint32_t decompress_thread_ct, TokenStream* tksp) {
  return InitTextStreamEx(fname, alloc_at_end, 0, kTokenStreamBlen - kDecompressChunkSize, decompress_thread_ct, &(tksp->txs));
}

HEADER_INLINE PglErr InitTokenStream(const char* fname, uint32_t decompress_thread_ct, TokenStream* tksp) {
  return InitTextStreamEx(fname, 0, 0, kTokenStreamBlen - kDecompressChunkSize, decompress_thread_ct, &(tksp->txs));
}

HEADER_INLINE void TokenStreamErrPrint(const char* file_descrip, const TokenStream* tksp) {
  PglErr reterr = TokenStreamErrcode(tksp);
  const char* errmsg = TokenStreamError(tksp);
  if (reterr != kPglRetMalformedInput) {
    TextErrPrint(file_descrip, errmsg, reterr);
  } else {
    logerrprintfww("Error: Pathologically long token in %s.\n", file_descrip);
  }
}

HEADER_INLINE BoolErr CleanupTokenStream2(const char* file_descrip, TokenStream* tksp, PglErr* reterrp) {
  if (unlikely(CleanupTokenStream(tksp, reterrp))) {
    logerrprintfww(kErrprintfFread, file_descrip, strerror(errno));
    return 1;
  }
  return 0;
}

HEADER_INLINE BoolErr CleanupTokenStream3(const char* file_descrip, TokenStream* tksp, PglErr* reterrp) {
  *reterrp = kPglRetSuccess;
  return CleanupTokenStream2(file_descrip, tksp, reterrp);
}

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_DECOMPRESS_H__
