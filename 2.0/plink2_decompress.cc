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

#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

const char kErrprintfDecompress[] = "Error: %s decompression failure: %s.\n";

PglErr InitTextStreamEx(const char* fname, uint32_t alloc_at_end, uint32_t enforced_max_line_blen, uint32_t max_line_blen, uint32_t decompress_thread_ct, TextStream* txsp) {
  const uint32_t dst_capacity = RoundUpPow2(max_line_blen + kDecompressChunkSize, kCacheline);
  if (unlikely(dst_capacity > bigstack_left())) {
    return kPglRetNomem;
  }
  char* dst;
  if (!alloc_at_end) {
    dst = S_CAST(char*, bigstack_alloc_raw(dst_capacity));
  } else {
    dst = S_CAST(char*, bigstack_end_alloc_raw(dst_capacity));
  }
  return TextStreamOpenEx(fname, enforced_max_line_blen, dst_capacity, decompress_thread_ct, nullptr, dst, txsp);
}

void TextErrPrint(const char* file_descrip, const char* errmsg, PglErr reterr) {
  assert(reterr != kPglRetSuccess);
  if (reterr == kPglRetOpenFail) {
    logerrprintfww(kErrprintfFopen, file_descrip, errmsg);
  } else if (reterr == kPglRetReadFail) {
    logerrprintfww(kErrprintfFread, file_descrip, errmsg);
  } else if (reterr == kPglRetDecompressFail) {
    logerrprintfww(kErrprintfDecompress, file_descrip, errmsg);
  } else if (reterr == kPglRetMalformedInput) {
    if (errmsg == kShortErrInteriorEmptyLine) {
      logerrprintfww("Error: Unexpected interior empty line in %s.\n", file_descrip);
    } else {
      assert(errmsg == kShortErrLongLine);
      logerrprintfww("Error: Pathologically long line in %s.\n", file_descrip);
    }
  } else if (reterr == kPglRetRewindFail) {
    // Not produced directly by TextStream, but it's inserted in between by
    // some consumers.
    logerrprintfww(kErrprintfRewind, file_descrip);
  }
}

#ifdef __cplusplus
}
#endif
