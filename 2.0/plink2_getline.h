#ifndef __PLINK2_GETLINE_H__
#define __PLINK2_GETLINE_H__

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


// Scanning one line at a time from a text file is one of the most common
// workflows in all of computing.
//
// Usually, text files are small; the obvious reason to choose text over binary
// is human-readability, after all, and humans can't read multi-gigabyte files
// in a reasonable amount of time.  As a consequence, the most commonly used C
// and C++ text-processing library functions sacrifice a substantial amount of
// performance in favor of ease-of-use.
//
// However, plink2 is frequently asked to load a multi-gigabyte text file and
// then do something very simple with it.  Often, the file is in the operating
// system's page cache, since the user or script is doing multiple things with
// the file and they're split across multiple invocations of plink2 and other
// programs.  In this setting, the usual "I/O cost > processing cost, it isn't
// worth worrying much about the latter" assumption is very, very wrong, and it
// is worth going to great lengths to keep baseline text-processing cost to a
// minimum.
//
// In addition, multi-gigabyte text files are practically guaranteed to
// compress well, and gzipped and bgzipped text files are widely used in
// bioinformatics practice.  Ordinarily, when sequentially processing a text
// file, there's little to gain from spawning a separate thread to issue
// file-read requests, since a modern operating system will recognize the
// access pattern and read-ahead from the disk on its own.  However, the
// operating system can't *decompress-ahead* for you; and when decompression
// has comparable latency to processing, decompress-ahead reduces runtime by up
// to 50%.
//
// Thus, this library provides a text reader that
// 1. allows the caller to treat gzipped and Zstd-compressed text files as if
//    they were uncompressed.  This is functionally identical to Zstd's
//    zlibWrapper, which was used for most of plink2's alpha testing period,
//    but I've decided to phase out zlibWrapper thanks to compilation headaches
//    and its static-linking requirement.
// 2. decompresses-ahead, potentially with multiple threads.
//    a. For now, multithreaded decompression can only kick in for bgzipped
//       files.  However, if a clear use case exists, it should be possible to
//       build a multithreaded Zstd decoder that isn't restricted to a Zstd
//       sub-format; see
//         https://github.com/facebook/zstd/issues/1702#issuecomment-515124700
//       If you have such a use case, I recommend responding to Yann Collet in
//       that GitHub issue.
//    b. Tabix-based seek support was considered and rejected, since the tabix
//       index only stores CHROM/POS, while plink2 also needs record numbers in
//       its most critical use case (.pvar loading).  A suitable index format
//       may be adopted or developed later, and if/when that's supported, there
//       will probably also be a way for callers which don't need record
//       numbers to exploit tabix indexes.  (In that event, plink2's Zstd
//       compressor will be modified to support seekable-Zstd output and index
//       generation.)
// 3. has line-reader functions that don't force the user to provide their own
//    buffer to put the line in.  Instead, they just return a (possibly const)
//    pointer to the beginning of the line and expose a pointer to the end of
//    the line.  This simultaneously saves memory and reduces overhead.
//    a. Since this reuses a single buffer, the string-view is invalidated when
//       the next line is read.
//    b. When the last line in the file is not terminated by '\n', this text
//       reader automatically appends '\n'.  Thus, while C library functions
//       that assume null-termination can't be used here (unless you're using a
//       wrapper function that replaces the terminating '\n' with '\0' after
//       the line-read function call, which is a totally valid thing to do),
//       plink2_string.h functions which either iterate to any end-of-line
//       character (ASCII code < 32 and unequal to 9=tab) or explicitly assume
//       '\n'-termination can be used, and these use essentially the same
//       optimizations as modern memchr implementations.
//    (A C++17 interface that returns std::string_view objects was considered,
//    but then rejected since std::string_view's design makes it much better
//    suited to be a function input-parameter type than a return type.  It is
//    easy enough to efficiently construct a string_view using the current
//    interface when it's time to call a function that accepts one.  The rest
//    of the time, there's no meaningful advantage over plain C pointers.)
// 4. can be used with either a single fixed-size memory buffer (this plays
//    well with plink2's memory allocation strategy), or dynamic resizing with
//    malloc()/realloc() calls.
//
// Two other readers are provided:
// - A decompress-ahead token reader.  This also shards the tokens, for the
//   common use case where the tokens don't need to be parsed in order (e.g.
//   --extract/--exclude).
// - A simpler single-threaded (no decompress-ahead) reader.

#ifdef STATIC_ZLIB
#  include "../zlib-1.2.11/zlib.h"
#else
#  include <zlib.h>
#  if !defined(ZLIB_VERNUM) || (ZLIB_VERNUM < 0x1240)
#    error "plink2_getline requires zlib 1.2.4 or later."
#  endif
#endif

#include "plink2_bgzf.h"
#include "plink2_zstfile.h"

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr GetFileType(const char* fname, FileCompressionType* ftype_ptr);

typedef struct TextRfileBaseStruct {
  // This *is* copied wholesale by the TextRstream move-constructor.

  // Positioned first so the compiler doesn't need to add an offset to access
  // it.
  char* consume_iter;

  char* consume_stop;  // should always point after the last loaded \n

  const char* errmsg;
  PglErr reterr;

  FileCompressionType file_type;
  FILE* ff;  // could use e.g. htslib for some network support later
  uint32_t dst_owned_by_consumer;
  uint32_t enforced_max_line_blen;
  // Note that when dst_owned_by_consumer isn't true, the reader thread can
  // alter the values below.  However, this can only happen when the consumer
  // thread is blocked waiting for the next line, or before it's attempted to
  // read the first line.  So there's no advantage to placing these in a
  // different cacheline from consume_iter (which is constantly altered by the
  // consumer).
  char* dst;
  uint32_t dst_len;
  uint32_t dst_capacity;
} TextRfileBase;

typedef struct GzRawDecompressStreamStruct {
  // Copied by TextRstream move-constructor.
  unsigned char* in;
  z_stream ds;
  uint32_t ds_initialized;
} GzRawDecompressStream;

// BgzfRawDecompressStream declared in plink2_bgzf.h.

typedef struct ZstRawDecompressStreamStruct {
  // Copied by TextRstream move-constructor.
  ZSTD_DStream* ds;
  ZSTD_inBuffer ib;
} ZstRawDecompressStream;

typedef union {
  GzRawDecompressStream gz;
  // Even in the single-threaded case, it's worth distinguishing bgzf from
  // generic .gz, since we can use libdeflate 100% of the time.
  BgzfRawDecompressStream bgzf;
  ZstRawDecompressStream zst;
} RawDecompressStream;

typedef struct textRFILEStruct {
  TextRfileBase base;
  RawDecompressStream rds;
} textRFILE;

void PreinitTextRfile(textRFILE* trfp);

// kDecompressChunkSize = 1 MiB currently declared in plink2_bgzf.h, may move
// somewhere closer to the base later.
CONSTI32(kTextRstreamBlenFast, 11 * kDecompressChunkSizeX);
CONSTI32(kTokenRstreamBlen, 11 * kDecompressChunkSizeX);
CONSTI32(kMaxTokenBlenX, 8 * kDecompressChunkSizeX);
static_assert(kMaxTokenBlenX >= kDecompressChunkSizeX, "kMaxTokenBlen too small.");
static_assert(kMaxTokenBlenX <= kTokenRstreamBlen, "kMaxTokenBlen can't be larger than kTokenRstreamBlen.");

// * Can return nomem, open-fail, or read-fail.
// * If dst == nullptr, this mallocs a buffer of size 2 * kDecompressChunkSize,
//   and it'll be realloced as necessary and freed by CleanupTextRfile().
//   The original value of dst_capacity doesn't matter in this case.
//   Otherwise, the buffer is owned by the caller, assumed to have size >=
//   dst_capacity, and never grown.
// * enforced_max_line_blen must be >= dst_capacity - kDecompressChunkSize.
//   It's the point at which long-line errors instead of out-of-memory errors
//   are reported.  It isn't permitted to be less than 1 MiB.
PglErr TextRfileOpenEx(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, char* dst, textRFILE* trfp);

HEADER_INLINE PglErr TextRfileOpen(const char* fname, textRFILE* trfp) {
  return TextRfileOpenEx(fname, kMaxLongLine, 0, nullptr, trfp);
}

extern const char kShortErrLongLine[];

PglErr TextRfileAdvance(textRFILE* trfp);

HEADER_INLINE PglErr TextRfileNextLine(textRFILE* trfp, char** line_startp) {
  if (trfp->base.consume_iter == trfp->base.consume_stop) {
    PglErr reterr = TextRfileAdvance(trfp);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
  }
  *line_startp = trfp->base.consume_iter;
  trfp->base.consume_iter = AdvPastDelim(trfp->base.consume_iter, '\n');
  return kPglRetSuccess;
}

HEADER_INLINE char* TextRfileLineEnd(textRFILE* trfp) {
  return trfp->base.consume_iter;
}

void TextRfileRewind(textRFILE* trfp);

HEADER_INLINE int32_t TextRfileIsOpen(const textRFILE* trfp) {
  return (trfp->base.ff != nullptr);
}

HEADER_INLINE int32_t TextRfileEof(const textRFILE* trfp) {
  return (trfp->base.reterr == kPglRetEof);
}

HEADER_INLINE const char* TextRfileError(const textRFILE* trfp) {
  return trfp->base.errmsg;
}

HEADER_INLINE PglErr TextRfileErrcode(const textRFILE* trfp) {
  if (trfp->base.reterr == kPglRetEof) {
    return kPglRetSuccess;
  }
  return trfp->base.reterr;
}

// Ok to pass reterrp == nullptr.
// Returns nonzero iff file-close fails, and either reterrp == nullptr or
// *reterrp == kPglRetSuccess; this is intended to be followed by logging of
// strerror(errno).  In the latter case, *reterrp is set to kPglRetReadFail.
BoolErr CleanupTextRfile(textRFILE* trfp, PglErr* reterrp);


// consumer -> reader message
// could add a "close current file and open another one" case
ENUM_U31_DEF_START()
  kTrsInterruptNone,
  kTrsInterruptRetarget,
  kTrsInterruptShutdown
ENUM_U31_DEF_END(TrsInterrupt);

typedef struct TextRstreamSyncStruct {
  // Mutex shared state, and everything guarded by the mutex.  Allocated to
  // different cacheline(s) than consume_stop.
#ifdef _WIN32
  CRITICAL_SECTION critical_section;
  HANDLE reader_progress_event;
  HANDLE consumer_progress_event;
#else
  pthread_mutex_t sync_mutex;
  pthread_cond_t reader_progress_condvar;
  pthread_cond_t consumer_progress_condvar;
  // bugfix (7 Mar 2018): need to avoid waiting on consumer_progress_condvar if
  // this is set.  (could also check an appropriate predicate)
  uint32_t consumer_progress_state;

  uint32_t sync_init_state;
#endif

  pthread_t read_thread;

  char* consume_tail;
  char* cur_circular_end;
  char* available_end;

  // Separate from the TextRfileBase instances of these values, since we don't
  // want to force the user to worry about these values changing at any moment.
  // Instead, the TextRfileBase instances are only updated during TextAdvance()
  // calls and the like.
  const char* errmsg;
  PglErr reterr;  // note that this is set to kPglRetEof once we reach eof

  uint32_t dst_reallocated;
  TrsInterrupt interrupt;
  const char* new_fname;
} TextRstreamSync;

typedef union {
  GzRawDecompressStream gz;
  BgzfRawMtDecompressStream bgzf;
  ZstRawDecompressStream zst;
} RawMtDecompressStream;

typedef struct TextRstreamStruct {
  TextRfileBase base;
  RawMtDecompressStream rds;
  uint32_t decompress_thread_ct;
  TextRstreamSync* syncp;
} TextRstream;

void PreinitTextRstream(TextRstream* trsp);

// * Can return nomem, open-fail, read-fail, or thread-create-fail.
// * Exactly one of fname and trfp must be nullptr.  If trfp is null, fname is
//   opened.  Otherwise, the returned stream is "move-constructed" from trfp.
//   When not move-constructing, enforced_max_line_blen, dst_capacity, and dst
//   are interpreted the same way as TextRfileOpenEx().
//   When move-constructing, enforced_max_line_blen and dst_capacity may be
//   smaller than what the textRFILE was opened with.
PglErr TextRstreamOpenEx(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, uint32_t decompress_thread_ct, textRFILE* trfp, char* dst, TextRstream* trsp);

HEADER_INLINE PglErr TextRstreamOpen(const char* fname, TextRstream* trsp) {
  return TextRstreamOpenEx(fname, kMaxLongLine, 0, NumCpu(nullptr), nullptr, nullptr, trsp);
}

HEADER_INLINE char* TextLineEnd(TextRstream* trsp) {
  return trsp->base.consume_iter;
}

HEADER_INLINE int32_t TextRstreamIsOpen(const TextRstream* trsp) {
  return (trsp->base.ff != nullptr);
}

HEADER_INLINE int32_t TextEof(const TextRstream* trsp) {
  return (trsp->base.reterr == kPglRetEof);
}

uint32_t TextDecompressThreadCt(const TextRstream* trsp);

PglErr TextAdvance(TextRstream* trsp);

HEADER_INLINE PglErr TextNextLine(TextRstream* trsp, char** line_startp) {
  if (trsp->base.consume_iter == trsp->base.consume_stop) {
    PglErr reterr = TextAdvance(trsp);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
  }
  *line_startp = trsp->base.consume_iter;
  trsp->base.consume_iter = AdvPastDelim(trsp->base.consume_iter, '\n');
  return kPglRetSuccess;
}

HEADER_INLINE PglErr TextNextLineK(TextRstream* trsp, const char** line_startp) {
  return TextNextLine(trsp, K_CAST(char**, line_startp));
}

HEADER_INLINE PglErr TextNextLineLstrip(TextRstream* trsp, char** line_startp) {
  if (trsp->base.consume_iter == trsp->base.consume_stop) {
    PglErr reterr = TextAdvance(trsp);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
  }
  *line_startp = FirstNonTspace(trsp->base.consume_iter);
  trsp->base.consume_iter = AdvPastDelim(*line_startp, '\n');
  return kPglRetSuccess;
}

HEADER_INLINE PglErr TextNextLineLstripK(TextRstream* trsp, const char** line_startp) {
  return TextNextLineLstrip(trsp, K_CAST(char**, line_startp));
}

// these do not update line_idx on eof.
PglErr TextNextNonemptyLineLstrip(TextRstream* trsp, uintptr_t* line_idx_ptr, char** line_startp);

HEADER_INLINE PglErr TextNextNonemptyLineLstripK(TextRstream* trsp, uintptr_t* line_idx_ptr, const char** line_startp) {
  return TextNextNonemptyLineLstrip(trsp, line_idx_ptr, K_CAST(char**, line_startp));
}

PglErr TextSkipNz(uintptr_t skip_ct, TextRstream* trsp);

HEADER_INLINE PglErr TextSkip(uintptr_t skip_ct, TextRstream* trsp) {
  if (skip_ct == 0) {
    return kPglRetSuccess;
  }
  return TextSkipNz(skip_ct, trsp);
}


// 'Unsafe' functions require line_iter to already point to the start of the
// line, and don't update trsp->base.consume_iter; they primarily wrap the
// TextAdvance() call.
HEADER_INLINE PglErr TextNextLineUnsafe(TextRstream* trsp, char** line_iterp) {
  if (*line_iterp != trsp->base.consume_stop) {
    return kPglRetSuccess;
  }
  trsp->base.consume_iter = *line_iterp;
  PglErr reterr = TextAdvance(trsp);
  // not unlikely() due to eof
  if (reterr) {
    return reterr;
  }
  *line_iterp = trsp->base.consume_iter;
  return kPglRetSuccess;
}


HEADER_INLINE PglErr TextNextLineLstripUnsafe(TextRstream* trsp, char** line_iterp) {
  char* line_iter = *line_iterp;
  if (line_iter == trsp->base.consume_stop) {
    trsp->base.consume_iter = line_iter;
    PglErr reterr = TextAdvance(trsp);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
    line_iter = trsp->base.consume_iter;
  }
  *line_iterp = FirstNonTspace(line_iter);
  return kPglRetSuccess;
}

PglErr TextNextNonemptyLineLstripUnsafe(TextRstream* trsp, uintptr_t* line_idx_ptr, char** line_iterp);


HEADER_INLINE uint32_t TextRstreamIsMt(const TextRstream* trsp) {
  // Only bgzf decoder is multithreaded for now.
  return (trsp->base.file_type == kFileBgzf);
}

PglErr TextRetarget(const char* new_fname, TextRstream* trsp);

HEADER_INLINE PglErr TextRewind(TextRstream* trsp) {
  return TextRetarget(nullptr, trsp);
}

HEADER_INLINE const char* TextRstreamError(const TextRstream* trsp) {
  return trsp->base.errmsg;
}

HEADER_INLINE PglErr TextRstreamErrcode(const TextRstream* trsp) {
  if (trsp->base.reterr == kPglRetEof) {
    return kPglRetSuccess;
  }
  return trsp->base.reterr;
}

// Ok to pass reterrp == nullptr.
// Returns nonzero iff file-close fails, and either reterrp == nullptr or
// *reterrp == kPglRetSuccess.  In the latter case, *reterrp is set to
// kPglRetReadFail.  (Note that this does *not* retrieve the existing
// trsp->reterr value; caller is responsible for checking TextRstreamErrcode()
// first when they care.)
BoolErr CleanupTextRstream(TextRstream* trsp, PglErr* reterrp);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_GETLINE_H__
