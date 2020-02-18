#ifndef __PLINK2_TEXT_H__
#define __PLINK2_TEXT_H__

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
#  include "../../zlib-1.2.11/zlib.h"
#else
#  include <zlib.h>
#  if !defined(ZLIB_VERNUM) || (ZLIB_VERNUM < 0x1240)
#    error "plink2_text requires zlib 1.2.4 or later."
#  endif
#endif

#include "plink2_bgzf.h"
#include "plink2_zstfile.h"

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr GetFileType(const char* fname, FileCompressionType* ftype_ptr);

typedef struct TextFileBaseStruct {
  // Not noncopyable, since this is copied wholesale by the TextStream
  // move-constructor.

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
} TextFileBase;

typedef struct GzRawDecompressStreamStruct {
  // Copied by TextStream move-constructor.
  unsigned char* in;
  z_stream ds;
  uint32_t ds_initialized;
} GzRawDecompressStream;

// BgzfRawDecompressStream declared in plink2_bgzf.h.

typedef struct ZstRawDecompressStreamStruct {
  // Copied by TextStream move-constructor.
  ZSTD_DStream* ds;
  ZSTD_inBuffer ib;
} ZstRawDecompressStream;

typedef union {
  GzRawDecompressStream gz;
  // Even in the single-threaded case, it's worth distinguishing bgzf from
  // generic .gz, since with bgzf we can use libdeflate.
  BgzfRawDecompressStream bgzf;
  ZstRawDecompressStream zst;
} RawDecompressStream;

typedef struct textFILEMainStruct {
  TextFileBase base;
  RawDecompressStream rds;
} textFILEMain;

typedef struct textFILEStruct {
#ifdef __cplusplus
  textFILEMain& GET_PRIVATE_m() { return m; }
  textFILEMain const& GET_PRIVATE_m() const { return m; }
 private:
#endif
  textFILEMain m;
} textFILE;

void PreinitTextFile(textFILE* txf_ptr);

// kDecompressChunkSize = 1 MiB currently declared in plink2_bgzf.h, may move
// somewhere closer to the base later.
CONSTI32(kTextStreamBlenFast, 11 * kDecompressChunkSize);
CONSTI32(kTokenStreamBlen, 12 * kDecompressChunkSize);
CONSTI32(kMaxTokenBlen, 8 * kDecompressChunkSize);
static_assert(kMaxTokenBlen >= kDecompressChunkSize, "kMaxTokenBlen too small.");
static_assert(kMaxTokenBlen + kDecompressChunkSize <= kTokenStreamBlen, "kMaxTokenBlen + kDecompressChunkSize can't be larger than kTokenStreamBlen.");

// max_line_blen and enforced_max_line_blen lower bound.
CONSTI32(kDecompressMinBlen, kDecompressChunkSize);

CONSTI32(kDecompressMinCapacity, kDecompressMinBlen + kDecompressChunkSize);

// * Can return nomem, open-fail, or read-fail.
// * If dst == nullptr, this mallocs a buffer of size 2 * kDecompressChunkSize,
//   and it'll be realloced as necessary and freed by CleanupTextFile().
//   The original value of dst_capacity doesn't matter in this case.
//   Otherwise, the buffer is owned by the caller, assumed to have size >=
//   dst_capacity, and never grown.
// * enforced_max_line_blen must be >= dst_capacity - kDecompressChunkSize.
//   It's the point at which long-line errors instead of out-of-memory errors
//   are reported.  It isn't permitted to be less than 1 MiB.
PglErr TextFileOpenEx(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, char* dst, textFILE* txf_ptr);

HEADER_INLINE PglErr TextFileOpen(const char* fname, textFILE* txf_ptr) {
  return TextFileOpenEx(fname, kMaxLongLine, 0, nullptr, txf_ptr);
}

extern const char kShortErrLongLine[];
extern const char kShortErrInteriorEmptyLine[];

PglErr TextFileAdvance(textFILE* txf_ptr);

HEADER_INLINE PglErr TextFileNextLine(textFILE* txf_ptr, char** line_startp) {
  TextFileBase* basep = &GET_PRIVATE(*txf_ptr, m).base;
  if (basep->consume_iter == basep->consume_stop) {
    PglErr reterr = TextFileAdvance(txf_ptr);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
  }
  *line_startp = basep->consume_iter;
  basep->consume_iter = AdvPastDelim(*line_startp, '\n');
  return kPglRetSuccess;
}

HEADER_INLINE char* TextFileLineEnd(textFILE* txf_ptr) {
  return GET_PRIVATE(*txf_ptr, m).base.consume_iter;
}

HEADER_INLINE PglErr TextFileNextLineLstrip(textFILE* txf_ptr, char** line_startp) {
  TextFileBase* basep = &GET_PRIVATE(*txf_ptr, m).base;
  if (basep->consume_iter == basep->consume_stop) {
    PglErr reterr = TextFileAdvance(txf_ptr);
    if (reterr) {
      return reterr;
    }
  }
  *line_startp = FirstNonTspace(basep->consume_iter);
  basep->consume_iter = AdvPastDelim(*line_startp, '\n');
  return kPglRetSuccess;
}

PglErr TextFileOnlyEmptyLinesLeft(textFILE* txf_ptr);

// API fix (30-31 Oct 2019): While we want to tolerate *trailing* empty lines
// (as well as lack of a final newline), since those arise from manual
// text-editing every once in a while, there's usually no good reason to
// tolerate *interior* empty lines like we were previously doing.  So the
// generic empty-line-skipping interface has been replaced with one that
// reports an "Unexpected empty line" error with a kPglRetMalformedInput error
// code unless the empty line(s) are at EOF.
HEADER_INLINE char* TextFileGet(textFILE* txf_ptr) {
  TextFileBase* basep = &GET_PRIVATE(*txf_ptr, m).base;
  if (basep->consume_iter == basep->consume_stop) {
    // not unlikely() due to eof
    if (TextFileAdvance(txf_ptr)) {
      return nullptr;
    }
  }
  char* line_start = FirstNonTspace(basep->consume_iter);
  basep->consume_iter = AdvPastDelim(line_start, '\n');
  if (!IsEolnKns(*line_start)) {
    return line_start;
  }
  TextFileOnlyEmptyLinesLeft(txf_ptr);
  return nullptr;
}

void TextFileRewind(textFILE* txf_ptr);

HEADER_INLINE int32_t TextFileIsOpen(const textFILE* txf_ptr) {
  return (GET_PRIVATE(*txf_ptr, m).base.ff != nullptr);
}

HEADER_INLINE int32_t TextFileEof(const textFILE* txf_ptr) {
  return (GET_PRIVATE(*txf_ptr, m).base.reterr == kPglRetEof);
}

HEADER_INLINE const char* TextFileError(const textFILE* txf_ptr) {
  return GET_PRIVATE(*txf_ptr, m).base.errmsg;
}

HEADER_INLINE PglErr TextFileErrcode(const textFILE* txf_ptr) {
  const PglErr reterr = GET_PRIVATE(*txf_ptr, m).base.reterr;
  if (reterr == kPglRetEof) {
    return kPglRetSuccess;
  }
  return reterr;
}

// Note that this does not assign to *reterrp at all in the usual eof case.
HEADER_INLINE BoolErr TextFileErrcode2(const textFILE* txf_ptr, PglErr* reterrp) {
  const PglErr reterr = GET_PRIVATE(*txf_ptr, m).base.reterr;
  if (reterr == kPglRetEof) {
    return 0;
  }
  *reterrp = reterr;
  return 1;
}

// Does not convert kPglRetEof -> kPglRetSuccess.
HEADER_INLINE PglErr TextFileRawErrcode(const textFILE* txf_ptr) {
  return GET_PRIVATE(*txf_ptr, m).base.reterr;
}

// Relevant when about to move-construct a TextStream; see plink2_glm.
HEADER_INLINE uint32_t TextFileLinebufLen(const textFILE* txf_ptr) {
  return GET_PRIVATE(*txf_ptr, m).base.dst_len;
}

// Ok to pass reterrp == nullptr.
// Returns nonzero iff file-close fails, and either reterrp == nullptr or
// *reterrp == kPglRetSuccess; this is intended to be followed by logging of
// strerror(errno).  In the latter case, *reterrp is set to kPglRetReadFail.
BoolErr CleanupTextFile(textFILE* txf_ptr, PglErr* reterrp);


// consumer -> reader message
// could add a "close current file and open another one" case
ENUM_U31_DEF_START()
  kTxsInterruptNone,
  kTxsInterruptRetarget,
  kTxsInterruptShutdown
ENUM_U31_DEF_END(TxsInterrupt);

typedef struct TextStreamSyncStruct {
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

  // Separate from the TextFileBase instances of these values, since we don't
  // want to force the user to worry about these values changing at any moment.
  // Instead, the TextFileBase instances are only updated during TextAdvance()
  // calls and the like.
  const char* errmsg;
  PglErr reterr;  // note that this is set to kPglRetEof once we reach eof

  uint32_t dst_reallocated;
  TxsInterrupt interrupt;
  const char* new_fname;
} TextStreamSync;

typedef union {
  GzRawDecompressStream gz;
  BgzfRawMtDecompressStream bgzf;
  ZstRawDecompressStream zst;
} RawMtDecompressStream;

typedef struct TextStreamMainStruct {
  TextFileBase base;
  RawMtDecompressStream rds;
  uint32_t decompress_thread_ct;
  TextStreamSync* syncp;
} TextStreamMain;

typedef struct TextStreamStruct {
#ifdef __cplusplus
  TextStreamMain& GET_PRIVATE_m() { return m; }
  TextStreamMain const& GET_PRIVATE_m() const { return m; }
 private:
#endif
  TextStreamMain m;
} TextStream;

void PreinitTextStream(TextStream* txs_ptr);

// * Can return nomem, open-fail, read-fail, or thread-create-fail.
// * Exactly one of fname and txf_ptr must be nullptr.  If txf_ptr is null,
//   fname is opened.  Otherwise, the returned stream is "move-constructed"
//   from txf_ptr.
//   When not move-constructing, enforced_max_line_blen, dst_capacity, and dst
//   are interpreted the same way as TextFileOpenEx().
//   When move-constructing, enforced_max_line_blen and dst_capacity may be
//   smaller than what the textFILE was opened with.
PglErr TextStreamOpenEx(const char* fname, uint32_t enforced_max_line_blen, uint32_t dst_capacity, uint32_t decompress_thread_ct, textFILE* txf_ptr, char* dst, TextStream* txs_ptr);

HEADER_INLINE PglErr TextStreamOpen(const char* fname, TextStream* txs_ptr) {
  return TextStreamOpenEx(fname, kMaxLongLine, 0, NumCpu(nullptr), nullptr, nullptr, txs_ptr);
}

// We drop 'Stream' from the function names outside of open/close, to emphasize
// that this is the default choice.
// (Originally this was named 'TextRstream', but then I realized that the write
// case doesn't really care about whether the input is text or binary.)
HEADER_INLINE char* TextLineEnd(TextStream* txs_ptr) {
  return GET_PRIVATE(*txs_ptr, m).base.consume_iter;
}

HEADER_INLINE int32_t TextIsOpen(const TextStream* txs_ptr) {
  return (GET_PRIVATE(*txs_ptr, m).base.ff != nullptr);
}

HEADER_INLINE int32_t TextEof(const TextStream* txs_ptr) {
  return (GET_PRIVATE(*txs_ptr, m).base.reterr == kPglRetEof);
}

uint32_t TextDecompressThreadCt(const TextStream* txs_ptr);

PglErr TextAdvance(TextStream* txs_ptr);

HEADER_INLINE PglErr TextNextLine(TextStream* txs_ptr, char** line_startp) {
  TextFileBase* basep = &GET_PRIVATE(*txs_ptr, m).base;
  if (basep->consume_iter == basep->consume_stop) {
    PglErr reterr = TextAdvance(txs_ptr);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
  }
  *line_startp = basep->consume_iter;
  basep->consume_iter = AdvPastDelim(basep->consume_iter, '\n');
  return kPglRetSuccess;
}

HEADER_INLINE PglErr TextNextLineK(TextStream* txs_ptr, const char** line_startp) {
  return TextNextLine(txs_ptr, K_CAST(char**, line_startp));
}

// plink2 functions strip leading whitespace by default (unless there's a clear
// reason not to, e.g. it isn't allowed in VCF files), since plink 1.x inserts
// it so often.
HEADER_INLINE PglErr TextNextLineLstrip(TextStream* txs_ptr, char** line_startp) {
  TextFileBase* basep = &GET_PRIVATE(*txs_ptr, m).base;
  if (basep->consume_iter == basep->consume_stop) {
    PglErr reterr = TextAdvance(txs_ptr);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
  }
  *line_startp = FirstNonTspace(basep->consume_iter);
  basep->consume_iter = AdvPastDelim(*line_startp, '\n');
  return kPglRetSuccess;
}

HEADER_INLINE PglErr TextNextLineLstripK(TextStream* txs_ptr, const char** line_startp) {
  return TextNextLineLstrip(txs_ptr, K_CAST(char**, line_startp));
}

// returns kPglRetEof if true, otherwise kPglRetMalformedInput and sets
// errmsg.
PglErr TextOnlyEmptyLinesLeft(TextStream* txs_ptr);

// This was previously named TextNextLineLstripNoempty, and returned a PglErr,
// but some compilers were (i) not inlining it and (ii) spamming spurious
// uninitialized-variable warnings as a side effect.  And it's plink2's most
// common use case.
// So we'll use a different interface here (and a shorter name).  User is
// expected to call TextStream[Raw]Errcode when nullptr is returned.
HEADER_INLINE char* TextGet(TextStream* txs_ptr) {
  TextFileBase* basep = &GET_PRIVATE(*txs_ptr, m).base;
  if (basep->consume_iter == basep->consume_stop) {
    // not unlikely() due to eof
    if (TextAdvance(txs_ptr)) {
      return nullptr;
    }
  }
  char* line_start = FirstNonTspace(basep->consume_iter);
  basep->consume_iter = AdvPastDelim(line_start, '\n');
  if (!IsEolnKns(*line_start)) {
    return line_start;
  }
  TextOnlyEmptyLinesLeft(txs_ptr);
  return nullptr;
}

PglErr TextSkipNz(uintptr_t skip_ct, TextStream* txs_ptr);

HEADER_INLINE PglErr TextSkip(uintptr_t skip_ct, TextStream* txs_ptr) {
  if (skip_ct == 0) {
    return kPglRetSuccess;
  }
  return TextSkipNz(skip_ct, txs_ptr);
}


// 'Unsafe' functions require line_iter to already point to the start of the
// line, and don't update txsp->base.consume_iter; they primarily wrap the
// TextAdvance() call.
HEADER_INLINE PglErr TextNextLineUnsafe(TextStream* txs_ptr, char** line_iterp) {
  TextFileBase* basep = &GET_PRIVATE(*txs_ptr, m).base;
  if (*line_iterp != basep->consume_stop) {
    return kPglRetSuccess;
  }
  basep->consume_iter = *line_iterp;
  PglErr reterr = TextAdvance(txs_ptr);
  // not unlikely() due to eof
  if (reterr) {
    return reterr;
  }
  *line_iterp = basep->consume_iter;
  return kPglRetSuccess;
}


HEADER_INLINE PglErr TextNextLineLstripUnsafe(TextStream* txs_ptr, char** line_iterp) {
  char* line_iter = *line_iterp;
  TextFileBase* basep = &GET_PRIVATE(*txs_ptr, m).base;
  if (line_iter == basep->consume_stop) {
    basep->consume_iter = line_iter;
    PglErr reterr = TextAdvance(txs_ptr);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
    line_iter = basep->consume_iter;
  }
  *line_iterp = FirstNonTspace(line_iter);
  return kPglRetSuccess;
}

HEADER_INLINE PglErr TextNextLineLstripUnsafeK(TextStream* txs_ptr, const char** line_iterp) {
  return TextNextLineLstripUnsafe(txs_ptr, K_CAST(char**, line_iterp));
}

HEADER_INLINE PglErr TextGetUnsafe(TextStream* txs_ptr, char** line_iterp) {
  char* line_iter = *line_iterp;
  TextFileBase* basep = &GET_PRIVATE(*txs_ptr, m).base;
  if (line_iter == basep->consume_stop) {
    basep->consume_iter = line_iter;
    PglErr reterr = TextAdvance(txs_ptr);
    // not unlikely() due to eof
    if (reterr) {
      return reterr;
    }
    line_iter = basep->consume_iter;
  }
  line_iter = FirstNonTspace(line_iter);
  if (!IsEolnKns(*line_iter)) {
    *line_iterp = line_iter;
    return kPglRetSuccess;
  }
  return TextOnlyEmptyLinesLeft(txs_ptr);
}

/*
HEADER_INLINE PglErr TextGetUnsafeK(TextStream* txs_ptr, const char** line_iterp) {
  return TextGetUnsafe(txs_ptr, K_CAST(char**, line_iterp));
}
*/

// Returns *zero* when it's time to stop iterating.  Designed to be the middle
// argument in a for (; ; ) loop.
HEADER_INLINE uint32_t TextGetUnsafe2(TextStream* txs_ptr, char** line_iterp) {
  char* line_iter = *line_iterp;
  TextFileBase* basep = &GET_PRIVATE(*txs_ptr, m).base;
  if (line_iter == basep->consume_stop) {
    basep->consume_iter = line_iter;
    // not unlikely() due to eof
    if (TextAdvance(txs_ptr)) {
      return 0;
    }
    line_iter = basep->consume_iter;
  }
  line_iter = FirstNonTspace(line_iter);
  if (!IsEolnKns(*line_iter)) {
    *line_iterp = line_iter;
    return 1;
  }
  TextOnlyEmptyLinesLeft(txs_ptr);
  return 0;
}

HEADER_INLINE uint32_t TextGetUnsafe2K(TextStream* txs_ptr, const char** line_iterp) {
  return TextGetUnsafe2(txs_ptr, K_CAST(char**, line_iterp));
}

HEADER_INLINE void TextSetPos(char* new_consume_iter, TextStream* txs_ptr) {
  GET_PRIVATE(*txs_ptr, m).base.consume_iter = new_consume_iter;
}


HEADER_INLINE uint32_t TextIsMt(const TextStream* txs_ptr) {
  // Only bgzf decoder is multithreaded for now.
  return (GET_PRIVATE(*txs_ptr, m).base.file_type == kFileBgzf);
}

PglErr TextRetarget(const char* new_fname, TextStream* txs_ptr);

HEADER_INLINE PglErr TextRewind(TextStream* txs_ptr) {
  return TextRetarget(nullptr, txs_ptr);
}

HEADER_INLINE const char* TextStreamError(const TextStream* txs_ptr) {
  return GET_PRIVATE(*txs_ptr, m).base.errmsg;
}

HEADER_INLINE PglErr TextStreamErrcode(const TextStream* txs_ptr) {
  const PglErr reterr = GET_PRIVATE(*txs_ptr, m).base.reterr;
  if (reterr == kPglRetEof) {
    return kPglRetSuccess;
  }
  return reterr;
}

// Note that this does not assign to *reterrp at all in the usual eof case.
HEADER_INLINE BoolErr TextStreamErrcode2(const TextStream* txs_ptr, PglErr* reterrp) {
  const PglErr reterr = GET_PRIVATE(*txs_ptr, m).base.reterr;
  if (reterr == kPglRetEof) {
    return 0;
  }
  *reterrp = reterr;
  return 1;
}

// Does not convert kPglRetEof -> kPglRetSuccess.
HEADER_INLINE PglErr TextStreamRawErrcode(const TextStream* txs_ptr) {
  return GET_PRIVATE(*txs_ptr, m).base.reterr;
}

// Ok to pass reterrp == nullptr.
// Returns nonzero iff file-close fails, and either reterrp == nullptr or
// *reterrp == kPglRetSuccess.  In the latter case, *reterrp is set to
// kPglRetReadFail.  (Note that this does *not* retrieve the existing
// txsp->reterr value; caller is responsible for checking TextStreamErrcode()
// first when they care.)
BoolErr CleanupTextStream(TextStream* txs_ptr, PglErr* reterrp);


// Low-level token-batch-reading interface, using the extra TextStream mode
// which cares about token rather than line endings.
typedef struct TokenStreamStruct {
  NONCOPYABLE(TokenStreamStruct);
  TextStream txs;
} TokenStream;

HEADER_INLINE void PreinitTokenStream(TokenStream* tksp) {
  PreinitTextStream(&tksp->txs);
}

// Note that shard_boundaries must have length (piece_ct + 1).
PglErr TksNext(TokenStream* tksp, uint32_t shard_ct, char** shard_boundaries);

HEADER_INLINE PglErr TokenStreamRetarget(const char* new_fname, TokenStream* tksp) {
  return TextRetarget(new_fname, &(tksp->txs));
}

HEADER_INLINE PglErr TokenRewind(TokenStream* tksp) {
  return TextRetarget(nullptr, &(tksp->txs));
}

HEADER_INLINE const char* TokenStreamError(const TokenStream* tksp) {
  return GET_PRIVATE(tksp->txs, m).base.errmsg;
}

HEADER_INLINE PglErr TokenStreamErrcode(const TokenStream* tksp) {
  const PglErr reterr = GET_PRIVATE(tksp->txs, m).base.reterr;
  if (reterr == kPglRetEof) {
    return kPglRetSuccess;
  }
  return reterr;
}

HEADER_INLINE BoolErr CleanupTokenStream(TokenStream* tksp, PglErr* reterrp) {
  return CleanupTextStream(&(tksp->txs), reterrp);
}

// Could create a slightly simpler interface for the one-token-at-a-time case,
// but I won't bother for now since it's kind of good for the relative cost of
// parallelizing token processing to be low.


#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_TEXT_H__
