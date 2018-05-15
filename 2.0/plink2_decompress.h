#ifndef __PLINK2_DECOMPRESS_H__
#define __PLINK2_DECOMPRESS_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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
// heavyweight dependence on zstd/zlibWrapper.
#include "plink2_cmdline.h"

// documentation on ZWRAP_USE_ZSTD is incorrect as of 11 Jan 2017, necessary to
// edit zstd_zlibwrapper.c or use compile flag.
#include "zstd/zlibWrapper/zstd_zlibwrapper.h"
#ifndef STATIC_ZLIB
#  if !defined(ZLIB_VERNUM) || ZLIB_VERNUM < 0x1240
#    error "plink2_decompress requires zlib 1.2.4 or later."
#  endif
#endif

#include "htslib/htslib/bgzf.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Also sets 128k read buffer.
PglErr gzopen_read_checked(const char* fname, gzFile* gzf_ptr);

// This sets loadbuf[loadbuf_size - 1] to ' ', just because.
// loadbuf_size assumed to be either exactly kMaxMediumLine (in which case
// any longer line is treated as pathological), or strictly larger (in which
// case we report an out-of-memory error when gzgets blows the buffer, unless
// loadbuf_size == kMaxLongLine, which is close to 2GiB).
PglErr GzopenAndSkipFirstLines(const char* fname, uint32_t lines_to_skip, uintptr_t loadbuf_size, char* loadbuf, gzFile* gzf_ptr);

// plink2_compress_stream interface should be used for writing .gz files.

HEADER_INLINE BoolErr gzclose_null(gzFile* gzf_ptr) {
  const int32_t ii = gzclose(*gzf_ptr);
  *gzf_ptr = nullptr;
  return (ii != Z_OK);
}

HEADER_INLINE void gzclose_cond(gzFile gz_infile) {
  if (gz_infile) {
    gzclose(gz_infile);
  }
}


// currently hardcoded to have maximum token length = kMaxMediumLine, buffer
// size = 2 * kMaxMediumLine * 2.
typedef struct GzTokenStreamStruct {
  NONCOPYABLE(GzTokenStreamStruct);
  gzFile gz_infile;
  char* buf_start;
  char* read_iter;
  char* buf_end;
} GzTokenStream;

void PreinitGzTokenStream(GzTokenStream* gtsp);

PglErr InitGzTokenStream(const char* fname, GzTokenStream* gtsp, char* buf_start);

// sets token_slen to 0xfffffffeU on read fail, 0xffffffffU on too-long token
// safe to null-terminate token between calls
char* AdvanceGzTokenStream(GzTokenStream* gtsp, uint32_t* token_slen_ptr);

// ok if already closed
BoolErr CloseGzTokenStream(GzTokenStream* gtsp);


PglErr IsBgzf(const char* fname, uint32_t* is_bgzf_ptr);

// While a separate non-compressing writer thread is practically useless (since
// the OS does its own scheduling, etc. of the writes anyway), and the OS is
// also usually good enough at sequential read-ahead, getting
// *decompress-ahead* functionality for 'free' is a pretty big deal.
//
// Possible todos: tabix index support (caller must declare the entire sequence
// of virtual offset ranges they want to read upfront), and optional
// multithreaded decompression of BGZF and seekable-Zstd files.  Note that
// these require us to stop using the Zstd zlib wrapper, so they're unlikely to
// be worth the time investment before the initial plink 2.0 beta release.

// (tested a few different values for this, 1 MiB appears to work well on the
// systems we care most about)
CONSTI32(kDecompressChunkSize, 1048576);
static_assert(!(kDecompressChunkSize % kCacheline), "kDecompressChunkSize must be a multiple of kCacheline.");

// consumer -> reader message
// could add a "close current file and open another one" case
ENUM_U31_DEF_START()
  kRlsInterruptNone,
  kRlsInterruptRetarget,
  kRlsInterruptShutdown
ENUM_U31_DEF_END(RlsInterrupt);

typedef struct ReadLineStreamSyncStruct {
  // Mutex shared state, and everything guarded by the mutex.  Allocated to
  // different cacheline(s) than consume_stop.
  NONCOPYABLE(ReadLineStreamSyncStruct);

#ifdef _WIN32
  CRITICAL_SECTION critical_section;
#else
  pthread_mutex_t sync_mutex;
  pthread_cond_t reader_progress_condvar;
  pthread_cond_t consumer_progress_condvar;
  // bugfix (7 Mar 2018): need to avoid waiting on consumer_progress_condvar if
  // this is set.  (could also check an appropriate predicate)
  uint32_t consumer_progress_state;
#endif

  char* consume_tail;
  char* cur_circular_end;
  char* available_end;
  PglErr reterr;  // note that this is set to kPglRetEof once we reach eof

  RlsInterrupt interrupt;
  const char* new_fname;
} ReadLineStreamSync;

// To minimize false (or true) sharing penalties, these values shouldn't change
// much; only the things they point to should be frequently changing.
typedef struct ReadLineStreamStruct {
  NONCOPYABLE(ReadLineStreamStruct);
  // Positioned first so the compiler doesn't need to add an offset to compare
  // to this.
  char* consume_stop;

#ifdef _WIN32
  // stored here since they're just pointers.
  HANDLE reader_progress_event;
  HANDLE consumer_progress_event;
#endif
  ReadLineStreamSync* syncp;

  gzFile gz_infile;
  BGZF* bgz_infile;
  uint32_t bgzf_decompress_thread_ct;  // this may change before rewind
  // This is aimed at (usually uncompressed) text files, so just use char*
  // instead of unsigned char*.
  char* buf;
  char* buf_end;

  pthread_t read_thread;
#ifndef _WIN32
  uint32_t sync_init_state;
#endif
  uint32_t enforced_max_line_blen;
} ReadLineStream;

void PreinitRLstream(ReadLineStream* rlsp);

CONSTI32(kRLstreamBlenLowerBound, 2 * kDecompressChunkSize);
static_assert(kRLstreamBlenLowerBound >= kMaxMediumLine, "max_line_blen lower limit too small.");
// Noticeable slowdown below this buffer size (11+1 = 12 MiB) on my dev
// machine.  (Additional buffer space beyond this point can still help a little
// bit.)
CONSTI32(kRLstreamBlenFast, 11 * kDecompressChunkSize);

// required_byte_ct can't be greater than kMaxLongLine.
// unstandardized_byte_ct cannot be bigstack_left() here, because the read
// stream requires a few additional allocations; use
// StandardizeLinebufSizemax() for that case instead.
HEADER_INLINE BoolErr StandardizeLinebufSize(uintptr_t unstandardized_byte_ct, uintptr_t required_byte_ct, uintptr_t* linebuf_sizep) {
  if (unstandardized_byte_ct >= kMaxLongLine) {
    *linebuf_sizep = kMaxLongLine;
    return 0;
  }
  if (unstandardized_byte_ct < RoundUpPow2(MAXV(kRLstreamBlenLowerBound, required_byte_ct), kCacheline)) {
    return 1;
  }
  *linebuf_sizep = RoundDownPow2(unstandardized_byte_ct, kCacheline);
  return 0;
}

// Must be synced with InitRLstreamRaw().
HEADER_INLINE BoolErr StandardizeLinebufSizemax(uintptr_t required_byte_ct, uintptr_t* linebuf_sizep) {
  const uintptr_t extra_alloc_byte_ct = RoundUpPow2(sizeof(ReadLineStreamSync) + 1, kCacheline) + kDecompressChunkSize;
  uintptr_t linebuf_size = bigstack_left();
  if (linebuf_size < extra_alloc_byte_ct) {
    return 1;
  }
  linebuf_size -= extra_alloc_byte_ct;
  if (linebuf_size >= kMaxLongLine) {
    *linebuf_sizep = kMaxLongLine;
    return 0;
  }
  if (linebuf_size < RoundUpPow2(MAXV(kRLstreamBlenLowerBound, required_byte_ct), kCacheline)) {
    return 1;
  }
  *linebuf_sizep = RoundDownPow2(linebuf_size, kCacheline);
  return 0;
}

// This supports two simple line iterators, one insane, and the other somewhat
// sane.  I'll discuss the insane one first.
// * consume_iter is effectively initialized to position -1, not 0, in the
//   file, and points to a '\n' byte.  This is necessary to make a simple loop
//   calling ReadNextLineFromRLstreamRaw() do the right thing.
// * On return from that function, consume_iter points to the beginning of a
//   line, that line is guaranteed to be '\n' terminated even if it's an
//   unterminated last line in the original file, and it is safe to mutate the
//   bytes on that line as long as you don't clobber the '\n' or insert another
//   one.  However, the line is NOT null-terminated.  You are expected to use
//   functions like strchrnul_n() and FirstPrespace() in place of the standard
//   C string library functions.
// * Since ReadNextLineFromRLstreamRaw starts with a rawmemchr operation, you
//   get the best performance by setting consume_iter to the last known
//   position in the previous line (pointing to the actual terminating '\n' is
//   fine), rather than leaving it in front.
// The sane interface gets rid of the roughest edges by returning a pointer to
// the end of the line (original position of '\n') as well, and automatically
// replacing the '\n' with a null character.  It also doesn't care if you
// insert a '\n' in the middle of the line or mutate the last character before
// continuing iteration, as long as you don't manually change line_last.

// ***** begin insane API function block ****

// if max_line_blen is a multiple of kCacheline and we're allocating on the
// bottom, total bigstack usage is max_line_blen + this value.
HEADER_CINLINE uintptr_t GetRLstreamExtraAlloc() {
  return RoundUpPow2(sizeof(ReadLineStreamSync) + 1, kCacheline) + kDecompressChunkSize;
}

// Does not print error message on fopen failure.
PglErr RlsOpenMaybeBgzf(const char* fname, uint32_t calc_thread_ct, ReadLineStream* rlsp);

// gz_infile assumed to already be opened.  Ok if e.g. header line has already
// been read from it.
// enforced_max_line_blen must be a multiple of kCacheline.
// The caller is responsible for ensuring max_line_blen >=
// kRLstreamBlenLowerBound, and preferably >= kRLstreamBlenFast, if at all
// possible.  The lower bound isn't hard--the stream will still work properly
// with a very small value--but there will be a lot of pointless
// synchronization in that case.
// name inspired by the Win32 API
PglErr InitRLstreamEx(uint32_t alloc_at_end, uint32_t enforced_max_line_blen, uint32_t max_line_blen, ReadLineStream* rlsp, char** consume_iterp);

HEADER_INLINE PglErr InitRLstreamRaw(const char* fname, uint32_t max_line_blen, ReadLineStream* rlsp, char** consume_iterp) {
  PglErr reterr = gzopen_read_checked(fname, &rlsp->gz_infile);
  if (reterr) {
    return reterr;
  }
  return InitRLstreamEx(0, kMaxLongLine, max_line_blen, rlsp, consume_iterp);
}

HEADER_INLINE PglErr InitRLstreamEndallocRaw(const char* fname, uint32_t max_line_blen, ReadLineStream* rlsp, char** consume_iterp) {
  PglErr reterr = gzopen_read_checked(fname, &rlsp->gz_infile);
  if (reterr) {
    return reterr;
  }
  return InitRLstreamEx(1, kMaxLongLine, max_line_blen, rlsp, consume_iterp);
}

HEADER_INLINE PglErr InitRLstreamFastsizeRaw(const char* fname, ReadLineStream* rlsp, char** consume_iterp) {
  PglErr reterr = gzopen_read_checked(fname, &rlsp->gz_infile);
  if (reterr) {
    return reterr;
  }
  return InitRLstreamEx(0, kRLstreamBlenFast, kRLstreamBlenFast, rlsp, consume_iterp);
}

HEADER_INLINE PglErr SizeAndInitRLstreamRaw(const char* fname, uintptr_t unstandardized_byte_ct, ReadLineStream* rlsp, char** consume_iterp) {
  // plink 1.9 immediately failed with an out-of-memory error if a "long line"
  // buffer would be smaller than kMaxMediumLine + 1 bytes, so may as well make
  // that the default lower bound.  (The precise value is currently irrelevant
  // since kRLstreamBlenLowerBound is larger and we take the maximum of the two
  // at compile time, but it's useful to distinguish "minimum acceptable
  // potentially-long-line buffer size" from "load/decompression block size
  // which generally has good performance".)
  uintptr_t linebuf_size;
  if (StandardizeLinebufSize(unstandardized_byte_ct, kMaxMediumLine + 1, &linebuf_size)) {
    return kPglRetNomem;
  }
  return InitRLstreamRaw(fname, linebuf_size, rlsp, consume_iterp);
}

// Easier to write correct code if consumer doesn't modify the buffer without
// very good reason.
HEADER_INLINE PglErr SizeAndInitRLstreamRawK(const char* fname, uintptr_t unstandardized_byte_ct, ReadLineStream* rlsp, const char** consume_iterp) {
  return SizeAndInitRLstreamRaw(fname, unstandardized_byte_ct, rlsp, K_CAST(char**, consume_iterp));
}

// When we want to enforce a tighter lower bound than kMaxMediumLine + 1.
HEADER_INLINE PglErr LboundedSizeAndInitRLstreamRaw(const char* fname, uintptr_t unstandardized_byte_ct, uintptr_t required_byte_ct, ReadLineStream* rlsp, char** consume_iterp) {
  uintptr_t linebuf_size;
  if (StandardizeLinebufSize(unstandardized_byte_ct, required_byte_ct, &linebuf_size)) {
    return kPglRetNomem;
  }
  return InitRLstreamRaw(fname, linebuf_size, rlsp, consume_iterp);
}

HEADER_INLINE PglErr SizemaxAndInitRLstreamRaw(const char* fname, ReadLineStream* rlsp, char** consume_iterp) {
  uintptr_t linebuf_size;
  if (StandardizeLinebufSizemax(kMaxMediumLine + 1, &linebuf_size)) {
    return kPglRetNomem;
  }
  return InitRLstreamRaw(fname, linebuf_size, rlsp, consume_iterp);
}

HEADER_INLINE PglErr LboundedSizemaxAndInitRLstreamRaw(const char* fname, uintptr_t required_byte_ct, ReadLineStream* rlsp, char** consume_iterp) {
  uintptr_t linebuf_size;
  if (StandardizeLinebufSizemax(required_byte_ct, &linebuf_size)) {
    return kPglRetNomem;
  }
  return InitRLstreamRaw(fname, linebuf_size, rlsp, consume_iterp);
}

PglErr AdvanceRLstream(ReadLineStream* rlsp, char** consume_iterp);

// consume_iter must be at the first unread byte.
HEADER_INLINE PglErr RlsPostlfNext(ReadLineStream* rlsp, char** consume_iterp) {
  if (*consume_iterp != rlsp->consume_stop) {
    return kPglRetSuccess;
  }
  return AdvanceRLstream(rlsp, consume_iterp);
}

// If read_iter was not previously advanced to the byte past the end of the
// most recently read line.
HEADER_INLINE PglErr RlsNext(ReadLineStream* rlsp, char** consume_iterp) {
  *consume_iterp = AdvPastDelim(*consume_iterp, '\n');
  return RlsPostlfNext(rlsp, consume_iterp);
}

// Also skips leading tabs/spaces.
HEADER_INLINE PglErr RlsNextLstrip(ReadLineStream* rlsp, char** consume_iterp) {
  *consume_iterp = AdvPastDelim(*consume_iterp, '\n');
  PglErr reterr = RlsPostlfNext(rlsp, consume_iterp);
  if (reterr) {
    return reterr;
  }
  *consume_iterp = FirstNonTspace(*consume_iterp);
  return kPglRetSuccess;
}

HEADER_INLINE PglErr RlsNextLstripK(ReadLineStream* rlsp, const char** consume_iterp) {
  return RlsNextLstrip(rlsp, K_CAST(char**, consume_iterp));
}

PglErr RlsNextNonemptyLstrip(ReadLineStream* rlsp, uintptr_t* line_idx_ptr, char** consume_iterp);

HEADER_INLINE PglErr RlsNextNonemptyLstripK(ReadLineStream* rlsp, uintptr_t* line_idx_ptr, const char** consume_iterp) {
  return RlsNextNonemptyLstrip(rlsp, line_idx_ptr, K_CAST(char**, consume_iterp));
}

// On success, consume_iter may be anywhere in the last skipped line; only
// guarantee is that a subsequent ReadNextLineFromRLstreamRaw() does what you'd
// expect.
// (The Nz function actually handles skip_ct == 0 properly, but we'd rather
// inline the skip_ct == 0 common-case check.)
PglErr RlsSkipNz(uintptr_t skip_ct, ReadLineStream* rlsp, char** consume_iterp);

HEADER_INLINE PglErr RlsSkip(uintptr_t skip_ct, ReadLineStream* rlsp, char** consume_iterp) {
  if (!skip_ct) {
    return kPglRetSuccess;
  }
  return RlsSkipNz(skip_ct, rlsp, consume_iterp);
}

HEADER_INLINE PglErr RlsSkipK(uintptr_t skip_ct, ReadLineStream* rlsp, const char** consume_iterp) {
  return RlsSkip(skip_ct, rlsp, K_CAST(char**, consume_iterp));
}

// If new_fname == nullptr, this just rewinds the current file.
// If it's non-null, the current file is closed and the specified file is
// opened for reading instead.  The filename is not copied, so it's unsafe to
// change the memory it points to before the first read call.
PglErr RetargetRLstreamRaw(const char* new_fname, ReadLineStream* rlsp, char** consume_iterp);

HEADER_INLINE PglErr RetargetRLstreamRawK(const char* new_fname, ReadLineStream* rlsp, const char** consume_iterp) {
  return RetargetRLstreamRaw(new_fname, rlsp, K_CAST(char**, consume_iterp));
}

HEADER_INLINE PglErr RewindRLstreamRaw(ReadLineStream* rlsp, char** consume_iterp) {
  return RetargetRLstreamRaw(nullptr, rlsp, consume_iterp);
}

HEADER_INLINE PglErr RewindRLstreamRawK(ReadLineStream* rlsp, const char** consume_iterp) {
  return RewindRLstreamRaw(rlsp, K_CAST(char**, consume_iterp));
}

// ***** end insane API function block *****


HEADER_INLINE PglErr InitRLstream(const char* fname, uintptr_t max_line_blen, ReadLineStream* rlsp, char** consume_iterp, char** line_lastp) {
  PglErr reterr = InitRLstreamRaw(fname, max_line_blen, rlsp, consume_iterp);
  *line_lastp = *consume_iterp;
  **line_lastp = '\0';
  return reterr;
}

HEADER_INLINE PglErr ReadNextLineFromRLstream(ReadLineStream* rlsp, char** consume_iterp, char** line_lastp) {
  *consume_iterp = &((*line_lastp)[1]);
  if (*consume_iterp == rlsp->consume_stop) {
    PglErr reterr = AdvanceRLstream(rlsp, consume_iterp);
    if (reterr) {
      return reterr;
    }
  }
  *line_lastp = S_CAST(char*, rawmemchr(*consume_iterp, '\n'));
  **line_lastp = '\0';
  return kPglRetSuccess;
}

HEADER_INLINE PglErr RewindRLstream(ReadLineStream* rlsp, char** consume_iterp, char** line_lastp) {
  PglErr reterr = RewindRLstreamRaw(rlsp, consume_iterp);
  *line_lastp = *consume_iterp;
  return reterr;
}

PglErr CleanupRLstream(ReadLineStream* rlsp);

void RLstreamErrPrint(const char* file_descrip, ReadLineStream* rlsp, PglErr* reterr_ptr);

HEADER_INLINE unsigned char* RLstreamMemStart(ReadLineStream* rlsp) {
  return R_CAST(unsigned char*, rlsp->syncp);
}

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_DECOMPRESS_H__
