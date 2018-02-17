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

#ifdef __cplusplus
namespace plink2 {
#endif

// Also sets 128k read buffer.
PglErr gzopen_read_checked(const char* fname, gzFile* gzf_ptr);

// This sets loadbuf[loadbuf_size - 1] to ' ', just because.
// loadbuf_size assumed to be either exactly kMaxMediumLine (in which case
// any longer line is treated as pathological), or strictly larger (in which
// case we report an out-of-memory error when gzgets blows the buffer, unless
// loadbuf_size == kMaxLongLine, which is close to 2GB).
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


/*
// While a separate non-compressing writer thread is practically useless (since
// the OS does its own scheduling, etc. of the writes anyway), a separate
// read-ahead thread is frequently useful since the typical read_ahead_kb = 64
// Linux kernel setting is sometimes smaller than we really want, and of course
// the ability to decompress-ahead speaks for itself.
//
// Possible todos: tabix index support (caller must declare the entire sequence
// of virtual offset ranges they want to read upfront), and optional
// multithreaded decompression of BGZF and seekable-Zstd files.  Note that
// these require us to stop using the Zstd zlib wrapper, so they're unlikely to
// be worth the time investment before the initial plink 2.0 beta release.

// Suppose the maximum line length is 1000 MB.  Then,
// * We open the file, and initialize a shared buffer of size 1001 MB, where
//   the read-head, available-head, and consume pointers start at position 0,
//   and the read-stop pointer is at the buffer end.
// * We launch a reader thread which owns the file object.  It decompresses the
//   first 1 MB into positions [0, 1 MB) in the shared buffer, and then checks
//   for a '\n'.  If there isn't one, it decompresses the next 1 MB into
//   positions [1 MB, 2 MB) in the buffer, and then rechecks, etc., and may
//   eventually throw a LongLine error.
// * Otherwise... suppose that the last '\n' is at position 3145000.  We had
//   read (and possibly decompressed) up to (but not including) position
//   3145728 in the buffer, so we set the available-head pointer to position
//   3145001, and set the line-available event.  Since this is past 1 MB, the
//   next line is not guaranteed to fit in the buffer, so while we can continue
//   reading forward, we need to be prepared to copy the new bytes back to the
//   beginning.
// * The consumer's initial RLstreamGet() call probably happens before the
//   first read operation has completed.  It sees that the consume pointer is
//   equal to available-head, so it waits on the line-available event.  Once
//   that's set, it returns.  (Note that no null-terminators are inserted; the
//   only properties that make this a 'ReadLineStream' instead of a generic
//   decompressor are (i) each successful return is guaranteed to end in '\n'
//   (one will be appended at eof if necessary), and (ii) the reader thread may
//   throw a LongLine error.)
// * The consumer sets the consume-done event once the consume pointer has
//   caught up to available-head, and waits on line-available.
// * The reader thread reads/decompresses the next 1 MB, checks for '\n',
//   reads an additional 1 MB if it's absent, etc., until it (i) finds another
//   '\n' or eof (in which case a '\n' is appended), or (ii) hits the read-stop
//   pointer.  Then it waits on the consume-done event; once that returns
//   (usually instant if decompressing),
//   - in case (ii), it memmoves bytes back if read-stop == end-of-buffer, and
//     sets read-stop := end-of-buffer otherwise.  Either way, it continues
//     with reading until '\n'/eof or a LongLine error.
//   - in case (i), the available-head pointer is set one byte past the last
//     '\n', the line-available event is set, and if eof we set at_eof and
//     return from the thread function.
//     If not eof, we check the old position of the available-head pointer.  If
//     it's before the 1 MB mark, we continue reading forward.  Otherwise, we
//     memcpy the trailing bytes (at or past the new available-head position)
//     to the beginning of the buffer, set the read-stop pointer to the old
//     available-head position rounded down to a cacheline boundary, and
//     continue reading from the end of the copied trailing bytes.

// (tested a few different values for this, 1 MB appears to work well on the
// systems we care most about)
CONSTU31(kDecompressChunkSize, 1048576);

// To minimize false (or true) sharing penalties, these values shouldn't change
// much; only the things they point to should be frequently changing.
typedef struct ReadLineStreamStruct {
  // ** These variables are owned by the reader thread.
  gzFile gz_infile;
  // This is aimed at (usually uncompressed) text files, so just use char*
  // instead of unsigned char*.
  char* buf;
  char* buf_end;

  // ** These variables are shared between the reader and the consumer.
  char* available_head;
#ifdef _WIN32
  HANDLE line_available_event;
  HANDLE consume_done_event;
#else
  pthread_mutex_t sync_mutex;
  // these two variables are guarded by sync_mutex
  pthread_cond_t line_available_condvar;
  pthread_cond_t consume_done_condvar;
#endif
  // 1 = rewind, 2 = exit immediately
  uint32_t rewind_or_interrupt;
  uint32_t at_eof;

  PglErr reterr;

  // ** These variables are owned by the consumer.
  pthread_t read_thread;
#ifndef _WIN32
  uint32_t sync_init_state;
#endif
} ReadLineStream;

void PreinitRLstream(ReadLineStream* rlsp);
*/

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_DECOMPRESS_H__
