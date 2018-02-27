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


#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr gzopen_read_checked(const char* fname, gzFile* gzf_ptr) {
  *gzf_ptr = gzopen(fname, FOPEN_RB);
  if (!(*gzf_ptr)) {
    logputs("\n");
    logerrprintfww(kErrprintfFopen, fname);
    return kPglRetOpenFail;
  }
  if (gzbuffer(*gzf_ptr, 131072)) {
    return kPglRetNomem;
  }
  return kPglRetSuccess;
}

PglErr GzopenAndSkipFirstLines(const char* fname, uint32_t lines_to_skip, uintptr_t loadbuf_size, char* loadbuf, gzFile* gzf_ptr) {
  PglErr reterr = gzopen_read_checked(fname, gzf_ptr);
  if (reterr) {
    return reterr;
  }
  loadbuf[loadbuf_size - 1] = ' ';
  for (uint32_t line_idx = 1; line_idx <= lines_to_skip; ++line_idx) {
    if (!gzgets(*gzf_ptr, loadbuf, loadbuf_size)) {
      if (gzeof(*gzf_ptr)) {
        logerrprintfww("Error: Fewer lines than expected in %s.\n", fname);
        return kPglRetInconsistentInput;
      }
      return kPglRetReadFail;
    }
    if (!loadbuf[loadbuf_size - 1]) {
      if ((loadbuf_size == kMaxMediumLine) || (loadbuf_size == kMaxLongLine)) {
        logerrprintfww("Error: Line %u of %s is pathologically long.\n", line_idx, fname);
        return kPglRetMalformedInput;
      }
      return kPglRetNomem;
    }
  }
  return kPglRetSuccess;
}


void PreinitGzTokenStream(GzTokenStream* gtsp) {
  gtsp->gz_infile = nullptr;
}

PglErr InitGzTokenStream(const char* fname, GzTokenStream* gtsp, char* buf_start) {
  PglErr reterr = gzopen_read_checked(fname, &(gtsp->gz_infile));
  if (reterr) {
    return reterr;
  }
  gtsp->buf_start = buf_start;
  gtsp->read_iter = &(buf_start[kMaxMediumLine]);
  gtsp->buf_end = gtsp->read_iter;
  gtsp->buf_end[0] = '0';  // force initial load
  return kPglRetSuccess;
}

char* AdvanceGzTokenStream(GzTokenStream* gtsp, uint32_t* token_slen_ptr) {
  char* token_start = gtsp->read_iter;
  char* buf_end = gtsp->buf_end;
 AdvanceGzTokenStream_restart:
  while (ctou32(*token_start) <= ' ') {
    ++token_start;
  }
  while (1) {
    if (token_start < buf_end) {
      char* token_end = &(token_start[1]);
      while (ctou32(*token_end) > ' ') {
        ++token_end;
      }
      const uint32_t token_slen = token_end - token_start;
      if (token_end < buf_end) {
        *token_slen_ptr = token_slen;
        gtsp->read_iter = &(token_end[1]);
        return token_start;
      }
      if (token_slen > kMaxMediumLine) {
        *token_slen_ptr = UINT32_MAX;
        return nullptr;
      }
      char* new_token_start = &(gtsp->buf_start[kMaxMediumLine - token_slen]);
      memcpy(new_token_start, token_start, token_slen);
      token_start = new_token_start;
    } else {
      token_start = &(gtsp->buf_start[kMaxMediumLine]);
    }
    char* load_start = &(gtsp->buf_start[kMaxMediumLine]);
    const int32_t bufsize = gzread(gtsp->gz_infile, load_start, kMaxMediumLine);
    if (bufsize < 0) {
      *token_slen_ptr = 0xfffffffeU;
      return nullptr;
    }
    buf_end = &(load_start[S_CAST(uint32_t, bufsize)]);
    buf_end[0] = ' ';
    buf_end[1] = '0';
    gtsp->buf_end = buf_end;
    if (!bufsize) {
      if (!gzeof(gtsp->gz_infile)) {
        *token_slen_ptr = 0xfffffffeU;
        return nullptr;
      }
      // bufsize == 0, eof
      if (token_start == load_start) {
        *token_slen_ptr = 0;
        return nullptr;
      }
      gtsp->read_iter = load_start;
      *token_slen_ptr = load_start - token_start;
      return token_start;
    }
    if (token_start == load_start) {
      goto AdvanceGzTokenStream_restart;
    }
  }
}

BoolErr CloseGzTokenStream(GzTokenStream* gtsp) {
  if (!gtsp->gz_infile) {
    return 0;
  }
  return gzclose_null(&(gtsp->gz_infile));
}


void PreinitRLstream(ReadLineStream* rlsp) {
  rlsp->gz_infile = nullptr;
#ifndef _WIN32
  rlsp->sync_init_state = 0;
#endif
}

THREAD_FUNC_DECL ReadLineStreamThread(void* arg) {
  ReadLineStream* context = R_CAST(ReadLineStream*, arg);
  ReadLineStreamSync* syncp = context->syncp;
  gzFile gz_infile = context->gz_infile;
  char* buf = context->buf;
  char* buf_end = context->buf_end;
  char* cur_block_start = buf;
  char* read_head = buf;

  // We can either be reading/decompressing into memory past the bytes passed
  // to the consumer, or we can be doing it before those bytes.
  // In the first case, read_stop is buf_end, but it gets changed to the
  // latest value of consume_tail when we return to the front of the buffer.
  // In the second case, read_stop is the position of the first passed byte.
  char* read_stop = buf_end;
#ifdef _WIN32
  CRITICAL_SECTION* critical_sectionp = &syncp->critical_section;
  HANDLE reader_progress_event = context->reader_progress_event;
  HANDLE consumer_progress_event = context->consumer_progress_event;
#else
  pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
  pthread_cond_t* reader_progress_condvarp = &syncp->reader_progress_condvar;
  pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
#endif
  while (1) {
    RlsInterrupt interrupt = kRlsInterruptNone;
    PglErr reterr;
    RlsInterrupt min_interrupt;
    while (1) {
      uintptr_t read_attempt_size = read_stop - read_head;
      if (!read_attempt_size) {
        if ((read_stop != buf_end) || (cur_block_start != buf)) {
          goto ReadLineStreamThread_SYNC_BEFORE_CONTINUING_READ;
        }
        goto ReadLineStreamThread_LONG_LINE;
      }
      if (read_attempt_size > kDecompressChunkSize) {
        read_attempt_size = kDecompressChunkSize;
      }
      int32_t bytes_read = gzread(gz_infile, read_head, read_attempt_size);
      char* cur_read_end = &(read_head[S_CAST(uint32_t, bytes_read)]);
      if (bytes_read < S_CAST(int32_t, S_CAST(uint32_t, read_attempt_size))) {
        if (!gzeof(gz_infile)) {
          goto ReadLineStreamThread_READ_FAIL;
        }
        read_head = cur_read_end;
        if (cur_block_start != read_head) {
          if (read_head[-1] != '\n') {
            // Append '\n' so consumer can always use rawmemchr(., '\n') to
            // find the end of the current line.
            *read_head++ = '\n';
          }
        }
        goto ReadLineStreamThread_EOF;
      }
      char* last_lf = memrchr_expected_far(read_head, '\n', read_attempt_size);
      if (last_lf) {
        char* next_available_end = &(last_lf[1]);
#ifdef _WIN32
        EnterCriticalSection(critical_sectionp);
#else
        pthread_mutex_lock(sync_mutexp);
#endif
        interrupt = syncp->interrupt;
        if (interrupt != kRlsInterruptNone) {
          goto ReadLineStreamThread_INTERRUPT;
        }
        char* latest_consume_tail = syncp->consume_tail;
        const uint32_t all_later_bytes_consumed = (latest_consume_tail <= cur_block_start);
        const uint32_t return_to_start = all_later_bytes_consumed && (latest_consume_tail >= &(buf[kDecompressChunkSize]));
        syncp->available_end = next_available_end;
        if (return_to_start) {
          syncp->cur_circular_end = next_available_end;
        }
#ifdef _WIN32
        LeaveCriticalSection(critical_sectionp);
        SetEvent(reader_progress_event);
#else
        pthread_cond_signal(reader_progress_condvarp);
        pthread_mutex_unlock(sync_mutexp);
#endif
        if (return_to_start) {
          // Best to return to the beginning of the buffer.
          // (Note that read_attempt_size is guaranteed to be
          // <= kDecompressChunkSize.)
          const uintptr_t trailing_byte_ct = cur_read_end - next_available_end;
          memcpy(buf, next_available_end, trailing_byte_ct);
          cur_block_start = buf;
          read_head = &(buf[trailing_byte_ct]);
          // May as well reduce false sharing risk.
          read_stop = R_CAST(char*, RoundDownPow2(R_CAST(uintptr_t, latest_consume_tail), kCacheline));
          continue;
        }
        if (all_later_bytes_consumed) {
          read_stop = buf_end;
        } else {
          read_stop = R_CAST(char*, RoundDownPow2(R_CAST(uintptr_t, latest_consume_tail), kCacheline));
        }
        cur_block_start = next_available_end;
      }
      read_head = cur_read_end;
    }
    while (0) {
    ReadLineStreamThread_READ_FAIL:
      min_interrupt = kRlsInterruptShutdown;
      reterr = kPglRetReadFail;
      break;
    ReadLineStreamThread_LONG_LINE:
      min_interrupt = kRlsInterruptShutdown;
      reterr = kPglRetLongLine;
      break;
    ReadLineStreamThread_EOF:
      min_interrupt = kRlsInterruptRewind;
      reterr = kPglRetSkipped;
      break;
    }
    // We need to wait for a message from the consumer before we can usefully
    // proceed.
    // More precisely:
    // * In the eof subcase, we're waiting for either a rewind or shutdown
    //   request.
    // * In the error subcase, we're just waiting for the shutdown request.

    // Pass the error code back.
#ifdef _WIN32
    EnterCriticalSection(critical_sectionp);
#else
    pthread_mutex_lock(sync_mutexp);
#endif
    syncp->reterr = reterr;
    interrupt = syncp->interrupt;
    if (interrupt >= min_interrupt) {
      // It's our lucky day: we don't need to wait again.
      goto ReadLineStreamThread_INTERRUPT;
    }
    if (reterr == kPglRetSkipped) {
      syncp->available_end = read_head;
    }
#ifdef _WIN32
    LeaveCriticalSection(critical_sectionp);
    SetEvent(reader_progress_event);
    while (1) {
      WaitForSingleObject(consumer_progress_event, INFINITE);
      EnterCriticalSection(critical_sectionp);
      interrupt = syncp->interrupt;
      if (interrupt >= min_interrupt) {
        break;
      }
      LeaveCriticalSection(critical_sectionp);
    }
#else
    pthread_cond_signal(reader_progress_condvarp);
    do {
      pthread_cond_wait(consumer_progress_condvarp, sync_mutexp);
      interrupt = syncp->interrupt;
    } while (interrupt < min_interrupt);
#endif
  ReadLineStreamThread_INTERRUPT:
    // must be in critical section here, or be holding the mutex.
    if (interrupt == kRlsInterruptRewind) {
      syncp->consume_tail = buf;
      syncp->cur_circular_end = nullptr;
      syncp->available_end = buf;
      syncp->interrupt = kRlsInterruptNone;
      syncp->reterr = kPglRetSuccess;
    }
#ifdef _WIN32
    LeaveCriticalSection(critical_sectionp);
#else
    pthread_mutex_unlock(sync_mutexp);
#endif
    if (interrupt == kRlsInterruptShutdown) {
      // possible todo: close the file here
      THREAD_RETURN;
    }
    assert(interrupt == kRlsInterruptRewind);
    if (gzrewind(gz_infile)) {
      goto ReadLineStreamThread_READ_FAIL;
    }
    cur_block_start = buf;
    read_head = buf;
    read_stop = buf_end;
    continue;
    {
      char* lbound;
    ReadLineStreamThread_SYNC_BEFORE_CONTINUING_READ:
      // We cannot continue reading forward.  Cases:
      // 1. read_stop == buf_end, cur_block_start != buf.  This means we're
      //    in the middle of reading/decompressing a really long line, and
      //    want to wait for consume_tail == cur_block_start so we can memmove
      //    all the bytes back and continue reading forward.
      // 2. read_stop == buf_end, cur_block_start == buf.  We failed with a
      //    long-line error here.
      // 3. read_stop < buf_end (usual case).  This means the consumer may
      //    not be done handling some bytes-in-front we handed off earlier.  We
      //    are waiting for consume_tail <= cur_block_start, which means all
      //    bytes in front have been consumed and we're free to continue
      //    reading forward.
      lbound = (read_stop == buf_end)? cur_block_start : buf;
      char* latest_consume_tail;
#ifdef _WIN32
      do {
        WaitForSingleObject(consumer_progress_event, INFINITE);
        EnterCriticalSection(critical_sectionp);
        interrupt = syncp->interrupt;
        if (interrupt != kRlsInterruptNone) {
          goto ReadLineStreamThread_INTERRUPT;
        }
        latest_consume_tail = syncp->consume_tail;
        LeaveCriticalSection(critical_sectionp);
      } while ((latest_consume_tail > cur_block_start) || (latest_consume_tail < lbound));
#else
      pthread_mutex_lock(sync_mutexp);
      interrupt = syncp->interrupt;
      if (interrupt != kRlsInterruptNone) {
        goto ReadLineStreamThread_INTERRUPT;
      }
      do {
        pthread_cond_wait(consumer_progress_condvarp, sync_mutexp);
        interrupt = syncp->interrupt;
        if (interrupt != kRlsInterruptNone) {
          goto ReadLineStreamThread_INTERRUPT;
        }
        latest_consume_tail = syncp->consume_tail;
      } while ((latest_consume_tail > cur_block_start) || (latest_consume_tail < lbound));
      pthread_mutex_unlock(sync_mutexp);
#endif
      if (read_stop == buf_end) {
        const uintptr_t cur_len = read_head - cur_block_start;
        memmove(buf, cur_block_start, cur_len);
        cur_block_start = buf;
        read_head = &(buf[cur_len]);
      } else {
        read_stop = buf_end;
      }
    }
  }
}

PglErr InitRLstreamRaw(const char* fname, uintptr_t max_line_blen, ReadLineStream* rlsp, char** consume_iterp) {
  PglErr reterr = kPglRetSuccess;
  {
    // To avoid a first-line special case, we start consume_iter at position
    // -1, pointing to a \n.  We make this safe by adding 1 to the previous
    // bigstack allocation size.  To simplify rewind, we also don't let the
    // reader thread overwrite this byte.
    ReadLineStreamSync* syncp = S_CAST(ReadLineStreamSync*, bigstack_alloc(sizeof(ReadLineStreamSync) + 1));
    if (!syncp) {
      goto InitRLstreamRaw_ret_NOMEM;
    }
    char* buf;
    if (bigstack_alloc_c(max_line_blen + kDecompressChunkSize, &buf)) {
      goto InitRLstreamRaw_ret_NOMEM;
    }
    reterr = gzopen_read_checked(fname, &rlsp->gz_infile);
    if (reterr) {
      goto InitRLstreamRaw_ret_1;
    }
    buf[-1] = '\n';
    rlsp->consume_stop = buf;
    rlsp->syncp = syncp;
    rlsp->buf = buf;
    rlsp->buf_end = R_CAST(char*, g_bigstack_base);
    syncp->consume_tail = buf;
    syncp->cur_circular_end = nullptr;
    syncp->available_end = buf;
    syncp->reterr = kPglRetSuccess;
    syncp->interrupt = kRlsInterruptNone;
#ifdef _WIN32
    // apparently this can raise a low-memory exception in older Windows
    // versions, but that's not really our problem.
    InitializeCriticalSection(&syncp->critical_section);

    rlsp->reader_progress_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (!rlsp->reader_progress_event) {
      DeleteCriticalSection(&syncp->critical_section);
      goto InitRLstreamRaw_ret_THREAD_CREATE_FAIL;
    }
    rlsp->consumer_progress_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (!rlsp->consumer_progress_event) {
      DeleteCriticalSection(&syncp->critical_section);
      CloseHandle(rlsp->reader_progress_event);
      goto InitRLstreamRaw_ret_THREAD_CREATE_FAIL;
    }
    rlsp->read_thread = R_CAST(HANDLE, _beginthreadex(nullptr, kDefaultThreadStack, ReadLineStreamThread, rlsp, 0, nullptr));
    if (!rlsp->read_thread) {
      DeleteCriticalSection(&syncp->critical_section);
      CloseHandle(rlsp->consumer_progress_event);
      CloseHandle(rlsp->reader_progress_event);
      goto InitRLstreamRaw_ret_THREAD_CREATE_FAIL;
    }
#else
    if (pthread_mutex_init(&syncp->sync_mutex, nullptr)) {
      goto InitRLstreamRaw_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 1;
    if (pthread_cond_init(&syncp->reader_progress_condvar, nullptr)) {
      goto InitRLstreamRaw_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 2;
    if (pthread_cond_init(&syncp->consumer_progress_condvar, nullptr)) {
      goto InitRLstreamRaw_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 3;
    if (pthread_create(&rlsp->read_thread, &g_smallstack_thread_attr, ReadLineStreamThread, rlsp)) {
      goto InitRLstreamRaw_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 4;
#endif
    *consume_iterp = &(buf[-1]);
  }
  while (0) {
  InitRLstreamRaw_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  InitRLstreamRaw_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 InitRLstreamRaw_ret_1:
  return reterr;
}

PglErr AdvanceRLstream(ReadLineStream* rlsp, char** consume_iterp) {
  char* consume_iter = *consume_iterp;
  ReadLineStreamSync* syncp = rlsp->syncp;
#ifdef _WIN32
  CRITICAL_SECTION* critical_sectionp = &syncp->critical_section;
  HANDLE consumer_progress_event = rlsp->consumer_progress_event;
  while (1) {
    EnterCriticalSection(critical_sectionp);
    const PglErr reterr = syncp->reterr;
    if (S_CAST(uint32_t, reterr) > S_CAST(uint32_t, kPglRetSkipped)) {
      LeaveCriticalSection(critical_sectionp);
      // No need to set consumer_progress event here, just let the cleanup
      // routine take care of that.
      return reterr;
    }
    char* available_end = syncp->available_end;
    if (consume_iter != available_end) {
      char* cur_circular_end = syncp->cur_circular_end;
      if (cur_circular_end) {
        if (cur_circular_end == consume_iter) {
          char* buf = rlsp->buf;
          syncp->consume_tail = buf;
          syncp->cur_circular_end = nullptr;
          LeaveCriticalSection(critical_sectionp);
          SetEvent(consumer_progress_event);
          *consume_iterp = buf;
          rlsp->consume_stop = available_end;
          return kPglRetSuccess;
        }
        rlsp->consume_stop = cur_circular_end;
      } else {
        rlsp->consume_stop = available_end;
      }
      syncp->consume_tail = consume_iter;
      LeaveCriticalSection(critical_sectionp);
      // We could set the consumer_progress event here, but it's not really
      // necessary.
      // SetEvent(consumer_progress_event);
      return kPglRetSuccess;
    }
    syncp->consume_tail = consume_iter;
    LeaveCriticalSection(critical_sectionp);
    // We've processed all the consume-ready bytes...
    if (reterr != kPglRetSuccess) {
      // ...and we're at eof.  Don't set consumer_progress event here; let that
      // wait until cleanup or rewind.
      return reterr;
    }
    // ...and there's probably more.
    SetEvent(consumer_progress_event);
    WaitForSingleObject(rlsp->reader_progress_event, INFINITE);
  }
#else
  pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
  pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
  pthread_cond_t* reader_progress_condvarp = &syncp->reader_progress_condvar;
  pthread_mutex_lock(sync_mutexp);
  while (1) {
    const PglErr reterr = syncp->reterr;
    if (S_CAST(uint32_t, reterr) > S_CAST(uint32_t, kPglRetSkipped)) {
      pthread_mutex_unlock(sync_mutexp);
      return reterr;
    }
    char* available_end = syncp->available_end;
    if (consume_iter != available_end) {
      char* cur_circular_end = syncp->cur_circular_end;
      if (cur_circular_end) {
        if (cur_circular_end == consume_iter) {
          char* buf = rlsp->buf;
          syncp->consume_tail = buf;
          syncp->cur_circular_end = nullptr;
          pthread_cond_signal(consumer_progress_condvarp);
          pthread_mutex_unlock(sync_mutexp);
          *consume_iterp = buf;
          rlsp->consume_stop = available_end;
          return kPglRetSuccess;
        }
        rlsp->consume_stop = cur_circular_end;
      } else {
        rlsp->consume_stop = available_end;
      }
      syncp->consume_tail = consume_iter;
      // pthread_cond_signal(consumer_progress_condvarp);
      pthread_mutex_unlock(sync_mutexp);
      return kPglRetSuccess;
    }
    syncp->consume_tail = consume_iter;
    // We've processed all the consume-ready bytes...
    if (reterr != kPglRetSuccess) {
      // ...and we're at eof.
      pthread_mutex_unlock(sync_mutexp);
      return reterr;
    }
    // ...and there's probably more.
    pthread_cond_signal(consumer_progress_condvarp);
    // no need for an explicit spurious-wakeup check, we'll check the progress
    // condition (available_end has advanced, or we have a read error) anyway
    // and get back here if it isn't satisfied
    pthread_cond_wait(reader_progress_condvarp, sync_mutexp);
  }
#endif
}

PglErr RewindRLstreamRaw(ReadLineStream* rlsp, char** consume_iterp) {
  ReadLineStreamSync* syncp = rlsp->syncp;
#ifdef _WIN32
  CRITICAL_SECTION* critical_sectionp = &syncp->critical_section;
  EnterCriticalSection(critical_sectionp);
  const PglErr reterr = syncp->reterr;
  if (reterr != kPglRetSuccess) {
    if (S_CAST(uint32_t, reterr) != S_CAST(uint32_t, kPglRetSkipped)) {
      LeaveCriticalSection(critical_sectionp);
      return reterr;
    }
    // clear eof
    syncp->reterr = kPglRetSuccess;
  }
  syncp->interrupt = kRlsInterruptRewind;
  LeaveCriticalSection(critical_sectionp);
  SetEvent(rlsp->consumer_progress_event);
#else
  pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
  pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
  pthread_mutex_lock(sync_mutexp);
  const PglErr reterr = syncp->reterr;
  if (reterr != kPglRetSuccess) {
    if (S_CAST(uint32_t, reterr) != S_CAST(uint32_t, kPglRetSkipped)) {
      pthread_mutex_unlock(sync_mutexp);
      return reterr;
    }
    // clear eof
    syncp->reterr = kPglRetSuccess;
  }
  syncp->interrupt = kRlsInterruptRewind;
  pthread_cond_signal(consumer_progress_condvarp);
  pthread_mutex_unlock(sync_mutexp);
#endif
  char* buf = rlsp->buf;
  *consume_iterp = &(buf[-1]);
  rlsp->consume_stop = buf;
  return kPglRetSuccess;
}

PglErr CleanupRLstream(ReadLineStream* rlsp) {
  PglErr reterr = kPglRetSuccess;
#ifdef _WIN32
  if (rlsp->read_thread) {
    ReadLineStreamSync* syncp = rlsp->syncp;
    CRITICAL_SECTION* critical_sectionp = &syncp->critical_section;
    EnterCriticalSection(critical_sectionp);
    syncp->interrupt = kRlsInterruptShutdown;
    LeaveCriticalSection(critical_sectionp);
    SetEvent(rlsp->consumer_progress_event);
    WaitForSingleObject(rlsp->read_thread, INFINITE);
    DeleteCriticalSection(critical_sectionp);
    rlsp->read_thread = nullptr;  // make it safe to call this multiple times
    CloseHandle(rlsp->consumer_progress_event);
    CloseHandle(rlsp->reader_progress_event);
  }
#else
  const uint32_t sync_init_state = rlsp->sync_init_state;
  if (sync_init_state) {
    ReadLineStreamSync* syncp = rlsp->syncp;
    pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
    pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
    if (sync_init_state == 4) {
      pthread_mutex_lock(sync_mutexp);
      syncp->interrupt = kRlsInterruptShutdown;
      pthread_cond_signal(consumer_progress_condvarp);
      pthread_mutex_unlock(sync_mutexp);
      pthread_join(rlsp->read_thread, nullptr);
    }
    pthread_mutex_destroy(sync_mutexp);
    if (sync_init_state > 1) {
      pthread_cond_destroy(&syncp->reader_progress_condvar);
      if (sync_init_state > 2) {
        pthread_cond_destroy(consumer_progress_condvarp);
      }
    }
    rlsp->sync_init_state = 0;  // make it safe to call this multiple times
  }
#endif
  if (rlsp->gz_infile) {
    if (gzclose_null(&rlsp->gz_infile)) {
      reterr = kPglRetReadFail;
    }
  }
  return reterr;
}

#ifdef __cplusplus
}
#endif
