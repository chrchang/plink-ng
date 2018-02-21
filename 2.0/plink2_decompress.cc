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


/*
void PreinitRLstream(ReadLineStream* rlsp) {
  rlsp->gz_infile = nullptr;
  rlsp->sync_init_state = 0;
}

THREAD_FUNC_DECL ReadLineStreamThread(void* arg) {
  ReadLineStream* context = R_CAST(ReadLineStream*, arg);
  ReadLineStreamSync* syncp = context->syncp;
  gzFile gz_infile = context->gz_infile;
  char* buf = context->buf;
  char* buf_end = &(context->buf_end);
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
  HANDLE reader_progress_event = syncp->reader_progress_event;
  HANDLE consumer_progress_event = syncp->consumer_progress_event;
#else
  pthread_mutex_t* sync_mutexp = &syncp->sync_mutex;
  pthread_cond_t* reader_progress_condvarp = &syncp->reader_progress_condvar;
  pthread_cond_t* consumer_progress_condvarp = &syncp->consumer_progress_condvar;
#endif
  while (1) {
    uint32_t unfinished_line = 0;
    uint32_t consumer_status = 0;
    while (1) {
      uintptr_t read_attempt_size = read_stop - read_head;
      if (!read_attempt_size) {
        if (cur_block_start != buf) {
          if (cur_block_start == read_head) {
          }
          goto ReadLineStreamThread_UNFINISHED_LINE;
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
#ifdef _WIN32
        EnterCriticalSection(critical_sectionp);
#else
        pthread_mutex_lock(sync_mutexp);
#endif
        consumer_status = syncp->consumer_status;
        if (consumer_status > 1) {
          goto ReadLineStreamThread_INTERRUPT;
        }
        syncp->consumer_status = 0;
        syncp->available_end = read_head;
        syncp->reterr = kPglRetSkipped;
#ifdef _WIN32
        LeaveCriticalSection(critical_sectionp);
        SetEvent(reader_progress_event);
#else
        pthread_cond_signal(reader_progress_condvarp);
        pthread_mutex_unlock(sync_mutexp);
#endif
        break;
      }
      char* last_lf = memrchr_expected_far(read_head, '\n', read_attempt_size);
      if (last_lf) {
        char* next_available_end = &(last_lf[1]);
#ifdef _WIN32
        EnterCriticalSection(critical_sectionp);
#else
        pthread_mutex_lock(sync_mutexp);
#endif
        consumer_status = syncp->consumer_status;
        if (consumer_status > 1) {
          goto ReadLineStreamThread_INTERRUPT;
        }
        syncp->consumer_status = 0;
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
          read_stop = latest_consume_tail;
          continue;
        }
        if (all_later_bytes_consumed) {
          read_stop = buf_end;
        } else {
          read_stop = latest_consume_tail;
        }
        cur_block_start = next_available_end;
      }
      read_head = cur_read_end;
    }
    while (0) {
    ReadLineStreamThread_READ_FAIL:
      context->reterr = kPglRetReadFail;
      break;
    ReadLineStreamThread_LONG_LINE:
      if (S_CAST(uintptr_t, buf_end - buf) == (S_CAST(uintptr_t, kMaxLongLine) + kDecompressChunkSize)) {
        context->reterr = kPglRetNomem;
      } else {
        context->reterr = kPglRetMalformedInput;
      }
      break;
    }
#ifdef _WIN32
    SetEvent(reader_progress_event);
#else
    pthread_mutex_lock(sync_mutexp);
    pthread_cond_signal(reader_progress_condvarp);
#endif
    while (0) {
    ReadLineStreamThread_UNFINISHED_LINE:
      // no reader progress in this case
      unfinished_line = 1;
#ifndef _WIN32
      pthread_mutex_lock(sync_mutexp);
#endif
    }
    // We need to wait for a message from the consumer before we can usefully
    // proceed.
#ifdef _WIN32
    WaitForSingleObject(consumer_progress_event, INFINITE);
    EnterCriticalSection(critical_sectionp);
    consumer_status = syncp->consumer_status;
    LeaveCriticalSection(critical_sectionp);
#else
    consumer_status = syncp->consumer_status;
    while (!consumer_status) {
      // spurious wakeup guard
      pthread_cond_wait(consumer_progress_condvarp, sync_mutexp);
      consumer_status = syncp->consumer_status;
    }
    pthread_mutex_unlock(sync_mutexp);
#endif
    if (consumer_status == 1) {
      if (!unfinished_line) {

      }
      const uintptr_t cur_len = read_head - cur_block_start;
      memmove(buf, cur_block_start, cur_len);
      cur_block_start = buf;
      read_head = &(buf[cur_len]);
      continue;
    }
    while (0) {
    ReadLineStreamThread_INTERRUPT:
#ifdef _WIN32
      LeaveCriticalSection(critical_sectionp);
#else
      pthread_mutex_unlock(sync_mutexp);
#endif
    }
    if (consumer_status == 3) {
      // todo: close the file here
      THREAD_RETURN;
    }
    // case 2: rewind
    if (gzrewind(gz_infile)) {
      goto ReadLineStreamThread_READ_FAIL;
    }
    // todo: set all pointers back to beginning
  }
}

PglErr InitRLstream(const char* fname, uintptr_t max_line_blen, ReadLineStream* rlsp, char** consume_iterp) {
  PglErr reterr = kPglRetSuccess;
  {
    // +1 since we may need to append an \n on the last line
    if (bigstack_alloc_c(max_line_blen + kDecompressChunkSize + 1, &rlsp->buf)) {
      goto InitRLstream_ret_NOMEM;
    }
    rlsp->buf_end = g_bigstack_base;
    rlsp->available_head = rlsp->buf;
    rlsp->rewind_or_interrupt = 0;
    rlsp->at_eof = 0;
    rlsp->reterr = kPglRetSuccess;
#ifdef _WIN32
    rlsp->reader_progress_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (!rlsp->reader_progress_event) {
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->consumer_progress_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (!rlsp->consumer_progress_event) {
      CloseHandle(rlsp->reader_progress_event);
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->read_thread = R_CAST(HANDLE, _beginthreadex(nullptr, kDefaultThreadStack, ReadLineStreamThread, rlsp, 0, nullptr));
    if (!rlsp->read_thread) {
      CloseHandle(rlsp->consumer_progress_event);
      CloseHandle(rlsp->reader_progress_event);
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
#else
    if (pthread_mutex_init(&rlsp->sync_mutex, nullptr)) {
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 1;
    if (pthread_cond_init(&rlsp->reader_progress_condvar, nullptr)) {
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 2;
    if (pthread_cond_init(&rlsp->consumer_progress_condvar, nullptr)) {
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 3;
    if (pthread_create(&rlsp->read_thread, &g_smallstack_thread_attr, ReadLineStreamThread, rlsp)) {
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 4;
#endif
    *consume_iterp = rlsp->buf;
  }
  while (0) {
  InitRLstream_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  InitRLstream_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  InitRLstream_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  return reterr;
}

PglErr AdvanceRLstream(ReadLineStream* rlsp, char** read_iterp) {
  char* read_iter = *read_iterp;
#ifdef _WIN32
  CRITICAL_SECTION* critical_sectionp = &rlsp->critical_section;
  HANDLE consumer_progress_event = rlsp->consumer_progress_event;
  while (1) {
    EnterCriticalSection(critical_sectionp);
    const PglErr reterr = rlsp->reterr;
    if (S_CAST(uint32_t, reterr) > S_CAST(uint32_t, kPglRetSkipped)) {
      LeaveCriticalSection(critical_sectionp);
      return reterr;
    }
    char* available_end = rlsp->available_end;
    if (read_iter != available_end) {
      char* cur_circular_end = rlsp->cur_circular_end;
      if (cur_circular_end) {
        if (cur_circular_end == read_iter) {
          char* buf = rlsp->buf;
          read_iter = buf;
          rlsp->consume_tail = buf;
          rlsp->cur_circular_end = nullptr;
          LeaveCriticalSection(critical_sectionp);
          SetEvent(consumer_progress_event);
          return kPglRetSuccess;
        }
        rlsp->consume_stop = cur_circular_end;
      } else {
        rlsp->consume_stop = available_end;
      }
      rlsp->consume_tail = read_iter;
      LeaveCriticalSection(critical_sectionp);
      SetEvent(consumer_progress_event);
      return kPglRetSuccess;
    }
    rlsp->consume_tail = read_iter;
    // We've processed all the consume-ready bytes...
    if (reterr != kPglRetSuccess) {
      // ...and we're at eof.
      LeaveCriticalSection(critical_sectionp);
      return reterr;
    }
    // ...and there's probably more.
    // Shouldn't be possible to get here twice.
    LeaveCriticalSection(critical_sectionp);
    SetEvent(consumer_progress_event);
    WaitForSingleObject(rlsp->reader_progress_event, INFINITE);
  }
#else
  pthread_mutex_t* sync_mutexp = &rlsp->sync_mutex;
  pthread_cond_t* consumer_progress_condvarp = &rlsp->consumer_progress_condvar;
  pthread_cond_t* reader_progress_condvarp = &rlsp->reader_progress_condvar;
  while (1) {
    pthread_mutex_lock(sync_mutexp);
    const PglErr reterr = rlsp->reterr;
    if (S_CAST(uint32_t, reterr) > S_CAST(uint32_t, kPglRetSkipped)) {
      pthread_mutex_unlock(sync_mutexp);
      return reterr;
    }
    char* available_end = rlsp->available_end;
    if (read_iter != available_end) {
      char* cur_circular_end = rlsp->cur_circular_end;
      if (cur_circular_end) {
        if (cur_circular_end == read_iter) {
          char* buf = rlsp->buf;
          read_iter = buf;
          rlsp->consume_tail = buf;
          rlsp->cur_circular_end = nullptr;
          pthread_cond_signal(consumer_progress_condvarp);
          pthread_mutex_unlock(sync_mutexp);
          return kPglRetSuccess;
        }
        rlsp->consume_stop = cur_circular_end;
      } else {
        rlsp->consume_stop = available_end;
      }
      rlsp->consume_tail = read_iter;
      pthread_cond_signal(consumer_progress_condvarp);
      pthread_mutex_unlock(sync_mutexp);
      return kPglRetSuccess;
    }
    rlsp->consume_tail = read_iter;
    // We've processed all the consume-ready bytes...
    if (reterr != kPglRetSuccess) {
      // ...and we're at eof.
      pthread_mutex_unlock(sync_mutexp);
      return reterr;
    }
    // ...and there's probably more.
    // Shouldn't be possible to get here twice.
    pthread_cond_signal(consumer_progress_condvarp);
    while (read_iter == rlsp->available_end) {
      // spurious wakeup guard
      pthread_cond_wait(reader_progress_condvarp, sync_mutexp);
    }
    pthread_mutex_unlock(sync_mutexp);
  }
#endif
}

PglErr RewindRLstream(ReadLineStream* rlsp, char** read_iterp) {
#ifdef _WIN32
#else
#endif
  ;;;
}

PglErr CleanupRLstream(ReadLineStream* rlsp) {
  PglErr reterr = kPglRetSuccess;
#ifdef _WIN32
  if (rlsp->read_thread) {
    if (!rlsp->at_eof) {
      rlsp->rewind_or_interrupt = 2;
    }
    CloseHandle(rlsp->consumer_progress_event);
    CloseHandle(rlsp->reader_progress_event);
  }
#else
  const uint32_t sync_init_state = rlsp->sync_init_state;
  if (sync_init_state) {
    if (sync_init_state == 4) {
      if (!rlsp->at_eof) {
        pthread_mutex_lock(&rlsp->sync_mutex);
        rlsp->rewind_or_interrupt = 2;

        while () {
          // Spurious wakeup guard.
          pthread_cond_wait(&rlsp->_, &rlsp->sync_mutex);
        }
        pthread_muteX_unlock(&rlsp->sync_mutex);
      }
    }
    pthread_mutex_destroy(&rlsp->sync_mutex);
    if (sync_init_state > 1) {
      pthread_cond_destroy(&rlsp->reader_progress_condvar);
      if (sync_init_state > 2) {
        pthread_cond_destroy(&rlsp->consumer_progress_condvar);
      }
    }
    rlsp->sync_init_state = 0;
  }
#endif
  if (rlsp->gz_infile) {
    if (gzclose_null(&rlsp->gz_infile)) {
      reterr = kPglRetReadFail;
    }
  }
  return reterr;
}
*/

#ifdef __cplusplus
}
#endif
