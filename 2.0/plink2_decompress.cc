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
  ReadStreamThread* context = R_CAST(ReadStreamThread*, arg);
  while (1) {
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
    rlsp->line_available_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (!rlsp->line_available_event) {
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->consume_done_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (!rlsp->consume_done_event) {
      CloseHandle(rlsp->line_available_event);
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->read_thread = R_CAST(HANDLE, _beginthreadex(nullptr, kDefaultThreadStack, ReadLineStreamThread, rlsp, 0, nullptr));
    if (!rlsp->read_thread) {
      CloseHandle(rlsp->consume_done_event);
      CloseHandle(rlsp->line_available_event);
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
#else
    if (pthread_mutex_init(&rlsp->sync_mutex, nullptr)) {
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 1;
    if (pthread_cond_init(&rlsp->line_available_condvar, nullptr)) {
      goto InitRLstream_ret_THREAD_CREATE_FAIL;
    }
    rlsp->sync_init_state = 2;
    if (pthread_cond_init(&rlsp->consume_done_condvar, nullptr)) {
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

PglErr RewindRLstream(ReadLineStream* rlsp) {
}

PglErr CleanupRLstream(ReadLineStream* rlsp) {
  PglErr reterr = kPglRetSuccess;
#ifdef _WIN32
  if (rlsp->read_thread) {
    if (!rlsp->at_eof) {
      rlsp->rewind_or_interrupt = 2;
    }
    CloseHandle(rlsp->consume_done_event);
    CloseHandle(rlsp->line_available_event);
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
      pthread_cond_destroy(&rlsp->line_available_condvar);
      if (sync_init_state > 2) {
        pthread_cond_destroy(&rlsp->consume_done_condvar);
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
