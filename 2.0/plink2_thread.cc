// This library is part of PLINK 2.00, copyright (C) 2005-2019 Shaun Purcell,
// Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.


#include "plink2_thread.h"

#ifndef _WIN32
#  include <unistd.h>  // sysconf()
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

void PreinitThreads(ThreadGroup* tgp) {
  tgp->shared.cb.is_last_block = 0;
  tgp->thread_func_ptr = nullptr;
  tgp->threads = nullptr;
  tgp->is_unjoined = 0;
  tgp->is_active = 0;
}

uint32_t NumCpu(int32_t* known_procs_ptr) {
#ifdef _WIN32
  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  const int32_t known_procs = sysinfo.dwNumberOfProcessors;
  uint32_t max_thread_ct = known_procs;
#else
  const int32_t known_procs = sysconf(_SC_NPROCESSORS_ONLN);
  uint32_t max_thread_ct = (known_procs == -1)? 1 : known_procs;
#endif
  if (known_procs_ptr) {
    *known_procs_ptr = known_procs;
  }
  if (max_thread_ct > kMaxThreads) {
    max_thread_ct = kMaxThreads;
  }
  return max_thread_ct;
}

BoolErr SetThreadCt(uint32_t thread_ct, ThreadGroup* tgp) {
  assert(!tgp->is_active);
  if (tgp->threads) {
    free(tgp->threads);
    tgp->threads = nullptr;
  }
  assert(thread_ct && (thread_ct <= kMaxThreads));
#ifdef _WIN32
  unsigned char* memptr = S_CAST(unsigned char*, malloc(thread_ct * (sizeof(ThreadGroupFuncArg) + 3 * sizeof(HANDLE))));
  if (unlikely(!memptr)) {
    return 1;
  }
  tgp->threads = R_CAST(HANDLE*, memptr);
  tgp->shared.cb.start_next_events = &(tgp->threads[thread_ct]);
  tgp->shared.cb.cur_block_done_events = &(tgp->shared.cb.start_next_events[thread_ct]);
  memset(tgp->threads, 0, thread_ct * (3 * sizeof(HANDLE)));
  memptr = &(memptr[thread_ct * (3 * sizeof(HANDLE))]);
#else
  unsigned char* memptr = S_CAST(unsigned char*, malloc(thread_ct * (sizeof(pthread_t) + sizeof(ThreadGroupFuncArg))));
  if (unlikely(!memptr)) {
    return 1;
  }
  tgp->threads = R_CAST(pthread_t*, memptr);
  // !is_active currently guarantees that the sync mutex/condvars are not
  // initialized.  Could change this later.
  tgp->shared.cb.active_ct = 0;
  tgp->sync_init_bits = 0;
  memptr = &(memptr[thread_ct * sizeof(pthread_t)]);
#endif
  tgp->thread_args = R_CAST(ThreadGroupFuncArg*, memptr);

  tgp->shared.cb.thread_ct = thread_ct;
  return 0;
}

// Note that thread_ct is permitted to be less than tgp->shared.cb.thread_ct,
// to support the SpawnThreads() error cases.
void JoinThreadsInternal(uint32_t thread_ct, ThreadGroup* tgp) {
#ifdef _WIN32
  if (!tgp->shared.cb.is_last_block) {
    WaitForMultipleObjects(thread_ct, tgp->shared.cb.cur_block_done_events, 1, INFINITE);
  } else {
    WaitForMultipleObjects(thread_ct, tgp->threads, 1, INFINITE);
    for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
      CloseHandle(tgp->threads[tidx]);
    }
    const uint32_t orig_thread_ct = tgp->shared.cb.thread_ct;
    for (uint32_t tidx = 0; tidx != orig_thread_ct; ++tidx) {
      CloseHandle(tgp->shared.cb.start_next_events[tidx]);
      CloseHandle(tgp->shared.cb.cur_block_done_events[tidx]);
    }
    memset(tgp->threads, 0, orig_thread_ct * (3 * sizeof(HANDLE)));
    tgp->is_active = 0;
  }
#else
  if (!tgp->shared.cb.is_last_block) {
    pthread_mutex_lock(&tgp->shared.cb.sync_mutex);
    while (tgp->shared.cb.active_ct) {
      pthread_cond_wait(&tgp->shared.cb.cur_block_done_condvar, &tgp->shared.cb.sync_mutex);
    }
    // keep mutex until next block loaded
  } else {
    for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
      pthread_join(tgp->threads[tidx], nullptr);
    }
    pthread_mutex_destroy(&tgp->shared.cb.sync_mutex);
    pthread_cond_destroy(&tgp->shared.cb.cur_block_done_condvar);
    pthread_cond_destroy(&tgp->shared.cb.start_next_condvar);
    tgp->is_active = 0;
    tgp->sync_init_bits = 0;
  }
#endif
  tgp->is_unjoined = 0;
}

#if defined(__cplusplus) && !defined(_WIN32)
Plink2ThreadStartup g_thread_startup;
#endif

BoolErr SpawnThreads(ThreadGroup* tgp) {
  ThreadGroupControlBlock* cbp = &(tgp->shared.cb);
  const uint32_t thread_ct = cbp->thread_ct;
  const uint32_t was_active = tgp->is_active;
  const uint32_t is_last_block = cbp->is_last_block;
  pthread_t* threads = tgp->threads;
  assert(threads != nullptr);
#ifdef _WIN32
  if (!was_active) {
    cbp->spawn_ct = 0;
    for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
      HANDLE cur_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
      if (unlikely(!cur_event)) {
        return 1;
      }
      cbp->start_next_events[tidx] = cur_event;
      cur_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
      if (unlikely(!cur_event)) {
        return 1;
      }
      cbp->cur_block_done_events[tidx] = cur_event;
    }
    for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
      ThreadGroupFuncArg* arg_slot = &(tgp->thread_args[tidx]);
      arg_slot->sharedp = &(tgp->shared);
      arg_slot->tidx = tidx;
      threads[tidx] = R_CAST(HANDLE, _beginthreadex(nullptr, kDefaultThreadStack, tgp->thread_func_ptr, arg_slot, 0, nullptr));
      if (unlikely(!threads[tidx])) {
        if (tidx) {
          JoinThreadsInternal(tidx, tgp);
        }
        if (!is_last_block) {
          for (uint32_t tidx2 = 0; tidx2 != tidx; ++tidx2) {
            TerminateThread(threads[tidx2], 0);
          }
          for (uint32_t tidx2 = 0; tidx2 != tidx; ++tidx2) {
            CloseHandle(threads[tidx2]);
            threads[tidx2] = nullptr;
          }
        }
        return 1;
      }
    }
    tgp->is_active = 1;
  } else {
    cbp->spawn_ct += 1;
    for (uintptr_t tidx = 0; tidx != thread_ct; ++tidx) {
      SetEvent(cbp->start_next_events[tidx]);
    }
  }
#else
  if (!is_last_block) {
    cbp->active_ct = thread_ct;
  }
  if (!was_active) {
    cbp->spawn_ct = 0;
    assert(!tgp->sync_init_bits);
    if (unlikely(pthread_mutex_init(&cbp->sync_mutex, nullptr))) {
      return 1;
    }
    tgp->sync_init_bits = 1;
    if (unlikely(pthread_cond_init(&cbp->cur_block_done_condvar, nullptr))) {
      return 1;
    }
    tgp->sync_init_bits |= 2;
    if (unlikely(pthread_cond_init(&cbp->start_next_condvar, nullptr))) {
      return 1;
    }
    tgp->sync_init_bits |= 4;
#  ifndef __cplusplus
    pthread_attr_t smallstack_thread_attr;
    if (unlikely(pthread_attr_init(&smallstack_thread_attr))) {
      return 1;
    }
    pthread_attr_setstacksize(&smallstack_thread_attr, kDefaultThreadStack);
#  endif
    for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
      ThreadGroupFuncArg* arg_slot = &(tgp->thread_args[tidx]);
      arg_slot->sharedp = &(tgp->shared);
      arg_slot->tidx = tidx;
      if (unlikely(pthread_create(&(threads[tidx]),
#  ifdef __cplusplus
                                  &g_thread_startup.smallstack_thread_attr,
#  else
                                  &smallstack_thread_attr,
#  endif
                                  tgp->thread_func_ptr, arg_slot))) {
        if (tidx) {
          if (is_last_block) {
            JoinThreadsInternal(tidx, tgp);
          } else {
            const uint32_t unstarted_thread_ct = thread_ct - tidx;
            pthread_mutex_lock(&cbp->sync_mutex);
            cbp->active_ct -= unstarted_thread_ct;
            while (cbp->active_ct) {
              pthread_cond_wait(&cbp->cur_block_done_condvar, &cbp->sync_mutex);
            }
            // not worth the trouble of demanding that all callers handle
            // pthread_create() failure cleanly
            for (uint32_t tidx2 = 0; tidx2 != tidx; ++tidx2) {
              pthread_cancel(threads[tidx2]);
            }
          }
        } else {
          cbp->active_ct = 0;
        }
#  ifndef __cplusplus
        pthread_attr_destroy(&smallstack_thread_attr);
#  endif
        return 1;
      }
    }
#  ifndef __cplusplus
    pthread_attr_destroy(&smallstack_thread_attr);
#  endif
    tgp->is_active = 1;
  } else {
    cbp->spawn_ct += 1;
    // still holding mutex
    pthread_mutex_unlock(&cbp->sync_mutex);
    pthread_cond_broadcast(&cbp->start_next_condvar);
  }
#endif
  tgp->is_unjoined = 1;
  return 0;
}

void JoinThreads(ThreadGroup* tgp) {
  JoinThreadsInternal(tgp->shared.cb.thread_ct, tgp);
}

void CleanupThreads(ThreadGroup* tgp) {
  if (tgp->threads) {
    ThreadGroupControlBlock* cbp = &(tgp->shared.cb);
    const uint32_t thread_ct = cbp->thread_ct;
    if (tgp->is_active) {
      if (tgp->is_unjoined) {
        JoinThreadsInternal(thread_ct, tgp);
      }
      if (!cbp->is_last_block) {
        cbp->is_last_block = 2;
        SpawnThreads(tgp);
        JoinThreadsInternal(thread_ct, tgp);
      }
    } else {
#ifdef _WIN32
      for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
        if (cbp->start_next_events[tidx]) {
          CloseHandle(cbp->start_next_events[tidx]);
        }
        if (cbp->cur_block_done_events[tidx]) {
          CloseHandle(cbp->cur_block_done_events[tidx]);
        }
      }
#else
      const uint32_t sync_init_bits = tgp->sync_init_bits;
      if (sync_init_bits & 1) {
        pthread_mutex_destroy(&cbp->sync_mutex);
        if (sync_init_bits & 2) {
          pthread_cond_destroy(&cbp->cur_block_done_condvar);
          if (sync_init_bits & 4) {
            pthread_cond_destroy(&cbp->start_next_condvar);
          }
        }
        tgp->sync_init_bits = 0;
      }
#endif
    }
#ifndef _WIN32
    assert(!tgp->shared.cb.active_ct);
#endif
    cbp->thread_ct = 0;
    free(tgp->threads);
    tgp->threads = nullptr;
  }
  tgp->shared.cb.is_last_block = 0;
  tgp->thread_func_ptr = nullptr;
}

#ifndef _WIN32
BoolErr THREAD_BLOCK_FINISH(ThreadGroupFuncArg* tgfap) {
  ThreadGroupControlBlock* cbp = &(tgfap->sharedp->cb);
  if (cbp->is_last_block) {
    return 1;
  }
  const uintptr_t initial_spawn_ct = cbp->spawn_ct;
  pthread_mutex_lock(&cbp->sync_mutex);
  if (!(--cbp->active_ct)) {
    pthread_cond_signal(&cbp->cur_block_done_condvar);
  }
  while (cbp->spawn_ct == initial_spawn_ct) {
    // spurious wakeup guard
    pthread_cond_wait(&cbp->start_next_condvar, &cbp->sync_mutex);
  }
  pthread_mutex_unlock(&cbp->sync_mutex);
  return (cbp->is_last_block == 2);
}
#endif

#ifdef __cplusplus
}  // namespace plink2
#endif
