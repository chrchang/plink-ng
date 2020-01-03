// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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

static inline ThreadGroupMain* GetTgp(ThreadGroup* tg_ptr) {
  return &GET_PRIVATE(*tg_ptr, m);
}

static inline ThreadGroupControlBlock* GetCbp(ThreadGroupShared* sharedp) {
  return &GET_PRIVATE(*sharedp, cb);
}

#ifdef _WIN32
void WaitForAllObjects(uint32_t ct, HANDLE* objs) {
  while (ct > 64) {
    WaitForMultipleObjects(64, objs, 1, INFINITE);
    objs = &(objs[64]);
    ct -= 64;
  }
  WaitForMultipleObjects(ct, objs, 1, INFINITE);
}
#endif

void PreinitThreads(ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = GetTgp(tg_ptr);
  GetCbp(&tgp->shared)->is_last_block = 0;
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

BoolErr SetThreadCt(uint32_t thread_ct, ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = GetTgp(tg_ptr);
  assert(!tgp->is_active);
  if (tgp->threads) {
    free(tgp->threads);
    tgp->threads = nullptr;
  }
  assert(thread_ct && (thread_ct <= kMaxThreads));
#ifdef _WIN32
  unsigned char* memptr = S_CAST(unsigned char*, malloc(thread_ct * (sizeof(ThreadGroupFuncArg) + sizeof(HANDLE))));
  if (unlikely(!memptr)) {
    return 1;
  }
  tgp->threads = R_CAST(HANDLE*, memptr);
  memset(tgp->threads, 0, thread_ct * sizeof(HANDLE));
  memptr = &(memptr[thread_ct * sizeof(HANDLE)]);
#else
  unsigned char* memptr = S_CAST(unsigned char*, malloc(thread_ct * (sizeof(pthread_t) + sizeof(ThreadGroupFuncArg))));
  if (unlikely(!memptr)) {
    return 1;
  }
  tgp->threads = R_CAST(pthread_t*, memptr);
  // !is_active currently guarantees that the sync events/mutex/condvars are
  // not initialized.  Could change this later.
  tgp->sync_init_state = 0;
  memptr = &(memptr[thread_ct * sizeof(pthread_t)]);
#endif
  ThreadGroupControlBlock* cbp = GetCbp(&tgp->shared);
  cbp->active_ct = 0;
  tgp->thread_args = R_CAST(ThreadGroupFuncArg*, memptr);

  cbp->thread_ct = thread_ct;
  return 0;
}

// Note that thread_ct is permitted to be less than tgp->shared.cb.thread_ct,
// to support the SpawnThreads() error cases.
void JoinThreadsInternal(uint32_t thread_ct, ThreadGroupMain* tgp) {
  assert(tgp->is_active);
  ThreadGroupControlBlock* cbp = GetCbp(&tgp->shared);
#ifdef _WIN32
  if (!cbp->is_last_block) {
    WaitForSingleObject(cbp->cur_block_done_event, INFINITE);
  } else {
    WaitForAllObjects(thread_ct, tgp->threads);
    for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
      CloseHandle(tgp->threads[tidx]);
    }
    CloseHandle(cbp->start_next_events[0]);
    CloseHandle(cbp->start_next_events[1]);
    CloseHandle(cbp->cur_block_done_event);
    memset(tgp->threads, 0, thread_ct * sizeof(HANDLE));
    tgp->is_active = 0;
  }
#else
  if (!cbp->is_last_block) {
    pthread_mutex_lock(&cbp->sync_mutex);
    while (cbp->active_ct) {
      pthread_cond_wait(&cbp->cur_block_done_condvar, &cbp->sync_mutex);
    }
    // keep mutex until next block loaded
  } else {
    for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
      pthread_join(tgp->threads[tidx], nullptr);
    }
    pthread_mutex_destroy(&cbp->sync_mutex);
    pthread_cond_destroy(&cbp->cur_block_done_condvar);
    pthread_cond_destroy(&cbp->start_next_condvar);
    tgp->is_active = 0;
    tgp->sync_init_state = 0;
  }
#endif
  tgp->is_unjoined = 0;
}

#if defined(__cplusplus) && !defined(_WIN32)
Plink2ThreadStartup g_thread_startup;
#endif

BoolErr SpawnThreads(ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = GetTgp(tg_ptr);
  ThreadGroupControlBlock* cbp = GetCbp(&tgp->shared);
  const uint32_t thread_ct = cbp->thread_ct;
  const uint32_t was_active = tgp->is_active;
  const uint32_t is_last_block = cbp->is_last_block;
  pthread_t* threads = tgp->threads;
  assert(threads != nullptr);
#ifdef _WIN32
  if (!was_active) {
    cbp->spawn_ct = 0;
    // manual-reset broadcast event
    HANDLE cur_event = CreateEvent(nullptr, TRUE, FALSE, nullptr);
    if (unlikely(!cur_event)) {
      return 1;
    }
    cbp->start_next_events[0] = cur_event;
    cur_event = CreateEvent(nullptr, TRUE, FALSE, nullptr);
    if (unlikely(!cur_event)) {
      CloseHandle(cbp->start_next_events[0]);
      return 1;
    }
    cbp->start_next_events[1] = cur_event;
    cur_event = CreateEvent(nullptr, FALSE, FALSE, nullptr);
    if (unlikely(!cur_event)) {
      CloseHandle(cbp->start_next_events[0]);
      CloseHandle(cbp->start_next_events[1]);
      return 1;
    }
    cbp->cur_block_done_event = cur_event;
    cbp->active_ct = thread_ct;
    for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
      ThreadGroupFuncArg* arg_slot = &(tgp->thread_args[tidx]);
      arg_slot->sharedp = &(tgp->shared);
      arg_slot->tidx = tidx;
      threads[tidx] = R_CAST(HANDLE, _beginthreadex(nullptr, kDefaultThreadStack, tgp->thread_func_ptr, arg_slot, 0, nullptr));
      if (unlikely(!threads[tidx])) {
        if (tidx) {
          if (!is_last_block) {
            // not sure the old TerminateThread() code ever worked properly in
            // the first place...
            // anyway, new contract makes clean error-shutdown easy
            if (!__atomic_sub_fetch(&cbp->active_ct, thread_ct - tidx, __ATOMIC_ACQ_REL)) {
              SetEvent(cbp->cur_block_done_event);
            }
            JoinThreadsInternal(tidx, tgp);
            cbp->is_last_block = 2;
            const uint32_t start_next_parity = cbp->spawn_ct & 1;
            // no need to reset start_prev_parity
            cbp->spawn_ct += 1;
            SetEvent(cbp->start_next_events[start_next_parity]);
          }
          JoinThreadsInternal(tidx, tgp);
        }
        return 1;
      }
    }
    tgp->is_active = 1;
  } else {
    cbp->spawn_ct += 1;
    cbp->active_ct = thread_ct;
    const uint32_t start_prev_parity = cbp->spawn_ct & 1;
    ResetEvent(cbp->start_next_events[start_prev_parity]);
    SetEvent(cbp->start_next_events[start_prev_parity ^ 1]);
  }
#else
  if (!is_last_block) {
    cbp->active_ct = thread_ct;
  }
  if (!was_active) {
    cbp->spawn_ct = 0;
    assert(!tgp->sync_init_state);
    if (unlikely(pthread_mutex_init(&cbp->sync_mutex, nullptr))) {
      return 1;
    }
    if (unlikely(pthread_cond_init(&cbp->cur_block_done_condvar, nullptr))) {
      tgp->sync_init_state = 1;
      return 1;
    }
    if (unlikely(pthread_cond_init(&cbp->start_next_condvar, nullptr))) {
      tgp->sync_init_state = 2;
      return 1;
    }
    tgp->sync_init_state = 3;
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
          if (!is_last_block) {
            JoinThreadsInternal(tidx, tgp);
            const uint32_t unstarted_thread_ct = thread_ct - tidx;
            cbp->active_ct -= unstarted_thread_ct;
            while (cbp->active_ct) {
              pthread_cond_wait(&cbp->cur_block_done_condvar, &cbp->sync_mutex);
            }
            cbp->is_last_block = 2;
            cbp->spawn_ct += 1;
            pthread_cond_broadcast(&cbp->start_next_condvar);
            pthread_mutex_unlock(&cbp->sync_mutex);
          }
          JoinThreadsInternal(tidx, tgp);
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
    pthread_cond_broadcast(&cbp->start_next_condvar);
    pthread_mutex_unlock(&cbp->sync_mutex);
  }
#endif
  tgp->is_unjoined = 1;
  return 0;
}

void JoinThreads(ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = GetTgp(tg_ptr);
  JoinThreadsInternal(GetCbp(&tgp->shared)->thread_ct, tgp);
}

void StopThreads(ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = GetTgp(tg_ptr);
  ThreadGroupControlBlock* cbp = GetCbp(&tgp->shared);
  cbp->is_last_block = 2;
  SpawnThreads(tg_ptr);
  JoinThreadsInternal(cbp->thread_ct, tgp);
}

void CleanupThreads(ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = GetTgp(tg_ptr);
  ThreadGroupControlBlock* cbp = GetCbp(&tgp->shared);
  if (tgp->threads) {
    const uint32_t thread_ct = cbp->thread_ct;
    if (tgp->is_active) {
      if (tgp->is_unjoined) {
        JoinThreadsInternal(thread_ct, tgp);
      }
      if (!cbp->is_last_block) {
        StopThreads(tg_ptr);
      }
#ifndef _WIN32
    } else {
      do {
        const uint32_t sync_init_state = tgp->sync_init_state;
        if (!sync_init_state) {
          break;
        }
        pthread_mutex_destroy(&cbp->sync_mutex);
        if (sync_init_state == 1) {
          break;
        }
        pthread_cond_destroy(&cbp->cur_block_done_condvar);
        if (sync_init_state == 2) {
          break;
        }
        pthread_cond_destroy(&cbp->start_next_condvar);
      } while (0);
      tgp->sync_init_state = 0;
#endif
    }
#ifndef _WIN32
    assert(!cbp->active_ct);
#endif
    cbp->thread_ct = 0;
    free(tgp->threads);
    tgp->threads = nullptr;
  }
  cbp->is_last_block = 0;
  tgp->thread_func_ptr = nullptr;
}

#ifndef _WIN32
BoolErr THREAD_BLOCK_FINISH(ThreadGroupFuncArg* tgfap) {
  ThreadGroupControlBlock* cbp = GetCbp(tgfap->sharedp);
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

void UpdateU64IfSmaller(uint64_t newval, uint64_t* oldval_ptr) {
  uint64_t oldval = *oldval_ptr;
  while (oldval > newval) {
    if (ATOMIC_COMPARE_EXCHANGE_N_U64(oldval_ptr, &oldval, newval, 1, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {
      break;
    }
  }
}

#ifdef __cplusplus
}  // namespace plink2
#endif
