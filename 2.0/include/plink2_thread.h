#ifndef __PLINK2_THREAD_H__
#define __PLINK2_THREAD_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
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
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Basic multithreading code.  Uses native Win32 API instead of pthreads
// emulation on Windows.
#include "plink2_base.h"

#ifdef _WIN32
#  include <process.h>
#else
#  include <pthread.h>
#endif

// Most thread functions should be of the form
//   THREAD_FUNC_DECL function_name(void* raw_arg) {
//     ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
//     uint32_t tidx = arg->tidx;
//     ...
//     do {
//       ... // process current block
//     } while (!THREAD_BLOCK_FINISH(arg));
//     THREAD_RETURN;
//   }
// or
//   THREAD_FUNC_DECL function_name(void* raw_arg) {
//     ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
//     uint32_t tidx = arg->tidx;
//     ...
//     // parallelized initialization code, must wait for other threads to also
//     // finish this before proceeding to main loop
//     while (!THREAD_BLOCK_FINISH(arg)) {
//       ... // process current block
//     }
//     THREAD_RETURN;
//   }
#ifdef _WIN32
#  define pthread_t HANDLE
#  define THREAD_FUNC_DECL unsigned __stdcall
#  define THREAD_FUNCPTR_T(func_ptr) unsigned (__stdcall *func_ptr)(void*)
  // #define THREAD_FUNCPP_T(func_pp) unsigned (__stdcall **func_pp)(void*)
#  define THREAD_RETURN return 0
#else
#  define THREAD_FUNC_DECL void*
#  define THREAD_FUNCPTR_T(func_ptr) void* (*func_ptr)(void*)
  // #define THREAD_FUNCPP_T(func_pp) void* (**func_pp)(void*)
#  define THREAD_RETURN return nullptr
#endif

#if (__GNUC__ == 4) && (__GNUC_MINOR__ < 7) && !defined(__clang__)
// todo: check if this is also needed for any clang versions we care about.
// (support was added in clang 3.1, I think?)
#  define __ATOMIC_RELAXED 0
#  define __ATOMIC_CONSUME 1
#  define __ATOMIC_ACQUIRE 2
#  define __ATOMIC_RELEASE 3
#  define __ATOMIC_ACQ_REL 4
#  define __ATOMIC_SEQ_CST 5
#  define __atomic_fetch_add(ptr, val, memorder) __sync_fetch_and_add((ptr), (val))
#  define __atomic_fetch_sub(ptr, val, memorder) __sync_fetch_and_sub((ptr), (val))
#  define __atomic_sub_fetch(ptr, val, memorder) __sync_sub_and_fetch((ptr), (val))

HEADER_INLINE uint32_t ATOMIC_COMPARE_EXCHANGE_N_U32(uint32_t* ptr, uint32_t* expected, uint32_t desired, __maybe_unused int weak, __maybe_unused int success_memorder, __maybe_unused int failure_memorder) {
  const uint32_t new_expected = __sync_val_compare_and_swap(ptr, *expected, desired);
  if (new_expected == (*expected)) {
    return 1;
  }
  *expected = new_expected;
  return 0;
}

HEADER_INLINE uint32_t ATOMIC_COMPARE_EXCHANGE_N_U64(uint64_t* ptr, uint64_t* expected, uint64_t desired, __maybe_unused int weak, __maybe_unused int success_memorder, __maybe_unused int failure_memorder) {
  const uint64_t new_expected = __sync_val_compare_and_swap(ptr, *expected, desired);
  if (new_expected == (*expected)) {
    return 1;
  }
  *expected = new_expected;
  return 0;
}
#else
#  define ATOMIC_COMPARE_EXCHANGE_N_U32 __atomic_compare_exchange_n
#  define ATOMIC_COMPARE_EXCHANGE_N_U64 __atomic_compare_exchange_n
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

#ifdef _WIN32
void WaitForAllObjects(uint32_t ct, HANDLE* objs);
#endif

#ifdef _WIN32
// This should be increased once the old-style threads code has been purged.
CONSTI32(kMaxThreads, 64);
#else
// currently assumed to be less than 2^16 (otherwise some multiply overflows
// are theoretically possible, at least in the 32-bit build)
CONSTI32(kMaxThreads, 512);
#endif

#ifdef __APPLE__
// cblas_dgemm may fail with 128k
CONSTI32(kDefaultThreadStack, 524288);
#else
// asserts didn't seem to work properly with a setting much smaller than this
CONSTI32(kDefaultThreadStack, 131072);
#endif

typedef struct ThreadGroupControlBlockStruct {
  // Neither thread-functions nor the thread-group owner should touch these
  // variables directly.
  uintptr_t spawn_ct;
#ifdef _WIN32
  HANDLE start_next_events[2];  // uses parity of spawn_ct
  HANDLE cur_block_done_event;
#else
  pthread_mutex_t sync_mutex;
  pthread_cond_t cur_block_done_condvar;
  pthread_cond_t start_next_condvar;
#endif
  uint32_t active_ct;

  // Thread-functions can safely read from these.
  uint32_t thread_ct;

  // 1 = process last block and exit; 2 = immediate termination requested
  uint32_t is_last_block;
} ThreadGroupControlBlock;

typedef struct ThreadGroupSharedStruct {
  void* context;
#ifdef __cplusplus
  ThreadGroupControlBlock& GET_PRIVATE_cb() { return cb; }
  ThreadGroupControlBlock const& GET_PRIVATE_cb() const { return cb; }
 private:
#endif
  ThreadGroupControlBlock cb;
} ThreadGroupShared;

typedef struct ThreadGroupFuncArgStruct {
  ThreadGroupShared* sharedp;
  uint32_t tidx;
} ThreadGroupFuncArg;

typedef struct ThreadGroupMainStruct {
  ThreadGroupShared shared;
  THREAD_FUNCPTR_T(thread_func_ptr);
  pthread_t* threads;
  ThreadGroupFuncArg* thread_args;
  // Generally favor uint16_t/uint32_t over unsigned char/uint8_t for isolated
  // bools, since in the latter case the compiler is fairly likely to generate
  // worse code due to aliasing paranoia; see e.g.
  //   https://travisdowns.github.io/blog/2019/08/26/vector-inc.html
  uint16_t is_unjoined;
  uint16_t is_active;

#ifndef _WIN32
  uint32_t sync_init_state;
#endif
} ThreadGroupMain;

typedef struct ThreadGroupStruct {
#ifdef __cplusplus
  ThreadGroupMain& GET_PRIVATE_m() { return m; }
  ThreadGroupMain const& GET_PRIVATE_m() const { return m; }
 private:
#endif
  ThreadGroupMain m;
} ThreadGroup;

void PreinitThreads(ThreadGroup* tg_ptr);

// Return value is clipped to 1..kMaxThreads.
// If known_procs_ptr is non-null, it's set to the raw unclipped value (which
// can theoretically be -1 if the sysconf call fails)
uint32_t NumCpu(int32_t* known_procs_ptr);

// Also allocates, returning 1 on failure.
BoolErr SetThreadCt(uint32_t thread_ct, ThreadGroup* tg_ptr);

HEADER_INLINE uint32_t GetThreadCt(const ThreadGroupShared* sharedp) {
  return GET_PRIVATE(*sharedp, cb).thread_ct;
}

HEADER_INLINE uint32_t GetThreadCtTg(const ThreadGroup* tg_ptr) {
  const ThreadGroupMain* tgp = &GET_PRIVATE(*tg_ptr, m);
  return GET_PRIVATE(tgp->shared, cb).thread_ct;
}

HEADER_INLINE void SetThreadFuncAndData(THREAD_FUNCPTR_T(start_routine), void* shared_context, ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = &GET_PRIVATE(*tg_ptr, m);
  assert(!tgp->is_active);
  tgp->shared.context = shared_context;
  GET_PRIVATE(tgp->shared, cb).is_last_block = 0;
  tgp->thread_func_ptr = start_routine;
}

// Equivalent to SetThreadFuncAndData() with unchanged
// start_routine/shared_context.  Ok to call this "unnecessarily".
HEADER_INLINE void ReinitThreads(ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = &GET_PRIVATE(*tg_ptr, m);
  assert(!tgp->is_active);
  GET_PRIVATE(tgp->shared, cb).is_last_block = 0;
}

HEADER_INLINE uint32_t ThreadsAreActive(ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = &GET_PRIVATE(*tg_ptr, m);
  return tgp->is_active;
}

// Technically unnecessary to call this, but it does save one sync cycle.
//
// Note that, if there's only one block of work-shards, this should be called
// before the first SpawnThreads() call.
HEADER_INLINE void DeclareLastThreadBlock(ThreadGroup* tg_ptr) {
  ThreadGroupMain* tgp = &GET_PRIVATE(*tg_ptr, m);
  assert(!tgp->is_unjoined);
  GET_PRIVATE(tgp->shared, cb).is_last_block = 1;
}

HEADER_INLINE uint32_t IsLastBlock(const ThreadGroup* tg_ptr) {
  const ThreadGroupMain* tgp = &GET_PRIVATE(*tg_ptr, m);
  return GET_PRIVATE(tgp->shared, cb).is_last_block;
}

#if defined(__cplusplus) && !defined(_WIN32)
class Plink2ThreadStartup {
public:
  pthread_attr_t smallstack_thread_attr;
  Plink2ThreadStartup() {
#  ifdef NDEBUG
    // we'll error out for another reason soon enough if there's insufficient
    // memory...
    pthread_attr_init(&smallstack_thread_attr);
#  else
    assert(!pthread_attr_init(&smallstack_thread_attr));
#  endif
    // if this fails due to kDefaultThreadStack being smaller than the system
    // page size, no need to error out
    pthread_attr_setstacksize(&smallstack_thread_attr, kDefaultThreadStack);
  }

  ~Plink2ThreadStartup() {
    pthread_attr_destroy(&smallstack_thread_attr);
  }
};

extern Plink2ThreadStartup g_thread_startup;
#endif

BoolErr SpawnThreads(ThreadGroup* tg_ptr);

void JoinThreads(ThreadGroup* tg_ptr);

// Assumes threads are joined.
void StopThreads(ThreadGroup* tg_ptr);

void CleanupThreads(ThreadGroup* tg_ptr);

#ifdef _WIN32
HEADER_INLINE BoolErr THREAD_BLOCK_FINISH(ThreadGroupFuncArg* tgfap) {
  ThreadGroupControlBlock* cbp = &(GET_PRIVATE(*tgfap->sharedp, cb));
  if (cbp->is_last_block) {
    return 1;
  }
  const uint32_t start_next_parity = cbp->spawn_ct & 1;
  if (!__atomic_sub_fetch(&cbp->active_ct, 1, __ATOMIC_ACQ_REL)) {
    SetEvent(cbp->cur_block_done_event);
  }
  WaitForSingleObject(cbp->start_next_events[start_next_parity], INFINITE);
  return (cbp->is_last_block == 2);
}
#else
BoolErr THREAD_BLOCK_FINISH(ThreadGroupFuncArg* tgfap);
#endif

// Convenience functions for potentially-small-and-frequent jobs where
// thread_ct == 0 corresponds to not wanting to launch threads at all; see
// MakeDupflagHtable in plink2_cmdline for a typical use case.
HEADER_INLINE BoolErr SetThreadCt0(uint32_t thread_ct, ThreadGroup* tg_ptr) {
  if (!thread_ct) {
    return 0;
  }
  return SetThreadCt(thread_ct, tg_ptr);
}

HEADER_INLINE void JoinThreads0(ThreadGroup* tg_ptr) {
  if (GET_PRIVATE(*tg_ptr, m).threads) {
    JoinThreads(tg_ptr);
  }
}

// This comes in handy a lot in multithreaded error-reporting code when
// deterministic behavior is desired.
void UpdateU64IfSmaller(uint64_t newval, uint64_t* oldval_ptr);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_THREAD_H__
