// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang, Benjamin Demaille.
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

#include "plink2_htable.h"

#include <errno.h>
#include <string.h>

#include "plink2_bits.h"
#include "plink2_memory.h"
#include "plink2_string.h"
#include "plink2_text.h"
#include "plink2_thread.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Hash32() is a C99 port of absl::Hash<std::string_view>'s portable 64-bit
// code path (the one used when neither CRC32 nor AES hardware acceleration is
// enabled), from Abseil LTS 20260817.0
// (absl/hash/internal/hash.h and hash.cc), followed by extraction of the top
// 32 bits.  Abseil is copyright 2018 The Abseil Authors, and licensed under
// the Apache License, Version 2.0:
//   https://www.apache.org/licenses/LICENSE-2.0
// Modifications from the original:
// - Translated to C99.
// - absl::Hash's per-process seed (the address of a global variable) is
//   replaced with a fixed seed of 0, so hash values are reproducible across
//   runs and builds.
// - The CRC32C (<= 32 bytes) and AES (> 32 bytes) code paths, and the 32-bit
//   CityHash path, are not included; all builds compute the same function.
//
// HashCombineContiguous(0, key, len) matches
// absl::hash_internal::MixingHashState::hash_with_seed(string_view, 0)
// exactly, as computed by an Abseil build without hardware acceleration.
//
// Seed 0 benchmarked slightly faster than a dense 64-bit seed (~1-2ns off
// builds and misses for 10M-entry tables on an Apple M2 Pro) and loses
// nothing in hash quality: <= 32-byte inputs start from
// kHashStaticRandomData-derived state regardless of seed, and longer inputs
// start from kHashStaticRandomData[0] ^ len.

static const uint64_t kHashMul = 0x79d5f9e0de1e8cf5LLU;

// Hexadecimal digits of pi's fractional part.
static const uint64_t kHashStaticRandomData[5] = {
  0x243f6a8885a308d3LLU, 0x13198a2e03707344LLU, 0xa4093822299f31d0LLU,
  0x082efa98ec4e6c89LLU, 0x452821e638d01377LLU
};

static inline uint64_t HashMix(uint64_t lhs, uint64_t rhs) {
  uint64_t high;
  const uint64_t low = multiply64to128(lhs, rhs, &high);
  return high ^ low;
}

static inline uint64_t HashRead8(const unsigned char* p) {
  uint64_t v;
  memcpy_k(&v, p, 8);
  return v;
}

static inline uint32_t HashRead4(const unsigned char* p) {
  uint32_t v;
  memcpy_k(&v, p, 4);
  return v;
}

// len must be <= 32.
static inline uint64_t HashPrecombineLengthMix(uint64_t state, uintptr_t len) {
  return state ^ HashRead8(&(R_CAST(const unsigned char*, kHashStaticRandomData)[len]));
}

static inline uint64_t HashMix32Bytes(const unsigned char* p, uint64_t state) {
  const uint64_t a = HashRead8(p);
  const uint64_t b = HashRead8(&(p[8]));
  const uint64_t c = HashRead8(&(p[16]));
  const uint64_t d = HashRead8(&(p[24]));
  const uint64_t cs0 = HashMix(a ^ kHashStaticRandomData[1], b ^ state);
  const uint64_t cs1 = HashMix(c ^ kHashStaticRandomData[2], d ^ state);
  return cs0 ^ cs1;
}

// len must be in (32, 1024].
static uint64_t HashLenGt32(uint64_t seed, const unsigned char* p, uintptr_t len) {
  uint64_t state = seed ^ kHashStaticRandomData[0] ^ len;
  const unsigned char* last32 = &(p[len - 32]);
  if (len <= 64) {
    return HashMix32Bytes(last32, HashMix32Bytes(p, state));
  }
  // Four independent dependency chains.
  uint64_t state1 = state;
  uint64_t state2 = state;
  uint64_t state3 = state;
  do {
    const uint64_t a = HashRead8(p);
    const uint64_t b = HashRead8(&(p[8]));
    const uint64_t c = HashRead8(&(p[16]));
    const uint64_t d = HashRead8(&(p[24]));
    const uint64_t e = HashRead8(&(p[32]));
    const uint64_t f = HashRead8(&(p[40]));
    const uint64_t g = HashRead8(&(p[48]));
    const uint64_t h = HashRead8(&(p[56]));
    state = HashMix(a ^ kHashStaticRandomData[1], b ^ state);
    state1 = HashMix(c ^ kHashStaticRandomData[2], d ^ state1);
    state2 = HashMix(e ^ kHashStaticRandomData[3], f ^ state2);
    state3 = HashMix(g ^ kHashStaticRandomData[4], h ^ state3);
    p = &(p[64]);
    len -= 64;
  } while (len > 64);
  state = (state ^ state1) ^ (state2 + state3);
  if (len > 32) {
    state = HashMix32Bytes(p, state);
  }
  return HashMix32Bytes(last32, state);
}

static uint64_t HashCombineContiguous(uint64_t state, const unsigned char* p, uintptr_t len) {
  if (len <= 8) {
    state = HashPrecombineLengthMix(state, len);
    uint64_t v;
    if (len >= 4) {
      // e.g. "ABCDE" is read as "BCDEABCD" on little-endian systems.
      v = (S_CAST(uint64_t, HashRead4(p)) << 32) | HashRead4(&(p[len - 4]));
    } else if (len) {
      // "A" is read as "AAA", "AB" as "ABB".
      v = (S_CAST(uint32_t, p[0]) << 16) | (S_CAST(uint32_t, p[len / 2]) << 8) | p[len - 1];
    } else {
      v = 0x57;
    }
    return HashMix(state ^ v, kHashMul);
  }
  if (len <= 16) {
    state = HashPrecombineLengthMix(state, len);
    return HashMix(state ^ HashRead8(p), kHashMul ^ HashRead8(&(p[len - 8])));
  }
  if (len <= 32) {
    // Two overlapping 16-byte halves, mixed in parallel.
    state = HashPrecombineLengthMix(state, len);
    const unsigned char* tail16 = &(p[len - 16]);
    const uint64_t m0 = HashMix(HashRead8(p) ^ kHashStaticRandomData[1], HashRead8(&(p[8])) ^ state);
    const uint64_t m1 = HashMix(HashRead8(tail16) ^ kHashStaticRandomData[3], HashRead8(&(tail16[8])) ^ state);
    return m0 ^ m1;
  }
  // absl processes inputs longer than 1024 bytes in 1024-byte pieces; this is
  // retained so that this is an exact port.
  while (len > 1024) {
    state = HashLenGt32(state, p, 1024);
    p = &(p[1024]);
    len -= 1024;
    if (len <= 32) {
      return HashCombineContiguous(state, p, len);
    }
  }
  return HashLenGt32(state, p, len);
}

uint32_t Hash32(const void* key, uint32_t len) {
  // Top 32 bits, since Hashceil() and other callers use high bits.
  return HashCombineContiguous(0, S_CAST(const unsigned char*, key), len) >> 32;
}

uint32_t IdHtableFindNoOverread(const char* cur_id, const char* const* item_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size) {
  // returns UINT32_MAX on failure
  for (uint32_t hashval = Hashceil(cur_id, cur_id_slen, id_htable_size); ; ) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if (cur_htable_idval == UINT32_MAX) {
      return UINT32_MAX;
    }
    // bit 31 may be set if the table was constructed with a dup-tolerant
    // function; look up/return the first-seen item in that case
    const uint32_t cur_item_idx = cur_htable_idval & 0x7fffffff;
    if (!strcmp(cur_id, item_ids[cur_item_idx])) {
      return cur_item_idx;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

uint32_t IdHtableFindNnt(const char* cur_id, const char* const* item_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size) {
  // returns UINT32_MAX on failure
  for (uint32_t hashval = Hashceil(cur_id, cur_id_slen, id_htable_size); ; ) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if (cur_htable_idval == UINT32_MAX) {
      return UINT32_MAX;
    }
    const uint32_t cur_item_idx = cur_htable_idval & 0x7fffffff;
    if (strequal_unsafe(item_ids[cur_item_idx], cur_id, cur_id_slen)) {
      return cur_item_idx;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

uint32_t IdHtableAdd(const char* cur_id, const char* const* item_ids, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t value, uint32_t* id_htable) {
  for (uint32_t hashval = Hashceil(cur_id, cur_id_slen, id_htable_size); ; ) {
    const uint32_t cur_htable_entry = id_htable[hashval];
    if (cur_htable_entry == UINT32_MAX) {
      id_htable[hashval] = value;
      return UINT32_MAX;
    }
    if (memequal(cur_id, item_ids[cur_htable_entry], cur_id_slen + 1)) {
      return cur_htable_entry;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

uint32_t IdHtableAddNnt(const char* cur_id, const char* const* item_ids, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t value, uint32_t* id_htable) {
  for (uint32_t hashval = Hashceil(cur_id, cur_id_slen, id_htable_size); ; ) {
    const uint32_t cur_htable_entry = id_htable[hashval];
    if (cur_htable_entry == UINT32_MAX) {
      id_htable[hashval] = value;
      return UINT32_MAX;
    }
    if (strequal_unsafe(item_ids[cur_htable_entry], cur_id, cur_id_slen)) {
      return cur_htable_entry;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

CONSTI32(kMaxDupflagThreads, 16);

typedef struct DupflagHtableMakerStruct {
  const uintptr_t* subset_mask;
  const char* const* item_ids;
  uint32_t item_ct;
  uint32_t id_htable_size;
  uint32_t item_uidx_starts[kMaxDupflagThreads];

  uint32_t* id_htable;
} DupflagHtableMaker;

void DupflagHtableMakerMain(uint32_t tidx, uint32_t thread_ct, DupflagHtableMaker* ctx) {
  const uint32_t id_htable_size = ctx->id_htable_size;
  const uintptr_t* subset_mask = ctx->subset_mask;
  const char* const* item_ids = ctx->item_ids;
  const uint32_t item_ct = ctx->item_ct;
  const uint32_t item_uidx_start = ctx->item_uidx_starts[tidx];
  const uint32_t item_idx_start = (item_ct * S_CAST(uint64_t, tidx)) / thread_ct;
  const uint32_t item_idx_end = (item_ct * (S_CAST(uint64_t, tidx) + 1)) / thread_ct;
  uint32_t* id_htable = ctx->id_htable;

  uintptr_t cur_bits;
  uintptr_t item_uidx_base;
  BitIter1Start(subset_mask, item_uidx_start, &item_uidx_base, &cur_bits);
  for (uint32_t item_idx = item_idx_start; item_idx != item_idx_end; ++item_idx) {
    const uintptr_t item_uidx = BitIter1(subset_mask, &item_uidx_base, &cur_bits);
    const char* sptr = item_ids[item_uidx];
    const uint32_t slen = strlen(sptr);
    for (uint32_t hashval = Hashceil(sptr, slen, id_htable_size); ; ) {
      uint32_t old_htable_entry = id_htable[hashval];
      if (old_htable_entry == UINT32_MAX) {
        if (ATOMIC_COMPARE_EXCHANGE_N_U32(&(id_htable[hashval]), &old_htable_entry, item_uidx, 0, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
          break;
        }
      }
      if (strequal_overread(sptr, item_ids[old_htable_entry & 0x7fffffff])) {
        if (!(old_htable_entry & 0x80000000U)) {
          // ok if multiple threads do this
          id_htable[hashval] = old_htable_entry | 0x80000000U;
        }
        break;
      }
      if (++hashval == id_htable_size) {
        hashval = 0;
      }
    }
  }
}

THREAD_FUNC_DECL DupflagHtableMakerThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  DupflagHtableMaker* ctx = S_CAST(DupflagHtableMaker*, arg->sharedp->context);

  // 1. Initialize id_htable with 1-bits in parallel.
  const uint32_t id_htable_size = ctx->id_htable_size;
  const uint32_t thread_ct = GetThreadCt(arg->sharedp) + 1;
  uint32_t* id_htable = ctx->id_htable;
  const uint32_t fill_start = RoundDownPow2((id_htable_size * S_CAST(uint64_t, tidx)) / thread_ct, kInt32PerCacheline);
  const uint32_t fill_end = RoundDownPow2((id_htable_size * (S_CAST(uint64_t, tidx) + 1)) / thread_ct, kInt32PerCacheline);
  SetAllU32Arr(fill_end - fill_start, &(id_htable[fill_start]));

  // 2. sync.Once
  if (THREAD_BLOCK_FINISH(arg)) {
    THREAD_RETURN;
  }

  // 3. Fill hash table in parallel, and then return.
  DupflagHtableMakerMain(tidx, thread_ct, ctx);
  THREAD_RETURN;
}

PglErr MakeDupflagHtable(const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t id_htable_size, uint32_t max_thread_ct, uint32_t* id_htable) {
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  DupflagHtableMaker ctx;
  {
    uint32_t thread_ct = item_ct / 65536;
    if (!thread_ct) {
      thread_ct = 1;
    } else {
      if (thread_ct > max_thread_ct) {
        thread_ct = max_thread_ct;
      }
      // bugfix (26 Nov 2019): forgot to force this
      if (thread_ct > kMaxDupflagThreads) {
        thread_ct = kMaxDupflagThreads;
      }
    }
    if (unlikely(SetThreadCt0(thread_ct - 1, &tg))) {
      goto MakeDupflagHtable_ret_NOMEM;
    }

    ctx.subset_mask = subset_mask;
    ctx.item_ids = item_ids;
    ctx.item_ct = item_ct;
    ctx.id_htable_size = id_htable_size;
    ctx.id_htable = id_htable;

    FillU32SubsetStarts(subset_mask, thread_ct, 0, item_ct, ctx.item_uidx_starts);
    if (thread_ct > 1) {
      SetThreadFuncAndData(DupflagHtableMakerThread, &ctx, &tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto MakeDupflagHtable_ret_THREAD_CREATE_FAIL;
      }
    }
    const uint32_t fill_start = RoundDownPow2((id_htable_size * S_CAST(uint64_t, thread_ct - 1)) / thread_ct, kInt32PerCacheline);
    SetAllU32Arr(id_htable_size - fill_start, &(id_htable[fill_start]));
    if (thread_ct > 1) {
      JoinThreads(&tg);
      DeclareLastThreadBlock(&tg);
      SpawnThreads(&tg);
    }
    DupflagHtableMakerMain(thread_ct - 1, thread_ct, &ctx);
    JoinThreads0(&tg);
  }
  while (0) {
  MakeDupflagHtable_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MakeDupflagHtable_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  return reterr;
}

CONSTI32(kMaxDupstoreThreads, 15);  // probably reduce this after switching to XXH3

typedef struct DupstoreHtableMakerStruct {
  const uintptr_t* subset_mask;
  const char* const* item_ids;
  uint32_t item_ct;
  uint32_t id_htable_size;
  uint32_t* id_htable;

  uint32_t item_uidx_start[2];
  uint32_t* hashes[2];
} DupstoreHtableMaker;

THREAD_FUNC_DECL DupstoreHtableMakerThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  DupstoreHtableMaker* ctx = S_CAST(DupstoreHtableMaker*, arg->sharedp->context);

  // 1. Initialize id_htable with 1-bits in parallel.
  const uint32_t id_htable_size = ctx->id_htable_size;
  const uint32_t thread_ct = GetThreadCt(arg->sharedp);
  uint32_t* id_htable = ctx->id_htable;
  // Add 1 to thread_ct in denominator, since unlike the DupflagHtableMaker
  // case, the parent thread has separate logic to help out here.
  const uint32_t fill_start = RoundDownPow2((id_htable_size * S_CAST(uint64_t, tidx)) / (thread_ct + 1), kInt32PerCacheline);
  uint32_t fill_end = RoundDownPow2((id_htable_size * (S_CAST(uint64_t, tidx) + 1)) / (thread_ct + 1), kInt32PerCacheline);
  SetAllU32Arr(fill_end - fill_start, &(id_htable[fill_start]));

  const uintptr_t* subset_mask = ctx->subset_mask;
  const char* const* item_ids = ctx->item_ids;
  uint32_t items_left = ctx->item_ct;
  const uint32_t idx_start_offset = tidx * kDupstoreBlockSize;
  uint32_t idx_stop_offset = idx_start_offset + kDupstoreBlockSize;
  const uint32_t is_last_thread = (tidx + 1 == thread_ct);
  uint32_t parity = 0;
  while (!THREAD_BLOCK_FINISH(arg)) {
    // 2. Compute up to kDupstoreBlockSize hashes.  (Parent thread is
    //    responsible for using them to update the hash table.)
    if (items_left < idx_stop_offset) {
      if (items_left <= idx_start_offset) {
        // Must be last iteration.  Don't need to update parity.
        // (Could load-balance this iteration differently, but it shouldn't
        // really matter.)
        continue;
      }
      idx_stop_offset = items_left;
    }
    uint32_t* hashes = ctx->hashes[parity];
    // bugfix (24 Jan 2020): IsSet(subset_mask, item_uidx) is not always true.
    uintptr_t item_uidx = FindNth1BitFrom(subset_mask, ctx->item_uidx_start[parity], idx_start_offset + 1);

    uintptr_t cur_bits;
    uintptr_t item_uidx_base;
    BitIter1Start(subset_mask, item_uidx, &item_uidx_base, &cur_bits);
    for (uint32_t item_idx = idx_start_offset; item_idx != idx_stop_offset; ++item_idx) {
      item_uidx = BitIter1(subset_mask, &item_uidx_base, &cur_bits);
      const char* sptr = item_ids[item_uidx];
      const uint32_t slen = strlen(sptr);
      hashes[item_idx] = Hashceil(sptr, slen, id_htable_size);
    }
    items_left -= thread_ct * kDupstoreBlockSize;  // final-iteration underflow ok
    parity = 1 - parity;
    if (is_last_thread) {
      ctx->item_uidx_start[parity] = item_uidx + 1;
    }
  }
  THREAD_RETURN;
}

uint32_t PopulateIdHtableMtDupstoreThreadCt(uint32_t max_thread_ct, uint32_t item_ct) {
  uint32_t thread_ct = item_ct / (2 * kDupstoreBlockSize);
  if (thread_ct >= max_thread_ct) {
    // parent thread is sufficiently busy
    thread_ct = max_thread_ct - 1;
  }
  if (!thread_ct) {
    return 1;
  }
  return MINV(thread_ct, kMaxDupstoreThreads);
}

// dup_ct assumed to be initialized to 0 when dup_ct_ptr != nullptr.
//
// This currently has totally separate code paths for the store_all_dups and
// !store_all_dups cases.  However, the table formats are nearly identical, so
// the code may re-converge, and it's reasonable to have just this API entry
// point.
PglErr PopulateIdHtableMt(unsigned char* arena_top, const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t store_all_dups, uint32_t id_htable_size, uint32_t thread_ct, unsigned char** arena_bottom_ptr, uint32_t* id_htable, uint32_t* dup_ct_ptr) {
  // Change from plink 1.9: if store_all_dups is false, we don't error out on
  // the first encountered duplicate ID; instead, we just flag it in the hash
  // table.  So if '.' is the only duplicate ID, and it never appears in a
  // variant ID list, plink2 never complains.
  //
  // When store_all_dups is true, additional linked lists are allocated past
  // the end of id_htable to track all raw indexes of duplicate names.
  if (!item_ct) {
    return kPglRetSuccess;
  }
  if (!store_all_dups) {
    return MakeDupflagHtable(subset_mask, item_ids, item_ct, id_htable_size, thread_ct, id_htable);
  }
  PglErr reterr = kPglRetSuccess;
  unsigned char* arena_bottom_mark = *arena_bottom_ptr;
  ThreadGroup tg;
  PreinitThreads(&tg);
  DupstoreHtableMaker ctx;
  {
    thread_ct = PopulateIdHtableMtDupstoreThreadCt(thread_ct, item_ct);
    uint32_t item_idx_stop = thread_ct * kDupstoreBlockSize;
    if (arena_end_alloc_u32(arena_bottom_mark, item_idx_stop, &arena_top, &ctx.hashes[0]) ||
        arena_end_alloc_u32(arena_bottom_mark, item_idx_stop, &arena_top, &ctx.hashes[1]) ||
        SetThreadCt(thread_ct, &tg)) {
      goto PopulateIdHtableMt_ret_NOMEM;
    }
    uintptr_t cur_arena_left = arena_top - arena_bottom_mark;
    // We only want to check for out-of-memory once per sync-iteration.  The
    // maximum number of 2x-uint32 entries that might be added per iteration is
    // thread_ct * kDupstoreBlockSize * 2 (*2 is because, when we find the
    // first duplicate of a string, we need to add two new entries instead of
    // just one).
    // htable_dup_base grows from the bottom of the arena, and during the main
    // loop, extra_alloc is the next htable_dup_base[] array index that'll be
    // written to (i.e. it's always an even number).  So if extra_alloc plus
    // twice the aforementioned limit isn't larger than cur_arena_left /
    // sizeof(int32_t), we can't run out of arena space.
    uint32_t extra_alloc_stop;
#ifdef __LP64__
    if (cur_arena_left >= 0x400000000LLU + item_idx_stop * 4 * sizeof(int32_t)) {
      // this can never be hit
      extra_alloc_stop = 0xfffffffe;
    } else {
#endif
      if (unlikely(cur_arena_left < item_idx_stop * 4 * sizeof(int32_t))) {
        goto PopulateIdHtableMt_ret_NOMEM;
      }
      extra_alloc_stop = (cur_arena_left / sizeof(int32_t)) - item_idx_stop * 4;
#ifdef __LP64__
    }
#endif

    ctx.subset_mask = subset_mask;
    ctx.item_ids = item_ids;
    ctx.item_ct = item_ct;
    ctx.id_htable_size = id_htable_size;
    ctx.id_htable = id_htable;
    SetThreadFuncAndData(DupstoreHtableMakerThread, &ctx, &tg);
    if (unlikely(SpawnThreads(&tg))) {
      goto PopulateIdHtableMt_ret_THREAD_CREATE_FAIL;
    }

    const uint32_t fill_start = RoundDownPow2((id_htable_size * S_CAST(uint64_t, thread_ct)) / (thread_ct + 1), kInt32PerCacheline);
    SetAllU32Arr(id_htable_size - fill_start, &(id_htable[fill_start]));

    const uint32_t item_uidx_start = AdvTo1Bit(subset_mask, 0);
    JoinThreads(&tg);
    ctx.item_uidx_start[0] = item_uidx_start;
    if (item_idx_stop >= item_ct) {
      DeclareLastThreadBlock(&tg);
      item_idx_stop = item_ct;
    }
    SpawnThreads(&tg);
    uint32_t items_left = item_ct;
    uint32_t* htable_dup_base = R_CAST(uint32_t*, arena_bottom_mark);
    uint32_t extra_alloc = 0;
    uint32_t prev_llidx = 0;
    uintptr_t cur_bits;
    uintptr_t item_uidx_base;
    BitIter1Start(subset_mask, item_uidx_start, &item_uidx_base, &cur_bits);
    uint32_t parity = 0;
    do {
      JoinThreads(&tg);
      if (extra_alloc > extra_alloc_stop) {
        goto PopulateIdHtableMt_ret_NOMEM;
      }
      if (item_idx_stop < items_left) {
        if (item_idx_stop * 2 >= items_left) {
          DeclareLastThreadBlock(&tg);
        }
        SpawnThreads(&tg);
      } else {
        item_idx_stop = items_left;
      }
      uint32_t* hashes = ctx.hashes[parity];
      for (uint32_t item_idx = 0; item_idx != item_idx_stop; ++item_idx) {
        const uintptr_t item_uidx = BitIter1(subset_mask, &item_uidx_base, &cur_bits);
        const uint32_t hashval_base = hashes[item_idx];
        uint32_t old_htable_entry = id_htable[hashval_base];
        if (old_htable_entry == UINT32_MAX) {
          id_htable[hashval_base] = item_uidx;
          continue;
        }
        const char* sptr = item_ids[item_uidx];
        for (uint32_t hashval = hashval_base; ; ) {
          const uint32_t was_dup = old_htable_entry >> 31;
          uint32_t prev_uidx;
          if (was_dup) {
            prev_llidx = old_htable_entry * 2;
            prev_uidx = htable_dup_base[prev_llidx];
          } else {
            prev_uidx = old_htable_entry;
          }
          if (strequal_overread(sptr, item_ids[prev_uidx])) {
            // point to linked list entry instead
            if (!was_dup) {
              htable_dup_base[extra_alloc] = old_htable_entry;
              htable_dup_base[extra_alloc + 1] = UINT32_MAX;
              prev_llidx = extra_alloc;
              extra_alloc += 2;
            }
            htable_dup_base[extra_alloc] = item_uidx;
            htable_dup_base[extra_alloc + 1] = prev_llidx;
            id_htable[hashval] = 0x80000000U | (extra_alloc >> 1);
            extra_alloc += 2;
            break;  // bugfix
          }
          if (++hashval == id_htable_size) {
            hashval = 0;
          }
          // duplicate this code since we want to avoid the item_ids[item_uidx]
          // read when possible
          old_htable_entry = id_htable[hashval];
          if (old_htable_entry == UINT32_MAX) {
            id_htable[hashval] = item_uidx;
            break;
          }
        }
      }
      items_left -= item_idx_stop;
      parity = 1 - parity;
    } while (items_left);
    // bugfix (15 Oct 2019): no threads to join here!
    if (extra_alloc) {
      // bugfix: forgot to align this
      *arena_bottom_ptr += RoundUpPow2(extra_alloc * sizeof(int32_t), kCacheline);
      if (dup_ct_ptr) {
        *dup_ct_ptr = extra_alloc / 2;
      }
    }
  }
  while (0) {
  PopulateIdHtableMt_ret_NOMEM:
    *arena_bottom_ptr = arena_bottom_mark;
    reterr = kPglRetNomem;
    break;
  PopulateIdHtableMt_ret_THREAD_CREATE_FAIL:
    // not currently possible for this to happen after *arena_bottom_ptr moved
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  return reterr;
}

// Similar to DupflagHtableMaker, but we don't need to flag duplicates, we just
// "error out" when we find one.
typedef struct NondupHtableMakerStruct {
  const uintptr_t* subset_mask;
  const char* const* item_ids;
  uint32_t item_ct;
  uint32_t id_htable_size;
  uint32_t item_uidx_starts[kMaxDupflagThreads];

  uint32_t* id_htable;

  uint32_t dup_found;
} NondupHtableMaker;

void NondupHtableMakerMain(uint32_t tidx, uint32_t thread_ct, NondupHtableMaker* ctx) {
  const uint32_t id_htable_size = ctx->id_htable_size;
  const uintptr_t* subset_mask = ctx->subset_mask;
  const char* const* item_ids = ctx->item_ids;
  const uint32_t item_ct = ctx->item_ct;
  const uint32_t item_uidx_start = ctx->item_uidx_starts[tidx];
  const uint32_t item_idx_start = (item_ct * S_CAST(uint64_t, tidx)) / thread_ct;
  const uint32_t item_idx_end = (item_ct * (S_CAST(uint64_t, tidx) + 1)) / thread_ct;
  uint32_t* id_htable = ctx->id_htable;

  uintptr_t cur_bits;
  uintptr_t item_uidx_base;
  BitIter1Start(subset_mask, item_uidx_start, &item_uidx_base, &cur_bits);
  for (uint32_t item_idx = item_idx_start; item_idx != item_idx_end; ) {
    const uint32_t item_idx_cur_end = MINV(item_idx_end, item_idx + 65536);
    for (; item_idx != item_idx_cur_end; ++item_idx) {
      const uintptr_t item_uidx = BitIter1(subset_mask, &item_uidx_base, &cur_bits);
      const char* sptr = item_ids[item_uidx];
      const uint32_t slen = strlen(sptr);
      for (uint32_t hashval = Hashceil(sptr, slen, id_htable_size); ; ) {
        uint32_t old_htable_entry = id_htable[hashval];
        if (old_htable_entry == UINT32_MAX) {
          if (ATOMIC_COMPARE_EXCHANGE_N_U32(&(id_htable[hashval]), &old_htable_entry, item_uidx, 0, __ATOMIC_ACQ_REL, __ATOMIC_ACQUIRE)) {
            break;
          }
        }
        if (strequal_overread(sptr, item_ids[old_htable_entry & 0x7fffffff])) {
          // no synchronization needed since this variable can't change in any
          // other way
          ctx->dup_found = 1;
          return;
        }
        if (++hashval == id_htable_size) {
          hashval = 0;
        }
      }
    }
    if (ctx->dup_found) {
      return;
    }
  }
}

THREAD_FUNC_DECL NondupHtableMakerThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  NondupHtableMaker* ctx = S_CAST(NondupHtableMaker*, arg->sharedp->context);

  // 1. Initialize id_htable with 1-bits in parallel.
  const uint32_t id_htable_size = ctx->id_htable_size;
  const uint32_t thread_ct = GetThreadCt(arg->sharedp) + 1;
  uint32_t* id_htable = ctx->id_htable;
  const uint32_t fill_start = RoundDownPow2((id_htable_size * S_CAST(uint64_t, tidx)) / thread_ct, kInt32PerCacheline);
  const uint32_t fill_end = RoundDownPow2((id_htable_size * (S_CAST(uint64_t, tidx) + 1)) / thread_ct, kInt32PerCacheline);
  SetAllU32Arr(fill_end - fill_start, &(id_htable[fill_start]));

  // 2. sync.Once
  if (THREAD_BLOCK_FINISH(arg)) {
    THREAD_RETURN;
  }

  // 3. Fill hash table in parallel, and then return.
  NondupHtableMakerMain(tidx, thread_ct, ctx);
  THREAD_RETURN;
}

PglErr MakeNondupHtable(const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t id_htable_size, uint32_t max_thread_ct, uint32_t* id_htable, uint32_t* dup_found_ptr) {
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  NondupHtableMaker ctx;
  {
    uint32_t thread_ct = item_ct / 65536;
    if (!thread_ct) {
      thread_ct = 1;
    } else {
      if (thread_ct > max_thread_ct) {
        thread_ct = max_thread_ct;
      }
      if (thread_ct > kMaxDupflagThreads) {
        thread_ct = kMaxDupflagThreads;
      }
    }
    if (unlikely(SetThreadCt0(thread_ct - 1, &tg))) {
      goto MakeNondupHtable_ret_NOMEM;
    }

    ctx.subset_mask = subset_mask;
    ctx.item_ids = item_ids;
    ctx.item_ct = item_ct;
    ctx.id_htable_size = id_htable_size;
    ctx.id_htable = id_htable;
    ctx.dup_found = 0;

    uint32_t item_uidx = AdvTo1Bit(subset_mask, 0);
    uint32_t item_idx = 0;
    ctx.item_uidx_starts[0] = item_uidx;
    for (uintptr_t tidx = 1; tidx != thread_ct; ++tidx) {
      const uint32_t item_idx_new = (item_ct * S_CAST(uint64_t, tidx)) / thread_ct;
      item_uidx = FindNth1BitFrom(subset_mask, item_uidx + 1, item_idx_new - item_idx);
      ctx.item_uidx_starts[tidx] = item_uidx;
      item_idx = item_idx_new;
    }

    if (thread_ct > 1) {
      SetThreadFuncAndData(NondupHtableMakerThread, &ctx, &tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto MakeNondupHtable_ret_THREAD_CREATE_FAIL;
      }
    }
    const uint32_t fill_start = RoundDownPow2((id_htable_size * S_CAST(uint64_t, thread_ct - 1)) / thread_ct, kInt32PerCacheline);
    SetAllU32Arr(id_htable_size - fill_start, &(id_htable[fill_start]));
    if (thread_ct > 1) {
      JoinThreads(&tg);
      DeclareLastThreadBlock(&tg);
      SpawnThreads(&tg);
    }
    NondupHtableMakerMain(thread_ct - 1, thread_ct, &ctx);
    JoinThreads0(&tg);
    *dup_found_ptr = ctx.dup_found;
  }
  while (0) {
  MakeNondupHtable_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MakeNondupHtable_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  return reterr;
}

PglErr AllocAndPopulateNondupHtableMt(unsigned char* arena_top, const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t max_thread_ct, unsigned char** arena_bottom_ptr, uint32_t** id_htable_ptr, uint32_t* id_htable_size_ptr, uint32_t* dup_found_ptr) {
  uint32_t id_htable_size = GetHtableFastSize(item_ct);
  const uintptr_t wkspace_left_rd = RoundDownPow2(S_CAST(uintptr_t, arena_top - (*arena_bottom_ptr)), kCacheline);
  if (wkspace_left_rd < id_htable_size * sizeof(int32_t)) {
    id_htable_size = wkspace_left_rd / sizeof(int32_t);
    const uint32_t min_htable_size = GetHtableMinSize(item_ct);
    if (id_htable_size < min_htable_size) {
      return kPglRetNomem;
    }
  }
  *id_htable_ptr = S_CAST(uint32_t*, arena_alloc_raw_rd(id_htable_size * sizeof(int32_t), arena_bottom_ptr));
  *id_htable_size_ptr = id_htable_size;
  return MakeNondupHtable(subset_mask, item_ids, item_ct, id_htable_size, max_thread_ct, *id_htable_ptr, dup_found_ptr);
}

void NondupIdLoadProcessTokens(const char* const* item_ids, const uint32_t* item_id_htable, const char* shard_start, const char* shard_end, uint32_t item_id_htable_size, uintptr_t* already_seen) {
  const char* shard_iter = shard_start;
  while (1) {
    shard_iter = FirstPostspaceBounded(shard_iter, shard_end);
    if (shard_iter == shard_end) {
      return;
    }
    const char* token_end = CurTokenEnd(shard_iter);
    const uint32_t item_uidx = IdHtableFindNnt(shard_iter, item_ids, item_id_htable, token_end - shard_iter, item_id_htable_size);
    shard_iter = token_end;
    if (item_uidx == UINT32_MAX) {
      continue;
    }
    SetBit(item_uidx, already_seen);
  }
}

CONSTI32(kMaxNondupIdLoadThreads, 8);

typedef struct NondupIdLoadCtxStruct {
  const char* const* item_ids;
  const uint32_t* item_id_htable;
  uintptr_t item_id_htable_size;

  char* shard_boundaries[kMaxNondupIdLoadThreads + 1];
  uintptr_t* already_seens[kMaxNondupIdLoadThreads];
} NondupIdLoadCtx;

THREAD_FUNC_DECL NondupIdLoadThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx_p1 = arg->tidx + 1;
  NondupIdLoadCtx* ctx = S_CAST(NondupIdLoadCtx*, arg->sharedp->context);

  const char* const* item_ids = ctx->item_ids;
  const uint32_t* item_id_htable = ctx->item_id_htable;
  const uintptr_t item_id_htable_size = ctx->item_id_htable_size;
  uintptr_t* already_seen = ctx->already_seens[tidx_p1];
  do {
    NondupIdLoadProcessTokens(item_ids, item_id_htable, ctx->shard_boundaries[tidx_p1], ctx->shard_boundaries[tidx_p1 + 1], item_id_htable_size, already_seen);
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr NondupIdLoad(unsigned char* arena_bottom, unsigned char* arena_top, const char* const* item_ids, const char* fname, uint32_t raw_item_ct, uint32_t item_ct, uint32_t max_thread_ct, uintptr_t* item_set, uint32_t* dup_found_ptr, char* errbuf) {
  PglErr reterr = kPglRetSuccess;
  TokenStream tks;
  ThreadGroup tg;
  PreinitTokenStream(&tks);
  PreinitThreads(&tg);
  NondupIdLoadCtx ctx;
  errbuf[0] = '\0';
  {
    const uint32_t calc_thread_ct_m1 = MINV(max_thread_ct, kMaxNondupIdLoadThreads) - 1;
    if (unlikely(SetThreadCt0(calc_thread_ct_m1, &tg))) {
      goto NondupIdLoad_ret_NOMEM;
    }
    if (!item_ct) {
      goto NondupIdLoad_ret_1;
    }
    uint32_t decompress_thread_ct = 1;
    if (max_thread_ct > calc_thread_ct_m1 + 2) {
      decompress_thread_ct = max_thread_ct - calc_thread_ct_m1 - 1;
    }
    const uint32_t raw_item_ctl = BitCtToWordCt(raw_item_ct);
    for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
      if (unlikely(arena_calloc_w(arena_top, raw_item_ctl, &arena_bottom, &(ctx.already_seens[tidx])))) {
        goto NondupIdLoad_ret_NOMEM;
      }
    }
    if (unlikely(S_CAST(uintptr_t, arena_top - arena_bottom) < kTokenStreamBlen)) {
      goto NondupIdLoad_ret_NOMEM;
    }
    char* dst = S_CAST(char*, arena_end_alloc_raw(kTokenStreamBlen, &arena_top));
    reterr = TextStreamOpenEx(fname, 0, kTokenStreamBlen, decompress_thread_ct, nullptr, dst, &(tks.txs));
    if (unlikely(reterr)) {
      goto NondupIdLoad_ret_TKSTREAM_FAIL;
    }
    // Allocate/initialize hash table last, so we can make the most appropriate
    // speed/size tradeoff.
    uint32_t* item_id_htable;
    uint32_t item_id_htable_size;
    reterr = AllocAndPopulateNondupHtableMt(arena_top, item_set, item_ids, item_ct, max_thread_ct, &arena_bottom, &item_id_htable, &item_id_htable_size, dup_found_ptr);
    if (reterr || (*dup_found_ptr)) {
      goto NondupIdLoad_ret_1;
    }
    if (calc_thread_ct_m1) {
      ctx.item_ids = item_ids;
      ctx.item_id_htable = item_id_htable;
      ctx.item_id_htable_size = item_id_htable_size;
      SetThreadFuncAndData(NondupIdLoadThread, &ctx, &tg);
    }
    ZeroWArr(raw_item_ctl, item_set);
    while (1) {
      reterr = TksNext(&tks, calc_thread_ct_m1 + 1, ctx.shard_boundaries);
      if (reterr) {
        break;
      }
      if (calc_thread_ct_m1) {
        if (unlikely(SpawnThreads(&tg))) {
          goto NondupIdLoad_ret_THREAD_CREATE_FAIL;
        }
      }
      NondupIdLoadProcessTokens(item_ids, item_id_htable, ctx.shard_boundaries[0], ctx.shard_boundaries[1], item_id_htable_size, item_set);
      JoinThreads0(&tg);
    }
    if (unlikely(reterr != kPglRetEof)) {
      goto NondupIdLoad_ret_TKSTREAM_FAIL;
    }
    reterr = kPglRetSuccess;
    for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
      BitvecOr(ctx.already_seens[tidx], raw_item_ctl, item_set);
    }
  }
  while (0) {
  NondupIdLoad_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  NondupIdLoad_ret_TKSTREAM_FAIL:
    // copy TextErrPrint() logic for now, will probably move to plink2_text
    {
      reterr = TokenStreamErrcode(&tks);
      const char* errmsg = TokenStreamError(&tks);
      if (reterr == kPglRetOpenFail) {
        snprintf(errbuf, 2 * kMaxMediumLine, "Error: Failed to open %s : %s.\n", fname, errmsg);
      } else if (reterr == kPglRetReadFail) {
        snprintf(errbuf, 2 * kMaxMediumLine, "Error: %s read failure: %s.\n", fname, errmsg);
      } else if (reterr == kPglRetDecompressFail) {
        snprintf(errbuf, 2 * kMaxMediumLine, "Error: %s decompression failure: %s.\n", fname, errmsg);
      } else if (reterr == kPglRetMalformedInput) {
        snprintf(errbuf, 2 * kMaxMediumLine, "Error: Pathologically long token in %s.\n", fname);
      }
    }
    break;
  NondupIdLoad_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
  }
 NondupIdLoad_ret_1:
  CleanupThreads(&tg);
  if (unlikely(CleanupTokenStream(&tks, &reterr))) {
    snprintf(errbuf, 2 * kMaxMediumLine, "Error: %s read failure: %s.\n", fname, strerror(errno));
  }
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
