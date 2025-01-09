#ifndef __PLINK2_HTABLE_H__
#define __PLINK2_HTABLE_H__

// This library is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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


// This implements a hash table over a const char* const* array of strings.

#include "plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// copy_subset() doesn't exist since a loop of the form
//   uintptr_t uidx_base = 0;
//   uintptr_t cur_bits = subset_mask[0];
//   for (uint32_t idx = 0; idx != subset_size; ++idx) {
//     const uintptr_t uidx = BitIter1(subset_mask, &uidx_base, &cur_bits);
//     *target_iter++ = source_arr[uidx];
//   }
// seems to compile better?
//
// similarly, copy_when_nonmissing() retired in favor of genoarr_to_nonmissing
// followed by a for loop

// Benchmark results justify replacing this with XXH3 when the latter is ready.
// Not considering it ready yet since the development XXH_INLINE_ALL setup
// isn't very friendly.
//
// (Can also try using _mm_aesenc_si128() in a similar manner to the Golang
// runtime when SSE4.2+ is available, and otherwise keeping Murmurhash.  See
// aeshashbody() in https://golang.org/src/runtime/asm_amd64.s .)
//
// Eventually want this to be a constexpr?  Seems painful to make that work,
// though.
uint32_t Hash32(const void* key, uint32_t len);

// see http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
// Note that this is a bit more vulnerable to adversarial input: modulo
// reduction requires lots of hash collisions (or near-collisions) or known
// hash table size, while this can be attacked without knowledge of the table
// size.  But it doesn't really matter; either way, you'd need to add something
// like a randomly generated salt if you care about defending.
HEADER_INLINE uint32_t Hashceil(const char* idstr, uint32_t idlen, uint32_t htable_size) {
  return (S_CAST(uint64_t, Hash32(idstr, idlen)) * htable_size) >> 32;
}

// In most cases, plink2 represents an array of strings in one of the following
// two ways:
//   (const char* strbox, uint32_t str_ct, uintptr_t max_str_blen):
//     null-terminated string #x starts at &(strbox[x * max_str_blen)), where
//     x is a 0-based index.
//   (const char* const* item_ids, uint32_t str_ct):
//     null-terminated string #x starts at item_ids[x].
// When we need to perform string -> string-index lookups into the array, we
// construct a hashmap as follows:
// - Allocate a uint32_t array of htable_size >= 2 * str_ct, and initialize all
//   entries to UINT32_MAX to mark them empty.
// - Iterate through the array, computing hashval := Hashceil(str, strlen(str),
//   htable_size) for each string, and then setting htable[hashval] :=
//   string-index whenever that htable entry is empty.  When there is a
//   conflict, use linear probing (increment hashval until an empty table cell
//   is found, wrapping around from (htable_size - 1) to 0).
//   - PLINK 1.9 used quadratic probing, but that's been scrapped since
//     benchmark results suggest that it has no meaningful advantage.  We
//     aren't dealing with adversarial input...
// The "StrboxHtable" functions below work with the const char* strbox
// string-array representation, while "IdHtable" functions work with the const
// char* const* item_ids representation.

// uintptr_t geqprime(uintptr_t floor);

// assumes ceil is odd and greater than 4.  Returns the first prime <= ceil.
// uintptr_t leqprime(uintptr_t ceil);

HEADER_INLINE uint32_t GetHtableMinSize(uintptr_t item_ct) {
  // Don't need primes any more.
  return item_ct * 2;
}

// load factor ~22% seems to yield the best speed/space tradeoff on my test
// machines
HEADER_INLINE uint32_t GetHtableFastSize(uint32_t item_ct) {
  if (item_ct < 954437176) {
    return (item_ct * 9LLU) / 2;
  }
  return 4294967290U;
}

// cur_id DOES need to be null-terminated
uint32_t IdHtableFind(const char* cur_id, const char* const* item_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size);

// item_ids overread must be ok.  In exchange, cur_id doesn't need to be
// null-terminated any more.
uint32_t IdHtableFindNnt(const char* cur_id, const char* const* item_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size);

HEADER_INLINE void HtableAddNondup(const char* cur_id, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t value, uint32_t* id_htable) {
  for (uint32_t hashval = Hashceil(cur_id, cur_id_slen, id_htable_size); ; ) {
    const uint32_t cur_htable_entry = id_htable[hashval];
    if (cur_htable_entry == UINT32_MAX) {
      id_htable[hashval] = value;
      return;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

// Assumes cur_id is null-terminated.
// item_ids overread must be ok.
// Returns string-index if cur_id is already in the table, UINT32_MAX if it was
// added.
uint32_t IdHtableAdd(const char* cur_id, const char* const* item_ids, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t value, uint32_t* id_htable);

// Does not require cur_id to be null-terminated.
uint32_t IdHtableAddNnt(const char* cur_id, const char* const* item_ids, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t value, uint32_t* id_htable);

CONSTI32(kDupstoreBlockSize, 65536);
CONSTI32(kDupstoreThreadWkspace, kDupstoreBlockSize * 2 * sizeof(int32_t));

uint32_t PopulateIdHtableMtDupstoreThreadCt(uint32_t max_thread_ct, uint32_t item_ct);

// store_all_dups currently must be true for dup_ct to be set, but this is easy
// to change
PglErr PopulateIdHtableMt(unsigned char* arena_top, const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t store_all_dups, uint32_t id_htable_size, uint32_t thread_ct, unsigned char** arena_bottom_ptr, uint32_t* id_htable, uint32_t* dup_ct_ptr);

PglErr MakeNondupHtable(const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t id_htable_size, uint32_t max_thread_ct, uint32_t* id_htable, uint32_t* dup_found_ptr);

PglErr AllocAndPopulateNondupHtableMt(unsigned char* arena_top, const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t max_thread_ct, unsigned char** arena_bottom_ptr, uint32_t** id_htable_ptr, uint32_t* id_htable_size_ptr, uint32_t* dup_found_ptr);

HEADER_INLINE PglErr CheckIdUniqueness(unsigned char* arena_bottom, unsigned char* arena_top, const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t max_thread_ct, uint32_t* dup_found_ptr) {
  uint32_t* id_htable;
  uint32_t id_htable_size;
  return AllocAndPopulateNondupHtableMt(arena_top, subset_mask, item_ids, item_ct, max_thread_ct, &arena_bottom, &id_htable, &id_htable_size, dup_found_ptr);
}

// item_set is in/out: it should be set to item_include on input, but it's
// overwritten with the loaded item set on successful return.
// errbuf assumed to have length at least 2 * kMaxMediumLine.
PglErr NondupIdLoad(unsigned char* arena_bottom, unsigned char* arena_top, const char* const* item_ids, const char* fname, uint32_t raw_item_ct, uint32_t item_ct, uint32_t max_thread_ct, uintptr_t* item_set, uint32_t* dup_found_ptr, char* errbuf);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_HTABLE_H__
