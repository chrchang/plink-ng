#ifndef __PLINK2_CMDLINE_H__
#define __PLINK2_CMDLINE_H__

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


// Cross-platform logging, command-line parsing, workspace
// initialization/allocation, basic multithreading code, and a few numeric
// constants.

#include <assert.h>
#include <stdarg.h>
#include <string.h>

#include "include/plink2_base.h"
#include "include/plink2_bits.h"
#include "include/plink2_memory.h"
#include "include/plink2_string.h"

#ifdef _WIN32
#  define NULL_STREAM_NAME "nul"
#else
#  define NULL_STREAM_NAME "/dev/null"
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

CONSTI32(kLogbufSize, 2 * kMaxMediumLine);

// must be at least 2 * kMaxMediumLine + 2 to support generic token loader.
CONSTI32(kTextbufSize, 2 * kMaxMediumLine + 256);

// when g_textbuf is used as a generic I/O buffer, this is a convenient
// power-of-2 size (must be <= kTextbufSize).
CONSTI32(kTextbufMainSize, 2 * kMaxMediumLine);

// "slen" is now used to indicate string lengths excluding terminating nulls,
// while "blen" includes the terminator.

// Maximum length of chromosome, variant, FID, IID, cluster, and set IDs (not
// including terminating null).  This value supports up to 8 IDs per line
// (maximum so far is 5, for e.g. --hom).
// Assumed by plink2_pvar to be a multiple of 16.
CONSTI32(kMaxIdSlen, 16000);
CONSTI32(kMaxIdBlen, kMaxIdSlen + 1);
// Don't see a better option than #define for this.
#define MAX_ID_SLEN_STR "16000"

// allow extensions like .model.trend.fisher.set.score.adjusted
CONSTI32(kMaxOutfnameExtBlen, 39);

// All global variables not initialized at compile time start with g_ (even if
// they're initialized very early and never changed afterwards, like
// g_one_char_strs).
extern char g_textbuf[];

extern const char* g_one_char_strs;

extern FILE* g_logfile;

// Mostly-safe log buffer (length kLogbufSize, currently 256k).  Good practice
// to use snprintf when writing an entire line to it in a single statement.
// Warning: Do NOT put allele codes or arbitrary-length lists in here.
extern char g_logbuf[];

extern uint32_t g_debug_on;
extern uint32_t g_log_failed;

// for --warning-errcode
extern uint32_t g_stderr_written_to;


// Warning: Do NOT include allele codes (unless they're guaranteed to be SNPs)
// in log strings; they can overflow the buffer.
void logputs_silent(const char* str);

void logputs(const char* str);

void logerrputs(const char* str);

void logputsb();

void logerrputsb();

HEADER_INLINE void logprintf(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vsnprintf(g_logbuf, kLogbufSize, fmt, args);
  logputsb();
}

HEADER_INLINE void logerrprintf(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vsnprintf(g_logbuf, kLogbufSize, fmt, args);
  logerrputsb();
}

// input for WordWrapB/logprintfww should have no intermediate '\n's.  If
// suffix_len is 0, there should be a terminating \n.
HEADER_INLINE void WordWrapB(uint32_t suffix_len) {
  WordWrap(suffix_len, g_logbuf);
}

HEADER_INLINE void logpreprintfww(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vsnprintf(g_logbuf, kLogbufSize, fmt, args);
  WordWrapB(0);
}

HEADER_INLINE void logprintfww(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vsnprintf(g_logbuf, kLogbufSize, fmt, args);
  WordWrapB(0);
  logputsb();
}

HEADER_INLINE void logerrprintfww(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vsnprintf(g_logbuf, kLogbufSize, fmt, args);
  WordWrapB(0);
  logerrputsb();
}

// 5 = length of "done." suffix, which is commonly used
HEADER_INLINE void logprintfww5(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vsnprintf(g_logbuf, kLogbufSize, fmt, args);
  WordWrapB(5);
  logputsb();
}

HEADER_INLINE void DebugPrintf(const char* fmt, ...) {
  if (g_debug_on) {
    va_list args;
    va_start(args, fmt);
    vsnprintf(g_logbuf, kLogbufSize, fmt, args);
    logputsb();
  }
}

HEADER_INLINE void DPrintf(const char* fmt, ...) {
  if (g_debug_on) {
    va_list args;
    va_start(args, fmt);
    vsnprintf(g_logbuf, kLogbufSize, fmt, args);
    logputsb();
  }
}

uint32_t FileExists(const char* fname);

// Returns kPglRetOpenFail if file doesn't exist, or kPglRetRewindFail if file
// is process-substitution/named-pipe.  Does not print an error message.
PglErr ForceNonFifo(const char* fname);

BoolErr fopen_checked(const char* fname, const char* mode, FILE** target_ptr);

HEADER_INLINE IntErr fputs_checked(const char* str, FILE* outfile) {
  fputs(str, outfile);
  return ferror_unlocked(outfile);
}

BoolErr fwrite_flush2(char* buf_flush, FILE* outfile, char** write_iter_ptr);

HEADER_INLINE BoolErr fwrite_uflush2(unsigned char* buf_flush, FILE* outfile, unsigned char** write_iter_ptr) {
  return fwrite_flush2(R_CAST(char*, buf_flush), outfile, R_CAST(char**, write_iter_ptr));
}

HEADER_INLINE BoolErr fwrite_ck(char* buf_flush, FILE* outfile, char** write_iter_ptr) {
  if ((*write_iter_ptr) < buf_flush) {
    return 0;
  }
  return fwrite_flush2(buf_flush, outfile, write_iter_ptr);
}

// fclose_null defined in plink2_base.h

BoolErr fclose_flush_null(char* buf_flush, char* write_iter, FILE** outfile_ptr);

HEADER_INLINE BoolErr fclose_uflush_null(unsigned char* buf_flush, unsigned char* write_iter, FILE** outfile_ptr) {
  return fclose_flush_null(R_CAST(char*, buf_flush), R_CAST(char*, write_iter), outfile_ptr);
}

// This should only be used when the file can only be open on error-early-exit;
// otherwise we care about fclose's error code.
HEADER_INLINE void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

HEADER_INLINE uint32_t ClipU32(uint32_t val, uint32_t lbound, uint32_t ubound) {
  if (val >= ubound) {
    return ubound;
  }
  return MAXV(val, lbound);
}

int32_t u32cmp(const void* aa, const void* bb);

int32_t double_cmp(const void* aa, const void* bb);

int32_t double_cmp_decr(const void* aa, const void* bb);

int32_t u64cmp(const void* aa, const void* bb);

#ifndef __cplusplus
int32_t u64cmp_decr(const void* aa, const void* bb);
#endif

HEADER_INLINE uint32_t U32ArrMax(const uint32_t* unsorted_arr, uintptr_t len) {
  const uint32_t* unsorted_arr_end = &(unsorted_arr[len]);
#ifndef __cplusplus
  const uint32_t* unsorted_arr_iter = unsorted_arr;
  uint32_t uimax = *unsorted_arr_iter++;
  while (unsorted_arr_iter < unsorted_arr_end) {
    const uint32_t cur_val = *unsorted_arr_iter++;
    if (cur_val > uimax) {
      uimax = cur_val;
    }
  }
  return uimax;
#else
  return *std::max_element(unsorted_arr, unsorted_arr_end);
#endif
}

HEADER_INLINE float FArrMax(const float* unsorted_arr, uintptr_t len) {
  const float* unsorted_arr_end = &(unsorted_arr[len]);
#if defined(__APPLE__) || !defined(__cplusplus)
  // std::max_element doesn't seem to be performant for floats/doubles on OS X
  const float* unsorted_arr_iter = unsorted_arr;
  float fmax = *unsorted_arr_iter++;
  while (unsorted_arr_iter < unsorted_arr_end) {
    const float cur_val = *unsorted_arr_iter++;
    if (cur_val > fmax) {
      fmax = cur_val;
    }
  }
  return fmax;
#else
  return *std::max_element(unsorted_arr, unsorted_arr_end);
#endif
}

HEADER_INLINE double DArrMax(const double* unsorted_arr, uintptr_t len) {
  const double* unsorted_arr_end = &(unsorted_arr[len]);
#if defined(__APPLE__) || !defined(__cplusplus)
  const double* unsorted_arr_iter = unsorted_arr;
  double dmax = *unsorted_arr_iter++;
  while (unsorted_arr_iter < unsorted_arr_end) {
    const double cur_val = *unsorted_arr_iter++;
    if (cur_val > dmax) {
      dmax = cur_val;
    }
  }
  return dmax;
#else
  return *std::max_element(unsorted_arr, unsorted_arr_end);
#endif
}

float DestructiveMedianF(uintptr_t len, float* unsorted_arr);

double DestructiveMedianD(uintptr_t len, double* unsorted_arr);


// This makes a temporary g_bigstack allocation.
// Must be safe to read up to (kBytesPerWord - 1) bytes past end of strbox.
// Results may technically vary between runs when duplicate elements are
// present; it's assumed that this doesn't matter because all duplicates will
// be handled in the same manner.
// id_map is new_to_old_idx.
BoolErr SortStrboxIndexed(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map);

// Offset of std::lower_bound, i.e. # of elements < needle.
// Requires arr_length > 0.
uintptr_t LowerBoundNonemptyU32(const uint32_t* sorted_u32_arr, uintptr_t arr_length, uint32_t needle);

HEADER_INLINE uintptr_t LowerBoundU32(const uint32_t* sorted_u32_arr, uintptr_t arr_length, uint32_t needle) {
  if (!arr_length) {
    return 0;
  }
  return LowerBoundNonemptyU32(sorted_u32_arr, arr_length, needle);
}

HEADER_INLINE uintptr_t LowerBoundConstrainedNonemptyU32(const uint32_t* sorted_u32_arr, uintptr_t start_idx, uintptr_t end_idx, uint32_t needle) {
  return start_idx + LowerBoundNonemptyU32(&(sorted_u32_arr[start_idx]), end_idx - start_idx, needle);
}

// Same as above, but optimized for case where result is probably close to the
// front, and tolerates end_idx (== arr_length) == 0.
uintptr_t Expsearch0U32(const uint32_t* sorted_u32_arr, uintptr_t end_idx, uint32_t needle);

HEADER_INLINE uintptr_t ExpsearchU32(const uint32_t* sorted_u32_arr, uintptr_t start_idx, uintptr_t end_idx, uint32_t needle) {
  return start_idx + Expsearch0U32(&(sorted_u32_arr[start_idx]), end_idx - start_idx, needle);
}

#ifdef __LP64__
uintptr_t Expsearch0U64(const uint64_t* sorted_u64_arr, uintptr_t end_idx, uint64_t needle);

HEADER_INLINE uintptr_t ExpsearchU64(const uint64_t* sorted_u64_arr, uintptr_t start_idx, uintptr_t end_idx, uint64_t needle) {
  return start_idx + Expsearch0U64(&(sorted_u64_arr[start_idx]), end_idx - start_idx, needle);
}

HEADER_INLINE uintptr_t ExpsearchW(const uintptr_t* sorted_w_arr, uintptr_t start_idx, uintptr_t end_idx, uintptr_t needle) {
  return ExpsearchU64(R_CAST(const uint64_t*, sorted_w_arr), start_idx, end_idx, needle);
}
#else
HEADER_INLINE uintptr_t ExpsearchW(const uintptr_t* sorted_w_arr, uintptr_t start_idx, uintptr_t end_idx, uintptr_t needle) {
  return ExpsearchU32(R_CAST(const uint32_t*, sorted_w_arr), start_idx, end_idx, needle);
}
#endif

uintptr_t LowerBoundNonemptyU64(const uint64_t* sorted_u64_arr, uintptr_t arr_length, uint64_t needle);

HEADER_INLINE uintptr_t LowerBoundNonemptyW(const uintptr_t* sorted_w_arr, uintptr_t arr_length, uintptr_t needle) {
#ifdef __LP64__
  return LowerBoundNonemptyU64(R_CAST(const uint64_t*, sorted_w_arr), arr_length, needle);
#else
  return LowerBoundNonemptyU32(R_CAST(const uint32_t*, sorted_w_arr), arr_length, needle);
#endif
}

HEADER_INLINE uintptr_t LowerBoundConstrainedNonemptyW(const uintptr_t* sorted_w_arr, uintptr_t start_idx, uintptr_t end_idx, uintptr_t needle) {
  return start_idx + LowerBoundNonemptyW(&(sorted_w_arr[start_idx]), end_idx - start_idx, needle);
}

uintptr_t LowerBoundNonemptyD(const double* sorted_dbl_arr, uintptr_t arr_length, double needle);

uintptr_t Expsearch0D(const double* sorted_dbl_arr, uintptr_t end_idx, double needle);

HEADER_INLINE uintptr_t ExpsearchD(const double* sorted_dbl_arr, uintptr_t start_idx, uintptr_t end_idx, double needle) {
  return start_idx + Expsearch0D(&(sorted_dbl_arr[start_idx]), end_idx - start_idx, needle);
}

uintptr_t UpperBoundNonemptyU32(const uint32_t* sorted_u32_arr, uintptr_t arr_length, uint32_t needle);

HEADER_INLINE uintptr_t UpperBoundConstrainedNonemptyU32(const uint32_t* sorted_u32_arr, uintptr_t start_idx, uintptr_t end_idx, uint32_t needle) {
  return start_idx + UpperBoundNonemptyU32(&(sorted_u32_arr[start_idx]), end_idx - start_idx, needle);
}

uintptr_t UpperBoundNonemptyU64(const uint64_t* sorted_u64_arr, uintptr_t arr_length, uint64_t needle);

HEADER_INLINE uintptr_t UpperBoundConstrainedNonemptyU64(const uint64_t* sorted_u64_arr, uintptr_t start_idx, uintptr_t end_idx, uint64_t needle) {
  return start_idx + UpperBoundNonemptyU64(&(sorted_u64_arr[start_idx]), end_idx - start_idx, needle);
}

HEADER_INLINE uintptr_t UpperBoundNonemptyW(const uintptr_t* sorted_w_arr, uintptr_t arr_length, uintptr_t needle) {
#ifdef __LP64__
  return UpperBoundNonemptyU64(R_CAST(const uint64_t*, sorted_w_arr), arr_length, needle);
#else
  return UpperBoundNonemptyU32(R_CAST(const uint32_t*, sorted_w_arr), arr_length, needle);
#endif
}

HEADER_INLINE uintptr_t UpperBoundConstrainedNonemptyW(const uintptr_t* sorted_w_arr, uintptr_t start_idx, uintptr_t end_idx, uintptr_t needle) {
  return start_idx + UpperBoundNonemptyW(&(sorted_w_arr[start_idx]), end_idx - start_idx, needle);
}

// Requires sorted_w_arr[start_idx] <= needle.
HEADER_INLINE uintptr_t LastLeqU32(const uint32_t* sorted_u32_arr, uintptr_t start_idx, uintptr_t end_idx, uint32_t needle) {
  return UpperBoundConstrainedNonemptyU32(sorted_u32_arr, start_idx, end_idx, needle) - 1;
}

HEADER_INLINE uintptr_t LastLeqU64(const uint64_t* sorted_u64_arr, uintptr_t start_idx, uintptr_t end_idx, uint64_t needle) {
  return UpperBoundConstrainedNonemptyU64(sorted_u64_arr, start_idx, end_idx, needle) - 1;
}

HEADER_INLINE uintptr_t LastLeqW(const uintptr_t* sorted_w_arr, uintptr_t start_idx, uintptr_t end_idx, uintptr_t needle) {
  return UpperBoundConstrainedNonemptyW(sorted_w_arr, start_idx, end_idx, needle) - 1;
}

// Requires sorted_w_arr[start_idx] <= needle, and needle < maxint.
HEADER_INLINE uintptr_t ExpsearchLastLeqU32(const uint32_t* sorted_u32_arr, uintptr_t start_idx, uintptr_t end_idx, uint32_t needle) {
  return ExpsearchU32(sorted_u32_arr, start_idx, end_idx, needle + 1) - 1;
}

HEADER_INLINE uintptr_t ExpsearchLastLeqW(const uintptr_t* sorted_w_arr, uintptr_t start_idx, uintptr_t end_idx, uintptr_t needle) {
  return ExpsearchW(sorted_w_arr, start_idx, end_idx, needle + 1) - 1;
}

uint32_t IdxToUidx(const uintptr_t* bitvec, const uint32_t* cumulative_popcounts, uint32_t widx_start, uint32_t widx_end, uint32_t idx);

uintptr_t IdxToUidxW(const uintptr_t* bitvec, const uintptr_t* cumulative_popcounts, uintptr_t widx_start, uintptr_t widx_end, uintptr_t idx);

uint32_t ExpsearchIdxToUidx(const uintptr_t* bitvec, const uint32_t* cumulative_popcounts, uint32_t widx_end, uint32_t idx, uint32_t* widx_startp);

uintptr_t ExpsearchIdxToUidxW(const uintptr_t* bitvec, const uintptr_t* cumulative_popcounts, uintptr_t widx_end, uintptr_t idx, uintptr_t* widx_startp);


HEADER_INLINE uint32_t IsCmdlineFlag(const char* param) {
  unsigned char ucc = param[1];
  return ((*param == '-') && ((ucc > '9') || ((ucc < '0') && (ucc != '.') && (ucc != '\0'))));
}

HEADER_INLINE CXXCONST_CP IsCmdlineFlagStart(const char* param) {
  unsigned char ucc = param[1];
  if ((*param == '-') && ((ucc > '9') || ((ucc < '0') && (ucc != '.') && (ucc != '\0')))) {
    return S_CAST(CXXCONST_CP, &(param[1 + (ucc == '-')]));
  }
  return nullptr;
}

#ifdef __cplusplus
HEADER_INLINE char* IsCmdlineFlagStart(char* param) {
  return const_cast<char*>(IsCmdlineFlagStart(const_cast<const char*>(param)));
}
#endif

uint32_t GetParamCt(const char* const* argvk, uint32_t argc, uint32_t flag_idx);

BoolErr EnforceParamCtRange(const char* flag_name, uint32_t param_ct, uint32_t min_ct, uint32_t max_ct);

// PglErr SortCmdlineFlags(uint32_t max_flag_blen, uint32_t flag_ct, char* flag_buf, uint32_t* flag_map);

BoolErr CleanupLogfile(uint32_t print_end_time);

CONSTI32(kNonBigstackMin, 67108864);

CONSTI32(kBigstackMinMib, 640);
CONSTI32(kBigstackDefaultMib, 2048);

static const double kLnPvalError = 9.0;

static const double kDblNormalMin = 2.2250738585072013e-308;

// probably time to flip arena_alloc and bigstack_alloc definitions...

PglErr InitLogfile(uint32_t always_stderr, char* outname, char* outname_end);

BoolErr CleanupLogfile(uint32_t print_end_time);

// manually managed, very large double-ended stack
extern unsigned char* g_bigstack_base;
extern unsigned char* g_bigstack_end;

uintptr_t GetDefaultAllocMib();

// caller is responsible for freeing bigstack_ua
PglErr InitBigstack(uintptr_t malloc_size_mib, uintptr_t* malloc_mib_final_ptr, unsigned char** bigstack_ua_ptr);


HEADER_INLINE uintptr_t bigstack_left() {
  return g_bigstack_end - g_bigstack_base;
}

HEADER_INLINE void* bigstack_alloc_raw(uintptr_t size) {
  // Assumes caller has already forced size to a multiple of
  // kCacheline, and verified that enough space is available.
  assert(!(size % kCacheline));
  unsigned char* alloc_ptr = g_bigstack_base;
  g_bigstack_base += size;
  return alloc_ptr;
}

HEADER_INLINE void* bigstack_alloc_raw_rd(uintptr_t size) {
  // Same as bigstack_alloc_raw(), except for rounding up size.
  unsigned char* alloc_ptr = g_bigstack_base;
  g_bigstack_base += RoundUpPow2(size, kCacheline);
  return alloc_ptr;
}

// Basic 64-byte-aligned allocation at bottom of stack.
// Note that --make-pgen switches gracefully to a less memory-hungry algorithm
// when it encounters an allocation failure with its default algorithm.  Since
// it can only happen once, unlikely() is still justified, but keep an eye on
// this.
HEADER_INLINE void* bigstack_alloc(uintptr_t size) {
  size = RoundUpPow2(size, kCacheline);
  if (unlikely(bigstack_left() < size)) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return bigstack_alloc_raw(size);
}

// When we only want vector-alignment.
HEADER_INLINE void* bigstack_allocv(uintptr_t size) {
  size = RoundUpPow2(size, kCacheline);
  if (unlikely(bigstack_left() < size)) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return bigstack_alloc_raw(size);
}

// Don't want uint64_t -> uintptr_t truncation to hide insufficient-memory
// scenarios in 32-bit builds.
HEADER_INLINE void* bigstack_alloc64(uint64_t size) {
  size = RoundUpPow2U64(size, kCacheline);
  if (unlikely(bigstack_left() < size)) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return bigstack_alloc_raw(size);
}


// Typesafe, return-0-iff-success interfaces.  (See also bigstack_calloc_...
// further below.)
// This interface deliberately not provided for bigstack_alloc_raw() and
// bigstack_alloc_raw_rd() (and I wouldn't blame a future maintainer for
// entirely eliminating those functions).
HEADER_INLINE BoolErr bigstack_alloc_c(uintptr_t ct, char** c_arr_ptr) {
  *c_arr_ptr = S_CAST(char*, bigstack_alloc(ct));
  return !(*c_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_d(uintptr_t ct, double** d_arr_ptr) {
  *d_arr_ptr = S_CAST(double*, bigstack_alloc(ct * sizeof(double)));
  return !(*d_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_f(uintptr_t ct, float** f_arr_ptr) {
  *f_arr_ptr = S_CAST(float*, bigstack_alloc(ct * sizeof(float)));
  return !(*f_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_i8(uintptr_t ct, int8_t** i8_arr_ptr) {
  *i8_arr_ptr = S_CAST(int8_t*, bigstack_alloc(ct * sizeof(int8_t)));
  return !(*i8_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_i16(uintptr_t ct, int16_t** i16_arr_ptr) {
  *i16_arr_ptr = S_CAST(int16_t*, bigstack_alloc(ct * sizeof(int16_t)));
  return !(*i16_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_i32(uintptr_t ct, int32_t** i32_arr_ptr) {
  *i32_arr_ptr = S_CAST(int32_t*, bigstack_alloc(ct * sizeof(int32_t)));
  return !(*i32_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = S_CAST(unsigned char*, bigstack_alloc(ct));
  return !(*uc_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_u16(uintptr_t ct, uint16_t** u16_arr_ptr) {
  *u16_arr_ptr = S_CAST(uint16_t*, bigstack_alloc(ct * sizeof(int16_t)));
  return !(*u16_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_u32(uintptr_t ct, uint32_t** u32_arr_ptr) {
  *u32_arr_ptr = S_CAST(uint32_t*, bigstack_alloc(ct * sizeof(int32_t)));
  return !(*u32_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_w(uintptr_t ct, uintptr_t** w_arr_ptr) {
  *w_arr_ptr = S_CAST(uintptr_t*, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*w_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_hw(uintptr_t ct, Halfword** hw_arr_ptr) {
  *hw_arr_ptr = S_CAST(Halfword*, bigstack_alloc(ct * sizeof(Halfword)));
  return !(*hw_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_i64(uintptr_t ct, int64_t** i64_arr_ptr) {
  *i64_arr_ptr = S_CAST(int64_t*, bigstack_alloc(ct * sizeof(int64_t)));
  return !(*i64_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_u64(uintptr_t ct, uint64_t** u64_arr_ptr) {
  *u64_arr_ptr = S_CAST(uint64_t*, bigstack_alloc(ct * sizeof(int64_t)));
  return !(*u64_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_v(uintptr_t ct, VecW** v_arr_ptr) {
  *v_arr_ptr = S_CAST(VecW*, bigstack_alloc(ct * sizeof(VecW)));
  return !(*v_arr_ptr);
}

// some versions of gcc give aliasing warnings if we use bigstack_alloc_w()
// for everything
// if sizeof(intptr_t) != sizeof(uintptr_t*), we're doomed anyway, so I won't
// bother with that static assert...
HEADER_INLINE BoolErr bigstack_alloc_wp(uintptr_t ct, uintptr_t*** ulp_arr_ptr) {
  *ulp_arr_ptr = S_CAST(uintptr_t**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*ulp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_cp(uintptr_t ct, char*** cp_arr_ptr) {
  *cp_arr_ptr = S_CAST(char**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*cp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_kcp(uintptr_t ct, const char*** kcp_arr_ptr) {
  *kcp_arr_ptr = S_CAST(const char**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*kcp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_kcpp(uintptr_t ct, const char**** kcpp_arr_ptr) {
  *kcpp_arr_ptr = S_CAST(const char***, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*kcpp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_i16p(uintptr_t ct, int16_t*** i16p_arr_ptr) {
  *i16p_arr_ptr = S_CAST(int16_t**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*i16p_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_i32p(uintptr_t ct, int32_t*** i32p_arr_ptr) {
  *i32p_arr_ptr = S_CAST(int32_t**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*i32p_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_ucp(uintptr_t ct, unsigned char*** ucp_arr_ptr) {
  *ucp_arr_ptr = S_CAST(unsigned char**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*ucp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_u16p(uintptr_t ct, uint16_t*** u16p_arr_ptr) {
  *u16p_arr_ptr = S_CAST(uint16_t**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*u16p_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_u32p(uintptr_t ct, uint32_t*** u32p_arr_ptr) {
  *u32p_arr_ptr = S_CAST(uint32_t**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*u32p_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_u64p(uintptr_t ct, uint64_t*** u64p_arr_ptr) {
  *u64p_arr_ptr = S_CAST(uint64_t**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*u64p_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_fp(uintptr_t ct, float*** fp_arr_ptr) {
  *fp_arr_ptr = S_CAST(float**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*fp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_dp(uintptr_t ct, double*** dp_arr_ptr) {
  *dp_arr_ptr = S_CAST(double**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*dp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_vp(uintptr_t ct, VecW*** vp_arr_ptr) {
  *vp_arr_ptr = S_CAST(VecW**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*vp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_wpp(uintptr_t ct, uintptr_t**** wpp_arr_ptr) {
  *wpp_arr_ptr = S_CAST(uintptr_t***, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*wpp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_cpp(uintptr_t ct, char**** cpp_arr_ptr) {
  *cpp_arr_ptr = S_CAST(char***, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*cpp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_u32pp(uintptr_t ct, uint32_t**** u32pp_arr_ptr) {
  *u32pp_arr_ptr = S_CAST(uint32_t***, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*u32pp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_vpp(uintptr_t ct, VecW**** vpp_arr_ptr) {
  *vpp_arr_ptr = S_CAST(VecW***, bigstack_alloc(ct * sizeof(VecW)));
  return !(*vpp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_allocv_w(uintptr_t ct, uintptr_t** w_arr_ptr) {
  *w_arr_ptr = S_CAST(uintptr_t*, bigstack_allocv(ct * sizeof(intptr_t)));
  return !(*w_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc64_c(uint64_t ct, char** c_arr_ptr) {
  *c_arr_ptr = S_CAST(char*, bigstack_alloc64(ct));
  return !(*c_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc64_d(uint64_t ct, double** d_arr_ptr) {
  *d_arr_ptr = S_CAST(double*, bigstack_alloc64(ct * sizeof(double)));
  return !(*d_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc64_f(uint64_t ct, float** f_arr_ptr) {
  *f_arr_ptr = S_CAST(float*, bigstack_alloc64(ct * sizeof(float)));
  return !(*f_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc64_uc(uint64_t ct, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = S_CAST(unsigned char*, bigstack_alloc64(ct));
  return !(*uc_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc64_u32(uint64_t ct, uint32_t** u32_arr_ptr) {
  *u32_arr_ptr = S_CAST(uint32_t*, bigstack_alloc64(ct * sizeof(int32_t)));
  return !(*u32_arr_ptr);
}

typedef struct LlStrStruct {
  NONCOPYABLE(LlStrStruct);
  struct LlStrStruct* next;
  char str[];
} LlStr;

// bigstack_end_alloc_llstr is usually better for small allocations like this.
HEADER_INLINE BoolErr bigstack_alloc_llstr(uintptr_t str_blen, LlStr** llstr_arr_ptr) {
  *llstr_arr_ptr = S_CAST(LlStr*, bigstack_alloc(str_blen + sizeof(LlStr)));
  return !(*llstr_arr_ptr);
}

BoolErr bigstack_calloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr);

BoolErr bigstack_calloc_d(uintptr_t ct, double** d_arr_ptr);

BoolErr bigstack_calloc_f(uintptr_t ct, float** f_arr_ptr);

BoolErr bigstack_calloc_u16(uintptr_t ct, uint16_t** u16_arr_ptr);

BoolErr bigstack_calloc_u32(uintptr_t ct, uint32_t** u32_arr_ptr);

BoolErr bigstack_calloc_w(uintptr_t ct, uintptr_t** w_arr_ptr);

BoolErr bigstack_calloc_u64(uintptr_t ct, uint64_t** u64_arr_ptr);

BoolErr bigstack_calloc_v(uintptr_t ct, VecW** v_arr_ptr);

BoolErr bigstack_calloc_wp(uintptr_t ct, uintptr_t*** wp_arr_ptr);

BoolErr bigstack_calloc_cp(uintptr_t ct, char*** cp_arr_ptr);

BoolErr bigstack_calloc_kcp(uintptr_t ct, const char*** kcp_arr_ptr);

BoolErr bigstack_calloc_cpp(uintptr_t ct, char**** cpp_arr_ptr);

#ifdef __LP64__
HEADER_INLINE BoolErr bigstack_calloc64_d(uint64_t ct, double** d_arr_ptr) {
  return bigstack_calloc_d(ct, d_arr_ptr);
}
#else
BoolErr bigstack_calloc64_d(uint64_t ct, double** d_arr_ptr);
#endif

#ifdef CPP11_TYPE_ENFORCEMENT

template <class T> BoolErr BigstackAllocX(uintptr_t ct, T** x_arr_ptr) {
  *x_arr_ptr = S_CAST(T*, bigstack_alloc(ct * sizeof(T)));
  return !(*x_arr_ptr);
}

// todo: define all BigstackAlloc functions in terms of ArenaAlloc; then these
// can be namespaced, and we only need ARENA_ALLOC_X and ARENA_ALLOC_STD_ARRAY
// macros
#  define BIGSTACK_ALLOC_X(tt, ct, pp) plink2::BigstackAllocX<tt>((ct), (pp))

#  define BIGSTACK_ALLOC_STD_ARRAY(tt, arr_size, len, pp) plink2::BigstackAllocX<std::array<tt, arr_size>>((len), (pp))

#else

#  define BIGSTACK_ALLOC_X(tt, ct, pp) (!((*(pp)) = S_CAST(tt*, bigstack_alloc((ct) * sizeof(tt)))))

#  define BIGSTACK_ALLOC_STD_ARRAY(tt, arr_size, len, pp) (!((*(pp)) = S_CAST(STD_ARRAY_PTR_TYPE(tt, arr_size), bigstack_alloc((len) * (arr_size * sizeof(tt))))))

#endif

HEADER_INLINE void BigstackReset(void* new_base) {
  g_bigstack_base = S_CAST(unsigned char*, new_base);
}

HEADER_INLINE void BigstackEndReset(void* new_end) {
  g_bigstack_end = S_CAST(unsigned char*, new_end);
}

HEADER_INLINE void BigstackDoubleReset(void* new_base, void* new_end) {
  BigstackReset(new_base);
  BigstackEndReset(new_end);
}

// assumes we've already been writing to wptr and have previously performed
// bounds-checking.
HEADER_INLINE void BigstackFinalizeW(__maybe_unused const uintptr_t* wptr, uintptr_t ct) {
  assert(wptr == R_CAST(const uintptr_t*, g_bigstack_base));
  g_bigstack_base += RoundUpPow2(ct * sizeof(intptr_t), kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}

HEADER_INLINE void BigstackFinalizeU32(__maybe_unused const uint32_t* u32ptr, uintptr_t ct) {
  assert(u32ptr == R_CAST(const uint32_t*, g_bigstack_base));
  g_bigstack_base += RoundUpPow2(ct * sizeof(int32_t), kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}

HEADER_INLINE void BigstackFinalizeU64(__maybe_unused const uint64_t* u64ptr, uintptr_t ct) {
  assert(u64ptr == R_CAST(const uint64_t*, g_bigstack_base));
  g_bigstack_base += RoundUpPow2(ct * sizeof(int64_t), kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}

HEADER_INLINE void BigstackFinalizeC(__maybe_unused const char* cptr, uintptr_t ct) {
  assert(cptr == R_CAST(const char*, g_bigstack_base));
  g_bigstack_base += RoundUpPow2(ct, kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}

HEADER_INLINE void BigstackFinalizeCp(__maybe_unused const char* const* cpptr, uintptr_t ct) {
  assert(cpptr == R_CAST(const char* const*, g_bigstack_base));
  g_bigstack_base += RoundUpPow2(ct * sizeof(intptr_t), kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}


HEADER_INLINE void BigstackBaseSet(const void* unaligned_base) {
  g_bigstack_base = R_CAST(unsigned char*, RoundUpPow2(R_CAST(uintptr_t, unaligned_base), kCacheline));
}

// When using BigstackBaseSet() after a loop where tmp_alloc_end is fixed, the
// latter should be initialized with BigstackEndRoundedDown().
// If tmp_alloc_end is not fixed, BigstackBaseSetChecked() should be called
// instead.
HEADER_INLINE BoolErr BigstackBaseSetChecked(const void* unaligned_base) {
  g_bigstack_base = R_CAST(unsigned char*, RoundUpPow2(R_CAST(uintptr_t, unaligned_base), kCacheline));
  return (g_bigstack_base > g_bigstack_end);
}

HEADER_INLINE unsigned char* BigstackEndRoundedDown() {
  return R_CAST(unsigned char*, RoundDownPow2(R_CAST(uintptr_t, g_bigstack_end), kCacheline));
}

HEADER_INLINE void BigstackShrinkTop(const void* rebase, uintptr_t new_size) {
  // could assert that this doesn't go in the wrong direction?
  g_bigstack_base = R_CAST(unsigned char*, RoundUpPow2(R_CAST(uintptr_t, rebase) + new_size, kCacheline));
}

HEADER_INLINE void BigstackEndSet(const void* unaligned_end) {
  g_bigstack_end = R_CAST(unsigned char*, RoundDownPow2(R_CAST(uintptr_t, unaligned_end), kEndAllocAlign));
}

// assumes size is divisible by kEndAllocAlign
// assumes enough space is available
HEADER_INLINE void* bigstack_end_alloc_raw(uintptr_t size) {
  assert(!(size % kEndAllocAlign));
  g_bigstack_end -= size;
  return g_bigstack_end;
}

HEADER_INLINE void* bigstack_end_alloc_raw_rd(uintptr_t size) {
  g_bigstack_end -= RoundUpPow2(size, kEndAllocAlign);
  return g_bigstack_end;
}

HEADER_INLINE void* bigstack_end_alloc_presized(uintptr_t size) {
  const uintptr_t cur_bigstack_left = bigstack_left();
  if (size > cur_bigstack_left) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return bigstack_end_alloc_raw(size);
}

HEADER_INLINE void* bigstack_end_alloc(uintptr_t size) {
  size = RoundUpPow2(size, kEndAllocAlign);
  return bigstack_end_alloc_presized(size);
}

HEADER_INLINE void* bigstack_end_alloc64_presized(uint64_t size) {
  const uintptr_t cur_bigstack_left = bigstack_left();
  if (size > cur_bigstack_left) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return bigstack_end_alloc_raw(size);
}

HEADER_INLINE void* bigstack_end_alloc64(uint64_t size) {
  size = RoundUpPow2U64(size, kEndAllocAlign);
  return bigstack_end_alloc64_presized(size);
}

HEADER_INLINE BoolErr bigstack_end_alloc_c(uintptr_t ct, char** c_arr_ptr) {
  *c_arr_ptr = S_CAST(char*, bigstack_end_alloc(ct));
  return !(*c_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_d(uintptr_t ct, double** d_arr_ptr) {
  *d_arr_ptr = S_CAST(double*, bigstack_end_alloc(ct * sizeof(double)));
  return !(*d_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_f(uintptr_t ct, float** f_arr_ptr) {
  *f_arr_ptr = S_CAST(float*, bigstack_end_alloc(ct * sizeof(float)));
  return !(*f_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_i32(uintptr_t ct, int32_t** i32_arr_ptr) {
  *i32_arr_ptr = S_CAST(int32_t*, bigstack_end_alloc(ct * sizeof(int32_t)));
  return !(*i32_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = S_CAST(unsigned char*, bigstack_end_alloc(ct));
  return !(*uc_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_u32(uintptr_t ct, uint32_t** u32_arr_ptr) {
  *u32_arr_ptr = S_CAST(uint32_t*, bigstack_end_alloc(ct * sizeof(int32_t)));
  return !(*u32_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_w(uintptr_t ct, uintptr_t** w_arr_ptr) {
  *w_arr_ptr = S_CAST(uintptr_t*, bigstack_end_alloc(ct * sizeof(intptr_t)));
  return !(*w_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_i64(uintptr_t ct, int64_t** i64_arr_ptr) {
  *i64_arr_ptr = S_CAST(int64_t*, bigstack_end_alloc(ct * sizeof(int64_t)));
  return !(*i64_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_u64(uintptr_t ct, uint64_t** u64_arr_ptr) {
  *u64_arr_ptr = S_CAST(uint64_t*, bigstack_end_alloc(ct * sizeof(int64_t)));
  return !(*u64_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_llstr(uintptr_t str_blen, LlStr** llstr_arr_ptr) {
  *llstr_arr_ptr = S_CAST(LlStr*, bigstack_end_alloc(str_blen + sizeof(LlStr)));
  return !(*llstr_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_wp(uintptr_t ct, uintptr_t*** wp_arr_ptr) {
  *wp_arr_ptr = S_CAST(uintptr_t**, bigstack_end_alloc(ct * sizeof(intptr_t)));
  return !(*wp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_cp(uintptr_t ct, char*** cp_arr_ptr) {
  *cp_arr_ptr = S_CAST(char**, bigstack_end_alloc(ct * sizeof(intptr_t)));
  return !(*cp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_kcp(uintptr_t ct, const char*** kcp_arr_ptr) {
  *kcp_arr_ptr = S_CAST(const char**, bigstack_end_alloc(ct * sizeof(intptr_t)));
  return !(*kcp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_ucp(uintptr_t ct, unsigned char*** ucp_arr_ptr) {
  *ucp_arr_ptr = S_CAST(unsigned char**, bigstack_end_alloc(ct * sizeof(intptr_t)));
  return !(*ucp_arr_ptr);
}

BoolErr bigstack_end_calloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr);

BoolErr bigstack_end_calloc_d(uintptr_t ct, double** d_arr_ptr);

BoolErr bigstack_end_calloc_f(uintptr_t ct, float** f_arr_ptr);

BoolErr bigstack_end_calloc_u32(uintptr_t ct, uint32_t** u32_arr_ptr);

BoolErr bigstack_end_calloc_w(uintptr_t ct, uintptr_t** w_arr_ptr);

BoolErr bigstack_end_calloc_u64(uintptr_t ct, uint64_t** u64_arr_ptr);

BoolErr bigstack_end_calloc_wp(uintptr_t ct, uintptr_t*** wp_arr_ptr);

BoolErr bigstack_end_calloc_cp(uintptr_t ct, char*** cp_arr_ptr);

// Not sure if this is needed.
#if (__GNUC__ >= 11) && !defined(__clang__)
#  pragma GCC push_options
#  pragma GCC optimize("-fno-ipa-modref")
#endif
HEADER_INLINE BoolErr bigstack_calloc_c(uintptr_t ct, char** c_arr_ptr) {
  return bigstack_calloc_uc(ct, R_CAST(unsigned char**, c_arr_ptr));
}

HEADER_INLINE BoolErr bigstack_calloc_i16(uintptr_t ct, int16_t** i16_arr_ptr) {
  return bigstack_calloc_u16(ct, R_CAST(uint16_t**, i16_arr_ptr));
}

HEADER_INLINE BoolErr bigstack_calloc_i32(uintptr_t ct, int32_t** i32_arr_ptr) {
  return bigstack_calloc_u32(ct, R_CAST(uint32_t**, i32_arr_ptr));
}

HEADER_INLINE BoolErr bigstack_calloc_i64(uintptr_t ct, int64_t** i64_arr_ptr) {
  return bigstack_calloc_u64(ct, R_CAST(uint64_t**, i64_arr_ptr));
}

HEADER_INLINE BoolErr bigstack_calloc_u32p(uintptr_t ct, uint32_t*** u32p_arr_ptr) {
  return bigstack_calloc_w(ct, R_CAST(uintptr_t**, u32p_arr_ptr));
}

HEADER_INLINE BoolErr bigstack_end_calloc_c(uintptr_t ct, char** c_arr_ptr) {
  return bigstack_end_calloc_uc(ct, R_CAST(unsigned char**, c_arr_ptr));
}

HEADER_INLINE BoolErr bigstack_end_calloc_i32(uintptr_t ct, int32_t** i32_arr_ptr) {
  return bigstack_end_calloc_u32(ct, R_CAST(uint32_t**, i32_arr_ptr));
}

HEADER_INLINE BoolErr bigstack_end_calloc_i64(uintptr_t ct, int64_t** i64_arr_ptr) {
  return bigstack_end_calloc_u64(ct, R_CAST(uint64_t**, i64_arr_ptr));
}

HEADER_INLINE BoolErr bigstack_end_calloc_kcp(uintptr_t ct, const char*** kcp_arr_ptr) {
  return bigstack_end_calloc_cp(ct, K_CAST(char***, kcp_arr_ptr));
}
#if (__GNUC__ >= 11) && !defined(__clang__)
#  pragma GCC pop_options
#endif

#ifdef __LP64__
HEADER_INLINE BoolErr bigstack_end_calloc64_w(uint64_t ct, uintptr_t** w_arr_ptr) {
  return bigstack_end_calloc_w(ct, w_arr_ptr);
}
#else
BoolErr bigstack_end_calloc64_w(uint64_t ct, uintptr_t** w_arr_ptr);
#endif

HEADER_INLINE void bigstack_end_clalign() {
  g_bigstack_end = BigstackEndRoundedDown();
}

HEADER_INLINE void* bigstack_end_clalloc(uintptr_t size) {
  // assumes g_bigstack_end is cacheline-aligned by e.g. a
  // bigstack_end_clalign() call.
  size = RoundUpPow2(size, kCacheline);
  return bigstack_end_alloc_presized(size);
}

HEADER_INLINE BoolErr bigstack_end_clalloc_c(uintptr_t ct, char** c_arr_ptr) {
  *c_arr_ptr = S_CAST(char*, bigstack_end_clalloc(ct));
  return !(*c_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_clalloc_d(uintptr_t ct, double** d_arr_ptr) {
  *d_arr_ptr = S_CAST(double*, bigstack_end_clalloc(ct * sizeof(double)));
  return !(*d_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_clalloc_v(uintptr_t ct, VecW** v_arr_ptr) {
  *v_arr_ptr = S_CAST(VecW*, bigstack_end_clalloc(ct * sizeof(VecW)));
  return !(*v_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_clalloc_cp(uintptr_t ct, char*** cp_arr_ptr) {
  *cp_arr_ptr = S_CAST(char**, bigstack_end_clalloc(ct * sizeof(intptr_t)));
  return !(*cp_arr_ptr);
}

// most of the non-global arena interface (necessary for some multithreaded
// code) now in include/plink2_memory; arena_alloc_raw() and
// arena_alloc_raw_rd() in include/plink2_base


BoolErr PushLlStr(const char* str, LlStr** ll_stack_ptr);

// Does not require null-termination.
// BoolErr push_llstr_counted(const char* str, uint32_t slen, LlStr** ll_stack_ptr);

void llstr_free_cond(LlStr* llstr_head);

typedef struct L32StrStruct {
  NONCOPYABLE(L32StrStruct);
  uint32_t len;
  char str[];
} L32Str;

// assumes multistr is nonempty
// BoolErr CountAndMeasureMultistrAlloc(const char* multistr, uintptr_t max_str_ct, const char** strptr_arr, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr);

// assumes multistr is nonempty
BoolErr CountAndMeasureMultistrReverseAlloc(const char* multistr, uintptr_t max_str_ct, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr, const char*** strptr_arrp);

// er, should benchmark this against hash table and InMultistr()
BoolErr MultistrToStrboxDedupArenaAlloc(unsigned char* arena_top, const char* multistr, unsigned char** arena_bottom_ptr, char** sorted_strbox_ptr, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr);

HEADER_INLINE BoolErr MultistrToStrboxDedupAlloc(const char* multistr, char** sorted_strbox_ptr, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr) {
  return MultistrToStrboxDedupArenaAlloc(g_bigstack_end, multistr, &g_bigstack_base, sorted_strbox_ptr, str_ct_ptr, max_blen_ptr);
}


void DivisionMagicNums(uint32_t divisor, uint64_t* multp, uint32_t* __restrict pre_shiftp, uint32_t* __restrict post_shiftp, uint32_t* __restrict incrp);

// ZeroU32Arr, ZeroWArr, ZeroU64Arr, SetAllWArr, SetAllU32Arr currently defined
// in plink2_base.h

HEADER_INLINE void ZeroVecArr(uintptr_t entry_ct, VecW* vvec) {
  memset(vvec, 0, entry_ct * kBytesPerVec);
}

HEADER_INLINE void SetAllU64Arr(uintptr_t entry_ct, uint64_t* u64arr) {
  // bugfix (1 Feb 2018): forgot to multiply by 2 in 32-bit case
  SetAllWArr(entry_ct * (sizeof(int64_t) / kBytesPerWord), R_CAST(uintptr_t*, u64arr));
}

HEADER_INLINE void ZeroI32Arr(uintptr_t entry_ct, int32_t* i32arr) {
  memset(i32arr, 0, entry_ct * sizeof(int32_t));
}

HEADER_INLINE void SetAllI32Arr(uintptr_t entry_ct, int32_t* i32arr) {
  for (uintptr_t ulii = 0; ulii != entry_ct; ulii++) {
    *i32arr++ = -1;
  }
}

HEADER_INLINE void ZeromovFArr(uintptr_t entry_ct, float** farr_ptr) {
  float* farr = *farr_ptr;
  for (uintptr_t ulii = 0; ulii != entry_ct; ulii++) {
    *farr++ = 0.0;
  }
  *farr_ptr = farr;
}

HEADER_INLINE void ZeromovDArr(uintptr_t entry_ct, double** darr_ptr) {
  double* darr = *darr_ptr;
  for (uintptr_t ulii = 0; ulii != entry_ct; ulii++) {
    *darr++ = 0.0;
  }
  *darr_ptr = darr;
}


// SetAllBits, IsSet, SetBit, ClearBit, AdvTo1Bit, AdvTo0Bit, AdvBoundedTo1Bit,
// FindLast1BitBefore, AllWordsAreZero, FillBitsNz, ClearBitsNz defined in
// plink2_bits.h

// Useful when we don't want to think about the signedness of a 32-bit int.
HEADER_INLINE void SetBitI(int32_t loc, uintptr_t* bitarr) {
  SetBit(S_CAST(uint32_t, loc), bitarr);
}

HEADER_INLINE void FlipBit(uintptr_t loc, uintptr_t* bitarr) {
  bitarr[loc / kBitsPerWord] ^= k1LU << (loc % kBitsPerWord);
}

// floor permitted to be -1, though not smaller than that.
int32_t FindLast1BitBeforeBounded(const uintptr_t* bitarr, uint32_t loc, int32_t floor);

// This can be made a tiny bit faster than memequal() in isolation, but the
// difference is almost certainly too small to justify additional i-cache
// pressure.
HEADER_INLINE uint32_t wordsequal(const uintptr_t* word_arr1, const uintptr_t* word_arr2, uintptr_t word_ct) {
  return memequal(word_arr1, word_arr2, word_ct * kBytesPerWord);
}


// BitvecAnd(), BitvecInvmask(), BitvecOr(), BitvecInvert(), BitvecXorCopy(),
// BitvecInvertCopy(), AlignedBitarrInvert(), and AlignedBitarrInvertCopy() in
// plink2_bits.h

void BitvecAndCopy(const uintptr_t* __restrict source1_bitvec, const uintptr_t* __restrict source2_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec);

void BitvecInvmaskCopy(const uintptr_t* __restrict source_bitvec, const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec);

void BitvecXor(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec);

void BitvecInvertAndMask(const uintptr_t* __restrict include_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

// void BitvecOrNot(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec);

// yes || (!no)
// 'AlignedBitarr' instead of Bitvec since this takes bit_ct instead of word_ct
// as the size argument, and zeroes trailing bits.
void AlignedBitarrOrnotCopy(const uintptr_t* __restrict argyes_bitvec, const uintptr_t* __restrict argno_bitvec, uintptr_t bit_ct, uintptr_t* __restrict target_bitvec);

void BitvecXor3Copy(const uintptr_t* __restrict src_bitvec1, const uintptr_t* __restrict src_bitvec2, const uintptr_t* __restrict src_bitvec3, uintptr_t word_ct, uintptr_t* __restrict target_bitvec);

// Address C-only incompatible-pointer-types-discards-qualifiers warning.
#ifdef __cplusplus
#  define TO_CONSTU32PCONSTP(u32_pp) (u32_pp)
#else
#  define TO_CONSTU32PCONSTP(u32_pp) ((const uint32_t* const*)(u32_pp))
#endif

HEADER_INLINE void U32VecAdd(const VecU32* arg_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] += arg_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecAdd(const uint32_t* arg_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecAdd(R_CAST(const VecU32*, arg_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

HEADER_INLINE void I32CastVecAdd(const int32_t* arg_i32arr, uintptr_t vec_ct, int32_t* main_i32arr) {
  U32VecAdd(R_CAST(const VecU32*, arg_i32arr), vec_ct, R_CAST(VecU32*, main_i32arr));
}

HEADER_INLINE void U32VecSub(const VecU32* arg_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] -= arg_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecSub(const uint32_t* arg_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecSub(R_CAST(const VecU32*, arg_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

HEADER_INLINE void U32VecMaskedAdd(const VecU32* mask_u32vec, const VecU32* arg_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] += mask_u32vec[vidx] & arg_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecMaskedAdd(const uint32_t* mask_u32arr, const uint32_t* arg_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecMaskedAdd(R_CAST(const VecU32*, mask_u32arr), R_CAST(const VecU32*, arg_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

HEADER_INLINE void U32VecInvmaskedAdd(const VecU32* invmask_u32vec, const VecU32* arg_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] += vecu32_and_notfirst(invmask_u32vec[vidx], arg_u32vec[vidx]);
  }
}

HEADER_INLINE void U32CastVecInvmaskedAdd(const uint32_t* invmask_u32arr, const uint32_t* arg_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecInvmaskedAdd(R_CAST(const VecU32*, invmask_u32arr), R_CAST(const VecU32*, arg_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}


HEADER_INLINE void U32VecAdd2(const VecU32* arg1_u32vec, const VecU32* arg2_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] += arg1_u32vec[vidx] + arg2_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecAdd2(const uint32_t* arg1_u32arr, const uint32_t* arg2_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecAdd2(R_CAST(const VecU32*, arg1_u32arr), R_CAST(const VecU32*, arg2_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

HEADER_INLINE void U32VecAssignAdd2(const VecU32* arg1_u32vec, const VecU32* arg2_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] = arg1_u32vec[vidx] + arg2_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecAssignAdd2(const uint32_t* arg1_u32arr, const uint32_t* arg2_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecAssignAdd2(R_CAST(const VecU32*, arg1_u32arr), R_CAST(const VecU32*, arg2_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

HEADER_INLINE void U32VecMaskedAdd2(const VecU32* mask_u32vec, const VecU32* arg1_u32vec, const VecU32* arg2_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] += mask_u32vec[vidx] & (arg1_u32vec[vidx] + arg2_u32vec[vidx]);
  }
}

HEADER_INLINE void U32CastVecMaskedAdd2(const uint32_t* mask_u32arr, const uint32_t* arg1_u32arr, const uint32_t* arg2_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecMaskedAdd2(R_CAST(const VecU32*, mask_u32arr), R_CAST(const VecU32*, arg1_u32arr), R_CAST(const VecU32*, arg2_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

HEADER_INLINE void U32VecAddSub(const VecU32* add_u32vec, const VecU32* sub_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] += add_u32vec[vidx] - sub_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecAddSub(const uint32_t* add_u32arr, const uint32_t* sub_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecAddSub(R_CAST(const VecU32*, add_u32arr), R_CAST(const VecU32*, sub_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}


HEADER_INLINE void U32VecAdd3(const VecU32* arg1_u32vec, const VecU32* arg2_u32vec, const VecU32* arg3_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] += arg1_u32vec[vidx] + arg2_u32vec[vidx] + arg3_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecAdd3(const uint32_t* arg1_u32arr, const uint32_t* arg2_u32arr, const uint32_t* arg3_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecAdd3(R_CAST(const VecU32*, arg1_u32arr), R_CAST(const VecU32*, arg2_u32arr), R_CAST(const VecU32*, arg3_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

HEADER_INLINE void U32VecAssignAdd3(const VecU32* arg1_u32vec, const VecU32* arg2_u32vec, const VecU32* arg3_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] = arg1_u32vec[vidx] + arg2_u32vec[vidx] + arg3_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecAssignAdd3(const uint32_t* arg1_u32arr, const uint32_t* arg2_u32arr, const uint32_t* arg3_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecAssignAdd3(R_CAST(const VecU32*, arg1_u32arr), R_CAST(const VecU32*, arg2_u32arr), R_CAST(const VecU32*, arg3_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

HEADER_INLINE void U32VecSub3(const VecU32* arg1_u32vec, const VecU32* arg2_u32vec, const VecU32* arg3_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] -= arg1_u32vec[vidx] + arg2_u32vec[vidx] + arg3_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecSub3(const uint32_t* arg1_u32arr, const uint32_t* arg2_u32arr, const uint32_t* arg3_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecSub3(R_CAST(const VecU32*, arg1_u32arr), R_CAST(const VecU32*, arg2_u32arr), R_CAST(const VecU32*, arg3_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

HEADER_INLINE void U32VecAssignAdd5(const VecU32* arg1_u32vec, const VecU32* arg2_u32vec, const VecU32* arg3_u32vec, const VecU32* arg4_u32vec, const VecU32* arg5_u32vec, uintptr_t vec_ct, VecU32* main_u32vec) {
  for (uintptr_t vidx = 0; vidx != vec_ct; ++vidx) {
    main_u32vec[vidx] = arg1_u32vec[vidx] + arg2_u32vec[vidx] + arg3_u32vec[vidx] + arg4_u32vec[vidx] + arg5_u32vec[vidx];
  }
}

HEADER_INLINE void U32CastVecAssignAdd5(const uint32_t* arg1_u32arr, const uint32_t* arg2_u32arr, const uint32_t* arg3_u32arr, const uint32_t* arg4_u32arr, const uint32_t* arg5_u32arr, uintptr_t vec_ct, uint32_t* main_u32arr) {
  U32VecAssignAdd5(R_CAST(const VecU32*, arg1_u32arr), R_CAST(const VecU32*, arg2_u32arr), R_CAST(const VecU32*, arg3_u32arr), R_CAST(const VecU32*, arg4_u32arr), R_CAST(const VecU32*, arg5_u32arr), vec_ct, R_CAST(VecU32*, main_u32arr));
}

// This usually makes a temporary g_bigstack allocation.
HEADER_INLINE BoolErr SortStrptrArrIndexed(uint32_t str_ct, uint32_t leave_first_alone, uint32_t overread_ok, uint32_t use_nsort, const char** strptrs, uint32_t* new_to_old_idx, uint32_t* old_to_new_idx) {
  const uint32_t str_sort_ct = str_ct - leave_first_alone;
  if ((str_sort_ct > 1) && (bigstack_left() < str_sort_ct * sizeof(StrSortIndexedDeref))) {
    return 1;
  }
  SortStrptrArrIndexed2(str_ct, leave_first_alone, overread_ok, use_nsort, strptrs, new_to_old_idx, old_to_new_idx, g_bigstack_base);
  return 0;
}

// basic linear scan
// returns -1 on failure to find, -2 if duplicate
int32_t GetVariantUidxWithoutHtable(const char* idstr, const char* const* variant_ids, const uintptr_t* variant_include, uint32_t variant_ct);

BoolErr HtableGoodSizeAlloc(uint32_t item_ct, uintptr_t bytes_avail, uint32_t** htable_ptr, uint32_t* htable_size_ptr);

// returned index in duplicate-pair case is unfiltered
// uint32_t populate_strbox_subset_htable(const uintptr_t* __restrict subset_mask, const char* strbox, uintptr_t raw_str_ct, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t str_htable_size, uint32_t* str_htable);

// assumes cur_id_slen < max_str_blen.
// requires cur_id to be null-terminated.
uint32_t StrboxHtableFind(const char* cur_id, const char* strbox, const uint32_t* id_htable, uintptr_t max_str_blen, uint32_t cur_id_slen, uint32_t id_htable_size);

uint32_t StrboxHtableFindNnt(const char* cur_id, const char* strbox, const uint32_t* id_htable, uintptr_t max_str_blen, uint32_t cur_id_slen, uint32_t id_htable_size);

// Assumes cur_id is null-terminated.
// Returns string-index if cur_id is already in the table, UINT32_MAX if it was
// added.
uint32_t StrboxHtableAdd(const char* cur_id, const char* strbox, uintptr_t max_str_blen, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t value, uint32_t* id_htable);

// useful for duplicate detection: returns 0 on no duplicates, a positive index
// of a duplicate pair if they're present
uint32_t PopulateStrboxHtable(const char* strbox, uint32_t str_ct, uintptr_t max_str_blen, uint32_t str_htable_size, uint32_t* str_htable);

uint32_t PopulateStrboxSubsetHtable(const char* strbox, const uintptr_t* subset_mask, uint32_t str_ct, uintptr_t max_str_blen, uint32_t allow_dups, uint32_t str_htable_size, uint32_t* str_htable);

// last variant_ids entry must be at least kMaxIdBlen bytes before end of
// bigstack
uint32_t VariantIdDupflagHtableFind(const char* idbuf, const char* const* variant_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen);

uint32_t VariantIdDupHtableFind(const char* idbuf, const char* const* variant_ids, const uint32_t* id_htable, const uint32_t* htable_dup_base, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen, uint32_t* llidx_ptr);


// This still perform a temporary bigstack allocation; 'noalloc' here just
// means that sorted_strbox and id_map must be allocated in advance.  (Overread
// must be safe.)
PglErr CopySortStrboxSubsetNoalloc(const uintptr_t* __restrict subset_mask, const char* __restrict orig_strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t collapse_idxs, uint32_t use_nsort, char* __restrict sorted_strbox, uint32_t* __restrict id_map);

PglErr CopySortStrboxSubset(const uintptr_t* __restrict subset_mask, const char* __restrict orig_strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t collapse_idxs, uint32_t use_nsort, char** sorted_strbox_ptr, uint32_t** id_map_ptr);


typedef struct RangeListStruct {
  NONCOPYABLE(RangeListStruct);
  char* names;
  unsigned char* starts_range;
  uint32_t name_ct;
  uint32_t name_max_blen;
} RangeList;

void InitRangeList(RangeList* range_list_ptr);

void CleanupRangeList(RangeList* range_list_ptr);

// bitarr assumed to be initialized (but not necessarily zero-initialized)
BoolErr NumericRangeListToBitarr(const RangeList* range_list_ptr, uint32_t bitarr_size, uint32_t offset, uint32_t ignore_overflow, uintptr_t* bitarr);

// returns UINT32_MAX on parsing error
uint32_t NumericRangeListMax(const RangeList* range_list_ptr);

PglErr StringRangeListToBitarr(const char* header_line, const RangeList* range_list_ptr, const char* __restrict sorted_ids, const uint32_t* __restrict id_map, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t* bitarr, int32_t* __restrict seen_idxs);

PglErr StringRangeListToBitarrAlloc(const char* header_line, const RangeList* range_list_ptr, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t** bitarr_ptr);


HEADER_INLINE uint32_t IntersectionIsEmpty(const uintptr_t* bitvec1, const uintptr_t* bitvec2, uintptr_t word_ct) {
#ifdef USE_SSE42
  const uintptr_t fullvec_ct = word_ct / kWordsPerVec;
#  ifdef USE_AVX2
  const __m256i* bitvvec1 = R_CAST(const __m256i*, bitvec1);
  const __m256i* bitvvec2 = R_CAST(const __m256i*, bitvec2);
  for (uintptr_t vidx = 0; vidx != fullvec_ct; ++vidx) {
    if (!_mm256_testz_si256(bitvvec1[vidx], bitvvec2[vidx])) {
      return 0;
    }
  }
#  else
  const __m128i* bitvvec1 = R_CAST(const __m128i*, bitvec1);
  const __m128i* bitvvec2 = R_CAST(const __m128i*, bitvec2);
  for (uintptr_t vidx = 0; vidx != fullvec_ct; ++vidx) {
    if (!_mm_testz_si128(bitvvec1[vidx], bitvvec2[vidx])) {
      return 0;
    }
  }
#  endif
  word_ct = word_ct & (kWordsPerVec - 1);
  bitvec1 = R_CAST(const uintptr_t*, &(bitvvec1[fullvec_ct]));
  bitvec2 = R_CAST(const uintptr_t*, &(bitvvec2[fullvec_ct]));
#endif
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    if (bitvec1[widx] & bitvec2[widx]) {
      return 0;
    }
  }
  return 1;
}

uint32_t IntersectionRangeIsEmpty(const uintptr_t* bitarr1, const uintptr_t* bitarr2, uintptr_t start_idx, uintptr_t end_idx);

// could use testz here, but it's more awkward
HEADER_INLINE uint32_t UnionIsFull(const uintptr_t* bitarr1, const uintptr_t* bitarr2, uintptr_t bit_ct) {
  const uintptr_t fullword_ct = bit_ct / kBitsPerWord;
  for (uintptr_t widx = 0; widx != fullword_ct; ++widx) {
    if ((bitarr1[widx] | bitarr2[widx]) != ~k0LU) {
      return 0;
    }
  }
  const uint32_t trailing_bit_ct = bit_ct % kBitsPerWord;
  if (trailing_bit_ct) {
    if ((bitarr1[fullword_ct] | bitarr2[fullword_ct]) != ((k1LU << trailing_bit_ct) - k1LU)) {
      return 0;
    }
  }
  return 1;
}

// PopcountWordsIntersect moved to plink2_bits

void PopcountWordsIntersect3val(const uintptr_t* __restrict bitvec1, const uintptr_t* __restrict bitvec2, uint32_t word_ct, uint32_t* __restrict popcount1_ptr, uint32_t* __restrict popcount2_ptr, uint32_t* __restrict popcount_intersect_ptr);

// uintptr_t count_11_longs(const uintptr_t* genovec, uintptr_t word_ct);
// (see CountNyp in pgenlib_misc for that functionality)

// Must be safe to read from bitarr[start_idx / kBitsPerWord].
uint32_t AllBitsAreZero(const uintptr_t* bitarr, uintptr_t start_idx, uintptr_t end_idx);

// does not assume relevant bits of target_bitarr are zero
HEADER_INLINE void CopyBits(uintptr_t src_word, uint64_t target_start_bitidx, uint32_t len, void* target_bitarr) {
  unsigned char* target_bitarr_uc = S_CAST(unsigned char*, target_bitarr);
  unsigned char* update_ptr = &(target_bitarr_uc[(target_start_bitidx / kBitsPerWord) * kBytesPerWord]);
  uintptr_t target_word;
  CopyFromUnalignedW(&target_word, update_ptr);
  const uint32_t bit_idx_start = target_start_bitidx % kBitsPerWord;
  const uintptr_t shifted_new_bits = src_word << bit_idx_start;
  const uint32_t bit_idx_end = bit_idx_start + len;
  if (bit_idx_end <= kBitsPerWord) {
    const uintptr_t invmask = bzhi_max((~k0LU) << bit_idx_start, bit_idx_end);
    target_word = (target_word & (~invmask)) | shifted_new_bits;
    CopyToUnalignedW(update_ptr, &target_word);
    return;
  }
  AppendW(bzhi(target_word, bit_idx_start) | shifted_new_bits, &update_ptr);

  CopyFromUnalignedW(&target_word, update_ptr);
  const uint32_t remainder = bit_idx_end - kBitsPerWord;
  target_word = (target_word & ((~k0LU) << remainder)) | (src_word >> (len - remainder));
  CopyToUnalignedW(update_ptr, &target_word);
}

// assumes len is positive, and relevant bits of target_bitarr are zero
void CopyBitarrRange(const uintptr_t* __restrict src_bitarr, uintptr_t src_start_bitidx, uintptr_t target_start_bitidx, uintptr_t len, uintptr_t* __restrict target_bitarr);


// vertical popcount support
// VcountScramble1() and friends are here since they apply to generic
// arrays; scramble_2_4_8_32() is more plink-specific

#ifdef __LP64__
#  ifdef USE_AVX2
HEADER_INLINE uint32_t VcountScramble1(uint32_t orig_idx) {
  // 1->4: 0 4 8 12 16 20 24 28 32 ... 252 1 5 9 ...
  // 4->8: 0 8 16 24 32 ... 248 4 12 20 ... 1 9 17 ...
  // 8->32: 0 32 ... 224 8 40 ... 232 ... 248 4 36 ... 252 1 33 ...
  return (orig_idx & (~255)) + ((orig_idx & 3) * 64) + ((orig_idx & 4) * 8) + (orig_idx & 24) + ((orig_idx & 224) / 32);
}
#  else
HEADER_INLINE uint32_t VcountScramble1(uint32_t orig_idx) {
  // 1->4: 0 4 8 12 16 20 24 28 32 ... 124 1 5 9 ...
  // 4->8: 0 8 16 24 32 ... 120 4 12 20 ... 1 9 17 ...
  // 8->32: 0 32 64 96 8 40 72 104 16 48 80 112 24 56 88 120 4 36 68 ... 1 33 ...
  return (orig_idx & (~127)) + ((orig_idx & 3) * 32) + ((orig_idx & 4) * 4) + ((orig_idx & 24) / 2) + ((orig_idx & 96) / 32);
}
#  endif
#else
// 1->4: 0 4 8 12 16 20 24 28 1 5 9 13 17 21 25 29 2 6 10 ...
// 4->8: 0 8 16 24 4 12 20 28 1 9 17 25 5 13 21 29 2 10 18 ...
// 8->32: 0 8 16 24 4 12 20 28 1 9 17 25 5 13 21 29 2 10 18 ...
HEADER_INLINE uint32_t VcountScramble1(uint32_t orig_idx) {
  return (orig_idx & (~31)) + ((orig_idx & 3) * 8) + (orig_idx & 4) + ((orig_idx & 24) / 8);
}
#endif

// change acc1 type to const void* if a VecW* use case arises
HEADER_INLINE void VcountIncr1To4(const uintptr_t* acc1, uint32_t acc1_vec_ct, VecW* acc4_iter) {
  const VecW m1x4 = VCONST_W(kMask1111);
  const VecW* acc1_iter = R_CAST(const VecW*, acc1);
  for (uint32_t vidx = 0; vidx != acc1_vec_ct; ++vidx) {
    VecW loader = *acc1_iter++;
    *acc4_iter = (*acc4_iter) + (loader & m1x4);
    ++acc4_iter;
    loader = vecw_srli(loader, 1);
    *acc4_iter = (*acc4_iter) + (loader & m1x4);
    ++acc4_iter;
    loader = vecw_srli(loader, 1);
    *acc4_iter = (*acc4_iter) + (loader & m1x4);
    ++acc4_iter;
    loader = vecw_srli(loader, 1);
    *acc4_iter = (*acc4_iter) + (loader & m1x4);
    ++acc4_iter;
  }
}

HEADER_INLINE void Vcount0Incr1To4(uint32_t acc1_vec_ct, uintptr_t* acc1, VecW* acc4_iter) {
  const VecW m1x4 = VCONST_W(kMask1111);
  VecW* acc1_iter = R_CAST(VecW*, acc1);
  for (uint32_t vidx = 0; vidx != acc1_vec_ct; ++vidx) {
    VecW loader = *acc1_iter;
    *acc1_iter++ = vecw_setzero();
    *acc4_iter = (*acc4_iter) + (loader & m1x4);
    ++acc4_iter;
    loader = vecw_srli(loader, 1);
    *acc4_iter = (*acc4_iter) + (loader & m1x4);
    ++acc4_iter;
    loader = vecw_srli(loader, 1);
    *acc4_iter = (*acc4_iter) + (loader & m1x4);
    ++acc4_iter;
    loader = vecw_srli(loader, 1);
    *acc4_iter = (*acc4_iter) + (loader & m1x4);
    ++acc4_iter;
  }
}

// er, should this just be the same function as unroll_incr_2_4 with extra
// parameters?...
HEADER_INLINE void VcountIncr4To8(const VecW* acc4_iter, uint32_t acc4_vec_ct, VecW* acc8_iter) {
  const VecW m4 = VCONST_W(kMask0F0F);
  for (uint32_t vidx = 0; vidx != acc4_vec_ct; ++vidx) {
    VecW loader = *acc4_iter++;
    *acc8_iter = (*acc8_iter) + (loader & m4);
    ++acc8_iter;
    loader = vecw_srli(loader, 4);
    *acc8_iter = (*acc8_iter) + (loader & m4);
    ++acc8_iter;
  }
}

HEADER_INLINE void Vcount0Incr4To8(uint32_t acc4_vec_ct, VecW* acc4_iter, VecW* acc8_iter) {
  const VecW m4 = VCONST_W(kMask0F0F);
  for (uint32_t vidx = 0; vidx != acc4_vec_ct; ++vidx) {
    VecW loader = *acc4_iter;
    *acc4_iter++ = vecw_setzero();
    *acc8_iter = (*acc8_iter) + (loader & m4);
    ++acc8_iter;
    loader = vecw_srli(loader, 4);
    *acc8_iter = (*acc8_iter) + (loader & m4);
    ++acc8_iter;
  }
}

HEADER_INLINE void VcountIncr8To32(const VecW* acc8_iter, uint32_t acc8_vec_ct, VecW* acc32_iter) {
  const VecW m8x32 = VCONST_W(kMask000000FF);
  for (uint32_t vidx = 0; vidx != acc8_vec_ct; ++vidx) {
    VecW loader = *acc8_iter++;
    *acc32_iter = (*acc32_iter) + (loader & m8x32);
    ++acc32_iter;
    loader = vecw_srli(loader, 8);
    *acc32_iter = (*acc32_iter) + (loader & m8x32);
    ++acc32_iter;
    loader = vecw_srli(loader, 8);
    *acc32_iter = (*acc32_iter) + (loader & m8x32);
    ++acc32_iter;
    loader = vecw_srli(loader, 8);
    *acc32_iter = (*acc32_iter) + (loader & m8x32);
    ++acc32_iter;
  }
}

HEADER_INLINE void Vcount0Incr8To32(uint32_t acc8_vec_ct, VecW* acc8_iter, VecW* acc32_iter) {
  const VecW m8x32 = VCONST_W(kMask000000FF);
  for (uint32_t vidx = 0; vidx != acc8_vec_ct; ++vidx) {
    VecW loader = *acc8_iter;
    *acc8_iter++ = vecw_setzero();
    *acc32_iter = (*acc32_iter) + (loader & m8x32);
    ++acc32_iter;
    loader = vecw_srli(loader, 8);
    *acc32_iter = (*acc32_iter) + (loader & m8x32);
    ++acc32_iter;
    loader = vecw_srli(loader, 8);
    *acc32_iter = (*acc32_iter) + (loader & m8x32);
    ++acc32_iter;
    loader = vecw_srli(loader, 8);
    *acc32_iter = (*acc32_iter) + (loader & m8x32);
    ++acc32_iter;
  }
}

void VerticalCounterUpdate(const uintptr_t* acc1, uint32_t acc1_vec_ct, uint32_t* rem15_and_255d15, VecW* acc4_8_32);


HEADER_INLINE uint32_t IdxToUidxBasic(const uintptr_t* bitvec, uint32_t idx) {
  return FindNth1BitFrom(bitvec, 0, idx + 1);
}

HEADER_INLINE void ComputeUidxStartPartition(const uintptr_t* variant_include, uint64_t variant_ct, uint32_t thread_ct, uint32_t first_uidx, uint32_t* variant_uidx_starts) {
  FillU32SubsetStarts(variant_include, thread_ct, first_uidx, variant_ct, variant_uidx_starts);
}

// These functions assume (bit_ct * (thread_ct - 1)) < 2^64.
// bit_ct must be positive, but can be smaller than thread_ct
// FillU32SubsetStarts moved to plink2_bits

void FillWStarts(uint32_t thread_ct, uintptr_t start, uint64_t bit_ct, uintptr_t* starts);

#ifdef __LP64__
void FillWSubsetStarts(const uintptr_t* subset, uint32_t thread_ct, uintptr_t start, uint64_t bit_ct, uintptr_t* starts);
#else
HEADER_INLINE void FillWSubsetStarts(const uintptr_t* subset, uint32_t thread_ct, uintptr_t start, uint64_t bit_ct, uintptr_t* starts) {
  return FillU32SubsetStarts(subset, thread_ct, start, bit_ct, R_CAST(uint32_t*, starts));
}
#endif

void InvertU32Arr(const uint32_t* new_idx_to_old_uidx, uint32_t raw_old_ct, uint32_t new_ct, uint32_t max_thread_ct, uint32_t* old_uidx_to_new_idx);


// Set multiplier to 0 to only count extra alleles, 1 to also count alt1 for
// those variants (useful for HWE), 2 to count both ref and alt1 for
// multiallelic variants.
// allele_idx_offsets == nullptr ok.
uintptr_t CountExtraAlleles(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t multiallelic_variant_ct_multiplier);

uint32_t MaxAlleleCtSubset(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct);

void ComputePartitionAligned(const uintptr_t* variant_include, uint32_t orig_thread_ct, uint32_t first_variant_uidx, uint32_t cur_variant_idx, uint32_t cur_variant_ct, uint32_t alignment, uint32_t* variant_uidx_starts, uint32_t* vidx_starts);

BoolErr ParseNextRange(const char* const* argvk, uint32_t param_ct, char range_delim, uint32_t* cur_param_idx_ptr, const char** cur_arg_pptr, const char** range_start_ptr, uint32_t* rs_len_ptr, const char** range_end_ptr, uint32_t* re_len_ptr);

PglErr ParseNameRanges(const char* const* argvk, const char* errstr_append, uint32_t param_ct, uint32_t require_posint, char range_delim, RangeList* range_list_ptr);


// Analytically finds all real roots of x^3 + ax^2 + bx + c, saving them in
// solutions[] (sorted from smallest to largest), and returning the count.
// Multiple roots are only returned/counted once.
uint32_t CubicRealRoots(double coef_a, double coef_b, double coef_c, STD_ARRAY_REF(double, 3) solutions);

HEADER_INLINE double PrevFloat64(double dxx) {
  uint64_t ullii;
  memcpy(&ullii, &dxx, 8);
  --ullii;
  memcpy(&dxx, &ullii, 8);
  return dxx;
}

HEADER_INLINE double NextFloat64(double dxx) {
  uint64_t ullii;
  memcpy(&ullii, &dxx, 8);
  ++ullii;
  memcpy(&dxx, &ullii, 8);
  return dxx;
}


// Pass in htable_dup_base_ptr == nullptr if just flagging duplicate IDs rather
// than tracking all their positions in item_ids.
// Otherwise, htable_dup_base entries are guaranteed to be filled in increasing
// order (briefly made that nondeterministic on 11 Oct 2019 and broke --rm-dup,
// not doing that again).
PglErr AllocAndPopulateIdHtableMt(const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uintptr_t fast_size_min_extra_bytes, uint32_t max_thread_ct, uint32_t** id_htable_ptr, uint32_t** htable_dup_base_ptr, uint32_t* id_htable_size_ptr, uint32_t* dup_ct_ptr);

typedef struct HelpCtrlStruct {
  NONCOPYABLE(HelpCtrlStruct);
  uint32_t iters_left;
  uint32_t param_ct;
  const char* const* argv;
  uintptr_t unmatched_ct;
  uintptr_t* all_match_arr;
  uintptr_t* prefix_match_arr;
  uintptr_t* perfect_match_arr;
  uint32_t* param_slens;
  uint32_t preprint_newline;
} HelpCtrl;

void HelpPrint(const char* cur_params, HelpCtrl* help_ctrl_ptr, uint32_t postprint_newline, const char* payload);


extern const char kErrstrReadCorrupted[];

// 'r' in front for (file) read.  We default to a generic "File appears to be
// corrupted" message when the C-library I/O function calls succeed (so errno
// is 0) but we determine that the logical contents of the file don't make
// sense.
HEADER_INLINE const char* rstrerror(int errnum) {
  if (errnum) {
    return strerror(errnum);
  }
  return kErrstrReadCorrupted;
}

// assumes logfile is open
void DispExitMsg(PglErr reterr);

BoolErr CheckExtraParam(const char* const* argvk, const char* permitted_modif, uint32_t* other_idx_ptr);

char ExtractCharParam(const char* ss);

PglErr CmdlineAllocString(const char* source, const char* flag_name, uint32_t max_slen, char** sbuf_ptr);

// AllocFname() and AllocFnamePrefix() now confirm file-existence.
PglErr AllocFname(const char* source, const char* flagname_p, char** fnbuf_ptr);

PglErr AllocFnamePrefix(const char* fname_prefix, const char* flattened_suffixes, const char* flagname_p, char** fnbuf_ptr);

PglErr AllocAndFlattenEx(const char* const* sources, const char* flagname_p, uint32_t param_ct, uint32_t max_blen, uint32_t check_file_existence, char** flattened_buf_ptr);

HEADER_INLINE PglErr AllocAndFlatten(const char* const* sources, const char* flagname_p, uint32_t param_ct, uint32_t max_blen, char** flattened_buf_ptr) {
  return AllocAndFlattenEx(sources, flagname_p, param_ct, max_blen, 0, flattened_buf_ptr);
}

HEADER_INLINE PglErr AllocAndFlattenFnames(const char* const* sources, const char* flagname_p, uint32_t param_ct, char** flattened_buf_ptr) {
  return AllocAndFlattenEx(sources, flagname_p, param_ct, kPglFnamesize, 1, flattened_buf_ptr);
}

typedef struct Plink2CmdlineMetaStruct {
  NONCOPYABLE(Plink2CmdlineMetaStruct);
  // need to be able to assign this to argv, so don't make it const char**
  char** subst_argv;

  char* script_buf;
  char* rerun_buf;
  char* flag_buf;
  uint32_t* flag_map;
} Plink2CmdlineMeta;

void PreinitPlink2CmdlineMeta(Plink2CmdlineMeta* pcmp);

// Handles --script, --rerun, --help, --version, and --silent.
// subst_argv, script_buf, and rerun_buf must be initialized to nullptr.
PglErr CmdlineParsePhase1(const char* ver_str, const char* ver_str2, const char* prog_name_str, const char* notestr_null_calc2, const char* cmdline_format_str, const char* errstr_append, uint32_t max_flag_blen, PglErr(* disp_help_fn)(const char* const*, uint32_t), int* argc_ptr, char*** argv_ptr, Plink2CmdlineMeta* pcmp, uint32_t* first_arg_idx_ptr, uint32_t* flag_ct_ptr);

// Assumes CmdlineParsePhase1() has completed, flag names have been copied to
// flag_buf/flag_map, aliases handled, and PROG_NAME_STR has been copied to
// outname (no null-terminator needed).  outname_end must be initialized to
// nullptr.
// This sorts the flag names so they're processed in a predictable order,
// handles --d/--out/--strict-extra-chr if present, initializes the log, and
// determines the number of processors the OS wants us to think the machine
// has.
PglErr CmdlineParsePhase2(const char* ver_str, const char* errstr_append, const char* const* argvk, uint32_t prog_name_str_slen, uint32_t max_flag_blen, int32_t argc, uint32_t flag_ct, Plink2CmdlineMeta* pcmp, char* outname, char** outname_end_ptr, char* range_delim_ptr, uint32_t* strict_extra_chr_ptr, int32_t* known_procs_ptr, uint32_t* max_thread_ct_ptr);

HEADER_INLINE void InvalidArg(const char* cur_arg) {
  logpreprintfww("Error: Unrecognized flag ('%s').\n", cur_arg);
}

PglErr CmdlineParsePhase3(uintptr_t max_default_mib, uintptr_t malloc_size_mib, uint32_t memory_require, Plink2CmdlineMeta* pcmp, unsigned char** bigstack_ua_ptr);

void CleanupPlink2CmdlineMeta(Plink2CmdlineMeta* pcmp);

ENUM_U31_DEF_START()
  kCmpExprTypeNull,

  kCmpExprTypeExists,

  kCmpExprTypeLe,
  kCmpExprTypeLeq,
  kCmpExprTypeNoteq,
  kCmpExprTypeEq,
  kCmpExprTypeGe,
  kCmpExprTypeGeq,

  kCmpExprTypeStrNoteq,
  kCmpExprTypeStrEq,
  // possible todo: kCmpExprTypeStrNotRegex, kCmpExprTypeStrRegex

  kCmpExprTypeAnd,
  kCmpExprTypeOr,

  kCmpExprTypeNot
ENUM_U31_DEF_END(CmpExprType);

HEADER_INLINE uint32_t CmpExprIsK(CmpExprType etype) {
  return (etype == kCmpExprTypeExists);
}

HEADER_INLINE uint32_t CmpExprIsKN(CmpExprType etype) {
  return (etype >= kCmpExprTypeLe) && (etype <= kCmpExprTypeGeq);
}

HEADER_INLINE uint32_t CmpExprIsKS(CmpExprType etype) {
  return (etype == kCmpExprTypeStrNoteq) || (etype == kCmpExprTypeStrEq);
}

HEADER_INLINE uint32_t CmpExprIsJct(CmpExprType etype) {
  return (etype >= kCmpExprTypeAnd);
}

struct CmpExprStruct;

// <key> existence check
typedef struct CmpExprKStruct {
  char* key;
} CmpExprK;

// <key> <operator> <number>
typedef struct CmpExprKNStruct {
  double value;
  char* key;
} CmpExprKN;

// <key> <operator> <usually-non-numeric string>
// can force numeric-string (in)equality comparison with quotes: key/str_value
// parser strips one pair of matched single-quotes if present
typedef struct CmpExprKSStruct {
  char* str_value;
  char* key;
} CmpExprKS;

// and, or, not
typedef struct CmpExprJctStruct {
  struct CmpExprStruct* children[2];
} CmpExprJct;

typedef union {
  CmpExprK k;
  CmpExprKN kn;
  CmpExprKS ks;
  CmpExprJct jct;
} CmpExprArgs;

typedef struct CmpExprStruct {
  CmpExprType etype;
  CmpExprArgs args;
} CmpExpr;

void InitCmpExpr(CmpExpr* cmp_expr_ptr);

void CleanupCmpExpr(CmpExpr* cmp_expr_ptr);

PglErr ValidateAndAllocCmpExpr(const char* const* sources, const char* flag_name, uint32_t param_ct, CmpExpr* cmp_expr_ptr);

PglErr SearchHeaderLine(const char* header_line_iter, const char* const* search_multistrs, const char* flag_nodash, uint32_t search_col_ct, uint32_t* found_col_ct_ptr, uint32_t* found_type_bitset_ptr, uint32_t* col_skips, uint32_t* col_types);

// col_descriptor is usually a pointer to argv[...][5] (first five characters
// are "cols=").  supported_ids is a multistr.
PglErr ParseColDescriptor(const char* col_descriptor_iter, const char* supported_ids, const char* cur_flag_name, uint32_t first_col_shifted, uint32_t default_cols_mask, uint32_t prohibit_empty, void* result_ptr);


// this is technically application-dependent, but let's keep this simple for
// now
// todo: recalibrate these numbers before each beta release
#ifndef __LP64__
  // 2047 seems to consistently fail on both macOS and Windows
#  ifdef _WIN32
CONSTI32(kMalloc32bitMibMax, 1728);
#  else
#    ifdef __APPLE__
CONSTI32(kMalloc32bitMibMax, 1888);
#    else
CONSTI32(kMalloc32bitMibMax, 2047);
#    endif
#  endif
#endif



#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_CMDLINE_H__
