#ifndef __PLINK2_CMDLINE_H__
#define __PLINK2_CMDLINE_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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

#include "plink2_base.h"

#include <math.h>
#include <stddef.h>

#ifndef _WIN32
  #include <sys/stat.h>
#endif

#ifdef _WIN32
  #include <process.h>
#else
  #include <pthread.h>
#endif

#ifdef __cplusplus
  #include <algorithm>
  #ifdef _WIN32
    // Windows C++11 <algorithm> resets these values :(
    #undef PRIu64
    #undef PRId64
    #define PRIu64 "I64u"
    #define PRId64 "I64d"
    #undef PRIuPTR
    #undef PRIdPTR
    #ifdef __LP64__
      #define PRIuPTR PRIu64
      #define PRIdPTR PRId64
    #else
      #if __cplusplus < 201103L
        #define PRIuPTR "lu"
        #define PRIdPTR "ld"
      #else
        #define PRIuPTR "u"
        #define PRIdPTR "d"
      #endif
    #endif
  #endif
#endif

#ifdef DYNAMIC_MKL
  #define USE_MKL
#endif

#ifdef USE_MKL
  #ifdef __APPLE__
    #error "plink2 cannot currently use MKL on OS X."
  #endif
  #ifdef LAPACK_ILP64
    #define MKL_ILP64
  #endif
  #ifdef DYNAMIC_MKL
    #include <mkl_service.h>
  #else
    #include "/opt/intel/mkl/include/mkl_service.h"
  #endif
  #define USE_MTBLAS
  #define BLAS_SET_NUM_THREADS mkl_set_num_threads
#else
  #ifdef USE_OPENBLAS
    #ifdef __cplusplus
extern "C" {
    #endif
      void openblas_set_num_threads(int num_threads);
    #ifdef __cplusplus
} // extern "C"
    #endif
    #define USE_MTBLAS
    #define BLAS_SET_NUM_THREADS openblas_set_num_threads
  #else
    #define BLAS_SET_NUM_THREADS(num)
  #endif
#endif

#ifdef _WIN32
  #define pthread_t HANDLE
  #define THREAD_FUNC_DECL unsigned __stdcall
  #define THREAD_FUNCPTR_T(func_ptr) unsigned (__stdcall *func_ptr)(void*)
  // #define THREAD_FUNCPP_T(func_pp) unsigned (__stdcall **func_pp)(void*)
  #define THREAD_RETURN return 0
  #define EOLN_STR "\r\n"
#else
  #define THREAD_FUNC_DECL void*
  #define THREAD_FUNCPTR_T(func_ptr) void* (*func_ptr)(void*)
  // #define THREAD_FUNCPP_T(func_pp) void* (**func_pp)(void*)
  #define THREAD_RETURN return nullptr
  #define EOLN_STR "\n"
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

#ifdef _GNU_SOURCE
// There was some recent (2016) discussion on the gcc mailing list on strlen()
// vs. rawmemchr(., 0), where it was claimed that rawmemchr(., 0) could be
// compiled to something slower than &(.[strlen(.)]), rather than being at
// least as good.  However, this didn't happen when I tried to benchmark this,
// so I'll stick to the function that returns the right type (except when
// rawmemchr itself has to be emulated).
HEADER_INLINE char* strnul(const char* ss) {
  return (char*)rawmemchr(ss, 0);
}
#else
  #ifdef __cplusplus
    #define rawmemchr(ss, cc) memchr((ss), (cc), (0x80000000U - plink2::kBytesPerVec))
  #else
    #define rawmemchr(ss, cc) memchr((ss), (cc), (0x80000000U - kBytesPerVec))
  #endif

HEADER_INLINE char* strnul(const char* ss) {
  return (char*)(&(ss[strlen(ss)]));
}

HEADER_INLINE char* strchrnul(const char* ss, int cc) {
  const char* strchr_result = strchr(ss, cc);
  if (strchr_result) {
    return (char*)strchr_result;
  }
  return strnul(ss);
}
#endif // !_GNU_SOURCE

#ifdef _WIN32
// if MAX_THREADS > 64, single WaitForMultipleObjects calls must be converted
// into loops
  CONSTU31(kMaxThreads, 64);
#else
// currently assumed to be less than 2^16 (otherwise some multiply overflows
// are theoretically possible, at least in the 32-bit build)
  CONSTU31(kMaxThreads, 512);
#endif

#ifdef __APPLE__
// cblas_dgemm may fail with 128k
CONSTU31(kDefaultThreadStack, 524288);
#else
// asserts didn't seem to work properly with a setting much smaller than this
CONSTU31(kDefaultThreadStack, 131072);
#endif

// generic maximum line byte length, currently also used as a default I/O
// buffer size.  .ped/.vcf/etc. lines can of course be longer.
CONSTU31(kMaxMediumLine, 131072);

CONSTU31(kLogbufSize, 2 * kMaxMediumLine);

// must be at least 2 * kMaxMediumLine + 2 to support generic token loader.
CONSTU31(kTextbufSize, 2 * kMaxMediumLine + 256);

// when g_textbuf is used as a generic I/O buffer, this is a convenient
// power-of-2 size (must be <= kTextbufSize).
CONSTU31(kTextbufMainSize, 2 * kMaxMediumLine);

// "slen" is now used to indicate string lengths excluding terminating nulls,
// while "blen" includes the terminator.

// Maximum length of chromosome, variant, FID, IID, cluster, and set IDs (not
// including terminating null).  This value supports up to 8 IDs per line
// (maximum so far is 5, for e.g. --hom).
// Assumed by plink2_pvar to be a multiple of 16.
CONSTU31(kMaxIdSlen, 16000);
CONSTU31(kMaxIdBlen, kMaxIdSlen + 1);
#define MAX_ID_SLEN_STR "16000"

// Maximum size of "dynamically" allocated line load buffer.  (This is the
// limit that applies to .vcf and similar files.)  Inconvenient to go higher
// since fgets() takes a int32_t size argument.
CONSTU31(kMaxLongLine, 0x7fffffc0);
static_assert(!(kMaxLongLine % kCacheline), "kMaxLongLine must be a multiple of kCacheline.");

// allow extensions like .model.trend.fisher.set.score.adjusted
CONSTU31(kMaxOutfnameExtBlen, 39);

#ifdef __LP64__
HEADER_CINLINE uint64_t round_up_pow2_ull(uint64_t val, uint64_t alignment) {
  return round_up_pow2(val, alignment);
}
#else
HEADER_CINLINE uint64_t round_up_pow2_ull(uint64_t val, uint64_t alignment) {
  return (val + alignment - 1) & (~(alignment - 1));
}
#endif

// file-scope string constants don't always have the g_ prefix, but multi-file
// constants are always tagged.
extern const char g_errstr_fopen[];
// extern const char g_cmdline_format_str[];

extern char g_textbuf[];

extern const char* g_one_char_strs;

// '.' missing genotype value is now taken for granted; this is in *addition*
// to it (default '0').
// (Yes, this might not belong in plink2_cmdline.)
extern const char* g_input_missing_geno_ptr;

extern const char* g_output_missing_geno_ptr; // now defaults to '.'

extern FILE* g_logfile;

// mostly-safe sprintf buffer (length kLogbufSize).  warning: do NOT put allele
// codes or arbitrary-length lists in here.
extern char g_logbuf[];

extern uint32_t g_debug_on;
extern uint32_t g_log_failed;

// for --warning-errcode
extern uint32_t g_stderr_written_to;


// warning: do NOT include allele codes (unless they're guaranteed to be SNPs)
// in log strings; they can overflow the buffer.
void logstr(const char* ss);

void logprint(const char* ss);

void logerrprint(const char* ss);

void logprintb();

void logerrprintb();

#define LOGPRINTF(...) sprintf(g_logbuf, __VA_ARGS__); logprintb();

#define LOGERRPRINTF(...) sprintf(g_logbuf, __VA_ARGS__); logerrprintb();

// input for wordwrap/LOGPRINTFWW should have no intermediate '\n's.  If
// suffix_len is 0, there should be a terminating \n.
void wordwrap(uint32_t suffix_len, char* ss);

void wordwrapb(uint32_t suffix_len);

#define LOGPREPRINTFWW(...) sprintf(g_logbuf, __VA_ARGS__); wordwrapb(0);

#define LOGPRINTFWW(...) sprintf(g_logbuf, __VA_ARGS__); wordwrapb(0); logprintb();

#define LOGERRPRINTFWW(...) sprintf(g_logbuf, __VA_ARGS__); wordwrapb(0); logerrprintb();

// 5 = length of "done." suffix, which is commonly used
#define LOGPRINTFWW5(...) sprintf(g_logbuf, __VA_ARGS__); wordwrapb(5); logprintb();

boolerr_t fopen_checked(const char* fname, const char* mode, FILE** target_ptr);

HEADER_INLINE interr_t putc_checked(int32_t ii, FILE* outfile) {
  putc_unlocked(ii, outfile);
  return ferror(outfile);
}

HEADER_INLINE interr_t fputs_checked(const char* ss, FILE* outfile) {
  fputs(ss, outfile);
  return ferror(outfile);
}

interr_t fwrite_flush2(char* buf_flush, FILE* outfile, char** write_iter_ptr);

HEADER_INLINE interr_t fwrite_ck(char* buf_flush, FILE* outfile, char** write_iter_ptr) {
  if ((*write_iter_ptr) < buf_flush) {
    return 0;
  }
  return fwrite_flush2(buf_flush, outfile, write_iter_ptr);
}

// fclose_null defined in plink2_base.h

boolerr_t fclose_flush_null(char* buf_flush, char* write_iter, FILE** outfile_ptr);

HEADER_INLINE void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

uint32_t int_slen(int32_t num);

HEADER_INLINE int32_t strequal_k(const char* s1, const char* k_s2, uint32_t s1_slen) {
  // any sane compiler should compute s2_slen at compile-time if k_s2 is a
  // constant string
  const uint32_t s2_slen = strlen(k_s2);
  return (s1_slen == s2_slen) && (!memcmp(s1, k_s2, s2_slen));
}

// Can use this when it's always safe to read first (1 + strlen(k_s2)) bytes of
// s1.
HEADER_INLINE int32_t strequal_k2(const char* s1, const char* k_s2) {
  const uint32_t s2_blen = 1 + strlen(k_s2);
  return !memcmp(s1, k_s2, s2_blen);
}

HEADER_CINLINE int32_t is_space_or_eoln(unsigned char ucc) {
  return (ucc <= 32);
}

// Assumes it's safe to read first 1 + strlen(s_const) bytes of s_read.
// Differs from strequal_k() since strings are not considered equal when
// s_read[strlen(s_const)] isn't a token-ender.
HEADER_INLINE int32_t tokequal_k(const char* s_read, const char* s_const) {
  const uint32_t s_const_slen = strlen(s_const);
  return (!memcmp(s_read, s_const, s_const_slen)) && is_space_or_eoln(s_read[s_const_slen]);
}

// s_prefix must be strictly contained.
HEADER_INLINE int32_t str_startswith(const char* s_read, const char* s_prefix_const, uint32_t s_read_slen) {
  const uint32_t s_const_slen = strlen(s_prefix_const);
  return (s_read_slen > s_const_slen) && (!memcmp(s_read, s_prefix_const, s_const_slen));
}

// permits s_read and s_prefix to be equal.
HEADER_INLINE int32_t str_startswithz(const char* s_read, const char* s_prefix_const, uint32_t s_read_slen) {
  const uint32_t s_const_slen = strlen(s_prefix_const);
  return (s_read_slen >= s_const_slen) && (!memcmp(s_read, s_prefix_const, s_const_slen));
}

// Can use this when it's always safe to read first strlen(s_prefix_const)
// bytes of s_read.
HEADER_INLINE int32_t str_startswith2(const char* s_read, const char* s_prefix_const) {
  const uint32_t s_const_slen = strlen(s_prefix_const);
  return !memcmp(s_read, s_prefix_const, s_const_slen);
}

// s_suffix must be strictly contained.
HEADER_INLINE int32_t str_endswith(const char* s_read, const char* s_suffix_const, uint32_t s_read_slen) {
  const uint32_t s_const_slen = strlen(s_suffix_const);
  return (s_read_slen > s_const_slen) && (!memcmp(&(s_read[s_read_slen - s_const_slen]), s_suffix_const, s_const_slen));
}

int32_t strcmp_casted(const void* s1, const void* s2);

int32_t strcmp_natural(const void* s1, const void* s2);

int32_t strcmp_deref(const void* s1, const void* s2);

int32_t strcmp_natural_deref(const void* s1, const void* s2);

int32_t double_cmp(const void* aa, const void* bb);

int32_t double_cmp_decr(const void* aa, const void* bb);

// requires all elements to be within 2^31 - 1 of each other
int32_t intcmp(const void* aa, const void* bb);

int32_t uint64cmp(const void* aa, const void* bb);

#ifndef __cplusplus
int32_t uint64cmp_decr(const void* aa, const void* bb);
#endif

HEADER_INLINE uint32_t get_uimax(uintptr_t len, const uint32_t* unsorted_arr) {
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

HEADER_INLINE float get_fmax(uintptr_t len, const float* unsorted_arr) {
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

HEADER_INLINE double get_dmax(uintptr_t len, const double* unsorted_arr) {
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

float destructive_get_fmedian(uintptr_t len, float* unsorted_arr);

double destructive_get_dmedian(uintptr_t len, double* unsorted_arr);

uintptr_t get_strboxsort_wentry_blen(uintptr_t max_str_blen);

#ifdef __cplusplus
typedef struct str_sort_deref_struct {
  char* strptr;
  bool operator<(const struct str_sort_deref_struct& rhs) const {
    return (strcmp(strptr, rhs.strptr) < 0);
  }
} str_sort_deref_t;

typedef struct str_nsort_deref_struct {
  char* strptr;
  bool operator<(const struct str_nsort_deref_struct& rhs) const {
    return (strcmp_natural(strptr, rhs.strptr) < 0);
  }
} str_nsort_deref_t;

HEADER_INLINE void strptr_arr_sort(uintptr_t ct, char** strptr_arr) {
  std::sort((str_sort_deref_t*)strptr_arr, &(((str_sort_deref_t*)strptr_arr)[ct]));
}

HEADER_INLINE void strptr_arr_nsort(uintptr_t ct, char** strptr_arr) {
  std::sort((str_nsort_deref_t*)strptr_arr, &(((str_nsort_deref_t*)strptr_arr)[ct]));
}

void sort_strbox_indexed2(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace);
#else
HEADER_INLINE void strptr_arr_sort(uintptr_t ct, char** strptr_arr) {
  qsort(strptr_arr, ct, sizeof(intptr_t), strcmp_deref);
}

HEADER_INLINE void strptr_arr_nsort(uintptr_t ct, char** strptr_arr) {
  qsort(strptr_arr, ct, sizeof(intptr_t), strcmp_natural_deref);
}

void sort_strbox_indexed2_fallback(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace);

HEADER_INLINE void sort_strbox_indexed2(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map, void* sort_wkspace) {
  sort_strbox_indexed2_fallback(str_ct, max_str_blen, use_nsort, strbox, id_map, sort_wkspace);
}
#endif

// This makes a temporary g_bigstack allocation.
boolerr_t sort_strbox_indexed(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map);

// Uses malloc instead of bigstack.
boolerr_t sort_strbox_indexed_malloc(uintptr_t str_ct, uintptr_t max_str_blen, char* strbox, uint32_t* id_map);

// Returns dedup'd strbox entry count.
uint32_t copy_and_dedup_sorted_strptrs_to_strbox(char** sorted_strptrs, uintptr_t str_ct, uintptr_t max_str_blen, char* strbox);


// note that this can be expected to have size 16 bytes, not 12, on 64-bit
// systems
typedef struct str_sort_indexed_deref_struct {
  const char* strptr;
  uint32_t orig_idx;
#ifdef __cplusplus
  bool operator<(const struct str_sort_indexed_deref_struct& rhs) const {
    return (strcmp(strptr, rhs.strptr) < 0);
  }
#endif
} str_sort_indexed_deref_t;

void strptr_arr_sort_main(uintptr_t str_ct, uint32_t use_nsort, str_sort_indexed_deref_t* wkspace_alias);

// This makes a temporary g_bigstack allocation.
boolerr_t strptr_arr_indexed_sort(char** unsorted_strptrs, uint32_t str_ct, uint32_t use_nsort, uint32_t* id_map);

/*
void qsort_ext2(void* main_arr, uintptr_t arr_length, uintptr_t item_length, int(* comparator_deref)(const void*, const void*), void* secondary_arr, uintptr_t secondary_item_len, unsigned char* proxy_arr, uintptr_t proxy_len);

// This makes a g_bigstack allocation, and returns -1 on alloc failure.
int32_t qsort_ext(void* main_arr, uintptr_t arr_length, uintptr_t item_length, int(* comparator_deref)(const void*, const void*), void* secondary_arr, uintptr_t secondary_item_len);
*/

uint32_t uint32arr_greater_than(const uint32_t* sorted_uint32_arr, uint32_t arr_length, uint32_t uii);

uintptr_t uint64arr_greater_than(const uint64_t* sorted_uint64_arr, uintptr_t arr_length, uint64_t ullii);

uintptr_t doublearr_greater_than(const double* sorted_dbl_arr, uintptr_t arr_length, double dxx);


uintptr_t uint64arr_geq(const uint64_t* sorted_uint64_arr, uintptr_t arr_length, uint64_t ullii);

HEADER_INLINE uint32_t is_flag(const char* param) {
  unsigned char ucc = param[1];
  return ((*param == '-') && ((ucc > '9') || ((ucc < '0') && (ucc != '.') && (ucc != '\0'))));
}

HEADER_INLINE char* is_flag_start(char* param) {
  unsigned char ucc = param[1];
  if ((*param == '-') && ((ucc > '9') || ((ucc < '0') && (ucc != '.') && (ucc != '\0')))) {
    return (ucc == '-')? (&(param[2])) : (&(param[1]));
  }
  return nullptr;
}

uint32_t param_count(char** argv, uint32_t argc, uint32_t flag_idx);

boolerr_t enforce_param_ct_range(const char* flag_name, uint32_t param_ct, uint32_t min_ct, uint32_t max_ct);

pglerr_t sort_cmdline_flags(uint32_t max_flag_blen, uint32_t flag_ct, char* flag_buf, uint32_t* flag_map);

boolerr_t cleanup_logfile(uint32_t print_end_time);

CONSTU31(kNonBigstackMin, 67108864);

CONSTU31(kBigstackMinMb, 640);
CONSTU31(kBigstackDefaultMb, 2048);

static const double kPi = 3.1415926535897932;
static const double kSqrt2 = 1.4142135623730951;
static const double kRecipE = 0.36787944117144233;
static const double kRecip2m53 = 0.00000000000000011102230246251565404236316680908203125;
static const double kRecip2m32 = 0.00000000023283064365386962890625;
static const double k2m64 = 18446744073709551616.0;

// floating point comparison-to-nonzero tolerance, currently 2^{-30}
static const double kEpsilon = 0.000000000931322574615478515625;
// less tolerant versions (2^{-35}, 2^{-44}) for some exact calculations
static const double kSmallishEpsilon = 0.00000000002910383045673370361328125;
static const double kSmallEpsilon = 0.00000000000005684341886080801486968994140625;

// 2^{-21}, must be >= sqrt(kSmallEpsilon)
static const double kBigEpsilon = 0.000000476837158203125;

// 2^{-83} bias to give exact tests maximum ability to determine tiny p-values.
// (~2^{-53} is necessary to take advantage of denormalized small numbers, then
// allow tail sum to be up to 2^30.)
static const double kExactTestBias = 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125;

// apparently these aren't always defined in limits.h
#ifndef DBL_MAX
  #define DBL_MAX 1.7976931348623157e308
#endif
#ifndef FLT_MAX
  #define FLT_MAX 3.40282347e38f
#endif

// probably time to flip arena_alloc and bigstack_alloc definitions...

// manually managed, very large double-ended stack
extern unsigned char* g_bigstack_base;
extern unsigned char* g_bigstack_end;

uintptr_t detect_mb();

uintptr_t get_default_alloc_mb();

// caller is responsible for freeing bigstack_ua
pglerr_t init_bigstack(uintptr_t malloc_size_mb, uintptr_t* malloc_mb_final_ptr, unsigned char** bigstack_ua_ptr);


HEADER_INLINE uintptr_t bigstack_left() {
  return (((uintptr_t)g_bigstack_end) - ((uintptr_t)g_bigstack_base));
}

HEADER_INLINE unsigned char* bigstack_alloc_raw(uintptr_t size) {
  // Assumes caller has already forced size to a multiple of
  // kCacheline, and verified that enough space is available.
  assert(!(size % kCacheline));
  unsigned char* alloc_ptr = g_bigstack_base;
  g_bigstack_base += size;
  return alloc_ptr;
}

HEADER_INLINE unsigned char* bigstack_alloc_raw_rd(uintptr_t size) {
  // Same as bigstack_alloc_raw(), except for rounding up size.
  unsigned char* alloc_ptr = g_bigstack_base;
  g_bigstack_base += round_up_pow2(size, kCacheline);
  return alloc_ptr;
}

// Basic 64-byte-aligned allocation at bottom of stack.
HEADER_INLINE unsigned char* bigstack_alloc(uintptr_t size) {
  size = round_up_pow2(size, kCacheline);
  if (bigstack_left() < size) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return bigstack_alloc_raw(size);
}


// Typesafe, return-0-iff-success interfaces.  (See also bigstack_calloc_...
// further below.)
HEADER_INLINE boolerr_t bigstack_alloc_c(uintptr_t ct, char** c_arr_ptr) {
  *c_arr_ptr = (char*)bigstack_alloc(ct);
  return !(*c_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_d(uintptr_t ct, double** d_arr_ptr) {
  *d_arr_ptr = (double*)bigstack_alloc(ct * sizeof(double));
  return !(*d_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_f(uintptr_t ct, float** f_arr_ptr) {
  *f_arr_ptr = (float*)bigstack_alloc(ct * sizeof(float));
  return !(*f_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_si(uintptr_t ct, int16_t** si_arr_ptr) {
  *si_arr_ptr = (int16_t*)bigstack_alloc(ct * sizeof(int16_t));
  return !(*si_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_i(uintptr_t ct, int32_t** i_arr_ptr) {
  *i_arr_ptr = (int32_t*)bigstack_alloc(ct * sizeof(int32_t));
  return !(*i_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = bigstack_alloc(ct);
  return !(*uc_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_usi(uintptr_t ct, uint16_t** usi_arr_ptr) {
  *usi_arr_ptr = (uint16_t*)bigstack_alloc(ct * sizeof(int16_t));
  return !(*usi_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_ui(uintptr_t ct, uint32_t** ui_arr_ptr) {
  *ui_arr_ptr = (uint32_t*)bigstack_alloc(ct * sizeof(int32_t));
  return !(*ui_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_ul(uintptr_t ct, uintptr_t** ul_arr_ptr) {
  *ul_arr_ptr = (uintptr_t*)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*ul_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_ll(uintptr_t ct, int64_t** ll_arr_ptr) {
  *ll_arr_ptr = (int64_t*)bigstack_alloc(ct * sizeof(int64_t));
  return !(*ll_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_ull(uintptr_t ct, uint64_t** ull_arr_ptr) {
  *ull_arr_ptr = (uint64_t*)bigstack_alloc(ct * sizeof(int64_t));
  return !(*ull_arr_ptr);
}

// some versions of gcc give aliasing warnings if we use bigstack_alloc_ul()
// for everything
// if sizeof(intptr_t) != sizeof(uintptr_t*), we're doomed anyway, so I won't
// bother with that static assert...
HEADER_INLINE boolerr_t bigstack_alloc_ulp(uintptr_t ct, uintptr_t*** ulp_arr_ptr) {
  *ulp_arr_ptr = (uintptr_t**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*ulp_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_cp(uintptr_t ct, char*** cp_arr_ptr) {
  *cp_arr_ptr = (char**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*cp_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_sip(uintptr_t ct, int16_t*** sip_arr_ptr) {
  *sip_arr_ptr = (int16_t**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*sip_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_ip(uintptr_t ct, int32_t*** ip_arr_ptr) {
  *ip_arr_ptr = (int32_t**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*ip_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_ucp(uintptr_t ct, unsigned char*** ucp_arr_ptr) {
  *ucp_arr_ptr = (unsigned char**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*ucp_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_usip(uintptr_t ct, uint16_t*** usip_arr_ptr) {
  *usip_arr_ptr = (uint16_t**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*usip_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_uip(uintptr_t ct, uint32_t*** uip_arr_ptr) {
  *uip_arr_ptr = (uint32_t**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*uip_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_dp(uintptr_t ct, double*** dp_arr_ptr) {
  *dp_arr_ptr = (double**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*dp_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_vp(uintptr_t ct, vul_t*** vp_arr_ptr) {
  *vp_arr_ptr = (vul_t**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*vp_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_thread(uintptr_t ct, pthread_t** thread_arr_ptr) {
  *thread_arr_ptr = (pthread_t*)bigstack_alloc(ct * sizeof(pthread_t));
  return !(*thread_arr_ptr);
}

HEADER_INLINE void bigstack_reset(void* new_base) {
  g_bigstack_base = (unsigned char*)new_base;
}

HEADER_INLINE void bigstack_end_reset(void* new_end) {
  g_bigstack_end = (unsigned char*)new_end;
}

HEADER_INLINE void bigstack_double_reset(void* new_base, void* new_end) {
  bigstack_reset(new_base);
  bigstack_end_reset(new_end);
}

// assumes we've already been writing to ulptr and have previously performed
// bounds-checking.
HEADER_INLINE void bigstack_finalize_ul(__maybe_unused const uintptr_t* ulptr, uintptr_t ct) {
  assert(ulptr == (const uintptr_t*)g_bigstack_base);
  g_bigstack_base += round_up_pow2(ct * sizeof(intptr_t), kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}

HEADER_INLINE void bigstack_finalize_ui(__maybe_unused const uint32_t* uiptr, uintptr_t ct) {
  assert(uiptr == (const uint32_t*)g_bigstack_base);
  g_bigstack_base += round_up_pow2(ct * sizeof(int32_t), kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}

HEADER_INLINE void bigstack_finalize_ull(__maybe_unused const uint64_t* ullptr, uintptr_t ct) {
  assert(ullptr == (const uint64_t*)g_bigstack_base);
  g_bigstack_base += round_up_pow2(ct * sizeof(int64_t), kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}

HEADER_INLINE void bigstack_finalize_c(__maybe_unused const char* cptr, uintptr_t ct) {
  assert(cptr == (const char*)g_bigstack_base);
  g_bigstack_base += round_up_pow2(ct, kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}

HEADER_INLINE void bigstack_finalize_cp(__maybe_unused char** cpptr, uintptr_t ct) {
  assert(cpptr == (char**)g_bigstack_base);
  g_bigstack_base += round_up_pow2(ct * sizeof(intptr_t), kCacheline);
  assert(g_bigstack_base <= g_bigstack_end);
}


HEADER_INLINE void bigstack_base_set(const void* unaligned_base) {
  g_bigstack_base = (unsigned char*)round_up_pow2((uintptr_t)unaligned_base, kCacheline);
}

HEADER_INLINE void bigstack_shrink_top(const void* rebase, uintptr_t new_size) {
  // could assert that this doesn't go in the wrong direction?
  g_bigstack_base = (unsigned char*)round_up_pow2(((uintptr_t)rebase) + new_size, kCacheline);
}

// simpler to have these allocations automatically AVX2-aligned when the time
// comes
CONSTU31(kEndAllocAlign, MAXV(kBytesPerVec, 16));

HEADER_INLINE void bigstack_end_set(const void* unaligned_end) {
  g_bigstack_end = (unsigned char*)round_down_pow2((uintptr_t)unaligned_end, kEndAllocAlign);
}

// assumes size is divisible by kEndAllocAlign
// assumes enough space is available
HEADER_INLINE unsigned char* bigstack_end_alloc_raw(uintptr_t size) {
  assert(!(size % kEndAllocAlign));
  g_bigstack_end -= size;
  return g_bigstack_end;
}

HEADER_INLINE unsigned char* bigstack_end_alloc_raw_rd(uintptr_t size) {
  g_bigstack_end -= round_up_pow2(size, kEndAllocAlign);
  return g_bigstack_end;
}

HEADER_INLINE unsigned char* bigstack_end_alloc_presized(uintptr_t size) {
  uintptr_t cur_bigstack_left = bigstack_left();
  if (size > cur_bigstack_left) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return bigstack_end_alloc_raw(size);
}

HEADER_INLINE unsigned char* bigstack_end_alloc(uintptr_t size) {
  size = round_up_pow2(size, kEndAllocAlign);
  return bigstack_end_alloc_presized(size);
}

HEADER_INLINE unsigned char* bigstack_end_aligned_alloc(uintptr_t size) {
  return bigstack_end_alloc(size);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_c(uintptr_t ct, char** c_arr_ptr) {
  *c_arr_ptr = (char*)bigstack_end_alloc(ct);
  return !(*c_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_d(uintptr_t ct, double** d_arr_ptr) {
  *d_arr_ptr = (double*)bigstack_end_alloc(ct * sizeof(double));
  return !(*d_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_f(uintptr_t ct, float** f_arr_ptr) {
  *f_arr_ptr = (float*)bigstack_end_alloc(ct * sizeof(float));
  return !(*f_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_i(uintptr_t ct, int32_t** i_arr_ptr) {
  *i_arr_ptr = (int32_t*)bigstack_end_alloc(ct * sizeof(int32_t));
  return !(*i_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = bigstack_end_alloc(ct);
  return !(*uc_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_ui(uintptr_t ct, uint32_t** ui_arr_ptr) {
  *ui_arr_ptr = (uint32_t*)bigstack_end_alloc(ct * sizeof(int32_t));
  return !(*ui_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_ul(uintptr_t ct, uintptr_t** ul_arr_ptr) {
  *ul_arr_ptr = (uintptr_t*)bigstack_end_alloc(ct * sizeof(intptr_t));
  return !(*ul_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_ll(uintptr_t ct, int64_t** ll_arr_ptr) {
  *ll_arr_ptr = (int64_t*)bigstack_end_alloc(ct * sizeof(int64_t));
  return !(*ll_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_ull(uintptr_t ct, uint64_t** ull_arr_ptr) {
  *ull_arr_ptr = (uint64_t*)bigstack_end_alloc(ct * sizeof(int64_t));
  return !(*ull_arr_ptr);
}

typedef struct ll_str_struct {
  struct ll_str_struct* next;
  char ss[];
} ll_str_t;

HEADER_INLINE boolerr_t bigstack_end_alloc_llstr(uintptr_t str_bytes, ll_str_t** llstr_arr_ptr) {
  *llstr_arr_ptr = (ll_str_t*)bigstack_end_alloc(str_bytes + sizeof(ll_str_t));
  return !(*llstr_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_cp(uintptr_t ct, char*** cp_arr_ptr) {
  *cp_arr_ptr = (char**)bigstack_end_alloc(ct * sizeof(intptr_t));
  return !(*cp_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_ucp(uintptr_t ct, unsigned char*** ucp_arr_ptr) {
  *ucp_arr_ptr = (unsigned char**)bigstack_end_alloc(ct * sizeof(intptr_t));
  return !(*ucp_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_thread(uintptr_t ct, pthread_t** thread_arr_ptr) {
  *thread_arr_ptr = (pthread_t*)bigstack_end_alloc(ct * sizeof(pthread_t));
  return !(*thread_arr_ptr);
}


// and here's the interface for a non-global arena (necessary for some
// multithreaded code).
HEADER_INLINE unsigned char* arena_alloc_raw(uintptr_t size, unsigned char** arena_bottom_ptr) {
  assert(!(size % kCacheline));
  unsigned char* alloc_ptr = *arena_bottom_ptr;
  *arena_bottom_ptr = &(alloc_ptr[size]);
  return alloc_ptr;
}

HEADER_INLINE unsigned char* arena_alloc_raw_rd(uintptr_t size, unsigned char** arena_bottom_ptr) {
  unsigned char* alloc_ptr = *arena_bottom_ptr;
  *arena_bottom_ptr = &(alloc_ptr[round_up_pow2(size, kCacheline)]);
  return alloc_ptr;
}

HEADER_INLINE unsigned char* arena_alloc(unsigned char* arena_top, uintptr_t size, unsigned char** arena_bottom_ptr) {
  size = round_up_pow2(size, kCacheline);
  if (((uintptr_t)(arena_top - (*arena_bottom_ptr))) < size) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return arena_alloc_raw(size, arena_bottom_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_c(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, char** c_arr_ptr) {
  *c_arr_ptr = (char*)arena_alloc(arena_top, ct, arena_bottom_ptr);
  return !(*c_arr_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_d(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, double** d_arr_ptr) {
  *d_arr_ptr = (double*)arena_alloc(arena_top, ct * sizeof(double), arena_bottom_ptr);
  return !(*d_arr_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_f(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, float** f_arr_ptr) {
  *f_arr_ptr = (float*)arena_alloc(arena_top, ct * sizeof(float), arena_bottom_ptr);
  return !(*f_arr_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_i(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, int32_t** i_arr_ptr) {
  *i_arr_ptr = (int32_t*)arena_alloc(arena_top, ct * sizeof(int32_t), arena_bottom_ptr);
  return !(*i_arr_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_uc(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = arena_alloc(arena_top, ct, arena_bottom_ptr);
  return !(*uc_arr_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_ui(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, uint32_t** ui_arr_ptr) {
  *ui_arr_ptr = (uint32_t*)arena_alloc(arena_top, ct * sizeof(int32_t), arena_bottom_ptr);
  return !(*ui_arr_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_ul(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, uintptr_t** ul_arr_ptr) {
  *ul_arr_ptr = (uintptr_t*)arena_alloc(arena_top, ct * sizeof(intptr_t), arena_bottom_ptr);
  return !(*ul_arr_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_ll(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, int64_t** ll_arr_ptr) {
  *ll_arr_ptr = (int64_t*)arena_alloc(arena_top, ct * sizeof(int64_t), arena_bottom_ptr);
  return !(*ll_arr_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_ull(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, uint64_t** ull_arr_ptr) {
  *ull_arr_ptr = (uint64_t*)arena_alloc(arena_top, ct * sizeof(int64_t), arena_bottom_ptr);
  return !(*ull_arr_ptr);
}

HEADER_INLINE boolerr_t arena_alloc_cp(unsigned char* arena_top, uintptr_t ct, unsigned char** arena_bottom_ptr, char*** cp_arr_ptr) {
  *cp_arr_ptr = (char**)arena_alloc(arena_top, ct * sizeof(intptr_t), arena_bottom_ptr);
  return !(*cp_arr_ptr);
}

HEADER_INLINE void arena_base_set(const void* unaligned_base, unsigned char** arena_bottom_ptr) {
  *arena_bottom_ptr = (unsigned char*)round_up_pow2((uintptr_t)unaligned_base, kCacheline);
}

HEADER_INLINE unsigned char* arena_end_alloc_raw(uintptr_t size, unsigned char** arena_top_ptr) {
  assert(!(size % kEndAllocAlign));
  unsigned char* alloc_ptr = *arena_top_ptr;
  alloc_ptr -= size;
  *arena_top_ptr = alloc_ptr;
  return alloc_ptr;
}

HEADER_INLINE unsigned char* arena_end_alloc(unsigned char* arena_bottom, uintptr_t size, unsigned char** arena_top_ptr) {
  size = round_up_pow2(size, kEndAllocAlign);
  if (((uintptr_t)((*arena_top_ptr) - arena_bottom)) < size) {
    g_failed_alloc_attempt_size = size;
    return nullptr;
  }
  return arena_end_alloc_raw(size, arena_top_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_c(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, char** c_arr_ptr) {
  *c_arr_ptr = (char*)arena_end_alloc(arena_bottom, ct, arena_top_ptr);
  return !(*c_arr_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_d(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, double** d_arr_ptr) {
  *d_arr_ptr = (double*)arena_end_alloc(arena_bottom, ct * sizeof(double), arena_top_ptr);
  return !(*d_arr_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_f(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, float** f_arr_ptr) {
  *f_arr_ptr = (float*)arena_end_alloc(arena_bottom, ct * sizeof(float), arena_top_ptr);
  return !(*f_arr_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_i(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, int32_t** i_arr_ptr) {
  *i_arr_ptr = (int32_t*)arena_end_alloc(arena_bottom, ct * sizeof(int32_t), arena_top_ptr);
  return !(*i_arr_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_uc(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = arena_end_alloc(arena_bottom, ct, arena_top_ptr);
  return !(*uc_arr_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_ui(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, uint32_t** ui_arr_ptr) {
  *ui_arr_ptr = (uint32_t*)arena_end_alloc(arena_bottom, ct * sizeof(int32_t), arena_top_ptr);
  return !(*ui_arr_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_ul(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, uintptr_t** ul_arr_ptr) {
  *ul_arr_ptr = (uintptr_t*)arena_end_alloc(arena_bottom, ct * sizeof(intptr_t), arena_top_ptr);
  return !(*ul_arr_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_ll(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, int64_t** ll_arr_ptr) {
  *ll_arr_ptr = (int64_t*)arena_end_alloc(arena_bottom, ct * sizeof(int64_t), arena_top_ptr);
  return !(*ll_arr_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_ull(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, uint64_t** ull_arr_ptr) {
  *ull_arr_ptr = (uint64_t*)arena_end_alloc(arena_bottom, ct * sizeof(int64_t), arena_top_ptr);
  return !(*ull_arr_ptr);
}

HEADER_INLINE boolerr_t arena_end_alloc_cp(unsigned char* arena_bottom, uintptr_t ct, unsigned char** arena_top_ptr, char*** cp_arr_ptr) {
  *cp_arr_ptr = (char**)arena_end_alloc(arena_bottom, ct * sizeof(intptr_t), arena_top_ptr);
  return !(*cp_arr_ptr);
}


boolerr_t push_llstr(const char* ss, ll_str_t** ll_stack_ptr);

typedef struct l32str_struct {
  uint32_t len;
  char ss[];
} l32str_t;


HEADER_CINLINE int32_t is_letter(unsigned char ucc) {
  return (((ucc & 192) == 64) && (((ucc - 1) & 31) < 26));
}

// if we need the digit value, better to use (unsigned char)cc - '0'...
HEADER_CINLINE int32_t is_digit(unsigned char ucc) {
  return (ucc <= '9') && (ucc >= '0');
}

HEADER_CINLINE int32_t is_not_digit(unsigned char ucc) {
  return (ucc > '9') || (ucc < '0');
}

HEADER_CINLINE int32_t is_not_nzdigit(unsigned char ucc) {
  return (ucc > '9') || (ucc <= '0');
}

// may as well treat all chars < 32, except tab, as eoln...
// kns = "known non-space" (where tab counts as a space)
/*
HEADER_INLINE int32_t is_eoln_kns(unsigned char ucc) {
  return (ucc < 32);
}
*/

HEADER_CINLINE int32_t is_eoln_kns(unsigned char ucc) {
  // could assert ucc is not a space/tab?
  return (ucc <= 32);
}

HEADER_CINLINE int32_t is_eoln_or_comment_kns(unsigned char ucc) {
  return (ucc < 32) || (ucc == '#');
}

HEADER_CINLINE int32_t no_more_tokens_kns(const char* sptr) {
  return ((!sptr) || is_eoln_kns(*sptr));
}

HEADER_INLINE char* skip_initial_spaces(char* sptr) {
  while ((*sptr == ' ') || (*sptr == '\t')) {
    ++sptr;
  }
  return sptr;
}

// assumes we are currently in a token -- UNSAFE OTHERWISE
HEADER_INLINE char* token_endnn(char* sptr) {
  // assert(((unsigned char)(*sptr)) > 32);
  while (!is_space_or_eoln(*(++sptr)));
  return sptr;
}

HEADER_INLINE char* next_prespace(char* sptr) {
  while (((unsigned char)(*(++sptr))) >= 32);
  return sptr;
}


// length-zero tokens and non-leading spaces are permitted in the
// comma-delimiter case
HEADER_INLINE char* comma_or_space_token_end(char* token_iter, uint32_t comma_delim) {
  if (comma_delim) {
    unsigned char ucc = (unsigned char)(*token_iter);
    while ((ucc >= ' ') && (ucc != ',')) {
      ucc = *(++token_iter);
    }
    return token_iter;
  }
  return token_endnn(token_iter);
}

HEADER_INLINE char* comma_or_space_next_token(char* token_end_iter, uint32_t comma_delim) {
  // assumes token_end_iter is non-null, returns nullptr if there are no more
  // tokens
  // assert(token_end_iter);
  if (comma_delim) {
    if ((*token_end_iter) != ',') {
      return nullptr;
    }
    return skip_initial_spaces(&(token_end_iter[1]));
  }
  char* ss = skip_initial_spaces(token_end_iter);
  return is_eoln_kns(*ss)? nullptr : ss;
}


// Returns whether uppercased ss matches nonempty fixed_str.  Assumes fixed_str
// contains nothing but letters and a null terminator.
// uint32_t match_upper(const char* ss, const char* fixed_str);

uint32_t match_upper_counted(const char* ss, const char* fixed_str, uint32_t ct);

HEADER_INLINE uint32_t match_upper_k(const char* ss, const char* fixed_str, uint32_t ss_slen) {
  const uint32_t fixed_slen = strlen(fixed_str);
  if (ss_slen != fixed_slen) {
    return 0;
  }
  return match_upper_counted(ss, fixed_str, fixed_slen);
}

HEADER_INLINE uint32_t match_upper_k2(const char* ss, const char* fixed_str) {
  return match_upper_counted(ss, fixed_str, strlen(fixed_str));
}

/*
void str_toupper(char* ss);

void buf_toupper(uint32_t slen, char* ss);

void strcpy_toupper(char* target, const char* source);
*/

uint32_t is_alphanumeric(const char* ss);

// scan_posint_capped(), scan_uint_capped(), scan_int_abs_bounded(),
// scan_int32(), scan_posint_defcap(), scan_uint_defcap(),
// scan_int_abs_defcap(), scan_uint_icap() in plink2_base

boolerr_t scan_posintptr(const char* ss, uintptr_t* valp);

#ifdef __LP64__
boolerr_t scanadv_posint_capped(uint64_t cap, char** ss_ptr, uint32_t* valp);

boolerr_t scanadv_uint_capped(uint64_t cap, char** ss_ptr, uint32_t* valp);
#else
boolerr_t scanadv_posint_capped32(uint32_t cap_div_10, uint32_t cap_mod_10, char** ss_ptr, uint32_t* valp);

boolerr_t scanadv_uint_capped32(uint32_t cap_div_10, uint32_t cap_mod_10, char** ss_ptr, uint32_t* valp);

HEADER_INLINE boolerr_t scanadv_posint_capped(uint32_t cap, char** ss_ptr, uint32_t* valp) {
 return scanadv_posint_capped32(cap / 10, cap % 10, ss_ptr, valp);
}

HEADER_INLINE boolerr_t scanadv_uint_capped(uint32_t cap, char** ss_ptr, uint32_t* valp) {
 return scanadv_uint_capped32(cap / 10, cap % 10, ss_ptr, valp);
}
#endif

HEADER_INLINE boolerr_t scanadv_uint_defcap(char** ss_ptr, uint32_t* valp) {
  return scanadv_uint_capped(0x7ffffffe, ss_ptr, valp);
}

// this has different semantics from scanadv_posint_capped, etc. since integer
// readers don't take much code (so it's fine to have a bunch of similar
// functions, optimized for slightly different use cases), but we only want one
// floating point reader
char* scanadv_double(char* ss, double* valp);

HEADER_INLINE boolerr_t scan_float(const char* ss, float* valp) {
  double dxx;
  // const_cast
  if (!scanadv_double((char*)((uintptr_t)ss), &dxx)) {
    return 1;
  }
  if (fabs(dxx) > 3.4028235677973362e38) {
    return 1;
  }
  *valp = (float)dxx;
  return 0;
}

// memcpya(), memseta() defined in plink2_base.h

HEADER_INLINE char* memcpyax(void* __restrict target, const void* __restrict source, uint32_t ct, char extra_char) {
  memcpy(target, source, ct);
  ((char*)target)[ct] = extra_char;
  return &(((char*)target)[ct + 1]);
}

HEADER_INLINE void memcpyx(void* __restrict target, const void* __restrict source, uint32_t ct, char extra_char) {
  memcpy(target, source, ct);
  ((char*)target)[ct] = extra_char;
}

HEADER_INLINE void memcpyl3(void* __restrict target, const void* __restrict source) {
  // when it's safe to clobber the fourth character, this is faster
  *((uint32_t*)target) = *((const uint32_t*)source);
}

HEADER_INLINE char* memcpyl3a(void* __restrict target, const void* __restrict source) {
  memcpyl3(target, source);
  return &(((char*)target)[3]);
}

HEADER_INLINE char* strcpya(char* __restrict target, const void* __restrict source) {
  uintptr_t slen = strlen((const char*)source);
  memcpy(target, source, slen);
  return &(target[slen]);
}

HEADER_INLINE char* strcpyax(char* __restrict target, const void* __restrict source, char extra_char) {
  uintptr_t slen = strlen((const char*)source);
  memcpy(target, source, slen);
  target[slen] = extra_char;
  return &(target[slen + 1]);
}

// MinGW support for stpcpy is a mess, so I'll use a different name
// ("strcpya0") which doesn't depend on MinGW knowing what it's doing
#if defined(_GNU_SOURCE) || defined(__APPLE__) || (_POSIX_C_SOURCE >= 200809L)
HEADER_INLINE char* strcpya0(char* __restrict target, const char* __restrict source) {
  return stpcpy(target, source);
}
#else
HEADER_INLINE char* strcpya0(char* __restrict target, const char* __restrict source) {
  uintptr_t slen = strlen(source);
  memcpy(target, source, slen + 1);
  return &(target[slen]);
}
#endif

HEADER_INLINE void append_binary_eoln(char** target_ptr) {
#ifdef _WIN32
  (*target_ptr)[0] = '\r';
  (*target_ptr)[1] = '\n';
  *target_ptr += 2;
#else
  **target_ptr = '\n';
  *target_ptr += 1;
#endif
}

HEADER_INLINE void decr_append_binary_eoln(char** target_ptr) {
#ifdef _WIN32
  (*target_ptr)[-1] = '\r';
  (*target_ptr)[0] = '\n';
  *target_ptr += 1;
#else
  (*target_ptr)[-1] = '\n';
#endif
}

void get_top_two_ui(const uint32_t* __restrict uint_arr, uintptr_t uia_size, uintptr_t* __restrict top_idx_ptr, uintptr_t* __restrict second_idx_ptr);

// safer than token_endnn(), since it handles length zero
// "se" = space/eoln treated as terminators
HEADER_INLINE uintptr_t strlen_se(const char* ss) {
  const char* ss2 = ss;
  while (!is_space_or_eoln(*ss2)) {
    ss2++;
  }
  return (uintptr_t)(ss2 - ss);
}

// ok if sptr is at end of current token
HEADER_INLINE char* next_token(char* sptr) {
  if (!sptr) {
    return nullptr;
  }
  unsigned char ucc = *sptr;
  while (ucc > 32) {
    ucc = *(++sptr);
  }
  while ((ucc == ' ') || (ucc == '\t')) {
    ucc = *(++sptr);
  }
  return (ucc > 32)? sptr : nullptr;
}

HEADER_INLINE char* next_token_mult(char* sptr, uint32_t ct) {
  // assert(ct);
  if (!sptr) {
    return nullptr;
  }
  unsigned char ucc = *sptr;
  do {
    while (ucc > 32) {
      ucc = *(++sptr);
    }
    while ((ucc == ' ') || (ucc == '\t')) {
      ucc = *(++sptr);
    }
    if (ucc <= 32) {
      return nullptr;
    }
  } while (--ct);
  return sptr;
}

HEADER_INLINE char* next_token_multz(char* sptr, uint32_t ct) {
  // tried replacing this with ternary operator, but that actually seemed to
  // slow things down a bit under gcc 4.2.1 (tail call optimization issue?).
  // todo: recheck this under newer gcc/clang.
  if (ct) {
    return next_token_mult(sptr, ct);
  }
  return sptr;
}

char* comma_or_space_next_token_mult(char* sptr, uint32_t ct, uint32_t comma_delim);

uint32_t count_tokens(const char* bufptr);

// uint32_t comma_or_space_count_tokens(const char* bufptr, uint32_t comma_delim);

HEADER_INLINE char* fw_strcpyn(const char* source, uint32_t min_width, uint32_t source_len, char* dest) {
  // right-justified strcpy with known source length
  if (source_len < min_width) {
    memcpy(memseta((unsigned char*)dest, ' ', min_width - source_len), source, source_len);
    return &(dest[min_width]);
  }
  return (char*)memcpya(dest, source, source_len);
}

HEADER_INLINE char* fw_strcpy(const char* source, uint32_t min_width, char* dest) {
  return fw_strcpyn(source, min_width, strlen(source), dest);
}

// uint32_t count_and_measure_multistr(const char* multistr, uintptr_t* max_blen_ptr);

boolerr_t count_and_measure_multistr_reverse_alloc(char* multistr, uintptr_t max_str_ct, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr, char*** strptr_arrp);

extern const uint16_t kDigitPair[];

char* uint32toa(uint32_t uii, char* start);

char* int32toa(int32_t ii, char* start);

char* uitoa_z5(uint32_t uii, char* start);

char* int64toa(int64_t llii, char* start);

char* uitoa_trunc4(uint32_t uii, char* start);

char* dtoa_g(double dxx, char* start);

// We try to avoid micromanaging floating point printing and just use %g
// everywhere, but occasionally we explicitly need more precision.
//
// dtoa_g_p8 provides generic 8-digit precision (instead of %g's 6-digit
// default), while print_dosage in plink2_common provides up to 3 places after
// the decimal point when dealing with dosages (which are internally
// represented as 32768ths).
// (may want to replace _p8 with _p10 for perfect int32 handling.)
char* dtoa_g_p8(double dxx, char* start);

HEADER_INLINE void trailing_zeroes_to_spaces(char* start) {
  --start;
  while (*start == '0') {
    *start-- = ' ';
  }
  if (*start == '.') {
    *start = ' ';
  }
}

HEADER_INLINE char* clip_trailing_zeroes(char* start) {
  char cc;
  do {
    cc = *(--start);
  } while (cc == '0');
  return &(start[(cc != '.')]);
}

// "prob" means that the number is guaranteed to be in [0, 1].
// no leading space is printed.  trailing zeroes (/decimal point) are erased
//   iff there is equality to ~13 decimal places.
char* dtoa_f_probp6_spaced(double dxx, char* start);

char* dtoa_f_probp6_clipped(double dxx, char* start);

// char* dtoa_f_p5_clipped(double dxx, char* start);

char* ftoa_g(float fxx, char* start);

HEADER_INLINE char* uint32toa_x(uint32_t uii, char extra_char, char* start) {
  char* penult = uint32toa(uii, start);
  *penult = extra_char;
  return &(penult[1]);
}

// fill_uint_zero, fill_ulong_zero, fill_ull_zero, fill_ulong_one currently
// defined in plink2_base.h

HEADER_INLINE void fill_vvec_zero(uintptr_t entry_ct, vul_t* vvec) {
  memset(vvec, 0, entry_ct * kBytesPerVec);
}

HEADER_INLINE void fill_ull_one(uintptr_t entry_ct, uint64_t* ullarr) {
  fill_ulong_one(entry_ct, (uintptr_t*)ullarr);
}

HEADER_INLINE void fill_int_zero(uintptr_t entry_ct, int32_t* iarr) {
  memset(iarr, 0, entry_ct * sizeof(int32_t));
}

HEADER_INLINE void fill_int_one(uintptr_t entry_ct, int32_t* iarr) {
  for (uintptr_t ulii = 0; ulii < entry_ct; ulii++) {
    *iarr++ = -1;
  }
}

HEADER_INLINE void fill_uint_one(uintptr_t entry_ct, uint32_t* uiarr) {
  for (uintptr_t ulii = 0; ulii < entry_ct; ulii++) {
    *uiarr++ = ~0U;
  }
}

HEADER_INLINE void fill_float_zero(uintptr_t entry_ct, float* farr) {
  for (uintptr_t ulii = 0; ulii < entry_ct; ulii++) {
    *farr++ = 0.0;
  }
}

HEADER_INLINE void fill_double_zero(uintptr_t entry_ct, double* darr) {
  for (uintptr_t ulii = 0; ulii < entry_ct; ulii++) {
    *darr++ = 0.0;
  }
}

HEADER_INLINE void fill_ptr_zero(uintptr_t entry_ct, void* pp) {
  memset(pp, 0, entry_ct * sizeof(intptr_t));
}


HEADER_INLINE void append_float_zero(uintptr_t entry_ct, float** farr_ptr) {
  float* farr = *farr_ptr;
  for (uintptr_t ulii = 0; ulii < entry_ct; ulii++) {
    *farr++ = 0.0;
  }
  *farr_ptr = farr;
}


// void magic_num(uint32_t divisor, uint64_t* multp, uint32_t* __restrict pre_shiftp, uint32_t* __restrict post_shiftp, uint32_t* __restrict incrp);


// fill_all_bits, IS_SET, SET_BIT, CLEAR_BIT, next_set_unsafe,
// next_set_unsafe_ck, next_unset_unsafe, next_unset_unsafe_ck,
// next_set, prev_set_unsafe, are_all_words_zero defined in plink2_base.h

// use this instead of IS_SET() for signed 32-bit integers
HEADER_INLINE uint32_t is_set(const uintptr_t* bitarr, uint32_t loc) {
  return (bitarr[loc / kBitsPerWord] >> (loc % kBitsPerWord)) & 1;
}

// useful for coercing int32_t loc to unsigned
HEADER_INLINE void set_bit(uint32_t loc, uintptr_t* bitarr) {
  bitarr[loc / kBitsPerWord] |= (k1LU << (loc % kBitsPerWord));
}

HEADER_INLINE void clear_bit(uint32_t loc, uintptr_t* bitarr) {
  bitarr[loc / kBitsPerWord] &= ~(k1LU << (loc % kBitsPerWord));
}

#define FLIP_BIT(idx, arr) ((arr)[(idx) / kBitsPerWord] ^= k1LU << ((idx) % kBitsPerWord))

// "_nz" added to names to make it obvious these require positive len
void fill_bits_nz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr);
void clear_bits_nz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr);

#ifdef __LP64__
uintptr_t next_set_ul_unsafe(const uintptr_t* bitarr, uintptr_t loc);
#else
HEADER_INLINE uintptr_t next_set_ul_unsafe(const uintptr_t* bitarr, uintptr_t loc) {
  return (uintptr_t)next_set_unsafe(bitarr, loc);
}
#endif

HEADER_INLINE void next_set_ul_unsafe_ck(const uintptr_t* __restrict bitarr, uintptr_t* __restrict loc_ptr) {
  if (!IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_set_ul_unsafe(bitarr, *loc_ptr);
  }
}

uint32_t next_unset(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil);

// floor permitted to be -1, though not smaller than that.
int32_t prev_set(const uintptr_t* bitarr, uint32_t loc, int32_t floor);

HEADER_INLINE uint32_t are_all_words_identical(const uintptr_t* word_arr1, const uintptr_t* word_arr2, uintptr_t word_ct) {
  for (uintptr_t widx = 0; widx < word_ct; ++widx) {
    if (word_arr1[widx] ^ word_arr2[widx]) {
      return 0;
    }
  }
  return 1;
}


boolerr_t bigstack_calloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr);

boolerr_t bigstack_calloc_d(uintptr_t ct, double** d_arr_ptr);

boolerr_t bigstack_calloc_f(uintptr_t ct, float** f_arr_ptr);

boolerr_t bigstack_calloc_usi(uintptr_t ct, uint16_t** usi_arr_ptr);

boolerr_t bigstack_calloc_ui(uintptr_t ct, uint32_t** ui_arr_ptr);

boolerr_t bigstack_calloc_ul(uintptr_t ct, uintptr_t** ul_arr_ptr);

boolerr_t bigstack_calloc_ull(uintptr_t ct, uint64_t** ull_arr_ptr);

HEADER_INLINE boolerr_t bigstack_calloc_c(uintptr_t ct, char** c_arr_ptr) {
  return bigstack_calloc_uc(ct, (unsigned char**)c_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_calloc_si(uintptr_t ct, int16_t** si_arr_ptr) {
  return bigstack_calloc_usi(ct, (uint16_t**)si_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_calloc_i(uintptr_t ct, int32_t** i_arr_ptr) {
  return bigstack_calloc_ui(ct, (uint32_t**)i_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_calloc_ll(uintptr_t ct, int64_t** ll_arr_ptr) {
  return bigstack_calloc_ull(ct, (uint64_t**)ll_arr_ptr);
}

boolerr_t bigstack_end_calloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr);

boolerr_t bigstack_end_calloc_d(uintptr_t ct, double** d_arr_ptr);

boolerr_t bigstack_end_calloc_f(uintptr_t ct, float** f_arr_ptr);

boolerr_t bigstack_end_calloc_ui(uintptr_t ct, uint32_t** ui_arr_ptr);

boolerr_t bigstack_end_calloc_ul(uintptr_t ct, uintptr_t** ul_arr_ptr);

boolerr_t bigstack_end_calloc_ull(uintptr_t ct, uint64_t** ull_arr_ptr);

HEADER_INLINE boolerr_t bigstack_end_calloc_c(uintptr_t ct, char** c_arr_ptr) {
  return bigstack_end_calloc_uc(ct, (unsigned char**)c_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_calloc_i(uintptr_t ct, int32_t** i_arr_ptr) {
  return bigstack_end_calloc_ui(ct, (uint32_t**)i_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_calloc_ll(uintptr_t ct, int64_t** ll_arr_ptr) {
  return bigstack_end_calloc_ull(ct, (uint64_t**)ll_arr_ptr);
}


// These ensure the trailing bits are zeroed out.
void bitarr_invert(uintptr_t bit_ct, uintptr_t* bitarr);

void bitarr_invert_copy(const uintptr_t* __restrict source_bitarr, uintptr_t bit_ct, uintptr_t* __restrict target_bitarr);

// bitvec_and(), bitvec_andnot() in plink2_base.h

void bitvec_and_copy(const uintptr_t* __restrict source1_bitvec, const uintptr_t* __restrict source2_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec);

void bitvec_andnot_copy(const uintptr_t* __restrict source_bitvec, const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec);

void bitvec_or(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec);

void bitvec_andnot2(const uintptr_t* __restrict include_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

// void bitvec_ornot(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec);


// basic linear scan
// returns -1 on failure to find, -2 if duplicate
int32_t get_variant_uidx_without_htable(const char* idstr, char** variant_ids, const uintptr_t* variant_include, uint32_t variant_ct);

// copy_subset() doesn't exist since a loop of the form
//   uint32_t uidx = 0;
//   for (uint32_t idx = 0; idx < subset_size; ++idx, ++uidx) {
//     next_set_unsafe_ck(subset_mask, &uidx);
//     *target_iter++ = source_arr[uidx];
//   }
// seems to compile better?

// void copy_when_nonmissing(const uintptr_t* loadbuf, const void* source, uintptr_t elem_size, uintptr_t unfiltered_sample_ct, uintptr_t missing_ct, void* dest);


// tried to replace this with a faster hash function, but turns out it's hard
// to meaningfully beat (and multithreading parts of hash table construction
// solved most of the initialization time issue, anyway)
// eventually want this to be a constexpr?  seems painful to make that work,
// though.
uint32_t murmurhash3_32(const void* key, uint32_t len);

// see http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
HEADER_INLINE uint32_t hashceil(const char* idstr, uint32_t idlen, uint32_t htable_size) {
  return (((uint64_t)murmurhash3_32(idstr, idlen)) * htable_size) >> 32;
}

// uintptr_t geqprime(uintptr_t floor);

// assumes ceil is odd and greater than 4.  Returns the first prime <= ceil.
// uintptr_t leqprime(uintptr_t ceil);

HEADER_INLINE uint32_t get_htable_min_size(uintptr_t item_ct) {
  if (item_ct > 6) {
    // Don't actually need primes any more.
    // return geqprime(item_ct * 2 + 1);
    return item_ct * 2;
  }
  return 13;
}

// load factor ~20% seems to yield the best speed/space tradeoff on my test
// machines
HEADER_INLINE uint32_t get_htable_fast_size(uint32_t item_ct) {
  // if (item_ct < 858993456) {
  //   return geqprime(round_up_pow2(item_ct * 5, 2) + 1);
  // }
  if (item_ct < 858993459) {
    return item_ct * 5;
  }
  return 4294967291U;
}

boolerr_t htable_good_size_alloc(uint32_t item_ct, uintptr_t bytes_avail, uint32_t** htable_ptr, uint32_t* htable_size_ptr);

// useful for duplicate detection: returns 0 on no duplicates, a positive index
// of a duplicate pair if they're present
uint32_t populate_strbox_htable(const char* strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t str_htable_size, uint32_t* str_htable);

// returned index in duplicate-pair case is unfiltered
// uint32_t populate_strbox_subset_htable(const uintptr_t* __restrict subset_mask, const char* strbox, uintptr_t raw_str_ct, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t str_htable_size, uint32_t* str_htable);

// cur_id DOES need to be null-terminated
uint32_t id_htable_find(const char* cur_id, char** item_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size);

// assumes cur_id_slen < max_str_blen.
// requires cur_id to be null-terminated.
uint32_t strbox_htable_find(const char* cur_id, const char* strbox, const uint32_t* id_htable, uintptr_t max_str_blen, uint32_t cur_id_slen, uint32_t id_htable_size);

// last variant_ids entry must be at least kMaxIdBlen bytes before end of
// bigstack
uint32_t variant_id_dupflag_htable_find(const char* idbuf, char** variant_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen);

uint32_t variant_id_dup_htable_find(const char* idbuf, char** variant_ids, const uint32_t* id_htable, const uint32_t* htable_dup_base, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen, uint32_t* llidx_ptr);

char* scan_for_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_blen);

// Collapses array of sorted IDs to remove duplicates, and writes pre-collapse
// positions to id_starts (so e.g. duplication count of any sample ID can be
// determined via subtraction) if it isn't nullptr.
// Returns id_ct of collapsed array.
uint32_t collapse_duplicate_ids(uintptr_t id_ct, uintptr_t max_id_blen, char* sorted_ids, uint32_t* id_starts);

pglerr_t copy_sort_strbox_subset_noalloc(const uintptr_t* __restrict subset_mask, const char* __restrict orig_strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t allow_dups, uint32_t collapse_idxs, uint32_t use_nsort, char* __restrict sorted_strbox, uint32_t* __restrict id_map);

pglerr_t copy_sort_strbox_subset(const uintptr_t* __restrict subset_mask, const char* __restrict orig_strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t allow_dups, uint32_t collapse_idxs, uint32_t use_nsort, char** sorted_strbox_ptr, uint32_t** id_map_ptr);


// returns position of string, or -1 if not found.
int32_t bsearch_str(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx);

// requires null-terminated string
int32_t bsearch_str_natural(const char* idbuf, const char* sorted_strbox, uintptr_t max_id_blen, uintptr_t end_idx);


// returns number of elements in sorted_strbox[] less than idbuf.
uintptr_t bsearch_str_lb(const char* idbuf, const char* sorted_strbox, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx);

// this is frequently preferable to bsearch_str(), since it's way too easy to
// forget to convert the sorted-stringbox index to the final index
// sample_id_map == nullptr is permitted; in this case id will be an index into
// the sorted array
HEADER_INLINE boolerr_t sorted_idbox_find(const char* idbuf, const char* sorted_idbox, const uint32_t* id_map, uintptr_t cur_id_slen, uintptr_t max_id_blen, uintptr_t end_idx, uint32_t* id_ptr) {
  const int32_t ii = bsearch_str(idbuf, sorted_idbox, cur_id_slen, max_id_blen, end_idx);
  if (ii == -1) {
    return 1;
  }
  *id_ptr = id_map? id_map[(uint32_t)ii] : ((uint32_t)ii);
  return 0;
}


typedef struct range_list_struct {
  char* names;
  unsigned char* starts_range;
  uint32_t name_ct;
  uint32_t name_max_blen;
} range_list_t;

void init_range_list(range_list_t* range_list_ptr);

void cleanup_range_list(range_list_t* range_list_ptr);

// bitarr assumed to be initialized (but not necessarily zero-initialized)
boolerr_t numeric_range_list_to_bitarr(const range_list_t* range_list_ptr, uint32_t bitarr_size, uint32_t offset, uint32_t ignore_overflow, uintptr_t* bitarr);

pglerr_t string_range_list_to_bitarr(char* header_line, const range_list_t* range_list_ptr, const char* __restrict sorted_ids, const uint32_t* __restrict id_map, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t* bitarr, int32_t* __restrict seen_idxs);

pglerr_t string_range_list_to_bitarr_alloc(char* header_line, const range_list_t* range_list_ptr, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t** bitarr_ptr);


HEADER_CINLINE uint32_t realnum(double dd) {
  return (dd == dd) && (dd != INFINITY) && (dd != -INFINITY);
}


#ifdef __LP64__
HEADER_INLINE uintptr_t popcount_longs_nzbase(const uintptr_t* bitvec, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t prefix_ct = 0;
  #ifdef USE_AVX2
  while (start_idx & 3) {
    if (end_idx == start_idx) {
      return prefix_ct;
    }
    prefix_ct += popcount_long(bitvec[start_idx++]);
  }
  #else
  if (start_idx & 1) {
    if (end_idx == start_idx) {
      return 0;
    }
    prefix_ct = popcount_long(bitvec[start_idx++]);
  }
  #endif // USE_AVX2
  return prefix_ct + popcount_longs(&(bitvec[start_idx]), end_idx - start_idx);
}
#else
HEADER_INLINE uintptr_t popcount_longs_nzbase(const uintptr_t* bitvec, uintptr_t start_idx, uintptr_t end_idx) {
  return popcount_longs(&(bitvec[start_idx]), end_idx - start_idx);
}
#endif

uintptr_t popcount_bit_idx(const uintptr_t* bitvec, uintptr_t start_idx, uintptr_t end_idx);

uintptr_t popcount_longs_intersect(const uintptr_t* __restrict bitvec1_iter, const uintptr_t* __restrict bitvec2_iter, uintptr_t word_ct);

void popcount_longs_intersect_3val(const uintptr_t* __restrict bitvec1, const uintptr_t* __restrict bitvec2, uint32_t word_ct, uint32_t* __restrict popcount1_ptr, uint32_t* __restrict popcount2_ptr, uint32_t* __restrict popcount_intersect_ptr);

// uintptr_t count_11_longs(const uintptr_t* genovec, uintptr_t word_ct);

uint32_t are_all_bits_zero(const uintptr_t* bitarr, uintptr_t start_idx, uintptr_t end_idx);

// assumes len is positive, and relevant bits of target_bitarr are zero
void copy_bitarr_range(const uintptr_t* __restrict src_bitarr, uintptr_t src_start_bitidx, uintptr_t target_start_bitidx, uintptr_t len, uintptr_t* __restrict target_bitarr);


// vertical popcount support
// scramble_1_4_8_32() and friends are here since they apply to generic
// arrays; scramble_2_4_8_32() is more plink-specific

#ifdef __LP64__
  #ifdef USE_AVX2
HEADER_CINLINE uint32_t scramble_1_4_8_32(uint32_t orig_idx) {
  // 1->4: 0 4 8 12 16 20 24 28 32 ... 252 1 5 9 ...
  // 4->8: 0 8 16 24 32 ... 248 4 12 20 ... 1 9 17 ...
  // 8->32: 0 32 ... 224 8 40 ... 232 ... 248 4 36 ... 252 1 33 ...
  return (orig_idx & (~255)) + ((orig_idx & 3) * 64) + ((orig_idx & 4) * 8) + (orig_idx & 24) + ((orig_idx & 224) / 32);
}
  #else
HEADER_CINLINE uint32_t scramble_1_4_8_32(uint32_t orig_idx) {
  // 1->4: 0 4 8 12 16 20 24 28 32 ... 124 1 5 9 ...
  // 4->8: 0 8 16 24 32 ... 120 4 12 20 ... 1 9 17 ...
  // 8->32: 0 32 64 96 8 40 72 104 16 48 80 112 24 56 88 120 4 36 68 ... 1 33 ...
  return (orig_idx & (~127)) + ((orig_idx & 3) * 32) + ((orig_idx & 4) * 4) + ((orig_idx & 24) / 2) + ((orig_idx & 96) / 32);
}
  #endif
#else
// 1->4: 0 4 8 12 16 20 24 28 1 5 9 13 17 21 25 29 2 6 10 ...
// 4->8: 0 8 16 24 4 12 20 28 1 9 17 25 5 13 21 29 2 10 18 ...
// 8->32: 0 8 16 24 4 12 20 28 1 9 17 25 5 13 21 29 2 10 18 ...
HEADER_CINLINE uint32_t scramble_1_4_8_32(uint32_t orig_idx) {
  return (orig_idx & (~31)) + ((orig_idx & 3) * 8) + (orig_idx & 4) + ((orig_idx & 24) / 8);
}
#endif

HEADER_INLINE void unroll_incr_1_4(const uintptr_t* acc1, uint32_t acc1_vec_ct, uintptr_t* acc4) {
  const vul_t m1x4 = VCONST_UL(kMask1111);
  const vul_t* acc1v_iter = (const vul_t*)acc1;
  vul_t* acc4v_iter = (vul_t*)acc4;
  for (uint32_t vidx = 0; vidx < acc1_vec_ct; ++vidx) {
    vul_t loader = *acc1v_iter++;
    *acc4v_iter = (*acc4v_iter) + (loader & m1x4);
    ++acc4v_iter;
    loader = vul_rshift(loader, 1);
    *acc4v_iter = (*acc4v_iter) + (loader & m1x4);
    ++acc4v_iter;
    loader = vul_rshift(loader, 1);
    *acc4v_iter = (*acc4v_iter) + (loader & m1x4);
    ++acc4v_iter;
    loader = vul_rshift(loader, 1);
    *acc4v_iter = (*acc4v_iter) + (loader & m1x4);
    ++acc4v_iter;
  }
}

HEADER_INLINE void unroll_zero_incr_1_4(uint32_t acc1_vec_ct, uintptr_t* acc1, uintptr_t* acc4) {
  const vul_t m1x4 = VCONST_UL(kMask1111);
  vul_t* acc1v_iter = (vul_t*)acc1;
  vul_t* acc4v_iter = (vul_t*)acc4;
  for (uint32_t vidx = 0; vidx < acc1_vec_ct; ++vidx) {
    vul_t loader = *acc1v_iter;
    *acc1v_iter++ = vul_setzero();
    *acc4v_iter = (*acc4v_iter) + (loader & m1x4);
    ++acc4v_iter;
    loader = vul_rshift(loader, 1);
    *acc4v_iter = (*acc4v_iter) + (loader & m1x4);
    ++acc4v_iter;
    loader = vul_rshift(loader, 1);
    *acc4v_iter = (*acc4v_iter) + (loader & m1x4);
    ++acc4v_iter;
    loader = vul_rshift(loader, 1);
    *acc4v_iter = (*acc4v_iter) + (loader & m1x4);
    ++acc4v_iter;
  }
}

// er, should this just be the same function as unroll_incr_2_4 with extra
// parameters?...
HEADER_INLINE void unroll_incr_4_8(const uintptr_t* acc4, uint32_t acc4_vec_ct, uintptr_t* acc8) {
  const vul_t m4 = VCONST_UL(kMask0F0F);
  const vul_t* acc4v_iter = (const vul_t*)acc4;
  vul_t* acc8v_iter = (vul_t*)acc8;
  for (uint32_t vidx = 0; vidx < acc4_vec_ct; ++vidx) {
    vul_t loader = *acc4v_iter++;
    *acc8v_iter = (*acc8v_iter) + (loader & m4);
    ++acc8v_iter;
    loader = vul_rshift(loader, 4);
    *acc8v_iter = (*acc8v_iter) + (loader & m4);
    ++acc8v_iter;
  }
}

HEADER_INLINE void unroll_zero_incr_4_8(uint32_t acc4_vec_ct, uintptr_t* acc4, uintptr_t* acc8) {
  const vul_t m4 = VCONST_UL(kMask0F0F);
  vul_t* acc4v_iter = (vul_t*)acc4;
  vul_t* acc8v_iter = (vul_t*)acc8;
  for (uint32_t vidx = 0; vidx < acc4_vec_ct; ++vidx) {
    vul_t loader = *acc4v_iter;
    *acc4v_iter++ = vul_setzero();
    *acc8v_iter = (*acc8v_iter) + (loader & m4);
    ++acc8v_iter;
    loader = vul_rshift(loader, 4);
    *acc8v_iter = (*acc8v_iter) + (loader & m4);
    ++acc8v_iter;
  }
}

HEADER_INLINE void unroll_incr_8_32(const uintptr_t* acc8, uint32_t acc8_vec_ct, uintptr_t* acc32) {
  const vul_t m8x32 = VCONST_UL(kMask000000FF);
  const vul_t* acc8v_iter = (const vul_t*)acc8;
  vul_t* acc32v_iter = (vul_t*)acc32;
  for (uint32_t vidx = 0; vidx < acc8_vec_ct; ++vidx) {
    vul_t loader = *acc8v_iter++;
    *acc32v_iter = (*acc32v_iter) + (loader & m8x32);
    ++acc32v_iter;
    loader = vul_rshift(loader, 8);
    *acc32v_iter = (*acc32v_iter) + (loader & m8x32);
    ++acc32v_iter;
    loader = vul_rshift(loader, 8);
    *acc32v_iter = (*acc32v_iter) + (loader & m8x32);
    ++acc32v_iter;
    loader = vul_rshift(loader, 8);
    *acc32v_iter = (*acc32v_iter) + (loader & m8x32);
    ++acc32v_iter;
  }
}

HEADER_INLINE void unroll_zero_incr_8_32(uint32_t acc8_vec_ct, uintptr_t* acc8, uintptr_t* acc32) {
  const vul_t m8x32 = VCONST_UL(kMask000000FF);
  vul_t* acc8v_iter = (vul_t*)acc8;
  vul_t* acc32v_iter = (vul_t*)acc32;
  for (uint32_t vidx = 0; vidx < acc8_vec_ct; ++vidx) {
    vul_t loader = *acc8v_iter;
    *acc8v_iter++ = vul_setzero();
    *acc32v_iter = (*acc32v_iter) + (loader & m8x32);
    ++acc32v_iter;
    loader = vul_rshift(loader, 8);
    *acc32v_iter = (*acc32v_iter) + (loader & m8x32);
    ++acc32v_iter;
    loader = vul_rshift(loader, 8);
    *acc32v_iter = (*acc32v_iter) + (loader & m8x32);
    ++acc32v_iter;
    loader = vul_rshift(loader, 8);
    *acc32v_iter = (*acc32v_iter) + (loader & m8x32);
    ++acc32v_iter;
  }
}


// advances forward_ct set bits; forward_ct must be positive.  (stays put if
// forward_ct == 1 and current bit is set.  may want to tweak this interface,
// easy to introduce off-by-one bugs...)
// In usual 64-bit case, also assumes bitvec is 16-byte aligned and the end of
// the trailing 16-byte block can be safely read from.
uintptr_t jump_forward_set_unsafe(const uintptr_t* bitvec, uintptr_t cur_pos, uintptr_t forward_ct);

// ...and here's the obvious tweaked interface.
HEADER_INLINE uint32_t idx_to_uidx_basic(const uintptr_t* bitvec, uint32_t idx) {
  return jump_forward_set_unsafe(bitvec, 0, idx + 1);
}

// variant_ct must be positive, but can be smaller than thread_ct
void compute_uidx_start_partition(const uintptr_t* variant_include, uint64_t variant_ct, uint32_t thread_ct, uint32_t first_uidx, uint32_t* variant_uidx_starts);

void compute_partition_aligned(const uintptr_t* variant_include, uint32_t orig_thread_ct, uint32_t first_variant_uidx, uint32_t cur_variant_idx, uint32_t cur_variant_ct, uint32_t alignment, uint32_t* variant_uidx_starts, uint32_t* vidx_starts);

#ifdef __arm__
  #error "Unaligned accesses in is_nan_str()."
#endif
// todo: check whether there's actually any point to the uint16_t type-pun
HEADER_INLINE uint32_t is_nan_str(const char* ss, uint32_t slen) {
  if ((slen > 3) || (slen == 1)) {
    return 0;
  }
  if (!slen) {
    return 1;
  }
  const uint32_t first_two_chars_code = ((const uint16_t*)ss)[0];
  // assumes little-endian
  if ((first_two_chars_code & 0xdfdf) != 0x414e) {
    return 0;
  }
  return (slen == 2) || ((((unsigned char)ss[2]) & 0xdf) == 78);
}

boolerr_t parse_next_range(char** argv, uint32_t param_ct, char range_delim, uint32_t* cur_param_idx_ptr, char** cur_arg_pptr, char** range_start_ptr, uint32_t* rs_len_ptr, char** range_end_ptr, uint32_t* re_len_ptr);

pglerr_t parse_name_ranges(char** argv, const char* errstr_append, uint32_t param_ct, uint32_t require_posint, char range_delim, range_list_t* range_list_ptr);


// Analytically finds all real roots of x^3 + ax^2 + bx + c, saving them in
// solutions[] (sorted from smallest to largest), and returning the count.
// Multiple roots are only returned/counted once.
uint32_t cubic_real_roots(double coef_a, double coef_b, double coef_c, double* solutions);


// For pure computations, where the launcher thread joins in as thread 0.
// threads[] is second rather than first parameter since, on Windows, we may
// need to call CloseHandle.
void join_threads(uint32_t ctp1, pthread_t* threads);

#ifndef _WIN32
extern pthread_attr_t g_smallstack_thread_attr;
#endif

boolerr_t spawn_threads(THREAD_FUNCPTR_T(start_routine), uintptr_t ct, pthread_t* threads);


// For double-buffering workloads where we don't want to respawn/join the
// threads on every block, and (unlike plink 1.9) the launcher thread does not
// participate.  (Function names end with "2z" instead of "2" since launched
// threads start with index 0 instead of 1.)
extern uintptr_t g_thread_spawn_ct;
extern uint32_t g_is_last_thread_block;

#ifdef _WIN32
extern HANDLE g_thread_start_next_event[];
extern HANDLE g_thread_cur_block_done_events[];

HEADER_INLINE void THREAD_BLOCK_FINISH(uintptr_t tidx) {
  SetEvent(g_thread_cur_block_done_events[tidx]);
  WaitForSingleObject(g_thread_start_next_event[tidx], INFINITE);
}
#else
void THREAD_BLOCK_FINISH(uintptr_t tidx);
#endif

void join_threads2z(uint32_t ct, uint32_t is_last_block, pthread_t* threads);

boolerr_t spawn_threads2z(THREAD_FUNCPTR_T(start_routine), uintptr_t ct, uint32_t is_last_block, pthread_t* threads);

// if a thread sets g_error_ret, and is_last_block wasn't true, caller should
// initialize globals to tell threads to stop, then call this function
void error_cleanup_threads2z(THREAD_FUNCPTR_T(start_routine), uintptr_t ct, pthread_t* threads);


// this interface simplifies error handling.  (todo: put most of these
// variables in a struct.)
typedef struct threads_state_struct {
  THREAD_FUNCPTR_T(thread_func_ptr);
  pthread_t* threads;
  uint32_t calc_thread_ct;
  uint32_t is_last_block;
  uint32_t is_unjoined;
} threads_state_t;

HEADER_INLINE void init_threads3z(threads_state_t* tsp) {
  tsp->thread_func_ptr = nullptr;
  tsp->threads = nullptr;
  tsp->is_last_block = 0;
  tsp->is_unjoined = 0;
}

HEADER_INLINE void reinit_threads3z(threads_state_t* tsp) {
  assert(!tsp->is_unjoined);
  tsp->is_last_block = 0;
}

HEADER_INLINE boolerr_t spawn_threads3z(uint32_t is_not_first_block, threads_state_t* tsp) {
  if (spawn_threads2z(tsp->thread_func_ptr, tsp->calc_thread_ct, tsp->is_last_block, tsp->threads)) {
    if (!is_not_first_block) {
      tsp->thread_func_ptr = nullptr;
    }
    return 1;
  }
  tsp->is_unjoined = 1;
  return 0;
}

HEADER_INLINE void join_threads3z(threads_state_t* tsp) {
  join_threads2z(tsp->calc_thread_ct, tsp->is_last_block, tsp->threads);
  tsp->is_unjoined = 0;
  if (tsp->is_last_block) {
    tsp->thread_func_ptr = nullptr;
  }
}

HEADER_INLINE void stop_threads3z(threads_state_t* tsp, uint32_t* cur_block_sizep) {
  assert(tsp->thread_func_ptr);
  if (tsp->is_unjoined) {
    join_threads2z(tsp->calc_thread_ct, tsp->is_last_block, tsp->threads);
  }
  tsp->is_unjoined = 0;
  if (!tsp->is_last_block) {
    tsp->is_last_block = 1;
    if (cur_block_sizep) {
      *cur_block_sizep = 0;
    }
    error_cleanup_threads2z(tsp->thread_func_ptr, tsp->calc_thread_ct, tsp->threads);
  }
  tsp->thread_func_ptr = nullptr;
}

HEADER_INLINE void threads3z_cleanup(threads_state_t* tsp, uint32_t* cur_block_sizep) {
  if (tsp->thread_func_ptr) {
    if (tsp->is_unjoined) {
      join_threads2z(tsp->calc_thread_ct, tsp->is_last_block, tsp->threads);
    }
    if (!tsp->is_last_block) {
      if (cur_block_sizep) {
        *cur_block_sizep = 0;
      }
      error_cleanup_threads2z(tsp->thread_func_ptr, tsp->calc_thread_ct, tsp->threads);
    }
  }
}


pglerr_t populate_id_htable_mt(const uintptr_t* subset_mask, char** item_ids, uintptr_t item_ct, uint32_t store_all_dups, uint32_t id_htable_size, uint32_t thread_ct, uint32_t* id_htable);

// pass in htable_dup_base_ptr == nullptr if not storing all duplicate IDs
pglerr_t alloc_and_populate_id_htable_mt(const uintptr_t* subset_mask, char** item_ids, uintptr_t item_ct, uint32_t max_thread_ct, uint32_t** id_htable_ptr, uint32_t** htable_dup_base_ptr, uint32_t* id_htable_size_ptr);


typedef struct help_ctrl_struct {
  uint32_t iters_left;
  uint32_t param_ct;
  char** argv;
  uintptr_t unmatched_ct;
  uintptr_t* all_match_arr;
  uintptr_t* prefix_match_arr;
  uintptr_t* perfect_match_arr;
  uint32_t* param_slens;
  uint32_t preprint_newline;
} help_ctrl_t;

void help_print(const char* cur_params, help_ctrl_t* help_ctrl_ptr, uint32_t postprint_newline, const char* payload);


extern const char errstr_nomem[];
extern const char errstr_write[];
extern const char errstr_read[];
extern const char errstr_thread_create[];

// assumes logfile is open
void disp_exit_msg(pglerr_t reterr);

boolerr_t check_extra_param(char** argv, const char* permitted_modif, uint32_t* other_idx_ptr);

char extract_char_param(const char* ss);

pglerr_t cmdline_alloc_string(const char* source, const char* flag_name, uint32_t max_slen, char** sbuf_ptr);

pglerr_t alloc_fname(const char* source, const char* flagname_p, uint32_t extra_size, char** fnbuf_ptr);

pglerr_t alloc_and_flatten(char** sources, uint32_t param_ct, uint32_t max_blen, char** flattened_buf_ptr);


typedef struct plink2_cmdline_meta_struct {
  char** subst_argv;
  char* script_buf;
  char* rerun_buf;
  char* flag_buf;
  uint32_t* flag_map;
} plink2_cmdline_meta_t;

void plink2_cmdline_meta_preinit(plink2_cmdline_meta_t* pcmp);

// Handles --script, --rerun, --help, --version, and --silent.
// subst_argv, script_buf, and rerun_buf must be initialized to nullptr.
pglerr_t cmdline_parse_phase1(const char* ver_str, const char* ver_str2, const char* prog_name_str, const char* notestr_null_calc2, const char* cmdline_format_str, const char* errstr_append, uint32_t max_flag_blen, pglerr_t(* disp_help_fn)(uint32_t, char**), int* argc_ptr, char*** argv_ptr, plink2_cmdline_meta_t* pcmp, uint32_t* first_arg_idx_ptr, uint32_t* flag_ct_ptr);

// Assumes cmdline_parse_phase1() has completed, flag names have been copied to
// flag_buf/flag_map, aliases handled, and PROG_NAME_STR has been copied to
// outname (no null-terminator needed).  outname_end must be initialized to
// nullptr.
// This sorts the flag names so they're processed in a predictable order,
// handles --out if present, initializes the log, and determines the number of
// processors the OS wants us to think the machine has.
pglerr_t cmdline_parse_phase2(const char* ver_str, const char* errstr_append, char** argv, uint32_t prog_name_str_slen, uint32_t max_flag_blen, int32_t argc, uint32_t flag_ct, plink2_cmdline_meta_t* pcmp, char* outname, char** outname_end_ptr, int32_t* known_procs_ptr, uint32_t* max_thread_ct_ptr);

HEADER_INLINE void invalid_arg(const char* cur_arg) {
  LOGPREPRINTFWW("Error: Unrecognized flag ('%s').\n", cur_arg);
}

pglerr_t cmdline_parse_phase3(uintptr_t max_default_mb, uintptr_t malloc_size_mb, uint32_t memory_require, plink2_cmdline_meta_t* pcmp, unsigned char** bigstack_ua_ptr);

void plink2_cmdline_meta_cleanup(plink2_cmdline_meta_t* pcmp);

// this is technically application-dependent, but let's keep this simple for
// now
#ifndef __LP64__
  // 2047 seems to consistently fail on both OS X and Windows
  #ifdef _WIN32
CONSTU31(kMalloc32bitMbMax, 1760);
  #else
    #ifdef __APPLE__
CONSTU31(kMalloc32bitMbMax, 1920);
    #else
CONSTU31(kMalloc32bitMbMax, 2047);
    #endif
  #endif
#endif



#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_CMDLINE_H__
