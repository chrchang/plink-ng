#ifndef __PLINK2_COMMON_H__
#define __PLINK2_COMMON_H__

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


// Resources needed across a variety of plink2 modules.  Now includes
// initialization code (init_chr_info_human, init_logfile, init_bigstack) to
// simplify inclusion in other programs.

#include "pgenlib_internal.h"

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

#ifdef __cplusplus
namespace plink2 {
#endif

#define PROG_NAME_STR "plink2"

// leave the door open to 32-bit dosages (or 31?  24?)
typedef uint16_t dosage_t;
typedef uint32_t dosage_prod_t;
#define kDosageMax (1U << (8 * sizeof(dosage_t) - 1))
CONSTU31(kDosageMid, kDosageMax / 2);
CONSTU31(kDosage4th, kDosageMax / 4);
static const double kRecipDosageMax = 0.000030517578125;
static const double kRecipDosageMid = 0.00006103515625;
static const float kRecipDosageMidf = 0.00006103515625;

// this is a bit arbitrary
CONSTU31(kMaxPhenoCt, 524287);
#define MAX_PHENO_CT_STR "524287"

// unnecessary to use e.g. (1LLU << 0), the FLAGSET64 macros should force the
// integer type to 64-bit.
FLAGSET64_DEF_START()
  kfMisc0,
  kfMiscAffection01 = (1 << 0),
  kfMiscAllowExtraChrs = (1 << 1),
  kfMiscRealRefAlleles = (1 << 2),
  kfMiscMajRef = (1 << 3),
  kfMiscMajRefForce = (1 << 4),
  kfMiscNonfounders = (1 << 5),
  kfMiscExtractRange = (1 << 6),
  kfMiscExcludeRange = (1 << 7),
  kfMiscKeepfileSid = (1 << 8),
  kfMiscRemovefileSid = (1 << 9),
  kfMiscKeepAutoconv = (1 << 10),
  kfMiscDoubleId = (1 << 11),
  kfMiscBiallelicOnly = (1 << 12),
  kfMiscBiallelicOnlyStrict = (1 << 13),
  kfMiscBiallelicOnlyList = (1 << 14),
  kfMiscExcludePvarFilterFail = (1 << 15),
  kfMiscVcfRequireGt = (1 << 16),
  kfMiscAutosomePar = (1 << 17),
  kfMiscAutosomeOnly = (1 << 18),
  kfMiscMergePar = (1 << 19),
  kfMiscAllowNoSamples = (1 << 20),
  kfMiscAllowNoVars = (1 << 21),
  kfMiscHweMidp = (1 << 22),
  kfMiscHweKeepFewhet = (1 << 23),
  kfMiscWriteSnplistZs = (1 << 24),
  kfMiscMafSucc = (1 << 25),
  kfMiscGenoDosage = (1 << 26),
  kfMiscGenoHhMissing = (1 << 27),
  kfMiscMindDosage = (1 << 28),
  kfMiscMindHhMissing = (1 << 29),
  kfMiscGenotypingRateDosage = (1 << 30),
  kfMiscSetMissingVarIds = (1LLU << 31),
  kfMiscChrOverrideCmdline = (1LLU << 32),
  kfMiscChrOverrideFile = (1LLU << 33),
  kfMiscNewVarIdOverflowMissing = (1LLU << 34),
  kfMiscNewVarIdOverflowTruncate = (1LLU << 35),
  kfMiscRequirePheno = (1LLU << 36),
  kfMiscRequireCovar = (1LLU << 37),
  kfMiscCatPhenoFamily = (1LLU << 38)
FLAGSET64_DEF_END(misc_flags_t);

FLAGSET64_DEF_START()
  kfExportf0,
  kfExportf01 = (1 << 0),
  kfExportf12 = (1 << 1),
  kfExportfSpaces = (1 << 2),
  kfExportfRefFirst = (1 << 3),
  kfExportf23 = (1 << 4),
  kfExportfA = (1 << 5),
  kfExportfATranspose = (1 << 6),
  kfExportfAD = (1 << 7),
  kfExportfBeagle = (1 << 8),
  kfExportfBeagleNomap = (1 << 9),
  kfExportfBgen11 = (1 << 10),
  kfExportfBgen12 = (1 << 11),
  kfExportfBgen13 = (1 << 12),
  kfExportfBimbam = (1 << 13),
  kfExportfBimbam1chr = (1 << 14),
  kfExportfFastphase = (1 << 15),
  kfExportfFastphase1chr = (1 << 16),
  kfExportfHaps = (1 << 17),
  kfExportfHapsLegend = (1 << 18),
  kfExportfHv = (1 << 19),
  kfExportfHv1chr = (1 << 20),
  kfExportfIndMajorBed = (1 << 21),
  kfExportfLgen = (1 << 22),
  kfExportfLgenRef = (1 << 23),
  kfExportfList = (1 << 24),
  kfExportfRlist = (1 << 25),
  kfExportfOxGen = (1 << 26),
  kfExportfPed = (1 << 27),
  kfExportfCompound = (1 << 28),
  kfExportfStructure = (1 << 29),
  kfExportfTranspose = (1 << 30),
  kfExportfVcf = (1U << 31),
  kfExportfTypemask = (2 * kfExportfVcf) - kfExportf23,
  kfExportfIncludeAlt = (1LLU << 32),
  kfExportfBgz = (1LLU << 33),
  kfExportfVcfDosageGp = (1LLU << 34),
  kfExportfVcfDosageDs = (1LLU << 35),
  kfExportfOmitNonmaleY = (1LLU << 36)
FLAGSET64_DEF_END(exportf_flags_t);

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

#ifndef _GNU_SOURCE
  #define rawmemchr(ss, cc) memchr((ss), (cc), (0x80000000U - kBytesPerVec))
#endif

#ifdef _WIN32
// if MAX_THREADS > 64, single WaitForMultipleObjects calls must be converted
// into loops
  CONSTU31(kMaxThreads, 64);
#else
// currently assumed to be less than 2^16 (otherwise some multiply overflows
// are theoretically possible, at least in the 32-bit build)
  CONSTU31(kMaxThreads, 512);
#endif

// asserts didn't seem to work properly with a setting much smaller than this
CONSTU31(kDefaultThreadStack, 131072);

// generic maximum line byte length.  .ped/.vcf/etc. lines can of course be
// longer
CONSTU31(kMaxMediumLine, 131072);

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
HEADER_INLINE uint64_t round_up_pow2_ull(uint64_t val, uint64_t alignment) {
  return round_up_pow2(val, alignment);
}
#else
HEADER_INLINE uint64_t round_up_pow2_ull(uint64_t val, uint64_t alignment) {
  uint64_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return (val + alignment_m1) & (~alignment_m1);
}
#endif

typedef struct aperm_struct {
  uint32_t min;
  uint32_t max;
  double alpha;
  double beta;
  double init_interval;
  double interval_slope;
} aperm_t;

// (2^31 - 1000001) / 2
CONSTU31(kApermMax, 1073241823);

// file-scope string constants don't always have the g_ prefix, but multi-file
// constants are always tagged.
extern const char g_errstr_fopen[];
// extern const char g_cmdline_format_str[];

extern char g_textbuf[];

extern const char* g_one_char_strs;

// '.' missing genotype value is now taken for granted; this is in *addition*
// to it (default '0').
extern const char* g_input_missing_geno_ptr;

extern const char* g_output_missing_geno_ptr; // now defaults to '.'

extern FILE* g_logfile;

// mostly-safe sprintf buffer.  warning: do NOT put allele codes or
// arbitrary-length lists in here.
extern char g_logbuf[];

extern uint32_t g_debug_on;
extern uint32_t g_log_failed;

// for --warning-errcode
extern uint32_t g_stderr_written_to;


typedef struct ll_str_struct {
  struct ll_str_struct* next;
  char ss[];
} ll_str_t;

boolerr_t push_llstr(const char* ss, ll_str_t** ll_stack_ptr);

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
// void wordwrap(uint32_t suffix_len, char* ss);

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

// fclose_null defined in pgenlib_internal.h

HEADER_INLINE void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

uint32_t int_slen(int32_t num);

// assumes it's safe to read first s_const_len bytes of s_read
int32_t strcmp_se(const char* s_read, const char* s_const, uint32_t s_const_len);

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

pglerr_t sort_cmdline_flags(uint32_t max_flag_len, uint32_t flag_ct, char* flag_buf, uint32_t* flag_map);

pglerr_t init_logfile(uint32_t always_stderr, char* outname, char* outname_end);

boolerr_t cleanup_logfile(uint32_t print_end_time);

CONSTU31(kNonBigstackMin, 67108864);

CONSTU31(kBigstackMinMb, 640);
CONSTU31(kBigstackDefaultMb, 2048);

static const double kPi = 3.1415926535897932;
static const double kRecipE = 0.36787944117144233;
static const double kRecip2m53 = 0.00000000000000011102230246251565404236316680908203125;
static const double kRecip2m32 = 0.00000000023283064365386962890625;
static const double k2m64 = 18446744073709551616.0;

// 2^{-44}
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

HEADER_INLINE boolerr_t bigstack_alloc_dosage(uintptr_t ct, dosage_t** dosage_arr_ptr) {
  *dosage_arr_ptr = (dosage_t*)bigstack_alloc(ct * sizeof(dosage_t));
  return !(*dosage_arr_ptr);
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

HEADER_INLINE boolerr_t bigstack_alloc_dosagep(uintptr_t ct, dosage_t*** dosagep_arr_ptr) {
  *dosagep_arr_ptr = (dosage_t**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*dosagep_arr_ptr);
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

HEADER_INLINE boolerr_t bigstack_end_alloc_dosage(uintptr_t ct, dosage_t** dosage_arr_ptr) {
  *dosage_arr_ptr = (dosage_t*)bigstack_end_alloc(ct * sizeof(dosage_t));
  return !(*dosage_arr_ptr);
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


typedef struct l32str_struct {
  uint32_t len;
  char ss[];
} l32str_t;


HEADER_INLINE int32_t is_letter(unsigned char ucc) {
  return (((ucc & 192) == 64) && (((ucc - 1) & 31) < 26));
}

// if we need the digit value, better to use (unsigned char)cc - '0'...
HEADER_INLINE int32_t is_digit(unsigned char ucc) {
  return (ucc <= '9') && (ucc >= '0');
}

HEADER_INLINE int32_t is_not_digit(unsigned char ucc) {
  return (ucc > '9') || (ucc < '0');
}

HEADER_INLINE int32_t is_not_nzdigit(unsigned char ucc) {
  return (ucc > '9') || (ucc <= '0');
}

// may as well treat all chars < 32, except tab, as eoln...
// kns = "known non-space" (where tab counts as a space)
/*
HEADER_INLINE int32_t is_eoln_kns(unsigned char ucc) {
  return (ucc < 32);
}
*/

HEADER_INLINE int32_t is_space_or_eoln(unsigned char ucc) {
  return (ucc <= 32);
}

HEADER_INLINE int32_t is_eoln_kns(unsigned char ucc) {
  // could assert ucc is not a space/tab?
  return (ucc <= 32);
}

HEADER_INLINE int32_t is_eoln_or_comment_kns(unsigned char ucc) {
  return (ucc < 32) || (ucc == '#');
}

HEADER_INLINE int32_t no_more_tokens_kns(const char* sptr) {
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

/*
void str_toupper(char* ss);

void buf_toupper(uint32_t slen, char* ss);

void strcpy_toupper(char* target, const char* source);
*/

uint32_t is_alphanumeric(const char* ss);

// scan_posint_capped(), scan_uint_capped(), scan_int_abs_bounded(),
// scan_int32(), scan_posint_defcap(), scan_uint_defcap(),
// scan_int_abs_defcap(), scan_uint_icap() in pgenlib_internal

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

// memcpya(), memseta() defined in pgenlib_internal.h

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

char* uint32toa(uint32_t uii, char* start);

char* int32toa(int32_t ii, char* start);

char* uitoa_z4(uint32_t uii, char* start);

char* int64toa(int64_t llii, char* start);

char* uitoa_trunc4(uint32_t uii, char* start);

char* dtoa_g(double dxx, char* start);

// We try to avoid micromanaging floating point printing and just use %g
// everywhere, but occasionally we explicitly need more precision.
//
// dtoa_g_p8 provides generic 8-digit precision (instead of %g's 6-digit
// default), while print_dosage provides up to 3 places after the decimal point
// when dealing with dosages (which are internally represented as 32768ths).
// (may want to replace _p8 with _p10 for perfect int32 handling.)
char* dtoa_g_p8(double dxx, char* start);

char* print_dosage(uint64_t dosage, char* start);

// char* dtoa_f_p5_clipped(double dxx, char* start);

char* ftoa_g(float fxx, char* start);

HEADER_INLINE char* uint32toa_x(uint32_t uii, char extra_char, char* start) {
  char* penult = uint32toa(uii, start);
  *penult = extra_char;
  return &(penult[1]);
}

// fill_uint_zero, fill_ulong_zero, fill_ull_zero, fill_ulong_one currently
// defined in pgenlib_internal.h

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
// next_set, prev_set_unsafe, are_all_words_zero defined in pgenlib_internal.h

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

// bitvec_and(), bitvec_andnot() in pgenlib_internal.h

void bitvec_and_copy(const uintptr_t* __restrict source1_bitvec, const uintptr_t* __restrict source2_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec);

void bitvec_andnot_copy(const uintptr_t* __restrict source_bitvec, const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec);

void bitvec_or(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec);

void bitvec_andnot2(const uintptr_t* __restrict include_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

void set_het_missing(uintptr_t word_ct, uintptr_t* genovec);

void genoarr_to_nonmissing(const uintptr_t* genoarr, uint32_t sample_ctl2, uintptr_t* nonmissing_bitarr);

uint32_t genoarr_count_missing_notsubset_unsafe(const uintptr_t* genoarr, const uintptr_t* exclude_mask, uint32_t sample_ct);

// dumb linear scan
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
// eventually want this to be a C++14 constexpr?
uint32_t murmurhash3_32(const void* key, uint32_t len);

// see http://lemire.me/blog/2016/06/27/a-fast-alternative-to-the-modulo-reduction/
HEADER_INLINE uint32_t hashceil(const char* idstr, uint32_t idlen, uint32_t htable_size) {
  return (((uint64_t)murmurhash3_32(idstr, idlen)) * htable_size) >> 32;
}

uintptr_t geqprime(uintptr_t floor);

// assumes ceil is odd and greater than 4.  Returns the first prime <= ceil.
uintptr_t leqprime(uintptr_t ceil);

HEADER_INLINE uint32_t get_htable_min_size(uintptr_t item_ct) {
  if (item_ct > 6) {
    return geqprime(item_ct * 2 + 1);
  }
  return 13;
}

// load factor ~20% seems to yield the best speed/space tradeoff on my test
// machines
HEADER_INLINE uint32_t get_htable_fast_size(uint32_t item_ct) {
  if (item_ct < 858993456) {
    return geqprime(round_up_pow2(item_ct * 5, 2) + 1);
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

uint32_t sid_col_required(const uintptr_t* sample_include, const char* sids, uint32_t sample_ct, uint32_t max_sid_blen, uint32_t maybe_modifier);

// forced SID '0' if sids == nullptr
// ok for sample_augid_map_ptr == nullptr
pglerr_t augid_init_alloc(const uintptr_t* sample_include, const char* sample_ids, const char* sids, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t** sample_augid_map_ptr, char** sample_augids_ptr, uintptr_t* max_sample_augid_blen_ptr);

HEADER_INLINE double get_nonmaj_freq(const double* cur_allele_freqs, uint32_t cur_allele_ct) {
  double tot_nonlast_freq = cur_allele_freqs[0];
  double max_freq = tot_nonlast_freq;
  const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
  for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct_m1; ++allele_idx) {
    const double cur_alt_freq = cur_allele_freqs[allele_idx];
    tot_nonlast_freq += cur_alt_freq;
    if (cur_alt_freq > max_freq) {
      max_freq = cur_alt_freq;
    }
  }
  const double nonmajor_freq = 1.0 - max_freq;
  return MINV(nonmajor_freq, tot_nonlast_freq);
}

HEADER_INLINE double get_allele_freq(const double* cur_allele_freqs, uint32_t allele_idx, uint32_t cur_allele_ct) {
  const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
  if (allele_idx < cur_allele_ct_m1) {
    return cur_allele_freqs[allele_idx];
  }
  double last_freq = 1.0 - cur_allele_freqs[0];
  for (uint32_t tmp_allele_idx = 1; tmp_allele_idx < cur_allele_ct_m1; ++tmp_allele_idx) {
    last_freq -= cur_allele_freqs[tmp_allele_idx];
  }
  // possible todo: force this to be nonnegative?
  return last_freq;
}


FLAGSET_DEF_START()
  kfXidMode0,
  
  kfXidModeFlagOneTokenOk = (1 << 0),
  kfXidModeFlagNeverFid = (1 << 1),
  kfXidModeFlagSid = (1 << 2),

  kfXidModeFidiid = 0,
  kfXidModeFidiidOrIid = kfXidModeFlagOneTokenOk,
  kfXidModeIid = (kfXidModeFlagOneTokenOk | kfXidModeFlagNeverFid),
  kfXidModeFidiidSid = kfXidModeFlagSid,
  kfXidModeIidSid = (kfXidModeFlagNeverFid | kfXidModeFlagSid)
FLAGSET_DEF_END(xid_mode_t);

// sample_xid_map allocated on bottom, to play well with --indiv-sort
pglerr_t sorted_xidbox_init_alloc(const uintptr_t* sample_include, const char* sample_ids, const char* sids, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, xid_mode_t xid_mode, uint32_t use_nsort, char** sorted_xidbox_ptr, uint32_t** xid_map_ptr, uintptr_t* max_xid_blen_ptr);

// returns 1 on missing token *or* if the sample ID is not present.  cases can
// be distinguished by checking whether *read_pp_new == nullptr: if it is, a
// missing-tokens error should probably be reported.
// sample_id_map == nullptr is permitted
// *read_pp is now set to point to the end of IID/SID instead of the beginning
// of the next token; this is a change from plink 1.9.
boolerr_t sorted_xidbox_read_find(const char* __restrict sorted_xidbox, const uint32_t* __restrict xid_map, uintptr_t max_xid_blen, uintptr_t end_idx, uint32_t comma_delim, xid_mode_t xid_mode, char** read_pp, uint32_t* sample_uidx_ptr, char* __restrict idbuf);

ENUM_U31_DEF_START()
  kSidDetectModeNotLoaded,
  kSidDetectModeLoaded,
  kSidDetectModeForce
ENUM_U31_DEF_END(sid_detect_mode_t);


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


HEADER_INLINE uint32_t realnum(double dd) {
  return (dd == dd) && (dd != INFINITY) && (dd != -INFINITY);
}

// note that this is no longer divisible by 64
CONSTU31(kMaxContigs, 65274);

// change chr_idx_t to uint32_t if (kMaxContigs + kChrOffsetCt) > 65536
typedef uint16_t chr_idx_t;

// get_htable_min_size(kChrRawEnd) (use constexpr once sufficient
// compiler support is available)
// (not get_htable_fast_size since, an overwhelming majority of the time, we'll
// have far fewer than 2^16 codes)
CONSTU31(kChrHtableSize, 130579);

// (note that n+1, n+2, n+3, and n+4 are reserved for X/Y/XY/MT)
CONSTU31(kMaxChrTextnum, 95);

// get_chr_code_raw() needs to be modified if this changes
CONSTU31(kMaxChrTextnumSlen, 2);

ENUM_U31_DEF_START()
  kChrOffsetX,
  kChrOffsetY,

  // old way of representing pseudo-autosomal regions.  clumsy since this
  // required changing chromosome order
  kChrOffsetXY,
  
  kChrOffsetMT,

  // plink 2.x pseudo-autosomal regions.
  kChrOffsetPAR1,
  kChrOffsetPAR2,
  kChrOffsetCt
ENUM_U31_DEF_END(xymt_offset_t);

CONSTU31(kChrRawX, kMaxContigs + kChrOffsetX);
CONSTU31(kChrRawY, kMaxContigs + kChrOffsetY);
CONSTU31(kChrRawXY, kMaxContigs + kChrOffsetXY);
CONSTU31(kChrRawMT, kMaxContigs + kChrOffsetMT);
CONSTU31(kChrRawPAR1, kMaxContigs + kChrOffsetPAR1);
CONSTU31(kChrRawPAR2, kMaxContigs + kChrOffsetPAR2);
CONSTU31(kChrRawEnd, kMaxContigs + kChrOffsetCt);

static_assert((!(kChrRawEnd % kBitsPerWord)), "kChrRawEnd expression must be updated.");
CONSTU31(kChrMaskWords, kChrRawEnd / kBitsPerWord);

#ifdef __LP64__
CONSTU31(kChrExcludeWords, 2);
#else
CONSTU31(kChrExcludeWords, 4);
#endif
static_assert(kChrExcludeWords * kBitsPerWord >= kMaxChrTextnum + 2 * kChrOffsetCt + 1, "kChrExcludeWords must be updated.");

ENUM_U31_DEF_START()
  kChrsetSourceDefault,
  kChrsetSourceCmdline,
  kChrsetSourceFile
ENUM_U31_DEF_END(chrset_source_t);

FLAGSET_DEF_START()
  kfChrOutput0,
  kfChrOutputPrefix = (1 << 0),
  kfChrOutputM = (1 << 1),
  kfChrOutputMT = (1 << 2),
  kfChrOutput0M = (1 << 3)
FLAGSET_DEF_END(chr_output_t);

typedef struct {
  // Main dynamic block intended to be allocated as a single aligned block of
  // memory on the heap freeable with vecaligned_free(), with chr_mask at the
  // base.

  uintptr_t* chr_mask; // which chromosomes aren't known to be absent?
  // This is a misnomer--it includes X and excludes MT.  Underlying concept is
  // "are some calls guaranteed to be homozygous (assuming >= 1 male)", which
  // is no longer true for MT since heteroplasmy is a thing.  (Well, the real
  // goal with MT is to enable dosage-based analysis, but until all pipelines
  // have adapted, diploid data handling loses slightly less information than
  // haploid.)
  uintptr_t* haploid_mask;

  // order of chromosomes in input files
  // currently tolerates out-of-order chromosomes, as long as all variants for
  // any given chromosome are together
  uint32_t* chr_file_order;
  
  // if the second chromosome in the dataset is chr5, chr_file_order[1] == 5,
  // the raw variant indexes for chr5 are in [chr_fo_vidx_start[1],
  // chr_fo_vidx_start[2]). and chr_idx_to_foidx[5] == 1.
  uint32_t* chr_fo_vidx_start;
  uint32_t* chr_idx_to_foidx;

  // --allow-extra-chr support
  char** nonstd_names;
  uint32_t* nonstd_id_htable;
  // end main dynamic block

  uint32_t chr_ct; // number of distinct chromosomes/contigs
  chrset_source_t chrset_source;

  uintptr_t chr_exclude[kChrExcludeWords];
  int32_t xymt_codes[kChrOffsetCt]; // X, Y, XY...; -2 = not in chromosome set
  uint32_t max_numeric_code;
  uint32_t max_code; // no longer identical to max_numeric_code, with PARs

  uint32_t autosome_ct;

  // yet more --allow-extra-chr support
  uint32_t zero_extra_chrs;
  uint32_t name_ct;
  ll_str_t* incl_excl_name_stack;
  uint32_t is_include_stack;
  chr_output_t output_encoding;
} chr_info_t;

extern const char g_xymt_log_names[][5];

pglerr_t init_chr_info(chr_info_t* cip);

void finalize_chrset(misc_flags_t misc_flags, chr_info_t* cip);

HEADER_INLINE pglerr_t init_chr_info_human(chr_info_t* cip) {
  // convenience wrapper
  if (init_chr_info(cip)) {
    return kPglRetNomem;
  }
  finalize_chrset(kfMisc0, cip);
  return kPglRetSuccess;
}

void forget_extra_chr_names(uint32_t reinitialize, chr_info_t* cip);

// in the usual case where the number of chromosomes/contigs is much less than
// kMaxContigs, this reduces chr_info's memory consumption and improves
// locality.
pglerr_t finalize_chr_info(chr_info_t* cip);

void cleanup_chr_info(chr_info_t* cip);

char* chr_name_write(const chr_info_t* cip, uint32_t chr_idx, char* buf);

uint32_t get_max_chr_slen(const chr_info_t* cip);

uint32_t haploid_chr_present(const chr_info_t* cip);

// any character <= ' ' is considered a terminator
// maps chrX -> kChrRawX, etc.
int32_t get_chr_code_raw(const char* sptr);

// requires chr_name to be null-terminated
// maps chrX -> xymt_codes[kChrOffsetX], etc.
// error codes:
//   -1 = --allow-extra-chr ok
//   -2 = total fail
int32_t get_chr_code(const char* chr_name, const chr_info_t* cip, uint32_t name_slen);

// when the chromosome name isn't null-terminated
// requires chr_name[name_slen] to be mutable
int32_t get_chr_code_counted(const chr_info_t* cip, uint32_t name_slen, char* chr_name);

HEADER_INLINE uint32_t get_variant_chr_fo_idx(const chr_info_t* cip, uintptr_t variant_uidx) {
  return uint32arr_greater_than(&(cip->chr_fo_vidx_start[1]), cip->chr_ct, variant_uidx + 1);
}

HEADER_INLINE uint32_t get_variant_chr(const chr_info_t* cip, uintptr_t variant_uidx) {
  return cip->chr_file_order[get_variant_chr_fo_idx(cip, variant_uidx)];
}

HEADER_INLINE uint32_t xymt_exists(const chr_info_t* cip, uint32_t xymt_offset, int32_t* xymt_code_ptr) {
  // too easy to forget is_set(chr_mask) check if we don't use this
  const int32_t xymt_code = cip->xymt_codes[xymt_offset];
  *xymt_code_ptr = xymt_code;
  return (xymt_code >= 0) && is_set(cip->chr_mask, xymt_code);
}

HEADER_INLINE void get_xymt_start_and_end(const chr_info_t* cip, uint32_t xymt_offset, uint32_t* xymt_start_ptr, uint32_t* xymt_end_ptr) {
  int32_t xymt_code;
  if (!xymt_exists(cip, xymt_offset, &xymt_code)) {
    *xymt_start_ptr = 0;
    *xymt_end_ptr = 0;
    return;
  }
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)xymt_code];
  *xymt_start_ptr = cip->chr_fo_vidx_start[chr_fo_idx];
  *xymt_end_ptr = cip->chr_fo_vidx_start[chr_fo_idx + 1];
}

HEADER_INLINE void get_xymt_code_start_and_end_unsafe(const chr_info_t* cip, uint32_t xymt_offset, int32_t* xymt_code_ptr, uint32_t* xymt_start_ptr, uint32_t* xymt_end_ptr) {
  // assumes xymt_exists was previously called, and is true
  const int32_t xymt_code = cip->xymt_codes[xymt_offset];
  *xymt_code_ptr = xymt_code;
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)xymt_code];
  *xymt_start_ptr = cip->chr_fo_vidx_start[chr_fo_idx];
  *xymt_end_ptr = cip->chr_fo_vidx_start[chr_fo_idx + 1];
}

// now assumes chr_name is null-terminated
pglerr_t try_to_add_chr_name(const char* chr_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t allow_extra_chrs, int32_t* chr_idx_ptr, chr_info_t* cip);

HEADER_INLINE pglerr_t get_or_add_chr_code(const char* chr_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t allow_extra_chrs, chr_info_t* cip, int32_t* chr_idx_ptr) {
  *chr_idx_ptr = get_chr_code(chr_name, cip, name_slen);
  if (*chr_idx_ptr >= 0) {
    return kPglRetSuccess;
  }
  return try_to_add_chr_name(chr_name, file_descrip, line_idx, name_slen, allow_extra_chrs, chr_idx_ptr, cip);
}

HEADER_INLINE pglerr_t get_or_add_chr_code_destructive(const char* file_descrip, uintptr_t line_idx, uint32_t allow_extra_chrs, char* chr_name, char* chr_name_end, chr_info_t* cip, int32_t* chr_idx_ptr) {
  *chr_name_end = '\0';
  return get_or_add_chr_code(chr_name, file_descrip, line_idx, (uintptr_t)(chr_name_end - chr_name), allow_extra_chrs, cip, chr_idx_ptr);
}

#ifdef __LP64__
HEADER_INLINE uintptr_t popcount_longs_nzbase(const uintptr_t* bitvec, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t prefix_ct = 0;
  #ifdef USE_AVX2
  while (start_idx & 3) {
    if (end_idx == start_idx) {
      return 0;
    }
    prefix_ct = popcount_long(bitvec[start_idx++]);
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

// uintptr_t count_11_longs(const uintptr_t* genovec, uintptr_t word_ct);

uint32_t are_all_bits_zero(const uintptr_t* bitarr, uintptr_t start_idx, uintptr_t end_idx);

// assumes len is positive, and relevant bits of target_bitarr are zero
void copy_bitarr_range(const uintptr_t* __restrict src_bitarr, uintptr_t src_start_bitidx, uintptr_t target_start_bitidx, uintptr_t len, uintptr_t* __restrict target_bitarr);

// zeroes out samples not in the mask
void interleaved_mask_zero(const uintptr_t* __restrict interleaved_mask, uintptr_t vec_ct, uintptr_t* __restrict genovec);

// sets samples in the mask to missing (0b11)
void interleaved_set_missing(const uintptr_t* __restrict interleaved_set, uintptr_t vec_ct, uintptr_t* __restrict genovec);

void set_male_het_missing(const uintptr_t* __restrict sex_male_interleaved, uint32_t vec_ct, uintptr_t* __restrict genovec);

// Clears each bit in bitarr which doesn't correspond to a genovec het.
// Assumes that either trailing bits of bitarr are already zero, or trailing
// bits of genovec are zero.
void mask_genovec_hets_unsafe(const uintptr_t* __restrict genovec, uint32_t raw_sample_ctl2, uintptr_t* __restrict bitarr);

// vertical popcount support
#ifdef __LP64__
static_assert(kBytesPerVec == 16, "scramble_2_4_8_32() assumes kBytesPerVec == 16.");
HEADER_INLINE uint32_t scramble_2_4_8_32(uint32_t orig_idx) {
  return (orig_idx & (~63)) + ((orig_idx & 1) * 32) + ((orig_idx & 2) * 8) + (orig_idx & 12) + ((orig_idx & 48) / 16);
}
#else
// 2->4: 0 2 4 6 8 10 12 14 1 3 5 ...
// 4->8: 0 4 8 12 2 6 10 14 1 5 9 ...
// 8->32: 0 4 8 12 2 6 10 14 1 5 9 13 3 7 11 15
HEADER_INLINE uint32_t scramble_2_4_8_32(uint32_t orig_idx) {
  return (orig_idx & (~15)) + ((orig_idx & 1) * 8) + ((orig_idx & 2) * 2) + ((orig_idx & 12) / 4);
}
#endif

// probable todo: switch to vul_t* parameters
HEADER_INLINE void unroll_incr_2_4(const uintptr_t* acc2, uint32_t acc2_vec_ct, uintptr_t* acc4) {
  const vul_t m2 = VCONST_UL(kMask3333);
  const vul_t* acc2v_iter = (const vul_t*)acc2;
  vul_t* acc4v_iter = (vul_t*)acc4;
  for (uint32_t vidx = 0; vidx < acc2_vec_ct; ++vidx) {
    vul_t loader = *acc2v_iter++;
    *acc4v_iter = (*acc4v_iter) + (loader & m2);
    ++acc4v_iter;
    loader = vul_rshift(loader, 2);
    *acc4v_iter = (*acc4v_iter) + (loader & m2);
    ++acc4v_iter;
  }
}

HEADER_INLINE void unroll_zero_incr_2_4(uint32_t acc2_vec_ct, uintptr_t* acc2, uintptr_t* acc4) {
  const vul_t m2 = VCONST_UL(kMask3333);
  vul_t* acc2v_iter = (vul_t*)acc2;
  vul_t* acc4v_iter = (vul_t*)acc4;
  for (uint32_t vidx = 0; vidx < acc2_vec_ct; ++vidx) {
    vul_t loader = *acc2v_iter;
    *acc2v_iter++ = vul_setzero();
    *acc4v_iter = (*acc4v_iter) + (loader & m2);
    ++acc4v_iter;
    loader = vul_rshift(loader, 2);
    *acc4v_iter = (*acc4v_iter) + (loader & m2);
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

#ifdef __LP64__
static_assert(kBytesPerVec == 16, "scramble_1_4_8_32() assumes kBytesPerVec == 16.");
HEADER_INLINE uint32_t scramble_1_4_8_32(uint32_t orig_idx) {
  // 1->4: 0 4 8 12 16 20 24 28 32 ... 124 1 5 9 ...
  // 4->8: 0 8 16 24 32 ... 120 4 12 20 ... 1 9 17 ...
  // 8->32: 0 32 64 96 8 40 72 104 16 48 80 112 24 56 88 120 4 36 68 ... 1 33 ...
  return (orig_idx & (~127)) + ((orig_idx & 3) * 32) + ((orig_idx & 4) * 4) + ((orig_idx & 24) / 2) + ((orig_idx & 96) / 32);
}
#else
// 1->4: 0 4 8 12 16 20 24 28 1 5 9 13 17 21 25 29 2 6 10 ...
// 4->8: 0 8 16 24 4 12 20 28 1 9 17 25 5 13 21 29 2 10 18 ...
// 8->32: 0 8 16 24 4 12 20 28 1 9 17 25 5 13 21 29 2 10 18 ...
HEADER_INLINE uint32_t scramble_1_4_8_32(uint32_t orig_idx) {
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


// uint32_t chr_window_max(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_pos, uint32_t chr_fo_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max);

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

HEADER_INLINE uint32_t count_chr_variants_unsafe(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t chr_idx) {
  assert(is_set(cip->chr_mask, chr_idx));
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
  const uint32_t min_idx = cip->chr_fo_vidx_start[chr_fo_idx];
  const uint32_t max_idx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  return popcount_bit_idx(variant_include, min_idx, max_idx);
}

HEADER_INLINE uint32_t chr_is_nonempty(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t chr_idx) {
  if (!is_set(cip->chr_mask, chr_idx)) {
    return 0;
  }
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
  const uint32_t min_idx = cip->chr_fo_vidx_start[chr_fo_idx];
  const uint32_t max_idx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  return !are_all_bits_zero(variant_include, min_idx, max_idx);
}

HEADER_INLINE uint32_t xymt_is_nonempty(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t xymt_offset) {
  const int32_t xymt_code = cip->xymt_codes[xymt_offset];
  if ((xymt_code < 0) || (!is_set(cip->chr_mask, xymt_code))) {
    return 0;
  }
  return chr_is_nonempty(variant_include, cip, (uint32_t)xymt_code);
}

// assumes there's at least one variant on specified chromosome
uint32_t not_only_xymt(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t raw_variant_ct, uint32_t xymt_offset);

uint32_t count_non_autosomal_variants(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t count_x, uint32_t count_mt);

pglerr_t conditional_allocate_non_autosomal_variants(const chr_info_t* cip, const char* calc_descrip, uint32_t raw_variant_ct, uintptr_t** variant_include_ptr, uint32_t* variant_ct_ptr);

void fill_subset_chr_fo_vidx_start(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t* subset_chr_fo_vidx_start);

HEADER_INLINE boolerr_t alloc_and_fill_subset_chr_fo_vidx_start(const uintptr_t* variant_include, const chr_info_t* cip, uint32_t** subset_chr_fo_vidx_start_ptr) {
  const uint32_t chr_ct = cip->chr_ct;
  if (bigstack_alloc_ui(chr_ct + 1, subset_chr_fo_vidx_start_ptr)) {
    return 1;
  }
  fill_subset_chr_fo_vidx_start(variant_include, cip, *subset_chr_fo_vidx_start_ptr);
  return 0;
}

// newval does not need to be null-terminated
// assumes *allele_ptr is not initialized
boolerr_t allele_set(const char* newval, uint32_t allele_slen, char** allele_ptr);

// *allele_ptr must be initialized; frees *allele_ptr if necessary
boolerr_t allele_reset(const char* newval, uint32_t allele_slen, char** allele_ptr);

void cleanup_allele_storage(uint32_t max_allele_slen, uintptr_t allele_storage_entry_ct, char** allele_storage);

CONSTU31(kMaxMissingPhenostrBlen, 32);
// might want g_input_missing_catname and/or g_output_missing_catname later,
// but let's start with the simplest implementation
extern char g_missing_catname[]; // default "NONE", not changeable for now

extern char g_output_missing_pheno[]; // default "NA"
extern char g_legacy_output_missing_pheno[]; // default "-9"

// don't care about kfUnsortedVarChrom
FLAGSET_DEF_START()
  kfUnsortedVar0,
  kfUnsortedVarBp = (1 << 0),
  kfUnsortedVarCm = (1 << 1),
  kfUnsortedVarSplitChr = (1 << 2)
FLAGSET_DEF_END(unsorted_var_t);

FLAGSET_DEF_START()
  kfFamCol0,
  kfFamCol1 = (1 << 0),
  kfFamCol34 = (1 << 1),
  kfFamCol5 = (1 << 2),
  kfFamCol6 = (1 << 3),
  kfFamCol13456 = (kfFamCol1 | kfFamCol34 | kfFamCol5 | kfFamCol6)
FLAGSET_DEF_END(fam_col_t);

HEADER_INLINE char sexchar(const uintptr_t* sex_nm, const uintptr_t* sex_male, uintptr_t sample_uidx) {
  if (is_set(sex_nm, sample_uidx)) {
    return '2' - is_set(sex_male, sample_uidx);
  }
  return '0';
}

// kPhenoDtypeCc and kPhenoDtypeQt currently can't change
// kPhenoDtypeOther currently used for --glm local covariates
ENUM_U31_DEF_START()
  kPhenoDtypeCc,
  kPhenoDtypeQt,
  kPhenoDtypeCat,
  kPhenoDtypeOther
ENUM_U31_DEF_END(pheno_dtype_t);

typedef union {
  uintptr_t* cc; // bitvector
  double* qt;
  uint32_t* cat; // always 0 for missing, nonmiss[] check unnecessary
} phenodata_t;

typedef struct {
  // * If categorical phenotype, [0] points to g_missing_catname, while [1],
  //   [2], etc. point to category names.  These are part of the same
  //   allocation as nonmiss, so no separate free is needed.
  //   Otherwise, this is nullptr.
  // * When .sample categorical variables are imported, 'P' is added in front
  //   of the integers.
  char** category_names;
  
  uintptr_t* nonmiss; // bitvector

  // essentially a tagged union; part of the same allocation as nonmiss
  phenodata_t data;
  pheno_dtype_t type_code;
  
  uint32_t nonnull_category_ct;
} pheno_col_t;

void init_pheno();


uint32_t is_categorical_phenostr(const char* phenostr);

uint32_t is_categorical_phenostr_nocsv(const char* phenostr);

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

// returns 0xffffffffU if none exists
uint32_t first_cc_or_qt_pheno_idx(const pheno_col_t* pheno_cols, uint32_t pheno_ct);

// "_covar" since this doesn't handle case/control
uint32_t is_const_covar(const pheno_col_t* covar_col, const uintptr_t* sample_include, uint32_t sample_ct);

uint32_t identify_remaining_cats(const uintptr_t* sample_include, const pheno_col_t* covar_col, uint32_t sample_ct, uintptr_t* cat_covar_wkspace);

// pheno_names is also allocated on the heap, but it can be handled with a
// simple free_cond().
void cleanup_pheno_cols(uint32_t pheno_ct, pheno_col_t* pheno_cols);

pglerr_t parse_chr_ranges(const char* flagname_p, const char* errstr_append, uint32_t param_ct, uint32_t allow_extra_chrs, uint32_t xymt_subtract, char range_delim, char** argv, chr_info_t* cip, uintptr_t* chr_mask);

pglerr_t parse_name_ranges(char** argv, const char* errstr_append, uint32_t param_ct, uint32_t require_posint, char range_delim, range_list_t* range_list_ptr);


// For pure computations, where the launcher thread joins in as thread 0.
// threads[] is second rather than first parameter since, on Windows, we may
// need to call CloseHandle.
void join_threads(uint32_t ctp1, pthread_t* threads);

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
extern pthread_attr_t g_smallstack_thread_attr;

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

// sample_ct not relevant if genovecs_ptr == nullptr
pglerr_t multithread_load_init(const uintptr_t* variant_include, uint32_t sample_ct, uint32_t variant_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t thread_xalloc_cacheline_ct, uintptr_t per_variant_xalloc_byte_ct, pgen_file_info_t* pgfip, uint32_t* calc_thread_ct_ptr, uintptr_t*** genovecs_ptr, uintptr_t*** dosage_present_ptr, dosage_t*** dosage_val_bufs_ptr, uint32_t* read_block_size_ptr, unsigned char** main_loadbufs, pthread_t** threads_ptr, pgen_reader_t*** pgr_pps, uint32_t** read_variant_uidx_starts_ptr);

pglerr_t write_sample_ids(const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* outname, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen);

#ifdef __cplusplus
} // namespace plink2
#endif
 
#endif // __PLINK2_COMMON_H__
