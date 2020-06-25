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


#include "plink2_cmdline.h"

#ifdef __APPLE__
  // needed for sysctl() call
#  include <sys/sysctl.h>
#endif

#include <sys/types.h>  // open()
#include <sys/stat.h>  // open()
#include <fcntl.h>  // open()
#include <time.h>  // time(), ctime()
#include <unistd.h>  // getcwd(), gethostname(), sysconf(), fstat()

#ifdef __cplusplus
namespace plink2 {
#endif

const char kErrprintfFopen[] = "Error: Failed to open %s : %s.\n";
const char kErrprintfFread[] = "Error: %s read failure: %s.\n";
const char kErrprintfRewind[] = "Error: %s could not be scanned twice. (Process-substitution/named-pipe input is not permitted in this use case.)\n";

char g_textbuf[kTextbufSize];

// now initialized by init_bigstack
// can't use std::array<char, 512>& since it's initially null
const char* g_one_char_strs = nullptr;
// If one-base indels become sufficiently common, might want to predefine
// g_two_char_strs[], and update allele string construction/destruction
// accordingly.  (Though that should either be programmatically initialized, or
// only cover a subset of the space; 192k is a lot to increase the binary image
// size for a single simple table.)

const char* g_input_missing_geno_ptr = nullptr;  // in addition to '.'
const char* g_output_missing_geno_ptr = nullptr;  // now '.'

FILE* g_logfile = nullptr;

char g_logbuf[kLogbufSize];

uint32_t g_debug_on = 0;
uint32_t g_log_failed = 0;
uint32_t g_stderr_written_to = 0;

void logputs_silent(const char* str) {
  if (!g_debug_on) {
    fputs(str, g_logfile);
    if (unlikely(ferror_unlocked(g_logfile))) {
      putchar('\n');
      fflush(stdout);
      fprintf(stderr, "Warning: Logging failure on:\n%s\nFurther logging will not be attempted in this run.\n", str);
      g_log_failed = 1;
    }
  } else {
    if (g_log_failed) {
      fflush(stdout);
      fputs(str, stderr);
    } else {
      fputs(str, g_logfile);
      if (unlikely(ferror_unlocked(g_logfile))) {
        putchar('\n');
        fflush(stdout);
        fprintf(stderr, "Error: Debug logging failure.  Dumping to stderr:\n%s", str);
        g_log_failed = 1;
        g_stderr_written_to = 1;
      } else {
        fflush(g_logfile);
      }
    }
  }
}

void logputs(const char* str) {
  logputs_silent(str);
  fputs(str, stdout);
}

void logerrputs(const char* str) {
  logputs_silent(str);
  fflush(stdout);
  fputs(str, stderr);
  g_stderr_written_to = 1;
}

void logputsb() {
  logputs_silent(g_logbuf);
  fputs(g_logbuf, stdout);
}

void logerrputsb() {
  logputs_silent(g_logbuf);
  fflush(stdout);
  fputs(g_logbuf, stderr);
  g_stderr_written_to = 1;
}

PglErr ForceNonFifo(const char* fname) {
  int32_t file_handle = open(fname, O_RDONLY);
  if (unlikely(file_handle < 0)) {
    return kPglRetOpenFail;
  }
  struct stat statbuf;
  if (unlikely(fstat(file_handle, &statbuf) < 0)) {
    close(file_handle);
    return kPglRetOpenFail;
  }
  if (unlikely(S_ISFIFO(statbuf.st_mode))) {
    close(file_handle);
    return kPglRetRewindFail;
  }
  close(file_handle);
  return kPglRetSuccess;
}

BoolErr fopen_checked(const char* fname, const char* mode, FILE** target_ptr) {
  /*
  if (!strcmp(mode, FOPEN_WB)) {
    // The fopen call may take a long time to return when overwriting in mode
    // "w" in some scenarios on some OSes (I've seen 15+ sec).  On the other
    // hand, sometimes it returns immediately and then writing takes less time.
    // Since it goes both ways, I won't second-guess the OSes for now;
    // presumably they're trying to minimize average amortized time.  But it
    // might be worth retesting in the future; perhaps some common plink
    // workflows violate normally-reasonable OS assumptions.
    if (access(fname, W_OK) != -1) {
      if (unlikely(unlink(fname))) {
        logputs("\n");
        logerrprintfww(kErrprintfFopen, fname, strerror(errno));
        return 1;
      }
    }
  }
  */
  *target_ptr = fopen(fname, mode);
  if (unlikely(!(*target_ptr))) {
    logputs("\n");
    logerrprintfww(kErrprintfFopen, fname, strerror(errno));
    return 1;
  }
  return 0;
}

BoolErr fwrite_flush2(char* buf_flush, FILE* outfile, char** write_iter_ptr) {
  char* buf = &(buf_flush[-S_CAST(int32_t, kMaxMediumLine)]);
  char* buf_end = *write_iter_ptr;
  *write_iter_ptr = buf;
  return fwrite_checked(buf, buf_end - buf, outfile);
}

BoolErr fclose_flush_null(char* buf_flush, char* write_iter, FILE** outfile_ptr) {
  char* buf = &(buf_flush[-S_CAST(int32_t, kMaxMediumLine)]);
  if (write_iter != buf) {
    if (unlikely(fwrite_checked(buf, write_iter - buf, *outfile_ptr))) {
      return 1;
    }
  }
  return fclose_null(outfile_ptr);
}


int32_t float_cmp(const void* aa, const void* bb) {
  const float fxx = *S_CAST(const float*, aa);
  const float fyy = *S_CAST(const float*, bb);
  if (fxx < fyy) {
    return -1;
  }
  return (fxx > fyy);
}

int32_t double_cmp(const void* aa, const void* bb) {
  const double dxx = *S_CAST(const double*, aa);
  const double dyy = *S_CAST(const double*, bb);
  if (dxx < dyy) {
    return -1;
  }
  return (dxx > dyy);
}

int32_t double_cmp_decr(const void* aa, const void* bb) {
  const double dxx = *S_CAST(const double*, aa);
  const double dyy = *S_CAST(const double*, bb);
  if (dxx > dyy) {
    return -1;
  }
  return (dxx < dyy);
}

int32_t u64cmp(const void* aa, const void* bb) {
  const uint64_t ullaa = *S_CAST(const uint64_t*, aa);
  const uint64_t ullbb = *S_CAST(const uint64_t*, bb);
  if (ullaa < ullbb) {
    return -1;
  }
  return (ullaa > ullbb);
}

#ifndef __cplusplus
int32_t u64cmp_decr(const void* aa, const void* bb) {
  const uint64_t ullaa = *S_CAST(const uint64_t*, aa);
  const uint64_t ullbb = *S_CAST(const uint64_t*, bb);
  if (ullaa > ullbb) {
    return -1;
  }
  return (ullaa < ullbb);
}
#endif

#ifdef __cplusplus
float DestructiveMedianF(uintptr_t len, float* unsorted_arr) {
  if (!len) {
    return 0.0;
  }
  const uintptr_t len_d2 = len / 2;
  std::nth_element(unsorted_arr, &(unsorted_arr[len_d2]), &(unsorted_arr[len]));
  const float median_upper = unsorted_arr[len_d2];
  if (len % 2) {
    return median_upper;
  }
  return (FArrMax(unsorted_arr, len_d2) + median_upper) * 0.5f;
}

double DestructiveMedianD(uintptr_t len, double* unsorted_arr) {
  if (!len) {
    return 0.0;
  }
  const uintptr_t len_d2 = len / 2;
  std::nth_element(unsorted_arr, &(unsorted_arr[len_d2]), &(unsorted_arr[len]));
  const double median_upper = unsorted_arr[len_d2];
  if (len % 2) {
    return median_upper;
  }
  return (DArrMax(unsorted_arr, len_d2) + median_upper) * 0.5;
}
#else
// these will probably be used in __cplusplus case too
float MedianF(const float* sorted_arr, uintptr_t len) {
  if (!len) {
    return 0.0;
  }
  if (len % 2) {
    return sorted_arr[len / 2];
  }
  return (sorted_arr[len / 2] + sorted_arr[(len / 2) - 1]) * 0.5f;
}

double MedianD(const double* sorted_arr, uintptr_t len) {
  if (!len) {
    return 0.0;
  }
  if (len % 2) {
    return sorted_arr[len / 2];
  }
  return (sorted_arr[len / 2] + sorted_arr[(len / 2) - 1]) * 0.5;
}

float DestructiveMedianF(uintptr_t len, float* unsorted_arr) {
  // no, I'm not gonna bother reimplementing introselect just for folks who
  // insist on compiling this as pure C instead of C++
  qsort(unsorted_arr, len, sizeof(float), float_cmp);
  return MedianF(unsorted_arr, len);
}

double DestructiveMedianD(uintptr_t len, double* unsorted_arr) {
  qsort(unsorted_arr, len, sizeof(double), double_cmp);
  return MedianD(unsorted_arr, len);
}
#endif


// Overread must be ok.
BoolErr SortStrboxIndexed(uintptr_t str_ct, uintptr_t max_str_blen, uint32_t use_nsort, char* strbox, uint32_t* id_map) {
  if (str_ct < 2) {
    return 0;
  }
  unsigned char* bigstack_mark = g_bigstack_base;
  const uintptr_t wkspace_entry_blen = GetStrboxsortWentryBlen(max_str_blen);
  unsigned char* sort_wkspace;
  if (unlikely(bigstack_alloc_uc(str_ct * wkspace_entry_blen, &sort_wkspace))) {
    return 1;
  }
  SortStrboxIndexed2(str_ct, max_str_blen, use_nsort, strbox, id_map, sort_wkspace);
  BigstackReset(bigstack_mark);
  return 0;
}

/*
BoolErr StrptrArrIndexedSort(const char* const* unsorted_strptrs, uint32_t str_ct, uint32_t overread_ok, uint32_t use_nsort, uint32_t* id_map) {
  if (str_ct < 2) {
    if (str_ct) {
      id_map[0] = 0;
    }
    return 0;
  }
  if (bigstack_left() < str_ct * sizeof(StrSortIndexedDeref)) {
    return 1;
  }
  StrSortIndexedDeref* wkspace_alias = (StrSortIndexedDeref*)g_bigstack_base;
  for (uint32_t str_idx = 0; str_idx != str_ct; ++str_idx) {
    wkspace_alias[str_idx].strptr = unsorted_strptrs[str_idx];
    wkspace_alias[str_idx].orig_idx = str_idx;
  }
  StrptrArrSortMain(str_ct, overread_ok, use_nsort, wkspace_alias);
  for (uint32_t str_idx = 0; str_idx != str_ct; ++str_idx) {
    id_map[str_idx] = wkspace_alias[str_idx].orig_idx;
  }
  BigstackReset(wkspace_alias);
  return 0;
}
*/


uint32_t CountSortedSmallerU32(const uint32_t* sorted_u32_arr, uint32_t arr_length, uint32_t needle) {
  // (strangely, this seems to be equal to or better than std::lower_bound with
  // -O2 optimization, but can become much slower with -O3?)

  // assumes arr_length is nonzero, and sorted_u32_arr is in nondecreasing
  // order.  (useful for searching variant_bps[].)
  // needle guaranteed to be larger than sorted_u32_arr[min_idx - 1] if it
  // exists, but NOT necessarily sorted_u32_arr[min_idx].
  int32_t min_idx = 0;
  // similarly, needle guaranteed to be no greater than
  // sorted_u32_arr[max_idx + 1] if it exists, but not necessarily
  // sorted_u32_arr[max_idx].  Signed integer since it could become -1, and
  // min_idx in turn is signed so comparisons are safe.
  int32_t max_idx = arr_length - 1;
  while (min_idx < max_idx) {
    const uint32_t mid_idx = (S_CAST(uint32_t, min_idx) + S_CAST(uint32_t, max_idx)) / 2;
    if (needle > sorted_u32_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  return min_idx + (needle > sorted_u32_arr[S_CAST(uint32_t, min_idx)]);
}

uint32_t ExpsearchU32(const uint32_t* sorted_u32_arr, uint32_t end_idx, uint32_t needle) {
  uint32_t next_incr = 1;
  uint32_t start_idx = 0;
  uint32_t cur_idx = 0;
  while (cur_idx < end_idx) {
    if (sorted_u32_arr[cur_idx] >= needle) {
      end_idx = cur_idx;
      break;
    }
    start_idx = cur_idx + 1;
    cur_idx += next_incr;
    next_incr *= 2;
  }
  while (start_idx < end_idx) {
    // this breaks if arr_length > 2^31
    const uint32_t mid_idx = (start_idx + end_idx) / 2;

    if (sorted_u32_arr[mid_idx] < needle) {
      start_idx = mid_idx + 1;
    } else {
      end_idx = mid_idx;
    }
  }
  return start_idx;
}

uintptr_t CountSortedSmallerU64(const uint64_t* sorted_u64_arr, uintptr_t arr_length, uint64_t needle) {
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  while (min_idx < max_idx) {
    const uintptr_t mid_idx = (S_CAST(uintptr_t, min_idx) + S_CAST(uintptr_t, max_idx)) / 2;
    if (needle > sorted_u64_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  return min_idx + (needle > sorted_u64_arr[S_CAST(uintptr_t, min_idx)]);
}

uintptr_t CountSortedSmallerD(const double* sorted_dbl_arr, uintptr_t arr_length, double needle) {
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  while (min_idx < max_idx) {
    const uintptr_t mid_idx = (S_CAST(uintptr_t, min_idx) + S_CAST(uintptr_t, max_idx)) / 2;
    if (needle > sorted_dbl_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  return min_idx + (needle > sorted_dbl_arr[S_CAST(uintptr_t, min_idx)]);
}

uintptr_t CountSortedLeqU64(const uint64_t* sorted_u64_arr, uintptr_t arr_length, uint64_t needle) {
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  while (min_idx < max_idx) {
    const uintptr_t mid_idx = (S_CAST(uintptr_t, min_idx) + S_CAST(uintptr_t, max_idx)) / 2;
    if (needle >= sorted_u64_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  return min_idx + (needle >= sorted_u64_arr[S_CAST(uintptr_t, min_idx)]);
}

uint32_t GetParamCt(const char* const* argvk, uint32_t argc, uint32_t flag_idx) {
  // Counts the number of optional arguments given to the flag at position
  // flag_idx, treating any nonnumeric argument beginning with "-" as optional.
  ++flag_idx;
  uint32_t cur_idx = flag_idx;
  while ((cur_idx < argc) && (!IsCmdlineFlag(argvk[cur_idx]))) {
    ++cur_idx;
  }
  return cur_idx - flag_idx;
}

BoolErr EnforceParamCtRange(const char* flag_name, uint32_t param_ct, uint32_t min_ct, uint32_t max_ct) {
  if (unlikely(param_ct > max_ct)) {
    if (max_ct > min_ct) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s accepts at most %u argument%s.\n", flag_name, max_ct, (max_ct == 1)? "" : "s");
    } else {
      snprintf(g_logbuf, kLogbufSize, "Error: %s only accepts %u argument%s.\n", flag_name, max_ct, (max_ct == 1)? "" : "s");
    }
    return 1;
  }
  if (likely(param_ct >= min_ct)) {
    return 0;
  }
  if (min_ct == 1) {
    snprintf(g_logbuf, kLogbufSize, "Error: Missing %s argument.\n", flag_name);
  } else {
    snprintf(g_logbuf, kLogbufSize, "Error: %s requires %s%u arguments.\n", flag_name, (min_ct < max_ct)? "at least " : "", min_ct);
  }
  return 1;
}

PglErr SortCmdlineFlags(uint32_t max_flag_blen, uint32_t flag_ct, char* flag_buf, uint32_t* flag_map) {
  // Assumes flag_ct is the number of flag (as opposed to value) arguments,
  // flag_buf[] points to a rectangular char* array (width max_flag_blen) of
  // flag names with leading dash(es) stripped, and flag_map[] maps flag_buf[]
  // entries to argv[] entries.
  // Lexicographically sorts flag_buf (updating flag_map in the process), and
  // then checks for duplicates.
  // Okay for flag_buf to contain entries with spaces (plink 1.9's alias
  // resolution takes advantage of this).
  assert(flag_ct);  // this must be skipped if there are no flags at all
  if (unlikely(SortStrboxIndexedMalloc(flag_ct, max_flag_blen, flag_buf, flag_map))) {
    return kPglRetNomem;
  }
  uint32_t prev_flag_len = strlen_se(flag_buf);
  char* prev_flag_ptr = flag_buf;
  for (uint32_t cur_flag_idx = 1; cur_flag_idx != flag_ct; ++cur_flag_idx) {
    char* cur_flag_ptr = &(prev_flag_ptr[max_flag_blen]);
    const uint32_t cur_flag_len = strlen_se(cur_flag_ptr);
    if (unlikely((prev_flag_len == cur_flag_len) && memequal(prev_flag_ptr, cur_flag_ptr, cur_flag_len))) {
      cur_flag_ptr[cur_flag_len] = '\0';  // just in case of aliases
      fflush(stdout);
      fprintf(stderr, "Error: Duplicate --%s flag.\n", cur_flag_ptr);
      // g_stderr_written_to = 1;
      return kPglRetInvalidCmdline;
    }
    prev_flag_ptr = cur_flag_ptr;
    prev_flag_len = cur_flag_len;
  }
  return kPglRetSuccess;
}

PglErr InitLogfile(uint32_t always_stderr, char* outname, char* outname_end) {
  snprintf(outname_end, kMaxOutfnameExtBlen, ".log");
  g_logfile = fopen(outname, "w");
  if (unlikely(!g_logfile)) {
    fflush(stdout);
    fprintf(stderr, "Error: Failed to open %s for logging: %s.\n", outname, strerror(errno));
    // g_stderr_written_to = 1;
    return kPglRetOpenFail;
  }
  fprintf(always_stderr? stderr : stdout, "Logging to %s.\n", outname);
  return kPglRetSuccess;
}

BoolErr CleanupLogfile(uint32_t print_end_time) {
  char* write_iter = strcpya_k(g_logbuf, "End time: ");
  time_t rawtime;
  time(&rawtime);
  write_iter = Stpcpy(write_iter, ctime(&rawtime));  // has trailing \n
  if (print_end_time) {
    fputs(g_logbuf, stdout);
  }
  BoolErr ret_boolerr = 0;
  if (g_logfile) {
    if (!g_log_failed) {
      logputs_silent("\n");
      logputs_silent(g_logbuf);
      if (unlikely(fclose(g_logfile))) {
        fflush(stdout);
        fprintf(stderr, "Error: Failed to finish writing to log: %s.\n", strerror(errno));
        ret_boolerr = 1;
      }
    } else {
      fclose(g_logfile);
    }
    g_logfile = nullptr;
  }
  return ret_boolerr;
}

// manually managed, very large stack
unsigned char* g_bigstack_base = nullptr;
unsigned char* g_bigstack_end = nullptr;

uintptr_t DetectMib() {
  int64_t llxx;
  // return zero if detection failed
  // see e.g. http://nadeausoftware.com/articles/2012/09/c_c_tip_how_get_physical_memory_size_system .
#ifdef __APPLE__
  int32_t mib[2];
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  llxx = 0;
  size_t sztmp = sizeof(int64_t);
  sysctl(&mib[0], 2, &llxx, &sztmp, nullptr, 0);
  llxx /= 1048576;
#else
#  ifdef _WIN32
  MEMORYSTATUSEX memstatus;
  memstatus.dwLength = sizeof(memstatus);
  GlobalMemoryStatusEx(&memstatus);
  llxx = memstatus.ullTotalPhys / 1048576;
#  else
  llxx = S_CAST(uint64_t, sysconf(_SC_PHYS_PAGES)) * S_CAST(size_t, sysconf(_SC_PAGESIZE)) / 1048576;
#  endif
#endif
  return llxx;
}

uintptr_t GetDefaultAllocMib() {
  const uintptr_t total_mib = DetectMib();
  if (!total_mib) {
    return kBigstackDefaultMib;
  }
  if (total_mib < (kBigstackMinMib * 2)) {
    return kBigstackMinMib;
  }
  return (total_mib / 2);
}

PglErr InitBigstack(uintptr_t malloc_size_mib, uintptr_t* malloc_mib_final_ptr, unsigned char** bigstack_ua_ptr) {
  // guarantee contiguous malloc space outside of main workspace
  unsigned char* bubble;

  // this is pointless with overcommit on, may want to conditionally compile
  // this out, and/or conditionally skip this
  if (unlikely(pgl_malloc(kNonBigstackMin, &bubble))) {
    return kPglRetNomem;
  }

  assert(malloc_size_mib >= kBigstackMinMib);
#ifndef __LP64__
  assert(malloc_size_mib <= 2047);
#endif
  // don't use pgl_malloc here since we don't automatically want to set
  // g_failed_alloc_attempt_size on failure
  unsigned char* bigstack_ua = S_CAST(unsigned char*, malloc(malloc_size_mib * 1048576 * sizeof(char)));
  // this is thwarted by overcommit, but still better than nothing...
  while (!bigstack_ua) {
    malloc_size_mib = (malloc_size_mib * 3) / 4;
    if (malloc_size_mib < kBigstackMinMib) {
      malloc_size_mib = kBigstackMinMib;
    }
    bigstack_ua = S_CAST(unsigned char*, malloc(malloc_size_mib * 1048576 * sizeof(char)));
    if (unlikely((!bigstack_ua) && (malloc_size_mib == kBigstackMinMib))) {
      // switch to "goto cleanup" pattern if any more exit points are needed
      g_failed_alloc_attempt_size = kBigstackMinMib * 1048576;
      free(bubble);
      return kPglRetNomem;
    }
  }
  // force 64-byte align to make cache line sensitivity work
  unsigned char* bigstack_initial_base = R_CAST(unsigned char*, RoundUpPow2(R_CAST(uintptr_t, bigstack_ua), kCacheline));
  g_bigstack_base = bigstack_initial_base;
  // last 576 bytes now reserved for g_one_char_strs + overread buffer
  g_bigstack_end = &(bigstack_initial_base[RoundDownPow2(malloc_size_mib * 1048576 - 576 - S_CAST(uintptr_t, bigstack_initial_base - bigstack_ua), kCacheline)]);
  free(bubble);
  uintptr_t* one_char_iter = R_CAST(uintptr_t*, g_bigstack_end);
#ifdef __LP64__
  // assumes little-endian
  uintptr_t cur_word = 0x3000200010000LLU;
  for (uint32_t uii = 0; uii != 64; ++uii) {
    *one_char_iter++ = cur_word;
    cur_word += 0x4000400040004LLU;
  }
#else
  uintptr_t cur_word = 0x10000;
  for (uint32_t uii = 0; uii != 128; ++uii) {
    *one_char_iter++ = cur_word;
    cur_word += 0x20002;
  }
#endif
  g_one_char_strs = R_CAST(const char*, g_bigstack_end);

  // plink2 doesn't actually need these here, but short programs using
  // plink2_common benefit from this
  g_input_missing_geno_ptr = &(g_one_char_strs[96]);
  g_output_missing_geno_ptr = &(g_one_char_strs[92]);

  *malloc_mib_final_ptr = malloc_size_mib;
  *bigstack_ua_ptr = bigstack_ua;
  return kPglRetSuccess;
}


BoolErr bigstack_calloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = S_CAST(unsigned char*, bigstack_alloc(ct));
  if (unlikely(!(*uc_arr_ptr))) {
    return 1;
  }
  memset(*uc_arr_ptr, 0, ct);
  return 0;
}

BoolErr bigstack_calloc_d(uintptr_t ct, double** d_arr_ptr) {
  *d_arr_ptr = S_CAST(double*, bigstack_alloc(ct * sizeof(double)));
  if (unlikely(!(*d_arr_ptr))) {
    return 1;
  }
  ZeroDArr(ct, *d_arr_ptr);
  return 0;
}

BoolErr bigstack_calloc_f(uintptr_t ct, float** f_arr_ptr) {
  *f_arr_ptr = S_CAST(float*, bigstack_alloc(ct * sizeof(float)));
  if (unlikely(!(*f_arr_ptr))) {
    return 1;
  }
  ZeroFArr(ct, *f_arr_ptr);
  return 0;
}

BoolErr bigstack_calloc_u16(uintptr_t ct, uint16_t** u16_arr_ptr) {
  *u16_arr_ptr = S_CAST(uint16_t*, bigstack_alloc(ct * sizeof(int16_t)));
  if (unlikely(!(*u16_arr_ptr))) {
    return 1;
  }
  memset(*u16_arr_ptr, 0, ct * sizeof(int16_t));
  return 0;
}

BoolErr bigstack_calloc_u32(uintptr_t ct, uint32_t** u32_arr_ptr) {
  *u32_arr_ptr = S_CAST(uint32_t*, bigstack_alloc(ct * sizeof(int32_t)));
  if (unlikely(!(*u32_arr_ptr))) {
    return 1;
  }
  ZeroU32Arr(ct, *u32_arr_ptr);
  return 0;
}

BoolErr bigstack_calloc_w(uintptr_t ct, uintptr_t** w_arr_ptr) {
  *w_arr_ptr = S_CAST(uintptr_t*, bigstack_alloc(ct * sizeof(intptr_t)));
  if (unlikely(!(*w_arr_ptr))) {
    return 1;
  }
  ZeroWArr(ct, *w_arr_ptr);
  return 0;
}

BoolErr bigstack_calloc_u64(uintptr_t ct, uint64_t** u64_arr_ptr) {
  *u64_arr_ptr = S_CAST(uint64_t*, bigstack_alloc(ct * sizeof(int64_t)));
  if (unlikely(!(*u64_arr_ptr))) {
    return 1;
  }
  ZeroU64Arr(ct, *u64_arr_ptr);
  return 0;
}

BoolErr bigstack_calloc_cp(uintptr_t ct, char*** cp_arr_ptr) {
  *cp_arr_ptr = S_CAST(char**, bigstack_alloc(ct * sizeof(intptr_t)));
  if (unlikely(!(*cp_arr_ptr))) {
    return 1;
  }
  ZeroPtrArr(ct, *cp_arr_ptr);
  return 0;
}

BoolErr bigstack_calloc_kcp(uintptr_t ct, const char*** kcp_arr_ptr) {
  *kcp_arr_ptr = S_CAST(const char**, bigstack_alloc(ct * sizeof(intptr_t)));
  if (unlikely(!(*kcp_arr_ptr))) {
    return 1;
  }
  ZeroPtrArr(ct, *kcp_arr_ptr);
  return 0;
}

BoolErr bigstack_end_calloc_uc(uintptr_t ct, unsigned char** uc_arr_ptr) {
  *uc_arr_ptr = S_CAST(unsigned char*, bigstack_end_alloc(ct));
  if (unlikely(!(*uc_arr_ptr))) {
    return 1;
  }
  memset(*uc_arr_ptr, 0, ct);
  return 0;
}

BoolErr bigstack_end_calloc_d(uintptr_t ct, double** d_arr_ptr) {
  *d_arr_ptr = S_CAST(double*, bigstack_end_alloc(ct * sizeof(double)));
  if (unlikely(!(*d_arr_ptr))) {
    return 1;
  }
  ZeroDArr(ct, *d_arr_ptr);
  return 0;
}

BoolErr bigstack_end_calloc_f(uintptr_t ct, float** f_arr_ptr) {
  *f_arr_ptr = S_CAST(float*, bigstack_end_alloc(ct * sizeof(float)));
  if (unlikely(!(*f_arr_ptr))) {
    return 1;
  }
  ZeroFArr(ct, *f_arr_ptr);
  return 0;
}

BoolErr bigstack_end_calloc_u32(uintptr_t ct, uint32_t** u32_arr_ptr) {
  *u32_arr_ptr = S_CAST(uint32_t*, bigstack_end_alloc(ct * sizeof(int32_t)));
  if (unlikely(!(*u32_arr_ptr))) {
    return 1;
  }
  ZeroU32Arr(ct, *u32_arr_ptr);
  return 0;
}

BoolErr bigstack_end_calloc_w(uintptr_t ct, uintptr_t** w_arr_ptr) {
  *w_arr_ptr = S_CAST(uintptr_t*, bigstack_end_alloc(ct * sizeof(intptr_t)));
  if (unlikely(!(*w_arr_ptr))) {
    return 1;
  }
  ZeroWArr(ct, *w_arr_ptr);
  return 0;
}

BoolErr bigstack_end_calloc_u64(uintptr_t ct, uint64_t** u64_arr_ptr) {
  *u64_arr_ptr = S_CAST(uint64_t*, bigstack_end_alloc(ct * sizeof(int64_t)));
  if (unlikely(!(*u64_arr_ptr))) {
    return 1;
  }
  ZeroU64Arr(ct, *u64_arr_ptr);
  return 0;
}

BoolErr bigstack_end_calloc_cp(uintptr_t ct, char*** cp_arr_ptr) {
  *cp_arr_ptr = S_CAST(char**, bigstack_end_alloc(ct * sizeof(intptr_t)));
  if (unlikely(!(*cp_arr_ptr))) {
    return 1;
  }
  ZeroPtrArr(ct, *cp_arr_ptr);
  return 0;
}


BoolErr PushLlStr(const char* str, LlStr** ll_stack_ptr) {
  uintptr_t blen = strlen(str) + 1;
  LlStr* new_llstr;
  if (unlikely(pgl_malloc(sizeof(LlStr) + blen, &new_llstr))) {
    return 1;
  }
  new_llstr->next = *ll_stack_ptr;
  memcpy(new_llstr->str, str, blen);
  *ll_stack_ptr = new_llstr;
  return 0;
}

/*
BoolErr push_llstr_counted(const char* str, uint32_t slen, LlStr** ll_stack_ptr) {
  LlStr* new_llstr;
  if (pgl_malloc(sizeof(LlStr) + slen + 1, &new_llstr)) {
    return 1;
  }
  new_llstr->next = *ll_stack_ptr;
  memcpy(new_llstr->str, str, slen);
  new_llstr->str[slen] = '\0';
  *ll_stack_ptr = new_llstr;
  return 0;
}
*/


BoolErr CountAndMeasureMultistrReverseAlloc(const char* multistr, uintptr_t max_str_ct, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr, const char*** strptr_arrp) {
  // assumes multistr is nonempty
  assert(multistr[0]);
  uintptr_t ct = 0;
  uintptr_t max_blen = *max_blen_ptr;
  const char** strptr_arr_iter = *strptr_arrp;
  do {
    if (unlikely(++ct > max_str_ct)) {
      return 1;
    }
    const uintptr_t blen = strlen(multistr) + 1;
    if (blen > max_blen) {
      max_blen = blen;
    }
    *(--strptr_arr_iter) = multistr;
    multistr = &(multistr[blen]);
  } while (*multistr);
  *str_ct_ptr = ct;
  *max_blen_ptr = max_blen;
  *strptr_arrp = strptr_arr_iter;
  return 0;
}

BoolErr MultistrToStrboxDedupArenaAlloc(unsigned char* arena_top, const char* multistr, unsigned char** arena_bottom_ptr, char** sorted_strbox_ptr, uint32_t* str_ct_ptr, uintptr_t* max_blen_ptr) {
  const char** strptr_arr = R_CAST(const char**, arena_top);
  uintptr_t max_str_blen = 0;
  uint32_t str_ct;
  if (unlikely(CountAndMeasureMultistrReverseAlloc(multistr, (arena_top - (*arena_bottom_ptr)) / sizeof(intptr_t), &str_ct, &max_str_blen, &strptr_arr))) {
    return 1;
  }
  const uintptr_t strbox_byte_ct = RoundUpPow2(str_ct * max_str_blen, kCacheline);
  if (unlikely((R_CAST(uintptr_t, strptr_arr) - R_CAST(uintptr_t, *arena_bottom_ptr)) < strbox_byte_ct)) {
    return 1;
  }
  // multistr can currently come from command line parser, which does not
  // guarantee safe overread
  StrptrArrSort(str_ct, strptr_arr);
  *sorted_strbox_ptr = R_CAST(char*, *arena_bottom_ptr);
  str_ct = CopyAndDedupSortedStrptrsToStrbox(strptr_arr, str_ct, max_str_blen, *sorted_strbox_ptr);
  *arena_bottom_ptr += RoundUpPow2(str_ct * max_str_blen, kCacheline);
  *str_ct_ptr = str_ct;
  *max_blen_ptr = max_str_blen;
  return 0;
}


void DivisionMagicNums(uint32_t divisor, uint64_t* multp, uint32_t* __restrict pre_shiftp, uint32_t* __restrict post_shiftp, uint32_t* __restrict incrp) {
  // Enables fast integer division by a constant not known until runtime.  See
  // http://ridiculousfish.com/blog/posts/labor-of-division-episode-iii.html .
  // Assumes divisor is not zero, of course.
  // May want to populate a struct instead.
  assert(divisor);
  if (!(divisor & (divisor - 1))) {
    // power of 2
    *multp = 1;
    *pre_shiftp = 0;
    *post_shiftp = ctzu32(divisor);
    *incrp = 0;
    return;
  }
  uint32_t quotient = 0x80000000U / divisor;
  uint32_t remainder = 0x80000000U - (quotient * divisor);
  const uint32_t ceil_log_2_d = 1 + bsru32(divisor);
  uint32_t down_multiplier = 0;
  uint32_t down_exponent = 0;
  uint32_t has_magic_down = 0;
  for (uint32_t exponent = 0; ; ++exponent) {
    if (remainder >= divisor - remainder) {
      quotient = quotient * 2 + 1;
      remainder = remainder * 2 - divisor;
    } else {
      quotient = quotient * 2;
      remainder = remainder * 2;
    }
    if ((exponent >= ceil_log_2_d) || (divisor - remainder) <= (1U << exponent)) {
      if (exponent < ceil_log_2_d) {
        *multp = quotient + 1;
        *pre_shiftp = 0;
        *post_shiftp = 32 + exponent;
        *incrp = 0;
        return;
      }
      break;
    }
    if ((!has_magic_down) && (remainder <= (1U << exponent))) {
      has_magic_down = 1;
      down_multiplier = quotient;
      down_exponent = exponent;
    }
  }
  if (divisor & 1) {
    *multp = down_multiplier;
    *pre_shiftp = 0;
    *post_shiftp = 32 + down_exponent;
    *incrp = 1;
    return;
  }
  *pre_shiftp = ctzu32(divisor);
  uint32_t dummy;
  DivisionMagicNums(divisor >> (*pre_shiftp), multp, &dummy, post_shiftp, incrp);
}


void FillBitsNz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr) {
  assert(end_idx > start_idx);
  uintptr_t maj_start = start_idx / kBitsPerWord;
  uintptr_t maj_end = end_idx / kBitsPerWord;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bitarr[maj_start] |= (k1LU << (end_idx % kBitsPerWord)) - (k1LU << (start_idx % kBitsPerWord));
  } else {
    bitarr[maj_start] |= ~((k1LU << (start_idx % kBitsPerWord)) - k1LU);
    SetAllWArr(maj_end - maj_start - 1, &(bitarr[maj_start + 1]));
    minor = end_idx % kBitsPerWord;
    if (minor) {
      bitarr[maj_end] |= (k1LU << minor) - k1LU;
    }
  }
}

void ClearBitsNz(uintptr_t start_idx, uintptr_t end_idx, uintptr_t* bitarr) {
  assert(end_idx > start_idx);
  uintptr_t maj_start = start_idx / kBitsPerWord;
  uintptr_t maj_end = end_idx / kBitsPerWord;
  uintptr_t minor;
  if (maj_start == maj_end) {
    bitarr[maj_start] &= ~((k1LU << (end_idx % kBitsPerWord)) - (k1LU << (start_idx % kBitsPerWord)));
  } else {
    bitarr[maj_start] = bzhi(bitarr[maj_start], start_idx % kBitsPerWord);
    ZeroWArr(maj_end - maj_start - 1, &(bitarr[maj_start + 1]));
    minor = end_idx % kBitsPerWord;
    if (minor) {
      bitarr[maj_end] &= ~((k1LU << minor) - k1LU);
    }
  }
}

// floor permitted to be -1, though not smaller than that.
int32_t FindLast1BitBeforeBounded(const uintptr_t* bitarr, uint32_t loc, int32_t floor) {
  const uintptr_t* bitarr_ptr = &(bitarr[loc / kBitsPerWord]);
  uint32_t remainder = loc % kBitsPerWord;
  uintptr_t ulii;
  if (remainder) {
    ulii = bzhi(*bitarr_ptr, remainder);
    if (ulii) {
      const uint32_t set_bit_loc = loc - remainder + bsrw(ulii);
      return MAXV(S_CAST(int32_t, set_bit_loc), floor);
    }
  }
  const uintptr_t* bitarr_last = &(bitarr[S_CAST(uint32_t, floor + 1) / kBitsPerWord]);
  do {
    if (bitarr_ptr <= bitarr_last) {
      return floor;
    }
    ulii = *(--bitarr_ptr);
  } while (!ulii);
  const uint32_t set_bit_loc = S_CAST(uintptr_t, bitarr_ptr - bitarr) * kBitsPerWord + bsrw(ulii);
  return MAXV(S_CAST(int32_t, set_bit_loc), floor);
}

void BitvecAndCopy(const uintptr_t* __restrict source1_bitvec, const uintptr_t* __restrict source2_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec) {
#ifdef __LP64__
  VecW* target_bitvvec = R_CAST(VecW*, target_bitvec);
  const VecW* source1_bitvvec = R_CAST(const VecW*, source1_bitvec);
  const VecW* source2_bitvvec = R_CAST(const VecW*, source2_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  for (uintptr_t ulii = 0; ulii != full_vec_ct; ++ulii) {
    target_bitvvec[ulii] = source1_bitvvec[ulii] & source2_bitvvec[ulii];
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    target_bitvec[base_idx] = source1_bitvec[base_idx] & source2_bitvec[base_idx];
    target_bitvec[base_idx + 1] = source1_bitvec[base_idx + 1] & source2_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    target_bitvec[word_ct - 1] = source1_bitvec[word_ct - 1] & source2_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    target_bitvec[widx] = source1_bitvec[widx] & source2_bitvec[widx];
  }
#endif
}

void BitvecInvmaskCopy(const uintptr_t* __restrict source_bitvec, const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* target_bitvec) {
  // target_bitvec := source_bitvec AND (~exclude_bitvec)
#ifdef __LP64__
  VecW* target_bitvvec = R_CAST(VecW*, target_bitvec);
  const VecW* source_bitvvec = R_CAST(const VecW*, source_bitvec);
  const VecW* exclude_bitvvec = R_CAST(const VecW*, exclude_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  for (uintptr_t ulii = 0; ulii != full_vec_ct; ++ulii) {
    target_bitvvec[ulii] = vecw_and_notfirst(exclude_bitvvec[ulii], source_bitvvec[ulii]);
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    target_bitvec[base_idx] = (~exclude_bitvec[base_idx]) & source_bitvec[base_idx];
    target_bitvec[base_idx + 1] = (~exclude_bitvec[base_idx + 1]) & source_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    target_bitvec[word_ct - 1] = source_bitvec[word_ct - 1] & (~exclude_bitvec[word_ct - 1]);
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    target_bitvec[widx] = source_bitvec[widx] & (~exclude_bitvec[widx]);
  }
#endif
}

void BitvecInvertCopy(const uintptr_t* __restrict source_bitvec, uintptr_t word_ct, uintptr_t* __restrict target_bitvec) {
#ifdef __LP64__
  const VecW* source_bitvvec_iter = R_CAST(const VecW*, source_bitvec);
  VecW* target_bitvvec_iter = R_CAST(VecW*, target_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  const VecW all1 = VCONST_W(~k0LU);
  // As of Apple clang 11, this manual unroll is no longer relevant.  todo:
  // check Linux performance, and remove all of these unrolls if perf is good
  // enough without them.
  if (full_vec_ct & 1) {
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
  }
  if (full_vec_ct & 2) {
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
    *target_bitvvec_iter++ = (*source_bitvvec_iter++) ^ all1;
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    target_bitvec[base_idx] = ~source_bitvec[base_idx];
    target_bitvec[base_idx + 1] = ~source_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    target_bitvec[word_ct - 1] = ~source_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    target_bitvec[widx] = ~source_bitvec[widx];
  }
#endif
}

/*
void BitvecXor(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec) {
  // main_bitvec := main_bitvec XOR arg_bitvec
#ifdef __LP64__
  VecW* main_bitvvec_iter = R_CAST(VecW*, main_bitvec);
  const VecW* arg_bitvvec_iter = R_CAST(const VecW*, arg_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    *main_bitvvec_iter++ ^= (*arg_bitvvec_iter++);
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter++ ^= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ ^= (*arg_bitvvec_iter++);
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter++ ^= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ ^= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ ^= (*arg_bitvvec_iter++);
    *main_bitvvec_iter++ ^= (*arg_bitvvec_iter++);
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] ^= arg_bitvec[base_idx];
    main_bitvec[base_idx + 1] ^= arg_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] ^= arg_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    main_bitvec[widx] ^= arg_bitvec[widx];
  }
#endif
}
*/

void BitvecInvertAndMask(const uintptr_t* __restrict include_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec) {
  // main_bitvec := (~main_bitvec) AND include_bitvec
  // this corresponds _mm_andnot() operand order
#ifdef __LP64__
  VecW* main_bitvvec_iter = R_CAST(VecW*, main_bitvec);
  const VecW* include_bitvvec_iter = R_CAST(const VecW*, include_bitvec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    *main_bitvvec_iter = vecw_and_notfirst(*main_bitvvec_iter, *include_bitvvec_iter++);
    ++main_bitvvec_iter;
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter = vecw_and_notfirst(*main_bitvvec_iter, *include_bitvvec_iter++);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = vecw_and_notfirst(*main_bitvvec_iter, *include_bitvvec_iter++);
    ++main_bitvvec_iter;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter = vecw_and_notfirst(*main_bitvvec_iter, *include_bitvvec_iter++);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = vecw_and_notfirst(*main_bitvvec_iter, *include_bitvvec_iter++);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = vecw_and_notfirst(*main_bitvvec_iter, *include_bitvvec_iter++);
    ++main_bitvvec_iter;
    *main_bitvvec_iter = vecw_and_notfirst(*main_bitvvec_iter, *include_bitvvec_iter++);
    ++main_bitvvec_iter;
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] = (~main_bitvec[base_idx]) & include_bitvec[base_idx];
    main_bitvec[base_idx + 1] = (~main_bitvec[base_idx + 1]) & include_bitvec[base_idx + 1];
  }
#  endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] = (~main_bitvec[word_ct - 1]) & include_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    main_bitvec[widx] = (~main_bitvec[widx]) & include_bitvec[widx];
  }
#endif
}

/*
void BitvecOrNot(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec) {
  // main_bitvec := main_bitvec OR (~arg_bitvec)
#ifdef __LP64__
  VecW* main_bitvvec_iter = (VecW*)main_bitvec;
  const VecW* arg_bitvvec_iter = (const VecW*)arg_bitvec;
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    // todo: verify the compiler isn't dumb here
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
  }
  if (full_vec_ct & 2) {
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
    *main_bitvvec_iter++ |= ~(*arg_bitvvec_iter++);
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    main_bitvec[base_idx] |= ~arg_bitvec[base_idx];
    main_bitvec[base_idx + 1] |= ~arg_bitvec[base_idx + 1]
  }
#  endif
  if (word_ct & 1) {
    main_bitvec[word_ct - 1] |= ~arg_bitvec[word_ct - 1];
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    main_bitvec[widx] |= ~arg_bitvec[widx];
  }
#endif
}
*/


int32_t GetVariantUidxWithoutHtable(const char* idstr, const char* const* variant_ids, const uintptr_t* variant_include, uint32_t variant_ct) {
  const uint32_t id_blen = strlen(idstr) + 1;
  uintptr_t widx = ~k0LU;
  int32_t retval = -1;
  for (uint32_t variant_idx = 0; variant_idx != variant_ct; ) {
    uintptr_t variant_include_bits;
    do {
      variant_include_bits = variant_include[++widx];
    } while (!variant_include_bits);
    const uintptr_t variant_uidx_base = widx * kBitsPerWord;
    const char* const* cur_variant_ids = &(variant_ids[variant_uidx_base]);
    variant_idx += PopcountWord(variant_include_bits);
    do {
      const uint32_t uidx_lowbits = ctzw(variant_include_bits);
      if (memequal(idstr, cur_variant_ids[uidx_lowbits], id_blen)) {
        if (retval != -1) {
          // duplicate
          return -2;
        }
        retval = S_CAST(int32_t, uidx_lowbits + variant_uidx_base);
      }
      variant_include_bits &= variant_include_bits - 1;
    } while (variant_include_bits);
  }
  return retval;
}


// MurmurHash3, from
// https://code.google.com/p/smhasher/source/browse/trunk/MurmurHash3.cpp
CSINLINE uint32_t rotl32(uint32_t x, int8_t r) {
  return (x << r) | (x >> (32 - r));
}

static inline uint32_t getblock32(const uint32_t* p, int i) {
  return p[i];
}

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche

CSINLINE2 uint32_t fmix32(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

uint32_t Hash32(const void* key, uint32_t len) {
  const uint8_t* data = S_CAST(const uint8_t*, key);
  const int32_t nblocks = len / 4;

  uint32_t h1 = 0;
  // uint32_t h1 = seed;

  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;

  //----------
  // body

  const uint32_t* blocks = R_CAST(const uint32_t*, data + nblocks*4);

  int32_t i;
  uint32_t k1;
  for(i = -nblocks; i; i++) {
      k1 = getblock32(blocks,i);

      k1 *= c1;
      k1 = rotl32(k1,15);
      k1 *= c2;

      h1 ^= k1;
      h1 = rotl32(h1,13);
      h1 = h1*5+0xe6546b64;
  }

  //----------
  // tail

  const uint8_t* tail = R_CAST(const uint8_t*, data + nblocks*4);

  k1 = 0;

  switch(len & 3) {
    case 3:
      k1 ^= tail[2] << 16;
      // fall through
    case 2:
      k1 ^= tail[1] << 8;
      // fall through
    case 1:
      k1 ^= tail[0];
      k1 *= c1;
      k1 = rotl32(k1,15);
      k1 *= c2;
      h1 ^= k1;
  }

  //----------
  // finalization

  h1 ^= len;

  return fmix32(h1);
}


/*
uint32_t is_composite6(uintptr_t num) {
  // assumes num is congruent to 1 or 5 mod 6.
  // can speed this up by ~50% by hardcoding avoidance of multiples of 5/7,
  // but this isn't currently a bottleneck so I'll keep this simple
  assert((num % 6 == 1) || (num % 6 == 5));
  uintptr_t divisor = 5;
  while (divisor * divisor <= num) {
    if (!(num % divisor)) {
      return 1;
    }
    divisor += 2;
    if (!(num % divisor)) {
      return 1;
    }
    divisor += 4;
  }
  return 0;
}

uintptr_t geqprime(uintptr_t floor) {
  // assumes floor is odd and greater than 1.  Returns 5 if floor = 3,
  // otherwise returns the first prime >= floor.
  assert((floor % 2) && (floor > 1));
  uintptr_t ulii = floor % 3;
  if (!ulii) {
    floor += 2;
  } else if (ulii == 1) {
    goto geqprime_1mod6;
  }
  while (is_composite6(floor)) {
    floor += 2;
  geqprime_1mod6:
    if (!is_composite6(floor)) {
      return floor;
    }
    floor += 4;
  }
  return floor;
}

uintptr_t leqprime(uintptr_t ceil) {
  // assumes ceil is odd and greater than 4.  Returns the first prime <= ceil.
  assert((ceil % 2) && (ceil > 4));
  uintptr_t ulii = ceil % 3;
  if (!ulii) {
    ceil -= 2;
  } else if (ulii == 2) {
    goto leqprime_5mod6;
  }
  while (is_composite6(ceil)) {
    ceil -= 2;
  leqprime_5mod6:
    if (!is_composite6(ceil)) {
      return ceil;
    }
    ceil -= 4;
  }
  return ceil;
}
*/

BoolErr HtableGoodSizeAlloc(uint32_t item_ct, uintptr_t bytes_avail, uint32_t** htable_ptr, uint32_t* htable_size_ptr) {
  bytes_avail &= (~(kCacheline - k1LU));
  uint32_t htable_size = GetHtableFastSize(item_ct);
  if (htable_size > bytes_avail / sizeof(int32_t)) {
    if (unlikely(!bytes_avail)) {
      return 1;
    }
    htable_size = bytes_avail / sizeof(int32_t);
    // htable_size = leqprime((bytes_avail / sizeof(int32_t)) - 1);
    if (unlikely(htable_size < item_ct * 2)) {
      return 1;
    }
  }
  *htable_ptr = S_CAST(uint32_t*, bigstack_alloc_raw_rd(htable_size * sizeof(int32_t)));
  *htable_size_ptr = htable_size;
  return 0;
}

uint32_t PopulateStrboxHtable(const char* strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t str_htable_size, uint32_t* str_htable) {
  // may want subset_mask parameter later
  SetAllU32Arr(str_htable_size, str_htable);
  const char* strbox_iter = strbox;
  for (uintptr_t str_idx = 0; str_idx != str_ct; ++str_idx) {
    const uint32_t slen = strlen(strbox_iter);
    // previously used quadratic probing, but turns out that that isn't
    // meaningfully better than linear probing.
    for (uint32_t hashval = Hashceil(strbox_iter, slen, str_htable_size); ; ) {
      const uint32_t cur_htable_entry = str_htable[hashval];
      if (cur_htable_entry == UINT32_MAX) {
        str_htable[hashval] = str_idx;
        break;
      }
      if (memequal(strbox_iter, &(strbox[cur_htable_entry * max_str_blen]), slen + 1)) {
        // guaranteed to be positive
        return str_idx;
      }
      if (++hashval == str_htable_size) {
        hashval = 0;
      }
    }
    strbox_iter = &(strbox_iter[max_str_blen]);
  }
  return 0;
}

// could merge this with non-subset case, but this isn't much code
/*
uint32_t populate_strbox_subset_htable(const uintptr_t* __restrict subset_mask, const char* strbox, uintptr_t raw_str_ct, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t str_htable_size, uint32_t* str_htable) {
  // may want subset_mask parameter later
  SetAllU32Arr(str_htable_size, str_htable);
  uintptr_t str_uidx = 0;
  for (uintptr_t str_idx = 0; str_idx != str_ct; ++str_idx, ++str_uidx) {
    MovWTo1Bit(subset_mask, &str_uidx);
    const char* cur_str = &(strbox[str_uidx * max_str_blen]);
    const uint32_t slen = strlen(cur_str);
    uint32_t hashval = Hashceil(cur_str, slen, str_htable_size);
    while (1) {
      const uint32_t cur_htable_entry = str_htable[hashval];
      if (cur_htable_entry == UINT32_MAX) {
        str_htable[hashval] = str_uidx;
        break;
      }
      if (memequal(cur_str, &(strbox[cur_htable_entry * max_str_blen]), slen + 1)) {
        // guaranteed to be positive
        return str_uidx;
      }
      if (++hashval == str_htable_size) {
        hashval = 0;
      }
    }
  }
  return 0;
}
*/

uint32_t IdHtableFind(const char* cur_id, const char* const* item_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size) {
  // returns UINT32_MAX on failure
  for (uint32_t hashval = Hashceil(cur_id, cur_id_slen, id_htable_size); ; ) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if ((cur_htable_idval == UINT32_MAX) || (!strcmp(cur_id, item_ids[cur_htable_idval]))) {
      return cur_htable_idval;
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
    if ((cur_htable_idval == UINT32_MAX) || (memequal(cur_id, item_ids[cur_htable_idval], cur_id_slen) && (!item_ids[cur_htable_idval][cur_id_slen]))) {
      return cur_htable_idval;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

// assumes cur_id_slen < max_str_blen.
// requires cur_id to be null-terminated.
uint32_t StrboxHtableFind(const char* cur_id, const char* strbox, const uint32_t* id_htable, uintptr_t max_str_blen, uint32_t cur_id_slen, uint32_t id_htable_size) {
  const uint32_t cur_id_blen = cur_id_slen + 1;
  for (uint32_t hashval = Hashceil(cur_id, cur_id_slen, id_htable_size); ; ) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if ((cur_htable_idval == UINT32_MAX) || memequal(cur_id, &(strbox[cur_htable_idval * max_str_blen]), cur_id_blen)) {
      return cur_htable_idval;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

uint32_t VariantIdDupflagHtableFind(const char* idbuf, const char* const* variant_ids, const uint32_t* id_htable, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen) {
  // assumes duplicate variant IDs are flagged, but full variant_uidx linked
  // lists are not stored
  // idbuf does not need to be null-terminated (note that this is currently
  // achieved in a way that forces variant_ids[] entries to not be too close
  // to the end of bigstack, otherwise memequal behavior is potentially
  // undefined)
  // returns UINT32_MAX on failure, value with bit 31 set on duplicate
  if (cur_id_slen > max_id_slen) {
    return UINT32_MAX;
  }
  for (uint32_t hashval = Hashceil(idbuf, cur_id_slen, id_htable_size); ; ) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    if ((cur_htable_idval == UINT32_MAX) || (memequal(idbuf, variant_ids[cur_htable_idval & 0x7fffffff], cur_id_slen) && (!variant_ids[cur_htable_idval & 0x7fffffff][cur_id_slen]))) {
      return cur_htable_idval;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}

uint32_t VariantIdDupHtableFind(const char* idbuf, const char* const* variant_ids, const uint32_t* id_htable, const uint32_t* htable_dup_base, uint32_t cur_id_slen, uint32_t id_htable_size, uint32_t max_id_slen, uint32_t* llidx_ptr) {
  // Permits duplicate entries.  Similar to plink 1.9
  // extract_exclude_process_token().
  // - Returns UINT32_MAX on failure (llidx currently unset in that case),
  //   otherwise returns the index of the first match (which will have the
  //   highest index, due to how the linked list is constructed)
  // - Sets llidx to UINT32_MAX if not a duplicate, otherwise it's the
  //   position in htable_dup_base[] of the next {variant_uidx, next_llidx}
  //   linked list entry.
  // - idbuf does not need to be null-terminated.
  if (cur_id_slen > max_id_slen) {
    return UINT32_MAX;
  }
  for (uint32_t hashval = Hashceil(idbuf, cur_id_slen, id_htable_size); ; ) {
    const uint32_t cur_htable_idval = id_htable[hashval];
    const uint32_t cur_dup = cur_htable_idval >> 31;
    uint32_t cur_llidx;
    uint32_t variant_uidx;
    if (cur_dup) {
      // UINT32_MAX empty-entry code has high bit set, so only need to check
      // here
      if (cur_htable_idval == UINT32_MAX) {
        return UINT32_MAX;
      }
      cur_llidx = cur_htable_idval << 1;
      variant_uidx = htable_dup_base[cur_llidx];
    } else {
      cur_llidx = UINT32_MAX;
      variant_uidx = cur_htable_idval;
    }
    const char* sptr = variant_ids[variant_uidx];
    if (memequal(idbuf, sptr, cur_id_slen) && (!sptr[cur_id_slen])) {
      *llidx_ptr = cur_llidx;
      return variant_uidx;
    }
    if (++hashval == id_htable_size) {
      hashval = 0;
    }
  }
}


PglErr CopySortStrboxSubsetNoalloc(const uintptr_t* __restrict subset_mask, const char* __restrict orig_strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t allow_dups, uint32_t collapse_idxs, uint32_t use_nsort, char* __restrict sorted_strbox, uint32_t* __restrict id_map) {
  // Stores a lexicographically sorted list of IDs in sorted_strbox and the raw
  // positions of the corresponding markers/samples in *id_map_ptr.  Does not
  // include excluded markers/samples in the list.
  // Assumes sorted_strbox and id_map have been allocated; use the
  // CopySortStrboxSubset() wrapper if they haven't been.
  // Note that this DOES still perform a "stack" allocation; 'Noalloc' just
  // means that sorted_strbox/id_map must already be allocated.
  if (!str_ct) {
    return kPglRetSuccess;
  }
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
#ifdef __cplusplus
    if (max_str_blen <= 60) {
      uintptr_t wkspace_entry_blen = (max_str_blen > 28)? sizeof(Strbuf60Ui) : sizeof(Strbuf28Ui);
      char* sort_wkspace;
      if (unlikely(bigstack_alloc_c(str_ct * wkspace_entry_blen, &sort_wkspace))) {
        goto CopySortStrboxSubsetNoalloc_ret_NOMEM;
      }
      const uint32_t wkspace_entry_blen_m4 = wkspace_entry_blen - 4;
      char* sort_wkspace_iter = sort_wkspace;
      uintptr_t str_uidx_base = 0;
      uintptr_t subset_mask_bits = subset_mask[0];
      for (uint32_t str_idx = 0; str_idx != str_ct; ++str_idx) {
        const uint32_t str_uidx = BitIter1(subset_mask, &str_uidx_base, &subset_mask_bits);
        strcpy(sort_wkspace_iter, &(orig_strbox[str_uidx * max_str_blen]));
        sort_wkspace_iter = &(sort_wkspace_iter[wkspace_entry_blen_m4]);
        if (collapse_idxs) {
          memcpy(sort_wkspace_iter, &str_idx, sizeof(int32_t));
        } else {
          memcpy(sort_wkspace_iter, &str_uidx, sizeof(int32_t));
        }
        sort_wkspace_iter = &(sort_wkspace_iter[sizeof(int32_t)]);
      }
      if (wkspace_entry_blen == 32) {
        SortStrbox32bFinish(str_ct, max_str_blen, use_nsort, R_CAST(Strbuf28Ui*, sort_wkspace), sorted_strbox, id_map);
      } else {
        SortStrbox64bFinish(str_ct, max_str_blen, use_nsort, R_CAST(Strbuf60Ui*, sort_wkspace), sorted_strbox, id_map);
      }
    } else {
#endif
      StrSortIndexedDerefOverread* sort_wkspace;
      if (unlikely(BIGSTACK_ALLOC_X(StrSortIndexedDerefOverread, str_ct, &sort_wkspace))) {
        goto CopySortStrboxSubsetNoalloc_ret_NOMEM;
      }
      uintptr_t str_uidx_base = 0;
      uintptr_t cur_bits = subset_mask[0];
      for (uint32_t str_idx = 0; str_idx != str_ct; ++str_idx) {
        uintptr_t str_uidx = BitIter1(subset_mask, &str_uidx_base, &cur_bits);
        sort_wkspace[str_idx].strptr = &(orig_strbox[str_uidx * max_str_blen]);
        if (collapse_idxs) {
          sort_wkspace[str_idx].orig_idx = str_idx;
        } else {
          sort_wkspace[str_idx].orig_idx = str_uidx;
        }
      }
      if (!use_nsort) {
        STD_SORT_PAR_UNSEQ(str_ct, strcmp_overread_deref, sort_wkspace);
      } else {
#ifdef __cplusplus
        StrNsortIndexedDeref* wkspace_alias = R_CAST(StrNsortIndexedDeref*, sort_wkspace);
        STD_SORT_PAR_UNSEQ(str_ct, nullptr, wkspace_alias);
#else
        STD_SORT_PAR_UNSEQ(str_ct, strcmp_natural_deref, sort_wkspace);
#endif
      }
      for (uintptr_t str_idx = 0; str_idx != str_ct; ++str_idx) {
        strcpy(&(sorted_strbox[str_idx * max_str_blen]), sort_wkspace[str_idx].strptr);
        id_map[str_idx] = sort_wkspace[str_idx].orig_idx;
      }
#ifdef __cplusplus
    }
#endif
    if (!allow_dups) {
      char* dup_id = ScanForDuplicateIds(sorted_strbox, str_ct, max_str_blen);
      if (unlikely(dup_id)) {
        TabsToSpaces(dup_id);
        logerrprintfww("Error: Duplicate ID '%s'.\n", dup_id);
        goto CopySortStrboxSubsetNoalloc_ret_MALFORMED_INPUT;
      }
    }
  }
  while (0) {
  CopySortStrboxSubsetNoalloc_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CopySortStrboxSubsetNoalloc_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr CopySortStrboxSubset(const uintptr_t* __restrict subset_mask, const char* __restrict orig_strbox, uintptr_t str_ct, uintptr_t max_str_blen, uint32_t allow_dups, uint32_t collapse_idxs, uint32_t use_nsort, char** sorted_strbox_ptr, uint32_t** id_map_ptr) {
  // id_map on bottom because --indiv-sort frees *sorted_strbox_ptr
  if (unlikely(
          bigstack_alloc_u32(str_ct, id_map_ptr) ||
          bigstack_alloc_c(str_ct * max_str_blen, sorted_strbox_ptr))) {
    return kPglRetNomem;
  }
  return CopySortStrboxSubsetNoalloc(subset_mask, orig_strbox, str_ct, max_str_blen, allow_dups, collapse_idxs, use_nsort, *sorted_strbox_ptr, *id_map_ptr);
}


void InitRangeList(RangeList* range_list_ptr) {
  range_list_ptr->names = nullptr;
  range_list_ptr->starts_range = nullptr;
  range_list_ptr->name_ct = 0;
  range_list_ptr->name_max_blen = 0;
}

void CleanupRangeList(RangeList* range_list_ptr) {
  free_cond(range_list_ptr->names);
  // starts_range now uses same allocation
}

BoolErr NumericRangeListToBitarr(const RangeList* range_list_ptr, uint32_t bitarr_size, uint32_t offset, uint32_t ignore_overflow, uintptr_t* bitarr) {
  // bitarr assumed to be initialized (but not necessarily zero-initialized)
  const char* names = range_list_ptr->names;
  const unsigned char* starts_range = range_list_ptr->starts_range;
  const uint32_t name_ct = range_list_ptr->name_ct;
  const uint32_t name_max_blen = range_list_ptr->name_max_blen;
  const uint32_t idx_max = bitarr_size + offset;
  for (uint32_t name_idx = 0; name_idx != name_ct; ++name_idx) {
    uint32_t idx1;
    if (ScanUintCapped(&(names[name_idx * name_max_blen]), idx_max, &idx1)) {
      if (likely(ignore_overflow)) {
        continue;
      }
      return 1;
    }
    if (unlikely(idx1 < offset)) {
      return 1;
    }
    if (starts_range[name_idx]) {
      ++name_idx;
      uint32_t idx2;
      if (ScanUintCapped(&(names[name_idx * name_max_blen]), idx_max, &idx2)) {
        if (unlikely(!ignore_overflow)) {
          return 1;
        }
        idx2 = idx_max - 1;
      }
      FillBitsNz(idx1 - offset, (idx2 - offset) + 1, bitarr);
    } else {
      SetBit(idx1 - offset, bitarr);
    }
  }
  return 0;
}

PglErr StringRangeListToBitarr(const char* header_line, const RangeList* range_list_ptr, const char* __restrict sorted_ids, const uint32_t* __restrict id_map, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t* bitarr, int32_t* __restrict seen_idxs) {
  // bitarr assumed to be zero-initialized
  // if fixed_len is zero, header_line is assumed to be a list of
  // space-delimited unequal-length names
  assert(token_ct);
  assert(!PopcountWords(bitarr, BitCtToWordCt(token_ct)));
  PglErr reterr = kPglRetSuccess;
  {
    const char* header_line_iter = header_line;
    const uintptr_t name_ct = range_list_ptr->name_ct;
    const uintptr_t max_id_blen = range_list_ptr->name_max_blen;
    for (uint32_t item_idx = 0; ; ) {
      const char* token_end = CommaOrTspaceTokenEnd(header_line_iter, comma_delim);
      const int32_t sorted_pos = bsearch_str(header_line_iter, sorted_ids, token_end - header_line_iter, max_id_blen, name_ct);
      if (sorted_pos != -1) {
        uint32_t cmdline_pos = sorted_pos;
        if (id_map) {
          cmdline_pos = id_map[S_CAST(uint32_t, sorted_pos)];
        }
        if (unlikely(seen_idxs[cmdline_pos] != -1)) {
          logerrprintfww("Error: --%s: '%s' appears multiple times in header line of %s.\n", range_list_flag, &(sorted_ids[S_CAST(uint32_t, sorted_pos) * max_id_blen]), file_descrip);
          goto StringRangeListToBitarr_ret_MALFORMED_INPUT;
        }
        seen_idxs[cmdline_pos] = item_idx;
        if (cmdline_pos && range_list_ptr->starts_range[cmdline_pos - 1]) {
          if (unlikely(seen_idxs[cmdline_pos - 1] == -1)) {
            uint32_t first_sorted_pos = cmdline_pos - 1;
            if (id_map) {
              for (first_sorted_pos = 0; ; ++first_sorted_pos) {
                if (id_map[first_sorted_pos] == cmdline_pos - 1) {
                  break;
                }
              }
            }
            snprintf(g_logbuf, kLogbufSize, "Error: --%s: '%s' appears before '%s' in header line of %s.\n", range_list_flag, &(sorted_ids[S_CAST(uint32_t, sorted_pos) * max_id_blen]), &(sorted_ids[first_sorted_pos * max_id_blen]), file_descrip);
            goto StringRangeListToBitarr_ret_INVALID_CMDLINE_WW;
          }
          FillBitsNz(seen_idxs[cmdline_pos - 1], item_idx + 1, bitarr);
        } else if (!(range_list_ptr->starts_range[cmdline_pos])) {
          SetBit(item_idx, bitarr);
        }
      }
      if (++item_idx == token_ct) {
        break;
      }
      if (fixed_len) {
        header_line_iter = &(header_line_iter[fixed_len]);
      } else {
        header_line_iter = FirstNonTspace(&(token_end[1]));
      }
    }
    for (uint32_t cmdline_pos = 0; cmdline_pos != name_ct; ++cmdline_pos) {
      if (unlikely(seen_idxs[cmdline_pos] == -1)) {
        uint32_t sorted_pos = cmdline_pos;
        if (id_map) {
          for (sorted_pos = 0; ; ++sorted_pos) {
            if (id_map[sorted_pos] == cmdline_pos) {
              break;
            }
          }
        }
        snprintf(g_logbuf, kLogbufSize, "Error: --%s: '%s' does not appear in header line of %s.\n", range_list_flag, &(sorted_ids[sorted_pos * max_id_blen]), file_descrip);
        goto StringRangeListToBitarr_ret_INVALID_CMDLINE_WW;
      }
    }
  }
  while (0) {
  StringRangeListToBitarr_ret_INVALID_CMDLINE_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInvalidCmdline;
    break;
  StringRangeListToBitarr_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  return reterr;
}

PglErr StringRangeListToBitarrAlloc(const char* header_line, const RangeList* range_list_ptr, const char* __restrict range_list_flag, const char* __restrict file_descrip, uint32_t token_ct, uint32_t fixed_len, uint32_t comma_delim, uintptr_t** bitarr_ptr) {
  // wrapper for StringRangeListToBitarr which allocates the bitfield and
  // temporary buffers on the heap
  uintptr_t token_ctl = BitCtToWordCt(token_ct);
  uintptr_t name_ct = range_list_ptr->name_ct;
  int32_t* seen_idxs;
  char* sorted_ids;
  uint32_t* id_map;
  if (unlikely(
          bigstack_calloc_w(token_ctl, bitarr_ptr) ||
          bigstack_alloc_i32(name_ct, &seen_idxs))) {
    return kPglRetNomem;
  }
  // kludge to use CopySortStrboxSubset()
  SetAllBits(name_ct, R_CAST(uintptr_t*, seen_idxs));
  if (unlikely(CopySortStrboxSubset(R_CAST(uintptr_t*, seen_idxs), range_list_ptr->names, name_ct, range_list_ptr->name_max_blen, 0, 0, 0, &sorted_ids, &id_map))) {
    return kPglRetNomem;
  }
  SetAllI32Arr(name_ct, seen_idxs);
  PglErr reterr = StringRangeListToBitarr(header_line, range_list_ptr, sorted_ids, id_map, range_list_flag, file_descrip, token_ct, fixed_len, comma_delim, *bitarr_ptr, seen_idxs);
  BigstackReset(seen_idxs);
  return reterr;
}


uintptr_t PopcountBitRange(const uintptr_t* bitvec, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t start_idxl = start_idx / kBitsPerWord;
  const uintptr_t start_idxlr = start_idx & (kBitsPerWord - 1);
  const uintptr_t end_idxl = end_idx / kBitsPerWord;
  const uintptr_t end_idxlr = end_idx & (kBitsPerWord - 1);
  uintptr_t ct = 0;
  if (start_idxl == end_idxl) {
    return PopcountWord(bitvec[start_idxl] & ((k1LU << end_idxlr) - (k1LU << start_idxlr)));
  }
  if (start_idxlr) {
    ct = PopcountWord(bitvec[start_idxl++] >> start_idxlr);
  }
  if (end_idxl > start_idxl) {
    ct += PopcountWordsNzbase(bitvec, start_idxl, end_idxl);
  }
  if (end_idxlr) {
    ct += PopcountWord(bzhi(bitvec[end_idxl], end_idxlr));
  }
  return ct;
}

#ifdef USE_SSE42
void PopcountWordsIntersect3val(const uintptr_t* __restrict bitvec1, const uintptr_t* __restrict bitvec2, uint32_t word_ct, uint32_t* __restrict popcount1_ptr, uint32_t* __restrict popcount2_ptr, uint32_t* __restrict popcount_intersect_ptr) {
  uint32_t ct1 = 0;
  uint32_t ct2 = 0;
  uint32_t ct3 = 0;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t word1 = bitvec1[widx];
    const uintptr_t word2 = bitvec2[widx];
    ct1 += PopcountWord(word1);
    ct2 += PopcountWord(word2);
    ct3 += PopcountWord(word1 & word2);
  }
  *popcount1_ptr = ct1;
  *popcount2_ptr = ct2;
  *popcount_intersect_ptr = ct3;
}
#else
static inline void PopcountVecsNoSse42Intersect3val(const VecW* __restrict vvec1_iter, const VecW* __restrict vvec2_iter, uint32_t vec_ct, uint32_t* __restrict popcount1_ptr, uint32_t* popcount2_ptr, uint32_t* __restrict popcount_intersect_ptr) {
  // ct must be a multiple of 3.
  assert(!(vec_ct % 3));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  VecW acc1 = vecw_setzero();
  VecW acc2 = vecw_setzero();
  VecW acc3 = vecw_setzero();
  uintptr_t cur_incr = 30;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 30) {
      if (!vec_ct) {
        *popcount1_ptr = HsumW(acc1);
        *popcount2_ptr = HsumW(acc2);
        *popcount_intersect_ptr = HsumW(acc3);
        return;
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc1 = vecw_setzero();
    VecW inner_acc2 = vecw_setzero();
    VecW inner_acc3 = vecw_setzero();
    const VecW* vvec1_stop = &(vvec1_iter[cur_incr]);
    do {
      VecW count1a = *vvec1_iter++;
      VecW count2a = *vvec2_iter++;
      VecW count3a = count1a & count2a;
      VecW count1b = *vvec1_iter++;
      VecW count2b = *vvec2_iter++;
      VecW count3b = count1b & count2b;
      VecW half1a = *vvec1_iter++;
      VecW half2a = *vvec2_iter++;
      const VecW half1b = vecw_srli(half1a, 1) & m1;
      const VecW half2b = vecw_srli(half2a, 1) & m1;
      half1a = half1a & m1;
      half2a = half2a & m1;

      count1a = count1a - (vecw_srli(count1a, 1) & m1);
      count2a = count2a - (vecw_srli(count2a, 1) & m1);
      count3a = count3a - (vecw_srli(count3a, 1) & m1);
      count1b = count1b - (vecw_srli(count1b, 1) & m1);
      count2b = count2b - (vecw_srli(count2b, 1) & m1);
      count3b = count3b - (vecw_srli(count3b, 1) & m1);
      count1a = count1a + half1a;
      count2a = count2a + half2a;
      count3a = count3a + (half1a & half2a);
      count1b = count1b + half1b;
      count2b = count2b + half2b;
      count3b = count3b + (half1b & half2b);

      count1a = (count1a & m2) + (vecw_srli(count1a, 2) & m2);
      count2a = (count2a & m2) + (vecw_srli(count2a, 2) & m2);
      count3a = (count3a & m2) + (vecw_srli(count3a, 2) & m2);
      count1a = count1a + (count1b & m2) + (vecw_srli(count1b, 2) & m2);
      count2a = count2a + (count2b & m2) + (vecw_srli(count2b, 2) & m2);
      count3a = count3a + (count3b & m2) + (vecw_srli(count3b, 2) & m2);
      inner_acc1 = inner_acc1 + (count1a & m4) + (vecw_srli(count1a, 4) & m4);
      inner_acc2 = inner_acc2 + (count2a & m4) + (vecw_srli(count2a, 4) & m4);
      inner_acc3 = inner_acc3 + (count3a & m4) + (vecw_srli(count3a, 4) & m4);
    } while (vvec1_iter < vvec1_stop);
    // too much register pressure to use prev_sad_result pattern?
    const VecW m0 = vecw_setzero();
    acc1 = acc1 + vecw_bytesum(inner_acc1, m0);
    acc2 = acc2 + vecw_bytesum(inner_acc2, m0);
    acc3 = acc3 + vecw_bytesum(inner_acc3, m0);
  }
}

void PopcountWordsIntersect3val(const uintptr_t* __restrict bitvec1, const uintptr_t* __restrict bitvec2, uint32_t word_ct, uint32_t* __restrict popcount1_ptr, uint32_t* __restrict popcount2_ptr, uint32_t* __restrict popcount_intersect_ptr) {
  const uint32_t trivec_ct = word_ct / (3 * kWordsPerVec);
  uint32_t ct1;
  uint32_t ct2;
  uint32_t ct3;
  PopcountVecsNoSse42Intersect3val(R_CAST(const VecW*, bitvec1), R_CAST(const VecW*, bitvec2), trivec_ct * 3, &ct1, &ct2, &ct3);
  const uint32_t words_consumed = trivec_ct * (3 * kWordsPerVec);
  bitvec1 = &(bitvec1[words_consumed]);
  bitvec2 = &(bitvec2[words_consumed]);
  const uint32_t remainder = word_ct - words_consumed;
  for (uint32_t widx = 0; widx != remainder; ++widx) {
    const uintptr_t word1 = bitvec1[widx];
    const uintptr_t word2 = bitvec2[widx];
    ct1 += PopcountWord(word1);
    ct2 += PopcountWord(word2);
    ct3 += PopcountWord(word1 & word2);
  }
  *popcount1_ptr = ct1;
  *popcount2_ptr = ct2;
  *popcount_intersect_ptr = ct3;
}
#endif

uint32_t AllBitsAreZero(const uintptr_t* bitarr, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t start_idxl = start_idx / kBitsPerWord;
  const uintptr_t start_idxlr = start_idx & (kBitsPerWord - 1);
  const uintptr_t end_idxl = end_idx / kBitsPerWord;
  const uintptr_t end_idxlr = end_idx & (kBitsPerWord - 1);
  if (start_idxl == end_idxl) {
    return !(bitarr[start_idxl] & ((k1LU << end_idxlr) - (k1LU << start_idxlr)));
  }
  if (start_idxlr) {
    if (bitarr[start_idxl++] >> start_idxlr) {
      return 0;
    }
  }
  for (; start_idxl != end_idxl; ++start_idxl) {
    if (bitarr[start_idxl]) {
      return 0;
    }
  }
  if (!end_idxlr) {
    return 1;
  }
  return !bzhi(bitarr[end_idxl], end_idxlr);
}

void CopyBitarrRange(const uintptr_t* __restrict src_bitarr, uintptr_t src_start_bitidx, uintptr_t target_start_bitidx, uintptr_t len, uintptr_t* __restrict target_bitarr) {
  // assumes len is positive, and relevant bits of target_bitarr are zero
  const uintptr_t* src_bitarr_iter = &(src_bitarr[src_start_bitidx / kBitsPerWord]);
  uint32_t src_rshift = src_start_bitidx % kBitsPerWord;
  uintptr_t* target_bitarr_iter = &(target_bitarr[target_start_bitidx / kBitsPerWord]);
  uint32_t target_initial_lshift = target_start_bitidx % kBitsPerWord;
  uintptr_t cur_src_word;
  if (target_initial_lshift) {
    const uint32_t initial_copy_bitct = kBitsPerWord - target_initial_lshift;
    if (len <= initial_copy_bitct) {
      goto CopyBitarrRange_last_partial_word;
    }
    cur_src_word = (*src_bitarr_iter) >> src_rshift;
    if (src_rshift >= target_initial_lshift) {
      ++src_bitarr_iter;
      cur_src_word |= (*src_bitarr_iter) << (kBitsPerWord - src_rshift);
    }
    *target_bitarr_iter++ |= cur_src_word << target_initial_lshift;
    src_rshift = (src_rshift + initial_copy_bitct) % kBitsPerWord;
    len -= initial_copy_bitct;
  }
  {
    const uintptr_t fullword_ct = len / kBitsPerWord;
    if (!src_rshift) {
      memcpy(target_bitarr_iter, src_bitarr_iter, fullword_ct * sizeof(intptr_t));
      target_bitarr_iter = &(target_bitarr_iter[fullword_ct]);
      src_bitarr_iter = &(src_bitarr_iter[fullword_ct]);
    } else {
      const uint32_t src_lshift = kBitsPerWord - src_rshift;
      cur_src_word = *src_bitarr_iter;
      for (uintptr_t widx = 0; widx != fullword_ct; ++widx) {
        const uintptr_t next_src_word = *(++src_bitarr_iter);
        *target_bitarr_iter++ = (cur_src_word >> src_rshift) | (next_src_word << src_lshift);
        cur_src_word = next_src_word;
      }
    }
  }
  len %= kBitsPerWord;
  if (len) {
    target_initial_lshift = 0;
  CopyBitarrRange_last_partial_word:
    cur_src_word = (*src_bitarr_iter) >> src_rshift;
    if (len + src_rshift > kBitsPerWord) {
      cur_src_word |= src_bitarr_iter[1] << (kBitsPerWord - src_rshift);
    }
    *target_bitarr_iter |= (cur_src_word & ((~k0LU) >> (kBitsPerWord - S_CAST(uint32_t, len)))) << target_initial_lshift;
  }
}

// advances forward_ct set bits; forward_ct must be positive.  (stays put if
// forward_ct == 1 and current bit is set.  may want to tweak this interface,
// easy to introduce off-by-one bugs...)
// In usual 64-bit case, also assumes bitvec is 16-byte aligned and the end of
// the trailing 16-byte block can be safely read from.
uintptr_t FindNth1BitFrom(const uintptr_t* bitvec, uintptr_t cur_pos, uintptr_t forward_ct) {
  assert(forward_ct);
  uintptr_t widx = cur_pos / kBitsPerWord;
  uintptr_t ulii = cur_pos % kBitsPerWord;
  const uintptr_t* bptr = &(bitvec[widx]);
  uintptr_t uljj;
  uintptr_t ulkk;
#ifdef __LP64__
  const VecW* vptr;
  assert(VecIsAligned(bitvec));
#endif
  if (ulii) {
    uljj = (*bptr) >> ulii;
    ulkk = PopcountWord(uljj);
    if (ulkk >= forward_ct) {
    JumpForwardSetUnsafe_finish:
      return widx * kBitsPerWord + ulii + WordBitIdxToUidx(uljj, forward_ct - 1);
    }
    forward_ct -= ulkk;
    ++widx;
    ++bptr;
  }
  ulii = 0;
#ifdef __LP64__
  while (widx & (kWordsPerVec - k1LU)) {
    uljj = *bptr;
    ulkk = PopcountWord(uljj);
    if (ulkk >= forward_ct) {
      goto JumpForwardSetUnsafe_finish;
    }
    forward_ct -= ulkk;
    ++widx;
    ++bptr;
  }
  vptr = R_CAST(const VecW*, bptr);
#ifdef USE_AVX2
  while (forward_ct > kBitsPerWord * (16 * kWordsPerVec)) {
    uljj = (forward_ct - 1) / (kBitsPerWord * kWordsPerVec);
    ulkk = PopcountVecsAvx2(vptr, uljj);
    vptr = &(vptr[uljj]);
    forward_ct -= ulkk;
  }
#else
  while (forward_ct > kBitsPerWord * (3 * kWordsPerVec)) {
    uljj = ((forward_ct - 1) / (kBitsPerWord * (3 * kWordsPerVec))) * 3;
    // yeah, yeah, this is suboptimal if we have SSE4.2
    ulkk = PopcountVecsNoAvx2(vptr, uljj);
    vptr = &(vptr[uljj]);
    forward_ct -= ulkk;
  }
#endif
  bptr = R_CAST(const uintptr_t*, vptr);
  while (forward_ct > kBitsPerWord) {
    forward_ct -= PopcountWord(*bptr++);
  }
#else
  while (forward_ct > kBitsPerWord) {
    uljj = (forward_ct - 1) / kBitsPerWord;
    ulkk = PopcountWords(bptr, uljj);
    bptr = &(bptr[uljj]);
    forward_ct -= ulkk;
  }
#endif
  for (; ; ++bptr) {
    uljj = *bptr;
    ulkk = PopcountWord(uljj);
    if (ulkk >= forward_ct) {
      widx = bptr - bitvec;
      goto JumpForwardSetUnsafe_finish;
    }
    forward_ct -= ulkk;
  }
}

void ComputeUidxStartPartition(const uintptr_t* variant_include, uint64_t variant_ct, uint32_t thread_ct, uint32_t first_variant_uidx, uint32_t* variant_uidx_starts) {
  assert(variant_ct);
  uint32_t cur_variant_uidx_start = AdvTo1Bit(variant_include, first_variant_uidx);
  uint32_t cur_variant_idx_start = 0;
  variant_uidx_starts[0] = cur_variant_uidx_start;
  for (uint32_t tidx = 1; tidx != thread_ct; ++tidx) {
    const uint32_t new_variant_idx_start = (tidx * variant_ct) / thread_ct;
    if (new_variant_idx_start != cur_variant_idx_start) {
      cur_variant_uidx_start = FindNth1BitFrom(variant_include, cur_variant_uidx_start + 1, new_variant_idx_start - cur_variant_idx_start);
      cur_variant_idx_start = new_variant_idx_start;
    }
    variant_uidx_starts[tidx] = cur_variant_uidx_start;
  }
}

// May want to have an multiallelic_set bitarray to accelerate this type of
// operation?  Probably only want to conditionally initialize it, and only
// after variant filtering is complete, though.
uintptr_t CountExtraAlleles(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t multiallelic_variant_ct_multiplier) {
  if (!allele_idx_offsets) {
    return 0;
  }
  const uint32_t variant_ct = PopcountBitRange(variant_include, variant_uidx_start, variant_uidx_end);
  const uint32_t const_subtract = 2 - multiallelic_variant_ct_multiplier;
  uintptr_t result = 0;
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, variant_uidx_start, &variant_uidx_base, &cur_bits);
  for (uint32_t uii = 0; uii != variant_ct; ++uii) {
    const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    const uintptr_t allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
    if (allele_ct != 2) {
      result += allele_ct - const_subtract;
    }
  }
  return result;
}

uint32_t MaxAlleleCtSubset(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct) {
  if (!allele_idx_offsets) {
    return 2;
  }
  if (raw_variant_ct == variant_ct) {
    return max_allele_ct;
  }
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  uintptr_t subset_max_allele_ct = 2;
  for (uint32_t uii = 0; uii != variant_ct; ++uii) {
    const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    const uint32_t allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
    if (allele_ct > subset_max_allele_ct) {
      if (allele_ct == max_allele_ct) {
        return max_allele_ct;
      }
      subset_max_allele_ct = allele_ct;
    }
  }
  return subset_max_allele_ct;
}

void ComputePartitionAligned(const uintptr_t* variant_include, uint32_t orig_thread_ct, uint32_t first_variant_uidx, uint32_t cur_variant_idx, uint32_t cur_variant_ct, uint32_t alignment, uint32_t* variant_uidx_starts, uint32_t* vidx_starts) {
  // Minimize size of the largest chunk, under the condition that all
  // intermediate variant_idx values are divisible by alignment.
  //
  // There are three possibilities:
  // 1. The straightforward solution one gets from rounding cur_variant_idx
  //    down, and (cur_variant_idx + variant_ct) up, is optimal.  Call this
  //    chunk size C.
  // 2. C - leading_idx_ct is ideal.
  // 3. C - trailing_idx_ct is ideal.

  assert(cur_variant_ct);
  const uint32_t variant_idx_end = cur_variant_idx + cur_variant_ct;
  const uint32_t log2_align = ctzu32(alignment);
  const uint32_t leading_idx_ct = cur_variant_idx % alignment;
  const uint32_t trailing_idx_ct = (-variant_idx_end) % alignment;
  const uint32_t block_ct = (cur_variant_ct + leading_idx_ct + trailing_idx_ct) >> log2_align;
  uint32_t cur_variant_uidx_start = AdvTo1Bit(variant_include, first_variant_uidx);
  variant_uidx_starts[0] = cur_variant_uidx_start;
  vidx_starts[0] = cur_variant_idx;
  const uint32_t thread_ct = MINV(orig_thread_ct, block_ct);
  if (thread_ct > 1) {
    const uint32_t std_blocks_per_thread = 1 + (block_ct - 1) / thread_ct;
    // Possibilities 2 and 3 are only live if
    //   block_ct == (std_blocks_per_thread - 1) * thread_ct + 1, or
    //   block_ct == (std_blocks_per_thread - 1) * thread_ct + 2,
    // In the first situation, the best solution is
    //    min(possibility 2, possibility 3).
    // In the second situation, the best solution is
    //    max(possibility 2, possibility 3).
    uint32_t central_variant_ct = std_blocks_per_thread * alignment;
    uint32_t first_block_variant_ct = central_variant_ct - leading_idx_ct;
    const uint32_t remainder_m1 = block_ct - (std_blocks_per_thread - 1) * thread_ct - 1;
    if (remainder_m1 <= 1) {
      central_variant_ct -= alignment;
      if ((!remainder_m1) && (leading_idx_ct < trailing_idx_ct)) {
        first_block_variant_ct -= alignment;
      }
    }
    cur_variant_uidx_start = FindNth1BitFrom(variant_include, cur_variant_uidx_start + 1, first_block_variant_ct);
    cur_variant_idx += first_block_variant_ct;
    variant_uidx_starts[1] = cur_variant_uidx_start;
    vidx_starts[1] = cur_variant_idx;
    for (uint32_t tidx = 2; tidx != thread_ct; ++tidx) {
      cur_variant_uidx_start = FindNth1BitFrom(variant_include, cur_variant_uidx_start + 1, central_variant_ct);
      cur_variant_idx += central_variant_ct;
      // bugfix (14 Nov 2017): this decrement was in the wrong place
      if (tidx == remainder_m1) {
        central_variant_ct -= alignment;
      }
      variant_uidx_starts[tidx] = cur_variant_uidx_start;
      vidx_starts[tidx] = cur_variant_idx;
    }
  }
  if (thread_ct < orig_thread_ct) {
    uint32_t last_vidx_ct = variant_idx_end - vidx_starts[thread_ct - 1];
    cur_variant_uidx_start = FindNth1BitFrom(variant_include, cur_variant_uidx_start + 1, last_vidx_ct);
    for (uint32_t tidx = thread_ct; tidx != orig_thread_ct; ++tidx) {
      variant_uidx_starts[tidx] = cur_variant_uidx_start;
    }
    for (uint32_t tidx = thread_ct; tidx != orig_thread_ct; ++tidx) {
      vidx_starts[tidx] = variant_idx_end;
    }
  }
  vidx_starts[orig_thread_ct] = variant_idx_end;
}

BoolErr ParseNextRange(const char* const* argvk, uint32_t param_ct, char range_delim, uint32_t* cur_param_idx_ptr, const char** cur_arg_pptr, const char** range_start_ptr, uint32_t* rs_len_ptr, const char** range_end_ptr, uint32_t* re_len_ptr) {
  // Starts reading from argv[cur_param_idx][cur_pos].  If a valid range is
  // next, range_start + rs_len + range_end + re_len are updated.  If only a
  // single item is next, range_end is set to nullptr and range_start + rs_len
  // are updated.  If there are no items left, range_start is set to nullptr.
  // If the input is not well-formed, -1 is returned instead of 0.
  uint32_t cur_param_idx = *cur_param_idx_ptr;
  if (cur_param_idx > param_ct) {
    *cur_arg_pptr = nullptr;
    return 0;
  }
  const char* cur_arg_ptr = *cur_arg_pptr;
  while (1) {
    char cc = *cur_arg_ptr;
    if (!cc) {
      *cur_param_idx_ptr = ++cur_param_idx;
      if (cur_param_idx > param_ct) {
        *range_start_ptr = nullptr;
        return 0;
      }
      cur_arg_ptr = argvk[cur_param_idx];
      cc = *cur_arg_ptr;
    }
    if (unlikely(cc == range_delim)) {
      return 1;
    }
    if (cc != ',') {
      break;
    }
    ++cur_arg_ptr;
  }
  *range_start_ptr = cur_arg_ptr;
  char cc;
  do {
    cc = *(++cur_arg_ptr);
    if ((!cc) || (cc == ',')) {
      *rs_len_ptr = cur_arg_ptr - (*range_start_ptr);
      *cur_arg_pptr = cur_arg_ptr;
      *range_end_ptr = nullptr;
      return 0;
    }
  } while (cc != range_delim);
  *rs_len_ptr = cur_arg_ptr - (*range_start_ptr);
  cc = *(++cur_arg_ptr);
  if (unlikely(((*rs_len_ptr) > kMaxIdSlen) || (!cc) || (cc == ',') || (cc == range_delim))) {
    return 1;
  }
  *range_end_ptr = cur_arg_ptr;
  do {
    cc = *(++cur_arg_ptr);
    if (unlikely(cc == range_delim)) {
      return 1;
    }
  } while (cc && (cc != ','));
  *re_len_ptr = cur_arg_ptr - (*range_end_ptr);
  if (unlikely((*re_len_ptr) > kMaxIdSlen)) {
    return 1;
  }
  *cur_arg_pptr = cur_arg_ptr;
  return 0;
}

PglErr ParseNameRanges(const char* const* argvk, const char* errstr_append, uint32_t param_ct, uint32_t require_posint, char range_delim, RangeList* range_list_ptr) {
  uint32_t name_ct = 0;
  uint32_t cur_param_idx = 1;
  uint32_t name_max_blen = 0;
  const char* cur_arg_ptr;
  const char* range_start;
  uint32_t rs_len;
  const char* range_end;
  uint32_t re_len;
  unsigned char* cur_name_starts_range;
  uint32_t last_val;
  uint32_t cur_val;
  // two passes.  first pass: count arguments, determine name_max_blen;
  // then allocate memory; then fill it.
  if (param_ct) {
    cur_arg_ptr = argvk[1];
    while (1) {
      if (unlikely(ParseNextRange(argvk, param_ct, range_delim, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len))) {
        logerrprintfww("Error: Invalid %s argument '%s'.\n", argvk[0], argvk[cur_param_idx]);
        logerrputs(errstr_append);
        return kPglRetInvalidCmdline;
      }
      if (!range_start) {
        break;
      }
      ++name_ct;
      if (rs_len > name_max_blen) {
        name_max_blen = rs_len;  // does NOT include trailing null yet
      }
      if (range_end) {
        ++name_ct;
        if (re_len > name_max_blen) {
          name_max_blen = re_len;
        }
      }
    }
  }
  if (unlikely(!name_ct)) {
    logerrprintf("Error: %s requires at least one value.\n%s", argvk[0], errstr_append);
    return kPglRetInvalidCmdline;
  }
  range_list_ptr->name_max_blen = ++name_max_blen;
  range_list_ptr->name_ct = name_ct;
  if (unlikely(pgl_malloc(name_ct * (S_CAST(uintptr_t, name_max_blen) + 1), &range_list_ptr->names))) {
    return kPglRetNomem;
  }
  range_list_ptr->starts_range = R_CAST(unsigned char*, &(range_list_ptr->names[name_ct * S_CAST(uintptr_t, name_max_blen)]));
  char* cur_name_str = range_list_ptr->names;
  cur_name_starts_range = range_list_ptr->starts_range;
  cur_param_idx = 1;
  cur_arg_ptr = argvk[1];
  while (1) {
    // second pass; this can't fail since we already validated
    ParseNextRange(argvk, param_ct, range_delim, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len);
    if (!range_start) {
      if (require_posint) {
        last_val = 0;
        for (cur_param_idx = 0; cur_param_idx != name_ct; ++cur_param_idx) {
          cur_name_str = &(range_list_ptr->names[cur_param_idx * S_CAST(uintptr_t, name_max_blen)]);
          const char* dup_check = cur_name_str;  // actually a numeric check
          do {
            if (unlikely(IsNotDigit(*dup_check))) {
              logerrprintfww("Error: Invalid %s argument '%s'.\n", argvk[0], cur_name_str);
              return kPglRetInvalidCmdline;
            }
          } while (*(++dup_check));
          if (unlikely(ScanPosintDefcapx(cur_name_str, &cur_val))) {
            logerrprintfww("Error: Invalid %s argument '%s'.\n", argvk[0], cur_name_str);
            return kPglRetInvalidCmdline;
          }
          if (range_list_ptr->starts_range[cur_param_idx]) {
            last_val = cur_val;
          } else {
            if (unlikely(cur_val <= last_val)) {
              logerrprintfww("Error: Invalid %s range '%s-%s'.\n", argvk[0], &(range_list_ptr->names[(cur_param_idx - 1) * name_max_blen]), cur_name_str);
              return kPglRetInvalidCmdline;
            }
            last_val = 0;
          }
        }
      }
      return kPglRetSuccess;
    }
    memcpyx(cur_name_str, range_start, rs_len, '\0');
    const char* dup_check = range_list_ptr->names;
    while (dup_check < cur_name_str) {
      if (unlikely(memequal(dup_check, cur_name_str, rs_len + 1))) {
        logerrprintfww("Error: Duplicate %s argument '%s'.\n", argvk[0], cur_name_str);
        return kPglRetInvalidCmdline;
      }
      dup_check = &(dup_check[name_max_blen]);
    }
    cur_name_str = &(cur_name_str[name_max_blen]);
    if (range_end) {
      *cur_name_starts_range++ = 1;
      memcpyx(cur_name_str, range_end, re_len, '\0');
      dup_check = range_list_ptr->names;
      while (dup_check < cur_name_str) {
        if (unlikely(memequal(dup_check, cur_name_str, rs_len + 1))) {
          logerrprintfww("Error: Duplicate %s argument '%s'.\n", argvk[0], cur_name_str);
          return kPglRetInvalidCmdline;
        }
        dup_check = &(dup_check[name_max_blen]);
      }
      cur_name_str = &(cur_name_str[name_max_blen]);
      *cur_name_starts_range++ = 0;
    } else {
      *cur_name_starts_range++ = 0;
    }
  }
}


uint32_t CubicRealRoots(double coef_a, double coef_b, double coef_c, STD_ARRAY_REF(double, 3) solutions) {
  // Additional research into numerical stability may be in order here...
  double a2 = coef_a * coef_a;
  double qq = (a2 - 3 * coef_b) * (1.0 / 9.0);
  double rr = (2 * a2 * coef_a - 9 * coef_a * coef_b + 27 * coef_c) * (1.0 / 54.0);
  double r2 = rr * rr;
  double q3 = qq * qq * qq;
  double adiv3 = coef_a * (1.0 / 3.0);
  double sq;
  double dxx;
  if (r2 < q3) {
    // three real roots
    sq = sqrt(qq);
    dxx = acos(rr / (qq * sq)) * (1.0 / 3.0);
    sq *= -2;
    solutions[0] = sq * cos(dxx) - adiv3;
    solutions[1] = sq * cos(dxx + (2.0 * kPi / 3.0)) - adiv3;
    solutions[2] = sq * cos(dxx - (2.0 * kPi / 3.0)) - adiv3;
    // now sort and check for within-epsilon equality
    if (solutions[0] > solutions[1]) {
      dxx = solutions[0];
      solutions[0] = solutions[1];
      if (dxx > solutions[2]) {
        solutions[1] = solutions[2];
        solutions[2] = dxx;
      } else {
        solutions[1] = dxx;
      }
      if (solutions[0] > solutions[1]) {
        dxx = solutions[0];
        solutions[0] = solutions[1];
        solutions[1] = dxx;
      }
    } else if (solutions[1] > solutions[2]) {
      dxx = solutions[1];
      solutions[1] = solutions[2];
      solutions[2] = dxx;
    }
    if (solutions[1] - solutions[0] < kEpsilon) {
      solutions[1] = solutions[2];
      return 2 - (solutions[1] - solutions[0] < kEpsilon);
    }
    return 3 - (solutions[2] - solutions[1] < kEpsilon);
  }
  dxx = -pow(fabs(rr) + sqrt(r2 - q3), 1.0 / 3.0);
  if (dxx == 0.0) {
    solutions[0] = -adiv3;
    return 1;
  }
  if (rr < 0.0) {
    dxx = -dxx;
  }
  sq = qq / dxx;
  solutions[0] = dxx + sq - adiv3;
  // use of regular epsilon here has actually burned us
  if (fabs(dxx - sq) >= (kEpsilon * 8)) {
    return 1;
  }
  if (dxx >= 0.0) {
    solutions[1] = solutions[0];
    solutions[0] = -dxx - adiv3;
  } else {
    solutions[1] = -dxx - adiv3;
  }
  return 2;
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
CONSTI32(kDupstoreBlockSize, 65536);
CONSTI32(kDupstoreThreadWkspace, kDupstoreBlockSize * 2 * sizeof(int32_t));

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
//
// It will probably be moved out of plink2_cmdline soon.
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
typedef struct IdUniquenessCheckerStruct {
  const uintptr_t* subset_mask;
  const char* const* item_ids;
  uint32_t item_ct;
  uint32_t id_htable_size;
  uint32_t item_uidx_starts[kMaxDupflagThreads];

  uint32_t* id_htable;

  uint32_t dup_found;
} IdUniquenessChecker;

void IdUniquenessCheckerMain(uint32_t tidx, uint32_t thread_ct, IdUniquenessChecker* ctx) {
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

THREAD_FUNC_DECL IdUniquenessCheckerThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  IdUniquenessChecker* ctx = S_CAST(IdUniquenessChecker*, arg->sharedp->context);

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
  IdUniquenessCheckerMain(tidx, thread_ct, ctx);
  THREAD_RETURN;
}

PglErr CheckIdUniqueness(unsigned char* arena_bottom, unsigned char* arena_top, const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uint32_t max_thread_ct, uint32_t* dup_found_ptr) {
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  IdUniquenessChecker ctx;
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
      goto CheckIdUniqueness_ret_NOMEM;
    }

    ctx.subset_mask = subset_mask;
    ctx.item_ids = item_ids;
    ctx.item_ct = item_ct;
    uint32_t id_htable_size = GetHtableFastSize(item_ct);
    {
      uintptr_t wkspace_left = S_CAST(uintptr_t, arena_top - arena_bottom);
      if (wkspace_left < id_htable_size * sizeof(int32_t)) {
        id_htable_size = wkspace_left / sizeof(int32_t);
        const uint32_t min_htable_size = GetHtableMinSize(item_ct);
        if (id_htable_size < min_htable_size) {
          goto CheckIdUniqueness_ret_NOMEM;
        }
      }
    }
    ctx.id_htable_size = id_htable_size;
    uint32_t* id_htable = R_CAST(uint32_t*, arena_bottom);
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
      SetThreadFuncAndData(IdUniquenessCheckerThread, &ctx, &tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto CheckIdUniqueness_ret_THREAD_CREATE_FAIL;
      }
    }
    const uint32_t fill_start = RoundDownPow2((id_htable_size * S_CAST(uint64_t, thread_ct - 1)) / thread_ct, kInt32PerCacheline);
    SetAllU32Arr(id_htable_size - fill_start, &(id_htable[fill_start]));
    if (thread_ct > 1) {
      JoinThreads(&tg);
      DeclareLastThreadBlock(&tg);
      SpawnThreads(&tg);
    }
    IdUniquenessCheckerMain(thread_ct - 1, thread_ct, &ctx);
    JoinThreads0(&tg);
    *dup_found_ptr = ctx.dup_found;
  }
  while (0) {
  CheckIdUniqueness_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CheckIdUniqueness_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  return reterr;
}

PglErr AllocAndPopulateIdHtableMt(const uintptr_t* subset_mask, const char* const* item_ids, uintptr_t item_ct, uintptr_t fast_size_min_extra_bytes, uint32_t max_thread_ct, uint32_t** id_htable_ptr, uint32_t** htable_dup_base_ptr, uint32_t* id_htable_size_ptr, uint32_t* dup_ct_ptr) {
  uint32_t id_htable_size = GetHtableFastSize(item_ct);
  // if store_all_dups, up to 8 bytes per variant in extra_alloc for duplicate
  //   tracking, as well as a small amount of per-thread temporary workspace
  const uint32_t store_all_dups = (htable_dup_base_ptr != nullptr);
  uintptr_t nonhtable_alloc = 0;
  if (store_all_dups) {
    const uint32_t actual_thread_ct = PopulateIdHtableMtDupstoreThreadCt(max_thread_ct, item_ct);
    // "2 *" in front of second term is temporary
    nonhtable_alloc = actual_thread_ct * kDupstoreThreadWkspace + 2 * RoundUpPow2(item_ct * 2 * sizeof(int32_t), kCacheline);
  }
  uintptr_t max_bytes = RoundDownPow2(bigstack_left(), kCacheline);
  // force max_bytes >= 5 so leqprime() doesn't fail
  // (not actually relevant any more, but whatever)
  if (nonhtable_alloc + (item_ct + 6) * sizeof(int32_t) > max_bytes) {
    return kPglRetNomem;
  }
  max_bytes -= nonhtable_alloc;
  if (id_htable_size * sizeof(int32_t) + fast_size_min_extra_bytes > max_bytes) {
    id_htable_size = max_bytes / sizeof(int32_t);
    // id_htable_size = leqprime((max_bytes / sizeof(int32_t)) - 1);
    const uint32_t min_htable_size = GetHtableMinSize(item_ct);
    if (id_htable_size < min_htable_size) {
      id_htable_size = min_htable_size;
    }
  }
  *id_htable_ptr = S_CAST(uint32_t*, bigstack_alloc_raw_rd(id_htable_size * sizeof(int32_t)));
  if (store_all_dups) {
    *htable_dup_base_ptr = &((*id_htable_ptr)[RoundUpPow2(id_htable_size, kInt32PerCacheline)]);
  }
  *id_htable_size_ptr = id_htable_size;
  return PopulateIdHtableMt(g_bigstack_end, subset_mask, item_ids, item_ct, store_all_dups, id_htable_size, max_thread_ct, &g_bigstack_base, *id_htable_ptr, dup_ct_ptr);
}


uint32_t Edit1Match(const char* s1, const char* s2, uint32_t len1, uint32_t len2) {
  // Permit one difference of the following forms (Damerau-Levenshtein distance
  // 1):
  // - inserted/deleted character
  // - replaced character
  // - adjacent pair of swapped characters
  // Returns 1 on approximate match, 0 on nonmatch.
  // With only one difference allowed, it is unnecessary to use a
  // fully-powered quadratic-time DP algorithm; instead, given the lengths of
  // the two strings, we select the appropriate linear-time loop.
  if (abs_i32(len1 - len2) > 1) {
    return 0;
  }
  uintptr_t mismatch_pos = FirstUnequal(s1, s2, len1);
  if (mismatch_pos == len1) {
    // bugfix (5 Nov 2019): this was backwards...
    return 1;
  }
  if (len1 != len2) {
    // At least one mismatch guaranteed.
    if (len1 < len2) {
      return memequal(&(s1[mismatch_pos]), &(s2[mismatch_pos + 1]), len1 - mismatch_pos);
    }
    return memequal(&(s1[mismatch_pos + 1]), &(s2[mismatch_pos]), len2 - mismatch_pos);
  }
  // Strings are equal length, and we have at least one mismatch.
  // Two ways to still have an approximate match:
  // 1. Remainder of strings are equal.
  // 2. Next character also differs, but due to transposition.  Remainder of
  //    strings are equal.
  if (s1[mismatch_pos + 1] != s2[mismatch_pos + 1]) {
    if ((s1[mismatch_pos] != s2[mismatch_pos + 1]) || (s1[mismatch_pos + 1] != s2[mismatch_pos])) {
      return 0;
    }
    ++mismatch_pos;
  }
  const uintptr_t post_mismatch_pos = mismatch_pos + 1;
  return memequal(&(s1[post_mismatch_pos]), &(s2[post_mismatch_pos]), len1 - post_mismatch_pos);
}

CONSTI32(kMaxEqualHelpParams, 64);

void HelpPrint(const char* cur_params, HelpCtrl* help_ctrl_ptr, uint32_t postprint_newline, const char* payload) {
  if (help_ctrl_ptr->param_ct) {
    const char* params_iter = cur_params;
    uint32_t cur_param_ct = 0;
    const char* cur_param_starts[kMaxEqualHelpParams];
    uint32_t cur_param_slens[kMaxEqualHelpParams];
    while (*params_iter) {
      cur_param_starts[cur_param_ct] = params_iter;
      const uint32_t cur_slen = strlen(params_iter);
      cur_param_slens[cur_param_ct++] = cur_slen;
      params_iter = &(params_iter[cur_slen + 1]);
    }
    if (help_ctrl_ptr->iters_left) {
      const uint32_t orig_unmatched_ct = help_ctrl_ptr->unmatched_ct;
      if (orig_unmatched_ct) {
        uint32_t arg_uidx = 0;
        if (help_ctrl_ptr->iters_left == 2) {
          // First pass: only exact matches.
          for (uint32_t arg_idx = 0; arg_idx != orig_unmatched_ct; ++arg_idx, ++arg_uidx) {
            arg_uidx = AdvTo0Bit(help_ctrl_ptr->all_match_arr, arg_uidx);
            const char* cur_arg = help_ctrl_ptr->argv[arg_uidx];
            const uint32_t cur_arg_slen = help_ctrl_ptr->param_slens[arg_uidx];
            for (uint32_t cur_param_idx = 0; cur_param_idx != cur_param_ct; ++cur_param_idx) {
              if ((cur_arg_slen == cur_param_slens[cur_param_idx]) && memequal(cur_arg, cur_param_starts[cur_param_idx], cur_arg_slen)) {
                SetBit(arg_uidx, help_ctrl_ptr->perfect_match_arr);
                SetBit(arg_uidx, help_ctrl_ptr->prefix_match_arr);
                SetBit(arg_uidx, help_ctrl_ptr->all_match_arr);
                help_ctrl_ptr->unmatched_ct -= 1;
                break;
              }
            }
          }
        } else {
          // Second pass: Prefix matches.
          for (uint32_t arg_idx = 0; arg_idx != orig_unmatched_ct; ++arg_idx, ++arg_uidx) {
            arg_uidx = AdvTo0Bit(help_ctrl_ptr->all_match_arr, arg_uidx);
            const char* cur_arg = help_ctrl_ptr->argv[arg_uidx];
            const uint32_t cur_arg_slen = help_ctrl_ptr->param_slens[arg_uidx];
            for (uint32_t cur_param_idx = 0; cur_param_idx != cur_param_ct; ++cur_param_idx) {
              if (cur_param_slens[cur_param_idx] > cur_arg_slen) {
                if (memequal(cur_arg, cur_param_starts[cur_param_idx], cur_arg_slen)) {
                  SetBit(arg_uidx, help_ctrl_ptr->prefix_match_arr);
                  SetBit(arg_uidx, help_ctrl_ptr->all_match_arr);
                  help_ctrl_ptr->unmatched_ct -= 1;
                  break;
                }
              }
            }
          }
        }
      }
    } else {
      uint32_t print_this = 0;
      for (uint32_t arg_uidx = 0; arg_uidx != help_ctrl_ptr->param_ct; ++arg_uidx) {
        if (IsSet(help_ctrl_ptr->prefix_match_arr, arg_uidx)) {
          // Current user-provided help argument has at least one prefix
          // match...
          if (!print_this) {
            const char* cur_arg = help_ctrl_ptr->argv[arg_uidx];
            const uint32_t cur_arg_slen = help_ctrl_ptr->param_slens[arg_uidx];
            if (IsSet(help_ctrl_ptr->perfect_match_arr, arg_uidx)) {
              // ...and at least one exact match.  So we only care about an
              // exact match.
              for (uint32_t cur_param_idx = 0; cur_param_idx != cur_param_ct; ++cur_param_idx) {
                if ((cur_arg_slen == cur_param_slens[cur_param_idx]) && memequal(cur_arg, cur_param_starts[cur_param_idx], cur_arg_slen)) {
                  print_this = 1;
                  break;
                }
              }
            } else {
              for (uint32_t cur_param_idx = 0; cur_param_idx != cur_param_ct; ++cur_param_idx) {
                if (cur_param_slens[cur_param_idx] > cur_arg_slen) {
                  if (memequal(cur_arg, cur_param_starts[cur_param_idx], cur_arg_slen)) {
                    print_this = 1;
                    break;
                  }
                }
              }
            }
          }
        } else {
          // Current user-provided help argument does not have an exact or
          // prefix match.  Print current help text if it's within
          // Damerau-Levenshtein distance 1 of any index term.
          for (uint32_t cur_param_idx = 0; cur_param_idx != cur_param_ct; ++cur_param_idx) {
            if (Edit1Match(cur_param_starts[cur_param_idx], help_ctrl_ptr->argv[arg_uidx], cur_param_slens[cur_param_idx], help_ctrl_ptr->param_slens[arg_uidx])) {
              print_this = 1;
              if (!IsSet(help_ctrl_ptr->all_match_arr, arg_uidx)) {
                SetBit(arg_uidx, help_ctrl_ptr->all_match_arr);
                help_ctrl_ptr->unmatched_ct -= 1;
              }
              break;
            }
          }
        }
      }
      if (print_this) {
        const uint32_t payload_slen = strlen(payload);
        const char* payload_end;
        if (payload[payload_slen - 2] == '\n') {
          payload_end = &(payload[payload_slen - 1]);
        } else {
          payload_end = &(payload[payload_slen]);
        }
        if (help_ctrl_ptr->preprint_newline) {
          putc_unlocked('\n', stdout);
        }
        help_ctrl_ptr->preprint_newline = postprint_newline;
        const char* payload_iter = payload;
        do {
          const char* line_end = AdvPastDelim(payload_iter, '\n');
          uint32_t line_slen = S_CAST(uint32_t, line_end - payload_iter);
          if (line_slen > 2) {
            payload_iter = &(payload_iter[2]);
            line_slen -= 2;
          }
          memcpyx(g_textbuf, payload_iter, line_slen, '\0');
          fputs(g_textbuf, stdout);
          payload_iter = line_end;
        } while (payload_iter < payload_end);
      }
    }
  } else {
    fputs(payload, stdout);
  }
}


void PreinitPlink2CmdlineMeta(Plink2CmdlineMeta* pcmp) {
  pcmp->subst_argv = nullptr;
  pcmp->script_buf = nullptr;
  pcmp->rerun_buf = nullptr;
  pcmp->flag_buf = nullptr;
  pcmp->flag_map = nullptr;
}

const char kErrstrNomem[] = "Error: Out of memory.  The --memory flag may be helpful.\n";
const char kErrstrWrite[] = "Error: File write failure: %s.\n";
const char kErrstrThreadCreate[] = "Error: Failed to create thread.\n";
const char kErrstrVarRecordTooLarge[] = "Error: Variant record size exceeds ~4 GiB limit.\n";
const char kErrstrReadCorrupted[] = "File appears to be corrupted";

// assumes logfile is open
void DispExitMsg(PglErr reterr) {
  if (reterr) {
    if (reterr == kPglRetNomem) {
      logputs("\n");
      logerrputs(kErrstrNomem);
      if (g_failed_alloc_attempt_size) {
        logerrprintf("Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
      }
      // kPglRetReadFail no longer gets a message here, for the same reason
      // kPglRetOpenFail doesn't: it's important to know which file we failed
      // to read.
    } else if (reterr == kPglRetWriteFail) {
      logputs("\n");
      if (errno) {
        logerrprintfww(kErrstrWrite, strerror(errno));
      } else {
        // Defensive.
        logerrputs("Error: File write failure: Untracked cause.\n");
      }
    } else if (reterr == kPglRetThreadCreateFail) {
      logputs("\n");
      logerrputs(kErrstrThreadCreate);
    } else if (reterr == kPglRetVarRecordTooLarge) {
      logputs("\n");
      logerrputs(kErrstrVarRecordTooLarge);
    } else if (reterr == kPglRetLongLine) {
      logputs("\n");
      logerrprintf("Error: Unhandled internal line-too-long message.\n");
    } else if (reterr == kPglRetEof) {
      logputs("\n");
      logerrprintf("Error: Unhandled internal EOF message.\n");
    }
  }
}

// useful when there's e.g. a filename and an optional modifier, and we want to
// permit either parmeter ordering
BoolErr CheckExtraParam(const char* const* argvk, const char* permitted_modif, uint32_t* other_idx_ptr) {
  const uint32_t idx_base = *other_idx_ptr;
  if (!strcmp(argvk[idx_base], permitted_modif)) {
    *other_idx_ptr = idx_base + 1;
  } else if (strcmp(argvk[idx_base + 1], permitted_modif)) {
    logerrprintf("Error: Invalid %s argument sequence.\n", argvk[0]);
    return 1;
  }
  return 0;
}

char ExtractCharParam(const char* ss) {
  // maps c, 'c', and "c" to c, and anything else to the null char.  This is
  // intended to support e.g. always using '#' to designate a # argument
  // without worrying about differences between shells.
  const char cc = ss[0];
  if (((cc == '\'') || (cc == '"')) && (ss[1]) && (ss[2] == cc) && (!ss[3])) {
    return ss[1];
  }
  if (cc && (!ss[1])) {
    return cc;
  }
  return '\0';
}

PglErr CmdlineAllocString(const char* source, const char* flag_name, uint32_t max_slen, char** sbuf_ptr) {
  const uint32_t slen = strlen(source);
  if (slen > max_slen) {
    logerrprintf("Error: %s argument too long.\n", flag_name);
    return kPglRetInvalidCmdline;
  }
  const uint32_t blen = slen + 1;
  if (pgl_malloc(blen, sbuf_ptr)) {
    return kPglRetNomem;
  }
  memcpy(*sbuf_ptr, source, blen);
  return kPglRetSuccess;
}

PglErr AllocFname(const char* source, const char* flagname_p, uint32_t extra_size, char** fnbuf_ptr) {
  const uint32_t blen = strlen(source) + 1;
  if (blen > (kPglFnamesize - extra_size)) {
    logerrprintf("Error: --%s filename too long.\n", flagname_p);
    return kPglRetOpenFail;
  }
  if (pgl_malloc(blen + extra_size, fnbuf_ptr)) {
    return kPglRetNomem;
  }
  memcpy(*fnbuf_ptr, source, blen);
  return kPglRetSuccess;
}

PglErr AllocAndFlatten(const char* const* sources, uint32_t param_ct, uint32_t max_blen, char** flattened_buf_ptr) {
  uintptr_t tot_blen = 1;
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
    const uint32_t cur_blen = 1 + strlen(sources[param_idx]);
    if (cur_blen > max_blen) {
      return kPglRetInvalidCmdline;
    }
    tot_blen += cur_blen;
  }
  char* buf_iter;
  if (pgl_malloc(tot_blen, &buf_iter)) {
    return kPglRetNomem;
  }
  *flattened_buf_ptr = buf_iter;
  for (uint32_t param_idx = 0; param_idx != param_ct; ++param_idx) {
    buf_iter = strcpyax(buf_iter, sources[param_idx], '\0');
  }
  *buf_iter = '\0';
  return kPglRetSuccess;
}

PglErr Rerun(const char* ver_str, const char* ver_str2, const char* prog_name_str, uint32_t rerun_argv_pos, uint32_t rerun_parameter_present, int32_t* argc_ptr, uint32_t* first_arg_idx_ptr, char*** argv_ptr, char*** subst_argv_ptr, char** rerun_buf_ptr) {
  // caller is responsible for freeing rerun_buf

  // ok, requiring zlib/zstd here is totally not worth it
  FILE* rerunfile = nullptr;

  char** subst_argv2 = nullptr;
  char** argv = *argv_ptr;
  char* rerun_fname = rerun_parameter_present? argv[rerun_argv_pos + 1] : g_textbuf;
  uintptr_t line_idx = 1;
  PglErr reterr = kPglRetSuccess;
  {
    if (!rerun_parameter_present) {
      char* write_iter = strcpya(rerun_fname, prog_name_str);
      snprintf(write_iter, kMaxOutfnameExtBlen, ".log");
    }
    rerunfile = fopen(rerun_fname, FOPEN_RB);
    if (unlikely(!rerunfile)) {
      goto Rerun_ret_OPEN_FAIL;
    }
    char* textbuf = g_textbuf;
    textbuf[kMaxMediumLine - 1] = ' ';
    if (unlikely(!fgets(textbuf, kMaxMediumLine, rerunfile))) {
      fputs(ver_str, stdout);
      fputs(ver_str2, stdout);
      fputs("Error: Empty log file for --rerun.\n", stderr);
      goto Rerun_ret_MALFORMED_INPUT;
    }
    if (unlikely(!textbuf[kMaxMediumLine - 1])) {
      goto Rerun_ret_LONG_LINE;
    }
    if (unlikely(!fgets(textbuf, kMaxMediumLine, rerunfile))) {
      fputs(ver_str, stdout);
      fputs(ver_str2, stdout);
      fputs("Error: Only one line in --rerun log file.\n", stderr);
      goto Rerun_ret_MALFORMED_INPUT;
    }
    ++line_idx;
    if (unlikely(!textbuf[kMaxMediumLine - 1])) {
      goto Rerun_ret_LONG_LINE;
    }
    // don't bother supporting "xx arguments: --aa bb --cc --dd" format
    while ((!StrStartsWithUnsafe(textbuf, "Options in effect:")) || (textbuf[strlen("Options in effect:")] >= ' ')) {
      ++line_idx;
      if (unlikely(!fgets(textbuf, kMaxMediumLine, rerunfile))) {
        fputs(ver_str, stdout);
        fputs(ver_str2, stdout);
        fputs("Error: Invalid log file for --rerun.\n", stderr);
        goto Rerun_ret_MALFORMED_INPUT;
      }
    }
    char* all_args_write_iter = textbuf;
    char* textbuf_limit = &(textbuf[kMaxMediumLine]);
    uint32_t loaded_arg_ct = 0;
    // We load each of the option lines in sequence into textbuf, always
    // overwriting the previous line's newline.  (Note that textbuf[] has
    // size > 2 * kMaxMediumLine; this lets us avoid additional
    // dynamic memory allocation as long as we impose the constraint that all
    // lines combined add up to less than kMaxMediumLine.)
    while (1) {
      all_args_write_iter[kMaxMediumLine - 1] = ' ';
      if (!fgets(all_args_write_iter, kMaxMediumLine, rerunfile)) {
        break;
      }
      ++line_idx;
      if (unlikely(!all_args_write_iter[kMaxMediumLine - 1])) {
        goto Rerun_ret_LONG_LINE;
      }
      char* arg_iter = FirstNonTspace(all_args_write_iter);
      if (IsEolnKns(*arg_iter)) {
        *all_args_write_iter = '\0';
        break;
      }
      char* token_end;
      do {
        token_end = CurTokenEnd(arg_iter);
        ++loaded_arg_ct;
        arg_iter = FirstNonTspace(token_end);
      } while (!IsEolnKns(*arg_iter));
      all_args_write_iter = token_end;
      if (unlikely(all_args_write_iter >= textbuf_limit)) {
        fputs(ver_str, stdout);
        fputs(ver_str2, stdout);
        fputs("Error: --rerun argument sequence too long.\n", stderr);
        goto Rerun_ret_MALFORMED_INPUT;
      }
    }
    fclose_null(&rerunfile);
    const uint32_t line_byte_ct = 1 + S_CAST(uintptr_t, all_args_write_iter - textbuf);
    char* rerun_buf;
    if (unlikely(pgl_malloc(line_byte_ct, &rerun_buf))) {
      goto Rerun_ret_NOMEM;
    }
    *rerun_buf_ptr = rerun_buf;
    memcpy(rerun_buf, textbuf, line_byte_ct);
    const uint32_t argc = *argc_ptr;
    const uint32_t first_arg_idx = *first_arg_idx_ptr;
    char* rerun_first_token = FirstNonTspace(rerun_buf);
    const char* arg_iter = rerun_first_token;
    // now use textbuf as a lame bitfield
    memset(textbuf, 1, loaded_arg_ct);
    uint32_t loaded_arg_idx = 0;
    uint32_t duplicate_arg_ct = 0;
    do {
      if (unlikely(NoMoreTokensKns(arg_iter))) {
        fputs(ver_str, stdout);
        fputs(ver_str2, stdout);
        fputs("Error: Line 2 of --rerun log file has fewer tokens than expected.\n", stderr);
        goto Rerun_ret_MALFORMED_INPUT;
      }
      const char* flagname_p = IsCmdlineFlagStart(arg_iter);
      if (flagname_p) {
        const uint32_t slen = strlen_se(flagname_p);
        uint32_t cmdline_arg_idx = first_arg_idx;
        for (; cmdline_arg_idx != argc; ++cmdline_arg_idx) {
          const char* later_flagname_p = IsCmdlineFlagStart(argv[cmdline_arg_idx]);
          if (later_flagname_p) {
            const uint32_t slen2 = strlen(later_flagname_p);
            if ((slen == slen2) && memequal(flagname_p, later_flagname_p, slen)) {
              cmdline_arg_idx = UINT32_MAX;
              break;
            }
          }
        }
        if (cmdline_arg_idx == UINT32_MAX) {
          // matching flag, override --rerun
          do {
            ++duplicate_arg_ct;
            textbuf[loaded_arg_idx++] = 0;
            if (loaded_arg_idx == loaded_arg_ct) {
              break;
            }
            arg_iter = NextToken(arg_iter);
          } while (!IsCmdlineFlag(arg_iter));
        } else {
          ++loaded_arg_idx;
          arg_iter = NextToken(arg_iter);
        }
      } else {
        ++loaded_arg_idx;
        arg_iter = NextToken(arg_iter);
      }
    } while (loaded_arg_idx < loaded_arg_ct);
    if (unlikely(pgl_malloc((argc + loaded_arg_ct - duplicate_arg_ct - rerun_parameter_present - 1 - first_arg_idx) * sizeof(intptr_t), &subst_argv2))) {
      goto Rerun_ret_NOMEM;
    }
    uint32_t new_arg_idx = rerun_argv_pos - first_arg_idx;
    memcpy(subst_argv2, &(argv[first_arg_idx]), new_arg_idx * sizeof(intptr_t));
    char* arg_nullterminate_iter = rerun_first_token;
    for (loaded_arg_idx = 0; loaded_arg_idx != loaded_arg_ct; ++loaded_arg_idx) {
      arg_nullterminate_iter = FirstNonTspace(arg_nullterminate_iter);
      char* token_end = CurTokenEnd(arg_nullterminate_iter);
      if (textbuf[loaded_arg_idx]) {
        subst_argv2[new_arg_idx++] = arg_nullterminate_iter;
        *token_end = '\0';
      }
      arg_nullterminate_iter = &(token_end[1]);
    }
    const uint32_t final_copy_start_idx = rerun_argv_pos + rerun_parameter_present + 1;
    memcpy(&(subst_argv2[new_arg_idx]), &(argv[final_copy_start_idx]), (argc - final_copy_start_idx) * sizeof(intptr_t));
    *first_arg_idx_ptr = 0;
    *argc_ptr = new_arg_idx + argc - final_copy_start_idx;
    if (*subst_argv_ptr) {
      free(*subst_argv_ptr);
    }
    *subst_argv_ptr = subst_argv2;
    *argv_ptr = subst_argv2;
    subst_argv2 = nullptr;
  }
  while (0) {
  Rerun_ret_NOMEM:
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
    reterr = kPglRetNomem;
    break;
  Rerun_ret_OPEN_FAIL:
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
    fprintf(stderr, kErrprintfFopen, rerun_fname, strerror(errno));
    reterr = kPglRetOpenFail;
    break;
  Rerun_ret_LONG_LINE:
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
    fprintf(stderr, "Error: Line %" PRIuPTR " of --rerun log file is pathologically long.\n", line_idx);
  Rerun_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  free_cond(subst_argv2);
  fclose_cond(rerunfile);
  return reterr;
}


// Handles --script, --rerun, --help, --version, and --silent.
// subst_argv, script_buf, and rerun_buf must be initialized to nullptr.
PglErr CmdlineParsePhase1(const char* ver_str, const char* ver_str2, const char* prog_name_str, const char* notestr_null_calc2, const char* cmdline_format_str, const char* errstr_append, uint32_t max_flag_blen, PglErr(* disp_help_fn)(const char* const*, uint32_t), int* argc_ptr, char*** argv_ptr, Plink2CmdlineMeta* pcmp, uint32_t* first_arg_idx_ptr, uint32_t* flag_ct_ptr) {
  FILE* scriptfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    int argc = *argc_ptr;
    const char* const* argvk = TO_CONSTCPCONSTP(*argv_ptr);
    char** subst_argv = nullptr;
    uint32_t first_arg_idx = 1;
    for (uint32_t arg_idx = 1; arg_idx != S_CAST(uint32_t, argc); ++arg_idx) {
      if ((!strcmp("-script", argvk[arg_idx])) || (!strcmp("--script", argvk[arg_idx]))) {
        const uint32_t param_ct = GetParamCt(argvk, argc, arg_idx);
        if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fputs(g_logbuf, stderr);
          fputs(errstr_append, stderr);
          goto CmdlineParsePhase1_ret_INVALID_CMDLINE;
        }
        for (uint32_t arg_idx2 = arg_idx + 2; arg_idx2 != S_CAST(uint32_t, argc); ++arg_idx2) {
          if (unlikely((!strcmp("-script", argvk[arg_idx2])) || (!strcmp("--script", argvk[arg_idx2])))) {
            fputs(ver_str, stdout);
            fputs(ver_str2, stdout);
            fputs("Error: Multiple --script flags.  Merge the files into one.\n", stderr);
            fputs(errstr_append, stderr);
            goto CmdlineParsePhase1_ret_INVALID_CMDLINE;
          }
        }
        // logging not yet active, so don't use fopen_checked()
        scriptfile = fopen(argvk[arg_idx + 1], FOPEN_RB);
        if (unlikely(!scriptfile)) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fprintf(stderr, kErrprintfFopen, argvk[arg_idx + 1], strerror(errno));
          goto CmdlineParsePhase1_ret_OPEN_FAIL;
        }
        if (unlikely(fseeko(scriptfile, 0, SEEK_END))) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fprintf(stderr, kErrprintfFread, argvk[arg_idx + 1], rstrerror(errno));
          goto CmdlineParsePhase1_ret_READ_FAIL;
        }
        int64_t fsize = ftello(scriptfile);
        if (unlikely(fsize < 0)) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fprintf(stderr, kErrprintfFread, argvk[arg_idx + 1], rstrerror(errno));
          goto CmdlineParsePhase1_ret_READ_FAIL;
        }
        if (unlikely(fsize > 0x7ffffffe)) {
          // could actually happen if user enters parameters in the wrong
          // order, so may as well catch it and print a somewhat informative
          // error message
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fputs("Error: --script file too large.", stderr);
          goto CmdlineParsePhase1_ret_INVALID_CMDLINE;
        }
        rewind(scriptfile);
        const uint32_t fsize_ui = fsize;
        if (unlikely(pgl_malloc(fsize_ui + 1, &pcmp->script_buf))) {
          goto CmdlineParsePhase1_ret_NOMEM;
        }
        char* script_buf = pcmp->script_buf;
        if (unlikely(!fread_unlocked(script_buf, fsize_ui, 1, scriptfile))) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          if (feof_unlocked(scriptfile)) {
            // shouldn't actually be possible, but let's handle this in as
            // uniform a manner as we can
            errno = 0;
          }
          fprintf(stderr, kErrprintfFread, argvk[arg_idx + 1], rstrerror(errno));
          goto CmdlineParsePhase1_ret_READ_FAIL;
        }
        script_buf[fsize_ui] = '\0';
        fclose_null(&scriptfile);
        uint32_t num_script_params = 0;
        char* script_buf_iter = script_buf;
        uint32_t char_code;
        do {
          uint32_t char_code_m1;
          do {
            char_code_m1 = ctou32(*script_buf_iter++) - 1;
          } while (char_code_m1 < 32);
          if (char_code_m1 == UINT32_MAX) {
            break;
          }
          ++num_script_params;
          do {
            char_code = ctou32(*script_buf_iter++);
          } while (char_code > 32);
        } while (char_code);
        if (unlikely(script_buf_iter != (&(script_buf[fsize_ui + 1])))) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fputs("Error: Null byte in --script file.\n", stderr);
          goto CmdlineParsePhase1_ret_INVALID_CMDLINE;
        }
        // probable todo: detect duplicate flags in the same manner as --rerun
        const uint32_t new_param_ct = num_script_params + argc - 3;
        if (unlikely(pgl_malloc(new_param_ct * sizeof(intptr_t), &subst_argv))) {
          goto CmdlineParsePhase1_ret_NOMEM;
        }
        pcmp->subst_argv = subst_argv;
        memcpy(subst_argv, &(argvk[1]), arg_idx * sizeof(intptr_t));
        const uint32_t load_param_idx_end = arg_idx + num_script_params;
        script_buf_iter = &(script_buf[-1]);
        for (uint32_t param_idx = arg_idx; param_idx != load_param_idx_end; ++param_idx) {
          while (ctou32(*(++script_buf_iter)) <= 32);
          subst_argv[param_idx] = script_buf_iter;
          while (ctou32(*(++script_buf_iter)) > 32);
          // could enforce some sort of length limit here
          *script_buf_iter = '\0';
        }
        memcpy(&(subst_argv[load_param_idx_end]), &(argvk[arg_idx + 2]), (argc - arg_idx - 2) * sizeof(intptr_t));
        argc = new_param_ct;
        *argc_ptr = argc;
        first_arg_idx = 0;
        argvk = TO_CONSTCPCONSTP(subst_argv);
        *argv_ptr = subst_argv;
        break;
      }
    }
    for (uint32_t arg_idx = first_arg_idx; arg_idx != S_CAST(uint32_t, argc); ++arg_idx) {
      if ((!strcmp("-rerun", argvk[arg_idx])) || (!strcmp("--rerun", argvk[arg_idx]))) {
        const uint32_t param_ct = GetParamCt(argvk, argc, arg_idx);
        if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 0, 1))) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          fputs(g_logbuf, stderr);
          fputs(errstr_append, stderr);
          goto CmdlineParsePhase1_ret_INVALID_CMDLINE;
        }
        for (uint32_t arg_idx2 = arg_idx + param_ct + 1; arg_idx2 != S_CAST(uint32_t, argc); ++arg_idx2) {
          if (unlikely((!strcmp("-rerun", argvk[arg_idx2])) || (!strcmp("--rerun", argvk[arg_idx2])))) {
            fputs(ver_str, stdout);
            fputs(ver_str2, stdout);
            fputs("Error: Duplicate --rerun flag.\n", stderr);
            goto CmdlineParsePhase1_ret_INVALID_CMDLINE;
          }
        }
        reterr = Rerun(ver_str, ver_str2, prog_name_str, arg_idx, param_ct, &argc, &first_arg_idx, argv_ptr, &pcmp->subst_argv, &pcmp->rerun_buf);
        if (unlikely(reterr)) {
          goto CmdlineParsePhase1_ret_1;
        }
        *argc_ptr = argc;
        subst_argv = pcmp->subst_argv;
        argvk = TO_CONSTCPCONSTP(*argv_ptr);
        break;
      }
    }
    if (unlikely((first_arg_idx < S_CAST(uint32_t, argc)) && (!IsCmdlineFlag(argvk[first_arg_idx])))) {
      fputs("Error: First argument must be a flag.\n", stderr);
      fputs(errstr_append, stderr);
      goto CmdlineParsePhase1_ret_INVALID_CMDLINE;
    }
    uint32_t flag_ct = 0;
    uint32_t version_present = 0;
    uint32_t silent_present = 0;
    for (uint32_t arg_idx = first_arg_idx; arg_idx != S_CAST(uint32_t, argc); ++arg_idx) {
      const char* flagname_p = IsCmdlineFlagStart(argvk[arg_idx]);
      if (flagname_p) {
        const uint32_t flagname_p_slen = strlen(flagname_p);
        if (strequal_k(flagname_p, "help", flagname_p_slen)) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          if ((!first_arg_idx) || (arg_idx != 1) || subst_argv) {
            fputs("--help present, ignoring other flags.\n", stdout);
          }
          if ((arg_idx == S_CAST(uint32_t, argc) - 1) && flag_ct) {
            // make "plink [valid flags/arguments] --help" work, and skip the
            // arguments
            const char** help_argv;
            if (unlikely(pgl_malloc(flag_ct * sizeof(intptr_t), &help_argv))) {
              goto CmdlineParsePhase1_ret_NOMEM2;
            }
            uint32_t arg_idx2 = 0;
            for (uint32_t flag_idx = 0; flag_idx != flag_ct; ++flag_idx) {
              while (!IsCmdlineFlagStart(argvk[++arg_idx2]));
              help_argv[flag_idx] = argvk[arg_idx2];
            }
            reterr = disp_help_fn(help_argv, flag_ct);
            free(help_argv);
          } else {
            reterr = disp_help_fn(&(argvk[arg_idx + 1]), argc - arg_idx - 1);
          }
          if (!reterr) {
            reterr = kPglRetHelp;
          }
          goto CmdlineParsePhase1_ret_1;
        }
        if (strequal_k(flagname_p, "h", flagname_p_slen) ||
            strequal_k(flagname_p, "?", flagname_p_slen)) {
          // these just act like the no-argument case
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          if ((!first_arg_idx) || (arg_idx != 1) || subst_argv) {
            printf("-%c present, ignoring other flags.\n", *flagname_p);
          }
          fputs(cmdline_format_str, stdout);
          fputs(notestr_null_calc2, stdout);
          reterr = kPglRetHelp;
          goto CmdlineParsePhase1_ret_1;
        }
        if (strequal_k(flagname_p, "version", flagname_p_slen)) {
          version_present = 1;
        } else if (strequal_k(flagname_p, "silent", flagname_p_slen)) {
          silent_present = 1;
        } else if (unlikely(flagname_p_slen >= max_flag_blen)) {
          fputs(ver_str, stdout);
          fputs(ver_str2, stdout);
          // shouldn't be possible for this to overflow the buffer...
          snprintf(g_logbuf, kLogbufSize, "Error: Unrecognized flag ('%s').\n", argvk[arg_idx]);
          WordWrapB(0);
          fputs(g_logbuf, stderr);
          fputs(errstr_append, stderr);
          goto CmdlineParsePhase1_ret_INVALID_CMDLINE;
        }
        ++flag_ct;
      }
    }
    if (version_present) {
      fputs(ver_str, stdout);
      putc_unlocked('\n', stdout);
      reterr = kPglRetHelp;
      goto CmdlineParsePhase1_ret_1;
    }
    if (silent_present) {
      if (unlikely(!freopen(NULL_STREAM_NAME, "w", stdout))) {
        fputs("Warning: --silent failed.", stderr);
        g_stderr_written_to = 1;
      }
    }
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
    *first_arg_idx_ptr = first_arg_idx;
    *flag_ct_ptr = flag_ct;
  }
  while (0) {
  CmdlineParsePhase1_ret_NOMEM:
    fputs(ver_str, stdout);
    fputs(ver_str2, stdout);
  CmdlineParsePhase1_ret_NOMEM2:
    fputs(kErrstrNomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
    }
    reterr = kPglRetNomem;
    break;
  CmdlineParsePhase1_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  CmdlineParsePhase1_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  CmdlineParsePhase1_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
 CmdlineParsePhase1_ret_1:
  fclose_cond(scriptfile);
  return reterr;
}

// Assumes CmdlineParsePhase1() has completed, flag names have been copied to
// flag_buf/flag_map, aliases handled, and PROG_NAME_STR has been copied to
// outname (no null-terminator needed).  outname_end must be initialized to
// nullptr.
// This sorts the flag names so they're processed in a predictable order,
// handles --out if present, initializes the log, and determines the number of
// processors the OS wants us to think the machine has.
PglErr CmdlineParsePhase2(const char* ver_str, const char* errstr_append, const char* const* argvk, uint32_t prog_name_str_slen, uint32_t max_flag_blen, int32_t argc, uint32_t flag_ct, Plink2CmdlineMeta* pcmp, char* outname, char** outname_end_ptr, char* range_delim_ptr, int32_t* known_procs_ptr, uint32_t* max_thread_ct_ptr) {
  PglErr reterr = kPglRetSuccess;
  {
    char* flag_buf = pcmp->flag_buf;
    uint32_t* flag_map = pcmp->flag_map;
    reterr = SortCmdlineFlags(max_flag_blen, flag_ct, flag_buf, flag_map);
    if (unlikely(reterr)) {
      if (reterr == kPglRetNomem) {
        goto CmdlineParsePhase2_ret_NOMEM_NOLOG;
      }
      goto CmdlineParsePhase2_ret_1;
    }

    *range_delim_ptr = '-';
    for (uint32_t cur_flag_idx = 0; cur_flag_idx != flag_ct; ++cur_flag_idx) {
      char* cur_flag = &(flag_buf[cur_flag_idx * max_flag_blen]);
      if (strequal_k_unsafe(cur_flag, "d")) {
        // Must be here, instead of in the main parse loop, for --covar-name +
        // --d to work.
        const uint32_t arg_idx = flag_map[cur_flag_idx];
        const uint32_t param_ct = GetParamCt(argvk, argc, arg_idx);
        if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
          fputs(g_logbuf, stderr);
          fputs(errstr_append, stderr);
          goto CmdlineParsePhase2_ret_INVALID_CMDLINE;
        }
        const char cc = ExtractCharParam(argvk[arg_idx + 1]);
        if (unlikely(!cc)) {
          fputs("Error: --d argument too long (must be a single character).\n", stderr);
          goto CmdlineParsePhase2_ret_INVALID_CMDLINE;
        }
        if ((cc == '-') || (cc == ',')) {
          fputs("Error: --d argument cannot be '-' or ','.\n", stderr);
          goto CmdlineParsePhase2_ret_INVALID_CMDLINE;
        }
        *range_delim_ptr = cc;
        // bugfix (31 Jul 2019): can't set *cur_flag = '\0' to mark the flag as
        // already-processed, since we haven't written the start of the .log
        // file yet.
        continue;
      }
      const int32_t memcmp_out_result = Memcmp("out", cur_flag, 4);
      if (!memcmp_out_result) {
        const uint32_t arg_idx = flag_map[cur_flag_idx];
        const uint32_t param_ct = GetParamCt(argvk, argc, arg_idx);
        if (unlikely(EnforceParamCtRange(argvk[arg_idx], param_ct, 1, 1))) {
          fputs(g_logbuf, stderr);
          fputs(errstr_append, stderr);
          goto CmdlineParsePhase2_ret_INVALID_CMDLINE;
        }
        if (unlikely(strlen(argvk[arg_idx + 1]) > (kPglFnamesize - kMaxOutfnameExtBlen))) {
          fflush(stdout);
          fputs("Error: --out argument too long.\n", stderr);
          goto CmdlineParsePhase2_ret_OPEN_FAIL;
        }
        const uint32_t slen = strlen(argvk[arg_idx + 1]);
        memcpy(outname, argvk[arg_idx + 1], slen + 1);
        *outname_end_ptr = &(outname[slen]);
      }
      if (memcmp_out_result <= 0) {
        break;
      }
    }
    if (unlikely(InitLogfile(0, outname, (*outname_end_ptr)? (*outname_end_ptr) : &(outname[prog_name_str_slen])))) {
      goto CmdlineParsePhase2_ret_OPEN_FAIL;
    }
    logputs_silent(ver_str);
    logputs_silent("\n");
    logputs("Options in effect:\n");
    for (uint32_t cur_flag_idx = 0; cur_flag_idx != flag_ct; ++cur_flag_idx) {
      logputs("  --");
      logputs(&(flag_buf[cur_flag_idx * max_flag_blen]));
      uint32_t arg_idx = flag_map[cur_flag_idx] + 1;
      while ((arg_idx < S_CAST(uint32_t, argc)) && (!IsCmdlineFlag(argvk[arg_idx]))) {
        logputs(" ");
        // Thought about special-casing argvk[arg_idx] == " ", so that
        // "--id-delim ' '" actually works with --rerun, but that adds too much
        // complexity to the rereader to be worth it.  Instead we just document
        // the incompatibility.
        logputs(argvk[arg_idx++]);
      }
      logputs("\n");
    }
    logputs("\n");

#ifdef _WIN32
    DWORD windows_dw = kTextbufSize;
    if (GetComputerName(g_textbuf, &windows_dw))
#else
    if (gethostname(g_textbuf, kTextbufSize) != -1)
#endif
    {
      logputs_silent("Hostname: ");
      logputs_silent(g_textbuf);
    }
    logputs_silent("\nWorking directory: ");
    if (unlikely(!getcwd(g_textbuf, kPglFnamesize))) {
      logputs_silent("\n");
      logerrprintfww("Error: Failed to get current working directory: %s.\n", strerror(errno));
      // debatable what error code applies here
      goto CmdlineParsePhase2_ret_OPEN_FAIL;
    }
    logputs_silent(g_textbuf);
    logputs_silent("\n");
    logputs("Start time: ");
    time_t rawtime;
    time(&rawtime);
    logputs(ctime(&rawtime));
    // ctime string always has a newline at the end
    logputs_silent("\n");

    *max_thread_ct_ptr = NumCpu(known_procs_ptr);
    // don't subtract 1 any more since, when max_thread_ct > 2, one of the
    // (virtual) cores will be dedicated to I/O and have lots of downtime.
    //
    // may make kMaxThreads a parameter later.
    if (*max_thread_ct_ptr > kMaxThreads) {
      *max_thread_ct_ptr = kMaxThreads;
    }
  }
  while (0) {
  CmdlineParsePhase2_ret_NOMEM_NOLOG:
    fputs(kErrstrNomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIuPTR "\n", g_failed_alloc_attempt_size);
    }
    reterr = kPglRetNomem;
    break;
  CmdlineParsePhase2_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  CmdlineParsePhase2_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
 CmdlineParsePhase2_ret_1:
  return reterr;
}

PglErr CmdlineParsePhase3(uintptr_t max_default_mib, uintptr_t malloc_size_mib, uint32_t memory_require, Plink2CmdlineMeta* pcmp, unsigned char** bigstack_ua_ptr) {
  PglErr reterr = kPglRetSuccess;
  {
    if (pcmp->subst_argv) {
      free(pcmp->subst_argv);
      pcmp->subst_argv = nullptr;
    }
    if (pcmp->script_buf) {
      free(pcmp->script_buf);
      pcmp->script_buf = nullptr;
    }
    if (pcmp->rerun_buf) {
      free(pcmp->rerun_buf);
      pcmp->rerun_buf = nullptr;
    }
    if (pcmp->flag_buf) {
      free(pcmp->flag_buf);
      pcmp->flag_buf = nullptr;
    }
    if (pcmp->flag_map) {
      free(pcmp->flag_map);
      pcmp->flag_map = nullptr;
    }

    uint64_t total_mib = DetectMib();
    if (!malloc_size_mib) {
      if (!total_mib) {
        malloc_size_mib = max_default_mib? max_default_mib : kBigstackDefaultMib;
      } else if (total_mib < (kBigstackMinMib * 2)) {
        malloc_size_mib = kBigstackMinMib;
      } else {
        malloc_size_mib = total_mib / 2;
        if (max_default_mib && (malloc_size_mib > max_default_mib)) {
          malloc_size_mib = max_default_mib;
        }
      }
    }
    assert(malloc_size_mib >= kBigstackMinMib);
#ifndef __LP64__
    if (malloc_size_mib > kMalloc32bitMibMax) {
      malloc_size_mib = kMalloc32bitMibMax;
    }
#endif
    if (total_mib) {
      snprintf(g_logbuf, kLogbufSize, "%" PRIu64 " MiB RAM detected; reserving %" PRIuPTR " MiB for main workspace.\n", total_mib, malloc_size_mib);
    } else {
      snprintf(g_logbuf, kLogbufSize, "Failed to determine total system memory.  Attempting to reserve %" PRIuPTR " MiB.\n", malloc_size_mib);
    }
    logputsb();
    uintptr_t malloc_mib_final;
    if (unlikely(InitBigstack(malloc_size_mib, &malloc_mib_final, bigstack_ua_ptr))) {
      goto CmdlineParsePhase3_ret_NOMEM;
    }
    if (malloc_size_mib != malloc_mib_final) {
      if (unlikely(memory_require)) {
        goto CmdlineParsePhase3_ret_NOMEM;
      }
      logprintf("Allocated %" PRIuPTR " MiB successfully, after larger attempt(s) failed.\n", malloc_mib_final);
    }
  }
  while (0) {
  CmdlineParsePhase3_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
  return reterr;
}

void CleanupPlink2CmdlineMeta(Plink2CmdlineMeta* pcmp) {
  free_cond(pcmp->subst_argv);
  free_cond(pcmp->script_buf);
  free_cond(pcmp->rerun_buf);
  free_cond(pcmp->flag_buf);
  free_cond(pcmp->flag_map);
}


void InitCmpExpr(CmpExpr* cmp_expr_ptr) {
  cmp_expr_ptr->pheno_name = nullptr;
}

void CleanupCmpExpr(CmpExpr* cmp_expr_ptr) {
  free_cond(cmp_expr_ptr->pheno_name);
}

// may want CXXCONST_CP treatment later
const char* ParseNextBinaryOp(const char* expr_str, uint32_t expr_slen, const char** op_start_ptr, CmpBinaryOp* binary_op_ptr) {
  // !=, <>: kCmpOperatorNoteq
  // <: kCmpOperatorLe
  // <=: kCmpOperatorLeq
  // =, ==: kCmpOperatorEq
  // >=: kCmpOperatorGeq
  // >: kCmpOperatorGe
  const char* next_eq = S_CAST(const char*, memchr(expr_str, '=', expr_slen));
  const char* next_lt = S_CAST(const char*, memchr(expr_str, '<', expr_slen));
  const char* next_gt = S_CAST(const char*, memchr(expr_str, '>', expr_slen));
  if (!next_eq) {
    if (!next_lt) {
      // may want to remove unlikely() later
      if (unlikely(!next_gt)) {
        return nullptr;
      }
      *op_start_ptr = next_gt;
      *binary_op_ptr = kCmpOperatorGe;
      return &(next_gt[1]);
    }
    if (next_gt == (&(next_lt[1]))) {
      *op_start_ptr = next_lt;
      *binary_op_ptr = kCmpOperatorNoteq;
      return &(next_lt[2]);
    }
    if ((!next_gt) || (next_gt > next_lt)) {
      *op_start_ptr = next_lt;
      *binary_op_ptr = kCmpOperatorLe;
      return &(next_lt[1]);
    }
    *op_start_ptr = next_gt;
    *binary_op_ptr = kCmpOperatorGe;
    return &(next_gt[1]);
  }
  if ((!next_lt) || (next_lt > next_eq)) {
    if ((!next_gt) || (next_gt > next_eq)) {
      if ((next_eq != expr_str) && (next_eq[-1] == '!')) {
        *op_start_ptr = &(next_eq[-1]);
        *binary_op_ptr = kCmpOperatorNoteq;
        return &(next_eq[1]);
      }
      *op_start_ptr = next_eq;
      *binary_op_ptr = kCmpOperatorEq;
      return (next_eq[1] == '=')? (&(next_eq[2])) : (&(next_eq[1]));
    }
    *op_start_ptr = next_gt;
    if (next_eq == (&(next_gt[1]))) {
      *binary_op_ptr = kCmpOperatorGeq;
      return &(next_gt[2]);
    }
    *binary_op_ptr = kCmpOperatorGe;
    return &(next_gt[1]);
  }
  if (next_gt == (&(next_lt[1]))) {
    *op_start_ptr = next_lt;
    *binary_op_ptr = kCmpOperatorNoteq;
    return &(next_lt[2]);
  }
  if ((!next_gt) || (next_gt > next_lt)) {
    *op_start_ptr = next_lt;
    if (next_eq == (&(next_lt[1]))) {
      *binary_op_ptr = kCmpOperatorLeq;
      return &(next_lt[2]);
    }
    *binary_op_ptr = kCmpOperatorLe;
    return &(next_lt[1]);
  }
  *op_start_ptr = next_gt;
  if (next_eq == (&(next_gt[1]))) {
    *binary_op_ptr = kCmpOperatorGeq;
    return &(next_gt[2]);
  }
  *binary_op_ptr = kCmpOperatorGe;
  return &(next_gt[1]);
}

PglErr ValidateAndAllocCmpExpr(const char* const* sources, const char* flag_name, uint32_t param_ct, CmpExpr* cmp_expr_ptr) {
  // Currently four use cases:
  //   [pheno/covar name] [operator] [pheno val]: regular comparison
  //   [pheno/covar name]: existence check
  //   [INFO key] [operator] [val]: regular comparison
  //   [INFO key]: existence check
  // Some key/value validation is deferred to LoadPvar()/KeepRemoveIf(),
  // since the requirements are different (e.g. no semicolons in anything
  // INFO-related, categorical phenotypes can be assumed to not start with a
  // valid number).
  // May support or/and, parentheses later, but need to be careful to not slow
  // down LoadPvar() too much in the no-INFO-filter case.
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely((param_ct != 1) && (param_ct != 3))) {
      goto ValidateAndAllocCmpExpr_ret_INVALID_EXPR_GENERIC;
    }
    const char* pheno_name_start = sources[0];
    const char* pheno_val_start;
    uint32_t pheno_name_slen;
    uint32_t pheno_val_slen;
    if (param_ct == 3) {
      pheno_name_slen = strlen(pheno_name_start);
      const char* op_str = sources[1];
      uint32_t op_slen = strlen(op_str);
      // ok to have single/double quotes around operator
      if (op_slen > 2) {
        const char cc = op_str[0];
        if (((cc == '\'') || (cc == '"')) && (op_str[op_slen - 1] == cc)) {
          ++op_str;
          op_slen -= 2;
        }
      }
      const char* op_start;
      const char* op_end = ParseNextBinaryOp(op_str, op_slen, &op_start, &cmp_expr_ptr->binary_op);
      if (unlikely((!op_end) || (*op_end) || (op_start != op_str))) {
        goto ValidateAndAllocCmpExpr_ret_INVALID_EXPR_GENERIC;
      }
      pheno_val_start = sources[2];
      pheno_val_slen = strlen(pheno_val_start);
    } else {
      // permit param_ct == 1 as long as tokens are unambiguous
      uint32_t expr_slen = strlen(pheno_name_start);
      const char* op_start;
      pheno_val_start = ParseNextBinaryOp(pheno_name_start, expr_slen, &op_start, &cmp_expr_ptr->binary_op);
      if (unlikely((!pheno_val_start) || (!(*pheno_val_start)) || (op_start == pheno_name_start))) {
        goto ValidateAndAllocCmpExpr_ret_INVALID_EXPR_GENERIC;
      }
      pheno_name_slen = op_start - pheno_name_start;
      // quasi-bugfix (13 Dec 2017): allow single argument to contain internal
      // spaces;
      //   --keep-if "PHENO1 > 1"
      // is more intuitive and consistent with usage of other command-line
      // tools than
      //   --keep-if PHENO1 '>' 1
      //
      // To prevent --rerun from breaking, if there's a space after the
      // operator, there must be a space before the operator as well, etc.
      if (*pheno_val_start == ' ') {
        if (unlikely(pheno_name_start[pheno_name_slen - 1] != ' ')) {
          goto ValidateAndAllocCmpExpr_ret_INVALID_EXPR_GENERIC;
        }
        do {
          ++pheno_val_start;
        } while (*pheno_val_start == ' ');
        do {
          --pheno_name_slen;
          if (unlikely(!pheno_name_slen)) {
            goto ValidateAndAllocCmpExpr_ret_INVALID_EXPR_GENERIC;
          }
        } while (pheno_name_start[pheno_name_slen - 1] == ' ');
      }
      pheno_val_slen = expr_slen - S_CAST(uintptr_t, pheno_val_start - pheno_name_start);
    }
    if (unlikely(memchr(pheno_name_start, ' ', pheno_name_slen) || memchr(pheno_val_start, ' ', pheno_val_slen))) {
      goto ValidateAndAllocCmpExpr_ret_INVALID_EXPR_GENERIC;
    }
    if (unlikely((pheno_name_slen > kMaxIdSlen) || (pheno_val_slen > kMaxIdSlen))) {
      logerrprintf("Error: ID too long in %s expression.\n", flag_name);
      goto ValidateAndAllocCmpExpr_ret_INVALID_CMDLINE;
    }
    if ((cmp_expr_ptr->binary_op != kCmpOperatorNoteq) && (cmp_expr_ptr->binary_op != kCmpOperatorEq)) {
      double dxx;
      if (unlikely(!ScantokDouble(pheno_val_start, &dxx))) {
        logerrprintfww("Error: Invalid %s value '%s' (finite number expected).\n", flag_name, pheno_val_start);
        goto ValidateAndAllocCmpExpr_ret_INVALID_CMDLINE;
      }
    }
    char* new_pheno_name_buf;
    if (unlikely(pgl_malloc(2 + pheno_name_slen + pheno_val_slen, &new_pheno_name_buf))) {
      goto ValidateAndAllocCmpExpr_ret_NOMEM;
    }
    memcpyx(new_pheno_name_buf, pheno_name_start, pheno_name_slen, '\0');
    // pheno_val_start guaranteed to be null-terminated for now
    memcpy(&(new_pheno_name_buf[pheno_name_slen + 1]), pheno_val_start, pheno_val_slen + 1);
    cmp_expr_ptr->pheno_name = new_pheno_name_buf;
  }
  while (0) {
  ValidateAndAllocCmpExpr_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ValidateAndAllocCmpExpr_ret_INVALID_EXPR_GENERIC:
    logerrprintf("Error: Invalid %s expression.\n", flag_name);
  ValidateAndAllocCmpExpr_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
  return reterr;
}

// See e.g. meta_analysis_open_and_read_header() in plink 1.9.
// Assumes at least one search term.
PglErr SearchHeaderLine(const char* header_line_iter, const char* const* search_multistrs, const char* flagname_p, uint32_t search_col_ct, uint32_t* found_col_ct_ptr, uint32_t* found_type_bitset_ptr, uint32_t* col_skips, uint32_t* col_types) {
  assert(search_col_ct <= 32);
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    uint32_t search_term_ct = 0;
    uintptr_t max_blen = 0;
    for (uint32_t search_col_idx = 0; search_col_idx != search_col_ct; ++search_col_idx) {
      const uint32_t cur_search_term_ct = CountAndMeasureMultistr(search_multistrs[search_col_idx], &max_blen);
      assert(cur_search_term_ct <= (1 << 26));
      search_term_ct += cur_search_term_ct;
    }
    char* merged_strbox;
    uint32_t* id_map;
    uint32_t* priority_vals;
    uint64_t* cols_and_types;
    if (unlikely(
            bigstack_alloc_c(search_term_ct * max_blen, &merged_strbox) ||
            bigstack_alloc_u32(search_term_ct, &id_map) ||
            bigstack_alloc_u32(search_col_ct, &priority_vals) ||
            bigstack_alloc_u64(search_col_ct, &cols_and_types))) {
      goto SearchHeaderLine_ret_NOMEM;
    }
    uint32_t search_term_idx = 0;
    for (uint32_t search_col_idx = 0; search_col_idx != search_col_ct; ++search_col_idx) {
      const char* multistr_iter = search_multistrs[search_col_idx];
      uint32_t priority_idx = 0;  // 0 = highest priority
      while (*multistr_iter) {
        const uint32_t cur_blen = strlen(multistr_iter) + 1;
        memcpy(&(merged_strbox[search_term_idx * max_blen]), multistr_iter, cur_blen);
        id_map[search_term_idx] = search_col_idx + (priority_idx * 32);
        ++search_term_idx;
        ++priority_idx;
        multistr_iter = &(multistr_iter[cur_blen]);
      }
    }
    SortStrboxIndexed(search_term_ct, max_blen, 0, merged_strbox, id_map);
    assert(search_term_ct);
    const char* duplicate_search_term = ScanForDuplicateIds(merged_strbox, search_term_ct, max_blen);
    if (unlikely(duplicate_search_term)) {
      logerrprintfww("Error: Duplicate term '%s' in --%s column search order.\n", duplicate_search_term, flagname_p);
      goto SearchHeaderLine_ret_INVALID_CMDLINE;
    }

    SetAllU32Arr(search_col_ct, priority_vals);
    SetAllU64Arr(search_col_ct, cols_and_types);
    for (uintptr_t col_idx = 0; ; ++col_idx) {
      const char* token_end = CurTokenEnd(header_line_iter);
      const uint32_t token_slen = token_end - header_line_iter;
      int32_t ii = bsearch_str(header_line_iter, merged_strbox, token_slen, max_blen, search_term_ct);
      if (ii != -1) {
        const uint32_t cur_map_idx = id_map[S_CAST(uint32_t, ii)];
        const uint32_t search_col_idx = cur_map_idx & 31;
        const uint32_t priority_idx = cur_map_idx >> 5;
        if (priority_vals[search_col_idx] >= priority_idx) {
          if (unlikely(priority_vals[search_col_idx] == priority_idx)) {
            logerrprintfww("Error: Duplicate column header '%s' in --%s file.\n", &(merged_strbox[max_blen * cur_map_idx]), flagname_p);
            goto SearchHeaderLine_ret_MALFORMED_INPUT;
          }
          priority_vals[search_col_idx] = priority_idx;
          cols_and_types[search_col_idx] = (S_CAST(uint64_t, col_idx) << 32) | search_col_idx;
        }
      }
      header_line_iter = FirstNonTspace(token_end);
      if (IsEolnKns(*header_line_iter)) {
        break;
      }
    }
    uint32_t found_type_bitset = 0;
    for (uint32_t search_col_idx = 0; search_col_idx != search_col_ct; ++search_col_idx) {
      if (priority_vals[search_col_idx] != UINT32_MAX) {
        found_type_bitset |= 1U << search_col_idx;
      }
    }
    const uint32_t found_col_ct = PopcountWord(found_type_bitset);
    *found_col_ct_ptr = found_col_ct;
    *found_type_bitset_ptr = found_type_bitset;
    if (found_col_ct) {
      STD_SORT(search_col_ct, u64cmp, cols_and_types);
      uint32_t prev_col_idx = cols_and_types[0] >> 32;
      col_skips[0] = prev_col_idx;
      col_types[0] = S_CAST(uint32_t, cols_and_types[0]);
      for (uint32_t found_col_idx = 1; found_col_idx != found_col_ct; ++found_col_idx) {
        const uint64_t cur_col_and_type = cols_and_types[found_col_idx];
        const uint32_t cur_col_idx = cur_col_and_type >> 32;
        col_skips[found_col_idx] = cur_col_idx - prev_col_idx;
        col_types[found_col_idx] = S_CAST(uint32_t, cur_col_and_type);
        prev_col_idx = cur_col_idx;
      }
    }
  }
  while (0) {
  SearchHeaderLine_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SearchHeaderLine_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  SearchHeaderLine_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr ParseColDescriptor(const char* col_descriptor_iter, const char* supported_ids, const char* cur_flag_name, uint32_t first_col_shifted, uint32_t default_cols_mask, uint32_t prohibit_empty, void* result_ptr) {
  // col_descriptor is usually a pointer to argv[...][5] (first five characters
  // are "cols=").  supported_ids is a multistr.
  // may need to switch first_col_shifted/default_cols_mask/result to uint64_t
  PglErr reterr = kPglRetSuccess;
  uint32_t* id_map = nullptr;
  {
    uint32_t max_id_blen = 0;
    uint32_t id_ct = 0;

    // work around strchr not returning a const char*?
    const char* supported_ids_iter = supported_ids;

    // can precompute this sorted index and avoid the dynamic
    // allocations/deallocations, but this is cheap enough that I'd rather make
    // it easier to extend functionality.
    do {
      const uint32_t blen = 1 + strlen(supported_ids_iter);
      if (blen > max_id_blen) {
        max_id_blen = blen;
      }
      ++id_ct;
      supported_ids_iter = &(supported_ids_iter[blen]);
    } while (*supported_ids_iter);
    // max_id_blen + 4 extra bytes at the end, to support a "maybe" search
    // (yes, this can also be precomputed)
    if (unlikely(pgl_malloc((max_id_blen + 4) * (id_ct + 1) + 1, &id_map))) {
      goto ParseColDescriptor_ret_NOMEM;
    }
    char* sorted_ids = R_CAST(char*, &(id_map[id_ct]));
    supported_ids_iter = supported_ids;
    for (uint32_t id_idx = 0; id_idx != id_ct; ++id_idx) {
      const uint32_t blen = strlen(supported_ids_iter) + 1;
      memcpy(&(sorted_ids[id_idx * max_id_blen]), supported_ids_iter, blen);
      id_map[id_idx] = id_idx;
      supported_ids_iter = &(supported_ids_iter[blen]);
    }
    if (unlikely(SortStrboxIndexedMalloc(id_ct, max_id_blen, sorted_ids, id_map))) {
      goto ParseColDescriptor_ret_NOMEM;
    }
    uint32_t result = *S_CAST(uint32_t*, result_ptr);
    // might not want to bother splitting this into two loops
    if ((col_descriptor_iter[0] == '+') || (col_descriptor_iter[0] == '-')) {
      result |= default_cols_mask;
      char* maybebuf = &(sorted_ids[max_id_blen * id_ct]);
      memcpy_k(maybebuf, "maybe", 5);
      while (1) {
        const char* id_start = &(col_descriptor_iter[1]);
        const char* tok_end = strchrnul(id_start, ',');
        const uint32_t slen = tok_end - id_start;
        int32_t alpha_idx = bsearch_str(id_start, sorted_ids, slen, max_id_blen, id_ct);
        if (unlikely(alpha_idx == -1)) {
          char* write_iter = strcpya_k(g_logbuf, "Error: Unrecognized ID '");
          write_iter = memcpya(write_iter, id_start, slen);
          write_iter = strcpya_k(write_iter, "' in --");
          write_iter = strcpya(write_iter, cur_flag_name);
          snprintf(write_iter, kLogbufSize - kMaxIdSlen - 128, " column set descriptor.\n");
          goto ParseColDescriptor_ret_INVALID_CMDLINE_WW;
        }
        uint32_t shift = id_map[S_CAST(uint32_t, alpha_idx)];
        if (col_descriptor_iter[0] == '+') {
          result |= first_col_shifted << shift;
        } else {
          if (result & (first_col_shifted << shift)) {
            result -= first_col_shifted << shift;
          } else if (slen + 5 < max_id_blen) {
            // special case: if default column set includes e.g. "maybesid",
            // and user types "-sid", that should work
            memcpy(&(maybebuf[5]), id_start, slen);
            alpha_idx = bsearch_str(maybebuf, sorted_ids, slen + 5, max_id_blen, id_ct);
            if (alpha_idx != -1) {
              shift = id_map[S_CAST(uint32_t, alpha_idx)];
              result &= ~(first_col_shifted << shift);
            }
          }
        }
        // bugfix (16 Oct 2017): forgot to switch from !tok_end to !(*tok_end)
        // after use of strchrnul().
        if (!(*tok_end)) {
          break;
        }
        col_descriptor_iter = &(tok_end[1]);
        if (unlikely((col_descriptor_iter[0] != '+') && (col_descriptor_iter[0] != '-'))) {
          goto ParseColDescriptor_ret_MIXED_SIGN;
        }
      }
    } else if (*col_descriptor_iter) {
      while (1) {
        const char* tok_end = strchrnul(col_descriptor_iter, ',');
        const uint32_t slen = tok_end - col_descriptor_iter;
        const int32_t alpha_idx = bsearch_str(col_descriptor_iter, sorted_ids, slen, max_id_blen, id_ct);
        if (unlikely(alpha_idx == -1)) {
          char* write_iter = strcpya_k(g_logbuf, "Error: Unrecognized ID '");
          write_iter = memcpya(write_iter, col_descriptor_iter, slen);
          write_iter = strcpya_k(write_iter, "' in --");
          write_iter = strcpya(write_iter, cur_flag_name);
          snprintf(write_iter, kLogbufSize - kMaxIdSlen - 128, " column set descriptor.\n");
          goto ParseColDescriptor_ret_INVALID_CMDLINE_WW;
        }
        uint32_t shift = id_map[S_CAST(uint32_t, alpha_idx)];
        result |= first_col_shifted << shift;
        if (!(*tok_end)) {
          break;
        }
        col_descriptor_iter = &(tok_end[1]);
        if (unlikely((col_descriptor_iter[0] == '+') || (col_descriptor_iter[0] == '-'))) {
          goto ParseColDescriptor_ret_MIXED_SIGN;
        }
      }
    }
    if (unlikely(prohibit_empty && (!(result & (first_col_shifted * (UINT32_MAX >> (32 - id_ct))))))) {
      char* write_iter = strcpya_k(g_logbuf, "Error: All columns excluded by --");
      write_iter = strcpya(write_iter, cur_flag_name);
      snprintf(write_iter, kLogbufSize - 128, " column set descriptor.\n");
      goto ParseColDescriptor_ret_INVALID_CMDLINE_WW;
    }
    *S_CAST(uint32_t*, result_ptr) = result;
  }
  while (0) {
  ParseColDescriptor_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ParseColDescriptor_ret_MIXED_SIGN:
    snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s column set descriptor (either all column set IDs must be preceded by +/-, or none of them can be).\n", cur_flag_name);
  ParseColDescriptor_ret_INVALID_CMDLINE_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInvalidCmdline;
    break;
  }
  free_cond(id_map);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
