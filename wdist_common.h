// Resources needed across all wdist modules.

#ifndef __WDIST_COMMON_H__
#define __WDIST_COMMON_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>

#ifdef __cplusplus
#include <algorithm>
#endif

#if __LP64__
#include <emmintrin.h>

typedef union {
  __m128i vi;
  __m128d vd;
  unsigned long u8[2];
  double d8[2];
  unsigned int u4[4];
} __uni16;
#endif

#include "zlib-1.2.7/zlib.h"

#define RET_SUCCESS 0
#define RET_NOMEM 1
#define RET_OPEN_FAIL 2
#define RET_INVALID_FORMAT 3
#define RET_CALC_NOT_YET_SUPPORTED 4
#define RET_INVALID_CMDLINE 5
#define RET_WRITE_FAIL 6
#define RET_READ_FAIL 7
#define RET_HELP 8
#define RET_THREAD_CREATE_FAIL 9
#define RET_ALLELE_MISMATCH 10
#define RET_NULL_CALC 11

#define CALC_RELATIONSHIP_COV 1
#define CALC_RELATIONSHIP_SQ 2
#define CALC_RELATIONSHIP_SQ0 4
#define CALC_RELATIONSHIP_TRI 6
#define CALC_RELATIONSHIP_SHAPEMASK 6
#define CALC_RELATIONSHIP_GZ 8
#define CALC_RELATIONSHIP_BIN 0x10
#define CALC_RELATIONSHIP_GRM 0x20
#define CALC_RELATIONSHIP_MASK 0x3f
#define CALC_IBC 0x40
#define CALC_DISTANCE_SQ 0x80
#define CALC_DISTANCE_SQ0 0x100
#define CALC_DISTANCE_TRI 0x180
#define CALC_DISTANCE_SHAPEMASK 0x180
#define CALC_DISTANCE_GZ 0x200
#define CALC_DISTANCE_BIN 0x400
#define CALC_DISTANCE_IBS 0x800
#define CALC_DISTANCE_1_MINUS_IBS 0x1000
#define CALC_DISTANCE_SNPS 0x2000
#define CALC_DISTANCE_FORMATMASK 0x3800
#define CALC_DISTANCE_MASK 0x3f80
#define CALC_PLINK_DISTANCE_MATRIX 0x4000
#define CALC_PLINK_IBS_MATRIX 0x8000
#define CALC_GDISTANCE_MASK 0xff80
#define CALC_LOAD_DISTANCES 0x10000
#define CALC_GROUPDIST 0x20000
#define CALC_REGRESS_DISTANCE 0x40000
#define CALC_UNRELATED_HERITABILITY 0x80000
#define CALC_UNRELATED_HERITABILITY_STRICT 0x100000
#define CALC_FREQ 0x200000
#define CALC_REL_CUTOFF 0x400000
#define CALC_WRITE_SNPLIST 0x800000
#define CALC_GENOME 0x1000000
#define CALC_REGRESS_REL 0x2000000
#define CALC_LD_PRUNE 0x4000000
#define CALC_LD_PRUNE_PAIRWISE 0x8000000
#define CALC_REGRESS_PCS 0x10000000
#define CALC_MAKE_BED 0x20000000
#define CALC_RECODE 0x40000000

#define MALLOC_DEFAULT_MB 2176

#define CACHELINE 64 // assumed number of bytes per cache line, for alignment
#define CACHELINE_DBL (CACHELINE / sizeof(double))

#define CACHEALIGN(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - 1)))
#define CACHEALIGN_DBL(val) ((val + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1)))

#define _FILE_OFFSET_BITS 64
#define MAX_THREADS 127
#define MAX_THREADS_P1 128

#ifdef __LP64__
#define BITCT 64
#else
#define BITCT 32
#endif

#define BITCT2 (BITCT / 2)

// size of generic text line load buffer.  .ped lines can of course be longer
#define MAXLINELEN 131072

#define FNAMESIZE 2048

// allow .mdist.bin.xxxxxxxxxx extension
#define MAX_POST_EXT 22

// number of different types of jackknife values to precompute (x^2, x, y, xy)
#define JACKKNIFE_VALS_DIST 5

// fit 4 pathologically long IDs plus a bit extra
extern char tbuf[];

extern const char errstr_fopen[];
extern const char errstr_append[];
extern const char errstr_thread_create[];

int fopen_checked(FILE** target_ptr, const char* fname, const char* mode);

static inline int fwrite_checked(const void* buf, size_t len, FILE* outfile) {
  if ((!len) || fwrite(buf, len, 1, outfile)) {
    return 0;
  }
  return -1;
}

static inline void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

static inline int fclose_null(FILE** fptr_ptr) {
  int ii;
  ii = fclose(*fptr_ptr);
  *fptr_ptr = NULL;
  return ii;
}

int gzopen_checked(gzFile* target_ptr, const char* fname, const char* mode);

static inline int gzwrite_checked(gzFile gz_outfile, const void* buf, size_t len) {
  if ((!len) || gzwrite(gz_outfile, buf, len)) {
    return 0;
  }
  return -1;
}

static inline int flexwrite_checked(FILE* outfile, gzFile gz_outfile, char* contents, unsigned long len) {
  if (outfile) {
    return fwrite_checked(contents, len, outfile);
  } else {
    return gzwrite_checked(gz_outfile, contents, len);
  }
}

// manually managed, very large stack
extern unsigned char* wkspace_base;
extern unsigned long wkspace_left;

unsigned char* wkspace_alloc(unsigned long size);

static inline int wkspace_alloc_c_checked(char** dc_ptr, unsigned long size) {
  *dc_ptr = (char*)wkspace_alloc(size);
  if (!(*dc_ptr)) {
    return 1;
  }
  return 0;
}

static inline int wkspace_alloc_d_checked(double** dp_ptr, unsigned long size) {
  *dp_ptr = (double*)wkspace_alloc(size);
  if (!(*dp_ptr)) {
    return 1;
  }
  return 0;
}

static inline int wkspace_alloc_i_checked(int** ip_ptr, unsigned long size) {
  *ip_ptr = (int*)wkspace_alloc(size);
  if (!(*ip_ptr)) {
    return 1;
  }
  return 0;
}

static inline int wkspace_alloc_uc_checked(unsigned char** ucp_ptr, unsigned long size) {
  *ucp_ptr = wkspace_alloc(size);
  if (!(*ucp_ptr)) {
    return 1;
  }
  return 0;
}

static inline int wkspace_alloc_ui_checked(unsigned int** uip_ptr, unsigned long size) {
  *uip_ptr = (unsigned int*)wkspace_alloc(size);
  if (!(*uip_ptr)) {
    return 1;
  }
  return 0;
}

static inline int wkspace_alloc_ul_checked(unsigned long** ulp_ptr, unsigned long size) {
  *ulp_ptr = (unsigned long*)wkspace_alloc(size);
  if (!(*ulp_ptr)) {
    return 1;
  }
  return 0;
}

static inline int wkspace_alloc_ull_checked(unsigned long long** ullp_ptr, unsigned long size) {
  *ullp_ptr = (unsigned long long*)wkspace_alloc(size);
  if (!(*ullp_ptr)) {
    return 1;
  }
  return 0;
}

void wkspace_reset(void* new_base);

static inline int no_more_items(char* sptr) {
  return ((!sptr) || (*sptr == '\n') || (*sptr == '\0'));
}

static inline char* skip_initial_spaces(char* sptr) {
  while ((*sptr == ' ') || (*sptr == '\t')) {
    sptr++;
  }
  return sptr;
}

char* item_end(char* sptr);

static inline int is_space_or_eoln(char cc) {
  return ((cc == ' ') || (cc == '\t') || (cc == '\n') || (cc == '\0'));
}

int strlen_se(char* ss);

int strcmp_se(char* s_read, char* s_const, int len);

char* next_item(char* sptr);

char* next_item_mult(char* sptr, unsigned int ct);

static inline void set_bit_noct(unsigned long* exclude_arr, int loc) {
  exclude_arr[loc / BITCT] |= (1LU << (loc % BITCT));
}

void set_bit(unsigned long* bit_arr, int loc, unsigned int* bit_set_ct_ptr);

void set_bit_sub(unsigned long* bit_arr, int loc, unsigned int* bit_unset_ct_ptr);

void clear_bit(unsigned long* exclude_arr, int loc, unsigned int* include_ct_ptr);

static inline int is_set(unsigned long* exclude_arr, int loc) {
  return (exclude_arr[loc / BITCT] >> (loc % BITCT)) & 1;
}

int next_non_set_unsafe(unsigned long* exclude_arr, unsigned int loc);

int next_set_unsafe(unsigned long* include_arr, unsigned int loc);

// These functions seem to optimize better than memset(arr, 0, x) under gcc.
static inline void fill_long_zero(long* larr, size_t size) {
  long* lptr = &(larr[size]);
  while (larr < lptr) {
    *larr++ = 0;
  }
}

static inline void fill_ulong_zero(unsigned long* ularr, size_t size) {
  unsigned long* ulptr = &(ularr[size]);
  while (ularr < ulptr) {
    *ularr++ = 0;
  }
}

static inline void fill_ulong_one(unsigned long* ularr, size_t size) {
  unsigned long* ulptr = &(ularr[size]);
  while (ularr < ulptr) {
    *ularr++ = ~0LU;
  }
}

static inline void fill_int_zero(int* iarr, size_t size) {
#if __LP64__
  fill_long_zero((long*)iarr, size >> 1);
  if (size & 1) {
    iarr[size - 1] = 0;
  }
#else
  fill_long_zero((long*)iarr, size);
#endif
}

static inline void fill_uint_zero(unsigned int* uiarr, size_t size) {
#if __LP64__
  fill_long_zero((long*)uiarr, size >> 1);
  if (size & 1) {
    uiarr[size - 1] = 0;
  }
#else
  fill_long_zero((long*)uiarr, size);
#endif
}

static inline void fill_int_one(unsigned int* iarr, size_t size) {
#if __LP64__
  fill_ulong_one((unsigned long*)iarr, size >> 1);
  if (size & 1) {
    iarr[size - 1] = -1;
  }
#else
  fill_ulong_one((unsigned long*)iarr, size);
#endif
}

static inline void fill_double_zero(double* darr, size_t size) {
  double* dptr = &(darr[size]);
  while (darr < dptr) {
    *darr++ = 0.0;
  }
}

static inline void free_cond(void* memptr) {
  if (memptr) {
    free(memptr);
  }
}

void triangle_fill(unsigned int* target_arr, int ct, int pieces, int parallel_idx, int parallel_tot, int start, int align);

int write_ids(char* outname, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, char* person_ids, unsigned int max_person_id_len);

int distance_d_write_ids(char* outname, char* outname_end, int calculation_type, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, char* person_ids, unsigned int max_person_id_len);

int distance_req(int calculation_type);

int double_cmp(const void* aa, const void* bb);

int double_cmp_deref(const void* aa, const void* bb);

int qsort_ext(char* main_arr, int arr_length, int item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, int secondary_item_len);

int distance_d_write(FILE** outfile_ptr, FILE** outfile2_ptr, FILE** outfile3_ptr, gzFile* gz_outfile_ptr, gzFile* gz_outfile2_ptr, gzFile* gz_outfile3_ptr, int calculation_type, char* outname, char* outname_end, double* dists, double half_marker_ct_recip, unsigned int indiv_ct, int first_indiv_idx, int end_indiv_idx, int parallel_idx, int parallel_tot, unsigned char* membuf);

void collapse_arr(char* item_arr, int fixed_item_len, unsigned long* exclude_arr, int exclude_arr_size);

void init_genrand(unsigned long s);

// unsigned long genrand_int32(void);

// void pick_d(unsigned char* cbuf, unsigned int ct, unsigned int dd);

void pick_d_small(unsigned char* tmp_cbuf, int* ibuf, unsigned int ct, unsigned int dd);

void print_pheno_stdev(double* pheno_d, unsigned int indiv_ct);

unsigned int set_default_jackknife_d(unsigned int ct);

int regress_distance(int calculation_type, double* dists_local, double* pheno_d_local, unsigned int unfiltered_indiv_ct, unsigned long* indiv_exclude, unsigned int indiv_ct_local, unsigned int thread_ct, unsigned long regress_iters, unsigned int regress_d);

#endif // __WDIST_COMMON_H__
