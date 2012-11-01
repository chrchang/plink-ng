// Resources needed across all wdist modules.

#ifndef __WDIST_COMMON_H__
#define __WDIST_COMMON_H__

#include <stdio.h>
#include <stdlib.h>

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

#define MALLOC_DEFAULT_MB 2176

#define CACHELINE 64 // assumed number of bytes per cache line, for alignment
#define CACHELINE_DBL (CACHELINE / sizeof(double))

#define CACHEALIGN(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - 1)))
#define CACHEALIGN_DBL(val) ((val + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1)))

extern const char errstr_fopen[];

int fopen_checked(FILE** target_ptr, const char* fname, const char* mode);

void fclose_cond(FILE* fptr);

// manually managed, very large stack
extern unsigned char* wkspace_base;
extern unsigned long wkspace_left;

unsigned char* wkspace_alloc(unsigned long size);

int wkspace_alloc_c_checked(char** dc_ptr, unsigned long size);

int wkspace_alloc_d_checked(double** dp_ptr, unsigned long size);

int wkspace_alloc_i_checked(int** ip_ptr, unsigned long size);

int wkspace_alloc_uc_checked(unsigned char** ucp_ptr, unsigned long size);

int wkspace_alloc_ui_checked(unsigned int** uip_ptr, unsigned long size);

int wkspace_alloc_ul_checked(unsigned long** ulp_ptr, unsigned long size);

int wkspace_alloc_ull_checked(unsigned long long** ullp_ptr, unsigned long size);

void wkspace_reset(void* new_base);

#endif
