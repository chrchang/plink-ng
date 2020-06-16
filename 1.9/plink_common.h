#ifndef __PLINK_COMMON_H__
#define __PLINK_COMMON_H__

// This file is part of PLINK 1.90, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Resources needed across all plink modules.

#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#ifndef __STDC_FORMAT_MACROS
#  define __STDC_FORMAT_MACROS 1
#endif
#include <inttypes.h>

// avoid compiler warning
#ifndef NDEBUG
  #define NDEBUG
#endif
#include <assert.h>

// Uncomment this to build this without CBLAS/CLAPACK.
// #define NOLAPACK

// Uncomment this to prevent all unstable features from being accessible from
// the command line.
// #define STABLE_BUILD

#define SPECIES_HUMAN 0
#define SPECIES_COW 1
#define SPECIES_DOG 2
#define SPECIES_HORSE 3
#define SPECIES_MOUSE 4
#define SPECIES_RICE 5
#define SPECIES_SHEEP 6
#define SPECIES_UNKNOWN 7
#define SPECIES_DEFAULT SPECIES_HUMAN

#define PROG_NAME_STR "plink"
#define PROG_NAME_CAPS "PLINK"

#ifdef _WIN32
  // needed for MEMORYSTATUSEX
  #ifndef _WIN64
    #define WINVER 0x0500
  #else
    #define __LP64__
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
#else // Unix
  #include <sys/stat.h>
#endif

#ifndef HAVE_NULLPTR
  #ifndef __cplusplus
    #define nullptr NULL
  #else
    #if __cplusplus <= 199711L
      #ifndef nullptr
        #define nullptr NULL
      #endif
    #endif
  #endif
#endif

#ifdef _WIN32
  #define fseeko fseeko64
  #define ftello ftello64
  #include <process.h>
#  undef PRId64
#  undef PRIu64
#  define PRId64 "I64d"
#  define PRIu64 "I64u"
  #define pthread_t HANDLE
  #define THREAD_RET_TYPE unsigned __stdcall
  #define THREAD_RETURN return 0
  #define EOLN_STR "\r\n"
  #define FOPEN_RB "rb"
  #define FOPEN_WB "wb"
  #ifdef _WIN64
    #define getc_unlocked _fgetc_nolock
    #define putc_unlocked _fputc_nolock
  #else
    #define getc_unlocked getc
    #define putc_unlocked putc
  #endif
  #if __cplusplus < 201103L
    #define uint64_t unsigned long long
    #define int64_t long long
  #endif
#else
  #include <pthread.h>
  #define THREAD_RET_TYPE void*
  #define THREAD_RETURN return nullptr
  #ifdef __cplusplus
    #ifndef PRId64
      #define PRId64 "lld"
    #endif
  #endif
  #define EOLN_STR "\n"
  #define FOPEN_RB "r"
  #define FOPEN_WB "w"
  #ifndef __APPLE__
    // argh
    // not sure what the right threshold actually is, but this works for now
    // (may break on gcc <3.0?  but that shouldn't matter anymore)
    // tried defining GCC_VERSION, but that didn't always work
    #if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8)
      #define uint64_t unsigned long long
      #define int64_t long long
    #endif
  #endif
#endif

#ifdef _WIN64
  #define __LP64__
  #define CTZLU __builtin_ctzll
  #define CLZLU __builtin_clzll
#else
  #define CTZLU __builtin_ctzl
  #define CLZLU __builtin_clzl
  #ifndef __LP64__
    // attempt to patch GCC 6 build failure
    #if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8)
      #ifndef uintptr_t
        #define uintptr_t unsigned long
      #endif
      #ifndef intptr_t
        #define intptr_t long
      #endif
    #endif
  #endif
#endif

#ifdef __cplusplus
  #include <algorithm>
#  ifdef _WIN32
// Windows C++11 <algorithm> resets these values :(
#    undef PRIu64
#    undef PRId64
#    define PRIu64 "I64u"
#    define PRId64 "I64d"
#    undef PRIuPTR
#    undef PRIdPTR
#    ifdef __LP64__
#      define PRIuPTR PRIu64
#      define PRIdPTR PRId64
#    else
#      if __cplusplus < 201103L
#        define PRIuPTR "lu"
#        define PRIdPTR "ld"
#      else
#        define PRIuPTR "u"
#        define PRIdPTR "d"
#      endif
#    endif
#  endif
  #define HEADER_INLINE inline
#else
  #define HEADER_INLINE static inline
#endif

// It would be useful to disable compilation on big-endian platforms, but I
// don't see a decent portable way to do this (see e.g. discussion at
// http://esr.ibiblio.org/?p=5095 ).

#ifdef __LP64__
  #ifndef __SSE2__
    // It's obviously possible to support this by writing 64-bit non-SSE2 code
    // shadowing each SSE2 intrinsic, but this almost certainly isn't worth the
    // development/testing effort until regular PLINK 2.0 development is
    // complete.  No researcher has ever asked me for this feature.
    #error "64-bit builds currently require SSE2.  Try producing a 32-bit build instead."
  #endif
  #include <emmintrin.h>

  #define VECFTYPE __m128
  #define VECITYPE __m128i
  #define VECDTYPE __m128d

  // useful because of its bitwise complement: ~ZEROLU is a word with all 1
  // bits, while ~0 is always 32 1 bits.
  #define ZEROLU 0LLU

  // mainly useful for bitshifts: (ONELU << 32) works in 64-bit builds, while
  // (1 << 32) is undefined.  also used to cast some numbers/expressions to
  // uintptr_t (e.g. multiplying an int constant by ONELU widens it to 64 bits
  // only in 64-bit builds; note that 1LU fails on Win64 while 1LLU doesn't do
  // the right thing for 32-bit builds).
  #define ONELU 1LLU

  #ifdef _WIN32 // i.e. Win64

    #ifndef PRIuPTR
      #define PRIuPTR PRIu64
    #endif
    #ifndef PRIdPTR
      #define PRIdPTR PRId64
    #endif
    #define PRIxPTR2 "016I64x"

  #else // not _WIN32

    #ifndef PRIuPTR
      #define PRIuPTR "lu"
    #endif
    #ifndef PRIdPTR
      #define PRIdPTR "ld"
    #endif
    #define PRIxPTR2 "016lx"

  #endif // Win64

  #define VEC_BYTES 16

#else // not __LP64__

  #define ZEROLU 0LU
  #define ONELU 1LU
#  if (__GNUC__ <= 4) && (__GNUC_MINOR__ < 8) && (__cplusplus < 201103L)
#    undef PRIuPTR
#    undef PRIdPTR
#    define PRIuPTR "lu"
#    define PRIdPTR "ld"
#  endif
  #define PRIxPTR2 "08lx"

  // todo: update code so this still works when reduced to 4
  #define VEC_BYTES 8

#endif // __LP64__

// use constexpr for these as soon as compiler support is available on all
// platforms
#define FIVEMASK ((~ZEROLU) / 3)
#define AAAAMASK (FIVEMASK * 2)

#define VEC_BYTES_M1 (VEC_BYTES - 1)
#define VEC_BITS (VEC_BYTES * 8)
#define VEC_BITS_M1 (VEC_BITS - 1)

#ifdef DYNAMIC_ZLIB
  #include <zlib.h>
  #if !defined(ZLIB_VERNUM) || ZLIB_VERNUM < 0x1240
    #error "zlib version 1.2.4 or later required."
  #endif
#else
  #include "../zlib-1.2.11/zlib.h"
#endif
#include "SFMT.h"

// 64MB of non-workspace memory guaranteed for now.
// Currently also serves as the maximum allele length.
#define NON_BIGSTACK_MIN 67108864

#define PI 3.1415926535897932
#define RECIP_2_32 0.00000000023283064365386962890625
#define RECIP_2_53 0.00000000000000011102230246251565404236316680908203125
// floating point comparison-to-nonzero tolerance, currently 2^{-30}
#define EPSILON 0.000000000931322574615478515625
// less tolerant versions (2^{-35}, 2^{-44}) for some exact calculations
#define SMALLISH_EPSILON 0.00000000002910383045673370361328125
#define SMALL_EPSILON 0.00000000000005684341886080801486968994140625
// at least sqrt(SMALL_EPSILON)
#define BIG_EPSILON 0.000000476837158203125
// 53-bit double precision limit
#define DOUBLE_PREC_LIMIT 0.00000000000000011102230246251565404236316680908203125
#define TWO_63 9223372036854775808.0
#define SQRT_HALF 0.70710678118654746

// 2^{-83} bias to give exact tests maximum ability to determine tiny p-values.
// (~2^{-53} is necessary to take advantage of denormalized small numbers, then
// allow tail sum to be up to 2^30.)
#define EXACT_TEST_BIAS 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125

// occasionally used as an infinity substitute that avoids the 32-bit Windows
// performance penalty
// can import from limits.h, we don't bother to include that for now
#ifndef DBL_MAX
  #define DBL_MAX 1.7976931348623157e308
#endif

// not quite the same as FLT_MAX since it's a double-precision constant
#define FLT_MAXD 3.4028234663852886e38

#define RET_SUCCESS 0
#define RET_NOMEM 1
#define RET_OPEN_FAIL 2
#define RET_INVALID_FORMAT 3
#define RET_CALC_NOT_YET_SUPPORTED 4
#define RET_INVALID_CMDLINE 5
#define RET_WRITE_FAIL 6
#define RET_READ_FAIL 7
#define RET_THREAD_CREATE_FAIL 8
#define RET_ALLELE_MISMATCH 9
#define RET_NULL_CALC 10
#define RET_ALL_SAMPLES_EXCLUDED 11
#define RET_ALL_MARKERS_EXCLUDED 12
#define RET_NETWORK 13
#define RET_DEGENERATE_DATA 14
#define LOAD_PHENO_LAST_COL 127

// for 2.0 -> 1.9 backports
#define RET_MALFORMED_INPUT RET_INVALID_FORMAT

#define MISC_AFFECTION_01 1LLU
#define MISC_NONFOUNDERS 2LLU
#define MISC_MAF_SUCC 4LLU
#define MISC_FREQ_COUNTS 8LLU
#define MISC_FREQ_CC 0x10LLU
#define MISC_FREQX 0x20LLU
#define MISC_KEEP_ALLELE_ORDER 0x40LLU
#define MISC_SET_HH_MISSING 0x80LLU
#define MISC_SET_MIXED_MT_MISSING 0x100LLU
#define MISC_KEEP_AUTOCONV 0x200LLU
#define MISC_LOAD_CLUSTER_KEEP_NA 0x400LLU
#define MISC_WRITE_CLUSTER_OMIT_UNASSIGNED 0x800LLU
#define MISC_ALLOW_EXTRA_CHROMS 0x1000LLU
#define MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING 0x2000LLU
#define MISC_MAKE_FOUNDERS_FIRST 0x4000LLU
#define MISC_LASSO_REPORT_ZEROES 0x8000LLU
#define MISC_LASSO_SELECT_COVARS 0x10000LLU
#define MISC_DOUBLE_ID 0x20000LLU
#define MISC_BIALLELIC_ONLY 0x40000LLU
#define MISC_BIALLELIC_ONLY_STRICT 0x80000LLU
#define MISC_BIALLELIC_ONLY_LIST 0x100000LLU
#define MISC_VCF_FILTER 0x200000LLU
#define MISC_GPLINK 0x400000LLU
#define MISC_SNPS_ONLY_JUST_ACGT 0x800000LLU
#define MISC_IMPUTE_SEX 0x1000000LLU
#define MISC_OXFORD_SNPID_CHR 0x2000000LLU
#define MISC_EXTRACT_RANGE 0x4000000LLU
#define MISC_EXCLUDE_RANGE 0x8000000LLU
#define MISC_MERGEX 0x10000000LLU
#define MISC_SET_ME_MISSING 0x20000000LLU
#define MISC_SEXCHECK_YCOUNT 0x40000000LLU
#define MISC_SEXCHECK_YONLY 0x80000000LLU
#define MISC_FAMILY_CLUSTERS 0x100000000LLU
#define MISC_FILL_MISSING_A2 0x200000000LLU
#define MISC_HET_SMALL_SAMPLE 0x400000000LLU
#define MISC_FST_CC 0x800000000LLU
#define MISC_SPLIT_MERGE_NOFAIL 0x1000000000LLU
#define MISC_REAL_REF_ALLELES 0x2000000000LLU
#define MISC_RPLUGIN_DEBUG 0x4000000000LLU
#define MISC_MISSING_GZ 0x8000000000LLU
#define MISC_FREQ_GZ 0x10000000000LLU
#define MISC_HET_GZ 0x20000000000LLU
#define MISC_ALLOW_NO_SAMPLES 0x40000000000LLU
#define MISC_ALLOW_NO_VARS 0x80000000000LLU
#define MISC_VCF_REQUIRE_GT 0x100000000000LLU

// assume for now that .bed must always be accompanied by both .bim and .fam
#define FILTER_ALL_REQ 1LLU
#define FILTER_BIM_REQ 2LLU
#define FILTER_FAM_REQ 4LLU

// ok with --dosage + --map, but not --dosage by itself
#define FILTER_DOSAGEMAP 8LLU

#define FILTER_NODOSAGE 0x10LLU
#define FILTER_NOCNV 0x20LLU
#define FILTER_EXCLUDE_MARKERNAME_SNP 0x40LLU
#define FILTER_BINARY_CASES 0x80LLU
#define FILTER_BINARY_CONTROLS 0x100LLU
#define FILTER_BINARY_FEMALES 0x200LLU
#define FILTER_BINARY_MALES 0x400LLU
#define FILTER_BINARY_FOUNDERS 0x800LLU
#define FILTER_BINARY_NONFOUNDERS 0x1000LLU
#define FILTER_MAKE_FOUNDERS 0x2000LLU
#define FILTER_PRUNE 0x4000LLU
#define FILTER_SNPS_ONLY 0x8000LLU
#define FILTER_TAIL_PHENO 0x10000LLU
#define FILTER_ZERO_CMS 0x20000LLU

#define CALC_RELATIONSHIP 1LLU
#define CALC_IBC 2LLU
#define CALC_DISTANCE 4LLU
#define CALC_PLINK1_DISTANCE_MATRIX 8LLU
#define CALC_PLINK1_IBS_MATRIX 0x10LLU
#define CALC_GDISTANCE_MASK 0x1cLLU
#define CALC_GROUPDIST 0x20LLU
#define CALC_REGRESS_DISTANCE 0x40LLU
#define CALC_UNRELATED_HERITABILITY 0x80LLU
#define CALC_FREQ 0x100LLU
#define CALC_REL_CUTOFF 0x200LLU
#define CALC_WRITE_SNPLIST 0x400LLU
#define CALC_LIST_23_INDELS 0x800LLU
#define CALC_GENOME 0x1000LLU
#define CALC_REGRESS_REL 0x2000LLU
#define CALC_LD_PRUNE 0x4000LLU
#define CALC_DFAM 0x8000LLU
#define CALC_TUCC 0x10000LLU
#define CALC_MAKE_BED 0x20000LLU
#define CALC_RECODE 0x40000LLU
#define CALC_MERGE 0x80000LLU
#define CALC_WRITE_COVAR 0x100000LLU
#define CALC_WRITE_CLUSTER 0x200000LLU
#define CALC_MODEL 0x400000LLU
#define CALC_HARDY 0x800000LLU
#define CALC_GXE 0x1000000LLU
#define CALC_IBS_TEST 0x2000000LLU
#define CALC_CLUSTER 0x4000000LLU
#define CALC_HOMOZYG 0x8000000LLU
#define CALC_NEIGHBOR 0x10000000LLU
#define CALC_GLM 0x20000000LLU
#define CALC_MISSING_REPORT 0x40000000LLU
#define CALC_CMH 0x80000000LLU
#define CALC_HOMOG 0x100000000LLU
#define CALC_LASSO 0x200000000LLU
#define CALC_LASSO_LAMBDA 0x400000000LLU
#define CALC_WRITE_SET 0x800000000LLU
#define CALC_LD 0x1000000000LLU
#define CALC_EPI 0x2000000000LLU
#define CALC_TESTMISS 0x4000000000LLU
#define CALC_TESTMISHAP 0x8000000000LLU
#define CALC_SEXCHECK 0x10000000000LLU
#define CALC_CLUMP 0x20000000000LLU
#define CALC_PCA 0x40000000000LLU
#define CALC_BLOCKS 0x80000000000LLU
#define CALC_SCORE 0x100000000000LLU
#define CALC_MENDEL 0x200000000000LLU
#define CALC_HET 0x400000000000LLU
#define CALC_FLIPSCAN 0x800000000000LLU
#define CALC_TDT 0x1000000000000LLU
#define CALC_MAKE_PERM_PHENO 0x2000000000000LLU
#define CALC_QFAM 0x4000000000000LLU
#define CALC_FST 0x8000000000000LLU
#define CALC_SHOW_TAGS 0x10000000000000LLU
#define CALC_MAKE_BIM 0x20000000000000LLU
#define CALC_MAKE_FAM 0x40000000000000LLU
#define CALC_WRITE_VAR_RANGES 0x80000000000000LLU
#define CALC_DUPVAR 0x100000000000000LLU
#define CALC_RPLUGIN 0x200000000000000LLU
#define CALC_ONLY_BIM (CALC_WRITE_SET | CALC_WRITE_SNPLIST | CALC_WRITE_VAR_RANGES | CALC_LIST_23_INDELS | CALC_MAKE_BIM | CALC_DUPVAR)
#define CALC_ONLY_FAM (CALC_MAKE_PERM_PHENO | CALC_WRITE_COVAR | CALC_MAKE_FAM)
// only room for 6 more basic commands before we need to switch from a single
// uint64_t to uintptr_t*/is_set()/etc.

// necessary to patch heterozygous haploids/female Y chromosome genotypes
// during loading?
#define XMHH_EXISTS 1
#define Y_FIX_NEEDED 2
#define NXMHH_EXISTS 4

#define ALLELE_RECODE 1
#define ALLELE_RECODE_MULTICHAR 2
#define ALLELE_RECODE_ACGT 4

// 0 = non-explicit error
#define VCF_HALF_CALL_ERROR 1
#define VCF_HALF_CALL_MISSING 2
#define VCF_HALF_CALL_HAPLOID 3
#define VCF_HALF_CALL_REFERENCE 4

#define M23_MALE 1
#define M23_FEMALE 2
#define M23_FORCE_MISSING_SEX 4
#define M23_SEX 7

#define MARKER_CMS_OPTIONAL 1
#define MARKER_CMS_FORCED 2

#define UNSORTED_CHROM 1
#define UNSORTED_BP 2
#define UNSORTED_CM 4
#define UNSORTED_SPLIT_CHROM 8

#define ALLOW_NO_SEX 1
#define MUST_HAVE_SEX 2

#define LGEN_REFERENCE 1
#define LGEN_ALLELE_COUNT 2

#define PHENO_ALL 1
#define PHENO_MERGE 2

#define FAM_COL_1 1
#define FAM_COL_34 2
#define FAM_COL_5 4
#define FAM_COL_6 8
#define FAM_COL_13456 15

#define COVAR_KEEP_PHENO_ON_MISSING_COV 1
#define COVAR_NAME 2
#define COVAR_NUMBER 4
#define COVAR_NO_CONST 8
#define COVAR_ALLOW_NONE 0x10

#define DISTANCE_SQ 1
#define DISTANCE_SQ0 2
#define DISTANCE_TRI 3
#define DISTANCE_SHAPEMASK 3
#define DISTANCE_GZ 4
#define DISTANCE_BIN 8
#define DISTANCE_BIN4 0x10
#define DISTANCE_IBS 0x20
#define DISTANCE_1_MINUS_IBS 0x40
#define DISTANCE_ALCT 0x80
#define DISTANCE_TYPEMASK 0xe0
#define DISTANCE_FLAT_MISSING 0x100
#define DISTANCE_CLUSTER 0x200
#define DISTANCE_WTS_NOHEADER 0x400

#define RECODE_01 1
#define RECODE_12 2
#define RECODE_TAB 4
#define RECODE_DELIMX 8
#define RECODE_23 0x10
#define RECODE_A 0x20
#define RECODE_A_TRANSPOSE 0x40
#define RECODE_AD 0x80
#define RECODE_BEAGLE 0x100
#define RECODE_BEAGLE_NOMAP 0x200
#define RECODE_BIMBAM 0x400
#define RECODE_BIMBAM_1CHR 0x800
#define RECODE_COMPOUND 0x1000
#define RECODE_FASTPHASE 0x2000
#define RECODE_FASTPHASE_1CHR 0x4000
#define RECODE_HV 0x8000
#define RECODE_HV_1CHR 0x10000
#define RECODE_LGEN 0x20000
#define RECODE_LGEN_REF 0x40000
#define RECODE_LIST 0x80000
#define RECODE_OXFORD 0x100000
#define RECODE_RLIST 0x200000
#define RECODE_STRUCTURE 0x400000
#define RECODE_TRANSPOSE 0x800000
#define RECODE_PED 0x1000000
#define RECODE_VCF 0x2000000
#define RECODE_TYPEMASK 0x3fffff0
#define RECODE_FID 0x4000000
#define RECODE_IID 0x8000000
#define RECODE_INCLUDE_ALT 0x10000000
#define RECODE_BGZ 0x20000000
#define RECODE_GEN_GZ 0x40000000
#define RECODE_OMIT_NONMALE_Y 0x80000000U

#define GENOME_OUTPUT_GZ 1
#define GENOME_REL_CHECK 2
#define GENOME_OUTPUT_FULL 4
#define GENOME_IBD_UNBOUNDED 8
#define GENOME_NUDGE 0x10
// separate flag to ensure behavior is unchanged under --unbounded
#define GENOME_FILTER_PI_HAT 0x20

#define WRITE_COVAR_PHENO 1
#define WRITE_COVAR_NO_PARENTS 2
#define WRITE_COVAR_NO_SEX 4
#define WRITE_COVAR_FEMALE_2 8
#define WRITE_COVAR_DUMMY 0x10
#define WRITE_COVAR_DUMMY_NO_ROUND 0x20

#define MERGE_MODE_MASK 7
#define MERGE_EQUAL_POS 8
#define MERGE_BINARY 16
#define MERGE_LIST 32

#define SAMPLE_SORT_NONE 1
#define SAMPLE_SORT_NATURAL 2
#define SAMPLE_SORT_ASCII 4
#define SAMPLE_SORT_FILE 8

#define HWE_MIDP 1
#define HWE_THRESH_MIDP 2
#define HWE_THRESH_ALL 4
#define HWE_GZ 8

#define MENDEL_FILTER 1
#define MENDEL_FILTER_VAR_FIRST 2
#define MENDEL_DUOS 4
#define MENDEL_MULTIGEN 8
#define MENDEL_SUMMARIES_ONLY 0x10

#define DUMMY_MISSING_GENO 1
#define DUMMY_MISSING_PHENO 2
#define DUMMY_SCALAR_PHENO 4
#define DUMMY_ACGT 8
#define DUMMY_1234 0x10
#define DUMMY_12 0x20

#define SIMULATE_QT 1
#define SIMULATE_TAGS 2
#define SIMULATE_HAPS 4
#define SIMULATE_ACGT 8
#define SIMULATE_1234 0x10
#define SIMULATE_12 0x20

#define MODEL_ASSOC 1
#define MODEL_FISHER 2
#define MODEL_FISHER_MIDP 4
#define MODEL_PERM 8
#define MODEL_MPERM 0x10
#define MODEL_GENEDROP 0x20
#define MODEL_PERM_COUNT 0x40
#define MODEL_ASSOC_COUNTS 0x80
#define MODEL_ASSOC_FDEPR 0x100
#define MODEL_DMASK 0x1a6
#define MODEL_QT_MEANS 0x200
#define MODEL_PDOM 0x400
#define MODEL_PREC 0x800
#define MODEL_PGEN 0x1000
#define MODEL_PTREND 0x2000
#define MODEL_TRENDONLY 0x4000
#define MODEL_PMASK (MODEL_PDOM | MODEL_PREC | MODEL_PGEN | MODEL_PTREND | MODEL_TRENDONLY)
#define MODEL_LIN 0x8000
#define MODEL_QMASK (MODEL_QT_MEANS | MODEL_LIN)
#define MODEL_SET_TEST 0x10000

#define GLM_LOGISTIC 1
#define GLM_PERM 2
#define GLM_MPERM 4
#define GLM_GENEDROP 8
#define GLM_PERM_COUNT 0x10
#define GLM_GENOTYPIC 0x20
#define GLM_HETHOM 0x40
#define GLM_DOMINANT 0x80
#define GLM_RECESSIVE 0x100
#define GLM_NO_SNP 0x200
#define GLM_HIDE_COVAR 0x400
#define GLM_SEX 0x800
#define GLM_NO_X_SEX 0x1000
#define GLM_INTERACTION 0x2000
#define GLM_STANDARD_BETA 0x4000
#define GLM_BETA 0x8000
#define GLM_TEST_ALL 0x10000
#define GLM_CONDITION_DOMINANT 0x20000
#define GLM_CONDITION_RECESSIVE 0x40000
#define GLM_SET_TEST 0x80000
#define GLM_NO_SNP_EXCL 0x831ea
#define GLM_INTERCEPT 0x100000

#define MPERM_DUMP_BEST 1
#define MPERM_DUMP_ALL 2

// (2^31 - 1000001) / 2
#define APERM_MAX 1073241823

#define ADJUST_GC 2
#define ADJUST_LOG10 4
#define ADJUST_QQ 8
#define ADJUST_LAMBDA 16

#define DUPVAR_REF 1
#define DUPVAR_IDS_ONLY 2
#define DUPVAR_SUPPRESS_FIRST 4

#define CNV_MAKE_MAP 1
#define CNV_MAKE_MAP_LONG 2
#define CNV_CHECK_NO_OVERLAP 4
#define CNV_DEL 8
#define CNV_DUP 0x10
#define CNV_WRITE_FREQ 0x20
#define CNV_UNIQUE 0x40
#define CNV_DROP_NO_SEGMENT 0x80
#define CNV_SAMPLE_PERM 0x100
#define CNV_ENRICHMENT_TEST 0x200
#define CNV_TEST 0x400
#define CNV_TEST_FORCE_1SIDED 0x800
#define CNV_TEST_FORCE_2SIDED 0x1000
#define CNV_TEST_REGION 0x2000
#define CNV_TRACK 0x4000
#define CNV_SEGLIST 0x8000
#define CNV_REPORT_REGIONS 0x10000
#define CNV_VERBOSE_REPORT_REGIONS 0x20000
#define CNV_WRITE 0x40000
#define CNV_EXCLUDE_OFF_BY_1 0x80000

#define CNV_INTERSECT 1
#define CNV_EXCLUDE 2
#define CNV_COUNT 4

#define CNV_OVERLAP 1
#define CNV_OVERLAP_REGION 2
#define CNV_OVERLAP_UNION 3
#define CNV_DISRUPT 4

#define CNV_FREQ_EXCLUDE_ABOVE 1
#define CNV_FREQ_EXCLUDE_BELOW 2
#define CNV_FREQ_EXCLUDE_EXACT 4
#define CNV_FREQ_INCLUDE_EXACT 8
#define CNV_FREQ_FILTER 15
#define CNV_FREQ_OVERLAP 16
#define CNV_FREQ_METHOD2 32

#define SEGMENT_GROUP 1

// default jackknife iterations
#define ITERS_DEFAULT 100000
#define MAX_PCS_DEFAULT 20

#define BIGSTACK_MIN_MB 64
#define BIGSTACK_DEFAULT_MB 2048

#ifdef __LP64__
  #define BITCT 64

  // unions generally shouldn't be used for reinterpret_cast's job (memcpy is
  // the right C-compatible way), but vectors are an exception to this rule.
  typedef union {
    VECFTYPE vf;
    VECITYPE vi;
    VECDTYPE vd;
    uintptr_t u8[VEC_BITS / BITCT];
    double d8[VEC_BYTES / sizeof(double)];
    float f4[VEC_BYTES / sizeof(float)];
    uint32_t u4[VEC_BYTES / sizeof(int32_t)];
  } __univec;
#else
  #define BITCT 32
#endif

#define BITCT2 (BITCT / 2)
#define BYTECT (BITCT / 8)
#define BYTECT4 (BITCT / 32)
#define VEC_WORDS (VEC_BITS / BITCT)
#define VEC_INT32 (VEC_BYTES / 4)

// assumed number of bytes per cache line, for alignment
#define CACHELINE 64

#define CACHELINE_BIT (CACHELINE * 8)
#define CACHELINE_INT32 (CACHELINE / 4)
#define CACHELINE_INT64 (CACHELINE / 8)
#define CACHELINE_WORD (CACHELINE / BYTECT)
#define CACHELINE_DBL (CACHELINE / 8)

// alignment must be a power of 2
HEADER_INLINE uintptr_t round_up_pow2(uintptr_t val, uintptr_t alignment) {
  uintptr_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return (val + alignment_m1) & (~alignment_m1);
}

#define BITCT_TO_VECCT(val) (((val) + (VEC_BITS - 1)) / VEC_BITS)
#define BITCT_TO_WORDCT(val) (((val) + (BITCT - 1)) / BITCT)
#define BITCT_TO_ALIGNED_WORDCT(val) (VEC_WORDS * BITCT_TO_VECCT(val))

#define QUATERCT_TO_VECCT(val) (((val) + ((VEC_BITS / 2) - 1)) / (VEC_BITS / 2))
#define QUATERCT_TO_WORDCT(val) (((val) + (BITCT2 - 1)) / BITCT2)
#define QUATERCT_TO_ALIGNED_WORDCT(val) (VEC_WORDS * QUATERCT_TO_VECCT(val))

// todo: get rid of (BITCT_TO_WORDCT(x) == QUATERCT_TO_VECCT(x)) and similar
// assumptions, in preparation for AVX2

#ifdef __LP64__
#define round_up_pow2_ull round_up_pow2
#else
HEADER_INLINE uint64_t round_up_pow2_ull(uint64_t val, uint64_t alignment) {
  uint64_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return (val + alignment_m1) & (~alignment_m1);
}
#endif

// 32-bit instead of word-length bitwise not here, when val can be assumed to
// be 32-bit.
// (note that the sizeof operator "returns" an uintptr_t, not a uint32_t; hence
// the lack of sizeof in the CACHELINE_INT32, etc. definitions.)
HEADER_INLINE uint32_t round_up_pow2_ui(uint32_t val, uint32_t alignment) {
  uint32_t alignment_m1 = alignment - 1;
  assert(!(alignment & alignment_m1));
  return (val + alignment_m1) & (~alignment_m1);
}

#define MAXV(aa, bb) (((bb) > (aa))? (bb) : (aa))
#define MINV(aa, bb) (((aa) > (bb))? (bb) : (aa))

#ifdef _WIN32
// if MAX_THREADS > 65, single WaitForMultipleObjects calls must be converted
// into loops
  #define MAX_THREADS 64
  #define MAX_THREADS_P1 65
#else
// shouldn't be larger than MODEL_BLOCKSIZE for now
  #define MAX_THREADS 512
  #define MAX_THREADS_P1 513
#endif

// defined as a macro since type of idx can vary; might want a debug
// compilation mode which performs type-checking, though
#define EXTRACT_2BIT_GENO(ulptr, idx) (((ulptr)[(idx) / BITCT2] >> (2 * ((idx) % BITCT2))) & 3)

// generic maximum line length.  .ped/.vcf/etc. lines can of course be longer
#define MAXLINELEN 131072

// must be at least 2 * MAXLINELEN + 2 to support generic token loader.
#define TEXTBUF_SIZE (2 * MAXLINELEN + 256)

// Maximum length of chromosome, variant, FID, IID, cluster, and set IDs (not
// including terminating null, that's what _P1 is for).  This value supports up
// to 8 IDs per line (maximum so far is 5, for e.g. --hom).
#define MAX_ID_SLEN 16000

#define MAX_ID_BLEN (MAX_ID_SLEN + 1)
#define MAX_ID_SLEN_STR "16000"

// Maximum size of "dynamically" allocated line load buffer.  (This is the
// limit that applies to .vcf and similar files.)  Inconvenient to go higher
// since fgets() takes a int32_t size argument.
#define MAXLINEBUFLEN 0x7fffffc0

// Default --perm-batch-size value in most contexts.  It may actually be better
// to *avoid* a power of two due to the need for transpositions involving this
// stride; see e.g. http://danluu.com/3c-conflict/ ; try 448 instead?  This
// should be tested during PLINK 2.0 development.
#define DEFAULT_PERM_BATCH_SIZE 512

// note that this is NOT foolproof: see e.g.
// http://insanecoding.blogspot.com/2007/11/pathmax-simply-isnt.html .  (This
// is why I haven't bothered with OS-based #ifdefs here.)  But it should be
// good enough in practice.
#define FNAMESIZE 4096

// allow extensions like .model.trend.fisher.set.score.adjusted
#define MAX_POST_EXT 39

// number of types of jackknife values to precompute (x^2, y^2, x, y, xy)
#define JACKKNIFE_VALS_DIST 5
#define JACKKNIFE_VALS_GROUPDIST 3

#ifdef __LP64__
  // number of snp-major .bed lines to read at once for distance calc if
  // exponent is nonzero.
  #define MULTIPLEX_DIST_EXP 64
  // number of snp-major .bed lines to read at once for relationship calc
  #define MULTIPLEX_REL 60
#else
  // N.B. 32-bit version not as carefully tested or optimized, but I'll try to
  // make sure it works properly
  #define MULTIPLEX_DIST_EXP 28
  #define MULTIPLEX_REL 30
#endif

// used to size a few tables
#define EXPECTED_MISSING_FREQ 0.05

// load markers in blocks to enable multithreading and, for quantitative
// phenotypes, PERMORY-style LD exploitation
#define MODEL_BLOCKSIZE 1024
#define MODEL_BLOCKKEEP 64

// string hash table constants, currently only relevant for merge operations
// and annotate()
// (dynamic sizing used for main marker name lookup)

// last prime before 2^19
// size chosen to be likely to fit in L3 cache
#define HASHSIZE 524287
#define HASHSIZE_S 524287

#ifdef __LP64__
#define HASHMEM 4194304
#else
#define HASHMEM 2097152
#endif

typedef struct {
  uint32_t min;
  uint32_t max;
  double alpha;
  double beta;
  double init_interval;
  double interval_slope;
} Aperm_info;

// Generic text I/O buffer: any function which reads from/writes to a text file
// or the console may clobber it.  Sized to fit two MAXLINELEN-length lines
// plus a bit extra.
extern char g_textbuf[];

extern const char g_one_char_strs[];
extern const char* g_missing_geno_ptr;
extern const char* g_output_missing_geno_ptr;

HEADER_INLINE const char* cond_replace(const char* ss, const char* match_str, const char* replace_str) {
  return (ss != match_str)? ss : replace_str;
}

uint32_t aligned_malloc(uintptr_t size, uintptr_t** aligned_pp);

void aligned_free(uintptr_t* aligned_pp);

HEADER_INLINE void aligned_free_cond(uintptr_t* aligned_ptr) {
  if (aligned_ptr) {
    aligned_free(aligned_ptr);
  }
}

HEADER_INLINE void aligned_free_null(uintptr_t** aligned_pp) {
  aligned_free(*aligned_pp);
  *aligned_pp = nullptr;
}

HEADER_INLINE void aligned_free_cond_null(uintptr_t** aligned_pp) {
  if (*aligned_pp) {
    aligned_free(*aligned_pp);
    *aligned_pp = nullptr;
  }
}

extern uintptr_t g_failed_alloc_attempt_size;

extern sfmt_t g_sfmt;

// file-scope string constants don't always have the g_ prefix, but multi-file
// constants are always tagged.
extern const char g_errstr_fopen[];
extern const char g_cmdline_format_str[];

extern FILE* g_logfile;

// mostly-safe sprintf buffer.  warning: do NOT put allele codes or
// arbitrary-length lists in here.
extern char g_logbuf[];

extern uint32_t g_debug_on;
extern uint32_t g_log_failed;

// should remove this global: multithreaded functions should use a file-local
// thread_ct which will occasionally be smaller due to job size.
extern uint32_t g_thread_ct;

typedef struct ll_str_struct {
  struct ll_str_struct* next;
  char ss[];
} Ll_str;

typedef struct ll_ctstr_entry_struct {
  struct ll_ctstr_entry_struct* next;
  uint32_t ct;
  char ss[];
} Ll_ctstr_entry;

typedef struct two_col_params_struct {
  uint32_t colx;
  uint32_t colid;
  uint32_t skip;
  char skipchar;
  char fname[];
} Two_col_params;

typedef struct range_list_struct {
  char* names;
  unsigned char* starts_range;
  uint32_t name_ct;
  uint32_t name_max_len;
} Range_list;

// Pushes a copy of ss (allocated via malloc) onto ll_stack.
uint32_t push_ll_str(const char* ss, Ll_str** ll_stack_ptr);

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

#ifdef STABLE_BUILD
  #define UNSTABLE(val) sptr = strcpya(&(g_logbuf[9]), val); goto main_unstable_disabled
#else
  #define UNSTABLE(val)
#endif

int32_t fopen_checked(const char* fname, const char* mode, FILE** target_ptr);

HEADER_INLINE int32_t putc_checked(int32_t ii, FILE* outfile) {
  putc_unlocked(ii, outfile);
  return ferror(outfile);
}

HEADER_INLINE int32_t fputs_checked(const char* ss, FILE* outfile) {
  fputs(ss, outfile);
  return ferror(outfile);
}

// This must be used for all fwrite() calls where len could be >= 2^31, since
// OS X raw fwrite() doesn't work in that case.
int32_t fwrite_checked(const void* buf, size_t len, FILE* outfile);

HEADER_INLINE int32_t fread_checked(char* buf, uintptr_t len, FILE* infile, uintptr_t* bytes_read_ptr) {
  *bytes_read_ptr = fread(buf, 1, len, infile);
  return ferror(infile);
}

HEADER_INLINE void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

HEADER_INLINE int32_t fclose_null(FILE** fptr_ptr) {
  int32_t ii = ferror(*fptr_ptr);
  int32_t jj = fclose(*fptr_ptr);
  *fptr_ptr = nullptr;
  return ii || jj;
}

// Also sets 128k read buffer.  Can return RET_OPEN_FAIL or RET_NOMEM.
int32_t gzopen_read_checked(const char* fname, gzFile* gzf_ptr);

// pigz interface should be used for writing .gz files.

HEADER_INLINE int32_t gzclose_null(gzFile* gzf_ptr) {
  int32_t ii = gzclose(*gzf_ptr);
  *gzf_ptr = nullptr;
  return (ii != Z_OK);
}

HEADER_INLINE void gzclose_cond(gzFile gz_infile) {
  if (gz_infile) {
    gzclose(gz_infile);
  }
}

HEADER_INLINE int32_t flexwrite_checked(const void* buf, size_t len, uint32_t output_gz, FILE* outfile, gzFile gz_outfile) {
  if (!output_gz) {
    return fwrite_checked(buf, len, outfile);
  } else {
    return (!gzwrite(gz_outfile, buf, len));
  }
}

HEADER_INLINE int32_t flexputc_checked(int32_t ii, uint32_t output_gz, FILE* outfile, gzFile gz_outfile) {
  if (!output_gz) {
    putc(ii, outfile);
    return ferror(outfile);
  } else {
    return (gzputc(gz_outfile, ii) == -1);
  }
}

HEADER_INLINE int32_t flexputs_checked(const char* ss, uint32_t output_gz, FILE* outfile, gzFile gz_outfile) {
  if (!output_gz) {
    return fputs_checked(ss, outfile);
  } else {
    return (gzputs(gz_outfile, ss) == -1);
  }
}

HEADER_INLINE int32_t flexclose_null(uint32_t output_gz, FILE** fptr_ptr, gzFile* gzf_ptr) {
  if (!output_gz) {
    return fclose_null(fptr_ptr);
  } else {
    return gzclose_null(gzf_ptr);
  }
}

// manually managed, very large double-ended stack
extern unsigned char* g_bigstack_base;
extern unsigned char* g_bigstack_end;

HEADER_INLINE uintptr_t bigstack_left() {
  return (((uintptr_t)g_bigstack_end) - ((uintptr_t)g_bigstack_base));
}

// Basic 64-byte-aligned allocation at bottom of stack.
unsigned char* bigstack_alloc(uintptr_t size);


// Typesafe, return-0-iff-success interfaces.  (See also bigstack_calloc_...
// further below.)
HEADER_INLINE int32_t bigstack_alloc_c(uintptr_t ct, char** cp_ptr) {
  *cp_ptr = (char*)bigstack_alloc(ct);
  return !(*cp_ptr);
}

HEADER_INLINE int32_t bigstack_alloc_d(uintptr_t ct, double** dp_ptr) {
  *dp_ptr = (double*)bigstack_alloc(ct * sizeof(double));
  return !(*dp_ptr);
}

HEADER_INLINE int32_t bigstack_alloc_f(uintptr_t ct, float** fp_ptr) {
  *fp_ptr = (float*)bigstack_alloc(ct * sizeof(float));
  return !(*fp_ptr);
}

HEADER_INLINE int32_t bigstack_alloc_i(uintptr_t ct, int32_t** ip_ptr) {
  *ip_ptr = (int32_t*)bigstack_alloc(ct * sizeof(int32_t));
  return !(*ip_ptr);
}

HEADER_INLINE int32_t bigstack_alloc_uc(uintptr_t ct, unsigned char** ucp_ptr) {
  *ucp_ptr = bigstack_alloc(ct);
  return !(*ucp_ptr);
}

HEADER_INLINE int32_t bigstack_alloc_ui(uintptr_t ct, uint32_t** uip_ptr) {
  *uip_ptr = (uint32_t*)bigstack_alloc(ct * sizeof(int32_t));
  return !(*uip_ptr);
}

HEADER_INLINE int32_t bigstack_alloc_ul(uintptr_t ct, uintptr_t** ulp_ptr) {
  *ulp_ptr = (uintptr_t*)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*ulp_ptr);
}

HEADER_INLINE int32_t bigstack_alloc_ll(uintptr_t ct, int64_t** llp_ptr) {
  *llp_ptr = (int64_t*)bigstack_alloc(ct * sizeof(int64_t));
  return !(*llp_ptr);
}

HEADER_INLINE int32_t bigstack_alloc_ull(uintptr_t ct, uint64_t** ullp_ptr) {
  *ullp_ptr = (uint64_t*)bigstack_alloc(ct * sizeof(int64_t));
  return !(*ullp_ptr);
}

HEADER_INLINE void bigstack_reset(const void* new_base) {
  g_bigstack_base = (unsigned char*)new_base;
}

HEADER_INLINE void bigstack_end_reset(const void* new_end) {
  g_bigstack_end = (unsigned char*)new_end;
}

HEADER_INLINE void bigstack_double_reset(const void* new_base, const void* new_end) {
  bigstack_reset(new_base);
  bigstack_end_reset(new_end);
}

void bigstack_shrink_top(const void* rebase, uintptr_t new_size);

#define END_ALLOC_CHUNK 16
#define END_ALLOC_CHUNK_M1 (END_ALLOC_CHUNK - 1)

HEADER_INLINE void bigstack_end_set(const void* unaligned_end) {
  g_bigstack_end = (unsigned char*)(((uintptr_t)unaligned_end) & (~(END_ALLOC_CHUNK_M1 * ONELU)));
}

// assumes size is divisible by END_ALLOC_CHUNK
// (no value in directly calling this with a constant size parameter: the
// compiler will properly optimize a bigstack_end_alloc() call)
unsigned char* bigstack_end_alloc_presized(uintptr_t size);

HEADER_INLINE unsigned char* bigstack_end_alloc(uintptr_t size) {
  // multiplication by ONELU is one way to widen an int to word-size.
  size = round_up_pow2(size, END_ALLOC_CHUNK);
  return bigstack_end_alloc_presized(size);
}

#define bigstack_end_aligned_alloc bigstack_end_alloc

HEADER_INLINE int32_t bigstack_end_alloc_c(uintptr_t ct, char** cp_ptr) {
  *cp_ptr = (char*)bigstack_end_alloc(ct);
  return !(*cp_ptr);
}

HEADER_INLINE int32_t bigstack_end_alloc_d(uintptr_t ct, double** dp_ptr) {
  *dp_ptr = (double*)bigstack_end_alloc(ct * sizeof(double));
  return !(*dp_ptr);
}

HEADER_INLINE int32_t bigstack_end_alloc_f(uintptr_t ct, float** fp_ptr) {
  *fp_ptr = (float*)bigstack_end_alloc(ct * sizeof(float));
  return !(*fp_ptr);
}

HEADER_INLINE int32_t bigstack_end_alloc_i(uintptr_t ct, int32_t** ip_ptr) {
  *ip_ptr = (int32_t*)bigstack_end_alloc(ct * sizeof(int32_t));
  return !(*ip_ptr);
}

HEADER_INLINE int32_t bigstack_end_alloc_uc(uintptr_t ct, unsigned char** ucp_ptr) {
  *ucp_ptr = bigstack_end_alloc(ct);
  return !(*ucp_ptr);
}

HEADER_INLINE int32_t bigstack_end_alloc_ui(uintptr_t ct, uint32_t** uip_ptr) {
  *uip_ptr = (uint32_t*)bigstack_end_alloc(ct * sizeof(int32_t));
  return !(*uip_ptr);
}

HEADER_INLINE int32_t bigstack_end_alloc_ul(uintptr_t ct, uintptr_t** ulp_ptr) {
  *ulp_ptr = (uintptr_t*)bigstack_end_alloc(ct * sizeof(intptr_t));
  return !(*ulp_ptr);
}

HEADER_INLINE int32_t bigstack_end_alloc_ll(uintptr_t ct, int64_t** llp_ptr) {
  *llp_ptr = (int64_t*)bigstack_end_alloc(ct * sizeof(int64_t));
  return !(*llp_ptr);
}

HEADER_INLINE int32_t bigstack_end_alloc_ull(uintptr_t ct, uint64_t** ullp_ptr) {
  *ullp_ptr = (uint64_t*)bigstack_end_alloc(ct * sizeof(int64_t));
  return !(*ullp_ptr);
}

HEADER_INLINE int32_t bigstack_end_alloc_llstr(uintptr_t str_bytes, Ll_str** llstrp_ptr) {
  *llstrp_ptr = (Ll_str*)bigstack_end_alloc(str_bytes + sizeof(Ll_str));
  return !(*llstrp_ptr);
}


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

// could assert ucc is not a space/tab
#define is_eoln_kns is_space_or_eoln

HEADER_INLINE int32_t is_eoln_char(unsigned char ucc) {
  return (ucc < 32) && (ucc != 9);
}

HEADER_INLINE int32_t is_eoln_or_comment_kns(unsigned char ucc) {
  return (ucc < 32) || (ucc == '#');
}

HEADER_INLINE int32_t no_more_tokens_kns(const char* sptr) {
  return ((!sptr) || is_eoln_kns(*sptr));
}

HEADER_INLINE char* skip_initial_spaces(char* sptr) {
  while ((*sptr == ' ') || (*sptr == '\t')) {
    sptr++;
  }
  return sptr;
}

/*
HEADER_INLINE int32_t is_space_or_eoln(unsigned char cc) {
  // ' ', \t, \n, \0, \r
#ifdef __LP64__
  return (ucc <= 32) && (0x100002601LLU & (1LLU << ucc));
#else
  return ((ucc <= 32) && ((ucc == ' ') || (0x2601LU & (ONELU << ucc))));
#endif
}
*/

// Returns whether uppercased ss matches nonempty fixed_str.  Assumes fixed_str
// contains nothing but letters and a null terminator.
uint32_t match_upper(const char* ss, const char* fixed_str);

uint32_t match_upper_counted(const char* ss, const char* fixed_str, uint32_t ct);

// Reads an integer in [1, cap].  Assumes first character is nonspace.  Has the
// overflow detection atoi() lacks.
#ifdef __LP64__
uint32_t scan_posint_capped(const char* ss, uint64_t cap, uint32_t* valp);

uint32_t scan_uint_capped(const char* ss, uint64_t cap, uint32_t* valp);

uint32_t scan_int_abs_bounded(const char* ss, uint64_t bound, int32_t* valp);
#else // not __LP64__
// Need to be more careful in 32-bit case due to overflow.
// A funny-looking div_10/mod_10 interface is used since the cap will usually
// be a constant, and we want the integer division/modulus to occur at compile
// time.
uint32_t scan_posint_capped32(const char* ss, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp);

uint32_t scan_uint_capped32(const char* ss, uint32_t cap_div_10, uint32_t cap_mod_10, uint32_t* valp);

uint32_t scan_int_abs_bounded32(const char* ss, uint32_t bound_div_10, uint32_t bound_mod_10, int32_t* valp);

  #define scan_posint_capped(aa, bb, cc) scan_posint_capped32((aa), (bb) / 10, (bb) % 10, (cc))

  #define scan_uint_capped(aa, bb, cc) scan_uint_capped32((aa), (bb) / 10, (bb) % 10, (cc))

  #define scan_int_abs_bounded(aa, bb, cc) scan_int_abs_bounded32((aa), (bb) / 10, (bb) % 10, (cc))
#endif

// intentionally rejects -2^31 for now
HEADER_INLINE uint32_t scan_int32(const char* ss, int32_t* valp) {
  return scan_int_abs_bounded(ss, 0x7fffffff, valp);
}

// default cap = 0x7ffffffe
HEADER_INLINE uint32_t scan_posint_defcap(const char* ss, uint32_t* valp) {
  return scan_posint_capped(ss, 0x7ffffffe, valp);
}

HEADER_INLINE uint32_t scan_uint_defcap(const char* ss, uint32_t* valp) {
  return scan_uint_capped(ss, 0x7ffffffe, valp);
}

HEADER_INLINE uint32_t scan_int_abs_defcap(const char* ss, int32_t* valp) {
  return scan_int_abs_bounded(ss, 0x7ffffffe, valp);
}

HEADER_INLINE uint32_t scan_uint_icap(const char* ss, uint32_t* valp) {
  return scan_uint_capped(ss, 0x7fffffff, valp);
}

uint32_t scan_posintptr(const char* ss, uintptr_t* valp);

HEADER_INLINE uint32_t scan_double(const char* ss, double* valp) {
  char* ss2;
  *valp = strtod(ss, &ss2);
  return (ss == ss2);
}

HEADER_INLINE uint32_t scan_float(const char* ss, float* valp) {
  char* ss2;
  *valp = strtof(ss, &ss2);
  return (ss == ss2);
}

// More restrictive parsing of command-line parameters.
HEADER_INLINE uint32_t scan_doublex(const char* ss, double* valp) {
  char* ss2;
  *valp = strtod(ss, &ss2);
  return (ss == ss2) || (!is_space_or_eoln(ss2[0]));
}

uint32_t scan_posint_cappedx(const char* ss, uint64_t cap, uint32_t* valp);

uint32_t scan_uint_cappedx(const char* ss, uint64_t cap, uint32_t* valp);

uint32_t scan_int_abs_boundedx(const char* ss, uint64_t bound, int32_t* valp);

HEADER_INLINE uint32_t scan_int32x(const char* ss, int32_t* valp) {
  return scan_int_abs_boundedx(ss, 0x7fffffff, valp);
}

HEADER_INLINE uint32_t scan_posint_defcapx(const char* ss, uint32_t* valp) {
  return scan_posint_cappedx(ss, 0x7ffffffe, valp);
}

HEADER_INLINE uint32_t scan_uint_defcapx(const char* ss, uint32_t* valp) {
  return scan_uint_cappedx(ss, 0x7ffffffe, valp);
}

uint32_t scan_posintptrx(const char* ss, uintptr_t* valp);

// __restrict isn't very important for newer x86 processors since loads/stores
// tend to be automatically reordered, but may as well use it properly in
// plink_common.
uint32_t scan_two_doubles(char* ss, double* __restrict val1p, double* __restrict val2p);

int32_t scan_token_ct_len(uintptr_t half_bufsize, FILE* infile, char* buf, uintptr_t* __restrict token_ct_ptr, uintptr_t* __restrict max_token_len_ptr);

int32_t read_tokens(uintptr_t half_bufsize, uintptr_t token_ct, uintptr_t max_token_len, FILE* infile, char* __restrict buf, char* __restrict token_name_buf);

HEADER_INLINE char* memseta(char* target, unsigned char val, uintptr_t ct) {
  memset(target, val, ct);
  return &(target[ct]);
}

HEADER_INLINE char* memcpya(char* __restrict target, const void* __restrict source, uintptr_t ct) {
  memcpy(target, source, ct);
  return &(target[ct]);
}

HEADER_INLINE char* memcpyb(char* __restrict target, const void* __restrict source, uint32_t ct) {
  // Same as memcpya, except the return value is one byte earlier.  Generally
  // used when source is a null-terminated string of known length and we want
  // to copy the null, but sometimes we need to append later.
  memcpy(target, source, ct);
  return &(target[ct - 1]);
}

HEADER_INLINE char* memcpyax(char* __restrict target, const void* __restrict source, uint32_t ct, char extra_char) {
  memcpy(target, source, ct);
  target[ct] = extra_char;
  return &(target[ct + 1]);
}

HEADER_INLINE void memcpyx(char* __restrict target, const void* __restrict source, uint32_t ct, char extra_char) {
  memcpy(target, source, ct);
  target[ct] = extra_char;
}

HEADER_INLINE void memcpyl3(char* __restrict target, const void* __restrict source) {
  // when it's safe to clobber the fourth character, this is faster
  *((uint32_t*)target) = *((const uint32_t*)source);
}

HEADER_INLINE char* memcpyl3a(char* __restrict target, const void* __restrict source) {
  memcpyl3(target, source);
  return &(target[3]);
}

// note that, unlike stpcpy(), this does not copy the null terminator
HEADER_INLINE char* strcpya(char* __restrict target, const void* __restrict source) {
  uintptr_t slen = strlen((char*)source);
  memcpy(target, source, slen);
  return &(target[slen]);
}

HEADER_INLINE char* strcpyax(char* __restrict target, const void* __restrict source, char extra_char) {
  uintptr_t slen = strlen((char*)source);
  memcpy(target, source, slen);
  target[slen] = extra_char;
  return &(target[slen + 1]);
}

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

HEADER_INLINE void fputs_w4(const char* ss, FILE* outfile) {
  // for efficient handling of width-4 allele columns; don't want to call
  // strlen() since that's redundant with fputs
  if (!ss[1]) {
    fputs("   ", outfile);
    putc(ss[0], outfile);
  } else {
    if (!ss[2]) {
      putc(' ', outfile);
      putc(' ', outfile);
    } else if (!ss[3]) {
      putc(' ', outfile);
    }
    fputs(ss, outfile);
  }
}

int32_t gzputs_w4(gzFile gz_outfile, const char* ss);

int32_t get_next_noncomment(FILE* fptr, char** lptr_ptr, uintptr_t* line_idx_ptr);

int32_t get_next_noncomment_excl(const uintptr_t* __restrict marker_exclude, FILE* fptr, char** lptr_ptr, uintptr_t* __restrict line_idx_ptr, uintptr_t* __restrict marker_uidx_ptr);

// assumes we are currently in a token -- UNSAFE OTHERWISE
HEADER_INLINE char* token_endnn(char* sptr) {
  while (!is_space_or_eoln(*(++sptr)));
  return sptr;
}

void get_top_two_ui(const uint32_t* __restrict uint_arr, uintptr_t uia_size, uintptr_t* __restrict top_idx_ptr, uintptr_t* __restrict second_idx_ptr);

uint32_t intlen(int32_t num);

// safer than token_endnn(), since it handles length zero
// "se" = stops at space or eoln character
HEADER_INLINE uintptr_t strlen_se(const char* ss) {
  const char* ss2 = ss;
  while (!is_space_or_eoln(*ss2)) {
    ss2++;
  }
  return (uintptr_t)(ss2 - ss);
}

int32_t strcmp_se(const char* s_read, const char* s_const, uint32_t s_const_len);

char* next_token(char* sptr);

char* next_token_mult(char* sptr, uint32_t ct);

HEADER_INLINE char* next_token_multz(char* sptr, uint32_t ct) {
  // tried replacing this with ternary operator, but that actually seemed to
  // slow things down a bit under gcc 4.2.1 (tail call optimization issue?).
  // todo: recheck this under newer gcc/clang.
  if (ct) {
    return next_token_mult(sptr, ct);
  } else {
    return sptr;
  }
}

uint32_t count_tokens(const char* bufptr);

HEADER_INLINE char* fw_strcpyn(uint32_t min_width, uint32_t source_len, const char* source, char* dest) {
  // right-justified strcpy with known source length
  if (source_len < min_width) {
    memcpy(memseta(dest, 32, min_width - source_len), source, source_len);
    return &(dest[min_width]);
  } else {
    return memcpya(dest, source, source_len);
  }
}

HEADER_INLINE char* fw_strcpy(uint32_t min_width, const char* source, char* dest) {
  return fw_strcpyn(min_width, strlen(source), source, dest);
}

uint32_t count_and_measure_multistr(const char* multistr, uintptr_t* max_slen_ptr);

char* uint32toa(uint32_t uii, char* start);

char* int32toa(int32_t ii, char* start);

// Write exactly four digits (padding with zeroes if necessary); useful for
// e.g. floating point encoders.  uii must not be >= 10^4.
char* uitoa_z4(uint32_t uii, char* start);

char* int64toa(int64_t llii, char* start);

// Minimum field width 4 (padding with spaces on left).
char* uint32toa_w4(uint32_t uii, char* start);

char* uint32toa_w6(uint32_t uii, char* start);

char* uint32toa_w7(uint32_t uii, char* start);

char* uint32toa_w8(uint32_t uii, char* start);

char* uint32toa_w10(uint32_t uii, char* start);

// These limited-precision converters are usually several times as fast as
// grisu2's descendants; and let's not even speak of sprintf.  (I'm guessing
// that the algorithm cannot be made round-trip-safe without throwing away its
// performance advantage, since we currently multiply by numbers like 1.0e256
// which don't have an exact representation.  But these functions are very,
// very good at what they do.)
char* dtoa_e(double dxx, char* start);

char* ftoa_e(float dxx, char* start);

char* dtoa_f_p2(double dxx, char* start);

char* dtoa_f_p3(double dxx, char* start);

char* dtoa_f_w9p6(double dxx, char* start);

char* dtoa_f_w7p4(double dxx, char* start);

HEADER_INLINE void trailing_zeroes_to_spaces(char* start) {
  // removes trailing zeroes
  start--;
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

char* dtoa_f_w9p6_spaced(double dxx, char* start);

char* dtoa_f_w9p6_clipped(double dxx, char* start);

char* dtoa_g(double dxx, char* start);

char* ftoa_g(float dxx, char* start);

HEADER_INLINE char* width_force(uint32_t min_width, char* startp, char* endp) {
  uintptr_t diff = (endp - startp);
  if (diff >= min_width) {
    return endp;
  } else {
    diff = min_width - diff;
    do {
      --endp;
      endp[diff] = *endp;
    } while (endp > startp);
    memset(startp, 32, diff);
    return &(startp[min_width]);
  }
}

// assumes min_width >= 5.
char* dtoa_g_wxp2(double dxx, uint32_t min_width, char* start);

// assumes min_width >= 5.
char* dtoa_g_wxp3(double dxx, uint32_t min_width, char* start);

// only requires min_width to be positive; less than 5 is ok
char* dtoa_g_wxp4(double dxx, uint32_t min_width, char* start);

// only requires min_width to be positive; less than 8 is ok
char* dtoa_g_wxp8(double dxx, uint32_t min_width, char* start);

HEADER_INLINE char* uint32toa_x(uint32_t uii, char extra_char, char* start) {
  char* penult = uint32toa(uii, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* int32toa_x(int32_t ii, char extra_char, char* start) {
  char* penult = int32toa(ii, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* uint32toa_w4x(uint32_t uii, char extra_char, char* start) {
  char* penult = uint32toa_w4(uii, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* uint32toa_w6x(uint32_t uii, char extra_char, char* start) {
  char* penult = uint32toa_w6(uii, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* uint32toa_w7x(uint32_t uii, char extra_char, char* start) {
  char* penult = uint32toa_w7(uii, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* uint32toa_w8x(uint32_t uii, char extra_char, char* start) {
  char* penult = uint32toa_w8(uii, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* uint32toa_w10x(uint32_t uii, char extra_char, char* start) {
  char* penult = uint32toa_w10(uii, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* dtoa_ex(double dxx, char extra_char, char* start) {
  char* penult = dtoa_e(dxx, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* ftoa_ex(float fxx, char extra_char, char* start) {
  char* penult = ftoa_e(fxx, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* dtoa_f_w9p6x(double dxx, char extra_char, char* start) {
  char* penult = dtoa_f_w9p6(dxx, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* dtoa_f_w7p4x(double dxx, char extra_char, char* start) {
  char* penult = dtoa_f_w7p4(dxx, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* dtoa_gx(double dxx, char extra_char, char* start) {
  char* penult = dtoa_g(dxx, start);
  *penult = extra_char;
  return &(penult[1]);
}

/*
HEADER_INLINE char* ftoa_gx(float dxx, char extra_char, char* start) {
  char* penult = ftoa_g(dxx, start);
  *penult = extra_char;
  return &(penult[1]);
}
*/

HEADER_INLINE char* dtoa_g_wxp3x(double dxx, uint32_t min_width, char extra_char, char* start) {
  char* penult = dtoa_g_wxp3(dxx, min_width, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* dtoa_g_wxp4x(double dxx, uint32_t min_width, char extra_char, char* start) {
  char* penult = dtoa_g_wxp4(dxx, min_width, start);
  *penult = extra_char;
  return &(penult[1]);
}

HEADER_INLINE char* dtoa_g_wxp8x(double dxx, uint32_t min_width, char extra_char, char* start) {
  char* penult = dtoa_g_wxp8(dxx, min_width, start);
  *penult = extra_char;
  return &(penult[1]);
}

char* chrom_print_human(uint32_t num, char* buf);

void magic_num(uint32_t divisor, uint64_t* multp, uint32_t* __restrict pre_shiftp, uint32_t* __restrict post_shiftp, uint32_t* __restrict incrp);

HEADER_INLINE uintptr_t tri_coord_no_diag(uintptr_t small_coord, uintptr_t big_coord) {
  // small_coord and big_coord are 0-based indices, small_coord < big_coord
  return ((big_coord * (big_coord - 1)) / 2) + small_coord;
}

HEADER_INLINE uint32_t tri_coord_no_diag_32(uint32_t small_coord, uint32_t big_coord) {
  return ((big_coord * (big_coord - 1)) / 2) + small_coord;
}

// let the compiler worry about the second argument's bit width here
#define SET_BIT(idx, arr) ((arr)[(idx) / BITCT] |= ONELU << ((idx) % BITCT))

#define SET_BIT_DBL(idx, arr) ((arr)[(idx) / BITCT2] |= ONELU << (2 * ((idx) % BITCT2)))

// useful for coercing int32_t loc to unsigned
HEADER_INLINE void set_bit(uint32_t loc, uintptr_t* bitarr) {
  bitarr[loc / BITCT] |= (ONELU << (loc % BITCT));
}

HEADER_INLINE void set_bit_ul(uintptr_t loc, uintptr_t* bitarr) {
  bitarr[loc / BITCT] |= (ONELU << (loc % BITCT));
}

// requires positive len
void fill_bits(uintptr_t loc_start, uintptr_t len, uintptr_t* bitarr);

// requires positive len
void clear_bits(uintptr_t loc_start, uintptr_t len, uintptr_t* bitarr);

#define CLEAR_BIT(idx, arr) ((arr)[(idx) / BITCT] &= ~(ONELU << ((idx) % BITCT)))

#define CLEAR_BIT_DBL(idx, arr) ((arr)[(idx) / BITCT2] &= ~(ONELU << (2 * ((idx) % BITCT2))))

HEADER_INLINE void clear_bit(uint32_t loc, uintptr_t* bitarr) {
  bitarr[loc / BITCT] &= ~(ONELU << (loc % BITCT));
}

HEADER_INLINE void clear_bit_ul(uintptr_t loc, uintptr_t* bitarr) {
  bitarr[loc / BITCT] &= ~(ONELU << (loc % BITCT));
}

#define IS_SET(arr, idx) (((arr)[(idx) / BITCT] >> ((idx) % BITCT)) & 1)

#define IS_SET_DBL(arr, idx) (((arr)[(idx) / BITCT2] >> (2 * ((idx) % BITCT2))) & 1)

// use this instead of IS_SET() for signed 32-bit integers
HEADER_INLINE uint32_t is_set(const uintptr_t* bitarr, uint32_t loc) {
  return (bitarr[loc / BITCT] >> (loc % BITCT)) & 1;
}

HEADER_INLINE uint32_t is_set_ul(const uintptr_t* bitarr, uintptr_t loc) {
  return (bitarr[loc / BITCT] >> (loc % BITCT)) & 1;
}

#define IS_NONNULL_AND_SET(arr, idx) ((arr) && IS_SET(arr, idx))

uint32_t next_unset_unsafe(const uintptr_t* bitarr, uint32_t loc);

HEADER_INLINE void next_unset_unsafe_ck(const uintptr_t* __restrict bitarr, uint32_t* __restrict loc_ptr) {
  if (IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_unset_unsafe(bitarr, *loc_ptr);
  }
}

#ifdef __LP64__
uintptr_t next_unset_ul_unsafe(const uintptr_t* bitarr, uintptr_t loc);
#else
HEADER_INLINE uintptr_t next_unset_ul_unsafe(const uintptr_t* bitarr, uintptr_t loc) {
  return (uintptr_t)next_unset_unsafe(bitarr, loc);
}
#endif

HEADER_INLINE void next_unset_ul_unsafe_ck(const uintptr_t* __restrict bitarr, uintptr_t* __restrict loc_ptr) {
  if (IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_unset_ul_unsafe(bitarr, *loc_ptr);
  }
}

uint32_t next_unset(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil);

HEADER_INLINE void next_unset_ck(const uintptr_t* __restrict bitarr, uint32_t ceil, uint32_t* __restrict loc_ptr) {
  if (IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_unset(bitarr, *loc_ptr, ceil);
  }
}

#ifdef __LP64__
uintptr_t next_unset_ul(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil);
#else
HEADER_INLINE uintptr_t next_unset_ul(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil) {
  return (uintptr_t)next_unset(bitarr, loc, ceil);
}
#endif

HEADER_INLINE void next_unset_ul_ck(const uintptr_t* __restrict bitarr, uintptr_t ceil, uintptr_t* __restrict loc_ptr) {
  if (IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_unset_ul(bitarr, *loc_ptr, ceil);
  }
}

uint32_t next_set_unsafe(const uintptr_t* bitarr, uint32_t loc);

HEADER_INLINE void next_set_unsafe_ck(const uintptr_t* __restrict bitarr, uint32_t* __restrict loc_ptr) {
  if (!IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_set_unsafe(bitarr, *loc_ptr);
  }
}

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

uint32_t next_set(const uintptr_t* bitarr, uint32_t loc, uint32_t ceil);

HEADER_INLINE void next_set_ck(const uintptr_t* __restrict bitarr, uint32_t ceil, uint32_t* __restrict loc_ptr) {
  if (!IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_set(bitarr, *loc_ptr, ceil);
  }
}

#ifdef __LP64__
uintptr_t next_set_ul(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil);
#else
HEADER_INLINE uintptr_t next_set_ul(const uintptr_t* bitarr, uintptr_t loc, uintptr_t ceil) {
  return (uintptr_t)next_set(bitarr, loc, ceil);
}
#endif

HEADER_INLINE void next_set_ul_ck(const uintptr_t* __restrict bitarr, uintptr_t ceil, uintptr_t* loc_ptr) {
  if (!IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = next_set_ul(bitarr, *loc_ptr, ceil);
  }
}

int32_t last_set_bit(const uintptr_t* bitarr, uint32_t word_ct);

// note different interface from last_set_bit()
// int32_t last_clear_bit(uintptr_t* bitarr, uint32_t ceil);

// unlike the next_[un]set family, this always returns a STRICTLY earlier
// position
uint32_t prev_unset_unsafe(const uintptr_t* bitarr, uint32_t loc);

// uint32_t prev_unset(uintptr_t* bitarr, uint32_t loc, uint32_t floor);

HEADER_INLINE void prev_unset_unsafe_ck(const uintptr_t* bitarr, uint32_t* loc_ptr) {
  *loc_ptr -= 1;
  if (IS_SET(bitarr, *loc_ptr)) {
    *loc_ptr = prev_unset_unsafe(bitarr, *loc_ptr);
  }
}

// These functions seem to optimize better than memset(arr, 0, x) under OS X
// <10.9's gcc, and they should be equivalent for later versions (looks like
// memcpy/memset were redone in gcc 4.3).
HEADER_INLINE void fill_ulong_zero(size_t size, uintptr_t* ularr) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *ularr++ = 0;
  }
}

#ifdef __LP64__
HEADER_INLINE void fill_ull_zero(size_t size, uint64_t* ullarr) {
  fill_ulong_zero(size, (uintptr_t*)ullarr);
}

// double v indicates that size is a vector count, not a word count.
HEADER_INLINE void fill_vvec_zero(size_t size, VECITYPE* vvec) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *vvec++ = _mm_setzero_si128();
  }
}
#else
HEADER_INLINE void fill_ull_zero(size_t size, uint64_t* ullarr) {
  fill_ulong_zero(size * 2, (uintptr_t*)ullarr);
}
#endif

HEADER_INLINE void fill_ulong_one(size_t size, uintptr_t* ularr) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *ularr++ = ~ZEROLU;
  }
}

#ifdef __LP64__
HEADER_INLINE void fill_ull_one(size_t size, uint64_t* ullarr) {
  fill_ulong_one(size, (uintptr_t*)ullarr);
}
#else
HEADER_INLINE void fill_ull_one(size_t size, uint64_t* ullarr) {
  fill_ulong_one(size * 2, (uintptr_t*)ullarr);
}
#endif

HEADER_INLINE void fill_int_zero(size_t size, int32_t* iarr) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *iarr++ = 0;
  }
}

HEADER_INLINE void fill_int_one(size_t size, int32_t* iarr) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *iarr++ = -1;
  }
}

HEADER_INLINE void fill_uint_zero(size_t size, uint32_t* uiarr) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *uiarr++ = 0;
  }
}

HEADER_INLINE void fill_uint_one(size_t size, uint32_t* uiarr) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *uiarr++ = ~0U;
  }
}

HEADER_INLINE void fill_float_zero(size_t size, float* farr) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *farr++ = 0.0;
  }
}

HEADER_INLINE void fill_double_zero(size_t size, double* darr) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *darr++ = 0.0;
  }
}


int32_t bigstack_calloc_uc(uintptr_t ct, unsigned char** ucp_ptr);

int32_t bigstack_calloc_d(uintptr_t ct, double** dp_ptr);

int32_t bigstack_calloc_f(uintptr_t ct, float** fp_ptr);

int32_t bigstack_calloc_ui(uintptr_t ct, uint32_t** uip_ptr);

int32_t bigstack_calloc_ul(uintptr_t ct, uintptr_t** ulp_ptr);

int32_t bigstack_calloc_ull(uintptr_t ct, uint64_t** ullp_ptr);

HEADER_INLINE int32_t bigstack_calloc_c(uintptr_t ct, char** cp_ptr) {
  return bigstack_calloc_uc(ct, (unsigned char**)cp_ptr);
}

HEADER_INLINE int32_t bigstack_calloc_i(uintptr_t ct, int32_t** ip_ptr) {
  return bigstack_calloc_ui(ct, (uint32_t**)ip_ptr);
}

HEADER_INLINE int32_t bigstack_calloc_ll(uintptr_t ct, int64_t** llp_ptr) {
  return bigstack_calloc_ull(ct, (uint64_t**)llp_ptr);
}

int32_t bigstack_end_calloc_uc(uintptr_t ct, unsigned char** ucp_ptr);

int32_t bigstack_end_calloc_d(uintptr_t ct, double** dp_ptr);

int32_t bigstack_end_calloc_f(uintptr_t ct, float** fp_ptr);

int32_t bigstack_end_calloc_ui(uintptr_t ct, uint32_t** uip_ptr);

int32_t bigstack_end_calloc_ul(uintptr_t ct, uintptr_t** ulp_ptr);

int32_t bigstack_end_calloc_ull(uintptr_t ct, uint64_t** ullp_ptr);

HEADER_INLINE int32_t bigstack_end_calloc_c(uintptr_t ct, char** cp_ptr) {
  return bigstack_end_calloc_uc(ct, (unsigned char**)cp_ptr);
}

HEADER_INLINE int32_t bigstack_end_calloc_i(uintptr_t ct, int32_t** ip_ptr) {
  return bigstack_end_calloc_ui(ct, (uint32_t**)ip_ptr);
}

HEADER_INLINE int32_t bigstack_end_calloc_ll(uintptr_t ct, int64_t** llp_ptr) {
  return bigstack_end_calloc_ull(ct, (uint64_t**)llp_ptr);
}


uint32_t murmurhash3_32(const void* key, uint32_t len);

HEADER_INLINE uint32_t hashval2(const char* idstr, uint32_t idlen) {
  return murmurhash3_32(idstr, idlen) % HASHSIZE;
}

uintptr_t geqprime(uintptr_t floor);

HEADER_INLINE uint32_t get_id_htable_size(uintptr_t item_ct) {
  if (item_ct < 32761) {
    return 65521;
  } else {
    return geqprime(item_ct * 2 + 1);
  }
}

int32_t populate_id_htable(uintptr_t unfiltered_ct, const uintptr_t* exclude_arr, uintptr_t item_ct, const char* item_ids, uintptr_t max_id_len, uint32_t store_dups, uint32_t id_htable_size, uint32_t* id_htable);

HEADER_INLINE int32_t alloc_and_populate_id_htable(uintptr_t unfiltered_ct, const uintptr_t* exclude_arr, uintptr_t item_ct, const char* item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t* id_htable_size_ptr, uint32_t** id_htable_ptr) {
  uint32_t id_htable_size = get_id_htable_size(item_ct);
  if (bigstack_alloc_ui(id_htable_size, id_htable_ptr)) {
    return RET_NOMEM;
  }
  *id_htable_size_ptr = id_htable_size;
  return populate_id_htable(unfiltered_ct, exclude_arr, item_ct, item_ids, max_id_len, allow_dups, id_htable_size, *id_htable_ptr);
}

uint32_t id_htable_find(const char* id_buf, uintptr_t cur_id_len, const uint32_t* id_htable, uint32_t id_htable_size, const char* item_ids, uintptr_t max_id_len);

void fill_idx_to_uidx(const uintptr_t* exclude_arr, uintptr_t unfiltered_item_ct, uintptr_t item_ct, uint32_t* idx_to_uidx);

void fill_idx_to_uidx_incl(const uintptr_t* include_arr, uintptr_t unfiltered_item_ct, uintptr_t item_ct, uint32_t* idx_to_uidx);

void fill_uidx_to_idx(const uintptr_t* exclude_arr, uint32_t unfiltered_item_ct, uint32_t item_ct, uint32_t* uidx_to_idx);

void fill_uidx_to_idx_incl(const uintptr_t* include_arr, uint32_t unfiltered_item_ct, uint32_t item_ct, uint32_t* uidx_to_idx);

void fill_midx_to_idx(const uintptr_t* exclude_arr_orig, const uintptr_t* exclude_arr, uint32_t item_ct, uint32_t* midx_to_idx);


// "quaterarr" refers to a packed group of base-4 (2-bit) elements, analogous
// to "bitarr".  (Based on "quaternary", not "quarter".)  "quatervec"
// indicates that vector-alignment is also required.
void fill_quatervec_55(uint32_t ct, uintptr_t* quatervec);

// Used to unpack e.g. unfiltered sex_male to a filtered quaterarr usable as a
// raw input bitmask.
// Assumes output_quaterarr is sized to a multiple of 16 bytes.
void quaterarr_collapse_init(const uintptr_t* __restrict unfiltered_bitarr, uint32_t unfiltered_ct, const uintptr_t* __restrict filter_bitarr, uint32_t filtered_ct, uintptr_t* __restrict output_quaterarr);

void quaterarr_collapse_init_exclude(const uintptr_t* __restrict unfiltered_bitarr, uint32_t unfiltered_ct, const uintptr_t* __restrict filter_exclude_bitarr, uint32_t filtered_ct, uintptr_t* __restrict output_quaterarr);

uint32_t alloc_collapsed_haploid_filters(const uintptr_t* __restrict sample_bitarr, const uintptr_t* __restrict sex_male, uint32_t unfiltered_sample_ct, uint32_t sample_ct, uint32_t hh_exists, uint32_t is_include, uintptr_t** sample_include_quatervec_ptr, uintptr_t** sample_male_include_quatervec_ptr);

HEADER_INLINE void free_cond(void* memptr) {
  if (memptr) {
    free(memptr);
  }
}

HEADER_INLINE uint32_t realnum(double dd) {
  return (dd == dd) && (dd != INFINITY) && (dd != -INFINITY);
}

HEADER_INLINE double get_maf(double allele_freq) {
  return (allele_freq <= 0.5)? allele_freq : (1.0 - allele_freq);
}

HEADER_INLINE int32_t filename_exists(const char* __restrict fname_append, char* fname, char* fname_end) {
#ifdef _WIN32
  DWORD file_attr;
  strcpy(fname_end, fname_append);
  file_attr = GetFileAttributes(fname);
  return (file_attr != 0xffffffffU);
#else
  struct stat st;
  strcpy(fname_end, fname_append);
  return (stat(fname, &st) == 0);
#endif
}

void sample_delim_convert(uintptr_t unfiltered_sample_ct, const uintptr_t* sample_exclude, uint32_t sample_ct, uintptr_t max_sample_id_len, char oldc, char newc, char* sample_ids);

void get_set_wrange_align(const uintptr_t* __restrict bitarr, uintptr_t word_ct, uintptr_t* __restrict firstw_ptr, uintptr_t* __restrict wlen_ptr);

// for hash tables where maximum ID string length is not known in advance.
uint32_t unklen_id_htable_find(const char* cur_id, const char* const* item_ids, const uint32_t* id_htable, uint32_t hashval, uint32_t id_htable_size);

// okay, time to provide O(c log c)-time instead of O(c^2)-time initialization
// (c = # of chromosomes/contigs).
#define MAX_POSSIBLE_CHROM 65280

// get_id_htable_size(MAX_POSSIBLE_CHROM) (use constexpr once sufficient
// compiler support is available)
#define CHROM_NAME_HTABLE_SIZE 130579

// assumes MAX_POSSIBLE_CHROM is a multiple of 64, otherwise add round-up
#define CHROM_MASK_WORDS (MAX_POSSIBLE_CHROM / BITCT)

// (note that n+1, n+2, n+3, and n+4 are reserved for X/Y/XY/MT)
#define MAX_CHROM_TEXTNUM 95

// get_chrom_code_raw() needs to be modified if this changes
#define MAX_CHROM_TEXTNUM_SLEN 2

#define X_OFFSET 0
#define Y_OFFSET 1
#define XY_OFFSET 2
#define MT_OFFSET 3
#define XYMT_OFFSET_CT 4

#define CHROM_X (MAX_POSSIBLE_CHROM + X_OFFSET)
#define CHROM_Y (MAX_POSSIBLE_CHROM + Y_OFFSET)
#define CHROM_XY (MAX_POSSIBLE_CHROM + XY_OFFSET)
#define CHROM_MT (MAX_POSSIBLE_CHROM + MT_OFFSET)

#ifdef __LP64__
  // dog requires 42 bits, and other species require less
  #define CHROM_MASK_INITIAL_WORDS 1
#else
  #define CHROM_MASK_INITIAL_WORDS 2
#endif

typedef struct {
  // Main dynamic block intended to be allocated as a single aligned block of
  // memory on the heap freeable with vecaligned_free(), with chrom_mask at the
  // base.

  uintptr_t* chrom_mask; // which chromosomes aren't known to be absent?
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
  uint32_t* chrom_file_order;

  // if the second chromosome in the dataset is chr5, chrom_file_order[1] == 5,
  // the raw variant indexes for chr5 are in [chrom_fo_vidx_start[1],
  // chrom_fo_vidx_start[2]). and chrom_idx_to_foidx[5] == 1.
  uint32_t* chrom_fo_vidx_start;
  uint32_t* chrom_idx_to_foidx;

  // --allow-extra-chr support
  char** nonstd_names;
  uint32_t* nonstd_id_htable;
  // end main dynamic block

  uint32_t chrom_ct; // number of distinct chromosomes/contigs
  uint32_t species;

  int32_t xymt_codes[XYMT_OFFSET_CT]; // x, y, xy, mt
  uint32_t max_code;

  uint32_t autosome_ct;

  // yet more --allow-extra-chr support
  uint32_t zero_extra_chroms;
  uint32_t name_ct;
  Ll_str* incl_excl_name_stack;
  uint32_t is_include_stack;
  uint32_t output_encoding;
} Chrom_info;

extern const char* g_species_singular;
extern const char* g_species_plural;

int32_t init_chrom_info(Chrom_info* chrom_info_ptr);

void init_species(uint32_t species_code, Chrom_info* chrom_info_ptr);

void init_default_chrom_mask(Chrom_info* chrom_info_ptr);

HEADER_INLINE int32_t init_chrom_info_human(Chrom_info* chrom_info_ptr) {
  // convenience wrapper
  if (init_chrom_info(chrom_info_ptr)) {
    return RET_NOMEM;
  }
  init_species(SPECIES_HUMAN, chrom_info_ptr);
  init_default_chrom_mask(chrom_info_ptr);
  return 0;
}

void forget_extra_chrom_names(uint32_t reinitialize, Chrom_info* chrom_info_ptr);

// in the usual case where the number of chromosomes/contigs is much less than
// MAX_POSSIBLE_CHROM, this reduces chrom_info's memory consumption and
// improves locality.
int32_t finalize_chrom_info(Chrom_info* chrom_info_ptr);

void cleanup_chrom_info(Chrom_info* chrom_info_ptr);

HEADER_INLINE const char* species_str(uintptr_t ct) {
  return (ct == ONELU)? g_species_singular : g_species_plural;
}

#define CHR_OUTPUT_PREFIX 1
#define CHR_OUTPUT_M 2
#define CHR_OUTPUT_MT 4
#define CHR_OUTPUT_0M 8

HEADER_INLINE uint32_t are_all_words_zero(const uintptr_t* word_arr, uintptr_t word_ct) {
  while (word_ct--) {
    if (*word_arr++) {
      return 0;
    }
  }
  return 1;
}

char* chrom_name_write(const Chrom_info* chrom_info_ptr, uint32_t chrom_idx, char* buf);

char* chrom_name_buf5w4write(const Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uint32_t* chrom_name_len_ptr, char* buf5);

uint32_t get_max_chrom_slen(const Chrom_info* chrom_info_ptr);

uint32_t haploid_chrom_present(const Chrom_info* chrom_info_ptr);

// does not require null-termination
// only handles 1-99, X, Y, XY, MT, and "chr" prefix
int32_t get_chrom_code_raw(const char* sptr);

// now requires null-termination
// now returns -1 when --allow-extra-chr may be ok, and -2 on total fail
int32_t get_chrom_code(const char* chrom_name, const Chrom_info* chrom_info_ptr, uint32_t name_slen);

// when the chromosome name isn't null-terminated, but we want to preserve the
// character there
// requires chrom_name[name_slen] to be mutable
int32_t get_chrom_code_counted(const Chrom_info* chrom_info_ptr, uint32_t name_slen, char* chrom_name);

// when it's okay to just replace the terminating space/tab with a \0
HEADER_INLINE int32_t get_chrom_code_destructive(const Chrom_info* chrom_info_ptr, char* chrom_name) {
  char* chrom_token_end = token_endnn(chrom_name);
  *chrom_token_end = '\0';
  return get_chrom_code(chrom_name, chrom_info_ptr, (uintptr_t)(chrom_token_end - chrom_name));
}

uint32_t get_variant_chrom_fo_idx(const Chrom_info* chrom_info_ptr, uintptr_t variant_uidx);

HEADER_INLINE uint32_t get_variant_chrom(const Chrom_info* chrom_info_ptr, uintptr_t variant_uidx) {
  return chrom_info_ptr->chrom_file_order[get_variant_chrom_fo_idx(chrom_info_ptr, variant_uidx)];
}


// these assume the chromosome is present in the dataset
HEADER_INLINE uint32_t get_chrom_start_vidx(const Chrom_info* chrom_info_ptr, uint32_t chrom_idx) {
  return chrom_info_ptr->chrom_fo_vidx_start[chrom_info_ptr->chrom_idx_to_foidx[chrom_idx]];
}

HEADER_INLINE uint32_t get_chrom_end_vidx(const Chrom_info* chrom_info_ptr, uint32_t chrom_idx) {
  return chrom_info_ptr->chrom_fo_vidx_start[chrom_info_ptr->chrom_idx_to_foidx[chrom_idx] + 1];
}

// now assumes chrom_name is null-terminated
int32_t try_to_add_chrom_name(const char* chrom_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t allow_extra_chroms, int32_t* chrom_idx_ptr, Chrom_info* chrom_info_ptr);

HEADER_INLINE int32_t get_or_add_chrom_code(const char* chrom_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t allow_extra_chroms, Chrom_info* chrom_info_ptr, int32_t* chrom_idx_ptr) {
  *chrom_idx_ptr = get_chrom_code(chrom_name, chrom_info_ptr, name_slen);
  if (*chrom_idx_ptr >= 0) {
    return 0;
  }
  return try_to_add_chrom_name(chrom_name, file_descrip, line_idx, name_slen, allow_extra_chroms, chrom_idx_ptr, chrom_info_ptr);
}

HEADER_INLINE int32_t get_or_add_chrom_code_destructive(const char* file_descrip, uintptr_t line_idx, uint32_t allow_extra_chroms, char* chrom_name, char* chrom_name_end, Chrom_info* chrom_info_ptr, int32_t* chrom_idx_ptr) {
  *chrom_name_end = '\0';
  return get_or_add_chrom_code(chrom_name, file_descrip, line_idx, (uintptr_t)(chrom_name_end - chrom_name), allow_extra_chroms, chrom_info_ptr, chrom_idx_ptr);
}

// newval does not need to be null-terminated
// assumes *allele_ptr is not initialized
// make last parameter const char** later
uint32_t allele_set(const char* newval, uint32_t allele_slen, char** allele_ptr);

// *allele_ptr must be initialized; frees *allele_ptr if necessary
uint32_t allele_reset(const char* newval, uint32_t allele_slen, char** allele_ptr);

void cleanup_allele_storage(uint32_t max_allele_slen, uintptr_t allele_storage_entry_ct, char** allele_storage);

// needed by fixed --merge-equal-pos implementation, which takes more liberties
// with allele_storage[]
void cleanup_allele_storage2(uintptr_t allele_storage_entry_ct, char** allele_storage);

// no need for this; code is simpler if we just create a copy of marker_exclude
// with all non-autosomal loci removed
/*
HEADER_INLINE uintptr_t next_autosomal_unsafe(uintptr_t* marker_exclude, uintptr_t marker_uidx, Chrom_info* chrom_info_ptr, uint32_t* chrom_end_ptr, uint32_t* chrom_fo_idx_ptr) {
  // assumes we are at an autosomal marker if marker_uidx < *chrom_end_ptr
  next_unset_ul_unsafe_ck(marker_exclude, &marker_uidx);
  if (marker_uidx < (*chrom_end_ptr)) {
    return marker_uidx;
  }
  uintptr_t* haploid_mask = chrom_info_ptr->haploid_mask;
  uint32_t chrom_idx;
  while (1) {
    do {
      *chrom_fo_idx_ptr += 1;
      *chrom_end_ptr = chrom_info_ptr->chrom_file_order_marker_idx[(*chrom_fo_idx_ptr) + 1];
    } while (marker_uidx >= (*chrom_end_ptr));
    chrom_idx = chrom_info_ptr->chrom_file_order[*chrom_fo_idx_ptr];
    if (!IS_SET(haploid_mask, chrom_idx)) {
      return marker_uidx;
    }
    marker_uidx = next_unset_ul_unsafe(marker_exclude, *chrom_end_ptr);
  }
}
*/

void refresh_chrom_info(const Chrom_info* chrom_info_ptr, uintptr_t marker_uidx, uint32_t* __restrict chrom_end_ptr, uint32_t* __restrict chrom_fo_idx_ptr, uint32_t* __restrict is_x_ptr, uint32_t* __restrict is_y_ptr, uint32_t* __restrict is_mt_ptr, uint32_t* __restrict is_haploid_ptr);

int32_t single_chrom_start(const Chrom_info* chrom_info_ptr, const uintptr_t* marker_exclude, uint32_t unfiltered_marker_ct);

double get_dmedian(const double* sorted_arr, uintptr_t len);

double destructive_get_dmedian(uintptr_t len, double* unsorted_arr);

int32_t strcmp_casted(const void* s1, const void* s2);

int32_t strcmp_natural(const void* s1, const void* s2);

int32_t strcmp_deref(const void* s1, const void* s2);

int32_t strcmp_natural_deref(const void* s1, const void* s2);

int32_t get_uidx_from_unsorted(const char* idstr, const uintptr_t* exclude_arr, uint32_t id_ct, const char* unsorted_ids, uintptr_t max_id_len);

// sorted_ids contents not changed, but not worth the trouble of returning a
// const char*
char* scan_for_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len);

char* scan_for_duplicate_or_overlap_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len, const char* sorted_nonoverlap_ids, uintptr_t nonoverlap_id_ct, uintptr_t max_nonoverlap_id_len);

int32_t eval_affection(const char* bufptr, double missing_phenod);

uint32_t triangle_divide(int64_t cur_prod, int32_t modif);

void triangle_fill(uint32_t ct, uint32_t pieces, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t start, uint32_t align, uint32_t* target_arr);

int32_t relationship_req(uint64_t calculation_type);

int32_t distance_req(const char* read_dists_fname, uint64_t calculation_type);

int32_t double_cmp(const void* aa, const void* bb);

int32_t double_cmp_decr(const void* aa, const void* bb);

int32_t double_cmp_deref(const void* aa, const void* bb);

int32_t char_cmp_deref(const void* aa, const void* bb);

int32_t intcmp(const void* aa, const void* bb);

int32_t uintcmp(const void* aa, const void* bb);

#ifndef __cplusplus
int32_t intcmp2(const void* aa, const void* bb);
#endif

int32_t intcmp3_decr(const void* aa, const void* bb);

#ifndef __cplusplus
int32_t llcmp(const void* aa, const void* bb);
#endif

void qsort_ext2(char* main_arr, uintptr_t arr_length, uintptr_t item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, uintptr_t secondary_item_len, char* proxy_arr, uintptr_t proxy_len);

int32_t qsort_ext(char* main_arr, uintptr_t arr_length, uintptr_t item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, intptr_t secondary_item_len);

int32_t sort_item_ids_noalloc(uintptr_t unfiltered_ct, const uintptr_t* exclude_arr, uintptr_t item_ct, const char* __restrict item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t collapse_idxs, int(* comparator_deref)(const void*, const void*), char* __restrict sorted_ids, uint32_t* id_map);

int32_t sort_item_ids(uintptr_t unfiltered_ct, const uintptr_t* exclude_arr, uintptr_t exclude_ct, const char* __restrict item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t collapse_idxs, int(* comparator_deref)(const void*, const void*), char** sorted_ids_ptr, uint32_t** id_map_ptr);

uint32_t uint32arr_greater_than(const uint32_t* sorted_uint32_arr, uint32_t arr_length, uint32_t uii);

uint32_t int32arr_greater_than(const int32_t* sorted_int32_arr, uint32_t arr_length, int32_t ii);

uintptr_t uint64arr_greater_than(const uint64_t* sorted_uint64_arr, uintptr_t arr_length, uint64_t ullii);

uintptr_t doublearr_greater_than(const double* sorted_dbl_arr, uintptr_t arr_length, double dxx);

uintptr_t nonincr_doublearr_leq_stride(const double* nonincr_dbl_arr, uintptr_t arr_length, uintptr_t stride, double dxx);

int32_t bsearch_str(const char* id_buf, uintptr_t cur_id_len, const char* lptr, uintptr_t max_id_len, uintptr_t end_idx);

HEADER_INLINE int32_t bsearch_str_nl(const char* id_buf, const char* lptr, uintptr_t max_id_len, intptr_t end_idx) {
  return bsearch_str(id_buf, strlen(id_buf), lptr, max_id_len, end_idx);
}

int32_t bsearch_str_natural(const char* id_buf, const char* lptr, uintptr_t max_id_len, uintptr_t end_idx);

uintptr_t bsearch_str_lb(const char* id_buf, uintptr_t cur_id_len, const char* lptr, uintptr_t max_id_len, uintptr_t end_idx);

uint32_t bsearch_read_fam_indiv(char* __restrict read_ptr, const char* __restrict lptr, uintptr_t max_id_len, uintptr_t filter_line_ct, char** read_pp_new, int32_t* retval_ptr, char* __restrict id_buf);

void bsearch_fam(const char* __restrict fam_id, const char* __restrict lptr, uintptr_t max_id_len, uint32_t filter_line_ct, uint32_t* __restrict first_idx_ptr, uint32_t* __restrict last_idx_ptr, char* __restrict id_buf);

// These ensure the trailing bits are zeroed out.
void bitarr_invert(uintptr_t bit_ct, uintptr_t* bitarr);

void bitarr_invert_copy(const uintptr_t* input_bitarr, uintptr_t bit_ct, uintptr_t* output_bitarr);


// "bitvec" indicates that word count is used instead of vector count.
void bitvec_and(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

void bitvec_andnot(const uintptr_t* __restrict exclude_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

void bitvec_andnot_reversed_args(const uintptr_t* __restrict include_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

void bitvec_or(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec);

void bitvec_ornot(const uintptr_t* __restrict inverted_or_bitvec, uintptr_t word_ct, uintptr_t* main_bitvec);

void bitvec_xor(const uintptr_t* __restrict arg_bitvec, uintptr_t word_ct, uintptr_t* __restrict main_bitvec);

HEADER_INLINE uint32_t popcount2_long(uintptr_t val) {
#ifdef __LP64__
  val = (val & 0x3333333333333333LLU) + ((val >> 2) & 0x3333333333333333LLU);
  return (((val + (val >> 4)) & 0x0f0f0f0f0f0f0f0fLLU) * 0x0101010101010101LLU) >> 56;
#else
  val = (val & 0x33333333) + ((val >> 2) & 0x33333333);
  return (((val + (val >> 4)) & 0x0f0f0f0f) * 0x01010101) >> 24;
#endif
}

HEADER_INLINE uint32_t popcount_long(uintptr_t val) {
  // the simple version, good enough for all non-time-critical stuff
  return popcount2_long(val - ((val >> 1) & FIVEMASK));
}

uint32_t is_monomorphic_a2(const uintptr_t* geno_arr, uint32_t sample_ct);

uint32_t is_monomorphic(const uintptr_t* geno_arr, uint32_t sample_ct);

// same as is_monomorphic, except it also flags the all-heterozygote case
uint32_t less_than_two_genotypes(const uintptr_t* geno_arr, uint32_t sample_ct);

// uint32_t has_three_genotypes(uintptr_t* lptr, uint32_t sample_ct);

uintptr_t popcount_longs(const uintptr_t* lptr, uintptr_t word_ct);

#ifdef __LP64__
HEADER_INLINE uintptr_t popcount_longs_nzbase(const uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
  uintptr_t prefix_ct = 0;
  if (start_idx & 1) {
    if (end_idx == start_idx) {
      return 0;
    }
    prefix_ct = popcount_long(lptr[start_idx++]);
  }
  return prefix_ct + popcount_longs(&(lptr[start_idx]), end_idx - start_idx);
}
#else
HEADER_INLINE uintptr_t popcount_longs_nzbase(const uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
  return popcount_longs(&(lptr[start_idx]), end_idx - start_idx);
}
#endif

uintptr_t popcount2_longs(const uintptr_t* lptr, uintptr_t word_ct);

#define popcount01_longs popcount2_longs

uintptr_t popcount_bit_idx(const uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx);

uint32_t chrom_window_max(const uint32_t* marker_pos, const uintptr_t* marker_exclude, const Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max);

uint32_t window_back(const uint32_t* __restrict marker_pos, const double* __restrict marker_cms, const uintptr_t* marker_exclude, uint32_t marker_uidx_min, uint32_t marker_uidx_start, uint32_t count_max, uint32_t bp_max, double cm_max, uint32_t* __restrict window_trail_ct_ptr);

uint32_t window_forward(const uint32_t* __restrict marker_pos, const double* __restrict marker_cms, const uintptr_t* marker_exclude, uint32_t marker_uidx_start, uint32_t marker_uidx_last, uint32_t count_max, uint32_t bp_max, double cm_max, uint32_t* __restrict window_lead_ct_ptr);

uintptr_t jump_forward_unset_unsafe(const uintptr_t* bitarr, uintptr_t cur_pos, uintptr_t forward_ct);

HEADER_INLINE uintptr_t popcount_chars(const uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
  return popcount_bit_idx(lptr, start_idx * 8, end_idx * 8);
}

// end_idx is in word, not bit, units
uintptr_t popcount_longs_exclude(const uintptr_t* __restrict lptr, const uintptr_t* __restrict exclude_arr, uintptr_t end_idx);

uintptr_t popcount_longs_intersect(const uintptr_t* __restrict lptr1, const uintptr_t* __restrict lptr2, uintptr_t word_ct);

void vertical_bitct_subtract(const uintptr_t* bitarr, uint32_t item_ct, uint32_t* sum_arr);

#ifdef __LP64__
void count_2freq_dbl_960b(const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end, const VECITYPE* __restrict mask1vp, const VECITYPE* __restrict mask2vp, uint32_t* __restrict ct1abp, uint32_t* __restrict ct1cp, uint32_t* __restrict ct2abp, uint32_t* __restrict ct2cp);

void count_3freq_1920b(const VECITYPE* geno_vvec, const VECITYPE* geno_vvec_end, const VECITYPE* __restrict maskvp, uint32_t* __restrict ctap, uint32_t* __restrict ctbp, uint32_t* __restrict ctcp);
#else
void count_2freq_dbl_24b(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict mask1p, const uintptr_t* __restrict mask2p, uint32_t* __restrict ct1abp, uint32_t* __restrict ct1cp, uint32_t* __restrict ct2abp, uint32_t* __restrict ct2cp);

void count_3freq_48b(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict maskp, uint32_t* __restrict ctap, uint32_t* __restrict ctbp, uint32_t* __restrict ctcp);
#endif

void genovec_set_freq(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict include_quatervec, uintptr_t sample_ctl2, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp);

void genovec_set_freq_x(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict include_quatervec, const uintptr_t* __restrict male_quatervec, uintptr_t sample_ctl2, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp);

void genovec_set_freq_y(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict include_quatervec, const uintptr_t* __restrict nonmale_quatervec, uintptr_t sample_ctl2, uint32_t* __restrict set_ctp, uint32_t* __restrict missing_ctp);

void genovec_3freq(const uintptr_t* __restrict geno_vec, const uintptr_t* __restrict include_quatervec, uintptr_t sample_ctl2, uint32_t* __restrict missing_ctp, uint32_t* __restrict het_ctp, uint32_t* __restrict homset_ctp);

uintptr_t count_01(const uintptr_t* quatervec, uintptr_t word_ct);

HEADER_INLINE void zero_trailing_bits(uintptr_t unfiltered_ct, uintptr_t* bitarr) {
  uintptr_t trail_ct = unfiltered_ct & (BITCT - 1);
  if (trail_ct) {
    bitarr[unfiltered_ct / BITCT] &= (ONELU << trail_ct) - ONELU;
  }
}

void fill_all_bits(uintptr_t ct, uintptr_t* bitarr);

uint32_t numeric_range_list_to_bitarr(const Range_list* range_list_ptr, uint32_t item_ct, uint32_t offset, uint32_t ignore_overflow, uintptr_t* bitarr);

int32_t string_range_list_to_bitarr(char* header_line, uint32_t item_ct, uint32_t fixed_len, const Range_list* range_list_ptr, const char* __restrict sorted_ids, const uint32_t* __restrict id_map, const char* __restrict range_list_flag, const char* __restrict file_descrip, uintptr_t* bitarr, int32_t* __restrict seen_idxs);

int32_t string_range_list_to_bitarr_alloc(char* header_line, uint32_t item_ct, uint32_t fixed_len, const Range_list* range_list_ptr, const char* __restrict range_list_flag, const char* __restrict file_descrip, uintptr_t** bitarr_ptr);

int32_t string_range_list_to_bitarr2(const char* __restrict sorted_ids, const uint32_t* id_map, uintptr_t item_ct, uintptr_t max_id_len, const Range_list* __restrict range_list_ptr, const char* __restrict range_list_flag, uintptr_t* bitarr_excl);

HEADER_INLINE uint32_t count_chrom_markers(const Chrom_info* chrom_info_ptr, const uintptr_t* marker_exclude, uint32_t chrom_idx) {
  if (!is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
    return 0;
  }
  const uint32_t chrom_fo_idx = chrom_info_ptr->chrom_idx_to_foidx[chrom_idx];
  const uint32_t min_idx = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx];
  const uint32_t max_idx = chrom_info_ptr->chrom_fo_vidx_start[chrom_fo_idx + 1];
  return (max_idx - min_idx) - ((uint32_t)popcount_bit_idx(marker_exclude, min_idx, max_idx));
}

uint32_t count_non_autosomal_markers(const Chrom_info* chrom_info_ptr, const uintptr_t* marker_exclude, uint32_t count_x, uint32_t count_mt);

int32_t conditional_allocate_non_autosomal_markers(const Chrom_info* chrom_info_ptr, uintptr_t unfiltered_marker_ct, const uintptr_t* marker_exclude_orig, uint32_t marker_ct, uint32_t count_x, uint32_t count_mt, const char* calc_descrip, uintptr_t** marker_exclude_ptr, uint32_t* newly_excluded_ct_ptr);

uint32_t get_max_chrom_size(const Chrom_info* chrom_info_ptr, const uintptr_t* marker_exclude, uint32_t* last_chrom_fo_idx_ptr);

void count_genders(const uintptr_t* __restrict sex_nm, const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sample_exclude, uintptr_t unfiltered_sample_ct, uint32_t* __restrict male_ct_ptr, uint32_t* __restrict female_ct_ptr, uint32_t* __restrict unk_ct_ptr);

void reverse_loadbuf(uintptr_t unfiltered_sample_ct, unsigned char* loadbuf);

// deprecated, try to just use copy_quaterarr_nonempty_subset()
void copy_quaterarr_nonempty_subset_excl(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_excl, uint32_t raw_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict output_quaterarr);

HEADER_INLINE uint32_t load_raw(uintptr_t unfiltered_sample_ct4, FILE* bedfile, uintptr_t* rawbuf) {
  // only use this if all accesses to the data involve
  // 1. some sort of mask, or
  // 2. explicit iteration from 0..(unfiltered_sample_ct-1).
  // otherwise improper trailing bits might cause a segfault, when we should
  // be ignoring them or just issuing a warning.
  return (fread(rawbuf, 1, unfiltered_sample_ct4, bedfile) < unfiltered_sample_ct4);
}

HEADER_INLINE uintptr_t get_final_mask(uint32_t sample_ct) {
  uint32_t uii = sample_ct % BITCT2;
  if (uii) {
    return (ONELU << (2 * uii)) - ONELU;
  } else {
    return ~ZEROLU;
  }
}

HEADER_INLINE uint32_t load_raw2(uintptr_t unfiltered_sample_ct4, uintptr_t unfiltered_sample_ctl2m1, uintptr_t final_mask, FILE* bedfile, uintptr_t* rawbuf) {
  if (fread(rawbuf, 1, unfiltered_sample_ct4, bedfile) < unfiltered_sample_ct4) {
    return 1;
  }
  rawbuf[unfiltered_sample_ctl2m1] &= final_mask;
  return 0;
}

uint32_t load_and_collapse(uint32_t unfiltered_sample_ct, uint32_t sample_ct, const uintptr_t* __restrict sample_exclude, uintptr_t final_mask, uint32_t do_reverse, FILE* bedfile, uintptr_t* __restrict rawbuf, uintptr_t* __restrict mainbuf);

// was "collapse_copy_quaterarr_incl", but this should be better way to think
// about it
void copy_quaterarr_nonempty_subset(const uintptr_t* __restrict raw_quaterarr, const uintptr_t* __restrict subset_mask, uint32_t raw_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict output_quaterarr);

/*
// in-place version of copy_quaterarr_subset (usually destroying original
// data).
// this doesn't seem to provide a meaningful advantage over
// copy_quaterarr_subset in practice, and the latter is more versatile without
// requiring much more memory.
void inplace_quaterarr_proper_subset(const uintptr_t* __restrict subset_mask, uint32_t orig_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict main_quaterarr);

HEADER_INLINE void inplace_quaterarr_subset(const uintptr_t* __restrict subset_mask, uint32_t orig_quaterarr_size, uint32_t subset_size, uintptr_t* __restrict main_quaterarr) {
  if (orig_quaterarr_size == subset_size) {
    return;
  }
  inplace_quaterarr_proper_subset(subset_mask, orig_quaterarr_size, subset_size, main_quaterarr);
}
*/

uint32_t load_and_collapse_incl(uint32_t unfiltered_sample_ct, uint32_t sample_ct, const uintptr_t* __restrict sample_include, uintptr_t final_mask, uint32_t do_reverse, FILE* bedfile, uintptr_t* __restrict rawbuf, uintptr_t* __restrict mainbuf);

// uint32_t load_and_collapse_incl_inplace(const uintptr_t* __restrict sample_include, uint32_t unfiltered_sample_ct, uint32_t sample_ct, uintptr_t final_mask, uint32_t do_reverse, FILE* bedfile, uintptr_t* __restrict mainbuf);

uint32_t load_and_split(uint32_t unfiltered_sample_ct, const uintptr_t* __restrict pheno_nm, const uintptr_t* __restrict pheno_c, FILE* bedfile, uintptr_t* __restrict rawbuf, uintptr_t* __restrict casebuf, uintptr_t* __restrict ctrlbuf);

void init_quaterarr_from_bitarr(const uintptr_t* __restrict bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict new_quaterarr);

void init_quaterarr_from_inverted_bitarr(const uintptr_t* __restrict inverted_bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict new_quaterarr);

void quatervec_01_init_invert(const uintptr_t* __restrict source_quatervec, uintptr_t entry_ct, uintptr_t* __restrict target_quatervec);

// target_vec := source_vec ANDNOT exclude_vec
// may write an extra word
void bitvec_andnot_copy(const uintptr_t* __restrict source_vec, const uintptr_t* __restrict exclude_vec, uintptr_t word_ct, uintptr_t* __restrict target_vec);

void apply_bitarr_mask_to_quaterarr_01(const uintptr_t* __restrict mask_bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* main_quaterarr);

void apply_bitarr_excl_to_quaterarr_01(const uintptr_t* __restrict excl_bitarr, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict main_quaterarr);

// excludes (excl_bitarr_1 & excl_bitarr_2).  (union can be excluded by calling
// apply_excl_to_quaterarr_01() twice.)
void apply_excl_intersect_to_quaterarr_01(const uintptr_t* __restrict excl_bitarr_1, const uintptr_t* __restrict excl_bitarr_2, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict main_quaterarr);

// initializes output_quatervec bits to 01 iff input_quatervec bits are 01,
// everything else zeroed out
void quatervec_copy_only_01(const uintptr_t* __restrict input_quatervec, uintptr_t unfiltered_sample_ct, uintptr_t* __restrict output_quatervec);

void quatervec_01_invert(uintptr_t unfiltered_sample_ct, uintptr_t* main_quatervec);

void vec_datamask(uintptr_t unfiltered_sample_ct, uint32_t matchval, uintptr_t* data_ptr, uintptr_t* mask_ptr, uintptr_t* result_ptr);

// void vec_rotate_plink1_to_plink2(uintptr_t* lptr, uint32_t word_ct);

void rotate_plink1_to_a2ct_and_copy(uintptr_t* loadbuf, uintptr_t* writebuf, uintptr_t word_ct);

void extract_collapsed_missing_bitfield(uintptr_t* lptr, uintptr_t unfiltered_sample_ct, uintptr_t* sample_include_quaterarr, uintptr_t sample_ct, uintptr_t* missing_bitfield);

void hh_reset(unsigned char* loadbuf, uintptr_t* sample_include_quaterarr, uintptr_t unfiltered_sample_ct);

void hh_reset_y(unsigned char* loadbuf, uintptr_t* sample_include_quaterarr, uintptr_t* sample_male_include_quaterarr, uintptr_t unfiltered_sample_ct);

HEADER_INLINE void haploid_fix(uint32_t hh_exists, uintptr_t* sample_include_quaterarr, uintptr_t* sample_male_include_quaterarr, uintptr_t sample_ct, uint32_t is_x, uint32_t is_y, unsigned char* loadbuf) {
  if (is_x) {
    if (hh_exists & XMHH_EXISTS) {
      hh_reset(loadbuf, sample_male_include_quaterarr, sample_ct);
    }
  } else if (is_y) {
    if (hh_exists & Y_FIX_NEEDED) {
      hh_reset_y(loadbuf, sample_include_quaterarr, sample_male_include_quaterarr, sample_ct);
    }
  } else if (hh_exists & NXMHH_EXISTS) {
    hh_reset(loadbuf, sample_include_quaterarr, sample_ct);
  }
}

uint32_t alloc_raw_haploid_filters(uint32_t unfiltered_sample_ct, uint32_t hh_exists, uint32_t is_include, uintptr_t* sample_bitarr, uintptr_t* sex_male, uintptr_t** sample_raw_include_quatervec_ptr, uintptr_t** sample_raw_male_quatervec_ptr);

void haploid_fix_multiple(uintptr_t* marker_exclude, uintptr_t marker_uidx_start, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uint32_t set_hh_missing, uint32_t set_mixed_mt_missing, uintptr_t* sample_raw_include2, uintptr_t* sample_raw_male_include2, uintptr_t unfiltered_sample_ct, uintptr_t byte_ct_per_marker, unsigned char* loadbuf);

void force_missing(unsigned char* loadbuf, uintptr_t* force_missing_include2, uintptr_t unfiltered_sample_ct);

HEADER_INLINE char sexchar(uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t sample_uidx) {
  if (is_set(sex_nm, sample_uidx)) {
    return '2' - is_set(sex_male, sample_uidx);
  } else {
    return '0';
  }
}

int32_t open_and_size_string_list(char* fname, FILE** infile_ptr, uintptr_t* list_len_ptr, uintptr_t* max_str_len_ptr);

int32_t load_string_list(FILE** infile_ptr, uintptr_t max_str_len, char* str_list);

int32_t open_and_skip_first_lines(FILE** infile_ptr, char* fname, char* loadbuf, uintptr_t loadbuf_size, uint32_t lines_to_skip);

int32_t load_to_first_token(FILE* infile, uintptr_t loadbuf_size, char comment_char, const char* file_descrip, char* loadbuf, char** bufptr_ptr, uintptr_t* line_idx_ptr);

int32_t open_and_load_to_first_token(FILE** infile_ptr, char* fname, uintptr_t loadbuf_size, char comment_char, const char* file_descrip, char* loadbuf, char** bufptr_ptr, uintptr_t* line_idx_ptr);

int32_t scan_max_strlen(char* fname, uint32_t colnum, uint32_t colnum2, uint32_t headerskip, char skipchar, uintptr_t* max_str_len_ptr, uintptr_t* max_str2_len_ptr);

int32_t scan_max_fam_indiv_strlen(char* fname, uint32_t colnum, uintptr_t* max_sample_id_len_ptr);

// void inplace_collapse_uint32(uint32_t* item_arr, uint32_t unfiltered_ct, uintptr_t* exclude_arr, uint32_t filtered_ct);

void inplace_collapse_uint32_incl(uint32_t* item_arr, uint32_t unfiltered_ct, uintptr_t* incl_arr, uint32_t filtered_ct);

char* alloc_and_init_collapsed_arr(char* item_arr, uintptr_t item_len, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t filtered_ct, uint32_t read_only);

char* alloc_and_init_collapsed_arr_incl(char* item_arr, uintptr_t item_len, uintptr_t unfiltered_ct, uintptr_t* include_arr, uintptr_t filtered_ct, uint32_t read_only);

void inplace_delta_collapse_arr(char* item_arr, uintptr_t item_len, uintptr_t filtered_ct_orig, uintptr_t filtered_ct_new, uintptr_t* exclude_orig, uintptr_t* exclude_new);

void inplace_delta_collapse_bitfield(uintptr_t* read_ptr, uint32_t filtered_ct_new, uintptr_t* exclude_orig, uintptr_t* exclude_new);

// deprecated, migrate to copy_bitarr_subset()
void copy_bitarr_subset_excl(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_excl, uint32_t raw_bitarr_size, uint32_t subset_size, uintptr_t* __restrict output_bitarr);

void copy_bitarr_subset(const uintptr_t* __restrict raw_bitarr, const uintptr_t* __restrict subset_mask, uint32_t raw_bitarr_size, uint32_t subset_size, uintptr_t* __restrict output_bitarr);

void uncollapse_copy_flip_include_arr(uintptr_t* collapsed_include_arr, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* output_exclude_arr);

void copy_when_nonmissing(uintptr_t* loadbuf, char* source, uintptr_t elem_size, uintptr_t unfiltered_sample_ct, uintptr_t missing_ct, char* dest);

uint32_t collapse_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len, uint32_t* id_starts);

HEADER_INLINE double rand_unif(void) {
  return (sfmt_genrand_uint32(&g_sfmt) + 0.5) * RECIP_2_32;
}

void range_list_init(Range_list* range_list_ptr);

void free_range_list(Range_list* range_list_ptr);

double normdist(double zz);

double rand_normal(double* secondval_ptr);

void init_sfmt64_from_sfmt32(sfmt_t* sfmt32, sfmt_t* sfmt64);

HEADER_INLINE void precompute_mods(uintptr_t sample_ct, uint32_t* precomputed_mods) {
  // sets precomputed_mods[n] = 2^32 mod (n-2)
  uintptr_t sample_idx;
  for (sample_idx = 2; sample_idx <= sample_ct; sample_idx++) {
    *precomputed_mods++ = (uint32_t)(0x100000000LLU % sample_idx);
  }
}

void generate_perm1_interleaved(uint32_t tot_ct, uint32_t set_ct, uintptr_t perm_idx, uintptr_t perm_ct, uintptr_t* perm_buf);

uint32_t cubic_real_roots(double coef_a, double coef_b, double coef_c, double* solutions);

void join_threads(pthread_t* threads, uint32_t ctp1);

#ifdef _WIN32
int32_t spawn_threads(pthread_t* threads, unsigned (__stdcall *start_routine)(void*), uintptr_t ct);
#else
int32_t spawn_threads(pthread_t* threads, void* (*start_routine)(void*), uintptr_t ct);
#endif

extern uintptr_t g_thread_spawn_ct;
extern uint32_t g_is_last_thread_block;

#ifdef _WIN32
extern HANDLE g_thread_start_next_event[];
extern HANDLE g_thread_cur_block_done_events[];

HEADER_INLINE void THREAD_BLOCK_FINISH(uintptr_t tidx) {
  SetEvent(g_thread_cur_block_done_events[tidx - 1]);
  WaitForSingleObject(g_thread_start_next_event[tidx - 1], INFINITE);
}

void join_threads2(pthread_t* threads, uint32_t ctp1, uint32_t is_last_block);

int32_t spawn_threads2(pthread_t* threads, unsigned (__stdcall *start_routine)(void*), uintptr_t ct, uint32_t is_last_block);
#else
void THREAD_BLOCK_FINISH(uintptr_t tidx);

void join_threads2(pthread_t* threads, uint32_t ctp1, uint32_t is_last_block);

int32_t spawn_threads2(pthread_t* threads, void* (*start_routine)(void*), uintptr_t ct, uint32_t is_last_block);
#endif

extern sfmt_t** g_sfmtp_arr;

uint32_t bigstack_init_sfmtp(uint32_t thread_ct);

#endif // __PLINK_COMMON_H__
