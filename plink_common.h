// Resources needed across all plink modules.

#ifndef __PLINK_COMMON_H__
#define __PLINK_COMMON_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <inttypes.h>

// Uncomment this to build this without CBLAS/CLAPACK.
// #define NOLAPACK

// Uncomment this to prevent all unstable features from being accessible from
// the command line.
// #define STABLE_BUILD

#define PROG_NAME_STR "plink"
#define PROG_NAME_CAPS "PLINK"

#ifdef _WIN32
  // needed for MEMORYSTATUSEX
  #ifndef _WIN64
    #define WINVER 0x0500
  #endif
#else // Unix
  #include <sys/stat.h>
#endif

#ifdef _WIN32
  #define PRId64 "I64d"
  #define PRIu64 "I64u"
  #define fseeko fseeko64
  #define ftello ftello64
  #include <windows.h>
  #include <process.h>
  #define pthread_t HANDLE
  #define THREAD_RET_TYPE unsigned __stdcall
  #define THREAD_RETURN return 0
  #define EOLN_STR "\r\n"
#else
  #include <pthread.h>
  #define THREAD_RET_TYPE void*
  #define THREAD_RETURN return NULL
  #ifdef __cplusplus
    #ifndef PRId64
      #define PRId64 "lld"
    #endif
  #endif
  #define EOLN_STR "\n"
#endif

#ifdef __APPLE__
  // fix OS X 10.9 build break; unfortunately, this doesn't work on some other
  // systems...
  typedef unsigned long long uint64_t;
  typedef long long int64_t;
#else
  #define uint64_t unsigned long long
  #define int64_t long long
#endif

#ifdef _WIN64
  #define __LP64__
  #define CTZLU __builtin_ctzll
  #define CLZLU __builtin_clzll
#else
  #define CTZLU __builtin_ctzl
  #define CLZLU __builtin_clzl
  #ifndef __LP64__
    #ifndef uintptr_t
      #define uintptr_t unsigned long
    #endif
    #ifndef intptr_t
      #define intptr_t long
    #endif
  #endif
#endif

#ifdef __cplusplus
  #include <algorithm>
#endif

#ifdef __LP64__
  #include <emmintrin.h>
  #define FIVEMASK 0x5555555555555555LLU
  typedef union {
    __m128 vf;
    __m128i vi;
    __m128d vd;
    uintptr_t u8[2];
    double d8[2];
    float f4[4];
    uint32_t u4[4];
  } __uni16;
  #define ZEROLU 0LLU
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

#else // not __LP64__

  #define FIVEMASK 0x55555555
  #define ZEROLU 0LU
  #define ONELU 1LU
  #ifndef PRIuPTR
    #define PRIuPTR "lu"
  #endif
  #ifndef PRIdPTR
    #define PRIdPTR "ld"
  #endif
  #define PRIxPTR2 "08lx"

#endif // __LP64__

#include "zlib-1.2.8/zlib.h"
#include "SFMT.h"

// 64MB of non-workspace memory guaranteed for now
#define NON_WKSPACE_MIN 67108864

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

// sometimes an infinity substitute that avoids the 32-bit Windows performance
// penalty
#define HUGE_DOUBLE 1.7976931348623157e308

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
#define RET_ALL_SAMPLES_EXCLUDED 12
#define RET_ALL_MARKERS_EXCLUDED 13
#define RET_NETWORK 14
#define LOAD_PHENO_LAST_COL 127

#define MISC_AFFECTION_01 1LLU
#define MISC_NONFOUNDERS 2LLU
#define MISC_MAF_SUCC 4LLU
#define MISC_FREQ_COUNTS 8LLU
#define MISC_FREQX 0x10LLU
#define MISC_KEEP_ALLELE_ORDER 0x20LLU
#define MISC_SET_HH_MISSING 0x40LLU
#define MISC_KEEP_AUTOCONV 0x80LLU
#define MISC_LOAD_CLUSTER_KEEP_NA 0x100LLU
#define MISC_WRITE_CLUSTER_OMIT_UNASSIGNED 0x200LLU
#define MISC_ALLOW_EXTRA_CHROMS 0x400LLU
#define MISC_MAKE_FOUNDERS_REQUIRE_2_MISSING 0x800LLU
#define MISC_MAKE_FOUNDERS_FIRST 0x1000LLU
#define MISC_LASSO_REPORT_ZEROES 0x2000LLU
#define MISC_LASSO_SELECT_COVARS 0x4000LLU
#define MISC_DOUBLE_ID 0x8000LLU
#define MISC_BIALLELIC_ONLY 0x10000LLU
#define MISC_BIALLELIC_ONLY_STRICT 0x20000LLU
#define MISC_BIALLELIC_ONLY_LIST 0x40000LLU
#define MISC_VCF_FILTER 0x80000LLU
#define MISC_GPLINK 0x100000LLU
#define MISC_SNPS_ONLY_NO_DI 0x200000LLU
#define MISC_IMPUTE_SEX 0x400000LLU
#define MISC_OXFORD_SNPID_CHR 0x800000LLU
#define MISC_EXTRACT_RANGE 0x1000000LLU
#define MISC_EXCLUDE_RANGE 0x2000000LLU
#define MISC_MERGEX 0x4000000LLU
#define MISC_SET_ME_MISSING 0x8000000LLU
#define MISC_SEXCHECK_YCOUNT 0x10000000LLU
#define MISC_SEXCHECK_YONLY 0x20000000LLU
#define MISC_FAMILY_CLUSTERS 0x40000000LLU
#define MISC_FILL_MISSING_A2 0x80000000LLU
#define MISC_HET_SMALL_SAMPLE 0x100000000LLU
#define MISC_FST_CC 0x200000000LLU
#define MISC_SPLIT_MERGE_NOFAIL 0x400000000LLU
#define MISC_REAL_REF_ALLELES 0x800000000LLU
#define MISC_RPLUGIN_DEBUG 0x1000000000LLU
#define MISC_MISSING_GZ 0x2000000000LLU
#define MISC_FREQ_GZ 0x4000000000LLU
#define MISC_HET_GZ 0x8000000000LLU

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
#define CALC_REGRESS_PCS 0x8000LLU
#define CALC_REGRESS_PCS_DISTANCE 0x10000LLU
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
#define CALC_WRITE_SET 0x400000000LLU
#define CALC_LD 0x800000000LLU
#define CALC_EPI 0x1000000000LLU
#define CALC_TESTMISS 0x2000000000LLU
#define CALC_TESTMISHAP 0x4000000000LLU
#define CALC_SEXCHECK 0x8000000000LLU
#define CALC_CLUMP 0x10000000000LLU
#define CALC_PCA 0x20000000000LLU
#define CALC_BLOCKS 0x40000000000LLU
#define CALC_SCORE 0x80000000000LLU
#define CALC_MENDEL 0x100000000000LLU
#define CALC_HET 0x200000000000LLU
#define CALC_FLIPSCAN 0x400000000000LLU
#define CALC_TDT 0x800000000000LLU
#define CALC_MAKE_PERM_PHENO 0x1000000000000LLU
#define CALC_QFAM 0x2000000000000LLU
#define CALC_FST 0x4000000000000LLU
#define CALC_SHOW_TAGS 0x8000000000000LLU
#define CALC_MAKE_BIM 0x10000000000000LLU
#define CALC_MAKE_FAM 0x20000000000000LLU
#define CALC_WRITE_VAR_RANGES 0x40000000000000LLU
#define CALC_DUPVAR 0x80000000000000LLU
#define CALC_RPLUGIN 0x100000000000000LLU
#define CALC_DFAM 0x200000000000000LLU
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
#define VCF_HALF_CALL_HAPLOID 2
#define VCF_HALF_CALL_MISSING 3

#define M23_MALE 1
#define M23_FEMALE 2
#define M23_FORCE_MISSING_SEX 4
#define M23_SEX 7

#define MARKER_CMS_OPTIONAL 1
#define MARKER_CMS_FORCED 2

#define UNSORTED_CHROM 1
#define UNSORTED_BP 2
// probably insert unsorted centimorgans later
#define UNSORTED_SPLIT_CHROM 4

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
#define RECODE_VCF 0x1000000
#define RECODE_TYPEMASK 0x1fffff0
#define RECODE_FID 0x2000000
#define RECODE_IID 0x4000000
#define RECODE_INCLUDE_ALT 0x8000000
#define RECODE_BGZ 0x10000000

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

#define REGRESS_PCS_NORMALIZE_PHENO 1
#define REGRESS_PCS_SEX_SPECIFIC 2
#define REGRESS_PCS_CLIP 4

#define HWE_MIDP 1
#define HWE_THRESH_MIDP 2
#define HWE_THRESH_ALL 4
#define HWE_GZ 8

#define MENDEL_FILTER 1
#define MENDEL_FILTER_VAR_FIRST 2
#define MENDEL_DUOS 4
#define MENDEL_MULTIGEN 8

#define DUMMY_MISSING_GENO 1
#define DUMMY_MISSING_PHENO 2
#define DUMMY_SCALAR_PHENO 4
#define DUMMY_ACGT 8
#define DUMMY_1234 16
#define DUMMY_12 32

#define SIMULATE_QT 1
#define SIMULATE_TAGS 2
#define SIMULATE_HAPS 4
#define SIMULATE_ACGT 8
#define SIMULATE_1234 16
#define SIMULATE_12 32

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

#define WKSPACE_MIN_MB 64
#define WKSPACE_DEFAULT_MB 2048

#define CACHELINE 64 // assumed number of bytes per cache line, for alignment
#define CACHELINE_INT32 (CACHELINE / 4)
#define CACHELINE_INT64 (CACHELINE / 8)
#define CACHELINE_WORD (CACHELINE / BYTECT)
#define CACHELINE_DBL (CACHELINE / 8)

#define CACHEALIGN(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - ONELU)))
#define CACHEALIGN_INT32(val) ((val + (CACHELINE_INT32 - 1)) & (~(CACHELINE_INT32 - ONELU)))
#define CACHEALIGN_WORD(val) ((val + (CACHELINE_WORD - 1)) & (~(CACHELINE_WORD - ONELU)))
#define CACHEALIGN_DBL(val) ((val + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - ONELU)))

// 32-bit instead of word-length bitwise not here, when val can be assumed to
// be 32-bit.
// (note that the sizeof operator "returns" an uintptr_t, not a uint32_t; hence
// the lack of sizeof in the CACHELINE_INT32, etc. definitions.)
#define CACHEALIGN32(val) ((val + (CACHELINE - 1)) & (~(CACHELINE - 1)))
#define CACHEALIGN32_INT32(val) ((val + (CACHELINE_INT32 - 1)) & (~(CACHELINE_INT32 - 1)))
#define CACHEALIGN32_WORD(val) ((val + (CACHELINE_WORD - 1)) & (~(CACHELINE_WORD - 1)))
#define CACHEALIGN32_DBL(val) ((val + (CACHELINE_DBL - 1)) & (~(CACHELINE_DBL - 1)))

#define MAXV(aa, bb) (((bb) > (aa))? (bb) : (aa))
#define MINV(aa, bb) (((aa) > (bb))? (bb) : (aa))

#define _FILE_OFFSET_BITS 64
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

#ifdef __LP64__
  #define BITCT 64
#else
  #define BITCT 32
#endif

#define BITCT2 (BITCT / 2)
#define BYTECT (BITCT / 8)

// generic maximum line length.  .ped/.vcf/etc. lines can of course be longer
#define MAXLINELEN 131072

// must be at least 2 * MAXLINELEN + 2 to support generic token loader.
#define TBUF_SIZE (2 * MAXLINELEN + 256)

// Maximum length of chromosome, variant, FID, IID, cluster, and set IDs (not
// including terminating null, that's what _P1 is for).  This value supports up
// to 8 IDs per line (maximum so far is 5, for e.g. --hom).
#define MAX_ID_LEN 16000

#define MAX_ID_LEN_P1 (MAX_ID_LEN + 1)
#define MAX_ID_LEN_STR "16000"

// maximum size of "dynamically" allocated line load buffer.  (this is the
// limit that applies to .vcf and similar files)
#define MAXLINEBUFLEN 0x7fffffc0

// Default --perm-batch-size value in most contexts.  It may actually be better
// to *avoid* a power of two due to the need for transpositions involving this
// stride; see e.g. http://danluu.com/3c-conflict/ .  This should be tested
// during PLINK 2.0 development.
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
  #define AAAAMASK 0xaaaaaaaaaaaaaaaaLLU
  // number of snp-major .bed lines to read at once for distance calc if
  // exponent is nonzero.
  #define MULTIPLEX_DIST_EXP 64
  // number of snp-major .bed lines to read at once for relationship calc
  #define MULTIPLEX_REL 60
#else
  // N.B. 32-bit version not as carefully tested or optimized, but I'll try to
  // make sure it works properly
  #define AAAAMASK 0xaaaaaaaaU
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
#define HASHMEM_S 4194304
#else
#define HASHMEM 2097152
#define HASHMEM_S 2097152
#endif

typedef struct {
  uint32_t min;
  uint32_t max;
  double alpha;
  double beta;
  double init_interval;
  double interval_slope;
} Aperm_info;

// fit 4 pathologically long IDs plus a bit extra
extern char tbuf[];

extern const char g_one_char_strs[];
extern const char* g_missing_geno_ptr;
extern const char* g_output_missing_geno_ptr;

static inline const char* cond_replace(const char* ss, const char* match_str, const char* replace_str) {
  if (ss != match_str) {
    return ss;
  } else {
    return replace_str;
  }
}

uint32_t aligned_malloc(uintptr_t** aligned_pp, uintptr_t size);

void aligned_free(uintptr_t* aligned_ptr);

static inline void aligned_free_cond(uintptr_t* aligned_ptr) {
  if (aligned_ptr) {
    aligned_free(aligned_ptr);
  }
}

static inline void aligned_free_null(uintptr_t** aligned_pp) {
  aligned_free(*aligned_pp);
  *aligned_pp = NULL;
}

static inline void aligned_free_cond_null(uintptr_t** aligned_pp) {
  if (*aligned_pp) {
    aligned_free(*aligned_pp);
    *aligned_pp = NULL;
  }
}

extern sfmt_t sfmt;

extern const char errstr_fopen[];
extern const char cmdline_format_str[];

extern FILE* logfile;
extern char logbuf[];
extern uint32_t g_debug_on;
extern uint32_t g_log_failed;

extern uintptr_t g_sample_ct;
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

typedef union {
  float ff;
  int32_t ii;
} __floatint32;

typedef union {
  double dd;
#ifdef __LP64__
  uintptr_t uu[1];
#else
  uintptr_t uu[2];
#endif
} __double_ulong;

uint32_t push_ll_str(Ll_str** ll_stack_ptr, const char* ss);

// warning: do NOT include allele codes (unless they're guaranteed to be SNPs)
// in log strings; they can overflow the buffer.
void logstr(const char* ss);

void logprint(const char* ss);

void logprintb();

#define LOGPRINTF(...) sprintf(logbuf, __VA_ARGS__); logprintb();

// input for wordwrap/LOGPRINTFWW should have no intermediate '\n's.  If
// suffix_len is 0, there should be a terminating \n.
void wordwrap(char* ss, uint32_t suffix_len);

#define LOGPREPRINTFWW(...) sprintf(logbuf, __VA_ARGS__); wordwrap(logbuf, 0);

#define LOGPRINTFWW(...) sprintf(logbuf, __VA_ARGS__); wordwrap(logbuf, 0); logprintb();

// 5 = length of "done." suffix, which is commonly used
#define LOGPRINTFWW5(...) sprintf(logbuf, __VA_ARGS__); wordwrap(logbuf, 5); logprintb();

#ifdef STABLE_BUILD
  #define UNSTABLE(val) sptr = strcpya(&(logbuf[9]), val); goto main_unstable_disabled
#else
  #define UNSTABLE(val)
#endif

int32_t fopen_checked(FILE** target_ptr, const char* fname, const char* mode);

static inline int32_t putc_checked(int32_t ii, FILE* outfile) {
  putc(ii, outfile);
  return ferror(outfile);
}

static inline int32_t fputs_checked(const char* ss, FILE* outfile) {
  fputs(ss, outfile);
  return ferror(outfile);
}

int32_t fwrite_checked(const void* buf, size_t len, FILE* outfile);

static inline int32_t fwrite_checkedz(const void* buf, size_t len, FILE* outfile) {
  if (len) {
    return fwrite_checked(buf, len, outfile);
  }
  return ferror(outfile);
}

static inline int32_t fread_checked(char* buf, uintptr_t len, FILE* infile, uintptr_t* bytes_read_ptr) {
  *bytes_read_ptr = fread(buf, 1, len, infile);
  return ferror(infile);
}

static inline void fclose_cond(FILE* fptr) {
  if (fptr) {
    fclose(fptr);
  }
}

static inline int32_t fclose_null(FILE** fptr_ptr) {
  int32_t ii = ferror(*fptr_ptr);
  int32_t jj = fclose(*fptr_ptr);
  *fptr_ptr = NULL;
  return ii || jj;
}

int32_t gzopen_checked(gzFile* target_ptr, const char* fname, const char* mode);

static inline int32_t gzclose_null(gzFile* gzf_ptr) {
  int32_t ii = gzclose(*gzf_ptr);
  *gzf_ptr = NULL;
  return (ii != Z_OK);
}

static inline void gzclose_cond(gzFile gz_infile) {
  if (gz_infile) {
    gzclose(gz_infile);
  }
}

static inline int32_t flexwrite_checked(const void* buf, size_t len, uint32_t output_gz, FILE* outfile, gzFile gz_outfile) {
  if (!output_gz) {
    return fwrite_checked(buf, len, outfile);
  } else {
    return (!gzwrite(gz_outfile, buf, len));
  }
}

static inline int32_t flexputc_checked(int32_t ii, uint32_t output_gz, FILE* outfile, gzFile gz_outfile) {
  if (!output_gz) {
    putc(ii, outfile);
    return ferror(outfile);
  } else {
    return (gzputc(gz_outfile, ii) == -1);
  }
}

static inline int32_t flexputs_checked(const char* ss, uint32_t output_gz, FILE* outfile, gzFile gz_outfile) {
  if (!output_gz) {
    return fputs_checked(ss, outfile);
  } else {
    return (gzputs(gz_outfile, ss) == -1);
  }
}

static inline int32_t flexclose_null(uint32_t output_gz, FILE** fptr_ptr, gzFile* gzf_ptr) {
  if (!output_gz) {
    return fclose_null(fptr_ptr);
  } else {
    return gzclose_null(gzf_ptr);
  }
}

static inline int32_t bed_suffix_conflict(uint64_t calculation_type, uint32_t recode_modifier) {
  return (calculation_type & CALC_MAKE_BED) || ((calculation_type & CALC_RECODE) && (recode_modifier & (RECODE_LGEN | RECODE_LGEN_REF | RECODE_RLIST)));
}

// manually managed, very large stack
extern unsigned char* wkspace_base;
extern uintptr_t wkspace_left;

unsigned char* wkspace_alloc(uintptr_t size);

static inline int32_t wkspace_alloc_c_checked(char** dc_ptr, uintptr_t size) {
  *dc_ptr = (char*)wkspace_alloc(size);
  return !(*dc_ptr);
}

static inline int32_t wkspace_alloc_d_checked(double** dp_ptr, uintptr_t size) {
  *dp_ptr = (double*)wkspace_alloc(size);
  return !(*dp_ptr);
}

static inline int32_t wkspace_alloc_f_checked(float** fp_ptr, uintptr_t size) {
  *fp_ptr = (float*)wkspace_alloc(size);
  return !(*fp_ptr);
}

static inline int32_t wkspace_alloc_i_checked(int32_t** ip_ptr, uintptr_t size) {
  *ip_ptr = (int32_t*)wkspace_alloc(size);
  return !(*ip_ptr);
}

static inline int32_t wkspace_alloc_uc_checked(unsigned char** ucp_ptr, uintptr_t size) {
  *ucp_ptr = wkspace_alloc(size);
  return !(*ucp_ptr);
}

static inline int32_t wkspace_alloc_ui_checked(uint32_t** uip_ptr, uintptr_t size) {
  *uip_ptr = (uint32_t*)wkspace_alloc(size);
  return !(*uip_ptr);
}

static inline int32_t wkspace_alloc_ul_checked(uintptr_t** ulp_ptr, uintptr_t size) {
  *ulp_ptr = (uintptr_t*)wkspace_alloc(size);
  return !(*ulp_ptr);
}

static inline int32_t wkspace_alloc_ll_checked(int64_t** llp_ptr, uintptr_t size) {
  *llp_ptr = (int64_t*)wkspace_alloc(size);
  return !(*llp_ptr);
}

static inline int32_t wkspace_alloc_ull_checked(uint64_t** ullp_ptr, uintptr_t size) {
  *ullp_ptr = (uint64_t*)wkspace_alloc(size);
  return !(*ullp_ptr);
}

void wkspace_reset(void* new_base);

void wkspace_shrink_top(void* rebase, uintptr_t new_size);

static inline unsigned char* top_alloc(uintptr_t* topsize_ptr, uintptr_t size) {
  uintptr_t ts = *topsize_ptr + ((size + 15) & (~(15 * ONELU)));
  if (ts > wkspace_left) {
    return NULL;
  } else {
    *topsize_ptr = ts;
    return &(wkspace_base[wkspace_left - ts]);
  }
}

static inline Ll_str* top_alloc_llstr(uintptr_t* topsize_ptr, uint32_t size) {
  return (Ll_str*)top_alloc(topsize_ptr, size + sizeof(Ll_str));
}

static inline int32_t is_letter(unsigned char ucc) {
  return (((ucc & 192) == 64) && (((ucc - 1) & 31) < 26));
}

// if we need the digit value, better to use (unsigned char)cc - '0'...
static inline int32_t is_digit(unsigned char ucc) {
  return (ucc <= '9') && (ucc >= '0');
}

static inline int32_t is_not_digit(unsigned char ucc) {
  return (ucc > '9') || (ucc < '0');
}

static inline int32_t is_not_nzdigit(unsigned char ucc) {
  return (ucc > '9') || (ucc <= '0');
}

// may as well treat all chars < 32, except tab, as eoln...
static inline int32_t is_eoln(unsigned char ucc) {
  return (ucc < 32) && (ucc != '\t');
}

// kns = "known non-space" (where tab counts as a space)
static inline int32_t is_eoln_kns(unsigned char ucc) {
  return (ucc < 32);
}

static inline int32_t is_eoln_or_comment(unsigned char ucc) {
  return (ucc < 32) || (ucc == '#');
}

static inline int32_t no_more_tokens(char* sptr) {
  return ((!sptr) || is_eoln(*sptr));
}

static inline int32_t no_more_tokens_kns(char* sptr) {
  return ((!sptr) || is_eoln_kns(*sptr));
}

static inline char* skip_initial_spaces(char* sptr) {
  while ((*sptr == ' ') || (*sptr == '\t')) {
    sptr++;
  }
  return sptr;
}

/*
static inline int32_t is_space_or_eoln(unsigned char cc) {
  // ' ', \t, \n, \0, \r
#ifdef __LP64__
  return (ucc <= 32) && (0x100002601LLU & (1LLU << ucc));
#else
  return ((ucc <= 32) && ((ucc == ' ') || (0x2601LU & (ONELU << ucc))));
#endif
}
*/
static inline int32_t is_space_or_eoln(unsigned char ucc) {
  return (ucc <= 32);
}

uint32_t match_upper(char* ss, const char* fixed_str);

uint32_t match_upper_nt(char* ss, const char* fixed_str, uint32_t ct);

uint32_t scan_posint_capped(char* ss, uint32_t* valp, uint32_t cap_div_10, uint32_t cap_mod_10);

uint32_t scan_uint_capped(char* ss, uint32_t* valp, uint32_t cap_div_10, uint32_t cap_mod_10);

uint32_t scan_int_abs_bounded(char* ss, int32_t* valp, uint32_t bound_div_10, uint32_t bound_mod_10);

// intentionally rejects -2^31 for now
static inline uint32_t scan_int32(char* ss, int32_t* valp) {
  return scan_int_abs_bounded(ss, valp, 0x7fffffff / 10, 0x7fffffff % 10);
}

// default cap = 0x7ffffffe
static inline uint32_t scan_posint_defcap(char* ss, uint32_t* valp) {
  return scan_posint_capped(ss, valp, 0x7ffffffe / 10, 0x7ffffffe % 10);
}

static inline uint32_t scan_uint_defcap(char* ss, uint32_t* valp) {
  return scan_uint_capped(ss, valp, 0x7ffffffe / 10, 0x7ffffffe % 10);
}

static inline uint32_t scan_int_abs_defcap(char* ss, int32_t* valp) {
  return scan_int_abs_bounded(ss, valp, 0x7ffffffe / 10, 0x7ffffffe % 10);
}

static inline uint32_t scan_uint_icap(char* ss, uint32_t* valp) {
  return scan_uint_capped(ss, valp, 0x7fffffff / 10, 0x7fffffff % 10);
}

uint32_t scan_posintptr(char* ss, uintptr_t* valp);

static inline char replace_if_zero(char cc, char replacement) {
  if (cc != '0') {
    return cc;
  } else {
    return replacement;
  }
}

static inline const char* replace_if_zstr(char* ss, const char* replacement) {
  if ((ss[0] != '0') || ss[1]) {
    return ss;
  }
  return replacement;
}

static inline uint32_t scan_double(char* ss, double* valp) {
  char* ss2;
  *valp = strtod(ss, &ss2);
  return (ss == ss2);
}

static inline uint32_t scan_float(char* ss, float* valp) {
  char* ss2;
  *valp = strtof(ss, &ss2);
  return (ss == ss2);
}

uint32_t scan_two_doubles(char* ss, double* val1p, double* val2p);

int32_t scan_token_ct_len(FILE* infile, char* buf, uintptr_t half_bufsize, uintptr_t* token_ct_ptr, uintptr_t* max_token_len_ptr);

int32_t read_tokens(FILE* infile, char* buf, uintptr_t half_bufsize, uintptr_t token_ct, uintptr_t max_token_len, char* token_name_buf);

static inline char* memseta(char* target, unsigned char val, uintptr_t ct) {
  memset(target, val, ct);
  return &(target[ct]);
}

static inline char* memcpya(char* target, const void* source, uintptr_t ct) {
  memcpy(target, source, ct);
  return &(target[ct]);
}

static inline char* memcpyb(char* target, const void* source, uint32_t ct) {
  // Same as memcpya, except the return value is one byte earlier.  Generally
  // used when source is a null-terminated string of known length and we want
  // to copy the null, but sometimes we need to append later.
  memcpy(target, source, ct);
  return &(target[ct - 1]);
}

static inline char* memcpyax(char* target, const void* source, uint32_t ct, char extra_char) {
  memcpy(target, source, ct);
  target[ct] = extra_char;
  return &(target[ct + 1]);
}

static inline void memcpyx(char* target, const void* source, uint32_t ct, char extra_char) {
  memcpy(target, source, ct);
  target[ct] = extra_char;
}

static inline void memcpyl3(char* target, const void* source) {
  // when it's safe to clobber the fourth character, this is faster
  *((uint32_t*)target) = *((uint32_t*)source);
}

static inline char* memcpyl3a(char* target, const void* source) {
  memcpyl3(target, source);
  return &(target[3]);
}

static inline char* strcpya(char* target, const void* source) {
  uintptr_t slen = strlen((char*)source);
  memcpy(target, source, slen);
  return &(target[slen]);
}

static inline char* strcpyax(char* target, const void* source, char extra_char) {
  uintptr_t slen = strlen((char*)source);
  memcpy(target, source, slen);
  target[slen] = extra_char;
  return &(target[slen + 1]);
}

static inline void append_binary_eoln(char** target_ptr) {
#ifdef _WIN32
  (*target_ptr)[0] = '\r';
  (*target_ptr)[1] = '\n';
  *target_ptr += 2;
#else
  **target_ptr = '\n';
  *target_ptr += 1;
#endif
}

static inline void fputs_w4(char* ss, FILE* outfile) {
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

int32_t get_next_noncomment_excl(FILE* fptr, char** lptr_ptr, uintptr_t* line_idx_ptr, uintptr_t* marker_exclude, uintptr_t* marker_uidx_ptr);

char* token_end(char* sptr);

// does not return NULL if token ends with null char
char* token_endl(char* sptr);

// token_endl without checking if sptr == NULL
// also assumes we are currently in a token -- UNSAFE OTHERWISE
static inline char* token_endnn(char* sptr) {
  while (!is_space_or_eoln(*(++sptr)));
  return sptr;
}

void get_top_two(uint32_t* uint_arr, uintptr_t uia_size, uintptr_t* top_idx_ptr, uintptr_t* second_idx_ptr);

static inline char* uint32_encode_5_hi_uchar(char* start, uint32_t uii) {
  // tried a few bit hacks here, but turns out nothing really beats this
  *start++ = (unsigned char)((uii >> 28) | 0x80);
  *start++ = (unsigned char)((uii >> 21) | 0x80);
  *start++ = (unsigned char)((uii >> 14) | 0x80);
  *start++ = (unsigned char)((uii >> 7) | 0x80);
  *start++ = (unsigned char)(uii | 0x80);
  return start;
}

static inline uint32_t uint32_decode_5_hi_uchar(char* start) {
  uint32_t uii = ((uint32_t)((unsigned char)(*start++))) << 28;
  uii |= (((uint32_t)((unsigned char)(*start++))) & 0x7f) << 21;
  uii |= (((uint32_t)((unsigned char)(*start++))) & 0x7f) << 14;
  uii |= (((uint32_t)((unsigned char)(*start++))) & 0x7f) << 7;
  uii |= ((uint32_t)((unsigned char)(*start))) & 0x7f;
  return uii;
}

int32_t intlen(int32_t num);

// safer than token_endnn(), since it handles length zero
static inline uintptr_t strlen_se(char* ss) {
  char* ss2 = ss;
  while (!is_space_or_eoln(*ss2)) {
    ss2++;
  }
  return (uintptr_t)(ss2 - ss);
}

int32_t strcmp_se(char* s_read, const char* s_const, uint32_t len);

char* next_token(char* sptr);

char* next_token_mult(char* sptr, uint32_t ct);

static inline char* next_token_multz(char* sptr, uint32_t ct) {
  if (!ct) {
    return sptr;
  } else {
    return next_token_mult(sptr, ct);
  }
}

uint32_t count_tokens(const char* bufptr);

static inline char* fw_strcpyn(uint32_t min_width, uint32_t source_len, const char* source, char* dest) {
  // right-justified strcpy with known source length
  if (source_len < min_width) {
    memcpy(memseta(dest, 32, min_width - source_len), source, source_len);
    return &(dest[min_width]);
  } else {
    return memcpya(dest, source, source_len);
  }
}

static inline char* fw_strcpy(uint32_t min_width, const char* source, char* dest) {
  return fw_strcpyn(min_width, strlen(source), source, dest);
}

uint32_t count_and_measure_multistr(const char* multistr, uintptr_t* max_slen_ptr);

char* uint32_write(char* start, uint32_t uii);

char* int32_write(char* start, int32_t ii);

void uint32_write4(char* start, uint32_t uii);

char* int64_write(char* start, int64_t llii);

char* uint32_writew4(char* start, uint32_t uii);

char* uint32_writew6(char* start, uint32_t uii);

char* uint32_writew7(char* start, uint32_t uii);

char* uint32_writew8(char* start, uint32_t uii);

char* uint32_writew10(char* start, uint32_t uii);

char* double_e_write(char* start, double dxx);

char* float_e_write(char* start, float dxx);

char* double_f_writew2(char* start, double dxx);

char* double_f_writew3(char* start, double dxx);

char* double_f_writew96(char* start, double dxx);

char* double_f_writew74(char* start, double dxx);

static inline void zeroes_to_spaces(char* start) {
  // removes trailing zeroes
  start--;
  while (*start == '0') {
    *start-- = ' ';
  }
  if (*start == '.') {
    *start = ' ';
  }
}

static inline char* clip_zeroes(char* start) {
  start--;
  while (*start == '0') {
    *(start--) = ' ';
  }
  if (*start == '.') {
    start--;
  }
  return &(start[1]);
}

char* double_f_writew96_spaced(char* start, double dxx);

char* double_f_writew96_clipped(char* start, double dxx);

char* double_g_write(char* start, double dxx);

char* float_g_write(char* start, float dxx);

static inline char* width_force(uint32_t min_width, char* startp, char* endp) {
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

char* double_g_writewx2(char* start, double dxx, uint32_t min_width);

char* double_g_writewx3(char* start, double dxx, uint32_t min_width);

char* double_g_writewx4(char* start, double dxx, uint32_t min_width);

char* double_g_writewx8(char* start, double dxx, uint32_t min_width);

static inline char* uint32_writex(char* start, uint32_t uii, char extra_char) {
  char* penult = uint32_write(start, uii);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* int32_writex(char* start, int32_t ii, char extra_char) {
  char* penult = int32_write(start, ii);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* uint32_writew4x(char* start, uint32_t uii, char extra_char) {
  char* penult = uint32_writew4(start, uii);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* uint32_writew6x(char* start, uint32_t uii, char extra_char) {
  char* penult = uint32_writew6(start, uii);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* uint32_writew7x(char* start, uint32_t uii, char extra_char) {
  char* penult = uint32_writew7(start, uii);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* uint32_writew8x(char* start, uint32_t uii, char extra_char) {
  char* penult = uint32_writew8(start, uii);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* uint32_writew10x(char* start, uint32_t uii, char extra_char) {
  char* penult = uint32_writew10(start, uii);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* double_e_writex(char* start, double dxx, char extra_char) {
  char* penult = double_e_write(start, dxx);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* float_e_writex(char* start, float dxx, char extra_char) {
  char* penult = float_e_write(start, dxx);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* double_f_writew96x(char* start, double dxx, char extra_char) {
  char* penult = double_f_writew96(start, dxx);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* double_f_writew74x(char* start, double dxx, char extra_char) {
  char* penult = double_f_writew74(start, dxx);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* double_g_writex(char* start, double dxx, char extra_char) {
  char* penult = double_g_write(start, dxx);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* float_g_writex(char* start, float dxx, char extra_char) {
  char* penult = float_g_write(start, dxx);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* double_g_writewx3x(char* start, double dxx, uint32_t min_width, char extra_char) {
  char* penult = double_g_writewx3(start, dxx, min_width);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* double_g_writewx4x(char* start, double dxx, uint32_t min_width, char extra_char) {
  char* penult = double_g_writewx4(start, dxx, min_width);
  *penult = extra_char;
  return &(penult[1]);
}

static inline char* double_g_writewx8x(char* start, double dxx, uint32_t min_width, char extra_char) {
  char* penult = double_g_writewx8(start, dxx, min_width);
  *penult = extra_char;
  return &(penult[1]);
}

static inline void read_next_terminate(char* target, char* source) {
  while (!is_space_or_eoln(*source)) {
    *target++ = *source++;
  }
  *target = '\0';
}

char* chrom_print_human(char* buf, uint32_t num);

uint32_t allele_set(char** allele_ptr, const char* newval, uint32_t slen);

uint32_t allele_reset(char** allele_ptr, const char* newval, uint32_t slen);

void magic_num(uint32_t divisor, uint64_t* multp, uint32_t* pre_shiftp, uint32_t* post_shiftp, uint32_t* incrp);

static inline uintptr_t tri_coord_no_diag(uintptr_t small_coord, uintptr_t big_coord) {
  // small_coord and big_coord are 0-based indices, small_coord < big_coord
  return ((big_coord * (big_coord - 1)) / 2) + small_coord;
}

static inline uint32_t tri_coord_no_diag_32(uint32_t small_coord, uint32_t big_coord) {
  return ((big_coord * (big_coord - 1)) / 2) + small_coord;
}

// let the compiler worry about the second argument's bit width here
#define SET_BIT(aa, bb) (aa[(bb) / BITCT] |= ONELU << ((bb) % BITCT))

#define SET_BIT_DBL(aa, bb) (aa[bb / BITCT2] |= ONELU << (2 * (bb % BITCT2)))

static inline void set_bit(uintptr_t* bit_arr, uint32_t loc) {
  bit_arr[loc / BITCT] |= (ONELU << (loc % BITCT));
}

static inline void set_bit_ul(uintptr_t* bit_arr, uintptr_t loc) {
  bit_arr[loc / BITCT] |= (ONELU << (loc % BITCT));
}

void fill_bits(uintptr_t* bit_arr, uintptr_t loc_start, uintptr_t len);

void clear_bits(uintptr_t* bit_arr, uintptr_t loc_start, uintptr_t len);

#define CLEAR_BIT(aa, bb) (aa[bb / BITCT] &= ~(ONELU << (bb % BITCT)))

#define CLEAR_BIT_DBL(aa, bb) (aa[bb / BITCT2] &= ~(ONELU << (2 * (bb % BITCT2))))

static inline void clear_bit(uintptr_t* bit_arr, uint32_t loc) {
  bit_arr[loc / BITCT] &= ~(ONELU << (loc % BITCT));
}

static inline void clear_bit_ul(uintptr_t* bit_arr, uintptr_t loc) {
  bit_arr[loc / BITCT] &= ~(ONELU << (loc % BITCT));
}

#define IS_SET(aa, bb) ((aa[bb / BITCT] >> (bb % BITCT)) & 1)

#define IS_SET_DBL(aa, bb) ((aa[bb / BITCT2] >> (2 * (bb % BITCT2))) & 1)

// use this instead of IS_SET() for signed 32-bit integers
static inline uint32_t is_set(uintptr_t* exclude_arr, uint32_t loc) {
  return (exclude_arr[loc / BITCT] >> (loc % BITCT)) & 1;
}

static inline uint32_t is_set_ul(uintptr_t* exclude_arr, uintptr_t loc) {
  return (exclude_arr[loc / BITCT] >> (loc % BITCT)) & 1;
}

#define IS_NONNULL_AND_SET(aa, bb) (aa && IS_SET(aa, bb))

uint32_t next_unset_unsafe(uintptr_t* bit_arr, uint32_t loc);

static inline void next_unset_unsafe_ck(uintptr_t* bit_arr, uint32_t* loc_ptr) {
  if (IS_SET(bit_arr, *loc_ptr)) {
    *loc_ptr = next_unset_unsafe(bit_arr, *loc_ptr);
  }
}

#ifdef __LP64__
uintptr_t next_unset_ul_unsafe(uintptr_t* bit_arr, uintptr_t loc);
#else
static inline uintptr_t next_unset_ul_unsafe(uintptr_t* bit_arr, uintptr_t loc) {
  return (uintptr_t)next_unset_unsafe(bit_arr, loc);
}
#endif

static inline void next_unset_ul_unsafe_ck(uintptr_t* bit_arr, uintptr_t* loc_ptr) {
  if (IS_SET(bit_arr, *loc_ptr)) {
    *loc_ptr = next_unset_ul_unsafe(bit_arr, *loc_ptr);
  }
}

uint32_t next_unset(uintptr_t* bit_arr, uint32_t loc, uint32_t ceil);

static inline void next_unset_ck(uintptr_t* bit_arr, uint32_t* loc_ptr, uint32_t ceil) {
  if (IS_SET(bit_arr, *loc_ptr)) {
    *loc_ptr = next_unset(bit_arr, *loc_ptr, ceil);
  }
}

#ifdef __LP64__
uintptr_t next_unset_ul(uintptr_t* bit_arr, uintptr_t loc, uintptr_t ceil);
#else
static inline uintptr_t next_unset_ul(uintptr_t* bit_arr, uintptr_t loc, uintptr_t ceil) {
  return (uintptr_t)next_unset(bit_arr, loc, ceil);
}
#endif

static inline void next_unset_ul_ck(uintptr_t* bit_arr, uintptr_t* loc_ptr, uintptr_t ceil) {
  if (IS_SET(bit_arr, *loc_ptr)) {
    *loc_ptr = next_unset_ul(bit_arr, *loc_ptr, ceil);
  }
}

uint32_t next_set_unsafe(uintptr_t* bit_arr, uint32_t loc);

static inline void next_set_unsafe_ck(uintptr_t* bit_arr, uint32_t* loc_ptr) {
  if (!IS_SET(bit_arr, *loc_ptr)) {
    *loc_ptr = next_set_unsafe(bit_arr, *loc_ptr);
  }
}

#ifdef __LP64__
uintptr_t next_set_ul_unsafe(uintptr_t* bit_arr, uintptr_t loc);
#else
static inline uintptr_t next_set_ul_unsafe(uintptr_t* bit_arr, uintptr_t loc) {
  return (uintptr_t)next_set_unsafe(bit_arr, loc);
}
#endif

static inline void next_set_ul_unsafe_ck(uintptr_t* bit_arr, uintptr_t* loc_ptr) {
  if (!IS_SET(bit_arr, *loc_ptr)) {
    *loc_ptr = next_set_ul_unsafe(bit_arr, *loc_ptr);
  }
}

uint32_t next_set(uintptr_t* bit_arr, uint32_t loc, uint32_t ceil);

static inline void next_set_ck(uintptr_t* bit_arr, uint32_t* loc_ptr, uint32_t ceil) {
  if (!IS_SET(bit_arr, *loc_ptr)) {
    *loc_ptr = next_set(bit_arr, *loc_ptr, ceil);
  }
}

#ifdef __LP64__
uintptr_t next_set_ul(uintptr_t* bit_arr, uintptr_t loc, uintptr_t ceil);
#else
static inline uintptr_t next_set_ul(uintptr_t* bit_arr, uintptr_t loc, uintptr_t ceil) {
  return (uintptr_t)next_set(bit_arr, loc, ceil);
}
#endif

static inline void next_set_ul_ck(uintptr_t* bit_arr, uintptr_t* loc_ptr, uintptr_t ceil) {
  if (!IS_SET(bit_arr, *loc_ptr)) {
    *loc_ptr = next_set_ul(bit_arr, *loc_ptr, ceil);
  }
}

int32_t last_set_bit(uintptr_t* bit_arr, uint32_t word_ct);

// note different interface from last_set_bit()
// int32_t last_clear_bit(uintptr_t* bit_arr, uint32_t ceil);

// unlike the next_[un]set family, this always returns a STRICTLY earlier
// position
uint32_t prev_unset_unsafe(uintptr_t* bit_arr, uint32_t loc);

// uint32_t prev_unset(uintptr_t* bit_arr, uint32_t loc, uint32_t floor);

static inline void prev_unset_unsafe_ck(uintptr_t* bit_arr, uint32_t* loc_ptr) {
  *loc_ptr -= 1;
  if (IS_SET(bit_arr, *loc_ptr)) {
    *loc_ptr = prev_unset_unsafe(bit_arr, *loc_ptr);
  }
}

// These functions seem to optimize better than memset(arr, 0, x) under gcc.
static inline void fill_long_zero(intptr_t* larr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *larr++ = 0;
  }
}

static inline void fill_ulong_zero(uintptr_t* ularr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *ularr++ = 0;
  }
}

#ifdef __LP64__
static inline void fill_ull_zero(uint64_t* ullarr, size_t size) {
  fill_ulong_zero((uintptr_t*)ullarr, size);
}
#else
static inline void fill_ull_zero(uint64_t* ullarr, size_t size) {
  fill_ulong_zero((uintptr_t*)ullarr, size * 2);
}
#endif

static inline void fill_long_one(intptr_t* larr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *larr++ = -1;
  }
}

static inline void fill_ulong_one(uintptr_t* ularr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *ularr++ = ~ZEROLU;
  }
}

#ifdef __LP64__
static inline void fill_ull_one(uint64_t* ullarr, size_t size) {
  fill_ulong_one((uintptr_t*)ullarr, size);
}
#else
static inline void fill_ull_one(uint64_t* ullarr, size_t size) {
  fill_ulong_one((uintptr_t*)ullarr, size * 2);
}
#endif

static inline void fill_int_zero(int32_t* iarr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *iarr++ = 0;
  }
}

static inline void fill_int_one(int32_t* iarr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *iarr++ = -1;
  }
}

static inline void fill_uint_zero(uint32_t* uiarr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *uiarr++ = 0;
  }
}

static inline void fill_uint_one(uint32_t* uiarr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *uiarr++ = ~0U;
  }
}

static inline void fill_float_zero(float* farr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *farr++ = 0.0;
  }
}

static inline void fill_double_zero(double* darr, size_t size) {
  size_t ulii;
  for (ulii = 0; ulii < size; ulii++) {
    *darr++ = 0.0;
  }
}

uint32_t murmurhash3_32(const void* key, uint32_t len);

static inline uint32_t hashval2(char* idstr, uint32_t idlen) {
  return murmurhash3_32(idstr, idlen) % HASHSIZE;
}

uintptr_t geqprime(uintptr_t floor);

static inline uint32_t get_id_htable_size(uintptr_t item_ct) {
  return (item_ct < 32761)? 65521 : geqprime(item_ct * 2 + 1);
}

int32_t populate_id_htable(uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t item_ct, const char* item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t* id_htable, uint32_t id_htable_size);

static inline int32_t alloc_and_populate_id_htable(uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t item_ct, const char* item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t** id_htable_ptr, uint32_t* id_htable_size_ptr) {
  uint32_t id_htable_size = get_id_htable_size(item_ct);
  if (wkspace_alloc_ui_checked(id_htable_ptr, id_htable_size * sizeof(int32_t))) {
    return RET_NOMEM;
  }
  *id_htable_size_ptr = id_htable_size;
  return populate_id_htable(unfiltered_ct, exclude_arr, item_ct, item_ids, max_id_len, allow_dups, *id_htable_ptr, id_htable_size);
}

uint32_t id_htable_find(const char* id_buf, uintptr_t cur_id_len, const uint32_t* id_htable, uint32_t id_htable_size, const char* item_ids, uintptr_t max_id_len);

void fill_idx_to_uidx(uintptr_t* exclude_arr, uintptr_t unfiltered_item_ct, uintptr_t item_ct, uint32_t* idx_to_uidx);

void fill_idx_to_uidx_incl(uintptr_t* include_arr, uintptr_t unfiltered_item_ct, uintptr_t item_ct, uint32_t* idx_to_uidx);

void fill_uidx_to_idx(uintptr_t* exclude_arr, uint32_t unfiltered_item_ct, uint32_t item_ct, uint32_t* uidx_to_idx);

void fill_uidx_to_idx_incl(uintptr_t* include_arr, uint32_t unfiltered_item_ct, uint32_t item_ct, uint32_t* uidx_to_idx);

void fill_midx_to_idx(uintptr_t* exclude_arr_orig, uintptr_t* exclude_arr, uint32_t item_ct, uint32_t* midx_to_idx);

void fill_vec_55(uintptr_t* vec, uint32_t ct);

void vec_collapse_init(uintptr_t* unfiltered_bitarr, uint32_t unfiltered_ct, uintptr_t* filter_bitarr, uint32_t filtered_ct, uintptr_t* output_vec);

void vec_collapse_init_exclude(uintptr_t* unfiltered_bitarr, uint32_t unfiltered_ct, uintptr_t* filter_exclude_bitarr, uint32_t filtered_ct, uintptr_t* output_vec);

uint32_t alloc_collapsed_haploid_filters(uint32_t unfiltered_sample_ct, uint32_t sample_ct, uint32_t hh_exists, uint32_t is_include, uintptr_t* sample_bitarr, uintptr_t* sex_male, uintptr_t** sample_include2_ptr, uintptr_t** sample_male_include2_ptr);

static inline void free_cond(void* memptr) {
  if (memptr) {
    free(memptr);
  }
}

static inline uint32_t realnum(double dd) {
  return (dd == dd) && (dd != INFINITY) && (dd != -INFINITY);
}

static inline double get_maf(double allele_freq) {
  if (allele_freq < 0.5) {
    return allele_freq;
  } else {
    return (1.0 - allele_freq);
  }
}

static inline int32_t filename_exists(char* fname, char* fname_end, const char* fname_append) {
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

void sample_delim_convert(uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uint32_t sample_ct, char* sample_ids, uintptr_t max_sample_id_len, char oldc, char newc);

void get_set_wrange_align(uintptr_t* bitfield, uintptr_t word_ct, uintptr_t* firstw_ptr, uintptr_t* wlen_ptr);

// Maximum accepted chromosome index is this minus 1.  Currently cannot exceed
// 2^14 due to SMALL_INTERVAL_BITS setting in plink_cnv.c...
#define MAX_POSSIBLE_CHROM 5120
// ...unless this is uncommented (it removes the entire CNV module).
// #define HIGH_MAX_CHROM

// assumes MAX_POSSIBLE_CHROM is a multiple of 64, otherwise add round-up
#define CHROM_MASK_WORDS (MAX_POSSIBLE_CHROM / BITCT)

#define MAX_CHROM_TEXTNUM 59
// usual PLINK 1.07 chromosome field length is 4.  So it's safe to increase
// MAX_CHROM_TEXTNUM to 9995, but 9996+ creates problems...
// (note that n+1, n+2, n+3, and n+4 are reserved for X/Y/XY/MT)
#define MAX_CHROM_TEXTNUM_LEN 2
#define CHROM_X MAX_POSSIBLE_CHROM
#define CHROM_Y (MAX_POSSIBLE_CHROM + 1)
#define CHROM_XY (MAX_POSSIBLE_CHROM + 2)
#define CHROM_MT (MAX_POSSIBLE_CHROM + 3)

#ifdef __LP64__
  // dog requires 42 bits, and other species require less
  #define CHROM_MASK_INITIAL_WORDS 1
#else
  #define CHROM_MASK_INITIAL_WORDS 2
#endif

typedef struct {
  // no point to dynamic allocation when MAX_POSSIBLE_CHROM is small and
  // there's only one copy of this

  // order of chromosomes in input files
  // currently tolerates out-of-order chromosomes, as long as all markers for
  // any given chromosome are together
  uint32_t chrom_file_order[MAX_POSSIBLE_CHROM];
  uint32_t chrom_ct; // length of chrom_file_order
  uint32_t chrom_file_order_marker_idx[MAX_POSSIBLE_CHROM + 1];

  // markers chrom_start[k] to (chrom_end[k] - 1) are part of chromosome k
  uint32_t chrom_start[MAX_POSSIBLE_CHROM];
  uint32_t chrom_end[MAX_POSSIBLE_CHROM];

  uintptr_t chrom_mask[CHROM_MASK_WORDS];

  uint32_t species;

  int32_t x_code;
  int32_t y_code;
  int32_t xy_code;
  int32_t mt_code;
  uint32_t max_code;

  uint32_t autosome_ct;
  // this is a misnomer--it includes X and excludes MT.  Underlying concept is
  // "are some calls guaranteed to be homozygous (assuming >= 1 male)", which
  // is no longer true for MT since heteroplasmy is a thing.
  uintptr_t haploid_mask[CHROM_MASK_WORDS];

  // --allow-extra-chroms support
  uint32_t zero_extra_chroms;
  uint32_t name_ct;
  Ll_str* incl_excl_name_stack;
  uint32_t is_include_stack;
  uint32_t output_encoding;
  char* nonstd_names[MAX_POSSIBLE_CHROM];
  uint32_t nonstd_name_order[MAX_POSSIBLE_CHROM];
} Chrom_info;

#define SPECIES_HUMAN 0
#define SPECIES_COW 1
#define SPECIES_DOG 2
#define SPECIES_HORSE 3
#define SPECIES_MOUSE 4
#define SPECIES_RICE 5
#define SPECIES_SHEEP 6
#define SPECIES_UNKNOWN 7

extern const char* g_species_singular;
extern const char* g_species_plural;

static inline const char* species_str(uintptr_t ct) {
  return (ct == ONELU)? g_species_singular : g_species_plural;
}

#define CHR_OUTPUT_PREFIX 1
#define CHR_OUTPUT_M 2
#define CHR_OUTPUT_MT 4
#define CHR_OUTPUT_0M 8

static inline uint32_t all_words_zero(uintptr_t* word_arr, uintptr_t word_ct) {
  while (word_ct--) {
    if (*word_arr++) {
      return 0;
    }
  }
  return 1;
}

char* chrom_name_write(char* buf, Chrom_info* chrom_info_ptr, uint32_t chrom_idx);

char* chrom_name_buf5w4write(char* buf5, Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uint32_t* chrom_name_len_ptr);

uint32_t get_max_chrom_len(Chrom_info* chrom_info_ptr);

void forget_extra_chrom_names(Chrom_info* chrom_info_ptr);

uint32_t haploid_chrom_present(Chrom_info* chrom_info_ptr);

int32_t get_chrom_code_raw(char* sptr);

int32_t get_chrom_code(Chrom_info* chrom_info_ptr, char* sptr);

int32_t get_chrom_code2(Chrom_info* chrom_info_ptr, char* sptr, uint32_t slen);

uint32_t get_marker_chrom_fo_idx(Chrom_info* chrom_info_ptr, uintptr_t marker_uidx);

static inline uint32_t get_marker_chrom(Chrom_info* chrom_info_ptr, uintptr_t marker_uidx) {
  return chrom_info_ptr->chrom_file_order[get_marker_chrom_fo_idx(chrom_info_ptr, marker_uidx)];
}

static inline int32_t chrom_exists(Chrom_info* chrom_info_ptr, uint32_t chrom_idx) {
  return is_set(chrom_info_ptr->chrom_mask, chrom_idx);
}

int32_t resolve_or_add_chrom_name(Chrom_info* chrom_info_ptr, char* bufptr, int32_t* chrom_idx_ptr, uintptr_t line_idx, const char* file_descrip);

// no need for this; code is simpler if we just create a copy of marker_exclude
// with all non-autosomal loci removed
/*
static inline uintptr_t next_autosomal_unsafe(uintptr_t* marker_exclude, uintptr_t marker_uidx, Chrom_info* chrom_info_ptr, uint32_t* chrom_end_ptr, uint32_t* chrom_fo_idx_ptr) {
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

void refresh_chrom_info(Chrom_info* chrom_info_ptr, uintptr_t marker_uidx, uint32_t* chrom_end_ptr, uint32_t* chrom_fo_idx_ptr, uint32_t* is_x_ptr, uint32_t* is_y_ptr, uint32_t* is_mt_ptr, uint32_t* is_haploid_ptr);

int32_t single_chrom_start(Chrom_info* chrom_info_ptr, uint32_t unfiltered_marker_ct, uintptr_t* marker_exclude);

double destructive_get_dmedian(double* unsorted_arr, uintptr_t len);

int32_t strcmp_casted(const void* s1, const void* s2);

int32_t strcmp_natural(const void* s1, const void* s2);

int32_t strcmp_deref(const void* s1, const void* s2);

int32_t strcmp_natural_deref(const void* s1, const void* s2);

int32_t get_uidx_from_unsorted(char* idstr, uintptr_t* exclude_arr, uint32_t id_ct, char* unsorted_ids, uintptr_t max_id_len);

char* scan_for_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len);

int32_t is_missing_pheno_cc(char* bufptr, double missing_phenod, uint32_t affection_01);

int32_t eval_affection(char* bufptr, double missing_phenod);

uint32_t triangle_divide(int64_t cur_prod, int32_t modif);

void triangle_fill(uint32_t* target_arr, uint32_t ct, uint32_t pieces, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t start, uint32_t align);

int32_t relationship_req(uint64_t calculation_type);

int32_t distance_req(uint64_t calculation_type, char* read_dists_fname);

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

int32_t qsort_ext(char* main_arr, intptr_t arr_length, intptr_t item_length, int(* comparator_deref)(const void*, const void*), char* secondary_arr, intptr_t secondary_item_len);

int32_t sort_item_ids_noalloc(char* sorted_ids, uint32_t* id_map, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t item_ct, char* item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t collapse_idxs, int(* comparator_deref)(const void*, const void*));

int32_t sort_item_ids(char** sorted_ids_ptr, uint32_t** id_map_ptr, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t exclude_ct, char* item_ids, uintptr_t max_id_len, uint32_t allow_dups, uint32_t collapse_idxs, int(* comparator_deref)(const void*, const void*));

uint32_t uint32arr_greater_than(uint32_t* sorted_uint32_arr, uint32_t arr_length, uint32_t uii);

uint32_t int32arr_greater_than(int32_t* sorted_int32_arr, uint32_t arr_length, int32_t ii);

uintptr_t uint64arr_greater_than(uint64_t* sorted_uint64_arr, uintptr_t arr_length, uint64_t ullii);

uintptr_t doublearr_greater_than(double* sorted_dbl_arr, uintptr_t arr_length, double dxx);

uintptr_t nonincr_doublearr_leq_stride(double* nonincr_dbl_arr, uintptr_t arr_length, uintptr_t stride, double dxx);

int32_t bsearch_str(const char* id_buf, uintptr_t cur_id_len, char* lptr, uintptr_t max_id_len, uintptr_t end_idx);

static inline int32_t bsearch_str_nl(const char* id_buf, char* lptr, uintptr_t max_id_len, intptr_t end_idx) {
  return bsearch_str(id_buf, strlen(id_buf), lptr, max_id_len, end_idx);
}

int32_t bsearch_str_natural(char* id_buf, char* lptr, uintptr_t max_id_len, uintptr_t end_idx);

uintptr_t bsearch_str_lb(const char* id_buf, uintptr_t cur_id_len, char* lptr, uintptr_t max_id_len, uintptr_t end_idx);

uint32_t bsearch_read_fam_indiv(char* id_buf, char* lptr, uintptr_t max_id_len, uintptr_t filter_line_ct, char* read_ptr, char** read_pp_new, int32_t* retval_ptr);

void bsearch_fam(char* id_buf, char* lptr, uintptr_t max_id_len, uint32_t filter_line_ct, char* fam_id, uint32_t* first_idx_ptr, uint32_t* last_idx_ptr);

void bitfield_invert(uintptr_t* bit_arr, uintptr_t bit_ct);

void bitfield_exclude_to_include(uintptr_t* exclude_arr, uintptr_t* include_arr, uintptr_t bit_ct);

void bitfield_and(uintptr_t* vv, uintptr_t* include_vec, uintptr_t word_ct);

void bitfield_andnot(uintptr_t* vv, uintptr_t* exclude_vec, uintptr_t word_ct);

void bitfield_andnot_reversed_args(uintptr_t* vv, uintptr_t* include_vec, uintptr_t word_ct);

void bitfield_or(uintptr_t* vv, uintptr_t* or_vec, uintptr_t word_ct);

void bitfield_ornot(uintptr_t* vv, uintptr_t* inverted_or_vec, uintptr_t word_ct);

void bitfield_xor(uintptr_t* bit_arr, uintptr_t* xor_arr, uintptr_t word_ct);

static inline uint32_t popcount2_long(uintptr_t val) {
#ifdef __LP64__
  val = (val & 0x3333333333333333LLU) + ((val >> 2) & 0x3333333333333333LLU);
  return (((val + (val >> 4)) & 0x0f0f0f0f0f0f0f0fLLU) * 0x0101010101010101LLU) >> 56;
#else
  val = (val & 0x33333333) + ((val >> 2) & 0x33333333);
  return (((val + (val >> 4)) & 0x0f0f0f0f) * 0x01010101) >> 24;
#endif
}

static inline uint32_t popcount_long(uintptr_t val) {
  // the simple version, good enough for all non-time-critical stuff
  return popcount2_long(val - ((val >> 1) & FIVEMASK));
}

uint32_t is_monomorphic_a2(uintptr_t* lptr, uint32_t sample_ct);

uint32_t is_monomorphic(uintptr_t* lptr, uint32_t sample_ct);

// same as is_monomorphic, except it also flags the all-heterozygote case
uint32_t less_than_two_genotypes(uintptr_t* lptr, uint32_t sample_ct);

// uint32_t has_three_genotypes(uintptr_t* lptr, uint32_t sample_ct);

uintptr_t popcount_longs(uintptr_t* lptr, uintptr_t word_ct);

#ifdef __LP64__
static inline uintptr_t popcount_longs_nzbase(uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
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
static inline uintptr_t popcount_longs_nzbase(uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
  return popcount_longs(&(lptr[start_idx]), end_idx - start_idx);
}
#endif

uintptr_t popcount2_longs(uintptr_t* lptr, uintptr_t word_ct);

#define popcount01_longs popcount2_longs

uintptr_t popcount_bit_idx(uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx);

uint32_t chrom_window_max(uint32_t* marker_pos, uintptr_t* marker_exclude, Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max);

uint32_t window_back(uint32_t* marker_pos, uintptr_t* marker_exclude, uint32_t marker_uidx_min, uint32_t marker_uidx, uint32_t count_max, uint32_t bp_max, uint32_t* window_trail_ct_ptr);

uint32_t window_forward(uint32_t* marker_pos, uintptr_t* marker_exclude, uint32_t marker_uidx_start, uint32_t marker_uidx_last, uint32_t count_max, uint32_t bp_max, uint32_t* window_lead_ct_ptr);

uintptr_t jump_forward_unset_unsafe(uintptr_t* bit_arr, uintptr_t cur_pos, uintptr_t forward_ct);

static inline uintptr_t popcount_chars(uintptr_t* lptr, uintptr_t start_idx, uintptr_t end_idx) {
  return popcount_bit_idx(lptr, start_idx * 8, end_idx * 8);
}

uintptr_t popcount_longs_exclude(uintptr_t* lptr, uintptr_t* exclude_arr, uintptr_t end_idx);

uintptr_t popcount_longs_intersect(uintptr_t* lptr1, uintptr_t* lptr2, uintptr_t word_ct);

void vertical_bitct_subtract(uintptr_t* bit_arr, uint32_t item_ct, uint32_t* sum_arr);

#ifdef __LP64__
void count_2freq_dbl_60v(__m128i* vptr, __m128i* vend, __m128i* mask1vp, __m128i* mask2vp, uint32_t* ct1abp, uint32_t* ct1cp, uint32_t* ct2abp, uint32_t* ct2cp);

void count_3freq_120v(__m128i* vptr, __m128i* vend, __m128i* maskvp, uint32_t* ctap, uint32_t* ctbp, uint32_t* ctcp);
#else
void count_2freq_dbl_6(uintptr_t* lptr, uintptr_t* mask1p, uintptr_t* mask2p, uint32_t* ct1abp, uint32_t* ct1cp, uint32_t* ct2abp, uint32_t* ct2cp);

void count_3freq_12(uintptr_t* lptr, uintptr_t* maskp, uint32_t* ctap, uint32_t* ctbp, uint32_t* ctcp);
#endif

void vec_set_freq(uintptr_t sample_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uint32_t* set_ctp, uint32_t* missing_ctp);

void vec_set_freq_x(uintptr_t sample_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* male_vec, uint32_t* set_ctp, uint32_t* missing_ctp);

void vec_set_freq_y(uintptr_t sample_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uintptr_t* nonmale_vec, uint32_t* set_ctp, uint32_t* missing_ctp);

void vec_3freq(uintptr_t sample_ctl2, uintptr_t* lptr, uintptr_t* include_vec, uint32_t* missing_ctp, uint32_t* het_ctp, uint32_t* homa2_ctp);

uintptr_t count_01(uintptr_t* lptr, uintptr_t word_ct);

static inline void zero_trailing_bits(uintptr_t* bitfield, uintptr_t unfiltered_ct) {
  uintptr_t trail_ct = unfiltered_ct & (BITCT - 1);
  if (trail_ct) {
    bitfield[unfiltered_ct / BITCT] &= (ONELU << trail_ct) - ONELU;
  }
}

void fill_all_bits(uintptr_t* bit_arr, uintptr_t ct);

uint32_t numeric_range_list_to_bitfield(Range_list* range_list_ptr, uint32_t item_ct, uintptr_t* bitfield, uint32_t offset, uint32_t ignore_overflow);

int32_t string_range_list_to_bitfield(char* header_line, uint32_t item_ct, uint32_t fixed_len, Range_list* range_list_ptr, char* sorted_ids, uint32_t* id_map, int32_t* seen_idxs, const char* range_list_flag, const char* file_descrip, uintptr_t* bitfield);

int32_t string_range_list_to_bitfield_alloc(char* header_line, uint32_t item_ct, uint32_t fixed_len, Range_list* range_list_ptr, uintptr_t** bitfield_ptr, const char* range_list_flag, const char* file_descrip);

int32_t string_range_list_to_bitfield2(char* sorted_ids, uint32_t* id_map, uintptr_t item_ct, uintptr_t max_id_len, Range_list* range_list_ptr, const char* range_list_flag, uintptr_t* bitfield_excl);

static inline uint32_t count_chrom_markers(Chrom_info* chrom_info_ptr, uint32_t chrom_idx, uintptr_t* marker_exclude) {
  uint32_t min_idx;
  uint32_t max_idx;
  if (!is_set(chrom_info_ptr->chrom_mask, chrom_idx)) {
    return 0;
  }
  min_idx = chrom_info_ptr->chrom_start[chrom_idx];
  max_idx = chrom_info_ptr->chrom_end[chrom_idx];
  return (max_idx - min_idx) - ((uint32_t)popcount_bit_idx(marker_exclude, min_idx, max_idx));
}

uint32_t count_non_autosomal_markers(Chrom_info* chrom_info_ptr, uintptr_t* marker_exclude, uint32_t count_x, uint32_t count_mt);

int32_t conditional_allocate_non_autosomal_markers(Chrom_info* chrom_info_ptr, uintptr_t unfiltered_marker_ct, uintptr_t* marker_exclude_orig, uint32_t marker_ct, uint32_t count_x, uint32_t count_mt, const char* calc_descrip, uintptr_t** marker_exclude_ptr, uint32_t* newly_excluded_ct_ptr);

uint32_t get_max_chrom_size(Chrom_info* chrom_info_ptr, uintptr_t* marker_exclude, uint32_t* last_chrom_fo_idx_ptr);

void count_genders(uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t unfiltered_sample_ct, uintptr_t* sample_exclude, uint32_t* male_ct_ptr, uint32_t* female_ct_ptr, uint32_t* unk_ct_ptr);

void reverse_loadbuf(unsigned char* loadbuf, uintptr_t unfiltered_sample_ct);

void collapse_copy_2bitarr(uintptr_t* rawbuf, uintptr_t* mainbuf, uint32_t unfiltered_sample_ct, uint32_t sample_ct, uintptr_t* sample_exclude);

static inline uint32_t load_raw(FILE* bedfile, uintptr_t* rawbuf, uintptr_t unfiltered_sample_ct4) {
  // only use this if all accesses to the data involve
  // 1. some sort of mask, or
  // 2. explicit iteration from 0..(unfiltered_sample_ct-1).
  // otherwise improper trailing bits might cause a segfault, when we should
  // be ignoring them or just issuing a warning.
  return (fread(rawbuf, 1, unfiltered_sample_ct4, bedfile) < unfiltered_sample_ct4);
}

static inline uintptr_t get_final_mask(uint32_t sample_ct) {
  uint32_t uii = sample_ct % BITCT2;
  if (uii) {
    return (ONELU << (2 * uii)) - ONELU;
  } else {
    return ~ZEROLU;
  }
}

static inline uint32_t load_raw2(FILE* bedfile, uintptr_t* rawbuf, uintptr_t unfiltered_sample_ct4, uintptr_t unfiltered_sample_ctl2m1, uintptr_t final_mask) {
  if (fread(rawbuf, 1, unfiltered_sample_ct4, bedfile) < unfiltered_sample_ct4) {
    return 1;
  }
  rawbuf[unfiltered_sample_ctl2m1] &= final_mask;
  return 0;
}

uint32_t load_and_collapse(FILE* bedfile, uintptr_t* rawbuf, uint32_t unfiltered_sample_ct, uintptr_t* mainbuf, uint32_t sample_ct, uintptr_t* sample_exclude, uintptr_t final_mask, uint32_t do_reverse);

void collapse_copy_2bitarr_incl(uintptr_t* rawbuf, uintptr_t* mainbuf, uint32_t unfiltered_sample_ct, uint32_t sample_ct, uintptr_t* sample_include);

uint32_t load_and_collapse_incl(FILE* bedfile, uintptr_t* rawbuf, uint32_t unfiltered_sample_ct, uintptr_t* mainbuf, uint32_t sample_ct, uintptr_t* sample_include, uintptr_t final_mask, uint32_t do_reverse);

uint32_t load_and_split(FILE* bedfile, uintptr_t* rawbuf, uint32_t unfiltered_sample_ct, uintptr_t* casebuf, uintptr_t* ctrlbuf, uintptr_t* pheno_nm, uintptr_t* pheno_c);

void vec_include_init(uintptr_t unfiltered_sample_ct, uintptr_t* new_include2, uintptr_t* old_include);

void exclude_to_vec_include(uintptr_t unfiltered_sample_ct, uintptr_t* include_vec, uintptr_t* exclude_arr);

void vec_init_invert(uintptr_t vec_entry_ct, uintptr_t* target_arr, uintptr_t* source_arr);

void vec_init_andnot(uintptr_t vec_entry_ct, uintptr_t* target_arr, uintptr_t* source_arr, uintptr_t* exclude_arr);

void vec_include_mask_in(uintptr_t unfiltered_sample_ct, uintptr_t* include_arr, uintptr_t* mask_arr);

void vec_include_mask_out(uintptr_t unfiltered_sample_ct, uintptr_t* include_arr, uintptr_t* mask_arr);

void vec_include_mask_out_intersect(uintptr_t unfiltered_sample_ct, uintptr_t* include_arr, uintptr_t* mask_arr, uintptr_t* mask2_arr);

void vec_init_01(uintptr_t unfiltered_sample_ct, uintptr_t* data_ptr, uintptr_t* result_ptr);

void vec_invert(uintptr_t unfiltered_sample_ct, uintptr_t* vec2);

void vec_datamask(uintptr_t unfiltered_sample_ct, uint32_t matchval, uintptr_t* data_ptr, uintptr_t* mask_ptr, uintptr_t* result_ptr);

// void vec_rotate_plink1_to_plink2(uintptr_t* lptr, uint32_t word_ct);

void rotate_plink1_to_plink2_and_copy(uintptr_t* loadbuf, uintptr_t* writebuf, uintptr_t word_ct);

void extract_collapsed_missing_bitfield(uintptr_t* lptr, uintptr_t unfiltered_sample_ct, uintptr_t* sample_include2, uintptr_t sample_ct, uintptr_t* missing_bitfield);

void hh_reset(unsigned char* loadbuf, uintptr_t* sample_include2, uintptr_t unfiltered_sample_ct);

void hh_reset_y(unsigned char* loadbuf, uintptr_t* sample_include2, uintptr_t* sample_male_include2, uintptr_t unfiltered_sample_ct);

static inline void haploid_fix(uint32_t hh_exists, uintptr_t* sample_include2, uintptr_t* sample_male_include2, uintptr_t sample_ct, uint32_t is_x, uint32_t is_y, unsigned char* loadbuf) {
  if (is_x) {
    if (hh_exists & XMHH_EXISTS) {
      hh_reset(loadbuf, sample_male_include2, sample_ct);
    }
  } else if (is_y) {
    if (hh_exists & Y_FIX_NEEDED) {
      hh_reset_y(loadbuf, sample_include2, sample_male_include2, sample_ct);
    }
  } else if (hh_exists & NXMHH_EXISTS) {
    hh_reset(loadbuf, sample_include2, sample_ct);
  }
}

uint32_t alloc_raw_haploid_filters(uint32_t unfiltered_sample_ct, uint32_t hh_exists, uint32_t is_include, uintptr_t* sample_bitarr, uintptr_t* sex_male, uintptr_t** sample_raw_include2_ptr, uintptr_t** sample_raw_male_include2_ptr);

void haploid_fix_multiple(uintptr_t* marker_exclude, uintptr_t marker_uidx_start, uintptr_t marker_ct, Chrom_info* chrom_info_ptr, uint32_t hh_exists, uintptr_t* sample_raw_include2, uintptr_t* sample_raw_male_include2, uintptr_t unfiltered_sample_ct, uintptr_t byte_ct_per_marker, unsigned char* loadbuf);

void force_missing(unsigned char* loadbuf, uintptr_t* force_missing_include2, uintptr_t unfiltered_sample_ct);

static inline char sexchar(uintptr_t* sex_nm, uintptr_t* sex_male, uintptr_t sample_uidx) {
  return is_set(sex_nm, sample_uidx)? (is_set(sex_male, sample_uidx)? '1' : '2') : '0';
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

void collapse_copy_bitarr(uint32_t orig_ct, uintptr_t* bit_arr, uintptr_t* exclude_arr, uint32_t filtered_ct, uintptr_t* output_arr);

void collapse_copy_bitarr_incl(uint32_t orig_ct, uintptr_t* bit_arr, uintptr_t* include_arr, uint32_t filtered_ct, uintptr_t* output_arr);

void uncollapse_copy_flip_include_arr(uintptr_t* collapsed_include_arr, uintptr_t unfiltered_ct, uintptr_t* exclude_arr, uintptr_t* output_exclude_arr);

void copy_when_nonmissing(uintptr_t* loadbuf, char* source, uintptr_t elem_size, uintptr_t unfiltered_sample_ct, uintptr_t missing_ct, char* dest);

uint32_t collapse_duplicate_ids(char* sorted_ids, uintptr_t id_ct, uintptr_t max_id_len, uint32_t* id_starts);

static inline double rand_unif(void) {
  return (sfmt_genrand_uint32(&sfmt) + 0.5) * RECIP_2_32;
}

void range_list_init(Range_list* range_list_ptr);

void free_range_list(Range_list* range_list_ptr);

double normdist(double zz);

double rand_normal(double* secondval_ptr);

void init_sfmt64_from_sfmt32(sfmt_t* sfmt32, sfmt_t* sfmt64);

static inline void precompute_mods(uintptr_t sample_ct, uint32_t* precomputed_mods) {
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

static inline void THREAD_BLOCK_FINISH(uintptr_t tidx) {
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

uint32_t wkspace_init_sfmtp(uint32_t thread_ct);

#endif // __PLINK_COMMON_H__
