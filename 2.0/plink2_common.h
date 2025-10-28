#ifndef __PLINK2_COMMON_H__
#define __PLINK2_COMMON_H__

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


// Resources needed across a variety of plink2 modules involving
// plink2-specific constructs (e.g. chr_info_t, variable-length sample IDs,
// 2-bit genotypes).  More generic library code has been moved to plink2_base
// and plink2_cmdline.

#include <assert.h>
#include <errno.h>
#include <string.h>

#include "include/pgenlib_misc.h"
#include "include/pgenlib_read.h"
#include "include/plink2_base.h"
#include "include/plink2_bits.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "plink2_cmdline.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Uncomment this if the default contig-count limit of ~65k is too low, and you
// want to rebuild this program to support ~980k instead.
// #define HIGH_CONTIG_BUILD

// Uncomment this if the default autosome-pair limit of 95 is too low, and you
// want to rebuild this program to support 995 instead.
// #define HIGH_AUTOSOME_NUM_BUILD

#define PROG_NAME_STR "plink2"

extern const double kSmallDoubles[4];

// Exclude 0x7fffffff, since the half-open interval containing it ends in
// 0x80000000U or higher.
CONSTI32(kMaxBp, 0x7ffffffe);

// leave the door semi-open to 32-bit dosages (or 31?  24?)
// note that u31tod()/u31tof() can't be used on 32-bit dosages
typedef uint16_t Dosage;
typedef int16_t SDosage;
typedef uint32_t DosageProd;
#define kDosageMax (1U << (8 * sizeof(Dosage) - 1))
CONSTI32(kDosageMid, kDosageMax / 2);
CONSTI32(kDosage4th, kDosageMax / 4);
CONSTI32(kDosageMissing, kDosageMax * 2 - 1);
static const double kRecipDosageMax = 0.000030517578125;
static const double kRecipDosageMid = 0.00006103515625;
static const double kRecipDosageMidSq = 0.0000000037252902984619140625;
static const float kRecipDosageMidf = 0.00006103515625;
CONSTI32(kDosagePerVec, kBytesPerVec / sizeof(Dosage));
CONSTI32(kDosagePerCacheline, kCacheline / sizeof(Dosage));

HEADER_CINLINE uintptr_t DosageCtToVecCt(uintptr_t val) {
  return DivUp(val, kDosagePerVec);
}

HEADER_CINLINE uintptr_t DosageCtToAlignedWordCt(uintptr_t val) {
  return kWordsPerVec * DosageCtToVecCt(val);
}

// dosage_int = 0..2 value in 16384ths
// returns distance from 0.5 or 1.5 in 16384ths, whichever is closer
HEADER_INLINE uint32_t BiallelicDosageHalfdist(uint32_t dosage_int) {
  const uint32_t dosage_int_rem = dosage_int & (kDosageMid - 1);
  return abs_i32(S_CAST(int32_t, dosage_int_rem) - kDosage4th);
}

HEADER_INLINE uint32_t HaploidDosageHalfdist(uint32_t dosage_int) {
  return abs_i32(S_CAST(int32_t, dosage_int) - kDosage4th);
}

HEADER_INLINE uint32_t DosageHetdist(uint32_t dosage_int) {
  return abs_i32(S_CAST(int32_t, dosage_int) - kDosageMid);
}

HEADER_INLINE uint32_t DosageHomdist(uint32_t dosage_int) {
  // todo: compare vs. branchless
  return (dosage_int > kDosageMid)? (kDosageMax - dosage_int) : dosage_int;
}

// on 0..32768 scale
HEADER_INLINE uint32_t DphaseHalfdist(uint32_t dphase_side) {
  return abs_i32(S_CAST(int32_t, dphase_side) - kDosageMid);
}

// this is a bit arbitrary
CONSTI32(kMaxPhenoCt, 524287);
#define MAX_PHENO_CT_STR "524287"

// maximum number of usable cluster computers, this is arbitrary though it
// shouldn't be larger than 2^32 - 1
// (actually, there's an overflow danger: <work units> * parallel_idx may not
// fit in a uint64 if parallel_tot is too high.)
CONSTI32(kParallelMax, 32768);

// unnecessary to use e.g. (1LLU << 0), the FLAGSET64 macros should force the
// integer type to 64-bit.
FLAGSET64_DEF_START()
  kfMisc0,
  kfMiscAffection01 = (1 << 0),
  kfMiscProhibitExtraChr = (1 << 1),
  kfMiscRealRefAlleles = (1 << 2),
  kfMiscMajRef = (1 << 3),
  kfMiscMajRefForce = (1 << 4),
  kfMiscNonfounders = (1 << 5),
  kfMiscMergePar = (1 << 6),
  kfMiscPhenoColNums = (1 << 7),
  kfMiscCovarColNums = (1 << 8),
  kfMiscHweMidp = (1 << 9),
  kfMiscHweKeepFewhet = (1 << 10),
  kfMiscWriteSnplistZs = (1 << 11),
  kfMiscWriteSnplistAllowDups = (1 << 12),
  kfMiscGenoDosage = (1 << 13),
  kfMiscGenoHhMissing = (1 << 14),
  kfMiscMindDosage = (1 << 15),
  kfMiscMindHhMissing = (1 << 16),
  kfMiscGenotypingRateDosage = (1 << 17),
  kfMiscSetMissingVarIds = (1 << 18),
  kfMiscChrOverrideCmdline = (1 << 19),
  kfMiscChrOverrideFile = (1 << 20),
  kfMiscNewVarIdOverflowMissing = (1 << 21),
  kfMiscNewVarIdOverflowTruncate = (1 << 22),
  kfMiscRequirePheno = (1 << 23),
  kfMiscRequireCovar = (1 << 24),
  kfMiscCatPhenoFamily = (1 << 25),
  kfMiscRefAlleleForce = (1 << 26),
  kfMiscAltAlleleForce = (1 << 27),
  kfMiscMergeX = (1 << 28),
  kfMiscNoIdHeader = (1 << 29),
  kfMiscNoIdHeaderIidOnly = (1 << 30),
  kfMiscBiallelicOnly = (1U << 31),
  kfMiscBiallelicOnlyStrict = (1LLU << 32),
  kfMiscBiallelicOnlyList = (1LLU << 33),
  kfMiscStrictSid0 = (1LLU << 34),
  kfMiscAllowBadFreqs = (1LLU << 35),
  kfMiscIidSid = (1LLU << 36),
  kfMiscPhenoIidOnly = (1LLU << 37),
  kfMiscCovarIidOnly = (1LLU << 38),
  kfMiscAllowBadLd = (1LLU << 39),
  kfMiscErrorOnFreqCalc = (1LLU << 40),
  kfMiscNoCategorical = (1LLU << 41),
  kfMiscAcFounders = (1LLU << 42),
  kfMiscMakeFoundersFirst = (1LLU << 43),
  kfMiscMakeFoundersNotfirst = (1LLU << 44),
  kfMiscMakeFoundersRequire2Missing = (1LLU << 45),
  kfMiscYNosexMissingStats = (1LLU << 46),
  kfMiscNeg9PhenoReallyMissing = (1LLU << 47),
  kfMiscAlt1Allele = (1LLU << 48),
  kfMiscSelectSidParentsOnly = (1LLU << 49)
FLAGSET64_DEF_END(MiscFlags);

FLAGSET64_DEF_START()
  kfExportf0,
  kfExportf01 = (1 << 0),
  kfExportf12 = (1 << 1),
  kfExportfSpaces = (1 << 2),
  kfExportfRefFirst = (1 << 3),
  kfExportf23 = (1 << 4),
  kfExportfA = (1 << 5),
  kfExportfAv = (1 << 6),
  kfExportfAD = (1 << 7),
  kfExportfBcf42 = (1 << 8),
  kfExportfBcf43 = (1 << 9),
  kfExportfBcf = kfExportfBcf42 | kfExportfBcf43,
  kfExportfBeagle = (1 << 10),
  kfExportfBeagleNomap = (1 << 11),
  kfExportfBgen11 = (1 << 12),
  kfExportfBgen12 = (1 << 13),
  kfExportfBgen13 = (1 << 14),
  kfExportfBimbam = (1 << 15),
  kfExportfBimbam1chr = (1 << 16),
  kfExportfEig = (1 << 17),
  kfExportfEigt = (1 << 18),
  kfExportfFastphase = (1 << 19),
  kfExportfFastphase1chr = (1 << 20),
  kfExportfHaps = (1 << 21),
  kfExportfHapsLegend = (1 << 22),
  kfExportfHv = (1 << 23),
  kfExportfHv1chr = (1 << 24),
  kfExportfIndMajorBed = (1 << 25),
  kfExportfLgen = (1 << 26),
  kfExportfLgenRef = (1 << 27),
  kfExportfList = (1 << 28),
  kfExportfRlist = (1 << 29),
  kfExportfOxGenV1 = (1 << 30),
  kfExportfOxGenV2 = (1U << 31),
  kfExportfOxGen = kfExportfOxGenV1 | kfExportfOxGenV2,
  kfExportfPed = (1LLU << 32),
  kfExportfCompound = (1LLU << 33),
  kfExportfPhylip = (1LLU << 34),
  kfExportfPhylipPhased = (1LLU << 35),
  kfExportfStructure = (1LLU << 36),
  kfExportfTped = (1LLU << 37),
  kfExportfVcf42 = (1LLU << 38),
  kfExportfVcf43 = (1LLU << 39),
  kfExportfVcf = kfExportfVcf42 | kfExportfVcf43,
  kfExportfTypemask = (2LLU * kfExportfVcf43) - kfExportf23,
  kfExportfIncludeAlt = (1LLU << 40),
  kfExportfBgz = (1LLU << 41),
  kfExportfOmitNonmaleY = (1LLU << 42),
  kfExportfSampleV2 = (1LLU << 43),
  kfExportfBgenOmitSampleIdBlock = (1LLU << 44),
  kfExportfPhylipUsedSites = (1LLU << 45)
FLAGSET64_DEF_END(ExportfFlags);

FLAGSET_DEF_START()
  kfInfo0,
  kfInfoPrFlagPresent = (1 << 0),
  kfInfoPrNonflagPresent = (1 << 1),
  kfInfoNonprPresent = (1 << 2),
  kfInfoAll = ((kfInfoNonprPresent * 2) - kfInfoPrFlagPresent),
  // this is set in the .bim case
  kfInfoPrNonrefDefault = (1 << 3),
FLAGSET_DEF_END(InfoFlags);

FLAGSET_DEF_START()
  kfPvarPsam0,
  kfPvarZs = (1 << 0),

  kfPvarColXheader = (1 << 1),
  kfPvarColVcfheader = (1 << 2),
  kfPvarColMaybequal = (1 << 3),
  kfPvarColQual = (1 << 4),
  kfPvarColMaybefilter = (1 << 5),
  kfPvarColFilter = (1 << 6),
  kfPvarColMaybeinfo = (1 << 7),
  kfPvarColInfo = (1 << 8),
  kfPvarColXinfo = (kfPvarColInfo * 2) - kfPvarColMaybeinfo,
  kfPvarColMaybecm = (1 << 9),
  kfPvarColCm = (1 << 10),
  kfPvarColDefault = (kfPvarColXheader | kfPvarColMaybequal | kfPvarColMaybefilter | kfPvarColMaybeinfo | kfPvarColMaybecm),
  kfPvarColAll = ((kfPvarColCm * 2) - kfPvarColXheader),
  kfPsamColMaybefid = (1 << 11),
  kfPsamColFid = (1 << 12),
  kfPsamColMaybesid = (1 << 13),
  kfPsamColSid = (1 << 14),
  kfPsamColMaybeparents = (1 << 15),
  kfPsamColParents = (1 << 16),
  kfPsamColSex = (1 << 17),
  kfPsamColPheno1 = (1 << 18),
  kfPsamColPhenos = (1 << 19),
  kfPsamColDefault = (kfPsamColMaybefid | kfPsamColMaybesid | kfPsamColMaybeparents | kfPsamColSex | kfPsamColPhenos),
  kfPsamColAll = ((kfPsamColPhenos * 2) - kfPsamColMaybefid)
FLAGSET_DEF_END(PvarPsamFlags);

// may want to rename FidPresent to FidMayBePresent
FLAGSET_DEF_START()
  kfSampleId0,
  kfSampleIdFidPresent = (1 << 0),
  kfSampleIdParentsPresent = (1 << 1),
  kfSampleIdNoIdHeader = (1 << 2),
  kfSampleIdNoIdHeaderIidOnly = (1 << 3),
  kfSampleIdStrictSid0 = (1 << 4)
FLAGSET_DEF_END(SampleIdFlags);

FLAGSET_DEF_START()
  kfLoadFilterLog0,
  // import-only
  kfLoadFilterLogImportMaxAlleles = (1 << 0),
  kfLoadFilterLogVcfRequireGt = (1 << 1),
  // merge-only
  kfLoadFilterLogMergeMaxAlleles = (1 << 2),
  // chromosome filters: always applied, but omitted from .pvar-load log if
  // import or merge already occurred
  kfLoadFilterLogAutosome = (1 << 3),
  kfLoadFilterLogAutosomePar = (1 << 4),
  kfLoadFilterLogChr = (1 << 5),
  kfLoadFilterLogNotChr = (1 << 6),
  kfLoadFilterLogImportMask = ((kfLoadFilterLogNotChr * 2) - kfLoadFilterLogMergeMaxAlleles - kfLoadFilterLogImportMaxAlleles),
  kfLoadFilterLogMergeMask = ((kfLoadFilterLogNotChr * 2) - kfLoadFilterLogMergeMaxAlleles),
  kfLoadFilterLogImportMergeAlreadyApplied = (1 << 7),
  kfLoadFilterLogImportMergeMask = ((kfLoadFilterLogImportMergeAlreadyApplied * 2) - kfLoadFilterLogImportMaxAlleles),
  // main .pvar load only
  kfLoadFilterLogExcludeIfInfo = (1 << 8),
  kfLoadFilterLogExcludePalindromicSnps = (1 << 9),
  kfLoadFilterLogExtractIfInfo = (1 << 10),
  kfLoadFilterLogMaxAlleles = (1 << 11),
  kfLoadFilterLogMinAlleles = (1 << 12),
  kfLoadFilterLogRequireInfo = (1 << 13),
  kfLoadFilterLogRequireNoInfo = (1 << 14),
  kfLoadFilterLogSnpsOnly = (1 << 15),
  kfLoadFilterLogVarFilter = (1 << 16),
  kfLoadFilterLogVarMinQual = (1 << 17)
FLAGSET_DEF_END(LoadFilterLogFlags);

// These structs are small enough and ownership of the pointed-to arrays is
// generally clear enough that the noncopyable annotation is just intended to
// be a tripwire.
typedef struct SampleIdInfoStruct {
  NONCOPYABLE(SampleIdInfoStruct);
  char* sample_ids;
  char* sids;
  uintptr_t max_sample_id_blen; // FID+IID only, does not count SID
  uintptr_t max_sid_blen;
  SampleIdFlags flags;
} SampleIdInfo;

typedef struct ParentalIdInfoStruct {
  NONCOPYABLE(ParentalIdInfoStruct);
  char* paternal_ids;
  char* maternal_ids;
  uintptr_t max_paternal_id_blen;
  uintptr_t max_maternal_id_blen;
} ParentalIdInfo;

typedef struct PedigreeIdInfoStruct {
  SampleIdInfo sii;
  ParentalIdInfo parental_id_info;
} PedigreeIdInfo;

typedef struct APermStruct {
  uint32_t min;
  uint32_t max;
  double alpha;
  double beta;
  double init_interval;
  double interval_slope;
} APerm;

// (2^31 - 1000001) / 2
CONSTI32(kApermMax, 1073241823);

typedef struct TwoColParamsStruct {
  NONCOPYABLE(TwoColParamsStruct);
  uint32_t colx;
  uint32_t colid;
  uint32_t skip_ct;
  char skipchar;
  char fname[];
} TwoColParams;

void InitPedigreeIdInfo(MiscFlags misc_flags, PedigreeIdInfo* piip);

// no CleanupPedigreeIdInfo function since LoadPsam() allocates the arrays in
// bigstack.

FLAGSET_DEF_START()
  kfFlip0,
  kfFlipPermissive = (1 << 0)
FLAGSET_DEF_END(FlipFlags);

typedef struct FlipInfoStruct {
  NONCOPYABLE(FlipInfoStruct);
  char* fname;
  char* subset_fname;
  FlipFlags flags;
} FlipInfo;

void InitFlip(FlipInfo* flip_info_ptr);

void CleanupFlip(FlipInfo* flip_info_ptr);


void AppendLoadFilterFlagnames(LoadFilterLogFlags load_filter_log_flags, char** write_iter_ptr);

HEADER_INLINE BoolErr bigstack_alloc_ac(uintptr_t ct, AlleleCode** allele_arr_ptr) {
  *allele_arr_ptr = S_CAST(AlleleCode*, bigstack_alloc(ct * sizeof(AlleleCode)));
  return !(*allele_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_kac(uintptr_t ct, const AlleleCode** allele_arr_ptr) {
  *allele_arr_ptr = S_CAST(const AlleleCode*, bigstack_alloc(ct * sizeof(AlleleCode)));
  return !(*allele_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_dosage(uintptr_t ct, Dosage** dosage_arr_ptr) {
  *dosage_arr_ptr = S_CAST(Dosage*, bigstack_alloc(ct * sizeof(Dosage)));
  return !(*dosage_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_dphase(uintptr_t ct, SDosage** dphase_arr_ptr) {
  *dphase_arr_ptr = S_CAST(SDosage*, bigstack_alloc(ct * sizeof(SDosage)));
  return !(*dphase_arr_ptr);
}


HEADER_INLINE BoolErr bigstack_alloc_acp(uintptr_t ct, AlleleCode*** acp_arr_ptr) {
  *acp_arr_ptr = S_CAST(AlleleCode**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*acp_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_dosagep(uintptr_t ct, Dosage*** dosagep_arr_ptr) {
  *dosagep_arr_ptr = S_CAST(Dosage**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*dosagep_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_dphasep(uintptr_t ct, SDosage*** dphasep_arr_ptr) {
  *dphasep_arr_ptr = S_CAST(SDosage**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*dphasep_arr_ptr);
}


HEADER_INLINE BoolErr bigstack_allocv_ac(uintptr_t ct, AlleleCode** allele_arr_ptr) {
  *allele_arr_ptr = S_CAST(AlleleCode*, bigstack_allocv(ct * sizeof(AlleleCode)));
  return !(*allele_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_allocv_dosage(uintptr_t ct, Dosage** dosage_arr_ptr) {
  *dosage_arr_ptr = S_CAST(Dosage*, bigstack_allocv(ct * sizeof(Dosage)));
  return !(*dosage_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_allocv_dphase(uintptr_t ct, SDosage** dphase_arr_ptr) {
  *dphase_arr_ptr = S_CAST(SDosage*, bigstack_allocv(ct * sizeof(SDosage)));
  return !(*dphase_arr_ptr);
}


HEADER_INLINE BoolErr bigstack_end_alloc_ac(uintptr_t ct, AlleleCode** allele_arr_ptr) {
  *allele_arr_ptr = S_CAST(AlleleCode*, bigstack_end_alloc(ct * sizeof(AlleleCode)));
  return !(*allele_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_dosage(uintptr_t ct, Dosage** dosage_arr_ptr) {
  *dosage_arr_ptr = S_CAST(Dosage*, bigstack_end_alloc(ct * sizeof(Dosage)));
  return !(*dosage_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_dphase(uintptr_t ct, SDosage** dphase_arr_ptr) {
  *dphase_arr_ptr = S_CAST(SDosage*, bigstack_end_alloc(ct * sizeof(SDosage)));
  return !(*dphase_arr_ptr);
}

BoolErr BigstackAllocPgv(uint32_t sample_ct, uint32_t multiallelic_needed, PgenGlobalFlags gflags, PgenVariant* pgvp);

// remainder must be in [1, 16383].
char* PrintDosageDecimal(uint32_t remainder, char* start);

// small_dosage must be in [0, kDosageMid * 10 - 1].
HEADER_INLINE char* PrintSmallDosage(uint32_t small_dosage, char* start) {
  *start++ = '0' + (small_dosage / kDosageMid);
  const uint32_t remainder = small_dosage % kDosageMid;
  if (!remainder) {
    return start;
  }
  return PrintDosageDecimal(remainder, start);
}

// 3 decimal places.  dosage on /kDosageMax rather than /kDosageMid scale
// (hence the extra 'd')
char* ddosagetoa(uint64_t dosage, char* start);

// remainder must be in [1, 32767].
char* PrintDdosageDecimal(uint32_t remainder, char* start);

// 5 decimal places.  Only used when it is important to be able to reconstruct
// the exact original value.
HEADER_INLINE char* ddosagetoa_full(uint64_t dosage, char* start) {
  start = u32toa(dosage / kDosageMax, start);
  const uint32_t remainder = dosage % kDosageMax;
  if (!remainder) {
    return start;
  }
  return PrintDdosageDecimal(remainder, start);
}

HEADER_INLINE void ZeroDosageArr(uintptr_t entry_ct, Dosage* dosage_arr) {
  memset(dosage_arr, 0, entry_ct * sizeof(Dosage));
}

HEADER_INLINE void ZeroDphaseArr(uintptr_t entry_ct, SDosage* dphase_arr) {
  memset(dphase_arr, 0, entry_ct * sizeof(SDosage));
}

HEADER_INLINE void SetAllDosageArr(uintptr_t entry_ct, Dosage* dosage_arr) {
  memset(dosage_arr, 255, entry_ct * sizeof(Dosage));
}

// Assumes dense_dosage is allocated up to vector boundary.
void PopulateDenseDosage(const uintptr_t* genoarr, const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, Dosage* dense_dosage);

void PopulateDenseDosageNonemptySubset(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* genoarr, const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t dosage_ct, Dosage* dense_dosage, uintptr_t* workspace);

void PopulateRescaledDosage(const uintptr_t* genoarr, const uintptr_t* dosage_present, const Dosage* dosage_main, double slope, double intercept, double missing_val, uint32_t sample_ct, uint32_t dosage_ct, double* expanded_dosages);

void PopulateRescaledDosageF(const uintptr_t* genoarr, const uintptr_t* dosage_present, const Dosage* dosage_main, float slope, float intercept, float missing_val, uint32_t sample_ct, uint32_t dosage_ct, float* expanded_dosages);

void PopulateDenseDphase(const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const uintptr_t* dosage_present, const Dosage* dense_dosage_vec, const uintptr_t* dphase_present, const SDosage* dphase_delta, uint32_t sample_ct, uint32_t phasepresent_ct, uint32_t dosage_ct, uint32_t dphase_ct, SDosage* dense_dphase_delta);

// assumes trailing bits of genoarr are zeroed out
HEADER_INLINE uint32_t AtLeastOneHetUnsafe(const uintptr_t* genoarr, uint32_t sample_ct) {
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  for (uint32_t uii = 0; uii != sample_ctl2; ++uii) {
    const uintptr_t geno_word = genoarr[uii];
    if (Word01(geno_word)) {
      return 1;
    }
  }
  return 0;
}

uint32_t AtLeastOneMultiallelicHet(const PgenVariant* pgvp, uint32_t sample_ct);

void SetHetMissing(uintptr_t word_ct, uintptr_t* genovec);

void SetHetMissingCleardosage(uintptr_t word_ct, uintptr_t* genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* dosagepresent, Dosage* dosage_main);

void SetHetMissingKeepdosage(uintptr_t word_ct, uintptr_t* genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* dosagepresent, Dosage* dosage_main);

void SetMissingRef(uintptr_t word_ct, uintptr_t* genovec);

void SetMissingRefY(const uintptr_t* __restrict sex_female_interleaved, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec);

void SetMissingRefDosage(const uintptr_t* __restrict dosagepresent, uintptr_t sample_ct, uintptr_t* __restrict genovec);

void SetMissingRefDosageY(const uintptr_t* __restrict sex_female, const uintptr_t* __restrict dosagepresent, uintptr_t sample_ct, uintptr_t* __restrict genovec);

uint32_t GenoarrCountMissingUnsafe(const uintptr_t* genoarr, uint32_t sample_ct);

uint32_t GenoarrCountMissingSubset(const uintptr_t* genoarr, const uintptr_t* interleaved_vec, uint32_t sample_ct);

uint32_t GenoarrCountMissingInvsubsetUnsafe(const uintptr_t* genoarr, const uintptr_t* exclude_mask, uint32_t sample_ct);

// See also DataFidColIsRequired(), etc. in plink2_data.h, which checks whether
// at least one remaining value is nonzero.
HEADER_INLINE uint32_t FidColIsRequired(const SampleIdInfo* siip, uint32_t maybe_modifier) {
  return (maybe_modifier & 2) || ((maybe_modifier & 1) && (siip->flags & kfSampleIdFidPresent));
}

HEADER_INLINE uint32_t SidColIsRequired(const char* sids, uint32_t maybe_modifier) {
  return (maybe_modifier & 2) || ((maybe_modifier & 1) && sids);
}

HEADER_INLINE uint32_t ParentalColsAreRequired(const PedigreeIdInfo* piip, uint32_t maybe_modifier) {
  return (maybe_modifier & 2) || ((maybe_modifier & 1) && (piip->sii.flags & kfSampleIdParentsPresent));
}

void CollapsedSampleFmtidInit(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t include_fid, uint32_t include_sid, uintptr_t max_sample_fmtid_blen, char* collapsed_sample_fmtids_iter);

HEADER_INLINE uintptr_t GetMaxSampleFmtidBlen(const SampleIdInfo* siip, uint32_t include_fid, uint32_t include_sid) {
  uintptr_t max_sid_blen = 0;
  if (include_sid) {
    max_sid_blen = siip->sids? siip->max_sid_blen : 2;
  }
  return siip->max_sample_id_blen + max_sid_blen - 2 * (!include_fid);
}

BoolErr CollapsedSampleFmtidInitAlloc(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t include_fid, uint32_t include_sid, char** collapsed_sample_fmtids_ptr, uintptr_t* max_sample_fmtid_blen_ptr);

// Assumes sample_ct > 0.
uint32_t OnlyOneFid(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct);

uint32_t GetMajIdxMulti(const double* cur_allele_freqs, uint32_t cur_allele_ct);

HEADER_INLINE uint32_t GetMajIdx(const double* cur_allele_freqs, uint32_t cur_allele_ct) {
  if (cur_allele_freqs[0] >= 0.5) {
    return 0;
  }
  if (cur_allele_ct == 2) {
    return 1;
  }
  return GetMajIdxMulti(cur_allele_freqs, cur_allele_ct);
}

HEADER_INLINE double GetNonmajFreq(const double* cur_allele_freqs, uint32_t allele_ct) {
  double tot_nonlast_freq = cur_allele_freqs[0];
  double max_freq = tot_nonlast_freq;
  const uint32_t allele_ct_m1 = allele_ct - 1;
  for (uint32_t allele_idx = 1; allele_idx != allele_ct_m1; ++allele_idx) {
    const double cur_alt_freq = cur_allele_freqs[allele_idx];
    tot_nonlast_freq += cur_alt_freq;
    if (cur_alt_freq > max_freq) {
      max_freq = cur_alt_freq;
    }
  }
  const double nonmajor_freq = 1.0 - max_freq;
  return MINV(nonmajor_freq, tot_nonlast_freq);
}

HEADER_INLINE double GetAlleleFreq(const double* cur_allele_freqs, uint32_t allele_idx, uint32_t allele_ct) {
  const uint32_t allele_ct_m1 = allele_ct - 1;
  if (allele_idx < allele_ct_m1) {
    return cur_allele_freqs[allele_idx];
  }
  double last_freq = 1.0 - cur_allele_freqs[0];
  for (uint32_t tmp_allele_idx = 1; tmp_allele_idx != allele_ct_m1; ++tmp_allele_idx) {
    last_freq -= cur_allele_freqs[tmp_allele_idx];
  }
  return MAXV(last_freq, 0.0);
}


// Functions with Xid in their name deal with both FID/IID (with a single tab
// separator) and FID/IID/SID (two tabs) sample IDs.  (Missing FID is
// represented as FID='0'.)

// With no header line, --keep/--remove and the like should interpret a single
// token as IID (treating FID as 0), and two tokens as FID/IID.  --no-id-header
// does not support IID/SID output, so we don't need to worry about supporting
// FidIidSidOrIidSid.
// With a header line, all four {FID present/absent, SID present/absent}
// combinations are allowed.
// bugfix (26 Dec 2020): XidRead() now skips the SID column when appropriate,
// and GetXidColCt() counts it.
FLAGSET_DEF_START()
  kfXidMode0,

  kfXidModeFlagOneCoreTokenOk = (1 << 0),
  kfXidModeFlagNeverFid = (1 << 1),
  kfXidModeFlagSid = (1 << 2),

  kfXidModeFlagSkipSid = (1 << 3),

  kfXidModeFidIid = 0,
  kfXidModeFidIidOrIid = kfXidModeFlagOneCoreTokenOk,
  kfXidModeIid = (kfXidModeFlagOneCoreTokenOk | kfXidModeFlagNeverFid),
  kfXidModeFidIidSid = kfXidModeFlagSid,
  kfXidModeIidSid = (kfXidModeFlagNeverFid | kfXidModeFlagSid),
  kfXidModeCoreMask = kfXidModeFlagSkipSid - 1
FLAGSET_DEF_END(XidMode);

// Assumes fixed-width.  Includes skipped SID column if there is one.
HEADER_INLINE uint32_t GetXidColCt(XidMode xid_mode) {
  const XidMode xid_mode_core = xid_mode & kfXidModeCoreMask;
  return 2 + ((xid_mode >> 3) & 1) + (xid_mode_core == kfXidModeFidIidSid) - (xid_mode_core == kfXidModeIid);
}

// sample_xid_map allocated on bottom, to play well with --indiv-sort
PglErr SortedXidboxInitAlloc(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, XidMode xid_mode, uint32_t use_nsort, char** sorted_xidbox_ptr, uint32_t** xid_map_ptr, uintptr_t* max_xid_blen_ptr);

PglErr SortedXidboxInitAllocEnd(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, XidMode xid_mode, uint32_t use_nsort, char** sorted_xidbox_ptr, uint32_t** xid_map_ptr, uintptr_t* max_xid_blen_ptr);

// returns slen for ID, or 0 on guaranteed mismatch (longer than max_xid_blen)
// or parse failure (*readpp set to nullptr in latter case).
uint32_t XidRead(uintptr_t max_xid_blen, uint32_t comma_delim, XidMode xid_mode, const char** read_pp, char* __restrict idbuf);

// returns 1 on missing token *or* if the sample ID is not present.  cases can
// be distinguished by checking whether *read_pp_new == nullptr: if it is, a
// missing-tokens error should probably be reported.
// sample_id_map == nullptr is permitted
// *read_pp is now set to point to the end of IID/SID instead of the beginning
// of the next token; this is a change from plink 1.9.
HEADER_INLINE BoolErr SortedXidboxReadFind(const char* __restrict sorted_xidbox, const uint32_t* __restrict xid_map, uintptr_t max_xid_blen, uintptr_t xid_ct, uint32_t comma_delim, XidMode xid_mode, const char** read_pp, uint32_t* sample_uidx_ptr, char* __restrict idbuf) {
  const uint32_t slen_final = XidRead(max_xid_blen, comma_delim, xid_mode, read_pp, idbuf);
  idbuf[slen_final] = '\0'; // needed for some error messages
  if (!slen_final) {
    return 1;
  }
  return SortedIdboxFind(idbuf, sorted_xidbox, xid_map, slen_final, max_xid_blen, xid_ct, sample_uidx_ptr);
}

// Matches a sample ID *prefix*.  Thus, if FID/IID/SID is loaded, but the input
// file contains just FID/IID, and there are some FID/IID pairs which
// correspond to multiple samples, this lets you iterate over all of them.
// (Caller is responsible for looking up xid_map[] to perform xid_idx ->
// sample_uidx conversions.)
BoolErr SortedXidboxReadMultifind(const char* __restrict sorted_xidbox, uintptr_t max_xid_blen, uintptr_t xid_ct, uint32_t comma_delim, XidMode xid_mode, const char** read_pp, uint32_t* __restrict xid_idx_start_ptr, uint32_t* __restrict xid_idx_end_ptr, char* __restrict idbuf);

FLAGSET_DEF_START()
  kfXidHeader0,

  kfXidHeaderFixedWidth = (1 << 0),
  kfXidHeaderIgnoreSid = (1 << 1),
  kfXidHeaderFixedWidthIgnoreSid = (kfXidHeaderFixedWidth | kfXidHeaderIgnoreSid)
FLAGSET_DEF_END(XidHeaderFlags);

// - May return kPglRetEof, or other TextStream errors.
// - line_startp can be nullptr.  If it isn't, it's set to the (lstripped)
//   beginning of the current line.
// - line_idx must be zero unless initial lines were skipped.
// - If no header line is present, xid_mode will be set to kfXidModeFidIid if
//   kfXidHeaderFixedWidth is set, and kfXidModeFidIidOrIid (which tolerates a
//   mix of single-token and multitoken lines) otherwise.
// - TSTREAM_FAIL errors are not printed, but other errors are.
PglErr LoadXidHeader(const char* flag_name, XidHeaderFlags xid_header_flags, uintptr_t* line_idx_ptr, TextStream* txsp, XidMode* xid_mode_ptr, char** line_startp);

PglErr OpenAndLoadXidHeader(const char* fname, const char* flag_name, XidHeaderFlags xid_header_flags, uint32_t max_line_blen, TextStream* txsp, XidMode* xid_mode_ptr, uintptr_t* line_idx_ptr, char** line_startp);

// header line expected to start with FID1, ID1, or IID1
PglErr LoadXidHeaderPair(const char* flag_name, uint32_t sid_over_fid, uintptr_t* line_idx_ptr, TextStream* txsp, XidMode* xid_mode_ptr, char** line_startp, char** line_iterp);

// sorted_xidbox must be created with kfXidModeFidIid or kfXidModeFidIidSid,
// and dst must initially be zeroed out.
// does not check for graph cycles, etc.
void MarkParents(const ParentalIdInfo* parental_id_infop, const char* sorted_xidbox, const uint32_t* xid_map, uint32_t sample_ct, uintptr_t max_xid_blen, uint32_t use_nsort, uintptr_t* dst, char* idbuf);

// Assumes no duplicates.
void InitXidHtable(const SampleIdInfo* siip, uint32_t sample_ct, uint32_t xid_htable_size, uint32_t* xid_htable, char* idbuf);

BoolErr LookupXidHtable(const char* line_start, const SampleIdInfo* siip, const uint32_t* xid_htable, uint32_t xid_htable_size, uint32_t fid_present, uint32_t sid_present, uint32_t* sample_idxp, char* idbuf);

PglErr CheckXidUniqueness(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* err_suffix_str, uint32_t sample_ct);

// could add MultifileIntersection mode if there's a user
FLAGSET_DEF_START()
  kfLoadSampleIds0,
  kfLoadSampleIdsMultifile = (1 << 0),
  kfLoadSampleIdsFamOnly = (1 << 1)
FLAGSET_DEF_END(LoadSampleIdsFlags);

// Allocates (at the bottom of bigstack) and returns results in loaded_bitarr.
PglErr LoadSampleIds(const char* fnames, const uintptr_t* sample_include, const SampleIdInfo* siip, const char* flag_name, uint32_t raw_sample_ct, uint32_t sample_ct, LoadSampleIdsFlags flags, uintptr_t** loaded_bitarrp, uint32_t* duplicate_ctp);

// AppendXid is inline since it's reasonable to call it in loops writing
// multi-GB files.  AppendSpacedXid is not since it's only called when printing
// error messages.
HEADER_INLINE char* AppendXid(const char* sample_ids, const char* sids, uint32_t write_fid, uint32_t write_sid, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t sample_uidx, char* write_iter) {
  const char* cur_sample_id = &(sample_ids[max_sample_id_blen * sample_uidx]);
  if (!write_fid) {
    cur_sample_id = AdvPastDelim(cur_sample_id, '\t');
  }
  write_iter = strcpya(write_iter, cur_sample_id);
  if (write_sid) {
    *write_iter++ = '\t';
    if (sids) {
      write_iter = strcpya(write_iter, &(sids[max_sid_blen * sample_uidx]));
    } else {
      *write_iter++ = '0';
    }
  }
  return write_iter;
}

char* AppendSpacedXid(const char* sample_ids, const char* sids, uint32_t write_fid, uint32_t write_sid, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t sample_uidx, char* write_iter);


extern const unsigned char g_char_to_sex[256];

HEADER_INLINE uint32_t CharToSex(char cc) {
  return g_char_to_sex[S_CAST(unsigned char, cc)];
}


#ifndef HIGH_CONTIG_BUILD
// note that this is no longer divisible by 64
#  ifndef HIGH_AUTOSOME_NUM_BUILD
CONSTI32(kMaxContigs, 65274);
#  else
CONSTI32(kMaxContigs, 64506);
#  endif
CONSTI32(kMaxChrCodeDigits, 5);

// change ChrIdx to uint32_t if (kMaxContigs + kChrOffsetCt) > 65536
typedef uint16_t ChrIdx;

// GetHtableMinSize(kChrRawEnd) (use constexpr once sufficient compiler support
// is available)
// (not GetHtableFastSize since, an overwhelming majority of the time, we'll
// have far fewer than 2^16 codes)
CONSTI32(kChrHtableSize, 130560);
#else
// Keep this number in sync with the top-of-file comment, and the error message
// in TryToAddChrName().
#  ifndef HIGH_AUTOSOME_NUM_BUILD
CONSTI32(kMaxContigs, 982778);
#  else
CONSTI32(kMaxContigs, 982010);
#  endif
CONSTI32(kMaxChrCodeDigits, 6);

typedef uint32_t ChrIdx;

CONSTI32(kChrHtableSize, 1965568);
#endif

static_assert((kMaxChrCodeDigits == 5) || (kMaxChrCodeDigits == 6), "u32toa_zchr() must be updated");
HEADER_INLINE char* u32toa_zchr(uint32_t uii, char* start) {
  if (kMaxChrCodeDigits == 5) {
    return u32toa_z5(uii, start);
  } else {
    return u32toa_z6(uii, start);
  }
}

HEADER_INLINE BoolErr bigstack_alloc_chridx(uintptr_t ct, ChrIdx** chridx_arr_ptr) {
  *chridx_arr_ptr = S_CAST(ChrIdx*, bigstack_alloc(ct * sizeof(ChrIdx)));
  return !(*chridx_arr_ptr);
}

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
ENUM_U31_DEF_END(XymtOffset);

CONSTI32(kChrRawX, kMaxContigs + kChrOffsetX);
CONSTI32(kChrRawY, kMaxContigs + kChrOffsetY);
CONSTI32(kChrRawXY, kMaxContigs + kChrOffsetXY);
CONSTI32(kChrRawMT, kMaxContigs + kChrOffsetMT);
CONSTI32(kChrRawPAR1, kMaxContigs + kChrOffsetPAR1);
CONSTI32(kChrRawPAR2, kMaxContigs + kChrOffsetPAR2);
CONSTI32(kChrRawEnd, kMaxContigs + kChrOffsetCt);

static_assert((!(kChrRawEnd % kBitsPerWord)), "kChrRawEnd expression must be updated.");
CONSTI32(kChrMaskWords, kChrRawEnd / kBitsPerWord);

#ifndef HIGH_AUTOSOME_NUM_BUILD
// (note that n+1, n+2, n+3, and n+4 are reserved for X/Y/XY/MT)
CONSTI32(kMaxChrTextnum, 95);

#  ifdef __LP64__
CONSTI32(kChrExcludeWords, 2);
#  else
CONSTI32(kChrExcludeWords, 4);
#  endif

// GetChrCodeRaw() needs to be modified if this changes
CONSTI32(kMaxChrTextnumSlen, 2);
#else  // HIGH_AUTOSOME_NUM_BUILD
CONSTI32(kMaxChrTextnum, 995);

#  ifdef __LP64__
CONSTI32(kChrExcludeWords, 16);
#  else
CONSTI32(kChrExcludeWords, 32);
#  endif

CONSTI32(kMaxChrTextnumSlen, 3);
#endif
static_assert(kChrExcludeWords * kBitsPerWord >= kMaxChrTextnum + 2 * kChrOffsetCt + 1, "kChrExcludeWords must be updated.");

// AnotherFile is for .pvar merge (chrSets must be consistent across all .pvar
// files in that case)
ENUM_U31_DEF_START()
  kChrsetSourceDefault,
  kChrsetSourceCmdline,
  kChrsetSourceFile,
  kChrsetSourceAnotherFile
ENUM_U31_DEF_END(ChrsetSource);

FLAGSET_DEF_START()
  kfChrOutput0,
  kfChrOutputPrefix = (1 << 0),
  kfChrOutputM = (1 << 1),
  kfChrOutputMT = (1 << 2),
  kfChrOutput0M = (1 << 3)
FLAGSET_DEF_END(ChrOutput);

typedef struct ChrInfoStruct {
  // Main dynamic block intended to be allocated as a single aligned block of
  // memory on the heap freeable with vecaligned_free(), with chr_mask at the
  // base.
  NONCOPYABLE(ChrInfoStruct);

  uintptr_t* chr_mask;  // which chromosomes aren't known to be absent?

  // This normally includes chrX.  However, library function implementations
  // cannot assume that: e.g. --glm now temporarily clears the chrX bit when
  // all samples are female.
  // As of alpha 2, it also includes MT again (like plink 1.07, and unlike 1.9
  // and 2.0a1), now that enough dosage functionality is in place.
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
  const char** nonstd_names;
  uint32_t* nonstd_id_htable;
  // end main dynamic block

  uint32_t chr_ct;  // number of distinct chromosomes/contigs
  ChrsetSource chrset_source;

  uintptr_t chr_exclude[kChrExcludeWords];

  // X, Y, XY...; UINT32_MAXM1 = not in chromosome set
  STD_ARRAY_DECL(uint32_t, kChrOffsetCt, xymt_codes);

  uint32_t max_numeric_code;
  uint32_t max_code;  // no longer identical to max_numeric_code, with PARs

  uint32_t autosome_ct;

  // yet more --allow-extra-chr support
  uint32_t zero_extra_chrs;
  uint32_t name_ct;
  LlStr* incl_excl_name_stack;
  uint32_t is_include_stack;
  ChrOutput output_encoding;
} ChrInfo;

extern const char g_xymt_log_names[kChrOffsetCt][5];

PglErr InitChrInfo(ChrInfo* cip);

void FinalizeChrset(LoadFilterLogFlags load_filter_log_flags, ChrInfo* cip);

HEADER_INLINE PglErr InitChrInfoHuman(ChrInfo* cip) {
  // convenience wrapper
  if (unlikely(InitChrInfo(cip))) {
    return kPglRetNomem;
  }
  FinalizeChrset(kfLoadFilterLog0, cip);
  return kPglRetSuccess;
}

HEADER_INLINE uint32_t IsHumanChrset(const ChrInfo* cip) {
  return (cip->max_code == 28) && (cip->max_numeric_code == 26) && (cip->autosome_ct == 22);
}

void ForgetExtraChrNames(uint32_t reinitialize, ChrInfo* cip);

// in the usual case where the number of chromosomes/contigs is much less than
// kMaxContigs, this reduces chr_info's memory consumption and improves
// locality.
PglErr FinalizeChrInfo(ChrInfo* cip);

void CleanupChrInfo(ChrInfo* cip);

// --export may need to override cip->output_encoding
char* ChrNameStdEx(const ChrInfo* cip, uint32_t chr_idx, ChrOutput output_encoding, char* buf);

HEADER_INLINE char* ChrNameStd(const ChrInfo* cip, uint32_t chr_idx, char* buf) {
  return ChrNameStdEx(cip, chr_idx, cip->output_encoding, buf);
}

// assumes chr_idx is valid
// note that chr_idx == 0 is always rendered as '0', never 'chr0'
char* chrtoa(const ChrInfo* cip, uint32_t chr_idx, char* buf);

uint32_t GetMaxChrSlen(const ChrInfo* cip);

uint32_t IsHaploidChrPresent(const ChrInfo* cip);

uint32_t IsAutosomalDiploidChrPresent(const ChrInfo* cip);

// any character <= ' ' is considered a terminator
// maps chrX -> kChrRawX, etc.
// now returns UINT32_MAXM1 on too many numeric digits
uint32_t GetChrCodeRaw(const char* str_iter);

// requires chr_name to be null-terminated
// maps chrX -> xymt_codes[kChrOffsetX], etc.
// error codes:
//   UINT32_MAX = --allow-extra-chr ok
//   UINT32_MAXM1 = total fail
uint32_t GetChrCode(const char* chr_name, const ChrInfo* cip, uint32_t name_slen);

// when the chromosome name isn't null-terminated
// requires chr_name[name_slen] to be mutable
uint32_t GetChrCodeCounted(const ChrInfo* cip, uint32_t name_slen, char* chr_name);

void ChrError(const char* chr_name, const char* file_descrip, const ChrInfo* cip, uintptr_t line_idx, uint32_t error_code);

HEADER_INLINE uint32_t GetVariantChrFoIdx(const ChrInfo* cip, uintptr_t variant_uidx) {
  return LowerBoundNonemptyU32(&(cip->chr_fo_vidx_start[1]), cip->chr_ct, variant_uidx + 1);
}

HEADER_INLINE uint32_t GetVariantChr(const ChrInfo* cip, uintptr_t variant_uidx) {
  return cip->chr_file_order[GetVariantChrFoIdx(cip, variant_uidx)];
}

HEADER_INLINE uint32_t XymtExists(const ChrInfo* cip, uint32_t xymt_offset, uint32_t* xymt_code_ptr) {
  // too easy to forget IsSet(chr_mask) check if we don't use this
  const uint32_t xymt_code = cip->xymt_codes[xymt_offset];
  *xymt_code_ptr = xymt_code;
  return (!IsI32Neg(xymt_code)) && IsSet(cip->chr_mask, xymt_code);
}

HEADER_INLINE void GetXymtStartAndEnd(const ChrInfo* cip, uint32_t xymt_offset, uint32_t* xymt_start_ptr, uint32_t* xymt_end_ptr) {
  uint32_t xymt_code;
  if (!XymtExists(cip, xymt_offset, &xymt_code)) {
    *xymt_start_ptr = 0;
    *xymt_end_ptr = 0;
    return;
  }
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[xymt_code];
  *xymt_start_ptr = cip->chr_fo_vidx_start[chr_fo_idx];
  *xymt_end_ptr = cip->chr_fo_vidx_start[chr_fo_idx + 1];
}

HEADER_INLINE void GetXymtCodeStartAndEndUnsafe(const ChrInfo* cip, uint32_t xymt_offset, uint32_t* xymt_code_ptr, uint32_t* xymt_start_ptr, uint32_t* xymt_end_ptr) {
  // assumes XymtExists was previously called, and is true
  const uint32_t xymt_code = cip->xymt_codes[xymt_offset];
  *xymt_code_ptr = xymt_code;
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[xymt_code];
  *xymt_start_ptr = cip->chr_fo_vidx_start[chr_fo_idx];
  *xymt_end_ptr = cip->chr_fo_vidx_start[chr_fo_idx + 1];
}

// now assumes chr_name is null-terminated
PglErr TryToAddChrName(const char* chr_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t prohibit_extra_chrs, uint32_t* chr_idx_ptr, ChrInfo* cip);

HEADER_INLINE PglErr GetOrAddChrCode(const char* chr_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t prohibit_extra_chrs, ChrInfo* cip, uint32_t* chr_idx_ptr) {
  *chr_idx_ptr = GetChrCode(chr_name, cip, name_slen);
  if (!IsI32Neg(*chr_idx_ptr)) {
    return kPglRetSuccess;
  }
  return TryToAddChrName(chr_name, file_descrip, line_idx, name_slen, prohibit_extra_chrs, chr_idx_ptr, cip);
}

HEADER_INLINE PglErr GetOrAddChrCodeDestructive(const char* file_descrip, uintptr_t line_idx, uint32_t prohibit_extra_chrs, char* chr_name, char* chr_name_end, ChrInfo* cip, uint32_t* chr_idx_ptr) {
  *chr_name_end = '\0';
  return GetOrAddChrCode(chr_name, file_descrip, line_idx, chr_name_end - chr_name, prohibit_extra_chrs, cip, chr_idx_ptr);
}

// Assumes sample_ct positive.  Does not require trailing bits to be clear.
uint32_t AllGenoEqual(const uintptr_t* genoarr, uint32_t sample_ct);

void InterleavedInvert(const uintptr_t* __restrict interleaved_set, uint32_t sample_ctv2, uintptr_t* __restrict genovec);

uint32_t InterleavedSubsetIs03(const uintptr_t* genovec, const uintptr_t* interleaved_set, uint32_t sample_ctv2);

uint32_t InterleavedSubsetIs23(const uintptr_t* genovec, const uintptr_t* interleaved_set, uint32_t sample_ctv2);

uint32_t InterleavedSubsetAllMissing(const uintptr_t* genovec, const uintptr_t* interleaved_set, uint32_t sample_ctv2);

// zeroes out samples not in the mask
void InterleavedMaskZero(const uintptr_t* __restrict interleaved_mask, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec);

// sets samples outside the mask to missing (0b11), trailing bits are dirtied
// void InterleavedMaskMissing(const uintptr_t* __restrict interleaved_set, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec);

// sets samples in the mask to missing (0b11)
void InterleavedSetMissing(const uintptr_t* __restrict interleaved_set, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec);

void InterleavedSetMissingCleardosage(const uintptr_t* __restrict orig_set, const uintptr_t* __restrict interleaved_set, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main);

void EraseDosages(const uintptr_t* __restrict erase_map, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main);

HEADER_INLINE void EraseDphases(const uintptr_t* __restrict erase_map, uint32_t* __restrict write_dphase_ct_ptr, uintptr_t* __restrict dphasepresent, SDosage* dphase_delta) {
  EraseDosages(erase_map, write_dphase_ct_ptr, dphasepresent, R_CAST(Dosage*, dphase_delta));
}

void SetMaleHetMissing(const uintptr_t* __restrict sex_male_interleaved, uint32_t geno_vec_ct, uintptr_t* __restrict genovec);

void EraseMaleDphases(const uintptr_t* __restrict sex_male, uint32_t* __restrict write_dphase_ct_ptr, uintptr_t* __restrict dphasepresent, SDosage* dphase_delta);

void SetMaleHetMissingCleardosage(const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_male_interleaved, uint32_t geno_vec_ct, uintptr_t* __restrict genovec, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main);

void SetMaleHetMissingKeepdosage(const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_male_interleaved, uint32_t geno_word_ct, uintptr_t* __restrict genovec, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main);

// Clears each bit in bitarr which doesn't correspond to a genovec het.
// Assumes that either trailing bits of bitarr are already zero, or trailing
// bits of genovec are zero.
HEADER_INLINE void MaskGenoarrHetsUnsafe(const uintptr_t* __restrict genoarr, uint32_t raw_sample_ctl2, uintptr_t* __restrict bitarr) {
  MaskWordsToHalfwordsInvmatch(genoarr, kMaskAAAA, raw_sample_ctl2, bitarr, bitarr);
}

void MaskGenoarrHetsMultiallelicUnsafe(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t raw_sample_ctl2, uintptr_t* __restrict bitarr);

// These require rare10_ct > 0.
void SetAltxyHetMissing(uintptr_t* __restrict genoarr, uint32_t* __restrict rare10_ct_ptr, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals);

void SetMaleAltxyHetMissing(const uintptr_t* sex_male, uintptr_t* __restrict genoarr, uint32_t* __restrict rare10_ct_ptr, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals);

// vertical popcount support
// VcountScramble1() and friends in plink2_cmdline
#ifdef __LP64__
#  ifdef USE_AVX2
// 2->4: 0 2 ... 126 1 3 ... 127
// 4->8: 0 4 ... 124 2 6 ... 126 1 5 ... 125 3 7 ... 127
// 8->32: 0 16 ... 112 4 20 ... 116 ... 124 2 18 ... 114 6 22 ... 118 ... 126 1 17 ...
HEADER_INLINE uint32_t VcountScramble2(uint32_t orig_idx) {
  return (orig_idx & (~127)) + ((orig_idx & 1) * 64) + ((orig_idx & 2) * 16) + ((orig_idx & 12) * 2) + ((orig_idx & 112) / 16);
}
#  else
HEADER_INLINE uint32_t VcountScramble2(uint32_t orig_idx) {
  return (orig_idx & (~63)) + ((orig_idx & 1) * 32) + ((orig_idx & 2) * 8) + (orig_idx & 12) + ((orig_idx & 48) / 16);
}
#  endif
#else
// 2->4: 0 2 4 6 8 10 12 14 1 3 5 ...
// 4->8: 0 4 8 12 2 6 10 14 1 5 9 ...
// 8->32: 0 4 8 12 2 6 10 14 1 5 9 13 3 7 11 15
HEADER_INLINE uint32_t VcountScramble2(uint32_t orig_idx) {
  return (orig_idx & (~15)) + ((orig_idx & 1) * 8) + ((orig_idx & 2) * 2) + ((orig_idx & 12) / 4);
}
#endif

HEADER_INLINE void VcountIncr2To4(const VecW* acc2_iter, uint32_t acc2_vec_ct, VecW* acc4_iter) {
  const VecW m2 = VCONST_W(kMask3333);
  for (uint32_t vidx = 0; vidx != acc2_vec_ct; ++vidx) {
    VecW loader = *acc2_iter++;
    *acc4_iter += loader & m2;
    ++acc4_iter;
    *acc4_iter += vecw_srli(loader, 2) & m2;
    ++acc4_iter;
  }
}

HEADER_INLINE void Vcount0Incr2To4(uint32_t acc2_vec_ct, VecW* acc2_iter, VecW* acc4_iter) {
  const VecW m2 = VCONST_W(kMask3333);
  for (uint32_t vidx = 0; vidx != acc2_vec_ct; ++vidx) {
    VecW loader = *acc2_iter;
    *acc2_iter++ = vecw_setzero();
    *acc4_iter += loader & m2;
    ++acc4_iter;
    *acc4_iter += vecw_srli(loader, 2) & m2;
    ++acc4_iter;
  }
}


// uint32_t chr_window_max(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_pos, uint32_t chr_fo_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max);

HEADER_INLINE uint32_t CountChrVariantsUnsafe(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t chr_idx) {
  assert(IsSet(cip->chr_mask, chr_idx));
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
  const uint32_t start_vidx = cip->chr_fo_vidx_start[chr_fo_idx];
  const uint32_t end_vidx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  return PopcountBitRange(variant_include, start_vidx, end_vidx);
}

HEADER_INLINE uint32_t ChrIsNonempty(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t chr_idx) {
  if (!IsSet(cip->chr_mask, chr_idx)) {
    return 0;
  }
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
  const uint32_t start_vidx = cip->chr_fo_vidx_start[chr_fo_idx];
  const uint32_t end_vidx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  return !AllBitsAreZero(variant_include, start_vidx, end_vidx);
}

HEADER_INLINE uint32_t XymtIsNonempty(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t xymt_offset) {
  const uint32_t xymt_code = cip->xymt_codes[xymt_offset];
  if (IsI32Neg(xymt_code) || (!IsSet(cip->chr_mask, xymt_code))) {
    return 0;
  }
  return ChrIsNonempty(variant_include, cip, xymt_code);
}

// assumes there's at least one variant on specified chromosome
uint32_t NotOnlyXymt(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t raw_variant_ct, uint32_t xymt_offset);

uint32_t CountNonAutosomalVariants(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t count_x, uint32_t count_mt);

void ExcludeNonAutosomalVariants(const ChrInfo* cip, uintptr_t* variant_include);

// If calc_descrip == nullptr, no logging occurs when some variants are
// removed, and no error occurs if all variants are removed.  (Might decouple
// these two later.)
PglErr ConditionalAllocateNonAutosomalVariants(const ChrInfo* cip, const char* calc_descrip, uint32_t raw_variant_ct, const uintptr_t** variant_include_ptr, uint32_t* variant_ct_ptr);

void FillSubsetChrFoVidxStart(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t* subset_chr_fo_vidx_start);

HEADER_INLINE BoolErr AllocAndFillSubsetChrFoVidxStart(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t** subset_chr_fo_vidx_start_ptr) {
  const uint32_t chr_ct = cip->chr_ct;
  if (unlikely(bigstack_alloc_u32(chr_ct + 1, subset_chr_fo_vidx_start_ptr))) {
    return 1;
  }
  FillSubsetChrFoVidxStart(variant_include, cip, *subset_chr_fo_vidx_start_ptr);
  return 0;
}

/*
// newval does not need to be null-terminated
// assumes *allele_ptr is not initialized
// (stop using these in main plink2 binary?)
BoolErr allele_set(const char* newval, uint32_t allele_slen, char** allele_ptr);

// *allele_ptr must be initialized; frees *allele_ptr if necessary
BoolErr allele_reset(const char* newval, uint32_t allele_slen, char** allele_ptr);

void cleanup_allele_storage(uint32_t max_allele_slen, uintptr_t allele_storage_entry_ct, const char** allele_storage);
*/

CONSTI32(kMaxMissingPhenostrBlen, 32);

// don't care about kfUnsortedVarChrom
FLAGSET_DEF_START()
  kfUnsortedVar0,
  kfUnsortedVarBp = (1 << 0),
  kfUnsortedVarCm = (1 << 1),
  kfUnsortedVarSplitChr = (1 << 2)
FLAGSET_DEF_END(UnsortedVar);

FLAGSET_DEF_START()
  kfFamCol0,
  kfFamCol1 = (1 << 0),
  kfFamCol34 = (1 << 1),
  kfFamCol5 = (1 << 2),
  kfFamCol6 = (1 << 3),
  kfFamCol13456 = (kfFamCol1 | kfFamCol34 | kfFamCol5 | kfFamCol6)
FLAGSET_DEF_END(FamCol);

HEADER_INLINE char Sexchar(const uintptr_t* sex_nm, const uintptr_t* sex_male, uintptr_t sample_uidx) {
  if (IsSet(sex_nm, sample_uidx)) {
    return '2' - IsSet(sex_male, sample_uidx);
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
ENUM_U31_DEF_END(PhenoDtype);

typedef union {
  uintptr_t* cc;  // bitvector
  double* qt;
  uint32_t* cat;  // always 0 for missing, nonmiss[] check unnecessary
} PhenoData;

typedef struct PhenoColStruct {
  // * If categorical phenotype, [0] points to g_missing_catname, while [1],
  //   [2], etc. point to category names.  These are part of the same
  //   allocation as nonmiss, so no separate free is needed.
  //   Otherwise, this is nullptr.
  // * When .sample non-V2 categorical variables are imported, 'C' is added in
  //   front of the integers.
  // * category_names[1..(n-1)] is guaranteed to be in natural-sorted order.
  // (MergePsams() has slightly different behavior since its pheno_cols
  // instance is only used for a one-time internal write.)
  const char** category_names;

  uintptr_t* nonmiss;  // bitvector

  // essentially a tagged union; part of the same allocation as nonmiss
  PhenoData data;
  PhenoDtype type_code;

  uint32_t nonnull_category_ct;
} PhenoCol;

uint32_t IsCategoricalPhenostr(const char* phenostr_iter);

uint32_t IsCategoricalPhenostrNocsv(const char* phenostr_iter);

// returns 0xffffffffU if none exists
uint32_t FirstCcOrQtPhenoIdx(const PhenoCol* pheno_cols, uint32_t pheno_ct);

// "_covar" since this doesn't handle case/control
uint32_t IsConstCovar(const PhenoCol* covar_col, const uintptr_t* sample_include, uint32_t raw_sample_ct);

// returns number of nonempty categories, null included
uint32_t IdentifyRemainingCats(const uintptr_t* sample_include, const PhenoCol* covar_col, uint32_t sample_ct, uintptr_t* observed_cat_bitarr);

// returns index of most common category
uint32_t IdentifyRemainingCatsAndMostCommon(const uintptr_t* sample_include, const PhenoCol* covar_col, uint32_t sample_ct, uintptr_t* observed_cat_bitarr, uint32_t* cat_obs_buf);

uint32_t GetCatSamples(const uintptr_t* sample_include_base, const PhenoCol* cat_pheno_col, uint32_t raw_sample_ctl, uint32_t sample_ct, uint32_t cat_uidx, uintptr_t* cur_cat_samples);

uint32_t RemoveExcludedCats(const uint32_t* data_cat, const uintptr_t* cat_keep_bitarr, uint32_t raw_sample_ct, uint32_t in_sample_ct, uintptr_t* sample_include);

// pheno_names is also allocated on the heap, but it can be handled with a
// simple free_cond().
void CleanupPhenoCols(uint32_t pheno_ct, PhenoCol* pheno_cols);

PglErr ParseChrRanges(const char* const* argvk, const char* flagname_p, const char* errstr_append, uint32_t param_ct, uint32_t prohibit_extra_chrs, uint32_t xymt_subtract, char range_delim, ChrInfo* cip, uintptr_t* chr_mask);


uint32_t MultiallelicVariantPresent(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t variant_ct);

uint32_t CountBiallelicVariants(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t variant_ct);

uintptr_t CountAlleles(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t variant_ct);

// This can be an overestimate (it won't return a value smaller than
// read_block_size).
uintptr_t GetMaxAltAlleleBlockSize(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t read_block_size);

// Mhc stands for "multiallelic hardcall" here, not major histocompatibility
// complex.  Though multiallelic hardcalls are of course especially relevant
// for analysis of major histocompatibility complex genomic data...
uintptr_t GetMhcWordCt(uintptr_t sample_ct);

// sample_ct not relevant if genovecs_ptr == nullptr
// only possible error is kPglRetNomem for now
// caller should reset pgfip->block_base to nullptr when it exits
PglErr PgenMtLoadInit(const uintptr_t* variant_include, uint32_t sample_ct, uint32_t variant_ct, uintptr_t bytes_avail, uintptr_t pgr_alloc_cacheline_ct, uintptr_t thread_xalloc_cacheline_ct, uintptr_t per_variant_xalloc_byte_ct, uintptr_t per_alt_allele_xalloc_byte_ct, PgenFileInfo* pgfip, uint32_t* calc_thread_ct_ptr, uintptr_t*** genovecs_ptr, uintptr_t*** mhc_ptr, uintptr_t*** phasepresent_ptr, uintptr_t*** phaseinfo_ptr, uintptr_t*** dosage_present_ptr, Dosage*** dosage_mains_ptr, uintptr_t*** dphase_present_ptr, SDosage*** dphase_delta_ptr, uint32_t* read_block_size_ptr, uintptr_t* max_alt_allele_block_size_ptr, STD_ARRAY_REF(unsigned char*, 2) main_loadbufs, PgenReader*** pgr_pps, uint32_t** read_variant_uidx_starts_ptr);

// Returns number of variants in current block.  Increases read_block_idx as
// necessary (to get to a nonempty block).
// reterr can be ReadFail but not MalformedInput, since this function just
// reads the raw bytes, it does not perform any unpacking.
uint32_t MultireadNonempty(const uintptr_t* variant_include, const ThreadGroup* tgp, uint32_t raw_variant_ct, uint32_t read_block_size, PgenFileInfo* pgfip, uint32_t* read_block_idxp, PglErr* reterrp);

// Assumes mhc != nullptr, and is vector-aligned.
void ExpandMhc(uint32_t sample_ct, uintptr_t* mhc, uintptr_t** patch_01_set_ptr, AlleleCode** patch_01_vals_ptr, uintptr_t** patch_10_set_ptr, AlleleCode** patch_10_vals_ptr);

HEADER_INLINE void SetPgvMhc(uint32_t sample_ct, uintptr_t* mhc, PgenVariant* pgvp) {
  ExpandMhc(sample_ct, mhc, &(pgvp->patch_01_set), &(pgvp->patch_01_vals), &(pgvp->patch_10_set), &(pgvp->patch_10_vals));
}

HEADER_INLINE void SetPgvThreadMhcNull(uint32_t sample_ct, uint32_t tidx, uintptr_t** thread_mhc, PgenVariant* pgvp) {
  if (!thread_mhc) {
    pgvp->patch_01_set = nullptr;
    pgvp->patch_01_vals = nullptr;
    pgvp->patch_10_set = nullptr;
    pgvp->patch_10_vals = nullptr;
  } else {
    ExpandMhc(sample_ct, thread_mhc[tidx], &(pgvp->patch_01_set), &(pgvp->patch_01_vals), &(pgvp->patch_10_set), &(pgvp->patch_10_vals));
  }
}


HEADER_INLINE BoolErr StoreStringAtBase(unsigned char* arena_top, const char* src, uintptr_t slen, unsigned char** arena_bottom_ptr, char** dst) {
  if (unlikely(slen >= S_CAST(uintptr_t, arena_top - (*arena_bottom_ptr)))) {
    return 1;
  }
  memcpyx(*arena_bottom_ptr, src, slen, '\0');
  *dst = R_CAST(char*, *arena_bottom_ptr);
  *arena_bottom_ptr += slen + 1;
  return 0;
}

HEADER_INLINE BoolErr StoreStringAtEnd(unsigned char* arena_bottom, const char* src, uintptr_t slen, unsigned char** arena_top_ptr, char** dst) {
  // minor todo: verify that the tiny amount of additional safety provided by
  // PtrWSubCk has no/negligible performance cost
  if (PtrWSubCk(arena_bottom, slen + 1, arena_top_ptr)) {
    return 1;
  }
  memcpyx(*arena_top_ptr, src, slen, '\0');
  *dst = R_CAST(char*, *arena_top_ptr);
  return 0;
}

HEADER_INLINE BoolErr StoreStringAtEndK(unsigned char* arena_bottom, const char* src, uintptr_t slen, unsigned char** arena_top_ptr, const char** dst) {
  if (PtrWSubCk(arena_bottom, slen + 1, arena_top_ptr)) {
    return 1;
  }
  memcpyx(*arena_top_ptr, src, slen, '\0');
  *dst = R_CAST(char*, *arena_top_ptr);
  return 0;
}

HEADER_INLINE BoolErr StoreStringAndPrecharAtEnd(unsigned char* arena_bottom, const char* src, unsigned char prechar, uintptr_t slen, unsigned char** arena_top_ptr, char** dst) {
  if (PtrWSubCk(arena_bottom, slen + 2, arena_top_ptr)) {
    return 1;
  }
  **arena_top_ptr = prechar;
  char* dst_write = 1 + R_CAST(char*, *arena_top_ptr);
  memcpyx(dst_write, src, slen, '\0');
  *dst = dst_write;
  return 0;
}

// Memory layout:
//   [       strset       ][packed strings]
//                                         ^
//                                         |
//                                   arena_bottom
// strset pointers do not remain valid after resize.

// returns 1 if resize needed, 2 if OOM
// uint32_t StrsetAdd(unsigned char* arena_top, const char* src, uint32_t slen, uint32_t strset_table_size, char** strset, uint32_t* str_ctp, unsigned char** arena_bottom_ptr);

// uint32_t StrsetAddEnd(unsigned char* arena_bottom, const char* src, uint32_t slen, uint32_t strset_table_size, char** strset, uint32_t* str_ctp, unsigned char** arena_top_ptr);

BoolErr StrsetAddResize(unsigned char* arena_top, const char* src, uint32_t slen, uint32_t strset_table_size_max, char** strset, uint32_t* strset_table_sizep, uint32_t* str_ctp, unsigned char** arena_bottom_ptr);

BoolErr StrsetAddEndResize(unsigned char* arena_bottom, const char* src, uint32_t slen, uint32_t strset_table_size_max, char*** strsetp, uint32_t* strset_table_sizep, uint32_t* str_ctp, unsigned char** arena_top_ptr);

// These use g_textbuf.
PglErr WriteSampleIdsOverride(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* outname, uint32_t sample_ct, SampleIdFlags override_flags);

HEADER_INLINE PglErr WriteSampleIds(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* outname, uint32_t sample_ct) {
  return WriteSampleIdsOverride(sample_include, siip, outname, sample_ct, siip->flags);
}

// read_realpath must be a buffer of size >= kPglFnamesize bytes
uint32_t RealpathIdentical(const char* outname, const char* read_realpath, char* write_realpath_buf);

// assumes rawval is in [1, 32767]
static_assert(kDosageMax == 32768, "PrintHaploidNonintDosage() needs to be updated.");
HEADER_INLINE char* PrintHaploidNonintDosage(uint32_t rawval, char* start) {
  // Instead of constant 5-digit precision, we print fewer digits whenever that
  // doesn't interfere with proper round-tripping.  I.e. we search for the
  // shortest string in
  //   ((n - 0.5)/32768, (n + 0.5)/32768).
  assert(rawval - 1 < 32767);
  *start++ = '0';
  return PrintDdosageDecimal(rawval, start);
}

HEADER_INLINE char* PrintHaploidDosage(uint32_t rawval, char* start) {
  *start++ = '0' + (rawval / kDosageMax);
  const uint32_t remainder = rawval % kDosageMax;
  if (!remainder) {
    return start;
  }
  return PrintDdosageDecimal(remainder, start);
}

char* PrintMultiallelicHcAsDs(uint32_t hc1, uint32_t hc2, uint32_t allele_ct, char* start);

char* PrintMultiallelicHcAsHaploidDs(uint32_t hc1, uint32_t hc2, uint32_t allele_ct, char* start);


HEADER_INLINE void OutnameZstSet(const char* ext, uint32_t output_zst, char* outname_end) {
  const uint32_t ext_slen = strlen(ext);
  assert(ext_slen < kMaxOutfnameExtBlen - 4);
  memcpy(outname_end, ext, ext_slen + 1);
  if (output_zst) {
    strcpy_k(&(outname_end[ext_slen]), ".zst");
  }
}

ENUM_U31_DEF_START()
  kVfilterExtract,
  kVfilterExtractIntersect,
  kVfilterExclude
ENUM_U31_DEF_END(VfilterType);

extern const char g_vft_names[3][18];

HEADER_INLINE BoolErr CleanupPgfi2(const char* file_descrip, PgenFileInfo* pgfip, PglErr* reterrp) {
  if (CleanupPgfi(pgfip, reterrp)) {
    logerrprintfww(kErrprintfFread, file_descrip, strerror(errno));
    return 1;
  }
  return 0;
}

HEADER_INLINE BoolErr CleanupPgr2(const char* file_descrip, PgenReader* pgrp, PglErr* reterrp) {
  if (CleanupPgr(pgrp, reterrp)) {
    logerrprintfww(kErrprintfFread, file_descrip, strerror(errno));
    return 1;
  }
  return 0;
}

void PgenErrPrintEx(const char* file_descrip, uint32_t prepend_lf, PglErr reterr, uint32_t variant_uidx);

HEADER_INLINE void PgenErrPrintN(PglErr reterr) {
  PgenErrPrintEx(".pgen file", 1, reterr, UINT32_MAX);
}

HEADER_INLINE void PgenErrPrintNV(PglErr reterr, uint32_t variant_uidx) {
  PgenErrPrintEx(".pgen file", 1, reterr, variant_uidx);
}

HEADER_INLINE void PgenErrPrint(PglErr reterr) {
  PgenErrPrintEx(".pgen file", 0, reterr, UINT32_MAX);
}

HEADER_INLINE void PgenErrPrintV(PglErr reterr, uint32_t variant_uidx) {
  PgenErrPrintEx(".pgen file", 0, reterr, variant_uidx);
}

// This should be called on SpgwFinish() failure when kPgenWriteAndCopy mode
// may have been used; otherwise it's unnecessary.
void PgenWriteFinishErrPrint(PglErr reterr, char* outname, char* outname_end);

// Given <outname>.tmp.pgen and <outname>.tmp.pgen.pgi, this generates
// <outname>.pgen and then deletes the two temporary files.
PglErr EmbedPgenIndex(char* outname, char* outname_end);

// Before 3 Jan 2023, PLINK 2 usually assumed that in INFO= and FORMAT= VCF
// header lines, the ID key-value pair appeared first, then Number, then Type.
// However, the 2022 Byrska-Bishop et al. 1000 Genomes hg38 callset includes a
// a header line where Type appears before Number, the VCF specification does
// not clearly forbid this, and parsing of these header lines takes negligible
// time in practically all real workflows, so it's time to support arbitrary
// key ordering.
//
// Assumptions:
// - The line is \n-terminated.
// - It's safe to read a few bytes past the end of the line.
CXXCONST_CP HkvlineValEnd(const char* val_start);

#ifdef __cplusplus
HEADER_INLINE char* HkvlineValEnd(char* val_start) {
  return const_cast<char*>(HkvlineValEnd(const_cast<const char*>(val_start)));
}
#endif

// returns 1 if key not found, 2 if line malformed
IntErr HkvlineFind2K(const char* hkvline_iter, const char* key, uint32_t key_slen, const char** val_ptr, uint32_t* val_slenp);

HEADER_INLINE IntErr HkvlineFind2(const char* hkvline_iter, const char* key, uint32_t key_slen, char** val_ptr, uint32_t* val_slenp) {
  return HkvlineFind2K(hkvline_iter, key, key_slen, K_CAST(const char**, val_ptr), val_slenp);
}

HEADER_INLINE IntErr HkvlineFindK(const char* hkvline_iter, const char* key, const char** val_ptr, uint32_t* val_slenp) {
  return HkvlineFind2K(hkvline_iter, key, strlen(key), val_ptr, val_slenp);
}

HEADER_INLINE IntErr HkvlineFind(const char* hkvline_iter, const char* key, char** val_ptr, uint32_t* val_slenp) {
  return HkvlineFind2(hkvline_iter, key, strlen(key), val_ptr, val_slenp);
}


// Input/output:
// - hkvline_iterp, initially pointing to the character after the '<'.  On exit
//   success, points to character after '>' if ID is only key; otherwise
//   beginning of the second key if ID is first key, or stays in place if ID is
//   positioned later.
//   Goal is to avoid reparsing the ID entry when we need to look at later
//   entries.
// Output:
// - idval_ptr and id_slenp, to be set to start of ID value and ID length,
//   respectively.
// Works on contig=/FILTER= header lines, as well as FORMAT=/INFO=.
HEADER_INLINE BoolErr HkvlineIdK(const char** hkvline_iterp, const char** idval_ptr, uint32_t* id_slenp) {
  const char* hkvline_iter = *hkvline_iterp;
  if (HkvlineFindK(hkvline_iter, "ID", idval_ptr, id_slenp)) {
    return 1;
  }
  if (*idval_ptr == &(hkvline_iter[1 + strlen("ID")])) {
    *hkvline_iterp = &((*idval_ptr)[1 + (*id_slenp)]);
  }
  return 0;
}

HEADER_INLINE BoolErr HkvlineId(char** hkvline_iterp, char** idval_ptr, uint32_t* id_slenp) {
  return HkvlineIdK(K_CAST(const char**, hkvline_iterp), K_CAST(const char**, idval_ptr), id_slenp);
}

// Convenience function for scraping Number and Type values simultaneously.
// (HkvlineFind works when only one of the two is needed.)
// Input:
// - hkvline_iter, as set by part 1
// Output:
// - numstr_ptr and num_slenp, to be set to start of Number value string and
//   string length, respectively
// - typestr_ptr and type_slenp, to be set to start of Type value and Type
//   length, respectively
BoolErr HkvlineNumTypeK(const char* hkvline_iter, const char** numstr_ptr, uint32_t* num_slenp, const char** typestr_ptr, uint32_t* type_slenp);

HEADER_INLINE BoolErr HkvlineNumType(char* hkvline_iter, char** numstr_ptr, uint32_t* num_slenp, char** typestr_ptr, uint32_t* type_slenp) {
  return HkvlineNumTypeK(hkvline_iter, K_CAST(const char**, numstr_ptr), num_slenp, K_CAST(const char**, typestr_ptr), type_slenp);
}

// In the unlikely event ID= isn't first, move it in front.
PglErr HkvlineForceIdFirst(uintptr_t workspace_size, char* hline_kv_start, void* workspace);


// 'maybe' bit should be at the bottom of shifted_flag_bits.
HEADER_INLINE uint32_t ProvrefCol(const uintptr_t* variant_include, const uintptr_t* nonref_flags, uint32_t shifted_flag_bits, uint32_t raw_variant_ct, uint32_t all_nonref) {
  if (shifted_flag_bits & 2) {
    return 1;
  }
  if (!(shifted_flag_bits & 1)) {
    return 0;
  }
  if (!nonref_flags) {
    return all_nonref;
  }
  return !IntersectionIsEmpty(variant_include, nonref_flags, BitCtToWordCt(raw_variant_ct));
}

PglErr NsortDedupAndWrite(const char* outname, uintptr_t str_ct, uint32_t max_slen, uint32_t output_zst, uint32_t max_thread_ct, const char** strptr_arr, uintptr_t* nwrite_ptr);

PglErr AllocAndFillVariantLastAlidxs(const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t max_thread_ct, const uintptr_t** variant_last_alidxsp, const uint32_t** variant_last_alidxs_cumulative_popcountsp);


// Assumes IsSet(variant_include, variant_uidx) is true.
// [start_vidx, end_vidx) must be initialized to search bounds.  (Could have
// wrapper that initializes them to variant_uidx's chromosome bounds.)
void GetBpWindow(const uintptr_t* variant_include, const uint32_t* variant_bps, uint32_t variant_uidx, uint32_t bp_radius, uint32_t* start_vidxp, uint32_t* end_vidxp);

PgenGlobalFlags GflagsVfilter(const uintptr_t* variant_include, const unsigned char* vrtypes, uint32_t raw_variant_ct, PgenGlobalFlags input_gflags);

uint32_t TriangleDivide(int64_t cur_prod_x2, int32_t modif);

void ParallelBounds(uint32_t ct, int32_t start, uint32_t parallel_idx, uint32_t parallel_tot, int32_t* __restrict bound_start_ptr, int32_t* __restrict bound_end_ptr);

// Unsafe since this may write up to 6 bytes past the end.
void AppendZerotabsUnsafe(uint32_t zerotab_ct, char** cswritep_ptr);

// Unsafe since this may write up to word boundary past the end.
PglErr InitAllelePermuteUnsafe(const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t max_thread_ct, AlleleCode* allele_permute);

HEADER_INLINE void PglLogputsb() {
  char* pgl_errbuf = PglReturnLog();
  if (pgl_errbuf) {
    WordWrapMultiline(pgl_errbuf);
    logputs(pgl_errbuf);
    PglResetLog();
  }
}

PglErr TrimAllelePermute(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const uintptr_t* allele_presents, uint32_t variant_ct, uint32_t max_thread_ct, AlleleCode* allele_permute);

HEADER_INLINE void UpdateEighash(const char* newstr, uint32_t slen, uint32_t* hash_ptr) {
  uint32_t strhash = 0;
  for (uint32_t uii = 0; uii != slen; ++uii) {
    strhash = (strhash * 23) + ctou32(newstr[uii]);
  }
  *hash_ptr = ((*hash_ptr) * 17) ^ strhash;
}

// '2' since plain LoadTokensNondup() should be based on dupflag hash table.
PglErr LoadTokensNondup2(const char* fname, const uintptr_t* variant_include, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* flagname_p, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t max_thread_ct, uintptr_t* loaded_variant_set);

PglErr LoadTokensNondupReindex(const char* fname, const uintptr_t* variant_include, const uint32_t* variant_include_cumulative_popcounts, const uint32_t* old_variant_uidx_to_new, const char* const* variant_ids, const char* flagname_p, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t max_thread_ct, uintptr_t* loaded_variant_set);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_COMMON_H__
