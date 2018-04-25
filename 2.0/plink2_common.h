#ifndef __PLINK2_COMMON_H__
#define __PLINK2_COMMON_H__

// This library is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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

#include "pgenlib_internal.h"
#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

#define PROG_NAME_STR "plink2"

// leave the door semi-open to 32-bit dosages (or 31?  24?)
// note that u31tod()/u31tof() can't be used on 32-bit dosages
typedef uint16_t Dosage;
typedef int16_t SDosage;
typedef uint32_t DosageProd;
#define kDosageMax (1U << (8 * sizeof(Dosage) - 1))
CONSTU31(kDosageMid, kDosageMax / 2);
CONSTU31(kDosage4th, kDosageMax / 4);
CONSTU31(kDosageMissing, kDosageMax * 2 - 1);
static const double kRecipDosageMax = 0.000030517578125;
static const double kRecipDosageMid = 0.00006103515625;
static const double kRecipDosageMidSq = 0.0000000037252902984619140625;
static const float kRecipDosageMidf = 0.00006103515625;
CONSTU31(kDosagePerVec, kBytesPerVec / sizeof(Dosage));
CONSTU31(kDosagePerCacheline, kCacheline / sizeof(Dosage));

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

// this is a bit arbitrary
CONSTU31(kMaxPhenoCt, 524287);
#define MAX_PHENO_CT_STR "524287"

// maximum number of usable cluster computers, this is arbitrary though it
// shouldn't be larger than 2^32 - 1
// (actually, there's an overflow danger: [work units] * parallel_idx may not
// fit in a uint64 if parallel_tot is too high.)
CONSTU31(kParallelMax, 32768);

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
  kfMiscExcludePvarFilterFail = (1 << 6),
  kfMiscAutosomePar = (1 << 7),
  kfMiscAutosomeOnly = (1 << 8),
  kfMiscMergePar = (1 << 9),
  kfMiscPhenoColNums = (1 << 10),
  kfMiscCovarColNums = (1 << 11),
  kfMiscHweMidp = (1 << 12),
  kfMiscHweKeepFewhet = (1 << 13),
  kfMiscWriteSnplistZs = (1 << 14),
  kfMiscMafSucc = (1 << 15),
  kfMiscGenoDosage = (1 << 16),
  kfMiscGenoHhMissing = (1 << 17),
  kfMiscMindDosage = (1 << 18),
  kfMiscMindHhMissing = (1 << 19),
  kfMiscGenotypingRateDosage = (1 << 20),
  kfMiscSetMissingVarIds = (1 << 21),
  kfMiscChrOverrideCmdline = (1 << 22),
  kfMiscChrOverrideFile = (1 << 23),
  kfMiscNewVarIdOverflowMissing = (1 << 24),
  kfMiscNewVarIdOverflowTruncate = (1 << 25),
  kfMiscRequirePheno = (1 << 26),
  kfMiscRequireCovar = (1 << 27),
  kfMiscCatPhenoFamily = (1 << 28),
  kfMiscRefAlleleForce = (1 << 29),
  kfMiscAlt1AlleleForce = (1 << 30),
  kfMiscRefFromFaForce = (1LLU << 31),
  kfMiscMergeX = (1LLU << 32),
  kfMiscNoIdHeader = (1LLU << 33),
  kfMiscNoIdHeaderIidOnly = (1LLU << 34),
  kfMiscBiallelicOnly = (1LLU << 35),
  kfMiscBiallelicOnlyStrict = (1LLU << 36),
  kfMiscBiallelicOnlyList = (1LLU << 37),
  kfMiscStrictSid0 = (1LLU << 38),
  kfMiscAllowBadFreqs = (1LLU << 39)
FLAGSET64_DEF_END(MiscFlags);

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
  kfExportfTypemask = (2LLU * kfExportfVcf) - kfExportf23,
  kfExportfIncludeAlt = (1LLU << 32),
  kfExportfBgz = (1LLU << 33),
  kfExportfOmitNonmaleY = (1LLU << 34)
FLAGSET64_DEF_END(ExportfFlags);

FLAGSET_DEF_START()
  kfSampleId0,
  kfSampleIdFidPresent = (1 << 0),
  kfSampleIdParentsPresent = (1 << 1),
  kfSampleIdNoIdHeader = (1 << 2),
  kfSampleIdNoIdHeaderIidOnly = (1 << 3),
  kfSampleIdStrictSid0 = (1 << 4)
FLAGSET_DEF_END(SampleIdFlags);

FLAGSET_DEF_START()
  kfInfo0,
  kfInfoPrFlagPresent = (1 << 0),
  kfInfoPrNonflagPresent = (1 << 1),
  kfInfoNonprPresent = (1 << 2),
  kfInfoAll = ((kfInfoNonprPresent * 2) - kfInfoPrFlagPresent),
  // this is set in the .bim case
  kfInfoPrNonrefDefault = (1 << 3),
FLAGSET_DEF_END(InfoFlags);

typedef struct SampleIdInfoStruct {
  char* sample_ids;
  char* sids;
  uintptr_t max_sample_id_blen;
  uintptr_t max_sid_blen;
  SampleIdFlags flags;
} SampleIdInfo;

typedef struct PaternalIdInfoStruct {
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
CONSTU31(kApermMax, 1073241823);

typedef struct TwoColParamsStruct {
  uint32_t colx;
  uint32_t colid;
  uint32_t skip_ct;
  char skipchar;
  char fname[];
} TwoColParams;

void InitPedigreeIdInfo(MiscFlags misc_flags, PedigreeIdInfo* piip);

// no CleanupPedigreeIdInfo function since LoadPsam() allocates the arrays in
// bigstack.


HEADER_INLINE BoolErr bigstack_alloc_dosage(uintptr_t ct, Dosage** dosage_arr_ptr) {
  *dosage_arr_ptr = S_CAST(Dosage*, bigstack_alloc(ct * sizeof(Dosage)));
  return !(*dosage_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_dphase(uintptr_t ct, SDosage** dphase_arr_ptr) {
  *dphase_arr_ptr = S_CAST(SDosage*, bigstack_alloc(ct * sizeof(SDosage)));
  return !(*dphase_arr_ptr);
}


HEADER_INLINE BoolErr bigstack_alloc_dosagep(uintptr_t ct, Dosage*** dosagep_arr_ptr) {
  *dosagep_arr_ptr = S_CAST(Dosage**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*dosagep_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_alloc_dphasep(uintptr_t ct, SDosage*** dphasep_arr_ptr) {
  *dphasep_arr_ptr = S_CAST(SDosage**, bigstack_alloc(ct * sizeof(intptr_t)));
  return !(*dphasep_arr_ptr);
}


HEADER_INLINE BoolErr bigstack_end_alloc_dosage(uintptr_t ct, Dosage** dosage_arr_ptr) {
  *dosage_arr_ptr = S_CAST(Dosage*, bigstack_end_alloc(ct * sizeof(Dosage)));
  return !(*dosage_arr_ptr);
}

HEADER_INLINE BoolErr bigstack_end_alloc_dphase(uintptr_t ct, SDosage** dphase_arr_ptr) {
  *dphase_arr_ptr = S_CAST(SDosage*, bigstack_end_alloc(ct * sizeof(SDosage)));
  return !(*dphase_arr_ptr);
}


char* dosagetoa(uint64_t dosage, char* start);

HEADER_INLINE void ZeroDosageArr(uintptr_t entry_ct, Dosage* dosage_arr) {
  memset(dosage_arr, 0, entry_ct * sizeof(Dosage));
}

HEADER_INLINE void ZeroDphaseArr(uintptr_t entry_ct, SDosage* dphase_arr) {
  memset(dphase_arr, 0, entry_ct * sizeof(SDosage));
}

HEADER_INLINE void SetAllDosageArr(uintptr_t entry_ct, Dosage* dosage_arr) {
  memset(dosage_arr, 255, entry_ct * sizeof(Dosage));
}

void PopulateDenseDosage(const uintptr_t* genovec, const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, Dosage* dense_dosage);

void SetHetMissing(uintptr_t word_ct, uintptr_t* genovec);

void SetHetMissingCleardosage(uintptr_t word_ct, uintptr_t* genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* dosagepresent, Dosage* dosage_main);

void SetHetMissingKeepdosage(uintptr_t word_ct, uintptr_t* genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* dosagepresent, Dosage* dosage_main);

void GenoarrToNonmissing(const uintptr_t* genoarr, uint32_t sample_ctl2, uintptr_t* nonmissing_bitarr);

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

HEADER_INLINE double GetNonmajFreq(const double* cur_allele_freqs, uint32_t cur_allele_ct) {
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

HEADER_INLINE double GetAlleleFreq(const double* cur_allele_freqs, uint32_t allele_idx, uint32_t cur_allele_ct) {
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


// Functions with Xid in their name deal with both FID/IID (with a single tab
// separator) and FID/IID/SID (two tabs) sample IDs.  (Missing FID is
// represented as FID='0'.)

// With no header line, --keep/--remove and the like should interpret a single
// token as IID (treating FID as 0), and two tokens as FID/IID.  --no-id-header
// does not support IID/SID output, so we don't need to worry about supporting
// FidIidSidOrIidSid.
// With a header line, all four {FID present/absent, SID present/absent}
// combinations are allowed.
FLAGSET_DEF_START()
  kfXidMode0,

  kfXidModeFlagOneTokenOk = (1 << 0),
  kfXidModeFlagNeverFid = (1 << 1),
  kfXidModeFlagSid = (1 << 2),

  kfXidModeFidIid = 0,
  kfXidModeFidIidOrIid = kfXidModeFlagOneTokenOk,
  kfXidModeIid = (kfXidModeFlagOneTokenOk | kfXidModeFlagNeverFid),
  kfXidModeFidIidSid = kfXidModeFlagSid,
  kfXidModeIidSid = (kfXidModeFlagNeverFid | kfXidModeFlagSid)
FLAGSET_DEF_END(XidMode);

// Assumes fixed-width.
HEADER_INLINE uint32_t GetXidColCt(XidMode xid_mode) {
  return 2 + (xid_mode == kfXidModeFidIidSid) - (xid_mode == kfXidModeIid);
}

// sample_xid_map allocated on bottom, to play well with --indiv-sort
PglErr SortedXidboxInitAlloc(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t allow_dups, XidMode xid_mode, uint32_t use_nsort, char** sorted_xidbox_ptr, uint32_t** xid_map_ptr, uintptr_t* max_xid_blen_ptr);

// returns slen for ID, or 0 on parse failure.
uint32_t XidRead(uintptr_t max_xid_blen, uint32_t comma_delim, XidMode xid_mode, const char** read_pp, char* __restrict idbuf);

// returns 1 on missing token *or* if the sample ID is not present.  cases can
// be distinguished by checking whether *read_pp_new == nullptr: if it is, a
// missing-tokens error should probably be reported.
// sample_id_map == nullptr is permitted
// *read_pp is now set to point to the end of IID/SID instead of the beginning
// of the next token; this is a change from plink 1.9.
HEADER_INLINE BoolErr SortedXidboxReadFind(const char* __restrict sorted_xidbox, const uint32_t* __restrict xid_map, uintptr_t max_xid_blen, uintptr_t xid_ct, uint32_t comma_delim, XidMode xid_mode, const char** read_pp, uint32_t* sample_uidx_ptr, char* __restrict idbuf) {
  const uint32_t slen_final = XidRead(max_xid_blen, comma_delim, xid_mode, read_pp, idbuf);
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

// May return kPglRetEof.
// loadbuf_iter_ptr can be nullptr.
// line_idx must be zero unless initial lines were skipped.
// If no header line is present, xid_mode will be set to kfXidModeFidIid if
// kfXidHeaderFixedWidth is set, and kfXidModeFidIidOrIid (which tolerates a
// mix of single-token and multitoken lines) otherwise.
PglErr LoadXidHeader(const char* flag_name, XidHeaderFlags xid_header_flags, char** line_iterp, uintptr_t* line_idx_ptr, char** linebuf_first_token_ptr, ReadLineStream* rlsp, XidMode* xid_mode_ptr);

PglErr OpenAndLoadXidHeader(const char* fname, const char* flag_name, XidHeaderFlags xid_header_flags, uintptr_t linebuf_size, char** line_iterp, uintptr_t* line_idx_ptr, char** linebuf_first_token_ptr, ReadLineStream* rlsp, XidMode* xid_mode_ptr);


// note that this is no longer divisible by 64
CONSTU31(kMaxContigs, 65274);
CONSTU31(kMaxChrCodeDigits, 5);

// change ChrIdx to uint32_t if (kMaxContigs + kChrOffsetCt) > 65536
typedef uint16_t ChrIdx;

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
ENUM_U31_DEF_END(XymtOffset);

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
ENUM_U31_DEF_END(ChrsetSource);

FLAGSET_DEF_START()
  kfChrOutput0,
  kfChrOutputPrefix = (1 << 0),
  kfChrOutputM = (1 << 1),
  kfChrOutputMT = (1 << 2),
  kfChrOutput0M = (1 << 3)
FLAGSET_DEF_END(ChrOutput);

typedef struct {
  // Main dynamic block intended to be allocated as a single aligned block of
  // memory on the heap freeable with vecaligned_free(), with chr_mask at the
  // base.

  uintptr_t* chr_mask;  // which chromosomes aren't known to be absent?

  // This includes chrX.  As of alpha 2, it also includes MT again (like plink
  // 1.07, and unlike 1.9 and 2.0a1), now that enough dosage functionality is
  // in place.
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
  uint32_t xymt_codes[kChrOffsetCt];

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

extern const char g_xymt_log_names[][5];

PglErr InitChrInfo(ChrInfo* cip);

void FinalizeChrset(MiscFlags misc_flags, ChrInfo* cip);

HEADER_INLINE PglErr InitChrInfoHuman(ChrInfo* cip) {
  // convenience wrapper
  if (InitChrInfo(cip)) {
    return kPglRetNomem;
  }
  FinalizeChrset(kfMisc0, cip);
  return kPglRetSuccess;
}

void ForgetExtraChrNames(uint32_t reinitialize, ChrInfo* cip);

// in the usual case where the number of chromosomes/contigs is much less than
// kMaxContigs, this reduces chr_info's memory consumption and improves
// locality.
PglErr FinalizeChrInfo(ChrInfo* cip);

void CleanupChrInfo(ChrInfo* cip);

// assumes chr_idx is valid
char* chrtoa(const ChrInfo* cip, uint32_t chr_idx, char* buf);

uint32_t GetMaxChrSlen(const ChrInfo* cip);

uint32_t IsHaploidChrPresent(const ChrInfo* cip);

// any character <= ' ' is considered a terminator
// maps chrX -> kChrRawX, etc.
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

HEADER_INLINE uint32_t GetVariantChrFoIdx(const ChrInfo* cip, uintptr_t variant_uidx) {
  return CountSortedSmallerU32(&(cip->chr_fo_vidx_start[1]), cip->chr_ct, variant_uidx + 1);
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
  // assumes xymt_exists was previously called, and is true
  const uint32_t xymt_code = cip->xymt_codes[xymt_offset];
  *xymt_code_ptr = xymt_code;
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[xymt_code];
  *xymt_start_ptr = cip->chr_fo_vidx_start[chr_fo_idx];
  *xymt_end_ptr = cip->chr_fo_vidx_start[chr_fo_idx + 1];
}

// now assumes chr_name is null-terminated
PglErr TryToAddChrName(const char* chr_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t allow_extra_chrs, uint32_t* chr_idx_ptr, ChrInfo* cip);

HEADER_INLINE PglErr GetOrAddChrCode(const char* chr_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t allow_extra_chrs, ChrInfo* cip, uint32_t* chr_idx_ptr) {
  *chr_idx_ptr = GetChrCode(chr_name, cip, name_slen);
  if (!IsI32Neg(*chr_idx_ptr)) {
    return kPglRetSuccess;
  }
  return TryToAddChrName(chr_name, file_descrip, line_idx, name_slen, allow_extra_chrs, chr_idx_ptr, cip);
}

HEADER_INLINE PglErr GetOrAddChrCodeDestructive(const char* file_descrip, uintptr_t line_idx, uint32_t allow_extra_chrs, char* chr_name, char* chr_name_end, ChrInfo* cip, uint32_t* chr_idx_ptr) {
  *chr_name_end = '\0';
  return GetOrAddChrCode(chr_name, file_descrip, line_idx, chr_name_end - chr_name, allow_extra_chrs, cip, chr_idx_ptr);
}

// zeroes out samples not in the mask
void InterleavedMaskZero(const uintptr_t* __restrict interleaved_mask, uintptr_t vec_ct, uintptr_t* __restrict genovec);

// sets samples in the mask to missing (0b11)
void InterleavedSetMissing(const uintptr_t* __restrict interleaved_set, uintptr_t vec_ct, uintptr_t* __restrict genovec);

void InterleavedSetMissingCleardosage(const uintptr_t* __restrict orig_set, const uintptr_t* __restrict interleaved_set, uintptr_t vec_ct, uintptr_t* __restrict genovec, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main);

void SetMaleHetMissing(const uintptr_t* __restrict sex_male_interleaved, uint32_t vec_ct, uintptr_t* __restrict genovec);

void EraseMaleDphases(const uintptr_t* __restrict sex_male, uint32_t* __restrict write_dphase_ct_ptr, uintptr_t* __restrict dphasepresent, SDosage* dphase_delta);

void SetMaleHetMissingCleardosage(const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_male_interleaved, uint32_t vec_ct, uintptr_t* __restrict genovec, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main);

void SetMaleHetMissingKeepdosage(const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_male_interleaved, uint32_t word_ct, uintptr_t* __restrict genovec, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main);

// Clears each bit in bitarr which doesn't correspond to a genovec het.
// Assumes that either trailing bits of bitarr are already zero, or trailing
// bits of genovec are zero.
void MaskGenovecHetsUnsafe(const uintptr_t* __restrict genovec, uint32_t raw_sample_ctl2, uintptr_t* __restrict bitarr);

// vertical popcount support
// VcountScramble1() and friends in plink2_cmdline
/*
#ifdef __LP64__
#  ifdef USE_AVX2
// 2->4: 0 2 ... 126 1 3 ... 127
// 4->8: 0 4 ... 124 2 6 ... 126 1 5 ... 125 3 7 ... 127
// 8->32: 0 16 ... 112 4 20 ... 116 ... 124 2 18 ... 114 6 22 ... 118 ... 126 1 17 ...
HEADER_CINLINE uint32_t VcountScramble2(uint32_t orig_idx) {
  return (orig_idx & (~127)) + ((orig_idx & 1) * 64) + ((orig_idx & 2) * 16) + ((orig_idx & 12) * 2) + ((orig_idx & 112) / 16);
}
#  else
HEADER_CINLINE uint32_t VcountScramble2(uint32_t orig_idx) {
  return (orig_idx & (~63)) + ((orig_idx & 1) * 32) + ((orig_idx & 2) * 8) + (orig_idx & 12) + ((orig_idx & 48) / 16);
}
#  endif
#else
// 2->4: 0 2 4 6 8 10 12 14 1 3 5 ...
// 4->8: 0 4 8 12 2 6 10 14 1 5 9 ...
// 8->32: 0 4 8 12 2 6 10 14 1 5 9 13 3 7 11 15
HEADER_CINLINE uint32_t VcountScramble2(uint32_t orig_idx) {
  return (orig_idx & (~15)) + ((orig_idx & 1) * 8) + ((orig_idx & 2) * 2) + ((orig_idx & 12) / 4);
}
#endif

// probable todo: switch to VecW* parameters
HEADER_INLINE void VcountIncr2To4(const uintptr_t* acc2, uint32_t acc2_vec_ct, uintptr_t* acc4) {
  const VecW m2 = VCONST_W(kMask3333);
  const VecW* acc2v_iter = R_CAST(const VecW*, acc2);
  VecW* acc4v_iter = R_CAST(VecW*, acc4);
  for (uint32_t vidx = 0; vidx < acc2_vec_ct; ++vidx) {
    VecW loader = *acc2v_iter++;
    *acc4v_iter = (*acc4v_iter) + (loader & m2);
    ++acc4v_iter;
    loader = vecw_srli(loader, 2);
    *acc4v_iter = (*acc4v_iter) + (loader & m2);
    ++acc4v_iter;
  }
}

HEADER_INLINE void Vcount0Incr2To4(uint32_t acc2_vec_ct, uintptr_t* acc2, uintptr_t* acc4) {
  const VecW m2 = VCONST_W(kMask3333);
  VecW* acc2v_iter = R_CAST(VecW*, acc2);
  VecW* acc4v_iter = R_CAST(VecW*, acc4);
  for (uint32_t vidx = 0; vidx < acc2_vec_ct; ++vidx) {
    VecW loader = *acc2v_iter;
    *acc2v_iter++ = vecw_setzero();
    *acc4v_iter = (*acc4v_iter) + (loader & m2);
    ++acc4v_iter;
    loader = vecw_srli(loader, 2);
    *acc4v_iter = (*acc4v_iter) + (loader & m2);
    ++acc4v_iter;
  }
}
*/


// uint32_t chr_window_max(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_pos, uint32_t chr_fo_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max);

HEADER_INLINE uint32_t CountChrVariantsUnsafe(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t chr_idx) {
  assert(IsSet(cip->chr_mask, chr_idx));
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
  const uint32_t min_idx = cip->chr_fo_vidx_start[chr_fo_idx];
  const uint32_t max_idx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  return PopcountBitRange(variant_include, min_idx, max_idx);
}

HEADER_INLINE uint32_t ChrIsNonempty(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t chr_idx) {
  if (!IsSet(cip->chr_mask, chr_idx)) {
    return 0;
  }
  const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
  const uint32_t min_idx = cip->chr_fo_vidx_start[chr_fo_idx];
  const uint32_t max_idx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  return !AllBitsAreZero(variant_include, min_idx, max_idx);
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

PglErr ConditionalAllocateNonAutosomalVariants(const ChrInfo* cip, const char* calc_descrip, uint32_t raw_variant_ct, const uintptr_t** variant_include_ptr, uint32_t* variant_ct_ptr);

void FillSubsetChrFoVidxStart(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t* subset_chr_fo_vidx_start);

HEADER_INLINE BoolErr AllocAndFillSubsetChrFoVidxStart(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t** subset_chr_fo_vidx_start_ptr) {
  const uint32_t chr_ct = cip->chr_ct;
  if (bigstack_alloc_u32(chr_ct + 1, subset_chr_fo_vidx_start_ptr)) {
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

CONSTU31(kMaxMissingPhenostrBlen, 32);
// might want g_input_missing_catname and/or g_output_missing_catname later,
// but let's start with the simplest implementation
extern char g_missing_catname[];  // default "NONE", not changeable for now

extern char g_output_missing_pheno[];  // default "NA"
extern char g_legacy_output_missing_pheno[];  // default "-9"

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

typedef struct {
  // * If categorical phenotype, [0] points to g_missing_catname, while [1],
  //   [2], etc. point to category names.  These are part of the same
  //   allocation as nonmiss, so no separate free is needed.
  //   Otherwise, this is nullptr.
  // * When .sample categorical variables are imported, 'P' is added in front
  //   of the integers.
  const char** category_names;

  uintptr_t* nonmiss;  // bitvector

  // essentially a tagged union; part of the same allocation as nonmiss
  PhenoData data;
  PhenoDtype type_code;

  uint32_t nonnull_category_ct;
} PhenoCol;

void InitPheno();


uint32_t IsCategoricalPhenostr(const char* phenostr_iter);

uint32_t IsCategoricalPhenostrNocsv(const char* phenostr_iter);

// returns 0xffffffffU if none exists
uint32_t FirstCcOrQtPhenoIdx(const PhenoCol* pheno_cols, uint32_t pheno_ct);

// "_covar" since this doesn't handle case/control
uint32_t IsConstCovar(const PhenoCol* covar_col, const uintptr_t* sample_include, uint32_t sample_ct);

uint32_t IdentifyRemainingCats(const uintptr_t* sample_include, const PhenoCol* covar_col, uint32_t sample_ct, uintptr_t* observed_cat_bitarr);

uint32_t GetCatSamples(const uintptr_t* sample_include_base, const PhenoCol* cat_pheno_col, uint32_t raw_sample_ctl, uint32_t sample_ct, uint32_t cat_uidx, uintptr_t* cur_cat_samples);

// pheno_names is also allocated on the heap, but it can be handled with a
// simple free_cond().
void CleanupPhenoCols(uint32_t pheno_ct, PhenoCol* pheno_cols);

PglErr ParseChrRanges(const char* const* argvk, const char* flagname_p, const char* errstr_append, uint32_t param_ct, uint32_t allow_extra_chrs, uint32_t xymt_subtract, char range_delim, ChrInfo* cip, uintptr_t* chr_mask);


// sample_ct not relevant if genovecs_ptr == nullptr
PglErr PgenMtLoadInit(const uintptr_t* variant_include, uint32_t sample_ct, uint32_t variant_ct, uintptr_t bytes_avail, uintptr_t pgr_alloc_cacheline_ct, uintptr_t thread_xalloc_cacheline_ct, uintptr_t per_variant_xalloc_byte_ct, PgenFileInfo* pgfip, uint32_t* calc_thread_ct_ptr, uintptr_t*** genovecs_ptr, uintptr_t*** phasepresent_ptr, uintptr_t*** phaseinfo_ptr, uintptr_t*** dosage_present_ptr, Dosage*** dosage_mains_ptr, uintptr_t*** dphase_present_ptr, SDosage*** dphase_delta_ptr, uint32_t* read_block_size_ptr, unsigned char** main_loadbufs, pthread_t** threads_ptr, PgenReader*** pgr_pps, uint32_t** read_variant_uidx_starts_ptr);


// These use g_textbuf.
PglErr WriteSampleIdsOverride(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* outname, uint32_t sample_ct, SampleIdFlags override_flags);
HEADER_INLINE PglErr WriteSampleIds(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* outname, uint32_t sample_ct) {
  return WriteSampleIdsOverride(sample_include, siip, outname, sample_ct, siip->flags);
}

uint32_t RealpathIdentical(const char* outname, const char* read_realpath, char* write_realpath_buf);


HEADER_INLINE void OutnameZstSet(const char* ext, uint32_t output_zst, char* outname_end) {
  const uint32_t ext_slen = strlen(ext);
  assert(ext_slen < kMaxOutfnameExtBlen - 4);
  memcpy(outname_end, ext, ext_slen + 1);
  if (output_zst) {
    memcpy(&(outname_end[ext_slen]), ".zst", 5);
  }
}

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_COMMON_H__
