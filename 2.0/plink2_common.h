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
typedef uint16_t dosage_t;
typedef uint32_t dosage_prod_t;
#define kDosageMax (1U << (8 * sizeof(dosage_t) - 1))
CONSTU31(kDosageMid, kDosageMax / 2);
CONSTU31(kDosage4th, kDosageMax / 4);
CONSTU31(kDosageMissing, kDosageMax * 2 - 1);
static const double kRecipDosageMax = 0.000030517578125;
static const double kRecipDosageMid = 0.00006103515625;
static const double kRecipDosageMidSq = 0.0000000037252902984619140625;
static const float kRecipDosageMidf = 0.00006103515625;
CONSTU31(kDosagePerVec, kBytesPerVec / sizeof(dosage_t));
CONSTU31(kDosagePerCacheline, kCacheline / sizeof(dosage_t));

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
  kfMiscKeepfileSid = (1 << 6),
  kfMiscRemovefileSid = (1 << 7),
  kfMiscKeepAutoconv = (1 << 8),
  kfMiscDoubleId = (1 << 9),
  kfMiscBiallelicOnly = (1 << 10),
  kfMiscBiallelicOnlyStrict = (1 << 11),
  kfMiscBiallelicOnlyList = (1 << 12),
  kfMiscExcludePvarFilterFail = (1 << 13),
  kfMiscVcfRequireGt = (1 << 14),
  kfMiscAutosomePar = (1 << 15),
  kfMiscAutosomeOnly = (1 << 16),
  kfMiscMergePar = (1 << 17),
  kfMiscAllowNoSamples = (1 << 18),
  kfMiscAllowNoVars = (1 << 19),
  kfMiscHweMidp = (1 << 20),
  kfMiscHweKeepFewhet = (1 << 21),
  kfMiscWriteSnplistZs = (1 << 22),
  kfMiscMafSucc = (1 << 23),
  kfMiscGenoDosage = (1 << 24),
  kfMiscGenoHhMissing = (1 << 25),
  kfMiscMindDosage = (1 << 26),
  kfMiscMindHhMissing = (1 << 27),
  kfMiscGenotypingRateDosage = (1 << 28),
  kfMiscSetMissingVarIds = (1 << 29),
  kfMiscChrOverrideCmdline = (1 << 30),
  kfMiscChrOverrideFile = (1LLU << 31),
  kfMiscNewVarIdOverflowMissing = (1LLU << 32),
  kfMiscNewVarIdOverflowTruncate = (1LLU << 33),
  kfMiscRequirePheno = (1LLU << 34),
  kfMiscRequireCovar = (1LLU << 35),
  kfMiscCatPhenoFamily = (1LLU << 36),
  kfMiscRefAlleleForce = (1LLU << 37),
  kfMiscAlt1AlleleForce = (1LLU << 38),
  kfMiscRefFromFaForce = (1LLU << 39),
  kfMiscMergeX = (1LLU << 40),
  kfMiscKeepFileStrsSid = (1LLU << 41),
  kfMiscPhenoColNums = (1LLU << 42),
  kfMiscCovarColNums = (1LLU << 43)
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
  kfExportfVcfDosageForce = (1LLU << 36),
  kfExportfOmitNonmaleY = (1LLU << 37)
FLAGSET64_DEF_END(exportf_flags_t);

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


typedef struct two_col_params_struct {
  uint32_t colx;
  uint32_t colid;
  uint32_t skip_ct;
  char skipchar;
  char fname[];
} two_col_params_t;


HEADER_INLINE boolerr_t bigstack_alloc_dosage(uintptr_t ct, dosage_t** dosage_arr_ptr) {
  *dosage_arr_ptr = (dosage_t*)bigstack_alloc(ct * sizeof(dosage_t));
  return !(*dosage_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_alloc_dosagep(uintptr_t ct, dosage_t*** dosagep_arr_ptr) {
  *dosagep_arr_ptr = (dosage_t**)bigstack_alloc(ct * sizeof(intptr_t));
  return !(*dosagep_arr_ptr);
}

HEADER_INLINE boolerr_t bigstack_end_alloc_dosage(uintptr_t ct, dosage_t** dosage_arr_ptr) {
  *dosage_arr_ptr = (dosage_t*)bigstack_end_alloc(ct * sizeof(dosage_t));
  return !(*dosage_arr_ptr);
}


char* print_dosage(uint64_t dosage, char* start);

HEADER_INLINE void fill_dosage_zero(uintptr_t entry_ct, dosage_t* dosage_arr) {
  memset(dosage_arr, 0, entry_ct * sizeof(dosage_t));
}

HEADER_INLINE void fill_dosage_one(uintptr_t entry_ct, dosage_t* dosage_arr) {
  memset(dosage_arr, 255, entry_ct * sizeof(dosage_t));
}

void populate_dense_dosage(const uintptr_t* genovec, const uintptr_t* dosage_present, const dosage_t* dosage_vals, uint32_t sample_ct, uint32_t dosage_ct, dosage_t* dense_dosage);

void set_het_missing(uintptr_t word_ct, uintptr_t* genovec);

void set_het_missing_cleardosage(uintptr_t word_ct, uintptr_t* genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* dosagepresent, dosage_t* dosage_vals);

void set_het_missing_keepdosage(uintptr_t word_ct, uintptr_t* genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* dosagepresent, dosage_t* dosage_vals);

void genoarr_to_nonmissing(const uintptr_t* genoarr, uint32_t sample_ctl2, uintptr_t* nonmissing_bitarr);

uint32_t genoarr_count_missing_notsubset_unsafe(const uintptr_t* genoarr, const uintptr_t* exclude_mask, uint32_t sample_ct);


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

// Assumes fixed-width.
HEADER_INLINE uint32_t get_xid_col_ct(xid_mode_t xid_mode) {
  if (xid_mode == kfXidModeIid) {
    return 1;
  } else {
    return 2 + (xid_mode == kfXidModeFidiidSid);
  }
}

// sample_xid_map allocated on bottom, to play well with --indiv-sort
pglerr_t sorted_xidbox_init_alloc(const uintptr_t* sample_include, const char* sample_ids, const char* sids, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t allow_dups, xid_mode_t xid_mode, uint32_t use_nsort, char** sorted_xidbox_ptr, uint32_t** xid_map_ptr, uintptr_t* max_xid_blen_ptr);

// returns slen for ID, or 0 on parse failure.
uint32_t xid_read(uintptr_t max_xid_blen, uint32_t comma_delim, xid_mode_t xid_mode, const char** read_pp, char* __restrict idbuf);

// returns 1 on missing token *or* if the sample ID is not present.  cases can
// be distinguished by checking whether *read_pp_new == nullptr: if it is, a
// missing-tokens error should probably be reported.
// sample_id_map == nullptr is permitted
// *read_pp is now set to point to the end of IID/SID instead of the beginning
// of the next token; this is a change from plink 1.9.
HEADER_INLINE boolerr_t sorted_xidbox_read_find(const char* __restrict sorted_xidbox, const uint32_t* __restrict xid_map, uintptr_t max_xid_blen, uintptr_t xid_ct, uint32_t comma_delim, xid_mode_t xid_mode, const char** read_pp, uint32_t* sample_uidx_ptr, char* __restrict idbuf) {
  const uint32_t slen_final = xid_read(max_xid_blen, comma_delim, xid_mode, read_pp, idbuf);
  if (!slen_final) {
    return 1;
  }
  return sorted_idbox_find(idbuf, sorted_xidbox, xid_map, slen_final, max_xid_blen, xid_ct, sample_uidx_ptr);
}

// Matches a sample ID *prefix*.  Thus, if FID/IID/SID is loaded, but the input
// file contains just FID/IID, and there are some FID/IID pairs which
// correspond to multiple samples, this lets you iterate over all of them.
// (Caller is responsible for looking up xid_map[] to perform xid_idx ->
// sample_uidx conversions.)
boolerr_t sorted_xidbox_read_multifind(const char* __restrict sorted_xidbox, uintptr_t max_xid_blen, uintptr_t xid_ct, uint32_t comma_delim, xid_mode_t xid_mode, const char** read_pp, uint32_t* __restrict xid_idx_start_ptr, uint32_t* __restrict xid_idx_end_ptr, char* __restrict idbuf);

ENUM_U31_DEF_START()
  kSidDetectModeNotLoaded,
  kSidDetectModeLoaded,
  kSidDetectModeForce
ENUM_U31_DEF_END(sid_detect_mode_t);

// may return kPglRetLongLine or kPglRetEmptyFile
// loadbuf_iter_ptr can be nullptr
// line_idx must be zero unless initial lines were skipped
// Follow this up with
//   if (xid_mode == kfXidModeFidiidOrIid) {
//     xid_mode = kfXidModeFidiid;
//   }
// when using this on a regular .tsv-like file.  (Otherwise, xid_read() will
// tolerate a mix of single-token and multitoken lines, where --double-id
// interpretation is applied to single-token lines.)
pglerr_t load_xid_header(const char* flag_name, sid_detect_mode_t sid_detect_mode, uintptr_t loadbuf_size, char* loadbuf, char** loadbuf_iter_ptr, uintptr_t* line_idx_ptr, char** loadbuf_first_token_ptr, gzFile* gz_infile_ptr, xid_mode_t* xid_mode_ptr);

// sets last character of loadbuf to ' '
pglerr_t open_and_load_xid_header(const char* fname, const char* flag_name, sid_detect_mode_t sid_detect_mode, uintptr_t loadbuf_size, char* loadbuf, char** loadbuf_iter_ptr, uintptr_t* line_idx_ptr, char** loadbuf_first_token_ptr, gzFile* gz_infile_ptr, xid_mode_t* xid_mode_ptr);


// note that this is no longer divisible by 64
CONSTU31(kMaxContigs, 65274);
CONSTU31(kMaxChrCodeDigits, 5);

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
  const char** nonstd_names;
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

// zeroes out samples not in the mask
void interleaved_mask_zero(const uintptr_t* __restrict interleaved_mask, uintptr_t vec_ct, uintptr_t* __restrict genovec);

// sets samples in the mask to missing (0b11)
void interleaved_set_missing(const uintptr_t* __restrict interleaved_set, uintptr_t vec_ct, uintptr_t* __restrict genovec);

void interleaved_set_missing_cleardosage(const uintptr_t* __restrict orig_set, const uintptr_t* __restrict interleaved_set, uintptr_t vec_ct, uintptr_t* __restrict genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, dosage_t* dosage_vals);

void set_male_het_missing(const uintptr_t* __restrict sex_male_interleaved, uint32_t vec_ct, uintptr_t* __restrict genovec);

void set_male_het_missing_cleardosage(const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_male_interleaved, uint32_t vec_ct, uintptr_t* __restrict genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, dosage_t* dosage_vals);

void set_male_het_missing_keepdosage(const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_male_interleaved, uint32_t word_ct, uintptr_t* __restrict genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, dosage_t* dosage_vals);

// Clears each bit in bitarr which doesn't correspond to a genovec het.
// Assumes that either trailing bits of bitarr are already zero, or trailing
// bits of genovec are zero.
void mask_genovec_hets_unsafe(const uintptr_t* __restrict genovec, uint32_t raw_sample_ctl2, uintptr_t* __restrict bitarr);

// vertical popcount support
// scramble_1_4_8_32() and friends in plink2_cmdline
#ifdef __LP64__
  #ifdef USE_AVX2
// 2->4: 0 2 ... 126 1 3 ... 127
// 4->8: 0 4 ... 124 2 6 ... 126 1 5 ... 125 3 7 ... 127
// 8->32: 0 16 ... 112 4 20 ... 116 ... 124 2 18 ... 114 6 22 ... 118 ... 126 1 17 ...
HEADER_CINLINE uint32_t scramble_2_4_8_32(uint32_t orig_idx) {
  return (orig_idx & (~127)) + ((orig_idx & 1) * 64) + ((orig_idx & 2) * 16) + ((orig_idx & 12) * 2) + ((orig_idx & 112) / 16);
}
  #else
HEADER_CINLINE uint32_t scramble_2_4_8_32(uint32_t orig_idx) {
  return (orig_idx & (~63)) + ((orig_idx & 1) * 32) + ((orig_idx & 2) * 8) + (orig_idx & 12) + ((orig_idx & 48) / 16);
}
  #endif
#else
// 2->4: 0 2 4 6 8 10 12 14 1 3 5 ...
// 4->8: 0 4 8 12 2 6 10 14 1 5 9 ...
// 8->32: 0 4 8 12 2 6 10 14 1 5 9 13 3 7 11 15
HEADER_CINLINE uint32_t scramble_2_4_8_32(uint32_t orig_idx) {
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


// uint32_t chr_window_max(const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_pos, uint32_t chr_fo_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max);

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

/*
// newval does not need to be null-terminated
// assumes *allele_ptr is not initialized
// (stop using these in main plink2 binary?)
boolerr_t allele_set(const char* newval, uint32_t allele_slen, char** allele_ptr);

// *allele_ptr must be initialized; frees *allele_ptr if necessary
boolerr_t allele_reset(const char* newval, uint32_t allele_slen, char** allele_ptr);

void cleanup_allele_storage(uint32_t max_allele_slen, uintptr_t allele_storage_entry_ct, const char** allele_storage);
*/

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
  const char** category_names;

  uintptr_t* nonmiss; // bitvector

  // essentially a tagged union; part of the same allocation as nonmiss
  phenodata_t data;
  pheno_dtype_t type_code;

  uint32_t nonnull_category_ct;
} pheno_col_t;

void init_pheno();


uint32_t is_categorical_phenostr(const char* phenostr);

uint32_t is_categorical_phenostr_nocsv(const char* phenostr);

// returns 0xffffffffU if none exists
uint32_t first_cc_or_qt_pheno_idx(const pheno_col_t* pheno_cols, uint32_t pheno_ct);

// "_covar" since this doesn't handle case/control
uint32_t is_const_covar(const pheno_col_t* covar_col, const uintptr_t* sample_include, uint32_t sample_ct);

uint32_t identify_remaining_cats(const uintptr_t* sample_include, const pheno_col_t* covar_col, uint32_t sample_ct, uintptr_t* observed_cat_bitarr);

uint32_t get_is_cat_include(const uintptr_t* sample_include_base, const pheno_col_t* cat_pheno_col, uint32_t raw_sample_ctl, uint32_t sample_ct, uint32_t cat_uidx, uintptr_t* is_cat_include);

// pheno_names is also allocated on the heap, but it can be handled with a
// simple free_cond().
void cleanup_pheno_cols(uint32_t pheno_ct, pheno_col_t* pheno_cols);

pglerr_t parse_chr_ranges(const char* const* argvc, const char* flagname_p, const char* errstr_append, uint32_t param_ct, uint32_t allow_extra_chrs, uint32_t xymt_subtract, char range_delim, chr_info_t* cip, uintptr_t* chr_mask);


// sample_ct not relevant if genovecs_ptr == nullptr
pglerr_t multithread_load_init(const uintptr_t* variant_include, uint32_t sample_ct, uint32_t variant_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t thread_xalloc_cacheline_ct, uintptr_t per_variant_xalloc_byte_ct, pgen_file_info_t* pgfip, uint32_t* calc_thread_ct_ptr, uintptr_t*** genovecs_ptr, uintptr_t*** dosage_present_ptr, dosage_t*** dosage_val_bufs_ptr, uint32_t* read_block_size_ptr, unsigned char** main_loadbufs, pthread_t** threads_ptr, pgen_reader_t*** pgr_pps, uint32_t** read_variant_uidx_starts_ptr);

pglerr_t write_sample_ids(const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* outname, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen);

uint32_t realpath_identical(const char* outname, const char* read_realpath, char* write_realpath_buf);

#ifdef __cplusplus
} // namespace plink2
#endif

#endif // __PLINK2_COMMON_H__
