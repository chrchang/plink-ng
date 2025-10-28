// This library is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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

#include "plink2_common.h"

#include <math.h>
#include <stdlib.h>
#include <unistd.h>

#include "include/plink2_htable.h"
#include "plink2_compress_stream.h"
#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

const double kSmallDoubles[4] = {0.0, 1.0, 2.0, 3.0};

void InitPedigreeIdInfo(MiscFlags misc_flags, PedigreeIdInfo* piip) {
  piip->sii.sample_ids = nullptr;
  piip->sii.sids = nullptr;
  piip->sii.max_sample_id_blen = 4;
  piip->sii.max_sid_blen = 0;
  piip->sii.flags = kfSampleId0;
  if (misc_flags & kfMiscNoIdHeader) {
    piip->sii.flags |= kfSampleIdNoIdHeader;
    if (misc_flags & kfMiscNoIdHeaderIidOnly) {
      piip->sii.flags |= kfSampleIdNoIdHeaderIidOnly;
    }
  }
  if (misc_flags & kfMiscStrictSid0) {
    piip->sii.flags |= kfSampleIdStrictSid0;
  }
  piip->parental_id_info.paternal_ids = nullptr;
  piip->parental_id_info.maternal_ids = nullptr;
  piip->parental_id_info.max_paternal_id_blen = 2;
  piip->parental_id_info.max_maternal_id_blen = 2;
}

void InitFlip(FlipInfo* flip_info_ptr) {
  flip_info_ptr->fname = nullptr;
  flip_info_ptr->subset_fname = nullptr;
  flip_info_ptr->flags = kfFlip0;
}

void CleanupFlip(FlipInfo* flip_info_ptr) {
  free_cond(flip_info_ptr->fname);
  free_cond(flip_info_ptr->subset_fname);
}

const char kLoadFilterLogFlagnames[18][25] = {
  "import-max-alleles",
  "vcf-require-gt",
  "merge-max-alleles",
  "autosome",
  "autosome-par",
  "chr",
  "not-chr",
  "",
  "exclude-if-info",
  "exclude-palindromic-snps",
  "extract-if-info",
  "max-alleles",
  "min-alleles",
  "require-info",
  "require-no-info",
  "snps-only",
  "var-filter",
  "var-min-qual"
};

void AppendLoadFilterFlagnames(LoadFilterLogFlags load_filter_log_flags, char** write_iter_ptr) {
  assert(load_filter_log_flags);
  char* write_iter = *write_iter_ptr;
  while (1) {
    write_iter = strcpya_k(write_iter, "--");
    const uint32_t flag_idx = ctzu32(load_filter_log_flags);
    write_iter = strcpya(write_iter, kLoadFilterLogFlagnames[flag_idx]);
    load_filter_log_flags = S_CAST(LoadFilterLogFlags, load_filter_log_flags & (load_filter_log_flags - 1));
    if (!load_filter_log_flags) {
      break;
    }
    write_iter = strcpya_k(write_iter, " + ");
  } while (load_filter_log_flags);
  *write_iter_ptr = write_iter;
}

BoolErr BigstackAllocPgv(uint32_t sample_ct, uint32_t multiallelic_needed, PgenGlobalFlags gflags, PgenVariant* pgvp) {
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  if (unlikely(bigstack_allocv_w(sample_ctl2, &(pgvp->genovec)))) {
    return 1;
  }
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  if (multiallelic_needed) {
    if (unlikely(bigstack_allocv_w(sample_ctl, &(pgvp->patch_01_set)) ||
                 bigstack_allocv_ac(sample_ct, &(pgvp->patch_01_vals)) ||
                 bigstack_allocv_w(sample_ctl, &(pgvp->patch_10_set)) ||
                 bigstack_allocv_ac(sample_ct * 2, &(pgvp->patch_10_vals)))) {
      return 1;
    }
  } else {
    // defensive
    pgvp->patch_01_set = nullptr;
    pgvp->patch_01_vals = nullptr;
    pgvp->patch_10_set = nullptr;
    pgvp->patch_10_vals = nullptr;
  }
  if (gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent)) {
    if (unlikely(bigstack_allocv_w(sample_ctl, &(pgvp->phasepresent)) ||
                 bigstack_allocv_w(sample_ctl, &(pgvp->phaseinfo)))) {
      return 1;
    }
  } else {
    pgvp->phasepresent = nullptr;
    pgvp->phaseinfo = nullptr;
  }
  if (gflags & kfPgenGlobalDosagePresent) {
    if (unlikely(bigstack_allocv_w(sample_ctl, &(pgvp->dosage_present)) ||
                 bigstack_allocv_dosage(sample_ct, &(pgvp->dosage_main)))) {
      return 1;
    }
    if (multiallelic_needed) {
      // todo
    }
    if (gflags & kfPgenGlobalDosagePhasePresent) {
      if (unlikely(bigstack_allocv_w(sample_ctl, &(pgvp->dphase_present)) ||
                   bigstack_allocv_dphase(sample_ct, &(pgvp->dphase_delta)))) {
        return 1;
      }
      if (multiallelic_needed) {
        // todo
      }
    }
  } else {
    pgvp->dosage_present = nullptr;
    pgvp->dosage_main = nullptr;
    pgvp->dphase_present = nullptr;
    pgvp->dphase_delta = nullptr;
    // todo: multiallelic-dosage buffers
  }
  return BigstackBaseSetChecked(g_bigstack_base);
}

static_assert(kDosageMid == 16384, "PrintDosageDecimal() needs to be updated.");
char* PrintDosageDecimal(uint32_t remainder, char* start) {
  // Instead of constant 5-digit precision, we print fewer digits whenever that
  // doesn't interfere with proper round-tripping.  I.e. we search for the
  // shortest string in
  //   ((n - 0.5)/16384, (n + 0.5)/16384).
  // E.g. 3277/16384 is 0.20001 when printed with 5-digit precision, but we'd
  // print that as 0.2 since that's still in (3276.5/16384, 3277.5/16384).
  *start++ = '.';
  // (remainder * 2) is in 32768ths
  // 32768 * 625 = 20480k, smallest common denominator with 10^4

  const uint32_t range_top_20480k = (remainder * 2 + 1) * 625;
  // this is technically checking a half-open rather than a fully-open
  // interval, but that's fine since we never hit the boundary points
  if ((range_top_20480k % 2048) < 1250) {
    // when this is true, the four-decimal-place approximation is in the range
    // which round-trips back to our original number.
    const uint32_t four_decimal_places = range_top_20480k / 2048;
    return u32toa_trunc4(four_decimal_places, start);
  }

  // we wish to print (100000 * remainder + 8192) / 16384, left-0-padded.  and
  // may as well banker's round too.
  //
  // banker's rounding yields a different result than regular rounding for n/64
  // when n is congruent to 1 mod 4:
  //   1/64 = .015625 -> print 0.01562
  //   3/64 = .046875 -> print 0.04688
  //   5/64 = .078125 -> print 0.07812
  const uint32_t five_decimal_places = ((3125 * remainder + 256) / 512) - ((remainder % 1024) == 256);
  const uint32_t first_decimal_place = five_decimal_places / 10000;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_four_digits = five_decimal_places - first_decimal_place * 10000;
  if (last_four_digits) {
    return u32toa_trunc4(last_four_digits, start);
  }
  return start;
}

static_assert(kDosageMax == 32768, "ddosagetoa() needs to be updated.");
char* ddosagetoa(uint64_t dosage, char* start) {
  // 3 digit precision seems like the best compromise between accuracy and
  // avoidance of rounding ugliness
  // (Rounding ugliness is not actually hidden for e.g. 1000 Genomes phase 1,
  // since there are lots of 0.05 and 0.1 dosages which all get rounded in the
  // same direction; oh well.)

  // +16 since we need to round .99951 up to 1
  const uint64_t dosage_p16 = dosage + 16;
  start = u32toa(dosage_p16 / kDosageMax, start);
  const uint32_t remainder_p16 = S_CAST(uint32_t, dosage_p16) % kDosageMax;
  if (remainder_p16 < 33) {
    return start;
  }
  // (1000 * remainder + 16384) / 32768
  //   1/16 = .0625 -> print 0.062
  //   3/16 = .1875 -> print 0.188
  //   5/16 = .3125 -> print 0.312
  // const uint32_t three_decimal_places = ((125 * remainder + 2048) / 4096) - ((remainder % 8192) == 2048);
  const uint32_t three_decimal_places = ((125 * remainder_p16 + 48) / 4096) - ((remainder_p16 % 8192) == 4048);
  // three_decimal_places guaranteed to be nonzero here
  *start++ = '.';
  const uint32_t first_decimal_place = three_decimal_places / 100;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_two_digits = three_decimal_places - first_decimal_place * 100;
  if (last_two_digits) {
    memcpy_k(start, &(kDigitPair[last_two_digits]), 2);
    return &(start[1 + (start[1] != '0')]);
  }
  return start;
}

static_assert(kDosageMax == 32768, "PrintDdosageDecimal() needs to be updated.");
char* PrintDdosageDecimal(uint32_t remainder, char* start) {
  // Instead of constant 5-digit precision, we print fewer digits whenever that
  // doesn't interfere with proper round-tripping.  I.e. we search for the
  // shortest string in
  //   ((n - 0.5)/32768, (n + 0.5)/32768).
  // E.g. 6554/32768 is 0.20001 when printed with 5-digit precision, but we'd
  // print that as 0.2 since that's still in (6553.5/32768, 6554.5/32768).
  // (This is more precise than necessary for most chromosomes, but it's
  // necessary for chrX.)
  *start++ = '.';
  // (remainder * 2) is in 65536ths
  // 65536 * 625 = 40960k, smallest common denominator with 10^4

  const uint32_t range_top_40960k = remainder * 1250 + 625;
  // this is technically checking a half-open rather than a fully-open
  // interval, but that's fine since we never hit the boundary points
  if ((range_top_40960k % 4096) < 1250) {
    // when this is true, the four-decimal-place approximation is in the range
    // which round-trips back to our original number.
    const uint32_t four_decimal_places = range_top_40960k / 4096;
    return u32toa_trunc4(four_decimal_places, start);
  }

  // we wish to print (100000 * remainder + 16384) / 32768, left-0-padded.  and
  // may as well banker's round too.
  //
  // banker's rounding yields a different result than regular rounding for n/64
  // when n is congruent to 1 mod 4:
  //   1/64 = .015625 -> print 0.01562
  //   3/64 = .046875 -> print 0.04688
  //   5/64 = .078125 -> print 0.07812
  const uint32_t five_decimal_places = ((3125 * remainder + 512) / 1024) - ((remainder % 2048) == 512);
  const uint32_t first_decimal_place = five_decimal_places / 10000;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_four_digits = five_decimal_places - first_decimal_place * 10000;
  if (last_four_digits) {
    return u32toa_trunc4(last_four_digits, start);
  }
  return start;
}

static_assert(sizeof(Dosage) == 2, "PopulateDenseDosage() needs to be updated.");
void PopulateDenseDosage(const uintptr_t* genoarr, const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t sample_ct, uint32_t dosage_ct, Dosage* dense_dosage) {
  GenoarrLookup256x2bx4(genoarr, kHcToDosage16, sample_ct, dense_dosage);
  // fill trailing bits with missing values so vector operations work
  const uint32_t trailing_entry_ct = (-sample_ct) % kDosagePerVec;
  if (trailing_entry_ct) {
    Dosage* dense_dosage_end = &(dense_dosage[sample_ct]);
    for (uint32_t uii = 0; uii != trailing_entry_ct; ++uii) {
      dense_dosage_end[uii] = kDosageMissing;
    }
  }
  if (dosage_ct) {
    uintptr_t sample_idx_base = 0;
    uintptr_t cur_bits = dosage_present[0];
    for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
      const uintptr_t sample_idx = BitIter1(dosage_present, &sample_idx_base, &cur_bits);
      dense_dosage[sample_idx] = dosage_main[dosage_idx];
    }
  }
}

// workspace must have at least NypCtToWordCt(sample_ct) words
static_assert(sizeof(Dosage) == 2, "PopulateDenseDosageNonemptySubset() needs to be updated.");
void PopulateDenseDosageNonemptySubset(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* genoarr, const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t dosage_ct, Dosage* dense_dosage, uintptr_t* workspace) {
  CopyNyparrNonemptySubset(genoarr, sample_include, raw_sample_ct, sample_ct, workspace);
  GenoarrLookup256x2bx4(workspace, kHcToDosage16, sample_ct, dense_dosage);
  // fill trailing bits with missing values so vector operations work
  const uint32_t trailing_entry_ct = (-sample_ct) % kDosagePerVec;
  if (trailing_entry_ct) {
    Dosage* dense_dosage_end = &(dense_dosage[sample_ct]);
    for (uint32_t uii = 0; uii != trailing_entry_ct; ++uii) {
      dense_dosage_end[uii] = kDosageMissing;
    }
  }
  if (dosage_ct) {
    uintptr_t widx = 0;
    uintptr_t cur_bits = dosage_present[0];
    for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
      const uintptr_t shifted_bit = BitIter1y(dosage_present, &widx, &cur_bits);
      const uintptr_t sample_include_word = sample_include[widx];
      if (sample_include_word & shifted_bit) {
        const uint32_t sample_idx = sample_include_cumulative_popcounts[widx] + PopcountWord(sample_include_word & (shifted_bit - 1));
        dense_dosage[sample_idx] = dosage_main[dosage_idx];
      }
    }
  }
}

void PopulateRescaledDosage(const uintptr_t* genoarr, const uintptr_t* dosage_present, const Dosage* dosage_main, double slope, double intercept, double missing_val, uint32_t sample_ct, uint32_t dosage_ct, double* expanded_dosages) {
  double lookup_vals[32] ALIGNV16;
  lookup_vals[0] = intercept;
  lookup_vals[2] = intercept + slope;
  lookup_vals[4] = intercept + 2 * slope;
  lookup_vals[6] = missing_val;
  InitLookup16x8bx2(lookup_vals);
  GenoarrLookup16x8bx2(genoarr, lookup_vals, sample_ct, expanded_dosages);
  if (!dosage_ct) {
    return;
  }
  slope *= kRecipDosageMid;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
    expanded_dosages[sample_uidx] = dosage_main[dosage_idx] * slope + intercept;
  }
}

void PopulateRescaledDosageF(const uintptr_t* genoarr, const uintptr_t* dosage_present, const Dosage* dosage_main, float slope, float intercept, float missing_val, uint32_t sample_ct, uint32_t dosage_ct, float* expanded_dosages) {
  float lookup_vals[32] ALIGNV16;
  lookup_vals[0] = intercept;
  lookup_vals[2] = intercept + slope;
  lookup_vals[4] = intercept + 2 * slope;
  lookup_vals[6] = missing_val;
  InitLookup16x4bx2(lookup_vals);
  GenoarrLookup16x4bx2(genoarr, lookup_vals, sample_ct, expanded_dosages);
  if (!dosage_ct) {
    return;
  }
  slope *= S_CAST(float, kRecipDosageMid);
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
    expanded_dosages[sample_uidx] = dosage_main[dosage_idx] * slope + intercept;
  }
}

void DosagePhaseinfoPatch(const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const uintptr_t* dosage_present, const Dosage* dosage_vec, const uintptr_t* dphase_present, uint32_t sample_ct, SDosage* dphase_delta) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  uintptr_t dosage_present_word = 0;
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    uintptr_t phasepresent_nodphase_word = phasepresent[widx];
    if (dphase_present) {
      phasepresent_nodphase_word &= ~dphase_present[widx];
    }
    if (phasepresent_nodphase_word) {
      const uintptr_t phaseinfo_word = phaseinfo[widx];
      if (dosage_present) {
        dosage_present_word = dosage_present[widx];
      }
      const uint32_t sample_idx_offset = widx * kBitsPerWord;
      do {
        const uint32_t sample_idx_lowbits = ctzw(phasepresent_nodphase_word);
        const uint32_t sample_idx = sample_idx_offset + sample_idx_lowbits;
        // should compile to blsi
        const uintptr_t shifted_bit = phasepresent_nodphase_word & (-phasepresent_nodphase_word);
        int32_t cur_diff = kDosageMid;
        if (dosage_present_word & shifted_bit) {
          cur_diff = DosageHomdist(dosage_vec[sample_idx]);
        }
        // todo: verify the compiler optimizes this well
        if (!(phaseinfo_word & shifted_bit)) {
          cur_diff = -cur_diff;
        }
        dphase_delta[sample_idx] = cur_diff;
        phasepresent_nodphase_word ^= shifted_bit;
      } while (phasepresent_nodphase_word);
    }
  }
}

void PopulateDenseDphase(const uintptr_t* phasepresent, const uintptr_t* phaseinfo, const uintptr_t* dosage_present, const Dosage* dense_dosage_vec, const uintptr_t* dphase_present, const SDosage* dphase_delta, uint32_t sample_ct, uint32_t phasepresent_ct, uint32_t dosage_ct, uint32_t dphase_ct, SDosage* dense_dphase_delta) {
  const uint32_t dosagev_ct = DivUp(sample_ct, kDosagePerVec);
  ZeroDphaseArr(dosagev_ct * kDosagePerVec, dense_dphase_delta);
  if (dphase_ct) {
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = dphase_present[0];
    for (uint32_t dphase_idx = 0; dphase_idx != dphase_ct; ++dphase_idx) {
      const uintptr_t sample_uidx = BitIter1(dphase_present, &sample_uidx_base, &cur_bits);
      const int32_t dphase_delta_val = dphase_delta[dphase_idx];
      dense_dphase_delta[sample_uidx] = dphase_delta_val;
    }
  } else {
    dphase_present = nullptr;
  }
  if (phasepresent_ct) {
    DosagePhaseinfoPatch(phasepresent, phaseinfo, dosage_ct? dosage_present : nullptr, dense_dosage_vec, dphase_present, sample_ct, dense_dphase_delta);
  }
}

static_assert(sizeof(AlleleCode) == 1, "AtLeastOneMultiallelicHet() needs to be updated.");
uint32_t AtLeastOneMultiallelicHet(const PgenVariant* pgvp, uint32_t sample_ct) {
  if (pgvp->patch_01_ct) {
    return 1;
  }
  {
    const uintptr_t* genoarr = pgvp->genovec;
    const uint32_t fullvec_ct = sample_ct / kBitsPerWordD2;
    for (uint32_t uii = 0; uii != fullvec_ct; ++uii) {
      const uintptr_t geno_word = genoarr[uii];
      if (Word01(geno_word)) {
        return 1;
      }
    }
    const uint32_t remainder = sample_ct % kBitsPerWordD2;
    if (remainder) {
      if (Word01(bzhi(genoarr[fullvec_ct], 2 * remainder))) {
        return 1;
      }
    }
  }
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  if (patch_10_ct) {
    const AlleleCode* patch_10_vals = pgvp->patch_10_vals;
#ifdef __LP64__
    const uint32_t fullvec_ct = patch_10_ct / (kBytesPerVec / 2);
    const VecU16* patch_10_valias = R_CAST(const VecU16*, patch_10_vals);
    const VecU16 m8 = vecu16_set1(0xff);
    for (uint32_t vidx = 0; vidx != fullvec_ct; ++vidx) {
      const VecU16 vv_orig = vecu16_loadu(&(patch_10_valias[vidx]));
      const VecU16 vv_hi = vecu16_srli(vv_orig, 8);
      const VecU16 vv_lo = vv_orig & m8;
      if (!vecu16s_are_equal(vv_hi, vv_lo)) {
        return 1;
      }
    }
    uint32_t uii = fullvec_ct * (kBytesPerVec / 2);
#else
    uint32_t uii = 0;
#endif
    for (; uii != patch_10_ct; ++uii) {
      if (patch_10_vals[2 * uii] != patch_10_vals[2 * uii + 1]) {
        return 1;
      }
    }
  }
  return 0;
}

void SetHetMissing(uintptr_t word_ct, uintptr_t* genovec) {
  // 01 -> 11, nothing else changes
#ifdef __LP64__
  const VecW m1 = VCONST_W(kMask5555);
  VecW* geno_vvec_iter = R_CAST(VecW*, genovec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    const VecW cur_geno_vword = *geno_vvec_iter;
    const VecW cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
  }
  if (full_vec_ct & 2) {
    VecW cur_geno_vword = *geno_vvec_iter;
    VecW cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    VecW cur_geno_vword = *geno_vvec_iter;
    VecW cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = cur_geno_vword | cur_geno_vword_low_lshifted;
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    uintptr_t geno_word = genovec[base_idx];
    genovec[base_idx] = geno_word | ((geno_word & kMask5555) << 1);
    geno_word = genovec[base_idx + 1];
    genovec[base_idx + 1] = geno_word | ((geno_word & kMask5555) << 1);
  }
#  endif
  if (word_ct & 1) {
    const uintptr_t geno_word = genovec[word_ct - 1];
    genovec[word_ct - 1] = geno_word | ((geno_word & kMask5555) << 1);
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    genovec[widx] = geno_word | ((geno_word & kMask5555) << 1);
  }
#endif
}

void SetHetMissingCleardosage(uintptr_t word_ct, uintptr_t* __restrict genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main) {
  const uint32_t orig_write_dosage_ct = *write_dosage_ct_ptr;
  if (orig_write_dosage_ct) {
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = dosagepresent[0];
    for (uint32_t dosage_read_idx = 0; dosage_read_idx != orig_write_dosage_ct; ++dosage_read_idx) {
      const uintptr_t sample_uidx = BitIter1(dosagepresent, &sample_uidx_base, &cur_bits);
      if (GetNyparrEntry(genovec, sample_uidx) == 1) {
        ClearBit(sample_uidx, dosagepresent);
        uint32_t dosage_write_idx = dosage_read_idx++;
        for (; dosage_read_idx == orig_write_dosage_ct; ++dosage_read_idx) {
          const uintptr_t sample_uidx2 = BitIter1(dosagepresent, &sample_uidx_base, &cur_bits);
          if (GetNyparrEntry(genovec, sample_uidx2) == 1) {
            ClearBit(sample_uidx2, dosagepresent);
          } else {
            dosage_main[dosage_write_idx++] = dosage_main[dosage_read_idx];
          }
        }
        *write_dosage_ct_ptr = dosage_write_idx;
        break;
      }
    }
  }
  SetHetMissing(word_ct, genovec);
}

void SetHetMissingKeepdosage(uintptr_t word_ct, uintptr_t* __restrict genovec, uint32_t* write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main) {
  // count # of 1.00000 dosages we need to insert, and then rewrite dosage_main
  // from back to front so we don't need temporary buffers.

  const uint32_t orig_dosage_ct = *write_dosage_ct_ptr;
  // can't assume dosagepresent is initialized in this case
  if (!orig_dosage_ct) {
    ZeroWArr(DivUp(word_ct, 2), dosagepresent);
  }
  Halfword* dosagepresent_alias = R_CAST(Halfword*, dosagepresent);
  uint32_t new_dosage_ct = 0;
  // todo: try vectorizing this loop
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t geno_hets = Word01(genovec[widx]);
    const uintptr_t dosagepresent_word = UnpackHalfwordToWord(dosagepresent_alias[widx]);
    new_dosage_ct += Popcount01Word(geno_hets & (~dosagepresent_word));
  }
  if (!new_dosage_ct) {
    SetHetMissing(word_ct, genovec);
    return;
  }
  uint32_t dosage_write_idx = orig_dosage_ct + new_dosage_ct;
  uint32_t dosage_read_idx = orig_dosage_ct;
  uint32_t widx = word_ct;
  *write_dosage_ct_ptr = dosage_write_idx;
  do {
    --widx;
    const uintptr_t geno_word = genovec[widx];
    const uintptr_t dosagepresent_word = UnpackHalfwordToWord(dosagepresent_alias[widx]);
    const uintptr_t geno_hets = Word01(geno_word);
    uintptr_t new_dosagepresent_word = dosagepresent_word | geno_hets;
    if (new_dosagepresent_word) {
      dosagepresent_alias[widx] = PackWordToHalfword(new_dosagepresent_word);
      do {
        const uint32_t top_set_bit = bsrw(new_dosagepresent_word);
        const uintptr_t cur_bit_word = k1LU << top_set_bit;
        Dosage cur_dosage = kDosageMid;
        if (cur_bit_word & dosagepresent_word) {
          cur_dosage = dosage_main[--dosage_read_idx];
        }
        dosage_main[--dosage_write_idx] = cur_dosage;
        new_dosagepresent_word -= cur_bit_word;
      } while (new_dosagepresent_word);
      genovec[widx] = geno_word | (geno_hets << 1);
    }
  } while (widx);
}

void SetMissingRef(uintptr_t word_ct, uintptr_t* genovec) {
  // 11 -> 00, nothing else changes
#ifdef __LP64__
  const VecW m1 = VCONST_W(kMask5555);
  VecW* geno_vvec_iter = R_CAST(VecW*, genovec);
  const uintptr_t full_vec_ct = word_ct / kWordsPerVec;
  if (full_vec_ct & 1) {
    const VecW cur_geno_vword = *geno_vvec_iter;
    const VecW cur_geno_vword_high_rshifted = vecw_srli(cur_geno_vword, 1) & m1;
    const VecW cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = vecw_and_notfirst(cur_geno_vword_high_rshifted | cur_geno_vword_low_lshifted, cur_geno_vword);
  }
  if (full_vec_ct & 2) {
    VecW cur_geno_vword = *geno_vvec_iter;
    VecW cur_geno_vword_high_rshifted = vecw_srli(cur_geno_vword, 1) & m1;
    VecW cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = vecw_and_notfirst(cur_geno_vword_high_rshifted | cur_geno_vword_low_lshifted, cur_geno_vword);
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_high_rshifted = vecw_srli(cur_geno_vword, 1) & m1;
    cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = vecw_and_notfirst(cur_geno_vword_high_rshifted | cur_geno_vword_low_lshifted, cur_geno_vword);
  }
  for (uintptr_t ulii = 3; ulii < full_vec_ct; ulii += 4) {
    VecW cur_geno_vword = *geno_vvec_iter;
    VecW cur_geno_vword_high_rshifted = vecw_srli(cur_geno_vword, 1) & m1;
    VecW cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = vecw_and_notfirst(cur_geno_vword_high_rshifted | cur_geno_vword_low_lshifted, cur_geno_vword);
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_high_rshifted = vecw_srli(cur_geno_vword, 1) & m1;
    cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = vecw_and_notfirst(cur_geno_vword_high_rshifted | cur_geno_vword_low_lshifted, cur_geno_vword);
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_high_rshifted = vecw_srli(cur_geno_vword, 1) & m1;
    cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = vecw_and_notfirst(cur_geno_vword_high_rshifted | cur_geno_vword_low_lshifted, cur_geno_vword);
    cur_geno_vword = *geno_vvec_iter;
    cur_geno_vword_high_rshifted = vecw_srli(cur_geno_vword, 1) & m1;
    cur_geno_vword_low_lshifted = vecw_slli(cur_geno_vword & m1, 1);
    *geno_vvec_iter++ = vecw_and_notfirst(cur_geno_vword_high_rshifted | cur_geno_vword_low_lshifted, cur_geno_vword);
  }
#  ifdef USE_AVX2
  if (word_ct & 2) {
    const uintptr_t base_idx = full_vec_ct * kWordsPerVec;
    uintptr_t geno_word = genovec[base_idx];
    uintptr_t both_set = geno_word & (geno_word >> 1) & kMask5555;
    genovec[base_idx] = geno_word ^ (both_set * 3);
    geno_word = genovec[base_idx + 1];
    both_set = geno_word & (geno_word >> 1) & kMask5555;
    genovec[base_idx + 1] = geno_word ^ (both_set * 3);
  }
#  endif
  if (word_ct & 1) {
    const uintptr_t geno_word = genovec[word_ct - 1];
    const uintptr_t both_set = geno_word & (geno_word >> 1) & kMask5555;
    genovec[word_ct - 1] = geno_word ^ (both_set * 3);
  }
#else
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    const uintptr_t both_set = geno_word & (geno_word >> 1) & kMask5555;
    genovec[widx] = geno_word ^ (both_set * 3);
  }
#endif
}

void SetMissingRefY(const uintptr_t* __restrict sex_female_interleaved, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec) {
  const uintptr_t twovec_ct = geno_vec_ct / 2;
#ifdef __LP64__
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m1inv = VCONST_W(kMaskAAAA);
  const VecW* sex_female_interleaved_iter = R_CAST(const VecW*, sex_female_interleaved);
  VecW* genovvec_iter = R_CAST(VecW*, genovec);
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const VecW set_vvec = *sex_female_interleaved_iter++;
    VecW cur_geno_vword = *genovvec_iter;
    const VecW nonfemale_first = vecw_and_notfirst(set_vvec, m1);
    VecW high_nonfemale_rshifted = vecw_srli(cur_geno_vword, 1) & nonfemale_first;
    VecW low_nonfemale_lshifted = vecw_slli(cur_geno_vword & nonfemale_first, 1);
    *genovvec_iter++ = vecw_and_notfirst(high_nonfemale_rshifted | low_nonfemale_lshifted, cur_geno_vword);
    const VecW nonfemale_second_lshifted = vecw_and_notfirst(set_vvec, m1inv);
    cur_geno_vword = *genovvec_iter;
    high_nonfemale_rshifted = vecw_srli(cur_geno_vword & nonfemale_second_lshifted, 1);
    low_nonfemale_lshifted = vecw_slli(cur_geno_vword, 1) & nonfemale_second_lshifted;
    *genovvec_iter++ = vecw_and_notfirst(high_nonfemale_rshifted | low_nonfemale_lshifted, cur_geno_vword);
  }
  if (geno_vec_ct & 1) {
    const VecW set_first = *sex_female_interleaved_iter;
    const VecW cur_geno_vword = *genovvec_iter;
    const VecW nonfemale_first = vecw_and_notfirst(set_first, m1);
    const VecW high_nonfemale_rshifted = vecw_srli(cur_geno_vword, 1) & nonfemale_first;
    const VecW low_nonfemale_lshifted = vecw_slli(cur_geno_vword & nonfemale_first, 1);
    *genovvec_iter = vecw_and_notfirst(high_nonfemale_rshifted | low_nonfemale_lshifted, cur_geno_vword);
  }
#else
  const uintptr_t* sex_female_interleaved_iter = sex_female_interleaved;
  uintptr_t* genovec_iter = genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const uintptr_t set_word = *sex_female_interleaved_iter++;
    uintptr_t geno_word = *genovec_iter;
    uintptr_t both_set_and_nonfemale = geno_word & (geno_word >> 1) & kMask5555 & (~set_word);
    *genovec_iter++ = geno_word ^ (both_set_and_nonfemale * 3);
    geno_word = *genovec_iter;
    both_set_and_nonfemale = geno_word & kMask5555 & ((geno_word & (~set_word)) >> 1);
    *genovec_iter++ = geno_word ^ (both_set_and_nonfemale * 3);
  }
  if (geno_vec_ct & 1) {
    const uintptr_t set_word = *sex_female_interleaved_iter;
    const uintptr_t geno_word = *genovec_iter;
    const uintptr_t both_set_and_nonfemale  = geno_word & (geno_word >> 1) & kMask5555 & (~set_word);
    *genovec_iter = geno_word ^ (both_set_and_nonfemale * 3);
  }
#endif
}

void SetMissingRefDosage(const uintptr_t* __restrict dosagepresent, uintptr_t sample_ct, uintptr_t* __restrict genovec) {
  const uintptr_t word_ct = NypCtToWordCt(sample_ct);
  const Halfword* dosagepresent_hw = DowncastKWToHW(dosagepresent);
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    const Halfword dosage_hw = dosagepresent_hw[widx];
    uintptr_t dosage_missing = geno_word & (geno_word >> 1) & kMask5555;
#ifdef USE_AVX2
    dosage_missing &= ~UnpackHalfwordToWord(dosage_hw);
#else
    if (!dosage_hw) {
      dosage_missing &= ~UnpackHalfwordToWord(dosage_hw);
    }
#endif
    genovec[widx] = geno_word ^ (dosage_missing * 3);
  }
}

void SetMissingRefDosageY(const uintptr_t* __restrict sex_female, const uintptr_t* __restrict dosagepresent, uintptr_t sample_ct, uintptr_t* __restrict genovec) {
  const uintptr_t word_ct = NypCtToWordCt(sample_ct);
  const Halfword* sex_female_hw = DowncastKWToHW(sex_female);
  const Halfword* dosagepresent_hw = DowncastKWToHW(dosagepresent);
  for (uintptr_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    const Halfword female_or_dosage_hw = sex_female_hw[widx] | dosagepresent_hw[widx];
    uintptr_t dosage_missing = geno_word & (geno_word >> 1) & kMask5555;
    dosage_missing &= ~UnpackHalfwordToWord(female_or_dosage_hw);
    genovec[widx] = geno_word ^ (dosage_missing * 3);
  }
}

// These functions may move to pgenlib_misc later.
uint32_t CountMissingVec6(const VecW* __restrict geno_vvec, uint32_t vec_ct) {
  assert(!(vec_ct % 6));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* geno_vvec_iter = geno_vvec;
  VecW acc_bothset = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        return HsumW(acc_bothset);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_bothset = vecw_setzero();
    const VecW* geno_vvec_stop = &(geno_vvec_iter[cur_incr]);
    do {
      VecW cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW bothset1 = m1 & cur_geno_vword & vecw_srli(cur_geno_vword, 1);

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      bothset1 = bothset1 + (m1 & cur_geno_vword & vecw_srli(cur_geno_vword, 1));

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      bothset1 = bothset1 + (m1 & cur_geno_vword & vecw_srli(cur_geno_vword, 1));
      bothset1 = (bothset1 & m2) + (vecw_srli(bothset1, 2) & m2);

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW bothset2 = m1 & cur_geno_vword & vecw_srli(cur_geno_vword, 1);

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      bothset2 = bothset2 + (m1 & cur_geno_vword & vecw_srli(cur_geno_vword, 1));

      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      bothset2 = bothset2 + (m1 & cur_geno_vword & vecw_srli(cur_geno_vword, 1));

      bothset1 = bothset1 + (bothset2 & m2) + (vecw_srli(bothset2, 2) & m2);
      // these now contain 4-bit values from 0-12

      inner_acc_bothset = inner_acc_bothset + (bothset1 & m4) + (vecw_srli(bothset1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_bothset = acc_bothset + vecw_bytesum(inner_acc_bothset, m0);
  }
}

uint32_t GenoarrCountMissingUnsafe(const uintptr_t* genoarr, uint32_t sample_ct) {
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uint32_t word_idx = sample_ctl2 - (sample_ctl2 % (6 * kWordsPerVec));
  uint32_t missing_ct = CountMissingVec6(R_CAST(const VecW*, genoarr), word_idx / kWordsPerVec);
  for (; word_idx != sample_ctl2; ++word_idx) {
    uintptr_t ww = genoarr[word_idx];
    // bugfix (19 Oct 2020): forgot to mask out high bits
    ww = ww & (ww >> 1) & kMask5555;
    missing_ct += Popcount01Word(ww);
  }
  return missing_ct;
}

// geno_vvec allowed to be unaligned.
uint32_t CountMissingMaskedVec6(const VecW* __restrict geno_vvec, const VecW* __restrict interleaved_mask_vvec, uint32_t vec_ct) {
  assert(!(vec_ct % 6));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* geno_vvec_iter = geno_vvec;
  const VecW* interleaved_mask_vvec_iter = interleaved_mask_vvec;
  VecW acc_bothset = vecw_setzero();
  uintptr_t cur_incr = 60;
  for (; ; vec_ct -= cur_incr) {
    if (vec_ct < 60) {
      if (!vec_ct) {
        return HsumW(acc_bothset);
      }
      cur_incr = vec_ct;
    }
    VecW inner_acc_bothset = vecw_setzero();
    const VecW* geno_vvec_stop = &(geno_vvec_iter[cur_incr]);
    do {
      VecW interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      VecW cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW cur_mask = interleaved_mask_vword & m1;
      VecW bothset1 = cur_mask & cur_geno_vword & vecw_srli(cur_geno_vword, 1);

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      bothset1 = bothset1 + (cur_mask & cur_geno_vword & vecw_srli(cur_geno_vword, 1));

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      bothset1 = bothset1 + (cur_mask & cur_geno_vword & vecw_srli(cur_geno_vword, 1));
      bothset1 = (bothset1 & m2) + (vecw_srli(bothset1, 2) & m2);

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      VecW bothset2 = cur_mask & cur_geno_vword & vecw_srli(cur_geno_vword, 1);

      interleaved_mask_vword = *interleaved_mask_vvec_iter++;
      cur_mask = interleaved_mask_vword & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      bothset2 = bothset2 + (cur_mask & cur_geno_vword & vecw_srli(cur_geno_vword, 1));

      cur_mask = vecw_srli(interleaved_mask_vword, 1) & m1;
      cur_geno_vword = vecw_loadu(geno_vvec_iter);
      ++geno_vvec_iter;
      bothset2 = bothset2 + (cur_mask & cur_geno_vword & vecw_srli(cur_geno_vword, 1));

      bothset1 = bothset1 + (bothset2 & m2) + (vecw_srli(bothset2, 2) & m2);
      // these now contain 4-bit values from 0-12

      inner_acc_bothset = inner_acc_bothset + (bothset1 & m4) + (vecw_srli(bothset1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    const VecW m0 = vecw_setzero();
    acc_bothset = acc_bothset + vecw_bytesum(inner_acc_bothset, m0);
  }
}

uint32_t GenoarrCountMissingSubset(const uintptr_t* genoarr, const uintptr_t* interleaved_vec, uint32_t sample_ct) {
  // See GenoarrCountSubsetFreqs().
  const uint32_t sample_ctv2 = NypCtToVecCt(sample_ct);
  uint32_t missing_ct;
#ifdef __LP64__
  uint32_t vec_idx = sample_ctv2 - (sample_ctv2 % 6);
  missing_ct = CountMissingMaskedVec6(R_CAST(const VecW*, genoarr), R_CAST(const VecW*, interleaved_vec), vec_idx);
  const uintptr_t* genoarr_iter = &(genoarr[kWordsPerVec * vec_idx]);
  const uintptr_t* interleaved_mask_iter = &(interleaved_vec[(kWordsPerVec / 2) * vec_idx]);
#  ifdef USE_AVX2
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  uintptr_t mask_base3 = 0;
  uintptr_t mask_base4 = 0;
  for (; vec_idx != sample_ctv2; ++vec_idx) {
    uintptr_t mask_word1;
    uintptr_t mask_word2;
    uintptr_t mask_word3;
    uintptr_t mask_word4;
    if (!(vec_idx % 2)) {
      mask_base1 = *interleaved_mask_iter++;
      mask_base2 = *interleaved_mask_iter++;
      mask_base3 = *interleaved_mask_iter++;
      mask_base4 = *interleaved_mask_iter++;
      mask_word1 = mask_base1 & kMask5555;
      mask_word2 = mask_base2 & kMask5555;
      mask_word3 = mask_base3 & kMask5555;
      mask_word4 = mask_base4 & kMask5555;
    } else {
      mask_word1 = (mask_base1 >> 1) & kMask5555;
      mask_word2 = (mask_base2 >> 1) & kMask5555;
      mask_word3 = (mask_base3 >> 1) & kMask5555;
      mask_word4 = (mask_base4 >> 1) & kMask5555;
    }
    for (uint32_t vechalf_idx = 0; ; ++vechalf_idx) {
      const uintptr_t cur_geno_word1 = *genoarr_iter++;
      const uintptr_t cur_geno_word2 = *genoarr_iter++;
      const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
      const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
      missing_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
      if (vechalf_idx) {
        break;
      }
      mask_word1 = mask_word3;
      mask_word2 = mask_word4;
    }
  }
#  else  // not USE_AVX2
  uintptr_t mask_base1 = 0;
  uintptr_t mask_base2 = 0;
  for (; vec_idx != sample_ctv2; ++vec_idx) {
    uintptr_t mask_word1;
    uintptr_t mask_word2;
    if (!(vec_idx % 2)) {
      mask_base1 = *interleaved_mask_iter++;
      mask_base2 = *interleaved_mask_iter++;
      mask_word1 = mask_base1 & kMask5555;
      mask_word2 = mask_base2 & kMask5555;
    } else {
      mask_word1 = (mask_base1 >> 1) & kMask5555;
      mask_word2 = (mask_base2 >> 1) & kMask5555;
    }
    const uintptr_t cur_geno_word1 = *genoarr_iter++;
    const uintptr_t cur_geno_word2 = *genoarr_iter++;
    const uintptr_t cur_geno_word1_high_masked = mask_word1 & (cur_geno_word1 >> 1);
    const uintptr_t cur_geno_word2_high_masked = mask_word2 & (cur_geno_word2 >> 1);
#    ifdef USE_SSE42
    missing_ct += PopcountWord(((cur_geno_word1 & cur_geno_word1_high_masked) << 1) | (cur_geno_word2 & cur_geno_word2_high_masked));
#    else
    missing_ct += NypsumWord((cur_geno_word1 & cur_geno_word1_high_masked) + (cur_geno_word2 & cur_geno_word2_high_masked));
#    endif
  }
#  endif  // not USE_AVX2
#else  // not __LP64__
  uint32_t word_idx = sample_ctv2 - (sample_ctv2 % 6);
  missing_ct = CountMissingMaskedVec6(R_CAST(const VecW*, genoarr), R_CAST(const VecW*, interleaved_vec), word_idx);
  const uintptr_t* interleaved_mask_iter = &(interleaved_vec[word_idx / 2]);
  uintptr_t mask_base = 0;
  for (; word_idx != sample_ctv2; ++word_idx) {
    uintptr_t mask_word;
    if (!(word_idx % 2)) {
      mask_base = *interleaved_mask_iter++;
      mask_word = mask_base & kMask5555;
    } else {
      mask_word = (mask_base >> 1) & kMask5555;
    }
    const uintptr_t cur_geno_word = genoarr[word_idx];
    const uintptr_t cur_geno_word_high_masked = mask_word & (cur_geno_word >> 1);
    missing_ct += Popcount01Word(cur_geno_word & cur_geno_word_high_masked);
  }
#endif
  return missing_ct;
}

// counts post-dosage missing, for which we don't have an interleaved mask
uint32_t GenoarrCountMissingInvsubsetUnsafe(const uintptr_t* genoarr, const uintptr_t* exclude_mask, uint32_t sample_ct) {
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uintptr_t* genoarr_iter = genoarr;
  const Halfword* exclude_alias_iter = R_CAST(const Halfword*, exclude_mask);
  uint32_t missing_ct = 0;
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    uintptr_t ww = *genoarr_iter++;
    ww = ww & (ww >> 1);
    const uint32_t include_mask = ~(*exclude_alias_iter++);
    missing_ct += Popcount01Word(ww & UnpackHalfwordToWord(include_mask));
  }
  return missing_ct;
}


void CollapsedSampleFmtidInit(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t include_fid, uint32_t include_sid, uintptr_t max_sample_fmtid_blen, char* collapsed_sample_fmtids_iter) {
  const char* sample_ids = siip->sample_ids;
  const char* sids = siip->sids;
  const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
  uintptr_t max_sid_blen = 0;
  if (include_sid) {
    max_sid_blen = sids? siip->max_sid_blen : 2;
  }
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
    if (!include_fid) {
      cur_sample_id = AdvPastDelim(cur_sample_id, '\t');
    }
    char* write_iter = strcpya(collapsed_sample_fmtids_iter, cur_sample_id);
    if (include_sid) {
      *write_iter++ = '\t';
      if (sids) {
        strcpy(write_iter, &(sids[sample_uidx * max_sid_blen]));
      } else {
        strcpy_k(write_iter, "0");
      }
    } else {
      *write_iter = '\0';
    }
    collapsed_sample_fmtids_iter = &(collapsed_sample_fmtids_iter[max_sample_fmtid_blen]);
  }
}

BoolErr CollapsedSampleFmtidInitAlloc(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t include_fid, uint32_t include_sid, char** collapsed_sample_fmtids_ptr, uintptr_t* max_sample_fmtid_blen_ptr) {
  const uintptr_t max_sample_fmtid_blen = GetMaxSampleFmtidBlen(siip, include_fid, include_sid);
  *max_sample_fmtid_blen_ptr = max_sample_fmtid_blen;
  // might add sample_fmtid_map later
  if (unlikely(bigstack_alloc_c(max_sample_fmtid_blen * sample_ct, collapsed_sample_fmtids_ptr))) {
    return 1;
  }
  CollapsedSampleFmtidInit(sample_include, siip, sample_ct, include_fid, include_sid, max_sample_fmtid_blen, *collapsed_sample_fmtids_ptr);
  return 0;
}

uint32_t OnlyOneFid(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct) {
  if (!(siip->flags & kfSampleIdFidPresent)) {
    return 1;
  }
  const char* sample_ids = siip->sample_ids;
  const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  const uintptr_t first_sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
  const char* fid_start = &(sample_ids[first_sample_uidx * max_sample_id_blen]);
  const uintptr_t fid_blen = AdvPastDelim(fid_start, '\t') - fid_start;
  for (uint32_t sample_idx = 1; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    const char* cur_fid_start = &(sample_ids[sample_uidx * max_sample_id_blen]);
    if (!memequal(fid_start, cur_fid_start, fid_blen)) {
      return 0;
    }
  }
  return 1;
}

uint32_t GetMajIdxMulti(const double* cur_allele_freqs, uint32_t cur_allele_ct) {
  assert(cur_allele_ct > 2);
  double max_freq = cur_allele_freqs[1];
  if (max_freq >= 0.5) {
    return 1;
  }
  double tot_nonlast_freq = cur_allele_freqs[0];
  uint32_t maj_allele_idx = 1;
  if (tot_nonlast_freq >= max_freq) {
    maj_allele_idx = 0;
    max_freq = tot_nonlast_freq;
  }
  tot_nonlast_freq += max_freq;
  const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
  for (uint32_t allele_idx = 2; allele_idx != cur_allele_ct_m1; ++allele_idx) {
    const double cur_freq = cur_allele_freqs[allele_idx];
    if (cur_freq > max_freq) {
      maj_allele_idx = allele_idx;
      max_freq = cur_freq;
    }
    tot_nonlast_freq += cur_freq;
  }
  if (max_freq + tot_nonlast_freq < 1.0 - kSmallEpsilon) {
    return cur_allele_ct_m1;
  }
  return maj_allele_idx;
}

// forced SID '0' if sids == nullptr
// ok for sample_augid_map_ptr == nullptr
PglErr AugidInitAlloc(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, uint32_t alloc_at_end, uint32_t** sample_augid_map_ptr, char** sample_augids_ptr, uintptr_t* max_sample_augid_blen_ptr) {
  const char* sample_ids = siip->sample_ids;
  const char* sids = siip->sids;
  const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
  uintptr_t max_sid_blen = siip->max_sid_blen;
  if (!sids) {
    max_sid_blen = 2;
  }
  const uintptr_t max_sample_augid_blen = max_sample_id_blen + max_sid_blen;
  *max_sample_augid_blen_ptr = max_sample_augid_blen;
  uint32_t* sample_augid_map = nullptr;
  if (sample_augid_map_ptr) {
    if (!alloc_at_end) {
      if (unlikely(bigstack_alloc_u32(sample_ct, sample_augid_map_ptr))) {
        return kPglRetNomem;
      }
    } else {
      if (unlikely(bigstack_end_alloc_u32(sample_ct, sample_augid_map_ptr))) {
        return kPglRetNomem;
      }
    }
    sample_augid_map = *sample_augid_map_ptr;
  }
  if (!alloc_at_end) {
    if (unlikely(bigstack_alloc_c(max_sample_augid_blen * sample_ct, sample_augids_ptr))) {
      return kPglRetNomem;
    }
  } else {
    if (unlikely(bigstack_end_alloc_c(max_sample_augid_blen * sample_ct, sample_augids_ptr))) {
      return kPglRetNomem;
    }
  }
  char* sample_augids_iter = *sample_augids_ptr;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    char* write_iter = strcpyax(sample_augids_iter, &(sample_ids[sample_uidx * max_sample_id_blen]), '\t');
    if (sids) {
      strcpy(write_iter, &(sids[sample_uidx * max_sid_blen]));
    } else {
      strcpy_k(write_iter, "0");
    }
    sample_augids_iter = &(sample_augids_iter[max_sample_augid_blen]);
    if (sample_augid_map) {
      sample_augid_map[sample_idx] = sample_uidx;
    }
  }
  return kPglRetSuccess;
}

PglErr SortedXidboxInitAlloc(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, XidMode xid_mode, uint32_t use_nsort, char** sorted_xidbox_ptr, uint32_t** xid_map_ptr, uintptr_t* max_xid_blen_ptr) {
  if (!(xid_mode & kfXidModeFlagSid)) {
    // two fields
    *max_xid_blen_ptr = siip->max_sample_id_blen;
    return CopySortStrboxSubset(sample_include, siip->sample_ids, sample_ct, siip->max_sample_id_blen, 0, use_nsort, sorted_xidbox_ptr, xid_map_ptr);
  }
  // three fields
  if (unlikely(AugidInitAlloc(sample_include, siip, sample_ct, 0, xid_map_ptr, sorted_xidbox_ptr, max_xid_blen_ptr))) {
    return kPglRetNomem;
  }
  if (unlikely(SortStrboxIndexed(sample_ct, *max_xid_blen_ptr, use_nsort, *sorted_xidbox_ptr, *xid_map_ptr))) {
    return kPglRetNomem;
  }
  return kPglRetSuccess;
}

PglErr SortedXidboxInitAllocEnd(const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t sample_ct, XidMode xid_mode, uint32_t use_nsort, char** sorted_xidbox_ptr, uint32_t** xid_map_ptr, uintptr_t* max_xid_blen_ptr) {
  if (!(xid_mode & kfXidModeFlagSid)) {
    // two fields
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    *max_xid_blen_ptr = max_sample_id_blen;
    if (unlikely(bigstack_end_alloc_u32(sample_ct, xid_map_ptr) ||
                 bigstack_end_alloc_c(sample_ct * max_sample_id_blen, sorted_xidbox_ptr))) {
      return kPglRetNomem;
    }
    return CopySortStrboxSubsetNoalloc(sample_include, siip->sample_ids, sample_ct, max_sample_id_blen, 0, use_nsort, *sorted_xidbox_ptr, *xid_map_ptr);
  }
  // three fields
  if (unlikely(AugidInitAlloc(sample_include, siip, sample_ct, 1, xid_map_ptr, sorted_xidbox_ptr, max_xid_blen_ptr))) {
    return kPglRetNomem;
  }
  if (unlikely(SortStrboxIndexed(sample_ct, *max_xid_blen_ptr, use_nsort, *sorted_xidbox_ptr, *xid_map_ptr))) {
    return kPglRetNomem;
  }
  return kPglRetSuccess;
}

uint32_t XidRead(uintptr_t max_xid_blen, uint32_t comma_delim, XidMode xid_mode, const char** read_pp, char* __restrict idbuf) {
  // Saves normalized tokens pointed to by *read_pp to idbuf, and advances
  // *read_pp.
  //
  // Input *read_pp must point to beginning of sample ID; this is a change from
  // plink 1.9.
  //
  // *read_pp is now set to point to the end of the last parsed token instead
  // of the beginning of the next; this is another change from plink 1.9.
  //
  // Returns 0 on missing token *or* guaranteed sample ID mismatch (longer than
  // max_xid_blen); cases can be distinguished by checking whether *read_pp ==
  // nullptr.  Otherwise, returns slen of the parsed ID saved to idbuf.  (idbuf
  // is not null-terminated.)
  //
  // Alpha 2 behavior change: If there's no FID column, it's now set to 0
  // instead of the IID value.
  const char* first_token_start = *read_pp;
  uintptr_t slen_fid = 0;  // now zero if no FID column
  uintptr_t blen_sid = 0;
  const char* sid_ptr = nullptr;
  const char* token_iter;
  const char* iid_ptr;
  uintptr_t slen_iid;
  if (comma_delim) {
    token_iter = first_token_start;
    unsigned char ucc = *token_iter;
    while (ucc != ',') {
      if (ucc < 32) {
        if (unlikely(!(xid_mode & kfXidModeFlagOneCoreTokenOk))) {
          *read_pp = nullptr;
          return 0;
        }
        goto XidRead_comma_single_core_token;
      }
      ucc = *(++token_iter);
    }
    if (xid_mode & kfXidModeFlagNeverFid) {
    XidRead_comma_single_core_token:
      iid_ptr = first_token_start;
      slen_iid = token_iter - first_token_start;
    } else {
      slen_fid = token_iter - first_token_start;
      do {
        ucc = *(++token_iter);
      } while ((ucc == ' ') || (ucc == '\t'));
      iid_ptr = token_iter;
      while ((ucc >= 32) && (ucc != ',')) {
        ucc = *(++token_iter);
      }
      slen_iid = token_iter - iid_ptr;
    }
    // token_iter now points to comma/eoln at end of IID
    if (xid_mode & (kfXidModeFlagSid | kfXidModeFlagSkipSid)) {
      if (*token_iter != ',') {
        *read_pp = nullptr;
        return 0;
      }
      do {
        ucc = *(++token_iter);
      } while ((ucc == ' ') || (ucc == '\t'));
      sid_ptr = token_iter;
      while ((ucc >= 32) && (ucc != ',')) {
        ucc = *(++token_iter);
      }
      if (xid_mode & kfXidModeFlagSid) {
        blen_sid = 1 + S_CAST(uintptr_t, token_iter - sid_ptr);
        if (token_iter == sid_ptr) {
          // special case: treat missing SID as '0'
          blen_sid = 2;
          sid_ptr = &(g_one_char_strs[96]);
        }
      }
    }
  } else {
    assert(!IsEolnKns(*first_token_start));
    token_iter = CurTokenEnd(first_token_start);
    if (xid_mode & kfXidModeFlagNeverFid) {
    XidRead_space_single_core_token:
      iid_ptr = first_token_start;
      slen_iid = token_iter - first_token_start;
    } else {
      slen_fid = token_iter - first_token_start;
      token_iter = FirstNonTspace(token_iter);
      if (IsEolnKns(*token_iter)) {
        if (unlikely(!(xid_mode & kfXidModeFlagOneCoreTokenOk))) {
          *read_pp = nullptr;
          return 0;
        }
        // need to backtrack
        token_iter = &(first_token_start[slen_fid]);
        // bugfix (26 Feb 2018): forgot to zero this
        slen_fid = 0;
        goto XidRead_space_single_core_token;
      }
      iid_ptr = token_iter;
      token_iter = CurTokenEnd(token_iter);
      slen_iid = token_iter - iid_ptr;
    }
    // token_iter now points to space/eoln at end of IID
    if (xid_mode & (kfXidModeFlagSid | kfXidModeFlagSkipSid)) {
      token_iter = FirstNonTspace(token_iter);
      if (unlikely(IsEolnKns(*token_iter))) {
        *read_pp = nullptr;
        return 0;
      }
      sid_ptr = token_iter;
      token_iter = CurTokenEnd(token_iter);
      if (xid_mode & kfXidModeFlagSid) {
        blen_sid = 1 + S_CAST(uintptr_t, token_iter - sid_ptr);
      }
    }
  }
  *read_pp = token_iter;
  const uintptr_t slen_final = slen_fid + (slen_fid == 0) + slen_iid + blen_sid + 1;
  if (slen_final >= max_xid_blen) {
    // avoid buffer overflow
    return 0;
  }
  char* idbuf_iter = idbuf;
  if (slen_fid) {
    idbuf_iter = memcpya(idbuf, first_token_start, slen_fid);
  } else {
    *idbuf_iter++ = '0';
  }
  *idbuf_iter++ = '\t';
  idbuf_iter = memcpya(idbuf_iter, iid_ptr, slen_iid);
  if (blen_sid) {
    *idbuf_iter++ = '\t';
    memcpy(idbuf_iter, sid_ptr, blen_sid - 1);
  }
  return slen_final;
}

BoolErr SortedXidboxReadMultifind(const char* __restrict sorted_xidbox, uintptr_t max_xid_blen, uintptr_t xid_ct, uint32_t comma_delim, XidMode xid_mode, const char** read_pp, uint32_t* __restrict xid_idx_start_ptr, uint32_t* __restrict xid_idx_end_ptr, char* __restrict idbuf) {
  // sorted_xidbox = packed, sorted list of ID strings to search over.
  const uint32_t slen_final = XidRead(max_xid_blen, comma_delim, xid_mode, read_pp, idbuf);
  if (!slen_final) {
    return 1;
  }
  const uint32_t lb_idx = bsearch_strbox_lb(idbuf, sorted_xidbox, slen_final, max_xid_blen, xid_ct);
  idbuf[slen_final] = ' ';
  const uint32_t ub_idx = ExpsearchStrLb(idbuf, sorted_xidbox, slen_final + 1, max_xid_blen, xid_ct, lb_idx);
  if (lb_idx == ub_idx) {
    return 1;
  }
  *xid_idx_start_ptr = lb_idx;
  *xid_idx_end_ptr = ub_idx;
  return 0;
}

PglErr LoadXidHeader(const char* flag_name, XidHeaderFlags xid_header_flags, uintptr_t* line_idx_ptr, TextStream* txsp, XidMode* xid_mode_ptr, char** line_startp) {
  // possible todo: support comma delimiter
  uintptr_t line_idx = *line_idx_ptr;
  char* line_iter;
  uint32_t is_header_line;
  do {
    ++line_idx;
    line_iter = TextGet(txsp);
    // eof may be ok
    if (!line_iter) {
      return TextStreamRawErrcode(txsp);
    }
    is_header_line = (line_iter[0] == '#');
  } while (is_header_line && (!tokequal_k(&(line_iter[1]), "FID")) && (!tokequal_k(&(line_iter[1]), "IID")));
  if (line_startp) {
    *line_startp = line_iter;
  }
  XidMode xid_mode = kfXidMode0;
  if (is_header_line) {
    // The following header leading columns are supported:
    // #FID IID
    // #FID IID SID (SID ignored if kfXidHeaderIgnoreSid set)
    // #IID
    // #IID SID
    char* linebuf_iter = &(line_iter[4]);
    if (line_iter[1] == 'I') {
      xid_mode = kfXidModeFlagNeverFid;
    } else {
      linebuf_iter = FirstNonTspace(linebuf_iter);
      if (unlikely(!tokequal_k(linebuf_iter, "IID"))) {
        logerrprintf("Error: No IID column on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
        return kPglRetMalformedInput;
      }
      linebuf_iter = &(linebuf_iter[3]);
    }
    linebuf_iter = FirstNonTspace(linebuf_iter);
    if (tokequal_k(linebuf_iter, "SID")) {
      if (!(xid_header_flags & kfXidHeaderIgnoreSid)) {
        xid_mode |= kfXidModeFlagSid;
      } else {
        xid_mode |= kfXidModeFlagSkipSid;
      }
      linebuf_iter = FirstNonTspace(&(linebuf_iter[3]));
    }
    if ((xid_mode & kfXidModeCoreMask) == kfXidModeFlagNeverFid) {
      xid_mode |= kfXidModeFlagOneCoreTokenOk;
    }
    line_iter = linebuf_iter;
  } else {
    xid_mode = (xid_header_flags & kfXidHeaderFixedWidth)? kfXidModeFidIid : kfXidModeFidIidOrIid;
  }
  *line_idx_ptr = line_idx;
  *xid_mode_ptr = xid_mode;
  return kPglRetSuccess;
}

PglErr OpenAndLoadXidHeader(const char* fname, const char* flag_name, XidHeaderFlags xid_header_flags, uint32_t max_line_blen, TextStream* txsp, XidMode* xid_mode_ptr, uintptr_t* line_idx_ptr, char** line_startp) {
  PglErr reterr = InitTextStream(fname, max_line_blen, 1, txsp);
  // remove unlikely() if any caller treats eof as ok
  if (unlikely(reterr)) {
    return reterr;
  }
  *line_idx_ptr = 0;
  return LoadXidHeader(flag_name, xid_header_flags, line_idx_ptr, txsp, xid_mode_ptr, line_startp);
}

PglErr LoadXidHeaderPair(const char* flag_name, uint32_t sid_over_fid, uintptr_t* line_idx_ptr, TextStream* txsp, XidMode* xid_mode_ptr, char** line_startp, char** line_iterp) {
  uintptr_t line_idx = *line_idx_ptr;
  char* line_iter;
  uint32_t is_header_line;
  do {
    ++line_idx;
    line_iter = TextGet(txsp);
    // eof may be ok
    if (!line_iter) {
      return TextStreamRawErrcode(txsp);
    }
    is_header_line = (line_iter[0] == '#');
  } while (is_header_line && (!tokequal_k(&(line_iter[1]), "FID1")) && (!tokequal_k(&(line_iter[1]), "ID1")) && (!tokequal_k(&(line_iter[1]), "IID1")));
  if (line_startp) {
    *line_startp = line_iter;
  }
  XidMode xid_mode = kfXidMode0;
  if (is_header_line) {
    // The following header leading columns are supported:
    // #FID1 {I}ID1 FID2 {I}ID2
    // #FID1 {I}ID1 SID1 FID2 {I}ID2 SID2
    // #{I}ID1 {I}ID2
    // #{I}ID1 SID1 {I}ID2 SID2
    char* linebuf_iter = &(line_iter[4]);
    if (line_iter[1] == 'I') {
      xid_mode = kfXidModeFlagNeverFid;
    } else {
      linebuf_iter = FirstNonTspace(linebuf_iter);
      if (tokequal_k(linebuf_iter, "ID1")) {
        linebuf_iter = &(linebuf_iter[3]);
      } else if (likely(tokequal_k(linebuf_iter, "IID1"))) {
        linebuf_iter = &(linebuf_iter[4]);
      } else {
        logerrprintf("Error: No [I]ID1 column on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
        return kPglRetMalformedInput;
      }
    }
    linebuf_iter = FirstNonTspace(linebuf_iter);
    if (tokequal_k(linebuf_iter, "SID1")) {
      xid_mode |= kfXidModeFlagSid;
      linebuf_iter = FirstNonTspace(&(linebuf_iter[4]));
    }
    // Now verify the expected ID2 token sequence follows.
    if (!(xid_mode & kfXidModeFlagNeverFid)) {
      if (unlikely(!tokequal_k(linebuf_iter, "FID2"))) {
        logerrprintf("Error: No FID2 column on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
        return kPglRetMalformedInput;
      }
      linebuf_iter = FirstNonTspace(&(linebuf_iter[4]));
    }
    if (tokequal_k(linebuf_iter, "ID2")) {
      linebuf_iter = &(linebuf_iter[3]);
    } else if (likely(tokequal_k(linebuf_iter, "IID2"))) {
      linebuf_iter = &(linebuf_iter[4]);
    } else {
      logerrprintf("Error: No [I]ID2 column on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
      return kPglRetMalformedInput;
    }
    if (xid_mode & kfXidModeFlagSid) {
      linebuf_iter = FirstNonTspace(linebuf_iter);
      if (unlikely(!tokequal_k(linebuf_iter, "SID2"))) {
        logerrprintf("Error: No SID2 column on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
        return kPglRetMalformedInput;
      }
      linebuf_iter = &(linebuf_iter[4]);
    }
    line_iter = linebuf_iter;
  } else {
    // Assume the file contains nothing but ID pairs.
    const uint32_t token_ct = CountTokens(line_iter);
    if (token_ct == 2) {
      xid_mode = kfXidModeIid;
    } else if (token_ct == 4) {
      xid_mode = sid_over_fid? kfXidModeIidSid : kfXidModeFidIid;
    } else if (likely(token_ct == 6)) {
      xid_mode = kfXidModeFidIidSid;
    } else {
      logerrprintf("Error: Unexpected token count on line %" PRIuPTR " of --%s file.\n", line_idx, flag_name);
      return kPglRetMalformedInput;
    }
  }
  if (line_iter) {
    *line_iterp = line_iter;
  }
  *line_idx_ptr = line_idx;
  *xid_mode_ptr = xid_mode;
  return kPglRetSuccess;
}

void MarkParents(const ParentalIdInfo* parental_id_infop, const char* sorted_xidbox, const uint32_t* xid_map, uint32_t sample_ct, uintptr_t max_xid_blen, uint32_t use_nsort, uintptr_t* dst, char* idbuf) {
  const char* paternal_ids = parental_id_infop->paternal_ids;
  const char* maternal_ids = parental_id_infop->maternal_ids;
  const uintptr_t max_paternal_id_blen = parental_id_infop->max_paternal_id_blen;
  const uintptr_t max_maternal_id_blen = parental_id_infop->max_maternal_id_blen;
  for (uint32_t xid_idx_start = 0; xid_idx_start != sample_ct; ) {
    const char* sorted_fam_xidbox = &(sorted_xidbox[xid_idx_start * max_xid_blen]);
    const char* fid_start = sorted_fam_xidbox;
    const char* fid_end = AdvToDelim(fid_start, '\t');
    const uint32_t fid_slen = fid_end - fid_start;
    memcpy(idbuf, fid_start, fid_slen);
    idbuf[fid_slen] = ' ';
    const uint32_t fid_blen = fid_slen + 1;
    idbuf[fid_blen] = '\0';
    uint32_t fam_xid_ct;
    if (use_nsort) {
      fam_xid_ct = ExpsearchNsortStrLb(idbuf, sorted_fam_xidbox, max_xid_blen, sample_ct - xid_idx_start, 1);
    } else {
      fam_xid_ct = ExpsearchStrLb(idbuf, sorted_fam_xidbox, fid_blen, max_xid_blen, sample_ct - xid_idx_start, 1);
    }
    idbuf[fid_slen] = '\t';
    const uint32_t* fam_xid_map = &(xid_map[xid_idx_start]);
    char* iid_write_start = &(idbuf[fid_blen]);
    for (uint32_t fam_xid_idx = 0; fam_xid_idx != fam_xid_ct; ++fam_xid_idx) {
      const uint32_t sample_uidx = fam_xid_map[fam_xid_idx];
      for (uint32_t parent_idx = 0; parent_idx != 2; ++parent_idx) {
        const char* parent_id;
        if (parent_idx == 0) {
          parent_id = &(paternal_ids[sample_uidx * max_paternal_id_blen]);
        } else {
          parent_id = &(maternal_ids[sample_uidx * max_maternal_id_blen]);
        }
        const uint32_t parent_id_slen = strlen(parent_id);
        const uintptr_t xid_slen = fid_blen + parent_id_slen;
        if (((parent_id_slen == 1) && (parent_id[0] == '0')) || (xid_slen >= max_xid_blen)) {
          continue;
        }
        memcpy(iid_write_start, parent_id, parent_id_slen + 1);
        uint32_t parent_fxid_idx_start;
        if (use_nsort) {
          parent_fxid_idx_start = bsearch_strbox_lb_natural(idbuf, sorted_fam_xidbox, xid_slen, max_xid_blen, fam_xid_ct);
        } else {
          parent_fxid_idx_start = bsearch_strbox_lb(idbuf, sorted_fam_xidbox, xid_slen, max_xid_blen, fam_xid_ct);
        }
        for (uint32_t uii = parent_fxid_idx_start; uii != fam_xid_ct; ++uii) {
          const char* cur_xid = &(sorted_fam_xidbox[uii * max_xid_blen]);
          if ((!memequal(idbuf, cur_xid, xid_slen)) || (cur_xid[xid_slen] > ' ')) {
            break;
          }
          SetBit(xid_map[uii], dst);
        }
      }
    }
    xid_idx_start += fam_xid_ct;
  }
}

// Assumes no duplicates.
void InitXidHtable(const SampleIdInfo* siip, uint32_t sample_ct, uint32_t xid_htable_size, uint32_t* xid_htable, char* idbuf) {
  SetAllU32Arr(xid_htable_size, xid_htable);
  const char* sample_ids_iter = siip->sample_ids;
  const char* sids_iter = siip->sids;
  const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
  const uintptr_t max_sid_blen = siip->max_sid_blen;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    uint32_t slen = strlen(sample_ids_iter);
    const char* cur_sample_id;
    if (!sids_iter) {
      cur_sample_id = sample_ids_iter;
    } else {
      char* id_iter = memcpyax(idbuf, sample_ids_iter, slen, '\t');
      const uint32_t sid_blen = 1 + strlen(sids_iter);
      memcpy(id_iter, sids_iter, sid_blen);
      slen += sid_blen;
      sids_iter = &(sids_iter[max_sid_blen]);
      cur_sample_id = idbuf;
    }
    HtableAddNondup(cur_sample_id, slen, xid_htable_size, sample_idx, xid_htable);
    sample_ids_iter = &(sample_ids_iter[max_sample_id_blen]);
  }
}

BoolErr LookupXidHtable(const char* line_start, const SampleIdInfo* siip, const uint32_t* xid_htable, uint32_t xid_htable_size, uint32_t fid_present, uint32_t sid_present, uint32_t* sample_idxp, char* idbuf) {
  const char* read_iid_start = line_start;
  const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
  char* write_id_iter = idbuf;
  if (fid_present) {
    const char* fid_end = CurTokenEnd(line_start);
    const uint32_t fid_slen = fid_end - line_start;
    // check this before copying, to prevent idbuf overflow
    if (fid_slen + 2 >= max_sample_id_blen) {
      *sample_idxp = UINT32_MAX;
      return 0;
    }
    write_id_iter = memcpya(write_id_iter, line_start, fid_slen);
    read_iid_start = FirstNonTspace(fid_end);
    if (unlikely(IsEolnKns(*read_iid_start))) {
      return 1;
    }
  } else {
    *write_id_iter++ = '0';
  }
  *write_id_iter++ = '\t';
  const char* read_iid_end = CurTokenEnd(read_iid_start);
  const uint32_t iid_slen = read_iid_end - read_iid_start;
  const uint32_t fid_iid_slen = iid_slen + S_CAST(uintptr_t, write_id_iter - idbuf);
  if (fid_iid_slen >= max_sample_id_blen) {
    *sample_idxp = UINT32_MAX;
    return 0;
  }
  write_id_iter = memcpya(write_id_iter, read_iid_start, iid_slen);

  const char* sample_ids = siip->sample_ids;
  const char* sids = siip->sids;
  if (!sids) {
    *write_id_iter = '\0';
    *sample_idxp = StrboxHtableFind(idbuf, sample_ids, xid_htable, max_sample_id_blen, fid_iid_slen, xid_htable_size);
    return 0;
  }

  *write_id_iter++ = '\t';
  char* write_sid_start = write_id_iter;
  uint32_t sid_slen = 1;
  const uintptr_t max_sid_blen = siip->max_sid_blen;
  if (sid_present) {
    const char* read_sid_start = FirstNonTspace(read_iid_end);
    if (unlikely(IsEolnKns(*read_sid_start))) {
      return 1;
    }
    const char* read_sid_end = CurTokenEnd(read_sid_start);
    sid_slen = read_sid_end - read_sid_start;
    if (sid_slen >= max_sid_blen) {
      *sample_idxp = UINT32_MAX;
      return 0;
    }
    write_id_iter = memcpya(write_id_iter, read_sid_start, sid_slen);
  } else {
    *write_id_iter++ = '0';
  }
  *write_id_iter = '\0';
  const uint32_t fid_iid_blen = fid_iid_slen + 1;
  const uint32_t sid_blen = sid_slen + 1;
  uint32_t hashval = Hashceil(idbuf, fid_iid_slen + sid_blen, xid_htable_size);
  write_sid_start[-1] = '\0';
  while (1) {
    const uint32_t sample_idx = xid_htable[hashval];
    if ((sample_idx == UINT32_MAX) ||
        (memequal(idbuf, &(sample_ids[sample_idx * max_sample_id_blen]), fid_iid_blen) &&
         memequal(write_sid_start, &(sids[sample_idx * max_sid_blen]), sid_blen))) {
      *sample_idxp = sample_idx;
      return 0;
    }
    if (++hashval == xid_htable_size) {
      hashval = 0;
    }
  }
}

PglErr CheckXidUniqueness(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* err_suffix_str, uint32_t sample_ct) {
  // PopulateStrboxHtable() isn't sufficient for this when SIDs are present.
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t xid_htable_size = GetHtableFastSize(sample_ct);
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    uint32_t* xid_htable;
    char* idbuf;
    if (unlikely(bigstack_alloc_u32(xid_htable_size, &xid_htable) ||
                 bigstack_alloc_c(max_sample_id_blen + max_sid_blen, &idbuf))) {
      goto CheckXidUniqueness_ret_NOMEM;
    }
    // Lots of overlap with InitXidHtable().
    SetAllU32Arr(xid_htable_size, xid_htable);
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
      const uint32_t sample_id_slen = strlen(cur_sample_id);
      uint32_t xid_slen = sample_id_slen;
      const char* cur_xid;
      if (!sids) {
        cur_xid = cur_sample_id;
      } else {
        char* id_iter = memcpyax(idbuf, cur_sample_id, sample_id_slen, '\t');
        const char* cur_sid = &(sids[sample_uidx * max_sid_blen]);
        const uint32_t sid_blen = 1 + strlen(cur_sid);
        memcpy(id_iter, cur_sid, sid_blen);
        xid_slen += sid_blen;
        cur_xid = idbuf;
      }
      for (uint32_t hashval = Hashceil(cur_xid, xid_slen, xid_htable_size); ; ) {
        const uint32_t cur_htable_entry = xid_htable[hashval];
        if (cur_htable_entry == UINT32_MAX) {
          xid_htable[hashval] = sample_uidx;
          break;
        }
        if (memequal(cur_sample_id, &(sample_ids[cur_htable_entry * max_sample_id_blen]), sample_id_slen + 1)) {
          if ((!sids) || memequal(&(cur_xid[sample_id_slen + 1]), &(sids[cur_htable_entry * max_sid_blen]), xid_slen - sample_id_slen)) {
            char* write_iter = strcpya_k(g_logbuf, "Error: Duplicate sample ID \"");
            const char* fid_end = AdvToDelim(cur_sample_id, '\t');
            const char* iid_start = &(fid_end[1]);
            write_iter = memcpyax(write_iter, cur_sample_id, fid_end - cur_sample_id, ' ');
            write_iter = strcpya(write_iter, iid_start);
            if (sids) {
              *write_iter++ = ' ';
              write_iter = memcpya(write_iter, &(cur_xid[sample_id_slen + 1]), xid_slen - sample_id_slen - 1);
            }
            *write_iter++ = '"';
            if (err_suffix_str) {
              write_iter = strcpya(write_iter, err_suffix_str);
            }
            strcpy_k(write_iter, ".\n");
            WordWrapB(0);
            logerrputsb();
            goto CheckXidUniqueness_ret_MALFORMED_INPUT;
          }
        }
        if (++hashval == xid_htable_size) {
          hashval = 0;
        }
      }
    }
  }
  while (0) {
  CheckXidUniqueness_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CheckXidUniqueness_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr LoadSampleIds(const char* fnames, const uintptr_t* sample_include, const SampleIdInfo* siip, const char* flag_name, uint32_t raw_sample_ct, uint32_t sample_ct, LoadSampleIdsFlags flags, uintptr_t** loaded_bitarrp, uint32_t* duplicate_ctp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  const char* fname_iter = fnames;
  const char* fname_txs = nullptr;
  TextStream txs;
  PreinitTextStream(&txs);
  uintptr_t line_idx;
  {
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* loaded_bitarr;
    if (unlikely(bigstack_calloc_w(raw_sample_ctl, &loaded_bitarr))) {
      goto LoadSampleIds_ret_NOMEM;
    }
    unsigned char* bigstack_mark2 = g_bigstack_base;

    const uint32_t families_only = flags & kfLoadSampleIdsFamOnly;
    const char* sample_ids = siip->sample_ids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    char* idbuf = nullptr;
    uint32_t* xid_map = nullptr;
    char* sorted_xidbox = nullptr;
    uintptr_t max_xid_blen = max_sample_id_blen - 1;
    if (families_only) {
      // only need to do this once
      if (unlikely(bigstack_alloc_u32(sample_ct, &xid_map) ||
                   bigstack_alloc_c(sample_ct * max_xid_blen, &sorted_xidbox))) {
        goto LoadSampleIds_ret_NOMEM;
      }
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
        const char* fidt_ptr = &(sample_ids[sample_uidx * max_sample_id_blen]);
        const char* fidt_end = AdvPastDelim(fidt_ptr, '\t');
        const uint32_t cur_fidt_slen = fidt_end - fidt_ptr;
        // include trailing tab, to simplify bsearch_strbox_lb() usage
        memcpyx(&(sorted_xidbox[sample_idx * max_xid_blen]), fidt_ptr, cur_fidt_slen, '\0');
        xid_map[sample_idx] = sample_uidx;
      }
      if (unlikely(SortStrboxIndexed(sample_ct, max_xid_blen, 0, sorted_xidbox, xid_map))) {
        goto LoadSampleIds_ret_NOMEM;
      }
    }

    unsigned char* bigstack_mark3 = nullptr;
    const char* fnames_iter = fnames;
    uintptr_t duplicate_ct = 0;
    do {
      if (!bigstack_mark3) {
        fname_txs = fnames_iter;
        reterr = SizeAndInitTextStream(fnames_iter, bigstack_left() - (bigstack_left() / 4), 1, &txs);
        if (unlikely(reterr)) {
          goto LoadSampleIds_ret_TSTREAM_FAIL;
        }
        bigstack_mark3 = g_bigstack_base;
      } else {
        reterr = TextRetarget(fnames_iter, &txs);
        if (unlikely(reterr)) {
          goto LoadSampleIds_ret_TSTREAM_FAIL;
        }
        fname_txs = fnames_iter;
      }
      char* line_start;
      XidMode xid_mode;
      line_idx = 0;
      uint32_t skip_header = 0;
      if (families_only) {
        skip_header = 1;
      } else {
        reterr = LoadXidHeader(flag_name, (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeader0 : kfXidHeaderIgnoreSid, &line_idx, &txs, &xid_mode, &line_start);
        if (reterr) {
          if (likely(reterr == kPglRetEof)) {
            reterr = kPglRetSuccess;
            goto LoadSampleIds_empty_file;
          }
          goto LoadSampleIds_ret_TSTREAM_XID_FAIL;
        }
        reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
        if (unlikely(reterr)) {
          goto LoadSampleIds_ret_1;
        }
        if (unlikely(bigstack_alloc_c(max_xid_blen, &idbuf))) {
          goto LoadSampleIds_ret_NOMEM;
        }
        if (*line_start == '#') {
          skip_header = 1;
        }
      }
      if (skip_header) {
        ++line_idx;
        line_start = TextGet(&txs);
      }
      for (; line_start; ++line_idx, line_start = TextGet(&txs)) {
        if (!families_only) {
          const char* linebuf_iter = line_start;
          uint32_t xid_idx_start;
          uint32_t xid_idx_end;
          if (!SortedXidboxReadMultifind(sorted_xidbox, max_xid_blen, sample_ct, 0, xid_mode, &linebuf_iter, &xid_idx_start, &xid_idx_end, idbuf)) {
            uint32_t sample_uidx = xid_map[xid_idx_start];
            if (IsSet(loaded_bitarr, sample_uidx)) {
              ++duplicate_ct;
            } else {
              for (uint32_t xid_idx = xid_idx_start; ; ) {
                SetBit(sample_uidx, loaded_bitarr);
                if (++xid_idx == xid_idx_end) {
                  break;
                }
                sample_uidx = xid_map[xid_idx];
              }
            }
          } else if (unlikely(!linebuf_iter)) {
            goto LoadSampleIds_ret_MISSING_TOKENS;
          }
        } else {
          char* token_end = CurTokenEnd(line_start);
          // bugfix (28 Oct 2018): \n was being clobbered and not replaced
          // const char orig_token_end_char = *token_end;
          *token_end = '\t';
          const uint32_t slen = 1 + S_CAST(uintptr_t, token_end - line_start);
          uint32_t lb_idx = bsearch_strbox_lb(line_start, sorted_xidbox, slen, max_xid_blen, sample_ct);
          *token_end = ' ';
          const uint32_t ub_idx = bsearch_strbox_lb(line_start, sorted_xidbox, slen, max_xid_blen, sample_ct);
          if (ub_idx != lb_idx) {
            uint32_t sample_uidx = xid_map[lb_idx];
            if (IsSet(loaded_bitarr, sample_uidx)) {
              ++duplicate_ct;
            } else {
              for (uint32_t xid_map_idx = lb_idx; ; ) {
                SetBit(sample_uidx, loaded_bitarr);
                if (++xid_map_idx == ub_idx) {
                  break;
                }
                sample_uidx = xid_map[xid_map_idx];
              }
            }
          }
          // *token_end = orig_token_end_char;
        }
      }
      if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
        goto LoadSampleIds_ret_TSTREAM_FAIL;
      }
    LoadSampleIds_empty_file:
      if (!(flags & kfLoadSampleIdsMultifile)) {
        break;
      }
      BigstackReset(bigstack_mark3);
      fnames_iter = strnul(fnames_iter);
      ++fnames_iter;
    } while (*fnames_iter);
    reterr = kPglRetSuccess;
    *loaded_bitarrp = loaded_bitarr;
    if (duplicate_ctp) {
      *duplicate_ctp = duplicate_ct;
    }
    bigstack_mark = bigstack_mark2;
  }
  while (0) {
  LoadSampleIds_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadSampleIds_ret_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&txs)) {
      break;
    }
  LoadSampleIds_ret_TSTREAM_FAIL:
    TextStreamErrPrint(fname_txs, &txs);
    break;
  LoadSampleIds_ret_MISSING_TOKENS:
    logerrprintf("Error: --%s: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", flag_name, line_idx, fname_iter);
    reterr = kPglRetMalformedInput;
    break;
  }
 LoadSampleIds_ret_1:
  if (fname_txs) {
    CleanupTextStream2(fname_txs, &txs, &reterr);
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

char* AppendSpacedXid(const char* sample_ids, const char* sids, uint32_t write_fid, uint32_t write_sid, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t sample_uidx, char* write_iter) {
  const char* fid = &(sample_ids[max_sample_id_blen * sample_uidx]);
  const char* fid_end = AdvToDelim(fid, '\t');
  if (write_fid) {
    write_iter = memcpyax(write_iter, fid, fid_end - fid, ' ');
  }
  const char* iid = &(fid_end[1]);
  write_iter = strcpya(write_iter, iid);
  if (write_sid) {
    *write_iter++ = ' ';
    if (sids) {
      write_iter = strcpya(write_iter, &(sids[max_sid_blen * sample_uidx]));
    } else {
      *write_iter++ = '0';
    }
  }
  return write_iter;
}

// accept 'M'/'F'/'m'/'f' since that's more readable without being any less
// efficient
const unsigned char g_char_to_sex[256] =
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const char g_xymt_log_names[kChrOffsetCt][5] = {"chrX", "chrY", "XY", "chrM", "PAR1", "PAR2"};

static_assert(!(kChrRawEnd % kBytesPerVec), "kChrRawEnd must be a multiple of kBytesPerVec.");
PglErr InitChrInfo(ChrInfo* cip) {
  // "constructor".  initializes with maximum capacity.  doesn't use bigstack.
  // chr_mask, haploid_mask: bits
  // chr_file_order, chr_idx_to_foidx: uint32s
  // chr_fo_vidx_start: uint32s, with an extra trailing element
  // nonstd_names: intptr_ts
  // nonstd_id_htable: kChrHtableSize uint32s

  // this assumes kChrRawEnd is divisible by kBytesPerVec
  const uintptr_t vecs_required = 2 * BitCtToVecCt(kChrRawEnd) + 3 * (kChrRawEnd / kInt32PerVec) + 1 + (kChrRawEnd / kWordsPerVec) + Int32CtToVecCt(kChrHtableSize);

  // needed for proper cleanup
  cip->name_ct = 0;
  cip->incl_excl_name_stack = nullptr;
  if (unlikely(vecaligned_malloc(vecs_required * kBytesPerVec, &(cip->chr_mask)))) {
    return kPglRetNomem;
  }
  uintptr_t* alloc_iter = &(cip->chr_mask[BitCtToVecCt(kChrRawEnd) * kWordsPerVec]);
  cip->haploid_mask = alloc_iter;
  alloc_iter = &(alloc_iter[BitCtToVecCt(kChrRawEnd) * kWordsPerVec]);
  cip->chr_file_order = R_CAST(uint32_t*, alloc_iter);
  alloc_iter = &(alloc_iter[(kChrRawEnd / kInt32PerVec) * kWordsPerVec]);
  cip->chr_fo_vidx_start = R_CAST(uint32_t*, alloc_iter);
  alloc_iter = &(alloc_iter[((kChrRawEnd / kInt32PerVec) + 1) * kWordsPerVec]);
  cip->chr_idx_to_foidx = R_CAST(uint32_t*, alloc_iter);
  alloc_iter = &(alloc_iter[(kChrRawEnd / kInt32PerVec) * kWordsPerVec]);
  cip->nonstd_names = R_CAST(const char**, alloc_iter);
  alloc_iter = &(alloc_iter[kChrRawEnd]);
  cip->nonstd_id_htable = R_CAST(uint32_t*, alloc_iter);
  // alloc_iter = &(alloc_iter[((kChrHtableSize + (kInt32PerVec - 1)) / kInt32PerVec) * kWordsPerVec]);
  // SetAllU32Arr(kChrHtableSize, cip->nonstd_id_htable);

  ZeroWArr(kChrMaskWords, cip->chr_mask);
  ZeroWArr(kChrExcludeWords, cip->chr_exclude);

  // This is a change from plink 1.x.  MT > M since the former matches Ensembl,
  // while the latter doesn't match any major resource.  no "chr" to reduce
  // file sizes and reduce the impact of this change.
  cip->output_encoding = kfChrOutputMT;

  cip->zero_extra_chrs = 0;
  cip->is_include_stack = 0;
  cip->chrset_source = kChrsetSourceDefault;
  cip->autosome_ct = 22;
  for (uint32_t xymt_idx = 0; xymt_idx != kChrOffsetCt; ++xymt_idx) {
    cip->xymt_codes[xymt_idx] = 23 + xymt_idx;
  }
  // Now includes MT again, as of alpha 2.
  cip->haploid_mask[0] = 0x5800000;
  ZeroWArr(kChrMaskWords - 1, &(cip->haploid_mask[1]));
  return kPglRetSuccess;
}

// explicit plink 1.07 species (now initialized by command line parser):
// human: 22, X, Y, XY, MT, PAR1, PAR2 (PAR1/PAR2 added, XY deprecated in plink
//   2.0)
// cow: 29, X, Y, MT
// dog: 38, X, Y, XY, MT, PAR1, PAR2 (updated 11 Aug 2023)
// horse: 31, X, Y
// mouse: 19, X, Y
// rice: 12
// sheep: 26, X, Y

// must be safe to call this twice.
void FinalizeChrset(LoadFilterLogFlags load_filter_log_flags, ChrInfo* cip) {
  uint32_t autosome_ct = cip->autosome_ct;
  uint32_t max_code = autosome_ct;
  for (uint32_t xymt_idx_p1 = kChrOffsetCt; xymt_idx_p1; --xymt_idx_p1) {
    if (!IsI32Neg(cip->xymt_codes[xymt_idx_p1 - 1])) {
      max_code = autosome_ct + xymt_idx_p1;
      break;
    }
  }

  // could initialize haploid_mask bits (after the first) here, instead of
  // earlier...

  cip->max_numeric_code = MINV(max_code, autosome_ct + 4);
  cip->max_code = max_code;
  uintptr_t* chr_mask = cip->chr_mask;
  uintptr_t last_chr_mask_word = chr_mask[kChrMaskWords - 1];
  STD_ARRAY_REF(uint32_t, kChrOffsetCt) xymt_codes = cip->xymt_codes;
  if (last_chr_mask_word) {
    // avoids repeating some work if this is called twice
    chr_mask[kChrMaskWords - 1] = 0;

    uint32_t xymt_include = last_chr_mask_word >> (kBitsPerWord - kChrOffsetCt);
    do {
      const uint32_t xymt_idx = ctzu32(xymt_include);
      const uint32_t cur_chr_code = xymt_codes[xymt_idx];
      if (!IsI32Neg(cur_chr_code)) {
        SetBit(cur_chr_code, chr_mask);
      }
      xymt_include &= xymt_include - 1;
    } while (xymt_include);
  } else if (AllWordsAreZero(chr_mask, kChrExcludeWords) && (!cip->is_include_stack)) {
    // init_default_chr_mask()
    SetAllBits(cip->autosome_ct + 1, chr_mask);
    for (uint32_t xymt_idx = 0; xymt_idx != kChrOffsetCt; ++xymt_idx) {
      const uint32_t cur_chr_code = cip->xymt_codes[xymt_idx];
      if (!IsI32Neg(cur_chr_code)) {
        SetBit(cur_chr_code, chr_mask);
      }
    }
  } else if (load_filter_log_flags & (kfLoadFilterLogAutosome | kfLoadFilterLogAutosomePar)) {
    FillBitsNz(1, cip->autosome_ct + 1, chr_mask);
    ClearBitsNz(cip->autosome_ct + 1, kChrExcludeWords * kBitsPerWord, chr_mask);
    if (load_filter_log_flags & kfLoadFilterLogAutosomePar) {
      uint32_t par_chr_code = cip->xymt_codes[kChrOffsetXY];
      if (!IsI32Neg(par_chr_code)) {
        SetBit(par_chr_code, chr_mask);
      }
      par_chr_code = cip->xymt_codes[kChrOffsetPAR1];
      if (!IsI32Neg(par_chr_code)) {
        SetBit(par_chr_code, chr_mask);
      }
      par_chr_code = cip->xymt_codes[kChrOffsetPAR2];
      if (!IsI32Neg(par_chr_code)) {
        SetBit(par_chr_code, chr_mask);
      }
    }
  }

  uintptr_t* chr_exclude = cip->chr_exclude;
  uintptr_t last_chr_exclude_word = chr_exclude[kChrExcludeWords - 1];
  uint32_t xymt_exclude = last_chr_exclude_word >> (kBitsPerWord - kChrOffsetCt);
  last_chr_exclude_word = bzhi(last_chr_exclude_word, kBitsPerWord - kChrOffsetCt);
  for (uint32_t widx = 0; widx != kChrExcludeWords - 1; ++widx) {
    chr_mask[widx] &= ~chr_exclude[widx];
  }
  chr_mask[kChrExcludeWords - 1] &= ~last_chr_exclude_word;
  if (xymt_exclude) {
    do {
      const uint32_t xymt_idx = ctzu32(xymt_exclude);
      const uint32_t cur_chr_code = xymt_codes[xymt_idx];
      if (!IsI32Neg(cur_chr_code)) {
        ClearBit(cur_chr_code, chr_mask);
      }
      xymt_exclude &= xymt_exclude - 1;
    } while (xymt_exclude);
  }
  SetAllU32Arr(max_code + 1, cip->chr_idx_to_foidx);
}

void ForgetExtraChrNames(uint32_t reinitialize, ChrInfo* cip) {
  const uint32_t name_ct = cip->name_ct;
  if (name_ct) {
    const char** nonstd_names = cip->nonstd_names;
    const uint32_t chr_idx_end = cip->max_code + 1 + name_ct;
    for (uint32_t chr_idx = cip->max_code + 1; chr_idx != chr_idx_end; ++chr_idx) {
      free_const(nonstd_names[chr_idx]);
      nonstd_names[chr_idx] = nullptr;
    }
    if (reinitialize) {
      // SetAllU32Arr(kChrHtableSize, cip->nonstd_id_htable);
      cip->name_ct = 0;
    }
  }
}

// not currently called.  might want to do so in the future.
PglErr FinalizeChrInfo(ChrInfo* cip) {
  const uint32_t chr_ct = cip->chr_ct;
  const uint32_t name_ct = cip->name_ct;
  const uint32_t chr_code_end = cip->max_code + 1 + name_ct;
  const uint32_t chr_code_bitvec_ct = BitCtToVecCt(chr_code_end);
  const uint32_t chr_ct_int32vec_ct = Int32CtToVecCt(chr_ct);
  const uint32_t chr_ct_p1_int32vec_ct = 1 + (chr_ct / kInt32PerVec);
  const uint32_t chr_code_end_int32vec_ct = Int32CtToVecCt(chr_code_end);
  const uint32_t chr_code_end_wordvec_ct = WordCtToVecCt(chr_code_end);
  uint32_t final_vecs_required = 2 * chr_code_bitvec_ct + chr_ct_int32vec_ct + chr_ct_p1_int32vec_ct + chr_code_end_int32vec_ct;
  if (name_ct) {
    final_vecs_required += chr_code_end_wordvec_ct + Int32CtToVecCt(kChrHtableSize);
  }
  uintptr_t* new_alloc;
  if (unlikely(vecaligned_malloc(final_vecs_required * kBytesPerVec, &new_alloc))) {
    return kPglRetNomem;
  }
  uintptr_t* old_alloc = cip->chr_mask;
  uintptr_t* new_alloc_iter = new_alloc;

  memcpy(new_alloc_iter, cip->chr_mask, chr_code_bitvec_ct * kBytesPerVec);
  cip->chr_mask = new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_code_bitvec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->haploid_mask, chr_code_bitvec_ct * kBytesPerVec);
  cip->haploid_mask = new_alloc_iter;
  new_alloc_iter = &(new_alloc_iter[chr_code_bitvec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->chr_file_order, chr_ct_int32vec_ct * kBytesPerVec);
  cip->chr_file_order = R_CAST(uint32_t*, new_alloc_iter);
  new_alloc_iter = &(new_alloc_iter[chr_ct_int32vec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->chr_fo_vidx_start, chr_ct_p1_int32vec_ct * kBytesPerVec);
  cip->chr_fo_vidx_start = R_CAST(uint32_t*, new_alloc_iter);
  new_alloc_iter = &(new_alloc_iter[chr_ct_p1_int32vec_ct * kWordsPerVec]);

  memcpy(new_alloc_iter, cip->chr_idx_to_foidx, chr_code_end_int32vec_ct * kBytesPerVec);
  cip->chr_idx_to_foidx = R_CAST(uint32_t*, new_alloc_iter);

  if (!name_ct) {
    cip->nonstd_names = nullptr;
    cip->nonstd_id_htable = nullptr;
  } else {
    new_alloc_iter = &(new_alloc_iter[chr_code_end_int32vec_ct * kWordsPerVec]);

    memcpy(new_alloc_iter, cip->nonstd_names, chr_code_end_wordvec_ct * kBytesPerVec);
    cip->nonstd_names = R_CAST(const char**, new_alloc_iter);
    new_alloc_iter = &(new_alloc_iter[chr_code_end_wordvec_ct * kWordsPerVec]);

    memcpy(new_alloc_iter, cip->nonstd_id_htable, kChrHtableSize * sizeof(int32_t));
    cip->nonstd_id_htable = R_CAST(uint32_t*, new_alloc_iter);
  }
  vecaligned_free(old_alloc);
  return kPglRetSuccess;
}

void CleanupChrInfo(ChrInfo* cip) {
  if (cip->chr_mask) {
    ForgetExtraChrNames(0, cip);
    vecaligned_free(cip->chr_mask);
    cip->chr_mask = nullptr;
  }
  LlStr* llstr_ptr = cip->incl_excl_name_stack;
  while (llstr_ptr) {
    LlStr* next_ptr = llstr_ptr->next;
    free(llstr_ptr);
    llstr_ptr = next_ptr;
  }
  cip->incl_excl_name_stack = nullptr;
}

char* ChrNameStdEx(const ChrInfo* cip, uint32_t chr_idx, ChrOutput output_encoding, char* buf) {
  if (chr_idx > cip->max_numeric_code) {
    // This is usually encoding-independent; no real numeric representation of
    // PAR1/PAR2 is defined.  However, since there'd otherwise be no way to
    // include the pseudoautosomal regions in e.g. ADMIXTURE, we render them
    // them as 25 (in the human case) when "--output-chr 26" is specified.
    if (output_encoding) {
      const uint32_t parx = 'P' + ('A' << 8) + ('R' << 16) + ('0' << 24) + ((chr_idx - cip->max_numeric_code) << 24);
      return memcpya(buf, &parx, 4);
    }
    return u32toa(cip->autosome_ct + (kChrOffsetXY + 1), buf);
  }
  if (output_encoding & (kfChrOutputPrefix | kfChrOutput0M)) {
    if (output_encoding == kfChrOutput0M) {
      // force two chars
      if (chr_idx <= cip->autosome_ct) {
        buf = memcpya_k(buf, &(kDigitPair[chr_idx]), 2);
      } else if (chr_idx == cip->xymt_codes[kChrOffsetY]) {
        buf = strcpya_k(buf, "XY");
      } else {
        *buf++ = '0';
        if (chr_idx == cip->xymt_codes[kChrOffsetX]) {
          *buf++ = 'X';
        } else {
          // assumes only X/Y/XY/MT defined
          *buf++ = (chr_idx == cip->xymt_codes[kChrOffsetY])? 'Y' : 'M';
        }
      }
      return buf;
    }
    buf = strcpya_k(buf, "chr");
  }
  if ((!(output_encoding & (kfChrOutputM | kfChrOutputMT))) || (chr_idx <= cip->autosome_ct)) {
    return u32toa(chr_idx, buf);
  }
  if (chr_idx == cip->xymt_codes[kChrOffsetX]) {
    *buf++ = 'X';
  } else if (chr_idx == cip->xymt_codes[kChrOffsetY]) {
    *buf++ = 'Y';
  } else if (chr_idx == cip->xymt_codes[kChrOffsetXY]) {
    buf = strcpya_k(buf, "XY");
  } else {
    *buf++ = 'M';
    if (output_encoding & kfChrOutputMT) {
      *buf++ = 'T';
    }
  }
  return buf;
}

char* chrtoa(const ChrInfo* cip, uint32_t chr_idx, char* buf) {
  // assumes chr_idx is valid
  if (!chr_idx) {
    // Note that this never has a 'chr' prefix.  Which is probably fine, it
    // isn't a real chromosome.
    *buf++ = '0';
    return buf;
  }
  if (chr_idx <= cip->max_code) {
    return ChrNameStd(cip, chr_idx, buf);
  }
  if (cip->zero_extra_chrs) {
    *buf++ = '0';
    return buf;
  }
  return strcpya(buf, cip->nonstd_names[chr_idx]);
}

uint32_t GetMaxChrSlen(const ChrInfo* cip) {
  // does not include trailing null
  // can be overestimate
  // if more functions start calling this, it should just be built into
  // LoadPvar() instead
  if (cip->zero_extra_chrs) {
    return 3 + kMaxChrTextnum;
  }
  const uint32_t chr_ct = cip->chr_ct;
  const uint32_t max_code = cip->max_code;
  uint32_t max_chr_slen = 3 + kMaxChrTextnum;
  for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
    if (!IsSet(cip->chr_mask, chr_idx)) {
      continue;
    }
    if (chr_idx > max_code) {
      const uint32_t name_slen = strlen(cip->nonstd_names[chr_idx]);
      if (name_slen > max_chr_slen) {
        max_chr_slen = name_slen;
      }
    }
  }
  return max_chr_slen;
}

uint32_t IsHaploidChrPresent(const ChrInfo* cip) {
  const uintptr_t* chr_mask = cip->chr_mask;
  const uintptr_t* haploid_mask = cip->haploid_mask;
  // since we don't load haploid vs. diploid info from ##contig header lines,
  // this is sufficient
  for (uint32_t widx = 0; widx != kChrExcludeWords; ++widx) {
    if (chr_mask[widx] & haploid_mask[widx]) {
      return 1;
    }
  }
  return 0;
}

uint32_t IsAutosomalDiploidChrPresent(const ChrInfo* cip) {
  const uintptr_t* haploid_mask = cip->haploid_mask;
  if (haploid_mask[0] & 1) {
    return 0;
  }
  const uint32_t max_code_p1 = cip->max_code + 1;
  const uint32_t name_ct = cip->name_ct;
  const uint32_t chr_code_end = max_code_p1 + name_ct;
  const uint32_t word_ct = BitCtToWordCt(chr_code_end);
  const uintptr_t* chr_mask = cip->chr_mask;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    if (chr_mask[widx] & (~haploid_mask[widx])) {
      return 1;
    }
  }
  return 0;
}

static inline uint32_t SingleCapLetterChrCode(uint32_t cap_letter) {
  if (cap_letter == 'X') {
    return kChrRawX;
  }
  if (cap_letter == 'Y') {
    return kChrRawY;
  }
  if (cap_letter == 'M') {
    return kChrRawMT;
  }
  return UINT32_MAX;
}

#ifndef HIGH_AUTOSOME_NUM_BUILD
static_assert(kMaxChrTextnumSlen == 2, "GetChrCodeRaw() must be updated.");
#else
static_assert(kMaxChrTextnumSlen == 3, "GetChrCodeRaw() must be updated.");
#endif
uint32_t GetChrCodeRaw(const char* str_iter) {
  // any character <= ' ' is considered a terminator
  // note that char arithmetic tends to be compiled to uint32 operations, so we
  // mostly work with ints here
  uint32_t first_char_code = ctou32(str_iter[0]);
  uint32_t first_char_toui;
  if (first_char_code < 58) {
  GetChrCodeRaw_digits:
    first_char_toui = first_char_code - '0';
    if (first_char_toui < 10) {
      const uint32_t second_char_code = ctou32(str_iter[1]);
      if (second_char_code <= ' ') {
        return first_char_toui;
      }
      const uint32_t third_char_code = ctou32(str_iter[2]);
#ifndef HIGH_AUTOSOME_NUM_BUILD
      if (third_char_code <= ' ') {
        const uint32_t second_char_toui = second_char_code - '0';
        if (second_char_toui < 10) {
          return first_char_toui * 10 + second_char_toui;
        }
        if (!first_char_toui) {
          // accept '0X', '0Y', '0M' emitted by Oxford software
          return SingleCapLetterChrCode(second_char_code & 0xdf);
        }
      } else if (((second_char_code - 48) < 10) && ((third_char_code - 48) < 10)) {
        // bugfix (23 Aug 2025): we don't want to return UINT32_MAXM1 on contig
        // names like chr1_KI270706v1_random .  We only want to do so for large
        // numbers.
        str_iter = &(str_iter[3]);
        while (1) {
          const uint32_t char_code = ctou32(*str_iter++);
          if (char_code <= ' ') {
            return UINT32_MAXM1;
          }
          if ((char_code - 48) >= 10) {
            break;
          }
        }
      }
#else
      const uint32_t second_char_toui = second_char_code - '0';
      if (second_char_toui < 10) {
        if (third_char_code <= ' ') {
          return first_char_toui * 10 + second_char_toui;
        }
        const uint32_t third_char_toui = third_char_code - '0';
        if (third_char_toui < 10) {
          const uint32_t fourth_char_code = ctou32(str_iter[3]);
          if (fourth_char_code <= ' ') {
            return first_char_toui * 100 + second_char_toui * 10 + third_char_toui;
          }
          if ((fourth_char_code - 48) < 10) {
            str_iter = &(str_iter[4]);
            while (1) {
              const uint32_t char_code = ctou32(*str_iter++);
              if (char_code <= ' ') {
                return UINT32_MAXM1;
              }
              if ((char_code - 48) >= 10) {
                break;
              }
            }
          }
        }
      } else {
        // accept '0X', '0Y', '0M' emitted by Oxford software
        if ((!first_char_toui) && (third_char_code <= ' ')) {
          return SingleCapLetterChrCode(second_char_code & 0xdf);
        }
      }
#endif
    }
    return UINT32_MAX;
  }
  first_char_code &= 0xdf;
  uint32_t second_char_code = ctou32(str_iter[1]);
  if (first_char_code == 'P') {
    // chrPAR1 *not* supported; has to be PAR1 by itself.
    // can't do uint16_t compare of multiple characters, since we could be
    // dealing with a length-1 null-terminated string; that IS faster when it's
    // safe, though
    if (((second_char_code & 0xdf) == 'A') && ((ctou32(str_iter[2]) & 0xdf) == 'R')) {
      const uint32_t par_idx_m1 = ctou32(str_iter[3]) - '1';
      if ((par_idx_m1 < 2) && (ctou32(str_iter[4]) <= ' ')) {
        return kChrRawPAR1 + par_idx_m1;
      }
    }
    return UINT32_MAX;
  }
  if (first_char_code == 'C') {
    if (((second_char_code & 0xdf) != 'H') || ((ctou32(str_iter[2]) & 0xdf) != 'R')) {
      return UINT32_MAX;
    }
    str_iter = &(str_iter[3]);
    first_char_code = ctou32(str_iter[0]);
    if (first_char_code < 58) {
      goto GetChrCodeRaw_digits;
    }
    first_char_code &= 0xdf;
    second_char_code = ctou32(str_iter[1]);
  }
  if (second_char_code <= ' ') {
    return SingleCapLetterChrCode(first_char_code);
  }
  if (ctou32(str_iter[2]) <= ' ') {
    second_char_code &= 0xdf;
    if ((first_char_code == 'X') && (second_char_code == 'Y')) {
      return kChrRawXY;
    } else if ((first_char_code == 'M') && (second_char_code == 'T')) {
      return kChrRawMT;
    }
  }
  return UINT32_MAX;
}

uint32_t GetChrCode(const char* chr_name, const ChrInfo* cip, uint32_t name_slen) {
  // requires chr_name to be null-terminated
  // in practice, name_slen will usually already be known, may as well avoid
  // redundant strlen() calls even though this uglifies the interface
  // does not perform exhaustive error-checking
  // UINT32_MAX = --allow-extra-chr ok, UINT32_MAXM1 = total fail
  const uint32_t chr_code_raw = GetChrCodeRaw(chr_name);
  const uint32_t autosome_ct = cip->autosome_ct;
  if (chr_code_raw <= autosome_ct) {
    return chr_code_raw;
  }
  if (chr_code_raw != UINT32_MAX) {
    if (chr_code_raw <= cip->max_numeric_code) {
      return cip->xymt_codes[chr_code_raw - autosome_ct - 1];
    }
    if ((chr_code_raw >= kMaxContigs) && (chr_code_raw != UINT32_MAXM1)) {
      return cip->xymt_codes[chr_code_raw - kMaxContigs];
    }
    return UINT32_MAXM1;
  }
  if (!cip->name_ct) {
    return UINT32_MAX;
  }
  // note that IdHtableFind returns UINT32_MAX if name not found
  // can't overread, nonstd_names not in main workspace
  return IdHtableFind(chr_name, cip->nonstd_names, cip->nonstd_id_htable, name_slen, kChrHtableSize);
}

uint32_t GetChrCodeCounted(const ChrInfo* cip, uint32_t name_slen, char* chr_name) {
  // When the chromosome name isn't null-terminated.
  // Yeah, probably want to revise this so that chr_name doesn't need to be
  // mutable here.  However, that currently requires new substitutes for both
  // GetChrCodeRaw AND IdHtableFind (IdHtableFindNnt doesn't work due to
  // overread); when either of those are needed for some other reason, that's a
  // good time to revise this function.
  char* s_end = &(chr_name[name_slen]);
  const char tmpc = *s_end;
  *s_end = '\0';
  const uint32_t chr_code = GetChrCode(chr_name, cip, name_slen);
  *s_end = tmpc;
  return chr_code;
}

void ChrError(const char* chr_name, const char* file_descrip, const ChrInfo* cip, uintptr_t line_idx, uint32_t error_code) {
  // assumes chr_name is null-terminated
  const uint32_t raw_code = GetChrCodeRaw(chr_name);
  logputs("\n");
  if (line_idx) {
    logerrprintfww("Error: Invalid chromosome code '%s' on line %" PRIuPTR " of %s.\n", chr_name, line_idx, file_descrip);
  } else {
    logerrprintfww("Error: Invalid chromosome code '%s' in %s.\n", chr_name, file_descrip);
  }
  if ((S_CAST(int32_t, raw_code) > S_CAST(int32_t, cip->max_numeric_code)) && ((raw_code <= kMaxChrTextnum + kChrOffsetCt) || (raw_code >= kMaxContigs))) {
    // numeric code or X/Y/MT/PAR
    if (cip->chrset_source == kChrsetSourceDefault) {
      logerrputs("(This is disallowed for humans.  Check if the problem is with your data, or if\nyou forgot to define a different chromosome set with e.g. --chr-set.)\n");
    } else if (cip->chrset_source == kChrsetSourceCmdline) {
      logerrputs("(This is disallowed by your command-line flags.)\n");
    } else {
      // kChrsetSourceFile
      logerrputs("(This is disallowed by the file's own ##chrSet header line.)\n");
    }
    // maybe want to print message(s) depending on whether chromosome set was
    // defined on the command line or by the input file?
  } else if (error_code == UINT32_MAX) {
    logerrputs("(Use --allow-extra-chr to force it to be accepted.)\n");
  }
}

PglErr TryToAddChrName(const char* chr_name, const char* file_descrip, uintptr_t line_idx, uint32_t name_slen, uint32_t prohibit_extra_chrs, uint32_t* chr_idx_ptr, ChrInfo* cip) {
  // assumes chr_name is either nonstandard (i.e. not "2", "chr2", "chrX",
  // etc.), or a rejected xymt.
  // requires chr_name to be null-terminated
  // assumes chr_idx currently has the return value of GetChrCode()
  if (unlikely(prohibit_extra_chrs || ((*chr_idx_ptr) == UINT32_MAXM1))) {
    ChrError(chr_name, file_descrip, cip, line_idx, *chr_idx_ptr);
    return kPglRetMalformedInput;
  }
  // When --allow-extra-chr was specified, we save all nonstandard but valid
  // chromosome names, even if they're all filtered out by the chromosome mask.

  // quasi-bugfix: remove redundant hash table check

  if (unlikely(chr_name[0] == '#')) {
    // redundant with some of the comment-skipping loaders, but this isn't
    // performance-critical
    logputs("\n");
    logerrputs("Error: Chromosome/contig names may not begin with '#'.\n");
    return kPglRetMalformedInput;
  }
  if (unlikely(name_slen > kMaxIdSlen)) {
    logputs("\n");
    if (line_idx) {
      logerrprintfww("Error: Line %" PRIuPTR " of %s has an excessively long chromosome/contig name. (The " PROG_NAME_STR " limit is " MAX_ID_SLEN_STR " characters.)\n", line_idx, file_descrip);
    } else {
      logerrprintfww("Error: Excessively long chromosome/contig name in %s. (The " PROG_NAME_STR " limit is " MAX_ID_SLEN_STR " characters.)\n", file_descrip);
    }
    return kPglRetMalformedInput;
  }
  const uint32_t max_code_p1 = cip->max_code + 1;
  const uint32_t name_ct = cip->name_ct;
  const uint32_t chr_code_end = max_code_p1 + name_ct;
  if (unlikely(chr_code_end == kMaxContigs)) {
    logputs("\n");
#ifdef HIGH_CONTIG_BUILD
    logerrputs("Error: Too many distinct nonstandard chromosome/contig names.\n");
#else
    logerrputs("Error: Too many distinct nonstandard chromosome/contig names.\nThe usual limit is about 65k.  You can raise this to ~980k by uncommenting\n\"#define HIGH_CONTIG_BUILD\" at the top of plink2_common.h and recompiling.  If\nthat still isn't enough, we should be able to help you if you make a post in\nthe plink2-users Google group.\n");
#endif
    return kPglRetMalformedInput;
  }
  if (!name_ct) {
    // lazy initialization
    SetAllU32Arr(kChrHtableSize, cip->nonstd_id_htable);
  }
  char* new_nonstd_name;
  if (unlikely(pgl_malloc(name_slen + 1, &new_nonstd_name))) {
    return kPglRetNomem;
  }
  LlStr* name_stack_ptr = cip->incl_excl_name_stack;
  uint32_t in_name_stack = 0;
  while (name_stack_ptr) {
    // there shouldn't be many of these, so sorting is unimportant
    if (!strcmp(chr_name, name_stack_ptr->str)) {
      in_name_stack = 1;
      break;
    }
    name_stack_ptr = name_stack_ptr->next;
  }
  if ((in_name_stack && cip->is_include_stack) || ((!in_name_stack) && (!cip->is_include_stack))) {
    SetBit(chr_code_end, cip->chr_mask);
    if (cip->haploid_mask[0] & 1) {
      SetBit(chr_code_end, cip->haploid_mask);
    }
  }
  memcpy(new_nonstd_name, chr_name, name_slen + 1);
  cip->nonstd_names[chr_code_end] = new_nonstd_name;
  *chr_idx_ptr = chr_code_end;
  cip->name_ct = name_ct + 1;
  uint32_t* id_htable = cip->nonstd_id_htable;
  HtableAddNondup(chr_name, name_slen, kChrHtableSize, chr_code_end, id_htable);
  return kPglRetSuccess;
}


/*
uintptr_t count_11_vecs(const VecW* geno_vvec, uintptr_t vec_ct) {
  // Counts number of aligned 11s in vptr[0..(vec_ct-1)].  Assumes vec_ct is a
  // multiple of 6 (0 ok).
  assert(!(vec_ct % 6));
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW m8 = VCONST_W(kMask00FF);
  const VecW* geno_vvec_iter = geno_vvec;
  const VecW* geno_vvec_end = &(geno_vvec[vec_ct]);
  uintptr_t tot = 0;

  while (1) {
    const VecW* geno_vvec_stop = &(geno_vvec_iter[60]);

    UniVec acc;
    acc.vw = vecw_setzero();

    if (geno_vvec_stop > geno_vvec_end) {
      if (geno_vvec_iter == geno_vvec_end) {
        return tot;
      }
      geno_vvec_stop = geno_vvec_end;
    }
    do {
      VecW cur_geno_vword = *geno_vvec_iter++;
      VecW count1 = cur_geno_vword & m1;
      count1 = count1 & vecw_srli(cur_geno_vword, 1);

      cur_geno_vword = *geno_vvec_iter++;
      VecW cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vecw_srli(cur_geno_vword, 1);
      count1 = count1 + cur_11;

      cur_geno_vword = *geno_vvec_iter++;
      cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vecw_srli(cur_geno_vword, 1);
      count1 = count1 + cur_11;
      count1 = (count1 & m2) + (vecw_srli(count1, 2) & m2);

      cur_geno_vword = *geno_vvec_iter++;
      VecW count2 = cur_geno_vword & m1;
      count2 = count2 & vecw_srli(cur_geno_vword, 1);

      cur_geno_vword = *geno_vvec_iter++;
      VecW cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vecw_srli(cur_geno_vword, 1);
      count2 = count2 + cur_11;

      cur_geno_vword = *geno_vvec_iter++;
      cur_11 = cur_geno_vword & m1;
      cur_11 = cur_11 & vecw_srli(cur_geno_vword, 1);
      count2 = count2 + cur_11;
      count1 = count1 + (count2 & m2) + (vecw_srli(count2, 2) & m2);

      acc.vw = acc.vw + (count1 & m4) + (vecw_srli(count1, 4) & m4);
    } while (geno_vvec_iter < geno_vvec_stop);
    acc.vw = (acc.vw & m8) + (vecw_srli(acc.vw, 8) & m8);
    tot += UniVecHsum16(acc);
  }
}

uintptr_t count_11_longs(const uintptr_t* genovec, uintptr_t word_ct) {
  uintptr_t tot = 0;
  if (word_ct >= (6 * kWordsPerVec)) {
    assert(IsVecAligned(genovec));
    const uintptr_t remainder = word_ct % (6 * kWordsPerVec);
    const uintptr_t main_block_word_ct = word_ct - remainder;
    tot = count_11_vecs((const VecW*)genovec, main_block_word_ct / kWordsPerVec);
    word_ct = remainder;
    genovec = &(genovec[main_block_word_ct]);
  }
  for (uintptr_t trailing_word_idx = 0; trailing_word_idx != word_ct; ++trailing_word_idx) {
    const uintptr_t cur_geno_word = genovec[trailing_word_idx];
    tot += Popcount01Word(Word11(cur_geno_word));
  }
}
*/

uint32_t AllGenoEqual(const uintptr_t* genoarr, uint32_t sample_ct) {
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uintptr_t match_word = (genoarr[0] & 3) * kMask5555;
  for (uint32_t widx = 0; widx != word_ct_m1; ++widx) {
    if (genoarr[widx] != match_word) {
      return 0;
    }
  }
  const uint32_t remainder = ModNz(sample_ct, kBitsPerWordD2);
  return !bzhi_max(genoarr[word_ct_m1] ^ match_word, 2 * remainder);
}

void InterleavedInvert(const uintptr_t* __restrict interleaved_set, uint32_t sample_ctv2, uintptr_t* __restrict genovec) {
  const uintptr_t twovec_ct = sample_ctv2 / 2;
#ifdef __LP64__
  const VecW m1 = VCONST_W(kMask5555);
  const VecW* interleaved_set_iter = R_CAST(const VecW*, interleaved_set);
  VecW* genovvec_iter = R_CAST(VecW*, genovec);
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const VecW set_vvec = *interleaved_set_iter++;
    const VecW set_first_lo = set_vvec & m1;
    const VecW set_second_hi = vecw_and_notfirst(m1, set_vvec);
    const VecW genovvec_first = *genovvec_iter;
    // Flip high bit iff (low geno bit is unset) & interleaved_set.
    *genovvec_iter = genovvec_first ^ vecw_slli(vecw_and_notfirst(genovvec_first, set_first_lo), 1);
    ++genovvec_iter;
    const VecW genovvec_second = *genovvec_iter;
    *genovvec_iter = genovvec_second ^ vecw_and_notfirst(vecw_slli(genovvec_second, 1), set_second_hi);
    ++genovvec_iter;
  }
  if (sample_ctv2 & 1) {
    const VecW set_vvec = *interleaved_set_iter;
    const VecW set_first_lo = set_vvec & m1;
    const VecW genovvec_first = *genovvec_iter;
    *genovvec_iter = genovvec_first ^ vecw_slli(vecw_and_notfirst(genovvec_first, set_first_lo), 1);
  }
#else
  const uintptr_t* interleaved_set_iter = interleaved_set;
  uintptr_t* genovec_iter = genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const uintptr_t set_word = *interleaved_set_iter++;
    const uintptr_t set_first_lo = set_word & kMask5555;
    const uintptr_t set_second_hi = (~kMask5555) & set_word;
    const uintptr_t genoword_first = *genovec_iter;
    *genovec_iter = genoword_first ^ (((~genoword_first) & set_first_lo) << 1);
    ++genovec_iter;
    const uintptr_t genoword_second = *genovec_iter;
    *genovec_iter = genoword_second ^ ((~(genoword_second << 1)) & set_second_hi);
    ++genovec_iter;
  }
  if (sample_ctv2 & 1) {
    const uintptr_t set_word = *interleaved_set_iter;
    const uintptr_t set_first_lo = set_word & kMask5555;
    const uintptr_t genoword_first = *genovec_iter;
    *genovec_iter = genoword_first ^ (((~genoword_first) & set_first_lo) << 1);
  }
#endif
}

uint32_t InterleavedSubsetIs03(const uintptr_t* genovec, const uintptr_t* interleaved_set, uint32_t sample_ctv2) {
  const uintptr_t twovec_ct = sample_ctv2 / 2;
#ifdef __LP64__
  const VecU16 m1 = VCONST_S(0x5555);
  const VecU16* interleaved_set_iter = R_CAST(const VecU16*, interleaved_set);
  const VecU16* genovvec_iter = R_CAST(const VecU16*, genovec);
  // This is expected to return true, so it's reasonable to bitwise-or all
  // intermediate failures and only check for a nonzero bit at the end.
  VecU16 result_vvec = vecu16_setzero();
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const VecU16 set_vvec = *interleaved_set_iter++;
    const VecU16 set_first_lo = set_vvec & m1;
    const VecU16 set_second_hi = vecu16_and_notfirst(m1, set_vvec);
    const VecU16 genovvec_first = *genovvec_iter++;
    const VecU16 genovvec_second = *genovvec_iter++;
    result_vvec |= (genovvec_first ^ vecu16_srli(genovvec_first, 1)) & set_first_lo;
    result_vvec |= (vecu16_slli(genovvec_second, 1) ^ genovvec_second) & set_second_hi;
  }
  if (sample_ctv2 & 1) {
    const VecU16 set_vvec = *interleaved_set_iter;
    const VecU16 set_first_lo = set_vvec & m1;
    const VecU16 genovvec_first = *genovvec_iter;
    result_vvec |= (genovvec_first ^ vecu16_srli(genovvec_first, 1)) & set_first_lo;
  }
  const VecU16 zero_vvec = vecu16_setzero();
  return vecu16s_are_equal(result_vvec, zero_vvec);
#else
  const uintptr_t* interleaved_set_iter = interleaved_set;
  const uintptr_t* genovec_iter = genovec;
  uintptr_t result_word = 0;
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const uintptr_t set_word = *interleaved_set_iter++;
    const uintptr_t set_first_lo = set_word & kMask5555;
    const uintptr_t set_second_hi = (~kMask5555) & set_word;
    const uintptr_t genoword_first = *genovec_iter++;
    const uintptr_t genoword_second = *genovec_iter++;
    result_word |= (genoword_first ^ (genoword_first >> 1)) & set_first_lo;
    result_word |= ((genoword_second << 1) ^ genoword_second) & set_second_hi;
  }
  if (sample_ctv2 & 1) {
    const uintptr_t set_word = *interleaved_set_iter;
    const uintptr_t set_first_lo = set_word & kMask5555;
    const uintptr_t genoword_first = *genovec_iter;
    result_word |= (genoword_first ^ (genoword_first >> 1)) & set_first_lo;
  }
  return (result_word == 0);
#endif
}

uint32_t InterleavedSubsetIs23(const uintptr_t* genovec, const uintptr_t* interleaved_set, uint32_t sample_ctv2) {
  const uintptr_t twovec_ct = sample_ctv2 / 2;
#ifdef __LP64__
  const VecU16 m1 = VCONST_S(0x5555);
  const VecU16* interleaved_set_iter = R_CAST(const VecU16*, interleaved_set);
  const VecU16* genovvec_iter = R_CAST(const VecU16*, genovec);
  // This is expected to return true, so it's reasonable to bitwise-or all
  // intermediate failures and only check for a nonzero bit at the end.
  VecU16 result_vvec = vecu16_setzero();
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const VecU16 set_vvec = *interleaved_set_iter++;
    const VecU16 set_first_hi = vecu16_slli(set_vvec & m1, 1);
    const VecU16 set_second_hi = vecu16_and_notfirst(m1, set_vvec);
    result_vvec |= vecu16_and_notfirst(*genovvec_iter++, set_first_hi);
    result_vvec |= vecu16_and_notfirst(*genovvec_iter++, set_second_hi);
  }
  if (sample_ctv2 & 1) {
    const VecU16 set_vvec = *interleaved_set_iter;
    const VecU16 set_first_hi = vecu16_slli(set_vvec & m1, 1);
    result_vvec |= vecu16_and_notfirst(*genovvec_iter++, set_first_hi);
  }
  const VecU16 zero_vvec = vecu16_setzero();
  return vecu16s_are_equal(result_vvec, zero_vvec);
#else
  const uintptr_t* interleaved_set_iter = interleaved_set;
  const uintptr_t* genovec_iter = genovec;
  uintptr_t result_word = 0;
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const uintptr_t set_word = *interleaved_set_iter++;
    const uintptr_t set_first_hi = (set_word & kMask5555) << 1;
    const uintptr_t set_second_hi = (~kMask5555) & set_word;
    result_word |= (~(*genovec_iter++)) & set_first_hi;
    result_word |= (~(*genovec_iter++)) & set_second_hi;
  }
  if (sample_ctv2 & 1) {
    const uintptr_t set_word = *interleaved_set_iter;
    const uintptr_t set_first_hi = (set_word & kMask5555) << 1;
    result_word |= (~(*genovec_iter)) & set_first_hi;
  }
  return (result_word == 0);
#endif
}

uint32_t InterleavedSubsetAllMissing(const uintptr_t* genovec, const uintptr_t* interleaved_set, uint32_t sample_ctv2) {
  const uintptr_t twovec_ct = sample_ctv2 / 2;
#ifdef __LP64__
  const VecU16 m1 = VCONST_S(0x5555);
  const VecU16* interleaved_set_iter = R_CAST(const VecU16*, interleaved_set);
  const VecU16* genovvec_iter = R_CAST(const VecU16*, genovec);
  // This is expected to return true, so it's reasonable to bitwise-or all
  // intermediate failures and only check for a nonzero bit at the end.
  VecU16 result_vvec = vecu16_setzero();
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const VecU16 set_vvec = *interleaved_set_iter++;
    const VecU16 set_first_lo = set_vvec & m1;
    const VecU16 set_second_hi = vecu16_and_notfirst(m1, set_vvec);
    const VecU16 set_first = set_first_lo | vecu16_slli(set_first_lo, 1);
    const VecU16 set_second = set_second_hi | vecu16_srli(set_second_hi, 1);
    result_vvec |= set_first & (set_first ^ (*genovvec_iter++));
    result_vvec |= set_second & (set_second ^ (*genovvec_iter++));
  }
  if (sample_ctv2 & 1) {
    const VecU16 set_vvec = *interleaved_set_iter;
    const VecU16 set_first_lo = set_vvec & m1;
    const VecU16 set_first = set_first_lo | vecu16_slli(set_first_lo, 1);
    result_vvec |= set_first & (set_first ^ (*genovvec_iter));
  }
  const VecU16 zero_vvec = vecu16_setzero();
  return vecu16s_are_equal(result_vvec, zero_vvec);
#else
  const uintptr_t* interleaved_set_iter = interleaved_set;
  const uintptr_t* genovec_iter = genovec;
  uintptr_t result_word = 0;
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const uintptr_t set_word = *interleaved_set_iter++;
    const uintptr_t set_first_lo = set_word & kMask5555;
    const uintptr_t set_second_hi = (~kMask5555) & set_word;
    const uintptr_t set_first = set_first_lo * 3;
    const uintptr_t set_second = set_second_hi | (set_second_hi >> 1);
    result_word |= set_first & (set_first ^ (*genovec_iter++));
    result_word |= set_second & (set_first ^ (*genovec_iter++));
  }
  if (sample_ctv2 & 1) {
    const uintptr_t set_word = *interleaved_set_iter;
    const uintptr_t set_first_lo = set_word & kMask5555;
    const uintptr_t set_first = set_first_lo * 3;
    result_word |= set_first & (set_first ^ (*genovec_iter));
  }
  return (result_word == 0);
#endif
}

void InterleavedMaskZero(const uintptr_t* __restrict interleaved_mask, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec) {
  const uintptr_t twovec_ct = geno_vec_ct / 2;
#ifdef __LP64__
  const VecW m1 = VCONST_W(kMask5555);
  const VecW* interleaved_mask_iter = R_CAST(const VecW*, interleaved_mask);
  VecW* genovvec_iter = R_CAST(VecW*, genovec);
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const VecW mask_vvec = *interleaved_mask_iter++;
    VecW mask_first = mask_vvec & m1;
    mask_first = mask_first | vecw_slli(mask_first, 1);
    VecW mask_second = vecw_and_notfirst(m1, mask_vvec);
    mask_second = mask_second | vecw_srli(mask_second, 1);
    *genovvec_iter = (*genovvec_iter) & mask_first;
    ++genovvec_iter;
    *genovvec_iter = (*genovvec_iter) & mask_second;
    ++genovvec_iter;
  }
  if (geno_vec_ct & 1) {
    VecW mask_first = *interleaved_mask_iter;
    mask_first = mask_first | vecw_slli(mask_first, 1);
    *genovvec_iter = (*genovvec_iter) & mask_first;
  }
#else
  const uintptr_t* interleaved_mask_iter = interleaved_mask;
  uintptr_t* genovec_iter = genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const uintptr_t mask_word = *interleaved_mask_iter++;
    *genovec_iter &= (mask_word & kMask5555) * 3;
    ++genovec_iter;
    *genovec_iter &= ((mask_word >> 1) & kMask5555) * 3;
    ++genovec_iter;
  }
  if (geno_vec_ct & 1) {
    const uintptr_t mask_word = *interleaved_mask_iter;
    *genovec_iter &= mask_word * 3;
  }
#endif
}

/*
// Set genovec entry missing iff interleaved_mask bit is NOT set.  Trailing
// genovec bits are dirtied.
void InterleavedMaskMissing(const uintptr_t* __restrict interleaved_mask, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec) {
  const uintptr_t twovec_ct = geno_vec_ct / 2;
#ifdef __LP64__
  const VecW m1 = VCONST_W(kMask5555);
  const VecW inv_m1 = VCONST_W(kMaskAAAA);
  const VecW* interleaved_mask_iter = R_CAST(const VecW*, interleaved_mask);
  VecW* genovvec_iter = R_CAST(VecW*, genovec);
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const VecW mask_vvec = *interleaved_mask_iter++;
    VecW set_first = vecw_and_notfirst(mask_vvec, m1);
    set_first = set_first | vecw_slli(set_first, 1);
    VecW set_second = vecw_and_notfirst(mask_vvec, inv_m1);
    set_second = set_second | vecw_srli(set_second, 1);
    *genovvec_iter = (*genovvec_iter) | set_first;
    ++genovvec_iter;
    *genovvec_iter = (*genovvec_iter) | set_second;
    ++genovvec_iter;
  }
  if (geno_vec_ct & 1) {
    VecW set_first = vecw_and_notfirst(*interleaved_mask_iter, m1);
    set_first = set_first | vecw_slli(set_first, 1);
    *genovvec_iter = (*genovvec_iter) | set_first;
  }
#else
  const uintptr_t* interleaved_mask_iter = interleaved_mask;
  uintptr_t* genovec_iter = genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const uintptr_t set_invword = *interleaved_mask_iter++;
    *genovec_iter |= ((~set_invword) & kMask5555) * 3;
    ++genovec_iter;
    *genovec_iter |= (((~set_invword) >> 1) & kMask5555) * 3;
    ++genovec_iter;
  }
  if (geno_vec_ct & 1) {
    const uintptr_t set_invword = *interleaved_mask_iter;
    *genovec_iter |= ((~set_invword) & kMask5555) * 3;
  }
#endif
}
*/

void InterleavedSetMissing(const uintptr_t* __restrict interleaved_set, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec) {
  const uintptr_t twovec_ct = geno_vec_ct / 2;
#ifdef __LP64__
  const VecW m1 = VCONST_W(kMask5555);
  const VecW* interleaved_set_iter = R_CAST(const VecW*, interleaved_set);
  VecW* genovvec_iter = R_CAST(VecW*, genovec);
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const VecW set_vvec = *interleaved_set_iter++;
    VecW set_first = set_vvec & m1;
    set_first = set_first | vecw_slli(set_first, 1);
    VecW set_second = vecw_and_notfirst(m1, set_vvec);
    set_second = set_second | vecw_srli(set_second, 1);
    *genovvec_iter = (*genovvec_iter) | set_first;
    ++genovvec_iter;
    *genovvec_iter = (*genovvec_iter) | set_second;
    ++genovvec_iter;
  }
  if (geno_vec_ct & 1) {
    VecW set_first = *interleaved_set_iter;
    set_first = set_first | vecw_slli(set_first, 1);
    *genovvec_iter = (*genovvec_iter) | set_first;
  }
#else
  const uintptr_t* interleaved_set_iter = interleaved_set;
  uintptr_t* genovec_iter = genovec;
  for (uintptr_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const uintptr_t set_word = *interleaved_set_iter++;
    *genovec_iter |= (set_word & kMask5555) * 3;
    ++genovec_iter;
    *genovec_iter |= ((set_word >> 1) & kMask5555) * 3;
    ++genovec_iter;
  }
  if (geno_vec_ct & 1) {
    const uintptr_t set_word = *interleaved_set_iter;
    *genovec_iter |= set_word * 3;
  }
#endif
}

void InterleavedSetMissingCleardosage(const uintptr_t* __restrict orig_set, const uintptr_t* __restrict interleaved_set, uintptr_t geno_vec_ct, uintptr_t* __restrict genovec, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main) {
  const uint32_t orig_write_dosage_ct = *write_dosage_ct_ptr;
  if (orig_write_dosage_ct) {
    uintptr_t sample_widx = 0;
    uintptr_t cur_bits = dosagepresent[0];
    for (uint32_t dosage_read_idx = 0; dosage_read_idx != orig_write_dosage_ct; ++dosage_read_idx) {
      const uintptr_t lowbit = BitIter1y(dosagepresent, &sample_widx, &cur_bits);
      if (orig_set[sample_widx] & lowbit) {
        dosagepresent[sample_widx] ^= lowbit;
        uint32_t dosage_write_idx = dosage_read_idx++;
        for (; dosage_read_idx == orig_write_dosage_ct; ++dosage_read_idx) {
          const uintptr_t lowbit2 = BitIter1y(dosagepresent, &sample_widx, &cur_bits);
          if (orig_set[sample_widx] & lowbit2) {
            dosagepresent[sample_widx] ^= lowbit2;
          } else {
            dosage_main[dosage_write_idx++] = dosage_main[dosage_read_idx];
          }
        }
        *write_dosage_ct_ptr = dosage_write_idx;
        break;
      }
    }
  }
  InterleavedSetMissing(interleaved_set, geno_vec_ct, genovec);
}

// may move this to pgenlib_misc
void EraseDosages(const uintptr_t* __restrict erase_map, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main) {
  const uint32_t orig_write_dosage_ct = *write_dosage_ct_ptr;
  if (!orig_write_dosage_ct) {
    return;
  }
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosagepresent[0];
  for (uint32_t dosage_read_idx = 0; dosage_read_idx != orig_write_dosage_ct; ++dosage_read_idx) {
    uintptr_t sample_uidx = BitIter1(dosagepresent, &sample_uidx_base, &cur_bits);
    if (IsSet(erase_map, sample_uidx)) {
      ClearBit(sample_uidx, dosagepresent);
      uint32_t dosage_write_idx = dosage_read_idx++;
      for (; dosage_read_idx != orig_write_dosage_ct; ++dosage_read_idx) {
        sample_uidx = BitIter1(dosagepresent, &sample_uidx_base, &cur_bits);
        if (IsSet(erase_map, sample_uidx)) {
          ClearBit(sample_uidx, dosagepresent);
        } else {
          dosage_main[dosage_write_idx++] = dosage_main[dosage_read_idx];
        }
      }
      *write_dosage_ct_ptr = dosage_write_idx;
      return;
    }
  }
}

void SetMaleHetMissing(const uintptr_t* __restrict sex_male_interleaved, uint32_t geno_vec_ct, uintptr_t* __restrict genovec) {
  const uint32_t twovec_ct = geno_vec_ct / 2;
#ifdef __LP64__
  const VecW m1 = VCONST_W(kMask5555);
  const VecW* sex_male_interleaved_iter = R_CAST(const VecW*, sex_male_interleaved);
  VecW* genovvec_iter = R_CAST(VecW*, genovec);
  for (uint32_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const VecW sex_male_vvec = *sex_male_interleaved_iter++;
    // we wish to bitwise-or with (sex_male_nypvec_01 & genovec) << 1
    const VecW sex_male_first = sex_male_vvec & m1;
    const VecW sex_male_second_shifted = vecw_and_notfirst(m1, sex_male_vvec);
    VecW cur_geno_vword = *genovvec_iter;

    const VecW missing_male_vword = sex_male_first & cur_geno_vword;

    *genovvec_iter++ = cur_geno_vword | vecw_slli(missing_male_vword, 1);
    cur_geno_vword = *genovvec_iter;
    *genovvec_iter++ = cur_geno_vword | (sex_male_second_shifted & vecw_slli(cur_geno_vword, 1));
  }
  if (geno_vec_ct & 1) {
    const VecW sex_male_first = (*sex_male_interleaved_iter) & m1;
    const VecW cur_geno_vword = *genovvec_iter;
    const VecW missing_male_vword = sex_male_first & cur_geno_vword;
    *genovvec_iter = cur_geno_vword | vecw_slli(missing_male_vword, 1);
  }
#else
  const uintptr_t* sex_male_interleaved_iter = sex_male_interleaved;
  uintptr_t* genovec_iter = genovec;
  for (uint32_t twovec_idx = 0; twovec_idx != twovec_ct; ++twovec_idx) {
    const uintptr_t sex_male_word = *sex_male_interleaved_iter++;
    uintptr_t cur_geno_word = *genovec_iter;
    *genovec_iter++ = cur_geno_word | ((sex_male_word & kMask5555 & cur_geno_word) << 1);
    cur_geno_word = *genovec_iter;
    *genovec_iter++ = cur_geno_word | (sex_male_word & kMaskAAAA & (cur_geno_word << 1));
  }
  if (geno_vec_ct & 1) {
    const uintptr_t sex_male_word = *sex_male_interleaved_iter;
    uintptr_t cur_geno_word = *genovec_iter;
    *genovec_iter = cur_geno_word | ((sex_male_word & kMask5555 & cur_geno_word) << 1);
  }
#endif
}

void EraseMaleHetDosages(const uintptr_t* __restrict sex_male, const uintptr_t* __restrict genoarr, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main) {
  const uint32_t orig_write_dosage_ct = *write_dosage_ct_ptr;
  if (!orig_write_dosage_ct) {
    return;
  }
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = dosagepresent[0];
  for (uint32_t dosage_read_idx = 0; dosage_read_idx != orig_write_dosage_ct; ++dosage_read_idx) {
    const uintptr_t sample_uidx = BitIter1(dosagepresent, &sample_uidx_base, &cur_bits);
    if (IsSet(sex_male, sample_uidx) && (GetNyparrEntry(genoarr, sample_uidx) == 1)) {
      ClearBit(sample_uidx, dosagepresent);
      uint32_t dosage_write_idx = dosage_read_idx++;
      for (; dosage_read_idx != orig_write_dosage_ct; ++dosage_read_idx) {
        const uintptr_t sample_uidx2 = BitIter1(dosagepresent, &sample_uidx_base, &cur_bits);
        if (IsSet(sex_male, sample_uidx2) && (GetNyparrEntry(genoarr, sample_uidx2) == 1)) {
          ClearBit(sample_uidx2, dosagepresent);
        } else {
          dosage_main[dosage_write_idx++] = dosage_main[dosage_read_idx];
        }
      }
      *write_dosage_ct_ptr = dosage_write_idx;
      return;
    }
  }
}

void EraseMaleDphases(const uintptr_t* __restrict sex_male, uint32_t* __restrict write_dphase_ct_ptr, uintptr_t* __restrict dphasepresent, SDosage* dphase_delta) {
  const uint32_t orig_write_dphase_ct = *write_dphase_ct_ptr;
  if (!orig_write_dphase_ct) {
    return;
  }
  uintptr_t sample_widx = 0;
  uintptr_t cur_bits = dphasepresent[0];
  for (uint32_t dphase_read_idx = 0; dphase_read_idx != orig_write_dphase_ct; ++dphase_read_idx) {
    const uintptr_t lowbit = BitIter1y(dphasepresent, &sample_widx, &cur_bits);
    if (sex_male[sample_widx] & lowbit) {
      dphasepresent[sample_widx] ^= lowbit;
      uint32_t dphase_write_idx = dphase_read_idx++;
      for (; dphase_read_idx == orig_write_dphase_ct; ++dphase_read_idx) {
        const uintptr_t lowbit2 = BitIter1y(dphasepresent, &sample_widx, &cur_bits);
        if (sex_male[sample_widx] & lowbit2) {
          dphasepresent[sample_widx] ^= lowbit2;
        } else {
          dphase_delta[dphase_write_idx++] = dphase_delta[dphase_read_idx];
        }
      }
      *write_dphase_ct_ptr = dphase_write_idx;
      return;
    }
  }
}

void SetMaleHetMissingCleardosage(const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_male_interleaved, uint32_t geno_vec_ct, uintptr_t* __restrict genovec, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main) {
  if (*write_dosage_ct_ptr) {
    EraseMaleHetDosages(sex_male, genovec, write_dosage_ct_ptr, dosagepresent, dosage_main);
  }
  SetMaleHetMissing(sex_male_interleaved, geno_vec_ct, genovec);
}

void SetMaleHetMissingKeepdosage(const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_male_interleaved, uint32_t geno_word_ct, uintptr_t* __restrict genovec, uint32_t* __restrict write_dosage_ct_ptr, uintptr_t* __restrict dosagepresent, Dosage* dosage_main) {
  // count # of 1.00000 dosages we need to insert, and then rewrite dosage_main
  // from back to front so we don't need temporary buffers.
  const uint32_t orig_dosage_ct = *write_dosage_ct_ptr;
  if (!orig_dosage_ct) {
    // can't assume dosagepresent is initialized in this case
    ZeroWArr(DivUp(geno_word_ct, 2), dosagepresent);
  }
  const Halfword* sex_male_alias = R_CAST(const Halfword*, sex_male);
  Halfword* dosagepresent_alias = R_CAST(Halfword*, dosagepresent);
  uint32_t new_dosage_ct = 0;
  // todo: try vectorizing this loop
  for (uintptr_t widx = 0; widx != geno_word_ct; ++widx) {
    const uintptr_t geno_hets = Word01(genovec[widx]);
    const uintptr_t male_nodosage_word = UnpackHalfwordToWord(sex_male_alias[widx] & (~dosagepresent_alias[widx]));
    new_dosage_ct += Popcount01Word(geno_hets & male_nodosage_word);
  }
  if (!new_dosage_ct) {
    SetMaleHetMissing(sex_male_interleaved, WordCtToVecCt(geno_word_ct), genovec);
    return;
  }
  uint32_t dosage_write_idx = orig_dosage_ct + new_dosage_ct;
  uint32_t dosage_read_idx = orig_dosage_ct;
  uint32_t widx = geno_word_ct;
  *write_dosage_ct_ptr = dosage_write_idx;
  do {
    --widx;
    const uintptr_t geno_word = genovec[widx];
    const uintptr_t male_geno_hets = geno_word & (~(geno_word >> 1)) & UnpackHalfwordToWord(sex_male_alias[widx]);
    const uintptr_t dosagepresent_word = UnpackHalfwordToWord(dosagepresent_alias[widx]);
    uintptr_t new_dosagepresent_word = dosagepresent_word | male_geno_hets;
    if (new_dosagepresent_word) {
      dosagepresent_alias[widx] = PackWordToHalfword(new_dosagepresent_word);
      do {
        const uint32_t top_set_bit = bsrw(new_dosagepresent_word);
        const uintptr_t cur_bit_word = k1LU << top_set_bit;
        Dosage cur_dosage = kDosageMid;
        if (cur_bit_word & dosagepresent_word) {
          cur_dosage = dosage_main[--dosage_read_idx];
        }
        dosage_main[--dosage_write_idx] = cur_dosage;
        new_dosagepresent_word -= cur_bit_word;
      } while (new_dosagepresent_word);
      genovec[widx] = geno_word | (male_geno_hets << 1);
    }
  } while (widx);
}

// Clears each bit in bitarr which doesn't correspond to a genovec het.
// Assumes that either trailing bits of bitarr are already zero, or trailing
// bits of genovec are zero.
//
// Similar to PgrDetectGenoarrHetsUnsafe().
void MaskGenoarrHetsMultiallelicUnsafe(const uintptr_t* __restrict genoarr, const uintptr_t* __restrict patch_10_set, const AlleleCode* __restrict patch_10_vals, uint32_t raw_sample_ctl2, uintptr_t* __restrict bitarr) {
  // Related to PgrDetectGenoarrHetsMultiallelic().
  const Halfword* patch_10_set_alias = R_CAST(const Halfword*, patch_10_set);
  const AlleleCode* patch_10_vals_iter = patch_10_vals;
  Halfword* bitarr_alias = R_CAST(Halfword*, bitarr);
  for (uint32_t widx = 0; widx != raw_sample_ctl2; ++widx) {
    const uintptr_t cur_word = genoarr[widx];
    uint32_t patch_10_hw = patch_10_set_alias[widx];
    uint32_t cur_hets = Pack01ToHalfword(cur_word);
    while (patch_10_hw) {
      const AlleleCode code1 = *patch_10_vals_iter++;
      const AlleleCode code2 = *patch_10_vals_iter++;
      if (code1 != code2) {
        cur_hets |= ctzu32(patch_10_hw);
      }
      patch_10_hw &= patch_10_hw - 1;
    }
    bitarr_alias[widx] &= cur_hets;
  }
}

void SetAltxyHetMissing(uintptr_t* __restrict genoarr, uint32_t* __restrict rare10_ct_ptr, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals) {
  const uint32_t orig_rare10_ct = *rare10_ct_ptr;
  uintptr_t sample_widx = 0;
  uintptr_t patch_10_bits = patch_10_set[0];
  uint32_t read_patch_10_idx = 0;
  while (1) {
    const uintptr_t lowbit = BitIter1y(patch_10_set, &sample_widx, &patch_10_bits);
    const AlleleCode lo_code = patch_10_vals[read_patch_10_idx * 2];
    const AlleleCode hi_code = patch_10_vals[read_patch_10_idx * 2 + 1];
    if (lo_code != hi_code) {
      // change genoarr nyp from 2 to 3, by setting the low bit.
      // (note that one of these genoarr operations does nothing due to
      // overflow/underflow; there may be a better way of doing this)
      genoarr[sample_widx * 2] |= lowbit * lowbit;
      genoarr[sample_widx * 2 + 1] |= (lowbit >> kBitsPerWordD2) * (lowbit >> kBitsPerWordD2);
      patch_10_set[sample_widx] ^= lowbit;
      break;
    }
    if (++read_patch_10_idx == orig_rare10_ct) {
      return;
    }
  }
  uint32_t write_patch_10_idx = read_patch_10_idx;
  while (++read_patch_10_idx != orig_rare10_ct) {
    const uintptr_t lowbit = BitIter1y(patch_10_set, &sample_widx, &patch_10_bits);
    const AlleleCode lo_code = patch_10_vals[read_patch_10_idx * 2];
    const AlleleCode hi_code = patch_10_vals[read_patch_10_idx * 2 + 1];
    if (lo_code != hi_code) {
      genoarr[sample_widx * 2] |= lowbit * lowbit;
      genoarr[sample_widx * 2 + 1] |= (lowbit >> kBitsPerWordD2) * (lowbit >> kBitsPerWordD2);
      patch_10_set[sample_widx] ^= lowbit;
    } else {
      patch_10_vals[write_patch_10_idx * 2] = lo_code;
      patch_10_vals[write_patch_10_idx * 2 + 1] = hi_code;
      ++write_patch_10_idx;
    }
  }
  *rare10_ct_ptr = write_patch_10_idx;
}

void SetMaleAltxyHetMissing(const uintptr_t* sex_male, uintptr_t* __restrict genoarr, uint32_t* __restrict rare10_ct_ptr, uintptr_t* __restrict patch_10_set, AlleleCode* __restrict patch_10_vals) {
  const uint32_t orig_rare10_ct = *rare10_ct_ptr;
  uintptr_t sample_widx = 0;
  uintptr_t patch_10_bits = patch_10_set[0];
  uint32_t read_patch_10_idx = 0;
  while (1) {
    const uintptr_t lowbit = BitIter1y(patch_10_set, &sample_widx, &patch_10_bits);
    if (sex_male[sample_widx] & lowbit) {
      const AlleleCode lo_code = patch_10_vals[read_patch_10_idx * 2];
      const AlleleCode hi_code = patch_10_vals[read_patch_10_idx * 2 + 1];
      if (lo_code != hi_code) {
        // change genoarr nyp from 2 to 3, by setting the low bit.
        // (note that one of these genoarr operations does nothing due to
        // overflow/underflow; there may be a better way of doing this)
        genoarr[sample_widx * 2] |= lowbit * lowbit;
        genoarr[sample_widx * 2 + 1] |= (lowbit >> kBitsPerWordD2) * (lowbit >> kBitsPerWordD2);
        patch_10_set[sample_widx] ^= lowbit;
        break;
      }
    }
    if (++read_patch_10_idx == orig_rare10_ct) {
      return;
    }
  }
  uint32_t write_patch_10_idx = read_patch_10_idx;
  while (++read_patch_10_idx != orig_rare10_ct) {
    const uintptr_t lowbit = BitIter1y(patch_10_set, &sample_widx, &patch_10_bits);
    const AlleleCode lo_code = patch_10_vals[read_patch_10_idx * 2];
    const AlleleCode hi_code = patch_10_vals[read_patch_10_idx * 2 + 1];
    if ((sex_male[sample_widx] & lowbit) && (lo_code != hi_code)) {
      genoarr[sample_widx * 2] |= lowbit * lowbit;
      genoarr[sample_widx * 2 + 1] |= (lowbit >> kBitsPerWordD2) * (lowbit >> kBitsPerWordD2);
      patch_10_set[sample_widx] ^= lowbit;
    } else {
      patch_10_vals[write_patch_10_idx * 2] = lo_code;
      patch_10_vals[write_patch_10_idx * 2 + 1] = hi_code;
      ++write_patch_10_idx;
    }
  }
  *rare10_ct_ptr = write_patch_10_idx;
}

/*
uint32_t chr_window_max(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bp, uint32_t chr_fo_idx, uint32_t ct_max, uint32_t bp_max, uint32_t cur_window_max) {
  if (cur_window_max >= ct_max) {
    return ct_max;
  }
  const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  uint32_t variant_uidx = AdvBoundedTo1Bit(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], chr_end);
  const uint32_t variant_ct = PopcountBitRange(variant_include, variant_uidx, chr_end);
  if (variant_ct <= cur_window_max) {
    return cur_window_max;
  }
  uint32_t window_idx_first = 0;
  uint32_t window_uidx_first = variant_uidx;
  uint32_t window_bp_first = variant_bp[variant_uidx];
  for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_uidx, ++variant_idx) {
    MovU32To1Bit(variant_include, &variant_uidx);
    uint32_t variant_bp_thresh = variant_bp[variant_uidx];
    if (variant_bp_thresh < bp_max) {
      variant_bp_thresh = 0;
    } else {
      variant_bp_thresh -= bp_max;
    }
    if (variant_bp_thresh > window_bp_first) {
      do {
        ++window_uidx_first;
        MovU32To1Bit(variant_include, &window_uidx_first);
        window_bp_first = variant_bp[window_uidx_first];
        ++window_idx_first;
      } while (variant_bp_thresh > window_bp_first);
    } else if (variant_idx - window_idx_first == cur_window_max) {
      if (++cur_window_max == ct_max) {
        return cur_window_max;
      }
    }
  }
  return cur_window_max;
}
*/

uint32_t NotOnlyXymt(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t raw_variant_ct, uint32_t xymt_offset) {
  const uint32_t xymt_code = cip->xymt_codes[xymt_offset];
  assert(!IsI32Neg(xymt_code));
  const uint32_t cur_chr_fo_idx = cip->chr_idx_to_foidx[xymt_code];
  const uint32_t chr_start = cip->chr_fo_vidx_start[cur_chr_fo_idx];
  if (chr_start) {
    const uint32_t first_uidx = AdvTo1Bit(variant_include, 0);
    if (first_uidx < chr_start) {
      return 1;
    }
  }
  const uint32_t chr_end = cip->chr_fo_vidx_start[cur_chr_fo_idx + 1];
  return (chr_end < raw_variant_ct) && (!AllBitsAreZero(variant_include, chr_end, raw_variant_ct));
}

uint32_t CountNonAutosomalVariants(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t count_x, uint32_t count_mt) {
  // for backward compatibility, unplaced markers are considered to be
  // autosomal here
  uint32_t ct = 0;
  if (count_x) {
    uint32_t x_code;
    if (XymtExists(cip, kChrOffsetX, &x_code)) {
      ct += CountChrVariantsUnsafe(variant_include, cip, x_code);
    }
  }
  uint32_t y_code;
  if (XymtExists(cip, kChrOffsetY, &y_code)) {
    ct += CountChrVariantsUnsafe(variant_include, cip, y_code);
  }
  if (count_mt) {
    uint32_t mt_code;
    if (XymtExists(cip, kChrOffsetMT, &mt_code)) {
      ct += CountChrVariantsUnsafe(variant_include, cip, mt_code);
    }
  }
  return ct;
}

void ExcludeNonAutosomalVariants(const ChrInfo* cip, uintptr_t* variant_include) {
  uint32_t x_code;
  if (XymtExists(cip, kChrOffsetX, &x_code)) {
    const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[x_code];
    ClearBitsNz(cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], variant_include);
  }
  uint32_t y_code;
  if (XymtExists(cip, kChrOffsetY, &y_code)) {
    const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[y_code];
    ClearBitsNz(cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], variant_include);
  }
  uint32_t mt_code;
  if (XymtExists(cip, kChrOffsetMT, &mt_code)) {
    const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[mt_code];
    ClearBitsNz(cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1], variant_include);
  }
}

PglErr ConditionalAllocateNonAutosomalVariants(const ChrInfo* cip, const char* calc_descrip, uint32_t raw_variant_ct, const uintptr_t** variant_include_ptr, uint32_t* variant_ct_ptr) {
  const uint32_t non_autosomal_variant_ct = CountNonAutosomalVariants(*variant_include_ptr, cip, 1, 1);
  if (!non_autosomal_variant_ct) {
    return kPglRetSuccess;
  }
  *variant_ct_ptr -= non_autosomal_variant_ct;
  if (calc_descrip) {
    logprintf("Excluding %u variant%s on non-autosomes from %s.\n", non_autosomal_variant_ct, (non_autosomal_variant_ct == 1)? "" : "s", calc_descrip);
    if (unlikely(!(*variant_ct_ptr))) {
      logerrprintf("Error: No variants remaining for %s.\n", calc_descrip);
      return kPglRetDegenerateData;
    }
  }
  const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
  uintptr_t* working_variant_include;
  if (unlikely(bigstack_alloc_w(raw_variant_ctl, &working_variant_include))) {
    return kPglRetNomem;
  }
  memcpy(working_variant_include, *variant_include_ptr, raw_variant_ctl * sizeof(intptr_t));
  ExcludeNonAutosomalVariants(cip, working_variant_include);
  *variant_include_ptr = working_variant_include;
  return kPglRetSuccess;
}

void FillSubsetChrFoVidxStart(const uintptr_t* variant_include, const ChrInfo* cip, uint32_t* subset_chr_fo_vidx_start) {
  const uint32_t chr_ct = cip->chr_ct;
  subset_chr_fo_vidx_start[0] = 0;
  uint32_t variant_uidx = 0;
  uint32_t variant_idx = 0;
  for (uint32_t chr_fo_idx = 1; chr_fo_idx <= chr_ct; ++chr_fo_idx) {
    const uint32_t chr_end_variant_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
    variant_idx += PopcountBitRange(variant_include, variant_uidx, chr_end_variant_uidx);
    subset_chr_fo_vidx_start[chr_fo_idx] = variant_idx;
    variant_uidx = chr_end_variant_uidx;
  }
}

/*
BoolErr allele_set(const char* newval, uint32_t allele_slen, char** allele_ptr) {
  char* newptr;
  if (allele_slen == 1) {
    // const_cast
    newptr = (char*)((uintptr_t)(&(g_one_char_strs[((unsigned char)(*newval)) * 2])));
  } else {
    char* new_alloc;
    if (pgl_malloc(allele_slen + 1, &new_alloc)) {
      return 1;
    }
    memcpyx(new_alloc, newval, allele_slen, '\0');
    newptr = new_alloc;
  }
  *allele_ptr = newptr;
  return 0;
}

BoolErr allele_reset(const char* newval, uint32_t allele_slen, char** allele_ptr) {
  char* newptr;
  if (allele_slen == 1) {
    // const_cast
    newptr = (char*)((uintptr_t)(&(g_one_char_strs[((unsigned char)(*newval)) * 2])));
  } else {
    char* new_alloc;
    if (pgl_malloc(allele_slen + 1, &new_alloc)) {
      return 1;
    }
    memcpyx(new_alloc, newval, allele_slen, '\0');
    newptr = new_alloc;
  }
  const uintptr_t bigstack_end_addr = (uintptr_t)g_bigstack_end;
  const uintptr_t maxdiff = ((uintptr_t)(&(g_one_char_strs[512]))) - bigstack_end_addr;
  // take advantage of unsigned wraparound
  if ((((uintptr_t)(*allele_ptr)) - bigstack_end_addr) >= maxdiff) {
    free(*allele_ptr);
  }
  *allele_ptr = newptr;
  return 0;
}

void cleanup_allele_storage(uint32_t max_allele_slen, uintptr_t allele_storage_entry_ct, const char** allele_storage) {
  // Now doesn't improperly free bigstack allocations (as long as they aren't
  // past g_bigstack_end), and doesn't need to be called at all most of the
  // time.

  // An alternative representation: have a separate bitarray which indicates
  // whether the allele_storage[] element should be interpreted as a heap
  // pointer or an in-place zero-terminated string (i.e. string length can be
  // up to 7 on 64-bit systems).  I expect that to be more efficient for new
  // datasets, but let's get the simple (and 1.9-codebase-compatible)
  // implementation working first, and then benchmark the fancier code later.
  if (allele_storage && (max_allele_slen > 1)) {
    const uintptr_t bigstack_end_addr = (uintptr_t)g_bigstack_end;
    const uintptr_t maxdiff = ((uintptr_t)(&(g_one_char_strs[512]))) - bigstack_end_addr;
    for (uintptr_t idx = 0; idx != allele_storage_entry_ct; ++idx) {
      const char* cur_entry = allele_storage[idx];
      assert(cur_entry);
      // take advantage of unsigned wraparound
      if ((((uintptr_t)cur_entry) - bigstack_end_addr) >= maxdiff) {
        free_const(cur_entry);
      }
    }
  }
}
*/

uint32_t IsCategoricalPhenostr(const char* phenostr_iter) {
  uint32_t first_char_code = ctou32(*phenostr_iter++);
  // allow leading +/-
  if ((first_char_code == 43) || (first_char_code == 45)) {
    first_char_code = ctou32(*phenostr_iter++);
  }
  if (((first_char_code - 48) < 10) || (first_char_code == 44) || (first_char_code < 32)) {
    // the last two conditions are for detecting CSV empty strings
    return 0;
  }
  if (first_char_code == 46) {
    // decimal point.  classify based on whether next character is a digit.
    const uint32_t second_char_code = ctou32(phenostr_iter[0]);
    return ((second_char_code - 48) >= 10);
  }
  // allow any capitalization of "NA"/"nan", but not "inf"
  if ((first_char_code & 0xdf) != 78) {
    return 1;
  }
  const uint32_t second_char_code = ctou32(phenostr_iter[0]);
  if ((second_char_code & 0xdf) != 65) {
    return 1;
  }
  const uint32_t third_char_code = ctou32(phenostr_iter[1]);
  if ((third_char_code & 0xdf) == 78) {
    return (ctou32(phenostr_iter[2]) > ' ');
  }
  return (third_char_code > 32);
}

uint32_t IsCategoricalPhenostrNocsv(const char* phenostr_iter) {
  uint32_t first_char_code = ctou32(*phenostr_iter++);
  // allow leading +/-
  if ((first_char_code == 43) || (first_char_code == 45)) {
    first_char_code = ctou32(*phenostr_iter++);
  }
  if ((first_char_code - 48) < 10) {
    return 0;
  }
  if (first_char_code == 46) {
    // decimal point.  classify based on whether next character is a digit.
    const uint32_t second_char_code = ctou32(phenostr_iter[0]);
    return ((second_char_code - 48) >= 10);
  }
  // allow any capitalization of "NA"/"nan", but not "inf"
  if ((first_char_code & 0xdf) != 78) {
    return 1;
  }
  const uint32_t second_char_code = ctou32(phenostr_iter[0]);
  if ((second_char_code & 0xdf) != 65) {
    return 1;
  }
  const uint32_t third_char_code = ctou32(phenostr_iter[1]);
  if ((third_char_code & 0xdf) == 78) {
    return (ctou32(phenostr_iter[2]) > ' ');
  }
  return (third_char_code > 32);
}

uint32_t FirstCcOrQtPhenoIdx(const PhenoCol* pheno_cols, uint32_t pheno_ct) {
  for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
    if (pheno_cols[pheno_idx].type_code < kPhenoDtypeCat) {
      return pheno_idx;
    }
  }
  return UINT32_MAX;
}

uint32_t IsConstCovar(const PhenoCol* covar_col, const uintptr_t* sample_include, uint32_t raw_sample_ct) {
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  const uintptr_t* nonmiss = covar_col->nonmiss;
  if (covar_col->type_code == kPhenoDtypeQt) {
    const double* covar_vals = covar_col->data.qt;
    for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
      uintptr_t cur_bits = sample_include[widx] & nonmiss[widx];
      if (cur_bits) {
        const double first_covar_val = covar_vals[widx * kBitsPerWord + ctzw(cur_bits)];
        cur_bits &= cur_bits - 1;
        while (1) {
          if (cur_bits) {
            const double* cur_covar_vals = &(covar_vals[widx * kBitsPerWord]);
            do {
              const uint32_t bit_idx = ctzw(cur_bits);
              if (cur_covar_vals[bit_idx] != first_covar_val) {
                return 0;
              }
              cur_bits &= cur_bits - 1;
            } while (cur_bits);
          }
          ++widx;
          if (widx == raw_sample_ctl) {
            return 1;
          }
          cur_bits = sample_include[widx] & nonmiss[widx];
        }
      }
    }
    return 1;
  }
  assert(covar_col->type_code == kPhenoDtypeCat);
  const uint32_t* covar_cats = covar_col->data.cat;
  for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
    uintptr_t cur_bits = sample_include[widx] & nonmiss[widx];
    if (cur_bits) {
      const uint32_t first_covar_cat = covar_cats[widx * kBitsPerWord + ctzw(cur_bits)];
      cur_bits &= cur_bits - 1;
      while (1) {
        if (cur_bits) {
          const uint32_t* cur_covar_cats = &(covar_cats[widx * kBitsPerWord]);
          do {
            const uint32_t bit_idx = ctzw(cur_bits);
            if (cur_covar_cats[bit_idx] != first_covar_cat) {
              return 0;
            }
            cur_bits &= cur_bits - 1;
          } while (cur_bits);
        }
        ++widx;
        if (widx == raw_sample_ctl) {
          return 1;
        }
        cur_bits = sample_include[widx] & nonmiss[widx];
      }
    }
  }
  return 1;
}

uint32_t IdentifyRemainingCats(const uintptr_t* sample_include, const PhenoCol* covar_col, uint32_t sample_ct, uintptr_t* observed_cat_bitarr) {
  // assumes covar_col->type_code == kPhenoTypeCat
  const uint32_t nonnull_cat_ct = covar_col->nonnull_category_ct;
  const uint32_t* covar_cats = covar_col->data.cat;
  const uint32_t word_ct = 1 + (nonnull_cat_ct / kBitsPerWord);
  ZeroWArr(word_ct, observed_cat_bitarr);
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    SetBit(covar_cats[sample_uidx], observed_cat_bitarr);
  }
  return PopcountWords(observed_cat_bitarr, word_ct);
}

// returns index of most common category
uint32_t IdentifyRemainingCatsAndMostCommon(const uintptr_t* sample_include, const PhenoCol* covar_col, uint32_t sample_ct, uintptr_t* observed_cat_bitarr, uint32_t* cat_obs_buf) {
  // assumes covar_col->type_code == kPhenoTypeCat
  const uint32_t nonnull_cat_ct = covar_col->nonnull_category_ct;
  const uint32_t* covar_cats = covar_col->data.cat;
  const uint32_t word_ct = 1 + (nonnull_cat_ct / kBitsPerWord);
  ZeroWArr(word_ct, observed_cat_bitarr);
  ZeroU32Arr(nonnull_cat_ct + 1, cat_obs_buf);
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    cat_obs_buf[covar_cats[sample_uidx]] += 1;
    SetBit(covar_cats[sample_uidx], observed_cat_bitarr);
  }
  if (cat_obs_buf[0]) {
    // don't actually need this for now since we're only calling this with the
    // missing category excluded, but let's be consistent.
    SetBit(0, observed_cat_bitarr);
  }
  uint32_t best_cat_idx = 0;
  uint32_t best_cat_obs_ct = 0;
  for (uint32_t cat_idx = 1; cat_idx <= nonnull_cat_ct; ++cat_idx) {
    const uint32_t cat_obs_ct = cat_obs_buf[cat_idx];
    if (cat_obs_ct) {
      SetBit(cat_idx, observed_cat_bitarr);
      if (cat_obs_ct > best_cat_obs_ct) {
        best_cat_idx = cat_idx;
        best_cat_obs_ct = cat_obs_ct;
      }
    }
  }
  return best_cat_idx;
}

uint32_t GetCatSamples(const uintptr_t* sample_include_base, const PhenoCol* cat_pheno_col, uint32_t raw_sample_ctl, uint32_t sample_ct, uint32_t cat_uidx, uintptr_t* cur_cat_samples) {
  ZeroWArr(raw_sample_ctl, cur_cat_samples);
  const uint32_t* cat_vals = cat_pheno_col->data.cat;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include_base[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include_base, &sample_uidx_base, &cur_bits);
    if (cat_vals[sample_uidx] == cat_uidx) {
      SetBit(sample_uidx, cur_cat_samples);
    }
  }
  return PopcountWords(cur_cat_samples, raw_sample_ctl);
}

uint32_t RemoveExcludedCats(const uint32_t* data_cat, const uintptr_t* cat_keep_bitarr, uint32_t raw_sample_ct, uint32_t in_sample_ct, uintptr_t* sample_include) {
  uintptr_t sample_include_base = 0;
  uintptr_t sample_include_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != in_sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_include_base, &sample_include_bits);
    if (!IsSet(cat_keep_bitarr, data_cat[sample_uidx])) {
      ClearBit(sample_uidx, sample_include);
    }
  }
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  return PopcountWords(sample_include, raw_sample_ctl);
}

void CleanupPhenoCols(uint32_t pheno_ct, PhenoCol* pheno_cols) {
  if (pheno_cols) {
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      vecaligned_free_cond(pheno_cols[pheno_idx].nonmiss);
    }
    free(pheno_cols);
  }
}

PglErr ParseChrRanges(const char* const* argvk, const char* flagname_p, const char* errstr_append, uint32_t param_ct, uint32_t prohibit_extra_chrs, uint32_t xymt_subtract, char range_delim, ChrInfo* cip, uintptr_t* chr_mask) {
  PglErr reterr = kPglRetSuccess;
  {
    const char* cur_arg_ptr = argvk[1];
    char* token_buf = g_textbuf;
    const char* range_end = nullptr;
    uint32_t cur_param_idx = 1;
    uint32_t rs_len = 0;
    uint32_t re_len = 0;
    while (1) {
      const char* range_start;
      if (unlikely(ParseNextRange(argvk, param_ct, range_delim, &cur_param_idx, &cur_arg_ptr, &range_start, &rs_len, &range_end, &re_len))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s parameter '%s'.\n", flagname_p, argvk[cur_param_idx]);
        goto ParseChrRanges_ret_INVALID_CMDLINE_WWA;
      }
      if (!range_start) {
        break;
      }
      // Could avoid this copy, but only at the cost of mutating argv.  Which
      // isn't actually problematic--it's not like we have multiple threads
      // parsing the command line--but it's not like this is a performance
      // bottleneck either.
      // ParseNextRange() prevents buffer overflow.
      memcpyx(token_buf, range_start, rs_len, '\0');
      uint32_t chr_code_start = GetChrCodeRaw(token_buf);
      if (IsI32Neg(chr_code_start)) {
        if (unlikely(prohibit_extra_chrs)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s chromosome code '%s'.\n", flagname_p, token_buf);
          goto ParseChrRanges_ret_INVALID_CMDLINE_WWA;
        }
        if (unlikely(range_end)) {
          goto ParseChrRanges_ret_INVALID_CMDLINE_NONSTD;
        }
        if (unlikely(PushLlStr(token_buf, &(cip->incl_excl_name_stack)))) {
          goto ParseChrRanges_ret_NOMEM;
        }
      } else {
        if (chr_code_start >= kMaxContigs) {
          chr_code_start -= xymt_subtract;
        }
        if (range_end) {
          memcpyx(token_buf, range_end, re_len, '\0');
          uint32_t chr_code_end = GetChrCodeRaw(token_buf);
          if (unlikely(IsI32Neg(chr_code_end))) {
            if (prohibit_extra_chrs) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s chromosome code '%s'.\n", flagname_p, range_end);
              goto ParseChrRanges_ret_INVALID_CMDLINE_WWA;
            }
            goto ParseChrRanges_ret_INVALID_CMDLINE_NONSTD;
          }
          if (unlikely(chr_code_end >= kMaxContigs)) {
            // prohibit stuff like "--chr par1-par2", "--chr x-y", "--chr x-26"
            snprintf(g_logbuf, kLogbufSize, "Error: --%s chromosome code '%s' cannot be the end of a range.\n", flagname_p, range_end);
            goto ParseChrRanges_ret_INVALID_CMDLINE_WWA;
          }
          if (unlikely(chr_code_end <= chr_code_start)) {
            snprintf(g_logbuf, kLogbufSize, "Error: --%s chromosome code '%s' is not greater than '%s'.\n", flagname_p, range_end, range_start);
            goto ParseChrRanges_ret_INVALID_CMDLINE_WWA;
          }
          FillBitsNz(chr_code_start, chr_code_end + 1, chr_mask);
        } else {
          SetBit(chr_code_start, chr_mask);
        }
      }
    }
    // no compelling reason to prohibit "--not-chr ,"
  }
  while (0) {
  ParseChrRanges_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ParseChrRanges_ret_INVALID_CMDLINE_NONSTD:
    logerrputs("Error: Chromosome ranges cannot include nonstandard names.\n");
    reterr = kPglRetInvalidCmdline;
    break;
  ParseChrRanges_ret_INVALID_CMDLINE_WWA:
    WordWrapB(0);
    logerrputsb();
    logerrputs(errstr_append);
    reterr = kPglRetInvalidCmdline;
    break;
  }
  return reterr;
}

uint32_t MultiallelicVariantPresent(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t variant_ct) {
  if (!allele_idx_offsets) {
    return 0;
  }
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
    const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    if (allele_idx_offsets[variant_uidx + 1] != allele_idx_offsets[variant_uidx] + 2) {
      return 1;
    }
  }
  return 0;
}

uint32_t CountBiallelicVariants(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t variant_ct) {
  if (!allele_idx_offsets) {
    return variant_ct;
  }
  uint32_t biallelic_variant_ct = 0;
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
    const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    biallelic_variant_ct += (allele_idx_offsets[variant_uidx + 1] == allele_idx_offsets[variant_uidx] + 2);
  }
  return biallelic_variant_ct;
}

uintptr_t CountAlleles(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t variant_ct) {
  if (!allele_idx_offsets) {
    return variant_ct * 2;
  }
  if (raw_variant_ct == variant_ct) {
    return allele_idx_offsets[raw_variant_ct];
  }
  uintptr_t allele_ct = 0;
  uint32_t variant_uidx = 0;
  while (variant_ct) {
    const uint32_t variant_uidx_start = AdvTo1Bit(variant_include, variant_uidx);
    variant_uidx = AdvBoundedTo0Bit(variant_include, variant_uidx_start + 1, raw_variant_ct);
    allele_ct += allele_idx_offsets[variant_uidx] - allele_idx_offsets[variant_uidx_start];
    variant_ct -= variant_uidx - variant_uidx_start;
  }
  return allele_ct;
}

/*
uintptr_t GetMaxAlleleBlockSize(const uintptr_t* __restrict variant_include, const uintptr_t* __restrict allele_idx_offsets, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t read_block_size) {
  uintptr_t max_allele_block_size = 2 * read_block_size;
  if (!allele_idx_offsets) {
    return max_allele_block_size;
  }
  uint32_t variant_idx = 0;
  uint32_t variant_uidx_end = 0;
  while (1) {
    uint32_t variant_uidx = variant_uidx_end;
    variant_uidx_end += read_block_size;
    if (variant_uidx_end > raw_variant_ct) {
      variant_uidx_end = raw_variant_ct;
      read_block_size = raw_variant_ct - variant_uidx;
    }
    uint32_t cur_variant_ct = PopcountBitRange(variant_include, variant_uidx, variant_uidx_end);
    if (cur_variant_ct) {
      const uint32_t omitted_variant_ct = read_block_size - cur_variant_ct;
      const uintptr_t block_size_ubound = allele_idx_offsets[variant_uidx_end] - allele_idx_offsets[variant_uidx] - (2 * omitted_variant_ct);
      variant_idx += cur_variant_ct;
      if (block_size_ubound > max_allele_block_size) {
        // subtract at beginning of each variant_include bitblock, add at end.
        uintptr_t cur_allele_block_size = 0;
        do {
          const uint32_t cur_bitblock_start = AdvTo1Bit(variant_include, variant_uidx);
          // deliberate underflow
          cur_allele_block_size -= allele_idx_offsets[cur_bitblock_start];
          if (variant_uidx_end - cur_bitblock_start == cur_variant_ct) {
            cur_allele_block_size += allele_idx_offsets[variant_uidx_end];
            break;
          }
          variant_uidx = AdvTo0Bit(variant_include, cur_bitblock_start);
          cur_allele_block_size += allele_idx_offsets[variant_uidx];
          cur_variant_ct -= variant_uidx - cur_bitblock_start;
        } while (cur_variant_ct);
        if (cur_allele_block_size > max_allele_block_size) {
          max_allele_block_size = cur_allele_block_size;
        }
      }
      if (variant_idx == variant_ct) {
        return max_allele_block_size;
      }
    }
  }
}
*/

uintptr_t GetMaxAltAlleleBlockSize(const uintptr_t* __restrict variant_include, const uintptr_t* __restrict allele_idx_offsets, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t read_block_size) {
  uintptr_t max_alt_allele_block_size = read_block_size;
  if (!allele_idx_offsets) {
    return max_alt_allele_block_size;
  }
  uint32_t variant_idx = 0;
  for (uint32_t variant_uidx_end = 0; ; ) {
    uint32_t variant_uidx = variant_uidx_end;
    variant_uidx_end += read_block_size;
    if (variant_uidx_end > raw_variant_ct) {
      variant_uidx_end = raw_variant_ct;
      read_block_size = raw_variant_ct - variant_uidx;
    }
    uint32_t cur_variant_ct = PopcountBitRange(variant_include, variant_uidx, variant_uidx_end);
    if (cur_variant_ct) {
      const uint32_t omitted_variant_ct = read_block_size - cur_variant_ct;
      const uintptr_t cur_ubound = allele_idx_offsets[variant_uidx_end] - allele_idx_offsets[variant_uidx] - omitted_variant_ct - read_block_size;
      variant_idx += cur_variant_ct;
      if (cur_ubound > max_alt_allele_block_size) {
        // subtract at beginning of each variant_include bitblock, add at end.
        // deliberate underflow
        uintptr_t cur_alt_allele_block_size = -S_CAST(uintptr_t, cur_variant_ct);
        do {
          const uint32_t cur_bitblock_start = AdvTo1Bit(variant_include, variant_uidx);
          // deliberate underflow
          cur_alt_allele_block_size -= allele_idx_offsets[cur_bitblock_start];
          if (variant_uidx_end - cur_bitblock_start == cur_variant_ct) {
            cur_alt_allele_block_size += allele_idx_offsets[variant_uidx_end];
            break;
          }
          variant_uidx = AdvTo0Bit(variant_include, cur_bitblock_start);
          cur_alt_allele_block_size += allele_idx_offsets[variant_uidx];
          cur_variant_ct -= variant_uidx - cur_bitblock_start;
        } while (cur_variant_ct);
        if (cur_alt_allele_block_size > max_alt_allele_block_size) {
          max_alt_allele_block_size = cur_alt_allele_block_size;
        }
      }
      if (variant_idx == variant_ct) {
        return max_alt_allele_block_size;
      }
    }
  }
}

uintptr_t GetMhcWordCt(uintptr_t sample_ct) {
  const uintptr_t sample_ctl = BitCtToWordCt(sample_ct);
  const uintptr_t patch_01_max_word_ct = sample_ctl + DivUp(sizeof(AlleleCode) * sample_ct, kBytesPerWord);
  const uintptr_t patch_10_max_word_ct = sample_ctl + DivUp(2 * sizeof(AlleleCode) * sample_ct, kBytesPerWord);
  return RoundUpPow2(patch_01_max_word_ct, kWordsPerVec) + RoundUpPow2(patch_10_max_word_ct, kWordsPerVec);
}

PglErr PgenMtLoadInit(const uintptr_t* variant_include, uint32_t sample_ct, uint32_t variant_ct, uintptr_t bytes_avail, uintptr_t pgr_alloc_cacheline_ct, uintptr_t thread_xalloc_cacheline_ct, uintptr_t per_variant_xalloc_byte_ct, uintptr_t per_alt_allele_xalloc_byte_ct, PgenFileInfo* pgfip, uint32_t* calc_thread_ct_ptr, uintptr_t*** genovecs_ptr, uintptr_t*** mhc_ptr, uintptr_t*** phasepresent_ptr, uintptr_t*** phaseinfo_ptr, uintptr_t*** dosage_present_ptr, Dosage*** dosage_mains_ptr, uintptr_t*** dphase_present_ptr, SDosage*** dphase_delta_ptr, uint32_t* read_block_size_ptr, uintptr_t* max_alt_allele_block_size_ptr, STD_ARRAY_REF(unsigned char*, 2) main_loadbufs, PgenReader*** pgr_pps, uint32_t** read_variant_uidx_starts_ptr) {
  uintptr_t cachelines_avail = bytes_avail / kCacheline;
  uint32_t read_block_size = kPglVblockSize;
  uint64_t multiread_cacheline_ct;
  for (; ; read_block_size /= 2) {
    multiread_cacheline_ct = PgfiMultireadGetCachelineReq(variant_include, pgfip, variant_ct, read_block_size);
    // limit each raw load buffer to 1/4 of remaining workspace
    // if there's an additional per-variant allocation, put it in the same bin
    // as the load buffers
    if ((multiread_cacheline_ct + (S_CAST(uint64_t, per_variant_xalloc_byte_ct) * read_block_size) / kCacheline) * 4 <= cachelines_avail) {
      break;
    }
    // lots of callers require read_block_size to be either raw_variant_ct or a
    // multiple of kBitsPerVec
#ifdef __LP64__
    if (read_block_size <= kBitsPerVec) {
      return kPglRetNomem;
    }
#else
    if (read_block_size <= kCacheline) {
      return kPglRetNomem;
    }
#endif
  }
#ifndef __LP64__
  if (multiread_cacheline_ct > (kMaxBytesPerIO / kCacheline)) {
    return kPglRetNomem;
  }
#endif
  main_loadbufs[0] = S_CAST(unsigned char*, bigstack_alloc_raw(multiread_cacheline_ct * kCacheline));
  main_loadbufs[1] = S_CAST(unsigned char*, bigstack_alloc_raw(multiread_cacheline_ct * kCacheline));
  pgfip->block_base = main_loadbufs[0];
  *read_block_size_ptr = read_block_size;
  cachelines_avail -= 2 * (multiread_cacheline_ct + (S_CAST(uint64_t, per_variant_xalloc_byte_ct) * read_block_size) / kCacheline);
  if (per_alt_allele_xalloc_byte_ct) {
    const uintptr_t max_alt_allele_block_size = GetMaxAltAlleleBlockSize(variant_include, pgfip->allele_idx_offsets, pgfip->raw_variant_ct, variant_ct, read_block_size);

    cachelines_avail -= DivUpU64(S_CAST(uint64_t, per_alt_allele_xalloc_byte_ct) * read_block_size, kCacheline);
    // we assume it's unnecessary to return this if
    // per_alt_allele_xalloc_byte_ct is zero
    *max_alt_allele_block_size_ptr = max_alt_allele_block_size;
  }
  // reduce calc_thread_ct if necessary
  uint32_t calc_thread_ct = *calc_thread_ct_ptr;
  if (calc_thread_ct > read_block_size) {
    calc_thread_ct = read_block_size;
    *calc_thread_ct_ptr = calc_thread_ct;
  }

  // pgr_pps, read_variant_uidx_starts_ptr, (*pgr_pps)[tidx], pgr_alloc;
  //   deliberately a slight overestimate
  const uintptr_t pgr_struct_alloc = RoundUpPow2(sizeof(PgenReader), kCacheline);
  uintptr_t thread_alloc_cacheline_ct = 1 + 1 + (pgr_struct_alloc / kCacheline) + pgr_alloc_cacheline_ct + thread_xalloc_cacheline_ct;

  const uint32_t sample_ctcl2 = NypCtToCachelineCt(sample_ct);
  const uint32_t sample_ctcl = BitCtToCachelineCt(sample_ct);

  // todo: multiallelic dosage
  const uintptr_t dosage_main_cl = DivUp(sample_ct, (kCacheline / sizeof(Dosage)));
  uintptr_t mhc_cl = 0;
  if (genovecs_ptr) {
    thread_alloc_cacheline_ct += 1 + sample_ctcl2;
    if (mhc_ptr) {
      mhc_cl = DivUp(GetMhcWordCt(sample_ct), kWordsPerCacheline);
      thread_alloc_cacheline_ct += 1 + mhc_cl;
    }
    if (phasepresent_ptr) {
      assert(phaseinfo_ptr);
      thread_alloc_cacheline_ct += 2 + sample_ctcl;
    }
    if (dosage_present_ptr) {
      assert(dosage_mains_ptr);
      thread_alloc_cacheline_ct += 2 + sample_ctcl + dosage_main_cl;
      if (dphase_present_ptr) {
        assert(dphase_delta_ptr);
        thread_alloc_cacheline_ct += 2 + sample_ctcl + dosage_main_cl;
      }
    }
  }
  if (thread_alloc_cacheline_ct * calc_thread_ct > cachelines_avail) {
    if (thread_alloc_cacheline_ct > cachelines_avail) {
      return kPglRetNomem;
    }
    calc_thread_ct = cachelines_avail / thread_alloc_cacheline_ct;
    *calc_thread_ct_ptr = calc_thread_ct;
  }

  const uint32_t array_of_ptrs_alloc = RoundUpPow2(calc_thread_ct * sizeof(intptr_t), kCacheline);
  *pgr_pps = S_CAST(PgenReader**, bigstack_alloc_raw(array_of_ptrs_alloc));
  *read_variant_uidx_starts_ptr = S_CAST(uint32_t*, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(int32_t)));
  for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
    (*pgr_pps)[tidx] = S_CAST(PgenReader*, bigstack_alloc_raw(pgr_struct_alloc));
    // PreinitPgr(g_pgr_ptrs[tidx]);
    unsigned char* pgr_alloc = S_CAST(unsigned char*, bigstack_alloc_raw(pgr_alloc_cacheline_ct * kCacheline));

    // shouldn't be possible for this to fail
    PgrInit(nullptr, 0, pgfip, (*pgr_pps)[tidx], pgr_alloc);
  }
  if (genovecs_ptr) {
    *genovecs_ptr = S_CAST(uintptr_t**, bigstack_alloc_raw(array_of_ptrs_alloc));
    if (mhc_ptr) {
      *mhc_ptr = S_CAST(uintptr_t**, bigstack_alloc_raw(array_of_ptrs_alloc));
    }
    if (phasepresent_ptr) {
      *phasepresent_ptr = S_CAST(uintptr_t**, bigstack_alloc_raw(array_of_ptrs_alloc));
      *phaseinfo_ptr = S_CAST(uintptr_t**, bigstack_alloc_raw(array_of_ptrs_alloc));
    }
    if (dosage_present_ptr) {
      *dosage_present_ptr = S_CAST(uintptr_t**, bigstack_alloc_raw(array_of_ptrs_alloc));
      *dosage_mains_ptr = S_CAST(Dosage**, bigstack_alloc_raw(array_of_ptrs_alloc));
      if (dphase_present_ptr) {
        *dphase_present_ptr = S_CAST(uintptr_t**, bigstack_alloc_raw(array_of_ptrs_alloc));
        *dphase_delta_ptr = S_CAST(SDosage**, bigstack_alloc_raw(array_of_ptrs_alloc));
      }
    }
    const uintptr_t genovec_alloc = sample_ctcl2 * kCacheline;
    const uintptr_t bitarray_alloc = sample_ctcl * kCacheline;
    const uintptr_t dosage_main_alloc = dosage_main_cl * kCacheline;
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      (*genovecs_ptr)[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(genovec_alloc));
      if (mhc_ptr) {
        (*mhc_ptr)[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(mhc_cl * kCacheline));
      }
      if (phasepresent_ptr) {
        (*phasepresent_ptr)[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitarray_alloc));
        (*phaseinfo_ptr)[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitarray_alloc));
      }
      if (dosage_present_ptr) {
        (*dosage_present_ptr)[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitarray_alloc));
        (*dosage_mains_ptr)[tidx] = S_CAST(Dosage*, bigstack_alloc_raw(dosage_main_alloc));
        if (dphase_present_ptr) {
          (*dphase_present_ptr)[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(bitarray_alloc));
          (*dphase_delta_ptr)[tidx] = S_CAST(SDosage*, bigstack_alloc_raw(dosage_main_alloc));
        }
      }
    }
  }
  return kPglRetSuccess;
}

uint32_t MultireadNonempty(const uintptr_t* variant_include, const ThreadGroup* tgp, uint32_t raw_variant_ct, uint32_t read_block_size, PgenFileInfo* pgfip, uint32_t* read_block_idxp, PglErr* reterrp) {
  if (IsLastBlock(tgp)) {
    return 0;
  }
  // This condition ensures the PopcountWords() calls are valid.  If it's
  // ever inconvenient, PopcountBitRange() can be used instead.
  assert((!(read_block_size % kBitsPerVec)) || (raw_variant_ct <= read_block_size));

  const uint32_t read_block_sizel = read_block_size / kBitsPerWord;
  uint32_t read_block_idx = *read_block_idxp;
  uint32_t offset = read_block_idx * read_block_size;
  uint32_t cur_read_block_size = read_block_size;
  uint32_t cur_block_write_ct;
  for (; ; ++read_block_idx, offset += read_block_size) {
    if (offset + read_block_size >= raw_variant_ct) {
      cur_read_block_size = raw_variant_ct - offset;
      // trailing bits of variant_include must not be set, otherwise this
      // miscounts
      cur_block_write_ct = PopcountWords(&(variant_include[read_block_idx * read_block_sizel]), BitCtToWordCt(cur_read_block_size));
      assert(cur_block_write_ct);  // otherwise, IsLastBlock should be set
      break;
    }
    cur_block_write_ct = PopcountWords(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
    if (cur_block_write_ct) {
      break;
    }
  }
  *read_block_idxp = read_block_idx;
  *reterrp = PgfiMultiread(variant_include, offset, offset + cur_read_block_size, cur_block_write_ct, pgfip);
  return cur_block_write_ct;
}

void ExpandMhc(uint32_t sample_ct, uintptr_t* mhc, uintptr_t** patch_01_set_ptr, AlleleCode** patch_01_vals_ptr, uintptr_t** patch_10_set_ptr, AlleleCode** patch_10_vals_ptr) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  *patch_01_set_ptr = mhc;
  AlleleCode* patch_01_vals = R_CAST(AlleleCode*, &(mhc[sample_ctl]));
  *patch_01_vals_ptr = patch_01_vals;
  AlleleCode* patch_01_vals_end = &(patch_01_vals[sample_ct]);
  AlignACToVec(&patch_01_vals_end);
  uintptr_t* patch_10_set = R_CAST(uintptr_t*, patch_01_vals_end);
  *patch_10_set_ptr = patch_10_set;
  *patch_10_vals_ptr = R_CAST(AlleleCode*, &(patch_10_set[sample_ctl]));
}

// Returns 2 if resize needed (str_ct == strset_table_size / 2), 1 if OOM
uint32_t StrsetAdd(unsigned char* arena_top, const char* src, uint32_t slen, uint32_t strset_table_size, char** strset, uint32_t* str_ctp, unsigned char** arena_bottom_ptr) {
  for (uint32_t hashval = Hashceil(src, slen, strset_table_size); ; ) {
    char* strset_entry = strset[hashval];
    if (!strset_entry) {
      const uint32_t str_ct = 1 + (*str_ctp);
      if (unlikely(str_ct * 2 > strset_table_size)) {
        return 2;
      }
      *str_ctp = str_ct;
      return S_CAST(uint32_t, StoreStringAtBase(arena_top, src, slen, arena_bottom_ptr, &(strset[hashval])));
    }
    if (strequal_unsafe(strset_entry, src, slen)) {
      return 0;
    }
    if (++hashval == strset_table_size) {
      hashval = 0;
    }
  }
}

uint32_t StrsetAddEnd(unsigned char* arena_bottom, const char* src, uint32_t slen, uint32_t strset_table_size, char** strset, uint32_t* str_ctp, unsigned char** arena_top_ptr) {
  for (uint32_t hashval = Hashceil(src, slen, strset_table_size); ; ) {
    char* strset_entry = strset[hashval];
    if (!strset_entry) {
      const uint32_t str_ct = 1 + (*str_ctp);
      if (unlikely(str_ct * 2 > strset_table_size)) {
        return 2;
      }
      *str_ctp = str_ct;
      return S_CAST(uint32_t, StoreStringAtEnd(arena_bottom, src, slen, arena_top_ptr, &(strset[hashval])));
    }
    if (strequal_unsafe(strset_entry, src, slen)) {
      return 0;
    }
    if (++hashval == strset_table_size) {
      hashval = 0;
    }
  }
}

void RepopulateStrset(char* str_iter, uint32_t str_ct, uint32_t strset_size, char** strset) {
  ZeroPtrArr(strset_size, strset);
  for (uint32_t str_idx = 0; str_idx != str_ct; ++str_idx) {
    char* str_end = strnul(str_iter);
    const uint32_t slen = str_end - str_iter;
    for (uint32_t hashval = Hashceil(str_iter, slen, strset_size); ; ) {
      char* strset_entry = strset[hashval];
      if (!strset_entry) {
        strset[hashval] = str_iter;
        break;
      }
      if (++hashval == strset_size) {
        hashval = 0;
      }
    }
    str_iter = &(str_end[1]);
  }
}

BoolErr StrsetAddResize(unsigned char* arena_top, const char* src, uint32_t slen, uint32_t strset_table_size_max, char** strset, uint32_t* strset_table_sizep, uint32_t* str_ctp, unsigned char** arena_bottom_ptr) {
  uint32_t strset_table_size = *strset_table_sizep;
  uint32_t retval = StrsetAdd(arena_top, src, slen, strset_table_size, strset, str_ctp, arena_bottom_ptr);
  if (!retval) {
    return 0;
  }
  if (unlikely((retval == 1) || (strset_table_size == strset_table_size_max))) {
    return 1;
  }
  const uint32_t new_table_size = MINV(strset_table_size * 2LLU, strset_table_size_max);
  const uint64_t bytes_needed = sizeof(intptr_t) * (new_table_size - strset_table_size);
  if (S_CAST(uintptr_t, arena_top - (*arena_bottom_ptr)) < bytes_needed) {
    return 1;
  }
  unsigned char* old_str_base = R_CAST(unsigned char*, &(strset[strset_table_size]));
  unsigned char* new_str_base = R_CAST(unsigned char*, &(strset[new_table_size]));
  memmove(new_str_base, old_str_base, (*arena_bottom_ptr) - old_str_base);
  *arena_bottom_ptr += bytes_needed;
  RepopulateStrset(R_CAST(char*, new_str_base), *str_ctp, new_table_size, strset);
  *strset_table_sizep = new_table_size;
  return (StrsetAdd(arena_top, src, slen, new_table_size, strset, str_ctp, arena_bottom_ptr) != 0);
}

BoolErr StrsetAddEndResize(unsigned char* arena_bottom, const char* src, uint32_t slen, uint32_t strset_table_size_max, char*** strsetp, uint32_t* strset_table_sizep, uint32_t* str_ctp, unsigned char** arena_top_ptr) {
  char** strset = *strsetp;
  uint32_t strset_table_size = *strset_table_sizep;
  uint32_t retval = StrsetAddEnd(arena_bottom, src, slen, strset_table_size, strset, str_ctp, arena_top_ptr);
  if (!retval) {
    return 0;
  }
  if (unlikely((retval == 1) || (strset_table_size == strset_table_size_max))) {
    return 1;
  }
  const uint32_t new_table_size = MINV(strset_table_size * 2LLU, strset_table_size_max);
  const uint64_t bytes_needed = sizeof(intptr_t) * (new_table_size - strset_table_size);
  unsigned char* old_str_base = *arena_top_ptr;
  if (unlikely(S_CAST(uintptr_t, old_str_base - arena_bottom) < bytes_needed)) {
    return 1;
  }
  unsigned char* new_str_base = old_str_base - bytes_needed;
  memmove(new_str_base, old_str_base, R_CAST(unsigned char*, strset) - old_str_base);
  *arena_top_ptr = new_str_base;
  strset -= new_table_size - strset_table_size;
  *strsetp = strset;
  RepopulateStrset(R_CAST(char*, new_str_base), *str_ctp, new_table_size, strset);
  *strset_table_sizep = new_table_size;
  return (StrsetAddEnd(arena_bottom, src, slen, new_table_size, strset, str_ctp, arena_top_ptr) != 0);
}

PglErr WriteSampleIdsOverride(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* outname, uint32_t sample_ct, SampleIdFlags override_flags) {
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto WriteSampleIdsOverride_ret_OPEN_FAIL;
    }
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    char* write_iter = g_textbuf;
    char* textbuf_flush = &(write_iter[kMaxMediumLine]);
    if (!(override_flags & kfSampleIdNoIdHeader)) {
      *write_iter++ = '#';
      if (override_flags & kfSampleIdFidPresent) {
        write_iter = strcpya_k(write_iter, "FID\t");
      } else {
        // every single row starts with "0\t", so this causes only IIDs to be
        // reported
        sample_ids = &(sample_ids[2]);
      }
      write_iter = strcpya_k(write_iter, "IID");
      if (sids) {
        write_iter = strcpya_k(write_iter, "\tSID");
      }
      AppendBinaryEoln(&write_iter);
    } else {
      if (override_flags & kfSampleIdNoIdHeaderIidOnly) {
        // We've previously verified that all FIDs are in fact '0'.
        sample_ids = &(sample_ids[2]);
      }
      sids = nullptr;
    }
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      write_iter = strcpya(write_iter, &(sample_ids[sample_uidx * max_sample_id_blen]));
      if (sids) {
        *write_iter++ = '\t';
        write_iter = strcpya(write_iter, &(sids[sample_uidx * max_sid_blen]));
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto WriteSampleIdsOverride_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto WriteSampleIdsOverride_ret_WRITE_FAIL;
    }
  }
  while (0) {
  WriteSampleIdsOverride_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  WriteSampleIdsOverride_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  return reterr;
}

uint32_t RealpathIdentical(const char* outname, const char* read_realpath, char* write_realpath_buf) {
#ifdef _WIN32
  const uint32_t fname_slen = GetFullPathName(outname, kPglFnamesize, write_realpath_buf, nullptr);
  return (fname_slen && (fname_slen <= kPglFnamesize) && memequal(read_realpath, write_realpath_buf, fname_slen + 1));
#else
  return (realpath(outname, write_realpath_buf) && strequal_overread(read_realpath, write_realpath_buf));
#endif
}

char* PrintMultiallelicHcAsDs(uint32_t hc1, uint32_t hc2, uint32_t allele_ct, char* start) {
  if (hc1 == kMissingAlleleCode) {
    *start++ = '.';
    return start;
  }
  for (uint32_t uii = 1; uii < hc1; ++uii) {
    start = strcpya_k(start, "0,");
  }
  if (hc1 == hc2) {
    if (hc1) {
      start = strcpya_k(start, "2,");
    }
    for (uint32_t uii = hc1 + 1; uii != allele_ct; ++uii) {
      start = strcpya_k(start, "0,");
    }
    return &(start[-1]);
  }
  if (hc1) {
    start = strcpya_k(start, "1,");
  }
  for (uint32_t uii = hc1 + 1; uii != hc2; ++uii) {
    start = strcpya_k(start, "0,");
  }
  start = strcpya_k(start, "1,");
  for (uint32_t uii = hc2 + 1; uii != allele_ct; ++uii) {
    start = strcpya_k(start, "0,");
  }
  return &(start[-1]);
}

char* PrintMultiallelicHcAsHaploidDs(uint32_t hc1, uint32_t hc2, uint32_t allele_ct, char* start) {
  if (hc1 == kMissingAlleleCode) {
    *start++ = '.';
    return start;
  }
  for (uint32_t uii = 1; uii < hc1; ++uii) {
    start = strcpya_k(start, "0,");
  }
  if (hc1 == hc2) {
    if (hc1) {
      start = strcpya_k(start, "1,");
    }
    for (uint32_t uii = hc1 + 1; uii != allele_ct; ++uii) {
      start = strcpya_k(start, "0,");
    }
    return &(start[-1]);
  }
  if (hc1) {
    start = strcpya_k(start, "0.5,");
  }
  for (uint32_t uii = hc1 + 1; uii != hc2; ++uii) {
    start = strcpya_k(start, "0,");
  }
  start = strcpya_k(start, "0.5,");
  for (uint32_t uii = hc2 + 1; uii != allele_ct; ++uii) {
    start = strcpya_k(start, "0,");
  }
  return &(start[-1]);
}

const char g_vft_names[3][18] = {"extract", "extract-intersect", "exclude"};

void PgenErrPrintEx(const char* file_descrip, uint32_t prepend_lf, PglErr reterr, uint32_t variant_uidx) {
  if (reterr == kPglRetReadFail) {
    if (prepend_lf) {
      logputs("\n");
    }
    logerrprintfww(kErrprintfFread, file_descrip, rstrerror(errno));
  } else if (reterr == kPglRetMalformedInput) {
    if (prepend_lf) {
      logputs("\n");
    }
    if (variant_uidx == UINT32_MAX) {
      logerrprintfww("Error: Failed to unpack variant in %s.\n", file_descrip);
    } else {
      logerrprintfww("Error: Failed to unpack (0-based) variant #%u in %s.\n", variant_uidx, file_descrip);
    }
    logerrputs("You can use --validate to check whether it is malformed.\n* If it is malformed, you probably need to either re-download the file, or\n  address an error in the command that generated the input .pgen.\n* If it appears to be valid, you have probably encountered a plink2 bug.  If\n  you report the error on GitHub or the plink2-users Google group (make sure to\n  include the full .log file in your report), we'll try to address it.\n");
  }
}

void PgenWriteFinishErrPrint(PglErr reterr, char* outname, char* outname_end) {
  // This mirrors include/pgenlib_write.cc PwcFinish().
  if (reterr == kPglRetOpenFail) {
    // Only way this can happen is if .pgen.tmp reopen fails.
    logputs("\n");
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen.tmp");
    logerrprintfww(kErrprintfFopen, outname, strerror(errno));
  } else if (reterr == kPglRetRewindFail) {
    logputs("\n");
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen.tmp");
    logerrprintfww("Error: %s is malformed.\n", outname);
  }
}

// Given <outname>.tmp.pgen and <outname>.tmp.pgen.pgi, this generates
// <outname>.pgen and then deletes the two temporary files.
// TODO: add a third pgenlib_write mode which takes care of this automatically.
// At the minimal cost of one extra malloc/free (for the base filename string),
// we avoid the need to create a temporary index file.
PglErr EmbedPgenIndex(char* outname, char* outname_end) {
  FILE* infile = nullptr;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    strcpy_k(outname_end, ".pgen");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto EmbedPgenIndex_ret_OPEN_FAIL;
    }
    strcpy_k(outname_end, ".tmp.pgen.pgi");
    if (unlikely(fopen_checked(outname, FOPEN_RB, &infile))) {
      goto EmbedPgenIndex_ret_OPEN_FAIL;
    }
    if (unlikely(fseeko(infile, 0, SEEK_END))) {
      goto EmbedPgenIndex_ret_READ_PGI_FAIL;
    }
    const uint64_t pgi_fsize = ftello(infile);
    if (unlikely(pgi_fsize < 20)) {
      goto EmbedPgenIndex_ret_INVALID_INDEX;
    }
    rewind(infile);
    uint32_t variant_ct;
    {
      char startbuf[12];
      if (unlikely(!fread(startbuf, 12, 1, infile))) {
        goto EmbedPgenIndex_ret_READ_PGI_FAIL;
      }
      if (unlikely(!memequal_sk(startbuf, "l\x1b\x30"))) {
        goto EmbedPgenIndex_ret_INVALID_INDEX;
      }
      memcpy(&variant_ct, &(startbuf[3]), 4);
      // May allow variant_ct == 0 later.
      if (unlikely((!variant_ct) || (variant_ct > kPglMaxVariantCt))) {
        goto EmbedPgenIndex_ret_INVALID_INDEX;
      }
      startbuf[2] = 0x10;
      if (unlikely(!fwrite_unlocked(startbuf, 12, 1, outfile))) {
        goto EmbedPgenIndex_ret_WRITE_FAIL;
      }
    }
    {
      const uint32_t vblock_ct = DivUp(variant_ct, kPglVblockSize);
      // The next (vblock_ct * 8) bytes of the .pgen.pgi file are a
      // random-access index, reporting the byte positions that variants 0,
      // 65536, 131072, etc. start at in the .pgen.
      // We must increment these positions by (pgi_fsize - 3), since the .pgen
      // without an embedded index has only 3 header bytes.
      if (unlikely(pgi_fsize < 12 + vblock_ct * 8)) {
        goto EmbedPgenIndex_ret_INVALID_INDEX;
      }
      const uint64_t offset_incr = pgi_fsize - 3;
      uint32_t remaining_vblock_ct = vblock_ct;
      uint64_t offset_buf[kPglFwriteBlockSize / 8];
      while (remaining_vblock_ct > kPglFwriteBlockSize / 8) {
        if (unlikely(!fread(offset_buf, kPglFwriteBlockSize, 1, infile))) {
          goto EmbedPgenIndex_ret_READ_PGI_FAIL;
        }
        for (uint32_t uii = 0; uii != kPglFwriteBlockSize / 8; ++uii) {
          offset_buf[uii] += offset_incr;
        }
        if (unlikely(!fwrite_unlocked(offset_buf, kPglFwriteBlockSize, 1, outfile))) {
          goto EmbedPgenIndex_ret_WRITE_FAIL;
        }
        remaining_vblock_ct -= kPglFwriteBlockSize / 8;
      }
      if (unlikely(!fread(offset_buf, remaining_vblock_ct * 8, 1, infile))) {
        goto EmbedPgenIndex_ret_READ_PGI_FAIL;
      }
      for (uint32_t uii = 0; uii != remaining_vblock_ct; ++uii) {
        offset_buf[uii] += offset_incr;
      }
      if (unlikely(!fwrite_unlocked(offset_buf, remaining_vblock_ct * 8, 1, outfile))) {
        goto EmbedPgenIndex_ret_WRITE_FAIL;
      }
    }
    while (1) {
      unsigned char copybuf[kPglFwriteBlockSize];
      uintptr_t nbyte = fread(copybuf, 1, kPglFwriteBlockSize, infile);
      if (!nbyte) {
        break;
      }
      if (unlikely(!fwrite_unlocked(copybuf, nbyte, 1, outfile))) {
        goto EmbedPgenIndex_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_null(&infile))) {
      goto EmbedPgenIndex_ret_READ_PGI_FAIL;
    }
    strcpy_k(outname_end, ".tmp.pgen");
    if (unlikely(fopen_checked(outname, FOPEN_RB, &infile))) {
      goto EmbedPgenIndex_ret_OPEN_FAIL;
    }
    unsigned char copybuf[kPglFwriteBlockSize + 3];
    uintptr_t nbyte = fread(copybuf, 1, kPglFwriteBlockSize + 3, infile);
    if (unlikely((nbyte <= 3) || (!memequal_sk(copybuf, "l\x1b\x20")))) {
      logerrprintfww("Error: %s does not appear to be a .pgen with an external index file.\n", outname);
      goto EmbedPgenIndex_ret_MALFORMED_INPUT;
    }
    nbyte -= 3;
    if (unlikely(!fwrite_unlocked(&(copybuf[3]), nbyte, 1, outfile))) {
      goto EmbedPgenIndex_ret_WRITE_FAIL;
    }
    while (1) {
      nbyte = fread(copybuf, 1, kPglFwriteBlockSize, infile);
      if (!nbyte) {
        break;
      }
      if (unlikely(!fwrite_unlocked(copybuf, nbyte, 1, outfile))) {
        goto EmbedPgenIndex_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_null(&infile))) {
      goto EmbedPgenIndex_ret_READ_PGEN_FAIL;
    }
    if (unlikely(fclose_null(&outfile))) {
      goto EmbedPgenIndex_ret_WRITE_FAIL;
    }
    if (unlikely(unlink(outname))) {
      logerrprintfww("Error: Failed to delete %s : %s.\n", outname, strerror(errno));
      goto EmbedPgenIndex_ret_WRITE_FAIL;
    }
    strcpy_k(outname_end, ".tmp.pgen.pgi");
    if (unlikely(unlink(outname))) {
      logerrprintfww("Error: Failed to delete %s : %s.\n", outname, strerror(errno));
      goto EmbedPgenIndex_ret_WRITE_FAIL;
    }
  }
  while (0) {
  EmbedPgenIndex_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  EmbedPgenIndex_ret_READ_PGI_FAIL:
    // strcpy_k(outname_end, ".tmp.pgen.pgi");
    logerrprintfww(kErrprintfFread, outname, rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  EmbedPgenIndex_ret_READ_PGEN_FAIL:
    // strcpy_k(outname_end, ".tmp.pgen");
    logerrprintfww(kErrprintfFread, outname, rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  EmbedPgenIndex_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  EmbedPgenIndex_ret_INVALID_INDEX:
    // strcpy_k(outname_end, ".tmp.pgen.pgi");
    logerrprintfww("Error: %s does not appear to be a .pgen index file.\n", outname);
  EmbedPgenIndex_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  fclose_cond(outfile);
  fclose_cond(infile);
  return reterr;
}

CXXCONST_CP HkvlineValEnd(const char* val_start) {
  if (val_start[0] != '"') {
    const char* val_end = strchrnul2_n(val_start, ',', '>');
    if (unlikely(*val_end == '\n')) {
      return nullptr;
    }
    return S_CAST(CXXCONST_CP, val_end);
  }
  const char* val_iter = &(val_start[1]);
  while (1) {
    val_iter = strchrnul2_n(val_iter, '\\', '"');
    if (*val_iter == '"') {
      break;
    }
    if (unlikely((*val_iter != '\\') || (val_iter[1] < 32))) {
      return nullptr;
    }
    val_iter = &(val_iter[2]);
  }
  const char* val_end = &(val_iter[1]);
  if (likely((*val_end == ',') || (*val_end == '>'))) {
    return S_CAST(CXXCONST_CP, val_end);
  }
  return nullptr;
}

IntErr HkvlineFind2K(const char* dict_iter, const char* key, uint32_t key_slen, const char** val_ptr, uint32_t* val_slenp) {
  while (1) {
    if (dict_iter[-1] == '>') {
      return 1;
    }
    const char* val_prestart = strchrnul_n(dict_iter, '=');
    if (unlikely(*val_prestart == '\n')) {
      return 2;
    }
    const char* val_start = &(val_prestart[1]);
    const char* val_end = HkvlineValEnd(val_start);
    if (unlikely(val_end == nullptr)) {
      return 2;
    }
    if ((S_CAST(uintptr_t, val_prestart - dict_iter) == key_slen) && memequal(dict_iter, key, key_slen)) {
      *val_ptr = val_start;
      *val_slenp = val_end - val_start;
      return 0;
    }
    dict_iter = &(val_end[1]);
  }
}

BoolErr HkvlineNumTypeK(const char* dict_iter, const char** numstr_ptr, uint32_t* num_slenp, const char** typestr_ptr, uint32_t* type_slenp) {
  uint32_t need_numstr = 1;
  uint32_t need_typestr = 1;
  while (1) {
    // Note that this is true on the first iteration when ID is the only key.
    if (unlikely(dict_iter[-1] == '>')) {
      return 1;
    }
    const char* val_prestart = strchrnul_n(dict_iter, '=');
    if (unlikely(*val_prestart == '\n')) {
      return 1;
    }
    const char* val_start = &(val_prestart[1]);
    const char* val_end = HkvlineValEnd(val_start);
    if (unlikely(val_end == nullptr)) {
      return 1;
    }
    const uint32_t key_slen = val_prestart - dict_iter;
    if (need_numstr && strequal_k(dict_iter, "Number", key_slen)) {
      *numstr_ptr = val_start;
      *num_slenp = val_end - val_start;
      if (!need_typestr) {
        return 0;
      }
      need_numstr = 0;
    } else if (need_typestr && strequal_k(dict_iter, "Type", key_slen)) {
      *typestr_ptr = val_start;
      *type_slenp = val_end - val_start;
      if (!need_numstr) {
        return 0;
      }
      need_typestr = 0;
    }
    ++dict_iter;
  }
}

PglErr HkvlineForceIdFirst(uintptr_t workspace_size, char* hline_kv_start, void* workspace) {
  if (memequal_sk(hline_kv_start, "ID=")) {
    return kPglRetSuccess;
  }
  char* idval;
  uint32_t id_slen;
  if (unlikely(HkvlineFind(hline_kv_start, "ID", &idval, &id_slen))) {
    return kPglRetMalformedInput;
  }
  const uintptr_t tmp_nbyte = id_slen + (3 * k1LU);
  if (tmp_nbyte > workspace_size) {
    return kPglRetNomem;
  }
  char* idkey_start = &(idval[-3]);
  char* workspace_c = S_CAST(char*, workspace);
  memcpy(workspace_c, idkey_start, tmp_nbyte);
  memmove(&(hline_kv_start[tmp_nbyte + 1]), hline_kv_start, &(idkey_start[-1]) - hline_kv_start);
  memcpy(hline_kv_start, workspace_c, tmp_nbyte);
  hline_kv_start[tmp_nbyte] = ',';
  return kPglRetSuccess;
}

PglErr NsortDedupAndWrite(const char* outname, uintptr_t str_ct, uint32_t max_slen, uint32_t output_zst, uint32_t max_thread_ct, const char** strptr_arr, uintptr_t* nwrite_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    StrptrArrNsort(str_ct, strptr_arr);
    const uintptr_t overflow_buf_size = kCompressStreamBlock + strlen(EOLN_STR) + max_slen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto NsortDedupAndWrite_ret_1;
    }
    const char* prev_id = nullptr;
    uintptr_t dup_ct = 0;
    uint32_t prev_slen = 0;
    for (uintptr_t ulii = 0; ulii != str_ct; ++ulii) {
      const char* cur_id = strptr_arr[ulii];
      const uint32_t cur_slen = strlen(cur_id);
      if ((cur_slen == prev_slen) && memequal(prev_id, cur_id, prev_slen)) {
        ++dup_ct;
        continue;
      }
      cswritep = memcpya(cswritep, cur_id, cur_slen);
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto NsortDedupAndWrite_ret_WRITE_FAIL;
      }
      prev_id = cur_id;
      prev_slen = cur_slen;
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto NsortDedupAndWrite_ret_WRITE_FAIL;
    }
    *nwrite_ptr = str_ct - dup_ct;
  }
  while (0) {
  NsortDedupAndWrite_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 NsortDedupAndWrite_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct VariantLastAlidxMakerStruct {
  const uintptr_t* allele_idx_offsets;
  uint32_t raw_variant_ct;

  // Bit k is set iff it's equal to (allele_idx_offsets[x] - 1) for some x.
  // This has the property that PopcountBitRange(0, allele_idx) = variant_uidx.
  uintptr_t* variant_last_alidxs;
  uint32_t* variant_last_alidxs_cumulative_popcounts;
} VariantLastAlidxMaker;

void VariantLastAlidxMakerMain(uint32_t tidx, uint32_t thread_ct, VariantLastAlidxMaker* ctx) {
  const uintptr_t* allele_idx_ends = &(ctx->allele_idx_offsets[1]);
  const uint32_t raw_variant_ct = ctx->raw_variant_ct;
  const uintptr_t raw_allele_ct = ctx->allele_idx_offsets[raw_variant_ct];
  const uintptr_t raw_allele_ctl = BitCtToWordCt(raw_allele_ct);
  const uintptr_t cacheline_ct = DivUp(raw_allele_ctl, kWordsPerCacheline);

  const uintptr_t widx_start = kWordsPerCacheline * ((cacheline_ct * S_CAST(uint64_t, tidx)) / thread_ct);
  const uintptr_t widx_stop = (tidx + 1 == thread_ct)? (raw_allele_ctl - 1) : (kWordsPerCacheline * ((cacheline_ct * (S_CAST(uint64_t, tidx) + 1)) / thread_ct));
  uintptr_t* alidxs = ctx->variant_last_alidxs;
  uint32_t* alidxs_cumulative_popcounts = ctx->variant_last_alidxs_cumulative_popcounts;

  uint32_t cur_vidx = 0;
  if (tidx) {
    // How many variants end before or at allele_idx == widx_start *
    // kBitsPerWord?
    cur_vidx = UpperBoundNonemptyW(allele_idx_ends, raw_variant_ct, widx_start * kBitsPerWord);
  }
  // This is >= widx_start * kBitsPerWord.
  uintptr_t cur_allele_idx = allele_idx_ends[cur_vidx] - 1;

  uintptr_t cur_allele_widx = cur_allele_idx / kBitsPerWord;
  for (uintptr_t widx = widx_start; widx != widx_stop; ++widx) {
    alidxs_cumulative_popcounts[widx] = cur_vidx;
    uintptr_t cur_word = 0;
    if (cur_allele_widx == widx) {
      uint32_t shift_ct = cur_allele_idx % kBitsPerWord;
      // fast path when multiallelic variants are rare.
      const uint32_t speculate_ct = (kBitsPerWord + 1 - shift_ct) / 2;
      const uint32_t speculative_vidx = cur_vidx + speculate_ct;
      if ((speculative_vidx < raw_variant_ct) && (allele_idx_ends[speculative_vidx] == cur_allele_idx + 1 + speculate_ct * 2)) {
        // next speculate_ct variants must all be biallelic
        cur_word = kMask5555 << shift_ct;
        cur_vidx = speculative_vidx;
        cur_allele_idx += speculate_ct * 2;
        ++cur_allele_widx;
      } else {
        while (1) {
          cur_word |= k1LU << shift_ct;
          ++cur_vidx;
          cur_allele_idx = allele_idx_ends[cur_vidx] - 1;
          cur_allele_widx = cur_allele_idx / kBitsPerWord;
          if (cur_allele_widx != widx) {
            break;
          }
          shift_ct = cur_allele_idx % kBitsPerWord;
        }
      }
    }
    alidxs[widx] = cur_word;
  }
  if (tidx + 1 == thread_ct) {
    alidxs_cumulative_popcounts[widx_stop] = cur_vidx;
    uintptr_t cur_word = 0;
    while (1) {
      cur_word |= k1LU << (cur_allele_idx % kBitsPerWord);
      ++cur_vidx;
      if (cur_vidx == raw_variant_ct) {
        break;
      }
      cur_allele_idx = allele_idx_ends[cur_vidx] - 1;
    }
    alidxs[widx_stop] = cur_word;
  }
}

THREAD_FUNC_DECL VariantLastAlidxMakerThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  VariantLastAlidxMaker* ctx = S_CAST(VariantLastAlidxMaker*, arg->sharedp->context);

  const uint32_t thread_ct = GetThreadCt(arg->sharedp) + 1;
  VariantLastAlidxMakerMain(tidx, thread_ct, ctx);
  THREAD_RETURN;
}

PglErr AllocAndFillVariantLastAlidxs(const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t max_thread_ct, const uintptr_t** variant_last_alidxsp, const uint32_t** variant_last_alidxs_cumulative_popcountsp) {
  if (!allele_idx_offsets) {
    *variant_last_alidxsp = nullptr;
    *variant_last_alidxs_cumulative_popcountsp = nullptr;
    return kPglRetSuccess;
  }
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  VariantLastAlidxMaker ctx;
  {
    const uintptr_t raw_allele_ct = allele_idx_offsets[raw_variant_ct];
    const uintptr_t raw_allele_ctl = BitCtToWordCt(raw_allele_ct);
    if (unlikely(bigstack_alloc_w(raw_allele_ctl, &ctx.variant_last_alidxs) ||
                 bigstack_alloc_u32(raw_allele_ctl, &ctx.variant_last_alidxs_cumulative_popcounts))) {
      goto AllocAndFillVariantLastAlidxs_ret_NOMEM;
    }
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.raw_variant_ct = raw_variant_ct;
    uint32_t thread_ct = raw_variant_ct / 65536;
    if (!thread_ct) {
      thread_ct = 1;
    } else if (thread_ct > max_thread_ct) {
      thread_ct = max_thread_ct;
    }
    if (unlikely(SetThreadCt0(thread_ct - 1, &tg))) {
      goto AllocAndFillVariantLastAlidxs_ret_NOMEM;
    }
    if (thread_ct > 1) {
      SetThreadFuncAndData(VariantLastAlidxMakerThread, &ctx, &tg);
      DeclareLastThreadBlock(&tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto AllocAndFillVariantLastAlidxs_ret_THREAD_CREATE_FAIL;
      }
    }
    VariantLastAlidxMakerMain(thread_ct - 1, thread_ct, &ctx);
    JoinThreads0(&tg);
    *variant_last_alidxsp = ctx.variant_last_alidxs;
    *variant_last_alidxs_cumulative_popcountsp = ctx.variant_last_alidxs_cumulative_popcounts;
  }
  while (0) {
  AllocAndFillVariantLastAlidxs_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  AllocAndFillVariantLastAlidxs_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  return reterr;
}

void GetBpWindow(const uintptr_t* variant_include, const uint32_t* variant_bps, uint32_t variant_uidx, uint32_t bp_radius, uint32_t* start_vidxp, uint32_t* end_vidxp) {
  assert(IsSet(variant_include, variant_uidx));
  const uint32_t center_bp = variant_bps[variant_uidx];
  {
    const uint32_t search_start_vidx = *start_vidxp;
    if (variant_uidx > search_start_vidx) {
      uint32_t first_bp = 0;
      if (center_bp > bp_radius) {
        first_bp = center_bp - bp_radius;
      }
      const uint32_t start_vidx = LowerBoundConstrainedNonemptyU32(variant_bps, search_start_vidx, variant_uidx, first_bp);
      *start_vidxp = AdvTo1Bit(variant_include, start_vidx);
    }
  }
  const uint32_t search_start_vidx = variant_uidx + 1;
  const uint32_t search_end_vidx = *end_vidxp;
  if (search_start_vidx < search_end_vidx) {
    const uint32_t end_vidx = LowerBoundConstrainedNonemptyU32(variant_bps, search_start_vidx, search_end_vidx, center_bp + bp_radius + 1);
    *end_vidxp = 1 + FindLast1BitBefore(variant_include, end_vidx);
  }
}

PgenGlobalFlags GflagsVfilter(const uintptr_t* variant_include, const unsigned char* vrtypes, uint32_t raw_variant_ct, PgenGlobalFlags input_gflags) {
  PgenGlobalFlags read_gflags = kfPgenGlobal0;
  if (!vrtypes) {
    return input_gflags;
  }
  const uintptr_t* vrtypes_alias = R_CAST(const uintptr_t*, vrtypes);
  const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
  uint32_t mask_multiply = ((input_gflags & kfPgenGlobalHardcallPhasePresent)? 0x10 : 0) + ((input_gflags & kfPgenGlobalDosagePresent)? 0x60 : 0) + ((input_gflags & kfPgenGlobalDosagePhasePresent)? 0x80 : 0);
  uintptr_t vrtypes_or = 0;
  // todo: try changing loop to be vec-based, use movemask to extract
  // information from vrtypes in 64-bit cases
  for (uint32_t widx = 0; widx != raw_variant_ctl; ++widx) {
    uintptr_t cur_variant_include_word = variant_include[widx];
    if (cur_variant_include_word) {
      // bugfix (20 Aug 2018): this needs to advance on every variant_include
      // word, not just the nonzero ones
      const uintptr_t* cur_vrtypes = &(vrtypes_alias[8 * widx]);
#ifdef __LP64__
      for (uint32_t vi_byte_idx = 0; vi_byte_idx != 8; ++vi_byte_idx) {
#  ifdef USE_AVX2
        // this doesn't seem to be much faster than non-AVX2 code on my Mac...
        // inverse-movemask shouldn't be better than regular movemask here
        const uintptr_t cur_mask = _pdep_u64(cur_variant_include_word, kMask0101);
#  else
        // this operation maps binary hgfedcba to h0000000g0000000f...
        //                                        ^       ^       ^
        //                                        |       |       |
        //                                       56      48      40
        // 1. (cur_variant_include_word & 0xfe) gives us hgfedcb0;
        //    necessary to avoid carryover.
        // 2. multiply by the number with bits 7, 14, 21, ..., 49 set, to
        //    get hgfedcbhgfedcbhgf...
        //        ^       ^       ^
        //        |       |       |
        //       56      48      40
        // 3. mask out all but bits 8, 16, 24, ..., 56
        // todo: test if this actually beats the per-character loop...
        const uintptr_t cur_mask = (((cur_variant_include_word & 0xfe) * 0x2040810204080LLU) & kMask0101) | (cur_variant_include_word & 1);
#  endif
        vrtypes_or |= cur_vrtypes[vi_byte_idx] & (cur_mask * mask_multiply);
        cur_variant_include_word >>= 8;
      }
#else
      for (uint32_t vi_hexa_idx = 0; vi_hexa_idx != 8; ++vi_hexa_idx) {
        // dcba -> d0000000c0000000b0000000a
        const uintptr_t cur_mask = ((cur_variant_include_word & 0xf) * 0x204081) & kMask0101;
        vrtypes_or |= cur_vrtypes[vi_hexa_idx] & (cur_mask * mask_multiply);
        cur_variant_include_word >>= 4;
      }
#endif
      if (vrtypes_or) {
        // bugfix (8 Oct 2017): forgot to multiply by kMask0101
        if (vrtypes_or & (0x10 * kMask0101)) {
          read_gflags |= kfPgenGlobalHardcallPhasePresent;
          mask_multiply -= 0x10;
        }
        if (vrtypes_or & (0x60 * kMask0101)) {
          read_gflags |= kfPgenGlobalDosagePresent;
          mask_multiply -= 0x60;
        }
        if (vrtypes_or & (0x80 * kMask0101)) {
          read_gflags |= kfPgenGlobalDosagePhasePresent;
          mask_multiply -= 0x80;
        }
        if (!mask_multiply) {
          return read_gflags;
        }
      }
    }
  }
  return read_gflags;
}

uint32_t TriangleDivide(int64_t cur_prod_x2, int32_t modif) {
  // return smallest integer vv for which (vv * (vv + modif)) is no smaller
  // than cur_prod_x2, and neither term in the product is negative.
  int64_t vv;
  if (cur_prod_x2 == 0) {
    if (modif < 0) {
      return -modif;
    }
    return 0;
  }
  vv = S_CAST(int64_t, sqrt(S_CAST(double, cur_prod_x2)));
  while ((vv - 1) * (vv + modif - 1) >= cur_prod_x2) {
    vv--;
  }
  while (vv * (vv + modif) < cur_prod_x2) {
    vv++;
  }
  return vv;
}

void ParallelBounds(uint32_t ct, int32_t start, uint32_t parallel_idx, uint32_t parallel_tot, int32_t* __restrict bound_start_ptr, int32_t* __restrict bound_end_ptr) {
  int32_t modif = 1 - start * 2;
  int64_t ct_tot = S_CAST(int64_t, ct) * (ct + modif);
  *bound_start_ptr = TriangleDivide((ct_tot * parallel_idx) / parallel_tot, modif);
  *bound_end_ptr = TriangleDivide((ct_tot * (parallel_idx + 1)) / parallel_tot, modif);
}

void AppendZerotabsUnsafe(uint32_t zerotab_ct, char** cswritep_ptr) {
  // (roughly same performance as creating a zero-tab constant buffer in
  // advance)
  const uint32_t word_ct = DivUp(zerotab_ct, kBytesPerWord / 2);
  const uintptr_t zerotab_word = 0x930 * kMask0001;
  char* cswritep = *cswritep_ptr;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    CAppendW(zerotab_word, &cswritep);
  }
  *cswritep_ptr = &((*cswritep_ptr)[2 * zerotab_ct]);
}

typedef struct InitAllelePermuteCtxStruct {
  const uintptr_t* allele_idx_offsets;
  uint32_t raw_variant_ct;

  AlleleCode* allele_permute;
} InitAllelePermuteCtx;

void InitAllelePermuteMain(uintptr_t tidx, uintptr_t thread_ct, InitAllelePermuteCtx* ctx) {
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uintptr_t* allele_idx_offset_ends = &(allele_idx_offsets[1]);
  const uint32_t raw_variant_ct = ctx->raw_variant_ct;
  const uint32_t variant_uidx_start = (raw_variant_ct * S_CAST(uint64_t, tidx)) / thread_ct;
  const uint32_t variant_uidx_end = (raw_variant_ct * (S_CAST(uint64_t, tidx) + 1)) / thread_ct;
  uintptr_t allele_idx_offset_start = allele_idx_offsets[variant_uidx_start];
  AlleleCode* allele_permute_iter = &(ctx->allele_permute[allele_idx_offset_start]);
  // could fill with zeroes first?
  for (uint32_t vidx = variant_uidx_start; vidx != variant_uidx_end; ++vidx) {
    const uintptr_t allele_idx_offset_end = allele_idx_offset_ends[vidx];
    const uint32_t allele_ct = allele_idx_offset_end - allele_idx_offset_start;
    for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
      *allele_permute_iter++ = aidx;
    }
    allele_idx_offset_start = allele_idx_offset_end;
  }
}

THREAD_FUNC_DECL InitAllelePermuteThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  InitAllelePermuteCtx* ctx = S_CAST(InitAllelePermuteCtx*, arg->sharedp->context);

  const uint32_t thread_ct = GetThreadCt(arg->sharedp) + 1;
  InitAllelePermuteMain(tidx, thread_ct, ctx);
  THREAD_RETURN;
}

PglErr InitAllelePermuteUnsafe(const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t max_thread_ct, AlleleCode* allele_permute) {
  if (!allele_idx_offsets) {
    uintptr_t* allele_permute_ul = R_CAST(uintptr_t*, allele_permute);
    const uintptr_t allele_code_range_size = k1LU << (8 * sizeof(AlleleCode));
    const uintptr_t fill_word = ((~k0LU) / ((allele_code_range_size - 1) * (allele_code_range_size + 1))) * allele_code_range_size;
    const uintptr_t word_ct = DivUp(raw_variant_ct, sizeof(intptr_t) / (sizeof(AlleleCode) * 2));
    for (uintptr_t widx = 0; widx != word_ct; ++widx) {
      allele_permute_ul[widx] = fill_word;
    }
    return kPglRetSuccess;
  }
  assert(allele_idx_offsets[raw_variant_ct] > raw_variant_ct * 2);
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  InitAllelePermuteCtx ctx;
  {
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.raw_variant_ct = raw_variant_ct;
    ctx.allele_permute = allele_permute;
    uint32_t thread_ct = raw_variant_ct / 65536;
    if (!thread_ct) {
      thread_ct = 1;
    } else if (thread_ct > max_thread_ct) {
      thread_ct = max_thread_ct;
    }
    if (unlikely(SetThreadCt0(thread_ct - 1, &tg))) {
      goto InitAllelePermuteUnsafe_ret_NOMEM;
    }
    if (thread_ct > 1) {
      SetThreadFuncAndData(InitAllelePermuteThread, &ctx, &tg);
      DeclareLastThreadBlock(&tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto InitAllelePermuteUnsafe_ret_THREAD_CREATE_FAIL;
      }
    }
    InitAllelePermuteMain(thread_ct - 1, thread_ct, &ctx);
    JoinThreads0(&tg);
  }
  while (0) {
  InitAllelePermuteUnsafe_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  InitAllelePermuteUnsafe_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  return reterr;
}

typedef struct TrimAllelePermuteCtxStruct {
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* allele_presents;

  AlleleCode* allele_permute;

  uint32_t* variant_uidx_starts;

  uint32_t variant_ct;
} TrimAllelePermuteCtx;

void TrimAllelePermuteMain(uint32_t tidx, uint32_t thread_ct, TrimAllelePermuteCtx* ctx) {
  const uintptr_t* variant_include = ctx->variant_include;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uintptr_t* allele_presents = ctx->allele_presents;
  const uint32_t variant_ct = ctx->variant_ct;
  const uint32_t variant_uidx_start = ctx->variant_uidx_starts[tidx];
  const uint32_t variant_idx_start = (variant_ct * S_CAST(uint64_t, tidx)) / thread_ct;
  const uint32_t variant_idx_end = (variant_ct * (S_CAST(uint64_t, tidx) + 1)) / thread_ct;
  AlleleCode* allele_permute = ctx->allele_permute;
  uintptr_t cur_allele_ct = 2;
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, variant_uidx_start, &variant_uidx_base, &cur_bits);
  for (uint32_t variant_idx = variant_idx_start; variant_idx != variant_idx_end; ++variant_idx) {
    const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    uintptr_t allele_idx_offset_base = 2 * variant_uidx;
    if (allele_idx_offsets) {
      allele_idx_offset_base = allele_idx_offsets[variant_uidx];
      cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
    }
    uintptr_t read_aidx = 1;
    for (; read_aidx != cur_allele_ct; ++read_aidx) {
      if (!IsSet(allele_presents, allele_idx_offset_base + read_aidx)) {
        break;
      }
    }
    if (read_aidx != cur_allele_ct) {
      AlleleCode* cur_allele_permute = &(allele_permute[allele_idx_offset_base]);
      uintptr_t write_aidx = read_aidx;
      ++read_aidx;
      for (; read_aidx != cur_allele_ct; ++read_aidx) {
        if (IsSet(allele_presents, allele_idx_offset_base + read_aidx)) {
          cur_allele_permute[write_aidx] = cur_allele_permute[read_aidx];
          ++write_aidx;
        }
      }
      for (; write_aidx != cur_allele_ct; ++write_aidx) {
        cur_allele_permute[write_aidx] = kMissingAlleleCode;
      }
    }
  }
}

THREAD_FUNC_DECL TrimAllelePermuteThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  const uint32_t thread_ct = GetThreadCt(arg->sharedp) + 1;
  TrimAllelePermuteCtx* ctx = S_CAST(TrimAllelePermuteCtx*, arg->sharedp->context);
  TrimAllelePermuteMain(tidx, thread_ct, ctx);
  THREAD_RETURN;
}

PglErr TrimAllelePermute(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const uintptr_t* allele_presents, uint32_t variant_ct, uint32_t max_thread_ct, AlleleCode* allele_permute) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  TrimAllelePermuteCtx ctx;
  {
    const uint32_t thread_ct = ClipU32(variant_ct / 65536, 1, max_thread_ct);
    if (unlikely(SetThreadCt0(thread_ct - 1, &tg) ||
                 bigstack_alloc_u32(thread_ct, &ctx.variant_uidx_starts))) {
      goto TrimAllelePermute_ret_NOMEM;
    }
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.allele_presents = allele_presents;
    ctx.allele_permute = allele_permute;
    ctx.variant_ct = variant_ct;
    FillU32SubsetStarts(variant_include, thread_ct, 0, variant_ct, ctx.variant_uidx_starts);

    if (thread_ct > 1) {
      SetThreadFuncAndData(TrimAllelePermuteThread, &ctx, &tg);
      DeclareLastThreadBlock(&tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto TrimAllelePermute_ret_THREAD_CREATE_FAIL;
      }
    }
    TrimAllelePermuteMain(thread_ct - 1, thread_ct, &ctx);
    JoinThreads0(&tg);
  }
  while (0) {
  TrimAllelePermute_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TrimAllelePermute_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  BigstackReset(bigstack_mark);
  return reterr;
}



// dupstore hash table (due to position in plink2 order of operations), but
// duplicated (in main dataset) ID not allowed.
//
// repeated entries in --flip file *are* allowed, those aren't dangerous in the
// same way that e.g. taking a '.' variant ID literally would be.
//
// returns UINT32_MAX on success, variant_uidx in duplicate-ID set on failure.
//
// could rename if any other operation wants this functionality, though this
// function probably doesn't belong in a more-central location (due to mismatch
// between hash table type and usage).
uint32_t LoadNondup2ProcessTokens(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* shard_start, const char* shard_end, uint32_t variant_id_htable_size, uint32_t max_variant_id_slen, uintptr_t* already_seen) {
  const char* shard_iter = shard_start;
  while (1) {
    shard_iter = FirstPostspaceBounded(shard_iter, shard_end);
    if (shard_iter == shard_end) {
      return UINT32_MAX;
    }
    const char* token_end = CurTokenEnd(shard_iter);
    uint32_t cur_llidx;
    uint32_t variant_uidx = VariantIdDupHtableFind(shard_iter, variant_ids, variant_id_htable, htable_dup_base, token_end - shard_iter, variant_id_htable_size, max_variant_id_slen, &cur_llidx);
    shard_iter = token_end;
    if (variant_uidx == UINT32_MAX) {
      continue;
    }
    if (IsSet(already_seen, variant_uidx)) {
      continue;
    }
    if (cur_llidx != UINT32_MAX) {
      return variant_uidx;
    }
    SetBit(variant_uidx, already_seen);
  }
}

CONSTI32(kMaxLoadTokensNondup2Threads, 8);

typedef struct LoadTokensNondup2CtxStruct {
  const char* const* variant_ids;
  const uint32_t* variant_id_htable;
  const uint32_t* htable_dup_base;
  uintptr_t variant_id_htable_size;
  uint32_t max_variant_id_slen;

  char* shard_boundaries[kMaxLoadTokensNondup2Threads + 1];
  uintptr_t* already_seens[kMaxLoadTokensNondup2Threads];
  uint32_t fail_variant_uidxs[kMaxLoadTokensNondup2Threads];
} LoadTokensNondup2Ctx;

THREAD_FUNC_DECL LoadTokensNondup2Thread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx_p1 = arg->tidx + 1;
  LoadTokensNondup2Ctx* ctx = S_CAST(LoadTokensNondup2Ctx*, arg->sharedp->context);

  const char* const* variant_ids = ctx->variant_ids;
  const uint32_t* variant_id_htable = ctx->variant_id_htable;
  const uint32_t* htable_dup_base = ctx->htable_dup_base;
  const uintptr_t variant_id_htable_size = ctx->variant_id_htable_size;
  const uint32_t max_variant_id_slen = ctx->max_variant_id_slen;
  uintptr_t* already_seen = ctx->already_seens[tidx_p1];
  do {
    ctx->fail_variant_uidxs[tidx_p1] = LoadNondup2ProcessTokens(variant_ids, variant_id_htable, htable_dup_base, ctx->shard_boundaries[tidx_p1], ctx->shard_boundaries[tidx_p1 + 1], variant_id_htable_size, max_variant_id_slen, already_seen);

  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr LoadTokensNondup2(const char* fname, const uintptr_t* variant_include, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* flagname_p, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t max_thread_ct, uintptr_t* loaded_variant_set) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  TokenStream tks;
  ThreadGroup tg;
  PreinitTokenStream(&tks);
  PreinitThreads(&tg);
  LoadTokensNondup2Ctx ctx;
  {
    const uint32_t calc_thread_ct_m1 = MINV(max_thread_ct, kMaxLoadTokensNondup2Threads) - 1;
    if (unlikely(SetThreadCt0(calc_thread_ct_m1, &tg))) {
      goto LoadTokensNondup2_ret_NOMEM;
    }
    uint32_t decompress_thread_ct = 1;
    if (max_thread_ct > calc_thread_ct_m1 + 2) {
      decompress_thread_ct = max_thread_ct - calc_thread_ct_m1 - 1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    ZeroWArr(raw_variant_ctl, loaded_variant_set);
    ctx.already_seens[0] = loaded_variant_set;
    for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
      if (unlikely(bigstack_calloc_w(raw_variant_ctl, &(ctx.already_seens[tidx])))) {
        goto LoadTokensNondup2_ret_NOMEM;
      }
    }
    if (calc_thread_ct_m1) {
      ctx.variant_ids = variant_ids;
      ctx.variant_id_htable = variant_id_htable;
      ctx.htable_dup_base = htable_dup_base;
      ctx.variant_id_htable_size = variant_id_htable_size;
      ctx.max_variant_id_slen = max_variant_id_slen;
      SetThreadFuncAndData(LoadTokensNondup2Thread, &ctx, &tg);
    }
    reterr = InitTokenStream(fname, decompress_thread_ct, &tks);
    if (unlikely(reterr)) {
      goto LoadTokensNondup2_ret_TKSTREAM_FAIL;
    }
    while (1) {
      reterr = TksNext(&tks, calc_thread_ct_m1 + 1, ctx.shard_boundaries);
      if (reterr) {
        break;
      }
      if (calc_thread_ct_m1) {
        if (unlikely(SpawnThreads(&tg))) {
          goto LoadTokensNondup2_ret_THREAD_CREATE_FAIL;
        }
      }
      ctx.fail_variant_uidxs[0] = LoadNondup2ProcessTokens(variant_ids, variant_id_htable, htable_dup_base, ctx.shard_boundaries[0], ctx.shard_boundaries[1], variant_id_htable_size, max_variant_id_slen, ctx.already_seens[0]);
      JoinThreads0(&tg);
      for (uint32_t tidx = 0; tidx <= calc_thread_ct_m1; ++tidx) {
        const uint32_t fail_variant_uidx = ctx.fail_variant_uidxs[tidx];
        if (unlikely(fail_variant_uidx != UINT32_MAX)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Variant '%s' in --%s file appears multiple times in dataset.\n", variant_ids[fail_variant_uidx], flagname_p);
          goto LoadTokensNondup2_ret_INCONSISTENT_INPUT_WW;
        }
      }
    }
    if (unlikely(reterr != kPglRetEof)) {
      goto LoadTokensNondup2_ret_TKSTREAM_FAIL;
    }
    reterr = kPglRetSuccess;
    for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
      BitvecOr(ctx.already_seens[tidx], raw_variant_ctl, loaded_variant_set);
    }
    // possible for some variants to have been excluded after hash table was
    // constructed
    BitvecAnd(variant_include, raw_variant_ctl, loaded_variant_set);
  }
  while (0) {
  LoadTokensNondup2_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadTokensNondup2_ret_TKSTREAM_FAIL:
    TokenStreamErrPrint(fname, &tks);
    break;
  LoadTokensNondup2_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  LoadTokensNondup2_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  CleanupTokenStream2(fname, &tks, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

// dupflag case, except we want (possibly post-permutation) variant_idx instead
// of variant_uidx
uint32_t LoadNondupReindexProcessTokens(const uintptr_t* variant_include, const uint32_t* variant_include_cumulative_popcounts, const uint32_t* old_variant_uidx_to_new, const char* const* variant_ids, const uint32_t* variant_id_htable, const char* shard_start, const char* shard_end, uint32_t variant_id_htable_size, uint32_t max_variant_id_slen, uintptr_t* already_seen) {
  const char* shard_iter = shard_start;
  while (1) {
    shard_iter = FirstPostspaceBounded(shard_iter, shard_end);
    if (shard_iter == shard_end) {
      return UINT32_MAX;
    }
    const char* token_end = CurTokenEnd(shard_iter);
    uint32_t variant_uidx = VariantIdDupflagHtableFind(shard_iter, variant_ids, variant_id_htable, token_end - shard_iter, variant_id_htable_size, max_variant_id_slen);
    shard_iter = token_end;
    if (variant_uidx >> 31) {
      if (unlikely(variant_uidx != UINT32_MAX)) {
        return variant_uidx & 0x7fffffff;
      }
      continue;
    }
    uint32_t variant_idx;
    if (old_variant_uidx_to_new) {
      variant_idx = old_variant_uidx_to_new[variant_uidx];
    } else {
      // could add a just-set-variant_uidx branch
      variant_idx = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx);
    }
    SetBit(variant_idx, already_seen);
  }
}

CONSTI32(kMaxLoadTokensNondupReindexThreads, 8);

typedef struct LoadTokensNondupReindexCtxStruct {
  const uintptr_t* variant_include;
  const uint32_t* variant_include_cumulative_popcounts;
  const uint32_t* old_variant_uidx_to_new;
  const char* const* variant_ids;
  const uint32_t* variant_id_htable;
  uintptr_t variant_id_htable_size;
  uint32_t max_variant_id_slen;

  char* shard_boundaries[kMaxLoadTokensNondupReindexThreads + 1];
  uintptr_t* already_seens[kMaxLoadTokensNondupReindexThreads];
  uint32_t fail_variant_uidxs[kMaxLoadTokensNondupReindexThreads];
} LoadTokensNondupReindexCtx;

THREAD_FUNC_DECL LoadTokensNondupReindexThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx_p1 = arg->tidx + 1;
  LoadTokensNondupReindexCtx* ctx = S_CAST(LoadTokensNondupReindexCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const uint32_t* variant_include_cumulative_popcounts = ctx->variant_include_cumulative_popcounts;
  const uint32_t* old_variant_uidx_to_new = ctx->old_variant_uidx_to_new;
  const char* const* variant_ids = ctx->variant_ids;
  const uint32_t* variant_id_htable = ctx->variant_id_htable;
  const uintptr_t variant_id_htable_size = ctx->variant_id_htable_size;
  const uint32_t max_variant_id_slen = ctx->max_variant_id_slen;
  uintptr_t* already_seen = ctx->already_seens[tidx_p1];
  do {
    ctx->fail_variant_uidxs[tidx_p1] = LoadNondupReindexProcessTokens(variant_include, variant_include_cumulative_popcounts, old_variant_uidx_to_new, variant_ids, variant_id_htable, ctx->shard_boundaries[tidx_p1], ctx->shard_boundaries[tidx_p1 + 1], variant_id_htable_size, max_variant_id_slen, already_seen);

  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr LoadTokensNondupReindex(const char* fname, const uintptr_t* variant_include, const uint32_t* variant_include_cumulative_popcounts, const uint32_t* old_variant_uidx_to_new, const char* const* variant_ids, const char* flagname_p, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t max_thread_ct, uintptr_t* loaded_variant_set) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  TokenStream tks;
  ThreadGroup tg;
  PreinitTokenStream(&tks);
  PreinitThreads(&tg);
  LoadTokensNondupReindexCtx ctx;
  {
    const uint32_t calc_thread_ct_m1 = MINV(max_thread_ct, kMaxLoadTokensNondupReindexThreads) - 1;
    if (unlikely(SetThreadCt0(calc_thread_ct_m1, &tg))) {
      goto LoadTokensNondupReindex_ret_NOMEM;
    }
    uint32_t decompress_thread_ct = 1;
    if (max_thread_ct > calc_thread_ct_m1 + 2) {
      decompress_thread_ct = max_thread_ct - calc_thread_ct_m1 - 1;
    }
    const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
    ZeroWArr(variant_ctl, loaded_variant_set);
    ctx.already_seens[0] = loaded_variant_set;
    for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
      if (unlikely(bigstack_calloc_w(variant_ctl, &(ctx.already_seens[tidx])))) {
        goto LoadTokensNondupReindex_ret_NOMEM;
      }
    }

    reterr = InitTokenStream(fname, decompress_thread_ct, &tks);
    if (unlikely(reterr)) {
      goto LoadTokensNondupReindex_ret_TKSTREAM_FAIL;
    }

    uint32_t* variant_id_htable;
    uint32_t variant_id_htable_size = 0;
    reterr = AllocAndPopulateIdHtableMt(variant_include, variant_ids, variant_ct, 0, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
    if (unlikely(reterr)) {
      goto LoadTokensNondupReindex_ret_1;
    }
    if (calc_thread_ct_m1) {
      ctx.variant_include = variant_include;
      ctx.variant_include_cumulative_popcounts = variant_include_cumulative_popcounts;
      ctx.old_variant_uidx_to_new = old_variant_uidx_to_new;
      ctx.variant_ids = variant_ids;
      ctx.variant_id_htable = variant_id_htable;
      ctx.variant_id_htable_size = variant_id_htable_size;
      ctx.max_variant_id_slen = max_variant_id_slen;
      SetThreadFuncAndData(LoadTokensNondupReindexThread, &ctx, &tg);
    }

    while (1) {
      reterr = TksNext(&tks, calc_thread_ct_m1 + 1, ctx.shard_boundaries);
      if (reterr) {
        break;
      }
      if (calc_thread_ct_m1) {
        if (unlikely(SpawnThreads(&tg))) {
          goto LoadTokensNondupReindex_ret_THREAD_CREATE_FAIL;
        }
      }
      ctx.fail_variant_uidxs[0] = LoadNondupReindexProcessTokens(variant_include, variant_include_cumulative_popcounts, old_variant_uidx_to_new, variant_ids, variant_id_htable, ctx.shard_boundaries[0], ctx.shard_boundaries[1], variant_id_htable_size, max_variant_id_slen, ctx.already_seens[0]);
      JoinThreads0(&tg);
      for (uint32_t tidx = 0; tidx <= calc_thread_ct_m1; ++tidx) {
        const uint32_t fail_variant_uidx = ctx.fail_variant_uidxs[tidx];
        if (unlikely(fail_variant_uidx != UINT32_MAX)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Variant '%s' in --%s file appears multiple times in dataset.\n", variant_ids[fail_variant_uidx], flagname_p);
          goto LoadTokensNondupReindex_ret_INCONSISTENT_INPUT_WW;
        }
      }
    }
    if (unlikely(reterr != kPglRetEof)) {
      goto LoadTokensNondupReindex_ret_TKSTREAM_FAIL;
    }
    reterr = kPglRetSuccess;
    for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
      BitvecOr(ctx.already_seens[tidx], variant_ctl, loaded_variant_set);
    }
  }
  while (0) {
  LoadTokensNondupReindex_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadTokensNondupReindex_ret_TKSTREAM_FAIL:
    TokenStreamErrPrint(fname, &tks);
    break;
  LoadTokensNondupReindex_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  LoadTokensNondupReindex_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 LoadTokensNondupReindex_ret_1:
  CleanupThreads(&tg);
  CleanupTokenStream2(fname, &tks, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
