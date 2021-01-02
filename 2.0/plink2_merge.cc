// This file is part of PLINK 2.00, copyright (C) 2005-2021 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "include/pgenlib_write.h"
#include "plink2_compress_stream.h"
#include "plink2_merge.h"
#include "plink2_psam.h"
#include "plink2_pvar.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitPmerge(PmergeInfo* pmerge_info_ptr) {
  pmerge_info_ptr->flags = kfPmerge0;
  pmerge_info_ptr->list_mode = kPmergeListModePfile;
  pmerge_info_ptr->merge_mode = kMergeModeNmMatch;
  pmerge_info_ptr->merge_pheno_mode = kMergePhenoModeNmMatch;
  pmerge_info_ptr->merge_xheader_mode = kMergeXheaderModeFirst;
  pmerge_info_ptr->merge_qual_mode = kMergeQualInfoModeNmFirst;
  pmerge_info_ptr->merge_filter_mode = kMergeFilterModeNmFirst;
  pmerge_info_ptr->merge_info_mode = kMergeQualInfoModeNmFirst;
  pmerge_info_ptr->merge_info_sort = kfSortNone;
  pmerge_info_ptr->max_allele_ct = 0;
  pmerge_info_ptr->pgen_fname = nullptr;
  pmerge_info_ptr->pvar_fname = nullptr;
  pmerge_info_ptr->psam_fname = nullptr;
  pmerge_info_ptr->list_fname = nullptr;
}

void CleanupPmerge(PmergeInfo* pmerge_info_ptr) {
  free_cond(pmerge_info_ptr->pgen_fname);
  free_cond(pmerge_info_ptr->pvar_fname);
  free_cond(pmerge_info_ptr->psam_fname);
  free_cond(pmerge_info_ptr->list_fname);
}

void InitPgenDiff(PgenDiffInfo* pgen_diff_info_ptr) {
  pgen_diff_info_ptr->flags = kfPgenDiff0;
  pgen_diff_info_ptr->dosage_hap_tol = kDosageMissing;
  pgen_diff_info_ptr->pgen_fname = nullptr;
  pgen_diff_info_ptr->pvar_fname = nullptr;
  pgen_diff_info_ptr->psam_fname = nullptr;
}

void CleanupPgenDiff(PgenDiffInfo* pgen_diff_info_ptr) {
  free_cond(pgen_diff_info_ptr->pgen_fname);
  free_cond(pgen_diff_info_ptr->pvar_fname);
  free_cond(pgen_diff_info_ptr->psam_fname);
}

typedef struct PmergeInputFilesetStruct {
  char* pgen_fname;
  char* pvar_fname;
  char* psam_fname;
} PmergeInputFileset;

PglErr Pmerge(__attribute__((unused)) const PmergeInfo* pmip, __attribute__((unused)) const char* sample_sort_fname, __attribute__((unused)) MiscFlags misc_flags, __attribute__((unused)) SortFlags sample_sort_flags, __attribute__((unused)) FamCol fam_cols, __attribute__((unused)) uint32_t max_thread_ct, __attribute__((unused)) char* pgenname, __attribute__((unused)) char* psamname, __attribute__((unused)) char* pvarname, __attribute__((unused)) char* outname, __attribute__((unused)) char* outname_end, __attribute__((unused)) ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    // 1. Construct/load fileset list.
    // 2. Merge .psam files.  Simplify merge-mode to 'first' if each sample
    //    appears in only one input fileset.
    // 3. Global .pvar scan, to determine merged header, chromosome sort order,
    //    and track first/last position in each fileset.
    // 4. If filesets cover disjoint positions, handle this as a concatenation
    //    job (or error out on --variant-inner-join).
    // 5. Otherwise, perform general-purpose incremental merge.
    // const PmergeFlags flags = pmip->flags;

    logerrputs("Error: --pmerge[-list] is under development.\n");
    reterr = kPglRetNotYetSupported;
    goto Pmerge_ret_1;
  }
  while (0) {

  }
 Pmerge_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

uint32_t DuplicateAllelePresent(const AlleleCode* remap, uint32_t remap_len, uint32_t merged_allele_ctl, uintptr_t* remap_seen) {
  ZeroWArr(merged_allele_ctl, remap_seen);
  for (uint32_t allele_idx = 0; allele_idx != remap_len; ++allele_idx) {
    const AlleleCode ac = remap[allele_idx];
    if (ac == kMissingAlleleCode) {
      continue;
    }
    if (IsSet(remap_seen, ac)) {
      return 1;
    }
    SetBit(ac, remap_seen);
  }
  return 0;
}

typedef struct PgenDiffGtEntryStruct {
  DoubleAlleleCode dac1;
  DoubleAlleleCode dac2;
} PgenDiffGtEntry;

static_assert(sizeof(Dosage) == 2, "PgenDiff must be updated.");
PglErr PgenDiff(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const PgenDiffInfo* pdip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t max_allele_ct1, uint32_t max_allele_slen, uint32_t max_thread_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  char* cswritep = nullptr;
  uintptr_t psam_line_idx = 0;
  uintptr_t pvar_line_idx = 0;
  TextStream psam_txs;
  PreinitTextStream(&psam_txs);
  TextStream pvar_txs;
  PreinitTextStream(&pvar_txs);
  PgenFileInfo pgfi2;
  PgenReader simple_pgr2;
  PreinitPgfi(&pgfi2);
  PreinitPgr(&simple_pgr2);
  CompressStreamState css;
  PreinitCstream(&css);
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* sample_include;
    if (unlikely(bigstack_calloc_w(raw_sample_ctl, &sample_include))) {
      goto PgenDiff_ret_NOMEM;
    }
    const uint32_t sex_needed = XymtIsNonempty(variant_include, cip, kChrOffsetX) || XymtIsNonempty(variant_include, cip, kChrOffsetY);
    uint32_t raw_sample_ct2 = 0;
    uint32_t* sample1_idx_to_2 = nullptr;
    uint32_t* sample_include_cumulative_popcounts;
    uint32_t* sample_idx_to_uidx;
    uintptr_t* sample_include2;
    uint32_t* sample_include2_cumulative_popcounts;
    uint32_t sample_ct;
    {
      uint32_t max_line_blen;
      if (unlikely(StandardizeMaxLineBlen(bigstack_left() / 4, &max_line_blen))) {
        goto PgenDiff_ret_NOMEM;
      }
      reterr = InitTextStreamEx(pdip->psam_fname, 1, kMaxLongLine, max_line_blen, MAXV(max_thread_ct - 1, 1), &psam_txs);
      if (unlikely(reterr)) {
        goto PgenDiff_ret_PSAM_TSTREAM_FAIL;
      }
      char* line_start;
      XidMode xid_mode;
      reterr = LoadXidHeader("pgen-diff .psam", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeaderFixedWidth : kfXidHeaderFixedWidthIgnoreSid, &psam_line_idx, &psam_txs, &xid_mode, &line_start);
      if (unlikely(reterr)) {
        if (reterr == kPglRetEof) {
          logerrputs("Error: Empty --pgen-diff .psam file.\n");
          goto PgenDiff_ret_MALFORMED_INPUT;
        }
        goto PgenDiff_ret_PSAM_TSTREAM_XID_FAIL;
      }
      char* sorted_xidbox;
      uint32_t* xid_map;
      uintptr_t max_xid_blen;
      reterr = SortedXidboxInitAllocEnd(orig_sample_include, siip, orig_sample_ct, 0, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
      if (unlikely(reterr)) {
        goto PgenDiff_ret_1;
      }
      uint32_t* sample1_uidx_to_2;
      char* idbuf;
      if (unlikely(bigstack_end_alloc_u32(raw_sample_ct, &sample1_uidx_to_2) ||
                   bigstack_end_alloc_c(max_xid_blen, &idbuf))) {
        goto PgenDiff_ret_NOMEM;
      }

      uint32_t postid_sex_col_idx = 0;
      if (*line_start == '#') {
        const uint32_t id_col_ct = GetXidColCt(xid_mode);
        char* token_end = CurTokenEnd(NextTokenMult0(line_start, id_col_ct - 1));
        for (uint32_t postid_col_idx = 1; ; ++postid_col_idx) {
          char* token_start = FirstNonTspace(token_end);
          if (IsEolnKns(*token_start)) {
            break;
          }
          token_end = CurTokenEnd(token_start);
          const uint32_t token_slen = token_end - token_start;
          if (strequal_k(token_start, "SEX", token_slen)) {
            if (unlikely(postid_sex_col_idx)) {
              logerrputs("Error: Multiple sex columns in --pgen-diff .psam file.\n");
              goto PgenDiff_ret_MALFORMED_INPUT;
            }
            postid_sex_col_idx = postid_col_idx;
          }
        }
        if (unlikely(!postid_sex_col_idx)) {
          logerrputs("Error: No sex column in --pgen-diff .psam file.\n");
          goto PgenDiff_ret_MALFORMED_INPUT;
        }
        ++psam_line_idx;
        line_start = TextGet(&psam_txs);
      } else {
        // default: FID IID PAT MAT SEX PHENO1
        postid_sex_col_idx = 3;
      }
      if (!sex_needed) {
        postid_sex_col_idx = 0;
      }
      const uintptr_t first_payload_line_idx = psam_line_idx;
      uint32_t matched_sample_ct = 0;
      for (; line_start; ++psam_line_idx, line_start = TextGet(&psam_txs)) {
        const char* token_iter = line_start;
        uint32_t sample_uidx;
        if (SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &token_iter, &sample_uidx, idbuf)) {
          if (unlikely(!token_iter)) {
            goto PgenDiff_ret_PSAM_MISSING_TOKENS;
          }
          continue;
        }
        if (unlikely(IsSet(sample_include, sample_uidx))) {
          TabsToSpaces(idbuf);
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID '%s' in --pgen-diff .psam file.\n", idbuf);
          goto PgenDiff_ret_MALFORMED_INPUT_WW;
        }
        SetBit(sample_uidx, sample_include);
        sample1_uidx_to_2[sample_uidx] = psam_line_idx - first_payload_line_idx;
        ++matched_sample_ct;
        if (postid_sex_col_idx) {
          token_iter = NextTokenMult(token_iter, postid_sex_col_idx);
          if (unlikely(!token_iter)) {
            goto PgenDiff_ret_PSAM_MISSING_TOKENS;
          }
          const uint32_t token_slen = strlen_se(token_iter);
          uint32_t cur_sex_code = 0;
          if (token_slen == 1) {
            const unsigned char sex_ucc = token_iter[0];
            const unsigned char sex_ucc_upcase = sex_ucc & 0xdfU;
            if ((sex_ucc == '1') || (sex_ucc_upcase == 'M')) {
              cur_sex_code = 1;
            } else if ((sex_ucc == '2') || (sex_ucc_upcase == 'F')) {
              cur_sex_code = 2;
            }
          }
          if (!cur_sex_code) {
            snprintf(g_logbuf, kLogbufSize, "Error: Missing sex code on line %" PRIuPTR " of --pgen-diff .psam file; this is prohibited when chrX or chrY is in the comparison.\n", psam_line_idx);
            goto PgenDiff_ret_INCONSISTENT_INPUT_WW;
          }
          if ((!IsSet(sex_nm, sample_uidx)) || (IsSet(sex_male, sample_uidx) != 2 - cur_sex_code)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Mismatching sex code on line %" PRIuPTR " of --pgen-diff .psam file; this is prohibited when chrX or chrY is in the comparison.\n", psam_line_idx);
            goto PgenDiff_ret_INCONSISTENT_INPUT_WW;
          }
        }
      }
      if (unlikely(TextStreamErrcode2(&psam_txs, &reterr))) {
        goto PgenDiff_ret_PSAM_TSTREAM_FAIL;
      }
      if (unlikely(CleanupTextStream2("--pgen-diff .psam file", &psam_txs, &reterr))) {
        goto PgenDiff_ret_1;
      }
      sample_ct = PopcountWords(sample_include, raw_sample_ctl);
      if (unlikely(!sample_ct)) {
        logerrputs("Error: No matching samples in --pgen-diff .psam file.\n");
        goto PgenDiff_ret_INCONSISTENT_INPUT;
      }
      const uintptr_t raw_sample_ct2_ul = psam_line_idx - first_payload_line_idx;
#ifdef __LP64__
      if (raw_sample_ct2_ul > 0x7ffffffe) {
        logerrputs("Error: Too many samples in --pgen-diff .psam file (max 2^31 - 2).\n");
        goto PgenDiff_ret_MALFORMED_INPUT;
      }
#endif
      raw_sample_ct2 = raw_sample_ct2_ul;
      const uint32_t raw_sample_ct2l = BitCtToWordCt(raw_sample_ct2);
      if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                   bigstack_alloc_u32(sample_ct, &sample_idx_to_uidx) ||
                   bigstack_calloc_w(raw_sample_ct2l, &sample_include2) ||
                   bigstack_alloc_u32(raw_sample_ct2l, &sample_include2_cumulative_popcounts))) {
        goto PgenDiff_ret_NOMEM;
      }
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = sample_include[0];
      uint32_t prev_uidx2 = 0; // ok for this to start at 0 instead of -1
      uint32_t is_increasing = 1;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uint32_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
        sample_idx_to_uidx[sample_idx] = sample_uidx;
        const uint32_t sample_uidx2 = sample1_uidx_to_2[sample_uidx];
        SetBit(sample_uidx2, sample_include2);
        if (is_increasing) {
          is_increasing = (sample_uidx2 >= prev_uidx2);
          prev_uidx2 = sample_uidx2;
        }
      }
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      FillCumulativePopcounts(sample_include2, raw_sample_ct2l, sample_include2_cumulative_popcounts);
      if (!is_increasing) {
        if (unlikely(bigstack_alloc_u32(sample_ct, &sample1_idx_to_2))) {
          goto PgenDiff_ret_NOMEM;
        }
        sample_uidx_base = 0;
        cur_bits = sample_include[0];
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          const uint32_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
          const uint32_t sample_uidx2 = sample1_uidx_to_2[sample_uidx];
          const uint32_t sample_idx2 = RawToSubsettedPos(sample_include2, sample_include2_cumulative_popcounts, sample_uidx2);
          sample1_idx_to_2[sample_idx] = sample_idx2;
        }
      }
      BigstackEndReset(bigstack_end_mark);
    }
    uintptr_t* allele_idx_offsets2 = nullptr;
    uint32_t raw_variant_ct2 = 0;
    uint32_t max_observed_pvar_line_blen = 0;
    uint32_t max_allele_slen2 = 1;
    uint32_t max_allele_ct2 = 2;
    reterr = LoadAlleleIdxOffsetsFromPvar(pdip->pvar_fname, "--pgen-diff .pvar file", max_thread_ct, &raw_variant_ct2, &max_allele_slen2, &max_observed_pvar_line_blen, &allele_idx_offsets2, &max_allele_ct2);
    if (unlikely(reterr)) {
      goto PgenDiff_ret_1;
    }
    if (max_allele_slen2 > max_allele_slen) {
      max_allele_slen = max_allele_slen2;
    }

    PgenHeaderCtrl header_ctrl;
    uintptr_t cur_alloc_cacheline_ct;
    reterr = PgfiInitPhase1(pdip->pgen_fname, raw_variant_ct2, raw_sample_ct2, 0, &header_ctrl, &pgfi2, &cur_alloc_cacheline_ct, g_logbuf);
    if (unlikely(reterr)) {
      WordWrapB(0);
      logerrputsb();
      goto PgenDiff_ret_1;
    }
    pgfi2.allele_idx_offsets = allele_idx_offsets2;
    pgfi2.max_allele_ct = max_allele_ct2;
    unsigned char* pgfi_alloc;
    if (unlikely(bigstack_alloc_uc(cur_alloc_cacheline_ct * kCacheline, &pgfi_alloc))) {
      goto PgenDiff_ret_NOMEM;
    }
    const uint32_t all_nonref2 = ((header_ctrl & 192) == 128);
    const uint32_t raw_variant_ctl2 = BitCtToWordCt(raw_variant_ct2);
    uintptr_t* nonref_flags2 = nullptr;
    if ((header_ctrl & 192) == 192) {
      if (unlikely(bigstack_alloc_w(raw_variant_ctl2, &nonref_flags2))) {
        goto PgenDiff_ret_NOMEM;
      }
      pgfi2.nonref_flags = nonref_flags2;
    }
    uintptr_t pgr_alloc_cacheline_ct;
    uint32_t max_vrec_width;
    reterr = PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct2, &max_vrec_width, &pgfi2, pgfi_alloc, &pgr_alloc_cacheline_ct, g_logbuf);
    if (unlikely(reterr)) {
      WordWrapB(0);
      logerrputsb();
      goto PgenDiff_ret_1;
    }
    if (unlikely((!allele_idx_offsets2) && (pgfi2.gflags & kfPgenGlobalMultiallelicHardcallFound))) {
      logerrputs("Error: --pgen-diff .pgen file contains multiallelic variants, while .pvar does\nnot.\n");
      goto PgenDiff_ret_INCONSISTENT_INPUT;
    }
    unsigned char* simple_pgr_alloc;
    if (unlikely(bigstack_alloc_uc((pgr_alloc_cacheline_ct + DivUp(max_vrec_width, kCacheline)) * kCacheline, &simple_pgr_alloc))) {
      goto PgenDiff_ret_NOMEM;
    }
    reterr = PgrInit(pdip->pgen_fname, max_vrec_width, &pgfi2, &simple_pgr2, simple_pgr_alloc);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, pdip->pgen_fname, strerror(errno));
      } else {
        assert(reterr == kPglRetReadFail);
        logerrprintfww(kErrprintfFread, pdip->pgen_fname, rstrerror(errno));
      }
      goto PgenDiff_ret_1;
    }
    PgrSampleSubsetIndex pssi1;
    PgrSampleSubsetIndex pssi2;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi1);
    PgrSetSampleSubsetIndex(sample_include2_cumulative_popcounts, &simple_pgr2, &pssi2);
    const PgenDiffFlags flags = pdip->flags;
    const uint32_t dosage_hap_tol = pdip->dosage_hap_tol;

    // these values are unimportant if dosage_hap_tol == kDosageMissing
    uint32_t dosage_sex_tols[2];
    dosage_sex_tols[0] = dosage_hap_tol / 2;
    dosage_sex_tols[1] = dosage_hap_tol;

    const uint32_t dosage_reported = (dosage_hap_tol != kDosageMissing);
    const uint32_t dosage_needed = ((PgrGetGflags(simple_pgrp) | PgrGetGflags(&simple_pgr2)) & kfPgenGlobalDosagePresent) && dosage_reported;
    const uint32_t max_merged_allele_ct = MINV(max_allele_ct1 + max_allele_ct2, kPglMaxAlleleCt);
    const uint32_t max_allele_htable_size = GetHtableFastSize(max_merged_allele_ct);
    PgenVariant pgv1;
    PgenVariant pgv2;
    if (unlikely(BigstackAllocPgv(sample_ct, allele_idx_offsets != nullptr, dosage_needed? kfPgenGlobalDosagePresent : kfPgenGlobal0, &pgv1) ||
                 BigstackAllocPgv(sample_ct, allele_idx_offsets2 != nullptr, dosage_needed? kfPgenGlobalDosagePresent : kfPgenGlobal0, &pgv2))) {
      goto PgenDiff_ret_NOMEM;
    }
    // separated to avoid spurious maybe-uninitialized warnings
    const char** cur_allele2s;
    const char** merged_alleles;
    uint32_t* merged_alleles_htable;
    uint32_t* diff_sample_idxs;
    AlleleCode* remap1;
    AlleleCode* remap2;
    AlleleCode* wide_codes1;
    AlleleCode* wide_codes2;
    uintptr_t* remap_seen;
    if (unlikely(bigstack_alloc_kcp(max_allele_ct2, &cur_allele2s) ||
                 bigstack_alloc_kcp(max_merged_allele_ct, &merged_alleles) ||
                 bigstack_alloc_u32(max_allele_htable_size, &merged_alleles_htable) ||
                 bigstack_alloc_u32(sample_ct, &diff_sample_idxs) ||
                 bigstack_alloc_ac(max_allele_ct1, &remap1) ||
                 bigstack_alloc_ac(max_allele_ct2, &remap2) ||
                 bigstack_alloc_ac(sample_ct * 2, &wide_codes1) ||
                 bigstack_alloc_ac(sample_ct * 2, &wide_codes2) ||
                 bigstack_alloc_w(BitCtToWordCt(max_merged_allele_ct), &remap_seen))) {
      goto PgenDiff_ret_NOMEM;
    }
    uintptr_t* sex_male_collapsed = nullptr;
    Dosage* biallelic_dosage1 = nullptr;
    Dosage* biallelic_dosage2 = nullptr;
    if (dosage_needed) {
      // allocations are automatically rounded up to vector boundary, so
      // PopulateDenseDosage() is safe
      if (unlikely(bigstack_alloc_dosage(sample_ct, &biallelic_dosage1) ||
                   bigstack_alloc_dosage(sample_ct, &biallelic_dosage2))) {
        goto PgenDiff_ret_NOMEM;
      }
      if (sex_needed) {
        const uintptr_t sample_ctl = BitCtToWordCt(sample_ct);
        if (unlikely(bigstack_alloc_w(sample_ctl, &sex_male_collapsed))) {
          goto PgenDiff_ret_NOMEM;
        }
        CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_male_collapsed);
      }
    }
    pgv1.patch_01_ct = 0;
    pgv1.patch_10_ct = 0;
    pgv1.dosage_ct = 0;
    pgv2.patch_01_ct = 0;
    pgv2.patch_10_ct = 0;
    pgv2.dosage_ct = 0;
    PgenDiffGtEntry* gt_entries = nullptr;
    Dosage* ds_entries = nullptr;
    if (!dosage_needed) {
      if (unlikely(BIGSTACK_ALLOC_X(PgenDiffGtEntry, sample_ct, &gt_entries))) {
        goto PgenDiff_ret_NOMEM;
      }
    } else {
      // don't define a struct for this, since entry length depends on number
      // of alleles
      if (unlikely(bigstack_alloc_dosage(sample_ct * (2 * k1LU) * (max_merged_allele_ct - 1), &ds_entries))) {
        goto PgenDiff_ret_NOMEM;
      }
    }

    const uint32_t output_zst = (flags / kfPgenDiffZs) & 1;
    OutnameZstSet(".pdiff", output_zst, outname_end);
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    const uintptr_t overflow_buf_size = kCompressStreamBlock + max_chr_blen + 4 * kMaxIdSlen + 128 + max_allele_slen + max_merged_allele_ct * 16;
    const uint32_t compress_thread_ct = MINV((max_thread_ct + 1) / 2, 6);
    reterr = InitCstreamAlloc(outname, 0, output_zst, compress_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto PgenDiff_ret_1;
    }
    *cswritep++ = '#';
    const uint32_t chr_col = flags & kfPgenDiffColChrom;

    // null-terminated
    char* chr_buf;
    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto PgenDiff_ret_NOMEM;
    }
    if (chr_col) {
      cswritep = strcpya_k(cswritep, "CHROM\t");
    }
    const uint32_t pos_col = flags & kfPgenDiffColPos;
    if (pos_col) {
      cswritep = strcpya_k(cswritep, "POS\t");
    }
    const uint32_t varid_col = flags & kfPgenDiffColId;
    if (varid_col) {
      cswritep = strcpya_k(cswritep, "ID\t");
    }
    const uint32_t ref_col = flags & kfPgenDiffColRef;
    if (ref_col) {
      cswritep = strcpya_k(cswritep, "REF\t");
    }
    const uint32_t alt_col = flags & kfPgenDiffColAlt;
    if (alt_col) {
      cswritep = strcpya_k(cswritep, "ALT\t");
    }
    const uint32_t fid_col = FidColIsRequired(siip, flags / kfPgenDiffColMaybefid);
    if (fid_col) {
      cswritep = strcpya_k(cswritep, "FID\t");
    }
    cswritep = strcpya_k(cswritep, "IID");
    const uint32_t sid_col = SidColIsRequired(siip->sids, flags / kfPgenDiffColMaybesid);
    if (sid_col) {
      cswritep = strcpya_k(cswritep, "\tSID");
    }
    const uint32_t geno_col = flags & kfPgenDiffColGeno;
    if (geno_col) {
      if (!dosage_reported) {
        cswritep = strcpya_k(cswritep, "\tGT1\tGT2");
      } else {
        cswritep = strcpya_k(cswritep, "\tDS1\tDS2");
      }
    }
    AppendBinaryEoln(&cswritep);

    reterr = InitTextStream(pdip->pvar_fname, MAXV(max_observed_pvar_line_blen, kTextStreamBlenFast), ClipU32(max_thread_ct - 1, 1, 3), &pvar_txs);
    if (unlikely(reterr)) {
      goto PgenDiff_ret_PVAR_REWIND_FAIL;
    }
    // Skip to header.
    char* line_start;
    do {
      line_start = TextGet(&pvar_txs);
      if (unlikely(!line_start)) {
        goto PgenDiff_ret_PVAR_REWIND_FAIL;
      }
      ++pvar_line_idx;
    } while ((*line_start == '#') && (!tokequal_k(line_start, "#CHROM")));
    uint32_t col_skips[4];
    uint32_t col_types[4];
    if (*line_start == '#') {
      // [-1] = #CHROM (must be first column)
      // [0] = POS
      // [1] = ID
      // [2] = REF
      // [3] = ALT
      // Can ignore other columns here.
      char* token_end = &(line_start[6]);
      uint32_t found_header_bitset = 0;
      uint32_t relevant_postchr_col_uidx = 0;
      for (uint32_t col_idx = 1; ; ++col_idx) {
        char* token_start = FirstNonTspace(token_end);
        if (IsEolnKns(*token_start)) {
          break;
        }
        token_end = CurTokenEnd(token_start);
        const uint32_t token_slen = token_end - token_start;
        uint32_t cur_col_type;
        if (strequal_k(token_start, "POS", token_slen)) {
          cur_col_type = 0;
        } else if (strequal_k(token_start, "ID", token_slen)) {
          cur_col_type = 1;
        } else if (strequal_k(token_start, "REF", token_slen)) {
          cur_col_type = 2;
        } else if (strequal_k(token_start, "ALT", token_slen)) {
          cur_col_type = 3;
        } else {
          continue;
        }
        const uint32_t cur_col_type_shifted = 1 << cur_col_type;
        if (unlikely(found_header_bitset & cur_col_type_shifted)) {
          // known token, so no overflow danger
          char* write_iter = strcpya_k(g_logbuf, "Error: Duplicate column header '");
          write_iter = memcpya(write_iter, token_start, token_slen);
          write_iter = strcpya_k(write_iter, "' on line ");
          write_iter = wtoa(pvar_line_idx, write_iter);
          write_iter = strcpya_k(write_iter, " of --pgen-diff file.\n");
          *write_iter = '\0';
          goto PgenDiff_ret_MALFORMED_INPUT_WW;
        }
        found_header_bitset |= cur_col_type_shifted;
        col_skips[relevant_postchr_col_uidx] = col_idx;
        col_types[relevant_postchr_col_uidx++] = cur_col_type;
      }
      if (unlikely(found_header_bitset != 0xf)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Missing column header(s) on line %" PRIuPTR " of --pgen-diff file. (POS, ID, REF, and ALT are required.)\n", pvar_line_idx);
        goto PgenDiff_ret_MALFORMED_INPUT_WW;
      }
      for (uint32_t rpc_col_idx = 3; rpc_col_idx; --rpc_col_idx) {
        col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
      }
      line_start = AdvPastDelim(token_end, '\n');
      ++pvar_line_idx;
    } else {
      char* fifth_token_start = NextTokenMult(line_start, 4);
      if (unlikely(!fifth_token_start)) {
        goto PgenDiff_ret_PVAR_MISSING_TOKENS;
      }
      col_skips[0] = 1;
      col_skips[2] = 1;
      col_skips[3] = 1;
      col_types[0] = 1;
      col_types[1] = 0;
      col_types[2] = 3;
      col_types[3] = 2;
      if (!NextToken(fifth_token_start)) {
        // #CHROM ID POS ALT REF
        col_skips[1] = 1;
      } else {
        // #CHROM ID CM POS ALT REF
        col_skips[1] = 2;
      }
    }

    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* already_seen;
    if (unlikely(bigstack_calloc_w(raw_variant_ctl, &already_seen))) {
      goto PgenDiff_ret_NOMEM;
    }
    const DoubleAlleleCode biallelic_dac[4] = {0, 1 << (8 * sizeof(AlleleCode)), 1 + (1 << (8 * sizeof(AlleleCode))), kMissingDoubleAlleleCode};
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    const uint32_t include_missing = (flags / kfPgenDiffIncludeMissing) & 1;
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    const char input_missing_geno_char = *g_input_missing_geno_ptr;
    const uintptr_t* nonref_flags1 = pgfip->nonref_flags;
    const uint32_t all_nonref1 = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags1);
    uint32_t chr_idx = UINT32_MAX;
    uint32_t chr_slen = 0;
    uint32_t cur_bp = 0; // just for .pvar-sorted sanity check
    uint32_t cur_included_bp = 0;
    uint32_t variant_uidx_start = 0;
    uint32_t variant_uidx_end = 0;
    uint32_t same_bp_variant_ct = 0;
    uint32_t chrom_end_variant_uidx = 0;
    uint32_t is_autosomal_diploid = 0;
    uint32_t is_x = 0;
    uint32_t is_y = 0;
    uint32_t dosage_cur_tol = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_uidx2 = raw_variant_ct2 / 100;
    uint64_t grand_diff_ct = 0;
    uint32_t cur_allele_ct1 = 2;
    uint32_t cur_allele_ct2 = 2;
    fputs("--pgen-diff: 0%", stdout);
    fflush(stdout);
    for (uint32_t variant_uidx2 = 0; variant_uidx2 != raw_variant_ct2; ++variant_uidx2, ++pvar_line_idx) {
      if (unlikely(!TextGetUnsafe2(&pvar_txs, &line_start))) {
        goto PgenDiff_ret_PVAR_REWIND_FAIL_N;
      }
      char* token_ptrs[4];
      uint32_t token_slens[4];
      char* last_token_end = TokenLex(line_start, col_types, col_skips, 4, token_ptrs, token_slens);
      if (unlikely(!last_token_end)) {
        goto PgenDiff_ret_PVAR_MISSING_TOKENS_N;
      }
      if (variant_uidx2 >= next_print_variant_uidx2) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_uidx2 * 100LLU) / raw_variant_ct2;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_uidx2 = (pct * S_CAST(uint64_t, raw_variant_ct2)) / 100;
      }
      {
        char* chrom_start = line_start;
        line_start = AdvPastDelim(line_start, '\n');

        char* chrom_end = CurTokenEnd(chrom_start);
        *chrom_end = '\0';
        const uint32_t new_chr_idx = GetChrCode(chrom_start, cip, chrom_end - chrom_start);
        if (IsI32Neg(new_chr_idx)) {
          continue;
        }
        if (new_chr_idx != chr_idx) {
          chr_idx = new_chr_idx;
          if (IsI32Neg(chr_idx) || !IsSet(cip->chr_mask, chr_idx)) {
            // Chromosome isn't loaded at all; skip it.
            chrom_end_variant_uidx = 0;
            continue;
          }
          const uint32_t is_haploid = IsSet(cip->haploid_mask, chr_idx);
          is_autosomal_diploid = 1 - is_haploid;
          is_x = (chr_idx == x_code);
          is_y = (chr_idx == y_code);
          dosage_cur_tol = dosage_sex_tols[is_haploid];
          cur_bp = 0;
          const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
          const uint32_t chrom_start_variant_uidx = cip->chr_fo_vidx_start[chr_fo_idx];
          chrom_end_variant_uidx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          variant_uidx_start = AdvBoundedTo1Bit(variant_include, chrom_start_variant_uidx, chrom_end_variant_uidx);
          if (variant_uidx_start == chrom_end_variant_uidx) {
            // All variants on this chromosome have been filtered out; skip it.
            chrom_end_variant_uidx = 0;
            continue;
          }
          char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
          *chr_name_end = '\0';
          chr_slen = chr_name_end - chr_buf;
          cur_included_bp = variant_bps[variant_uidx_start];
          const uint32_t search_start = variant_uidx_start + 1;
          variant_uidx_end = search_start + ExpsearchU32(&(variant_bps[search_start]), chrom_end_variant_uidx - search_start, cur_included_bp + 1);
          variant_uidx_end = 1 + FindLast1BitBefore(variant_include, variant_uidx_end);
          same_bp_variant_ct = PopcountBitRange(variant_include, variant_uidx_start, variant_uidx_end);
        } else if (!chrom_end_variant_uidx) {
          continue;
        }
      }

      {
        int32_t new_bp;
        if (unlikely(ScanIntAbsDefcap(token_ptrs[0], &new_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of --pgen-diff .pvar file.\n", pvar_line_idx);
          goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
        }
        if (new_bp < 0) {
          continue;
        }
        const uint32_t new_bp_u32 = new_bp;
        if (new_bp_u32 < cur_bp) {
          // Could also verify that .pvar has no split chromosomes.
          logputs("\n");
          logerrputs("Error: --pgen-diff .pvar file is unsorted.\n");
          goto PgenDiff_ret_MALFORMED_INPUT;
        }
        cur_bp = new_bp_u32;
        if (new_bp_u32 < cur_included_bp) {
          continue;
        }
        if (new_bp_u32 > cur_included_bp) {
          variant_uidx_start = variant_uidx_end + ExpsearchU32(&(variant_bps[variant_uidx_end]), chrom_end_variant_uidx - variant_uidx_end, new_bp_u32);
          variant_uidx_start = AdvBoundedTo1Bit(variant_include, variant_uidx_start, chrom_end_variant_uidx);
          if (variant_uidx_start == chrom_end_variant_uidx) {
            variant_uidx_end = variant_uidx_start;
            cur_included_bp = UINT32_MAX;
            continue;
          }
          cur_included_bp = variant_bps[variant_uidx_start];
          const uint32_t search_start = variant_uidx_start + 1;
          variant_uidx_end = search_start + ExpsearchU32(&(variant_bps[search_start]), chrom_end_variant_uidx - search_start, cur_included_bp + 1);
          variant_uidx_end = 1 + FindLast1BitBefore(variant_include, variant_uidx_end);
          same_bp_variant_ct = PopcountBitRange(variant_include, variant_uidx_start, variant_uidx_end);
          if (new_bp_u32 < cur_included_bp) {
            continue;
          }
        }
      }

      uint32_t variant_uidx = UINT32_MAX;
      {
        // Note that this is inefficient if we're dealing with a large pile of
        // unmapped variants with CHROM and POS set to 0.  We can construct a
        // variant ID hash table when same_bp_variant_ct is large if it's
        // important to speed up that use case.
        char* cur_variant_id = token_ptrs[1];
        const uint32_t cur_variant_id_slen = token_slens[1];
        if (cur_variant_id_slen > kMaxIdSlen) {
          snprintf(g_logbuf, kLogbufSize, "Error: Variant ID too long on line %" PRIuPTR " of --pgen-diff .pvar file (max " MAX_ID_SLEN_STR " characters).\n", pvar_line_idx);
          goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
        }
        cur_variant_id[cur_variant_id_slen] = '\0';
        const uint32_t cur_variant_id_blen = cur_variant_id_slen + 1;
        uintptr_t variant_uidx_base;
        uintptr_t cur_bits;
        BitIter1Start(variant_include, variant_uidx_start, &variant_uidx_base, &cur_bits);
        for (uint32_t same_bp_variant_idx = 0; same_bp_variant_idx != same_bp_variant_ct; ++same_bp_variant_idx) {
          const uint32_t scan_variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          if (!memequal(cur_variant_id, variant_ids[scan_variant_uidx], cur_variant_id_blen)) {
            continue;
          }
          if (unlikely(variant_uidx != UINT32_MAX)) {
            snprintf(g_logbuf, kLogbufSize, "Error: --pgen-diff: Multiple instances of variant ID '%s' at position %s:%u in main fileset.\n", cur_variant_id, chr_buf, cur_included_bp);
            goto PgenDiff_ret_INCONSISTENT_INPUT_WW_N;
          }
          variant_uidx = scan_variant_uidx;
        }
        if (variant_uidx == UINT32_MAX) {
          continue;
        }
        if (unlikely(IsSet(already_seen, variant_uidx))) {
          snprintf(g_logbuf, kLogbufSize, "Error: --pgen-diff: Multiple instances of variant ID '%s' at position %s:%u in .pvar file.\n", cur_variant_id, chr_buf, cur_included_bp);
          goto PgenDiff_ret_INCONSISTENT_INPUT_WW_N;
        }
        SetBit(variant_uidx, already_seen);
      }

      uintptr_t allele_idx_offset_base1 = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base1 = allele_idx_offsets[variant_uidx];
        cur_allele_ct1 = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base1;
      }
      const char* const* cur_allele1s = &(allele_storage[allele_idx_offset_base1]);
      uintptr_t allele_idx_offset_base2 = variant_uidx2 * 2;
      if (allele_idx_offsets2) {
        allele_idx_offset_base2 = allele_idx_offsets2[variant_uidx2];
        cur_allele_ct2 = allele_idx_offsets2[variant_uidx2 + 1] - allele_idx_offset_base2;
      }

      if (!dosage_needed) {
        if (cur_allele_ct1 == 2) {
          reterr = PgrGet(sample_include, pssi1, sample_ct, variant_uidx, simple_pgrp, pgv1.genovec);
        } else {
          reterr = PgrGetM(sample_include, pssi1, sample_ct, variant_uidx, simple_pgrp, &pgv1);
        }
        if (unlikely(reterr)) {
          goto PgenDiff_ret_PGR_FAIL;
        }
        if (cur_allele_ct2 == 2) {
          reterr = PgrGet(sample_include2, pssi2, sample_ct, variant_uidx2, &simple_pgr2, pgv2.genovec);
        } else {
          reterr = PgrGetM(sample_include2, pssi2, sample_ct, variant_uidx2, &simple_pgr2, &pgv2);
        }
      } else {
        if (cur_allele_ct1 == 2) {
          reterr = PgrGetD(sample_include, pssi1, sample_ct, variant_uidx, simple_pgrp, pgv1.genovec, pgv1.dosage_present, pgv1.dosage_main, &pgv1.dosage_ct);
        } else {
          reterr = PgrGetMD(sample_include, pssi1, sample_ct, variant_uidx, simple_pgrp, &pgv1);
        }
        if (unlikely(reterr)) {
          goto PgenDiff_ret_PGR_FAIL;
        }
        if (cur_allele_ct2 == 2) {
          reterr = PgrGetD(sample_include2, pssi2, sample_ct, variant_uidx2, &simple_pgr2, pgv2.genovec, pgv2.dosage_present, pgv2.dosage_main, &pgv2.dosage_ct);
        } else {
          reterr = PgrGetMD(sample_include2, pssi2, sample_ct, variant_uidx2, &simple_pgr2, &pgv2);
        }
      }
      if (unlikely(reterr)) {
        goto PgenDiff_ret_PGR_FAIL;
      }

      uintptr_t* genovec1 = pgv1.genovec;
      uintptr_t* genovec2 = pgv2.genovec;
      uint32_t merged_allele_ct = 0;
      {
        char* ref_allele2 = token_ptrs[2];
        const uint32_t ref_allele2_slen = token_slens[2];
        ref_allele2[ref_allele2_slen] = '\0';
        cur_allele2s[0] = ref_allele2;

        const uint32_t last_allele2_idx = cur_allele_ct2 - 1;
        char* alt_iter = token_ptrs[3];
        alt_iter[token_slens[3]] = '\0';
        for (uint32_t allele2_idx = 1; allele2_idx != last_allele2_idx; ++allele2_idx) {
          char* alt_end = AdvToDelim(alt_iter, ',');
          *alt_end = '\0';
          cur_allele2s[allele2_idx] = alt_iter;
          alt_iter = &(alt_end[1]);
        }
        cur_allele2s[last_allele2_idx] = alt_iter;

        // Missing allele codes are only permitted in biallelic variants, and
        // they must not appear in the actual genotype calls.
        uint32_t missing1_idx = UINT32_MAX;
        for (uint32_t allele1_idx = 0; allele1_idx != cur_allele_ct1; ++allele1_idx) {
          const char* cur_allele1 = cur_allele1s[allele1_idx];
          if (memequal(cur_allele1, ".", 2)) {
            if (unlikely(missing1_idx != UINT32_MAX)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Multiple missing alleles for variant '%s' at position %s:%u.\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
              goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
            }
            missing1_idx = allele1_idx;
          }
        }
        ZeroTrailingNyps(sample_ct, genovec1);
        ZeroTrailingNyps(sample_ct, genovec2);
        if (missing1_idx != UINT32_MAX) {
          if (unlikely(cur_allele_ct1 > 2)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Missing allele in multiallelic variant '%s' at position %s:%u.\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
            goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
          }
          STD_ARRAY_DECL(uint32_t, 4, genocounts);
          GenoarrCountFreqsUnsafe(genovec1, sample_ct, genocounts);
          if (unlikely(genocounts[1] || genocounts[2 * missing1_idx] || pgv1.patch_01_ct || pgv1.patch_10_ct || pgv1.dosage_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Missing allele for variant '%s' at position %s:%u is present in the .pgen.\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
            goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
          }
          remap1[missing1_idx] = kMissingAlleleCode;
        }
        uint32_t missing2_idx = UINT32_MAX;
        for (uint32_t allele2_idx = 0; allele2_idx != cur_allele_ct2; ++allele2_idx) {
          const char* cur_allele2 = cur_allele2s[allele2_idx];
          if (((cur_allele2[0] == '.') || (cur_allele2[0] == input_missing_geno_char)) && (cur_allele2[1] == '\0')) {
            if (unlikely(missing2_idx != UINT32_MAX)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Multiple missing alleles on line %" PRIuPTR " of --pgen-diff .pvar file.\n", pvar_line_idx);
              goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
            }
            missing2_idx = allele2_idx;
          }
        }
        if (missing2_idx != UINT32_MAX) {
          if (unlikely(cur_allele_ct2 > 2)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Missing allele in multiallelic variant on line %" PRIuPTR " of --pgen-diff .pvar file.\n", pvar_line_idx);
            goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
          }
          STD_ARRAY_DECL(uint32_t, 4, genocounts);
          GenoarrCountFreqsUnsafe(genovec2, sample_ct, genocounts);
          if (unlikely(genocounts[1] || genocounts[2 * missing2_idx] || pgv2.patch_01_ct || pgv2.patch_10_ct || pgv2.dosage_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Missing allele on line %" PRIuPTR " of --pgen-diff .pvar file is present in the .pgen.\n", pvar_line_idx);
            goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
          }
          remap2[missing2_idx] = kMissingAlleleCode;
        }
        const uint32_t provisional_ref1 = all_nonref1 || (nonref_flags1 && IsSet(nonref_flags1, variant_uidx));
        if ((!provisional_ref1) && (missing1_idx == 0)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Missing REF allele for variant '%s' at position %s:%u is not flagged as provisional.\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
          goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
        }
        const uint32_t provisional_ref2 = all_nonref2 || (nonref_flags2 && IsSet(nonref_flags2, variant_uidx2));
        if ((!provisional_ref2) && (missing2_idx == 0)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Missing REF allele on line %" PRIuPTR " of --pgen-diff .pvar file is not flagged as provisional.\n", pvar_line_idx);
          goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
        }
        uint32_t allele1_alt_start = (missing1_idx == 0);
        uint32_t allele2_alt_start = (missing2_idx == 0);
        if ((!provisional_ref1) || ((missing1_idx != 0) && provisional_ref2)) {
          if (unlikely((!provisional_ref2) && (!strequal_overread(cur_allele1s[0], cur_allele2s[0])))) {
            snprintf(g_logbuf, kLogbufSize, "Error: REF allele on line %" PRIuPTR " of --pgen-diff .pvar file conflicts with loaded REF allele, and neither are flagged as provisional.\n", pvar_line_idx);
            goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
          }
          merged_alleles[0] = cur_allele1s[0];
          remap1[0] = 0;
          allele1_alt_start = 1 + (missing1_idx == 1);
        } else if (missing2_idx != 0) {
          merged_alleles[0] = cur_allele2s[0];
          remap2[0] = 0;
          allele2_alt_start = 1 + (missing2_idx == 1);
        } else {
          // Both REFs missing; may as well treat main dataset ALT1 as REF.
          merged_alleles[0] = cur_allele1s[1];
          remap1[1] = 0;
          allele1_alt_start = 2;
        }
        const uint32_t cur_max_merged_alleles = 1 + cur_allele_ct1 + cur_allele_ct2 - allele1_alt_start - allele2_alt_start;
        const uint32_t cur_htable_size = (cur_max_merged_alleles * 9) / 2;
        SetAllU32Arr(cur_htable_size, merged_alleles_htable);
        // See e.g PopulateStrboxHtable().
        uint32_t hashval = Hashceil(merged_alleles[0], strlen(merged_alleles[0]), cur_htable_size);
        merged_alleles_htable[hashval] = 0;
        ++merged_allele_ct;

        for (uint32_t allele1_idx = allele1_alt_start; allele1_idx != cur_allele_ct1; ++allele1_idx) {
          if (allele1_idx == missing1_idx) {
            continue;
          }
          const char* cur_allele1 = cur_allele1s[allele1_idx];
          const uint32_t cur_allele1_slen = strlen(cur_allele1);
          for (hashval = Hashceil(cur_allele1, cur_allele1_slen, cur_htable_size); ; ) {
            const uint32_t cur_htable_entry = merged_alleles_htable[hashval];
            if (cur_htable_entry == UINT32_MAX) {
              if (unlikely(merged_allele_ct == kPglMaxAlleleCt)) {
                logerrprintfww("Error: Too many alleles across --pgen-diff and main filesets for variant '%s' at position %s:%u. (This " PROG_NAME_STR " build is limited to " PGL_MAX_ALLELE_CT_STR ".)\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
                reterr = kPglRetNotYetSupported;
                goto PgenDiff_ret_1;
              }
              remap1[allele1_idx] = merged_allele_ct;
              merged_alleles_htable[hashval] = merged_allele_ct;
              merged_alleles[merged_allele_ct] = cur_allele1;
              ++merged_allele_ct;
              break;
            }
            if (memequal(cur_allele1, merged_alleles[cur_htable_entry], cur_allele1_slen + 1)) {
              remap1[allele1_idx] = cur_htable_entry;
              break;
            }
            if (++hashval == cur_htable_size) {
              hashval = 0;
            }
          }
        }
        for (uint32_t allele2_idx = allele2_alt_start; allele2_idx != cur_allele_ct2; ++allele2_idx) {
          if (allele2_idx == missing2_idx) {
            continue;
          }
          const char* cur_allele2 = cur_allele2s[allele2_idx];
          const uint32_t cur_allele2_slen = strlen(cur_allele2);
          for (hashval = Hashceil(cur_allele2, cur_allele2_slen, cur_htable_size); ; ) {
            const uint32_t cur_htable_entry = merged_alleles_htable[hashval];
            if (cur_htable_entry == UINT32_MAX) {
              if (unlikely(merged_allele_ct == kPglMaxAlleleCt)) {
                logerrprintfww("Error: Too many alleles across --pgen-diff and main filesets for variant '%s' at position %s:%u. (This " PROG_NAME_STR " build is limited to " PGL_MAX_ALLELE_CT_STR ".)\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
                reterr = kPglRetNotYetSupported;
                goto PgenDiff_ret_1;
              }
              remap2[allele2_idx] = merged_allele_ct;
              merged_alleles_htable[hashval] = merged_allele_ct;
              merged_alleles[merged_allele_ct] = cur_allele2;
              ++merged_allele_ct;
              break;
            }
            if (memequal(cur_allele2, merged_alleles[cur_htable_entry], cur_allele2_slen + 1)) {
              remap2[allele2_idx] = cur_htable_entry;
              break;
            }
            if (++hashval == cur_htable_size) {
              hashval = 0;
            }
          }
        }
        const uint32_t merged_allele_ctl = BitCtToWordCt(merged_allele_ct);
        if (unlikely(DuplicateAllelePresent(remap1, cur_allele_ct1, merged_allele_ctl, remap_seen))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate allele code in variant '%s' at position %s:%u.\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
          goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
        }
        if (unlikely(DuplicateAllelePresent(remap2, cur_allele_ct2, merged_allele_ctl, remap_seen))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate allele code on line %" PRIuPTR " of --pgen-diff .pvar file.\n", pvar_line_idx);
          goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
        }
        if (merged_allele_ct == 1) {
          if (!include_missing) {
            // Only possible difference is missing vs. not-missing.
            continue;
          }
          merged_alleles[1] = cur_allele1s[missing1_idx];
          merged_allele_ct = 2;
        }
      }
      uint32_t diff_ct = 0;
      if (merged_allele_ct == 2) {
        // Note that these two conditions aren't necessarily synonymous, due to
        // missing alleles.
        if ((remap1[0] == 1) || (remap1[1] == 0)) {
          GenovecInvertUnsafe(sample_ct, genovec1);
          ZeroTrailingNyps(sample_ct, genovec1);
          if (dosage_needed) {
            BiallelicDosage16Invert(pgv1.dosage_ct, pgv1.dosage_main);
          }
        }
        if ((remap2[0] == 1) || (remap2[1] == 0)) {
          GenovecInvertUnsafe(sample_ct, genovec2);
          ZeroTrailingNyps(sample_ct, genovec2);
          if (dosage_needed) {
            BiallelicDosage16Invert(pgv2.dosage_ct, pgv2.dosage_main);
          }
        }
        if (!dosage_needed) {
          if (!sample1_idx_to_2) {
            // Optimize common case where no reshuffling is needed.
            for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
              uintptr_t geno_word1 = genovec1[widx];
              uintptr_t geno_word2 = genovec2[widx];
              uintptr_t diff_word = geno_word1 ^ geno_word2;
              if (!diff_word) {
                continue;
              }
              if (!include_missing) {
                const uintptr_t missing_word1 = geno_word1 & (geno_word1 >> 1) & kMask5555;
                const uintptr_t missing_word2 = geno_word2 & (geno_word2 >> 1) & kMask5555;
                diff_word &= ~((missing_word1 | missing_word2) * 3);
                if (!diff_word) {
                  continue;
                }
              }
              const uint32_t sample_idx_base = widx * kBitsPerWordD2;
              do {
                const uint32_t rshift = ctzw(diff_word) & (~1);
                const uintptr_t geno1 = (geno_word1 >> rshift) & 3;
                const uintptr_t geno2 = (geno_word2 >> rshift) & 3;
                diff_sample_idxs[diff_ct] = sample_idx_base + (rshift / 2);
                PgenDiffGtEntry* gt_entryp = &(gt_entries[diff_ct]);
                gt_entryp->dac1 = biallelic_dac[geno1];
                gt_entryp->dac2 = biallelic_dac[geno2];
                ++diff_ct;
                diff_word &= ~((3 * k1LU) << rshift);
              } while (diff_word);
            }
          } else {
            // biallelic, !dosage_needed, must reshuffle
            // possible todo: parallelize this
            const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
            uint32_t loop_len = kBitsPerWordD2;
            for (uint32_t widx = 0; ; ++widx) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                loop_len = ModNz(sample_ct, kBitsPerWordD2);
              }
              const uint32_t sample_idx_base = widx * kBitsPerWordD2;
              const uint32_t* cur_sample1_idx_to_2 = &(sample1_idx_to_2[sample_idx_base]);
              uintptr_t geno_word1 = genovec1[widx];
              if (include_missing) {
                for (uint32_t uii = 0; uii != loop_len; ++uii) {
                  const uint32_t sample_idx2 = cur_sample1_idx_to_2[uii];
                  const uintptr_t geno1 = geno_word1 & 3;
                  const uintptr_t geno2 = GetNyparrEntry(genovec2, sample_idx2);
                  geno_word1 >>= 2;
                  if (geno1 == geno2) {
                    continue;
                  }
                  diff_sample_idxs[diff_ct] = sample_idx_base + uii;
                  PgenDiffGtEntry* gt_entryp = &(gt_entries[diff_ct]);
                  gt_entryp->dac1 = biallelic_dac[geno1];
                  gt_entryp->dac2 = biallelic_dac[geno2];
                  ++diff_ct;
                }
              } else {
                for (uint32_t uii = 0; uii != loop_len; ++uii) {
                  const uint32_t sample_idx2 = cur_sample1_idx_to_2[uii];
                  const uintptr_t geno1 = geno_word1 & 3;
                  const uintptr_t geno2 = GetNyparrEntry(genovec2, sample_idx2);
                  geno_word1 >>= 2;
                  if ((geno1 == geno2) || (geno1 == 3) || (geno2 == 3)) {
                    continue;
                  }
                  diff_sample_idxs[diff_ct] = sample_idx_base + uii;
                  PgenDiffGtEntry* gt_entryp = &(gt_entries[diff_ct]);
                  gt_entryp->dac1 = biallelic_dac[geno1];
                  gt_entryp->dac2 = biallelic_dac[geno2];
                  ++diff_ct;
                }
              }
            }
          }
        } else {
          // biallelic, dosage_needed
          PopulateDenseDosage(genovec1, pgv1.dosage_present, pgv1.dosage_main, sample_ct, pgv1.dosage_ct, biallelic_dosage1);
          PopulateDenseDosage(genovec2, pgv2.dosage_present, pgv2.dosage_main, sample_ct, pgv2.dosage_ct, biallelic_dosage2);
          if (!sample1_idx_to_2) {
            const uintptr_t dosage_blen = sample_ct * sizeof(Dosage);
            uintptr_t sample_idx = 0;
            if (!is_x) {
              // don't need to special-case chrY here since we just skip
              // nonmales in the reporting step
              while (1) {
                const uintptr_t diff_byte_offset = FirstUnequalFrom(biallelic_dosage1, biallelic_dosage2, sample_idx * sizeof(Dosage), dosage_blen);
                if (diff_byte_offset == dosage_blen) {
                  break;
                }
                sample_idx = diff_byte_offset / sizeof(Dosage);
                const uint32_t d1 = biallelic_dosage1[sample_idx];
                const uint32_t d2 = biallelic_dosage2[sample_idx];
                if (((d1 != kDosageMissing) && (d2 != kDosageMissing)) || include_missing) {
                  // Since this doesn't separate out the kDosageMissing cases,
                  // it needs to be updated if Dosage is widened (since
                  // kDosageMissing would then be adjacent to 0 in uint32_t
                  // space).
                  if (abs_i32(d1 - d2) > dosage_cur_tol) {
                    diff_sample_idxs[diff_ct] = sample_idx;
                    ds_entries[diff_ct * 2] = d1;
                    ds_entries[diff_ct * 2 + 1] = d2;
                    ++diff_ct;
                  }
                }
                ++sample_idx;
              }
            } else {
              // is_x
              while (1) {
                const uintptr_t diff_byte_offset = FirstUnequalFrom(biallelic_dosage1, biallelic_dosage2, sample_idx * sizeof(Dosage), dosage_blen);
                if (diff_byte_offset == dosage_blen) {
                  break;
                }
                sample_idx = diff_byte_offset / sizeof(Dosage);
                const uint32_t d1 = biallelic_dosage1[sample_idx];
                const uint32_t d2 = biallelic_dosage2[sample_idx];
                if (((d1 != kDosageMissing) && (d2 != kDosageMissing)) || include_missing) {
                  const uint32_t is_male = IsSet(sex_male_collapsed, sample_idx);
                  if (abs_i32(d1 - d2) > dosage_sex_tols[is_male]) {
                    diff_sample_idxs[diff_ct] = sample_idx;
                    ds_entries[diff_ct * 2] = d1;
                    ds_entries[diff_ct * 2 + 1] = d2;
                    ++diff_ct;
                  }
                }
                ++sample_idx;
              }
            }
          } else {
            if (!is_x) {
              for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
                const uint32_t d1 = biallelic_dosage1[sample_idx];
                const uint32_t sample_idx2 = sample1_idx_to_2[sample_idx];
                const uint32_t d2 = biallelic_dosage2[sample_idx2];
                if (((d1 != kDosageMissing) && (d2 != kDosageMissing)) || include_missing) {
                  if (abs_i32(d1 - d2) > dosage_cur_tol) {
                    diff_sample_idxs[diff_ct] = sample_idx;
                    ds_entries[diff_ct * 2] = d1;
                    ds_entries[diff_ct * 2 + 1] = d2;
                    ++diff_ct;
                  }
                }
              }
            } else {
              for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
                const uint32_t d1 = biallelic_dosage1[sample_idx];
                const uint32_t sample_idx2 = sample1_idx_to_2[sample_idx];
                const uint32_t d2 = biallelic_dosage2[sample_idx2];
                if (((d1 != kDosageMissing) && (d2 != kDosageMissing)) || include_missing) {
                  const uint32_t is_male = IsSet(sex_male_collapsed, sample_idx);
                  if (abs_i32(d1 - d2) > dosage_sex_tols[is_male]) {
                    diff_sample_idxs[diff_ct] = sample_idx;
                    ds_entries[diff_ct * 2] = d1;
                    ds_entries[diff_ct * 2 + 1] = d2;
                    ++diff_ct;
                  }
                }
              }
            }
          }
        }
      } else {
        if (!dosage_needed) {
          PglMultiallelicSparseToDense(genovec1, pgv1.patch_01_set, pgv1.patch_01_vals, pgv1.patch_10_set, pgv1.patch_10_vals, remap1, sample_ct, pgv1.patch_01_ct, pgv1.patch_10_ct, nullptr, wide_codes1);
          PglMultiallelicSparseToDense(genovec2, pgv2.patch_01_set, pgv2.patch_01_vals, pgv2.patch_10_set, pgv2.patch_10_vals, remap2, sample_ct, pgv2.patch_01_ct, pgv2.patch_10_ct, nullptr, wide_codes2);
          const DoubleAlleleCode* wc1_alias = R_CAST(DoubleAlleleCode*, wide_codes1);
          const DoubleAlleleCode* wc2_alias = R_CAST(DoubleAlleleCode*, wide_codes2);
          if (!sample1_idx_to_2) {
            const uintptr_t wide_codes_blen = sample_ct * 2 * sizeof(AlleleCode);
            uintptr_t sample_idx = 0;
            while (1) {
              const uintptr_t diff_byte_offset = FirstUnequalFrom(wc1_alias, wc2_alias, sample_idx * sizeof(DoubleAlleleCode), wide_codes_blen);
              if (diff_byte_offset == wide_codes_blen) {
                break;
              }
              sample_idx = diff_byte_offset / (2 * sizeof(AlleleCode));
              const DoubleAlleleCode dac1 = wc1_alias[sample_idx];
              const DoubleAlleleCode dac2 = wc2_alias[sample_idx];
              if (((dac1 != kMissingDoubleAlleleCode) && (dac2 != kMissingDoubleAlleleCode)) || include_missing) {
                diff_sample_idxs[diff_ct] = sample_idx;
                PgenDiffGtEntry* gt_entryp = &(gt_entries[diff_ct]);
                gt_entryp->dac1 = dac1;
                gt_entryp->dac2 = dac2;
                ++diff_ct;
              }
              ++sample_idx;
            }
          } else {
            for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
              const DoubleAlleleCode dac1 = wc1_alias[sample_idx];
              const uint32_t sample_idx2 = sample1_idx_to_2[sample_idx];
              const DoubleAlleleCode dac2 = wc2_alias[sample_idx2];
              if ((dac1 == dac2) || (((dac1 == kMissingAlleleCode) || (dac2 == kMissingAlleleCode)) && (!include_missing))) {
                continue;
              }
              diff_sample_idxs[diff_ct] = sample_idx;
              PgenDiffGtEntry* gt_entryp = &(gt_entries[diff_ct]);
              gt_entryp->dac1 = dac1;
              gt_entryp->dac2 = dac2;
              ++diff_ct;
            }
          }
        } else {
          // multiallelic, dosage_needed
          logerrputs("Error: --pgen-diff multiallelic-variant dosage support is under development.\n");
          reterr = kPglRetNotYetSupported;
          goto PgenDiff_ret_1;
        }
      }
      if (!diff_ct) {
        continue;
      }
      grand_diff_ct += diff_ct;
      const uint32_t merged_allele_ct_m1 = merged_allele_ct - 1;
      uint32_t cur_autosomal_diploid = is_autosomal_diploid;
      for (uint32_t diff_idx = 0; diff_idx != diff_ct; ++diff_idx) {
        const uint32_t sample_uidx = sample_idx_to_uidx[diff_sample_idxs[diff_idx]];
        if (is_x) {
          cur_autosomal_diploid = !IsSet(sex_male, sample_uidx);
        } else if (is_y) {
          if (!IsSet(sex_male, sample_uidx)) {
            --grand_diff_ct;
            continue;
          }
        }
        if (chr_col) {
          cswritep = memcpyax(cswritep, chr_buf, chr_slen, '\t');
        }
        if (pos_col) {
          cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
        }
        if (varid_col) {
          cswritep = strcpyax(cswritep, variant_ids[variant_uidx], '\t');
        }
        if (ref_col) {
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto PgenDiff_ret_WRITE_FAIL;
          }
          cswritep = strcpyax(cswritep, merged_alleles[0], '\t');
        }
        if (alt_col) {
          for (uint32_t allele_idx = 1; allele_idx != merged_allele_ct; ++allele_idx) {
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto PgenDiff_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, merged_alleles[allele_idx], ',');
          }
          cswritep[-1] = '\t';
        }
        const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
        if (!fid_col) {
          cur_sample_id = AdvPastDelim(cur_sample_id, '\t');
        }
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_sample_id);
        if (sid_col) {
          *cswritep++ = '\t';
          if (sids) {
            cswritep = strcpya(cswritep, &(sids[sample_uidx * max_sid_blen]));
          } else {
            *cswritep++ = '0';
          }
        }
        if (geno_col) {
          if (!dosage_needed) {
            const PgenDiffGtEntry entry = gt_entries[diff_idx];
            if (!dosage_reported) {
              *cswritep++ = '\t';
              if (entry.dac1 == kMissingDoubleAlleleCode) {
                if (cur_autosomal_diploid) {
                  cswritep = strcpya_k(cswritep, "./.");
                } else {
                  *cswritep++ = '.';
                }
              } else {
                const AlleleCode gt1_low = entry.dac1; // truncate
                const AlleleCode gt1_high = entry.dac1 >> (8 * sizeof(AlleleCode));
                cswritep = u32toa(gt1_low, cswritep);
                if (cur_autosomal_diploid || (gt1_low != gt1_high)) {
                  *cswritep++ = '/';
                  cswritep = u32toa(gt1_high, cswritep);
                }
              }
              *cswritep++ = '\t';
              if (entry.dac2 == kMissingDoubleAlleleCode) {
                if (cur_autosomal_diploid) {
                  cswritep = strcpya_k(cswritep, "./.");
                } else {
                  *cswritep++ = '.';
                }
              } else {
                const AlleleCode gt2_low = entry.dac2; // truncate
                const AlleleCode gt2_high = entry.dac2 >> (8 * sizeof(AlleleCode));
                cswritep = u32toa(gt2_low, cswritep);
                if (cur_autosomal_diploid || (gt2_low != gt2_high)) {
                  *cswritep++ = '/';
                  cswritep = u32toa(gt2_high, cswritep);
                }
              }
            } else {
              // dosage_reported
              *cswritep++ = '\t';
              const AlleleCode gt1_low = entry.dac1; // truncate
              const AlleleCode gt1_high = entry.dac1 >> (8 * sizeof(AlleleCode));
              if (cur_autosomal_diploid) {
                cswritep = PrintMultiallelicHcAsDs(gt1_low, gt1_high, merged_allele_ct, cswritep);
              } else {
                cswritep = PrintMultiallelicHcAsHaploidDs(gt1_low, gt1_high, merged_allele_ct, cswritep);
              }
              *cswritep++ = '\t';
              const AlleleCode gt2_low = entry.dac2; // truncate
              const AlleleCode gt2_high = entry.dac2 >> (8 * sizeof(AlleleCode));
              if (cur_autosomal_diploid) {
                cswritep = PrintMultiallelicHcAsDs(gt2_low, gt2_high, merged_allele_ct, cswritep);
              } else {
                cswritep = PrintMultiallelicHcAsHaploidDs(gt2_low, gt2_high, merged_allele_ct, cswritep);
              }
            }
          } else {
            // dosage_needed
            const Dosage* cur_ds_entry = &(ds_entries[diff_idx * (2 * k1LU) * merged_allele_ct_m1]);
            *cswritep++ = '\t';
            for (uint32_t uii = 0; uii != 2; ++uii) {
              if (cur_ds_entry[0] == kDosageMissing) {
                cswritep = strcpya_k(cswritep, ".\t");
              } else {
                if (cur_autosomal_diploid) {
                  for (uint32_t alt_idx = 0; alt_idx != merged_allele_ct_m1; ++alt_idx) {
                    cswritep = PrintSmallDosage(cur_ds_entry[alt_idx], cswritep);
                    *cswritep++ = ',';
                  }
                } else {
                  for (uint32_t alt_idx = 0; alt_idx != merged_allele_ct_m1; ++alt_idx) {
                    cswritep = PrintHaploidDosage(cur_ds_entry[alt_idx], cswritep);
                    *cswritep++ = ',';
                  }
                }
                cswritep[-1] = '\t';
              }
              cur_ds_entry = &(cur_ds_entry[merged_allele_ct_m1]);
            }
            --cswritep;
          }
        }
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto PgenDiff_ret_WRITE_FAIL;
        }
      }
    }
    // could verify we're at .pvar EOF

    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto PgenDiff_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    const uint32_t matched_variant_ct = PopcountWords(already_seen, raw_variant_ctl);
    logprintfww("--pgen-diff: %u sample%s and %u variant%s compared, %" PRIu64 " difference%s reported to %s .\n", sample_ct, (sample_ct == 1)? "" : "s", matched_variant_ct, (matched_variant_ct == 1)? "" : "s", grand_diff_ct, (grand_diff_ct == 1)? "" : "s", outname);

  }
  while (0) {
  PgenDiff_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  PgenDiff_ret_PVAR_REWIND_FAIL_N:
    logputs("\n");
  PgenDiff_ret_PVAR_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, "--pgen-diff .pvar file");
    reterr = kPglRetRewindFail;
    break;
  PgenDiff_ret_PVAR_MISSING_TOKENS_N:
    logputs("\n");
  PgenDiff_ret_PVAR_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of --pgen-diff .pvar file has fewer tokens than expected.\n", pvar_line_idx);
    reterr = kPglRetMalformedInput;
    break;
  PgenDiff_ret_PSAM_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&psam_txs)) {
      break;
    }
  PgenDiff_ret_PSAM_TSTREAM_FAIL:
    TextStreamErrPrint("--pgen-diff .psam file", &psam_txs);
    break;
  PgenDiff_ret_PSAM_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of --pgen-diff .psam file has fewer tokens than expected.\n", psam_line_idx);
    reterr = kPglRetMalformedInput;
    break;
  PgenDiff_ret_MALFORMED_INPUT_WW_N:
    logputs("\n");
  PgenDiff_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  PgenDiff_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  PgenDiff_ret_INCONSISTENT_INPUT_WW_N:
    logputs("\n");
  PgenDiff_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  PgenDiff_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  PgenDiff_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  PgenDiff_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 PgenDiff_ret_1:
  CswriteCloseCond(&css, cswritep);
  CleanupPgr2("--pgen-diff .pgen file", &simple_pgr2, &reterr);
  CleanupPgfi2("--pgen-diff .pgen file", &pgfi2, &reterr);
  CleanupTextStream2("--pgen-diff .pvar file", &pvar_txs, &reterr);
  CleanupTextStream2("--pgen-diff .psam file", &psam_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
