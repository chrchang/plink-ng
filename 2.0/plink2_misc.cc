// This file is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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


#include "include/plink2_stats.h"
#include "plink2_compress_stream.h"
#include "plink2_data.h"
#include "plink2_misc.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitUpdateSex(UpdateSexInfo* update_sex_info_ptr) {
  update_sex_info_ptr->flags = kfUpdateSex0;
  update_sex_info_ptr->col_num = 0;
  update_sex_info_ptr->fname = nullptr;
}

void CleanupUpdateSex(UpdateSexInfo* update_sex_info_ptr) {
  free_cond(update_sex_info_ptr->fname);
}

void InitSdiff(SdiffInfo* sdiff_info_ptr) {
  sdiff_info_ptr->flags = kfSdiff0;
  sdiff_info_ptr->dosage_hap_tol = kDosageMissing;
  sdiff_info_ptr->fname_id_delim = '\0';
  sdiff_info_ptr->other_id_ct = 0;
  sdiff_info_ptr->first_id_or_fname = nullptr;
  sdiff_info_ptr->other_ids_flattened = nullptr;
}

void CleanupSdiff(SdiffInfo* sdiff_info_ptr) {
  free_cond(sdiff_info_ptr->first_id_or_fname);
  free_cond(sdiff_info_ptr->other_ids_flattened);
}

void InitFst(FstInfo* fst_info_ptr) {
  fst_info_ptr->flags = kfFst0;
  fst_info_ptr->blocksize = 0;
  fst_info_ptr->other_id_ct = 0;
  fst_info_ptr->pheno_name = nullptr;
  fst_info_ptr->first_id_or_fname = nullptr;
  fst_info_ptr->other_ids_flattened = nullptr;
}

void CleanupFst(FstInfo* fst_info_ptr) {
  free_cond(fst_info_ptr->pheno_name);
  free_cond(fst_info_ptr->first_id_or_fname);
  free_cond(fst_info_ptr->other_ids_flattened);
}

PglErr UpdateVarBps(const ChrInfo* cip, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const TwoColParams* params, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uint32_t htable_size, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* __restrict variant_bps, uint32_t* __restrict variant_ct_ptr, UnsortedVar* vpos_sortstatusp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    uintptr_t* already_seen;
    if (unlikely(
            bigstack_calloc_w(BitCtToWordCt(raw_variant_ct), &already_seen))) {
      goto UpdateVarBps_ret_NOMEM;
    }
    // This could be pointed at a file containing allele codes, so don't limit
    // line length to minimum value.
    reterr = SizeAndInitTextStream(params->fname, bigstack_left(), MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto UpdateVarBps_ret_TSTREAM_FAIL;
    }
    reterr = TextSkip(params->skip_ct, &txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        snprintf(g_logbuf, kLogbufSize, "Error: Fewer lines than expected in %s.\n", params->fname);
        goto UpdateVarBps_ret_INCONSISTENT_INPUT_WW;
      }
      goto UpdateVarBps_ret_TSTREAM_FAIL;
    }
    line_idx = params->skip_ct;
    const uint32_t colid_first = (params->colid < params->colx);
    uint32_t variant_ct = *variant_ct_ptr;
    uint32_t colmin;
    uint32_t coldiff;
    if (colid_first) {
      colmin = params->colid - 1;
      coldiff = params->colx - params->colid;
    } else {
      colmin = params->colx - 1;
      coldiff = params->colid - params->colx;
    }
    const char skipchar = params->skipchar;
    uintptr_t miss_ct = 0;
    uint32_t hit_ct = 0;
    while (1) {
      ++line_idx;
      const char* line_start = TextGet(&txs);
      if (!line_start) {
        if (likely(!TextStreamErrcode2(&txs, &reterr))) {
          break;
        }
        goto UpdateVarBps_ret_TSTREAM_FAIL;
      }
      char cc = *line_start;
      if (cc == skipchar) {
        continue;
      }
      const char* colid_ptr;
      const char* colbp_ptr;
      if (colid_first) {
        colid_ptr = NextTokenMult0(line_start, colmin);
        colbp_ptr = NextTokenMult(colid_ptr, coldiff);
        if (unlikely(!colbp_ptr)) {
          goto UpdateVarBps_ret_MISSING_TOKENS;
        }
      } else {
        colbp_ptr = NextTokenMult0(line_start, colmin);
        colid_ptr = NextTokenMult(colbp_ptr, coldiff);
        if (unlikely(!colid_ptr)) {
          goto UpdateVarBps_ret_MISSING_TOKENS;
        }
      }
      const uint32_t varid_slen = strlen_se(colid_ptr);
      uint32_t cur_llidx;
      uint32_t variant_uidx = VariantIdDupHtableFind(colid_ptr, variant_ids, variant_id_htable, htable_dup_base, varid_slen, htable_size, max_variant_id_slen, &cur_llidx);
      if (variant_uidx == UINT32_MAX) {
        ++miss_ct;
        continue;
      }
      const char* cur_var_id = variant_ids[variant_uidx];
      if (unlikely(cur_llidx != UINT32_MAX)) {
        // we could check if some copies have been filtered out after hash
        // table construction?
        snprintf(g_logbuf, kLogbufSize, "Error: --update-map variant ID '%s' appears multiple times in dataset.\n", cur_var_id);
        goto UpdateVarBps_ret_INCONSISTENT_INPUT_WW;
      }
      if (unlikely(IsSet(already_seen, variant_uidx))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --update-map file.\n", cur_var_id);
        goto UpdateVarBps_ret_INCONSISTENT_INPUT_WW;
      }
      SetBit(variant_uidx, already_seen);
      if (!IsSet(variant_include, variant_uidx)) {
        continue;
      }
      int32_t bp_coord;
      if (ScanIntAbsDefcap(colbp_ptr, &bp_coord)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of --update-map file.\n", line_idx);
        goto UpdateVarBps_ret_MALFORMED_INPUT;
      }
      if (bp_coord < 0) {
        ClearBit(variant_uidx, variant_include);
        --variant_ct;
      } else {
        variant_bps[variant_uidx] = bp_coord;
      }
      ++hit_ct;
    }
    if (miss_ct) {
      snprintf(g_logbuf, kLogbufSize, "--update-map: %u value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      snprintf(g_logbuf, kLogbufSize, "--update-map: %u value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logputsb();

    UnsortedVar vpos_sortstatus = (*vpos_sortstatusp) & (~kfUnsortedVarBp);
    if (!(vpos_sortstatus & kfUnsortedVarSplitChr)) {
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t last_bp = 0;
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        if (variant_uidx >= chr_end) {
          do {
            ++chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          } while (variant_uidx >= chr_end);
          last_bp = 0;
        }
        const uint32_t cur_bp = variant_bps[variant_uidx];
        if (last_bp > cur_bp) {
          vpos_sortstatus |= kfUnsortedVarBp;
          if (!((*vpos_sortstatusp) & kfUnsortedVarBp)) {
            logerrputs("Warning: Base-pair positions are now unsorted!\n");
          }
          break;
        }
        last_bp = cur_bp;
      }
      if (((*vpos_sortstatusp) & kfUnsortedVarBp) && (!(vpos_sortstatus & kfUnsortedVarBp))) {
        logputs("Base-pair positions are now sorted.\n");
      }
      *vpos_sortstatusp = vpos_sortstatus;
    }
  }
  while (0) {
  UpdateVarBps_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  UpdateVarBps_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--update-map file", &txs);
    break;
  UpdateVarBps_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --update-map file has fewer tokens than expected.\n", line_idx);
  UpdateVarBps_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  UpdateVarBps_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
  CleanupTextStream2("--update-map file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr UpdateVarNames(const uintptr_t* variant_include, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const TwoColParams* params, uint32_t raw_variant_ct, uint32_t htable_size, uint32_t max_thread_ct, char** variant_ids, uint32_t* max_variant_id_slen_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    const uint32_t orig_max_variant_id_slen = *max_variant_id_slen_ptr;
    uint32_t max_variant_id_slen = orig_max_variant_id_slen;
    char** variant_ids_copy;
    uintptr_t* already_seen;
    if (unlikely(
            bigstack_alloc_cp(raw_variant_ct, &variant_ids_copy) ||
            bigstack_calloc_w(BitCtToWordCt(raw_variant_ct), &already_seen))) {
      goto UpdateVarNames_ret_NOMEM;
    }
    memcpy(variant_ids_copy, variant_ids, raw_variant_ct * sizeof(intptr_t));
    // This could be pointed at a file containing allele codes, so don't limit
    // line length to minimum value.
    // On the other hand, new variant IDs are allocated off the end of
    // bigstack, and that could result in a lot of memory pressure.
    reterr = SizeAndInitTextStream(params->fname, bigstack_left() / 4, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto UpdateVarNames_ret_TSTREAM_FAIL;
    }
    reterr = TextSkip(params->skip_ct, &txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        snprintf(g_logbuf, kLogbufSize, "Error: Fewer lines than expected in %s.\n", params->fname);
        goto UpdateVarNames_ret_INCONSISTENT_INPUT_WW;
      }
      goto UpdateVarNames_ret_TSTREAM_FAIL;
    }
    line_idx = params->skip_ct;
    const uint32_t colold_first = (params->colid < params->colx);
    uint32_t colmin;
    uint32_t coldiff;
    if (colold_first) {
      colmin = params->colid - 1;
      coldiff = params->colx - params->colid;
    } else {
      colmin = params->colx - 1;
      coldiff = params->colid - params->colx;
    }
    const char skipchar = params->skipchar;
    char* alloc_base = R_CAST(char*, g_bigstack_base);
    char* alloc_end = R_CAST(char*, g_bigstack_end);
    uintptr_t miss_ct = 0;
    uint32_t hit_ct = 0;
    while (1) {
      ++line_idx;
      const char* line_start = TextGet(&txs);
      if (!line_start) {
        if (likely(TextStreamErrcode2(&txs, &reterr))) {
          break;
        }
        goto UpdateVarNames_ret_TSTREAM_FAIL;
      }
      char cc = *line_start;
      if (cc == skipchar) {
        continue;
      }
      const char* colold_ptr;
      const char* colnew_ptr;
      if (colold_first) {
        colold_ptr = NextTokenMult0(line_start, colmin);
        colnew_ptr = NextTokenMult(colold_ptr, coldiff);
        if (unlikely(!colnew_ptr)) {
          goto UpdateVarNames_ret_MISSING_TOKENS;
        }
      } else {
        colnew_ptr = NextTokenMult0(line_start, colmin);
        colold_ptr = NextTokenMult(colnew_ptr, coldiff);
        if (unlikely(!colold_ptr)) {
          goto UpdateVarNames_ret_MISSING_TOKENS;
        }
      }
      const uint32_t colold_slen = strlen_se(colold_ptr);
      uint32_t cur_llidx;
      uint32_t variant_uidx = VariantIdDupHtableFind(colold_ptr, TO_CONSTCPCONSTP(variant_ids_copy), variant_id_htable, htable_dup_base, colold_slen, htable_size, orig_max_variant_id_slen, &cur_llidx);
      if (variant_uidx == UINT32_MAX) {
        ++miss_ct;
        continue;
      }
      const char* cur_var_id = variant_ids_copy[variant_uidx];
      if (unlikely(cur_llidx != UINT32_MAX)) {
        // we could check if some copies have been filtered out after hash
        // table construction?
        snprintf(g_logbuf, kLogbufSize, "Error: --update-name variant ID '%s' appears multiple times in dataset.\n", cur_var_id);
        goto UpdateVarNames_ret_INCONSISTENT_INPUT_WW;
      }
      if (!IsSet(variant_include, variant_uidx)) {
        continue;
      }
      if (unlikely(IsSet(already_seen, variant_uidx))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --update-name file.\n", cur_var_id);
        goto UpdateVarNames_ret_INCONSISTENT_INPUT_WW;
      }
      SetBit(variant_uidx, already_seen);
      ++hit_ct;
      const uint32_t colnew_slen = strlen_se(colnew_ptr);
      if (colnew_slen <= colold_slen) {
        const uint32_t colold_blen = colold_slen + 1;
        if (unlikely(S_CAST(uintptr_t, alloc_end - alloc_base) < colold_blen)) {
          goto UpdateVarNames_ret_NOMEM;
        }
        memcpy(alloc_base, cur_var_id, colold_blen);
        variant_ids_copy[variant_uidx] = alloc_base;
        alloc_base = &(alloc_base[colold_blen]);
        memcpyx(variant_ids[variant_uidx], colnew_ptr, colnew_slen, '\0');
      } else {
        if (colnew_slen > max_variant_id_slen) {
          max_variant_id_slen = colnew_slen;
        }
        const uint32_t colnew_blen = colnew_slen + 1;
        if (unlikely(S_CAST(uintptr_t, alloc_end - alloc_base) < colnew_blen)) {
          goto UpdateVarNames_ret_NOMEM;
        }
        alloc_end -= colnew_blen;
        memcpyx(alloc_end, colnew_ptr, colnew_slen, '\0');
        variant_ids[variant_uidx] = alloc_end;
      }
    }
    BigstackEndSet(alloc_end);
    if (miss_ct) {
      snprintf(g_logbuf, kLogbufSize, "--update-name: %u value%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      snprintf(g_logbuf, kLogbufSize, "--update-name: %u value%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logputsb();
    *max_variant_id_slen_ptr = max_variant_id_slen;
  }
  while (0) {
  UpdateVarNames_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  UpdateVarNames_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--update-name file", &txs);
    break;
  UpdateVarNames_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --update-name file has fewer tokens than expected.\n", line_idx);
  UpdateVarNames_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
  CleanupTextStream2("--update-name file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr UpdateVarAlleles(const char* fname, const uintptr_t* variant_include, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const uintptr_t* allele_idx_offsets, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uint32_t htable_size, uint32_t max_thread_ct, char** allele_storage_mutable, uint32_t* max_allele_slen_ptr, char* outname, char* outname_end) {
  // probable todos:
  // - add '3col' modifier to support three-column input file where second
  //   column has comma-delimited old alleles and third column has
  //   comma-delimited new alleles
  // - add 'allow-absent' modifier to do partial-replace instead of writing to
  //   .allele.no.snp when some old allele entries match loaded alleles and
  //   some don't
  // - add 'strict-missing' modifier to remove missing code -> remaining allele
  //   behavior
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  FILE* errfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    uintptr_t* already_seen;
    if (unlikely(bigstack_calloc_w(BitCtToWordCt(raw_variant_ct), &already_seen))) {
      goto UpdateVarAlleles_ret_NOMEM;
    }
    reterr = SizeAndInitTextStream(fname, bigstack_left() / 4, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto UpdateVarAlleles_ret_TSTREAM_FAIL;
    }
    const char* std_input_missing_geno = &(g_one_char_strs[92]);
    const char input_missing_geno_char = *g_input_missing_geno_ptr;
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = g_bigstack_end;
    uintptr_t miss_ct = 0;
    uintptr_t err_ct = 0;
    uint32_t hit_ct = 0;
    uint32_t max_allele_slen = *max_allele_slen_ptr;
    uint32_t cur_allele_ct = 2;
    const char* line_iter = TextLineEnd(&txs);
    ++line_idx;
    for (; TextGetUnsafe2K(&txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      const char* varid_start = line_iter;
      const char* varid_end = CurTokenEnd(varid_start);
      uint32_t cur_llidx;
      uint32_t variant_uidx = VariantIdDupHtableFind(varid_start, variant_ids, variant_id_htable, htable_dup_base, varid_end - varid_start, htable_size, max_variant_id_slen, &cur_llidx);
      if (variant_uidx == UINT32_MAX) {
        line_iter = varid_end;
        ++miss_ct;
        continue;
      }
      const char* varid = variant_ids[variant_uidx];
      if (unlikely(cur_llidx != UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --update-alleles variant ID '%s' appears multiple times in dataset.\n", varid);
        goto UpdateVarAlleles_ret_INCONSISTENT_INPUT_WW;
      }
      if (unlikely(IsSet(already_seen, variant_uidx))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --update-alleles file.\n", varid);
        goto UpdateVarAlleles_ret_INCONSISTENT_INPUT_WW;
      }
      SetBit(variant_uidx, already_seen);
      if (!IsSet(variant_include, variant_uidx)) {
        line_iter = varid_end;
        ++miss_ct;
        continue;
      }
      // change these to arrays once new input format is supported
      const char* old_allele1_start = FirstNonTspace(varid_end);
      const char* old_allele1_end = FirstSpaceOrEoln(old_allele1_start);
      const char* old_allele2_start = FirstNonTspace(old_allele1_end);
      const char* old_allele2_end = FirstSpaceOrEoln(old_allele2_start);
      const char* new_allele1_start = FirstNonTspace(old_allele2_end);
      const char* new_allele1_end = FirstSpaceOrEoln(new_allele1_start);
      const char* new_allele2_start = FirstNonTspace(new_allele1_end);
      if (IsEolnKns(*new_allele2_start)) {
        goto UpdateVarAlleles_ret_MISSING_TOKENS;
      }
      const char* new_allele2_end = FirstSpaceOrEoln(new_allele2_start);
      line_iter = new_allele2_end;
      const uint32_t old_slen1 = old_allele1_end - old_allele1_start;
      const uint32_t old_slen2 = old_allele2_end - old_allele2_start;
      // standardize missing allele codes before proceeding
      if ((old_slen1 == 1) && (old_allele1_start[0] == input_missing_geno_char)) {
        old_allele1_start = std_input_missing_geno;
      }
      if ((old_slen2 == 1) && (old_allele2_start[0] == input_missing_geno_char)) {
        old_allele2_start = std_input_missing_geno;
      }
      if ((old_slen1 == old_slen2) && memequal(old_allele1_start, old_allele2_start, old_slen1)) {
        goto UpdateVarAlleles_ret_DUPLICATE_INPUT_ALLELE_CODE;
      }
      const uint32_t new_slen1 = new_allele1_end - new_allele1_start;
      const uint32_t new_slen2 = new_allele2_end - new_allele2_start;
      if ((new_slen1 == 1) && (new_allele1_start[0] == input_missing_geno_char)) {
        new_allele1_start = std_input_missing_geno;
      }
      if ((new_slen2 == 1) && (new_allele2_start[0] == input_missing_geno_char)) {
        new_allele2_start = std_input_missing_geno;
      }
      if ((new_slen1 == new_slen2) && memequal(new_allele1_start, new_allele2_start, new_slen1)) {
        goto UpdateVarAlleles_ret_DUPLICATE_INPUT_ALLELE_CODE;
      }
      if (memchr(new_allele1_start, ',', new_slen1) || memchr(new_allele2_start, ',', new_slen2)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Comma-containing new allele code on line %" PRIuPTR " of --update-alleles file.\n", line_idx);
        goto UpdateVarAlleles_ret_MALFORMED_INPUT_WW;
      }
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      char** cur_alleles = &(allele_storage_mutable[allele_idx_offset_base]);
      uint32_t old_allele1_match = UINT32_MAX;
      uint32_t old_allele2_match = UINT32_MAX;
      for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
        const char* old_allele = cur_alleles[allele_idx];
        if (memequal(old_allele1_start, old_allele, old_slen1) && (!old_allele[old_slen1])) {
          if (old_allele1_match != UINT32_MAX) {
            snprintf(g_logbuf, kLogbufSize, "Error: Duplicate allele code in variant '%s'.\n", varid);
            goto UpdateVarAlleles_ret_MALFORMED_INPUT_WW;
          }
          old_allele1_match = allele_idx;
        } else if (memequal(old_allele2_start, old_allele, old_slen2) && (!old_allele[old_slen2])) {
          if (old_allele2_match != UINT32_MAX) {
            snprintf(g_logbuf, kLogbufSize, "Error: Duplicate allele code in variant '%s'.\n", varid);
            goto UpdateVarAlleles_ret_MALFORMED_INPUT_WW;
          }
          old_allele2_match = allele_idx;
        }
      }
      if (old_allele1_match == UINT32_MAX) {
        // Allow one missing allele code which doesn't match a loaded missing
        // allele code when it doesn't introduce ambiguity.
        // I.e. since we're only changing two alleles for now, the variant must
        // be biallelic, and the other old_allele must match one of the loaded
        // alleles.
        if ((old_allele2_match == UINT32_MAX) || (cur_allele_ct != 2) || (old_slen1 != 1) || (old_allele1_start[0] != '.')) {
        UpdateVarAlleles_errfile:
          if (!err_ct) {
            strcpy(outname_end, ".allele.no.snp");
            if (fopen_checked(outname, FOPEN_WB, &errfile)) {
              goto UpdateVarAlleles_ret_OPEN_FAIL;
            }
          }
          // variant ID, expected allele 1, expected allele 2
          // not optimized for now, could explicitly manage a write buffer if
          // it ever matters
          fputs(varid, errfile);
          putc_unlocked('\t', errfile);
          fwrite(old_allele1_start, old_slen1, 1, errfile);
          putc_unlocked('\t', errfile);
          fwrite(old_allele2_start, old_slen2, 1, errfile);
#ifdef _WIN32
          putc_unlocked('\r', errfile);
#endif
          if (unlikely(putc_checked('\n', errfile))) {
            goto UpdateVarAlleles_ret_WRITE_FAIL;
          }
          ++err_ct;
          continue;
        }
        old_allele1_match = 1 - old_allele2_match;
      }
      if (old_allele2_match == UINT32_MAX) {
        if ((cur_allele_ct != 2) || (old_slen2 != 1) || (old_allele2_start[0] != '.')) {
          goto UpdateVarAlleles_errfile;
        }
        old_allele2_match = 1 - old_allele1_match;
      }
      if (new_slen1 == 1) {
        cur_alleles[old_allele1_match] = K_CAST(char*, &(g_one_char_strs[ctou32(new_allele1_start[0]) * 2]));
      } else {
        // reuse old storage if we can, allocate when we must
        if (old_slen1 < new_slen1) {
          if (PtrWSubCk(tmp_alloc_base, new_slen1 + 1, &tmp_alloc_end)) {
            goto UpdateVarAlleles_ret_NOMEM;
          }
          cur_alleles[old_allele1_match] = R_CAST(char*, tmp_alloc_end);
        }
        memcpyx(cur_alleles[old_allele1_match], new_allele1_start, new_slen1, '\0');
      }
      if (new_slen2 == 1) {
        cur_alleles[old_allele2_match] = K_CAST(char*, &(g_one_char_strs[ctou32(new_allele2_start[0]) * 2]));
      } else {
        if (old_slen2 < new_slen2) {
          if (PtrWSubCk(tmp_alloc_base, new_slen2 + 1, &tmp_alloc_end)) {
            goto UpdateVarAlleles_ret_NOMEM;
          }
          cur_alleles[old_allele2_match] = R_CAST(char*, tmp_alloc_end);
        }
        memcpyx(cur_alleles[old_allele2_match], new_allele2_start, new_slen2, '\0');
      }
      if (cur_allele_ct > 2) {
        // verify we aren't creating a duplicate allele
        // todo: create a tiny hash table instead (upfront), this will be
        // especially important if we're updating more than two alleles
        for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
          if ((allele_idx == old_allele1_match) || (allele_idx == old_allele2_match)) {
            continue;
          }
          const char* cur_allele = cur_alleles[allele_idx];
          if ((memequal(new_allele1_start, cur_allele, new_slen1) && (!cur_allele[new_slen1])) ||
              (memequal(new_allele2_start, cur_allele, new_slen2) && (!cur_allele[new_slen2]))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Duplicate allele code in variant '%s'.\n", varid);
            goto UpdateVarAlleles_ret_MALFORMED_INPUT_WW;
          }
        }
      }
      ++hit_ct;
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto UpdateVarAlleles_ret_TSTREAM_FAIL;
    }
    *max_allele_slen_ptr = max_allele_slen;
    if (miss_ct) {
      snprintf(g_logbuf, kLogbufSize, "--update-alleles: %u variant%s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      snprintf(g_logbuf, kLogbufSize, "--update-alleles: %u variant%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logputsb();
    if (err_ct) {
      logprintfww("%" PRIuPTR " update failure%s logged to %s .\n", err_ct, (err_ct == 1)? "" : "s", outname);
    }
    BigstackEndSet(tmp_alloc_end);
  }
  while (0) {
  UpdateVarAlleles_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  UpdateVarAlleles_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  UpdateVarAlleles_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  UpdateVarAlleles_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--update-alleles file", &txs);
    break;
  UpdateVarAlleles_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of --update-alleles file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  UpdateVarAlleles_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  UpdateVarAlleles_ret_DUPLICATE_INPUT_ALLELE_CODE:
    snprintf(g_logbuf, kLogbufSize, "Error: Duplicate allele code on line %" PRIuPTR " of --update-alleles file.\n", line_idx);
  UpdateVarAlleles_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  }
  fclose_cond(errfile);
  CleanupTextStream2("--update-alleles file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr RecoverVarIds(const char* fname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const char* missing_varid, uint32_t raw_variant_ct, uint32_t variant_ct, RecoverVarIdsFlags flags, uint32_t max_thread_ct, char** variant_ids, uint32_t* max_variant_id_slen_ptr, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* alloc_end = R_CAST(char*, g_bigstack_end);
  uintptr_t line_idx = 0;
  FILE* errfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    // Significant overlap with UpdateVarNames().
    const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    char** variant_ids_copy;
    uintptr_t* chr_already_seen;
    uintptr_t* already_seen;
    uintptr_t* conflict_bitarr;
    char** orig_alt_starts;
    if (unlikely(
            bigstack_alloc_cp(raw_variant_ct, &variant_ids_copy) ||
            bigstack_calloc_w(BitCtToWordCt(chr_code_end), &chr_already_seen) ||
            bigstack_calloc_w(raw_variant_ctl, &already_seen) ||
            bigstack_calloc_w(raw_variant_ctl, &conflict_bitarr) ||
            bigstack_alloc_cp(kPglMaxAltAlleleCt, &orig_alt_starts))) {
      goto RecoverVarIds_ret_NOMEM;
    }
    // need this to be able to write original variant IDs in conflict case
    memcpy(variant_ids_copy, variant_ids, raw_variant_ct * sizeof(intptr_t));

    reterr = SizeAndInitTextStream(fname, bigstack_left() / 4, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto RecoverVarIds_ret_TSTREAM_FAIL;
    }
    char* line_start;
    do {
      ++line_idx;
      line_start = TextGet(&txs);
      if (unlikely(!line_start)) {
        if (!TextStreamErrcode2(&txs, &reterr)) {
          logerrputs("Error: Empty --recover-var-ids file.\n");
          goto RecoverVarIds_ret_1;
        }
        goto RecoverVarIds_ret_TSTREAM_FAIL;
      }
    } while ((*line_start == '#') && (!tokequal_k(line_start, "#CHROM")));
    // simplified copy of plink2_pvar code, probably want to write a library
    // function for this
    uint32_t col_skips[4];
    uint32_t col_types[4];
    uint32_t strict_allele_order = 1;
    const uint32_t is_bim = !(line_start[0] == '#');
    char* line_iter = nullptr;
    if (!is_bim) {
      // parse header
      // [-1] = #CHROM (must be first column)
      // [0] = POS
      // [1] = ID
      // [2] = REF
      // [3] = ALT
      char* token_end = &(line_start[6]);
      uint32_t found_header_bitset = 0;
      uint32_t relevant_postchr_col_ct = 0;
      char* linebuf_iter;
      for (uint32_t col_idx = 1; ; ++col_idx) {
        linebuf_iter = FirstNonTspace(token_end);
        if (IsEolnKns(*linebuf_iter)) {
          break;
        }
        token_end = CurTokenEnd(linebuf_iter);
        const uint32_t token_slen = token_end - linebuf_iter;
        uint32_t cur_col_type;
        if (token_slen == 3) {
          if (memequal_k(linebuf_iter, "POS", 3)) {
            cur_col_type = 0;
          } else if (memequal_k(linebuf_iter, "REF", 3)) {
            cur_col_type = 2;
          } else if (memequal_k(linebuf_iter, "ALT", 3)) {
            cur_col_type = 3;
          } else {
            continue;
          }
        } else if (token_slen == 2) {
          if (memequal_k(linebuf_iter, "ID", 2)) {
            cur_col_type = 1;
          } else {
            continue;
          }
        } else {
          continue;
        }
        const uint32_t cur_col_type_shifted = 1 << cur_col_type;
        if (unlikely(found_header_bitset & cur_col_type_shifted)) {
          // known token, so no overflow danger
          char* write_iter = strcpya_k(g_logbuf, "Error: Duplicate column header '");
          write_iter = memcpya(write_iter, linebuf_iter, token_slen);
          write_iter = strcpya_k(write_iter, "' on line ");
          write_iter = wtoa(line_idx, write_iter);
          strcpy_k(write_iter, " of --recover-var-ids file.\n");
          goto RecoverVarIds_ret_MALFORMED_INPUT_WW;
        }
        found_header_bitset |= cur_col_type_shifted;
        col_skips[relevant_postchr_col_ct] = col_idx;
        col_types[relevant_postchr_col_ct++] = cur_col_type;
      }
      if (unlikely(relevant_postchr_col_ct != 4)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Missing column header(s) on line %" PRIuPTR " of --recover-var-ids file. (POS, ID, REF, and ALT are required.)\n", line_idx);
        goto RecoverVarIds_ret_MALFORMED_INPUT_WW;
      }
      for (uint32_t rpc_col_idx = 3; rpc_col_idx; --rpc_col_idx) {
        col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
      }
      line_iter = AdvPastDelim(linebuf_iter, '\n');
    } else {
      // .bim.  Interpret as #CHROM ID POS ALT REF if there are exactly 5
      // columns, otherwise #CHROM ID CM POS ALT REF.
      col_skips[0] = 1;
      col_skips[2] = 1;
      col_skips[3] = 1;
      col_types[0] = 1;
      col_types[1] = 0;
      col_types[2] = 3;
      col_types[3] = 2;
      char* linebuf_iter = NextTokenMult(line_start, 4);
      if (unlikely(!linebuf_iter)) {
        goto RecoverVarIds_ret_MISSING_TOKENS;
      }
      linebuf_iter = NextToken(linebuf_iter);
      if (!linebuf_iter) {
        col_skips[1] = 1;
      } else {
        col_skips[1] = 2;
      }
      strict_allele_order = (flags / kfRecoverVarIdsStrictBimOrder) & 1;
      line_iter = line_start;
    }
    // Only nonzero if we're actually replacing existing IDs with the missing
    // ID in the conflict case.
    uint32_t missing_varid_blen = 0;
    if (flags & kfRecoverVarIdsForce) {
      if (!missing_varid) {
        missing_varid = &(g_one_char_strs[92]);  // '.'
      }
      missing_varid_blen = 1 + strlen(missing_varid);
    }
    const uint32_t is_rigid = (flags / kfRecoverVarIdsRigid) & 1;
    uint32_t max_variant_id_slen = *max_variant_id_slen_ptr;
    char* alloc_base = R_CAST(char*, g_bigstack_base);
    uint32_t is_unsorted = 0;
    uint32_t prev_chr = UINT32_MAX;
    uint32_t prev_chr_vidx_end = 0;
    uint32_t prev_chr_vidx_start = 0;
    uint32_t prev_bp = 0;
    uint32_t prev_vidx_start = 0;
    uint32_t prev_vidx_end = 0;
    uint32_t rm_dup_warning = 0;
    uintptr_t record_uidx = ~k0LU;  // deliberate overflow
    uint32_t cur_allele_ct = 2;
    for (; TextGetUnsafe2(&txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n')) {
      if (!(line_idx % 1000000)) {
        printf("\r--recover-var-ids: %" PRIuPTR "m lines scanned.", line_idx / 1000000);
        fflush(stdout);
      }
      ++line_idx;
      line_start = line_iter;
      if (unlikely(line_start[0] == '#')) {
        putc_unlocked('\n', stdout);
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --recover-var-ids file starts with a '#'. (This is only permitted before the first nonheader line, and if a #CHROM header line is present it must denote the end of the header block.)\n", line_idx);
        goto RecoverVarIds_ret_MALFORMED_INPUT_WW;
      }
      ++record_uidx;
      line_iter = CurTokenEnd(line_start);
      uint32_t cur_chr_code = GetChrCodeCounted(cip, line_iter - line_start, line_start);
      if (IsI32Neg(cur_chr_code)) {
        continue;
      }
      // Could add some special handling of chrX/PAR1/PAR2(/XY?) later, if
      // that proves to be a pain point; this does operate on VCF input,
      // after all.  It's a bit annoying to implement, though, so I'm
      // excluding it from the first version.
      if (cur_chr_code == prev_chr) {
        if (!prev_chr_vidx_end) {
          continue;
        }
      } else {
        prev_chr = cur_chr_code;
        if (!is_unsorted) {
          if (IsSet(chr_already_seen, cur_chr_code)) {
            is_unsorted = 1;
          } else {
            SetBit(cur_chr_code, chr_already_seen);
            prev_bp = 0;
          }
        }
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          prev_chr_vidx_end = 0;
          continue;
        }
        const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[cur_chr_code];
        prev_chr_vidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
        prev_chr_vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        prev_vidx_start = prev_chr_vidx_start;
        prev_vidx_end = prev_chr_vidx_start;
      }
      char* token_ptrs[4];
      uint32_t token_slens[4];
      line_iter = TokenLex(line_iter, col_types, col_skips, 4, token_ptrs, token_slens);
      if (unlikely(!line_iter)) {
        putc_unlocked('\n', stdout);
        goto RecoverVarIds_ret_MISSING_TOKENS;
      }
      // Usually, everything is sorted by position, so we use exponential
      // search over binary search unless we've detected unsorted data.
      int32_t cur_bp_signed;
      if (unlikely(ScanIntAbsDefcap(token_ptrs[0], &cur_bp_signed))) {
        putc_unlocked('\n', stdout);
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of --recover-var-ids file.\n", line_idx);
        goto RecoverVarIds_ret_MALFORMED_INPUT_WW;
      }
      if (cur_bp_signed < 0) {
        continue;
      }
      const uint32_t cur_bp = cur_bp_signed;
      if (cur_bp < prev_bp) {
        is_unsorted = 1;
      }
      if (!is_unsorted) {
        if (prev_vidx_end == prev_chr_vidx_end) {
          continue;
        }
        if (is_rigid || (cur_bp > prev_bp)) {
          const uint32_t arr_length = prev_chr_vidx_end - prev_vidx_end;
          const uint32_t incr = ExpsearchU32(&(variant_bps[prev_vidx_end]), arr_length, cur_bp);
          if (incr == arr_length) {
            prev_vidx_end = prev_chr_vidx_end;
            continue;
          }
          prev_vidx_start = incr + prev_vidx_end;
        }
        prev_bp = cur_bp;
      } else {
        const uint32_t arr_length = prev_chr_vidx_end - prev_chr_vidx_start;
        const uint32_t incr = CountSortedSmallerU32(&(variant_bps[prev_chr_vidx_start]), arr_length, cur_bp);
        if (incr == arr_length) {
          continue;
        }
        prev_vidx_start = incr + prev_chr_vidx_start;
      }
      prev_vidx_start = AdvBoundedTo1Bit(variant_include, prev_vidx_start, prev_chr_vidx_end);

      // variant_bps[prev_vidx_start] is the first element >= cur_bp that
      // hasn't been filtered out.
      if (variant_bps[prev_vidx_start] != cur_bp) {
        prev_vidx_end = prev_vidx_start;
        continue;
      }
      if (is_rigid) {
        if (unlikely(record_uidx >= raw_variant_ct)) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: \"--recover-var-ids rigid\" file does not match the main dataset.\n");
          goto RecoverVarIds_ret_INCONSISTENT_INPUT;
        }
        if (!IsSet(variant_include, record_uidx)) {
          continue;
        }
        if (record_uidx != prev_vidx_start) {
          putc_unlocked('\n', stdout);
          logerrputs("Error: \"--recover-var-ids rigid\" file does not match the main dataset.\n");
          goto RecoverVarIds_ret_INCONSISTENT_INPUT;
        }
      }
      char* orig_alt_start = token_ptrs[3];
      const uint32_t orig_alt_slen = token_slens[3];
      const uint32_t extra_orig_alt_ct = CountByte(orig_alt_start, ',', orig_alt_slen);
      if (unlikely(extra_orig_alt_ct && is_bim)) {
        // probably want to enforce this in main .pvar loader too...
        putc_unlocked('\n', stdout);
        snprintf(g_logbuf, kLogbufSize, "Error: Multiple ALT alleles on line %" PRIuPTR " of headerless --recover-var-ids file.\n", line_idx);
        goto RecoverVarIds_ret_MALFORMED_INPUT_WW;
      }
      // Null-terminate all allele codes so we can use strequal_overread();
      // and make sure to reset the final null terminator so that unsafe
      // line_iter advancement works.
      const char line_iter_char = *line_iter;
      char* orig_ref_start = token_ptrs[2];
      orig_ref_start[token_slens[2]] = '\0';
      {
        char* orig_alt_iter = orig_alt_start;
        for (uint32_t alt_idx = 0; alt_idx != extra_orig_alt_ct; ++alt_idx) {
          orig_alt_starts[alt_idx] = orig_alt_iter;
          orig_alt_iter = AdvToDelim(orig_alt_iter, ',');
          *orig_alt_iter++ = '\0';
        }
        orig_alt_starts[extra_orig_alt_ct] = orig_alt_iter;
      }
      char* orig_alt_end = &(orig_alt_start[orig_alt_slen]);
      *orig_alt_end = '\0';
      const uint32_t orig_alt_ct = extra_orig_alt_ct + 1;

      uint32_t cur_vidx = prev_vidx_start;
      // Multiple variants in the current dataset may have this position;
      // check them all for allele-code concordance.  If there are multiple
      // matches, allow this but print a warning suggesting --rm-dup.  In the
      // usual subcase where there are multiple instances in the original
      // file, we'll error out instead unless 'force' was specified.
      uint32_t already_matched = 0;
      goto RecoverVarIds_first_iter;
      for (; (!is_rigid) && (cur_vidx != prev_chr_vidx_end) && (variant_bps[cur_vidx] == cur_bp); cur_vidx = AdvBoundedTo1Bit(variant_include, cur_vidx + 1, prev_chr_vidx_end)) {
      RecoverVarIds_first_iter: ;
        uintptr_t allele_idx_offset_base = cur_vidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[cur_vidx];
          cur_allele_ct = allele_idx_offsets[cur_vidx + 1] - allele_idx_offset_base;
        }
        if (cur_allele_ct != orig_alt_ct + 1) {
          continue;
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        const char* ref_allele = cur_alleles[0];
        if (strict_allele_order) {
          if (!strequal_overread(ref_allele, orig_ref_start)) {
            continue;
          }
          const char* const* cur_alt_alleles = &(cur_alleles[1]);
          uint32_t alt_idx = 0;
          for (; alt_idx != orig_alt_ct; ++alt_idx) {
            if (!strequal_overread(cur_alt_alleles[alt_idx], orig_alt_starts[alt_idx])) {
              break;
            }
          }
          if (alt_idx != orig_alt_ct) {
            continue;
          }
        } else {
          // cur_allele_ct == 2 guaranteed.
          const char* other_allele = cur_alleles[1];
          if (!strequal_overread(ref_allele, orig_ref_start)) {
            if ((!strequal_overread(ref_allele, orig_alt_start)) || (!strequal_overread(other_allele, orig_ref_start))) {
              continue;
            }
          } else if (!strequal_overread(other_allele, orig_alt_start)) {
            continue;
          }
        }
        char* old_id = variant_ids[cur_vidx];
        const uint32_t old_id_slen = strlen(old_id);
        const char* new_id = token_ptrs[1];
        const uint32_t new_id_slen = token_slens[1];
        if (IsSet(already_seen, cur_vidx)) {
          if ((new_id_slen != old_id_slen) || (!memequal(new_id, old_id, old_id_slen))) {
            SetBit(cur_vidx, conflict_bitarr);
            rm_dup_warning = 1;
          }
          continue;
        }
        // Okay, chr/pos/alleles match, and there's no conflict with a
        // previous entry.  Update the variant ID.
        SetBit(cur_vidx, already_seen);
        if (new_id_slen <= old_id_slen) {
          const uint32_t old_id_blen = old_id_slen + 1;
          if (unlikely(S_CAST(uintptr_t, alloc_end - alloc_base) < old_id_blen)) {
            goto RecoverVarIds_ret_NOMEM;
          }
          memcpy(alloc_base, old_id, old_id_blen);
          variant_ids_copy[cur_vidx] = alloc_base;
          alloc_base = &(alloc_base[old_id_blen]);
          memcpyx(old_id, new_id, new_id_slen, '\0');
        } else {
          if (new_id_slen > max_variant_id_slen) {
            max_variant_id_slen = new_id_slen;
          }
          const uint32_t new_id_blen = new_id_slen + 1;
          if (unlikely(S_CAST(uintptr_t, alloc_end - alloc_base) < new_id_blen)) {
            goto RecoverVarIds_ret_NOMEM;
          }
          alloc_end -= new_id_blen;
          memcpyx(alloc_end, new_id, new_id_slen, '\0');
          variant_ids[cur_vidx] = alloc_end;
        }
        rm_dup_warning |= already_matched;
        already_matched = 1;
      }
      *line_iter = line_iter_char;
      prev_vidx_end = cur_vidx;
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto RecoverVarIds_ret_TSTREAM_FAIL;
    }
    putc_unlocked('\r', stdout);
    logprintf("--recover-var-ids: %" PRIuPTR " line%s scanned.\n", line_idx, (line_idx == 1)? "" : "s");
    const uint32_t conflict_ct = PopcountWords(conflict_bitarr, raw_variant_ctl);
    if (conflict_ct) {
      strcpy(outname_end, ".recoverid.dup");
      if (fopen_checked(outname, FOPEN_WB, &errfile)) {
        goto RecoverVarIds_ret_OPEN_FAIL;
      }
      // Simply a list of affected variant IDs, one per line, sorted
      // in order of original-file appearance.
      // The reported variant IDs will usually be duplicated, since
      // --set-all-var-ids makes them a function of CHROM/POS/REF/ALT, and
      // that's what old-variant-ID recovery is based on.  (But one could e.g.
      // set the IDs to zero-based indexes before calling --recover-var-ids, if
      // they wanted every .recoverid.dup entry to be unambiguous.)
      char* textbuf = g_textbuf;
      char* textbuf_flush = &(textbuf[kMaxMediumLine]);
      char* write_iter = textbuf;
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = conflict_bitarr[0];
      for (uint32_t conflict_idx = 0; conflict_idx != conflict_ct; ++conflict_idx) {
        const uint32_t variant_uidx = BitIter1(conflict_bitarr, &variant_uidx_base, &cur_bits);
        write_iter = strcpya(write_iter, variant_ids_copy[variant_uidx]);
        AppendBinaryEoln(&write_iter);
        if (unlikely(fwrite_ck(textbuf_flush, errfile, &write_iter))) {
          goto RecoverVarIds_ret_WRITE_FAIL;
        }
        if (missing_varid_blen) {
          char* cur_id = variant_ids[variant_uidx];
          const uint32_t cur_id_slen = strlen(cur_id);
          if (missing_varid_blen <= cur_id_slen + 1) {
            memcpy(cur_id, missing_varid, missing_varid_blen);
          } else {
            if (unlikely(S_CAST(uintptr_t, alloc_end - alloc_base) < missing_varid_blen)) {
              goto RecoverVarIds_ret_NOMEM;
            }
            alloc_end -= missing_varid_blen;
            memcpy(alloc_end, missing_varid, missing_varid_blen);
            variant_ids[variant_uidx] = alloc_end;
          }
        }
      }
      if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &errfile))) {
        goto RecoverVarIds_ret_WRITE_FAIL;
      }
      if (!missing_varid_blen) {
        // not really unlikely, conditional on being in this branch...
        snprintf(g_logbuf, kLogbufSize, "Error: %u variant%s had conflicting matching-position-and-alleles records in the --recover-var-ids file; affected ID%s been written to %s . (Add the 'force' modifier if you just want to set %s to the missing code.)\n", conflict_ct, (conflict_ct == 1)? "" : "s", (conflict_ct == 1)? "" : "s", outname, (conflict_ct == 1)? "it" : "them");
        goto RecoverVarIds_ret_INCONSISTENT_INPUT_WW;
      }
      logerrprintfww("Warning: %u variant ID%s set to '%s' by \"--recover-var-ids force\". Previous ID%s been written to %s .\n", conflict_ct, (conflict_ct == 1)? "" : "s", missing_varid, (conflict_ct == 1)? " has" : "s have", outname);
    }
    const uint32_t update_ct = PopcountWords(already_seen, raw_variant_ctl);
    const uint32_t unupdated_ct = variant_ct - update_ct;
    if (unupdated_ct && (!(flags & kfRecoverVarIdsPartial))) {
      snprintf(g_logbuf, kLogbufSize, "Error: %u/%u variant%s had no matching-position-and-alleles records in the --recover-var-ids file. (Add the 'partial' modifier when this is expected.)\n", unupdated_ct, variant_ct, (unupdated_ct == 1)? "" : "s");
      goto RecoverVarIds_ret_INCONSISTENT_INPUT_WW;
    }
    logprintf("--recover-var-ids: %u/%u ID%s updated.\n", update_ct, variant_ct, (update_ct == 1)? "" : "s");
    if (rm_dup_warning) {
      logerrputs("Warning: Some variants had identical positions and alleles, and therefore were\nassigned the same ID.  Consider deduplicating your data with --rm-dup.\n");
    }
    *max_variant_id_slen_ptr = max_variant_id_slen;
  }
  while (0) {
  RecoverVarIds_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  RecoverVarIds_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  RecoverVarIds_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  RecoverVarIds_ret_TSTREAM_FAIL:
    putc_unlocked('\n', stdout);
    TextStreamErrPrint("--recover-var-ids file", &txs);
    break;
  RecoverVarIds_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --recover-var-ids file has fewer tokens than expected.\n", line_idx);
  RecoverVarIds_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  RecoverVarIds_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  RecoverVarIds_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 RecoverVarIds_ret_1:
  CleanupTextStream2("--recover-var-ids file", &txs, &reterr);
  fclose_cond(errfile);

  // doesn't hurt to ensure all replaced variant_ids[] entries are valid in the
  // error case...
  BigstackEndSet(alloc_end);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr Plink1ClusterImport(const char* within_fname, const char* catpheno_name, const char* family_missing_catname, const uintptr_t* sample_include, const char* sample_ids, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t mwithin_val, uint32_t max_thread_ct, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream within_txs;
  PreinitTextStream(&within_txs);
  {
    if (!sample_ct) {
      goto Plink1ClusterImport_ret_1;
    }
    const char catpheno_name_default[] = "CATPHENO";
    uint32_t catpheno_name_blen;
    if (!catpheno_name) {
      catpheno_name = catpheno_name_default;
      catpheno_name_blen = 9;
    } else {
      catpheno_name_blen = 1 + strlen(catpheno_name);
    }
    const uintptr_t old_max_pheno_name_blen = *max_pheno_name_blen_ptr;
    const uint32_t old_pheno_ct = *pheno_ct_ptr;
    const char* old_pheno_names = *pheno_names_ptr;
    uintptr_t new_max_pheno_name_blen;
    if (old_pheno_names && (catpheno_name_blen <= old_max_pheno_name_blen)) {
      new_max_pheno_name_blen = old_max_pheno_name_blen;
      for (uint32_t pheno_idx = 0; pheno_idx != old_pheno_ct; ++pheno_idx) {
        if (unlikely(memequal(catpheno_name, &(old_pheno_names[pheno_idx * old_max_pheno_name_blen]), catpheno_name_blen))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Cannot create a new categorical phenotype named '%s', since another phenotype of the same name already exists.\n", catpheno_name);
          goto Plink1ClusterImport_ret_INCONSISTENT_INPUT_WW;
        }
      }
    } else {
      new_max_pheno_name_blen = catpheno_name_blen;
    }
    const uint32_t new_pheno_ct = old_pheno_ct + 1;
    uintptr_t new_pheno_names_byte_ct = new_pheno_ct * new_max_pheno_name_blen;
    char* pheno_names;
    if (unlikely(pgl_malloc(new_pheno_names_byte_ct, &pheno_names))) {
      goto Plink1ClusterImport_ret_NOMEM;
    }
    if (old_pheno_names && (old_max_pheno_name_blen == new_max_pheno_name_blen)) {
      memcpy(pheno_names, old_pheno_names, old_pheno_ct * new_max_pheno_name_blen);
    } else {
      for (uint32_t pheno_idx = 0; pheno_idx != old_pheno_ct; ++pheno_idx) {
        strcpy(&(pheno_names[pheno_idx * new_max_pheno_name_blen]), &(old_pheno_names[pheno_idx * old_max_pheno_name_blen]));
      }
    }
    memcpy(&(pheno_names[old_pheno_ct * new_max_pheno_name_blen]), catpheno_name, catpheno_name_blen);
    free_cond(old_pheno_names);
    *pheno_names_ptr = pheno_names;

    PhenoCol* new_pheno_cols = S_CAST(PhenoCol*, realloc(*pheno_cols_ptr, new_pheno_ct * sizeof(PhenoCol)));
    if (unlikely(!new_pheno_cols)) {
      goto Plink1ClusterImport_ret_NOMEM;
    }
    *pheno_cols_ptr = new_pheno_cols;
    *pheno_ct_ptr = new_pheno_ct;
    *max_pheno_name_blen_ptr = new_max_pheno_name_blen;
    new_pheno_cols[old_pheno_ct].nonmiss = nullptr;
    new_pheno_cols[old_pheno_ct].type_code = S_CAST(PhenoDtype, kPhenoDtypeCat);

    const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
    uintptr_t* cat_nm = nullptr;
    uint32_t* cat_idxs = nullptr;
    if (!within_fname) {
      if (unlikely(
              bigstack_alloc_w(raw_sample_ctaw, &cat_nm) ||
              bigstack_calloc_u32(raw_sample_ct, &cat_idxs))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      memcpy(cat_nm, sample_include, raw_sample_ctaw * sizeof(intptr_t));
    }
    uint32_t* cat_htable;
    uint32_t cat_htable_size;
    if (unlikely(HtableGoodSizeAlloc(sample_ct + 2, bigstack_left() / 4, &cat_htable, &cat_htable_size))) {
      goto Plink1ClusterImport_ret_NOMEM;
    }
    SetAllU32Arr(cat_htable_size, cat_htable);
    const char* missing_catname = g_missing_catname;
    const uintptr_t data_vec_ct = Int32CtToVecCt(raw_sample_ct);
    const uint32_t missing_catname_slen = strlen(missing_catname);
    const uint32_t missing_catname_hval = Hashceil(missing_catname, missing_catname_slen, cat_htable_size);
    if (within_fname) {
      uintptr_t* already_seen;
      uint32_t* sorted_cat_idxs;
      char* idbuf;
      const char** cur_cat_names;
      if (unlikely(
              bigstack_calloc_w(BitCtToWordCt(sample_ct), &already_seen) ||
              bigstack_calloc_u32(sample_ct, &sorted_cat_idxs) ||
              bigstack_alloc_c(max_sample_id_blen, &idbuf) ||
              bigstack_alloc_kcp(sample_ct + 2, &cur_cat_names))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      cat_htable[missing_catname_hval] = 0;
      cur_cat_names[0] = missing_catname;
      const char na_str[] = "NA";
      uint32_t na_hashval = Hashceil(na_str, 2, cat_htable_size);
      if (na_hashval == missing_catname_hval) {
        if (++na_hashval == cat_htable_size) {
          na_hashval = 0;
        }
      }
      cat_htable[na_hashval] = sample_ct + 1;
      cur_cat_names[sample_ct + 1] = na_str;

      uint32_t* id_map;
      char* sorted_idbox;
      if (unlikely(CopySortStrboxSubset(sample_include, sample_ids, sample_ct, max_sample_id_blen, 1, 0, 0, &sorted_idbox, &id_map))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      reterr = SizeAndInitTextStream(within_fname, bigstack_left() - (bigstack_left() / 4), MAXV(max_thread_ct - 1, 1), &within_txs);
      if (unlikely(reterr)) {
        goto Plink1ClusterImport_ret_TSTREAM_FAIL;
      }
      char* cat_name_write_start = R_CAST(char*, g_bigstack_base);
      char* cat_name_iter = cat_name_write_start;
      char* cat_name_write_max = R_CAST(char*, g_bigstack_end);

      uint32_t nonnull_cat_ct = 0;
      uintptr_t miss_ct = 0;
      uintptr_t duplicate_ct = 0;
      for (char* line_iter = &(TextLineEnd(&within_txs)[-1]); ; line_iter = AdvToDelim(line_iter, '\n')) {
        // need this since we may clobber \n
      Plink1ClusterImport_LINE_ITER_ALREADY_ADVANCED:
        ++line_iter;
        ++line_idx;
        if (!TextGetUnsafe2(&within_txs, &line_iter)) {
          if (likely(!TextStreamErrcode2(&within_txs, &reterr))) {
            break;
          }
          goto Plink1ClusterImport_ret_TSTREAM_FAIL;
        }
        char* fid_start = line_iter;
        char* fid_end = CurTokenEnd(fid_start);
        char* iid_start = FirstNonTspace(fid_end);
        if (unlikely(IsEolnKns(*iid_start))) {
          goto Plink1ClusterImport_ret_MISSING_TOKENS;
        }
        char* iid_end = CurTokenEnd(iid_start);
        const uint32_t fid_slen = fid_end - fid_start;
        const uint32_t iid_slen = iid_end - iid_start;
        const uint32_t id_blen = fid_slen + iid_slen + 2;
        if (id_blen > max_sample_id_blen) {
          ++miss_ct;
          line_iter = iid_end;
          continue;
        }
        char* idbuf_iter = memcpyax(idbuf, fid_start, fid_slen, '\t');
        idbuf_iter = memcpya(idbuf_iter, iid_start, iid_slen);
        *idbuf_iter = '\0';
        uint32_t lb_idx = bsearch_str_lb(idbuf, sorted_idbox, id_blen, max_sample_id_blen, sample_ct);
        *idbuf_iter = ' ';
        const uint32_t ub_idx = bsearch_str_lb(idbuf, sorted_idbox, id_blen, max_sample_id_blen, sample_ct);
        if (ub_idx == lb_idx) {
          ++miss_ct;
          line_iter = iid_end;
          continue;
        }
        char* main_token_start = NextTokenMult(iid_end, mwithin_val);
        if (unlikely(!main_token_start)) {
          goto Plink1ClusterImport_ret_MISSING_TOKENS;
        }
        char* main_token_end = CurTokenEnd(main_token_start);
        line_iter = AdvToDelim(line_iter, '\n');
        *main_token_end = '\0';
        const uint32_t main_token_slen = main_token_end - main_token_start;
        if (unlikely(main_token_slen > kMaxIdSlen)) {
          logerrputs("Error: Category names are limited to " MAX_ID_SLEN_STR " characters.\n");
          goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
        }
        uint32_t hashval = Hashceil(main_token_start, main_token_slen, cat_htable_size);
        const uint32_t main_token_blen = main_token_slen + 1;
        uint32_t cur_htable_entry;
        while (1) {
          cur_htable_entry = cat_htable[hashval];
          if (cur_htable_entry == UINT32_MAX) {
            if (unlikely(main_token_blen > S_CAST(uintptr_t, cat_name_write_max - cat_name_iter))) {
              goto Plink1ClusterImport_ret_NOMEM;
            }
            char* cat_name_start = cat_name_iter;
            cat_name_iter = memcpya(cat_name_iter, main_token_start, main_token_blen);
            cur_cat_names[++nonnull_cat_ct] = cat_name_start;
            cur_htable_entry = nonnull_cat_ct;
            cat_htable[hashval] = cur_htable_entry;
            break;
          }
          if (memequal(main_token_start, cur_cat_names[cur_htable_entry], main_token_blen)) {
            break;
          }
          if (++hashval == cat_htable_size) {
            hashval = 0;
          }
        }
        // permit duplicates if category is identical
        if (IsSet(already_seen, lb_idx)) {
          const uint32_t existing_cat_idx = sorted_cat_idxs[lb_idx];
          if (unlikely(existing_cat_idx != cur_htable_entry)) {
            idbuf[fid_slen] = ' ';
            logpreprintfww("Error: Duplicate sample ID '%s' with conflicting category assignments in --within file.\n", idbuf);
            goto Plink1ClusterImport_ret_MALFORMED_INPUT_2;
          }
          ++duplicate_ct;
        } else {
          SetBit(lb_idx, already_seen);
          for (; lb_idx != ub_idx; ++lb_idx) {
            sorted_cat_idxs[lb_idx] = cur_htable_entry;
          }
        }
        goto Plink1ClusterImport_LINE_ITER_ALREADY_ADVANCED;
      }
      if (unlikely(!nonnull_cat_ct)) {
        logerrputs("Error: All --within categories are null.\n");
        goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
      }
      double dxx;
      const uint32_t prepend_c = (ScanadvDouble(cur_cat_names[1], &dxx) != nullptr);
      if (prepend_c) {
        for (uint32_t catname_idx = 2; catname_idx <= nonnull_cat_ct; ++catname_idx) {
          if (unlikely(!ScanadvDouble(cur_cat_names[catname_idx], &dxx))) {
            logerrputs("Error: Either all non-null --within categories must be numeric, or none can be.\n");
            goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
          }
        }
        logputs("Note: Prepending 'C' to all --within category names.\n");
      } else {
        for (uint32_t catname_idx = 2; catname_idx <= nonnull_cat_ct; ++catname_idx) {
          if (unlikely(ScanadvDouble(cur_cat_names[catname_idx], &dxx))) {
            logerrputs("Error: Either all non-null --within categories must be numeric, or none can be.\n");
            goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
          }
        }
      }
      // see end of e.g. LoadPsam()
      const uintptr_t catname_vec_ct = WordCtToVecCt(nonnull_cat_ct + 1);
      const uintptr_t total_catname_blen = (prepend_c * nonnull_cat_ct) + S_CAST(uintptr_t, cat_name_iter - cat_name_write_start);
      const uintptr_t catname_storage_vec_ct = DivUp(total_catname_blen, kBytesPerVec);
      if (unlikely(vecaligned_malloc((raw_sample_ctaw * kWordsPerVec + data_vec_ct + catname_vec_ct + catname_storage_vec_ct) * kBytesPerVec, &(new_pheno_cols[old_pheno_ct].nonmiss)))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      new_pheno_cols[old_pheno_ct].nonnull_category_ct = nonnull_cat_ct;
      uintptr_t* catdata_iter = new_pheno_cols[old_pheno_ct].nonmiss;
      cat_nm = catdata_iter;
      ZeroWArr(raw_sample_ctaw, cat_nm);
      catdata_iter = &(catdata_iter[raw_sample_ctaw]);

      cat_idxs = R_CAST(uint32_t*, catdata_iter);
      ZeroU32Arr(raw_sample_ct, cat_idxs);
      new_pheno_cols[old_pheno_ct].data.cat = cat_idxs;
      catdata_iter = &(catdata_iter[data_vec_ct * kWordsPerVec]);

      for (uint32_t sorted_sample_idx = 0; sorted_sample_idx != sample_ct; ++sorted_sample_idx) {
        const uint32_t cur_sample_uidx = id_map[sorted_sample_idx];
        uint32_t cur_cat_idx = sorted_cat_idxs[sorted_sample_idx];
        if (cur_cat_idx > sample_ct) {
          cur_cat_idx = 0;
        }
        if (cur_cat_idx) {
          SetBit(cur_sample_uidx, cat_nm);
        }
        cat_idxs[cur_sample_uidx] = cur_cat_idx;
      }

      const char** cur_name_ptrs = R_CAST(const char**, catdata_iter);
      new_pheno_cols[old_pheno_ct].category_names = cur_name_ptrs;
      *cur_name_ptrs++ = missing_catname;
      char* name_storage_iter = R_CAST(char*, &(catdata_iter[catname_vec_ct * kWordsPerVec]));
      cat_name_iter = cat_name_write_start;
      for (uint32_t uii = 0; uii != nonnull_cat_ct; ++uii) {
        char* cur_name_start = name_storage_iter;
        if (prepend_c) {
          *name_storage_iter++ = 'C';
        }
        const uint32_t cur_catname_blen = 1 + strlen(cat_name_iter);
        name_storage_iter = memcpya(name_storage_iter, cat_name_iter, cur_catname_blen);
        *cur_name_ptrs++ = cur_name_start;
        cat_name_iter = &(cat_name_iter[cur_catname_blen]);
      }

      if (duplicate_ct) {
        logprintfww("Note: %" PRIuPTR " duplicate sample ID%s) in --within file.\n", duplicate_ct, (duplicate_ct == 1)? " (with a consistent category assignment" : "s (with consistent category assignments");
      }
      if (miss_ct) {
        snprintf(g_logbuf, kLogbufSize, "--within: %u non-null categories present, %" PRIuPTR " sample ID%s skipped.\n", nonnull_cat_ct, miss_ct, (miss_ct == 1)? "" : "s");
        WordWrapB(0);
      } else {
        snprintf(g_logbuf, kLogbufSize, "--within: %u non-null categories present.\n", nonnull_cat_ct);
      }
      logputsb();
    } else {
      // --family
      cat_htable[missing_catname_hval] = 0xfffffffdU;
      uintptr_t total_catname_blen = 0;  // does not need to include 'NONE'
      uint32_t family_missing_catname_slen = 0;
      if (family_missing_catname) {
        family_missing_catname_slen = strlen(family_missing_catname);
        uint32_t family_missing_catname_hval = Hashceil(family_missing_catname, family_missing_catname_slen, cat_htable_size);
        if (cat_htable[family_missing_catname_hval] == UINT32_MAX) {
          cat_htable[family_missing_catname_hval] = UINT32_MAXM1;
        } else if ((missing_catname_slen != family_missing_catname_slen) || (!memequal(family_missing_catname, missing_catname, missing_catname_slen))) {
          if (++family_missing_catname_hval == cat_htable_size) {
            family_missing_catname_hval = 0;
          }
          cat_htable[family_missing_catname_hval] = UINT32_MAXM1;
        }
      }
      // guaranteed to have enough space
      uint32_t* cat_idx_m1_to_first_sample_uidx = R_CAST(uint32_t*, g_bigstack_base);
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = sample_include[0];
      uint32_t nonnull_cat_ct = 0;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
        const char* cur_fid = &(sample_ids[sample_uidx * max_sample_id_blen]);
        const char* cur_fid_end = AdvToDelim(cur_fid, '\t');
        const uint32_t slen = cur_fid_end - cur_fid;
        const uint32_t blen = slen + 1;
        for (uint32_t hashval = Hashceil(cur_fid, slen, cat_htable_size); ; ) {
          const uint32_t cur_htable_entry = cat_htable[hashval];
          if (cur_htable_entry >= 0xfffffffdU) {
            if (cur_htable_entry == UINT32_MAX) {
              cat_htable[hashval] = sample_uidx;
              total_catname_blen += blen;
              cat_idx_m1_to_first_sample_uidx[nonnull_cat_ct] = sample_uidx;
              cat_idxs[sample_uidx] = ++nonnull_cat_ct;
              break;
            } else if (cur_htable_entry == UINT32_MAXM1) {
              if ((slen == family_missing_catname_slen) && memequal(cur_fid, family_missing_catname, family_missing_catname_slen)) {
                ClearBit(sample_uidx, cat_nm);
                cat_idxs[sample_uidx] = 0;
                break;
              }
            } else {
              if ((slen == missing_catname_slen) && memequal(cur_fid, missing_catname, missing_catname_slen)) {
                ClearBit(sample_uidx, cat_nm);
                cat_idxs[sample_uidx] = 0;
                break;
              }
            }
          } else {
            if (memequal(cur_fid, &(sample_ids[cur_htable_entry * max_sample_id_blen]), blen)) {
              cat_idxs[sample_uidx] = cat_idxs[cur_htable_entry];
              break;
            }
          }
          if (++hashval == cat_htable_size) {
            hashval = 0;
          }
        }
      }
      if (unlikely(!nonnull_cat_ct)) {
        logerrputs("Error: All --family FIDs are null.\n");
        goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
      }
      // add 'P' prefixes?
      double dxx;
      const uint32_t prepend_c = (ScanadvDouble(&(sample_ids[cat_idx_m1_to_first_sample_uidx[0] * max_sample_id_blen]), &dxx) != nullptr);
      if (prepend_c) {
        for (uint32_t uii = 1; uii != nonnull_cat_ct; ++uii) {
          if (unlikely(!ScanadvDouble(&(sample_ids[cat_idx_m1_to_first_sample_uidx[uii] * max_sample_id_blen]), &dxx))) {
            logerrputs("Error: Either all non-null --family FIDs must be numeric, or none can be.\n");
            goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
          }
        }
        logputs("Note: Prepending 'C' to all --family category names.\n");
        total_catname_blen += nonnull_cat_ct;
      } else {
        for (uint32_t uii = 1; uii != nonnull_cat_ct; ++uii) {
          if (unlikely(ScanadvDouble(&(sample_ids[cat_idx_m1_to_first_sample_uidx[uii] * max_sample_id_blen]), &dxx))) {
            logerrputs("Error: Either all non-null --family FIDs must be numeric, or none can be.\n");
            goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
          }
        }
      }
      // see end of e.g. LoadPsam()
      const uintptr_t catname_vec_ct = WordCtToVecCt(nonnull_cat_ct + 1);
      const uintptr_t catname_storage_vec_ct = DivUp(total_catname_blen, kBytesPerVec);
      if (unlikely(vecaligned_malloc((raw_sample_ctaw * kWordsPerVec + data_vec_ct + catname_vec_ct + catname_storage_vec_ct) * kBytesPerVec, &(new_pheno_cols[old_pheno_ct].nonmiss)))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      new_pheno_cols[old_pheno_ct].nonnull_category_ct = nonnull_cat_ct;
      uintptr_t* catdata_iter = new_pheno_cols[old_pheno_ct].nonmiss;
      memcpy(catdata_iter, cat_nm, raw_sample_ctaw * sizeof(intptr_t));
      catdata_iter = &(catdata_iter[raw_sample_ctaw]);

      new_pheno_cols[old_pheno_ct].data.cat = R_CAST(uint32_t*, catdata_iter);
      memcpy(catdata_iter, cat_idxs, data_vec_ct * kBytesPerVec);
      catdata_iter = &(catdata_iter[data_vec_ct * kWordsPerVec]);

      const char** cur_name_ptrs = R_CAST(const char**, catdata_iter);
      new_pheno_cols[old_pheno_ct].category_names = cur_name_ptrs;
      *cur_name_ptrs++ = missing_catname;
      char* name_storage_iter = R_CAST(char*, &(catdata_iter[catname_vec_ct * kWordsPerVec]));
      for (uint32_t uii = 0; uii != nonnull_cat_ct; ++uii) {
        char* cur_name_start = name_storage_iter;
        if (prepend_c) {
          *name_storage_iter++ = 'C';
        }
        const char* cur_fid = &(sample_ids[cat_idx_m1_to_first_sample_uidx[uii] * max_sample_id_blen]);
        const char* cur_fid_end = AdvToDelim(cur_fid, '\t');
        name_storage_iter = memcpyax(name_storage_iter, cur_fid, cur_fid_end - cur_fid, '\0');
        *cur_name_ptrs++ = cur_name_start;
      }
      logprintf("--family: %u non-null categor%s present.\n", nonnull_cat_ct, (nonnull_cat_ct == 1)? "y" : "ies");
    }
  }
  while (0) {
  Plink1ClusterImport_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Plink1ClusterImport_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--within file", &within_txs);
    break;
  Plink1ClusterImport_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --within file has fewer tokens than expected.\n", line_idx);
  Plink1ClusterImport_ret_MALFORMED_INPUT_2:
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  Plink1ClusterImport_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  Plink1ClusterImport_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 Plink1ClusterImport_ret_1:
  CleanupTextStream2("--within file", &within_txs, &reterr);
  BigstackReset(bigstack_mark);
  if (reterr) {
    if (*pheno_names_ptr) {
      free(*pheno_names_ptr);
      *pheno_names_ptr = nullptr;
    }
    CleanupPhenoCols(*pheno_ct_ptr, *pheno_cols_ptr);
    *pheno_ct_ptr = 0;
    *pheno_cols_ptr = nullptr;
  }
  return reterr;
}

PglErr PrescanSampleIds(const char* fname, SampleIdInfo* siip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    reterr = ForceNonFifo(fname);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, "--update-ids file", strerror(errno));
      } else {
        logerrprintfww(kErrprintfRewind, "--update-ids file");
      }
      goto PrescanSampleIds_ret_1;
    }
    reterr = InitTextStream(fname, kTextStreamBlenFast, 1, &txs);
    if (unlikely(reterr)) {
      goto PrescanSampleIds_ret_TSTREAM_FAIL;
    }
    uint32_t is_header_line;
    const char* line_iter;
    do {
      ++line_idx;
      line_iter = TextGet(&txs);
      if (!line_iter) {
        reterr = TextStreamRawErrcode(&txs);
        if (likely(reterr == kPglRetEof)) {
          // permit empty file here, but signal this to caller
          // (reterr == kPglRetSuccess when file is nonempty, so we can detect
          // rewind-fail)
          goto PrescanSampleIds_ret_1;
        }
        goto PrescanSampleIds_ret_TSTREAM_FAIL;
      }
      is_header_line = (*line_iter == '#');
    } while (is_header_line && (!tokequal_k(&(line_iter[1]), "OLD_FID")) && (!tokequal_k(&(line_iter[1]), "OLD_IID")));
    uint32_t old_fid_present = 0;
    uint32_t old_sid_present = 0;
    uint32_t new_fid_present = 0;
    uint32_t new_sid_present = 0;
    if (is_header_line) {
      if (line_iter[5] == 'F') {
        old_fid_present = 1;
        line_iter = FirstNonTspace(&(line_iter[8]));
        if (unlikely(!tokequal_k(line_iter, "OLD_IID"))) {
          logerrputs("Error: Invalid --update-ids file (second header column must be OLD_IID when\nfirst is #OLD_FID).\n");
          goto PrescanSampleIds_ret_MALFORMED_INPUT;
        }
      }
      line_iter = FirstNonTspace(CurTokenEnd(line_iter));
      if (tokequal_k(line_iter, "OLD_SID")) {
        old_sid_present = 1;
        line_iter = FirstNonTspace(&(line_iter[7]));
      }
      if (tokequal_k(line_iter, "NEW_FID")) {
        new_fid_present = 1;
        line_iter = FirstNonTspace(&(line_iter[7]));
      }
      if (unlikely(!tokequal_k(line_iter, "NEW_IID"))) {
        logerrputs("Error: Invalid --update-ids file (no NEW_IID column in expected position).\n");
        goto PrescanSampleIds_ret_MALFORMED_INPUT;
      }
      line_iter = FirstNonTspace(&(line_iter[7]));
      if (tokequal_k(line_iter, "NEW_SID")) {
        new_sid_present = 1;
        line_iter = FirstNonTspace(&(line_iter[7]));
      }
      if (unlikely(!IsSpaceOrEoln(*line_iter))) {
        // Unlike --update-parents/--update-sex, there's no clear benefit to
        // tolerating extra columns, so let's keep this simple.
        logerrputs("Error: Invalid --update-ids file main header line (only permitted columns are\nOLD_FID, OLD_IID, OLD_SID, NEW_FID, NEW_IID, and NEW_SID, in that order).\n");
        goto PrescanSampleIds_ret_MALFORMED_INPUT;
      }
    } else {
      const uint32_t token_ct = CountTokens(line_iter);
      if (token_ct == 4) {
        old_fid_present = 1;
        new_fid_present = 1;
      } else if (unlikely(token_ct != 2)) {
        logerrputs("Error: Invalid --update-ids file (with no #OLD_FID or #OLD_IID header line, 2\nor 4 columns expected).\n");
        goto PrescanSampleIds_ret_MALFORMED_INPUT;
      }
    }
    const uint32_t initial_skip_ct = 1 + old_fid_present + old_sid_present;
    uintptr_t max_sample_id_blen = 0;
    uintptr_t max_sid_blen = 0;
    if (is_header_line) {
      line_iter = AdvPastDelim(line_iter, '\n');
      ++line_idx;
    }
    for (; TextGetUnsafe2K(&txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      line_iter = NextTokenMult(line_iter, initial_skip_ct);
      if (unlikely(!line_iter)) {
        goto PrescanSampleIds_ret_MISSING_TOKENS;
      }
      const char* id_start = line_iter;
      line_iter = CurTokenEnd(id_start);
      uintptr_t cur_sample_id_blen = 1 + S_CAST(uintptr_t, line_iter - id_start);
      if (new_fid_present) {
        const char* iid_start = FirstNonTspace(line_iter);
        if (unlikely(IsEolnKns(*iid_start))) {
          goto PrescanSampleIds_ret_MISSING_TOKENS;
        }
        line_iter = CurTokenEnd(iid_start);
        cur_sample_id_blen += 1 + S_CAST(uintptr_t, line_iter - iid_start);
      }
      if (cur_sample_id_blen > max_sample_id_blen) {
        max_sample_id_blen = cur_sample_id_blen;
      }
      if (new_sid_present) {
        const char* sid_start = FirstNonTspace(line_iter);
        if (unlikely(IsEolnKns(*sid_start))) {
          goto PrescanSampleIds_ret_MISSING_TOKENS;
        }
        line_iter = CurTokenEnd(sid_start);
        const uintptr_t sid_slen = line_iter - sid_start;
        if (sid_slen >= max_sid_blen) {
          max_sid_blen = sid_slen + 1;
        }
      }
      if (unlikely(!IsEolnKns(*line_iter))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --update-ids file has more tokens than expected.\n", line_idx);
        goto PrescanSampleIds_ret_MALFORMED_INPUT_WW;
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto PrescanSampleIds_ret_TSTREAM_FAIL;
    }
    reterr = kPglRetSuccess;
    siip->max_sample_id_blen = max_sample_id_blen;
    if (new_sid_present) {
      siip->max_sid_blen = max_sid_blen;
    }
  }
  while (0) {
  PrescanSampleIds_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--update-ids file", &txs);
    break;
  PrescanSampleIds_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --update-ids file has fewer tokens than expected.\n", line_idx);
  PrescanSampleIds_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  PrescanSampleIds_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 PrescanSampleIds_ret_1:
  CleanupTextStream2("--update-ids file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr PrescanParentalIds(const char* fname, uint32_t max_thread_ct, ParentalIdInfo* parental_id_infop) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    reterr = ForceNonFifo(fname);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, "--update-parents file", strerror(errno));
      } else {
        logerrprintfww(kErrprintfRewind, "--update-parents file");
      }
      goto PrescanParentalIds_ret_1;
    }
    // permit very long lines since this can be pointed at .ped files
    // possible minor todo: could save longest line length for later reference
    reterr = SizeAndInitTextStream(fname, bigstack_left() - (bigstack_left() / 4), MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto PrescanParentalIds_ret_TSTREAM_FAIL;
    }

    uint32_t is_header_line;
    const char* line_iter;
    do {
      ++line_idx;
      line_iter = TextGet(&txs);
      if (!line_iter) {
        reterr = TextStreamRawErrcode(&txs);
        if (likely(reterr == kPglRetEof)) {
          // permit empty file here, but signal this to caller
          // (reterr == kPglRetSuccess when file is nonempty, so we can detect
          // rewind-fail)
          goto PrescanParentalIds_ret_1;
        }
        goto PrescanParentalIds_ret_TSTREAM_FAIL;
      }
      is_header_line = (*line_iter == '#');
    } while (is_header_line && (!tokequal_k(&(line_iter[1]), "FID")) && (!tokequal_k(&(line_iter[1]), "IID")));
    uintptr_t max_paternal_id_blen = 0;
    uintptr_t max_maternal_id_blen = 0;
    uint32_t pat_col_idx;
    if (is_header_line) {
      // Search for 'PAT' column.  Require all-caps.
      pat_col_idx = 0;
      line_iter = &(line_iter[4]);
      while (1) {
        ++pat_col_idx;
        const char* token_start = FirstNonTspace(line_iter);
        if (unlikely(IsEolnKns(*token_start))) {
          logerrputs("Error: Invalid --update-parents file (no PAT column).\n");
          goto PrescanParentalIds_ret_MALFORMED_INPUT;
        }
        line_iter = CurTokenEnd(token_start);
        if (strequal_k(token_start, "PAT", line_iter - token_start)) {
          break;
        }
      }
      // Require immediately-following column to be 'MAT'.
      const char* token_start = FirstNonTspace(line_iter);
      if (unlikely(IsEolnKns(*token_start))) {
        logerrputs("Error: Invalid --update-parents file (no MAT column immediately following PAT\ncolumn).\n");
        goto PrescanParentalIds_ret_MALFORMED_INPUT;
      }
      line_iter = CurTokenEnd(token_start);
      if (unlikely(!strequal_k(token_start, "MAT", line_iter - token_start))) {
        logerrputs("Error: Invalid --update-parents file (no MAT column immediately following PAT\ncolumn).\n");
        goto PrescanParentalIds_ret_MALFORMED_INPUT;
      }
      line_iter = AdvPastDelim(line_iter, '\n');
      ++line_idx;
    } else {
      const uint32_t token_ct = CountTokens(line_iter);
      if (unlikely(token_ct < 3)) {
        logerrputs("Error: Invalid --update-parents file (3+ columns expected).\n");
        goto PrescanParentalIds_ret_MALFORMED_INPUT;
      }
      pat_col_idx = (token_ct == 3)? 1 : 2;
    }
    for (; TextGetUnsafe2K(&txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      line_iter = NextTokenMult(line_iter, pat_col_idx);
      if (unlikely(!line_iter)) {
        goto PrescanParentalIds_ret_MISSING_TOKENS;
      }
      const char* pat_start = line_iter;
      line_iter = CurTokenEnd(pat_start);
      const uintptr_t cur_paternal_id_slen = line_iter - pat_start;
      if (cur_paternal_id_slen >= max_paternal_id_blen) {
        max_paternal_id_blen = cur_paternal_id_slen + 1;
      }
      const char* mat_start = FirstNonTspace(line_iter);
      if (unlikely(IsEolnKns(*mat_start))) {
        goto PrescanParentalIds_ret_MISSING_TOKENS;
      }
      line_iter = CurTokenEnd(mat_start);
      uintptr_t cur_maternal_id_slen = line_iter - mat_start;
      if (cur_maternal_id_slen >= max_maternal_id_blen) {
        max_maternal_id_blen = cur_maternal_id_slen + 1;
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto PrescanParentalIds_ret_TSTREAM_FAIL;
    }
    reterr = kPglRetSuccess;
    parental_id_infop->max_paternal_id_blen = max_paternal_id_blen;
    parental_id_infop->max_maternal_id_blen = max_maternal_id_blen;
  }
  while (0) {
  PrescanParentalIds_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--update-parents file", &txs);
    break;
  PrescanParentalIds_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --update-parents file has fewer tokens than expected.\n", line_idx);
    WordWrapB(0);
    logerrputsb();
  PrescanParentalIds_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 PrescanParentalIds_ret_1:
  CleanupTextStream2("--update-parents file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr UpdateSampleIds(const char* fname, const uintptr_t* sample_include, uint32_t raw_sample_ct, uintptr_t sample_ct, SampleIdInfo* siip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    if (!sample_ct) {
      goto UpdateSampleIds_ret_1;
    }
    // probable todo: deduplicate shared code with PrescanSampleIds
    reterr = InitTextStream(fname, kTextStreamBlenFast, 1, &txs);
    if (unlikely(reterr)) {
      goto UpdateSampleIds_ret_TSTREAM_FAIL;
    }
    const char* line_iter;
    uint32_t is_header_line;
    do {
      ++line_idx;
      line_iter = TextGet(&txs);
      if (unlikely(!line_iter)) {
        reterr = TextStreamRawErrcode(&txs);
        // This function is no longer called when the original file is empty.
        goto UpdateSampleIds_ret_TSTREAM_FAIL;
      }
      is_header_line = (*line_iter == '#');
    } while (is_header_line && (!tokequal_k(&(line_iter[1]), "OLD_FID")) && (!tokequal_k(&(line_iter[1]), "OLD_IID")));
    uint32_t old_fid_present = 0;
    uint32_t old_sid_present = 0;
    uint32_t new_fid_present = 0;
    uint32_t new_sid_present = 0;
    if (is_header_line) {
      if (line_iter[5] == 'F') {
        old_fid_present = 1;
        line_iter = FirstNonTspace(&(line_iter[8]));
        if (unlikely(!tokequal_k(line_iter, "OLD_IID"))) {
          goto UpdateSampleIds_ret_REWIND_FAIL;
        }
      }
      line_iter = FirstNonTspace(CurTokenEnd(line_iter));
      if (tokequal_k(line_iter, "OLD_SID")) {
        old_sid_present = 1;
        line_iter = FirstNonTspace(&(line_iter[7]));
      }
      if (tokequal_k(line_iter, "NEW_FID")) {
        new_fid_present = 1;
        line_iter = FirstNonTspace(&(line_iter[7]));
      }
      if (unlikely(!tokequal_k(line_iter, "NEW_IID"))) {
        goto UpdateSampleIds_ret_REWIND_FAIL;
      }
      line_iter = FirstNonTspace(&(line_iter[7]));
      if (tokequal_k(line_iter, "NEW_SID")) {
        new_sid_present = 1;
        line_iter = FirstNonTspace(&(line_iter[7]));
      }
      if (unlikely(!IsSpaceOrEoln(*line_iter))) {
        goto UpdateSampleIds_ret_REWIND_FAIL;
      }
    } else {
      const uint32_t token_ct = CountTokens(line_iter);
      if (token_ct == 4) {
        old_fid_present = 1;
        new_fid_present = 1;
      } else if (unlikely(token_ct != 2)) {
        goto UpdateSampleIds_ret_REWIND_FAIL;
      }
    }
    XidMode xid_mode;
    if (old_fid_present) {
      xid_mode = old_sid_present? kfXidModeFidIidSid : kfXidModeFidIid;
    } else {
      xid_mode = old_sid_present? kfXidModeIidSid : kfXidModeIid;
    }
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    uint32_t* xid_map = nullptr;
    char* sorted_xidbox = nullptr;
    const uint32_t allow_dups = siip->sids && (!old_sid_present);
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, allow_dups, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto UpdateSampleIds_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* already_seen;
    char* idbuf;
    if (unlikely(
            bigstack_calloc_w(raw_sample_ctl, &already_seen) ||
            bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto UpdateSampleIds_ret_NOMEM;
    }
    uint32_t hit_ct = 0;
    uintptr_t miss_ct = 0;
    if (is_header_line) {
      line_iter = AdvPastDelim(line_iter, '\n');
      ++line_idx;
    }
    for (; TextGetUnsafe2K(&txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      const char* linebuf_iter = line_iter;
      uint32_t xid_idx;
      uint32_t xid_idx_end;
      if (SortedXidboxReadMultifind(sorted_xidbox, max_xid_blen, sample_ct, 0, xid_mode, &linebuf_iter, &xid_idx, &xid_idx_end, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto UpdateSampleIds_ret_REWIND_FAIL;
        }
        ++miss_ct;
        continue;
      }
      uint32_t sample_uidx = xid_map[xid_idx];
      if (unlikely(IsSet(already_seen, sample_uidx))) {
        // possible todo: tolerate duplicates if payload matches, like
        // --update-sex does
        logerrprintfww("Error: Sample ID on line %" PRIuPTR " of --update-ids file duplicates one earlier in the file.\n", line_idx);
        goto UpdateSampleIds_ret_MALFORMED_INPUT;
      }
      SetBit(sample_uidx, already_seen);
      char* new_sample_id = &(siip->sample_ids[sample_uidx * max_sample_id_blen]);
      char* new_sample_id_iter = new_sample_id;
      if (new_fid_present) {
        const char* sample_id_start = FirstNonTspace(linebuf_iter);
        linebuf_iter = CurTokenEnd(sample_id_start);
        new_sample_id_iter = memcpya(new_sample_id_iter, sample_id_start, linebuf_iter - sample_id_start);
      } else {
        *new_sample_id_iter++ = '0';
      }
      *new_sample_id_iter++ = '\t';
      const char* iid_start = FirstNonTspace(linebuf_iter);
      linebuf_iter = CurTokenEnd(iid_start);
      new_sample_id_iter = memcpyax(new_sample_id_iter, iid_start, linebuf_iter - iid_start, '\0');
      char* new_sid = nullptr;
      uint32_t new_sid_blen = 0;
      if (new_sid_present) {
        new_sid = &(siip->sids[sample_uidx * max_sid_blen]);
        const char* sid_start = FirstNonTspace(linebuf_iter);
        linebuf_iter = CurTokenEnd(sid_start);
        const uint32_t new_sid_slen = linebuf_iter - sid_start;
        memcpyx(new_sid, sid_start, new_sid_slen, '\0');
        new_sid_blen = new_sid_slen + 1;
      }
      line_iter = linebuf_iter;
      hit_ct += xid_idx_end - xid_idx;
      if (++xid_idx != xid_idx_end) {
        uintptr_t new_sample_id_blen = new_sample_id_iter - new_sample_id;
        do {
          sample_uidx = xid_map[xid_idx];
          memcpy(&(siip->sample_ids[sample_uidx * max_sample_id_blen]), new_sample_id, new_sample_id_blen);
          if (new_sid) {
            // now possible since NEW_SID may be present without OLD_SID
            memcpy(&(siip->sids[sample_uidx * max_sid_blen]), new_sid, new_sid_blen);
          }
        } while (++xid_idx != xid_idx_end);
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto UpdateSampleIds_ret_TSTREAM_FAIL;
    }
    reterr = kPglRetSuccess;
    if (miss_ct) {
      snprintf(g_logbuf, kLogbufSize, "--update-ids: %u sample%s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      snprintf(g_logbuf, kLogbufSize, "--update-ids: %u sample%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logputsb();
  }
  while (0) {
  UpdateSampleIds_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  UpdateSampleIds_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--update-ids file", &txs);
    break;
  UpdateSampleIds_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, "--update-ids file");
    reterr = kPglRetRewindFail;
    break;
  UpdateSampleIds_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 UpdateSampleIds_ret_1:
  CleanupTextStream2("--update-ids file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr UpdateSampleParents(const char* fname, const SampleIdInfo* siip, const uintptr_t* sample_include, uint32_t raw_sample_ct, uintptr_t sample_ct, uint32_t max_thread_ct, ParentalIdInfo* parental_id_infop, uintptr_t* founder_info) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    if (!sample_ct) {
      goto UpdateSampleParents_ret_1;
    }
    // This repeats code in PrescanParentalIds(); probably want to pull this
    // into its own function.
    reterr = SizeAndInitTextStream(fname, bigstack_left() - (bigstack_left() / 4), MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto UpdateSampleParents_ret_TSTREAM_FAIL;
    }
    const char* line_iter;
    uint32_t is_header_line;
    do {
      ++line_idx;
      line_iter = TextGet(&txs);
      if (unlikely(!line_iter)) {
        reterr = TextStreamRawErrcode(&txs);
        goto UpdateSampleParents_ret_TSTREAM_FAIL;
      }
      is_header_line = (*line_iter == '#');
    } while (is_header_line && (!tokequal_k(&(line_iter[1]), "FID")) && (!tokequal_k(&(line_iter[1]), "IID")));
    uint32_t postid_pat_col_idx = 1;
    XidMode xid_mode;
    if (is_header_line) {
      if (line_iter[1] == 'I') {
        xid_mode = kfXidModeFlagNeverFid;
        line_iter = &(line_iter[4]);
      } else {
        xid_mode = kfXidModeFidIid;
        line_iter = FirstNonTspace(&(line_iter[4]));
        if (unlikely(!tokequal_k(line_iter, "IID"))) {
          logerrprintf("Error: No IID column on line %" PRIuPTR " of --update-parents file.\n", line_idx);
          goto UpdateSampleParents_ret_MALFORMED_INPUT;
        }
        line_iter = &(line_iter[3]);
      }
      line_iter = FirstNonTspace(line_iter);
      if (tokequal_k(line_iter, "SID")) {
        xid_mode |= kfXidModeFlagSid;
        line_iter = &(line_iter[3]);
      }
      // Search for 'PAT' column.  Require all-caps.
      while (1) {
        const char* token_start = FirstNonTspace(line_iter);
        if (unlikely(IsEolnKns(*token_start))) {
          // previously validated
          goto UpdateSampleParents_ret_REWIND_FAIL;
        }
        line_iter = CurTokenEnd(token_start);
        if (strequal_k(token_start, "PAT", line_iter - token_start)) {
          break;
        }
        ++postid_pat_col_idx;
      }
    } else {
      const uint32_t token_ct = CountTokens(line_iter);
      if (unlikely(token_ct < 3)) {
        goto UpdateSampleParents_ret_REWIND_FAIL;
      }
      xid_mode = (token_ct == 3)? kfXidModeIid : kfXidModeFidIid;
    }

    uint32_t* xid_map = nullptr;
    char* sorted_xidbox = nullptr;
    const uint32_t allow_dups = siip->sids && (!(xid_mode & kfXidModeFlagSid));
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, allow_dups, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto UpdateSampleParents_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* already_seen;
    char* idbuf;
    if (unlikely(
            bigstack_calloc_w(raw_sample_ctl, &already_seen) ||
            bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto UpdateSampleParents_ret_NOMEM;
    }
    const uintptr_t max_paternal_id_blen = parental_id_infop->max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = parental_id_infop->max_maternal_id_blen;
    char* paternal_ids = parental_id_infop->paternal_ids;
    char* maternal_ids = parental_id_infop->maternal_ids;
    uint32_t hit_ct = 0;
    uintptr_t miss_ct = 0;
    if (is_header_line) {
      line_iter = AdvPastDelim(line_iter, '\n');
      ++line_idx;
    }
    for (; TextGetUnsafe2K(&txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      const char* linebuf_iter = line_iter;
      uint32_t xid_idx;
      uint32_t xid_idx_end;
      if (SortedXidboxReadMultifind(sorted_xidbox, max_xid_blen, sample_ct, 0, xid_mode, &linebuf_iter, &xid_idx, &xid_idx_end, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto UpdateSampleParents_ret_REWIND_FAIL;
        }
        ++miss_ct;
        continue;
      }
      uint32_t sample_uidx = xid_map[xid_idx];
      if (unlikely(IsSet(already_seen, sample_uidx))) {
        // possible todo: tolerate duplicates if payload matches, like
        // --update-sex does
        logerrprintfww("Error: Sample ID on line %" PRIuPTR " of --update-parents file duplicates one earlier in the file.\n", line_idx);
        goto UpdateSampleParents_ret_MALFORMED_INPUT;
      }
      SetBit(sample_uidx, already_seen);
      const char* pat_start = NextTokenMult(linebuf_iter, postid_pat_col_idx);
      const char* pat_end = CurTokenEnd(pat_start);
      const char* mat_start = FirstNonTspace(pat_end);
      line_iter = CurTokenEnd(mat_start);
      const uint32_t plen = pat_end - pat_start;
      const uint32_t mlen = line_iter - mat_start;
      const uint32_t is_founder = (plen == 1) && (*pat_start == '0') && (mlen == 1) && (*mat_start == '0');
      hit_ct += xid_idx_end - xid_idx;
      while (1) {
        memcpyx(&(paternal_ids[sample_uidx * max_paternal_id_blen]), pat_start, plen, '\0');
        memcpyx(&(maternal_ids[sample_uidx * max_maternal_id_blen]), mat_start, mlen, '\0');
        AssignBit(sample_uidx, is_founder, founder_info);
        if (++xid_idx == xid_idx_end) {
          break;
        }
        sample_uidx = xid_map[xid_idx];
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto UpdateSampleParents_ret_REWIND_FAIL;
    }
    reterr = kPglRetSuccess;
    if (miss_ct) {
      snprintf(g_logbuf, kLogbufSize, "--update-parents: %u sample%s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      snprintf(g_logbuf, kLogbufSize, "--update-parents: %u sample%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logputsb();
  }
  while (0) {
  UpdateSampleParents_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  UpdateSampleParents_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--update-parents file", &txs);
    break;
  UpdateSampleParents_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, "--update-parents file");
    reterr = kPglRetRewindFail;
    break;
  UpdateSampleParents_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 UpdateSampleParents_ret_1:
  CleanupTextStream2("--update-parents file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr UpdateSampleSexes(const uintptr_t* sample_include, const SampleIdInfo* siip, const UpdateSexInfo* update_sex_info_ptr, uint32_t raw_sample_ct, uintptr_t sample_ct, uint32_t max_thread_ct, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    if (!sample_ct) {
      goto UpdateSampleSexes_ret_1;
    }
    // permit very long lines since this can be pointed at .ped files
    reterr = SizeAndInitTextStream(update_sex_info_ptr->fname, bigstack_left() - (bigstack_left() / 4), MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto UpdateSampleSexes_ret_1;
    }

    // (Much of this boilerplate is shared with e.g. KeepFcol(); it probably
    // belongs in its own function.)
    char* line_start;
    char* header_sample_id_end;
    XidMode xid_mode;
    reterr = LoadXidHeader("update-sex", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeaderFixedWidth : kfXidHeaderFixedWidthIgnoreSid, &line_idx, &txs, &xid_mode, &line_start, &header_sample_id_end);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: Empty --update-sex file.\n");
        goto UpdateSampleSexes_ret_MALFORMED_INPUT;
      }
      goto UpdateSampleSexes_ret_TSTREAM_XID_FAIL;
    }
    const uint32_t id_col_ct = GetXidColCt(xid_mode);
    uint32_t col_num = update_sex_info_ptr->col_num;
    uint32_t postid_col_idx = 0;
    if ((*line_start == '#') && (!col_num)) {
      // search for 'SEX' column (any capitalization)
      const char* token_end = header_sample_id_end;
      while (1) {
        ++postid_col_idx;
        const char* linebuf_iter = FirstNonTspace(token_end);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          logerrputs("Error: No 'SEX' column in --update-sex file, and no column number specified.\n");
          goto UpdateSampleSexes_ret_MALFORMED_INPUT;
        }
        token_end = CurTokenEnd(linebuf_iter);
        if (MatchUpperKLen(linebuf_iter, "SEX", token_end - linebuf_iter)) {
          break;
        }
      }
    } else {
      if (!col_num) {
        if (unlikely(id_col_ct == 3)) {
          logerrputs("Error: You must use 'col-num=' to specify the position of the sex column in the\n--update-sex file.\n");
          goto UpdateSampleSexes_ret_MALFORMED_INPUT;
        }
        col_num = 3;
      }
      if (unlikely(id_col_ct >= col_num)) {
        logerrputs("Error: --update-sex 'col-num=' argument too small (it refers to a sample ID\ncolumn).\n");
        goto UpdateSampleSexes_ret_MALFORMED_INPUT;
      }
      postid_col_idx = col_num - id_col_ct;
    }

    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uint32_t* xid_map = nullptr;
    char* sorted_xidbox = nullptr;
    const uint32_t allow_dups = siip->sids && (!(xid_mode & kfXidModeFlagSid));
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, allow_dups, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto UpdateSampleSexes_ret_1;
    }
    uintptr_t* already_seen;
    char* idbuf;
    if (unlikely(
            bigstack_calloc_w(raw_sample_ctl, &already_seen) ||
            bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto UpdateSampleSexes_ret_NOMEM;
    }

    const uint32_t male0 = (update_sex_info_ptr->flags / kfUpdateSexMale0) & 1;
    uint32_t hit_ct = 0;
    uintptr_t miss_ct = 0;
    uintptr_t duplicate_ct = 0;
    if (*line_start == '#') {
      ++line_idx;
      line_start = TextGet(&txs);
    }
    for (; line_start; ++line_idx, line_start = TextGet(&txs)) {
      const char* linebuf_iter = line_start;
      uint32_t xid_idx_start;
      uint32_t xid_idx_end;
      if (SortedXidboxReadMultifind(sorted_xidbox, max_xid_blen, sample_ct, 0, xid_mode, &linebuf_iter, &xid_idx_start, &xid_idx_end, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto UpdateSampleSexes_ret_MISSING_TOKENS;
        }
        ++miss_ct;
        continue;
      }
      const char* sex_start = NextTokenMult(linebuf_iter, postid_col_idx);
      if (unlikely(!sex_start)) {
        goto UpdateSampleSexes_ret_MISSING_TOKENS;
      }
      uint32_t sexval = ctou32(*sex_start);
      const uint32_t ujj = sexval & 0xdfU;
      sexval -= 48;
      if (sexval > 2) {
        if (ujj == 77) {
          // 'M'/'m'
          sexval = 1;
        } else if (ujj == 70) {
          // 'F'/'f'
          sexval = 2;
        } else if (unlikely((!male0) && (sexval != 30))) {
          // allow 'N' = missing to make 1/2/NA work
          // don't permit 'n' for now
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid sex value on line %" PRIuPTR " of --update-sex file. (Acceptable values: 1/M/m = male, 2/F/f = female, 0/N = missing.)\n", line_idx);
          goto UpdateSampleSexes_ret_MALFORMED_INPUT_WW;
        } else {
          // with 'male0', everything else is treated as missing
          sexval = 0;
        }
      } else if (male0) {
        if (unlikely(sexval == 2)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid sex value on line %" PRIuPTR " of --update-sex file. ('2' is prohibited when the 'male0' modifier is present.)\n", line_idx);
          goto UpdateSampleSexes_ret_MALFORMED_INPUT_WW;
        }
        ++sexval;
      }
      uint32_t sample_uidx = xid_map[xid_idx_start];
      if (IsSet(already_seen, sample_uidx)) {
        // permit duplicates iff sex value is identical
        const uint32_t old_sexval = IsSet(sex_nm, sample_uidx) * (2 - IsSet(sex_male, sample_uidx));
        if (unlikely(sexval != old_sexval)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Sample ID on line %" PRIuPTR " of --update-sex file duplicates one earlier in the file, and sex values don't match.\n", line_idx);
          goto UpdateSampleSexes_ret_MALFORMED_INPUT_WW;
        }
        ++duplicate_ct;
        continue;
      }
      SetBit(sample_uidx, already_seen);
      hit_ct += xid_idx_end - xid_idx_start;
      for (uint32_t xid_idx = xid_idx_start; ; sample_uidx = xid_map[xid_idx]) {
        if (sexval) {
          SetBit(sample_uidx, sex_nm);
          if (sexval == 1) {
            SetBit(sample_uidx, sex_male);
          } else {
            ClearBit(sample_uidx, sex_male);
          }
        } else {
          ClearBit(sample_uidx, sex_nm);
          ClearBit(sample_uidx, sex_male);
        }
        if (++xid_idx == xid_idx_end) {
          break;
        }
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto UpdateSampleSexes_ret_TSTREAM_FAIL;
    }
    if (duplicate_ct) {
      logprintfww("Note: %" PRIuPTR " duplicate sample ID%s) in --update-sex file.\n", duplicate_ct, (duplicate_ct == 1)? " (with a consistent sex assignment" : "s (with consistent sex assignments");
    }
    if (miss_ct) {
      snprintf(g_logbuf, kLogbufSize, "--update-sex: %u sample%s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      snprintf(g_logbuf, kLogbufSize, "--update-sex: %u sample%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logputsb();
  }
  while (0) {
  UpdateSampleSexes_ret_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&txs)) {
      break;
    }
  UpdateSampleSexes_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--update-sex file", &txs);
    break;
  UpdateSampleSexes_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  UpdateSampleSexes_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  UpdateSampleSexes_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  UpdateSampleSexes_ret_MISSING_TOKENS:
    logerrprintf("Error: Line %" PRIuPTR " of --update-sex file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  }
 UpdateSampleSexes_ret_1:
  CleanupTextStream2("--update-sex file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr SplitCatPheno(const char* split_cat_phenonames_flattened, const uintptr_t* sample_include, uint32_t raw_sample_ct, PhenoTransformFlags pheno_transform_flags, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr, PhenoCol** covar_cols_ptr, char** covar_names_ptr, uint32_t* covar_ct_ptr, uintptr_t* max_covar_name_blen_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  const char* doomed_pheno_names = nullptr;
  PhenoCol* doomed_pheno_cols = nullptr;
  uint32_t doomed_pheno_ct = 0;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t is_omit = ((pheno_transform_flags & (kfPhenoTransformSplitCatOmitMost | kfPhenoTransformSplitCatOmitLast)) != 0);
    const uint32_t omit_most = (pheno_transform_flags / kfPhenoTransformSplitCatOmitMost) & 1;
    uint32_t qt_12 = 0;
    uint32_t at_least_one_cat_pheno_processed = 0;
    for (uint32_t is_covar = 0; is_covar != 2; ++is_covar) {
      PhenoCol** xpheno_cols_ptr;
      char** xpheno_names_ptr;
      uint32_t* xpheno_ct_ptr;
      uintptr_t* max_xpheno_name_blen_ptr;
      if (!is_covar) {
        xpheno_cols_ptr = pheno_cols_ptr;
        xpheno_names_ptr = pheno_names_ptr;
        xpheno_ct_ptr = pheno_ct_ptr;
        max_xpheno_name_blen_ptr = max_pheno_name_blen_ptr;
      } else {
        if (!split_cat_phenonames_flattened) {
          break;
        }
        BigstackReset(bigstack_mark);
        xpheno_cols_ptr = covar_cols_ptr;
        xpheno_names_ptr = covar_names_ptr;
        xpheno_ct_ptr = covar_ct_ptr;
        max_xpheno_name_blen_ptr = max_covar_name_blen_ptr;
        qt_12 = !(pheno_transform_flags & kfPhenoTransformSplitCatCovar01);
      }
      const uint32_t old_pheno_ct = *xpheno_ct_ptr;
      if (!old_pheno_ct) {
        continue;
      }
      const uint32_t old_pheno_ctl = BitCtToWordCt(old_pheno_ct);
      const uintptr_t old_max_pheno_name_blen = *max_xpheno_name_blen_ptr;
      PhenoCol* old_pheno_cols = *xpheno_cols_ptr;
      const char* old_pheno_names = *xpheno_names_ptr;
      uintptr_t* phenos_to_split;
      if (unlikely(bigstack_calloc_w(old_pheno_ctl, &phenos_to_split))) {
        goto SplitCatPheno_ret_NOMEM;
      }
      if (!split_cat_phenonames_flattened) {
        for (uint32_t pheno_idx = 0; pheno_idx != old_pheno_ct; ++pheno_idx) {
          const PhenoCol* cur_pheno_col = &(old_pheno_cols[pheno_idx]);
          if (cur_pheno_col->type_code == kPhenoDtypeCat) {
            if (unlikely(strchr(&(old_pheno_names[pheno_idx * old_max_pheno_name_blen]), '='))) {
              logerrputs("Error: --split-cat-pheno cannot be used on phenotypes containing the '='\ncharacter.\n");
              goto SplitCatPheno_ret_INCONSISTENT_INPUT;
            }
            SetBit(pheno_idx, phenos_to_split);
          }
        }
      } else {
        uint32_t* id_htable;
        uint32_t id_htable_size;
        if (unlikely(HtableGoodSizeAlloc(old_pheno_ct, bigstack_left(), &id_htable, &id_htable_size))) {
          goto SplitCatPheno_ret_NOMEM;
        }
        // shouldn't be possible for this to fail
        PopulateStrboxHtable(old_pheno_names, old_pheno_ct, old_max_pheno_name_blen, id_htable_size, id_htable);
        const char* split_cat_phenonames_iter = split_cat_phenonames_flattened;
        do {
          const uint32_t cur_phenoname_slen = strlen(split_cat_phenonames_iter);
          if (cur_phenoname_slen < old_max_pheno_name_blen) {
            uint32_t pheno_idx = StrboxHtableFind(split_cat_phenonames_iter, old_pheno_names, id_htable, old_max_pheno_name_blen, cur_phenoname_slen, id_htable_size);
            if (pheno_idx != UINT32_MAX) {
              if (unlikely(old_pheno_cols[pheno_idx].type_code != kPhenoDtypeCat)) {
                snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a categorical %s.\n", split_cat_phenonames_iter, is_covar? "covariate" : "phenotype");
                goto SplitCatPheno_ret_INCONSISTENT_INPUT_WW;
              }
              SetBit(pheno_idx, phenos_to_split);
            }
          }
          split_cat_phenonames_iter = &(split_cat_phenonames_iter[cur_phenoname_slen + 1]);
        } while (*split_cat_phenonames_iter);
        BigstackReset(id_htable);
      }
      const uint32_t split_pheno_ct = PopcountWords(phenos_to_split, old_pheno_ctl);
      if (!split_pheno_ct) {
        continue;
      }
      at_least_one_cat_pheno_processed = 1;
      // first pass: determine new memory allocation sizes
      const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
      uintptr_t new_max_pheno_name_blen = old_max_pheno_name_blen;

      // excludes null; also excludes one other category if omit-most/omit-last
      uint32_t* observed_cat_cts;

      uintptr_t** observed_cats;
      uintptr_t* sample_include_intersect;
      if (unlikely(
              bigstack_alloc_u32(split_pheno_ct, &observed_cat_cts) ||
              bigstack_alloc_wp(split_pheno_ct, &observed_cats) ||
              bigstack_alloc_w(raw_sample_ctl, &sample_include_intersect))) {
        goto SplitCatPheno_ret_NOMEM;
      }
      uint32_t* omitted_cat_uidxs = nullptr;
      if (is_omit) {
        if (unlikely(bigstack_alloc_u32(split_pheno_ct, &omitted_cat_uidxs))) {
          goto SplitCatPheno_ret_NOMEM;
        }
      }
      uint32_t* cat_obs_buf = nullptr;
      uintptr_t create_pheno_ct = 0;
      uintptr_t split_pheno_uidx_base = 0;
      uintptr_t phenos_to_split_bits = phenos_to_split[0];
      uint32_t max_cat_uidx = 0;
      for (uint32_t split_pheno_idx = 0; split_pheno_idx != split_pheno_ct; ++split_pheno_idx) {
        const uintptr_t split_pheno_uidx = BitIter1(phenos_to_split, &split_pheno_uidx_base, &phenos_to_split_bits);
        const PhenoCol* cur_pheno_col = &(old_pheno_cols[split_pheno_uidx]);
        BitvecAndCopy(sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, sample_include_intersect);
        const uint32_t cur_cat_ct = cur_pheno_col->nonnull_category_ct + 1;
        const uint32_t cur_cat_ctl = BitCtToWordCt(cur_cat_ct);
        uintptr_t* cur_observed_cats;
        if (unlikely(bigstack_alloc_w(cur_cat_ctl, &cur_observed_cats))) {
          goto SplitCatPheno_ret_NOMEM;
        }
        if (omit_most) {
          if (unlikely(bigstack_alloc_u32(cur_cat_ct, &cat_obs_buf))) {
            goto SplitCatPheno_ret_NOMEM;
          }
        }
        observed_cats[split_pheno_idx] = cur_observed_cats;
        const uint32_t cur_nmiss_ct = PopcountWords(sample_include_intersect, raw_sample_ctl);
        uint32_t cur_observed_cat_ct;
        if (!omit_most) {
          cur_observed_cat_ct = IdentifyRemainingCats(sample_include_intersect, cur_pheno_col, cur_nmiss_ct, cur_observed_cats);
          if (is_omit && cur_observed_cat_ct) {
            const uint32_t last_cat_uidx = FindLast1BitBefore(cur_observed_cats, cur_cat_ct + 1);
            omitted_cat_uidxs[split_pheno_idx] = last_cat_uidx;
            ClearBit(last_cat_uidx, cur_observed_cats);
            --cur_observed_cat_ct;
          }
        } else {
          const uint32_t largest_cat_uidx = IdentifyRemainingCatsAndMostCommon(sample_include_intersect, cur_pheno_col, cur_nmiss_ct, cur_observed_cats, cat_obs_buf);
          cur_observed_cat_ct = PopcountWords(cur_observed_cats, cur_cat_ctl);
          if (cur_observed_cat_ct) {
            omitted_cat_uidxs[split_pheno_idx] = largest_cat_uidx;
            ClearBit(largest_cat_uidx, cur_observed_cats);
            --cur_observed_cat_ct;
          }
        }
        if (cur_observed_cat_ct) {
          // old phenotype name, '=' character, null terminator
          const uintptr_t blen_base = strlen(&(old_pheno_names[split_pheno_uidx * old_max_pheno_name_blen])) + 2;
          const char* const* cat_names = cur_pheno_col->category_names;
          uintptr_t cat_uidx_base = 0;
          uintptr_t cur_observed_cats_bits = cur_observed_cats[0];
          uint32_t cat_uidx = 0;
          for (uint32_t cat_idx = 0; cat_idx != cur_observed_cat_ct; ++cat_idx) {
            cat_uidx = BitIter1(cur_observed_cats, &cat_uidx_base, &cur_observed_cats_bits);
            const char* cur_cat_name = cat_names[cat_uidx];
            const uint32_t cur_slen = strlen(cur_cat_name);
            if (unlikely(memchr(cur_cat_name, '=', cur_slen))) {
              logerrputs("Error: --split-cat-pheno category names may not contain the '=' character.\n");
              goto SplitCatPheno_ret_INCONSISTENT_INPUT;
            }
            const uintptr_t total_blen = cur_slen + blen_base;
            if (total_blen > new_max_pheno_name_blen) {
              new_max_pheno_name_blen = total_blen;
            }
          }
          if (cat_uidx > max_cat_uidx) {
            max_cat_uidx = cat_uidx;
          }
        }
        observed_cat_cts[split_pheno_idx] = cur_observed_cat_ct;
        create_pheno_ct += cur_observed_cat_ct;
        if (cat_obs_buf) {
          BigstackReset(cat_obs_buf);
        }
      }
      if (unlikely(new_max_pheno_name_blen > kMaxIdBlen)) {
        logerrputs("Error: New --split-cat-pheno phenotype name too long.  Shorten your phenotype\nor your category names.\n");
        goto SplitCatPheno_ret_INCONSISTENT_INPUT;
      }
      const uint32_t copy_pheno_ct = old_pheno_ct - split_pheno_ct;
      // before new_pheno_ct variable definition due to potential integer
      // overflow
      if (unlikely(create_pheno_ct + copy_pheno_ct > kMaxPhenoCt)) {
        logerrputs("Error: --split-cat-pheno would create too many phenotypes (" PROG_NAME_STR " is limited to\n" MAX_PHENO_CT_STR ").\n");
        goto SplitCatPheno_ret_INCONSISTENT_INPUT;
      }
      const uint32_t new_pheno_ct = create_pheno_ct + copy_pheno_ct;
      uintptr_t** write_data_ptrs;
      if (unlikely(bigstack_alloc_wp(max_cat_uidx + 1, &write_data_ptrs))) {
        goto SplitCatPheno_ret_NOMEM;
      }
      const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
      uint32_t new_data_word_ct = raw_sample_ctaw;
      if (is_covar) {
        new_data_word_ct = DblCtToVecCt(raw_sample_ct) * kWordsPerVec;
      }
      uintptr_t* omit_dummy = nullptr;
      if (is_omit) {
        if (unlikely(bigstack_alloc_w(new_data_word_ct, &omit_dummy))) {
          goto SplitCatPheno_ret_NOMEM;
        }
      }

      // second pass: allocate memory and actually create the new phenotypes
      char* new_pheno_names;
      if (unlikely(pgl_malloc(new_pheno_ct * new_max_pheno_name_blen, &new_pheno_names))) {
        goto SplitCatPheno_ret_NOMEM;
      }
      doomed_pheno_names = old_pheno_names;
      *xpheno_names_ptr = new_pheno_names;
      PhenoCol* new_pheno_cols;
      if (unlikely(pgl_malloc(new_pheno_ct * sizeof(PhenoCol), &new_pheno_cols))) {
        goto SplitCatPheno_ret_NOMEM;
      }
      doomed_pheno_cols = old_pheno_cols;
      doomed_pheno_ct = old_pheno_ct;
      *xpheno_cols_ptr = new_pheno_cols;
      *xpheno_ct_ptr = new_pheno_ct;
      *max_xpheno_name_blen_ptr = new_max_pheno_name_blen;
      uintptr_t pheno_read_idx_base = 0;
      uintptr_t phenos_to_split_inv_bits = ~phenos_to_split[0];
      for (uint32_t pheno_write_idx = 0; pheno_write_idx != copy_pheno_ct; ++pheno_write_idx) {
        const uintptr_t pheno_read_idx = BitIter0(phenos_to_split, &pheno_read_idx_base, &phenos_to_split_inv_bits);

        // manually move this data
        memcpy(&(new_pheno_cols[pheno_write_idx]), &(doomed_pheno_cols[pheno_read_idx]), sizeof(PhenoCol));
        doomed_pheno_cols[pheno_read_idx].nonmiss = nullptr;

        strcpy(&(new_pheno_names[pheno_write_idx * new_max_pheno_name_blen]), &(old_pheno_names[pheno_read_idx * old_max_pheno_name_blen]));
      }
      for (uint32_t pheno_write_idx = copy_pheno_ct; pheno_write_idx != new_pheno_ct; ++pheno_write_idx) {
        new_pheno_cols[pheno_write_idx].nonmiss = nullptr;
      }

      const uintptr_t new_pheno_bytes_req = (raw_sample_ctaw + new_data_word_ct) * sizeof(intptr_t);
      uint32_t pheno_write_idx = copy_pheno_ct;
      pheno_read_idx_base = 0;
      phenos_to_split_bits = phenos_to_split[0];
      for (uint32_t split_pheno_idx = 0; split_pheno_idx != split_pheno_ct; ++split_pheno_idx) {
        const uintptr_t pheno_read_idx = BitIter1(phenos_to_split, &pheno_read_idx_base, &phenos_to_split_bits);
        const uint32_t cur_pheno_write_ct = observed_cat_cts[split_pheno_idx];
        if (!cur_pheno_write_ct) {
          continue;
        }
        const uintptr_t* cur_observed_cats = observed_cats[split_pheno_idx];
        const PhenoCol* old_pheno_col = &(old_pheno_cols[pheno_read_idx]);
        BitvecAndCopy(sample_include, old_pheno_col->nonmiss, raw_sample_ctaw, sample_include_intersect);
        const char* old_pheno_name = &(old_pheno_names[pheno_read_idx * old_max_pheno_name_blen]);
        const uint32_t old_pheno_name_slen = strlen(old_pheno_name);
        const char* const* old_cat_names = old_pheno_col->category_names;
        uintptr_t orig_cat_idx_base;
        uintptr_t cur_observed_cats_bits;
        BitIter1Start(cur_observed_cats, 1, &orig_cat_idx_base, &cur_observed_cats_bits);
        for (uint32_t uii = 0; uii != cur_pheno_write_ct; ++uii, ++pheno_write_idx) {
          const uintptr_t orig_cat_idx = BitIter1(cur_observed_cats, &orig_cat_idx_base, &cur_observed_cats_bits);
          uintptr_t* new_pheno_data_iter;
          if (unlikely(vecaligned_malloc(new_pheno_bytes_req, &new_pheno_data_iter))) {
            goto SplitCatPheno_ret_NOMEM;
          }
          char* new_phenoname_write_iter = memcpyax(&(new_pheno_names[pheno_write_idx * new_max_pheno_name_blen]), old_pheno_name, old_pheno_name_slen, '=');
          strcpy(new_phenoname_write_iter, old_cat_names[orig_cat_idx]);
          PhenoCol* pheno_write_col = &(new_pheno_cols[pheno_write_idx]);
          pheno_write_col->nonmiss = new_pheno_data_iter;
          pheno_write_col->category_names = nullptr;
          pheno_write_col->type_code = S_CAST(PhenoDtype, is_covar);
          pheno_write_col->nonnull_category_ct = 0;
          memcpy(new_pheno_data_iter, sample_include_intersect, raw_sample_ctaw * sizeof(intptr_t));
          new_pheno_data_iter = &(new_pheno_data_iter[raw_sample_ctaw]);
          write_data_ptrs[orig_cat_idx] = new_pheno_data_iter;
          // assigning to one element of a union and reading from another with
          // a different type is undefined behavior in C++11
          if (!is_covar) {
            pheno_write_col->data.cc = new_pheno_data_iter;
            ZeroWArr(raw_sample_ctaw, new_pheno_data_iter);
          } else {
            double* pheno_qt = R_CAST(double*, new_pheno_data_iter);
            pheno_write_col->data.qt = pheno_qt;
            if (qt_12) {
              for (uint32_t ujj = 0; ujj != raw_sample_ct; ++ujj) {
                pheno_qt[ujj] = 1.0;
              }
            } else {
              ZeroDArr(raw_sample_ct, pheno_qt);
            }
          }
        }
        if (is_omit) {
          write_data_ptrs[omitted_cat_uidxs[split_pheno_idx]] = omit_dummy;
        }

        const uint32_t cur_nmiss_ct = PopcountWords(sample_include_intersect, raw_sample_ctl);
        const uint32_t* cur_cats = old_pheno_col->data.cat;
        uintptr_t sample_uidx_base = 0;
        uintptr_t sample_include_intersect_bits = sample_include_intersect[0];
        if (!is_covar) {
          for (uint32_t sample_idx = 0; sample_idx != cur_nmiss_ct; ++sample_idx) {
            const uintptr_t sample_uidx = BitIter1(sample_include_intersect, &sample_uidx_base, &sample_include_intersect_bits);
            SetBit(sample_uidx, write_data_ptrs[cur_cats[sample_uidx]]);
          }
        } else {
          double** write_qt_ptrs = R_CAST(double**, write_data_ptrs);
          const double write_val = u31tod(1 + qt_12);
          for (uint32_t sample_idx = 0; sample_idx != cur_nmiss_ct; ++sample_idx) {
            const uintptr_t sample_uidx = BitIter1(sample_include_intersect, &sample_uidx_base, &sample_include_intersect_bits);
            write_qt_ptrs[cur_cats[sample_uidx]][sample_uidx] = write_val;
          }
        }
      }

      // if any preexisting phenotype names contain a single copy of the '='
      // character, verify that we didn't create any duplicate IDs
      for (uint32_t pheno_idx = 0; pheno_idx != copy_pheno_ct; ++pheno_idx) {
        const char* first_eq = strchr(&(new_pheno_names[pheno_idx * new_max_pheno_name_blen]), '=');
        if (first_eq && (!strchr(&(first_eq[1]), '='))) {
          uint32_t* id_htable;
          uint32_t id_htable_size;
          if (unlikely(HtableGoodSizeAlloc(new_pheno_ct, bigstack_left(), &id_htable, &id_htable_size))) {
            goto SplitCatPheno_ret_NOMEM;
          }
          uint32_t duplicate_idx = PopulateStrboxHtable(new_pheno_names, new_pheno_ct, new_max_pheno_name_blen, id_htable_size, id_htable);
          if (unlikely(duplicate_idx)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Duplicate %s '%s' created by --split-cat-pheno.\n", is_covar? "covariate" : "phenotype", &(new_pheno_names[duplicate_idx * new_max_pheno_name_blen]));
            goto SplitCatPheno_ret_INCONSISTENT_INPUT_WW;
          }
          break;
        }
      }

      free_const(doomed_pheno_names);
      doomed_pheno_names = nullptr;
      CleanupPhenoCols(doomed_pheno_ct, doomed_pheno_cols);
      doomed_pheno_cols = nullptr;

      logprintfww("--split-cat-pheno: %u categorical %s%s converted to %" PRIuPTR " %s%s.\n", split_pheno_ct, is_covar? "covariate" : "phenotype", (split_pheno_ct == 1)? "" : "s", create_pheno_ct, is_covar? "covariate" : "phenotype", (create_pheno_ct == 1)? "" : "s");
    }
    if (!at_least_one_cat_pheno_processed) {
      logerrprintf("Warning: No categorical phenotypes%s processed by --split-cat-pheno.\n", (split_cat_phenonames_flattened && (!(*covar_ct_ptr)))? "/covariates" : "");
    }
  }
  while (0) {
  SplitCatPheno_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SplitCatPheno_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  SplitCatPheno_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  BigstackReset(bigstack_mark);
  free_cond(doomed_pheno_names);
  CleanupPhenoCols(doomed_pheno_ct, doomed_pheno_cols);
  if (reterr) {
    if (*pheno_names_ptr) {
      free(*pheno_names_ptr);
      *pheno_names_ptr = nullptr;
    }
    CleanupPhenoCols(*pheno_ct_ptr, *pheno_cols_ptr);
    *pheno_ct_ptr = 0;
    *pheno_cols_ptr = nullptr;
    if (*covar_names_ptr) {
      free(*covar_names_ptr);
      *covar_names_ptr = nullptr;
    }
    CleanupPhenoCols(*covar_ct_ptr, *covar_cols_ptr);
    *covar_ct_ptr = 0;
    *covar_cols_ptr = nullptr;
  }
  return reterr;
}

PglErr PhenoVarianceStandardize(const char* vstd_flattened, const uintptr_t* sample_include, const char* pheno_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uint32_t is_covar_flag, PhenoCol* pheno_cols) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    if (!pheno_ct) {
      goto PhenoVarianceStandardize_ret_SKIP;
    }
    const uint32_t pheno_ctl = BitCtToWordCt(pheno_ct);
    uintptr_t* phenos_to_transform;
    if (unlikely(bigstack_calloc_w(pheno_ctl, &phenos_to_transform))) {
      goto PhenoVarianceStandardize_ret_NOMEM;
    }
    if (!vstd_flattened) {
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
        if (cur_pheno_col->type_code == kPhenoDtypeQt) {
          SetBit(pheno_idx, phenos_to_transform);
        }
      }
    } else {
      uint32_t* id_htable;
      uint32_t id_htable_size;
      if (unlikely(HtableGoodSizeAlloc(pheno_ct, bigstack_left(), &id_htable, &id_htable_size))) {
        goto PhenoVarianceStandardize_ret_NOMEM;
      }
      PopulateStrboxHtable(pheno_names, pheno_ct, max_pheno_name_blen, id_htable_size, id_htable);
      const char* vstd_phenonames_iter = vstd_flattened;
      do {
        const uint32_t cur_phenoname_slen = strlen(vstd_phenonames_iter);
        if (cur_phenoname_slen < max_pheno_name_blen) {
          uint32_t pheno_idx = StrboxHtableFind(vstd_phenonames_iter, pheno_names, id_htable, max_pheno_name_blen, cur_phenoname_slen, id_htable_size);
          if (pheno_idx != UINT32_MAX) {
            if (unlikely(pheno_cols[pheno_idx].type_code != kPhenoDtypeQt)) {
              snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a quantitative %s.\n", vstd_phenonames_iter, is_covar? "covariate" : "phenotype");
              goto PhenoVarianceStandardize_ret_INCONSISTENT_INPUT_WW;
            }
            SetBit(pheno_idx, phenos_to_transform);
          }
        }
        vstd_phenonames_iter = &(vstd_phenonames_iter[cur_phenoname_slen + 1]);
      } while (*vstd_phenonames_iter);
      BigstackReset(id_htable);
    }
    const uint32_t pheno_transform_ct = PopcountWords(phenos_to_transform, pheno_ctl);
    if (!pheno_transform_ct) {
      goto PhenoVarianceStandardize_ret_SKIP;
    }
    double* shifted_pheno_qt;
    if (unlikely(bigstack_alloc_d(raw_sample_ct, &shifted_pheno_qt))) {
      goto PhenoVarianceStandardize_ret_NOMEM;
    }
    const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
    uintptr_t pheno_uidx_base = 0;
    uintptr_t phenos_to_transform_bits = phenos_to_transform[0];
    for (uint32_t pheno_transform_idx = 0; pheno_transform_idx != pheno_transform_ct; ++pheno_transform_idx) {
      const uintptr_t pheno_uidx = BitIter1(phenos_to_transform, &pheno_uidx_base, &phenos_to_transform_bits);
      PhenoCol* cur_pheno_col = &(pheno_cols[pheno_uidx]);
      uintptr_t* pheno_nm = cur_pheno_col->nonmiss;
      BitvecAnd(sample_include, raw_sample_ctaw, pheno_nm);
      const uint32_t cur_sample_ct = PopcountWords(pheno_nm, raw_sample_ctaw);
      if (cur_sample_ct < 2) {
        if (cur_sample_ct) {
          logerrprintfww("Warning: Exactly one value present for %s '%s'; standardizing to missing.\n", is_covar? "covariate" : "quantitative phenotype", &(pheno_names[pheno_uidx * max_pheno_name_blen]));
          ZeroWArr(raw_sample_ctaw, pheno_nm);
        }
        continue;
      }
      double* pheno_qt = cur_pheno_col->data.qt;
      double shifted_pheno_sum = 0.0;
      double shifted_pheno_ssq = 0.0;
      const uint32_t first_sample_uidx = AdvTo1Bit(pheno_nm, 0);
      shifted_pheno_qt[first_sample_uidx] = 0.0;
      const double pheno_shift = pheno_qt[first_sample_uidx];
      uintptr_t sample_uidx_base;
      uintptr_t pheno_nm_bits;
      BitIter1Start(pheno_nm, first_sample_uidx + 1, &sample_uidx_base, &pheno_nm_bits);
      for (uint32_t sample_idx = 1; sample_idx != cur_sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(pheno_nm, &sample_uidx_base, &pheno_nm_bits);
        const double cur_shifted_pheno_val = pheno_qt[sample_uidx] - pheno_shift;
        shifted_pheno_sum += cur_shifted_pheno_val;
        shifted_pheno_ssq += cur_shifted_pheno_val * cur_shifted_pheno_val;
        shifted_pheno_qt[sample_uidx] = cur_shifted_pheno_val;
      }
      const double cur_shifted_mean = shifted_pheno_sum / u31tod(cur_sample_ct);
      const double variance_numer = shifted_pheno_ssq - shifted_pheno_sum * cur_shifted_mean;
      if (!(variance_numer > 0.0)) {
        logerrprintfww("Warning: %s '%s' is constant; standardizing to all-missing.\n", is_covar? "Covariate" : "Quantitative phenotype", &(pheno_names[pheno_uidx * max_pheno_name_blen]));
        ZeroWArr(raw_sample_ctaw, pheno_nm);
        continue;
      }
      const double cur_stdev_recip = sqrt(u31tod(cur_sample_ct - 1) / variance_numer);
      BitIter1Start(pheno_nm, first_sample_uidx, &sample_uidx_base, &pheno_nm_bits);
      for (uint32_t sample_idx = 0; sample_idx != cur_sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(pheno_nm, &sample_uidx_base, &pheno_nm_bits);
        pheno_qt[sample_uidx] = (shifted_pheno_qt[sample_uidx] - cur_shifted_mean) * cur_stdev_recip;
      }
    }
    // could reduce the reported number when all values were originally missing
    logprintf("--%svariance-standardize: %u %s%s transformed.\n", is_covar_flag? "covar-" : "", pheno_transform_ct, is_covar? "covariate" : "phenotype", (pheno_transform_ct == 1)? "" : "s");
  }
  while (0) {
  PhenoVarianceStandardize_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  PhenoVarianceStandardize_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  PhenoVarianceStandardize_ret_SKIP:
    logprintf("--%svariance-standardize: No %s affected.\n", is_covar_flag? "covar-" : "", is_covar? "covariates" : "quantitative phenotypes");
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct DblIndexStruct {
  double dxx;
  uint32_t uii;
#ifdef __cplusplus
  bool operator<(const struct DblIndexStruct& rhs) const {
    return dxx < rhs.dxx;
  }
#endif
} DblIndex;

PglErr PhenoQuantileNormalize(const char* quantnorm_flattened, const uintptr_t* sample_include, const char* pheno_names, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uint32_t is_subset_flag, PhenoCol* pheno_cols) {
  unsigned char* bigstack_mark = g_bigstack_base;
  const char* flag_prefix = is_subset_flag? (is_covar? "covar-" : "pheno-") : "";
  PglErr reterr = kPglRetSuccess;
  {
    if ((!pheno_ct) || (!sample_ct)) {
      goto PhenoQuantileNormalize_ret_SKIP;
    }
    // this boilerplate probably belongs in its own function
    const uint32_t pheno_ctl = BitCtToWordCt(pheno_ct);
    uintptr_t* phenos_to_transform;
    if (unlikely(bigstack_calloc_w(pheno_ctl, &phenos_to_transform))) {
      goto PhenoQuantileNormalize_ret_NOMEM;
    }
    if (!quantnorm_flattened) {
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
        if (cur_pheno_col->type_code == kPhenoDtypeQt) {
          SetBit(pheno_idx, phenos_to_transform);
        }
      }
    } else {
      uint32_t* id_htable;
      uint32_t id_htable_size;
      if (unlikely(HtableGoodSizeAlloc(pheno_ct, bigstack_left(), &id_htable, &id_htable_size))) {
        goto PhenoQuantileNormalize_ret_NOMEM;
      }
      PopulateStrboxHtable(pheno_names, pheno_ct, max_pheno_name_blen, id_htable_size, id_htable);
      const char* quantnorm_phenonames_iter = quantnorm_flattened;
      do {
        const uint32_t cur_phenoname_slen = strlen(quantnorm_phenonames_iter);
        if (cur_phenoname_slen < max_pheno_name_blen) {
          uint32_t pheno_idx = StrboxHtableFind(quantnorm_phenonames_iter, pheno_names, id_htable, max_pheno_name_blen, cur_phenoname_slen, id_htable_size);
          if (pheno_idx != UINT32_MAX) {
            if (unlikely(pheno_cols[pheno_idx].type_code != kPhenoDtypeQt)) {
              snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a quantitative %s.\n", quantnorm_phenonames_iter, is_covar? "covariate" : "phenotype");
              goto PhenoQuantileNormalize_ret_INCONSISTENT_INPUT_WW;
            }
            SetBit(pheno_idx, phenos_to_transform);
          }
        }
        quantnorm_phenonames_iter = &(quantnorm_phenonames_iter[cur_phenoname_slen + 1]);
      } while (*quantnorm_phenonames_iter);
      BigstackReset(id_htable);
    }
    const uint32_t pheno_transform_ct = PopcountWords(phenos_to_transform, pheno_ctl);
    if (!pheno_transform_ct) {
      goto PhenoQuantileNormalize_ret_SKIP;
    }
    DblIndex* tagged_raw_pheno_vals;
    if (unlikely(BIGSTACK_ALLOC_X(DblIndex, sample_ct, &tagged_raw_pheno_vals))) {
      goto PhenoQuantileNormalize_ret_NOMEM;
    }
    const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
    uintptr_t pheno_uidx_base = 0;
    uintptr_t phenos_to_transform_bits = phenos_to_transform[0];
    for (uint32_t pheno_transform_idx = 0; pheno_transform_idx != pheno_transform_ct; ++pheno_transform_idx) {
      const uintptr_t pheno_uidx = BitIter1(phenos_to_transform, &pheno_uidx_base, &phenos_to_transform_bits);
      PhenoCol* cur_pheno_col = &(pheno_cols[pheno_uidx]);
      uintptr_t* pheno_nm = cur_pheno_col->nonmiss;
      BitvecAnd(sample_include, raw_sample_ctaw, pheno_nm);
      const uint32_t cur_sample_ct = PopcountWords(pheno_nm, raw_sample_ctaw);
      if (!cur_sample_ct) {
        continue;
      }
      double* pheno_qt = cur_pheno_col->data.qt;
      uintptr_t sample_uidx_base = 0;
      uintptr_t pheno_nm_bits = pheno_nm[0];
      for (uint32_t sample_idx = 0; sample_idx != cur_sample_ct; ++sample_idx) {
        // bugfix (1 Sep 2017): this needs to iterate over pheno_nm, not
        // sample_include
        const uintptr_t sample_uidx = BitIter1(pheno_nm, &sample_uidx_base, &pheno_nm_bits);
        tagged_raw_pheno_vals[sample_idx].dxx = pheno_qt[sample_uidx];
        tagged_raw_pheno_vals[sample_idx].uii = sample_uidx;
      }
      STD_SORT_PAR_UNSEQ(cur_sample_ct, double_cmp, tagged_raw_pheno_vals);
      const double sample_ct_x2_recip = 1.0 / S_CAST(double, 2 * cur_sample_ct);
      for (uint32_t sample_idx_start = 0; sample_idx_start != cur_sample_ct; ) {
        const double cur_raw_pheno = tagged_raw_pheno_vals[sample_idx_start].dxx;
        uint32_t sample_idx_end = sample_idx_start + 1;
        for (; sample_idx_end != cur_sample_ct; ++sample_idx_end) {
          if (tagged_raw_pheno_vals[sample_idx_end].dxx != cur_raw_pheno) {
            break;
          }
        }
        const double cur_zscore = QuantileToZscore(S_CAST(double, sample_idx_start + sample_idx_end) * sample_ct_x2_recip);
        for (; sample_idx_start != sample_idx_end; ++sample_idx_start) {
          pheno_qt[tagged_raw_pheno_vals[sample_idx_start].uii] = cur_zscore;
        }
      }
    }
    logprintf("--%squantile-normalize: %u %s%s transformed.\n", flag_prefix, pheno_transform_ct, is_covar? "covariate" : "phenotype", (pheno_transform_ct == 1)? "" : "s");
  }
  while (0) {
  PhenoQuantileNormalize_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  PhenoQuantileNormalize_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  PhenoQuantileNormalize_ret_SKIP:
    logprintf("--%squantile-normalize: No %s affected.\n", flag_prefix, is_covar? "covariates" : "quantitative phenotypes");
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}


PglErr ProcessBoundaryToken(const char* tok_start, const char* tok_end, const char* token_source_str, uint32_t max_boundary_ct, PglErr err_type, double* prev_boundary_ptr, uint32_t* boundary_ct_ptr, double** freq_bounds_ptr, uint64_t** ddosage_bounds_ptr) {
  double cur_boundary;
  const char* scan_end = ScanadvDouble(tok_start, &cur_boundary);
  if (unlikely((!scan_end) || (scan_end != tok_end))) {
    logerrprintf("Error: Invalid token in %s.\n", token_source_str);
    return err_type;
  }
  if (unlikely(cur_boundary <= (*prev_boundary_ptr))) {
    logerrputs("Error: Invalid bin boundary sequence (must be strictly increasing, and start\nwith a positive number).\n");
    return err_type;
  }
  uint32_t boundary_ct = *boundary_ct_ptr;
  if (unlikely(boundary_ct == max_boundary_ct)) {
#ifdef __LP64__
    if (max_boundary_ct == 0x40000000) {
      logerrputs("Error: Too many bin boundaries.\n");
      return err_type;
    }
#endif
    return kPglRetNomem;
  }
  if (freq_bounds_ptr) {
    if (unlikely(cur_boundary > 1.0)) {
      logerrputs("Error: --freq bin boundary too large (must be <= 1).\n");
      return err_type;
    }
    // strictly-greater-than comparisons
    (*freq_bounds_ptr)[boundary_ct] = cur_boundary * (1 - kSmallEpsilon);
  } else {
    // max 2^31 - 3 variants
    if (unlikely(cur_boundary > 4294967290.0)) {
      logerrputs("Error: --freq counts bin boundary too large.\n");
      return err_type;
    }
    // due to the use of strictly-greater-than for comparison, we round
    // exact multiples of 1/32768 down
    const int64_t int_part = S_CAST(int64_t, cur_boundary);
    const double cur_boundary_frac_part = cur_boundary - int_part;
    const int64_t int_part_scaled = int_part * kDosageMax;
    if (cur_boundary_frac_part == 0.0) {
      (*ddosage_bounds_ptr)[boundary_ct] = int_part_scaled - 1;
    } else {
      (*ddosage_bounds_ptr)[boundary_ct] = int_part_scaled + S_CAST(int64_t, cur_boundary_frac_part * (kDosageMax * (1 - kSmallEpsilon)));
    }
  }
  *prev_boundary_ptr = cur_boundary;
  *boundary_ct_ptr = boundary_ct + 1;
  return kPglRetSuccess;
}

PglErr InitHistogramFromFileOrCommalist(const char* binstr, uint32_t is_fname, double** freq_bounds_ptr, uint64_t** ddosage_bounds_ptr, uint32_t* boundary_ct_ptr, uint32_t** histogram_ptr) {
  unsigned char* bigstack_end_mark = g_bigstack_end;
  TokenStream tks;
  PreinitTokenStream(&tks);
  uint32_t max_boundary_ct = 0;
  PglErr reterr = kPglRetSuccess;
  {
    if (is_fname) {
      // we want to accept >100000 numbers on a single line.  this will reject
      // "000...{a million more zeroes}...1"; pretty sure that's okay.
      reterr = InitTokenStreamEx(binstr, 1, 1, &tks);
      if (unlikely(reterr)) {
        goto InitHistogramFromFileOrCommalist_ret_TKSTREAM_FAIL;
      }
    }
    uintptr_t ulii = RoundDownPow2(bigstack_left(), kCacheline);
    if (unlikely(ulii < 2 * kCacheline)) {
      goto InitHistogramFromFileOrCommalist_ret_NOMEM;
    }
    // 12 = 8 bytes for boundary value + 4 bytes for histogram entry
    ulii = (ulii - 2 * kCacheline) / 12;
#ifdef __LP64__
    max_boundary_ct = MINV(ulii, 0x40000000);
#else
    max_boundary_ct = ulii;
#endif
    if (freq_bounds_ptr) {
      *freq_bounds_ptr = R_CAST(double*, g_bigstack_base);
    } else {
      *ddosage_bounds_ptr = R_CAST(uint64_t*, g_bigstack_base);
    }
    uint32_t boundary_ct = 0;
    double prev_boundary = 0.0;
    if (is_fname) {
      while (1) {
        char* shard_boundaries[2];
        reterr = TksNext(&tks, 1, shard_boundaries);
        if (reterr) {
          break;
        }
        char* shard_iter = shard_boundaries[0];
        char* shard_end = shard_boundaries[1];
        while (1) {
          shard_iter = FirstPostspaceBounded(shard_iter, shard_end);
          if (shard_iter == shard_end) {
            break;
          }
          char* token_end = CurTokenEnd(shard_iter);
          reterr = ProcessBoundaryToken(shard_iter, token_end, binstr, max_boundary_ct, kPglRetMalformedInput, &prev_boundary, &boundary_ct, freq_bounds_ptr, ddosage_bounds_ptr);
          if (unlikely(reterr)) {
            goto InitHistogramFromFileOrCommalist_ret_1;
          }
          shard_iter = token_end;
        }
      }
      if (unlikely(reterr != kPglRetEof)) {
        goto InitHistogramFromFileOrCommalist_ret_TKSTREAM_FAIL;
      }
      if (CleanupTokenStream3("--freq {ref|alt1}bins-file= file", &tks, &reterr)) {
        goto InitHistogramFromFileOrCommalist_ret_1;
      }
    } else {
      const char* binstr_iter = binstr;
      while (1) {
        const char* tok_end = strchrnul(binstr_iter, ',');
        reterr = ProcessBoundaryToken(binstr_iter, tok_end, "--freq {ref,alt1}bins= list", max_boundary_ct, kPglRetInvalidCmdline, &prev_boundary, &boundary_ct, freq_bounds_ptr, ddosage_bounds_ptr);
        if (unlikely(reterr)) {
          goto InitHistogramFromFileOrCommalist_ret_1;
        }
        if (!(*tok_end)) {
          break;
        }
        binstr_iter = &(tok_end[1]);
      }
    }
    *boundary_ct_ptr = boundary_ct;
    g_bigstack_base += RoundUpPow2(boundary_ct * (8 * k1LU), kCacheline);
    *histogram_ptr = S_CAST(uint32_t*, bigstack_alloc_raw_rd((boundary_ct + 1) * sizeof(int32_t)));
    ZeroU32Arr(boundary_ct + 1, *histogram_ptr);
  }
  while (0) {
  InitHistogramFromFileOrCommalist_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  InitHistogramFromFileOrCommalist_ret_TKSTREAM_FAIL:
    TokenStreamErrPrint("--freq {ref|alt1}bins-file= file", &tks);
    break;
  }
 InitHistogramFromFileOrCommalist_ret_1:
  CleanupTokenStream2("--freq {ref|alt1}bins-file= file", &tks, &reterr);
  BigstackEndReset(bigstack_end_mark);
  return reterr;
}

PglErr WriteAlleleFreqs(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint64_t* founder_allele_ddosages, const double* imp_r2_vals, const char* ref_binstr, const char* alt1_binstr, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, FreqRptFlags freq_rpt_flags, uint32_t max_thread_ct, uint32_t nonfounders, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t counts = (freq_rpt_flags / kfAlleleFreqCounts) & 1;
    if (counts) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".acount");
    } else {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".afreq");
    }
    if (!(freq_rpt_flags & kfAlleleFreqBinsOnly)) {
      const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
      const uintptr_t overflow_buf_size = kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + max_allele_ct * (24 * k1LU) + 2 * max_allele_slen;
      const uint32_t output_zst = freq_rpt_flags & kfAlleleFreqZs;
      if (output_zst) {
        snprintf(&(outname_end[6 + counts]), kMaxOutfnameExtBlen - 7, ".zst");
      }
      reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
      if (unlikely(reterr)) {
        goto WriteAlleleFreqs_ret_1;
      }
      *cswritep++ = '#';
      const uint32_t chr_col = freq_rpt_flags & kfAlleleFreqColChrom;

      // includes trailing tab
      char* chr_buf;
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto WriteAlleleFreqs_ret_NOMEM;
      }
      if (chr_col) {
        cswritep = strcpya_k(cswritep, "CHROM\t");
      }
      if (freq_rpt_flags & kfAlleleFreqColPos) {
        cswritep = strcpya_k(cswritep, "POS\t");
      } else {
        variant_bps = nullptr;
      }
      cswritep = strcpya_k(cswritep, "ID");
      const uint32_t ref_col = freq_rpt_flags & kfAlleleFreqColRef;
      if (ref_col) {
        cswritep = strcpya_k(cswritep, "\tREF");
      }
      const uint32_t alt1_col = freq_rpt_flags & kfAlleleFreqColAlt1;
      if (alt1_col) {
        cswritep = strcpya_k(cswritep, "\tALT1");
      }
      const uint32_t alt_col = freq_rpt_flags & kfAlleleFreqColAlt;
      if (alt_col) {
        cswritep = strcpya_k(cswritep, "\tALT");
      }
      const uint32_t reffreq_col = freq_rpt_flags & kfAlleleFreqColReffreq;
      if (reffreq_col) {
        cswritep = strcpya_k(cswritep, "\tREF_");
        if (counts) {
          cswritep = strcpya_k(cswritep, "CT");
        } else {
          cswritep = strcpya_k(cswritep, "FREQ");
        }
      }
      const uint32_t alt1freq_col = freq_rpt_flags & kfAlleleFreqColAlt1freq;
      if (alt1freq_col) {
        cswritep = strcpya_k(cswritep, "\tALT1_");
        if (counts) {
          cswritep = strcpya_k(cswritep, "CT");
        } else {
          cswritep = strcpya_k(cswritep, "FREQ");
        }
      }
      const uint32_t freq_col = freq_rpt_flags & (kfAlleleFreqColFreq | kfAlleleFreqColAltfreq);
      const uint32_t commalist_exclude_ref = (freq_rpt_flags & (kfAlleleFreqColAltfreq | kfAlleleFreqColAlteq | kfAlleleFreqColAlteqz | kfAlleleFreqColAltnumeq))? 1 : 0;
      const uint32_t eq_col = freq_rpt_flags & (kfAlleleFreqColEq | kfAlleleFreqColEqz | kfAlleleFreqColAlteq | kfAlleleFreqColAlteqz | kfAlleleFreqColNumeq | kfAlleleFreqColAltnumeq);
      const uint32_t eq_includez = freq_rpt_flags & (kfAlleleFreqColEqz | kfAlleleFreqColAlteqz);
      const uint32_t eq_num = freq_rpt_flags & (kfAlleleFreqColNumeq | kfAlleleFreqColAltnumeq);
      if (freq_col || eq_col) {
        *cswritep++ = '\t';
        if (commalist_exclude_ref) {
          cswritep = strcpya_k(cswritep, "ALT_");
        }
        if (eq_num) {
          cswritep = strcpya_k(cswritep, "NUM_");
        }
        if (counts) {
          cswritep = strcpya_k(cswritep, "CTS");
        } else {
          cswritep = strcpya_k(cswritep, "FREQS");
        }
      }
      const uint32_t imp_r2_col = freq_rpt_flags & (kfAlleleFreqColMachR2 | kfAlleleFreqColMinimac3R2);
      if (imp_r2_col) {
        // These two columns are currently mutually exclusive.
        if (freq_rpt_flags & kfAlleleFreqColMachR2) {
          cswritep = strcpya_k(cswritep, "\tMACH_R2");
        } else {
          cswritep = strcpya_k(cswritep, "\tMINIMAC3_R2");
        }
      }
      const uint32_t nobs_col = freq_rpt_flags & kfAlleleFreqColNobs;
      if (nobs_col) {
        cswritep = strcpya_k(cswritep, "\tOBS_CT");
      }
      AppendBinaryEoln(&cswritep);

      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t suppress_imp_r2 = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uint32_t cur_allele_ct = 2;
      printf("--freq%s%s: 0%%", output_zst? " zs" : "", counts? " counts" : "");
      fflush(stdout);
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        if (variant_uidx >= chr_end) {
          do {
            ++chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          } while (variant_uidx >= chr_end);
          const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
          char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
          suppress_imp_r2 = (chr_idx == x_code) || (chr_idx == mt_code);
          *chr_name_end = '\t';
          chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
        }
        if (chr_col) {
          cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
        }
        if (variant_bps) {
          cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
        }
        cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        if (ref_col) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, cur_alleles[0]);
        }
        if (alt1_col) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, cur_alleles[1]);
        }
        if (alt_col) {
          *cswritep++ = '\t';
          for (uint32_t allele_idx = 1; allele_idx != cur_allele_ct; ++allele_idx) {
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto WriteAlleleFreqs_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
          }
          --cswritep;
        }
        const uint64_t* cur_allele_ddosages = &(founder_allele_ddosages[allele_idx_offset_base]);
        uint64_t tot_allele_ddosage = cur_allele_ddosages[0];
        for (uint32_t allele_idx = 1; allele_idx != cur_allele_ct; ++allele_idx) {
          tot_allele_ddosage += cur_allele_ddosages[allele_idx];
        }
        double tot_allele_ddosage_recip = 0.0;
        if (!counts) {
          tot_allele_ddosage_recip = 1.0 / u63tod(tot_allele_ddosage);
        }
        if (reffreq_col) {
          *cswritep++ = '\t';
          const uint64_t ref_ddosage = cur_allele_ddosages[0];
          if (counts) {
            // update (22 Nov 2019): we want --read-freq to be able to exactly
            // replicate the original frequency.  dosagetoa() is not good
            // enough for this purpose.
            cswritep = ddosagetoa_full(ref_ddosage, cswritep);
          } else {
            cswritep = dtoa_g(u63tod(ref_ddosage) * tot_allele_ddosage_recip, cswritep);
          }
        }
        if (alt1freq_col) {
          *cswritep++ = '\t';
          const uint64_t alt1_ddosage = cur_allele_ddosages[1];
          if (counts) {
            cswritep = ddosagetoa_full(alt1_ddosage, cswritep);
          } else {
            cswritep = dtoa_g(u63tod(alt1_ddosage) * tot_allele_ddosage_recip, cswritep);
          }
        }
        if (freq_col) {
          *cswritep++ = '\t';
          for (uint32_t allele_idx = commalist_exclude_ref; allele_idx != cur_allele_ct; ++allele_idx) {
            const uint64_t cur_allele_ddosage = cur_allele_ddosages[allele_idx];
            if (counts) {
              cswritep = ddosagetoa_full(cur_allele_ddosage, cswritep);
            } else {
              cswritep = dtoa_g(u63tod(cur_allele_ddosage) * tot_allele_ddosage_recip, cswritep);
            }
            *cswritep++ = ',';
          }
          --cswritep;
        } else if (eq_col) {
          *cswritep++ = '\t';
          uint32_t at_least_one_entry = 0;
          for (uint32_t allele_idx = commalist_exclude_ref; allele_idx != cur_allele_ct; ++allele_idx) {
            const uint64_t cur_allele_ddosage = cur_allele_ddosages[allele_idx];
            if (eq_includez || cur_allele_ddosage) {
              if (eq_num) {
                cswritep = u32toa(allele_idx, cswritep);
              } else {
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto WriteAlleleFreqs_ret_WRITE_FAIL;
                }
                const char* cur_allele = cur_alleles[allele_idx];
                const char* cur_allele_end_or_eq = strchrnul(cur_allele, '=');
                if (unlikely(*cur_allele_end_or_eq == '=')) {
                  logerrputs("Error: --freq's 'eq', 'eqz', 'alteq', and 'alteqz' columns cannot be requested\nwhen an allele code contains a '='.\n");
                  goto WriteAlleleFreqs_ret_INCONSISTENT_INPUT;
                }
                cswritep = memcpya(cswritep, cur_allele, cur_allele_end_or_eq - cur_allele);
              }
              *cswritep++ = '=';
              if (counts) {
                cswritep = ddosagetoa_full(cur_allele_ddosage, cswritep);
              } else {
                cswritep = dtoa_g(u63tod(cur_allele_ddosage) * tot_allele_ddosage_recip, cswritep);
              }
              *cswritep++ = ',';
              at_least_one_entry = 1;
            }
          }
          if (at_least_one_entry) {
            --cswritep;
          } else {
            *cswritep++ = '.';
          }
        }
        if (imp_r2_col) {
          *cswritep++ = '\t';
          if (!suppress_imp_r2) {
            cswritep = dtoa_g(imp_r2_vals[variant_uidx], cswritep);
          } else {
            cswritep = strcpya_k(cswritep, "NA");
          }
        }
        if (nobs_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(tot_allele_ddosage / kDosageMax, cswritep);
        }
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto WriteAlleleFreqs_ret_WRITE_FAIL;
        }
        if (variant_idx >= next_print_variant_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (variant_idx * 100LLU) / variant_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
        }
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto WriteAlleleFreqs_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      logprintfww("--freq%s%s: Allele %s (%s) written to %s .\n", output_zst? " zs" : "", counts? " counts" : "", counts? "counts" : "frequencies", nonfounders? "all samples" : "founders only", outname);
    }

    if (ref_binstr || alt1_binstr) {
      BigstackReset(bigstack_mark);
      double* ref_freq_bounds = nullptr;
      uint64_t* ref_ddosage_bounds = nullptr;
      uint32_t* ref_histogram = nullptr;
      uint32_t ref_boundary_ct = 0;
      if (ref_binstr) {
        reterr = InitHistogramFromFileOrCommalist(ref_binstr, (freq_rpt_flags / kfAlleleFreqBinsRefFname) & 1, counts? nullptr : (&ref_freq_bounds), counts? (&ref_ddosage_bounds) : nullptr, &ref_boundary_ct, &ref_histogram);
        if (unlikely(reterr)) {
          goto WriteAlleleFreqs_ret_1;
        }
      }
      double* alt1_freq_bounds = nullptr;
      uint64_t* alt1_ddosage_bounds = nullptr;
      uint32_t* alt1_histogram = nullptr;
      uint32_t alt1_boundary_ct = 0;
      if (alt1_binstr) {
        reterr = InitHistogramFromFileOrCommalist(alt1_binstr, (freq_rpt_flags / kfAlleleFreqBinsAlt1Fname) & 1, counts? nullptr : (&alt1_freq_bounds), counts? (&alt1_ddosage_bounds) : nullptr, &alt1_boundary_ct, &alt1_histogram);
        if (unlikely(reterr)) {
          goto WriteAlleleFreqs_ret_1;
        }
      }
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      if (!counts) {
        uint32_t cur_allele_ct = 2;
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          uintptr_t allele_idx_offset_base = variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[variant_uidx];
            cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
          }
          const uint64_t* cur_allele_ddosages = &(founder_allele_ddosages[allele_idx_offset_base]);
          const uint64_t ref_allele_ddosage = cur_allele_ddosages[0];
          const uint64_t alt1_allele_ddosage = cur_allele_ddosages[1];
          uint64_t tot_allele_ddosage = ref_allele_ddosage + alt1_allele_ddosage;
          for (uint32_t allele_idx = 2; allele_idx != cur_allele_ct; ++allele_idx) {
            tot_allele_ddosage += cur_allele_ddosages[allele_idx];
          }
          const double tot_allele_ddosage_recip = 1.0 / u63tod(tot_allele_ddosage);
          if (ref_histogram) {
            ref_histogram[CountSortedSmallerD(ref_freq_bounds, ref_boundary_ct, ref_allele_ddosage * tot_allele_ddosage_recip)] += 1;
          }
          if (alt1_histogram) {
            alt1_histogram[CountSortedSmallerD(alt1_freq_bounds, alt1_boundary_ct, alt1_allele_ddosage * tot_allele_ddosage_recip)] += 1;
          }
        }
      } else {
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          uintptr_t allele_idx_offset_base = variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          }
          const uint64_t* cur_allele_ddosages = &(founder_allele_ddosages[allele_idx_offset_base]);
          if (ref_histogram) {
            ref_histogram[CountSortedSmallerU64(ref_ddosage_bounds, ref_boundary_ct, cur_allele_ddosages[0])] += 1;
          }
          if (alt1_histogram) {
            alt1_histogram[CountSortedSmallerU64(alt1_ddosage_bounds, alt1_boundary_ct, cur_allele_ddosages[1])] += 1;
          }
        }
      }
      for (uint32_t is_alt1 = 0; is_alt1 != 2; ++is_alt1) {
        const uint32_t* cur_histogram = is_alt1? alt1_histogram : ref_histogram;
        if (!cur_histogram) {
          continue;
        }
        char* outname_end2 = &(outname_end[6 + counts]);
        if (!is_alt1) {
          outname_end2 = strcpya_k(outname_end2, ".ref");
        } else {
          outname_end2 = strcpya_k(outname_end2, ".alt1");
        }
        snprintf(outname_end2, kMaxOutfnameExtBlen - 12, ".bins");
        if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
          goto WriteAlleleFreqs_ret_OPEN_FAIL;
        }
        char* textbuf = g_textbuf;
        char* textbuf_flush = &(textbuf[kMaxMediumLine]);
        char* write_iter = strcpya_k(textbuf, "#BIN_START\tOBS_CT" EOLN_STR);
        const uint32_t cur_boundary_ct = is_alt1? alt1_boundary_ct : ref_boundary_ct;
        if (!counts) {
          const double* cur_freq_bounds = is_alt1? alt1_freq_bounds : ref_freq_bounds;
          for (uint32_t bin_idx = 0; bin_idx <= cur_boundary_ct; ++bin_idx) {
            if (!bin_idx) {
              *write_iter++ = '0';
            } else {
              write_iter = dtoa_g(cur_freq_bounds[bin_idx - 1] * (1.0 / (1 - kSmallEpsilon)), write_iter);
            }
            *write_iter++ = '\t';
            write_iter = u32toa(cur_histogram[bin_idx], write_iter);
            AppendBinaryEoln(&write_iter);
            if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
              goto WriteAlleleFreqs_ret_WRITE_FAIL;
            }
          }
        } else {
          const uint64_t* cur_ddosage_bounds = is_alt1? alt1_ddosage_bounds : ref_ddosage_bounds;
          for (uint32_t bin_idx = 0; bin_idx <= cur_boundary_ct; ++bin_idx) {
            if (!bin_idx) {
              *write_iter++ = '0';
            } else {
              write_iter = ddosagetoa(1 + cur_ddosage_bounds[bin_idx - 1], write_iter);
            }
            *write_iter++ = '\t';
            write_iter = u32toa(cur_histogram[bin_idx], write_iter);
            AppendBinaryEoln(&write_iter);
            if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
              goto WriteAlleleFreqs_ret_WRITE_FAIL;
            }
          }
        }
        if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
          goto WriteAlleleFreqs_ret_WRITE_FAIL;
        }
        const uint32_t cur_is_file = freq_rpt_flags & (is_alt1? kfAlleleFreqBinsAlt1Fname : kfAlleleFreqBinsRefFname);
        logprintfww("--freq%s %sbins%s=: Histogram written to %s .\n", counts? " counts" : "", is_alt1? "alt1" : "ref", cur_is_file? "-file" : "", outname);
      }
    }
  }
  while (0) {
  WriteAlleleFreqs_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WriteAlleleFreqs_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  WriteAlleleFreqs_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  WriteAlleleFreqs_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 WriteAlleleFreqs_ret_1:
  CswriteCloseCond(&css, cswritep);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr WriteGenoCounts(const uintptr_t* sample_include, __attribute__((unused)) const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts), uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t variant_ct, uint32_t x_start, uint32_t max_allele_slen, GenoCountsFlags geno_counts_flags, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    char* chr_buf;
    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
      goto WriteGenoCounts_ret_NOMEM;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t max_allele_ct = PgrGetMaxAlleleCt(simple_pgrp);
    uint32_t* sample_include_cumulative_popcounts = nullptr;
    uintptr_t* sex_male_collapsed = nullptr;  // chrX
    uint32_t* sex_male_cumulative_popcounts = nullptr;  // chrY
    PgenVariant pgv;
    pgv.genovec = nullptr;
    pgv.patch_01_set = nullptr;
    pgv.patch_01_vals = nullptr;
    pgv.patch_10_set = nullptr;
    pgv.patch_10_vals = nullptr;
    uint32_t* diploid_pair_cts = nullptr;
    uint32_t* hap_cts = nullptr;
    const uint32_t more_counts_needed = (max_allele_ct > 2) && (geno_counts_flags & (kfGenoCountsColRefalt1 | kfGenoCountsColRefalt | kfGenoCountsColHomalt1 | kfGenoCountsColAltxy | kfGenoCountsColXy | kfGenoCountsColHapalt1 | kfGenoCountsColHapalt | kfGenoCountsColHap | kfGenoCountsColNumeq));
    if (more_counts_needed) {
      if (unlikely(
              bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
              bigstack_alloc_w(sample_ctl, &sex_male_collapsed) ||
              bigstack_alloc_u32(raw_sample_ctl, &sex_male_cumulative_popcounts) ||
              bigstack_alloc_w(NypCtToWordCt(raw_sample_ct), &pgv.genovec) ||
              bigstack_alloc_w(sample_ctl, &pgv.patch_01_set) ||
              bigstack_alloc_ac(sample_ct, &pgv.patch_01_vals) ||
              bigstack_alloc_w(sample_ctl, &pgv.patch_10_set) ||
              bigstack_alloc_ac(2 * sample_ct, &pgv.patch_10_vals) ||
              bigstack_alloc_u32(max_allele_ct * max_allele_ct, &diploid_pair_cts) ||
              bigstack_alloc_u32(max_allele_ct, &hap_cts))) {
        goto WriteGenoCounts_ret_NOMEM;
      }
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_male_collapsed);
      FillCumulativePopcounts(sex_male, raw_sample_ctl, sex_male_cumulative_popcounts);
    }
    // Will need to remove quadratic dependency on max_allele_ct if
    // sizeof(AlleleCode) > 1.
    // The actual number of numeq and {alt}xy entries is roughly half of
    // allele_ct^2, so 12 * allele_ct^2 allows 20 bytes per entry when both
    // fields are present simultaneously, which is more than enough.  (The
    // surplus and the 512 constant are enough to cover the linearly-growing
    // fields.)
    const uintptr_t overflow_buf_size = kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + max_allele_ct * max_allele_ct * 20 + 2 * max_allele_slen;
    const uint32_t output_zst = geno_counts_flags & kfGenoCountsZs;
    OutnameZstSet(".gcount", output_zst, outname_end);
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WriteGenoCounts_ret_1;
    }
    *cswritep++ = '#';
    const uint32_t chr_col = geno_counts_flags & kfGenoCountsColChrom;

    // includes trailing tab
    if (chr_col) {
      cswritep = strcpya_k(cswritep, "CHROM\t");
    }
    if (geno_counts_flags & kfGenoCountsColPos) {
      cswritep = strcpya_k(cswritep, "POS\t");
    } else {
      variant_bps = nullptr;
    }
    cswritep = strcpya_k(cswritep, "ID");
    const uint32_t ref_col = geno_counts_flags & kfGenoCountsColRef;
    if (ref_col) {
      cswritep = strcpya_k(cswritep, "\tREF");
    }
    const uint32_t alt1_col = geno_counts_flags & kfGenoCountsColAlt1;
    if (alt1_col) {
      cswritep = strcpya_k(cswritep, "\tALT1");
    }
    const uint32_t alt_col = geno_counts_flags & kfGenoCountsColAlt;
    if (alt_col) {
      cswritep = strcpya_k(cswritep, "\tALT");
    }
    const uint32_t homref_col = geno_counts_flags & kfGenoCountsColHomref;
    if (homref_col) {
      cswritep = strcpya_k(cswritep, "\tHOM_REF_CT");
    }
    const uint32_t refalt1_col = geno_counts_flags & kfGenoCountsColRefalt1;
    if (refalt1_col) {
      cswritep = strcpya_k(cswritep, "\tHET_REF_ALT1_CT");
    }
    const uint32_t refalt_col = geno_counts_flags & kfGenoCountsColRefalt;
    if (refalt_col) {
      cswritep = strcpya_k(cswritep, "\tHET_REF_ALT_CTS");
    }
    const uint32_t homalt1_col = geno_counts_flags & kfGenoCountsColHomalt1;
    if (homalt1_col) {
      cswritep = strcpya_k(cswritep, "\tHOM_ALT1_CT");
    }
    const uint32_t xy_col = geno_counts_flags & (kfGenoCountsColAltxy | kfGenoCountsColXy);
    const uint32_t xy_col_altonly = (geno_counts_flags / kfGenoCountsColAltxy) & 1;
    if (xy_col) {
      *cswritep++ = '\t';
      if (xy_col_altonly) {
        cswritep = strcpya_k(cswritep, "TWO_ALT_GENO_CTS");
      } else {
        cswritep = strcpya_k(cswritep, "DIPLOID_GENO_CTS");
      }
    }
    const uint32_t hapref_col = geno_counts_flags & kfGenoCountsColHapref;
    if (hapref_col) {
      cswritep = strcpya_k(cswritep, "\tHAP_REF_CT");
    }
    const uint32_t hapalt1_col = geno_counts_flags & kfGenoCountsColHapalt1;
    if (hapalt1_col) {
      cswritep = strcpya_k(cswritep, "\tHAP_ALT1_CT");
    }
    const uint32_t hap_col = geno_counts_flags & (kfGenoCountsColHapalt | kfGenoCountsColHap);
    const uint32_t hap_col_altonly = (geno_counts_flags / kfGenoCountsColHapalt) & 1;
    if (hap_col) {
      if (hap_col_altonly) {
        cswritep = strcpya_k(cswritep, "\tHAP_ALT_CTS");
      } else {
        cswritep = strcpya_k(cswritep, "\tHAP_CTS");
      }
    }
    const uint32_t numeq_col = geno_counts_flags & kfGenoCountsColNumeq;
    if (numeq_col) {
      cswritep = strcpya_k(cswritep, "\tGENO_NUM_CTS");
    }
    const uint32_t missing_col = geno_counts_flags & kfGenoCountsColMissing;
    if (missing_col) {
      cswritep = strcpya_k(cswritep, "\tMISSING_CT");
    }
    const uint32_t nobs_col = geno_counts_flags & kfGenoCountsColNobs;
    if (nobs_col) {
      cswritep = strcpya_k(cswritep, "\tOBS_CT");
    }
    AppendBinaryEoln(&cswritep);

    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    const uintptr_t* cur_sample_include = nullptr;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t is_autosomal_diploid = 0;
    uint32_t is_x = 0;
    uint32_t nobs_base = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t homref_ct = 0;
    uint32_t het_ct = 0;
    uint32_t homalt1_ct = 0;
    uint32_t hapref_ct = 0;
    uint32_t hapalt1_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    printf("--geno-counts%s: 0%%", output_zst? " zs" : "");
    fflush(stdout);
    PgrSampleSubsetIndex pssi;
    PgrClearSampleSubsetIndex(simple_pgrp, &pssi);
    uint32_t allele_ct = 2;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
        *chr_name_end = '\t';
        chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
        is_autosomal_diploid = !IsSet(cip->haploid_mask, chr_idx);
        nobs_base = sample_ct;
        is_x = (chr_idx == x_code);
        cur_sample_include = sample_include;
        const uint32_t* cur_cumulative_popcounts = sample_include_cumulative_popcounts;
        if (!is_autosomal_diploid) {
          if (chr_idx == y_code) {
            cur_sample_include = sex_male;
            cur_cumulative_popcounts = sex_male_cumulative_popcounts;
            nobs_base = male_ct;
          }
        }
        PgrSetSampleSubsetIndex(cur_cumulative_popcounts, simple_pgrp, &pssi);
        homref_ct = 0;
        het_ct = 0;
        homalt1_ct = 0;
        hapref_ct = 0;
        hapalt1_ct = 0;
      }
      if (chr_col) {
        cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
      }
      if (variant_bps) {
        cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
      }
      cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (ref_col) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[0]);
      }
      if (alt1_col) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cur_alleles[1]);
      }
      if (alt_col) {
        *cswritep++ = '\t';
        for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto WriteGenoCounts_ret_WRITE_FAIL;
          }
          cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
        }
        --cswritep;
      }
      STD_ARRAY_KREF(uint32_t, 3) cur_raw_geno_cts = raw_geno_cts[variant_uidx];
      uint32_t missing_ct;
      if ((allele_ct == 2) || (!more_counts_needed)) {
        if (is_autosomal_diploid) {
          homref_ct = cur_raw_geno_cts[0];
          het_ct = cur_raw_geno_cts[1];
          homalt1_ct = cur_raw_geno_cts[2];
          missing_ct = nobs_base - homref_ct - het_ct - homalt1_ct;
        } else {
          if (is_x) {
            // bugfix (8 Oct 2017): don't set these in haploid case
            homref_ct = cur_raw_geno_cts[0];
            het_ct = cur_raw_geno_cts[1];
            homalt1_ct = cur_raw_geno_cts[2];
            if (x_male_geno_cts) {
              STD_ARRAY_KREF(uint32_t, 3) cur_male_geno_cts = x_male_geno_cts[variant_uidx - x_start];
              hapref_ct = cur_male_geno_cts[0];
              homref_ct -= hapref_ct;
              het_ct -= cur_male_geno_cts[1];
              hapalt1_ct = cur_male_geno_cts[2];
              homalt1_ct -= hapalt1_ct;
            }
            missing_ct = nobs_base - homref_ct - het_ct - homalt1_ct - hapref_ct - hapalt1_ct;
          } else {
            // chrY or other pure-haploid; hethap treated as missing
            hapref_ct = cur_raw_geno_cts[0];
            hapalt1_ct = cur_raw_geno_cts[2];
            missing_ct = nobs_base - hapref_ct - hapalt1_ct;
          }
        }
        if (homref_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(homref_ct, cswritep);
        }
        if (refalt1_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(het_ct, cswritep);
        }
        if (refalt_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(het_ct, cswritep);
        }
        if (homalt1_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(homalt1_ct, cswritep);
        }
        if (xy_col_altonly) {
          *cswritep++ = '\t';
          cswritep = u32toa(homalt1_ct, cswritep);
        } else if (xy_col) {
          *cswritep++ = '\t';
          cswritep = u32toa_x(homref_ct, ',', cswritep);
          cswritep = u32toa_x(het_ct, ',', cswritep);
          cswritep = u32toa(homalt1_ct, cswritep);
        }
        if (hapref_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(hapref_ct, cswritep);
        }
        if (hapalt1_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(hapalt1_ct, cswritep);
        }
        if (hap_col) {
          *cswritep++ = '\t';
          if (!hap_col_altonly) {
            cswritep = u32toa_x(hapref_ct, ',', cswritep);
          }
          cswritep = u32toa(hapalt1_ct, cswritep);
        }
        if (numeq_col) {
          *cswritep++ = '\t';
          if (homref_ct) {
            cswritep = strcpya_k(cswritep, "0/0=");
            cswritep = u32toa_x(homref_ct, ',', cswritep);
          }
          if (het_ct) {
            cswritep = strcpya_k(cswritep, "0/1=");
            cswritep = u32toa_x(het_ct, ',', cswritep);
          }
          if (homalt1_ct) {
            cswritep = strcpya_k(cswritep, "1/1=");
            cswritep = u32toa_x(homalt1_ct, ',', cswritep);
          }
          if (hapref_ct) {
            cswritep = strcpya_k(cswritep, "0=");
            cswritep = u32toa_x(hapref_ct, ',', cswritep);
          }
          if (hapalt1_ct) {
            cswritep = strcpya_k(cswritep, "1=");
            cswritep = u32toa_x(hapalt1_ct, ',', cswritep);
          }
          if (missing_ct != nobs_base) {
            --cswritep;
          } else {
            *cswritep++ = '.';
          }
        }
      } else {
        reterr = PgrGetM(cur_sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
        if (unlikely(reterr)) {
          goto WriteGenoCounts_ret_PGR_FAIL;
        }
        // Usually don't care about contents of genovec, patch_01_set, and
        // patch_10_set, but chrX is an exception.
        // Need to have an alternate strategy once sizeof(AlleleCode) > 1.
        ZeroU32Arr(allele_ct * allele_ct, diploid_pair_cts);
        diploid_pair_cts[0] = cur_raw_geno_cts[0];
        // lo_idx * allele_ct + hi_idx
        diploid_pair_cts[1] = cur_raw_geno_cts[1] - pgv.patch_01_ct;
        for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
          diploid_pair_cts[pgv.patch_01_vals[uii]] += 1;
        }
        const uintptr_t allele_ctp1 = allele_ct + 1;
        diploid_pair_cts[allele_ctp1] = cur_raw_geno_cts[2] - pgv.patch_10_ct;
        for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
          const uintptr_t lo_code = pgv.patch_10_vals[2 * uii];
          const uintptr_t hi_code = pgv.patch_10_vals[2 * uii + 1];
          diploid_pair_cts[lo_code * allele_ct + hi_code] += 1;
        }
        if (is_autosomal_diploid) {
          homref_ct = cur_raw_geno_cts[0];
          homalt1_ct = diploid_pair_cts[allele_ctp1];
          missing_ct = nobs_base - cur_raw_geno_cts[0] - cur_raw_geno_cts[1] - cur_raw_geno_cts[2];
          hap_cts[0] = 0;
          hap_cts[1] = 0;
        } else {
          if (is_x) {
            missing_ct = nobs_base - cur_raw_geno_cts[0] - cur_raw_geno_cts[1] - cur_raw_geno_cts[2];
            ZeroU32Arr(allele_ct, hap_cts);
            if (x_male_geno_cts) {
              STD_ARRAY_KREF(uint32_t, 3) cur_male_geno_cts = x_male_geno_cts[variant_uidx - x_start];
              hap_cts[0] = cur_male_geno_cts[0];
              diploid_pair_cts[0] -= cur_male_geno_cts[0];
              // possible todo: try making two disjoint PgrGetM calls instead.
              // don't expect that to be better, though.
              uintptr_t sample_widx = 0;
              uintptr_t cur_patch_bits = pgv.patch_01_set[0];
              uint32_t male_patch_01_ct = 0;
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t lowbit = BitIter1y(pgv.patch_01_set, &sample_widx, &cur_patch_bits);
                if (sex_male_collapsed[sample_widx] & lowbit) {
                  diploid_pair_cts[pgv.patch_01_vals[uii]] -= 1;
                  ++male_patch_01_ct;
                }
              }
              missing_ct += male_patch_01_ct;
              diploid_pair_cts[1] -= cur_male_geno_cts[1] - male_patch_01_ct;
              sample_widx = 0;
              cur_patch_bits = pgv.patch_10_set[0];
              uint32_t* hap_cts_offset1 = &(hap_cts[1]);
              uint32_t male_patch_10_ct = 0;
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t lowbit = BitIter1y(pgv.patch_10_set, &sample_widx, &cur_patch_bits);
                if (sex_male_collapsed[sample_widx] & lowbit) {
                  const uintptr_t lo_code = pgv.patch_10_vals[2 * uii];
                  const uintptr_t hi_code = pgv.patch_10_vals[2 * uii + 1];
                  diploid_pair_cts[lo_code * allele_ct + hi_code] -= 1;
                  ++male_patch_10_ct;
                  if (lo_code == hi_code) {
                    hap_cts_offset1[lo_code] += 1;
                  }
                }
              }
              missing_ct += male_patch_10_ct;
              hap_cts[1] = cur_male_geno_cts[2] - male_patch_10_ct;
              diploid_pair_cts[allele_ctp1] -= hap_cts[1];
              for (uint32_t aidx = 2; aidx != allele_ct; ++aidx) {
                // subtract male rarehoms
                missing_ct -= hap_cts[aidx];
              }
            }
            homref_ct = diploid_pair_cts[0];
            homalt1_ct = diploid_pair_cts[allele_ctp1];
          } else {
            // chrY or other pure-haploid; hethap treated as missing
            uint32_t nonmissing_ct = 0;
            for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
              const uint32_t cur_ct = diploid_pair_cts[aidx * allele_ctp1];
              hap_cts[aidx] = cur_ct;
              nonmissing_ct += cur_ct;
            }
            missing_ct = nobs_base - nonmissing_ct;
          }
        }
        if (homref_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(homref_ct, cswritep);
        }
        if (refalt1_col) {
          if (is_autosomal_diploid || is_x) {
            *cswritep++ = '\t';
            cswritep = u32toa(diploid_pair_cts[1], cswritep);
          } else {
            cswritep = strcpya_k(cswritep, "\t0");
          }
        }
        if (refalt_col) {
          *cswritep++ = '\t';
          if (is_autosomal_diploid || is_x) {
            for (uint32_t aidx = 1; aidx != allele_ct; ++aidx) {
              cswritep = u32toa_x(diploid_pair_cts[aidx], ',', cswritep);
            }
            --cswritep;
          } else {
            // repeat "0,"
            cswritep = u16setsa(cswritep, 0x2c30, allele_ct - 2);
            *cswritep++ = '0';
          }
        }
        if (homalt1_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(homalt1_ct, cswritep);
        }
        if (xy_col) {
          *cswritep++ = '\t';
          if (is_autosomal_diploid || is_x) {
            for (uint32_t aidx_hi = xy_col_altonly; aidx_hi != allele_ct; ++aidx_hi) {
              for (uint32_t aidx_lo = xy_col_altonly; aidx_lo <= aidx_hi; ++aidx_lo) {
                cswritep = u32toa_x(diploid_pair_cts[aidx_lo * allele_ct + aidx_hi], ',', cswritep);
              }
            }
            --cswritep;
          } else {
            const uint32_t triangle_base = allele_ct - xy_col_altonly;
            cswritep = u16setsa(cswritep, 0x2c30, (triangle_base * (triangle_base + 1) / 2) - 1);
            *cswritep++ = '0';
          }
        }
        if (hapref_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(hap_cts[0], cswritep);
        }
        if (hapalt1_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(hap_cts[1], cswritep);
        }
        if (hap_col) {
          *cswritep++ = '\t';
          if (!is_autosomal_diploid) {
            for (uint32_t aidx = hap_col_altonly; aidx != allele_ct; ++aidx) {
              cswritep = u32toa_x(hap_cts[aidx], ',', cswritep);
            }
            --cswritep;
          } else {
            cswritep = u16setsa(cswritep, 0x2c30, allele_ct - hap_col_altonly - 1);
            *cswritep++ = '0';
          }
        }
        if (numeq_col) {
          *cswritep++ = '\t';
          if (is_autosomal_diploid || is_x) {
            for (uint32_t aidx_hi = 0; aidx_hi != allele_ct; ++aidx_hi) {
              for (uint32_t aidx_lo = 0; aidx_lo <= aidx_hi; ++aidx_lo) {
                const uint32_t cur_ct = diploid_pair_cts[aidx_lo * allele_ct + aidx_hi];
                if (cur_ct) {
                  cswritep = u32toa_x(aidx_lo, '/', cswritep);
                  cswritep = u32toa_x(aidx_hi, '=', cswritep);
                  cswritep = u32toa_x(cur_ct, ',', cswritep);
                }
              }
            }
          }
          if (!is_autosomal_diploid) {
            for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
              const uint32_t cur_ct = hap_cts[aidx];
              if (cur_ct) {
                cswritep = u32toa_x(aidx, '=', cswritep);
                cswritep = u32toa_x(cur_ct, ',', cswritep);
              }
            }
          }
          if (missing_ct != nobs_base) {
            --cswritep;
          } else {
            *cswritep++ = '.';
          }
        }
      }
      if (missing_col) {
        *cswritep++ = '\t';
        cswritep = u32toa(missing_ct, cswritep);
      }
      if (nobs_col) {
        *cswritep++ = '\t';
        cswritep = u32toa(nobs_base - missing_ct, cswritep);
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto WriteGenoCounts_ret_WRITE_FAIL;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto WriteGenoCounts_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    logprintfww("--geno-counts%s: Genotype counts written to %s .\n", output_zst? " zs" : "", outname);
  }
  while (0) {
  WriteGenoCounts_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WriteGenoCounts_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  WriteGenoCounts_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 WriteGenoCounts_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr WriteMissingnessReports(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* sample_missing_hc_cts, const uint32_t* sample_missing_dosage_cts, const uint32_t* sample_hethap_cts, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint32_t* variant_missing_hc_cts, const uint32_t* variant_missing_dosage_cts, const uint32_t* variant_hethap_cts, uint32_t sample_ct, uint32_t male_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t variant_ct, uintptr_t max_allele_slen, uint32_t first_hap_uidx, MissingRptFlags missing_rpt_flags, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t output_zst = missing_rpt_flags & kfMissingRptZs;
    if (!(missing_rpt_flags & kfMissingRptVariantOnly)) {
      const uintptr_t overflow_buf_size = kCompressStreamBlock + kMaxMediumLine + pheno_ct * 2;
      OutnameZstSet(".smiss", output_zst, outname_end);
      reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
      if (unlikely(reterr)) {
        goto WriteMissingnessReports_ret_1;
      }
      *cswritep++ = '#';
      const uint32_t scol_fid = FidColIsRequired(siip, missing_rpt_flags / kfMissingRptScolMaybefid);
      if (scol_fid) {
        cswritep = strcpya_k(cswritep, "FID\t");
      }
      cswritep = strcpya_k(cswritep, "IID");
      const char* sample_ids = siip->sample_ids;
      const char* sids = siip->sids;
      const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
      const uintptr_t max_sid_blen = siip->max_sid_blen;
      const uint32_t scol_sid = SidColIsRequired(sids, missing_rpt_flags / kfMissingRptScolMaybesid);
      if (scol_sid) {
        cswritep = strcpya_k(cswritep, "\tSID");
      }
      const uint32_t scol_empty_pheno = (missing_rpt_flags & kfMissingRptScolMisspheno1) && (!pheno_ct);
      if (scol_empty_pheno) {
        cswritep = strcpya_k(cswritep, "\tMISS_PHENO1");
      }
      const uint32_t scol_phenos = (missing_rpt_flags & (kfMissingRptScolMisspheno1 | kfMissingRptScolMissphenos)) && pheno_ct;
      if (scol_phenos) {
        if (!(missing_rpt_flags & kfMissingRptScolMissphenos)) {
          pheno_ct = 1;
        }
        for (uintptr_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, &(pheno_names[pheno_idx * max_pheno_name_blen]));
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto WriteMissingnessReports_ret_WRITE_FAIL;
          }
        }
      }
      const uint32_t scol_nmiss_dosage = (missing_rpt_flags / kfMissingRptScolNmissDosage) & 1;
      if (scol_nmiss_dosage) {
        cswritep = strcpya_k(cswritep, "\tMISSING_DOSAGE_CT");
      }
      const uint32_t scol_nmiss = (missing_rpt_flags / kfMissingRptScolNmiss) & 1;
      if (scol_nmiss) {
        cswritep = strcpya_k(cswritep, "\tMISSING_CT");
      }
      const uint32_t scol_nmiss_hh = (missing_rpt_flags / kfMissingRptScolNmissHh) & 1;
      if (scol_nmiss_hh) {
        cswritep = strcpya_k(cswritep, "\tMISSING_AND_HETHAP_CT");
      }
      const uint32_t scol_hethap = (missing_rpt_flags / kfMissingRptScolHethap) & 1;
      if (scol_hethap) {
        cswritep = strcpya_k(cswritep, "\tHETHAP_CT");
      }
      const uint32_t scol_nobs = (missing_rpt_flags / kfMissingRptScolNobs) & 1;
      if (scol_nobs) {
        cswritep = strcpya_k(cswritep, "\tOBS_CT");
      }
      const uint32_t scol_fmiss_dosage = (missing_rpt_flags / kfMissingRptScolFmissDosage) & 1;
      if (scol_fmiss_dosage) {
        cswritep = strcpya_k(cswritep, "\tF_MISS_DOSAGE");
      }
      const uint32_t scol_fmiss = (missing_rpt_flags / kfMissingRptScolFmiss) & 1;
      if (scol_fmiss) {
        cswritep = strcpya_k(cswritep, "\tF_MISS");
      }
      const uint32_t scol_fmiss_hh = (missing_rpt_flags / kfMissingRptScolFmissHh) & 1;
      if (scol_fmiss_hh) {
        cswritep = strcpya_k(cswritep, "\tF_MISS_AND_HETHAP");
      }
      AppendBinaryEoln(&cswritep);
      uint32_t variant_ct_y = 0;
      uint32_t y_code;
      if (XymtExists(cip, kChrOffsetY, &y_code)) {
        variant_ct_y = CountChrVariantsUnsafe(variant_include, cip, y_code);
      }
      const uint32_t variant_ct_nony = variant_ct - variant_ct_y;
      char nobs_strs[2][16];
      uint32_t nobs_slens[2];
      double variant_ct_recips[2];
      char* write_iter = nobs_strs[0];
      *write_iter++ = '\t';
      write_iter = u32toa(variant_ct_nony, write_iter);
      nobs_slens[0] = write_iter - nobs_strs[0];
      variant_ct_recips[0] = 1.0 / u31tod(variant_ct_nony);
      write_iter = nobs_strs[1];
      *write_iter++ = '\t';
      write_iter = u32toa(variant_ct, write_iter);
      nobs_slens[1] = write_iter - nobs_strs[1];
      variant_ct_recips[1] = 1.0 / u31tod(variant_ct);
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
        const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
        if (!scol_fid) {
          cur_sample_id = AdvPastDelim(cur_sample_id, '\t');
        }
        cswritep = strcpya(cswritep, cur_sample_id);
        if (scol_sid) {
          *cswritep++ = '\t';
          if (sids) {
            cswritep = strcpya(cswritep, &(sids[sample_uidx * max_sid_blen]));
          } else {
            *cswritep++ = '0';
          }
        }
        if (scol_phenos) {
          for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
            *cswritep++ = '\t';
            // 'Y' - 'N' == 11
            *cswritep++ = 'Y' - 11 * IsSet(pheno_cols[pheno_idx].nonmiss, sample_uidx);
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto WriteMissingnessReports_ret_WRITE_FAIL;
            }
          }
        } else {
          if (scol_empty_pheno) {
            cswritep = strcpya_k(cswritep, "\tY");
          }
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto WriteMissingnessReports_ret_WRITE_FAIL;
          }
        }
        if (scol_nmiss_dosage) {
          *cswritep++ = '\t';
          cswritep = u32toa(sample_missing_dosage_cts[sample_uidx], cswritep);
        }
        const uint32_t cur_missing_hc_base = sample_missing_hc_cts[sample_uidx];
        if (scol_nmiss) {
          *cswritep++ = '\t';
          cswritep = u32toa(cur_missing_hc_base, cswritep);
        }
        if (scol_nmiss_hh) {
          *cswritep++ = '\t';
          cswritep = u32toa(cur_missing_hc_base + sample_hethap_cts[sample_uidx], cswritep);
        }
        if (scol_hethap) {
          *cswritep++ = '\t';
          cswritep = u32toa(sample_hethap_cts[sample_uidx], cswritep);
        }
        const uint32_t is_male = IsSet(sex_male, sample_uidx);
        if (scol_nobs) {
          cswritep = memcpya(cswritep, nobs_strs[is_male], nobs_slens[is_male]);
        }
        const double cur_variant_ct_recip = variant_ct_recips[is_male];
        if (scol_fmiss_dosage) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(u31tod(sample_missing_dosage_cts[sample_uidx]) * cur_variant_ct_recip, cswritep);
        }
        if (scol_fmiss) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(u31tod(cur_missing_hc_base) * cur_variant_ct_recip, cswritep);
        }
        if (scol_fmiss_hh) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(u31tod(cur_missing_hc_base + sample_hethap_cts[sample_uidx]) * cur_variant_ct_recip, cswritep);
        }
        AppendBinaryEoln(&cswritep);
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto WriteMissingnessReports_ret_WRITE_FAIL;
      }
      BigstackReset(bigstack_mark);
      logprintfww("--missing: Sample missing data report written to %s .\n", outname);
    }
    if (!(missing_rpt_flags & kfMissingRptSampleOnly)) {
      const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
      char* chr_buf;  // includes trailing tab
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto WriteMissingnessReports_ret_NOMEM;
      }
      const uintptr_t overflow_buf_size = kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + 2 * max_allele_slen;
      OutnameZstSet(".vmiss", output_zst, outname_end);
      reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
      if (unlikely(reterr)) {
        goto WriteMissingnessReports_ret_1;
      }
      *cswritep++ = '#';
      const uint32_t chr_col = missing_rpt_flags & kfMissingRptVcolChrom;

      if (chr_col) {
        cswritep = strcpya_k(cswritep, "CHROM\t");
      }
      if (missing_rpt_flags & kfMissingRptVcolPos) {
        cswritep = strcpya_k(cswritep, "POS\t");
      } else {
        variant_bps = nullptr;
      }
      cswritep = strcpya_k(cswritep, "ID");
      const uint32_t ref_col = missing_rpt_flags & kfMissingRptVcolRef;
      if (ref_col) {
        cswritep = strcpya_k(cswritep, "\tREF");
      }
      const uint32_t alt1_col = missing_rpt_flags & kfMissingRptVcolAlt1;
      if (alt1_col) {
        cswritep = strcpya_k(cswritep, "\tALT1");
      }
      const uint32_t alt_col = missing_rpt_flags & kfMissingRptVcolAlt;
      if (alt_col) {
        cswritep = strcpya_k(cswritep, "\tALT");
      }
      const uint32_t nmiss_dosage_col = missing_rpt_flags & kfMissingRptVcolNmissDosage;
      if (nmiss_dosage_col) {
        cswritep = strcpya_k(cswritep, "\tMISSING_DOSAGE_CT");
      }
      const uint32_t nmiss_col = (missing_rpt_flags / kfMissingRptVcolNmiss) & 1;
      if (nmiss_col) {
        cswritep = strcpya_k(cswritep, "\tMISSING_CT");
      }
      const uint32_t nmiss_hh_col = (missing_rpt_flags / kfMissingRptVcolNmissHh) & 1;
      if (nmiss_hh_col) {
        cswritep = strcpya_k(cswritep, "\tMISSING_AND_HETHAP_CT");
      }
      const uint32_t hethap_col = (missing_rpt_flags / kfMissingRptVcolHethap) & 1;
      if (hethap_col) {
        cswritep = strcpya_k(cswritep, "\tHETHAP_CT");
      }
      const uint32_t nobs_col = (missing_rpt_flags / kfMissingRptVcolNobs) & 1;
      if (nobs_col) {
        cswritep = strcpya_k(cswritep, "\tOBS_CT");
      }
      const uint32_t fmiss_dosage_col = missing_rpt_flags & kfMissingRptVcolFmissDosage;
      if (fmiss_dosage_col) {
        cswritep = strcpya_k(cswritep, "\tF_MISS_DOSAGE");
      }
      const uint32_t fmiss_col = (missing_rpt_flags / kfMissingRptVcolFmiss) & 1;
      if (fmiss_col) {
        cswritep = strcpya_k(cswritep, "\tF_MISS");
      }
      const uint32_t fmiss_hh_col = (missing_rpt_flags / kfMissingRptVcolFmissHh) & 1;
      if (fmiss_hh_col) {
        cswritep = strcpya_k(cswritep, "\tF_MISS_AND_HETHAP");
      }
      const uint32_t fhethap_col = (missing_rpt_flags / kfMissingRptVcolFhethap) & 1;
      if (fhethap_col) {
        cswritep = strcpya_k(cswritep, "\tF_HETHAP");
      }
      AppendBinaryEoln(&cswritep);
      char nobs_str[16];
      nobs_str[0] = '\t';
      const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      uint32_t nobs_slen = 0;
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uint32_t is_y = 2;
      double nobs_recip = 0.0;
      fputs("--missing variant report: 0%", stdout);
      fflush(stdout);
      uint32_t cur_allele_ct = 2;
      uint32_t cur_missing_hc_ct = 0;
      uint32_t cur_hethap_ct = 0;
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        if (variant_uidx >= chr_end) {
          do {
            ++chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          } while (variant_uidx >= chr_end);
          const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
          char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
          *chr_name_end = '\t';
          chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
          const uint32_t new_is_y = (chr_idx == y_code);
          if (new_is_y != is_y) {
            is_y = new_is_y;
            const uint32_t cur_nobs = is_y? male_ct : sample_ct;
            nobs_recip = 1.0 / u31tod(cur_nobs);
            char* nobs_str_end = u32toa(cur_nobs, &(nobs_str[1]));
            nobs_slen = nobs_str_end - nobs_str;
          }
        }
        if (chr_col) {
          cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
        }
        if (variant_bps) {
          cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
        }
        cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        if (ref_col) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, cur_alleles[0]);
        }
        if (alt1_col) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, cur_alleles[1]);
        }
        if (alt_col) {
          *cswritep++ = '\t';
          for (uint32_t allele_idx = 1; allele_idx != cur_allele_ct; ++allele_idx) {
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto WriteMissingnessReports_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
          }
          --cswritep;
        }
        if (nmiss_dosage_col) {
          *cswritep++ = '\t';
          cswritep = u32toa(variant_missing_dosage_cts[variant_uidx], cswritep);
        }
        if (variant_missing_hc_cts) {
          cur_missing_hc_ct = variant_missing_hc_cts[variant_uidx];
          cur_hethap_ct = 0;
          if (variant_uidx >= first_hap_uidx) {
            cur_hethap_ct = variant_hethap_cts[variant_uidx - first_hap_uidx];
          }
          if (nmiss_col) {
            *cswritep++ = '\t';
            cswritep = u32toa(cur_missing_hc_ct, cswritep);
          }
          if (nmiss_hh_col) {
            *cswritep++ = '\t';
            cswritep = u32toa(cur_missing_hc_ct + cur_hethap_ct, cswritep);
          }
          if (hethap_col) {
            *cswritep++ = '\t';
            cswritep = u32toa(cur_hethap_ct, cswritep);
          }
        }
        if (nobs_col) {
          cswritep = memcpya(cswritep, nobs_str, nobs_slen);
        }
        if (fmiss_dosage_col) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(u31tod(variant_missing_dosage_cts[variant_uidx]) * nobs_recip, cswritep);
        }
        if (fmiss_col) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(u31tod(cur_missing_hc_ct) * nobs_recip, cswritep);
        }
        if (fmiss_hh_col) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(u31tod(cur_missing_hc_ct + cur_hethap_ct) * nobs_recip, cswritep);
        }
        if (fhethap_col) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(u31tod(cur_hethap_ct) * nobs_recip, cswritep);
        }
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto WriteMissingnessReports_ret_WRITE_FAIL;
        }
        if (variant_idx >= next_print_variant_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (variant_idx * 100LLU) / variant_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
        }
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto WriteMissingnessReports_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      logprintfww("--missing: Variant missing data report written to %s .\n", outname);
    }
  }
  while (0) {
  WriteMissingnessReports_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WriteMissingnessReports_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 WriteMissingnessReports_ret_1:
  BigstackReset(bigstack_mark);
  CswriteCloseCond(&css, cswritep);
  return reterr;
}

PglErr GetMultiallelicMarginalCounts(const uintptr_t* founder_info, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_geno_cts), uint32_t raw_sample_ct, uint32_t autosomal_variant_ct, uint32_t autosomal_xallele_ct, uint32_t hwe_x_ct, uint32_t x_xallele_ct, PgenReader* simple_pgrp, STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts), STD_ARRAY_PTR_DECL(uint32_t, 2, autosomal_xgeno_cts), STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts)) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uint32_t* cumulative_popcounts;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &cumulative_popcounts))) {
      goto GetMultiallelicMarginalCounts_ret_NOMEM;
    }
    FillCumulativePopcounts(founder_info, raw_sample_ctl, cumulative_popcounts);
    const uint32_t founder_ct = cumulative_popcounts[raw_sample_ctl - 1] + PopcountWord(founder_info[raw_sample_ctl - 1]);
    const uint32_t founder_ctl2 = NypCtToWordCt(founder_ct);
    const uint32_t founder_ctl = BitCtToWordCt(founder_ct);
    const uint32_t max_allele_ct = PgrGetMaxAlleleCt(simple_pgrp);
    PgenVariant pgv;
    uint32_t* one_cts;
    uint32_t* two_cts;
    if (unlikely(
            bigstack_alloc_w(founder_ctl2, &(pgv.genovec)) ||
            bigstack_alloc_w(founder_ctl, &(pgv.patch_01_set)) ||
            bigstack_alloc_ac(founder_ct, &(pgv.patch_01_vals)) ||
            bigstack_alloc_w(founder_ctl, &(pgv.patch_10_set)) ||
            bigstack_alloc_ac(2 * founder_ct, &(pgv.patch_10_vals)) ||
            bigstack_alloc_u32(max_allele_ct, &one_cts) ||
            bigstack_alloc_u32(max_allele_ct, &two_cts))) {
      goto GetMultiallelicMarginalCounts_ret_NOMEM;
    }
    if (autosomal_xallele_ct) {
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cumulative_popcounts, simple_pgrp, &pssi);
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      uintptr_t xgeno_idx = 0;
      for (uint32_t variant_idx = 0; variant_idx != autosomal_variant_ct; ++variant_idx) {
        uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        if (variant_uidx >= chr_end) {
          uint32_t chr_idx;
          do {
            ++chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            chr_idx = cip->chr_file_order[chr_fo_idx];
          } while ((variant_uidx >= chr_end) || IsSet(cip->haploid_mask, chr_idx));
          BitIter1Start(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], &variant_uidx_base, &cur_bits);
          variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        }
        const uint32_t allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
        if (allele_ct > 2) {
          reterr = PgrGetM(founder_info, pssi, founder_ct, variant_uidx, simple_pgrp, &pgv);
          if (unlikely(reterr)) {
            goto GetMultiallelicMarginalCounts_ret_PGR_FAIL;
          }
          ZeroTrailingNyps(founder_ct, pgv.genovec);
          ZeroU32Arr(allele_ct, one_cts);
          ZeroU32Arr(allele_ct, two_cts);
          // const uint32_t hom_ref_ct = hwe_geno_cts[variant_uidx][0];
          const uint32_t het_ref_ct = hwe_geno_cts[variant_uidx][1];
          const uint32_t altxy_ct = hwe_geno_cts[variant_uidx][2];
          // two_cts[0] = hom_ref_ct;
          // one_cts[0] = het_ref_ct;
          one_cts[1] = het_ref_ct - pgv.patch_01_ct;
          two_cts[1] = altxy_ct - pgv.patch_10_ct;
          for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
            one_cts[pgv.patch_01_vals[uii]] += 1;
          }
          for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
            const uintptr_t lo_code = pgv.patch_10_vals[2 * uii];
            const uintptr_t hi_code = pgv.patch_10_vals[2 * uii + 1];
            if (lo_code == hi_code) {
              two_cts[lo_code] += 1;
            } else {
              one_cts[lo_code] += 1;
              one_cts[hi_code] += 1;
            }
          }
          for (uint32_t aidx = 1; aidx != allele_ct; ++aidx) {
            // This is currently inverted, to match how hwe_geno_cts represents
            // REF allele counts.
            autosomal_xgeno_cts[xgeno_idx][0] = two_cts[aidx];
            autosomal_xgeno_cts[xgeno_idx][1] = one_cts[aidx];
            ++xgeno_idx;
          }
        }
      }
    }
    if (x_xallele_ct) {
      // Related to multiallelic-chrX --geno-counts implementation.
      uintptr_t* founder_knownsex;
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &founder_knownsex))) {
        goto GetMultiallelicMarginalCounts_ret_NOMEM;
      }
      BitvecAndCopy(founder_info, sex_nm, raw_sample_ctl, founder_knownsex);
      FillCumulativePopcounts(founder_knownsex, raw_sample_ctl, cumulative_popcounts);
      const uint32_t founder_x_ct = cumulative_popcounts[raw_sample_ctl - 1] + PopcountWord(founder_knownsex[raw_sample_ctl - 1]);
      const uint32_t founder_x_ctaw = BitCtToAlignedWordCt(founder_x_ct);
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cumulative_popcounts, simple_pgrp, &pssi);
      uintptr_t* founder_male_collapsed;
      uintptr_t* founder_male_interleaved_vec;
      if (unlikely(
              bigstack_alloc_w(founder_x_ctaw, &founder_male_collapsed) ||
              bigstack_alloc_w(founder_x_ctaw, &founder_male_interleaved_vec))) {
        goto GetMultiallelicMarginalCounts_ret_NOMEM;
      }
      CopyBitarrSubset(sex_male, founder_knownsex, founder_x_ct, founder_male_collapsed);
      const uint32_t founder_x_ctl = BitCtToWordCt(founder_x_ct);
      const uint32_t founder_male_ct = PopcountWords(founder_male_collapsed, founder_ctl);
      ZeroWArr(founder_x_ctaw - founder_x_ctl, &(founder_male_collapsed[founder_x_ctl]));
      FillInterleavedMaskVec(founder_male_collapsed, founder_x_ctaw / kWordsPerVec, founder_male_interleaved_vec);
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
      const uint32_t x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
      uintptr_t variant_uidx_base;
      uintptr_t cur_bits;
      BitIter1Start(variant_include, x_start, &variant_uidx_base, &cur_bits);
      uintptr_t xgeno_idx = 0;
      for (uint32_t variant_idx = 0; variant_idx != hwe_x_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        const uint32_t allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
        if (allele_ct > 2) {
          reterr = PgrGetM(founder_knownsex, pssi, founder_x_ct, variant_uidx, simple_pgrp, &pgv);
          if (unlikely(reterr)) {
            goto GetMultiallelicMarginalCounts_ret_PGR_FAIL;
          }
          ZeroTrailingNyps(founder_x_ct, pgv.genovec);
          ZeroU32Arr(allele_ct, one_cts);
          ZeroU32Arr(allele_ct, two_cts);
          STD_ARRAY_DECL(uint32_t, 4, genocounts);
          GenoarrCountFreqsUnsafe(pgv.genovec, founder_x_ct, genocounts);
          // const uint32_t hom_ref_ct = genocounts[0];
          const uint32_t het_ref_ct = genocounts[1];
          const uint32_t altxy_ct = genocounts[2];
          // two_cts[0] = hom_ref_ct;
          // one_cts[0] = het_ref_ct;
          one_cts[1] = het_ref_ct - pgv.patch_01_ct;
          two_cts[1] = altxy_ct - pgv.patch_10_ct;
          for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
            one_cts[pgv.patch_01_vals[uii]] += 1;
          }
          for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
            const uintptr_t lo_code = pgv.patch_10_vals[2 * uii];
            const uintptr_t hi_code = pgv.patch_10_vals[2 * uii + 1];
            if (lo_code == hi_code) {
              two_cts[lo_code] += 1;
            } else {
              one_cts[lo_code] += 1;
              one_cts[hi_code] += 1;
            }
          }
          uintptr_t knownsex_xgeno_idx = xgeno_idx;
          for (uint32_t aidx = 1; aidx != allele_ct; ++aidx) {
            x_knownsex_xgeno_cts[knownsex_xgeno_idx][0] = two_cts[aidx];
            x_knownsex_xgeno_cts[knownsex_xgeno_idx][1] = one_cts[aidx];
            ++knownsex_xgeno_idx;
          }
          if (x_male_xgeno_cts) {
            ZeroU32Arr(allele_ct, one_cts);
            ZeroU32Arr(allele_ct, two_cts);
            GenoarrCountSubsetFreqs(pgv.genovec, founder_male_interleaved_vec, founder_x_ct, founder_male_ct, genocounts);
            const uint32_t male_het_ref_ct = genocounts[1];
            const uint32_t male_altxy_ct = genocounts[2];
            one_cts[1] = male_het_ref_ct;
            two_cts[1] = male_altxy_ct;
            uintptr_t sample_widx = 0;
            uintptr_t cur_patch_bits = pgv.patch_01_set[0];
            uint32_t male_patch_01_ct = 0;
            for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
              const uintptr_t lowbit = BitIter1y(pgv.patch_01_set, &sample_widx, &cur_patch_bits);
              if (founder_male_collapsed[sample_widx] & lowbit) {
                ++male_patch_01_ct;
                one_cts[pgv.patch_01_vals[uii]] += 1;
              }
            }
            one_cts[1] -= male_patch_01_ct;
            sample_widx = 0;
            cur_patch_bits = pgv.patch_10_set[0];
            uint32_t male_patch_10_ct = 0;
            for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
              const uintptr_t lowbit = BitIter1y(pgv.patch_10_set, &sample_widx, &cur_patch_bits);
              if (founder_male_collapsed[sample_widx] & lowbit) {
                ++male_patch_10_ct;
                const uintptr_t lo_code = pgv.patch_10_vals[2 * uii];
                const uintptr_t hi_code = pgv.patch_10_vals[2 * uii + 1];
                if (lo_code == hi_code) {
                  two_cts[lo_code] += 1;
                } else {
                  one_cts[lo_code] += 1;
                  one_cts[hi_code] += 1;
                }
              }
            }
            two_cts[1] -= male_patch_10_ct;
            for (uint32_t aidx = 1; aidx != allele_ct; ++aidx) {
              x_male_xgeno_cts[xgeno_idx][0] = two_cts[aidx];
              x_male_xgeno_cts[xgeno_idx][1] = one_cts[aidx];
              ++xgeno_idx;
            }
          } else {
            xgeno_idx = knownsex_xgeno_idx;
          }
        }
      }
    }
  }
  while (0) {
  GetMultiallelicMarginalCounts_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GetMultiallelicMarginalCounts_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct ComputeHweXPvalsCtxStruct {
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts);
  const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts);
  const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts);
  const STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts);
  const STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts);
  uint32_t* variant_uidx_starts;
  uintptr_t* extra_aidx_starts;
  uint32_t x_start;
  uint32_t hwe_x_start;
  uint32_t hwe_midp;
  uint32_t hwe_x_ct;

  double* hwe_x_pvals;
} ComputeHweXPvalsCtx;

void ComputeHweXPvalsMain(uintptr_t tidx, uintptr_t thread_ct, ComputeHweXPvalsCtx* ctx) {
  const uintptr_t* variant_include = ctx->variant_include;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts) = ctx->founder_raw_geno_cts;
  const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts) = ctx->founder_x_male_geno_cts;
  const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts) = ctx->founder_x_nosex_geno_cts;
  const STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts) = ctx->x_knownsex_xgeno_cts;
  const STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts) = ctx->x_male_xgeno_cts;
  const uint32_t x_start = ctx->x_start;
  const uint32_t hwe_x_ct = ctx->hwe_x_ct;
  const uint32_t hwe_midp = ctx->hwe_midp;

  // this needs to be aligned with ComputeUidxStartPartition()
  const uint32_t variant_idx_end = (hwe_x_ct * (S_CAST(uint64_t, tidx) + 1)) / thread_ct;
  uint32_t variant_idx = (hwe_x_ct * S_CAST(uint64_t, tidx)) / thread_ct;

  uintptr_t xgeno_idx = ctx->extra_aidx_starts[tidx];
  double* hwe_x_pvals_iter = &(ctx->hwe_x_pvals[variant_idx + xgeno_idx]);
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, ctx->variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
  uint32_t pct = 0;
  uint32_t next_print_variant_idx = variant_idx_end;
  if (!tidx) {
    next_print_variant_idx = variant_idx_end / 100;
  }
  uint32_t male_1copy_ct = 0;
  uint32_t male_hethap_ct = 0;
  uint32_t male_0copy_ct = 0;
  for (; variant_idx != variant_idx_end; ++variant_idx) {
    const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    STD_ARRAY_KREF(uint32_t, 3) cur_raw_geno_cts = founder_raw_geno_cts[variant_uidx];
    uint32_t female_2copy_ct = cur_raw_geno_cts[0];
    uint32_t female_1copy_ct = cur_raw_geno_cts[1];
    uint32_t female_0copy_ct = cur_raw_geno_cts[2];
    if (founder_x_male_geno_cts) {
      STD_ARRAY_KREF(uint32_t, 3) cur_male_geno_cts = founder_x_male_geno_cts[variant_uidx - x_start];
      male_1copy_ct = cur_male_geno_cts[0];
      female_2copy_ct -= male_1copy_ct;
      male_hethap_ct = cur_male_geno_cts[1];
      female_1copy_ct -= male_hethap_ct;
      male_0copy_ct = cur_male_geno_cts[2];
      female_0copy_ct -= male_0copy_ct;
    }
    if (founder_x_nosex_geno_cts) {
      STD_ARRAY_KREF(uint32_t, 3) cur_nosex_geno_cts = founder_x_nosex_geno_cts[variant_uidx - x_start];
      female_2copy_ct -= cur_nosex_geno_cts[0];
      female_1copy_ct -= cur_nosex_geno_cts[1];
      female_0copy_ct -= cur_nosex_geno_cts[2];
    }
    *hwe_x_pvals_iter++ = HweXchrP(female_1copy_ct, female_2copy_ct, female_0copy_ct, male_1copy_ct, male_0copy_ct, hwe_midp);
    if (allele_idx_offsets) {
      const uint32_t allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
      if (allele_ct != 2) {
        const uint32_t female_obs_ct = female_2copy_ct + female_1copy_ct + female_0copy_ct;
        const uint32_t male_obs_ct = male_1copy_ct + male_hethap_ct + male_0copy_ct;
        for (uint32_t aidx = 1; aidx != allele_ct; ++aidx) {
          female_2copy_ct = x_knownsex_xgeno_cts[xgeno_idx][0];
          female_1copy_ct = x_knownsex_xgeno_cts[xgeno_idx][1];
          if (x_male_xgeno_cts) {
            male_1copy_ct = x_male_xgeno_cts[xgeno_idx][0];
            female_2copy_ct -= male_1copy_ct;
            male_hethap_ct = x_male_xgeno_cts[xgeno_idx][1];
            female_1copy_ct -= male_hethap_ct;
            male_0copy_ct = male_obs_ct - male_1copy_ct - male_hethap_ct;
          }
          female_0copy_ct = female_obs_ct - female_2copy_ct - female_1copy_ct;
          *hwe_x_pvals_iter++ = HweXchrP(female_1copy_ct, female_2copy_ct, female_0copy_ct, male_1copy_ct, male_0copy_ct, hwe_midp);
          ++xgeno_idx;
        }
      }
    }
    if (variant_idx >= next_print_variant_idx) {
      // only possible for tidx == 0
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      pct = (variant_idx * 100LLU) / variant_idx_end;
      printf("\b\b%u%%", pct++);
      fflush(stdout);
      next_print_variant_idx = (pct * S_CAST(uint64_t, variant_idx_end)) / 100;
    }
  }
  if (pct > 10) {
    putc_unlocked('\b', stdout);
  }
}

THREAD_FUNC_DECL ComputeHweXPvalsThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  ComputeHweXPvalsCtx* ctx = S_CAST(ComputeHweXPvalsCtx*, arg->sharedp->context);
  ComputeHweXPvalsMain(arg->tidx, GetThreadCt(arg->sharedp) + 1, ctx);
  THREAD_RETURN;
}

PglErr ComputeHweXPvals(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts), uint32_t x_start, uint32_t hwe_x_ct, uintptr_t x_xallele_ct, uint32_t hwe_midp, uint32_t calc_thread_ct, double** hwe_x_pvals_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  ComputeHweXPvalsCtx ctx;
  {
    assert(hwe_x_ct);
    if (unlikely(bigstack_alloc_d(hwe_x_ct + x_xallele_ct, hwe_x_pvals_ptr))) {
      goto ComputeHweXPvals_ret_NOMEM;
    }
    bigstack_mark = g_bigstack_base;
    ctx.hwe_x_pvals = *hwe_x_pvals_ptr;

    if (calc_thread_ct > hwe_x_ct) {
      calc_thread_ct = hwe_x_ct;
    }
    if (unlikely(
            SetThreadCt0(calc_thread_ct - 1, &tg) ||
            bigstack_alloc_u32(calc_thread_ct, &ctx.variant_uidx_starts) ||
            bigstack_alloc_w(calc_thread_ct, &ctx.extra_aidx_starts))) {
      goto ComputeHweXPvals_ret_NOMEM;
    }
    // possible todo: extra-allele-based load balancer
    ComputeUidxStartPartition(variant_include, hwe_x_ct, calc_thread_ct, x_start, ctx.variant_uidx_starts);
    ctx.extra_aidx_starts[0] = 0;
    uintptr_t extra_aidx = 0;
    uint32_t prev_variant_uidx = ctx.variant_uidx_starts[0];
    for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
      const uint32_t cur_variant_uidx = ctx.variant_uidx_starts[tidx];
      extra_aidx += CountExtraAlleles(variant_include, allele_idx_offsets, prev_variant_uidx, cur_variant_uidx, 1);
      ctx.extra_aidx_starts[tidx] = extra_aidx;
      prev_variant_uidx = cur_variant_uidx;
    }
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.founder_raw_geno_cts = founder_raw_geno_cts;
    ctx.founder_x_male_geno_cts = founder_x_male_geno_cts;
    ctx.founder_x_nosex_geno_cts = founder_x_nosex_geno_cts;
    ctx.x_knownsex_xgeno_cts = x_knownsex_xgeno_cts;
    ctx.x_male_xgeno_cts = x_male_xgeno_cts;
    ctx.x_start = x_start;
    ctx.hwe_x_ct = hwe_x_ct;
    ctx.hwe_midp = hwe_midp;
    logprintf("Computing chrX Hardy-Weinberg %sp-values... ", hwe_midp? "mid" : "");
    fputs("0%", stdout);
    fflush(stdout);
    if (calc_thread_ct > 1) {
      SetThreadFuncAndData(ComputeHweXPvalsThread, &ctx, &tg);
      DeclareLastThreadBlock(&tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto ComputeHweXPvals_ret_THREAD_CREATE_FAIL;
      }
    }
    ComputeHweXPvalsMain(calc_thread_ct - 1, calc_thread_ct, &ctx);
    JoinThreads0(&tg);
    fputs("\b\b", stdout);
    logputs("done.\n");
  }
  while (0) {
  ComputeHweXPvals_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ComputeHweXPvals_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr HardyReport(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, autosomal_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_nosex_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts), const double* hwe_x_pvals, uint32_t variant_ct, uint32_t hwe_x_ct, uint32_t max_allele_slen, double output_min_ln, HardyFlags hardy_flags, uint32_t max_thread_ct, uint32_t nonfounders, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
    const uint32_t chr_code_endl = BitCtToWordCt(chr_code_end);
    const uintptr_t overflow_buf_size = RoundUpPow2(kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + 2 * max_allele_slen, kCacheline);
    const uint32_t output_zst = hardy_flags & kfHardyZs;
    const double output_min_p = (output_min_ln < kLnNormalMin)? 0 : exp(output_min_ln);
    uintptr_t overflow_buf_alloc = overflow_buf_size;
    if (output_zst) {
      overflow_buf_alloc += CstreamWkspaceReq(overflow_buf_size);
    }
    char* overflow_buf;
    uintptr_t* chr_skips;
    if (unlikely(
            bigstack_alloc_c(overflow_buf_alloc, &overflow_buf) ||
            bigstack_alloc_w(chr_code_endl, &chr_skips))) {
      goto HardyReport_ret_NOMEM;
    }
    // skip chrX, chrY, chrM here
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    memcpy(chr_skips, cip->haploid_mask, chr_code_endl * sizeof(intptr_t));
    if (!IsI32Neg(mt_code)) {
      SetBit(mt_code, chr_skips);
    }
    const uint32_t chr_skip_ct = PopcountWords(chr_skips, chr_code_endl);
    uint32_t variant_skip_ct = 0;
    uintptr_t chr_uidx_base = 0;
    uintptr_t chr_skips_bits = chr_skips[0];
    for (uint32_t chr_skip_idx = 0; chr_skip_idx != chr_skip_ct; ++chr_skip_idx) {
      const uintptr_t chr_uidx = BitIter1(chr_skips, &chr_uidx_base, &chr_skips_bits);
      if (IsSet(cip->chr_mask, chr_uidx)) {
        const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_uidx];
        variant_skip_ct += PopcountBitRange(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1]);
      }
    }
    if (variant_skip_ct - hwe_x_ct) {
      logprintf("--hardy: Skipping %u haploid variant%s.\n", variant_skip_ct - hwe_x_ct, (variant_skip_ct - hwe_x_ct == 1)? "" : "s");
    }
    variant_ct -= variant_skip_ct;
    const uint32_t midp = (hardy_flags / kfHardyMidp) & 1;
    const uint32_t redundant = (hardy_flags / kfHardyRedundant) & 1;
    const uint32_t chr_col = hardy_flags & kfHardyColChrom;
    const uint32_t ref_col = hardy_flags & kfHardyColRef;
    const uint32_t alt1_col = hardy_flags & kfHardyColAlt1;
    const uint32_t alt_col = hardy_flags & kfHardyColAlt;
    const uint32_t ax_col = hardy_flags & kfHardyColAx;
    const uint32_t gcounts = hardy_flags & (kfHardyColGcounts | kfHardyColGcount1col);
    const uint32_t gcount_1col = hardy_flags & kfHardyColGcount1col;
    const char gcount_delim = gcount_1col? ',' : '\t';
    const uint32_t hetfreq_cols = hardy_flags & kfHardyColHetfreq;
    const uint32_t p_col = hardy_flags & kfHardyColP;
    if (variant_ct) {
      OutnameZstSet(".hardy", output_zst, outname_end);
      reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
      if (unlikely(reterr)) {
        goto HardyReport_ret_1;
      }
      cswritep = overflow_buf;
      *cswritep++ = '#';

      // includes trailing tab
      char* chr_buf = nullptr;
      if (chr_col) {
        if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
          goto HardyReport_ret_NOMEM;
        }
        cswritep = strcpya_k(cswritep, "CHROM\t");
      }
      if (hardy_flags & kfHardyColPos) {
        cswritep = strcpya_k(cswritep, "POS\t");
      } else {
        variant_bps = nullptr;
      }
      cswritep = strcpya_k(cswritep, "ID");
      if (ref_col) {
        cswritep = strcpya_k(cswritep, "\tREF");
      }
      if (alt1_col) {
        cswritep = strcpya_k(cswritep, "\tALT1");
      }
      if (alt_col) {
        cswritep = strcpya_k(cswritep, "\tALT");
      }
      cswritep = strcpya_k(cswritep, "\tA1");
      if (ax_col) {
        cswritep = strcpya_k(cswritep, "\tAX");
      }
      if (gcounts) {
        if (gcount_1col) {
          cswritep = strcpya_k(cswritep, "\tGCOUNTS");
        } else {
          cswritep = strcpya_k(cswritep, "\tHOM_A1_CT\tHET_A1_CT\tTWO_AX_CT");
        }
      }
      if (hetfreq_cols) {
        cswritep = strcpya_k(cswritep, "\tO(HET_A1)\tE(HET_A1)");
      }
      if (p_col) {
        *cswritep++ = '\t';
        if (midp) {
          cswritep = strcpya_k(cswritep, "MIDP");
        } else {
          *cswritep++ = 'P';
        }
      }
      AppendBinaryEoln(&cswritep);
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      printf("--hardy%s%s: 0%%", output_zst? " zs" : "", midp? " midp" : "");
      fflush(stdout);
      uintptr_t xgeno_idx = 0;
      uint32_t allele_ct = 2;
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        // bugfix (15 May 2018): this needs to happen even if we aren't
        // printing #CHROM column
        if (variant_uidx >= chr_end) {
          uint32_t chr_idx;
          do {
            ++chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            chr_idx = cip->chr_file_order[chr_fo_idx];
          } while ((variant_uidx >= chr_end) || IsSet(chr_skips, chr_idx));
          BitIter1Start(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], &variant_uidx_base, &variant_include_bits);
          variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
          if (chr_col) {
            char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
            *chr_name_end = '\t';
            chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
          }
        }
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        const uint32_t print_allele_ct = ((allele_ct == 2) && (!redundant))? 1 : allele_ct;
        STD_ARRAY_KREF(uint32_t, 3) cur_geno_cts = hwe_geno_cts[variant_uidx];
        uint32_t het_a1_ct = cur_geno_cts[1];
        const uint32_t nonmissing_ct = cur_geno_cts[0] + het_a1_ct + cur_geno_cts[2];
        for (uint32_t a1_idx = 0; a1_idx != print_allele_ct; ++a1_idx) {
          if (chr_col) {
            cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
          }
          if (variant_bps) {
            cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
          }
          cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
          if (ref_col) {
            *cswritep++ = '\t';
            cswritep = strcpya(cswritep, cur_alleles[0]);
          }
          if (alt1_col) {
            *cswritep++ = '\t';
            cswritep = strcpya(cswritep, cur_alleles[1]);
          }
          if (alt_col) {
            *cswritep++ = '\t';
            for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto HardyReport_ret_WRITE_FAIL;
              }
              cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
            }
            --cswritep;
          }
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, cur_alleles[a1_idx]);
          if (ax_col) {
            *cswritep++ = '\t';
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if (allele_idx == a1_idx) {
                continue;
              }
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto HardyReport_ret_WRITE_FAIL;
              }
              cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
            }
            --cswritep;
          }
          uint32_t hom_a1_ct;
          uint32_t two_ax_ct;
          if (!a1_idx) {
            hom_a1_ct = cur_geno_cts[0];
            two_ax_ct = cur_geno_cts[2];
          } else if (allele_ct == 2) {
            // special case, don't read from autosomal_xgeno_cts
            hom_a1_ct = cur_geno_cts[2];
            two_ax_ct = cur_geno_cts[0];
          } else {
            STD_ARRAY_KREF(uint32_t, 2) xgeno_cts = autosomal_xgeno_cts[xgeno_idx];
            hom_a1_ct = xgeno_cts[0];
            het_a1_ct = xgeno_cts[1];
            two_ax_ct = nonmissing_ct - hom_a1_ct - het_a1_ct;
            ++xgeno_idx;
          }
          if (gcounts) {
            *cswritep++ = '\t';
            cswritep = u32toa_x(hom_a1_ct, gcount_delim, cswritep);
            cswritep = u32toa_x(het_a1_ct, gcount_delim, cswritep);
            cswritep = u32toa(two_ax_ct, cswritep);
          }
          if (hetfreq_cols) {
            *cswritep++ = '\t';
            const double nonmissing_ct_recip = 1.0 / u31tod(nonmissing_ct);
            cswritep = dtoa_g(u31tod(het_a1_ct) * nonmissing_ct_recip, cswritep);
            *cswritep++ = '\t';
            const double dbl_maj_freq = (hom_a1_ct * 2 + het_a1_ct) * nonmissing_ct_recip;
            // (1.0 - maj_freq) is vulnerable to catastrophic cancellation when
            // maj_freq is 1
            if (hom_a1_ct == nonmissing_ct) {
              *cswritep++ = '0';
            } else {
              const double expected_het_freq = dbl_maj_freq * (1.0 - dbl_maj_freq * 0.5);
              cswritep = dtoa_g(expected_het_freq, cswritep);
            }
          }
          if (p_col) {
            // possible todo: multithread this
            *cswritep++ = '\t';
            const double hwe_p = HweP(het_a1_ct, hom_a1_ct, two_ax_ct, midp);
            cswritep = dtoa_g(MAXV(hwe_p, output_min_p), cswritep);
          }
          AppendBinaryEoln(&cswritep);
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto HardyReport_ret_WRITE_FAIL;
          }
        }
        if (variant_idx >= next_print_variant_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (variant_idx * 100LLU) / variant_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
        }
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto HardyReport_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      logprintfww("--hardy%s%s: Autosomal Hardy-Weinberg report (%s) written to %s .\n", output_zst? " zs" : "", midp? " midp" : "", nonfounders? "all samples" : "founders only", outname);
    }
    if (hwe_x_ct) {
      BigstackReset(chr_skips);
      OutnameZstSet(".hardy.x", output_zst, outname_end);
      reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
      if (unlikely(reterr)) {
        goto HardyReport_ret_1;
      }
      cswritep = overflow_buf;
      *cswritep++ = '#';

      // includes trailing tab
      char x_name_buf[8];
      uint32_t x_name_blen = 0;
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      if (chr_col) {
        cswritep = strcpya_k(cswritep, "CHROM\t");
        char* write_iter = chrtoa(cip, x_code, x_name_buf);
        *write_iter++ = '\t';
        x_name_blen = write_iter - x_name_buf;
      }
      if (hardy_flags & kfHardyColPos) {
        cswritep = strcpya_k(cswritep, "POS\t");
      } else {
        variant_bps = nullptr;
      }
      cswritep = strcpya_k(cswritep, "ID");
      if (ref_col) {
        cswritep = strcpya_k(cswritep, "\tREF");
      }
      if (alt1_col) {
        cswritep = strcpya_k(cswritep, "\tALT1");
      }
      if (alt_col) {
        cswritep = strcpya_k(cswritep, "\tALT");
      }
      cswritep = strcpya_k(cswritep, "\tA1");
      if (ax_col) {
        cswritep = strcpya_k(cswritep, "\tAX");
      }
      if (gcounts) {
        if (gcount_1col) {
          cswritep = strcpya_k(cswritep, "\tGCOUNTS");
        } else {
          cswritep = strcpya_k(cswritep, "\tFEMALE_HOM_A1_CT\tFEMALE_HET_A1_CT\tFEMALE_TWO_AX_CT\tMALE_A1_CT\tMALE_AX_CT");
        }
      }
      if (hetfreq_cols) {
        cswritep = strcpya_k(cswritep, "\tO(FEMALE_HET_A1)\tE(FEMALE_HET_A1)");
      }
      const uint32_t sexaf_cols = hardy_flags & kfHardyColSexaf;
      if (sexaf_cols) {
        cswritep = strcpya_k(cswritep, "\tFEMALE_A1_FREQ\tMALE_A1_FREQ");
      }
      const uint32_t femalep_col = hardy_flags & kfHardyColFemalep;
      if (femalep_col) {
        cswritep = strcpya_k(cswritep, "\tFEMALE_ONLY_");
        if (midp) {
          cswritep = strcpya_k(cswritep, "MIDP");
        } else {
          *cswritep++ = 'P';
        }
      }
      if (p_col) {
        *cswritep++ = '\t';
        if (midp) {
          cswritep = strcpya_k(cswritep, "MIDP");
        } else {
          *cswritep++ = 'P';
        }
      }
      AppendBinaryEoln(&cswritep);
      fputs("--hardy: Writing chrX results...", stdout);
      fflush(stdout);
      const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
      const uint32_t x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
      uintptr_t variant_uidx_base;
      uintptr_t variant_include_bits;
      BitIter1Start(variant_include, x_start, &variant_uidx_base, &variant_include_bits);
      uintptr_t xgeno_idx = 0;
      uint32_t allele_ct = 2;
      uint32_t male_a1_ct = 0;
      uint32_t male_ax_ct = 0;
      for (uint32_t variant_idx = 0; variant_idx != hwe_x_ct; ++variant_idx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        const uint32_t print_allele_ct = ((allele_ct == 2) && (!redundant))? 1 : allele_ct;
        STD_ARRAY_KREF(uint32_t, 3) cur_geno_cts = hwe_geno_cts[variant_uidx];
        uint32_t female_obs_ct = cur_geno_cts[0] + cur_geno_cts[1] + cur_geno_cts[2];
        if (hwe_x_nosex_geno_cts) {
          STD_ARRAY_KREF(uint32_t, 3) cur_nosex_geno_cts = hwe_x_nosex_geno_cts[variant_uidx - x_start];
          female_obs_ct -= cur_nosex_geno_cts[0] + cur_nosex_geno_cts[1] + cur_nosex_geno_cts[2];
        }
        uint32_t male_obs_ct = 0;
        if (hwe_x_male_geno_cts) {
          STD_ARRAY_KREF(uint32_t, 3) cur_male_geno_cts = hwe_x_male_geno_cts[variant_uidx - x_start];
          male_obs_ct = cur_male_geno_cts[0] + cur_male_geno_cts[1] + cur_male_geno_cts[2];
          female_obs_ct -= male_obs_ct;
        }
        for (uint32_t a1_idx = 0; a1_idx != print_allele_ct; ++a1_idx) {
          cswritep = memcpya(cswritep, x_name_buf, x_name_blen);
          if (variant_bps) {
            cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
          }
          cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
          if (ref_col) {
            *cswritep++ = '\t';
            cswritep = strcpya(cswritep, cur_alleles[0]);
          }
          if (alt1_col) {
            *cswritep++ = '\t';
            cswritep = strcpya(cswritep, cur_alleles[1]);
          }
          if (alt_col) {
            *cswritep++ = '\t';
            for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto HardyReport_ret_WRITE_FAIL;
              }
              cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
            }
            --cswritep;
          }
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, cur_alleles[a1_idx]);
          if (ax_col) {
            *cswritep++ = '\t';
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if (allele_idx == a1_idx) {
                continue;
              }
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto HardyReport_ret_WRITE_FAIL;
              }
              cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
            }
            --cswritep;
          }
          uint32_t female_hom_a1_ct;
          uint32_t female_het_a1_ct;
          uint32_t female_two_ax_ct;
          if ((!a1_idx) || (allele_ct == 2)) {
            female_hom_a1_ct = cur_geno_cts[2 * a1_idx];
            female_het_a1_ct = cur_geno_cts[1];
            female_two_ax_ct = cur_geno_cts[2 - 2 * a1_idx];
            if (hwe_x_male_geno_cts) {
              STD_ARRAY_KREF(uint32_t, 3) cur_male_geno_cts = hwe_x_male_geno_cts[variant_uidx - x_start];
              male_a1_ct = cur_male_geno_cts[2 * a1_idx];
              female_hom_a1_ct -= male_a1_ct;
              female_het_a1_ct -= cur_male_geno_cts[1];
              male_ax_ct = cur_male_geno_cts[2 - 2 * a1_idx];
              female_two_ax_ct -= male_ax_ct;
            }
            if (hwe_x_nosex_geno_cts) {
              STD_ARRAY_KREF(uint32_t, 3) cur_nosex_geno_cts = hwe_x_nosex_geno_cts[variant_uidx - x_start];
              female_hom_a1_ct -= cur_nosex_geno_cts[2 * a1_idx];
              female_het_a1_ct -= cur_nosex_geno_cts[1];
              female_two_ax_ct -= cur_nosex_geno_cts[2 - 2 * a1_idx];
            }
          } else {
            STD_ARRAY_KREF(uint32_t, 2) cur_knownsex_xgeno_cts = x_knownsex_xgeno_cts[xgeno_idx];
            female_hom_a1_ct = cur_knownsex_xgeno_cts[0];
            female_het_a1_ct = cur_knownsex_xgeno_cts[1];
            if (x_male_xgeno_cts) {
              STD_ARRAY_KREF(uint32_t, 2) cur_male_xgeno_cts = x_male_xgeno_cts[xgeno_idx];
              male_a1_ct = cur_male_xgeno_cts[0];
              female_hom_a1_ct -= male_a1_ct;
              const uint32_t male_hethap_ct = cur_male_xgeno_cts[1];
              female_het_a1_ct -= male_hethap_ct;
              male_ax_ct = male_obs_ct - male_a1_ct - male_hethap_ct;
            }
            female_two_ax_ct = female_obs_ct - female_hom_a1_ct - female_het_a1_ct;
            // Correct to increment this before looking up hwe_x_pvals[] (and
            // to not increment on a1_idx == 0).
            ++xgeno_idx;
          }
          if (gcounts) {
            *cswritep++ = '\t';
            cswritep = u32toa_x(female_hom_a1_ct, gcount_delim, cswritep);
            cswritep = u32toa_x(female_het_a1_ct, gcount_delim, cswritep);
            cswritep = u32toa_x(female_two_ax_ct, gcount_delim, cswritep);
            cswritep = u32toa_x(male_a1_ct, gcount_delim, cswritep);
            cswritep = u32toa(male_ax_ct, cswritep);
          }
          if (hetfreq_cols || sexaf_cols) {
            const double female_obs_ct_recip = 1.0 / u31tod(female_obs_ct);
            const double dbl_a1_freq = (female_hom_a1_ct * 2 + female_het_a1_ct) * female_obs_ct_recip;
            if (hetfreq_cols) {
              *cswritep++ = '\t';
              cswritep = dtoa_g(u31tod(female_het_a1_ct) * female_obs_ct_recip, cswritep);
              *cswritep++ = '\t';
              // (1.0 - a1_freq) is vulnerable to catastrophic cancellation
              // when actual a1_freq is 1.
              if (female_hom_a1_ct == female_obs_ct) {
                *cswritep++ = '0';
              } else {
                const double expected_het_freq = dbl_a1_freq * (1.0 - 0.5 * dbl_a1_freq);
                cswritep = dtoa_g(expected_het_freq, cswritep);
              }
            }
            if (sexaf_cols) {
              *cswritep++ = '\t';
              cswritep = dtoa_g(dbl_a1_freq * 0.5, cswritep);
              *cswritep++ = '\t';
              const double male_a1_freq = u31tod(male_a1_ct) / u31tod(male_a1_ct + male_ax_ct);
              cswritep = dtoa_g(male_a1_freq, cswritep);
            }
          }
          if (femalep_col) {
            *cswritep++ = '\t';
            const double female_hwe_p = HweP(female_het_a1_ct, female_hom_a1_ct, female_two_ax_ct, midp);
            cswritep = dtoa_g(MAXV(female_hwe_p, output_min_p), cswritep);
          }
          if (p_col) {
            *cswritep++ = '\t';
            // bugfix (27 Jun 2020): forgot to correct this for multiallelic
            // variants
            cswritep = dtoa_g(MAXV(hwe_x_pvals[variant_idx + xgeno_idx], output_min_p), cswritep);
          }
          AppendBinaryEoln(&cswritep);
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto HardyReport_ret_WRITE_FAIL;
          }
        }
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto HardyReport_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      logprintfww("--hardy%s%s: chrX Hardy-Weinberg report (%s) written to %s .\n", output_zst? " zs" : "", midp? " midp" : "", nonfounders? "all samples" : "founders only", outname);
    }
  }
  while (0) {
  HardyReport_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  HardyReport_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 HardyReport_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

ENUM_U31_DEF_START()
  kSubstCodeNonsnpsymb,
  kSubstCodeSymbolic,
  kSubstCodeTs,
  kSubstCodeTv,
  kSubstCodeWeirdSnp,

  kSubstCodeCt
ENUM_U31_DEF_END(SubstCode);

SubstCode GetSubstCode(const char* ref, const char* alt) {
  if (ref[1] || alt[1]) {
    return (alt[0] == '<')? kSubstCodeSymbolic : kSubstCodeNonsnpsymb;
  }
  const unsigned char ref_upcase = ref[0] & 0xdf;
  const unsigned char alt_upcase = alt[0] & 0xdf;
  const unsigned char uchar_xor = ref_upcase ^ alt_upcase;
  if ((ref_upcase == 'A') || (ref_upcase == 'G')) {
    // 1 ^ 7 = 6
    if (uchar_xor == 6) {
      return kSubstCodeTs;
    }
    // bcftools sample-stats treats all SNP non-transitions are
    // transversions.  Don't replicate that.
    if ((alt_upcase == 'C') || (alt_upcase == 'T')) {
      return kSubstCodeTv;
    }
  } else if ((ref_upcase == 'C') || (ref_upcase == 'T')) {
    // 3 ^ 20 = 23
    if (uchar_xor == 23) {
      return kSubstCodeTs;
    }
    if ((alt_upcase == 'A') || (alt_upcase == 'G')) {
      return kSubstCodeTv;
    }
  }
  return kSubstCodeWeirdSnp;
}

// Assumes trailing bits have been zeroed or filled.
// Returns UINT32_MAX if no singleton.
uint32_t GetSingletonIdx(const uintptr_t* genovec, uint32_t sample_ct) {
  const uint32_t word_ct = NypCtToWordCt(sample_ct);
  uint32_t singleton_idx = UINT32_MAX;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    uintptr_t geno_word = genovec[widx];
    if (!geno_word) {
      continue;
    }
    geno_word = (geno_word ^ (geno_word >> 1)) & kMask5555;
    if (!geno_word) {
      continue;
    }
    if ((singleton_idx != UINT32_MAX) || (geno_word & (geno_word - 1))) {
      return UINT32_MAX;
    }
    singleton_idx = widx * kBitsPerWordD2 + (ctzw(geno_word) / 2);
  }
  return singleton_idx;
}

// If the ALT allele in a chrX biallelic variant is observed in one female and
// one male, bcftools counts it as a singleton for the female.  Default to
// replicating this weird behavior.
void UpdateSampleDiploidSingletonCountX(const uintptr_t* sex_male, const uintptr_t* genovec, uint32_t sample_ct, uint32_t* diploid_singleton_cts) {
  const uint32_t word_ct = NypCtToWordCt(sample_ct);
  const Halfword* sex_male_alias = R_CAST(const Halfword*, sex_male);
  uint32_t singleton_idx = UINT32_MAX;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    uintptr_t geno_word = genovec[widx];
    if (!geno_word) {
      continue;
    }
    geno_word = (geno_word ^ (geno_word >> 1)) & kMask5555;
    if (!geno_word) {
      continue;
    }
    geno_word = geno_word & (~UnpackHalfwordToWord(sex_male_alias[widx]));
    if (!geno_word) {
      continue;
    }
    if ((singleton_idx != UINT32_MAX) || (geno_word & (geno_word - 1))) {
      return;
    }
    singleton_idx = widx * kBitsPerWordD2 + ctzw(geno_word) / 2;
  }
  if (singleton_idx != UINT32_MAX) {
    diploid_singleton_cts[singleton_idx] += 1;
  }
}

void UpdateSampleSingletonCountX(const uintptr_t* sex_male, const uintptr_t* genovec, uint32_t sample_ct, uint32_t* singleton_cts) {
  const uint32_t word_ct = NypCtToWordCt(sample_ct);
  const Halfword* sex_male_alias = R_CAST(const Halfword*, sex_male);
  uint32_t singleton_idx = UINT32_MAX;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    if (!geno_word) {
      continue;
    }
    const uintptr_t geno2_word = (~geno_word) & (geno_word >> 1) & kMask5555;
    if (geno2_word) {
      if ((singleton_idx != UINT32_MAX) || (geno2_word & (geno2_word - 1))) {
        return;
      }
      singleton_idx = widx * kBitsPerWordD2 + ctzw(geno2_word) / 2;
    }
    uintptr_t geno1_word = geno_word & (~(geno_word >> 1)) & kMask5555;
    if (geno1_word) {
      geno1_word = geno1_word & (~UnpackHalfwordToWord(sex_male_alias[widx]));
      if (geno1_word) {
        if ((singleton_idx != UINT32_MAX) || (geno1_word & (geno1_word - 1))) {
          return;
        }
        singleton_idx = widx * kBitsPerWordD2 + ctzw(geno1_word) / 2;
      }
    }
  }
  if (singleton_idx != UINT32_MAX) {
    singleton_cts[singleton_idx] += 1;
  }
}

void UpdateSampleSingletonCountY(const uintptr_t* sex_male, const uintptr_t* genovec, uint32_t sample_ct, uint32_t* singleton_cts) {
  const uint32_t word_ct = NypCtToWordCt(sample_ct);
  const Halfword* sex_male_alias = R_CAST(const Halfword*, sex_male);
  uint32_t singleton_idx = UINT32_MAX;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    if (!geno_word) {
      continue;
    }
    uintptr_t geno2_word = (~geno_word) & (geno_word >> 1) & kMask5555;
    if (geno2_word) {
      geno2_word = geno2_word & UnpackHalfwordToWord(sex_male_alias[widx]);
      if (geno2_word) {
        if ((singleton_idx != UINT32_MAX) || (geno2_word & (geno2_word - 1))) {
          return;
        }
        singleton_idx = widx * kBitsPerWordD2 + ctzw(geno2_word) / 2;
      }
    }
  }
  if (singleton_idx != UINT32_MAX) {
    singleton_cts[singleton_idx] += 1;
  }
}

uint32_t GetSingletonIdxSparse(const uintptr_t* raregeno, const uint32_t* difflist_sample_ids, uint32_t difflist_len) {
  const uint32_t word_ct = NypCtToWordCt(difflist_len);
  uint32_t singleton_idx = UINT32_MAX;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    uintptr_t raregeno_word = raregeno[widx];
    if (!raregeno_word) {
      continue;
    }
    raregeno_word = (raregeno_word ^ (raregeno_word >> 1)) & kMask5555;
    if (!raregeno_word) {
      continue;
    }
    if ((singleton_idx != UINT32_MAX) || (raregeno_word & (raregeno_word - 1))) {
      return UINT32_MAX;
    }
    singleton_idx = difflist_sample_ids[widx * kBitsPerWordD2 + (ctzw(raregeno_word) / 2)];
  }
  return singleton_idx;
}

void UpdateSampleDiploidSingletonCountSparseX(const uintptr_t* sex_male, const uintptr_t* raregeno, const uint32_t* difflist_sample_ids, uint32_t difflist_len, uint32_t* diploid_singleton_cts) {
  const uint32_t word_ct = NypCtToWordCt(difflist_len);
  uint32_t singleton_idx = UINT32_MAX;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    uintptr_t raregeno_word = raregeno[widx];
    if (!raregeno_word) {
      continue;
    }
    raregeno_word = (raregeno_word ^ (raregeno_word >> 1)) & kMask5555;
    if (!raregeno_word) {
      continue;
    }
    const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
    do {
      const uint32_t difflist_idx_lowbits = ctzw(raregeno_word) / 2;
      const uint32_t sample_idx = cur_difflist_sample_ids[difflist_idx_lowbits];
      if (!IsSet(sex_male, sample_idx)) {
        if (singleton_idx != UINT32_MAX) {
          return;
        }
        singleton_idx = sample_idx;
      }
      raregeno_word = raregeno_word & (raregeno_word - 1);
    } while (raregeno_word);
  }
  if (singleton_idx != UINT32_MAX) {
    diploid_singleton_cts[singleton_idx] += 1;
  }
}

void UpdateSampleSingletonCountSparseX(const uintptr_t* sex_male, const uintptr_t* raregeno, const uint32_t* difflist_sample_ids, uint32_t difflist_len, uint32_t* singleton_cts) {
  const uint32_t word_ct = NypCtToWordCt(difflist_len);
  uint32_t singleton_idx = UINT32_MAX;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t raregeno_word = raregeno[widx];
    if (!raregeno_word) {
      continue;
    }
    const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
    const uintptr_t geno2_word = (~raregeno_word) & (raregeno_word >> 1) & kMask5555;
    if (geno2_word) {
      if ((singleton_idx != UINT32_MAX) || (geno2_word & (geno2_word - 1))) {
        return;
      }
      singleton_idx = cur_difflist_sample_ids[ctzw(geno2_word) / 2];
    }
    uintptr_t geno1_word = raregeno_word & (~(raregeno_word >> 1)) & kMask5555;
    while (geno1_word) {
      const uint32_t difflist_idx_lowbits = ctzw(geno1_word) / 2;
      const uint32_t sample_idx = cur_difflist_sample_ids[difflist_idx_lowbits];
      if (!IsSet(sex_male, sample_idx)) {
        if (singleton_idx != UINT32_MAX) {
          return;
        }
        singleton_idx = sample_idx;
      }
      geno1_word = geno1_word & (geno1_word - 1);
    }
  }
  if (singleton_idx != UINT32_MAX) {
    singleton_cts[singleton_idx] += 1;
  }
}

void UpdateSampleSingletonCountSparseY(const uintptr_t* sex_male, const uintptr_t* raregeno, const uint32_t* difflist_sample_ids, uint32_t difflist_len, uint32_t* singleton_cts) {
  const uint32_t word_ct = NypCtToWordCt(difflist_len);
  uint32_t singleton_idx = UINT32_MAX;
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    uintptr_t raregeno_word = raregeno[widx];
    if (!raregeno_word) {
      continue;
    }
    raregeno_word = (~raregeno_word) & (raregeno_word >> 1) & kMask5555;
    if (!raregeno_word) {
      continue;
    }
    const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
    do {
      const uint32_t difflist_idx_lowbits = ctzw(raregeno_word) / 2;
      const uint32_t sample_idx = cur_difflist_sample_ids[difflist_idx_lowbits];
      if (IsSet(sex_male, sample_idx)) {
        if (singleton_idx != UINT32_MAX) {
          return;
        }
        singleton_idx = sample_idx;
      }
      raregeno_word = raregeno_word & (raregeno_word - 1);
    } while (raregeno_word);
  }
  if (singleton_idx != UINT32_MAX) {
    singleton_cts[singleton_idx] += 1;
  }
}

void UpdateDenseSampleCounts3(const uintptr_t* genovec, uint32_t acc2_vec_ct, VecW* acc2_0, uint16_t* remainders) {
  const VecW m1 = VCONST_W(kMask5555);
  const VecW* genovvec = R_CAST(const VecW*, genovec);
  VecW* acc2_2 = &(acc2_0[acc2_vec_ct * 23]);
  VecW* acc2_1 = &(acc2_2[acc2_vec_ct * 23]);
  for (uint32_t vidx = 0; vidx != acc2_vec_ct; ++vidx) {
    // Tried iterating over word-indexes instead (since LEA instruction can
    // only multiply by 1/2/4/8, not 16 or 32), but that doesn't help, it just
    // makes the code messier.
    const VecW vv = genovvec[vidx];
    const VecW vv_hi = vecw_srli(vv, 1);
    const VecW vv_lo_0 = vecw_and_notfirst(vv, m1);
    acc2_0[vidx] += vecw_and_notfirst(vv_hi, vv_lo_0);
    acc2_2[vidx] += vv_lo_0 & vv_hi;
    acc2_1[vidx] += vecw_and_notfirst(vv_hi, vv) & m1;
  }
  remainders[0] -= 1;
  if (!remainders[0]) {
    VecW* acc4_0 = &(acc2_0[acc2_vec_ct]);
    Vcount0Incr2To4(acc2_vec_ct, acc2_0, acc4_0);
    VecW* acc4_2 = &(acc2_2[acc2_vec_ct]);
    Vcount0Incr2To4(acc2_vec_ct, acc2_2, acc4_2);
    VecW* acc4_1 = &(acc2_1[acc2_vec_ct]);
    Vcount0Incr2To4(acc2_vec_ct, acc2_1, acc4_1);
    remainders[1] -= 1;
    if (!remainders[1]) {
      const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
      VecW* acc8_0 = &(acc4_0[acc4_vec_ct]);
      Vcount0Incr4To8(acc4_vec_ct, acc4_0, acc8_0);
      VecW* acc8_2 = &(acc4_2[acc4_vec_ct]);
      Vcount0Incr4To8(acc4_vec_ct, acc4_2, acc8_2);
      VecW* acc8_1 = &(acc4_1[acc4_vec_ct]);
      Vcount0Incr4To8(acc4_vec_ct, acc4_1, acc8_1);
      remainders[2] -= 1;
      if (!remainders[2]) {
        const uint32_t acc8_vec_ct = acc4_vec_ct * 2;
        VecW* acc32_0 = &(acc8_0[acc8_vec_ct]);
        Vcount0Incr8To32(acc8_vec_ct, acc8_0, acc32_0);
        VecW* acc32_2 = &(acc8_2[acc8_vec_ct]);
        Vcount0Incr8To32(acc8_vec_ct, acc8_2, acc32_2);
        VecW* acc32_1 = &(acc8_1[acc8_vec_ct]);
        Vcount0Incr8To32(acc8_vec_ct, acc8_1, acc32_1);
        remainders[2] = 17;
      }
      remainders[1] = 5;
    }
    remainders[0] = 3;
  }
}

void UpdateDenseSampleCounts2(const uintptr_t* genovec, uint32_t acc2_vec_ct, VecW* acc2_0, uint16_t* remainders) {
  const VecW m1 = VCONST_W(kMask5555);
  const VecW* genovvec = R_CAST(const VecW*, genovec);
  VecW* acc2_2 = &(acc2_0[acc2_vec_ct * 23]);
  for (uint32_t vidx = 0; vidx != acc2_vec_ct; ++vidx) {
    const VecW vv = genovvec[vidx];
    const VecW vv_hi = vecw_srli(vv, 1);
    const VecW vv_lo_0 = vecw_and_notfirst(vv, m1);
    acc2_0[vidx] += vecw_and_notfirst(vv_hi, vv_lo_0);
    acc2_2[vidx] += vv_lo_0 & vv_hi;
  }
  remainders[0] -= 1;
  if (!remainders[0]) {
    VecW* acc4_0 = &(acc2_0[acc2_vec_ct]);
    Vcount0Incr2To4(acc2_vec_ct, acc2_0, acc4_0);
    VecW* acc4_2 = &(acc2_2[acc2_vec_ct]);
    Vcount0Incr2To4(acc2_vec_ct, acc2_2, acc4_2);
    remainders[1] -= 1;
    if (!remainders[1]) {
      const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
      VecW* acc8_0 = &(acc4_0[acc4_vec_ct]);
      Vcount0Incr4To8(acc4_vec_ct, acc4_0, acc8_0);
      VecW* acc8_2 = &(acc4_2[acc4_vec_ct]);
      Vcount0Incr4To8(acc4_vec_ct, acc4_2, acc8_2);
      remainders[2] -= 1;
      if (!remainders[2]) {
        const uint32_t acc8_vec_ct = acc4_vec_ct * 2;
        VecW* acc32_0 = &(acc8_0[acc8_vec_ct]);
        Vcount0Incr8To32(acc8_vec_ct, acc8_0, acc32_0);
        VecW* acc32_2 = &(acc8_2[acc8_vec_ct]);
        Vcount0Incr8To32(acc8_vec_ct, acc8_2, acc32_2);
        remainders[2] = 17;
      }
      remainders[1] = 5;
    }
    remainders[0] = 3;
  }
}

// This shares a lot with LoadSampleMissingCts() in plink2_filter.  Might want
// to merge this into a single driver function at some point, similar to how
// LoadAlleleAndGenoCounts() reduces duplication of variant-stat calculation.
typedef struct SampleCountsCtxStruct {
  const uintptr_t* variant_include;
  const ChrInfo* cip;
  const uintptr_t* allele_idx_offsets;
  const char* const* allele_storage;
  const uintptr_t* sample_include;
  uint32_t* sample_include_cumulative_popcounts;
  uintptr_t* sex_male_collapsed;
  uint32_t sample_ct;
  uint32_t male_ct;
  uint32_t y_nonmale_needed;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** raregenos;
  uint32_t** difflist_sample_id_bufs;
  uint32_t* read_variant_uidx_starts;
  uint32_t cur_block_size;

  // only kPglRetMalformedInput possible, no atomic ops needed
  PglErr reterr;

  // top-level: length calc_thread_ct array
  // second level: length-16 arrays, corresponding to the 16 chr_type x
  //               subst_code possibilities
  // bottom level: acc2 partial counts, then acc4, then acc8, then acc32 for
  //               genotype=0.  Then the same for genotype=2, and finally (in
  //               autosomal-diploid/chrX cases) genotype=1.
  VecW*** thread_dense_counts;

  // top-level: length calc_thread_ct array
  // second level: length-32 arrays, corresponding to chr_type x common_geno x
  //               subst_code possibilities
  // bottom level: three vector-aligned length-sample_ct arrays, corresponding
  //               to the three ((common_geno ^ rare_geno) - 1) possibilities.
  uint32_t*** thread_sparse_counts;
  uint32_t** thread_sparse_common0_cts;

  uint32_t** thread_diploid_singleton_cts;
  uint32_t** thread_singleton_cts;

  uint16_t** thread_alt_subst_codes;
  // Multiallelic complications: geno=2 entries can correspond to het
  // altx/alty, and subst_code can vary between alt alleles.  We track the
  // necessary adjustments with the het_rarealt_cts, het2alt_cts,
  // hom_rarealt_cts, and hap_rarealt_cts arrays below.
  // - het2alt_cts has 10 vector-aligned length-sample_ct arrays.  The first
  //   four correspond to ALT1 subst_codes, and the last four correspond to the
  //   actually-observed-ALT subst_codes.  The sum of the last four arrays is
  //   always twice that of the first four arrays, since each allele is counted
  //   separately.
  // - het_rarealt_cts, hom_rarealt_cts and hap_rarealt_cts have 5
  //   vector-aligned lenght-sample_ct arrays, corresponding to subst_code
  //   frequency deltas.
  uint32_t** thread_het2alt_cts;
  int32_t** thread_het_rarealt_cts;
  int32_t** thread_hom_rarealt_cts;
  int32_t** thread_hap_rarealt_cts;
} SampleCountsCtx;

THREAD_FUNC_DECL SampleCountsThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  SampleCountsCtx* ctx = S_CAST(SampleCountsCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const ChrInfo* cip = ctx->cip;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const char* const* allele_storage = ctx->allele_storage;
  const uintptr_t* sample_include = ctx->sample_include;
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  const uintptr_t* sex_male_collapsed = ctx->sex_male_collapsed;
  const uintptr_t sample_ct = ctx->sample_ct;
  const uint32_t male_ct = ctx->male_ct;
  const uint32_t skip_y = (!male_ct) && (!ctx->y_nonmale_needed);
  uintptr_t* genovec = ctx->genovecs[tidx];
  PgenVariant pgv;
  pgv.genovec = genovec;
  SetPgvThreadMhcNull(sample_ct, tidx, ctx->thread_read_mhc, &pgv);
  uintptr_t* raregeno = ctx->raregenos[tidx];
  uint32_t* difflist_sample_ids = ctx->difflist_sample_id_bufs[tidx];

  // todo: tune this threshold
  const uint32_t max_simple_difflist_len = sample_ct / 32;

  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t acc2_vec_ct = NypCtToVecCt(sample_ct);
  const uintptr_t dense_counts_vstride = acc2_vec_ct * 23;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  VecW** dense_counts = ctx->thread_dense_counts[tidx];
  // 4 chr_types, 3 remainders per variant_type
  // uint16_t instead of uint8_t to avoid aliasing paranoia
  uint16_t dense_remainders[12 * kSubstCodeCt];
  for (uint32_t dense_vtype = 0; dense_vtype != 4 * kSubstCodeCt; ++dense_vtype) {
    if (dense_counts[dense_vtype]) {
      ZeroVecArr(dense_counts_vstride * (2 + (dense_vtype < 2 * kSubstCodeCt)), dense_counts[dense_vtype]);
      dense_remainders[3 * dense_vtype] = 3;
      dense_remainders[3 * dense_vtype + 1] = 5;
      dense_remainders[3 * dense_vtype + 2] = 17;
    }
  }

  const uintptr_t sample_ct_i32av = RoundUpPow2(sample_ct, kInt32PerVec);
  uint32_t** sparse_counts = ctx->thread_sparse_counts[tidx];
  // 4 chr_types * 2 common_geno
  for (uint32_t uii = 0; uii != 8 * kSubstCodeCt; ++uii) {
    if (sparse_counts[uii]) {
      ZeroU32Arr(3 * sample_ct_i32av, sparse_counts[uii]);
    }
  }
  uint32_t* sparse_common0_cts = ctx->thread_sparse_common0_cts[tidx];
  ZeroU32Arr(4 * kSubstCodeCt, sparse_common0_cts);

  uint32_t* diploid_singleton_cts = nullptr;
  if (ctx->thread_diploid_singleton_cts) {
    diploid_singleton_cts = ctx->thread_diploid_singleton_cts[tidx];
    ZeroU32Arr(sample_ct_i32av, diploid_singleton_cts);
  }
  uint32_t* singleton_cts = nullptr;
  if (ctx->thread_singleton_cts) {
    singleton_cts = ctx->thread_singleton_cts[tidx];
    ZeroU32Arr(sample_ct_i32av, singleton_cts);
  }

  uint16_t* alt_subst_codes = nullptr;
  uint32_t* het2alt_minus_cts = nullptr;
  uint32_t* het2alt_plus_cts = nullptr;
  int32_t* het_rarealt_cts = nullptr;
  int32_t* hom_rarealt_cts = nullptr;
  int32_t* hap_rarealt_cts = nullptr;
  if (ctx->thread_alt_subst_codes) {
    alt_subst_codes = ctx->thread_alt_subst_codes[tidx];

    het2alt_minus_cts = ctx->thread_het2alt_cts[tidx];
    ZeroU32Arr(2 * kSubstCodeCt * sample_ct_i32av, het2alt_minus_cts);
    het2alt_plus_cts = &(het2alt_minus_cts[kSubstCodeCt * sample_ct_i32av]);
    het_rarealt_cts = ctx->thread_het_rarealt_cts[tidx];
    ZeroI32Arr(kSubstCodeCt * sample_ct_i32av, het_rarealt_cts);
    hom_rarealt_cts = ctx->thread_hom_rarealt_cts[tidx];
    ZeroI32Arr(kSubstCodeCt * sample_ct_i32av, hom_rarealt_cts);
    if (ctx->thread_hap_rarealt_cts) {
      hap_rarealt_cts = ctx->thread_hap_rarealt_cts[tidx];
      ZeroI32Arr(kSubstCodeCt * sample_ct_i32av, hap_rarealt_cts);
    }
  }

  uint32_t cur_allele_ct = 2;
  do {
    const uint32_t cur_block_size = ctx->cur_block_size;
    const uint32_t cur_idx_ct = (((tidx + 1) * cur_block_size) / calc_thread_ct) - ((tidx * cur_block_size) / calc_thread_ct);
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
    uint32_t chr_end = 0;
    uint32_t chr_type = 0;
    uint32_t is_haploid = 0;
    uint32_t is_diploid_x = 0;
    uint32_t is_y = 0;
    for (uint32_t cur_idx = 0; cur_idx != cur_idx_ct; ++cur_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        is_haploid = IsSet(cip->haploid_mask, chr_idx);
        is_diploid_x = 0;
        is_y = (chr_idx == y_code);
        if (chr_idx == x_code) {
          is_diploid_x = !IsSet(cip->haploid_mask, 0);
        }
        chr_type = 0;
        if (is_haploid) {
          if (is_diploid_x) {
            if (!male_ct) {
              // all-unknown-sex override
              is_haploid = 0;
              is_diploid_x = 0;
            } else {
              chr_type = 1;
            }
          } else if (is_y) {
            chr_type = 2;
          } else {
            chr_type = 3;
          }
        }
      }
      if (is_y && skip_y) {
        continue;
      }
      uintptr_t allele_idx_offset_base;
      if (!allele_idx_offsets) {
        allele_idx_offset_base = variant_uidx * 2;
      } else {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (cur_allele_ct == 2) {
        // Finally, a scenario where exposing this capability really pays off.
        uint32_t difflist_common_geno;
        uint32_t difflist_len;
        PglErr reterr = PgrGetDifflistOrGenovec(sample_include, pssi, sample_ct, max_simple_difflist_len, variant_uidx, pgrp, genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
        if (unlikely(reterr)) {
          ctx->reterr = reterr;
          break;
        }
        if (difflist_common_geno == 2) {
          // Don't bother with sparse optimization here, should be rare and it
          // complicates singleton-counting code
          PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, sample_ct, difflist_len, genovec);
          difflist_common_geno = UINT32_MAX;
        }
        uint32_t subst_code = GetSubstCode(cur_alleles[0], cur_alleles[1]);
        const uint32_t dense_vtype = chr_type * kSubstCodeCt + subst_code;
        if (difflist_common_geno == UINT32_MAX) {
          if (chr_type < 2) {
            UpdateDenseSampleCounts3(genovec, acc2_vec_ct, dense_counts[dense_vtype], &(dense_remainders[dense_vtype * 3]));
          } else {
            UpdateDenseSampleCounts2(genovec, acc2_vec_ct, dense_counts[dense_vtype], &(dense_remainders[dense_vtype * 3]));
          }
          if (diploid_singleton_cts || singleton_cts) {
            ZeroTrailingNyps(sample_ct, genovec);
            if (is_haploid) {
              if (is_diploid_x) {
                if (diploid_singleton_cts) {
                  UpdateSampleDiploidSingletonCountX(sex_male_collapsed, genovec, sample_ct, diploid_singleton_cts);
                }
                if (singleton_cts) {
                  UpdateSampleSingletonCountX(sex_male_collapsed, genovec, sample_ct, singleton_cts);
                }
                continue;
              }
              if (!singleton_cts) {
                continue;
              }
              if (is_y) {
                UpdateSampleSingletonCountY(sex_male_collapsed, genovec, sample_ct, singleton_cts);
                continue;
              }
              SetHetMissing(sample_ctl2, genovec);
              const uint32_t singleton_idx = GetSingletonIdx(genovec, sample_ct);
              if (singleton_idx != UINT32_MAX) {
                singleton_cts[singleton_idx] += 1;
              }
              continue;
            }
            const uint32_t singleton_idx = GetSingletonIdx(genovec, sample_ct);
            if (singleton_idx != UINT32_MAX) {
              if (diploid_singleton_cts) {
                diploid_singleton_cts[singleton_idx] += 1;
              }
              if (singleton_cts) {
                singleton_cts[singleton_idx] += 1;
              }
            }
          }
          continue;
        }
        uint32_t vtrans_type = dense_vtype + chr_type * kSubstCodeCt;
        if (!difflist_common_geno) {
          sparse_common0_cts[dense_vtype] += 1;
        } else {
          vtrans_type += kSubstCodeCt;
        }
        if (!difflist_len) {
          continue;
        }
        uint32_t* cur_sparse_cts = sparse_counts[vtrans_type];
        const uint32_t word_ct_m1 = (difflist_len - 1) / kBitsPerWordD2;
        // difflist_common_geno_word == 0 if difflist_common_geno == 0,
        // ~k0LU if difflist_common_geno == 3.
        const uintptr_t difflist_common_geno_word = -S_CAST(uintptr_t, difflist_common_geno & 1);
        uint32_t loop_len = kBitsPerWordD2;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= word_ct_m1) {
            if (widx > word_ct_m1) {
              break;
            }
            loop_len = ModNz(difflist_len, kBitsPerWordD2);
          }
          uintptr_t raregeno_xor_word = raregeno[widx] ^ difflist_common_geno_word;
          const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
          for (uint32_t uii = 0; uii != loop_len; ++uii) {
            const uintptr_t sample_idx = cur_difflist_sample_ids[uii];
            const uint32_t cur_geno_xor = raregeno_xor_word & 3;
            cur_sparse_cts[(cur_geno_xor - 1) * sample_ct_i32av + sample_idx] += 1;
            raregeno_xor_word = raregeno_xor_word >> 2;
          }
        }
        if (diploid_singleton_cts || singleton_cts) {
          ZeroTrailingNyps(difflist_len, raregeno);
          if (is_haploid) {
            if (is_diploid_x) {
              if (diploid_singleton_cts) {
                UpdateSampleDiploidSingletonCountSparseX(sex_male_collapsed, raregeno, difflist_sample_ids, difflist_len, diploid_singleton_cts);
              }
              if (singleton_cts) {
                UpdateSampleSingletonCountSparseX(sex_male_collapsed, raregeno, difflist_sample_ids, difflist_len, singleton_cts);
              }
              continue;
            }
            if (!singleton_cts) {
              continue;
            }
            if (is_y) {
              UpdateSampleSingletonCountSparseY(sex_male_collapsed, raregeno, difflist_sample_ids, difflist_len, singleton_cts);
              continue;
            }
            const uint32_t difflist_word_ct = NypCtToWordCt(difflist_len);
            SetHetMissing(difflist_word_ct, raregeno);
            const uint32_t singleton_idx = GetSingletonIdxSparse(raregeno, difflist_sample_ids, difflist_len);
            if (singleton_idx != UINT32_MAX) {
              singleton_cts[singleton_idx] += 1;
            }
            continue;
          }
          const uint32_t singleton_idx = GetSingletonIdxSparse(raregeno, difflist_sample_ids, difflist_len);
          if (singleton_idx != UINT32_MAX) {
            if (diploid_singleton_cts) {
              diploid_singleton_cts[singleton_idx] += 1;
            }
            if (singleton_cts) {
              singleton_cts[singleton_idx] += 1;
            }
          }
        }
        continue;
      }
      // Multiallelic case.
      PglErr reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, pgrp, &pgv);
      if (unlikely(reterr)) {
        ctx->reterr = reterr;
        break;
      }
      const uint32_t subst_code1 = GetSubstCode(cur_alleles[0], cur_alleles[1]);
      alt_subst_codes[1] = subst_code1;
      uint32_t subst_codes_vary = 0;
      for (uint32_t allele_idx = 2; allele_idx != cur_allele_ct; ++allele_idx) {
        const uint32_t cur_subst_code = GetSubstCode(cur_alleles[0], cur_alleles[allele_idx]);
        alt_subst_codes[allele_idx] = cur_subst_code;
        subst_codes_vary |= (cur_subst_code != subst_code1);
      }
      const uint32_t dense_vtype = chr_type * kSubstCodeCt + subst_code1;
      if (chr_type < 2) {
        UpdateDenseSampleCounts3(genovec, acc2_vec_ct, dense_counts[dense_vtype], &(dense_remainders[dense_vtype * 3]));
      } else {
        UpdateDenseSampleCounts2(genovec, acc2_vec_ct, dense_counts[dense_vtype], &(dense_remainders[dense_vtype * 3]));
      }
      const uint32_t patch_10_ct = pgv.patch_10_ct;
      const uint32_t subst_code1_offset = subst_code1 * sample_ct_i32av;
      if (!subst_codes_vary) {
        if (patch_10_ct) {
          const uintptr_t* patch_10_set = pgv.patch_10_set;
          const AlleleCode* patch_10_vals = pgv.patch_10_vals;
          uintptr_t sample_idx_base = 0;
          uintptr_t sample_idx_bits = patch_10_set[0];
          if (!chr_type) {
            for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
              const AlleleCode a0 = patch_10_vals[uii * 2];
              const AlleleCode a1 = patch_10_vals[uii * 2 + 1];
              if (a0 != a1) {
                het2alt_minus_cts[subst_code1_offset + sample_idx] += 1;
                het2alt_plus_cts[subst_code1_offset + sample_idx] += 2;
              }
            }
          } else if (chr_type == 1) {
            for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
              const uint32_t is_male = IsSet(sex_male_collapsed, sample_idx);
              const AlleleCode a0 = patch_10_vals[uii * 2];
              const AlleleCode a1 = patch_10_vals[uii * 2 + 1];
              if (a0 != a1) {
                if (!is_male) {
                  het2alt_minus_cts[subst_code1_offset + sample_idx] += 1;
                  het2alt_plus_cts[subst_code1_offset + sample_idx] += 2;
                } else {
                  // het-haploid, treat as missing
                  hap_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
                }
              }
            }
          } else {
            for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
              const AlleleCode a0 = patch_10_vals[uii * 2];
              const AlleleCode a1 = patch_10_vals[uii * 2 + 1];
              if (a0 != a1) {
                hap_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
              }
            }
          }
        }
      } else {
        if (chr_type < 2) {
          const uint32_t patch_01_ct = pgv.patch_01_ct;
          if (patch_01_ct) {
            const uintptr_t* patch_01_set = pgv.patch_01_set;
            const AlleleCode* patch_01_vals = pgv.patch_01_vals;
            uintptr_t sample_idx_base = 0;
            uintptr_t sample_idx_bits = patch_01_set[0];
            for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &sample_idx_bits);
              const uint32_t cur_subst_code = alt_subst_codes[patch_01_vals[uii]];
              if (subst_code1 != cur_subst_code) {
                het_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
                het_rarealt_cts[cur_subst_code * sample_ct_i32av + sample_idx] += 1;
              }
            }
          }
        }
        if (patch_10_ct) {
          const uintptr_t* patch_10_set = pgv.patch_10_set;
          const AlleleCode* patch_10_vals = pgv.patch_10_vals;
          uintptr_t sample_idx_base = 0;
          uintptr_t sample_idx_bits = patch_10_set[0];
          if (!chr_type) {
            for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
              const AlleleCode a0 = patch_10_vals[uii * 2];
              const AlleleCode a1 = patch_10_vals[uii * 2 + 1];
              const uint32_t new_code0 = alt_subst_codes[a0];
              if (a0 != a1) {
                const uint32_t new_code1 = alt_subst_codes[a1];
                het2alt_minus_cts[subst_code1_offset + sample_idx] += 1;
                het2alt_plus_cts[new_code0 * sample_ct_i32av + sample_idx] += 1;
                het2alt_plus_cts[new_code1 * sample_ct_i32av + sample_idx] += 1;
              } else if (new_code0 != subst_code1) {
                hom_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
                hom_rarealt_cts[new_code0 * sample_ct_i32av + sample_idx] += 1;
              }
            }
          } else if (chr_type == 1) {
            // ugh
            for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
              const uint32_t is_male = IsSet(sex_male_collapsed, sample_idx);
              const AlleleCode a0 = patch_10_vals[uii * 2];
              const AlleleCode a1 = patch_10_vals[uii * 2 + 1];
              const uint32_t new_code0 = alt_subst_codes[a0];
              if (a0 != a1) {
                if (!is_male) {
                  const uint32_t new_code1 = alt_subst_codes[a1];
                  het2alt_minus_cts[subst_code1_offset + sample_idx] += 1;
                  het2alt_plus_cts[new_code0 * sample_ct_i32av + sample_idx] += 1;
                  het2alt_plus_cts[new_code1 * sample_ct_i32av + sample_idx] += 1;
                } else {
                  hap_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
                }
              } else if (new_code0 != subst_code1) {
                if (!is_male) {
                  hom_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
                  hom_rarealt_cts[new_code0 * sample_ct_i32av + sample_idx] += 1;
                } else {
                  hap_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
                  hap_rarealt_cts[new_code0 * sample_ct_i32av + sample_idx] += 1;
                }
              }
            }
          } else {
            for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
              const AlleleCode a0 = patch_10_vals[uii * 2];
              const AlleleCode a1 = patch_10_vals[uii * 2 + 1];
              const uint32_t new_code0 = alt_subst_codes[a0];
              if (a0 != a1) {
                hap_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
              } else if (new_code0 != subst_code1) {
                hap_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
                hap_rarealt_cts[new_code0 * sample_ct_i32av + sample_idx] += 1;
              }
            }
          }
        }
      }
      if (diploid_singleton_cts || singleton_cts) {
        // this logic is completely unchanged, basic dense case should fall
        // through to this?
        ZeroTrailingNyps(sample_ct, genovec);
        if (is_haploid) {
          if (is_diploid_x) {
            if (diploid_singleton_cts) {
              UpdateSampleDiploidSingletonCountX(sex_male_collapsed, genovec, sample_ct, diploid_singleton_cts);
            }
            if (singleton_cts) {
              UpdateSampleSingletonCountX(sex_male_collapsed, genovec, sample_ct, singleton_cts);
            }
            continue;
          }
          if (!singleton_cts) {
            continue;
          }
          if (is_y) {
            UpdateSampleSingletonCountY(sex_male_collapsed, genovec, sample_ct, singleton_cts);
            continue;
          }
          SetHetMissing(sample_ctl2, genovec);
          const uint32_t singleton_idx = GetSingletonIdx(genovec, sample_ct);
          if (singleton_idx != UINT32_MAX) {
            singleton_cts[singleton_idx] += 1;
          }
          continue;
        }
        const uint32_t singleton_idx = GetSingletonIdx(genovec, sample_ct);
        if (singleton_idx != UINT32_MAX) {
          if (diploid_singleton_cts) {
            diploid_singleton_cts[singleton_idx] += 1;
          }
          if (singleton_cts) {
            singleton_cts[singleton_idx] += 1;
          }
        }
      }
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  const uintptr_t acc4_vec_ct = acc2_vec_ct * 2;
  const uintptr_t acc8_vec_ct = acc2_vec_ct * 4;
  for (uint32_t dense_vtype = 0; dense_vtype != 4 * kSubstCodeCt; ++dense_vtype) {
    VecW* acc2_0 = dense_counts[dense_vtype];
    if (acc2_0) {
      VecW* acc4_0 = &(acc2_0[acc2_vec_ct]);
      VcountIncr2To4(acc2_0, acc2_vec_ct, acc4_0);
      VecW* acc8_0 = &(acc4_0[acc4_vec_ct]);
      VcountIncr4To8(acc4_0, acc4_vec_ct, acc8_0);
      VcountIncr8To32(acc8_0, acc8_vec_ct, &(acc8_0[acc8_vec_ct]));

      VecW* acc2_2 = &(acc2_0[dense_counts_vstride]);
      VecW* acc4_2 = &(acc2_2[acc2_vec_ct]);
      VcountIncr2To4(acc2_2, acc2_vec_ct, acc4_2);
      VecW* acc8_2 = &(acc4_2[acc4_vec_ct]);
      VcountIncr4To8(acc4_2, acc4_vec_ct, acc8_2);
      VcountIncr8To32(acc8_2, acc8_vec_ct, &(acc8_2[acc8_vec_ct]));

      if (dense_vtype < 2 * kSubstCodeCt) {
        VecW* acc2_1 = &(acc2_2[dense_counts_vstride]);
        VecW* acc4_1 = &(acc2_1[acc2_vec_ct]);
        VcountIncr2To4(acc2_1, acc2_vec_ct, acc4_1);
        VecW* acc8_1 = &(acc4_1[acc4_vec_ct]);
        VcountIncr4To8(acc4_1, acc4_vec_ct, acc8_1);
        VcountIncr8To32(acc8_1, acc8_vec_ct, &(acc8_1[acc8_vec_ct]));
      }
    }
  }
  THREAD_RETURN;
}

void Unscramble2(uint32_t sample_ct, uint32_t* dst, uint32_t* unscramble_buf) {
  const uint32_t acc2_vec_ct = NypCtToWordCt(sample_ct);
  memcpy(unscramble_buf, dst, acc2_vec_ct * 16 * kBytesPerVec);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uint32_t scrambled_idx = VcountScramble2(sample_idx);
    dst[sample_idx] = unscramble_buf[scrambled_idx];
  }
}

// diploid-only if haploid_counts_subst == nullptr
void IncrSubstType(const SampleCountsCtx* ctx, const uint32_t* const* diploid_counts_subst, const uint32_t* const* haploid_counts_subst, uint32_t sample_ct_i32v, uint32_t subst_code, uint32_t* dst) {
  {
    const uint32_t* dip_src1 = diploid_counts_subst[1];
    const uint32_t* dip_src2 = diploid_counts_subst[2];
    U32CastVecAdd2(dip_src1, dip_src2, sample_ct_i32v, dst);
    if (haploid_counts_subst) {
      const uint32_t* hap_src = haploid_counts_subst[1];
      if (hap_src) {
        U32CastVecAdd(hap_src, sample_ct_i32v, dst);
      }
    }
  }
  if (ctx->thread_het2alt_cts) {
    const uintptr_t sample_ct_i32av = sample_ct_i32v * kInt32PerVec;
    const uint32_t* sub_src = &(ctx->thread_het2alt_cts[0][subst_code * sample_ct_i32av]);
    const uint32_t* add_src = &(sub_src[kSubstCodeCt * sample_ct_i32av]);
    U32CastVecAddSub(add_src, sub_src, sample_ct_i32v, dst);
  }
}

uint32_t* AllocAndCountSubstType(const SampleCountsCtx* ctx, const uint32_t* const* diploid_counts_subst, const uint32_t* const* haploid_counts_subst, uint32_t sample_ct_i32v, uint32_t subst_code) {
  uint32_t* dst = S_CAST(uint32_t*, bigstack_alloc(sample_ct_i32v * sizeof(VecU32)));
  if (unlikely(!dst)) {
    return nullptr;
  }
  ZeroU32Arr(sample_ct_i32v * kInt32PerVec, dst);
  IncrSubstType(ctx, diploid_counts_subst, haploid_counts_subst, sample_ct_i32v, subst_code, dst);
  return dst;
}

ENUM_U31_DEF_START()
  kSampleCountHom,
  kSampleCountHomref,
  kSampleCountHomalt,
  kSampleCountHomaltSnp,
  kSampleCountHet,
  kSampleCountRefalt,
  kSampleCountHet2alt,

  kSampleCountHetSnp,
  kSampleCountDiploidTs,
  kSampleCountTs,
  kSampleCountDiploidTv,
  kSampleCountTv,
  kSampleCountDiploidNonsnpsymb,
  kSampleCountNonsnpsymb,
  kSampleCountSymbolic,
  kSampleCountNonsnp,

  kSampleCountDiploidSingle,
  kSampleCountSingle,
  kSampleCountHaprefWithFemaleY,
  kSampleCountHapref,
  kSampleCountHapaltWithFemaleY,
  kSampleCountHapalt,
  kSampleCountMissingWithFemaleY,
  kSampleCountMissing,

  kSampleCountTypeCt
ENUM_U31_DEF_END(SampleCountType);

PglErr SampleCounts(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, SampleCountsFlags flags, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  const char sample_counts_headers[kSampleCountTypeCt][30] = {
    "HOM_CT",
    "HOM_REF_CT",
    "HOM_ALT_CT",
    "HOM_ALT_SNP_CT",
    "HET_CT",
    "HET_REF_ALT_CT",
    "HET_2ALT_CT",

    "HET_SNP_CT",
    "DIPLOID_TRANSITION_CT",
    "TRANSITION_CT",
    "DIPLOID_TRANSVERSION_CT",
    "TRANSVERSION_CT",
    "DIPLOID_NONSNP_NONSYMBOLIC_CT",
    "NONSNP_NONSYMBOLIC_CT",
    "SYMBOLIC_CT",
    "NONSNP_CT",

    "DIPLOID_SINGLETON_CT",
    "SINGLETON_CT",
    "HAP_REF_INCL_FEMALE_Y_CT",
    "HAP_REF_CT",
    "HAP_ALT_INCL_FEMALE_Y_CT",
    "HAP_ALT_CT",
    "MISSING_INCL_FEMALE_Y_CT",
    "MISSING_CT"
  };
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  ThreadGroup tg;
  PreinitThreads(&tg);
  SampleCountsCtx ctx;
  {
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    if (unlikely(
            bigstack_alloc_u32(raw_sample_ctl, &ctx.sample_include_cumulative_popcounts))) {
      goto SampleCounts_ret_NOMEM;
    }
    FillCumulativePopcounts(sample_include, raw_sample_ctl, ctx.sample_include_cumulative_popcounts);

    // Biallelic variants are handled as follows:
    // 1. Define four chromosome types (autosomal-diploid, chrX-diploid, chrY,
    //    and generic-haploid), and four substitution types (non-SNP,
    //    SNP-transition, SNP-transversion, non-ACGT 'SNP').  We compute
    //    per-sample 0/1/2 genotype counts for each (chromosome type,
    //    substitution type) pair.  This is accelerated by distinguishing the
    //    dense- and sparse-storage cases, and also parallelized across
    //    variants.
    // 2. Given these raw genotype counts, it is straightforward to build all
    //    of the supported columns, except for SINGLETON_CT.  (MISSING_CT can
    //    be obtained by subtraction from the total number of variants,
    //    excluding chrY for nonmales.)
    // 3. SINGLETON_CT gets its own logic in the main loop.
    uint32_t chr_types = 0;
    uint32_t diploid_chr_type_ct = 0;
    if (IsAutosomalDiploidChrPresent(cip)) {
      chr_types = 1;
      ++diploid_chr_type_ct;
    }
    uint32_t haploid_chr_type_ct = 0;
    uint32_t y_variant_ct = 0;
    uint32_t y_nonmale_needed = 0;
    {
      uint32_t y_code;
      if (XymtExists(cip, kChrOffsetY, &y_code)) {
        y_variant_ct = CountChrVariantsUnsafe(variant_include, cip, y_code);
      }
      uintptr_t haploid_mask_lowbit = cip->haploid_mask[0] & 1;
      if (male_ct) {
        uint32_t x_code;
        if ((!haploid_mask_lowbit) && XymtExists(cip, kChrOffsetX, &x_code)) {
          chr_types |= 2;
          ++diploid_chr_type_ct;
        }
      }
      if (y_variant_ct) {
        y_nonmale_needed = ((flags & (kfSampleCountsColHaprefWithFemaleY | kfSampleCountsColHapaltWithFemaleY | kfSampleCountsColMissingWithFemaleY)) != 0);
        if (male_ct || y_nonmale_needed) {
          chr_types |= 4;
          ++haploid_chr_type_ct;
        }
      }
      uint32_t mt_code;
      if (haploid_mask_lowbit || XymtExists(cip, kChrOffsetMT, &mt_code)) {
        chr_types |= 8;
        ++haploid_chr_type_ct;
      }
    }
    assert(diploid_chr_type_ct + haploid_chr_type_ct);
    // probable todo: for sample_ct > raw_sample_ct / 2 or so, may be better to
    // perform un-subsetted computation and then subset the results?  Benchmark
    // this.
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    if (chr_types & 6) {
      if (unlikely(
              bigstack_alloc_w(sample_ctl, &ctx.sex_male_collapsed))) {
        goto SampleCounts_ret_NOMEM;
      }
      CopyBitarrSubset(sex_male, sample_include, sample_ct, ctx.sex_male_collapsed);
    } else {
      ctx.sex_male_collapsed = nullptr;
    }

    // don't subtract 1 after load-balancing is improved?
    uint32_t calc_thread_ct = (max_thread_ct > 4)? (max_thread_ct - 1) : max_thread_ct;
    if (unlikely(
             bigstack_alloc_wp(calc_thread_ct, &ctx.raregenos) ||
             bigstack_alloc_u32p(calc_thread_ct, &ctx.difflist_sample_id_bufs) ||
             bigstack_alloc_vpp(calc_thread_ct, &ctx.thread_dense_counts) ||
             bigstack_alloc_u32pp(calc_thread_ct, &ctx.thread_sparse_counts) ||
             bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_sparse_common0_cts))) {
      goto SampleCounts_ret_NOMEM;
    }
    const uint32_t max_returned_difflist_len = 2 * (raw_sample_ct / kPglMaxDifflistLenDivisor);

    const uintptr_t raregeno_vec_ct = DivUp(max_returned_difflist_len, kNypsPerVec);
    const uintptr_t difflist_sample_id_vec_ct = DivUp(max_returned_difflist_len, kInt32PerVec);

    const uintptr_t acc2_vec_ct = NypCtToVecCt(sample_ct);
    const uintptr_t dense_counts_vstride = acc2_vec_ct * 23;
    // - need (0, 1, 2) genotype counts in diploid-chromosome case, only need
    //   (0, 1) for haploid
    const uintptr_t dense_counts_bottom_vec_ct = kSubstCodeCt * (dense_counts_vstride * (3 * diploid_chr_type_ct + 2 * haploid_chr_type_ct));
    const uintptr_t dense_counts_middle_vec_ct = DivUp(4 * kSubstCodeCt, kWordsPerVec);

    const uintptr_t sample_ct_i32v = DivUp(sample_ct, kInt32PerVec);
    // - x6 for 2 common_geno possibilities (0, 3), 3 common_geno -> rare_geno
    //   transitions per common_geno
    const uintptr_t sparse_counts_bottom_vec_ct = 6 * kSubstCodeCt * sample_ct_i32v * (diploid_chr_type_ct + haploid_chr_type_ct);
    const uintptr_t sample_ct_i32av = sample_ct_i32v * kInt32PerVec;
    // - x8 for 4 chromosome types, 2 common_geno possibilities
    const uintptr_t sparse_counts_middle_vec_ct = DivUp(8 * kSubstCodeCt, kWordsPerVec);
    const uintptr_t sparse_common0_vec_ct = DivUp(4 * kSubstCodeCt, kInt32PerVec);

    uintptr_t thread_xalloc_vec_ct = raregeno_vec_ct + difflist_sample_id_vec_ct + dense_counts_bottom_vec_ct + dense_counts_middle_vec_ct + sparse_counts_bottom_vec_ct + sparse_counts_middle_vec_ct + sparse_common0_vec_ct;

    ctx.thread_diploid_singleton_cts = nullptr;
    if (flags & kfSampleCountsColDiploidSingle) {
      if (unlikely(
              bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_diploid_singleton_cts))) {
        goto SampleCounts_ret_NOMEM;
      }
      thread_xalloc_vec_ct += sample_ct_i32v;
    }
    ctx.thread_singleton_cts = nullptr;
    if (flags & kfSampleCountsColSingle) {
      if (unlikely(
              bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_singleton_cts))) {
        goto SampleCounts_ret_NOMEM;
      }
      thread_xalloc_vec_ct += sample_ct_i32v;
    }
    // We support multiallelic variants by also tracking a few het_2alt,
    // hom_rarealt, and hap_rarealt counts.
    ctx.thread_alt_subst_codes = nullptr;
    ctx.thread_het2alt_cts = nullptr;
    ctx.thread_het_rarealt_cts = nullptr;
    ctx.thread_hom_rarealt_cts = nullptr;
    ctx.thread_hap_rarealt_cts = nullptr;
    const uint32_t mhc_needed = (max_allele_ct > 2);
    uintptr_t alt_subst_codes_vec_ct = 0;
    uintptr_t het2alt_vec_ct = 0;
    uintptr_t het_rarealt_vec_ct = 0;
    uintptr_t hom_rarealt_vec_ct = 0;
    uintptr_t hap_rarealt_vec_ct = 0;
    if (mhc_needed) {
      if (unlikely(
              bigstack_alloc_u16p(calc_thread_ct, &ctx.thread_alt_subst_codes) ||
              bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_het2alt_cts) ||
              bigstack_alloc_i32p(calc_thread_ct, &ctx.thread_het_rarealt_cts) ||
              bigstack_alloc_i32p(calc_thread_ct, &ctx.thread_hom_rarealt_cts))) {
        goto SampleCounts_ret_NOMEM;
      }
      alt_subst_codes_vec_ct = DivUp(kPglMaxAltAlleleCt + 1, kInt16PerVec);
      het2alt_vec_ct = sample_ct_i32v * 2 * kSubstCodeCt;
      het_rarealt_vec_ct = sample_ct_i32v * kSubstCodeCt;
      hom_rarealt_vec_ct = sample_ct_i32v * kSubstCodeCt;
      thread_xalloc_vec_ct += alt_subst_codes_vec_ct + het2alt_vec_ct + het_rarealt_vec_ct + hom_rarealt_vec_ct;
      if (chr_types & 14) {
        if (unlikely(
                bigstack_alloc_i32p(calc_thread_ct, &ctx.thread_hap_rarealt_cts))) {
          goto SampleCounts_ret_NOMEM;
        }
        hap_rarealt_vec_ct = sample_ct_i32v * kSubstCodeCt;
        thread_xalloc_vec_ct += hap_rarealt_vec_ct;
      }
    }

    const uintptr_t thread_xalloc_cacheline_ct = DivUp(thread_xalloc_vec_ct, kVecsPerCacheline);
    // Cacheline-alignment is relevant for per-thread workspaces.  We normally
    // get that automatically by allocating from the bottom of the arena, but
    // in this case it's more convenient to allocate them at the end, which is
    // normally only vector-aligned.  So we force the end to be
    // cacheline-aligned here.
    g_bigstack_end = R_CAST(unsigned char*, RoundDownPow2(R_CAST(uintptr_t, g_bigstack_end), kCacheline));

    unsigned char* bigstack_mark2 = g_bigstack_base;
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    ctx.thread_read_mhc = nullptr;
    uint32_t read_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, raw_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &ctx.genovecs, mhc_needed? (&ctx.thread_read_mhc) : nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto SampleCounts_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto SampleCounts_ret_NOMEM;
    }
    ctx.variant_include = variant_include;
    ctx.cip = cip;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.allele_storage = allele_storage;
    ctx.sample_include = sample_include;
    ctx.sample_ct = sample_ct;
    ctx.male_ct = male_ct;
    ctx.y_nonmale_needed = y_nonmale_needed;
    ctx.reterr = kPglRetSuccess;
    unsigned char* bigstack_end_mark2 = nullptr;
    assert(bigstack_left() >= thread_xalloc_cacheline_ct * kCacheline * calc_thread_ct);
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      unsigned char* cur_alloc = S_CAST(unsigned char*, bigstack_end_alloc_raw(thread_xalloc_cacheline_ct * kCacheline));
      ctx.raregenos[tidx] = R_CAST(uintptr_t*, cur_alloc);
      cur_alloc = &(cur_alloc[raregeno_vec_ct * kBytesPerVec]);
      ctx.difflist_sample_id_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
      cur_alloc = &(cur_alloc[difflist_sample_id_vec_ct * kBytesPerVec]);
      if (mhc_needed) {
        ctx.thread_alt_subst_codes[tidx] = R_CAST(uint16_t*, cur_alloc);
        cur_alloc = &(cur_alloc[alt_subst_codes_vec_ct * kBytesPerVec]);
      }
      if (!tidx) {
        bigstack_end_mark2 = cur_alloc;
      }
      ctx.thread_dense_counts[tidx] = R_CAST(VecW**, cur_alloc);
      cur_alloc = &(cur_alloc[dense_counts_middle_vec_ct * kBytesPerVec]);
      ZeroPtrArr(4 * kSubstCodeCt, ctx.thread_dense_counts[tidx]);
      for (uint32_t cur_chr_type = 0; cur_chr_type != 4; ++cur_chr_type) {
        if (!((chr_types >> cur_chr_type) & 1)) {
          continue;
        }
        const uintptr_t record_byte_ct = (dense_counts_vstride * (2 + (cur_chr_type < 2))) * kBytesPerVec;
        for (uint32_t subst_code = 0; subst_code != kSubstCodeCt; ++subst_code) {
          ctx.thread_dense_counts[tidx][cur_chr_type * kSubstCodeCt + subst_code] = R_CAST(VecW*, cur_alloc);
          cur_alloc = &(cur_alloc[record_byte_ct]);
        }
      }
      ctx.thread_sparse_counts[tidx] = R_CAST(uint32_t**, cur_alloc);
      cur_alloc = &(cur_alloc[sparse_counts_middle_vec_ct * kBytesPerVec]);
      ZeroPtrArr(8 * kSubstCodeCt, ctx.thread_sparse_counts[tidx]);
      for (uint32_t cur_chr_type = 0; cur_chr_type != 4; ++cur_chr_type) {
        if (!((chr_types >> cur_chr_type) & 1)) {
          continue;
        }
        const uintptr_t record_byte_ct = sample_ct_i32v * 3 * kBytesPerVec;
        for (uint32_t common_geno_and_subst_code = 0; common_geno_and_subst_code != 2 * kSubstCodeCt; ++common_geno_and_subst_code) {
          ctx.thread_sparse_counts[tidx][cur_chr_type * 2 * kSubstCodeCt + common_geno_and_subst_code] = R_CAST(uint32_t*, cur_alloc);
          cur_alloc = &(cur_alloc[record_byte_ct]);
        }
      }
      ctx.thread_sparse_common0_cts[tidx] = R_CAST(uint32_t*, cur_alloc);
      cur_alloc = &(cur_alloc[sparse_common0_vec_ct * kBytesPerVec]);

      if (ctx.thread_diploid_singleton_cts) {
        ctx.thread_diploid_singleton_cts[tidx] = R_CAST(uint32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[sample_ct_i32v * kBytesPerVec]);
      }
      if (ctx.thread_singleton_cts) {
        ctx.thread_singleton_cts[tidx] = R_CAST(uint32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[sample_ct_i32v * kBytesPerVec]);
      }
      if (mhc_needed) {
        ctx.thread_het2alt_cts[tidx] = R_CAST(uint32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[het2alt_vec_ct * kBytesPerVec]);
        ctx.thread_het_rarealt_cts[tidx] = R_CAST(int32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[het_rarealt_vec_ct * kBytesPerVec]);
        ctx.thread_hom_rarealt_cts[tidx] = R_CAST(int32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[hom_rarealt_vec_ct * kBytesPerVec]);
        if (ctx.thread_hap_rarealt_cts) {
          ctx.thread_hap_rarealt_cts[tidx] = R_CAST(int32_t*, cur_alloc);
          cur_alloc = &(cur_alloc[hap_rarealt_vec_ct * kBytesPerVec]);
        }
      }
      if (!tidx) {
        assert(S_CAST(uintptr_t, cur_alloc - g_bigstack_end) == thread_xalloc_vec_ct * kBytesPerVec);
      }
    }
    SetThreadFuncAndData(SampleCountsThread, &ctx, &tg);

    logputs("Calculating sample counts... ");
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;

    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_size = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto SampleCounts_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = ctx.reterr;
        if (unlikely(reterr)) {
          goto SampleCounts_ret_PGR_FAIL;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_size = cur_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_size, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_size == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto SampleCounts_ret_THREAD_CREATE_FAIL;
        }
      }

      parity = 1 - parity;
      if (variant_idx == variant_ct) {
        break;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }

      ++read_block_idx;
      variant_idx += cur_block_size;
      // crucially, this is independent of the PgenReader block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    const uintptr_t acc32_voffset = acc2_vec_ct * 7;
    const uintptr_t acc32_vec_ct = acc2_vec_ct * 16;
    for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
      // any reasonable compiler should vectorize the inner loops
      // could also multithread this (unscramble first, then split by sample
      // ranges), but computation time here should be negligible
      for (uint32_t dense_vtype = 0; dense_vtype != 4 * kSubstCodeCt; ++dense_vtype) {
        VecW* dst_acc2_0 = ctx.thread_dense_counts[0][dense_vtype];
        if (dst_acc2_0) {
          const VecW* src_acc2_0 = ctx.thread_dense_counts[tidx][dense_vtype];
          const VecU32* src_acc32_0 = R_CAST(const VecU32*, &(src_acc2_0[acc32_voffset]));
          VecU32* dst_acc32_0 = R_CAST(VecU32*, &(dst_acc2_0[acc32_voffset]));
          U32VecAdd(src_acc32_0, acc32_vec_ct, dst_acc32_0);
          const VecU32* src_acc32_2 = &(src_acc32_0[acc2_vec_ct * 23]);
          VecU32* dst_acc32_2 = &(dst_acc32_0[acc2_vec_ct * 23]);
          U32VecAdd(src_acc32_2, acc32_vec_ct, dst_acc32_2);
          if (dense_vtype < 2 * kSubstCodeCt) {
            const VecU32* src_acc32_1 = &(src_acc32_2[acc2_vec_ct * 23]);
            VecU32* dst_acc32_1 = &(dst_acc32_2[acc2_vec_ct * 23]);
            U32VecAdd(src_acc32_1, acc32_vec_ct, dst_acc32_1);
          }
        }
      }
      for (uint32_t sparse_vtype = 0; sparse_vtype != 8 * kSubstCodeCt; ++sparse_vtype) {
        uint32_t* dst = ctx.thread_sparse_counts[0][sparse_vtype];
        if (dst) {
          const uint32_t* src = ctx.thread_sparse_counts[tidx][sparse_vtype];
          U32CastVecAdd(src, sample_ct_i32v * 3, dst);
        }
      }
      {
        const uint32_t* src = ctx.thread_sparse_common0_cts[tidx];
        uint32_t* dst = ctx.thread_sparse_common0_cts[0];
        for (uint32_t dense_vtype = 0; dense_vtype != 4 * kSubstCodeCt; ++dense_vtype) {
          dst[dense_vtype] += src[dense_vtype];
        }
      }
      if (ctx.thread_diploid_singleton_cts) {
        const uint32_t* src = ctx.thread_diploid_singleton_cts[tidx];
        uint32_t* dst = ctx.thread_diploid_singleton_cts[0];
        U32CastVecAdd(src, sample_ct_i32v, dst);
      }
      if (ctx.thread_singleton_cts) {
        const uint32_t* src = ctx.thread_singleton_cts[tidx];
        uint32_t* dst = ctx.thread_singleton_cts[0];
        U32CastVecAdd(src, sample_ct_i32v, dst);
      }
      if (mhc_needed) {
        {
          uint32_t* dst = ctx.thread_het2alt_cts[0];
          const uint32_t* src = ctx.thread_het2alt_cts[tidx];
          U32CastVecAdd(src, sample_ct_i32v * 2 * kSubstCodeCt, dst);
        }

        int32_t* dst = ctx.thread_het_rarealt_cts[0];
        const int32_t* src = ctx.thread_het_rarealt_cts[tidx];
        I32CastVecAdd(src, sample_ct_i32v * kSubstCodeCt, dst);
        dst = ctx.thread_hom_rarealt_cts[0];
        src = ctx.thread_hom_rarealt_cts[tidx];
        I32CastVecAdd(src, sample_ct_i32v * kSubstCodeCt, dst);

        if (ctx.thread_hap_rarealt_cts) {
          dst = ctx.thread_hap_rarealt_cts[0];
          src = ctx.thread_hap_rarealt_cts[tidx];
          I32CastVecAdd(src, sample_ct_i32v * kSubstCodeCt, dst);
        }
      }
    }
    BigstackDoubleReset(bigstack_mark2, bigstack_end_mark2);
    // Now we're ready to derive our final values from the raw count-matrices.
    // 1. Unscramble the final dense raw counts.
    uint32_t* unscramble_buf;
    if (unlikely(
            bigstack_alloc_u32(acc32_vec_ct * kInt32PerVec, &unscramble_buf))) {
      goto SampleCounts_ret_NOMEM;
    }
    uint32_t* male_u32_mask = nullptr;
    uint32_t* haploid_counts[kSubstCodeCt][2];
    uint32_t* y_nonmale_counts[2];
    uint32_t* diploid_counts[kSubstCodeCt][3];
    {
      uint32_t* unscrambled_dense_counts[4 * kSubstCodeCt][3];
      for (uint32_t cur_chr_type = 0; cur_chr_type != 4; ++cur_chr_type) {
        const uint32_t vtype_stop = kSubstCodeCt * (cur_chr_type + 1);
        if (!((chr_types >> cur_chr_type) & 1)) {
          for (uint32_t vtype = kSubstCodeCt * cur_chr_type; vtype != vtype_stop; ++vtype) {
            unscrambled_dense_counts[vtype][0] = nullptr;
            unscrambled_dense_counts[vtype][1] = nullptr;
            unscrambled_dense_counts[vtype][2] = nullptr;
          }
          continue;
        }
        for (uint32_t vtype = kSubstCodeCt * cur_chr_type; vtype != vtype_stop; ++vtype) {
          VecW* dst_acc2_0 = ctx.thread_dense_counts[0][vtype];
          uint32_t* dst_0 = R_CAST(uint32_t*, &(dst_acc2_0[acc32_voffset]));
          Unscramble2(sample_ct, dst_0, unscramble_buf);
          unscrambled_dense_counts[vtype][0] = dst_0;

          uint32_t* dst_2 = &(dst_0[acc2_vec_ct * 23 * kInt32PerVec]);
          Unscramble2(sample_ct, dst_2, unscramble_buf);
          unscrambled_dense_counts[vtype][2] = dst_2;

          if (cur_chr_type < 2) {
            uint32_t* dst_1 = &(dst_2[acc2_vec_ct * 23 * kInt32PerVec]);
            Unscramble2(sample_ct, dst_1, unscramble_buf);
            unscrambled_dense_counts[vtype][1] = dst_1;
          } else {
            unscrambled_dense_counts[vtype][1] = nullptr;
          }
        }
      }

      BigstackReset(bigstack_mark2);
      // 2. Patch in the sparse counts.
      for (uint32_t cur_chr_type = 0; cur_chr_type != 4; ++cur_chr_type) {
        if (!((chr_types >> cur_chr_type) & 1)) {
          continue;
        }
        for (uint32_t subst_code = 0; subst_code != kSubstCodeCt; ++subst_code) {
          const uint32_t dense_vtype = kSubstCodeCt * cur_chr_type + subst_code;
          uint32_t* dense0 = unscrambled_dense_counts[dense_vtype][0];
          uint32_t* dense1 = unscrambled_dense_counts[dense_vtype][1];
          uint32_t* dense2 = unscrambled_dense_counts[dense_vtype][2];

          const uint32_t base0 = ctx.thread_sparse_common0_cts[0][dense_vtype];
          uint32_t sparse_vtype = dense_vtype + kSubstCodeCt * cur_chr_type;
          // 0->1
          const uint32_t* src = ctx.thread_sparse_counts[0][sparse_vtype];
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            dense0[sample_idx] += base0 - src[sample_idx];
          }
          if (dense1) {
            U32CastVecAdd(src, sample_ct_i32v, dense1);
          }
          // 0->2
          src = &(src[sample_ct_i32av]);
          U32CastVecSub(src, sample_ct_i32v, dense0);
          U32CastVecAdd(src, sample_ct_i32v, dense2);
          // 0->3
          src = &(src[sample_ct_i32av]);
          U32CastVecSub(src, sample_ct_i32v, dense0);

          sparse_vtype += kSubstCodeCt;
          src = ctx.thread_sparse_counts[0][sparse_vtype];
          // 3->2
          U32CastVecAdd(src, sample_ct_i32v, dense2);
          // 3->1
          src = &(src[sample_ct_i32av]);
          if (dense1) {
            U32CastVecAdd(src, sample_ct_i32v, dense1);
          }
          // 3->0
          src = &(src[sample_ct_i32av]);
          U32CastVecAdd(src, sample_ct_i32v, dense0);
        }
      }

      // 3. Allocate the reported-counts arrays, and fill them in a sensible
      //    order.

      // We start with
      //   haploid_counts := chrX/chrY male + generic-haploid
      //   y_nonmale_counts := chrY nonmale 0+2
      //   diploid_counts := autosomal-diploid + chrX nonmale
      // It's safe to move and alter the unscrambled_dense_counts arrays in the
      // process, since we only refer to haploid_counts/y_nonmale_counts/
      // diploid_counts afterward, but we have to be careful about timing when
      // we do that.
      for (uint32_t alt_ct = 0; alt_ct != 2; ++alt_ct) {
        y_nonmale_counts[alt_ct] = nullptr;
        for (uint32_t subst_code = 0; subst_code != kSubstCodeCt; ++subst_code) {
          haploid_counts[subst_code][alt_ct] = nullptr;
        }
      }
      if (chr_types & 14) {
        if (chr_types & 6) {
          // possible todo: not important here, but should be interesting to
          // benchmark vector-based ways of doing this (including a simple
          // 16- or 256-element lookup table)
          if (unlikely(
                  bigstack_end_calloc_u32(sample_ct_i32av, &male_u32_mask))) {
            goto SampleCounts_ret_NOMEM;
          }
          const uintptr_t* sex_male_collapsed = ctx.sex_male_collapsed;
          uintptr_t sample_idx_base = 0;
          uintptr_t cur_bits = sex_male_collapsed[0];
          for (uint32_t male_idx = 0; male_idx != male_ct; ++male_idx) {
            const uint32_t sample_idx = BitIter1(sex_male_collapsed, &sample_idx_base, &cur_bits);
            male_u32_mask[sample_idx] = UINT32_MAX;
          }
        }
        uint32_t* female_y_dst = nullptr;
        for (uint32_t alt_ct = 0; alt_ct != 2; ++alt_ct) {
          if (y_nonmale_needed) {
            if (unlikely(bigstack_end_calloc_u32(sample_ct_i32av, &female_y_dst))) {
              goto SampleCounts_ret_NOMEM;
            }
            y_nonmale_counts[alt_ct] = female_y_dst;
          }
          for (uint32_t subst_code = 0; subst_code != kSubstCodeCt; ++subst_code) {
            uint32_t* dst;
            if (chr_types & 8) {
              dst = unscrambled_dense_counts[3 * kSubstCodeCt + subst_code][alt_ct * 2];
            } else {
              if (unlikely(
                      bigstack_end_calloc_u32(sample_ct_i32av, &dst))) {
                goto SampleCounts_ret_NOMEM;
              }
            }
            haploid_counts[subst_code][alt_ct] = dst;
            if (chr_types & 6) {
              if (chr_types & 4) {
                const uint32_t* src2 = unscrambled_dense_counts[2 * kSubstCodeCt + subst_code][alt_ct * 2];
                if (female_y_dst) {
                  U32CastVecInvmaskedAdd(male_u32_mask, src2, sample_ct_i32v, female_y_dst);
                }
                if (chr_types & 2) {
                  const uint32_t* src1 = unscrambled_dense_counts[kSubstCodeCt + subst_code][alt_ct * 2];
                  U32CastVecMaskedAdd2(male_u32_mask, src1, src2, sample_ct_i32v, dst);
                } else {
                  U32CastVecMaskedAdd(male_u32_mask, src2, sample_ct_i32v, dst);
                }
              } else {
                const uint32_t* src1 = unscrambled_dense_counts[kSubstCodeCt + subst_code][alt_ct * 2];
                U32CastVecMaskedAdd(male_u32_mask, src1, sample_ct_i32v, dst);
              }
            }
          }
        }
      }
      for (uint32_t subst_code = 0; subst_code != kSubstCodeCt; ++subst_code) {
        for (uint32_t alt_ct = 0; alt_ct != 3; ++alt_ct) {
          uint32_t* dst;
          if (chr_types & 1) {
            dst = unscrambled_dense_counts[subst_code][alt_ct];
          } else {
            if (unlikely(
                    bigstack_end_calloc_u32(sample_ct_i32av, &dst))) {
              goto SampleCounts_ret_NOMEM;
            }
          }

          diploid_counts[subst_code][alt_ct] = dst;
          if (chr_types & 2) {
            // chrX nonmales
            const uint32_t* src = unscrambled_dense_counts[subst_code + kSubstCodeCt][alt_ct];
            U32CastVecInvmaskedAdd(male_u32_mask, src, sample_ct_i32v, dst);
          }
        }
      }
      if (mhc_needed) {
        // het_rarealt: modifies diploid alt_ct=1 subst_code distribution
        // hom_rarealt: modifies diploid alt_ct=2 subst_code distribution
        // hap_rarealt: modifies haploid alt_ct=1 subst_code distribution
        for (uint32_t subst_code = 0; subst_code != kSubstCodeCt; ++subst_code) {
          // let's make sure overflow behavior is defined, so if there's a bug
          // earlier on it's less likely to be a nightmare
          const uint32_t* src = R_CAST(uint32_t*, &(ctx.thread_het_rarealt_cts[0][subst_code * sample_ct_i32av]));
          uint32_t* dst = diploid_counts[subst_code][1];
          U32CastVecAdd(src, sample_ct_i32v, dst);
        }

        for (uint32_t subst_code = 0; subst_code != kSubstCodeCt; ++subst_code) {
          const uint32_t* src = R_CAST(uint32_t*, &(ctx.thread_hom_rarealt_cts[0][subst_code * sample_ct_i32av]));
          uint32_t* dst = diploid_counts[subst_code][2];
          U32CastVecAdd(src, sample_ct_i32v, dst);
        }

        if (ctx.thread_hap_rarealt_cts) {
          for (uint32_t subst_code = 0; subst_code != kSubstCodeCt; ++subst_code) {
            const uint32_t* src = R_CAST(uint32_t*, &(ctx.thread_hap_rarealt_cts[0][subst_code * sample_ct_i32av]));
            uint32_t* dst = haploid_counts[subst_code][1];
            U32CastVecAdd(src, sample_ct_i32v, dst);
          }
        }

        // het2alt: handled separately
      }
    }
    // Main dependency graph:
    //   het2alt ALT1= pieces
    //            |
    //            v
    //   het2alt reported sum + refalt -> het
    //                   |                 |
    //                   |                  ---> missing
    //                   v                 |
    //                 homalt + homref -> hom
    uint32_t* het2alt_reported_vals = nullptr;
    if (mhc_needed) {
      if (unlikely(
              bigstack_alloc_u32(sample_ct_i32av, &het2alt_reported_vals))) {
        goto SampleCounts_ret_NOMEM;
      }
      const uint32_t* src0 = ctx.thread_het2alt_cts[0];
      const uint32_t* src1 = &(src0[sample_ct_i32av]);
      const uint32_t* src2 = &(src1[sample_ct_i32av]);
      const uint32_t* src3 = &(src2[sample_ct_i32av]);
      const uint32_t* src4 = &(src3[sample_ct_i32av]);
      U32CastVecAssignAdd5(src0, src1, src2, src3, src4, sample_ct_i32v, het2alt_reported_vals);
    }
    uint32_t* refalt_vals;
    uint32_t* homalt_vals;
    uint32_t* homref_vals;
    uint32_t* hom_vals;
    if (unlikely(
            bigstack_alloc_u32(sample_ct_i32av, &refalt_vals) ||
            bigstack_alloc_u32(sample_ct_i32av, &homalt_vals) ||
            bigstack_alloc_u32(sample_ct_i32av, &homref_vals) ||
            bigstack_alloc_u32(sample_ct_i32av, &hom_vals))) {
      goto SampleCounts_ret_NOMEM;
    }
    {
      const uint32_t* src0 = diploid_counts[0][1];
      const uint32_t* src1 = diploid_counts[1][1];
      const uint32_t* src2 = diploid_counts[2][1];
      const uint32_t* src3 = diploid_counts[3][1];
      const uint32_t* src4 = diploid_counts[4][1];
      U32CastVecAssignAdd5(src0, src1, src2, src3, src4, sample_ct_i32v, refalt_vals);
    }
    uint32_t* het_vals;
    if (!mhc_needed) {
      het_vals = refalt_vals;
    } else {
      if (unlikely(
              bigstack_alloc_u32(sample_ct_i32av, &het_vals))) {
        goto SampleCounts_ret_NOMEM;
      }
      U32CastVecAssignAdd2(het2alt_reported_vals, refalt_vals, sample_ct_i32v, het_vals);
    }
    {
      const uint32_t* src0 = diploid_counts[0][2];
      const uint32_t* src1 = diploid_counts[1][2];
      const uint32_t* src2 = diploid_counts[2][2];
      const uint32_t* src3 = diploid_counts[3][2];
      const uint32_t* src4 = diploid_counts[4][2];
      U32CastVecAssignAdd5(src0, src1, src2, src3, src4, sample_ct_i32v, homalt_vals);
      if (het2alt_reported_vals) {
        U32CastVecSub(het2alt_reported_vals, sample_ct_i32v, homalt_vals);
      }

      src0 = diploid_counts[0][0];
      src1 = diploid_counts[1][0];
      src2 = diploid_counts[2][0];
      src3 = diploid_counts[3][0];
      src4 = diploid_counts[4][0];
      U32CastVecAssignAdd5(src0, src1, src2, src3, src4, sample_ct_i32v, homref_vals);
      U32CastVecAssignAdd2(homref_vals, homalt_vals, sample_ct_i32v, hom_vals);
    }
    uint32_t* hapref_vals;
    uint32_t* hapalt_vals;
    if (unlikely(
            bigstack_alloc_u32(sample_ct_i32av, &hapref_vals) ||
            bigstack_alloc_u32(sample_ct_i32av, &hapalt_vals))) {
      goto SampleCounts_ret_NOMEM;
    }
    {
      if (chr_types & 14) {
        const uint32_t* src0 = haploid_counts[0][0];
        const uint32_t* src1 = haploid_counts[1][0];
        const uint32_t* src2 = haploid_counts[2][0];
        const uint32_t* src3 = haploid_counts[3][0];
        const uint32_t* src4 = haploid_counts[4][0];
        U32CastVecAssignAdd5(src0, src1, src2, src3, src4, sample_ct_i32v, hapref_vals);
        src0 = haploid_counts[0][1];
        src1 = haploid_counts[1][1];
        src2 = haploid_counts[2][1];
        src3 = haploid_counts[3][1];
        src4 = haploid_counts[4][1];
        U32CastVecAssignAdd5(src0, src1, src2, src3, src4, sample_ct_i32v, hapalt_vals);
      } else {
        ZeroU32Arr(sample_ct, hapref_vals);
        ZeroU32Arr(sample_ct, hapalt_vals);
      }
    }

    uint32_t* final_counts[kSampleCountTypeCt];
    ZeroPtrArr(kSampleCountTypeCt, final_counts);
    if (flags & kfSampleCountsColHom) {
      final_counts[kSampleCountHom] = hom_vals;
    }
    if (flags & kfSampleCountsColHomref) {
      final_counts[kSampleCountHomref] = homref_vals;
    }
    if (flags & kfSampleCountsColHomalt) {
      final_counts[kSampleCountHomalt] = homalt_vals;
    }
    if (flags & kfSampleCountsColHomaltSnp) {
      uint32_t* dst;
      if (unlikely(
              bigstack_alloc_u32(sample_ct_i32av, &dst))) {
        goto SampleCounts_ret_NOMEM;
      }
      final_counts[kSampleCountHomaltSnp] = dst;
      const uint32_t* src1 = diploid_counts[kSubstCodeTs][2];
      const uint32_t* src2 = diploid_counts[kSubstCodeTv][2];
      const uint32_t* src3 = diploid_counts[kSubstCodeWeirdSnp][2];
      U32CastVecAssignAdd3(src1, src2, src3, sample_ct_i32v, dst);
      if (mhc_needed) {
        const uint32_t* src_base = ctx.thread_het2alt_cts[0];
        src1 = &(src_base[sample_ct_i32av * kSubstCodeTs]);
        src2 = &(src_base[sample_ct_i32av * kSubstCodeTv]);
        src3 = &(src_base[sample_ct_i32av * kSubstCodeWeirdSnp]);
        U32CastVecSub3(src1, src2, src3, sample_ct_i32v, dst);
      }
    }
    if (flags & kfSampleCountsColHet) {
      final_counts[kSampleCountHet] = het_vals;
    }
    if (flags & kfSampleCountsColRefalt) {
      final_counts[kSampleCountRefalt] = refalt_vals;
    }
    if (flags & kfSampleCountsColHet2alt) {
      if (het2alt_reported_vals) {
        final_counts[kSampleCountHet2alt] = het2alt_reported_vals;
      } else {
        if (unlikely(
                bigstack_calloc_u32(sample_ct_i32av, &(final_counts[kSampleCountHet2alt])))) {
          goto SampleCounts_ret_NOMEM;
        }
      }
    }

    if (flags & kfSampleCountsColHetSnp) {
      uint32_t* dst;
      if (unlikely(
              bigstack_alloc_u32(sample_ct_i32av, &dst))) {
        goto SampleCounts_ret_NOMEM;
      }
      final_counts[kSampleCountHetSnp] = dst;
      const uint32_t* src1 = diploid_counts[kSubstCodeTs][1];
      const uint32_t* src2 = diploid_counts[kSubstCodeTv][1];
      const uint32_t* src3 = diploid_counts[kSubstCodeWeirdSnp][1];
      U32CastVecAssignAdd3(src1, src2, src3, sample_ct_i32v, dst);
      if (mhc_needed) {
        const uint32_t* src_base = &(ctx.thread_het2alt_cts[0][kSubstCodeCt * sample_ct_i32av]);
        src1 = &(src_base[kSubstCodeTs * sample_ct_i32av]);
        src2 = &(src_base[kSubstCodeTv * sample_ct_i32av]);
        src3 = &(src_base[kSubstCodeWeirdSnp * sample_ct_i32av]);
        U32CastVecAdd3(src1, src2, src3, sample_ct_i32v, dst);
      }
    }
    if (flags & kfSampleCountsColDiploidTs) {
      final_counts[kSampleCountDiploidTs] = AllocAndCountSubstType(&ctx, TO_CONSTU32PCONSTP(diploid_counts[kSubstCodeTs]), nullptr, sample_ct_i32v, kSubstCodeTs);
      if (unlikely(!final_counts[kSampleCountDiploidTs])) {
        goto SampleCounts_ret_NOMEM;
      }
    }
    if (flags & kfSampleCountsColTs) {
      final_counts[kSampleCountTs] = AllocAndCountSubstType(&ctx, TO_CONSTU32PCONSTP(diploid_counts[kSubstCodeTs]), TO_CONSTU32PCONSTP(haploid_counts[kSubstCodeTs]), sample_ct_i32v, kSubstCodeTs);
      if (unlikely(!final_counts[kSampleCountTs])) {
        goto SampleCounts_ret_NOMEM;
      }
    }
    if (flags & kfSampleCountsColDiploidTv) {
      final_counts[kSampleCountDiploidTv] = AllocAndCountSubstType(&ctx, TO_CONSTU32PCONSTP(diploid_counts[kSubstCodeTv]), nullptr, sample_ct_i32v, kSubstCodeTv);
      if (unlikely(!final_counts[kSampleCountDiploidTv])) {
        goto SampleCounts_ret_NOMEM;
      }
    }
    if (flags & kfSampleCountsColTv) {
      final_counts[kSampleCountTv] = AllocAndCountSubstType(&ctx, TO_CONSTU32PCONSTP(diploid_counts[kSubstCodeTv]), TO_CONSTU32PCONSTP(haploid_counts[kSubstCodeTv]), sample_ct_i32v, kSubstCodeTv);
      if (unlikely(!final_counts[kSampleCountTv])) {
        goto SampleCounts_ret_NOMEM;
      }
    }
    if (flags & kfSampleCountsColDiploidNonsnpsymb) {
      final_counts[kSampleCountDiploidNonsnpsymb] = AllocAndCountSubstType(&ctx, TO_CONSTU32PCONSTP(diploid_counts[kSubstCodeNonsnpsymb]), nullptr, sample_ct_i32v, kSubstCodeNonsnpsymb);
      if (unlikely(!final_counts[kSampleCountDiploidNonsnpsymb])) {
        goto SampleCounts_ret_NOMEM;
      }
    }
    if (flags & kfSampleCountsColNonsnpsymb) {
      final_counts[kSampleCountNonsnpsymb] = AllocAndCountSubstType(&ctx, TO_CONSTU32PCONSTP(diploid_counts[kSubstCodeNonsnpsymb]), TO_CONSTU32PCONSTP(haploid_counts[kSubstCodeNonsnpsymb]), sample_ct_i32v, kSubstCodeNonsnpsymb);
      if (unlikely(!final_counts[kSampleCountNonsnpsymb])) {
        goto SampleCounts_ret_NOMEM;
      }
    }
    if (flags & kfSampleCountsColSymbolic) {
      final_counts[kSampleCountSymbolic] = AllocAndCountSubstType(&ctx, TO_CONSTU32PCONSTP(diploid_counts[kSubstCodeSymbolic]), TO_CONSTU32PCONSTP(haploid_counts[kSubstCodeSymbolic]), sample_ct_i32v, kSubstCodeSymbolic);
      if (unlikely(!final_counts[kSampleCountSymbolic])) {
        goto SampleCounts_ret_NOMEM;
      }
    }
    if (flags & kfSampleCountsColNonsnp) {
      final_counts[kSampleCountNonsnp] = AllocAndCountSubstType(&ctx, TO_CONSTU32PCONSTP(diploid_counts[kSubstCodeNonsnpsymb]), TO_CONSTU32PCONSTP(haploid_counts[kSubstCodeNonsnpsymb]), sample_ct_i32v, kSubstCodeNonsnpsymb);
      if (unlikely(!final_counts[kSampleCountNonsnp])) {
        goto SampleCounts_ret_NOMEM;
      }
      IncrSubstType(&ctx, TO_CONSTU32PCONSTP(diploid_counts[kSubstCodeSymbolic]), TO_CONSTU32PCONSTP(haploid_counts[kSubstCodeSymbolic]), sample_ct_i32v, kSubstCodeSymbolic, final_counts[kSampleCountNonsnp]);
    }

    if (flags & kfSampleCountsColDiploidSingle) {
      final_counts[kSampleCountDiploidSingle] = ctx.thread_diploid_singleton_cts[0];
    }
    if (flags & kfSampleCountsColSingle) {
      final_counts[kSampleCountSingle] = ctx.thread_singleton_cts[0];
    }
    if (flags & kfSampleCountsColHaprefWithFemaleY) {
      uint32_t* dst;
      if (unlikely(
              bigstack_alloc_u32(sample_ct_i32av, &dst))) {
        goto SampleCounts_ret_NOMEM;
      }
      final_counts[kSampleCountHaprefWithFemaleY] = dst;
      memcpy(dst, hapref_vals, sample_ct * sizeof(int32_t));
      if (chr_types & 4) {
        const uint32_t* src = y_nonmale_counts[0];
        U32CastVecAdd(src, sample_ct_i32v, dst);
      }
    }
    if (flags & kfSampleCountsColHapref) {
      final_counts[kSampleCountHapref] = hapref_vals;
    }
    if (flags & kfSampleCountsColHapaltWithFemaleY) {
      uint32_t* dst;
      if (unlikely(
              bigstack_alloc_u32(sample_ct_i32av, &dst))) {
        goto SampleCounts_ret_NOMEM;
      }
      final_counts[kSampleCountHapaltWithFemaleY] = dst;
      memcpy(dst, hapalt_vals, sample_ct * sizeof(int32_t));
      if (chr_types & 4) {
        const uint32_t* src = y_nonmale_counts[1];
        U32CastVecAdd(src, sample_ct_i32v, dst);
      }
    }
    if (flags & kfSampleCountsColHapalt) {
      final_counts[kSampleCountHapalt] = hapalt_vals;
    }
    if (flags & kfSampleCountsColMissingWithFemaleY) {
      uint32_t* dst;
      if (unlikely(
              bigstack_alloc_u32(sample_ct_i32av, &dst))) {
        goto SampleCounts_ret_NOMEM;
      }
      final_counts[kSampleCountMissingWithFemaleY] = dst;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        // todo: verify this vectorizes well
        dst[sample_idx] = variant_ct - hom_vals[sample_idx] - het_vals[sample_idx] - hapref_vals[sample_idx] - hapalt_vals[sample_idx];
      }
      if (chr_types & 4) {
        const uint32_t* src0 = y_nonmale_counts[0];
        const uint32_t* src1 = y_nonmale_counts[1];
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          dst[sample_idx] -= src0[sample_idx] + src1[sample_idx];
        }
      }
    }
    if (flags & kfSampleCountsColMissing) {
      uint32_t* dst;
      if (unlikely(
              bigstack_alloc_u32(sample_ct_i32av, &dst))) {
        goto SampleCounts_ret_NOMEM;
      }
      final_counts[kSampleCountMissing] = dst;
      const uint32_t nony_variant_ct = variant_ct - y_variant_ct;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        dst[sample_idx] = nony_variant_ct - hom_vals[sample_idx] - het_vals[sample_idx] - hapref_vals[sample_idx] - hapalt_vals[sample_idx];
      }
      if (chr_types & 4) {
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          dst[sample_idx] += male_u32_mask[sample_idx] & y_variant_ct;
        }
      }
    }

    // Ready to write to disk.
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    BigstackEndReset(bigstack_end_mark2);
    const uintptr_t overflow_buf_size = RoundUpPow2(kCompressStreamBlock + 3 * kMaxIdSlen + 512, kCacheline);
    const uint32_t output_zst = flags & kfSampleCountsZs;
    uintptr_t overflow_buf_alloc = overflow_buf_size;
    if (output_zst) {
      overflow_buf_alloc += CstreamWkspaceReq(overflow_buf_size);
    }
    char* overflow_buf;
    if (unlikely(
            bigstack_alloc_c(overflow_buf_alloc, &overflow_buf))) {
      goto SampleCounts_ret_NOMEM;
    }
    OutnameZstSet(".scount", output_zst, outname_end);
    reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
    if (unlikely(reterr)) {
      goto SampleCounts_ret_1;
    }
    cswritep = overflow_buf;
    *cswritep++ = '#';
    const uint32_t col_fid = FidColIsRequired(siip, flags / kfSampleCountsColMaybefid);
    if (col_fid) {
      cswritep = strcpya_k(cswritep, "FID\t");
    }
    cswritep = strcpya_k(cswritep, "IID");
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    const uint32_t col_sid = SidColIsRequired(sids, flags / kfSampleCountsColMaybesid);
    if (col_sid) {
      cswritep = strcpya_k(cswritep, "\tSID");
    }
    const uint32_t col_sex = (flags / kfSampleCountsColSex) & 1;
    if (col_sex) {
      cswritep = strcpya_k(cswritep, "\tSEX");
    }
    uint32_t type_bitarr;
    memcpy(&type_bitarr, &flags, 4);
    type_bitarr = type_bitarr / kfSampleCountsColHom;
    {
      uint32_t types_left = type_bitarr;
      while (types_left) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, sample_counts_headers[ctzu32(types_left)]);
        types_left = types_left & (types_left - 1);
      }
    }
    AppendBinaryEoln(&cswritep);
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
      if (!col_fid) {
        cur_sample_id = AdvPastDelim(cur_sample_id, '\t');
      }
      cswritep = strcpya(cswritep, cur_sample_id);
      if (col_sid) {
        *cswritep++ = '\t';
        if (sids) {
          cswritep = strcpya(cswritep, &(sids[sample_uidx * max_sid_blen]));
        } else {
          *cswritep++ = '0';
        }
      }
      if (col_sex) {
        *cswritep++ = '\t';
        if (IsSet(sex_nm, sample_uidx)) {
          *cswritep++ = '2' - IsSet(sex_male, sample_uidx);
        } else {
          cswritep = strcpya_k(cswritep, "NA");
        }
      }
      uint32_t types_left = type_bitarr;
      while (types_left) {
        const uint32_t type_idx = ctzu32(types_left);
        *cswritep++ = '\t';
        cswritep = u32toa(final_counts[type_idx][sample_idx], cswritep);
        types_left = types_left & (types_left - 1);
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto SampleCounts_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto SampleCounts_ret_WRITE_FAIL;
    }
    logprintfww("--sample-counts: Results written to %s .\n", outname);
  }
  while (0) {
  SampleCounts_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SampleCounts_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  SampleCounts_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  SampleCounts_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 SampleCounts_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&css, cswritep);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

static const Dosage kGenoToDosage[4] = {0, 16384, 32768, 65535};

typedef struct SdiffCountsStruct {
  uint32_t missing_ct; // never includes halfmiss
  uint32_t ibsmiss_ct; // can determine ibs2_ct from this
  uint32_t ibsx_cts[2];
  uint32_t halfmiss_ct;
  uint32_t diff_ct; // never includes halfmiss
} SdiffCounts;

// Simplest case.  May multithread this, add biallelic hardcall and other
// optimizations, etc. later.
PglErr SdiffCountsOnly(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts, const uint32_t* __restrict id_pairs, const uintptr_t* __restrict pair_sex_male, const uintptr_t* __restrict variant_include, const ChrInfo* cip, const uintptr_t* __restrict allele_idx_offsets, const SdiffInfo* sdip, uint32_t sample_ct, uint32_t variant_ct, uintptr_t id_pair_ct, PgenReader* simple_pgrp, SdiffCounts* sdiff_counts) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t dosage_hap_tol = sdip->dosage_hap_tol;
    const uint32_t dosage_needed = (PgrGetGflags(simple_pgrp) & kfPgenGlobalDosagePresent) && (dosage_hap_tol != kDosageMissing);

    // values unimportant if dosage_hap_tol == kDosageMissing
    const uint32_t dosage_dip_tol = dosage_hap_tol / 2;
    uint32_t dosage_sex_tols[2]; // for chrX
    dosage_sex_tols[0] = dosage_dip_tol;
    dosage_sex_tols[1] = dosage_hap_tol;

    const uint32_t ibs_needed = ((sdip->flags & kfSdiffCountsIbsNeeded) != kfSdiff0);
    PgenVariant pgv;
    if (unlikely(BigstackAllocPgv(sample_ct, allele_idx_offsets != nullptr, dosage_needed? kfPgenGlobalDosagePresent : kfPgenGlobal0, &pgv))) {
      goto SdiffCountsOnly_ret_NOMEM;
    }
    Dosage* dosage_buf = nullptr;
    if (dosage_needed) {
      // todo: multidosage
      if (unlikely(bigstack_alloc_dosage(sample_ct, &dosage_buf))) {
        goto SdiffCountsOnly_ret_NOMEM;
      }
    }
    AlleleCode* allele_code_buf = nullptr;
    if (allele_idx_offsets != nullptr) {
      if (unlikely(bigstack_alloc_ac(2 * sample_ct, &allele_code_buf))) {
        goto SdiffCountsOnly_ret_NOMEM;
      }
    }
    pgv.patch_01_ct = 0;
    pgv.patch_10_ct = 0;
    pgv.dosage_ct = 0;
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t is_autosomal_diploid = 0;
    uint32_t is_x = 0;
    uint32_t is_y = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t allele_ct = 2;
    fputs("--sample-diff counts-only: 0%", stdout);
    fflush(stdout);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        is_autosomal_diploid = !IsSet(cip->haploid_mask, chr_idx);
        is_x = (chr_idx == x_code);
        is_y = (chr_idx == y_code);
      }
      if (allele_idx_offsets) {
        allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
      }
      if (allele_ct == 2) {
        if (!dosage_needed) {
          reterr = PgrGet(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec);
        } else {
          reterr = PgrGetD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
        }
      } else {
        if (!dosage_needed) {
          reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
        } else {
          reterr = PgrGetMD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
        }
        if ((!pgv.patch_01_ct) && (!pgv.patch_10_ct)) {
          // todo: also check for multidosage
          allele_ct = 2;
        }
      }
      if (unlikely(reterr)) {
        goto SdiffCountsOnly_ret_PGR_FAIL;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }
      if (allele_ct == 2) {
        if (!pgv.dosage_ct) {
          if (AllGenoEqual(pgv.genovec, sample_ct)) {
            // don't need to do anything!
            continue;
          }
          const uint32_t* id_pair_iter = id_pairs;
          if (is_autosomal_diploid) {
            for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
              const uint32_t sample_idx1 = *id_pair_iter++;
              const uint32_t sample_idx2 = *id_pair_iter++;
              const uintptr_t hc1 = GetNyparrEntry(pgv.genovec, sample_idx1);
              const uintptr_t hc2 = GetNyparrEntry(pgv.genovec, sample_idx2);
              if (hc1 == hc2) {
                if (hc1 == 3) {
                  sdiff_counts[pair_idx].missing_ct += 1;
                  sdiff_counts[pair_idx].ibsmiss_ct += 1;
                }
                continue;
              }
              SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
              if ((hc1 == 3) || (hc2 == 3)) {
                cur_counts->halfmiss_ct += 1;
                cur_counts->ibsmiss_ct += 1;
                continue;
              }
              cur_counts->diff_ct += 1;
              cur_counts->ibsx_cts[(hc1 | hc2) & 1] += 1;
            }
          } else if (is_x) {
            for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
              const uint32_t sample_idx1 = *id_pair_iter++;
              const uint32_t sample_idx2 = *id_pair_iter++;
              const uintptr_t hc1 = GetNyparrEntry(pgv.genovec, sample_idx1);
              const uintptr_t hc2 = GetNyparrEntry(pgv.genovec, sample_idx2);
              if (hc1 == hc2) {
                if (hc1 == 3) {
                  sdiff_counts[pair_idx].missing_ct += 1;
                  if (!IsSet(pair_sex_male, pair_idx)) {
                    sdiff_counts[pair_idx].ibsmiss_ct += 1;
                  }
                }
                continue;
              }
              SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
              const uint32_t is_male = IsSet(pair_sex_male, pair_idx);
              if ((hc1 == 3) || (hc2 == 3)) {
                cur_counts->halfmiss_ct += 1;
                cur_counts->ibsmiss_ct += !is_male;
                continue;
              }
              cur_counts->diff_ct += 1;
              if (!is_male) {
                cur_counts->ibsx_cts[(hc1 | hc2) & 1] += 1;
              }
            }
          } else {
            for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
              if (is_y && (!IsSet(pair_sex_male, pair_idx))) {
                continue;
              }
              const uint32_t sample_idx1 = id_pairs[2 * pair_idx];
              const uint32_t sample_idx2 = id_pairs[2 * pair_idx + 1];
              const uintptr_t hc1 = GetNyparrEntry(pgv.genovec, sample_idx1);
              const uintptr_t hc2 = GetNyparrEntry(pgv.genovec, sample_idx2);
              if (hc1 == hc2) {
                if (hc1 == 3) {
                  sdiff_counts[pair_idx].missing_ct += 1;
                }
                continue;
              }
              SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
              if ((hc1 == 3) || (hc2 == 3)) {
                cur_counts->halfmiss_ct += 1;
                continue;
              }
              cur_counts->diff_ct += 1;
            }
          }
        } else {
          // dosages present
          const Dosage* dosage_read_iter = pgv.dosage_main;
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            if (!IsSet(pgv.dosage_present, sample_idx)) {
              dosage_buf[sample_idx] = kGenoToDosage[GetNyparrEntry(pgv.genovec, sample_idx)];
            } else {
              dosage_buf[sample_idx] = *dosage_read_iter++;
            }
          }
          const uint32_t* id_pair_iter = id_pairs;
          if (is_autosomal_diploid) {
            for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
              const uint32_t sample_idx1 = *id_pair_iter++;
              const uint32_t sample_idx2 = *id_pair_iter++;
              const uint32_t dosage1 = dosage_buf[sample_idx1];
              const uint32_t dosage2 = dosage_buf[sample_idx2];
              SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
              if (dosage1 == kDosageMissing) {
                cur_counts->ibsmiss_ct += 1;
                if (dosage2 == kDosageMissing) {
                  cur_counts->missing_ct += 1;
                } else {
                  cur_counts->halfmiss_ct += 1;
                }
                continue;
              }
              if (dosage2 == kDosageMissing) {
                cur_counts->halfmiss_ct += 1;
                cur_counts->ibsmiss_ct += 1;
                continue;
              }
              if (abs_i32(dosage1 - dosage2) > dosage_dip_tol) {
                cur_counts->diff_ct += 1;
              }
              if (!ibs_needed) {
                continue;
              }
              const uintptr_t hc1 = GetNyparrEntry(pgv.genovec, sample_idx1);
              const uintptr_t hc2 = GetNyparrEntry(pgv.genovec, sample_idx2);
              if (hc1 == hc2) {
                continue;
              }
              if ((hc1 == 3) || (hc2 == 3)) {
                cur_counts->ibsmiss_ct += 1;
              }
              cur_counts->ibsx_cts[(hc1 | hc2) & 1] += 1;
            }
          } else if (is_x) {
            for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
              const uint32_t sample_idx1 = *id_pair_iter++;
              const uint32_t sample_idx2 = *id_pair_iter++;
              const uint32_t dosage1 = dosage_buf[sample_idx1];
              const uint32_t dosage2 = dosage_buf[sample_idx2];
              SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
              if (dosage1 == kDosageMissing) {
                cur_counts->ibsmiss_ct += 1;
                if (dosage2 == kDosageMissing) {
                  cur_counts->missing_ct += 1;
                } else {
                  cur_counts->halfmiss_ct += 1;
                }
                continue;
              }
              if (dosage2 == kDosageMissing) {
                cur_counts->halfmiss_ct += 1;
                cur_counts->ibsmiss_ct += 1;
                continue;
              }
              const uint32_t is_male = IsSet(pair_sex_male, pair_idx);
              if (abs_i32(dosage1 - dosage2) > dosage_sex_tols[is_male]) {
                cur_counts->diff_ct += 1;
              }
              if ((!ibs_needed) || is_male) {
                continue;
              }
              const uintptr_t hc1 = GetNyparrEntry(pgv.genovec, sample_idx1);
              const uintptr_t hc2 = GetNyparrEntry(pgv.genovec, sample_idx2);
              if (hc1 == hc2) {
                continue;
              }
              if ((hc1 == 3) || (hc2 == 3)) {
                cur_counts->ibsmiss_ct += 1;
              }
              cur_counts->ibsx_cts[(hc1 | hc2) & 1] += 1;
            }
          } else {
            for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
              if (is_y && (!IsSet(pair_sex_male, pair_idx))) {
                continue;
              }
              const uint32_t sample_idx1 = *id_pair_iter++;
              const uint32_t sample_idx2 = *id_pair_iter++;
              const uint32_t dosage1 = dosage_buf[sample_idx1];
              const uint32_t dosage2 = dosage_buf[sample_idx2];
              SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
              if (dosage1 == kDosageMissing) {
                cur_counts->ibsmiss_ct += 1;
                if (dosage2 == kDosageMissing) {
                  cur_counts->missing_ct += 1;
                } else {
                  cur_counts->halfmiss_ct += 1;
                }
                continue;
              }
              if (dosage2 == kDosageMissing) {
                cur_counts->halfmiss_ct += 1;
                cur_counts->ibsmiss_ct += 1;
                continue;
              }
              if (abs_i32(dosage1 - dosage2) > dosage_hap_tol) {
                cur_counts->diff_ct += 1;
              }
            }
          }
        }
      } else {
        // multiallelic
        if (!pgv.dosage_ct) {
          PglMultiallelicSparseToDenseMiss(&pgv, sample_ct, allele_code_buf);
          const uint32_t* id_pair_iter = id_pairs;
          if (is_autosomal_diploid) {
            for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
              const uint32_t sample_idx1 = *id_pair_iter++;
              const uint32_t sample_idx2 = *id_pair_iter++;
              const AlleleCode ac11 = allele_code_buf[2 * sample_idx1];
              const AlleleCode ac12 = allele_code_buf[2 * sample_idx1 + 1];
              const AlleleCode ac21 = allele_code_buf[2 * sample_idx2];
              const AlleleCode ac22 = allele_code_buf[2 * sample_idx2 + 1];
              if ((ac11 == ac21) && (ac12 == ac22)) {
                if (ac11 == kMissingAlleleCode) {
                  sdiff_counts[pair_idx].missing_ct += 1;
                  sdiff_counts[pair_idx].ibsmiss_ct += 1;
                }
                continue;
              }
              SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
              if ((ac11 == kMissingAlleleCode) || (ac21 == kMissingAlleleCode)) {
                cur_counts->halfmiss_ct += 1;
                cur_counts->ibsmiss_ct += 1;
                continue;
              }
              cur_counts->diff_ct += 1;
              if (ibs_needed) {
                cur_counts->ibsx_cts[(ac11 == ac21) || (ac11 == ac22) || (ac12 == ac21) || (ac12 == ac22)] += 1;
              }
            }
          } else if (is_x) {
            for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
              const uint32_t sample_idx1 = *id_pair_iter++;
              const uint32_t sample_idx2 = *id_pair_iter++;
              const AlleleCode ac11 = allele_code_buf[2 * sample_idx1];
              const AlleleCode ac12 = allele_code_buf[2 * sample_idx1 + 1];
              const AlleleCode ac21 = allele_code_buf[2 * sample_idx2];
              const AlleleCode ac22 = allele_code_buf[2 * sample_idx2 + 1];
              if ((ac11 == ac21) && (ac12 == ac22)) {
                if (ac11 == kMissingAlleleCode) {
                  sdiff_counts[pair_idx].missing_ct += 1;
                  if (!IsSet(pair_sex_male, pair_idx)) {
                    sdiff_counts[pair_idx].ibsmiss_ct += 1;
                  }
                }
                continue;
              }
              SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
              const uint32_t is_male = IsSet(pair_sex_male, pair_idx);
              if ((ac11 == kMissingAlleleCode) || (ac21 == kMissingAlleleCode)) {
                cur_counts->halfmiss_ct += 1;
                cur_counts->ibsmiss_ct += !is_male;
                continue;
              }
              cur_counts->diff_ct += 1;
              if (ibs_needed && (!is_male)) {
                cur_counts->ibsx_cts[(ac11 == ac21) || (ac11 == ac22) || (ac12 == ac21) || (ac12 == ac22)] += 1;
              }
            }
          } else {
            for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
              if (is_y && (!IsSet(pair_sex_male, pair_idx))) {
                continue;
              }
              const uint32_t sample_idx1 = id_pairs[2 * pair_idx];
              const uint32_t sample_idx2 = id_pairs[2 * pair_idx + 1];
              const AlleleCode ac11 = allele_code_buf[2 * sample_idx1];
              const AlleleCode ac12 = allele_code_buf[2 * sample_idx1 + 1];
              const AlleleCode ac21 = allele_code_buf[2 * sample_idx2];
              const AlleleCode ac22 = allele_code_buf[2 * sample_idx2 + 1];
              if ((ac11 == ac21) && (ac12 == ac22)) {
                if (ac11 == kMissingAlleleCode) {
                  sdiff_counts[pair_idx].missing_ct += 1;
                }
                continue;
              }
              SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
              if ((ac11 == kMissingAlleleCode) || (ac21 == kMissingAlleleCode)) {
                cur_counts->halfmiss_ct += 1;
              }
              cur_counts->diff_ct += 1;
            }
          }
        } else {
          // multiallelic-dosage; todo
          exit(S_CAST(int32_t, kPglRetInternalError));
        }
      }
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
  }
  while (0) {
  SdiffCountsOnly_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SdiffCountsOnly_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

void AppendSdiffHeaderLine(SdiffFlags flags, uint32_t col_fid, uint32_t col_sid, uint32_t dosage_needed, char** cswritepp) {
  char* cswritep = *cswritepp;
  *cswritep++ = '#';
  if (flags & kfSdiffColChrom) {
    cswritep = strcpya_k(cswritep, "CHROM\t");
  }
  if (flags & kfSdiffColPos) {
    cswritep = strcpya_k(cswritep, "POS\t");
  }
  cswritep = strcpya_k(cswritep, "ID");
  if (flags & kfSdiffColRef) {
    cswritep = strcpya_k(cswritep, "\tREF");
  }
  if (flags & kfSdiffColAlt) {
    cswritep = strcpya_k(cswritep, "\tALT");
  }
  if (flags & kfSdiffColId) {
    if (col_fid) {
      cswritep = strcpya_k(cswritep, "\tFID1");
    }
    cswritep = strcpya_k(cswritep, "\tIID1");
    if (col_sid) {
      cswritep = strcpya_k(cswritep, "\tSID1");
    }
    if (col_fid) {
      cswritep = strcpya_k(cswritep, "\tFID2");
    }
    cswritep = strcpya_k(cswritep, "\tIID2");
    if (col_sid) {
      cswritep = strcpya_k(cswritep, "\tSID2");
    }
  }
  if (flags & kfSdiffColGeno) {
    if (!dosage_needed) {
      cswritep = strcpya_k(cswritep, "\tGT1\tGT2");
    } else {
      cswritep = strcpya_k(cswritep, "\tDS1\tDS2");
    }
  }
  *cswritepp = cswritep;
  AppendBinaryEoln(cswritepp);
}

typedef struct SdiffWriteContextStruct {
  char* collapsed_sample_fmtids;
  char* chr_buf;
  const uint32_t* variant_bps;
  const char* const* variant_ids;
  SdiffFlags flags;
  uintptr_t max_sample_fmtid_blen;
  uint32_t chr_buf_blen;
} SdiffWriteContext;

BoolErr AppendSdiffPregenoFields(const SdiffWriteContext* swcp, const char* const* cur_alleles, uint32_t sample_idx1, uint32_t sample_idx2, uint32_t variant_uidx, uint32_t allele_ct, CompressStreamState* cssp, char** cswritepp) {
  char* cswritep = *cswritepp;
  if (swcp->chr_buf) {
    cswritep = memcpya(cswritep, swcp->chr_buf, swcp->chr_buf_blen);
  }
  if (swcp->variant_bps) {
    cswritep = u32toa_x(swcp->variant_bps[variant_uidx], '\t', cswritep);
  }
  cswritep = strcpya(cswritep, swcp->variant_ids[variant_uidx]);
  if (swcp->flags & kfSdiffColRef) {
    *cswritep++ = '\t';
    const char* cur_allele = cur_alleles[0];
    const uint32_t allele_slen = strlen(cur_allele);
    if (unlikely(CsputsStd(cur_allele, allele_slen, cssp, &cswritep))) {
      // might not need this assignment, but play it safe for now
      *cswritepp = cswritep;
      return 1;
    }
  }
  if (swcp->flags & kfSdiffColAlt) {
    *cswritep++ = '\t';
    for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
      const char* cur_allele = cur_alleles[allele_idx];
      const uint32_t allele_slen = strlen(cur_allele);
      if (unlikely(CsputsStd(cur_allele, allele_slen, cssp, &cswritep))) {
        *cswritepp = cswritep;
        return 1;
      }
      *cswritep++ = ',';
    }
    --cswritep;
  }
  if (swcp->collapsed_sample_fmtids) {
    *cswritep++ = '\t';
    cswritep = strcpyax(cswritep, &(swcp->collapsed_sample_fmtids[sample_idx1 * swcp->max_sample_fmtid_blen]), '\t');
    cswritep = strcpya(cswritep, &(swcp->collapsed_sample_fmtids[sample_idx2 * swcp->max_sample_fmtid_blen]));
  }
  *cswritepp = cswritep;
  return 0;
}

PglErr SdiffMainBatch(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts,  const SampleIdInfo* siip, const uint32_t* __restrict id_pairs, const uintptr_t* __restrict pair_sex_male, const uintptr_t* __restrict variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* __restrict allele_idx_offsets, const char* const* allele_storage, const SdiffInfo* sdip, uint32_t sample_ct, uint32_t variant_ct, uintptr_t id_pair_ct, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end, SdiffCounts* sdiff_counts) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char** cswritep_arr = nullptr;
  CompressStreamState* css_arr = nullptr;
  const SdiffFlags flags = sdip->flags;
  const uint32_t is_pairwise = (flags / kfSdiffPairwise) & 1;
  const uintptr_t file_ct = is_pairwise? id_pair_ct : 1;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t col_fid = FidColIsRequired(siip, flags / kfSdiffColMaybefid);
    const uint32_t col_sid = SidColIsRequired(siip->sids, flags / kfSdiffColMaybesid);
    SdiffWriteContext swc;
    if (unlikely(CollapsedSampleFmtidInitAlloc(sample_include, siip, sample_ct, col_fid, col_sid, &swc.collapsed_sample_fmtids, &swc.max_sample_fmtid_blen))) {
      goto SdiffMainBatch_ret_NOMEM;
    }
    swc.chr_buf = nullptr;
    if (flags & kfSdiffColChrom) {
      const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
      if (unlikely(bigstack_alloc_c(max_chr_blen, &swc.chr_buf))) {
        goto SdiffMainBatch_ret_NOMEM;
      }
    }
    if (!(flags & kfSdiffColPos)) {
      swc.variant_bps = nullptr;
    } else {
      swc.variant_bps = variant_bps;
    }
    swc.variant_ids = variant_ids;
    if (unlikely(
            bigstack_calloc_cp(file_ct, &cswritep_arr) ||
            BIGSTACK_ALLOC_X(CompressStreamState, file_ct, &css_arr))) {
      goto SdiffMainBatch_ret_NOMEM;
    }
    swc.flags = flags;
    swc.chr_buf_blen = 0;
    for (uintptr_t fidx = 0; fidx != id_pair_ct; ++fidx) {
      PreinitCstream(&(css_arr[fidx]));
    }
    const uint32_t output_zst = (flags / kfSdiffZs) & 1;
    // ".sdiff" + terminating null = 7 bytes
    const uintptr_t fname_extrachar_limit = kPglFnamesize - 7 - (output_zst * 4) - S_CAST(uintptr_t, outname_end - outname);
    const char fname_id_delim = sdip->fname_id_delim;
    const uint32_t dosage_hap_tol = sdip->dosage_hap_tol;
    const uint32_t dosage_reported = (dosage_hap_tol != kDosageMissing);
    const uint32_t dosage_needed = (PgrGetGflags(simple_pgrp) & kfPgenGlobalDosagePresent) && dosage_reported;
    if (is_pairwise) {
      const uint32_t* id_pair_iter = id_pairs;
      const uint32_t stream_thread_ct = (max_thread_ct > (2 * id_pair_ct))? ((max_thread_ct - 1) / id_pair_ct) : 1;
      *outname_end = '.';
      for (uintptr_t fidx = 0; fidx != id_pair_ct; ++fidx) {
        const uint32_t sample_idx1 = *id_pair_iter++;
        const uint32_t sample_idx2 = *id_pair_iter++;
        const char* sample_fmtid1 = &(swc.collapsed_sample_fmtids[sample_idx1 * swc.max_sample_fmtid_blen]);
        const char* sample_fmtid2 = &(swc.collapsed_sample_fmtids[sample_idx2 * swc.max_sample_fmtid_blen]);
        const uint32_t slen1 = strlen(sample_fmtid1);
        const uint32_t slen2 = strlen(sample_fmtid2);
        if (unlikely(slen1 + slen2 + 2 > fname_extrachar_limit)) {
          logputs("\n");
          logerrputs("Error: Sample ID and/or --out argument too long for --sample-diff pairwise\nmode.\n");
          goto SdiffMainBatch_ret_INCONSISTENT_INPUT;
        }
        char* fname_iter = &(outname_end[1]);
        // note that this copies tab(s); replace them with fname_id_delim
        fname_iter = memcpya(fname_iter, sample_fmtid1, slen1);
        char* tab_iter = &(outname_end[2]);
        tab_iter = S_CAST(char*, memchr(tab_iter, '\t', fname_iter - tab_iter));
        if (tab_iter) {
          *tab_iter++ = fname_id_delim;
          tab_iter = S_CAST(char*, memchr(tab_iter, '\t', fname_iter - tab_iter));
          if (tab_iter) {
            *tab_iter = fname_id_delim;
          }
        }
        *fname_iter++ = '.';
        tab_iter = &(fname_iter[1]);
        fname_iter = memcpya(fname_iter, sample_fmtid2, slen2);
        tab_iter = S_CAST(char*, memchr(tab_iter, '\t', fname_iter - tab_iter));
        if (tab_iter) {
          *tab_iter++ = fname_id_delim;
          tab_iter = S_CAST(char*, memchr(tab_iter, '\t', fname_iter - tab_iter));
          if (tab_iter) {
            *tab_iter = fname_id_delim;
          }
        }
        fname_iter = strcpya_k(fname_iter, ".sdiff");
        if (output_zst) {
          fname_iter = strcpya_k(fname_iter, ".zst");
        }
        *fname_iter = '\0';
        reterr = InitCstreamAlloc(outname, 0, output_zst, stream_thread_ct, kCompressStreamBlock + kMaxMediumLine, &(css_arr[fidx]), &(cswritep_arr[fidx]));
        if (unlikely(reterr)) {
          goto SdiffMainBatch_ret_1;
        }
        AppendSdiffHeaderLine(flags, col_fid, col_sid, dosage_reported, &(cswritep_arr[fidx]));
      }
    } else {
      char* fname_iter = outname_end;
      const uint32_t stream_thread_ct = (max_thread_ct > 1)? (max_thread_ct - 1) : 1;
      if (flags & kfSdiffOneBase) {
        const uint32_t base_sample_idx = id_pairs[0];
        const char* sample_fmtid = &(swc.collapsed_sample_fmtids[base_sample_idx * swc.max_sample_fmtid_blen]);
        const uint32_t base_slen = strlen(sample_fmtid);
        if (unlikely(base_slen >= fname_extrachar_limit)) {
          logputs("\n");
          logerrputs("Error: Sample ID and/or --out argument too long for --sample-diff base= mode.\n");
          goto SdiffMainBatch_ret_INCONSISTENT_INPUT;
        }
        *fname_iter++ = '.';
        char* tab_iter = &(fname_iter[1]);
        fname_iter = memcpya(fname_iter, sample_fmtid, base_slen);
        tab_iter = S_CAST(char*, memchr(tab_iter, '\t', fname_iter - tab_iter));
        if (tab_iter) {
          *tab_iter++ = fname_id_delim;
          tab_iter = S_CAST(char*, memchr(tab_iter, '\t', fname_iter - tab_iter));
          if (tab_iter) {
            *tab_iter = fname_id_delim;
          }
        }
      }
      fname_iter = strcpya_k(fname_iter, ".sdiff");
      if (output_zst) {
        fname_iter = strcpya_k(fname_iter, ".zst");
      }
      *fname_iter = '\0';
      reterr = InitCstreamAlloc(outname, 0, output_zst, stream_thread_ct, kCompressStreamBlock + kMaxMediumLine, &(css_arr[0]), &(cswritep_arr[0]));
      if (unlikely(reterr)) {
        goto SdiffMainBatch_ret_1;
      }
      AppendSdiffHeaderLine(flags, col_fid, col_sid, dosage_reported, &(cswritep_arr[0]));
    }
    if (!(flags & kfSdiffColId)) {
      swc.collapsed_sample_fmtids = nullptr;
    }

    // values unimportant if dosage_hap_tol == kDosageMissing
    const uint32_t dosage_dip_tol = dosage_hap_tol / 2;
    uint32_t dosage_sex_tols[2]; // for chrX
    dosage_sex_tols[0] = dosage_dip_tol;
    dosage_sex_tols[1] = dosage_hap_tol;

    const uint32_t ibs_needed = ((sdip->flags & kfSdiffCountsIbsNeeded) != kfSdiff0);
    PgenVariant pgv;
    if (unlikely(BigstackAllocPgv(sample_ct, allele_idx_offsets != nullptr, dosage_needed? kfPgenGlobalDosagePresent : kfPgenGlobal0, &pgv))) {
      goto SdiffMainBatch_ret_NOMEM;
    }
    Dosage* dosage_buf = nullptr;
    if (dosage_needed) {
      // todo: multidosage
      if (unlikely(bigstack_alloc_dosage(sample_ct, &dosage_buf))) {
        goto SdiffMainBatch_ret_NOMEM;
      }
    }
    AlleleCode* allele_code_buf = nullptr;
    if (allele_idx_offsets != nullptr) {
      if (unlikely(bigstack_alloc_ac(2 * sample_ct, &allele_code_buf))) {
        goto SdiffMainBatch_ret_NOMEM;
      }
    }
    pgv.patch_01_ct = 0;
    pgv.patch_10_ct = 0;
    pgv.dosage_ct = 0;
    // assumes little-endian
    // diploid GT: \t0/0 \t0/1 \t1/1 \t./.
    // diploid DS: \t0 \t1 \t2 \t.
    // haploid GT: \t0 \t0/1 \t1 \t.
    // haploid DS: \t0 \t0.5 \t1 \t.
    // '!' = ascii 33, should never appear in output
    const uint32_t gt_text[8] = {0x302f3009, 0x312f3009, 0x312f3109, 0x2e2f2e09, 0x21213009, 0x312f3009, 0x21213109, 0x21212e09};
    const uint32_t ds_text[8] = {0x21213009, 0x21213109, 0x21213209, 0x21212e09, 0x21213009, 0x352e3009, 0x21213109, 0x21212e09};
    const uint32_t* hc_text = dosage_reported? ds_text : gt_text;
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    const uint32_t include_missing = (flags / kfSdiffIncludeMissing) & 1;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t is_autosomal_diploid = 0;
    uint32_t is_x = 0;
    uint32_t is_y = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t allele_ct = 2;
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        is_autosomal_diploid = !IsSet(cip->haploid_mask, chr_idx);
        is_x = (chr_idx == x_code);
        is_y = (chr_idx == y_code);
        if (swc.chr_buf) {
          char* chr_name_end = chrtoa(cip, chr_idx, swc.chr_buf);
          *chr_name_end = '\t';
          swc.chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - swc.chr_buf);
        }
      }
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      if (allele_ct == 2) {
        if (!dosage_needed) {
          reterr = PgrGet(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec);
        } else {
          reterr = PgrGetD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &pgv.dosage_ct);
        }
      } else {
        if (!dosage_needed) {
          reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
        } else {
          reterr = PgrGetMD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
        }
        if ((!pgv.patch_01_ct) && (!pgv.patch_10_ct)) {
          // todo: also check for multidosage
          allele_ct = 2;
        }
      }
      if (unlikely(reterr)) {
        goto SdiffMainBatch_ret_PGR_FAIL;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }
      if (allele_ct == 2) {
        if (!pgv.dosage_ct) {
          if (AllGenoEqual(pgv.genovec, sample_ct)) {
            // don't need to do anything!
            continue;
          }
          const uint32_t* id_pair_iter = id_pairs;
          for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
            const uint32_t sample_idx1 = *id_pair_iter++;
            const uint32_t sample_idx2 = *id_pair_iter++;
            if (is_y && (!IsSet(pair_sex_male, pair_idx))) {
              continue;
            }
            const uintptr_t hc1 = GetNyparrEntry(pgv.genovec, sample_idx1);
            const uintptr_t hc2 = GetNyparrEntry(pgv.genovec, sample_idx2);
            if (hc1 == hc2) {
              if (hc1 == 3) {
                sdiff_counts[pair_idx].missing_ct += 1;
                sdiff_counts[pair_idx].ibsmiss_ct += is_autosomal_diploid || (is_x && (!IsSet(pair_sex_male, pair_idx)));
              }
              continue;
            }
            SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
            const uint32_t is_diploid_pair = is_autosomal_diploid || (is_x && (!IsSet(pair_sex_male, pair_idx)));
            if ((hc1 == 3) || (hc2 == 3)) {
              cur_counts->halfmiss_ct += 1;
              cur_counts->ibsmiss_ct += is_diploid_pair;
              if (!include_missing) {
                continue;
              }
            } else {
              cur_counts->diff_ct += 1;
              if (is_diploid_pair) {
                cur_counts->ibsx_cts[(hc1 | hc2) & 1] += 1;
              }
            }
            const uintptr_t fidx = is_pairwise? pair_idx : 0;
            char* cswritep = cswritep_arr[fidx];
            const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
            if (unlikely(AppendSdiffPregenoFields(&swc, cur_alleles, sample_idx1, sample_idx2, variant_uidx, allele_ct, &(css_arr[fidx]), &cswritep))) {
              // might not need this assignment, but play it safe for now
              cswritep_arr[fidx] = cswritep;
              goto SdiffMainBatch_ret_WRITE_FAIL;
            }
            if (flags & kfSdiffColGeno) {
              if (is_diploid_pair) {
                cswritep = memcpya(cswritep, &(hc_text[hc1]), 4);
                cswritep -= 2 * dosage_reported;
                cswritep = memcpya(cswritep, &(hc_text[hc2]), 4);
                cswritep -= 2 * dosage_reported;
              } else {
                cswritep = memcpya(cswritep, &(hc_text[hc1 + 4]), 4);
                cswritep -= 2 * (hc1 != 1);
                cswritep = memcpya(cswritep, &(hc_text[hc2 + 4]), 4);
                cswritep -= 2 * (hc2 != 1);
              }
            }
            AppendBinaryEoln(&cswritep);
            cswritep_arr[fidx] = cswritep;
            if (unlikely(Cswrite(&(css_arr[fidx]), &(cswritep_arr[fidx])))) {
              cswritep_arr[fidx] = cswritep;
              goto SdiffMainBatch_ret_WRITE_FAIL;
            }
          }
        } else {
          // dosages present
          const Dosage* dosage_read_iter = pgv.dosage_main;
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            if (!IsSet(pgv.dosage_present, sample_idx)) {
              dosage_buf[sample_idx] = kGenoToDosage[GetNyparrEntry(pgv.genovec, sample_idx)];
            } else {
              dosage_buf[sample_idx] = *dosage_read_iter++;
            }
          }
          const uint32_t* id_pair_iter = id_pairs;
          for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
            const uint32_t sample_idx1 = *id_pair_iter++;
            const uint32_t sample_idx2 = *id_pair_iter++;
            if (is_y && (!IsSet(pair_sex_male, pair_idx))) {
              continue;
            }
            const uint32_t dosage1 = dosage_buf[sample_idx1];
            const uint32_t dosage2 = dosage_buf[sample_idx2];
            if (dosage1 == dosage2) {
              if (dosage1 == kDosageMissing) {
                sdiff_counts[pair_idx].missing_ct += 1;
                sdiff_counts[pair_idx].ibsmiss_ct += is_autosomal_diploid || (is_x && (!IsSet(pair_sex_male, pair_idx)));
              }
              continue;
            }
            SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
            const uint32_t is_diploid_pair = is_autosomal_diploid || (is_x && (!IsSet(pair_sex_male, pair_idx)));
            if ((dosage1 == kDosageMissing) || (dosage2 == kDosageMissing)) {
              cur_counts->halfmiss_ct += 1;
              cur_counts->ibsmiss_ct += is_diploid_pair;
              if (!include_missing) {
                continue;
              }
            } else {
              if (ibs_needed && is_diploid_pair) {
                const uintptr_t hc1 = GetNyparrEntry(pgv.genovec, sample_idx1);
                const uintptr_t hc2 = GetNyparrEntry(pgv.genovec, sample_idx2);
                if (hc1 != hc2) {
                  if ((hc1 == 3) || (hc2 == 3)) {
                    cur_counts->ibsmiss_ct += 1;
                  }
                  cur_counts->ibsx_cts[(hc1 | hc2) & 1] += 1;
                }
              }
              if (abs_i32(dosage1 - dosage2) <= dosage_sex_tols[!is_diploid_pair]) {
                continue;
              }
              cur_counts->diff_ct += 1;
            }
            const uintptr_t fidx = is_pairwise? pair_idx : 0;
            char* cswritep = cswritep_arr[fidx];
            const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
            if (unlikely(AppendSdiffPregenoFields(&swc, cur_alleles, sample_idx1, sample_idx2, variant_uidx, allele_ct, &(css_arr[fidx]), &cswritep))) {
              cswritep_arr[fidx] = cswritep;
              goto SdiffMainBatch_ret_WRITE_FAIL;
            }
            if (flags & kfSdiffColGeno) {
              *cswritep++ = '\t';
              if (dosage1 == kDosageMissing) {
                *cswritep++ = '.';
              } else if (is_diploid_pair) {
                cswritep = PrintSmallDosage(dosage1, cswritep);
              } else if (!(dosage1 % kDosageMax)) {
                *cswritep++ = '0' + (dosage1 / 32768);
              } else {
                cswritep = PrintHaploidNonintDosage(dosage1, cswritep);
              }
              *cswritep++ = '\t';
              if (dosage2 == kDosageMissing) {
                *cswritep++ = '.';
              } else if (is_diploid_pair) {
                cswritep = PrintSmallDosage(dosage2, cswritep);
              } else if (!(dosage2 % kDosageMax)) {
                *cswritep++ = '0' + (dosage2 / 32768);
              } else {
                cswritep = PrintHaploidNonintDosage(dosage2, cswritep);
              }
            }
            AppendBinaryEoln(&cswritep);
            cswritep_arr[fidx] = cswritep;
            if (unlikely(Cswrite(&(css_arr[fidx]), &(cswritep_arr[fidx])))) {
              cswritep_arr[fidx] = cswritep;
              goto SdiffMainBatch_ret_WRITE_FAIL;
            }
          }
        }
      } else {
        // multiallelic
        if (!pgv.dosage_ct) {
          PglMultiallelicSparseToDenseMiss(&pgv, sample_ct, allele_code_buf);
          const uint32_t* id_pair_iter = id_pairs;
          for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
            const uint32_t sample_idx1 = *id_pair_iter++;
            const uint32_t sample_idx2 = *id_pair_iter++;
            if (is_y && (!IsSet(pair_sex_male, pair_idx))) {
              continue;
            }
            const AlleleCode ac11 = allele_code_buf[2 * sample_idx1];
            const AlleleCode ac12 = allele_code_buf[2 * sample_idx1 + 1];
            const AlleleCode ac21 = allele_code_buf[2 * sample_idx2];
            const AlleleCode ac22 = allele_code_buf[2 * sample_idx2 + 1];
            if ((ac11 == ac21) && (ac12 == ac22)) {
              if (ac11 == kMissingAlleleCode) {
                sdiff_counts[pair_idx].missing_ct += 1;
                sdiff_counts[pair_idx].ibsmiss_ct += is_autosomal_diploid || (is_x && (!IsSet(pair_sex_male, pair_idx)));
              }
              continue;
            }
            SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
            const uint32_t is_diploid_pair = is_autosomal_diploid || (is_x && (!IsSet(pair_sex_male, pair_idx)));
            if ((ac11 == kMissingAlleleCode) || (ac21 == kMissingAlleleCode)) {
              cur_counts->halfmiss_ct += 1;
              cur_counts->ibsmiss_ct += is_diploid_pair;
              if (!include_missing) {
                continue;
              }
            } else {
              cur_counts->diff_ct += 1;
              if (ibs_needed && is_diploid_pair) {
                cur_counts->ibsx_cts[(ac11 == ac21) || (ac11 == ac22) || (ac12 == ac21) || (ac12 == ac22)] += 1;
              }
            }
            const uintptr_t fidx = is_pairwise? pair_idx : 0;
            char* cswritep = cswritep_arr[fidx];
            const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
            if (unlikely(AppendSdiffPregenoFields(&swc, cur_alleles, sample_idx1, sample_idx2, variant_uidx, allele_ct, &(css_arr[fidx]), &cswritep))) {
              // might not need this assignment, but play it safe for now
              cswritep_arr[fidx] = cswritep;
              goto SdiffMainBatch_ret_WRITE_FAIL;
            }
            if (flags & kfSdiffColGeno) {
              *cswritep++ = '\t';
              if (is_diploid_pair) {
                if (!dosage_reported) {
                  if (ac11 == kMissingAlleleCode) {
                    cswritep = strcpya_k(cswritep, "./.");
                  } else {
                    cswritep = u32toa_x(ac11, '/', cswritep);
                    cswritep = u32toa(ac12, cswritep);
                  }
                  *cswritep++ = '\t';
                  if (ac21 == kMissingAlleleCode) {
                    cswritep = strcpya_k(cswritep, "./.");
                  } else {
                    cswritep = u32toa_x(ac21, '/', cswritep);
                    cswritep = u32toa(ac22, cswritep);
                  }
                } else {
                  cswritep = PrintMultiallelicHcAsDs(ac11, ac12, allele_ct, cswritep);
                  *cswritep++ = '\t';
                  cswritep = PrintMultiallelicHcAsDs(ac21, ac22, allele_ct, cswritep);
                }
              } else {
                if (!dosage_reported) {
                  if (ac11 == kMissingAlleleCode) {
                    *cswritep++ = '.';
                  } else {
                    cswritep = u32toa(ac11, cswritep);
                    if (ac11 != ac12) {
                      *cswritep++ = '/';
                      cswritep = u32toa(ac12, cswritep);
                    }
                  }
                  *cswritep++ = '\t';
                  if (ac21 == kMissingAlleleCode) {
                    *cswritep++ = '.';
                  } else {
                    cswritep = u32toa(ac21, cswritep);
                    if (ac21 != ac22) {
                      *cswritep++ = '/';
                      cswritep = u32toa(ac22, cswritep);
                    }
                  }
                } else {
                  cswritep = PrintMultiallelicHcAsHaploidDs(ac11, ac12, allele_ct, cswritep);
                  *cswritep++ = '\t';
                  cswritep = PrintMultiallelicHcAsHaploidDs(ac21, ac22, allele_ct, cswritep);
                }
              }
            }
            AppendBinaryEoln(&cswritep);
            cswritep_arr[fidx] = cswritep;
            if (unlikely(Cswrite(&(css_arr[fidx]), &(cswritep_arr[fidx])))) {
              cswritep_arr[fidx] = cswritep;
              goto SdiffMainBatch_ret_WRITE_FAIL;
            }
          }
        } else {
          // multiallelic-dosage; todo
          exit(S_CAST(int32_t, kPglRetInternalError));
        }
      }
    }
    for (uintptr_t fidx = 0; fidx != file_ct; ++fidx) {
      if (unlikely(CswriteCloseNull(&(css_arr[fidx]), cswritep_arr[fidx]))) {
        goto SdiffMainBatch_ret_WRITE_FAIL;
      }
    }
    if (pct > 10) {
      fputs("\b \b", stdout);
    }
    fputs("\b\b", stdout);
  }
  while (0) {
  SdiffMainBatch_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SdiffMainBatch_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  SdiffMainBatch_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  SdiffMainBatch_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 SdiffMainBatch_ret_1:
  if (css_arr) {
    for (uintptr_t fidx = 0; fidx != file_ct; ++fidx) {
      CswriteCloseCond(&(css_arr[fidx]), cswritep_arr[fidx]);
    }
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

CONSTI32(kSdiffBatchMax, (kMaxOpenFiles - 32) & (~(kBitsPerWord - 1)));
static_assert((kSdiffBatchMax % kBitsPerWord) == 0, "kSdiffBatchMax must be a multiple of kBitsPerWord.");
static_assert(kSdiffBatchMax > 0, "kSdiffBatchMax must be positive.");

PglErr Sdiff(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const SdiffInfo* sdip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t variant_ct, uint32_t iid_sid, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  TextStream txs; // file=
  PreinitTextStream(&txs);
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    // Determine xid_mode.
    const SdiffFlags flags = sdip->flags;
    const uint32_t other_id_ct = sdip->other_id_ct;
    uintptr_t line_idx = 0;
    char* line_start = nullptr;
    char* line_iter = nullptr;
    XidMode xid_mode;
    if (other_id_ct) {
      // todo: make this its own function if we use the same logic elsewhere
      // For sanity's sake, require all sample IDs provided on command line
      // to contain the same number of delimiters.
      const char* first_delim_ptr = strchr(sdip->first_id_or_fname, '\t');
      if (!first_delim_ptr) {
        xid_mode = kfXidModeIid;
      } else {
        const char* second_delim_ptr = strchr(&(first_delim_ptr[1]), '\t');
        if (!second_delim_ptr) {
          if (iid_sid) {
            xid_mode = kfXidModeIidSid;
          } else {
            xid_mode = kfXidModeFidIid;
          }
        } else if (unlikely(strchr(&(second_delim_ptr[1]), '\t'))) {
          logerrputs("Error: Too many instances of id-delim= character in --sample-diff sample ID.\n");
          goto Sdiff_ret_INVALID_CMDLINE;
        } else {
          xid_mode = kfXidModeFidIidSid;
        }
      }
    } else {
      reterr = InitTextStream(sdip->first_id_or_fname, kTextStreamBlenFast, MAXV(max_thread_ct - 1, 1), &txs);
      if (unlikely(reterr)) {
        goto Sdiff_ret_TSTREAM_FAIL;
      }
      reterr = LoadXidHeaderPair("sample-diff", iid_sid, &line_idx, &txs, &xid_mode, &line_start, &line_iter);
      if (unlikely(reterr)) {
        if (reterr == kPglRetEof) {
          logerrputs("Error: Empty --sample-diff file.\n");
          reterr = kPglRetInconsistentInput;
        }
        goto Sdiff_ret_TSTREAM_XID_FAIL;
      }
    }
    unsigned char* bigstack_mark2 = g_bigstack_base;
    // May as well have --strict-sid0 or lack of it apply to base=/ids= too.
    const uint32_t skip_sid1 = (xid_mode & kfXidModeFlagSid) && (!siip->sids) && (!(siip->flags & kfSampleIdStrictSid0));
    char* sorted_xidbox;
    uint32_t* xid_map;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(orig_sample_include, siip, orig_sample_ct, 0, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto Sdiff_ret_1;
    }
    char* idbuf;
    if (unlikely(bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto Sdiff_ret_NOMEM;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    // Now we can parse the sample IDs from the file or command line properly.
    uintptr_t* sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    uint32_t* id_pairs;
    uintptr_t id_pair_ct;
    uint32_t sample_ct;
    if (other_id_ct) {
      sample_ct = other_id_ct + 1;
      uint32_t* sample_idxs;
      if (unlikely(bigstack_end_alloc_u32(sample_ct, &sample_idxs))) {
        goto Sdiff_ret_NOMEM;
      }
      const char* first_id = sdip->first_id_or_fname;
      const char* dummy_iter = first_id;
      if (unlikely(SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &dummy_iter, &(sample_idxs[0]), idbuf))) {
        logerrprintfww("Error: --sample-diff sample ID '%s' not found.\n", first_id);
        goto Sdiff_ret_INCONSISTENT_INPUT;
      }
      const char* other_ids_iter = sdip->other_ids_flattened;
      for (uint32_t id_idx = 1; id_idx != sample_ct; ++id_idx) {
        const char* other_id = other_ids_iter;
        if (unlikely(SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &other_ids_iter, &(sample_idxs[id_idx]), idbuf))) {
          logerrprintfww("Error: --sample-diff sample ID '%s' not found.\n", other_id);
          goto Sdiff_ret_INCONSISTENT_INPUT;
        }
        ++other_ids_iter;
      }
      BigstackReset(bigstack_mark2);
      if (flags & kfSdiffOneBase) {
        id_pair_ct = other_id_ct;
      } else {
        id_pair_ct = (S_CAST(uintptr_t, other_id_ct) * sample_ct) / 2;
      }
      if (unlikely(
              bigstack_calloc_w(raw_sample_ctl, &sample_include) ||
              bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
              bigstack_alloc_u32(id_pair_ct * 2, &id_pairs))) {
        goto Sdiff_ret_NOMEM;
      }
      for (uint32_t id_idx = 0; id_idx != sample_ct; ++id_idx) {
        const uint32_t sample_uidx = sample_idxs[id_idx];
        if (IsSet(sample_include, sample_uidx)) {
          logerrputs("Error: Duplicate ID in --sample-diff list.\n");
          goto Sdiff_ret_INVALID_CMDLINE;
        }
        SetBit(sample_uidx, sample_include);
      }
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      UidxsToIdxs(sample_include, sample_include_cumulative_popcounts, sample_ct, sample_idxs);
      if (flags & kfSdiffOneBase) {
        const uint32_t sample_idx_base = sample_idxs[0];
        for (uint32_t id_idx = 1; id_idx != sample_ct; ++id_idx) {
          id_pairs[2 * id_idx - 2] = sample_idx_base;
          id_pairs[2 * id_idx - 1] = sample_idxs[id_idx];
        }
      } else {
        uint32_t* id_pairs_iter = id_pairs;
        for (uint32_t id_idx2 = 1; id_idx2 != sample_ct; ++id_idx2) {
          const uint32_t sample_idx2 = sample_idxs[id_idx2];
          for (uint32_t id_idx1 = 0; id_idx1 != id_idx2; ++id_idx1) {
            *id_pairs_iter++ = sample_idxs[id_idx1];
            *id_pairs_iter++ = sample_idx2;
          }
        }
      }
      BigstackEndReset(bigstack_end_mark);
    } else {
      uint32_t* sample_pair_uidxs_end = R_CAST(uint32_t*, bigstack_end_mark);
      uint32_t* sample_pair_uidxs_iter = sample_pair_uidxs_end;
      uint32_t* sample_pair_uidxs_oom = R_CAST(uint32_t*, g_bigstack_base);
      if (*line_start == '#') {
        line_iter = AdvPastDelim(line_iter, '\n');
        ++line_idx;
      } else {
        line_iter = line_start;
      }
      for (; TextGetUnsafe2(&txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
        const char* linebuf_iter = line_iter;
        uint32_t sample_uidx1;
        if (unlikely(SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &linebuf_iter, &sample_uidx1, idbuf))) {
          TabsToSpaces(idbuf);
          logerrprintfww("Error: --sample-diff sample ID '%s' (on line %" PRIuPTR " of file) not found.\n", idbuf, line_idx);
          goto Sdiff_ret_INCONSISTENT_INPUT;
        }
        linebuf_iter = FirstNonTspace(linebuf_iter);
        if (skip_sid1) {
          linebuf_iter = FirstNonTspace(FirstSpaceOrEoln(linebuf_iter));
        }
        uint32_t sample_uidx2;
        if (unlikely(SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &linebuf_iter, &sample_uidx2, idbuf))) {
          TabsToSpaces(idbuf);
          logerrprintfww("Error: --sample-diff sample ID '%s' (on line %" PRIuPTR " of file) not found.\n", idbuf, line_idx);
          goto Sdiff_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(sample_uidx1 == sample_uidx2)) {
          TabsToSpaces(idbuf);
          logerrprintfww("Error: Duplicate sample ID '%s' on line %" PRIuPTR " of --sample-diff file.\n", idbuf, line_idx);
          goto Sdiff_ret_MALFORMED_INPUT;
        }
        if (unlikely(sample_pair_uidxs_iter == sample_pair_uidxs_oom)) {
          goto Sdiff_ret_NOMEM;
        }
        *(--sample_pair_uidxs_iter) = sample_uidx1;
        *(--sample_pair_uidxs_iter) = sample_uidx2;
        line_iter = K_CAST(char*, linebuf_iter);
      }
      if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
        goto Sdiff_ret_TSTREAM_FAIL;
      }
      const uintptr_t id_ct = sample_pair_uidxs_end - sample_pair_uidxs_iter;
      if (unlikely(!id_ct)) {
        logerrputs("Error: Empty --sample-diff file.\n");
        goto Sdiff_ret_MALFORMED_INPUT;
      }
      BigstackEndSet(sample_pair_uidxs_iter);
      BigstackReset(bigstack_mark2);
      if (unlikely(
              bigstack_calloc_w(raw_sample_ctl, &sample_include) ||
              bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
              bigstack_alloc_u32(id_ct, &id_pairs))) {
        goto Sdiff_ret_NOMEM;
      }
      sample_pair_uidxs_iter = sample_pair_uidxs_end;
      for (uintptr_t ulii = 0; ulii != id_ct; ++ulii) {
        const uint32_t sample_uidx = *(--sample_pair_uidxs_iter);
        SetBit(sample_uidx, sample_include);
        id_pairs[ulii] = sample_uidx;
      }
      BigstackEndReset(bigstack_end_mark);
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      sample_ct = sample_include_cumulative_popcounts[raw_sample_ctl - 1] + PopcountWord(sample_include[raw_sample_ctl - 1]);
      // Don't strictly enforce pair-uniqueness for now, but do perform a basic
      // sanity check.
      if (id_ct > S_CAST(uintptr_t, sample_ct) * (sample_ct - 1)) {
        logerrputs("Error: Duplicate ID pairs in --sample-diff file.\n");
        goto Sdiff_ret_MALFORMED_INPUT;
      }
      UidxsToIdxs(sample_include, sample_include_cumulative_popcounts, id_ct, id_pairs);
      id_pair_ct = id_ct / 2;
    }

    SdiffCounts* sdiff_counts;
    if (unlikely(BIGSTACK_ALLOC_X(SdiffCounts, id_pair_ct, &sdiff_counts))) {
      goto Sdiff_ret_NOMEM;
    }
    memset(sdiff_counts, 0, id_pair_ct * sizeof(SdiffCounts));
    uintptr_t* pair_sex_male = nullptr;
    uint32_t x_ct = 0;
    uint32_t y_ct = 0;
    {
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      if ((!IsI32Neg(x_code)) && IsSet(cip->chr_mask, x_code)) {
        x_ct = CountChrVariantsUnsafe(variant_include, cip, x_code);
      }
      const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
      if ((!IsI32Neg(y_code)) && IsSet(cip->chr_mask, y_code)) {
        y_ct = CountChrVariantsUnsafe(variant_include, cip, y_code);
      }
    }
    if (x_ct || y_ct) {
      const uintptr_t id_pair_ctl = BitCtToWordCt(id_pair_ct);
      const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
      uintptr_t* sex_nm_collapsed;
      uintptr_t* sex_male_collapsed;
      if (unlikely(
              bigstack_alloc_w(id_pair_ctl, &pair_sex_male) ||
              bigstack_alloc_w(sample_ctl, &sex_nm_collapsed) ||
              bigstack_alloc_w(sample_ctl, &sex_male_collapsed))) {
        goto Sdiff_ret_NOMEM;
      }
      CopyBitarrSubset(sex_nm, sample_include, sample_ct, sex_nm_collapsed);
      CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_male_collapsed);
      const uint32_t* id_pairs_iter = id_pairs;
      const uintptr_t word_ct_m1 = id_pair_ctl - 1;
      uint32_t loop_len = kBitsPerWord;
      for (uintptr_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          loop_len = ModNz(id_pair_ct, kBitsPerWord);
        }
        uintptr_t cur_word = 0;
        for (uint32_t uii = 0; uii != loop_len; ++uii) {
          const uint32_t sample_idx1 = *id_pairs_iter++;
          const uint32_t sample_idx2 = *id_pairs_iter++;
          const uint32_t sex_nm1 = IsSet(sex_nm_collapsed, sample_idx1);
          const uint32_t sex_nm_xor = sex_nm1 ^ IsSet(sex_nm_collapsed, sample_idx2);
          const uintptr_t sex_male1 = IsSet(sex_male_collapsed, sample_idx1);
          const uintptr_t sex_male_xor = sex_male1 ^ IsSet(sex_male_collapsed, sample_idx2);
          if (sex_nm_xor) {
            // Exactly one sex is missing.
            cur_word |= sex_male_xor << uii;
          } else if (unlikely(!sex_nm1)) {
            logerrputs("Error: --sample-diff requires at least one sample in each pair to have known\nsex when chrX or chrY is present.\n");
            goto Sdiff_ret_INCONSISTENT_INPUT;
          } else if (unlikely(sex_male_xor)) {
            logerrputs("Error: --sample-diff cannot perform male-female comparisons when chrX or chrY\nis present.  (Consider \"--not-chr X,Y\".)\n");
            goto Sdiff_ret_INCONSISTENT_INPUT;
          } else {
            cur_word |= sex_male1 << uii;
          }
        }
        pair_sex_male[widx] = cur_word;
      }
      BigstackReset(sex_nm_collapsed);
    }
    // Main loops.
    if (flags & kfSdiffCountsOnly) {
      reterr = SdiffCountsOnly(sample_include, sample_include_cumulative_popcounts, id_pairs, pair_sex_male, variant_include, cip, allele_idx_offsets, sdip, sample_ct, variant_ct, id_pair_ct, simple_pgrp, sdiff_counts);
      if (unlikely(reterr)) {
        goto Sdiff_ret_1;
      }
    } else {
      const uint32_t is_pairwise = (flags / kfSdiffPairwise) & 1;
      uintptr_t batch_ct = 1;
      if (is_pairwise && (id_pair_ct > kSdiffBatchMax)) {
        batch_ct = 1 + ((id_pair_ct - 1) / kSdiffBatchMax);
      }
      for (uintptr_t batch_idx = 0; batch_idx != batch_ct; ++batch_idx) {
        const uintptr_t id_pair_idx_start = batch_idx * kSdiffBatchMax;
        uintptr_t id_pair_idx_end = id_pair_ct;
        if (batch_ct == 1) {
          printf("--sample-diff%s: 0%%", is_pairwise? " pairwise" : "");
        } else {
          if (batch_idx != batch_ct - 1) {
            id_pair_idx_end = id_pair_idx_start + kSdiffBatchMax;
          }
          printf("\r--sample-diff pairwise batch %" PRIuPTR "/%" PRIuPTR ": 0%%", batch_idx + 1, batch_ct);
        }
        fflush(stdout);
        reterr = SdiffMainBatch(sample_include, sample_include_cumulative_popcounts, siip, &(id_pairs[2 * id_pair_idx_start]), &(pair_sex_male[id_pair_idx_start / kBitsPerWord]), variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, sdip, sample_ct, variant_ct, id_pair_idx_end - id_pair_idx_start, max_thread_ct, simple_pgrp, outname, outname_end, &(sdiff_counts[id_pair_idx_start]));
        if (unlikely(reterr)) {
          goto Sdiff_ret_1;
        }
      }
      fputs("done.\n", stdout);
      if (is_pairwise) {
        *outname_end = '\0';
        logprintfww("--sample-diff pairwise: Discordances written to %s.[ID1].[ID2].sdiff%s (%" PRIuPTR " file%s).\n", outname, (flags & kfSdiffZs)? ".zst" : "", id_pair_ct, (id_pair_ct == 1)? "" : "s");
      } else {
        logprintfww("--sample-diff: Discordances written to %s .\n", outname);
      }
    }

    // Final count-summary.
    const uint32_t col_fid = FidColIsRequired(siip, flags / kfSdiffCountsColMaybefid);
    const uint32_t col_sid = SidColIsRequired(siip->sids, flags / kfSdiffCountsColMaybesid);
    char* collapsed_sample_fmtids;
    uintptr_t max_sample_fmtid_blen;
    if (unlikely(CollapsedSampleFmtidInitAlloc(sample_include, siip, sample_ct, col_fid, col_sid, &collapsed_sample_fmtids, &max_sample_fmtid_blen))) {
      goto Sdiff_ret_NOMEM;
    }
    uint32_t autosomal_diploid_ct;
    uint32_t nonsex_haploid_ct;
    if (!IsSet(cip->chr_mask, 0)) {
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
      const uint32_t chr_idx_end = cip->max_code + 1 + cip->name_ct;
      nonsex_haploid_ct = 0;
      for (uint32_t chr_idx = 1; ; ++chr_idx) {
        // usually just MT
        chr_idx = AdvBoundedTo1Bit(cip->haploid_mask, chr_idx, chr_idx_end);
        if (chr_idx == chr_idx_end) {
          break;
        }
        if ((!IsSet(cip->chr_mask, chr_idx)) || (chr_idx == x_code) || (chr_idx == y_code)) {
          continue;
        }
        nonsex_haploid_ct += CountChrVariantsUnsafe(variant_include, cip, chr_idx);
      }
      autosomal_diploid_ct = variant_ct - x_ct - y_ct - nonsex_haploid_ct;
    } else {
      autosomal_diploid_ct = 0;
      nonsex_haploid_ct = variant_ct - x_ct - y_ct;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".sdiff.summary");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto Sdiff_ret_OPEN_FAIL;
    }
    char* textbuf = g_textbuf;
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);
    char* write_iter = textbuf;
    *write_iter++ = '#';
    if (col_fid) {
      write_iter = strcpya_k(write_iter, "FID1\t");
    }
    write_iter = strcpya_k(write_iter, "IID1");
    if (col_sid) {
      write_iter = strcpya_k(write_iter, "\tSID1");
    }
    if (col_fid) {
      write_iter = strcpya_k(write_iter, "\tFID2");
    }
    write_iter = strcpya_k(write_iter, "\tIID2");
    if (col_sid) {
      write_iter = strcpya_k(write_iter, "\tSID2");
    }
    if (flags & kfSdiffCountsColNobs) {
      write_iter = strcpya_k(write_iter, "\tOBS_CT");
    }
    if (flags & kfSdiffCountsColNobsIbs) {
      write_iter = strcpya_k(write_iter, "\tIBS_OBS_CT");
    }
    if (flags & kfSdiffCountsColIbs0) {
      write_iter = strcpya_k(write_iter, "\tIBS0_CT");
    }
    if (flags & kfSdiffCountsColIbs1) {
      write_iter = strcpya_k(write_iter, "\tIBS1_CT");
    }
    if (flags & kfSdiffCountsColIbs2) {
      write_iter = strcpya_k(write_iter, "\tIBS2_CT");
    }
    if (flags & kfSdiffCountsColHalfmiss) {
      write_iter = strcpya_k(write_iter, "\tHALFMISS_CT");
    }
    if (flags & kfSdiffCountsColDiff) {
      write_iter = strcpya_k(write_iter, "\tDIFF_CT");
    }
    AppendBinaryEoln(&write_iter);
    const uint32_t obs_ct_base = autosomal_diploid_ct + x_ct + nonsex_haploid_ct;
    const uint32_t* id_pair_iter = id_pairs;
    for (uintptr_t pair_idx = 0; pair_idx != id_pair_ct; ++pair_idx) {
      const uint32_t sample_idx1 = *id_pair_iter++;
      const uint32_t sample_idx2 = *id_pair_iter++;
      write_iter = strcpyax(write_iter, &(collapsed_sample_fmtids[max_sample_fmtid_blen * sample_idx1]), '\t');
      write_iter = strcpya(write_iter, &(collapsed_sample_fmtids[max_sample_fmtid_blen * sample_idx2]));
      const SdiffCounts* cur_counts = &(sdiff_counts[pair_idx]);
      if (flags & kfSdiffCountsColNobs) {
        uint32_t obs_ct = obs_ct_base;
        if (y_ct && IsSet(pair_sex_male, pair_idx)) {
          obs_ct += y_ct;
        }
        if (!(flags & kfSdiffIncludeMissing)) {
          obs_ct -= cur_counts->missing_ct + cur_counts->halfmiss_ct;
        }
        *write_iter++ = '\t';
        write_iter = u32toa(obs_ct, write_iter);
      }
      if (flags & kfSdiffCountsIbsNeeded) {
        uint32_t ibs_obs_ct = autosomal_diploid_ct - cur_counts->ibsmiss_ct;
        if (x_ct && (!IsSet(pair_sex_male, pair_idx))) {
          ibs_obs_ct += x_ct;
        }
        if (flags & kfSdiffCountsColNobsIbs) {
          *write_iter++ = '\t';
          write_iter = u32toa(ibs_obs_ct, write_iter);
        }
        if (flags & kfSdiffCountsColIbs0) {
          *write_iter++ = '\t';
          write_iter = u32toa(cur_counts->ibsx_cts[0], write_iter);
        }
        if (flags & kfSdiffCountsColIbs1) {
          *write_iter++ = '\t';
          write_iter = u32toa(cur_counts->ibsx_cts[1], write_iter);
        }
        if (flags & kfSdiffCountsColIbs2) {
          *write_iter++ = '\t';
          write_iter = u32toa(ibs_obs_ct - cur_counts->ibsx_cts[0] - cur_counts->ibsx_cts[1], write_iter);
        }
      }
      if (flags & kfSdiffCountsColHalfmiss) {
        *write_iter++ = '\t';
        write_iter = u32toa(cur_counts->halfmiss_ct, write_iter);
      }
      if (flags & kfSdiffCountsColDiff) {
        uint32_t diff_ct = cur_counts->diff_ct;
        if (flags & kfSdiffIncludeMissing) {
          diff_ct += cur_counts->halfmiss_ct;
        }
        *write_iter++ = '\t';
        write_iter = u32toa(diff_ct, write_iter);
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto Sdiff_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto Sdiff_ret_WRITE_FAIL;
    }
    logprintfww("--sample-diff: Discordance count summary written to %s .\n", outname);
  }
  while (0) {
  Sdiff_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Sdiff_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  Sdiff_ret_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&txs)) {
      break;
    }
  Sdiff_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--sample-diff file", &txs);
    break;
  Sdiff_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  Sdiff_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  Sdiff_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  Sdiff_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 Sdiff_ret_1:
  fclose_cond(outfile);
  CleanupTextStream2("--sample-diff file", &txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

PglErr WriteSnplist(const uintptr_t* variant_include, const char* const* variant_ids, uint32_t variant_ct, uint32_t output_zst, uint32_t allow_dups, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    if (!allow_dups) {
      uint32_t dup_found;
      reterr = CheckIdUniqueness(g_bigstack_base, g_bigstack_end, variant_include, variant_ids, variant_ct, max_thread_ct, &dup_found);
      if (unlikely(reterr)) {
        goto WriteSnplist_ret_1;
      }
      if (dup_found) {
        // I estimate that this corresponds to a pipeline bug (or, at minimum,
        // inefficiency) >90% of the time.
        // --extract[-intersect] could be modified in a similar manner, but my
        // working hypothesis is that interception of suspect --write-snplist
        // operations is enough to stop most duplicate-related --extract
        // misuses.
        logerrputs("Error: --write-snplist normally shouldn't be used with duplicate variant IDs.\n(--set-all-var-ids helps with ID deduplication, and --rm-dup addresses actual\nduplicate data.)  However, if you've already accounted for them, you can\nsuppress this error with the 'allow-dups' modifier.\n");
        goto WriteSnplist_ret_INCONSISTENT_INPUT;
      }
    }
    OutnameZstSet(".snplist", output_zst, outname_end);
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, kCompressStreamBlock + kMaxIdSlen + 2, &css, &cswritep);
    if (unlikely(reterr)) {
      goto WriteSnplist_ret_1;
    }
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto WriteSnplist_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto WriteSnplist_ret_WRITE_FAIL;
    }
    logprintfww("--write-snplist%s%s: Variant IDs written to %s .\n", output_zst? " zs" : "", allow_dups? " allow-dups" : "", outname);
  }
  while (0) {
  WriteSnplist_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  WriteSnplist_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 WriteSnplist_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

// similar to write_psam().
PglErr WriteCovar(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, const uint32_t* new_sample_idx_to_old, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, WriteCovarFlags write_covar_flags, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    snprintf(outname_end, kMaxOutfnameExtBlen, ".cov");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto WriteCovar_ret_OPEN_FAIL;
    }
    const char* output_missing_pheno = g_output_missing_pheno;
    const uint32_t omp_slen = strlen(output_missing_pheno);

    char* textbuf = g_textbuf;
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);

    const uint32_t write_fid = FidColIsRequired(&piip->sii, write_covar_flags / kfWriteCovarColMaybefid);
    const char* sample_ids = piip->sii.sample_ids;
    const char* sids = piip->sii.sids;
    const char* paternal_ids = piip->parental_id_info.paternal_ids;
    const char* maternal_ids = piip->parental_id_info.maternal_ids;
    const uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    const uintptr_t max_sid_blen = piip->sii.max_sid_blen;
    const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
    const uint32_t write_sid = SidColIsRequired(sids, write_covar_flags / kfWriteCovarColMaybesid);
    const uint32_t write_parents = ParentalColsAreRequired(piip, write_covar_flags / kfWriteCovarColMaybeparents);
    const uint32_t write_sex = (write_covar_flags / kfWriteCovarColSex) & 1;
    const uint32_t write_empty_pheno = (write_covar_flags & kfWriteCovarColPheno1) && (!pheno_ct);
    const uint32_t write_phenos = (write_covar_flags & (kfWriteCovarColPheno1 | kfWriteCovarColPhenos)) && pheno_ct;
    if (write_phenos && (!(write_covar_flags & kfWriteCovarColPhenos))) {
      pheno_ct = 1;
    }
    char* write_iter = textbuf;
    *write_iter++ = '#';
    if (write_fid) {
      write_iter = strcpya_k(write_iter, "FID\t");
    }
    write_iter = strcpya_k(write_iter, "IID");
    if (write_sid) {
      write_iter = strcpya_k(write_iter, "\tSID");
    }
    if (write_parents) {
      write_iter = strcpya_k(write_iter, "\tPAT\tMAT");
    }
    if (write_sex) {
      write_iter = strcpya_k(write_iter, "\tSEX");
    }
    if (write_phenos || write_empty_pheno || write_sex) {
      // verify that no names are duplicated
      uint32_t* covar_name_htable;
      uint32_t covar_name_htable_size;
      if (unlikely(HtableGoodSizeAlloc(covar_ct + write_sex, bigstack_left(), &covar_name_htable, &covar_name_htable_size))) {
        goto WriteCovar_ret_NOMEM;
      }
      // shouldn't be possible for this to fail
      PopulateStrboxHtable(covar_names, covar_ct, max_covar_name_blen, covar_name_htable_size, covar_name_htable);
      uint32_t max_xcovar_name_blen = max_covar_name_blen;
      if (write_sex) {
        // add "SEX"
        for (uint32_t hashval = Hashceil("SEX", 3, covar_name_htable_size); ; ) {
          const uint32_t cur_htable_entry = covar_name_htable[hashval];
          if (cur_htable_entry == UINT32_MAX) {
            covar_name_htable[hashval] = covar_ct;
            break;
          }
          if (unlikely(strequal_k_unsafe(&(covar_names[cur_htable_entry * max_covar_name_blen]), "SEX"))) {
            logerrputs("Error: .cov file cannot have both a regular SEX column and a covariate named\n'SEX'.  Exclude or rename one of these columns.\n");
            goto WriteCovar_ret_INCONSISTENT_INPUT;
          }
          if (++hashval == covar_name_htable_size) {
            hashval = 0;
          }
        }
        if (max_xcovar_name_blen < 4) {
          max_xcovar_name_blen = 4;
        }
      }
      if (write_phenos) {
        const char* pheno_name_iter = pheno_names;
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          const uint32_t cur_pheno_name_slen = strlen(pheno_name_iter);
          if (cur_pheno_name_slen < max_xcovar_name_blen) {
            const uint32_t cur_pheno_name_blen = cur_pheno_name_slen + 1;
            for (uint32_t hashval = Hashceil(pheno_name_iter, cur_pheno_name_slen, covar_name_htable_size); ; ) {
              const uint32_t cur_htable_idval = covar_name_htable[hashval];
              if (cur_htable_idval >= covar_ct) {
                if (cur_htable_idval == UINT32_MAX) {
                  break;
                }
                if (unlikely(strequal_k_unsafe(pheno_name_iter, "SEX"))) {
                  logerrputs(write_sex? "Error: .cov file cannot have both a regular SEX column and a phenotype named\n'SEX'.  Exclude or rename one of these columns.\n" : "Error: .cov file cannot have a phenotype and a covariate with the same name.\n");
                  goto WriteCovar_ret_INCONSISTENT_INPUT;
                }
              } else {
                if (unlikely(memequal(pheno_name_iter, &(covar_names[cur_htable_idval * max_covar_name_blen]), cur_pheno_name_blen))) {
                  logerrputs("Error: .cov file cannot have a phenotype and a covariate with the same name.\n");
                  goto WriteCovar_ret_INCONSISTENT_INPUT;
                }
              }
              if (++hashval == covar_name_htable_size) {
                hashval = 0;
              }
            }
          }
          write_iter = memcpya(write_iter, pheno_name_iter, cur_pheno_name_slen);
          pheno_name_iter = &(pheno_name_iter[max_pheno_name_blen]);
          if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
            goto WriteCovar_ret_WRITE_FAIL;
          }
        }
      } else if (write_empty_pheno) {
        if (max_covar_name_blen > 6) {
          for (uint32_t hashval = Hashceil("PHENO1", 6, covar_name_htable_size); ; ) {
            const uint32_t cur_htable_idval = covar_name_htable[hashval];
            if (cur_htable_idval >= covar_ct) {
              if (cur_htable_idval == UINT32_MAX) {
                break;
              }
            } else {
              if (unlikely(strequal_k_unsafe(&(covar_names[cur_htable_idval * max_covar_name_blen]), "PHENO1"))) {
                logerrputs("Error: .cov file cannot have a phenotype and a covariate with the same name.\n");
                goto WriteCovar_ret_INCONSISTENT_INPUT;
              }
            }
            if (++hashval == covar_name_htable_size) {
              hashval = 0;
            }
          }
        }
        write_iter = strcpya_k(write_iter, "\tPHENO1");
      }
    }
    for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
      *write_iter++ = '\t';
      const char* cur_covar_name = &(covar_names[covar_idx * max_covar_name_blen]);
      const uint32_t cur_covar_name_slen = strlen(cur_covar_name);
      write_iter = memcpya(write_iter, cur_covar_name, cur_covar_name_slen);
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto WriteCovar_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);

    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    uint32_t sample_uidx2 = 0;
    // not really necessary to make sample_uidx increment dependent on
    // new_sample_idx_to_old == nullptr
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      uintptr_t sample_uidx;
      if (!new_sample_idx_to_old) {
        sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      } else {
        do {
          sample_uidx = new_sample_idx_to_old[sample_uidx2++];
        } while (!IsSet(sample_include, sample_uidx));
      }
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
      if (write_parents) {
        *write_iter++ = '\t';
        write_iter = strcpyax(write_iter, &(paternal_ids[max_paternal_id_blen * sample_uidx]), '\t');
        write_iter = strcpya(write_iter, &(maternal_ids[max_maternal_id_blen * sample_uidx]));
      }
      if (write_sex) {
        *write_iter++ = '\t';
        if (IsSet(sex_nm, sample_uidx)) {
          *write_iter++ = '2' - IsSet(sex_male, sample_uidx);
        } else {
          // this is better than '0' since it allows the raw column to be used
          // as --covar input
          // (can't do this for .fam export, though: not worth the
          // compatibility issues)
          write_iter = strcpya_k(write_iter, "NA");
        }
      }
      if (write_phenos) {
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          write_iter = AppendPhenoStr(&(pheno_cols[pheno_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
          if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
            goto WriteCovar_ret_WRITE_FAIL;
          }
        }
      } else {
        if (write_empty_pheno) {
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
        }
        if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
          goto WriteCovar_ret_WRITE_FAIL;
        }
      }
      for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
        *write_iter++ = '\t';
        write_iter = AppendPhenoStr(&(covar_cols[covar_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
        if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
          goto WriteCovar_ret_WRITE_FAIL;
        }
      }
      AppendBinaryEoln(&write_iter);
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto WriteCovar_ret_WRITE_FAIL;
    }
    logprintfww("Covariates written to %s .\n", outname);
  }
  while (0) {
  WriteCovar_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  WriteCovar_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  WriteCovar_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  WriteCovar_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

// This is almost a subset of SampleCounts().
typedef struct HetCtxStruct {
  const uintptr_t* autosomal_variant_include;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  const uintptr_t* founder_info_collapsed;  // only non-null if small-sample
  const uintptr_t* founder_info_interleaved_vec;
  const double* allele_freqs;  // nullptr if small-sample
  uint32_t sample_ct;
  uint32_t founder_ct;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** raregenos;
  uint32_t** difflist_sample_id_bufs;
  uint32_t** allele_nobs_bufs;
  VecW** scrambled_ohet_bufs;
  uint32_t* read_variant_uidx_starts;
  uint32_t cur_block_size;

  // only kPglRetMalformedInput possible, no atomic ops needed
  PglErr reterr;

  double* thread_ehet_base;
  uint32_t* thread_nobs_base;
  uint32_t** thread_ohets;
  double** thread_ehet_incrs;
  int32_t** thread_nobs_incrs;
} HetCtx;

THREAD_FUNC_DECL HetThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  HetCtx* ctx = S_CAST(HetCtx*, arg->sharedp->context);

  const uintptr_t* autosomal_variant_include = ctx->autosomal_variant_include;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const double* allele_freqs = ctx->allele_freqs;
  const uintptr_t* sample_include = ctx->sample_include;
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  const uintptr_t* founder_info_collapsed = ctx->founder_info_collapsed;
  const uintptr_t* founder_info_interleaved_vec = ctx->founder_info_interleaved_vec;
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t founder_ct = ctx->founder_ct;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  uintptr_t* genovec = ctx->genovecs[tidx];
  PgenVariant pgv;
  pgv.genovec = genovec;
  SetPgvThreadMhcNull(sample_ct, tidx, ctx->thread_read_mhc, &pgv);
  uintptr_t* raregeno = ctx->raregenos[tidx];
  uint32_t* difflist_sample_ids = ctx->difflist_sample_id_bufs[tidx];
  uint32_t* allele_nobs = nullptr;
  if (ctx->thread_read_mhc) {
    allele_nobs = ctx->allele_nobs_bufs[tidx];
  }

  // todo: tune this threshold
  const uint32_t max_simple_difflist_len = sample_ct / 32;

  double ehet_base = 0.0;
  uint32_t nobs_base = 0;
  uint32_t* ohets = ctx->thread_ohets[tidx];
  double* ehet_incrs = ctx->thread_ehet_incrs[tidx];
  int32_t* nobs_incrs = ctx->thread_nobs_incrs[tidx];
  ZeroU32Arr(RoundUpPow2(sample_ct, kInt32PerVec), ohets);
  ZeroDArr(sample_ct, ehet_incrs);
  ZeroI32Arr(RoundUpPow2(sample_ct, kInt32PerVec), nobs_incrs);
  const uint32_t acc2_vec_ct = NypCtToVecCt(sample_ct);
  const uintptr_t dense_counts_vstride = acc2_vec_ct * 23;
  VecW* scrambled_ohets = ctx->scrambled_ohet_bufs[tidx];
  uint32_t dense_ct_rem3 = 3;
  uint32_t dense_ct_rem15d3 = 5;
  uint32_t dense_ct_rem255d15 = 17;
  ZeroVecArr(dense_counts_vstride, scrambled_ohets);

  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);

  double ehet = 0.0;
  uint32_t cur_allele_ct = 2;
  do {
    const uint32_t cur_block_size = ctx->cur_block_size;
    const uint32_t cur_idx_ct = (((tidx + 1) * cur_block_size) / calc_thread_ct) - ((tidx * cur_block_size) / calc_thread_ct);
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(autosomal_variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
    for (uint32_t cur_idx = 0; cur_idx != cur_idx_ct; ++cur_idx) {
      const uint32_t variant_uidx = BitIter1(autosomal_variant_include, &variant_uidx_base, &cur_bits);
      uintptr_t allele_idx_offset_base;
      if (!allele_idx_offsets) {
        allele_idx_offset_base = variant_uidx * 2;
      } else {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      if (cur_allele_ct == 2) {
        if (allele_freqs) {
          const double ref_freq = allele_freqs[allele_idx_offset_base - variant_uidx];
          ehet = 2 * ref_freq * (1 - ref_freq);
          if (ehet < kSmallishEpsilon) {
            continue;
          }
        }
        uint32_t difflist_common_geno;
        uint32_t difflist_len;
        PglErr reterr = PgrGetDifflistOrGenovec(sample_include, pssi, sample_ct, max_simple_difflist_len, variant_uidx, pgrp, genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
        if (unlikely(reterr)) {
          ctx->reterr = reterr;
          break;
        }
        if (difflist_common_geno != UINT32_MAX) {
          ZeroTrailingNyps(difflist_len, raregeno);
          const uint32_t word_ct = NypCtToWordCt(difflist_len);
          if (!allele_freqs) {
            STD_ARRAY_DECL(uint32_t, 4, genocounts);
            genocounts[0] = 0;
            genocounts[1] = 0;
            genocounts[2] = 0;
            genocounts[3] = 0;
            if (word_ct) {
              const uint32_t word_ct_m1 = word_ct - 1;
              uint32_t loop_len = kBitsPerWordD2;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= word_ct_m1) {
                  if (widx > word_ct_m1) {
                    break;
                  }
                  loop_len = ModNz(difflist_len, kBitsPerWordD2);
                }
                const uint32_t difflist_idx_base = widx * kBitsPerWordD2;
                const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[difflist_idx_base]);
                uintptr_t raregeno_word = raregeno[widx];
                for (uint32_t uii = 0; uii != loop_len; ++uii) {
                  const uint32_t sample_idx = cur_difflist_sample_ids[uii];
                  genocounts[raregeno_word & 3] += IsSet(founder_info_collapsed, sample_idx);
                  raregeno_word >>= 2;
                }
              }
            }
            genocounts[difflist_common_geno] += founder_ct - difflist_len;
            const uint32_t numer1 = 2 * genocounts[0] + genocounts[1];
            const uint32_t numer2 = genocounts[1] + 2 * genocounts[2];
            if ((!numer1) || (!numer2)) {
              continue;
            }
            const uint32_t denom = numer1 + numer2;
            ehet = 2 * S_CAST(double, numer1) * S_CAST(double, numer2) / ((S_CAST(double, denom)) * S_CAST(double, denom - 1));
          }
          if (difflist_common_geno != 3) {
            ehet_base += ehet;
            ++nobs_base;
            for (uint32_t widx = 0; widx != word_ct; ++widx) {
              const uintptr_t raregeno_word = raregeno[widx];
              const uint32_t difflist_idx_base = widx * kBitsPerWordD2;
              const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[difflist_idx_base]);
              uintptr_t het_word = (raregeno_word & kMask5555) & (~(raregeno_word >> 1));
              while (het_word) {
                const uint32_t difflist_idx_lowbits = ctzw(het_word) / 2;
                const uint32_t sample_idx = cur_difflist_sample_ids[difflist_idx_lowbits];
                ohets[sample_idx] += 1;
                het_word &= het_word - 1;
              }
              uintptr_t missing_word = (raregeno_word & kMask5555) & (raregeno_word >> 1);
              while (missing_word) {
                const uint32_t difflist_idx_lowbits = ctzw(missing_word) / 2;
                const uint32_t sample_idx = cur_difflist_sample_ids[difflist_idx_lowbits];
                ehet_incrs[sample_idx] -= ehet;
                nobs_incrs[sample_idx] -= 1;
                missing_word &= missing_word - 1;
              }
            }
          } else {
            if (!word_ct) {
              continue;
            }
            const uint32_t word_ct_m1 = word_ct - 1;
            uint32_t loop_len = kBitsPerWordD2;
            for (uint32_t widx = 0; ; ++widx) {
              if (widx >= word_ct_m1) {
                if (widx > word_ct_m1) {
                  break;
                }
                loop_len = ModNz(difflist_len, kBitsPerWordD2);
              }
              uintptr_t raregeno_word = raregeno[widx];
              const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
              for (uint32_t uii = 0; uii != loop_len; ++uii) {
                const uint32_t sample_idx = cur_difflist_sample_ids[uii];
                ohets[sample_idx] += raregeno_word & 1;
                ehet_incrs[sample_idx] += ehet;
                nobs_incrs[sample_idx] += 1;
                raregeno_word >>= 2;
              }
            }
          }
          continue;
        }
        ZeroTrailingNyps(sample_ct, genovec);
        if (!allele_freqs) {
          STD_ARRAY_DECL(uint32_t, 4, genocounts);
          GenoarrCountSubsetFreqs(genovec, founder_info_interleaved_vec, sample_ct, founder_ct, genocounts);
          const uint32_t numer1 = 2 * genocounts[0] + genocounts[1];
          const uint32_t numer2 = genocounts[1] + 2 * genocounts[2];
          if ((!numer1) || (!numer2)) {
            continue;
          }
          const uint32_t denom = numer1 + numer2;
          ehet = 2 * S_CAST(double, numer1) * S_CAST(double, numer2) / ((S_CAST(double, denom)) * S_CAST(double, denom - 1));
        }
      } else {
        // multiallelic
        if (allele_freqs) {
          const double* cur_allele_freqs = &(allele_freqs[allele_idx_offset_base - variant_uidx]);
          const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
          double freq_sum = 0.0;
          double freq_ssq = 0.0;
          for (uint32_t aidx = 0; aidx != cur_allele_ct_m1; ++aidx) {
            const double cur_freq = cur_allele_freqs[aidx];
            freq_sum += cur_freq;
            freq_ssq += cur_freq * cur_freq;
          }
          const double last_allele_freq = 1.0 - freq_sum;
          ehet = 1.0 - freq_ssq - last_allele_freq * last_allele_freq;
          if (ehet < kSmallishEpsilon) {
            continue;
          }
        }
        PglErr reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, pgrp, &pgv);
        if (unlikely(reterr)) {
          ctx->reterr = reterr;
          break;
        }
        const uint32_t patch_10_ct = pgv.patch_10_ct;
        const AlleleCode* patch_10_vals = pgv.patch_10_vals;
        if (!allele_freqs) {
          STD_ARRAY_DECL(uint32_t, 4, genocounts);
          GenoarrCountSubsetFreqs(genovec, founder_info_interleaved_vec, sample_ct, founder_ct, genocounts);
          allele_nobs[0] = 2 * genocounts[0] + genocounts[1];
          allele_nobs[1] = genocounts[1] + 2 * genocounts[2];
          const uint32_t denom = allele_nobs[0] + allele_nobs[1];
          if (!denom) {
            continue;
          }
          ZeroU32Arr(cur_allele_ct - 2, &(allele_nobs[2]));
          const uint32_t patch_01_ct = pgv.patch_01_ct;
          if (patch_01_ct) {
            const uintptr_t* patch_01_set = pgv.patch_01_set;
            const AlleleCode* patch_01_vals = pgv.patch_01_vals;
            uintptr_t sample_idx_base = 0;
            uintptr_t sample_idx_bits = patch_01_set[0];
            for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &sample_idx_bits);
              allele_nobs[patch_01_vals[uii]] += IsSet(founder_info_collapsed, sample_idx);
            }
            allele_nobs[1] -= PopcountWordsIntersect(founder_info_collapsed, patch_01_set, sample_ctl);
          }
          if (patch_10_ct) {
            const uintptr_t* patch_10_set = pgv.patch_10_set;
            uintptr_t sample_idx_base = 0;
            uintptr_t sample_idx_bits = patch_10_set[0];
            for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
              const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
              if (IsSet(founder_info_collapsed, sample_idx)) {
                allele_nobs[patch_10_vals[2 * uii]] += 1;
                allele_nobs[patch_10_vals[2 * uii + 1]] += 1;
              }
            }
            allele_nobs[1] -= 2 * PopcountWordsIntersect(founder_info_collapsed, patch_10_set, sample_ctl);
          }

          uint64_t allele_nobs_ssq = 0;
          for (uint32_t uii = 0; uii != cur_allele_ct; ++uii) {
            const uintptr_t cur_allele_nobs = allele_nobs[uii];
            allele_nobs_ssq += S_CAST(uint64_t, cur_allele_nobs) * cur_allele_nobs;
          }
          const double denom_d = S_CAST(double, denom);
          ehet = (1.0 - (S_CAST(double, allele_nobs_ssq) / (denom_d * denom_d))) * (denom_d / S_CAST(double, denom - 1));
        }
        if (patch_10_ct) {
          // For every altx/alty genotype where x and y are different, change
          // the genovec entry from 2 to 1.
          const uintptr_t* patch_10_set = pgv.patch_10_set;
          uintptr_t sample_idx_base = 0;
          uintptr_t sample_idx_bits = patch_10_set[0];
          for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
            const AlleleCode a0 = patch_10_vals[uii * 2];
            const AlleleCode a1 = patch_10_vals[uii * 2 + 1];
            if (a0 != a1) {
              const uintptr_t widx = sample_idx / kBitsPerWordD2;
              const uint32_t bit_shift_ct = 2 * (sample_idx % kBitsPerWordD2);
              genovec[widx] ^= (3 * k1LU) << bit_shift_ct;
            }
          }
        }
      }
      ehet_base += ehet;
      ++nobs_base;

      // See e.g. UpdateDenseSampleCounts2().
      const VecW m1 = VCONST_W(kMask5555);
      const VecW* genovvec = R_CAST(const VecW*, genovec);
      for (uint32_t vidx = 0; vidx != acc2_vec_ct; ++vidx) {
        const VecW vv = genovvec[vidx];
        const VecW vv_shifted = vecw_srli(vv, 1);
        const VecW vv_lo_only = vv & m1;
        scrambled_ohets[vidx] += vecw_and_notfirst(vv_shifted, vv_lo_only);
      }
      --dense_ct_rem3;
      if (!dense_ct_rem3) {
        VecW* acc4 = &(scrambled_ohets[acc2_vec_ct]);
        Vcount0Incr2To4(acc2_vec_ct, scrambled_ohets, acc4);
        --dense_ct_rem15d3;
        if (!dense_ct_rem15d3) {
          const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
          VecW* acc8 = &(acc4[acc4_vec_ct]);
          Vcount0Incr4To8(acc4_vec_ct, acc4, acc8);
          --dense_ct_rem255d15;
          if (!dense_ct_rem255d15) {
            const uint32_t acc8_vec_ct = acc4_vec_ct * 2;
            VecW* acc32 = &(acc8[acc8_vec_ct]);
            Vcount0Incr8To32(acc8_vec_ct, acc8, acc32);
            dense_ct_rem255d15 = 17;
          }
          dense_ct_rem15d3 = 5;
        }
        dense_ct_rem3 = 3;
      }

      for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
        const uintptr_t geno_word = genovec[widx];
        uintptr_t missing_word = geno_word & (geno_word >> 1) & kMask5555;
        if (!missing_word) {
          continue;
        }
        uint32_t sample_idx_base = widx * kBitsPerWordD2;
        double* cur_ehet_incrs = &(ehet_incrs[sample_idx_base]);
        int32_t* cur_nobs_incrs = &(nobs_incrs[sample_idx_base]);
        do {
          const uint32_t sample_idx_lowbits = ctzw(missing_word) / 2;
          cur_ehet_incrs[sample_idx_lowbits] -= ehet;
          cur_nobs_incrs[sample_idx_lowbits] -= 1;
          missing_word &= missing_word - 1;
        } while (missing_word);
      }
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  ctx->thread_ehet_base[tidx] = ehet_base;
  ctx->thread_nobs_base[tidx] = nobs_base;

  VecW* acc4 = &(scrambled_ohets[acc2_vec_ct]);
  VcountIncr2To4(scrambled_ohets, acc2_vec_ct, acc4);
  const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
  VecW* acc8 = &(acc4[acc4_vec_ct]);
  VcountIncr4To8(acc4, acc4_vec_ct, acc8);
  const uint32_t acc8_vec_ct = acc4_vec_ct * 2;
  VecW* acc32 = &(acc8[acc8_vec_ct]);
  VcountIncr8To32(acc8, acc8_vec_ct, acc32);
  uint32_t* acc32_alias = R_CAST(uint32_t*, acc32);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uint32_t scrambled_idx = VcountScramble2(sample_idx);
    ohets[sample_idx] += acc32_alias[scrambled_idx];
  }
  THREAD_RETURN;
}

PglErr HetReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const double* allele_freqs, const uintptr_t* founder_info, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_allele_ct, HetFlags flags, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  ThreadGroup tg;
  PreinitThreads(&tg);
  HetCtx ctx;
  {
    if (IsSet(cip->haploid_mask, 0)) {
      logerrputs("Error: --het cannot be used on haploid genomes.\n");
      goto HetReport_ret_INCONSISTENT_INPUT;
    }
    const uint32_t small_sample = (flags / kfHetSmallSample) & 1;
    if (small_sample && (!founder_ct)) {
      logerrputs("Error: '--het small-sample' requires founders.  (PLINK 1.9 --make-founders may\ncome in handy here.\n)");
      goto HetReport_ret_INCONSISTENT_INPUT;
    }

    ctx.autosomal_variant_include = orig_variant_include;
    uint32_t autosomal_variant_ct = orig_variant_ct;
    reterr = ConditionalAllocateNonAutosomalVariants(cip, "--het", raw_variant_ct, &ctx.autosomal_variant_include, &autosomal_variant_ct);
    if (!autosomal_variant_ct) {
      goto HetReport_ret_NO_VARIATION;
    }

    const uint32_t output_zst = flags & kfHetZs;
    OutnameZstSet(".het", output_zst, outname_end);
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, 2 * kCompressStreamBlock, &css, &cswritep);
    if (unlikely(reterr)) {
      goto HetReport_ret_1;
    }
    *cswritep++ = '#';
    const uint32_t col_fid = FidColIsRequired(siip, flags / kfHetColMaybefid);
    if (col_fid) {
      cswritep = strcpya_k(cswritep, "FID\t");
    }
    cswritep = strcpya_k(cswritep, "IID");
    const uint32_t col_sid = SidColIsRequired(siip->sids, flags / kfHetColMaybesid);
    if (col_sid) {
      cswritep = strcpya_k(cswritep, "\tSID");
    }
    const uint32_t col_hom = (flags / kfHetColHom) & 1;
    if (col_hom) {
      cswritep = strcpya_k(cswritep, "\tO(HOM)\tE(HOM)");
    }
    const uint32_t col_het = (flags / kfHetColHet) & 1;
    if (col_het) {
      cswritep = strcpya_k(cswritep, "\tO(HET)\tE(HET)");
    }
    const uint32_t col_nobs = (flags / kfHetColNobs) & 1;
    if (col_nobs) {
      cswritep = strcpya_k(cswritep, "\tOBS_CT");
    }
    const uint32_t col_f = (flags / kfHetColF) & 1;
    if (col_f) {
      cswritep = strcpya_k(cswritep, "\tF");
    }
    AppendBinaryEoln(&cswritep);
    if (unlikely(Cswrite(&css, &cswritep))) {
      goto HetReport_ret_WRITE_FAIL;
    }

    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.sample_include = sample_include;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uint32_t* sample_include_cumulative_popcounts;
    if (unlikely(
            bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts))) {
      goto HetReport_ret_NOMEM;
    }
    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
    ctx.founder_info_collapsed = nullptr;
    ctx.founder_info_interleaved_vec = nullptr;
    ctx.allele_freqs = allele_freqs;
    if (small_sample) {
      const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
      uintptr_t* founder_info_collapsed;
      uintptr_t* founder_info_interleaved_vec;
      if (unlikely(
              bigstack_alloc_w(sample_ctaw, &founder_info_collapsed) ||
              bigstack_alloc_w(sample_ctaw, &founder_info_interleaved_vec))) {
        goto HetReport_ret_NOMEM;
      }
      CopyBitarrSubset(founder_info, sample_include, sample_ct, founder_info_collapsed);
      const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
      ZeroWArr(sample_ctaw - sample_ctl, &(founder_info_collapsed[sample_ctl]));
      FillInterleavedMaskVec(founder_info_collapsed, sample_ctaw / kWordsPerVec, founder_info_interleaved_vec);
      ctx.founder_info_collapsed = founder_info_collapsed;
      ctx.founder_info_interleaved_vec = founder_info_interleaved_vec;
      ctx.allele_freqs = nullptr;
    }
    ctx.sample_ct = sample_ct;
    ctx.founder_ct = founder_ct;

    uint32_t calc_thread_ct = max_thread_ct;
    if (unlikely(
            bigstack_alloc_wp(calc_thread_ct, &ctx.raregenos) ||
            bigstack_alloc_u32p(calc_thread_ct, &ctx.difflist_sample_id_bufs) ||
            bigstack_alloc_vp(calc_thread_ct, &ctx.scrambled_ohet_bufs) ||
            bigstack_alloc_d(calc_thread_ct, &ctx.thread_ehet_base) ||
            bigstack_alloc_u32(calc_thread_ct, &ctx.thread_nobs_base) ||
            bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_ohets) ||
            bigstack_alloc_dp(calc_thread_ct, &ctx.thread_ehet_incrs) ||
            bigstack_alloc_i32p(calc_thread_ct, &ctx.thread_nobs_incrs))) {
      goto HetReport_ret_NOMEM;
    }
    const uint32_t max_returned_difflist_len = 2 * (raw_sample_ct / kPglMaxDifflistLenDivisor);
    const uintptr_t raregeno_vec_ct = DivUp(max_returned_difflist_len, kNypsPerVec);
    const uintptr_t difflist_sample_id_vec_ct = DivUp(max_returned_difflist_len, kInt32PerVec);
    const uint32_t mhc_needed = (max_allele_ct > 2);
    uintptr_t allele_nobs_vec_ct = 0;
    ctx.allele_nobs_bufs = nullptr;
    if (mhc_needed) {
      if (unlikely(
              bigstack_alloc_u32p(calc_thread_ct, &ctx.allele_nobs_bufs))) {
        goto HetReport_ret_NOMEM;
      }
      allele_nobs_vec_ct = DivUp(max_allele_ct, kInt32PerVec);
    }

    const uintptr_t acc2_vec_ct = NypCtToVecCt(sample_ct);
    const uintptr_t scrambled_ohet_vec_ct = acc2_vec_ct * 23;
    const uintptr_t sample_ct_i32v = DivUp(sample_ct, kInt32PerVec);
    const uintptr_t sample_ct_dv = DivUp(sample_ct * sizeof(double), kBytesPerVec);

    const uintptr_t thread_xalloc_vec_ct = raregeno_vec_ct + difflist_sample_id_vec_ct + allele_nobs_vec_ct + scrambled_ohet_vec_ct + 2 * sample_ct_i32v + sample_ct_dv;
    const uintptr_t thread_xalloc_cacheline_ct = DivUp(thread_xalloc_vec_ct, kVecsPerCacheline);
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    ctx.thread_read_mhc = nullptr;
    uint32_t read_block_size;
    if (unlikely(PgenMtLoadInit(ctx.autosomal_variant_include, raw_sample_ct, autosomal_variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &ctx.genovecs, mhc_needed? (&ctx.thread_read_mhc) : nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto HetReport_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto HetReport_ret_NOMEM;
    }
    ctx.reterr = kPglRetSuccess;
    // assert(bigstack_left() >= thread_xalloc_cacheline_ct * kCacheline * calc_thread_ct);
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      unsigned char* cur_alloc = S_CAST(unsigned char*, bigstack_alloc_raw(thread_xalloc_cacheline_ct * kCacheline));
      ctx.raregenos[tidx] = R_CAST(uintptr_t*, cur_alloc);
      cur_alloc = &(cur_alloc[raregeno_vec_ct * kBytesPerVec]);
      ctx.difflist_sample_id_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
      cur_alloc = &(cur_alloc[difflist_sample_id_vec_ct * kBytesPerVec]);
      if (mhc_needed) {
        ctx.allele_nobs_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[allele_nobs_vec_ct * kBytesPerVec]);
      }
      ctx.scrambled_ohet_bufs[tidx] = R_CAST(VecW*, cur_alloc);
      cur_alloc = &(cur_alloc[scrambled_ohet_vec_ct * kBytesPerVec]);
      ctx.thread_ohets[tidx] = R_CAST(uint32_t*, cur_alloc);
      cur_alloc = &(cur_alloc[sample_ct_i32v * kBytesPerVec]);
      ctx.thread_ehet_incrs[tidx] = R_CAST(double*, cur_alloc);
      cur_alloc = &(cur_alloc[sample_ct_dv * kBytesPerVec]);
      ctx.thread_nobs_incrs[tidx] = R_CAST(int32_t*, cur_alloc);
      // cur_alloc = &(cur_alloc[sample_ct_i32v * kBytesPerVec]);
    }
    SetThreadFuncAndData(HetThread, &ctx, &tg);

    logputs("--het: ");
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;

    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t next_print_variant_idx = autosomal_variant_ct / 100;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_size = MultireadNonempty(ctx.autosomal_variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto HetReport_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = ctx.reterr;
        if (unlikely(reterr)) {
          goto HetReport_ret_PGR_FAIL;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_size = cur_block_size;
        ComputeUidxStartPartition(ctx.autosomal_variant_include, cur_block_size, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_size == autosomal_variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto HetReport_ret_THREAD_CREATE_FAIL;
        }
      }

      parity = 1 - parity;
      if (variant_idx == autosomal_variant_ct) {
        break;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / autosomal_variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, autosomal_variant_ct)) / 100;
      }

      ++read_block_idx;
      variant_idx += cur_block_size;
      pgfip->block_base = main_loadbufs[parity];
    }
    double ehet_base = ctx.thread_ehet_base[0];
    uint32_t nobs_base = ctx.thread_nobs_base[0];
    uint32_t* ohets = ctx.thread_ohets[0];
    double* ehet_incrs = ctx.thread_ehet_incrs[0];
    int32_t* nobs_incrs = ctx.thread_nobs_incrs[0];
    for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
      ehet_base += ctx.thread_ehet_base[tidx];
      nobs_base += ctx.thread_nobs_base[tidx];
      U32CastVecAdd(ctx.thread_ohets[tidx], sample_ct_i32v, ohets);
      const double* src = ctx.thread_ehet_incrs[tidx];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        ehet_incrs[sample_idx] += src[sample_idx];
      }
      I32CastVecAdd(ctx.thread_nobs_incrs[tidx], sample_ct_i32v, nobs_incrs);
    }
    // Ready to write results.
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      const char* fid_or_iid_start = &(sample_ids[sample_uidx * max_sample_id_blen]);
      if (!col_fid) {
        fid_or_iid_start = AdvPastDelim(fid_or_iid_start, '\t');
      }
      cswritep = strcpya(cswritep, fid_or_iid_start);

      if (col_sid) {
        *cswritep++ = '\t';
        if (sids) {
          cswritep = strcpya(cswritep, &(sids[sample_uidx * max_sid_blen]));
        } else {
          *cswritep++ = '0';
        }
      }
      const uint32_t ohet = ohets[sample_idx];
      const double ehet = ehet_base + ehet_incrs[sample_idx];
      const uint32_t nobs = nobs_base + nobs_incrs[sample_idx];
      if (col_hom) {
        *cswritep++ = '\t';
        cswritep = u32toa_x(nobs - ohet, '\t', cswritep);
        cswritep = dtoa_g(u31tod(nobs) - ehet, cswritep);
      }
      if (col_het) {
        *cswritep++ = '\t';
        cswritep = u32toa_x(ohet, '\t', cswritep);
        cswritep = dtoa_g(ehet, cswritep);
      }
      if (col_nobs) {
        *cswritep++ = '\t';
        cswritep = u32toa(nobs, cswritep);
      }
      if (col_f) {
        *cswritep++ = '\t';
        cswritep = dtoa_g(1.0 - ohet / ehet, cswritep);
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto HetReport_ret_WRITE_FAIL;
      }
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto HetReport_ret_WRITE_FAIL;
    }
    logprintfww("--het: Results written to %s .\n", outname);
  }
  while (0) {
  HetReport_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  HetReport_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  HetReport_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  HetReport_ret_NO_VARIATION:
    logerrputs("Error: --het requires at least one polymorphic autosomal variant.\n");
  HetReport_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  HetReport_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 HetReport_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct FstCtxStruct {
  const uintptr_t* autosomal_variant_include;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  const uint32_t* sample_idx_to_pop_idx;
  uint32_t sample_ct;
  uint32_t pop_ct;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** raregenos;
  uint32_t** difflist_sample_id_bufs;
  uint32_t** pop_allele_nobs_bufs;  // major dimension = pop, minor = allele
  uint32_t** pop_het_ct_bufs;  // nullptr if method=hudson
  uint32_t* read_variant_uidx_starts;
  uint32_t cur_block_size;  // not for jackknife

  // only kPglRetMalformedInput possible, no atomic ops needed
  PglErr reterr;

  // worker threads compute intermediate stats; parent thread computes final
  // Fst estimates and executes block-jackknife.
  uint32_t* variant_obs_cts[2];
  double* variant_[2];

  double* thread_ehet_base;
  uint32_t* thread_nobs_base;
  uint32_t** thread_ohets;
  double** thread_ehet_incrs;
  int32_t** thread_nobs_incrs;
} FstCtx;

/*
PglErr FstReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const FstInfo* fst_infop, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end) {
  // TODO
  return kPglRetSuccess;
}
*/

#ifdef __cplusplus
}  // namespace plink2
#endif
