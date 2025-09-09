// This file is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
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

#include "plink2_misc.h"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_stats.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "plink2_cmdline.h"
#include "plink2_compress_stream.h"
#include "plink2_decompress.h"
#include "plink2_data.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitUpdateAlleles(UpdateAllelesInfo* update_alleles_info_ptr) {
  update_alleles_info_ptr->flags = kfUpdateAlleles0;
  update_alleles_info_ptr->fname = nullptr;
}

void CleanupUpdateAlleles(UpdateAllelesInfo* update_alleles_info_ptr) {
  free_cond(update_alleles_info_ptr->fname);
}

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
  fst_info_ptr->pheno_name = nullptr;
  fst_info_ptr->first_id_or_fname = nullptr;
  fst_info_ptr->other_ids_flattened = nullptr;
}

void CleanupFst(FstInfo* fst_info_ptr) {
  free_cond(fst_info_ptr->pheno_name);
  free_cond(fst_info_ptr->first_id_or_fname);
  free_cond(fst_info_ptr->other_ids_flattened);
}

void InitCheckSex(CheckSexInfo* check_sex_info_ptr) {
  check_sex_info_ptr->flags = kfCheckSex0;
  check_sex_info_ptr->max_female_xf = -1.0;
  check_sex_info_ptr->min_male_xf = -1.0;
  check_sex_info_ptr->max_female_ycount = UINT32_MAX;
  check_sex_info_ptr->min_male_ycount = UINT32_MAX;
  check_sex_info_ptr->max_female_yrate = -1.0;
  check_sex_info_ptr->min_male_yrate = -1.0;
}

PglErr FlipAlleles(const uintptr_t* variant_include, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const uintptr_t* allele_idx_offsets, const FlipInfo* flip_info_ptr, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t max_thread_ct, char** allele_storage_mutable) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    if (!variant_ct) {
      goto FlipAlleles_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* flip_include;
    if (unlikely(bigstack_alloc_w(raw_variant_ctl, &flip_include))) {
      goto FlipAlleles_ret_NOMEM;
    }
    reterr = LoadTokensNondup2(flip_info_ptr->fname, variant_include, variant_ids, variant_id_htable, htable_dup_base, "flip", raw_variant_ct, max_variant_id_slen, variant_id_htable_size, max_thread_ct, flip_include);
    if (unlikely(reterr)) {
      goto FlipAlleles_ret_1;
    }

    const uint32_t flip_ct = PopcountWords(flip_include, raw_variant_ctl);
    const uint32_t permissive = (flip_info_ptr->flags / kfFlipPermissive) & 1;
    char* dash_ptr = K_CAST(char*, &(g_one_char_strs[90]));
    // dash = ASCII 45.  relevant allele codes are '.' (46), 'A' (65), 'C'
    // (67), 'G' (71), 'N' (78), 'T' (84), and lowercase letters (+32).
    // When table entry is nonzero, flipped allele code is
    // &(dash_ptr[table entry]).
    uint8_t dash_ptrdiff_table[144] = {
      0, 0, 2, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 78, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 44, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 66, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 40, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 78, 0, 0, 0, 52, 0, 0, 0, 0, 0, 0, 0, 44, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 66, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 40, 0};
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = flip_include[0];
    uint32_t cur_allele_ct = 2;
    uint32_t nonsnp_ct = 0;
    // possible todo: multithread this part
    for (uint32_t flip_idx = 0; flip_idx != flip_ct; ++flip_idx) {
      const uint32_t variant_uidx = BitIter1(flip_include, &variant_uidx_base, &cur_bits);
      uintptr_t allele_idx_offset_base;
      if (!allele_idx_offsets) {
        allele_idx_offset_base = 2 * variant_uidx;
      } else {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      char** cur_alleles = &(allele_storage_mutable[allele_idx_offset_base]);
      uint32_t non_acgtn_found = 0;
      for (uint32_t aidx = 0; aidx != cur_allele_ct; ++aidx) {
        const uintptr_t ptr_diff = cur_alleles[aidx] - dash_ptr;
        if ((ptr_diff > 144) || (!dash_ptrdiff_table[ptr_diff])) {
          non_acgtn_found = 1;
          break;
        }
      }
      if (non_acgtn_found) {
        if (unlikely(!permissive)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Variant '%s' in --flip file is not a SNP. (Add the 'permissive' modifier to skip instead of erroring out.)\n", variant_ids[variant_uidx]);
          goto FlipAlleles_ret_INCONSISTENT_INPUT_WW;
        }
        ++nonsnp_ct;
        continue;
      }
      for (uint32_t aidx = 0; aidx != cur_allele_ct; ++aidx) {
        const uintptr_t ptr_diff = cur_alleles[aidx] - dash_ptr;
        cur_alleles[aidx] = &(dash_ptr[dash_ptrdiff_table[ptr_diff]]);
      }
    }
    if (!nonsnp_ct) {
      logprintf("--flip: %u variant%s flipped.\n", flip_ct, (flip_ct == 1)? "" : "s");
    } else {
      logprintf("--flip: %u variant%s flipped, %u non-SNP%s skipped.\n", flip_ct - nonsnp_ct, (flip_ct - nonsnp_ct == 1)? "" : "s", nonsnp_ct, (nonsnp_ct == 1)? "" : "s");
    }
  }
  while (0) {
  FlipAlleles_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  FlipAlleles_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 FlipAlleles_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr UpdateVarBps(const ChrInfo* cip, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const TwoColParams* params, uint32_t sort_vars_in_cmd, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uint32_t htable_size, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* __restrict variant_bps, uint32_t* __restrict variant_ct_ptr, UnsortedVar* vpos_sortstatusp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    uintptr_t* already_seen;
    if (unlikely(bigstack_calloc_w(BitCtToWordCt(raw_variant_ct), &already_seen))) {
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
            if (sort_vars_in_cmd) {
              logerrputs("Note: Base-pair positions are now unsorted.\n");
            } else {
              logerrputs("Warning: Base-pair positions are now unsorted!\n");
            }
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
    if (unlikely(bigstack_alloc_cp(raw_variant_ct, &variant_ids_copy) ||
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
        // bugfix (7 Jun 2023): this condition was flipped
        if (likely(!TextStreamErrcode2(&txs, &reterr))) {
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

static_assert(kPglMaxAlleleCt <= 65535, "UpdateVarAlleles() must be updated.");
PglErr UpdateVarAlleles(const uintptr_t* variant_include, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const uintptr_t* allele_idx_offsets, const UpdateAllelesInfo* update_alleles_info_ptr, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uint32_t htable_size, uint32_t max_allele_ct, char input_missing_geno_char, uint32_t max_thread_ct, char** allele_storage_mutable, uint32_t* max_allele_slen_ptr, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  // Set this to 65536 even though kPglMaxAlleleCt is only 255 as of this
  // writing, since VCF variants with >255 alleles are out there and we don't
  // want to crash when partial-match is a possibility.
  const uint32_t max_impl_input_allele_ct = 65536;
  uintptr_t line_idx = 0;
  FILE* errfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    const uint32_t max_allele_ctl = BitCtToWordCt(max_allele_ct);
    char** input_allele_starts;
    uint32_t* input_allele_slens;
    const char** orig_alleles_sorted;
    uint32_t* sorted_to_orig_idx;
    uint16_t* orig_to_input_idx;
    uintptr_t* alleles_seen;
    uintptr_t* variants_seen;
    unsigned char* sort_wkspace;
    if (unlikely(bigstack_alloc_cp(max_impl_input_allele_ct, &input_allele_starts) ||
                 bigstack_alloc_u32(max_impl_input_allele_ct, &input_allele_slens) ||
                 bigstack_alloc_kcp(max_allele_ct, &orig_alleles_sorted) ||
                 bigstack_alloc_u32(max_allele_ct, &sorted_to_orig_idx) ||
                 bigstack_alloc_u16(max_allele_ct, &orig_to_input_idx) ||
                 bigstack_alloc_w(max_allele_ctl, &alleles_seen) ||
                 bigstack_calloc_w(BitCtToWordCt(raw_variant_ct), &variants_seen) ||
                 bigstack_alloc_uc(max_allele_ct * sizeof(StrSortIndexedDeref), &sort_wkspace))) {
      goto UpdateVarAlleles_ret_NOMEM;
    }
    reterr = SizeAndInitTextStream(update_alleles_info_ptr->fname, bigstack_left() / 4, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto UpdateVarAlleles_ret_TSTREAM_FAIL;
    }
    const uint32_t allow_mismatch = (update_alleles_info_ptr->flags / kfUpdateAllelesAllowMismatch) & 1;
    const uint32_t strict_missing = (update_alleles_info_ptr->flags / kfUpdateAllelesStrictMissing) & 1;
    const char* std_input_missing_geno = &(g_one_char_strs[92]);
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = g_bigstack_end;
    uintptr_t miss_ct = 0;
    uintptr_t err_ct = 0;
    uint32_t hit_ct = 0;
    uint32_t is_3col = 2; // 0 = orig format, 1 = 3col, 2 = not determined yet
    uint32_t max_allele_slen = *max_allele_slen_ptr;
    uint32_t cur_input_allele_ct = 2;
    uint32_t cur_allele_ct = 2;
    char* line_iter = TextLineEnd(&txs);
    ++line_idx;
    for (; TextGetUnsafe2(&txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      char* varid_start = line_iter;
      char* varid_end = CurTokenEnd(varid_start);
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
      if (unlikely(IsSet(variants_seen, variant_uidx))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --update-alleles file.\n", varid);
        goto UpdateVarAlleles_ret_INCONSISTENT_INPUT_WW;
      }
      SetBit(variant_uidx, variants_seen);
      if (!IsSet(variant_include, variant_uidx)) {
        line_iter = varid_end;
        ++miss_ct;
        continue;
      }

      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      char** orig_alleles = &(allele_storage_mutable[allele_idx_offset_base]);

      char* col2_start = FirstNonTspace(varid_end);
      char* col2_end = FirstSpaceOrEoln(col2_start);
      char* col3_start = FirstNonTspace(col2_end);
      if (unlikely(IsEolnKns(*col3_start))) {
        goto UpdateVarAlleles_ret_MISSING_TOKENS;
      }
      char* col3_end = FirstSpaceOrEoln(col3_start);
      if (is_3col == 2) {
        is_3col = (NextToken(col3_end) == nullptr);
      }
      char* new_alleles_start = col3_start;
      if (is_3col) {
        cur_input_allele_ct = 0;
        for (char* col2_iter = col2_start; ; ) {
          char* input_allele_end = AdvToDelimOrEnd(col2_iter, col2_end, ',');
          input_allele_starts[cur_input_allele_ct] = col2_iter;
          input_allele_slens[cur_input_allele_ct] = input_allele_end - col2_iter;
          ++cur_input_allele_ct;
          // Null-terminate to set up bsearch_strptr_overread().
          *input_allele_end = '\0';
          if (input_allele_end == col2_end) {
            break;
          }
          if (unlikely(cur_input_allele_ct == max_impl_input_allele_ct)) {
            goto UpdateVarAlleles_ret_TOO_MANY_ALLELES;
          }
          col2_iter = &(input_allele_end[1]);
        }
        if ((cur_input_allele_ct > cur_allele_ct) && (!allow_mismatch)) {
        UpdateVarAlleles_errfile:
          if (!err_ct) {
            strcpy_k(outname_end, ".allele.no.snp");
            if (fopen_checked(outname, FOPEN_WB, &errfile)) {
              goto UpdateVarAlleles_ret_OPEN_FAIL;
            }
          }
          // variant ID, expected allele 1, remaining comma-separated expected
          // alleles
          // not optimized for now, could explicitly manage a write buffer if
          // it ever matters
          fputs(varid, errfile);
          putc_unlocked('\t', errfile);
          fwrite(input_allele_starts[0], input_allele_slens[0], 1, errfile);
          putc_unlocked('\t', errfile);
          if (cur_input_allele_ct == 1) {
            putc_unlocked('.', errfile);
          } else {
            for (uint32_t aidx = 1; ; ) {
              fwrite(input_allele_starts[aidx], input_allele_slens[aidx], 1, errfile);
              ++aidx;
              if (aidx == cur_input_allele_ct) {
                break;
              }
              putc_unlocked(',', errfile);
            }
          }
#ifdef _WIN32
          putc_unlocked('\r', errfile);
#endif
          if (unlikely(putc_checked('\n', errfile))) {
            goto UpdateVarAlleles_ret_WRITE_FAIL;
          }
          ++err_ct;
          continue;
        }
      } else {
        // cur_input_allele_ct == 2
        new_alleles_start = FirstNonTspace(col3_end);
        input_allele_starts[0] = col2_start;
        input_allele_starts[1] = col3_start;
        input_allele_slens[0] = col2_end - col2_start;
        input_allele_slens[1] = col3_end - col3_start;
        *col2_end = '\0';
        *col3_end = '\0';
      }
      uint32_t wildcard_idx = UINT32_MAX;
      // safe to assume preexisting allele codes are distinct since
      // CheckAlleleUniqueness() was run.
      if (cur_allele_ct == 2) {
        if (strcmp_overread_lt(orig_alleles[0], orig_alleles[1])) {
          sorted_to_orig_idx[0] = 0;
          sorted_to_orig_idx[1] = 1;
        } else {
          sorted_to_orig_idx[0] = 1;
          sorted_to_orig_idx[1] = 0;
        }
        orig_alleles_sorted[0] = orig_alleles[sorted_to_orig_idx[0]];
        orig_alleles_sorted[1] = orig_alleles[sorted_to_orig_idx[1]];
        if (!strict_missing) {
          for (uint32_t uii = 0; uii != 2; ++uii) {
            if (memequal(orig_alleles_sorted[uii], std_input_missing_geno, 2)) {
              wildcard_idx = sorted_to_orig_idx[uii];
              break;
            }
          }
        }
      } else {
        memcpy(orig_alleles_sorted, orig_alleles, cur_allele_ct * sizeof(intptr_t));
        SortStrptrArrIndexed2(cur_allele_ct, 0, 1, 0, orig_alleles_sorted, sorted_to_orig_idx, nullptr, sort_wkspace);
      }
      ZeroWArr(max_allele_ctl, alleles_seen);
      for (uint32_t old_idx = 0; old_idx != cur_input_allele_ct; ++old_idx) {
        const char* old_allele = input_allele_starts[old_idx];
        uint32_t old_slen = input_allele_slens[old_idx];
        // standardize missing allele codes
        if ((old_slen == 1) && (old_allele[0] == input_missing_geno_char)) {
          old_allele = std_input_missing_geno;
        }
        const int32_t sorted_orig_idx = bsearch_strptr_overread(old_allele, orig_alleles_sorted, cur_allele_ct);
        if (sorted_orig_idx != -1) {
          const uint32_t orig_idx = sorted_to_orig_idx[sorted_orig_idx];
          if (IsSet(alleles_seen, orig_idx)) {
            goto UpdateVarAlleles_ret_DUPLICATE_INPUT_ALLELE_CODE;
          }
          SetBit(orig_idx, alleles_seen);
          orig_to_input_idx[orig_idx] = old_idx;
        }
      }
      if ((wildcard_idx != UINT32_MAX) && (alleles_seen[0] == 2 - wildcard_idx) && (cur_input_allele_ct == 2)) {
        // Biallelic variant with one missing allele code, only non-missing
        // allele has been matched, exactly two input alleles.
        // Match the other old input allele with the original missing allele.
        alleles_seen[0] = 3;
        orig_to_input_idx[wildcard_idx] = 1 - orig_to_input_idx[1 - wildcard_idx];
      }
      const uint32_t seen_ct = PopcountWords(alleles_seen, max_allele_ctl);
      if ((!seen_ct) || ((!allow_mismatch) && (seen_ct < cur_input_allele_ct))) {
        goto UpdateVarAlleles_errfile;
      }

      if (is_3col) {
        uint32_t new_input_allele_idx = 0;
        for (char* col3_iter = col3_start; ; ) {
          char* input_allele_end = AdvToDelimOrEnd(col3_iter, col3_end, ',');
          input_allele_starts[new_input_allele_idx] = col3_iter;
          input_allele_slens[new_input_allele_idx] = input_allele_end - col3_iter;
          // No need to null-terminate this time.  (And if we did, we'd need to
          // handle the case where *col3_end was '\n'.)
          ++new_input_allele_idx;
          if (input_allele_end == col3_end) {
            if (unlikely(new_input_allele_idx != cur_input_allele_ct)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Too few new alleles on line %" PRIuPTR " of --update-alleles file.\n", line_idx);
              goto UpdateVarAlleles_ret_MALFORMED_INPUT_WW;
            }
            break;
          }
          if (unlikely(new_input_allele_idx == cur_input_allele_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Too many new alleles on line %" PRIuPTR " of --update-alleles file.\n", line_idx);
            goto UpdateVarAlleles_ret_MALFORMED_INPUT_WW;
          }
          col3_iter = &(input_allele_end[1]);
        }
      } else {
        char* col4_end = FirstSpaceOrEoln(new_alleles_start);
        char* col5_start = FirstNonTspace(col4_end);
        if (unlikely(IsEolnKns(*col5_start))) {
          goto UpdateVarAlleles_ret_MISSING_TOKENS;
        }
        char* col5_end = FirstSpaceOrEoln(col5_start);
        input_allele_starts[0] = new_alleles_start;
        input_allele_starts[1] = col5_start;
        input_allele_slens[0] = col4_end - new_alleles_start;
        input_allele_slens[1] = col5_end - col5_start;
      }

      uintptr_t aidx_base = 0;
      uintptr_t cur_bits = alleles_seen[0];
      for (uint32_t seen_idx = 0; seen_idx != seen_ct; ++seen_idx) {
        const uint32_t aidx = BitIter1(alleles_seen, &aidx_base, &cur_bits);
        const uint32_t new_idx = orig_to_input_idx[aidx];
        const char* new_allele_start = input_allele_starts[new_idx];
        const uint32_t new_slen = input_allele_slens[new_idx];
        if (new_slen == 1) {
          char cc = new_allele_start[0];
          if (cc == input_missing_geno_char) {
            cc = '.';
          }
          orig_alleles[aidx] = K_CAST(char*, &(g_one_char_strs[ctou32(cc) * 2]));
        } else {
          // reuse old storage if we can, allocate when we must
          if (strlen(orig_alleles[aidx]) < new_slen) {
            if (PtrWSubCk(tmp_alloc_base, new_slen + 1, &tmp_alloc_end)) {
              goto UpdateVarAlleles_ret_NOMEM;
            }
            orig_alleles[aidx] = R_CAST(char*, tmp_alloc_end);
          }
          memcpyx(orig_alleles[aidx], new_allele_start, new_slen, '\0');
        }
      }

      // Confirm there are no duplicate allele codes after the change.
      if (cur_allele_ct == 2) {
        if (unlikely(strequal_overread(orig_alleles[0], orig_alleles[1]))) {
          goto UpdateVarAlleles_ret_DUPLICATE_ALLELE_CODE_AFTER_UPDATE;
        }
      } else {
        memcpy(orig_alleles_sorted, orig_alleles, cur_allele_ct * sizeof(intptr_t));
        SortStrptrArrIndexed2(cur_allele_ct, 0, 1, 0, orig_alleles_sorted, sorted_to_orig_idx, nullptr, sort_wkspace);
        const char* prev_allele = orig_alleles_sorted[0];
        for (uint32_t aidx = 1; aidx != cur_allele_ct; ++aidx) {
          const char* cur_allele = orig_alleles_sorted[aidx];
          if (unlikely(strequal_overread(prev_allele, cur_allele))) {
            goto UpdateVarAlleles_ret_DUPLICATE_ALLELE_CODE_AFTER_UPDATE;
          }
          prev_allele = cur_allele;
        }
        // For multiallelic variants, also confirm missing allele code is not
        // present.
        if (unlikely(bsearch_strptr_overread(std_input_missing_geno, orig_alleles_sorted, cur_allele_ct) != -1)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --update-alleles file results in a missing allele code in a multiallelic variant.\n", line_idx);
          goto UpdateVarAlleles_ret_INCONSISTENT_INPUT_WW;
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
  UpdateVarAlleles_ret_TOO_MANY_ALLELES:
    logerrprintfww("Error: Line %" PRIuPTR " of --update-alleles file has too many alleles (this " PROG_NAME_STR " build is limited to %u).\n", line_idx, max_impl_input_allele_ct);
    reterr = kPglRetNotYetSupported;
    break;
  UpdateVarAlleles_ret_DUPLICATE_ALLELE_CODE_AFTER_UPDATE:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --update-alleles file results in a duplicated allele code.\n", line_idx);
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
    char** orig_alt_starts = nullptr; // spurious g++ 4.8 warning
    if (unlikely(bigstack_alloc_cp(raw_variant_ct, &variant_ids_copy) ||
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
          goto RecoverVarIds_ret_MALFORMED_INPUT;
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
          if (memequal_sk(linebuf_iter, "POS")) {
            cur_col_type = 0;
          } else if (memequal_sk(linebuf_iter, "REF")) {
            cur_col_type = 2;
          } else if (memequal_sk(linebuf_iter, "ALT")) {
            cur_col_type = 3;
          } else {
            continue;
          }
        } else if (token_slen == 2) {
          if (memequal_sk(linebuf_iter, "ID")) {
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
      --line_idx;
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
      if ((!(line_idx % 1000000)) && line_idx) {
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
        if (is_rigid || (cur_bp > prev_bp)) {
          const uint32_t arr_length = prev_chr_vidx_end - prev_vidx_end;
          const uint32_t incr = Expsearch0U32(&(variant_bps[prev_vidx_end]), arr_length, cur_bp);
          if (incr == arr_length) {
            prev_vidx_end = prev_chr_vidx_end;
            continue;
          }
          prev_vidx_start = incr + prev_vidx_end;
        }
        prev_bp = cur_bp;
      } else {
        const uint32_t arr_length = prev_chr_vidx_end - prev_chr_vidx_start;
        const uint32_t incr = LowerBoundNonemptyU32(&(variant_bps[prev_chr_vidx_start]), arr_length, cur_bp);
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
      strcpy_k(outname_end, ".recoverid.dup");
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
  RecoverVarIds_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  RecoverVarIds_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  RecoverVarIds_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  CleanupTextStream2("--recover-var-ids file", &txs, &reterr);
  fclose_cond(errfile);

  // doesn't hurt to ensure all replaced variant_ids[] entries are valid in the
  // error case...
  BigstackEndSet(alloc_end);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr Plink1ClusterImport(const char* within_fname, const char* catpheno_name, const char* family_missing_catname, const uintptr_t* sample_include, const char* sample_ids, const char* missing_catname, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t mwithin_val, uint32_t max_thread_ct, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr) {
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
      if (unlikely(bigstack_alloc_w(raw_sample_ctaw, &cat_nm) ||
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
    const uintptr_t data_vec_ct = Int32CtToVecCt(raw_sample_ct);
    const uint32_t missing_catname_slen = strlen(missing_catname);
    const uint32_t missing_catname_hval = Hashceil(missing_catname, missing_catname_slen, cat_htable_size);
    if (within_fname) {
      uintptr_t* already_seen;
      uint32_t* sorted_cat_idxs;
      char* idbuf;
      const char** cur_cat_names;
      if (unlikely(bigstack_calloc_w(BitCtToWordCt(sample_ct), &already_seen) ||
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
      if (unlikely(CopySortStrboxSubset(sample_include, sample_ids, sample_ct, max_sample_id_blen, 0, 0, &sorted_idbox, &id_map))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      reterr = SizeAndInitTextStream(within_fname, bigstack_left() - (bigstack_left() / 4), MAXV(max_thread_ct - 1, 1), &within_txs);
      if (unlikely(reterr)) {
        goto Plink1ClusterImport_ret_TSTREAM_FAIL;
      }
      char* cat_name_write_start = R_CAST(char*, g_bigstack_base);
      char* cat_name_iter = cat_name_write_start;
      char* cat_name_write_max = R_CAST(char*, BigstackEndRoundedDown());

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
        uint32_t lb_idx = bsearch_strbox_lb(idbuf, sorted_idbox, id_blen, max_sample_id_blen, sample_ct);
        *idbuf_iter = ' ';
        const uint32_t ub_idx = bsearch_strbox_lb(idbuf, sorted_idbox, id_blen, max_sample_id_blen, sample_ct);
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
        uint32_t cur_cat_idx = IdHtableAdd(main_token_start, cur_cat_names, main_token_slen, cat_htable_size, nonnull_cat_ct + 1, cat_htable);
        if (cur_cat_idx == UINT32_MAX) {
          if (unlikely(main_token_slen >= S_CAST(uintptr_t, cat_name_write_max - cat_name_iter))) {
            goto Plink1ClusterImport_ret_NOMEM;
          }
          char* cat_name_start = cat_name_iter;
          cat_name_iter = memcpya(cat_name_iter, main_token_start, main_token_slen + 1);
          // bugfix (5 May 2021): dropped this increment in mid-Jan refactor
          ++nonnull_cat_ct;
          cur_cat_idx = nonnull_cat_ct;
          cur_cat_names[cur_cat_idx] = cat_name_start;
        }
        // permit duplicates if category is identical
        if (IsSet(already_seen, lb_idx)) {
          const uint32_t existing_cat_idx = sorted_cat_idxs[lb_idx];
          if (unlikely(existing_cat_idx != cur_cat_idx)) {
            idbuf[fid_slen] = ' ';
            logpreprintfww("Error: Duplicate sample ID \"%s\" with conflicting category assignments in --within file.\n", idbuf);
            goto Plink1ClusterImport_ret_MALFORMED_INPUT_2;
          }
          ++duplicate_ct;
        } else {
          SetBit(lb_idx, already_seen);
          for (; lb_idx != ub_idx; ++lb_idx) {
            sorted_cat_idxs[lb_idx] = cur_cat_idx;
          }
        }
        goto Plink1ClusterImport_LINE_ITER_ALREADY_ADVANCED;
      }
      if (unlikely(!nonnull_cat_ct)) {
        if (line_idx == miss_ct + 1) {
          // could be fancy and only print the last part of this message if
          // all FIDs are 0
          logerrputs("Error: No sample IDs in --within file are in the main dataset.  Compare it with\nyour .fam/.psam file; note that if the latter has no FID column at all, that\ncorresponds to FIDs of '0', not FIDs equal to the IIDs.\n");
        } else {
          logerrputs("Error: All --within categories are null.\n");
        }
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
      BigstackBaseSet(cat_name_write_start);
      uint32_t* old_cat_idxs_to_new;
      if (unlikely(bigstack_alloc_u32(nonnull_cat_ct + 1, &old_cat_idxs_to_new))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
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
      if (unlikely(SortStrptrArrIndexed(nonnull_cat_ct + 1, 1, 0, 1, new_pheno_cols[old_pheno_ct].category_names, nullptr, old_cat_idxs_to_new))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      for (uint32_t sorted_sample_idx = 0; sorted_sample_idx != sample_ct; ++sorted_sample_idx) {
        const uint32_t cur_sample_uidx = id_map[sorted_sample_idx];
        uint32_t cur_cat_idx = sorted_cat_idxs[sorted_sample_idx];
        if (cur_cat_idx > sample_ct) {
          cur_cat_idx = 0;
        }
        if (cur_cat_idx) {
          SetBit(cur_sample_uidx, cat_nm);
        }
        cat_idxs[cur_sample_uidx] = old_cat_idxs_to_new[cur_cat_idx];
      }

      if (duplicate_ct) {
        logprintfww("Note: %" PRIuPTR " duplicate sample ID%s) in --within file.\n", duplicate_ct, (duplicate_ct == 1)? " (with a consistent category assignment" : "s (with consistent category assignments");
      }
      if (miss_ct) {
        snprintf(g_logbuf, kLogbufSize, "--within: %u non-null categor%s present, %" PRIuPTR " sample ID%s skipped.\n", nonnull_cat_ct, (nonnull_cat_ct == 1)? "y" : "ies", miss_ct, (miss_ct == 1)? "" : "s");
        WordWrapB(0);
      } else {
        snprintf(g_logbuf, kLogbufSize, "--within: %u non-null categor%s present.\n", nonnull_cat_ct, (nonnull_cat_ct == 1)? "y" : "ies");
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
      // guaranteed to have enough space, otherwise HtableGoodSizeAlloc would
      // have failed
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
      // add 'C' prefixes?
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
      BigstackFinalizeU32(cat_idx_m1_to_first_sample_uidx, nonnull_cat_ct);
      uint32_t* old_cat_idxs_to_new;
      if (unlikely(bigstack_alloc_u32(nonnull_cat_ct + 1, &old_cat_idxs_to_new))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      const uintptr_t catname_vec_ct = WordCtToVecCt(nonnull_cat_ct + 1);
      const uintptr_t catname_storage_vec_ct = DivUp(total_catname_blen, kBytesPerVec);
      if (unlikely(vecaligned_malloc((raw_sample_ctaw * kWordsPerVec + data_vec_ct + catname_vec_ct + catname_storage_vec_ct) * kBytesPerVec, &(new_pheno_cols[old_pheno_ct].nonmiss)))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      new_pheno_cols[old_pheno_ct].nonnull_category_ct = nonnull_cat_ct;
      uintptr_t* catdata_iter = new_pheno_cols[old_pheno_ct].nonmiss;
      memcpy(catdata_iter, cat_nm, raw_sample_ctaw * sizeof(intptr_t));
      catdata_iter = &(catdata_iter[raw_sample_ctaw]);

      uint32_t* cat_idx_dst = R_CAST(uint32_t*, catdata_iter);
      new_pheno_cols[old_pheno_ct].data.cat = cat_idx_dst;
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
      if (unlikely(SortStrptrArrIndexed(nonnull_cat_ct + 1, 1, 0, 1, new_pheno_cols[old_pheno_ct].category_names, nullptr, old_cat_idxs_to_new))) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      sample_uidx_base = 0;
      cur_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
        cat_idx_dst[sample_uidx] = old_cat_idxs_to_new[cat_idxs[sample_uidx]];
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

PglErr AlleleAlphanumUpdate(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, uint32_t variant_ct, AlleleAlphanumFlags flags, __attribute__((unused)) uint32_t max_thread_ct, char** allele_storage_mutable) {
  PglErr reterr = kPglRetSuccess;
  {
    // TODO: trivial to parallelize.
    const uint32_t is_acgt_to_1234 = (flags / kfAlleleAlphanum1234) & 1;
    const uint32_t allow_multichar = (flags / kfAlleleAlphanumMultichar) & 1;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t allele_ct = 2;
    if (is_acgt_to_1234) {
      // Single-character-allele conversion: operate directly on byte-offset
      // from g_one_char_strs.  Missing code permitted.
      const unsigned char acgtm_ptr_to_1234m[512] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '.'*2, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

        0, 0, '1'*2, 0, 0, 0, '2'*2, 0, 0, 0, 0, 0, 0, 0, '3'*2, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, '4'*2, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, '1'*2, 0, 0, 0, '2'*2, 0, 0, 0, 0, 0, 0, 0, '3'*2, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, '4'*2, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      // Multi-character conversion: ACGT->1234 conversion table if
      // allow_multichar, guaranteed-failure otherwise.  Missing allele code
      // doesn't come into play here.
      unsigned char acgt_to_1234[256];
      memset(acgt_to_1234, 0, 256);
      if (allow_multichar) {
        acgt_to_1234['A'] = '1';
        acgt_to_1234['C'] = '2';
        acgt_to_1234['G'] = '3';
        acgt_to_1234['T'] = '4';
        acgt_to_1234['a'] = '1';
        acgt_to_1234['c'] = '2';
        acgt_to_1234['g'] = '3';
        acgt_to_1234['t'] = '4';
      }
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        char** cur_alleles = &(allele_storage_mutable[allele_idx_offset_base]);
        for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
          char* cur_allele = cur_alleles[aidx];
          const uintptr_t one_char_offset = S_CAST(uintptr_t, cur_allele - K_CAST(char*, g_one_char_strs));
          if (one_char_offset < 512) {
            const uint32_t one_char_offset_conv = acgtm_ptr_to_1234m[one_char_offset];
            if (unlikely(!one_char_offset_conv)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --allele1234: Variant '%s' has an allele code that isn't in {A, C, G, T, a, c, g, t, <missing>}.\n", variant_ids[variant_uidx]);
              goto AlleleAlphanumUpdate_ret_INCONSISTENT_INPUT_WW;
            }
            cur_alleles[aidx] = &(K_CAST(char*, g_one_char_strs)[one_char_offset_conv]);
          } else {
            char* cur_allele_iter = cur_allele;
            unsigned char ucc = *cur_allele_iter;
            do {
              const unsigned char ucc_conv = acgt_to_1234[ucc];
              if (unlikely(!ucc_conv)) {
                if (allow_multichar) {
                  snprintf(g_logbuf, kLogbufSize, "Error: --allele1234: Variant '%s' has a multi-character allele code containing a non-ACGT character. (You can use --snps-only to exclude all variants with multi-character allele codes.)\n", variant_ids[variant_uidx]);
                } else {
                  snprintf(g_logbuf, kLogbufSize, "Error: --allele1234: Variant '%s' has a multi-character allele code. (Add the 'multichar' modifier if such allele codes are expected.)\n", variant_ids[variant_uidx]);
                }
                goto AlleleAlphanumUpdate_ret_INCONSISTENT_INPUT_WW;
              }
              *cur_allele_iter = ucc_conv;
              ucc = *(++cur_allele_iter);
            } while (ucc);
          }
        }
      }
    } else {
      const unsigned char acgtm_ptr_from_1234m[512] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '.'*2, 0, 0, 0,
        0, 0, 65*2, 0, 67*2, 0, 71*2, 0, 84*2, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
      unsigned char acgt_from_1234[256];
      memset(acgt_from_1234, 0, 256);
      if (allow_multichar) {
        acgt_from_1234['1'] = 'A';
        acgt_from_1234['2'] = 'C';
        acgt_from_1234['3'] = 'G';
        acgt_from_1234['4'] = 'T';
      }
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        char** cur_alleles = &(allele_storage_mutable[allele_idx_offset_base]);
        for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
          char* cur_allele = cur_alleles[aidx];
          const uintptr_t one_char_offset = S_CAST(uintptr_t, cur_allele - K_CAST(char*, g_one_char_strs));
          if (one_char_offset < 512) {
            const uint32_t one_char_offset_conv = acgtm_ptr_from_1234m[one_char_offset];
            if (unlikely(!one_char_offset_conv)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --alleleACGT: Variant '%s' has an allele code that isn't in {1, 2, 3, 4, <missing>}.\n", variant_ids[variant_uidx]);
              goto AlleleAlphanumUpdate_ret_INCONSISTENT_INPUT_WW;
            }
            cur_alleles[aidx] = &(K_CAST(char*, g_one_char_strs)[one_char_offset_conv]);
          } else {
            char* cur_allele_iter = cur_allele;
            unsigned char ucc = *cur_allele_iter;
            do {
              const unsigned char ucc_conv = acgt_from_1234[ucc];
              if (unlikely(!ucc_conv)) {
                if (allow_multichar) {
                  snprintf(g_logbuf, kLogbufSize, "Error: --alleleACGT: Variant '%s' has a multi-character allele code containing a non-1234 character. (You can use --snps-only to exclude all variants with multi-character allele codes.)\n", variant_ids[variant_uidx]);
                } else {
                  snprintf(g_logbuf, kLogbufSize, "Error: --alleleACGT: Variant '%s' has a multi-character allele code. (Add the 'multichar' modifier if such allele codes are expected.)\n", variant_ids[variant_uidx]);
                }
                goto AlleleAlphanumUpdate_ret_INCONSISTENT_INPUT_WW;
              }
              *cur_allele_iter = ucc_conv;
              ucc = *(++cur_allele_iter);
            } while (ucc);
          }
        }
      }
    }
    logprintf("--allele%s: %u variant%s updated.\n", is_acgt_to_1234? "1234" : "ACGT", variant_ct, (variant_ct == 1)? "" : "s");
  }
  while (0) {
  AlleleAlphanumUpdate_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
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
      } else {
        // bugfix (23 Sep 2020): forgot to add 2 for implicit FID=0
        cur_sample_id_blen += 2;
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
    // don't use LoadXidHeader since we can ignore sample ID on this pass
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
      line_iter = FirstSpaceOrEoln(token_start);
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
    uintptr_t max_paternal_id_blen = 0;
    uintptr_t max_maternal_id_blen = 0;
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
    if (new_fid_present) {
      // bugfix (14 Jan 2021)
      siip->flags |= kfSampleIdFidPresent;
    }
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    uint32_t* xid_map = nullptr;
    char* sorted_xidbox = nullptr;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto UpdateSampleIds_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* already_seen;
    char* idbuf;
    if (unlikely(bigstack_calloc_w(raw_sample_ctl, &already_seen) ||
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
    if (miss_ct) {
      snprintf(g_logbuf, kLogbufSize, "--update-ids: %u sample%s updated, %" PRIuPTR " ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      snprintf(g_logbuf, kLogbufSize, "--update-ids: %u sample%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logputsb();
    reterr = CheckXidUniqueness(sample_include, siip, " from --update-ids", sample_ct);
    if (unlikely(reterr)) {
      goto UpdateSampleIds_ret_1;
    }
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
    // permit very long lines since this can be pointed at .ped files
    reterr = SizeAndInitTextStream(fname, bigstack_left() - (bigstack_left() / 4), MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto UpdateSampleParents_ret_TSTREAM_FAIL;
    }

    const char* line_iter;
    XidMode xid_mode;
    {
      char* line_start;
      reterr = LoadXidHeader("update-parents", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeader0 : kfXidHeaderIgnoreSid, &line_idx, &txs, &xid_mode, &line_start);
      if (unlikely(reterr)) {
        if (reterr == kPglRetEof) {
          goto UpdateSampleParents_ret_REWIND_FAIL;
        }
        goto UpdateSampleParents_ret_TSTREAM_XID_FAIL;
      }
      line_iter = line_start;
    }
    uint32_t postid_pat_col_idx = 1;
    if (*line_iter == '#') {
      const uint32_t id_col_ct = GetXidColCt(xid_mode);
      const char* token_end = CurTokenEnd(NextTokenMult0(line_iter, id_col_ct - 1));
      while (1) {
        const char* token_start = FirstNonTspace(token_end);
        if (unlikely(IsEolnKns(*token_start))) {
          goto UpdateSampleParents_ret_REWIND_FAIL;
        }
        token_end = CurTokenEnd(token_start);
        if (strequal_k(token_start, "PAT", token_end - token_start)) {
          break;
        }
        ++postid_pat_col_idx;
      }
      line_iter = AdvPastDelim(line_iter, '\n');
      ++line_idx;
    } else {
      const uint32_t token_ct = CountTokens(line_iter);
      if (unlikely(token_ct < 3)) {
        goto UpdateSampleParents_ret_REWIND_FAIL;
      }
      xid_mode = (token_ct == 3)? kfXidModeIid : kfXidModeFidIid;
    }

    uint32_t* xid_map = nullptr;
    char* sorted_xidbox = nullptr;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto UpdateSampleParents_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* already_seen;
    char* idbuf;
    if (unlikely(bigstack_calloc_w(raw_sample_ctl, &already_seen) ||
                 bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto UpdateSampleParents_ret_NOMEM;
    }
    const uintptr_t max_paternal_id_blen = parental_id_infop->max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = parental_id_infop->max_maternal_id_blen;
    char* paternal_ids = parental_id_infop->paternal_ids;
    char* maternal_ids = parental_id_infop->maternal_ids;
    uint32_t hit_ct = 0;
    uintptr_t miss_ct = 0;
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
  UpdateSampleParents_ret_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&txs)) {
      break;
    }
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
      goto UpdateSampleSexes_ret_TSTREAM_FAIL;
    }

    // (Much of this boilerplate is shared with e.g. KeepColMatch(); it
    // probably belongs in its own function.)
    char* line_start;
    XidMode xid_mode;
    reterr = LoadXidHeader("update-sex", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeaderFixedWidth : kfXidHeaderFixedWidthIgnoreSid, &line_idx, &txs, &xid_mode, &line_start);
    if (reterr) {
      if (likely(reterr == kPglRetEof)) {
        reterr = kPglRetSuccess;
        logputs("--update-sex: 0 samples updated.\n");
        goto UpdateSampleSexes_ret_1;
      }
      goto UpdateSampleSexes_ret_TSTREAM_XID_FAIL;
    }
    const uint32_t id_col_ct = GetXidColCt(xid_mode);
    uint32_t col_num = update_sex_info_ptr->col_num;
    uint32_t postid_col_idx = 0;
    if ((*line_start == '#') && (!col_num)) {
      // search for 'SEX' column (any capitalization)
      const char* token_end = CurTokenEnd(NextTokenMult0(line_start, id_col_ct - 1));
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
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto UpdateSampleSexes_ret_1;
    }
    uintptr_t* already_seen;
    char* idbuf;
    if (unlikely(bigstack_calloc_w(raw_sample_ctl, &already_seen) ||
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
        } else if (unlikely((!male0) && (sexval != 30) && (ujj != 85))) {
          // allow 'N' = missing to make 1/2/NA work
          // allow 'U'/'u' since this is actually being used by Illumina
          // GenCall and Affymetrix APT
          // don't permit 'n' for now
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid sex value on line %" PRIuPTR " of --update-sex file. (Acceptable values: 1/M/m = male, 2/F/f = female, 0/N/U = missing.)\n", line_idx);
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
      if (unlikely(bigstack_alloc_u32(split_pheno_ct, &observed_cat_cts) ||
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
        const char* tok_end = Strchrnul(binstr_iter, ',');
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

PglErr WriteAlleleFreqs(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* nonref_flags, const uint64_t* founder_allele_ddosages, const double* imp_r2_vals, const char* ref_binstr, const char* alt1_binstr, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, PgenGlobalFlags gflags, FreqRptFlags freq_rpt_flags, uint32_t max_thread_ct, uint32_t nonfounders, char* outname, char* outname_end) {
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
      const uint32_t all_nonref = (gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
      const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags, freq_rpt_flags / kfAlleleFreqColMaybeprovref, raw_variant_ct, all_nonref);
      if (provref_col) {
        cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
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
      uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
        if (provref_col) {
          *cswritep++ = '\t';
          *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
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
                const char* cur_allele_end_or_eq = Strchrnul(cur_allele, '=');
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
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
            ref_histogram[LowerBoundNonemptyD(ref_freq_bounds, ref_boundary_ct, ref_allele_ddosage * tot_allele_ddosage_recip)] += 1;
          }
          if (alt1_histogram) {
            alt1_histogram[LowerBoundNonemptyD(alt1_freq_bounds, alt1_boundary_ct, alt1_allele_ddosage * tot_allele_ddosage_recip)] += 1;
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
            ref_histogram[LowerBoundNonemptyU64(ref_ddosage_bounds, ref_boundary_ct, cur_allele_ddosages[0])] += 1;
          }
          if (alt1_histogram) {
            alt1_histogram[LowerBoundNonemptyU64(alt1_ddosage_bounds, alt1_boundary_ct, cur_allele_ddosages[1])] += 1;
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

PglErr WriteGenoCounts(const uintptr_t* sample_include, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* nonref_flags, const STD_ARRAY_PTR_DECL(uint32_t, 3, raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, x_male_geno_cts), uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t male_ct, uint32_t nosex_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t x_start, uint32_t max_allele_slen, PgenGlobalFlags gflags, GenoCountsFlags geno_counts_flags, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
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
    const uintptr_t* sex_nonfemale = sex_male;
    uint32_t* sex_nonfemale_cumulative_popcounts = nullptr;  // chrY
    PgenVariant pgv;
    pgv.genovec = nullptr;
    pgv.patch_01_set = nullptr;
    pgv.patch_01_vals = nullptr;
    pgv.patch_10_set = nullptr;
    pgv.patch_10_vals = nullptr;
    uint32_t* diploid_pair_cts = nullptr;
    uint32_t* hap_cts = nullptr;
    const uint32_t more_counts_needed = (max_allele_ct > 2) && (geno_counts_flags & (kfGenoCountsColRefalt1 | kfGenoCountsColRefalt | kfGenoCountsColHomalt1 | kfGenoCountsColAltxy | kfGenoCountsColXy | kfGenoCountsColHapalt1 | kfGenoCountsColHapalt | kfGenoCountsColHap | kfGenoCountsColNumeq));
    const uint32_t nonfemale_ct = male_ct + nosex_ct;
    if (more_counts_needed) {
      if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                   bigstack_alloc_w(sample_ctl, &sex_male_collapsed) ||
                   bigstack_alloc_u32(raw_sample_ctl, &sex_nonfemale_cumulative_popcounts) ||
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
      if (nosex_ct) {
        uintptr_t* nonfemale_tmp;
        if (unlikely(bigstack_alloc_w(raw_sample_ctl, &nonfemale_tmp))) {
          goto WriteGenoCounts_ret_NOMEM;
        }
        BitvecXor3Copy(sample_include, sex_male, sex_nm, raw_sample_ct, nonfemale_tmp);
        sex_nonfemale = nonfemale_tmp;
      }
      FillCumulativePopcounts(sex_nonfemale, raw_sample_ctl, sex_nonfemale_cumulative_popcounts);
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
    const uint32_t all_nonref = (gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags, geno_counts_flags / kfGenoCountsColMaybeprovref, raw_variant_ct, all_nonref);
    if (provref_col) {
      cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
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
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uintptr_t* cur_sample_include = nullptr;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t is_autosomal_diploid = 0;
    uint32_t is_x = 0;
    uint32_t is_mt = 0;
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
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
        const uint32_t is_y = (chr_idx == y_code);
        is_mt = (chr_idx == mt_code);
        cur_sample_include = sample_include;
        const uint32_t* cur_cumulative_popcounts = sample_include_cumulative_popcounts;
        if (!is_autosomal_diploid) {
          if (is_y) {
            cur_sample_include = sex_nonfemale;
            cur_cumulative_popcounts = sex_nonfemale_cumulative_popcounts;
            nobs_base = nonfemale_ct;
          }
        } else if (hap_cts) {
          hap_cts[0] = 0;
          hap_cts[1] = 0;
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
      if (provref_col) {
        *cswritep++ = '\t';
        *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
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
            hapref_ct = cur_raw_geno_cts[0];
            hapalt1_ct = cur_raw_geno_cts[2];
            // treat hethap as missing in chrY or other pure-haploid case
            missing_ct = nobs_base - hapref_ct - hapalt1_ct;
            if (is_mt) {
              // behavior update (6 Aug 2025): mixed no longer treated as
              // missing
              het_ct = cur_raw_geno_cts[1];
              missing_ct -= het_ct;
            }
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
        reterr = PgrGetM(cur_sample_include, pssi, nobs_base, variant_uidx, simple_pgrp, &pgv);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, variant_uidx);
          goto WriteGenoCounts_ret_1;
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
          assert(hap_cts[0] == 0);
          assert(hap_cts[1] == 0);
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
          } else if (is_mt) {
            // behavior update (6 Aug 2025): mixed no longer treated as missing
            missing_ct = nobs_base - cur_raw_geno_cts[0] - cur_raw_geno_cts[1] - cur_raw_geno_cts[2];
            for (uintptr_t aidx = 0; aidx != allele_ct; ++aidx) {
              hap_cts[aidx] = diploid_pair_cts[aidx * allele_ctp1];
            }
          } else {
            // treat hethap as missing in chrY or other pure-haploid case
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
          } else if (is_mt) {
            for (uint32_t aidx_hi = xy_col_altonly; aidx_hi != allele_ct; ++aidx_hi) {
              for (uint32_t aidx_lo = xy_col_altonly; aidx_lo < aidx_hi; ++aidx_lo) {
                cswritep = u32toa_x(diploid_pair_cts[aidx_lo * allele_ct + aidx_hi], ',', cswritep);
              }
              cswritep = strcpya_k(cswritep, "0,");
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
          if (is_autosomal_diploid || is_x || is_mt) {
            for (uint32_t aidx_hi = is_mt; aidx_hi != allele_ct; ++aidx_hi) {
              for (uint32_t aidx_lo = 0; aidx_lo <= aidx_hi - is_mt; ++aidx_lo) {
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
  WriteGenoCounts_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 WriteGenoCounts_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr WriteMissingnessReports(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uint32_t* sample_missing_hc_cts, const uint32_t* sample_missing_dosage_cts, const uint32_t* sample_hethap_cts, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* nonref_flags, const uint32_t* variant_missing_hc_cts, const uint32_t* variant_missing_dosage_cts, const uint32_t* variant_hethap_cts, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uintptr_t max_allele_slen, uint32_t y_nosex_missing_stats, PgenGlobalFlags gflags, uint32_t first_hap_uidx, MissingRptFlags missing_rpt_flags, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uintptr_t* chry_missingstat_include = sex_male;
    if (y_nosex_missing_stats) {
      uintptr_t* nonfemale_tmp;
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &nonfemale_tmp))) {
        goto WriteMissingnessReports_ret_NOMEM;
      }
      BitvecXor3Copy(sample_include, sex_male, sex_nm, raw_sample_ctl, nonfemale_tmp);
      chry_missingstat_include = nonfemale_tmp;
    }

    const uint32_t chry_missingstat_sample_ct = PopcountWords(chry_missingstat_include, raw_sample_ctl);
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
      {
        char* write_iter = nobs_strs[0];
        *write_iter++ = '\t';
        write_iter = u32toa(variant_ct_nony, write_iter);
        nobs_slens[0] = write_iter - nobs_strs[0];
        variant_ct_recips[0] = 1.0 / u31tod(variant_ct_nony);
        write_iter = nobs_strs[1];
        *write_iter++ = '\t';
        write_iter = u32toa(variant_ct, write_iter);
        nobs_slens[1] = write_iter - nobs_strs[1];
      }
      variant_ct_recips[1] = 1.0 / u31tod(variant_ct);
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
        cswritep = AppendXid(sample_ids, sids, scol_fid, scol_sid, max_sample_id_blen, max_sid_blen, sample_uidx, cswritep);
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
        const uint32_t chry_included = IsSet(chry_missingstat_include, sample_uidx);
        if (scol_nobs) {
          cswritep = memcpya(cswritep, nobs_strs[chry_included], nobs_slens[chry_included]);
        }
        const double cur_variant_ct_recip = variant_ct_recips[chry_included];
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
      const uint32_t alt_col = missing_rpt_flags & kfMissingRptVcolAlt;
      if (alt_col) {
        cswritep = strcpya_k(cswritep, "\tALT");
      }
      const uint32_t all_nonref = (gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
      const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags, missing_rpt_flags / kfMissingRptVcolMaybeprovref, raw_variant_ct, all_nonref);
      if (provref_col) {
        cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
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
      uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
            const uint32_t cur_nobs = is_y? chry_missingstat_sample_ct : sample_ct;
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
        if (provref_col) {
          *cswritep++ = '\t';
          *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
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
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
    if (unlikely(bigstack_alloc_w(founder_ctl2, &(pgv.genovec)) ||
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
            PgenErrPrintNV(reterr, variant_uidx);
            goto GetMultiallelicMarginalCounts_ret_1;
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
      if (unlikely(bigstack_alloc_w(founder_x_ctaw, &founder_male_collapsed) ||
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
            PgenErrPrintNV(reterr, variant_uidx);
            goto GetMultiallelicMarginalCounts_ret_1;
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
  }
 GetMultiallelicMarginalCounts_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct ComputeHweXLnPvalsCtxStruct {
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

  double* hwe_x_ln_pvals;
} ComputeHweXLnPvalsCtx;

void ComputeHweXLnPvalsMain(uintptr_t tidx, uintptr_t thread_ct, ComputeHweXLnPvalsCtx* ctx) {
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
  double* hwe_x_ln_pvals_iter = &(ctx->hwe_x_ln_pvals[variant_idx + xgeno_idx]);
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, ctx->variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
  uint32_t pct = 0;
  uint32_t next_print_variant_idx = variant_idx_end;
  if (!tidx) {
    next_print_variant_idx = (variant_idx_end + 99) / 100;
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
    *hwe_x_ln_pvals_iter++ = HweXchrLnP(female_1copy_ct, female_2copy_ct, female_0copy_ct, male_1copy_ct, male_0copy_ct, hwe_midp);
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
          *hwe_x_ln_pvals_iter++ = HweXchrLnP(female_1copy_ct, female_2copy_ct, female_0copy_ct, male_1copy_ct, male_0copy_ct, hwe_midp);
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
      next_print_variant_idx = (pct * S_CAST(uint64_t, variant_idx_end) + 99) / 100;
    }
  }
  if (pct > 10) {
    putc_unlocked('\b', stdout);
  }
}

THREAD_FUNC_DECL ComputeHweXLnPvalsThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  ComputeHweXLnPvalsCtx* ctx = S_CAST(ComputeHweXLnPvalsCtx*, arg->sharedp->context);
  ComputeHweXLnPvalsMain(arg->tidx, GetThreadCt(arg->sharedp) + 1, ctx);
  THREAD_RETURN;
}

PglErr ComputeHweXLnPvals(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts), uint32_t x_start, uint32_t hwe_x_ct, uintptr_t x_xallele_ct, uint32_t hwe_midp, uint32_t calc_thread_ct, double** hwe_x_ln_pvals_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  ComputeHweXLnPvalsCtx ctx;
  {
    assert(hwe_x_ct);
    if (unlikely(bigstack_alloc_d(hwe_x_ct + x_xallele_ct, hwe_x_ln_pvals_ptr))) {
      goto ComputeHweXLnPvals_ret_NOMEM;
    }
    bigstack_mark = g_bigstack_base;
    ctx.hwe_x_ln_pvals = *hwe_x_ln_pvals_ptr;

    if (calc_thread_ct > hwe_x_ct) {
      calc_thread_ct = hwe_x_ct;
    }
    if (unlikely(SetThreadCt0(calc_thread_ct - 1, &tg) ||
                 bigstack_alloc_u32(calc_thread_ct, &ctx.variant_uidx_starts) ||
                 bigstack_alloc_w(calc_thread_ct, &ctx.extra_aidx_starts))) {
      goto ComputeHweXLnPvals_ret_NOMEM;
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
      SetThreadFuncAndData(ComputeHweXLnPvalsThread, &ctx, &tg);
      DeclareLastThreadBlock(&tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto ComputeHweXLnPvals_ret_THREAD_CREATE_FAIL;
      }
    }
    ComputeHweXLnPvalsMain(calc_thread_ct - 1, calc_thread_ct, &ctx);
    JoinThreads0(&tg);
    fputs("\b\b", stdout);
    logputs("done.\n");
  }
  while (0) {
  ComputeHweXLnPvals_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ComputeHweXLnPvals_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr HardyReport(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* nonref_flags, const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, autosomal_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_nosex_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts), const double* hwe_x_ln_pvals, uint32_t variant_ct, uint32_t hwe_x_ct, uint32_t max_allele_slen, PgenGlobalFlags gflags, double output_min_ln, HardyFlags hardy_flags, uint32_t max_thread_ct, uint32_t nonfounders, char* outname, char* outname_end) {
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
    uintptr_t overflow_buf_alloc = overflow_buf_size;
    if (output_zst) {
      overflow_buf_alloc += CstreamWkspaceReq(overflow_buf_size);
    }
    char* overflow_buf;
    uintptr_t* chr_skips;
    if (unlikely(bigstack_alloc_c(overflow_buf_alloc, &overflow_buf) ||
                 bigstack_alloc_w(chr_code_endl, &chr_skips))) {
      goto HardyReport_ret_NOMEM;
    }
    // skip chrX, chrY, chrM here
    memcpy(chr_skips, cip->haploid_mask, chr_code_endl * sizeof(intptr_t));
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
    const uint32_t all_nonref = (gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    // customize this, since the standard calculation does not exclude haploid
    // chromosomes
    uint32_t provref_col = 0;
    if (ref_col) {
      if (hardy_flags & kfHardyColProvref) {
        provref_col = 1;
      } else if (!(hardy_flags & kfHardyColMaybeprovref)) {
        // do nothing
      } else if (!nonref_flags) {
        provref_col = all_nonref;
      } else {
        const uint32_t chr_ct = cip->chr_ct;
        for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
          const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
          if ((!IsSet(cip->chr_mask, chr_idx)) || IsSet(chr_skips, chr_idx)) {
            continue;
          }
          if (!IntersectionRangeIsEmpty(variant_include, nonref_flags, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1])) {
            provref_col = 1;
            break;
          }
        }
      }
    }
    const uint32_t report_neglog10p = (hardy_flags / kfHardyLog10) & 1;
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
      if (provref_col) {
        cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
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
        if (report_neglog10p) {
          cswritep = strcpya_k(cswritep, "NEG_LOG10_");
        }
        if (midp) {
          cswritep = strcpya_k(cswritep, "MID");
        }
        *cswritep++ = 'P';
      }
      AppendBinaryEoln(&cswritep);
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
          if (provref_col) {
            *cswritep++ = '\t';
            *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
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
            const double hwe_ln_p = HweLnP(het_a1_ct, hom_a1_ct, two_ax_ct, midp);
            if (report_neglog10p) {
              const double reported_val = (-kRecipLn10) * hwe_ln_p;
              cswritep = dtoa_g(reported_val, cswritep);
            } else {
              cswritep = lntoa_g(MAXV(hwe_ln_p, output_min_ln), cswritep);
            }
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
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
      const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
      const uint32_t x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
      if (ref_col && ((hardy_flags & (kfHardyColMaybeprovref | kfHardyColProvref)) == kfHardyColMaybeprovref) && nonref_flags) {
        provref_col = !IntersectionRangeIsEmpty(variant_include, nonref_flags, x_start, cip->chr_fo_vidx_start[x_chr_fo_idx + 1]);
      }
      if (provref_col) {
        cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
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
        if (report_neglog10p) {
          cswritep = strcpya_k(cswritep, "NEG_LOG10_");
        }
        if (midp) {
          cswritep = strcpya_k(cswritep, "MID");
        }
        *cswritep++ = 'P';
      }
      if (p_col) {
        *cswritep++ = '\t';
        if (report_neglog10p) {
          cswritep = strcpya_k(cswritep, "NEG_LOG10_");
        }
        if (midp) {
          cswritep = strcpya_k(cswritep, "MID");
        }
        *cswritep++ = 'P';
      }
      AppendBinaryEoln(&cswritep);
      fputs("--hardy: Writing chrX results...", stdout);
      fflush(stdout);
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
          if (provref_col) {
            *cswritep++ = '\t';
            *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
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
            // Correct to increment this before looking up hwe_x_ln_pvals[]
            // (and to not increment on a1_idx == 0).
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
            const double female_hwe_ln_p = HweLnP(female_het_a1_ct, female_hom_a1_ct, female_two_ax_ct, midp);
            if (report_neglog10p) {
              const double reported_val = (-kRecipLn10) * female_hwe_ln_p;
              cswritep = dtoa_g(reported_val, cswritep);
            } else {
              cswritep = lntoa_g(MAXV(female_hwe_ln_p, output_min_ln), cswritep);
            }
          }
          if (p_col) {
            *cswritep++ = '\t';
            // bugfix (27 Jun 2020): forgot to correct this for multiallelic
            // variants
            const double ln_pval = hwe_x_ln_pvals[variant_idx + xgeno_idx];
            if (report_neglog10p) {
              const double reported_val = (-kRecipLn10) * ln_pval;
              cswritep = dtoa_g(reported_val, cswritep);
            } else {
              cswritep = lntoa_g(MAXV(ln_pval, output_min_ln), cswritep);
            }
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
  uint32_t max_difflist_len;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** raregenos;
  uint32_t** difflist_sample_id_bufs;
  uint32_t* read_variant_uidx_starts;
  uint32_t cur_block_size;

  uint64_t err_info;

  // top-level: length calc_thread_ct array
  // second level: length-20 arrays, corresponding to the 20 chr_type x
  //               subst_code possibilities
  // bottom level: acc2 partial counts, then acc4, then acc8, then acc32 for
  //               genotype=0.  Then the same for genotype=2, and finally (in
  //               autosomal-diploid/chrX cases) genotype=1.
  VecW*** thread_dense_counts;

  // top-level: length calc_thread_ct array
  // second level: length-40 arrays, corresponding to chr_type x common_geno x
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
  //   five correspond to ALT1 subst_codes, and the last five correspond to the
  //   actually-observed-ALT subst_codes.  The sum of the last five arrays is
  //   always twice that of the first five arrays, since each allele is counted
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
  const uint32_t max_difflist_len = ctx->max_difflist_len;
  const uint32_t skip_y = (!male_ct) && (!ctx->y_nonmale_needed);
  uintptr_t* genovec = ctx->genovecs[tidx];
  PgenVariant pgv;
  pgv.genovec = genovec;
  SetPgvThreadMhcNull(sample_ct, tidx, ctx->thread_read_mhc, &pgv);
  uintptr_t* raregeno = ctx->raregenos[tidx];
  uint32_t* difflist_sample_ids = ctx->difflist_sample_id_bufs[tidx];

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
  uint64_t new_err_info = 0;
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
        const PglErr reterr = PgrGetDifflistOrGenovec(sample_include, pssi, sample_ct, max_difflist_len, variant_uidx, pgrp, genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto SampleCountsThread_err;
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
      const PglErr reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, pgrp, &pgv);
      if (unlikely(reterr)) {
        new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
        goto SampleCountsThread_err;
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
      } else { // subst_codes_vary
        if (chr_type < 2) {
          const uint32_t patch_01_ct = pgv.patch_01_ct;
          if (patch_01_ct) {
            const uintptr_t* patch_01_set = pgv.patch_01_set;
            const AlleleCode* patch_01_vals = pgv.patch_01_vals;
            uintptr_t sample_idx_base = 0;
            uintptr_t sample_idx_bits = patch_01_set[0];
            if (!chr_type) {
              for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &sample_idx_bits);
                const uint32_t cur_subst_code = alt_subst_codes[patch_01_vals[uii]];
                if (subst_code1 != cur_subst_code) {
                  het_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
                  het_rarealt_cts[cur_subst_code * sample_ct_i32av + sample_idx] += 1;
                }
              }
            } else {
              for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &sample_idx_bits);
                // bugfix (9 Jan 2023): in male case, we will adjust hethap
                // to missing later under the assumption that this genotype
                // is ref/alt1.  So, subst_code doesn't matter.
                const uint32_t is_male = IsSet(sex_male_collapsed, sample_idx);
                if (!is_male) {
                  const uint32_t cur_subst_code = alt_subst_codes[patch_01_vals[uii]];
                  if (subst_code1 != cur_subst_code) {
                    het_rarealt_cts[subst_code1_offset + sample_idx] -= 1;
                    het_rarealt_cts[cur_subst_code * sample_ct_i32av + sample_idx] += 1;
                  }
                }
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
  {
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
  }
  while (0) {
  SampleCountsThread_err:
    UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    THREAD_BLOCK_FINISH(arg);
    break;
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
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &ctx.sample_include_cumulative_popcounts))) {
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
      // bugfix (28 Mar 2021): if chrX is present, no autosomal diploid
      // chromosomes are present, and there are no males, we need to set the
      // low bit of chr_types and increment diploid_chr_type_ct.
      uint32_t x_code;
      if ((!haploid_mask_lowbit) && XymtExists(cip, kChrOffsetX, &x_code)) {
        if (male_ct) {
          chr_types |= 2;
          ++diploid_chr_type_ct;
        } else if (!chr_types) {
          chr_types = 1;
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
      if (unlikely(bigstack_alloc_w(sample_ctl, &ctx.sex_male_collapsed))) {
        goto SampleCounts_ret_NOMEM;
      }
      CopyBitarrSubset(sex_male, sample_include, sample_ct, ctx.sex_male_collapsed);
    } else {
      ctx.sex_male_collapsed = nullptr;
    }

    // don't subtract 1 after load-balancing is improved?
    uint32_t calc_thread_ct = (max_thread_ct > 4)? (max_thread_ct - 1) : max_thread_ct;
    if (unlikely(bigstack_alloc_wp(calc_thread_ct, &ctx.raregenos) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.difflist_sample_id_bufs) ||
                 bigstack_alloc_vpp(calc_thread_ct, &ctx.thread_dense_counts) ||
                 bigstack_alloc_u32pp(calc_thread_ct, &ctx.thread_sparse_counts) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_sparse_common0_cts))) {
      goto SampleCounts_ret_NOMEM;
    }
    // todo: tune this threshold
    const uint32_t max_difflist_len = sample_ct / 32;
    ctx.max_difflist_len = max_difflist_len;

    const uintptr_t raregeno_vec_ct = DivUp(max_difflist_len, kNypsPerVec);
    const uintptr_t difflist_sample_id_vec_ct = DivUp(max_difflist_len, kInt32PerVec);

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
      if (unlikely(bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_diploid_singleton_cts))) {
        goto SampleCounts_ret_NOMEM;
      }
      thread_xalloc_vec_ct += sample_ct_i32v;
    }
    ctx.thread_singleton_cts = nullptr;
    if (flags & kfSampleCountsColSingle) {
      if (unlikely(bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_singleton_cts))) {
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
      if (unlikely(bigstack_alloc_u16p(calc_thread_ct, &ctx.thread_alt_subst_codes) ||
                   bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_het2alt_cts) ||
                   bigstack_alloc_i32p(calc_thread_ct, &ctx.thread_het_rarealt_cts) ||
                   bigstack_alloc_i32p(calc_thread_ct, &ctx.thread_hom_rarealt_cts))) {
        goto SampleCounts_ret_NOMEM;
      }
      alt_subst_codes_vec_ct = DivUp(kPglMaxAlleleCt, kInt16PerVec);
      het2alt_vec_ct = sample_ct_i32v * 2 * kSubstCodeCt;
      het_rarealt_vec_ct = sample_ct_i32v * kSubstCodeCt;
      hom_rarealt_vec_ct = sample_ct_i32v * kSubstCodeCt;
      thread_xalloc_vec_ct += alt_subst_codes_vec_ct + het2alt_vec_ct + het_rarealt_vec_ct + hom_rarealt_vec_ct;
      if (chr_types & 14) {
        if (unlikely(bigstack_alloc_i32p(calc_thread_ct, &ctx.thread_hap_rarealt_cts))) {
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
    ctx.err_info = (~0LLU) << 32;
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
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_size = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto SampleCounts_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, ctx.err_info);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, ctx.err_info >> 32);
          goto SampleCounts_ret_1;
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
    if (unlikely(bigstack_alloc_u32(acc32_vec_ct * kInt32PerVec, &unscramble_buf))) {
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
          if (unlikely(bigstack_end_calloc_u32(sample_ct_i32av, &male_u32_mask))) {
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
              if (unlikely(bigstack_end_calloc_u32(sample_ct_i32av, &dst))) {
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
            if (unlikely(bigstack_end_calloc_u32(sample_ct_i32av, &dst))) {
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
      if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &het2alt_reported_vals))) {
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
    if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &refalt_vals) ||
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
      if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &het_vals))) {
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
    if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &hapref_vals) ||
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
      if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &dst))) {
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
        if (unlikely(bigstack_calloc_u32(sample_ct_i32av, &(final_counts[kSampleCountHet2alt])))) {
          goto SampleCounts_ret_NOMEM;
        }
      }
    }

    if (flags & kfSampleCountsColHetSnp) {
      uint32_t* dst;
      if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &dst))) {
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
      if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &dst))) {
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
      if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &dst))) {
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
      if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &dst))) {
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
      if (unlikely(bigstack_alloc_u32(sample_ct_i32av, &dst))) {
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
    if (unlikely(bigstack_alloc_c(overflow_buf_alloc, &overflow_buf))) {
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
      cswritep = AppendXid(sample_ids, sids, col_fid, col_sid, max_sample_id_blen, max_sid_blen, sample_uidx, cswritep);
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
  pgfip->block_base = nullptr;
  return reterr;
}

static const Dosage kGenoToDosage[4] = {0, kDosageMid, kDosageMax, kDosageMissing};

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
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
        PgenErrPrintNV(reterr, variant_uidx);
        goto SdiffCountsOnly_ret_1;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
              // don't need pair_sex_nonfemale here since we prohibit
              // unknown-sex pairs when chrY is present
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
    PgenErrPrintN(reterr);
    break;
  }
 SdiffCountsOnly_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

void AppendSdiffHeaderLine(SdiffFlags flags, uint32_t col_provref, uint32_t col_fid, uint32_t col_sid, uint32_t dosage_reported, char** cswritepp) {
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
  if (col_provref) {
    cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
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
    if (!dosage_reported) {
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

BoolErr AppendSdiffPregenoFields(const SdiffWriteContext* swcp, const char* const* cur_alleles, uint32_t sample_idx1, uint32_t sample_idx2, uint32_t variant_uidx, uint32_t allele_ct, uint32_t col_provref, uint32_t provref, CompressStreamState* cssp, char** cswritepp) {
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
  if (col_provref) {
    *cswritep++ = '\t';
    *cswritep++ = provref? 'Y' : 'N';
  }
  if (swcp->collapsed_sample_fmtids) {
    *cswritep++ = '\t';
    cswritep = strcpyax(cswritep, &(swcp->collapsed_sample_fmtids[sample_idx1 * swcp->max_sample_fmtid_blen]), '\t');
    cswritep = strcpya(cswritep, &(swcp->collapsed_sample_fmtids[sample_idx2 * swcp->max_sample_fmtid_blen]));
  }
  *cswritepp = cswritep;
  return 0;
}

PglErr SdiffMainBatch(const uintptr_t* __restrict sample_include, const uint32_t* __restrict sample_include_cumulative_popcounts,  const SampleIdInfo* siip, const uint32_t* __restrict id_pairs, const uintptr_t* __restrict pair_sex_male, const uintptr_t* __restrict variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* __restrict allele_idx_offsets, const char* const* allele_storage, const SdiffInfo* sdip, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uintptr_t id_pair_ct, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end, SdiffCounts* sdiff_counts) {
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
    if (unlikely(bigstack_calloc_cp(file_ct, &cswritep_arr) ||
                 BIGSTACK_ALLOC_X(CompressStreamState, file_ct, &css_arr))) {
      goto SdiffMainBatch_ret_NOMEM;
    }
    swc.flags = flags;
    swc.chr_buf_blen = 0;
    for (uintptr_t fidx = 0; fidx != id_pair_ct; ++fidx) {
      PreinitCstream(&(css_arr[fidx]));
    }
    const uintptr_t* nonref_flags = PgrGetNonrefFlags(simple_pgrp);
    const uint32_t all_nonref = (PgrGetGflags(simple_pgrp) & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t col_provref = (flags & kfSdiffColRef) && ProvrefCol(variant_include, nonref_flags, flags / kfSdiffColMaybeprovref, raw_variant_ct, all_nonref);
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
        AppendSdiffHeaderLine(flags, col_provref, col_fid, col_sid, dosage_reported, &(cswritep_arr[fidx]));
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
      AppendSdiffHeaderLine(flags, col_provref, col_fid, col_sid, dosage_reported, &(cswritep_arr[0]));
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
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
      }
      if (unlikely(reterr)) {
        PgenErrPrintNV(reterr, variant_uidx);
        goto SdiffMainBatch_ret_1;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
      }
      const uint32_t provref = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)));
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
            if (unlikely(AppendSdiffPregenoFields(&swc, cur_alleles, sample_idx1, sample_idx2, variant_uidx, allele_ct, col_provref, provref, &(css_arr[fidx]), &cswritep))) {
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
            if (unlikely(AppendSdiffPregenoFields(&swc, cur_alleles, sample_idx1, sample_idx2, variant_uidx, allele_ct, col_provref, provref, &(css_arr[fidx]), &cswritep))) {
              cswritep_arr[fidx] = cswritep;
              goto SdiffMainBatch_ret_WRITE_FAIL;
            }
            if (flags & kfSdiffColGeno) {
              *cswritep++ = '\t';
              if (dosage1 == kDosageMissing) {
                *cswritep++ = '.';
              } else if (is_diploid_pair) {
                cswritep = PrintSmallDosage(dosage1, cswritep);
              } else {
                cswritep = PrintHaploidDosage(dosage1, cswritep);
              }
              *cswritep++ = '\t';
              if (dosage2 == kDosageMissing) {
                *cswritep++ = '.';
              } else if (is_diploid_pair) {
                cswritep = PrintSmallDosage(dosage2, cswritep);
              } else {
                cswritep = PrintHaploidDosage(dosage2, cswritep);
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
            if (unlikely(AppendSdiffPregenoFields(&swc, cur_alleles, sample_idx1, sample_idx2, variant_uidx, allele_ct, col_provref, provref, &(css_arr[fidx]), &cswritep))) {
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

PglErr Sdiff(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const SdiffInfo* sdip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t iid_sid, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  TextStream txs; // file=
  PreinitTextStream(&txs);
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    // Determine xid_mode.
    const SdiffFlags flags = sdip->flags;
    const uint32_t other_id_ct = sdip->other_id_ct;
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
        } else {
          if (unlikely(strchr(&(second_delim_ptr[1]), '\t'))) {
            logerrputs("Error: Too many instances of id-delim= character in --sample-diff sample ID.\n");
            goto Sdiff_ret_INVALID_CMDLINE;
          }
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
    if ((xid_mode & kfXidModeFlagSid) && (!siip->sids) && (!(siip->flags & kfSampleIdStrictSid0))) {
      xid_mode ^= kfXidModeFlagSid | kfXidModeFlagSkipSid;
    }
    char* sorted_xidbox;
    uint32_t* xid_map;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(orig_sample_include, siip, orig_sample_ct, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
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
      if (unlikely(bigstack_calloc_w(raw_sample_ctl, &sample_include) ||
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
          if (!idbuf[0]) {
            logerrprintfww("Error: --sample-diff sample ID (on line %" PRIuPTR " of file) not found.\n", line_idx);
          } else {
            TabsToSpaces(idbuf);
            logerrprintfww("Error: --sample-diff sample ID '%s' (on line %" PRIuPTR " of file) not found.\n", idbuf, line_idx);
          }
          goto Sdiff_ret_INCONSISTENT_INPUT;
        }
        linebuf_iter = FirstNonTspace(linebuf_iter);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto Sdiff_ret_MISSING_TOKENS;
        }
        uint32_t sample_uidx2;
        if (unlikely(SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &linebuf_iter, &sample_uidx2, idbuf))) {
          if (!idbuf[0]) {
            logerrprintfww("Error: --sample-diff sample ID (on line %" PRIuPTR " of file) not found.\n", line_idx);
          } else {
            TabsToSpaces(idbuf);
            logerrprintfww("Error: --sample-diff sample ID '%s' (on line %" PRIuPTR " of file) not found.\n", idbuf, line_idx);
          }
          goto Sdiff_ret_INCONSISTENT_INPUT;
        }
        if (unlikely(sample_uidx1 == sample_uidx2)) {
          TabsToSpaces(idbuf);
          logerrprintfww("Error: Duplicate sample ID \"%s\" on line %" PRIuPTR " of --sample-diff file.\n", idbuf, line_idx);
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
      if (unlikely(bigstack_calloc_w(raw_sample_ctl, &sample_include) ||
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
      if (unlikely(bigstack_alloc_w(id_pair_ctl, &pair_sex_male) ||
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
        reterr = SdiffMainBatch(sample_include, sample_include_cumulative_popcounts, siip, &(id_pairs[2 * id_pair_idx_start]), &(pair_sex_male[id_pair_idx_start / kBitsPerWord]), variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, sdip, sample_ct, raw_variant_ct, variant_ct, id_pair_idx_end - id_pair_idx_start, max_thread_ct, simple_pgrp, outname, outname_end, &(sdiff_counts[id_pair_idx_start]));
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
  Sdiff_ret_MISSING_TOKENS:
    logerrprintf("Error: Line %" PRIuPTR " of --sample-diff file has fewer tokens than expected.\n", line_idx);
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
PglErr WriteCovar(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, const uint32_t* new_sample_idx_to_old, const char* output_missing_pheno, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, WriteCovarFlags write_covar_flags, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    snprintf(outname_end, kMaxOutfnameExtBlen, ".cov");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto WriteCovar_ret_OPEN_FAIL;
    }
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
        if (unlikely(StrboxHtableAdd("SEX", covar_names, max_covar_name_blen, strlen("SEX"), covar_name_htable_size, covar_ct, covar_name_htable) != UINT32_MAX)) {
          logerrputs("Error: .cov file cannot have both a regular SEX column and a covariate named\n'SEX'.  Exclude or rename one of these columns.\n");
          goto WriteCovar_ret_INCONSISTENT_INPUT;
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
            // can't just use StrboxHtableFind() since "SEX" may not be stored
            // in covar_names[]
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
      write_iter = AppendXid(sample_ids, sids, write_fid, write_sid, max_sample_id_blen, max_sid_blen, sample_uidx, write_iter);
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
  const uintptr_t* variant_subset;

  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  const uintptr_t* founder_info_collapsed;  // only non-null if small-sample
  const uintptr_t* founder_info_interleaved_vec;
  const double* allele_freqs;  // nullptr if small-sample
  uint32_t sample_ct;
  uint32_t founder_ct;
  uint32_t max_difflist_len;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** raregenos;
  uint32_t** difflist_sample_id_bufs;
  uint32_t** allele_nobs_bufs;
  VecW** scrambled_ohet_bufs;
  uint32_t* read_variant_uidx_starts;
  uint32_t cur_block_size;

  uint64_t err_info;

  double* thread_ehet_base;
  uint32_t* thread_nobs_base;
  uint32_t* thread_monomorphic_ct;
  uint32_t** thread_ohets;
  double** thread_ehet_incrs;
  int32_t** thread_nobs_incrs;
} HetCtx;

THREAD_FUNC_DECL HetThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  HetCtx* ctx = S_CAST(HetCtx*, arg->sharedp->context);

  const uintptr_t* variant_subset = ctx->variant_subset;
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
  const uint32_t max_difflist_len = ctx->max_difflist_len;

  double ehet_base = 0.0;
  uint32_t nobs_base = 0;
  uint32_t monomorphic_ct = 0;
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
  uint64_t new_err_info = 0;
  do {
    const uint32_t cur_block_size = ctx->cur_block_size;
    const uint32_t cur_idx_ct = (((tidx + 1) * cur_block_size) / calc_thread_ct) - ((tidx * cur_block_size) / calc_thread_ct);
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_subset, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
    for (uint32_t cur_idx = 0; cur_idx != cur_idx_ct; ++cur_idx) {
      const uint32_t variant_uidx = BitIter1(variant_subset, &variant_uidx_base, &cur_bits);
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
          if (ehet < k2m35) {
            ++monomorphic_ct;
            continue;
          }
        }
        uint32_t difflist_common_geno;
        uint32_t difflist_len;
        const PglErr reterr = PgrGetDifflistOrGenovec(sample_include, pssi, sample_ct, max_difflist_len, variant_uidx, pgrp, genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto HetThread_err;
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
              ++monomorphic_ct;
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
              ++monomorphic_ct;
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
            ++monomorphic_ct;
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
          if (ehet < k2m35) {
            ++monomorphic_ct;
            continue;
          }
        }
        const PglErr reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, pgrp, &pgv);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto HetThread_err;
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
            ++monomorphic_ct;
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
  {
    ctx->thread_ehet_base[tidx] = ehet_base;
    ctx->thread_nobs_base[tidx] = nobs_base;
    ctx->thread_monomorphic_ct[tidx] = monomorphic_ct;

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
  }
  while (0) {
  HetThread_err:
    UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    THREAD_BLOCK_FINISH(arg);
    break;
  }
  THREAD_RETURN;
}

PglErr HetCalcMain(const uintptr_t* sample_include, const uintptr_t* variant_subset, const uintptr_t* allele_idx_offsets, const double* allele_freqs, const uintptr_t* founder_info, const char* calcstr, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t small_sample, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uint32_t** ohets_ptr, double** ehet_incrs_ptr, int32_t** nobs_incrs_ptr, double* ehet_base_ptr, int32_t* nobs_base_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  HetCtx ctx;
  {
    // return values
    if (unlikely(bigstack_alloc_u32(sample_ct, ohets_ptr) ||
                 bigstack_alloc_d(sample_ct, ehet_incrs_ptr) ||
                 bigstack_alloc_i32(sample_ct, nobs_incrs_ptr))) {
      goto HetCalcMain_ret_NOMEM;
    }
    bigstack_mark = g_bigstack_base;

    ctx.variant_subset = variant_subset;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.sample_include = sample_include;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uint32_t* sample_include_cumulative_popcounts;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts))) {
      goto HetCalcMain_ret_NOMEM;
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
      if (unlikely(bigstack_alloc_w(sample_ctaw, &founder_info_collapsed) ||
                   bigstack_alloc_w(sample_ctaw, &founder_info_interleaved_vec))) {
        goto HetCalcMain_ret_NOMEM;
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
    if (unlikely(bigstack_alloc_wp(calc_thread_ct, &ctx.raregenos) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.difflist_sample_id_bufs) ||
                 bigstack_alloc_vp(calc_thread_ct, &ctx.scrambled_ohet_bufs) ||
                 bigstack_alloc_d(calc_thread_ct, &ctx.thread_ehet_base) ||
                 bigstack_alloc_u32(calc_thread_ct, &ctx.thread_nobs_base) ||
                 bigstack_alloc_u32(calc_thread_ct, &ctx.thread_monomorphic_ct) ||
                 bigstack_alloc_u32p(calc_thread_ct, &ctx.thread_ohets) ||
                 bigstack_alloc_dp(calc_thread_ct, &ctx.thread_ehet_incrs) ||
                 bigstack_alloc_i32p(calc_thread_ct, &ctx.thread_nobs_incrs))) {
      goto HetCalcMain_ret_NOMEM;
    }
    // todo: tune this threshold
    const uint32_t max_difflist_len = sample_ct / 32;
    ctx.max_difflist_len = max_difflist_len;
    const uintptr_t raregeno_vec_ct = DivUp(max_difflist_len, kNypsPerVec);
    const uintptr_t difflist_sample_id_vec_ct = DivUp(max_difflist_len, kInt32PerVec);
    const uint32_t mhc_needed = (max_allele_ct > 2);
    uintptr_t allele_nobs_vec_ct = 0;
    ctx.allele_nobs_bufs = nullptr;
    if (mhc_needed) {
      if (unlikely(bigstack_alloc_u32p(calc_thread_ct, &ctx.allele_nobs_bufs))) {
        goto HetCalcMain_ret_NOMEM;
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
    if (unlikely(PgenMtLoadInit(ctx.variant_subset, raw_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &ctx.genovecs, mhc_needed? (&ctx.thread_read_mhc) : nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto HetCalcMain_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto HetCalcMain_ret_NOMEM;
    }
    ctx.err_info = (~0LLU) << 32;
    assert(bigstack_left() >= thread_xalloc_cacheline_ct * kCacheline * calc_thread_ct);
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
      if (!tidx) {
        ctx.thread_ohets[0] = *ohets_ptr;
        ctx.thread_ehet_incrs[0] = *ehet_incrs_ptr;
        ctx.thread_nobs_incrs[0] = *nobs_incrs_ptr;
      } else {
        ctx.thread_ohets[tidx] = R_CAST(uint32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[sample_ct_i32v * kBytesPerVec]);
        ctx.thread_ehet_incrs[tidx] = R_CAST(double*, cur_alloc);
        cur_alloc = &(cur_alloc[sample_ct_dv * kBytesPerVec]);
        ctx.thread_nobs_incrs[tidx] = R_CAST(int32_t*, cur_alloc);
        cur_alloc = &(cur_alloc[sample_ct_i32v * kBytesPerVec]);
        assert(cur_alloc <= g_bigstack_base);
      }
    }
    SetThreadFuncAndData(HetThread, &ctx, &tg);

    logprintf("%s: ", calcstr);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;

    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_size = MultireadNonempty(ctx.variant_subset, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto HetCalcMain_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, ctx.err_info);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, ctx.err_info >> 32);
          goto HetCalcMain_ret_1;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_size = cur_block_size;
        ComputeUidxStartPartition(ctx.variant_subset, cur_block_size, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_size == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto HetCalcMain_ret_THREAD_CREATE_FAIL;
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
      }

      ++read_block_idx;
      variant_idx += cur_block_size;
      pgfip->block_base = main_loadbufs[parity];
    }
    double ehet_base = ctx.thread_ehet_base[0];
    uint32_t nobs_base = ctx.thread_nobs_base[0];
    uint32_t monomorphic_ct = ctx.thread_monomorphic_ct[0];
    uint32_t* ohets = *ohets_ptr;
    double* ehet_incrs = *ehet_incrs_ptr;
    int32_t* nobs_incrs = *nobs_incrs_ptr;
    for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
      ehet_base += ctx.thread_ehet_base[tidx];
      nobs_base += ctx.thread_nobs_base[tidx];
      monomorphic_ct += ctx.thread_monomorphic_ct[tidx];
      U32CastVecAdd(ctx.thread_ohets[tidx], sample_ct_i32v, ohets);
      const double* src = ctx.thread_ehet_incrs[tidx];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        ehet_incrs[sample_idx] += src[sample_idx];
      }
      I32CastVecAdd(ctx.thread_nobs_incrs[tidx], sample_ct_i32v, nobs_incrs);
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    if (monomorphic_ct) {
      logerrprintfww("Warning: %u variant%s skipped because %s monomorphic.%s\n", monomorphic_ct, (monomorphic_ct == 1)? "" : "s", (monomorphic_ct == 1)? "it was" : "they were", small_sample? "" : " You may want to use --read-freq to provide more accurate allele frequency estimates.");
    }
    *ehet_base_ptr = ehet_base;
    *nobs_base_ptr = nobs_base;
  }
  while (0) {
  HetCalcMain_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  HetCalcMain_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  HetCalcMain_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 HetCalcMain_ret_1:
  CleanupThreads(&tg);
  BigstackReset(bigstack_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

PglErr HetReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const double* allele_freqs, const uintptr_t* founder_info, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t founder_ct, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_allele_ct, HetFlags flags, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    if (unlikely(IsSet(cip->haploid_mask, 0))) {
      logerrputs("Error: --het cannot be used on haploid genomes.\n");
      goto HetReport_ret_INCONSISTENT_INPUT;
    }
    const uint32_t small_sample = (flags / kfHetSmallSample) & 1;
    if (unlikely(small_sample && (!founder_ct))) {
      logerrputs("Error: '--het small-sample' requires founders.  (--make-founders may come in\nhandy here.\n)");
      goto HetReport_ret_INCONSISTENT_INPUT;
    }

    const uintptr_t* autosomal_variant_include = orig_variant_include;
    uint32_t autosomal_variant_ct = orig_variant_ct;
    reterr = ConditionalAllocateNonAutosomalVariants(cip, "--het", raw_variant_ct, &autosomal_variant_include, &autosomal_variant_ct);
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

    uint32_t* ohets;
    double* ehet_incrs;
    int32_t* nobs_incrs;
    double ehet_base;
    int32_t nobs_base;
    reterr = HetCalcMain(sample_include, autosomal_variant_include, allele_idx_offsets, allele_freqs, founder_info, "--het", raw_sample_ct, sample_ct, founder_ct, raw_variant_ct, autosomal_variant_ct, max_allele_ct, small_sample, max_thread_ct, pgr_alloc_cacheline_ct, pgfip, &ohets, &ehet_incrs, &nobs_incrs, &ehet_base, &nobs_base);
    if (unlikely(reterr)) {
      goto HetReport_ret_1;
    }
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      cswritep = AppendXid(sample_ids, sids, col_fid, col_sid, max_sample_id_blen, max_sid_blen, sample_uidx, cswritep);
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
        cswritep = dtoa_g(1.0 - u31tod(ohet) / ehet, cswritep);
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
  HetReport_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  HetReport_ret_NO_VARIATION:
    logerrputs("Error: --het requires at least one polymorphic autosomal variant.\n");
  HetReport_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 HetReport_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr CheckOrImputeSex(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const double* allele_freqs, const CheckSexInfo* csip, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t max_allele_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t* sex_nm, uintptr_t* sex_male, PgenFileInfo* pgfip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    const CheckSexFlags flags = csip->flags;
    const char* flagstr = (flags & kfCheckSexImpute)? "--impute-sex" : "--check-sex";
    if (unlikely(IsSet(cip->haploid_mask, 0))) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s cannot be used on haploid genomes.\n", flagstr);
      goto CheckOrImputeSex_ret_INCONSISTENT_INPUT_2;
    }

    double* xfs = nullptr;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    double max_female_xf = csip->max_female_xf;
    double min_male_xf = csip->min_male_xf;
    uint32_t used_variant_ct_x = 0;
    uint32_t x_code;
    if ((flags & kfCheckSexUseX) && XymtExists(cip, kChrOffsetX, &x_code)) {
      used_variant_ct_x = CountChrVariantsUnsafe(orig_variant_include, cip, x_code);
      if (used_variant_ct_x) {
        if (max_female_xf == -1.0) {
          max_female_xf = PrevFloat64(min_male_xf);
        } else if (min_male_xf == -1.0) {
          min_male_xf = NextFloat64(max_female_xf);
        }
        uintptr_t* variant_include_x;
        if (unlikely(bigstack_alloc_d(sample_ct, &xfs) ||
                     bigstack_alloc_w(raw_variant_ctl, &variant_include_x))) {
          goto CheckOrImputeSex_ret_NOMEM;
        }
        const uint32_t x_fo_idx = cip->chr_idx_to_foidx[x_code];
        const uint32_t x_start = cip->chr_fo_vidx_start[x_fo_idx];
        const uint32_t x_end = cip->chr_fo_vidx_start[x_fo_idx + 1];
        const uint32_t x_start_widx = x_start / kBitsPerWord;
        const uint32_t x_end_widx = (x_end + kBitsPerWord - 1) / kBitsPerWord;
        memcpy(&(variant_include_x[x_start_widx]), &(orig_variant_include[x_start_widx]), (x_end_widx - x_start_widx) * sizeof(intptr_t));
        if (x_start) {
          ClearBitsNz(0, x_start, variant_include_x);
        }
        // bugfix (17 Dec 2024): trailing bits must also be zeroed out
        const uint32_t raw_variant_ct_rounded_up = RoundUpPow2(raw_variant_ct, kBitsPerWord);
        if (x_end < raw_variant_ct_rounded_up) {
          ClearBitsNz(x_end, raw_variant_ct_rounded_up, variant_include_x);
        }
        // Don't actually need nobs.
        uint32_t* ohets;
        double* ehet_incrs;
        int32_t* nobs_incrs;
        double ehet_base;
        int32_t nobs_base;
        reterr = HetCalcMain(sample_include, variant_include_x, allele_idx_offsets, allele_freqs, nullptr, (flags & kfCheckSexImpute)? "--impute-sex chrX" : "--check-sex chrX", raw_sample_ct, sample_ct, sample_ct, raw_variant_ct, used_variant_ct_x, max_allele_ct, 0, max_thread_ct, pgr_alloc_cacheline_ct, pgfip, &ohets, &ehet_incrs, &nobs_incrs, &ehet_base, &nobs_base);
        if (unlikely(reterr)) {
          goto CheckOrImputeSex_ret_1;
        }
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          const uint32_t ohet = ohets[sample_idx];
          const double ehet = ehet_base + ehet_incrs[sample_idx];
          xfs[sample_idx] = 1.0 - u31tod(ohet) / ehet;
        }
        BigstackReset(variant_include_x);
      }
    }
    uint32_t* y_valid_geno_cts = nullptr;
    uint32_t max_female_ycount = csip->max_female_ycount;
    uint32_t min_male_ycount = csip->min_male_ycount;
    uint32_t used_variant_ct_y = 0;
    uint32_t y_code = 0;
    if ((flags & kfCheckSexUseY) && XymtExists(cip, kChrOffsetY, &y_code)) {
      used_variant_ct_y = CountChrVariantsUnsafe(orig_variant_include, cip, y_code);
      if (used_variant_ct_y) {
        if (max_female_ycount == UINT32_MAX) {
          if (min_male_ycount != UINT32_MAX) {
            max_female_ycount = min_male_ycount - 1;
          } else {
            // No ycount condition specified.  yrate condition will be merged
            // correctly if we initialize max_female_ycount=UINT32_MAX,
            // min_male_ycount=0.
            min_male_ycount = 0;
          }
        } else if (min_male_ycount == UINT32_MAX) {
          min_male_ycount = max_female_ycount + 1;
        }
        {
          const double max_female_yrate = csip->max_female_yrate;
          const double min_male_yrate = csip->min_male_yrate;
          if (max_female_yrate != -1.0) {
            const uint32_t uii = S_CAST(int32_t, max_female_yrate * used_variant_ct_y);
            if (uii < max_female_ycount) {
              max_female_ycount = uii;
            }
          }
          if (min_male_yrate != -1.0) {
            const uint32_t uii = used_variant_ct_y - S_CAST(int32_t, (1.0 - min_male_yrate) * used_variant_ct_y);
            if (uii > min_male_ycount) {
              min_male_ycount = uii;
            }
          }
        }
        uintptr_t* variant_include_y;
        uint32_t* sample_missing_hc_cts;
        uint32_t* sample_hethap_cts;
        if (unlikely(bigstack_alloc_u32(sample_ct, &y_valid_geno_cts) ||
                     bigstack_alloc_w(raw_variant_ctl, &variant_include_y) ||
                     bigstack_alloc_u32(raw_sample_ct, &sample_missing_hc_cts) ||
                     bigstack_alloc_u32(raw_sample_ct, &sample_hethap_cts))) {
          goto CheckOrImputeSex_ret_NOMEM;
        }
        const uint32_t y_fo_idx = cip->chr_idx_to_foidx[y_code];
        const uint32_t y_start = cip->chr_fo_vidx_start[y_fo_idx];
        const uint32_t y_end = cip->chr_fo_vidx_start[y_fo_idx + 1];
        const uint32_t y_start_widx = y_start / kBitsPerWord;
        const uint32_t y_end_widx = (y_end + kBitsPerWord - 1) / kBitsPerWord;
        memcpy(&(variant_include_y[y_start_widx]), &(orig_variant_include[y_start_widx]), (y_end_widx - y_start_widx) * sizeof(intptr_t));
        if (y_start) {
          ClearBitsNz(0, y_start, variant_include_y);
        }
        const uint32_t raw_variant_ct_rounded_up = RoundUpPow2(raw_variant_ct, kBitsPerWord);
        if (y_end < raw_variant_ct_rounded_up) {
          ClearBitsNz(y_end, raw_variant_ct_rounded_up, variant_include_y);
        }
        logprintf("%s: ", flagstr);
        reterr = LoadSampleMissingCts(sample_include, sample_include, sample_include, variant_include_y, cip, "chrY valid genotype call", raw_variant_ct, used_variant_ct_y, raw_sample_ct, 0, max_thread_ct, pgr_alloc_cacheline_ct, pgfip, sample_missing_hc_cts, nullptr, sample_hethap_cts);
        if (unlikely(reterr)) {
          goto CheckOrImputeSex_ret_1;
        }
        uintptr_t sample_uidx_base = 0;
        uintptr_t cur_bits = sample_include[0];
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          const uint32_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
          y_valid_geno_cts[sample_idx] = used_variant_ct_y - sample_missing_hc_cts[sample_uidx] - sample_hethap_cts[sample_uidx];
        }
        BigstackReset(variant_include_y);
      }
    }
    if (unlikely((!used_variant_ct_x) && (!used_variant_ct_y))) {
      snprintf(g_logbuf, kLogbufSize, "Error: %s: No usable variants.\n", flagstr);
      goto CheckOrImputeSex_ret_INCONSISTENT_INPUT_2;
    }
    uintptr_t* imputed_sex_nm;
    uintptr_t* imputed_sex_male;
    if (unlikely(bigstack_calloc_w(raw_sample_ctl, &imputed_sex_nm) ||
                 bigstack_calloc_w(raw_sample_ctl, &imputed_sex_male))) {
      goto CheckOrImputeSex_ret_NOMEM;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".sexcheck");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto CheckOrImputeSex_ret_OPEN_FAIL;
    }
    char* textbuf_flush = &(g_textbuf[kMaxMediumLine]);
    char* write_iter = g_textbuf;
    *write_iter++ = '#';
    const uint32_t col_fid = FidColIsRequired(siip, flags / kfCheckSexColMaybefid);
    if (col_fid) {
      write_iter = strcpya_k(write_iter, "FID\t");
    }
    write_iter = strcpya_k(write_iter, "IID");
    const uint32_t col_sid = SidColIsRequired(siip->sids, flags / kfCheckSexColMaybesid);
    if (col_sid) {
      write_iter = strcpya_k(write_iter, "\tSID");
    }
    const uint32_t col_pedsex = (flags / kfCheckSexColPedsex) & 1;
    if (col_pedsex) {
      write_iter = strcpya_k(write_iter, "\tPEDSEX");
    }
    write_iter = strcpya_k(write_iter, "\tSNPSEX");
    const uint32_t col_status = (flags / kfCheckSexColStatus) & 1;
    if (col_status) {
      write_iter = strcpya_k(write_iter, "\tSTATUS");
    }
    const uint32_t col_xf = used_variant_ct_x && (flags & kfCheckSexColXF);
    if (col_xf) {
      write_iter = strcpya_k(write_iter, "\tF");
    }
    const uint32_t col_ycount = used_variant_ct_y && (flags & kfCheckSexColYcount);
    if (col_ycount) {
      write_iter = strcpya_k(write_iter, "\tYCOUNT");
    }
    const uint32_t col_yrate = used_variant_ct_y && (flags & kfCheckSexColYrate);
    if (col_yrate) {
      write_iter = strcpya_k(write_iter, "\tYRATE");
    }
    const uint32_t col_yobs = used_variant_ct_y && (flags & kfCheckSexColYobs);
    if (col_yobs) {
      write_iter = strcpya_k(write_iter, "\tYOBS");
    }
    AppendBinaryEoln(&write_iter);

    const double y_denom = used_variant_ct_y? (1.0 / u31tod(used_variant_ct_y)) : 0.0;
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    double cur_xf = 0.0;
    uint32_t cur_ycount = 0;
    uint32_t problem_ct = 0;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      write_iter = AppendXid(sample_ids, sids, col_fid, col_sid, max_sample_id_blen, max_sid_blen, sample_uidx, write_iter);
      const uint32_t ped_nm = IsSet(sex_nm, sample_uidx);
      const uint32_t ped_male = IsSet(sex_male, sample_uidx);
      if (col_pedsex) {
        *write_iter++ = '\t';
        if (ped_nm) {
          *write_iter++ = '2' - ped_male;
        } else {
          write_iter = strcpya_k(write_iter, "NA");
        }
      }
      uint32_t not_male = 0;
      uint32_t not_female = 0;
      if (xfs) {
        cur_xf = xfs[sample_idx];
        if (cur_xf <= max_female_xf) {
          not_male = 1;
        } else if (cur_xf >= min_male_xf) {
          not_female = 1;
        } else {
          not_male = 1;
          not_female = 1;
        }
      }
      if (y_valid_geno_cts) {
        cur_ycount = y_valid_geno_cts[sample_idx];
        if (cur_ycount <= max_female_ycount) {
          not_male = 1;
        } else if (cur_ycount >= min_male_ycount) {
          not_female = 1;
        } else {
          not_male = 1;
          not_female = 1;
        }
      }
      *write_iter++ = '\t';
      uint32_t is_problem;
      if (not_male) {
        if (not_female) {
          write_iter = strcpya_k(write_iter, "NA");
          is_problem = 1;
        } else {
          *write_iter++ = '2';
          SetBit(sample_uidx, imputed_sex_nm);
          is_problem = (!ped_nm) || ped_male;
        }
      } else {
        *write_iter++ = '1';
        SetBit(sample_uidx, imputed_sex_nm);
        SetBit(sample_uidx, imputed_sex_male);
        is_problem = (!ped_nm) || (!ped_male);
      }
      problem_ct += is_problem;
      if (col_status) {
        if (is_problem) {
          write_iter = strcpya_k(write_iter, "\tPROBLEM");
        } else {
          write_iter = strcpya_k(write_iter, "\tOK");
        }
      }
      if (col_xf) {
        *write_iter++ = '\t';
        write_iter = dtoa_g(cur_xf, write_iter);
      }
      if (col_ycount) {
        *write_iter++ = '\t';
        write_iter = u32toa(cur_ycount, write_iter);
      }
      if (col_yrate) {
        *write_iter++ = '\t';
        write_iter = dtoa_g(u31tod(cur_ycount) * y_denom, write_iter);
      }
      if (col_yobs) {
        *write_iter++ = '\t';
        write_iter = u32toa(used_variant_ct_y, write_iter);
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto CheckOrImputeSex_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto CheckOrImputeSex_ret_WRITE_FAIL;
    }
    write_iter = strcpya(g_logbuf, flagstr);
    write_iter = strcpya_k(write_iter, ": ");
    if (used_variant_ct_x) {
      write_iter = u32toa(used_variant_ct_x, write_iter);
      write_iter = strcpya_k(write_iter, " chrX variant");
      if (used_variant_ct_x != 1) {
        *write_iter++ = 's';
      }
      if (used_variant_ct_y) {
        write_iter = strcpya_k(write_iter, " and ");
        write_iter = u32toa(used_variant_ct_y, write_iter);
        write_iter = strcpya_k(write_iter, " variant");
        if (used_variant_ct_y != 1) {
          *write_iter++ = 's';
        }
      }
    } else {
      write_iter = u32toa(used_variant_ct_y, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (used_variant_ct_y != 1) {
        *write_iter++ = 's';
      }
    }
    write_iter = strcpya_k(write_iter, " scanned, ");
    if (flags & kfCheckSexImpute) {
      memcpy(sex_nm, imputed_sex_nm, raw_sample_ctl * sizeof(intptr_t));
      memcpy(sex_male, imputed_sex_male, raw_sample_ctl * sizeof(intptr_t));
      const uint32_t imputed_nm_ct = PopcountWords(imputed_sex_nm, raw_sample_ctl);
      const uint32_t imputed_male_ct = PopcountWords(imputed_sex_male, raw_sample_ctl);
      write_iter = u32toa_x(imputed_nm_ct, '/', write_iter);
      write_iter = u32toa(sample_ct, write_iter);
      write_iter = strcpya_k(write_iter, " sex");
      if (sample_ct != 1) {
        write_iter = strcpya_k(write_iter, "es");
      }
      write_iter = strcpya_k(write_iter, " imputed (");
      write_iter = u32toa(imputed_nm_ct - imputed_male_ct, write_iter);
      write_iter = strcpya_k(write_iter, " female, ");
      write_iter = u32toa(imputed_male_ct, write_iter);
      write_iter = strcpya_k(write_iter, " male)");
    } else {
      write_iter = u32toa(problem_ct, write_iter);
      write_iter = strcpya_k(write_iter, " problem");
      if (problem_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " detected");
    }
    write_iter = strcpya_k(write_iter, ". Report written to ");
    write_iter = memcpya(write_iter, outname, &(outname_end[strlen(".sexcheck")]) - outname);
    strcpy_k(write_iter, " .\n");
    WordWrapB(0);
    logputsb();
  }
  while (0) {
  CheckOrImputeSex_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CheckOrImputeSex_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  CheckOrImputeSex_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  CheckOrImputeSex_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 CheckOrImputeSex_ret_1:
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct FstCtxStruct {
  const uintptr_t* cur_variant_include;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  // sex_male_collapsed == nullptr unless chrX /w males
  const uintptr_t* sex_male_collapsed;
  const uint32_t* sample_to_pop_idx;
  uint32_t* diploid_pop_sizes;
  const uint32_t* haploid_pop_sizes;
  uint32_t sample_ct;
  uint32_t pop_ct;
  uint32_t max_difflist_len;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** raregenos;
  uint32_t** difflist_sample_id_bufs;
  uint32_t** pop_geno_bufs;
  uint32_t* read_variant_uidx_starts;
  uint32_t cur_block_variant_ct;  // .pgen iteration, not jackknife

  uint64_t err_info;

  // worker threads compute intermediate stats; parent thread computes final
  // Fst estimates (to avoid O(pop_ct^2) space problem) and executes
  // block-jackknife.
  uint32_t* pop_nm_sample_cts[2];
  uint32_t* pop_allele_obs_cts[2];  // major dimension = pop, minor = allele

  // Hudson only
  double* half_within[2];

  // W-C only
  uint32_t* pop_allele_het_cts[2];
} FstCtx;

THREAD_FUNC_DECL FstThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  FstCtx* ctx = S_CAST(FstCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->cur_variant_include;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uintptr_t* sample_include = ctx->sample_include;
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  // sex_male nullptr iff not chrX /w males
  const uintptr_t* sex_male = ctx->sex_male_collapsed;
  const uint32_t* sample_to_pop_idx = ctx->sample_to_pop_idx;
  const uint32_t* diploid_pop_sizes = ctx->diploid_pop_sizes;
  const uint32_t* haploid_pop_sizes = ctx->haploid_pop_sizes;
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
  const uint32_t haploid_present = (sex_male != nullptr);
  const uintptr_t pop_ct = ctx->pop_ct;
  const uintptr_t pop_ct_x4 = pop_ct * 4;
  uintptr_t* genovec = ctx->genovecs[tidx];
  PgenVariant pgv;
  pgv.genovec = genovec;
  SetPgvThreadMhcNull(sample_ct, tidx, ctx->thread_read_mhc, &pgv);
  uintptr_t* raregeno = ctx->raregenos[tidx];
  uint32_t* difflist_sample_ids = ctx->difflist_sample_id_bufs[tidx];
  uint32_t* pop_geno_buf = ctx->pop_geno_bufs[tidx];
  const uintptr_t pop_geno_buf_size = pop_ct_x4 << haploid_present;
  const uint32_t max_difflist_len = ctx->max_difflist_len;

  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);

  uint32_t allele_ct = 2;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    uint32_t cur_idx_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    uint32_t* pop_nm_sample_cts_iter;
    uint32_t* pop_allele_obs_cts_iter;
    double* half_within_iter;
    uint32_t* pop_allele_het_cts_iter;
    {
      const uint32_t cur_block_variant_ct = ctx->cur_block_variant_ct;
      const uint32_t first_variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
      cur_idx_ct = (((tidx + 1) * cur_block_variant_ct) / calc_thread_ct) - first_variant_bidx;
      const uint32_t variant_uidx_start = ctx->read_variant_uidx_starts[tidx];
      BitIter1Start(variant_include, variant_uidx_start, &variant_uidx_base, &variant_include_bits);
      const uintptr_t allele_bidx = first_variant_bidx * 2 + CountExtraAlleles(variant_include, allele_idx_offsets, ctx->read_variant_uidx_starts[0], variant_uidx_start, 0);
      // non-null iff Weir-Cockerham method, or OBS_CT .fst.var col requested
      pop_nm_sample_cts_iter = ctx->pop_nm_sample_cts[parity];
      if (pop_nm_sample_cts_iter) {
        pop_nm_sample_cts_iter = &(pop_nm_sample_cts_iter[pop_ct * first_variant_bidx]);
      }
      pop_allele_obs_cts_iter = ctx->pop_allele_obs_cts[parity];
      pop_allele_obs_cts_iter = &(pop_allele_obs_cts_iter[pop_ct * allele_bidx]);
      // non-null iff Hudson
      half_within_iter = ctx->half_within[parity];
      if (half_within_iter) {
        half_within_iter = &(half_within_iter[pop_ct * first_variant_bidx]);
      }

      // non-null iff W-C
      pop_allele_het_cts_iter = ctx->pop_allele_het_cts[parity];
      if (pop_allele_het_cts_iter) {
        pop_allele_het_cts_iter = &(pop_allele_het_cts_iter[pop_ct * allele_bidx]);
      }
    }
    for (uint32_t cur_idx = 0; cur_idx != cur_idx_ct; ++cur_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      if (allele_idx_offsets) {
        allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
      }
      ZeroU32Arr(pop_geno_buf_size, pop_geno_buf);
      uint32_t difflist_common_geno = UINT32_MAX;
      if (allele_ct == 2) {
        uint32_t difflist_len;
        const PglErr reterr = PgrGetDifflistOrGenovec(sample_include, pssi, sample_ct, max_difflist_len, variant_uidx, pgrp, genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto FstThread_err;
        }
        if (difflist_common_geno != UINT32_MAX) {
          const uint32_t word_ct = NypCtToWordCt(difflist_len);
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
              if (!haploid_present) {
                for (uint32_t uii = 0; uii != loop_len; ++uii) {
                  const uint32_t sample_idx = cur_difflist_sample_ids[uii];
                  const uintptr_t cur_geno = raregeno_word & 3;
                  const uintptr_t pop_idx = sample_to_pop_idx[sample_idx];
                  pop_geno_buf[pop_idx * 4 + cur_geno] += 1;
                  raregeno_word >>= 2;
                }
              } else {
                for (uint32_t uii = 0; uii != loop_len; ++uii) {
                  const uint32_t sample_idx = cur_difflist_sample_ids[uii];
                  uintptr_t cur_geno = raregeno_word & 3;
                  const uintptr_t pop_idx = sample_to_pop_idx[sample_idx];
                  const uintptr_t male_mask = -IsSet(sex_male, sample_idx);
                  pop_geno_buf[pop_idx * 4 + cur_geno + (male_mask & pop_ct_x4)] += 1;
                  raregeno_word >>= 2;
                }
              }
            }
          }
          if (difflist_common_geno != 3) {
            const uint32_t other_hom_geno = 2 - difflist_common_geno;
            for (uintptr_t pop_idx = 0; pop_idx != pop_ct; ++pop_idx) {
              const uint32_t diploid_pop_size = diploid_pop_sizes[pop_idx];
              uint32_t* cur_pop_geno_buf = &(pop_geno_buf[pop_idx * 4]);
              cur_pop_geno_buf[difflist_common_geno] = diploid_pop_size - cur_pop_geno_buf[1] - cur_pop_geno_buf[other_hom_geno] - cur_pop_geno_buf[3];
            }
            if (haploid_present) {
              uint32_t* haploid_pop_geno_buf = &(pop_geno_buf[pop_ct_x4]);
              for (uintptr_t pop_idx = 0; pop_idx != pop_ct; ++pop_idx) {
                const uint32_t haploid_pop_size = haploid_pop_sizes[pop_idx];
                uint32_t* cur_pop_geno_buf = &(haploid_pop_geno_buf[pop_idx * 4]);
                cur_pop_geno_buf[difflist_common_geno] = haploid_pop_size - cur_pop_geno_buf[1] - cur_pop_geno_buf[other_hom_geno] - cur_pop_geno_buf[3];
              }
            }
          }
        }
      } else {
        const PglErr reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, pgrp, &pgv);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto FstThread_err;
        }
      }
      if (difflist_common_geno == UINT32_MAX) {
        uint32_t loop_len = kBitsPerWordD2;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            loop_len = ModNz(sample_ct, kBitsPerWordD2);
          }
          const uint32_t* cur_pop_idxs = &(sample_to_pop_idx[widx * kBitsPerWordD2]);
          uintptr_t geno_word = genovec[widx];
          if (!haploid_present) {
            for (uint32_t uii = 0; uii != loop_len; ++uii) {
              const uintptr_t cur_geno = geno_word & 3;
              const uintptr_t pop_idx = cur_pop_idxs[uii];
              pop_geno_buf[pop_idx * 4 + cur_geno] += 1;
              geno_word >>= 2;
            }
          } else {
            uintptr_t male_hw = R_CAST(const Halfword*, sex_male)[widx];
            for (uint32_t uii = 0; uii != loop_len; ++uii) {
              const uintptr_t cur_geno = geno_word & 3;
              const uintptr_t pop_idx = cur_pop_idxs[uii];
              const uintptr_t male_mask = -(male_hw & k1LU);
              pop_geno_buf[pop_idx * 4 + cur_geno + (male_mask & pop_ct_x4)] += 1;
              geno_word >>= 2;
              male_hw >>= 1;
            }
          }
        }
      }
      if (pop_nm_sample_cts_iter) {
        const uint32_t* pop_geno_buf_iter = pop_geno_buf;
        for (uintptr_t pop_idx = 0; pop_idx != pop_ct; ++pop_idx) {
          pop_nm_sample_cts_iter[pop_idx] = pop_geno_buf_iter[0] + pop_geno_buf_iter[1] + pop_geno_buf_iter[2];
          pop_geno_buf_iter = &(pop_geno_buf_iter[4]);
        }
        if (haploid_present) {
          for (uintptr_t pop_idx = 0; pop_idx != pop_ct; ++pop_idx) {
            pop_nm_sample_cts_iter[pop_idx] += pop_geno_buf_iter[0] + pop_geno_buf_iter[2];
            pop_geno_buf_iter = &(pop_geno_buf_iter[4]);
          }
        }
      }
      ZeroU32Arr(pop_ct * allele_ct, pop_allele_obs_cts_iter);
      if (pop_allele_het_cts_iter) {
        ZeroU32Arr(pop_ct * allele_ct, pop_allele_het_cts_iter);
      }
      if (allele_ct != 2) {
        if (pgv.patch_01_ct) {
          const uint32_t patch_01_ct = pgv.patch_01_ct;
          const uintptr_t* patch_01_set = pgv.patch_01_set;
          const AlleleCode* patch_01_vals = pgv.patch_01_vals;
          uintptr_t sample_idx_base = 0;
          uintptr_t sample_idx_bits = patch_01_set[0];
          for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
            const uint32_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &sample_idx_bits);
            if (sex_male && IsSet(sex_male, sample_idx)) {
              // already counted as missing
              continue;
            }
            const uintptr_t pop_idx = sample_to_pop_idx[sample_idx];
            pop_geno_buf[pop_idx * 4 + 1] -= 1;
            pop_allele_obs_cts_iter[pop_idx] += 1;
            const AlleleCode ac = patch_01_vals[uii];
            pop_allele_obs_cts_iter[ac * pop_ct + pop_idx] += 1;
            if (pop_allele_het_cts_iter) {
              pop_allele_het_cts_iter[pop_idx] += 1;
              pop_allele_het_cts_iter[ac * pop_ct + pop_idx] += 1;
            }
          }
        }
        if (pgv.patch_10_ct) {
          const uint32_t patch_10_ct = pgv.patch_10_ct;
          const uintptr_t* patch_10_set = pgv.patch_10_set;
          const AlleleCode* patch_10_vals = pgv.patch_10_vals;
          uintptr_t sample_idx_base = 0;
          uintptr_t sample_idx_bits = patch_10_set[0];
          if (!haploid_present) {
            for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
              const uint32_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
              const uintptr_t pop_idx = sample_to_pop_idx[sample_idx];
              pop_geno_buf[pop_idx * 4 + 2] -= 1;
              const AlleleCode ac0 = patch_10_vals[uii * 2];
              const AlleleCode ac1 = patch_10_vals[uii * 2 + 1];
              pop_allele_obs_cts_iter[pop_idx + ac0 * pop_ct] += 1;
              pop_allele_obs_cts_iter[pop_idx + ac1 * pop_ct] += 1;
              if (pop_allele_het_cts_iter && (ac0 != ac1)) {
                pop_allele_het_cts_iter[pop_idx + ac0 * pop_ct] += 1;
                pop_allele_het_cts_iter[pop_idx + ac1 * pop_ct] += 1;
              }
            }
          } else {
            for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
              const uint32_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &sample_idx_bits);
              const uintptr_t pop_idx = sample_to_pop_idx[sample_idx];
              const uintptr_t is_male = IsSet(sex_male, sample_idx);
              pop_geno_buf[pop_idx * 4 + 2 + ((-is_male) & pop_ct_x4)] -= 1;
              const AlleleCode ac0 = patch_10_vals[uii * 2];
              const AlleleCode ac1 = patch_10_vals[uii * 2 + 1];
              if (ac0 == ac1) {
                pop_allele_obs_cts_iter[pop_idx + ac0 * pop_ct] += 2 - is_male;
              } else {
                if (is_male) {
                  if (pop_nm_sample_cts_iter) {
                    pop_nm_sample_cts_iter[pop_idx] -= 1;
                  }
                  continue;
                }
                pop_allele_obs_cts_iter[pop_idx + ac0 * pop_ct] += 1;
                pop_allele_obs_cts_iter[pop_idx + ac1 * pop_ct] += 1;
                // guaranteed to be performing Hudson computation if we're
                // here, so no need to check pop_allele_het_cts_iter
              }
            }
          }
        }
      }
      for (uintptr_t pop_idx = 0; pop_idx != pop_ct; ++pop_idx) {
        const uint32_t hom_ref_ct = pop_geno_buf[pop_idx * 4];
        const uint32_t het_ref_alt1_ct = pop_geno_buf[pop_idx * 4 + 1];
        const uint32_t hom_alt1_ct = pop_geno_buf[pop_idx * 4 + 2];
        pop_allele_obs_cts_iter[pop_idx] += hom_ref_ct * 2 + het_ref_alt1_ct;
        pop_allele_obs_cts_iter[pop_idx + pop_ct] += hom_alt1_ct * 2 + het_ref_alt1_ct;
        if (pop_allele_het_cts_iter) {
          pop_allele_het_cts_iter[pop_idx] += het_ref_alt1_ct;
          pop_allele_het_cts_iter[pop_idx + pop_ct] += het_ref_alt1_ct;
        }
      }
      if (haploid_present) {
        const uint32_t* haploid_pop_geno_buf_iter = &(pop_geno_buf[pop_ct_x4]);
        for (uintptr_t pop_idx = 0; pop_idx != pop_ct; ++pop_idx) {
          const uint32_t hap_ref_ct = haploid_pop_geno_buf_iter[0];
          const uint32_t hap_alt1_ct = haploid_pop_geno_buf_iter[2];
          pop_allele_obs_cts_iter[pop_idx] += hap_ref_ct;
          pop_allele_obs_cts_iter[pop_idx + pop_ct] += hap_alt1_ct;
          haploid_pop_geno_buf_iter = &(haploid_pop_geno_buf_iter[4]);
        }
      }
      if (half_within_iter) {
        // Hudson
        for (uintptr_t pop_idx = 0; pop_idx != pop_ct; ++pop_idx) {
          uintptr_t n_hap = 0;
          uint64_t ssq = 0;
          for (uintptr_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            const uintptr_t allele_obs_ct = pop_allele_obs_cts_iter[allele_idx * pop_ct + pop_idx];
            n_hap += allele_obs_ct;
            ssq += S_CAST(uint64_t, allele_obs_ct) * allele_obs_ct;
          }
          const uint64_t n_pairs_x2 = S_CAST(uint64_t, n_hap) * (n_hap - 1);
          const uint64_t n_same = (ssq - n_hap) / 2;
          const uint64_t n_diff = n_pairs_x2 / 2 - n_same;
          *half_within_iter++ = u63tod(n_diff) / u63tod(n_pairs_x2);
        }
      } else {
        // Weir-Cockerham
        pop_allele_het_cts_iter = &(pop_allele_het_cts_iter[pop_ct * allele_ct]);
      }
      if (pop_nm_sample_cts_iter) {
        pop_nm_sample_cts_iter = &(pop_nm_sample_cts_iter[pop_ct]);
      }
      pop_allele_obs_cts_iter = &(pop_allele_obs_cts_iter[pop_ct * allele_ct]);
    }
    parity = 1 - parity;
    while (0) {
    FstThread_err:
      UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

const char g_cc_cat_names[2][8] = {"CONTROL", "CASE"};

CONSTI32(kFstReportVariantsBatchMax, kMaxOpenFiles - 12);

PglErr FstReport(const uintptr_t* orig_sample_include, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const FstInfo* fst_infop, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_allele_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uint32_t max_variant_file_ct = 0;
  FILE* s_outfile = nullptr;
  char** v_cswritep_arr = nullptr;
  CompressStreamState* v_css_arr = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs; // file=
  PreinitTextStream(&txs);
  ThreadGroup tg;
  PreinitThreads(&tg);
  FstCtx ctx;
  {
    if (IsSet(cip->haploid_mask, 0)) {
      logerrputs("Error: --fst cannot be used on haploid genomes.\n");
      goto FstReport_ret_INCONSISTENT_INPUT;
    }
    const char* pheno_name = fst_infop->pheno_name;
    const PhenoCol* pheno_col;
    {
      const uint32_t pheno_blen = strlen(pheno_name) + 1;
      if (pheno_blen > max_pheno_name_blen) {
        goto FstReport_ret_PHENO_NOT_FOUND;
      }
      for (uintptr_t pheno_idx = 0; ; ++pheno_idx) {
        if (pheno_idx == pheno_ct) {
          goto FstReport_ret_PHENO_NOT_FOUND;
        }
        if (memequal(pheno_name, &(pheno_names[pheno_idx * max_pheno_name_blen]), pheno_blen)) {
          pheno_col = &(pheno_cols[pheno_idx]);
          break;
        }
      }
    }
    if (pheno_col->type_code == kPhenoDtypeQt) {
      logerrprintfww("Error: --fst phenotype '%s' is quantitative (binary or categorical required).\n", pheno_name);
      goto FstReport_ret_INCONSISTENT_INPUT;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    if (pheno_col->type_code != kPhenoDtypeCat) {
      assert(pheno_col->type_code == kPhenoDtypeCc);
      // this may belong in plink2_cmdline
      PhenoCol* synthetic_pheno_col;
      uint32_t* cat_tmp;
      const char** category_names;
      if (unlikely(BIGSTACK_ALLOC_X(PhenoCol, 1, &synthetic_pheno_col) ||
                   bigstack_end_calloc_u32(raw_sample_ct, &cat_tmp) ||
                   bigstack_end_alloc_kcp(3, &category_names))) {
        goto FstReport_ret_NOMEM;
      }
      const uintptr_t* raw_pheno_nm = pheno_col->nonmiss;
      const uintptr_t* pheno_cc = pheno_col->data.cc;
      for (uint32_t sample_uidx = 0; sample_uidx != raw_sample_ct; ++sample_uidx) {
        if (!IsSet(raw_pheno_nm, sample_uidx)) {
          continue;
        }
        // 'CASE' is lexicographically before 'CONTROL'
        cat_tmp[sample_uidx] = 2 - IsSet(pheno_cc, sample_uidx);
      }
      category_names[0] = nullptr;
      category_names[1] = g_cc_cat_names[1];
      category_names[2] = g_cc_cat_names[0];
      synthetic_pheno_col->category_names = category_names;
      synthetic_pheno_col->nonmiss = pheno_col->nonmiss;
      synthetic_pheno_col->data.cat = cat_tmp;
      synthetic_pheno_col->type_code = kPhenoDtypeCat;
      synthetic_pheno_col->nonnull_category_ct = 2;
      pheno_col = synthetic_pheno_col;
    }
    ctx.allele_idx_offsets = allele_idx_offsets;
    uintptr_t* sample_include;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &sample_include))) {
      goto FstReport_ret_NOMEM;
    }
    BitvecAndCopy(orig_sample_include, pheno_col->nonmiss, raw_sample_ctl, sample_include);
    uint32_t sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    const FstFlags flags = fst_infop->flags;
    const uint32_t is_wc = (flags / kfFstMethodWc) & 1;
    uintptr_t pop_pair_ct = 0;
    uint32_t* pop_pairs = nullptr;
    const char** pop_names;
    uintptr_t pop_ct;
    {
      const uint32_t nonnull_category_ct = pheno_col->nonnull_category_ct;
      const uint32_t raw_cat_ctl = 1 + (nonnull_category_ct / kBitsPerWord);
      uintptr_t* cats_seen;
      if (unlikely(bigstack_end_alloc_w(raw_cat_ctl, &cats_seen))) {
        goto FstReport_ret_NOMEM;
      }
      pop_ct = IdentifyRemainingCats(sample_include, pheno_col, sample_ct, cats_seen);
      assert(!(cats_seen[0] & 1));
      if (pop_ct < 2) {
        logerrputs("Error: --fst requires two or more populations.\n");
        goto FstReport_ret_INCONSISTENT_INPUT;
      }
      const uint32_t* old_cats = pheno_col->data.cat;
      if (flags & (kfFstOneBasePop | kfFstExplicitPopIds)) {
        const char** nonnull_category_names = &(pheno_col->category_names[1]);
        const char* first_id = fst_infop->first_id_or_fname;
        const uint32_t first_cat_idx = 1 + bsearch_strptr_natural(first_id, nonnull_category_names, nonnull_category_ct);
        if (unlikely(!IsSet(cats_seen, first_cat_idx))) {
          logerrprintfww("Error: --fst phenotype '%s' does not have a nonempty population with ID '%s'.\n", pheno_name, first_id);
          goto FstReport_ret_INCONSISTENT_INPUT;
        }
        const char* other_ids = fst_infop->other_ids_flattened;
        if (!other_ids) {
          // base= with single argument
          const uintptr_t other_pop_ct = pop_ct - 1;
          if (unlikely(bigstack_alloc_u32(2 * other_pop_ct, &pop_pairs))) {
            goto FstReport_ret_NOMEM;
          }
          pop_pair_ct = other_pop_ct;
          uint32_t* pop_pairs_iter = pop_pairs;
          uintptr_t cats_seen_base = 0;
          uintptr_t cats_seen_bits = cats_seen[0];
          uintptr_t ulii = 0;
          for (; ulii != other_pop_ct; ++ulii) {
            const uint32_t other_cat_idx = BitIter1(cats_seen, &cats_seen_base, &cats_seen_bits);
            if (other_cat_idx == first_cat_idx) {
              break;
            }
            *pop_pairs_iter++ = other_cat_idx;
            *pop_pairs_iter++ = first_cat_idx;
          }
          for (; ulii != other_pop_ct; ++ulii) {
            const uint32_t other_cat_idx = BitIter1(cats_seen, &cats_seen_base, &cats_seen_bits);
            *pop_pairs_iter++ = first_cat_idx;
            *pop_pairs_iter++ = other_cat_idx;
          }
        } else {
          // either base= or id= with multiple arguments
          uintptr_t* pop_include;
          if (unlikely(bigstack_end_calloc_w(raw_cat_ctl, &pop_include))) {
            goto FstReport_ret_NOMEM;
          }
          SetBit(first_cat_idx, pop_include);
          const char* other_ids_iter = other_ids;
          do {
            const uint32_t cat_idx = 1 + bsearch_strptr_natural(other_ids_iter, nonnull_category_names, nonnull_category_ct);
            if (unlikely(!IsSet(cats_seen, cat_idx))) {
              logerrprintfww("Error: --fst phenotype '%s' does not have a nonempty population with ID '%s'.\n", pheno_name, other_ids_iter);
              goto FstReport_ret_INCONSISTENT_INPUT;
            }
            if (unlikely(IsSet(pop_include, cat_idx))) {
              logerrprintfww("Error: Duplicate --fst population ID '%s'.\n", other_ids_iter);
              goto FstReport_ret_INVALID_CMDLINE;
            }
            SetBit(cat_idx, pop_include);
            other_ids_iter = strnul(other_ids_iter);
            ++other_ids_iter;
          } while (*other_ids_iter);
          const uintptr_t named_pop_ct = PopcountWords(pop_include, raw_cat_ctl);
          if (named_pop_ct != pop_ct) {
            // Safe to ignore unnamed populations.
            pop_ct = named_pop_ct;
            memcpy(cats_seen, pop_include, raw_cat_ctl * sizeof(intptr_t));
            sample_ct = RemoveExcludedCats(old_cats, cats_seen, raw_sample_ct, sample_ct, sample_include);
          }
          if (flags & kfFstOneBasePop) {
            pop_pair_ct = named_pop_ct - 1;
            if (unlikely(bigstack_alloc_u32(2 * pop_pair_ct, &pop_pairs))) {
              goto FstReport_ret_NOMEM;
            }
            uint32_t* pop_pairs_iter = pop_pairs;
            uint32_t uii = 0;
            uintptr_t cats_seen_base = 0;
            uintptr_t cats_seen_bits = cats_seen[0];
            for (; uii != named_pop_ct; ++uii) {
              uint32_t other_cat_idx = BitIter1(cats_seen, &cats_seen_base, &cats_seen_bits);
              if (first_cat_idx == other_cat_idx) {
                while (1) {
                  ++uii;
                  if (uii == named_pop_ct) {
                    break;
                  }
                  other_cat_idx = BitIter1(cats_seen, &cats_seen_base, &cats_seen_bits);
                  *pop_pairs_iter++ = first_cat_idx;
                  *pop_pairs_iter++ = other_cat_idx;
                }
                break;
              }
              *pop_pairs_iter++ = other_cat_idx;
              *pop_pairs_iter++ = first_cat_idx;
            }
          } else {
            // ids=
            SetBit(first_cat_idx, cats_seen);
            const uint64_t named_pop_ct_m1 = named_pop_ct - k1LU;
            const uint64_t pop_pair_ct_u64 = (named_pop_ct * named_pop_ct_m1) / 2;
            if (unlikely(bigstack_alloc64_u32(2 * pop_pair_ct_u64, &pop_pairs))) {
              goto FstReport_ret_NOMEM;
            }
            pop_pair_ct = pop_pair_ct_u64;
            uint32_t* pop_pairs_iter = pop_pairs;
            uintptr_t cats_seen_base1 = 0;
            uintptr_t cats_seen_bits1 = cats_seen[0];
            for (uintptr_t ulii = 0; ulii != named_pop_ct_m1; ++ulii) {
              const uint32_t cat_idx1 = BitIter1(cats_seen, &cats_seen_base1, &cats_seen_bits1);
              uintptr_t cats_seen_base2 = cats_seen_base1;
              uintptr_t cats_seen_bits2 = cats_seen_bits1;
              for (uintptr_t uljj = ulii; uljj != named_pop_ct_m1; ++uljj) {
                const uint32_t cat_idx2 = BitIter1(cats_seen, &cats_seen_base2, &cats_seen_bits2);
                *pop_pairs_iter++ = cat_idx1;
                *pop_pairs_iter++ = cat_idx2;
              }
            }
          }
        }
      } else if (flags & kfFstPopPairFile) {
        const uint64_t possible_pair_ct_u64 = S_CAST(uint64_t, nonnull_category_ct) * (nonnull_category_ct + 1LLU) / 2;
        const uint64_t possible_pair_ctl_u64 = DivUpU64(possible_pair_ct_u64, kBitsPerWord);
        uintptr_t* selected_pairs;
        uintptr_t* pop_include;
        if (unlikely(bigstack_end_calloc64_w(possible_pair_ctl_u64, &selected_pairs) ||
                     bigstack_end_calloc_w(raw_cat_ctl, &pop_include))) {
          goto FstReport_ret_NOMEM;
        }
        unsigned char* bigstack_end_mark2 = g_bigstack_end;
        const char* in_fname = fst_infop->first_id_or_fname;
        reterr = InitTextStreamEx(in_fname, 1, kMaxLongLine, kTextStreamBlenFast, 1, &txs);
        if (unlikely(reterr)) {
          goto FstReport_ret_TSTREAM_FAIL;
        }
        const char** nonnull_category_names = &(pheno_col->category_names[1]);
        uintptr_t duplicate_ct = 0;
        uintptr_t line_idx = 0;
        while (1) {
          ++line_idx;
          char* first_token_start = TextGet(&txs);
          if (!first_token_start) {
            if (likely(!TextStreamErrcode2(&txs, &reterr))) {
              break;
            }
            goto FstReport_ret_TSTREAM_FAIL;
          }
          char* first_token_end = CurTokenEnd(first_token_start);
          char* second_token_start = FirstNonTspace(first_token_end);
          if (unlikely(IsEolnKns(*second_token_start))) {
            logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, in_fname);
            goto FstReport_ret_MALFORMED_INPUT;
          }
          char* second_token_end = CurTokenEnd(second_token_start);
          *first_token_end = '\0';
          const uintptr_t cat_idx1 = 1 + bsearch_strptr_natural(first_token_start, nonnull_category_names, nonnull_category_ct);
          if (unlikely(!IsSet(cats_seen, cat_idx1))) {
            logerrprintfww("Error: --fst phenotype '%s' does not have a nonempty population with ID '%s'.\n", pheno_name, first_token_start);
            goto FstReport_ret_INCONSISTENT_INPUT;
          }
          *second_token_end = '\0';
          const uintptr_t cat_idx2 = 1 + bsearch_strptr_natural(second_token_start, nonnull_category_names, nonnull_category_ct);
          if (unlikely(!IsSet(cats_seen, cat_idx2))) {
            logerrprintfww("Error: --fst phenotype '%s' does not have a nonempty population with ID '%s'.\n", pheno_name, second_token_start);
            goto FstReport_ret_INCONSISTENT_INPUT;
          }
          uintptr_t bit_idx;
          if (cat_idx1 < cat_idx2) {
            bit_idx = ((cat_idx2 * (cat_idx2 - 1)) / 2) + cat_idx1;
          } else if (likely(cat_idx2 < cat_idx1)) {
            bit_idx = ((cat_idx1 * (cat_idx1 - 1)) / 2) + cat_idx2;
          } else {
            logerrprintfww("Error: Both populations on line %" PRIuPTR " of %s are the same.\n", line_idx, in_fname);
            goto FstReport_ret_MALFORMED_INPUT;
          }
          SetBit(cat_idx1, pop_include);
          SetBit(cat_idx2, pop_include);
          duplicate_ct += IsSet(selected_pairs, bit_idx);
          SetBit(bit_idx, selected_pairs);
        }
        if (unlikely(CleanupTextStream2("--fst file", &txs, &reterr))) {
          goto FstReport_ret_1;
        }
        const uintptr_t named_pop_ct = PopcountWords(pop_include, raw_cat_ctl);
        if (named_pop_ct != pop_ct) {
          pop_ct = named_pop_ct;
          memcpy(cats_seen, pop_include, raw_cat_ctl * sizeof(intptr_t));
          sample_ct = RemoveExcludedCats(old_cats, cats_seen, raw_sample_ct, sample_ct, sample_include);
        }
        BigstackEndReset(bigstack_end_mark2);
        pop_pair_ct = PopcountWords(selected_pairs, possible_pair_ctl_u64);
        if (unlikely(!pop_pair_ct)) {
          logerrprintfww("Error: %s is empty.\n", in_fname);
          goto FstReport_ret_INCONSISTENT_INPUT;
        }
        logprintfww("--fst: %" PRIuPTR "%s population pair%s loaded from %s.\n", pop_pair_ct, duplicate_ct? " distinct" : "", (pop_pair_ct == 1)? "" : "s", in_fname);
        if (duplicate_ct) {
          logerrprintf("Warning: %" PRIuPTR " duplicate pair%s in --fst file.\n", duplicate_ct, (duplicate_ct == 1)? "" : "s");
        }
        if (unlikely(bigstack_alloc_u32(2 * pop_pair_ct, &pop_pairs))) {
          goto FstReport_ret_NOMEM;
        }
        uint32_t* pop_pairs_iter = pop_pairs;
        uintptr_t selected_pairs_base = 0;
        uintptr_t selected_pairs_bits = selected_pairs[0];
        uintptr_t cur_tri = 0;
        uintptr_t idx_hi = 1;
        for (uintptr_t ulii = 0; ulii != pop_pair_ct; ++ulii) {
          const uintptr_t selected_pairs_uidx = BitIter1(selected_pairs, &selected_pairs_base, &selected_pairs_bits);
          uintptr_t idx_lo = selected_pairs_uidx - cur_tri;
          while (idx_lo >= idx_hi) {
            idx_lo -= idx_hi;
            cur_tri += idx_hi;
            ++idx_hi;
          }
          *pop_pairs_iter++ = idx_lo;
          *pop_pairs_iter++ = idx_hi;
        }
      }
      uint32_t* sample_include_cumulative_popcounts;
      uint32_t* sample_to_pop_idx;
      uint32_t* diploid_pop_sizes;
      uint32_t* old_cat_idx_to_new;
      if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                   bigstack_alloc_u32(sample_ct, &sample_to_pop_idx) ||
                   bigstack_calloc_u32(pop_ct, &diploid_pop_sizes) ||
                   bigstack_alloc_kcp(pop_ct, &pop_names) ||
                   bigstack_end_alloc_u32(1 + nonnull_category_ct, &old_cat_idx_to_new))) {
        goto FstReport_ret_NOMEM;
      }
      ctx.sample_include = sample_include;
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
      uintptr_t cat_uidx_base = 0;
      uintptr_t cats_seen_bits = cats_seen[0];
      for (uintptr_t cat_idx = 0; cat_idx != pop_ct; ++cat_idx) {
        const uintptr_t cat_uidx = BitIter1(cats_seen, &cat_uidx_base, &cats_seen_bits);
        old_cat_idx_to_new[cat_uidx] = cat_idx;
        pop_names[cat_idx] = pheno_col->category_names[cat_uidx];
      }
      uintptr_t sample_uidx_base = 0;
      uintptr_t sample_include_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
        const uint32_t pop_idx = old_cat_idx_to_new[old_cats[sample_uidx]];
        sample_to_pop_idx[sample_idx] = pop_idx;
        diploid_pop_sizes[pop_idx] += 1;
      }
      ctx.sample_to_pop_idx = sample_to_pop_idx;
      ctx.diploid_pop_sizes = diploid_pop_sizes;
      ctx.sample_ct = sample_ct;
      ctx.pop_ct = pop_ct;
      if (!(flags & (kfFstOneBasePop | kfFstExplicitPopIds | kfFstPopPairFile))) {
        const uint64_t pop_pair_ct_u64 = (S_CAST(uint64_t, pop_ct) * (pop_ct - k1LU)) / 2;
        if (unlikely(bigstack_alloc64_u32(2 * pop_pair_ct_u64, &pop_pairs))) {
          goto FstReport_ret_NOMEM;
        }
        pop_pair_ct = pop_pair_ct_u64;
        uint32_t* pop_pairs_iter = pop_pairs;
        const uint32_t pop_ct_m1 = pop_ct - 1;
        for (uint32_t first_pop_idx = 0; first_pop_idx != pop_ct_m1; ++first_pop_idx) {
          for (uint32_t second_pop_idx = first_pop_idx + 1; second_pop_idx != pop_ct; ++second_pop_idx) {
            *pop_pairs_iter++ = first_pop_idx;
            *pop_pairs_iter++ = second_pop_idx;
          }
        }
      } else {
        const uintptr_t pop_pair_ct_x2 = pop_pair_ct * 2;
        for (uintptr_t ulii = 0; ulii != pop_pair_ct_x2; ++ulii) {
          pop_pairs[ulii] = old_cat_idx_to_new[pop_pairs[ulii]];
        }
      }
      BigstackEndReset(bigstack_end_mark);
    }
    logprintf("--fst: Analyzing %u samples across %" PRIuPTR " populations.\n", sample_ct, pop_ct);
    if (unlikely(bigstack_alloc_wp(max_thread_ct, &ctx.raregenos) ||
                 bigstack_alloc_u32p(max_thread_ct, &ctx.difflist_sample_id_bufs) ||
                 bigstack_alloc_u32p(max_thread_ct, &ctx.pop_geno_bufs))) {
      goto FstReport_ret_NOMEM;
    }
    const uint32_t report_variants = (flags / kfFstReportVariants) & 1;
    if (report_variants) {
      max_variant_file_ct = MINV(pop_pair_ct, kFstReportVariantsBatchMax);
      if (unlikely(bigstack_calloc_cp(max_variant_file_ct, &v_cswritep_arr) ||
                   BIGSTACK_ALLOC_X(CompressStreamState, max_variant_file_ct, &v_css_arr))) {
        goto FstReport_ret_NOMEM;
      }
      for (uintptr_t fidx = 0; fidx != max_variant_file_ct; ++fidx) {
        PreinitCstream(&(v_css_arr[fidx]));
      }
    }
    const uint32_t v_output_zst = (flags / kfFstZs) & 1;
    const uint32_t chrom_vcol = (flags / kfFstVcolChrom) & 1;
    char* chr_buf = nullptr;
    if (chrom_vcol) {
      uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto FstReport_ret_NOMEM;
      }
    }
    const uint32_t pos_vcol = (flags / kfFstVcolPos) & 1;
    const uint32_t ref_vcol = (flags / kfFstVcolRef) & 1;
    const uint32_t alt_vcol = (flags / kfFstVcolAlt) & 1;
    const uint32_t alleles_needed = ref_vcol || alt_vcol;
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t all_nonref = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t nobs_vcol = (flags / kfFstVcolNobs) & 1;
    const uint32_t nallele_vcol = (flags / kfFstVcolNallele) & 1;
    const uint32_t fstfrac_vcol = (flags / kfFstVcolFstfrac) & 1;
    const uint32_t fst_vcol = (flags / kfFstVcolFst) & 1;
    const uint32_t jackknife_blocksize = fst_infop->blocksize;
    // todo: tune this threshold
    const uint32_t max_difflist_len = sample_ct / 32;
    ctx.max_difflist_len = max_difflist_len;
    const uintptr_t raregeno_vec_ct = DivUp(max_difflist_len, kNypsPerVec);
    const uintptr_t difflist_sample_id_vec_ct = DivUp(max_difflist_len, kInt32PerVec);
    const uint32_t mhc_needed = (max_allele_ct > 2);

    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    unsigned char* bigstack_mark2 = g_bigstack_base;
    uint32_t allele_ct = 2;
    for (uint32_t is_x = 0; is_x != 2; ++is_x) {
      const uintptr_t* cur_variant_include = orig_variant_include;
      char* outname_end2 = outname_end;
      uint32_t cur_variant_ct;
      if (!is_x) {
        cur_variant_ct = orig_variant_ct;
        reterr = ConditionalAllocateNonAutosomalVariants(cip, nullptr, raw_variant_ct, &cur_variant_include, &cur_variant_ct);
        if (unlikely(reterr)) {
          goto FstReport_ret_1;
        }
        ctx.sex_male_collapsed = nullptr;
        ctx.haploid_pop_sizes = nullptr;
      } else {
        if (is_wc) {
          break;
        }
        outname_end2 = strcpya_k(outname_end2, ".x");
        cur_variant_ct = 0;
        uint32_t x_code;
        if (XymtExists(cip, kChrOffsetX, &x_code)) {
          const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[x_code];
          const uint32_t start_idx = cip->chr_fo_vidx_start[chr_fo_idx];
          const uint32_t end_idx = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          cur_variant_ct = PopcountBitRange(orig_variant_include, start_idx, end_idx);
          if (cur_variant_ct && (cur_variant_ct != orig_variant_ct)) {
            // may want this in a plink2_common function
            uintptr_t* tmp_variant_include;
            if (unlikely(bigstack_alloc_w(raw_variant_ctl, &tmp_variant_include))) {
              goto FstReport_ret_NOMEM;
            }
            const uint32_t start_widx = start_idx / kBitsPerWord;
            const uint32_t end_widx = DivUp(end_idx, kBitsPerWord);
            ZeroWArr(start_widx, tmp_variant_include);
            memcpy(&(tmp_variant_include[start_widx]), &(orig_variant_include[start_widx]), (end_widx - start_widx) * sizeof(intptr_t));
            ZeroWArr(raw_variant_ctl - end_widx, &(tmp_variant_include[end_widx]));
            const uint32_t start_remainder = start_idx % kBitsPerWord;
            tmp_variant_include[start_widx] &= (~k0LU) << start_remainder;
            const uint32_t end_remainder = end_idx % kBitsPerWord;
            if (end_remainder) {
              tmp_variant_include[end_widx - 1] &= (k1LU << end_remainder) - 1;
            }
            cur_variant_include = tmp_variant_include;
          }
          const uint32_t haploid_present = !IntersectionIsEmpty(sample_include, sex_male, raw_sample_ctl);
          if (haploid_present) {
            const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
            uintptr_t* sex_male_collapsed;
            uint32_t* haploid_pop_sizes;
            if (unlikely(bigstack_alloc_w(sample_ctl, &sex_male_collapsed) ||
                         bigstack_calloc_u32(pop_ct, &haploid_pop_sizes))) {
              goto FstReport_ret_NOMEM;
            }
            CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_male_collapsed);
            const uint32_t* sample_to_pop_idx = ctx.sample_to_pop_idx;
            for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
              const uint32_t pop_idx = sample_to_pop_idx[sample_idx];
              haploid_pop_sizes[pop_idx] += IsSet(sex_male_collapsed, sample_idx);
            }
            // This assumes is_x happens after !is_x.
            for (uintptr_t pop_idx = 0; pop_idx != pop_ct; ++pop_idx) {
              ctx.diploid_pop_sizes[pop_idx] -= haploid_pop_sizes[pop_idx];
            }
            ctx.sex_male_collapsed = sex_male_collapsed;
            ctx.haploid_pop_sizes = haploid_pop_sizes;
          }
        }
      }
      if (!cur_variant_ct) {
        continue;
      }
      const uint32_t provref_vcol = ref_vcol && ProvrefCol(cur_variant_include, nonref_flags, flags / kfFstVcolMaybeprovref, raw_variant_ct, all_nonref);
      ctx.cur_variant_include = cur_variant_include;
      const uintptr_t pop_geno_vec_ct = DivUp(pop_ct * (4 << (ctx.sex_male_collapsed != nullptr)), kInt32PerVec);
      const uintptr_t thread_xalloc_vec_ct = raregeno_vec_ct + difflist_sample_id_vec_ct + pop_geno_vec_ct;
      const uintptr_t thread_xalloc_cacheline_ct = DivUp(thread_xalloc_vec_ct, kVecsPerCacheline);
      const uintptr_t per_allele_xalloc_byte_ct = pop_ct * sizeof(int32_t) * 2 * (1 + is_wc);
      uintptr_t per_variant_xalloc_byte_ct = per_allele_xalloc_byte_ct;
      if (nobs_vcol || is_wc) {
        per_variant_xalloc_byte_ct += pop_ct * sizeof(int32_t) * 2;
      } else {
        ctx.pop_nm_sample_cts[0] = nullptr;
        ctx.pop_nm_sample_cts[1] = nullptr;
      }
      if (!is_wc) {
        // half_within
        per_variant_xalloc_byte_ct += pop_ct * sizeof(double) * 2;

        ctx.pop_allele_het_cts[0] = nullptr;
        ctx.pop_allele_het_cts[1] = nullptr;
      } else {
        // pop_allele_het_cts already included in (1 + is_wc) term above
        ctx.half_within[0] = nullptr;
        ctx.half_within[1] = nullptr;
      }
      uint32_t calc_thread_ct = max_thread_ct;

      strcpy_k(outname_end2, ".fst.summary");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &s_outfile))) {
        goto FstReport_ret_OPEN_FAIL;
      }
      char* s_write_iter = g_textbuf;
      char* s_textbuf_flush = &(s_write_iter[kMaxMediumLine]);
      s_write_iter = strcpya_k(s_write_iter, "#POP1\tPOP2\t");
      const uint32_t nobs_scol = (flags / kfFstColNobs) & 1;
      if (nobs_scol) {
        s_write_iter = strcpya_k(s_write_iter, "OBS_CT\t");
      }
      if (!is_wc) {
        s_write_iter = strcpya_k(s_write_iter, "HUDSON_FST");
      } else {
        s_write_iter = strcpya_k(s_write_iter, "WC_FST");
      }
      if (jackknife_blocksize) {
        s_write_iter = strcpya_k(s_write_iter, "\tSE");
      }
      AppendBinaryEoln(&s_write_iter);

      uintptr_t pop_pair_batch_size = pop_pair_ct;
      uintptr_t pass_ct = 1;
      if (report_variants && (pop_pair_ct > kFstReportVariantsBatchMax)) {
        pop_pair_batch_size = kFstReportVariantsBatchMax;
        pass_ct = DivUp(pop_pair_ct, kFstReportVariantsBatchMax);
      }
      double* fst_numer_sums;
      double* fst_denom_sums;
      uint32_t* fst_nobs;
      if (unlikely(bigstack_alloc_d(pop_pair_batch_size, &fst_numer_sums) ||
                   bigstack_alloc_d(pop_pair_batch_size, &fst_denom_sums) ||
                   bigstack_alloc_u32(pop_pair_batch_size, &fst_nobs))) {
        goto FstReport_ret_NOMEM;
      }
      uint32_t jackknife_block_max = 0;
      if (jackknife_blocksize) {
        if (cur_variant_ct <= jackknife_blocksize) {
          logerrprintf("Warning: Too few %s variants for --fst blocksize=%u to be useful.\n", is_x? "chrX" : "autosomal", jackknife_blocksize);
        }
        jackknife_block_max = DivUp(cur_variant_ct, jackknife_blocksize);
      }
      // <output prefix, including .x>.<popID1>.<popID2>.fst.var[.zst]
      const uint32_t pop_name_capacity = kPglFnamesize - 11 - (4 * v_output_zst) - S_CAST(uintptr_t, outname_end2 - outname);
      unsigned char* bigstack_mark3 = g_bigstack_base;
      for (uintptr_t pass_idx = 0; pass_idx != pass_ct; ++pass_idx) {
        const uintptr_t pop_pair_start_idx = pass_idx * pop_pair_batch_size;
        if (pass_idx == pass_ct - 1) {
          pop_pair_batch_size = pop_pair_ct - pop_pair_start_idx;
        }
        const uint32_t* cur_pop_pairs = &(pop_pairs[2 * pop_pair_start_idx]);
        if (report_variants) {
          if (pop_pair_batch_size > pop_pair_ct - pop_pair_start_idx) {
            pop_pair_batch_size = pop_pair_ct - pop_pair_start_idx;
          }
          for (uintptr_t fidx = 0; fidx != pop_pair_batch_size; ++fidx) {
            {
              const uint32_t pop_idx1 = cur_pop_pairs[2 * fidx];
              const uint32_t pop_idx2 = cur_pop_pairs[2 * fidx + 1];
              char* fname_iter = outname_end2;
              const uint32_t slen1 = strlen(pop_names[pop_idx1]);
              const uint32_t slen2 = strlen(pop_names[pop_idx2]);
              if (unlikely(slen1 + slen2 > pop_name_capacity)) {
                logerrputs("Error: Population name and/or --out argument too long.\n");
                goto FstReport_ret_INCONSISTENT_INPUT;
              }
              *fname_iter++ = '.';
              fname_iter = memcpyax(fname_iter, pop_names[pop_idx1], slen1, '.');
              fname_iter = memcpya(fname_iter, pop_names[pop_idx2], slen2);
              OutnameZstSet(".fst.var", v_output_zst, fname_iter);
            }
            reterr = InitCstreamAlloc(outname, 0, v_output_zst, 1, 2 * kCompressStreamBlock, &(v_css_arr[fidx]), &(v_cswritep_arr[fidx]));
            if (unlikely(reterr)) {
              goto FstReport_ret_1;
            }
            char* cswritep = v_cswritep_arr[fidx];
            *cswritep++ = '#';
            if (chrom_vcol) {
              cswritep = strcpya_k(cswritep, "CHROM\t");
            }
            if (pos_vcol) {
              cswritep = strcpya_k(cswritep, "POS\t");
            }
            cswritep = strcpya_k(cswritep, "ID");
            if (ref_vcol) {
              cswritep = strcpya_k(cswritep, "\tREF");
            }
            if (alt_vcol) {
              cswritep = strcpya_k(cswritep, "\tALT");
            }
            if (provref_vcol) {
              cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
            }
            if (nobs_vcol) {
              cswritep = strcpya_k(cswritep, "\tOBS_CT");
            }
            if (nallele_vcol) {
              cswritep = strcpya_k(cswritep, "\tPOP1_ALLELE_CT\tPOP2_ALLELE_CT");
            }
            if (fstfrac_vcol) {
              cswritep = strcpya_k(cswritep, "\tFST_NUMER\tFST_DENOM");
            }
            if (fst_vcol) {
              if (!is_wc) {
                cswritep = strcpya_k(cswritep, "\tHUDSON_FST");
              } else {
                cswritep = strcpya_k(cswritep, "\tWC_FST");
              }
            }
            AppendBinaryEoln(&cswritep);
            v_cswritep_arr[fidx] = cswritep;
          }
        }
        ZeroDArr(pop_pair_batch_size, fst_numer_sums);
        ZeroDArr(pop_pair_batch_size, fst_denom_sums);
        ZeroU32Arr(pop_pair_batch_size, fst_nobs);
        double* jackknife_fst_numers = nullptr;
        double* jackknife_fst_denoms = nullptr;
        uint32_t* jackknife_next_block_starts = nullptr;
        // In the jackknife case, fst_{numer,denom}_sums are only for the
        // current block while the calculation is in progress; we save them off
        // to jackknife_fst_{numers,denoms} when each block is complete.
        if (jackknife_blocksize) {
          // leading dimension = pop-pair index
          // strictly speaking, pop_pair_batch_size should be reduced if this
          // would be too large, but that should never come up with any sane
          // blocksize= choice
          if (unlikely(bigstack_alloc_d(jackknife_block_max * pop_pair_batch_size, &jackknife_fst_numers) ||
                       bigstack_alloc_d(jackknife_block_max * pop_pair_batch_size, &jackknife_fst_denoms) ||
                       bigstack_alloc_u32(pop_pair_batch_size, &jackknife_next_block_starts))) {
            goto FstReport_ret_NOMEM;
          }
          for (uintptr_t ulii = 0; ulii != pop_pair_batch_size; ++ulii) {
            jackknife_next_block_starts[ulii] = jackknife_blocksize;
          }
        }

        uintptr_t bytes_avail = bigstack_left();
        // defend against adverse rounding
        if (unlikely(bytes_avail < 6 * kCacheline)) {
          goto FstReport_ret_NOMEM;
        }
        bytes_avail -= 6 * kCacheline;

        STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
        ctx.thread_read_mhc = nullptr;
        uint32_t read_block_size;
        uintptr_t max_alt_allele_block_size;
        if (unlikely(PgenMtLoadInit(cur_variant_include, sample_ct, cur_variant_ct, bytes_avail, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, per_allele_xalloc_byte_ct, pgfip, &calc_thread_ct, &ctx.genovecs, mhc_needed? (&ctx.thread_read_mhc) : nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, &max_alt_allele_block_size, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
          goto FstReport_ret_NOMEM;
        }
        if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
          goto FstReport_ret_NOMEM;
        }
        ctx.err_info = (~0LLU) << 32;
        for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
          unsigned char* cur_alloc = S_CAST(unsigned char*, bigstack_alloc_raw(thread_xalloc_cacheline_ct * kCacheline));
          ctx.raregenos[tidx] = R_CAST(uintptr_t*, cur_alloc);
          cur_alloc = &(cur_alloc[raregeno_vec_ct * kBytesPerVec]);
          ctx.difflist_sample_id_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
          cur_alloc = &(cur_alloc[difflist_sample_id_vec_ct * kBytesPerVec]);
          ctx.pop_geno_bufs[tidx] = R_CAST(uint32_t*, cur_alloc);
          cur_alloc = &(cur_alloc[pop_geno_vec_ct * kBytesPerVec]);
        }
        if (nobs_vcol || is_wc) {
          const uintptr_t geno_obs_alloc = RoundUpPow2(pop_ct * sizeof(int32_t) * read_block_size, kCacheline);
          ctx.pop_nm_sample_cts[0] = S_CAST(uint32_t*, bigstack_alloc_raw(geno_obs_alloc));
          ctx.pop_nm_sample_cts[1] = S_CAST(uint32_t*, bigstack_alloc_raw(geno_obs_alloc));
        }
        const uintptr_t max_allele_block_size = read_block_size + max_alt_allele_block_size;
        const uintptr_t allele_pop_alloc = RoundUpPow2(pop_ct * max_allele_block_size * sizeof(int32_t), kCacheline);
        ctx.pop_allele_obs_cts[0] = S_CAST(uint32_t*, bigstack_alloc_raw(allele_pop_alloc));
        ctx.pop_allele_obs_cts[1] = S_CAST(uint32_t*, bigstack_alloc_raw(allele_pop_alloc));
        if (!is_wc) {
          const uintptr_t half_within_alloc = RoundUpPow2(pop_ct * sizeof(double) * read_block_size, kCacheline);
          ctx.half_within[0] = S_CAST(double*, bigstack_alloc_raw(half_within_alloc));
          ctx.half_within[1] = S_CAST(double*, bigstack_alloc_raw(half_within_alloc));
        } else {
          ctx.pop_allele_het_cts[0] = S_CAST(uint32_t*, bigstack_alloc_raw(allele_pop_alloc));
          ctx.pop_allele_het_cts[1] = S_CAST(uint32_t*, bigstack_alloc_raw(allele_pop_alloc));
        }
        SetThreadFuncAndData(FstThread, &ctx, &tg);

        if (pop_pair_ct == pop_pair_batch_size) {
          logprintf("%s --fst: ", is_x? "chrX" : "Autosomal");
        } else {
          logprintf("%s --fst pass %u/%u: ", is_x? "chrX" : "Autosomal", pass_idx + 1, pass_ct);
        }
        fputs("0%", stdout);
        fflush(stdout);
        uint32_t pct = 0;

        uintptr_t write_variant_uidx_base = 0;
        uintptr_t write_variant_bits = cur_variant_include[0];
        uint32_t parity = 0;
        uint32_t read_block_idx = 0;
        uint32_t chr_fo_idx = UINT32_MAX;
        uint32_t chr_end = 0;
        uint32_t chr_buf_blen = 0;
        uint32_t prev_block_variant_ct = 0;
        uint32_t next_print_variant_idx = (cur_variant_ct + 99) / 100;
        for (uint32_t variant_idx = 0; ; ) {
          const uint32_t cur_block_variant_ct = MultireadNonempty(cur_variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
          if (unlikely(reterr)) {
            goto FstReport_ret_PGR_FAIL;
          }
          if (variant_idx) {
            JoinThreads(&tg);
            reterr = S_CAST(PglErr, ctx.err_info);
            if (unlikely(reterr)) {
              PgenErrPrintNV(reterr, ctx.err_info >> 32);
              goto FstReport_ret_1;
            }
          }
          if (!IsLastBlock(&tg)) {
            ctx.cur_block_variant_ct = cur_block_variant_ct;
            ComputeUidxStartPartition(cur_variant_include, cur_block_variant_ct, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
            PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
            if (variant_idx + cur_block_variant_ct == cur_variant_ct) {
              DeclareLastThreadBlock(&tg);
            }
            if (unlikely(SpawnThreads(&tg))) {
              goto FstReport_ret_THREAD_CREATE_FAIL;
            }
          }

          parity = 1 - parity;
          if (variant_idx) {
            // process *previous* block results
            const uint32_t* pop_nm_sample_cts_iter = ctx.pop_nm_sample_cts[parity];
            const uint32_t* pop_allele_obs_cts_iter = ctx.pop_allele_obs_cts[parity];
            const double* half_within_iter = ctx.half_within[parity];
            const uint32_t* pop_allele_het_cts_iter = ctx.pop_allele_het_cts[parity];
            for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx) {
              const uint32_t write_variant_uidx = BitIter1(cur_variant_include, &write_variant_uidx_base, &write_variant_bits);
              if (chr_buf && (write_variant_uidx >= chr_end)) {
                do {
                  ++chr_fo_idx;
                  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
                } while (write_variant_uidx >= chr_end);
                const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
                char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
                *chr_name_end = '\t';
                chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
              }
              uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
              if (allele_idx_offsets) {
                allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
                allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offset_base;
              }
              for (uintptr_t pop_pair_bidx = 0; pop_pair_bidx != pop_pair_batch_size; ++pop_pair_bidx) {
                const uint32_t pop_idx1 = cur_pop_pairs[2 * pop_pair_bidx];
                const uint32_t pop_idx2 = cur_pop_pairs[2 * pop_pair_bidx + 1];
                // See allel/stats/fst.py under
                // https://github.com/cggh/scikit-allel for a more readable
                // form of these computations.
                double fst_numer = 0.0;
                double fst_denom = 0.0;
                uintptr_t pop1_n_hap = 0;
                uintptr_t pop2_n_hap = 0;
                if (!is_wc) {
                  // Hudson
                  uint64_t n_same = 0;
                  for (uintptr_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    const uintptr_t ct1 = pop_allele_obs_cts_iter[allele_idx * pop_ct + pop_idx1];
                    const uintptr_t ct2 = pop_allele_obs_cts_iter[allele_idx * pop_ct + pop_idx2];
                    pop1_n_hap += ct1;
                    pop2_n_hap += ct2;
                    n_same += S_CAST(uint64_t, ct1) * ct2;
                  }
                  const uint64_t n_pairs = S_CAST(uint64_t, pop1_n_hap) * pop2_n_hap;
                  const uint64_t n_diff = n_pairs - n_same;
                  if (n_diff) {
                    const double within = half_within_iter[pop_idx1] + half_within_iter[pop_idx2];
                    fst_denom = u63tod(n_diff) / u63tod(n_pairs);
                    fst_numer = fst_denom - within;
                  }
                } else {
                  // Weir-Cockerham
                  const uintptr_t pop1_size = pop_nm_sample_cts_iter[pop_idx1];
                  const uintptr_t pop2_size = pop_nm_sample_cts_iter[pop_idx2];
                  pop1_n_hap = pop1_size * 2;
                  pop2_n_hap = pop2_size * 2;
                  const uintptr_t n_total = pop1_size + pop2_size;
                  const double n_total_d = u31tod(n_total);
                  const double n_total_recip = 1.0 / n_total_d;
                  const double n_bar = n_total_d * 0.5;
                  const double n_bar_m1_recip = 1.0 / (n_bar - 1);
                  const double n_bar_div_n_c = n_bar / (n_total_d - u63tod(pop1_size * pop1_size + pop2_size * pop2_size) * n_total_recip);
                  const double pop1_size_d = u31tod(pop1_size);
                  const double pop2_size_d = u31tod(pop2_size);
                  const double pop1_size_recip = 1.0 / pop1_size_d;
                  const double pop2_size_recip = 1.0 / pop2_size_d;
                  double a_sum = 0.0;
                  double b_sum = 0.0;
                  double c_sum = 0.0;
                  const uintptr_t wc_allele_ct = (allele_ct == 2)? 1 : allele_ct;
                  for (uintptr_t allele_idx = 0; allele_idx != wc_allele_ct; ++allele_idx) {
                    const uint32_t pop1_cur_allele_ct = pop_allele_obs_cts_iter[allele_idx * pop_ct + pop_idx1];
                    const uint32_t pop2_cur_allele_ct = pop_allele_obs_cts_iter[allele_idx * pop_ct + pop_idx2];
                    const uint32_t total_allele_ct = pop1_cur_allele_ct + pop2_cur_allele_ct;
                    if (!total_allele_ct) {
                      continue;
                    }
                    if (total_allele_ct == 2 * n_total) {
                      // guaranteed to have no variation; don't want to deal
                      // with rounding errors, etc.
                      break;
                    }
                    const double pop1_cur_allele_ct_d = u31tod(pop1_cur_allele_ct);
                    const double pop2_cur_allele_ct_d = u31tod(pop2_cur_allele_ct);
                    const double pop1_afreq = pop1_cur_allele_ct_d * pop1_size_recip * 0.5;
                    const double pop2_afreq = pop2_cur_allele_ct_d * pop2_size_recip * 0.5;
                    const double p_bar = (pop1_cur_allele_ct_d + pop2_cur_allele_ct_d) * 0.5 * n_total_recip;
                    const double s_squared_pop1_term = pop1_afreq - p_bar;
                    const double s_squared_pop2_term = pop2_afreq - p_bar;
                    const double s_squared = (pop1_size_d * s_squared_pop1_term * s_squared_pop1_term + pop2_size_d * s_squared_pop2_term * s_squared_pop2_term) * n_total_recip * 2;
                    const double h_bar = u31tod(pop_allele_het_cts_iter[allele_idx * pop_ct + pop_idx1] + pop_allele_het_cts_iter[allele_idx * pop_ct + pop_idx2]) * n_total_recip;
                    const double p_bar_times_1_minus_p_bar = p_bar * (1 - p_bar);
                    a_sum += n_bar_div_n_c * (s_squared - (p_bar_times_1_minus_p_bar - 0.5 * s_squared - 0.25 * h_bar) * n_bar_m1_recip);
                    b_sum += n_bar * n_bar_m1_recip * (p_bar_times_1_minus_p_bar - 0.5 * s_squared - (0.5 - 0.5 * n_total_recip) * h_bar);
                    c_sum += h_bar * 0.5;
                  }
                  fst_numer = a_sum;
                  fst_denom = a_sum + b_sum + c_sum;
                }
                if ((fst_denom != 0.0) && (fst_numer == fst_numer)) {
                  fst_numer_sums[pop_pair_bidx] += fst_numer;
                  fst_denom_sums[pop_pair_bidx] += fst_denom;
                  const uint32_t new_nobs = fst_nobs[pop_pair_bidx] + 1;
                  fst_nobs[pop_pair_bidx] = new_nobs;
                  if (jackknife_next_block_starts && (new_nobs == jackknife_next_block_starts[pop_pair_bidx])) {
                    const uintptr_t idx_2d = pop_pair_bidx * jackknife_block_max + (new_nobs / jackknife_blocksize) - 1;
                    const double cur_block_fst_numer_sum = fst_numer_sums[pop_pair_bidx];
                    const double cur_block_fst_denom_sum = fst_denom_sums[pop_pair_bidx];
                    jackknife_fst_numers[idx_2d] = cur_block_fst_numer_sum;
                    jackknife_fst_denoms[idx_2d] = cur_block_fst_denom_sum;
                    fst_numer_sums[pop_pair_bidx] = 0.0;
                    fst_denom_sums[pop_pair_bidx] = 0.0;
                    jackknife_next_block_starts[pop_pair_bidx] = new_nobs + jackknife_blocksize;
                  }
                }
                if (!report_variants) {
                  continue;
                }
                CompressStreamState* cssp = &(v_css_arr[pop_pair_bidx]);
                char* cswritep = v_cswritep_arr[pop_pair_bidx];
                // first several columns are practically identical to e.g.
                // AppendSdiffPregenoFields(); probably want a library
                // function for this
                if (chrom_vcol) {
                  cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
                }
                if (pos_vcol) {
                  cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
                }
                cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
                if (alleles_needed) {
                  const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
                  if (ref_vcol) {
                    *cswritep++ = '\t';
                    const char* cur_allele = cur_alleles[0];
                    const uint32_t allele_slen = strlen(cur_allele);
                    if (unlikely(CsputsStd(cur_allele, allele_slen, cssp, &cswritep))) {
                      v_cswritep_arr[pop_pair_bidx] = cswritep;
                      goto FstReport_ret_WRITE_FAIL;
                    }
                  }
                  if (alt_vcol) {
                    *cswritep++ = '\t';
                    for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
                      const char* cur_allele = cur_alleles[allele_idx];
                      const uint32_t allele_slen = strlen(cur_allele);
                      if (unlikely(CsputsStd(cur_allele, allele_slen, cssp, &cswritep))) {
                        v_cswritep_arr[pop_pair_bidx] = cswritep;
                        goto FstReport_ret_WRITE_FAIL;
                      }
                      *cswritep++ = ',';
                    }
                    --cswritep;
                  }
                  if (provref_vcol) {
                    *cswritep++ = '\t';
                    *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, write_variant_uidx)))? 'Y' : 'N';
                  }
                }
                if (nobs_vcol) {
                  *cswritep++ = '\t';
                  cswritep = u32toa(pop_nm_sample_cts_iter[pop_idx1] + pop_nm_sample_cts_iter[pop_idx2], cswritep);
                }
                if (nallele_vcol) {
                  *cswritep++ = '\t';
                  cswritep = u32toa_x(pop1_n_hap, '\t', cswritep);
                  cswritep = u32toa(pop2_n_hap, cswritep);
                }
                if (fstfrac_vcol) {
                  *cswritep++ = '\t';
                  cswritep = dtoa_g(fst_numer, cswritep);
                  *cswritep++ = '\t';
                  cswritep = dtoa_g(fst_denom, cswritep);
                }
                if (fst_vcol) {
                  *cswritep++ = '\t';
                  cswritep = dtoa_g(fst_numer / fst_denom, cswritep);
                }
                AppendBinaryEoln(&cswritep);
                if (unlikely(Cswrite(cssp, &cswritep))) {
                  goto FstReport_ret_WRITE_FAIL;
                }
                v_cswritep_arr[pop_pair_bidx] = cswritep;
              }
              if (pop_nm_sample_cts_iter) {
                pop_nm_sample_cts_iter = &(pop_nm_sample_cts_iter[pop_ct]);
              }
              pop_allele_obs_cts_iter = &(pop_allele_obs_cts_iter[pop_ct * allele_ct]);
              if (half_within_iter) {
                half_within_iter = &(half_within_iter[pop_ct]);
              } else {
                pop_allele_het_cts_iter = &(pop_allele_het_cts_iter[pop_ct * allele_ct]);
              }
            }
          }
          if (variant_idx == cur_variant_ct) {
            break;
          }
          if (variant_idx >= next_print_variant_idx) {
            if (pct > 10) {
              putc_unlocked('\b', stdout);
            }
            pct = (variant_idx * 100LLU) / cur_variant_ct;
            printf("\b\b%u%%", pct++);
            fflush(stdout);
            next_print_variant_idx = (pct * S_CAST(uint64_t, cur_variant_ct) + 99) / 100;
          }
          ++read_block_idx;
          prev_block_variant_ct = cur_block_variant_ct;
          variant_idx += cur_block_variant_ct;
          pgfip->block_base = main_loadbufs[parity];
        }
        for (uintptr_t pop_pair_bidx = 0; pop_pair_bidx != pop_pair_batch_size; ++pop_pair_bidx) {
          const uint32_t pop_idx1 = cur_pop_pairs[2 * pop_pair_bidx];
          const uint32_t pop_idx2 = cur_pop_pairs[2 * pop_pair_bidx + 1];
          const uint32_t nobs = fst_nobs[pop_pair_bidx];
          double fst_numer_sum = fst_numer_sums[pop_pair_bidx];
          double fst_denom_sum = fst_denom_sums[pop_pair_bidx];
          double fst_se = 0.0;
          if (jackknife_fst_numers) {
            double* block_fst_numers = &(jackknife_fst_numers[pop_pair_bidx * jackknife_block_max]);
            double* block_fst_denoms = &(jackknife_fst_denoms[pop_pair_bidx * jackknife_block_max]);
            if (nobs > jackknife_blocksize) {
              const uint32_t n_block = DivUp(nobs, jackknife_blocksize);
              uint32_t prev_block_ct = n_block;
              uint32_t last_block_size = jackknife_blocksize;
              if (nobs != n_block * jackknife_blocksize) {
                --prev_block_ct;
                block_fst_numers[prev_block_ct] = fst_numer_sum;
                block_fst_denoms[prev_block_ct] = fst_denom_sum;
                last_block_size = nobs - (prev_block_ct * jackknife_blocksize);
              }
              for (uint32_t block_idx = 0; block_idx != prev_block_ct; ++block_idx) {
                fst_numer_sum += block_fst_numers[block_idx];
                fst_denom_sum += block_fst_denoms[block_idx];
              }
              const double nobs_d = u31tod(nobs);
              const double nobs_recip = 1.0 / nobs_d;
              // https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf
              const double theta_hat = fst_numer_sum / fst_denom_sum;
              // first pass: compute theta_jack
              double theta_jack = 0.0;
              double cur_block_size_d = u31tod(jackknife_blocksize);
              for (uint32_t block_idx = 0; ; ++block_idx) {
                if (block_idx >= prev_block_ct) {
                  if (block_idx > prev_block_ct) {
                    break;
                  }
                  cur_block_size_d = u31tod(last_block_size);
                }
                // we can avoid computing some of these values twice, but this
                // isn't a bottleneck
                const double theta_with_j_removed = (fst_numer_sum - block_fst_numers[block_idx]) / (fst_denom_sum - block_fst_denoms[block_idx]);
                theta_jack += (theta_hat - theta_with_j_removed) + cur_block_size_d * theta_with_j_removed * nobs_recip;
              }
              // second pass: compute variance estimate
              double main_sum = 0.0;
              double hh = nobs_d / u31tod(jackknife_blocksize);
              for (uint32_t block_idx = 0; ; ++block_idx) {
                if (block_idx >= prev_block_ct) {
                  if (block_idx > prev_block_ct) {
                    break;
                  }
                  hh = nobs_d / u31tod(last_block_size);
                }
                const double h_minus_1 = hh - 1.0;
                const double theta_with_j_removed = (fst_numer_sum - block_fst_numers[block_idx]) / (fst_denom_sum - block_fst_denoms[block_idx]);
                const double tau_j = hh * theta_hat - h_minus_1 * theta_with_j_removed;
                const double tau_j_minus_theta_jack = tau_j - theta_jack;
                main_sum += tau_j_minus_theta_jack * tau_j_minus_theta_jack / h_minus_1;
              }
              const double variance_estimate = main_sum / u31tod(n_block);
              fst_se = sqrt(variance_estimate);
            } else {
              if (nobs == jackknife_blocksize) {
                fst_numer_sum += block_fst_numers[0];
                fst_denom_sum += block_fst_denoms[0];
              }
              fst_se = 0.0 / 0.0;
            }
          }
          s_write_iter = strcpyax(s_write_iter, pop_names[pop_idx1], '\t');
          s_write_iter = strcpyax(s_write_iter, pop_names[pop_idx2], '\t');
          if (nobs_scol) {
            s_write_iter = u32toa_x(nobs, '\t', s_write_iter);
          }
          s_write_iter = dtoa_g(fst_numer_sum / fst_denom_sum, s_write_iter);
          if (jackknife_blocksize) {
            *s_write_iter++ = '\t';
            s_write_iter = dtoa_g(fst_se, s_write_iter);
          }
          AppendBinaryEoln(&s_write_iter);
          if (unlikely(fwrite_ck(s_textbuf_flush, s_outfile, &s_write_iter))) {
            goto FstReport_ret_WRITE_FAIL;
          }
        }
        if (report_variants) {
          for (uintptr_t fidx = 0; fidx != pop_pair_batch_size; ++fidx) {
            if (unlikely(CswriteCloseNull(&(v_css_arr[fidx]), v_cswritep_arr[fidx]))) {
              goto FstReport_ret_WRITE_FAIL;
            }
          }
        }
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        fputs("\b\b", stdout);
        logputs("done.\n");
        BigstackReset(bigstack_mark3);
      }
      if (report_variants) {
        logprintf("%s --fst: %u .fst.var%s file%s written.\n", is_x? "chrX" : "Autosomal", pop_pair_ct, v_output_zst? ".zst" : "", (pop_pair_ct == 1)? "" : "s");
        strcpy_k(outname_end2, ".fst.summary");
      }
      if (unlikely(fclose_flush_null(s_textbuf_flush, s_write_iter, &s_outfile))) {
        goto FstReport_ret_WRITE_FAIL;
      }
      logprintfww("%s --fst: Summary written to %s .\n", is_x? "chrX" : "Autosomal", outname);
      BigstackReset(bigstack_mark2);
    }
  }
  while (0) {
  FstReport_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  FstReport_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  FstReport_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  FstReport_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  FstReport_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  FstReport_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--fst file", &txs);
    break;
  FstReport_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  FstReport_ret_PHENO_NOT_FOUND:
    logerrprintfww("Error: --fst phenotype '%s' not loaded.\n", fst_infop->pheno_name);
  FstReport_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  FstReport_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 FstReport_ret_1:
  if (v_css_arr) {
    for (uintptr_t fidx = 0; fidx != max_variant_file_ct; ++fidx) {
      CswriteCloseCond(&(v_css_arr[fidx]), v_cswritep_arr[fidx]);
    }
  }
  CleanupThreads(&tg);
  fclose_cond(s_outfile);
  CleanupTextStream2("--fst file", &txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

typedef struct AlleleUniquenessCheckerStruct {
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const char* const* allele_storage;

  uint32_t* variant_uidx_starts;
  const char*** allele_sortbufs;

  uint32_t variant_ct;

  uint32_t dup_vidx; // UINT32_MAX if no variant found with duplicate alleles
} AlleleUniquenessChecker;

void CheckAlleleUniquenessMain(uint32_t tidx, uint32_t thread_ct, AlleleUniquenessChecker* ctx) {
  const uintptr_t* variant_include = ctx->variant_include;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const char* const* allele_storage = ctx->allele_storage;
  const char** allele_sortbuf = ctx->allele_sortbufs? ctx->allele_sortbufs[tidx]  : nullptr;
  const uint32_t variant_ct = ctx->variant_ct;
  const uint32_t variant_uidx_start = ctx->variant_uidx_starts[tidx];
  const uint32_t variant_idx_start = (variant_ct * S_CAST(uint64_t, tidx)) / thread_ct;
  const uint32_t variant_idx_end = (variant_ct * (S_CAST(uint64_t, tidx) + 1)) / thread_ct;
  uint32_t cur_allele_ct = 2;
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
    const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
    if (cur_allele_ct == 2) {
      if (unlikely(strequal_overread(cur_alleles[0], cur_alleles[1]))) {
        UpdateU32IfSmaller(variant_uidx, &(ctx->dup_vidx));
        return;
      }
    } else {
      memcpy(allele_sortbuf, cur_alleles, cur_allele_ct * sizeof(intptr_t));
      StrptrArrSortOverread(cur_allele_ct, allele_sortbuf);
      const char* prev_allele = allele_sortbuf[0];
      for (uint32_t aidx = 1; aidx != cur_allele_ct; ++aidx) {
        const char* cur_allele = allele_sortbuf[aidx];
        if (unlikely(strequal_overread(prev_allele, cur_allele))) {
          UpdateU32IfSmaller(variant_uidx, &(ctx->dup_vidx));
          return;
        }
        cur_allele = prev_allele;
      }
    }
  }
}

THREAD_FUNC_DECL CheckAlleleUniquenessThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  const uint32_t thread_ct = GetThreadCt(arg->sharedp) + 1;
  AlleleUniquenessChecker* ctx = S_CAST(AlleleUniquenessChecker*, arg->sharedp->context);
  CheckAlleleUniquenessMain(tidx, thread_ct, ctx);
  THREAD_RETURN;
}

PglErr CheckAlleleUniqueness(const uintptr_t* variant_include, const ChrInfo* cip, const ChrIdx* chr_idxs, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_thread_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  AlleleUniquenessChecker ctx;
  {
    const uint32_t thread_ct = ClipU32(variant_ct / 65536, 1, max_thread_ct);
    if (unlikely(SetThreadCt0(thread_ct - 1, &tg) ||
                 bigstack_alloc_u32(thread_ct, &ctx.variant_uidx_starts))) {
      goto CheckAlleleUniqueness_ret_NOMEM;
    }
    if (max_allele_ct == 2) {
      ctx.allele_sortbufs = nullptr;
    } else {
      if (unlikely(bigstack_alloc_kcpp(thread_ct, &ctx.allele_sortbufs))) {
        goto CheckAlleleUniqueness_ret_NOMEM;
      }
      for (uint32_t tidx = 0; tidx != thread_ct; ++tidx) {
        if (unlikely(bigstack_alloc_kcp(max_allele_ct, &(ctx.allele_sortbufs[tidx])))) {
          goto CheckAlleleUniqueness_ret_NOMEM;
        }
      }
    }
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.allele_storage = allele_storage;
    ctx.variant_ct = variant_ct;
    ctx.dup_vidx = UINT32_MAX;
    FillU32SubsetStarts(variant_include, thread_ct, 0, variant_ct, ctx.variant_uidx_starts);

    if (thread_ct > 1) {
      SetThreadFuncAndData(CheckAlleleUniquenessThread, &ctx, &tg);
      DeclareLastThreadBlock(&tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto CheckAlleleUniqueness_ret_THREAD_CREATE_FAIL;
      }
    }
    CheckAlleleUniquenessMain(thread_ct - 1, thread_ct, &ctx);
    JoinThreads0(&tg);
    const uint32_t dup_vidx = ctx.dup_vidx;
    if (unlikely(dup_vidx != UINT32_MAX)) {
      // Not a big deal, but best to make this still work when split
      // chromosomes are present.
      const uint32_t chr_idx = chr_idxs? chr_idxs[dup_vidx] : GetVariantChr(cip, dup_vidx);
      char* write_iter = strcpya_k(g_logbuf, "Error: Duplicate allele code in variant '");
      write_iter = strcpya(write_iter, variant_ids[dup_vidx]);
      write_iter = strcpya_k(write_iter, "' at position ");
      write_iter = chrtoa(cip, chr_idx, write_iter);
      *write_iter++ = ':';
      write_iter = u32toa(variant_bps[dup_vidx], write_iter);
      strcpy_k(write_iter, ".\n");
      WordWrapB(0);
      logerrputsb();
      reterr = kPglRetInconsistentInput;
    }
  }
  while (0) {
  CheckAlleleUniqueness_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CheckAlleleUniqueness_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
