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
#include "plink2_glm.h"
#include "plink2_matrix.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitGlm(GlmInfo* glm_info_ptr) {
  glm_info_ptr->flags = kfGlm0;
  glm_info_ptr->cols = kfGlmCol0;
  glm_info_ptr->mperm_ct = 0;
  glm_info_ptr->local_cat_ct = 0;
  glm_info_ptr->local_header_line_ct = 0;
  glm_info_ptr->local_chrom_col = 0;
  glm_info_ptr->local_bp_col = 0;
  glm_info_ptr->local_first_covar_col = 0;
  glm_info_ptr->max_corr = 0.999;
  glm_info_ptr->condition_varname = nullptr;
  glm_info_ptr->condition_list_fname = nullptr;
  InitRangeList(&(glm_info_ptr->parameters_range_list));
  InitRangeList(&(glm_info_ptr->tests_range_list));
}

void CleanupGlm(GlmInfo* glm_info_ptr) {
  free_cond(glm_info_ptr->condition_varname);
  free_cond(glm_info_ptr->condition_list_fname);
  CleanupRangeList(&(glm_info_ptr->parameters_range_list));
  CleanupRangeList(&(glm_info_ptr->tests_range_list));
}


// refugees from plink2_stats.h; important to remove its plink2_matrix.h
// dependency

// outer_buf = constraint_ct
// inner_buf = constraint_ct x constraint_ct
// tmphxs_buf and h_transpose_buf are constraint_ct x predictor_ct
// mi_buf only needs to be of length 2 * constraint_ct
BoolErr LinearHypothesisChisqF(const float* coef, const float* constraints_con_major, const float* cov_matrix, uint32_t constraint_ct, uint32_t predictor_ct, uint32_t cov_stride, double* chisq_ptr, float* tmphxs_buf, float* h_transpose_buf, float* inner_buf, double* half_inverted_buf, MatrixInvertBuf1* mi_buf, double* dbl_2d_buf, float* outer_buf) {
  ColMajorFvectorMatrixMultiplyStrided(coef, constraints_con_major, predictor_ct, predictor_ct, constraint_ct, outer_buf);
  // h-transpose does not have a special stride
  FmatrixTransposeCopy(constraints_con_major, constraint_ct, predictor_ct, predictor_ct, h_transpose_buf);
  ColMajorFmatrixMultiplyStrided(h_transpose_buf, cov_matrix, constraint_ct, constraint_ct, predictor_ct, cov_stride, predictor_ct, constraint_ct, tmphxs_buf);
  // tmp[][] is now predictor-major
  ColMajorFmatrixMultiplyStrided(tmphxs_buf, constraints_con_major, constraint_ct, constraint_ct, constraint_ct, predictor_ct, predictor_ct, constraint_ct, inner_buf);

  if (InvertFmatrixFirstHalf(constraint_ct, constraint_ct, inner_buf, half_inverted_buf, mi_buf, dbl_2d_buf)) {
    return 1;
  }
  InvertFmatrixSecondHalf(constraint_ct, constraint_ct, half_inverted_buf, inner_buf, mi_buf, dbl_2d_buf);
  double result = 0.0;
  const float* inner_iter = inner_buf;
  if (constraint_ct > kDotprodFThresh) {
    for (uint32_t constraint_idx = 0; constraint_idx != constraint_ct; ++constraint_idx) {
      result += S_CAST(double, DotprodF(inner_iter, outer_buf, constraint_ct) * outer_buf[constraint_idx]);
      inner_iter = &(inner_iter[constraint_ct]);
    }
  } else {
    for (uint32_t constraint_idx = 0; constraint_idx != constraint_ct; ++constraint_idx) {
      result += S_CAST(double, DotprodFShort(inner_iter, outer_buf, constraint_ct) * outer_buf[constraint_idx]);
      inner_iter = &(inner_iter[constraint_ct]);
    }
  }
  if (result < 0.0) {
    // guard against floating point error
    result = 0.0;
  }
  *chisq_ptr = result;
  return 0;
}

BoolErr LinearHypothesisChisq(const double* coef, const double* constraints_con_major, const double* cov_matrix, uintptr_t constraint_ct, uintptr_t predictor_ct, double* chisq_ptr, double* tmphxs_buf, double* h_transpose_buf, double* inner_buf, MatrixInvertBuf1* mi_buf, double* outer_buf) {
  // See PLINK model.cpp Model::linearHypothesis().
  //
  // outer_buf = constraint_ct
  // inner_buf = constraint_ct x constraint_ct
  // tmphxs_buf and h_transpose_buf are constraint_ct x predictor_ct
  // mi_buf only needs to be of length 2 * constraint_ct
  //
  // Since no PLINK function ever calls this with nonzero h[] values, this just
  // takes a df (constraint_ct) parameter for now; it's trivial to switch to
  // the more general interface later.
  ColMajorVectorMatrixMultiplyStrided(coef, constraints_con_major, predictor_ct, predictor_ct, constraint_ct, outer_buf);
  MatrixTransposeCopy(constraints_con_major, constraint_ct, predictor_ct, h_transpose_buf);
  ColMajorMatrixMultiply(h_transpose_buf, cov_matrix, constraint_ct, predictor_ct, predictor_ct, tmphxs_buf);
  // tmp[][] is now predictor-major
  ColMajorMatrixMultiply(tmphxs_buf, constraints_con_major, constraint_ct, constraint_ct, predictor_ct, inner_buf);

  // don't need H-transpose any more, so we can use h_transpose_buf for matrix
  // inversion
  if (InvertMatrix(constraint_ct, inner_buf, mi_buf, h_transpose_buf)) {
    return 1;
  }
  double result = 0.0;
  const double* inner_iter = inner_buf;
  if (constraint_ct > kDotprodDThresh) {
    for (uintptr_t constraint_idx = 0; constraint_idx != constraint_ct; ++constraint_idx) {
      result += DotprodD(inner_iter, outer_buf, constraint_ct) * outer_buf[constraint_idx];
      inner_iter = &(inner_iter[constraint_ct]);
    }
  } else {
    for (uintptr_t constraint_idx = 0; constraint_idx != constraint_ct; ++constraint_idx) {
      result += DotprodDShort(inner_iter, outer_buf, constraint_ct) * outer_buf[constraint_idx];
      inner_iter = &(inner_iter[constraint_ct]);
    }
  }
  if (result < 0.0) {
    // guard against floating point error
    result = 0.0;
  }
  *chisq_ptr = result;
  return 0;
}


PglErr GlmLocalOpen(const char* local_covar_fname, const char* local_pvar_fname, const char* local_psam_fname, const char* sample_ids, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const GlmInfo* glm_info_ptr, uint32_t raw_sample_ct, uintptr_t max_sample_id_blen, uint32_t raw_variant_ct, const uintptr_t** sample_include_ptr, const uintptr_t** sex_nm_ptr, const uintptr_t** sex_male_ptr, const uintptr_t** variant_include_ptr, uint32_t* sample_ct_ptr, uint32_t* variant_ct_ptr, TextStream* local_covar_txsp, uint32_t** local_sample_uidx_order_ptr, uintptr_t** local_variant_include_ptr, uint32_t* local_sample_ct_ptr, uint32_t* local_variant_ctl_ptr, uint32_t* local_covar_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  const char* txs_fname = nullptr;
  textFILE local_covar_txf;
  TextStream txs;
  PreinitTextFile(&local_covar_txf);
  PreinitTextStream(&txs);
  {
    // 1. read .psam/.fam file, update sample_ct, initialize
    //    local_sample_uidx_order (use OpenAndLoadXidHeader()?)
    uint32_t max_line_blen;
    if (unlikely(StandardizeMaxLineBlen(bigstack_left() / 4, &max_line_blen))) {
      goto GlmLocalOpen_ret_NOMEM;
    }
    txs_fname = local_psam_fname;
    reterr = InitTextStreamEx(local_psam_fname, 1, kMaxLongLine, max_line_blen, 1, &txs);
    if (unlikely(reterr)) {
      goto GlmLocalOpen_ret_TSTREAM_FAIL;
    }
    uint32_t is_header_line;
    char* line_start;
    do {
      ++line_idx;
      line_start = TextGet(&txs);
      if (unlikely(!line_start)) {
        if (!TextStreamErrcode2(&txs, &reterr)) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", local_psam_fname);
          goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
        }
        goto GlmLocalOpen_ret_TSTREAM_FAIL;
      }
      is_header_line = (line_start[0] == '#');
    } while (is_header_line && (!tokequal_k(&(line_start[1]), "FID")) && (!tokequal_k(&(line_start[1]), "IID")));
    XidMode xid_mode = kfXidModeFidIid;
    if (is_header_line) {
      if (line_start[1] == 'I') {
        xid_mode = kfXidModeIid;
      }
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    char* sorted_sample_idbox;
    uint32_t* sample_id_map;
    uintptr_t* new_sample_include;
    char* idbuf;
    if (unlikely(
            bigstack_end_alloc_c(orig_sample_ct * max_sample_id_blen, &sorted_sample_idbox) ||
            bigstack_end_alloc_u32(orig_sample_ct, &sample_id_map) ||
            bigstack_end_calloc_w(raw_sample_ctl, &new_sample_include) ||
            bigstack_end_alloc_c(max_sample_id_blen, &idbuf))) {
      goto GlmLocalOpen_ret_NOMEM;
    }
    // (don't permit duplicate FID+IID for now, but maybe we'll want to use
    // xid interface later?)
    reterr = CopySortStrboxSubsetNoalloc(*sample_include_ptr, sample_ids, orig_sample_ct, max_sample_id_blen, 0, 0, 0, sorted_sample_idbox, sample_id_map);
    if (unlikely(reterr)) {
      goto GlmLocalOpen_ret_1;
    }
    uint32_t* local_sample_uidx_order = R_CAST(uint32_t*, g_bigstack_base);
    uintptr_t max_local_sample_ct = RoundDownPow2(bigstack_left(), kCacheline) / sizeof(int32_t);
#ifdef __LP64__
    if (max_local_sample_ct > kMaxLongLine / 2) {
      max_local_sample_ct = kMaxLongLine / 2;
    }
#endif
    uintptr_t local_sample_ct = 0;
    if (is_header_line) {
      ++line_idx;
      line_start = TextGet(&txs);
    }
    for (; line_start; ++local_sample_ct, ++line_idx, line_start = TextGet(&txs)) {
      if (unlikely(line_start[0] == '#')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, local_psam_fname);
        goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
      }
      if (unlikely(local_sample_ct == max_local_sample_ct)) {
#ifdef __LP64__
        if (local_sample_ct == kMaxLongLine / 2) {
          snprintf(g_logbuf, kLogbufSize, "Error: Too many samples in %s.\n", local_psam_fname);
          goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
        }
#endif
        goto GlmLocalOpen_ret_NOMEM;
      }
      const char* read_ptr = line_start;
      uint32_t sample_uidx;
      if (!SortedXidboxReadFind(sorted_sample_idbox, sample_id_map, max_sample_id_blen, orig_sample_ct, 0, xid_mode, &read_ptr, &sample_uidx, idbuf)) {
        if (unlikely(IsSet(new_sample_include, sample_uidx))) {
          char* first_tab = AdvToDelim(idbuf, '\t');
          *first_tab = ' ';
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate ID '%s' in %s.\n", idbuf, local_psam_fname);
          goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
        }
        SetBit(sample_uidx, new_sample_include);
        local_sample_uidx_order[local_sample_ct] = sample_uidx;
      } else {
        if (unlikely(!read_ptr)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Fewer tokens than expected on line %" PRIuPTR " of %s.\n", line_idx, local_psam_fname);
          goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
        }
        local_sample_uidx_order[local_sample_ct] = UINT32_MAX;
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto GlmLocalOpen_ret_TSTREAM_FAIL;
    }
    BigstackFinalizeU32(local_sample_uidx_order, local_sample_ct);
    *local_sample_uidx_order_ptr = local_sample_uidx_order;
    *local_sample_ct_ptr = local_sample_ct;
    const uint32_t new_sample_ct = PopcountWords(new_sample_include, raw_sample_ctl);
    assert(new_sample_ct <= orig_sample_ct);
    if (new_sample_ct < orig_sample_ct) {
      uintptr_t* sample_include_copy;
      uintptr_t* sex_nm_copy;
      uintptr_t* sex_male_copy;
      if (unlikely(
              bigstack_alloc_w(raw_sample_ctl, &sample_include_copy) ||
              bigstack_alloc_w(raw_sample_ctl, &sex_nm_copy) ||
              bigstack_alloc_w(raw_sample_ctl, &sex_male_copy))) {
        goto GlmLocalOpen_ret_NOMEM;
      }
      memcpy(sample_include_copy, new_sample_include, raw_sample_ctl * sizeof(intptr_t));
      BitvecAndCopy(sample_include_copy, *sex_nm_ptr, raw_sample_ctl, sex_nm_copy);
      *sex_nm_ptr = sex_nm_copy;
      BitvecAndCopy(sample_include_copy, *sex_male_ptr, raw_sample_ctl, sex_male_copy);
      *sex_male_ptr = sex_male_copy;
      *sample_include_ptr = sample_include_copy;
      *sample_ct_ptr = new_sample_ct;
    }
    BigstackEndReset(TextStreamMemStart(&txs));

    // 2. if local-pvar= specified, read .pvar/.bim file, update variant_ct,
    //    initialize local_variant_ct/local_variant_include.
    if (local_pvar_fname) {
      reterr = TextRetarget(local_pvar_fname, &txs);
      if (unlikely(reterr)) {
        goto GlmLocalOpen_ret_TSTREAM_FAIL;
      }
      txs_fname = local_pvar_fname;
      line_idx = 0;
      do {
        ++line_idx;
        line_start = TextGet(&txs);
        if (unlikely(!line_start)) {
          if (!TextStreamErrcode2(&txs, &reterr)) {
            snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", local_pvar_fname);
            goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
          }
          goto GlmLocalOpen_ret_TSTREAM_FAIL;
        }
        is_header_line = (line_start[0] == '#');
      } while (is_header_line && (!tokequal_k(line_start, "#CHROM")));
      uint32_t col_skips[2];
      uint32_t col_types[2];
      // uint32_t relevant_postchr_col_ct = 2;
      if (is_header_line) {
        // parse header
        // [-1] = #CHROM (must be first column)
        // [0] = POS
        // [1] = ID
        // don't care about the rest
        const char* token_end = &(line_start[6]);
        uint32_t found_header_bitset = 0;
        uint32_t relevant_postchr_col_ct = 0;
        const char* linebuf_iter;
        for (uint32_t col_idx = 1; ; ++col_idx) {
          linebuf_iter = FirstNonTspace(token_end);
          if (IsEolnKns(*linebuf_iter)) {
            break;
          }
          token_end = CurTokenEnd(linebuf_iter);
          const uint32_t token_slen = token_end - linebuf_iter;
          uint32_t cur_col_type;
          if (strequal_k(linebuf_iter, "POS", token_slen)) {
            cur_col_type = 0;
          } else if (strequal_k(linebuf_iter, "ID", token_slen)) {
            cur_col_type = 1;
          } else {
            continue;
          }
          const uint32_t cur_col_type_shifted = 1 << cur_col_type;
          if (unlikely(found_header_bitset & cur_col_type_shifted)) {
            // Known column type, so length is limited and we don't have to
            // worry about buffer overflow.
            char* write_iter = strcpya_k(g_logbuf, "Error: Duplicate column header '");
            write_iter = memcpya(write_iter, linebuf_iter, token_slen);
            write_iter = strcpya_k(write_iter, "' on line ");
            write_iter = wtoa(line_idx, write_iter);
            write_iter = strcpya_k(write_iter, " of ");
            write_iter = strcpya(write_iter, local_pvar_fname);
            memcpy_k(write_iter, ".\n\0", 4);
            goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
          }
          found_header_bitset |= cur_col_type_shifted;
          col_skips[relevant_postchr_col_ct] = col_idx;
          col_types[relevant_postchr_col_ct++] = cur_col_type;
        }
        if (unlikely(found_header_bitset != 3)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Missing column header(s) on line %" PRIuPTR " of %s. (POS and ID are required.)\n", line_idx, local_pvar_fname);
          goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
        }
        for (uint32_t rpc_col_idx = relevant_postchr_col_ct - 1; rpc_col_idx; --rpc_col_idx) {
          col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
        }
      } else {
        col_types[0] = 1;
        col_types[1] = 0;
        col_skips[0] = 1;
        // CM column may be omitted
        const char* linebuf_iter = NextTokenMult(line_start, 4);
        if (unlikely(!linebuf_iter)) {
          goto GlmLocalOpen_ret_MISSING_TOKENS_PVAR;
        }
        linebuf_iter = NextToken(linebuf_iter);
        if (!linebuf_iter) {
          // #CHROM ID POS ALT REF
          col_skips[1] = 1;
        } else {
          // #CHROM ID CM POS ALT REF
          col_skips[1] = 2;
        }
      }
      const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
      uintptr_t* new_variant_include;
      if (unlikely(bigstack_end_calloc_w(raw_variant_ctl, &new_variant_include))) {
        goto GlmLocalOpen_ret_NOMEM;
      }
      uint32_t max_local_variant_ct = 0x7ffffffd;
      if (bigstack_left() < (0x80000000U / CHAR_BIT)) {
        max_local_variant_ct = RoundDownPow2(bigstack_left(), kCacheline) * CHAR_BIT;
      }
      const uintptr_t* orig_variant_include = *variant_include_ptr;
      uintptr_t* local_variant_include = R_CAST(uintptr_t*, g_bigstack_base);
      *local_variant_include_ptr = local_variant_include;
      uint32_t local_variant_ct = 0;
      uint32_t new_variant_ct = 0;
      uint32_t prev_variant_uidx = AdvTo1Bit(orig_variant_include, 0);
      uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, prev_variant_uidx);
      uint32_t prev_chr_code = cip->chr_file_order[chr_fo_idx];
      uint32_t prev_bp = variant_bps[prev_variant_uidx];
      uint32_t chr_end = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[prev_chr_code] + 1];
      if (is_header_line) {
        ++line_idx;
        line_start = TextGet(&txs);
      }
      for (; line_start; ++local_variant_ct, ++line_idx, line_start = TextGet(&txs)) {
        if (unlikely(line_start[0] == '#')) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #CHROM header line is present it must denote the end of the header block.)\n", line_idx, local_pvar_fname);
          goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
        }
        if (unlikely(local_variant_ct == max_local_variant_ct)) {
          if (max_local_variant_ct == 0x7ffffffd) {
            snprintf(g_logbuf, kLogbufSize, "Error: Too many samples in %s.\n", local_pvar_fname);
            goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
          }
          goto GlmLocalOpen_ret_NOMEM;
        }
        if (!(local_variant_ct % kBitsPerWord)) {
          local_variant_include[local_variant_ct / kBitsPerWord] = 0;
        }
        char* line_iter = CurTokenEnd(line_start);
        // #CHROM
        if (unlikely((*line_iter != '\t') && (*line_iter != ' '))) {
          goto GlmLocalOpen_ret_MISSING_TOKENS_PVAR;
        }
        const uint32_t cur_chr_code = GetChrCodeCounted(cip, line_iter - line_start, line_start);
        if (IsI32Neg(cur_chr_code)) {
          continue;
        }
        const uint32_t cur_chr_code_u = cur_chr_code;
        if (cur_chr_code_u != prev_chr_code) {
          uint32_t first_variant_uidx_in_chr = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[cur_chr_code_u]];
          if (first_variant_uidx_in_chr < prev_variant_uidx) {
            if (unlikely(new_variant_ct)) {
              // not worth the trouble of handling this
              snprintf(g_logbuf, kLogbufSize, "Error: Chromosome order in %s is different from chromosome order in main dataset.\n", local_pvar_fname);
              goto GlmLocalOpen_ret_INCONSISTENT_INPUT_WW;
            }
            continue;
          }
          prev_variant_uidx = AdvBoundedTo1Bit(orig_variant_include, first_variant_uidx_in_chr, raw_variant_ct);
          if (prev_variant_uidx == raw_variant_ct) {
            break;
          }
          chr_fo_idx = GetVariantChrFoIdx(cip, prev_variant_uidx);
          prev_chr_code = cip->chr_file_order[chr_fo_idx];
          prev_bp = variant_bps[prev_variant_uidx];
          chr_end = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[prev_chr_code] + 1];
          if (cur_chr_code_u != prev_chr_code) {
            continue;
          }
        }
        char* token_ptrs[2];
        uint32_t token_slens[2];
        if (unlikely(!TokenLex(line_iter, col_types, col_skips, 2, token_ptrs, token_slens))) {
          goto GlmLocalOpen_ret_MISSING_TOKENS_PVAR;
        }
        // POS
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(token_ptrs[0], &cur_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, local_pvar_fname);
          goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
        }
        if (cur_bp < S_CAST(int32_t, prev_bp)) {
          continue;
        }
        const uint32_t cur_bp_u = cur_bp;
        if (cur_bp_u > prev_bp) {
          do {
            prev_variant_uidx = AdvBoundedTo1Bit(orig_variant_include, prev_variant_uidx + 1, raw_variant_ct);
          } while ((prev_variant_uidx < chr_end) && (cur_bp_u > variant_bps[prev_variant_uidx]));
          if (prev_variant_uidx >= chr_end) {
            goto GlmLocalOpen_skip_variant_and_update_chr;
          }
          prev_bp = variant_bps[prev_variant_uidx];
        }
        if (cur_bp_u == prev_bp) {
          // ID
          // note that if there are two same-position variants which appear in
          // a different order in the main dataset and the local-pvar file, one
          // will be skipped.  (probably want to add a warning in this case.)
          const uint32_t variant_id_slen = token_slens[1];
          char* cur_variant_id = token_ptrs[1];
          cur_variant_id[variant_id_slen] = '\0';
          const uint32_t variant_id_blen = variant_id_slen + 1;
          do {
            const char* loaded_variant_id = variant_ids[prev_variant_uidx];
            if (memequal(cur_variant_id, loaded_variant_id, variant_id_blen)) {
              if (unlikely(IsSet(new_variant_include, prev_variant_uidx))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate ID (with duplicate CHROM/POS) '%s' in %s.\n", cur_variant_id, local_pvar_fname);
                goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
              }
              SetBit(prev_variant_uidx, new_variant_include);
              ++new_variant_ct;
              SetBit(local_variant_ct, local_variant_include);
              break;
            }
            prev_variant_uidx = AdvBoundedTo1Bit(orig_variant_include, prev_variant_uidx + 1, raw_variant_ct);
            if (prev_variant_uidx >= chr_end) {
              goto GlmLocalOpen_skip_variant_and_update_chr;
            }
            prev_bp = variant_bps[prev_variant_uidx];
          } while (cur_bp_u == prev_bp);
        }
        continue;
      GlmLocalOpen_skip_variant_and_update_chr:
        if (prev_variant_uidx == raw_variant_ct) {
          continue;
        }
        chr_fo_idx = GetVariantChrFoIdx(cip, prev_variant_uidx);
        prev_chr_code = cip->chr_file_order[chr_fo_idx];
        prev_bp = variant_bps[prev_variant_uidx];
        chr_end = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[prev_chr_code] + 1];
      }
      if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
        goto GlmLocalOpen_ret_TSTREAM_FAIL;
      }
      BigstackFinalizeW(local_variant_include, BitCtToWordCt(local_variant_ct));
      *local_variant_ctl_ptr = BitCtToWordCt(local_variant_ct);
      assert(new_variant_ct <= *variant_ct_ptr);
      if (new_variant_ct < *variant_ct_ptr) {
        uintptr_t* variant_include_copy;
        if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include_copy))) {
          goto GlmLocalOpen_ret_NOMEM;
        }
        memcpy(variant_include_copy, new_variant_include, raw_variant_ctl * sizeof(intptr_t));
        *variant_include_ptr = variant_include_copy;
        *variant_ct_ptr = new_variant_ct;
      }
    }

    // 3. If not local-cats[0]=, scan first line of local-covar= file to
    //    determine covariate count.
    if (CleanupTextStream2(local_pvar_fname, &txs, &reterr)) {
      goto GlmLocalOpen_ret_1;
    }
    txs_fname = nullptr;
    BigstackEndReset(bigstack_end_mark);
    uintptr_t loadbuf_size = bigstack_left();
    if (loadbuf_size > kMaxLongLine) {
      loadbuf_size = kMaxLongLine;
    } else if (unlikely(loadbuf_size < kDecompressMinCapacity)) {
      // this shouldn't be possible
      goto GlmLocalOpen_ret_NOMEM;
    }
    // don't formally allocate
    char* dst = R_CAST(char*, g_bigstack_base);
    reterr = TextFileOpenEx(local_covar_fname, kMaxLongLine, loadbuf_size, dst, &local_covar_txf);
    if (unlikely(reterr)) {
      goto GlmLocalOpen_ret_RFILE_FAIL;
    }
    const uint32_t local_sample_or_hap_ct = local_sample_ct << ((glm_info_ptr->flags / kfGlmLocalHaps) & 1);
    const uint32_t skip_ct = glm_info_ptr->local_first_covar_col? (glm_info_ptr->local_first_covar_col - 1) : 0;
    uint32_t local_covar_ct;
    if (!glm_info_ptr->local_cat_ct) {
      const uint32_t first_nonheader_line_idx = glm_info_ptr->local_header_line_ct + 1;
      for (line_idx = 1; line_idx <= first_nonheader_line_idx; ++line_idx) {
        line_start = TextFileGet(&local_covar_txf);
        if (unlikely(!line_start)) {
          if (!TextFileErrcode2(&local_covar_txf, &reterr)) {
            snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", local_covar_fname);
            goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
          }
          goto GlmLocalOpen_ret_RFILE_FAIL;
        }
      }
      const uint32_t token_ct = CountTokens(line_start);
      local_covar_ct = (token_ct - skip_ct) / local_sample_or_hap_ct;
      if (unlikely((token_ct < skip_ct) || (local_covar_ct * local_sample_or_hap_ct + skip_ct != token_ct))) {
        if (!skip_ct) {
          snprintf(g_logbuf, kLogbufSize, "Error: Unexpected token count on line %u of %s (%u, %smultiple of %u expected).\n", first_nonheader_line_idx, local_covar_fname, token_ct, (token_ct == local_sample_or_hap_ct)? "larger " : "", local_sample_or_hap_ct);
        } else {
          snprintf(g_logbuf, kLogbufSize, "Error: Unexpected token count on line %u of %s (%u, %u more than %smultiple of %u expected).\n", first_nonheader_line_idx, local_covar_fname, token_ct, skip_ct, (token_ct == local_sample_or_hap_ct)? "larger " : "", local_sample_or_hap_ct);
        }
        goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
      }
      if (glm_info_ptr->flags & kfGlmLocalOmitLast) {
        if (unlikely(local_covar_ct == 1)) {
          logerrputs("Error: --glm 'local-omit-last' modifier cannot be used when there is only one\nlocal covariate.\n");
          goto GlmLocalOpen_ret_INCONSISTENT_INPUT;
        }
        logprintf("--glm local-covar=: %u local covariates present, %u used.\n", local_covar_ct, local_covar_ct - 1);
        --local_covar_ct;
      } else {
        logprintf("--glm local-covar=: %u local covariate%s present.\n", local_covar_ct, (local_covar_ct == 1)? "" : "s");
      }
    } else {
      local_covar_ct = glm_info_ptr->local_cat_ct - 1;
      if (unlikely(local_covar_ct * S_CAST(uint64_t, local_sample_or_hap_ct) > (kMaxLongLine / 2) - skip_ct)) {
        snprintf(g_logbuf, kLogbufSize, "Error: <# samples/haplotypes in %s> * <# categories - 1> too large (limited to %u).\n", local_covar_fname, kMaxLongLine / 2 - skip_ct);
        goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
      }
    }
    *local_covar_ct_ptr = local_covar_ct;

    // 4. Initialize final local-covar reader.
    // Quasi-bugfix (6 Jun 2020): even in local-cats case, a header line may
    // have a floating point number per (sample, covariate) tuple.
    // Permit 24 characters per floating point number instead of 16, since some
    // tools dump 15-17 significant digits.
    uint64_t enforced_max_line_blen = local_sample_or_hap_ct * 24LLU;
    if (!glm_info_ptr->local_cat_ct) {
      enforced_max_line_blen *= local_covar_ct + ((glm_info_ptr->flags / kfGlmLocalOmitLast) & 1);
    } else {
      // Separate ID per (sample, category) possible on header line here.
      enforced_max_line_blen *= glm_info_ptr->local_cat_ct;
    }
    // \r\n
    enforced_max_line_blen += 2;
    if (unlikely(enforced_max_line_blen > kMaxLongLine)) {
      logerrputs("Error: Too many samples/covariates for --glm local-covar=.\n");
      goto GlmLocalOpen_ret_MALFORMED_INPUT;
    }
    if (enforced_max_line_blen < kDecompressMinBlen) {
      enforced_max_line_blen = kDecompressMinBlen;
    }
    uint32_t dst_capacity = enforced_max_line_blen + kDecompressChunkSize;
    if (TextFileLinebufLen(&local_covar_txf) > dst_capacity) {
      dst_capacity = TextFileLinebufLen(&local_covar_txf);
    }
    dst_capacity = RoundUpPow2(dst_capacity, kCacheline);
    if (bigstack_left() < dst_capacity) {
      goto GlmLocalOpen_ret_NOMEM;
    }
    g_bigstack_base = R_CAST(unsigned char*, &(dst[dst_capacity]));
    reterr = TextStreamOpenEx(nullptr, enforced_max_line_blen, dst_capacity, 1, &local_covar_txf, nullptr, local_covar_txsp);
    if (unlikely(reterr)) {
      TextStreamErrPrint(local_covar_fname, local_covar_txsp);
      goto GlmLocalOpen_ret_1;
    }
    bigstack_mark = g_bigstack_base;
  }
  while (0) {
  GlmLocalOpen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GlmLocalOpen_ret_TSTREAM_FAIL:
    TextStreamErrPrint(txs_fname, &txs);
    break;
  GlmLocalOpen_ret_RFILE_FAIL:
    TextFileErrPrint(local_covar_fname, &local_covar_txf);
    break;
  GlmLocalOpen_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  GlmLocalOpen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  GlmLocalOpen_ret_MISSING_TOKENS_PVAR:
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, local_pvar_fname);
    reterr = kPglRetMalformedInput;
    break;
  GlmLocalOpen_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  GlmLocalOpen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 GlmLocalOpen_ret_1:
  if (txs_fname) {
    CleanupTextStream2(txs_fname, &txs, &reterr);
  }
  CleanupTextFile2(local_covar_fname, &local_covar_txf, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

uint32_t CheckForSeparatedQtCovar(const uintptr_t* pheno_cc, const uintptr_t* cur_sample_include, const double* covar_vals, uint32_t sample_ct, double* covar_val_keep_ptr) {
  const uint32_t first_sample_uidx = AdvTo1Bit(cur_sample_include, 0);
  double cur_covar_val = covar_vals[first_sample_uidx];
  double ctrl_min;
  double ctrl_max;
  double case_min;
  double case_max;
  if (IsSet(pheno_cc, first_sample_uidx)) {
    ctrl_min = DBL_MAX;
    ctrl_max = -DBL_MAX;
    case_min = cur_covar_val;
    case_max = cur_covar_val;
  } else {
    ctrl_min = cur_covar_val;
    ctrl_max = cur_covar_val;
    case_min = DBL_MAX;
    case_max = -DBL_MAX;
  }
  uintptr_t sample_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(cur_sample_include, first_sample_uidx + 1, &sample_uidx_base, &cur_bits);
  for (uint32_t sample_idx = 1; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(cur_sample_include, &sample_uidx_base, &cur_bits);
    cur_covar_val = covar_vals[sample_uidx];
    if (IsSet(pheno_cc, sample_uidx)) {
      if (cur_covar_val < case_min) {
        if ((ctrl_min < case_max) && (cur_covar_val < ctrl_max)) {
          return 0;
        }
        case_min = cur_covar_val;
      } else if (cur_covar_val > case_max) {
        if ((ctrl_max > case_min) && (cur_covar_val > ctrl_min)) {
          return 0;
        }
        case_max = cur_covar_val;
      }
    } else {
      if (cur_covar_val < ctrl_min) {
        if ((case_min < ctrl_max) && (cur_covar_val < case_max)) {
          return 0;
        }
        ctrl_min = cur_covar_val;
      } else if (cur_covar_val > ctrl_max) {
        if ((case_max > ctrl_min) && (cur_covar_val > case_min)) {
          return 0;
        }
        ctrl_max = cur_covar_val;
      }
    }
  }
  if ((case_min > ctrl_max) || (ctrl_min > case_max)) {
    return 2;
  }
  if (covar_val_keep_ptr) {
    *covar_val_keep_ptr = (case_min == ctrl_max)? case_min : case_max;
  }
  return 1;
}

uint32_t CheckForSeparatedCatCovar(const uintptr_t* pheno_cc, const uintptr_t* cur_sample_include, const PhenoCol* covar_col, uint32_t sample_ct, uint32_t* case_and_ctrl_cat_ct_ptr, uintptr_t* cat_covar_wkspace) {
  const uint32_t first_sample_uidx = AdvTo1Bit(cur_sample_include, 0);
  const uint32_t nonnull_cat_ct = covar_col->nonnull_category_ct;
  const uint32_t cur_word_ct = 1 + nonnull_cat_ct / kBitsPerWordD2;
  ZeroWArr(cur_word_ct, cat_covar_wkspace);
  const uint32_t* covar_cats = covar_col->data.cat;
  // If no remaining categories have both cases and controls, we have complete
  // separation.
  // If some do and some do not, we have quasi-complete separation, and when
  // performing ordinary logistic regression, we must remove samples in the
  // all-case and/or all-control categories.
  uintptr_t sample_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(cur_sample_include, first_sample_uidx, &sample_uidx_base, &cur_bits);
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(cur_sample_include, &sample_uidx_base, &cur_bits);
    const uint32_t cur_cat_idx = covar_cats[sample_uidx];
    // Odd bits represent presence of a case, even bits represent presence of a
    // control.
    const uint32_t is_case = IsSet(pheno_cc, sample_uidx);
    SetBit(cur_cat_idx * 2 + is_case, cat_covar_wkspace);
  }
  uint32_t case_and_ctrl_cat_ct = 0;
  uint32_t pheno_by_cat_ct = 0;
  for (uint32_t widx = 0; widx != cur_word_ct; ++widx) {
    const uintptr_t cur_word = cat_covar_wkspace[widx];
    case_and_ctrl_cat_ct += Popcount01Word(Word11(cur_word));
    pheno_by_cat_ct += PopcountWord(cur_word);
  }
  if (!case_and_ctrl_cat_ct) {
    return 2;
  }
  if (case_and_ctrl_cat_ct == pheno_by_cat_ct * 2) {
    // all categories contain both cases and controls.
    return 0;
  }
  if (case_and_ctrl_cat_ct_ptr) {
    *case_and_ctrl_cat_ct_ptr = case_and_ctrl_cat_ct;
  }
  return 1;
}

// Returns 0 if not separated, 1 if quasi-separated, 2 if fully separated.
uint32_t CheckForSeparatedCovar(const uintptr_t* pheno_cc, const uintptr_t* cur_sample_include, const PhenoCol* covar_col, uint32_t sample_ct, uintptr_t* cat_covar_wkspace) {
  if (sample_ct < 2) {
    return 2;
  }
  if (covar_col->type_code == kPhenoDtypeOther) {
    // local covariate, no precomputation possible
    return 0;
  }
  if (covar_col->type_code == kPhenoDtypeQt) {
    return CheckForSeparatedQtCovar(pheno_cc, cur_sample_include, covar_col->data.qt, sample_ct, nullptr);
  }
  return CheckForSeparatedCatCovar(pheno_cc, cur_sample_include, covar_col, sample_ct, nullptr, cat_covar_wkspace);
}

// Returns 1 if phenotype fully separates the covariate, and the non-Firth
// logistic regression should be aborted.
// Otherwise, if we have quasi-separation, use Stata approach of throwing out
// the covariate, and keeping only samples which have the one covariate value
// present in both cases and controls.
// Note that covar_ct is not a parameter; caller is responsible for
// re-popcounting covar_include.
BoolErr CheckForAndHandleSeparatedCovar(const uintptr_t* pheno_cc, const PhenoCol* covar_cols, uint32_t raw_sample_ctl, uint32_t covar_uidx, uintptr_t* cur_sample_include, uintptr_t* covar_include, uint32_t* sample_ct_ptr, uintptr_t* cat_covar_wkspace) {
  uint32_t sample_ct = *sample_ct_ptr;
  if (sample_ct < 2) {
    return 1;
  }
  const PhenoCol* cur_covar_col = &(covar_cols[covar_uidx]);
  if (cur_covar_col->type_code == kPhenoDtypeOther) {
    return 0;
  }
  const uint32_t first_sample_uidx = AdvTo1Bit(cur_sample_include, 0);
  if (cur_covar_col->type_code == kPhenoDtypeQt) {
    const double* covar_vals = cur_covar_col->data.qt;
    double covar_val_keep;
    const uint32_t separation_type = CheckForSeparatedQtCovar(pheno_cc, cur_sample_include, covar_vals, sample_ct, &covar_val_keep);
    if (separation_type != 1) {
      return separation_type / 2;
    }
    uintptr_t sample_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(cur_sample_include, first_sample_uidx, &sample_uidx_base, &cur_bits);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(cur_sample_include, &sample_uidx_base, &cur_bits);
      if (covar_vals[sample_uidx] != covar_val_keep) {
        ClearBit(sample_uidx, cur_sample_include);
      }
    }
    *sample_ct_ptr = PopcountWords(cur_sample_include, raw_sample_ctl);
    ClearBit(covar_uidx, covar_include);
    return 0;
  }
  uint32_t case_and_ctrl_cat_ct;
  const uint32_t separation_type = CheckForSeparatedCatCovar(pheno_cc, cur_sample_include, cur_covar_col, sample_ct, &case_and_ctrl_cat_ct, cat_covar_wkspace);
  if (separation_type != 1) {
    return separation_type / 2;
  }
  // more than one category contains both cases and controls (so we don't need
  // to remove the categorical covariate), but at least one does not, so we
  // still have to prune some samples.
  const uint32_t nonnull_cat_ct = cur_covar_col->nonnull_category_ct;
  const uint32_t cur_word_ct = 1 + nonnull_cat_ct / kBitsPerWordD2;
  for (uint32_t widx = 0; widx != cur_word_ct; ++widx) {
    const uintptr_t cur_word = cat_covar_wkspace[widx];
    cat_covar_wkspace[widx] = Word11(cur_word);
  }
  uintptr_t sample_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(cur_sample_include, first_sample_uidx, &sample_uidx_base, &cur_bits);
  const uint32_t* covar_cats = cur_covar_col->data.cat;
  // bugfix (7 Jul 2018): this needs to check sample_idx, not sample_uidx
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(cur_sample_include, &sample_uidx_base, &cur_bits);
    if (!IsSet(cat_covar_wkspace, covar_cats[sample_uidx] * 2)) {
      ClearBit(sample_uidx, cur_sample_include);
    }
  }
  *sample_ct_ptr = PopcountWords(cur_sample_include, raw_sample_ctl);
  if (case_and_ctrl_cat_ct == 1) {
    ClearBit(covar_uidx, covar_include);
  }
  return 0;
}

BoolErr GlmDetermineCovars(const uintptr_t* pheno_cc, const uintptr_t* initial_covar_include, const PhenoCol* covar_cols, uint32_t raw_sample_ct, uint32_t raw_covar_ctl, uint32_t initial_covar_ct, uint32_t covar_max_nonnull_cat_ct, uint32_t is_sometimes_firth, uint32_t is_always_firth, uintptr_t* cur_sample_include, uintptr_t* covar_include, uint32_t* sample_ct_ptr, uint32_t* covar_ct_ptr, uint32_t* extra_cat_ct_ptr, uint16_t* separation_found_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  BoolErr reterr = 0;
  {
    // bugfix (11 Mar 2020): if no covariates provided, covar_include ==
    // nullptr is possible, and covar_include[0] dereference below is optimized
    // out by some but not all compilers.
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    if (!covar_include) {
      *sample_ct_ptr = PopcountWords(cur_sample_include, raw_sample_ctl);
      goto GlmDetermineCovars_ret_NOCOVAR;
    }

    uintptr_t* sample_include_backup;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &sample_include_backup))) {
      goto GlmDetermineCovars_ret_NOMEM;
    }
    memcpy(sample_include_backup, cur_sample_include, raw_sample_ctl * sizeof(intptr_t));
    memcpy(covar_include, initial_covar_include, raw_covar_ctl * sizeof(intptr_t));
    // bugfix (28 Apr 2018): initial_covar_include can contain an extra sex set
    // bit at the end.
    if (PopcountWords(initial_covar_include, raw_covar_ctl) > initial_covar_ct) {
      uintptr_t last_word = covar_include[raw_covar_ctl - 1];
      const uint32_t highest_set_bit_idx = bsrw(last_word);
      covar_include[raw_covar_ctl - 1] -= k1LU << highest_set_bit_idx;
    }

    // 1. Determine samples for which all phenotype/covariate values are
    //    present, then provisionally remove the covariates which are constant
    //    over that set in linear case, or produce separation in logistic case
    uintptr_t covar_uidx_base = 0;
    uintptr_t cur_bits = covar_include[0];
    for (uint32_t covar_idx = 0; covar_idx != initial_covar_ct; ++covar_idx) {
      const uintptr_t covar_uidx = BitIter1(covar_include, &covar_uidx_base, &cur_bits);
      if (covar_cols[covar_uidx].nonmiss) {
        BitvecAnd(covar_cols[covar_uidx].nonmiss, raw_sample_ctl, cur_sample_include);
      }
    }
    uint32_t prev_sample_ct = PopcountWords(cur_sample_include, raw_sample_ctl);
    covar_uidx_base = 0;
    cur_bits = initial_covar_include[0];
    for (uint32_t covar_idx = 0; covar_idx != initial_covar_ct; ++covar_idx) {
      const uintptr_t covar_uidx = BitIter1(initial_covar_include, &covar_uidx_base, &cur_bits);
      if ((covar_cols[covar_uidx].type_code != kPhenoDtypeOther) && IsConstCovar(&(covar_cols[covar_uidx]), cur_sample_include, prev_sample_ct)) {
        ClearBit(covar_uidx, covar_include);
      }
    }
    uint32_t covar_ct = PopcountWords(covar_include, raw_covar_ctl);
    if (covar_ct != initial_covar_ct) {
      // 2. If any covariates were removed, redetermine the samples for which
      //    all phenotype/covariate values are present.  If there are more than
      //    before, add back any provisionally removed covariates which don't
      //    have any missing values in the new sample set, and are now
      //    nonconstant.
      //    Categorical covariates should behave just like they had been
      //    pre-split into n-1 0/1 indicator variables.
      memcpy(cur_sample_include, sample_include_backup, raw_sample_ctl * sizeof(intptr_t));
      covar_uidx_base = 0;
      cur_bits = covar_include[0];
      for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
        const uintptr_t covar_uidx = BitIter1(covar_include, &covar_uidx_base, &cur_bits);
        if (covar_cols[covar_uidx].nonmiss) {
          BitvecAnd(covar_cols[covar_uidx].nonmiss, raw_sample_ctl, cur_sample_include);
        }
      }
      uint32_t new_sample_ct = PopcountWords(cur_sample_include, raw_sample_ctl);
      if (new_sample_ct > prev_sample_ct) {
        prev_sample_ct = new_sample_ct;
        covar_uidx_base = 0;
        cur_bits = initial_covar_include[0];
        for (uint32_t covar_idx = 0; covar_idx != initial_covar_ct; ++covar_idx) {
          const uintptr_t covar_uidx = BitIter1(initial_covar_include, &covar_uidx_base, &cur_bits);
          if (!IsSet(covar_include, covar_uidx)) {
            const PhenoCol* cur_covar_col = &(covar_cols[covar_uidx]);
            if (PopcountWordsIntersect(cur_sample_include, cur_covar_col->nonmiss, raw_sample_ctl) == prev_sample_ct) {
              if (!IsConstCovar(cur_covar_col, cur_sample_include, prev_sample_ct)) {
                SetBit(covar_uidx, covar_include);
              }
            }
          }
        }
        covar_ct = PopcountWords(covar_include, raw_covar_ctl);
      }
    }
    if (!covar_ct) {
      *sample_ct_ptr = prev_sample_ct;
      goto GlmDetermineCovars_ret_NOCOVAR;
    }
    if (prev_sample_ct < 3) {
      goto GlmDetermineCovars_ret_SKIP;
    }
    // 3. if quantitative trait, remove samples corresponding to single-element
    //    categories or constant-except-for-one-sample regular covariates.
    uint32_t sample_ct = prev_sample_ct;
    uint32_t extra_cat_ct = 0;
    BigstackReset(sample_include_backup);
    if (!pheno_cc) {
      uintptr_t* cat_one_obs = nullptr;
      uintptr_t* cat_two_or_more_obs = nullptr;
      uint32_t* cat_first_sample_uidxs = nullptr;
      if (covar_max_nonnull_cat_ct) {
        const uint32_t max_cat_ctl = 1 + (covar_max_nonnull_cat_ct / kBitsPerWord);
        if (unlikely(
                bigstack_alloc_w(max_cat_ctl, &cat_one_obs) ||
                bigstack_alloc_w(max_cat_ctl, &cat_two_or_more_obs) ||
                bigstack_alloc_u32(covar_max_nonnull_cat_ct + 1, &cat_first_sample_uidxs))) {
          goto GlmDetermineCovars_ret_NOMEM;
        }
      }
      do {
        prev_sample_ct = sample_ct;
        covar_uidx_base = 0;
        cur_bits = covar_include[0];
        extra_cat_ct = 0;
        for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
          const uintptr_t covar_uidx = BitIter1(covar_include, &covar_uidx_base, &cur_bits);
          const PhenoCol* cur_covar_col = &(covar_cols[covar_uidx]);
          if (cur_covar_col->type_code == kPhenoDtypeOther) {
            continue;
          }
          // bugfix (31 Aug 2018): sample_uidx_base/cur_sample_include_bits
          // must be reset on each iteration
          uintptr_t sample_uidx_base = 0;
          uintptr_t cur_sample_include_bits = cur_sample_include[0];
          if (cur_covar_col->type_code == kPhenoDtypeQt) {
            const double* pheno_vals = cur_covar_col->data.qt;
            const uint32_t first_sample_uidx = BitIter1(cur_sample_include, &sample_uidx_base, &cur_sample_include_bits);
            double common_pheno_val = pheno_vals[first_sample_uidx];
            const uint32_t second_sample_uidx = BitIter1(cur_sample_include, &sample_uidx_base, &cur_sample_include_bits);
            const double second_pheno_val = pheno_vals[second_sample_uidx];
            uint32_t sample_idx = 2;
            uint32_t sample_uidx_remove;
            if (second_pheno_val != common_pheno_val) {
              sample_uidx_remove = second_sample_uidx;
              const uint32_t third_sample_uidx = BitIter1(cur_sample_include, &sample_uidx_base, &cur_sample_include_bits);
              const double third_pheno_val = pheno_vals[third_sample_uidx];
              if (third_pheno_val == second_pheno_val) {
                common_pheno_val = second_pheno_val;
                sample_uidx_remove = first_sample_uidx;
              } else if (third_pheno_val != common_pheno_val) {
                continue;
              }
              sample_idx = 3;
            } else {
              sample_uidx_remove = UINT32_MAX;
            }
            // minor bugfix (8 Jul 2018): forgot to increment sample_uidx on
            // first loop iteration here
            // sample_ct >= 3 guaranteed
            for (; sample_idx != sample_ct; ++sample_idx) {
              const uint32_t sample_uidx = BitIter1(cur_sample_include, &sample_uidx_base, &cur_sample_include_bits);
              if (pheno_vals[sample_uidx] != common_pheno_val) {
                if (sample_uidx_remove == UINT32_MAX) {
                  sample_uidx_remove = sample_uidx;
                } else {
                  break;
                }
              }
            }
            if (sample_idx == sample_ct) {
              if (sample_uidx_remove != UINT32_MAX) {
                if (--sample_ct == 2) {
                  goto GlmDetermineCovars_ret_SKIP;
                }
                ClearBit(sample_uidx_remove, cur_sample_include);
              } else {
                // constant covariate, remove it
                ClearBit(covar_uidx, covar_include);
              }
            }
          } else {
            const uint32_t cur_nonnull_cat_ct = cur_covar_col->nonnull_category_ct;
            const uint32_t cur_cat_ctl = 1 + (cur_nonnull_cat_ct / kBitsPerWord);
            ZeroWArr(cur_cat_ctl, cat_one_obs);
            ZeroWArr(cur_cat_ctl, cat_two_or_more_obs);
            const uint32_t* pheno_vals = cur_covar_col->data.cat;
            for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
              const uint32_t sample_uidx = BitIter1(cur_sample_include, &sample_uidx_base, &cur_sample_include_bits);
              const uint32_t cur_cat_idx = pheno_vals[sample_uidx];
              if (!IsSet(cat_two_or_more_obs, cur_cat_idx)) {
                if (IsSet(cat_one_obs, cur_cat_idx)) {
                  SetBit(cur_cat_idx, cat_two_or_more_obs);
                } else {
                  SetBit(cur_cat_idx, cat_one_obs);
                  cat_first_sample_uidxs[cur_cat_idx] = sample_uidx;
                }
              }
            }
            for (uint32_t widx = 0; widx != cur_cat_ctl; ++widx) {
              uintptr_t cur_word = cat_one_obs[widx] & (~cat_two_or_more_obs[widx]);
              if (cur_word) {
                const uint32_t* cat_first_sample_uidxs_iter = &(cat_first_sample_uidxs[widx * kBitsPerWord]);
                do {
                  const uint32_t cat_idx_lowbits = ctzw(cur_word);

                  // this could be moved out of the loop, but unlikely to be
                  // profitable...
                  --sample_ct;

                  ClearBit(cat_first_sample_uidxs_iter[cat_idx_lowbits], cur_sample_include);
                  cur_word &= cur_word - 1;
                } while (cur_word);
              }
            }
            if (sample_ct < 3) {
              goto GlmDetermineCovars_ret_SKIP;
            }
            uint32_t remaining_cat_ct = PopcountWords(cat_two_or_more_obs, cur_cat_ctl);
            if (remaining_cat_ct <= 1) {
              // now a constant covariate, remove it
              ClearBit(covar_uidx, covar_include);
            } else {
              extra_cat_ct += remaining_cat_ct - 2;
            }
          }
        }
        covar_ct = PopcountWords(covar_include, raw_covar_ctl);
      } while (sample_ct < prev_sample_ct);
    } else {
      uintptr_t* cat_covar_wkspace;
      if (unlikely(bigstack_alloc_w(1 + (covar_max_nonnull_cat_ct / kBitsPerWordD2), &cat_covar_wkspace))) {
        goto GlmDetermineCovars_ret_NOMEM;
      }
      if (!is_always_firth) {
        covar_uidx_base = 0;
        cur_bits = covar_include[0];
        if (!is_sometimes_firth) {
          do {
            prev_sample_ct = sample_ct;
            for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
              const uintptr_t covar_uidx = BitIter1(covar_include, &covar_uidx_base, &cur_bits);
              if (CheckForAndHandleSeparatedCovar(pheno_cc, covar_cols, raw_sample_ctl, covar_uidx, cur_sample_include, covar_include, &sample_ct, cat_covar_wkspace)) {
                *separation_found_ptr = 1;
                goto GlmDetermineCovars_ret_SKIP;
              }
            }
            covar_ct = PopcountWords(covar_include, raw_covar_ctl);
          } while (sample_ct < prev_sample_ct);
        } else {
          for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
            const uintptr_t covar_uidx = BitIter1(covar_include, &covar_uidx_base, &cur_bits);
            if (CheckForSeparatedCovar(pheno_cc, cur_sample_include, &(covar_cols[covar_uidx]), sample_ct, cat_covar_wkspace)) {
              *separation_found_ptr = 1;
              break;
            }
          }
        }
      }

      // now count extra categories
      covar_uidx_base = 0;
      cur_bits = covar_include[0];
      extra_cat_ct = 0;
      for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
        const uintptr_t covar_uidx = BitIter1(covar_include, &covar_uidx_base, &cur_bits);
        const PhenoCol* cur_covar_col = &(covar_cols[covar_uidx]);
        if (cur_covar_col->type_code == kPhenoDtypeCat) {
          const uint32_t remaining_cat_ct = IdentifyRemainingCats(cur_sample_include, cur_covar_col, sample_ct, cat_covar_wkspace);
          if (remaining_cat_ct > 2) {
            extra_cat_ct += remaining_cat_ct - 2;
          }
        }
      }
    }
    *sample_ct_ptr = sample_ct;
    *covar_ct_ptr = covar_ct;
    *extra_cat_ct_ptr = extra_cat_ct;
  }
  while (0) {
  GlmDetermineCovars_ret_NOMEM:
    reterr = 1;
    break;
  GlmDetermineCovars_ret_SKIP:
    *sample_ct_ptr = 0;
  GlmDetermineCovars_ret_NOCOVAR:
    *covar_ct_ptr = 0;
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

uint32_t CollapseParamOrTestSubset(const uintptr_t* covar_include, const uintptr_t* raw_params_or_tests, uint32_t domdev_present, uint32_t raw_covar_ct, uint32_t covar_ct, uint32_t add_interactions, uintptr_t* new_params_or_tests) {
  const uint32_t first_covar_pred_idx = 2 + domdev_present;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  uint32_t first_interaction_pred_read_idx = 0;
  uint32_t first_interaction_pred_write_idx = 0;
  if (add_interactions) {
    first_interaction_pred_read_idx = first_covar_pred_idx + raw_covar_ct;
    first_interaction_pred_write_idx = first_covar_pred_idx + covar_ct;
  }
  const uint32_t write_idx_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
  const uint32_t write_idx_ctl = BitCtToWordCt(write_idx_ct);
  ZeroWArr(write_idx_ctl, new_params_or_tests);
  // intercept, additive, domdev
  new_params_or_tests[0] = raw_params_or_tests[0] & (3 + 4 * domdev_present);
  // bugfix (26 Jun 2019): covar_include == nullptr when covar_ct == 0
  if (covar_ct) {
    uintptr_t covar_uidx_base = 0;
    uintptr_t cur_bits = covar_include[0];
    for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
      const uintptr_t covar_uidx = BitIter1(covar_include, &covar_uidx_base, &cur_bits);
      if (IsSet(raw_params_or_tests, first_covar_pred_idx + covar_uidx)) {
        SetBit(first_covar_pred_idx + covar_idx, new_params_or_tests);
      }
      if (add_interactions) {
        if (IsSet(raw_params_or_tests, first_interaction_pred_read_idx + domdev_present_p1 * covar_uidx)) {
          SetBit(first_interaction_pred_write_idx + domdev_present_p1 * covar_idx, new_params_or_tests);
        }
        if (domdev_present) {
          if (IsSet(raw_params_or_tests, first_interaction_pred_read_idx + 2 * covar_uidx + 1)) {
            SetBit(first_interaction_pred_write_idx + 2 * covar_idx + 1, new_params_or_tests);
          }
        }
      }
    }
  }
  return PopcountWords(new_params_or_tests, write_idx_ctl);
}


// may want to try to identify a specific linear dependency in rank-deficient
// case
ENUM_U31_DEF_START()
  kGlmErrcodeNone,
  kGlmErrcodeSampleCtLtePredictorCt,
  kGlmErrcodeConstOmittedAllele,
  kGlmErrcodeConstAllele,
  kGlmErrcodeCorrTooHigh, // 2 predictor args
  kGlmErrcodeVifInfinite,
  kGlmErrcodeVifTooHigh, // 1 predictor arg
  kGlmErrcodeSeparation, // 1 allele arg
  kGlmErrcodeRankDeficient,
  kGlmErrcodeLogisticConvergeFail,
  kGlmErrcodeFirthConvergeFail,
  kGlmErrcodeInvalidResult,
  // no codes for logistic-unfinished and firth-unfinished for now since we
  // still report results there

  // only during initial scan
  kGlmErrcodeUnstableScale
ENUM_U31_DEF_END(GlmErrcode);

#if __cplusplus >= 201103L
// see IntErr in plink2_base.h
struct GlmErr {
  GlmErr() {}

  GlmErr(uint64_t source) : value_(source) {}

  explicit operator uint64_t() const {
    return value_;
  }

  explicit operator uint32_t() const {
    return static_cast<uint32_t>(value_);
  }

  explicit operator bool() const {
    return (value_ != 0);
  }

private:
  uint64_t value_;
};
#else
typedef uint64_t GlmErr;
#endif

// The error code, along with up to two predictor/allele-index arguments, must
// fit in 8 bytes for now.  Since kMaxPhenoCt is larger than 2^16, we need a
// custom encoding.
// Current encoding has GlmErrcode in bits 0..7, the first argument in bits
// 8..31, the second argument in bits 32..55, and the top 8 bits are guaranteed
// to be zero.
static inline GlmErr SetGlmErr0(GlmErrcode errcode) {
  return errcode;
}

static inline GlmErr SetGlmErr1(GlmErrcode errcode, uint32_t arg1) {
  return errcode | (arg1 << 8);
}

static inline GlmErr SetGlmErr2(GlmErrcode errcode, uint32_t arg1, uint32_t arg2) {
  return S_CAST(uint64_t, errcode) | (arg1 << 8) | (S_CAST(uint64_t, arg2) << 32);
}

static inline GlmErrcode GetGlmErrCode(GlmErr glm_err) {
  return S_CAST(GlmErrcode, S_CAST(uint64_t, glm_err) & 255);
}

static inline uint32_t GetGlmErrArg1(GlmErr glm_err) {
  return S_CAST(uint32_t, glm_err) >> 8;
}

static inline uint32_t GetGlmErrArg2(GlmErr glm_err) {
  return S_CAST(uint64_t, glm_err) >> 32;
}

static const char kGlmErrcodeStrs[][24] = {"", "SAMPLE_CT<=PREDICTOR_CT", "CONST_OMITTED_ALLELE", "CONST_ALLELE", "CORR_TOO_HIGH", "VIF_INFINITE", "VIF_TOO_HIGH", "SEPARATION", "RANK_DEFICIENT", "LOGISTIC_CONVERGE_FAIL", "FIRTH_CONVERGE_FAIL", "INVALID_RESULT"};

char* AppendGlmErrstr(GlmErr glm_err, char* write_iter) {
  // todo: support predictor args
  GlmErrcode errcode = GetGlmErrCode(glm_err);
  write_iter = strcpya(write_iter, kGlmErrcodeStrs[errcode]);
  if (errcode == kGlmErrcodeSeparation) {
    *write_iter++ = ',';
    const uint32_t allele_idx = GetGlmErrArg1(glm_err);
    if (!allele_idx) {
      write_iter = strcpya_k(write_iter, "REF");
    } else {
      write_iter = strcpya_k(write_iter, "ALT");
      write_iter = u32toa(allele_idx, write_iter);
    }
  }
  return write_iter;
}

// * first_predictor_idx should be 1 if first term is constant-1 intercept, 0
//   otherwise
// * dbl_2d_buf[n] expected to be sum of predictor row n on input, is destroyed
// * if corr_buf is not nullptr, lower left is filled with uninverted
//   correlation matrix on exit
// * lower left of inverse_corr_buf is filled on return
GlmErr CheckMaxCorrAndVif(const double* predictor_dotprods, uint32_t first_predictor_idx, uint32_t predictor_ct, uintptr_t sample_ct, double max_corr, double vif_thresh, double* dbl_2d_buf, double* corr_buf, double* inverse_corr_buf, MatrixInvertBuf1* matrix_invert_buf1) {
  // we have dot products, now determine
  //   (dotprod - sum(a)mean(b)) / (N-1)
  // to get small-sample covariance
  const uintptr_t relevant_predictor_ct = predictor_ct - first_predictor_idx;
  if (relevant_predictor_ct == 1) {
    // bugfix (31 Jul 2019): precomputed images are wrong if these aren't
    // initialized
    if (corr_buf) {
      corr_buf[0] = 1.0;
      corr_buf[1] = 1.0;
    }
    inverse_corr_buf[0] = 1.0;
    return 0;
  }
  const uintptr_t relevant_predictor_ct_p1 = relevant_predictor_ct + 1;
  const double sample_ct_recip = 1.0 / u31tod(sample_ct);
  const double sample_ct_m1_d = u31tod(sample_ct - 1);
  const double sample_ct_m1_recip = 1.0 / sample_ct_m1_d;
  for (uintptr_t pred_idx1 = 0; pred_idx1 != relevant_predictor_ct; ++pred_idx1) {
    double* sample_corr_row = &(inverse_corr_buf[pred_idx1 * relevant_predictor_ct]);
    const uintptr_t input_pred_idx1 = pred_idx1 + first_predictor_idx;
    const double* predictor_dotprods_row = &(predictor_dotprods[input_pred_idx1 * predictor_ct]);
    const double covar1_mean_adj = dbl_2d_buf[input_pred_idx1] * sample_ct_recip;
    for (uintptr_t pred_idx2 = 0; pred_idx2 <= pred_idx1; ++pred_idx2) {
      const uintptr_t input_pred_idx2 = pred_idx2 + first_predictor_idx;
      sample_corr_row[pred_idx2] = (predictor_dotprods_row[input_pred_idx2] - covar1_mean_adj * dbl_2d_buf[input_pred_idx2]) * sample_ct_m1_recip;
    }
  }
  // now use dbl_2d_buf to store inverse-sqrts, to get to correlation matrix
  for (uintptr_t pred_idx = 0; pred_idx != relevant_predictor_ct; ++pred_idx) {
    dbl_2d_buf[pred_idx] = 1.0 / sqrt(inverse_corr_buf[pred_idx * relevant_predictor_ct_p1]);
  }
  // invert_symmdef_matrix only cares about bottom left of inverse_corr_buf[]
  for (uintptr_t pred_idx1 = 1; pred_idx1 != relevant_predictor_ct; ++pred_idx1) {
    const double inverse_stdev1 = dbl_2d_buf[pred_idx1];
    double* corr_row_iter = &(inverse_corr_buf[pred_idx1 * relevant_predictor_ct]);
    const double* inverse_stdev2_iter = dbl_2d_buf;
    for (uintptr_t pred_idx2 = 0; pred_idx2 != pred_idx1; ++pred_idx2) {
      const double cur_corr = (*corr_row_iter) * inverse_stdev1 * (*inverse_stdev2_iter++);
      // bugfix (14 Sep 2017): need to take absolute value here
      if (fabs(cur_corr) > max_corr) {
        return SetGlmErr2(kGlmErrcodeCorrTooHigh, pred_idx2, pred_idx1);
      }
      *corr_row_iter++ = cur_corr;
    }
  }
  for (uintptr_t pred_idx = 0; pred_idx != relevant_predictor_ct; ++pred_idx) {
    inverse_corr_buf[pred_idx * relevant_predictor_ct_p1] = 1.0;
  }
  if (corr_buf) {
    memcpy(corr_buf, inverse_corr_buf, relevant_predictor_ct * relevant_predictor_ct * sizeof(double));
    memcpy(&(corr_buf[relevant_predictor_ct * relevant_predictor_ct]), dbl_2d_buf, relevant_predictor_ct * sizeof(double));
  }
  if (InvertSymmdefMatrixChecked(relevant_predictor_ct, inverse_corr_buf, matrix_invert_buf1, dbl_2d_buf)) {
    return SetGlmErr0(kGlmErrcodeVifInfinite);
  }
  // VIFs = diagonal elements of inverse correlation matrix
  for (uintptr_t pred_idx = 0; pred_idx != relevant_predictor_ct; ++pred_idx) {
    if (inverse_corr_buf[pred_idx * relevant_predictor_ct_p1] > vif_thresh) {
      return SetGlmErr1(kGlmErrcodeVifTooHigh, pred_idx);
    }
  }
  return 0;
}

void PrintPrescanErrmsg(const char* domain_str, const char* pheno_name, const char** covar_names, GlmErr glm_err, uint32_t local_covar_ct, uint32_t skip_modifier, uint32_t is_batch) {
  const GlmErrcode errcode = GetGlmErrCode(glm_err);
  const char* msg_start = skip_modifier? "Note: Skipping " : "Error: Cannot proceed with ";
  if (errcode == kGlmErrcodeUnstableScale) {
    snprintf(g_logbuf, kLogbufSize, "%s%s--glm regression on phenotype '%s'%s, since genotype/covariate scales vary too widely for numerical stability of the current implementation. Try rescaling your covariates with e.g. --covar-variance-standardize.\n", msg_start, domain_str, pheno_name, is_batch? ", and other(s) with identical missingness patterns" : "");
  } else if (errcode == kGlmErrcodeVifInfinite) {
    snprintf(g_logbuf, kLogbufSize, "%s%s--glm regression on phenotype '%s'%s, since covariate correlation matrix could not be inverted (VIF_INFINITE).%s\n", msg_start, domain_str, pheno_name, is_batch? ", and other(s) with identical missingness patterns" : "", domain_str[0]? "" : " You may want to remove redundant covariates and try again.");
  } else {
    const char* covar_name1 = covar_names[GetGlmErrArg1(glm_err) + local_covar_ct];
    if (errcode == kGlmErrcodeVifTooHigh) {
      snprintf(g_logbuf, kLogbufSize, "%s%s--glm regression on phenotype '%s'%s, since variance inflation factor for covariate '%s' is too high (VIF_TOO_HIGH).%s\n", msg_start, domain_str, pheno_name, is_batch? ", and other(s) with identical missingness patterns" : "", covar_name1, domain_str[0]? "" : " You may want to remove redundant covariates and try again.");
    } else {
      const char* covar_name2 = covar_names[GetGlmErrArg2(glm_err) + local_covar_ct];
      snprintf(g_logbuf, kLogbufSize, "%s%s--glm regression on phenotype '%s'%s, since correlation between covariates '%s' and '%s' is too high (CORR_TOO_HIGH).%s\n", msg_start, domain_str, pheno_name, is_batch? ", and other(s) with identical missingness patterns" : "", covar_name1, covar_name2, domain_str[0]? "" : " You may want to remove redundant covariates and try again.");
    }
  }
  WordWrapB(0);
  if (skip_modifier) {
    logputsb();
  } else {
    logerrputsb();
  }
}

// no-missing-genotype optimizations:
// * most of the inter-predictor correlation matrix can be initialized once
//   from an image, doesn't even need to be refreshed; same goes for
//   inv_corr_sqrts
// * inverse of that matrix can be precomputed, and used in rank 1 inverse
//   update (rank 2 for domdev case)
//
// Other notes:
// * dbl_2d_buf does not need to store row sums, since we can get that from
//   predictor_dotprods (first actual, though "irrelevant", predictor is
//   all-1).
// * predictor_ct includes intercept.
// * geno_pred_ct must be 1 or 2.
// * ainv_b_buf[] must have size at least 4 * nongeno_pred_ct.
GlmErr CheckMaxCorrAndVifNm(const double* predictor_dotprods, const double* corr_inv, uint32_t predictor_ct, uint32_t geno_pred_ct, double sample_ct_recip, double sample_ct_m1_recip, double max_corr, double vif_thresh, double* __restrict semicomputed_corr_matrix, double* __restrict semicomputed_inv_corr_sqrts, double* __restrict corr_row_buf, double* __restrict inverse_corr_diag, double* __restrict ainv_b_buf) {
  // we have dot products, now determine
  //   (dotprod - sum(a)mean(b)) / (N-1)
  // to get small-sample covariance
  if (predictor_ct == 2) {
    return 0;
  }
  const uintptr_t relevant_predictor_ct = predictor_ct - 1;
  const uintptr_t relevant_predictor_ct_p1 = predictor_ct;
  // predictor_dotprods[] *rows* 1 (and 2, if geno_pred_ct == 2) are filled,
  // rather than columns
  {
    const double covar1_mean_adj = predictor_dotprods[predictor_ct] * sample_ct_recip;
    semicomputed_corr_matrix[0] = (predictor_dotprods[predictor_ct + 1] - covar1_mean_adj * predictor_dotprods[predictor_ct]) * sample_ct_m1_recip;
  }
  for (uintptr_t pred_idx1 = 1; pred_idx1 != relevant_predictor_ct; ++pred_idx1) {
    double* sample_corr_row = &(semicomputed_corr_matrix[pred_idx1 * relevant_predictor_ct]);
    const uintptr_t input_pred_idx1 = pred_idx1 + 1;
    const double* predictor_dotprods_row = &(predictor_dotprods[input_pred_idx1 * predictor_ct]);
    const double covar1_mean_adj = predictor_dotprods_row[0] * sample_ct_recip;
    uintptr_t pred_idx2 = 0;
    for (; pred_idx2 != geno_pred_ct; ++pred_idx2) {
      const uintptr_t input_pred_idx2 = pred_idx2 + 1;
      sample_corr_row[pred_idx2] = (predictor_dotprods[input_pred_idx2 * predictor_ct + input_pred_idx1] - covar1_mean_adj * predictor_dotprods[input_pred_idx2 * predictor_ct]) * sample_ct_m1_recip;
    }
    // document whether pred_idx2 guaranteed to be <= input_pred_idx1 if this
    // code is revisited
    for (; pred_idx2 < input_pred_idx1; ++pred_idx2) {
      const uintptr_t input_pred_idx2 = pred_idx2 + 1;
      sample_corr_row[pred_idx2] = (predictor_dotprods_row[input_pred_idx2] - covar1_mean_adj * predictor_dotprods[input_pred_idx2 * predictor_ct]) * sample_ct_m1_recip;
    }
  }
  // assumes semicomputed_inv_corr_sqrts[geno_pred_ct..] is pre-initialized
  semicomputed_inv_corr_sqrts[0] = 1.0 / sqrt(semicomputed_corr_matrix[0]);
  const uintptr_t nongeno_pred_ct = relevant_predictor_ct - geno_pred_ct;
  double inverse_stdev1 = semicomputed_inv_corr_sqrts[0];
  const double* nongeno_inverse_stdevs = &(semicomputed_inv_corr_sqrts[geno_pred_ct]);
  const double* corr_col = &(semicomputed_corr_matrix[geno_pred_ct * relevant_predictor_ct]);
  for (uintptr_t nongeno_pred_idx = 0; nongeno_pred_idx != nongeno_pred_ct; ++nongeno_pred_idx) {
    const double inverse_stdev2 = nongeno_inverse_stdevs[nongeno_pred_idx];
    const double cur_corr = inverse_stdev1 * inverse_stdev2 * corr_col[nongeno_pred_idx * relevant_predictor_ct];
    if (fabs(cur_corr) > max_corr) {
      return SetGlmErr2(kGlmErrcodeCorrTooHigh, 0, geno_pred_ct + nongeno_pred_idx);
    }
    corr_row_buf[nongeno_pred_idx] = cur_corr;
  }
  if (geno_pred_ct == 1) {
    if (InvertRank1SymmDiag(corr_inv, corr_row_buf, nongeno_pred_ct, 1.0, inverse_corr_diag, ainv_b_buf)) {
      return SetGlmErr0(kGlmErrcodeVifInfinite);
    }
  } else {
    inverse_stdev1 = 1.0 / sqrt(semicomputed_corr_matrix[relevant_predictor_ct_p1]);
    corr_col = &(semicomputed_corr_matrix[geno_pred_ct * relevant_predictor_ct + 1]);
    double* corr_row2 = &(corr_row_buf[nongeno_pred_ct]);
    for (uintptr_t nongeno_pred_idx = 0; nongeno_pred_idx != nongeno_pred_ct; ++nongeno_pred_idx) {
      const double inverse_stdev2 = nongeno_inverse_stdevs[nongeno_pred_idx];
      const double cur_corr = inverse_stdev1 * inverse_stdev2 * corr_col[nongeno_pred_idx * relevant_predictor_ct];
      if (fabs(cur_corr) > max_corr) {
        return SetGlmErr2(kGlmErrcodeCorrTooHigh, 1, 2 + nongeno_pred_idx);
      }
      corr_row2[nongeno_pred_idx] = cur_corr;
    }
    const double inverse_stdev2 = semicomputed_inv_corr_sqrts[0];
    const double cur_corr = inverse_stdev1 * inverse_stdev2 * semicomputed_corr_matrix[relevant_predictor_ct];
    if (fabs(cur_corr) > max_corr) {
      return SetGlmErr2(kGlmErrcodeCorrTooHigh, 0, 1);
    }
    // do we want special handling of nongeno_pred_ct == 0?
    if (InvertRank2SymmDiag(corr_inv, corr_row_buf, nongeno_pred_ct, 1.0, cur_corr, 1.0, inverse_corr_diag, ainv_b_buf, &(ainv_b_buf[2 * nongeno_pred_ct]))) {
      return SetGlmErr0(kGlmErrcodeVifInfinite);
    }
  }
  // VIFs = diagonal elements of inverse correlation matrix
  for (uintptr_t pred_idx = 0; pred_idx != relevant_predictor_ct; ++pred_idx) {
    if (inverse_corr_diag[pred_idx] > vif_thresh) {
      return SetGlmErr1(kGlmErrcodeVifTooHigh, pred_idx);
    }
  }
  return 0;
}

// Only called by GlmLogisticThread(), so there are the following differences
// from CheckMaxCorrAndVif():
// * predictor_dotprods is not precomputed; we start with predictors_pmaj
//   instead.
// * predictors_pmaj already has the intercept stripped off, so we don't need
//   relevant_predictor_ct := predictor_ct - 1, etc.
// * sample_stride parameter added, since predictors_pmaj has vector-aligned
//   rather than packed rows.
// * dbl_2d_buf not assumed to be filled with row sums, we compute them here.
//   (probably want to modify CheckMaxCorrAndVif() to do the same.)
// This now uses double-precision arithmetic since matrix inversion is too
// inconsistent if we stick to single-precision.
GlmErr CheckMaxCorrAndVifF(const float* predictors_pmaj, uint32_t predictor_ct, uint32_t sample_ct, uint32_t sample_stride, double max_corr, double vif_thresh, float* predictor_dotprod_buf, double* dbl_2d_buf, double* inverse_corr_buf, MatrixInvertBuf1* inv_1d_buf) {
  MultiplySelfTransposeStridedF(predictors_pmaj, predictor_ct, sample_ct, sample_stride, predictor_dotprod_buf);
  for (uintptr_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    const float* predictor_row = &(predictors_pmaj[pred_idx * sample_stride]);
    double row_sum = 0.0;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      row_sum += S_CAST(double, predictor_row[sample_idx]);
    }
    dbl_2d_buf[pred_idx] = row_sum;
  }
  const uint32_t predictor_ct_p1 = predictor_ct + 1;
  const double sample_ct_recip = 1.0 / u31tod(sample_ct);
  const double sample_ct_m1_d = u31tod(sample_ct - 1);
  const double sample_ct_m1_recip = 1.0 / sample_ct_m1_d;
  for (uint32_t pred_idx1 = 0; pred_idx1 != predictor_ct; ++pred_idx1) {
    double* sample_corr_row = &(inverse_corr_buf[pred_idx1 * predictor_ct]);
    const float* predictor_dotprod_row = &(predictor_dotprod_buf[pred_idx1 * predictor_ct]);
    const double covar1_mean_adj = dbl_2d_buf[pred_idx1] * sample_ct_recip;
    for (uint32_t pred_idx2 = 0; pred_idx2 <= pred_idx1; ++pred_idx2) {
      sample_corr_row[pred_idx2] = (S_CAST(double, predictor_dotprod_row[pred_idx2]) - covar1_mean_adj * dbl_2d_buf[pred_idx2]) * sample_ct_m1_recip;
    }
  }
  // now use dbl_2d_buf to store inverse-sqrts, to get to correlation matrix
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    dbl_2d_buf[pred_idx] = 1.0 / sqrt(inverse_corr_buf[pred_idx * predictor_ct_p1]);
  }
  // invert_symmdef_matrix only cares about bottom left of inverse_corr_buf[]
  for (uint32_t pred_idx1 = 1; pred_idx1 != predictor_ct; ++pred_idx1) {
    const double inverse_stdev1 = dbl_2d_buf[pred_idx1];
    double* corr_row_iter = &(inverse_corr_buf[pred_idx1 * predictor_ct]);
    const double* inverse_stdev2_iter = dbl_2d_buf;
    for (uintptr_t pred_idx2 = 0; pred_idx2 != pred_idx1; ++pred_idx2) {
      const double cur_corr = (*corr_row_iter) * inverse_stdev1 * (*inverse_stdev2_iter++);
      if (fabs(cur_corr) > max_corr) {
        return SetGlmErr2(kGlmErrcodeCorrTooHigh, pred_idx2, pred_idx1);
      }
      *corr_row_iter++ = cur_corr;
    }
  }
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    inverse_corr_buf[pred_idx * predictor_ct_p1] = 1.0;
  }
  if (InvertSymmdefMatrixChecked(predictor_ct, inverse_corr_buf, inv_1d_buf, dbl_2d_buf)) {
    return SetGlmErr0(kGlmErrcodeVifInfinite);
  }

  // VIFs = diagonal elements of inverse correlation matrix
  for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
    if (inverse_corr_buf[pred_idx * predictor_ct_p1] > vif_thresh) {
      return SetGlmErr1(kGlmErrcodeVifTooHigh, pred_idx);
    }
  }
  return 0;
}

PglErr GlmFillAndTestCovars(const uintptr_t* sample_include, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, double* covar_dotprod, double* corr_buf, double* inverse_corr_buf, double* covars_cmaj, const char** cur_covar_names, GlmErr* glm_err_ptr) {
  *glm_err_ptr = 0;
  if (covar_ct == local_covar_ct) {
    // bugfix (5 Mar 2018): need to copy local-covar names
    for (uintptr_t local_covar_read_idx = 0; local_covar_read_idx != covar_ct; ++local_covar_read_idx) {
      cur_covar_names[local_covar_read_idx] = &(covar_names[local_covar_read_idx * max_covar_name_blen]);
    }
    return kPglRetSuccess;
  }
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  const uint32_t covar_max_cat_ctl = 1 + (covar_max_nonnull_cat_ct / kBitsPerWord);
  MatrixInvertBuf1* matrix_invert_buf1;
  uintptr_t* cat_covar_wkspace;
  double* dbl_2d_buf;
  if (unlikely(
          BIGSTACK_ALLOC_X(MatrixInvertBuf1, kMatrixInvertBuf1CheckedAlloc * new_nonlocal_covar_ct, &matrix_invert_buf1) ||
          bigstack_alloc_w(covar_max_cat_ctl, &cat_covar_wkspace) ||
          bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &dbl_2d_buf))) {
    return kPglRetNomem;
  }
  uint32_t* cat_obs_buf;
  // bugfix (11 May 2020): we were previously only initializing cat_obs_buf
  // when extra_cat_ct > 0; this resulted in segfaults when every categorical
  // covariate had only two categories.
  if (unlikely(
          bigstack_alloc_u32(covar_max_nonnull_cat_ct + 1, &cat_obs_buf))) {
    return kPglRetNomem;
  }
  unsigned char* alloc_base = g_bigstack_base;
  unsigned char* new_covar_name_alloc = g_bigstack_end;
  const uint32_t first_sample_uidx = AdvTo1Bit(sample_include, 0);
  uintptr_t covar_read_uidx_base = 0;
  uintptr_t covar_include_bits = covar_include[0];
  const char** cur_covar_names_iter = cur_covar_names;
  double* covar_write_iter = covars_cmaj;
  double* sum_iter = dbl_2d_buf;
  // Quasi-bugfix (22 Jan 2019): Main loop is too numerically unstable,
  // especially in the single-precision logistic/Firth regression case, when a
  // covariate is on a very different scale from the main genotype column (or
  // other covariates).  We now check sum-of-squares and (rescaled) sample
  // variance of each fixed covariate column, and error out when the range of
  // either value exceeds ~6 orders of magnitude.  (Conveniently, this catches
  // the fairly common "year of birth" problem.)
  // Initial ssq and variance values correspond to genotype columns with a MAF
  // range of [0.01, 0.5].
  double max_ssq = 1.5 * u31tod(sample_ct);
  double min_ssq_minus_sqmean = 0.0198 * u31tod(sample_ct - 1);
  for (uintptr_t covar_read_idx = 0; covar_read_idx != covar_ct; ++covar_read_idx) {
    const uintptr_t covar_read_uidx = BitIter1(covar_include, &covar_read_uidx_base, &covar_include_bits);
    const PhenoCol* cur_covar_col = &(covar_cols[covar_read_uidx]);
    const char* covar_name_base = &(covar_names[covar_read_uidx * max_covar_name_blen]);
    if (cur_covar_col->type_code == kPhenoDtypeOther) {
      // local covariate
      *cur_covar_names_iter++ = covar_name_base;
    } else if (cur_covar_col->type_code == kPhenoDtypeQt) {
      *cur_covar_names_iter++ = covar_name_base;
      const double* covar_vals = cur_covar_col->data.qt;
      uintptr_t sample_uidx_base;
      uintptr_t sample_include_bits;
      BitIter1Start(sample_include, first_sample_uidx, &sample_uidx_base, &sample_include_bits);
      double covar_sum = 0.0;
      double covar_ssq = 0.0;
      for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
        const double cur_covar_val = covar_vals[sample_uidx];
        covar_sum += cur_covar_val;
        covar_ssq += cur_covar_val * cur_covar_val;
        *covar_write_iter++ = cur_covar_val;
      }
      *sum_iter++ = covar_sum;
      if (covar_ssq > max_ssq) {
        max_ssq = covar_ssq;
      }
      covar_ssq -= covar_sum * covar_sum / u31tod(sample_ct);
      if (covar_ssq < min_ssq_minus_sqmean) {
        min_ssq_minus_sqmean = covar_ssq;
      }
    } else {
      // this is equivalent to "--split-cat-pheno omit-most covar-01"
      const uint32_t cur_cat_ct = cur_covar_col->nonnull_category_ct + 1;
      const uint32_t cur_cat_ctl = BitCtToWordCt(cur_cat_ct);
      const uint32_t largest_cat_uidx = IdentifyRemainingCatsAndMostCommon(sample_include, cur_covar_col, sample_ct, cat_covar_wkspace, cat_obs_buf);
      const uint32_t remaining_cat_ct = PopcountWords(cat_covar_wkspace, cur_cat_ctl);
      assert(remaining_cat_ct >= 2);
      ClearBit(largest_cat_uidx, cat_covar_wkspace);
      const uint32_t* covar_vals = cur_covar_col->data.cat;
      const char* const* cur_category_names = cur_covar_col->category_names;
      const uint32_t covar_name_base_slen = strlen(covar_name_base);
      uintptr_t cat_uidx_base;
      uintptr_t cat_covar_wkspace_bits;
      BitIter1Start(cat_covar_wkspace, 1, &cat_uidx_base, &cat_covar_wkspace_bits);
      for (uint32_t cat_idx = 1; cat_idx != remaining_cat_ct; ++cat_idx) {
        const uintptr_t cat_uidx = BitIter1(cat_covar_wkspace, &cat_uidx_base, &cat_covar_wkspace_bits);
        const char* catname = cur_category_names[cat_uidx];
        const uint32_t catname_slen = strlen(catname);
        new_covar_name_alloc -= covar_name_base_slen + catname_slen + 2;
        if (unlikely(new_covar_name_alloc < alloc_base)) {
          return kPglRetNomem;
        }
        char* new_covar_name_write = memcpyax(new_covar_name_alloc, covar_name_base, covar_name_base_slen, '=');
        memcpy(new_covar_name_write, catname, catname_slen + 1);
        *cur_covar_names_iter++ = R_CAST(const char*, new_covar_name_alloc);

        uintptr_t sample_uidx_base;
        uintptr_t sample_include_bits;
        BitIter1Start(sample_include, first_sample_uidx, &sample_uidx_base, &sample_include_bits);
        uint32_t cur_cat_obs_ct = 0;
        for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
          const uint32_t cur_sample_is_in_cat = (covar_vals[sample_uidx] == cat_uidx);
          cur_cat_obs_ct += cur_sample_is_in_cat;
          *covar_write_iter++ = u31tod(cur_sample_is_in_cat);
        }
        double covar_ssq = u31tod(cur_cat_obs_ct);
        *sum_iter++ = covar_ssq;
        if (covar_ssq > max_ssq) {
          max_ssq = covar_ssq;
        }
        covar_ssq -= covar_ssq * covar_ssq / u31tod(sample_ct);
        if (covar_ssq < min_ssq_minus_sqmean) {
          min_ssq_minus_sqmean = covar_ssq;
        }
      }
    }
  }
  BigstackEndSet(new_covar_name_alloc);
  if (min_ssq_minus_sqmean * 1048576.0 < max_ssq) {
    // probable todo: automatically variance-standardize, while keeping track
    // of the linear transformations so we can translate results back to
    // original units in the final output.
    *glm_err_ptr = SetGlmErr0(kGlmErrcodeUnstableScale);
    return kPglRetSkipped;
  }
  assert(covar_write_iter == &(covars_cmaj[new_nonlocal_covar_ct * sample_ct]));
  MultiplySelfTranspose(covars_cmaj, new_nonlocal_covar_ct, sample_ct, covar_dotprod);
  // intentionally ignore error code, since all callers check glm_err
  *glm_err_ptr = CheckMaxCorrAndVif(covar_dotprod, 0, new_nonlocal_covar_ct, sample_ct, max_corr, vif_thresh, dbl_2d_buf, corr_buf, inverse_corr_buf,  matrix_invert_buf1);
  BigstackReset(matrix_invert_buf1);
  return kPglRetSuccess;
}

// Useful precomputed values for linear and logistic regression, for variants
// with no missing genotypes.
typedef struct {
  double* xtx_image;  // (covar_ct + domdev_present + 2)^2, genotype cols empty
  double* covarx_dotprod_inv;  // (covar_ct + 1) x (covar_ct + 1), reflected
  double* corr_inv;  // covar_ct x covar_ct, reflected
  double* corr_image;  // (covar_ct + domdev_present_p1)^2, genotype cols empty
  double* corr_inv_sqrts;  // covar_ct x 1
  double* xt_y_image;  // (covar_ct + domdev_present + 2) x 1
} RegressionNmPrecomp;

BoolErr InitNmPrecomp(const double* covars_cmaj, const double* covar_dotprod, const double* corr_buf, const double* corr_inv_tri, uint32_t sample_ct, uint32_t is_qt, uintptr_t new_covar_ct, uintptr_t xtx_state, RegressionNmPrecomp* nm_precomp) {
  uintptr_t stride = new_covar_ct + 1 + xtx_state;
  double* xtx_image = nm_precomp->xtx_image;
  ZeroDArr(stride * stride, xtx_image);
  xtx_image[0] = u31tod(sample_ct);
  for (uintptr_t covar_idx = 0; covar_idx != new_covar_ct; ++covar_idx) {
    const double* cur_covar = &(covars_cmaj[covar_idx * sample_ct]);
    double dxx = 0.0;
    for (uint32_t uii = 0; uii != sample_ct; ++uii) {
      dxx += cur_covar[uii];
    }
    double* xtx_image_row = &(xtx_image[(covar_idx + 1 + xtx_state) * stride]);
    xtx_image_row[0] = dxx;
    memcpy(&(xtx_image_row[1 + xtx_state]), &(covar_dotprod[covar_idx * new_covar_ct]), (covar_idx + 1) * sizeof(double));
  }
  if (is_qt) {
    // also save precomputed inverse of covar dotprods, with intercept included
    // (not currently needed in logistic case)
    double* covarx_dotprod_inv = nm_precomp->covarx_dotprod_inv;
    covarx_dotprod_inv[0] = xtx_image[0];
    const uintptr_t new_covar_ct_p1 = new_covar_ct + 1;
    for (uintptr_t row_idx = 1; row_idx != new_covar_ct_p1; ++row_idx) {
      covarx_dotprod_inv[row_idx * new_covar_ct_p1] = xtx_image[(row_idx + xtx_state) * stride];
      memcpy(&(covarx_dotprod_inv[row_idx * new_covar_ct_p1 + 1]), &(covar_dotprod[(row_idx - 1) * new_covar_ct]), row_idx * sizeof(double));
    }
    // this makes assumptions about amount of space past corr_inv...
    if (InvertSymmdefMatrixChecked(new_covar_ct_p1, covarx_dotprod_inv, R_CAST(MatrixInvertBuf1*, &(nm_precomp->corr_inv[new_covar_ct_p1 * MAXV(new_covar_ct_p1, 3)])), nm_precomp->corr_inv)) {
      return 1;
    }
    ReflectMatrix(new_covar_ct_p1, covarx_dotprod_inv);
    // xt_y_image fill no longer happens in this function (since that's
    // inappropriate in the QT-batch case)
  }

  memcpy(nm_precomp->corr_inv, corr_inv_tri, new_covar_ct * new_covar_ct * sizeof(double));
  ReflectMatrix(new_covar_ct, nm_precomp->corr_inv);
  --stride;
  double* corr_image = nm_precomp->corr_image;
  // also store uninverted correlation matrix
  ZeroDArr(stride * stride, corr_image);
  for (uintptr_t orig_row_idx = 0; orig_row_idx != new_covar_ct; ++orig_row_idx) {
    memcpy(&(corr_image[(orig_row_idx + xtx_state) * stride + xtx_state]), &(corr_buf[orig_row_idx * new_covar_ct]), (orig_row_idx + 1) * sizeof(double));
  }
  // and store inverse-sqrts after the correlation matrix
  memcpy(nm_precomp->corr_inv_sqrts, &(corr_buf[new_covar_ct * new_covar_ct]), new_covar_ct * sizeof(double));
  return 0;
}

// xtx_state 0: either interactions or local covariates present, no xtx_image
// xtx_state 1: only additive effect
// xtx_state 2: additive and domdev effects
// Note that the immediate return value only corresponds to out-of-memory.
// Other errors are indicated by glm_err_ptr on a *false* return
// value, due to how the caller is expected to handle them.
BoolErr GlmAllocFillAndTestCovarsQt(const uintptr_t* sample_include, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, RegressionNmPrecomp** nm_precomp_ptr, double** covars_cmaj_d_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  if (unlikely(
          bigstack_alloc_kcp(new_covar_ct, cur_covar_names_ptr) ||
          bigstack_alloc_d(new_nonlocal_covar_ct * sample_ct, covars_cmaj_d_ptr))) {
    return 1;
  }
  double* corr_buf = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  *nm_precomp_ptr = nullptr;
  if (xtx_state) {
    assert(!local_covar_ct);
    if (unlikely(BIGSTACK_ALLOC_X(RegressionNmPrecomp, 1, nm_precomp_ptr))) {
      return 1;
    }
    // x^2 + (2*xtx_state + 2)x + (5*xtx_state - 1)
    // 2x^2 + (2*xtx_state + 4)x + (5*xtx_state)
    // 3x^2 + (2*xtx_state + 4)x + (5*xtx_state)
    // 4x^2 + (4*xtx_state + 4)x + (8*xtx_state - 2)
    // 4x^2 + (4*xtx_state + 5)x + (8*xtx_state - 2)
    //
    // xt_y_image: (1 + xtx_state + x) elements per phenotype in subbatch;
    //   now allocated later
    if (unlikely(
            bigstack_alloc_d(8 * xtx_state - 2 +
                             new_covar_ct * (5 + 4 * xtx_state + 4 * new_covar_ct),
                             &((*nm_precomp_ptr)->xtx_image)) ||
            bigstack_alloc_d(new_covar_ct * (new_covar_ct + 1), &corr_buf))) {
      return 1;
    }
    (*nm_precomp_ptr)->covarx_dotprod_inv = &((*nm_precomp_ptr)->xtx_image[(new_covar_ct + xtx_state + 1) * (new_covar_ct + xtx_state + 1)]);
    (*nm_precomp_ptr)->corr_inv = &((*nm_precomp_ptr)->covarx_dotprod_inv[(new_covar_ct + 1) * (new_covar_ct + 1)]);
    (*nm_precomp_ptr)->corr_image = &((*nm_precomp_ptr)->corr_inv[new_covar_ct * new_covar_ct]);
    (*nm_precomp_ptr)->corr_inv_sqrts = &((*nm_precomp_ptr)->corr_image[(new_covar_ct + xtx_state) * (new_covar_ct + xtx_state)]);
    (*nm_precomp_ptr)->xt_y_image = nullptr;  // defensive
    bigstack_mark = R_CAST(unsigned char*, corr_buf);
  }
  double* covar_dotprod;
  double* inverse_corr_buf;
  if (unlikely(
          bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &covar_dotprod) ||
          bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &inverse_corr_buf))) {
    return 1;
  }
  PglErr reterr = GlmFillAndTestCovars(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, max_corr, vif_thresh, covar_dotprod, corr_buf, inverse_corr_buf, *covars_cmaj_d_ptr, *cur_covar_names_ptr, glm_err_ptr);
  if (unlikely(reterr)) {
    return (reterr == kPglRetNomem);
  }
  if (xtx_state) {
    if (InitNmPrecomp(*covars_cmaj_d_ptr, covar_dotprod, corr_buf, inverse_corr_buf, sample_ct, 1, new_covar_ct, xtx_state, *nm_precomp_ptr)) {
      *glm_err_ptr = SetGlmErr0(kGlmErrcodeVifInfinite);
      return 0;  // not out-of-memory error
    }
  }
  BigstackReset(bigstack_mark);
  return 0;
}

void FillPhenoAndXtY(const uintptr_t* sample_include, const double* __restrict pheno_qt, const double* __restrict covars_cmaj_d, uintptr_t sample_ct, uintptr_t domdev_present_p1, uintptr_t covar_ct, double* xt_y_image, double* __restrict pheno_d) {
  double* pheno_d_iter = pheno_d;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  double pheno_sum = 0.0;
  for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    const double cur_pheno = pheno_qt[sample_uidx];
    *pheno_d_iter++ = cur_pheno;
    pheno_sum += cur_pheno;
  }
  if (!xt_y_image) {
    return;
  }
  xt_y_image[0] = pheno_sum;
  ZeroDArr(domdev_present_p1, &(xt_y_image[1]));
  ColMajorVectorMatrixMultiplyStrided(pheno_d, covars_cmaj_d, sample_ct, sample_ct, covar_ct, &(xt_y_image[1 + domdev_present_p1]));
}

BoolErr GlmAllocFillAndTestPhenoCovarsQt(const uintptr_t* sample_include, const double* pheno_qt, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, double** pheno_d_ptr, RegressionNmPrecomp** nm_precomp_ptr, double** covars_cmaj_d_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  if (unlikely(GlmAllocFillAndTestCovarsQt(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, max_corr, vif_thresh, xtx_state, nm_precomp_ptr, covars_cmaj_d_ptr, cur_covar_names_ptr, glm_err_ptr))) {
    return 1;
  }
  if (*glm_err_ptr) {
    // this is a bit messy
    return 0;
  }
  if (unlikely(
          bigstack_alloc_d(sample_ct, pheno_d_ptr))) {
    return 1;
  }
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  double* xt_y_image = nullptr;
  if (xtx_state) {
    if (unlikely(
            bigstack_alloc_d(1 + xtx_state + new_covar_ct, &xt_y_image))) {
      return 1;
    }
    (*nm_precomp_ptr)->xt_y_image = xt_y_image;
  }
  FillPhenoAndXtY(sample_include, pheno_qt, *covars_cmaj_d_ptr, sample_ct, xtx_state, new_covar_ct, xt_y_image, *pheno_d_ptr);
  return 0;
}

static const float kSmallFloats[4] = {0.0, 1.0, 2.0, 3.0};

BoolErr GlmAllocFillAndTestPhenoCovarsCc(const uintptr_t* sample_include, const uintptr_t* pheno_cc, const uintptr_t* covar_include, const PhenoCol* covar_cols, const char* covar_names, uintptr_t sample_ct, uintptr_t covar_ct, uint32_t local_covar_ct, uint32_t covar_max_nonnull_cat_ct, uintptr_t extra_cat_ct, uintptr_t max_covar_name_blen, double max_corr, double vif_thresh, uintptr_t xtx_state, uintptr_t** pheno_cc_collapsed_ptr, uintptr_t** gcount_case_interleaved_vec_ptr, float** pheno_f_ptr, RegressionNmPrecomp** nm_precomp_ptr, float** covars_cmaj_f_ptr, const char*** cur_covar_names_ptr, GlmErr* glm_err_ptr) {
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  const uintptr_t new_covar_ct = covar_ct + extra_cat_ct;
  const uintptr_t new_nonlocal_covar_ct = new_covar_ct - local_covar_ct;
  const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
  if (unlikely(
          bigstack_alloc_w(sample_ctv * kWordsPerVec, pheno_cc_collapsed_ptr) ||
          bigstack_alloc_f(sample_ctav, pheno_f_ptr) ||
          bigstack_alloc_f(new_nonlocal_covar_ct * sample_ctav, covars_cmaj_f_ptr) ||
          bigstack_alloc_kcp(new_covar_ct, cur_covar_names_ptr))) {
    return 1;
  }
  double* corr_buf = nullptr;
  unsigned char* bigstack_mark = g_bigstack_base;
  *nm_precomp_ptr = nullptr;
  if (xtx_state) {
    assert(!local_covar_ct);
    if (unlikely(BIGSTACK_ALLOC_X(RegressionNmPrecomp, 1, nm_precomp_ptr))) {
      return 1;
    }
    // x^2 + (2*xtx_state + 2)x + (5*xtx_state - 1)
    // 2x^2 + (2*xtx_state + 2)x + (5*xtx_state - 1)
    // 3x^2 + (4*xtx_state + 2)x + (8*xtx_state - 3)
    // 3x^2 + (4*xtx_state + 3)x + (8*xtx_state - 3)
    if (unlikely(
            bigstack_alloc_d(8 * xtx_state - 3 +
                             new_covar_ct * (3 + 4 * xtx_state + 3 * new_covar_ct),
                             &((*nm_precomp_ptr)->xtx_image)) ||
            bigstack_alloc_d(new_covar_ct * (new_covar_ct + 1), &corr_buf))) {
      return 1;
    }
    (*nm_precomp_ptr)->covarx_dotprod_inv = nullptr;
    (*nm_precomp_ptr)->corr_inv = &((*nm_precomp_ptr)->xtx_image[(new_covar_ct + xtx_state + 1) * (new_covar_ct + xtx_state + 1)]);
    (*nm_precomp_ptr)->corr_image = &((*nm_precomp_ptr)->corr_inv[new_covar_ct * new_covar_ct]);
    (*nm_precomp_ptr)->corr_inv_sqrts = &((*nm_precomp_ptr)->corr_image[(new_covar_ct + xtx_state) * (new_covar_ct + xtx_state)]);
    (*nm_precomp_ptr)->xt_y_image = nullptr;
    bigstack_mark = R_CAST(unsigned char*, corr_buf);
  }
  double* covars_cmaj_d;
  double* covar_dotprod;
  double* inverse_corr_buf;
  if (unlikely(
          bigstack_alloc_d(new_nonlocal_covar_ct * sample_ct, &covars_cmaj_d) ||
          bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &covar_dotprod) ||
          bigstack_alloc_d(new_nonlocal_covar_ct * new_nonlocal_covar_ct, &inverse_corr_buf))) {
    return 1;
  }
  uintptr_t* pheno_cc_collapsed = *pheno_cc_collapsed_ptr;
  CopyBitarrSubset(pheno_cc, sample_include, sample_ct, pheno_cc_collapsed);
  float* pheno_f_iter = *pheno_f_ptr;
  for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    // can use the bitvector equivalent of GenoarrLookup...(), but this isn't
    // in a critical loop so I'll postpone writing those functions for now
    *pheno_f_iter++ = kSmallFloats[IsSet(pheno_cc_collapsed, sample_idx)];
  }
  const uint32_t sample_remv = sample_ctav - sample_ct;
  ZeroFArr(sample_remv, pheno_f_iter);
  PglErr reterr = GlmFillAndTestCovars(sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, max_corr, vif_thresh, covar_dotprod, corr_buf, inverse_corr_buf, covars_cmaj_d, *cur_covar_names_ptr, glm_err_ptr);
  if (unlikely(reterr)) {
    return (reterr == kPglRetNomem);
  }
  double* covar_read_iter = covars_cmaj_d;
  float* covar_write_iter = *covars_cmaj_f_ptr;
  for (uintptr_t covar_idx = 0; covar_idx != new_nonlocal_covar_ct; ++covar_idx) {
    for (uintptr_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      *covar_write_iter++ = S_CAST(float, *covar_read_iter++);
    }
    ZeroFArr(sample_remv, covar_write_iter);
    covar_write_iter = &(covar_write_iter[sample_remv]);
  }
  if (xtx_state) {
    // error-out should be impossible
    InitNmPrecomp(covars_cmaj_d, covar_dotprod, corr_buf, inverse_corr_buf, sample_ct, 0, new_covar_ct, xtx_state, *nm_precomp_ptr);
  }
  BigstackReset(bigstack_mark);
  if (gcount_case_interleaved_vec_ptr) {
    if (unlikely(bigstack_alloc_w(sample_ctv * kWordsPerVec, gcount_case_interleaved_vec_ptr))) {
      return 1;
    }
    ZeroTrailingWords(BitCtToWordCt(sample_ct), pheno_cc_collapsed);
    FillInterleavedMaskVec(pheno_cc_collapsed, sample_ctv, *gcount_case_interleaved_vec_ptr);
  }
  return 0;
}

uint32_t DosageIsConstant(uint64_t dosage_sum, uint64_t dosage_ssq, uint32_t nm_sample_ct) {
  // Dosages are all identical iff dosage_sum * dosage_sum ==
  // dosage_ssq * nm_sample_ct.
  // Unfortunately, uint64 isn't enough for this comparison.
  const uint64_t dosage_sum_hi = dosage_sum >> 32;
  const uint64_t dosage_sum_lo = S_CAST(uint32_t, dosage_sum);
  const uint64_t dosage_ssq_hi = dosage_ssq >> 32;
  const uint64_t dosage_ssq_lo = S_CAST(uint32_t, dosage_ssq);

  // (a * 2^32 + b)^2 = a^2 * 2^64 + 2ab * 2^32 + b^2
  const uint64_t lhs_ab = dosage_sum_hi * dosage_sum_lo;
  const uint64_t lhs_b2 = dosage_sum_lo * dosage_sum_lo;
  const uint64_t lhs_lo = lhs_b2 + (lhs_ab << 33);
  const uint64_t lhs_hi = (lhs_lo < lhs_b2) + (lhs_ab >> 31) + dosage_sum_hi * dosage_sum_hi;

  // (a * 2^32 + b) * c = ac * 2^32 + bc
  const uint64_t rhs_ac = dosage_ssq_hi * nm_sample_ct;
  const uint64_t rhs_bc = dosage_ssq_lo * nm_sample_ct;
  const uint64_t rhs_lo = rhs_bc + (rhs_ac << 32);
  const uint64_t rhs_hi = (rhs_lo < rhs_bc) + (rhs_ac >> 32);
  return (lhs_hi == rhs_hi) && (lhs_lo == rhs_lo);
}

static const float kSmallFloatPairs[32] = PAIR_TABLE16(0.0, 1.0, 2.0, 3.0);

static const float kSmallInvFloatPairs[32] = PAIR_TABLE16(2.0, 1.0, 0.0, 3.0);

static const float kSmallInvFloats[4] = {2.0, 1.0, 0.0, 3.0};

uint32_t GenoarrToFloatsRemoveMissing(const uintptr_t* genoarr, const float* __restrict table, uint32_t sample_ct, float* __restrict dst) {
  assert(sample_ct);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t subgroup_len = kBitsPerWordD2;
  float* dst_iter = dst;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        return S_CAST(uint32_t, dst_iter - dst);
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      const uintptr_t cur_geno = geno_word & 3;
      if (cur_geno < 3) {
        *dst_iter++ = table[cur_geno];
      }
      geno_word >>= 2;
    }
  }
}

// #####
// The following code is based on the winning submission of Pascal Pons in the
// "GWASSpeedup" contest run in April 2013 by Babbage Analytics & Innovation
// and TopCoder, who have donated the results to be used in PLINK.  See:
//   Hill A, Loh PR, Bharadwaj RB, Pons P, Shang J, Guinan E, Lakhani K,
//   Kilty I, Jelinsky SA (2017) Stepwise Distributed Open Innovation Contests
//   for Software Development - Acceleration of Genome-Wide Association
//   Analysis.  Gigascience, 6.
// #####

#ifdef __LP64__
// The two instances of fmath_exp_ps() are C ports of Shigeo Mitsunari's fast
// math library functions posted at https://github.com/herumi/fmath .  License
// is http://opensource.org/licenses/BSD-3-Clause .
// (I tried porting fmath_log_ps, but it turns out that Firth regression needs
// double-precision log accuracy; logf() actually interferes with convergence.)

// programmatically generated by:
// typedef union {
//   float f4;
//   uint32_t u4;
// } __uni4;
//
// __uni4 u4;
// int32_t ii;
// for (ii = 0; ii < 1024; ii++) {
//   u4.f4 = pow(2.0f, ((float)ii) / 1024.0);
//   printf("0x%08x", u4.u4 & 0x7fffff);
//   if (ii % 4 != 3) {
//     printf(", ");
//   } else {
//     printf(",\n");
//   }
// }

const uint32_t kFloatExpLookupInt[]
#ifdef FVEC_32
  __attribute__((aligned(32)))
#else
  __attribute__((aligned(16)))
#endif
  = {
0x00000000, 0x00001630, 0x00002c64, 0x0000429c,
0x000058d8, 0x00006f17, 0x0000855b, 0x00009ba2,
0x0000b1ed, 0x0000c83c, 0x0000de8f, 0x0000f4e6,
0x00010b41, 0x0001219f, 0x00013802, 0x00014e68,
0x000164d2, 0x00017b40, 0x000191b2, 0x0001a828,
0x0001bea1, 0x0001d51f, 0x0001eba1, 0x00020226,
0x000218af, 0x00022f3c, 0x000245ce, 0x00025c63,
0x000272fc, 0x00028998, 0x0002a039, 0x0002b6de,
0x0002cd87, 0x0002e433, 0x0002fae4, 0x00031198,
0x00032850, 0x00033f0d, 0x000355cd, 0x00036c91,
0x00038359, 0x00039a25, 0x0003b0f5, 0x0003c7c9,
0x0003dea1, 0x0003f57d, 0x00040c5d, 0x00042341,
0x00043a29, 0x00045115, 0x00046804, 0x00047ef8,
0x000495f0, 0x0004aceb, 0x0004c3eb, 0x0004daef,
0x0004f1f6, 0x00050902, 0x00052012, 0x00053725,
0x00054e3d, 0x00056558, 0x00057c78, 0x0005939c,
0x0005aac3, 0x0005c1ef, 0x0005d91f, 0x0005f052,
0x0006078a, 0x00061ec6, 0x00063606, 0x00064d4a,
0x00066491, 0x00067bdd, 0x0006932d, 0x0006aa81,
0x0006c1d9, 0x0006d935, 0x0006f095, 0x000707f9,
0x00071f62, 0x000736ce, 0x00074e3e, 0x000765b3,
0x00077d2b, 0x000794a8, 0x0007ac28, 0x0007c3ad,
0x0007db35, 0x0007f2c2, 0x00080a53, 0x000821e8,
0x00083981, 0x0008511e, 0x000868c0, 0x00088065,
0x0008980f, 0x0008afbc, 0x0008c76e, 0x0008df23,
0x0008f6dd, 0x00090e9b, 0x0009265d, 0x00093e24,
0x000955ee, 0x00096dbc, 0x0009858f, 0x00099d66,
0x0009b541, 0x0009cd20, 0x0009e503, 0x0009fcea,
0x000a14d5, 0x000a2cc5, 0x000a44b9, 0x000a5cb1,
0x000a74ad, 0x000a8cad, 0x000aa4b1, 0x000abcba,
0x000ad4c6, 0x000aecd7, 0x000b04ec, 0x000b1d05,
0x000b3523, 0x000b4d44, 0x000b656a, 0x000b7d94,
0x000b95c2, 0x000badf4, 0x000bc62b, 0x000bde65,
0x000bf6a4, 0x000c0ee7, 0x000c272f, 0x000c3f7a,
0x000c57ca, 0x000c701e, 0x000c8876, 0x000ca0d2,
0x000cb933, 0x000cd198, 0x000cea01, 0x000d026e,
0x000d1adf, 0x000d3355, 0x000d4bcf, 0x000d644d,
0x000d7cd0, 0x000d9556, 0x000dade1, 0x000dc671,
0x000ddf04, 0x000df79c, 0x000e1038, 0x000e28d8,
0x000e417d, 0x000e5a25, 0x000e72d3, 0x000e8b84,
0x000ea43a, 0x000ebcf3, 0x000ed5b2, 0x000eee74,
0x000f073b, 0x000f2006, 0x000f38d5, 0x000f51a9,
0x000f6a81, 0x000f835d, 0x000f9c3e, 0x000fb523,
0x000fce0c, 0x000fe6fa, 0x000fffec, 0x001018e2,
0x001031dc, 0x00104adb, 0x001063de, 0x00107ce6,
0x001095f2, 0x0010af02, 0x0010c816, 0x0010e12f,
0x0010fa4d, 0x0011136e, 0x00112c94, 0x001145be,
0x00115eed, 0x00117820, 0x00119158, 0x0011aa93,
0x0011c3d3, 0x0011dd18, 0x0011f661, 0x00120fae,
0x00122900, 0x00124256, 0x00125bb0, 0x0012750f,
0x00128e72, 0x0012a7da, 0x0012c146, 0x0012dab7,
0x0012f42c, 0x00130da5, 0x00132723, 0x001340a5,
0x00135a2b, 0x001373b6, 0x00138d46, 0x0013a6d9,
0x0013c072, 0x0013da0e, 0x0013f3af, 0x00140d55,
0x001426ff, 0x001440ae, 0x00145a60, 0x00147418,
0x00148dd4, 0x0014a794, 0x0014c159, 0x0014db22,
0x0014f4f0, 0x00150ec2, 0x00152898, 0x00154274,
0x00155c53, 0x00157637, 0x00159020, 0x0015aa0d,
0x0015c3ff, 0x0015ddf5, 0x0015f7ef, 0x001611ee,
0x00162bf2, 0x001645fa, 0x00166006, 0x00167a18,
0x0016942d, 0x0016ae47, 0x0016c866, 0x0016e289,
0x0016fcb1, 0x001716dd, 0x0017310e, 0x00174b43,
0x0017657d, 0x00177fbc, 0x001799ff, 0x0017b446,
0x0017ce92, 0x0017e8e3, 0x00180338, 0x00181d92,
0x001837f0, 0x00185253, 0x00186cbb, 0x00188727,
0x0018a197, 0x0018bc0d, 0x0018d686, 0x0018f105,
0x00190b88, 0x0019260f, 0x0019409c, 0x00195b2c,
0x001975c2, 0x0019905c, 0x0019aafa, 0x0019c59e,
0x0019e046, 0x0019faf2, 0x001a15a3, 0x001a3059,
0x001a4b13, 0x001a65d2, 0x001a8096, 0x001a9b5e,
0x001ab62b, 0x001ad0fd, 0x001aebd3, 0x001b06ae,
0x001b218d, 0x001b3c71, 0x001b575a, 0x001b7248,
0x001b8d3a, 0x001ba831, 0x001bc32c, 0x001bde2c,
0x001bf931, 0x001c143b, 0x001c2f49, 0x001c4a5c,
0x001c6573, 0x001c8090, 0x001c9bb1, 0x001cb6d6,
0x001cd201, 0x001ced30, 0x001d0864, 0x001d239c,
0x001d3eda, 0x001d5a1c, 0x001d7562, 0x001d90ae,
0x001dabfe, 0x001dc753, 0x001de2ad, 0x001dfe0b,
0x001e196e, 0x001e34d6, 0x001e5043, 0x001e6bb4,
0x001e872a, 0x001ea2a5, 0x001ebe25, 0x001ed9a9,
0x001ef532, 0x001f10c0, 0x001f2c53, 0x001f47eb,
0x001f6387, 0x001f7f28, 0x001f9ace, 0x001fb679,
0x001fd228, 0x001feddc, 0x00200996, 0x00202553,
0x00204116, 0x00205cde, 0x002078aa, 0x0020947b,
0x0020b051, 0x0020cc2c, 0x0020e80b, 0x002103f0,
0x00211fd9, 0x00213bc7, 0x002157ba, 0x002173b2,
0x00218faf, 0x0021abb0, 0x0021c7b7, 0x0021e3c2,
0x0021ffd2, 0x00221be7, 0x00223801, 0x0022541f,
0x00227043, 0x00228c6b, 0x0022a899, 0x0022c4cb,
0x0022e102, 0x0022fd3e, 0x0023197f, 0x002335c5,
0x0023520f, 0x00236e5f, 0x00238ab3, 0x0023a70d,
0x0023c36b, 0x0023dfce, 0x0023fc37, 0x002418a4,
0x00243516, 0x0024518d, 0x00246e08, 0x00248a89,
0x0024a70f, 0x0024c39a, 0x0024e029, 0x0024fcbe,
0x00251958, 0x002535f6, 0x00255299, 0x00256f42,
0x00258bef, 0x0025a8a2, 0x0025c559, 0x0025e215,
0x0025fed7, 0x00261b9d, 0x00263868, 0x00265538,
0x0026720e, 0x00268ee8, 0x0026abc7, 0x0026c8ac,
0x0026e595, 0x00270283, 0x00271f76, 0x00273c6f,
0x0027596c, 0x0027766e, 0x00279376, 0x0027b082,
0x0027cd94, 0x0027eaaa, 0x002807c6, 0x002824e6,
0x0028420c, 0x00285f37, 0x00287c66, 0x0028999b,
0x0028b6d5, 0x0028d414, 0x0028f158, 0x00290ea1,
0x00292bef, 0x00294942, 0x0029669b, 0x002983f8,
0x0029a15b, 0x0029bec2, 0x0029dc2f, 0x0029f9a1,
0x002a1718, 0x002a3494, 0x002a5215, 0x002a6f9b,
0x002a8d26, 0x002aaab7, 0x002ac84c, 0x002ae5e7,
0x002b0387, 0x002b212c, 0x002b3ed6, 0x002b5c85,
0x002b7a3a, 0x002b97f3, 0x002bb5b2, 0x002bd376,
0x002bf13f, 0x002c0f0d, 0x002c2ce0, 0x002c4ab9,
0x002c6897, 0x002c867a, 0x002ca462, 0x002cc24f,
0x002ce041, 0x002cfe39, 0x002d1c36, 0x002d3a38,
0x002d583f, 0x002d764b, 0x002d945d, 0x002db274,
0x002dd090, 0x002deeb1, 0x002e0cd8, 0x002e2b03,
0x002e4934, 0x002e676b, 0x002e85a6, 0x002ea3e7,
0x002ec22d, 0x002ee078, 0x002efec8, 0x002f1d1e,
0x002f3b79, 0x002f59d9, 0x002f783e, 0x002f96a9,
0x002fb519, 0x002fd38e, 0x002ff209, 0x00301089,
0x00302f0e, 0x00304d98, 0x00306c28, 0x00308abd,
0x0030a957, 0x0030c7f7, 0x0030e69c, 0x00310546,
0x003123f6, 0x003142aa, 0x00316165, 0x00318024,
0x00319ee9, 0x0031bdb3, 0x0031dc83, 0x0031fb57,
0x00321a32, 0x00323911, 0x003257f6, 0x003276e0,
0x003295d0, 0x0032b4c5, 0x0032d3bf, 0x0032f2bf,
0x003311c4, 0x003330cf, 0x00334fde, 0x00336ef4,
0x00338e0e, 0x0033ad2e, 0x0033cc54, 0x0033eb7e,
0x00340aaf, 0x003429e4, 0x0034491f, 0x00346860,
0x003487a6, 0x0034a6f1, 0x0034c642, 0x0034e598,
0x003504f3, 0x00352454, 0x003543bb, 0x00356327,
0x00358298, 0x0035a20f, 0x0035c18b, 0x0035e10d,
0x00360094, 0x00362020, 0x00363fb2, 0x00365f4a,
0x00367ee7, 0x00369e89, 0x0036be31, 0x0036dddf,
0x0036fd92, 0x00371d4a, 0x00373d08, 0x00375ccc,
0x00377c95, 0x00379c63, 0x0037bc37, 0x0037dc11,
0x0037fbf0, 0x00381bd4, 0x00383bbe, 0x00385bae,
0x00387ba3, 0x00389b9e, 0x0038bb9e, 0x0038dba4,
0x0038fbaf, 0x00391bc0, 0x00393bd7, 0x00395bf3,
0x00397c14, 0x00399c3b, 0x0039bc68, 0x0039dc9a,
0x0039fcd2, 0x003a1d10, 0x003a3d53, 0x003a5d9b,
0x003a7dea, 0x003a9e3e, 0x003abe97, 0x003adef6,
0x003aff5b, 0x003b1fc5, 0x003b4035, 0x003b60aa,
0x003b8126, 0x003ba1a6, 0x003bc22d, 0x003be2b9,
0x003c034a, 0x003c23e2, 0x003c447f, 0x003c6521,
0x003c85ca, 0x003ca678, 0x003cc72b, 0x003ce7e5,
0x003d08a4, 0x003d2968, 0x003d4a33, 0x003d6b03,
0x003d8bd8, 0x003dacb4, 0x003dcd95, 0x003dee7c,
0x003e0f68, 0x003e305a, 0x003e5152, 0x003e7250,
0x003e9353, 0x003eb45c, 0x003ed56b, 0x003ef67f,
0x003f179a, 0x003f38ba, 0x003f59df, 0x003f7b0b,
0x003f9c3c, 0x003fbd73, 0x003fdeb0, 0x003ffff2,
0x0040213b, 0x00404289, 0x004063dc, 0x00408536,
0x0040a695, 0x0040c7fb, 0x0040e966, 0x00410ad6,
0x00412c4d, 0x00414dc9, 0x00416f4b, 0x004190d3,
0x0041b261, 0x0041d3f5, 0x0041f58e, 0x0042172d,
0x004238d2, 0x00425a7d, 0x00427c2e, 0x00429de4,
0x0042bfa1, 0x0042e163, 0x0043032b, 0x004324f9,
0x004346cd, 0x004368a7, 0x00438a86, 0x0043ac6b,
0x0043ce57, 0x0043f048, 0x0044123f, 0x0044343c,
0x0044563f, 0x00447848, 0x00449a56, 0x0044bc6b,
0x0044de85, 0x004500a5, 0x004522cc, 0x004544f8,
0x0045672a, 0x00458962, 0x0045aba0, 0x0045cde4,
0x0045f02e, 0x0046127e, 0x004634d3, 0x0046572f,
0x00467991, 0x00469bf8, 0x0046be66, 0x0046e0d9,
0x00470353, 0x004725d2, 0x00474858, 0x00476ae3,
0x00478d75, 0x0047b00c, 0x0047d2aa, 0x0047f54d,
0x004817f7, 0x00483aa6, 0x00485d5b, 0x00488017,
0x0048a2d8, 0x0048c5a0, 0x0048e86d, 0x00490b41,
0x00492e1b, 0x004950fa, 0x004973e0, 0x004996cc,
0x0049b9be, 0x0049dcb5, 0x0049ffb3, 0x004a22b7,
0x004a45c1, 0x004a68d1, 0x004a8be8, 0x004aaf04,
0x004ad226, 0x004af54f, 0x004b187d, 0x004b3bb2,
0x004b5eed, 0x004b822e, 0x004ba575, 0x004bc8c2,
0x004bec15, 0x004c0f6e, 0x004c32ce, 0x004c5633,
0x004c799f, 0x004c9d11, 0x004cc089, 0x004ce407,
0x004d078c, 0x004d2b16, 0x004d4ea7, 0x004d723d,
0x004d95da, 0x004db97e, 0x004ddd27, 0x004e00d6,
0x004e248c, 0x004e4848, 0x004e6c0a, 0x004e8fd2,
0x004eb3a1, 0x004ed775, 0x004efb50, 0x004f1f31,
0x004f4319, 0x004f6706, 0x004f8afa, 0x004faef4,
0x004fd2f4, 0x004ff6fb, 0x00501b08, 0x00503f1b,
0x00506334, 0x00508753, 0x0050ab79, 0x0050cfa5,
0x0050f3d7, 0x00511810, 0x00513c4f, 0x00516094,
0x005184df, 0x0051a931, 0x0051cd89, 0x0051f1e7,
0x0052164c, 0x00523ab7, 0x00525f28, 0x005283a0,
0x0052a81e, 0x0052cca2, 0x0052f12c, 0x005315bd,
0x00533a54, 0x00535ef2, 0x00538396, 0x0053a840,
0x0053ccf1, 0x0053f1a8, 0x00541665, 0x00543b29,
0x00545ff3, 0x005484c3, 0x0054a99a, 0x0054ce77,
0x0054f35b, 0x00551845, 0x00553d35, 0x0055622c,
0x00558729, 0x0055ac2d, 0x0055d137, 0x0055f647,
0x00561b5e, 0x0056407b, 0x0056659f, 0x00568ac9,
0x0056affa, 0x0056d531, 0x0056fa6e, 0x00571fb2,
0x005744fd, 0x00576a4e, 0x00578fa5, 0x0057b503,
0x0057da67, 0x0057ffd2, 0x00582543, 0x00584abb,
0x00587039, 0x005895be, 0x0058bb49, 0x0058e0db,
0x00590673, 0x00592c12, 0x005951b8, 0x00597763,
0x00599d16, 0x0059c2cf, 0x0059e88e, 0x005a0e54,
0x005a3421, 0x005a59f4, 0x005a7fcd, 0x005aa5ae,
0x005acb94, 0x005af182, 0x005b1776, 0x005b3d70,
0x005b6371, 0x005b8979, 0x005baf87, 0x005bd59c,
0x005bfbb8, 0x005c21da, 0x005c4802, 0x005c6e32,
0x005c9468, 0x005cbaa4, 0x005ce0e7, 0x005d0731,
0x005d2d82, 0x005d53d9, 0x005d7a36, 0x005da09b,
0x005dc706, 0x005ded77, 0x005e13f0, 0x005e3a6f,
0x005e60f5, 0x005e8781, 0x005eae14, 0x005ed4ae,
0x005efb4e, 0x005f21f5, 0x005f48a3, 0x005f6f58,
0x005f9613, 0x005fbcd5, 0x005fe39e, 0x00600a6d,
0x00603143, 0x00605820, 0x00607f03, 0x0060a5ee,
0x0060ccdf, 0x0060f3d7, 0x00611ad5, 0x006141db,
0x006168e7, 0x00618ffa, 0x0061b713, 0x0061de34,
0x0062055b, 0x00622c89, 0x006253be, 0x00627af9,
0x0062a23c, 0x0062c985, 0x0062f0d5, 0x0063182c,
0x00633f89, 0x006366ee, 0x00638e59, 0x0063b5cb,
0x0063dd44, 0x006404c4, 0x00642c4b, 0x006453d8,
0x00647b6d, 0x0064a308, 0x0064caaa, 0x0064f253,
0x00651a03, 0x006541b9, 0x00656977, 0x0065913c,
0x0065b907, 0x0065e0d9, 0x006608b2, 0x00663092,
0x00665879, 0x00668067, 0x0066a85c, 0x0066d058,
0x0066f85b, 0x00672064, 0x00674875, 0x0067708c,
0x006798ab, 0x0067c0d0, 0x0067e8fd, 0x00681130,
0x0068396a, 0x006861ac, 0x006889f4, 0x0068b243,
0x0068da99, 0x006902f7, 0x00692b5b, 0x006953c6,
0x00697c38, 0x0069a4b1, 0x0069cd32, 0x0069f5b9,
0x006a1e47, 0x006a46dd, 0x006a6f79, 0x006a981c,
0x006ac0c7, 0x006ae978, 0x006b1231, 0x006b3af1,
0x006b63b7, 0x006b8c85, 0x006bb55a, 0x006bde36,
0x006c0719, 0x006c3003, 0x006c58f4, 0x006c81ec,
0x006caaec, 0x006cd3f2, 0x006cfd00, 0x006d2614,
0x006d4f30, 0x006d7853, 0x006da17d, 0x006dcaae,
0x006df3e7, 0x006e1d26, 0x006e466d, 0x006e6fbb,
0x006e9910, 0x006ec26c, 0x006eebcf, 0x006f1539,
0x006f3eab, 0x006f6824, 0x006f91a4, 0x006fbb2b,
0x006fe4ba, 0x00700e4f, 0x007037ec, 0x00706190,
0x00708b3b, 0x0070b4ee, 0x0070dea8, 0x00710868,
0x00713231, 0x00715c00, 0x007185d7, 0x0071afb5,
0x0071d99a, 0x00720386, 0x00722d7a, 0x00725775,
0x00728177, 0x0072ab81, 0x0072d592, 0x0072ffaa,
0x007329c9, 0x007353f0, 0x00737e1e, 0x0073a853,
0x0073d290, 0x0073fcd4, 0x0074271f, 0x00745172,
0x00747bcc, 0x0074a62d, 0x0074d096, 0x0074fb06,
0x0075257d, 0x00754ffc, 0x00757a82, 0x0075a50f,
0x0075cfa4, 0x0075fa40, 0x007624e4, 0x00764f8f,
0x00767a41, 0x0076a4fb, 0x0076cfbc, 0x0076fa85,
0x00772555, 0x0077502d, 0x00777b0b, 0x0077a5f2,
0x0077d0df, 0x0077fbd5, 0x007826d1, 0x007851d5,
0x00787ce1, 0x0078a7f4, 0x0078d30e, 0x0078fe30,
0x0079295a, 0x0079548b, 0x00797fc3, 0x0079ab03,
0x0079d64a, 0x007a0199, 0x007a2cf0, 0x007a584d,
0x007a83b3, 0x007aaf20, 0x007ada94, 0x007b0610,
0x007b3194, 0x007b5d1f, 0x007b88b2, 0x007bb44c,
0x007bdfed, 0x007c0b97, 0x007c3748, 0x007c6300,
0x007c8ec0, 0x007cba88, 0x007ce657, 0x007d122e,
0x007d3e0c, 0x007d69f2, 0x007d95e0, 0x007dc1d5,
0x007dedd2, 0x007e19d6, 0x007e45e2, 0x007e71f6,
0x007e9e11, 0x007eca34, 0x007ef65f, 0x007f2291,
0x007f4ecb, 0x007f7b0d, 0x007fa756, 0x007fd3a7
};

#  ifdef FVEC_32
static inline VecF fmath_exp_ps(VecF xxv) {
  __m256 xx = R_CAST(__m256, xxv);
  const __m256i mask7ff = {0x7fffffff7fffffffLLU, 0x7fffffff7fffffffLLU, 0x7fffffff7fffffffLLU, 0x7fffffff7fffffffLLU};
  // 88
  const __m256i max_x = {0x42b0000042b00000LLU, 0x42b0000042b00000LLU, 0x42b0000042b00000LLU, 0x42b0000042b00000LLU};
  // -88
  // more sensible 0xc2b00000... not used here due to "narrowing conversion"
  // warning
  const __m256i min_x = {-0x3d4fffff3d500000LL, -0x3d4fffff3d500000LL, -0x3d4fffff3d500000LL, -0x3d4fffff3d500000LL};
  // 2^10 / log(2)
  const __m256i const_aa = {0x44b8aa3b44b8aa3bLLU, 0x44b8aa3b44b8aa3bLLU, 0x44b8aa3b44b8aa3bLLU, 0x44b8aa3b44b8aa3bLLU};
  // log(2) / 2^10
  const __m256i const_bb = {0x3a3172183a317218LLU, 0x3a3172183a317218LLU, 0x3a3172183a317218LLU, 0x3a3172183a317218LLU};
  const __m256i f1 = {0x3f8000003f800000LLU, 0x3f8000003f800000LLU, 0x3f8000003f800000LLU, 0x3f8000003f800000LLU};
  const __m256i mask_s = {0x3ff000003ffLLU, 0x3ff000003ffLLU, 0x3ff000003ffLLU, 0x3ff000003ffLLU};
  const __m256i i127s = {0x1fc000001fc00LLU, 0x1fc000001fc00LLU, 0x1fc000001fc00LLU, 0x1fc000001fc00LLU};
  const __m256i limit = _mm256_castps_si256(_mm256_and_ps(xx, R_CAST(__m256, mask7ff)));
  const int32_t over = _mm256_movemask_epi8(_mm256_cmpgt_epi32(limit, max_x));
  if (over) {
    xx = _mm256_min_ps(xx, R_CAST(__m256, max_x));
    xx = _mm256_max_ps(xx, R_CAST(__m256, min_x));
  }
  const __m256i rr = _mm256_cvtps_epi32(_mm256_mul_ps(xx, R_CAST(__m256, const_aa)));
  __m256 tt = _mm256_fnmadd_ps(_mm256_cvtepi32_ps(rr), R_CAST(__m256, const_bb), xx);
  tt = _mm256_add_ps(R_CAST(__m256, f1), tt);
  const __m256i v8 = _mm256_and_si256(rr, mask_s);
  __m256i u8 = _mm256_add_epi32(rr, i127s);
  u8 = _mm256_srli_epi32(u8, 10);
  u8 = _mm256_slli_epi32(u8, 23);
  __m256i ti = _mm256_i32gather_epi32(R_CAST(const int*, kFloatExpLookupInt), v8, 4);
  __m256 t0 = _mm256_castsi256_ps(ti);
  t0 = _mm256_or_ps(t0, _mm256_castsi256_ps(u8));
  return R_CAST(VecF, _mm256_mul_ps(tt, t0));
}

// For equivalent "normal" C/C++ code, see the non-__LP64__ versions of these
// functions.

// N.B. This requires all mm[] rows to be zero-padded at the end, and there
// can't be nan values at the end of vect[].  (The other way around works too.)
//
// This is currently a bit faster than sgemm and sgemv on my Mac, so it isn't
// appropriate to throw out this code yet.
static inline void MultMatrixDxnVectN(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uintptr_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  uint32_t row_idx = 0;
  __m256 s1;
  __m256 s2;
  __m256 s3;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    // Handle 4 rows at a time in this loop, regardless of vector size.
    for (; row_idx < row_ctm3; row_idx += 4) {
      s1 = _mm256_setzero_ps();
      s2 = _mm256_setzero_ps();
      s3 = _mm256_setzero_ps();
      __m256 s4 = _mm256_setzero_ps();
      for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
        const float* mm_ptr = &(mm[row_idx * col_ctav + col_idx]);
        const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
        __m256 a1 = _mm256_load_ps(mm_ptr);
        __m256 a2 = _mm256_load_ps(&(mm_ptr[col_ctav]));
        __m256 a3 = _mm256_load_ps(&(mm_ptr[2 * col_ctav]));
        __m256 a4 = _mm256_load_ps(&(mm_ptr[3 * col_ctav]));
        s1 = _mm256_fmadd_ps(a1, vv, s1);
        s2 = _mm256_fmadd_ps(a2, vv, s2);
        s3 = _mm256_fmadd_ps(a3, vv, s3);
        s4 = _mm256_fmadd_ps(a4, vv, s4);
      }
      *dest++ = VecFHsum(R_CAST(VecF, s1));
      *dest++ = VecFHsum(R_CAST(VecF, s2));
      *dest++ = VecFHsum(R_CAST(VecF, s3));
      *dest++ = VecFHsum(R_CAST(VecF, s4));
    }
  }
  s1 = _mm256_setzero_ps();
  s2 = _mm256_setzero_ps();
  s3 = _mm256_setzero_ps();
  switch (row_ct % 4) {
  case 3:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
      const float* mm_ptr = &(mm[row_idx * col_ctav + col_idx]);
      const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
      __m256 a1 = _mm256_load_ps(mm_ptr);
      __m256 a2 = _mm256_load_ps(&(mm_ptr[col_ctav]));
      __m256 a3 = _mm256_load_ps(&(mm_ptr[2 * col_ctav]));
      s1 = _mm256_fmadd_ps(a1, vv, s1);
      s2 = _mm256_fmadd_ps(a2, vv, s2);
      s3 = _mm256_fmadd_ps(a3, vv, s3);
    }
    *dest++ = VecFHsum(R_CAST(VecF, s1));
    *dest++ = VecFHsum(R_CAST(VecF, s2));
    *dest = VecFHsum(R_CAST(VecF, s3));
    break;
  case 2:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
      const float* mm_ptr = &(mm[row_idx * col_ctav + col_idx]);
      const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
      __m256 a1 = _mm256_load_ps(mm_ptr);
      __m256 a2 = _mm256_load_ps(&(mm_ptr[col_ctav]));
      s1 = _mm256_fmadd_ps(a1, vv, s1);
      s2 = _mm256_fmadd_ps(a2, vv, s2);
    }
    *dest++ = VecFHsum(R_CAST(VecF, s1));
    *dest = VecFHsum(R_CAST(VecF, s2));
    break;
  case 1:
    for (uint32_t col_idx = 0; col_idx < col_ct; col_idx += kFloatPerFVec) {
      const __m256 vv = _mm256_load_ps(&(vect[col_idx]));
      __m256 a1 = _mm256_load_ps(&(mm[row_idx * col_ctav + col_idx]));
      s1 = _mm256_fmadd_ps(a1, vv, s1);
    }
    *dest = VecFHsum(R_CAST(VecF, s1));
    break;
  }
}

#  else  // !FVEC_32
const float* const kFloatExpLookup = R_CAST(const float*, kFloatExpLookupInt);

static inline VecF fmath_exp_ps(VecF xxv) {
  __m128 xx = xxv;
  const __m128i mask7ff = {0x7fffffff7fffffffLLU, 0x7fffffff7fffffffLLU};

  // 88
  const __m128i max_x = {0x42b0000042b00000LLU, 0x42b0000042b00000LLU};
  // -88
  // more sensible 0xc2b00000... not used here due to "narrowing conversion"
  // warning
  const __m128i min_x = {-0x3d4fffff3d500000LL, -0x3d4fffff3d500000LL};
  // 2^10 / log(2)
  const __m128i const_aa = {0x44b8aa3b44b8aa3bLLU, 0x44b8aa3b44b8aa3bLLU};
  // log(2) / 2^10
  const __m128i const_bb = {0x3a3172183a317218LLU, 0x3a3172183a317218LLU};

  const __m128i f1 = {0x3f8000003f800000LLU, 0x3f8000003f800000LLU};
  const __m128i mask_s = {0x3ff000003ffLLU, 0x3ff000003ffLLU};
  const __m128i i127s = {0x1fc000001fc00LLU, 0x1fc000001fc00LLU};
  const __m128i limit = _mm_castps_si128(_mm_and_ps(xx, R_CAST(__m128, mask7ff)));
  const int32_t over = _mm_movemask_epi8(_mm_cmpgt_epi32(limit, max_x));
  if (over) {
    xx = _mm_min_ps(xx, R_CAST(__m128, max_x));
    xx = _mm_max_ps(xx, R_CAST(__m128, min_x));
  }
  const __m128i rr = _mm_cvtps_epi32(_mm_mul_ps(xx, R_CAST(__m128, const_aa)));
  __m128 tt = _mm_sub_ps(xx, _mm_mul_ps(_mm_cvtepi32_ps(rr), R_CAST(__m128, const_bb)));
  tt = _mm_add_ps(tt, R_CAST(__m128, f1));
  const __m128i v4 = _mm_and_si128(rr, mask_s);
  __m128i u4 = _mm_add_epi32(rr, i127s);
  u4 = _mm_srli_epi32(u4, 10);
  u4 = _mm_slli_epi32(u4, 23);
  const uint32_t v0 = _mm_cvtsi128_si32(v4);
  // uint32_t v1 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(2)));
  // uint32_t v2 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(4)));
  // uint32_t v3 = ((int32_t)(uint16_t)__builtin_ia32_vec_ext_v8hi((__v8hi)(__m128i)(v4), (int32_t)(6)));
  // make this work with LLVM
  const uint32_t v1 = _mm_extract_epi16(R_CAST(__m128i, v4), 2);
  const uint32_t v2 = _mm_extract_epi16(R_CAST(__m128i, v4), 4);
  const uint32_t v3 = _mm_extract_epi16(R_CAST(__m128i, v4), 6);

  __m128 t0 = _mm_set_ss(kFloatExpLookup[v0]);
  __m128 t1 = _mm_set_ss(kFloatExpLookup[v1]);
  const __m128 t2 = _mm_set_ss(kFloatExpLookup[v2]);
  const __m128 t3 = _mm_set_ss(kFloatExpLookup[v3]);
  t1 = _mm_movelh_ps(t1, t3);
  t1 = _mm_castsi128_ps(_mm_slli_epi64(_mm_castps_si128(t1), 32));
  t0 = _mm_movelh_ps(t0, t2);
  t0 = _mm_or_ps(t0, t1);
  t0 = _mm_or_ps(t0, _mm_castsi128_ps(u4));
  tt = _mm_mul_ps(tt, t0);
  return R_CAST(VecF, tt);
}

static inline void MultMatrixDxnVectN(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uint32_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  ColMajorFvectorMatrixMultiplyStrided(vect, mm, col_ct, col_ctav, row_ct, dest);
}

#  endif  // !FVEC_32

static inline void LogisticSse(uint32_t nn, float* vect) {
  const VecF zero = vecf_setzero();
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF aa = *R_CAST(VecF*, &(vect[uii]));
    aa = zero - aa;
    // tried substituting in vexpf() here on OS X; it was slower without being
    // significantly more accurate.
    aa = fmath_exp_ps(aa);
    aa = aa + one;
    aa = one / aa;
    *R_CAST(VecF*, &(vect[uii])) = aa;
  }
}

static inline void ComputeVAndPMinusY(const float* yy, uint32_t nn, float* pp, float* vv) {
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF ptmp = *R_CAST(VecF*, &(pp[uii]));
    VecF one_minus_ptmp = one - ptmp;
    *R_CAST(VecF*, &(vv[uii])) = ptmp * one_minus_ptmp;
    VecF ytmp = *R_CAST(const VecF*, &(yy[uii]));
    *R_CAST(VecF*, &(pp[uii])) = ptmp - ytmp;
  }
}

static inline void ComputeV(const float* pp, uint32_t nn, float* vv) {
  const VecF one = VCONST_F(1.0);
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF ptmp = *R_CAST(const VecF*, &(pp[uii]));
    VecF one_minus_ptmp = one - ptmp;
    *R_CAST(VecF*, &(vv[uii])) = ptmp * one_minus_ptmp;
  }
}

static inline float TripleProduct(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  VecF sum = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    VecF aa = *R_CAST(const VecF*, &(v1[uii]));
    VecF bb = *R_CAST(const VecF*, &(v2[uii]));
    VecF cc = *R_CAST(const VecF*, &(v3[uii]));
    sum = sum + aa * bb * cc;
  }
  return VecFHsum(sum);
}

static inline void ComputeTwoDiagTripleProduct(const float* aa, const float* bb, const float* vv, uint32_t nn, float* __restrict raa_ptr, float* __restrict rab_ptr, float* __restrict rbb_ptr) {
  VecF saa = vecf_setzero();
  VecF sab = vecf_setzero();
  VecF sbb = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    const VecF vtmp = *R_CAST(const VecF*, &(vv[uii]));
    const VecF atmp = *R_CAST(const VecF*, &(aa[uii]));
    const VecF btmp = *R_CAST(const VecF*, &(bb[uii]));
    const VecF av = atmp * vtmp;
    const VecF bv = btmp * vtmp;
    saa = saa + atmp * av;
    sab = sab + atmp * bv;
    sbb = sbb + btmp * bv;
  }
  *raa_ptr = VecFHsum(saa);
  *rab_ptr = VecFHsum(sab);
  *rbb_ptr = VecFHsum(sbb);
}

static inline void ComputeThreeTripleProduct(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  VecF s1 = vecf_setzero();
  VecF s2 = vecf_setzero();
  VecF s3 = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    const VecF a1tmp = *R_CAST(const VecF*, &(a1[uii]));
    const VecF a2tmp = *R_CAST(const VecF*, &(a2[uii]));
    const VecF a3tmp = *R_CAST(const VecF*, &(a3[uii]));
    const VecF vtmp = *R_CAST(const VecF*, &(vv[uii]));
    VecF btmp = *R_CAST(const VecF*, &(bb[uii]));
    btmp = btmp * vtmp;
    s1 = s1 + a1tmp * btmp;
    s2 = s2 + a2tmp * btmp;
    s3 = s3 + a3tmp * btmp;
  }
  *r1_ptr = VecFHsum(s1);
  *r2_ptr = VecFHsum(s2);
  *r3_ptr = VecFHsum(s3);
}

static inline void ComputeTwoPlusOneTripleProduct(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* __restrict r1_ptr, float* __restrict r2_ptr, float* __restrict r3_ptr) {
  VecF s1 = vecf_setzero();
  VecF s2 = vecf_setzero();
  VecF s3 = vecf_setzero();
  for (uint32_t uii = 0; uii < nn; uii += kFloatPerFVec) {
    const VecF a1tmp = *R_CAST(const VecF*, &(a1[uii]));
    const VecF a2tmp = *R_CAST(const VecF*, &(a2[uii]));
    const VecF btmp = *R_CAST(const VecF*, &(bb[uii]));
    const VecF vtmp = *R_CAST(const VecF*, &(vv[uii]));
    const VecF bv = btmp * vtmp;
    s1 = s1 + btmp * bv;
    s2 = s2 + a1tmp * bv;
    s3 = s3 + a2tmp * bv;
  }
  *r1_ptr = VecFHsum(s1);
  *r2_ptr = VecFHsum(s2);
  *r3_ptr = VecFHsum(s3);
}
#else  // no __LP64__ (and hence, unsafe to assume presence of SSE2)
static inline void LogisticSse(uint32_t nn, float* vect) {
  // We use explicit static_cast<float> instead of e.g. 1.0f because
  // handling of the latter is actually implementation-specific; see
  //   http://nullprogram.com/blog/2018/05/01/
  // In particular, that blog post claims that
  //   int float_compare() {
  //     float x = 1.3f;
  //     return x == 1.3f;
  //   }
  // returns 0 under gcc and 1 under clang (with -std=c99 -m32, which is one of
  // plink2's compilation settings)????!!!!!!!
  // Unless the author is outright mistaken, this suggests that use of the f
  // suffix should be considered a bug ~100% of the time.
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vect[uii] = S_CAST(float, 1.0) / (1 + expf(-vect[uii]));
  }
}

static inline void ComputeVAndPMinusY(const float* yy, uint32_t nn, float* pp, float* vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vv[uii] = pp[uii] * (S_CAST(float, 1.0) - pp[uii]);
    pp[uii] -= yy[uii];
  }
}

static inline void ComputeV(const float* pp, uint32_t nn, float* vv) {
  for (uint32_t uii = 0; uii != nn; ++uii) {
    vv[uii] = pp[uii] * (S_CAST(float, 1.0) - pp[uii]);
  }
}

static inline void MultMatrixDxnVectN(const float* mm, const float* vect, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uint32_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  ColMajorFvectorMatrixMultiplyStrided(vect, mm, col_ct, col_ctav, row_ct, dest);
}

static inline float TripleProduct(const float* v1, const float* v2, const float* v3, uint32_t nn) {
  float fxx = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    fxx += (*v1++) * (*v2++) * (*v3++);
  }
  return fxx;
}

static inline void ComputeTwoDiagTripleProduct(const float* aa, const float* bb, const float* vv, uint32_t nn, float* raa_ptr, float* rab_ptr, float* rbb_ptr) {
  float raa = 0.0;
  float rab = 0.0;
  float rbb = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = (*aa++);
    const float fyy = (*bb++);
    float fzz = (*vv++);
    raa += fxx * fxx * fzz;
    fzz *= fyy;
    rab += fxx * fzz;
    rbb += fyy * fzz;
  }
  *raa_ptr = raa;
  *rab_ptr = rab;
  *rbb_ptr = rbb;
}

static inline void ComputeThreeTripleProduct(const float* bb, const float* a1, const float* a2, const float* a3, const float* vv, uint32_t nn, float* r1_ptr, float* r2_ptr, float* r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = (*bb++) * (*vv++);
    r1 += (*a1++) * fxx;
    r2 += (*a2++) * fxx;
    r3 += (*a3++) * fxx;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}

static inline void ComputeTwoPlusOneTripleProduct(const float* bb, const float* a1, const float* a2, const float* vv, uint32_t nn, float* r1_ptr, float* r2_ptr, float* r3_ptr) {
  float r1 = 0.0;
  float r2 = 0.0;
  float r3 = 0.0;
  for (uint32_t uii = 0; uii != nn; ++uii) {
    const float fxx = (*bb++);
    const float fyy = fxx * (*vv++);
    r1 += fxx * fyy;
    r2 += (*a1++) * fyy;
    r3 += (*a2++) * fyy;
  }
  *r1_ptr = r1;
  *r2_ptr = r2;
  *r3_ptr = r3;
}
#endif
double ComputeLoglik(const float* yy, const float* pp, uint32_t sample_ct) {
  // possible todo: look for a high-precision way to accelerate this.
  double loglik = 0.0;
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const double new_pi = S_CAST(double, pp[sample_idx]);
    loglik += (yy[sample_idx] != S_CAST(float, 0.0))? log(new_pi) : log1p(-new_pi);
  }
  return loglik;
}

// M V M^T
// This is the biggest logistic/Firth regression bottleneck.
// Tried to replace this with sqrt(v) followed by ssyrk, but that was slower.
// Also tried to take advantage of the first row of M being constant-1, and
// managed to fail.
void ComputeHessian(const float* mm, const float* vv, uint32_t col_ct, uint32_t row_ct, float* dest) {
  const uintptr_t col_ctav = RoundUpPow2(col_ct, kFloatPerFVec);
  const uintptr_t row_ctav = RoundUpPow2(row_ct, kFloatPerFVec);
  const uintptr_t row_ctavp1 = row_ctav + 1;
  if (row_ct > 3) {
    const uint32_t row_ctm3 = row_ct - 3;
    for (uint32_t row_idx = 0; row_idx < row_ctm3; row_idx += 3) {
      const float* mm_cur = &(mm[row_idx * col_ctav]);
      ComputeTwoDiagTripleProduct(mm_cur, &(mm_cur[col_ctav]), vv, col_ct, &(dest[row_idx * row_ctavp1]), &(dest[(row_idx + 1) * row_ctavp1 - 1]), &(dest[(row_idx + 1) * row_ctavp1]));
      ComputeTwoPlusOneTripleProduct(&(mm_cur[2 * col_ctav]), &(mm_cur[col_ctav]), mm_cur, vv, col_ct, &(dest[(row_idx + 2) * row_ctavp1]), &(dest[(row_idx + 2) * row_ctavp1 - 1]), &(dest[(row_idx + 2) * row_ctavp1 - 2]));
      for (uint32_t row_idx2 = row_idx + 3; row_idx2 != row_ct; ++row_idx2) {
        ComputeThreeTripleProduct(&(mm[row_idx2 * col_ctav]), mm_cur, &(mm_cur[col_ctav]), &(mm_cur[2 * col_ctav]), vv, col_ct, &(dest[row_idx2 * row_ctav + row_idx]), &(dest[row_idx2 * row_ctav + row_idx + 1]), &(dest[row_idx2 * row_ctav + row_idx + 2]));
      }
    }
  }
  switch (row_ct % 3) {
  case 0:
    ComputeTwoPlusOneTripleProduct(&(mm[(row_ct - 3) * col_ctav]), &(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 3) * row_ctavp1]), &(dest[(row_ct - 2) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1 - 2]));
    // fall through
  case 2:
    ComputeTwoDiagTripleProduct(&(mm[(row_ct - 2) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct, &(dest[(row_ct - 2) * row_ctavp1]), &(dest[(row_ct - 1) * row_ctavp1 - 1]), &(dest[(row_ct - 1) * row_ctavp1]));
    break;
  case 1:
    dest[(row_ct - 1) * row_ctavp1] = TripleProduct(&(mm[(row_ct - 1) * col_ctav]), &(mm[(row_ct - 1) * col_ctav]), vv, col_ct);
  }
}

void CholeskyDecomposition(const float* aa, uint32_t predictor_ct, float* ll) {
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  const uintptr_t predictor_ctavp1 = predictor_ctav + 1;
  for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
    float fxx = aa[row_idx * predictor_ctavp1];
    float* ll_row_iter = &(ll[row_idx * predictor_ctav]);
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      const float fyy = (*ll_row_iter++);
      fxx -= fyy * fyy;
    }
    float fyy;
    if (fxx >= S_CAST(float, 0.0)) {
      fyy = sqrtf(fxx);
    } else {
      fyy = S_CAST(float, 1e-6);
    }
    ll[row_idx * predictor_ctavp1] = fyy;
    fyy = S_CAST(float, 1.0) / fyy;  // now 1.0 / L[j][j]
    for (uint32_t row_idx2 = row_idx + 1; row_idx2 != predictor_ct; ++row_idx2) {
      float fxx2 = aa[row_idx2 * predictor_ctav + row_idx];
      float* ll_row_iter2 = &(ll[row_idx * predictor_ctav]);
      float* ll_row_iter3 = &(ll[row_idx2 * predictor_ctav]);
      for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
        fxx2 -= (*ll_row_iter2++) * (*ll_row_iter3++);
      }
      ll[row_idx2 * predictor_ctav + row_idx] = fxx2 * fyy;
    }
  }
}

void SolveLinearSystem(const float* ll, const float* yy, uint32_t predictor_ct, float* xx) {
  // Finds x such that y = L(L^T)x, via forward and backward substitution
  //
  // might want to use this in NOLAPACK case only, since we can now produce
  // 32-bit Linux builds with statically linked LAPACK
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  for (uint32_t row_idx = 0; row_idx != predictor_ct; ++row_idx) {
    float fxx = yy[row_idx];
    const float* ll_row_iter = &(ll[row_idx * predictor_ctav]);
    float* xx_iter = xx;
    for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
      fxx -= (*ll_row_iter++) * (*xx_iter++);
    }
    *xx_iter = fxx / (*ll_row_iter);
  }
  for (uint32_t col_idx = predictor_ct; col_idx; ) {
    float fxx = xx[--col_idx];
    float* xx_iter = &(xx[predictor_ct - 1]);
    for (uint32_t row_idx = predictor_ct - 1; row_idx > col_idx; --row_idx) {
      fxx -= ll[row_idx * predictor_ctav + col_idx] * (*xx_iter--);
    }
    *xx_iter = fxx / ll[col_idx * (predictor_ctav + 1)];
  }
}

BoolErr LogisticRegression(const float* yy, const float* xx, uint32_t sample_ct, uint32_t predictor_ct, float* coef, uint32_t* is_unfinished_ptr, float* ll, float* pp, float* vv, float* hh, float* grad, float* dcoef) {
  // Similar to first part of logistic.cpp fitLM(), but incorporates changes
  // from Pascal Pons et al.'s TopCoder code.
  //
  // Preallocated buffers (initial contents irrelevant):
  // vv    = sample variance buffer
  // hh    = hessian matrix buffer, predictor_ct^2, rows vector-aligned
  // grad  = gradient buffer Y[] (length predictor_ct)
  // dcoef = current coefficient change buffer (length predictor_ct)
  //
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         vector-aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype; trailing elements must be zeroed out
  //
  // Input/output:
  // coef  = starting point, overwritten with logistic regression betas.  Must
  //         be vector-16-byte waligned.
  //
  // Outputs:
  // ll    = cholesky decomposition matrix, predictor_ct^2, rows vector-aligned
  // pp    = final likelihoods minus Y[] (not currently used by callers)
  //
  // Returns 0 on success, 1 on convergence failure.
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  float min_delta_coef = 1e9;

  ZeroFArr(predictor_ct * predictor_ctav, ll);
  for (uint32_t iteration = 0; ; ++iteration) {
    // P[i] = \sum_j X[i][j] * coef[j];
    ColMajorFmatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);
    // Suppose categorical covariates are represented as
    // categorical-covariate-major uint16_t* kk, indicating for each sample
    // which raw covariate index is 1 (with one "fallow" covariate index at the
    // end).
    // Then the above expression becomes
    // P[i] = \sum_j^{regular} X[i][j] * coef[j] +
    //        \sum_j^{cats} coef[K[i][j]]

    // P[i] = 1 / (1 + exp(-P[i]));
    LogisticSse(sample_ct, pp);

    // V[i] = P[i] * (1 - P[i]);
    // P[i] -= Y[i];
    ComputeVAndPMinusY(yy, sample_ct, pp, vv);

    // Possible categorical optimizations:
    // 1. skip terms between different categories of the same covariate
    // 2. all same-category terms within the same covariate can be handled with
    //    a single loop over the samples
    // 3. terms involving a regular covariate and a categorical covariate can
    //    be handled with one loop over the samples; covers all categories
    // 4. similarly, one loop over the samples is enough to update all category
    //    pairs between two categorical covariates
    ComputeHessian(xx, vv, sample_ct, predictor_ct, hh);

    // grad = X^T P
    // Separate categorical loop also possible here
    MultMatrixDxnVectN(xx, pp, sample_ct, predictor_ct, grad);

    CholeskyDecomposition(hh, predictor_ct, ll);

    // ZeroFArr(predictor_ct, dcoef);
    SolveLinearSystem(ll, grad, predictor_ct, dcoef);

    float delta_coef = 0.0;
    for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
      const float cur_dcoef = dcoef[pred_idx];
      delta_coef += fabsf(cur_dcoef);
      coef[pred_idx] -= cur_dcoef;
    }
    if (delta_coef < min_delta_coef) {
      min_delta_coef = delta_coef;
    }
    if (delta_coef != delta_coef) {
      return 1;
    }
    if (iteration > 3) {
      if (((delta_coef > S_CAST(float, 20.0)) && (delta_coef > 2 * min_delta_coef)) || ((iteration > 6) && (fabsf(S_CAST(float, 1.0) - delta_coef) < S_CAST(float, 1e-3)))) {
        return 1;
      }
      if (iteration > 13) {
        // If fabsf(any coefficient) > 8e3, this is almost certainly a form of
        // convergence failure that didn't get caught by the
        // (delta_coef > 20.0) check due to a precision quirk.  (8e3 threshold
        // ~= 1e-4 * 2^23, since floats have 23 bits of precision)
        for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
          if (fabsf(coef[pred_idx]) > S_CAST(float, 8e3)) {
            return 1;
          }
        }
        *is_unfinished_ptr = 1;
        return 0;
      }
    }
    // Pons reported that 1.1e-3 was dangerous, so I agree with the decision to
    // tighten this threshold from 1e-3 to 1e-4.
    if (delta_coef < S_CAST(float, 1e-4)) {
      // Be more conservative in throwing out results when we don't hit the
      // iteration limit.
      for (uint32_t pred_idx = 0; pred_idx != predictor_ct; ++pred_idx) {
        if (fabsf(coef[pred_idx]) > S_CAST(float, 6e4)) {
          return 1;
        }
      }
      return 0;
    }
  }
}

#ifdef __LP64__
// tmpNxK, interpreted as column-major, is sample_ct x predictor_ct
// X, interpreted as column-major, is also sample_ct x predictor_ct
// Hdiag[i] = V[i] (\sum_j tmpNxK[i][j] X[i][j])
void FirthComputeWeights(const float* yy, const float* xx, const float* pp, const float* vv, const float* tmpnxk, uint32_t predictor_ct, __maybe_unused uint32_t sample_ct, uint32_t sample_ctav, float* ww) {
  const VecF half = VCONST_F(0.5);
  for (uint32_t sample_offset = 0; sample_offset < sample_ctav; sample_offset += kFloatPerFVec) {
    VecF dotprods = vecf_setzero();
    const float* xx_row = &(xx[sample_offset]);
    const float* tmpnxk_row = &(tmpnxk[sample_offset]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const VecF cur_xx = *R_CAST(const VecF*, &(xx_row[pred_uidx * sample_ctav]));
      const VecF cur_tmpnxk = *R_CAST(const VecF*, &(tmpnxk_row[pred_uidx * sample_ctav]));
      dotprods = dotprods + cur_xx * cur_tmpnxk;
    }
    // Can handle categorical covariates in a separate loop here, and load the
    // dotprods increment into a union, etc.
    const VecF cur_weights = *R_CAST(const VecF*, &(vv[sample_offset]));
    const VecF cur_pis = *R_CAST(const VecF*, &(pp[sample_offset]));
    const VecF cur_yy = *R_CAST(const VecF*, &(yy[sample_offset]));
    const VecF half_minus_cur_pis = half - cur_pis;
    const VecF yy_minus_cur_pis = cur_yy - cur_pis;
    const VecF second_term = half_minus_cur_pis * (cur_weights * dotprods);
    const VecF cur_wws = yy_minus_cur_pis + second_term;
    *R_CAST(VecF*, &(ww[sample_offset])) = cur_wws;
  }
}
#else
void FirthComputeWeights(const float* yy, const float* xx, const float* pp, const float* vv, const float* tmpnxk, uint32_t predictor_ct, uint32_t sample_ct, uint32_t sample_ctav, float* ww) {
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    float dotprod = 0.0;
    const float* xx_row = &(xx[sample_idx]);
    const float* tmpnxk_row = &(tmpnxk[sample_idx]);
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      dotprod += xx_row[pred_uidx * sample_ctav] * tmpnxk_row[pred_uidx * sample_ctav];
    }
    const float cur_weight = vv[sample_idx];
    const float cur_pi = pp[sample_idx];
    ww[sample_idx] = (yy[sample_idx] - cur_pi) + (S_CAST(float, 0.5) - cur_pi) * cur_weight * dotprod;
  }
}
#endif

BoolErr FirthRegression(const float* yy, const float* xx, uint32_t sample_ct, uint32_t predictor_ct, float* coef, uint32_t* is_unfinished_ptr, float* hh, double* half_inverted_buf, MatrixInvertBuf1* inv_1d_buf, double* dbl_2d_buf, float* pp, float* vv, float* grad, float* dcoef, float* ww, float* tmpnxk_buf) {
  // This is a port of Georg Heinze's logistf R function, adapted to use many
  // of plink 1.9's optimizations; see
  //   http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/fllogistf/
  //
  // Preallocated buffers (initial contents irrelevant):
  // half_inverted_buf, inv_1d_buf, dbl_2d_buf = for matrix inversion
  // pp    = likelihoods minus Y[] (not currently used by callers)
  // vv    = sample variance buffer
  // grad  = gradient buffer (length predictor_ct)
  // dcoef = current coefficient change buffer (length predictor_ct)
  // ww    = Firth-adjusted scores, sample_ct
  //
  // Inputs:
  // xx    = covariate (and usually genotype) matrix, covariate-major, rows are
  //         vector-aligned, trailing row elements must be zeroed out
  // yy    = case/control phenotype
  //
  // Input/output:
  // coef  = starting point, overwritten with logistic regression betas.  Must
  //         be vector-aligned.
  //
  // Outputs:
  // hh    = variance-covariance matrix buffer, predictor_ct^2, rows
  //         vector-aligned.  (spends some time as pre-inversion Hessian matrix
  //         too)
  //
  // Returns 0 on success, 1 on convergence failure.
  // is_unfinished assumed to be initialized to 0, and is set to 1 if we hit
  // the iteration limit without satisfying other convergence criteria.
  const uintptr_t predictor_ctav = RoundUpPow2(predictor_ct, kFloatPerFVec);
  const uintptr_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  uint32_t is_last_iter = 0;

  // pull these out of the start of the loop, since they happen again in the
  // log-likelihood update
  // P[i] = \sum_j coef[j] * X[i][j];
  // categorical optimization possible here
  ColMajorFmatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);
  // P[i] = 1 / (1 + exp(-P[i]));
  LogisticSse(sample_ct, pp);
  // V[i] = P[i] * (1 - P[i]);
  ComputeV(pp, sample_ct, vv);
  // P[i] -= Y[i] NOT done here

  // hessian = X diag(V) X'
  // note that only lower triangle is filled here
  ComputeHessian(xx, vv, sample_ct, predictor_ct, hh);

  // we shouldn't need to compute the log directly, since underflow <->
  // regression failure, right?  check this.
  if (InvertSymmdefFmatrixFirstHalf(predictor_ct, predictor_ctav, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf)) {
    return 1;
  }
  double dethh = HalfSymmInvertedDet(half_inverted_buf, inv_1d_buf, predictor_ct);
  double loglik = ComputeLoglik(yy, pp, sample_ct);
  // printf("loglik: %g\n", loglik);
  loglik += 0.5 * log(dethh);

  // bugfix (4 Nov 2017): grad[] trailing elements must be zeroed out
  ZeroFArr(predictor_ctav - predictor_ct, &(grad[predictor_ct]));

  // start with 80% of most logistf convergence defaults (some reduction is
  // appropriate to be consistent with single-precision arithmetic).  (Update,
  // 9 Apr 2020: max_iter now matches logistf default since we've been shown a
  // concrete example where it matters in a way that isn't masked by the
  // single- vs. double-precision difference.)
  //
  // see also the hs_bail condition: if we ever try all five halfsteps, when
  // dcoef_max and grad_max aren't that far from the normal convergence
  // conditions, it's probably pointless to continue with single-precision
  // arithmetic.  (possible todo: use a fully-double-precision routine to
  // finish the job when that happens.)
  const uint32_t max_iter_m1 = 24;
  const float gconv = S_CAST(float, 0.0001);
  const float xconv = S_CAST(float, 0.0001);
  const double lconv = 0.0001;
  uint32_t hs_bail = 0;
  for (uint32_t iter_idx = 0; ; ++iter_idx) {
    InvertSymmdefFmatrixSecondHalf(predictor_ct, predictor_ctav, half_inverted_buf, hh, inv_1d_buf, dbl_2d_buf);
    if (is_last_iter) {
      return 0;
    }
    // bugfix (13 Oct 2017): trailing elements of hh[] rows can't be arbitrary
    // for later MultMatrixDxnVectN() call
    ReflectFmatrix0(predictor_ct, predictor_ctav, hh);

    // categorical optimization possible here
    ColMajorFmatrixMultiplyStrided(xx, hh, sample_ct, sample_ctav, predictor_ct, predictor_ctav, predictor_ct, sample_ctav, tmpnxk_buf);

    FirthComputeWeights(yy, xx, pp, vv, tmpnxk_buf, predictor_ct, sample_ct, sample_ctav, ww);

    // trailing elements of ww can't be nan for MultMatrixDxnVectN()
    ZeroFArr(sample_ctav - sample_ct, &(ww[sample_ct]));

    // gradient (Ustar in logistf) = X' W
    // categorical optimization possible here
    MultMatrixDxnVectN(xx, ww, sample_ct, predictor_ct, grad);
    float grad_max = 0.0;
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const float abs_grad_cur = fabsf(grad[pred_uidx]);
      if (abs_grad_cur > grad_max) {
        grad_max = abs_grad_cur;
      }
    }

    // dcoef := hh * grad (note that hh is inverted already)
    MultMatrixDxnVectN(hh, grad, predictor_ct, predictor_ct, dcoef);

    float dcoef_max = 0.0;
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      const float abs_dcoef_cur = fabsf(dcoef[pred_uidx]);
      if (abs_dcoef_cur > dcoef_max) {
        dcoef_max = abs_dcoef_cur;
      }
    }
    const float maxstep = 5.0;
    if (dcoef_max > maxstep) {
      const float scaling_factor = maxstep / dcoef_max;
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        dcoef[pred_uidx] *= scaling_factor;
      }
      dcoef_max = maxstep;
    }
    for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
      coef[pred_uidx] += dcoef[pred_uidx];
    }
    const uint32_t delta_and_grad_converged = (dcoef_max <= xconv) && (grad_max < gconv);
    const double loglik_old = loglik;
    double loglik_thresh = loglik_old;
    if (delta_and_grad_converged) {
      // on the last iteration, we would frequently try all 5 halfsteps when
      // the log-likelihood change was effectively random due to floating point
      // error.  detect this and exit the loop earlier.
      loglik_thresh -= 0.999999 * lconv;
    }

    uint32_t maxhs = 5;
    uint32_t halfstep_idx = 1;
    while (1) {
      // categorical optimization possible here
      ColMajorFmatrixVectorMultiplyStrided(xx, coef, sample_ct, sample_ctav, predictor_ct, pp);

      LogisticSse(sample_ct, pp);
      loglik = ComputeLoglik(yy, pp, sample_ct);
      ComputeV(pp, sample_ct, vv);
      ComputeHessian(xx, vv, sample_ct, predictor_ct, hh);
      if (InvertSymmdefFmatrixFirstHalf(predictor_ct, predictor_ctav, hh, half_inverted_buf, inv_1d_buf, dbl_2d_buf)) {
        return 1;
      }
      dethh = HalfSymmInvertedDet(half_inverted_buf, inv_1d_buf, predictor_ct);
      loglik += 0.5 * log(dethh);
      if (halfstep_idx > maxhs) {
        break;
      }
      if (loglik >= loglik_thresh) {
        if (loglik >= loglik_old) {
          break;
        }
        maxhs = halfstep_idx;
      } else if (halfstep_idx == maxhs) {
        if ((dcoef_max < S_CAST(float, 0.001)) && (grad_max < S_CAST(float, 0.05)) && (loglik >= loglik_old - lconv)) {
          // we've converged as much as we can with single-precision
          // arithmetic, and now we're flailing around.  don't even take the
          // 2^{-maxhs} step, undo it all and bail.
          // (0.001 and 0.05 constants can obviously be tuned; they were chosen
          // based on a test 500k sample/5 covariate regression.)
          --halfstep_idx;
          --maxhs;
          hs_bail = 1;
        }
      }
      const float multiplier = exp2f(-u31tof(halfstep_idx));
      for (uint32_t pred_uidx = 0; pred_uidx != predictor_ct; ++pred_uidx) {
        coef[pred_uidx] -= dcoef[pred_uidx] * multiplier;
      }
      ++halfstep_idx;
    }
    // printf("%.9g %.9g %g %g\n", loglik, loglik_old, dcoef_max, grad_max);
    const double loglik_change = loglik - loglik_old;
    if ((fabs(loglik_change) <= lconv) && (delta_and_grad_converged || hs_bail)) {
      is_last_iter = 1;
    } else if (iter_idx == max_iter_m1) {
      is_last_iter = 1;
      *is_unfinished_ptr = 1;
    }
  }
}

uintptr_t GetLogisticWorkspaceSize(uint32_t sample_ct, uint32_t biallelic_predictor_ct, uint32_t max_extra_allele_ct, uint32_t constraint_ct, uint32_t xmain_ct, uint32_t gcount_cc, uint32_t is_sometimes_firth) {
  // sample_ctav * max_predictor_ct < 2^31, and sample_ct >=
  // biallelic_predictor_ct, so no overflows?
  // could round everything up to multiples of 16 instead of 64
  const uint32_t max_predictor_ct = biallelic_predictor_ct + max_extra_allele_ct;
  const uint32_t sample_ctav = RoundUpPow2(sample_ct, kFloatPerFVec);
  const uint32_t max_predictor_ctav = RoundUpPow2(max_predictor_ct, kFloatPerFVec);
  // sample_nm, pheno_cc_nm, tmp_nm = sample_ctl words
  uintptr_t workspace_size = 3 * RoundUpPow2(BitCtToWordCt(sample_ct) * sizeof(intptr_t), kCacheline);

  // yy = sample_ctav floats
  workspace_size += RoundUpPow2(sample_ctav * sizeof(float), kCacheline);

  // xx = (max_predictor_ct + main_mutated + main_omitted) * sample_ctav floats
  workspace_size += RoundUpPow2((max_predictor_ct + xmain_ct) * sample_ctav * sizeof(float), kCacheline);

  // hh = max_predictor_ct * max_predictor_ctav floats
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ctav * sizeof(float), kCacheline);

  // pp, vv = sample_ctav floats
  workspace_size += 2 * RoundUpPow2(sample_ctav * sizeof(float), kCacheline);

  // coef, grad, dcoef = max_predictor_ctav floats
  workspace_size += 3 * RoundUpPow2(max_predictor_ctav * sizeof(float), kCacheline);

  // ll = max_predictor_ct * max_predictor_ctav floats
  // (technically not needed in pure-Firth case)
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ctav * sizeof(float), kCacheline);

  // semicomputed_biallelic_xtx
  workspace_size += RoundUpPow2(biallelic_predictor_ct * biallelic_predictor_ct * sizeof(double), kCacheline);

  // semicomputed_biallelic_corr_matrix
  workspace_size += RoundUpPow2((biallelic_predictor_ct - 1) * (biallelic_predictor_ct - 1) * sizeof(double), kCacheline);

  // semicomputed_biallelic_inv_corr_sqrts
  workspace_size += RoundUpPow2(biallelic_predictor_ct * sizeof(double), kCacheline);

  // inv_1d_buf
  workspace_size += RoundUpPow2(max_predictor_ct * kMatrixInvertBuf1CheckedAlloc, kCacheline);

  // dbl_2d_buf = max_predictor_ct * max_predictor_ctav floats, or VIF/Firth
  // dbl in practice, the latter value is never smaller due to the max(x, 7)
  workspace_size += RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);

  // a1_dosages, a1_case_dosages
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(double) * 2, kCacheline);

  // machr2_dosage_sums, machr2_dosage_ssqs
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(uint64_t) * 2, kCacheline);

  if (gcount_cc && max_extra_allele_ct) {
    // case_one_cts, case_two_cts
    workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(int32_t) * 2, kCacheline);
  }

  // predictor_dotprod_buf
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ct * sizeof(float), kCacheline);

  const uintptr_t other_2d_byte_ct = max_predictor_ct * MAXV(max_predictor_ct, 3) * sizeof(double);
  // inverse_corr_buf/half_inverted_buf
  workspace_size += RoundUpPow2(other_2d_byte_ct, kCacheline);

  if (is_sometimes_firth) {
    // ww = sample_ctav floats
    workspace_size += RoundUpPow2(sample_ctav * sizeof(float), kCacheline);

    // tmpnxk_buf = max_predictor_ct * sample_ctav floats
    workspace_size += RoundUpPow2(max_predictor_ct * sample_ctav * sizeof(float), kCacheline);
  }
  if (constraint_ct) {
    // tmphxs_buf, h_transpose_buf = constraint_ct * max_predictor_ctav floats
    workspace_size += 2 * RoundUpPow2(constraint_ct * max_predictor_ctav * sizeof(float), kCacheline);

    // inner_buf = constraint_ct * constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * constraint_ct * sizeof(float), kCacheline);

    // outer_buf = constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * sizeof(float), kCacheline);

    // constraints_con_major = constraint_ct * max_predictor_ct
    workspace_size += RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(float), kCacheline);
  }
  return workspace_size;
}


static const double kSmallDoubles[4] = {0.0, 1.0, 2.0, 3.0};

// could split this into per-variant and per-tested-allele parts, but only 12
// bytes (sample_obs_ct, multiallelic mach_r2) can go into the former; probably
// unimportant
typedef struct {
  // double beta;
  //   odds ratio = exp(beta)
  // double se;
  //   zval = beta / se
  //   width of asymptotic CI (beta units) = ci_zt * se
  //   T-statistic = zval
  //   pval = ZscoreToP(zval)

  uint32_t sample_obs_ct;

  uint32_t allele_obs_ct;
  double a1_dosage;

  uint16_t firth_fallback;
  uint16_t is_unfinished;
  uint32_t case_allele_obs_ct;
  double a1_case_dosage;

  double mach_r2;

  // case hom-ref, case ref-alt, case alt-alt, ctrl hom-ref, ...
  STD_ARRAY_DECL(uint32_t, 6, geno_hardcall_cts);
} LogisticAuxResult;

typedef struct {
  // double beta;
  // double se;
  //   zval = beta / se
  //   width of asymptotic CI = ci_zt * se
  //   T-statistic = zval
  //   pval = TstatToP(zval, sample_obs_ct - predictor_ct)

  uint32_t sample_obs_ct;

  uint32_t allele_obs_ct;
  double a1_dosage;

  double mach_r2;
} LinearAuxResult;

typedef struct GlmCtxStruct {
  const uintptr_t* variant_include;
  const ChrInfo* cip;
  const uintptr_t* allele_idx_offsets;
  const AlleleCode* omitted_alleles;
  uint32_t* subset_chr_fo_vidx_start;
  uintptr_t* sample_include;
  const uintptr_t* sample_include_x;
  const uintptr_t* sample_include_y;
  uint32_t* sample_include_cumulative_popcounts;
  uint32_t* sample_include_x_cumulative_popcounts;
  uint32_t* sample_include_y_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed;
  uintptr_t* parameter_subset;
  uintptr_t* parameter_subset_x;
  uintptr_t* parameter_subset_y;
  uintptr_t* joint_test_params;
  uintptr_t* joint_test_params_x;
  uintptr_t* joint_test_params_y;
  uint32_t variant_ct;
  uint32_t sample_ct;
  uint32_t sample_ct_x;
  uint32_t sample_ct_y;
  uint32_t max_extra_allele_ct;
  uint32_t covar_ct;
  uint32_t local_covar_ct;
  uint32_t covar_ct_x;
  uint32_t covar_ct_y;
  uint32_t constraint_ct;
  uint32_t constraint_ct_x;
  uint32_t constraint_ct_y;
  uint32_t tests_flag;
  GlmFlags glm_flags;
  uint32_t is_xchr_model_1;
  double max_corr;
  double vif_thresh;
  uintptr_t max_reported_test_ct;

  RegressionNmPrecomp* nm_precomp;
  RegressionNmPrecomp* nm_precomp_x;
  RegressionNmPrecomp* nm_precomp_y;

  uint32_t cur_block_variant_ct;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** thread_mhc;
  uintptr_t** dosage_presents;
  Dosage** dosage_mains;
  uint32_t* read_variant_uidx_starts;

  unsigned char** workspace_bufs;

  double* block_beta_se;

  uint64_t err_info;
} GlmCtx;

typedef struct GlmLogisticCtxStruct {
  GlmCtx *common;

  uintptr_t* pheno_cc;
  uintptr_t* pheno_x_cc;
  uintptr_t* pheno_y_cc;
  uintptr_t* gcount_case_interleaved_vec;
  uintptr_t* gcount_case_interleaved_vec_x;
  uintptr_t* gcount_case_interleaved_vec_y;
  const float* pheno_f;
  float* pheno_x_f;
  float* pheno_y_f;
  const float* covars_cmaj_f;
  float* covars_cmaj_x_f;
  float* covars_cmaj_y_f;
  uint16_t separation_found;
  uint16_t separation_found_x;
  uint16_t separation_found_y;
  float* local_covars_vcmaj_f[2];
  LogisticAuxResult* block_aux;
} GlmLogisticCtx;

THREAD_FUNC_DECL GlmLogisticThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GlmLogisticCtx* ctx = S_CAST(GlmLogisticCtx*, arg->sharedp->context);
  GlmCtx* common = ctx->common;

  PgenReader* pgrp = common->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = common->genovecs[tidx];
  pgv.dosage_present = nullptr;
  pgv.dosage_main = nullptr;
  if (common->dosage_presents) {
    pgv.dosage_present = common->dosage_presents[tidx];
    pgv.dosage_main = common->dosage_mains[tidx];
  }
  unsigned char* workspace_buf = common->workspace_bufs[tidx];
  const uintptr_t* variant_include = common->variant_include;
  const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
  const AlleleCode* omitted_alleles = common->omitted_alleles;
  const uintptr_t* sex_male_collapsed = common->sex_male_collapsed;
  const ChrInfo* cip = common->cip;
  const uint32_t* subset_chr_fo_vidx_start = common->subset_chr_fo_vidx_start;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const GlmFlags glm_flags = common->glm_flags;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t is_sometimes_firth = !(glm_flags & kfGlmNoFirth);
  const uint32_t is_always_firth = (glm_flags / kfGlmFirth) & 1;
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t joint_genotypic = (glm_flags / kfGlmGenotypic) & 1;
  const uint32_t joint_hethom = (glm_flags / kfGlmHethom) & 1;
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uint32_t domdev_present = joint_genotypic || joint_hethom;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t reported_pred_uidx_start = 1 - include_intercept;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const uintptr_t max_reported_test_ct = common->max_reported_test_ct;
  const uintptr_t local_covar_ct = common->local_covar_ct;
  const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
  // bugfix (20 Mar 2020): Also need to exclude dominant/recessive.
  const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!model_dominant) && (!model_recessive) && (!common->tests_flag) && (!add_interactions);
  uintptr_t max_sample_ct = MAXV(common->sample_ct, common->sample_ct_x);
  if (max_sample_ct < common->sample_ct_y) {
    max_sample_ct = common->sample_ct_y;
  }
  SetPgvThreadMhcNull(max_sample_ct, tidx, common->thread_mhc, &pgv);
  pgv.patch_01_ct = 0;
  pgv.patch_10_ct = 0;
  pgv.multidosage_sample_ct = 0;
  uint32_t variant_idx_offset = 0;
  uint32_t allele_ct = 2;
  uint32_t omitted_allele_idx = 0;
  uint32_t extra_regression_ct = 0;
  double main_dosage_sum = 0.0;
  double main_dosage_ssq = 0.0;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);

    double* beta_se_iter = common->block_beta_se;
    uintptr_t allele_bidx = variant_bidx;
    if (max_extra_allele_ct) {
      allele_bidx = variant_bidx + CountExtraAlleles(variant_include, allele_idx_offsets, common->read_variant_uidx_starts[0], common->read_variant_uidx_starts[tidx], 0);
    }
    if (beta_se_multiallelic_fused) {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * variant_bidx]);
    } else {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * allele_bidx]);
    }

    LogisticAuxResult* block_aux_iter = &(ctx->block_aux[allele_bidx]);
    const float* local_covars_iter = nullptr;
    if (local_covar_ct) {
      // &(nullptr[0]) is okay in C++, but undefined in C
      local_covars_iter = &(ctx->local_covars_vcmaj_f[parity][variant_bidx * max_sample_ct * local_covar_ct]);
    }
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = CountSortedSmallerU32(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      uint32_t cur_variant_bidx_end = subset_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
      if (cur_variant_bidx_end > variant_bidx_end) {
        cur_variant_bidx_end = variant_bidx_end;
      }
      const uint32_t is_x = (chr_idx == x_code);
      const uint32_t is_y = (chr_idx == y_code);
      const uint32_t is_nonx_haploid = (!is_x) && IsSet(cip->haploid_mask, chr_idx);
      const uintptr_t* cur_sample_include;
      const uint32_t* cur_sample_include_cumulative_popcounts;
      const uintptr_t* cur_pheno_cc;
      const uintptr_t* cur_gcount_case_interleaved_vec;
      const float* cur_pheno;
      const RegressionNmPrecomp* nm_precomp;
      const float* cur_covars_cmaj;
      const uintptr_t* cur_parameter_subset;
      const uintptr_t* cur_joint_test_params;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      uint32_t cur_constraint_ct;
      uint32_t cur_is_always_firth;
      if (is_y && common->sample_include_y) {
        cur_sample_include = common->sample_include_y;
        cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
        cur_pheno_cc = ctx->pheno_y_cc;
        cur_gcount_case_interleaved_vec = ctx->gcount_case_interleaved_vec_y;
        cur_pheno = ctx->pheno_y_f;
        nm_precomp = common->nm_precomp_y;
        cur_covars_cmaj = ctx->covars_cmaj_y_f;
        cur_parameter_subset = common->parameter_subset_y;
        cur_joint_test_params = common->joint_test_params_y;
        cur_sample_ct = common->sample_ct_y;
        cur_covar_ct = common->covar_ct_y;
        cur_constraint_ct = common->constraint_ct_y;
        cur_is_always_firth = is_always_firth || ctx->separation_found_y;
      } else if (is_x && common->sample_include_x) {
        cur_sample_include = common->sample_include_x;
        cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
        cur_pheno_cc = ctx->pheno_x_cc;
        cur_gcount_case_interleaved_vec = ctx->gcount_case_interleaved_vec_x;
        cur_pheno = ctx->pheno_x_f;
        nm_precomp = common->nm_precomp_x;
        cur_covars_cmaj = ctx->covars_cmaj_x_f;
        cur_parameter_subset = common->parameter_subset_x;
        cur_joint_test_params = common->joint_test_params_x;
        cur_sample_ct = common->sample_ct_x;
        cur_covar_ct = common->covar_ct_x;
        cur_constraint_ct = common->constraint_ct_x;
        cur_is_always_firth = is_always_firth || ctx->separation_found_x;
      } else {
        cur_sample_include = common->sample_include;
        cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
        cur_pheno_cc = ctx->pheno_cc;
        cur_gcount_case_interleaved_vec = ctx->gcount_case_interleaved_vec;
        cur_pheno = ctx->pheno_f;
        nm_precomp = common->nm_precomp;
        cur_covars_cmaj = ctx->covars_cmaj_f;
        cur_parameter_subset = common->parameter_subset;
        cur_joint_test_params = common->joint_test_params;
        cur_sample_ct = common->sample_ct;
        cur_covar_ct = common->covar_ct;
        cur_constraint_ct = common->constraint_ct;
        cur_is_always_firth = is_always_firth || ctx->separation_found;
      }
      const uint32_t sample_ctl = BitCtToWordCt(cur_sample_ct);
      const uint32_t sample_ctav = RoundUpPow2(cur_sample_ct, kFloatPerFVec);
      const uint32_t cur_case_ct = PopcountWords(cur_pheno_cc, sample_ctl);
      const uint32_t cur_biallelic_predictor_ct_base = 2 + domdev_present + cur_covar_ct * (1 + add_interactions * domdev_present_p1);
      uint32_t cur_biallelic_predictor_ct = cur_biallelic_predictor_ct_base;
      uint32_t literal_covar_ct = cur_covar_ct;
      if (cur_parameter_subset) {
        cur_biallelic_predictor_ct = PopcountWords(cur_parameter_subset, BitCtToWordCt(cur_biallelic_predictor_ct_base));
        literal_covar_ct = PopcountBitRange(cur_parameter_subset, 2 + domdev_present, 2 + domdev_present + cur_covar_ct);
      }
      const uint32_t max_predictor_ct = cur_biallelic_predictor_ct + max_extra_allele_ct;
      const uint32_t max_predictor_ctav = RoundUpPow2(max_predictor_ct, kFloatPerFVec);
      uint32_t reported_pred_uidx_biallelic_end;
      if (hide_covar) {
        if (!cur_parameter_subset) {
          reported_pred_uidx_biallelic_end = 2 + domdev_present;
        } else {
          reported_pred_uidx_biallelic_end = 1 + IsSet(cur_parameter_subset, 1) + domdev_present;
        }
      } else {
        reported_pred_uidx_biallelic_end = cur_biallelic_predictor_ct;
      }
      // nm_predictors_pmaj_buf may require up to two extra columns omitted
      // from the main regression.
      // 1. In the multiallelic dominant/recessive/hethom cases, the original
      //    genotype column does not appear in the regression, and we'd rather
      //    not reconstruct it from genovec, etc. when we need to swap it out
      //    for another allele, so we keep the original genotype in an extra
      //    column.
      //    To reduce code bloat, we now handle the biallelic cases in the same
      //    way; this is one of the more peripheral code paths so adding more
      //    complexity to speed it up is less justifiable.
      // 2. If --parameters excludes the main (possibly dominant/recessive)
      //    genotype column but does care about an interaction, we want a copy
      //    of what the main genotype column's contents would have been to
      //    refer to.
      const uint32_t main_omitted = cur_parameter_subset && (!IsSet(cur_parameter_subset, 1));
      const uint32_t main_mutated = model_dominant || model_recessive || joint_hethom;
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* pheno_cc_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* tmp_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      float* nm_pheno_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
      float* nm_predictors_pmaj_buf = S_CAST(float*, arena_alloc_raw_rd((max_predictor_ct + main_mutated + main_omitted) * sample_ctav * sizeof(float), &workspace_iter));
      float* coef_return = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ctav * sizeof(float), &workspace_iter));
      float* hh_return = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ctav * sizeof(float), &workspace_iter));
      float* pp_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
      float* sample_variance_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
      float* gradient_buf = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ctav * sizeof(float), &workspace_iter));
      float* dcoef_buf = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ctav * sizeof(float), &workspace_iter));
      float* cholesky_decomp_return = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ctav * sizeof(float), &workspace_iter));

      double* semicomputed_biallelic_xtx = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));
      // currently overallocates
      double* semicomputed_biallelic_corr_matrix = S_CAST(double*, arena_alloc_raw_rd((cur_biallelic_predictor_ct - 1) * (cur_biallelic_predictor_ct - 1) * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_inv_corr_sqrts = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));

      MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, arena_alloc_raw_rd(max_predictor_ct * kMatrixInvertBuf1CheckedAlloc, &workspace_iter));
      const uintptr_t dbl_2d_byte_ct = RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);
      double* dbl_2d_buf = S_CAST(double*, arena_alloc_raw(dbl_2d_byte_ct, &workspace_iter));
      double* a1_dosages = S_CAST(double*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(double) * 2, &workspace_iter));
      double* a1_case_dosages = &(a1_dosages[max_extra_allele_ct + 2]);
      uint64_t* machr2_dosage_sums = S_CAST(uint64_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(uint64_t) * 2, &workspace_iter));
      uint64_t* machr2_dosage_ssqs = &(machr2_dosage_sums[max_extra_allele_ct + 2]);
      uint32_t* case_one_cts = nullptr;
      uint32_t* case_two_cts = nullptr;
      if (cur_gcount_case_interleaved_vec && max_extra_allele_ct) {
        case_one_cts = S_CAST(uint32_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(int32_t) * 2, &workspace_iter));
        case_two_cts = &(case_one_cts[max_extra_allele_ct + 2]);
      }
      float* predictor_dotprod_buf = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ct * sizeof(float), &workspace_iter));
      const uintptr_t other_2d_byte_ct = RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 3) * sizeof(double), kCacheline);
      double* inverse_corr_buf = S_CAST(double*, arena_alloc_raw(other_2d_byte_ct, &workspace_iter));

      // these could use the same memory, but not a big deal, use the less
      // bug-prone approach for now
      // Firth-only
      float* score_buf = nullptr;
      float* tmpnxk_buf = nullptr;
      if (is_sometimes_firth) {
        score_buf = S_CAST(float*, arena_alloc_raw_rd(sample_ctav * sizeof(float), &workspace_iter));
        tmpnxk_buf = S_CAST(float*, arena_alloc_raw_rd(max_predictor_ct * sample_ctav * sizeof(float), &workspace_iter));
      }

      // joint test only
      float* tmphxs_buf = nullptr;
      float* h_transpose_buf = nullptr;
      float* inner_buf = nullptr;
      float* outer_buf = nullptr;
      float* cur_constraints_con_major = nullptr;
      if (cur_constraint_ct) {
        tmphxs_buf = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ctav * sizeof(float), &workspace_iter));
        h_transpose_buf = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ctav * sizeof(float), &workspace_iter));
        inner_buf = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * cur_constraint_ct * sizeof(float), &workspace_iter));
        outer_buf = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * sizeof(float), &workspace_iter));
        // bugfix (27 Jan 2019): forgot sizeof(float) here
        cur_constraints_con_major = S_CAST(float*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(float), &workspace_iter));
        ZeroFArr(cur_constraint_ct * max_predictor_ct, cur_constraints_con_major);
        const uint32_t first_joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
        cur_constraints_con_major[first_joint_test_idx] = 1.0;
        // Rest of this matrix must be updated later, since cur_predictor_ct
        // changes at multiallelic variants.
      }
      assert(S_CAST(uintptr_t, workspace_iter - workspace_buf) == GetLogisticWorkspaceSize(cur_sample_ct, cur_biallelic_predictor_ct, max_extra_allele_ct, cur_constraint_ct, main_mutated + main_omitted, cur_gcount_case_interleaved_vec != nullptr, is_sometimes_firth));
      const double cur_sample_ct_recip = 1.0 / u31tod(cur_sample_ct);
      const double cur_sample_ct_m1_recip = 1.0 / u31tod(cur_sample_ct - 1);
      const double* corr_inv = nullptr;
      if (nm_precomp) {
        memcpy(semicomputed_biallelic_xtx, nm_precomp->xtx_image, cur_biallelic_predictor_ct * cur_biallelic_predictor_ct * sizeof(double));
        corr_inv = nm_precomp->corr_inv;
        const uintptr_t nongeno_pred_ct = cur_biallelic_predictor_ct - domdev_present - 2;
        const uintptr_t nonintercept_biallelic_pred_ct = cur_biallelic_predictor_ct - 1;
        memcpy(semicomputed_biallelic_corr_matrix, nm_precomp->corr_image, nonintercept_biallelic_pred_ct * nonintercept_biallelic_pred_ct * sizeof(double));
        memcpy(&(semicomputed_biallelic_inv_corr_sqrts[domdev_present_p1]), nm_precomp->corr_inv_sqrts, nongeno_pred_ct * sizeof(double));
      }
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      // when this is set, the last fully-processed variant had no missing
      // genotypes, and if the current variant also has no missing genotypes we
      // may be able to skip reinitialization of most of
      // nm_predictors_pmaj_buf.
      // (todo: do we want to track prev_biallelic_nm?)
      uint32_t prev_nm = 0;

      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        if (allele_idx_offsets) {
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
          if (!beta_se_multiallelic_fused) {
            extra_regression_ct = allele_ct - 2;
          }
        }
        const uint32_t allele_ct_m2 = allele_ct - 2;
        const uint32_t expected_predictor_ct = cur_biallelic_predictor_ct + allele_ct_m2;
        PglErr reterr;
        if (!allele_ct_m2) {
          reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
        } else {
          reterr = PgrGetMD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, &pgv);
          // todo: proper multiallelic dosage support
        }
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto GlmLogisticThread_err;
        }
        ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
        GenoarrCountFreqsUnsafe(pgv.genovec, cur_sample_ct, genocounts);
        uint32_t missing_ct = genocounts[3];
        if (!missing_ct) {
          SetAllBits(cur_sample_ct, sample_nm);
        } else {
          GenoarrToNonmissing(pgv.genovec, cur_sample_ct, sample_nm);
          if (pgv.dosage_ct) {
            BitvecOr(pgv.dosage_present, sample_ctl, sample_nm);
            missing_ct = cur_sample_ct - PopcountWords(sample_nm, sample_ctl);
          }
        }
        if (omitted_alleles) {
          omitted_allele_idx = omitted_alleles[variant_uidx];
        }
        // Once sizeof(AlleleCode) > 1, we probably want to allocate this from
        // g_bigstack instead of the thread stack.
        uintptr_t const_alleles[DivUp(kPglMaxAltAlleleCt + 1, kBitsPerWord)];
        const uint32_t allele_ctl = DivUp(allele_ct, kBitsPerWord);
        ZeroWArr(allele_ctl, const_alleles);
        const uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
        const uint32_t nm_sample_ctl = BitCtToWordCt(nm_sample_ct);
        const uint32_t nm_sample_ctav = RoundUpPow2(nm_sample_ct, kFloatPerFVec);
        const uint32_t nm_sample_ct_rem = nm_sample_ctav - nm_sample_ct;
        // first predictor column: intercept
        if (!prev_nm) {
          for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
            nm_predictors_pmaj_buf[sample_idx] = 1.0;
          }
          ZeroFArr(nm_sample_ct_rem, &(nm_predictors_pmaj_buf[nm_sample_ct]));
        }
        // second predictor column: genotype
        float* genotype_vals = &(nm_predictors_pmaj_buf[nm_sample_ctav]);
        if (main_mutated || main_omitted) {
          genotype_vals = &(nm_predictors_pmaj_buf[expected_predictor_ct * nm_sample_ctav]);
        }
        CopyBitarrSubset(cur_pheno_cc, sample_nm, nm_sample_ct, pheno_cc_nm);
        const uint32_t nm_case_ct = PopcountWords(pheno_cc_nm, nm_sample_ctl);
        float* multi_start = nullptr;
        if (!allele_ct_m2) {
          if (omitted_allele_idx) {
            GenovecInvertUnsafe(cur_sample_ct, pgv.genovec);
            // ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
            }
            const uint32_t uii = genocounts[0];
            genocounts[0] = genocounts[2];
            genocounts[2] = uii;
          }
          uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
          uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
          if (!missing_ct) {
            GenoarrLookup16x4bx2(pgv.genovec, kSmallFloatPairs, nm_sample_ct, genotype_vals);
            if (pgv.dosage_ct) {
              uintptr_t sample_idx_base = 0;
              uintptr_t dosage_present_bits = pgv.dosage_present[0];
              for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                const uintptr_t sample_idx = BitIter1(pgv.dosage_present, &sample_idx_base, &dosage_present_bits);
                const uint32_t dosage_val = pgv.dosage_main[dosage_idx];
                // 32768 -> 2, 16384 -> 1, 0 -> 0
                genotype_vals[sample_idx] = kRecipDosageMidf * u31tof(dosage_val);
                dosage_sum += dosage_val;
                dosage_ssq += dosage_val * dosage_val;
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_idx);
                if (cur_geno && (cur_geno != 3)) {
                  const uintptr_t prev_val = cur_geno * kDosageMid;
                  dosage_sum -= prev_val;
                  dosage_ssq -= prev_val * prev_val;
                }
              }
            }
          } else {
            if (!pgv.dosage_ct) {
              GenoarrToFloatsRemoveMissing(pgv.genovec, kSmallFloats, cur_sample_ct, genotype_vals);
            } else {
              uintptr_t sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              uint32_t dosage_idx = 0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_midx);
                float cur_val;
                if (IsSet(pgv.dosage_present, sample_midx)) {
                  const uint32_t dosage_val = pgv.dosage_main[dosage_idx++];
                  cur_val = kRecipDosageMidf * u31tof(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                } else {
                  // cur_geno != 3 guaranteed
                  cur_val = kSmallFloats[cur_geno];
                }
                genotype_vals[sample_idx] = cur_val;
              }
            }
          }
          // Check for constant genotype column.
          // (Technically, we should recheck later in the chrX no-sex-covariate
          // --xchr-model 1 corner case.)
          if (!pgv.dosage_ct) {
            if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
              // bugfix (28 Mar 2020): didn't set the bit that actually
              // mattered last week...
              const_alleles[0] = 3;
            }
          } else if (pgv.dosage_ct == nm_sample_ct) {
            if (DosageIsConstant(dosage_sum, dosage_ssq, nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          }
          machr2_dosage_sums[1 - omitted_allele_idx] = dosage_sum;
          machr2_dosage_ssqs[1 - omitted_allele_idx] = dosage_ssq;
          machr2_dosage_sums[omitted_allele_idx] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[omitted_allele_idx] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
          if (cur_gcount_case_interleaved_vec) {
            // gcountcc
            STD_ARRAY_REF(uint32_t, 6) cur_geno_hardcall_cts = block_aux_iter->geno_hardcall_cts;
            GenoarrCountSubsetFreqs(pgv.genovec, cur_gcount_case_interleaved_vec, cur_sample_ct, cur_case_ct, R_CAST(STD_ARRAY_REF(uint32_t, 4), cur_geno_hardcall_cts));
            for (uint32_t geno_hardcall_idx = 0; geno_hardcall_idx != 3; ++geno_hardcall_idx) {
              cur_geno_hardcall_cts[3 + geno_hardcall_idx] = genocounts[geno_hardcall_idx] - cur_geno_hardcall_cts[geno_hardcall_idx];
            }
          }
        } else {
          // multiallelic.
          // Update (18 Mar 2020): If some but not all alleles have constant
          // dosages, we remove just those alleles from the regressions;
          // trim-alts is not necessary to see what's going on with the other
          // alleles.  To reduce parsing complexity, the number of output lines
          // is not affected by this; the ones corresponding to the constant
          // alleles have NA values.

          // dosage_ct == 0 temporarily guaranteed if we reach here.
          assert(!pgv.dosage_ct);
          multi_start = &(nm_predictors_pmaj_buf[(expected_predictor_ct - allele_ct_m2) * nm_sample_ctav]);
          ZeroU64Arr(allele_ct, machr2_dosage_sums);
          ZeroU64Arr(allele_ct, machr2_dosage_ssqs);
          // postpone multiply for now, since no multiallelic dosages
          // Use sums as ones[] and ssqs as twos[] for rarealts; transform to
          // actual sums/ssqs later.
          machr2_dosage_sums[0] = genocounts[1];
          machr2_dosage_ssqs[0] = genocounts[0];
          if (omitted_allele_idx) {
            // Main genotype column starts as REF.
            if (!missing_ct) {
              GenoarrLookup16x4bx2(pgv.genovec, kSmallInvFloatPairs, nm_sample_ct, genotype_vals);
            } else {
              GenoarrToFloatsRemoveMissing(pgv.genovec, kSmallInvFloats, cur_sample_ct, genotype_vals);
            }
          }
          uint32_t rare_allele_ct = allele_ct_m2;
          float* alt1_start = nullptr;
          float* rarealt_start = multi_start;
          if (omitted_allele_idx != 1) {
            if (omitted_allele_idx) {
              alt1_start = multi_start;
              ZeroFArr(nm_sample_ct_rem, &(alt1_start[nm_sample_ct]));
              rarealt_start = &(rarealt_start[nm_sample_ctav]);
              --rare_allele_ct;
            } else {
              alt1_start = genotype_vals;
            }
            if (!missing_ct) {
              GenoarrLookup16x4bx2(pgv.genovec, kSmallFloatPairs, nm_sample_ct, alt1_start);
            } else {
              GenoarrToFloatsRemoveMissing(pgv.genovec, kSmallFloats, cur_sample_ct, alt1_start);
            }
          }
          ZeroFArr(rare_allele_ct * nm_sample_ctav, rarealt_start);
          if (pgv.patch_01_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_01_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_01_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ctav + sample_idx] = 1.0;
                alt1_start[sample_idx] = 0.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ctav + sample_idx] = 1.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                alt1_start[sample_idx] = 0.0;
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                machr2_dosage_sums[allele_code] += 1;
                if (allele_code == omitted_allele_idx) {
                  continue;
                }
                const uint32_t cur_col = allele_code - 2 - (allele_code > omitted_allele_idx);
                rarealt_start[cur_col * nm_sample_ctav + sample_idx] = 1.0;
              }
            }
          }
          uintptr_t alt1_het_ct = genocounts[1] - pgv.patch_01_ct;
          if (pgv.patch_10_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_10_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_10_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] = 2.0;
                  alt1_start[sample_idx] = 0.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ctav + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] += S_CAST(float, 1.0);
                    alt1_start[sample_idx] = 0.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] = 2.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ctav + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ctav + sample_idx] += S_CAST(float, 1.0);
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  machr2_dosage_ssqs[ac0] += 1;
                  alt1_start[sample_idx] = 0.0;
                  if (ac0 != omitted_allele_idx) {
                    const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                    rarealt_start[ac0_col * nm_sample_ctav + sample_idx] = 2.0;
                  }
                } else {
                  machr2_dosage_sums[ac1] += 1;
                  if (ac1 != omitted_allele_idx) {
                    const uint32_t ac1_col = ac1 - 2 - (ac1 > omitted_allele_idx);
                    rarealt_start[ac1_col * nm_sample_ctav + sample_idx] = 1.0;
                  }
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    machr2_dosage_sums[ac0] += 1;
                    alt1_start[sample_idx] = 0.0;
                    if (ac0 != omitted_allele_idx) {
                      const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                      rarealt_start[ac0_col * nm_sample_ctav + sample_idx] += S_CAST(float, 1.0);
                    }
                  }
                }
              }
            }
          }
          machr2_dosage_sums[1] = alt1_het_ct;
          machr2_dosage_ssqs[1] = genocounts[2] - pgv.patch_10_ct;
          if (cur_gcount_case_interleaved_vec) {
            // gcountcc.  Need case-specific one_cts and two_cts for each
            // allele.
            STD_ARRAY_DECL(uint32_t, 4, case_hardcall_cts);
            GenoarrCountSubsetFreqs(pgv.genovec, cur_gcount_case_interleaved_vec, cur_sample_ct, cur_case_ct, case_hardcall_cts);
            ZeroU32Arr(allele_ct, case_one_cts);
            ZeroU32Arr(allele_ct, case_two_cts);
            uint32_t case_alt1_het_ct = case_hardcall_cts[1];
            case_one_cts[0] = case_alt1_het_ct;
            case_two_cts[0] = case_hardcall_cts[0];
            if (pgv.patch_01_ct) {
              uintptr_t sample_widx = 0;
              uintptr_t cur_bits = pgv.patch_01_set[0];
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t lowbit = BitIter1y(pgv.patch_01_set, &sample_widx, &cur_bits);
                if (cur_pheno_cc[sample_widx] & lowbit) {
                  const uint32_t allele_code = pgv.patch_01_vals[uii];
                  case_one_cts[allele_code] += 1;
                }
              }
              for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
                case_alt1_het_ct -= case_one_cts[allele_idx];
              }
            }
            uint32_t case_alt1_hom_ct = case_hardcall_cts[2];
            if (pgv.patch_10_ct) {
              uintptr_t sample_widx = 0;
              uintptr_t cur_bits = pgv.patch_10_set[0];
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t lowbit = BitIter1y(pgv.patch_10_set, &sample_widx, &cur_bits);
                if (cur_pheno_cc[sample_widx] & lowbit) {
                  const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                  const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                  --case_alt1_hom_ct;
                  if (ac0 == ac1) {
                    case_two_cts[ac0] += 1;
                  } else {
                    case_one_cts[ac1] += 1;
                    if (ac0 == 1) {
                      ++case_alt1_het_ct;
                    } else {
                      case_one_cts[ac0] += 1;
                    }
                  }
                }
              }
            }
            case_one_cts[1] = case_alt1_het_ct;
            case_two_cts[1] = case_alt1_hom_ct;
            uint32_t nonomitted_allele_idx = 0;
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if (allele_idx == omitted_allele_idx) {
                continue;
              }
              const uint32_t one_ct = machr2_dosage_sums[allele_idx];
              const uint32_t two_ct = machr2_dosage_ssqs[allele_idx];
              const uint32_t case_one_ct = case_one_cts[allele_idx];
              const uint32_t case_two_ct = case_two_cts[allele_idx];
              STD_ARRAY_REF(uint32_t, 6) dst = block_aux_iter[nonomitted_allele_idx].geno_hardcall_cts;
              dst[0] = nm_case_ct - case_one_ct - case_two_ct;
              dst[1] = case_one_ct;
              dst[2] = case_two_ct;
              dst[3] = nm_sample_ct - one_ct - two_ct - dst[0];
              dst[4] = one_ct - case_one_ct;
              dst[5] = two_ct - case_two_ct;
              ++nonomitted_allele_idx;
            }
          }
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            const uintptr_t one_ct = machr2_dosage_sums[allele_idx];
            const uintptr_t two_ct = machr2_dosage_ssqs[allele_idx];
            machr2_dosage_sums[allele_idx] = (one_ct + 2 * two_ct) * 0x4000LLU;
            machr2_dosage_ssqs[allele_idx] = (one_ct + 4LLU * two_ct) * 0x10000000LLU;
            if ((one_ct == nm_sample_ct) || (two_ct == nm_sample_ct) || ((!one_ct) && (!two_ct))) {
              SetBit(allele_idx, const_alleles);
            }
          }
        }
        ZeroFArr(nm_sample_ct_rem, &(genotype_vals[nm_sample_ct]));
        // usually need to save some of {sample_obs_ct, allele_obs_ct,
        // a1_dosage, case_allele_obs_ct, a1_case_dosage, mach_r2 even for
        // skipped variants
        // compute them all for now, could conditionally skip later
        uint32_t allele_obs_ct = nm_sample_ct * 2;
        uint32_t case_allele_obs_ct = nm_case_ct * 2;
        if (!is_x) {
          if (is_nonx_haploid) {
            allele_obs_ct = nm_sample_ct;
            case_allele_obs_ct = nm_case_ct;
            // everything is on 0..1 scale, not 0..2
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              genotype_vals[sample_idx] *= S_CAST(float, 0.5);
            }
            const uint32_t high_ct = nm_sample_ct * allele_ct_m2;
            for (uint32_t uii = 0; uii != high_ct; ++uii) {
              multi_start[uii] *= S_CAST(float, 0.5);
            }
          }
        } else {
          CopyBitarrSubset(sex_male_collapsed, sample_nm, nm_sample_ct, tmp_nm);
          const uintptr_t* male_nm = tmp_nm;
          const uint32_t nm_male_ct = PopcountWords(male_nm, nm_sample_ctl);
          if (is_xchr_model_1) {
            // special case: multiply male values by 0.5
            uintptr_t sample_idx_base = 0;
            uintptr_t male_nm_bits = male_nm[0];
            for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
              const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
              genotype_vals[sample_idx] *= S_CAST(float, 0.5);
              // could insert multiallelic loop here isntead, but I'm guessing
              // that's worse due to locality of writes?
            }
            for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
              float* cur_start = &(multi_start[extra_allele_idx * nm_sample_ctav]);
              sample_idx_base = 0;
              male_nm_bits = male_nm[0];
              for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
                const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
                cur_start[sample_idx] *= S_CAST(float, 0.5);
              }
            }
            allele_obs_ct -= nm_male_ct;
            case_allele_obs_ct -= PopcountWordsIntersect(pheno_cc_nm, male_nm, nm_sample_ctl);
          }
        }
        const double mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, allele_ct);
        uint32_t nonomitted_allele_idx = 0;
        for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
          if (allele_idx == omitted_allele_idx) {
            continue;
          }

          float* geno_col = genotype_vals;
          if (allele_idx > (!omitted_allele_idx)) {
            geno_col = &(nm_predictors_pmaj_buf[(expected_predictor_ct - (allele_ct - allele_idx) + (allele_idx < omitted_allele_idx)) * nm_sample_ctav]);
          }
          double a1_dosage = u63tod(machr2_dosage_sums[allele_idx]) * kRecipDosageMid;
          if (is_xchr_model_1) {
            // ugh.
            a1_dosage = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              a1_dosage += S_CAST(double, geno_col[sample_idx]);
            }
          } else {
            if (is_nonx_haploid) {
              a1_dosage *= 0.5;
            }
          }
          a1_dosages[allele_idx] = a1_dosage;

          // todo: shortcut if gcountcc computed and no dosages
          double a1_case_dosage = 0.0;
          uintptr_t sample_idx_base = 0;
          uintptr_t pheno_cc_nm_bits = pheno_cc_nm[0];
          for (uint32_t uii = 0; uii != nm_case_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(pheno_cc_nm, &sample_idx_base, &pheno_cc_nm_bits);
            a1_case_dosage += S_CAST(double, geno_col[sample_idx]);
          }
          a1_case_dosages[allele_idx] = a1_case_dosage;
          block_aux_iter[nonomitted_allele_idx].sample_obs_ct = nm_sample_ct;
          block_aux_iter[nonomitted_allele_idx].allele_obs_ct = allele_obs_ct;
          if (!allele_ct_m2) {
            // Need main_dosage_sum and main_dosage_ssq for now (probably move
            // this computation in-place later).
            if (is_xchr_model_1) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = 0.0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const double cur_dosage = S_CAST(double, geno_col[sample_idx]);
                main_dosage_ssq += cur_dosage * cur_dosage;
              }
            } else {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = u63tod(machr2_dosage_ssqs[allele_idx]) * kRecipDosageMidSq;
              if (is_nonx_haploid) {
                main_dosage_ssq *= 0.25;
              }
            }
          }
          block_aux_iter[nonomitted_allele_idx].a1_dosage = a1_dosage;

          // bugfix (4 Sep 2018): forgot to save this
          block_aux_iter[nonomitted_allele_idx].case_allele_obs_ct = case_allele_obs_ct;

          block_aux_iter[nonomitted_allele_idx].a1_case_dosage = a1_case_dosage;
          block_aux_iter[nonomitted_allele_idx].firth_fallback = 0;
          block_aux_iter[nonomitted_allele_idx].is_unfinished = 0;
          block_aux_iter[nonomitted_allele_idx].mach_r2 = mach_r2;
          ++nonomitted_allele_idx;
        }
        // Now free to skip the actual regression if there are too few samples,
        // or omitted allele corresponds to a zero-variance genotype column.
        // If another allele has zero variance but the omitted allele does not,
        // we now salvage as many alleles as we can.
        GlmErr glm_err = 0;
        if (nm_sample_ct <= expected_predictor_ct) {
          // reasonable for this to override CONST_ALLELE
          glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
        } else if (IsSet(const_alleles, omitted_allele_idx)) {
          glm_err = SetGlmErr0(kGlmErrcodeConstOmittedAllele);
        }
        if (glm_err) {
          if (missing_ct) {
            // covariates have not been copied yet, so we can't usually change
            // prev_nm from 0 to 1 when missing_ct == 0 (and there's little
            // reason to optimize the zero-covariate case)
            prev_nm = 0;
          }
          uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
          if (allele_ct_m2 && (beta_se_multiallelic_fused || (!hide_covar))) {
            reported_ct += allele_ct_m2;
          }
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            for (uint32_t uii = 0; uii != reported_ct; ++uii) {
              memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
              beta_se_iter[uii * 2 + 1] = -9.0;
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        } else {
          {
            double omitted_dosage = u63tod(allele_obs_ct);
            double omitted_case_dosage = u63tod(case_allele_obs_ct);
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if (allele_idx == omitted_allele_idx) {
                continue;
              }
              omitted_dosage -= a1_dosages[allele_idx];
              omitted_case_dosage -= a1_case_dosages[allele_idx];
            }
            a1_dosages[omitted_allele_idx] = omitted_dosage;
            a1_case_dosages[omitted_allele_idx] = omitted_case_dosage;
          }
          uint32_t parameter_uidx = 2 + domdev_present;
          float* nm_predictors_pmaj_istart = nullptr;
          // only need to do this part once per variant in multiallelic case
          float* nm_predictors_pmaj_iter = &(nm_predictors_pmaj_buf[nm_sample_ctav * (parameter_uidx - main_omitted)]);
          if (missing_ct || (!prev_nm)) {
            // fill phenotype
            uintptr_t sample_midx_base = 0;
            uintptr_t sample_nm_bits = sample_nm[0];
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
              nm_pheno_buf[sample_idx] = cur_pheno[sample_midx];
            }
            // bugfix (13 Oct 2017): must guarantee trailing phenotype values
            // are valid (exact contents don't matter since they are multiplied
            // by zero, but they can't be nan)
            ZeroFArr(nm_sample_ct_rem, &(nm_pheno_buf[nm_sample_ct]));

            // fill covariates
            for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx, ++parameter_uidx) {
              // strictly speaking, we don't need cur_covars_cmaj to be
              // vector-aligned
              if (cur_parameter_subset && (!IsSet(cur_parameter_subset, parameter_uidx))) {
                continue;
              }
              const float* cur_covar_col;
              if (covar_idx < local_covar_ct) {
                cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
              } else {
                cur_covar_col = &(cur_covars_cmaj[(covar_idx - local_covar_ct) * sample_ctav]);
              }
              sample_midx_base = 0;
              sample_nm_bits = sample_nm[0];
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                *nm_predictors_pmaj_iter++ = cur_covar_col[sample_midx];
              }
              ZeromovFArr(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
            }
            nm_predictors_pmaj_istart = nm_predictors_pmaj_iter;
            prev_nm = !missing_ct;
          } else {
            // bugfix (15 Aug 2018): this was not handling --parameters
            // correctly when a covariate was only needed as part of an
            // interaction
            parameter_uidx += cur_covar_ct;
            nm_predictors_pmaj_istart = &(nm_predictors_pmaj_iter[literal_covar_ct * nm_sample_ctav]);
          }
          const uint32_t const_allele_ct = PopcountWords(const_alleles, allele_ctl);
          if (const_allele_ct) {
            // Must delete constant-allele columns from nm_predictors_pmaj, and
            // shift later columns back.
            float* read_iter = genotype_vals;
            float* write_iter = genotype_vals;
            for (uint32_t read_allele_idx = 0; read_allele_idx != allele_ct; ++read_allele_idx) {
              if (read_allele_idx == omitted_allele_idx) {
                continue;
              }
              if (!IsSet(const_alleles, read_allele_idx)) {
                if (write_iter != read_iter) {
                  memcpy(write_iter, read_iter, nm_sample_ctav * sizeof(float));
                }
                if (write_iter == genotype_vals) {
                  write_iter = multi_start;
                } else {
                  write_iter = &(write_iter[nm_sample_ctav]);
                }
              }
              if (read_iter == genotype_vals) {
                read_iter = multi_start;
              } else {
                read_iter = &(read_iter[nm_sample_ctav]);
              }
            }
          }
          const uint32_t cur_predictor_ct = expected_predictor_ct - const_allele_ct;
          const uint32_t cur_predictor_ctav = RoundUpPow2(cur_predictor_ct, kFloatPerFVec);
          const uint32_t cur_predictor_ctavp1 = cur_predictor_ctav + 1;
          uint32_t nonconst_extra_regression_idx = UINT32_MAX;  // deliberate overflow
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            float* main_vals = &(nm_predictors_pmaj_buf[nm_sample_ctav]);
            float* domdev_vals = nullptr;
            uint32_t is_unfinished = 0;
            if (extra_regression_ct) {
              if (IsSet(const_alleles, extra_regression_idx + (extra_regression_idx >= omitted_allele_idx))) {
                glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
                goto GlmLogisticThread_skip_regression;
              }
              ++nonconst_extra_regression_idx;
              if (nonconst_extra_regression_idx) {
                float* swap_target = &(multi_start[(nonconst_extra_regression_idx - 1) * nm_sample_ctav]);
                for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
                  float fxx = genotype_vals[uii];
                  genotype_vals[uii] = swap_target[uii];
                  swap_target[uii] = fxx;
                }
              }
            }
            if (main_omitted) {
              // if main_mutated, this will be filled below
              // if not, this aliases genotype_vals
              main_vals = &(nm_predictors_pmaj_buf[(cur_predictor_ct + main_mutated) * nm_sample_ctav]);
            } else if (joint_genotypic || joint_hethom) {
              // in hethom case, do this before clobbering genotype data
              domdev_vals = &(main_vals[nm_sample_ctav]);
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                float cur_genotype_val = genotype_vals[sample_idx];
                if (cur_genotype_val > S_CAST(float, 1.0)) {
                  cur_genotype_val = S_CAST(float, 2.0) - cur_genotype_val;
                }
                domdev_vals[sample_idx] = cur_genotype_val;
              }
              ZeroFArr(nm_sample_ct_rem, &(domdev_vals[nm_sample_ct]));
            }
            if (model_dominant) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                float cur_genotype_val = genotype_vals[sample_idx];
                // 0..1..1
                if (cur_genotype_val > S_CAST(float, 1.0)) {
                  cur_genotype_val = 1.0;
                }
                main_vals[sample_idx] = cur_genotype_val;
              }
            } else if (model_recessive || joint_hethom) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                float cur_genotype_val = genotype_vals[sample_idx];
                // 0..0..1
                if (cur_genotype_val < S_CAST(float, 1.0)) {
                  cur_genotype_val = 0.0;
                } else {
                  cur_genotype_val -= S_CAST(float, 1.0);
                }
                main_vals[sample_idx] = cur_genotype_val;
              }
            }

            // fill interaction terms
            if (add_interactions) {
              nm_predictors_pmaj_iter = nm_predictors_pmaj_istart;
              for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
                const float* cur_covar_col;
                if (covar_idx < local_covar_ct) {
                  cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                } else {
                  cur_covar_col = &(cur_covars_cmaj[covar_idx * sample_ctav]);
                }
                if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                  uintptr_t sample_midx_base = 0;
                  uintptr_t sample_nm_bits = sample_nm[0];
                  for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                    const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                    *nm_predictors_pmaj_iter++ = main_vals[sample_idx] * cur_covar_col[sample_midx];
                  }
                  ZeromovFArr(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
                }
                ++parameter_uidx;
                if (domdev_present) {
                  if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                    uintptr_t sample_midx_base = 0;
                    uintptr_t sample_nm_bits = sample_nm[0];
                    for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                      const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                      *nm_predictors_pmaj_iter++ = domdev_vals[sample_idx] * cur_covar_col[sample_midx];
                    }
                    ZeromovFArr(nm_sample_ct_rem, &nm_predictors_pmaj_iter);
                  }
                  ++parameter_uidx;
                }
              }
            }
            if (corr_inv && prev_nm && (!allele_ct_m2)) {
              uintptr_t start_pred_idx = 0;
              if (!(model_dominant || model_recessive || joint_hethom)) {
                start_pred_idx = domdev_present + 2;
                semicomputed_biallelic_xtx[cur_predictor_ct] = main_dosage_sum;
                semicomputed_biallelic_xtx[cur_predictor_ct + 1] = main_dosage_ssq;
              }
              if (cur_predictor_ct > start_pred_idx) {
                ColMajorFvectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[nm_sample_ctav]), &(nm_predictors_pmaj_buf[start_pred_idx * nm_sample_ctav]), nm_sample_ct, nm_sample_ctav, cur_predictor_ct - start_pred_idx, &(predictor_dotprod_buf[start_pred_idx]));
                for (uint32_t uii = start_pred_idx; uii != cur_predictor_ct; ++uii) {
                  semicomputed_biallelic_xtx[cur_predictor_ct + uii] = S_CAST(double, predictor_dotprod_buf[uii]);
                }
              }
              if (domdev_present) {
                ColMajorFvectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[2 * nm_sample_ctav]), nm_predictors_pmaj_buf, nm_sample_ct, nm_sample_ctav, cur_predictor_ct, predictor_dotprod_buf);
                for (uint32_t uii = 0; uii != cur_predictor_ct; ++uii) {
                  semicomputed_biallelic_xtx[2 * cur_predictor_ct + uii] = S_CAST(double, predictor_dotprod_buf[uii]);
                }
                semicomputed_biallelic_xtx[cur_predictor_ct + 2] = semicomputed_biallelic_xtx[2 * cur_predictor_ct + 1];
              }
              glm_err = CheckMaxCorrAndVifNm(semicomputed_biallelic_xtx, corr_inv, cur_predictor_ct, domdev_present_p1, cur_sample_ct_recip, cur_sample_ct_m1_recip, max_corr, vif_thresh, semicomputed_biallelic_corr_matrix, semicomputed_biallelic_inv_corr_sqrts, dbl_2d_buf, &(dbl_2d_buf[2 * cur_predictor_ct]), &(dbl_2d_buf[3 * cur_predictor_ct]));
              if (glm_err) {
                goto GlmLogisticThread_skip_regression;
              }
            } else {
              glm_err = CheckMaxCorrAndVifF(&(nm_predictors_pmaj_buf[nm_sample_ctav]), cur_predictor_ct - 1, nm_sample_ct, nm_sample_ctav, max_corr, vif_thresh, predictor_dotprod_buf, dbl_2d_buf, inverse_corr_buf, inv_1d_buf);
              if (glm_err) {
                goto GlmLogisticThread_skip_regression;
              }
            }
            ZeroFArr(cur_predictor_ctav, coef_return);
            if (!cur_is_always_firth) {
              // Does any genotype column have zero case or zero control
              // dosage?  If yes, faster to skip logistic regression than
              // wait for convergence failure.
              for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                if (IsSet(const_alleles, allele_idx)) {
                  continue;
                }
                const double tot_dosage = a1_dosages[allele_idx];
                const double case_dosage = a1_case_dosages[allele_idx];
                if ((case_dosage == 0.0) || (case_dosage == tot_dosage)) {
                  if (is_sometimes_firth) {
                    goto GlmLogisticThread_firth_fallback;
                  }
                  glm_err = SetGlmErr1(kGlmErrcodeSeparation, allele_idx);
                  goto GlmLogisticThread_skip_regression;
                }
              }
              if (LogisticRegression(nm_pheno_buf, nm_predictors_pmaj_buf, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, cholesky_decomp_return, pp_buf, sample_variance_buf, hh_return, gradient_buf, dcoef_buf)) {
                if (is_sometimes_firth) {
                  ZeroFArr(cur_predictor_ctav, coef_return);
                  goto GlmLogisticThread_firth_fallback;
                }
                glm_err = SetGlmErr0(kGlmErrcodeLogisticConvergeFail);
                goto GlmLogisticThread_skip_regression;
              }
              // unlike FirthRegression(), hh_return isn't inverted yet, do
              // that here
              for (uint32_t pred_uidx = 0; pred_uidx != cur_predictor_ct; ++pred_uidx) {
                float* hh_inv_row = &(hh_return[pred_uidx * cur_predictor_ctav]);
                // ZeroFArr(cur_predictor_ct, gradient_buf);
                // gradient_buf[pred_uidx] = 1.0;
                // (y is gradient_buf, x is dcoef_buf)
                // SolveLinearSystem(cholesky_decomp_return, gradient_buf, cur_predictor_ct, hh_inv_row);
                // that works, but doesn't exploit the sparsity of y

                // hh_return does now have vector-aligned rows
                ZeroFArr(pred_uidx, hh_inv_row);

                float fxx = 1.0;
                for (uint32_t row_idx = pred_uidx; row_idx != cur_predictor_ct; ++row_idx) {
                  const float* ll_row = &(cholesky_decomp_return[row_idx * cur_predictor_ctav]);
                  for (uint32_t col_idx = pred_uidx; col_idx != row_idx; ++col_idx) {
                    fxx -= ll_row[col_idx] * hh_inv_row[col_idx];
                  }
                  hh_inv_row[row_idx] = fxx / ll_row[row_idx];
                  fxx = 0.0;
                }
                for (uint32_t col_idx = cur_predictor_ct; col_idx; ) {
                  fxx = hh_inv_row[--col_idx];
                  float* hh_inv_row_iter = &(hh_inv_row[cur_predictor_ct - 1]);
                  for (uint32_t row_idx = cur_predictor_ct - 1; row_idx > col_idx; --row_idx) {
                    fxx -= cholesky_decomp_return[row_idx * cur_predictor_ctav + col_idx] * (*hh_inv_row_iter--);
                  }
                  *hh_inv_row_iter = fxx / cholesky_decomp_return[col_idx * cur_predictor_ctavp1];
                }
              }
            } else {
              if (!is_always_firth) {
              GlmLogisticThread_firth_fallback:
                block_aux_iter[extra_regression_idx].firth_fallback = 1;
                if (allele_ct_m2 && beta_se_multiallelic_fused) {
                  for (uint32_t uii = 1; uii != allele_ct - 1; ++uii) {
                    block_aux_iter[uii].firth_fallback = 1;
                  }
                }
              }
              if (FirthRegression(nm_pheno_buf, nm_predictors_pmaj_buf, nm_sample_ct, cur_predictor_ct, coef_return, &is_unfinished, hh_return, inverse_corr_buf, inv_1d_buf, dbl_2d_buf, pp_buf, sample_variance_buf, gradient_buf, dcoef_buf, score_buf, tmpnxk_buf)) {
                glm_err = SetGlmErr0(kGlmErrcodeFirthConvergeFail);
                goto GlmLogisticThread_skip_regression;
              }
            }
            // validParameters() check
            for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
              const float hh_inv_diag_element = hh_return[pred_uidx * cur_predictor_ctavp1];
              if ((hh_inv_diag_element < S_CAST(float, 1e-20)) || (!isfinite_f(hh_inv_diag_element))) {
                glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                goto GlmLogisticThread_skip_regression;
              }
              // use sample_variance_buf[] to store diagonal square roots
              sample_variance_buf[pred_uidx] = sqrtf(hh_inv_diag_element);
            }
            sample_variance_buf[0] = sqrtf(hh_return[0]);
            for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
              const float cur_hh_inv_diag_sqrt = S_CAST(float, 0.99999) * sample_variance_buf[pred_uidx];
              const float* hh_inv_row_iter = &(hh_return[pred_uidx * cur_predictor_ctav]);
              const float* hh_inv_diag_sqrts_iter = sample_variance_buf;
              for (uint32_t pred_uidx2 = 0; pred_uidx2 != pred_uidx; ++pred_uidx2) {
                if ((*hh_inv_row_iter++) > cur_hh_inv_diag_sqrt * (*hh_inv_diag_sqrts_iter++)) {
                  glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                  goto GlmLogisticThread_skip_regression;
                }
              }
            }
            if (is_unfinished) {
              block_aux_iter[extra_regression_idx].is_unfinished = 1;
              if (allele_ct_m2 && beta_se_multiallelic_fused) {
                for (uint32_t uii = 1; uii != allele_ct - 1; ++uii) {
                  block_aux_iter[uii].is_unfinished = 1;
                }
              }
            }
            {
              double* beta_se_iter2 = beta_se_iter;
              for (uint32_t pred_uidx = reported_pred_uidx_start; pred_uidx != reported_pred_uidx_biallelic_end; ++pred_uidx) {
                // In the multiallelic-fused case, if the first allele is
                // constant, this writes the beta/se values for the first
                // nonconstant, non-omitted allele where the results for the
                // first allele belong.  We correct that at the end of this
                // block.
                *beta_se_iter2++ = S_CAST(double, coef_return[pred_uidx]);
                *beta_se_iter2++ = S_CAST(double, sample_variance_buf[pred_uidx]);
              }
              if (cur_constraint_ct) {
                *beta_se_iter2++ = 0.0;

                uint32_t joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 1.0;
                }
                double chisq;
                if (!LinearHypothesisChisqF(coef_return, cur_constraints_con_major, hh_return, cur_constraint_ct, cur_predictor_ct, cur_predictor_ctav, &chisq, tmphxs_buf, h_transpose_buf, inner_buf, inverse_corr_buf, inv_1d_buf, dbl_2d_buf, outer_buf)) {
                  *beta_se_iter2++ = chisq;
                } else {
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeRankDeficient);
                  memcpy(&(beta_se_iter2[-1]), &glm_err2, 8);
                  *beta_se_iter2++ = -9.0;
                }
                // next test may have different alt allele count
                joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 0.0;
                }
              }
              if (!const_allele_ct) {
                if (beta_se_multiallelic_fused || (!hide_covar)) {
                  for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
                    *beta_se_iter2++ = S_CAST(double, coef_return[cur_biallelic_predictor_ct + extra_allele_idx]);
                    *beta_se_iter2++ = S_CAST(double, sample_variance_buf[cur_biallelic_predictor_ct + extra_allele_idx]);
                  }
                }
              } else if (!beta_se_multiallelic_fused) {
                if (!hide_covar) {
                  // Need to insert some {CONST_ALLELE, -9} entries.
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                  const uint32_t cur_raw_allele_idx = extra_regression_idx + (extra_regression_idx >= omitted_allele_idx);
                  uint32_t extra_read_allele_idx = 0;
                  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    if ((allele_idx == omitted_allele_idx) || (allele_idx == cur_raw_allele_idx)) {
                      continue;
                    }
                    if (IsSet(const_alleles, allele_idx)) {
                      memcpy(beta_se_iter2, &glm_err2, 8);
                      beta_se_iter2[1] = -9.0;
                      beta_se_iter2 = &(beta_se_iter2[2]);
                    } else {
                      *beta_se_iter2++ = S_CAST(double, coef_return[cur_biallelic_predictor_ct + extra_read_allele_idx]);
                      *beta_se_iter2++ = S_CAST(double, sample_variance_buf[cur_biallelic_predictor_ct + extra_read_allele_idx]);
                      ++extra_read_allele_idx;
                    }
                  }
                }
              } else {
                const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                // Special-case first nonconst allele since it's positioned
                // discontinuously, and its BETA/SE may already be correctly
                // filled.
                uint32_t allele_idx = omitted_allele_idx? 0 : 1;
                if (IsSet(const_alleles, allele_idx)) {
                  memcpy(&(beta_se_iter[2 * include_intercept]), &glm_err2, 8);
                  beta_se_iter[2 * include_intercept + 1] = -9.0;
                  allele_idx = AdvTo0Bit(const_alleles, 1);
                  if (allele_idx == omitted_allele_idx) {
                    allele_idx = AdvTo0Bit(const_alleles, omitted_allele_idx + 1);
                  }
                  const uint32_t skip_ct = allele_idx - 1 - (allele_idx > omitted_allele_idx);
                  for (uint32_t uii = 0; uii != skip_ct; ++uii) {
                    memcpy(beta_se_iter2, &glm_err2, 8);
                    beta_se_iter2[1] = -9.0;
                    beta_se_iter2 = &(beta_se_iter2[2]);
                  }
                  *beta_se_iter2++ = S_CAST(double, coef_return[1]);
                  *beta_se_iter2++ = S_CAST(double, sample_variance_buf[1]);
                }
                ++allele_idx;
                uint32_t nonconst_allele_idx_m1 = 0;
                for (; allele_idx != allele_ct; ++allele_idx) {
                  if (allele_idx == omitted_allele_idx) {
                    continue;
                  }
                  if (!IsSet(const_alleles, allele_idx)) {
                    *beta_se_iter2++ = S_CAST(double, coef_return[cur_biallelic_predictor_ct + nonconst_allele_idx_m1]);
                    *beta_se_iter2++ = S_CAST(double, sample_variance_buf[cur_biallelic_predictor_ct + nonconst_allele_idx_m1]);
                    ++nonconst_allele_idx_m1;
                  } else {
                    memcpy(beta_se_iter2, &glm_err2, 8);
                    beta_se_iter2[1] = -9.0;
                    beta_se_iter2 = &(beta_se_iter2[2]);
                  }
                }
              }
            }
            while (0) {
            GlmLogisticThread_skip_regression:
              {
                uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
                if (allele_ct_m2 && (beta_se_multiallelic_fused || (!hide_covar))) {
                  reported_ct += allele_ct_m2;
                }
                for (uint32_t uii = 0; uii != reported_ct; ++uii) {
                  memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
                  beta_se_iter[uii * 2 + 1] = -9.0;
                }
              }
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        }
        block_aux_iter = &(block_aux_iter[allele_ct - 1]);
        if (local_covars_iter) {
          local_covars_iter = &(local_covars_iter[local_covar_ct * max_sample_ct]);
        }
      }
    }
    while (0) {
    GlmLogisticThread_err:
      UpdateU64IfSmaller(new_err_info, &common->err_info);
      break;
    }
    parity = 1 - parity;
    variant_idx_offset += cur_block_variant_ct;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

uint32_t GetBiallelicReportedTestCt(const uintptr_t* parameter_subset, GlmFlags glm_flags, uint32_t covar_ct, uint32_t tests_flag) {
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
  const uint32_t joint_test = domdev_present || tests_flag;

  if (hide_covar) {
    if (!parameter_subset) {
      return 1 + include_intercept + domdev_present + joint_test;
    }
    return include_intercept + domdev_present + joint_test + IsSet(parameter_subset, 1);
  }

  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t biallelic_predictor_ct_base = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
  uint32_t biallelic_predictor_ct = biallelic_predictor_ct_base;
  if (parameter_subset) {
    biallelic_predictor_ct = PopcountWords(parameter_subset, BitCtToWordCt(biallelic_predictor_ct_base));
  }
  return biallelic_predictor_ct + joint_test + include_intercept - 1;
}

BoolErr AllocAndInitReportedTestNames(const uintptr_t* parameter_subset, const char* const* covar_names, GlmFlags glm_flags, uint32_t covar_ct, uint32_t user_constraint_ct, const char*** cur_test_names_ptr) {
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t is_hethom = (glm_flags / kfGlmHethom) & 1;
  const uint32_t domdev_present = (glm_flags & kfGlmGenotypic) || is_hethom;
  char main_effect[4];
  if (model_dominant) {
    memcpy(main_effect, "DOMx", 4);
  } else if (model_recessive) {
    memcpy(main_effect, "RECx", 4);
  } else if (is_hethom) {
    memcpy(main_effect, "HOMx", 4);
  } else {
    memcpy(main_effect, "ADDx", 4);
  }
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t include_main_effect = (!parameter_subset) || IsSet(parameter_subset, 1);
  char domdev_str[8];
  uint32_t domdev_slen = 7;
  if (!is_hethom) {
    snprintf(domdev_str, 8, "DOMDEVx");
  } else {
    snprintf(domdev_str, 8, "HETx");
    domdev_slen = 4;
  }

  const uint32_t joint_test = domdev_present || user_constraint_ct;

  if (glm_flags & kfGlmHideCovar) {
    const uint32_t biallelic_reported_test_ct = include_intercept + include_main_effect + domdev_present + joint_test;
    char* test_name_buf_iter;
    if (unlikely(
            bigstack_alloc_kcp(biallelic_reported_test_ct, cur_test_names_ptr) ||
            bigstack_alloc_c(64, &test_name_buf_iter))) {
      return 1;
    }
    const char** cur_test_names = *cur_test_names_ptr;
    uint32_t write_idx = 0;
    if (include_intercept) {
      char* iter_next = memcpya_k(test_name_buf_iter, "INTERCEPT", 10);
      cur_test_names[write_idx++] = test_name_buf_iter;
      test_name_buf_iter = iter_next;
    }
    if (include_main_effect) {
      char* iter_next = memcpyax_k(test_name_buf_iter, main_effect, 3, '\0');
      cur_test_names[write_idx++] = test_name_buf_iter;
      test_name_buf_iter = iter_next;
    }
    if (domdev_present) {
      char* iter_next = memcpyax(test_name_buf_iter, domdev_str, domdev_slen - 1, '\0');
      cur_test_names[write_idx++] = test_name_buf_iter;
      test_name_buf_iter = iter_next;
    }
    if (joint_test) {
      if (user_constraint_ct) {
        char* iter_next = strcpya_k(test_name_buf_iter, "USER_");
        iter_next = u32toa(user_constraint_ct, iter_next);
        snprintf(iter_next, 64 - 10 - 4 - 8 - 15, "DF");
      } else {
        snprintf(test_name_buf_iter, 64 - 10 - 4 - 8, "GENO_2DF");
      }
      cur_test_names[write_idx++] = test_name_buf_iter;
    }
    assert(write_idx == biallelic_reported_test_ct);
    return 0;
  }
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  uint32_t biallelic_predictor_ct_base = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
  uint32_t biallelic_predictor_ct = biallelic_predictor_ct_base;
  if (parameter_subset) {
    biallelic_predictor_ct = PopcountWords(parameter_subset, BitCtToWordCt(biallelic_predictor_ct_base));
  }
  const uint32_t biallelic_reported_test_ct = biallelic_predictor_ct + joint_test + include_intercept - 1;
  uintptr_t test_name_buf_alloc = 64;
  if (add_interactions) {
    // don't bother optimizing this for parameter_subset case for now
    uintptr_t covar_name_total_blen = covar_ct;
    for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
      covar_name_total_blen += strlen(covar_names[covar_idx]);
    }
    // ADDx[covar name], etc.
    test_name_buf_alloc += 4 * covar_ct + covar_name_total_blen;
    if (is_hethom) {
      // HETx
      test_name_buf_alloc += 4 * covar_ct + covar_name_total_blen;
    } else if (domdev_present) {
      // DOMDEVx
      test_name_buf_alloc += 7 * covar_ct + covar_name_total_blen;
    }
  }
  char* test_name_buf_iter;
  if (unlikely(
          bigstack_alloc_kcp(biallelic_reported_test_ct, cur_test_names_ptr) ||
          bigstack_alloc_c(test_name_buf_alloc, &test_name_buf_iter))) {
    return 1;
  }
  const char** cur_test_names = *cur_test_names_ptr;
  uint32_t write_idx = 0;
  if (include_intercept) {
    char* iter_next = memcpya_k(test_name_buf_iter, "INTERCEPT", 10);
    cur_test_names[write_idx++] = test_name_buf_iter;
    test_name_buf_iter = iter_next;
  }
  if (include_main_effect) {
    char* iter_next = memcpyax_k(test_name_buf_iter, main_effect, 3, '\0');
    cur_test_names[write_idx++] = test_name_buf_iter;
    test_name_buf_iter = iter_next;
  }
  if (domdev_present) {
    char* iter_next = memcpyax(test_name_buf_iter, domdev_str, domdev_slen - 1, '\0');
    cur_test_names[write_idx++] = test_name_buf_iter;
    test_name_buf_iter = iter_next;
  }
  uint32_t pred_uidx = 2 + domdev_present;
  for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx, ++pred_uidx) {
    if (parameter_subset && (!IsSet(parameter_subset, pred_uidx))) {
      continue;
    }
    // just point to the existing string, its lifetime is sufficient
    cur_test_names[write_idx++] = covar_names[covar_idx];
  }
  if (add_interactions) {
    for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
      const char* cur_covar_name = covar_names[covar_idx];
      if ((!parameter_subset) || IsSet(parameter_subset, pred_uidx)) {
        char* iter_next = memcpya(test_name_buf_iter, main_effect, 4);
        iter_next = strcpyax(iter_next, cur_covar_name, '\0');
        cur_test_names[write_idx++] = test_name_buf_iter;
        test_name_buf_iter = iter_next;
      }
      ++pred_uidx;
      if (domdev_present) {
        if ((!parameter_subset) || IsSet(parameter_subset, pred_uidx)) {
          char* iter_next = memcpya(test_name_buf_iter, domdev_str, domdev_slen);
          iter_next = strcpyax(iter_next, cur_covar_name, '\0');
          cur_test_names[write_idx++] = test_name_buf_iter;
          test_name_buf_iter = iter_next;
        }
        ++pred_uidx;
      }
    }
  }
  if (joint_test) {
    if (user_constraint_ct) {
      char* iter_next = strcpya_k(test_name_buf_iter, "USER_");
      iter_next = u32toa(user_constraint_ct, iter_next);
      snprintf(iter_next, 64 - 10 - 4 - 8 - 15, "DF");
    } else {
      snprintf(test_name_buf_iter, 64 - 10 - 4 - 8, "GENO_2DF");
    }
    cur_test_names[write_idx++] = test_name_buf_iter;
  }
  assert(write_idx == biallelic_reported_test_ct);
  return 0;
}

typedef struct LocalCovarCoeffparseCtxStruct {
  uint32_t* sample_idx_order;
  uint32_t max_sample_ct;
  uint32_t cur_sample_ct;
  uint32_t tokens_per_sample;
  uint32_t local_covar_ct;
  uint32_t omit_last;
  uint32_t local_haps;
  uint32_t local_cat_ct;
  uint32_t local_cats_1based;
} LocalCovarCoeffparseCtx;

// Processes a single line's worth of payload.
PglErr LoadLocalCovarCoeffs(const LocalCovarCoeffparseCtx* ctx, const char* local_line_iter, uint32_t local_line_idx, float* local_covars_vcmaj_f_iter, double* local_covars_vcmaj_d_iter) {
  const uint32_t* sample_idx_order = ctx->sample_idx_order;
  const uint32_t max_sample_ct = ctx->max_sample_ct;
  const uint32_t cur_sample_ct = ctx->cur_sample_ct;
  const uint32_t tokens_per_sample = ctx->tokens_per_sample;
  const uint32_t local_covar_ct = ctx->local_covar_ct;
  const uint32_t omit_last = ctx->omit_last;
  const uint32_t local_haps = ctx->local_haps;
  const uint32_t local_cat_ct = ctx->local_cat_ct;
  const uint32_t local_cats_1based = ctx->local_cats_1based;
  const uint32_t max_cat_idx = local_cat_ct + local_cats_1based - 1;
  uint32_t sample_idx = 0;
  for (uint32_t local_sample_idx = 0; sample_idx != cur_sample_ct; ++local_sample_idx) {
    const uint32_t cur_sample_idx = sample_idx_order[local_sample_idx];
    if (cur_sample_idx == UINT32_MAX) {
      local_line_iter = NextTokenMult(local_line_iter, tokens_per_sample);
      if (unlikely(!local_line_iter)) {
        logputs("\n");
        logerrprintfww("Error: Fewer tokens than expected on line %u of --glm local-covar= file.\n", local_line_idx);
        return kPglRetMalformedInput;
      }
      continue;
    }
    if (local_cat_ct) {
      uint32_t cat_idx;
      if (unlikely(ScanmovUintCapped(max_cat_idx, &local_line_iter, &cat_idx) || (cat_idx < local_cats_1based))) {
        logputs("\n");
        logerrprintf("Error: Invalid category index on line %u of --glm local-covar= file.\n", local_line_idx);
        return kPglRetMalformedInput;
      }
      cat_idx -= local_cats_1based;
      local_line_iter = FirstNonTspace(FirstSpaceOrEoln(local_line_iter));
      if (!local_haps) {
        if (cat_idx != max_cat_idx) {
          const uint32_t offset = cat_idx * max_sample_ct + cur_sample_idx;
          if (local_covars_vcmaj_f_iter) {
            local_covars_vcmaj_f_iter[offset] = 1.0;
          } else {
            local_covars_vcmaj_d_iter[offset] = 1.0;
          }
        }
      } else {
        if (cat_idx != max_cat_idx) {
          const uint32_t offset = cat_idx * max_sample_ct + cur_sample_idx;
          if (local_covars_vcmaj_f_iter) {
            local_covars_vcmaj_f_iter[offset] = 0.5;
          } else {
            local_covars_vcmaj_d_iter[offset] = 0.5;
          }
        }
        if (unlikely(ScanmovUintCapped(max_cat_idx, &local_line_iter, &cat_idx) || (cat_idx < local_cats_1based))) {
          logputs("\n");
          logerrprintf("Error: Invalid category index on line %u of --glm local-covar= file.\n", local_line_idx);
          return kPglRetMalformedInput;
        }
        cat_idx -= local_cats_1based;
        local_line_iter = FirstNonTspace(FirstSpaceOrEoln(local_line_iter));
        if (cat_idx != max_cat_idx) {
          const uint32_t offset = cat_idx * max_sample_ct + cur_sample_idx;
          if (local_covars_vcmaj_f_iter) {
            local_covars_vcmaj_f_iter[offset] += S_CAST(float, 0.5);
          } else {
            local_covars_vcmaj_d_iter[offset] += 0.5;
          }
        }
      }
    } else {
      if (local_covars_vcmaj_f_iter) {
        float* local_covars_f_iter2 = &(local_covars_vcmaj_f_iter[cur_sample_idx]);
        for (uint32_t covar_idx = 0; covar_idx != local_covar_ct; ++covar_idx) {
          double dxx;
          local_line_iter = ScantokDouble(local_line_iter, &dxx);
          if (unlikely((!local_line_iter) || (fabs(dxx) > 3.4028235677973362e38))) {
            logputs("\n");
            logerrprintf("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
            return kPglRetMalformedInput;
          }
          *local_covars_f_iter2 = S_CAST(float, dxx);
          local_covars_f_iter2 = &(local_covars_f_iter2[max_sample_ct]);
          local_line_iter = FirstNonTspace(local_line_iter);
        }
      } else {
        double* local_covars_d_iter2 = &(local_covars_vcmaj_d_iter[cur_sample_idx]);
        for (uint32_t covar_idx = 0; covar_idx != local_covar_ct; ++covar_idx) {
          double dxx;
          local_line_iter = ScantokDouble(local_line_iter, &dxx);
          if (unlikely(!local_line_iter)) {
            logputs("\n");
            logerrprintf("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
            return kPglRetMalformedInput;
          }
          *local_covars_d_iter2 = dxx;
          local_covars_d_iter2 = &(local_covars_d_iter2[max_sample_ct]);
          local_line_iter = FirstNonTspace(local_line_iter);
        }
      }
      if (omit_last) {
        local_line_iter = FirstNonTspace(FirstSpaceOrEoln(local_line_iter));
      }
      if (local_haps) {
        if (local_covars_vcmaj_f_iter) {
          float* local_covars_f_iter2 = &(local_covars_vcmaj_f_iter[cur_sample_idx]);
          for (uint32_t covar_idx = 0; covar_idx != local_covar_ct; ++covar_idx) {
            double dxx;
            local_line_iter = ScantokDouble(local_line_iter, &dxx);
            if (unlikely((!local_line_iter) || (fabs(dxx) > 3.4028235677973362e38))) {
              logputs("\n");
              logerrprintf("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
              return kPglRetMalformedInput;
            }
            *local_covars_f_iter2 = S_CAST(float, (S_CAST(double, *local_covars_f_iter2) + dxx) * 0.5);
            local_covars_f_iter2 = &(local_covars_f_iter2[max_sample_ct]);
            local_line_iter = FirstNonTspace(local_line_iter);
          }
        } else {
          double* local_covars_d_iter2 = &(local_covars_vcmaj_d_iter[cur_sample_idx]);
          for (uint32_t covar_idx = 0; covar_idx != local_covar_ct; ++covar_idx) {
            double dxx;
            local_line_iter = ScantokDouble(local_line_iter, &dxx);
            if (unlikely(!local_line_iter)) {
              logputs("\n");
              logerrprintf("Error: Invalid or missing token on line %u of --glm local-covar= file.\n", local_line_idx);
              return kPglRetMalformedInput;
            }
            // may as well defend against overflow
            *local_covars_d_iter2 = (*local_covars_d_iter2) * 0.5 + dxx * 0.5;
            local_covars_d_iter2 = &(local_covars_d_iter2[max_sample_ct]);
            local_line_iter = FirstNonTspace(local_line_iter);
          }
        }
        if (omit_last) {
          local_line_iter = FirstNonTspace(FirstSpaceOrEoln(local_line_iter));
        }
      }
    }
    ++sample_idx;
  }
  return kPglRetSuccess;
}

PglErr ReadLocalCovarBlock(const GlmCtx* common, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t cur_block_variant_ct, uint32_t local_sample_ct, uint32_t local_cat_ct, TextStream* local_covar_txsp, uint32_t* local_line_idx_ptr, uint32_t* local_xy_ptr, float* local_covars_vcmaj_f_iter, double* local_covars_vcmaj_d_iter, uint32_t* local_sample_idx_order) {
  const ChrInfo* cip = common->cip;
  const uintptr_t* variant_include = common->variant_include;
  const uint32_t sample_ct = common->sample_ct;
  const uint32_t sample_ct_x = common->sample_ct_x;
  const uint32_t sample_ct_y = common->sample_ct_y;
  const uint32_t local_covar_ct = common->local_covar_ct;
  const GlmFlags flags = common->glm_flags;
  const uint32_t omit_last = (flags / kfGlmLocalOmitLast) & 1;
  const uint32_t local_haps = (flags / kfGlmLocalHaps) & 1;

  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
  if (max_sample_ct < sample_ct_y) {
    max_sample_ct = sample_ct_y;
  }
  LocalCovarCoeffparseCtx coeffparse_ctx;
  coeffparse_ctx.sample_idx_order = local_sample_idx_order;
  coeffparse_ctx.max_sample_ct = max_sample_ct;
  // cur_sample_ct filled a bit later
  coeffparse_ctx.tokens_per_sample = (local_cat_ct? 1 : (local_covar_ct + omit_last)) << local_haps;
  coeffparse_ctx.local_covar_ct = local_covar_ct;
  coeffparse_ctx.omit_last = omit_last;
  coeffparse_ctx.local_haps = local_haps;
  coeffparse_ctx.local_cat_ct = local_cat_ct;
  coeffparse_ctx.local_cats_1based = (flags / kfGlmLocalCats1based) & 1;
  uint32_t variant_bidx = 0;
  if (local_cat_ct) {
    // assert(local_covar_ct == local_cat_ct - 1);
    if (local_covars_vcmaj_f_iter) {
      ZeroFArr(local_covar_ct * max_sample_ct * S_CAST(uintptr_t, cur_block_variant_ct), local_covars_vcmaj_f_iter);
    } else {
      ZeroDArr(local_covar_ct * max_sample_ct * S_CAST(uintptr_t, cur_block_variant_ct), local_covars_vcmaj_d_iter);
    }
  }
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, variant_uidx_start, &variant_uidx_base, &cur_bits);
  uint32_t local_line_idx = *local_line_idx_ptr;
  while (variant_bidx < cur_block_variant_ct) {
    const uint32_t variant_uidx1 = BitIter1NoAdv(variant_include, &variant_uidx_base, &cur_bits);
    const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx1);
    const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
    const uint32_t chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    uint32_t cur_variant_bidx_end = cur_block_variant_ct;
    if (chr_end < variant_uidx_end) {
      cur_variant_bidx_end = variant_bidx + PopcountBitRange(variant_include, variant_uidx1, chr_end);
      assert(cur_variant_bidx_end <= cur_block_variant_ct);
    }
    const uint32_t is_x = (chr_idx == x_code);
    const uint32_t is_y = (chr_idx == y_code);
    const uintptr_t* cur_sample_include;
    const uint32_t* cur_sample_include_cumulative_popcounts;
    if (is_y && common->sample_include_y) {
      cur_sample_include = common->sample_include_y;
      cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
      coeffparse_ctx.cur_sample_ct = sample_ct_y;
    } else if (is_x && common->sample_include_x) {
      cur_sample_include = common->sample_include_x;
      cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
      coeffparse_ctx.cur_sample_ct = sample_ct_x;
    } else {
      cur_sample_include = common->sample_include;
      cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
      coeffparse_ctx.cur_sample_ct = sample_ct;
    }
    const uint32_t new_local_xy = is_x + 2 * is_y;
    if (new_local_xy != *local_xy_ptr) {
      for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
        const uint32_t cur_uidx = local_sample_uidx_order[uii];
        uint32_t cur_idx = UINT32_MAX;
        if ((cur_uidx != UINT32_MAX) && IsSet(cur_sample_include, cur_uidx)) {
          cur_idx = RawToSubsettedPos(cur_sample_include, cur_sample_include_cumulative_popcounts, cur_uidx);
        }
        local_sample_idx_order[uii] = cur_idx;
      }
      *local_xy_ptr = new_local_xy;
    }
    for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
      BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (!IsSet(local_variant_include, local_line_idx)) {
        uint32_t local_line_idx_target_m1 = AdvTo1Bit(local_variant_include, local_line_idx);
        PglErr reterr = TextSkipNz(local_line_idx_target_m1 - local_line_idx, local_covar_txsp);
        if (unlikely(reterr)) {
          if (reterr == kPglRetEof) {
            logputs("\n");
            logerrputs("Error: --glm local-covar= file has fewer lines than local-pvar= file.\n");
            return kPglRetInconsistentInput;
          }
          TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
          return reterr;
        }
        local_line_idx = local_line_idx_target_m1;
      }
      ++local_line_idx;
      PglErr reterr = kPglRetSuccess;
      char* local_covar_line_start = TextGet(local_covar_txsp);
      if (unlikely(!local_covar_line_start)) {
        if (!TextStreamErrcode2(local_covar_txsp, &reterr)) {
          logputs("\n");
          logerrputs("Error: --glm local-covar= file has fewer lines than local-pvar= file.\n");
          return kPglRetInconsistentInput;
        }
        TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
        return reterr;
      }
      reterr = LoadLocalCovarCoeffs(&coeffparse_ctx, local_covar_line_start, local_line_idx, local_covars_vcmaj_f_iter, local_covars_vcmaj_d_iter);
      if (unlikely(reterr)) {
        return reterr;
      }
      if (local_covars_vcmaj_f_iter) {
        local_covars_vcmaj_f_iter += max_sample_ct * local_covar_ct;
      } else {
        local_covars_vcmaj_d_iter += max_sample_ct * local_covar_ct;
      }
    }
  }
  *local_line_idx_ptr = local_line_idx;
  return kPglRetSuccess;
}

static inline void ZeroLocalCovarRows(uint32_t row_ct, uintptr_t row_width, uint32_t* variant_bidx_ptr, float** local_covars_vcmaj_f_iterp, double** local_covars_vcmaj_d_iterp) {
  *variant_bidx_ptr += row_ct;
  const uintptr_t elem_ct = row_ct * row_width;
  if (*local_covars_vcmaj_f_iterp) {
    ZeroFArr(elem_ct, *local_covars_vcmaj_f_iterp);
    *local_covars_vcmaj_f_iterp += elem_ct;
  } else {
    ZeroDArr(elem_ct, *local_covars_vcmaj_d_iterp);
    *local_covars_vcmaj_d_iterp += elem_ct;
  }
}

void DuplicateLocalCovarRow(uint32_t row_ct, uintptr_t row_width, uint32_t* variant_bidx_ptr, float** local_covars_vcmaj_f_iterp, double** local_covars_vcmaj_d_iterp) {
  if (*local_covars_vcmaj_f_iterp) {
    float* local_covars_vcmaj_f_iter = *local_covars_vcmaj_f_iterp;
    // Could also copy one row at a time, but I'd expect doubling to be
    // slightly more efficient.
    for (uint32_t row_idx = 1; row_idx < row_ct; row_idx *= 2) {
      uint32_t row_copy_ct = row_ct - row_idx;
      if (row_copy_ct > row_idx) {
        row_copy_ct = row_idx;
      }
      memcpy(&(local_covars_vcmaj_f_iter[row_idx * row_width]), local_covars_vcmaj_f_iter, row_copy_ct * row_width * sizeof(float));
    }
    local_covars_vcmaj_f_iter = &(local_covars_vcmaj_f_iter[row_ct * row_width]);
  } else {
    double* local_covars_vcmaj_d_iter = *local_covars_vcmaj_d_iterp;
    for (uint32_t row_idx = 1; row_idx < row_ct; row_idx *= 2) {
      uint32_t row_copy_ct = row_ct - row_idx;
      if (row_copy_ct > row_idx) {
        row_copy_ct = row_idx;
      }
      memcpy(&(local_covars_vcmaj_d_iter[row_idx * row_width]), local_covars_vcmaj_d_iter, row_copy_ct * row_width * sizeof(double));
    }
    *local_covars_vcmaj_d_iterp = &(local_covars_vcmaj_d_iter[row_ct * row_width]);
  }
  *variant_bidx_ptr += row_ct;
}


PglErr ReadRfmix2Block(const GlmCtx* common, const uint32_t* variant_bps, const uint32_t* local_sample_uidx_order, const float* prev_local_covar_row_f, const double* prev_local_covar_row_d, uint32_t variant_uidx_start, uint32_t variant_uidx_end, uint32_t cur_block_variant_ct, uint32_t local_sample_ct, uint32_t local_cat_ct, uint32_t local_chrom_col, uint32_t local_bp_col, uint32_t local_first_covar_col, TextStream* local_covar_txsp, const char** local_line_iterp, uint32_t* local_line_idx_ptr, uint32_t* local_prev_chr_code_ptr, uint32_t* local_chr_code_ptr, uint32_t* local_bp_ptr, uint32_t* local_skip_chr_ptr, float* local_covars_vcmaj_f_iter, double* local_covars_vcmaj_d_iter, uint32_t* local_sample_idx_order) {
  const ChrInfo* cip = common->cip;
  const uintptr_t* variant_include = common->variant_include;
  const uint32_t sample_ct = common->sample_ct;
  const uint32_t sample_ct_x = common->sample_ct_x;
  const uint32_t sample_ct_y = common->sample_ct_y;
  const uint32_t local_covar_ct = common->local_covar_ct;
  const GlmFlags flags = common->glm_flags;
  const uint32_t omit_last = (flags / kfGlmLocalOmitLast) & 1;
  const uint32_t local_haps = (flags / kfGlmLocalHaps) & 1;
  // There are several complications here:
  // 1. Until we've seen the beginning of the next line, we don't know the end
  //    of the (possibly empty) variant range the current set of local
  //    covariates applies to.
  // 2. The aforementioned variant range may go past the end of the current
  //    variant block.
  // 3. There may be variants which aren't contained in any local-covar
  //    interval (either before the first local-covar position on that
  //    chromosome, or on a chromosome absent from the local-covar file).
  // To address (1), we structure the main loop such that the first part
  // duplicates row contents and advances the matrix-filling iterator when
  // necessary, the second part speculatively parses local covariates to the
  // next matrix row, and the last part reads the position fields on the next
  // line.
  // To manage (2), we track local_prev_chr_code, local_chr_code, local_bp, and
  // prev_local_covar_row across function calls.  When local_prev_chr_code !=
  // local_chr_code on function restart, we know we need to replicate
  // prev_local_covar_row to the end of that chromosome, exiting early if the
  // end of the chromosome is past variant_uidx_end; otherwise, we need to
  // replicate it up to the first variant with position >= local_bp, possibly
  // exiting early.  This is no different from within-variant-block iteration,
  // so the main loop is written such that function reentry isn't
  // special-cased outside of needing to know prev_local_covar_row.
  // For (3), we zero-fill the affected local covariate rows.
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
  if (max_sample_ct < sample_ct_y) {
    max_sample_ct = sample_ct_y;
  }

  LocalCovarCoeffparseCtx coeffparse_ctx;
  coeffparse_ctx.sample_idx_order = local_sample_idx_order;
  coeffparse_ctx.max_sample_ct = max_sample_ct;
  // cur_sample_ct filled a bit later
  coeffparse_ctx.tokens_per_sample = (local_cat_ct? 1 : (local_covar_ct + omit_last)) << local_haps;
  coeffparse_ctx.local_covar_ct = local_covar_ct;
  coeffparse_ctx.omit_last = omit_last;
  coeffparse_ctx.local_haps = local_haps;
  coeffparse_ctx.local_cat_ct = local_cat_ct;
  coeffparse_ctx.local_cats_1based = (flags / kfGlmLocalCats1based) & 1;

  const char* local_line_iter = *local_line_iterp;
  uint32_t local_line_idx = *local_line_idx_ptr;
  uint32_t first_skip;
  uint32_t second_skip;
  uint32_t last_skip;
  if (local_chrom_col < local_bp_col) {
    first_skip = local_chrom_col - 1;
    second_skip = local_bp_col - local_chrom_col;
    last_skip = local_first_covar_col - local_bp_col;
  } else {
    first_skip = local_bp_col - 1;
    second_skip = local_chrom_col - local_bp_col;
    last_skip = local_first_covar_col - local_chrom_col;
  }
  const uintptr_t row_width = local_covar_ct * max_sample_ct;
  uint32_t local_prev_chr_code = *local_prev_chr_code_ptr;
  uint32_t local_chr_code = *local_chr_code_ptr;
  uint32_t local_bp = *local_bp_ptr;
  uint32_t local_skip_chr = *local_skip_chr_ptr;
  // might be inaccurate at file-start and EOF, but that's okay since these
  // variables are reinitialized before use in the first case, and irrelevant
  // in the second
  uint32_t is_x = (local_chr_code == x_code);
  uint32_t is_y = (local_chr_code == y_code);
  uint32_t local_xy = is_x + 2 * is_y;
  const uintptr_t* cur_sample_include;
  const uint32_t* cur_sample_include_cumulative_popcounts;
  if (is_y && common->sample_include_y) {
    cur_sample_include = common->sample_include_y;
    cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
    coeffparse_ctx.cur_sample_ct = sample_ct_y;
  } else if (is_x && common->sample_include_x) {
    cur_sample_include = common->sample_include_x;
    cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
    coeffparse_ctx.cur_sample_ct = sample_ct_x;
  } else {
    cur_sample_include = common->sample_include;
    cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
    coeffparse_ctx.cur_sample_ct = sample_ct;
  }
  // Necessary for cross-block duplication to work.
  if (prev_local_covar_row_f) {
    memcpy(local_covars_vcmaj_f_iter, prev_local_covar_row_f, row_width * sizeof(float));
  } else if (prev_local_covar_row_d) {
    memcpy(local_covars_vcmaj_d_iter, prev_local_covar_row_d, row_width * sizeof(double));
  }

  uint32_t variant_bidx = 0;
  uint32_t variant_uidx = AdvTo1Bit(variant_include, variant_uidx_start);
  uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
  uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
  uint32_t variant_uidx_chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
  if (variant_uidx_chr_end > variant_uidx_end) {
    variant_uidx_chr_end = variant_uidx_end;
  }
  uint32_t variant_bp = (chr_idx == local_prev_chr_code)? variant_bps[variant_uidx] : UINT32_MAX;
  while (1) {
    // Loop invariants:
    //   local_chr_code and local_bp are from the current line, if we're not at
    //     the beginning or end of the file.  In the latter two cases,
    //     local_chr_code is UINT32_MAX.
    //   local_prev_chr_code is from the previous relevant line, or UINT32_MAX
    //     at the beginning of the file.
    //   variant_uidx = next variant we need to fill local covariates for.
    //   chr_fo_idx and chr_idx correspond to variant_uidx.
    //   variant_uidx_chr_end is min(end of chr_fo_idx, variant_uidx_end).
    //   variant_bp is UINT32_MAX if we've already iterated through all
    //     variants on local_prev_chr_code (or we're at the beginning of the
    //     file); otherwise it also corresponds to variant_uidx.
    //   local_skip_chr is set iff either all variants on local_chr_code were
    //     filtered out, or local_chr_code == UINT32_MAX.
    //   If local_skip_chr isn't true, variant_uidx is not before the start of
    //     local_prev_chr_code, or after the end of local_chr_code.
    if (!local_skip_chr) {
      // Part 1.
      const uint32_t backfill_stop_bp = (local_chr_code == chr_idx)? local_bp : UINT32_MAX;
      if (backfill_stop_bp > variant_bp) {
        // At least one variant in [prev pos, current pos) (or [prev pos, end
        // of chromosome) if we're moving on to a new chromosome).  Count how
        // many there are, and fill that many lines of the matrix.
        const uint32_t next_variant_uidx = variant_uidx + ExpsearchU32(&(variant_bps[variant_uidx]), variant_uidx_chr_end - variant_uidx, backfill_stop_bp);
        const uint32_t row_ct = PopcountBitRange(variant_include, variant_uidx, next_variant_uidx);
        DuplicateLocalCovarRow(row_ct, row_width, &variant_bidx, &local_covars_vcmaj_f_iter, &local_covars_vcmaj_d_iter);
        if (variant_bidx == cur_block_variant_ct) {
          break;
        }
        variant_uidx = AdvTo1Bit(variant_include, next_variant_uidx);
        if (variant_uidx >= variant_uidx_chr_end) {
          // probable todo: helper function for this
          chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
          chr_idx = cip->chr_file_order[chr_fo_idx];
          variant_uidx_chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          if (variant_uidx_chr_end > variant_uidx_end) {
            variant_uidx_chr_end = variant_uidx_end;
          }
          variant_bp = UINT32_MAX;
        }
      }
      if (local_chr_code != local_prev_chr_code) {
        const uint32_t local_chr_fo_idx = cip->chr_idx_to_foidx[local_chr_code];
        if (local_chr_code != chr_idx) {
          // Some variants are on chromosomes entirely absent from the
          // local-covar file.  Zero-fill the covariate values for these
          // variants; NA results will be reported for them.
          const uint32_t local_chr_start_vidx = cip->chr_fo_vidx_start[local_chr_fo_idx];
          uint32_t row_ct;
          if (local_chr_start_vidx >= variant_uidx_chr_end) {
            row_ct = cur_block_variant_ct - variant_bidx;
          } else {
            // Verify that we aren't going backwards.  (This condition works
            // since we don't enter the block in the local_skip_chr case.)
            if (unlikely(local_chr_start_vidx < variant_uidx)) {
              logputs("\n");
              logerrputs("Error: --glm local-covar= file has a different chromosome order than the main\ndataset.\n");
              return kPglRetInconsistentInput;
            }
            row_ct = PopcountBitRange(variant_include, variant_uidx, local_chr_start_vidx);
          }
          ZeroLocalCovarRows(row_ct, row_width, &variant_bidx, &local_covars_vcmaj_f_iter, &local_covars_vcmaj_d_iter);
          if (variant_bidx == cur_block_variant_ct) {
            break;
          }
          variant_uidx = AdvTo1Bit(variant_include, local_chr_start_vidx);
          chr_fo_idx = local_chr_fo_idx;
          chr_idx = local_chr_code;
          variant_uidx_chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          if (variant_uidx_chr_end > variant_uidx_end) {
            variant_uidx_chr_end = variant_uidx_end;
          }
        }
        variant_bp = variant_bps[variant_uidx];
      }
      if (local_cat_ct) {
        if (local_covars_vcmaj_f_iter) {
          ZeroFArr(local_covar_ct * max_sample_ct, local_covars_vcmaj_f_iter);
        } else {
          ZeroDArr(local_covar_ct * max_sample_ct, local_covars_vcmaj_d_iter);
        }
      }
      // Part 2.
      PglErr reterr = LoadLocalCovarCoeffs(&coeffparse_ctx, local_line_iter, local_line_idx, local_covars_vcmaj_f_iter, local_covars_vcmaj_d_iter);
      if (unlikely(reterr)) {
        return reterr;
      }
    }
    // Part 3.

    // Not const for now, due to limitation of GetChrCodeCounted().
    char* line_start = TextGet(local_covar_txsp);
    if (!line_start) {
      PglErr reterr = TextStreamErrcode(local_covar_txsp);
      if (unlikely(reterr)) {
        return reterr;
      }
      // EOF.
      local_line_iter = nullptr;
      local_chr_code = UINT32_MAX;
      local_skip_chr = 1;
      if (local_prev_chr_code == chr_idx) {
        uint32_t row_ct = cur_block_variant_ct - variant_bidx;
        if (variant_uidx_chr_end < variant_uidx_end) {
          row_ct = PopcountBitRange(variant_include, variant_uidx, variant_uidx_chr_end);
        }
        if (row_ct) {
          DuplicateLocalCovarRow(row_ct, row_width, &variant_bidx, &local_covars_vcmaj_f_iter, &local_covars_vcmaj_d_iter);
          if (variant_bidx == cur_block_variant_ct) {
            break;
          }
        }
        variant_uidx = AdvTo1Bit(variant_include, variant_uidx_chr_end);
        chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
        chr_idx = cip->chr_file_order[chr_fo_idx];
        // no need to update variant_uidx_chr_end or variant_bp
      }
      const uintptr_t remaining_variant_ct = cur_block_variant_ct - variant_bidx;
      ZeroLocalCovarRows(remaining_variant_ct, row_width, &variant_bidx, &local_covars_vcmaj_f_iter, &local_covars_vcmaj_d_iter);
      break;
    }
    ++local_line_idx;
    char* tok1_start = NextTokenMult0(line_start, first_skip);
    if (unlikely(!tok1_start)) {
      logputs("\n");
      logerrprintf("Error: Line %u of --glm local-covar= file has fewer tokens than expected.\n", local_line_idx);
      return kPglRetMalformedInput;
    }
    char* tok1_end = CurTokenEnd(tok1_start);
    char* tok2_start = NextTokenMult(tok1_end, second_skip);
    if (unlikely(!tok2_start)) {
      logputs("\n");
      logerrprintf("Error: Line %u of --glm local-covar= file has fewer tokens than expected.\n", local_line_idx);
      return kPglRetMalformedInput;
    }
    char* tok2_end = CurTokenEnd(tok2_start);
    local_line_iter = NextTokenMult(tok2_end, last_skip);
    if (unlikely(!local_line_iter)) {
      logputs("\n");
      logerrprintf("Error: Line %u of --glm local-covar= file has fewer tokens than expected.\n", local_line_idx);
      return kPglRetMalformedInput;
    }
    char* chr_code_start = tok1_start;
    char* chr_code_end = tok1_end;
    char* bp_col_start = tok2_start;
    if (local_chrom_col > local_bp_col) {
      chr_code_start = tok2_start;
      chr_code_end = tok2_end;
      bp_col_start = tok1_start;
    }
    const uint32_t next_chr_code = GetChrCodeCounted(cip, chr_code_end - chr_code_start, chr_code_start);
    if (next_chr_code != local_chr_code) {
      local_skip_chr = IsI32Neg(next_chr_code) || (!IsSet(cip->chr_mask, next_chr_code));
      if (local_skip_chr) {
        continue;
      }
      // Still possible for all variants in this chromosome to have been
      // filtered out.
      const uint32_t new_chr_fo_idx = cip->chr_idx_to_foidx[next_chr_code];
      if (local_chr_code != UINT32_MAX) {
        // May as well sanity-check this.
        const uint32_t old_chr_fo_idx = cip->chr_idx_to_foidx[local_chr_code];
        if (unlikely(old_chr_fo_idx > new_chr_fo_idx)) {
          logputs("\n");
          logerrputs("Error: --glm local-covar= file has a different chromosome order than the main\ndataset.\n");
          return kPglRetInconsistentInput;
        }
      }
      const uint32_t chr_start_vidx = cip->chr_fo_vidx_start[new_chr_fo_idx];
      const uint32_t chr_end_vidx = cip->chr_fo_vidx_start[new_chr_fo_idx + 1];
      local_skip_chr = !PopcountBitRange(variant_include, chr_start_vidx, chr_end_vidx);
      if (local_skip_chr) {
        continue;
      }
      // May as well ensure that local_prev_chr_code is either UINT32_MAX,
      // or corresponds to a not-totally-excluded chromosome, so failures
      // of the sanity check above are more likely to be meaningful.
      local_prev_chr_code = local_chr_code;
      local_chr_code = next_chr_code;
      local_bp = UINT32_MAX;
      is_x = (local_chr_code == x_code);
      is_y = (local_chr_code == y_code);
      const uint32_t new_local_xy = is_x + 2 * is_y;
      if (new_local_xy != local_xy) {
        local_xy = new_local_xy;
        if (is_y && common->sample_include_y) {
          cur_sample_include = common->sample_include_y;
          cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
          coeffparse_ctx.cur_sample_ct = sample_ct_y;
        } else if (is_x && common->sample_include_x) {
          cur_sample_include = common->sample_include_x;
          cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
          coeffparse_ctx.cur_sample_ct = sample_ct_x;
        } else {
          cur_sample_include = common->sample_include;
          cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
          coeffparse_ctx.cur_sample_ct = sample_ct;
        }
        for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
          const uint32_t cur_uidx = local_sample_uidx_order[uii];
          uint32_t cur_idx = UINT32_MAX;
          if ((cur_uidx != UINT32_MAX) && IsSet(cur_sample_include, cur_uidx)) {
            cur_idx = RawToSubsettedPos(cur_sample_include, cur_sample_include_cumulative_popcounts, cur_uidx);
          }
          local_sample_idx_order[uii] = cur_idx;
        }
      }
    } else if (local_skip_chr) {
      continue;
    }
    uint32_t next_bp;
    if (unlikely(ScanPosintDefcap(bp_col_start, &next_bp))) {
      logputs("\n");
      logerrprintf("Error: Line %u of --glm local-covar= file has fewer tokens than expected.\n", local_line_idx);
      return kPglRetMalformedInput;
    }
    // could conditionally prohibit duplicate positions
    if (unlikely(S_CAST(int32_t, next_bp) < S_CAST(int32_t, local_bp))) {
      logputs("\n");
      logerrputs("Error: Positions in --glm local-covar= file are not sorted.\n");
      return kPglRetMalformedInput;
    }
    local_bp = next_bp;
  }
  *local_line_iterp = local_line_iter;
  *local_line_idx_ptr = local_line_idx;
  *local_prev_chr_code_ptr = local_prev_chr_code;
  *local_chr_code_ptr = local_chr_code;
  *local_bp_ptr = local_bp;
  *local_skip_chr_ptr = local_skip_chr;
  return kPglRetSuccess;
}

// only pass the parameters which aren't also needed by the compute threads,
// for now
// valid_variants and valid_alleles are a bit redundant, may want to remove the
// former later, but let's make that decision during/after permutation test
// implementation
PglErr GlmLogistic(const char* cur_pheno_name, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const char* outname, uint32_t raw_variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLogisticCtx* ctx, TextStream* local_covar_txsp, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, double* orig_permstat, uintptr_t* valid_allele_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  ThreadGroup tg;
  PreinitCstream(&css);
  PreinitThreads(&tg);
  {
    GlmCtx* common = ctx->common;
    const uintptr_t* variant_include = common->variant_include;
    const ChrInfo* cip = common->cip;
    const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
    const AlleleCode* omitted_alleles = common->omitted_alleles;

    const uint32_t sample_ct = common->sample_ct;
    const uint32_t sample_ct_x = common->sample_ct_x;
    const uint32_t sample_ct_y = common->sample_ct_y;
    const uint32_t covar_ct = common->covar_ct;
    const uintptr_t local_covar_ct = common->local_covar_ct;
    const uint32_t covar_ct_x = common->covar_ct_x;
    const uint32_t covar_ct_y = common->covar_ct_y;

    uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
    if (max_sample_ct < sample_ct_y) {
      max_sample_ct = sample_ct_y;
    }
    // obvious todo: wrap these in structs
    uint32_t* local_sample_idx_order = nullptr;
    uint32_t local_line_idx = 0;
    uint32_t local_xy = 0;  // 1 = chrX, 2 = chrY

    const char* local_line_iter = nullptr;
    uint32_t local_prev_chr_code = UINT32_MAX;
    uint32_t local_chr_code = UINT32_MAX;
    uint32_t local_bp = UINT32_MAX;
    uint32_t local_skip_chr = 1;
    if (local_covar_ct) {
      reterr = TextRewind(local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLogistic_ret_TSTREAM_FAIL;
      }
      local_line_idx = glm_info_ptr->local_header_line_ct;
      reterr = TextSkip(local_line_idx, local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLogistic_ret_TSTREAM_FAIL;
      }
      if (unlikely(bigstack_alloc_u32(local_sample_ct, &local_sample_idx_order))) {
        goto GlmLogistic_ret_NOMEM;
      }
      for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
        const uint32_t cur_uidx = local_sample_uidx_order[uii];
        uint32_t cur_idx = UINT32_MAX;
        if ((cur_uidx != UINT32_MAX) && IsSet(common->sample_include, cur_uidx)) {
          cur_idx = RawToSubsettedPos(common->sample_include, common->sample_include_cumulative_popcounts, cur_uidx);
        }
        local_sample_idx_order[uii] = cur_idx;
      }
    }

    const uint32_t variant_ct = common->variant_ct;

    const GlmFlags glm_flags = glm_info_ptr->flags;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    // forced-singlethreaded
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto GlmLogistic_ret_1;
    }
    const uint32_t report_neglog10p = (glm_flags / kfGlmLog10) & 1;
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t domdev_present_p1 = domdev_present + 1;

    const uint32_t constraint_ct = common->constraint_ct;
    const uint32_t constraint_ct_x = common->constraint_ct_x;
    const uint32_t constraint_ct_y = common->constraint_ct_y;

    const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
    uint32_t biallelic_predictor_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
    uint32_t biallelic_predictor_ct_x = 2 + domdev_present + covar_ct_x * (1 + add_interactions * domdev_present_p1);
    uint32_t biallelic_predictor_ct_y = 2 + domdev_present + covar_ct_y * (1 + add_interactions * domdev_present_p1);
    const uintptr_t* parameter_subset = common->parameter_subset;
    const uintptr_t* parameter_subset_x = common->parameter_subset_x;
    const uintptr_t* parameter_subset_y = common->parameter_subset_y;
    if (parameter_subset) {
      biallelic_predictor_ct = PopcountWords(parameter_subset, BitCtToWordCt(biallelic_predictor_ct));
      if (sample_ct_x) {
        biallelic_predictor_ct_x = PopcountWords(parameter_subset_x, BitCtToWordCt(biallelic_predictor_ct_x));
      } else {
        biallelic_predictor_ct_x = 0;
      }
      if (sample_ct_y) {
        biallelic_predictor_ct_y = PopcountWords(parameter_subset_y, BitCtToWordCt(biallelic_predictor_ct_x));
      } else {
        biallelic_predictor_ct_y = 0;
      }
    }
    uint32_t biallelic_reported_test_ct = GetBiallelicReportedTestCt(parameter_subset, glm_flags, covar_ct, common->tests_flag);
    uintptr_t max_reported_test_ct = biallelic_reported_test_ct;
    uint32_t biallelic_reported_test_ct_x = 0;
    if (sample_ct_x) {
      biallelic_reported_test_ct_x = GetBiallelicReportedTestCt(parameter_subset_x, glm_flags, covar_ct_x, common->tests_flag);
      if (biallelic_reported_test_ct_x > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_x;
      }
    }
    uint32_t biallelic_reported_test_ct_y = 0;
    if (sample_ct_y) {
      biallelic_reported_test_ct_y = GetBiallelicReportedTestCt(parameter_subset_y, glm_flags, covar_ct_y, common->tests_flag);
      if (biallelic_reported_test_ct_y > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_y;
      }
    }
    const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
    const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
    const GlmColFlags glm_cols = glm_info_ptr->cols;
    const uint32_t test_col = glm_cols & kfGlmColTest;
    if (unlikely((!test_col) && (max_reported_test_ct > 1))) {
      // this is okay in plain multiallelic case due to A1 column
      logerrputs("Error: --glm's 'test' column cannot be omitted when results for multiple\npredictors are reported.  (Did you forget 'hide-covar'?)\n");
      goto GlmLogistic_ret_INCONSISTENT_INPUT;
    }
    const uint32_t main_mutated = ((glm_flags & (kfGlmDominant | kfGlmRecessive | kfGlmHethom)) != kfGlm0);
    // if 'fused', one row per variant
    // otherwise, one row per tested allele
    const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!main_mutated) && (!common->tests_flag) && (!add_interactions);
    if (beta_se_multiallelic_fused || (!hide_covar)) {
      max_reported_test_ct += max_extra_allele_ct;
    }
    common->max_reported_test_ct = max_reported_test_ct;

    const uint32_t is_sometimes_firth = !(glm_flags & kfGlmNoFirth);
    const uint32_t is_always_firth = (glm_flags / kfGlmFirth) & 1;

    uint32_t x_code = UINT32_MAXM1;
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    if (sample_ct_x) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetX, &x_code, &x_start, &x_end);
    }
    uint32_t y_code = UINT32_MAXM1;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    if (sample_ct_y) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetY, &y_code, &y_start, &y_end);
    }
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t chr_col = glm_cols & kfGlmColChrom;

    // includes trailing tab
    char* chr_buf = nullptr;
    if (chr_col) {
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto GlmLogistic_ret_NOMEM;
      }
    }

    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > variant_ct) {
      calc_thread_ct = variant_ct;
    }

    const uint32_t main_omitted = parameter_subset && (!IsSet(parameter_subset, 1));
    const uint32_t xmain_ct = main_mutated + main_omitted;
    const uint32_t gcount_cc_col = glm_cols & kfGlmColGcountcc;
    // workflow is similar to --make-bed
    uintptr_t workspace_alloc = GetLogisticWorkspaceSize(sample_ct, biallelic_predictor_ct, max_extra_allele_ct, constraint_ct, xmain_ct, gcount_cc_col, is_sometimes_firth);
    if (sample_ct_x) {
      const uintptr_t workspace_alloc_x = GetLogisticWorkspaceSize(sample_ct_x, biallelic_predictor_ct_x, max_extra_allele_ct, constraint_ct_x, xmain_ct, gcount_cc_col, is_sometimes_firth);
      if (workspace_alloc_x > workspace_alloc) {
        workspace_alloc = workspace_alloc_x;
      }
    }
    if (sample_ct_y) {
      const uintptr_t workspace_alloc_y = GetLogisticWorkspaceSize(sample_ct_y, biallelic_predictor_ct_y, max_extra_allele_ct, constraint_ct_y, xmain_ct, gcount_cc_col, is_sometimes_firth);
      if (workspace_alloc_y > workspace_alloc) {
        workspace_alloc = workspace_alloc_y;
      }
    }
    // +1 is for top-level common->workspace_bufs
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;

    uintptr_t per_variant_xalloc_byte_ct = max_sample_ct * local_covar_ct * sizeof(float);
    uintptr_t per_alt_allele_xalloc_byte_ct = sizeof(LogisticAuxResult);
    if (beta_se_multiallelic_fused) {
      per_variant_xalloc_byte_ct += 2 * max_reported_test_ct * sizeof(double);
    } else {
      per_alt_allele_xalloc_byte_ct += 2 * max_reported_test_ct * sizeof(double);
    }
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    common->thread_mhc = nullptr;
    common->dosage_presents = nullptr;
    common->dosage_mains = nullptr;
    uint32_t read_block_size;
    uintptr_t max_alt_allele_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, max_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, per_alt_allele_xalloc_byte_ct, pgfip, &calc_thread_ct, &common->genovecs, max_extra_allele_ct? (&common->thread_mhc) : nullptr, nullptr, nullptr, dosage_is_present? (&common->dosage_presents) : nullptr, dosage_is_present? (&common->dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, &max_alt_allele_block_size, main_loadbufs, &common->pgr_ptrs, &common->read_variant_uidx_starts))) {
      goto GlmLogistic_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto GlmLogistic_ret_NOMEM;
    }
    LogisticAuxResult* logistic_block_aux_bufs[2];
    double* block_beta_se_bufs[2];

    for (uint32_t uii = 0; uii != 2; ++uii) {
      if (unlikely(BIGSTACK_ALLOC_X(LogisticAuxResult, max_alt_allele_block_size, &(logistic_block_aux_bufs[uii])))) {
        goto GlmLogistic_ret_NOMEM;
      }
      if (beta_se_multiallelic_fused) {
        if (unlikely(bigstack_alloc_d(read_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii])))) {
          goto GlmLogistic_ret_NOMEM;
        }
      } else {
        if (unlikely(bigstack_alloc_d(max_alt_allele_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii])))) {
          goto GlmLogistic_ret_NOMEM;
        }
      }
      if (local_covar_ct) {
        // bugfix (18 May 2018): don't want sizeof(float) here
        if (unlikely(bigstack_alloc_f(read_block_size * max_sample_ct * local_covar_ct, &(ctx->local_covars_vcmaj_f[uii])))) {
          goto GlmLogistic_ret_NOMEM;
        }
      } else {
        ctx->local_covars_vcmaj_f[uii] = nullptr;
      }
    }

    if (max_sample_ct > 2000000) {
      // may eventually want a large-matrix double-precision fallback, but that
      // can probably wait till 2020 or later
      logerrputs("Warning: --glm logistic regression is unreliable on more than ~2 million\nsamples, since it uses single-precision arithmetic.\n");
    }
    common->workspace_bufs = S_CAST(unsigned char**, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t)));
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      common->workspace_bufs[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(workspace_alloc));
    }
    common->err_info = (~0LLU) << 32;
    SetThreadFuncAndData(GlmLogisticThread, ctx, &tg);

    const uint32_t ref_col = glm_cols & kfGlmColRef;
    const uint32_t alt1_col = glm_cols & kfGlmColAlt1;
    const uint32_t alt_col = glm_cols & kfGlmColAlt;
    const uint32_t ax_col = glm_cols & kfGlmColAx;
    const uint32_t a1_ct_col = glm_cols & kfGlmColA1count;
    const uint32_t tot_allele_col = glm_cols & kfGlmColTotallele;
    const uint32_t a1_ct_cc_col = glm_cols & kfGlmColA1countcc;
    const uint32_t tot_allele_cc_col = glm_cols & kfGlmColTotallelecc;
    const uint32_t a1_freq_col = glm_cols & kfGlmColA1freq;
    const uint32_t a1_freq_cc_col = glm_cols & kfGlmColA1freqcc;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t firth_yn_col = (glm_cols & kfGlmColFirthYn) && is_sometimes_firth && (!is_always_firth);
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t orbeta_col = glm_cols & (kfGlmColBeta | kfGlmColOrbeta);
    const uint32_t report_beta_instead_of_odds_ratio = glm_cols & kfGlmColBeta;
    const uint32_t se_col = glm_cols & kfGlmColSe;
    const uint32_t ci_col = (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t z_col = glm_cols & kfGlmColTz;
    const uint32_t p_col = glm_cols & kfGlmColP;
    const uint32_t err_col = glm_cols & kfGlmColErr;
    *cswritep++ = '#';
    if (chr_col) {
      cswritep = strcpya_k(cswritep, "CHROM\t");
    }
    if (variant_bps) {
      cswritep = strcpya_k(cswritep, "POS\t");
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
    if (a1_ct_col) {
      cswritep = strcpya_k(cswritep, "\tA1_CT");
    }
    if (tot_allele_col) {
      cswritep = strcpya_k(cswritep, "\tALLELE_CT");
    }
    if (a1_ct_cc_col) {
      cswritep = strcpya_k(cswritep, "\tA1_CASE_CT\tA1_CTRL_CT");
    }
    if (tot_allele_cc_col) {
      cswritep = strcpya_k(cswritep, "\tCASE_ALLELE_CT\tCTRL_ALLELE_CT");
    }
    if (gcount_cc_col) {
      cswritep = strcpya_k(cswritep, "\tCASE_NON_A1_CT\tCASE_HET_A1_CT\tCASE_HOM_A1_CT\tCTRL_NON_A1_CT\tCTRL_HET_A1_CT\tCTRL_HOM_A1_CT");
    }
    if (a1_freq_col) {
      cswritep = strcpya_k(cswritep, "\tA1_FREQ");
    }
    if (a1_freq_cc_col) {
      cswritep = strcpya_k(cswritep, "\tA1_CASE_FREQ\tA1_CTRL_FREQ");
    }
    if (mach_r2_col) {
      cswritep = strcpya_k(cswritep, "\tMACH_R2");
    }
    if (firth_yn_col) {
      cswritep = strcpya_k(cswritep, "\tFIRTH?");
    }
    if (test_col) {
      cswritep = strcpya_k(cswritep, "\tTEST");
    }
    if (nobs_col) {
      cswritep = strcpya_k(cswritep, "\tOBS_CT");
    }
    if (orbeta_col) {
      if (report_beta_instead_of_odds_ratio) {
        cswritep = strcpya_k(cswritep, "\tBETA");
      } else {
        cswritep = strcpya_k(cswritep, "\tOR");
      }
    }
    if (se_col) {
      if (report_beta_instead_of_odds_ratio) {
        cswritep = strcpya_k(cswritep, "\tSE");
      } else {
        cswritep = strcpya_k(cswritep, "\tLOG(OR)_SE");
      }
    }
    double ci_zt = 0.0;
    if (ci_col) {
      cswritep = strcpya_k(cswritep, "\tL");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      cswritep = strcpya_k(cswritep, "\tU");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      ci_zt = QuantileToZscore((ci_size + 1.0) * 0.5);
    }
    if (z_col) {
      if (!constraint_ct) {
        cswritep = strcpya_k(cswritep, "\tZ_STAT");
      } else {
        // F-statistic for joint tests.
        cswritep = strcpya_k(cswritep, "\tZ_OR_F_STAT");
      }
    }
    if (p_col) {
      if (report_neglog10p) {
        cswritep = strcpya_k(cswritep, "\tLOG10_P");
      } else {
        cswritep = strcpya_k(cswritep, "\tP");
      }
    }
    if (err_col) {
      cswritep = strcpya_k(cswritep, "\tERRCODE");
    }
    AppendBinaryEoln(&cswritep);

    // Main workflow:
    // 1. Set n=0, load/skip block 0
    //
    // 2. Spawn threads processing block n
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Load/skip block n unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8, Write results for last block
    uintptr_t write_variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t suppress_mach_r2 = 0;

    uint32_t cur_biallelic_reported_test_ct = 0;
    uint32_t primary_reported_test_idx = include_intercept;
    uint32_t cur_constraint_ct = 0;

    const char* const* cur_test_names = nullptr;
    uint32_t prev_block_variant_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t allele_ct = 2;
    uint32_t omitted_allele_idx = 0;
    uintptr_t valid_allele_ct = 0;
    logprintfww5("--glm %s regression on phenotype '%s': ", is_always_firth? "Firth" : (is_sometimes_firth? "logistic-Firth hybrid" : "logistic"), cur_pheno_name);
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_variant_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto GlmLogistic_ret_PGR_FAIL;
      }
      if (local_covar_ct && cur_block_variant_ct) {
        const uint32_t uidx_start = read_block_idx * read_block_size;
        const uint32_t uidx_end = MINV(raw_variant_ct, uidx_start + read_block_size);
        if (local_variant_include) {
          reterr = ReadLocalCovarBlock(common, local_sample_uidx_order, local_variant_include, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, local_covar_txsp, &local_line_idx, &local_xy, ctx->local_covars_vcmaj_f[parity], nullptr, local_sample_idx_order);
        } else {
          float* prev_local_covar_row_f = nullptr;
          if (variant_idx) {
            prev_local_covar_row_f = &(ctx->local_covars_vcmaj_f[1 - parity][S_CAST(uintptr_t, read_block_size - 1) * max_sample_ct * local_covar_ct]);
          }
          reterr = ReadRfmix2Block(common, variant_bps, local_sample_uidx_order, prev_local_covar_row_f, nullptr, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, glm_info_ptr->local_chrom_col, glm_info_ptr->local_bp_col, glm_info_ptr->local_first_covar_col, local_covar_txsp, &local_line_iter, &local_line_idx, &local_prev_chr_code, &local_chr_code, &local_bp, &local_skip_chr, ctx->local_covars_vcmaj_f[parity], nullptr, local_sample_idx_order);
        }
        if (unlikely(reterr)) {
          goto GlmLogistic_ret_1;
        }
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, common->err_info);
        if (unlikely(reterr)) {
          goto GlmLogistic_ret_PGR_FAIL;
        }
      }
      if (!IsLastBlock(&tg)) {
        common->cur_block_variant_ct = cur_block_variant_ct;
        const uint32_t uidx_start = read_block_idx * read_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_variant_ct, calc_thread_ct, uidx_start, common->read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, common->pgr_ptrs);
        ctx->block_aux = logistic_block_aux_bufs[parity];
        common->block_beta_se = block_beta_se_bufs[parity];
        if (variant_idx + cur_block_variant_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto GlmLogistic_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const double* beta_se_iter = block_beta_se_bufs[parity];
        const LogisticAuxResult* cur_block_aux = logistic_block_aux_bufs[parity];
        uintptr_t allele_bidx = 0;
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            if ((chr_idx == x_code) && sample_ct_x) {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct_x;
              cur_constraint_ct = constraint_ct_x;
              cur_test_names = test_names_x;
            } else if ((chr_idx == y_code) && sample_ct_y) {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct_y;
              cur_constraint_ct = constraint_ct_y;
              cur_test_names = test_names_y;
            } else {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct;
              cur_constraint_ct = constraint_ct;
              cur_test_names = test_names;
            }
            suppress_mach_r2 = (chr_idx == x_code) || (chr_idx == mt_code);
            if (cur_constraint_ct) {
              // bugfix (17 May 2018): this was using reported_test_ct instead
              // of cur_reported_test_ct.
              primary_reported_test_idx = cur_biallelic_reported_test_ct - 1;
            }
            if (chr_col) {
              char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
              *chr_name_end = '\t';
              chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
            }
          }
          uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offsets[write_variant_uidx];
          }
          const uint32_t allele_ct_m1 = allele_ct - 1;
          const uint32_t extra_allele_ct = allele_ct - 2;
          if (omitted_alleles) {
            omitted_allele_idx = omitted_alleles[write_variant_uidx];
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          uint32_t variant_is_valid = 0;
          uint32_t a1_allele_idx = 0;
          for (uint32_t nonomitted_allele_idx = 0; nonomitted_allele_idx != allele_ct_m1; ++nonomitted_allele_idx, ++a1_allele_idx) {
            if (beta_se_multiallelic_fused) {
              if (!nonomitted_allele_idx) {
                primary_reported_test_idx = include_intercept;
              } else {
                primary_reported_test_idx = cur_biallelic_reported_test_ct + nonomitted_allele_idx - 1;
              }
            }
            if (nonomitted_allele_idx == omitted_allele_idx) {
              ++a1_allele_idx;
            }
            const double primary_beta = beta_se_iter[primary_reported_test_idx * 2];
            const double primary_se = beta_se_iter[primary_reported_test_idx * 2 + 1];
            const uint32_t allele_is_valid = (primary_se != -9.0);
            variant_is_valid |= allele_is_valid;
            {
              const LogisticAuxResult* auxp = &(cur_block_aux[allele_bidx]);
              if (ln_pfilter <= 0.0) {
                if (!allele_is_valid) {
                  goto GlmLogistic_allele_iterate;
                }
                double permstat;
                double primary_ln_pval;
                if (!cur_constraint_ct) {
                  permstat = fabs(primary_beta / primary_se);
                  // could precompute a tstat threshold instead
                  primary_ln_pval = ZscoreToLnP(permstat);
                } else {
                  // cur_constraint_ct may be different on chrX/chrY than it is
                  // on autosomes, so just have permstat be -log(pval) to be
                  // safe
                  primary_ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                  permstat = -primary_ln_pval;
                }
                if (primary_ln_pval > ln_pfilter) {
                  if (orig_ln_pvals) {
                    orig_ln_pvals[valid_allele_ct] = primary_ln_pval;
                  }
                  if (orig_permstat) {
                    orig_permstat[valid_allele_ct] = permstat;
                  }
                  goto GlmLogistic_allele_iterate;
                }
              }
              uint32_t inner_reported_test_ct = cur_biallelic_reported_test_ct;
              if (extra_allele_ct) {
                if (beta_se_multiallelic_fused) {
                  // in fused case, we're only performing a single multiple
                  // regression, so list all additive results together,
                  // possibly with intercept before.
                  if (!nonomitted_allele_idx) {
                    inner_reported_test_ct = 1 + include_intercept;
                  } else if (nonomitted_allele_idx == extra_allele_ct) {
                    inner_reported_test_ct -= include_intercept;
                  } else {
                    inner_reported_test_ct = 1;
                  }
                } else if (!hide_covar) {
                  inner_reported_test_ct += extra_allele_ct;
                }
              }
              // possible todo: make number-to-string operations, strlen(),
              // etc. happen only once per variant.
              for (uint32_t allele_test_idx = 0; allele_test_idx != inner_reported_test_ct; ++allele_test_idx) {
                uint32_t test_idx = allele_test_idx;
                if (beta_se_multiallelic_fused && nonomitted_allele_idx) {
                  if (!allele_test_idx) {
                    test_idx = primary_reported_test_idx;
                  } else {
                    // bugfix (26 Jun 2019): only correct to add 1 here in
                    // include_intercept case
                    test_idx += include_intercept;
                  }
                }
                if (chr_col) {
                  cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
                }
                if (variant_bps) {
                  cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
                }
                cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
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
                  for (uint32_t tmp_allele_idx = 1; tmp_allele_idx != allele_ct; ++tmp_allele_idx) {
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLogistic_ret_WRITE_FAIL;
                    }
                    cswritep = strcpyax(cswritep, cur_alleles[tmp_allele_idx], ',');
                  }
                  --cswritep;
                }
                *cswritep++ = '\t';
                const uint32_t multi_a1 = extra_allele_ct && beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx);
                if (multi_a1) {
                  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    if (allele_idx == omitted_allele_idx) {
                      continue;
                    }
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLogistic_ret_WRITE_FAIL;
                    }
                    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                  }
                  --cswritep;
                } else {
                  cswritep = strcpya(cswritep, cur_alleles[a1_allele_idx]);
                }
                if (ax_col) {
                  *cswritep++ = '\t';
                  if (beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx)) {
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLogistic_ret_WRITE_FAIL;
                    }
                    cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
                  } else {
                    for (uint32_t tmp_allele_idx = 0; tmp_allele_idx != allele_ct; ++tmp_allele_idx) {
                      if (tmp_allele_idx == a1_allele_idx) {
                        continue;
                      }
                      if (unlikely(Cswrite(&css, &cswritep))) {
                        goto GlmLogistic_ret_WRITE_FAIL;
                      }
                      cswritep = strcpyax(cswritep, cur_alleles[tmp_allele_idx], ',');
                    }
                    --cswritep;
                  }
                }
                if (a1_ct_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_dosage, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (tot_allele_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa(auxp->allele_obs_ct, cswritep);
                }
                if (a1_ct_cc_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_case_dosage, cswritep);
                    *cswritep++ = '\t';
                    cswritep = dtoa_g(auxp->a1_dosage - auxp->a1_case_dosage, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA\tNA");
                  }
                }
                if (tot_allele_cc_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa_x(auxp->case_allele_obs_ct, '\t', cswritep);
                  cswritep = u32toa(auxp->allele_obs_ct - auxp->case_allele_obs_ct, cswritep);
                }
                if (gcount_cc_col) {
                  if (!multi_a1) {
                    STD_ARRAY_KREF(uint32_t, 6) cur_geno_hardcall_cts = auxp->geno_hardcall_cts;
                    for (uint32_t uii = 0; uii != 6; ++uii) {
                      *cswritep++ = '\t';
                      cswritep = u32toa(cur_geno_hardcall_cts[uii], cswritep);
                    }
                  } else {
                    cswritep = strcpya_k(cswritep, "\tNA\tNA\tNA\tNA\tNA\tNA");
                  }
                }
                if (a1_freq_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_dosage / S_CAST(double, auxp->allele_obs_ct), cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (a1_freq_cc_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_case_dosage / S_CAST(double, auxp->case_allele_obs_ct), cswritep);
                    *cswritep++ = '\t';
                    cswritep = dtoa_g((auxp->a1_dosage - auxp->a1_case_dosage) / S_CAST(double, auxp->allele_obs_ct - auxp->case_allele_obs_ct), cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA\tNA");
                  }
                }
                if (mach_r2_col) {
                  *cswritep++ = '\t';
                  if (!suppress_mach_r2) {
                    cswritep = dtoa_g(auxp->mach_r2, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (firth_yn_col) {
                  *cswritep++ = '\t';
                  // 'Y' - 'N' = 11
                  *cswritep++ = 'N' + 11 * auxp->firth_fallback;
                }
                if (test_col) {
                  *cswritep++ = '\t';
                  if (test_idx < cur_biallelic_reported_test_ct) {
                    cswritep = strcpya(cswritep, cur_test_names[test_idx]);
                  } else {
                    // always use basic dosage for untested alleles
                    cswritep = strcpya_k(cswritep, "ADD");
                    if (!beta_se_multiallelic_fused) {
                      // extra alt allele covariate.
                      uint32_t test_xallele_idx = test_idx - cur_biallelic_reported_test_ct;
                      // now we have the 0-based relative position in a list
                      // with the omitted_allele_idx and a1_allele_idx removed.
                      // correct this to the absolute index.  (there may be a
                      // cleaner way to do this with nonomitted_allele_idx?)
                      if (omitted_allele_idx < a1_allele_idx) {
                        test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                      }
                      test_xallele_idx = test_xallele_idx + (test_xallele_idx >= a1_allele_idx);
                      if (a1_allele_idx < omitted_allele_idx) {
                        test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                      }
                      if (!test_xallele_idx) {
                        cswritep = strcpya_k(cswritep, "_REF");
                      } else {
                        cswritep = strcpya_k(cswritep, "_ALT");
                        cswritep = u32toa(test_xallele_idx, cswritep);
                      }
                    }
                  }
                }
                if (nobs_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa(auxp->sample_obs_ct, cswritep);
                }
                double ln_pval = kLnPvalError;
                double permstat = 0.0;
                uint32_t test_is_valid;
                if ((!cur_constraint_ct) || (test_idx != primary_reported_test_idx)) {
                  double beta = beta_se_iter[2 * test_idx];
                  double se = beta_se_iter[2 * test_idx + 1];
                  test_is_valid = (se != -9.0);
                  if (test_is_valid) {
                    permstat = beta / se;
                    ln_pval = ZscoreToLnP(permstat);
                  }
                  if (orbeta_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      if (report_beta_instead_of_odds_ratio) {
                        cswritep = dtoa_g(beta, cswritep);
                      } else {
                        cswritep = lntoa_g(beta, cswritep);
                      }
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  if (se_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(se, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  if (ci_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      const double ci_halfwidth = ci_zt * se;
                      if (report_beta_instead_of_odds_ratio) {
                        cswritep = dtoa_g(beta - ci_halfwidth, cswritep);
                        *cswritep++ = '\t';
                        cswritep = dtoa_g(beta + ci_halfwidth, cswritep);
                      } else {
                        cswritep = lntoa_g(beta - ci_halfwidth, cswritep);
                        *cswritep++ = '\t';
                        cswritep = lntoa_g(beta + ci_halfwidth, cswritep);
                      }
                    } else {
                      cswritep = strcpya_k(cswritep, "NA\tNA");
                    }
                  }
                  if (z_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(permstat, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                } else {
                  // joint test: use F-test instead of Wald test
                  test_is_valid = allele_is_valid;
                  if (orbeta_col) {
                    cswritep = strcpya_k(cswritep, "\tNA");
                  }
                  if (se_col) {
                    cswritep = strcpya_k(cswritep, "\tNA");
                  }
                  if (ci_col) {
                    cswritep = strcpya_k(cswritep, "\tNA\tNA");
                  }
                  if (z_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(primary_se / u31tod(cur_constraint_ct), cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  // could avoid recomputing
                  if (test_is_valid) {
                    ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                    permstat = -ln_pval;
                  }
                }
                if (p_col) {
                  *cswritep++ = '\t';
                  if (test_is_valid) {
                    if (report_neglog10p) {
                      double reported_val = (-kRecipLn10) * ln_pval;
                      cswritep = dtoa_g(reported_val, cswritep);
                    } else {
                      double reported_ln = MAXV(ln_pval, output_min_ln);
                      cswritep = lntoa_g(reported_ln, cswritep);
                    }
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (err_col) {
                  *cswritep++ = '\t';
                  if (test_is_valid) {
                    if (!auxp->is_unfinished) {
                      *cswritep++ = '.';
                    } else {
                      cswritep = strcpya_k(cswritep, "UNFINISHED");
                    }
                  } else {
                    uint64_t glm_errcode;
                    memcpy(&glm_errcode, &(beta_se_iter[2 * test_idx]), 8);
                    cswritep = AppendGlmErrstr(glm_errcode, cswritep);
                  }
                }
                AppendBinaryEoln(&cswritep);
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto GlmLogistic_ret_WRITE_FAIL;
                }
                if ((test_idx == primary_reported_test_idx) && allele_is_valid) {
                  if (orig_ln_pvals) {
                    orig_ln_pvals[valid_allele_ct] = ln_pval;
                  }
                  if (orig_permstat) {
                    orig_permstat[valid_allele_ct] = permstat;
                  }
                }
              }
            }
          GlmLogistic_allele_iterate:
            ++allele_bidx;
            valid_allele_ct += allele_is_valid;
            if (valid_alleles && allele_is_valid) {
              SetBit(allele_idx_offset_base + a1_allele_idx, valid_alleles);
            }
            if (!beta_se_multiallelic_fused) {
              beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
            }
          }
          if (beta_se_multiallelic_fused) {
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
          if ((!variant_is_valid) && valid_alleles) {
            ClearBit(write_variant_uidx, valid_variants);
          }
        }
      }
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
      prev_block_variant_ct = cur_block_variant_ct;
      variant_idx += cur_block_variant_ct;
      // crucially, this is independent of the PgenReader block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto GlmLogistic_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    logprintf("Results written to %s .\n", outname);
    *valid_allele_ct_ptr = valid_allele_ct;
  }
  while (0) {
  GlmLogistic_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GlmLogistic_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
    break;
  GlmLogistic_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  GlmLogistic_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GlmLogistic_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  GlmLogistic_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 GlmLogistic_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

static const double kSmallDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, 3.0);

static const double kSmallInvDoublePairs[32] ALIGNV16 = PAIR_TABLE16(2.0, 1.0, 0.0, 3.0);

static const double kSmallInvDoubles[4] = {2.0, 1.0, 0.0, 3.0};

uint32_t GenoarrToDoublesRemoveMissing(const uintptr_t* genoarr, const double* __restrict table, uint32_t sample_ct, double* __restrict dst) {
  assert(sample_ct);
  const uint32_t sample_ctl2m1 = (sample_ct - 1) / kBitsPerWordD2;
  uint32_t subgroup_len = kBitsPerWordD2;
  double* dst_iter = dst;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2m1) {
      if (widx > sample_ctl2m1) {
        return dst_iter - dst;
      }
      subgroup_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genoarr[widx];
    for (uint32_t uii = 0; uii != subgroup_len; ++uii) {
      const uintptr_t cur_geno = geno_word & 3;
      if (cur_geno < 3) {
        // *dst_iter++ = u31tod(cur_geno);
        *dst_iter++ = table[cur_geno];
      }
      geno_word >>= 2;
    }
  }
}

uintptr_t GetLinearWorkspaceSize(uint32_t sample_ct, uint32_t biallelic_predictor_ct, uint32_t max_extra_allele_ct, uint32_t constraint_ct, uint32_t xmain_ct) {
  // sample_ct * max_predictor_ct < 2^31, and max_predictor_ct < sqrt(2^31), so
  // no overflows

  // sample_nm, tmp_nm = sample_ctl words
  uintptr_t workspace_size = 2 * RoundUpPow2(BitCtToWordCt(sample_ct) * sizeof(intptr_t), kCacheline);

  // nm_pheno_buf = sample_ct doubles
  workspace_size += RoundUpPow2(sample_ct * sizeof(double), kCacheline);

  const uint32_t max_predictor_ct = biallelic_predictor_ct + max_extra_allele_ct;
  // predictors_pmaj = (max_predictor_ct + main_mutated + main_omitted) * sample_ct
  // doubles
  workspace_size += RoundUpPow2((max_predictor_ct + xmain_ct) * sample_ct * sizeof(double), kCacheline);

  // xtx_inv = max_predictor_ct * max_predictor_ct doubles
  workspace_size += RoundUpPow2(max_predictor_ct * max_predictor_ct * sizeof(double), kCacheline);

  // dbl_2d_buf = max_predictor_ct * max(max_predictor_ct, 7) doubles
  workspace_size += RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);

  // inverse_corr_buf = (max_predictor_ct - 1) * max(max_predictor_ct - 1, 4) doubles
  workspace_size += RoundUpPow2((max_predictor_ct - 1) * MAXV((max_predictor_ct - 1), 4) * sizeof(double), kCacheline);

  // semicomputed_biallelic_corr_matrix = (max_predictor_ct - 1)^2 doubles
  workspace_size += RoundUpPow2((biallelic_predictor_ct - 1) * (biallelic_predictor_ct - 1) * sizeof(double), kCacheline);

  // semicomputed_biallelic_inv_corr_sqrts = biallelic_predictor_ct doubles
  workspace_size += RoundUpPow2(biallelic_predictor_ct * sizeof(double), kCacheline);

  // fitted_coefs, xt_y = max_predictor_ct doubles
  workspace_size += 2 * RoundUpPow2(max_predictor_ct * sizeof(double), kCacheline);

  // inv_1d_buf
  workspace_size += RoundUpPow2(MAXV(max_predictor_ct, constraint_ct) * kMatrixInvertBuf1CheckedAlloc, kCacheline);

  // machr2_dosage_sums, machr2_dosage_ssqs
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(uint64_t) * 2, kCacheline);

  if (constraint_ct) {
    // tmphxs_buf, h_transpose_buf = constraint_ct * max_predictor_ct doubles
    workspace_size += 2 * RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(double), kCacheline);

    // inner_buf = constraint_ct * constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * constraint_ct * sizeof(double), kCacheline);

    // cur_constraints_con_major = constraint_ct * max_predictor_ct doubles
    workspace_size += RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(double), kCacheline);
  }
  return workspace_size;
}

typedef struct GlmLinearCtxStruct {
  GlmCtx* common;

  double* pheno_d;
  double* pheno_x_d;
  double* pheno_y_d;
  const double* covars_cmaj_d;
  double* covars_cmaj_x_d;
  double* covars_cmaj_y_d;
  double* local_covars_vcmaj_d[2];
  LinearAuxResult* block_aux;

  uint32_t subbatch_size;
} GlmLinearCtx;

// possible todo: delete this, and GlmLinear(), if GlmLinearBatchThread is good
// enough at the same job.
THREAD_FUNC_DECL GlmLinearThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GlmLinearCtx* ctx = S_CAST(GlmLinearCtx*, arg->sharedp->context);
  GlmCtx* common = ctx->common;

  PgenReader* pgrp = common->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = common->genovecs[tidx];
  pgv.dosage_present = nullptr;
  pgv.dosage_main = nullptr;
  if (common->dosage_presents) {
    pgv.dosage_present = common->dosage_presents[tidx];
    pgv.dosage_main = common->dosage_mains[tidx];
  }
  unsigned char* workspace_buf = common->workspace_bufs[tidx];
  const uintptr_t* variant_include = common->variant_include;
  const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
  const AlleleCode* omitted_alleles = common->omitted_alleles;
  const uintptr_t* sex_male_collapsed = common->sex_male_collapsed;
  const ChrInfo* cip = common->cip;
  const uint32_t* subset_chr_fo_vidx_start = common->subset_chr_fo_vidx_start;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const GlmFlags glm_flags = common->glm_flags;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t joint_genotypic = (glm_flags / kfGlmGenotypic) & 1;
  const uint32_t joint_hethom = (glm_flags / kfGlmHethom) & 1;
  const uint32_t domdev_present = joint_genotypic || joint_hethom;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t reported_pred_uidx_start = 1 - include_intercept;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uintptr_t max_reported_test_ct = common->max_reported_test_ct;
  const uintptr_t local_covar_ct = common->local_covar_ct;
  const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
  const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!model_dominant) && (!model_recessive) && (!common->tests_flag) && (!add_interactions);
  uintptr_t max_sample_ct = MAXV(common->sample_ct, common->sample_ct_x);
  if (max_sample_ct < common->sample_ct_y) {
    max_sample_ct = common->sample_ct_y;
  }
  pgv.patch_01_set = nullptr;
  pgv.patch_01_vals = nullptr;
  pgv.patch_10_set = nullptr;
  pgv.patch_10_vals = nullptr;
  if (common->thread_mhc) {
    const uint32_t max_sample_ctl = BitCtToWordCt(max_sample_ct);
    pgv.patch_01_set = common->thread_mhc[tidx];
    pgv.patch_01_vals = R_CAST(AlleleCode*, &(pgv.patch_01_set[max_sample_ctl]));
    AlleleCode* patch_01_vals_end = &(pgv.patch_01_vals[max_sample_ct]);
    VecAlignUp(&patch_01_vals_end);
    pgv.patch_10_set = R_CAST(uintptr_t*, patch_01_vals_end);
    pgv.patch_10_vals = R_CAST(AlleleCode*, &(pgv.patch_10_set[max_sample_ctl]));
  }
  pgv.patch_01_ct = 0;
  pgv.patch_10_ct = 0;
  pgv.multidosage_sample_ct = 0;
  uint32_t variant_idx_offset = 0;
  uint32_t allele_ct = 2;
  uint32_t omitted_allele_idx = 0;
  uint32_t extra_regression_ct = 0;
  double main_dosage_sum = 0.0;
  double main_dosage_ssq = 0.0;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);

    uintptr_t allele_bidx = variant_bidx;
    if (max_extra_allele_ct) {
      allele_bidx = variant_bidx + CountExtraAlleles(variant_include, allele_idx_offsets, common->read_variant_uidx_starts[0], common->read_variant_uidx_starts[tidx], 0);
    }
    double* beta_se_iter = common->block_beta_se;
    if (beta_se_multiallelic_fused) {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * variant_bidx]);
    } else {
      beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct * allele_bidx]);
    }

    LinearAuxResult* block_aux_iter = &(ctx->block_aux[allele_bidx]);
    const double* local_covars_iter = nullptr;
    if (local_covar_ct) {
      // &(nullptr[0]) is okay in C++, but undefined in C
      local_covars_iter = &(ctx->local_covars_vcmaj_d[parity][variant_bidx * max_sample_ct * local_covar_ct]);
    }
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = CountSortedSmallerU32(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      uint32_t cur_variant_bidx_end = subset_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
      if (cur_variant_bidx_end > variant_bidx_end) {
        cur_variant_bidx_end = variant_bidx_end;
      }
      const uint32_t is_x = (chr_idx == x_code);
      const uint32_t is_y = (chr_idx == y_code);
      const uint32_t is_nonx_haploid = (!is_x) && IsSet(cip->haploid_mask, chr_idx);
      const uintptr_t* cur_sample_include;
      const uint32_t* cur_sample_include_cumulative_popcounts;
      const double* cur_pheno;
      const RegressionNmPrecomp* nm_precomp;
      const double* cur_covars_cmaj;
      const uintptr_t* cur_parameter_subset;
      const uintptr_t* cur_joint_test_params;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      uint32_t cur_constraint_ct;
      if (is_y && common->sample_include_y) {
        cur_sample_include = common->sample_include_y;
        cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
        cur_pheno = ctx->pheno_y_d;
        nm_precomp = common->nm_precomp_y;
        cur_covars_cmaj = ctx->covars_cmaj_y_d;
        cur_parameter_subset = common->parameter_subset_y;
        cur_joint_test_params = common->joint_test_params_y;
        cur_sample_ct = common->sample_ct_y;
        cur_covar_ct = common->covar_ct_y;
        cur_constraint_ct = common->constraint_ct_y;
      } else if (is_x && common->sample_include_x) {
        cur_sample_include = common->sample_include_x;
        cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
        cur_pheno = ctx->pheno_x_d;
        nm_precomp = common->nm_precomp_x;
        cur_covars_cmaj = ctx->covars_cmaj_x_d;
        cur_parameter_subset = common->parameter_subset_x;
        cur_joint_test_params = common->joint_test_params_x;
        cur_sample_ct = common->sample_ct_x;
        cur_covar_ct = common->covar_ct_x;
        cur_constraint_ct = common->constraint_ct_x;
      } else {
        cur_sample_include = common->sample_include;
        cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
        cur_pheno = ctx->pheno_d;
        nm_precomp = common->nm_precomp;
        cur_covars_cmaj = ctx->covars_cmaj_d;
        cur_parameter_subset = common->parameter_subset;
        cur_joint_test_params = common->joint_test_params;
        cur_sample_ct = common->sample_ct;
        cur_covar_ct = common->covar_ct;
        cur_constraint_ct = common->constraint_ct;
      }
      const uint32_t sample_ctl = BitCtToWordCt(cur_sample_ct);
      const uint32_t sample_ctl2 = NypCtToWordCt(cur_sample_ct);
      const uint32_t cur_biallelic_predictor_ct_base = 2 + domdev_present + cur_covar_ct * (1 + add_interactions * domdev_present_p1);
      uint32_t cur_biallelic_predictor_ct = cur_biallelic_predictor_ct_base;
      uint32_t literal_covar_ct = cur_covar_ct;
      if (cur_parameter_subset) {
        cur_biallelic_predictor_ct = PopcountWords(cur_parameter_subset, BitCtToWordCt(cur_biallelic_predictor_ct_base));
        literal_covar_ct = PopcountBitRange(cur_parameter_subset, 2 + domdev_present, 2 + domdev_present + cur_covar_ct);
      }
      const uint32_t max_predictor_ct = cur_biallelic_predictor_ct + max_extra_allele_ct;
      uint32_t reported_pred_uidx_biallelic_end;
      if (hide_covar) {
        if (!cur_parameter_subset) {
          reported_pred_uidx_biallelic_end = 2 + domdev_present;
        } else {
          reported_pred_uidx_biallelic_end = IsSet(cur_parameter_subset, 1) + domdev_present_p1;
        }
      } else {
        reported_pred_uidx_biallelic_end = cur_biallelic_predictor_ct;
      }
      // nm_predictors_pmaj_buf may require up to two extra columns omitted
      // from the main regression.
      // 1. In the multiallelic dominant/recessive/hethom cases, the original
      //    genotype column does not appear in the regression, and we'd rather
      //    not reconstruct it from genovec, etc. when we need to swap it out
      //    for another allele, so we keep the original genotype in an extra
      //    column.
      //    To reduce code bloat, we now handle the biallelic cases in the same
      //    way; this is one of the more peripheral code paths so adding more
      //    complexity to speed it up is less justifiable.
      // 2. If --parameters excludes the main (possibly dominant/recessive)
      //    genotype column but does care about an interaction, we want a copy
      //    of what the main genotype column's contents would have been to
      //    refer to.
      const uint32_t main_omitted = cur_parameter_subset && (!IsSet(cur_parameter_subset, 1));
      const uint32_t main_mutated = model_dominant || model_recessive || joint_hethom;
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* tmp_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      double* nm_pheno_buf = S_CAST(double*, arena_alloc_raw_rd(cur_sample_ct * sizeof(double), &workspace_iter));
      double* nm_predictors_pmaj_buf = S_CAST(double*, arena_alloc_raw_rd((max_predictor_ct + main_mutated + main_omitted) * cur_sample_ct * sizeof(double), &workspace_iter));
      double* xtx_inv = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ct * sizeof(double), &workspace_iter));
      double* fitted_coefs = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * sizeof(double), &workspace_iter));
      double* xt_y = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_corr_matrix = S_CAST(double*, arena_alloc_raw_rd((cur_biallelic_predictor_ct - 1) * (cur_biallelic_predictor_ct - 1) * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_inv_corr_sqrts = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));
      MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, arena_alloc_raw_rd(MAXV(max_predictor_ct, cur_constraint_ct) * kMatrixInvertBuf1CheckedAlloc, &workspace_iter));
      double* dbl_2d_buf = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), &workspace_iter));
      uint64_t* machr2_dosage_sums = S_CAST(uint64_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(uint64_t) * 2, &workspace_iter));
      uint64_t* machr2_dosage_ssqs = &(machr2_dosage_sums[max_extra_allele_ct + 2]);

      // could technically have this overlap fitted_coefs/xt_y, but that sets
      // the stage for future bugs
      double* inverse_corr_buf = S_CAST(double*, arena_alloc_raw_rd((max_predictor_ct - 1) * MAXV((max_predictor_ct - 1), 4) * sizeof(double), &workspace_iter));

      // joint test only
      double* tmphxs_buf = nullptr;
      double* h_transpose_buf = nullptr;
      double* inner_buf = nullptr;
      double* cur_constraints_con_major = nullptr;
      if (cur_constraint_ct) {
        tmphxs_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        h_transpose_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        inner_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * cur_constraint_ct * sizeof(double), &workspace_iter));
        cur_constraints_con_major = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        ZeroDArr(cur_constraint_ct * max_predictor_ct, cur_constraints_con_major);
        const uint32_t first_joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
        cur_constraints_con_major[first_joint_test_idx] = 1.0;
        // Rest of this matrix must be updated later, since cur_predictor_ct
        // changes at multiallelic variants.
      }
      assert(S_CAST(uintptr_t, workspace_iter - workspace_buf) == GetLinearWorkspaceSize(cur_sample_ct, cur_biallelic_predictor_ct, max_extra_allele_ct, cur_constraint_ct, main_mutated + main_omitted));
      const double pheno_ssq_base = DotprodD(cur_pheno, cur_pheno, cur_sample_ct);
      const double cur_sample_ct_recip = 1.0 / u31tod(cur_sample_ct);
      const double cur_sample_ct_m1_recip = 1.0 / u31tod(cur_sample_ct - 1);
      const uint32_t sparse_optimization_eligible = (!is_x) && nm_precomp;
      double geno_d_lookup[2];
      if (sparse_optimization_eligible) {
        geno_d_lookup[1] = 1.0;
        if (is_nonx_haploid) {
          geno_d_lookup[0] = 0.5;
        } else if (model_recessive || joint_hethom) {
          geno_d_lookup[0] = 0.0;
        } else {
          geno_d_lookup[0] = 1.0;
          if (!model_dominant) {
            geno_d_lookup[1] = 2.0;
          }
        }
      }
      const double* xtx_image = nullptr;
      const double* covarx_dotprod_inv = nullptr;
      const double* corr_inv = nullptr;
      const double* xt_y_image = nullptr;
      if (nm_precomp) {
        xtx_image = nm_precomp->xtx_image;
        covarx_dotprod_inv = nm_precomp->covarx_dotprod_inv;
        corr_inv = nm_precomp->corr_inv;
        const uintptr_t nongeno_pred_ct = cur_biallelic_predictor_ct - domdev_present - 2;
        const uintptr_t nonintercept_biallelic_pred_ct = cur_biallelic_predictor_ct - 1;
        memcpy(semicomputed_biallelic_corr_matrix, nm_precomp->corr_image, nonintercept_biallelic_pred_ct * nonintercept_biallelic_pred_ct * sizeof(double));
        memcpy(&(semicomputed_biallelic_inv_corr_sqrts[domdev_present_p1]), nm_precomp->corr_inv_sqrts, nongeno_pred_ct * sizeof(double));
        xt_y_image = nm_precomp->xt_y_image;
      }
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      // when this is set, the last fully-processed variant had no missing
      // genotypes, and if the current variant also has no missing genotypes we
      // may be able to skip reinitialization of most of
      // nm_predictors_pmaj_buf.
      uint32_t prev_nm = 0;

      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        if (allele_idx_offsets) {
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
          if (!beta_se_multiallelic_fused) {
            extra_regression_ct = allele_ct - 2;
          }
        }
        const uint32_t allele_ct_m2 = allele_ct - 2;
        const uint32_t expected_predictor_ct = cur_biallelic_predictor_ct + allele_ct_m2;
        PglErr reterr;
        if (!allele_ct_m2) {
          reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
        } else {
          reterr = PgrGetMD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, &pgv);
          // todo: proper multiallelic dosage support
        }
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto GlmLinearThread_err;
        }
        ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
        GenoarrCountFreqsUnsafe(pgv.genovec, cur_sample_ct, genocounts);
        uint32_t missing_ct = genocounts[3];
        if (!missing_ct) {
          SetAllBits(cur_sample_ct, sample_nm);
        } else {
          GenoarrToNonmissing(pgv.genovec, cur_sample_ct, sample_nm);
          if (pgv.dosage_ct) {
            BitvecOr(pgv.dosage_present, sample_ctl, sample_nm);
            missing_ct = cur_sample_ct - PopcountWords(sample_nm, sample_ctl);
          }
        }
        if (omitted_alleles) {
          omitted_allele_idx = omitted_alleles[variant_uidx];
        }
        uintptr_t const_alleles[DivUp(kPglMaxAltAlleleCt + 1, kBitsPerWord)];
        const uint32_t allele_ctl = DivUp(allele_ct, kBitsPerWord);
        ZeroWArr(allele_ctl, const_alleles);
        const uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
        const uint32_t nm_sample_ctl = BitCtToWordCt(nm_sample_ct);
        // first predictor column: intercept
        if (!prev_nm) {
          for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
            nm_predictors_pmaj_buf[sample_idx] = 1.0;
          }
        }
        // second predictor column: main genotype
        double* genotype_vals = &(nm_predictors_pmaj_buf[nm_sample_ct]);
        if (main_mutated || main_omitted) {
          genotype_vals = &(nm_predictors_pmaj_buf[expected_predictor_ct * nm_sample_ct]);
        }
        double cur_pheno_ssq = pheno_ssq_base;
        uintptr_t sample_midx_base = 0;
        uintptr_t sample_nm_inv_bits = ~sample_nm[0];
        for (uint32_t missing_idx = 0; missing_idx != missing_ct; ++missing_idx) {
          const uintptr_t sample_midx = BitIter0(sample_nm, &sample_midx_base, &sample_nm_inv_bits);
          cur_pheno_ssq -= cur_pheno[sample_midx] * cur_pheno[sample_midx];
        }
        uint32_t sparse_optimization = 0;
        double* multi_start = nullptr;
        if (!allele_ct_m2) {
          // When prev_nm is set and missing_ct is zero, we don't need to call
          // MultiplySelfTranspose() on nm_predictors_pmaj_buf to get all the
          // predictor x predictor dot products; instead we patch in the
          // genotype x predictor dot products that may change, copying the
          // rest from xtx_image.
          // As a side effect, it is no longer strictly necessary to fill the
          // genotype row of nm_predictors_pmaj_buf.  sparse_optimization
          // indicates that plink 1.9's QT --assoc sparse dot product algorithm
          // will be used instead.
          // probable todos: allow a few dosages to be present, cover chrX
          // case.
          if (omitted_allele_idx) {
            GenovecInvertUnsafe(cur_sample_ct, pgv.genovec);
            ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
            }
            const uint32_t uii = genocounts[0];
            genocounts[0] = genocounts[2];
            genocounts[2] = uii;
          }
          uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
          uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
          if (!missing_ct) {
            // originally had genocounts[0] > 0.875 * nm_sample_ct threshold,
            // but then tried this on high-MAF data and it was still
            // substantially faster
            sparse_optimization = sparse_optimization_eligible && (!pgv.dosage_ct) && prev_nm;
            if (!sparse_optimization) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, genotype_vals);
              if (pgv.dosage_ct) {
                uintptr_t sample_idx_base = 0;
                uintptr_t dosage_present_bits = pgv.dosage_present[0];
                for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                  const uintptr_t sample_idx = BitIter1(pgv.dosage_present, &sample_idx_base, &dosage_present_bits);
                  const uintptr_t dosage_val = pgv.dosage_main[dosage_idx];
                  // 32768 -> 2, 16384 -> 1, 0 -> 0
                  genotype_vals[sample_idx] = kRecipDosageMid * swtod(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_idx);
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                }
              }
            }
          } else {
            if (!pgv.dosage_ct) {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, genotype_vals);
            } else {
              sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              uint32_t dosage_idx = 0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_midx);
                double cur_val;
                if (IsSet(pgv.dosage_present, sample_midx)) {
                  const uintptr_t dosage_val = pgv.dosage_main[dosage_idx++];
                  cur_val = kRecipDosageMid * swtod(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                } else {
                  cur_val = kSmallDoubles[cur_geno];
                }
                genotype_vals[sample_idx] = cur_val;
              }
            }
          }
          // Check for constant genotype column.
          // (Technically, we should recheck later in the chrX no-sex-covariate
          // --xchr-model 1 corner case.)
          if (!pgv.dosage_ct) {
            if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          } else if (pgv.dosage_ct == nm_sample_ct) {
            if (DosageIsConstant(dosage_sum, dosage_ssq, nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          }
          machr2_dosage_sums[1 - omitted_allele_idx] = dosage_sum;
          machr2_dosage_ssqs[1 - omitted_allele_idx] = dosage_ssq;
          machr2_dosage_sums[omitted_allele_idx] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[omitted_allele_idx] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
        } else {
          // multiallelic.
          // If some but not all alleles have constant dosages, we remove just
          // those alleles from the regressions; trim-alts is not necessary to
          // see what's going on with the other alleles.  To reduce parsing
          // complexity, the number of output lines is not affected by this;
          // the ones corresponding to the constant alleles have NA values.
          // Punt on sparse_optimization for now; may be worth revisiting after
          // multiallelic dosage implemented.
          // dosage_ct == 0 temporarily guaranteed if we reach here.
          assert(!pgv.dosage_ct);
          multi_start = &(nm_predictors_pmaj_buf[(expected_predictor_ct - allele_ct_m2) * nm_sample_ct]);
          ZeroU64Arr(allele_ct, machr2_dosage_sums);
          ZeroU64Arr(allele_ct, machr2_dosage_ssqs);
          // postpone multiply for now, since no multiallelic dosages
          machr2_dosage_sums[0] = genocounts[1] + 2 * genocounts[0];
          machr2_dosage_ssqs[0] = genocounts[1] + 4LLU * genocounts[0];
          if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
            SetBit(0, const_alleles);
          }
          if (omitted_allele_idx) {
            // Main genotype column starts as REF.
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallInvDoublePairs, nm_sample_ct, genotype_vals);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallInvDoubles, cur_sample_ct, genotype_vals);
            }
          }
          uint32_t rare_allele_ct = allele_ct_m2;
          double* alt1_start = nullptr;
          double* rarealt_start = multi_start;
          if (omitted_allele_idx != 1) {
            if (omitted_allele_idx) {
              alt1_start = multi_start;
              rarealt_start = &(rarealt_start[nm_sample_ct]);
              --rare_allele_ct;
            } else {
              alt1_start = genotype_vals;
            }
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, alt1_start);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, alt1_start);
            }
          }
          ZeroDArr(rare_allele_ct * nm_sample_ct, rarealt_start);
          // Use sums as ones[] and ssqs as twos[] for rarealts; transform to
          // actual sums/ssqs later.
          if (pgv.patch_01_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_01_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_01_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ct + sample_idx] = 1.0;
                alt1_start[sample_idx] = 0.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ct + sample_idx] = 1.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                alt1_start[sample_idx] = 0.0;
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                machr2_dosage_sums[allele_code] += 1;
                if (allele_code == omitted_allele_idx) {
                  continue;
                }
                const uint32_t cur_col = allele_code - 2 - (allele_code > omitted_allele_idx);
                rarealt_start[cur_col * nm_sample_ct + sample_idx] = 1.0;
              }
            }
          }
          uintptr_t alt1_het_ct = genocounts[1] - pgv.patch_01_ct;
          if (pgv.patch_10_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_10_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_10_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] = 2.0;
                  alt1_start[sample_idx] = 0.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ct + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] += 1.0;
                    alt1_start[sample_idx] = 0.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] = 2.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ct + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] += 1.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  machr2_dosage_ssqs[ac0] += 1;
                  alt1_start[sample_idx] = 0.0;
                  if (ac0 != omitted_allele_idx) {
                    const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                    rarealt_start[ac0_col * nm_sample_ct + sample_idx] = 2.0;
                  }
                } else {
                  machr2_dosage_sums[ac1] += 1;
                  if (ac1 != omitted_allele_idx) {
                    const uint32_t ac1_col = ac1 - 2 - (ac1 > omitted_allele_idx);
                    rarealt_start[ac1_col * nm_sample_ct + sample_idx] = 1.0;
                  }
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    machr2_dosage_sums[ac0] += 1;
                    alt1_start[sample_idx] = 0.0;
                    if (ac0 != omitted_allele_idx) {
                      const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                      rarealt_start[ac0_col * nm_sample_ct + sample_idx] += 1.0;
                    }
                  }
                }
              }
            }
          }
          for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
            const uintptr_t one_ct = machr2_dosage_sums[allele_idx];
            const uintptr_t two_ct = machr2_dosage_ssqs[allele_idx];
            machr2_dosage_sums[allele_idx] = one_ct + 2 * two_ct;
            machr2_dosage_ssqs[allele_idx] = one_ct + 4LLU * two_ct;
            if ((one_ct == nm_sample_ct) || (two_ct == nm_sample_ct) || ((!one_ct) && (!two_ct))) {
              SetBit(allele_idx, const_alleles);
            }
          }
          const uintptr_t alt1_hom_ct = genocounts[2] - pgv.patch_10_ct;
          machr2_dosage_sums[1] = alt1_het_ct + 2 * alt1_hom_ct;
          machr2_dosage_ssqs[1] = alt1_het_ct + 4LLU * alt1_hom_ct;
          if ((alt1_het_ct == nm_sample_ct) || (alt1_hom_ct == nm_sample_ct) || ((!alt1_het_ct) && (!alt1_hom_ct))) {
            SetBit(1, const_alleles);
          }
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            machr2_dosage_sums[allele_idx] *= 0x4000LLU;
            machr2_dosage_ssqs[allele_idx] *= 0x10000000LLU;
          }
        }
        // usually need to save some of {sample_obs_ct, allele_obs_ct,
        // a1_dosage, mach_r2 even for skipped variants
        // compute them all for now, could conditionally skip later
        uint32_t allele_obs_ct = nm_sample_ct * 2;
        if (!is_x) {
          if (is_nonx_haploid) {
            allele_obs_ct = nm_sample_ct;
            // everything is on 0..1 scale, not 0..2
            if (!sparse_optimization) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                genotype_vals[sample_idx] *= 0.5;
              }
              const uint32_t high_ct = nm_sample_ct * allele_ct_m2;
              for (uint32_t uii = 0; uii != high_ct; ++uii) {
                multi_start[uii] *= 0.5;
              }
            }
          }
        } else {
          CopyBitarrSubset(sex_male_collapsed, sample_nm, nm_sample_ct, tmp_nm);
          const uintptr_t* male_nm = tmp_nm;
          const uint32_t nm_male_ct = PopcountWords(male_nm, nm_sample_ctl);
          if (is_xchr_model_1) {
            // special case: multiply male values by 0.5
            uintptr_t sample_idx_base = 0;
            uintptr_t male_nm_bits = male_nm[0];
            for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
              const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
              genotype_vals[sample_idx] *= 0.5;
              // could insert multiallelic loop here instead, but I'm guessing
              // that's worse due to locality of writes?
            }
            for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
              double* cur_start = &(multi_start[extra_allele_idx * nm_sample_ct]);
              sample_idx_base = 0;
              male_nm_bits = male_nm[0];
              for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
                const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
                cur_start[sample_idx] *= 0.5;
              }
            }
            allele_obs_ct -= nm_male_ct;
          }
        }
        const double mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, allele_ct);
        for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
          if (allele_idx == omitted_allele_idx) {
            continue;
          }
          block_aux_iter->sample_obs_ct = nm_sample_ct;
          block_aux_iter->allele_obs_ct = allele_obs_ct;
          double a1_dosage = u63tod(machr2_dosage_sums[allele_idx]) * kRecipDosageMid;
          if (is_xchr_model_1) {
            // ugh.
            double* geno_col = genotype_vals;
            if (allele_idx > (!omitted_allele_idx)) {
              geno_col = &(nm_predictors_pmaj_buf[(expected_predictor_ct - (allele_ct - allele_idx) + (allele_idx < omitted_allele_idx)) * nm_sample_ct]);
            }
            a1_dosage = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              a1_dosage += geno_col[sample_idx];
            }
            if (!allele_ct_m2) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = 0.0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const double cur_dosage = geno_col[sample_idx];
                main_dosage_ssq += cur_dosage * cur_dosage;
              }
            }
          } else {
            if (is_nonx_haploid) {
              a1_dosage *= 0.5;
            }
            if (!allele_ct_m2) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = u63tod(machr2_dosage_ssqs[allele_idx]) * kRecipDosageMidSq;
              if (is_nonx_haploid) {
                main_dosage_ssq *= 0.25;
              }
            }
          }
          block_aux_iter->a1_dosage = a1_dosage;
          block_aux_iter->mach_r2 = mach_r2;
          ++block_aux_iter;
        }
        // now free to skip the actual regression if there are too few samples,
        // or there's a zero-variance genotype column
        GlmErr glm_err = 0;
        if (nm_sample_ct <= expected_predictor_ct) {
          glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
        } else if (IsSet(const_alleles, omitted_allele_idx)) {
          glm_err = SetGlmErr0(kGlmErrcodeConstOmittedAllele);
        }
        if (glm_err) {
          if (missing_ct) {
            // covariates have not been copied yet, so we can't usually change
            // prev_nm from 0 to 1 when missing_ct == 0 (and there's little
            // reason to optimize the zero-covariate case).
            prev_nm = 0;
          }
          uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
          if (beta_se_multiallelic_fused || (!hide_covar)) {
            reported_ct += allele_ct_m2;
          }
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            for (uint32_t uii = 0; uii != reported_ct; ++uii) {
              memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
              beta_se_iter[uii * 2 + 1] = -9.0;
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        } else {
          uint32_t parameter_uidx = 2 + domdev_present;
          double* nm_predictors_pmaj_istart = nullptr;
          if (!sparse_optimization) {
            // only need to do this part once per variant in multiallelic case
            double* nm_predictors_pmaj_iter = &(nm_predictors_pmaj_buf[nm_sample_ct * (parameter_uidx - main_omitted)]);
            if (missing_ct || (!prev_nm)) {
              // fill phenotype
              sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                nm_pheno_buf[sample_idx] = cur_pheno[sample_midx];
              }

              // fill covariates
              for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx, ++parameter_uidx) {
                // strictly speaking, we don't need cur_covars_cmaj to be
                // vector-aligned
                if (cur_parameter_subset && (!IsSet(cur_parameter_subset, parameter_uidx))) {
                  continue;
                }
                const double* cur_covar_col;
                if (covar_idx < local_covar_ct) {
                  cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                } else {
                  cur_covar_col = &(cur_covars_cmaj[(covar_idx - local_covar_ct) * cur_sample_ct]);
                }
                sample_midx_base = 0;
                sample_nm_bits = sample_nm[0];
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                  *nm_predictors_pmaj_iter++ = cur_covar_col[sample_midx];
                }
              }
              nm_predictors_pmaj_istart = nm_predictors_pmaj_iter;
              prev_nm = !missing_ct;
            } else {
              // bugfix (15 Aug 2018): this was not handling --parameters
              // correctly when a covariate was only needed as part of an
              // interaction
              parameter_uidx += cur_covar_ct;
              nm_predictors_pmaj_istart = &(nm_predictors_pmaj_iter[literal_covar_ct * nm_sample_ct]);
            }
          }
          const uint32_t const_allele_ct = PopcountWords(const_alleles, allele_ctl);
          if (const_allele_ct) {
            // Must delete constant-allele columns from nm_predictors_pmaj, and
            // shift later columns back.
            double* read_iter = genotype_vals;
            double* write_iter = genotype_vals;
            for (uint32_t read_allele_idx = 0; read_allele_idx != allele_ct; ++read_allele_idx) {
              if (read_allele_idx == omitted_allele_idx) {
                continue;
              }
              if (!IsSet(const_alleles, read_allele_idx)) {
                if (write_iter != read_iter) {
                  memcpy(write_iter, read_iter, nm_sample_ct * sizeof(double));
                }
                if (write_iter == genotype_vals) {
                  write_iter = multi_start;
                } else {
                  write_iter = &(write_iter[nm_sample_ct]);
                }
              }
              if (read_iter == genotype_vals) {
                read_iter = multi_start;
              } else {
                read_iter = &(read_iter[nm_sample_ct]);
              }
            }
          }
          const uint32_t cur_predictor_ct = expected_predictor_ct - const_allele_ct;
          uint32_t nonconst_extra_regression_idx = UINT32_MAX;  // deliberate overflow
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            if (extra_regression_ct) {
              if (IsSet(const_alleles, extra_regression_idx + (extra_regression_idx >= omitted_allele_idx))) {
                glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
                goto GlmLinearThread_skip_regression;
              }
              ++nonconst_extra_regression_idx;
              if (nonconst_extra_regression_idx) {
                double* swap_target = &(multi_start[(nonconst_extra_regression_idx - 1) * nm_sample_ct]);
                for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
                  double dxx = genotype_vals[uii];
                  genotype_vals[uii] = swap_target[uii];
                  swap_target[uii] = dxx;
                }
              }
            }
            if (!sparse_optimization) {
              double* main_vals = &(nm_predictors_pmaj_buf[nm_sample_ct]);
              double* domdev_vals = nullptr;
              if (main_omitted) {
                // if main_mutated, this will be filled below
                // if not, this aliases genotype_vals
                main_vals = &(nm_predictors_pmaj_buf[(cur_predictor_ct + main_mutated) * nm_sample_ct]);
              } else if (joint_genotypic || joint_hethom) {
                domdev_vals = &(main_vals[nm_sample_ct]);
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 2.0 - cur_genotype_val;
                  }
                  domdev_vals[sample_idx] = cur_genotype_val;
                }
              }
              if (model_dominant) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..1..1
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 1.0;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              } else if (model_recessive || joint_hethom) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..0..1
                  if (cur_genotype_val < 1.0) {
                    cur_genotype_val = 0.0;
                  } else {
                    cur_genotype_val -= 1.0;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              }

              // fill interaction terms
              if (add_interactions) {
                double* nm_predictors_pmaj_iter = nm_predictors_pmaj_istart;
                for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
                  const double* cur_covar_col;
                  if (covar_idx < local_covar_ct) {
                    cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                  } else {
                    cur_covar_col = &(cur_covars_cmaj[covar_idx * cur_sample_ct]);
                  }
                  if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                    sample_midx_base = 0;
                    uintptr_t sample_nm_bits = sample_nm[0];
                    for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                      const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                      *nm_predictors_pmaj_iter++ = main_vals[sample_idx] * cur_covar_col[sample_midx];
                    }
                  }
                  ++parameter_uidx;
                  if (domdev_present) {
                    if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                      sample_midx_base = 0;
                      uintptr_t sample_nm_bits = sample_nm[0];
                      for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                        const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                        *nm_predictors_pmaj_iter++ = domdev_vals[sample_idx] * cur_covar_col[sample_midx];
                      }
                    }
                    ++parameter_uidx;
                  }
                }
              }
            }

            // bugfix (12 Sep 2017): forgot to implement per-variant VIF and
            // max-corr checks
            if (xtx_image && prev_nm && (!allele_ct_m2)) {
              // only need to fill in additive and possibly domdev dot
              // products
              memcpy(xtx_inv, xtx_image, cur_predictor_ct * cur_predictor_ct * sizeof(double));
              memcpy(xt_y, xt_y_image, cur_predictor_ct * sizeof(double));
              if (sparse_optimization) {
                // currently does not handle chrX
                double geno_pheno_prod = 0.0;
                double domdev_pheno_prod = 0.0;
                double domdev_geno_prod = 0.0;
                double* geno_dotprod_row = &(xtx_inv[cur_predictor_ct]);
                double* domdev_dotprod_row = &(xtx_inv[2 * cur_predictor_ct]);
                for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
                  uintptr_t geno_word = pgv.genovec[widx];
                  if (geno_word) {
                    const uint32_t sample_idx_base = widx * kBitsPerWordD2;
                    do {
                      const uint32_t lowest_set_bit = ctzw(geno_word);
                      // since there are no missing values, we have a het if
                      // (lowest_set_bit & 1) is zero, and a hom-alt when
                      // it's one.
                      const uint32_t sample_idx = sample_idx_base + (lowest_set_bit / 2);
                      const double geno_d = geno_d_lookup[lowest_set_bit & 1];
                      const double cur_pheno_val = nm_pheno_buf[sample_idx];
                      geno_pheno_prod += geno_d * cur_pheno_val;
                      for (uintptr_t pred_idx = domdev_present + 2; pred_idx != cur_predictor_ct; ++pred_idx) {
                        geno_dotprod_row[pred_idx] += geno_d * nm_predictors_pmaj_buf[pred_idx * nm_sample_ct + sample_idx];
                      }
                      // can have a separate categorical loop here

                      if (domdev_present && (!(lowest_set_bit & 1))) {
                        // domdev = 1
                        domdev_pheno_prod += cur_pheno_val;
                        domdev_geno_prod += geno_d;
                        for (uintptr_t pred_idx = 3; pred_idx != cur_predictor_ct; ++pred_idx) {
                          domdev_dotprod_row[pred_idx] += nm_predictors_pmaj_buf[pred_idx * nm_sample_ct + sample_idx];
                        }
                        // categorical optimization possible here
                      }
                      geno_word &= geno_word - 1;
                    } while (geno_word);
                  }
                }
                xt_y[1] = geno_pheno_prod;
                const double het_ctd = u31tod(genocounts[1]);
                const double homalt_ctd = u31tod(genocounts[2]);
                xtx_inv[cur_predictor_ct] = het_ctd * geno_d_lookup[0] + homalt_ctd * geno_d_lookup[1];
                xtx_inv[cur_predictor_ct + 1] = het_ctd * geno_d_lookup[0] * geno_d_lookup[0] + homalt_ctd * geno_d_lookup[1] * geno_d_lookup[1];
                if (domdev_present) {
                  xt_y[2] = domdev_pheno_prod;
                  xtx_inv[cur_predictor_ct + 2] = domdev_geno_prod;
                  xtx_inv[2 * cur_predictor_ct] = het_ctd;
                  xtx_inv[2 * cur_predictor_ct + 2] = het_ctd;
                }
              } else {
                // !sparse_optimization
                xt_y[1] = DotprodD(&(nm_predictors_pmaj_buf[nm_sample_ct]), nm_pheno_buf, nm_sample_ct);
                uintptr_t start_pred_idx = 0;
                if (!(model_dominant || model_recessive || joint_hethom)) {
                  start_pred_idx = domdev_present + 2;
                  xtx_inv[cur_predictor_ct] = main_dosage_sum;
                  xtx_inv[cur_predictor_ct + 1] = main_dosage_ssq;
                }
                if (cur_predictor_ct > start_pred_idx) {
                  // categorical optimization possible here
                  ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[nm_sample_ct]), &(nm_predictors_pmaj_buf[start_pred_idx * nm_sample_ct]), nm_sample_ct, nm_sample_ct, cur_predictor_ct - start_pred_idx, &(xtx_inv[cur_predictor_ct + start_pred_idx]));
                }
                if (domdev_present) {
                  xt_y[2] = DotprodD(&(nm_predictors_pmaj_buf[2 * nm_sample_ct]), nm_pheno_buf, nm_sample_ct);
                  // categorical optimization possible here
                  ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[2 * nm_sample_ct]), nm_predictors_pmaj_buf, nm_sample_ct, nm_sample_ct, cur_predictor_ct, &(xtx_inv[2 * cur_predictor_ct]));
                  xtx_inv[cur_predictor_ct + 2] = xtx_inv[2 * cur_predictor_ct + 1];
                }
              }
              glm_err = CheckMaxCorrAndVifNm(xtx_inv, corr_inv, cur_predictor_ct, domdev_present_p1, cur_sample_ct_recip, cur_sample_ct_m1_recip, max_corr, vif_thresh, semicomputed_biallelic_corr_matrix, semicomputed_biallelic_inv_corr_sqrts, dbl_2d_buf, &(dbl_2d_buf[2 * cur_predictor_ct]), &(dbl_2d_buf[3 * cur_predictor_ct]));
              if (glm_err) {
                goto GlmLinearThread_skip_regression;
              }
              const double geno_ssq = xtx_inv[1 + cur_predictor_ct];
              if (!domdev_present) {
                xtx_inv[1 + cur_predictor_ct] = xtx_inv[cur_predictor_ct];
                if (InvertRank1Symm(covarx_dotprod_inv, &(xtx_inv[1 + cur_predictor_ct]), cur_predictor_ct - 1, 1, geno_ssq, dbl_2d_buf, inverse_corr_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                  goto GlmLinearThread_skip_regression;
                }
              } else {
                const double domdev_geno_prod = xtx_inv[2 + cur_predictor_ct];
                const double domdev_ssq = xtx_inv[2 + 2 * cur_predictor_ct];
                xtx_inv[2 + cur_predictor_ct] = xtx_inv[cur_predictor_ct];
                xtx_inv[2 + 2 * cur_predictor_ct] = xtx_inv[2 * cur_predictor_ct];
                if (InvertRank2Symm(covarx_dotprod_inv, &(xtx_inv[2 + cur_predictor_ct]), cur_predictor_ct - 2, cur_predictor_ct, 1, geno_ssq, domdev_geno_prod, domdev_ssq, dbl_2d_buf, inverse_corr_buf, &(inverse_corr_buf[2 * (cur_predictor_ct - 2)]))) {
                  glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                  goto GlmLinearThread_skip_regression;
                }
              }
              // need to make sure xtx_inv remains reflected in NOLAPACK case
              memcpy(xtx_inv, dbl_2d_buf, cur_predictor_ct * cur_predictor_ct * sizeof(double));
              ReflectMatrix(cur_predictor_ct, xtx_inv);
              ColMajorVectorMatrixMultiplyStrided(xt_y, xtx_inv, cur_predictor_ct, cur_predictor_ct, cur_predictor_ct, fitted_coefs);
            } else {
              // generic case
              // major categorical optimization possible here, some
              // multiallelic optimizations possible
              MultiplySelfTranspose(nm_predictors_pmaj_buf, cur_predictor_ct, nm_sample_ct, xtx_inv);

              for (uint32_t pred_idx = 1; pred_idx != cur_predictor_ct; ++pred_idx) {
                dbl_2d_buf[pred_idx] = xtx_inv[pred_idx * cur_predictor_ct];
              }
              glm_err = CheckMaxCorrAndVif(xtx_inv, 1, cur_predictor_ct, nm_sample_ct, max_corr, vif_thresh, dbl_2d_buf, nullptr, inverse_corr_buf, inv_1d_buf);
              if (glm_err) {
                goto GlmLinearThread_skip_regression;
              }
              if (LinearRegressionInv(nm_pheno_buf, nm_predictors_pmaj_buf, cur_predictor_ct, nm_sample_ct, 1, xtx_inv, fitted_coefs, xt_y, inv_1d_buf, dbl_2d_buf)) {
                glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                goto GlmLinearThread_skip_regression;
              }
            }
            {
              // RSS = y^T y - y^T X (X^T X)^{-1} X^T y
              //     = cur_pheno_ssq - xt_y * fitted_coefs
              // s^2 = RSS / df
              // possible todo: improve numerical stability of this computation
              // in non-mean-centered phenotype case
              const double sigma = (cur_pheno_ssq - DotprodxD(xt_y, fitted_coefs, cur_predictor_ct)) / u31tod(nm_sample_ct - cur_predictor_ct);
              for (uint32_t uii = 0; uii != cur_predictor_ct; ++uii) {
                double* s_iter = &(xtx_inv[uii * cur_predictor_ct]);
#ifdef NOLAPACK
                for (uint32_t ujj = 0; ujj != cur_predictor_ct; ++ujj) {
                  s_iter[ujj] *= sigma;
                }
#else
                for (uint32_t ujj = 0; ujj <= uii; ++ujj) {
                  s_iter[ujj] *= sigma;
                }
#endif
              }
              // validParameters() check
              for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
                const double xtx_inv_diag_element = xtx_inv[pred_uidx * (cur_predictor_ct + 1)];
                if (xtx_inv_diag_element < 1e-20) {
                  glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                  goto GlmLinearThread_skip_regression;
                }
                // use dbl_2d_buf[] to store diagonal square roots
                dbl_2d_buf[pred_uidx] = sqrt(xtx_inv_diag_element);
              }
              dbl_2d_buf[0] = sqrt(xtx_inv[0]);
              for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
                const double cur_xtx_inv_diag_sqrt = 0.99999 * dbl_2d_buf[pred_uidx];
                const double* xtx_inv_row = &(xtx_inv[pred_uidx * cur_predictor_ct]);
                for (uint32_t pred_uidx2 = 0; pred_uidx2 != pred_uidx; ++pred_uidx2) {
                  if (xtx_inv_row[pred_uidx2] > cur_xtx_inv_diag_sqrt * dbl_2d_buf[pred_uidx2]) {
                    glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                    goto GlmLinearThread_skip_regression;
                  }
                }
              }
              double* beta_se_iter2 = beta_se_iter;
              for (uint32_t pred_uidx = reported_pred_uidx_start; pred_uidx != reported_pred_uidx_biallelic_end; ++pred_uidx) {
                // In the multiallelic-fused case, if the first allele is
                // constant, this writes the beta/se values for the first
                // nonconstant, non-omitted allele where the results for the
                // first allele belong.  We correct that below.
                *beta_se_iter2++ = fitted_coefs[pred_uidx];
                *beta_se_iter2++ = dbl_2d_buf[pred_uidx];
              }
              // move this up since dbl_2d_buf may be clobbered
              if (cur_constraint_ct) {
                const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeRankDeficient);
                memcpy(beta_se_iter2, &glm_err2, 8);
                ++beta_se_iter2;
                *beta_se_iter2++ = -9.0;
              }
              if (!const_allele_ct) {
                if (beta_se_multiallelic_fused || (!hide_covar)) {
                  for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
                    beta_se_iter2[2 * extra_allele_idx] = fitted_coefs[cur_biallelic_predictor_ct + extra_allele_idx];
                    beta_se_iter2[2 * extra_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + extra_allele_idx];
                  }
                }
              } else if (!beta_se_multiallelic_fused) {
                if (!hide_covar) {
                  // Need to insert some {CONST_ALLELE, -9} entries.
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                  const uint32_t cur_raw_allele_idx = extra_regression_idx + (extra_regression_idx >= omitted_allele_idx);
                  uint32_t extra_read_allele_idx = 0;
                  uint32_t extra_write_allele_idx = 0;
                  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    if ((allele_idx == omitted_allele_idx) || (allele_idx == cur_raw_allele_idx)) {
                      continue;
                    }
                    if (IsSet(const_alleles, allele_idx)) {
                      memcpy(&(beta_se_iter2[2 * extra_write_allele_idx]), &glm_err2, 8);
                      beta_se_iter2[2 * extra_write_allele_idx + 1] = -9.0;
                    } else {
                      beta_se_iter2[2 * extra_write_allele_idx] = fitted_coefs[cur_biallelic_predictor_ct + extra_read_allele_idx];
                      beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + extra_read_allele_idx];
                      ++extra_read_allele_idx;
                    }
                    ++extra_write_allele_idx;
                  }
                }
              } else {
                const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                // Special-case first nonconst allele since it's positioned
                // discontinuously, and its BETA/SE may already be correctly
                // filled.
                uint32_t allele_idx = omitted_allele_idx? 0 : 1;
                uint32_t extra_write_allele_idx = 0;
                if (IsSet(const_alleles, allele_idx)) {
                  memcpy(&(beta_se_iter[2 * include_intercept]), &glm_err2, 8);
                  beta_se_iter[2 * include_intercept + 1] = -9.0;
                  allele_idx = AdvTo0Bit(const_alleles, 1);
                  if (allele_idx == omitted_allele_idx) {
                    allele_idx = AdvTo0Bit(const_alleles, omitted_allele_idx + 1);
                  }
                  extra_write_allele_idx = allele_idx - 1 - (allele_idx > omitted_allele_idx);
                  for (uint32_t uii = 0; uii != extra_write_allele_idx; ++uii) {
                    memcpy(&(beta_se_iter2[2 * uii]), &glm_err2, 8);
                    beta_se_iter2[2 * uii + 1] = -9.0;
                  }
                  beta_se_iter2[2 * extra_write_allele_idx] = fitted_coefs[1];
                  beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[1];
                  ++extra_write_allele_idx;
                }
                ++allele_idx;
                uint32_t nonconst_allele_idx_m1 = 0;
                for (; allele_idx != allele_ct; ++allele_idx) {
                  if (allele_idx == omitted_allele_idx) {
                    continue;
                  }
                  if (!IsSet(const_alleles, allele_idx)) {
                    beta_se_iter2[2 * extra_write_allele_idx] = fitted_coefs[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                    beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                    ++nonconst_allele_idx_m1;
                  } else {
                    memcpy(&(beta_se_iter2[2 * extra_write_allele_idx]), &glm_err2, 8);
                    beta_se_iter2[2 * extra_write_allele_idx + 1] = -9.0;
                  }
                  ++extra_write_allele_idx;
                }
              }
              if (cur_constraint_ct) {
                uint32_t joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 1.0;
                }
#ifndef NOLAPACK
                // xtx_inv upper triangle was not filled
                ReflectMatrix(cur_predictor_ct, xtx_inv);
#endif
                double chisq;
                if (!LinearHypothesisChisq(fitted_coefs, cur_constraints_con_major, xtx_inv, cur_constraint_ct, cur_predictor_ct, &chisq, tmphxs_buf, h_transpose_buf, inner_buf, inv_1d_buf, dbl_2d_buf)) {
                  beta_se_iter2[-1] = chisq;
                }
                // next test may have different alt allele count
                joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                  cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 0.0;
                }
              }
            }
            while (0) {
            GlmLinearThread_skip_regression:
              {
                uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
                if (beta_se_multiallelic_fused || (!hide_covar)) {
                  reported_ct += allele_ct_m2;
                }
                for (uint32_t uii = 0; uii != reported_ct; ++uii) {
                  memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
                  beta_se_iter[uii * 2 + 1] = -9.0;
                }
              }
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        }
        // bugfix (1 Apr 2019): this needs to execute when regression is
        // skipped
        if (local_covars_iter) {
          local_covars_iter = &(local_covars_iter[local_covar_ct * max_sample_ct]);
        }
      }
    }
    while (0) {
    GlmLinearThread_err:
      UpdateU64IfSmaller(new_err_info, &common->err_info);
      break;
    }
    parity = 1 - parity;
    variant_idx_offset += cur_block_variant_ct;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr GlmLinear(const char* cur_pheno_name, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, const char* outname, uint32_t raw_variant_ct, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLinearCtx* ctx, TextStream* local_covar_txsp, uintptr_t* valid_variants, uintptr_t* valid_alleles, double* orig_ln_pvals, uintptr_t* valid_allele_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  ThreadGroup tg;
  PreinitCstream(&css);
  PreinitThreads(&tg);
  {
    GlmCtx* common = ctx->common;
    const uintptr_t* variant_include = common->variant_include;
    const ChrInfo* cip = common->cip;
    const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
    const AlleleCode* omitted_alleles = common->omitted_alleles;

    const uint32_t sample_ct = common->sample_ct;
    const uint32_t sample_ct_x = common->sample_ct_x;
    const uint32_t sample_ct_y = common->sample_ct_y;
    const uint32_t covar_ct = common->covar_ct;
    const uintptr_t local_covar_ct = common->local_covar_ct;
    const uint32_t covar_ct_x = common->covar_ct_x;
    const uint32_t covar_ct_y = common->covar_ct_y;

    uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
    if (max_sample_ct < sample_ct_y) {
      max_sample_ct = sample_ct_y;
    }
    uint32_t* local_sample_idx_order = nullptr;
    uint32_t local_line_idx = 0;
    uint32_t local_xy = 0;  // 1 = chrX, 2 = chrY

    const char* local_line_iter = nullptr;
    uint32_t local_prev_chr_code = UINT32_MAX;
    uint32_t local_chr_code = UINT32_MAX;
    uint32_t local_bp = UINT32_MAX;
    uint32_t local_skip_chr = 1;
    if (local_covar_ct) {
      reterr = TextRewind(local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLinear_ret_TSTREAM_FAIL;
      }
      local_line_idx = glm_info_ptr->local_header_line_ct;
      reterr = TextSkip(local_line_idx, local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLinear_ret_TSTREAM_FAIL;
      }
      if (unlikely(bigstack_alloc_u32(local_sample_ct, &local_sample_idx_order))) {
        goto GlmLinear_ret_NOMEM;
      }
      for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
        const uint32_t cur_uidx = local_sample_uidx_order[uii];
        uint32_t cur_idx = UINT32_MAX;
        if ((cur_uidx != UINT32_MAX) && IsSet(common->sample_include, cur_uidx)) {
          cur_idx = RawToSubsettedPos(common->sample_include, common->sample_include_cumulative_popcounts, cur_uidx);
        }
        local_sample_idx_order[uii] = cur_idx;
      }
    }

    const uint32_t variant_ct = common->variant_ct;

    const GlmFlags glm_flags = glm_info_ptr->flags;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    // forced-singlethreaded
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto GlmLinear_ret_1;
    }
    const uint32_t report_neglog10p = (glm_flags / kfGlmLog10) & 1;
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t domdev_present_p1 = domdev_present + 1;

    const uint32_t constraint_ct = common->constraint_ct;
    const uint32_t constraint_ct_x = common->constraint_ct_x;
    const uint32_t constraint_ct_y = common->constraint_ct_y;

    const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
    uint32_t biallelic_predictor_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
    uint32_t biallelic_predictor_ct_x = 2 + covar_ct_x * (1 + add_interactions);
    uint32_t biallelic_predictor_ct_y = 2 + covar_ct_y * (1 + add_interactions);
    const uintptr_t* parameter_subset = common->parameter_subset;
    const uintptr_t* parameter_subset_x = common->parameter_subset_x;
    const uintptr_t* parameter_subset_y = common->parameter_subset_y;
    if (parameter_subset) {
      biallelic_predictor_ct = PopcountWords(parameter_subset, BitCtToWordCt(biallelic_predictor_ct));
      if (sample_ct_x) {
        biallelic_predictor_ct_x = PopcountWords(parameter_subset_x, BitCtToWordCt(biallelic_predictor_ct_x));
      } else {
        biallelic_predictor_ct_x = 0;
      }
      if (sample_ct_y) {
        // bugfix (7 Feb 2018): had biallelic_predictor_ct_x on right side
        // here, oops
        biallelic_predictor_ct_y = PopcountWords(parameter_subset_y, BitCtToWordCt(biallelic_predictor_ct_y));
      } else {
        biallelic_predictor_ct_y = 0;
      }
    }
    uint32_t biallelic_reported_test_ct = GetBiallelicReportedTestCt(parameter_subset, glm_flags, covar_ct, common->tests_flag);
    uintptr_t max_reported_test_ct = biallelic_reported_test_ct;
    uint32_t biallelic_reported_test_ct_x = 0;
    if (sample_ct_x) {
      biallelic_reported_test_ct_x = GetBiallelicReportedTestCt(parameter_subset_x, glm_flags, covar_ct_x, common->tests_flag);
      if (biallelic_reported_test_ct_x > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_x;
      }
    }
    uint32_t biallelic_reported_test_ct_y = 0;
    if (sample_ct_y) {
      biallelic_reported_test_ct_y = GetBiallelicReportedTestCt(parameter_subset_y, glm_flags, covar_ct_y, common->tests_flag);
      if (biallelic_reported_test_ct_y > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_y;
      }
    }
    const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
    const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
    const GlmColFlags glm_cols = glm_info_ptr->cols;
    const uint32_t test_col = glm_cols & kfGlmColTest;
    if (unlikely((!test_col) && (max_reported_test_ct > 1))) {
      // this is okay in plain multiallelic case due to A1 column
      logerrputs("Error: --glm's 'test' column cannot be omitted when results for multiple\npredictors are reported.  (Did you forget 'hide-covar'?)\n");
      goto GlmLinear_ret_INCONSISTENT_INPUT;
    }
    const uint32_t main_mutated = ((glm_flags & (kfGlmDominant | kfGlmRecessive | kfGlmHethom)) != kfGlm0);
    // bugfix (4 Mar 2019): forgot to update this for --tests
    const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!main_mutated) && (!common->tests_flag) && (!add_interactions);
    if (beta_se_multiallelic_fused || (!hide_covar)) {
      max_reported_test_ct += max_extra_allele_ct;
    }
    common->max_reported_test_ct = max_reported_test_ct;

    uint32_t x_code = UINT32_MAXM1;
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    if (sample_ct_x) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetX, &x_code, &x_start, &x_end);
    }
    uint32_t y_code = UINT32_MAXM1;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    if (sample_ct_y) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetY, &y_code, &y_start, &y_end);
    }
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t chr_col = glm_cols & kfGlmColChrom;

    // includes trailing tab
    char* chr_buf = nullptr;
    if (chr_col) {
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto GlmLinear_ret_NOMEM;
      }
    }

    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > variant_ct) {
      calc_thread_ct = variant_ct;
    }

    const uint32_t main_omitted = (parameter_subset && (!IsSet(parameter_subset, 1)));
    const uint32_t xmain_ct = main_mutated + main_omitted;
    uintptr_t workspace_alloc = GetLinearWorkspaceSize(sample_ct, biallelic_predictor_ct, max_extra_allele_ct, constraint_ct, xmain_ct);
    if (sample_ct_x) {
      const uintptr_t workspace_alloc_x = GetLinearWorkspaceSize(sample_ct_x, biallelic_predictor_ct_x, max_extra_allele_ct, constraint_ct_x, xmain_ct);
      if (workspace_alloc_x > workspace_alloc) {
        workspace_alloc = workspace_alloc_x;
      }
    }
    if (sample_ct_y) {
      const uintptr_t workspace_alloc_y = GetLinearWorkspaceSize(sample_ct_y, biallelic_predictor_ct_y, max_extra_allele_ct, constraint_ct_y, xmain_ct);
      if (workspace_alloc_y > workspace_alloc) {
        workspace_alloc = workspace_alloc_y;
      }
    }
    // +1 is for top-level common->workspace_bufs
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;

    uintptr_t per_variant_xalloc_byte_ct = max_sample_ct * local_covar_ct * sizeof(double);
    uintptr_t per_alt_allele_xalloc_byte_ct = sizeof(LinearAuxResult);
    if (beta_se_multiallelic_fused) {
      per_variant_xalloc_byte_ct += 2 * max_reported_test_ct * sizeof(double);
    } else {
      per_alt_allele_xalloc_byte_ct += 2 * max_reported_test_ct * sizeof(double);
    }
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    common->thread_mhc = nullptr;
    common->dosage_presents = nullptr;
    common->dosage_mains = nullptr;
    uint32_t read_block_size;
    uintptr_t max_alt_allele_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, max_sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, per_alt_allele_xalloc_byte_ct, pgfip, &calc_thread_ct, &common->genovecs, max_extra_allele_ct? (&common->thread_mhc) : nullptr, nullptr, nullptr, dosage_is_present? (&common->dosage_presents) : nullptr, dosage_is_present? (&common->dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, &max_alt_allele_block_size, main_loadbufs, &common->pgr_ptrs, &common->read_variant_uidx_starts))) {
      goto GlmLinear_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto GlmLinear_ret_NOMEM;
    }
    LinearAuxResult* linear_block_aux_bufs[2];
    double* block_beta_se_bufs[2];

    for (uint32_t uii = 0; uii != 2; ++uii) {
      if (unlikely(BIGSTACK_ALLOC_X(LinearAuxResult, max_alt_allele_block_size, &(linear_block_aux_bufs[uii])))) {
        // shouldn't be possible for these to fail?
        goto GlmLinear_ret_NOMEM;
      }
      if (beta_se_multiallelic_fused) {
        if (unlikely(bigstack_alloc_d(read_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii])))) {
          goto GlmLinear_ret_NOMEM;
        }
      } else {
        if (unlikely(bigstack_alloc_d(max_alt_allele_block_size * 2 * max_reported_test_ct, &(block_beta_se_bufs[uii])))) {
          goto GlmLinear_ret_NOMEM;
        }
      }
      if (local_covar_ct) {
        // bugfix (5 Mar 2018): don't want sizeof(double) here
        if (unlikely(bigstack_alloc_d(read_block_size * max_sample_ct * local_covar_ct, &(ctx->local_covars_vcmaj_d[uii])))) {
          goto GlmLinear_ret_NOMEM;
        }
      } else {
        ctx->local_covars_vcmaj_d[uii] = nullptr;
      }
    }

    common->workspace_bufs = S_CAST(unsigned char**, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t)));
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      common->workspace_bufs[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(workspace_alloc));
    }
    common->err_info = (~0LLU) << 32;
    SetThreadFuncAndData(GlmLinearThread, ctx, &tg);

    const uint32_t ref_col = glm_cols & kfGlmColRef;
    const uint32_t alt1_col = glm_cols & kfGlmColAlt1;
    const uint32_t alt_col = glm_cols & kfGlmColAlt;
    const uint32_t ax_col = glm_cols & kfGlmColAx;
    const uint32_t a1_ct_col = glm_cols & kfGlmColA1count;
    const uint32_t tot_allele_col = glm_cols & kfGlmColTotallele;
    const uint32_t a1_freq_col = glm_cols & kfGlmColA1freq;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t beta_col = glm_cols & (kfGlmColBeta | kfGlmColOrbeta);
    const uint32_t se_col = glm_cols & kfGlmColSe;
    const uint32_t ci_col = (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t t_col = glm_cols & kfGlmColTz;
    const uint32_t p_col = glm_cols & kfGlmColP;
    const uint32_t err_col = glm_cols & kfGlmColErr;
    *cswritep++ = '#';
    if (chr_col) {
      cswritep = strcpya_k(cswritep, "CHROM\t");
    }
    if (variant_bps) {
      cswritep = strcpya_k(cswritep, "POS\t");
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
    if (a1_ct_col) {
      cswritep = strcpya_k(cswritep, "\tA1_CT");
    }
    if (tot_allele_col) {
      cswritep = strcpya_k(cswritep, "\tALLELE_CT");
    }
    if (a1_freq_col) {
      cswritep = strcpya_k(cswritep, "\tA1_FREQ");
    }
    if (mach_r2_col) {
      cswritep = strcpya_k(cswritep, "\tMACH_R2");
    }
    if (test_col) {
      cswritep = strcpya_k(cswritep, "\tTEST");
    }
    if (nobs_col) {
      cswritep = strcpya_k(cswritep, "\tOBS_CT");
    }
    if (beta_col) {
      cswritep = strcpya_k(cswritep, "\tBETA");
    }
    if (se_col) {
      cswritep = strcpya_k(cswritep, "\tSE");
    }
    double ci_zt = 0.0;
    if (ci_col) {
      cswritep = strcpya_k(cswritep, "\tL");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      cswritep = strcpya_k(cswritep, "\tU");
      cswritep = dtoa_g(ci_size * 100, cswritep);
      ci_zt = QuantileToZscore((ci_size + 1.0) * 0.5);
    }
    if (t_col) {
      if (!constraint_ct) {
        cswritep = strcpya_k(cswritep, "\tT_STAT");
      } else {
        // F-statistic for joint tests.
        cswritep = strcpya_k(cswritep, "\tT_OR_F_STAT");
      }
    }
    if (p_col) {
      if (report_neglog10p) {
        cswritep = strcpya_k(cswritep, "\tLOG10_P");
      } else {
        cswritep = strcpya_k(cswritep, "\tP");
      }
    }
    if (err_col) {
      cswritep = strcpya_k(cswritep, "\tERRCODE");
    }
    AppendBinaryEoln(&cswritep);

    // Main workflow:
    // 1. Set n=0, load/skip block 0
    //
    // 2. Spawn threads processing block n
    // 3. If n>0, write results for block (n-1)
    // 4. Increment n by 1
    // 5. Load/skip block n unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8, Write results for last block
    uintptr_t write_variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t suppress_mach_r2 = 0;

    uint32_t cur_biallelic_reported_test_ct = 0;
    uint32_t primary_reported_test_idx = include_intercept;
    uint32_t cur_biallelic_predictor_ct = 0;
    uint32_t cur_constraint_ct = 0;

    const char* const* cur_test_names = nullptr;
    uint32_t prev_block_variant_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t allele_ct = 2;
    uint32_t omitted_allele_idx = 0;
    uintptr_t valid_allele_ct = 0;
    logprintfww5("--glm linear regression on phenotype '%s': ", cur_pheno_name);
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_variant_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto GlmLinear_ret_PGR_FAIL;
      }
      if (local_covar_ct && cur_block_variant_ct) {
        const uint32_t uidx_start = read_block_idx * read_block_size;
        const uint32_t uidx_end = MINV(raw_variant_ct, uidx_start + read_block_size);
        if (local_variant_include) {
          reterr = ReadLocalCovarBlock(common, local_sample_uidx_order, local_variant_include, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, local_covar_txsp, &local_line_idx, &local_xy, nullptr, ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
        } else {
          double* prev_local_covar_row_d = nullptr;
          if (variant_idx) {
            prev_local_covar_row_d = &(ctx->local_covars_vcmaj_d[1 - parity][S_CAST(uintptr_t, read_block_size - 1) * max_sample_ct * local_covar_ct]);
          }
          reterr = ReadRfmix2Block(common, variant_bps, local_sample_uidx_order, nullptr, prev_local_covar_row_d, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, glm_info_ptr->local_chrom_col, glm_info_ptr->local_bp_col, glm_info_ptr->local_first_covar_col, local_covar_txsp, &local_line_iter, &local_line_idx, &local_prev_chr_code, &local_chr_code, &local_bp, &local_skip_chr, nullptr, ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
          /*
          for (uint32_t uii = 0; uii < max_sample_ct; ++uii) {
            printf("%g ", ctx->local_covars_vcmaj_d[parity][uii]);
          }
          printf("\n");
          exit(1);
          */
        }
        if (unlikely(reterr)) {
          goto GlmLinear_ret_1;
        }
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, common->err_info);
        if (unlikely(reterr)) {
          goto GlmLinear_ret_PGR_FAIL;
        }
      }
      if (!IsLastBlock(&tg)) {
        common->cur_block_variant_ct = cur_block_variant_ct;
        const uint32_t uidx_start = read_block_idx * read_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_variant_ct, calc_thread_ct, uidx_start, common->read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, common->pgr_ptrs);
        ctx->block_aux = linear_block_aux_bufs[parity];
        common->block_beta_se = block_beta_se_bufs[parity];
        if (variant_idx + cur_block_variant_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto GlmLinear_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const double* beta_se_iter = block_beta_se_bufs[parity];
        const LinearAuxResult* cur_block_aux = linear_block_aux_bufs[parity];
        uintptr_t allele_bidx = 0;
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            if ((chr_idx == x_code) && sample_ct_x) {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct_x;
              cur_biallelic_predictor_ct = biallelic_predictor_ct_x;
              cur_constraint_ct = constraint_ct_x;
              cur_test_names = test_names_x;
            } else if ((chr_idx == y_code) && sample_ct_y) {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct_y;
              cur_biallelic_predictor_ct = biallelic_predictor_ct_y;
              cur_constraint_ct = constraint_ct_y;
              cur_test_names = test_names_y;
            } else {
              cur_biallelic_reported_test_ct = biallelic_reported_test_ct;
              cur_biallelic_predictor_ct = biallelic_predictor_ct;
              cur_constraint_ct = constraint_ct;
              cur_test_names = test_names;
            }
            suppress_mach_r2 = (chr_idx == x_code) || (chr_idx == mt_code);
            if (cur_constraint_ct) {
              primary_reported_test_idx = cur_biallelic_reported_test_ct - 1;
            }
            if (chr_col) {
              char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
              *chr_name_end = '\t';
              chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
            }
          }
          uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offsets[write_variant_uidx];
          }
          const uint32_t allele_ct_m1 = allele_ct - 1;
          const uint32_t extra_allele_ct = allele_ct - 2;
          if (omitted_alleles) {
            omitted_allele_idx = omitted_alleles[write_variant_uidx];
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          uint32_t variant_is_valid = 0;
          uint32_t a1_allele_idx = 0;
          for (uint32_t nonomitted_allele_idx = 0; nonomitted_allele_idx != allele_ct_m1; ++nonomitted_allele_idx, ++a1_allele_idx) {
            if (beta_se_multiallelic_fused) {
              if (!nonomitted_allele_idx) {
                primary_reported_test_idx = include_intercept;
              } else {
                primary_reported_test_idx = cur_biallelic_reported_test_ct + nonomitted_allele_idx - 1;
              }
            }
            if (nonomitted_allele_idx == omitted_allele_idx) {
              ++a1_allele_idx;
            }
            const double primary_beta = beta_se_iter[primary_reported_test_idx * 2];
            const double primary_se = beta_se_iter[primary_reported_test_idx * 2 + 1];
            const uint32_t allele_is_valid = (primary_se != -9.0);
            variant_is_valid |= allele_is_valid;
            {
              const LinearAuxResult* auxp = &(cur_block_aux[allele_bidx]);
              if (ln_pfilter <= 0.0) {
                if (!allele_is_valid) {
                  goto GlmLinear_allele_iterate;
                }
                double primary_ln_pval;
                if (!cur_constraint_ct) {
                  if (primary_beta == 0.0) {
                    primary_ln_pval = 0.0;
                  } else if (primary_se == 0.0) {
                    primary_ln_pval = -DBL_MAX;
                  } else {
                    const double primary_tstat = primary_beta / primary_se;
                    primary_ln_pval = TstatToLnP(primary_tstat, auxp->sample_obs_ct - cur_biallelic_predictor_ct - extra_allele_ct);
                  }
                } else {
                  primary_ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                }
                if (primary_ln_pval > ln_pfilter) {
                  if (orig_ln_pvals) {
                    orig_ln_pvals[valid_allele_ct] = primary_ln_pval;
                  }
                  goto GlmLinear_allele_iterate;
                }
              }
              uint32_t inner_reported_test_ct = cur_biallelic_reported_test_ct;
              if (extra_allele_ct) {
                if (beta_se_multiallelic_fused) {
                  if (!nonomitted_allele_idx) {
                    inner_reported_test_ct = 1 + include_intercept;
                  } else if (nonomitted_allele_idx == extra_allele_ct) {
                    inner_reported_test_ct -= include_intercept;
                  } else {
                    inner_reported_test_ct = 1;
                  }
                } else if (!hide_covar) {
                  inner_reported_test_ct += extra_allele_ct;
                }
              }
              // possible todo: make number-to-string operations, strlen(),
              // etc. happen only once per variant.
              for (uint32_t allele_test_idx = 0; allele_test_idx != inner_reported_test_ct; ++allele_test_idx) {
                uint32_t test_idx = allele_test_idx;
                if (beta_se_multiallelic_fused && nonomitted_allele_idx) {
                  if (!allele_test_idx) {
                    test_idx = primary_reported_test_idx;
                  } else {
                    // bugfix (26 Jun 2019): only correct to add 1 here in
                    // include_intercept case
                    test_idx += include_intercept;
                  }
                }
                if (chr_col) {
                  cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
                }
                if (variant_bps) {
                  cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
                }
                cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
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
                      goto GlmLinear_ret_WRITE_FAIL;
                    }
                    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                  }
                  --cswritep;
                }
                *cswritep++ = '\t';
                const uint32_t multi_a1 = extra_allele_ct && beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx);
                if (multi_a1) {
                  for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                    if (allele_idx == omitted_allele_idx) {
                      continue;
                    }
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLinear_ret_WRITE_FAIL;
                    }
                    cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                  }
                  --cswritep;
                } else {
                  cswritep = strcpya(cswritep, cur_alleles[a1_allele_idx]);
                }
                if (ax_col) {
                  *cswritep++ = '\t';
                  if (beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx)) {
                    if (unlikely(Cswrite(&css, &cswritep))) {
                      goto GlmLinear_ret_WRITE_FAIL;
                    }
                    cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
                  } else {
                    for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                      if (allele_idx == a1_allele_idx) {
                        continue;
                      }
                      if (unlikely(Cswrite(&css, &cswritep))) {
                        goto GlmLinear_ret_WRITE_FAIL;
                      }
                      cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                    }
                    --cswritep;
                  }
                }
                if (a1_ct_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_dosage, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (tot_allele_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa(auxp->allele_obs_ct, cswritep);
                }
                if (a1_freq_col) {
                  *cswritep++ = '\t';
                  if (!multi_a1) {
                    cswritep = dtoa_g(auxp->a1_dosage / S_CAST(double, auxp->allele_obs_ct), cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (mach_r2_col) {
                  *cswritep++ = '\t';
                  if (!suppress_mach_r2) {
                    cswritep = dtoa_g(auxp->mach_r2, cswritep);
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (test_col) {
                  *cswritep++ = '\t';
                  if (test_idx < cur_biallelic_reported_test_ct) {
                    cswritep = strcpya(cswritep, cur_test_names[test_idx]);
                  } else {
                    // always use basic dosage for untested alleles
                    cswritep = strcpya_k(cswritep, "ADD");
                    if (!beta_se_multiallelic_fused) {
                      // extra alt allele covariate.
                      uint32_t test_xallele_idx = test_idx - cur_biallelic_reported_test_ct;
                      if (omitted_allele_idx < a1_allele_idx) {
                        test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                      }
                      test_xallele_idx = test_xallele_idx + (test_xallele_idx >= a1_allele_idx);
                      if (a1_allele_idx < omitted_allele_idx) {
                        test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                      }
                      if (!test_xallele_idx) {
                        cswritep = strcpya_k(cswritep, "_REF");
                      } else {
                        cswritep = strcpya_k(cswritep, "_ALT");
                        cswritep = u32toa(test_xallele_idx, cswritep);
                      }
                    }
                  }
                }
                if (nobs_col) {
                  *cswritep++ = '\t';
                  cswritep = u32toa(auxp->sample_obs_ct, cswritep);
                }
                double ln_pval = kLnPvalError;
                double tstat = 0.0;
                uint32_t test_is_valid;
                if ((!cur_constraint_ct) || (test_idx != primary_reported_test_idx)) {
                  double beta = beta_se_iter[2 * test_idx];
                  double se = beta_se_iter[2 * test_idx + 1];
                  test_is_valid = (se != -9.0);
                  if (test_is_valid) {
                    tstat = beta / se;
                    ln_pval = TstatToLnP(tstat, auxp->sample_obs_ct - cur_biallelic_predictor_ct - extra_allele_ct);
                  }
                  if (beta_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(beta, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  if (se_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(se, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  if (ci_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      const double ci_halfwidth = ci_zt * se;
                      cswritep = dtoa_g(beta - ci_halfwidth, cswritep);
                      *cswritep++ = '\t';
                      cswritep = dtoa_g(beta + ci_halfwidth, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA\tNA");
                    }
                  }
                  if (t_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(tstat, cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                } else {
                  // joint test
                  test_is_valid = allele_is_valid;
                  if (beta_col) {
                    cswritep = strcpya_k(cswritep, "\tNA");
                  }
                  if (se_col) {
                    cswritep = strcpya_k(cswritep, "\tNA");
                  }
                  if (ci_col) {
                    cswritep = strcpya_k(cswritep, "\tNA\tNA");
                  }
                  if (t_col) {
                    *cswritep++ = '\t';
                    if (test_is_valid) {
                      cswritep = dtoa_g(primary_se / u31tod(cur_constraint_ct), cswritep);
                    } else {
                      cswritep = strcpya_k(cswritep, "NA");
                    }
                  }
                  // could avoid recomputing
                  if (test_is_valid) {
                    ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                  }
                }
                if (p_col) {
                  *cswritep++ = '\t';
                  if (test_is_valid) {
                    if (report_neglog10p) {
                      const double reported_val = (-kRecipLn10) * ln_pval;
                      cswritep = dtoa_g(reported_val, cswritep);
                    } else {
                      const double reported_ln = MAXV(ln_pval, output_min_ln);
                      cswritep = lntoa_g(reported_ln, cswritep);
                    }
                  } else {
                    cswritep = strcpya_k(cswritep, "NA");
                  }
                }
                if (err_col) {
                  *cswritep++ = '\t';
                  if (test_is_valid) {
                    *cswritep++ = '.';
                  } else {
                    uint64_t glm_errcode;
                    memcpy(&glm_errcode, &(beta_se_iter[2 * test_idx]), 8);
                    cswritep = AppendGlmErrstr(glm_errcode, cswritep);
                  }
                }
                AppendBinaryEoln(&cswritep);
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto GlmLinear_ret_WRITE_FAIL;
                }
                if ((test_idx == primary_reported_test_idx) && allele_is_valid) {
                  if (orig_ln_pvals) {
                    orig_ln_pvals[valid_allele_ct] = ln_pval;
                  }
                }
              }
            }
          GlmLinear_allele_iterate:
            ++allele_bidx;
            valid_allele_ct += allele_is_valid;
            if (valid_alleles && allele_is_valid) {
              SetBit(allele_idx_offset_base + a1_allele_idx, valid_alleles);
            }
            if (!beta_se_multiallelic_fused) {
              beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
            }
          }
          if (beta_se_multiallelic_fused) {
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
          if ((!variant_is_valid) && valid_alleles) {
            ClearBit(write_variant_uidx, valid_variants);
          }
        }
      }
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
      prev_block_variant_ct = cur_block_variant_ct;
      variant_idx += cur_block_variant_ct;
      // crucially, this is independent of the PgenReader block_base
      // pointers
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto GlmLinear_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    logprintf("Results written to %s .\n", outname);
    *valid_allele_ct_ptr = valid_allele_ct;
  }
  while (0) {
  GlmLinear_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GlmLinear_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
    break;
  GlmLinear_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  GlmLinear_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GlmLinear_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  GlmLinear_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 GlmLinear_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

uintptr_t GetLinearSubbatchWorkspaceSize(uint32_t sample_ct, uint32_t subbatch_size, uint32_t biallelic_predictor_ct, uint32_t max_extra_allele_ct, uint32_t constraint_ct, uint32_t xmain_ct) {
  // sample_ct * max_predictor_ct < 2^31, sample_ct * subbatch_size < 2^31,
  // subbatch_size <= 240, and max_predictor_ct < sqrt(2^31), so no overflows

  // sample_nm, tmp_nm = sample_ctl words
  uintptr_t workspace_size = 2 * RoundUpPow2(BitCtToWordCt(sample_ct) * sizeof(intptr_t), kCacheline);

  // nm_pheno_buf = sample_ct * subbatch_size doubles
  workspace_size += RoundUpPow2(sample_ct * subbatch_size * sizeof(double), kCacheline);

  const uint32_t max_predictor_ct = biallelic_predictor_ct + max_extra_allele_ct;
  // predictors_pmaj = (max_predictor_ct + main_mutated + main_omitted) *
  //                   sample_ct doubles
  workspace_size += RoundUpPow2((max_predictor_ct + xmain_ct) * sample_ct * sizeof(double), kCacheline);

  // xtx_inv, xtx_inv2 = max_predictor_ct * max_predictor_ct doubles
  workspace_size += 2 * RoundUpPow2(max_predictor_ct * max_predictor_ct * sizeof(double), kCacheline);

  // dbl_2d_buf = max_predictor_ct * max(max_predictor_ct, 7) doubles
  workspace_size += RoundUpPow2(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), kCacheline);

  // inverse_corr_buf = (max_predictor_ct - 1) * max(max_predictor_ct - 1, 4) doubles
  workspace_size += RoundUpPow2((max_predictor_ct - 1) * MAXV((max_predictor_ct - 1), 4) * sizeof(double), kCacheline);

  // semicomputed_biallelic_corr_matrix = (max_predictor_ct - 1)^2 doubles
  workspace_size += RoundUpPow2((biallelic_predictor_ct - 1) * (biallelic_predictor_ct - 1) * sizeof(double), kCacheline);

  // semicomputed_biallelic_inv_corr_sqrts = biallelic_predictor_ct doubles
  workspace_size += RoundUpPow2(biallelic_predictor_ct * sizeof(double), kCacheline);

  // fitted_coefs, xt_y = max_predictor_ct * subbatch_size doubles
  workspace_size += 2 * RoundUpPow2(max_predictor_ct * subbatch_size * sizeof(double), kCacheline);

  // inv_1d_buf
  workspace_size += RoundUpPow2(MAXV(max_predictor_ct, constraint_ct) * kMatrixInvertBuf1CheckedAlloc, kCacheline);

  // machr2_dosage_sums, machr2_dosage_ssqs
  workspace_size += RoundUpPow2((2 + max_extra_allele_ct) * sizeof(uint64_t) * 2, kCacheline);

  // pheno_ssq_bases, geno_pheno_prods, domdev_pheno_prods: subbatch_size
  workspace_size += 3 * RoundUpPow2(subbatch_size * sizeof(double), kCacheline);

  if (constraint_ct) {
    // tmphxs_buf, h_transpose_buf = constraint_ct * max_predictor_ct doubles
    workspace_size += 2 * RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(double), kCacheline);

    // inner_buf = constraint_ct * constraint_ct
    workspace_size += RoundUpPow2(constraint_ct * constraint_ct * sizeof(double), kCacheline);

    // cur_constraints_con_major = constraint_ct * max_predictor_ct doubles
    workspace_size += RoundUpPow2(constraint_ct * max_predictor_ct * sizeof(double), kCacheline);
  }
  return workspace_size;
}

THREAD_FUNC_DECL GlmLinearSubbatchThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  GlmLinearCtx* ctx = S_CAST(GlmLinearCtx*, arg->sharedp->context);
  GlmCtx* common = ctx->common;

  PgenReader* pgrp = common->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = common->genovecs[tidx];
  pgv.dosage_present = nullptr;
  pgv.dosage_main = nullptr;
  if (common->dosage_presents) {
    pgv.dosage_present = common->dosage_presents[tidx];
    pgv.dosage_main = common->dosage_mains[tidx];
  }
  unsigned char* workspace_buf = common->workspace_bufs[tidx];
  const uintptr_t* variant_include = common->variant_include;
  const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
  const AlleleCode* omitted_alleles = common->omitted_alleles;
  const uintptr_t* sex_male_collapsed = common->sex_male_collapsed;
  const ChrInfo* cip = common->cip;
  const uint32_t* subset_chr_fo_vidx_start = common->subset_chr_fo_vidx_start;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const GlmFlags glm_flags = common->glm_flags;
  const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
  const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
  const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t joint_genotypic = (glm_flags / kfGlmGenotypic) & 1;
  const uint32_t joint_hethom = (glm_flags / kfGlmHethom) & 1;
  const uint32_t domdev_present = joint_genotypic || joint_hethom;
  const uint32_t domdev_present_p1 = domdev_present + 1;
  const uint32_t reported_pred_uidx_start = 1 - include_intercept;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = common->is_xchr_model_1;
  const double max_corr = common->max_corr;
  const double vif_thresh = common->vif_thresh;
  const uintptr_t max_reported_test_ct = common->max_reported_test_ct;
  const uintptr_t local_covar_ct = common->local_covar_ct;
  const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
  const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!model_dominant) && (!model_recessive) && (!common->tests_flag) && (!add_interactions);
  const uint32_t subbatch_size = ctx->subbatch_size;
  uintptr_t max_sample_ct = MAXV(common->sample_ct, common->sample_ct_x);
  if (max_sample_ct < common->sample_ct_y) {
    max_sample_ct = common->sample_ct_y;
  }
  pgv.patch_01_set = nullptr;
  pgv.patch_01_vals = nullptr;
  pgv.patch_10_set = nullptr;
  pgv.patch_10_vals = nullptr;
  if (common->thread_mhc) {
    const uint32_t max_sample_ctl = BitCtToWordCt(max_sample_ct);
    pgv.patch_01_set = common->thread_mhc[tidx];
    pgv.patch_01_vals = R_CAST(AlleleCode*, &(pgv.patch_01_set[max_sample_ctl]));
    AlleleCode* patch_01_vals_end = &(pgv.patch_01_vals[max_sample_ct]);
    VecAlignUp(&patch_01_vals_end);
    pgv.patch_10_set = R_CAST(uintptr_t*, patch_01_vals_end);
    pgv.patch_10_vals = R_CAST(AlleleCode*, &(pgv.patch_10_set[max_sample_ctl]));
  }
  pgv.patch_01_ct = 0;
  pgv.patch_10_ct = 0;
  pgv.multidosage_sample_ct = 0;
  uint32_t variant_idx_offset = 0;
  uint32_t allele_ct = 2;
  uint32_t omitted_allele_idx = 0;
  uint32_t extra_regression_ct = 0;
  double main_dosage_sum = 0.0;
  double main_dosage_ssq = 0.0;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_variant_ct = common->cur_block_variant_ct;
    uint32_t variant_bidx = (tidx * cur_block_variant_ct) / calc_thread_ct;
    const uint32_t variant_bidx_end = ((tidx + 1) * cur_block_variant_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, common->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);

    uintptr_t allele_bidx = variant_bidx;
    if (max_extra_allele_ct) {
      allele_bidx = variant_bidx + CountExtraAlleles(variant_include, allele_idx_offsets, common->read_variant_uidx_starts[0], common->read_variant_uidx_starts[tidx], 0);
    }
    double* beta_se_iter = common->block_beta_se;
    if (beta_se_multiallelic_fused) {
      beta_se_iter = &(beta_se_iter[(2 * k1LU * subbatch_size) * max_reported_test_ct * variant_bidx]);
    } else {
      beta_se_iter = &(beta_se_iter[(2 * k1LU * subbatch_size) * max_reported_test_ct * allele_bidx]);
    }

    LinearAuxResult* block_aux_iter = &(ctx->block_aux[allele_bidx]);
    const double* local_covars_iter = nullptr;
    if (local_covar_ct) {
      // &(nullptr[0]) is okay in C++, but undefined in C
      local_covars_iter = &(ctx->local_covars_vcmaj_d[parity][variant_bidx * max_sample_ct * local_covar_ct]);
    }
    while (variant_bidx < variant_bidx_end) {
      const uint32_t variant_idx = variant_bidx + variant_idx_offset;
      const uint32_t chr_fo_idx = CountSortedSmallerU32(&(subset_chr_fo_vidx_start[1]), cip->chr_ct, variant_idx + 1);
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      uint32_t cur_variant_bidx_end = subset_chr_fo_vidx_start[chr_fo_idx + 1] - variant_idx_offset;
      if (cur_variant_bidx_end > variant_bidx_end) {
        cur_variant_bidx_end = variant_bidx_end;
      }
      const uint32_t is_x = (chr_idx == x_code);
      const uint32_t is_y = (chr_idx == y_code);
      const uint32_t is_nonx_haploid = (!is_x) && IsSet(cip->haploid_mask, chr_idx);
      const uintptr_t* cur_sample_include;
      const uint32_t* cur_sample_include_cumulative_popcounts;
      const double* cur_pheno_pmaj;
      const RegressionNmPrecomp* nm_precomp;
      const double* cur_covars_cmaj;
      const uintptr_t* cur_parameter_subset;
      const uintptr_t* cur_joint_test_params;
      uint32_t cur_sample_ct;
      uint32_t cur_covar_ct;
      uint32_t cur_constraint_ct;
      if (is_y && common->sample_include_y) {
        cur_sample_include = common->sample_include_y;
        cur_sample_include_cumulative_popcounts = common->sample_include_y_cumulative_popcounts;
        cur_pheno_pmaj = ctx->pheno_y_d;
        nm_precomp = common->nm_precomp_y;
        cur_covars_cmaj = ctx->covars_cmaj_y_d;
        cur_parameter_subset = common->parameter_subset_y;
        cur_joint_test_params = common->joint_test_params_y;
        cur_sample_ct = common->sample_ct_y;
        cur_covar_ct = common->covar_ct_y;
        cur_constraint_ct = common->constraint_ct_y;
      } else if (is_x && common->sample_include_x) {
        cur_sample_include = common->sample_include_x;
        cur_sample_include_cumulative_popcounts = common->sample_include_x_cumulative_popcounts;
        cur_pheno_pmaj = ctx->pheno_x_d;
        nm_precomp = common->nm_precomp_x;
        cur_covars_cmaj = ctx->covars_cmaj_x_d;
        cur_parameter_subset = common->parameter_subset_x;
        cur_joint_test_params = common->joint_test_params_x;
        cur_sample_ct = common->sample_ct_x;
        cur_covar_ct = common->covar_ct_x;
        cur_constraint_ct = common->constraint_ct_x;
      } else {
        cur_sample_include = common->sample_include;
        cur_sample_include_cumulative_popcounts = common->sample_include_cumulative_popcounts;
        cur_pheno_pmaj = ctx->pheno_d;
        nm_precomp = common->nm_precomp;
        cur_covars_cmaj = ctx->covars_cmaj_d;
        cur_parameter_subset = common->parameter_subset;
        cur_joint_test_params = common->joint_test_params;
        cur_sample_ct = common->sample_ct;
        cur_covar_ct = common->covar_ct;
        cur_constraint_ct = common->constraint_ct;
      }
      const uint32_t sample_ctl = BitCtToWordCt(cur_sample_ct);
      const uint32_t sample_ctl2 = NypCtToWordCt(cur_sample_ct);
      const uint32_t cur_biallelic_predictor_ct_base = 2 + domdev_present + cur_covar_ct * (1 + add_interactions * domdev_present_p1);
      uint32_t cur_biallelic_predictor_ct = cur_biallelic_predictor_ct_base;
      uint32_t literal_covar_ct = cur_covar_ct;
      if (cur_parameter_subset) {
        cur_biallelic_predictor_ct = PopcountWords(cur_parameter_subset, BitCtToWordCt(cur_biallelic_predictor_ct_base));
        literal_covar_ct = PopcountBitRange(cur_parameter_subset, 2 + domdev_present, 2 + domdev_present + cur_covar_ct);
      }
      const uint32_t max_predictor_ct = cur_biallelic_predictor_ct + max_extra_allele_ct;
      uint32_t reported_pred_uidx_biallelic_end;
      if (hide_covar) {
        if (!cur_parameter_subset) {
          reported_pred_uidx_biallelic_end = 2 + domdev_present;
        } else {
          reported_pred_uidx_biallelic_end = IsSet(cur_parameter_subset, 1) + domdev_present_p1;
        }
      } else {
        reported_pred_uidx_biallelic_end = cur_biallelic_predictor_ct;
      }
      // nm_predictors_pmaj_buf may require up to two extra columns omitted
      // from the main regression.
      // 1. In the multiallelic dominant/recessive/hethom cases, the original
      //    genotype column does not appear in the regression, and we'd rather
      //    not reconstruct it from genovec, etc. when we need to swap it out
      //    for another allele, so we keep the original genotype in an extra
      //    column.
      //    To reduce code bloat, we now handle the biallelic cases in the same
      //    way; this is one of the more peripheral code paths so adding more
      //    complexity to speed it up is less justifiable.
      // 2. If --parameters excludes the main (possibly dominant/recessive)
      //    genotype column but does care about an interaction, we want a copy
      //    of what the main genotype column's contents would have been to
      //    refer to.
      const uint32_t main_omitted = cur_parameter_subset && (!IsSet(cur_parameter_subset, 1));
      const uint32_t main_mutated = model_dominant || model_recessive || joint_hethom;
      unsigned char* workspace_iter = workspace_buf;
      uintptr_t* sample_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      uintptr_t* tmp_nm = S_CAST(uintptr_t*, arena_alloc_raw_rd(sample_ctl * sizeof(intptr_t), &workspace_iter));
      double* nm_pheno_buf = S_CAST(double*, arena_alloc_raw_rd(cur_sample_ct * subbatch_size * sizeof(double), &workspace_iter));
      double* nm_predictors_pmaj_buf = S_CAST(double*, arena_alloc_raw_rd((max_predictor_ct + main_mutated + main_omitted) * cur_sample_ct * sizeof(double), &workspace_iter));
      double* xtx_inv = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ct * sizeof(double), &workspace_iter));
      double* fitted_coefs = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * subbatch_size * sizeof(double), &workspace_iter));
      double* xt_y = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * subbatch_size * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_corr_matrix = S_CAST(double*, arena_alloc_raw_rd((cur_biallelic_predictor_ct - 1) * (cur_biallelic_predictor_ct - 1) * sizeof(double), &workspace_iter));
      double* semicomputed_biallelic_inv_corr_sqrts = S_CAST(double*, arena_alloc_raw_rd(cur_biallelic_predictor_ct * sizeof(double), &workspace_iter));
      MatrixInvertBuf1* inv_1d_buf = S_CAST(MatrixInvertBuf1*, arena_alloc_raw_rd(MAXV(max_predictor_ct, cur_constraint_ct) * kMatrixInvertBuf1CheckedAlloc, &workspace_iter));
      double* dbl_2d_buf = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * MAXV(max_predictor_ct, 7) * sizeof(double), &workspace_iter));
      uint64_t* machr2_dosage_sums = S_CAST(uint64_t*, arena_alloc_raw_rd((max_extra_allele_ct + 2) * sizeof(uint64_t) * 2, &workspace_iter));
      uint64_t* machr2_dosage_ssqs = &(machr2_dosage_sums[max_extra_allele_ct + 2]);
      double* pheno_ssq_bases = S_CAST(double*, arena_alloc_raw_rd(subbatch_size * sizeof(double), &workspace_iter));
      double* geno_pheno_prods = S_CAST(double*, arena_alloc_raw_rd(subbatch_size * sizeof(double), &workspace_iter));
      double* domdev_pheno_prods = S_CAST(double*, arena_alloc_raw_rd(subbatch_size * sizeof(double), &workspace_iter));
      double* xtx_inv2 = S_CAST(double*, arena_alloc_raw_rd(max_predictor_ct * max_predictor_ct * sizeof(double), &workspace_iter));

      // could technically have this overlap fitted_coefs/xt_y, but that sets
      // the stage for future bugs
      double* inverse_corr_buf = S_CAST(double*, arena_alloc_raw_rd((max_predictor_ct - 1) * MAXV((max_predictor_ct - 1), 4) * sizeof(double), &workspace_iter));

      // joint test only
      double* tmphxs_buf = nullptr;
      double* h_transpose_buf = nullptr;
      double* inner_buf = nullptr;
      double* cur_constraints_con_major = nullptr;
      if (cur_constraint_ct) {
        tmphxs_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        h_transpose_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        inner_buf = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * cur_constraint_ct * sizeof(double), &workspace_iter));
        cur_constraints_con_major = S_CAST(double*, arena_alloc_raw_rd(cur_constraint_ct * max_predictor_ct * sizeof(double), &workspace_iter));
        ZeroDArr(cur_constraint_ct * max_predictor_ct, cur_constraints_con_major);
        const uint32_t first_joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
        cur_constraints_con_major[first_joint_test_idx] = 1.0;
        // Rest of this matrix must be updated later, since cur_predictor_ct
        // changes at multiallelic variants.
      }
      assert(S_CAST(uintptr_t, workspace_iter - workspace_buf) == GetLinearSubbatchWorkspaceSize(cur_sample_ct, subbatch_size, cur_biallelic_predictor_ct, max_extra_allele_ct, cur_constraint_ct, main_mutated + main_omitted));
      for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
        const double* cur_pheno = &(cur_pheno_pmaj[pheno_idx * cur_sample_ct]);
        pheno_ssq_bases[pheno_idx] = DotprodD(cur_pheno, cur_pheno, cur_sample_ct);
      }
      const double cur_sample_ct_recip = 1.0 / u31tod(cur_sample_ct);
      const double cur_sample_ct_m1_recip = 1.0 / u31tod(cur_sample_ct - 1);
      const uint32_t sparse_optimization_eligible = (!is_x) && nm_precomp;
      double geno_d_lookup[2];
      if (sparse_optimization_eligible) {
        geno_d_lookup[1] = 1.0;
        if (is_nonx_haploid) {
          geno_d_lookup[0] = 0.5;
        } else if (model_recessive || joint_hethom) {
          geno_d_lookup[0] = 0.0;
        } else {
          geno_d_lookup[0] = 1.0;
          if (!model_dominant) {
            geno_d_lookup[1] = 2.0;
          }
        }
      }
      const double* xtx_image = nullptr;
      const double* covarx_dotprod_inv = nullptr;
      const double* corr_inv = nullptr;
      const double* xt_y_image = nullptr;
      if (nm_precomp) {
        xtx_image = nm_precomp->xtx_image;
        covarx_dotprod_inv = nm_precomp->covarx_dotprod_inv;
        corr_inv = nm_precomp->corr_inv;
        const uintptr_t nongeno_pred_ct = cur_biallelic_predictor_ct - domdev_present - 2;
        const uintptr_t nonintercept_biallelic_pred_ct = cur_biallelic_predictor_ct - 1;
        memcpy(semicomputed_biallelic_corr_matrix, nm_precomp->corr_image, nonintercept_biallelic_pred_ct * nonintercept_biallelic_pred_ct * sizeof(double));
        memcpy(&(semicomputed_biallelic_inv_corr_sqrts[domdev_present_p1]), nm_precomp->corr_inv_sqrts, nongeno_pred_ct * sizeof(double));
        xt_y_image = nm_precomp->xt_y_image;
      }
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(cur_sample_include_cumulative_popcounts, pgrp, &pssi);
      // when this is set, the last fully-processed variant had no missing
      // genotypes, and if the current variant also has no missing genotypes we
      // may be able to skip reinitialization of most of
      // nm_predictors_pmaj_buf.
      uint32_t prev_nm = 0;

      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      for (; variant_bidx != cur_variant_bidx_end; ++variant_bidx) {
        const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        if (allele_idx_offsets) {
          allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
          if (!beta_se_multiallelic_fused) {
            extra_regression_ct = allele_ct - 2;
          }
        }
        const uint32_t allele_ct_m2 = allele_ct - 2;
        const uint32_t expected_predictor_ct = cur_biallelic_predictor_ct + allele_ct_m2;
        PglErr reterr;
        if (!allele_ct_m2) {
          reterr = PgrGetD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
        } else {
          reterr = PgrGetMD(cur_sample_include, pssi, cur_sample_ct, variant_uidx, pgrp, &pgv);
          // todo: proper multiallelic dosage support
        }
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto GlmLinearSubbatchThread_err;
        }
        ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
        GenoarrCountFreqsUnsafe(pgv.genovec, cur_sample_ct, genocounts);
        uint32_t missing_ct = genocounts[3];
        if (!missing_ct) {
          SetAllBits(cur_sample_ct, sample_nm);
        } else {
          GenoarrToNonmissing(pgv.genovec, cur_sample_ct, sample_nm);
          if (pgv.dosage_ct) {
            BitvecOr(pgv.dosage_present, sample_ctl, sample_nm);
            missing_ct = cur_sample_ct - PopcountWords(sample_nm, sample_ctl);
          }
        }
        if (omitted_alleles) {
          omitted_allele_idx = omitted_alleles[variant_uidx];
        }
        uintptr_t const_alleles[DivUp(kPglMaxAltAlleleCt + 1, kBitsPerWord)];
        const uint32_t allele_ctl = DivUp(allele_ct, kBitsPerWord);
        ZeroWArr(allele_ctl, const_alleles);
        const uint32_t nm_sample_ct = cur_sample_ct - missing_ct;
        const uint32_t nm_sample_ctl = BitCtToWordCt(nm_sample_ct);
        // first predictor column: intercept
        if (!prev_nm) {
          for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
            nm_predictors_pmaj_buf[sample_idx] = 1.0;
          }
        }
        // second predictor column: main genotype
        double* genotype_vals = &(nm_predictors_pmaj_buf[nm_sample_ct]);
        if (main_mutated || main_omitted) {
          genotype_vals = &(nm_predictors_pmaj_buf[expected_predictor_ct * nm_sample_ct]);
        }
        uint32_t sparse_optimization = 0;
        double* multi_start = nullptr;
        if (!allele_ct_m2) {
          // When prev_nm is set and missing_ct is zero, we don't need to call
          // MultiplySelfTranspose() on nm_predictors_pmaj_buf to get all the
          // predictor x predictor dot products; instead we patch in the
          // genotype x predictor dot products that may change, copying the
          // rest from xtx_image.
          // As a side effect, it is no longer strictly necessary to fill the
          // genotype row of nm_predictors_pmaj_buf.  sparse_optimization
          // indicates that plink 1.9's QT --assoc sparse dot product algorithm
          // will be used instead.
          // probable todos: allow a few dosages to be present, cover chrX
          // case.
          if (omitted_allele_idx) {
            GenovecInvertUnsafe(cur_sample_ct, pgv.genovec);
            ZeroTrailingNyps(cur_sample_ct, pgv.genovec);
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
            }
            const uint32_t uii = genocounts[0];
            genocounts[0] = genocounts[2];
            genocounts[2] = uii;
          }
          uint64_t dosage_sum = (genocounts[1] + 2 * genocounts[2]) * 0x4000LLU;
          uint64_t dosage_ssq = (genocounts[1] + 4LLU * genocounts[2]) * 0x10000000LLU;
          if (!missing_ct) {
            // originally had genocounts[0] > 0.875 * nm_sample_ct threshold,
            // but then tried this on high-MAF data and it was still
            // substantially faster
            sparse_optimization = sparse_optimization_eligible && (!pgv.dosage_ct) && prev_nm;
            if (!sparse_optimization) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, genotype_vals);
              if (pgv.dosage_ct) {
                uintptr_t sample_idx_base = 0;
                uintptr_t dosage_present_bits = pgv.dosage_present[0];
                for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
                  const uintptr_t sample_idx = BitIter1(pgv.dosage_present, &sample_idx_base, &dosage_present_bits);
                  const uintptr_t dosage_val = pgv.dosage_main[dosage_idx];
                  // 32768 -> 2, 16384 -> 1, 0 -> 0
                  genotype_vals[sample_idx] = kRecipDosageMid * swtod(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_idx);
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                }
              }
            }
          } else {
            if (!pgv.dosage_ct) {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, genotype_vals);
            } else {
              uintptr_t sample_midx_base = 0;
              uintptr_t sample_nm_bits = sample_nm[0];
              uint32_t dosage_idx = 0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                const uintptr_t cur_geno = GetNyparrEntry(pgv.genovec, sample_midx);
                double cur_val;
                if (IsSet(pgv.dosage_present, sample_midx)) {
                  const uintptr_t dosage_val = pgv.dosage_main[dosage_idx++];
                  cur_val = kRecipDosageMid * swtod(dosage_val);
                  dosage_sum += dosage_val;
                  dosage_ssq += dosage_val * dosage_val;
                  if (cur_geno && (cur_geno != 3)) {
                    const uintptr_t prev_val = cur_geno * kDosageMid;
                    dosage_sum -= prev_val;
                    dosage_ssq -= prev_val * prev_val;
                  }
                } else {
                  cur_val = kSmallDoubles[cur_geno];
                }
                genotype_vals[sample_idx] = cur_val;
              }
            }
          }
          // Check for constant genotype column.
          // (Technically, we should recheck later in the chrX no-sex-covariate
          // --xchr-model 1 corner case.)
          if (!pgv.dosage_ct) {
            if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          } else if (pgv.dosage_ct == nm_sample_ct) {
            if (DosageIsConstant(dosage_sum, dosage_ssq, nm_sample_ct)) {
              const_alleles[0] = 3;
            }
          }
          machr2_dosage_sums[1 - omitted_allele_idx] = dosage_sum;
          machr2_dosage_ssqs[1 - omitted_allele_idx] = dosage_ssq;
          machr2_dosage_sums[omitted_allele_idx] = kDosageMax * S_CAST(uint64_t, nm_sample_ct) - dosage_sum;
          machr2_dosage_ssqs[omitted_allele_idx] = kDosageMax * (kDosageMax * S_CAST(uint64_t, nm_sample_ct) - 2 * dosage_sum) + dosage_ssq;
        } else {
          // multiallelic.
          // If some but not all alleles have constant dosages, we remove just
          // those alleles from the regressions; trim-alts is not necessary to
          // see what's going on with the other alleles.  To reduce parsing
          // complexity, the number of output lines is not affected by this;
          // the ones corresponding to the constant alleles have NA values.
          // Punt on sparse_optimization for now; may be worth revisiting after
          // multiallelic dosage implemented.
          // dosage_ct == 0 temporarily guaranteed if we reach here.
          assert(!pgv.dosage_ct);
          multi_start = &(nm_predictors_pmaj_buf[(expected_predictor_ct - allele_ct_m2) * nm_sample_ct]);
          ZeroU64Arr(allele_ct, machr2_dosage_sums);
          ZeroU64Arr(allele_ct, machr2_dosage_ssqs);
          // postpone multiply for now, since no multiallelic dosages
          machr2_dosage_sums[0] = genocounts[1] + 2 * genocounts[0];
          machr2_dosage_ssqs[0] = genocounts[1] + 4LLU * genocounts[0];
          if ((genocounts[0] == nm_sample_ct) || (genocounts[1] == nm_sample_ct) || (genocounts[2] == nm_sample_ct)) {
            SetBit(0, const_alleles);
          }
          if (omitted_allele_idx) {
            // Main genotype column starts as REF.
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallInvDoublePairs, nm_sample_ct, genotype_vals);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallInvDoubles, cur_sample_ct, genotype_vals);
            }
          }
          uint32_t rare_allele_ct = allele_ct_m2;
          double* alt1_start = nullptr;
          double* rarealt_start = multi_start;
          if (omitted_allele_idx != 1) {
            if (omitted_allele_idx) {
              alt1_start = multi_start;
              rarealt_start = &(rarealt_start[nm_sample_ct]);
              --rare_allele_ct;
            } else {
              alt1_start = genotype_vals;
            }
            if (!missing_ct) {
              GenoarrLookup16x8bx2(pgv.genovec, kSmallDoublePairs, nm_sample_ct, alt1_start);
            } else {
              GenoarrToDoublesRemoveMissing(pgv.genovec, kSmallDoubles, cur_sample_ct, alt1_start);
            }
          }
          ZeroDArr(rare_allele_ct * nm_sample_ct, rarealt_start);
          // Use sums as ones[] and ssqs as twos[] for rarealts; transform to
          // actual sums/ssqs later.
          if (pgv.patch_01_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_01_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_01_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ct + sample_idx] = 1.0;
                alt1_start[sample_idx] = 0.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                rarealt_start[(allele_code - 2) * nm_sample_ct + sample_idx] = 1.0;
                machr2_dosage_sums[allele_code] += 1;
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                alt1_start[sample_idx] = 0.0;
                const uint32_t allele_code = pgv.patch_01_vals[uii];
                machr2_dosage_sums[allele_code] += 1;
                if (allele_code == omitted_allele_idx) {
                  continue;
                }
                const uint32_t cur_col = allele_code - 2 - (allele_code > omitted_allele_idx);
                rarealt_start[cur_col * nm_sample_ct + sample_idx] = 1.0;
              }
            }
          }
          uintptr_t alt1_het_ct = genocounts[1] - pgv.patch_01_ct;
          if (pgv.patch_10_ct) {
            const uintptr_t* patch_set_nm = pgv.patch_10_set;
            if (missing_ct) {
              CopyBitarrSubset(pgv.patch_10_set, sample_nm, nm_sample_ct, tmp_nm);
              patch_set_nm = tmp_nm;
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t cur_bits = patch_set_nm[0];
            if (!omitted_allele_idx) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] = 2.0;
                  alt1_start[sample_idx] = 0.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ct + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] += 1.0;
                    alt1_start[sample_idx] = 0.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else if (omitted_allele_idx == 1) {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const AlleleCode ac0 = pgv.patch_10_vals[2 * uii];
                const AlleleCode ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] = 2.0;
                  machr2_dosage_ssqs[ac0] += 1;
                } else {
                  rarealt_start[(ac1 - 2) * nm_sample_ct + sample_idx] = 1.0;
                  machr2_dosage_sums[ac1] += 1;
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                  } else {
                    rarealt_start[(ac0 - 2) * nm_sample_ct + sample_idx] += 1.0;
                    machr2_dosage_sums[ac0] += 1;
                  }
                }
              }
            } else {
              for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                const uintptr_t sample_idx = BitIter1(patch_set_nm, &sample_idx_base, &cur_bits);
                const uint32_t ac0 = pgv.patch_10_vals[2 * uii];
                const uint32_t ac1 = pgv.patch_10_vals[2 * uii + 1];
                if (ac0 == ac1) {
                  machr2_dosage_ssqs[ac0] += 1;
                  alt1_start[sample_idx] = 0.0;
                  if (ac0 != omitted_allele_idx) {
                    const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                    rarealt_start[ac0_col * nm_sample_ct + sample_idx] = 2.0;
                  }
                } else {
                  machr2_dosage_sums[ac1] += 1;
                  if (ac1 != omitted_allele_idx) {
                    const uint32_t ac1_col = ac1 - 2 - (ac1 > omitted_allele_idx);
                    rarealt_start[ac1_col * nm_sample_ct + sample_idx] = 1.0;
                  }
                  if (ac0 == 1) {
                    ++alt1_het_ct;
                    alt1_start[sample_idx] = 1.0;
                  } else {
                    machr2_dosage_sums[ac0] += 1;
                    alt1_start[sample_idx] = 0.0;
                    if (ac0 != omitted_allele_idx) {
                      const uint32_t ac0_col = ac0 - 2 - (ac0 > omitted_allele_idx);
                      rarealt_start[ac0_col * nm_sample_ct + sample_idx] += 1.0;
                    }
                  }
                }
              }
            }
          }
          for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
            const uintptr_t one_ct = machr2_dosage_sums[allele_idx];
            const uintptr_t two_ct = machr2_dosage_ssqs[allele_idx];
            machr2_dosage_sums[allele_idx] = one_ct + 2 * two_ct;
            machr2_dosage_ssqs[allele_idx] = one_ct + 4LLU * two_ct;
            if ((one_ct == nm_sample_ct) || (two_ct == nm_sample_ct) || ((!one_ct) && (!two_ct))) {
              SetBit(allele_idx, const_alleles);
            }
          }
          const uintptr_t alt1_hom_ct = genocounts[2] - pgv.patch_10_ct;
          machr2_dosage_sums[1] = alt1_het_ct + 2 * alt1_hom_ct;
          machr2_dosage_ssqs[1] = alt1_het_ct + 4LLU * alt1_hom_ct;
          if ((alt1_het_ct == nm_sample_ct) || (alt1_hom_ct == nm_sample_ct) || ((!alt1_het_ct) && (!alt1_hom_ct))) {
            SetBit(1, const_alleles);
          }
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            machr2_dosage_sums[allele_idx] *= 0x4000LLU;
            machr2_dosage_ssqs[allele_idx] *= 0x10000000LLU;
          }
        }
        // usually need to save some of {sample_obs_ct, allele_obs_ct,
        // a1_dosage, mach_r2 even for skipped variants
        // compute them all for now, could conditionally skip later
        uint32_t allele_obs_ct = nm_sample_ct * 2;
        if (!is_x) {
          if (is_nonx_haploid) {
            allele_obs_ct = nm_sample_ct;
            // everything is on 0..1 scale, not 0..2
            if (!sparse_optimization) {
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                genotype_vals[sample_idx] *= 0.5;
              }
              const uint32_t high_ct = nm_sample_ct * allele_ct_m2;
              for (uint32_t uii = 0; uii != high_ct; ++uii) {
                multi_start[uii] *= 0.5;
              }
            }
          }
        } else {
          CopyBitarrSubset(sex_male_collapsed, sample_nm, nm_sample_ct, tmp_nm);
          const uintptr_t* male_nm = tmp_nm;
          const uint32_t nm_male_ct = PopcountWords(male_nm, nm_sample_ctl);
          if (is_xchr_model_1) {
            // special case: multiply male values by 0.5
            uintptr_t sample_idx_base = 0;
            uintptr_t male_nm_bits = male_nm[0];
            for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
              const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
              genotype_vals[sample_idx] *= 0.5;
              // could insert multiallelic loop here instead, but I'm guessing
              // that's worse due to locality of writes?
            }
            for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
              double* cur_start = &(multi_start[extra_allele_idx * nm_sample_ct]);
              sample_idx_base = 0;
              male_nm_bits = male_nm[0];
              for (uint32_t male_idx = 0; male_idx != nm_male_ct; ++male_idx) {
                const uintptr_t sample_idx = BitIter1(male_nm, &sample_idx_base, &male_nm_bits);
                cur_start[sample_idx] *= 0.5;
              }
            }
            allele_obs_ct -= nm_male_ct;
          }
        }
        const double mach_r2 = MultiallelicDiploidMachR2(machr2_dosage_sums, machr2_dosage_ssqs, nm_sample_ct, allele_ct);
        for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
          if (allele_idx == omitted_allele_idx) {
            continue;
          }
          block_aux_iter->sample_obs_ct = nm_sample_ct;
          block_aux_iter->allele_obs_ct = allele_obs_ct;
          double a1_dosage = u63tod(machr2_dosage_sums[allele_idx]) * kRecipDosageMid;
          if (is_xchr_model_1) {
            // ugh.
            double* geno_col = genotype_vals;
            if (allele_idx > (!omitted_allele_idx)) {
              geno_col = &(nm_predictors_pmaj_buf[(expected_predictor_ct - (allele_ct - allele_idx) + (allele_idx < omitted_allele_idx)) * nm_sample_ct]);
            }
            a1_dosage = 0.0;
            for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
              a1_dosage += geno_col[sample_idx];
            }
            if (!allele_ct_m2) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = 0.0;
              for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                const double cur_dosage = geno_col[sample_idx];
                main_dosage_ssq += cur_dosage * cur_dosage;
              }
            }
          } else {
            if (is_nonx_haploid) {
              a1_dosage *= 0.5;
            }
            if (!allele_ct_m2) {
              main_dosage_sum = a1_dosage;
              main_dosage_ssq = u63tod(machr2_dosage_ssqs[allele_idx]) * kRecipDosageMidSq;
              if (is_nonx_haploid) {
                main_dosage_ssq *= 0.25;
              }
            }
          }
          block_aux_iter->a1_dosage = a1_dosage;
          block_aux_iter->mach_r2 = mach_r2;
          ++block_aux_iter;
        }
        // now free to skip the actual regression if there are too few samples,
        // or there's a zero-variance genotype column
        GlmErr glm_err = 0;
        if (nm_sample_ct <= expected_predictor_ct) {
          glm_err = SetGlmErr0(kGlmErrcodeSampleCtLtePredictorCt);
        } else if (IsSet(const_alleles, omitted_allele_idx)) {
          glm_err = SetGlmErr0(kGlmErrcodeConstOmittedAllele);
        }
        if (glm_err) {
          if (missing_ct) {
            // covariates have not been copied yet, so we can't usually change
            // prev_nm from 0 to 1 when missing_ct == 0 (and there's little
            // reason to optimize the zero-covariate case).
            prev_nm = 0;
          }
          uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
          if (beta_se_multiallelic_fused || (!hide_covar)) {
            reported_ct += allele_ct_m2;
          }
          const uintptr_t result_row_ct = (extra_regression_ct + 1) * subbatch_size;
          for (uintptr_t ulii = 0; ulii != result_row_ct; ++ulii) {
            for (uint32_t uii = 0; uii != reported_ct; ++uii) {
              memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
              beta_se_iter[uii * 2 + 1] = -9.0;
            }
            beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
          }
        } else {
          uint32_t parameter_uidx = 2 + domdev_present;
          double* nm_predictors_pmaj_istart = nullptr;
          if (!sparse_optimization) {
            // only need to do this part once per variant in multiallelic case
            double* nm_predictors_pmaj_iter = &(nm_predictors_pmaj_buf[nm_sample_ct * (parameter_uidx - main_omitted)]);
            if (missing_ct || (!prev_nm)) {
              // fill phenotypes
              for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                const double* cur_pheno = &(cur_pheno_pmaj[pheno_idx * cur_sample_ct]);
                double* nm_pheno_col = &(nm_pheno_buf[pheno_idx * nm_sample_ct]);
                uintptr_t sample_midx_base = 0;
                uintptr_t sample_nm_bits = sample_nm[0];
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                  nm_pheno_col[sample_idx] = cur_pheno[sample_midx];
                }
              }

              // fill covariates
              for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx, ++parameter_uidx) {
                // strictly speaking, we don't need cur_covars_cmaj to be
                // vector-aligned
                if (cur_parameter_subset && (!IsSet(cur_parameter_subset, parameter_uidx))) {
                  continue;
                }
                const double* cur_covar_col;
                if (covar_idx < local_covar_ct) {
                  cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                } else {
                  cur_covar_col = &(cur_covars_cmaj[(covar_idx - local_covar_ct) * cur_sample_ct]);
                }
                uintptr_t sample_midx_base = 0;
                uintptr_t sample_nm_bits = sample_nm[0];
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                  *nm_predictors_pmaj_iter++ = cur_covar_col[sample_midx];
                }
              }
              nm_predictors_pmaj_istart = nm_predictors_pmaj_iter;
              prev_nm = !missing_ct;
            } else {
              // bugfix (15 Aug 2018): this was not handling --parameters
              // correctly when a covariate was only needed as part of an
              // interaction
              parameter_uidx += cur_covar_ct;
              nm_predictors_pmaj_istart = &(nm_predictors_pmaj_iter[literal_covar_ct * nm_sample_ct]);
            }
          }
          const uint32_t const_allele_ct = PopcountWords(const_alleles, allele_ctl);
          if (const_allele_ct) {
            // Must delete constant-allele columns from nm_predictors_pmaj, and
            // shift later columns back.
            double* read_iter = genotype_vals;
            double* write_iter = genotype_vals;
            for (uint32_t read_allele_idx = 0; read_allele_idx != allele_ct; ++read_allele_idx) {
              if (read_allele_idx == omitted_allele_idx) {
                continue;
              }
              if (!IsSet(const_alleles, read_allele_idx)) {
                if (write_iter != read_iter) {
                  memcpy(write_iter, read_iter, nm_sample_ct * sizeof(double));
                }
                if (write_iter == genotype_vals) {
                  write_iter = multi_start;
                } else {
                  write_iter = &(write_iter[nm_sample_ct]);
                }
              }
              if (read_iter == genotype_vals) {
                read_iter = multi_start;
              } else {
                read_iter = &(read_iter[nm_sample_ct]);
              }
            }
          }
          const uint32_t cur_predictor_ct = expected_predictor_ct - const_allele_ct;
          uint32_t nonconst_extra_regression_idx = UINT32_MAX;  // deliberate overflow
          for (uint32_t extra_regression_idx = 0; extra_regression_idx <= extra_regression_ct; ++extra_regression_idx) {
            if (extra_regression_ct) {
              if (IsSet(const_alleles, extra_regression_idx + (extra_regression_idx >= omitted_allele_idx))) {
                glm_err = SetGlmErr0(kGlmErrcodeConstAllele);
                goto GlmLinearSubbatchThread_skip_regression;
              }
              ++nonconst_extra_regression_idx;
              if (nonconst_extra_regression_idx) {
                double* swap_target = &(multi_start[(nonconst_extra_regression_idx - 1) * nm_sample_ct]);
                for (uint32_t uii = 0; uii != nm_sample_ct; ++uii) {
                  double dxx = genotype_vals[uii];
                  genotype_vals[uii] = swap_target[uii];
                  swap_target[uii] = dxx;
                }
              }
            }
            if (!sparse_optimization) {
              double* main_vals = &(nm_predictors_pmaj_buf[nm_sample_ct]);
              double* domdev_vals = nullptr;
              if (main_omitted) {
                // if main_mutated, this will be filled below
                // if not, this aliases genotype_vals
                main_vals = &(nm_predictors_pmaj_buf[(cur_predictor_ct + main_mutated) * nm_sample_ct]);
              } else if (joint_genotypic || joint_hethom) {
                domdev_vals = &(main_vals[nm_sample_ct]);
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 2.0 - cur_genotype_val;
                  }
                  domdev_vals[sample_idx] = cur_genotype_val;
                }
              }
              if (model_dominant) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..1..1
                  if (cur_genotype_val > 1.0) {
                    cur_genotype_val = 1.0;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              } else if (model_recessive || joint_hethom) {
                for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                  double cur_genotype_val = genotype_vals[sample_idx];
                  // 0..0..1
                  if (cur_genotype_val < 1.0) {
                    cur_genotype_val = 0.0;
                  } else {
                    cur_genotype_val -= 1.0;
                  }
                  main_vals[sample_idx] = cur_genotype_val;
                }
              }

              // fill interaction terms
              if (add_interactions) {
                double* nm_predictors_pmaj_iter = nm_predictors_pmaj_istart;
                for (uint32_t covar_idx = 0; covar_idx != cur_covar_ct; ++covar_idx) {
                  const double* cur_covar_col;
                  if (covar_idx < local_covar_ct) {
                    cur_covar_col = &(local_covars_iter[covar_idx * max_sample_ct]);
                  } else {
                    cur_covar_col = &(cur_covars_cmaj[covar_idx * cur_sample_ct]);
                  }
                  if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                    uintptr_t sample_midx_base = 0;
                    uintptr_t sample_nm_bits = sample_nm[0];
                    for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                      const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                      *nm_predictors_pmaj_iter++ = main_vals[sample_idx] * cur_covar_col[sample_midx];
                    }
                  }
                  ++parameter_uidx;
                  if (domdev_present) {
                    if ((!cur_parameter_subset) || IsSet(cur_parameter_subset, parameter_uidx)) {
                      uintptr_t sample_midx_base = 0;
                      uintptr_t sample_nm_bits = sample_nm[0];
                      for (uint32_t sample_idx = 0; sample_idx != nm_sample_ct; ++sample_idx) {
                        const uintptr_t sample_midx = BitIter1(sample_nm, &sample_midx_base, &sample_nm_bits);
                        *nm_predictors_pmaj_iter++ = domdev_vals[sample_idx] * cur_covar_col[sample_midx];
                      }
                    }
                    ++parameter_uidx;
                  }
                }
              }
            }

            // bugfix (12 Sep 2017): forgot to implement per-variant VIF and
            // max-corr checks
            if (xtx_image && prev_nm && (!allele_ct_m2)) {
              // only need to fill in additive and possibly domdev dot
              // products
              memcpy(xtx_inv, xtx_image, cur_predictor_ct * cur_predictor_ct * sizeof(double));
              memcpy(xt_y, xt_y_image, cur_predictor_ct * subbatch_size * sizeof(double));
              if (sparse_optimization) {
                // currently does not handle chrX
                ZeroDArr(subbatch_size, geno_pheno_prods);
                if (domdev_pheno_prods) {
                  ZeroDArr(subbatch_size, domdev_pheno_prods);
                }
                double domdev_geno_prod = 0.0;
                double* geno_dotprod_row = &(xtx_inv[cur_predictor_ct]);
                double* domdev_dotprod_row = &(xtx_inv[2 * cur_predictor_ct]);
                for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
                  uintptr_t geno_word = pgv.genovec[widx];
                  if (geno_word) {
                    const uint32_t sample_idx_base = widx * kBitsPerWordD2;
                    do {
                      const uint32_t lowest_set_bit = ctzw(geno_word);
                      // since there are no missing values, we have a het if
                      // (lowest_set_bit & 1) is zero, and a hom-alt when
                      // it's one.
                      const uint32_t sample_idx = sample_idx_base + (lowest_set_bit / 2);
                      const double geno_d = geno_d_lookup[lowest_set_bit & 1];
                      for (uintptr_t pred_idx = domdev_present + 2; pred_idx != cur_predictor_ct; ++pred_idx) {
                        geno_dotprod_row[pred_idx] += geno_d * nm_predictors_pmaj_buf[pred_idx * nm_sample_ct + sample_idx];
                      }
                      // can have a separate categorical loop here

                      if (domdev_present && (!(lowest_set_bit & 1))) {
                        // domdev = 1
                        for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                          const double cur_pheno_val = nm_pheno_buf[pheno_idx * nm_sample_ct + sample_idx];
                          geno_pheno_prods[pheno_idx] += geno_d * cur_pheno_val;
                          domdev_pheno_prods[pheno_idx] += cur_pheno_val;
                        }
                        domdev_geno_prod += geno_d;
                        for (uintptr_t pred_idx = 3; pred_idx != cur_predictor_ct; ++pred_idx) {
                          domdev_dotprod_row[pred_idx] += nm_predictors_pmaj_buf[pred_idx * nm_sample_ct + sample_idx];
                        }
                        // categorical optimization possible here
                      } else {
                        for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                          geno_pheno_prods[pheno_idx] += geno_d * nm_pheno_buf[pheno_idx * nm_sample_ct + sample_idx];
                        }
                      }
                      geno_word &= geno_word - 1;
                    } while (geno_word);
                  }
                }
                const double het_ctd = u31tod(genocounts[1]);
                const double homalt_ctd = u31tod(genocounts[2]);
                xtx_inv[cur_predictor_ct] = het_ctd * geno_d_lookup[0] + homalt_ctd * geno_d_lookup[1];
                xtx_inv[cur_predictor_ct + 1] = het_ctd * geno_d_lookup[0] * geno_d_lookup[0] + homalt_ctd * geno_d_lookup[1] * geno_d_lookup[1];
                if (domdev_present) {
                  xtx_inv[cur_predictor_ct + 2] = domdev_geno_prod;
                  xtx_inv[2 * cur_predictor_ct] = het_ctd;
                  xtx_inv[2 * cur_predictor_ct + 2] = het_ctd;
                }
                for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                  double* cur_xt_y = &(xt_y[cur_predictor_ct * pheno_idx]);
                  cur_xt_y[1] = geno_pheno_prods[pheno_idx];
                  if (domdev_present) {
                    cur_xt_y[2] = domdev_pheno_prods[pheno_idx];
                  }
                }
              } else {
                // !sparse_optimization
                uintptr_t start_pred_idx = 0;
                if (!(model_dominant || model_recessive || joint_hethom)) {
                  start_pred_idx = domdev_present + 2;
                  xtx_inv[cur_predictor_ct] = main_dosage_sum;
                  xtx_inv[cur_predictor_ct + 1] = main_dosage_ssq;
                }
                if (cur_predictor_ct > start_pred_idx) {
                  // categorical optimization possible here
                  ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[nm_sample_ct]), &(nm_predictors_pmaj_buf[start_pred_idx * nm_sample_ct]), nm_sample_ct, nm_sample_ct, cur_predictor_ct - start_pred_idx, &(xtx_inv[cur_predictor_ct + start_pred_idx]));
                }
                RowMajorMatrixMultiplyStrided(nm_pheno_buf, &(nm_predictors_pmaj_buf[nm_sample_ct]), subbatch_size, nm_sample_ct, 1, 1, nm_sample_ct, cur_predictor_ct, &(xt_y[1]));
                if (domdev_present) {
                  RowMajorMatrixMultiplyStrided(nm_pheno_buf, &(nm_predictors_pmaj_buf[2 * nm_sample_ct]), subbatch_size, nm_sample_ct, 1, 1, nm_sample_ct, cur_predictor_ct, &(xt_y[2]));
                  // categorical optimization possible here
                  ColMajorVectorMatrixMultiplyStrided(&(nm_predictors_pmaj_buf[2 * nm_sample_ct]), nm_predictors_pmaj_buf, nm_sample_ct, nm_sample_ct, cur_predictor_ct, &(xtx_inv[2 * cur_predictor_ct]));
                  xtx_inv[cur_predictor_ct + 2] = xtx_inv[2 * cur_predictor_ct + 1];
                }
              }
              glm_err = CheckMaxCorrAndVifNm(xtx_inv, corr_inv, cur_predictor_ct, domdev_present_p1, cur_sample_ct_recip, cur_sample_ct_m1_recip, max_corr, vif_thresh, semicomputed_biallelic_corr_matrix, semicomputed_biallelic_inv_corr_sqrts, dbl_2d_buf, &(dbl_2d_buf[2 * cur_predictor_ct]), &(dbl_2d_buf[3 * cur_predictor_ct]));
              if (glm_err) {
                goto GlmLinearSubbatchThread_skip_regression;
              }
              const double geno_ssq = xtx_inv[1 + cur_predictor_ct];
              if (!domdev_present) {
                xtx_inv[1 + cur_predictor_ct] = xtx_inv[cur_predictor_ct];
                if (InvertRank1Symm(covarx_dotprod_inv, &(xtx_inv[1 + cur_predictor_ct]), cur_predictor_ct - 1, 1, geno_ssq, dbl_2d_buf, inverse_corr_buf)) {
                  glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                  goto GlmLinearSubbatchThread_skip_regression;
                }
              } else {
                const double domdev_geno_prod = xtx_inv[2 + cur_predictor_ct];
                const double domdev_ssq = xtx_inv[2 + 2 * cur_predictor_ct];
                xtx_inv[2 + cur_predictor_ct] = xtx_inv[cur_predictor_ct];
                xtx_inv[2 + 2 * cur_predictor_ct] = xtx_inv[2 * cur_predictor_ct];
                if (InvertRank2Symm(covarx_dotprod_inv, &(xtx_inv[2 + cur_predictor_ct]), cur_predictor_ct - 2, cur_predictor_ct, 1, geno_ssq, domdev_geno_prod, domdev_ssq, dbl_2d_buf, inverse_corr_buf, &(inverse_corr_buf[2 * (cur_predictor_ct - 2)]))) {
                  glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                  goto GlmLinearSubbatchThread_skip_regression;
                }
              }
              // need to make sure xtx_inv remains reflected in NOLAPACK case
              memcpy(xtx_inv, dbl_2d_buf, cur_predictor_ct * cur_predictor_ct * sizeof(double));
              ReflectMatrix(cur_predictor_ct, xtx_inv);
              RowMajorMatrixMultiply(xt_y, xtx_inv, subbatch_size, cur_predictor_ct, cur_predictor_ct, fitted_coefs);
            } else {
              // generic case
              // major categorical optimization possible here, some
              // multiallelic optimizations possible
              MultiplySelfTranspose(nm_predictors_pmaj_buf, cur_predictor_ct, nm_sample_ct, xtx_inv);

              for (uint32_t pred_idx = 1; pred_idx != cur_predictor_ct; ++pred_idx) {
                dbl_2d_buf[pred_idx] = xtx_inv[pred_idx * cur_predictor_ct];
              }
              glm_err = CheckMaxCorrAndVif(xtx_inv, 1, cur_predictor_ct, nm_sample_ct, max_corr, vif_thresh, dbl_2d_buf, nullptr, inverse_corr_buf, inv_1d_buf);
              if (glm_err) {
                goto GlmLinearSubbatchThread_skip_regression;
              }
              if (LinearRegressionInv(nm_pheno_buf, nm_predictors_pmaj_buf, cur_predictor_ct, nm_sample_ct, subbatch_size, xtx_inv, fitted_coefs, xt_y, inv_1d_buf, dbl_2d_buf)) {
                glm_err = SetGlmErr0(kGlmErrcodeRankDeficient);
                goto GlmLinearSubbatchThread_skip_regression;
              }
            }
            // RSS = y^T y - y^T X (X^T X)^{-1} X^T y
            //     = cur_pheno_ssq - xt_y * fitted_coefs
            // s^2 = RSS / df
            // possible todo: improve numerical stability of this computation
            // in non-mean-centered phenotype case

            for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
              {
                double tmp_pheno_ssq = pheno_ssq_bases[pheno_idx];
                if (missing_ct) {
                  const double* cur_pheno = &(cur_pheno_pmaj[pheno_idx * cur_sample_ct]);
                  uintptr_t sample_midx_base = 0;
                  uintptr_t sample_nm_inv_bits = ~sample_nm[0];
                  for (uint32_t missing_idx = 0; missing_idx != missing_ct; ++missing_idx) {
                    const uintptr_t sample_midx = BitIter0(sample_nm, &sample_midx_base, &sample_nm_inv_bits);
                    tmp_pheno_ssq -= cur_pheno[sample_midx] * cur_pheno[sample_midx];
                  }
                }
                const double* tmp_fitted_coefs = &(fitted_coefs[pheno_idx * cur_predictor_ct]);
                const double sigma = (tmp_pheno_ssq - DotprodxD(&(xt_y[pheno_idx * cur_predictor_ct]), tmp_fitted_coefs, cur_predictor_ct)) / u31tod(nm_sample_ct - cur_predictor_ct);
                for (uint32_t uii = 0; uii != cur_predictor_ct; ++uii) {
                  const double* s_iter = &(xtx_inv[uii * cur_predictor_ct]);
                  double* s_iter2 = &(xtx_inv2[uii * cur_predictor_ct]);
#ifdef NOLAPACK
                  for (uint32_t ujj = 0; ujj != cur_predictor_ct; ++ujj) {
                    s_iter2[ujj] = s_iter[ujj] * sigma;
                  }
#else
                  for (uint32_t ujj = 0; ujj <= uii; ++ujj) {
                    s_iter2[ujj] = s_iter[ujj] * sigma;
                  }
#endif
                }
                // validParameters() check
                for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
                  const double xtx_inv2_diag_element = xtx_inv2[pred_uidx * (cur_predictor_ct + 1)];
                  if (xtx_inv2_diag_element < 1e-20) {
                    glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                    goto GlmLinearSubbatchThread_skip_regression_for_one_pheno;
                  }
                  // use dbl_2d_buf[] to store diagonal square roots
                  dbl_2d_buf[pred_uidx] = sqrt(xtx_inv2_diag_element);
                }
                dbl_2d_buf[0] = sqrt(xtx_inv2[0]);
                for (uint32_t pred_uidx = 1; pred_uidx != cur_predictor_ct; ++pred_uidx) {
                  const double cur_xtx_inv2_diag_sqrt = 0.99999 * dbl_2d_buf[pred_uidx];
                  const double* xtx_inv2_row = &(xtx_inv2[pred_uidx * cur_predictor_ct]);
                  for (uint32_t pred_uidx2 = 0; pred_uidx2 != pred_uidx; ++pred_uidx2) {
                    if (xtx_inv2_row[pred_uidx2] > cur_xtx_inv2_diag_sqrt * dbl_2d_buf[pred_uidx2]) {
                      glm_err = SetGlmErr0(kGlmErrcodeInvalidResult);
                      goto GlmLinearSubbatchThread_skip_regression_for_one_pheno;
                    }
                  }
                }
                double* beta_se_iter2 = beta_se_iter;
                for (uint32_t pred_uidx = reported_pred_uidx_start; pred_uidx != reported_pred_uidx_biallelic_end; ++pred_uidx) {
                  // In the multiallelic-fused case, if the first allele is
                  // constant, this writes the beta/se values for the first
                  // nonconstant, non-omitted allele where the results for the
                  // first allele belong.  We correct that below.
                  *beta_se_iter2++ = tmp_fitted_coefs[pred_uidx];
                  *beta_se_iter2++ = dbl_2d_buf[pred_uidx];
                }
                // move this up since dbl_2d_buf may be clobbered
                if (cur_constraint_ct) {
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeRankDeficient);
                  memcpy(beta_se_iter2, &glm_err2, 8);
                  ++beta_se_iter2;
                  *beta_se_iter2++ = -9.0;
                }
                if (!const_allele_ct) {
                  if (beta_se_multiallelic_fused || (!hide_covar)) {
                    for (uint32_t extra_allele_idx = 0; extra_allele_idx != allele_ct_m2; ++extra_allele_idx) {
                      beta_se_iter2[2 * extra_allele_idx] = tmp_fitted_coefs[cur_biallelic_predictor_ct + extra_allele_idx];
                      beta_se_iter2[2 * extra_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + extra_allele_idx];
                    }
                  }
                } else if (!beta_se_multiallelic_fused) {
                  if (!hide_covar) {
                    // Need to insert some {CONST_ALLELE, -9} entries.
                    const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                    const uint32_t cur_raw_allele_idx = extra_regression_idx + (extra_regression_idx >= omitted_allele_idx);
                    uint32_t extra_read_allele_idx = 0;
                    uint32_t extra_write_allele_idx = 0;
                    for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                      if ((allele_idx == omitted_allele_idx) || (allele_idx == cur_raw_allele_idx)) {
                        continue;
                      }
                      if (IsSet(const_alleles, allele_idx)) {
                        memcpy(&(beta_se_iter2[2 * extra_write_allele_idx]), &glm_err2, 8);
                        beta_se_iter2[2 * extra_write_allele_idx + 1] = -9.0;
                      } else {
                        beta_se_iter2[2 * extra_write_allele_idx] = tmp_fitted_coefs[cur_biallelic_predictor_ct + extra_read_allele_idx];
                        beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + extra_read_allele_idx];
                        ++extra_read_allele_idx;
                      }
                      ++extra_write_allele_idx;
                    }
                  }
                } else {
                  const GlmErr glm_err2 = SetGlmErr0(kGlmErrcodeConstAllele);
                  // Special-case first nonconst allele since it's positioned
                  // discontinuously, and its BETA/SE may already be correctly
                  // filled.
                  uint32_t allele_idx = omitted_allele_idx? 0 : 1;
                  uint32_t extra_write_allele_idx = 0;
                  if (IsSet(const_alleles, allele_idx)) {
                    memcpy(&(beta_se_iter[2 * include_intercept]), &glm_err2, 8);
                    beta_se_iter[2 * include_intercept + 1] = -9.0;
                    allele_idx = AdvTo0Bit(const_alleles, 1);
                    if (allele_idx == omitted_allele_idx) {
                      allele_idx = AdvTo0Bit(const_alleles, omitted_allele_idx + 1);
                    }
                    extra_write_allele_idx = allele_idx - 1 - (allele_idx > omitted_allele_idx);
                    for (uint32_t uii = 0; uii != extra_write_allele_idx; ++uii) {
                      memcpy(&(beta_se_iter2[2 * uii]), &glm_err2, 8);
                      beta_se_iter2[2 * uii + 1] = -9.0;
                    }
                    beta_se_iter2[2 * extra_write_allele_idx] = tmp_fitted_coefs[1];
                    beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[1];
                    ++extra_write_allele_idx;
                  }
                  ++allele_idx;
                  uint32_t nonconst_allele_idx_m1 = 0;
                  for (; allele_idx != allele_ct; ++allele_idx) {
                    if (allele_idx == omitted_allele_idx) {
                      continue;
                    }
                    if (!IsSet(const_alleles, allele_idx)) {
                      beta_se_iter2[2 * extra_write_allele_idx] = tmp_fitted_coefs[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                      beta_se_iter2[2 * extra_write_allele_idx + 1] = dbl_2d_buf[cur_biallelic_predictor_ct + nonconst_allele_idx_m1];
                      ++nonconst_allele_idx_m1;
                    } else {
                      memcpy(&(beta_se_iter2[2 * extra_write_allele_idx]), &glm_err2, 8);
                      beta_se_iter2[2 * extra_write_allele_idx + 1] = -9.0;
                    }
                    ++extra_write_allele_idx;
                  }
                }
                if (cur_constraint_ct) {
                  // probable todo: cur_constraints_con_major manipulation can
                  // go outside the pheno_idx loop
                  uint32_t joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                  for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                    joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                    cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 1.0;
                  }
#ifndef NOLAPACK
                  // xtx_inv2 upper triangle was not filled
                  ReflectMatrix(cur_predictor_ct, xtx_inv2);
#endif
                  double chisq;
                  if (!LinearHypothesisChisq(tmp_fitted_coefs, cur_constraints_con_major, xtx_inv2, cur_constraint_ct, cur_predictor_ct, &chisq, tmphxs_buf, h_transpose_buf, inner_buf, inv_1d_buf, dbl_2d_buf)) {
                    beta_se_iter2[-1] = chisq;
                  }
                  // next test may have different alt allele count
                  joint_test_idx = AdvTo1Bit(cur_joint_test_params, 0);
                  for (uint32_t uii = 1; uii != cur_constraint_ct; ++uii) {
                    joint_test_idx = AdvTo1Bit(cur_joint_test_params, joint_test_idx + 1);
                    cur_constraints_con_major[uii * cur_predictor_ct + joint_test_idx] = 0.0;
                  }
                }
              }
              while (0) {
              GlmLinearSubbatchThread_skip_regression_for_one_pheno:
                {
                  uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
                  if (beta_se_multiallelic_fused || (!hide_covar)) {
                    reported_ct += allele_ct_m2;
                  }
                  for (uint32_t uii = 0; uii != reported_ct; ++uii) {
                    memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
                    beta_se_iter[uii * 2 + 1] = -9.0;
                  }
                }
              }
              beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
            }
            while (0) {
            GlmLinearSubbatchThread_skip_regression:
              {
                uint32_t reported_ct = reported_pred_uidx_biallelic_end + (cur_constraint_ct != 0) - reported_pred_uidx_start;
                if (beta_se_multiallelic_fused || (!hide_covar)) {
                  reported_ct += allele_ct_m2;
                }
                for (uint32_t pheno_idx = 0; pheno_idx != subbatch_size; ++pheno_idx) {
                  for (uint32_t uii = 0; uii != reported_ct; ++uii) {
                    memcpy(&(beta_se_iter[uii * 2]), &glm_err, 8);
                    beta_se_iter[uii * 2 + 1] = -9.0;
                  }
                  beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
                }
              }
            }
          }
        }
        if (local_covars_iter) {
          local_covars_iter = &(local_covars_iter[local_covar_ct * max_sample_ct]);
        }
      }
    }
    while (0) {
    GlmLinearSubbatchThread_err:
      UpdateU64IfSmaller(new_err_info, &common->err_info);
      break;
    }
    parity = 1 - parity;
    variant_idx_offset += cur_block_variant_ct;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

CONSTI32(kMaxLinearSubbatchSize, 240);
static_assert(kMaxLinearSubbatchSize + 12 <= kMaxOpenFiles, "kMaxLinearSubbatchSize can't be too close to or larger than kMaxOpenFiles.");

PglErr GlmLinearBatch(const uintptr_t* pheno_batch, const PhenoCol* pheno_cols, const char* pheno_names, const char* const* test_names, const char* const* test_names_x, const char* const* test_names_y, const uint32_t* variant_bps, const char* const* variant_ids, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const uint32_t* local_sample_uidx_order, const uintptr_t* local_variant_include, uint32_t raw_variant_ct, uint32_t completed_pheno_ct, uint32_t batch_size, uintptr_t max_pheno_name_blen, uint32_t max_chr_blen, double ci_size, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, uintptr_t overflow_buf_size, uint32_t local_sample_ct, PgenFileInfo* pgfip, GlmLinearCtx* ctx, TextStream* local_covar_txsp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char** cswritep_arr = nullptr;
  CompressStreamState* css_arr = nullptr;
  uint32_t subbatch_size = 0;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  {
    GlmCtx* common = ctx->common;
    const uintptr_t* variant_include = common->variant_include;
    const ChrInfo* cip = common->cip;
    const uintptr_t* allele_idx_offsets = common->allele_idx_offsets;
    const AlleleCode* omitted_alleles = common->omitted_alleles;

    const uint32_t sample_ct = common->sample_ct;
    const uint32_t sample_ct_x = common->sample_ct_x;
    const uint32_t sample_ct_y = common->sample_ct_y;
    const uint32_t covar_ct = common->covar_ct;
    const uintptr_t local_covar_ct = common->local_covar_ct;
    const uint32_t covar_ct_x = common->covar_ct_x;
    const uint32_t covar_ct_y = common->covar_ct_y;

    uint32_t max_sample_ct = MAXV(sample_ct, sample_ct_x);
    if (max_sample_ct < sample_ct_y) {
      max_sample_ct = sample_ct_y;
    }
    const uintptr_t* sample_include = common->sample_include;
    uint32_t* local_sample_idx_order = nullptr;
    uint32_t local_line_idx = 0;
    uint32_t local_xy = 0;  // 1 = chrX, 2 = chrY

    const char* local_line_iter = nullptr;
    uint32_t local_prev_chr_code = UINT32_MAX;
    uint32_t local_chr_code = UINT32_MAX;
    uint32_t local_bp = UINT32_MAX;
    uint32_t local_skip_chr = 1;
    if (local_covar_ct) {
      reterr = TextRewind(local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLinearBatch_ret_TSTREAM_FAIL;
      }
      local_line_idx = glm_info_ptr->local_header_line_ct;
      reterr = TextSkip(local_line_idx, local_covar_txsp);
      if (unlikely(reterr)) {
        goto GlmLinearBatch_ret_TSTREAM_FAIL;
      }
      if (unlikely(bigstack_alloc_u32(local_sample_ct, &local_sample_idx_order))) {
        goto GlmLinearBatch_ret_NOMEM;
      }
      for (uint32_t uii = 0; uii != local_sample_ct; ++uii) {
        const uint32_t cur_uidx = local_sample_uidx_order[uii];
        uint32_t cur_idx = UINT32_MAX;
        if ((cur_uidx != UINT32_MAX) && IsSet(sample_include, cur_uidx)) {
          cur_idx = RawToSubsettedPos(sample_include, common->sample_include_cumulative_popcounts, cur_uidx);
        }
        local_sample_idx_order[uii] = cur_idx;
      }
    }

    const uint32_t variant_ct = common->variant_ct;

    const GlmFlags glm_flags = glm_info_ptr->flags;
    const uint32_t report_neglog10p = (glm_flags / kfGlmLog10) & 1;
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t domdev_present_p1 = domdev_present + 1;

    const uint32_t constraint_ct = common->constraint_ct;
    const uint32_t constraint_ct_x = common->constraint_ct_x;
    const uint32_t constraint_ct_y = common->constraint_ct_y;

    const uint32_t max_extra_allele_ct = common->max_extra_allele_ct;
    uint32_t biallelic_predictor_ct = 2 + domdev_present + covar_ct * (1 + add_interactions * domdev_present_p1);
    uint32_t biallelic_predictor_ct_x = 2 + covar_ct_x * (1 + add_interactions);
    uint32_t biallelic_predictor_ct_y = 2 + covar_ct_y * (1 + add_interactions);
    const uintptr_t* parameter_subset = common->parameter_subset;
    const uintptr_t* parameter_subset_x = common->parameter_subset_x;
    const uintptr_t* parameter_subset_y = common->parameter_subset_y;
    if (parameter_subset) {
      biallelic_predictor_ct = PopcountWords(parameter_subset, BitCtToWordCt(biallelic_predictor_ct));
      if (sample_ct_x) {
        biallelic_predictor_ct_x = PopcountWords(parameter_subset_x, BitCtToWordCt(biallelic_predictor_ct_x));
      } else {
        biallelic_predictor_ct_x = 0;
      }
      if (sample_ct_y) {
        biallelic_predictor_ct_y = PopcountWords(parameter_subset_y, BitCtToWordCt(biallelic_predictor_ct_y));
      } else {
        biallelic_predictor_ct_y = 0;
      }
    }
    uint32_t biallelic_reported_test_ct = GetBiallelicReportedTestCt(parameter_subset, glm_flags, covar_ct, common->tests_flag);
    uintptr_t max_reported_test_ct = biallelic_reported_test_ct;
    uint32_t biallelic_reported_test_ct_x = 0;
    if (sample_ct_x) {
      biallelic_reported_test_ct_x = GetBiallelicReportedTestCt(parameter_subset_x, glm_flags, covar_ct_x, common->tests_flag);
      if (biallelic_reported_test_ct_x > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_x;
      }
    }
    uint32_t biallelic_reported_test_ct_y = 0;
    if (sample_ct_y) {
      biallelic_reported_test_ct_y = GetBiallelicReportedTestCt(parameter_subset_y, glm_flags, covar_ct_y, common->tests_flag);
      if (biallelic_reported_test_ct_y > max_reported_test_ct) {
        max_reported_test_ct = biallelic_reported_test_ct_y;
      }
    }
    const uint32_t hide_covar = (glm_flags / kfGlmHideCovar) & 1;
    const uint32_t include_intercept = (glm_flags / kfGlmIntercept) & 1;
    const GlmColFlags glm_cols = glm_info_ptr->cols;
    const uint32_t test_col = glm_cols & kfGlmColTest;
    if (unlikely((!test_col) && (max_reported_test_ct > 1))) {
      logerrputs("Error: --glm's 'test' column cannot be omitted when results for multiple\npredictors are reported.  (Did you forget 'hide-covar'?)\n");
      goto GlmLinearBatch_ret_INCONSISTENT_INPUT;
    }
    const uint32_t main_mutated = ((glm_flags & (kfGlmDominant | kfGlmRecessive | kfGlmHethom)) != kfGlm0);
    const uint32_t beta_se_multiallelic_fused = (!domdev_present) && (!main_mutated) && (!common->tests_flag) && (!add_interactions);
    if (beta_se_multiallelic_fused || (!hide_covar)) {
      max_reported_test_ct += max_extra_allele_ct;
    }
    common->max_reported_test_ct = max_reported_test_ct;

    uint32_t x_code = UINT32_MAXM1;
    uint32_t x_start = 0;
    uint32_t x_end = 0;
    if (sample_ct_x) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetX, &x_code, &x_start, &x_end);
    }
    uint32_t y_code = UINT32_MAXM1;
    uint32_t y_start = 0;
    uint32_t y_end = 0;
    if (sample_ct_y) {
      GetXymtCodeStartAndEndUnsafe(cip, kChrOffsetY, &y_code, &y_start, &y_end);
    }
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t chr_col = glm_cols & kfGlmColChrom;

    char* chr_buf = nullptr;
    if (chr_col) {
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto GlmLinearBatch_ret_NOMEM;
      }
    }

    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > variant_ct) {
      calc_thread_ct = variant_ct;
    }

    subbatch_size = MINV(kMaxLinearSubbatchSize, batch_size);
    if (sample_ct * S_CAST(uint64_t, subbatch_size) > 0x7fffffff) {
      subbatch_size = 0x7fffffff / sample_ct;
    }
    if (unlikely(
            bigstack_calloc_cp(subbatch_size, &cswritep_arr) ||
            BIGSTACK_ALLOC_X(CompressStreamState, subbatch_size, &css_arr))) {
      goto GlmLinearBatch_ret_NOMEM;
    }
    for (uint32_t fidx = 0; fidx != subbatch_size; ++fidx) {
      PreinitCstream(&(css_arr[fidx]));
    }
    // Subbatch-size-dependent memory allocations:
    //   nm_pheno_buf
    //   fitted_coefs
    //   xt_y
    //   block_beta_se_bufs
    //   css_arr

    const uint32_t main_omitted = (parameter_subset && (!IsSet(parameter_subset, 1)));
    const uint32_t xmain_ct = main_mutated + main_omitted;
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    common->thread_mhc = nullptr;
    common->dosage_presents = nullptr;
    common->dosage_mains = nullptr;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    // This cannot be less than what InitCstreamAlloc() actually allocates.
    uintptr_t cstream_alloc_size = RoundUpPow2(overflow_buf_size, kCacheline);
    if (output_zst) {
      cstream_alloc_size += RoundUpPow2(CstreamWkspaceReq(overflow_buf_size), kCacheline);
    }
    unsigned char* bigstack_mark2 = g_bigstack_base;
    uintptr_t workspace_alloc;
    uint32_t read_block_size;
    uintptr_t max_alt_allele_block_size;
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    for (; ; subbatch_size = (subbatch_size + 1) / 2) {
      // may permit size 1 here and get rid of GlmLinear() later
      if (subbatch_size == 1) {
        goto GlmLinearBatch_ret_NOMEM;
      }
      BigstackReset(bigstack_mark2);
      if (bigstack_alloc_d(sample_ct * subbatch_size, &ctx->pheno_d)) {
        continue;
      }
      if (common->nm_precomp) {
        if (bigstack_alloc_d((2 + domdev_present + covar_ct) * subbatch_size, &(common->nm_precomp->xt_y_image))) {
          continue;
        }
      }
      if (sample_ct_x) {
        if (bigstack_alloc_d(sample_ct_x * subbatch_size, &ctx->pheno_x_d)) {
          continue;
        }
        if (common->nm_precomp_x) {
          // domdev_present can't be set here.
          if (bigstack_alloc_d((2 + covar_ct_x) * subbatch_size, &(common->nm_precomp_x->xt_y_image))) {
            continue;
          }
        }
      }
      if (sample_ct_y) {
        if (bigstack_alloc_d(sample_ct_y * subbatch_size, &ctx->pheno_y_d)) {
          continue;
        }
        if (common->nm_precomp_y) {
          if (bigstack_alloc_d((2 + covar_ct_y) * subbatch_size, &(common->nm_precomp_y->xt_y_image))) {
            continue;
          }
        }
      }
      workspace_alloc = GetLinearSubbatchWorkspaceSize(sample_ct, subbatch_size, biallelic_predictor_ct, max_extra_allele_ct, constraint_ct, xmain_ct);
      if (sample_ct_x) {
        const uintptr_t workspace_alloc_x = GetLinearSubbatchWorkspaceSize(sample_ct_x, subbatch_size, biallelic_predictor_ct_x, max_extra_allele_ct, constraint_ct_x, xmain_ct);
        if (workspace_alloc_x > workspace_alloc) {
          workspace_alloc = workspace_alloc_x;
        }
      }
      if (sample_ct_y) {
        const uintptr_t workspace_alloc_y = GetLinearSubbatchWorkspaceSize(sample_ct_y, subbatch_size, biallelic_predictor_ct_y, max_extra_allele_ct, constraint_ct_y, xmain_ct);
        if (workspace_alloc_y > workspace_alloc) {
          workspace_alloc = workspace_alloc_y;
        }
      }
      uintptr_t thread_xalloc_cacheline_ct = (workspace_alloc / kCacheline) + 1;
      uintptr_t per_variant_xalloc_byte_ct = max_sample_ct * local_covar_ct * sizeof(double);
      uintptr_t per_alt_allele_xalloc_byte_ct = sizeof(LinearAuxResult);
      if (beta_se_multiallelic_fused) {
        per_variant_xalloc_byte_ct += 2 * max_reported_test_ct * subbatch_size * sizeof(double);
      } else {
        per_alt_allele_xalloc_byte_ct += 2 * max_reported_test_ct * subbatch_size * sizeof(double);
      }

      uintptr_t bytes_avail = bigstack_left();
      if (bytes_avail < cstream_alloc_size * subbatch_size) {
        continue;
      }
      bytes_avail -= cstream_alloc_size * subbatch_size;
      if (!PgenMtLoadInit(variant_include, max_sample_ct, variant_ct, bytes_avail, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, per_alt_allele_xalloc_byte_ct, pgfip, &calc_thread_ct, &common->genovecs, max_extra_allele_ct? (&common->thread_mhc) : nullptr, nullptr, nullptr, dosage_is_present? (&common->dosage_presents) : nullptr, dosage_is_present? (&common->dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, &max_alt_allele_block_size, main_loadbufs, &common->pgr_ptrs, &common->read_variant_uidx_starts)) {
        break;
      }
    }
    g_failed_alloc_attempt_size = 0;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto GlmLinearBatch_ret_NOMEM;
    }
    LinearAuxResult* linear_block_aux_bufs[2];
    double* block_beta_se_bufs[2];

    for (uint32_t uii = 0; uii != 2; ++uii) {
      if (unlikely(BIGSTACK_ALLOC_X(LinearAuxResult, max_alt_allele_block_size, &(linear_block_aux_bufs[uii])))) {
        // shouldn't be possible for these to fail?
        goto GlmLinearBatch_ret_NOMEM;
      }
      if (beta_se_multiallelic_fused) {
        if (unlikely(bigstack_alloc_d(read_block_size * (2 * k1LU) * max_reported_test_ct * subbatch_size, &(block_beta_se_bufs[uii])))) {
          goto GlmLinearBatch_ret_NOMEM;
        }
      } else {
        if (unlikely(bigstack_alloc_d(max_alt_allele_block_size * 2 * max_reported_test_ct * subbatch_size, &(block_beta_se_bufs[uii])))) {
          goto GlmLinearBatch_ret_NOMEM;
        }
      }
      if (local_covar_ct) {
        if (unlikely(bigstack_alloc_d(read_block_size * max_sample_ct * local_covar_ct, &(ctx->local_covars_vcmaj_d[uii])))) {
          goto GlmLinearBatch_ret_NOMEM;
        }
      } else {
        ctx->local_covars_vcmaj_d[uii] = nullptr;
      }
    }

    common->workspace_bufs = S_CAST(unsigned char**, bigstack_alloc_raw_rd(calc_thread_ct * sizeof(intptr_t)));
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      common->workspace_bufs[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(workspace_alloc));
    }
    common->err_info = (~0LLU) << 32;
    SetThreadFuncAndData(GlmLinearSubbatchThread, ctx, &tg);

    const uint32_t subbatch_ct = 1 + (batch_size - 1) / subbatch_size;

    const uint32_t ref_col = glm_cols & kfGlmColRef;
    const uint32_t alt1_col = glm_cols & kfGlmColAlt1;
    const uint32_t alt_col = glm_cols & kfGlmColAlt;
    const uint32_t ax_col = glm_cols & kfGlmColAx;
    const uint32_t a1_ct_col = glm_cols & kfGlmColA1count;
    const uint32_t tot_allele_col = glm_cols & kfGlmColTotallele;
    const uint32_t a1_freq_col = glm_cols & kfGlmColA1freq;
    const uint32_t mach_r2_col = glm_cols & kfGlmColMachR2;
    const uint32_t nobs_col = glm_cols & kfGlmColNobs;
    const uint32_t beta_col = glm_cols & (kfGlmColBeta | kfGlmColOrbeta);
    const uint32_t se_col = glm_cols & kfGlmColSe;
    const uint32_t ci_col = (ci_size != 0.0) && (glm_cols & kfGlmColCi);
    const uint32_t t_col = glm_cols & kfGlmColTz;
    const uint32_t p_col = glm_cols & kfGlmColP;
    const uint32_t err_col = glm_cols & kfGlmColErr;
    double ci_zt = 0.0;
    if (ci_col) {
      ci_zt = QuantileToZscore((ci_size + 1.0) * 0.5);
    }
    for (uint32_t subbatch_idx = 0; subbatch_idx != subbatch_ct; ++subbatch_idx) {
      const uint32_t pheno_uidx_start = IdxToUidxBasic(pheno_batch, subbatch_idx * subbatch_size);
      if (subbatch_idx == subbatch_ct - 1) {
        subbatch_size = batch_size - subbatch_idx * subbatch_size;
      }
      ctx->subbatch_size = subbatch_size;
      uint32_t pheno_uidx = pheno_uidx_start;
      for (uint32_t fidx = 0; fidx != subbatch_size; ++fidx, ++pheno_uidx) {
        pheno_uidx = AdvTo1Bit(pheno_batch, pheno_uidx);
        const char* cur_pheno_name = &(pheno_names[pheno_uidx * max_pheno_name_blen]);
        char* outname_end2 = strcpya(&(outname_end[1]), cur_pheno_name);
        outname_end2 = strcpya_k(outname_end2, ".glm.linear");
        if (output_zst) {
          snprintf(outname_end2, 22, ".zst");
        } else {
          *outname_end2 = '\0';
        }

        // forced-singlethreaded
        char* cswritep;
        reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &(css_arr[fidx]), &cswritep);
        if (unlikely(reterr)) {
          goto GlmLinearBatch_ret_1;
        }
        *cswritep++ = '#';
        if (chr_col) {
          cswritep = strcpya_k(cswritep, "CHROM\t");
        }
        if (variant_bps) {
          cswritep = strcpya_k(cswritep, "POS\t");
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
        if (a1_ct_col) {
          cswritep = strcpya_k(cswritep, "\tA1_CT");
        }
        if (tot_allele_col) {
          cswritep = strcpya_k(cswritep, "\tALLELE_CT");
        }
        if (a1_freq_col) {
          cswritep = strcpya_k(cswritep, "\tA1_FREQ");
        }
        if (mach_r2_col) {
          cswritep = strcpya_k(cswritep, "\tMACH_R2");
        }
        if (test_col) {
          cswritep = strcpya_k(cswritep, "\tTEST");
        }
        if (nobs_col) {
          cswritep = strcpya_k(cswritep, "\tOBS_CT");
        }
        if (beta_col) {
          cswritep = strcpya_k(cswritep, "\tBETA");
        }
        if (se_col) {
          cswritep = strcpya_k(cswritep, "\tSE");
        }
        if (ci_col) {
          cswritep = strcpya_k(cswritep, "\tL");
          cswritep = dtoa_g(ci_size * 100, cswritep);
          cswritep = strcpya_k(cswritep, "\tU");
          cswritep = dtoa_g(ci_size * 100, cswritep);
        }
        if (t_col) {
          if (!constraint_ct) {
            cswritep = strcpya_k(cswritep, "\tT_STAT");
          } else {
            // F-statistic for joint tests.
            cswritep = strcpya_k(cswritep, "\tT_OR_F_STAT");
          }
        }
        if (p_col) {
          if (report_neglog10p) {
            cswritep = strcpya_k(cswritep, "\tLOG10_P");
          } else {
            cswritep = strcpya_k(cswritep, "\tP");
          }
        }
        if (err_col) {
          cswritep = strcpya_k(cswritep, "\tERRCODE");
        }
        AppendBinaryEoln(&cswritep);
        cswritep_arr[fidx] = cswritep;

        FillPhenoAndXtY(sample_include, pheno_cols[pheno_uidx].data.qt, ctx->covars_cmaj_d, sample_ct, domdev_present_p1, covar_ct, common->nm_precomp? (&(common->nm_precomp->xt_y_image[fidx * (1 + domdev_present_p1 + covar_ct)])) : nullptr, &(ctx->pheno_d[fidx * sample_ct]));
        if (sample_ct_x) {
          FillPhenoAndXtY(common->sample_include_x, pheno_cols[pheno_uidx].data.qt, ctx->covars_cmaj_x_d, sample_ct_x, domdev_present_p1, covar_ct_x, common->nm_precomp_x? (&(common->nm_precomp_x->xt_y_image[fidx * (1 + domdev_present_p1 + covar_ct_x)])) : nullptr, &(ctx->pheno_x_d[fidx * sample_ct_x]));
        }
        if (sample_ct_y) {
          FillPhenoAndXtY(common->sample_include_y, pheno_cols[pheno_uidx].data.qt, ctx->covars_cmaj_y_d, sample_ct_y, domdev_present_p1, covar_ct_y, common->nm_precomp_y? (&(common->nm_precomp_y->xt_y_image[fidx * (1 + domdev_present_p1 + covar_ct_y)])) : nullptr, &(ctx->pheno_y_d[fidx * sample_ct_y]));
        }
      }

      // Main workflow:
      // 1. Set n=0, load/skip block 0
      //
      // 2. Spawn threads processing block n
      // 3. If n>0, write results for block (n-1)
      // 4. Increment n by 1
      // 5. Load/skip block n unless eof
      // 6. Join threads
      // 7. Goto step 2 unless eof
      //
      // 8, Write results for last block
      uintptr_t write_variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      uint32_t parity = 0;
      uint32_t read_block_idx = 0;
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t suppress_mach_r2 = 0;

      uint32_t cur_biallelic_reported_test_ct = 0;
      uint32_t primary_reported_test_idx = include_intercept;
      uint32_t cur_biallelic_predictor_ct = 0;
      uint32_t cur_constraint_ct = 0;

      const char* const* cur_test_names = nullptr;
      uint32_t prev_block_variant_ct = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uint32_t allele_ct = 2;
      uint32_t omitted_allele_idx = 0;
      if (subbatch_size > 1) {
        logprintfww5("--glm linear regression on quantitative phenotypes #%u-%u: ", completed_pheno_ct + 1, completed_pheno_ct + subbatch_size);
      } else {
        logprintfww5("--glm linear regression on phenotype '%s': ", &(pheno_names[pheno_uidx_start * max_pheno_name_blen]));
      }
      fputs("0%", stdout);
      fflush(stdout);
      // bugfix (12 May 2019): need to reinitialize this
      pgfip->block_base = main_loadbufs[0];
      ReinitThreads(&tg);
      for (uint32_t variant_idx = 0; ; ) {
        const uint32_t cur_block_variant_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
        if (unlikely(reterr)) {
          goto GlmLinearBatch_ret_PGR_FAIL;
        }
        if (local_covar_ct && cur_block_variant_ct) {
          const uint32_t uidx_start = read_block_idx * read_block_size;
          const uint32_t uidx_end = MINV(raw_variant_ct, uidx_start + read_block_size);
          if (local_variant_include) {
            reterr = ReadLocalCovarBlock(common, local_sample_uidx_order, local_variant_include, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, local_covar_txsp, &local_line_idx, &local_xy, nullptr, ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
          } else {
            double* prev_local_covar_row_d = nullptr;
            if (variant_idx) {
              prev_local_covar_row_d = &(ctx->local_covars_vcmaj_d[1 - parity][S_CAST(uintptr_t, read_block_size - 1) * max_sample_ct * local_covar_ct]);
            }
            reterr = ReadRfmix2Block(common, variant_bps, local_sample_uidx_order, nullptr, prev_local_covar_row_d, uidx_start, uidx_end, cur_block_variant_ct, local_sample_ct, glm_info_ptr->local_cat_ct, glm_info_ptr->local_chrom_col, glm_info_ptr->local_bp_col, glm_info_ptr->local_first_covar_col, local_covar_txsp, &local_line_iter, &local_line_idx, &local_prev_chr_code, &local_chr_code, &local_bp, &local_skip_chr, nullptr, ctx->local_covars_vcmaj_d[parity], local_sample_idx_order);
          }
          if (unlikely(reterr)) {
            goto GlmLinearBatch_ret_1;
          }
        }
        if (variant_idx) {
          JoinThreads(&tg);
          reterr = S_CAST(PglErr, common->err_info);
          if (unlikely(reterr)) {
            goto GlmLinearBatch_ret_PGR_FAIL;
          }
        }
        if (!IsLastBlock(&tg)) {
          common->cur_block_variant_ct = cur_block_variant_ct;
          const uint32_t uidx_start = read_block_idx * read_block_size;
          ComputeUidxStartPartition(variant_include, cur_block_variant_ct, calc_thread_ct, uidx_start, common->read_variant_uidx_starts);
          PgrCopyBaseAndOffset(pgfip, calc_thread_ct, common->pgr_ptrs);
          ctx->block_aux = linear_block_aux_bufs[parity];
          common->block_beta_se = block_beta_se_bufs[parity];
          if (variant_idx + cur_block_variant_ct == variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto GlmLinearBatch_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (variant_idx) {
          // write *previous* block results
          const double* beta_se_iter = block_beta_se_bufs[parity];
          const LinearAuxResult* cur_block_aux = linear_block_aux_bufs[parity];
          uintptr_t allele_bidx = 0;
          for (uint32_t variant_bidx = 0; variant_bidx != prev_block_variant_ct; ++variant_bidx) {
            const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
            if (write_variant_uidx >= chr_end) {
              do {
                ++chr_fo_idx;
                chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
              } while (write_variant_uidx >= chr_end);
              const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
              if ((chr_idx == x_code) && sample_ct_x) {
                cur_biallelic_reported_test_ct = biallelic_reported_test_ct_x;
                cur_biallelic_predictor_ct = biallelic_predictor_ct_x;
                cur_constraint_ct = constraint_ct_x;
                cur_test_names = test_names_x;
              } else if ((chr_idx == y_code) && sample_ct_y) {
                cur_biallelic_reported_test_ct = biallelic_reported_test_ct_y;
                cur_biallelic_predictor_ct = biallelic_predictor_ct_y;
                cur_constraint_ct = constraint_ct_y;
                cur_test_names = test_names_y;
              } else {
                cur_biallelic_reported_test_ct = biallelic_reported_test_ct;
                cur_biallelic_predictor_ct = biallelic_predictor_ct;
                cur_constraint_ct = constraint_ct;
                cur_test_names = test_names;
              }
              suppress_mach_r2 = (chr_idx == x_code) || (chr_idx == mt_code);
              if (cur_constraint_ct) {
                primary_reported_test_idx = cur_biallelic_reported_test_ct - 1;
              }
              if (chr_col) {
                char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
                *chr_name_end = '\t';
                chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
              }
            }
            uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
            if (allele_idx_offsets) {
              allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
              allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offsets[write_variant_uidx];
            }
            const uint32_t allele_ct_m1 = allele_ct - 1;
            const uint32_t extra_allele_ct = allele_ct - 2;
            if (omitted_alleles) {
              omitted_allele_idx = omitted_alleles[write_variant_uidx];
            }
            const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
            for (uint32_t fidx = 0; fidx != subbatch_size; ++fidx) {
              char* cswritep = cswritep_arr[fidx];
              uint32_t a1_allele_idx = 0;
              uintptr_t allele_bidx_tmp = allele_bidx;
              // bugfix (20 Mar 2020): In multiallelic unfused case, fidx is
              // the *inner* index, allele index is on the outside.
              const double* beta_se_iter2 = beta_se_iter;
              for (uint32_t nonomitted_allele_idx = 0; nonomitted_allele_idx != allele_ct_m1; ++nonomitted_allele_idx, ++a1_allele_idx) {
                if (beta_se_multiallelic_fused) {
                  if (!nonomitted_allele_idx) {
                    primary_reported_test_idx = include_intercept;
                  } else {
                    primary_reported_test_idx = cur_biallelic_reported_test_ct + nonomitted_allele_idx - 1;
                  }
                }
                if (nonomitted_allele_idx == omitted_allele_idx) {
                  ++a1_allele_idx;
                }
                const double primary_beta = beta_se_iter2[primary_reported_test_idx * 2];
                const double primary_se = beta_se_iter2[primary_reported_test_idx * 2 + 1];
                const uint32_t allele_is_valid = (primary_se != -9.0);
                {
                  const LinearAuxResult* auxp = &(cur_block_aux[allele_bidx_tmp]);
                  if (ln_pfilter <= 0.0) {
                    if (!allele_is_valid) {
                      goto GlmLinearBatch_allele_iterate;
                    }
                    double primary_ln_pval;
                    if (!cur_constraint_ct) {
                      if (primary_beta == 0.0) {
                        primary_ln_pval = 0.0;
                      } else if (primary_se == 0.0) {
                        primary_ln_pval = -DBL_MAX;
                      } else {
                        const double primary_tstat = primary_beta / primary_se;
                        primary_ln_pval = TstatToLnP(primary_tstat, auxp->sample_obs_ct - cur_biallelic_predictor_ct - extra_allele_ct);
                      }
                    } else {
                      primary_ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                    }
                    if (primary_ln_pval > ln_pfilter) {
                      goto GlmLinearBatch_allele_iterate;
                    }
                  }
                  uint32_t inner_reported_test_ct = cur_biallelic_reported_test_ct;
                  if (extra_allele_ct) {
                    if (beta_se_multiallelic_fused) {
                      if (!nonomitted_allele_idx) {
                        inner_reported_test_ct = 1 + include_intercept;
                      } else if (nonomitted_allele_idx == extra_allele_ct) {
                        inner_reported_test_ct -= include_intercept;
                      } else {
                        inner_reported_test_ct = 1;
                      }
                    } else if (!hide_covar) {
                      inner_reported_test_ct += extra_allele_ct;
                    }
                  }
                  // possible todo: make number-to-string operations, strlen(),
                  // etc. happen only once per variant.
                  for (uint32_t allele_test_idx = 0; allele_test_idx != inner_reported_test_ct; ++allele_test_idx) {
                    uint32_t test_idx = allele_test_idx;
                    if (beta_se_multiallelic_fused && nonomitted_allele_idx) {
                      if (!allele_test_idx) {
                        test_idx = primary_reported_test_idx;
                      } else {
                        ++test_idx;
                      }
                    }
                    if (chr_col) {
                      cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
                    }
                    if (variant_bps) {
                      cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
                    }
                    cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
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
                        if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                          // might not need this assignment, but play it safe
                          // for now
                          cswritep_arr[fidx] = cswritep;
                          goto GlmLinearBatch_ret_WRITE_FAIL;
                        }
                        cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                      }
                      --cswritep;
                    }
                    *cswritep++ = '\t';
                    const uint32_t multi_a1 = extra_allele_ct && beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx);
                    if (multi_a1) {
                      for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                        if (allele_idx == omitted_allele_idx) {
                          continue;
                        }
                        if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                          cswritep_arr[fidx] = cswritep;
                          goto GlmLinearBatch_ret_WRITE_FAIL;
                        }
                        cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                      }
                      --cswritep;
                    } else {
                      cswritep = strcpya(cswritep, cur_alleles[a1_allele_idx]);
                    }
                    if (ax_col) {
                      *cswritep++ = '\t';
                      if (beta_se_multiallelic_fused && (test_idx != primary_reported_test_idx)) {
                        if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                          cswritep_arr[fidx] = cswritep;
                          goto GlmLinearBatch_ret_WRITE_FAIL;
                        }
                        cswritep = strcpya(cswritep, cur_alleles[omitted_allele_idx]);
                      } else {
                        for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
                          if (allele_idx == a1_allele_idx) {
                            continue;
                          }
                          if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                            cswritep_arr[fidx] = cswritep;
                            goto GlmLinearBatch_ret_WRITE_FAIL;
                          }
                          cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
                        }
                        --cswritep;
                      }
                    }
                    if (a1_ct_col) {
                      *cswritep++ = '\t';
                      if (!multi_a1) {
                        cswritep = dtoa_g(auxp->a1_dosage, cswritep);
                      } else {
                        cswritep = strcpya_k(cswritep, "NA");
                      }
                    }
                    if (tot_allele_col) {
                      *cswritep++ = '\t';
                      cswritep = u32toa(auxp->allele_obs_ct, cswritep);
                    }
                    if (a1_freq_col) {
                      *cswritep++ = '\t';
                      if (!multi_a1) {
                        cswritep = dtoa_g(auxp->a1_dosage / S_CAST(double, auxp->allele_obs_ct), cswritep);
                      } else {
                        cswritep = strcpya_k(cswritep, "NA");
                      }
                    }
                    if (mach_r2_col) {
                      *cswritep++ = '\t';
                      if (!suppress_mach_r2) {
                        cswritep = dtoa_g(auxp->mach_r2, cswritep);
                      } else {
                        cswritep = strcpya_k(cswritep, "NA");
                      }
                    }
                    if (test_col) {
                      *cswritep++ = '\t';
                      if (test_idx < cur_biallelic_reported_test_ct) {
                        cswritep = strcpya(cswritep, cur_test_names[test_idx]);
                      } else {
                        // always use basic dosage for untested alleles
                        cswritep = strcpya_k(cswritep, "ADD");
                        if (!beta_se_multiallelic_fused) {
                          // extra alt allele covariate.
                          uint32_t test_xallele_idx = test_idx - cur_biallelic_reported_test_ct;
                          if (omitted_allele_idx < a1_allele_idx) {
                            test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                          }
                          test_xallele_idx = test_xallele_idx + (test_xallele_idx >= a1_allele_idx);
                          if (a1_allele_idx < omitted_allele_idx) {
                            test_xallele_idx = test_xallele_idx + (test_xallele_idx >= omitted_allele_idx);
                          }
                          if (!test_xallele_idx) {
                            cswritep = strcpya_k(cswritep, "_REF");
                          } else {
                            cswritep = strcpya_k(cswritep, "_ALT");
                            cswritep = u32toa(test_xallele_idx, cswritep);
                          }
                        }
                      }
                    }
                    if (nobs_col) {
                      *cswritep++ = '\t';
                      cswritep = u32toa(auxp->sample_obs_ct, cswritep);
                    }
                    double ln_pval = kLnPvalError;
                    double tstat = 0.0;
                    uint32_t test_is_valid;
                    if ((!cur_constraint_ct) || (test_idx != primary_reported_test_idx)) {
                      double beta = beta_se_iter2[2 * test_idx];
                      double se = beta_se_iter2[2 * test_idx + 1];
                      test_is_valid = (se != -9.0);
                      if (test_is_valid) {
                        tstat = beta / se;
                        ln_pval = TstatToLnP(tstat, auxp->sample_obs_ct - cur_biallelic_predictor_ct - extra_allele_ct);
                      }
                      if (beta_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          cswritep = dtoa_g(beta, cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA");
                        }
                      }
                      if (se_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          cswritep = dtoa_g(se, cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA");
                        }
                      }
                      if (ci_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          const double ci_halfwidth = ci_zt * se;
                          cswritep = dtoa_g(beta - ci_halfwidth, cswritep);
                          *cswritep++ = '\t';
                          cswritep = dtoa_g(beta + ci_halfwidth, cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA\tNA");
                        }
                      }
                      if (t_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          cswritep = dtoa_g(tstat, cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA");
                        }
                      }
                    } else {
                      // joint test
                      test_is_valid = allele_is_valid;
                      if (beta_col) {
                        cswritep = strcpya_k(cswritep, "\tNA");
                      }
                      if (se_col) {
                        cswritep = strcpya_k(cswritep, "\tNA");
                      }
                      if (ci_col) {
                        cswritep = strcpya_k(cswritep, "\tNA\tNA");
                      }
                      if (t_col) {
                        *cswritep++ = '\t';
                        if (test_is_valid) {
                          cswritep = dtoa_g(primary_se / u31tod(cur_constraint_ct), cswritep);
                        } else {
                          cswritep = strcpya_k(cswritep, "NA");
                        }
                      }
                      // could avoid recomputing
                      if (test_is_valid) {
                        ln_pval = FstatToLnP(primary_se / u31tod(cur_constraint_ct), cur_constraint_ct, auxp->sample_obs_ct);
                      }
                    }
                    if (p_col) {
                      *cswritep++ = '\t';
                      if (test_is_valid) {
                        if (report_neglog10p) {
                          const double reported_val = (-kRecipLn10) * ln_pval;
                          cswritep = dtoa_g(reported_val, cswritep);
                        } else {
                          const double reported_ln = MAXV(ln_pval, output_min_ln);
                          cswritep = lntoa_g(reported_ln, cswritep);
                        }
                      } else {
                        cswritep = strcpya_k(cswritep, "NA");
                      }
                    }
                    if (err_col) {
                      *cswritep++ = '\t';
                      if (test_is_valid) {
                        *cswritep++ = '.';
                      } else {
                        uint64_t glm_errcode;
                        memcpy(&glm_errcode, &(beta_se_iter2[2 * test_idx]), 8);
                        cswritep = AppendGlmErrstr(glm_errcode, cswritep);
                      }
                    }
                    AppendBinaryEoln(&cswritep);
                    if (unlikely(Cswrite(&(css_arr[fidx]), &cswritep))) {
                      cswritep_arr[fidx] = cswritep;
                      goto GlmLinearBatch_ret_WRITE_FAIL;
                    }
                  }
                }
              GlmLinearBatch_allele_iterate:
                ++allele_bidx_tmp;
                if (!beta_se_multiallelic_fused) {
                  beta_se_iter2 = &(beta_se_iter2[subbatch_size * (2 * k1LU) * max_reported_test_ct]);
                }
              }  // for nonomitted_allele_idx
              beta_se_iter = &(beta_se_iter[2 * max_reported_test_ct]);
              cswritep_arr[fidx] = cswritep;
            }  // for fidx
            if (!beta_se_multiallelic_fused) {
              beta_se_iter = &(beta_se_iter[extra_allele_ct * subbatch_size * (2 * k1LU) * max_reported_test_ct]);
            }
            allele_bidx += allele_ct_m1;
          }  // for variant_bidx
        }
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
        prev_block_variant_ct = cur_block_variant_ct;
        variant_idx += cur_block_variant_ct;
        // crucially, this is independent of the PgenReader block_base
        // pointers
        pgfip->block_base = main_loadbufs[parity];
      }
      for (uintptr_t fidx = 0; fidx != subbatch_size; ++fidx) {
        if (unlikely(CswriteCloseNull(&(css_arr[fidx]), cswritep_arr[fidx]))) {
          goto GlmLinearBatch_ret_WRITE_FAIL;
        }
      }
      if (pct > 10) {
        putc_unlocked('\b', stdout);
      }
      fputs("\b\b", stdout);
      logputs("done.\n");
      // bugfix (12 May 2019): added batch_size instead of subbatch_size here
      completed_pheno_ct += subbatch_size;
    }
    outname_end[1] = '\0';
    logprintfww("Results written to %s<phenotype name>.glm.linear%s .\n", outname, output_zst? ".zst" : "");
  }
  while (0) {
  GlmLinearBatch_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GlmLinearBatch_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--glm local-covar= file", local_covar_txsp);
    break;
  GlmLinearBatch_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  GlmLinearBatch_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GlmLinearBatch_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  GlmLinearBatch_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 GlmLinearBatch_ret_1:
  CleanupThreads(&tg);
  if (css_arr) {
    for (uintptr_t fidx = 0; fidx != subbatch_size; ++fidx) {
      CswriteCloseCond(&(css_arr[fidx]), cswritep_arr[fidx]);
    }
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

static const double kSexMaleToCovarD[2] = {2.0, 1.0};

void SexInteractionReshuffle(uint32_t first_interaction_pred_uidx, uint32_t raw_covar_ct, uint32_t domdev_present, uint32_t biallelic_raw_predictor_ctl, uintptr_t* __restrict parameters_or_tests, uintptr_t* __restrict parameter_subset_reshuffle_buf) {
  ZeroWArr(biallelic_raw_predictor_ctl, parameter_subset_reshuffle_buf);
  CopyBitarrRange(parameters_or_tests, 0, 0, first_interaction_pred_uidx - 1, parameter_subset_reshuffle_buf);
  const uint32_t raw_interaction_ct = raw_covar_ct * (domdev_present + 1);
  CopyBitarrRange(parameters_or_tests, first_interaction_pred_uidx - 1, first_interaction_pred_uidx, raw_interaction_ct, parameter_subset_reshuffle_buf);
  const uint32_t first_sex_parameter_idx = first_interaction_pred_uidx - 1 + raw_interaction_ct;
  if (IsSet(parameters_or_tests, first_sex_parameter_idx)) {
    SetBit(first_interaction_pred_uidx - 1, parameter_subset_reshuffle_buf);
  }
  if (IsSet(parameters_or_tests, first_sex_parameter_idx + 1)) {
    SetBit(first_sex_parameter_idx + 1, parameter_subset_reshuffle_buf);
  }
  if (domdev_present && IsSet(parameters_or_tests, first_sex_parameter_idx + 2)) {
    SetBit(first_sex_parameter_idx + 2, parameter_subset_reshuffle_buf);
  }
  memcpy(parameters_or_tests, parameter_subset_reshuffle_buf, biallelic_raw_predictor_ctl * sizeof(intptr_t));
}

PglErr GlmMain(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const AdjustInfo* adjust_info_ptr, const APerm* aperm_ptr, const char* local_covar_fname, const char* local_pvar_fname, const char* local_psam_fname, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t orig_covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t xchr_model, double ci_size, double vif_thresh, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  TextStream local_covar_txs;
  TokenStream tks;
  PreinitTextStream(&local_covar_txs);
  PreinitTokenStream(&tks);
  GlmCtx common;
  GlmLogisticCtx logistic_ctx;
  GlmLinearCtx linear_ctx;
  logistic_ctx.common = &common;
  linear_ctx.common = &common;
  {
    if (unlikely(!pheno_ct)) {
      logerrputs("Error: No phenotypes loaded.\n");
      goto GlmMain_ret_INCONSISTENT_INPUT;
    }
    if (unlikely(orig_sample_ct < 2)) {
      logerrputs("Error: --glm requires at least two samples.\n");
      goto GlmMain_ret_DEGENERATE_DATA;
    }
    assert(orig_variant_ct);
    // common linear/logistic initialization
    const GlmFlags glm_flags = glm_info_ptr->flags;
    const uintptr_t* early_variant_include = orig_variant_include;
    uint32_t* local_sample_uidx_order = nullptr;
    uintptr_t* local_variant_include = nullptr;
    uint32_t variant_ct = orig_variant_ct;
    uint32_t local_sample_ct = 0;
    uint32_t local_variant_ctl = 0;
    uint32_t local_covar_ct = 0;
    if (local_covar_fname) {
      reterr = GlmLocalOpen(local_covar_fname, local_pvar_fname, local_psam_fname, siip->sample_ids, cip, variant_bps, variant_ids, glm_info_ptr, raw_sample_ct, siip->max_sample_id_blen, raw_variant_ct, &orig_sample_include, &sex_nm, &sex_male, &early_variant_include, &orig_sample_ct, &variant_ct, &local_covar_txs, &local_sample_uidx_order, &local_variant_include, &local_sample_ct, &local_variant_ctl, &local_covar_ct);
      if (unlikely(reterr)) {
        goto GlmMain_ret_1;
      }
    }

    common.glm_flags = glm_flags;
    common.dosage_presents = nullptr;
    common.dosage_mains = nullptr;
    const uint32_t output_zst = (glm_flags / kfGlmZs) & 1;
    const uint32_t perm_adapt = (glm_flags / kfGlmPerm) & 1;
    const uint32_t perms_total = perm_adapt? aperm_ptr->max : glm_info_ptr->mperm_ct;
    // <output prefix>.<pheno name>.glm.logistic.hybrid{,.perm,.mperm}[.zst]
    uint32_t pheno_name_blen_capacity = kPglFnamesize - 21 - (4 * output_zst) - S_CAST(uintptr_t, outname_end - outname);
    if (perms_total) {
      pheno_name_blen_capacity -= 6 - perm_adapt;
    }
    if (unlikely(max_pheno_name_blen > pheno_name_blen_capacity)) {
      logerrputs("Error: Phenotype name and/or --out parameter too long.\n");
      goto GlmMain_ret_INCONSISTENT_INPUT;
    }
    *outname_end = '.';
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;

    // synthetic categorical covariate name could be ~twice max ID length?
    const uintptr_t overflow_buf_size = kCompressStreamBlock + 2 * kMaxIdSlen + max_chr_blen + kMaxIdSlen + 1024 + 2 * max_allele_slen;

    uintptr_t* cur_sample_include;
    if (unlikely(
            bigstack_alloc_w(raw_sample_ctl, &cur_sample_include) ||
            bigstack_alloc_u32(raw_sample_ctl, &common.sample_include_cumulative_popcounts))) {
      goto GlmMain_ret_NOMEM;
    }
    common.sample_include = cur_sample_include;
    common.cip = cip;
    common.allele_idx_offsets = allele_idx_offsets;

    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uint32_t max_variant_ct = variant_ct;

    uint32_t x_start;
    uint32_t x_end;
    GetXymtStartAndEnd(cip, kChrOffsetX, &x_start, &x_end);
    uint32_t y_start;
    uint32_t y_end;
    GetXymtStartAndEnd(cip, kChrOffsetY, &y_start, &y_end);

    uintptr_t* sex_male_collapsed_buf = nullptr;
    uint32_t x_code;
    uint32_t variant_ct_x = 0;
    uint32_t variant_ct_y = 0;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t sex_nm_ct = PopcountWords(sex_nm, raw_sample_ctl);
    const uint32_t male_ct = PopcountWords(sex_male, raw_sample_ctl);
    uint32_t add_sex_covar = !(glm_flags & kfGlmNoXSex);
    if (add_sex_covar && ((!male_ct) || (male_ct == sex_nm_ct))) {
      add_sex_covar = 0;
    }
    uintptr_t* cur_sample_include_y_buf = nullptr;
    if (domdev_present || (glm_flags & (kfGlmDominant | kfGlmRecessive))) {
      // dominant/recessive/genotypic/hethom suppress all chromosomes which
      // aren't fully diploid.  (could throw in a hack to permit chrX if
      // all samples are female?  i.e. synthesize a ChrInfo where
      // xymt_codes[0] is UINT32_MAXM1 and haploid_mask X bit is cleared)
      uintptr_t* variant_include_nohap = nullptr;
      const uint32_t chr_ct = cip->chr_ct;
      uint32_t removed_variant_ct = 0;
      for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        if (IsSet(cip->haploid_mask, chr_idx)) {
          const uint32_t variant_uidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
          const uint32_t variant_uidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          const uint32_t cur_chr_variant_ct = PopcountBitRange(early_variant_include, variant_uidx_start, variant_uidx_end);
          if (cur_chr_variant_ct) {
            if (!removed_variant_ct) {
              // no main-loop logic for excluding all haploid chromosomes, so
              // make a full copy of early_variant_include and throw away our
              // reference to the original
              if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include_nohap))) {
                goto GlmMain_ret_NOMEM;
              }
              memcpy(variant_include_nohap, early_variant_include, raw_variant_ctl * sizeof(intptr_t));
            }
            ClearBitsNz(variant_uidx_start, variant_uidx_end, variant_include_nohap);
            removed_variant_ct += cur_chr_variant_ct;
          }
        }
      }
      if (removed_variant_ct) {
        if (unlikely(variant_ct == removed_variant_ct)) {
          logerrputs("Error: No variants remaining for --glm ('dominant', 'recessive', 'genotypic',\nand 'hethom' only operate on diploid data).\n");
          goto GlmMain_ret_DEGENERATE_DATA;
        }
        variant_ct -= removed_variant_ct;
        early_variant_include = variant_include_nohap;
        max_variant_ct = variant_ct;
      }
    } else {
      if (XymtExists(cip, kChrOffsetX, &x_code)) {
        variant_ct_x = CountChrVariantsUnsafe(early_variant_include, cip, x_code);
        // --xchr-model 0 now only suppresses chrX.
        if (xchr_model) {
          if (variant_ct_x) {
            if (unlikely(bigstack_alloc_w(BitCtToWordCt(orig_sample_ct), &sex_male_collapsed_buf))) {
              goto GlmMain_ret_NOMEM;
            }
          }
        } else {
          max_variant_ct -= variant_ct_x;
          if (unlikely(!max_variant_ct)) {
            logerrputs("Error: No variants remaining for --glm, due to --xchr-model 0.\n");
            goto GlmMain_ret_DEGENERATE_DATA;
          }
        }
      }
      uint32_t y_code;
      if (XymtExists(cip, kChrOffsetY, &y_code)) {
        variant_ct_y = CountChrVariantsUnsafe(early_variant_include, cip, y_code);
        if (variant_ct_y) {
          if (!male_ct) {
            logputs("--glm: Skipping chrY since there are no males.\n");
            max_variant_ct -= variant_ct_y;
            if (unlikely(!max_variant_ct)) {
              logerrputs("Error: No variants remaining for --glm.\n");
              goto GlmMain_ret_DEGENERATE_DATA;
            }
          } else if (male_ct < orig_sample_ct) {
            // may as well check for only-chrY special case
            if (max_variant_ct != variant_ct_y) {
              if (unlikely(bigstack_alloc_w(raw_sample_ctl, &cur_sample_include_y_buf))) {
                // covar_include_y allocation postponed since raw_covar_ct not
                // yet known
                goto GlmMain_ret_NOMEM;
              }
            } else {
              orig_sample_include = sex_male;
              orig_sample_ct = male_ct;
            }
          }
        }
      }
    }
    if (add_sex_covar && (!variant_ct_x) && (!(glm_flags & kfGlmSex))) {
      add_sex_covar = 0;
    }
    common.sex_male_collapsed = sex_male_collapsed_buf;
    common.omitted_alleles = (glm_flags & kfGlmOmitRef)? nullptr : maj_alleles;
    uint32_t raw_covar_ct = orig_covar_ct + local_covar_ct;
    if (glm_info_ptr->condition_varname || glm_info_ptr->condition_list_fname || local_covar_ct || add_sex_covar) {
      uint32_t condition_ct = 0;
      PhenoCol* new_covar_cols;
      char* new_covar_names;
      uintptr_t new_max_covar_name_blen = max_covar_name_blen;
      if (add_sex_covar && (new_max_covar_name_blen < 4)) {
        new_max_covar_name_blen = 4;
      }
      if (local_covar_ct && (new_max_covar_name_blen < 6 + UintSlen(local_covar_ct + 1))) {
        new_max_covar_name_blen = 6 + UintSlen(local_covar_ct + 1);
      }
      if (glm_info_ptr->condition_varname || glm_info_ptr->condition_list_fname) {
        assert(g_bigstack_end == bigstack_end_mark);
        const uint32_t condition_multiallelic = (glm_flags / kfGlmConditionMultiallelic) & 1;
        if (condition_multiallelic) {
          logerrputs("Error: --condition[-list] 'multiallelic' implementation is under development.\n");
          reterr = kPglRetNotYetSupported;
          goto GlmMain_ret_1;
        }
        // reserve space for condition-list worst case (roughly sqrt(2^31)),
        // since that's relatively small
        const uint32_t condition_ct_max = 46338;
        uint32_t* condition_uidxs;
        if (unlikely(bigstack_end_alloc_u32(condition_ct_max, &condition_uidxs))) {
          goto GlmMain_ret_NOMEM;
        }
        if (glm_info_ptr->condition_varname) {
          int32_t ii = GetVariantUidxWithoutHtable(glm_info_ptr->condition_varname, variant_ids, orig_variant_include, orig_variant_ct);
          if (ii >= 0) {
            condition_uidxs[0] = ii;
            condition_ct = 1;
            const uint32_t condition_blen = strlen(glm_info_ptr->condition_varname) + 1;
            // drop "CSNP" column name for sanity's sake
            if (new_max_covar_name_blen < condition_blen) {
              new_max_covar_name_blen = condition_blen;
            }
          } else {
            if (unlikely(ii == -2)) {
              logerrprintfww("Error: Duplicate --condition variant ID '%s'.\n", glm_info_ptr->condition_varname);
              goto GlmMain_ret_INVALID_CMDLINE;
            }
            logerrprintfww("Warning: --condition variant ID '%s' not found.\n", glm_info_ptr->condition_varname);
          }
        } else {
          // 1. (re)construct variant ID hash table
          uintptr_t* already_seen;
          if (unlikely(bigstack_calloc_w(raw_variant_ctl, &already_seen))) {
            goto GlmMain_ret_NOMEM;
          }
          reterr = InitTokenStream(glm_info_ptr->condition_list_fname, MAXV(max_thread_ct - 1, 1), &tks);
          if (unlikely(reterr)) {
            goto GlmMain_ret_TKSTREAM_FAIL;
          }
          uint32_t* variant_id_htable = nullptr;
          uint32_t variant_id_htable_size;
          reterr = AllocAndPopulateIdHtableMt(orig_variant_include, variant_ids, orig_variant_ct, 0, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
          if (unlikely(reterr)) {
            goto GlmMain_ret_1;
          }

          // 2. iterate through --condition-list file, make sure no IDs are
          //    duplicate in loaded fileset, warn about duplicates in
          //    --condition-list file
          uintptr_t skip_ct = 0;
          uintptr_t duplicate_ct = 0;
          while (1) {
            char* shard_boundaries[2];
            reterr = TksNext(&tks, 1, shard_boundaries);
            if (reterr) {
              break;
            }
            const char* shard_iter = shard_boundaries[0];
            const char* shard_end = shard_boundaries[1];
            while (1) {
              shard_iter = FirstPostspaceBounded(shard_iter, shard_end);
              if (shard_iter == shard_end) {
                break;
              }
              const char* token_end = CurTokenEnd(shard_iter);
              const uint32_t token_slen = token_end - shard_iter;
              uint32_t cur_variant_uidx = VariantIdDupflagHtableFind(shard_iter, variant_ids, variant_id_htable, token_slen, variant_id_htable_size, max_variant_id_slen);
              if (cur_variant_uidx >> 31) {
                if (unlikely(cur_variant_uidx != UINT32_MAX)) {
                  logerrprintfww("Error: --condition-list variant ID '%s' appears multiple times.\n", variant_ids[cur_variant_uidx & 0x7fffffff]);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                ++skip_ct;
              } else if (IsSet(already_seen, cur_variant_uidx)) {
                ++duplicate_ct;
              } else {
                if (unlikely(condition_ct == condition_ct_max)) {
                  logerrputs("Error: Too many --condition-list variant IDs.\n");
                  goto GlmMain_ret_MALFORMED_INPUT;
                }
                SetBit(cur_variant_uidx, already_seen);
                condition_uidxs[condition_ct++] = cur_variant_uidx;
                if (new_max_covar_name_blen <= token_slen) {
                  new_max_covar_name_blen = token_slen + 1;
                }
              }
              shard_iter = token_end;
            }
          }
          if (unlikely(reterr != kPglRetEof)) {
            goto GlmMain_ret_TKSTREAM_FAIL;
          }
          if (CleanupTokenStream3("--condition-list file", &tks, &reterr)) {
            goto GlmMain_ret_1;
          }
          if (skip_ct || duplicate_ct) {
            if (skip_ct && duplicate_ct) {
              logerrprintfww("Warning: %" PRIuPTR " --condition-list variant ID%s not found, and %" PRIuPTR " duplicate ID%s present.\n", skip_ct, (skip_ct == 1)? "" : "s", duplicate_ct, (duplicate_ct == 1)? "" : "s");
            } else if (skip_ct) {
              logerrprintf("Warning: %" PRIuPTR " --condition-list variant ID%s not found.\n", skip_ct, (skip_ct == 1)? "" : "s");
            } else {
              logerrprintf("Warning: %" PRIuPTR " duplicate --condition-list variant ID%s present.\n", duplicate_ct, (duplicate_ct == 1)? "" : "s");
            }
          }
          logprintf("--condition-list: %u variant ID%s loaded.\n", condition_ct, (condition_ct == 1)? "" : "s");

          // free hash table, duplicate tracker, TokenStream
          BigstackReset(already_seen);
        }
        raw_covar_ct += condition_ct;
        if (unlikely(
                BIGSTACK_ALLOC_X(PhenoCol, raw_covar_ct + add_sex_covar, &new_covar_cols) ||
                bigstack_alloc_c((raw_covar_ct + add_sex_covar) * new_max_covar_name_blen, &new_covar_names))) {
          goto GlmMain_ret_NOMEM;
        }
        if (condition_ct) {
          BigstackEndSet(condition_uidxs);
          uintptr_t* genovec;
          uintptr_t* dosage_present;
          Dosage* dosage_main;
          if (unlikely(
                  bigstack_end_alloc_w(NypCtToWordCt(raw_sample_ct), &genovec) ||
                  bigstack_end_alloc_w(raw_sample_ctl, &dosage_present) ||
                  bigstack_end_alloc_dosage(raw_sample_ct, &dosage_main))) {
            goto GlmMain_ret_NOMEM;
          }
          PgrSampleSubsetIndex null_pssi;
          PgrClearSampleSubsetIndex(simple_pgrp, &null_pssi);
          uint32_t allele_ct = 2;
          for (uint32_t condition_idx = 0; condition_idx != condition_ct; ++condition_idx) {
            const uint32_t cur_variant_uidx = condition_uidxs[condition_idx];
            if (allele_idx_offsets) {
              allele_ct = allele_idx_offsets[cur_variant_uidx + 1] - allele_idx_offsets[cur_variant_uidx];
              if ((allele_ct != 2) && (!condition_multiallelic)) {
                logerrputs("Error: --condition[-list] includes a multiallelic variant, but 'multiallelic'\nmodifier was not specified.\n");
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
            }
            uint32_t dosage_ct;
            reterr = PgrGetD(nullptr, null_pssi, raw_sample_ct, cur_variant_uidx, simple_pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
            if (unlikely(reterr)) {
              goto GlmMain_ret_PGR_FAIL;
            }
            // alpha 2 update: default to major allele, respect omit-ref
            if (common.omitted_alleles && common.omitted_alleles[cur_variant_uidx]) {
              GenovecInvertUnsafe(raw_sample_ct, genovec);
              if (dosage_ct) {
                BiallelicDosage16Invert(dosage_ct, dosage_main);
              }
            }
            PhenoCol* cur_covar_col = &(new_covar_cols[local_covar_ct + condition_idx]);
            uintptr_t* cur_nonmiss;
            double* cur_covar_vals;
            if (unlikely(
                    bigstack_alloc_w(raw_sample_ctl, &cur_nonmiss) ||
                    bigstack_alloc_d(raw_sample_ct, &cur_covar_vals))) {
              goto GlmMain_ret_NOMEM;
            }
            cur_covar_col->category_names = nullptr;
            cur_covar_col->nonmiss = cur_nonmiss;
            cur_covar_col->data.qt = cur_covar_vals;
            cur_covar_col->type_code = kPhenoDtypeQt;
            cur_covar_col->nonnull_category_ct = 0;
            GenoarrToNonmissing(genovec, raw_sample_ct, cur_nonmiss);
            GenoarrLookup16x8bx2(genovec, kSmallDoublePairs, raw_sample_ct, cur_covar_vals);
            if (dosage_ct) {
              uintptr_t sample_uidx_base = 0;
              uintptr_t cur_bits = dosage_present[0];
              for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
                const uintptr_t sample_uidx = BitIter1(dosage_present, &sample_uidx_base, &cur_bits);
                cur_covar_vals[sample_uidx] = kRecipDosageMid * u31tod(dosage_main[dosage_idx]);
              }
              BitvecOr(dosage_present, raw_sample_ctl, cur_nonmiss);
            }
            if (glm_flags & kfGlmConditionDominant) {
              for (uint32_t sample_uidx = 0; sample_uidx != raw_sample_ct; ++sample_uidx) {
                if (cur_covar_vals[sample_uidx] > 1.0) {
                  cur_covar_vals[sample_uidx] = 1.0;
                }
              }
            } else if (glm_flags & kfGlmConditionRecessive) {
              for (uint32_t sample_uidx = 0; sample_uidx != raw_sample_ct; ++sample_uidx) {
                double dxx = cur_covar_vals[sample_uidx];
                if (dxx <= 1.0) {
                  dxx = 0;
                } else {
                  dxx -= 1.0;
                }
                cur_covar_vals[sample_uidx] = dxx;
              }
            }
            // quasi-bugfix (21 Feb 2018): also should respect
            // "--xchr-model 1", and divide by 2 in haploid case.
            const uint32_t chr_idx = GetVariantChr(cip, cur_variant_uidx);
            if (IsSet(cip->haploid_mask, chr_idx)) {
              if (chr_idx == cip->xymt_codes[kChrOffsetX]) {
                if (xchr_model == 1) {
                  if (unlikely(glm_flags & (kfGlmConditionDominant | kfGlmConditionRecessive))) {
                    logerrputs("Error: --condition[-list] 'dominant'/'recessive' cannot be used with a chrX\nvariant when \"--xchr-model 1\" is in effect.\n");
                    goto GlmMain_ret_INCONSISTENT_INPUT;
                  }
                  uintptr_t sample_uidx_base = 0;
                  uintptr_t cur_bits = sex_male[0];
                  for (uint32_t male_idx = 0; male_idx != male_ct; ++male_idx) {
                    const uintptr_t sample_uidx = BitIter1(sex_male, &sample_uidx_base, &cur_bits);
                    cur_covar_vals[sample_uidx] *= 0.5;
                  }
                }
              } else {
                if (unlikely(glm_flags & (kfGlmConditionDominant | kfGlmConditionRecessive))) {
                  logerrputs("Error: --condition[-list] 'dominant'/'recessive' cannot be used with haploid\nvariants.\n");
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                for (uint32_t sample_uidx = 0; sample_uidx != raw_sample_ct; ++sample_uidx) {
                  cur_covar_vals[sample_uidx] *= 0.5;
                }
              }
            }
            strcpy(&(new_covar_names[(local_covar_ct + condition_idx) * new_max_covar_name_blen]), variant_ids[cur_variant_uidx]);
          }
          BigstackEndReset(bigstack_end_mark);
        }
      } else {
        if (unlikely(
                BIGSTACK_ALLOC_X(PhenoCol, raw_covar_ct + add_sex_covar, &new_covar_cols) ||
                bigstack_alloc_c((raw_covar_ct + add_sex_covar) * new_max_covar_name_blen, &new_covar_names))) {
          goto GlmMain_ret_NOMEM;
        }
      }
      memcpy(&(new_covar_cols[condition_ct + local_covar_ct]), covar_cols, orig_covar_ct * sizeof(PhenoCol));
      const char* covar_names_read_iter = covar_names;
      // bugfix (11 May 2017): local covar names come before, not after,
      //   --condition[-list] covar names
      char* covar_names_write_iter = new_covar_names;
      for (uint32_t local_covar_idx = 0; local_covar_idx != local_covar_ct; ++local_covar_idx) {
        memcpy_k(covar_names_write_iter, "LOCAL", 5);
        char* name_end = u32toa(local_covar_idx + 1, &(covar_names_write_iter[5]));
        *name_end = '\0';
        new_covar_cols[local_covar_idx].type_code = kPhenoDtypeOther;
        new_covar_cols[local_covar_idx].nonmiss = nullptr;
        covar_names_write_iter = &(covar_names_write_iter[new_max_covar_name_blen]);
      }
      covar_names_write_iter = &(covar_names_write_iter[condition_ct * new_max_covar_name_blen]);
      for (uint32_t old_covar_idx = 0; old_covar_idx != orig_covar_ct; ++old_covar_idx) {
        strcpy(covar_names_write_iter, covar_names_read_iter);
        covar_names_read_iter = &(covar_names_read_iter[max_covar_name_blen]);
        covar_names_write_iter = &(covar_names_write_iter[new_max_covar_name_blen]);
      }
      if (add_sex_covar) {
        PhenoCol* new_sex_col = &(new_covar_cols[raw_covar_ct++]);
        double* sex_covar_vals;
        if (unlikely(bigstack_alloc_d(raw_sample_ct, &sex_covar_vals))) {
          goto GlmMain_ret_NOMEM;
        }
        uintptr_t sample_uidx_base = 0;
        uintptr_t cur_bits = sex_nm[0];
        for (uint32_t sample_idx = 0; sample_idx != orig_sample_ct; ++sample_idx) {
          const uintptr_t sample_uidx = BitIter1(sex_nm, &sample_uidx_base, &cur_bits);
          // 1/2 instead of 1/0 coding; user shouldn't have to worry about
          // signs changing when they use --sex instead of using the sex column
          // from a .bim/.psam file
          sex_covar_vals[sample_uidx] = kSexMaleToCovarD[IsSet(sex_male, sample_uidx)];
        }
        new_sex_col->category_names = nullptr;
        new_sex_col->nonmiss = K_CAST(uintptr_t*, sex_nm);
        new_sex_col->data.qt = sex_covar_vals;
        new_sex_col->type_code = kPhenoDtypeQt;
        new_sex_col->nonnull_category_ct = 0;
        strcpy_k(covar_names_write_iter, "SEX");
      }
      covar_cols = new_covar_cols;
      covar_names = new_covar_names;
      max_covar_name_blen = new_max_covar_name_blen;
    }
    const uint32_t raw_covar_ctl = BitCtToWordCt(raw_covar_ct);
    uintptr_t* initial_covar_include = nullptr;
    uintptr_t* covar_include = nullptr;
    uintptr_t* cur_sample_include_x_buf = nullptr;
    uintptr_t* covar_include_x = nullptr;
    uint32_t covar_max_nonnull_cat_ct = 0;
    if (raw_covar_ctl) {
      if (unlikely(
              bigstack_alloc_w(raw_covar_ctl, &initial_covar_include) ||
              bigstack_alloc_w(raw_covar_ctl, &covar_include))) {
        goto GlmMain_ret_NOMEM;
      }
      ZeroWArr(raw_covar_ctl, initial_covar_include);
      for (uint32_t covar_uidx = 0; covar_uidx != raw_covar_ct; ++covar_uidx) {
        const PhenoCol* cur_covar_col = &(covar_cols[covar_uidx]);
        if (cur_covar_col->type_code != kPhenoDtypeOther) {
          if (!IsConstCovar(cur_covar_col, orig_sample_include, orig_sample_ct)) {
            SetBit(covar_uidx, initial_covar_include);
            if (cur_covar_col->type_code == kPhenoDtypeCat) {
              if (cur_covar_col->nonnull_category_ct > covar_max_nonnull_cat_ct) {
                covar_max_nonnull_cat_ct = cur_covar_col->nonnull_category_ct;
              }
            }
          } else {
            logerrprintf("Warning: Excluding constant covariate '%s' from --glm.\n", &(covar_names[covar_uidx * max_covar_name_blen]));
          }
        } else {
          // local covariate, always include
          SetBit(covar_uidx, initial_covar_include);
        }
      }
      if (unlikely(covar_max_nonnull_cat_ct && (glm_info_ptr->parameters_range_list.name_ct || glm_info_ptr->tests_range_list.name_ct))) {
        // todo: permit this, and automatically expand a single parameter index
        // referring to a categorical covariate into the appropriate range of
        // final predictor indices
        logerrputs("Error: --parameters/--tests cannot currently be used directly with categorical\ncovariates; expand them into binary covariates with --split-cat-pheno first.\n");
        goto GlmMain_ret_INCONSISTENT_INPUT;
      }
    }
    const uint32_t domdev_present_p1 = domdev_present + 1;
    const uint32_t add_interactions = (glm_flags / kfGlmInteraction) & 1;
    // multiallelic case: one more predictor per extra allele
    const uint32_t biallelic_raw_predictor_ct = 2 + domdev_present + raw_covar_ct * (1 + add_interactions * domdev_present_p1);

    const uint32_t max_extra_allele_ct = MaxAlleleCtSubset(early_variant_include, allele_idx_offsets, raw_variant_ct, variant_ct, PgrGetMaxAlleleCt(simple_pgrp)) - 2;
    if (unlikely(biallelic_raw_predictor_ct + max_extra_allele_ct > 46340)) {
      logerrputs("Error: Too many predictors for --glm.\n");
      if ((biallelic_raw_predictor_ct > 46000) && (max_extra_allele_ct < 23170)) {
        logerrputs("Try reducing the number of covariates.\n");
      } else if ((max_extra_allele_ct > 46000) && (biallelic_raw_predictor_ct < 23170)) {
        logerrputs("Try removing/splitting the variants with the most alternate alleles.\n");
      } else {
        logerrputs("Try reducing the number of covariates, and removing/splitting the variants with\nthe most alternate alleles.\n");
      }
      reterr = kPglRetNotYetSupported;
      goto GlmMain_ret_1;
    }
    const uint32_t biallelic_raw_predictor_ctl = BitCtToWordCt(biallelic_raw_predictor_ct);
    const uint32_t first_covar_pred_uidx = 2 + domdev_present;
    uint32_t first_interaction_pred_uidx = 0;
    if (add_interactions) {
      first_interaction_pred_uidx = first_covar_pred_uidx + raw_covar_ct;
    }

    uintptr_t* raw_parameter_subset = nullptr;
    common.parameter_subset = nullptr;
    common.parameter_subset_x = nullptr;
    common.parameter_subset_y = nullptr;
    common.vif_thresh = vif_thresh;
    common.max_corr = glm_info_ptr->max_corr;
    // bugfix (20 Feb 2018): common.is_xchr_model_1 initialization was either
    // accidentally deleted, or I forgot to add it in the first place...
    common.is_xchr_model_1 = (xchr_model == 1);
    common.tests_flag = glm_info_ptr->tests_range_list.name_ct || (glm_flags & kfGlmTestsAll);
    const uint32_t joint_test = domdev_present || common.tests_flag;
    if (glm_info_ptr->parameters_range_list.name_ct) {
      if (unlikely(
              bigstack_calloc_w(biallelic_raw_predictor_ctl, &raw_parameter_subset) ||
              bigstack_alloc_w(biallelic_raw_predictor_ctl, &common.parameter_subset) ||
              bigstack_alloc_w(biallelic_raw_predictor_ctl, &common.parameter_subset_x) ||
              bigstack_alloc_w(biallelic_raw_predictor_ctl, &common.parameter_subset_y))) {
        goto GlmMain_ret_NOMEM;
      }
      raw_parameter_subset[0] = 1;  // intercept (index 0) always included
      NumericRangeListToBitarr(&(glm_info_ptr->parameters_range_list), biallelic_raw_predictor_ct, 0, 1, raw_parameter_subset);
      if (unlikely(domdev_present && ((raw_parameter_subset[0] & 7) != 7))) {
        // this breaks the joint test
        logerrputs("Error: --parameters cannot exclude 1 or 2 when the 'genotypic' or 'hethom'\nmodifier is present.\n");
        goto GlmMain_ret_INVALID_CMDLINE;
      }
      if (unlikely((glm_flags & kfGlmHideCovar) && (!joint_test) && (!(raw_parameter_subset[0] & 2)))) {
        logerrputs("Error: 'hide-covar' modifier suppresses all output due to --parameters setting.\n");
        goto GlmMain_ret_INVALID_CMDLINE;
      }
    }
    uintptr_t* raw_joint_test_params = nullptr;
    uintptr_t* joint_test_params_buf = nullptr;
    common.joint_test_params = nullptr;
    common.joint_test_params_x = nullptr;
    common.joint_test_params_y = nullptr;
    common.constraint_ct = 0;
    common.constraint_ct_x = 0;
    common.constraint_ct_y = 0;
    if (joint_test) {
      if (unlikely(
              bigstack_calloc_w(biallelic_raw_predictor_ctl, &raw_joint_test_params) ||
              bigstack_alloc_w(biallelic_raw_predictor_ctl, &common.joint_test_params) ||
              bigstack_alloc_w(biallelic_raw_predictor_ctl, &common.joint_test_params_x))) {
        goto GlmMain_ret_NOMEM;
      }

      // includes intercept
      uint32_t raw_param_ct = biallelic_raw_predictor_ct;
      if (raw_parameter_subset) {
        if (unlikely(bigstack_alloc_w(biallelic_raw_predictor_ctl, &joint_test_params_buf))) {
          goto GlmMain_ret_NOMEM;
        }
        raw_param_ct = PopcountWords(raw_parameter_subset, biallelic_raw_predictor_ctl);
      }

      if (!common.tests_flag) {
        // 1, 2
        raw_joint_test_params[0] = 6;
      } else {
        uintptr_t* tests_buf;
        if (unlikely(
                bigstack_alloc_w(biallelic_raw_predictor_ctl, &common.joint_test_params_y) ||
                bigstack_calloc_w(biallelic_raw_predictor_ctl, &tests_buf))) {
          goto GlmMain_ret_NOMEM;
        }
        if (glm_info_ptr->tests_range_list.name_ct) {
          if (unlikely(NumericRangeListToBitarr(&(glm_info_ptr->tests_range_list), raw_param_ct, 0, 0, tests_buf))) {
            logerrputs("Error: Invalid --tests expression.\n");
            goto GlmMain_ret_INVALID_CMDLINE;
          }
        } else {
          FillBitsNz(1, raw_param_ct, tests_buf);
        }
        if (raw_parameter_subset) {
          ExpandBytearr(tests_buf, raw_parameter_subset, biallelic_raw_predictor_ctl, raw_param_ct, 0, raw_joint_test_params);
        } else {
          memcpy(raw_joint_test_params, tests_buf, biallelic_raw_predictor_ctl * sizeof(intptr_t));
        }
        BigstackReset(tests_buf);
      }
    }
    if (raw_parameter_subset) {
      if (add_sex_covar && first_interaction_pred_uidx) {
        // special case: when add_sex_covar is true, the added sex covariate is
        // simply the last covariate, with predictor index
        // (first_interaction_pred_uidx - 1).  This lines up with --parameters
        // when interactions are not requested; but when they are, we have a
        // small reshuffle to do.
        uintptr_t* parameter_subset_reshuffle_buf;
        if (unlikely(bigstack_alloc_w(biallelic_raw_predictor_ctl, &parameter_subset_reshuffle_buf))) {
          goto GlmMain_ret_NOMEM;
        }
        SexInteractionReshuffle(first_interaction_pred_uidx, raw_covar_ct, domdev_present, biallelic_raw_predictor_ctl, raw_parameter_subset, parameter_subset_reshuffle_buf);
        if (common.tests_flag) {
          SexInteractionReshuffle(first_interaction_pred_uidx, raw_covar_ct, domdev_present, biallelic_raw_predictor_ctl, raw_joint_test_params, parameter_subset_reshuffle_buf);
        }
        BigstackReset(parameter_subset_reshuffle_buf);
      }
      // if there were any constant covariates, exclude them from
      // raw_parameter_subset... but only after --tests has been processed, so
      // we handle indexes properly.
      // note that, if appended sex covariate is present at all, it is always
      // nonconstant.
      if (initial_covar_include) {
        const uint32_t nonconst_covar_ct = PopcountWords(initial_covar_include, raw_covar_ctl);
        const uint32_t removed_covar_ct = raw_covar_ct - nonconst_covar_ct;
        uintptr_t covar_uidx_base = 0;
        uintptr_t cur_inv_bits = ~initial_covar_include[0];
        for (uint32_t removed_covar_idx = 0; removed_covar_idx != removed_covar_ct; ++removed_covar_idx) {
          const uintptr_t covar_uidx = BitIter0(initial_covar_include, &covar_uidx_base, &cur_inv_bits);
          ClearBit(first_covar_pred_uidx + covar_uidx, raw_parameter_subset);
          if (first_interaction_pred_uidx) {
            const uint32_t geno_interaction_uidx = first_interaction_pred_uidx + covar_uidx * domdev_present_p1;
            ClearBit(geno_interaction_uidx, raw_parameter_subset);
            if (domdev_present) {
              ClearBit(geno_interaction_uidx + 1, raw_parameter_subset);
            }
          }
        }
        // if any loaded nonconstant covariates aren't referenced in
        // raw_parameter_subset, remove them from initial_covar_include
        covar_uidx_base = 0;
        uintptr_t cur_bits = initial_covar_include[0];
        for (uint32_t nonconst_covar_idx = 0; nonconst_covar_idx != nonconst_covar_ct; ++nonconst_covar_idx) {
          const uintptr_t covar_uidx = BitIter1(initial_covar_include, &covar_uidx_base, &cur_bits);
          uint32_t cur_covar_is_referenced = IsSet(raw_parameter_subset, first_covar_pred_uidx + covar_uidx);
          if (add_interactions) {
            cur_covar_is_referenced = cur_covar_is_referenced || IsSet(raw_parameter_subset, first_interaction_pred_uidx + covar_uidx * domdev_present_p1);
            if (domdev_present) {
              cur_covar_is_referenced = cur_covar_is_referenced || IsSet(raw_parameter_subset, first_interaction_pred_uidx + covar_uidx * 2 + 1);
            }
          }
          if (!cur_covar_is_referenced) {
            ClearBit(covar_uidx, initial_covar_include);
          }
        }
        // May as well remove constant covariates from joint test.
        if (common.tests_flag) {
          BitvecAnd(raw_parameter_subset, biallelic_raw_predictor_ctl, raw_joint_test_params);
        }
      }
      // if your regression doesn't involve genotype data, you should be using
      // e.g. R, not plink...
      if (unlikely(
              (!(raw_parameter_subset[0] & 2)) &&
              ((!domdev_present) || (!(raw_parameter_subset[0] & 4))) &&
              ((!add_interactions) || (!PopcountBitRange(raw_parameter_subset, first_interaction_pred_uidx, biallelic_raw_predictor_ct))))) {
        logerrputs("Error: --parameters must retain at least one dosage-dependent variable.\n");
        goto GlmMain_ret_INCONSISTENT_INPUT;
      }
    }
    // computation of these counts moved here, since --parameters can reduce
    // the number of relevant covariates
    uint32_t initial_nonx_covar_ct = 0;
    if (initial_covar_include) {
      initial_nonx_covar_ct = PopcountWords(initial_covar_include, raw_covar_ctl);
    }
    uint32_t initial_y_covar_ct = 0;
    uintptr_t* covar_include_y = nullptr;
    if (!initial_nonx_covar_ct) {
      // BigstackReset(initial_covar_include);  // not ok with parameters
      initial_covar_include = nullptr;
      covar_include = nullptr;
    } else {
      initial_y_covar_ct = initial_nonx_covar_ct - (cur_sample_include_y_buf && add_sex_covar && IsSet(initial_covar_include, raw_covar_ct - 1));
      if (add_sex_covar && (!(glm_flags & kfGlmSex))) {
        // may as well verify there's at least one non-x/non-y variant
        // (if only chrX and chrY present, don't allocate
        // cur_sample_include_x_buf, just make chrX the baseline instead)
        if (IsSet(initial_covar_include, raw_covar_ct - 1) && (variant_ct != variant_ct_x + variant_ct_y)) {
          if (unlikely(
                  bigstack_alloc_w(raw_sample_ctl, &cur_sample_include_x_buf) ||
                  bigstack_alloc_w(raw_covar_ctl, &covar_include_x))) {
            goto GlmMain_ret_NOMEM;
          }
          --initial_nonx_covar_ct;
        }
      }
      if (cur_sample_include_y_buf) {
        if (unlikely(bigstack_alloc_w(raw_covar_ctl, &covar_include_y))) {
          goto GlmMain_ret_NOMEM;
        }
      }
    }

    const uint32_t report_adjust = (adjust_info_ptr->flags & kfAdjustColAll);
    const uint32_t is_sometimes_firth = !(glm_flags & kfGlmNoFirth);
    const uint32_t is_always_firth = glm_flags & kfGlmFirth;
    const uint32_t skip_modifier = (glm_flags / kfGlmSkip) & 1;
    const uint32_t glm_pos_col = glm_info_ptr->cols & kfGlmColPos;
    const uint32_t gcount_cc_col = glm_info_ptr->cols & kfGlmColGcountcc;
    const uint32_t xtx_state = (add_interactions || local_covar_ct)? 0 : domdev_present_p1;

    const uintptr_t raw_allele_ct = allele_idx_offsets? allele_idx_offsets[raw_variant_ct] : (2 * raw_variant_ct);
    const uintptr_t raw_allele_ctl = BitCtToWordCt(raw_allele_ct);

    common.max_extra_allele_ct = max_extra_allele_ct;

    const uint32_t pheno_ctl = BitCtToWordCt(pheno_ct);
    uintptr_t* pheno_include;
    if (unlikely(bigstack_alloc_w(pheno_ctl, &pheno_include))) {
      goto GlmMain_ret_NOMEM;
    }
    SetAllBits(pheno_ct, pheno_include);

    uintptr_t* valid_variants = nullptr;
    uintptr_t* valid_alleles = nullptr;
    unsigned char* bigstack_mark2 = g_bigstack_base;
    if (report_adjust || perms_total) {
      if (unlikely(
              bigstack_alloc_w(raw_variant_ctl, &valid_variants) ||
              bigstack_alloc_w(raw_allele_ctl, &valid_alleles))) {
        goto GlmMain_ret_NOMEM;
      }
      bigstack_mark2 = g_bigstack_base;
    } else if (pheno_ct > 1) {
      // When there are multiple quantitative phenotypes with the same
      // missingness pattern, they can be processed more efficiently together.
      uintptr_t* pheno_batch;
      uint32_t* pheno_nm_hashes;
      uintptr_t* pheno_nonmiss_tmp; // might be able to move this later
      if (unlikely(
              bigstack_alloc_w(pheno_ctl, &pheno_batch) ||
              bigstack_alloc_u32(pheno_ct, &pheno_nm_hashes) ||
              bigstack_alloc_w(raw_sample_ctl, &pheno_nonmiss_tmp))) {
        goto GlmMain_ret_NOMEM;
      }
      uint32_t completed_pheno_ct = 0;
      for (uint32_t pheno_uidx = 0; pheno_uidx != pheno_ct; ++pheno_uidx) {
        const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_uidx]);
        const PhenoDtype dtype_code = cur_pheno_col->type_code;
        if (dtype_code != kPhenoDtypeQt) {
          ClearBit(pheno_uidx, pheno_include);
          continue;
        }
        BitvecAndCopy(orig_sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include);
        const uint32_t sample_ct = PopcountWords(cur_sample_include, raw_sample_ctl);
        if (IsConstCovar(cur_pheno_col, cur_sample_include, sample_ct)) {
          const char* cur_pheno_name = &(pheno_names[pheno_uidx * max_pheno_name_blen]);
          if (unlikely(!skip_modifier)) {
            logerrprintfww("Error: --glm quantitative phenotype '%s' is constant.\n", cur_pheno_name);
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logprintfww("--glm: Skipping constant quantitative phenotype '%s'.\n", cur_pheno_name);
          ClearBit(pheno_uidx, pheno_include);
          continue;
        }
        // May want to use 64-bit XXH3 hashes instead.  But let's get this
        // working first.
        pheno_nm_hashes[pheno_uidx] = Hash32(cur_sample_include, raw_sample_ctl * sizeof(intptr_t));
      }
      unsigned char* bigstack_mark3 = g_bigstack_base;

      for (uint32_t pheno_uidx = 0; pheno_uidx != pheno_ct; ++pheno_uidx) {
        if (!IsSet(pheno_include, pheno_uidx)) {
          continue;
        }
        ZeroWArr(pheno_ctl, pheno_batch);
        SetBit(pheno_uidx, pheno_batch);
        const uint32_t cur_hash = pheno_nm_hashes[pheno_uidx];
        const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_uidx]);
        BitvecAndCopy(orig_sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include);
        // possible todo: switch to O(n log n)-expected-time algorithm here.
        // But unlikely to matter since this cost should be minuscule compared
        // to the actual regressions, and in the
        // all-phenotype-values-are-present case this is O(n) anyway.
        for (uint32_t pheno_uidx2 = pheno_uidx + 1; pheno_uidx2 != pheno_ct; ++pheno_uidx2) {
          if ((!IsSet(pheno_include, pheno_uidx2)) || (cur_hash != pheno_nm_hashes[pheno_uidx2])) {
            continue;
          }
          const PhenoCol* candidate_pheno_col = &(pheno_cols[pheno_uidx2]);
          BitvecAndCopy(orig_sample_include, candidate_pheno_col->nonmiss, raw_sample_ctl, pheno_nonmiss_tmp);
          if (memequal(cur_sample_include, pheno_nonmiss_tmp, raw_sample_ctl * sizeof(intptr_t))) {
            SetBit(pheno_uidx2, pheno_batch);
          }
        }
        uint32_t batch_size = PopcountWords(pheno_batch, pheno_ctl);
        if (batch_size == 1) {
          continue;
        }
        // probable todo: pull out lots of shared code with the non-batch
        // phenotype loop.
        uint32_t sample_ct = PopcountWords(cur_sample_include, raw_sample_ctl);
        uint32_t covar_ct = 0;
        uint32_t extra_cat_ct = 0;
        BigstackDoubleReset(bigstack_mark3, bigstack_end_mark);
        if (initial_nonx_covar_ct) {
          uint16_t dummy = 0;
          if (unlikely(GlmDetermineCovars(nullptr, initial_covar_include, covar_cols, raw_sample_ct, raw_covar_ctl, initial_nonx_covar_ct, covar_max_nonnull_cat_ct, 0, 0, cur_sample_include, covar_include, &sample_ct, &covar_ct, &extra_cat_ct, &dummy))) {
            goto GlmMain_ret_NOMEM;
          }
        }
        uint32_t biallelic_predictor_ct = 2 + domdev_present + (covar_ct + extra_cat_ct) * (1 + add_interactions * domdev_present_p1);
        if (raw_parameter_subset) {
          biallelic_predictor_ct = CollapseParamOrTestSubset(covar_include, raw_parameter_subset, domdev_present, raw_covar_ct, covar_ct, add_interactions, common.parameter_subset);
        }
        if (raw_joint_test_params) {
          common.constraint_ct = CollapseParamOrTestSubset(covar_include, raw_joint_test_params, domdev_present, raw_covar_ct, covar_ct, add_interactions, common.joint_test_params);
          if (raw_parameter_subset) {
            memcpy(joint_test_params_buf, common.joint_test_params, biallelic_raw_predictor_ctl * sizeof(intptr_t));
            ZeroWArr(biallelic_raw_predictor_ctl, common.joint_test_params);
            CopyBitarrSubset(joint_test_params_buf, common.parameter_subset, biallelic_predictor_ct, common.joint_test_params);
          }
        }
        const char* first_pheno_name = &(pheno_names[pheno_uidx * max_pheno_name_blen]);
        if (sample_ct <= biallelic_predictor_ct) {
          if (unlikely(!skip_modifier)) {
            logerrprintfww("Error: # samples <= # predictor columns for --glm phenotype '%s'.\n", first_pheno_name);
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logprintfww("Note: Skipping --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since # samples <= # predictor columns.\n", first_pheno_name);
          BitvecInvmask(pheno_batch, pheno_ctl, pheno_include);
          continue;
        }
#ifdef __LP64__
        if (RoundUpPow2(sample_ct, 4) * S_CAST(uint64_t, biallelic_predictor_ct + max_extra_allele_ct) > 0x7fffffff) {
          // todo: remove this constraint in LAPACK_ILP64 case?
          if (unlikely(!skip_modifier)) {
            logerrprintfww("Error: Too many samples or predictors for --glm regression on phenotype '%s' (internal matrices currently limited to ~2^31 entries).\n", first_pheno_name);
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logerrprintfww("Warning: Skipping --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since there are too many samples or predictors (internal matrices currently limited to ~2^31 entries).\n", first_pheno_name);
          BitvecInvmask(pheno_batch, pheno_ctl, pheno_include);
          continue;
        }
#endif
        if (common.tests_flag && (!common.constraint_ct)) {
          if (unlikely(!skip_modifier)) {
            logerrprintfww("Error: --tests predictor(s) are constant for all remaining samples for --glm phenotype '%s'.\n", first_pheno_name);
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logprintfww("Note: Skipping --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since --tests predictor(s) are constant for all remaining samples.\n", first_pheno_name);
          BitvecInvmask(pheno_batch, pheno_ctl, pheno_include);
          continue;
        }
        if (covar_ct < initial_nonx_covar_ct) {
          uintptr_t covar_uidx_base = 0;
          uintptr_t cur_bits = initial_covar_include[0];
          for (uint32_t covar_idx = 0; covar_idx != initial_nonx_covar_ct; ++covar_idx) {
            const uintptr_t covar_uidx = BitIter1(initial_covar_include, &covar_uidx_base, &cur_bits);
            if (!IsSet(covar_include, covar_uidx)) {
              logerrprintfww("Warning: %sot including covariate '%s' in --glm regression on phenotype '%s', and other(s) with identical missingness patterns.\n", cur_sample_include_x_buf? (cur_sample_include_y_buf? "Outside of chrX, n" : "Outside of chrX and chrY, n") : (cur_sample_include_y_buf? "Outside of chrY, n" : "N"), &(covar_names[covar_uidx * max_covar_name_blen]), first_pheno_name);
            }
          }
        }

        // cur_sample_include_x == nullptr: chrX uses same samples and
        //   covariates as the rest of the genome.  sample_ct_x always zero to
        //   force most chrX-specific initialization to be skipped (exception:
        //   sex_male_collapsed, needed for allele count/freq reporting)
        // cur_sample_include_x non-null: if sample_ct_x == 0, we skip the
        //   entire chromosome.  otherwise, we have different covariates than
        //   the rest of the genome.
        uintptr_t* cur_sample_include_x = cur_sample_include_x_buf;
        uint32_t sample_ct_x = 0;
        uint32_t covar_ct_x = 0;
        uint32_t extra_cat_ct_x = 0;
        uint32_t biallelic_predictor_ct_x = 0;
        uint32_t x_samples_are_different = 0;
        if (cur_sample_include_x) {
          BitvecAndCopy(orig_sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include_x);
          uint16_t dummy = 0;
          if (unlikely(GlmDetermineCovars(nullptr, initial_covar_include, covar_cols, raw_sample_ct, raw_covar_ctl, initial_nonx_covar_ct + 1, covar_max_nonnull_cat_ct, 0, 0, cur_sample_include_x, covar_include_x, &sample_ct_x, &covar_ct_x, &extra_cat_ct_x, &dummy))) {
            goto GlmMain_ret_NOMEM;
          }
          x_samples_are_different = (sample_ct_x != sample_ct) || (!wordsequal(cur_sample_include, cur_sample_include_x, raw_sample_ctl));
          if ((!x_samples_are_different) && (covar_ct == covar_ct_x) && wordsequal(covar_include, covar_include_x, raw_covar_ctl)) {
            logprintfww("Note: chrX samples and covariate(s) in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, are the same as that for the rest of the genome.\n", first_pheno_name);
            sample_ct_x = 0;
            cur_sample_include_x = nullptr;
          } else {
            if (!sample_ct_x) {
              if (unlikely(!skip_modifier)) {
                logerrprintfww("Error: --glm regression on phenotype '%s' is degenerate on chrX.\n", first_pheno_name);
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', and other(s) with identical missingness patterns.\n", first_pheno_name);
            } else {
              biallelic_predictor_ct_x = 2 + domdev_present + (covar_ct_x + extra_cat_ct_x) * (1 + add_interactions * domdev_present_p1);
              if (raw_parameter_subset) {
                biallelic_predictor_ct_x = CollapseParamOrTestSubset(covar_include, raw_parameter_subset, domdev_present, raw_covar_ct, covar_ct_x, add_interactions, common.parameter_subset_x);
              }
              if (raw_joint_test_params) {
                common.constraint_ct_x = CollapseParamOrTestSubset(covar_include, raw_joint_test_params, domdev_present, raw_covar_ct, covar_ct_x, add_interactions, common.joint_test_params_x);
                if (raw_parameter_subset) {
                  memcpy(joint_test_params_buf, common.joint_test_params_x, biallelic_raw_predictor_ctl * sizeof(intptr_t));
                  ZeroWArr(biallelic_raw_predictor_ctl, common.joint_test_params_x);
                  CopyBitarrSubset(joint_test_params_buf, common.parameter_subset, biallelic_predictor_ct_x, common.joint_test_params_x);
                }
              }
              if (sample_ct_x <= biallelic_predictor_ct_x) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: # samples <= # predictor columns for --glm phenotype '%s' on chrX.\n", first_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since # remaining samples <= # predictor columns.\n", first_pheno_name);
                sample_ct_x = 0;
#ifdef __LP64__
              } else if (RoundUpPow2(sample_ct_x, 4) * S_CAST(uint64_t, biallelic_predictor_ct_x + max_extra_allele_ct) > 0x7fffffff) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: Too many samples or predictors for --glm regression on phenotype '%s' on chrX (internal matrices currently limited to ~2^31 entries).\n", first_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logerrprintfww("Warning: Skipping chrX in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since there are too many samples or predictors (internal matrices currently limited to ~2^31 entries).\n", first_pheno_name);
                sample_ct_x = 0;
#endif
              } else {
                for (uint32_t pheno_uidx2 = pheno_uidx; pheno_uidx2 != pheno_ct; ++pheno_uidx2) {
                  if (!IsSet(pheno_batch, pheno_uidx2)) {
                    continue;
                  }
                  const PhenoCol* candidate_pheno_col = &(pheno_cols[pheno_uidx2]);
                  if (IsConstCovar(candidate_pheno_col, cur_sample_include_x, sample_ct_x)) {
                    // Punt to single-phenotype-at-a-time handler.
                    ClearBit(pheno_uidx2, pheno_batch);
                  }
                }
                // possible todo: reset first_pheno_name here
              }
              if (sample_ct_x && common.tests_flag && (!common.constraint_ct_x)) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: --tests predictor(s) are constant for all remaining samples on chrX for --glm phenotype '%s'.\n", first_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since --tests predictor(s) are constant for all remaining samples.\n", first_pheno_name);
                sample_ct_x = 0;
              }
              if (sample_ct_x && (covar_ct_x < initial_nonx_covar_ct + 1)) {
                uintptr_t covar_uidx_base = 0;
                uintptr_t cur_bits = initial_covar_include[0];
                for (uint32_t covar_idx = 0; covar_idx != covar_ct_x; ++covar_idx) {
                  const uintptr_t covar_uidx = BitIter1(initial_covar_include, &covar_uidx_base, &cur_bits);
                  if (!IsSet(covar_include_x, covar_uidx)) {
                    logerrprintfww("Warning: On chrX, not including covariate '%s' in --glm regression on phenotype '%s', and other(s) with identical missingness patterns.\n", &(covar_names[covar_uidx * max_covar_name_blen]), first_pheno_name);
                  }
                }
              }
            }
          }
        }

        uintptr_t* cur_sample_include_y = cur_sample_include_y_buf;
        uint32_t sample_ct_y = 0;
        uint32_t covar_ct_y = 0;
        uint32_t extra_cat_ct_y = 0;
        uint32_t biallelic_predictor_ct_y = 0;
        uint32_t y_samples_are_different = 0;
        if (cur_sample_include_y) {
          BitvecAndCopy(orig_sample_include, sex_male, raw_sample_ctl, cur_sample_include_y);
          BitvecAnd(cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include_y);
          uint16_t dummy = 0;
          if (unlikely(GlmDetermineCovars(nullptr, initial_covar_include, covar_cols, raw_sample_ct, raw_covar_ctl, initial_y_covar_ct, covar_max_nonnull_cat_ct, 0, 0, cur_sample_include_y, covar_include_y, &sample_ct_y, &covar_ct_y, &extra_cat_ct_y, &dummy))) {
            goto GlmMain_ret_NOMEM;
          }
          y_samples_are_different = (sample_ct_y != sample_ct) || (!wordsequal(cur_sample_include, cur_sample_include_y, raw_sample_ctl));
          if ((!y_samples_are_different) && (covar_ct == covar_ct_y) && wordsequal(covar_include, covar_include_y, raw_covar_ctl)) {
            logprintfww("Note: chrY samples and covariate(s) in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, are the same as that for the rest of the genome.\n", first_pheno_name);
            sample_ct_y = 0;
            cur_sample_include_y = nullptr;
          } else {
            if (!sample_ct_y) {
              if (unlikely(!skip_modifier)) {
                logerrprintfww("Error: --glm regression on phenotype '%s' is degenerate on chrY.\n", first_pheno_name);
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', and other(s) with identical missingness patterns.\n", first_pheno_name);
            } else {
              biallelic_predictor_ct_y = 2 + domdev_present + (covar_ct_y + extra_cat_ct_y) * (1 + add_interactions * domdev_present_p1);
              if (raw_parameter_subset) {
                biallelic_predictor_ct_y = CollapseParamOrTestSubset(covar_include, raw_parameter_subset, domdev_present, raw_covar_ct, covar_ct_y, add_interactions, common.parameter_subset_y);
              }
              if (raw_joint_test_params) {
                assert(common.tests_flag);
                common.constraint_ct_y = CollapseParamOrTestSubset(covar_include, raw_joint_test_params, domdev_present, raw_covar_ct, covar_ct_y, add_interactions, common.joint_test_params_y);
                if (raw_parameter_subset) {
                  memcpy(joint_test_params_buf, common.joint_test_params_y, biallelic_raw_predictor_ctl * sizeof(intptr_t));
                  ZeroWArr(biallelic_raw_predictor_ctl, common.joint_test_params_y);
                  CopyBitarrSubset(joint_test_params_buf, common.parameter_subset, biallelic_predictor_ct_y, common.joint_test_params_y);
                }
              }
              if (sample_ct_y <= biallelic_predictor_ct_y) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: # samples <= # predictor columns for --glm phenotype '%s' on chrY.\n", first_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since # remaining samples <= # predictor columns.\n", first_pheno_name);
                sample_ct_y = 0;
#ifdef __LP64__
              } else if (RoundUpPow2(sample_ct_y, 4) * S_CAST(uint64_t, biallelic_predictor_ct_y + max_extra_allele_ct) > 0x7fffffff) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: Too many samples or predictors for --glm regression on phenotype '%s' on chrY (internal matrices currently limited to ~2^31 entries).\n", first_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logerrprintfww("Warning: Skipping chrY in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since there are too many samples or predictors (internal matrices currently limited to ~2^31 entries).\n", first_pheno_name);
                sample_ct_y = 0;
#endif
              } else {
                for (uint32_t pheno_uidx2 = pheno_uidx; pheno_uidx2 != pheno_ct; ++pheno_uidx2) {
                  if (!IsSet(pheno_batch, pheno_uidx2)) {
                    continue;
                  }
                  const PhenoCol* candidate_pheno_col = &(pheno_cols[pheno_uidx2]);
                  if (IsConstCovar(candidate_pheno_col, cur_sample_include_y, sample_ct_y)) {
                    // Punt to single-phenotype-at-a-time handler.
                    ClearBit(pheno_uidx2, pheno_batch);
                  }
                }
              }
              if (sample_ct_y && common.tests_flag && (!common.constraint_ct_y)) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: --tests predictor(s) are constant for all remaining samples on chrY for --glm phenotype '%s'.\n", first_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since --tests predictor(s) are constant for all remaining samples.\n", first_pheno_name);
                sample_ct_y = 0;
              }
              if (sample_ct_y && (covar_ct_y < initial_y_covar_ct)) {
                uintptr_t covar_uidx_base = 0;
                uintptr_t cur_bits = initial_covar_include[0];
                for (uint32_t covar_idx = 0; covar_idx != covar_ct_y; ++covar_idx) {
                  const uintptr_t covar_uidx = BitIter1(initial_covar_include, &covar_uidx_base, &cur_bits);
                  if (!IsSet(covar_include_y, covar_uidx)) {
                    logerrprintfww("Warning: On chrY, not including covariate '%s' in --glm regression on phenotype '%s', and other(s) with identical missingness patterns.\n", &(covar_names[covar_uidx * max_covar_name_blen]), first_pheno_name);
                  }
                }
              }
            }
          }
        }

        batch_size = PopcountWords(pheno_batch, pheno_ctl);
        if (batch_size <= 1) {
          continue;
        }

        // Expand categorical covariates and perform VIF and correlation checks
        // here.
        double* covars_cmaj_d = nullptr;
        const char** cur_covar_names = nullptr;
        GlmErr glm_err;
        if (unlikely(GlmAllocFillAndTestCovarsQt(cur_sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &common.nm_precomp, &covars_cmaj_d, &cur_covar_names, &glm_err))) {
          goto GlmMain_ret_NOMEM;
        }
        if (glm_err) {
          PrintPrescanErrmsg("", first_pheno_name, cur_covar_names, glm_err, local_covar_ct, skip_modifier, 1);
          if (unlikely(!skip_modifier)) {
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          BitvecInvmask(pheno_batch, pheno_ctl, pheno_include);
          continue;
        }
        const char** cur_covar_names_x = nullptr;
        common.nm_precomp_x = nullptr;
        linear_ctx.covars_cmaj_x_d = nullptr;
        if (sample_ct_x) {
          if (unlikely(GlmAllocFillAndTestCovarsQt(cur_sample_include_x, covar_include_x, covar_cols, covar_names, sample_ct_x, covar_ct_x, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_x, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &common.nm_precomp_x, &linear_ctx.covars_cmaj_x_d, &cur_covar_names_x, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
          if (glm_err) {
            PrintPrescanErrmsg("chrX in ", first_pheno_name, cur_covar_names_x, glm_err, local_covar_ct, skip_modifier, 1);
            if (unlikely(!skip_modifier)) {
              goto GlmMain_ret_INCONSISTENT_INPUT;
            }
            sample_ct_x = 0;
          }
        }
        const char** cur_covar_names_y = nullptr;
        common.nm_precomp_y = nullptr;
        linear_ctx.covars_cmaj_y_d = nullptr;
        if (sample_ct_y) {
          if (unlikely(GlmAllocFillAndTestCovarsQt(cur_sample_include_y, covar_include_y, covar_cols, covar_names, sample_ct_y, covar_ct_y, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_y, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &common.nm_precomp_y, &linear_ctx.covars_cmaj_y_d, &cur_covar_names_y, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
          if (glm_err) {
            PrintPrescanErrmsg("chrY in ", first_pheno_name, cur_covar_names_y, glm_err, local_covar_ct, skip_modifier, 1);
            if (unlikely(!skip_modifier)) {
              goto GlmMain_ret_INCONSISTENT_INPUT;
            }
            sample_ct_y = 0;
          }
        }

        const char** cur_test_names = nullptr;
        const char** cur_test_names_x = nullptr;
        const char** cur_test_names_y = nullptr;
        if (unlikely(AllocAndInitReportedTestNames(common.parameter_subset, cur_covar_names, glm_flags, covar_ct + extra_cat_ct, common.tests_flag? common.constraint_ct : 0, &cur_test_names))) {
          goto GlmMain_ret_NOMEM;
        }
        if (sample_ct_x) {
          if (unlikely(AllocAndInitReportedTestNames(common.parameter_subset_x, cur_covar_names_x, glm_flags, covar_ct_x + extra_cat_ct_x, common.tests_flag? common.constraint_ct_x : 0, &cur_test_names_x))) {
            goto GlmMain_ret_NOMEM;
          }
        }
        if (sample_ct_y) {
          if (unlikely(AllocAndInitReportedTestNames(common.parameter_subset_y, cur_covar_names_y, glm_flags, covar_ct_y + extra_cat_ct_y, common.tests_flag? common.constraint_ct_y : 0, &cur_test_names_y))) {
            goto GlmMain_ret_NOMEM;
          }
        }

        const uintptr_t* cur_variant_include = early_variant_include;
        const uintptr_t* cur_local_variant_include = local_variant_include;
        const uint32_t skip_x = variant_ct_x && ((!xchr_model) || (cur_sample_include_x && (!sample_ct_x)));
        const uint32_t skip_y = variant_ct_y && ((!male_ct) || (cur_sample_include_y && (!sample_ct_y)));
        uint32_t cur_variant_ct = variant_ct;
        if (skip_x || skip_y) {
          uintptr_t* tmp_variant_include;
          if (unlikely(bigstack_alloc_w(raw_variant_ctl, &tmp_variant_include))) {
            goto GlmMain_ret_NOMEM;
          }
          memcpy(tmp_variant_include, early_variant_include, raw_variant_ctl * sizeof(intptr_t));
          uintptr_t* tmp_local_variant_include = nullptr;
          if (local_variant_include) {
            if (unlikely(bigstack_alloc_w(local_variant_ctl, &tmp_local_variant_include))) {
              goto GlmMain_ret_NOMEM;
            }
            memcpy(tmp_local_variant_include, local_variant_include, local_variant_ctl * sizeof(intptr_t));
          }
          if (skip_x) {
            if (local_variant_include) {
              const uint32_t variant_ct_before_x = PopcountBitRange(early_variant_include, 0, x_start);
              uint32_t local_uidx_first = IdxToUidxBasic(local_variant_include, variant_ct_before_x);
              uint32_t local_uidx_last = FindNth1BitFrom(local_variant_include, local_uidx_first, variant_ct_x);
              ClearBitsNz(local_uidx_first, local_uidx_last + 1, tmp_local_variant_include);
            }
            ClearBitsNz(x_start, x_end, tmp_variant_include);
            cur_variant_ct -= variant_ct_x;
          }
          if (skip_y) {
            if (local_variant_include) {
              const uint32_t variant_ct_before_y = PopcountBitRange(early_variant_include, 0, y_start);
              uint32_t local_uidx_first = IdxToUidxBasic(local_variant_include, variant_ct_before_y);
              uint32_t local_uidx_last = FindNth1BitFrom(local_variant_include, local_uidx_first, variant_ct_y);
              ClearBitsNz(local_uidx_first, local_uidx_last + 1, tmp_local_variant_include);
            }
            ClearBitsNz(y_start, y_end, tmp_variant_include);
            cur_variant_ct -= variant_ct_y;
          }
          cur_variant_include = tmp_variant_include;
          cur_local_variant_include = tmp_local_variant_include;
        }
        if (sex_male_collapsed_buf && (!skip_x)) {
          if (!cur_sample_include_x) {
            CopyBitarrSubset(sex_male, cur_sample_include, sample_ct, sex_male_collapsed_buf);
          } else {
            CopyBitarrSubset(sex_male, cur_sample_include_x, sample_ct_x, sex_male_collapsed_buf);
          }
        }
        FillCumulativePopcounts(cur_sample_include, raw_sample_ctl, common.sample_include_cumulative_popcounts);
        common.sample_ct = sample_ct;
        common.sample_ct_x = sample_ct_x;
        common.covar_ct = covar_ct + extra_cat_ct;
        common.local_covar_ct = local_covar_ct;
        if (sample_ct_x) {
          if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &common.sample_include_x_cumulative_popcounts))) {
            goto GlmMain_ret_NOMEM;
          }
          FillCumulativePopcounts(cur_sample_include_x, raw_sample_ctl, common.sample_include_x_cumulative_popcounts);
          common.sample_include_x = cur_sample_include_x;
          common.covar_ct_x = covar_ct_x + extra_cat_ct_x;
          // common.male_ct = PopcountWordsIntersect(cur_sample_include_x, sex_male, raw_sample_ctl);
        } else {
          // technically only need this if variant_ct_x && (!skip_x)
          // common.male_ct = PopcountWordsIntersect(cur_sample_include, sex_male, raw_sample_ctl);

          // defensive
          common.sample_include_x = nullptr;
          common.sample_include_x_cumulative_popcounts = nullptr;
          common.covar_ct_x = 0;
        }
        common.sample_ct_y = sample_ct_y;
        if (sample_ct_y) {
          if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &common.sample_include_y_cumulative_popcounts))) {
            goto GlmMain_ret_NOMEM;
          }
          FillCumulativePopcounts(cur_sample_include_y, raw_sample_ctl, common.sample_include_y_cumulative_popcounts);
          common.sample_include_y = cur_sample_include_y;
          common.covar_ct_y = covar_ct_y + extra_cat_ct_y;
        } else {
          common.sample_include_y = nullptr;
          common.sample_include_y_cumulative_popcounts = nullptr;
          common.covar_ct_y = 0;
        }

        if (unlikely(AllocAndFillSubsetChrFoVidxStart(cur_variant_include, cip, &common.subset_chr_fo_vidx_start))) {
          goto GlmMain_ret_NOMEM;
        }
        common.variant_include = cur_variant_include;
        common.variant_ct = cur_variant_ct;
        linear_ctx.covars_cmaj_d = covars_cmaj_d;

        if (glm_flags & kfGlmPhenoIds) {
          // Possible todo: have file-copy library function which uses
          // sendfile() with Linux kernel 2.6.33+, etc.
          uint32_t pheno_uidx2 = 0;
          for (uint32_t pheno_idx = 0; pheno_idx != batch_size; ++pheno_idx, ++pheno_uidx2) {
            pheno_uidx2 = AdvTo1Bit(pheno_batch, pheno_uidx2);
            char* outname_end2 = strcpya(&(outname_end[1]), &(pheno_names[pheno_uidx2 * max_pheno_name_blen]));
            outname_end2 = strcpya_k(outname_end2, ".glm.linear");
            strcpy_k(outname_end2, ".id");
            reterr = WriteSampleIds(cur_sample_include, siip, outname, sample_ct);
            if (unlikely(reterr)) {
              goto GlmMain_ret_1;
            }
            if (sample_ct_x && x_samples_are_different) {
              strcpy_k(outname_end2, ".x.id");
              reterr = WriteSampleIds(cur_sample_include_x, siip, outname, sample_ct_x);
              if (unlikely(reterr)) {
                goto GlmMain_ret_1;
              }
            }
            if (sample_ct_y && y_samples_are_different) {
              strcpy_k(outname_end2, ".y.id");
              reterr = WriteSampleIds(cur_sample_include_y, siip, outname, sample_ct_y);
              if (unlikely(reterr)) {
                goto GlmMain_ret_1;
              }
            }
          }
        }

        reterr = GlmLinearBatch(pheno_batch, pheno_cols, pheno_names, cur_test_names, cur_test_names_x, cur_test_names_y, glm_pos_col? variant_bps : nullptr, variant_ids, allele_storage, glm_info_ptr, local_sample_uidx_order, cur_local_variant_include, raw_variant_ct, completed_pheno_ct, batch_size, max_pheno_name_blen, max_chr_blen, ci_size, ln_pfilter, output_min_ln, max_thread_ct, pgr_alloc_cacheline_ct, overflow_buf_size, local_sample_ct, pgfip, &linear_ctx, &local_covar_txs, outname, outname_end);
        if (unlikely(reterr)) {
          goto GlmMain_ret_1;
        }
        completed_pheno_ct += batch_size;
        BitvecInvmask(pheno_batch, pheno_ctl, pheno_include);
      }
      for (uint32_t pheno_uidx = 0; pheno_uidx != pheno_ct; ++pheno_uidx) {
        const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_uidx]);
        const PhenoDtype dtype_code = cur_pheno_col->type_code;
        if (dtype_code != kPhenoDtypeQt) {
          SetBit(pheno_uidx, pheno_include);
        }
      }
    }

    for (uint32_t pheno_uidx = 0; pheno_uidx != pheno_ct; ++pheno_uidx) {
      if (!IsSet(pheno_include, pheno_uidx)) {
        continue;
      }
      const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_uidx]);
      const PhenoDtype dtype_code = cur_pheno_col->type_code;
      const char* cur_pheno_name = &(pheno_names[pheno_uidx * max_pheno_name_blen]);
      if (dtype_code == kPhenoDtypeCat) {
        // todo: check if there are only two categories after linear-style
        // covariate QC, and automatically use ordinary logistic regression in
        // that case?  (need to indicate which category is treated as 'case'
        // and which is 'control'...)
        // longer-term todo: multinomial logistic regression?
        logprintfww("--glm: Skipping categorical phenotype '%s'.\n", cur_pheno_name);
        continue;
      }

      BitvecAndCopy(orig_sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include);
      const uint32_t is_logistic = (dtype_code == kPhenoDtypeCc);
      uint32_t sample_ct = PopcountWords(cur_sample_include, raw_sample_ctl);
      if (is_logistic) {
        const uint32_t initial_case_ct = PopcountWordsIntersect(cur_sample_include, cur_pheno_col->data.cc, raw_sample_ctl);
        if ((!initial_case_ct) || (initial_case_ct == sample_ct)) {
          if (unlikely(!skip_modifier)) {
            logerrprintfww("Error: All samples for --glm phenotype '%s' are %s.\n", cur_pheno_name, initial_case_ct? "cases" : "controls");
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logprintfww("--glm: Skipping case/control phenotype '%s' since all samples are %s.\n", cur_pheno_name, initial_case_ct? "cases" : "controls");
          continue;
        }
      } else {
        if (IsConstCovar(cur_pheno_col, cur_sample_include, sample_ct)) {
          if (unlikely(!skip_modifier)) {
            logerrprintfww("Error: --glm quantitative phenotype '%s' is constant.\n", cur_pheno_name);
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logprintfww("--glm: Skipping constant quantitative phenotype '%s'.\n", cur_pheno_name);
          continue;
        }
      }
      uint32_t covar_ct = 0;
      uint32_t extra_cat_ct = 0;
      logistic_ctx.separation_found = 0;
      BigstackDoubleReset(bigstack_mark2, bigstack_end_mark);
      if (initial_nonx_covar_ct) {
        if (unlikely(GlmDetermineCovars(is_logistic? cur_pheno_col->data.cc : nullptr, initial_covar_include, covar_cols, raw_sample_ct, raw_covar_ctl, initial_nonx_covar_ct, covar_max_nonnull_cat_ct, is_sometimes_firth, is_always_firth, cur_sample_include, covar_include, &sample_ct, &covar_ct, &extra_cat_ct, &logistic_ctx.separation_found))) {
          goto GlmMain_ret_NOMEM;
        }
      }
      uint32_t biallelic_predictor_ct = 2 + domdev_present + (covar_ct + extra_cat_ct) * (1 + add_interactions * domdev_present_p1);
      if (raw_parameter_subset) {
        biallelic_predictor_ct = CollapseParamOrTestSubset(covar_include, raw_parameter_subset, domdev_present, raw_covar_ct, covar_ct, add_interactions, common.parameter_subset);
      }
      if (raw_joint_test_params) {
        common.constraint_ct = CollapseParamOrTestSubset(covar_include, raw_joint_test_params, domdev_present, raw_covar_ct, covar_ct, add_interactions, common.joint_test_params);
        if (raw_parameter_subset) {
          // bugfix (2 Feb 2019): forgot to properly take --parameters into
          // account
          memcpy(joint_test_params_buf, common.joint_test_params, biallelic_raw_predictor_ctl * sizeof(intptr_t));
          ZeroWArr(biallelic_raw_predictor_ctl, common.joint_test_params);
          CopyBitarrSubset(joint_test_params_buf, common.parameter_subset, biallelic_predictor_ct, common.joint_test_params);
        }
      }
      if (sample_ct <= biallelic_predictor_ct) {
        if (unlikely(!skip_modifier)) {
          if ((!is_sometimes_firth) && logistic_ctx.separation_found) {
            logerrprintfww("Error: (Quasi-)separated covariate(s) were present for --glm phenotype '%s'. Try removing inappropriate covariates, and/or using Firth logistic regression.\n", cur_pheno_name);
          } else {
            logerrprintfww("Error: # samples <= # predictor columns for --glm phenotype '%s'.\n", cur_pheno_name);
          }
          goto GlmMain_ret_INCONSISTENT_INPUT;
        }
        if ((!is_sometimes_firth) && logistic_ctx.separation_found) {
          logerrprintfww("Warning: Skipping --glm regression on phenotype '%s' since (quasi-)separated covariate(s) were present. Try removing inappropriate covariates, and/or using Firth logistic regression.\n", cur_pheno_name);
        } else {
          logprintfww("Note: Skipping --glm regression on phenotype '%s' since # samples <= # predictor columns.\n", cur_pheno_name);
        }
        continue;
      }
#ifdef __LP64__
      if (RoundUpPow2(sample_ct, 4) * S_CAST(uint64_t, biallelic_predictor_ct + max_extra_allele_ct) > 0x7fffffff) {
        // todo: remove this constraint in LAPACK_ILP64 case?
        if (unlikely(!skip_modifier)) {
          logerrprintfww("Error: Too many samples or predictors for --glm regression on phenotype '%s' (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
          goto GlmMain_ret_INCONSISTENT_INPUT;
        }
        logprintfww("Note: Skipping --glm regression on phenotype '%s' since there are too many samples or predictors (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
        continue;
      }
#endif
      uint32_t case_ct = 0;
      if (is_logistic) {
        case_ct = PopcountWordsIntersect(cur_sample_include, cur_pheno_col->data.cc, raw_sample_ctl);
        if ((!case_ct) || (case_ct == sample_ct)) {
          if (unlikely(!skip_modifier)) {
            logerrprintfww("Error: All remaining samples for --glm phenotype '%s' are %s.\n", cur_pheno_name, case_ct? "cases" : "controls");
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logprintfww("--glm: Skipping case/control phenotype '%s' since all remaining samples are %s.\n", cur_pheno_name, case_ct? "cases" : "controls");
          // without any e.g. cases in the dataset, every single covariate
          // should fail the separation check, so covar_ct should be zero here
          assert(!covar_ct);
          continue;
        }
        // quasi-bugfix (4 Jun 2018): "one in ten" rule of thumb applies to
        // minimum of case and control counts, not total sample count
        if (MINV(case_ct, sample_ct - case_ct) < 10 * biallelic_predictor_ct) {
          if (case_ct * 2 < sample_ct) {
            logerrprintfww("Warning: --glm remaining case count is less than 10x predictor count for phenotype '%s'.\n", cur_pheno_name);
          } else {
            logerrprintfww("Warning: --glm remaining control count is less than 10x predictor count for phenotype '%s'.\n", cur_pheno_name);
          }
        }
      } else {
        // verify phenotype is still nonconstant
        if (IsConstCovar(cur_pheno_col, cur_sample_include, sample_ct)) {
          if (unlikely(!skip_modifier)) {
            logerrprintfww("Error: --glm quantitative phenotype '%s' is constant for all remaining samples.\n", cur_pheno_name);
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logprintfww("--glm: Skipping quantitative phenotype '%s' since phenotype is constant for all remaining samples.\n", cur_pheno_name);
          continue;
        }
      }
      if (common.tests_flag && (!common.constraint_ct)) {
        if (unlikely(!skip_modifier)) {
          logerrprintfww("Error: --tests predictor(s) are constant for all remaining samples for --glm phenotype '%s'.\n", cur_pheno_name);
          goto GlmMain_ret_INCONSISTENT_INPUT;
        }
        logprintfww("Note: Skipping --glm regression on phenotype '%s' since --tests predictor(s) are constant for all remaining samples.\n", cur_pheno_name);
        continue;
      }
      if (covar_ct < initial_nonx_covar_ct) {
        uintptr_t covar_uidx_base = 0;
        uintptr_t cur_bits = initial_covar_include[0];
        for (uint32_t covar_idx = 0; covar_idx != initial_nonx_covar_ct; ++covar_idx) {
          const uintptr_t covar_uidx = BitIter1(initial_covar_include, &covar_uidx_base, &cur_bits);
          if (!IsSet(covar_include, covar_uidx)) {
            logerrprintfww("Warning: %sot including covariate '%s' in --glm regression on phenotype '%s'.\n", cur_sample_include_x_buf? (cur_sample_include_y_buf? "Outside of chrX, n" : "Outside of chrX and chrY, n") : (cur_sample_include_y_buf? "Outside of chrY, n" : "N"), &(covar_names[covar_uidx * max_covar_name_blen]), cur_pheno_name);
          }
        }
      }

      // cur_sample_include_x == nullptr: chrX uses same samples and covariates
      //   as the rest of the genome.  sample_ct_x always zero to force most
      //   chrX-specific initialization to be skipped (exception:
      //   sex_male_collapsed, needed for allele count/freq reporting)
      // cur_sample_include_x non-null: if sample_ct_x == 0, we skip the entire
      //   chromosome.  otherwise, we have different covariates than the rest
      //   of the genome.
      uintptr_t* cur_sample_include_x = cur_sample_include_x_buf;
      uint32_t sample_ct_x = 0;
      uint32_t covar_ct_x = 0;
      uint32_t extra_cat_ct_x = 0;
      uint32_t biallelic_predictor_ct_x = 0;
      uint32_t x_samples_are_different = 0;
      if (cur_sample_include_x) {
        BitvecAndCopy(orig_sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include_x);
        logistic_ctx.separation_found_x = 0;
        if (unlikely(GlmDetermineCovars(is_logistic? cur_pheno_col->data.cc : nullptr, initial_covar_include, covar_cols, raw_sample_ct, raw_covar_ctl, initial_nonx_covar_ct + 1, covar_max_nonnull_cat_ct, is_sometimes_firth, is_always_firth, cur_sample_include_x, covar_include_x, &sample_ct_x, &covar_ct_x, &extra_cat_ct_x, &logistic_ctx.separation_found_x))) {
          goto GlmMain_ret_NOMEM;
        }
        x_samples_are_different = (sample_ct_x != sample_ct) || (!wordsequal(cur_sample_include, cur_sample_include_x, raw_sample_ctl));
        if ((!x_samples_are_different) && (covar_ct == covar_ct_x) && wordsequal(covar_include, covar_include_x, raw_covar_ctl)) {
          logprintfww("Note: chrX samples and covariate(s) in --glm regression on phenotype '%s' are the same as that for the rest of the genome.\n", cur_pheno_name);
          sample_ct_x = 0;
          cur_sample_include_x = nullptr;
        } else {
          if (!sample_ct_x) {
            if (unlikely(!skip_modifier)) {
              logerrprintfww("Error: --glm regression on phenotype '%s' is degenerate on chrX.\n", cur_pheno_name);
              goto GlmMain_ret_INCONSISTENT_INPUT;
            }
            logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s'.\n", cur_pheno_name);
          } else {
            biallelic_predictor_ct_x = 2 + domdev_present + (covar_ct_x + extra_cat_ct_x) * (1 + add_interactions * domdev_present_p1);
            if (raw_parameter_subset) {
              biallelic_predictor_ct_x = CollapseParamOrTestSubset(covar_include, raw_parameter_subset, domdev_present, raw_covar_ct, covar_ct_x, add_interactions, common.parameter_subset_x);
            }
            if (raw_joint_test_params) {
              common.constraint_ct_x = CollapseParamOrTestSubset(covar_include, raw_joint_test_params, domdev_present, raw_covar_ct, covar_ct_x, add_interactions, common.joint_test_params_x);
              if (raw_parameter_subset) {
                memcpy(joint_test_params_buf, common.joint_test_params_x, biallelic_raw_predictor_ctl * sizeof(intptr_t));
                ZeroWArr(biallelic_raw_predictor_ctl, common.joint_test_params_x);
                CopyBitarrSubset(joint_test_params_buf, common.parameter_subset, biallelic_predictor_ct_x, common.joint_test_params_x);
              }
            }
            if (sample_ct_x <= biallelic_predictor_ct_x) {
              if (unlikely(!skip_modifier)) {
                if ((!is_sometimes_firth) && logistic_ctx.separation_found_x) {
                  logerrprintfww("Error: (Quasi-)separated covariate(s) were present for --glm phenotype '%s' on chrX. Try removing inappropriate covariates, and/or using Firth logistic regression.\n", cur_pheno_name);
                } else {
                  logerrprintfww("Error: # samples <= # predictor columns for --glm phenotype '%s' on chrX.\n", cur_pheno_name);
                }
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              if ((!is_sometimes_firth) && logistic_ctx.separation_found_x) {
                logerrprintfww("Warning: Skipping --glm regression on phenotype '%s' on chrX, since (quasi-)separated covariate(s) were present. Try removing inappropriate covariates, and/or using Firth logistic regression.\n", cur_pheno_name);
              } else {
                logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', since # remaining samples <= # predictor columns.\n", cur_pheno_name);
              }
              sample_ct_x = 0;
#ifdef __LP64__
            } else if (RoundUpPow2(sample_ct_x, 4) * S_CAST(uint64_t, biallelic_predictor_ct_x + max_extra_allele_ct) > 0x7fffffff) {
              if (unlikely(!skip_modifier)) {
                logerrprintfww("Error: Too many samples or predictors for --glm regression on phenotype '%s' on chrX (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              logerrprintfww("Warning: Skipping chrX in --glm regression on phenotype '%s', since there are too many samples or predictors (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
              sample_ct_x = 0;
#endif
            } else if (is_logistic) {
              const uint32_t case_ct_x = PopcountWordsIntersect(cur_sample_include_x, cur_pheno_col->data.cc, raw_sample_ctl);
              if ((!case_ct_x) || (case_ct_x == sample_ct_x)) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: All remaining samples on chrX for --glm phenotype '%s' are %s.\n", cur_pheno_name, case_ct? "cases" : "controls");
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', since all remaining samples are %s.\n", cur_pheno_name, case_ct_x? "cases" : "controls");
                sample_ct_x = 0;
              }
            } else {
              if (IsConstCovar(cur_pheno_col, cur_sample_include_x, sample_ct_x)) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: --glm quantitative phenotype '%s' is constant for all remaining samples on chrX.\n", cur_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', since phenotype is constant for all remaining samples.\n", cur_pheno_name);
                sample_ct_x = 0;
              }
            }
            if (sample_ct_x && common.tests_flag && (!common.constraint_ct_x)) {
              if (unlikely(!skip_modifier)) {
                logerrprintfww("Error: --tests predictor(s) are constant for all remaining samples on chrX for --glm phenotype '%s'.\n", cur_pheno_name);
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', since --tests predictor(s) are constant for all remaining samples.\n", cur_pheno_name);
              sample_ct_x = 0;
            }
            if (sample_ct_x && (covar_ct_x < initial_nonx_covar_ct + 1)) {
              uintptr_t covar_uidx_base = 0;
              uintptr_t cur_bits = initial_covar_include[0];
              for (uint32_t covar_idx = 0; covar_idx != covar_ct_x; ++covar_idx) {
                const uintptr_t covar_uidx = BitIter1(initial_covar_include, &covar_uidx_base, &cur_bits);
                if (!IsSet(covar_include_x, covar_uidx)) {
                  logerrprintfww("Warning: On chrX, not including covariate '%s' in --glm regression on phenotype '%s'.\n", &(covar_names[covar_uidx * max_covar_name_blen]), cur_pheno_name);
                }
              }
            }
          }
        }
      }

      uintptr_t* cur_sample_include_y = cur_sample_include_y_buf;
      uint32_t sample_ct_y = 0;
      uint32_t covar_ct_y = 0;
      uint32_t extra_cat_ct_y = 0;
      uint32_t biallelic_predictor_ct_y = 0;
      uint32_t y_samples_are_different = 0;
      if (cur_sample_include_y) {
        BitvecAndCopy(orig_sample_include, sex_male, raw_sample_ctl, cur_sample_include_y);
        BitvecAnd(cur_pheno_col->nonmiss, raw_sample_ctl, cur_sample_include_y);
        logistic_ctx.separation_found_y = 0;
        if (unlikely(GlmDetermineCovars(is_logistic? cur_pheno_col->data.cc : nullptr, initial_covar_include, covar_cols, raw_sample_ct, raw_covar_ctl, initial_y_covar_ct, covar_max_nonnull_cat_ct, is_sometimes_firth, is_always_firth, cur_sample_include_y, covar_include_y, &sample_ct_y, &covar_ct_y, &extra_cat_ct_y, &logistic_ctx.separation_found_y))) {
          goto GlmMain_ret_NOMEM;
        }
        y_samples_are_different = (sample_ct_y != sample_ct) || (!wordsequal(cur_sample_include, cur_sample_include_y, raw_sample_ctl));
        if ((!y_samples_are_different) && (covar_ct == covar_ct_y) && wordsequal(covar_include, covar_include_y, raw_covar_ctl)) {
          logprintfww("Note: chrY samples and covariate(s) in --glm regression on phenotype '%s' are the same as that for the rest of the genome.\n", cur_pheno_name);
          sample_ct_y = 0;
          cur_sample_include_y = nullptr;
        } else {
          if (!sample_ct_y) {
            if (unlikely(!skip_modifier)) {
              logerrprintfww("Error: --glm regression on phenotype '%s' is degenerate on chrY.\n", cur_pheno_name);
              goto GlmMain_ret_INCONSISTENT_INPUT;
            }
            logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s'.\n", cur_pheno_name);
          } else {
            biallelic_predictor_ct_y = 2 + domdev_present + (covar_ct_y + extra_cat_ct_y) * (1 + add_interactions * domdev_present_p1);
            if (raw_parameter_subset) {
              biallelic_predictor_ct_y = CollapseParamOrTestSubset(covar_include, raw_parameter_subset, domdev_present, raw_covar_ct, covar_ct_y, add_interactions, common.parameter_subset_y);
            }
            if (raw_joint_test_params) {
              assert(common.tests_flag);
              common.constraint_ct_y = CollapseParamOrTestSubset(covar_include, raw_joint_test_params, domdev_present, raw_covar_ct, covar_ct_y, add_interactions, common.joint_test_params_y);
              if (raw_parameter_subset) {
                memcpy(joint_test_params_buf, common.joint_test_params_y, biallelic_raw_predictor_ctl * sizeof(intptr_t));
                ZeroWArr(biallelic_raw_predictor_ctl, common.joint_test_params_y);
                CopyBitarrSubset(joint_test_params_buf, common.parameter_subset, biallelic_predictor_ct_y, common.joint_test_params_y);
              }
            }
            if (sample_ct_y <= biallelic_predictor_ct_y) {
              if (unlikely(!skip_modifier)) {
                if ((!is_sometimes_firth) && logistic_ctx.separation_found_x) {
                  logerrprintfww("Error: (Quasi-)separated covariate(s) were present for --glm phenotype '%s' on chrY. Try removing inappropriate covariates, and/or using Firth logistic regression.\n", cur_pheno_name);
                } else {
                  logerrprintfww("Error: # samples <= # predictor columns for --glm phenotype '%s' on chrY.\n", cur_pheno_name);
                }
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              if ((!is_sometimes_firth) && logistic_ctx.separation_found_x) {
                logerrprintfww("Warning: Skipping --glm regression on phenotype '%s' on chrY, since (quasi-)separated covariate(s) were present. Try removing inappropriate covariates, and/or using Firth logistic regression.\n", cur_pheno_name);
              } else {
                logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', since # remaining samples <= # predictor columns.\n", cur_pheno_name);
              }
              sample_ct_y = 0;
#ifdef __LP64__
            } else if (RoundUpPow2(sample_ct_y, 4) * S_CAST(uint64_t, biallelic_predictor_ct_y + max_extra_allele_ct) > 0x7fffffff) {
              if (unlikely(!skip_modifier)) {
                logerrprintfww("Error: Too many samples or predictors for --glm regression on phenotype '%s' on chrY (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              logerrprintfww("Warning: Skipping chrY in --glm regression on phenotype '%s', since there are too many samples or predictors (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
              sample_ct_y = 0;
#endif
            } else if (is_logistic) {
              const uint32_t case_ct_y = PopcountWordsIntersect(cur_sample_include_y, cur_pheno_col->data.cc, raw_sample_ctl);
              if ((!case_ct_y) || (case_ct_y == sample_ct_y)) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: All remaining samples on chrY for --glm phenotype '%s' are %s.\n", cur_pheno_name, case_ct? "cases" : "controls");
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', since all remaining samples are %s.\n", cur_pheno_name, case_ct_y? "cases" : "controls");
                sample_ct_y = 0;
              }
            } else {
              if (IsConstCovar(cur_pheno_col, cur_sample_include_y, sample_ct_y)) {
                if (unlikely(!skip_modifier)) {
                  logerrprintfww("Error: --glm quantitative phenotype '%s' is constant for all remaining samples on chrY.\n", cur_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', since phenotype is constant for all remaining samples.\n", cur_pheno_name);
                sample_ct_y = 0;
              }
            }
            if (sample_ct_y && common.tests_flag && (!common.constraint_ct_y)) {
              if (unlikely(!skip_modifier)) {
                logerrprintfww("Error: --tests predictor(s) are constant for all remaining samples on chrY for --glm phenotype '%s'.\n", cur_pheno_name);
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', since --tests predictor(s) are constant for all remaining samples.\n", cur_pheno_name);
              sample_ct_y = 0;
            }
            if (sample_ct_y && (covar_ct_y < initial_y_covar_ct)) {
              uintptr_t covar_uidx_base = 0;
              uintptr_t cur_bits = initial_covar_include[0];
              for (uint32_t covar_idx = 0; covar_idx != covar_ct_y; ++covar_idx) {
                const uintptr_t covar_uidx = BitIter1(initial_covar_include, &covar_uidx_base, &cur_bits);
                if (!IsSet(covar_include_y, covar_uidx)) {
                  logerrprintfww("Warning: On chrY, not including covariate '%s' in --glm regression on phenotype '%s'.\n", &(covar_names[covar_uidx * max_covar_name_blen]), cur_pheno_name);
                }
              }
            }
          }
        }
      }

      // Expand categorical covariates and perform VIF and correlation checks
      // here.
      double* covars_cmaj_d = nullptr;
      float* pheno_f = nullptr;
      float* covars_cmaj_f = nullptr;
      const char** cur_covar_names = nullptr;
      GlmErr glm_err;
      logistic_ctx.pheno_cc = nullptr;
      logistic_ctx.gcount_case_interleaved_vec = nullptr;
      logistic_ctx.pheno_f = nullptr;
      logistic_ctx.covars_cmaj_f = nullptr;
      linear_ctx.pheno_d = nullptr;
      linear_ctx.covars_cmaj_d = nullptr;
      if (is_logistic) {
        if (unlikely(GlmAllocFillAndTestPhenoCovarsCc(cur_sample_include, cur_pheno_col->data.cc, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &logistic_ctx.pheno_cc, gcount_cc_col? (&logistic_ctx.gcount_case_interleaved_vec) : nullptr, &pheno_f, &common.nm_precomp, &covars_cmaj_f, &cur_covar_names, &glm_err))) {
          goto GlmMain_ret_NOMEM;
        }
      } else {
        if (unlikely(GlmAllocFillAndTestPhenoCovarsQt(cur_sample_include, cur_pheno_col->data.qt, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &linear_ctx.pheno_d, &common.nm_precomp, &covars_cmaj_d, &cur_covar_names, &glm_err))) {
          goto GlmMain_ret_NOMEM;
        }
      }
      if (glm_err) {
        PrintPrescanErrmsg("", cur_pheno_name, cur_covar_names, glm_err, local_covar_ct, skip_modifier, 0);
        if (unlikely(!skip_modifier)) {
          goto GlmMain_ret_INCONSISTENT_INPUT;
        }
        continue;
      }
      const char** cur_covar_names_x = nullptr;
      common.nm_precomp_x = nullptr;
      logistic_ctx.pheno_x_cc = nullptr;
      logistic_ctx.gcount_case_interleaved_vec_x = nullptr;
      logistic_ctx.pheno_x_f = nullptr;
      logistic_ctx.covars_cmaj_x_f = nullptr;
      linear_ctx.pheno_x_d = nullptr;
      linear_ctx.covars_cmaj_x_d = nullptr;
      if (sample_ct_x) {
        if (is_logistic) {
          // bugfix (28 Apr 2018): common.nm_precomp -> common.nm_precomp_x
          if (unlikely(GlmAllocFillAndTestPhenoCovarsCc(cur_sample_include_x, cur_pheno_col->data.cc, covar_include_x, covar_cols, covar_names, sample_ct_x, covar_ct_x, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_x, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &logistic_ctx.pheno_x_cc, gcount_cc_col? (&logistic_ctx.gcount_case_interleaved_vec_x) : nullptr, &logistic_ctx.pheno_x_f, &common.nm_precomp_x, &logistic_ctx.covars_cmaj_x_f, &cur_covar_names_x, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
        } else {
          if (unlikely(GlmAllocFillAndTestPhenoCovarsQt(cur_sample_include_x, cur_pheno_col->data.qt, covar_include_x, covar_cols, covar_names, sample_ct_x, covar_ct_x, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_x, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &linear_ctx.pheno_x_d, &common.nm_precomp_x, &linear_ctx.covars_cmaj_x_d, &cur_covar_names_x, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
        }
        if (glm_err) {
          PrintPrescanErrmsg("chrX in ", cur_pheno_name, cur_covar_names_x, glm_err, local_covar_ct, skip_modifier, 0);
          if (unlikely(!skip_modifier)) {
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          sample_ct_x = 0;
        }
      }
      const char** cur_covar_names_y = nullptr;
      common.nm_precomp_y = nullptr;
      logistic_ctx.pheno_y_cc = nullptr;
      logistic_ctx.gcount_case_interleaved_vec_y = nullptr;
      logistic_ctx.pheno_y_f = nullptr;
      logistic_ctx.covars_cmaj_y_f = nullptr;
      linear_ctx.pheno_y_d = nullptr;
      linear_ctx.covars_cmaj_y_d = nullptr;
      if (sample_ct_y) {
        if (is_logistic) {
          if (unlikely(GlmAllocFillAndTestPhenoCovarsCc(cur_sample_include_y, cur_pheno_col->data.cc, covar_include_y, covar_cols, covar_names, sample_ct_y, covar_ct_y, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_y, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &logistic_ctx.pheno_y_cc, gcount_cc_col? (&logistic_ctx.gcount_case_interleaved_vec_y) : nullptr, &logistic_ctx.pheno_y_f, &common.nm_precomp_y, &logistic_ctx.covars_cmaj_y_f, &cur_covar_names_y, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
        } else {
          if (unlikely(GlmAllocFillAndTestPhenoCovarsQt(cur_sample_include_y, cur_pheno_col->data.qt, covar_include_y, covar_cols, covar_names, sample_ct_y, covar_ct_y, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_y, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &linear_ctx.pheno_y_d, &common.nm_precomp_y, &linear_ctx.covars_cmaj_y_d, &cur_covar_names_y, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
        }
        if (glm_err) {
          PrintPrescanErrmsg("chrY in ", cur_pheno_name, cur_covar_names_y, glm_err, local_covar_ct, skip_modifier, 0);
          if (unlikely(!skip_modifier)) {
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          sample_ct_y = 0;
        }
      }
      const char** cur_test_names = nullptr;
      const char** cur_test_names_x = nullptr;
      const char** cur_test_names_y = nullptr;
      if (unlikely(AllocAndInitReportedTestNames(common.parameter_subset, cur_covar_names, glm_flags, covar_ct + extra_cat_ct, common.tests_flag? common.constraint_ct : 0, &cur_test_names))) {
        goto GlmMain_ret_NOMEM;
      }
      if (sample_ct_x) {
        if (unlikely(AllocAndInitReportedTestNames(common.parameter_subset_x, cur_covar_names_x, glm_flags, covar_ct_x + extra_cat_ct_x, common.tests_flag? common.constraint_ct_x : 0, &cur_test_names_x))) {
          goto GlmMain_ret_NOMEM;
        }
      }
      if (sample_ct_y) {
        if (unlikely(AllocAndInitReportedTestNames(common.parameter_subset_y, cur_covar_names_y, glm_flags, covar_ct_y + extra_cat_ct_y, common.tests_flag? common.constraint_ct_y : 0, &cur_test_names_y))) {
          goto GlmMain_ret_NOMEM;
        }
      }

      // okay, we know what variants we're running the regression on, and we've
      // done much of the necessary covariate preprocessing.  now prepare to
      // launch GlmLogistic()/GlmLinear().

      const uintptr_t* cur_variant_include = early_variant_include;
      const uintptr_t* cur_local_variant_include = local_variant_include;
      const uint32_t skip_x = variant_ct_x && ((!xchr_model) || (cur_sample_include_x && (!sample_ct_x)));
      const uint32_t skip_y = variant_ct_y && ((!male_ct) || (cur_sample_include_y && (!sample_ct_y)));
      uint32_t cur_variant_ct = variant_ct;
      if (skip_x || skip_y) {
        uintptr_t* tmp_variant_include;
        if (unlikely(bigstack_alloc_w(raw_variant_ctl, &tmp_variant_include))) {
          goto GlmMain_ret_NOMEM;
        }
        memcpy(tmp_variant_include, early_variant_include, raw_variant_ctl * sizeof(intptr_t));
        uintptr_t* tmp_local_variant_include = nullptr;
        if (local_variant_include) {
          if (unlikely(bigstack_alloc_w(local_variant_ctl, &tmp_local_variant_include))) {
            goto GlmMain_ret_NOMEM;
          }
          memcpy(tmp_local_variant_include, local_variant_include, local_variant_ctl * sizeof(intptr_t));
        }
        if (skip_x) {
          if (local_variant_include) {
            const uint32_t variant_ct_before_x = PopcountBitRange(early_variant_include, 0, x_start);
            uint32_t local_uidx_first = IdxToUidxBasic(local_variant_include, variant_ct_before_x);
            uint32_t local_uidx_last = FindNth1BitFrom(local_variant_include, local_uidx_first, variant_ct_x);
            ClearBitsNz(local_uidx_first, local_uidx_last + 1, tmp_local_variant_include);
          }
          ClearBitsNz(x_start, x_end, tmp_variant_include);
          cur_variant_ct -= variant_ct_x;
        }
        if (skip_y) {
          if (local_variant_include) {
            const uint32_t variant_ct_before_y = PopcountBitRange(early_variant_include, 0, y_start);
            uint32_t local_uidx_first = IdxToUidxBasic(local_variant_include, variant_ct_before_y);
            uint32_t local_uidx_last = FindNth1BitFrom(local_variant_include, local_uidx_first, variant_ct_y);
            ClearBitsNz(local_uidx_first, local_uidx_last + 1, tmp_local_variant_include);
          }
          ClearBitsNz(y_start, y_end, tmp_variant_include);
          cur_variant_ct -= variant_ct_y;
        }
        cur_variant_include = tmp_variant_include;
        cur_local_variant_include = tmp_local_variant_include;
      }
      if (sex_male_collapsed_buf && (!skip_x)) {
        if (!cur_sample_include_x) {
          CopyBitarrSubset(sex_male, cur_sample_include, sample_ct, sex_male_collapsed_buf);
        } else {
          CopyBitarrSubset(sex_male, cur_sample_include_x, sample_ct_x, sex_male_collapsed_buf);
        }
      }
      // todo: if permutation test, also keep whatever statistic is most
      // appropriate for that
      FillCumulativePopcounts(cur_sample_include, raw_sample_ctl, common.sample_include_cumulative_popcounts);
      common.sample_ct = sample_ct;
      common.sample_ct_x = sample_ct_x;
      common.covar_ct = covar_ct + extra_cat_ct;
      common.local_covar_ct = local_covar_ct;
      if (sample_ct_x) {
        if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &common.sample_include_x_cumulative_popcounts))) {
          goto GlmMain_ret_NOMEM;
        }
        FillCumulativePopcounts(cur_sample_include_x, raw_sample_ctl, common.sample_include_x_cumulative_popcounts);
        common.sample_include_x = cur_sample_include_x;
        common.covar_ct_x = covar_ct_x + extra_cat_ct_x;
        // common.male_ct = PopcountWordsIntersect(cur_sample_include_x, sex_male, raw_sample_ctl);
      } else {
        // technically only need this if variant_ct_x && (!skip_x)
        // common.male_ct = PopcountWordsIntersect(cur_sample_include, sex_male, raw_sample_ctl);

        // defensive
        common.sample_include_x = nullptr;
        common.sample_include_x_cumulative_popcounts = nullptr;
        common.covar_ct_x = 0;
      }
      common.sample_ct_y = sample_ct_y;
      if (sample_ct_y) {
        if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &common.sample_include_y_cumulative_popcounts))) {
          goto GlmMain_ret_NOMEM;
        }
        FillCumulativePopcounts(cur_sample_include_y, raw_sample_ctl, common.sample_include_y_cumulative_popcounts);
        common.sample_include_y = cur_sample_include_y;
        common.covar_ct_y = covar_ct_y + extra_cat_ct_y;
      } else {
        common.sample_include_y = nullptr;
        common.sample_include_y_cumulative_popcounts = nullptr;
        common.covar_ct_y = 0;
      }

      double* orig_ln_pvals = nullptr;
      double* orig_permstat = nullptr;
      if (report_adjust || perms_total) {
        const uintptr_t cur_allele_ct = cur_variant_ct + CountExtraAlleles(cur_variant_include, allele_idx_offsets, 0, raw_variant_ct, 0);
        if (unlikely(bigstack_alloc_d(cur_allele_ct, &orig_ln_pvals))) {
          goto GlmMain_ret_NOMEM;
        }
        memcpy(valid_variants, cur_variant_include, raw_variant_ctl * sizeof(intptr_t));
        ZeroWArr(raw_allele_ctl, valid_alleles);
        if (perms_total) {
          if (is_logistic) {
            if (unlikely(bigstack_alloc_d(cur_allele_ct, &orig_permstat))) {
              goto GlmMain_ret_NOMEM;
            }
          } else {
            // unfortunately, signs are reversed
            orig_permstat = orig_ln_pvals;
          }
        }
      }

      if (unlikely(AllocAndFillSubsetChrFoVidxStart(cur_variant_include, cip, &common.subset_chr_fo_vidx_start))) {
        goto GlmMain_ret_NOMEM;
      }
      common.variant_include = cur_variant_include;
      common.variant_ct = cur_variant_ct;
      // this is safe, see pheno_name_blen_capacity check above
      char* outname_end2 = strcpya(&(outname_end[1]), cur_pheno_name);
      if (is_logistic) {
        logistic_ctx.pheno_f = pheno_f;
        logistic_ctx.covars_cmaj_f = covars_cmaj_f;
        if (is_always_firth) {
          outname_end2 = strcpya_k(outname_end2, ".glm.firth");
        } else if (is_sometimes_firth) {
          outname_end2 = strcpya_k(outname_end2, ".glm.logistic.hybrid");
        } else {
          outname_end2 = strcpya_k(outname_end2, ".glm.logistic");
        }
      } else {
        linear_ctx.covars_cmaj_d = covars_cmaj_d;
        outname_end2 = strcpya_k(outname_end2, ".glm.linear");
      }
      // write IDs
      if (glm_flags & kfGlmPhenoIds) {
        snprintf(outname_end2, 22, ".id");
        reterr = WriteSampleIds(cur_sample_include, siip, outname, sample_ct);
        if (unlikely(reterr)) {
          goto GlmMain_ret_1;
        }
        if (sample_ct_x && x_samples_are_different) {
          // quasi-bugfix (7 Jan 2017): use ".x.id" suffix instead of ".id.x",
          // since the last part of the file extension should indicate format
          snprintf(outname_end2, 22, ".x.id");
          reterr = WriteSampleIds(cur_sample_include_x, siip, outname, sample_ct_x);
          if (unlikely(reterr)) {
            goto GlmMain_ret_1;
          }
        }
        if (sample_ct_y && y_samples_are_different) {
          snprintf(outname_end2, 22, ".y.id");
          reterr = WriteSampleIds(cur_sample_include_y, siip, outname, sample_ct_y);
          if (unlikely(reterr)) {
            goto GlmMain_ret_1;
          }
        }
      }

      if (output_zst) {
        snprintf(outname_end2, 22, ".zst");
      } else {
        *outname_end2 = '\0';
      }

      uintptr_t valid_allele_ct = 0;
      if (is_logistic) {
        reterr = GlmLogistic(cur_pheno_name, cur_test_names, cur_test_names_x, cur_test_names_y, glm_pos_col? variant_bps : nullptr, variant_ids, allele_storage, glm_info_ptr, local_sample_uidx_order, cur_local_variant_include, outname, raw_variant_ct, max_chr_blen, ci_size, ln_pfilter, output_min_ln, max_thread_ct, pgr_alloc_cacheline_ct, overflow_buf_size, local_sample_ct, pgfip, &logistic_ctx, &local_covar_txs, valid_variants, valid_alleles, orig_ln_pvals, orig_permstat, &valid_allele_ct);
      } else {
        reterr = GlmLinear(cur_pheno_name, cur_test_names, cur_test_names_x, cur_test_names_y, glm_pos_col? variant_bps : nullptr, variant_ids, allele_storage, glm_info_ptr, local_sample_uidx_order, cur_local_variant_include, outname, raw_variant_ct, max_chr_blen, ci_size, ln_pfilter, output_min_ln, max_thread_ct, pgr_alloc_cacheline_ct, overflow_buf_size, local_sample_ct, pgfip, &linear_ctx, &local_covar_txs, valid_variants, valid_alleles, orig_ln_pvals, &valid_allele_ct);
      }
      if (unlikely(reterr)) {
        goto GlmMain_ret_1;
      }
      if (report_adjust) {
        reterr = Multcomp(valid_variants, cip, nullptr, variant_bps, variant_ids, valid_alleles, allele_idx_offsets, allele_storage, nullptr, adjust_info_ptr, orig_ln_pvals, nullptr, valid_allele_ct, max_allele_slen, ln_pfilter, output_min_ln, joint_test, max_thread_ct, outname, outname_end2);
        if (unlikely(reterr)) {
          goto GlmMain_ret_1;
        }
      }
      if (perms_total) {
        // todo
        logerrputs("Error: --glm permutation tests are under development.\n");
        reterr = kPglRetNotYetSupported;
        goto GlmMain_ret_1;
        // note that orig_permstat signs are now flipped for linear case
      }
    }
  }
  while (0) {
  GlmMain_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GlmMain_ret_TKSTREAM_FAIL:
    TokenStreamErrPrint("--condition-list file", &tks);
    break;
  GlmMain_ret_PGR_FAIL:
    PgenErrPrint(reterr);
    break;
  GlmMain_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  GlmMain_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  GlmMain_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  GlmMain_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 GlmMain_ret_1:
  CleanupTokenStream2("--condition-list file", &tks, &reterr);
  CleanupTextStream2(local_covar_fname, &local_covar_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
