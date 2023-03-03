// This file is part of PLINK 2.00, copyright (C) 2005-2023 Shaun Purcell,
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

#include "plink2_adjust.h"
#include "plink2_compress_stream.h"
#include "plink2_glm.h"
#include "plink2_glm_linear.h"
#include "plink2_glm_logistic.h"

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

void InitGwasSsf(GwasSsfInfo* gwas_ssf_info_ptr) {
  gwas_ssf_info_ptr->flags = kfGwasSsf0;
  gwas_ssf_info_ptr->a1freq_lower_limit = 0.0;
  gwas_ssf_info_ptr->rsid_mode = kGwasSsfRsidMode0;
  gwas_ssf_info_ptr->fname = nullptr;
  gwas_ssf_info_ptr->list_fname = nullptr;
}

void CleanupGwasSsf(GwasSsfInfo* gwas_ssf_info_ptr) {
  free_cond(gwas_ssf_info_ptr->fname);
  free_cond(gwas_ssf_info_ptr->list_fname);
}

// [-1] = #CHROM (must be first column)
// [0] = POS
// [1] = ID
// [2] = REF
// [3] = ALT
// [4] = A1
// [5] = OMITTED (optional, required to keep multiallelic variants)
// [6] = A1_FREQ
// [7] = TEST
// [8] = OBS_CT
// [9] = BETA/OR
// [10] = SE/LOG(OR)_SE
// [11] = L## (optional)
// [12] = U## (optional)
// [13] = P/LOG10_P
ENUM_U31_DEF_START()
  kGwasSsfColPos = 0,
  kGwasSsfColId,
  kGwasSsfColRef,
  kGwasSsfColAlt,
  kGwasSsfColA1,
  kGwasSsfColOmitted,
  kGwasSsfColA1Freq,
  kGwasSsfColTest,
  kGwasSsfColObsCt,
  kGwasSsfColBetaOr,
  kGwasSsfColSe,
  kGwasSsfColCiLower,
  kGwasSsfColCiUpper,
  kGwasSsfColP,
  kGwasSsfColNull
ENUM_U31_DEF_END(GwasSsfColidx);

FLAGSET_DEF_START()
  kfGwasSsfColset0,
  kfGwasSsfColsetPos = (1 << kGwasSsfColPos),
  kfGwasSsfColsetId = (1 << kGwasSsfColId),
  kfGwasSsfColsetRef = (1 << kGwasSsfColRef),
  kfGwasSsfColsetAlt = (1 << kGwasSsfColAlt),
  kfGwasSsfColsetA1 = (1 << kGwasSsfColA1),
  kfGwasSsfColsetOmitted = (1 << kGwasSsfColOmitted),
  kfGwasSsfColsetA1Freq = (1 << kGwasSsfColA1Freq),
  kfGwasSsfColsetTest = (1 << kGwasSsfColTest),
  kfGwasSsfColsetObsCt = (1 << kGwasSsfColObsCt),
  kfGwasSsfColsetBetaOr = (1 << kGwasSsfColBetaOr),
  kfGwasSsfColsetSe = (1 << kGwasSsfColSe),
  kfGwasSsfColsetCiLower = (1 << kGwasSsfColCiLower),
  kfGwasSsfColsetCiUpper = (1 << kGwasSsfColCiUpper),
  kfGwasSsfColsetP = (1 << kGwasSsfColP),

  kfGwasSsfColsetRequired = kfGwasSsfColsetPos | kfGwasSsfColsetId | kfGwasSsfColsetRef | kfGwasSsfColsetAlt | kfGwasSsfColsetA1 | kfGwasSsfColsetA1Freq | kfGwasSsfColsetTest | kfGwasSsfColsetObsCt | kfGwasSsfColsetBetaOr | kfGwasSsfColsetSe | kfGwasSsfColsetP
FLAGSET_DEF_END(GwasSsfColFlags);

// If rsid_mode == kGwasSsfRsidModeInfer and force_rsid is false, we initially
// proceed under the assumption that rsID is absent.  Then, if/when we
// encounter a rsID (this should happen quickly if rsIDs are present at all),
// we return kPglRetRetry, and the caller retries with force_rsid true.
PglErr GwasSsfInternal(const GwasSsfInfo* gsip, const char* in_fname, const char* out_fname, uint32_t force_rsid, uint32_t max_line_blen, uint32_t max_thread_ct, TextStream* txsp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  CompressStreamState css;
  PreinitCstream(&css);
  {
    char* line_iter = TextLineEnd(txsp);
    ++line_idx;
    if (!TextGetUnsafe2(txsp, &line_iter)) {
      if (unlikely(TextStreamErrcode2(txsp, &reterr))) {
        snprintf(g_logbuf, kLogbufSize, "Error: --gwas-ssf: %s is empty.\n", in_fname);
        goto GwasSsfInternal_ret_MALFORMED_INPUT_WW;
      }
      goto GwasSsfInternal_ret_TSTREAM_FAIL;
    }
    uint32_t col_skips[kGwasSsfColNull];
    GwasSsfColidx col_types[kGwasSsfColNull];
    char* first_token_end = FirstSpaceOrEoln(line_iter);
    const uint32_t first_token_slen = first_token_end - line_iter;
    GwasSsfColFlags header_cols = kfGwasSsfColset0;
    uint32_t relevant_postchr_col_ct = 0;
    uint32_t is_odds_ratio = 0;
    if (unlikely(!strequal_k(line_iter, "#CHROM", first_token_slen))) {
      if (strequal_k(line_iter, "#POS", first_token_slen) ||
          strequal_k(line_iter, "#ID", first_token_slen)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --gwas-ssf: %s does not have a #CHROM column.\n", in_fname);
        goto GwasSsfInternal_ret_MALFORMED_INPUT_WW;
      }
      if (strequal_k(line_iter, "CHR", first_token_slen)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --gwas-ssf: %s appears to be a PLINK 1.x --linear/--logistic output file, which is missing the required A1_FREQ output column. Redo the regression with PLINK 2 --glm.\n", in_fname);
        goto GwasSsfInternal_ret_MALFORMED_INPUT_WW;
      }
      snprintf(g_logbuf, kLogbufSize, "Error: --gwas-ssf: %s does not appear to be a PLINK 2 --glm output file.\n", in_fname);
      goto GwasSsfInternal_ret_MALFORMED_INPUT_WW;
    }

    const char* token_end = first_token_end;
    for (uint32_t col_idx = 1; ; ++col_idx) {
      const char* linebuf_iter = FirstNonTspace(token_end);
      if (IsEolnKns(*linebuf_iter)) {
        break;
      }
      token_end = CurTokenEnd(linebuf_iter);
      const uint32_t token_slen = token_end - linebuf_iter;
      GwasSsfColidx cur_colidx = kGwasSsfColNull;
      if (token_slen <= 3) {
        if (token_slen == 3) {
          if (memequal_sk(linebuf_iter, "POS")) {
            cur_colidx = kGwasSsfColPos;
          } else if (memequal_sk(linebuf_iter, "REF")) {
            cur_colidx = kGwasSsfColRef;
          } else if (memequal_sk(linebuf_iter, "ALT")) {
            cur_colidx = kGwasSsfColAlt;
          } else if (IsDigit(linebuf_iter[1]) && IsDigit(linebuf_iter[2])) {
            if (linebuf_iter[0] == 'L') {
              cur_colidx = kGwasSsfColCiLower;
            } else if (linebuf_iter[0] == 'U') {
              cur_colidx = kGwasSsfColCiUpper;
            }
          }
        } else if (token_slen == 2) {
          if (memequal_sk(linebuf_iter, "ID")) {
            cur_colidx = kGwasSsfColId;
          } else if (memequal_sk(linebuf_iter, "A1")) {
            cur_colidx = kGwasSsfColA1;
          } else if (memequal_sk(linebuf_iter, "OR")) {
            cur_colidx = kGwasSsfColBetaOr;
            is_odds_ratio = 1;
          } else if (memequal_sk(linebuf_iter, "SE")) {
            cur_colidx = kGwasSsfColSe;
          }
        } else if (strequal_k(linebuf_iter, "P", token_slen)) {
          cur_colidx = kGwasSsfColP;
        }
      } else if (strequal_k(linebuf_iter, "OMITTED", token_slen)) {
        cur_colidx = kGwasSsfColOmitted;
      } else if (strequal_k(linebuf_iter, "A1_FREQ", token_slen)) {
        cur_colidx = kGwasSsfColA1Freq;
      } else if (strequal_k(linebuf_iter, "TEST", token_slen)) {
        cur_colidx = kGwasSsfColTest;
      } else if (strequal_k(linebuf_iter, "OBS_CT", token_slen)) {
        cur_colidx = kGwasSsfColObsCt;
      } else if (strequal_k(linebuf_iter, "BETA", token_slen)) {
        cur_colidx = kGwasSsfColBetaOr;
      } else if (strequal_k(linebuf_iter, "LOG(OR)_SE", token_slen)) {
        cur_colidx = kGwasSsfColSe;
      } else if (strequal_k(linebuf_iter, "LOG10_P", token_slen)) {
        cur_colidx = kGwasSsfColP;
      }
      if (cur_colidx != kGwasSsfColNull) {
        const GwasSsfColFlags cur_colset = S_CAST(GwasSsfColFlags, 1U << cur_colidx);
        if (unlikely(header_cols & cur_colset)) {
          snprintf(g_logbuf, kLogbufSize, "Error: --gwas-ssf: Conflicting columns in header line of %s .\n", in_fname);
          goto GwasSsfInternal_ret_MALFORMED_INPUT_WW;
        }
        header_cols |= cur_colset;
        col_skips[relevant_postchr_col_ct] = col_idx;
        col_types[relevant_postchr_col_ct++] = cur_colidx;
      }
    }
    for (uint32_t rpc_col_idx = relevant_postchr_col_ct - 1; rpc_col_idx; --rpc_col_idx) {
      col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
    }
    if (unlikely((header_cols & kfGwasSsfColsetRequired) != kfGwasSsfColsetRequired)) {
      snprintf(g_logbuf, kLogbufSize, "Error: --gwas-ssf: %s does not have all required input columns. (CHROM, POS, ID, REF, ALT, A1, A1_FREQ, TEST, OBS_CT, BETA/OR, SE/LOG(OR)_SE, and P are required.)\n", in_fname);
      goto GwasSsfInternal_ret_MALFORMED_INPUT_WW;
    }

    const GwasSsfRsidMode rsid_mode = force_rsid? kGwasSsfRsidModeYes : gsip->rsid_mode;
    const uint32_t rsid_col = (rsid_mode == kGwasSsfRsidModeYes);
    const uint32_t output_zst = (gsip->flags / kfGwasSsfZs) & 1;

    // Output line length won't be more than ~twice input line length.
    // (It can be close to 2x because of the variant_id output column.)
    const uintptr_t overflow_buf_size = kCompressStreamBlock + 512 + (2 * k1LU) * max_line_blen;
    reterr = InitCstreamAlloc(out_fname, 0, output_zst, MAXV(1, max_thread_ct - 1), overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto GwasSsfInternal_ret_1;
    }
    cswritep = strcpya_k(cswritep, "chromosome\tbase_pair_location\teffect_allele\tother_allele\t");
    if (!is_odds_ratio) {
      cswritep = strcpya_k(cswritep, "beta");
    } else {
      cswritep = strcpya_k(cswritep, "odds_ratio");
    }
    cswritep = strcpya_k(cswritep, "\tstandard_error\teffect_allele_frequency\tp_value\tvariant_id");
    if (rsid_col) {
      cswritep = strcpya_k(cswritep, "\trsid");
    }
    if (header_cols & kfGwasSsfColsetCiUpper) {
      cswritep = strcpya_k(cswritep, "\tci_upper");
    }
    if (header_cols & kfGwasSsfColsetCiLower) {
      cswritep = strcpya_k(cswritep, "\tci_lower");
    }
    cswritep = strcpya_k(cswritep, "\tn\tref_allele" EOLN_STR);

    const uint8_t acgt_bool_table[256] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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
    // X, XY, PAR1, PAR2 -> 23
    // Y -> 24
    // MT -> 25, not 26
    const unsigned char gwas_ssf_chr_remap[6] = {23, 24, 23, 25, 23, 23};
    const double a1freq_lower_limit = gsip->a1freq_lower_limit;
    char a1freq_lower_limit_str[16]; // > kMaxDoubleGSlen
    uint32_t a1freq_lower_limit_slen = 0;
    if (a1freq_lower_limit > 0.0) {
      char* a1freq_lower_limit_end = dtoa_g(a1freq_lower_limit, a1freq_lower_limit_str);
      a1freq_lower_limit_slen = a1freq_lower_limit_end - a1freq_lower_limit_str;
    }
    ++line_idx;
    line_iter = AdvPastDelim(line_iter, '\n');
    for (; TextGetUnsafe2(txsp, &line_iter); ++line_idx) {
      if (IsEolnKns(*line_iter)) {
        goto GwasSsfInternal_ret_MISSING_TOKENS;
      }
      const uint32_t chr_code_raw = GetChrCodeRaw(line_iter);
      if ((!chr_code_raw) || ((chr_code_raw > 26) && (chr_code_raw < kMaxContigs))) {
        line_iter = AdvPastDelim(line_iter, '\n');
        continue;
      }

      char* token_ptrs[14];
      uint32_t token_slens[14];
      line_iter = TokenLex(line_iter, R_CAST(uint32_t*, col_types), col_skips, relevant_postchr_col_ct, token_ptrs, token_slens);
      if (unlikely(!line_iter)) {
        goto GwasSsfInternal_ret_MISSING_TOKENS;
      }
      line_iter = AdvPastDelim(line_iter, '\n');
      // skip if TEST isn't ADD, or result is NA, or REF allele has non-ACGT
      // character, or OMITTED column is absent and variant is multiallelic
      if (!strequal_k(token_ptrs[7], "ADD", token_slens[7])) {
        continue;
      }
      if ((token_ptrs[13][0] & 0xdf) == 'N') {
        continue;
      }
      const char* ref_allele = token_ptrs[2];
      const uint32_t ref_allele_slen = token_slens[2];
      {
        uint32_t uii = 0;
        for (; uii != ref_allele_slen; ++uii) {
          if (!acgt_bool_table[S_CAST(unsigned char, ref_allele[uii])]) {
            break;
          }
        }
        if (uii != ref_allele_slen) {
          continue;
        }
      }
      if (!(header_cols & kfGwasSsfColsetOmitted)) {
        if (memchr(token_ptrs[3], ',', token_slens[3]) != nullptr) {
          continue;
        }
      }

      uint32_t gwas_ssf_chr_code;
      if (chr_code_raw < kMaxContigs) {
        if (chr_code_raw < 25) {
          gwas_ssf_chr_code = chr_code_raw;
        } else if (chr_code_raw == 25) {
          gwas_ssf_chr_code = 23;
        } else {
          gwas_ssf_chr_code = 25;
        }
      } else {
        gwas_ssf_chr_code = gwas_ssf_chr_remap[chr_code_raw - kMaxContigs];
      }
      cswritep = u32toa_x(gwas_ssf_chr_code, '\t', cswritep);

      // POS
      const char* pos_str = token_ptrs[0];
      const uint32_t pos_slen = token_slens[0];
      cswritep = memcpyax(cswritep, pos_str, pos_slen, '\t');

      // A1
      const char* a1_str = token_ptrs[4];
      const uint32_t a1_slen = token_slens[4];
      cswritep = memcpyax(cswritep, a1_str, a1_slen, '\t');

      // OMITTED
      const char* alt_allele = token_ptrs[3];
      const uint32_t alt_allele_slen = token_slens[3];
      if (header_cols & kfGwasSsfColsetOmitted) {
        cswritep = memcpya(cswritep, token_ptrs[5], token_slens[5]);
      } else {
        // we already verified variant is biallelic.
        // check for REF/ALT match, ALT first (much more likely), error out if
        // neither.
        if ((a1_slen == alt_allele_slen) && memequal(a1_str, alt_allele, a1_slen)) {
          cswritep = memcpya(cswritep, ref_allele, ref_allele_slen);
        } else if (likely((a1_slen == ref_allele_slen) && memequal(a1_str, ref_allele, a1_slen))) {
          cswritep = memcpya(cswritep, alt_allele, alt_allele_slen);
        } else {
          snprintf(g_logbuf, kLogbufSize, "Error: A1 allele on line %" PRIuPTR " of %s matches neither REF nor ALT.\n", line_idx, in_fname);
          goto GwasSsfInternal_ret_MALFORMED_INPUT_WW;
        }
      }
      *cswritep++ = '\t';

      // BETA/OR
      cswritep = memcpyax(cswritep, token_ptrs[9], token_slens[9], '\t');

      // SE
      cswritep = memcpyax(cswritep, token_ptrs[10], token_slens[10], '\t');

      // A1_FREQ
      const char* a1_freq_str = token_ptrs[6];
      if (a1freq_lower_limit == 0.0) {
        cswritep = memcpya(cswritep, a1_freq_str, token_slens[6]);
      } else {
        double a1_freq;
        if (unlikely(!ScanadvDouble(a1_freq_str, &a1_freq))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid A1_FREQ on line %" PRIuPTR " of %s .\n", line_idx, in_fname);
          goto GwasSsfInternal_ret_MALFORMED_INPUT_WW;
        }
        if (a1_freq >= a1freq_lower_limit) {
          cswritep = memcpya(cswritep, a1_freq_str, token_slens[6]);
        } else {
          cswritep = memcpya(cswritep, a1freq_lower_limit_str, a1freq_lower_limit_slen);
        }
      }
      *cswritep++ = '\t';

      // P
      cswritep = memcpyax(cswritep, token_ptrs[13], token_slens[13], '\t');

      // C_P_R_A
      cswritep = u32toa_x(gwas_ssf_chr_code, '_', cswritep);
      cswritep = memcpyax(cswritep, pos_str, pos_slen, '_');
      cswritep = memcpyax(cswritep, ref_allele, ref_allele_slen, '_');
      cswritep = memcpyax(cswritep, alt_allele, alt_allele_slen, '\t');

      if (rsid_mode != kGwasSsfRsidModeNo) {
        const char* variant_id = token_ptrs[1];
        const uint32_t variant_id_slen = token_slens[1];
        uint32_t is_rsid = 0;
        if (StrStartsWith(variant_id, "rs", variant_id_slen)) {
          uint32_t uii = strlen("rs");
          for (; uii != variant_id_slen; ++uii) {
            if (!IsDigit(variant_id[uii])) {
              break;
            }
          }
          is_rsid = (uii == variant_id_slen);
        }
        if (rsid_col) {
          if (is_rsid) {
            cswritep = memcpyax(cswritep, variant_id, variant_id_slen, '\t');
          } else {
            cswritep = strcpya_k(cswritep, "#NA\t");
          }
        } else if (is_rsid) {
          reterr = kPglRetRetry;
          goto GwasSsfInternal_ret_1;
        }
      }

      if (header_cols & kfGwasSsfColsetCiUpper) {
        cswritep = memcpyax(cswritep, token_ptrs[12], token_slens[12], '\t');
      }
      if (header_cols & kfGwasSsfColsetCiLower) {
        cswritep = memcpyax(cswritep, token_ptrs[11], token_slens[11], '\t');
      }

      // OBS_CT
      cswritep = memcpyax(cswritep, token_ptrs[8], token_slens[8], '\t');

      cswritep = memcpya(cswritep, ref_allele, ref_allele_slen);

      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(&css, &cswritep))) {
        goto GwasSsfInternal_ret_WRITE_FAIL;
      }
    }
    if (unlikely(TextStreamErrcode2(txsp, &reterr))) {
      goto GwasSsfInternal_ret_TSTREAM_FAIL;
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto GwasSsfInternal_ret_WRITE_FAIL;
    }
    reterr = kPglRetSuccess;
  }
  while (0) {
  GwasSsfInternal_ret_TSTREAM_FAIL:
    TextStreamErrPrint(in_fname, txsp);
    break;
  GwasSsfInternal_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  GwasSsfInternal_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, in_fname);
  GwasSsfInternal_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  }
 GwasSsfInternal_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr GwasSsfOneFile(const GwasSsfInfo* gsip, const char* in_fname, uint32_t max_thread_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* out_fname = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    uint32_t shared_fname_slen = strlen(in_fname);
    if (StrEndsWith(in_fname, ".zst", shared_fname_slen)) {
      shared_fname_slen -= 4;
    }
    // ".ssf.tsv.zst"
    if (unlikely(bigstack_alloc_c(shared_fname_slen + 13, &out_fname))) {
      goto GwasSsfOneFile_ret_NOMEM;
    }
    const uint32_t output_zst = (gsip->flags / kfGwasSsfZs) & 1;
    char* out_fname_iter = memcpya(out_fname, in_fname, shared_fname_slen);
    out_fname_iter = strcpya_k(out_fname_iter, ".ssf.tsv");
    if (output_zst) {
      out_fname_iter = strcpya_k(out_fname_iter, ".zst");
    }
    *out_fname_iter = '\0';

    uint32_t max_line_blen;
    if (unlikely(StandardizeMaxLineBlen(bigstack_left() / 4, &max_line_blen))) {
      goto GwasSsfOneFile_ret_NOMEM;
    }
    // Invert overflow_buf_size expression so we can assume output lines aren't
    // longer than kMaxLongLine.
    if (max_line_blen > (kMaxLongLine - 512) / 2) {
      max_line_blen = (kMaxLongLine - 512) / 2;
    }
    reterr = InitTextStream(in_fname, max_line_blen, output_zst? 1 : MAXV(1, max_thread_ct - 1), &txs);
    if (unlikely(reterr)) {
      goto GwasSsfOneFile_ret_TSTREAM_FAIL;
    }

    reterr = GwasSsfInternal(gsip, in_fname, out_fname, 0, max_line_blen, max_thread_ct, &txs);
    if (reterr) {
      if (unlikely(reterr != kPglRetRetry)) {
        goto GwasSsfOneFile_ret_1;
      }
      reterr = TextRewind(&txs);
      if (unlikely(reterr)) {
        goto GwasSsfOneFile_ret_TSTREAM_FAIL;
      }
      reterr = GwasSsfInternal(gsip, in_fname, out_fname, 1, max_line_blen, max_thread_ct, &txs);
      if (unlikely(reterr)) {
        goto GwasSsfOneFile_ret_1;
      }
    }
    logprintfww("--gwas-ssf: %s written.\n", out_fname);
  }
  while (0) {
  GwasSsfOneFile_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  GwasSsfOneFile_ret_TSTREAM_FAIL:
    TextStreamErrPrint(in_fname, &txs);
    break;
  }
 GwasSsfOneFile_ret_1:
  CleanupTextStream2(in_fname, &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr GwasSsfStandalone(const GwasSsfInfo* gsip, uint32_t max_thread_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  TextStream list_txs;
  PreinitTextStream(&list_txs);
  {
    const GwasSsfFlags flags = gsip->flags;
    uintptr_t file_ct = 0;
    if (gsip->fname) {
      reterr = GwasSsfOneFile(gsip, gsip->fname, max_thread_ct);
      if (unlikely(reterr)) {
        goto GwasSsfStandalone_ret_1;
      }
      ++file_ct;
    }
    if (gsip->list_fname) {
      reterr = InitTextStream(gsip->list_fname, kTextStreamBlenFast, 1, &list_txs);
      if (unlikely(reterr)) {
        goto GwasSsfStandalone_ret_TSTREAM_FAIL;
      }
      while (1) {
        char* first_token_start = TextGet(&list_txs);
        if (!first_token_start) {
          break;
        }
        if (IsSpaceOrEoln(*first_token_start)) {
          continue;
        }
        char* first_token_end = CurTokenEnd(first_token_start);
        *first_token_end = '\0';
        reterr = GwasSsfOneFile(gsip, first_token_start, max_thread_ct);
        if (unlikely(reterr)) {
          goto GwasSsfStandalone_ret_1;
        }
        ++file_ct;
      }
      if (unlikely(TextStreamErrcode2(&list_txs, &reterr))) {
        goto GwasSsfStandalone_ret_TSTREAM_FAIL;
      }
    }
    logprintf("--gwas-ssf file=/file-list=: %" PRIuPTR " file%s processed.\n", file_ct, (file_ct == 1)? "" : "s");
  }
  while (0) {
  GwasSsfStandalone_ret_TSTREAM_FAIL:
    TextStreamErrPrint(gsip->list_fname, &list_txs);
    break;
  }
 GwasSsfStandalone_ret_1:
  CleanupTextStream2(gsip->list_fname, &list_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr GlmLocalOpen(const char* local_covar_fname, const char* local_pvar_fname, const char* local_psam_fname, const SampleIdInfo* siip, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const GlmInfo* glm_info_ptr, uint32_t raw_sample_ct, uint32_t raw_variant_ct, const uintptr_t** sample_include_ptr, const uintptr_t** sex_nm_ptr, const uintptr_t** sex_male_ptr, const uintptr_t** variant_include_ptr, uint32_t* sample_ct_ptr, uint32_t* variant_ct_ptr, TextStream* local_covar_txsp, uint32_t** local_sample_uidx_order_ptr, uintptr_t** local_variant_include_ptr, uint32_t* local_sample_ct_ptr, uint32_t* local_variant_ctl_ptr, uint32_t* local_covar_ct_ptr) {
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
    //    local_sample_uidx_order
    uint32_t max_line_blen;
    if (unlikely(StandardizeMaxLineBlen(bigstack_left() / 4, &max_line_blen))) {
      goto GlmLocalOpen_ret_NOMEM;
    }
    txs_fname = local_psam_fname;
    reterr = InitTextStreamEx(local_psam_fname, 1, kMaxLongLine, max_line_blen, 1, &txs);
    if (unlikely(reterr)) {
      goto GlmLocalOpen_ret_TSTREAM_FAIL;
    }
    XidMode xid_mode;
    char* line_start;
    reterr = LoadXidHeader("glm local-psam=", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeader0 : kfXidHeaderIgnoreSid, &line_idx, &txs, &xid_mode, &line_start);
    if (unlikely(reterr)) {
      goto GlmLocalOpen_ret_TSTREAM_XID_FAIL;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    char* sorted_xidbox;
    uint32_t* xid_map;
    uintptr_t max_xid_blen;
    // IID -> multiple IID+SID not supported for now, may want to add this
    reterr = SortedXidboxInitAllocEnd(*sample_include_ptr, siip, orig_sample_ct, 0, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto GlmLocalOpen_ret_1;
    }
    uintptr_t* new_sample_include;
    char* idbuf;
    if (unlikely(bigstack_end_calloc_w(raw_sample_ctl, &new_sample_include) ||
                 bigstack_end_alloc_c(max_xid_blen, &idbuf))) {
      goto GlmLocalOpen_ret_NOMEM;
    }
    uint32_t* local_sample_uidx_order = R_CAST(uint32_t*, g_bigstack_base);
    uintptr_t max_local_sample_ct = RoundDownPow2(bigstack_left(), kCacheline) / sizeof(int32_t);
#ifdef __LP64__
    if (max_local_sample_ct > kMaxLongLine / 2) {
      max_local_sample_ct = kMaxLongLine / 2;
    }
#endif
    uintptr_t local_sample_ct = 0;
    if (line_start[0] == '#') {
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
      if (!SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &read_ptr, &sample_uidx, idbuf)) {
        if (unlikely(IsSet(new_sample_include, sample_uidx))) {
          TabsToSpaces(idbuf);
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
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &sample_include_copy) ||
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
      uint32_t is_header_line;
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
      uint32_t max_local_variant_ct = kPglMaxVariantCt;
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
          if (max_local_variant_ct == kPglMaxVariantCt) {
            snprintf(g_logbuf, kLogbufSize, "Error: Too many variants in %s.\n", local_pvar_fname);
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
          const uint32_t first_variant_uidx_in_chr = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[cur_chr_code_u]];
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
          // bugfix (16 Mar 2022): if a variant filter was applied to the main
          // dataset, the local-pvar entry corresponding to prev_variant_uidx
          // may be on a later line, even if we fail to find an ID match for
          // the current local-pvar ID.
          // So we no longer advance prev_variant_uidx/prev_bp here under any
          // circumstances.  We just scan through the same-bp variant IDs,
          // looking for a match, and then advance to the next local-pvar line
          // with prev_variant_uidx/prev_bp unchanged.
          uint32_t scan_variant_uidx = prev_variant_uidx;
          uint32_t scan_bp;
          do {
            const char* loaded_variant_id = variant_ids[scan_variant_uidx];
            if (memequal(cur_variant_id, loaded_variant_id, variant_id_blen)) {
              if (unlikely(IsSet(new_variant_include, scan_variant_uidx))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate ID (with duplicate CHROM/POS) '%s' in %s.\n", cur_variant_id, local_pvar_fname);
                goto GlmLocalOpen_ret_MALFORMED_INPUT_WW;
              }
              SetBit(scan_variant_uidx, new_variant_include);
              ++new_variant_ct;
              SetBit(local_variant_ct, local_variant_include);
              break;
            }
            scan_variant_uidx = AdvBoundedTo1Bit(orig_variant_include, scan_variant_uidx + 1, raw_variant_ct);
            if (scan_variant_uidx >= chr_end) {
              break;
            }
            scan_bp = variant_bps[scan_variant_uidx];
          } while (cur_bp_u == scan_bp);
        }
        continue;
      GlmLocalOpen_skip_variant_and_update_chr:
        if (prev_variant_uidx == raw_variant_ct) {
          break;
        }
        chr_fo_idx = GetVariantChrFoIdx(cip, prev_variant_uidx);
        prev_chr_code = cip->chr_file_order[chr_fo_idx];
        prev_bp = variant_bps[prev_variant_uidx];
        chr_end = cip->chr_fo_vidx_start[cip->chr_idx_to_foidx[prev_chr_code] + 1];
      }
      if (!line_start) {
        if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
          goto GlmLocalOpen_ret_TSTREAM_FAIL;
        }
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
  GlmLocalOpen_ret_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&txs)) {
      break;
    }
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
  if (case_and_ctrl_cat_ct * 2 == pheno_by_cat_ct) {
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
      if ((covar_cols[covar_uidx].type_code != kPhenoDtypeOther) && IsConstCovar(&(covar_cols[covar_uidx]), cur_sample_include, raw_sample_ct)) {
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
              if (!IsConstCovar(cur_covar_col, cur_sample_include, raw_sample_ct)) {
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
        if (unlikely(bigstack_alloc_w(max_cat_ctl, &cat_one_obs) ||
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

void PrintPrescanErrmsg(const char* domain_str, const char* pheno_name, const char** covar_names, GlmErr glm_err, uint32_t local_covar_ct, uint32_t skip_invalid_pheno, uint32_t is_batch) {
  const GlmErrcode errcode = GetGlmErrCode(glm_err);
  const char* msg_start = skip_invalid_pheno? "Note: Skipping " : "Error: Cannot proceed with ";
  if (errcode == kGlmErrcodeUnstableScale) {
    snprintf(g_logbuf, kLogbufSize, "%s%s--glm regression on phenotype '%s'%s, since genotype/covariate scales vary too widely for numerical stability of the current implementation. Try rescaling your covariates with e.g. --covar-variance-standardize.\n", msg_start, domain_str, pheno_name, is_batch? ", and other(s) with identical missingness patterns" : "");
  } else if (errcode == kGlmErrcodeVifInfinite) {
    snprintf(g_logbuf, kLogbufSize, "%s%s--glm regression on phenotype '%s'%s, since covariate correlation matrix could not be inverted (VIF_INFINITE).%s\n", msg_start, domain_str, pheno_name, is_batch? ", and other(s) with identical missingness patterns" : "", domain_str[0]? "" : " You may want to remove redundant covariates and try again.");
  } else if (errcode == kGlmErrcodeLogisticConvergeFail) {
    snprintf(g_logbuf, kLogbufSize, "%s%s--glm regression on phenotype '%s', since covariate-only logistic regression failed to converge, and Firth-fallback was disabled.\n", msg_start, domain_str, pheno_name);
  } else if (errcode == kGlmErrcodeFirthConvergeFail) {
    snprintf(g_logbuf, kLogbufSize, "%s%s--glm regression on phenotype '%s', since covariate-only Firth regression failed to converge.\n", msg_start, domain_str, pheno_name);
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
  if (skip_invalid_pheno) {
    logputsb();
  } else {
    logerrputsb();
  }
}

BoolErr AllocAndInitReportedTestNames(const uintptr_t* parameter_subset, const char* const* covar_names, GlmFlags glm_flags, uint32_t covar_ct, uint32_t user_constraint_ct, const char*** cur_test_names_ptr) {
  const uint32_t model_dominant = (glm_flags / kfGlmDominant) & 1;
  const uint32_t model_recessive = (glm_flags / kfGlmRecessive) & 1;
  const uint32_t model_hetonly = (glm_flags / kfGlmHetonly) & 1;
  const uint32_t is_hethom = (glm_flags / kfGlmHethom) & 1;
  const uint32_t domdev_present = (glm_flags & kfGlmGenotypic) || is_hethom;
  char main_effect[4];
  if (model_dominant) {
    memcpy(main_effect, "DOMx", 4);
  } else if (model_recessive) {
    memcpy(main_effect, "RECx", 4);
  } else if (model_hetonly) {
    memcpy(main_effect, "HETx", 4);
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
    if (unlikely(bigstack_alloc_kcp(biallelic_reported_test_ct, cur_test_names_ptr) ||
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
  if (unlikely(bigstack_alloc_kcp(biallelic_reported_test_ct, cur_test_names_ptr) ||
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
  // bugfix (16 Aug 2021): sex + interaction?
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

static const double kSexMaleToCovarD[2] = {2.0, 1.0};

void SexInteractionReshuffle(uint32_t first_interaction_pred_uidx, uint32_t raw_covar_ct, uint32_t domdev_present, uint32_t biallelic_raw_predictor_ctl, uintptr_t* __restrict parameters_or_tests, uintptr_t* __restrict parameter_subset_reshuffle_buf) {
  ZeroWArr(biallelic_raw_predictor_ctl, parameter_subset_reshuffle_buf);
  CopyBitarrRange(parameters_or_tests, 0, 0, first_interaction_pred_uidx - 1, parameter_subset_reshuffle_buf);
  // bugfix (16 Aug 2021): raw_covar_ct includes sex
  const uint32_t raw_nonsex_interaction_ct = (raw_covar_ct - 1) * (domdev_present + 1);
  CopyBitarrRange(parameters_or_tests, first_interaction_pred_uidx - 1, first_interaction_pred_uidx, raw_nonsex_interaction_ct, parameter_subset_reshuffle_buf);
  const uint32_t first_sex_parameter_idx = first_interaction_pred_uidx - 1 + raw_nonsex_interaction_ct;
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

PglErr GlmMain(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const char* const* allele_storage, const GlmInfo* glm_info_ptr, const AdjustInfo* adjust_info_ptr, const APerm* aperm_ptr, const char* local_covar_fname, const char* local_pvar_fname, const char* local_psam_fname, const GwasSsfInfo* gsip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t orig_covar_ct, uintptr_t max_covar_name_blen, uint32_t raw_variant_ct, uint32_t orig_variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t xchr_model, double ci_size, double vif_thresh, double ln_pfilter, double output_min_ln, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;

  // We may temporarily clear the cip->haploid_mask bit for chrX if all samples
  // are female.  (This is only checked on function entry, it is not rechecked
  // on a per-phenotype basis.)  Track this here so we can reverse it on
  // function exit.
  uint32_t x_fully_diploid = 0;

  LlStr* gwas_ssf_ll = nullptr;
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
      reterr = GlmLocalOpen(local_covar_fname, local_pvar_fname, local_psam_fname, siip, cip, variant_bps, variant_ids, glm_info_ptr, raw_sample_ct, raw_variant_ct, &orig_sample_include, &sex_nm, &sex_male, &early_variant_include, &orig_sample_ct, &variant_ct, &local_covar_txs, &local_sample_uidx_order, &local_variant_include, &local_sample_ct, &local_variant_ctl, &local_covar_ct);
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
      logerrputs("Error: Phenotype name and/or --out argument too long.\n");
      goto GlmMain_ret_INCONSISTENT_INPUT;
    }
    *outname_end = '.';
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;

    // synthetic categorical covariate name could be ~twice max ID length?
    const uintptr_t overflow_buf_size = kCompressStreamBlock + 2 * kMaxIdSlen + max_chr_blen + kMaxIdSlen + 1024 + 2 * max_allele_slen;

    uintptr_t* cur_sample_include;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &cur_sample_include) ||
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
    uint32_t variant_ct_y = 0;
    const uint32_t domdev_present = (glm_flags & (kfGlmGenotypic | kfGlmHethom))? 1 : 0;
    const uint32_t sex_nm_ct = PopcountWords(sex_nm, raw_sample_ctl);
    const uint32_t male_ct = PopcountWords(sex_male, raw_sample_ctl);
    uint32_t add_sex_covar = !(glm_flags & kfGlmNoXSex);
    if (add_sex_covar && ((!male_ct) || (male_ct == sex_nm_ct))) {
      add_sex_covar = 0;
    }
    uint32_t variant_ct_x = 0;
    {
      uint32_t x_code;
      if (XymtExists(cip, kChrOffsetX, &x_code)) {
        variant_ct_x = CountChrVariantsUnsafe(early_variant_include, cip, x_code);
        x_fully_diploid = (!male_ct) && (sex_nm_ct == orig_sample_ct) && variant_ct_x && xchr_model;
        if (x_fully_diploid) {
          ClearBit(x_code, cip->haploid_mask);
        }
      }
    }
    uintptr_t* cur_sample_include_y_buf = nullptr;
    if (domdev_present || (glm_flags & (kfGlmDominant | kfGlmRecessive | kfGlmHetonly))) {
      // dominant/recessive/hetonly/genotypic/hethom suppress all chromosomes
      // which aren't fully diploid.

      xchr_model = 0;
      // update (18 Sep 2021): chrX is no longer suppressed if all samples are
      // female.
      if (x_fully_diploid) {
        xchr_model = 2;
        logputs("--glm: Including chrX, despite presence of a diploid-only modifier\n('dominant', 'recessive', 'hetonly', 'genotypic', 'hethom'), since all samples\nare female.\n");
      } else {
        // bugfix (20 Sep 2021): need this to avoid double-subtraction.
        variant_ct_x = 0;
      }
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
          logerrputs("Error: No variants remaining for --glm ('dominant', 'recessive', 'hetonly',\n'genotypic', and 'hethom' only operate on diploid data).\n");
          goto GlmMain_ret_DEGENERATE_DATA;
        }
        variant_ct -= removed_variant_ct;
        early_variant_include = variant_include_nohap;
        max_variant_ct = variant_ct;
      }
    } else {
      if (variant_ct_x) {
        // --xchr-model 0 now only suppresses chrX.
        if (xchr_model) {
          if (unlikely(bigstack_alloc_w(BitCtToWordCt(orig_sample_ct), &sex_male_collapsed_buf))) {
            goto GlmMain_ret_NOMEM;
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
    if (unlikely((!raw_covar_ct) && (!(glm_flags & kfGlmAllowNoCovars)))) {
      // now possible due to --not-covar
      logerrputs("Error: --glm invoked with no covariates, and 'allow-no-covars' was not\nspecified.\n");
      goto GlmMain_ret_INCONSISTENT_INPUT;
    }
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
            // TODO: this count will be different with a multiallelic variant
            logputs("--glm: One --condition covariate added.\n");
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
        if (unlikely(BIGSTACK_ALLOC_X(PhenoCol, raw_covar_ct + add_sex_covar, &new_covar_cols) ||
                     bigstack_alloc_c((raw_covar_ct + add_sex_covar) * new_max_covar_name_blen, &new_covar_names))) {
          goto GlmMain_ret_NOMEM;
        }
        if (condition_ct) {
          BigstackEndSet(condition_uidxs);
          uintptr_t* genovec;
          uintptr_t* dosage_present;
          Dosage* dosage_main;
          if (unlikely(bigstack_end_alloc_w(NypCtToWordCt(raw_sample_ct), &genovec) ||
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
              PgenErrPrintV(reterr, cur_variant_uidx);
              goto GlmMain_ret_1;
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
            if (unlikely(bigstack_alloc_w(raw_sample_ctl, &cur_nonmiss) ||
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
                    // this is technically allowed when all samples are female,
                    // but unimportant to mention that in the error message.
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
        if (unlikely(BIGSTACK_ALLOC_X(PhenoCol, raw_covar_ct + add_sex_covar, &new_covar_cols) ||
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
      if (unlikely(bigstack_alloc_w(raw_covar_ctl, &initial_covar_include) ||
                   bigstack_alloc_w(raw_covar_ctl, &covar_include))) {
        goto GlmMain_ret_NOMEM;
      }
      ZeroWArr(raw_covar_ctl, initial_covar_include);
      for (uint32_t covar_uidx = 0; covar_uidx != raw_covar_ct; ++covar_uidx) {
        const PhenoCol* cur_covar_col = &(covar_cols[covar_uidx]);
        if (cur_covar_col->type_code != kPhenoDtypeOther) {
          if (!IsConstCovar(cur_covar_col, orig_sample_include, raw_sample_ct)) {
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
      if (unlikely(bigstack_calloc_w(biallelic_raw_predictor_ctl, &raw_parameter_subset) ||
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
      if (unlikely(bigstack_calloc_w(biallelic_raw_predictor_ctl, &raw_joint_test_params) ||
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
        if (unlikely(bigstack_alloc_w(biallelic_raw_predictor_ctl, &common.joint_test_params_y) ||
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
      if (unlikely((!(raw_parameter_subset[0] & 2)) &&
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
          if (unlikely(bigstack_alloc_w(raw_sample_ctl, &cur_sample_include_x_buf) ||
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
    const uint32_t skip_invalid_pheno = (glm_flags / kfGlmSkipInvalidPheno) & 1;
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

    LlStr** gwas_ssf_ll_ptr = IsGwasSsf(gsip)? (&gwas_ssf_ll) : nullptr;
    uintptr_t* valid_variants = nullptr;
    uintptr_t* valid_alleles = nullptr;
    unsigned char* bigstack_mark2 = g_bigstack_base;
    if (report_adjust || perms_total) {
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, &valid_variants) ||
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
      if (unlikely(bigstack_alloc_w(pheno_ctl, &pheno_batch) ||
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
        if (IsConstCovar(cur_pheno_col, cur_sample_include, raw_sample_ct)) {
          const char* cur_pheno_name = &(pheno_names[pheno_uidx * max_pheno_name_blen]);
          if (unlikely(!skip_invalid_pheno)) {
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
          if (unlikely(!skip_invalid_pheno)) {
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
          if (unlikely(!skip_invalid_pheno)) {
            logerrprintfww("Error: Too many samples or predictors for --glm regression on phenotype '%s' (internal matrices currently limited to ~2^31 entries).\n", first_pheno_name);
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logerrprintfww("Warning: Skipping --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since there are too many samples or predictors (internal matrices currently limited to ~2^31 entries).\n", first_pheno_name);
          BitvecInvmask(pheno_batch, pheno_ctl, pheno_include);
          continue;
        }
#endif
        if (common.tests_flag && (!common.constraint_ct)) {
          if (unlikely(!skip_invalid_pheno)) {
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
              if (unlikely(!skip_invalid_pheno)) {
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
                if (unlikely(!skip_invalid_pheno)) {
                  logerrprintfww("Error: # samples <= # predictor columns for --glm phenotype '%s' on chrX.\n", first_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since # remaining samples <= # predictor columns.\n", first_pheno_name);
                sample_ct_x = 0;
#ifdef __LP64__
              } else if (RoundUpPow2(sample_ct_x, 4) * S_CAST(uint64_t, biallelic_predictor_ct_x + max_extra_allele_ct) > 0x7fffffff) {
                if (unlikely(!skip_invalid_pheno)) {
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
                  if (IsConstCovar(candidate_pheno_col, cur_sample_include_x, raw_sample_ct)) {
                    // Punt to single-phenotype-at-a-time handler.
                    ClearBit(pheno_uidx2, pheno_batch);
                  }
                }
                // possible todo: reset first_pheno_name here
              }
              if (sample_ct_x && common.tests_flag && (!common.constraint_ct_x)) {
                if (unlikely(!skip_invalid_pheno)) {
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
              if (unlikely(!skip_invalid_pheno)) {
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
                if (unlikely(!skip_invalid_pheno)) {
                  logerrprintfww("Error: # samples <= # predictor columns for --glm phenotype '%s' on chrY.\n", first_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', and other(s) with identical missingness patterns, since # remaining samples <= # predictor columns.\n", first_pheno_name);
                sample_ct_y = 0;
#ifdef __LP64__
              } else if (RoundUpPow2(sample_ct_y, 4) * S_CAST(uint64_t, biallelic_predictor_ct_y + max_extra_allele_ct) > 0x7fffffff) {
                if (unlikely(!skip_invalid_pheno)) {
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
                  if (IsConstCovar(candidate_pheno_col, cur_sample_include_y, raw_sample_ct)) {
                    // Punt to single-phenotype-at-a-time handler.
                    ClearBit(pheno_uidx2, pheno_batch);
                  }
                }
              }
              if (sample_ct_y && common.tests_flag && (!common.constraint_ct_y)) {
                if (unlikely(!skip_invalid_pheno)) {
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
        const char** cur_covar_names = nullptr;
        GlmErr glm_err;
        {
          double* covars_cmaj_d = nullptr;
          if (unlikely(GlmAllocFillAndTestCovarsQt(cur_sample_include, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &common.nm_precomp, &covars_cmaj_d, &cur_covar_names, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
          linear_ctx.covars_cmaj_d = covars_cmaj_d;
          if (glm_err) {
            PrintPrescanErrmsg("", first_pheno_name, cur_covar_names, glm_err, local_covar_ct, skip_invalid_pheno, 1);
            if (unlikely(!skip_invalid_pheno)) {
              goto GlmMain_ret_INCONSISTENT_INPUT;
            }
            BitvecInvmask(pheno_batch, pheno_ctl, pheno_include);
            continue;
          }
        }
        const char** cur_covar_names_x = nullptr;
        common.nm_precomp_x = nullptr;
        linear_ctx.covars_cmaj_x_d = nullptr;
        if (sample_ct_x) {
          double* covars_cmaj_d = nullptr;
          if (unlikely(GlmAllocFillAndTestCovarsQt(cur_sample_include_x, covar_include_x, covar_cols, covar_names, sample_ct_x, covar_ct_x, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_x, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &common.nm_precomp_x, &covars_cmaj_d, &cur_covar_names_x, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
          linear_ctx.covars_cmaj_x_d = covars_cmaj_d;
          if (glm_err) {
            PrintPrescanErrmsg("chrX in ", first_pheno_name, cur_covar_names_x, glm_err, local_covar_ct, skip_invalid_pheno, 1);
            if (unlikely(!skip_invalid_pheno)) {
              goto GlmMain_ret_INCONSISTENT_INPUT;
            }
            sample_ct_x = 0;
          }
        }
        const char** cur_covar_names_y = nullptr;
        common.nm_precomp_y = nullptr;
        linear_ctx.covars_cmaj_y_d = nullptr;
        if (sample_ct_y) {
          double* covars_cmaj_d = nullptr;
          if (unlikely(GlmAllocFillAndTestCovarsQt(cur_sample_include_y, covar_include_y, covar_cols, covar_names, sample_ct_y, covar_ct_y, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_y, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &common.nm_precomp_y, &covars_cmaj_d, &cur_covar_names_y, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
          linear_ctx.covars_cmaj_y_d = covars_cmaj_d;
          if (glm_err) {
            PrintPrescanErrmsg("chrY in ", first_pheno_name, cur_covar_names_y, glm_err, local_covar_ct, skip_invalid_pheno, 1);
            if (unlikely(!skip_invalid_pheno)) {
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
          uint32_t* cumulative_popcounts;
          if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &cumulative_popcounts))) {
            goto GlmMain_ret_NOMEM;
          }
          FillCumulativePopcounts(cur_sample_include_x, raw_sample_ctl, cumulative_popcounts);
          common.sample_include_x_cumulative_popcounts = cumulative_popcounts;
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
          uint32_t* cumulative_popcounts;
          if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &cumulative_popcounts))) {
            goto GlmMain_ret_NOMEM;
          }
          FillCumulativePopcounts(cur_sample_include_y, raw_sample_ctl, cumulative_popcounts);
          common.sample_include_y_cumulative_popcounts = cumulative_popcounts;
          common.sample_include_y = cur_sample_include_y;
          common.covar_ct_y = covar_ct_y + extra_cat_ct_y;
        } else {
          common.sample_include_y = nullptr;
          common.sample_include_y_cumulative_popcounts = nullptr;
          common.covar_ct_y = 0;
        }

        uint32_t* subset_chr_fo_vidx_start;
        if (unlikely(AllocAndFillSubsetChrFoVidxStart(cur_variant_include, cip, &subset_chr_fo_vidx_start))) {
          goto GlmMain_ret_NOMEM;
        }
        common.subset_chr_fo_vidx_start = subset_chr_fo_vidx_start;
        common.variant_include = cur_variant_include;
        common.variant_ct = cur_variant_ct;

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

        reterr = GlmLinearBatch(pheno_batch, pheno_cols, pheno_names, cur_test_names, cur_test_names_x, cur_test_names_y, glm_pos_col? variant_bps : nullptr, variant_ids, allele_storage, glm_info_ptr, local_sample_uidx_order, cur_local_variant_include, raw_variant_ct, completed_pheno_ct, batch_size, max_pheno_name_blen, max_chr_blen, ci_size, ln_pfilter, output_min_ln, max_thread_ct, pgr_alloc_cacheline_ct, overflow_buf_size, local_sample_ct, pgfip, &linear_ctx, &local_covar_txs, gwas_ssf_ll_ptr, outname, outname_end);
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
          if (unlikely(!skip_invalid_pheno)) {
            logerrprintfww("Error: All samples for --glm phenotype '%s' are %s.\n", cur_pheno_name, initial_case_ct? "cases" : "controls");
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logprintfww("--glm: Skipping case/control phenotype '%s' since all samples are %s.\n", cur_pheno_name, initial_case_ct? "cases" : "controls");
          continue;
        }
      } else {
        if (IsConstCovar(cur_pheno_col, cur_sample_include, raw_sample_ct)) {
          if (unlikely(!skip_invalid_pheno)) {
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
        if (unlikely(!skip_invalid_pheno)) {
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
        if (unlikely(!skip_invalid_pheno)) {
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
          if (unlikely(!skip_invalid_pheno)) {
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
        if (IsConstCovar(cur_pheno_col, cur_sample_include, raw_sample_ct)) {
          if (unlikely(!skip_invalid_pheno)) {
            logerrprintfww("Error: --glm quantitative phenotype '%s' is constant for all remaining samples.\n", cur_pheno_name);
            goto GlmMain_ret_INCONSISTENT_INPUT;
          }
          logprintfww("--glm: Skipping quantitative phenotype '%s' since phenotype is constant for all remaining samples.\n", cur_pheno_name);
          continue;
        }
      }
      if (common.tests_flag && (!common.constraint_ct)) {
        if (unlikely(!skip_invalid_pheno)) {
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
            if (unlikely(!skip_invalid_pheno)) {
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
                // bugfix (14 Sep 2020): forgot to use parameter_subset_x
                // instead of parameter_subset here
                CopyBitarrSubset(joint_test_params_buf, common.parameter_subset_x, biallelic_predictor_ct_x, common.joint_test_params_x);
              }
            }
            if (sample_ct_x <= biallelic_predictor_ct_x) {
              if (unlikely(!skip_invalid_pheno)) {
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
              if (unlikely(!skip_invalid_pheno)) {
                logerrprintfww("Error: Too many samples or predictors for --glm regression on phenotype '%s' on chrX (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              logerrprintfww("Warning: Skipping chrX in --glm regression on phenotype '%s', since there are too many samples or predictors (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
              sample_ct_x = 0;
#endif
            } else if (is_logistic) {
              const uint32_t case_ct_x = PopcountWordsIntersect(cur_sample_include_x, cur_pheno_col->data.cc, raw_sample_ctl);
              if ((!case_ct_x) || (case_ct_x == sample_ct_x)) {
                if (unlikely(!skip_invalid_pheno)) {
                  logerrprintfww("Error: All remaining samples on chrX for --glm phenotype '%s' are %s.\n", cur_pheno_name, case_ct? "cases" : "controls");
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', since all remaining samples are %s.\n", cur_pheno_name, case_ct_x? "cases" : "controls");
                sample_ct_x = 0;
              }
            } else {
              if (IsConstCovar(cur_pheno_col, cur_sample_include_x, raw_sample_ct)) {
                if (unlikely(!skip_invalid_pheno)) {
                  logerrprintfww("Error: --glm quantitative phenotype '%s' is constant for all remaining samples on chrX.\n", cur_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrX in --glm regression on phenotype '%s', since phenotype is constant for all remaining samples.\n", cur_pheno_name);
                sample_ct_x = 0;
              }
            }
            if (sample_ct_x && common.tests_flag && (!common.constraint_ct_x)) {
              if (unlikely(!skip_invalid_pheno)) {
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
            if (unlikely(!skip_invalid_pheno)) {
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
                CopyBitarrSubset(joint_test_params_buf, common.parameter_subset_y, biallelic_predictor_ct_y, common.joint_test_params_y);
              }
            }
            if (sample_ct_y <= biallelic_predictor_ct_y) {
              if (unlikely(!skip_invalid_pheno)) {
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
              if (unlikely(!skip_invalid_pheno)) {
                logerrprintfww("Error: Too many samples or predictors for --glm regression on phenotype '%s' on chrY (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
                goto GlmMain_ret_INCONSISTENT_INPUT;
              }
              logerrprintfww("Warning: Skipping chrY in --glm regression on phenotype '%s', since there are too many samples or predictors (internal matrices currently limited to ~2^31 entries).\n", cur_pheno_name);
              sample_ct_y = 0;
#endif
            } else if (is_logistic) {
              const uint32_t case_ct_y = PopcountWordsIntersect(cur_sample_include_y, cur_pheno_col->data.cc, raw_sample_ctl);
              if ((!case_ct_y) || (case_ct_y == sample_ct_y)) {
                if (unlikely(!skip_invalid_pheno)) {
                  logerrprintfww("Error: All remaining samples on chrY for --glm phenotype '%s' are %s.\n", cur_pheno_name, case_ct? "cases" : "controls");
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', since all remaining samples are %s.\n", cur_pheno_name, case_ct_y? "cases" : "controls");
                sample_ct_y = 0;
              }
            } else {
              if (IsConstCovar(cur_pheno_col, cur_sample_include_y, raw_sample_ct)) {
                if (unlikely(!skip_invalid_pheno)) {
                  logerrprintfww("Error: --glm quantitative phenotype '%s' is constant for all remaining samples on chrY.\n", cur_pheno_name);
                  goto GlmMain_ret_INCONSISTENT_INPUT;
                }
                logprintfww("Note: Skipping chrY in --glm regression on phenotype '%s', since phenotype is constant for all remaining samples.\n", cur_pheno_name);
                sample_ct_y = 0;
              }
            }
            if (sample_ct_y && common.tests_flag && (!common.constraint_ct_y)) {
              if (unlikely(!skip_invalid_pheno)) {
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
      const char** cur_covar_names = nullptr;
      GlmErr glm_err;
      logistic_ctx.pheno_cc = nullptr;
      logistic_ctx.gcount_case_interleaved_vec = nullptr;
      logistic_ctx.cc_residualize = nullptr;
      if (is_logistic) {
        linear_ctx.pheno_d = nullptr;
        linear_ctx.covars_cmaj_d = nullptr;
        float* pheno_f = nullptr;
        double* pheno_d = nullptr;
        float* covars_cmaj_f = nullptr;
        double* covars_cmaj_d = nullptr;
        if (unlikely(GlmAllocFillAndTestPhenoCovarsCc(cur_sample_include, cur_pheno_col->data.cc, covar_include, covar_cols, covar_names, sample_ct, domdev_present_p1, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, glm_flags, &logistic_ctx.pheno_cc, gcount_cc_col? (&logistic_ctx.gcount_case_interleaved_vec) : nullptr, &pheno_f, &pheno_d, &common.nm_precomp, &covars_cmaj_f, &covars_cmaj_d, &logistic_ctx.cc_residualize, &cur_covar_names, &glm_err))) {
          goto GlmMain_ret_NOMEM;
        }
        logistic_ctx.pheno_f = pheno_f;
        logistic_ctx.pheno_d = pheno_d;
        logistic_ctx.covars_cmaj_f = covars_cmaj_f;
        logistic_ctx.covars_cmaj_d = covars_cmaj_d;
      } else {
        logistic_ctx.pheno_f = nullptr;
        logistic_ctx.covars_cmaj_f = nullptr;
        logistic_ctx.pheno_d = nullptr;
        logistic_ctx.covars_cmaj_d = nullptr;
        double* pheno_d = nullptr;
        double* covars_cmaj_d = nullptr;
        if (unlikely(GlmAllocFillAndTestPhenoCovarsQt(cur_sample_include, cur_pheno_col->data.qt, covar_include, covar_cols, covar_names, sample_ct, covar_ct, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &pheno_d, &common.nm_precomp, &covars_cmaj_d, &cur_covar_names, &glm_err))) {
          goto GlmMain_ret_NOMEM;
        }
        linear_ctx.pheno_d = pheno_d;
        linear_ctx.covars_cmaj_d = covars_cmaj_d;
      }
      if (glm_err) {
        PrintPrescanErrmsg("", cur_pheno_name, cur_covar_names, glm_err, local_covar_ct, skip_invalid_pheno, 0);
        if (unlikely(!skip_invalid_pheno)) {
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
      logistic_ctx.pheno_x_d = nullptr;
      logistic_ctx.covars_cmaj_x_d = nullptr;
      logistic_ctx.cc_residualize_x = nullptr;
      linear_ctx.pheno_x_d = nullptr;
      linear_ctx.covars_cmaj_x_d = nullptr;
      if (sample_ct_x) {
        if (is_logistic) {
          float* pheno_f = nullptr;
          double* pheno_d = nullptr;
          float* covars_cmaj_f = nullptr;
          double* covars_cmaj_d = nullptr;
          if (unlikely(GlmAllocFillAndTestPhenoCovarsCc(cur_sample_include_x, cur_pheno_col->data.cc, covar_include_x, covar_cols, covar_names, sample_ct_x, domdev_present_p1, covar_ct_x, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_x, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, glm_flags, &logistic_ctx.pheno_x_cc, gcount_cc_col? (&logistic_ctx.gcount_case_interleaved_vec_x) : nullptr, &pheno_f, &pheno_d, &common.nm_precomp_x, &covars_cmaj_f, &covars_cmaj_d, &logistic_ctx.cc_residualize_x, &cur_covar_names_x, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
          logistic_ctx.pheno_x_f = pheno_f;
          logistic_ctx.pheno_x_d = pheno_d;
          logistic_ctx.covars_cmaj_x_f = covars_cmaj_f;
          logistic_ctx.covars_cmaj_x_d = covars_cmaj_d;
        } else {
          double* pheno_d = nullptr;
          double* covars_cmaj_d = nullptr;
          if (unlikely(GlmAllocFillAndTestPhenoCovarsQt(cur_sample_include_x, cur_pheno_col->data.qt, covar_include_x, covar_cols, covar_names, sample_ct_x, covar_ct_x, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_x, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &pheno_d, &common.nm_precomp_x, &covars_cmaj_d, &cur_covar_names_x, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
          linear_ctx.pheno_x_d = pheno_d;
          linear_ctx.covars_cmaj_x_d = covars_cmaj_d;
        }
        if (glm_err) {
          PrintPrescanErrmsg("chrX in ", cur_pheno_name, cur_covar_names_x, glm_err, local_covar_ct, skip_invalid_pheno, 0);
          if (unlikely(!skip_invalid_pheno)) {
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
      logistic_ctx.pheno_y_d = nullptr;
      logistic_ctx.covars_cmaj_y_d = nullptr;
      logistic_ctx.cc_residualize_y = nullptr;
      linear_ctx.pheno_y_d = nullptr;
      linear_ctx.covars_cmaj_y_d = nullptr;
      if (sample_ct_y) {
        if (is_logistic) {
          float* pheno_f = nullptr;
          double* pheno_d = nullptr;
          float* covars_cmaj_f = nullptr;
          double* covars_cmaj_d = nullptr;
          if (unlikely(GlmAllocFillAndTestPhenoCovarsCc(cur_sample_include_y, cur_pheno_col->data.cc, covar_include_y, covar_cols, covar_names, sample_ct_y, domdev_present_p1, covar_ct_y, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_y, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, glm_flags, &logistic_ctx.pheno_y_cc, gcount_cc_col? (&logistic_ctx.gcount_case_interleaved_vec_y) : nullptr, &pheno_f, &pheno_d, &common.nm_precomp_y, &covars_cmaj_f, &covars_cmaj_d, &logistic_ctx.cc_residualize_y, &cur_covar_names_y, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
          logistic_ctx.pheno_y_f = pheno_f;
          logistic_ctx.pheno_y_d = pheno_d;
          logistic_ctx.covars_cmaj_y_f = covars_cmaj_f;
          logistic_ctx.covars_cmaj_y_d = covars_cmaj_d;
        } else {
          double* pheno_d = nullptr;
          double* covars_cmaj_d = nullptr;
          if (unlikely(GlmAllocFillAndTestPhenoCovarsQt(cur_sample_include_y, cur_pheno_col->data.qt, covar_include_y, covar_cols, covar_names, sample_ct_y, covar_ct_y, local_covar_ct, covar_max_nonnull_cat_ct, extra_cat_ct_y, max_covar_name_blen, common.max_corr, vif_thresh, xtx_state, &pheno_d, &common.nm_precomp_y, &covars_cmaj_d, &cur_covar_names_y, &glm_err))) {
            goto GlmMain_ret_NOMEM;
          }
          linear_ctx.pheno_y_d = pheno_d;
          linear_ctx.covars_cmaj_y_d = covars_cmaj_d;
        }
        if (glm_err) {
          PrintPrescanErrmsg("chrY in ", cur_pheno_name, cur_covar_names_y, glm_err, local_covar_ct, skip_invalid_pheno, 0);
          if (unlikely(!skip_invalid_pheno)) {
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
        uint32_t* cumulative_popcounts;
        if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &cumulative_popcounts))) {
          goto GlmMain_ret_NOMEM;
        }
        FillCumulativePopcounts(cur_sample_include_x, raw_sample_ctl, cumulative_popcounts);
        common.sample_include_x_cumulative_popcounts = cumulative_popcounts;
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
        uint32_t* cumulative_popcounts;
        if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &cumulative_popcounts))) {
          goto GlmMain_ret_NOMEM;
        }
        FillCumulativePopcounts(cur_sample_include_y, raw_sample_ctl, cumulative_popcounts);
        common.sample_include_y_cumulative_popcounts = cumulative_popcounts;
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

      uint32_t* subset_chr_fo_vidx_start;
      if (unlikely(AllocAndFillSubsetChrFoVidxStart(cur_variant_include, cip, &subset_chr_fo_vidx_start))) {
        goto GlmMain_ret_NOMEM;
      }
      common.subset_chr_fo_vidx_start = subset_chr_fo_vidx_start;
      common.variant_include = cur_variant_include;
      common.variant_ct = cur_variant_ct;
      // this is safe, see pheno_name_blen_capacity check above
      char* outname_end2 = strcpya(&(outname_end[1]), cur_pheno_name);
      if (is_logistic) {
        if (is_always_firth) {
          outname_end2 = strcpya_k(outname_end2, ".glm.firth");
        } else if (is_sometimes_firth) {
          outname_end2 = strcpya_k(outname_end2, ".glm.logistic.hybrid");
        } else {
          outname_end2 = strcpya_k(outname_end2, ".glm.logistic");
        }
      } else {
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
        reterr = GlmLogistic(cur_pheno_name, cur_test_names, cur_test_names_x, cur_test_names_y, glm_pos_col? variant_bps : nullptr, variant_ids, allele_storage, glm_info_ptr, local_sample_uidx_order, cur_local_variant_include, outname, raw_variant_ct, max_chr_blen, ci_size, ln_pfilter, output_min_ln, max_thread_ct, pgr_alloc_cacheline_ct, overflow_buf_size, local_sample_ct, pgfip, &logistic_ctx, &local_covar_txs, gwas_ssf_ll_ptr, valid_variants, valid_alleles, orig_ln_pvals, orig_permstat, &valid_allele_ct);
      } else {
        reterr = GlmLinear(cur_pheno_name, cur_test_names, cur_test_names_x, cur_test_names_y, glm_pos_col? variant_bps : nullptr, variant_ids, allele_storage, glm_info_ptr, local_sample_uidx_order, cur_local_variant_include, outname, raw_variant_ct, max_chr_blen, ci_size, ln_pfilter, output_min_ln, max_thread_ct, pgr_alloc_cacheline_ct, overflow_buf_size, local_sample_ct, pgfip, &linear_ctx, &local_covar_txs, gwas_ssf_ll_ptr, valid_variants, valid_alleles, orig_ln_pvals, &valid_allele_ct);
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
    if (gwas_ssf_ll_ptr) {
      const uint32_t delete_orig_glm = gsip->flags & kfGwasSsfDeleteOrigGlm;
      while (gwas_ssf_ll) {
        const char* in_fname = gwas_ssf_ll->str;
        reterr = GwasSsfOneFile(gsip, in_fname, max_thread_ct);
        if (unlikely(reterr)) {
          goto GlmMain_ret_1;
        }
        // Don't need to free linked list nodes yet, but we may need to delete
        // files here.
        if (delete_orig_glm) {
          if (unlikely(unlink(in_fname))) {
            logerrprintfww("Error: --gwas-ssf delete-orig-glm: Failed to delete %s .\n", in_fname);
            goto GlmMain_ret_WRITE_FAIL;
          }
        }
        gwas_ssf_ll = gwas_ssf_ll->next;
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
  GlmMain_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
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
  llstr_free_cond(gwas_ssf_ll);
  CleanupTokenStream2("--condition-list file", &tks, &reterr);
  CleanupTextStream2(local_covar_fname, &local_covar_txs, &reterr);
  if (x_fully_diploid) {
    SetBit(cip->xymt_codes[kChrOffsetX], cip->haploid_mask);
  }
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

/*
void LogisticTest() {
  LogisticTestInternal();
}
*/

#ifdef __cplusplus
}  // namespace plink2
#endif
