// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
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

#include "plink2_psam.h"

#ifdef __cplusplus
namespace plink2 {
#endif

typedef struct PsamInfoLlStruct {
  NONCOPYABLE(PsamInfoLlStruct);
  // vardata[] starts with 8-byte phenotypes (we don't want to parse the same
  // same numeric string twice), followed by NON-null-terminated sample_id, and
  // then non-terminated paternal and maternal IDs.
  struct PsamInfoLlStruct* next;
  uint32_t sample_id_slen;
  uint32_t sid_slen;
  uint32_t paternal_id_slen;
  uint32_t maternal_id_slen;
  uint32_t sex_code;  // 0 = unknown, 1 = male, 2 = female
  unsigned char vardata[];
} PsamInfoLl;

typedef struct CatnameLl2Struct {
  NONCOPYABLE(CatnameLl2Struct);
  struct CatnameLl2Struct* htable_next;
  struct CatnameLl2Struct* pheno_next;
  uint32_t cat_idx;  // 0 == "NONE", etc.
  char str[];
} CatnameLl2;

PglErr LoadPsam(const char* psamname, const RangeList* pheno_range_list_ptr, FamCol fam_cols, uint32_t pheno_ct_max, int32_t missing_pheno, uint32_t affection_01, uint32_t max_thread_ct, PedigreeIdInfo* piip, uintptr_t** sample_include_ptr, uintptr_t** founder_info_ptr, uintptr_t** sex_nm_ptr, uintptr_t** sex_male_ptr, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* raw_sample_ct_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr) {
  // outparameter pointers assumed to be initialized to nullptr
  //
  // pheno_ct_max should default to something like 0x7fffffff, not UINT32_MAX
  //
  // FidPresent flag and max_{sample,sid,paternal,maternal}_id_blen are in/out,
  // for interoperation with --update-ids and --update-parents.
  //
  // permanent allocations are at stack end, not base, to work better with
  // VariantIdDupflagHtableFind()

  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;

  PhenoCol* pheno_cols = nullptr;
  uintptr_t line_idx = 0;
  uint32_t pheno_ct = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream psam_txs;
  PreinitTextStream(&psam_txs);
  {
    reterr = SizeAndInitTextStream(psamname, bigstack_left() / 4, MAXV(max_thread_ct - 1, 1), &psam_txs);
    if (unlikely(reterr)) {
      goto LoadPsam_ret_TSTREAM_FAIL;
    }
    const char* line_iter;
    do {
      ++line_idx;
      line_iter = TextGet(&psam_txs);
      if (unlikely(!line_iter)) {
        if (!TextStreamErrcode2(&psam_txs, &reterr)) {
          logerrprintfww("Error: No samples in %s.\n", psamname);
          goto LoadPsam_ret_MALFORMED_INPUT;
        }
        goto LoadPsam_ret_TSTREAM_FAIL;
      }
    } while ((line_iter[0] == '#') && (!tokequal_k(&(line_iter[1]), "FID")) && (!tokequal_k(&(line_iter[1]), "IID")));
    const uint32_t pheno_name_subset = pheno_range_list_ptr && pheno_range_list_ptr->names;
    uint32_t* col_skips = nullptr;
    uint32_t* col_types = nullptr;
    uint32_t psam_cols_mask = 0;
    LlStr* pheno_names_reverse_ll = nullptr;
    uintptr_t max_pheno_name_blen = 0;
    uint32_t relevant_postfid_col_ct = 0;
    g_bigstack_end -= kMaxIdSlen;
    unsigned char* tmp_bigstack_end = g_bigstack_end;
    unsigned char* bigstack_mark2;
    uint32_t fid_present = 1;
    if (line_iter[0] == '#') {
      // parse header
      // [-1] = #FID (if present, must be first column)
      // [0] = IID (could also be first column)
      // [1] = SID
      // [2] = PAT
      // [3] = MAT
      // [4] = SEX
      // [5+] = phenotypes
      relevant_postfid_col_ct = CountTokens(line_iter);
      if (relevant_postfid_col_ct > pheno_ct_max + 5) {
        relevant_postfid_col_ct = pheno_ct_max + 5;
      }
      if (unlikely(
              bigstack_alloc_u32(relevant_postfid_col_ct, &col_skips) ||
              bigstack_alloc_u32(relevant_postfid_col_ct, &col_types))) {
        goto LoadPsam_ret_NOMEM;
      }
      bigstack_mark2 = g_bigstack_base;
      uint32_t rpf_col_idx = 0;
      if (line_iter[1] == 'I') {
        col_skips[0] = 0;
        col_types[0] = 0;
        ++rpf_col_idx;
        psam_cols_mask = 1;
        fid_present = 0;
      }
      uint32_t in_interval = 0;
      char* cmdline_pheno_sorted_ids = nullptr;
      uint32_t* cmdline_pheno_id_map = nullptr;
      uintptr_t max_cmdline_pheno_id_blen = 0;
      uintptr_t cmdline_pheno_name_ct = 0;
      if (pheno_name_subset) {
        max_cmdline_pheno_id_blen = pheno_range_list_ptr->name_max_blen;
        cmdline_pheno_name_ct = pheno_range_list_ptr->name_ct;
        uintptr_t* dummy_bitarr;
        // don't bother freeing these before LoadPsam() is done
        if (unlikely(
                bigstack_alloc_c(cmdline_pheno_name_ct * max_cmdline_pheno_id_blen, &cmdline_pheno_sorted_ids) ||
                bigstack_alloc_u32(cmdline_pheno_name_ct, &cmdline_pheno_id_map) ||
                bigstack_alloc_w(BitCtToWordCt(cmdline_pheno_name_ct), &dummy_bitarr))) {
          goto LoadPsam_ret_NOMEM;
        }
        SetAllBits(cmdline_pheno_name_ct, dummy_bitarr);
        reterr = CopySortStrboxSubsetNoalloc(dummy_bitarr, pheno_range_list_ptr->names, cmdline_pheno_name_ct, max_cmdline_pheno_id_blen, 0, 0, 0, cmdline_pheno_sorted_ids, cmdline_pheno_id_map);
        if (unlikely(reterr)) {
          goto LoadPsam_ret_1;
        }
        BigstackReset(dummy_bitarr);
      }
      const char* token_end = &(line_iter[4]);
      unsigned char* ll_alloc_base = g_bigstack_base;
      const char* linebuf_iter;
      for (uint32_t col_idx = 1; ; ++col_idx) {
        linebuf_iter = FirstNonTspace(token_end);
        if (IsEolnKns(*linebuf_iter)) {
          break;
        }
        token_end = CurTokenEnd(linebuf_iter);
        const uint32_t token_slen = token_end - linebuf_iter;
        if (token_slen == 3) {
          uint32_t cur_col_type = UINT32_MAX;
          if (memequal_k(linebuf_iter, "IID", 3)) {
            cur_col_type = 0;
          } else if (memequal_k(linebuf_iter, "SID", 3)) {
            cur_col_type = 1;
          } else if (memequal_k(linebuf_iter, "PAT", 3)) {
            cur_col_type = 2;
          } else if (memequal_k(linebuf_iter, "MAT", 3)) {
            cur_col_type = 3;
          } else if (memequal_k(linebuf_iter, "SEX", 3)) {
            cur_col_type = 4;
          } else if (unlikely(memequal_k(linebuf_iter, "FID", 3))) {
            snprintf(g_logbuf, kLogbufSize, "Error: 'FID' column header on line %" PRIuPTR " of %s is not at the beginning.\n", line_idx, psamname);
            goto LoadPsam_ret_MALFORMED_INPUT_WW;
          }
          if (cur_col_type != UINT32_MAX) {
            const uint32_t cur_col_type_shifted = 1 << cur_col_type;
            if (unlikely(psam_cols_mask & cur_col_type_shifted)) {
              // known token, so no overflow danger
              char* write_iter = strcpya_k(g_logbuf, "Error: Duplicate column header '");
              write_iter = memcpya_k(write_iter, linebuf_iter, 4);
              write_iter = strcpya_k(write_iter, "' on line ");
              write_iter = wtoa(line_idx, write_iter);
              write_iter = strcpya_k(write_iter, " of ");
              write_iter = strcpya(write_iter, psamname);
              memcpy_k(write_iter, ".\n\0", 4);
              goto LoadPsam_ret_MALFORMED_INPUT_WW;
            }
            psam_cols_mask |= cur_col_type_shifted;
            col_skips[rpf_col_idx] = col_idx;
            col_types[rpf_col_idx++] = cur_col_type;
            continue;
          }
        }
        if (pheno_ct < pheno_ct_max) {
          if (pheno_name_subset) {
            uint32_t cmdline_pos;
            if (!SortedIdboxFind(linebuf_iter, cmdline_pheno_sorted_ids, cmdline_pheno_id_map, token_slen, max_cmdline_pheno_id_blen, cmdline_pheno_name_ct, &cmdline_pos)) {
              // similar to string_range_list_to_bitarr()
              if (pheno_range_list_ptr->starts_range[cmdline_pos]) {
                if (unlikely(in_interval)) {
                  logerrputs("Error: Overlapping --pheno-name ranges.\n");
                  goto LoadPsam_ret_INCONSISTENT_INPUT;
                }
                in_interval = 1;
              } else if (cmdline_pos && pheno_range_list_ptr->starts_range[cmdline_pos - 1]) {
                if (unlikely(!in_interval)) {
                  snprintf(g_logbuf, kLogbufSize, "Error: --pheno-name range is inconsistent with %s.\n", psamname);
                  goto LoadPsam_ret_INCONSISTENT_INPUT_WW;
                }
                in_interval = 0;
              }
            } else if (!in_interval) {
              continue;
            }
          }
          const uint32_t tok_blen = token_slen + 1;
          LlStr* ll_str_new = R_CAST(LlStr*, ll_alloc_base);
          // just word-aligned, not cacheline-aligned
          ll_alloc_base += RoundUpPow2(tok_blen + sizeof(LlStr), kBytesPerWord);
          if (unlikely(ll_alloc_base > tmp_bigstack_end)) {
            goto LoadPsam_ret_NOMEM;
          }
          ll_str_new->next = pheno_names_reverse_ll;
          memcpyx(ll_str_new->str, linebuf_iter, token_slen, '\0');
          if (tok_blen > max_pheno_name_blen) {
            max_pheno_name_blen = tok_blen;
          }
          pheno_names_reverse_ll = ll_str_new;
          col_skips[rpf_col_idx] = col_idx;
          col_types[rpf_col_idx++] = pheno_ct + 5;
          ++pheno_ct;
        }
      }
      if (unlikely(max_pheno_name_blen > kMaxIdBlen)) {
        logerrputs("Error: Phenotype/covariate names are limited to " MAX_ID_SLEN_STR " characters.\n");
        goto LoadPsam_ret_MALFORMED_INPUT;
      }
      BigstackBaseSet(ll_alloc_base);
      if (unlikely(!(psam_cols_mask & 1))) {
        snprintf(g_logbuf, kLogbufSize, "Error: No IID column on line %" PRIuPTR " of %s.\n", line_idx, psamname);
        goto LoadPsam_ret_MALFORMED_INPUT_WW;
      }
      if (unlikely(in_interval)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --pheno-name range is inconsistent with %s.\n", psamname);
        goto LoadPsam_ret_INCONSISTENT_INPUT_WW;
      }
      relevant_postfid_col_ct = rpf_col_idx;
      for (rpf_col_idx = relevant_postfid_col_ct - 1; rpf_col_idx; --rpf_col_idx) {
        col_skips[rpf_col_idx] -= col_skips[rpf_col_idx - 1];
      }

      line_iter = AdvPastDelim(linebuf_iter, '\n');
      ++line_idx;
    } else {
      if (unlikely(pheno_name_subset)) {
        logerrputs("Error: --pheno-name requires a --pheno or .psam file with a header.\n");
        goto LoadPsam_ret_INCONSISTENT_INPUT;
      }

      pheno_ct = (fam_cols & kfFamCol6) && pheno_ct_max;
      fid_present = (fam_cols / kfFamCol1) & 1;
      relevant_postfid_col_ct = fid_present + ((fam_cols / (kfFamCol34 / 2)) & 2) + ((fam_cols / kfFamCol5) & 1) + pheno_ct;
      // these small allocations can't fail, since kMaxMediumLine <
      // linebuf_size <= 1/3 of remaining space
      col_skips = S_CAST(uint32_t*, bigstack_alloc_raw_rd(relevant_postfid_col_ct * sizeof(int32_t)));
      col_types = S_CAST(uint32_t*, bigstack_alloc_raw_rd(relevant_postfid_col_ct * sizeof(int32_t)));
      bigstack_mark2 = g_bigstack_base;
      col_skips[0] = fam_cols & 1;  // assumes kfFamCol1 == 1
      col_types[0] = 0;
      // psam_cols_mask = 1;  // may need this later
      uint32_t rpf_col_idx = 1;
      if (fam_cols & kfFamCol34) {
        col_skips[rpf_col_idx] = 1;
        col_types[rpf_col_idx++] = 2;
        col_skips[rpf_col_idx] = 1;
        col_types[rpf_col_idx++] = 3;
        psam_cols_mask |= 12;
      }
      if (fam_cols & kfFamCol5) {
        col_skips[rpf_col_idx] = 1;
        col_types[rpf_col_idx++] = 4;
        psam_cols_mask |= 0x10;
      }
      if (pheno_ct) {
        col_skips[rpf_col_idx] = 1;
        // col_types[rpf_col_idx++] = 6;
        col_types[rpf_col_idx] = 5;
        LlStr* ll_str_new = S_CAST(LlStr*, bigstack_alloc_raw_rd(7 + sizeof(LlStr)));
        ll_str_new->next = pheno_names_reverse_ll;
        strcpy_k(ll_str_new->str, "PHENO1");
        max_pheno_name_blen = 7;
        pheno_names_reverse_ll = ll_str_new;
      }
    }
    if (psam_cols_mask & 12) {
      piip->sii.flags |= kfSampleIdParentsPresent;
    }
    if (pheno_ct) {
      char* pheno_names;
      if (unlikely(pgl_malloc(pheno_ct * max_pheno_name_blen, &pheno_names))) {
        goto LoadPsam_ret_NOMEM;
      }
      *pheno_names_ptr = pheno_names;
      for (uint32_t pheno_idx = pheno_ct; pheno_idx; ) {
        --pheno_idx;
        strcpy(&(pheno_names[pheno_idx * max_pheno_name_blen]), pheno_names_reverse_ll->str);
        pheno_names_reverse_ll = pheno_names_reverse_ll->next;
      }
      if (pheno_ct > 1) {
        if (unlikely(pheno_ct > kMaxPhenoCt)) {
          // yeah, yeah, this will never come up
          logerrputs("Error: " PROG_NAME_STR " does not support more than " MAX_PHENO_CT_STR " phenotypes.\n");
          goto LoadPsam_ret_MALFORMED_INPUT;
        }
        // verify there are no duplicates
        uint32_t tmp_htable_size;
        uint32_t* htable_tmp;
        if (unlikely(HtableGoodSizeAlloc(pheno_ct, bigstack_left(), &htable_tmp, &tmp_htable_size))) {
          goto LoadPsam_ret_NOMEM;
        }
        const uint32_t duplicate_idx = PopulateStrboxHtable(pheno_names, pheno_ct, max_pheno_name_blen, tmp_htable_size, htable_tmp);
        if (unlikely(duplicate_idx)) {
          const char* duplicate_pheno_name = &(pheno_names[duplicate_idx * max_pheno_name_blen]);
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate phenotype/covariate name '%s' on line %" PRIuPTR " of %s.\n", duplicate_pheno_name, line_idx, psamname);
          goto LoadPsam_ret_MALFORMED_INPUT_WW;
        }
      }
      // free pheno_names_reverse_ll
      BigstackReset(bigstack_mark2);
    }

    // make sure to error out properly in sample_ct == 0 case
    PsamInfoLl* psam_info_reverse_ll = nullptr;
    const uint32_t sids_present = (psam_cols_mask / 2) & 1;
    const uint32_t paternal_ids_present = psam_cols_mask & 4;
    const uint32_t maternal_ids_present = psam_cols_mask & 8;
    const uint32_t sex_present = psam_cols_mask & 0x10;
    const uint32_t col_type_end = 5 + pheno_ct;
    const uint32_t pheno_ctl = BitCtToWordCt(pheno_ct);
    const double missing_phenod = missing_pheno? S_CAST(double, missing_pheno) : HUGE_VAL;
    const double pheno_ctrld = u31tod(1 - affection_01);
    const double pheno_cased = pheno_ctrld + 1.0;
    uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    uintptr_t max_sid_blen = piip->sii.max_sid_blen;
    uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
    uint32_t raw_sample_ct = 0;
    uint32_t categorical_pheno_ct = 0;

    const char** token_ptrs;
    uint32_t* token_slens;
    uintptr_t* categorical_phenos;
    uintptr_t* quantitative_phenos;
    if (unlikely(
            bigstack_alloc_kcp(col_type_end, &token_ptrs) ||
            bigstack_alloc_u32(col_type_end, &token_slens) ||
            bigstack_calloc_w(pheno_ctl, &categorical_phenos) ||
            bigstack_calloc_w(pheno_ctl, &quantitative_phenos))) {
      goto LoadPsam_ret_NOMEM;
    }
    const char* missing_catname = g_missing_catname;
    const uint32_t missing_catname_blen = strlen(missing_catname) + 1;
    const uint32_t missing_catname_hval = Hashceil(missing_catname, missing_catname_blen - 1, kCatHtableSize);
    unsigned char* tmp_bigstack_base = g_bigstack_base;
    CatnameLl2** catname_htable = nullptr;
    CatnameLl2** pheno_catname_last = nullptr;
    uintptr_t* total_catname_blens = nullptr;
    for (; TextGetUnsafe2K(&psam_txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      if (unlikely(line_iter[0] == '#')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, psamname);
        goto LoadPsam_ret_MALFORMED_INPUT_WW;
      }
      if (unlikely(raw_sample_ct == 0x7ffffffe)) {
        logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
        goto LoadPsam_ret_MALFORMED_INPUT;
      }
      const char* line_start = line_iter;
      const char* iid_ptr;
      uint32_t iid_slen;
      if (relevant_postfid_col_ct) {
        // bugfix (3 Jul 2018): didn't handle no-other-columns case correctly
        line_iter = TokenLexK0(line_start, col_types, col_skips, relevant_postfid_col_ct, token_ptrs, token_slens);
        if (unlikely(!line_iter)) {
          goto LoadPsam_ret_MISSING_TOKENS;
        }
        iid_ptr = token_ptrs[0];
        iid_slen = token_slens[0];
      } else {
        iid_ptr = line_start;
        iid_slen = CurTokenEnd(line_start) - line_start;
      }
      uint32_t fid_slen = 2;
      if (fid_present) {
        fid_slen = CurTokenEnd(line_start) - line_start;
      }
      const uint32_t sid_slen = sids_present? token_slens[1] : 0;
      const uint32_t paternal_id_slen = paternal_ids_present? token_slens[2] : 1;
      const uint32_t maternal_id_slen = maternal_ids_present? token_slens[3] : 1;
      // phenotypes
      if (!raw_sample_ct) {
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          if (IsCategoricalPhenostrNocsv(token_ptrs[pheno_idx + 5])) {
            SetBit(pheno_idx, categorical_phenos);
          }
        }
        categorical_pheno_ct = PopcountWords(categorical_phenos, pheno_ctl);
        if (categorical_pheno_ct) {
          // initialize hash table
          const uint32_t cat_ul_byte_ct = categorical_pheno_ct * sizeof(intptr_t);
          const uint32_t htable_byte_ct = kCatHtableSize * sizeof(uintptr_t);
          const uintptr_t entry_byte_ct = RoundUpPow2(offsetof(CatnameLl2, str) + missing_catname_blen, sizeof(intptr_t));
          if (unlikely(S_CAST(uintptr_t, tmp_bigstack_end - tmp_bigstack_base) < htable_byte_ct + categorical_pheno_ct * entry_byte_ct + 2 * cat_ul_byte_ct)) {
            goto LoadPsam_ret_NOMEM;
          }
          pheno_catname_last = R_CAST(CatnameLl2**, tmp_bigstack_base);
          tmp_bigstack_base += cat_ul_byte_ct;
          total_catname_blens = R_CAST(uintptr_t*, tmp_bigstack_base);
          tmp_bigstack_base += cat_ul_byte_ct;
          ZeroWArr(categorical_pheno_ct, total_catname_blens);
          catname_htable = R_CAST(CatnameLl2**, tmp_bigstack_base);
          tmp_bigstack_base += htable_byte_ct;
          for (uint32_t uii = 0; uii != kCatHtableSize; ++uii) {
            catname_htable[uii] = nullptr;
          }
          uint32_t cur_hval = missing_catname_hval;
          for (uint32_t cat_pheno_idx = 0; cat_pheno_idx != categorical_pheno_ct; ++cat_pheno_idx) {
            CatnameLl2* new_entry = R_CAST(CatnameLl2*, tmp_bigstack_base);
            tmp_bigstack_base += entry_byte_ct;
            pheno_catname_last[cat_pheno_idx] = new_entry;
            new_entry->cat_idx = 0;
            new_entry->htable_next = nullptr;
            new_entry->pheno_next = nullptr;
            memcpy(new_entry->str, missing_catname, missing_catname_blen);
            catname_htable[cur_hval++] = new_entry;
            if (cur_hval == kCatHtableSize) {
              cur_hval = 0;
            }
          }
        }
      }
      // 1 extra byte for tab between FID and IID; this gets absorbed into
      // the "+ sizeof(intptr_t)" at the end, since that would normally be
      // "+ (sizeof(intptr_t) - 1)"
      // bugfix: pheno_ct * sizeof(intptr_t) -> pheno_ct * 8
      const uint32_t alloc_byte_ct = sizeof(PsamInfoLl) + sizeof(intptr_t) + RoundDownPow2(fid_slen + iid_slen + sid_slen + paternal_id_slen + maternal_id_slen + pheno_ct * 8, sizeof(intptr_t));
      PsamInfoLl* new_psam_info = R_CAST(PsamInfoLl*, tmp_bigstack_base);
      tmp_bigstack_base += alloc_byte_ct;
      if (unlikely(tmp_bigstack_base > tmp_bigstack_end)) {
        goto LoadPsam_ret_NOMEM;
      }
      new_psam_info->next = psam_info_reverse_ll;
      char* sample_id_storage = R_CAST(char*, &(new_psam_info->vardata[pheno_ct * 8]));
      char* ss_iter = sample_id_storage;
      if (fid_present) {
        ss_iter = memcpya(ss_iter, line_start, fid_slen);
      } else {
        *ss_iter++ = '0';
      }
      *ss_iter++ = '\t';
      psam_info_reverse_ll = new_psam_info;
      if (unlikely((iid_slen == 1) && (iid_ptr[0] == '0'))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid IID '0' on line %" PRIuPTR " of %s.\n", line_idx, psamname);
        goto LoadPsam_ret_MALFORMED_INPUT_WW;
      }
      ss_iter = memcpya(ss_iter, iid_ptr, iid_slen);
      const uint32_t sample_id_slen = ss_iter - sample_id_storage;
      if (sample_id_slen >= max_sample_id_blen) {
        max_sample_id_blen = sample_id_slen + 1;
      }
      new_psam_info->sample_id_slen = sample_id_slen;
      if (sids_present) {
        ss_iter = memcpya(ss_iter, token_ptrs[1], sid_slen);
        if (sid_slen >= max_sid_blen) {
          max_sid_blen = sid_slen + 1;
        }
        // bugfix (3 Oct 2019): forgot this
        new_psam_info->sid_slen = sid_slen;
      }
      if (paternal_ids_present) {
        if (paternal_id_slen >= max_paternal_id_blen) {
          max_paternal_id_blen = paternal_id_slen + 1;
        }
        ss_iter = memcpya(ss_iter, token_ptrs[2], paternal_id_slen);
      } else {
        *ss_iter++ = '0';
      }
      new_psam_info->paternal_id_slen = paternal_id_slen;
      if (maternal_ids_present) {
        if (maternal_id_slen >= max_maternal_id_blen) {
          max_maternal_id_blen = maternal_id_slen + 1;
        }
        ss_iter = memcpya(ss_iter, token_ptrs[3], maternal_id_slen);
      } else {
        *ss_iter++ = '0';
      }
      new_psam_info->maternal_id_slen = maternal_id_slen;
      uint32_t cur_sex_code = 0;
      // accept 'M'/'F'/'m'/'f' since that's more readable without being any
      // less efficient
      // don't accept "male"/"female", that's overkill
      if (sex_present && (token_slens[4] == 1)) {
        const unsigned char sex_ucc = token_ptrs[4][0];
        const unsigned char sex_ucc_upcase = sex_ucc & 0xdfU;
        if ((sex_ucc == '1') || (sex_ucc_upcase == 'M')) {
          cur_sex_code = 1;
        } else if ((sex_ucc == '2') || (sex_ucc_upcase == 'F')) {
          cur_sex_code = 2;
        }
      }
      new_psam_info->sex_code = cur_sex_code;
      // phenotypes
      unsigned char* pheno_data = new_psam_info->vardata;
      uint32_t cat_pheno_idx = 0;
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        const uint32_t col_type_idx = pheno_idx + 5;
        const char* cur_phenostr = token_ptrs[col_type_idx];
        double dxx;
        const char* cur_phenostr_end = ScanadvDouble(cur_phenostr, &dxx);
        if (!cur_phenostr_end) {
          // possible todo: defend against out-of-range numbers like 1e1000;
          // right now they get treated as categorical variables...
          const uint32_t slen = token_slens[col_type_idx];
          if (IsNanStr(cur_phenostr, slen)) {
            dxx = missing_phenod;
          } else {
            if (unlikely(!IsSet(categorical_phenos, pheno_idx))) {
              assert(psam_info_reverse_ll->next);
              const uint32_t is_second_relevant_line = !(psam_info_reverse_ll->next->next);
              logerrprintfww("Error: '%s' entry on line %" PRIuPTR " of %s is categorical, while %s not.\n", &((*pheno_names_ptr)[pheno_idx * max_pheno_name_blen]), line_idx, psamname, is_second_relevant_line? "an earlier entry is" : "earlier entries are");
              goto LoadPsam_ret_INCOMPATIBLE_PHENOSTRS;
            }
            if (unlikely(slen > kMaxIdSlen)) {
              logerrputs("Error: Categorical phenotypes are limited to " MAX_ID_SLEN_STR " characters.\n");
              goto LoadPsam_ret_MALFORMED_INPUT;
            }
            uint32_t hashval = Hashceil(cur_phenostr, slen, kCatHtableSize) + cat_pheno_idx;
            if (hashval >= kCatHtableSize) {
              hashval -= kCatHtableSize;
            }
            uintptr_t htable_idx = 0;
            for (CatnameLl2** cur_entry_ptr = &(catname_htable[hashval]); ; ) {
              CatnameLl2* cur_entry = *cur_entry_ptr;
              if (!cur_entry) {
                const uint32_t entry_byte_ct = RoundUpPow2(offsetof(CatnameLl2, str) + slen + 1, sizeof(intptr_t));
                htable_idx = pheno_catname_last[cat_pheno_idx]->cat_idx + 1;
                CatnameLl2* new_entry = R_CAST(CatnameLl2*, tmp_bigstack_base);
                tmp_bigstack_base += entry_byte_ct;
                if (unlikely(tmp_bigstack_base > tmp_bigstack_end)) {
                  goto LoadPsam_ret_NOMEM;
                }
                new_entry->htable_next = nullptr;
                new_entry->pheno_next = nullptr;
                pheno_catname_last[cat_pheno_idx]->pheno_next = new_entry;
                pheno_catname_last[cat_pheno_idx] = new_entry;
                *cur_entry_ptr = new_entry;
                new_entry->cat_idx = htable_idx;
                memcpyx(new_entry->str, cur_phenostr, slen, '\0');
                total_catname_blens[cat_pheno_idx] += slen + 1;
                break;
              }
              // safe since we guarantee kMaxIdSlen spare bytes at the end
              // of bigstack
              if (memequal(cur_entry->str, cur_phenostr, slen) && (!cur_entry->str[slen])) {
                htable_idx = cur_entry->cat_idx;
                break;
              }
              cur_entry_ptr = &(cur_entry->htable_next);
            }
            // don't bother writing top 4 bytes in 32-bit build
            memcpy(&(pheno_data[pheno_idx * 8]), &htable_idx, sizeof(intptr_t));
            ++cat_pheno_idx;
            continue;
          }
        } else {
          if (unlikely(!IsSpaceOrEoln(*cur_phenostr_end))) {
            cur_phenostr_end = CurTokenEnd(cur_phenostr_end);
            *K_CAST(char*, cur_phenostr_end) = '\0';
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid numeric token '%s' on line %" PRIuPTR " of %s.\n", cur_phenostr, line_idx, psamname);
            goto LoadPsam_ret_MALFORMED_INPUT_WW;
          }
        }
        if (unlikely(IsSet(categorical_phenos, pheno_idx))) {
          assert(psam_info_reverse_ll->next);
          const uint32_t is_second_relevant_line = !(psam_info_reverse_ll->next->next);
          logerrprintfww("Error: '%s' entry on line %" PRIuPTR " of %s is numeric/'NA', while %s categorical.\n", &((*pheno_names_ptr)[pheno_idx * max_pheno_name_blen]), line_idx, psamname, is_second_relevant_line? "an earlier entry is" : "earlier entries are");
          goto LoadPsam_ret_INCOMPATIBLE_PHENOSTRS;
        }
        if (!IsSet(quantitative_phenos, pheno_idx)) {
          if ((dxx != missing_phenod) && (dxx != pheno_ctrld) && (dxx != pheno_cased) && (dxx != 0.0)) {
            SetBit(pheno_idx, quantitative_phenos);
          }
        }
        memcpy(&(pheno_data[pheno_idx * 8]), &dxx, sizeof(double));
      }
      ++raw_sample_ct;
    }
    if (unlikely(TextStreamErrcode2(&psam_txs, &reterr))) {
      goto LoadPsam_ret_TSTREAM_FAIL;
    }
    reterr = kPglRetSuccess;
    if (unlikely((max_sample_id_blen > 2 * kMaxIdBlen) || (max_paternal_id_blen > kMaxIdBlen) || (max_maternal_id_blen > kMaxIdBlen))) {
      logerrputs("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto LoadPsam_ret_MALFORMED_INPUT;
    }
    if (unlikely(max_sid_blen > kMaxIdBlen)) {
      logerrputs("Error: SIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto LoadPsam_ret_MALFORMED_INPUT;
    }
    if (unlikely(!raw_sample_ct)) {
      logerrprintfww("Error: No samples in %s.\n", psamname);
      goto LoadPsam_ret_MALFORMED_INPUT;
    }
    const uintptr_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    BigstackBaseSet(tmp_bigstack_base);

    if (pheno_ct) {
      if (unlikely(pgl_malloc(pheno_ct * sizeof(PhenoCol), &pheno_cols))) {
        goto LoadPsam_ret_NOMEM;
      }
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        // ensure cleanup works if initialization fails in the middle
        pheno_cols[pheno_idx].nonmiss = nullptr;
      }
      uint32_t cat_pheno_idx = 0;
      PhenoCol* pheno_cols_iter = pheno_cols;
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        const uintptr_t nonmiss_vec_ct = BitCtToVecCt(raw_sample_ct);
        const uint32_t is_categorical = IsSet(categorical_phenos, pheno_idx);
        const uint32_t is_qt = IsSet(quantitative_phenos, pheno_idx);
        uintptr_t data_vec_ct = 0;
        uintptr_t catname_vec_ct = 0;
        uintptr_t catname_storage_vec_ct = 0;
        uint32_t nonnull_catname_ct = 0;
        if (!is_categorical) {
          pheno_cols_iter->category_names = nullptr;
          pheno_cols_iter->type_code = S_CAST(PhenoDtype, is_qt);
          pheno_cols_iter->nonnull_category_ct = 0;
          if (is_qt) {
            data_vec_ct = DblCtToVecCt(raw_sample_ct);
          } else {
            data_vec_ct = nonmiss_vec_ct;
          }
        } else {
          nonnull_catname_ct = pheno_catname_last[cat_pheno_idx]->cat_idx;
          data_vec_ct = Int32CtToVecCt(raw_sample_ct);
          catname_vec_ct = WordCtToVecCt(nonnull_catname_ct + 1);
          catname_storage_vec_ct = DivUp(total_catname_blens[cat_pheno_idx], kBytesPerVec);
          pheno_cols_iter->type_code = kPhenoDtypeCat;
          pheno_cols_iter->nonnull_category_ct = nonnull_catname_ct;
        }
        // pheno_cols_iter->nonmiss = nullptr;
        uintptr_t* new_pheno_data_iter;
        if (unlikely(vecaligned_malloc((nonmiss_vec_ct + data_vec_ct + catname_vec_ct + catname_storage_vec_ct) * kBytesPerVec, &new_pheno_data_iter))) {
          goto LoadPsam_ret_NOMEM;
        }
        pheno_cols_iter->nonmiss = new_pheno_data_iter;
        ZeroWArr(nonmiss_vec_ct * kWordsPerVec, new_pheno_data_iter);
        new_pheno_data_iter = &(new_pheno_data_iter[nonmiss_vec_ct * kWordsPerVec]);
        if (is_categorical) {
          pheno_cols_iter->data.cat = R_CAST(uint32_t*, new_pheno_data_iter);
          new_pheno_data_iter = &(new_pheno_data_iter[data_vec_ct * kWordsPerVec]);
          const char** cur_name_ptrs = R_CAST(const char**, new_pheno_data_iter);
          pheno_cols_iter->category_names = cur_name_ptrs;
          *cur_name_ptrs++ = missing_catname;
          char* name_storage_iter = R_CAST(char*, &(new_pheno_data_iter[catname_vec_ct * kWordsPerVec]));
          uint32_t cur_hval = missing_catname_hval + cat_pheno_idx;
          if (cur_hval >= kCatHtableSize) {
            cur_hval -= kCatHtableSize;
          }
          // make this point to the "NONE" entry for the current phenotype,
          // which starts the linked list
          CatnameLl2* catname_entry_ptr = catname_htable[cur_hval];

          for (uint32_t catname_idx = 0; catname_idx != nonnull_catname_ct; ++catname_idx) {
            catname_entry_ptr = catname_entry_ptr->pheno_next;
            char* cur_name_start = name_storage_iter;
            name_storage_iter = strcpyax(name_storage_iter, catname_entry_ptr->str, '\0');
            *cur_name_ptrs++ = cur_name_start;
          }
          ++cat_pheno_idx;
        } else if (!is_qt) {
          pheno_cols_iter->data.cc = new_pheno_data_iter;
          ZeroWArr(nonmiss_vec_ct * kWordsPerVec, new_pheno_data_iter);
        } else {
          pheno_cols_iter->data.qt = R_CAST(double*, new_pheno_data_iter);
        }
        ++pheno_cols_iter;
      }
    }
    // real allocations start here
    // could make these cacheline-aligned?
    g_bigstack_end = bigstack_end_mark;
    const uint32_t aligned_wct = BitCtToAlignedWordCt(raw_sample_ct);
    if (unlikely(
            bigstack_end_alloc_c(raw_sample_ct * max_sample_id_blen, &(piip->sii.sample_ids)) ||
            bigstack_end_alloc_c(raw_sample_ct * max_paternal_id_blen, &(piip->parental_id_info.paternal_ids)) ||
            bigstack_end_alloc_c(raw_sample_ct * max_maternal_id_blen, &(piip->parental_id_info.maternal_ids)) ||
            bigstack_end_alloc_w(raw_sample_ctl, sample_include_ptr) ||
            bigstack_end_calloc_w(aligned_wct, founder_info_ptr) ||
            bigstack_end_calloc_w(aligned_wct, sex_nm_ptr) ||
            bigstack_end_calloc_w(aligned_wct, sex_male_ptr))) {
      goto LoadPsam_ret_NOMEM;
    }
    if (max_sid_blen) {
      if (unlikely(bigstack_end_alloc_c(raw_sample_ct * max_sid_blen, &(piip->sii.sids)))) {
        goto LoadPsam_ret_NOMEM;
      }
      if (!sids_present) {
        // Possible with --update-ids.
        char* sids_iter = piip->sii.sids;
        for (uint32_t sample_idx = 0; sample_idx != raw_sample_ct; ++sample_idx) {
          strcpy_k(sids_iter, "0");
          sids_iter = &(sids_iter[max_sid_blen]);
        }
      }
    }
    bigstack_end_mark = g_bigstack_end;
    SetAllBits(raw_sample_ct, *sample_include_ptr);
    // make FillInterleavedMaskVec() work by default
    ZeroTrailingWords(raw_sample_ctl, *sample_include_ptr);
    ZeroTrailingWords(raw_sample_ctl, *founder_info_ptr);
    ZeroTrailingWords(raw_sample_ctl, *sex_male_ptr);
    *raw_sample_ct_ptr = raw_sample_ct;
    piip->sii.max_sample_id_blen = max_sample_id_blen;
    piip->sii.max_sid_blen = max_sid_blen;
    piip->parental_id_info.max_paternal_id_blen = max_paternal_id_blen;
    piip->parental_id_info.max_maternal_id_blen = max_maternal_id_blen;
    *max_pheno_name_blen_ptr = max_pheno_name_blen;
    char* sample_ids = piip->sii.sample_ids;
    char* sids = piip->sii.sids;
    char* paternal_ids = piip->parental_id_info.paternal_ids;
    char* maternal_ids = piip->parental_id_info.maternal_ids;
    uintptr_t* founder_info = *founder_info_ptr;
    uintptr_t* sex_nm = *sex_nm_ptr;
    uintptr_t* sex_male = *sex_male_ptr;
    uint32_t sample_uidx = raw_sample_ct;
    while (sample_uidx) {
      --sample_uidx;
      unsigned char* cur_vardata = psam_info_reverse_ll->vardata;
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        if (IsSet(categorical_phenos, pheno_idx)) {
          uint32_t cur_cat;
          memcpy(&cur_cat, &(cur_vardata[pheno_idx * 8]), sizeof(int32_t));
          pheno_cols[pheno_idx].data.cat[sample_uidx] = cur_cat;
          if (cur_cat) {
            SetBit(sample_uidx, pheno_cols[pheno_idx].nonmiss);
          }
        } else {
          double dxx;
          memcpy(&dxx, &(cur_vardata[pheno_idx * 8]), sizeof(double));
          if (IsSet(quantitative_phenos, pheno_idx)) {
            if (dxx != missing_phenod) {
              SetBit(sample_uidx, pheno_cols[pheno_idx].nonmiss);
              pheno_cols[pheno_idx].data.qt[sample_uidx] = dxx;
            }
          } else {
            if (dxx == pheno_cased) {
              SetBit(sample_uidx, pheno_cols[pheno_idx].data.cc);
              SetBit(sample_uidx, pheno_cols[pheno_idx].nonmiss);
            } else if (dxx == pheno_ctrld) {
              SetBit(sample_uidx, pheno_cols[pheno_idx].nonmiss);
            }
          }
        }
      }
      const uint32_t sample_id_slen = psam_info_reverse_ll->sample_id_slen;
      const uint32_t paternal_id_slen = psam_info_reverse_ll->paternal_id_slen;
      const uint32_t maternal_id_slen = psam_info_reverse_ll->maternal_id_slen;
      const uint32_t sex_code = psam_info_reverse_ll->sex_code;
      char* cur_sample_id = R_CAST(char*, &(cur_vardata[pheno_ct * 8]));
      memcpyx(&(sample_ids[sample_uidx * max_sample_id_blen]), cur_sample_id, sample_id_slen, '\0');
      char* cur_paternal_id = &(cur_sample_id[sample_id_slen]);
      if (sids_present) {
        const uint32_t sid_slen = psam_info_reverse_ll->sid_slen;
        char* cur_sid = cur_paternal_id;
        memcpyx(&(sids[sample_uidx * max_sid_blen]), cur_sid, sid_slen, '\0');
        cur_paternal_id = &(cur_sid[sid_slen]);
      }
      memcpyx(&(paternal_ids[sample_uidx * max_paternal_id_blen]), cur_paternal_id, paternal_id_slen, '\0');
      char* cur_maternal_id = &(cur_paternal_id[paternal_id_slen]);
      if ((paternal_id_slen == 1) && (maternal_id_slen == 1) && (cur_paternal_id[0] == '0') && (cur_maternal_id[0] == '0')) {
        SetBit(sample_uidx, founder_info);
      }
      memcpyx(&(maternal_ids[sample_uidx * max_maternal_id_blen]), cur_maternal_id, maternal_id_slen, '\0');
      if (sex_code) {
        SetBit(sample_uidx, sex_nm);
        if (sex_code == 1) {
          SetBit(sample_uidx, sex_male);
        }
      }
      psam_info_reverse_ll = psam_info_reverse_ll->next;
    }
    if (fid_present) {
      piip->sii.flags |= kfSampleIdFidPresent;
    }
    // special case: if there's exactly one phenotype and it's all-missing,
    // discard it
    if ((pheno_ct == 1) && (!PopcountWords(pheno_cols[0].nonmiss, raw_sample_ctl))) {
      free(*pheno_names_ptr);
      *pheno_names_ptr = nullptr;
      CleanupPhenoCols(1, pheno_cols);
      *pheno_cols_ptr = nullptr;
      *pheno_ct_ptr = 0;
    } else {
      *pheno_cols_ptr = pheno_cols;
      *pheno_ct_ptr = pheno_ct;
    }
  }
  while (0) {
  LoadPsam_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadPsam_ret_TSTREAM_FAIL:
    TextStreamErrPrint(psamname, &psam_txs);
    break;
  LoadPsam_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  LoadPsam_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  LoadPsam_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  LoadPsam_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  LoadPsam_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, psamname);
    reterr = kPglRetMalformedInput;
    break;
  LoadPsam_ret_INCOMPATIBLE_PHENOSTRS:
    logerrputs("(Case/control and quantitative phenotypes must all be numeric/'NA'.\nCategorical phenotypes cannot be 'NA'--use e.g. 'NONE' to represent missing\ncategorical values instead--or start with a number.)\n");
    reterr = kPglRetMalformedInput;
    break;
  }
 LoadPsam_ret_1:
  CleanupTextStream2(psamname, &psam_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  if (reterr) {
    if (*pheno_names_ptr) {
      free(*pheno_names_ptr);
      *pheno_names_ptr = nullptr;
    }
    CleanupPhenoCols(pheno_ct, pheno_cols);
    *pheno_ct_ptr = 0;
    *pheno_cols_ptr = nullptr;
  }
  return reterr;
}


typedef struct PhenoInfoLlStruct {
  NONCOPYABLE(PhenoInfoLlStruct);
  // for categorical phenotypes, phenodata entry should be reinterpreted as
  // uint32_t
  struct PhenoInfoLlStruct* next;
  uint32_t sample_uidx;
  double phenodata[];
} PhenoInfoLl;

// also for loading covariates.  set affection_01 to 2 to prohibit case/control
// and make unnamed variables start with "COVAR" instead of "PHENO"
PglErr LoadPhenos(const char* pheno_fname, const RangeList* pheno_range_list_ptr, const uintptr_t* sample_include, const char* sample_ids, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, int32_t missing_pheno, uint32_t affection_01, uint32_t iid_only, uint32_t numeric_ranges, uint32_t max_thread_ct, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* pheno_names = nullptr;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream pheno_txs;
  PreinitTextStream(&pheno_txs);
  {
    if (!sample_ct) {
      goto LoadPhenos_ret_1;
    }
    reterr = SizeAndInitTextStream(pheno_fname, bigstack_left() / 4, MAXV(max_thread_ct - 1, 1), &pheno_txs);
    if (unlikely(reterr)) {
      goto LoadPhenos_ret_TSTREAM_FAIL;
    }
    const char* line_iter;
    do {
      ++line_idx;
      line_iter = TextGet(&pheno_txs);
      if (!line_iter) {
        if (likely(!TextStreamErrcode2(&pheno_txs, &reterr))) {
          goto LoadPhenos_ret_1;
        }
        goto LoadPhenos_ret_TSTREAM_FAIL;
      }
    } while ((line_iter[0] == '#') && (((line_iter[1] != 'F') && (line_iter[1] != 'I')) || (!memequal_k(&(line_iter[2]), "ID", 2)) || ((ctou32(line_iter[4]) > ' ') && (line_iter[4] != ','))));
    if (line_iter[0] == '#') {
      ++line_iter;
    }
    const uint32_t old_pheno_ct = *pheno_ct_ptr;
    const uintptr_t old_max_pheno_name_blen = *max_pheno_name_blen_ptr;
    uintptr_t max_pheno_name_blen = old_max_pheno_name_blen;
    uint32_t comma_delim = 0;
    XidMode xid_mode = (line_iter[0] == 'I')? kfXidModeIid : kfXidModeFidIid;
    uint32_t* col_types = nullptr;
    uint32_t* col_skips = nullptr;
    uint32_t new_pheno_ct;
    uint32_t final_pheno_ct;
    uintptr_t final_pheno_names_byte_ct;
    if (((line_iter[0] == 'F') || (xid_mode != kfXidModeFidIid)) && memequal_k(&(line_iter[1]), "ID", 2) && ((ctou32(line_iter[3]) <= 32) || (line_iter[3] == ','))) {
      // treat this as a header line
      // autodetect CSV vs. space/tab-delimited
      // (note that we don't permit CSVs without header lines)
      comma_delim = (line_iter[3] == ',');
      const char* linebuf_iter = FirstNonTspace(&(line_iter[3 + comma_delim]));
      if (unlikely(IsEolnKns(*linebuf_iter))) {
        goto LoadPhenos_ret_MISSING_TOKENS;
      }
      const char* iid_end;
      if (xid_mode == kfXidModeFidIid) {
        if (iid_only) {
          snprintf(g_logbuf, kLogbufSize, "Error: \"--%s iid-only\" file has a FID column.\n", (affection_01 == 2)? "covar" : "pheno");
          goto LoadPhenos_ret_INCONSISTENT_INPUT_2;
        }
        iid_end = CommaOrTspaceTokenEnd(linebuf_iter, comma_delim);
        const uintptr_t token_slen = iid_end - linebuf_iter;
        if (unlikely(!strequal_k(linebuf_iter, "IID", token_slen))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Second column header in %s must be 'IID'.\n", pheno_fname);
          goto LoadPhenos_ret_MALFORMED_INPUT_WW;
        }
        linebuf_iter = CommaOrTspaceFirstToken(iid_end, comma_delim);
      } else {
        iid_end = &(line_iter[3]);
      }
      uint32_t pheno_col_ct = 0;
      const char* pheno_start = linebuf_iter;
      while (linebuf_iter) {
        const char* token_end = CommaOrTspaceTokenEnd(linebuf_iter, comma_delim);
        const uintptr_t token_slen = token_end - linebuf_iter;
        if (max_pheno_name_blen <= token_slen) {
          max_pheno_name_blen = token_slen + 1;
        }
        ++pheno_col_ct;
        linebuf_iter = CommaOrTspaceFirstToken(token_end, comma_delim);
      }
      if (unlikely(max_pheno_name_blen > kMaxIdBlen)) {
        logerrputs("Error: Phenotype/covariate names are limited to " MAX_ID_SLEN_STR " characters.\n");
        goto LoadPhenos_ret_MALFORMED_INPUT;
      }
      if (pheno_range_list_ptr->names && pheno_col_ct) {
        // bugfix (20 Oct 2017): forgot to make error message different in
        // --covar case
        uintptr_t* bitarr;
        if (numeric_ranges) {
          const uint32_t leading_col_ct = GetXidColCt(xid_mode);
          if (unlikely(bigstack_calloc_w(BitCtToWordCt(pheno_col_ct), &bitarr))) {
            goto LoadPhenos_ret_NOMEM;
          }
          if (unlikely(NumericRangeListToBitarr(pheno_range_list_ptr, pheno_col_ct, leading_col_ct + 1, 0, bitarr))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s-col-nums argument for %s.\n", (affection_01 == 2)? "covar" : "pheno", pheno_fname);
            goto LoadPhenos_ret_INCONSISTENT_INPUT_WW;
          }
        } else {
          reterr = StringRangeListToBitarrAlloc(pheno_start, pheno_range_list_ptr, (affection_01 == 2)? "covar-name" : "pheno-name", pheno_fname, pheno_col_ct, 0, comma_delim, &bitarr);
          if (unlikely(reterr)) {
            goto LoadPhenos_ret_1;
          }
        }
        new_pheno_ct = PopcountWords(bitarr, BitCtToWordCt(pheno_col_ct));
        if (unlikely(
                bigstack_alloc_u32(new_pheno_ct, &col_types) ||
                bigstack_alloc_u32(new_pheno_ct, &col_skips))) {
          goto LoadPhenos_ret_NOMEM;
        }
        uintptr_t col_uidx_base = 0;
        uintptr_t cur_bits = bitarr[0];
        int32_t prev_col_uidx = -1;
        for (uint32_t col_idx = 0; col_idx != new_pheno_ct; ++col_idx) {
          const uint32_t col_uidx = BitIter1(bitarr, &col_uidx_base, &cur_bits);
          col_types[col_idx] = col_idx;
          col_skips[col_idx] = col_uidx - prev_col_uidx;
          prev_col_uidx = col_uidx;
        }
      } else {
        // usual case, load all phenotypes
        new_pheno_ct = pheno_col_ct;
        if (unlikely(
                bigstack_alloc_u32(new_pheno_ct, &col_types) ||
                bigstack_alloc_u32(new_pheno_ct, &col_skips))) {
          goto LoadPhenos_ret_NOMEM;
        }
        for (uint32_t col_idx = 0; col_idx != pheno_col_ct; ++col_idx) {
          col_types[col_idx] = col_idx;
          col_skips[col_idx] = 1;
        }
      }
      final_pheno_ct = new_pheno_ct + old_pheno_ct;
      final_pheno_names_byte_ct = final_pheno_ct * max_pheno_name_blen;
      if (unlikely(pgl_malloc(final_pheno_names_byte_ct, &pheno_names))) {
        goto LoadPhenos_ret_NOMEM;
      }
      linebuf_iter = iid_end;
      char* pheno_names_iter = &(pheno_names[old_pheno_ct * max_pheno_name_blen]);
      for (uint32_t new_pheno_idx = 0; new_pheno_idx != new_pheno_ct; ++new_pheno_idx) {
        linebuf_iter = CommaOrTspaceNextTokenMult(linebuf_iter, col_skips[new_pheno_idx], comma_delim);
        const char* token_end = CommaOrTspaceTokenEnd(linebuf_iter, comma_delim);
        const uint32_t name_slen = token_end - linebuf_iter;
        if (unlikely(IsReservedPhenoName(linebuf_iter, name_slen))) {
          char* write_iter = strcpya_k(g_logbuf, "Error: '");
          // length verified to be <= kMaxIdSlen
          write_iter = memcpya(write_iter, linebuf_iter, name_slen);
          snprintf(write_iter, kLogbufSize - kMaxIdSlen - 16, "' cannot be used as a phenotype/covariate name.\n");
          goto LoadPhenos_ret_MALFORMED_INPUT_2;
        }
        memcpyx(pheno_names_iter, linebuf_iter, name_slen, '\0');
        pheno_names_iter = &(pheno_names_iter[max_pheno_name_blen]);
        linebuf_iter = token_end;
      }
      line_iter = AdvPastDelim(linebuf_iter, '\n');
      ++line_idx;
    } else {
      // no header line
      xid_mode = iid_only? kfXidModeIid : kfXidModeFidIid;
      // don't support comma delimiter here, since we don't have guaranteed
      // leading FID/IID to distinguish it
      const uint32_t col_ct = CountTokens(line_iter);
      if (unlikely(col_ct < 3 - iid_only)) {
        // todo: tolerate col_ct == 2 with --allow-no-phenos
        goto LoadPhenos_ret_MISSING_TOKENS;
      }
      if (pheno_range_list_ptr->names) {
        if (unlikely(!numeric_ranges)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Header line expected in %s, due to --pheno-name/--covar-name. (This line must start with '#FID', 'FID', '#IID', or 'IID'.)\n", pheno_fname);
          goto LoadPhenos_ret_INCONSISTENT_INPUT_WW;
        }
        const uint32_t bitarr_size = col_ct + iid_only - 2;
        uintptr_t* bitarr;
        if (unlikely(bigstack_calloc_w(BitCtToWordCt(bitarr_size), &bitarr))) {
          goto LoadPhenos_ret_NOMEM;
        }
        if (unlikely(NumericRangeListToBitarr(pheno_range_list_ptr, bitarr_size, 3 - iid_only, 0, bitarr))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid --%s-col-nums argument for %s.\n", (affection_01 == 2)? "covar" : "pheno", pheno_fname);
          goto LoadPhenos_ret_INCONSISTENT_INPUT_WW;
        }
        // this boilerplate may belong in its own function
        new_pheno_ct = PopcountWords(bitarr, BitCtToWordCt(col_ct));
        if (unlikely(
                bigstack_alloc_u32(new_pheno_ct, &col_types) ||
                bigstack_alloc_u32(new_pheno_ct, &col_skips))) {
          goto LoadPhenos_ret_NOMEM;
        }
        uintptr_t col_uidx_base = 0;
        uintptr_t cur_bits = bitarr[0];
        int32_t prev_col_uidx = -1;
        for (uint32_t col_idx = 0; col_idx != new_pheno_ct; ++col_idx) {
          const uint32_t col_uidx = BitIter1(bitarr, &col_uidx_base, &cur_bits);
          col_types[col_idx] = col_idx;
          col_skips[col_idx] = col_uidx - prev_col_uidx;
          prev_col_uidx = col_uidx;
        }
      } else {
        new_pheno_ct = col_ct + iid_only - 2;
        // bugfix (11 Mar 2018): forgot to initialize col_types here
        if (unlikely(
                bigstack_alloc_u32(new_pheno_ct, &col_types) ||
                bigstack_alloc_u32(new_pheno_ct, &col_skips))) {
          goto LoadPhenos_ret_NOMEM;
        }
        for (uint32_t col_idx = 0; col_idx != new_pheno_ct; ++col_idx) {
          col_types[col_idx] = col_idx;
          col_skips[col_idx] = 1;
        }
      }
      final_pheno_ct = new_pheno_ct + old_pheno_ct;
      // bugfix (29 Mar 2018): don't subtract 1 from final_pheno_ct here since
      // names are 1-based
      const uintptr_t max_new_name_blen = 6 + UintSlen(final_pheno_ct);
      if (max_new_name_blen > max_pheno_name_blen) {
        max_pheno_name_blen = max_new_name_blen;
      }
      final_pheno_names_byte_ct = final_pheno_ct * max_pheno_name_blen;
      if (unlikely(pgl_malloc(final_pheno_names_byte_ct, &pheno_names))) {
        goto LoadPhenos_ret_NOMEM;
      }
      const char* default_prefix = (affection_01 == 2)? "COVAR" : "PHENO";
      for (uint32_t pheno_idx = old_pheno_ct; pheno_idx != final_pheno_ct; ) {
        char* write_iter = memcpya_k(&(pheno_names[pheno_idx * max_pheno_name_blen]), default_prefix, 5);
        ++pheno_idx;  // 1-based default names, not 0-based
        write_iter = u32toa(pheno_idx, write_iter);
        *write_iter = '\0';
      }
    }
    if (unlikely(final_pheno_ct > kMaxPhenoCt)) {
      // yeah, yeah, this will never come up
      logerrputs("Error: " PROG_NAME_STR " does not support more than " MAX_PHENO_CT_STR " phenotypes.\n");
      goto LoadPhenos_ret_INCONSISTENT_INPUT;
    }
    for (uint32_t old_pheno_idx = 0; old_pheno_idx != old_pheno_ct; ++old_pheno_idx) {
      strcpy(&(pheno_names[old_pheno_idx * max_pheno_name_blen]), &((*pheno_names_ptr)[old_pheno_idx * old_max_pheno_name_blen]));
    }

    uint32_t tmp_htable_size;
    uint32_t* htable_tmp;
    if (unlikely(HtableGoodSizeAlloc(final_pheno_ct, bigstack_left(), &htable_tmp, &tmp_htable_size))) {
      goto LoadPhenos_ret_NOMEM;
    }
    // possible todo: implement something like --pheno-merge, allow conditional
    // duplication in that case
    const uint32_t duplicate_idx = PopulateStrboxHtable(pheno_names, final_pheno_ct, max_pheno_name_blen, tmp_htable_size, htable_tmp);
    if (unlikely(duplicate_idx)) {
      const char* duplicate_pheno_name = &(pheno_names[duplicate_idx * max_pheno_name_blen]);
      snprintf(g_logbuf, kLogbufSize, "Error: Duplicate phenotype/covariate ID '%s'.\n", duplicate_pheno_name);
      goto LoadPhenos_ret_MALFORMED_INPUT_WW;
    }
    BigstackReset(htable_tmp);

    PhenoCol* new_pheno_cols = S_CAST(PhenoCol*, realloc(*pheno_cols_ptr, final_pheno_ct * sizeof(PhenoCol)));
    if (unlikely(!new_pheno_cols)) {
      goto LoadPhenos_ret_NOMEM;
    }
    // ensure cleanup works if initialization fails in the middle
    for (uint32_t pheno_idx = old_pheno_ct; pheno_idx != final_pheno_ct; ++pheno_idx) {
      new_pheno_cols[pheno_idx].nonmiss = nullptr;
    }
    *pheno_ct_ptr = final_pheno_ct;
    *pheno_cols_ptr = new_pheno_cols;

    // switch to hash table?
    char* sorted_sample_ids;
    uint32_t* sample_id_map;
    // todo: permit duplicates if SIDs are defined
    reterr = CopySortStrboxSubset(sample_include, sample_ids, sample_ct, max_sample_id_blen, 0, 0, 0, &sorted_sample_ids, &sample_id_map);
    if (unlikely(reterr)) {
      goto LoadPhenos_ret_1;
    }
    const uintptr_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    char* id_buf;
    uintptr_t* already_seen;
    if (unlikely(
            bigstack_alloc_c(max_sample_id_blen, &id_buf) ||
            bigstack_calloc_w(raw_sample_ctl, &already_seen))) {
      goto LoadPhenos_ret_NOMEM;
    }

    PhenoInfoLl* pheno_info_reverse_ll = nullptr;
    const uint32_t pheno_info_alloc_byte_ct = sizeof(PhenoInfoLl) + new_pheno_ct * sizeof(double);
    const uint32_t new_pheno_ctl = BitCtToWordCt(new_pheno_ct);
    const double missing_phenod = missing_pheno? S_CAST(double, missing_pheno) : HUGE_VAL;

    // affection_01 can be 2, so don't use u31tod()
    const double pheno_ctrld = S_CAST(int32_t, 1 - affection_01);
    const double pheno_cased = pheno_ctrld + 1.0;
    uint32_t categorical_pheno_ct = 0;
    const char** token_ptrs;
    uint32_t* token_slens;
    uintptr_t* categorical_phenos;
    uintptr_t* quantitative_phenos;
    if (unlikely(
            bigstack_alloc_kcp(new_pheno_ct, &token_ptrs) ||
            bigstack_alloc_u32(new_pheno_ct, &token_slens) ||
            bigstack_calloc_w(new_pheno_ctl, &categorical_phenos) ||
            bigstack_calloc_w(new_pheno_ctl, &quantitative_phenos))) {
      goto LoadPhenos_ret_NOMEM;
    }
    const char* missing_catname = g_missing_catname;
    const uint32_t missing_catname_blen = strlen(missing_catname) + 1;
    const uint32_t missing_catname_hval = Hashceil(missing_catname, missing_catname_blen - 1, kCatHtableSize);
    unsigned char* bigstack_base_copy = g_bigstack_base;
    unsigned char* tmp_bigstack_end = g_bigstack_end;
    CatnameLl2** catname_htable = nullptr;
    CatnameLl2** pheno_catname_last = nullptr;
    uintptr_t* total_catname_blens = nullptr;
    for (; TextGetUnsafe2K(&pheno_txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      if (unlikely(line_iter[0] == '#')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, pheno_fname);
        goto LoadPhenos_ret_MALFORMED_INPUT_WW;
      }
      uint32_t sample_uidx;
      if (SortedXidboxReadFind(sorted_sample_ids, sample_id_map, max_sample_id_blen, sample_ct, comma_delim, xid_mode, &line_iter, &sample_uidx, id_buf)) {
        if (unlikely(!line_iter)) {
          goto LoadPhenos_ret_MISSING_TOKENS;
        }
        continue;
      }
      if (unlikely(IsSet(already_seen, sample_uidx))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID in %s.\n", pheno_fname);
        goto LoadPhenos_ret_MALFORMED_INPUT_WW;
      }
      SetBit(sample_uidx, already_seen);
      if (!comma_delim) {
        line_iter = TokenLexK(line_iter, col_types, col_skips, new_pheno_ct, token_ptrs, token_slens);
      } else {
        line_iter = CsvLexK(line_iter, col_types, col_skips, new_pheno_ct, token_ptrs, token_slens);
      }
      if (unlikely(!line_iter)) {
        goto LoadPhenos_ret_MISSING_TOKENS;
      }
      if (!pheno_info_reverse_ll) {
        // first relevant line, detect categorical phenotypes
        for (uint32_t new_pheno_idx = 0; new_pheno_idx != new_pheno_ct; ++new_pheno_idx) {
          if (IsCategoricalPhenostr(token_ptrs[new_pheno_idx])) {
            SetBit(new_pheno_idx, categorical_phenos);
          } else if (affection_01 == 2) {
            SetBit(new_pheno_idx, quantitative_phenos);
          }
        }
        categorical_pheno_ct = PopcountWords(categorical_phenos, new_pheno_ctl);
        if (categorical_pheno_ct) {
          // initialize hash table
          const uint32_t cat_ul_byte_ct = categorical_pheno_ct * sizeof(intptr_t);
          const uint32_t htable_byte_ct = kCatHtableSize * sizeof(uintptr_t);
          const uintptr_t entry_byte_ct = RoundUpPow2(offsetof(CatnameLl2, str) + missing_catname_blen, sizeof(intptr_t));

          if (unlikely(S_CAST(uintptr_t, tmp_bigstack_end - bigstack_base_copy) < htable_byte_ct + categorical_pheno_ct * entry_byte_ct + 2 * cat_ul_byte_ct)) {
            goto LoadPhenos_ret_NOMEM;
          }
          tmp_bigstack_end -= cat_ul_byte_ct;
          total_catname_blens = R_CAST(uintptr_t*, tmp_bigstack_end);
          tmp_bigstack_end -= cat_ul_byte_ct;
          pheno_catname_last = R_CAST(CatnameLl2**, tmp_bigstack_end);
          ZeroWArr(categorical_pheno_ct, total_catname_blens);
          tmp_bigstack_end -= htable_byte_ct;
          catname_htable = R_CAST(CatnameLl2**, tmp_bigstack_end);
          ZeroPtrArr(kCatHtableSize, catname_htable);
          uint32_t cur_hval = missing_catname_hval;
          for (uint32_t cat_pheno_idx = 0; cat_pheno_idx != categorical_pheno_ct; ++cat_pheno_idx) {
            tmp_bigstack_end -= entry_byte_ct;
            CatnameLl2* new_entry = R_CAST(CatnameLl2*, tmp_bigstack_end);
            pheno_catname_last[cat_pheno_idx] = new_entry;
            new_entry->cat_idx = 0;
            new_entry->htable_next = nullptr;
            new_entry->pheno_next = nullptr;
            memcpy(new_entry->str, missing_catname, missing_catname_blen);
            catname_htable[cur_hval++] = new_entry;
            if (cur_hval == kCatHtableSize) {
              cur_hval = 0;
            }
          }
        }
      }
      if (unlikely(S_CAST(uintptr_t, tmp_bigstack_end - bigstack_base_copy) < pheno_info_alloc_byte_ct)) {
        goto LoadPhenos_ret_NOMEM;
      }
      tmp_bigstack_end -= pheno_info_alloc_byte_ct;

      PhenoInfoLl* new_pheno_info = R_CAST(PhenoInfoLl*, tmp_bigstack_end);
      new_pheno_info->next = pheno_info_reverse_ll;
      new_pheno_info->sample_uidx = sample_uidx;
      double* pheno_data = new_pheno_info->phenodata;
      uint32_t cat_pheno_idx = 0;
      for (uint32_t new_pheno_idx = 0; new_pheno_idx != new_pheno_ct; ++new_pheno_idx) {
        const char* cur_phenostr = token_ptrs[new_pheno_idx];
        double dxx;
        const char* cur_phenostr_end = ScanadvDouble(cur_phenostr, &dxx);
        if (!cur_phenostr_end) {
          const uint32_t slen = token_slens[new_pheno_idx];
          if (IsNanStr(cur_phenostr, slen)) {
            // note that, in CSVs, empty string is interpreted as a missing
            // non-categorical phenotype; explicit "NONE" is needed to denote
            // a missing category
            dxx = missing_phenod;
          } else {
            if (unlikely(!IsSet(categorical_phenos, new_pheno_idx))) {
              assert(pheno_info_reverse_ll);
              const uint32_t is_second_relevant_line = !(pheno_info_reverse_ll->next);
              logerrprintfww("Error: '%s' entry on line %" PRIuPTR " of %s is categorical, while %s not.\n", &(pheno_names[(old_pheno_ct + new_pheno_idx) * max_pheno_name_blen]), line_idx, pheno_fname, is_second_relevant_line? "an earlier entry is" : "earlier entries are");
              goto LoadPhenos_ret_INCOMPATIBLE_PHENOSTRS;
            }
            uint32_t hashval;
            hashval = Hashceil(cur_phenostr, slen, kCatHtableSize) + cat_pheno_idx;
            if (hashval >= kCatHtableSize) {
              hashval -= kCatHtableSize;
            }
            uintptr_t htable_idx = 0;
            CatnameLl2** cur_entry_ptr = &(catname_htable[hashval]);
            while (1) {
              CatnameLl2* cur_entry = *cur_entry_ptr;
              if (!cur_entry) {
                const uint32_t entry_byte_ct = RoundUpPow2(offsetof(CatnameLl2, str) + slen + 1, sizeof(intptr_t));
                if (unlikely(S_CAST(uintptr_t, tmp_bigstack_end - bigstack_base_copy) < entry_byte_ct)) {
                  goto LoadPhenos_ret_NOMEM;
                }
                tmp_bigstack_end -= entry_byte_ct;
                htable_idx = pheno_catname_last[cat_pheno_idx]->cat_idx + 1;
                CatnameLl2* new_entry = R_CAST(CatnameLl2*, tmp_bigstack_end);
                new_entry->htable_next = nullptr;
                new_entry->pheno_next = nullptr;
                pheno_catname_last[cat_pheno_idx]->pheno_next = new_entry;
                pheno_catname_last[cat_pheno_idx] = new_entry;
                *cur_entry_ptr = new_entry;
                new_entry->cat_idx = htable_idx;
                memcpyx(new_entry->str, cur_phenostr, slen, '\0');
                total_catname_blens[cat_pheno_idx] += slen + 1;
                break;
              }
              // safe since hash table entries are in the middle of bigstack
              if (memequal(cur_entry->str, cur_phenostr, slen) && (!cur_entry->str[slen])) {
                htable_idx = cur_entry->cat_idx;
                break;
              }
              cur_entry_ptr = &(cur_entry->htable_next);
            }
            // don't bother writing top 4 bytes in 32-bit build
            memcpy(&(pheno_data[new_pheno_idx]), &htable_idx, sizeof(intptr_t));
            ++cat_pheno_idx;
            continue;
          }
        } else {
          if (unlikely(!IsSpaceOrEoln(*cur_phenostr_end))) {
            cur_phenostr_end = CurTokenEnd(cur_phenostr_end);
            *K_CAST(char*, cur_phenostr_end) = '\0';
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid numeric token '%s' on line %" PRIuPTR " of %s.\n", cur_phenostr, line_idx, pheno_fname);
            goto LoadPhenos_ret_MALFORMED_INPUT_WW;
          }
        }
        if (unlikely(IsSet(categorical_phenos, new_pheno_idx))) {
          assert(pheno_info_reverse_ll);
          const uint32_t is_second_relevant_line = !(pheno_info_reverse_ll->next);
          logerrprintfww("Error: '%s' entry on line %" PRIuPTR " of %s is numeric/'NA', while %s categorical.\n", &(pheno_names[(old_pheno_ct + new_pheno_idx) * max_pheno_name_blen]), line_idx, pheno_fname, is_second_relevant_line? "an earlier entry is" : "earlier entries are");
          goto LoadPhenos_ret_INCOMPATIBLE_PHENOSTRS;
        }
        if (!IsSet(quantitative_phenos, new_pheno_idx)) {
          if ((dxx != missing_phenod) && (dxx != pheno_ctrld) && (dxx != pheno_cased) && (dxx != 0.0)) {
            SetBit(new_pheno_idx, quantitative_phenos);
          }
        }
        pheno_data[new_pheno_idx] = dxx;
      }
      pheno_info_reverse_ll = new_pheno_info;
    }
    if (unlikely(TextStreamErrcode2(&pheno_txs, &reterr))) {
      goto LoadPhenos_ret_TSTREAM_FAIL;
    }
    reterr = kPglRetSuccess;
    if (unlikely(!pheno_info_reverse_ll)) {
      if (line_idx == 1) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", pheno_fname);
        goto LoadPhenos_ret_MALFORMED_INPUT_WW;
      }
      // could make this a warning, and automatically delete phenotypes?
      logerrprintf("Error: No entries in %s correspond to loaded sample IDs.\n", pheno_fname);
      goto LoadPhenos_ret_INCONSISTENT_INPUT;
    }
    if (new_pheno_ct) {
      const uintptr_t nonmiss_vec_ct = BitCtToVecCt(raw_sample_ct);
      uint32_t cat_pheno_idx = 0;
      PhenoCol* pheno_cols_iter = &(new_pheno_cols[old_pheno_ct]);
      for (uint32_t new_pheno_idx = 0; new_pheno_idx != new_pheno_ct; ++new_pheno_idx) {
        const uint32_t is_categorical = IsSet(categorical_phenos, new_pheno_idx);
        const uint32_t is_qt = IsSet(quantitative_phenos, new_pheno_idx);
        uintptr_t data_vec_ct = 0;
        uintptr_t catname_vec_ct = 0;
        uintptr_t catname_storage_vec_ct = 0;
        uint32_t nonnull_catname_ct = 0;
        if (!is_categorical) {
          pheno_cols_iter->category_names = nullptr;
          pheno_cols_iter->type_code = S_CAST(PhenoDtype, is_qt);
          pheno_cols_iter->nonnull_category_ct = 0;
          if (is_qt) {
            data_vec_ct = DblCtToVecCt(raw_sample_ct);
          } else {
            data_vec_ct = nonmiss_vec_ct;
          }
        } else {
          nonnull_catname_ct = pheno_catname_last[cat_pheno_idx]->cat_idx;
          data_vec_ct = Int32CtToVecCt(raw_sample_ct);
          catname_vec_ct = WordCtToVecCt(nonnull_catname_ct + 1);
          catname_storage_vec_ct = DivUp(total_catname_blens[cat_pheno_idx], kBytesPerVec);
          pheno_cols_iter->type_code = kPhenoDtypeCat;
          pheno_cols_iter->nonnull_category_ct = nonnull_catname_ct;
        }
        // pheno_cols_iter->nonmiss = nullptr;
        uintptr_t* new_pheno_data_iter;
        if (unlikely(vecaligned_malloc((nonmiss_vec_ct + data_vec_ct + catname_vec_ct + catname_storage_vec_ct) * kBytesPerVec, &new_pheno_data_iter))) {
          goto LoadPhenos_ret_NOMEM;
        }
        pheno_cols_iter->nonmiss = new_pheno_data_iter;
        ZeroWArr(nonmiss_vec_ct * kWordsPerVec, new_pheno_data_iter);
        new_pheno_data_iter = &(new_pheno_data_iter[nonmiss_vec_ct * kWordsPerVec]);
        if (is_categorical) {
          // allow nonmiss[] to be ignored in categorical case
          ZeroWArr(data_vec_ct, new_pheno_data_iter);
          pheno_cols_iter->data.cat = R_CAST(uint32_t*, new_pheno_data_iter);
          new_pheno_data_iter = &(new_pheno_data_iter[data_vec_ct * kWordsPerVec]);
          const char** cur_name_ptrs = R_CAST(const char**, new_pheno_data_iter);
          pheno_cols_iter->category_names = cur_name_ptrs;
          *cur_name_ptrs++ = g_missing_catname;
          char* name_storage_iter = R_CAST(char*, &(new_pheno_data_iter[catname_vec_ct * kWordsPerVec]));
          uint32_t cur_hval = missing_catname_hval + cat_pheno_idx;
          if (cur_hval >= kCatHtableSize) {
            cur_hval -= kCatHtableSize;
          }
          // make this point to the "NONE" entry for the current phenotype,
          // which starts the linked list
          CatnameLl2* catname_entry_ptr = catname_htable[cur_hval];

          for (uint32_t catname_idx = 0; catname_idx != nonnull_catname_ct; ++catname_idx) {
            catname_entry_ptr = catname_entry_ptr->pheno_next;
            char* cur_name_start = name_storage_iter;
            name_storage_iter = strcpyax(name_storage_iter, catname_entry_ptr->str, '\0');
            *cur_name_ptrs++ = cur_name_start;
          }
          ++cat_pheno_idx;
        } else if (!is_qt) {
          pheno_cols_iter->data.cc = new_pheno_data_iter;
          ZeroWArr(nonmiss_vec_ct * kWordsPerVec, new_pheno_data_iter);
        } else {
          pheno_cols_iter->data.qt = R_CAST(double*, new_pheno_data_iter);
        }
        ++pheno_cols_iter;
      }
      while (pheno_info_reverse_ll) {
        const uint32_t sample_uidx = pheno_info_reverse_ll->sample_uidx;
        double* pheno_data = pheno_info_reverse_ll->phenodata;
        pheno_cols_iter = &(new_pheno_cols[old_pheno_ct]);
        for (uint32_t new_pheno_idx = 0; new_pheno_idx != new_pheno_ct; ++new_pheno_idx) {
          if (IsSet(categorical_phenos, new_pheno_idx)) {
            uint32_t cur_cat;
            memcpy(&cur_cat, &(pheno_data[new_pheno_idx]), sizeof(int32_t));
            pheno_cols_iter->data.cat[sample_uidx] = cur_cat;
            if (cur_cat) {
              SetBit(sample_uidx, pheno_cols_iter->nonmiss);
            }
          } else {
            double dxx = pheno_data[new_pheno_idx];
            // bugfix (6 May 2017): forgot to accept 0 as missing value for
            // case/control
            if (IsSet(quantitative_phenos, new_pheno_idx)) {
              if (dxx != missing_phenod) {
                SetBit(sample_uidx, pheno_cols_iter->nonmiss);
                pheno_cols_iter->data.qt[sample_uidx] = dxx;
              }
            } else {
              if (dxx == pheno_cased) {
                SetBit(sample_uidx, pheno_cols_iter->data.cc);
                SetBit(sample_uidx, pheno_cols_iter->nonmiss);
              } else if (dxx == pheno_ctrld) {
                SetBit(sample_uidx, pheno_cols_iter->nonmiss);
              }
            }
          }
          ++pheno_cols_iter;
        }
        pheno_info_reverse_ll = pheno_info_reverse_ll->next;
      }
    }
    *pheno_names_ptr = pheno_names;
    *max_pheno_name_blen_ptr = max_pheno_name_blen;
  }
  while (0) {
  LoadPhenos_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadPhenos_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pheno_fname, &pheno_txs);
    break;
  LoadPhenos_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
  LoadPhenos_ret_MALFORMED_INPUT_2:
    logerrputsb();
  LoadPhenos_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  LoadPhenos_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, pheno_fname);
    reterr = kPglRetMalformedInput;
    break;
  LoadPhenos_ret_INCOMPATIBLE_PHENOSTRS:
    logerrputs("(Case/control and quantitative phenotypes must all be numeric/'NA'.\nCategorical phenotypes cannot be 'NA'--use e.g. 'NONE' to represent missing\ncategorical values instead--or start with a number.)\n");
    reterr = kPglRetMalformedInput;
    break;
  LoadPhenos_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
  LoadPhenos_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
  LoadPhenos_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 LoadPhenos_ret_1:
  CleanupTextStream2(pheno_fname, &pheno_txs, &reterr);
  BigstackReset(bigstack_mark);
  if (reterr) {
    free_cond(pheno_names);
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

#ifdef __cplusplus
}  // namespace plink2
#endif
