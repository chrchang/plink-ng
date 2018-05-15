// This file is part of PLINK 2.00, copyright (C) 2005-2018 Shaun Purcell,
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


#include "plink2_compress_stream.h"
#include "plink2_data.h"
#include "plink2_misc.h"
#include "plink2_stats.h"

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

PglErr UpdateVarNames(const uintptr_t* variant_include, const uint32_t* variant_id_htable, const TwoColParams* params, uint32_t raw_variant_ct, uint32_t htable_size, char** variant_ids, uint32_t* max_variant_id_slen_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  ReadLineStream rls;
  PreinitRLstream(&rls);
  {
    const uint32_t orig_max_variant_id_slen = *max_variant_id_slen_ptr;
    uint32_t max_variant_id_slen = orig_max_variant_id_slen;
    char** variant_ids_copy;
    uintptr_t* already_seen;
    if (bigstack_alloc_cp(raw_variant_ct, &variant_ids_copy) ||
        bigstack_calloc_w(BitCtToWordCt(raw_variant_ct), &already_seen)) {
      goto UpdateVarNames_ret_NOMEM;
    }
    memcpy(variant_ids_copy, variant_ids, raw_variant_ct * sizeof(intptr_t));
    // This could be pointed at a file containing allele codes, so don't limit
    // line length to minimum value.
    // On the other hand, new variant IDs are allocated off the end of
    // bigstack, and that could result in a lot of memory pressure.
    const char* line_iter;
    reterr = SizeAndInitRLstreamRawK(params->fname, bigstack_left() / 4, &rls, &line_iter);
    if (reterr) {
      goto UpdateVarNames_ret_1;
    }
    reterr = RlsSkipK(params->skip_ct, &rls, &line_iter);
    if (reterr) {
      if (reterr == kPglRetEof) {
        snprintf(g_logbuf, kLogbufSize, "Error: Fewer lines than expected in %s.\n", params->fname);
        goto UpdateVarNames_ret_INCONSISTENT_INPUT_WW;
      }
      goto UpdateVarNames_ret_READ_RLSTREAM;
    }
    line_idx = params->skip_ct;
    // ok, this should be a parameter instead...
    const uint32_t* htable_dup_base = &(variant_id_htable[RoundUpPow2(htable_size, kInt32PerCacheline)]);
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
      reterr = RlsNextLstripK(&rls, &line_iter);
      if (reterr) {
        if (reterr == kPglRetEof) {
          reterr = kPglRetSuccess;
          break;
        }
        goto UpdateVarNames_ret_READ_RLSTREAM;
      }
      const char* linebuf_first_token = line_iter;
      char cc = *linebuf_first_token;
      if (IsEolnKns(cc) || (cc == skipchar)) {
        continue;
      }
      const char* colold_ptr;
      const char* colnew_ptr;
      if (colold_first) {
        colold_ptr = NextTokenMult0(linebuf_first_token, colmin);
        colnew_ptr = NextTokenMult(colold_ptr, coldiff);
        if (!colnew_ptr) {
          goto UpdateVarNames_ret_MISSING_TOKENS;
        }
      } else {
        colnew_ptr = NextTokenMult0(linebuf_first_token, colmin);
        colold_ptr = NextTokenMult(colnew_ptr, coldiff);
        if (!colold_ptr) {
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
      if (cur_llidx != UINT32_MAX) {
        // we could check if some copies have been filtered out after hash
        // table construction?
        snprintf(g_logbuf, kLogbufSize, "Error: --update-name variant ID '%s' appears multiple times in dataset.\n", cur_var_id);
        goto UpdateVarNames_ret_INCONSISTENT_INPUT_WW;
      }
      if (!IsSet(variant_include, variant_uidx)) {
        continue;
      }
      if (IsSet(already_seen, variant_uidx)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --update-name file.\n", cur_var_id);
        goto UpdateVarNames_ret_INCONSISTENT_INPUT_WW;
      }
      SetBit(variant_uidx, already_seen);
      ++hit_ct;
      const uint32_t colnew_slen = strlen_se(colnew_ptr);
      if (colnew_slen <= colold_slen) {
        const uint32_t colold_blen = colold_slen + 1;
        if (S_CAST(uintptr_t, alloc_end - alloc_base) < colold_blen) {
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
        if (S_CAST(uintptr_t, alloc_end - alloc_base) < colnew_blen) {
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
  UpdateVarNames_ret_READ_RLSTREAM:
    RLstreamErrPrint("--update-name file", &rls, &reterr);
    break;
  UpdateVarNames_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --update-name file has fewer tokens than expected.\n", line_idx);
  UpdateVarNames_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 UpdateVarNames_ret_1:
  CleanupRLstream(&rls);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr Plink1ClusterImport(const char* within_fname, const char* catpheno_name, const char* family_missing_catname, const uintptr_t* sample_include, const char* sample_ids, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t mwithin_val, PhenoCol** pheno_cols_ptr, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  ReadLineStream within_rls;
  PreinitRLstream(&within_rls);
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
      for (uint32_t pheno_idx = 0; pheno_idx < old_pheno_ct; ++pheno_idx) {
        if (!memcmp(catpheno_name, &(old_pheno_names[pheno_idx * old_max_pheno_name_blen]), catpheno_name_blen)) {
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
    if (pgl_malloc(new_pheno_names_byte_ct, &pheno_names)) {
      goto Plink1ClusterImport_ret_NOMEM;
    }
    if (old_pheno_names && (old_max_pheno_name_blen == new_max_pheno_name_blen)) {
      memcpy(pheno_names, old_pheno_names, old_pheno_ct * new_max_pheno_name_blen);
    } else {
      for (uint32_t pheno_idx = 0; pheno_idx < old_pheno_ct; ++pheno_idx) {
        strcpy(&(pheno_names[pheno_idx * new_max_pheno_name_blen]), &(old_pheno_names[pheno_idx * old_max_pheno_name_blen]));
      }
    }
    memcpy(&(pheno_names[old_pheno_ct * new_max_pheno_name_blen]), catpheno_name, catpheno_name_blen);
    free_cond(old_pheno_names);
    *pheno_names_ptr = pheno_names;

    PhenoCol* new_pheno_cols = S_CAST(PhenoCol*, realloc(*pheno_cols_ptr, new_pheno_ct * sizeof(PhenoCol)));
    if (!new_pheno_cols) {
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
      if (bigstack_alloc_w(raw_sample_ctaw, &cat_nm) ||
          bigstack_calloc_u32(raw_sample_ct, &cat_idxs)) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      memcpy(cat_nm, sample_include, raw_sample_ctaw * sizeof(intptr_t));
    }
    uint32_t* cat_htable;
    uint32_t cat_htable_size;
    if (HtableGoodSizeAlloc(sample_ct + 2, bigstack_left() / 4, &cat_htable, &cat_htable_size)) {
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
      if (bigstack_calloc_w(BitCtToWordCt(sample_ct), &already_seen) ||
          bigstack_calloc_u32(sample_ct, &sorted_cat_idxs) ||
          bigstack_alloc_c(max_sample_id_blen, &idbuf) ||
          bigstack_alloc_kcp(sample_ct + 2, &cur_cat_names)) {
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
      if (CopySortStrboxSubset(sample_include, sample_ids, sample_ct, max_sample_id_blen, 1, 0, 0, &sorted_idbox, &id_map)) {
        goto Plink1ClusterImport_ret_NOMEM;
      }
      char* line_iter;
      reterr = SizeAndInitRLstreamRaw(within_fname, bigstack_left() - (bigstack_left() / 4), &within_rls, &line_iter);
      if (reterr) {
        goto Plink1ClusterImport_ret_1;
      }
      char* cat_name_write_start = R_CAST(char*, g_bigstack_base);
      char* cat_name_iter = cat_name_write_start;
      char* cat_name_write_max = R_CAST(char*, g_bigstack_end);

      uint32_t nonnull_cat_ct = 0;
      uintptr_t miss_ct = 0;
      uintptr_t duplicate_ct = 0;
      while (1) {
        line_iter = AdvToDelim(line_iter, '\n');
      Plink1ClusterImport_LINE_ITER_ALREADY_ADVANCED:
        ++line_iter;
        ++line_idx;
        reterr = RlsPostlfNext(&within_rls, &line_iter);
        if (reterr) {
          if (reterr == kPglRetEof) {
            reterr = kPglRetSuccess;
            break;
          }
          goto Plink1ClusterImport_ret_READ_RLSTREAM;
        }
        line_iter = FirstNonTspace(line_iter);
        if (IsEolnKns(*line_iter)) {
          continue;
        }
        char* fid_start = line_iter;
        char* fid_end = CurTokenEnd(fid_start);
        char* iid_start = FirstNonTspace(fid_end);
        if (IsEolnKns(*iid_start)) {
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
        if (!main_token_start) {
          goto Plink1ClusterImport_ret_MISSING_TOKENS;
        }
        char* main_token_end = CurTokenEnd(main_token_start);
        line_iter = AdvToDelim(line_iter, '\n');
        *main_token_end = '\0';
        const uint32_t main_token_slen = main_token_end - main_token_start;
        if (main_token_slen > kMaxIdSlen) {
          logerrputs("Error: Category names are limited to " MAX_ID_SLEN_STR " characters.\n");
          goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
        }
        uint32_t hashval = Hashceil(main_token_start, main_token_slen, cat_htable_size);
        const uint32_t main_token_blen = main_token_slen + 1;
        uint32_t cur_htable_entry;
        while (1) {
          cur_htable_entry = cat_htable[hashval];
          if (cur_htable_entry == UINT32_MAX) {
            if (main_token_blen > S_CAST(uintptr_t, cat_name_write_max - cat_name_iter)) {
              goto Plink1ClusterImport_ret_NOMEM;
            }
            char* cat_name_start = cat_name_iter;
            cat_name_iter = memcpya(cat_name_iter, main_token_start, main_token_blen);
            cur_cat_names[++nonnull_cat_ct] = cat_name_start;
            cur_htable_entry = nonnull_cat_ct;
            cat_htable[hashval] = cur_htable_entry;
            break;
          }
          if (!memcmp(main_token_start, cur_cat_names[cur_htable_entry], main_token_blen)) {
            break;
          }
          if (++hashval == cat_htable_size) {
            hashval = 0;
          }
        }
        // permit duplicates if category is identical
        if (IsSet(already_seen, lb_idx)) {
          const uint32_t existing_cat_idx = sorted_cat_idxs[lb_idx];
          if (existing_cat_idx != cur_htable_entry) {
            idbuf[fid_slen] = ' ';
            logpreprintfww("Error: Duplicate sample ID '%s' with conflicting category assignments in --within file.\n", idbuf);
            goto Plink1ClusterImport_ret_MALFORMED_INPUT_2;
          }
          ++duplicate_ct;
        } else {
          SetBit(lb_idx, already_seen);
          for (; lb_idx < ub_idx; ++lb_idx) {
            sorted_cat_idxs[lb_idx] = cur_htable_entry;
          }
        }
        goto Plink1ClusterImport_LINE_ITER_ALREADY_ADVANCED;
      }
      if (!nonnull_cat_ct) {
        logerrputs("Error: All --within categories are null.\n");
        goto Plink1ClusterImport_ret_INCONSISTENT_INPUT_WW;
      }
      double dxx;
      const uint32_t prepend_c = (ScanadvDouble(cur_cat_names[1], &dxx) != nullptr);
      if (prepend_c) {
        for (uint32_t catname_idx = 2; catname_idx <= nonnull_cat_ct; ++catname_idx) {
          if (!ScanadvDouble(cur_cat_names[catname_idx], &dxx)) {
            logerrputs("Error: Either all non-null --within categories must be numeric, or none can be.\n");
            goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
          }
        }
        logputs("Note: Prepending 'C' to all --within category names.\n");
      } else {
        for (uint32_t catname_idx = 2; catname_idx <= nonnull_cat_ct; ++catname_idx) {
          if (ScanadvDouble(cur_cat_names[catname_idx], &dxx)) {
            logerrputs("Error: Either all non-null --within categories must be numeric, or none can be.\n");
            goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
          }
        }
      }
      // see end of e.g. LoadPsam()
      const uintptr_t catname_vec_ct = WordCtToVecCt(nonnull_cat_ct + 1);
      const uintptr_t total_catname_blen = (prepend_c * nonnull_cat_ct) + S_CAST(uintptr_t, cat_name_iter - cat_name_write_start);
      const uintptr_t catname_storage_vec_ct = DivUp(total_catname_blen, kBytesPerVec);
      if (vecaligned_malloc((raw_sample_ctaw * kWordsPerVec + data_vec_ct + catname_vec_ct + catname_storage_vec_ct) * kBytesPerVec, &(new_pheno_cols[old_pheno_ct].nonmiss))) {
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

      for (uint32_t sorted_sample_idx = 0; sorted_sample_idx < sample_ct; ++sorted_sample_idx) {
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
      for (uint32_t uii = 0; uii < nonnull_cat_ct; ++uii) {
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
        } else if ((missing_catname_slen != family_missing_catname_slen) || memcmp(family_missing_catname, missing_catname, missing_catname_slen)) {
          if (++family_missing_catname_hval == cat_htable_size) {
            family_missing_catname_hval = 0;
          }
          cat_htable[family_missing_catname_hval] = UINT32_MAXM1;
        }
      }
      // guaranteed to have enough space
      uint32_t* cat_idx_m1_to_first_sample_uidx = R_CAST(uint32_t*, g_bigstack_base);
      uint32_t sample_uidx = 0;
      uint32_t nonnull_cat_ct = 0;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
        MovU32To1Bit(sample_include, &sample_uidx);
        const char* cur_fid = &(sample_ids[sample_uidx * max_sample_id_blen]);
        const char* cur_fid_end = AdvToDelim(cur_fid, '\t');
        const uint32_t slen = cur_fid_end - cur_fid;
        uint32_t hashval = Hashceil(cur_fid, slen, cat_htable_size);
        const uint32_t blen = slen + 1;
        while (1) {
          const uint32_t cur_htable_entry = cat_htable[hashval];
          if (cur_htable_entry >= 0xfffffffdU) {
            if (cur_htable_entry == UINT32_MAX) {
              cat_htable[hashval] = sample_uidx;
              total_catname_blen += blen;
              cat_idx_m1_to_first_sample_uidx[nonnull_cat_ct] = sample_uidx;
              cat_idxs[sample_uidx] = ++nonnull_cat_ct;
              break;
            } else if (cur_htable_entry == UINT32_MAXM1) {
              if ((slen == family_missing_catname_slen) && (!memcmp(cur_fid, family_missing_catname, family_missing_catname_slen))) {
                ClearBit(sample_uidx, cat_nm);
                cat_idxs[sample_uidx] = 0;
                break;
              }
            } else {
              if ((slen == missing_catname_slen) && (!memcmp(cur_fid, missing_catname, missing_catname_slen))) {
                ClearBit(sample_uidx, cat_nm);
                cat_idxs[sample_uidx] = 0;
                break;
              }
            }
          } else {
            if (!memcmp(cur_fid, &(sample_ids[cur_htable_entry * max_sample_id_blen]), blen)) {
              cat_idxs[sample_uidx] = cat_idxs[cur_htable_entry];
              break;
            }
          }
          if (++hashval == cat_htable_size) {
            hashval = 0;
          }
        }
      }
      if (!nonnull_cat_ct) {
        logerrputs("Error: All --family FIDs are null.\n");
        goto Plink1ClusterImport_ret_INCONSISTENT_INPUT_WW;
      }
      // add 'P' prefixes?
      double dxx;
      const uint32_t prepend_c = (ScanadvDouble(&(sample_ids[cat_idx_m1_to_first_sample_uidx[0] * max_sample_id_blen]), &dxx) != nullptr);
      if (prepend_c) {
        for (uint32_t uii = 1; uii < nonnull_cat_ct; ++uii) {
          if (!ScanadvDouble(&(sample_ids[cat_idx_m1_to_first_sample_uidx[uii] * max_sample_id_blen]), &dxx)) {
            logerrputs("Error: Either all non-null --family FIDs must be numeric, or none can be.\n");
            goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
          }
        }
        logputs("Note: Prepending 'C' to all --family category names.\n");
        total_catname_blen += nonnull_cat_ct;
      } else {
        for (uint32_t uii = 1; uii < nonnull_cat_ct; ++uii) {
          if (ScanadvDouble(&(sample_ids[cat_idx_m1_to_first_sample_uidx[uii] * max_sample_id_blen]), &dxx)) {
            logerrputs("Error: Either all non-null --family FIDs must be numeric, or none can be.\n");
            goto Plink1ClusterImport_ret_INCONSISTENT_INPUT;
          }
        }
      }
      // see end of e.g. LoadPsam()
      const uintptr_t catname_vec_ct = WordCtToVecCt(nonnull_cat_ct + 1);
      const uintptr_t catname_storage_vec_ct = DivUp(total_catname_blen, kBytesPerVec);
      if (vecaligned_malloc((raw_sample_ctaw * kWordsPerVec + data_vec_ct + catname_vec_ct + catname_storage_vec_ct) * kBytesPerVec, &(new_pheno_cols[old_pheno_ct].nonmiss))) {
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
      for (uint32_t uii = 0; uii < nonnull_cat_ct; ++uii) {
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
  Plink1ClusterImport_ret_READ_RLSTREAM:
    RLstreamErrPrint("--within file", &within_rls, &reterr);
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
  CleanupRLstream(&within_rls);
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

PglErr UpdateSampleSexes(const uintptr_t* sample_include, const SampleIdInfo* siip, const UpdateSexInfo* update_sex_info_ptr, uint32_t raw_sample_ct, uintptr_t sample_ct, uintptr_t* sex_nm, uintptr_t* sex_male) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  ReadLineStream rls;
  PreinitRLstream(&rls);
  {
    if (!sample_ct) {
      goto UpdateSampleSexes_ret_1;
    }
    // permit very long lines since this can be pointed at .ped files
    char* line_iter;
    reterr = SizeAndInitRLstreamRaw(update_sex_info_ptr->fname, bigstack_left() - (bigstack_left() / 4), &rls, &line_iter);
    if (reterr) {
      goto UpdateSampleSexes_ret_1;
    }

    // (Much of this boilerplate is shared with e.g. KeepFcol(); it probably
    // belongs in its own function.)
    char* linebuf_first_token;
    XidMode xid_mode;
    reterr = LoadXidHeader("update-sex", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeaderFixedWidth : kfXidHeaderFixedWidthIgnoreSid, &line_iter, &line_idx, &linebuf_first_token, &rls, &xid_mode);
    if (reterr) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: Empty --update-sex file.\n");
        goto UpdateSampleSexes_ret_MALFORMED_INPUT;
      }
      goto UpdateSampleSexes_ret_READ_RLSTREAM;
    }
    char* header_sample_id_end = line_iter;
    const uint32_t id_col_ct = GetXidColCt(xid_mode);
    uint32_t col_num = update_sex_info_ptr->col_num;
    uint32_t postid_col_idx = 0;
    if ((*linebuf_first_token == '#') && (!col_num)) {
      // search for 'SEX' column (any capitalization)
      const char* token_end = header_sample_id_end;
      while (1) {
        ++postid_col_idx;
        const char* linebuf_iter = FirstNonTspace(token_end);
        if (IsEolnKns(*linebuf_iter)) {
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
        if (id_col_ct == 3) {
          logerrputs("Error: You must use 'col-num=' to specify the position of the sex column in the\n--update-sex file.\n");
          goto UpdateSampleSexes_ret_MALFORMED_INPUT;
        }
        col_num = 3;
      }
      if (id_col_ct >= col_num) {
        logerrputs("Error: --update-sex 'col-num=' parameter too small (it refers to a sample ID\ncolumn).\n");
        goto UpdateSampleSexes_ret_MALFORMED_INPUT;
      }
      postid_col_idx = col_num - id_col_ct;
    }
    if (*linebuf_first_token == '#') {
      // advance to next line
      *linebuf_first_token = '\0';
    }

    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uint32_t* xid_map = nullptr;
    char* sorted_xidbox = nullptr;
    const uint32_t allow_dups = siip->sids && (!(xid_mode & kfXidModeFlagSid));
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, allow_dups, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (reterr) {
      goto UpdateSampleSexes_ret_1;
    }
    uintptr_t* already_seen;
    char* idbuf;
    if (bigstack_calloc_w(raw_sample_ctl, &already_seen) ||
        bigstack_alloc_c(max_xid_blen, &idbuf)) {
      goto UpdateSampleSexes_ret_NOMEM;
    }

    const uint32_t male0 = (update_sex_info_ptr->flags / kfUpdateSexMale0) & 1;
    uint32_t hit_ct = 0;
    uintptr_t miss_ct = 0;
    uintptr_t duplicate_ct = 0;
    while (1) {
      if (!IsEolnKns(*linebuf_first_token)) {
        const char* linebuf_iter = linebuf_first_token;
        uint32_t xid_idx_start;
        uint32_t xid_idx_end;
        if (!SortedXidboxReadMultifind(sorted_xidbox, max_xid_blen, sample_ct, 0, xid_mode, &linebuf_iter, &xid_idx_start, &xid_idx_end, idbuf)) {
          const char* sex_start = NextTokenMult(linebuf_iter, postid_col_idx);
          if (!sex_start) {
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
            } else if ((!male0) && (sexval != 30)) {
              // allow 'N' = missing to make 1/2/NA work
              // don't permit 'n' for now
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid sex value on line %" PRIuPTR " of --update-sex file. (Acceptable values: 1/M/m = male, 2/F/f = female, 0/N = missing.)\n", line_idx);
              goto UpdateSampleSexes_ret_MALFORMED_INPUT_WW;
            } else {
              // with 'male0', everything else is treated as missing
              sexval = 0;
            }
          } else if (male0) {
            if (sexval == 2) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid sex value on line %" PRIuPTR " of --update-sex file. ('2' is prohibited when the 'male0' modifier is present.)\n", line_idx);
              goto UpdateSampleSexes_ret_MALFORMED_INPUT_WW;
            }
            ++sexval;
          }
          uint32_t sample_uidx = xid_map[xid_idx_start];
          if (IsSet(already_seen, sample_uidx)) {
            // permit duplicates iff sex value is identical
            const uint32_t old_sexval = IsSet(sex_nm, sample_uidx) * (2 - IsSet(sex_male, sample_uidx));
            if (sexval != old_sexval) {
              snprintf(g_logbuf, kLogbufSize, "Error: Sample ID on line %" PRIuPTR " of --update-sex file duplicates one earlier in the file, and sex values don't match.\n", line_idx);
              goto UpdateSampleSexes_ret_MALFORMED_INPUT_WW;
            }
            ++duplicate_ct;
            continue;
          }
          SetBit(sample_uidx, already_seen);
          while (1) {
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
            ++hit_ct;
            if (++xid_idx_start == xid_idx_end) {
              break;
            }
            sample_uidx = xid_map[xid_idx_start];
          }
        } else if (!linebuf_iter) {
          goto UpdateSampleSexes_ret_MISSING_TOKENS;
        }
      }
      ++line_idx;
      reterr = RlsNextLstrip(&rls, &line_iter);
      if (reterr) {
        if (reterr == kPglRetEof) {
          reterr = kPglRetSuccess;
          break;
        }
        goto UpdateSampleSexes_ret_READ_RLSTREAM;
      }
      linebuf_first_token = line_iter;
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
  UpdateSampleSexes_ret_READ_RLSTREAM:
    RLstreamErrPrint("--update-sex file", &rls, &reterr);
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
  CleanupRLstream(&rls);
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
    const uint32_t omit_last = (pheno_transform_flags / kfPhenoTransformSplitCatOmitLast) & 1;
    uint32_t qt_12 = 0;
    uint32_t at_least_one_cat_pheno_processed = 0;
    for (uint32_t is_covar = 0; is_covar < 2; ++is_covar) {
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
      if (bigstack_calloc_w(old_pheno_ctl, &phenos_to_split)) {
        goto SplitCatPheno_ret_NOMEM;
      }
      if (!split_cat_phenonames_flattened) {
        for (uint32_t pheno_idx = 0; pheno_idx < old_pheno_ct; ++pheno_idx) {
          const PhenoCol* cur_pheno_col = &(old_pheno_cols[pheno_idx]);
          if (cur_pheno_col->type_code == kPhenoDtypeCat) {
            if (strchr(&(old_pheno_names[pheno_idx * old_max_pheno_name_blen]), '=')) {
              logerrputs("Error: --split-cat-pheno cannot be used on phenotypes containing the '='\ncharacter.\n");
              goto SplitCatPheno_ret_INCONSISTENT_INPUT;
            }
            SetBit(pheno_idx, phenos_to_split);
          }
        }
      } else {
        uint32_t* id_htable;
        uint32_t id_htable_size;
        if (HtableGoodSizeAlloc(old_pheno_ct, bigstack_left(), &id_htable, &id_htable_size)) {
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
              if (old_pheno_cols[pheno_idx].type_code != kPhenoDtypeCat) {
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
      uint32_t* observed_cat_cts;  // excludes null, excludes last if omit-last
      uintptr_t** observed_cats;
      uintptr_t* sample_include_intersect;
      if (bigstack_alloc_u32(split_pheno_ct, &observed_cat_cts) ||
          bigstack_alloc_wp(split_pheno_ct, &observed_cats) ||
          bigstack_alloc_w(raw_sample_ctl, &sample_include_intersect)) {
        goto SplitCatPheno_ret_NOMEM;
      }
      uintptr_t create_pheno_ct = 0;
      uint32_t split_pheno_uidx = 0;
      uint32_t max_cat_uidx_p1 = 0;
      for (uint32_t split_pheno_idx = 0; split_pheno_idx < split_pheno_ct; ++split_pheno_idx, ++split_pheno_uidx) {
        MovU32To1Bit(phenos_to_split, &split_pheno_uidx);
        const PhenoCol* cur_pheno_col = &(old_pheno_cols[split_pheno_uidx]);
        BitvecAndCopy(sample_include, cur_pheno_col->nonmiss, raw_sample_ctl, sample_include_intersect);
        const uint32_t cur_cat_ct = cur_pheno_col->nonnull_category_ct + 1;
        const uint32_t cur_cat_ctl = BitCtToWordCt(cur_cat_ct);
        uintptr_t* cur_observed_cats;
        if (bigstack_alloc_w(cur_cat_ctl, &cur_observed_cats)) {
          goto SplitCatPheno_ret_NOMEM;
        }
        observed_cats[split_pheno_idx] = cur_observed_cats;
        const uint32_t cur_nmiss_ct = PopcountWords(sample_include_intersect, raw_sample_ctl);
        uint32_t cur_observed_cat_ct = IdentifyRemainingCats(sample_include_intersect, cur_pheno_col, cur_nmiss_ct, cur_observed_cats);
        if (cur_observed_cat_ct > omit_last) {
          cur_observed_cat_ct -= omit_last;
          // old phenotype name, '=' character, null terminator
          const uintptr_t blen_base = strlen(&(old_pheno_names[split_pheno_uidx * old_max_pheno_name_blen])) + 2;
          const char* const* cat_names = cur_pheno_col->category_names;
          uint32_t cat_uidx = 0;
          for (uint32_t cat_idx = 0; cat_idx < cur_observed_cat_ct; ++cat_idx, ++cat_uidx) {
            MovU32To1Bit(cur_observed_cats, &cat_uidx);
            const char* cur_cat_name = cat_names[cat_uidx];
            const uint32_t cur_slen = strlen(cur_cat_name);
            if (memchr(cur_cat_name, '=', cur_slen)) {
              logerrputs("Error: --split-cat-pheno category names may not contain the '=' character.\n");
              goto SplitCatPheno_ret_INCONSISTENT_INPUT;
            }
            const uintptr_t total_blen = cur_slen + blen_base;
            if (total_blen > new_max_pheno_name_blen) {
              new_max_pheno_name_blen = total_blen;
            }
          }
          if (cat_uidx > max_cat_uidx_p1) {
            max_cat_uidx_p1 = cat_uidx;
          }
        } else {
          cur_observed_cat_ct = 0;
        }
        observed_cat_cts[split_pheno_idx] = cur_observed_cat_ct;
        create_pheno_ct += cur_observed_cat_ct;
      }
      if (new_max_pheno_name_blen > kMaxIdBlen) {
        logerrputs("Error: New --split-cat-pheno phenotype name too long.  Shorten your phenotype\nor your category names.\n");
        goto SplitCatPheno_ret_INCONSISTENT_INPUT;
      }
      const uint32_t copy_pheno_ct = old_pheno_ct - split_pheno_ct;
      // before new_pheno_ct variable definition due to potential integer
      // overflow
      if (create_pheno_ct + copy_pheno_ct > kMaxPhenoCt) {
        logerrputs("Error: --split-cat-pheno would create too many phenotypes (" PROG_NAME_STR " is limited to\n" MAX_PHENO_CT_STR ").\n");
        goto SplitCatPheno_ret_INCONSISTENT_INPUT;
      }
      const uint32_t new_pheno_ct = create_pheno_ct + copy_pheno_ct;
      uintptr_t** write_data_ptrs;
      if (bigstack_alloc_wp(max_cat_uidx_p1, &write_data_ptrs)) {
        goto SplitCatPheno_ret_NOMEM;
      }
      const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
      uint32_t new_data_word_ct = raw_sample_ctaw;
      if (is_covar) {
        new_data_word_ct = DblCtToVecCt(raw_sample_ct) * kWordsPerVec;
      }
      uintptr_t* omit_last_dummy = nullptr;
      if (omit_last) {
        if (bigstack_alloc_w(new_data_word_ct, &omit_last_dummy)) {
          goto SplitCatPheno_ret_NOMEM;
        }
      }

      // second pass: allocate memory and actually create the new phenotypes
      char* new_pheno_names;
      if (pgl_malloc(new_pheno_ct * new_max_pheno_name_blen, &new_pheno_names)) {
        goto SplitCatPheno_ret_NOMEM;
      }
      doomed_pheno_names = old_pheno_names;
      *xpheno_names_ptr = new_pheno_names;
      PhenoCol* new_pheno_cols;
      if (pgl_malloc(new_pheno_ct * sizeof(PhenoCol), &new_pheno_cols)) {
        goto SplitCatPheno_ret_NOMEM;
      }
      doomed_pheno_cols = old_pheno_cols;
      doomed_pheno_ct = old_pheno_ct;
      *xpheno_cols_ptr = new_pheno_cols;
      *xpheno_ct_ptr = new_pheno_ct;
      *max_xpheno_name_blen_ptr = new_max_pheno_name_blen;
      uint32_t pheno_read_idx = 0;
      for (uint32_t pheno_write_idx = 0; pheno_write_idx < copy_pheno_ct; ++pheno_write_idx, ++pheno_read_idx) {
        MovU32To0Bit(phenos_to_split, &pheno_read_idx);

        // manually move this data
        memcpy(&(new_pheno_cols[pheno_write_idx]), &(doomed_pheno_cols[pheno_read_idx]), sizeof(PhenoCol));
        doomed_pheno_cols[pheno_read_idx].nonmiss = nullptr;

        strcpy(&(new_pheno_names[pheno_write_idx * new_max_pheno_name_blen]), &(old_pheno_names[pheno_read_idx * old_max_pheno_name_blen]));
      }
      for (uint32_t pheno_write_idx = copy_pheno_ct; pheno_write_idx < new_pheno_ct; ++pheno_write_idx) {
        new_pheno_cols[pheno_write_idx].nonmiss = nullptr;
      }

      const uintptr_t new_pheno_bytes_req = (raw_sample_ctaw + new_data_word_ct) * sizeof(intptr_t);
      uint32_t pheno_write_idx = copy_pheno_ct;
      pheno_read_idx = 0;
      for (uint32_t split_pheno_idx = 0; split_pheno_idx < split_pheno_ct; ++split_pheno_idx, ++pheno_read_idx) {
        MovU32To1Bit(phenos_to_split, &pheno_read_idx);
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
        uint32_t orig_cat_idx = 1;
        for (uint32_t uii = 0; uii < cur_pheno_write_ct; ++uii, ++orig_cat_idx, ++pheno_write_idx) {
          MovU32To1Bit(cur_observed_cats, &orig_cat_idx);
          uintptr_t* new_pheno_data_iter;
          if (vecaligned_malloc(new_pheno_bytes_req, &new_pheno_data_iter)) {
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
              for (uint32_t ujj = 0; ujj < raw_sample_ct; ++ujj) {
                pheno_qt[ujj] = 1.0;
              }
            } else {
              ZeroDArr(raw_sample_ct, pheno_qt);
            }
          }
        }
        if (omit_last) {
          MovU32To1Bit(cur_observed_cats, &orig_cat_idx);
          write_data_ptrs[orig_cat_idx] = omit_last_dummy;
        }

        const uint32_t cur_nmiss_ct = PopcountWords(sample_include_intersect, raw_sample_ctl);
        const uint32_t* cur_cats = old_pheno_col->data.cat;
        uint32_t sample_uidx = 0;
        if (!is_covar) {
          for (uint32_t sample_idx = 0; sample_idx < cur_nmiss_ct; ++sample_idx, ++sample_uidx) {
            MovU32To1Bit(sample_include_intersect, &sample_uidx);
            SetBit(sample_uidx, write_data_ptrs[cur_cats[sample_uidx]]);
          }
        } else {
          double** write_qt_ptrs = R_CAST(double**, write_data_ptrs);
          const double write_val = u31tod(1 + qt_12);
          for (uint32_t sample_idx = 0; sample_idx < cur_nmiss_ct; ++sample_idx, ++sample_uidx) {
            MovU32To1Bit(sample_include_intersect, &sample_uidx);
            write_qt_ptrs[cur_cats[sample_uidx]][sample_uidx] = write_val;
          }
        }
      }

      // if any preexisting phenotype names contain a single copy of the '='
      // character, verify that we didn't create any duplicate IDs
      for (uint32_t pheno_idx = 0; pheno_idx < copy_pheno_ct; ++pheno_idx) {
        const char* first_eq = strchr(&(new_pheno_names[pheno_idx * new_max_pheno_name_blen]), '=');
        if (first_eq && (!strchr(&(first_eq[1]), '='))) {
          uint32_t* id_htable;
          uint32_t id_htable_size;
          if (HtableGoodSizeAlloc(new_pheno_ct, bigstack_left(), &id_htable, &id_htable_size)) {
            goto SplitCatPheno_ret_NOMEM;
          }
          uint32_t duplicate_idx = PopulateStrboxHtable(new_pheno_names, new_pheno_ct, new_max_pheno_name_blen, id_htable_size, id_htable);
          if (duplicate_idx) {
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
    if (bigstack_calloc_w(pheno_ctl, &phenos_to_transform)) {
      goto PhenoVarianceStandardize_ret_NOMEM;
    }
    if (!vstd_flattened) {
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
        const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
        if (cur_pheno_col->type_code == kPhenoDtypeQt) {
          SetBit(pheno_idx, phenos_to_transform);
        }
      }
    } else {
      uint32_t* id_htable;
      uint32_t id_htable_size;
      if (HtableGoodSizeAlloc(pheno_ct, bigstack_left(), &id_htable, &id_htable_size)) {
        goto PhenoVarianceStandardize_ret_NOMEM;
      }
      PopulateStrboxHtable(pheno_names, pheno_ct, max_pheno_name_blen, id_htable_size, id_htable);
      const char* vstd_phenonames_iter = vstd_flattened;
      do {
        const uint32_t cur_phenoname_slen = strlen(vstd_phenonames_iter);
        if (cur_phenoname_slen < max_pheno_name_blen) {
          uint32_t pheno_idx = StrboxHtableFind(vstd_phenonames_iter, pheno_names, id_htable, max_pheno_name_blen, cur_phenoname_slen, id_htable_size);
          if (pheno_idx != UINT32_MAX) {
            if (pheno_cols[pheno_idx].type_code != kPhenoDtypeQt) {
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
    if (bigstack_alloc_d(raw_sample_ct, &shifted_pheno_qt)) {
      goto PhenoVarianceStandardize_ret_NOMEM;
    }
    const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
    uint32_t pheno_uidx = 0;
    for (uint32_t pheno_transform_idx = 0; pheno_transform_idx < pheno_transform_ct; ++pheno_transform_idx, ++pheno_uidx) {
      MovU32To1Bit(phenos_to_transform, &pheno_uidx);
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
      uint32_t sample_uidx = first_sample_uidx;
      shifted_pheno_qt[sample_uidx] = 0.0;
      const double pheno_shift = pheno_qt[sample_uidx++];
      for (uint32_t sample_idx = 1; sample_idx < cur_sample_ct; ++sample_idx, ++sample_uidx) {
        MovU32To1Bit(pheno_nm, &sample_uidx);
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
      sample_uidx = first_sample_uidx;
      for (uint32_t sample_idx = 0; sample_idx < cur_sample_ct; ++sample_idx, ++sample_uidx) {
        MovU32To1Bit(pheno_nm, &sample_uidx);
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
    if (bigstack_calloc_w(pheno_ctl, &phenos_to_transform)) {
      goto PhenoQuantileNormalize_ret_NOMEM;
    }
    if (!quantnorm_flattened) {
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
        const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
        if (cur_pheno_col->type_code == kPhenoDtypeQt) {
          SetBit(pheno_idx, phenos_to_transform);
        }
      }
    } else {
      uint32_t* id_htable;
      uint32_t id_htable_size;
      if (HtableGoodSizeAlloc(pheno_ct, bigstack_left(), &id_htable, &id_htable_size)) {
        goto PhenoQuantileNormalize_ret_NOMEM;
      }
      PopulateStrboxHtable(pheno_names, pheno_ct, max_pheno_name_blen, id_htable_size, id_htable);
      const char* quantnorm_phenonames_iter = quantnorm_flattened;
      do {
        const uint32_t cur_phenoname_slen = strlen(quantnorm_phenonames_iter);
        if (cur_phenoname_slen < max_pheno_name_blen) {
          uint32_t pheno_idx = StrboxHtableFind(quantnorm_phenonames_iter, pheno_names, id_htable, max_pheno_name_blen, cur_phenoname_slen, id_htable_size);
          if (pheno_idx != UINT32_MAX) {
            if (pheno_cols[pheno_idx].type_code != kPhenoDtypeQt) {
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
    if (BIGSTACK_ALLOC_X(DblIndex, sample_ct, &tagged_raw_pheno_vals)) {
      goto PhenoQuantileNormalize_ret_NOMEM;
    }
    const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
    uint32_t pheno_uidx = 0;
    for (uint32_t pheno_transform_idx = 0; pheno_transform_idx < pheno_transform_ct; ++pheno_transform_idx, ++pheno_uidx) {
      MovU32To1Bit(phenos_to_transform, &pheno_uidx);
      PhenoCol* cur_pheno_col = &(pheno_cols[pheno_uidx]);
      uintptr_t* pheno_nm = cur_pheno_col->nonmiss;
      BitvecAnd(sample_include, raw_sample_ctaw, pheno_nm);
      const uint32_t cur_sample_ct = PopcountWords(pheno_nm, raw_sample_ctaw);
      if (!cur_sample_ct) {
        continue;
      }
      double* pheno_qt = cur_pheno_col->data.qt;
      uint32_t sample_uidx = 0;
      for (uint32_t sample_idx = 0; sample_idx < cur_sample_ct; ++sample_idx, ++sample_uidx) {
        // bugfix (1 Sep 2017): this needs to iterate over pheno_nm, not
        // sample_include
        MovU32To1Bit(pheno_nm, &sample_uidx);
        tagged_raw_pheno_vals[sample_idx].dxx = pheno_qt[sample_uidx];
        tagged_raw_pheno_vals[sample_idx].uii = sample_uidx;
      }
      STD_SORT(cur_sample_ct, double_cmp, tagged_raw_pheno_vals);
      const double sample_ct_x2_recip = 1.0 / S_CAST(double, 2 * cur_sample_ct);
      for (uint32_t sample_idx_start = 0; sample_idx_start < cur_sample_ct;) {
        const double cur_raw_pheno = tagged_raw_pheno_vals[sample_idx_start].dxx;
        uint32_t sample_idx_end = sample_idx_start + 1;
        for (; sample_idx_end < cur_sample_ct; ++sample_idx_end) {
          if (tagged_raw_pheno_vals[sample_idx_end].dxx != cur_raw_pheno) {
            break;
          }
        }
        const double cur_zscore = QuantileToZscore(S_CAST(double, sample_idx_start + sample_idx_end) * sample_ct_x2_recip);
        for (; sample_idx_start < sample_idx_end; ++sample_idx_start) {
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


PglErr ProcessBoundaryToken(const char* tok_start, const char* tok_end, const char* token_source_str, uint32_t max_boundary_ct, PglErr err_type, double* prev_boundary_ptr, uint32_t* boundary_ct_ptr, double** freq_bounds_ptr, uint64_t** dosage_bounds_ptr) {
  double cur_boundary;
  const char* scan_end = ScanadvDouble(tok_start, &cur_boundary);
  if ((!scan_end) || (scan_end != tok_end)) {
    logerrprintf("Error: Invalid token in %s.\n", token_source_str);
    return err_type;
  }
  if (cur_boundary <= (*prev_boundary_ptr)) {
    logerrputs("Error: Invalid bin boundary sequence (must be strictly increasing, and start\nwith a positive number).\n");
    return err_type;
  }
  uint32_t boundary_ct = *boundary_ct_ptr;
  if (boundary_ct == max_boundary_ct) {
#ifdef __LP64__
    if (max_boundary_ct == 0x40000000) {
      logerrputs("Error: Too many bin boundaries.\n");
      return err_type;
    }
#endif
    return kPglRetNomem;
  }
  if (freq_bounds_ptr) {
    if (cur_boundary > 1.0) {
      logerrputs("Error: --freq bin boundary too large (must be <= 1).\n");
      return err_type;
    }
    // strictly-greater-than comparisons
    (*freq_bounds_ptr)[boundary_ct] = cur_boundary * (1 - kSmallEpsilon);
  } else {
    // max 2^31 - 3 variants
    if (cur_boundary > 4294967290.0) {
      logerrputs("Error: --freq counts bin boundary too large.\n");
      return err_type;
    }
    // due to the use of strictly-greater-than for comparison, we round
    // exact multiples of 1/32768 down
    const int64_t int_part = S_CAST(int64_t, cur_boundary);
    const double cur_boundary_frac_part = cur_boundary - int_part;
    const int64_t int_part_scaled = int_part * kDosageMax;
    if (cur_boundary_frac_part == 0.0) {
      (*dosage_bounds_ptr)[boundary_ct] = int_part_scaled - 1;
    } else {
      (*dosage_bounds_ptr)[boundary_ct] = int_part_scaled + S_CAST(int64_t, cur_boundary_frac_part * (kDosageMax * (1 - kSmallEpsilon)));
    }
  }
  *prev_boundary_ptr = cur_boundary;
  *boundary_ct_ptr = boundary_ct + 1;
  return kPglRetSuccess;
}

PglErr InitHistogramFromFileOrCommalist(const char* binstr, uint32_t is_fname, double** freq_bounds_ptr, uint64_t** dosage_bounds_ptr, uint32_t* boundary_ct_ptr, uint32_t** histogram_ptr) {
  GzTokenStream gts;
  PreinitGzTokenStream(&gts);
  uint32_t max_boundary_ct = 0;
  PglErr reterr = kPglRetSuccess;
  {
    uintptr_t ulii = RoundDownPow2(bigstack_left(), kCacheline);
    if (ulii < 2 * kCacheline) {
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
      *dosage_bounds_ptr = R_CAST(uint64_t*, g_bigstack_base);
    }
    uint32_t boundary_ct = 0;
    double prev_boundary = 0.0;
    if (is_fname) {
      // we want to accept >100000 numbers on a single line.  this will reject
      // "000...{a million more zeroes}...1"; pretty sure that's okay.
      reterr = InitGzTokenStream(binstr, &gts, g_textbuf);
      if (reterr) {
        goto InitHistogramFromFileOrCommalist_ret_1;
      }
      uint32_t token_slen;
      while (1) {
        const char* token_start = AdvanceGzTokenStream(&gts, &token_slen);
        if (!token_start) {
          break;
        }
        reterr = ProcessBoundaryToken(token_start, &(token_start[token_slen]), binstr, max_boundary_ct, kPglRetMalformedInput, &prev_boundary, &boundary_ct, freq_bounds_ptr, dosage_bounds_ptr);
        if (reterr) {
          goto InitHistogramFromFileOrCommalist_ret_1;
        }
      }
      if (token_slen) {
        if (token_slen == UINT32_MAX) {
          snprintf(g_logbuf, kLogbufSize, "Error: Excessively long token in %s.\n", binstr);
          goto InitHistogramFromFileOrCommalist_ret_MALFORMED_INPUT_2;
        } else {
          goto InitHistogramFromFileOrCommalist_ret_READ_FAIL;
        }
      }
    } else {
      const char* binstr_iter = binstr;
      while (1) {
        const char* tok_end = strchrnul(binstr_iter, ',');
        reterr = ProcessBoundaryToken(binstr_iter, tok_end, "--freq {ref,alt1}bins= list", max_boundary_ct, kPglRetInvalidCmdline, &prev_boundary, &boundary_ct, freq_bounds_ptr, dosage_bounds_ptr);
        if (reterr) {
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
  InitHistogramFromFileOrCommalist_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  InitHistogramFromFileOrCommalist_ret_MALFORMED_INPUT_2:
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  }
 InitHistogramFromFileOrCommalist_ret_1:
  CloseGzTokenStream(&gts);
  return reterr;
}

PglErr WriteAlleleFreqs(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uint64_t* founder_allele_dosages, const double* mach_r2_vals, const char* ref_binstr, const char* alt1_binstr, uint32_t variant_ct, uint32_t max_alt_allele_ct, uint32_t max_allele_slen, FreqRptFlags freq_rpt_flags, uint32_t max_thread_ct, uint32_t nonfounders, char* outname, char* outname_end) {
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
      const uintptr_t overflow_buf_size = kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + max_alt_allele_ct * (24 * k1LU) + 2 * max_allele_slen;
      const uint32_t output_zst = freq_rpt_flags & kfAlleleFreqZs;
      if (output_zst) {
        snprintf(&(outname_end[6 + counts]), kMaxOutfnameExtBlen - 7, ".zst");
      }
      reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
      if (reterr) {
        goto WriteAlleleFreqs_ret_1;
      }
      *cswritep++ = '#';
      const uint32_t chr_col = freq_rpt_flags & kfAlleleFreqColChrom;

      // includes trailing tab
      char* chr_buf;
      if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
        goto WriteAlleleFreqs_ret_NOMEM;
      }
      if (chr_col) {
        cswritep = strcpya(cswritep, "CHROM\t");
      }
      if (freq_rpt_flags & kfAlleleFreqColPos) {
        cswritep = strcpya(cswritep, "POS\t");
      } else {
        variant_bps = nullptr;
      }
      cswritep = strcpya(cswritep, "ID");
      const uint32_t ref_col = freq_rpt_flags & kfAlleleFreqColRef;
      if (ref_col) {
        cswritep = strcpya(cswritep, "\tREF");
      }
      const uint32_t alt1_col = freq_rpt_flags & kfAlleleFreqColAlt1;
      if (alt1_col) {
        cswritep = strcpya(cswritep, "\tALT1");
      }
      const uint32_t alt_col = freq_rpt_flags & kfAlleleFreqColAlt;
      if (alt_col) {
        cswritep = strcpya(cswritep, "\tALT");
      }
      const uint32_t reffreq_col = freq_rpt_flags & kfAlleleFreqColReffreq;
      if (reffreq_col) {
        cswritep = strcpya(cswritep, "\tREF_");
        if (counts) {
          cswritep = strcpya(cswritep, "CT");
        } else {
          cswritep = strcpya(cswritep, "FREQ");
        }
      }
      const uint32_t alt1freq_col = freq_rpt_flags & kfAlleleFreqColAlt1freq;
      if (alt1freq_col) {
        cswritep = strcpya(cswritep, "\tALT1_");
        if (counts) {
          cswritep = strcpya(cswritep, "CT");
        } else {
          cswritep = strcpya(cswritep, "FREQ");
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
          cswritep = strcpya(cswritep, "ALT_");
        }
        if (eq_num) {
          cswritep = strcpya(cswritep, "NUM_");
        }
        if (counts) {
          cswritep = memcpyl3a(cswritep, "CTS");
        } else {
          cswritep = strcpya(cswritep, "FREQS");
        }
      }
      const uint32_t mach_r2_col = freq_rpt_flags & kfAlleleFreqColMachR2;
      if (mach_r2_col) {
        cswritep = strcpya(cswritep, "\tMACH_R2");
      }
      const uint32_t nobs_col = freq_rpt_flags & kfAlleleFreqColNobs;
      if (nobs_col) {
        cswritep = strcpya(cswritep, "\tOBS_CT");
      }
      AppendBinaryEoln(&cswritep);

      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
      uint32_t variant_uidx = 0;
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t suppress_mach_r2 = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      uint32_t cur_allele_ct = 2;
      printf("--freq%s%s: 0%%", output_zst? " zs" : "", counts? " counts" : "");
      fflush(stdout);
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
        MovU32To1Bit(variant_include, &variant_uidx);
        if (variant_uidx >= chr_end) {
          do {
            ++chr_fo_idx;
            chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
          } while (variant_uidx >= chr_end);
          const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
          char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
          suppress_mach_r2 = (chr_idx == x_code) || (chr_idx == mt_code);
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
          for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
            if (Cswrite(&css, &cswritep)) {
              goto WriteAlleleFreqs_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
          }
          --cswritep;
        }
        const uint64_t* cur_allele_dosages = &(founder_allele_dosages[allele_idx_offset_base]);
        uint64_t tot_allele_dosage = cur_allele_dosages[0];
        for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
          tot_allele_dosage += cur_allele_dosages[allele_idx];
        }
        double tot_allele_dosage_recip = 0.0;
        if (!counts) {
          tot_allele_dosage_recip = 1.0 / u63tod(tot_allele_dosage);
        }
        if (reffreq_col) {
          *cswritep++ = '\t';
          const uint64_t ref_dosage = cur_allele_dosages[0];
          if (counts) {
            cswritep = dosagetoa(ref_dosage, cswritep);
          } else {
            cswritep = dtoa_g(u63tod(ref_dosage) * tot_allele_dosage_recip, cswritep);
          }
        }
        if (alt1freq_col) {
          *cswritep++ = '\t';
          const uint64_t alt1_dosage = cur_allele_dosages[1];
          if (counts) {
            cswritep = dosagetoa(alt1_dosage, cswritep);
          } else {
            cswritep = dtoa_g(u63tod(alt1_dosage) * tot_allele_dosage_recip, cswritep);
          }
        }
        if (freq_col) {
          *cswritep++ = '\t';
          for (uint32_t allele_idx = commalist_exclude_ref; allele_idx < cur_allele_ct; ++allele_idx) {
            const uint64_t cur_allele_dosage = cur_allele_dosages[allele_idx];
            if (counts) {
              cswritep = dosagetoa(cur_allele_dosage, cswritep);
            } else {
              cswritep = dtoa_g(u63tod(cur_allele_dosage) * tot_allele_dosage_recip, cswritep);
            }
            *cswritep++ = ',';
          }
          --cswritep;
        } else if (eq_col) {
          *cswritep++ = '\t';
          uint32_t at_least_one_entry = 0;
          for (uint32_t allele_idx = commalist_exclude_ref; allele_idx < cur_allele_ct; ++allele_idx) {
            const uint64_t cur_allele_dosage = cur_allele_dosages[allele_idx];
            if (eq_includez || cur_allele_dosage) {
              if (eq_num) {
                cswritep = u32toa(allele_idx, cswritep);
              } else {
                if (Cswrite(&css, &cswritep)) {
                  goto WriteAlleleFreqs_ret_WRITE_FAIL;
                }
                const char* cur_allele = cur_alleles[allele_idx];
                const char* cur_allele_end_or_eq = strchrnul(cur_allele, '=');
                if (*cur_allele_end_or_eq == '=') {
                  logerrputs("Error: --freq's 'eq', 'eqz', 'alteq', and 'alteqz' columns cannot be requested\nwhen an allele code contains a '='.\n");
                  goto WriteAlleleFreqs_ret_INCONSISTENT_INPUT;
                }
                cswritep = memcpya(cswritep, cur_allele, cur_allele_end_or_eq - cur_allele);
              }
              *cswritep++ = '=';
              if (counts) {
                cswritep = dosagetoa(cur_allele_dosage, cswritep);
              } else {
                cswritep = dtoa_g(u63tod(cur_allele_dosage) * tot_allele_dosage_recip, cswritep);
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
        if (mach_r2_col) {
          *cswritep++ = '\t';
          if (!suppress_mach_r2) {
            cswritep = dtoa_g(mach_r2_vals[variant_uidx], cswritep);
          } else {
            cswritep = strcpya(cswritep, "NA");
          }
        }
        if (nobs_col) {
          *cswritep++ = '\t';
          cswritep = dosagetoa(tot_allele_dosage, cswritep);
        }
        AppendBinaryEoln(&cswritep);
        if (Cswrite(&css, &cswritep)) {
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
      if (CswriteCloseNull(&css, cswritep)) {
        goto WriteAlleleFreqs_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      logprintfww("--freq%s%s: Allele %s (%s) written to %s .\n", output_zst? " zs" : "", counts? " counts" : "", counts? "counts" : "frequencies", nonfounders? "all samples" : "founders only", outname);
    }

    if (ref_binstr || alt1_binstr) {
      BigstackReset(bigstack_mark);
      double* ref_freq_bounds = nullptr;
      uint64_t* ref_dosage_bounds = nullptr;
      uint32_t* ref_histogram = nullptr;
      uint32_t ref_boundary_ct = 0;
      if (ref_binstr) {
        reterr = InitHistogramFromFileOrCommalist(ref_binstr, (freq_rpt_flags / kfAlleleFreqBinsRefFname) & 1, counts? nullptr : (&ref_freq_bounds), counts? (&ref_dosage_bounds) : nullptr, &ref_boundary_ct, &ref_histogram);
        if (reterr) {
          goto WriteAlleleFreqs_ret_1;
        }
      }
      double* alt1_freq_bounds = nullptr;
      uint64_t* alt1_dosage_bounds = nullptr;
      uint32_t* alt1_histogram = nullptr;
      uint32_t alt1_boundary_ct = 0;
      if (alt1_binstr) {
        reterr = InitHistogramFromFileOrCommalist(alt1_binstr, (freq_rpt_flags / kfAlleleFreqBinsAlt1Fname) & 1, counts? nullptr : (&alt1_freq_bounds), counts? (&alt1_dosage_bounds) : nullptr, &alt1_boundary_ct, &alt1_histogram);
        if (reterr) {
          goto WriteAlleleFreqs_ret_1;
        }
      }
      uint32_t variant_uidx = 0;
      if (!counts) {
        uint32_t cur_allele_ct = 2;
        for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
          MovU32To1Bit(variant_include, &variant_uidx);
          uintptr_t allele_idx_offset_base = variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[variant_uidx];
            cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
          }
          const uint64_t* cur_allele_dosages = &(founder_allele_dosages[allele_idx_offset_base]);
          const uint64_t ref_allele_dosage = cur_allele_dosages[0];
          const uint64_t alt1_allele_dosage = cur_allele_dosages[1];
          uint64_t tot_allele_dosage = ref_allele_dosage + alt1_allele_dosage;
          for (uint32_t allele_idx = 2; allele_idx < cur_allele_ct; ++allele_idx) {
            tot_allele_dosage += cur_allele_dosages[allele_idx];
          }
          const double tot_allele_dosage_recip = 1.0 / u63tod(tot_allele_dosage);
          if (ref_histogram) {
            ref_histogram[CountSortedSmallerD(ref_freq_bounds, ref_boundary_ct, ref_allele_dosage * tot_allele_dosage_recip)] += 1;
          }
          if (alt1_histogram) {
            alt1_histogram[CountSortedSmallerD(alt1_freq_bounds, alt1_boundary_ct, alt1_allele_dosage * tot_allele_dosage_recip)] += 1;
          }
        }
      } else {
        for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
          MovU32To1Bit(variant_include, &variant_uidx);
          uintptr_t allele_idx_offset_base = variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          }
          const uint64_t* cur_allele_dosages = &(founder_allele_dosages[allele_idx_offset_base]);
          if (ref_histogram) {
            ref_histogram[CountSortedSmallerU64(ref_dosage_bounds, ref_boundary_ct, cur_allele_dosages[0])] += 1;
          }
          if (alt1_histogram) {
            alt1_histogram[CountSortedSmallerU64(alt1_dosage_bounds, alt1_boundary_ct, cur_allele_dosages[1])] += 1;
          }
        }
      }
      for (uint32_t is_alt1 = 0; is_alt1 < 2; ++is_alt1) {
        const uint32_t* cur_histogram = is_alt1? alt1_histogram : ref_histogram;
        if (!cur_histogram) {
          continue;
        }
        char* outname_end2 = &(outname_end[6 + counts]);
        if (!is_alt1) {
          outname_end2 = strcpya(outname_end2, ".ref");
        } else {
          outname_end2 = strcpya(outname_end2, ".alt1");
        }
        snprintf(outname_end2, kMaxOutfnameExtBlen - 12, ".bins");
        if (fopen_checked(outname, FOPEN_WB, &outfile)) {
          goto WriteAlleleFreqs_ret_OPEN_FAIL;
        }
        char* textbuf = g_textbuf;
        char* textbuf_flush = &(textbuf[kMaxMediumLine]);
        char* write_iter = strcpya(textbuf, "#BIN_START\tOBS_CT" EOLN_STR);
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
            if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
              goto WriteAlleleFreqs_ret_WRITE_FAIL;
            }
          }
        } else {
          const uint64_t* cur_dosage_bounds = is_alt1? alt1_dosage_bounds : ref_dosage_bounds;
          for (uint32_t bin_idx = 0; bin_idx <= cur_boundary_ct; ++bin_idx) {
            if (!bin_idx) {
              *write_iter++ = '0';
            } else {
              write_iter = dosagetoa(1 + cur_dosage_bounds[bin_idx - 1], write_iter);
            }
            *write_iter++ = '\t';
            write_iter = u32toa(cur_histogram[bin_idx], write_iter);
            AppendBinaryEoln(&write_iter);
            if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
              goto WriteAlleleFreqs_ret_WRITE_FAIL;
            }
          }
        }
        if (fclose_flush_null(textbuf_flush, write_iter, &outfile)) {
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
    if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
      goto WriteGenoCounts_ret_NOMEM;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t raw_sample_ctv = BitCtToVecCt(raw_sample_ct);
    uintptr_t* sample_include_interleaved_vec = nullptr;
    uint32_t* sample_include_cumulative_popcounts = nullptr;
    uintptr_t* genovec = nullptr;
    uintptr_t* sex_male_interleaved_vec = nullptr;
    uint32_t* sex_male_cumulative_popcounts = nullptr;
    const uint32_t more_counts_needed = (simple_pgrp->fi.max_alt_allele_ct > 1) && (geno_counts_flags & (kfGenoCountsColRefalt1 | kfGenoCountsColRefalt | kfGenoCountsColHomalt1 | kfGenoCountsColAltxy | kfGenoCountsColXy | kfGenoCountsColHapalt1 | kfGenoCountsColHapalt | kfGenoCountsColHap | kfGenoCountsColNumeq));
    if (more_counts_needed) {
      if (bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &sample_include_interleaved_vec) ||
          bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
          bigstack_alloc_w(QuaterCtToWordCt(raw_sample_ct), &genovec) ||
          bigstack_alloc_w(raw_sample_ctv * kWordsPerVec, &sex_male_interleaved_vec) ||
          bigstack_alloc_u32(raw_sample_ctl, &sex_male_cumulative_popcounts)) {
        goto WriteGenoCounts_ret_NOMEM;
      }
      FillInterleavedMaskVec(sample_include, raw_sample_ctv, sample_include_interleaved_vec);
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      FillInterleavedMaskVec(sex_male, raw_sample_ctv, sex_male_interleaved_vec);
      FillCumulativePopcounts(sex_male, raw_sample_ctl, sex_male_cumulative_popcounts);
      // currently clear in front of every chromosome
      // PgrClearLdCache(simple_pgrp);
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + simple_pgrp->fi.max_alt_allele_ct * (24 * k1LU) + 2 * max_allele_slen;
    const uint32_t output_zst = geno_counts_flags & kfGenoCountsZs;
    OutnameZstSet(".gcount", output_zst, outname_end);
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
    if (reterr) {
      goto WriteGenoCounts_ret_1;
    }
    *cswritep++ = '#';
    const uint32_t chr_col = geno_counts_flags & kfGenoCountsColChrom;

    // includes trailing tab
    if (chr_col) {
      cswritep = strcpya(cswritep, "CHROM\t");
    }
    if (geno_counts_flags & kfGenoCountsColPos) {
      cswritep = strcpya(cswritep, "POS\t");
    } else {
      variant_bps = nullptr;
    }
    cswritep = strcpya(cswritep, "ID");
    const uint32_t ref_col = geno_counts_flags & kfGenoCountsColRef;
    if (ref_col) {
      cswritep = strcpya(cswritep, "\tREF");
    }
    const uint32_t alt1_col = geno_counts_flags & kfGenoCountsColAlt1;
    if (alt1_col) {
      cswritep = strcpya(cswritep, "\tALT1");
    }
    const uint32_t alt_col = geno_counts_flags & kfGenoCountsColAlt;
    if (alt_col) {
      cswritep = strcpya(cswritep, "\tALT");
    }
    const uint32_t homref_col = geno_counts_flags & kfGenoCountsColHomref;
    if (homref_col) {
      cswritep = strcpya(cswritep, "\tHOM_REF_CT");
    }
    const uint32_t refalt1_col = geno_counts_flags & kfGenoCountsColRefalt1;
    if (refalt1_col) {
      cswritep = strcpya(cswritep, "\tHET_REF_ALT1_CT");
    }
    const uint32_t refalt_col = geno_counts_flags & kfGenoCountsColRefalt;
    if (refalt_col) {
      cswritep = strcpya(cswritep, "\tHET_REF_ALT_CTS");
    }
    const uint32_t homalt1_col = geno_counts_flags & kfGenoCountsColHomalt1;
    if (homalt1_col) {
      cswritep = strcpya(cswritep, "\tHOM_ALT1_CT");
    }
    const uint32_t xy_col = geno_counts_flags & (kfGenoCountsColAltxy | kfGenoCountsColXy);
    const uint32_t xy_col_altonly = (geno_counts_flags / kfGenoCountsColAltxy) & 1;
    if (xy_col) {
      *cswritep++ = '\t';
      if (xy_col_altonly) {
        cswritep = strcpya(cswritep, "NONREF_");
      }
      cswritep = strcpya(cswritep, "DIPLOID_GENO_CTS");
    }
    const uint32_t hapref_col = geno_counts_flags & kfGenoCountsColHapref;
    if (hapref_col) {
      cswritep = strcpya(cswritep, "\tHAP_REF_CT");
    }
    const uint32_t hapalt1_col = geno_counts_flags & kfGenoCountsColHapalt1;
    if (hapalt1_col) {
      cswritep = strcpya(cswritep, "\tHAP_ALT1_CT");
    }
    const uint32_t hap_col = geno_counts_flags & (kfGenoCountsColHapalt | kfGenoCountsColHap);
    const uint32_t hap_col_altonly = (geno_counts_flags / kfGenoCountsColHapalt) & 1;
    if (hap_col) {
      if (hap_col_altonly) {
        cswritep = strcpya(cswritep, "\tHAP_ALT_CTS");
      } else {
        cswritep = strcpya(cswritep, "\tHAP_CTS");
      }
    }
    const uint32_t numeq_col = geno_counts_flags & kfGenoCountsColNumeq;
    if (numeq_col) {
      cswritep = strcpya(cswritep, "\tGENO_NUM_CTS");
    }
    const uint32_t missing_col = geno_counts_flags & kfGenoCountsColMissing;
    if (missing_col) {
      cswritep = strcpya(cswritep, "\tMISSING_CT");
    }
    const uint32_t nobs_col = geno_counts_flags & kfGenoCountsColNobs;
    if (nobs_col) {
      cswritep = strcpya(cswritep, "\tOBS_CT");
    }
    AppendBinaryEoln(&cswritep);

    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    const uintptr_t* cur_sample_include = nullptr;
    uint32_t is_autosomal_diploid = 0;
    uint32_t is_x = 0;
    uint32_t nobs_base = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t variant_uidx = 0;
    uint32_t homref_ct = 0;
    uint32_t het_ct = 0;
    uint32_t homalt1_ct = 0;
    uint32_t hapref_ct = 0;
    uint32_t hapalt1_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    printf("--geno-counts%s: 0%%", output_zst? " zs" : "");
    fflush(stdout);
    uint32_t cur_allele_ct = 2;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      MovU32To1Bit(variant_include, &variant_uidx);
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
        PgrClearLdCache(simple_pgrp);
        if (!is_autosomal_diploid) {
          if (chr_idx == y_code) {
            cur_sample_include = sex_male;
            nobs_base = male_ct;
          }
        }
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
        for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
          if (Cswrite(&css, &cswritep)) {
            goto WriteGenoCounts_ret_WRITE_FAIL;
          }
          cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
        }
        --cswritep;
      }
      uint32_t missing_ct;
      if ((cur_allele_ct == 2) || (!more_counts_needed)) {
        STD_ARRAY_KREF(uint32_t, 3) cur_raw_geno_cts = raw_geno_cts[variant_uidx];
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
            // chrY or other pure-haploid
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
            cswritep = strcpya(cswritep, "0/0=");
            cswritep = u32toa_x(homref_ct, ',', cswritep);
          }
          if (het_ct) {
            cswritep = strcpya(cswritep, "0/1=");
            cswritep = u32toa_x(het_ct, ',', cswritep);
          }
          if (homalt1_ct) {
            cswritep = strcpya(cswritep, "1/1=");
            cswritep = u32toa_x(homalt1_ct, ',', cswritep);
          }
          if (hapref_ct) {
            cswritep = strcpya(cswritep, "0=");
            cswritep = u32toa_x(hapref_ct, ',', cswritep);
          }
          if (hapalt1_ct) {
            cswritep = strcpya(cswritep, "1=");
            cswritep = u32toa_x(hapalt1_ct, ',', cswritep);
          }
          if (missing_ct != nobs_base) {
            --cswritep;
          } else {
            *cswritep++ = '.';
          }
        }
      } else {
        // todo; need PgrGetM()
        exit(63);
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
      if (Cswrite(&css, &cswritep)) {
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
    if (CswriteCloseNull(&css, cswritep)) {
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
      if (reterr) {
        goto WriteMissingnessReports_ret_1;
      }
      *cswritep++ = '#';
      const uint32_t scol_fid = FidColIsRequired(siip, missing_rpt_flags / kfMissingRptScolMaybefid);
      if (scol_fid) {
        cswritep = strcpya(cswritep, "FID\t");
      }
      cswritep = memcpyl3a(cswritep, "IID");
      const char* sample_ids = siip->sample_ids;
      const char* sids = siip->sids;
      const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
      const uintptr_t max_sid_blen = siip->max_sid_blen;
      const uint32_t scol_sid = SidColIsRequired(sids, missing_rpt_flags / kfMissingRptScolMaybesid);
      if (scol_sid) {
        cswritep = strcpya(cswritep, "\tSID");
      }
      const uint32_t scol_empty_pheno = (missing_rpt_flags & kfMissingRptScolMisspheno1) && (!pheno_ct);
      if (scol_empty_pheno) {
        cswritep = strcpya(cswritep, "\tMISS_PHENO1");
      }
      const uint32_t scol_phenos = (missing_rpt_flags & (kfMissingRptScolMisspheno1 | kfMissingRptScolMissphenos)) && pheno_ct;
      if (scol_phenos) {
        if (!(missing_rpt_flags & kfMissingRptScolMissphenos)) {
          pheno_ct = 1;
        }
        for (uintptr_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, &(pheno_names[pheno_idx * max_pheno_name_blen]));
          if (Cswrite(&css, &cswritep)) {
            goto WriteMissingnessReports_ret_WRITE_FAIL;
          }
        }
      }
      const uint32_t scol_nmiss_dosage = (missing_rpt_flags / kfMissingRptScolNmissDosage) & 1;
      if (scol_nmiss_dosage) {
        cswritep = strcpya(cswritep, "\tMISSING_DOSAGE_CT");
      }
      const uint32_t scol_nmiss = (missing_rpt_flags / kfMissingRptScolNmiss) & 1;
      if (scol_nmiss) {
        cswritep = strcpya(cswritep, "\tMISSING_CT");
      }
      const uint32_t scol_nmiss_hh = (missing_rpt_flags / kfMissingRptScolNmissHh) & 1;
      if (scol_nmiss_hh) {
        cswritep = strcpya(cswritep, "\tMISSING_AND_HETHAP_CT");
      }
      const uint32_t scol_hethap = (missing_rpt_flags / kfMissingRptScolHethap) & 1;
      if (scol_hethap) {
        cswritep = strcpya(cswritep, "\tHETHAP_CT");
      }
      const uint32_t scol_nobs = (missing_rpt_flags / kfMissingRptScolNobs) & 1;
      if (scol_nobs) {
        cswritep = strcpya(cswritep, "\tOBS_CT");
      }
      const uint32_t scol_fmiss_dosage = (missing_rpt_flags / kfMissingRptScolFmissDosage) & 1;
      if (scol_fmiss_dosage) {
        cswritep = strcpya(cswritep, "\tF_MISS_DOSAGE");
      }
      const uint32_t scol_fmiss = (missing_rpt_flags / kfMissingRptScolFmiss) & 1;
      if (scol_fmiss) {
        cswritep = strcpya(cswritep, "\tF_MISS");
      }
      const uint32_t scol_fmiss_hh = (missing_rpt_flags / kfMissingRptScolFmissHh) & 1;
      if (scol_fmiss_hh) {
        cswritep = strcpya(cswritep, "\tF_MISS_AND_HETHAP");
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
      uintptr_t sample_uidx = 0;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
        MovWTo1Bit(sample_include, &sample_uidx);
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
          for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
            *cswritep++ = '\t';
            // 'Y' - 'N' == 11
            *cswritep++ = 'Y' - 11 * IsSet(pheno_cols[pheno_idx].nonmiss, sample_uidx);
            if (Cswrite(&css, &cswritep)) {
              goto WriteMissingnessReports_ret_WRITE_FAIL;
            }
          }
        } else {
          if (scol_empty_pheno) {
            cswritep = strcpya(cswritep, "\tY");
          }
          if (Cswrite(&css, &cswritep)) {
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
      if (CswriteCloseNull(&css, cswritep)) {
        goto WriteMissingnessReports_ret_WRITE_FAIL;
      }
      BigstackReset(bigstack_mark);
      logprintfww("--missing: Sample missing data report written to %s .\n", outname);
    }
    if (!(missing_rpt_flags & kfMissingRptSampleOnly)) {
      const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
      char* chr_buf;  // includes trailing tab
      if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
        goto WriteMissingnessReports_ret_NOMEM;
      }
      const uintptr_t overflow_buf_size = kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + 2 * max_allele_slen;
      OutnameZstSet(".vmiss", output_zst, outname_end);
      reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, overflow_buf_size, &css, &cswritep);
      if (reterr) {
        goto WriteMissingnessReports_ret_1;
      }
      *cswritep++ = '#';
      const uint32_t chr_col = missing_rpt_flags & kfMissingRptVcolChrom;

      if (chr_col) {
        cswritep = strcpya(cswritep, "CHROM\t");
      }
      if (missing_rpt_flags & kfMissingRptVcolPos) {
        cswritep = strcpya(cswritep, "POS\t");
      } else {
        variant_bps = nullptr;
      }
      cswritep = strcpya(cswritep, "ID");
      const uint32_t ref_col = missing_rpt_flags & kfMissingRptVcolRef;
      if (ref_col) {
        cswritep = strcpya(cswritep, "\tREF");
      }
      const uint32_t alt1_col = missing_rpt_flags & kfMissingRptVcolAlt1;
      if (alt1_col) {
        cswritep = strcpya(cswritep, "\tALT1");
      }
      const uint32_t alt_col = missing_rpt_flags & kfMissingRptVcolAlt;
      if (alt_col) {
        cswritep = strcpya(cswritep, "\tALT");
      }
      const uint32_t nmiss_dosage_col = missing_rpt_flags & kfMissingRptVcolNmissDosage;
      if (nmiss_dosage_col) {
        cswritep = strcpya(cswritep, "\tMISSING_DOSAGE_CT");
      }
      const uint32_t nmiss_col = (missing_rpt_flags / kfMissingRptVcolNmiss) & 1;
      if (nmiss_col) {
        cswritep = strcpya(cswritep, "\tMISSING_CT");
      }
      const uint32_t nmiss_hh_col = (missing_rpt_flags / kfMissingRptVcolNmissHh) & 1;
      if (nmiss_hh_col) {
        cswritep = strcpya(cswritep, "\tMISSING_AND_HETHAP_CT");
      }
      const uint32_t hethap_col = (missing_rpt_flags / kfMissingRptVcolHethap) & 1;
      if (hethap_col) {
        cswritep = strcpya(cswritep, "\tHETHAP_CT");
      }
      const uint32_t nobs_col = (missing_rpt_flags / kfMissingRptVcolNobs) & 1;
      if (nobs_col) {
        cswritep = strcpya(cswritep, "\tOBS_CT");
      }
      const uint32_t fmiss_dosage_col = missing_rpt_flags & kfMissingRptVcolFmissDosage;
      if (fmiss_dosage_col) {
        cswritep = strcpya(cswritep, "\tF_MISS_DOSAGE");
      }
      const uint32_t fmiss_col = (missing_rpt_flags / kfMissingRptVcolFmiss) & 1;
      if (fmiss_col) {
        cswritep = strcpya(cswritep, "\tF_MISS");
      }
      const uint32_t fmiss_hh_col = (missing_rpt_flags / kfMissingRptVcolFmissHh) & 1;
      if (fmiss_hh_col) {
        cswritep = strcpya(cswritep, "\tF_MISS_AND_HETHAP");
      }
      const uint32_t fhethap_col = (missing_rpt_flags / kfMissingRptVcolFhethap) & 1;
      if (fhethap_col) {
        cswritep = strcpya(cswritep, "\tF_HETHAP");
      }
      AppendBinaryEoln(&cswritep);
      char nobs_str[16];
      nobs_str[0] = '\t';
      const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
      uint32_t nobs_slen = 0;
      uint32_t variant_uidx = 0;
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
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
        MovU32To1Bit(variant_include, &variant_uidx);
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
          for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
            if (Cswrite(&css, &cswritep)) {
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
        if (Cswrite(&css, &cswritep)) {
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
      if (CswriteCloseNull(&css, cswritep)) {
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

PglErr HardyMaj(__maybe_unused const uintptr_t* founder_info, __maybe_unused const uintptr_t* sex_nm, __maybe_unused const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, __maybe_unused uint32_t raw_sample_ct, uint32_t variant_ct, __maybe_unused PgenReader* simple_pgrp, STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_male_geno_cts), STD_ARRAY_PTR_DECL(uint32_t, 3, hwe_x_nosex_geno_cts)) {
  PglErr reterr = kPglRetSuccess;
  {
    uint32_t x_start = UINT32_MAX;
    uint32_t x_end = UINT32_MAX;
    if (hwe_x_male_geno_cts || hwe_x_nosex_geno_cts) {
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
      x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
      x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
    }
    uint32_t x_thresh = x_start;
    uint32_t is_x = 0;
    uint32_t cur_allele_ct = 2;
    uint32_t variant_uidx = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      MovU32To1Bit(variant_include, &variant_uidx);
      if (variant_uidx >= x_thresh) {
        is_x = (variant_uidx < x_end);
        if (!is_x) {
          x_thresh = UINT32_MAX;
        }
      }
      // could also iterate on nonzero maj_alleles[] bytes
      const uint32_t maj_idx = maj_alleles[variant_uidx];
      if (!maj_idx) {
        continue;
      }
      if (allele_idx_offsets) {
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
      }
      if (cur_allele_ct == 2) {
        STD_ARRAY_REF(uint32_t, 3) cur_geno_cts = hwe_geno_cts[variant_uidx];
        uint32_t tmp_ct = cur_geno_cts[0];
        cur_geno_cts[0] = cur_geno_cts[2];
        cur_geno_cts[2] = tmp_ct;
        if (is_x) {
          const uint32_t offset = variant_uidx - x_start;
          if (hwe_x_male_geno_cts) {
            STD_ARRAY_REF(uint32_t, 3) cur_x_male_geno_cts = hwe_x_male_geno_cts[offset];
            tmp_ct = cur_x_male_geno_cts[0];
            cur_x_male_geno_cts[0] = cur_x_male_geno_cts[2];
            cur_x_male_geno_cts[2] = tmp_ct;
          }
          if (hwe_x_nosex_geno_cts) {
            STD_ARRAY_REF(uint32_t, 3) cur_x_nosex_geno_cts = hwe_x_nosex_geno_cts[offset];
            tmp_ct = cur_x_nosex_geno_cts[0];
            cur_x_nosex_geno_cts[0] = cur_x_nosex_geno_cts[2];
            cur_x_nosex_geno_cts[2] = tmp_ct;
          }
        }
      } else {
        logerrputs("Error: --hardy/--hwe multiallelic support is under development.\n");
        exit(63);
      }
    }
  }
  while (0) {
  }
  return reterr;
}

// multithread globals
static const uintptr_t* g_variant_include = nullptr;

static const STD_ARRAY_PTR_DECL(uint32_t, 3, g_founder_raw_geno_cts) = nullptr;
static const STD_ARRAY_PTR_DECL(uint32_t, 3, g_founder_x_male_geno_cts) = nullptr;
static const STD_ARRAY_PTR_DECL(uint32_t, 3, g_founder_x_nosex_geno_cts) = nullptr;

static uint32_t* g_variant_uidx_starts = nullptr;
static double* g_hwe_x_pvals = nullptr;
static uint32_t g_x_start = 0;
static uint32_t g_hwe_x_ct = 0;
static uint32_t g_calc_thread_ct = 0;
static uint32_t g_hwe_midp = 0;

THREAD_FUNC_DECL ComputeHweXPvalsThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uintptr_t* variant_include = g_variant_include;
  const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts) = g_founder_raw_geno_cts;
  const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts) = g_founder_x_male_geno_cts;
  const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts) = g_founder_x_nosex_geno_cts;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t x_start = g_x_start;
  const uint32_t hwe_x_ct = g_hwe_x_ct;
  const uint32_t hwe_midp = g_hwe_midp;

  // this needs to be aligned with ComputeUidxStartPartition()
  const uint32_t variant_idx_end = (hwe_x_ct * (S_CAST(uint64_t, tidx) + 1)) / calc_thread_ct;
  uint32_t variant_idx = (hwe_x_ct * S_CAST(uint64_t, tidx)) / calc_thread_ct;

  double* hwe_x_pvals_iter = &(g_hwe_x_pvals[variant_idx]);
  uint32_t variant_uidx = g_variant_uidx_starts[tidx];
  uint32_t pct = 0;
  uint32_t next_print_variant_idx = variant_idx_end;
  if (!tidx) {
    next_print_variant_idx = variant_idx_end / 100;
  }
  uint32_t male_maj_ct = 0;
  uint32_t male_nonmaj_ct = 0;
  for (; variant_idx < variant_idx_end; ++variant_idx, ++variant_uidx) {
    MovU32To1Bit(variant_include, &variant_uidx);
    STD_ARRAY_KREF(uint32_t, 3) cur_raw_geno_cts = founder_raw_geno_cts[variant_uidx];
    uint32_t female_hommaj_ct = cur_raw_geno_cts[0];
    uint32_t female_het_maj_ct = cur_raw_geno_cts[1];
    uint32_t female_two_nonmaj_ct = cur_raw_geno_cts[2];
    if (founder_x_male_geno_cts) {
      STD_ARRAY_KREF(uint32_t, 3) cur_male_geno_cts = founder_x_male_geno_cts[variant_uidx - x_start];
      male_maj_ct = cur_male_geno_cts[0];
      female_hommaj_ct -= male_maj_ct;
      female_het_maj_ct -= cur_male_geno_cts[1];
      male_nonmaj_ct = cur_male_geno_cts[2];
      female_two_nonmaj_ct -= male_nonmaj_ct;
    }
    if (founder_x_nosex_geno_cts) {
      STD_ARRAY_KREF(uint32_t, 3) cur_nosex_geno_cts = founder_x_nosex_geno_cts[variant_uidx - x_start];
      female_hommaj_ct -= cur_nosex_geno_cts[0];
      female_het_maj_ct -= cur_nosex_geno_cts[1];
      female_two_nonmaj_ct -= cur_nosex_geno_cts[2];
    }
    *hwe_x_pvals_iter++ = HweXchrP(female_het_maj_ct, female_hommaj_ct, female_two_nonmaj_ct, male_maj_ct, male_nonmaj_ct, hwe_midp);
    if (variant_idx >= next_print_variant_idx) {
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
  THREAD_RETURN;
}

PglErr ComputeHweXPvals(const uintptr_t* variant_include, const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), uint32_t x_start, uint32_t hwe_x_ct, uint32_t hwe_midp, uint32_t calc_thread_ct, double** hwe_x_pvals_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    assert(hwe_x_ct);
    if (bigstack_alloc_d(hwe_x_ct, hwe_x_pvals_ptr)) {
      goto ComputeHweXPvals_ret_NOMEM;
    }
    bigstack_mark = g_bigstack_base;
    g_hwe_x_pvals = *hwe_x_pvals_ptr;

    if (calc_thread_ct > hwe_x_ct) {
      calc_thread_ct = hwe_x_ct;
    }
    pthread_t* threads;
    if (bigstack_alloc_thread(calc_thread_ct, &threads) ||
        bigstack_alloc_u32(calc_thread_ct, &g_variant_uidx_starts)) {
      goto ComputeHweXPvals_ret_NOMEM;
    }
    ComputeUidxStartPartition(variant_include, hwe_x_ct, calc_thread_ct, x_start, g_variant_uidx_starts);
    g_variant_include = variant_include;
    g_founder_raw_geno_cts = founder_raw_geno_cts;
    g_founder_x_male_geno_cts = founder_x_male_geno_cts;
    g_founder_x_nosex_geno_cts = founder_x_nosex_geno_cts;
    g_calc_thread_ct = calc_thread_ct;
    g_x_start = x_start;
    g_hwe_x_ct = hwe_x_ct;
    g_hwe_midp = hwe_midp;
    logprintf("Computing chrX Hardy-Weinberg %sp-values... ", hwe_midp? "mid" : "");
    fputs("0%", stdout);
    fflush(stdout);
    if (SpawnThreads(ComputeHweXPvalsThread, calc_thread_ct, threads)) {
      goto ComputeHweXPvals_ret_THREAD_CREATE_FAIL;
    }
    ComputeHweXPvalsThread(S_CAST(void*, 0));
    JoinThreads(calc_thread_ct, threads);
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
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr HardyReport(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), const double* hwe_x_pvals, uint32_t variant_ct, uint32_t hwe_x_ct, uint32_t max_allele_slen, double output_min_ln, HardyFlags hardy_flags, uint32_t max_thread_ct, uint32_t nonfounders, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    if (cip->haploid_mask[0] & 1) {
      logerrputs("Error: --hardy is pointless on an all-haploid genome.\n");
      goto HardyReport_ret_INCONSISTENT_INPUT;
    }
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    const uint32_t chr_code_end = cip->max_code + 1 + cip->name_ct;
    const uint32_t chr_code_endl = BitCtToWordCt(chr_code_end);
    const uintptr_t overflow_buf_size = RoundUpPow2(kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 512 + 2 * max_allele_slen, kCacheline);
    const uint32_t output_zst = hardy_flags & kfHardyZs;
    const double output_min_p = (output_min_ln < kLnDenormalMin)? 0 : exp(output_min_ln);
    uintptr_t overflow_buf_alloc = overflow_buf_size;
    if (output_zst) {
      overflow_buf_alloc += CstreamWkspaceReq(overflow_buf_size);
    }
    char* overflow_buf;
    uintptr_t* chr_skips;
    if (bigstack_alloc_c(overflow_buf_alloc, &overflow_buf) ||
        bigstack_alloc_w(chr_code_endl, &chr_skips)) {
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
    uint32_t chr_uidx = 0;
    for (uint32_t chr_skip_idx = 0; chr_skip_idx < chr_skip_ct; ++chr_skip_idx, ++chr_uidx) {
      MovU32To1Bit(chr_skips, &chr_uidx);
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
    const uint32_t chr_col = hardy_flags & kfHardyColChrom;
    const uint32_t ref_col = hardy_flags & kfHardyColRef;
    const uint32_t alt1_col = hardy_flags & kfHardyColAlt1;
    const uint32_t alt_col = hardy_flags & kfHardyColAlt;
    const uint32_t maj_col = hardy_flags & kfHardyColMaj;
    const uint32_t nonmaj_col = hardy_flags & kfHardyColNonmaj;
    const uint32_t gcounts = hardy_flags & (kfHardyColGcounts | kfHardyColGcount1col);
    const uint32_t gcount_1col = hardy_flags & kfHardyColGcount1col;
    const char gcount_delim = gcount_1col? ',' : '\t';
    const uint32_t hetfreq_cols = hardy_flags & kfHardyColHetfreq;
    const uint32_t p_col = hardy_flags & kfHardyColP;
    if (variant_ct) {
      OutnameZstSet(".hardy", output_zst, outname_end);
      reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
      if (reterr) {
        goto HardyReport_ret_1;
      }
      cswritep = overflow_buf;
      *cswritep++ = '#';

      // includes trailing tab
      char* chr_buf = nullptr;
      if (chr_col) {
        if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
          goto HardyReport_ret_NOMEM;
        }
        cswritep = strcpya(cswritep, "CHROM\t");
      }
      if (hardy_flags & kfHardyColPos) {
        cswritep = strcpya(cswritep, "POS\t");
      } else {
        variant_bps = nullptr;
      }
      cswritep = strcpya(cswritep, "ID");
      if (ref_col) {
        cswritep = strcpya(cswritep, "\tREF");
      }
      if (alt1_col) {
        cswritep = strcpya(cswritep, "\tALT1");
      }
      if (alt_col) {
        cswritep = strcpya(cswritep, "\tALT");
      }
      if (maj_col) {
        cswritep = strcpya(cswritep, "\tMAJ");
      }
      if (nonmaj_col) {
        cswritep = strcpya(cswritep, "\tNONMAJ");
      }
      if (gcounts) {
        if (gcount_1col) {
          cswritep = strcpya(cswritep, "\tGCOUNTS");
        } else {
          cswritep = strcpya(cswritep, "\tHOM_MAJ_CT\tHET_MAJ_CT\tTWO_NONMAJ_CT");
        }
      }
      if (hetfreq_cols) {
        cswritep = strcpya(cswritep, "\tO(HET_MAJ)\tE(HET_MAJ)");
      }
      if (p_col) {
        *cswritep++ = '\t';
        if (midp) {
          cswritep = strcpya(cswritep, "MIDP");
        } else {
          *cswritep++ = 'P';
        }
      }
      AppendBinaryEoln(&cswritep);
      uint32_t variant_uidx = 0;
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      uint32_t pct = 0;
      uint32_t next_print_variant_idx = variant_ct / 100;
      printf("--hardy%s%s: 0%%", output_zst? " zs" : "", midp? " midp" : "");
      fflush(stdout);
      uint32_t cur_allele_ct = 2;
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
        MovU32To1Bit(variant_include, &variant_uidx);
        if (chr_col) {
          if (variant_uidx >= chr_end) {
            uint32_t chr_idx;
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
              chr_idx = cip->chr_file_order[chr_fo_idx];
            } while ((variant_uidx >= chr_end) || IsSet(chr_skips, chr_idx));
            variant_uidx = AdvTo1Bit(variant_include, cip->chr_fo_vidx_start[chr_fo_idx]);
            char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
            *chr_name_end = '\t';
            chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
          }
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
          for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
            if (Cswrite(&css, &cswritep)) {
              goto HardyReport_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
          }
          --cswritep;
        }
        const uint32_t maj_idx = maj_alleles[variant_uidx];
        if (maj_col) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, cur_alleles[maj_idx]);
        }
        if (nonmaj_col) {
          *cswritep++ = '\t';
          for (uint32_t allele_idx = 0; allele_idx < cur_allele_ct; ++allele_idx) {
            if (allele_idx == maj_idx) {
              continue;
            }
            if (Cswrite(&css, &cswritep)) {
              goto HardyReport_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
          }
          --cswritep;
        }
        STD_ARRAY_KREF(uint32_t, 3) cur_geno_cts = founder_raw_geno_cts[variant_uidx];
        if (gcounts) {
          *cswritep++ = '\t';
          cswritep = u32toa_x(cur_geno_cts[0], gcount_delim, cswritep);
          cswritep = u32toa_x(cur_geno_cts[1], gcount_delim, cswritep);
          cswritep = u32toa(cur_geno_cts[2], cswritep);
        }
        if (hetfreq_cols) {
          *cswritep++ = '\t';
          const uint32_t tot_obs = cur_geno_cts[0] + cur_geno_cts[1] + cur_geno_cts[2];
          const double tot_obs_recip = 1.0 / u31tod(tot_obs);
          cswritep = dtoa_g(u31tod(cur_geno_cts[1]) * tot_obs_recip, cswritep);
          *cswritep++ = '\t';
          const double dbl_maj_freq = (cur_geno_cts[0] * 2 + cur_geno_cts[1]) * tot_obs_recip;
          const double expected_het_freq = dbl_maj_freq * (1.0 - dbl_maj_freq * 0.5);
          cswritep = dtoa_g(expected_het_freq, cswritep);
        }
        if (p_col) {
          // possible todo: multithread this
          *cswritep++ = '\t';
          const double hwe_p = HweP(cur_geno_cts[1], cur_geno_cts[0], cur_geno_cts[2], midp);
          cswritep = dtoa_g(MAXV(hwe_p, output_min_p), cswritep);
        }
        AppendBinaryEoln(&cswritep);
        if (Cswrite(&css, &cswritep)) {
          goto HardyReport_ret_WRITE_FAIL;
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
      if (CswriteCloseNull(&css, cswritep)) {
        goto HardyReport_ret_WRITE_FAIL;
      }
      putc_unlocked('\r', stdout);
      logprintfww("--hardy%s%s: Autosomal Hardy-Weinberg report (%s) written to %s .\n", output_zst? " zs" : "", midp? " midp" : "", nonfounders? "all samples" : "founders only", outname);
    }
    if (hwe_x_ct) {
      BigstackReset(chr_skips);
      OutnameZstSet(".hardy.x", output_zst, outname_end);
      reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
      if (reterr) {
        goto HardyReport_ret_1;
      }
      cswritep = overflow_buf;
      *cswritep++ = '#';

      // includes trailing tab
      char x_name_buf[8];
      uint32_t x_name_blen = 0;
      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      if (chr_col) {
        cswritep = strcpya(cswritep, "CHROM\t");
        char* write_iter = chrtoa(cip, x_code, x_name_buf);
        *write_iter++ = '\t';
        x_name_blen = write_iter - x_name_buf;
      }
      if (hardy_flags & kfHardyColPos) {
        cswritep = strcpya(cswritep, "POS\t");
      } else {
        variant_bps = nullptr;
      }
      cswritep = strcpya(cswritep, "ID");
      if (ref_col) {
        cswritep = strcpya(cswritep, "\tREF");
      }
      if (alt1_col) {
        cswritep = strcpya(cswritep, "\tALT1");
      }
      if (alt_col) {
        cswritep = strcpya(cswritep, "\tALT");
      }
      if (maj_col) {
        cswritep = strcpya(cswritep, "\tMAJ");
      }
      if (nonmaj_col) {
        cswritep = strcpya(cswritep, "\tNONMAJ");
      }
      if (gcounts) {
        if (gcount_1col) {
          cswritep = strcpya(cswritep, "\tGCOUNTS");
        } else {
          cswritep = strcpya(cswritep, "\tFEMALE_HOM_MAJ_CT\tFEMALE_HET_MAJ_CT\tFEMALE_TWO_NONMAJ_CT\tMALE_MAJ_CT\tMALE_NONMAJ_CT");
        }
      }
      if (hetfreq_cols) {
        cswritep = strcpya(cswritep, "\tO(FEMALE_HET_MAJ)\tE(FEMALE_HET_MAJ)");
      }
      const uint32_t sexaf_cols = hardy_flags & kfHardyColSexaf;
      if (sexaf_cols) {
        cswritep = strcpya(cswritep, "\tFEMALE_MAJ_FREQ\tMALE_MAJ_FREQ");
      }
      const uint32_t femalep_col = hardy_flags & kfHardyColFemalep;
      if (femalep_col) {
        cswritep = strcpya(cswritep, "\tFEMALE_ONLY_");
        if (midp) {
          cswritep = strcpya(cswritep, "MIDP");
        } else {
          *cswritep++ = 'P';
        }
      }
      if (p_col) {
        *cswritep++ = '\t';
        if (midp) {
          cswritep = strcpya(cswritep, "MIDP");
        } else {
          *cswritep++ = 'P';
        }
      }
      AppendBinaryEoln(&cswritep);
      fputs("--hardy: Writing chrX results...", stdout);
      fflush(stdout);
      const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
      const uint32_t x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
      uint32_t variant_uidx = x_start;
      uint32_t cur_allele_ct = 2;
      uint32_t male_maj_ct = 0;
      uint32_t male_nonmaj_ct = 0;
      for (uint32_t variant_idx = 0; variant_idx < hwe_x_ct; ++variant_idx, ++variant_uidx) {
        MovU32To1Bit(variant_include, &variant_uidx);
        cswritep = memcpya(cswritep, x_name_buf, x_name_blen);
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
          for (uint32_t allele_idx = 1; allele_idx < cur_allele_ct; ++allele_idx) {
            if (Cswrite(&css, &cswritep)) {
              goto HardyReport_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
          }
          --cswritep;
        }
        const uint32_t maj_idx = maj_alleles[variant_uidx];
        if (maj_col) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, cur_alleles[maj_idx]);
        }
        if (nonmaj_col) {
          *cswritep++ = '\t';
          for (uint32_t allele_idx = 0; allele_idx < cur_allele_ct; ++allele_idx) {
            if (allele_idx == maj_idx) {
              continue;
            }
            if (Cswrite(&css, &cswritep)) {
              goto HardyReport_ret_WRITE_FAIL;
            }
            cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
          }
          --cswritep;
        }
        STD_ARRAY_KREF(uint32_t, 3) cur_geno_cts = founder_raw_geno_cts[variant_uidx];
        uint32_t female_hommaj_ct = cur_geno_cts[0];
        uint32_t female_het_maj_ct = cur_geno_cts[1];
        uint32_t female_two_nonmaj_ct = cur_geno_cts[2];
        if (founder_x_male_geno_cts) {
          STD_ARRAY_KREF(uint32_t, 3) cur_male_geno_cts = founder_x_male_geno_cts[variant_uidx - x_start];
          male_maj_ct = cur_male_geno_cts[0];
          female_hommaj_ct -= male_maj_ct;
          female_het_maj_ct -= cur_male_geno_cts[1];
          male_nonmaj_ct = cur_male_geno_cts[2];
          female_two_nonmaj_ct -= male_nonmaj_ct;
        }
        if (founder_x_nosex_geno_cts) {
          STD_ARRAY_KREF(uint32_t, 3) cur_nosex_geno_cts = founder_x_nosex_geno_cts[variant_uidx - x_start];
          female_hommaj_ct -= cur_nosex_geno_cts[0];
          female_het_maj_ct -= cur_nosex_geno_cts[1];
          female_two_nonmaj_ct -= cur_nosex_geno_cts[2];
        }
        if (gcounts) {
          *cswritep++ = '\t';
          cswritep = u32toa_x(female_hommaj_ct, gcount_delim, cswritep);
          cswritep = u32toa_x(female_het_maj_ct, gcount_delim, cswritep);
          cswritep = u32toa_x(female_two_nonmaj_ct, gcount_delim, cswritep);
          cswritep = u32toa_x(male_maj_ct, gcount_delim, cswritep);
          cswritep = u32toa(male_nonmaj_ct, cswritep);
        }
        if (hetfreq_cols || sexaf_cols) {
          const uint32_t tot_female_obs = female_hommaj_ct + female_het_maj_ct + female_two_nonmaj_ct;
          const double tot_female_obs_recip = 1.0 / u31tod(tot_female_obs);
          const double dbl_ref_freq = (female_hommaj_ct * 2 + female_het_maj_ct) * tot_female_obs_recip;
          const double ref_freq = dbl_ref_freq * 0.5;
          if (hetfreq_cols) {
            *cswritep++ = '\t';
            cswritep = dtoa_g(u31tod(female_het_maj_ct) * tot_female_obs_recip, cswritep);
            *cswritep++ = '\t';
            const double expected_het_freq = dbl_ref_freq * (1.0 - ref_freq);
            cswritep = dtoa_g(expected_het_freq, cswritep);
          }
          if (sexaf_cols) {
            *cswritep++ = '\t';
            cswritep = dtoa_g(ref_freq, cswritep);
            *cswritep++ = '\t';
            const double male_ref_freq = u31tod(male_maj_ct) / u31tod(male_maj_ct + male_nonmaj_ct);
            cswritep = dtoa_g(male_ref_freq, cswritep);
          }
        }
        if (femalep_col) {
          *cswritep++ = '\t';
          const double female_hwe_p = HweP(female_het_maj_ct, female_hommaj_ct, female_two_nonmaj_ct, midp);
          cswritep = dtoa_g(MAXV(female_hwe_p, output_min_p), cswritep);
        }
        if (p_col) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(MAXV(hwe_x_pvals[variant_idx], output_min_p), cswritep);
        }
        AppendBinaryEoln(&cswritep);
        if (Cswrite(&css, &cswritep)) {
          goto HardyReport_ret_WRITE_FAIL;
        }
      }
      if (CswriteCloseNull(&css, cswritep)) {
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
  HardyReport_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 HardyReport_ret_1:
  CswriteCloseCond(&css, cswritep);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr WriteSnplist(const uintptr_t* variant_include, const char* const* variant_ids, uint32_t variant_ct, uint32_t output_zst, uint32_t max_thread_ct, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char* cswritep = nullptr;
  CompressStreamState css;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    OutnameZstSet(".snplist", output_zst, outname_end);
    reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, kCompressStreamBlock + kMaxIdSlen + 2, &css, &cswritep);
    if (reterr) {
      goto WriteSnplist_ret_1;
    }
    uint32_t variant_uidx = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      MovU32To1Bit(variant_include, &variant_uidx);
      cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
      AppendBinaryEoln(&cswritep);
      if (Cswrite(&css, &cswritep)) {
        goto WriteSnplist_ret_WRITE_FAIL;
      }
    }
    if (CswriteCloseNull(&css, cswritep)) {
      goto WriteSnplist_ret_WRITE_FAIL;
    }
    logprintfww("--write-snplist%s: Variant IDs written to %s .\n", output_zst? " zs" : "", outname);
  }
  while (0) {
  WriteSnplist_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
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
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
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
      write_iter = strcpya(write_iter, "FID\t");
    }
    write_iter = memcpyl3a(write_iter, "IID");
    if (write_sid) {
      write_iter = strcpya(write_iter, "\tSID");
    }
    if (write_parents) {
      write_iter = strcpya(write_iter, "\tPAT\tMAT");
    }
    if (write_sex) {
      write_iter = strcpya(write_iter, "\tSEX");
    }
    if (write_phenos || write_empty_pheno || write_sex) {
      // verify that no names are duplicated
      uint32_t* covar_name_htable;
      uint32_t covar_name_htable_size;
      if (HtableGoodSizeAlloc(covar_ct + write_sex, bigstack_left(), &covar_name_htable, &covar_name_htable_size)) {
        goto WriteCovar_ret_NOMEM;
      }
      // shouldn't be possible for this to fail
      PopulateStrboxHtable(covar_names, covar_ct, max_covar_name_blen, covar_name_htable_size, covar_name_htable);
      uint32_t max_xcovar_name_blen = max_covar_name_blen;
      if (write_sex) {
        // add "SEX"
        uint32_t hashval = Hashceil("SEX", 3, covar_name_htable_size);
        while (1) {
          const uint32_t cur_htable_entry = covar_name_htable[hashval];
          if (cur_htable_entry == UINT32_MAX) {
            covar_name_htable[hashval] = covar_ct;
            break;
          }
          if (strequal_k_unsafe(&(covar_names[cur_htable_entry * max_covar_name_blen]), "SEX")) {
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
        for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          const uint32_t cur_pheno_name_slen = strlen(pheno_name_iter);
          if (cur_pheno_name_slen < max_xcovar_name_blen) {
            const uint32_t cur_pheno_name_blen = cur_pheno_name_slen + 1;
            uint32_t hashval = Hashceil(pheno_name_iter, cur_pheno_name_slen, covar_name_htable_size);
            while (1) {
              uint32_t cur_htable_idval = covar_name_htable[hashval];
              if (cur_htable_idval >= covar_ct) {
                if (cur_htable_idval == UINT32_MAX) {
                  break;
                }
                if (strequal_k_unsafe(pheno_name_iter, "SEX")) {
                  logerrputs(write_sex? "Error: .cov file cannot have both a regular SEX column and a phenotype named\n'SEX'.  Exclude or rename one of these columns.\n" : "Error: .cov file cannot have a phenotype and a covariate with the same name.\n");
                  goto WriteCovar_ret_INCONSISTENT_INPUT;
                }
              } else {
                if (!memcmp(pheno_name_iter, &(covar_names[cur_htable_idval * max_covar_name_blen]), cur_pheno_name_blen)) {
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
          if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
            goto WriteCovar_ret_WRITE_FAIL;
          }
        }
      } else if (write_empty_pheno) {
        if (max_covar_name_blen > 6) {
          uint32_t hashval = Hashceil("PHENO1", 6, covar_name_htable_size);
          while (1) {
            uint32_t cur_htable_idval = covar_name_htable[hashval];
            if (cur_htable_idval >= covar_ct) {
              if (cur_htable_idval == UINT32_MAX) {
                break;
              }
            } else {
              if (strequal_k_unsafe(&(covar_names[cur_htable_idval * max_covar_name_blen]), "PHENO1")) {
                logerrputs("Error: .cov file cannot have a phenotype and a covariate with the same name.\n");
                goto WriteCovar_ret_INCONSISTENT_INPUT;
              }
            }
            if (++hashval == covar_name_htable_size) {
              hashval = 0;
            }
          }
        }
        write_iter = strcpya(write_iter, "\tPHENO1");
      }
    }
    for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx) {
      *write_iter++ = '\t';
      const char* cur_covar_name = &(covar_names[covar_idx * max_covar_name_blen]);
      const uint32_t cur_covar_name_slen = strlen(cur_covar_name);
      write_iter = memcpya(write_iter, cur_covar_name, cur_covar_name_slen);
      if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
        goto WriteCovar_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);

    uintptr_t sample_uidx = 0;
    uint32_t sample_uidx2 = 0;
    // not really necessary to make sample_uidx increment dependent on
    // new_sample_idx_to_old == nullptr
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      if (!new_sample_idx_to_old) {
        MovWTo1Bit(sample_include, &sample_uidx);
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
          write_iter = strcpya(write_iter, "NA");
        }
      }
      if (write_phenos) {
        for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
          *write_iter++ = '\t';
          write_iter = AppendPhenoStr(&(pheno_cols[pheno_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
          if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
            goto WriteCovar_ret_WRITE_FAIL;
          }
        }
      } else {
        if (write_empty_pheno) {
          *write_iter++ = '\t';
          write_iter = memcpya(write_iter, output_missing_pheno, omp_slen);
        }
        if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
          goto WriteCovar_ret_WRITE_FAIL;
        }
      }
      for (uint32_t covar_idx = 0; covar_idx < covar_ct; ++covar_idx) {
        *write_iter++ = '\t';
        write_iter = AppendPhenoStr(&(covar_cols[covar_idx]), output_missing_pheno, omp_slen, sample_uidx, write_iter);
        if (fwrite_ck(textbuf_flush, outfile, &write_iter)) {
          goto WriteCovar_ret_WRITE_FAIL;
        }
      }
      AppendBinaryEoln(&write_iter);
    }
    if (fclose_flush_null(textbuf_flush, write_iter, &outfile)) {
      goto WriteCovar_ret_WRITE_FAIL;
    }
    logprintfww("Covariates written to %s.\n", outname);
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

#ifdef __cplusplus
}  // namespace plink2
#endif
