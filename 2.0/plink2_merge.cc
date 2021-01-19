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
#include "plink2_data.h"
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
  pmerge_info_ptr->merge_parents_mode = kMergePhenoModeNmMatch;
  pmerge_info_ptr->merge_sex_mode = kMergePhenoModeNmMatch;
  pmerge_info_ptr->merge_pheno_mode = kMergePhenoModeNmMatch;
  pmerge_info_ptr->merge_xheader_mode = kMergeXheaderModeFirst;
  pmerge_info_ptr->merge_qual_mode = kMergeQualInfoModeNmFirst;
  pmerge_info_ptr->merge_filter_mode = kMergeFilterModeNmFirst;
  pmerge_info_ptr->merge_info_mode = kMergeQualInfoModeNmFirst;
  pmerge_info_ptr->merge_pheno_sort = kSortNone;
  pmerge_info_ptr->merge_info_sort = kSortNone;
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

typedef struct PmergeInputFilesetLlStruct {
  struct PmergeInputFilesetLlStruct* next;
  char* pgen_fname;
  char* pvar_fname;
  char* psam_fname;
  char* pgen_locked_fname;
  uint32_t sample_ct;
  uint32_t variant_ct;
  uint32_t max_line_blen;
  uint32_t max_single_pos_ct;
  uintptr_t max_single_pos_blen;
  uintptr_t grand_allele_ct;
  ChrIdx first_chr_idx;
  ChrIdx last_chr_idx;
  uint32_t first_pos;
  uint32_t last_pos;
} PmergeInputFilesetLl;

// Allocates at end of bigstack.
PmergeInputFilesetLl* AllocFilesetLlEntry(PmergeInputFilesetLl*** filesets_endpp) {
  const uintptr_t alloc_size = RoundUpPow2(sizeof(PmergeInputFilesetLl), kEndAllocAlign);
  if (unlikely(S_CAST(uintptr_t, g_bigstack_end - g_bigstack_base) < alloc_size)) {
    return nullptr;
  }
  g_bigstack_end -= alloc_size;
  PmergeInputFilesetLl* new_entry = R_CAST(PmergeInputFilesetLl*, g_bigstack_end);
  new_entry->next = nullptr;
  **filesets_endpp = new_entry;
  *filesets_endpp = &(new_entry->next);
  return new_entry;
}

PglErr LoadPmergeList(const char* list_fname, PmergeListMode mode, uint32_t main_fileset_present, PmergeInputFilesetLl*** filesets_endpp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    reterr = InitTextStream(list_fname, kTextStreamBlenFast, 1, &txs);
    if (unlikely(reterr)) {
      goto LoadPmergeList_ret_TSTREAM_FAIL;
    }
    const char* pgen_suffix;
    const char* pvar_suffix;
    const char* psam_suffix;
    switch (mode) {
    case kPmergeListModeBfile:
      pgen_suffix = ".bed";
      pvar_suffix = ".bim";
      psam_suffix = ".fam";
      break;
    case kPmergeListModeBpfile:
      pgen_suffix = ".pgen";
      pvar_suffix = ".bim";
      psam_suffix = ".fam";
      break;
    case kPmergeListModePfile:
      pgen_suffix = ".pgen";
      pvar_suffix = ".pvar";
      psam_suffix = ".psam";
      break;
    default:
      assert(mode == kPmergeListModePfileVzs);
      pgen_suffix = ".pgen";
      pvar_suffix = ".pvar.zst";
      psam_suffix = ".psam";
      break;
    }
    const uint32_t pgen_blen = strlen(pgen_suffix) + 1;
    const uint32_t pvar_blen = strlen(pvar_suffix) + 1;
    const uint32_t psam_blen = strlen(psam_suffix) + 1;
    const uint32_t max_single_token_slen = kPglFnamesize - 1 - MAXV(pgen_blen, pvar_blen);
    while (1) {
      const char* first_token_start = TextGet(&txs);
      if (!first_token_start) {
        break;
      }
      ++line_idx;
      PmergeInputFilesetLl* cur_entry = AllocFilesetLlEntry(filesets_endpp);
      if (unlikely(!cur_entry)) {
        goto LoadPmergeList_ret_NOMEM;
      }
      const char* first_token_end = CurTokenEnd(first_token_start);
      const uint32_t first_token_slen = first_token_end - first_token_start;
      const char* second_token_start = FirstNonTspace(first_token_end);
      if (IsEolnKns(*second_token_start)) {
        // Expand single token, using the --pmerge-list mode.
        if (unlikely(first_token_slen > max_single_token_slen)) {
          logerrprintf("Error: Entry on line %" PRIuPTR " of --pmerge-list file is too long.\n", line_idx);
          goto LoadPmergeList_ret_MALFORMED_INPUT;
        }
        char* fname_iter;
        if (unlikely(bigstack_end_alloc_c(first_token_slen * 3 + pgen_blen + pvar_blen + psam_blen, &fname_iter))) {
          goto LoadPmergeList_ret_NOMEM;
        }
        cur_entry->pgen_fname = fname_iter;
        fname_iter = memcpya(fname_iter, first_token_start, first_token_slen);
        fname_iter = memcpya(fname_iter, pgen_suffix, pgen_blen);
        cur_entry->pvar_fname = fname_iter;
        fname_iter = memcpya(fname_iter, first_token_start, first_token_slen);
        fname_iter = memcpya(fname_iter, pvar_suffix, pvar_blen);
        cur_entry->psam_fname = fname_iter;
        fname_iter = memcpya(fname_iter, first_token_start, first_token_slen);
        memcpy(fname_iter, psam_suffix, psam_blen);
        cur_entry->pgen_locked_fname = nullptr;
      } else {
        const char* second_token_end = CurTokenEnd(second_token_start);
        const char* third_token_start = FirstNonTspace(second_token_end);
        if (unlikely(IsEolnKns(*third_token_start))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --pmerge-list has exactly two tokens; 1 or 3 expected.\n", line_idx);
          goto LoadPmergeList_ret_MALFORMED_INPUT_WW;
        }
        const char* third_token_end = CurTokenEnd(third_token_start);
        const char* fourth_token_start = FirstNonTspace(third_token_end);
        if (unlikely(!IsEolnKns(*fourth_token_start))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --pmerge-list has more than 3 tokens.\n", line_idx);
          goto LoadPmergeList_ret_MALFORMED_INPUT_WW;
        }
        const uint32_t second_token_slen = second_token_end - second_token_start;
        const uint32_t third_token_slen = third_token_end - third_token_start;
        if (unlikely((first_token_slen >= kPglFnamesize) ||
                     (second_token_slen >= kPglFnamesize) ||
                     (third_token_slen >= kPglFnamesize))) {
          logerrprintf("Error: Filename on line %" PRIuPTR " of --pmerge-list file is too long.\n", line_idx);
          goto LoadPmergeList_ret_MALFORMED_INPUT;
        }
        char* fname_iter;
        if (unlikely(bigstack_end_alloc_c(first_token_slen + second_token_slen + third_token_slen + 3, &fname_iter))) {
          goto LoadPmergeList_ret_NOMEM;
        }
        cur_entry->pgen_fname = fname_iter;
        fname_iter = memcpyax(fname_iter, first_token_start, first_token_slen, '\0');
        cur_entry->pvar_fname = fname_iter;
        fname_iter = memcpyax(fname_iter, second_token_start, third_token_slen, '\0');
        cur_entry->psam_fname = fname_iter;
        memcpyx(fname_iter, third_token_start, third_token_slen, '\0');
        cur_entry->pgen_locked_fname = nullptr;
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto LoadPmergeList_ret_TSTREAM_FAIL;
    }
    if (main_fileset_present + line_idx < 2) {
      logerrputs("Error: --pmerge-list requires at least two filesets to be specified.\n");
      goto LoadPmergeList_ret_INCONSISTENT_INPUT;
    }
    logprintf("--pmerge-list: %" PRIuPTR " filesets specified%s.\n", main_fileset_present + line_idx, main_fileset_present? " (including main fileset)" : "");
  }
  while (0) {
  LoadPmergeList_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadPmergeList_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--pmerge-list file", &txs);
    break;
  LoadPmergeList_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  LoadPmergeList_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  LoadPmergeList_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

// Permanent allocations are at end of bigstack.
PglErr MergePsams(const PmergeInfo* pmip, const char* sample_sort_fname, MiscFlags misc_flags, SortMode sample_sort_mode, FamCol fam_cols, int32_t missing_pheno, uint32_t max_thread_ct, char* outname, char* outname_end, PmergeInputFilesetLl* filesets, SampleIdInfo* siip, uint32_t* sample_ctp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  const char* cur_fname = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    // First pass: determine sample IDs and phenotype set.
    // Intermission: sort sample and phenotype IDs.
    // Second pass: populate data structure, then call WritePsam().

    // Allocate a single linebuf that's used in all first-pass loads, so that
    // we can have the pheno strset grow from bigstack_end down in an
    // unfragmented manner.
    uintptr_t linebuf_capacity = MINV(kMaxLongLine, bigstack_left() / 4) + kDecompressChunkSize;
    char* linebuf;
    if (unlikely(bigstack_end_alloc_c(linebuf_capacity, &linebuf))) {
      goto MergePsams_ret_NOMEM;
    }
    // If --sample-inner-join is specified, we append raw null-terminated
    // sample-ID strings to the bottom of bigstack while scanning the first
    // file, construct the usual non-resizable string-hash-table data structure
    // when we're done.  Each ID then has a 0-based index, and we use
    // sample_include/cur_sample_include to track which sample IDs are present
    // in all files from that point on.
    // Otherwise, we add sample IDs to the resizable sample_id_strset data
    // structure and wait till we've processed all files before ordering the
    // elements.
    // Analogous thing happens for phenotype names if --pheno-inner-join is
    // specified.  Phenotype-name data structure grows down from the top of
    // bigstack.
    const PmergeFlags flags = pmip->flags;
    char* first_sample_ids_start = nullptr;
    const char** first_sample_ids = nullptr;
    uint32_t* first_sample_ids_htable = nullptr;
    uintptr_t* sample_include = nullptr;
    uintptr_t* cur_sample_include = nullptr;
    char** sample_id_strset = nullptr;
    uint32_t first_sample_id_ct = 0;
    uint32_t sample_id_table_size = 512;
    if (flags & kfPmergeSampleInnerJoin) {
      first_sample_ids_start = R_CAST(char*, g_bigstack_base);
      // - arena_bottom tracks the current append point
      // - first_sample_ids, first_sample_ids_htable, sample_include,
      //   cur_sample_include allocated after we're done scanning first file
      // - sample_id_table_size set to first_sample_ids_htable size at that
      //   point
    } else {
      if (unlikely(bigstack_calloc_cp(sample_id_table_size, &sample_id_strset))) {
        goto MergePsams_ret_NOMEM;
      }
    }
    const char** first_pheno_names = nullptr;
    uint32_t* first_pheno_names_htable = nullptr;
    uintptr_t* pheno_include = nullptr;
    uintptr_t* cur_pheno_include = nullptr;
    char** pheno_strset = nullptr;
    uint32_t first_pheno_ct = 0;
    uint32_t pheno_names_table_size = 512;
    if (!(flags & kfPmergePhenoInnerJoin)) {
      if (unlikely(bigstack_end_calloc_cp(pheno_names_table_size, &pheno_strset))) {
        goto MergePsams_ret_NOMEM;
      }
    }
    unsigned char* arena_bottom = g_bigstack_base;
    unsigned char* arena_top = g_bigstack_end;

    const uint32_t decompress_thread_ct = MAXV(max_thread_ct - 1, 1);
    // possible todo: track seen realpaths, skip duplicates
    PmergeInputFilesetLl* filesets_iter = filesets;
    char* idbuf = g_textbuf;
    uint32_t max_line_blen = 0;
    uint32_t max_sample_id_blen_m2 = 2;
    uintptr_t max_sid_blen = 0;
    uintptr_t max_paternal_id_blen = 2;
    uintptr_t max_maternal_id_blen = 2;
    uintptr_t max_pheno_name_blen = 0;
    uint32_t sample_ct = 0;
    uint32_t pheno_ct = 0;
    do {
      cur_fname = filesets_iter->psam_fname;
      reterr = TextStreamOpenEx(cur_fname, kMaxLongLine, linebuf_capacity, decompress_thread_ct, nullptr, linebuf, &txs);
      if (unlikely(reterr)) {
        goto MergePsams_ret_TSTREAM_FAIL;
      }
      // Worth optimizing this more than most text-reading loops, since we may
      // be processing a LOT of files.
      const char* line_start = TextLineEnd(&txs);
      for (line_idx = 1; ; ++line_idx) {
        if (unlikely(!TextGetUnsafe2K(&txs, &line_start))) {
          if (TextStreamErrcode2(&txs, &reterr)) {
            goto MergePsams_ret_TSTREAM_FAIL;
          }
          logerrprintfww("Error: No samples in %s.\n", cur_fname);
          goto MergePsams_ret_MALFORMED_INPUT;
        }
        if ((line_start[0] != '#') || tokequal_k(&(line_start[1]), "FID") || tokequal_k(&(line_start[1]), "IID")) {
          break;
        }
        const char* line_end = AdvPastDelim(line_start, '\n');
        const uint32_t line_blen = line_end - line_start;
        if (max_line_blen < line_blen) {
          max_line_blen = line_blen;
        }
        line_start = line_end;
      }
      uint32_t sid_present = 0;
      uint32_t postid_pat_col_idx = 0;
      uint32_t postid_mat_col_idx = 0;
      uint32_t fid_present;
      if (line_start[0] == '#') {
        const char* iid_end = &(line_start[4]);
        fid_present = (line_start[1] == 'F');
        if (fid_present) {
          const char* iid_start = FirstNonTspace(iid_end);
          iid_end = FirstPrechar(iid_start, 33);
          if (unlikely(!strequal_k(iid_start, "IID", iid_end - iid_start))) {
            goto MergePsams_ret_MALFORMED_INPUT;
          }
        }
        const char* token_start = FirstNonTspace(iid_end);
        if (tokequal_k(token_start, "SID")) {
          sid_present = 1;
          token_start = FirstNonTspace(&(token_start[3]));
        }
        const char* token_end;
        for (uint32_t postid_col_idx = 1; !IsEolnKns(*token_start); token_start = FirstNonTspace(token_end), ++postid_col_idx) {
          token_end = CurTokenEnd(token_start);
          const uint32_t token_slen = token_end - token_start;
          if (token_slen == 3) {
            if (unlikely(memequal_k(token_start, "FID", 3))) {
              if (fid_present) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate FID column in %s.\n", cur_fname);
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Improperly positioned FID column in %s (must be first).\n", cur_fname);
              }
              goto MergePsams_ret_MALFORMED_INPUT_WW;
            }
            if (unlikely(memequal_k(token_start, "IID", 3))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Duplicate IID column in %s.\n", cur_fname);
              goto MergePsams_ret_MALFORMED_INPUT_WW;
            }
            if (unlikely(memequal_k(token_start, "SID", 3))) {
              if (sid_present) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate SID column in %s.\n", cur_fname);
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Improperly positioned SID column in %s (must immediately follow IID).\n", cur_fname);
              }
              goto MergePsams_ret_MALFORMED_INPUT_WW;
            }
            if (memequal_k(token_start, "PAT", 3)) {
              if (unlikely(postid_pat_col_idx)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate PAT column in %s.\n", cur_fname);
                goto MergePsams_ret_MALFORMED_INPUT_WW;
              }
              postid_pat_col_idx = postid_col_idx;
              continue;
            }
            if (memequal_k(token_start, "MAT", 3)) {
              if (unlikely(postid_mat_col_idx)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate MAT column in %s.\n", cur_fname);
                goto MergePsams_ret_MALFORMED_INPUT_WW;
              }
              postid_mat_col_idx = postid_col_idx;
              continue;
            }
            if (memequal_k(token_start, "SEX", 3)) {
              continue;
            }
          } else if (token_slen > kMaxIdBlen) {
            logerrputs("Error: Phenotype names are limited to " MAX_ID_SLEN_STR " characters.\n");
            goto MergePsams_ret_MALFORMED_INPUT;
          }
          if (pheno_strset) {
            if (token_slen >= max_pheno_name_blen) {
              max_pheno_name_blen = token_slen + 1;
            }
            if (unlikely(StrsetAddEndResize(arena_bottom, token_start, token_slen, kMaxPhenoCt * 2, &pheno_strset, &pheno_names_table_size, &pheno_ct, &arena_top))) {
              if (pheno_ct == kMaxPhenoCt) {
                logerrputs("Error: " PROG_NAME_STR " does not support more than " MAX_PHENO_CT_STR " phenotypes.\n");
                goto MergePsams_ret_INCONSISTENT_INPUT;
              }
              goto MergePsams_ret_NOMEM;
            }
          } else {
            // max_pheno_name_blen calculation deferred since any name may not
            // appear in a later file
            if (!first_pheno_names_htable) {
              // StoreStringAtEnd(), without dst assignment
              if (unlikely(PtrWSubCk(arena_bottom, token_slen + 1, &arena_top))) {
                goto MergePsams_ret_NOMEM;
              }
              memcpyx(arena_top, token_start, token_slen, '\0');
              ++first_pheno_ct;
            } else if (pheno_ct) {
              const uint32_t pheno_idx = IdHtableFindNnt(token_start, first_pheno_names, first_pheno_names_htable, token_slen, pheno_names_table_size);
              if (pheno_idx != UINT32_MAX) {
                if (unlikely(IsSet(cur_pheno_include, pheno_idx))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Duplicate phenotype name '%s' in %s.\n", first_pheno_names[pheno_idx], cur_fname);
                  goto MergePsams_ret_MALFORMED_INPUT_WW;
                }
                SetBit(pheno_idx, cur_pheno_include);
              }
            }
          }
        }
        if ((!postid_pat_col_idx) != (!postid_mat_col_idx)) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s has a '%cAT' column without a '%cAT' column; either both or neither must be present.\n", cur_fname, postid_pat_col_idx? 'P' : 'M', postid_pat_col_idx? 'M' : 'P');
          goto MergePsams_ret_MALFORMED_INPUT_WW;
        }
        const char* line_end = AdvPastDelim(token_start, '\n');
        const uint32_t line_blen = line_end - line_start;
        if (max_line_blen < line_blen) {
          max_line_blen = line_blen;
        }
        line_start = line_end;
        ++line_idx;
        if (unlikely(!TextGetUnsafe2K(&txs, &line_start))) {
          if (TextStreamErrcode2(&txs, &reterr)) {
            goto MergePsams_ret_TSTREAM_FAIL;
          }
          logerrprintfww("Error: No samples in %s.\n", cur_fname);
          goto MergePsams_ret_MALFORMED_INPUT;
        }
      } else {
        // .fam
        fid_present = (fam_cols / kfFamCol1) & 1;
        if (fam_cols & kfFamCol34) {
          postid_pat_col_idx = 1;
          postid_mat_col_idx = 2;
        }
        if (fam_cols & kfFamCol6) {
          if (pheno_strset) {
            if (unlikely(StrsetAddEndResize(arena_bottom, "PHENO1", strlen("PHENO1"), kMaxPhenoCt * 2, &pheno_strset, &pheno_names_table_size, &pheno_ct, &arena_top))) {
              if (pheno_ct == kMaxPhenoCt) {
                logerrputs("Error: " PROG_NAME_STR " does not support more than " MAX_PHENO_CT_STR " phenotypes.\n");
                goto MergePsams_ret_INCONSISTENT_INPUT;
              }
              goto MergePsams_ret_NOMEM;
            }
            if (strlen("PHENO1") >= max_pheno_name_blen) {
              max_pheno_name_blen = strlen("PHENO1") + 1;
            }
          } else {
            if (!first_pheno_names_htable) {
              if (unlikely(PtrWSubCk(arena_bottom, strlen("PHENO1") + 1, &arena_top))) {
                goto MergePsams_ret_NOMEM;
              }
              strcpy_k(R_CAST(char*, arena_top), "PHENO1");
              first_pheno_ct = 1;
            } else if (pheno_ct) {
              const uint32_t pheno_idx = IdHtableFind("PHENO1", TO_CONSTCPCONSTP(first_pheno_names), first_pheno_names_htable, strlen("PHENO1"), pheno_names_table_size);
              if (pheno_idx != UINT32_MAX) {
                SetBit(pheno_idx, cur_pheno_include);
              }
            }
          }
        }
      }
      uint32_t first_parent_postid_col_idx = postid_pat_col_idx;
      uint32_t parent_col_skip = postid_mat_col_idx - postid_pat_col_idx;
      uintptr_t first_parent_max_blen = max_paternal_id_blen;
      uintptr_t second_parent_max_blen = max_maternal_id_blen;
      if (postid_mat_col_idx < postid_pat_col_idx) {
        first_parent_postid_col_idx = postid_mat_col_idx;
        parent_col_skip = -parent_col_skip;
        first_parent_max_blen = max_maternal_id_blen;
        second_parent_max_blen = max_paternal_id_blen;
      }
      const char* fid_start = &(g_one_char_strs[96]);
      uint32_t fid_slen = 1;
      const char* sid_start = &(g_one_char_strs[96]);
      uint32_t sid_slen = 1;
      while (1) {
        const char* token_start = line_start;
        if (fid_present) {
          fid_start = line_start;
          const char* fid_end = CurTokenEnd(fid_start);
          fid_slen = fid_end - fid_start;
          token_start = FirstNonTspace(fid_end);
          if (unlikely(IsEolnKns(*token_start))) {
            goto MergePsams_ret_MISSING_TOKENS;
          }
        }
        const char* iid_start = token_start;
        const char* iid_end = CurTokenEnd(iid_start);
        const char* token_end = iid_end;
        const uint32_t iid_slen = iid_end - iid_start;
        if (sid_present) {
          sid_start = FirstNonTspace(iid_end);
          token_end = CurTokenEnd(token_start);
          sid_slen = token_end - sid_start;
          if ((sid_slen > 1) || (sid_start[0] != '0')) {
            if (sid_slen >= max_sid_blen) {
              if (unlikely(sid_slen > kMaxIdBlen)) {
                logerrputs("Error: SIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
                goto MergePsams_ret_MALFORMED_INPUT;
              }
              max_sid_blen = sid_slen + 1;
            }
          }
        }
        if (fid_slen + iid_slen > max_sample_id_blen_m2) {
          max_sample_id_blen_m2 = fid_slen + iid_slen;
          if (unlikely(max_sample_id_blen_m2 > 2 * kMaxIdSlen)) {
            logerrputs("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
            goto MergePsams_ret_MALFORMED_INPUT;
          }
        }
        char* id_iter = memcpyax(idbuf, fid_start, fid_slen, '\t');
        id_iter = memcpyax(id_iter, iid_start, iid_slen, '\t');
        id_iter = memcpya(id_iter, sid_start, sid_slen);
        *id_iter = '\0';
        const uint32_t id_slen = id_iter - idbuf;
        if (sample_id_strset) {
          if (unlikely(StrsetAddResize(arena_top, idbuf, id_slen, 2U * 0x7ffffffe, sample_id_strset, &sample_id_table_size, &sample_ct, &arena_bottom))) {
            if (sample_ct == 0x7ffffffe) {
              logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
              goto MergePsams_ret_INCONSISTENT_INPUT;
            }
            goto MergePsams_ret_NOMEM;
          }
        } else {
          if (!first_sample_ids_htable) {
            if (unlikely(id_slen >= S_CAST(uintptr_t, arena_top - arena_bottom))) {
              goto MergePsams_ret_NOMEM;
            }
            arena_bottom = memcpyua(arena_bottom, idbuf, id_slen + 1);
            if (unlikely(first_sample_id_ct == 0x7ffffffe)) {
              logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
              goto MergePsams_ret_INCONSISTENT_INPUT;
            }
            ++first_sample_id_ct;
          } else {
            const uint32_t sample_idx = IdHtableFind(idbuf, first_sample_ids, first_sample_ids_htable, id_slen, sample_id_table_size);
            if (sample_idx != UINT32_MAX) {
              if (unlikely(IsSet(cur_sample_include, sample_idx))) {
                TabsToSpaces(idbuf);
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID \"%s\" in %s.\n", idbuf, cur_fname);
                goto MergePsams_ret_MALFORMED_INPUT_WW;
              }
              SetBit(sample_idx, cur_sample_include);
            }
          }
        }
        if (first_parent_postid_col_idx) {
          const char* first_parent_start = NextTokenMult(token_end, first_parent_postid_col_idx);
          if (unlikely(!first_parent_start)) {
            goto MergePsams_ret_MISSING_TOKENS;
          }
          const char* first_parent_end = CurTokenEnd(first_parent_start);
          const uintptr_t first_parent_slen = first_parent_end - first_parent_start;
          if (first_parent_slen >= first_parent_max_blen) {
            if (unlikely(first_parent_slen > kMaxIdSlen)) {
              logerrputs("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
              goto MergePsams_ret_MALFORMED_INPUT;
            }
            first_parent_max_blen = first_parent_slen + 1;
          }
          const char* second_parent_start = NextTokenMult(first_parent_end, parent_col_skip);
          if (unlikely(!second_parent_start)) {
            goto MergePsams_ret_MISSING_TOKENS;
          }
          token_end = CurTokenEnd(second_parent_start);
          const uintptr_t second_parent_slen = token_end - second_parent_start;
          if (second_parent_slen >= second_parent_max_blen) {
            if (unlikely(second_parent_slen > kMaxIdSlen)) {
              logerrputs("Error: FIDs and IIDs are limited to " MAX_ID_SLEN_STR " characters.\n");
              goto MergePsams_ret_MALFORMED_INPUT;
            }
            second_parent_max_blen = second_parent_slen + 1;
          }
        }
        const char* line_end = AdvPastDelim(token_end, '\n');
        const uint32_t line_blen = line_end - line_start;
        if (max_line_blen < line_blen) {
          max_line_blen = line_blen;
        }
        line_start = line_end;
        ++line_idx;
        if (!TextGetUnsafe2K(&txs, &line_start)) {
          break;
        }
        if (unlikely(line_start[0] == '#')) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx, cur_fname);
          goto MergePsams_ret_MALFORMED_INPUT_WW;
        }
      }
      if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
        goto MergePsams_ret_TSTREAM_FAIL;
      }
      if (unlikely(CleanupTextStream2(cur_fname, &txs, &reterr))) {
        goto MergePsams_ret_1;
      }
      if (postid_pat_col_idx) {
        if (postid_mat_col_idx > postid_pat_col_idx) {
          max_paternal_id_blen = first_parent_max_blen;
          max_maternal_id_blen = second_parent_max_blen;
        } else {
          max_paternal_id_blen = second_parent_max_blen;
          max_maternal_id_blen = first_parent_max_blen;
        }
      }
      if (first_sample_ids_start) {
        const uint32_t first_sample_ctl = BitCtToWordCt(first_sample_id_ct);
        if (filesets_iter == filesets) {
          ArenaBaseSet(arena_bottom, &arena_bottom);
          sample_id_table_size = GetHtableFastSize(first_sample_id_ct);
          if (unlikely((arena_bottom > arena_top) ||
                       arena_alloc_kcp(arena_top, first_sample_id_ct, &arena_bottom, &first_sample_ids) ||
                       arena_alloc_u32(arena_top, sample_id_table_size, &arena_bottom, &first_sample_ids_htable) ||
                       arena_alloc_w(arena_top, first_sample_ctl, &arena_bottom, &sample_include) ||
                       arena_alloc_w(arena_top, first_sample_ctl, &arena_bottom, &cur_sample_include))) {
            goto MergePsams_ret_NOMEM;
          }
          SetAllU32Arr(sample_id_table_size, first_sample_ids_htable);
          SetAllBits(first_sample_id_ct, sample_include);
          const char* first_sample_ids_iter = first_sample_ids_start;
          for (uint32_t sample_idx = 0; sample_idx != first_sample_id_ct; ++sample_idx) {
            first_sample_ids[sample_idx] = first_sample_ids_iter;
            const uint32_t slen = strlen(first_sample_ids_iter);
            if (unlikely(IdHtableAdd(first_sample_ids_iter, first_sample_ids, slen, sample_id_table_size, sample_idx, first_sample_ids_htable) != UINT32_MAX)) {
              char* mutable_id = K_CAST(char*, first_sample_ids_iter);
              TabsToSpaces(mutable_id);
              snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID \"%s\" in %s.\n", mutable_id, cur_fname);
              goto MergePsams_ret_MALFORMED_INPUT_WW;
            }
            first_sample_ids_iter = &(first_sample_ids_iter[slen + 1]);
          }
        } else {
          BitvecAnd(cur_sample_include, first_sample_ctl, sample_include);
          if (unlikely(AllWordsAreZero(sample_include, first_sample_ctl))) {
            snprintf(g_logbuf, kLogbufSize, "Error: No common samples in --pmerge%s --sample-inner-join job.\n", pmip->list_fname? "-list" : "");
            goto MergePsams_ret_MALFORMED_INPUT_2;
          }
        }
        ZeroWArr(first_sample_ctl, cur_sample_include);
      }
      if (!pheno_strset) {
        const uint32_t first_pheno_ctl = BitCtToWordCt(first_pheno_ct);
        if (filesets_iter == filesets) {
          // must handle pheno_ct == 0
          // (only really need to set first_pheno_names_htable to non-null in
          // that case, since that's how we currently recognize we aren't
          // processing the first .psam file)
          const char* pheno_names_iter = R_CAST(char*, arena_top);
          ArenaEndSet(arena_top, &arena_top);
          pheno_names_table_size = GetHtableFastSize(first_pheno_ct);
          if (unlikely(arena_end_alloc_kcp(arena_bottom, first_pheno_ct, &arena_top, &first_pheno_names) ||
                       arena_end_alloc_u32(arena_bottom, pheno_names_table_size, &arena_top, &first_pheno_names_htable) ||
                       arena_end_alloc_w(arena_bottom, first_pheno_ctl, &arena_top, &pheno_include) ||
                       arena_end_alloc_w(arena_bottom, first_pheno_ctl, &arena_top, &cur_pheno_include))) {
            goto MergePsams_ret_NOMEM;
          }
          SetAllU32Arr(pheno_names_table_size, first_pheno_names_htable);
          SetAllBits(first_pheno_ct, pheno_include);
          ZeroWArr(first_pheno_ctl, cur_pheno_include);
          pheno_ct = first_pheno_ct;
          for (uint32_t pheno_idx = 0; pheno_idx != first_pheno_ct; ++pheno_idx) {
            first_pheno_names[pheno_idx] = pheno_names_iter;
            const uint32_t slen = strlen(pheno_names_iter);
            if (unlikely(IdHtableAdd(pheno_names_iter, first_pheno_names, slen, pheno_names_table_size, pheno_idx, first_pheno_names_htable) != UINT32_MAX)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Duplicate phenotype name '%s' in %s.\n", pheno_names_iter, cur_fname);
              goto MergePsams_ret_MALFORMED_INPUT_WW;
            }
            pheno_names_iter = &(pheno_names_iter[slen + 1]);
          }
        } else if (pheno_ct) {
          BitvecAnd(cur_pheno_include, first_pheno_ctl, pheno_include);
          pheno_ct = PopcountWords(pheno_include, first_pheno_ctl);
          ZeroWArr(first_pheno_ctl, cur_pheno_include);
        }
      }
      filesets_iter = filesets_iter->next;
    } while (filesets_iter);
    // Intermission:
    // 1. Move final set of phenotype ID strings to bottom, compute
    //    max_pheno_name_blen if --pheno-inner-join.
    // 2. Compute max_sample_id_blen and max_sid_blen if --sample-inner-join.
    // 3. Free all end-of-bigstack allocations.
    // 4. Construct SampleIdInfo at end of bigstack, to be returned.
    // 5. Construct pheno_names, hash tables.
    char* pheno_names = nullptr;
    uint32_t* pheno_names_htable = nullptr;
    uint32_t* sample_id_htable;
    uintptr_t max_sample_id_blen;
    uint32_t sample_ctl;
    {
      char* pheno_names_tmp_start = R_CAST(char*, arena_bottom);
      char* pheno_names_tmp_end;
      if (pheno_strset) {
        const uintptr_t byte_ct = R_CAST(unsigned char*, pheno_strset) - arena_top;
        memmove(pheno_names_tmp_start, arena_top, byte_ct);
        pheno_names_tmp_end = &(pheno_names_tmp_start[byte_ct]);
      } else {
        max_pheno_name_blen = 0;
        char* pheno_names_write_iter = pheno_names_tmp_start;
        if (pheno_ct) {
          if (S_CAST(uintptr_t, first_pheno_names[0] - pheno_names_write_iter) < kMaxIdBlen) {
            // ensure memcpy is safe
            goto MergePsams_ret_NOMEM;
          }
          uintptr_t pheno_uidx_base = 0;
          uintptr_t cur_bits = pheno_include[0];
          for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
            const uint32_t pheno_uidx = BitIter1(pheno_include, &pheno_uidx_base, &cur_bits);
            const char* cur_pheno_name = first_pheno_names[pheno_uidx];
            const uint32_t blen = 1 + strlen(cur_pheno_name);
            if (blen > max_pheno_name_blen) {
              max_pheno_name_blen = blen;
            }
            pheno_names_write_iter = memcpya(pheno_names_write_iter, cur_pheno_name, blen);
          }
        }
        pheno_names_tmp_end = pheno_names_write_iter;
      }
      BigstackBaseSet(pheno_names_tmp_end);

      if (sample_id_strset) {
        max_sample_id_blen = max_sample_id_blen_m2 + 2;
      } else {
        const uint32_t first_sample_ctl = BitCtToWordCt(first_sample_id_ct);
        sample_ct = PopcountWords(sample_include, first_sample_ctl);
        max_sample_id_blen = 4;
        max_sid_blen = 0;
        uintptr_t sample_uidx_base = 0;
        uintptr_t cur_bits = sample_include[0];
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          const uint32_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
          const char* cur_sample_id = first_sample_ids[sample_uidx];
          const char* iid_end = AdvToDelim(AdvPastDelim(cur_sample_id, '\t'), '\t');
          const uint32_t sample_id_blen_m1 = iid_end - cur_sample_id;
          if (max_sample_id_blen <= sample_id_blen_m1) {
            max_sample_id_blen = sample_id_blen_m1 + 1;
          }
          const char* sid_start = &(iid_end[1]);
          const uint32_t sid_slen = strlen(sid_start);
          if ((sid_slen > 1) || (sid_start[0] != '0')) {
            if (max_sid_blen <= sid_slen) {
              max_sid_blen = sid_slen + 1;
            }
          }
        }
      }
      // defensive
      linebuf = nullptr;
      first_pheno_names = nullptr;
      first_pheno_names_htable = nullptr;
      pheno_include = nullptr;
      cur_pheno_include = nullptr;
      pheno_strset = nullptr;

      sample_ctl = BitCtToWordCt(sample_ct);
      siip->flags = kfSampleIdFidPresent | kfSampleIdParentsPresent;
      if (misc_flags & kfMiscStrictSid0) {
        // affects --indiv-sort
        siip->flags |= kfSampleIdStrictSid0;
      }
      siip->max_sample_id_blen = max_sample_id_blen;
      siip->max_sid_blen = max_sid_blen;
      siip->sids = nullptr;
      if (max_sid_blen) {
        if (unlikely(bigstack_end_alloc_c(max_sid_blen * sample_ct, &(siip->sids)))) {
          goto MergePsams_ret_NOMEM;
        }
      }
      // place this below optional siip->sids, so that
      // BigstackEndSet(siip->sample_ids) is always correct on success
      if (unlikely(bigstack_end_alloc_c(max_sample_id_blen * sample_ct, &(siip->sample_ids)))) {
        goto MergePsams_ret_NOMEM;
      }
      if (pheno_ct) {
        pheno_names_table_size = GetHtableFastSize(pheno_ct);
        if (unlikely(bigstack_end_alloc_c(max_pheno_name_blen * pheno_ct, &pheno_names) ||
                     bigstack_end_alloc_u32(pheno_names_table_size, &pheno_names_htable))) {
          goto MergePsams_ret_NOMEM;
        }
      }

      char* sample_ids_write_iter = siip->sample_ids;
      if (sample_id_strset) {
        const char* sample_ids_read_iter = R_CAST(char*, &(sample_id_strset[sample_id_table_size]));
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          const char* iid_end = AdvToDelim(AdvPastDelim(sample_ids_read_iter, '\t'), '\t');
          memcpyx(sample_ids_write_iter, sample_ids_read_iter, iid_end - sample_ids_read_iter, '\0');
          sample_ids_write_iter = &(sample_ids_write_iter[max_sample_id_blen]);
          const char* sid_end = strnul(iid_end);
          if (max_sid_blen) {
            const uint32_t sid_blen = sid_end - iid_end;
            const char* sid_start = &(iid_end[1]);
            memcpy(&(siip->sids[sample_idx * max_sid_blen]), sid_start, sid_blen);
          }
          sample_ids_read_iter = &(sid_end[1]);
        }
      } else {
        uintptr_t sample_uidx_base = 0;
        uintptr_t cur_bits = sample_include[0];
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          const uint32_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
          const char* cur_sample_id = first_sample_ids[sample_uidx];
          const char* iid_end = AdvToDelim(AdvPastDelim(cur_sample_id, '\t'), '\t');
          memcpyx(sample_ids_write_iter, cur_sample_id, iid_end - cur_sample_id, '\0');
          sample_ids_write_iter = &(sample_ids_write_iter[max_sample_id_blen]);
          if (max_sid_blen) {
            const char* sid_start = &(iid_end[1]);
            strcpy(&(siip->sids[sample_idx * max_sid_blen]), sid_start);
          }
        }
      }

      // defensive
      first_sample_ids_start = nullptr;
      first_sample_ids = nullptr;
      first_sample_ids_htable = nullptr;
      sample_include = nullptr;
      cur_sample_include = nullptr;
      sample_id_strset = nullptr;
      g_bigstack_base = bigstack_mark;

      if (pheno_ct) {
        const char* pheno_names_read_iter = pheno_names_tmp_start;
        // phenotypes are stored in reverse
        char* pheno_names_write_iter = &(pheno_names[pheno_ct * max_pheno_name_blen]);
        for (uint32_t pheno_idx = pheno_ct; pheno_idx; ) {
          --pheno_idx;
          pheno_names_write_iter -= max_pheno_name_blen;
          const uint32_t slen = strlen(pheno_names_read_iter);
          memcpy(pheno_names_write_iter, pheno_names_read_iter, slen + 1);
          pheno_names_read_iter = &(pheno_names_read_iter[slen + 1]);
        }
        const SortMode pheno_sort_mode = pmip->merge_pheno_sort;
        if (pheno_sort_mode == kSortAscii) {
          qsort(pheno_names, pheno_ct, max_pheno_name_blen, strcmp_overread_casted);
        } else if (pheno_sort_mode == kSortNatural) {
          qsort(pheno_names, pheno_ct, max_pheno_name_blen, strcmp_natural);
        }

        SetAllU32Arr(pheno_names_table_size, pheno_names_htable);
        pheno_names_read_iter = pheno_names;
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          const uint32_t slen = strlen(pheno_names_read_iter);
          HtableAddNondup(pheno_names_read_iter, slen, pheno_names_table_size, pheno_idx, pheno_names_htable);
          pheno_names_read_iter = &(pheno_names_read_iter[max_pheno_name_blen]);
        }
      }

      if (unlikely(bigstack_alloc_w(sample_ctl, &sample_include))) {
        goto MergePsams_ret_NOMEM;
      }
      SetAllBits(sample_ct, sample_include);
      char* sample_ids = siip->sample_ids;
      char* sids = siip->sids;
      if (sample_sort_mode != kSortNone) {
        if (sample_sort_mode == kSortFile) {
          // yes, this is a bit circuitous
          unsigned char* bigstack_end_mark2 = g_bigstack_end;
          uint32_t* new_sample_idx_to_old;
          char* sample_ids_tmp;
          if (unlikely(bigstack_end_alloc_u32(sample_ct, &new_sample_idx_to_old) ||
                       bigstack_end_alloc_c(sample_ct * max_sample_id_blen, &sample_ids_tmp))) {
            goto MergePsams_ret_NOMEM;
          }
          char* sids_tmp = nullptr;
          if (max_sid_blen) {
            if (unlikely(bigstack_end_alloc_c(sample_ct * max_sid_blen, &sids_tmp))) {
              goto MergePsams_ret_NOMEM;
            }
          }
          reterr = SampleSortFileMap(sample_include, siip, sample_sort_fname, sample_ct, sample_ct, &new_sample_idx_to_old);
          if (unlikely(reterr)) {
            goto MergePsams_ret_1;
          }
          char* new_sample_ids_iter = sample_ids_tmp;
          char* new_sids_iter = sids_tmp;
          for (uint32_t new_idx = 0; new_idx != sample_ct; ++new_idx) {
            const uint32_t old_idx = new_sample_idx_to_old[new_idx];
            strcpy(new_sample_ids_iter, &(sample_ids[old_idx * max_sample_id_blen]));
            new_sample_ids_iter = &(new_sample_ids_iter[max_sample_id_blen]);
            if (sids) {
              strcpy(new_sids_iter, &(sids[old_idx * max_sid_blen]));
              new_sids_iter = &(new_sids_iter[max_sid_blen]);
            }
          }
          memcpy(sample_ids, sample_ids_tmp, sample_ct * max_sample_id_blen);
          if (sids) {
            memcpy(sids, sids_tmp, sample_ct * max_sid_blen);
          }
          BigstackEndReset(bigstack_end_mark2);
        } else if (sample_sort_mode == kSortAscii) {
          qsort(siip->sample_ids, sample_ct, max_sample_id_blen, strcmp_overread_casted);
        } else {
          // natural-sort, even if sample_sort_mode == kSort0
          qsort(siip->sample_ids, sample_ct, max_sample_id_blen, strcmp_natural);
        }
      }

      sample_id_table_size = GetHtableFastSize(sample_ct);
      if (unlikely(bigstack_alloc_u32(sample_id_table_size, &sample_id_htable))) {
        goto MergePsams_ret_NOMEM;
      }
      SetAllU32Arr(sample_id_table_size, sample_id_htable);

      const char* sample_ids_iter = siip->sample_ids;
      const char* sids_iter = siip->sids;
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
        HtableAddNondup(cur_sample_id, slen, sample_id_table_size, sample_idx, sample_id_htable);
        sample_ids_iter = &(sample_ids_iter[max_sample_id_blen]);
      }
    }
    linebuf_capacity = MAXV(max_line_blen, kDecompressMinBlen) + kDecompressChunkSize;
    // max_{p,m}aternal_id_blen may be overestimated in --sample-inner-join
    // case, but that's harmless since we free this strbox soon enough
    ParentalIdInfo parental_id_info;
    uintptr_t* sex_nm;
    uintptr_t* sex_male;
    if (unlikely(bigstack_alloc_c(sample_ct * max_paternal_id_blen, &parental_id_info.paternal_ids) ||
                 bigstack_alloc_c(sample_ct * max_maternal_id_blen, &parental_id_info.maternal_ids) ||
                 bigstack_calloc_w(sample_ctl, &sex_nm) ||
                 bigstack_calloc_w(sample_ctl, &sex_male))) {
      goto MergePsams_ret_NOMEM;
    }
    {
      char* paternal_ids_iter = parental_id_info.paternal_ids;
      char* maternal_ids_iter = parental_id_info.maternal_ids;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        strcpy_k(paternal_ids_iter, "0");
        paternal_ids_iter = &(paternal_ids_iter[max_paternal_id_blen]);
        strcpy_k(maternal_ids_iter, "0");
        maternal_ids_iter = &(maternal_ids_iter[max_maternal_id_blen]);
      }
    }
    parental_id_info.max_paternal_id_blen = max_paternal_id_blen;
    parental_id_info.max_maternal_id_blen = max_maternal_id_blen;

    // unlike core pheno_cols, this just lives on bigstack
    PhenoCol* pheno_cols = nullptr;

    if (pheno_ct) {
      if (unlikely(BIGSTACK_ALLOC_X(PhenoCol, pheno_ct, &pheno_cols))) {
        goto MergePsams_ret_NOMEM;
      }
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        // we don't really distinguish between case/control and quantitative
        // phenotypes here (just use data.qt in both cases), and it should be
        // ok to overallocate a bit in the categorical-phenotype case
        if (unlikely(bigstack_calloc_w(sample_ctl, &(pheno_cols[pheno_idx].nonmiss)) ||
                     bigstack_alloc_d(sample_ct, &(pheno_cols[pheno_idx].data.qt)))) {
          goto MergePsams_ret_NOMEM;
        }
        // We use kPhenoDtypeCc instead of kPhenoTypeOther to indicate "type
        // not yet known", since that lets us leave it unchanged if all values
        // are missing.
        pheno_cols[pheno_idx].type_code = kPhenoDtypeCc;
        // defensive, nonnull_category_ct is ignored by WritePsam() except in
        // the 'SEX' phenotype case which can't happen here
        pheno_cols[pheno_idx].category_names = nullptr;
        pheno_cols[pheno_idx].nonnull_category_ct = 0;
      }
    }
    unsigned char* bigstack_mark2 = g_bigstack_base;
    const uint32_t kNonphenoColCt = 5;
    const uint32_t max_col_ct = pheno_ct + kNonphenoColCt;
    // note that ZeroWArr(pheno_ctl, pheno_include) is safe when pheno_ct == 0
    const uint32_t pheno_ctl = BitCtToWordCt(pheno_ct);
    const double missing_phenod = missing_pheno? S_CAST(double, missing_pheno) : HUGE_VAL;
    uint32_t* col_skips;
    uint32_t* col_types;
    const char** token_ptrs;
    uint32_t* token_slens;
    if (unlikely(bigstack_alloc_u32(max_col_ct, &col_skips) ||
                 bigstack_alloc_u32(max_col_ct, &col_types) ||
                 bigstack_alloc_kcp(max_col_ct, &token_ptrs) ||
                 bigstack_alloc_u32(max_col_ct, &token_slens) ||
                 bigstack_alloc_w(pheno_ctl, &pheno_include) ||
                 bigstack_alloc_c(linebuf_capacity, &linebuf))) {
      goto MergePsams_ret_NOMEM;
    }
    // sample-major, PAT before MAT
    const MergePhenoMode merge_parents_mode = pmip->merge_parents_mode;
    uintptr_t* parents_locked = nullptr;
    if (merge_parents_mode != kMergePhenoModeNmFirst) {
      if (unlikely(bigstack_calloc_w(BitCtToWordCt(2 * sample_ct), &parents_locked))) {
        goto MergePsams_ret_NOMEM;
      }
    }
    const MergePhenoMode merge_sex_mode = pmip->merge_sex_mode;
    uintptr_t* sex_locked = nullptr;
    if (merge_sex_mode != kMergePhenoModeNmFirst) {
      if (unlikely(bigstack_calloc_w(sample_ctl, &sex_locked))) {
        goto MergePsams_ret_NOMEM;
      }
    }
    // phenotype-major
    const MergePhenoMode merge_pheno_mode = pmip->merge_pheno_mode;
    uintptr_t* pheno_locked = nullptr;
    if (pheno_ct && (merge_pheno_mode != kMergePhenoModeNmFirst)) {
      // no need to check for overflow in 32-bit case since we would have
      // already run out of memory on data.qt allocations
      const uintptr_t bit_ct = S_CAST(uintptr_t, pheno_ct) * sample_ct;
      if (unlikely(bigstack_calloc_w(BitCtToWordCt(bit_ct), &pheno_locked))) {
        goto MergePsams_ret_NOMEM;
      }
    }
    // When --1 isn't specified, we need to be careful when interpreting
    // zeroes: if the phenotype winds up being quantitative, zero is a regular
    // value, but if it's binary, it's a missing value.
    // Except in the trivial "--merge-pheno-mode first" case, we treat it as a
    // missing value during the main pass, but set the relevant pheno_zero_seen
    // bit if that phenotype cell isn't locked yet.
    // We then correct the final phenotype values at the end, when we know
    // whether each column will appear to be binary or quantitative in the
    // merged .psam.
    // Note that it is possible for a quantitative column to degrade to a
    // binary column due to e.g. conflicts knocking out all values not in {-9,
    // 0, 1, 2}.  There isn't much we can do about that here, though.
    const uint32_t affection_01 = (misc_flags / kfMiscAffection01) & 1;
    uintptr_t* pheno_zero_seen = nullptr;
    if (pheno_ct && (!affection_01) && (merge_pheno_mode != kMergePhenoModeFirst)) {
      const uintptr_t bit_ct = S_CAST(uintptr_t, pheno_ct) * sample_ct;
      if (unlikely(bigstack_calloc_w(BitCtToWordCt(bit_ct), &pheno_zero_seen))) {
        goto MergePsams_ret_NOMEM;
      }
    }

    if (pheno_ct) {
      logprintf("--pmerge%s: %u sample%s and %u phenotype%s present.\n", pmip->list_fname? "-list" : "", sample_ct, (sample_ct == 1)? "" : "s", pheno_ct, (pheno_ct == 1)? "" : "s");
    } else {
      logprintf(g_logbuf, kLogbufSize, "--pmerge%s: %u sample%s present.\n", pmip->list_fname? "-list" : "", sample_ct, (sample_ct == 1)? "" : "s");
    }
    //  [category_names] [category_names_htable]            [strings]
    // ^                                        ^          ^         ^
    // |                                        |          |         |
    // g_bigstack_base                    arena_bottom arena_top g_bigstack_end
    //
    // category-name strings are allocated at the end of bigstack, and need to
    // remain allocated until WritePsam() returns.
    // category_names_htable is repositioned, resized, and rebuilt whenever
    // category_names would otherwise overflow.
    // No other allocations allowed past this point, until WritePsam() call.
    const char** category_names = nullptr;
    uint32_t* category_names_htable = nullptr;
    uint32_t category_names_ct = 0;
    // category_names_capacity == category_names_htable_size / 2
    uint32_t category_names_htable_size = 0;
    arena_top = g_bigstack_end;
    if (pheno_ct) {
      category_names_htable_size = 512;
      if (unlikely(bigstack_alloc_kcp(category_names_htable_size / 2, &category_names) ||
                   bigstack_alloc_u32(category_names_htable_size, &category_names_htable))) {
        goto MergePsams_ret_NOMEM;
      }
      SetAllU32Arr(category_names_htable_size, category_names_htable);
      // copy into arena so overread is safe
      const uint32_t missing_catname_slen = strlen(g_missing_catname);
      if (unlikely(StoreStringAtEndK(g_bigstack_base, g_missing_catname, missing_catname_slen, &arena_top, &(category_names[0])))) {
        goto MergePsams_ret_NOMEM;
      }
      const uint32_t hashval = Hashceil(g_missing_catname, missing_catname_slen, category_names_htable_size);
      category_names_htable[0] = hashval;
      category_names_ct = 1;
    }
    arena_bottom = g_bigstack_base;
    filesets_iter = filesets;
    do {
      cur_fname = filesets_iter->psam_fname;
      reterr = TextStreamOpenEx(cur_fname, kMaxLongLine, linebuf_capacity, decompress_thread_ct, nullptr, linebuf, &txs);
      if (unlikely(reterr)) {
        goto MergePsams_ret_TSTREAM_REWIND_FAIL;
      }
      const char* line_start = TextLineEnd(&txs);
      for (line_idx = 1; ; ++line_idx) {
        if (unlikely(!TextGetUnsafe2K(&txs, &line_start))) {
          reterr = TextStreamRawErrcode(&txs);
          goto MergePsams_ret_TSTREAM_REWIND_FAIL;
        }
        if ((line_start[0] != '#') || tokequal_k(&(line_start[1]), "FID") || tokequal_k(&(line_start[1]), "IID")) {
          break;
        }
        line_start = AdvPastDelim(line_start, '\n');
      }
      ZeroWArr(pheno_ctl, pheno_include);
      uint32_t sid_present = 0;
      uint32_t parents_present = 0;
      uint32_t sex_present = 0;
      uint32_t relevant_postfid_col_ct = 1;
      uint32_t fid_present;
      col_types[0] = 0;
      if (line_start[0] == '#') {
        fid_present = (line_start[1] == 'F');
        if (!fid_present) {
          col_skips[0] = 0;
        }
        const char* token_end = &(line_start[4]);
        for (uint32_t col_idx = 1; ; ++col_idx) {
          const char* token_start = FirstNonTspace(token_end);
          if (IsEolnKns(*token_start)) {
            break;
          }
          token_end = CurTokenEnd(token_start);
          const uint32_t token_slen = token_end - token_start;
          uint32_t cur_col_type = UINT32_MAX;
          if (token_slen == 3) {
            if (memequal_k(token_start, "IID", 3)) {
              cur_col_type = 0;
            } else if (memequal_k(token_start, "SID", 3)) {
              sid_present = 1;
              cur_col_type = 1;
            } else if (memequal_k(token_start, "PAT", 3)) {
              cur_col_type = 2;
              parents_present = 1;
            } else if (memequal_k(token_start, "MAT", 3)) {
              cur_col_type = 3;
            } else if (memequal_k(token_start, "SEX", 3)) {
              cur_col_type = 4;
              sex_present = 1;
            }
          }
          if (cur_col_type == UINT32_MAX) {
            if (token_slen < max_pheno_name_blen) {
              cur_col_type = StrboxHtableFindNnt(token_start, pheno_names, pheno_names_htable, max_pheno_name_blen, token_slen, pheno_names_table_size);
            }
            if (cur_col_type == UINT32_MAX) {
              continue;
            }
            SetBit(cur_col_type, pheno_include);
            cur_col_type += kNonphenoColCt;
          }
          col_skips[relevant_postfid_col_ct] = col_idx;
          col_types[relevant_postfid_col_ct++] = cur_col_type;
        }
        for (uint32_t rpf_col_idx = relevant_postfid_col_ct - 1; rpf_col_idx; --rpf_col_idx) {
          col_skips[rpf_col_idx] -= col_skips[rpf_col_idx - 1];
        }
        line_start = AdvPastDelim(token_end, '\n');
        ++line_idx;
      } else {
        // .fam
        fid_present = (fam_cols / kfFamCol1) & 1;
        col_skips[0] = fid_present;
        if (fam_cols & kfFamCol34) {
          col_skips[relevant_postfid_col_ct] = 1;
          col_types[relevant_postfid_col_ct++] = 2;
          col_skips[relevant_postfid_col_ct] = 1;
          col_types[relevant_postfid_col_ct++] = 3;
          parents_present = 1;
        }
        if (fam_cols & kfFamCol5) {
          col_skips[relevant_postfid_col_ct] = 1;
          col_types[relevant_postfid_col_ct++] = 4;
          sex_present = 1;
        }
        if (fam_cols & kfFamCol6) {
          col_skips[relevant_postfid_col_ct] = 1;
          const uint32_t pheno_idx = StrboxHtableFind("PHENO1", pheno_names, pheno_names_htable, max_pheno_name_blen, strlen("PHENO1"), pheno_names_table_size);
          col_types[relevant_postfid_col_ct++] = kNonphenoColCt + pheno_idx;
          SetBit(pheno_idx, pheno_include);
        }
      }

      const uint32_t cur_pheno_ct = PopcountWords(pheno_include, pheno_ctl);
      const char* sample_ids = siip->sample_ids;
      const char* sids = siip->sids;
      uint32_t cur_sample_ct = 0;
      uint32_t sid_slen = 1;
      const char* line_iter;
      for (; TextGetUnsafe2K(&txs, &line_start); line_start = AdvPastDelim(line_iter, '\n'), ++line_idx) {
        line_iter = TokenLexK0(line_start, col_types, col_skips, relevant_postfid_col_ct, token_ptrs, token_slens);
        if (unlikely(!line_iter)) {
          goto MergePsams_ret_MISSING_TOKENS;
        }
        char* id_iter = idbuf;
        if (fid_present) {
          const char* fid_end = CurTokenEnd(line_start);
          id_iter = memcpya(id_iter, line_start, fid_end - line_start);
        } else {
          *id_iter++ = '0';
        }
        *id_iter++ = '\t';
        id_iter = memcpya(id_iter, token_ptrs[0], token_slens[0]);
        uint32_t sample_idx;
        if (!sids) {
          *id_iter = '\0';
          const uint32_t id_slen = id_iter - idbuf;
          sample_idx = StrboxHtableFind(idbuf, sample_ids, sample_id_htable, max_sample_id_blen, id_slen, sample_id_table_size);
        } else {
          *id_iter++ = '\t';
          char* sid_start = id_iter;
          if (sid_present) {
            sid_slen = token_slens[1];
            id_iter = memcpya(id_iter, token_ptrs[1], sid_slen);
          } else {
            *id_iter++ = '0';
          }
          *id_iter = '\0';
          // probably want to move this to a plink2_common function
          uint32_t hashval = Hashceil(idbuf, id_iter - idbuf, sample_id_table_size);
          const uint32_t fid_iid_blen = sid_start - idbuf;
          sid_start[-1] = '\0';
          while (1) {
            sample_idx = sample_id_htable[hashval];
            if (sample_idx == UINT32_MAX) {
              break;
            }
            if (memequal(idbuf, &(sample_ids[sample_idx * max_sample_id_blen]), fid_iid_blen) && memequal(sid_start, &(sids[sample_idx * max_sid_blen]), sid_slen + 1)) {
              break;
            }
            if (++hashval == sample_id_table_size) {
              hashval = 0;
            }
          }
        }
        if (sample_idx == UINT32_MAX) {
          continue;
        }
        ++cur_sample_ct;
        if (parents_present) {
          char* dad_id_dst = &(parental_id_info.paternal_ids[sample_idx * max_paternal_id_blen]);
          char* mom_id_dst = &(parental_id_info.maternal_ids[sample_idx * max_maternal_id_blen]);
          if (parents_locked) {
            if (!IsSet(parents_locked, sample_idx * 2)) {
              const char* dad_id = token_ptrs[2];
              const uint32_t dad_slen = token_slens[2];
              if (merge_parents_mode == kMergePhenoModeNmMatch) {
                if (((dad_slen != 1) || (dad_id[0] != '0')) && (!strequal_unsafe(dad_id_dst, dad_id, dad_slen))) {
                  if (memequal_k(dad_id_dst, "0", 2)) {
                    memcpyx(dad_id_dst, dad_id, dad_slen, '\0');
                  } else {
                    strcpy_k(dad_id_dst, "0");
                    SetBit(sample_idx * 2, parents_locked);
                  }
                }
              } else {
                memcpyx(dad_id_dst, dad_id, dad_slen, '\0');
                SetBit(sample_idx * 2, parents_locked);
              }
            }
            if (!IsSet(parents_locked, sample_idx * 2 + 1)) {
              const char* mom_id = token_ptrs[3];
              const uint32_t mom_slen = token_slens[3];
              if (merge_parents_mode == kMergePhenoModeNmMatch) {
                if (((mom_slen != 1) || (mom_id[0] != '0')) && (!strequal_unsafe(mom_id_dst, mom_id, mom_slen))) {
                  if (memequal_k(mom_id_dst, "0", 2)) {
                    memcpyx(mom_id_dst, mom_id, mom_slen, '\0');
                  } else {
                    strcpy_k(mom_id_dst, "0");
                    SetBit(sample_idx * 2, parents_locked);
                  }
                }
              } else {
                memcpyx(mom_id_dst, mom_id, mom_slen, '\0');
                SetBit(sample_idx * 2 + 1, parents_locked);
              }
            }
          } else {
            // nm-first, skip if not missing
            if (memequal_k(dad_id_dst, "0", 2)) {
              memcpyx(dad_id_dst, token_ptrs[2], token_slens[2], '\0');
            }
            if (memequal_k(mom_id_dst, "0", 2)) {
              memcpyx(mom_id_dst, token_ptrs[3], token_slens[3], '\0');
            }
          }
        }
        if (sex_present) {
          if (!sex_locked) {
            if (!IsSet(sex_nm, sample_idx)) {
              // nm-first and previously missing
              const uint32_t cur_sex_code = (token_slens[4] == 1)? CharToSex(token_ptrs[4][0]) : 0;
              if (cur_sex_code) {
                SetBit(sample_idx, sex_nm);
                AssignBit(sample_idx, cur_sex_code & 1, sex_male);
              }
            }
          } else {
            if (!IsSet(sex_locked, sample_idx)) {
              const uint32_t cur_sex_code = (token_slens[4] == 1)? CharToSex(token_ptrs[4][0]) : 0;
              if (merge_sex_mode == kMergePhenoModeNmMatch) {
                if (cur_sex_code) {
                  if (!IsSet(sex_nm, sample_idx)) {
                    SetBit(sample_idx, sex_nm);
                    AssignBit(sample_idx, cur_sex_code & 1, sex_male);
                  } else {
                    if (IsSet(sex_male, sample_idx) != (cur_sex_code & 1)) {
                      ClearBit(sample_idx, sex_nm);
                      ClearBit(sample_idx, sex_male);
                      SetBit(sample_idx, sex_locked);
                    }
                  }
                }
              } else {
                // first
                AssignBit(sample_idx, (cur_sex_code + 1) / 2, sex_nm);
                AssignBit(sample_idx, cur_sex_code & 1, sex_male);
                SetBit(sample_idx, sex_locked);
              }
            }
          }
        }
        if (cur_pheno_ct) {
          uintptr_t pheno_uidx_base = 0;
          uintptr_t cur_bits = pheno_include[0];
          for (uint32_t cur_pheno_idx = 0; cur_pheno_idx != cur_pheno_ct; ++cur_pheno_idx) {
            const uintptr_t pheno_uidx = BitIter1(pheno_include, &pheno_uidx_base, &cur_bits);
            PhenoCol* cur_pheno_col = &(pheno_cols[pheno_uidx]);
            if (((!pheno_locked) && IsSet(cur_pheno_col->nonmiss, sample_idx)) || (pheno_locked && IsSet(pheno_locked, sample_idx + pheno_uidx * sample_ct))) {
              continue;
            }
            const uint32_t col_type_idx = pheno_uidx + kNonphenoColCt;
            const char* cur_phenostr = token_ptrs[col_type_idx];
            double dxx;
            const char* cur_phenostr_end = ScanadvDouble(cur_phenostr, &dxx);
            if (!cur_phenostr_end) {
              const uint32_t slen = token_slens[col_type_idx];
              if (IsNanStr(cur_phenostr, slen)) {
                dxx = missing_phenod;
              } else {
                // categorical case
                if (cur_pheno_col->type_code != kPhenoDtypeCat) {
                  if (unlikely(cur_pheno_col->type_code == kPhenoDtypeQt)) {
                    snprintf(g_logbuf, kLogbufSize, "Error: '%s' entry on line %" PRIuPTR " of %s is categorical, while an earlier entry is numeric/'NA'.\n", &(pheno_names[pheno_uidx * max_pheno_name_blen]), line_idx, cur_fname);
                    goto MergePsams_ret_INCONSISTENT_INPUT_WW;
                  }
                  cur_pheno_col->type_code = kPhenoDtypeCat;
                  ZeroU32Arr(sample_ct, cur_pheno_col->data.cat);
                }
                uint32_t catname_idx = IdHtableFindNnt(cur_phenostr, category_names, category_names_htable, slen, category_names_htable_size);
                if (catname_idx == UINT32_MAX) {
                  if (category_names_ct * 2 == category_names_htable_size) {
                    // resize
                    if (unlikely(S_CAST(uintptr_t, arena_top - arena_bottom) < category_names_ct * (2 * sizeof(int32_t) + sizeof(intptr_t)))) {
                      goto MergePsams_ret_NOMEM;
                    }
                    uint32_t next_capacity = category_names_ct * 2;
                    if (next_capacity >= 0x80000000U) {
                      if (unlikely(next_capacity != 0x80000000U)) {
                        logerrputs("Error: This implementation of --pmerge[-list] is limited to ~2^31 distinct\ncategory names.\n");
                        reterr = kPglRetNotYetSupported;
                        goto MergePsams_ret_1;
                      }
                      next_capacity = 0x7fffffff;
                    }
                    category_names_htable = R_CAST(uint32_t*, &(category_names[next_capacity]));
                    category_names_htable_size = 2 * next_capacity;
                    SetAllU32Arr(category_names_htable_size, category_names_htable);
                    arena_bottom = R_CAST(unsigned char*, &(category_names_htable[category_names_htable_size]));
                    for (uint32_t category_idx = 0; category_idx != category_names_ct; ++category_idx) {
                      const char* cur_cat_name = category_names[category_idx];
                      HtableAddNondup(cur_cat_name, strlen(cur_cat_name), category_names_htable_size, category_idx, category_names_htable);
                    }
                  }
                  catname_idx = category_names_ct;
                  if (unlikely(StoreStringAtEndK(arena_bottom, cur_phenostr, slen, &arena_top, &(category_names[catname_idx])))) {
                    goto MergePsams_ret_NOMEM;
                  }
                  ++category_names_ct;
                  HtableAddNondup(cur_phenostr, slen, category_names_htable_size, catname_idx, category_names_htable);
                }
                if (!pheno_locked) {
                  // nm-first and previously missing
                  if (catname_idx) {
                    SetBit(sample_idx, cur_pheno_col->nonmiss);
                    cur_pheno_col->data.cat[sample_idx] = catname_idx;
                  }
                } else {
                  if (merge_pheno_mode == kMergePhenoModeNmMatch) {
                    if (catname_idx) {
                      const uint32_t prev_catname_idx = cur_pheno_col->data.cat[sample_idx];
                      if (prev_catname_idx != catname_idx) {
                        if (!prev_catname_idx) {
                          SetBit(sample_idx, cur_pheno_col->nonmiss);
                          cur_pheno_col->data.cat[sample_idx] = catname_idx;
                        } else {
                          ClearBit(sample_idx, cur_pheno_col->nonmiss);
                          cur_pheno_col->data.cat[sample_idx] = 0;
                          SetBit(sample_idx + pheno_uidx * sample_ct, pheno_locked);
                        }
                      }
                    }
                  } else {
                    // first
                    if (catname_idx) {
                      SetBit(sample_idx, cur_pheno_col->nonmiss);
                      cur_pheno_col->data.cat[sample_idx] = catname_idx;
                    }
                    SetBit(sample_idx + pheno_uidx * sample_ct, pheno_locked);
                  }
                }
                continue;
              }
            } else {
              if (unlikely(!IsSpaceOrEoln(*cur_phenostr_end))) {
                cur_phenostr_end = CurTokenEnd(cur_phenostr_end);
                *K_CAST(char*, cur_phenostr_end) = '\0';
                snprintf(g_logbuf, kLogbufSize, "Error: Invalid numeric token '%s' on line %" PRIuPTR " of %s.\n", cur_phenostr, line_idx, cur_fname);
                goto MergePsams_ret_MALFORMED_INPUT_WW;
              }
            }
            if (unlikely(cur_pheno_col->type_code == kPhenoDtypeCat)) {
              snprintf(g_logbuf, kLogbufSize, "Error: '%s' entry on line %" PRIuPTR " of %s is numeric/'NA', while an earlier entry is categorical.\n", &(pheno_names[pheno_uidx * max_pheno_name_blen]), line_idx, cur_fname);
              goto MergePsams_ret_INCONSISTENT_INPUT_WW;
            }
            cur_pheno_col->type_code = kPhenoDtypeQt;
            if (!pheno_locked) {
              // nm-first
              if (dxx != missing_phenod) {
                if ((dxx == 0.0) && pheno_zero_seen) {
                  SetBit(sample_idx + pheno_uidx * sample_ct, pheno_zero_seen);
                } else {
                  cur_pheno_col->data.qt[sample_idx] = dxx;
                  SetBit(sample_idx, cur_pheno_col->nonmiss);
                }
              }
            } else {
              if (merge_pheno_mode == kMergePhenoModeNmMatch) {
                if (dxx != missing_phenod) {
                  if ((dxx == 0.0) && pheno_zero_seen) {
                    SetBit(sample_idx + pheno_uidx * sample_ct, pheno_zero_seen);
                  } else if (!IsSet(cur_pheno_col->nonmiss, sample_idx)) {
                    SetBit(sample_idx, cur_pheno_col->nonmiss);
                    cur_pheno_col->data.qt[sample_idx] = dxx;
                  } else {
                    if (dxx != cur_pheno_col->data.qt[sample_idx]) {
                      ClearBit(sample_idx, cur_pheno_col->nonmiss);
                      SetBit(sample_idx + pheno_uidx * sample_ct, pheno_locked);
                    }
                  }
                }
              } else {
                if (dxx != missing_phenod) {
                  SetBit(sample_idx, cur_pheno_col->nonmiss);
                  cur_pheno_col->data.qt[sample_idx] = dxx;
                }
                SetBit(sample_idx + pheno_uidx * sample_ct, pheno_locked);
              }
            }
          }
        }
      }
      if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
        goto MergePsams_ret_TSTREAM_FAIL;
      }
      if (unlikely(CleanupTextStream2(cur_fname, &txs, &reterr))) {
        goto MergePsams_ret_1;
      }
      filesets_iter->sample_ct = cur_sample_ct;
      filesets_iter = filesets_iter->next;
    } while (filesets_iter);
    BigstackReset(bigstack_mark2);
    BigstackEndSet(arena_top);
    if (pheno_ct) {
      if (category_names) {
        const char** category_names_final = R_CAST(const char**, g_bigstack_base);
        BigstackBaseSet(&(category_names_final[category_names_ct]));
        memmove(category_names_final, category_names, category_names_ct * sizeof(intptr_t));
        category_names = category_names_final;
      }
      for (uintptr_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        PhenoCol* pheno_col = &(pheno_cols[pheno_idx]);
        if (pheno_col->type_code == kPhenoDtypeQt) {
          uintptr_t* nonmiss = pheno_col->nonmiss;
          double* qt = pheno_col->data.qt;
          const uint32_t sample_nm_ct = PopcountWords(nonmiss, sample_ctl);
          uintptr_t sample_uidx_base = 0;
          uintptr_t cur_bits = nonmiss[0];
          uint32_t sample_idx = 0;
          if (!affection_01) {
            if (merge_pheno_mode == kMergePhenoModeFirst) {
              // 1. Check for values outside {0, 1, 2, missing}.
              // 2. If none exist, convert zeroes to missing.
              for (; sample_idx != sample_nm_ct; ++sample_idx) {
                const uint32_t sample_uidx = BitIter1(nonmiss, &sample_uidx_base, &cur_bits);
                const double dxx = qt[sample_uidx];
                if ((dxx != 0.0) && (dxx != 1.0) && (dxx != 2.0)) {
                  break;
                }
              }
              if (!(sample_idx != sample_nm_ct)) {
                // Could use explicit word-based iteration instead.
                sample_uidx_base = 0;
                cur_bits = nonmiss[0];
                for (sample_idx = 0; sample_idx != sample_nm_ct; ++sample_idx) {
                  const uint32_t sample_uidx = BitIter1(nonmiss, &sample_uidx_base, &cur_bits);
                  if (qt[sample_uidx] == 0.0) {
                    ClearBit(sample_uidx, nonmiss);
                  }
                }
              }
            } else {
              // 1. Check for values outside {1, 2, missing}.
              // 2. If at least one exists, patch in zeroes.  (In nm-match
              //    case, when a nonmissing value is in the cell but a zero was
              //    also observed, convert to missing.)
              for (; sample_idx != sample_nm_ct; ++sample_idx) {
                const uint32_t sample_uidx = BitIter1(nonmiss, &sample_uidx_base, &cur_bits);
                const double dxx = qt[sample_uidx];
                if ((dxx != 1.0) && (dxx != 2.0)) {
                  break;
                }
              }
              if (sample_idx != sample_nm_ct) {
                const uintptr_t bit_offset = pheno_idx * sample_ct;
                const uint32_t zero_ct = PopcountBitRange(pheno_zero_seen, bit_offset, bit_offset + sample_ct);
                if (zero_ct) {
                  uintptr_t zero_seen_uidx_base;
                  BitIter1Start(pheno_zero_seen, bit_offset, &zero_seen_uidx_base, &cur_bits);
                  if (merge_pheno_mode == kMergePhenoModeNmMatch) {
                    for (uint32_t zero_idx = 0; zero_idx != zero_ct; ++zero_idx) {
                      const uintptr_t sample_uidx = BitIter1(pheno_zero_seen, &zero_seen_uidx_base, &cur_bits) - bit_offset;
                      if (IsSet(nonmiss, sample_uidx)) {
                        ClearBit(sample_uidx, nonmiss);
                      } else {
                        SetBit(sample_uidx, nonmiss);
                        qt[sample_uidx] = 0.0;
                      }
                    }
                  } else {
                    // nm-first
                    for (uint32_t zero_idx = 0; zero_idx != zero_ct; ++zero_idx) {
                      const uintptr_t sample_uidx = BitIter1(pheno_zero_seen, &zero_seen_uidx_base, &cur_bits) - bit_offset;
                      SetBit(sample_uidx, nonmiss);
                      qt[sample_uidx] = 0.0;
                    }
                  }
                }
              }
            }
          } else {
            // 1. Check for values outside {0, 1, missing}.
            // 2. If none exist, convert to {1, 2, missing}.
            for (; sample_idx != sample_nm_ct; ++sample_idx) {
              const uint32_t sample_uidx = BitIter1(nonmiss, &sample_uidx_base, &cur_bits);
              const double dxx = qt[sample_uidx];
              if ((dxx != 0.0) && (dxx != 1.0)) {
                break;
              }
            }
            if (!(sample_idx != sample_nm_ct)) {
              sample_uidx_base = 0;
              cur_bits = nonmiss[0];
              for (sample_idx = 0; sample_idx != sample_nm_ct; ++sample_idx) {
                const uint32_t sample_uidx = BitIter1(nonmiss, &sample_uidx_base, &cur_bits);
                qt[sample_uidx] += 1.0;
              }
            }
          }
        } else if (pheno_col->type_code == kPhenoDtypeCat) {
          pheno_col->category_names = category_names;
        }
      }
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    reterr = WritePsam(outname, sample_include, siip, &parental_id_info, sex_nm, sex_male, pheno_cols, pheno_names, nullptr, "NA", sample_ct, pheno_ct, max_pheno_name_blen, kfPsamColDefault);
    if (unlikely(reterr)) {
      goto MergePsams_ret_1;
    }
    logprintfww("--pmerge%s: Merged .psam written to %s .\n", pmip->list_fname? "-list" : "", outname);
    *sample_ctp = sample_ct;
    ArenaEndSet(siip->sample_ids, &bigstack_end_mark);
  }
  while (0) {
  MergePsams_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  MergePsams_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind(cur_fname, &txs, &reterr);
    break;
  MergePsams_ret_TSTREAM_FAIL:
    TextStreamErrPrint(cur_fname, &txs);
    break;
  MergePsams_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, cur_fname);
  MergePsams_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
  MergePsams_ret_MALFORMED_INPUT_2:
    logerrputsb();
  MergePsams_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  MergePsams_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  MergePsams_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 MergePsams_ret_1:
  CleanupTextStream2(cur_fname, &txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

// supports topological sort of chromosome order
typedef struct ChrGraphEdgeLlStruct {
  struct ChrGraphEdgeLlStruct* next;
  ChrIdx dst_idx;
} ChrGraphEdgeLl;

/*
// cip->chr_file_order is filled with the final chromosome sort order.
// info_keys, pointed-to InfoVtype entries, and info_keys_htable are allocated
// at the end of bigstack.
PglErr ScanPvarsAndMergeHeader(const PmergeInfo* pmip, ImportFlags import_flags, uint32_t max_thread_ct, char* outname, char* outname_end, PmergeInputFilesetLl* filesets, ChrInfo* cip, const char* const** info_keys_ptr, uint32_t* info_key_ctp, uint32_t** info_keys_htablep, uint32_t* info_keys_htable_sizep) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  const char* cur_fname = nullptr;
  PglErr reterr = kPglRetSuccess;
  char* cswritep = nullptr;
  uintptr_t line_idx = 0;
  CompressStreamState css;
  PreinitCstream(&css);
  TextStream txs;
  PreinitTextStream(&txs);
  {
    const uintptr_t linebuf_capacity = MINV(kMaxLongLine, bigstack_left() / 4) + kDecompressChunkSize;
    char* linebuf;
    if (unlikely(bigstack_alloc_c(linebuf_capacity, &linebuf))) {
      goto ScanPvarsAndMergeHeader_ret_NOMEM;
    }

    const MergeXheaderMode merge_xheader_mode = pmip->merge_xheader_mode;
    // Each entry is an xheader key, followed by a null terminator, followed by
    // the null-terminated remainder of the header line.  The latter includes
    // the original punctuation character ('=', ',', or '>') ending the key,
    // but not the final newline.
    // Header line classes defined in the VCF v4.3 specification, and the
    // chrSet class defined in the .pvar specification, are recognized.  I.e.
    // when the header line starts with "##INFO=", "##FILTER=", "##FORMAT=",
    // "##ALT=", "##contig=", "##META=", "##SAMPLE=", or "##PEDIGREE=", the key
    // includes everything up to (and not including) the ',' ending the ID,
    // skipping the initial "##".  Otherwise, the key includes everything up to
    // and not including the '=' (also skipping the initial "##").
    const char** xheader_entries = nullptr;
    uint32_t* xheader_entries_htable = nullptr;
    uint32_t xheader_entry_ct = 0;
    // xheader_entry_capacity == xheader_entry_htable_size / 2
    uint32_t xheader_entry_htable_size = 0;
    if (merge_xheader_mode != kMergeXheaderModeErase) {
      xheader_entry_htable_size = 2048;
      if (unlikely(bigstack_alloc_kcp())) {
      }
    }
  }
  while (0) {
  ScanPvarsAndMergeHeader_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
  CleanupTextStream2(cur_fname, &txs, &reterr);
  CswriteCloseCond(&css, cswritep);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}
*/

PglErr Pmerge(const PmergeInfo* pmip, const char* sample_sort_fname, MiscFlags misc_flags, __attribute__((unused)) ImportFlags import_flags, SortMode sample_sort_mode, FamCol fam_cols, int32_t missing_pheno, uint32_t max_thread_ct, char* pgenname, char* psamname, char* pvarname, char* outname, char* outname_end, __attribute__((unused)) ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    // 1. Construct/load fileset list.
    // 2. Merge .psam files.
    // 3. Global .pvar scan, to determine merged header, chromosome sort order,
    //    and track first/last position in each fileset.
    // 4. If filesets cover disjoint positions, handle this as a concatenation
    //    job (or error out on --variant-inner-join).
    // 5. Otherwise, perform general-purpose incremental merge.
    PmergeInputFilesetLl* filesets = nullptr;
    PmergeInputFilesetLl** filesets_endp = &filesets;
    if (pgenname[0]) {
      PmergeInputFilesetLl* cur_entry = AllocFilesetLlEntry(&filesets_endp);
      if (unlikely(!cur_entry)) {
        goto Pmerge_ret_NOMEM;
      }
      const uint32_t pgen_fname_blen = strlen(pgenname) + 1;
      const uint32_t pvar_fname_blen = strlen(pvarname) + 1;
      const uint32_t psam_fname_blen = strlen(psamname) + 1;
      char* fname_iter;
      if (unlikely(bigstack_end_alloc_c(pgen_fname_blen + pvar_fname_blen + psam_fname_blen, &fname_iter))) {
        goto Pmerge_ret_NOMEM;
      }
      cur_entry->pgen_fname = fname_iter;
      fname_iter = memcpya(fname_iter, pgenname, pgen_fname_blen);
      cur_entry->pvar_fname = fname_iter;
      fname_iter = memcpya(fname_iter, pvarname, pvar_fname_blen);
      cur_entry->psam_fname = fname_iter;
      memcpy(fname_iter, psamname, psam_fname_blen);
      cur_entry->pgen_locked_fname = nullptr;
    }
    if (!pmip->list_fname) {
      PmergeInputFilesetLl* cur_entry = AllocFilesetLlEntry(&filesets_endp);
      if (unlikely(!cur_entry)) {
        goto Pmerge_ret_NOMEM;
      }
      cur_entry->pgen_fname = pmip->pgen_fname;
      cur_entry->pvar_fname = pmip->pvar_fname;
      cur_entry->psam_fname = pmip->psam_fname;
      cur_entry->pgen_locked_fname = nullptr;
    } else {
      reterr = LoadPmergeList(pmip->list_fname, pmip->list_mode, pgenname[0] != '\0', &filesets_endp);
      if (unlikely(reterr)) {
        goto Pmerge_ret_1;
      }
    }

    SampleIdInfo sii;
    uint32_t sample_ct;
    reterr = MergePsams(pmip, sample_sort_fname, misc_flags, sample_sort_mode, fam_cols, missing_pheno, max_thread_ct, outname, outname_end, filesets, &sii, &sample_ct);
    if (unlikely(reterr)) {
      goto Pmerge_ret_1;
    }
    /*
    const char* const* info_keys = nullptr;
    uint32_t* info_keys_htable = nullptr;
    uint32_t info_key_ct = 0;
    uint32_t info_keys_htable_size = 0;
    reterr = ScanPvarsAndMergeHeader(pmip, import_flags, max_thread_ct, outname, outname_end, filesets, cip, &info_keys, &info_key_ct, &info_keys_htable, &info_keys_htable_size);
    if (unlikely(reterr)) {
      goto Pmerge_ret_1;
    }
    */
    logerrputs("Error: --pmerge[-list] is under development.\n");
    reterr = kPglRetNotYetSupported;
    goto Pmerge_ret_1;
  }
  while (0) {
  Pmerge_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 Pmerge_ret_1:
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
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

static_assert(sizeof(Dosage) == 2, "PgenDiff() must be updated.");
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
          const uint32_t cur_sex_code = CharToSex(token_iter[0]);
          if ((token_slen != 1) || (!cur_sex_code)) {
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

        // We verify during .pvar load that, if a missing allele code is
        // present, it's in a biallelic variant, and the other allele code is
        // nonmissing.
        // We verify below that the missing allele code does not appear in the
        // actual genotype calls.
        uint32_t missing1_idx = UINT32_MAX;
        for (uint32_t allele1_idx = 0; allele1_idx != cur_allele_ct1; ++allele1_idx) {
          const char* cur_allele1 = cur_allele1s[allele1_idx];
          if (memequal_k(cur_allele1, ".", 2)) {
            assert(missing1_idx == UINT32_MAX);
            missing1_idx = allele1_idx;
          }
        }
        ZeroTrailingNyps(sample_ct, genovec1);
        ZeroTrailingNyps(sample_ct, genovec2);
        if (missing1_idx != UINT32_MAX) {
          assert(cur_allele_ct1 == 2);
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
        // See e.g. PopulateStrboxHtable().
        uint32_t hashval = Hashceil(merged_alleles[0], strlen(merged_alleles[0]), cur_htable_size);
        merged_alleles_htable[hashval] = 0;
        ++merged_allele_ct;

        for (uint32_t allele1_idx = allele1_alt_start; allele1_idx != cur_allele_ct1; ++allele1_idx) {
          if (allele1_idx == missing1_idx) {
            continue;
          }
          const char* cur_allele1 = cur_allele1s[allele1_idx];
          const uint32_t cur_allele1_slen = strlen(cur_allele1);
          const uint32_t cur_idx = IdHtableAdd(cur_allele1, merged_alleles, cur_allele1_slen, cur_htable_size, merged_allele_ct, merged_alleles_htable);
          if (cur_idx == UINT32_MAX) {
            if (unlikely(merged_allele_ct == kPglMaxAlleleCt)) {
              logerrprintfww("Error: Too many alleles across --pgen-diff and main filesets for variant '%s' at position %s:%u. (This " PROG_NAME_STR " build is limited to " PGL_MAX_ALLELE_CT_STR ".)\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
              reterr = kPglRetNotYetSupported;
              goto PgenDiff_ret_1;
            }
            remap1[allele1_idx] = merged_allele_ct;
            merged_alleles[merged_allele_ct] = cur_allele1;
            ++merged_allele_ct;
          } else {
            remap1[allele1_idx] = cur_idx;
          }
        }
        for (uint32_t allele2_idx = allele2_alt_start; allele2_idx != cur_allele_ct2; ++allele2_idx) {
          if (allele2_idx == missing2_idx) {
            continue;
          }
          const char* cur_allele2 = cur_allele2s[allele2_idx];
          const uint32_t cur_allele2_slen = strlen(cur_allele2);
          const uint32_t cur_idx = IdHtableAdd(cur_allele2, merged_alleles, cur_allele2_slen, cur_htable_size, merged_allele_ct, merged_alleles_htable);
          if (cur_idx == UINT32_MAX) {
            if (unlikely(merged_allele_ct == kPglMaxAlleleCt)) {
              logerrprintfww("Error: Too many alleles across --pgen-diff and main filesets for variant '%s' at position %s:%u. (This " PROG_NAME_STR " build is limited to " PGL_MAX_ALLELE_CT_STR ".)\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
              reterr = kPglRetNotYetSupported;
              goto PgenDiff_ret_1;
            }
            remap2[allele2_idx] = merged_allele_ct;
            merged_alleles[merged_allele_ct] = cur_allele2;
            ++merged_allele_ct;
          } else {
            remap2[allele2_idx] = cur_idx;
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
