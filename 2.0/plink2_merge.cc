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

#include "plink2_merge.h"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "include/pgenlib_misc.h"
#include "include/pgenlib_write.h"
#include "include/plink2_bgzf.h"
#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_memory.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "plink2_cmdline.h"
#include "plink2_compress_stream.h"
#include "plink2_decompress.h"
#include "plink2_data.h"
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
  pmerge_info_ptr->merge_qual_mode = kMergeQualModeMin;
  pmerge_info_ptr->merge_filter_mode = kMergeFilterModeNmFirst;
  pmerge_info_ptr->merge_info_mode = kMergeInfoCmModeNmFirst;
  pmerge_info_ptr->merge_cm_mode = kMergeInfoCmModeNmFirst;
  pmerge_info_ptr->merge_pheno_sort = kSortNone;
  pmerge_info_ptr->merge_info_sort = kSortNone;
  pmerge_info_ptr->max_allele_ct = 0;
  pmerge_info_ptr->pgen_fname = nullptr;
  pmerge_info_ptr->pvar_fname = nullptr;
  pmerge_info_ptr->psam_fname = nullptr;
  pmerge_info_ptr->list_fname = nullptr;
  pmerge_info_ptr->list_base_dir = nullptr;
}

void CleanupPmerge(PmergeInfo* pmerge_info_ptr) {
  free_cond(pmerge_info_ptr->pgen_fname);
  free_cond(pmerge_info_ptr->pvar_fname);
  free_cond(pmerge_info_ptr->psam_fname);
  free_cond(pmerge_info_ptr->list_fname);
  free_cond(pmerge_info_ptr->list_base_dir);
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
  // Heap-allocated iff temporary
  char* pgen_fname;
  char* pvar_fname;
  char* psam_fname;
  char* pgen_locked_fname;

  uint32_t read_sample_ct;
  uint32_t write_sample_ct;
  uint32_t read_variant_ct;
  // accounts for --chr, negative bp, same-ID + same-position
  uint32_t write_variant_ct;
  // Also accounts for --merge-max-alleles.  Separate from write_variant_ct
  // since, for large --pmerge-list jobs, we may need to write a tombstone for
  // variants doomed to fail the --merge-max-alleles filter.
  // We compute this upfront to make PmergeConcat() more efficient, since
  // that's the most common use case.
  uint32_t write_nondoomed_variant_ct;
  uint32_t max_pvar_line_blen;
  uint32_t max_single_pos_ct;
  uintptr_t max_single_pos_blen; // includes chrom/pos, could remove later
  char* first_varid;  // heap-allocated
  char* last_varid;  // heap-allocated
  uint32_t first_pos;
  uint32_t last_pos;
  ChrIdx first_chr_idx;
  ChrIdx last_chr_idx;
  uint32_t read_max_nonpass_filter_ct;
  AlleleCode read_max_allele_ct;
  AlleleCode write_nondoomed_max_allele_ct;
  unsigned char nm_qual_exists;
  unsigned char nm_filter_exists;
  unsigned char nm_info_exists;
  unsigned char pvar_info_pr_exists;
  unsigned char nz_cm_present;

  // initialized by ScanPgenHeaders()
  unsigned char nonref_flags_storage;
  unsigned char vrtype_8bit_needed;

  // Note that this struct is allocated on bigstack, rather than the heap, iff
  // the fileset is not temporary.
  unsigned char is_temporary;
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

PglErr LoadPmergeList(const char* list_fname, const char* list_base_dir, PmergeListMode mode, uint32_t main_fileset_present, PmergeInputFilesetLl*** filesets_endpp, uintptr_t* fileset_ctp) {
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
    uint32_t base_dir_slen = 0;
    if (list_base_dir) {
      base_dir_slen = strlen(list_base_dir);
    }
    const uint32_t max_single_token_slen = kPglFnamesize - 1 - MAXV(pgen_blen, pvar_blen) - base_dir_slen;
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
      cur_entry->is_temporary = 0;
      cur_entry->pgen_locked_fname = nullptr;
      cur_entry->first_varid = nullptr;
      cur_entry->last_varid = nullptr;
      const char* first_token_end = CurTokenEnd(first_token_start);
      const uint32_t first_token_slen = first_token_end - first_token_start;
      const char* second_token_start = FirstNonTspace(first_token_end);
      if (IsEolnKns(*second_token_start)) {
        // Expand single token, using the --pmerge-list mode.
        if (unlikely(first_token_slen > max_single_token_slen)) {
          logerrprintf("Error: Entry on line %" PRIuPTR " of --pmerge-list file is too long.\n", line_idx);
          goto LoadPmergeList_ret_INCONSISTENT_INPUT;
        }
        char* fname_iter;
        if (unlikely(bigstack_end_alloc_c((first_token_slen + base_dir_slen) * 3 + pgen_blen + pvar_blen + psam_blen, &fname_iter))) {
          goto LoadPmergeList_ret_NOMEM;
        }
        cur_entry->pgen_fname = fname_iter;
        if (base_dir_slen) {
          fname_iter = memcpya(fname_iter, list_base_dir, base_dir_slen);
        }
        fname_iter = memcpya(fname_iter, first_token_start, first_token_slen);
        fname_iter = memcpya(fname_iter, pgen_suffix, pgen_blen);
        cur_entry->pvar_fname = fname_iter;
        if (base_dir_slen) {
          fname_iter = memcpya(fname_iter, list_base_dir, base_dir_slen);
        }
        fname_iter = memcpya(fname_iter, first_token_start, first_token_slen);
        fname_iter = memcpya(fname_iter, pvar_suffix, pvar_blen);
        cur_entry->psam_fname = fname_iter;
        if (base_dir_slen) {
          fname_iter = memcpya(fname_iter, list_base_dir, base_dir_slen);
        }
        fname_iter = memcpya(fname_iter, first_token_start, first_token_slen);
        memcpy(fname_iter, psam_suffix, psam_blen);
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
          goto LoadPmergeList_ret_INCONSISTENT_INPUT;
        }
        char* fname_iter;
        if (unlikely(bigstack_end_alloc_c(base_dir_slen * 3 + first_token_slen + second_token_slen + third_token_slen + 3, &fname_iter))) {
          goto LoadPmergeList_ret_NOMEM;
        }
        cur_entry->pgen_fname = fname_iter;
        if (base_dir_slen) {
          fname_iter = memcpya(fname_iter, list_base_dir, base_dir_slen);
        }
        fname_iter = memcpyax(fname_iter, first_token_start, first_token_slen, '\0');
        cur_entry->pvar_fname = fname_iter;
        if (base_dir_slen) {
          fname_iter = memcpya(fname_iter, list_base_dir, base_dir_slen);
        }
        fname_iter = memcpyax(fname_iter, second_token_start, second_token_slen, '\0');
        cur_entry->psam_fname = fname_iter;
        if (base_dir_slen) {
          fname_iter = memcpya(fname_iter, list_base_dir, base_dir_slen);
        }
        memcpyx(fname_iter, third_token_start, third_token_slen, '\0');
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto LoadPmergeList_ret_TSTREAM_FAIL;
    }
    const uintptr_t fileset_ct = main_fileset_present + line_idx;
    // update (26 Jul 2025): permit just one fileset, since it makes sense for
    // e.g. --merge-sids
    if (unlikely(fileset_ct == 0)) {
      logerrputs("Error: --pmerge-list requires at least one fileset to be specified.\n");
      goto LoadPmergeList_ret_INCONSISTENT_INPUT;
    }
    *fileset_ctp = fileset_ct;
    logprintf("--pmerge-list: %" PRIuPTR " fileset%s specified%s.\n", fileset_ct, (fileset_ct == 1)? "" : "s", main_fileset_present? " (including main fileset)" : "");
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
    reterr = kPglRetMalformedInput;
    break;
  LoadPmergeList_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  CleanupTextStream2(list_fname, &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

// Permanent allocations are at end of bigstack.
PglErr MergePsams(const PmergeInfo* pmip, const char* sample_sort_fname, const char* missing_catname, MiscFlags misc_flags, SortMode sample_sort_mode, FamCol fam_cols, int32_t missing_pheno, uint32_t max_thread_ct, char* outname, char* outname_end, PmergeInputFilesetLl* filesets, SampleIdInfo* siip, uint32_t* sample_ctp, uint32_t* linebuf_capacityp) {
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
    if (flags & kfPmergeSids) {
      logerrputs("Error: --merge-sids is under development.\n");
      reterr = kPglRetNotYetSupported;
      goto MergePsams_ret_1;
    }
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
            if (unlikely(memequal_sk(token_start, "FID"))) {
              if (fid_present) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate FID column in %s.\n", cur_fname);
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Improperly positioned FID column in %s (must be first).\n", cur_fname);
              }
              goto MergePsams_ret_MALFORMED_INPUT_WW;
            }
            if (unlikely(memequal_sk(token_start, "IID"))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Duplicate IID column in %s.\n", cur_fname);
              goto MergePsams_ret_MALFORMED_INPUT_WW;
            }
            if (unlikely(memequal_sk(token_start, "SID"))) {
              if (sid_present) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate SID column in %s.\n", cur_fname);
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: Improperly positioned SID column in %s (must immediately follow IID).\n", cur_fname);
              }
              goto MergePsams_ret_MALFORMED_INPUT_WW;
            }
            if (memequal_sk(token_start, "PAT")) {
              if (unlikely(postid_pat_col_idx)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate PAT column in %s.\n", cur_fname);
                goto MergePsams_ret_MALFORMED_INPUT_WW;
              }
              postid_pat_col_idx = postid_col_idx;
              continue;
            }
            if (memequal_sk(token_start, "MAT")) {
              if (unlikely(postid_mat_col_idx)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Duplicate MAT column in %s.\n", cur_fname);
                goto MergePsams_ret_MALFORMED_INPUT_WW;
              }
              postid_mat_col_idx = postid_col_idx;
              continue;
            }
            if (memequal_sk(token_start, "SEX")) {
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
          token_end = CurTokenEnd(sid_start);
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
            if (sample_ct == kPglMaxSampleCt) {
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
            if (unlikely(first_sample_id_ct == kPglMaxSampleCt)) {
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
    // 2. (TODO: --merge-sids, once general-case merge is implemented.)
    // 3. Compute max_sample_id_blen and max_sid_blen if --sample-inner-join.
    // 4. Free all end-of-bigstack allocations.
    // 5. Construct SampleIdInfo at end of bigstack, to be returned.
    // 6. Construct pheno_names, hash tables.
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
      if ((sample_sort_mode != kSortNone) && (sample_ct > 1)) {
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
      InitXidHtable(siip, sample_ct, sample_id_table_size, sample_id_htable, idbuf);
    }
    linebuf_capacity = MAXV(max_line_blen, kDecompressMinBlen) + kDecompressChunkSize;
    *linebuf_capacityp = linebuf_capacity;
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
    const uint32_t kNonphenoPostIdColCt = 3;
    const uint32_t max_col_ct = pheno_ct + kNonphenoPostIdColCt;
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
    // bugfix (22 Mar 2022): We've only checked for duplicate sample IDs in the
    // --sample-inner-join case.  The current implementation does not support
    // them, so we must error out in the non-sample-inner-join case too.
    cur_sample_include = nullptr;
    if (!(flags & kfPmergeSampleInnerJoin)) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &cur_sample_include))) {
        goto MergePsams_ret_NOMEM;
      }
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
      logprintf("--pmerge%s: %u sample%s present.\n", pmip->list_fname? "-list" : "", sample_ct, (sample_ct == 1)? "" : "s");
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
      const uint32_t missing_catname_slen = strlen(missing_catname);
      if (unlikely(StoreStringAtEndK(g_bigstack_base, missing_catname, missing_catname_slen, &arena_top, &(category_names[0])))) {
        goto MergePsams_ret_NOMEM;
      }
      const uint32_t hashval = Hashceil(missing_catname, missing_catname_slen, category_names_htable_size);
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
      if (cur_sample_include) {
        ZeroWArr(sample_ctl, cur_sample_include);
      }
      ZeroWArr(pheno_ctl, pheno_include);
      uint32_t sid_present = 0;
      uint32_t parents_present = 0;
      uint32_t sex_present = 0;
      uint32_t relevant_postid_col_ct = 0;
      uint32_t fid_present;
      col_types[0] = 0;
      if (line_start[0] == '#') {
        fid_present = (line_start[1] == 'F');
        uint32_t col_idx = 1;
        const char* token_iter = FirstNonTspace(&(line_start[4]));
        if (fid_present) {
          token_iter = FirstNonTspace(CurTokenEnd(token_iter));
          ++col_idx;
        }
        sid_present = tokequal_k(token_iter, "SID");
        if (sid_present) {
          token_iter = FirstNonTspace(CurTokenEnd(token_iter));
          ++col_idx;
        }
        for (; ; ++col_idx) {
          if (IsEolnKns(*token_iter)) {
            break;
          }
          const char* token_end = CurTokenEnd(token_iter);
          const uint32_t token_slen = token_end - token_iter;
          uint32_t cur_col_type = UINT32_MAX;
          if (token_slen == 3) {
            if (memequal_sk(token_iter, "PAT")) {
              cur_col_type = 0;
              parents_present = 1;
            } else if (memequal_sk(token_iter, "MAT")) {
              cur_col_type = 1;
            } else if (memequal_sk(token_iter, "SEX")) {
              cur_col_type = 2;
              sex_present = 1;
            }
          }
          const char* cur_token_start = token_iter;
          token_iter = FirstNonTspace(token_end);
          if (cur_col_type == UINT32_MAX) {
            if (token_slen < max_pheno_name_blen) {
              cur_col_type = StrboxHtableFindNnt(cur_token_start, pheno_names, pheno_names_htable, max_pheno_name_blen, token_slen, pheno_names_table_size);
            }
            if (cur_col_type == UINT32_MAX) {
              continue;
            }
            SetBit(cur_col_type, pheno_include);
            cur_col_type += kNonphenoPostIdColCt;
          }
          col_skips[relevant_postid_col_ct] = col_idx;
          col_types[relevant_postid_col_ct++] = cur_col_type;
        }
        if (relevant_postid_col_ct) {
          for (uint32_t rp_col_idx = relevant_postid_col_ct - 1; rp_col_idx; --rp_col_idx) {
            col_skips[rp_col_idx] -= col_skips[rp_col_idx - 1];
          }
        }
        line_start = AdvPastDelim(token_iter, '\n');
        ++line_idx;
      } else {
        // .fam
        fid_present = (fam_cols / kfFamCol1) & 1;
        uint32_t prev_col_idx = 0;
        uint32_t col_idx = fid_present + 1;
        if (fam_cols & kfFamCol34) {
          col_skips[relevant_postid_col_ct] = col_idx - prev_col_idx;
          col_types[relevant_postid_col_ct++] = 0;
          col_skips[relevant_postid_col_ct] = 1;
          col_types[relevant_postid_col_ct++] = 1;
          parents_present = 1;
          prev_col_idx = col_idx + 1;
          col_idx += 2;
        }
        if (fam_cols & kfFamCol5) {
          col_skips[relevant_postid_col_ct] = col_idx - prev_col_idx;
          col_types[relevant_postid_col_ct++] = 2;
          sex_present = 1;
          prev_col_idx = col_idx;
          ++col_idx;
        }
        if (fam_cols & kfFamCol6) {
          col_skips[relevant_postid_col_ct] = col_idx - prev_col_idx;
          const uint32_t pheno_idx = StrboxHtableFind("PHENO1", pheno_names, pheno_names_htable, max_pheno_name_blen, strlen("PHENO1"), pheno_names_table_size);
          col_types[relevant_postid_col_ct++] = kNonphenoPostIdColCt + pheno_idx;
          SetBit(pheno_idx, pheno_include);
        }
      }

      const uintptr_t line_idx_body_start = line_idx;
      const uint32_t cur_pheno_ct = PopcountWords(pheno_include, pheno_ctl);
      uint32_t write_sample_ct = 0;
      const char* line_iter;
      for (; TextGetUnsafe2K(&txs, &line_start); line_start = AdvPastDelim(line_iter, '\n'), ++line_idx) {
        line_iter = line_start;
        if (relevant_postid_col_ct) {
          line_iter = TokenLexK(line_start, col_types, col_skips, relevant_postid_col_ct, token_ptrs, token_slens);
          if (unlikely(!line_iter)) {
            goto MergePsams_ret_MISSING_TOKENS;
          }
        }
        uint32_t sample_idx;
        if (unlikely(LookupXidHtable(line_start, siip, sample_id_htable, sample_id_table_size, fid_present, sid_present, &sample_idx, idbuf))) {
          goto MergePsams_ret_REWIND_FAIL;
        }
        if (sample_idx == UINT32_MAX) {
          continue;
        }
        if (cur_sample_include) {
          if (IsSet(cur_sample_include, sample_idx)) {
            TabsToSpaces(idbuf);
            snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID \"%s\" in %s.\n", idbuf, cur_fname);
            goto MergePsams_ret_MALFORMED_INPUT_WW;
          }
          SetBit(sample_idx, cur_sample_include);
        }
        ++write_sample_ct;
        if (parents_present) {
          char* dad_id_dst = &(parental_id_info.paternal_ids[sample_idx * max_paternal_id_blen]);
          char* mom_id_dst = &(parental_id_info.maternal_ids[sample_idx * max_maternal_id_blen]);
          if (parents_locked) {
            if (!IsSet(parents_locked, sample_idx * 2)) {
              const char* dad_id = token_ptrs[0];
              const uint32_t dad_slen = token_slens[0];
              if (merge_parents_mode == kMergePhenoModeNmMatch) {
                if (((dad_slen != 1) || (dad_id[0] != '0')) && (!strequal_unsafe(dad_id_dst, dad_id, dad_slen))) {
                  if (strequal_k_unsafe(dad_id_dst, "0")) {
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
              const char* mom_id = token_ptrs[1];
              const uint32_t mom_slen = token_slens[1];
              if (merge_parents_mode == kMergePhenoModeNmMatch) {
                if (((mom_slen != 1) || (mom_id[0] != '0')) && (!strequal_unsafe(mom_id_dst, mom_id, mom_slen))) {
                  if (strequal_k_unsafe(mom_id_dst, "0")) {
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
            if (strequal_k_unsafe(dad_id_dst, "0")) {
              memcpyx(dad_id_dst, token_ptrs[0], token_slens[0], '\0');
            }
            if (strequal_k_unsafe(mom_id_dst, "0")) {
              memcpyx(mom_id_dst, token_ptrs[1], token_slens[1], '\0');
            }
          }
        }
        if (sex_present) {
          if (!sex_locked) {
            if (!IsSet(sex_nm, sample_idx)) {
              // nm-first and previously missing
              const uint32_t cur_sex_code = (token_slens[2] == 1)? CharToSex(token_ptrs[2][0]) : 0;
              if (cur_sex_code) {
                SetBit(sample_idx, sex_nm);
                AssignBit(sample_idx, cur_sex_code & 1, sex_male);
              }
            }
          } else {
            if (!IsSet(sex_locked, sample_idx)) {
              const uint32_t cur_sex_code = (token_slens[2] == 1)? CharToSex(token_ptrs[2][0]) : 0;
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
            const uint32_t col_type_idx = pheno_uidx + kNonphenoPostIdColCt;
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
      filesets_iter->read_sample_ct = line_idx - line_idx_body_start;
      filesets_iter->write_sample_ct = write_sample_ct;
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
    reterr = WritePsam(outname, sample_include, siip, &parental_id_info, sex_nm, sex_male, pheno_cols, pheno_names, nullptr, "NA", sample_ct, pheno_ct, max_pheno_name_blen, kfPsamColDefault, 0);
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
  MergePsams_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, cur_fname);
    reterr = kPglRetRewindFail;
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

// This executes before ScanPvarsAndMergeHeader(), so we know whether to write
// an INFO/PR header line even when it doesn't appear in any input .pvar.
PglErr ScanPgenHeaders(uint32_t is_list, MiscFlags misc_flags, PmergeInputFilesetLl* filesets) {
  const char* read_pgen_fname = nullptr;
  PglErr reterr = kPglRetSuccess;
  PgenFileInfo pgfi;
  PreinitPgfi(&pgfi);
  {
    const uint32_t real_ref_alleles = (misc_flags / kfMiscRealRefAlleles) & 1;
    PmergeInputFilesetLl* filesets_iter = filesets;
    do {
      read_pgen_fname = filesets_iter->pgen_fname;
      PgenHeaderCtrl header_ctrl;
      uintptr_t cur_alloc_cacheline_ct;  // unused
      reterr = PgfiInitPhase1(read_pgen_fname, nullptr, UINT32_MAX, filesets_iter->read_sample_ct, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, g_logbuf);
      if (unlikely(reterr)) {
        if (reterr == kPglRetSampleMajorBed) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s is a sample-major .bed file; this is not supported by --pmerge%s. Retry after converting it to a .pgen.\n", read_pgen_fname, is_list? "-list" : "");
          goto ScanPgenHeaders_ret_INCONSISTENT_INPUT_WW;
        }
        WordWrapB(0);
        logerrputsb();
        goto ScanPgenHeaders_ret_1;
      }
      uint32_t vrtype_8bit_needed = 0;
      uint32_t nonref_flags_storage;
      if (pgfi.const_vrtype == kPglVrtypePlink1) {
        nonref_flags_storage = 2 - real_ref_alleles;
      } else {
        nonref_flags_storage = header_ctrl >> 6;
        if (pgfi.const_vrtype == UINT32_MAX) {
          if (((header_ctrl & 12) == 4) || ((header_ctrl & 15) > 9)) {
            vrtype_8bit_needed = 1;
          }
        } else if (pgfi.const_vrtype > 15) {
          vrtype_8bit_needed = 1;
        }
      }
      filesets_iter->nonref_flags_storage = nonref_flags_storage;
      filesets_iter->vrtype_8bit_needed = vrtype_8bit_needed;
      if (unlikely(CleanupPgfi2(read_pgen_fname, &pgfi, &reterr))) {
        goto ScanPgenHeaders_ret_1;
      }
      filesets_iter = filesets_iter->next;
    } while (filesets_iter);
  }
  while (0) {
  ScanPgenHeaders_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 ScanPgenHeaders_ret_1:
  CleanupPgfi2(read_pgen_fname, &pgfi, &reterr);
  return reterr;
}

typedef struct RescanOnePosRecordStruct {
  uint32_t rec_blen;
  AlleleCode allele_ct;
  char variant_id[];  // null-terminated, followed by null-terminated REF, ALT
} RescanOnePosRecord;

typedef struct RescanOnePosContextStruct {
  RescanOnePosRecord* first_record;
  uint32_t write_allele_ct_max;
  uint32_t sort_vars_ascii;
  uint32_t multiallelics_already_joined;
  char input_missing_geno_char;
  uint32_t write_doomed_variant_ct;
  uint32_t write_nondoomed_max_allele_ct;
  uint32_t first_chr_idx;
  uint32_t first_bp;
  char** first_varid_ptr;
  const char* cur_fname;
  const ChrInfo* cip;
} RescanOnePosContext;

PglErr RescanOnePos(unsigned char* arena_top, uint32_t batch_size, uint32_t prev_chr_code, uint32_t prev_bp, unsigned char* arena_bottom, RescanOnePosContext* ctxp, uint32_t* nonwrite_variant_ctp) {
  char* first_varid;
  if (batch_size == 1) {
    // Fast path for common case.
    RescanOnePosRecord* first_record = ctxp->first_record;
    const uint32_t prev_allele_ct = first_record->allele_ct;
    if (prev_allele_ct != 2) {
      if (prev_allele_ct > ctxp->write_allele_ct_max) {
        ctxp->write_doomed_variant_ct += 1;
      } else if (prev_allele_ct > ctxp->write_nondoomed_max_allele_ct) {
        ctxp->write_nondoomed_max_allele_ct = prev_allele_ct;
      }
    }
    if (ctxp->first_bp != UINT32_MAX) {
      return kPglRetSuccess;
    }
    first_varid = first_record->variant_id;
  } else {
    if (!batch_size) {
      return kPglRetSuccess;
    }
    // Sort variant IDs in ASCII or natural order depending on
    // --sort-vars setting, then count number of alleles in merged
    // variant when duplicates are present.
    const uintptr_t bytes_to_round_up = (-R_CAST(uintptr_t, arena_bottom)) % sizeof(intptr_t);
    const uintptr_t extra_bytes_needed = bytes_to_round_up + batch_size * sizeof(intptr_t);
    if (unlikely(S_CAST(uintptr_t, arena_top - arena_bottom) < extra_bytes_needed)) {
      return kPglRetNomem;
    }
    char** sorted_variant_ids = R_CAST(char**, &(arena_bottom[bytes_to_round_up]));
    arena_bottom = &(arena_bottom[extra_bytes_needed]);
    RescanOnePosRecord* record_iter = ctxp->first_record;
    for (uint32_t uii = 0; uii != batch_size; ++uii) {
      const uintptr_t rec_blen = record_iter->rec_blen;
      sorted_variant_ids[uii] = record_iter->variant_id;
      record_iter = R_CAST(RescanOnePosRecord*, &(R_CAST(unsigned char*, record_iter)[rec_blen]));
    }
    if (ctxp->sort_vars_ascii) {
      StrptrArrSortOverread(batch_size, K_CAST(const char**, sorted_variant_ids));
    } else {
      StrptrArrNsort(batch_size, K_CAST(const char**, sorted_variant_ids));
    }
    char* cur_variant_id = sorted_variant_ids[0];
    uint32_t variant_idx_start = 0;
    uint32_t next_allele_ct = container_of(cur_variant_id, RescanOnePosRecord, variant_id)->allele_ct;
    uintptr_t allele_ct_limit = next_allele_ct;
    for (uint32_t variant_idx_end = 1; ; ++variant_idx_end) {
      if (!(variant_idx_end == batch_size)) {
        char* next_variant_id = sorted_variant_ids[variant_idx_end];
        next_allele_ct = container_of(next_variant_id, RescanOnePosRecord, variant_id)->allele_ct;
        if (strequal_overread(cur_variant_id, next_variant_id)) {
          allele_ct_limit += next_allele_ct;
          continue;
        }
      }
      uint32_t merged_allele_ct = allele_ct_limit;
      const uint32_t extra_read_variant_ct = variant_idx_end - variant_idx_start - 1;
      if (extra_read_variant_ct) {
        *nonwrite_variant_ctp += extra_read_variant_ct;

        // Count the number of distinct alleles across this group of
        // same-position, same-ID variants by creating a strptr array and
        // sorting it.  (Could try a char** hash table instead if this is ever
        // a bottleneck, but that seems exceedingly unlikely.)
        unsigned char* arena_bottom_tmp = arena_bottom;
        const char** alleles;
        if (unlikely(arena_alloc_kcp(arena_top, allele_ct_limit, &arena_bottom_tmp, &alleles))) {
          return kPglRetNomem;
        }
        const uint32_t variant_id_blen = strlen(cur_variant_id) + 1;
        const char input_missing_geno_char = ctxp->input_missing_geno_char;
        const char** alleles_iter = alleles;
        uintptr_t missing_allele_ct = 0;
        for (uint32_t variant_idx = variant_idx_start; variant_idx != variant_idx_end; ++variant_idx) {
          char* tmp_variant_id = sorted_variant_ids[variant_idx];
          char* ref_allele = &(tmp_variant_id[variant_id_blen]);
          if (((ref_allele[0] == '.') || (ref_allele[0] == input_missing_geno_char)) && (ref_allele[1] == '\0')) {
            ++missing_allele_ct;
          } else {
            *alleles_iter++ = ref_allele;
          }

          const uint32_t extra_alt_ct = container_of(tmp_variant_id, RescanOnePosRecord, variant_id)->allele_ct - 2;
          char* alt_iter = &(strnul(ref_allele)[1]);
          for (uint32_t extra_alt_idx = 0; extra_alt_idx != extra_alt_ct; ++extra_alt_idx) {
            char* cur_alt_end = AdvToDelim(alt_iter, ',');
            *cur_alt_end = '\0';
            // Missing allele prohibited here.
            *alleles_iter++ = alt_iter;
            alt_iter = &(cur_alt_end[1]);
          }
          if (((alt_iter[0] == '.') || (alt_iter[0] == input_missing_geno_char)) && (alt_iter[1] == '\0')) {
            ++missing_allele_ct;
          } else {
            *alleles_iter++ = alt_iter;
          }
        }
        allele_ct_limit -= missing_allele_ct;
        StrptrArrSortOverread(allele_ct_limit, alleles);
        merged_allele_ct = 0;
        if (allele_ct_limit) {
          merged_allele_ct = 1;
          for (uintptr_t ulii = 1; ulii != allele_ct_limit; ++ulii) {
            if (!strequal_overread(alleles[ulii], alleles[ulii - 1])) {
              ++merged_allele_ct;
            }
          }
        }
        // bugfix (8 Sep 2021): incorrect to perform this check if variant is
        // being filtered out by --merge-max-alleles setting
        if ((!ctxp->multiallelics_already_joined) && (merged_allele_ct <= ctxp->write_allele_ct_max)) {
          const uint32_t merged_variant_ct = variant_idx_end - variant_idx_start;
          if ((allele_ct_limit == merged_variant_ct * 2) && (merged_allele_ct == merged_variant_ct + 1) && (!missing_allele_ct)) {
            // If all REF alleles are equal, sound the alarm.
            char* first_ref_allele = &(sorted_variant_ids[variant_idx_start][variant_id_blen]);
            const uint32_t ref_blen = strlen(first_ref_allele) + 1;
            uint32_t variant_idx = variant_idx_start + 1;
            for (; variant_idx != variant_idx_end; ++variant_idx) {
              char* ref_allele = &(sorted_variant_ids[variant_idx][variant_id_blen]);
              if (!memequal(ref_allele, first_ref_allele, ref_blen)) {
                break;
              }
            }
            if (unlikely(variant_idx == variant_idx_end)) {
              const uint32_t variant_id_slen = variant_id_blen - 1;
              char* write_iter = strcpya_k(g_logbuf, "Error: The biallelic variants with ID '");
              write_iter = memcpya(write_iter, cur_variant_id, variant_id_slen);
              write_iter = strcpya_k(write_iter, "' at position ");
              write_iter = chrtoa(ctxp->cip, prev_chr_code, write_iter);
              *write_iter++ = ':';
              write_iter = u32toa(prev_bp, write_iter);
              write_iter = strcpya_k(write_iter, " in ");
              write_iter = strcpya(write_iter, ctxp->cur_fname);
              write_iter = strcpya_k(write_iter, " appear to be the components of a 'split' multiallelic variant; if so, it must be 'joined' (with e.g. \"bcftools norm -m\") before a correct merge can occur. Alternatively, you can keep the variants separate by assigning them distinct IDs; unless you have very long indels, adding --set-all-var-ids");
              if (strequal_k(cur_variant_id, ".", variant_id_slen)) {
                write_iter = strcpya(write_iter, " or --set-missing-var-ids");
              }
              strcpy_k(write_iter, " to your merge command is a simple way to do this.\n");
              WordWrapB(0);
              logerrputsb();
              return kPglRetInconsistentInput;
            }
          }
        }
      }
      if (merged_allele_ct > 2) {
        if (merged_allele_ct > ctxp->write_allele_ct_max) {
          if (unlikely(ctxp->write_allele_ct_max == kPglMaxAlleleCt)) {
            char* write_iter = strcpya_k(g_logbuf, "Error: Too many alleles for variant '");
            write_iter = strcpya(write_iter, cur_variant_id);
            write_iter = strcpya_k(write_iter, "' at position ");
            write_iter = chrtoa(ctxp->cip, prev_chr_code, write_iter);
            *write_iter++ = ':';
            write_iter = u32toa(prev_bp, write_iter);
            write_iter = strcpya_k(write_iter, " in ");
            write_iter = strcpya(write_iter, ctxp->cur_fname);
            strcpy_k(write_iter, ". (This " PROG_NAME_STR " build is limited to " PGL_MAX_ALLELE_CT_STR ".)\n");
            WordWrapB(0);
            logerrputsb();
            return kPglRetNotYetSupported;
          }
          ctxp->write_doomed_variant_ct += 1;
        } else if (merged_allele_ct > ctxp->write_nondoomed_max_allele_ct) {
          ctxp->write_nondoomed_max_allele_ct = merged_allele_ct;
        }
      }
      if (variant_idx_end == batch_size) {
        break;
      }
      cur_variant_id = sorted_variant_ids[variant_idx_end];
      variant_idx_start = variant_idx_end;
      allele_ct_limit = next_allele_ct;
    }

    if (ctxp->first_bp != UINT32_MAX) {
      return kPglRetSuccess;
    }
    first_varid = sorted_variant_ids[0];
  }
  ctxp->first_chr_idx = prev_chr_code;
  ctxp->first_bp = prev_bp;
  const uint32_t first_id_blen = strlen(first_varid) + 1;
  if (unlikely(pgl_malloc(first_id_blen, ctxp->first_varid_ptr))) {
    return kPglRetNomem;
  }
  memcpy(*(ctxp->first_varid_ptr), first_varid, first_id_blen);
  // Only need to natural-sort first and last positions.  We manually reset
  // this before any possibly-last call.
  ctxp->sort_vars_ascii = 1;
  return kPglRetSuccess;
}

BoolErr ScrapeLastVarid(const RescanOnePosContext* ctxp, unsigned char* arena_bottom, uint32_t batch_size, char** last_varid_ptr) {
  char* last_varid;
  if (batch_size == 1) {
    last_varid = ctxp->first_record->variant_id;
  } else {
    // See middle of RescanOnePos().  We look up the last element of the
    // sorted_variant_ids array it created.
    const uintptr_t bytes_to_round_up = (-R_CAST(uintptr_t, arena_bottom)) % sizeof(intptr_t);
    char** sorted_variant_ids = R_CAST(char**, &(arena_bottom[bytes_to_round_up]));
    last_varid = sorted_variant_ids[batch_size - 1];
  }
  const uint32_t last_id_blen = strlen(last_varid) + 1;
  free_cond(*last_varid_ptr);
  if (unlikely(pgl_malloc(last_id_blen, last_varid_ptr))) {
    return 1;
  }
  memcpy(*last_varid_ptr, last_varid, last_id_blen);
  return 0;
}

// cip->chr_file_order is filled with the final chromosome sort order.
// info_keys, pointed-to InfoVtype entries, and info_keys_htable are allocated
// at the end of bigstack.
static_assert(kCompressStreamBlock <= kDecompressChunkSize, "ScanPvarsAndMergeHeader() needs to be updated.");
PglErr ScanPvarsAndMergeHeader(const PmergeInfo* pmip, const char* missing_varid_match, MiscFlags misc_flags, LoadFilterLogFlags load_filter_log_merge_flags, uint32_t missing_varid_match_slen, char input_missing_geno_char, uint32_t max_thread_ct, SortMode sort_vars_mode, VaridTemplate* varid_templatep, VaridTemplate* varid_multi_templatep, VaridTemplate* varid_multi_nonsnp_templatep, char* outname, char* outname_end, PmergeInputFilesetLl** filesets_ptr, ChrInfo* cip, uintptr_t* fileset_ctp, uint32_t* info_has_g_keyp, const char* const** info_keys_ptr, uint32_t* info_key_ctp, uint32_t** info_keys_htablep, uint32_t* info_keys_htable_sizep, uint32_t* info_conflict_presentp) {
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
    char* chr_output_name_buf = nullptr;
    char* variant_id_buf = nullptr;
    if (varid_templatep) {
      // Don't know max_chr_blen yet.
      if (unlikely(bigstack_end_alloc_c(kMaxIdBlen, &chr_output_name_buf))) {
        goto ScanPvarsAndMergeHeader_ret_NOMEM;
      }
      varid_templatep->chr_output_name_buf = chr_output_name_buf;
      uint32_t max_base_len = varid_templatep->base_len;
      if (varid_multi_templatep) {
        varid_multi_templatep->chr_output_name_buf = chr_output_name_buf;
        if (varid_multi_templatep->base_len > max_base_len) {
          max_base_len = varid_multi_templatep->base_len;
        }
      }
      if (varid_multi_nonsnp_templatep) {
        varid_multi_nonsnp_templatep->chr_output_name_buf = chr_output_name_buf;
        if (varid_multi_nonsnp_templatep->base_len > max_base_len) {
          max_base_len = varid_multi_nonsnp_templatep->base_len;
        }
      }
      // strlen(CHROM) <= kMaxIdSlen + 3 (don't want to worry about "chr")
      // strlen(POS) <= 10
      // also need null-terminator
      if (unlikely(bigstack_end_alloc_c(kMaxIdSlen + 14 + 2 * varid_templatep->new_id_max_allele_slen, &variant_id_buf))) {
        goto ScanPvarsAndMergeHeader_ret_NOMEM;
      }
    }
    // We represent the chromosome-ordering graph as follows:
    // - chr_outedges[] is indexed by chr_idx.  When each chromosome is first
    //   seen, chr_outedges[x] is set to an empty bitarray allocated off the
    //   end of the arena.  An x -> y edge is tracked by setting bit y in
    //   chr_outedges[x].
    // - chr_inedge_cts[] is also indexed by chr_idx, and tracks how many
    //   in-edges each chromosome has.
    // TODO: this is unsuitable for files with a huge number of contigs; linked
    // lists should work better in practice
    uintptr_t* chr_present;
    uintptr_t** chr_outedges;
    uint32_t* chr_inedge_cts;
    if (unlikely(bigstack_end_calloc_w(BitCtToWordCt(kMaxContigs), &chr_present) ||
                 bigstack_end_calloc_wp(kMaxContigs, &chr_outedges) ||
                 bigstack_end_calloc_u32(kMaxContigs, &chr_inedge_cts))) {
      goto ScanPvarsAndMergeHeader_ret_NOMEM;
    }

    const uintptr_t linebuf_capacity = MINV(kMaxLongLine, bigstack_left() / 4) + kDecompressChunkSize;
    char* linebuf;
    if (unlikely(bigstack_alloc_c(linebuf_capacity, &linebuf))) {
      goto ScanPvarsAndMergeHeader_ret_NOMEM;
    }

    const MergeXheaderMode merge_xheader_mode = pmip->merge_xheader_mode;
    // Each entry is an xheader key, followed by a null terminator, followed by
    // either the null-terminated remainder of the header line or a single
    // ASCII code 1 to indicate a mismatch.  The latter includes the original
    // punctuation character ('=', ',', or '>') ending the key, but not the
    // final newline.
    // For '##' header lines where the first '=' character is followed by a
    // '<', the key is everything from the second '#' up to the first comma in
    // the '<' (or '>' if there is none) in the '<' expression; otherwise, the
    // key is everything from the second '#' up to the '='.
    // If there is no '=' character, or the header line starts with only one
    // '#', the key is the entire line (excluding the newline) and the value is
    // empty.
    char** xheader_entries = nullptr;
    uint32_t* xheader_entries_htable = nullptr;
    uint32_t xheader_entry_ct = 0;
    // xheader_entry_capacity == xheader_entry_htable_size / 2
    uint32_t xheader_entry_htable_size = 0;
    if (merge_xheader_mode != kMergeXheaderModeErase) {
      xheader_entry_htable_size = 512;
      if (unlikely(bigstack_alloc_cp(xheader_entry_htable_size / 2, &xheader_entries) ||
                   bigstack_alloc_u32(xheader_entry_htable_size, &xheader_entries_htable))) {
        goto ScanPvarsAndMergeHeader_ret_NOMEM;
      }
      SetAllU32Arr(xheader_entry_htable_size, xheader_entries_htable);
    }
    unsigned char* arena_bottom = g_bigstack_base;
    unsigned char* arena_top = g_bigstack_end;

    const uintptr_t fileset_ct = *fileset_ctp;
    const uint32_t decompress_thread_ct = MAXV(max_thread_ct - 1, 1);
    RescanOnePosContext ctx;
    ctx.write_allele_ct_max = pmip->max_allele_ct? pmip->max_allele_ct : kPglMaxAlleleCt;
    ctx.sort_vars_ascii = (sort_vars_mode == kSortAscii);
    ctx.multiallelics_already_joined = (pmip->flags / kfPmergeMultiallelicsAlreadyJoined) & 1;
    ctx.input_missing_geno_char = input_missing_geno_char;
    ctx.cip = cip;
    // For each .pvar, need to determine:
    // - variant_ct
    // - write_variant_ct (same as variant_ct unless chromosome filter or
    //   negative POS)
    // - write_nondoomed_variant_ct (--merge-max-alleles applied)
    // - max_line_blen
    // - max_single_pos_ct
    // - max_single_pos_blen
    // - first and last (chr_idx, pos, varid)
    // - nm_{qual,filter,info}_present, nz_cm_present
    // - write_max_allele_ct
    // possible todo: track seen realpaths, skip duplicates
    // possible todo: after merging headers, split workspace into equal-size
    //                pieces and optimistically scan .pvar bodies in parallel.
    //                Any file that cannot be processed with the smaller
    //                per-thread workspace is retried afterward with the
    //                workspace merged again.  cip and the chromosome-ordering
    //                graph would be shared and protected by a mutex; the
    //                string representation of the last-seen chromosome code
    //                would be cached to reduce the need for cip
    //                lookups/updates, and the chromosome-ordering graph would
    //                be heap-allocated using a sparse representation (unless
    //                we don't care about supporting higher compile-time
    //                kMaxContigs values, then we can just allocate a ~500MB
    //                matrix).
    // Lots of overlap with LoadPvar().
    PmergeInputFilesetLl** filesets_iterp = filesets_ptr;
    // Chromosome set must be either defined on the command line, or there must
    // be equivalent chrSet header lines in *all* .pvar files.
    ChrsetSource orig_chrset_source = cip->chrset_source;
    uintptr_t null_fileset_ct = 0;
    uint32_t info_pr_exists = 0;
    uint32_t info_pr_nonflag_exists = 0;
    uint32_t max_xheader_line_blen = 0;
    uint32_t at_least_one_info_exists = 0;
    uint32_t max_allele_overflow_slen = 0;
    uint32_t neg_bp_seen = 0;
    uintptr_t info_conflict_ct = 0;
    for (uintptr_t fileset_idx = 0; fileset_idx != fileset_ct; ++fileset_idx) {
      PmergeInputFilesetLl* cur_fileset = *filesets_iterp;
      cur_fname = cur_fileset->pvar_fname;
      reterr = TextStreamOpenEx(cur_fname, kMaxLongLine, linebuf_capacity, decompress_thread_ct, nullptr, linebuf, &txs);
      if (unlikely(reterr)) {
        goto ScanPvarsAndMergeHeader_ret_TSTREAM_FAIL;
      }
      uint32_t max_line_blen = 0;
      char* line_start = TextLineEnd(&txs);
      uint32_t info_pr_exists_here = 0;
      uint32_t chrset_seen_in_this_file = 0;
      for (line_idx = 1; ; ++line_idx) {
        if (unlikely(!TextGetUnsafe2(&txs, &line_start))) {
          if (TextStreamErrcode2(&txs, &reterr)) {
            goto ScanPvarsAndMergeHeader_ret_TSTREAM_FAIL;
          }
          logerrprintf("Error: No variants in %s.\n", cur_fname);
          goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT;
        }
        if ((line_start[0] != '#') || tokequal_k(line_start, "#CHROM")) {
          break;
        }
        char* key_start = &(line_start[1]);
        // Contents of ##chrSet and ##INFO=<ID=PR,...> header lines matter even
        // when we're not propagating the header.
        if (StrStartsWithUnsafe(key_start, "#chrSet=<")) {
          if (unlikely(chrset_seen_in_this_file)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Multiple ##chrSet header lines in %s.\n", cur_fname);
            goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
          }
          reterr = ReadChrsetHeaderLine(&(key_start[strlen("#chrSet=<")]), cur_fname, misc_flags, line_idx, cip);
          if (unlikely(reterr)) {
            goto ScanPvarsAndMergeHeader_ret_1;
          }
          if (cip->chrset_source != kChrsetSourceCmdline) {
            if (unlikely(fileset_idx && (cip->chrset_source == kChrsetSourceDefault))) {
              goto ScanPvarsAndMergeHeader_ret_INCONSISTENTLY_PRESENT_CHRSET;
            }
            orig_chrset_source = kChrsetSourceAnotherFile;
            cip->chrset_source = kChrsetSourceCmdline;
          }
          chrset_seen_in_this_file = 1;
        }
        if (StrStartsWithUnsafe(key_start, "#INFO=<")) {
          char* info_read_iter = &(key_start[strlen("#INFO=<")]);
          char* idval;
          uint32_t id_slen;
          if (unlikely(HkvlineId(&info_read_iter, &idval, &id_slen))) {
            snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s is malformed.\n", line_idx, cur_fname);
            goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
          }
          if (strequal_k(idval, "PR", id_slen)) {
            char* typestr;
            uint32_t type_slen;
            if (unlikely(HkvlineFind(info_read_iter, "Type", &typestr, &type_slen))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s is malformed.\n", line_idx, cur_fname);
              goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
            }
            if (strequal_k(typestr, "Flag", type_slen)) {
              if (unlikely(info_pr_nonflag_exists)) {
                // Other conflicting keys can be skipped under
                // "--merge-info-mode nm-match", but we don't want to remove
                // the INFO/PR key if it's present.
                logerrputs("Error: Inconsistent INFO/PR header lines across --pmerge[-list] files.\n");
                goto ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT;
              }
              info_pr_exists = 1;
              info_pr_exists_here = 1;
              line_start = AdvPastDelim(&(typestr[type_slen]), '\n');
              continue;
            }
            if (unlikely(info_pr_exists)) {
              logerrputs("Error: Inconsistent INFO/PR header lines across --pmerge[-list] files.\n");
              goto ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT;
            }
            info_pr_nonflag_exists = 1;
          }
          if (!xheader_entries_htable) {
            line_start = AdvPastDelim(info_read_iter, '\n');
            continue;
          }
        } else {
          // if the "pvar file" was actually a VCF, suppress the same lines
          // we'd suppress when importing with --vcf.
          if ((!xheader_entries_htable) ||
              StrStartsWithUnsafe(key_start, "#fileformat=") ||
              StrStartsWithUnsafe(key_start, "#fileDate=") ||
              StrStartsWithUnsafe(key_start, "#source=") ||
              StrStartsWithUnsafe(key_start, "#FORMAT=")) {
            line_start = AdvPastDelim(line_start, '\n');
            continue;
          }
        }
        char* value_start;
        char* value_end;
        if (line_start[1] != '#') {
          value_start = AdvToDelim(key_start, '\n');
          value_end = value_start;
        } else {
          value_start = strchrnul_n(&(key_start[1]), '=');
          if (*value_start == '\n') {
            value_end = value_start;
          } else {
            if (value_start[1] == '<') {
              char* hline_kv_start = &(value_start[2]);
              if (unlikely(HkvlineForceIdFirst(arena_top - arena_bottom, &(value_start[2]), arena_bottom))) {
                snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s is malformed.\n", line_idx, cur_fname);
                goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
              }
              value_start = strchrnul2_n(hline_kv_start, ',', '>');
              if (unlikely(*value_start == '\n')) {
                snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s is malformed (value starts with '<', but there is no closing '>').\n", line_idx, cur_fname);
                goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
              }
            }
            value_end = AdvToDelim(&(value_start[1]), '\n');
          }
        }
        const uint32_t key_slen = value_start - key_start;
        if (key_slen) {
          const uint32_t value_slen = value_end - value_start;
          // might want to move this logic into plink2_common
          for (uint32_t hashval = Hashceil(key_start, key_slen, xheader_entry_htable_size); ; ) {
            const uint32_t cur_htable_entry = xheader_entries_htable[hashval];
            if (cur_htable_entry == UINT32_MAX) {
              if (unlikely(PtrWSubCk(arena_bottom, key_slen + value_slen + 2, &arena_top))) {
                goto ScanPvarsAndMergeHeader_ret_NOMEM;
              }
              char* entry_start = R_CAST(char*, arena_top);
              char* entry_iter = memcpyax(entry_start, key_start, key_slen, '\0');
              memcpyx(entry_iter, value_start, value_slen, '\0');
              if (xheader_entry_ct * 2 < xheader_entry_htable_size) {
                xheader_entries[xheader_entry_ct] = entry_start;
                xheader_entries_htable[hashval] = xheader_entry_ct;
                ++xheader_entry_ct;
              } else {
                // resize
                uint32_t next_xheader_table_size;
                if (xheader_entry_htable_size < 0x80000000U) {
                  next_xheader_table_size = xheader_entry_htable_size * 2;
                } else if (likely(xheader_entry_htable_size == 0x80000000U)) {
                  next_xheader_table_size = 0xfffffffaU;
                } else {
                  logerrputs("Error: --pmerge[-list] is limited to 2^31 - 3 header lines.\n");
                  reterr = kPglRetNotYetSupported;
                  goto ScanPvarsAndMergeHeader_ret_1;
                }
                if (unlikely(S_CAST(uintptr_t, arena_top - arena_bottom) < (next_xheader_table_size - xheader_entry_htable_size) * ((sizeof(intptr_t) / 2) + sizeof(int32_t)))) {
                  goto ScanPvarsAndMergeHeader_ret_NOMEM;
                }
                xheader_entry_htable_size = next_xheader_table_size;
                xheader_entries_htable = R_CAST(uint32_t*, &(xheader_entries[xheader_entry_htable_size / 2]));
                arena_bottom = R_CAST(unsigned char*, &(xheader_entries_htable[xheader_entry_htable_size]));
                SetAllU32Arr(xheader_entry_htable_size, xheader_entries_htable);
                xheader_entries[xheader_entry_ct] = entry_start;
                ++xheader_entry_ct;
                for (uint32_t entry_idx = 0; entry_idx != xheader_entry_ct; ++entry_idx) {
                  const char* cur_entry = xheader_entries[entry_idx];
                  HtableAddNondup(cur_entry, strlen(cur_entry), xheader_entry_htable_size, entry_idx, xheader_entries_htable);
                }
              }
              break;
            }
            if (strequal_unsafe(xheader_entries[cur_htable_entry], key_start, key_slen)) {
              char* entry_value_start = &(xheader_entries[cur_htable_entry][key_slen + 1]);
              if (merge_xheader_mode == kMergeXheaderModeMatch) {
                // entry_value_start[0] == 1 marks a conflict.
                if (entry_value_start[0] != 1) {
                  if (!strequal_unsafe(entry_value_start, value_start, value_slen)) {
                    info_conflict_ct += StrStartsWithUnsafe(key_start, "#INFO=<");
                    entry_value_start[0] = 1;
                  }
                }
              } else if (StrStartsWithUnsafe(key_start, "#INFO=<")) {
                // Check for possible Number= incompatibility.
                char* existing_read_iter = &(entry_value_start[1]);
                char* numstr;
                uint32_t num_slen;
                int32_t existing_info_vtype_num;
                if (unlikely(HkvlineFind(existing_read_iter, "Number", &numstr, &num_slen) ||
                             FillInfoVtypeNum(numstr, num_slen, info_has_g_keyp, &existing_info_vtype_num))) {
                  char* write_iter = strcpya_k(g_logbuf, "Error: --pmerge");
                  if (pmip->list_fname) {
                    write_iter = strcpya_k(write_iter, "-list");
                  }
                  write_iter = strcpya_k(write_iter, ": First INFO/");
                  write_iter = memcpya(write_iter, &(key_start[strlen("#INFO=<ID=")]), key_slen - strlen("#INFO=<ID="));
                  strcpy_k(write_iter, " header line is malformed or unrecognized.\n");
                  goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
                }
                char* cur_read_iter = &(value_start[1]);
                int32_t suppressed_info_vtype_num;
                if (unlikely(HkvlineFind(cur_read_iter, "Number", &numstr, &num_slen) ||
                             FillInfoVtypeNum(numstr, num_slen, info_has_g_keyp, &suppressed_info_vtype_num))) {
                  snprintf(g_logbuf, kLogbufSize, "Error: Header line %" PRIuPTR " of %s is malformed.\n", line_idx, cur_fname);
                  goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
                }
                if ((existing_info_vtype_num != suppressed_info_vtype_num) &&
                    (existing_info_vtype_num != kInfoVtypeUnknown) &&
                    (suppressed_info_vtype_num != kInfoVtypeUnknown)) {
                  char* write_iter = strcpya_k(g_logbuf, "Warning: INFO/");
                  write_iter = memcpya(write_iter, &(key_start[strlen("#INFO=<ID=")]), key_slen - strlen("#INFO=<ID="));
                  write_iter = strcpya_k(write_iter, " header entry on line ");
                  write_iter = wtoa(line_idx, write_iter);
                  write_iter = strcpya_k(write_iter, " of ");
                  write_iter = strcpya(write_iter, cur_fname);
                  strcpy_k(write_iter, " has a different Number= value than an earlier file; this key may be merged incorrectly, if the merge completes at all.\n");
                  WordWrapB(0);
                  logerrputsb();
                }
                // possible todo: also warn on most type mismatches
              }
              break;
            }
            if (++hashval == xheader_entry_htable_size) {
              hashval = 0;
            }
          }
        }
        char* line_end = &(value_end[1]);
        uint32_t line_blen = line_end - line_start;
        if (max_line_blen < line_blen) {
          max_line_blen = line_blen;
        }
        line_start = line_end;
      }
      if (unlikely((orig_chrset_source == kChrsetSourceAnotherFile) && (!chrset_seen_in_this_file))) {
        goto ScanPvarsAndMergeHeader_ret_INCONSISTENTLY_PRESENT_CHRSET;
      }
      if (!fileset_idx) {
        FinalizeChrset(load_filter_log_merge_flags, cip);
      }
      if (max_xheader_line_blen < max_line_blen) {
        max_xheader_line_blen = max_line_blen;
      }
      // [-1] = #CHROM (must be first column)
      // [0] = POS
      // [1] = ID
      // [2] = REF
      // [3] = ALT
      // [4] = QUAL
      // [5] = FILTER
      // [6] = INFO
      // [7] = CM (usually absent)
      uint32_t col_skips[8];
      uint32_t col_types[8];
      uint32_t no_multiallelic_allowed = 0;
      uint32_t check_qual = 0;
      uint32_t check_filter = 0;
      uint32_t info_col_present = 0;
      uint32_t check_info = 0;
      uint32_t check_cm = 0;
      uint32_t relevant_postchr_col_ct;
      if (line_start[0] == '#') {
        char* token_end = &(line_start[6]);
        uint32_t found_header_bitset = 0;
        relevant_postchr_col_ct = 0;
        char* token_start;
        for (uint32_t col_idx = 1; ; ++col_idx) {
          token_start = FirstNonTspace(token_end);
          if (IsEolnKns(*token_start)) {
            break;
          }
          token_end = CurTokenEnd(token_start);
          const uint32_t token_slen = token_end - token_start;
          uint32_t cur_col_type;
          if (token_slen <= 3) {
            if (token_slen == 3) {
              if (memequal_sk(token_start, "POS")) {
                cur_col_type = 0;
              } else if (memequal_sk(token_start, "REF")) {
                cur_col_type = 2;
              } else if (memequal_sk(token_start, "ALT")) {
                cur_col_type = 3;
              } else {
                continue;
              }
            } else if (token_slen == 2) {
              if (memequal_sk(token_start, "ID")) {
                cur_col_type = 1;
              } else if (memequal_sk(token_start, "CM")) {
                cur_col_type = 7;
                check_cm = 1;
              } else {
                continue;
              }
            } else {
              continue;
            }
          } else if (strequal_k(token_start, "QUAL", token_slen)) {
            if (pmip->merge_qual_mode == kMergeQualModeErase) {
              continue;
            }
            cur_col_type = 4;
            check_qual = 1;
          } else if (strequal_k(token_start, "INFO", token_slen)) {
            info_col_present = 1;
            if (pmip->merge_info_mode == kMergeInfoCmModeErase) {
              continue;
            }
            cur_col_type = 6;
            check_info = 1;
          } else if (token_slen == 6) {
            if (memequal_sk(token_start, "FILTER")) {
              if (pmip->merge_filter_mode == kMergeFilterModeErase) {
                continue;
              }
              cur_col_type = 5;
              check_filter = 1;
            } else if (memequal_sk(token_start, "FORMAT")) {
              break;
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
            write_iter = memcpya(write_iter, token_start, token_slen);
            write_iter = strcpya_k(write_iter, "' on line ");
            write_iter = wtoa(line_idx, write_iter);
            write_iter = strcpya_k(write_iter, " of ");
            write_iter = strcpya(write_iter, cur_fname);
            memcpy_k(write_iter, ".\n\0", 4);
            goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
          }
          found_header_bitset |= cur_col_type_shifted;
          col_skips[relevant_postchr_col_ct] = col_idx;
          col_types[relevant_postchr_col_ct++] = cur_col_type;
        }
        if (unlikely((found_header_bitset & 0x0f) != 0x0f)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Missing column header(s) on line %" PRIuPTR " of %s. (POS, ID, REF, and ALT are required.)\n", line_idx, cur_fname);
          goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
        }
        for (uint32_t rpc_col_idx = relevant_postchr_col_ct - 1; rpc_col_idx; --rpc_col_idx) {
          col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
        }
        // skip this line in main loop
        char* line_end = AdvPastDelim(token_start, '\n');
        const uint32_t line_blen = line_end - line_start;
        if (max_line_blen < line_blen) {
          max_line_blen = line_blen;
        }
        line_start = line_end;
        ++line_idx;
      } else {
        no_multiallelic_allowed = 1;
        col_skips[0] = 1;
        col_skips[1] = 1;
        col_skips[2] = 1;
        col_skips[3] = 1;
        col_types[0] = 1;
        const char* fifth_col_start = NextTokenMult(line_start, 4);
        if (unlikely(!fifth_col_start)) {
          goto ScanPvarsAndMergeHeader_ret_MISSING_TOKENS;
        }
        const char* sixth_col_start = NextToken(fifth_col_start);
        if (!sixth_col_start) {
          relevant_postchr_col_ct = 4;
          col_types[1] = 0;
          col_types[2] = 3;
          col_types[3] = 2;
        } else {
          relevant_postchr_col_ct = 5;
          col_skips[4] = 1;
          col_types[1] = 7;
          col_types[2] = 0;
          col_types[3] = 3;
          col_types[4] = 2;
          check_cm = 1;
        }
      }
      // In order to perform 'concatenation' with only one more pass through
      // each .pvar, without making that yield a different result than
      // general-purpose merge, we want to track (variant ID, alleles) for
      // each variant in the current group of same-position variant(s).  (It is
      // not necessary to distinguish REF/ALT here.)  This is necessary to
      // compute write_nondoomed_variant_ct and write_nondoomed_max_allele_ct
      // accurately, both of which must be known before the .pgen writer can be
      // constructed.
      //
      // We store this as a sequence of records growing up from
      // arena_bottom_mark, structured as follows:
      //   4 byte uint32_t, storing record length in bytes
      //   AlleleCode storing extra_alt_ct
      //   null-terminated variant ID
      //   null-terminated REF
      //   null-terminated ALT, internally still comma-separated
      unsigned char* arena_bottom_mark = arena_bottom;
      const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
      const uint32_t filter_count_needed = (pmip->merge_filter_mode == kMergeFilterModeNonpassUnion) || (pmip->merge_filter_mode == kMergeFilterModeNmMatch);
      const uintptr_t line_idx_body_start = line_idx;
      uint32_t nonwrite_variant_ct = 0;
      ctx.cur_fname = cur_fname;
      ctx.first_record = R_CAST(RescanOnePosRecord*, arena_bottom_mark);
      ctx.write_doomed_variant_ct = 0;
      ctx.write_nondoomed_max_allele_ct = 2;
      // ctx.first_chr_idx = 0;
      ctx.first_bp = UINT32_MAX;
      ctx.first_varid_ptr = &(cur_fileset->first_varid);
      uint32_t cur_single_pos_ct = 0;
      uint32_t max_single_pos_ct = 1;
      uintptr_t cur_single_pos_blen = 0;
      uintptr_t max_single_pos_blen = 0;
      uint32_t read_max_allele_ct = 2;
      uint32_t read_max_nonpass_filter_ct = 0;
      uint32_t prev_chr_code = UINT32_MAX;
      int32_t prev_bp = 0;
      cur_fileset->nm_qual_exists = 0;
      cur_fileset->nm_filter_exists = 0;
      cur_fileset->nm_info_exists = 0;
      cur_fileset->nz_cm_present = 0;

      info_pr_exists_here = info_pr_exists_here && info_col_present;
      if (unlikely((cur_fileset->nonref_flags_storage == 0) && info_pr_exists_here)) {
        // Provisional-vs.-not REF status is directly relevant during merge.
        // It's of limited importance for most other operations, so
        // Plink2Core() only prints a warning.
        snprintf(g_logbuf, kLogbufSize, "Error: %s indicates that provisional-REF information (required by --pmerge%s) is stored in the companion .pvar, but that .pvar lacks an INFO/PR header line and/or an INFO column.\n", cur_fileset->pgen_fname, pmip->list_fname? "-list" : "");
        goto ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT_WW;
      }
      cur_fileset->pvar_info_pr_exists = info_pr_exists_here;
      for (; TextGetUnsafe2(&txs, &line_start); ++line_idx) {
        if (unlikely(line_start[0] == '#')) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s starts with a '#'. (This is only permitted before the first nonheader line, and if a #CHROM header line is present it must denote the end of the header block.)\n", line_idx, cur_fname);
          goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
        }
        char* first_token_end = CurTokenEnd(line_start);
        if (unlikely(*first_token_end == '\n')) {
          goto ScanPvarsAndMergeHeader_ret_MISSING_TOKENS;
        }
        uint32_t cur_chr_code;
        reterr = GetOrAddChrCodeDestructive(cur_fname, line_idx, prohibit_extra_chr, line_start, first_token_end, cip, &cur_chr_code);
        if (unlikely(reterr)) {
          goto ScanPvarsAndMergeHeader_ret_1;
        }
        if (cur_chr_code != prev_chr_code) {
          SetBit(cur_chr_code, chr_present);
          if (chr_output_name_buf) {
            char* chr_name_end = chrtoa(cip, cur_chr_code, chr_output_name_buf);
            const uint32_t chr_slen = chr_name_end - chr_output_name_buf;
            const int32_t chr_slen_delta = chr_slen - varid_templatep->chr_slen;
            varid_templatep->chr_slen = chr_slen;
            varid_templatep->base_len += chr_slen_delta;
            if (varid_multi_templatep) {
              varid_multi_templatep->chr_slen = chr_slen;
              varid_multi_templatep->base_len += chr_slen_delta;
            }
            if (varid_multi_nonsnp_templatep) {
              varid_multi_nonsnp_templatep->chr_slen = chr_slen;
              varid_multi_nonsnp_templatep->base_len += chr_slen_delta;
            }
          }
          if (prev_chr_code != UINT32_MAX) {
            if (cur_single_pos_ct) {
              // Rescan now, before clobbering prev_chr_code/prev_bp.
              if (max_single_pos_ct < cur_single_pos_ct) {
                max_single_pos_ct = cur_single_pos_ct;
              }
              if (max_single_pos_blen < cur_single_pos_blen) {
                max_single_pos_blen = cur_single_pos_blen;
              }
              // This could be the last included variant in the entire file
              // (e.g. all remaining POS values could be -1), in which case we
              // need to save off last_varid now.
              ctx.sort_vars_ascii = (sort_vars_mode == kSortAscii);
              reterr = RescanOnePos(arena_top, cur_single_pos_ct, prev_chr_code, prev_bp, arena_bottom, &ctx, &nonwrite_variant_ct);
              if (unlikely(reterr)) {
                goto ScanPvarsAndMergeHeader_ret_1;
              }
              ctx.sort_vars_ascii = 1;
              if (unlikely(ScrapeLastVarid(&ctx, arena_bottom, cur_single_pos_ct, &cur_fileset->last_varid))) {
                goto ScanPvarsAndMergeHeader_ret_NOMEM;
              }
              cur_fileset->last_chr_idx = prev_chr_code;
              cur_fileset->last_pos = prev_bp;
            }
            arena_bottom = arena_bottom_mark;
            cur_single_pos_ct = 0;
            cur_single_pos_blen = 0;
            // Add prev_chr_code -> cur_chr_code graph edge.
            if (!chr_outedges[prev_chr_code]) {
              ArenaEndSet(arena_top, &arena_top);
              if (unlikely(arena_end_alloc_w(arena_bottom, BitCtToWordCt(kMaxContigs), &arena_top, &(chr_outedges[prev_chr_code])))) {
                goto ScanPvarsAndMergeHeader_ret_NOMEM;
              }
              ZeroWArr(BitCtToWordCt(kMaxContigs), chr_outedges[prev_chr_code]);
            }
            if (!IsSet(chr_outedges[prev_chr_code], cur_chr_code)) {
              SetBit(cur_chr_code, chr_outedges[prev_chr_code]);
              chr_inedge_cts[cur_chr_code] += 1;
            }
          }
          prev_chr_code = cur_chr_code;
          // no explicit split-chr check needed here, we'll error out anyway
          // during topological sort
          prev_bp = -1;
        }

        *first_token_end = '\t';
        char* token_ptrs[8];
        uint32_t token_slens[8];
        char* line_iter = TokenLex(first_token_end, col_types, col_skips, relevant_postchr_col_ct, token_ptrs, token_slens);
        if (unlikely(!line_iter)) {
          goto ScanPvarsAndMergeHeader_ret_MISSING_TOKENS;
        }
        const char* alt_start = token_ptrs[3];
        const uint32_t alt_slen = token_slens[3];
        const uint32_t extra_alt_ct = CountByte(alt_start, ',', alt_slen);
        if (unlikely(extra_alt_ct >= kPglMaxAltAlleleCt)) {
          logerrprintfww("Error: Too many ALT alleles on line %" PRIuPTR " of %s. (This " PROG_NAME_STR " build is limited to " PGL_MAX_ALT_ALLELE_CT_STR ".)\n", line_idx, cur_fname);
          reterr = kPglRetNotYetSupported;
          goto ScanPvarsAndMergeHeader_ret_1;
        }

        char* line_end = AdvPastDelim(line_iter, '\n');
        uint32_t line_blen = line_end - line_start;
        if (max_line_blen < line_blen) {
          max_line_blen = line_blen;
        }
        line_start = line_end;

        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          ++nonwrite_variant_ct;
          continue;
        }
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(token_ptrs[0], &cur_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid POS on line %" PRIuPTR " of %s.\n", line_idx, cur_fname);
          goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
        }
        const uint32_t ref_slen = token_slens[2];
        char* variant_id = token_ptrs[1];
        uint32_t id_slen = token_slens[1];
        if (varid_templatep && ((!missing_varid_match_slen) || ((id_slen == missing_varid_match_slen) && memequal(variant_id, missing_varid_match, missing_varid_match_slen)))) {
          const VaridTemplate* cur_varid_templatep = varid_templatep;
          if (extra_alt_ct && (varid_multi_templatep || varid_multi_nonsnp_templatep)) {
            if (varid_multi_templatep) {
              cur_varid_templatep = varid_multi_templatep;
            }
            if (varid_multi_nonsnp_templatep) {
              if ((ref_slen > 1) || (token_slens[3] != 2 * extra_alt_ct + 1)) {
                cur_varid_templatep = varid_multi_nonsnp_templatep;
              }
            }
          }
          variant_id = variant_id_buf;
          char* variant_id_end = VaridTemplateWrite(cur_varid_templatep, token_ptrs[2], alt_start, cur_bp, ref_slen, extra_alt_ct, alt_slen, &max_allele_overflow_slen, variant_id);
          // error message fix (11 Oct 2024)
          if (max_allele_overflow_slen && (!(misc_flags & (kfMiscNewVarIdOverflowMissing | kfMiscNewVarIdOverflowTruncate)))) {
            ctx.multiallelics_already_joined = 1;
          }
          id_slen = variant_id_end - variant_id;
          // bugfix (11 Oct 2024): need to update line_blen value used in
          // cur_single_pos_blen accounting.
          line_blen += id_slen - token_slens[1];
        }
        if (unlikely(id_slen > kMaxIdSlen)) {
          logerrputs("Error: Variant IDs are limited to " MAX_ID_SLEN_STR " characters.\n");
          goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT;
        }
        if (cur_bp <= prev_bp) {
          if (cur_bp < 0) {
            neg_bp_seen = 1;
            ++nonwrite_variant_ct;
            continue;
          }
          if (unlikely(cur_bp < prev_bp)) {
            snprintf(g_logbuf, kLogbufSize, "Error: %s is not position-sorted. Retry --pmerge[-list] after using --make-pgen/--make-bed + --sort-vars to sort your data.\n", cur_fname);
            goto ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT_WW;
          }
          // same position as previous included variant
          ++cur_single_pos_ct;
          cur_single_pos_blen += line_blen;
        } else {
          if (max_single_pos_ct < cur_single_pos_ct) {
            max_single_pos_ct = cur_single_pos_ct;
          }
          if (max_single_pos_blen < cur_single_pos_blen) {
            max_single_pos_blen = cur_single_pos_blen;
          }
          reterr = RescanOnePos(arena_top, cur_single_pos_ct, prev_chr_code, prev_bp, arena_bottom, &ctx, &nonwrite_variant_ct);
          if (unlikely(reterr)) {
            goto ScanPvarsAndMergeHeader_ret_1;
          }
          arena_bottom = arena_bottom_mark;
          cur_single_pos_ct = 1;
          cur_single_pos_blen = line_blen;
          prev_bp = cur_bp;
        }
        variant_id[id_slen] = '\0';
        const uint32_t id_blen = id_slen + 1;
        const uint32_t rec_blen = sizeof(int32_t) + sizeof(AlleleCode) + id_blen + ref_slen + alt_slen + 2;
        if (S_CAST(uintptr_t, arena_top - arena_bottom) < rec_blen) {
          goto ScanPvarsAndMergeHeader_ret_NOMEM;
        }
        RescanOnePosRecord* cur_record = R_CAST(RescanOnePosRecord*, arena_bottom);
        arena_bottom = &(arena_bottom[rec_blen]);
        cur_record->rec_blen = rec_blen;
        const uint32_t cur_allele_ct = extra_alt_ct + 2;
        cur_record->allele_ct = cur_allele_ct;
        if (read_max_allele_ct < cur_allele_ct) {
          read_max_allele_ct = cur_allele_ct;
        }
        char* write_iter = memcpya(cur_record->variant_id, variant_id, id_blen);
        write_iter = memcpyax(write_iter, token_ptrs[2], ref_slen, '\0');
        memcpyx(write_iter, alt_start, alt_slen, '\0');

        if (check_qual) {
          const char* qual_token = token_ptrs[4];
          if ((qual_token[0] != '.') || (qual_token[1] > ' ')) {
            cur_fileset->nm_qual_exists = 1;
            // possible todo: update col_types and col_skips to remove this
            // column.
            check_qual = 0;
          }
        }
        if (check_filter) {
          const char* filter_token = token_ptrs[5];
          const uint32_t filter_slen = token_slens[5];
          if ((filter_slen > 1) || (filter_token[0] != '.')) {
            cur_fileset->nm_filter_exists = 1;
            if (filter_count_needed) {
              if (!strequal_k(filter_token, "PASS", filter_slen)) {
                const uint32_t cur_filter_ct_m1 = CountByte(filter_token, ';', filter_slen);
                if (cur_filter_ct_m1 >= read_max_nonpass_filter_ct) {
                  read_max_nonpass_filter_ct = cur_filter_ct_m1 + 1;
                }
              }
            } else {
              check_filter = 0;
            }
          }
        }
        if (check_info) {
          const char* info_token = token_ptrs[6];
          const uint32_t info_slen = token_slens[6];
          if ((info_slen > 1) || (info_token[0] != '.')) {
            cur_fileset->nm_info_exists = 1;
            at_least_one_info_exists = 1;
            check_info = 0;
          }
        }
        if (check_cm) {
          const char* cm_token = token_ptrs[7];
          if ((cm_token[0] != '0') || (cm_token[1] > ' ')) {
            double cur_cm;
            if (unlikely(!ScantokDouble(cm_token, &cur_cm))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, cur_fname);
              goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
            }
            if (cur_cm != 0.0) {
              cur_fileset->nz_cm_present = 1;
              check_cm = 0;
            }
          }
        }
      }
      if (unlikely(CleanupTextStream2(cur_fname, &txs, &reterr))) {
        goto ScanPvarsAndMergeHeader_ret_1;
      }
      const uintptr_t read_variant_ct = line_idx - line_idx_body_start;
      if (unlikely(read_variant_ct > kPglMaxVariantCt)) {
        logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
        goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT;
      }
      if (cur_single_pos_ct || (read_variant_ct != nonwrite_variant_ct)) {
        cur_fileset->read_variant_ct = read_variant_ct;
        cur_fileset->max_pvar_line_blen = max_line_blen;
        if (max_single_pos_ct < cur_single_pos_ct) {
          max_single_pos_ct = cur_single_pos_ct;
        }
        cur_fileset->max_single_pos_ct = max_single_pos_ct;
        if (cur_single_pos_ct) {
          ctx.sort_vars_ascii = (sort_vars_mode == kSortAscii);
          reterr = RescanOnePos(arena_top, cur_single_pos_ct, prev_chr_code, prev_bp, arena_bottom, &ctx, &nonwrite_variant_ct);
          if (unlikely(reterr)) {
            goto ScanPvarsAndMergeHeader_ret_1;
          }
          if (unlikely(ScrapeLastVarid(&ctx, arena_bottom, cur_single_pos_ct, &cur_fileset->last_varid))) {
            goto ScanPvarsAndMergeHeader_ret_NOMEM;
          }
          cur_fileset->last_chr_idx = prev_chr_code;
          cur_fileset->last_pos = prev_bp;
        }
        cur_fileset->first_chr_idx = ctx.first_chr_idx;
        cur_fileset->first_pos = ctx.first_bp;
        if (max_single_pos_blen < cur_single_pos_blen) {
          max_single_pos_blen = cur_single_pos_blen;
        }
        cur_fileset->max_single_pos_blen = max_single_pos_blen;
        const uint32_t write_variant_ct = read_variant_ct - nonwrite_variant_ct;
        cur_fileset->write_variant_ct = write_variant_ct;
        cur_fileset->write_nondoomed_variant_ct = write_variant_ct - ctx.write_doomed_variant_ct;
        if (unlikely(no_multiallelic_allowed && (read_max_allele_ct > 2))) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s contains multiallelic variant(s), despite having no #CHROM header line. Add that header line to make it obvious that this isn't a valid .bim.\n", cur_fname);
          goto ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW;
        }
        cur_fileset->read_max_allele_ct = read_max_allele_ct;
        cur_fileset->write_nondoomed_max_allele_ct = ctx.write_nondoomed_max_allele_ct;
        cur_fileset->read_max_nonpass_filter_ct = read_max_nonpass_filter_ct;
        arena_bottom = arena_bottom_mark;
        filesets_iterp = &((*filesets_iterp)->next);
      } else {
        PmergeInputFilesetLl** next_filesets_iterp = &((*filesets_iterp)->next);
        *filesets_iterp = cur_fileset->next;
        filesets_iterp = next_filesets_iterp;
        ++null_fileset_ct;
      }
    }
    if (unlikely(null_fileset_ct && ((null_fileset_ct == fileset_ct) || (pmip->flags & kfPmergeVariantInnerJoin)))) {
      char* write_iter = strcpya_k(g_logbuf, "Error: --pmerge[-list]: All variants excluded by ");
      if (load_filter_log_merge_flags) {
        AppendLoadFilterFlagnames(load_filter_log_merge_flags, &write_iter);
      }
      if (neg_bp_seen) {
        if (load_filter_log_merge_flags) {
          write_iter = strcpya_k(write_iter, " and/or ");
        }
        write_iter = strcpya_k(write_iter, "negative bp coordinates");
      }
      strcpy_k(write_iter, ".\n");
      goto ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT_WW;
    }
    cip->chrset_source = orig_chrset_source;

    BigstackEndSet(arena_top);
    if (unlikely(BigstackBaseSetChecked(arena_bottom))) {
      goto ScanPvarsAndMergeHeader_ret_NOMEM;
    }
    // defensive
    arena_bottom = nullptr;
    arena_top = nullptr;

    if (xheader_entries && (!info_pr_exists)) {
      // If INFO/PR is not present in any input .pvar/.bim, but nonref_flags
      // is stored in at least one input .pgen, and we're already writing the
      // INFO column for other reasons, we *do* want to write INFO/PR.
      PmergeInputFilesetLl* filesets_iter = *filesets_ptr;
      do {
        if (filesets_iter->nonref_flags_storage == 3) {
          if (unlikely(info_pr_nonflag_exists)) {
            logerrputs("Error: Conflicting INFO/PR definitions.  Either fix all REF alleles so that the\n\"provisional reference\" flag is no longer needed, or remove/rename the other\nuse of the INFO/PR key.\n");
            goto ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT;
          }
          info_pr_exists = 1;
          break;
        }
        filesets_iter = filesets_iter->next;
      } while (filesets_iter != nullptr);
    }
    // Note that, under "--merge-info-mode erase", we may need to scrape
    // INFO/PR from .pvar files to be saved later to the .pgen.
    const uint32_t write_info_pr_to_pvar = info_pr_exists && (pmip->merge_info_mode != kMergeInfoCmModeErase);
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    // kfImportKeepAutoconvVzs cannot be set here.
    const uint32_t output_zst = (pmip->flags / kfPmergeOutputVzs) & 1;
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    if (xheader_entry_ct || write_info_pr_to_pvar) {
      const uintptr_t overflow_buf_size = max_xheader_line_blen + kCompressStreamBlock;
      unsigned char* compress_wkspace = nullptr;
      uint32_t compress_thread_ct = 1;
      if (output_zst) {
        const uintptr_t compress_wkspace_req = CstreamWkspaceReq(overflow_buf_size);
        if (overflow_buf_size + compress_wkspace_req <= linebuf_capacity) {
          compress_wkspace = R_CAST(unsigned char*, &(linebuf[overflow_buf_size]));
        } else {
          if (unlikely(bigstack_alloc_uc(compress_wkspace_req, &compress_wkspace))) {
            goto ScanPvarsAndMergeHeader_ret_NOMEM;
          }
        }
        compress_thread_ct = decompress_thread_ct;
      }
      reterr = InitCstream(outname, 0, output_zst, compress_thread_ct, overflow_buf_size, linebuf, compress_wkspace, &css);
      if (unlikely(reterr)) {
        goto ScanPvarsAndMergeHeader_ret_1;
      }
      cswritep = linebuf;
      uint32_t info_key_ct = 0;
      for (uint32_t xheader_entry_idx = 0; xheader_entry_idx != xheader_entry_ct; ++xheader_entry_idx) {
        const char* key = xheader_entries[xheader_entry_idx];
        const uint32_t key_slen = strlen(key);
        const char* value = &(key[key_slen + 1]);
        if (value[0] == 1) {
          // conflict that we're skipping
          continue;
        }
        if (at_least_one_info_exists) {
          if (StrStartsWithUnsafe(key, "#INFO=<")) {
            ++info_key_ct;
          }
        }
        *cswritep++ = '#';
        cswritep = memcpya(cswritep, key, key_slen);
        cswritep = strcpya(cswritep, value);
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto ScanPvarsAndMergeHeader_ret_WRITE_FAIL;
        }
      }
      if (write_info_pr_to_pvar) {
        cswritep = strcpya_k(cswritep, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">" EOLN_STR);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto ScanPvarsAndMergeHeader_ret_WRITE_FAIL;
        }
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto ScanPvarsAndMergeHeader_ret_WRITE_FAIL;
      }

      BigstackReset(bigstack_mark);
      if (at_least_one_info_exists) {
        // Export INFO key hash table for future use, positioning it at the
        // bottom of bigstack.
        // Similar to ParseInfoHeader() in plink2_data.cc.
        char** xheader_entries_next = S_CAST(char**, bigstack_end_alloc_raw_rd(xheader_entry_ct * sizeof(intptr_t)));
        memmove(xheader_entries_next, xheader_entries, xheader_entry_ct * sizeof(intptr_t));
        xheader_entries = xheader_entries_next;
        const char** info_keys;
        if (unlikely(bigstack_alloc_kcp(info_key_ct + write_info_pr_to_pvar, &info_keys))) {
          goto ScanPvarsAndMergeHeader_ret_NOMEM;
        }
        arena_bottom = g_bigstack_base;
        arena_top = BigstackEndRoundedDown();

        uint32_t xheader_entry_idx = 0;
        for (uint32_t info_key_idx = 0; info_key_idx != info_key_ct; ++info_key_idx) {
          const char* key;
          const char* value;
          uint32_t key_slen;
          do {
            key = xheader_entries[xheader_entry_idx++];
            key_slen = strlen(key);
            value = &(key[key_slen + 1]);
          } while ((value[0] == 1) || (!StrStartsWithUnsafe(key, "#INFO=<")));
          const char* info_key = &(key[strlen("#INFO=<ID=")]);
          const uint32_t info_key_slen = key_slen - strlen("#INFO=<ID=");
          if (info_key_slen > kMaxInfoKeySlen) {
            logerrputs("Error: " PROG_NAME_STR " does not support INFO keys longer than " MAX_INFO_KEY_SLEN_STR " characters.\n");
            reterr = kPglRetNotYetSupported;
            goto ScanPvarsAndMergeHeader_ret_1;
          }
          if (unlikely(value[0] != ',')) {
            goto ScanPvarsAndMergeHeader_ret_MALFORMED_INFO_HEADER_LINE;
          }
          const uintptr_t entry_byte_ct = RoundUpPow2(offsetof(InfoVtype, key) + 1 + info_key_slen, sizeof(intptr_t));
          if (unlikely(S_CAST(uintptr_t, arena_top - arena_bottom) < entry_byte_ct)) {
            goto ScanPvarsAndMergeHeader_ret_NOMEM;
          }
          InfoVtype* new_entry = R_CAST(InfoVtype*, arena_bottom);
          arena_bottom = &(arena_bottom[entry_byte_ct]);
          memcpy(new_entry->key, info_key, info_key_slen + 1);
          info_keys[info_key_idx] = new_entry->key;
          char* numstr;
          uint32_t num_slen;
          if (unlikely(HkvlineFind(&(value[1]), "Number", &numstr, &num_slen) ||
                       FillInfoVtypeNum(numstr, num_slen, info_has_g_keyp, &(new_entry->num)))) {
            goto ScanPvarsAndMergeHeader_ret_MALFORMED_INFO_HEADER_LINE;
          }
        }
        if (write_info_pr_to_pvar) {
          const uintptr_t entry_byte_ct = RoundUpPow2(offsetof(InfoVtype, key) + 1 + strlen("PR"), sizeof(intptr_t));
          if (unlikely(S_CAST(uintptr_t, arena_top - arena_bottom) < entry_byte_ct)) {
            goto ScanPvarsAndMergeHeader_ret_NOMEM;
          }
          InfoVtype* new_entry = R_CAST(InfoVtype*, arena_bottom);
          arena_bottom = &(arena_bottom[entry_byte_ct]);
          strcpy_k(new_entry->key, "PR");
          info_keys[info_key_ct] = new_entry->key;
          new_entry->num = 0;
          ++info_key_ct;
        }
        BigstackBaseSet(arena_bottom);
        arena_bottom = nullptr;
        arena_top = nullptr;
        if (pmip->merge_info_sort == kSortAscii) {
          StrptrArrSortOverread(info_key_ct, info_keys);
        } else if (pmip->merge_info_sort == kSortNatural) {
          StrptrArrNsort(info_key_ct, info_keys);
        }

        const uint32_t info_key_ctl = BitCtToWordCt(info_key_ct);
        uintptr_t* dummy_include;
        if (unlikely(bigstack_end_alloc_w(info_key_ctl, &dummy_include))) {
          goto ScanPvarsAndMergeHeader_ret_NOMEM;
        }
        SetAllBits(info_key_ct, dummy_include);
        const uintptr_t bytes_left = bigstack_left();
        reterr = AllocAndPopulateIdHtableMt(dummy_include, info_keys, info_key_ct, bytes_left - (bytes_left / 64), 1, info_keys_htablep, nullptr, info_keys_htable_sizep, nullptr);
        if (unlikely(reterr)) {
          goto ScanPvarsAndMergeHeader_ret_1;
        }
        *info_keys_ptr = info_keys;
        *info_key_ctp = info_key_ct;

        bigstack_mark = g_bigstack_base;
      }
    } else {
      // bugfix (2 Mar 2021): if there's an existing file with the
      // eventual-output-.pvar name, truncate it down to zero, since we append
      // later.
      FILE* pvar_stub_file = fopen(outname, FOPEN_WB);
      if (unlikely(!pvar_stub_file)) {
        goto ScanPvarsAndMergeHeader_ret_OPEN_FAIL;
      }
      if (unlikely(fclose(pvar_stub_file))) {
        goto ScanPvarsAndMergeHeader_ret_WRITE_FAIL;
      }
    }
    const uint32_t name_ct = cip->name_ct;
    const uint32_t autosome_code_end = cip->autosome_ct + 1;
    const uint32_t name_code_start = autosome_code_end + kChrOffsetCt;
    // When there are multiple chromosomes with no remaining in-edges, we
    // prioritize as follows:
    // 1. Smallest autosomal chromosome index.
    // 2. PAR1 < chrX < PAR2 < chrY < XY < chrM
    // 3. Contig names in natural-sort or ASCII-sort increasing order,
    //    depending on --sort-vars setting.
    // See also SortChr() in plink2_data.cc.  This is a bit simpler since we
    // don't need the sort_idxs to be dense.
    const uint32_t max_code_p1 = cip->max_code + 1;
    const uint32_t chr_code_end = max_code_p1 + name_ct;
    uint32_t* chr_idx_to_sort_idx;
    if (unlikely(bigstack_end_alloc_u32(chr_code_end, &chr_idx_to_sort_idx))) {
      goto ScanPvarsAndMergeHeader_ret_NOMEM;
    }
    for (uint32_t chr_code = 0; chr_code != autosome_code_end; ++chr_code) {
      chr_idx_to_sort_idx[chr_code] = chr_code;
    }
    // bugfix (20 Apr 2025): forgot to add 1
    const uint32_t xymt_ct = max_code_p1 - autosome_code_end;
    const uint32_t xymt_idx_to_chr_sort_offset[kChrOffsetCt] = {1, 3, 4, 5, 0, 2};
    for (uint32_t xymt_idx = 0; xymt_idx != xymt_ct; ++xymt_idx) {
      chr_idx_to_sort_idx[autosome_code_end + xymt_idx] = autosome_code_end + xymt_idx_to_chr_sort_offset[xymt_idx];
    }
    if (name_ct) {
      unsigned char* bigstack_end_mark2 = g_bigstack_end;
      StrSortIndexedDeref* nonstd_sort_buf = S_CAST(StrSortIndexedDeref*, bigstack_end_alloc(name_ct * sizeof(StrSortIndexedDeref)));
      if (unlikely(!nonstd_sort_buf)) {
        goto ScanPvarsAndMergeHeader_ret_NOMEM;
      }
      const char** nonstd_names = cip->nonstd_names;
      for (uint32_t name_idx = 0; name_idx != name_ct; ++name_idx) {
        // bugfix (17 Feb 2022): first non-null nonstd_names[] entry is at
        // [max_code_p1], not [0]
        nonstd_sort_buf[name_idx].strptr = nonstd_names[name_idx + max_code_p1];
        nonstd_sort_buf[name_idx].orig_idx = name_idx;
      }
      // nonstd_names are not allocated in main workspace, so can't overread.
      StrptrArrSortMain(name_ct, 0, (sort_vars_mode != kSortAscii), nonstd_sort_buf);
      for (uint32_t name_idx = 0; name_idx != name_ct; ++name_idx) {
        chr_idx_to_sort_idx[max_code_p1 + nonstd_sort_buf[name_idx].orig_idx] = name_code_start + name_idx;
      }
      BigstackEndReset(bigstack_end_mark2);
    }
    const uint32_t sort_code_end = name_code_start + name_ct;
    const uint32_t sort_code_endl = BitCtToWordCt(sort_code_end);
    uint32_t* chr_sort_idx_to_idx;
    uintptr_t* no_incoming_set;
    if (unlikely(bigstack_end_alloc_u32(sort_code_end, &chr_sort_idx_to_idx) ||
                 bigstack_end_calloc_w(sort_code_endl, &no_incoming_set))) {
      goto ScanPvarsAndMergeHeader_ret_NOMEM;
    }
    const uint32_t chr_code_endl = BitCtToWordCt(chr_code_end);
    const uint32_t chr_ct = PopcountWords(chr_present, chr_code_endl);
    {
      uintptr_t chr_code_base = 0;
      uintptr_t cur_bits = chr_present[0];
      for (uint32_t uii = 0; uii != chr_ct; ++uii) {
        uint32_t chr_code = BitIter1(chr_present, &chr_code_base, &cur_bits);
        const uint32_t chr_sort_idx = chr_idx_to_sort_idx[chr_code];
        chr_sort_idx_to_idx[chr_sort_idx] = chr_code;
        if (!chr_inedge_cts[chr_code]) {
          SetBit(chr_sort_idx, no_incoming_set);
        }
      }
    }
    uint32_t chr_fo_idx = 0;
    while (1) {
      uint32_t chr_sort_idx;
      {
        uint32_t widx = 0;
        for (; widx != sort_code_endl; ++widx) {
          const uintptr_t cur_bits = no_incoming_set[widx];
          if (cur_bits) {
            chr_sort_idx = widx * kBitsPerWord + ctzw(cur_bits);
            // clear this bit
            no_incoming_set[widx] = cur_bits & (cur_bits - 1);
            break;
          }
        }
        if (widx == sort_code_endl) {
          break;
        }
      }
      const uint32_t chr_idx = chr_sort_idx_to_idx[chr_sort_idx];
      cip->chr_file_order[chr_fo_idx] = chr_idx;
      cip->chr_idx_to_foidx[chr_idx] = chr_fo_idx;
      ++chr_fo_idx;
      const uintptr_t* outedges = chr_outedges[chr_idx];
      if (outedges) {
        for (uint32_t widx = 0; widx != chr_code_endl; ++widx) {
          uintptr_t cur_bits = outedges[widx];
          if (cur_bits) {
            const uint32_t other_chr_idx_base = widx * kBitsPerWord;
            do {
              const uint32_t other_chr_idx = other_chr_idx_base + ctzw(cur_bits);
              chr_inedge_cts[other_chr_idx] -= 1;
              if (!chr_inedge_cts[other_chr_idx]) {
                const uint32_t other_chr_sort_idx = chr_idx_to_sort_idx[other_chr_idx];
                SetBit(other_chr_sort_idx, no_incoming_set);
              }
              cur_bits &= cur_bits - 1;
            } while (cur_bits);
          }
        }
      }
    }
    if (unlikely(chr_fo_idx != chr_ct)) {
      logerrputs("Error: Chromosomes are not in a consistent order.  Retry --pmerge[-list] after\nusing --make-pgen/--make-bed + --sort-vars to sort your variants in a\nconsistent manner.\n");
      goto ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT;
    }
    if (unlikely(max_allele_overflow_slen && (!(misc_flags & (kfMiscNewVarIdOverflowMissing | kfMiscNewVarIdOverflowTruncate))))) {
      logputs("\n");
      logerrprintf("Error: Allele code(s) too long for --set-%s-var-ids.\n", (misc_flags & kfMiscSetMissingVarIds)? "missing" : "all");
      // Not a precise bound, but in practice this should print the more useful
      // message >99% of the time.
      if (max_allele_overflow_slen < kMaxIdSlen / 2) {
        logerrprintfww("The longest observed allele code across these datasets has length %u. If you're fine with the corresponding ID length, rerun with \"--new-id-max-allele-len %u\" added to your command line.\n", max_allele_overflow_slen, max_allele_overflow_slen);
      } else {
        logerrprintfww("The longest observed allele code across these datasets has length %u. We recommend deciding on a length-limit, and then adding \"--new-id-max-allele-len <limit> missing\" to your command line to cause variants with longer allele codes to be assigned '.' IDs. (You can then process just those variants with another script, if necessary.)\n", max_allele_overflow_slen);
      }
      goto ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT;
    }
    cip->chr_ct = chr_ct;
    // edge case: if no header meta lines, but we're writing INFO/PR, no need
    // to log "headers merged"
    logprintf("--pmerge%s: %" PRIuPTR " .pvar files scanned%s.\n", pmip->list_fname? "-list" : "", fileset_ct, (xheader_entry_ct || (info_pr_exists && xheader_entries))? ", headers merged" : "");
    if (null_fileset_ct) {
      logprintfww("Note: Ignoring %" PRIuPTR " .pgen file%s intersect the chromosome filter.\n", null_fileset_ct, (null_fileset_ct == 1)? " since it doesn't" : "s since they don't");
      *fileset_ctp -= null_fileset_ct;
    }
    if (info_conflict_ct) {
      logprintfww("Note: Skipping %" PRIuPTR " INFO key%s with mismatching header lines.\n", info_conflict_ct, (info_conflict_ct == 1)? "" : "s");
    }
    *info_conflict_presentp = !!info_conflict_ct;
  }
  while (0) {
  ScanPvarsAndMergeHeader_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ScanPvarsAndMergeHeader_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ScanPvarsAndMergeHeader_ret_TSTREAM_FAIL:
    TextStreamErrPrint(cur_fname, &txs);
    break;
  ScanPvarsAndMergeHeader_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, cur_fname);
  ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  ScanPvarsAndMergeHeader_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ScanPvarsAndMergeHeader_ret_INCONSISTENTLY_PRESENT_CHRSET:
    logerrputs("Error: ##chrSet header line present in some, but not all, --pmerge[-list] input\nfiles.  This is only permitted when the chromosome set is also defined on the\ncommand line.\n");
    reterr = kPglRetInconsistentInput;
    break;
  ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  ScanPvarsAndMergeHeader_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ScanPvarsAndMergeHeader_ret_MALFORMED_INFO_HEADER_LINE:
    logerrputs("Error: Malformed or unrecognized INFO header line.\n");
    reterr = kPglRetMalformedInput;
    break;
  ScanPvarsAndMergeHeader_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 ScanPvarsAndMergeHeader_ret_1:
  CleanupTextStream2(cur_fname, &txs, &reterr);
  CswriteCloseCond(&css, cswritep);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

typedef struct PmergeFilesetSorterStruct {
  PmergeInputFilesetLl* pp;
  uint64_t first_coord;
  uint64_t last_coord;
#ifdef __cplusplus
  bool operator<(const struct PmergeFilesetSorterStruct& rhs) const {
    if (first_coord != rhs.first_coord) {
      return (first_coord < rhs.first_coord);
    }
    return (strcmp(pp->first_varid, rhs.pp->first_varid) < 0);
  }
#endif
} PmergeFilesetSorter;

typedef struct PmergeFilesetNsorterStruct {
  PmergeInputFilesetLl* pp;
  uint64_t first_coord;
  uint64_t last_coord;
#ifdef __cplusplus
  bool operator<(const struct PmergeFilesetNsorterStruct& rhs) const {
    if (first_coord != rhs.first_coord) {
      return (first_coord < rhs.first_coord);
    }
    return (strcmp_natural_uncasted(pp->first_varid, rhs.pp->first_varid) < 0);
  }
#endif
} PmergeFilesetNsorter;

#ifndef __cplusplus
int32_t FilesetAsciiCmp(const void* aa, const void* bb) {
  const PmergeFilesetSorter* pfs1 = S_CAST(const PmergeFilesetSorter*, aa);
  const PmergeFilesetSorter* pfs2 = S_CAST(const PmergeFilesetSorter*, bb);
  const uint64_t first_coord1 = pfs1->first_coord;
  const uint64_t first_coord2 = pfs2->first_coord;
  if (first_coord1 != first_coord2) {
    return (first_coord1 < first_coord2)? -1 : 1;
  }
  return strcmp(pfs1->pp->first_varid, pfs2->pp->first_varid);
}

int32_t FilesetNaturalCmp(const void* aa, const void* bb) {
  const PmergeFilesetNsorter* pfs1 = S_CAST(const PmergeFilesetNsorter*, aa);
  const PmergeFilesetNsorter* pfs2 = S_CAST(const PmergeFilesetNsorter*, bb);
  const uint64_t first_coord1 = pfs1->first_coord;
  const uint64_t first_coord2 = pfs2->first_coord;
  if (first_coord1 != first_coord2) {
    return (first_coord1 < first_coord2)? -1 : 1;
  }
  return strcmp_natural_uncasted(pfs1->pp->first_varid, pfs2->pp->first_varid);
}
#endif

// Determines whether there's no overlap between the positional ranges covered
// by the filesets.  If so, the filesets linked list is sorted.
PglErr DetectConcatJob(const uint32_t* chr_idx_to_foidx, uintptr_t fileset_ct, SortMode sort_vars_mode, PmergeInputFilesetLl** filesets_ptr, uint32_t* is_concat_jobp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    PmergeFilesetSorter* sorted_filesets;
    if (unlikely(BIGSTACK_ALLOC_X(PmergeFilesetSorter, fileset_ct, &sorted_filesets))) {
      goto DetectConcatJob_ret_NOMEM;
    }
    PmergeInputFilesetLl* filesets_iter = *filesets_ptr;
    for (uintptr_t fileset_idx = 0; fileset_idx != fileset_ct; ++fileset_idx, filesets_iter = filesets_iter->next) {
      sorted_filesets[fileset_idx].pp = filesets_iter;
      sorted_filesets[fileset_idx].first_coord = (S_CAST(uint64_t, chr_idx_to_foidx[filesets_iter->first_chr_idx]) << 32) | filesets_iter->first_pos;
      sorted_filesets[fileset_idx].last_coord = (S_CAST(uint64_t, chr_idx_to_foidx[filesets_iter->last_chr_idx]) << 32) | filesets_iter->last_pos;
    }
    const uint32_t sort_ascii = (sort_vars_mode == kSortAscii);
    if (sort_ascii) {
      STD_SORT(fileset_ct, FilesetAsciiCmp, sorted_filesets);
    } else {
      STD_SORT(fileset_ct, FilesetNaturalCmp, R_CAST(PmergeFilesetNsorter*, sorted_filesets));
    }
    uintptr_t prev_last_coord = sorted_filesets[0].last_coord;
    uintptr_t fileset_idx = 1;
    for (; fileset_idx != fileset_ct; ++fileset_idx) {
      const uintptr_t cur_first_coord = sorted_filesets[fileset_idx].first_coord;
      if (cur_first_coord <= prev_last_coord) {
        if (cur_first_coord < prev_last_coord) {
          break;
        }
        if (sort_ascii) {
          if (strcmp(sorted_filesets[fileset_idx - 1].pp->last_varid, sorted_filesets[fileset_idx].pp->first_varid) >= 0) {
            break;
          }
        } else {
          if (strcmp_natural(sorted_filesets[fileset_idx - 1].pp->last_varid, sorted_filesets[fileset_idx].pp->first_varid) >= 0) {
            break;
          }
        }
        continue;
      }
      prev_last_coord = sorted_filesets[fileset_idx].last_coord;
    }
    if (fileset_idx == fileset_ct) {
      PmergeInputFilesetLl** filesets_iterp = filesets_ptr;
      for (fileset_idx = 0; fileset_idx != fileset_ct; ++fileset_idx) {
        *filesets_iterp = sorted_filesets[fileset_idx].pp;
        filesets_iterp = &((*filesets_iterp)->next);
      }
      *filesets_iterp = nullptr;
      logputs("Concatenation job detected.\n");
      *is_concat_jobp = 1;
    }
  }
  while (0) {
  DetectConcatJob_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

void CleanupFilesetLl(PmergeInputFilesetLl* filesets_iter, PglErr* reterrp) {
  PglErr reterr = *reterrp;
  while (filesets_iter != nullptr) {
    free_cond(filesets_iter->first_varid);
    free_cond(filesets_iter->last_varid);
    const uint32_t is_temporary = filesets_iter->is_temporary;
    if (is_temporary) {
      if (filesets_iter->pgen_fname) {
        if (reterr == kPglRetSuccess) {
          if (unlikely(unlink(filesets_iter->pgen_fname))) {
            reterr = kPglRetWriteFail;
          }
        }
        free(filesets_iter->pgen_fname);
      }
      if (filesets_iter->pvar_fname) {
        if (reterr == kPglRetSuccess) {
          if (unlikely(unlink(filesets_iter->pvar_fname))) {
            reterr = kPglRetWriteFail;
          }
        }
        free(filesets_iter->pvar_fname);
      }
      if (filesets_iter->psam_fname) {
        if (reterr == kPglRetSuccess) {
          if (unlikely(unlink(filesets_iter->psam_fname))) {
            reterr = kPglRetWriteFail;
          }
        }
        free(filesets_iter->psam_fname);
      }
      if (filesets_iter->pgen_locked_fname) {
        if (reterr == kPglRetSuccess) {
          if (unlikely(unlink(filesets_iter->pgen_locked_fname))) {
            reterr = kPglRetWriteFail;
          }
        }
        free(filesets_iter->pgen_locked_fname);
      }
    }
    PmergeInputFilesetLl* cur_node = filesets_iter;
    filesets_iter = filesets_iter->next;
    if (is_temporary) {
      free(cur_node);
    }
  }
  *reterrp = reterr;
}

PglErr ScrapeSampleOrder(const char* psam_fname, const SampleIdInfo* siip, const uint32_t* xid_htable, uint32_t read_sample_ct, uint32_t xid_ct, uint32_t xid_htable_size, FamCol fam_cols, uint32_t merge_sids, uint32_t psam_linebuf_capacity, uint32_t max_thread_ct, uint32_t* old_sample_idx_to_new, uint32_t* sample_idx_increasingp, uint32_t* cur_write_sample_ctp, uintptr_t* read_sample_include, uintptr_t* sample_span) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    const uint32_t decompress_thread_ct = MAXV(max_thread_ct - 1, 1);
    char* linebuf;
    if (unlikely(bigstack_alloc_c(psam_linebuf_capacity, &linebuf))) {
      goto ScrapeSampleOrder_ret_NOMEM;
    }
    // Lots of overlap with second half of MergePsams().
    reterr = TextStreamOpenEx(psam_fname, kMaxLongLine, psam_linebuf_capacity, decompress_thread_ct, nullptr, linebuf, &txs);
    if (unlikely(reterr)) {
      goto ScrapeSampleOrder_ret_TSTREAM_REWIND_FAIL_N;
    }
    const char* line_start = TextLineEnd(&txs);
    while (1) {
      if (unlikely(!TextGetUnsafe2K(&txs, &line_start))) {
        reterr = TextStreamRawErrcode(&txs);
        goto ScrapeSampleOrder_ret_TSTREAM_REWIND_FAIL_N;
      }
      if ((line_start[0] != '#') || tokequal_k(&(line_start[1]), "FID") || tokequal_k(&(line_start[1]), "IID")) {
        break;
      }
      line_start = AdvPastDelim(line_start, '\n');
    }
    uint32_t sid_present = 0;
    uint32_t fid_present;
    if (line_start[0] == '#') {
      fid_present = (line_start[1] == 'F');
      const char* postiid_token_start = FirstNonTspace(&(line_start[4]));
      if (fid_present) {
        postiid_token_start = FirstNonTspace(CurTokenEnd(postiid_token_start));
      }
      sid_present = tokequal_k(postiid_token_start, "SID");
      line_start = AdvPastDelim(postiid_token_start, '\n');
    } else {
      // .fam
      fid_present = (fam_cols / kfFamCol1) & 1;
    }
    const uint32_t read_sample_ctl = BitCtToWordCt(read_sample_ct);
    ZeroWArr(read_sample_ctl, read_sample_include);
    const uint32_t xid_ctl = BitCtToWordCt(xid_ct);
    ZeroWArr(xid_ctl, sample_span);
    uint32_t* old_sample_idx_to_new_iter = old_sample_idx_to_new;
    uint32_t prev_write_sample_idx = 0;
    uint32_t sample_idx_increasing = 1;
    for (uint32_t read_sample_idx = 0; read_sample_idx != read_sample_ct; ++read_sample_idx, line_start = AdvPastDelim(line_start, '\n')) {
      if (unlikely(!TextGetUnsafe2K(&txs, &line_start))) {
        reterr = TextStreamRawErrcode(&txs);
        goto ScrapeSampleOrder_ret_TSTREAM_REWIND_FAIL_N;
      }
      uint32_t write_sample_idx;
      if (unlikely(LookupXidHtable(line_start, siip, xid_htable, xid_htable_size, fid_present, sid_present, &write_sample_idx, g_textbuf))) {
        goto ScrapeSampleOrder_ret_REWIND_FAIL_N;
      }
      if (write_sample_idx == UINT32_MAX) {
        if (!merge_sids) {
          continue;
        }
        // in --merge-sids case, initial lookup could fail
        if (unlikely(LookupXidHtable(line_start, siip, xid_htable, xid_htable_size, fid_present, 0, &write_sample_idx, g_textbuf) || (write_sample_idx == UINT32_MAX))) {
          goto ScrapeSampleOrder_ret_REWIND_FAIL_N;
        }
      }
      SetBit(read_sample_idx, read_sample_include);
      SetBit(write_sample_idx, sample_span);
      *old_sample_idx_to_new_iter++ = write_sample_idx;
      if (write_sample_idx < prev_write_sample_idx) {
        sample_idx_increasing = 0;
      }
      prev_write_sample_idx = write_sample_idx;
    }
    const uint32_t cur_write_sample_ct = PopcountWords(read_sample_include, read_sample_ctl);
    *cur_write_sample_ctp = cur_write_sample_ct;
    *sample_idx_increasingp = sample_idx_increasing;
  }
  while (0) {
  ScrapeSampleOrder_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ScrapeSampleOrder_ret_TSTREAM_REWIND_FAIL_N:
    logputs("\n");
    TextStreamErrPrintRewind(psam_fname, &txs, &reterr);
    break;
  ScrapeSampleOrder_ret_REWIND_FAIL_N:
    logputs("\n");
    logerrprintfww(kErrprintfRewind, psam_fname);
    reterr = kPglRetRewindFail;
    break;
  }
  CleanupTextStream2(psam_fname, &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct PvariantMergeContextStruct {
  char* cswritep;
  CompressStreamState css;

  unsigned char write_qual;
  unsigned char write_filter;
  unsigned char write_info;
  unsigned char write_cm;
  MergeQualMode merge_qual_mode;
  MergeFilterMode merge_filter_mode;
  MergeInfoCmMode merge_info_mode;
  MergeInfoCmMode merge_cm_mode;
  SortMode merge_info_sort;
  uint32_t max_allele_ct;

  uintptr_t read_max_allele_ct;
  char** cur_fields_buf;
  uint32_t* cur_slens_buf;
  AlleleCode* allele_cts_buf;
  uintptr_t* is_pr_buf;

  // used for both alleles and FILTER
  const char** tmp_strs;
  uint32_t* tmp_htable;

  uintptr_t* filter_set;

  // concatenation: allele_remap results used immediately for .pgen
  // general case: allele_remap results from .pvar merging pass thrown away,
  //   we re-scrape when merging .pgen
  AlleleCode* allele_remap; // variant-major

  const char* const* info_keys;
  const uint32_t* info_keys_htable;
  uint32_t info_key_ct;
  uint32_t info_keys_htable_size;
  uint32_t info_conflict_present;
  char* info_missing_str;
  char* locked_missing_semicolon_str;
  char* missing_semicolon_str;
  char* locked_missing_comma_str;
  char* missing_comma_str;
  char* pr_str;
  uintptr_t* info_field_set;
  char** basic_info_fields;
  char*** ar_info_fields; // types A and R have multiple subfields
  uint32_t* read_info_field_order;
  uint32_t pr_kidx;

  const char** fnames;
  uintptr_t* line_idx_body_starts;

  // low bit = writing to final, next bit = reading from tmp
  unsigned char tmp_status;

  char* chr_buf;
  uint32_t chr_slen;
  uint32_t cur_bp;
} PvariantMergeContext;

typedef struct PvariantPosMergeContextStruct {
  PvariantMergeContext pmc;

  SortMode sort_vars_mode;
  AlleleCode max_written_allele_ct;
  uint32_t max_written_single_pos_ct;
  uintptr_t max_written_line_blen;
  uintptr_t max_written_pos_blen;

  uintptr_t* write_allele_idx_offsets;
  uintptr_t* write_nonref_flags;
  uint32_t write_variant_idx;
  uint32_t next_print_variant_idx;
  uint32_t write_variant_ct;
} PvariantPosMergeContext;

void PreinitPvariantPosMergeContext(PvariantPosMergeContext* ppmcp) {
  ppmcp->pmc.cswritep = nullptr;
  PreinitCstream(&ppmcp->pmc.css);
}

// this should be more than enough
CONSTI32(kMaxFilterCt, 65535);
#define MAX_FILTER_CT_STR "65535"

PglErr InitPvariantPosMergeContext(const PmergeInfo* pmip, const char* out_fname, const char* const* info_keys, const uint32_t* info_keys_htable, const char** fnames, uintptr_t* line_idx_body_starts, uint32_t write_qual, uint32_t write_filter, uint32_t write_info, uint32_t write_cm, uint32_t info_key_ct, uint32_t info_keys_htable_size, uint32_t info_conflict_present, SortMode sort_vars_mode, uint32_t is_tmp_src, uint32_t is_final_dst, uint32_t output_zst, uintptr_t overflow_buf_size, uint32_t read_max_allele_ct, uint32_t write_max_allele_ct, uintptr_t max_single_pos_ct, uint32_t read_max_nonpass_filter_ct, PvariantPosMergeContext* ppmcp) {
  PvariantMergeContext* pmcp = &ppmcp->pmc;
  PglErr reterr = InitCstreamAlloc(out_fname, is_final_dst, output_zst, 1, overflow_buf_size, &pmcp->css, &pmcp->cswritep);
  if (unlikely(reterr)) {
    return reterr;
  }

  pmcp->write_qual = write_qual;
  pmcp->write_filter = write_filter;
  pmcp->write_info = write_info;
  pmcp->write_cm = write_cm;
  pmcp->merge_qual_mode = pmip->merge_qual_mode;
  MergeFilterMode merge_filter_mode = pmip->merge_filter_mode;
  if (!read_max_nonpass_filter_ct) {
    if ((merge_filter_mode == kMergeFilterModeNonpassUnion) || (merge_filter_mode == kMergeFilterModeNmMatch)) {
      // only need to distinguish PASS vs. missing in this common case
      merge_filter_mode = kMergeFilterModeNmFirst;
    }
  }
  pmcp->merge_filter_mode = merge_filter_mode;
  pmcp->merge_info_mode = pmip->merge_info_mode;
  pmcp->merge_cm_mode = pmip->merge_info_mode;
  pmcp->merge_info_sort = pmip->merge_info_sort;
  const uint32_t max_allele_ct = pmip->max_allele_ct;
  pmcp->max_allele_ct = max_allele_ct;

  char* cswritep = strcpya_k(pmcp->cswritep, "#CHROM\tPOS\tID\tREF\tALT");
  if (write_qual) {
    cswritep = strcpya_k(cswritep, "\tQUAL");
  }
  if (write_filter) {
    cswritep = strcpya_k(cswritep, "\tFILTER");
  }
  if (write_info) {
    cswritep = strcpya_k(cswritep, "\tINFO");
  }
  if (write_cm) {
    cswritep = strcpya_k(cswritep, "\tCM");
  }
  AppendBinaryEoln(&cswritep);
  pmcp->cswritep = cswritep;

  if (unlikely(bigstack_alloc_ac(max_single_pos_ct * read_max_allele_ct, &pmcp->allele_remap))) {
    return kPglRetNomem;
  }
  pmcp->cur_fields_buf = nullptr;
  pmcp->cur_slens_buf = nullptr;
  pmcp->allele_cts_buf = nullptr;
  pmcp->is_pr_buf = nullptr;
  pmcp->tmp_strs = nullptr;
  pmcp->tmp_htable = nullptr;
  pmcp->filter_set = nullptr;
  pmcp->read_max_allele_ct = read_max_allele_ct;
  if (max_single_pos_ct > 1) {
    // Need space for 1 more than the allele-count-limit, since we need to be
    // able to write tombstones.
    const uint32_t merged_allele_ct_max_p1 = max_allele_ct? (max_allele_ct + 1) : (kPglMaxAlleleCt + 1);

    const uint32_t filter_ct_max = MINV(1 + max_single_pos_ct * read_max_nonpass_filter_ct, kMaxFilterCt);
    const uint32_t str_ct_max = MAXV(merged_allele_ct_max_p1, filter_ct_max);
    const uint32_t htable_size = GetHtableMinSize(str_ct_max);
    if (unlikely(bigstack_alloc_cp(max_single_pos_ct, &pmcp->cur_fields_buf) ||
                 bigstack_alloc_u32(max_single_pos_ct, &pmcp->cur_slens_buf) ||
                 bigstack_alloc_ac(max_single_pos_ct, &pmcp->allele_cts_buf) ||
                 bigstack_alloc_w(BitCtToWordCt(max_single_pos_ct), &pmcp->is_pr_buf) ||
                 bigstack_alloc_kcp(str_ct_max, &pmcp->tmp_strs) ||
                 bigstack_alloc_u32(htable_size, &pmcp->tmp_htable))) {
      return kPglRetNomem;
    }
    if (merge_filter_mode == kMergeFilterModeNmMatch) {
      if (unlikely(bigstack_alloc_w(BitCtToWordCt(filter_ct_max), &pmcp->filter_set))) {
        return kPglRetNomem;
      }
    }
  }

  pmcp->info_keys = info_keys;
  pmcp->info_keys_htable = info_keys_htable;
  pmcp->info_key_ct = info_key_ct;
  pmcp->info_keys_htable_size = info_keys_htable_size;
  pmcp->info_conflict_present = info_conflict_present;
  pmcp->info_missing_str = nullptr;
  pmcp->locked_missing_semicolon_str = nullptr;
  pmcp->missing_semicolon_str = nullptr;
  pmcp->locked_missing_comma_str = nullptr;
  pmcp->missing_comma_str = nullptr;
  pmcp->pr_str = nullptr;
  pmcp->info_field_set = nullptr;
  pmcp->basic_info_fields = nullptr;
  pmcp->ar_info_fields = nullptr;
  pmcp->read_info_field_order = nullptr;
  pmcp->pr_kidx = UINT32_MAX;
  if (info_key_ct) {
    char*** ar_info_fields = nullptr;
    if (unlikely(bigstack_alloc_w(BitCtToWordCt(info_key_ct), &pmcp->info_field_set) ||
                 bigstack_calloc_cpp(info_key_ct, &ar_info_fields) ||
                 bigstack_alloc_u32(info_key_ct, &pmcp->read_info_field_order))) {
      return kPglRetNomem;
    }
    char** ar_info_alloc_iter = R_CAST(char**, g_bigstack_base);
    char** ar_info_alloc_limit = R_CAST(char**, BigstackEndRoundedDown());
    uint32_t ar_info_ct = 0;
    int32_t max_num = 1;
    for (uint32_t info_key_idx = 0; info_key_idx != info_key_ct; ++info_key_idx) {
      const int32_t info_vtype = const_container_of(info_keys[info_key_idx], InfoVtype, key)->num;
      if (!IsInfoVtypeARSkip(info_vtype)) {
        if (max_num < info_vtype) {
          max_num = info_vtype;
        }
        continue;
      }
      ++ar_info_ct;
      const uintptr_t subfield_ct = write_max_allele_ct + kInfoVtypeR - info_vtype;
      if (unlikely(S_CAST(uintptr_t, ar_info_alloc_limit - ar_info_alloc_iter) < subfield_ct)) {
        return kPglRetNomem;
      }
      ar_info_fields[info_key_idx] = ar_info_alloc_iter;
      ar_info_alloc_iter = &(ar_info_alloc_iter[subfield_ct]);
    }
    if (ar_info_ct) {
      BigstackBaseSet(ar_info_alloc_iter);
      pmcp->ar_info_fields = ar_info_fields;
    } else {
      BigstackReset(ar_info_fields);
    }
    if (ar_info_ct != info_key_ct) {
      if (unlikely(bigstack_calloc_cp(info_key_ct, &pmcp->basic_info_fields))) {
        return kPglRetNomem;
      }
    }
    if (unlikely(bigstack_alloc_c(2 * max_num, &pmcp->info_missing_str) ||
                 bigstack_alloc_c(strlen("PR;==;=.;=,.,"), &pmcp->pr_str))) {
      return kPglRetNomem;
    }
    u16setsa(pmcp->info_missing_str, 0x2e2c, max_num);
    pmcp->info_missing_str[0] = '=';
    pmcp->locked_missing_semicolon_str = &(pmcp->pr_str[3]);
    pmcp->missing_semicolon_str = &(pmcp->pr_str[6]);
    pmcp->locked_missing_comma_str = &(pmcp->pr_str[9]);
    pmcp->missing_comma_str = &(pmcp->pr_str[11]);

    pmcp->pr_kidx = IdHtableFind("PR", info_keys, info_keys_htable, strlen("PR"), info_keys_htable_size);
  }
  pmcp->fnames = fnames;
  pmcp->line_idx_body_starts = line_idx_body_starts;

  pmcp->tmp_status = is_final_dst + 2 * is_tmp_src;

  ppmcp->sort_vars_mode = sort_vars_mode;
  ppmcp->max_written_allele_ct = 2;
  ppmcp->max_written_single_pos_ct = 1;
  ppmcp->max_written_line_blen = 0;
  ppmcp->max_written_pos_blen = 0;
  ppmcp->write_allele_idx_offsets = nullptr;
  ppmcp->write_nonref_flags = nullptr;
  ppmcp->write_variant_idx = 0;
  ppmcp->next_print_variant_idx = 10000;
  ppmcp->write_variant_ct = 0;
  return kPglRetSuccess;
}

BoolErr FlushPvarCounted(PvariantMergeContext* pmcp, char** cswritepp, uint64_t* cur_line_blenp) {
  if (S_CAST(uintptr_t, (*cswritepp) - pmcp->css.overflow_buf) < kCompressStreamBlock + 1) {
    return 0;
  }
  *cur_line_blenp += (*cswritepp) - pmcp->cswritep;
  BoolErr reterr;
  if (IsUncompressedCstream(&pmcp->css)) {
    reterr = ForceUncompressedCswrite(&pmcp->css, cswritepp);
  } else {
    reterr = ForceCompressedCswrite(&pmcp->css, cswritepp);
  }
  pmcp->cswritep = *cswritepp;
  return reterr;
}

static_assert(sizeof(AlleleCode) == 1, "RenderFinalInfoFromSingleTmp() needs to be updated.");
PglErr RenderFinalInfoFromSingleTmp(uint32_t read_info_slen, PvariantMergeContext* pmcp, char* read_info_start, char** cswritepp) {
  // Since the input is from a tmp file, we can assume keys are already sorted
  // and all values are valid.  We only need to convert '=' locked-missing
  // values to regular missing values.
  //   Flag (0) -> can't happen, there's no notion of a missing value (or a
  //               value at all, for that matter), there are only missing keys
  //   positive number -> that number of '.'; we do not resolve each entry
  //                      separately since we generally want to do the simplest
  //                      thing that might work.  A/R is different since we
  //                      actually know a raw string comparison won't work when
  //                      allele reordering is involved.  (G is punted since
  //                      chrX is a nightmare.)
  //   unknown/G (-1) -> single '.'
  //   A/R -> each allele entry resolved separately
  char* cswritep = *cswritepp;
  *cswritep++ = '\t';
  const char* const* info_keys = pmcp->info_keys;
  const uint32_t* info_keys_htable = pmcp->info_keys_htable;
  const uint32_t info_keys_htable_size = pmcp->info_keys_htable_size;
  char* read_info_end = &(read_info_start[read_info_slen + 1]);
  read_info_end[-1] = ';';
  char* read_info_iter = read_info_start;
  do {
    char* key_end = S_CAST(char*, rawmemchr2(read_info_iter, '=', ';'));
    const uint32_t key_slen = key_end - read_info_iter;
    if (*key_end == ';') {
      // flag
      cswritep = memcpya(cswritep, read_info_iter, key_slen + 1);
      read_info_iter = &(key_end[1]);
      continue;
    }
    const uint32_t kidx = IdHtableFindNnt(read_info_iter, info_keys, info_keys_htable, key_slen, info_keys_htable_size);
    assert(kidx && (kidx != UINT32_MAX));
    const int32_t knum = const_container_of(info_keys[kidx], InfoVtype, key)->num;
    char* value_iter = &(key_end[1]);
    read_info_iter = AdvPastDelim(value_iter, ';');
    if (!IsInfoVtypeARSkip(knum)) {
      if (value_iter[0] != '=') {
        cswritep = memcpya(cswritep, key_end, read_info_iter - key_end);
      } else {
        cswritep = strcpya_k(cswritep, "=.");
        if (knum > 1) {
          // ",."
          cswritep = u16setsa(cswritep, 0x2e2c, knum - 1);
        }
        *cswritep++ = ';';
      }
    } else {
      char* field_end = &(read_info_iter[-1]);
      *field_end = '=';
      while (1) {
        char* next_eq = AdvToDelim(value_iter, '=');
        cswritep = memcpya(cswritep, value_iter, next_eq - value_iter);
        if (next_eq == field_end) {
          break;
        }
        *cswritep++ = '.';
        value_iter = &(next_eq[1]);
      }
      *cswritep++ = ';';
    }
  } while (read_info_iter != read_info_end);
  *cswritepp = &(cswritep[-1]);
  return kPglRetSuccess;
}

// Need to sort by key, but don't need to worry about reading/writing
// locked-missing values or managing Number=A/R/G subfields.
PglErr RenderTmpInfoFromSingleUnsorted(const char* read_pvar_fname, uintptr_t line_idx, uint32_t read_info_slen, uint32_t lone_pgen_pr, PvariantMergeContext* pmcp, char* read_info_start, char** cswritepp, uint32_t* is_pr_ptr) {
  const char* const* info_keys = pmcp->info_keys;
  const uint32_t* info_keys_htable = pmcp->info_keys_htable;
  const uint32_t info_key_ct = pmcp->info_key_ct;
  const uint32_t info_keys_htable_size = pmcp->info_keys_htable_size;
  const uint32_t info_key_ctl = BitCtToWordCt(info_key_ct);
  uintptr_t* info_field_set = pmcp->info_field_set;
  ZeroWArr(info_key_ctl, info_field_set);
  char** basic_info_fields = pmcp->basic_info_fields;
  char* read_info_end = &(read_info_start[read_info_slen + 1]);
  read_info_end[-1] = ';';
  char* read_info_iter = read_info_start;
  do {
    char* key_end = S_CAST(char*, rawmemchr2(read_info_iter, '=', ';'));
    const uint32_t key_slen = key_end - read_info_iter;
    const uint32_t kidx = IdHtableFindNnt(read_info_iter, info_keys, info_keys_htable, key_slen, info_keys_htable_size);
    if (kidx == UINT32_MAX) {
      // Skip this key.  This can happen on e.g. a header conflict.
      if (*key_end != ';') {
        key_end = AdvToDelim(&(key_end[1]), ';');
      }
      read_info_iter = &(key_end[1]);
      continue;
    }
    if (unlikely(IsSet(info_field_set, kidx))) {
      *key_end = '\0';
      logputs("\n");
      logerrprintfww("Error: Duplicate INFO key '%s' on line %" PRIuPTR " of %s.\n", read_info_iter, line_idx, read_pvar_fname);
      return kPglRetMalformedInput;
    }
    SetBit(kidx, info_field_set);
    const int32_t knum = const_container_of(info_keys[kidx], InfoVtype, key)->num;
    if (unlikely((knum == 0) != (*key_end == ';'))) {
      logputs("\n");
      logerrprintfww("Error: Line %" PRIuPTR " of %s has a malformed INFO field (%s immediate ';').\n", line_idx, read_pvar_fname, (knum != 0)? "non-flag key followed by" : "flag key followed by '=' instead of");
      // InconsistentInput instead of MalformedInput error code since this
      // may be caused by inconsistent headers across files
      return kPglRetInconsistentInput;
    }
    // just leave everything semicolon-terminated
    basic_info_fields[kidx] = read_info_iter;
    if (knum == 0) {
      read_info_iter = &(key_end[1]);
    } else {
      char* value_end = AdvToDelim(&(key_end[1]), ';');
      read_info_iter = &(value_end[1]);
    }
  } while (read_info_iter != read_info_end);
  if (lone_pgen_pr) {
    const uint32_t pr_kidx = pmcp->pr_kidx;
    SetBit(pr_kidx, info_field_set);
    basic_info_fields[pr_kidx] = pmcp->pr_str;
  }
  const uint32_t observed_key_ct = PopcountWords(info_field_set, info_key_ctl);
  if (!observed_key_ct) {
    *cswritepp = strcpya_k(*cswritepp, "\t.");
    return kPglRetSuccess;
  }
  char* cswritep = *cswritepp;
  *cswritep++ = '\t';
  uintptr_t kidx_base = 0;
  uintptr_t cur_bits = info_field_set[0];
  for (uint32_t uii = 0; uii != observed_key_ct; ++uii) {
    const uint32_t kidx = BitIter1(info_field_set, &kidx_base, &cur_bits);
    char* kv_start = basic_info_fields[kidx];
    char* kv_end = AdvPastDelim(kv_start, ';');
    cswritep = memcpya(cswritep, kv_start, kv_end - kv_start);
  }
  if (is_pr_ptr) {
    *is_pr_ptr = IsSet(info_field_set, pmcp->pr_kidx);
  }
  *cswritepp = &(cswritep[-1]);
  return kPglRetSuccess;
}

// [0] = REF
// [1] = ALT
// [2] = QUAL
// [3] = FILTER
// [4] = INFO
// [5] = CM (null-terminated string always present, default "0")
typedef struct SamePosPvarRecordStruct {
  uint32_t other_field_offsets[6];
  // High 32 bits = filename index, low bits = variant index
  uint64_t secondary_key;
  AlleleCode allele_ct;
  // low bit = .pgen PR flag; next bit = does .pvar have a INFO/PR field?
  unsigned char pgen_pr_status;
  char variant_id[];
} SamePosPvarRecord;

// returns allele_ct == 0 if allele is filtered out, otherwise at least 2
// if non-null, *is_pr_ptr must be initialized to 0
// record fields assumed to be null-terminated, but this may change in the
// middle of the function
PglErr MergePvariant(uintptr_t merge_rec_ct, PvariantMergeContext* pmcp, SamePosPvarRecord** same_id_records, uint32_t* is_pr_ptr, uint32_t* allele_ct_ptr, uint64_t* cur_line_blenp) {
  char* cswritep = memcpyax(pmcp->cswritep, pmcp->chr_buf, pmcp->chr_slen, '\t');
  cswritep = u32toa_x(pmcp->cur_bp, '\t', cswritep);
  char* cur_variant_id = same_id_records[0]->variant_id;
  cswritep = strcpyax(cswritep, cur_variant_id, '\t');
  const uint32_t max_allele_ct = pmcp->max_allele_ct;
  const uint32_t tmp_status = pmcp->tmp_status;
  AlleleCode* allele_remap = pmcp->allele_remap;
  if (merge_rec_ct == 1) {
    SamePosPvarRecord* cur_record = same_id_records[0];
    const uint32_t read_allele_ct = cur_record->allele_ct;
    if (max_allele_ct && (read_allele_ct > max_allele_ct)) {
      if (tmp_status & 1) {
        // cswritep not saved, so this rolls back the 3-column-write above
        *allele_ct_ptr = 0;
        return kPglRetSuccess;
      }
      // Need to write a tombstone.
      *allele_ct_ptr = max_allele_ct + 1;
      const uint32_t* other_field_offsets = cur_record->other_field_offsets;
      const uint32_t ref_offset = other_field_offsets[0];
      const uint32_t alt_offset = other_field_offsets[1];
      cswritep = memcpya(cswritep, &(cur_variant_id[ref_offset]), alt_offset - ref_offset);
      cswritep[-1] = '\t';
      char* alt_start = &(cur_variant_id[alt_offset]);
      if (read_allele_ct == max_allele_ct + 1) {
        const uint32_t qual_offset = other_field_offsets[2];
        cswritep = memcpya(cswritep, alt_start, qual_offset - alt_offset - 1);
      } else {
        char* alt_trunc = AdvToNthDelim(alt_start, max_allele_ct, ',');
        cswritep = memcpya(cswritep, alt_start, alt_trunc - alt_start);
      }
      // doesn't matter if we write a few extra fields
      cswritep = strcpya_k(cswritep, "\t.\t.\t.\t0" EOLN_STR);
      *cur_line_blenp = cswritep - pmcp->cswritep;
      if (unlikely(Cswrite(&pmcp->css, &cswritep))) {
        return kPglRetWriteFail;
      }
      pmcp->cswritep = cswritep;
      return kPglRetSuccess;
    }
    *allele_ct_ptr = read_allele_ct;
    const uint32_t* other_field_offsets = cur_record->other_field_offsets;
    const uint32_t ref_offset = other_field_offsets[0];
    const uint32_t alt_offset = other_field_offsets[1];
    const uint32_t qual_offset = other_field_offsets[2];
    {
      const char* ref_and_alt_str = &(cur_variant_id[ref_offset]);
      const uint32_t ref_and_alt_slen = qual_offset - 1 - ref_offset;
      cswritep = memcpya(cswritep, ref_and_alt_str, ref_and_alt_slen);
      cswritep[S_CAST(int32_t, alt_offset - qual_offset)] = '\t';
      // only need to initialize first two entries
      // possible todo: detect missing allele codes here so we can validate
      allele_remap[0] = 0;
      allele_remap[1] = 1;
    }
    if (pmcp->write_qual) {
      const uint32_t read_qual_blen = other_field_offsets[3] - qual_offset;
      if (!read_qual_blen) {
        cswritep = strcpya_k(cswritep, "\t.");
      } else {
        char* read_qual_start = &(cur_variant_id[qual_offset]);
        if ((tmp_status != 3) || (!strequal_k_unsafe(read_qual_start, "="))) {
          *cswritep++ = '\t';
          cswritep = memcpya(cswritep, read_qual_start, read_qual_blen - 1);
        } else {
          cswritep = strcpya_k(cswritep, "\t.");
        }
      }
    }
    if (pmcp->write_filter) {
      const uint32_t filter_offset = other_field_offsets[3];
      const uint32_t read_filter_blen = other_field_offsets[4] - filter_offset;
      if (!read_filter_blen) {
        cswritep = strcpya_k(cswritep, "\t.");
      } else {
        char* read_filter_start = &(cur_variant_id[filter_offset]);
        if ((tmp_status != 3) || (!strequal_k_unsafe(read_filter_start, "="))) {
          *cswritep++ = '\t';
          cswritep = memcpya(cswritep, read_filter_start, read_filter_blen - 1);
        } else {
          cswritep = strcpya_k(cswritep, "\t.");
        }
      }
    }
    if (pmcp->write_info) {
      const uint32_t info_offset = other_field_offsets[4];
      const uint32_t read_info_blen = other_field_offsets[5] - info_offset;
      const uint32_t lone_pgen_pr = (cur_record->pgen_pr_status == 1);
      if (is_pr_ptr && (cur_record->pgen_pr_status & 1)) {
        *is_pr_ptr = 1;
        is_pr_ptr = nullptr;
      }
      char* read_info_start = &(cur_variant_id[info_offset]);
      if ((!read_info_blen) || ((read_info_blen == 1) && (read_info_start[0] == '.'))) {
        if (!lone_pgen_pr) {
          cswritep = strcpya_k(cswritep, "\t.");
        } else {
          cswritep = strcpya_k(cswritep, "\tPR");
        }
      } else {
        const uint32_t read_info_slen = read_info_blen - 1;
        if ((pmcp->merge_info_sort == kSortNone) || (tmp_status & 2)) {
          if ((pmcp->merge_info_mode == kMergeInfoCmModeNmFirst) || (tmp_status != 3)) {
            *cswritep++ = '\t';
            if ((!pmcp->info_conflict_present) || (tmp_status & 2)) {
              cswritep = memcpya(cswritep, read_info_start, read_info_slen);
            } else {
              // Line may contain key(s) that were filtered out during .pvar
              // header merge due to a conflict.  Remove the corresponding
              // INFO keys from every line.  This incidentally also removes
              // INFO keys that were (improperly) not declared at all in the
              // header.
              // We don't promise to remove undeclared INFO keys in the general
              // case, since this code path is substantially more expensive
              // than the memcpya() above.  Any user who does want to guarantee
              // removal of such keys (and detection of duplicate keys) can use
              // --merge-info-sort.
              const char* const* info_keys = pmcp->info_keys;
              const uint32_t* info_keys_htable = pmcp->info_keys_htable;
              const uint32_t info_keys_htable_size = pmcp->info_keys_htable_size;
              char* read_info_iter = read_info_start;
              char* read_info_end = &(read_info_start[read_info_slen]);
              char* initial_write_ptr = cswritep;
              *read_info_end++ = ';';
              do {
                char* key_end = S_CAST(char*, rawmemchr2(read_info_iter, '=', ';'));
                const uint32_t key_slen = key_end - read_info_iter;
                const uint32_t kidx = IdHtableFindNnt(read_info_iter, info_keys, info_keys_htable, key_slen, info_keys_htable_size);
                char* value_end = key_end;
                if (*value_end != ';') {
                  value_end = AdvToDelim(&(value_end[1]), ';');
                }
                ++value_end;
                if (kidx != UINT32_MAX) {
                  cswritep = memcpya(cswritep, read_info_iter, value_end - read_info_iter);
                }
                read_info_iter = value_end;
              } while (read_info_iter != read_info_end);
              read_info_end[-1] = '\0';
              if (cswritep == initial_write_ptr) {
                *cswritep++ = '.';
              } else {
                // strip trailing semicolon
                --cswritep;
              }
            }
            if (lone_pgen_pr) {
              cswritep = strcpya_k(cswritep, ";PR");
            } else if (is_pr_ptr != nullptr) {
              *is_pr_ptr = PrInInfo(read_info_slen, read_info_start);
            }
          } else {
            // convert from tmp to final.
            // lone_pgen_pr is impossible here, since we wrote an INFO column
            // to the tmp file, and when we do we'll always insert INFO/PR when
            // relevant.  is_pr_ptr != nullptr is also impossible, for the
            // converse reason.
            PglErr reterr = RenderFinalInfoFromSingleTmp(read_info_slen, pmcp, read_info_start, &cswritep);
            if (unlikely(reterr)) {
              return reterr;
            }
          }
        } else {
          // merge_info_sort != kSortNone, reading from original.
          // So we need to sort by key, but we don't need to worry about
          // reading/writing locked-missing values.
          const uint32_t file_idx = cur_record->secondary_key >> 32;
          const uint32_t variant_uidx = S_CAST(uint32_t, cur_record->secondary_key);
          PglErr reterr = RenderTmpInfoFromSingleUnsorted(pmcp->fnames[file_idx], pmcp->line_idx_body_starts[file_idx] + variant_uidx, read_info_slen, lone_pgen_pr, pmcp, read_info_start, &cswritep, is_pr_ptr);
          if (unlikely(reterr)) {
            return reterr;
          }
        }
      }
    } else if (is_pr_ptr) {
      const uint32_t info_offset = other_field_offsets[4];
      const uint32_t read_info_blen = other_field_offsets[5] - info_offset;
      if (read_info_blen >= 2) {
        char* read_info_start = &(cur_variant_id[info_offset]);
        *is_pr_ptr = PrInInfo(read_info_blen - 1, read_info_start);
      }
    }
    if (pmcp->write_cm) {
      char* cm_start = &(cur_variant_id[other_field_offsets[5]]);
      if ((tmp_status != 3) || (!strequal_k_unsafe(cm_start, "="))) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, cm_start);
      } else {
        cswritep = strcpya_k(cswritep, "\t0");
      }
    }
  } else {
    // 1. Merge alleles.  Fast exit path if any source variant has more than
    //    max_allele_ct alleles on its own.
    // 2. Merge QUAL.
    // 3. Merge FILTER.
    // 4. Merge INFO.
    // 5. Merge CM.
    const uintptr_t read_max_allele_ct = pmcp->read_max_allele_ct;
    AlleleCode* allele_cts = pmcp->allele_cts_buf;
    uint32_t is_pr = 0;
    uint32_t merged_info_allele_ct;
    {
      char** cur_refs = pmcp->cur_fields_buf;
      uintptr_t* is_pr_bitarr = pmcp->is_pr_buf;
      const uintptr_t merge_rec_ctl = BitCtToWordCt(merge_rec_ct);
      ZeroWArr(merge_rec_ctl, is_pr_bitarr);
      // First pass: fill cur_refs, allele_cts, is_pr_bitarr; null-terminate
      // any extra ALT alleles
      for (uintptr_t rec_idx = 0; rec_idx != merge_rec_ct; ++rec_idx) {
        SamePosPvarRecord* cur_record = same_id_records[rec_idx];
        const uint32_t* other_field_offsets = cur_record->other_field_offsets;
        cur_variant_id = cur_record->variant_id;
        cur_refs[rec_idx] = &(cur_variant_id[other_field_offsets[0]]);
        const AlleleCode allele_ct = cur_record->allele_ct;
        allele_cts[rec_idx] = cur_record->allele_ct;
        if (allele_ct > 2) {
          const uint32_t extra_alt_ct = allele_ct - 2;
          char* alt_iter = &(cur_variant_id[other_field_offsets[1]]);
          for (uint32_t uii = 0; uii != extra_alt_ct; ++uii) {
            char* cur_alt_end = AdvToDelim(alt_iter, ',');
            *cur_alt_end = '\0';
            alt_iter = &(cur_alt_end[1]);
          }
        }
        const uint32_t pgen_pr_status = cur_record->pgen_pr_status;
        if (pgen_pr_status & 1) {
          SetBit(rec_idx, is_pr_bitarr);
        } else if (pgen_pr_status) {
          const uint32_t info_offset = other_field_offsets[4];
          const uint32_t info_blen = other_field_offsets[5] - info_offset;
          char* info_start = &(cur_record->variant_id[info_offset]);
          AssignBit(rec_idx, PrInInfo(info_blen - 1, info_start), is_pr_bitarr);
        }
      }
      const uint32_t allele_ct_capacity = max_allele_ct? (max_allele_ct + 1) : (kPglMaxAlleleCt + 1);
      uint32_t allele_ct_limit = 0;
      for (uintptr_t rec_idx = 0; rec_idx != merge_rec_ct; ++rec_idx) {
        allele_ct_limit += allele_cts[rec_idx];
        if (allele_ct_limit >= allele_ct_capacity) {
          allele_ct_limit = allele_ct_capacity;
          break;
        }
      }
      const uint32_t htable_size = GetHtableMinSize(allele_ct_limit);
      uint32_t* merged_alleles_htable = pmcp->tmp_htable;
      SetAllU32Arr(htable_size, merged_alleles_htable);
      const uintptr_t pr_ct = PopcountWords(is_pr_bitarr, merge_rec_ctl);
      const char** merged_alleles = pmcp->tmp_strs;
      uint32_t merged_allele_ct = 0;
      if (pr_ct != merge_rec_ct) {
        // REF allele known.  Verify that we have agreement.
        const uintptr_t known_ref_ct = merge_rec_ct - pr_ct;
        uintptr_t rec_idx_base = 0;
        uintptr_t inv_bits = ~is_pr_bitarr[0];
        uintptr_t rec_idx = BitIter0(is_pr_bitarr, &rec_idx_base, &inv_bits);
        char* ref_allele = cur_refs[rec_idx];
        const uint32_t ref_allele_slen = strlen(ref_allele);
        allele_remap[rec_idx * read_max_allele_ct] = 0;
        for (uintptr_t ulii = 1; ulii != known_ref_ct; ++ulii) {
          rec_idx = BitIter0(is_pr_bitarr, &rec_idx_base, &inv_bits);
          allele_remap[rec_idx * read_max_allele_ct] = 0;
          char* cur_ref = cur_refs[rec_idx];
          if (unlikely(!strequal_unsafe(cur_ref, ref_allele, ref_allele_slen))) {
            logputs("\n");
            logerrprintfww("Error: Conflicting REF alleles for variant '%s' at %s:%u.\n", cur_variant_id, pmcp->chr_buf, pmcp->cur_bp);
            return kPglRetInconsistentInput;
          }
        }
        merged_alleles[0] = ref_allele;
        HtableAddNondup(ref_allele, ref_allele_slen, htable_size, 0, merged_alleles_htable);
        ++merged_allele_ct;
      } else {
        is_pr = 1;
        if (is_pr_ptr) {
          *is_pr_ptr = 1;
        }
      }
      for (uintptr_t rec_idx = 0; rec_idx != merge_rec_ct; ++rec_idx) {
        char* allele_iter = cur_refs[rec_idx];
        uint32_t allele_idx = 0;
        if (!IsSet(is_pr_bitarr, rec_idx)) {
          allele_idx = 1;
          allele_iter = &(strnul(allele_iter)[1]);
        }
        const uint32_t allele_ct = allele_cts[rec_idx];
        AlleleCode* cur_allele_remap = &(allele_remap[rec_idx * read_max_allele_ct]);
        for (; allele_idx != allele_ct; ++allele_idx) {
          const uint32_t cur_allele_slen = strlen(allele_iter);
          uint32_t merged_allele_idx;
          if ((allele_iter[0] != '.') || (cur_allele_slen != 1)) {
            merged_allele_idx = IdHtableAdd(allele_iter, merged_alleles, cur_allele_slen, htable_size, merged_allele_ct, merged_alleles_htable);
            if (merged_allele_idx == UINT32_MAX) {
              merged_alleles[merged_allele_ct] = allele_iter;
              merged_allele_idx = merged_allele_ct;
              ++merged_allele_ct;
              if (merged_allele_ct == allele_ct_limit) {
                // Too many alleles.
                if (unlikely(allele_ct_limit == kPglMaxAlleleCt + 1)) {
                  logputs("\n");
                  logerrprintfww("Error: Too many alleles for variant '%s' at %s:%u. (This " PROG_NAME_STR " build is limited to " PGL_MAX_ALLELE_CT_STR ".)\n", cur_variant_id, pmcp->chr_buf, pmcp->cur_bp);
                  return kPglRetNotYetSupported;
                }
                if (tmp_status & 1) {
                  *allele_ct_ptr = 0;
                  return kPglRetSuccess;
                }
                // Write tombstone.
                *allele_ct_ptr = allele_ct_limit;
                cswritep = strcpyax(cswritep, merged_alleles[0], '\t');
                for (uint32_t aidx = 1; aidx != allele_ct_limit; ++aidx) {
                  cswritep = strcpyax(cswritep, merged_alleles[aidx], ',');
                }
                cswritep[-1] = '\t';
                cswritep = strcpya_k(cswritep, ".\t.\t.\t0" EOLN_STR);
                *cur_line_blenp = cswritep - pmcp->cswritep;
                if (unlikely(Cswrite(&pmcp->css, &cswritep))) {
                  return kPglRetWriteFail;
                }
                pmcp->cswritep = cswritep;
                return kPglRetSuccess;
              }
            }
          } else {
            if (unlikely(allele_ct != 2)) {
              SamePosPvarRecord* cur_record = same_id_records[rec_idx];
              const uint32_t file_idx = cur_record->secondary_key >> 32;
              const uint32_t variant_uidx = S_CAST(uint32_t, cur_record->secondary_key);
              logputs("\n");
              logerrprintfww("Error: Missing allele in multiallelic variant on line %" PRIuPTR " of %s.\n", pmcp->line_idx_body_starts[file_idx] + variant_uidx, pmcp->fnames[file_idx]);
              return kPglRetMalformedInput;
            }
            merged_allele_idx = kMissingAlleleCode;
          }
          cur_allele_remap[allele_idx] = merged_allele_idx;
          allele_iter = &(allele_iter[cur_allele_slen + 1]);
        }
      }
      // When merging .ped-derived datasets, it's possible for *all* alleles to
      // be missing.
      merged_info_allele_ct = merged_allele_ct;
      for (; merged_allele_ct < 2; ++merged_allele_ct) {
        merged_alleles[merged_allele_ct] = &(g_one_char_strs[92]);
      }
      cswritep = strcpyax(cswritep, merged_alleles[0], '\t');
      for (uint32_t aidx = 1; aidx != merged_allele_ct; ++aidx) {
        cswritep = strcpyax(cswritep, merged_alleles[aidx], ',');
      }
      --cswritep;
      *allele_ct_ptr = merged_allele_ct;
    }
    if (pmcp->write_qual) {
      const MergeQualMode merge_qual_mode = pmcp->merge_qual_mode;
      *cswritep++ = '\t';
      if ((merge_qual_mode == kMergeQualModeFirst) || (merge_qual_mode == kMergeQualModeNmFirst)) {
        uintptr_t rec_idx = 0;
        for (; rec_idx != merge_rec_ct; ++rec_idx) {
          SamePosPvarRecord* cur_record = same_id_records[rec_idx];
          const uint32_t* other_field_offsets = cur_record->other_field_offsets;
          char* tmp_variant_id = cur_record->variant_id;
          const uint32_t qual_start = other_field_offsets[2];
          const uint32_t qual_blen = other_field_offsets[3] - qual_start;
          if (qual_blen) {
            if ((!((qual_blen == 2) && (tmp_variant_id[qual_start] == '.'))) || (merge_qual_mode == kMergeQualModeFirst)) {
              cswritep = memcpya(cswritep, &(tmp_variant_id[qual_start]), qual_blen - 1);
              break;
            }
          }
        }
        if (rec_idx == merge_rec_ct) {
          *cswritep++ = '.';
        }
      } else {
        char** cur_qual_strs = pmcp->cur_fields_buf;
        uintptr_t nm_qual_ct = 0;
        for (uintptr_t rec_idx = 0; rec_idx != merge_rec_ct; ++rec_idx) {
          SamePosPvarRecord* cur_record = same_id_records[rec_idx];
          const uint32_t* other_field_offsets = cur_record->other_field_offsets;
          char* tmp_variant_id = cur_record->variant_id;
          const uint32_t qual_start = other_field_offsets[2];
          const uint32_t qual_blen = other_field_offsets[3] - qual_start;
          if (qual_blen && (!((qual_blen == 2) && (tmp_variant_id[qual_start] == '.')))) {
            cur_qual_strs[nm_qual_ct] = &(tmp_variant_id[qual_start]);
            ++nm_qual_ct;
          }
        }
        if ((!nm_qual_ct) || ((merge_qual_mode == kMergeQualModeMin) && (nm_qual_ct != merge_rec_ct))) {
          *cswritep++ = '.';
        } else if (nm_qual_ct == 1) {
          cswritep = strcpya(cswritep, cur_qual_strs[0]);
        } else {
          // Must parse and compare QUAL values.
          float main_qual = 0.0;
          uintptr_t nm_rec_idx = 0;
          for (; nm_rec_idx != nm_qual_ct; ++nm_rec_idx) {
            float cur_qual;
            if (unlikely(ScanFloat(cur_qual_strs[nm_rec_idx], &cur_qual))) {
              // Note that we only validate QUAL when we have to.

              // Yes, line number/filename would be better, but to get those
              // we'd need to determine the original variant index from the
              // nonmissing index, which isn't really worth the trouble here.
              logputs("\n");
              logerrprintfww("Error: Invalid QUAL value for variant '%s' at %s:%u.\n", cur_variant_id, pmcp->chr_buf, pmcp->cur_bp);
              return kPglRetMalformedInput;
            }
            if (!nm_rec_idx) {
              main_qual = cur_qual;
              continue;
            }
            if (merge_qual_mode == kMergeQualModeMin) {
              if (cur_qual < main_qual) {
                main_qual = cur_qual;
              }
            } else if (cur_qual != main_qual) {
              break;
            }
          }
          if (nm_rec_idx == nm_qual_ct) {
            cswritep = ftoa_g(main_qual, cswritep);
          } else {
            *cswritep++ = '.';
          }
        }
      }
    }
    if (pmcp->write_filter) {
      const MergeFilterMode merge_filter_mode = pmcp->merge_filter_mode;
      *cswritep++ = '\t';
      if ((merge_filter_mode == kMergeFilterModeFirst) || (merge_filter_mode == kMergeFilterModeNmFirst)) {
        uintptr_t rec_idx = 0;
        for (; rec_idx != merge_rec_ct; ++rec_idx) {
          SamePosPvarRecord* cur_record = same_id_records[rec_idx];
          const uint32_t* other_field_offsets = cur_record->other_field_offsets;
          char* tmp_variant_id = cur_record->variant_id;
          const uint32_t filter_start = other_field_offsets[3];
          const uint32_t filter_blen = other_field_offsets[4] - filter_start;
          if (filter_blen) {
            if ((!((filter_blen == 2) && (tmp_variant_id[filter_start] == '.'))) || (merge_filter_mode == kMergeFilterModeFirst)) {
              cswritep = memcpya(cswritep, &(tmp_variant_id[filter_start]), filter_blen - 1);
              break;
            }
          }
        }
        if (rec_idx == merge_rec_ct) {
          *cswritep++ = '.';
        }
      } else {
        // First pass: count number of tokens in non-PASS entries so we can
        // size the hash table properly, distinguish between PASS and missing
        // if there are no non-PASS entries.
        char** cur_filter_strs = pmcp->cur_fields_buf;
        uintptr_t nm_merge_rec_ct = 0;
        uintptr_t nonpass_subtoken_ct = 0;
        uint32_t pass_exists = 0;
        for (uintptr_t rec_idx = 0; rec_idx != merge_rec_ct; ++rec_idx) {
          SamePosPvarRecord* cur_record = same_id_records[rec_idx];
          const uint32_t* other_field_offsets = cur_record->other_field_offsets;
          char* tmp_variant_id = cur_record->variant_id;
          const uint32_t filter_start = other_field_offsets[3];
          const uint32_t filter_blen = other_field_offsets[4] - filter_start;
          if (filter_blen && (!((filter_blen == 2) && (tmp_variant_id[filter_start] == '.')))) {
            char* filter_str = &(tmp_variant_id[filter_start]);
            // bugfix (7 Jan 2023)
            if ((filter_blen == 5) && memequal_sk(filter_str, "PASS")) {
              pass_exists = 1;
            } else {
              nonpass_subtoken_ct += 1 + CountByte(filter_str, ';', filter_blen - 1);
              cur_filter_strs[nm_merge_rec_ct] = filter_str;
              ++nm_merge_rec_ct;
            }
          }
        }
        if (!nm_merge_rec_ct) {
          if (pass_exists) {
            cswritep = strcpya_k(cswritep, "PASS");
          } else {
            *cswritep++ = '.';
          }
        } else if (nm_merge_rec_ct == 1) {
          cswritep = strcpya(cswritep, cur_filter_strs[0]);
        } else {
          const uint32_t htable_size = GetHtableMinSize(MINV(nonpass_subtoken_ct + 1, kMaxFilterCt));
          uint32_t* filter_htable = pmcp->tmp_htable;
          SetAllU32Arr(htable_size, filter_htable);
          const char** subtokens = pmcp->tmp_strs;
          subtokens[0] = "PASS";
          HtableAddNondup(subtokens[0], strlen("PASS"), htable_size, 0, filter_htable);
          uint32_t cur_filter_ct = 1;

          const uintptr_t first_loop_stop = (merge_filter_mode == kMergeFilterModeNmMatch)? 1 : nm_merge_rec_ct;
          for (uintptr_t nm_rec_idx = 0; nm_rec_idx != first_loop_stop; ++nm_rec_idx) {
            char* filter_subtoken_iter = cur_filter_strs[nm_rec_idx];
            while (1) {
              char* subtoken_end = Strchrnul(filter_subtoken_iter, ';');
              const uint32_t is_last = (*subtoken_end == '\0');
              if (!is_last) {
                *subtoken_end = '\0';
              }
              if (IdHtableAdd(filter_subtoken_iter, subtokens, subtoken_end - filter_subtoken_iter, htable_size, cur_filter_ct, filter_htable) == UINT32_MAX) {
                if (unlikely(cur_filter_ct == kMaxFilterCt + 1)) {
                  logputs("\n");
                  logerrprintfww("Error: Too many distinct FILTERs for variant '%s' at %s:%u. (This " PROG_NAME_STR " build is limited to " MAX_FILTER_CT_STR ".)\n", cur_variant_id, pmcp->chr_buf, pmcp->cur_bp);
                  return kPglRetNotYetSupported;
                }
                subtokens[cur_filter_ct] = filter_subtoken_iter;
                ++cur_filter_ct;
              }
              if (is_last) {
                break;
              }
              filter_subtoken_iter = &(subtoken_end[1]);
            }
          }
          uint32_t mismatch_found = 0;
          if (merge_filter_mode == kMergeFilterModeNmMatch) {
            uintptr_t* filter_set = pmcp->filter_set;
            const uint32_t filter_ctl = BitCtToWordCt(cur_filter_ct);
            for (uintptr_t nm_rec_idx = 1; nm_rec_idx != nm_merge_rec_ct; ++nm_rec_idx) {
              filter_set[0] = 1;
              ZeroWArr(filter_ctl - 1, &(filter_set[1]));
              char* filter_subtoken_iter = cur_filter_strs[nm_rec_idx];
              while (1) {
                char* subtoken_end = Strchrnul(filter_subtoken_iter, ';');
                const uint32_t is_last = (*subtoken_end == '\0');
                if (!is_last) {
                  *subtoken_end = '\0';
                }
                const uint32_t filter_idx = IdHtableFind(filter_subtoken_iter, subtokens, filter_htable, subtoken_end - filter_subtoken_iter, htable_size);
                if (filter_idx == UINT32_MAX) {
                  mismatch_found = 1;
                  break;
                }
                SetBit(filter_idx, filter_set);
                if (is_last) {
                  break;
                }
                filter_subtoken_iter = &(subtoken_end[1]);
              }
              if (mismatch_found) {
                break;
              }
              if (PopcountWords(filter_set, filter_ctl) != cur_filter_ct) {
                mismatch_found = 1;
                break;
              }
            }
          }
          if (!mismatch_found) {
            if (unlikely((subtokens[1][0] == '\0') && (cur_filter_ct == 2))) {
              // We usually don't explicitly check for this type of malformed
              // input, but if the empty string is the *only* non-PASS FILTER
              // key, we must error out since we'd otherwise generate an
              // invalid .pvar.
              logputs("\n");
              logerrprintfww("Error: Length-0 FILTER key present for variant '%s' at %s:%u.\n", cur_variant_id, pmcp->chr_buf, pmcp->cur_bp);
              return kPglRetMalformedInput;
            }
            for (uint32_t uii = 1; uii != cur_filter_ct; ++uii) {
              cswritep = strcpyax(cswritep, subtokens[uii], ';');
            }
            --cswritep;
          }
        }
      }
    }
    if (pmcp->write_info) {
      *cswritep++ = '\t';
      // Lots of overlap with RenderTmpInfoFromSingleUnsorted().
      const char* const* info_keys = pmcp->info_keys;
      const uint32_t* info_keys_htable = pmcp->info_keys_htable;
      const uint32_t info_key_ct = pmcp->info_key_ct;
      const uint32_t info_keys_htable_size = pmcp->info_keys_htable_size;
      const uint32_t info_key_ctl = BitCtToWordCt(info_key_ct);
      const char* info_missing_str = pmcp->info_missing_str;
      uintptr_t* info_field_set = pmcp->info_field_set;
      ZeroWArr(info_key_ctl, info_field_set);
      char** basic_info_fields = pmcp->basic_info_fields;
      char*** ar_info_fields = pmcp->ar_info_fields;
      uint32_t* read_info_field_order_iter = (pmcp->merge_info_sort == kSortNone)? pmcp->read_info_field_order : nullptr;
      const MergeInfoCmMode merge_info_mode = pmcp->merge_info_mode;
      char* locked_missing_semicolon_str;
      char* locked_missing_comma_str;
      if (tmp_status == 1) {
        locked_missing_semicolon_str = pmcp->missing_semicolon_str;
        locked_missing_comma_str = pmcp->missing_comma_str;
      } else {
        locked_missing_semicolon_str = pmcp->locked_missing_semicolon_str;
        locked_missing_comma_str = pmcp->locked_missing_comma_str;
      }
      uintptr_t rec_idx = 0;
      for (; rec_idx != merge_rec_ct; ++rec_idx) {
        SamePosPvarRecord* cur_record = same_id_records[rec_idx];
        const uint32_t* other_field_offsets = cur_record->other_field_offsets;
        char* tmp_variant_id = cur_record->variant_id;
        const uint32_t info_start = other_field_offsets[4];
        const uint32_t info_blen = other_field_offsets[5] - info_start;
        if ((!info_blen) || ((info_blen == 2) && (tmp_variant_id[info_start] == '.'))) {
          continue;
        }
        char* read_info_iter = &(tmp_variant_id[info_start]);
        char* read_info_end = &(read_info_iter[info_blen]);
        read_info_end[-1] = ';';

        uint32_t read_allele_ct = allele_cts[rec_idx];
        AlleleCode* cur_allele_remap = &(allele_remap[rec_idx * read_max_allele_ct]);
        if (read_allele_ct == 2) {
          if (cur_allele_remap[0] == kMissingAlleleCode) {
            read_allele_ct = 0;
          } else if (cur_allele_remap[1] == kMissingAlleleCode) {
            read_allele_ct = 1;
          }
        }

        do {
          char* key_end = S_CAST(char*, rawmemchr2(read_info_iter, '=', ';'));
          const uint32_t key_slen = key_end - read_info_iter;
          const uint32_t kidx = IdHtableFindNnt(read_info_iter, info_keys, info_keys_htable, key_slen, info_keys_htable_size);
          if (kidx == UINT32_MAX) {
            // Skip this key.  This can happen on e.g. a header conflict.
            if (*key_end != ';') {
              key_end = AdvToDelim(&(key_end[1]), ';');
            }
            read_info_iter = &(key_end[1]);
            continue;
          }
          // possible todo: error out on duplicated INFO key in a single
          // input variant, like we already do in
          // RenderTmpInfoFromSingleUnsorted().
          const int32_t knum = const_container_of(info_keys[kidx], InfoVtype, key)->num;
          if (unlikely((knum == 0) != (*key_end == ';'))) {
            const uint32_t file_idx = cur_record->secondary_key >> 32;
            const uint32_t variant_uidx = S_CAST(uint32_t, cur_record->secondary_key);
            logputs("\n");
            logerrprintfww("Error: Line %" PRIuPTR " of %s has a malformed INFO field (%s immediate ';').\n", pmcp->line_idx_body_starts[file_idx] + variant_uidx, pmcp->fnames[file_idx], (knum != 0)? "non-flag key followed by" : "flag key followed by '=' instead of");
            // InconsistentInput instead of MalformedInput error code since
            // this may be caused by inconsistent headers across files
            return kPglRetInconsistentInput;
          }
          if (!IsInfoVtypeARSkip(knum)) {
            // Treat payload as a single string, including '=' and ';'.  Do not
            // null-terminate, semicolon is enough.
            if (knum == 0) {
              if (!IsSet(info_field_set, kidx)) {
                SetBit(kidx, info_field_set);
                basic_info_fields[kidx] = key_end;
                if (read_info_field_order_iter) {
                  *read_info_field_order_iter++ = kidx;
                }
              }
              read_info_iter = &(key_end[1]);
            } else {
              char* value_end = AdvToDelim(&(key_end[1]), ';');
              if (!IsSet(info_field_set, kidx)) {
                SetBit(kidx, info_field_set);
                if (read_info_field_order_iter) {
                  *read_info_field_order_iter++ = kidx;
                }
                if (merge_info_mode == kMergeInfoCmModeFirst) {
                  basic_info_fields[kidx] = key_end;
                } else {
                  const uint32_t value_slen = value_end - key_end;
                  basic_info_fields[kidx] = nullptr;
                  if (knum <= 1) {
                    // Regular missing string is "=."
                    if ((key_end[1] != '.') || (value_slen != 2)) {
                      basic_info_fields[kidx] = key_end;
                    }
                  } else {
                    // Missing string looks like "=.,.,.", where number of dots
                    // is knum.
                    if ((value_slen != S_CAST(uint32_t, knum) * 2) || (!memequal(key_end, info_missing_str, value_slen))) {
                      basic_info_fields[kidx] = key_end;
                    }
                  }
                }
              } else if (merge_info_mode != kMergeInfoCmModeFirst) {
                uint32_t is_nonmissing = 0;
                // Locked-missing string is "==" (or "=." when we don't need to
                // generate a temporary .pvar).  This is the only way for
                // basic_info_fields[kidx] == locked_missing_semicolon_str to
                // be true.
                const uint32_t value_slen = value_end - key_end;
                if ((!basic_info_fields[kidx]) || ((merge_info_mode == kMergeInfoCmModeNmMatch) && (!memequal_k(basic_info_fields[kidx], locked_missing_semicolon_str, 3)))) {
                  if (knum <= 1) {
                    is_nonmissing = (key_end[1] != '.') || (value_slen != 2);
                  } else {
                    is_nonmissing = (value_slen != S_CAST(uint32_t, knum) * 2) || (!memequal(key_end, info_missing_str, value_slen));
                  }
                }
                if (is_nonmissing) {
                  if (!basic_info_fields[kidx]) {
                    basic_info_fields[kidx] = key_end;
                  } else {
                    // set to locked-missing if unequal
                    if (!memequal(key_end, basic_info_fields[kidx], value_slen + 1)) {
                      // This may need to be postprocessed in the knum>1 case.
                      basic_info_fields[kidx] = locked_missing_semicolon_str;
                    }
                  }
                }
              }
              read_info_iter = &(value_end[1]);
            }
          } else {
            // Treat per-allele payloads separately.
            const uint32_t first_aidx = knum - kInfoVtypeR;
            if (unlikely(first_aidx >= read_allele_ct)) {
              const uint32_t file_idx = cur_record->secondary_key >> 32;
              const uint32_t variant_uidx = S_CAST(uint32_t, cur_record->secondary_key);
              logputs("\n");
              logerrprintfww("Error: Invalid INFO key '%s' on line %" PRIuPTR " of %s (this variant has no %s allele, and key has Number=%c).\n", info_keys[kidx], pmcp->line_idx_body_starts[file_idx] + variant_uidx, pmcp->fnames[file_idx], read_allele_ct? "ALT" : "REF", first_aidx? 'A' : 'R');
              return kPglRetInconsistentInput;
            }
            char** cur_ar_info_fields = ar_info_fields[kidx];
            const uint32_t subfield_ct = merged_info_allele_ct - first_aidx;
            if (!IsSet(info_field_set, kidx)) {
              SetBit(kidx, info_field_set);
              ZeroPtrArr(subfield_ct, cur_ar_info_fields);
              if (read_info_field_order_iter) {
                *read_info_field_order_iter++ = kidx;
              }
            }
            read_info_iter = &(key_end[1]);
            char* value_end = AdvToDelim(read_info_iter, ';');
            if (unlikely(CountByte(read_info_iter, ',', value_end - read_info_iter) != read_allele_ct - first_aidx - 1)) {
              const uint32_t file_idx = cur_record->secondary_key >> 32;
              const uint32_t variant_uidx = S_CAST(uint32_t, cur_record->secondary_key);
              logputs("\n");
              logerrprintfww("Error: Invalid INFO entry for key '%s' on line %" PRIuPTR " of %s (incorrect number of values).\n", info_keys[kidx], pmcp->line_idx_body_starts[file_idx] + variant_uidx, pmcp->fnames[file_idx]);
              return kPglRetInconsistentInput;
            }
            *value_end = ',';
            uint32_t read_aidx = first_aidx;
            if (merge_info_mode == kMergeInfoCmModeNmFirst) {
              do {
                const uint32_t write_aidx = cur_allele_remap[read_aidx] - first_aidx;
                if ((!cur_ar_info_fields[write_aidx]) && (!memequal_sk(read_info_iter, ".,"))) {
                  cur_ar_info_fields[write_aidx] = read_info_iter;
                }
                read_info_iter = AdvPastDelim(read_info_iter, ',');
                ++read_aidx;
              } while (read_aidx != read_allele_ct);
            } else if (merge_info_mode == kMergeInfoCmModeFirst) {
              do {
                const uint32_t write_aidx = cur_allele_remap[read_aidx] - first_aidx;
                cur_ar_info_fields[write_aidx] = read_info_iter;
                read_info_iter = AdvPastDelim(read_info_iter, ',');
                ++read_aidx;
              } while (read_aidx != read_allele_ct);
            } else {
              // NmMatch
              do {
                const uint32_t write_aidx = cur_allele_remap[read_aidx] - first_aidx;
                char* next_subtoken_start = AdvPastDelim(read_info_iter, ',');
                if (!memequal_sk(read_info_iter, ".,")) {
                  if (!cur_ar_info_fields[write_aidx]) {
                    cur_ar_info_fields[write_aidx] = read_info_iter;
                  } else if (!memequal_k(cur_ar_info_fields[write_aidx], locked_missing_comma_str, 2)) {
                    if (!memequal(read_info_iter, cur_ar_info_fields[write_aidx], next_subtoken_start - read_info_iter)) {
                      cur_ar_info_fields[write_aidx] = locked_missing_comma_str;
                    }
                  }
                }
                read_info_iter = next_subtoken_start;
                ++read_aidx;
              } while (read_aidx != read_allele_ct);
            }
          }
        } while (read_info_iter != read_info_end);
      }
      if (is_pr) {
        const uint32_t pr_kidx = pmcp->pr_kidx;
        if (!IsSet(info_field_set, pr_kidx)) {
          SetBit(pr_kidx, info_field_set);
          basic_info_fields[pr_kidx] = pmcp->pr_str;
          if (read_info_field_order_iter) {
            *read_info_field_order_iter++ = pr_kidx;
          }
        }
      }
      const uint32_t observed_key_ct = PopcountWords(info_field_set, info_key_ctl);
      if (!observed_key_ct) {
        *cswritep++ = '.';
      } else {
        uint32_t* read_info_field_order = pmcp->read_info_field_order;
        if (!read_info_field_order_iter) {
          // We can now fill in the sorted field order.
          uintptr_t kidx_base = 0;
          uintptr_t cur_bits = info_field_set[0];
          for (uint32_t uii = 0; uii != observed_key_ct; ++uii) {
            const uint32_t kidx = BitIter1(info_field_set, &kidx_base, &cur_bits);
            read_info_field_order[uii] = kidx;
          }
        }
        if ((merge_info_mode != kMergeInfoCmModeNmMatch) || (tmp_status != 3)) {
          for (uint32_t uii = 0; uii != observed_key_ct; ++uii) {
            const uint32_t kidx = read_info_field_order[uii];
            const char* cur_key = info_keys[kidx];
            const int32_t knum = const_container_of(cur_key, InfoVtype, key)->num;
            cswritep = strcpya(cswritep, cur_key);
            if (!IsInfoVtypeARSkip(knum)) {
              char* v_start = basic_info_fields[kidx];
              if (!v_start) {
                // Missing value.
                if (knum <= 1) {
                  cswritep = strcpya_k(cswritep, "=.;");
                } else {
                  // Need to expand to "=.,.,.;"
                  *cswritep++ = '=';
                  cswritep = u16setsa(cswritep, 0x2c2e, knum);
                  cswritep[-1] = ';';
                }
              } else if ((tmp_status != 1) || (knum <= 1) || (!memequal_k(v_start, locked_missing_semicolon_str, 3))) {
                char* v_end = AdvPastDelim(v_start, ';');
                cswritep = memcpya(cswritep, v_start, v_end - v_start);
              } else {
                // Need to expand to "=.,.,.;"
                *cswritep++ = '=';
                cswritep = u16setsa(cswritep, 0x2c2e, knum);
                cswritep[-1] = ';';
              }
            } else {
              *cswritep++ = '=';
              char** cur_ar_info_fields = ar_info_fields[kidx];
              const uint32_t first_aidx = knum - kInfoVtypeR;
              for (uint32_t write_idx = 0; write_idx != merged_info_allele_ct - first_aidx; ++write_idx) {
                char* v_start = cur_ar_info_fields[write_idx];
                if (!v_start) {
                  cswritep = strcpya_k(cswritep, ".,");
                } else {
                  char* v_end = AdvPastDelim(v_start, ',');
                  cswritep = memcpya(cswritep, v_start, v_end - v_start);
                }
              }
              cswritep[-1] = ';';
            }
          }
        } else {
          // Necessary to convert all intermediate locked-missing
          // representations to final form.
          for (uint32_t uii = 0; uii != observed_key_ct; ++uii) {
            const uint32_t kidx = read_info_field_order[uii];
            const char* cur_key = info_keys[kidx];
            const int32_t knum = const_container_of(cur_key, InfoVtype, key)->num;
            cswritep = strcpya(cswritep, cur_key);
            if (!IsInfoVtypeARSkip(knum)) {
              char* v_start = basic_info_fields[kidx];
              if (v_start && (!memequal_k(v_start, locked_missing_semicolon_str, 3))) {
                char* v_end = AdvPastDelim(v_start, ';');
                cswritep = memcpya(cswritep, v_start, v_end - v_start);
              } else if (knum <= 1) {
                cswritep = strcpya_k(cswritep, "=.;");
              } else {
                *cswritep++ = '=';
                cswritep = u16setsa(cswritep, 0x2c2e, knum);
                cswritep[-1] = ';';
              }
            } else {
              *cswritep++ = '=';
              char** cur_ar_info_fields = ar_info_fields[kidx];
              const uint32_t first_aidx = knum - kInfoVtypeR;
              for (uint32_t write_idx = 0; write_idx != merged_info_allele_ct - first_aidx; ++write_idx) {
                char* v_start = cur_ar_info_fields[write_idx];
                if (v_start && (!memequal_k(v_start, locked_missing_comma_str, 2))) {
                  char* v_end = AdvPastDelim(v_start, ',');
                  cswritep = memcpya(cswritep, v_start, v_end - v_start);
                } else {
                  cswritep = strcpya_k(cswritep, ".,");
                }
              }
              cswritep[-1] = ';';
            }
          }
        }
        --cswritep;
      }
    }
    if (pmcp->write_cm) {
      // Very similar to QUAL.  Differences:
      // - No 'min' option.
      // - Double-precision float parsing required.
      // - '0' is the missing value, not '.'.
      const MergeInfoCmMode merge_cm_mode = pmcp->merge_cm_mode;
      *cswritep++ = '\t';
      if ((merge_cm_mode == kMergeInfoCmModeFirst) || (merge_cm_mode == kMergeInfoCmModeNmFirst)) {
        uintptr_t rec_idx = 0;
        for (; rec_idx != merge_rec_ct; ++rec_idx) {
          SamePosPvarRecord* cur_record = same_id_records[rec_idx];
          const uint32_t* other_field_offsets = cur_record->other_field_offsets;
          char* tmp_variant_id = cur_record->variant_id;
          char* cm_str = &(tmp_variant_id[other_field_offsets[5]]);
          if (cm_str[0] != '\0') {
            if ((!strequal_k_unsafe(cm_str, "0")) || (merge_cm_mode == kMergeInfoCmModeFirst)) {
              cswritep = strcpya(cswritep, cm_str);
              break;
            }
          }
        }
        if (rec_idx == merge_rec_ct) {
          *cswritep++ = '0';
        }
      } else {
        char** cur_cm_strs = pmcp->cur_fields_buf;
        uintptr_t nz_cm_ct = 0;
        for (uintptr_t rec_idx = 0; rec_idx != merge_rec_ct; ++rec_idx) {
          SamePosPvarRecord* cur_record = same_id_records[rec_idx];
          const uint32_t* other_field_offsets = cur_record->other_field_offsets;
          char* tmp_variant_id = cur_record->variant_id;
          char* cm_str = &(tmp_variant_id[other_field_offsets[5]]);
          if ((cm_str[0] != '\0') && (!strequal_k_unsafe(cm_str, "0"))) {
            cur_cm_strs[nz_cm_ct] = cm_str;
            ++nz_cm_ct;
          }
        }
        if (!nz_cm_ct) {
          *cswritep++ = '0';
        } else if (nz_cm_ct == 1) {
          cswritep = strcpya(cswritep, cur_cm_strs[0]);
        } else {
          // Parse and check equality of CM values.
          double main_cm = 0.0;
          uintptr_t nm_rec_idx = 0;
          for (; nm_rec_idx != nz_cm_ct; ++nm_rec_idx) {
            double cur_cm;
            if (unlikely(!ScantokDouble(cur_cm_strs[nm_rec_idx], &cur_cm))) {
              logputs("\n");
              logerrprintfww("Error: Invalid CM value for variant '%s' at %s:%u.\n", cur_variant_id, pmcp->chr_buf, pmcp->cur_bp);
              return kPglRetMalformedInput;
            }
            if (!nm_rec_idx) {
              main_cm = cur_cm;
              continue;
            }
            if (cur_cm != main_cm) {
              break;
            }
          }
          if (nm_rec_idx == nz_cm_ct) {
            cswritep = dtoa_g(main_cm, cswritep);
          } else {
            *cswritep++ = '.';
          }
        }
      }
    }
  }
  AppendBinaryEoln(&cswritep);
  *cur_line_blenp = cswritep - pmcp->cswritep;
  if (unlikely(Cswrite(&pmcp->css, &cswritep))) {
    return kPglRetWriteFail;
  }
  pmcp->cswritep = cswritep;
  return kPglRetSuccess;
}

// We don't want to crash in an uncontrolled manner if the .pgen contains an
// instance of an allele that the .pvar claims shouldn't exist.
// If normalize is true, this also flips the variant when it's REF that's
// missing in the .pvar.
PglErr ValidateBiallelicVariantWithMissingCode(const AlleleCode* cur_allele_remap, uint32_t sample_ct, uint32_t normalize, PgenVariant* pgvp) {
  uintptr_t* genovec = pgvp->genovec;
  ZeroTrailingNyps(sample_ct, genovec);
  STD_ARRAY_DECL(uint32_t, 4, genocounts);
  GenoarrCountFreqsUnsafe(genovec, sample_ct, genocounts);
  if (unlikely(genocounts[1] || (genocounts[0] && (cur_allele_remap[0] == kMissingAlleleCode)) || (genocounts[2] && (cur_allele_remap[1] == kMissingAlleleCode)) || pgvp->dosage_ct)) {
    logputs("\n");
    logerrputs("Error: Missing allele in .pvar is present in the .pgen.\n");
    return kPglRetInconsistentInput;
  }
  if (normalize && ((cur_allele_remap[0] == 1) || (cur_allele_remap[1] == 0))) {
    // only possible when [0] == missing
    GenovecInvertUnsafe(sample_ct, genovec);
    if (pgvp->phasepresent_ct) {
      BitvecInvert(BitCtToWordCt(sample_ct), pgvp->phaseinfo);
    }
    if (pgvp->dosage_ct) {
      BiallelicDosage16Invert(pgvp->dosage_ct, pgvp->dosage_main);
      if (pgvp->dphase_ct) {
        BiallelicDphase16Invert(pgvp->dphase_ct, pgvp->dphase_delta);
      }
    }
  }
  return kPglRetSuccess;
}

void Update8bitDenseFromSparse(const uintptr_t* subset_mask, const void* __restrict src, uint32_t subset_size, void* __restrict dst) {
  const unsigned char* src_uc = S_CAST(const unsigned char*, src);
  unsigned char* dst_uc = S_CAST(unsigned char*, dst);
  uintptr_t uidx_base = 0;
  uintptr_t cur_bits = subset_mask[0];
  for (uint32_t uii = 0; uii != subset_size; ++uii) {
    const uintptr_t cur_uidx = BitIter1(subset_mask, &uidx_base, &cur_bits);
    dst_uc[cur_uidx] = src_uc[uii];
  }
}

void Update16bitDenseFromSparse(const uintptr_t* subset_mask, const void* __restrict src, uint32_t subset_size, void* __restrict dst) {
  const uint16_t* src_u16 = S_CAST(const uint16_t*, src);
  uint16_t* dst_u16 = S_CAST(uint16_t*, dst);
  uintptr_t uidx_base = 0;
  uintptr_t cur_bits = subset_mask[0];
  for (uint32_t uii = 0; uii != subset_size; ++uii) {
    const uintptr_t cur_uidx = BitIter1(subset_mask, &uidx_base, &cur_bits);
    dst_u16[cur_uidx] = src_u16[uii];
  }
}

void CopyAndPermuteBitarr(const uintptr_t* __restrict read_bitarr, const uint32_t* __restrict old_sample_idx_to_new, uint32_t write_sample_ctl, uint32_t read_ct, uintptr_t* __restrict write_bitarr) {
  ZeroWArr(write_sample_ctl, write_bitarr);
  uintptr_t read_uidx_base = 0;
  uintptr_t cur_bits = read_bitarr[0];
  for (uint32_t read_idx = 0; read_idx != read_ct; ++read_idx) {
    const uint32_t read_uidx = BitIter1(read_bitarr, &read_uidx_base, &cur_bits);
    const uint32_t write_sample_idx = old_sample_idx_to_new[read_uidx];
    SetBit(write_sample_idx, write_bitarr);
  }
}

// Subset of CopyAndPermute8bit().
void PermuteUpdate8bitDenseFromSparse(const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uint32_t* __restrict old_sample_idx_to_new, uint32_t sample_ct, uint32_t val_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  ZeroWArr(sample_ctl, dst_subset);
  const unsigned char* src_vals_uc = S_CAST(const unsigned char*, src_vals);
  unsigned char* dst_vals_uc = S_CAST(unsigned char*, dst_vals);

  uintptr_t cur_bits = src_subset[0];
  uintptr_t old_sample_idx_base = 0;
  for (uint32_t old_val_idx = 0; old_val_idx != val_ct; ++old_val_idx) {
    const uint32_t old_sample_idx = BitIter1(src_subset, &old_sample_idx_base, &cur_bits);
    const uint32_t new_sample_idx = old_sample_idx_to_new[old_sample_idx];
    SetBit(new_sample_idx, dst_subset);
    dst_vals_uc[new_sample_idx] = src_vals_uc[old_val_idx];
  }
}

// ok for dst_subset to be nullptr.
void PermuteUpdate16bitDenseFromSparse(const uintptr_t* __restrict src_subset, const void* __restrict src_vals, const uint32_t* __restrict old_sample_idx_to_new, uint32_t sample_ct, uint32_t val_ct, uintptr_t* __restrict dst_subset, void* __restrict dst_vals) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint16_t* src_vals_u16 = S_CAST(const uint16_t*, src_vals);
  uint16_t* dst_vals_u16 = S_CAST(uint16_t*, dst_vals);

  uintptr_t cur_bits = src_subset[0];
  uintptr_t old_sample_idx_base = 0;
  if (dst_subset) {
    // bugfix (15 May 2025): forgot to put this behind dst_subset != nullptr
    ZeroWArr(sample_ctl, dst_subset);
    for (uint32_t old_val_idx = 0; old_val_idx != val_ct; ++old_val_idx) {
      const uint32_t old_sample_idx = BitIter1(src_subset, &old_sample_idx_base, &cur_bits);
      const uint32_t new_sample_idx = old_sample_idx_to_new[old_sample_idx];
      SetBit(new_sample_idx, dst_subset);
      dst_vals_u16[new_sample_idx] = src_vals_u16[old_val_idx];
    }
  } else {
    for (uint32_t old_val_idx = 0; old_val_idx != val_ct; ++old_val_idx) {
      const uint32_t old_sample_idx = BitIter1(src_subset, &old_sample_idx_base, &cur_bits);
      const uint32_t new_sample_idx = old_sample_idx_to_new[old_sample_idx];
      dst_vals_u16[new_sample_idx] = src_vals_u16[old_val_idx];
    }
  }
}

void PermuteUpdateGenovec(const uintptr_t* __restrict r_genovec, const uint32_t* __restrict clobber_sample_idx_to_new, uint32_t clobber_sample_ct, uintptr_t* __restrict genovec) {
  // possible todo: add most-common-geno optimization when clobber_sample_ct is
  // sufficiently large.  That should be unlikely, though.
  const uint32_t read_genoword_ct_m1 = (clobber_sample_ct - 1) / kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    // We know the existing genovec bits at the positions we're clobbering are
    // all set, so xor operations make sense here.
    uintptr_t geno_word_xor;
    if (widx >= read_genoword_ct_m1) {
      if (widx > read_genoword_ct_m1) {
        break;
      }
      geno_word_xor = bzhi_max(~(r_genovec[widx]), 2 * ModNz(clobber_sample_ct, kBitsPerWordD2));
    } else {
      geno_word_xor = ~(r_genovec[widx]);
    }
    if (!geno_word_xor) {
      continue;
    }
    const uint32_t* cur_old_sample_idx_to_new = &(clobber_sample_idx_to_new[widx * kBitsPerWordD2]);
    // possible todo: branch between this and tighter loop, depending on number
    // of missing genotypes
    do {
      const uint32_t bit_read_shift_ct = ctzw(geno_word_xor) & (kBitsPerWord - 2);
      const uintptr_t cur_geno_xor = (geno_word_xor >> bit_read_shift_ct) & 3;
      const uint32_t new_sample_idx = cur_old_sample_idx_to_new[bit_read_shift_ct / 2];
      const uint32_t new_word_idx = new_sample_idx / kBitsPerWordD2;
      const uint32_t bit_write_shift_ct = 2 * (new_sample_idx % kBitsPerWordD2);
      genovec[new_word_idx] ^= cur_geno_xor << bit_write_shift_ct;
      geno_word_xor &= (~(3 * k1LU)) << bit_read_shift_ct;
    } while (geno_word_xor);
  }
}

void PermuteUpdateHphase(const uintptr_t* __restrict r_phasepresent, const uintptr_t* __restrict r_phaseinfo, const uint32_t* __restrict old_sample_idx_to_new, uint32_t phasepresent_ct, uintptr_t* __restrict phasepresent, uintptr_t* __restrict phaseinfo) {
  uintptr_t read_widx = 0;
  uintptr_t cur_bits = r_phasepresent[0];
  for (uint32_t phasepresent_idx = 0; phasepresent_idx != phasepresent_ct; ++phasepresent_idx) {
    const uintptr_t read_sample_idx_lowbit_idx = BitIter1x(r_phasepresent, &read_widx, &cur_bits);
    const uintptr_t is_phaseinfo = (r_phaseinfo[read_widx] >> read_sample_idx_lowbit_idx) & 1;
    const uint32_t read_sample_idx = read_widx * kBitsPerWord + read_sample_idx_lowbit_idx;
    const uint32_t write_sample_idx = old_sample_idx_to_new[read_sample_idx];
    const uint32_t write_widx = write_sample_idx / kBitsPerWord;
    const uintptr_t write_bit = k1LU << (write_sample_idx % kBitsPerWord);
    phasepresent[write_widx] |= write_bit;
    phaseinfo[write_widx] |= write_bit & (-is_phaseinfo);
  }
}

void CopyAndPermuteHphase(const uintptr_t* __restrict read_phasepresent, const uintptr_t* __restrict read_phaseinfo, const uint32_t* __restrict old_sample_idx_to_new, uint32_t write_sample_ctl, uint32_t phasepresent_ct, uintptr_t* __restrict write_phasepresent, uintptr_t* __restrict write_phaseinfo) {
  ZeroWArr(write_sample_ctl, write_phasepresent);
  ZeroWArr(write_sample_ctl, write_phaseinfo);
  PermuteUpdateHphase(read_phasepresent, read_phaseinfo, old_sample_idx_to_new, phasepresent_ct, write_phasepresent, write_phaseinfo);
}

void Compare8bitDense(const uintptr_t* __restrict present1, const void* __restrict dense1, const uintptr_t* __restrict present2, const void* __restrict dense2, uint32_t word_ct, uintptr_t* __restrict compare_mask) {
  const unsigned char* dense1_uc = S_CAST(const unsigned char*, dense1);
  const unsigned char* dense2_uc = S_CAST(const unsigned char*, dense2);
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    uintptr_t compare_word = compare_mask[widx];
    if (!compare_word) {
      continue;
    }
    const uintptr_t present1_word = present1[widx];
    const uintptr_t present2_word = present2[widx];
    uintptr_t byte_compare_word = compare_word & present1_word & present2_word;
    compare_word &= ~(present1_word ^ present2_word);
    if (byte_compare_word) {
      // possible todo: when byte_compare_word has enough set bits, use
      // movemask instead
      const unsigned char* cur_dense1_uc = &(dense1_uc[widx * kBitsPerWord]);
      const unsigned char* cur_dense2_uc = &(dense2_uc[widx * kBitsPerWord]);
      do {
        const uint32_t sample_idx_lowbits = ctzw(byte_compare_word);
        if (cur_dense1_uc[sample_idx_lowbits] != cur_dense2_uc[sample_idx_lowbits]) {
          // (byte_compare_word & (-byte_compare_word)) also isolates low bit
          compare_word -= k1LU << sample_idx_lowbits;
        }
        byte_compare_word &= byte_compare_word - 1;
      } while (byte_compare_word);
    }
    compare_mask[widx] = compare_word;
  }
}

void Compare16bitDense(const uintptr_t* __restrict present1, const void* __restrict dense1, const uintptr_t* __restrict present2, const void* __restrict dense2, uint32_t word_ct, uintptr_t* __restrict compare_mask) {
  const uint16_t* dense1_u16 = S_CAST(const uint16_t*, dense1);
  const uint16_t* dense2_u16 = S_CAST(const uint16_t*, dense2);
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    uintptr_t compare_word = compare_mask[widx];
    if (!compare_word) {
      continue;
    }
    const uintptr_t present1_word = present1[widx];
    const uintptr_t present2_word = present2[widx];
    uintptr_t byte_compare_word = compare_word & present1_word & present2_word;
    compare_word &= ~(present1_word ^ present2_word);
    if (byte_compare_word) {
      const uint16_t* cur_dense1_u16 = &(dense1_u16[widx * kBitsPerWord]);
      const uint16_t* cur_dense2_u16 = &(dense2_u16[widx * kBitsPerWord]);
      do {
        const uint32_t sample_idx_lowbits = ctzw(byte_compare_word);
        if (cur_dense1_u16[sample_idx_lowbits] != cur_dense2_u16[sample_idx_lowbits]) {
          compare_word -= k1LU << sample_idx_lowbits;
        }
        byte_compare_word &= byte_compare_word - 1;
      } while (byte_compare_word);
    }
    compare_mask[widx] = compare_word;
  }
}

void Update8bitDense(const uintptr_t* __restrict update_mask, const uintptr_t* __restrict src_present, const void* __restrict src_dense, uint32_t word_ct, uintptr_t* __restrict dst_present, void* __restrict dst_dense) {
  const unsigned char* src_dense_uc = S_CAST(const unsigned char*, src_dense);
  unsigned char* dst_dense_uc = S_CAST(unsigned char*, dst_dense);
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t update_word = update_mask[widx];
    if (!update_word) {
      continue;
    }
    uintptr_t src_present_word = src_present[widx] & update_word;
    if (!src_present_word) {
      continue;
    }
    dst_present[widx] |= src_present_word;
    const unsigned char* cur_src_dense = &(src_dense_uc[widx * kBitsPerWord]);
    unsigned char* cur_dst_dense = &(dst_dense_uc[widx * kBitsPerWord]);
    do {
      const uint32_t sample_idx_lowbits = ctzw(src_present_word);
      cur_dst_dense[sample_idx_lowbits] = cur_src_dense[sample_idx_lowbits];
      src_present_word &= src_present_word - 1;
    } while (src_present_word);
  }
}

void Update16bitDense(const uintptr_t* __restrict update_mask, const uintptr_t* __restrict src_present, const void* __restrict src_dense, uint32_t word_ct, uintptr_t* __restrict dst_present, void* __restrict dst_dense) {
  const uint16_t* src_dense_u16 = S_CAST(const uint16_t*, src_dense);
  uint16_t* dst_dense_u16 = S_CAST(uint16_t*, dst_dense);
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t update_word = update_mask[widx];
    if (!update_word) {
      continue;
    }
    uintptr_t src_present_word = src_present[widx] & update_word;
    if (!src_present_word) {
      continue;
    }
    dst_present[widx] |= src_present_word;
    const uint16_t* cur_src_dense = &(src_dense_u16[widx * kBitsPerWord]);
    uint16_t* cur_dst_dense = &(dst_dense_u16[widx * kBitsPerWord]);
    do {
      const uint32_t sample_idx_lowbits = ctzw(src_present_word);
      cur_dst_dense[sample_idx_lowbits] = cur_src_dense[sample_idx_lowbits];
      src_present_word &= src_present_word - 1;
    } while (src_present_word);
  }
}

typedef struct MergeReaderStruct {
  PgenReader pgr;
  uintptr_t* sample_include;
  // yes, this breaks the standard "local variable only" usage rule; we don't
  // want to clear the LD-cache on every single variant.  Not dangerous here
  // since we never change the sample-subset.
  PgrSampleSubsetIndex pssi;

  uintptr_t* sample_span;
  uint32_t* old_sample_idx_to_new;
  uint32_t sample_idx_increasing;  // =2 in simplest case

  uint32_t sample_ct;
} MergeReader;

typedef struct MergeWriterStruct {
  STPgenWriter spgw;
  // Main write buffers.
  // probable todo: define PgenVariant-based write functions, and replace these
  // fields with a third PgenVariant buffer.
  // todo: document which buffers must be allocated up to next vector boundary.
  // Not urgent since we're allocating up to that point by default, but this
  // would make it easier to patch in debug-allocator logic to test for
  // out-of-bounds memory writes.
  uintptr_t* genovec;
  uintptr_t* patch_01_set;
  AlleleCode* patch_01_vals;
  uintptr_t* patch_10_set;
  AlleleCode* patch_10_vals;
  uintptr_t* phasepresent;
  uintptr_t* phaseinfo;
  uintptr_t* dosage_present;
  uint16_t* dosage_main;
  uintptr_t* dphase_present;
  int16_t* dphase_delta;

  PgenVariant pgv_readbuf;

  uintptr_t* unlocked_set;
  // "--merge-mode nm-match" only
  uintptr_t* unlocked_missing_set;

  // "--merge-mode nm-match"-specific temporary buffers
  uintptr_t* clobber_sample_span;
  uintptr_t* unlocked_nonmissing_sample_span;

  // Buffers supporting efficient clobber when only some of the currently-read
  // samples qualify.
  uintptr_t* unlocked_sample_span;
  uint32_t* clobber_sample_idx_to_new;
  uintptr_t* mask_buf;
  PgenVariant pgv_midbuf;
  uintptr_t* phaseinfo_xor;
  AlleleCode* wide_codes;

  MergeMode merge_mode;
} MergeWriter;

PglErr MergePgenVariantNoTmpLocked(SamePosPvarRecord** same_id_records, const AlleleCode* master_allele_remap, uintptr_t merge_rec_ct, uint32_t write_allele_ct, uint32_t allele_remap_stride, MergeReader** mrp_arr, MergeWriter* mwp) {
  PglErr reterr = kPglRetSuccess;
  {
    STPgenWriter* spgwp = &(mwp->spgw);
    PgenVariant* pgvp = &(mwp->pgv_readbuf);
    const uint32_t write_sample_ct = SpgwGetSampleCt(spgwp);
    const uint32_t write_sample_ctl = BitCtToWordCt(write_sample_ct);

    // Only used in multi-file case.
    uint32_t hphase_exists = 0;
    uint32_t dosage_exists = 0;
    uint32_t dphase_exists = 0;
    // "simple" = either non-allele-permuting, or {[0] = kMissingAlleleCode,
    // [1] = 0}.
    // (These are the only possibilities when merge_rec_ct == 1.)
    uint32_t simple_first_allele_remap = 1;
    if (merge_rec_ct > 1) {
      const uint32_t read_allele_ct = same_id_records[0]->allele_ct;
      for (uint32_t allele_idx = 0; allele_idx != read_allele_ct; ++allele_idx) {
        if (master_allele_remap[allele_idx] != allele_idx) {
          if (read_allele_ct == 2) {
            // Check for permitted kMissingAlleleCode cases.
            if (!allele_idx) {
              if (((master_allele_remap[0] == kMissingAlleleCode) && ((master_allele_remap[1] <= 1) || (master_allele_remap[1] == kMissingAlleleCode))) || ((master_allele_remap[0] == 1) && (master_allele_remap[1] == kMissingAlleleCode))) {
                break;
              }
            } else {
              if (master_allele_remap[1] == kMissingAlleleCode) {
                break;
              }
            }
          }
          simple_first_allele_remap = 0;
          break;
        }
      }
      uint32_t vrtype_or = 0;
      for (uintptr_t rec_idx = 0; rec_idx != merge_rec_ct; ++rec_idx) {
        const uint32_t file_idx = same_id_records[rec_idx]->secondary_key >> 32;
        const uint32_t read_variant_uidx = S_CAST(uint32_t, same_id_records[0]->secondary_key);
        MergeReader* cur_mrp = mrp_arr[file_idx];
        PgenReader* pgrp = &(cur_mrp->pgr);
        vrtype_or |= PgrGetVrtype(pgrp, read_variant_uidx);
      }
      hphase_exists = (vrtype_or / 0x10) & 1;
      dosage_exists = !!(vrtype_or & 0x60);
      dphase_exists = (vrtype_or / 0x80) & 1;
      assert((!dphase_exists) || dosage_exists);
      if (dosage_exists && (write_allele_ct > 2)) {
        logerrputs("Error: --pmerge[-list] multiallelic-variant dosage support is under development.\n");
        reterr = kPglRetNotYetSupported;
        goto MergePgenVariantNoTmpLocked_ret_1;
      }
    }
    const MergeMode merge_mode = mwp->merge_mode;
    if (simple_first_allele_remap) {
      const uint32_t file_idx = same_id_records[0]->secondary_key >> 32;
      const uint32_t read_variant_uidx = S_CAST(uint32_t, same_id_records[0]->secondary_key);
      MergeReader* cur_mrp = mrp_arr[file_idx];
      const uint32_t read_sample_ct = cur_mrp->sample_ct;
      const uint32_t read_allele_ct = same_id_records[0]->allele_ct;
      PgenReader* pgrp = &(cur_mrp->pgr);
      const uint32_t vrtype = PgrGetVrtype(pgrp, read_variant_uidx);
      const uint32_t read_phase_present = !!(vrtype & 0x90);
      const uint32_t read_dosage_present = !!(vrtype & 0x60);
      const uintptr_t* sample_include = cur_mrp->sample_include;
      if ((read_allele_ct == 2) && (!read_dosage_present)) {
        pgvp->patch_01_ct = 0;
        pgvp->patch_10_ct = 0;
        pgvp->dosage_ct = 0;
        if (!read_phase_present) {
          pgvp->phasepresent_ct = 0;
          reterr = PgrGet(sample_include, cur_mrp->pssi, read_sample_ct, read_variant_uidx, pgrp, pgvp->genovec);
        } else {
          reterr = PgrGetP(sample_include, cur_mrp->pssi, read_sample_ct, read_variant_uidx, pgrp, pgvp->genovec, pgvp->phasepresent, pgvp->phaseinfo, &pgvp->phasepresent_ct);
        }
      } else {
        reterr = PgrGetMDp(sample_include, cur_mrp->pssi, read_sample_ct, read_variant_uidx, pgrp, pgvp);
      }
      if (unlikely(reterr)) {
        PgenErrPrintNV(reterr, read_variant_uidx);
        goto MergePgenVariantNoTmpLocked_ret_1;
      }
      if ((master_allele_remap[0] == kMissingAlleleCode) || (master_allele_remap[1] == kMissingAlleleCode)) {
        reterr = ValidateBiallelicVariantWithMissingCode(master_allele_remap, read_sample_ct, 1, pgvp);
        if (unlikely(reterr)) {
          goto MergePgenVariantNoTmpLocked_ret_1;
        }
      }
      uintptr_t* read_genovec = pgvp->genovec;
      ZeroTrailingNyps(read_sample_ct, read_genovec);
      const uint32_t sample_idx_increasing = cur_mrp->sample_idx_increasing;
      const uint32_t write_biallelic = (!pgvp->patch_01_ct) && (!pgvp->patch_10_ct);
      uint32_t unlocked_ct = 0;
      if (sample_idx_increasing == 2) {
        // read_sample_ct == write_sample_ct
        if ((merge_rec_ct != 1) && (merge_mode != kMergeModeFirst)) {
          // If there are any missing genotypes/dosages, those are still
          // unlocked.
          uintptr_t* missingness = (merge_mode == kMergeModeNmFirst)? mwp->unlocked_set : mwp->unlocked_missing_set;
          GenoarrToMissingnessUnsafe(read_genovec, write_sample_ct, missingness);
          if (pgvp->dosage_ct) {
            BitvecInvmask(pgvp->dosage_present, write_sample_ctl, missingness);
          }
          if (merge_mode == kMergeModeNmFirst) {
            unlocked_ct = PopcountWords(missingness, write_sample_ctl);
          } else {
            // kMergeModeNmMatch
            SetAllBits(write_sample_ct, mwp->unlocked_set);
            unlocked_ct = write_sample_ct;
          }
        }
        if (unlocked_ct == 0) {
          if (write_biallelic) {
            // Simplest cases: write directly from pgvp, no need to use
            // mwp.genovec, etc.
            // We explicitly branch on phasepresent_ct since otherwise we'd
            // need to pass phasepresent == nullptr.
            if (pgvp->dosage_ct == 0) {
              if (!pgvp->phasepresent_ct) {
                reterr = SpgwAppendBiallelicGenovec(read_genovec, spgwp);
              } else {
                reterr = SpgwAppendBiallelicGenovecHphase(read_genovec, pgvp->phasepresent, pgvp->phaseinfo, spgwp);
              }
            } else {
              if ((!pgvp->phasepresent_ct) && (!pgvp->dphase_ct)) {
                reterr = SpgwAppendBiallelicGenovecDosage16(read_genovec, pgvp->dosage_present, pgvp->dosage_main, pgvp->dosage_ct, spgwp);
              } else {
                if (!pgvp->phasepresent_ct) {
                  ZeroWArr(write_sample_ctl, pgvp->phasepresent);
                }
                reterr = SpgwAppendBiallelicGenovecDphase16(read_genovec, pgvp->phasepresent, pgvp->phaseinfo, pgvp->dosage_present, pgvp->dphase_present, pgvp->dosage_main, pgvp->dphase_delta, pgvp->dosage_ct, pgvp->dphase_ct, spgwp);
              }
            }
          } else {
            // Multiallelic.
            ZeroTrailingNyps(write_sample_ct, read_genovec);
            assert(!pgvp->dosage_ct);  // not yet supported
            if (!pgvp->phasepresent_ct) {
              reterr = SpgwAppendMultiallelicSparse(read_genovec, pgvp->patch_01_set, pgvp->patch_01_vals, pgvp->patch_10_set, pgvp->patch_10_vals, write_allele_ct, pgvp->patch_01_ct, pgvp->patch_10_ct, spgwp);
            } else {
              reterr = SpgwAppendMultiallelicGenovecHphase(read_genovec, pgvp->patch_01_set, pgvp->patch_01_vals, pgvp->patch_10_set, pgvp->patch_10_vals, pgvp->phasepresent, pgvp->phaseinfo, write_allele_ct, pgvp->patch_01_ct, pgvp->patch_10_ct, spgwp);
            }
          }
          goto MergePgenVariantNoTmpLocked_ret_1;
        }
        // Copy to write buffers.
        // bugfix (2 Mar 2021): did not copy enough bytes here
        const uint32_t write_sample_ctl2 = NypCtToWordCt(write_sample_ct);
        memcpy(mwp->genovec, read_genovec, write_sample_ctl2 * sizeof(intptr_t));
        const uint32_t write_sample_ctb = write_sample_ctl * sizeof(intptr_t);
        if (write_allele_ct > 2) {
          if (pgvp->patch_01_ct) {
            memcpy(mwp->patch_01_set, pgvp->patch_01_set, write_sample_ctb);
            Update8bitDenseFromSparse(pgvp->patch_01_set, pgvp->patch_01_vals, pgvp->patch_01_ct, mwp->patch_01_vals);
          } else {
            ZeroWArr(write_sample_ctl, mwp->patch_01_set);
          }
          if (pgvp->patch_10_ct) {
            memcpy(mwp->patch_10_set, pgvp->patch_10_set, write_sample_ctb);
            Update16bitDenseFromSparse(pgvp->patch_10_set, pgvp->patch_10_vals, pgvp->patch_10_ct, mwp->patch_10_vals);
          } else {
            ZeroWArr(write_sample_ctl, mwp->patch_10_set);
          }
        }
        if (hphase_exists || dphase_exists) {
          if (pgvp->phasepresent_ct) {
            memcpy(mwp->phasepresent, pgvp->phasepresent, write_sample_ctb);
            memcpy(mwp->phaseinfo, pgvp->phaseinfo, write_sample_ctb);
          } else {
            ZeroWArr(write_sample_ctl, mwp->phasepresent);
          }
        }
        if (dosage_exists) {
          if (pgvp->dosage_ct) {
            memcpy(mwp->dosage_present, pgvp->dosage_present, write_sample_ctb);
            Update16bitDenseFromSparse(pgvp->dosage_present, pgvp->dosage_main, pgvp->dosage_ct, mwp->dosage_main);
          } else {
            ZeroWArr(write_sample_ctl, mwp->dosage_present);
          }
          if (dphase_exists) {
            if (pgvp->dphase_ct) {
              memcpy(mwp->dphase_present, pgvp->dphase_present, write_sample_ctb);
              Update16bitDenseFromSparse(pgvp->dphase_present, pgvp->dphase_delta, pgvp->dosage_ct, mwp->dphase_delta);
            } else {
              ZeroWArr(write_sample_ctl, mwp->dphase_present);
            }
          }
        }
      } else {
        // Necessary to expand and/or permute.
        // We'd rather not iterate over every single genotype.  So we exploit
        // sparsity, efficiently prefilling all relevant genotypes in the
        // output array with the most common value, and then patching only the
        // exceptions.

        // (strictly speaking, we don't need to recompute genocounts if we
        // already did so in the kMissingAlleleCode cases above, but that event
        // should be rare enough that it's not worth deduplicating.)
        const uintptr_t most_common_geno_word = MostCommonGenoUnsafe(read_genovec, read_sample_ct) * kMask5555;

        const uint32_t* old_sample_idx_to_new = cur_mrp->old_sample_idx_to_new;
        const uintptr_t* sample_span = cur_mrp->sample_span;
        // This part is very similar to GenovecPermute().  (We could just call
        // it in the read_sample_ct == write_sample_ct case.)
        uintptr_t* genovec = mwp->genovec;
        const uint32_t write_sample_ctl2 = NypCtToWordCt(write_sample_ct);
        if (read_sample_ct == write_sample_ct) {
          for (uint32_t hwidx = 0; hwidx != write_sample_ctl2; ++hwidx) {
            genovec[hwidx] = most_common_geno_word;
          }
        } else {
          const Halfword* sample_span_hw = DowncastKWToHW(sample_span);
          for (uint32_t hwidx = 0; hwidx != write_sample_ctl2; ++hwidx) {
            // Initialize each genotype outside the single source's sample-span
            // to missing (both bits set).
            const Halfword inv_sample_span_hw = ~sample_span_hw[hwidx];
            genovec[hwidx] = most_common_geno_word | (UnpackHalfwordToWord(inv_sample_span_hw) * 3);
          }
        }
        ZeroTrailingNyps(write_sample_ct, genovec);
        const uint32_t read_genoword_ct_m1 = (read_sample_ct - 1) / kBitsPerWordD2;
        for (uint32_t widx = 0; ; ++widx) {
          uintptr_t geno_word_xor;
          if (widx >= read_genoword_ct_m1) {
            if (widx > read_genoword_ct_m1) {
              break;
            }
            geno_word_xor = bzhi_max(most_common_geno_word ^ read_genovec[widx], 2 * ModNz(read_sample_ct, kBitsPerWordD2));
          } else {
            geno_word_xor = most_common_geno_word ^ read_genovec[widx];
          }
          if (!geno_word_xor) {
            continue;
          }
          const uint32_t* cur_old_sample_idx_to_new = &(old_sample_idx_to_new[widx * kBitsPerWordD2]);
          do {
            const uint32_t bit_read_shift_ct = ctzw(geno_word_xor) & (kBitsPerWord - 2);
            const uintptr_t cur_geno_xor = (geno_word_xor >> bit_read_shift_ct) & 3;
            const uint32_t new_sample_idx = cur_old_sample_idx_to_new[bit_read_shift_ct / 2];
            const uint32_t new_word_idx = new_sample_idx / kBitsPerWordD2;
            const uint32_t bit_write_shift_ct = 2 * (new_sample_idx % kBitsPerWordD2);
            // Value has been preset to most_common_geno, so if we xor it with
            // (most_common_geno ^ actual_geno), the result is actual_geno.
            genovec[new_word_idx] ^= cur_geno_xor << bit_write_shift_ct;
            geno_word_xor &= (~(3 * k1LU)) << bit_read_shift_ct;
          } while (geno_word_xor);
        }
        // Resolve dosage_present first, so we know unlocked_ct status and can
        // select correct (sparse vs. dense) conversion for patch_01_vals,
        // patch_10_vals, dosage_main, and dphase_delta.
        const uint32_t dosage_ct = pgvp->dosage_ct;
        if (dosage_ct) {
          assert(write_biallelic);
          if (sample_idx_increasing) {
            ExpandBytearr(pgvp->dosage_present, sample_span, write_sample_ctl, read_sample_ct, 0, mwp->dosage_present);
          } else {
            CopyAndPermuteBitarr(pgvp->dosage_present, old_sample_idx_to_new, write_sample_ctl, dosage_ct, mwp->dosage_present);
          }
        }
        if (merge_rec_ct != 1) {
          if (merge_mode == kMergeModeFirst) {
            if (read_sample_ct != write_sample_ct) {
              unlocked_ct = write_sample_ct - read_sample_ct;
              BitvecInvertCopy(sample_span, write_sample_ctl, mwp->unlocked_set);
              ZeroTrailingBits(write_sample_ct, mwp->unlocked_set);
            }
          } else {
            // If there are any missing genotypes/dosages, those are still
            // unlocked.
            uintptr_t* missingness = (merge_mode == kMergeModeNmFirst)? mwp->unlocked_set : mwp->unlocked_missing_set;
            GenoarrToMissingnessUnsafe(mwp->genovec, write_sample_ct, missingness);
            if (dosage_ct) {
              BitvecInvmask(mwp->dosage_present, write_sample_ctl, missingness);
            }
            if (merge_mode == kMergeModeNmFirst) {
              unlocked_ct = PopcountWords(missingness, write_sample_ctl);
            } else {
              // kMergeModeNmMatch
              SetAllBits(write_sample_ct, mwp->unlocked_set);
              unlocked_ct = write_sample_ct;
            }
          }
        }
        uint32_t patch_01_ct = 0;
        uint32_t patch_10_ct = 0;
        if (!write_biallelic) {
          patch_01_ct = pgvp->patch_01_ct;
          if (patch_01_ct) {
            if (sample_idx_increasing) {
              ExpandBytearr(pgvp->patch_01_set, sample_span, write_sample_ctl, read_sample_ct, 0, mwp->patch_01_set);
              if (!unlocked_ct) {
                memcpy(mwp->patch_01_vals, pgvp->patch_01_vals, patch_01_ct);
              } else {
                Update8bitDenseFromSparse(mwp->patch_01_set, pgvp->patch_01_vals, pgvp->patch_01_ct, mwp->patch_01_vals);
              }
            } else {
              if (!unlocked_ct) {
                CopyAndPermute8bit(nullptr, pgvp->patch_01_set, pgvp->patch_01_vals, old_sample_idx_to_new, read_sample_ct, patch_01_ct, mwp->patch_01_set, mwp->patch_01_vals);
              } else {
                PermuteUpdate8bitDenseFromSparse(pgvp->patch_01_set, pgvp->patch_01_vals, old_sample_idx_to_new, read_sample_ct, patch_01_ct, mwp->patch_01_set, mwp->patch_01_vals);
              }
            }
          }
          patch_10_ct = pgvp->patch_10_ct;
          if (patch_10_ct) {
            if (sample_idx_increasing) {
              ExpandBytearr(pgvp->patch_10_set, sample_span, write_sample_ctl, read_sample_ct, 0, mwp->patch_10_set);
              if (!unlocked_ct) {
                memcpy(mwp->patch_10_vals, pgvp->patch_10_vals, patch_10_ct * 2);
              } else {
                Update16bitDenseFromSparse(mwp->patch_10_set, pgvp->patch_10_vals, patch_10_ct, mwp->patch_10_vals);
              }
            } else {
              if (!unlocked_ct) {
                CopyAndPermute16bit(nullptr, pgvp->patch_10_set, pgvp->patch_10_vals, old_sample_idx_to_new, read_sample_ct, patch_10_ct, mwp->patch_10_set, mwp->patch_10_vals);
              } else {
                PermuteUpdate16bitDenseFromSparse(pgvp->patch_10_set, pgvp->patch_10_vals, old_sample_idx_to_new, read_sample_ct, patch_10_ct, mwp->patch_10_set, mwp->patch_10_vals);
              }
            }
          }
        }
        const uint32_t phasepresent_ct = pgvp->phasepresent_ct;
        if (phasepresent_ct) {
          if (sample_idx_increasing) {
            ExpandBytearr(pgvp->phasepresent, sample_span, write_sample_ctl, read_sample_ct, 0, mwp->phasepresent);
            ExpandBytearr(pgvp->phaseinfo, sample_span, write_sample_ctl, read_sample_ct, 0, mwp->phaseinfo);
          } else {
            CopyAndPermuteHphase(pgvp->phasepresent, pgvp->phaseinfo, old_sample_idx_to_new, write_sample_ctl, phasepresent_ct, mwp->phasepresent, mwp->phaseinfo);
          }
        }
        uint32_t dphase_ct = 0;
        if (dosage_ct) {
          assert(write_biallelic);
          // mwp->dosage_present already correct
          if (sample_idx_increasing) {
            if (!unlocked_ct) {
              memcpy(mwp->dosage_main, pgvp->dosage_main, dosage_ct * 2);
            } else {
              Update16bitDenseFromSparse(mwp->dosage_present, pgvp->dosage_main, dosage_ct, mwp->dosage_main);
            }
          } else {
            if (!unlocked_ct) {
              CopyAndPermute16bit(nullptr, pgvp->dosage_present, pgvp->dosage_main, old_sample_idx_to_new, read_sample_ct, dosage_ct, mwp->dosage_present, mwp->dosage_main);
            } else {
              PermuteUpdate16bitDenseFromSparse(pgvp->dosage_present, pgvp->dosage_main, old_sample_idx_to_new, read_sample_ct, dosage_ct, nullptr, mwp->dosage_main);
            }
          }
          dphase_ct = pgvp->dphase_ct;
          if (dphase_ct) {
            if (sample_idx_increasing) {
              ExpandBytearr(pgvp->dphase_present, sample_span, write_sample_ctl, read_sample_ct, 0, mwp->dphase_present);
              if (!unlocked_ct) {
                memcpy(mwp->dphase_delta, pgvp->dphase_delta, dphase_ct * 2);
              } else {
                Update16bitDenseFromSparse(mwp->dphase_present, pgvp->dphase_delta, dosage_ct, mwp->dphase_delta);
              }
            } else {
              if (!unlocked_ct) {
                CopyAndPermute16bit(nullptr, pgvp->dphase_present, pgvp->dphase_delta, old_sample_idx_to_new, read_sample_ct, dphase_ct, mwp->dphase_present, mwp->dphase_delta);
              } else {
                PermuteUpdate16bitDenseFromSparse(pgvp->dphase_present, pgvp->dphase_delta, old_sample_idx_to_new, read_sample_ct, dphase_ct, mwp->dphase_present, mwp->dphase_delta);
              }
            }
          }
        }
        if (unlocked_ct == 0) {
          if (!dosage_ct) {
            if (write_biallelic) {
              if (!phasepresent_ct) {
                reterr = SpgwAppendBiallelicGenovec(genovec, spgwp);
              } else {
                reterr = SpgwAppendBiallelicGenovecHphase(genovec, mwp->phasepresent, mwp->phaseinfo, spgwp);
              }
            } else {
              if (!phasepresent_ct) {
                reterr = SpgwAppendMultiallelicSparse(genovec, mwp->patch_01_set, mwp->patch_01_vals, mwp->patch_10_set, mwp->patch_10_vals, write_allele_ct, patch_01_ct, patch_10_ct, spgwp);
              } else {
                reterr = SpgwAppendMultiallelicGenovecHphase(genovec, mwp->patch_01_set, mwp->patch_01_vals, mwp->patch_10_set, mwp->patch_10_vals, mwp->phasepresent, mwp->phaseinfo, write_allele_ct, patch_01_ct, patch_10_ct, spgwp);
              }
            }
          } else {
            if ((!phasepresent_ct) && (!dphase_ct)) {
              reterr = SpgwAppendBiallelicGenovecDosage16(genovec, mwp->dosage_present, mwp->dosage_main, dosage_ct, spgwp);
            } else {
              if (!phasepresent_ct) {
                ZeroWArr(write_sample_ctl, mwp->phasepresent);
              }
              reterr = SpgwAppendBiallelicGenovecDphase16(genovec, mwp->phasepresent, mwp->phaseinfo, mwp->dosage_present, mwp->dphase_present, mwp->dosage_main, mwp->dphase_delta, dosage_ct, dphase_ct, spgwp);
            }
          }
          goto MergePgenVariantNoTmpLocked_ret_1;
        }
        // unlocked_ct > 0, some write buffers already copied to.  Only need
        // to initialize write buffers which don't happen to have been touched
        // by processing of the first file.
        if (write_allele_ct > 2) {
          if (!patch_01_ct) {
            ZeroWArr(write_sample_ctl, mwp->patch_01_set);
          }
          if (!patch_10_ct) {
            ZeroWArr(write_sample_ctl, mwp->patch_10_set);
          }
        }
        if ((hphase_exists || dphase_exists) && (!phasepresent_ct)) {
          ZeroWArr(write_sample_ctl, mwp->phasepresent);
        }
        if (dosage_exists && (!dosage_ct)) {
          ZeroWArr(write_sample_ctl, mwp->dosage_present);
          if (dphase_exists && (!dphase_ct)) {
            ZeroWArr(write_sample_ctl, mwp->dphase_present);
          }
        }
      }
    } else {
      SetAllBits(write_sample_ct * 2, mwp->genovec);
      SetAllBits(write_sample_ct, mwp->unlocked_set);
      if (mwp->unlocked_missing_set) {
        SetAllBits(write_sample_ct, mwp->unlocked_missing_set);
      }
      if (write_allele_ct > 2) {
        ZeroWArr(write_sample_ctl, mwp->patch_01_set);
        ZeroWArr(write_sample_ctl, mwp->patch_10_set);
      }
      if (hphase_exists || dphase_exists) {
        ZeroWArr(write_sample_ctl, mwp->phasepresent);
      }
      if (dosage_exists) {
        ZeroWArr(write_sample_ctl, mwp->dosage_present);
        if (dphase_exists) {
          ZeroWArr(write_sample_ctl, mwp->dphase_present);
        }
      }
    }

    // Multiple variants to merge.
    // 1. Missing-initialize write buffer and unlocked_set.  Dense
    //    representation of patch_01_vals, patch_10_vals, dosage_vals, and
    //    dphase_deltas.
    // 2. Incrementally update write buffer.
    // 3. Pack patch_01_vals, patch_10_vals, dosage_vals, and dphase_deltas,
    //    and flush.
    uintptr_t* genovec = mwp->genovec;
    uintptr_t* unlocked_set = mwp->unlocked_set;
    uintptr_t* unlocked_missing_set = mwp->unlocked_missing_set;
    uintptr_t* patch_01_set = nullptr;
    AlleleCode* patch_01_dense = nullptr;
    uintptr_t* patch_10_set = nullptr;
    AlleleCode* patch_10_dense = nullptr;
    if (write_allele_ct > 2) {
      patch_01_set = mwp->patch_01_set;
      patch_01_dense = mwp->patch_01_vals;
      patch_10_set = mwp->patch_10_set;
      patch_10_dense = mwp->patch_10_vals;
    }
    uintptr_t* phasepresent = nullptr;
    uintptr_t* phaseinfo = nullptr;
    if (hphase_exists || dphase_exists) {
      phasepresent = mwp->phasepresent;
      phaseinfo = mwp->phaseinfo;
    }
    uintptr_t* dosage_present = nullptr;
    uint16_t* dosage_dense = nullptr;
    uintptr_t* dphase_present = nullptr;
    int16_t* dphase_dense = nullptr;
    if (dosage_exists) {
      dosage_present = mwp->dosage_present;
      dosage_dense = mwp->dosage_main;
      if (dphase_exists) {
        dphase_present = mwp->dphase_present;
        dphase_dense = mwp->dphase_delta;
      }
    }
    const uint32_t write_sample_ctl2 = NypCtToWordCt(write_sample_ct);
    for (uintptr_t rec_idx = simple_first_allele_remap; rec_idx != merge_rec_ct; ++rec_idx) {
      const uint32_t file_idx = same_id_records[rec_idx]->secondary_key >> 32;
      MergeReader* cur_mrp = mrp_arr[file_idx];
      uintptr_t* sample_span = cur_mrp->sample_span;
      // Note that there should be a fast path for read_sample_ct == 1 in the
      // generic version of this function.  But probably not worth the effort
      // here; we'll only get here on a large merge of single-sample filesets
      // when the user has specified --merge-mode nm-first.
      const uint32_t read_sample_ct = cur_mrp->sample_ct;
      uintptr_t* unlocked_sample_span = mwp->unlocked_sample_span;
      BitvecAndCopy(sample_span, unlocked_set, write_sample_ctl, unlocked_sample_span);
      uint32_t unlocked_sample_ct = PopcountWords(unlocked_sample_span, write_sample_ctl);
      if (!unlocked_sample_ct) {
        // No point in loading this entry, all genotypes for these samples are
        // already locked.
        continue;
      }
      const uintptr_t* sample_include = cur_mrp->sample_include;
      const uint32_t read_variant_uidx = S_CAST(uint32_t, same_id_records[rec_idx]->secondary_key);
      // not const since there's an allele-rotation case where we change it
      uint32_t read_allele_ct = same_id_records[rec_idx]->allele_ct;
      PgenReader* pgrp = &(cur_mrp->pgr);
      const uint32_t vrtype = PgrGetVrtype(pgrp, read_variant_uidx);
      const uint32_t read_hphase_present = (vrtype / 0x10) & 1;
      const uint32_t read_dosage_present = !!(vrtype & 0x60);
      if ((read_allele_ct == 2) && (!read_dosage_present)) {
        pgvp->patch_01_ct = 0;
        pgvp->patch_10_ct = 0;
        pgvp->dosage_ct = 0;
        if (!read_hphase_present) {
          pgvp->phasepresent_ct = 0;
          reterr = PgrGet(sample_include, cur_mrp->pssi, read_sample_ct, read_variant_uidx, pgrp, pgvp->genovec);
        } else {
          reterr = PgrGetP(sample_include, cur_mrp->pssi, read_sample_ct, read_variant_uidx, pgrp, pgvp->genovec, pgvp->phasepresent, pgvp->phaseinfo, &pgvp->phasepresent_ct);
        }
      } else {
        reterr = PgrGetMDp(sample_include, cur_mrp->pssi, read_sample_ct, read_variant_uidx, pgrp, pgvp);
      }
      // bugfix (9 Aug 2022)
      ZeroTrailingNyps(read_sample_ct, pgvp->genovec);
      if (unlikely(reterr)) {
        PgenErrPrintNV(reterr, read_variant_uidx);
        goto MergePgenVariantNoTmpLocked_ret_1;
      }

      // Identify which samples we can blindly clobber the previous entry for,
      // and then do it.
      // For "--merge-mode nm-first" or "--merge-mode first", this set covers
      // every unlocked sample.
      // For "--merge-mode nm-match", this set covers only unlocked-and-missing
      // samples.  We check concordance of unlocked-and-nonmissing entries in a
      // subsequent loop.
      uint32_t clobber_sample_ct = unlocked_sample_ct;
      uintptr_t* clobber_sample_span = unlocked_sample_span;
      const AlleleCode* cur_allele_remap = &(master_allele_remap[rec_idx * allele_remap_stride]);
      if (unlocked_missing_set) {
        // If allele rotation is needed, do it on all of pgvp in this case.
        // (In the "--merge-mode nm-first" and "--merge-mode first" cases, we
        // wait until after we've subsetted down to the samples we'll be
        // clobbering.)
        uint32_t read_allele_idx = 0;
        for (; read_allele_idx != read_allele_ct; ++read_allele_idx) {
          if (cur_allele_remap[read_allele_idx] != read_allele_idx) {
            break;
          }
        }
        if (read_allele_idx != read_allele_ct) {
          // Allele rotation may be needed.
          if ((cur_allele_remap[0] == kMissingAlleleCode) || (cur_allele_remap[1] == kMissingAlleleCode)) {
            reterr = ValidateBiallelicVariantWithMissingCode(cur_allele_remap, read_sample_ct, (write_allele_ct == 2), pgvp);
            if (unlikely(reterr)) {
              goto MergePgenVariantNoTmpLocked_ret_1;
            }
          } else {
            const uint32_t read_sample_ctl = BitCtToWordCt(read_sample_ct);
            if (write_allele_ct == 2) {
              // Only possibility is cur_allele_remap[0] == 1,
              // cur_allele_remap[1] == 0.
              GenovecInvertUnsafe(read_sample_ct, pgvp->genovec);
              ZeroTrailingNyps(read_sample_ct, pgvp->genovec);
              if (pgvp->phasepresent_ct) {
                BitvecInvert(read_sample_ct, pgvp->phaseinfo);
              }
              if (pgvp->dosage_ct) {
                BiallelicDosage16Invert(pgvp->dosage_ct, pgvp->dosage_main);
                if (pgvp->dphase_ct) {
                  BiallelicDphase16Invert(pgvp->dphase_ct, pgvp->dphase_delta);
                }
              }
            } else {
              uintptr_t* phaseinfo_xor = nullptr;
              if (pgvp->phasepresent_ct) {
                phaseinfo_xor = mwp->phaseinfo_xor;
              }
              AlleleCode* wide_codes = mwp->wide_codes;
              PglMultiallelicSparseToDense(pgvp->genovec, pgvp->patch_01_set, pgvp->patch_01_vals, pgvp->patch_10_set, pgvp->patch_10_vals, cur_allele_remap, read_sample_ct, pgvp->patch_01_ct, pgvp->patch_10_ct, phaseinfo_xor, wide_codes);
              if (phaseinfo_xor) {
                BitvecXor(phaseinfo_xor, read_sample_ctl, pgvp->phaseinfo);
              }
              PglMultiallelicDenseToSparse(wide_codes, read_sample_ct, pgvp->genovec, pgvp->patch_01_set, pgvp->patch_01_vals, pgvp->patch_10_set, pgvp->patch_10_vals, &(pgvp->patch_01_ct), &(pgvp->patch_10_ct));
              if (pgvp->patch_01_ct || pgvp->patch_10_ct) {
                // Only need to distinguish between 2 vs. greater-than-2 for
                // now.
                read_allele_ct = 3;
              }
            }
          }
        }
        cur_allele_remap = nullptr;
        clobber_sample_span = mwp->clobber_sample_span;
        BitvecAndCopy(unlocked_missing_set, sample_span, write_sample_ctl, clobber_sample_span);
        clobber_sample_ct = PopcountWords(clobber_sample_span, write_sample_ctl);
      }
      if ((merge_mode == kMergeModeNmMatch) && (clobber_sample_ct != read_sample_ct)) {
        // Need to scan some samples for conflicts (clearing unlocked_set bit
        // whenever one exists), and maybe clobber some others.
        // Number of samples involved is likely to be a large fraction of
        // write_sample_ct, so current logic is just optimized for the dense
        // case.  Can add a separate codepath for the sparse case later.
        // 1a. Expand+permute all fields from pgvp to compare_pgvp.  (Although
        //     clobber_pgvp can be modified after pointing to pgvp above, that
        //     can only happen when clobber_sample_ct == read_sample_ct.)
        // 1b. Perform the following word-based comparisons, incrementally
        //     updating compare_mask:
        //       (2->1)(compare_pgvp->genovec ^ genovec) & compare_mask
        //       if r_patch_01_ct:
        //         (compare_pgvp->patch_01_set ^ patch_01_set) & compare_mask
        //         iterate over compare_pgvp->patch_01_vals
        //       r_patch_10_ct is similar
        //       if r_phasepresent_ct:
        //         (compare_pgvp->phasepresent ^ phasepresent) & compare_mask
        //         compare_pgvp->phasepresent & phasepresent & compare_mask &
        //           (compare_pgvp->phaseinfo ^ phaseinfo)
        //       dosage, dphase identical to patch_10
        //     This is interleaved with (1a).
        // 2. Clear unlocked_set bits.
        PgenVariant* compare_pgvp = &(mwp->pgv_midbuf);
        BitvecInvmaskCopy(unlocked_sample_span, clobber_sample_span, write_sample_ctl, mwp->unlocked_nonmissing_sample_span);
        uintptr_t* compare_mask = mwp->mask_buf;
        memcpy(compare_mask, mwp->unlocked_nonmissing_sample_span, write_sample_ctl * sizeof(intptr_t));
        uintptr_t* orig_genovec = pgvp->genovec;
        uintptr_t* r_genovec = compare_pgvp->genovec;
        ZeroWArr(write_sample_ctl2, r_genovec);
        const uint32_t* old_sample_idx_to_new = cur_mrp->old_sample_idx_to_new;
        const uint32_t read_sample_ctl2 = NypCtToWordCt(read_sample_ct);
        for (uint32_t widx = 0; widx != read_sample_ctl2; ++widx) {
          uintptr_t geno_word = orig_genovec[widx];
          if (!geno_word) {
            continue;
          }
          const uint32_t* cur_old_sample_idx_to_new = &(old_sample_idx_to_new[widx * kBitsPerWordD2]);
          do {
            const uint32_t bit_read_shift_ct = ctzw(geno_word) & (kBitsPerWord - 2);
            const uintptr_t cur_geno = (geno_word >> bit_read_shift_ct) & 3;
            const uint32_t new_sample_idx = cur_old_sample_idx_to_new[bit_read_shift_ct / 2];
            const uint32_t new_word_idx = new_sample_idx / kBitsPerWordD2;
            const uint32_t bit_write_shift_ct = 2 * (new_sample_idx % kBitsPerWordD2);
            r_genovec[new_word_idx] |= cur_geno << bit_write_shift_ct;
            geno_word &= (~(3 * k1LU)) << bit_read_shift_ct;
          } while (geno_word);
        }
        Halfword* compare_mask_hwalias = R_CAST(Halfword*, compare_mask);
        for (uint32_t widx = 0; widx != write_sample_ctl2; ++widx) {
          Halfword compare_mask_hw = compare_mask_hwalias[widx];
          if (!compare_mask_hw) {
            continue;
          }
          const uintptr_t diff_bits = r_genovec[widx] ^ genovec[widx];
          const Halfword diff_hw = PackWordToHalfwordMask5555(diff_bits | (diff_bits >> 1));
          compare_mask_hwalias[widx] = compare_mask_hw & (~diff_hw);
        }
        // possible todo: check if performance gain from switching to single
        // compare+clobber loops is significant
        if (clobber_sample_ct) {
          const Halfword* clobber_span_hwalias = R_CAST(const Halfword*, clobber_sample_span);
          for (uint32_t widx = 0; widx != write_sample_ctl2; ++widx) {
            const Halfword clobber_span_hw = clobber_span_hwalias[widx];
            if (!clobber_span_hw) {
              continue;
            }
            genovec[widx] ^= (UnpackHalfwordToWord(clobber_span_hw) * 3) & (~r_genovec[widx]);
          }
        }
        if (write_allele_ct > 2) {
          if (pgvp->patch_01_ct) {
            uintptr_t* r_patch_01_set = compare_pgvp->patch_01_set;
            AlleleCode* r_patch_01_dense = compare_pgvp->patch_01_vals;
            ZeroWArr(write_sample_ctl, r_patch_01_set);
            PermuteUpdate8bitDenseFromSparse(pgvp->patch_01_set, pgvp->patch_01_vals, old_sample_idx_to_new, read_sample_ct, pgvp->patch_01_ct, r_patch_01_set, r_patch_01_dense);
            Compare8bitDense(r_patch_01_set, r_patch_01_dense, patch_01_set, patch_01_dense, write_sample_ctl, compare_mask);
            if (clobber_sample_ct) {
              Update8bitDense(clobber_sample_span, r_patch_01_set, r_patch_01_dense, write_sample_ctl, patch_01_set, patch_01_dense);
            }
          } else {
            // Anything which is currently marked as multiallelic ref-altx, and
            // is in the current comparison set, must be a mismatch since there
            // isn't a single multiallelic ref-altx in the current source.
            BitvecInvmask(patch_01_set, write_sample_ctl, compare_mask);
          }
          if (pgvp->patch_10_ct) {
            uintptr_t* r_patch_10_set = compare_pgvp->patch_10_set;
            AlleleCode* r_patch_10_dense = compare_pgvp->patch_10_vals;
            ZeroWArr(write_sample_ctl, r_patch_10_set);
            PermuteUpdate16bitDenseFromSparse(pgvp->patch_10_set, pgvp->patch_10_vals, old_sample_idx_to_new, read_sample_ct, pgvp->patch_10_ct, r_patch_10_set, r_patch_10_dense);
            Compare16bitDense(r_patch_10_set, r_patch_10_dense, patch_10_set, patch_10_dense, write_sample_ctl, compare_mask);
            if (clobber_sample_ct) {
              Update16bitDense(clobber_sample_span, r_patch_10_set, r_patch_10_dense, write_sample_ctl, patch_10_set, patch_10_dense);
            }
          } else {
            BitvecInvmask(patch_10_set, write_sample_ctl, compare_mask);
          }
        }
        if (phasepresent) {
          if (pgvp->phasepresent_ct) {
            uintptr_t* r_phasepresent = compare_pgvp->phasepresent;
            uintptr_t* r_phaseinfo = compare_pgvp->phaseinfo;
            CopyAndPermuteHphase(pgvp->phasepresent, pgvp->phaseinfo, old_sample_idx_to_new, write_sample_ctl, pgvp->phasepresent_ct, r_phasepresent, r_phaseinfo);
            for (uint32_t widx = 0; widx != write_sample_ctl; ++widx) {
              uintptr_t compare_word = compare_mask[widx];
              if (!compare_word) {
                continue;
              }
              const uintptr_t r_phasepresent_word = r_phasepresent[widx];
              const uintptr_t w_phasepresent_word = phasepresent[widx];
              compare_word &= ~((r_phasepresent_word ^ w_phasepresent_word) | (r_phasepresent_word & w_phasepresent_word & (r_phaseinfo[widx] ^ phaseinfo[widx])));
              compare_mask[widx] = compare_word;
            }
            if (clobber_sample_ct) {
              for (uint32_t widx = 0; widx != write_sample_ctl; ++widx) {
                const uintptr_t update_word = clobber_sample_span[widx];
                if (!update_word) {
                  continue;
                }
                const uintptr_t src_present_word = r_phasepresent[widx] & update_word;
                if (!src_present_word) {
                  continue;
                }
                phasepresent[widx] |= src_present_word;
                phaseinfo[widx] |= src_present_word & r_phaseinfo[widx];
              }
            }
          } else {
            BitvecInvmask(phasepresent, write_sample_ctl, compare_mask);
          }
        }
        if (dosage_present) {
          if (pgvp->dosage_ct) {
            uintptr_t* r_dosage_present = compare_pgvp->dosage_present;
            uint16_t* r_dosage_dense = compare_pgvp->dosage_main;
            ZeroWArr(write_sample_ctl, r_dosage_present);
            PermuteUpdate16bitDenseFromSparse(pgvp->dosage_present, pgvp->dosage_main, old_sample_idx_to_new, read_sample_ct, pgvp->dosage_ct, r_dosage_present, r_dosage_dense);
            Compare16bitDense(r_dosage_present, r_dosage_dense, dosage_present, dosage_dense, write_sample_ctl, compare_mask);
            if (clobber_sample_ct) {
              Update16bitDense(clobber_sample_span, r_dosage_present, r_dosage_dense, write_sample_ctl, dosage_present, dosage_dense);
            }
          } else {
            BitvecInvmask(dosage_present, write_sample_ctl, compare_mask);
          }
          if (dphase_present) {
            if (pgvp->dphase_ct) {
              uintptr_t* r_dphase_present = compare_pgvp->dphase_present;
              int16_t* r_dphase_dense = compare_pgvp->dphase_delta;
              ZeroWArr(write_sample_ctl, r_dphase_present);
              PermuteUpdate16bitDenseFromSparse(pgvp->dphase_present, pgvp->dphase_delta, old_sample_idx_to_new, read_sample_ct, pgvp->dphase_ct, r_dphase_present, r_dphase_dense);
              Compare16bitDense(r_dphase_present, r_dphase_dense, dphase_present, dphase_dense, write_sample_ctl, compare_mask);
              if (clobber_sample_ct) {
                Update16bitDense(clobber_sample_span, r_dphase_present, r_dphase_dense, write_sample_ctl, dosage_present, dosage_dense);
              }
            } else {
              BitvecInvmask(dphase_present, write_sample_ctl, compare_mask);
            }
          }
        }
        // unlocked_set ^= unlocked_nonmissing_sample_span ^ compare_mask
        BitvecXor(mwp->unlocked_nonmissing_sample_span, write_sample_ctl, unlocked_set);
        BitvecXor(compare_mask, write_sample_ctl, unlocked_set);
      } else if (clobber_sample_ct) {
        // possible todo: only use this branch on sufficiently small
        // clobber_sample_ct
        const uint32_t clobber_sample_ctl = BitCtToWordCt(clobber_sample_ct);
        PgenVariant* clobber_pgvp = pgvp;
        uint32_t* clobber_sample_idx_to_new = cur_mrp->old_sample_idx_to_new;
        if (read_sample_ct != clobber_sample_ct) {
          // merge_mode is either 'first' or 'nm-first'.
          // Similar to subsetting logic in MakePgenThread().
          uintptr_t* clobber_subset_of_read = mwp->mask_buf;
          CopyBitarrSubset(clobber_sample_span, sample_span, read_sample_ct, clobber_subset_of_read);
          clobber_sample_idx_to_new = mwp->clobber_sample_idx_to_new;
          const uint32_t* old_sample_idx_to_new = cur_mrp->old_sample_idx_to_new;
          const uint32_t sample_idx_increasing = cur_mrp->sample_idx_increasing;
          uintptr_t read_idx_base = 0;
          uintptr_t cur_bits = clobber_subset_of_read[0];
          if (sample_idx_increasing < 2) {
            for (uint32_t clobber_idx = 0; clobber_idx != clobber_sample_ct; ++clobber_idx) {
              const uint32_t read_idx = BitIter1(clobber_subset_of_read, &read_idx_base, &cur_bits);
              clobber_sample_idx_to_new[clobber_idx] = old_sample_idx_to_new[read_idx];
            }
          } else {
            for (uint32_t clobber_idx = 0; clobber_idx != clobber_sample_ct; ++clobber_idx) {
              const uint32_t read_idx = BitIter1(clobber_subset_of_read, &read_idx_base, &cur_bits);
              clobber_sample_idx_to_new[clobber_idx] = read_idx;
            }
          }
          clobber_pgvp = &(mwp->pgv_midbuf);
          CopyNyparrNonemptySubset(pgvp->genovec, clobber_subset_of_read, read_sample_ct, clobber_sample_ct, clobber_pgvp->genovec);
          if (read_allele_ct > 2) {
            clobber_pgvp->patch_01_ct = Copy1bit8Subset(pgvp->patch_01_set, pgvp->patch_01_vals, clobber_subset_of_read, pgvp->patch_01_ct, clobber_sample_ct, clobber_pgvp->patch_01_set, clobber_pgvp->patch_01_vals);
            clobber_pgvp->patch_10_ct = Copy1bit16Subset(pgvp->patch_10_set, pgvp->patch_10_vals, clobber_subset_of_read, pgvp->patch_10_ct, clobber_sample_ct, clobber_pgvp->patch_10_set, clobber_pgvp->patch_10_vals);
          }
          clobber_pgvp->phasepresent_ct = 0;
          if (pgvp->phasepresent_ct) {
            CopyBitarrSubset(pgvp->phasepresent, clobber_subset_of_read, clobber_sample_ct, clobber_pgvp->phasepresent);
            clobber_pgvp->phasepresent_ct = PopcountWords(clobber_pgvp->phasepresent, clobber_sample_ctl);
            if (clobber_pgvp->phasepresent_ct) {
              CopyBitarrSubset(pgvp->phaseinfo, clobber_subset_of_read, clobber_sample_ct, clobber_pgvp->phaseinfo);
            }
          }
          clobber_pgvp->dosage_ct = 0;
          if (pgvp->dosage_ct) {
            clobber_pgvp->dosage_ct = Copy1bit16Subset(pgvp->dosage_present, pgvp->dosage_main, clobber_subset_of_read, pgvp->dosage_ct, clobber_sample_ct, clobber_pgvp->dosage_present, clobber_pgvp->dosage_main);
            clobber_pgvp->dphase_ct = Copy1bit16Subset(pgvp->dphase_present, pgvp->dphase_delta, clobber_subset_of_read, pgvp->dphase_ct, clobber_sample_ct, clobber_pgvp->dphase_present, clobber_pgvp->dphase_delta);
          }
        }
        if (cur_allele_remap && ((cur_allele_remap[0] == kMissingAlleleCode) || (cur_allele_remap[1] == kMissingAlleleCode))) {
          reterr = ValidateBiallelicVariantWithMissingCode(cur_allele_remap, clobber_sample_ct, (write_allele_ct == 2), clobber_pgvp);
          if (unlikely(reterr)) {
            goto MergePgenVariantNoTmpLocked_ret_1;
          }
        }
        uintptr_t* r_genovec = clobber_pgvp->genovec;
        const uint32_t r_phasepresent_ct = clobber_pgvp->phasepresent_ct;
        if (write_allele_ct == 2) {
          const uint32_t r_dosage_ct = clobber_pgvp->dosage_ct;
          uint32_t r_dphase_ct = r_dosage_ct? clobber_pgvp->dphase_ct : 0;
          if (cur_allele_remap && (cur_allele_remap[0] == 1) && (cur_allele_remap[1] == 0)) {
            // This is the only remaining case where clobber_pgvp is still
            // flipped relative to the write buffers; we already flipped the
            // kMissingAlleleCode cases.
            GenovecInvertUnsafe(clobber_sample_ct, r_genovec);
            ZeroTrailingNyps(clobber_sample_ct, r_genovec);
            if (r_phasepresent_ct) {
              BitvecInvert(clobber_sample_ctl, clobber_pgvp->phaseinfo);
            }
            if (r_dosage_ct) {
              BiallelicDosage16Invert(r_dosage_ct, clobber_pgvp->dosage_main);
              if (r_dphase_ct) {
                BiallelicDphase16Invert(r_dphase_ct, clobber_pgvp->dphase_delta);
              }
            }
          }
          PermuteUpdateGenovec(r_genovec, clobber_sample_idx_to_new, clobber_sample_ct, genovec);
          if (r_phasepresent_ct) {
            PermuteUpdateHphase(clobber_pgvp->phasepresent, clobber_pgvp->phaseinfo, clobber_sample_idx_to_new, r_phasepresent_ct, phasepresent, phaseinfo);
          }
          if (r_dosage_ct) {
            PermuteUpdate16bitDenseFromSparse(clobber_pgvp->dosage_present, clobber_pgvp->dosage_main, clobber_sample_idx_to_new, clobber_sample_ct, r_dosage_ct, dosage_present, dosage_dense);
            if (r_dphase_ct) {
              PermuteUpdate16bitDenseFromSparse(clobber_pgvp->dphase_present, clobber_pgvp->dphase_delta, clobber_sample_idx_to_new, clobber_sample_ct, r_dosage_ct, dphase_present, dphase_dense);
            }
          }
        } else {
          // Multiallelic variant.
          uintptr_t* r_patch_01_set = clobber_pgvp->patch_01_set;
          AlleleCode* r_patch_01_vals = clobber_pgvp->patch_01_vals;
          uintptr_t* r_patch_10_set = clobber_pgvp->patch_10_set;
          AlleleCode* r_patch_10_vals = clobber_pgvp->patch_10_vals;
          uint32_t r_patch_01_ct = 0;
          uint32_t r_patch_10_ct = 0;
          if (read_allele_ct > 2) {
            r_patch_01_ct = clobber_pgvp->patch_01_ct;
            r_patch_10_ct = clobber_pgvp->patch_10_ct;
          }
          if (cur_allele_remap) {
            uint32_t read_allele_idx = 0;
            for (; read_allele_idx != read_allele_ct; ++read_allele_idx) {
              if (cur_allele_remap[read_allele_idx] != read_allele_idx) {
                break;
              }
            }
            if (read_allele_idx != read_allele_ct) {
              // Allele rotation needed.
              uintptr_t* phaseinfo_xor = nullptr;
              if (r_phasepresent_ct) {
                phaseinfo_xor = mwp->phaseinfo_xor;
              }
              AlleleCode* wide_codes = mwp->wide_codes;
              PglMultiallelicSparseToDense(r_genovec, r_patch_01_set, r_patch_01_vals, r_patch_10_set, r_patch_10_vals, cur_allele_remap, clobber_sample_ct, r_patch_01_ct, r_patch_10_ct, phaseinfo_xor, wide_codes);
              if (phaseinfo_xor) {
                BitvecXor(phaseinfo_xor, clobber_sample_ctl, clobber_pgvp->phaseinfo);
              }
              PglMultiallelicDenseToSparse(wide_codes, clobber_sample_ct, r_genovec, r_patch_01_set, r_patch_01_vals, r_patch_10_set, r_patch_10_vals, &r_patch_01_ct, &r_patch_10_ct);
            }
          }
          PermuteUpdateGenovec(r_genovec, clobber_sample_idx_to_new, clobber_sample_ct, genovec);
          if (r_patch_01_ct) {
            PermuteUpdate8bitDenseFromSparse(r_patch_01_set, r_patch_01_vals, clobber_sample_idx_to_new, clobber_sample_ct, r_patch_01_ct, patch_01_set, patch_01_dense);
          }
          if (r_patch_10_ct) {
            PermuteUpdate16bitDenseFromSparse(r_patch_10_set, r_patch_10_vals, clobber_sample_idx_to_new, clobber_sample_ct, r_patch_10_ct, patch_10_set, patch_10_dense);
          }
          if (r_phasepresent_ct) {
            PermuteUpdateHphase(clobber_pgvp->phasepresent, clobber_pgvp->phaseinfo, clobber_sample_idx_to_new, r_phasepresent_ct, phasepresent, phaseinfo);
          }
        }
      }
      if (clobber_sample_ct && (rec_idx + 1 != merge_rec_ct)) {
        // Update unlocked_set for modes 'first'/'nm-first', or
        // unlocked_missing_set for mode 'nm-match'.
        //   'first': clear all bits intersecting clobber_sample_span.
        //   'nm-first', 'nm-match': clear all bits corresponding to nonmissing
        //                           genotypes.
        // Skip this on last input file, since this only affects eligibility
        // for future clobbers, which don't exist.
        uintptr_t* untouched_set = unlocked_missing_set? unlocked_missing_set : unlocked_set;
        if (merge_mode == kMergeModeFirst) {
          BitvecInvmask(clobber_sample_span, write_sample_ctl, untouched_set);
        } else {
          GenoarrToMissingnessUnsafe(genovec, write_sample_ct, untouched_set);
          if (dosage_present) {
            // Possible for dosage to be present when hardcall is missing.
            BitvecInvmask(dosage_present, write_sample_ctl, untouched_set);
          }
        }
      }
    }
    if (merge_mode == kMergeModeNmMatch) {
      const Halfword* unlocked_set_hwalias = R_CAST(Halfword*, unlocked_set);
      for (uint32_t widx = 0; widx != write_sample_ctl2; ++widx) {
        const Halfword locked_hw = ~unlocked_set_hwalias[widx];
        if (!locked_hw) {
          continue;
        }
        genovec[widx] |= UnpackHalfwordToWord(locked_hw) * 3;
      }
      ZeroTrailingNyps(write_sample_ct, genovec);
      if (write_allele_ct > 2) {
        BitvecAnd(unlocked_set, write_sample_ctl, patch_01_set);
        BitvecAnd(unlocked_set, write_sample_ctl, patch_10_set);
      }
      if (phasepresent) {
        BitvecAnd(unlocked_set, write_sample_ctl, phasepresent);
      }
      if (dosage_present) {
        BitvecAnd(unlocked_set, write_sample_ctl, dosage_present);
        if (dphase_present) {
          BitvecAnd(unlocked_set, write_sample_ctl, dphase_present);
        }
      }
    }
    // Collapse dense representations to sparse when the latter is expected by
    // the writer.
    uint32_t patch_01_ct = 0;
    uint32_t patch_10_ct = 0;
    if (write_allele_ct > 2) {
      patch_01_ct = Dense8bitToSparse(patch_01_set, write_sample_ctl, patch_01_dense);
      patch_10_ct = Dense16bitToSparse(patch_10_set, write_sample_ctl, patch_10_dense);
    }
    // rename for clarity
    AlleleCode* patch_01_vals = patch_01_dense;
    AlleleCode* patch_10_vals = patch_10_dense;

    uint32_t dosage_ct = 0;
    uint32_t dphase_ct = 0;
    if (dosage_present) {
      dosage_ct = Dense16bitToSparse(dosage_present, write_sample_ctl, dosage_dense);
      if (dphase_present) {
        dphase_ct = Dense16bitToSparse(dphase_present, write_sample_ctl, dphase_dense);
      }
    }
    uint16_t* dosage_main = dosage_dense;
    int16_t* dphase_delta = dphase_dense;

    uint32_t phasepresent_ct = 0;
    if (phasepresent) {
      phasepresent_ct = PopcountWords(phasepresent, write_sample_ctl);
    }
    if (write_allele_ct == 2) {
      if (!dosage_ct) {
        if (!phasepresent_ct) {
          reterr = SpgwAppendBiallelicGenovec(genovec, spgwp);
        } else {
          reterr = SpgwAppendBiallelicGenovecHphase(genovec, phasepresent, phaseinfo, spgwp);
        }
      } else {
        if ((!phasepresent_ct) && (!dphase_ct)) {
          reterr = SpgwAppendBiallelicGenovecDosage16(genovec, dosage_present, dosage_main, dosage_ct, spgwp);
        } else {
          reterr = SpgwAppendBiallelicGenovecDphase16(genovec, phasepresent, phaseinfo, dosage_present, dphase_present, dosage_main, dphase_delta, dosage_ct, dphase_ct, spgwp);
        }
      }
    } else {
      assert(!dosage_ct); // not yet supported
      if (!phasepresent_ct) {
        reterr = SpgwAppendMultiallelicSparse(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, write_allele_ct, patch_01_ct, patch_10_ct, spgwp);
      } else {
        reterr = SpgwAppendMultiallelicGenovecHphase(genovec, patch_01_set, patch_01_vals, patch_10_set, patch_10_vals, phasepresent, phaseinfo, write_allele_ct, patch_01_ct, patch_10_ct, spgwp);
      }
    }
    if (unlikely(reterr)) {
      goto MergePgenVariantNoTmpLocked_ret_1;
    }
  }
 MergePgenVariantNoTmpLocked_ret_1:
  return reterr;
}

// TODO: --merge-mode 'first'/'nm-match' + intermediate-file support

typedef struct SamePosPvarRecordAsorterStruct {
  SamePosPvarRecord* pp;
#ifdef __cplusplus
  bool operator<(const struct SamePosPvarRecordAsorterStruct& rhs) const {
    const int32_t strcmp_result = strcmp_overread(pp->variant_id, rhs.pp->variant_id);
    if (strcmp_result) {
      return strcmp_result < 0;
    }
    return pp->secondary_key < rhs.pp->secondary_key;
  }
#endif
} SamePosPvarRecordAsorter;

typedef struct SamePosPvarRecordNsorterStruct {
  SamePosPvarRecord* pp;
#ifdef __cplusplus
  bool operator<(const struct SamePosPvarRecordNsorterStruct& rhs) const {
    const int32_t strcmp_result = strcmp_natural_uncasted(pp->variant_id, rhs.pp->variant_id);
    if (strcmp_result) {
      return strcmp_result < 0;
    }
    return pp->secondary_key < rhs.pp->secondary_key;
  }
#endif
} SamePosPvarRecordNsorter;

#ifndef __cplusplus
int32_t SamePosPvarRecordAcmp(const void* r1, const void* r2) {
  const SamePosPvarRecord* rec1 = S_CAST(const SamePosPvarRecordAsorter*, r1)->pp;
  const SamePosPvarRecord* rec2 = S_CAST(const SamePosPvarRecordAsorter*, r2)->pp;
  const int32_t strcmp_result = strcmp_overread(rec1->variant_id, rec2->variant_id);
  if (strcmp_result) {
    return strcmp_result;
  }
  if (rec1->secondary_key < rec2->secondary_key) {
    return -1;
  }
  return (rec1->secondary_key > rec2->secondary_key);
}

int32_t SamePosPvarRecordNcmp(const void* r1, const void* r2) {
  const SamePosPvarRecord* rec1 = S_CAST(const SamePosPvarRecordNsorter*, r1)->pp;
  const SamePosPvarRecord* rec2 = S_CAST(const SamePosPvarRecordNsorter*, r2)->pp;
  const int32_t strcmp_result = strcmp_natural_uncasted(rec1->variant_id, rec2->variant_id);
  if (strcmp_result) {
    return strcmp_result;
  }
  if (rec1->secondary_key < rec2->secondary_key) {
    return -1;
  }
  return (rec1->secondary_key > rec2->secondary_key);
}
#endif

PglErr ConcatPvariantPos(int32_t cur_bp, uintptr_t variant_ct, PvariantPosMergeContext* ppmcp, SamePosPvarRecord** same_pos_records, MergeReader* mrp, MergeWriter* mwp) {
  if (!variant_ct) {
    return kPglRetSuccess;
  }
  if (ppmcp->sort_vars_mode == kSortAscii) {
    SamePosPvarRecordAsorter* asorter = R_CAST(SamePosPvarRecordAsorter*, same_pos_records);
    STD_SORT(variant_ct, SamePosPvarRecordAcmp, asorter);
  } else {
    SamePosPvarRecordNsorter* nsorter = R_CAST(SamePosPvarRecordNsorter*, same_pos_records);
    STD_SORT(variant_ct, SamePosPvarRecordNcmp, nsorter);
  }
  ppmcp->pmc.cur_bp = cur_bp;
  uintptr_t* write_nonref_flags = ppmcp->write_nonref_flags;
  uint32_t write_variant_idx = ppmcp->write_variant_idx;
  uint32_t next_print_variant_idx = ppmcp->next_print_variant_idx;
  uintptr_t* write_allele_idx_offsets = ppmcp->write_allele_idx_offsets;

  for (uintptr_t rec_idx_start = 0; ; ) {
    char* cur_variant_id = same_pos_records[rec_idx_start]->variant_id;
    const uint32_t cur_variant_id_slen = strlen(cur_variant_id);
    uintptr_t rec_idx_end = rec_idx_start + 1;
    for (; rec_idx_end != variant_ct; ++rec_idx_end) {
      if (!strequal_unsafe(same_pos_records[rec_idx_end]->variant_id, cur_variant_id, cur_variant_id_slen)) {
        break;
      }
    }
    SamePosPvarRecord** same_id_records = &(same_pos_records[rec_idx_start]);
    const uintptr_t merge_rec_ct = rec_idx_end - rec_idx_start;
    uint32_t is_pr = 0;
    uint32_t allele_ct;
    uint64_t cur_line_blen;
    PglErr reterr = MergePvariant(merge_rec_ct, &ppmcp->pmc, same_id_records, write_nonref_flags? (&is_pr) : nullptr, &allele_ct, &cur_line_blen);
    if (unlikely(reterr)) {
      return reterr;
    }
    // bugfix (8 Sep 2021): forgot that allele_ct == 0 indicates variant skip
    if (allele_ct) {
      if (unlikely(cur_line_blen > kMaxLongLine)) {
        logerrprintfww("Error: Merged .pvar entry for variant '%s' at %s:%d is too long for this " PROG_NAME_STR " build.\n", cur_variant_id, ppmcp->pmc.chr_buf, cur_bp);
        return kPglRetNotYetSupported;
      }
      if (write_allele_idx_offsets) {
        write_allele_idx_offsets[write_variant_idx + 1] = write_allele_idx_offsets[write_variant_idx] + allele_ct;
      }
      if (write_nonref_flags) {
        AssignBit(write_variant_idx, is_pr, write_nonref_flags);
      }

      reterr = MergePgenVariantNoTmpLocked(same_id_records, ppmcp->pmc.allele_remap, merge_rec_ct, allele_ct, ppmcp->pmc.read_max_allele_ct, &mrp, mwp);
      if (unlikely(reterr)) {
        return reterr;
      }

      ++write_variant_idx;
      if (write_variant_idx == next_print_variant_idx) {
        printf("\rConcatenating... %u/%u variants complete.", write_variant_idx, ppmcp->write_variant_ct);
        fflush(stdout);
        next_print_variant_idx += 10000;
        ppmcp->next_print_variant_idx = next_print_variant_idx;
      }
    }
    if (rec_idx_end == variant_ct) {
      ppmcp->write_variant_idx = write_variant_idx;
      return kPglRetSuccess;
    }
    rec_idx_start = rec_idx_end;
  }
}

void CleanupPvariantPosMergeContext(PvariantPosMergeContext* ppmcp) {
  CswriteCloseCond(&ppmcp->pmc.css, ppmcp->pmc.cswritep);
}

// This can actually deviate from pure concatenation: same-position variants
// are reordered by ID, and same-position same-ID variants are merged.  The
// distinction from the general case is that we never need to have more than
// one .pvar + .pgen open for reading at a time.
PglErr PmergeConcat(const PmergeInfo* pmip, const SampleIdInfo* siip, const ChrInfo* cip, const PmergeInputFilesetLl* filesets, const char* missing_varid_match, const char* const* info_keys, const uint32_t* info_keys_htable, uint32_t sample_ct, FamCol fam_cols, uintptr_t fileset_ct, uint32_t psam_linebuf_capacity, uint32_t missing_varid_match_slen, uint32_t info_key_ct, uint32_t info_keys_htable_size, uint32_t info_conflict_present, char input_missing_geno_char, uint32_t max_thread_ct, SortMode sort_vars_mode, VaridTemplate* varid_templatep, VaridTemplate* varid_multi_templatep, VaridTemplate* varid_multi_nonsnp_templatep, char* outname, char* outname_end) {
  // Don't need to reset bigstack at function end, since Pmerge() will do it.
  const char* read_pgen_fname = nullptr;
  const char* read_pvar_fname = nullptr;
  PglErr reterr = kPglRetSuccess;
  PvariantPosMergeContext ppmc;
  PreinitPvariantPosMergeContext(&ppmc);
  uintptr_t pvar_line_idx = 0;
  TextStream pvar_txs;
  PreinitTextStream(&pvar_txs);
  PgenFileInfo pgfi;
  MergeReader mr;
  MergeWriter mw;
  PreinitPgfi(&pgfi);
  PreinitPgr(&mr.pgr);
  PreinitSpgw(&mw.spgw);
  {
    // 1. Scan .pgen headers, to determine appropriate write_gflags.
    // 2. Initialize single-threaded .pgen writer.  (Possible todo:
    //    support multithreaded writer when sufficient memory is available.)
    // 3. Iterate through filesets:
    //    a. Scan .psam, save sample subset/order.
    //    b. Scan through .pvar and .pgen simultaneously.
    uintptr_t write_variant_ct = 0;
    uint32_t write_qual = 0;
    uint32_t write_filter = 0;
    uint32_t write_info = 0;
    uint32_t write_cm = 0;
    uintptr_t overflow_buf_size = 0;
    uint32_t vrtype_8bit_needed = 0;
    // 1 = all known, 2 = all provisional-REF, 3 = enough evidence for mixed
    uint32_t nonref_flags_storage = 0;

    uint32_t read_max_allele_ct = 2;
    uint32_t write_max_allele_ct = 2;
    uint32_t max_single_pos_ct = 1;
    uintptr_t max_single_pos_blen = 0;
    uint32_t read_max_nonpass_filter_ct = 0;
    const PmergeInputFilesetLl* filesets_iter = filesets;
    for (uintptr_t fileset_idx = 0; fileset_idx != fileset_ct; ++fileset_idx) {
      write_variant_ct += filesets_iter->write_nondoomed_variant_ct;
      write_qual |= filesets_iter->nm_qual_exists;
      write_filter |= filesets_iter->nm_filter_exists;
      write_info |= filesets_iter->nm_info_exists;
      write_cm |= filesets_iter->nz_cm_present;
      if (overflow_buf_size < filesets_iter->max_pvar_line_blen) {
        overflow_buf_size = filesets_iter->max_pvar_line_blen;
      }
      if (read_max_allele_ct < filesets_iter->read_max_allele_ct) {
        read_max_allele_ct = filesets_iter->read_max_allele_ct;
      }
      if (write_max_allele_ct < filesets_iter->write_nondoomed_max_allele_ct) {
        write_max_allele_ct = filesets_iter->write_nondoomed_max_allele_ct;
      }
      if (max_single_pos_ct < filesets_iter->max_single_pos_ct) {
        max_single_pos_ct = filesets_iter->max_single_pos_ct;
      }
      if (max_single_pos_blen < filesets_iter->max_single_pos_blen) {
        max_single_pos_blen = filesets_iter->max_single_pos_blen;
      }
      if (read_max_nonpass_filter_ct < filesets_iter->read_max_nonpass_filter_ct) {
        read_max_nonpass_filter_ct = filesets_iter->read_max_nonpass_filter_ct;
      }
      vrtype_8bit_needed |= filesets_iter->vrtype_8bit_needed;
      const uint32_t cur_nonref_flags_storage = filesets_iter->nonref_flags_storage;
      if (!cur_nonref_flags_storage) {
        nonref_flags_storage = 3;
      } else {
        nonref_flags_storage |= cur_nonref_flags_storage;
      }
      filesets_iter = filesets_iter->next;
    }
    if (unlikely(write_variant_ct > kPglMaxVariantCt)) {
      logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
      goto PmergeConcat_ret_INCONSISTENT_INPUT;
    }
    if (unlikely(!write_variant_ct)) {
      logerrputs("Error: All variants filtered out by --merge-max-alleles.\n");
      goto PmergeConcat_ret_INCONSISTENT_INPUT;
    }

    // a few extra bytes for miscellaneous delimiters
    overflow_buf_size += 32;
    if (info_key_ct) {
      // We'd rather not be forced to perform a bunch of write-buffer flushes
      // at the end of each INFO Number=A or =R field, or when expanding a
      // Number=<positive #> field when there's a conflict.
      // Each INFO Number=A requires up to 2 * (n-2) extra bytes, and each INFO
      // Number=R requires up to 2 * (n-1).
      // Each INFO Number=k requires up to 2 * (k-1) extra bytes.
      uint32_t info_ra_cts[2];
      info_ra_cts[0] = 0; // R
      info_ra_cts[1] = 0; // A
      uintptr_t num_m1_sum = 0;
      for (uint32_t info_key_idx = 0; info_key_idx != info_key_ct; ++info_key_idx) {
        const int32_t info_vtype = const_container_of(info_keys[info_key_idx], InfoVtype, key)->num;
        if (!IsInfoVtypeARSkip(info_vtype)) {
          if (info_vtype > 1) {
            num_m1_sum += info_vtype - 1;
          }
          continue;
        }
        info_ra_cts[info_vtype - kInfoVtypeR] += 1;
      }
      const uintptr_t max_extra_cost = 2 * (info_ra_cts[1] * (write_max_allele_ct - 2) + info_ra_cts[0] * (write_max_allele_ct - 1) + num_m1_sum);
      overflow_buf_size += max_extra_cost;
    }
    if (overflow_buf_size < kCompressStreamBlock) {
      overflow_buf_size = kCompressStreamBlock;
    }
    overflow_buf_size += kCompressStreamBlock;
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    const uint32_t pvar_zst = (pmip->flags / kfPmergeOutputVzs) & 1;
    if (pvar_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    uintptr_t line_idx_body_start = 0;
    reterr = InitPvariantPosMergeContext(pmip, outname, info_keys, info_keys_htable, &read_pvar_fname, &line_idx_body_start, write_qual, write_filter, write_info, write_cm, info_key_ct, info_keys_htable_size, info_conflict_present, sort_vars_mode, 0, 1, pvar_zst, overflow_buf_size, read_max_allele_ct, write_max_allele_ct, max_single_pos_ct, read_max_nonpass_filter_ct, &ppmc);
    if (unlikely(reterr)) {
      goto PmergeConcat_ret_1;
    }

    ppmc.write_variant_ct = write_variant_ct;
    if (write_max_allele_ct > 2) {
      // bugfix (13 Apr 2021): forgot +1
      if (bigstack_alloc_w(write_variant_ct + 1, &ppmc.write_allele_idx_offsets)) {
        goto PmergeConcat_ret_NOMEM;
      }
      ppmc.write_allele_idx_offsets[0] = 0;
    }
    const uint32_t write_variant_ctl = BitCtToWordCt(write_variant_ct);
    if (nonref_flags_storage == 3) {
      if (bigstack_calloc_w(write_variant_ctl, &ppmc.write_nonref_flags)) {
        goto PmergeConcat_ret_NOMEM;
      }
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    const PgenGlobalFlags write_gflags = vrtype_8bit_needed? (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePresent | kfPgenGlobalDosagePhasePresent) : kfPgenGlobal0;
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    // bugfix (1 Jun 2022): in multiallelic case, contents of
    // write_allele_idx_offsets have not been initialized, and this triggered
    // an assert-failure in debug builds or a write_max_allele_ct
    // miscalculation otherwise (practically always a harmless overestimate)
    reterr = SpgwInitPhase1(outname, nullptr, ppmc.write_nonref_flags, write_variant_ct, sample_ct, write_max_allele_ct, kPgenWriteBackwardSeek, write_gflags, nonref_flags_storage, &mw.spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, outname, strerror(errno));
      }
      goto PmergeConcat_ret_1;
    }
    const uint32_t sample_id_htable_size = GetHtableMinSize(sample_ct);
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    unsigned char* spgw_alloc;
    uint32_t* sample_id_htable;
    uint32_t* old_sample_idx_to_new_buf;
    if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc) ||
                 bigstack_alloc_u32(sample_id_htable_size, &sample_id_htable) ||
                 bigstack_alloc_u32(sample_ct, &old_sample_idx_to_new_buf) ||
                 bigstack_alloc_w(sample_ctl2, &mw.genovec))) {
      goto PmergeConcat_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &mw.spgw, spgw_alloc);
    mw.patch_01_set = nullptr;
    mw.patch_01_vals = nullptr;
    mw.patch_10_set = nullptr;
    mw.patch_10_vals = nullptr;
    if (write_max_allele_ct > 2) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &mw.patch_01_set) ||
                   bigstack_alloc_ac(sample_ct, &mw.patch_01_vals) ||
                   bigstack_alloc_w(sample_ctl, &mw.patch_10_set) ||
                   bigstack_alloc_ac(2 * sample_ct, &mw.patch_10_vals) ||
                   bigstack_alloc_ac(2 * sample_ct, &mw.wide_codes))) {
        goto PmergeConcat_ret_NOMEM;
      }
    }
    mw.phasepresent = nullptr;
    mw.phaseinfo = nullptr;
    mw.dosage_present = nullptr;
    mw.dosage_main = nullptr;
    mw.dphase_present = nullptr;
    mw.dphase_delta = nullptr;
    if (vrtype_8bit_needed) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &mw.phasepresent) ||
                   bigstack_alloc_w(sample_ctl, &mw.phaseinfo) ||
                   bigstack_alloc_w(sample_ctl, &mw.dosage_present) ||
                   bigstack_alloc_dosage(sample_ct, &mw.dosage_main) ||
                   bigstack_alloc_w(sample_ctl, &mw.dphase_present) ||
                   bigstack_alloc_dphase(sample_ct, &mw.dphase_delta) ||
                   bigstack_alloc_w(sample_ctl, &mw.phaseinfo_xor))) {
        goto PmergeConcat_ret_NOMEM;
      }
    }
    // pgv_readbuf reinitialized for each file we're reading from
    if (unlikely(bigstack_alloc_w(sample_ctl, &mw.unlocked_set) ||
                 bigstack_alloc_w(sample_ctl, &mw.unlocked_sample_span) ||
                 bigstack_alloc_u32(sample_ct, &mw.clobber_sample_idx_to_new) ||
                 bigstack_alloc_w(sample_ctl, &mw.mask_buf) ||
                 BigstackAllocPgv(sample_ct, write_max_allele_ct > 2, write_gflags, &mw.pgv_midbuf))) {
      goto PmergeConcat_ret_NOMEM;
    }
    mw.unlocked_missing_set = nullptr;
    mw.clobber_sample_span = nullptr;
    mw.unlocked_nonmissing_sample_span = nullptr;
    if (pmip->merge_mode == kMergeModeNmMatch) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &mw.unlocked_missing_set) ||
                   bigstack_alloc_w(sample_ctl, &mw.clobber_sample_span) ||
                   bigstack_alloc_w(sample_ctl, &mw.unlocked_nonmissing_sample_span))) {
        goto PmergeConcat_ret_NOMEM;
      }
    }
    mw.merge_mode = pmip->merge_mode;

    InitXidHtable(siip, sample_ct, sample_id_htable_size, sample_id_htable, g_textbuf);

    unsigned char* bigstack_mark = g_bigstack_base;
    filesets_iter = filesets;
    logputs("Concatenating... ");
    printf("0/%" PRIuPTR " variant%s complete.", write_variant_ct, (write_variant_ct == 1)? "" : "s");
    fflush(stdout);
    for (uintptr_t fileset_idx = 0; fileset_idx != fileset_ct; ++fileset_idx, filesets_iter = filesets_iter->next) {
      if (!filesets_iter->write_nondoomed_variant_ct) {
        continue;
      }
      BigstackReset(bigstack_mark);
      const uint32_t read_sample_ct = filesets_iter->read_sample_ct;
      const uint32_t read_sample_ctl = BitCtToWordCt(read_sample_ct);
      uint32_t* read_cumulative_popcounts;
      if (unlikely(bigstack_alloc_w(read_sample_ctl, &mr.sample_include) ||
                   bigstack_alloc_u32(read_sample_ctl, &read_cumulative_popcounts) ||
                   bigstack_alloc_w(sample_ctl, &mr.sample_span))) {
        goto PmergeConcat_ret_NOMEM;
      }
      mr.old_sample_idx_to_new = old_sample_idx_to_new_buf;
      uint32_t cur_write_sample_ct;
      reterr = ScrapeSampleOrder(filesets_iter->psam_fname, siip, sample_id_htable, read_sample_ct, sample_ct, sample_id_htable_size, fam_cols, 0, psam_linebuf_capacity, max_thread_ct, mr.old_sample_idx_to_new, &mr.sample_idx_increasing, &cur_write_sample_ct, mr.sample_include, mr.sample_span);
      if (unlikely(reterr)) {
        goto PmergeConcat_ret_1;
      }
      FillCumulativePopcounts(mr.sample_include, read_sample_ctl, read_cumulative_popcounts);
      if (mr.sample_idx_increasing && (cur_write_sample_ct == sample_ct)) {
        mr.sample_idx_increasing = 2;
      }
      mr.sample_ct = cur_write_sample_ct;

      read_pgen_fname = filesets_iter->pgen_fname;
      const uint32_t read_variant_ct = filesets_iter->read_variant_ct;
      PgenHeaderCtrl header_ctrl;
      uintptr_t cur_alloc_cacheline_ct;
      reterr = PgfiInitPhase1(read_pgen_fname, nullptr, read_variant_ct, read_sample_ct, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, g_logbuf);
      if (unlikely(reterr)) {
        if (reterr == kPglRetInconsistentInput) {
          // .pgen was not checked for consistency with .pvar on the first
          // pass.
          logputs("\n");
          WordWrapB(0);
          logerrputsb();
          goto PmergeConcat_ret_1;
        }
        goto PmergeConcat_ret_PGEN_REWIND_FAIL_N;
      }
      unsigned char* pgfi_alloc;
      if (unlikely(bigstack_alloc_uc(cur_alloc_cacheline_ct * kCacheline, &pgfi_alloc))) {
        goto PmergeConcat_ret_NOMEM;
      }
      if ((header_ctrl & 192) == 192) {
        if (unlikely(bigstack_alloc_w(BitCtToWordCt(read_variant_ct), &pgfi.nonref_flags))) {
          goto PmergeConcat_ret_NOMEM;
        }
      }
      read_pvar_fname = filesets_iter->pvar_fname;
      if (filesets_iter->read_max_allele_ct > 2) {
        if (bigstack_alloc_w(read_variant_ct + 1, &pgfi.allele_idx_offsets)) {
          goto PmergeConcat_ret_NOMEM;
        }
        pgfi.allele_idx_offsets[0] = 0;
        pgfi.max_allele_ct = filesets_iter->read_max_allele_ct;
      }
      uint32_t max_vrec_width;
      reterr = PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, read_variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, &cur_alloc_cacheline_ct, g_logbuf);
      if (unlikely(reterr)) {
        logputs("\n");
        WordWrapB(0);
        logerrputsb();
        goto PmergeConcat_ret_1;
      }
      unsigned char* pgr_alloc;
      if (unlikely(bigstack_alloc_uc(cur_alloc_cacheline_ct * kCacheline, &pgr_alloc))) {
        goto PmergeConcat_ret_NOMEM;
      }
      reterr = PgrInit(read_pgen_fname, max_vrec_width, &pgfi, &mr.pgr, pgr_alloc);
      if (unlikely(reterr)) {
        goto PmergeConcat_ret_PGEN_REWIND_FAIL_N;
      }
      PgrSetSampleSubsetIndex(read_cumulative_popcounts, &mr.pgr, &mr.pssi);
      // Must check write_max_allele_ct instead of just whether the input file
      // has multiallelic variants, since we may need to rotate a biallelic
      // variant into a "multiallelic variant" in these buffers.
      if (unlikely(BigstackAllocPgv(read_sample_ct, write_max_allele_ct > 2, write_gflags, &mw.pgv_readbuf))) {
        goto PmergeConcat_ret_NOMEM;
      }
      reterr = InitTextStream(read_pvar_fname, MAXV(filesets_iter->max_pvar_line_blen, kDecompressMinBlen), 1, &pvar_txs);
      if (unlikely(reterr)) {
        goto PmergeConcat_ret_PVAR_TSTREAM_REWIND_FAIL_N;
      }
      char* line_start = TextLineEnd(&pvar_txs);
      for (pvar_line_idx = 1; ; ++pvar_line_idx) {
        if (unlikely(!TextGetUnsafe2(&pvar_txs, &line_start))) {
          reterr = TextStreamRawErrcode(&pvar_txs);
          goto PmergeConcat_ret_PVAR_TSTREAM_REWIND_FAIL_N;
        }
        if ((line_start[0] != '#') || tokequal_k(line_start, "#CHROM")) {
          break;
        }
        line_start = AdvPastDelim(line_start, '\n');
      }
      uint32_t col_skips[8];
      uint32_t col_types[8];
      const uint32_t read_qual = filesets_iter->nm_qual_exists;
      const uint32_t read_filter = filesets_iter->nm_filter_exists;
      const uint32_t read_info_pr = filesets_iter->pvar_info_pr_exists;
      const uint32_t read_info = read_info_pr | filesets_iter->nm_info_exists;
      const uint32_t pgen_pr_status_base = 2 * read_info_pr + (filesets_iter->nonref_flags_storage == 2);
      const uint32_t read_cm = filesets_iter->nz_cm_present;
      uint32_t relevant_postchr_col_ct;
      if (line_start[0] == '#') {
        char* token_end = &(line_start[6]);
        relevant_postchr_col_ct = 0;
        for (uint32_t col_idx = 1; ; ++col_idx) {
          char* token_start = FirstNonTspace(token_end);
          if (IsEolnKns(*token_start)) {
            break;
          }
          token_end = CurTokenEnd(token_start);
          const uint32_t token_slen = token_end - token_start;
          uint32_t cur_col_type;
          if (token_slen <= 3) {
            if (token_slen == 3) {
              if (memequal_sk(token_start, "POS")) {
                cur_col_type = 0;
              } else if (memequal_sk(token_start, "REF")) {
                cur_col_type = 2;
              } else if (memequal_sk(token_start, "ALT")) {
                cur_col_type = 3;
              } else {
                continue;
              }
            } else if (token_slen == 2) {
              if (memequal_sk(token_start, "ID")) {
                cur_col_type = 1;
              } else if (memequal_sk(token_start, "CM")) {
                if (!read_cm) {
                  continue;
                }
                cur_col_type = 7;
              } else {
                continue;
              }
            } else {
              continue;
            }
          } else if (strequal_k(token_start, "QUAL", token_slen)) {
            if (!read_qual) {
              continue;
            }
            cur_col_type = 4;
          } else if (strequal_k(token_start, "INFO", token_slen)) {
            if (!read_info) {
              continue;
            }
            cur_col_type = 6;
          } else if (token_slen == 6) {
            if (memequal_sk(token_start, "FILTER")) {
              if (!read_filter) {
                continue;
              }
              cur_col_type = 5;
            } else if (memequal_sk(token_start, "FORMAT")) {
              break;
            } else {
              continue;
            }
          } else {
            continue;
          }
          col_skips[relevant_postchr_col_ct] = col_idx;
          col_types[relevant_postchr_col_ct++] = cur_col_type;
        }
        for (uint32_t rpc_col_idx = relevant_postchr_col_ct - 1; rpc_col_idx; --rpc_col_idx) {
          col_skips[rpc_col_idx] -= col_skips[rpc_col_idx - 1];
        }
        line_start = AdvPastDelim(token_end, '\n');
        ++pvar_line_idx;
      } else {
        col_skips[0] = 1;
        col_skips[1] = 1;
        col_skips[2] = 1;
        col_skips[3] = 1;
        col_types[0] = 1;
        if (!read_cm) {
          relevant_postchr_col_ct = 4;
          col_types[1] = 0;
          col_types[2] = 3;
          col_types[3] = 2;
          const char* sixth_col_start = NextTokenMult(line_start, 5);
          if (sixth_col_start) {
            // bugfix (25 Mar 2021)
            col_skips[1] = 2;
          }
        } else {
          relevant_postchr_col_ct = 5;
          col_skips[4] = 1;
          col_types[1] = 7;
          col_types[2] = 0;
          col_types[3] = 3;
          col_types[4] = 2;
        }
      }

      max_single_pos_ct = filesets_iter->max_single_pos_ct;
      const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
      char* cur_pos_readbuf;
      SamePosPvarRecord** same_pos_records;
      if (unlikely(bigstack_alloc_c(max_chr_blen, &ppmc.pmc.chr_buf) ||
                   bigstack_alloc_c(max_single_pos_blen + (sizeof(SamePosPvarRecord) + 1) * max_single_pos_ct, &cur_pos_readbuf) ||
                   BIGSTACK_ALLOC_X(SamePosPvarRecord*, max_single_pos_ct, &same_pos_records))) {
        goto PmergeConcat_ret_NOMEM;
      }
      if (varid_templatep) {
        varid_templatep->chr_output_name_buf = ppmc.pmc.chr_buf;
        if (varid_multi_templatep) {
          varid_multi_templatep->chr_output_name_buf = ppmc.pmc.chr_buf;
        }
        if (varid_multi_nonsnp_templatep) {
          varid_multi_nonsnp_templatep->chr_output_name_buf = ppmc.pmc.chr_buf;
        }
      }
      line_idx_body_start = pvar_line_idx;
      char* cur_pos_readbuf_iter = cur_pos_readbuf;
      uint32_t cur_single_pos_ct = 0;
      uint32_t prev_chr_idx = UINT32_MAX;
      int32_t prev_bp = 0;
      if (g_debug_on && (fileset_idx == 2)) {
        logprintf("Starting chr3.\n");
      }
      for (uint32_t read_variant_idx = 0; read_variant_idx != read_variant_ct; ++read_variant_idx) {
        if (unlikely(!TextGetUnsafe2(&pvar_txs, &line_start))) {
          reterr = TextStreamRawErrcode(&pvar_txs);
          goto PmergeConcat_ret_PVAR_TSTREAM_REWIND_FAIL_N;
        }
        char* chr_token_end = CurTokenEnd(line_start);
        const uint32_t chr_idx = GetChrCodeCounted(cip, chr_token_end - line_start, line_start);
        assert(chr_idx < UINT32_MAXM1);
        if (!IsSet(cip->chr_mask, chr_idx)) {
          if (pgfi.allele_idx_offsets) {
            // No need to actually parse ALT: we only need
            // (pgfi.allele_idx_offsets[x+1] - pgfi.allele_idx_offsets[x]) to
            // account for extra ALT alleles for variants we're properly
            // reading.  We may still process the *biallelic* portion of
            // variants we aren't explicitly reading (due to how LD-compression
            // works), so we defensively set the allele count to 2.
            pgfi.allele_idx_offsets[read_variant_idx + 1] = pgfi.allele_idx_offsets[read_variant_idx] + 2;
          }
          line_start = AdvPastDelim(chr_token_end, '\n');
          continue;
        }
        if (chr_idx != prev_chr_idx) {
          reterr = ConcatPvariantPos(prev_bp, cur_single_pos_ct, &ppmc, same_pos_records, &mr, &mw);
          if (unlikely(reterr)) {
            goto PmergeConcat_ret_N;
          }
          cur_pos_readbuf_iter = cur_pos_readbuf;
          cur_single_pos_ct = 0;
          char* chr_name_end = chrtoa(cip, chr_idx, ppmc.pmc.chr_buf);
          *chr_name_end = '\0';
          const uint32_t chr_slen = chr_name_end - ppmc.pmc.chr_buf;
          ppmc.pmc.chr_slen = chr_slen;
          if (varid_templatep) {
            const int32_t chr_slen_delta = chr_slen - varid_templatep->chr_slen;
            varid_templatep->chr_slen = chr_slen;
            varid_templatep->base_len += chr_slen_delta;
            if (varid_multi_templatep) {
              varid_multi_templatep->chr_slen = chr_slen;
              varid_multi_templatep->base_len += chr_slen_delta;
            }
            if (varid_multi_nonsnp_templatep) {
              varid_multi_nonsnp_templatep->chr_slen = chr_slen;
              varid_multi_nonsnp_templatep->base_len += chr_slen_delta;
            }
          }
          prev_chr_idx = chr_idx;
          prev_bp = -1;
        }
        char* token_ptrs[8];
        uint32_t token_slens[8];
        char* line_iter = TokenLex(chr_token_end, col_types, col_skips, relevant_postchr_col_ct, token_ptrs, token_slens);
        if (unlikely(!line_iter)) {
          goto PmergeConcat_ret_PVAR_REWIND_FAIL_N;
        }
        line_start = AdvPastDelim(line_iter, '\n');
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(token_ptrs[0], &cur_bp))) {
          goto PmergeConcat_ret_PVAR_REWIND_FAIL_N;
        }
        if (cur_bp < 0) {
          if (pgfi.allele_idx_offsets) {
            pgfi.allele_idx_offsets[read_variant_idx + 1] = pgfi.allele_idx_offsets[read_variant_idx] + 2;
          }
          continue;
        }
        if (cur_bp > prev_bp) {
          reterr = ConcatPvariantPos(prev_bp, cur_single_pos_ct, &ppmc, same_pos_records, &mr, &mw);
          if (unlikely(reterr)) {
            goto PmergeConcat_ret_N;
          }
          cur_pos_readbuf_iter = cur_pos_readbuf;
          cur_single_pos_ct = 0;
          prev_bp = cur_bp;
        }
        SamePosPvarRecord* cur_record = R_CAST(SamePosPvarRecord*, cur_pos_readbuf_iter);
        uint32_t* other_field_offsets = cur_record->other_field_offsets;
        cur_record->secondary_key = read_variant_idx;
        uint32_t variant_id_slen = token_slens[1];
        char* cur_variant_id_start = cur_record->variant_id;
        char* ref_start = token_ptrs[2];
        const uint32_t ref_slen = token_slens[2];
        char* alt_start = token_ptrs[3];
        const uint32_t alt_slen = token_slens[3];
        const uint32_t extra_alt_ct = CountByte(alt_start, ',', alt_slen);
        if (varid_templatep && ((!missing_varid_match_slen) || ((variant_id_slen == missing_varid_match_slen) && memequal(token_ptrs[1], missing_varid_match, missing_varid_match_slen)))) {
          const VaridTemplate* cur_varid_templatep = varid_templatep;
          if (extra_alt_ct && (varid_multi_templatep || varid_multi_nonsnp_templatep)) {
            if (varid_multi_templatep) {
              cur_varid_templatep = varid_multi_templatep;
            }
            if (varid_multi_nonsnp_templatep) {
              if ((ref_slen > 1) || (alt_slen != 2 * extra_alt_ct + 1)) {
                cur_varid_templatep = varid_multi_nonsnp_templatep;
              }
            }
          }
          // already scanned, shouldn't be possible for this to fail
          uint32_t discard;
          cur_pos_readbuf_iter = VaridTemplateWrite(cur_varid_templatep, ref_start, alt_start, cur_bp, ref_slen, extra_alt_ct, alt_slen, &discard, cur_variant_id_start);
          variant_id_slen = cur_pos_readbuf_iter - cur_variant_id_start;
        } else {
          cur_pos_readbuf_iter = memcpya(cur_variant_id_start, token_ptrs[1], variant_id_slen);
        }
        *cur_pos_readbuf_iter++ = '\0';
        other_field_offsets[0] = variant_id_slen + 1;

        if ((ref_start[0] != input_missing_geno_char) || (ref_slen != 1)) {
          cur_pos_readbuf_iter = memcpya(cur_pos_readbuf_iter, ref_start, ref_slen);
        } else {
          *cur_pos_readbuf_iter++ = '.';
        }
        *cur_pos_readbuf_iter++ = '\0';
        other_field_offsets[1] = cur_pos_readbuf_iter - cur_variant_id_start;

        if ((alt_start[0] != input_missing_geno_char) || (alt_slen != 1)) {
          cur_pos_readbuf_iter = memcpya(cur_pos_readbuf_iter, alt_start, alt_slen);
        } else {
          *cur_pos_readbuf_iter++ = '.';
        }
        *cur_pos_readbuf_iter++ = '\0';
        other_field_offsets[2] = cur_pos_readbuf_iter - cur_variant_id_start;
        cur_record->allele_ct = 2 + extra_alt_ct;
        if (pgfi.allele_idx_offsets) {
          pgfi.allele_idx_offsets[read_variant_idx + 1] = pgfi.allele_idx_offsets[read_variant_idx] + 2 + extra_alt_ct;
        }

        uint32_t pgen_pr_status = pgen_pr_status_base;
        if (pgfi.nonref_flags) {
          pgen_pr_status |= IsSet(pgfi.nonref_flags, read_variant_idx);
        }
        cur_record->pgen_pr_status = pgen_pr_status;

        if (read_qual) {
          cur_pos_readbuf_iter = memcpyax(cur_pos_readbuf_iter, token_ptrs[4], token_slens[4], '\0');
        }
        other_field_offsets[3] = cur_pos_readbuf_iter - cur_variant_id_start;

        if (read_filter) {
          cur_pos_readbuf_iter = memcpyax(cur_pos_readbuf_iter, token_ptrs[5], token_slens[5], '\0');
        }
        other_field_offsets[4] = cur_pos_readbuf_iter - cur_variant_id_start;

        if (read_info) {
          cur_pos_readbuf_iter = memcpyax(cur_pos_readbuf_iter, token_ptrs[6], token_slens[6], '\0');
        }
        other_field_offsets[5] = cur_pos_readbuf_iter - cur_variant_id_start;

        if (read_cm) {
          cur_pos_readbuf_iter = memcpya(cur_pos_readbuf_iter, token_ptrs[7], token_slens[7]);
        }
        *cur_pos_readbuf_iter++ = '\0';
        // could align up to 8-byte boundary?
        assert(cur_pos_readbuf_iter <= R_CAST(char*, same_pos_records));

        same_pos_records[cur_single_pos_ct] = cur_record;
        ++cur_single_pos_ct;
      }
      reterr = ConcatPvariantPos(prev_bp, cur_single_pos_ct, &ppmc, same_pos_records, &mr, &mw);
      if (unlikely(reterr)) {
        goto PmergeConcat_ret_N;
      }
      // bugfix (14 Apr 2021): forgot to close .pgen
      if (unlikely(CleanupTextStream2(read_pvar_fname, &pvar_txs, &reterr) ||
                   CleanupPgr2(read_pgen_fname, &mr.pgr, &reterr) ||
                   CleanupPgfi2(read_pgen_fname, &pgfi, &reterr))) {
        goto PmergeConcat_ret_N;
      }
    }
    reterr = SpgwFinish(&mw.spgw);
    if (unlikely(reterr)) {
      goto PmergeConcat_ret_1;
    }
    if (unlikely(CswriteCloseNull(&ppmc.pmc.css, ppmc.pmc.cswritep))) {
      goto PmergeConcat_ret_WRITE_FAIL_N;
    }
    fputs("\rConcatenating... ", stdout);
    logprintf("%" PRIuPTR "/%" PRIuPTR " variant%s complete.\n", write_variant_ct, write_variant_ct, (write_variant_ct == 1)? "" : "s");
    *outname_end = '\0';
    logprintfww("Results written to %s.pgen + %s.pvar%s .\n", outname, outname, pvar_zst? ".zst" : "");
  }
  while (0) {
  PmergeConcat_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  PmergeConcat_ret_PVAR_TSTREAM_REWIND_FAIL_N:
    logputs("\n");
    TextStreamErrPrintRewind(read_pvar_fname, &pvar_txs, &reterr);
    break;
  PmergeConcat_ret_PVAR_REWIND_FAIL_N:
    logputs("\n");
    logerrprintfww(kErrprintfRewind, read_pvar_fname);
    reterr = kPglRetRewindFail;
    break;
  PmergeConcat_ret_PGEN_REWIND_FAIL_N:
    logputs("\n");
    logerrprintfww(kErrprintfRewind, read_pgen_fname);
    reterr = kPglRetRewindFail;
    break;
  PmergeConcat_ret_WRITE_FAIL_N:
    logputs("\n");
    reterr = kPglRetWriteFail;
    break;
  PmergeConcat_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  PmergeConcat_ret_N:
    logputs("\n");
    break;
  }
 PmergeConcat_ret_1:
  CleanupPvariantPosMergeContext(&ppmc);
  CleanupTextStream2(read_pvar_fname, &pvar_txs, &reterr);
  CleanupSpgw(&mw.spgw, &reterr);
  CleanupPgr2(read_pgen_fname, &mr.pgr, &reterr);
  CleanupPgfi2(read_pgen_fname, &pgfi, &reterr);
  return reterr;
}

// Performs one pass of a possibly-multipass merge.
// *input_filesets_ptr is advanced as files are processed, and temporary files
// are deleted while this happens.
// next_filesets is left as nullptr on the last pass.
PglErr PmergePass(__attribute__((unused)) const PmergeInfo* pmip, __attribute__((unused)) const SampleIdInfo* siip, __attribute__((unused)) const ChrInfo* cip, __attribute__((unused)) const char* const* info_keys, __attribute__((unused)) const uint32_t* info_keys_htable, __attribute__((unused)) uint32_t pass_idx, __attribute__((unused)) uint32_t sample_ct, __attribute__((unused)) FamCol fam_cols, uintptr_t fileset_ct, __attribute__((unused)) uint32_t psam_linebuf_capacity, __attribute__((unused)) uint32_t info_key_ct, __attribute__((unused)) uint32_t info_keys_htable_size, __attribute__((unused)) uint32_t info_conflict_present, __attribute__((unused)) uint32_t max_thread_ct, __attribute__((unused)) SortMode sort_vars_mode, __attribute__((unused)) char* outname, __attribute__((unused)) char* outname_end, __attribute__((unused)) PmergeInputFilesetLl** input_filesets_ptr, __attribute__((unused)) PmergeInputFilesetLl** next_filesets_ptr, __attribute__((unused)) uint32_t* next_fileset_ctp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    // 1. Determine how many input filesets from the front of the list can be
    //    processed in the current loop iteration.  It's bounded by
    //    kMaxOpenFiles and workspace memory.  The workspace memory computation
    //    is annoying.
    //    - If there's only one input fileset left, just move that entry to the
    //      tail of next_filesets, and exit.
    //    - If all input filesets can be merged at once, write to the final
    //      output filenames, and leave next_filesets as nullptr.
    //    - If we don't obviously have enough memory to even merge the next two
    //      filesets... could either error out, or speculatively continue with
    //      two depending on how much we might overestimate memory usage.
    // 2. Incremental-merge them.
    // 3. Return to step 1 if input fileset(s) remain.
    uintptr_t input_filesets_remaining = fileset_ct;
    PmergeInputFilesetLl** next_filesets_end_ptr = next_filesets_ptr;
    do {
      logerrputs("Error: Non-concatenating --pmerge[-list] is under development.\n");
      reterr = kPglRetNotYetSupported;
      goto PmergePass_ret_1;
    } while (input_filesets_remaining > 1);
    if (input_filesets_remaining == 1) {
      *next_filesets_end_ptr = *input_filesets_ptr;
      *input_filesets_ptr = nullptr;
    }
  }
  while (0) {
  }
 PmergePass_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr Pmerge(const PmergeInfo* pmip, const char* sample_sort_fname, const char* missing_catname, const char* varid_template_str, const char* varid_multi_template_str, const char* varid_multi_nonsnp_template_str, const char* missing_varid_match, MiscFlags misc_flags, LoadFilterLogFlags load_filter_log_merge_flags, SortMode sample_sort_mode, FamCol fam_cols, int32_t missing_pheno, uint32_t new_variant_id_max_allele_slen, char input_missing_geno_char, uint32_t max_thread_ct, SortMode sort_vars_mode, char* pgenname, char* psamname, char* pvarname, char* outname, char* outname_end, ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;

  // nodes and filenames are heap-allocated, not just first_varid/last_varid
  PmergeInputFilesetLl* input_filesets = nullptr;
  PmergeInputFilesetLl* next_filesets = nullptr;

  PglErr reterr = kPglRetSuccess;
  {
    // 1. Construct/load fileset list.
    // 2. Merge .psam files.
    // 3. Global .pvar scan, to determine merged header, chromosome sort order,
    //    and track first/last position in each fileset.
    // 4. If filesets cover disjoint positions, handle this as a concatenation
    //    job (or error out on --variant-inner-join).
    // 5. Otherwise, perform general-purpose incremental merge.
    uintptr_t fileset_ct = 2;
    {
      PmergeInputFilesetLl** filesets_endp = &input_filesets;
      if (pgenname[0]) {
        PmergeInputFilesetLl* cur_entry = AllocFilesetLlEntry(&filesets_endp);
        if (unlikely(!cur_entry)) {
          goto Pmerge_ret_NOMEM;
        }
        cur_entry->is_temporary = 0;
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
        cur_entry->first_varid = nullptr;
        cur_entry->last_varid = nullptr;
      }
      if (!pmip->list_fname) {
        PmergeInputFilesetLl* cur_entry = AllocFilesetLlEntry(&filesets_endp);
        if (unlikely(!cur_entry)) {
          goto Pmerge_ret_NOMEM;
        }
        cur_entry->is_temporary = 0;
        cur_entry->pgen_fname = pmip->pgen_fname;
        cur_entry->pvar_fname = pmip->pvar_fname;
        cur_entry->psam_fname = pmip->psam_fname;
        cur_entry->pgen_locked_fname = nullptr;
        cur_entry->first_varid = nullptr;
        cur_entry->last_varid = nullptr;
      } else {
        reterr = LoadPmergeList(pmip->list_fname, pmip->list_base_dir, pmip->list_mode, pgenname[0] != '\0', &filesets_endp, &fileset_ct);
        if (unlikely(reterr)) {
          goto Pmerge_ret_1;
        }
      }
    }

    SampleIdInfo sii;
    uint32_t sample_ct = 0;
    uint32_t psam_linebuf_capacity = 0;
    reterr = MergePsams(pmip, sample_sort_fname, missing_catname, misc_flags, sample_sort_mode, fam_cols, missing_pheno, max_thread_ct, outname, outname_end, input_filesets, &sii, &sample_ct, &psam_linebuf_capacity);
    if (unlikely(reterr)) {
      goto Pmerge_ret_1;
    }
    reterr = ScanPgenHeaders(!!(pmip->list_fname), misc_flags, input_filesets);
    if (unlikely(reterr)) {
      goto Pmerge_ret_1;
    }

    VaridTemplate* varid_templatep = nullptr;
    VaridTemplate* varid_multi_templatep = nullptr;
    VaridTemplate* varid_multi_nonsnp_templatep = nullptr;
    uint32_t missing_varid_match_slen = 0;
    if (varid_template_str) {
      if (unlikely(VaridInitAll(g_bigstack_end, varid_template_str, varid_multi_template_str, varid_multi_nonsnp_template_str, misc_flags, new_variant_id_max_allele_slen, &g_bigstack_base, &missing_varid_match, nullptr, &varid_templatep, &varid_multi_templatep, &varid_multi_nonsnp_templatep, nullptr, &missing_varid_match_slen))) {
        goto Pmerge_ret_NOMEM;
      }
    }

    uint32_t info_has_g_key = 0;
    const char* const* info_keys = nullptr;
    uint32_t* info_keys_htable = nullptr;
    uint32_t info_key_ct = 0;
    uint32_t info_keys_htable_size = 0;
    uint32_t info_conflict_present;
    reterr = ScanPvarsAndMergeHeader(pmip, missing_varid_match, misc_flags, load_filter_log_merge_flags, missing_varid_match_slen, input_missing_geno_char, max_thread_ct, sort_vars_mode, varid_templatep, varid_multi_templatep, varid_multi_nonsnp_templatep, outname, outname_end, &input_filesets, cip, &fileset_ct, &info_has_g_key, &info_keys, &info_key_ct, &info_keys_htable, &info_keys_htable_size, &info_conflict_present);
    if (unlikely(reterr)) {
      goto Pmerge_ret_1;
    }
    uint32_t is_concat_job = 0;
    if (!(pmip->flags & (kfPmergeVariantInnerJoin | kfPmergeSids))) {
      reterr = DetectConcatJob(cip->chr_idx_to_foidx, fileset_ct, sort_vars_mode, &input_filesets, &is_concat_job);
      if (unlikely(reterr)) {
        goto Pmerge_ret_1;
      }
    }
    if (is_concat_job) {
      reterr = PmergeConcat(pmip, &sii, cip, input_filesets, missing_varid_match, info_keys, info_keys_htable, sample_ct, fam_cols, fileset_ct, psam_linebuf_capacity, missing_varid_match_slen, info_key_ct, info_keys_htable_size, info_conflict_present, input_missing_geno_char, max_thread_ct, sort_vars_mode, varid_templatep, varid_multi_templatep, varid_multi_nonsnp_templatep, outname, outname_end);
    } else {
      for (uint32_t pass_idx = 1; ; ++pass_idx) {
        uint32_t next_fileset_ct;
        reterr = PmergePass(pmip, &sii, cip, info_keys, info_keys_htable, pass_idx, sample_ct, fam_cols, fileset_ct, psam_linebuf_capacity, info_key_ct, info_keys_htable_size, info_conflict_present, max_thread_ct, sort_vars_mode, outname, outname_end, &input_filesets, &next_filesets, &next_fileset_ct);
        if (unlikely(reterr)) {
          goto Pmerge_ret_1;
        }
        if (!next_filesets) {
          break;
        }
        assert(input_filesets == nullptr);
        fileset_ct = next_fileset_ct;
        input_filesets = next_filesets;
        next_filesets = nullptr;
      }
    }
    const uint32_t outname_slen = outname_end - outname;
    memcpy(pgenname, outname, outname_slen);
    strcpy_k(&(pgenname[outname_slen]), ".pgen");
    memcpy(pvarname, outname, outname_slen);
    const uint32_t output_zst = (pmip->flags / kfPmergeOutputVzs) & 1;
    if (output_zst) {
      strcpy_k(&(pvarname[outname_slen]), ".pvar.zst");
    } else {
      strcpy_k(&(pvarname[outname_slen]), ".pvar");
    }
    memcpy(psamname, outname, outname_slen);
    strcpy_k(&(psamname[outname_slen]), ".psam");
    if ((pmip->flags & kfPmergeVariantInnerJoin) || pmip->max_allele_ct) {
      // Some input variants that pass the chromosome filter might still be
      // excluded from the merged dataset.
      ForgetExtraChrNames(1, cip);
    }
  }
  while (0) {
  Pmerge_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 Pmerge_ret_1:
  CleanupFilesetLl(input_filesets, &reterr);
  CleanupFilesetLl(next_filesets, &reterr);
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
PglErr PgenDiff(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const PgenDiffInfo* pdip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t max_allele_ct1, uint32_t max_allele_slen, char input_missing_geno_char, uint32_t max_thread_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
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
      reterr = SortedXidboxInitAllocEnd(orig_sample_include, siip, orig_sample_ct, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
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
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID \"%s\" in --pgen-diff .psam file.\n", idbuf);
          goto PgenDiff_ret_MALFORMED_INPUT_WW;
        }
        SetBit(sample_uidx, sample_include);
        sample1_uidx_to_2[sample_uidx] = psam_line_idx - first_payload_line_idx;
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
      if (raw_sample_ct2_ul > kPglMaxSampleCt) {
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
    reterr = PgfiInitPhase1(pdip->pgen_fname, nullptr, raw_variant_ct2, raw_sample_ct2, &header_ctrl, &pgfi2, &cur_alloc_cacheline_ct, g_logbuf);
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
    const uintptr_t* nonref_flags1 = pgfip->nonref_flags;
    const uint32_t all_nonref1 = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags1);
    const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags1, flags / kfPgenDiffColMaybeprovref, raw_variant_ct, all_nonref1);
    if (provref_col) {
      cswritep = strcpya_k(cswritep, "PROVISIONAL_REF?\t");
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
    uint32_t next_print_variant_uidx2 = (raw_variant_ct2 + 99) / 100;
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
        next_print_variant_uidx2 = (pct * S_CAST(uint64_t, raw_variant_ct2) + 99) / 100;
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
          variant_uidx_end = ExpsearchU32(variant_bps, search_start, chrom_end_variant_uidx, cur_included_bp + 1);
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
          variant_uidx_start = ExpsearchU32(variant_bps, variant_uidx_end, chrom_end_variant_uidx, new_bp_u32);
          variant_uidx_start = AdvBoundedTo1Bit(variant_include, variant_uidx_start, chrom_end_variant_uidx);
          if (variant_uidx_start == chrom_end_variant_uidx) {
            variant_uidx_end = variant_uidx_start;
            cur_included_bp = UINT32_MAX;
            continue;
          }
          cur_included_bp = variant_bps[variant_uidx_start];
          const uint32_t search_start = variant_uidx_start + 1;
          variant_uidx_end = ExpsearchU32(variant_bps, search_start, chrom_end_variant_uidx, cur_included_bp + 1);
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
          PgenErrPrintNV(reterr, variant_uidx);
          goto PgenDiff_ret_1;
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
          PgenErrPrintNV(reterr, variant_uidx);
          goto PgenDiff_ret_1;
        }
        if (cur_allele_ct2 == 2) {
          reterr = PgrGetD(sample_include2, pssi2, sample_ct, variant_uidx2, &simple_pgr2, pgv2.genovec, pgv2.dosage_present, pgv2.dosage_main, &pgv2.dosage_ct);
        } else {
          reterr = PgrGetMD(sample_include2, pssi2, sample_ct, variant_uidx2, &simple_pgr2, &pgv2);
        }
      }
      if (unlikely(reterr)) {
        PgenErrPrintNV(reterr, variant_uidx2);
        goto PgenDiff_ret_1;
      }

      uintptr_t* genovec1 = pgv1.genovec;
      uintptr_t* genovec2 = pgv2.genovec;
      uint32_t merged_allele_ct = 0;
      uint32_t merged_provref = 0;
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
        // present, it's in a biallelic variant.  bugfix (3 Feb 2021): forgot
        // that biallelic variants can have *both* allele codes missing
        // (consider a .ped-derived variant with only missing calls).
        // We verify below that the missing allele code does not appear in the
        // actual genotype calls.
        uint32_t missing1_state = 0;
        for (uint32_t allele1_idx = 0; allele1_idx != cur_allele_ct1; ++allele1_idx) {
          const char* cur_allele1 = cur_allele1s[allele1_idx];
          if (strequal_k_unsafe(cur_allele1, ".")) {
            missing1_state |= allele1_idx + 1;
          }
        }
        ZeroTrailingNyps(sample_ct, genovec1);
        ZeroTrailingNyps(sample_ct, genovec2);
        if (missing1_state) {
          assert(cur_allele_ct1 == 2);
          STD_ARRAY_DECL(uint32_t, 4, genocounts);
          GenoarrCountFreqsUnsafe(genovec1, sample_ct, genocounts);
          if (unlikely(genocounts[1] || ((missing1_state & 1) && genocounts[0]) || ((missing1_state & 2) && genocounts[2]) || pgv1.dosage_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Missing allele for variant '%s' at position %s:%u is present in the .pgen.\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
            goto PgenDiff_ret_INCONSISTENT_INPUT_WW_N;
          }
          if (missing1_state & 1) {
            remap1[0] = kMissingAlleleCode;
          }
          if (missing1_state & 2) {
            remap1[1] = kMissingAlleleCode;
          }
        }
        uint32_t missing2_state = 0;
        for (uint32_t allele2_idx = 0; allele2_idx != cur_allele_ct2; ++allele2_idx) {
          const char* cur_allele2 = cur_allele2s[allele2_idx];
          if (((cur_allele2[0] == '.') || (cur_allele2[0] == input_missing_geno_char)) && (cur_allele2[1] == '\0')) {
            if (unlikely(cur_allele_ct2 > 2)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Missing allele in multiallelic variant on line %" PRIuPTR " of --pgen-diff .pvar file.\n", pvar_line_idx);
              goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
            }
            missing2_state |= allele2_idx + 1;
          }
        }
        if (missing2_state) {
          STD_ARRAY_DECL(uint32_t, 4, genocounts);
          GenoarrCountFreqsUnsafe(genovec2, sample_ct, genocounts);
          if (unlikely(genocounts[1] || ((missing2_state & 1) && genocounts[0]) || ((missing2_state & 2) && genocounts[2]) || ((cur_allele_ct2 > 2) && (pgv2.patch_01_ct || pgv2.patch_10_ct)) || pgv2.dosage_ct)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Missing allele on line %" PRIuPTR " of --pgen-diff .pvar file is present in the .pgen.\n", pvar_line_idx);
            goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
          }
          if (missing2_state & 1) {
            remap2[0] = kMissingAlleleCode;
          }
          if (missing2_state & 2) {
            remap2[1] = kMissingAlleleCode;
          }
        }
        const uint32_t provisional_ref1 = all_nonref1 || (nonref_flags1 && IsSet(nonref_flags1, variant_uidx));
        if ((!provisional_ref1) && (missing1_state & 1)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Missing REF allele for variant '%s' at position %s:%u is not flagged as provisional.\n", variant_ids[variant_uidx], chr_buf, cur_included_bp);
          goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
        }
        const uint32_t provisional_ref2 = all_nonref2 || (nonref_flags2 && IsSet(nonref_flags2, variant_uidx2));
        if ((!provisional_ref2) && (missing2_state & 1)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Missing REF allele on line %" PRIuPTR " of --pgen-diff .pvar file is not flagged as provisional.\n", pvar_line_idx);
          goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
        }
        if ((missing1_state == 3) && (missing2_state == 3)) {
          // Both variants are all-missing, no differences possible.
          continue;
        }
        // Initialize these to the index of the first nonmissing allele.
        uint32_t allele1_alt_start = (missing1_state & 1) + (missing1_state == 3);
        uint32_t allele2_alt_start = (missing2_state & 1) + (missing2_state == 3);
        merged_provref = provisional_ref1;
        if ((!provisional_ref1) || ((!(missing1_state & 1)) && provisional_ref2)) {
          // - If REF1 is not provisional, it's always the merged REF allele.
          // - If REF1 is provisional, it's still the merged REF allele if REF2
          //   is also provisional and REF1 is nonmissing.
          if (unlikely((!provisional_ref2) && (!strequal_overread(cur_allele1s[0], cur_allele2s[0])))) {
            snprintf(g_logbuf, kLogbufSize, "Error: REF allele on line %" PRIuPTR " of --pgen-diff .pvar file conflicts with loaded REF allele, and neither are flagged as provisional.\n", pvar_line_idx);
            goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
          }
          merged_alleles[0] = cur_allele1s[0];
          remap1[0] = 0;
          allele1_alt_start = 1 + ((missing1_state >> 1) & 1);
        } else {
          // REF1 is not the merged REF allele.
          if (!(missing2_state & 1)) {
            // If REF2 can be the merged REF allele (i.e. it isn't missing),
            // make it so.
            merged_alleles[0] = cur_allele2s[0];
            remap2[0] = 0;
            allele2_alt_start = 1 + ((missing2_state >> 1) & 1);
            merged_provref = provisional_ref2;
          } else {
            // Both REFs missing.  Treat main dataset ALT1 as REF if it isn't
            // missing.
            if (!(missing1_state & 2)) {
              merged_alleles[0] = cur_allele1s[1];
              remap1[1] = 0;
              allele1_alt_start = 2;
            } else {
              // Other dataset's ALT1 guaranteed to be nonmissing if we get
              // here.
              merged_alleles[0] = cur_allele2s[1];
              remap2[1] = 0;
              allele2_alt_start = 2;
            }
          }
        }
        const uint32_t cur_max_merged_alleles = 1 + cur_allele_ct1 + cur_allele_ct2 - allele1_alt_start - allele2_alt_start;
        const uint32_t cur_htable_size = (cur_max_merged_alleles * 9) / 2;
        SetAllU32Arr(cur_htable_size, merged_alleles_htable);
        // See e.g. PopulateStrboxHtable().
        uint32_t hashval = Hashceil(merged_alleles[0], strlen(merged_alleles[0]), cur_htable_size);
        merged_alleles_htable[hashval] = 0;
        ++merged_allele_ct;

        for (uint32_t allele1_idx = allele1_alt_start; allele1_idx != cur_allele_ct1; ++allele1_idx) {
          if ((allele1_idx < 2) && ((missing1_state >> allele1_idx) & 1)) {
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
          if ((allele2_idx < 2) && ((missing2_state >> allele2_idx) & 1)) {
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
        if (unlikely(DuplicateAllelePresent(remap2, cur_allele_ct2, merged_allele_ctl, remap_seen))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate allele code on line %" PRIuPTR " of --pgen-diff .pvar file.\n", pvar_line_idx);
          goto PgenDiff_ret_MALFORMED_INPUT_WW_N;
        }
        if (merged_allele_ct == 1) {
          if (!include_missing) {
            // Only possible difference is missing vs. not-missing.
            continue;
          }
          assert(missing1_state);
          merged_alleles[1] = cur_allele1s[ctzu32(missing1_state)];
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
          // unknown-sex prohibited, so sex_male is fine
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
        if (provref_col) {
          *cswritep++ = merged_provref? 'Y' : 'N';
          *cswritep++ = '\t';
        }
        const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
        if (!fid_col) {
          cur_sample_id = AdvPastDelim(cur_sample_id, '\t');
        }
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
