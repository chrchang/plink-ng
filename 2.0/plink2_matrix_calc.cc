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

#include "plink2_matrix_calc.h"

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "plink2_compress_stream.h"
#include "plink2_decompress.h"
#include "plink2_matrix.h"
#include "plink2_random.h"

#ifdef USE_CUDA
#  include "cuda/plink2_matrix_cuda.h"
#endif

#include <unistd.h>  // unlink()

#ifdef __cplusplus
namespace plink2 {
#endif

void InitScore(ScoreInfo* score_info_ptr) {
  score_info_ptr->flags = kfScore0;
  score_info_ptr->varid_col_p1 = 1;
  score_info_ptr->allele_col_p1 = 0;  // defensive
  score_info_ptr->input_fname = nullptr;
  InitRangeList(&(score_info_ptr->input_col_idx_range_list));

  score_info_ptr->qsr_range_fname = nullptr;
  score_info_ptr->qsr_data_fname = nullptr;
  score_info_ptr->qsr_varid_col_p1 = 1;
  score_info_ptr->qsr_val_col_p1 = 0;  // defensive
}

void CleanupScore(ScoreInfo* score_info_ptr) {
  free_cond(score_info_ptr->input_fname);
  CleanupRangeList(&(score_info_ptr->input_col_idx_range_list));

  free_cond(score_info_ptr->qsr_range_fname);
  free_cond(score_info_ptr->qsr_data_fname);
}

void InitPhenoSvd(PhenoSvdInfo* pheno_svd_info_ptr) {
  pheno_svd_info_ptr->flags = kfPhenoSvd0;
  pheno_svd_info_ptr->ct = 0;
  pheno_svd_info_ptr->min_variance_explained = 0.0;
}


// Cost function for thread load-balancing: (max - min) * ((max + min)/2)
// i.e. we don't waste any time processing entries in the irrelevant
// upper-right triangle
// set align to 1 for no alignment
void TriangleFill(uint32_t ct, uint32_t piece_ct, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t start, uint32_t align, uint32_t* target_arr) {
  int32_t modif = 1 - start * 2;
  int64_t cur_prod_x2;
  int32_t lbound;
  int32_t ubound;
  uint32_t uii;
  uint32_t align_m1;
  ParallelBounds(ct, start, parallel_idx, parallel_tot, &lbound, &ubound);
  // x(x+1)/2 is divisible by y iff (x % (2y)) is 0 or (2y - 1).
  align *= 2;
  align_m1 = align - 1;
  target_arr[0] = lbound;
  target_arr[piece_ct] = ubound;
  cur_prod_x2 = S_CAST(int64_t, lbound) * (lbound + modif);
  const int64_t ct_tr = (S_CAST(int64_t, ubound) * (ubound + modif) - cur_prod_x2) / piece_ct;
  for (uint32_t piece_idx = 1; piece_idx != piece_ct; ++piece_idx) {
    cur_prod_x2 += ct_tr;
    lbound = TriangleDivide(cur_prod_x2, modif);
    uii = (lbound - S_CAST(int32_t, start)) & align_m1;
    if ((uii) && (uii != align_m1)) {
      lbound = start + ((lbound - S_CAST(int32_t, start)) | align_m1);
    }
    // lack of this check caused a nasty bug earlier
    if (S_CAST(uint32_t, lbound) > ct) {
      lbound = ct;
    }
    target_arr[piece_idx] = lbound;
  }
}

// Cost function for thread load-balancing: (max - min) * max
// i.e. when we can't avoid processing entries in the irrelevant upper-right
// triangle
void TriangleFill2(uint32_t ct, uint32_t piece_ct, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t start, uint32_t* target_arr) {
  int32_t lbound_s;
  int32_t ubound_s;
  ParallelBounds(ct, start, parallel_idx, parallel_tot, &lbound_s, &ubound_s);
  uint32_t lbound = lbound_s;
  const uint32_t ubound = ubound_s;
  target_arr[0] = lbound;
  target_arr[piece_ct] = ubound;

  uint32_t remaining_row_ct = ubound - lbound;
  uint32_t remaining_piece_ct = piece_ct;

  for (uint32_t piece_idx = 1; piece_idx != piece_ct; ++piece_idx) {
    // Start by assigning an equal-row-count piece, rounding up (since later
    // pieces are wider).  Then compare current-piece cost with a lower bound
    // of average remaining-piece cost, and append rows until current-piece
    // cost exceeds that average.
    uint32_t candidate_piece_size = (remaining_row_ct + remaining_piece_ct - 1) / remaining_piece_ct;
    remaining_row_ct -= candidate_piece_size;
    remaining_piece_ct -= 1;
    uint32_t candidate_boundary = lbound + candidate_piece_size;

    // Consider lbound == 0, ubound == 8, piece_ct == 3.  The ideal solution
    // would be:
    //
    // row 0: 0 0 0 0
    // row 1: 0 0 0 0
    // row 2: 0 0 0 0
    // row 3: 0 0 0 0
    // row 4: 1 1 1 1 1 1
    // row 5: 1 1 1 1 1 1
    // row 6: 2 2 2 2 2 2 2 2
    // row 7: 2 2 2 2 2 2 2 2
    //
    // The initial candidate_piece_size is 3.
    //
    // With no waste, the remaining cost would be (8-3) * (8+3+1)/2 = 30.
    // The true remaining cost is at least 34, because there are 5 remaining
    // rows, the upper-right-triangle waste is minimized if they're split as
    // evenly as possible between the remaining threads, and that {3, 2} split
    // results in waste of 3 for the 3-row piece and 1 for the 2-row piece.
    //
    // So current cost is 9, average remaining cost is 17.  If we incremented
    // candidate_piece_size, the current cost would become 16, and the
    // average remaining cost would drop to ~15 (actually 14, but we use a
    // slight overestimate in the current comparison); the latter pair of
    // numbers is closer, so we perform the increment.
    //
    // On the next iteration, we clearly don't want to increment
    // candidate_piece_size any more.
    //
    //
    // Then, when determining the boundary between threads 1 and 2,
    // candidate_piece_size starts at 2.  Current cost is 6*2=12, average
    // remaining cost is 16.  If we incremented, candidate_piece_size, the
    // current cost would become 21, and the average remaining cost would drop
    // to ~8.  |16-12| is smaller, so thread 1 gets just 2 rows.


    // If the remaining rows were split as evenly as possible, what's the
    // smaller chunk size?
    uint32_t remaining_piece_lower_size = remaining_row_ct / remaining_piece_ct;
    // What's the waste associated with that smaller chunk size?
    uint64_t remaining_piece_lower_waste = (S_CAST(uint64_t, remaining_piece_lower_size) * (remaining_piece_lower_size - 1)) >> 1;
    // How many threads would be assigned the higher chunk size?
    uint32_t remaining_piece_higher_ct = remaining_row_ct - remaining_piece_lower_size * remaining_piece_ct;
    // Note that the waste associated with the larger chunk size is
    // (remaining_piece_lower_waste + remaining_piece_lower_size).
    const uint64_t min_waste = remaining_piece_lower_waste * remaining_piece_ct + remaining_piece_higher_ct * remaining_piece_lower_size;
    uint64_t remaining_cost = ((remaining_row_ct * (ubound + 1LLU + candidate_boundary)) >> 1) + min_waste;
    while (1) {
      const uint64_t candidate_piece_cost = S_CAST(uint64_t, candidate_piece_size) * candidate_boundary;
      const uint64_t next_candidate_piece_cost = candidate_piece_cost + candidate_piece_size + candidate_boundary + 1;
      // printf("lbound: %u  ubound: %u  candidate_piece_cost: %" PRIu64 "  remaining_piece_ct: %u  remaining_cost: %" PRIu64 "\n", lbound, ubound, candidate_piece_cost, remaining_piece_ct, remaining_cost);
      if ((candidate_piece_cost + next_candidate_piece_cost) * remaining_piece_ct > 2 * (remaining_cost - candidate_boundary)) {
        break;
      }
      // There are circumstances where we can make larger jumps here, but the
      // overall computational cost here is low enough; let's keep this
      // simpler.
      ++candidate_piece_size;
      --remaining_row_ct;
      ++candidate_boundary;
      // Incremental update of non-waste portion of remaining_cost.
      remaining_cost -= candidate_boundary;
      // Incremental update of waste portion of remaining_cost.
      if (remaining_piece_higher_ct) {
        --remaining_piece_higher_ct;
        remaining_cost -= remaining_piece_lower_size;
      } else {
        --remaining_piece_lower_size;
        remaining_piece_lower_waste -= remaining_piece_lower_size;
        remaining_cost -= remaining_piece_lower_size;
        remaining_piece_higher_ct = remaining_piece_ct - 1;
      }
    }
    lbound += candidate_piece_size;
    target_arr[piece_idx] = lbound;
  }
}

// Returns 0 if cells_avail is insufficient.
uint32_t CountTrianglePasses(uintptr_t start_idx, uintptr_t end_idx, uintptr_t is_no_diag, uintptr_t cells_avail) {
  start_idx -= is_no_diag;
  end_idx -= is_no_diag;
  if (cells_avail < end_idx) {
    return 0;
  }
  cells_avail *= 2;  // don't want to worry about /2 in triangular numbers
  const uint64_t end_tri = S_CAST(uint64_t, end_idx) * (end_idx + 1);
  uint64_t start_tri = S_CAST(uint64_t, start_idx) * (start_idx + 1);
  for (uint32_t pass_ct = 1; ; ++pass_ct) {
    const uint64_t delta_tri = end_tri - start_tri;
    if (delta_tri <= cells_avail) {
      return pass_ct;
    }
    const uint64_t next_target = start_tri + cells_avail;
    start_idx = S_CAST(int64_t, sqrt(u63tod(next_target)));
    start_tri = S_CAST(uint64_t, start_idx) * (start_idx + 1);
    if (start_tri > next_target) {
      --start_idx;
      start_tri = S_CAST(uint64_t, start_idx) * (start_idx + 1);
      assert(start_tri <= next_target);
    }
  }
}

uint64_t NextTrianglePass(uintptr_t start_idx, uintptr_t grand_end_idx, uintptr_t is_no_diag, uintptr_t cells_avail) {
  cells_avail *= 2;
  start_idx -= is_no_diag;
  grand_end_idx -= is_no_diag;
  const uint64_t end_tri = S_CAST(uint64_t, grand_end_idx) * (grand_end_idx + 1);
  uint64_t start_tri = S_CAST(uint64_t, start_idx) * (start_idx + 1);
  const uint64_t delta_tri = end_tri - start_tri;
  if (delta_tri <= cells_avail) {
    return grand_end_idx + is_no_diag;
  }
  const uint64_t next_target = start_tri + cells_avail;
  start_idx = S_CAST(int64_t, sqrt(u63tod(next_target)));
  start_tri = S_CAST(uint64_t, start_idx) * (start_idx + 1);
  return start_idx + is_no_diag - (start_tri > next_target);
}

void TriangleLoadBalance(uint32_t piece_ct, uintptr_t start_idx, uintptr_t end_idx, uint32_t is_no_diag, uint32_t* target_arr) {
  target_arr[0] = start_idx;
  target_arr[piece_ct] = end_idx;
  start_idx -= is_no_diag;
  end_idx -= is_no_diag;
  const uint64_t end_tri = S_CAST(uint64_t, end_idx) * (end_idx + 1);
  uint64_t cur_target = S_CAST(uint64_t, start_idx) * (start_idx + 1);
  const uint64_t std_size = (end_tri - cur_target) / piece_ct;
  for (uint32_t piece_idx = 1; piece_idx != piece_ct; ++piece_idx) {
    // don't use cur_target = start_tri + (piece_idx * delta_tri) / piece_ct
    // because of potential overflow
    cur_target += std_size;
    start_idx = S_CAST(int64_t, sqrt(u63tod(cur_target)));
    const uint64_t start_tri = S_CAST(uint64_t, start_idx) * (start_idx + 1);
    if (start_tri > cur_target) {
      --start_idx;
    }
    target_arr[piece_idx] = start_idx + is_no_diag;
  }
}

PglErr KinshipPruneDestructive(uintptr_t* kinship_table, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  PglErr reterr = kPglRetSuccess;
  {
    const uintptr_t orig_sample_ct = *sample_ct_ptr;
    const uintptr_t orig_sample_ctl = BitCtToWordCt(orig_sample_ct);
    uintptr_t* sample_include_collapsed_nz;
    uintptr_t* sample_remove_collapsed;
    uint32_t* vertex_degree;
    if (unlikely(bigstack_calloc_w(orig_sample_ctl, &sample_include_collapsed_nz) ||
                 bigstack_calloc_w(orig_sample_ctl, &sample_remove_collapsed) ||
                 bigstack_alloc_u32(orig_sample_ct, &vertex_degree))) {
      goto KinshipPruneDestructive_ret_NOMEM;
    }
    // 1. count the number of constraints for each remaining sample
    uint32_t degree_1_vertex_ct = 0;
    for (uint32_t sample_idx = 0; sample_idx != orig_sample_ct; ++sample_idx) {
      const uintptr_t woffset = sample_idx * orig_sample_ctl;
      const uintptr_t* read_iter1 = &(kinship_table[woffset]);
      // don't currently guarantee vector-alignment of kinship_table rows, so
      // can't use PopcountWords().  (change this?)
      uint32_t cur_degree = 0;
      for (uint32_t widx = 0; widx != orig_sample_ctl; ++widx) {
        const uintptr_t cur_word = *read_iter1++;
        cur_degree += PopcountWord(cur_word);
      }
      if (cur_degree) {
        vertex_degree[sample_idx] = cur_degree;
        degree_1_vertex_ct += (cur_degree == 1);
        SetBit(sample_idx, sample_include_collapsed_nz);
      }
    }
    uint32_t cur_sample_nz_ct = PopcountWords(sample_include_collapsed_nz, orig_sample_ctl);
    // 2. as long as edges remain,
    //    a. remove partner of first degree-one vertex, if such a vertex exists
    //    b. otherwise, remove first maximal-degree vertex
    //    (similar to plink 1.9 rel_cutoff_batch(), but data structure is not
    //    triangular since more speed is needed)
    while (cur_sample_nz_ct) {
      uint32_t prune_uidx;
      uint32_t cur_degree;
      if (degree_1_vertex_ct) {
        uint32_t degree_1_vertex_uidx = 0;
        while (1) {
          // sparse
          degree_1_vertex_uidx = AdvTo1Bit(sample_include_collapsed_nz, degree_1_vertex_uidx);
          if (vertex_degree[degree_1_vertex_uidx] == 1) {
            break;
          }
          ++degree_1_vertex_uidx;
        }
        // find partner
        prune_uidx = AdvTo1Bit(&(kinship_table[degree_1_vertex_uidx * orig_sample_ctl]), 0);
        cur_degree = vertex_degree[prune_uidx];
      } else {
        uint32_t sample_uidx = AdvTo1Bit(sample_include_collapsed_nz, 0);
        cur_degree = vertex_degree[sample_uidx];
        prune_uidx = sample_uidx;
        for (uint32_t sample_idx = 1; sample_idx != cur_sample_nz_ct; ++sample_idx) {
          // sparse
          sample_uidx = AdvTo1Bit(sample_include_collapsed_nz, sample_uidx + 1);
          const uint32_t new_degree = vertex_degree[sample_uidx];
          if (new_degree > cur_degree) {
            cur_degree = new_degree;
            prune_uidx = sample_uidx;
          }
        }
      }
      // remove row/column
      uintptr_t* cur_kinship_col = &(kinship_table[prune_uidx / kBitsPerWord]);
      const uintptr_t kinship_col_mask = ~(k1LU << (prune_uidx % kBitsPerWord));
      uintptr_t* cur_kinship_row = &(kinship_table[prune_uidx * orig_sample_ctl]);
      uint32_t sample_uidx = 0;
      for (uint32_t partner_idx = 0; partner_idx != cur_degree; ++partner_idx, ++sample_uidx) {
        // sparse
        sample_uidx = AdvTo1Bit(cur_kinship_row, sample_uidx);
        const uint32_t new_degree = vertex_degree[sample_uidx] - 1;
        if (!new_degree) {
          ClearBit(sample_uidx, sample_include_collapsed_nz);
          --degree_1_vertex_ct;
          --cur_sample_nz_ct;
          // unnecessary to write to kinship_table[] or vertex_degree[]
        } else {
          cur_kinship_col[sample_uidx * orig_sample_ctl] &= kinship_col_mask;
          degree_1_vertex_ct += (new_degree == 1);
          vertex_degree[sample_uidx] = new_degree;
        }
      }
      if (vertex_degree[prune_uidx] == 1) {
        --degree_1_vertex_ct;
      }
      sample_remove_collapsed[prune_uidx / kBitsPerWord] |= ~kinship_col_mask;
      sample_include_collapsed_nz[prune_uidx / kBitsPerWord] &= kinship_col_mask;
      // unnecessary to update current kinship_table[] row
      --cur_sample_nz_ct;
    }
    uint32_t sample_ct = orig_sample_ct;
    uintptr_t sample_widx = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != orig_sample_ct; ++sample_idx) {
      const uintptr_t lowbit = BitIter1y(sample_include, &sample_widx, &cur_bits);
      if (IsSet(sample_remove_collapsed, sample_idx)) {
        sample_include[sample_widx] ^= lowbit;
        --sample_ct;
      }
    }
    *sample_ct_ptr = sample_ct;
  }
  while (0) {
  KinshipPruneDestructive_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
  return reterr;
}

PglErr KingCutoffBatchBinary(const SampleIdInfo* siip, uint32_t raw_sample_ct, double king_cutoff, uintptr_t* sample_include, char* king_cutoff_fprefix, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* binfile = nullptr;
  char* fprefix_end = &(king_cutoff_fprefix[strlen(king_cutoff_fprefix)]);
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    uint32_t sample_ct = *sample_ct_ptr;
    const uint32_t orig_sample_ctl = BitCtToWordCt(sample_ct);
    uintptr_t* kinship_table;
    uint32_t* sample_uidx_to_king_uidx;
    if (unlikely(bigstack_calloc_w(sample_ct * orig_sample_ctl, &kinship_table) ||
                 bigstack_alloc_u32(raw_sample_ct, &sample_uidx_to_king_uidx))) {
      goto KingCutoffBatchBinary_ret_NOMEM;
    }

    snprintf(fprefix_end, 9, ".king.id");
    reterr = InitTextStream(king_cutoff_fprefix, kTextStreamBlenFast, 1, &txs);
    if (unlikely(reterr)) {
      goto KingCutoffBatchBinary_ret_TSTREAM_FAIL;
    }
    // bugfix (18 Aug 2018): this missed some xid_mode possibilities
    // todo: try to simplify this interface, it's bordering on incomprehensible
    char* line_start;
    XidMode xid_mode;
    reterr = LoadXidHeader("king-cutoff", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeader0 : kfXidHeaderIgnoreSid, &line_idx, &txs, &xid_mode, &line_start);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: Empty --king-cutoff ID file.\n");
        goto KingCutoffBatchBinary_ret_MALFORMED_INPUT;
      }
      goto KingCutoffBatchBinary_ret_TSTREAM_XID_FAIL;
    }

    uint32_t* xid_map;  // IDs not collapsed
    char* sorted_xidbox;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto KingCutoffBatchBinary_ret_1;
    }
    char* idbuf;
    if (unlikely(bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto KingCutoffBatchBinary_ret_NOMEM;
    }
    SetAllU32Arr(raw_sample_ct, sample_uidx_to_king_uidx);
    uintptr_t king_id_ct = 0;
    if (*line_start == '#') {
      ++line_idx;
      line_start = TextGet(&txs);
    }
    for (; line_start; ++line_idx, line_start = TextGet(&txs)) {
      const char* linebuf_iter = line_start;
      uint32_t sample_uidx;
      if (SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, sample_ct, 0, xid_mode, &linebuf_iter, &sample_uidx, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto KingCutoffBatchBinary_ret_MISSING_TOKENS;
        }
        continue;
      }
      if (unlikely(sample_uidx_to_king_uidx[sample_uidx] != UINT32_MAX)) {
        char* first_tab = AdvToDelim(idbuf, '\t');
        char* second_tab = strchr(&(first_tab[1]), '\t');
        *first_tab = ' ';
        if (second_tab) {
          *second_tab = ' ';
        }
        snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID \"%s\" in %s .\n", idbuf, king_cutoff_fprefix);
        goto KingCutoffBatchBinary_ret_MALFORMED_INPUT_WW;
      }
      sample_uidx_to_king_uidx[sample_uidx] = king_id_ct;
      ++king_id_ct;
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto KingCutoffBatchBinary_ret_TSTREAM_FAIL;
    }

    BigstackReset(TextStreamMemStart(&txs));
    if (unlikely(CleanupTextStream2(king_cutoff_fprefix, &txs, &reterr))) {
      goto KingCutoffBatchBinary_ret_1;
    }
    uintptr_t* king_include;
    uint32_t* king_uidx_to_sample_idx;
    if (unlikely(bigstack_calloc_w(BitCtToWordCt(king_id_ct), &king_include) ||
                 bigstack_alloc_u32(king_id_ct, &king_uidx_to_sample_idx))) {
      goto KingCutoffBatchBinary_ret_NOMEM;
    }
    uintptr_t sample_uidx_base = 0;
    uintptr_t sample_include_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
      const uint32_t king_uidx = sample_uidx_to_king_uidx[sample_uidx];
      if (king_uidx != UINT32_MAX) {
        SetBit(king_uidx, king_include);
        king_uidx_to_sample_idx[king_uidx] = sample_idx;
      }
    }
    snprintf(fprefix_end, 10, ".king.bin");
    if (unlikely(fopen_checked(king_cutoff_fprefix, FOPEN_RB, &binfile))) {
      goto KingCutoffBatchBinary_ret_OPEN_FAIL;
    }
    if (unlikely(fseeko(binfile, 0, SEEK_END))) {
      goto KingCutoffBatchBinary_ret_READ_FAIL;
    }
    const uint64_t fsize = ftello(binfile);
    const uint64_t fsize_double_expected = (king_id_ct * (S_CAST(uint64_t, king_id_ct) - 1) * (sizeof(double) / 2));
    const uint32_t is_double = (fsize == fsize_double_expected);
    rewind(binfile);
    const uint32_t first_king_uidx = AdvBoundedTo1Bit(king_include, 0, king_id_ct);
    uintptr_t king_uidx = AdvBoundedTo1Bit(king_include, first_king_uidx + 1, king_id_ct);
    if (king_uidx > 1) {
      if (fseeko(binfile, king_uidx * (S_CAST(uint64_t, king_uidx) - 1) * (2 + (2 * is_double)), SEEK_SET)) {
        goto KingCutoffBatchBinary_ret_READ_FAIL;
      }
    }
    uintptr_t constraint_ct = 0;
    if (is_double) {
      // fread limit
      assert(king_id_ct <= ((kMaxBytesPerIO / sizeof(double)) + 1));
      double* king_drow;
      if (unlikely(bigstack_alloc_d(king_id_ct - 1, &king_drow))) {
        goto KingCutoffBatchBinary_ret_NOMEM;
      }
      for (uint32_t king_idx = 1; king_uidx != king_id_ct; ++king_idx, ++king_uidx) {
        if (!IsSet(king_include, king_uidx)) {
          king_uidx = AdvBoundedTo1Bit(king_include, king_uidx + 1, king_id_ct);
          if (king_uidx == king_id_ct) {
            break;
          }
          if (unlikely(fseeko(binfile, S_CAST(uint64_t, king_uidx) * (king_uidx - 1) * (sizeof(double) / 2), SEEK_SET))) {
            goto KingCutoffBatchBinary_ret_READ_FAIL;
          }
        }
        if (unlikely(!fread_unlocked(king_drow, king_uidx * sizeof(double), 1, binfile))) {
          goto KingCutoffBatchBinary_ret_READ_FAIL;
        }
        const uintptr_t sample_idx = king_uidx_to_sample_idx[king_uidx];
        uintptr_t* kinship_table_row = &(kinship_table[sample_idx * orig_sample_ctl]);
        uintptr_t* kinship_table_col = &(kinship_table[sample_idx / kBitsPerWord]);
        const uintptr_t kinship_new_bit = k1LU << (sample_idx % kBitsPerWord);
        uintptr_t king_uidx2_base;
        uintptr_t king_include_bits;
        BitIter1Start(king_include, first_king_uidx, &king_uidx2_base, &king_include_bits);
        for (uint32_t king_idx2 = 0; king_idx2 != king_idx; ++king_idx2) {
          const uintptr_t king_uidx2 = BitIter1(king_include, &king_uidx2_base, &king_include_bits);
          if (king_drow[king_uidx2] > king_cutoff) {
            const uintptr_t sample_idx2 = king_uidx_to_sample_idx[king_uidx2];
            SetBit(sample_idx2, kinship_table_row);
            kinship_table_col[sample_idx2 * orig_sample_ctl] |= kinship_new_bit;
            ++constraint_ct;
          }
        }
      }
    } else {
      if (unlikely(fsize != (fsize_double_expected / 2))) {
        const uint64_t fsize_double_square = king_id_ct * S_CAST(uint64_t, king_id_ct) * sizeof(double);
        if ((fsize == fsize_double_square) || (fsize == fsize_double_square / 2)) {
          logerrputs("Error: --king-cutoff currently requires a *triangular* .bin file; the provided\nfile appears to be square.\n");
        } else {
          logerrprintfww("Error: Invalid --king-cutoff .bin file size (expected %" PRIu64 " or %" PRIu64 " bytes).\n", fsize_double_expected / 2, fsize_double_expected);
        }
        goto KingCutoffBatchBinary_ret_MALFORMED_INPUT;
      }
      assert(king_id_ct <= ((0x7ffff000 / sizeof(float)) + 1));
      const float king_cutoff_f = S_CAST(float, king_cutoff);
      float* king_frow;
      if (unlikely(bigstack_alloc_f(king_id_ct - 1, &king_frow))) {
        goto KingCutoffBatchBinary_ret_NOMEM;
      }
      for (uint32_t king_idx = 1; king_uidx != king_id_ct; ++king_idx, ++king_uidx) {
        if (!IsSet(king_include, king_uidx)) {
          king_uidx = AdvBoundedTo1Bit(king_include, king_uidx + 1, king_id_ct);
          if (king_uidx == king_id_ct) {
            break;
          }
          if (unlikely(fseeko(binfile, S_CAST(uint64_t, king_uidx) * (king_uidx - 1) * (sizeof(float) / 2), SEEK_SET))) {
            goto KingCutoffBatchBinary_ret_READ_FAIL;
          }
        }
        if (unlikely(!fread_unlocked(king_frow, king_uidx * sizeof(float), 1, binfile))) {
          goto KingCutoffBatchBinary_ret_READ_FAIL;
        }
        const uintptr_t sample_idx = king_uidx_to_sample_idx[king_uidx];
        uintptr_t* kinship_table_row = &(kinship_table[sample_idx * orig_sample_ctl]);
        uintptr_t* kinship_table_col = &(kinship_table[sample_idx / kBitsPerWord]);
        const uintptr_t kinship_new_bit = k1LU << (sample_idx % kBitsPerWord);
        uintptr_t king_uidx2_base;
        uintptr_t king_include_bits;
        BitIter1Start(king_include, first_king_uidx, &king_uidx2_base, &king_include_bits);
        for (uint32_t king_idx2 = 0; king_idx2 != king_idx; ++king_idx2) {
          const uintptr_t king_uidx2 = BitIter1(king_include, &king_uidx2_base, &king_include_bits);
          if (king_frow[king_uidx2] > king_cutoff_f) {
            const uintptr_t sample_idx2 = king_uidx_to_sample_idx[king_uidx2];
            SetBit(sample_idx2, kinship_table_row);
            kinship_table_col[sample_idx2 * orig_sample_ctl] |= kinship_new_bit;
            ++constraint_ct;
          }
        }
      }
    }
    logprintf("--king-cutoff: %" PRIuPTR " constraint%s loaded.\n", constraint_ct, (constraint_ct == 1)? "" : "s");
    BigstackReset(sample_uidx_to_king_uidx);
    if (unlikely(KinshipPruneDestructive(kinship_table, sample_include, sample_ct_ptr))) {
      goto KingCutoffBatchBinary_ret_NOMEM;
    }
  }
  while (0) {
  KingCutoffBatchBinary_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  KingCutoffBatchBinary_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  KingCutoffBatchBinary_ret_READ_FAIL:
    if (feof_unlocked(binfile)) {
      errno = 0;
    }
    logerrprintfww(kErrprintfFread, king_cutoff_fprefix, rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  KingCutoffBatchBinary_ret_MISSING_TOKENS:
    logerrprintfww("Error: Fewer tokens than expected on line %" PRIuPTR " of %s .\n", line_idx, king_cutoff_fprefix);
    reterr = kPglRetMalformedInput;
    break;
  KingCutoffBatchBinary_ret_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&txs)) {
      break;
    }
  KingCutoffBatchBinary_ret_TSTREAM_FAIL:
    TextStreamErrPrint(king_cutoff_fprefix, &txs);
    break;
  KingCutoffBatchBinary_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  KingCutoffBatchBinary_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 KingCutoffBatchBinary_ret_1:
  fclose_cond(binfile);
  if (CleanupTextStream(&txs, &reterr)) {
    snprintf(fprefix_end, 9, ".king.id");
    logerrprintfww(kErrprintfFread, king_cutoff_fprefix, rstrerror(errno));
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr KingCutoffBatchTable(const SampleIdInfo* siip, const char* kin0_fname, uint32_t raw_sample_ct, double king_cutoff, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    // defend from ScanadvDouble() rounding error
    king_cutoff *= 1.0 + kSmallEpsilon;

    const uint32_t orig_sample_ct = *sample_ct_ptr;
    const uint32_t orig_sample_ctl = BitCtToWordCt(orig_sample_ct);
    uintptr_t* kinship_table;
    if (unlikely(bigstack_calloc_w(orig_sample_ct * orig_sample_ctl, &kinship_table))) {
      goto KingCutoffBatchTable_ret_NOMEM;
    }

    // possible todo: deduplicate with CalcKingTableSubset() .kin0 header-line
    // scanner.
    uint32_t kinship_skip = 0;
    XidMode xid_mode = siip->sids? kfXidModeFidIidSid : kfXidModeIidSid;
    reterr = InitTextStream(kin0_fname, kTextStreamBlenFast, 1, &txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: Empty --king-cutoff-table file.\n");
        goto KingCutoffBatchTable_ret_MALFORMED_INPUT;
      }
      goto KingCutoffBatchTable_ret_TSTREAM_FAIL;
    }
    ++line_idx;
    {
      // Parse header line.
      const char* linebuf_iter = TextGet(&txs);
      if (unlikely(!linebuf_iter)) {
        if (!TextStreamErrcode2(&txs, &reterr)) {
          logerrputs("Error: Empty --king-cutoff-table file.\n");
          goto KingCutoffBatchTable_ret_MALFORMED_INPUT;
        }
        goto KingCutoffBatchTable_ret_TSTREAM_FAIL;
      }
      const char* token_end = CurTokenEnd(linebuf_iter);
      uint32_t token_slen = token_end - linebuf_iter;
      // Make this work with both KING- and plink2-generated .kin0 files.
      const uint32_t fid_present = strequal_k(linebuf_iter, "#FID1", token_slen) || strequal_k(linebuf_iter, "FID", token_slen);
      if (fid_present) {
        linebuf_iter = FirstNonTspace(token_end);
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
        xid_mode = kfXidModeFidIid;
      } else {
        if (unlikely(*linebuf_iter != '#')) {
          goto KingCutoffBatchTable_ret_INVALID_HEADER;
        }
        ++linebuf_iter;
        --token_slen;
        xid_mode = kfXidModeIid;
      }
      if (unlikely((!strequal_k(linebuf_iter, "ID1", token_slen)) && (!strequal_k(linebuf_iter, "IID1", token_slen)))) {
        goto KingCutoffBatchTable_ret_INVALID_HEADER;
      }
      linebuf_iter = FirstNonTspace(token_end);
      token_end = CurTokenEnd(linebuf_iter);
      token_slen = token_end - linebuf_iter;
      if (strequal_k(linebuf_iter, "SID1", token_slen)) {
        if (siip->sids) {
          xid_mode = fid_present? kfXidModeFidIidSid : kfXidModeIidSid;
        } else {
          xid_mode |= kfXidModeFlagSkipSid;
        }
        linebuf_iter = FirstNonTspace(token_end);
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
      }
      if (fid_present) {
        if (unlikely(!strequal_k(linebuf_iter, "FID2", token_slen))) {
          goto KingCutoffBatchTable_ret_INVALID_HEADER;
        }
        linebuf_iter = FirstNonTspace(token_end);
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
      }
      if (unlikely((!strequal_k(linebuf_iter, "ID2", token_slen)) && (!strequal_k(linebuf_iter, "IID2", token_slen)))) {
        goto KingCutoffBatchTable_ret_INVALID_HEADER;
      }
      if (xid_mode & (kfXidModeFlagSid | kfXidModeFlagSkipSid)) {
        linebuf_iter = FirstNonTspace(token_end);
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
        if (unlikely(!strequal_k(linebuf_iter, "SID2", token_slen))) {
          goto KingCutoffBatchTable_ret_INVALID_HEADER;
        }
      }
      while (1) {
        linebuf_iter = FirstNonTspace(token_end);
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
        if (unlikely(!token_slen)) {
          logerrputs("Error: No kinship-coefficient column in --king-cutoff-table file.\n");
          goto KingCutoffBatchTable_ret_INCONSISTENT_INPUT;
        }
        if (strequal_k(linebuf_iter, "KINSHIP", token_slen) || strequal_k(linebuf_iter, "Kinship", token_slen)) {
          break;
        }
        ++kinship_skip;
      }
    }

    uint32_t* xid_cmap;
    char* sorted_xidbox;
    uintptr_t max_xid_blen;
    {
      uint32_t* xid_umap;  // IDs not collapsed
      reterr = SortedXidboxInitAlloc(sample_include, siip, orig_sample_ct, xid_mode, 0, &sorted_xidbox, &xid_umap, &max_xid_blen);
      if (unlikely(reterr)) {
        goto KingCutoffBatchTable_ret_1;
      }
      xid_cmap = xid_umap;
      if (orig_sample_ct != raw_sample_ct) {
        // Convert to collapsed IDs in-place.
        const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
        uint32_t* sample_include_cumulative_popcounts;
        if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts))) {
          goto KingCutoffBatchTable_ret_NOMEM;
        }
        FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
        UidxsToIdxs(sample_include, sample_include_cumulative_popcounts, orig_sample_ct, xid_cmap);
        BigstackReset(sample_include_cumulative_popcounts);
      }
    }
    char* idbuf;
    if (unlikely(bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto KingCutoffBatchTable_ret_NOMEM;
    }

    // possible todo: deduplicate with KingTableSubsetLoad().
    ++line_idx;
    uintptr_t constraint_ct = 0;
    for (char* line_iter = TextLineEnd(&txs); TextGetUnsafe2(&txs, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      const char* linebuf_iter = line_iter;
      uint32_t sample_idx1;
      if (SortedXidboxReadFind(sorted_xidbox, xid_cmap, max_xid_blen, orig_sample_ct, 0, xid_mode, &linebuf_iter, &sample_idx1, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto KingCutoffBatchTable_ret_MISSING_TOKENS;
        }
        line_iter = K_CAST(char*, linebuf_iter);
        continue;
      }
      linebuf_iter = FirstNonTspace(linebuf_iter);
      uint32_t sample_idx2;
      if (SortedXidboxReadFind(sorted_xidbox, xid_cmap, max_xid_blen, orig_sample_ct, 0, xid_mode, &linebuf_iter, &sample_idx2, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto KingCutoffBatchTable_ret_MISSING_TOKENS;
        }
        line_iter = K_CAST(char*, linebuf_iter);
        continue;
      }
      if (unlikely(sample_idx1 == sample_idx2)) {
        // could technically be due to unloaded SID, so use inconsistent-input
        // error code
        snprintf(g_logbuf, kLogbufSize, "Error: Identical sample IDs on line %" PRIuPTR " of --king-cutoff-table file.\n", line_idx);
        goto KingCutoffBatchTable_ret_INCONSISTENT_INPUT_WW;
      }
      linebuf_iter = FirstNonTspace(linebuf_iter);
      linebuf_iter = NextTokenMult0(linebuf_iter, kinship_skip);
      if (unlikely(!linebuf_iter)) {
        goto KingCutoffBatchTable_ret_MISSING_TOKENS;
      }
      double cur_kinship;
      const char* kinship_end = ScanadvDouble(linebuf_iter, &cur_kinship);
      if (!kinship_end) {
        line_iter = K_CAST(char*, linebuf_iter);
        continue;
      }
      if (unlikely(!IsSpaceOrEoln(*kinship_end))) {
        kinship_end = CurTokenEnd(kinship_end);
        *K_CAST(char*, kinship_end) = '\0';
        logerrprintfww("Error: Invalid numeric token '%s' on line %" PRIuPTR " of --king-cutoff-table file.\n", linebuf_iter, line_idx);
        goto KingCutoffBatchTable_ret_MALFORMED_INPUT;
      }
      if (cur_kinship > king_cutoff) {
        SetBit(sample_idx2, &(kinship_table[sample_idx1 * orig_sample_ctl]));
        SetBit(sample_idx1, &(kinship_table[sample_idx2 * orig_sample_ctl]));
        ++constraint_ct;
      }
    }

    logprintf("--king-cutoff-table: %" PRIuPTR " constraint%s loaded.\n", constraint_ct, (constraint_ct == 1)? "" : "s");
    if (unlikely(KinshipPruneDestructive(kinship_table, sample_include, sample_ct_ptr))) {
      goto KingCutoffBatchTable_ret_NOMEM;
    }
  }
  while (0) {
  KingCutoffBatchTable_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  KingCutoffBatchTable_ret_MISSING_TOKENS:
    logerrprintfww("Error: Fewer tokens than expected on line %" PRIuPTR " of %s .\n", line_idx, kin0_fname);
    reterr = kPglRetMalformedInput;
    break;
  KingCutoffBatchTable_ret_TSTREAM_FAIL:
    TextStreamErrPrint(kin0_fname, &txs);
    break;
  KingCutoffBatchTable_ret_INVALID_HEADER:
    logerrputs("Error: Invalid header line in --king-cutoff-table file.\n");
    reterr = kPglRetMalformedInput;
    break;
  KingCutoffBatchTable_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  KingCutoffBatchTable_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  KingCutoffBatchTable_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 KingCutoffBatchTable_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

CONSTI32(kKingOffsetIbs0, 0);
CONSTI32(kKingOffsetHethet, 1);
CONSTI32(kKingOffsetHet2Hom1, 2);
CONSTI32(kKingOffsetHet1Hom2, 3);
CONSTI32(kKingOffsetHomhom, 4);

typedef struct CalcKingSparseCtxStruct {
  const uintptr_t* variant_include_orig;
  uintptr_t* sample_include;
  uint32_t* sample_include_cumulative_popcounts;
  uint32_t row_start_idx;
  uint32_t row_end_idx;
  uint32_t homhom_needed;

  uint32_t max_sparse_ct;

  uint32_t read_block_size;  // guaranteed to be power of 2

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uint32_t* read_variant_uidx_starts;

  // this has length >= 3 * max_sparse_ct
  uint32_t** thread_idx_bufs;

  uint32_t cur_block_size;

  uint32_t** thread_singleton_het_cts;
  uint32_t** thread_singleton_hom_cts;
  uint32_t** thread_singleton_missing_cts;
  uint32_t* thread_skip_cts;

  // single global copy
  uint32_t* king_counts;

  uintptr_t** thread_sparse_excludes[2];

  uint64_t err_info;
} CalcKingSparseCtx;

THREAD_FUNC_DECL CalcKingSparseThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcKingSparseCtx* ctx = S_CAST(CalcKingSparseCtx*, arg->sharedp->context);

  const uintptr_t* variant_include_orig = ctx->variant_include_orig;
  const uintptr_t* sample_include = ctx->sample_include;

  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  uintptr_t* genovec = ctx->genovecs[tidx];
  uint32_t row_start_idx = ctx->row_start_idx;
  const uint64_t tri_start = ((row_start_idx - 1) * S_CAST(uint64_t, row_start_idx)) / 2;
  if (row_start_idx == 1) {
    row_start_idx = 0;
  }
  const uint32_t sample_ct = ctx->row_end_idx;
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t remainder = sample_ct % kBitsPerWordD2;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uint32_t homhom_needed = ctx->homhom_needed;
  const uintptr_t homhom_needed_p4 = homhom_needed + 4;
  const uint32_t max_sparse_ct = ctx->max_sparse_ct;
  const uint32_t read_block_size_mask = ctx->read_block_size - 1;
  const uint32_t read_block_sizel = ctx->read_block_size / kBitsPerWord;

  uint32_t* idx_bufs[4];
  idx_bufs[0] = nullptr;
  idx_bufs[1] = ctx->thread_idx_bufs[tidx];
  idx_bufs[2] = &(idx_bufs[1][max_sparse_ct]);
  idx_bufs[3] = &(idx_bufs[2][max_sparse_ct]);
  const uint32_t min_common_ct = sample_ct - max_sparse_ct;

  uint32_t* singleton_het_cts = ctx->thread_singleton_het_cts[tidx];
  uint32_t* singleton_hom_cts = ctx->thread_singleton_hom_cts[tidx];
  uint32_t* singleton_missing_cts = ctx->thread_singleton_missing_cts[tidx];
  ZeroU32Arr(sample_ct, singleton_het_cts);
  ZeroU32Arr(sample_ct, singleton_hom_cts);
  ZeroU32Arr(sample_ct, singleton_missing_cts);
  uint32_t skip_ct = 0;

  uint32_t* king_counts = ctx->king_counts;
  uint64_t new_err_info = 0;
  {
    // This matrix can be huge, so we multithread zero-initialization.
    const uint64_t entry_ct = homhom_needed_p4 * (((sample_ct - 1) * S_CAST(uint64_t, sample_ct)) / 2 - tri_start);
    const uintptr_t fill_start = RoundDownPow2((tidx * entry_ct) / calc_thread_ct, kInt32PerCacheline);
    uintptr_t fill_end = entry_ct;
    if (tidx + 1 != calc_thread_ct) {
      fill_end = RoundDownPow2(((tidx + 1) * entry_ct) / calc_thread_ct, kInt32PerCacheline);
    }
    ZeroU32Arr(fill_end - fill_start, &(king_counts[fill_start]));
  }
  uint32_t parity = 0;
  // sync.Once before main loop; we need the other threads to be done with
  // their zero-initialization jobs before we can proceed.
  while (!THREAD_BLOCK_FINISH(arg)) {
    const uint32_t cur_block_size = ctx->cur_block_size;
    const uint32_t idx_end = ((tidx + 1) * cur_block_size) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include_orig, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);
    uintptr_t* sparse_exclude = ctx->thread_sparse_excludes[parity][tidx];
    ZeroWArr(read_block_sizel, sparse_exclude);
    // probable todo: better load-balancing
    for (uint32_t cur_idx = (tidx * cur_block_size) / calc_thread_ct; cur_idx != idx_end; ++cur_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include_orig, &variant_uidx_base, &variant_include_bits);
      // tried DifflistOrGenovec, difference was negligible.  Not really worth
      // considering it when calculation is inherently >O(mn).
      const PglErr reterr = PgrGet(sample_include, pssi, sample_ct, variant_uidx, pgrp, genovec);
      if (unlikely(reterr)) {
        new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
        goto CalcKingSparseThread_err;
      }
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      ZeroTrailingNyps(sample_ct, genovec);
      GenoarrCountFreqsUnsafe(genovec, sample_ct, genocounts);
      uintptr_t mask_word;
      uintptr_t common_idx;
      if (genocounts[0] >= min_common_ct) {
        common_idx = 0;
        mask_word = 0;
      } else if (genocounts[2] >= min_common_ct) {
        common_idx = 2;
        mask_word = kMaskAAAA;
      } else if (genocounts[3] >= min_common_ct) {
        common_idx = 3;
        mask_word = ~k0LU;
        ++skip_ct;
      } else {
        if ((!homhom_needed) && ((genocounts[0] + genocounts[3] == sample_ct) || (genocounts[2] + genocounts[3] == sample_ct))) {
          SetBit(variant_uidx & read_block_size_mask, sparse_exclude);
          ++skip_ct;
        }
        continue;
      }
      SetBit(variant_uidx & read_block_size_mask, sparse_exclude);
      if (genocounts[common_idx] == sample_ct) {
        continue;
      }
      if (remainder) {
        genovec[sample_ctl2 - 1] |= mask_word << (2 * remainder);
      }
      uint32_t* idx_buf_iters[4];
      memcpy(idx_buf_iters, idx_bufs, 4 * sizeof(intptr_t));
      for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
        uintptr_t xor_word = genovec[widx] ^ mask_word;
        if (xor_word) {
          const uint32_t offset_base = widx * kBitsPerWordD2;
          do {
            const uint32_t shift_ct = ctzw(xor_word) & (~1);
            const uint32_t cur_xor = (xor_word >> shift_ct) & 3;
            *(idx_buf_iters[cur_xor])++ = offset_base + (shift_ct / 2);
            xor_word &= ~((3 * k1LU) << shift_ct);
          } while (xor_word);
        }
      }
      // We do two things here.
      // 1. Update singleton_{het,hom,missing}_cts for every observed rare
      //    genotype.  This is enough for correct accounting for any pair
      //    involving only one (or none) of these rare genotypes, and the
      //    arrays are small enough that each thread can keep its own copy
      //    (they're added up at the end).
      // 2. For each pair of rare genotypes, atomically correct the main
      //    king_counts[] array.  This is messy (9x2 cases) but conceptually
      //    straightforward.
      const uint32_t* het_idxs = idx_bufs[common_idx ^ 1];
      const uint32_t het_ct = genocounts[1];
      if (common_idx != 3) {
        const uint32_t* other_hom_idxs = idx_bufs[2];
        const uint32_t* missing_idxs = idx_bufs[common_idx ^ 3];
        const uint32_t other_hom_ct = idx_buf_iters[2] - other_hom_idxs;
        const uint32_t missing_ct = genocounts[3];
        for (uint32_t uii = 0; uii != het_ct; ++uii) {
          // We want to iterate over one row at a time, for better
          // memory-access locality.  So the outer loop must correspond to the
          // larger sample-index.
          const uintptr_t sample_idx_hi = het_idxs[uii];
          singleton_het_cts[sample_idx_hi] += 1;
          if (sample_idx_hi < row_start_idx) {
            continue;
          }
          const uintptr_t tri_base = (S_CAST(uint64_t, sample_idx_hi) * (sample_idx_hi - 1)) / 2 - tri_start;
          for (uint32_t ujj = 0; ujj != uii; ++ujj) {
            const uintptr_t sample_idx_lo = het_idxs[ujj];
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHethet]), 1, __ATOMIC_RELAXED);
            __atomic_fetch_sub(&(king_counts_ptr[kKingOffsetHet2Hom1]), 1, __ATOMIC_RELAXED);
            __atomic_fetch_sub(&(king_counts_ptr[kKingOffsetHet1Hom2]), 1, __ATOMIC_RELAXED);
            if (homhom_needed) {
              __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHomhom]), 1, __ATOMIC_RELAXED);
            }
          }
          for (uint32_t ujj = 0; ujj != other_hom_ct; ++ujj) {
            const uintptr_t sample_idx_lo = other_hom_idxs[ujj];
            if (sample_idx_lo > sample_idx_hi) {
              break;
            }
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_sub(&(king_counts_ptr[kKingOffsetIbs0]), 1, __ATOMIC_RELAXED);
          }
          for (uint32_t ujj = 0; ujj != missing_ct; ++ujj) {
            const uintptr_t sample_idx_lo = missing_idxs[ujj];
            if (sample_idx_lo > sample_idx_hi) {
              break;
            }
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_sub(&(king_counts_ptr[kKingOffsetHet2Hom1]), 1, __ATOMIC_RELAXED);
            if (homhom_needed) {
              __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHomhom]), 1, __ATOMIC_RELAXED);
            }
          }
        }
        for (uint32_t uii = 0; uii != other_hom_ct; ++uii) {
          const uintptr_t sample_idx_hi = other_hom_idxs[uii];
          singleton_hom_cts[sample_idx_hi] += 1;
          if (sample_idx_hi < row_start_idx) {
            continue;
          }
          const uintptr_t tri_base = (S_CAST(uint64_t, sample_idx_hi) * (sample_idx_hi - 1)) / 2 - tri_start;
          for (uint32_t ujj = 0; ujj != uii; ++ujj) {
            const uintptr_t sample_idx_lo = other_hom_idxs[ujj];
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_sub(&(king_counts_ptr[kKingOffsetIbs0]), 2, __ATOMIC_RELAXED);
          }
          for (uint32_t ujj = 0; ujj != het_ct; ++ujj) {
            const uintptr_t sample_idx_lo = het_idxs[ujj];
            if (sample_idx_lo > sample_idx_hi) {
              break;
            }
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_sub(&(king_counts_ptr[kKingOffsetIbs0]), 1, __ATOMIC_RELAXED);
          }
          for (uint32_t ujj = 0; ujj != missing_ct; ++ujj) {
            const uintptr_t sample_idx_lo = missing_idxs[ujj];
            if (sample_idx_lo > sample_idx_hi) {
              break;
            }
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_sub(&(king_counts_ptr[kKingOffsetIbs0]), 1, __ATOMIC_RELAXED);
            if (homhom_needed) {
              __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHomhom]), 1, __ATOMIC_RELAXED);
            }
          }
        }
        for (uint32_t uii = 0; uii != missing_ct; ++uii) {
          const uintptr_t sample_idx_hi = missing_idxs[uii];
          singleton_missing_cts[sample_idx_hi] += 1;
          if (sample_idx_hi < row_start_idx) {
            continue;
          }
          const uintptr_t tri_base = (S_CAST(uint64_t, sample_idx_hi) * (sample_idx_hi - 1)) / 2 - tri_start;
          if (homhom_needed) {
            for (uint32_t ujj = 0; ujj != uii; ++ujj) {
              const uintptr_t sample_idx_lo = missing_idxs[ujj];
              const uintptr_t tri_coord = tri_base + sample_idx_lo;
              // bugfix (12 Nov 2019): added 4 twice
              uint32_t* king_counts_ptr = &(king_counts[tri_coord * 5]);
              __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHomhom]), 1, __ATOMIC_RELAXED);
            }
          }
          for (uint32_t ujj = 0; ujj != het_ct; ++ujj) {
            const uintptr_t sample_idx_lo = het_idxs[ujj];
            if (sample_idx_lo > sample_idx_hi) {
              break;
            }
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_sub(&(king_counts_ptr[kKingOffsetHet1Hom2]), 1, __ATOMIC_RELAXED);
            if (homhom_needed) {
              __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHomhom]), 1, __ATOMIC_RELAXED);
            }
          }
          for (uint32_t ujj = 0; ujj != other_hom_ct; ++ujj) {
            const uintptr_t sample_idx_lo = other_hom_idxs[ujj];
            if (sample_idx_lo > sample_idx_hi) {
              break;
            }
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_sub(&(king_counts_ptr[kKingOffsetIbs0]), 1, __ATOMIC_RELAXED);
            if (homhom_needed) {
              __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHomhom]), 1, __ATOMIC_RELAXED);
            }
          }
        }
      } else {
        // merge hom0 and hom2 cases.
        for (uint32_t hom_geno = 0; hom_geno != 4; hom_geno += 2) {
          const uint32_t* cur_hom_idxs = idx_bufs[3 - hom_geno];
          const uint32_t* opp_hom_idxs = idx_bufs[1 + hom_geno];
          const uint32_t cur_hom_ct = genocounts[hom_geno];
          const uint32_t opp_hom_ct = genocounts[2 - hom_geno];
          for (uint32_t uii = 0; uii != cur_hom_ct; ++uii) {
            const uintptr_t sample_idx_hi = cur_hom_idxs[uii];
            if (sample_idx_hi < row_start_idx) {
              continue;
            }
            const uintptr_t tri_base = (S_CAST(uint64_t, sample_idx_hi) * (sample_idx_hi - 1)) / 2 - tri_start;
            if (homhom_needed) {
              for (uint32_t ujj = 0; ujj != uii; ++ujj) {
                const uintptr_t sample_idx_lo = cur_hom_idxs[ujj];
                const uintptr_t tri_coord = tri_base + sample_idx_lo;
                uint32_t* king_counts_ptr = &(king_counts[tri_coord * 5]);
                __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHomhom]), 1, __ATOMIC_RELAXED);
              }
            }
            for (uint32_t ujj = 0; ujj != het_ct; ++ujj) {
              const uintptr_t sample_idx_lo = het_idxs[ujj];
              if (sample_idx_lo > sample_idx_hi) {
                break;
              }
              const uintptr_t tri_coord = tri_base + sample_idx_lo;
              uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
              __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHet1Hom2]), 1, __ATOMIC_RELAXED);
            }
            for (uint32_t ujj = 0; ujj != opp_hom_ct; ++ujj) {
              const uintptr_t sample_idx_lo = opp_hom_idxs[ujj];
              if (sample_idx_lo > sample_idx_hi) {
                break;
              }
              const uintptr_t tri_coord = tri_base + sample_idx_lo;
              uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
              __atomic_fetch_add(&(king_counts_ptr[kKingOffsetIbs0]), 1, __ATOMIC_RELAXED);
              if (homhom_needed) {
                __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHomhom]), 1, __ATOMIC_RELAXED);
              }
            }
          }
        }
        const uint32_t* hom0_idxs = idx_bufs[3];
        const uint32_t* hom2_idxs = idx_bufs[1];
        const uint32_t hom0_ct = genocounts[0];
        const uint32_t hom2_ct = genocounts[2];
        for (uint32_t uii = 0; uii != het_ct; ++uii) {
          const uintptr_t sample_idx_hi = het_idxs[uii];
          if (sample_idx_hi < row_start_idx) {
            continue;
          }
          const uintptr_t tri_base = (S_CAST(uint64_t, sample_idx_hi) * (sample_idx_hi - 1)) / 2 - tri_start;
          for (uint32_t ujj = 0; ujj != uii; ++ujj) {
            const uintptr_t sample_idx_lo = het_idxs[ujj];
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHethet]), 1, __ATOMIC_RELAXED);
          }
          for (uint32_t ujj = 0; ujj != hom0_ct; ++ujj) {
            const uintptr_t sample_idx_lo = hom0_idxs[ujj];
            if (sample_idx_lo > sample_idx_hi) {
              break;
            }
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHet2Hom1]), 1, __ATOMIC_RELAXED);
          }
          for (uint32_t ujj = 0; ujj != hom2_ct; ++ujj) {
            const uintptr_t sample_idx_lo = hom2_idxs[ujj];
            if (sample_idx_lo > sample_idx_hi) {
              break;
            }
            const uintptr_t tri_coord = tri_base + sample_idx_lo;
            uint32_t* king_counts_ptr = &(king_counts[tri_coord * homhom_needed_p4]);
            __atomic_fetch_add(&(king_counts_ptr[kKingOffsetHet2Hom1]), 1, __ATOMIC_RELAXED);
          }
        }
      }
    }
    parity = 1 - parity;
  }
  {
    ctx->thread_skip_cts[tidx] = skip_ct;
  }
  while (0) {
  CalcKingSparseThread_err:
    UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    THREAD_BLOCK_FINISH(arg);
    break;
  }
  THREAD_RETURN;
}

#ifdef USE_SSE42
CONSTI32(kKingMultiplex, 1024);
CONSTI32(kKingMultiplexWords, kKingMultiplex / kBitsPerWord);
void IncrKing(const uintptr_t* smaj_hom, const uintptr_t* smaj_ref2het, uint32_t start_idx, uint32_t end_idx, uint32_t* king_counts_iter) {
  // Tried adding another level of blocking, but couldn't get it to make a
  // difference.
  for (uint32_t second_idx = start_idx; second_idx != end_idx; ++second_idx) {
    // technically overflows for huge sample_ct
    const uint32_t second_offset = second_idx * kKingMultiplexWords;
    const uintptr_t* second_hom = &(smaj_hom[second_offset]);
    const uintptr_t* second_ref2het = &(smaj_ref2het[second_offset]);
    const uintptr_t* first_hom_iter = smaj_hom;
    const uintptr_t* first_ref2het_iter = smaj_ref2het;
    while (first_hom_iter < second_hom) {
      uint32_t acc_ibs0 = 0;
      uint32_t acc_hethet = 0;
      uint32_t acc_het2hom1 = 0;
      uint32_t acc_het1hom2 = 0;
      for (uint32_t widx = 0; widx != kKingMultiplexWords; ++widx) {
        const uintptr_t hom1 = first_hom_iter[widx];
        const uintptr_t hom2 = second_hom[widx];
        const uintptr_t ref2het1 = first_ref2het_iter[widx];
        const uintptr_t ref2het2 = second_ref2het[widx];
        const uintptr_t homhom = hom1 & hom2;
        const uintptr_t het1 = ref2het1 & (~hom1);
        const uintptr_t het2 = ref2het2 & (~hom2);
        acc_ibs0 += PopcountWord((ref2het1 ^ ref2het2) & homhom);
        acc_hethet += PopcountWord(het1 & het2);
        acc_het2hom1 += PopcountWord(hom1 & het2);
        acc_het1hom2 += PopcountWord(hom2 & het1);
      }
      king_counts_iter[kKingOffsetIbs0] += acc_ibs0;
      king_counts_iter[kKingOffsetHethet] += acc_hethet;
      king_counts_iter[kKingOffsetHet2Hom1] += acc_het2hom1;
      king_counts_iter[kKingOffsetHet1Hom2] += acc_het1hom2;
      king_counts_iter = &(king_counts_iter[4]);

      first_hom_iter = &(first_hom_iter[kKingMultiplexWords]);
      first_ref2het_iter = &(first_ref2het_iter[kKingMultiplexWords]);
    }
  }
}

void IncrKingHomhom(const uintptr_t* smaj_hom, const uintptr_t* smaj_ref2het, uint32_t start_idx, uint32_t end_idx, uint32_t* king_counts_iter) {
  for (uint32_t second_idx = start_idx; second_idx != end_idx; ++second_idx) {
    // technically overflows for huge sample_ct
    const uint32_t second_offset = second_idx * kKingMultiplexWords;
    const uintptr_t* second_hom = &(smaj_hom[second_offset]);
    const uintptr_t* second_ref2het = &(smaj_ref2het[second_offset]);
    const uintptr_t* first_hom_iter = smaj_hom;
    const uintptr_t* first_ref2het_iter = smaj_ref2het;
    while (first_hom_iter < second_hom) {
      uint32_t acc_homhom = 0;
      uint32_t acc_ibs0 = 0;
      uint32_t acc_hethet = 0;
      uint32_t acc_het2hom1 = 0;
      uint32_t acc_het1hom2 = 0;
      for (uint32_t widx = 0; widx != kKingMultiplexWords; ++widx) {
        const uintptr_t hom1 = first_hom_iter[widx];
        const uintptr_t hom2 = second_hom[widx];
        const uintptr_t ref2het1 = first_ref2het_iter[widx];
        const uintptr_t ref2het2 = second_ref2het[widx];
        const uintptr_t homhom = hom1 & hom2;
        const uintptr_t het1 = ref2het1 & (~hom1);
        const uintptr_t het2 = ref2het2 & (~hom2);
        acc_homhom += PopcountWord(homhom);
        acc_ibs0 += PopcountWord((ref2het1 ^ ref2het2) & homhom);
        acc_hethet += PopcountWord(het1 & het2);
        acc_het2hom1 += PopcountWord(hom1 & het2);
        acc_het1hom2 += PopcountWord(hom2 & het1);
      }
      king_counts_iter[kKingOffsetIbs0] += acc_ibs0;
      king_counts_iter[kKingOffsetHethet] += acc_hethet;
      king_counts_iter[kKingOffsetHet2Hom1] += acc_het2hom1;
      king_counts_iter[kKingOffsetHet1Hom2] += acc_het1hom2;
      king_counts_iter[kKingOffsetHomhom] += acc_homhom;
      king_counts_iter = &(king_counts_iter[5]);

      first_hom_iter = &(first_hom_iter[kKingMultiplexWords]);
      first_ref2het_iter = &(first_ref2het_iter[kKingMultiplexWords]);
    }
  }
}
#else  // !USE_SSE42
#  ifdef __LP64__
CONSTI32(kKingMultiplex, 1536);
#  else
CONSTI32(kKingMultiplex, 960);
#  endif
static_assert(kKingMultiplex % (3 * kBitsPerVec) == 0, "Invalid kKingMultiplex value.");
CONSTI32(kKingMultiplexWords, kKingMultiplex / kBitsPerWord);
CONSTI32(kKingMultiplexVecs, kKingMultiplex / kBitsPerVec);
// expensive PopcountWord().  Use Lauradoux/Walisch accumulators, since
// Harley-Seal requires too many variables.
void IncrKing(const uintptr_t* smaj_hom, const uintptr_t* smaj_ref2het, uint32_t start_idx, uint32_t end_idx, uint32_t* king_counts_iter) {
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  for (uint32_t second_idx = start_idx; second_idx != end_idx; ++second_idx) {
    // technically overflows for huge sample_ct
    const uint32_t second_offset = second_idx * kKingMultiplexWords;
    const VecW* second_hom = R_CAST(const VecW*, &(smaj_hom[second_offset]));
    const VecW* second_ref2het = R_CAST(const VecW*, &(smaj_ref2het[second_offset]));
    const VecW* first_hom_iter = R_CAST(const VecW*, smaj_hom);
    const VecW* first_ref2het_iter = R_CAST(const VecW*, smaj_ref2het);
    while (first_hom_iter < second_hom) {
      UniVec acc_ibs0;
      UniVec acc_hethet;
      UniVec acc_het2hom1;
      UniVec acc_het1hom2;
      acc_ibs0.vw = vecw_setzero();
      acc_hethet.vw = vecw_setzero();
      acc_het2hom1.vw = vecw_setzero();
      acc_het1hom2.vw = vecw_setzero();
      for (uint32_t vec_idx = 0; vec_idx < kKingMultiplexVecs; vec_idx += 3) {
        VecW hom1 = first_hom_iter[vec_idx];
        VecW hom2 = second_hom[vec_idx];
        VecW ref2het1 = first_ref2het_iter[vec_idx];
        VecW ref2het2 = second_ref2het[vec_idx];
        VecW het1 = vecw_and_notfirst(hom1, ref2het1);
        VecW het2 = vecw_and_notfirst(hom2, ref2het2);
        VecW agg_ibs0 = (ref2het1 ^ ref2het2) & (hom1 & hom2);
        VecW agg_hethet = het1 & het2;
        VecW agg_het2hom1 = hom1 & het2;
        VecW agg_het1hom2 = hom2 & het1;
        agg_ibs0 = agg_ibs0 - (vecw_srli(agg_ibs0, 1) & m1);
        agg_hethet = agg_hethet - (vecw_srli(agg_hethet, 1) & m1);
        agg_het2hom1 = agg_het2hom1 - (vecw_srli(agg_het2hom1, 1) & m1);
        agg_het1hom2 = agg_het1hom2 - (vecw_srli(agg_het1hom2, 1) & m1);
        agg_ibs0 = (agg_ibs0 & m2) + (vecw_srli(agg_ibs0, 2) & m2);
        agg_hethet = (agg_hethet & m2) + (vecw_srli(agg_hethet, 2) & m2);
        agg_het2hom1 = (agg_het2hom1 & m2) + (vecw_srli(agg_het2hom1, 2) & m2);
        agg_het1hom2 = (agg_het1hom2 & m2) + (vecw_srli(agg_het1hom2, 2) & m2);

        for (uint32_t offset = 1; offset != 3; ++offset) {
          hom1 = first_hom_iter[vec_idx + offset];
          hom2 = second_hom[vec_idx + offset];
          ref2het1 = first_ref2het_iter[vec_idx + offset];
          ref2het2 = second_ref2het[vec_idx + offset];
          het1 = vecw_and_notfirst(hom1, ref2het1);
          het2 = vecw_and_notfirst(hom2, ref2het2);
          VecW cur_ibs0 = (ref2het1 ^ ref2het2) & (hom1 & hom2);
          VecW cur_hethet = het1 & het2;
          VecW cur_het2hom1 = hom1 & het2;
          VecW cur_het1hom2 = hom2 & het1;
          cur_ibs0 = cur_ibs0 - (vecw_srli(cur_ibs0, 1) & m1);
          cur_hethet = cur_hethet - (vecw_srli(cur_hethet, 1) & m1);
          cur_het2hom1 = cur_het2hom1 - (vecw_srli(cur_het2hom1, 1) & m1);
          cur_het1hom2 = cur_het1hom2 - (vecw_srli(cur_het1hom2, 1) & m1);
          agg_ibs0 += (cur_ibs0 & m2) + (vecw_srli(cur_ibs0, 2) & m2);
          agg_hethet += (cur_hethet & m2) + (vecw_srli(cur_hethet, 2) & m2);
          agg_het2hom1 += (cur_het2hom1 & m2) + (vecw_srli(cur_het2hom1, 2) & m2);
          agg_het1hom2 += (cur_het1hom2 & m2) + (vecw_srli(cur_het1hom2, 2) & m2);
        }
        acc_ibs0.vw = acc_ibs0.vw + (agg_ibs0 & m4) + (vecw_srli(agg_ibs0, 4) & m4);
        acc_hethet.vw = acc_hethet.vw + (agg_hethet & m4) + (vecw_srli(agg_hethet, 4) & m4);
        acc_het2hom1.vw = acc_het2hom1.vw + (agg_het2hom1 & m4) + (vecw_srli(agg_het2hom1, 4) & m4);
        acc_het1hom2.vw = acc_het1hom2.vw + (agg_het1hom2 & m4) + (vecw_srli(agg_het1hom2, 4) & m4);
      }
      const VecW m8 = VCONST_W(kMask00FF);
      acc_ibs0.vw = (acc_ibs0.vw & m8) + (vecw_srli(acc_ibs0.vw, 8) & m8);
      acc_hethet.vw = (acc_hethet.vw & m8) + (vecw_srli(acc_hethet.vw, 8) & m8);
      acc_het2hom1.vw = (acc_het2hom1.vw & m8) + (vecw_srli(acc_het2hom1.vw, 8) & m8);
      acc_het1hom2.vw = (acc_het1hom2.vw & m8) + (vecw_srli(acc_het1hom2.vw, 8) & m8);
      king_counts_iter[kKingOffsetIbs0] += UniVecHsum16(acc_ibs0);
      king_counts_iter[kKingOffsetHethet] += UniVecHsum16(acc_hethet);
      king_counts_iter[kKingOffsetHet2Hom1] += UniVecHsum16(acc_het2hom1);
      king_counts_iter[kKingOffsetHet1Hom2] += UniVecHsum16(acc_het1hom2);
      king_counts_iter = &(king_counts_iter[4]);

      first_hom_iter = &(first_hom_iter[kKingMultiplexVecs]);
      first_ref2het_iter = &(first_ref2het_iter[kKingMultiplexVecs]);
    }
  }
}

void IncrKingHomhom(const uintptr_t* smaj_hom, const uintptr_t* smaj_ref2het, uint32_t start_idx, uint32_t end_idx, uint32_t* king_counts_iter) {
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  for (uint32_t second_idx = start_idx; second_idx != end_idx; ++second_idx) {
    // technically overflows for huge sample_ct
    const uint32_t second_offset = second_idx * kKingMultiplexWords;
    const VecW* second_hom = R_CAST(const VecW*, &(smaj_hom[second_offset]));
    const VecW* second_ref2het = R_CAST(const VecW*, &(smaj_ref2het[second_offset]));
    const VecW* first_hom_iter = R_CAST(const VecW*, smaj_hom);
    const VecW* first_ref2het_iter = R_CAST(const VecW*, smaj_ref2het);
    while (first_hom_iter < second_hom) {
      UniVec acc_homhom;
      UniVec acc_ibs0;
      UniVec acc_hethet;
      UniVec acc_het2hom1;
      UniVec acc_het1hom2;
      acc_homhom.vw = vecw_setzero();
      acc_ibs0.vw = vecw_setzero();
      acc_hethet.vw = vecw_setzero();
      acc_het2hom1.vw = vecw_setzero();
      acc_het1hom2.vw = vecw_setzero();
      for (uint32_t vec_idx = 0; vec_idx < kKingMultiplexVecs; vec_idx += 3) {
        VecW hom1 = first_hom_iter[vec_idx];
        VecW hom2 = second_hom[vec_idx];
        VecW ref2het1 = first_ref2het_iter[vec_idx];
        VecW ref2het2 = second_ref2het[vec_idx];
        VecW agg_homhom = hom1 & hom2;
        VecW het1 = vecw_and_notfirst(hom1, ref2het1);
        VecW het2 = vecw_and_notfirst(hom2, ref2het2);
        VecW agg_ibs0 = (ref2het1 ^ ref2het2) & agg_homhom;
        VecW agg_hethet = het1 & het2;
        VecW agg_het2hom1 = hom1 & het2;
        VecW agg_het1hom2 = hom2 & het1;
        agg_homhom = agg_homhom - (vecw_srli(agg_homhom, 1) & m1);
        agg_ibs0 = agg_ibs0 - (vecw_srli(agg_ibs0, 1) & m1);
        agg_hethet = agg_hethet - (vecw_srli(agg_hethet, 1) & m1);
        agg_het2hom1 = agg_het2hom1 - (vecw_srli(agg_het2hom1, 1) & m1);
        agg_het1hom2 = agg_het1hom2 - (vecw_srli(agg_het1hom2, 1) & m1);
        agg_homhom = (agg_homhom & m2) + (vecw_srli(agg_homhom, 2) & m2);
        agg_ibs0 = (agg_ibs0 & m2) + (vecw_srli(agg_ibs0, 2) & m2);
        agg_hethet = (agg_hethet & m2) + (vecw_srli(agg_hethet, 2) & m2);
        agg_het2hom1 = (agg_het2hom1 & m2) + (vecw_srli(agg_het2hom1, 2) & m2);
        agg_het1hom2 = (agg_het1hom2 & m2) + (vecw_srli(agg_het1hom2, 2) & m2);

        for (uint32_t offset = 1; offset != 3; ++offset) {
          hom1 = first_hom_iter[vec_idx + offset];
          hom2 = second_hom[vec_idx + offset];
          ref2het1 = first_ref2het_iter[vec_idx + offset];
          ref2het2 = second_ref2het[vec_idx + offset];
          VecW cur_homhom = hom1 & hom2;
          het1 = vecw_and_notfirst(hom1, ref2het1);
          het2 = vecw_and_notfirst(hom2, ref2het2);
          VecW cur_ibs0 = (ref2het1 ^ ref2het2) & cur_homhom;
          VecW cur_hethet = het1 & het2;
          VecW cur_het2hom1 = hom1 & het2;
          VecW cur_het1hom2 = hom2 & het1;
          cur_homhom = cur_homhom - (vecw_srli(cur_homhom, 1) & m1);
          cur_ibs0 = cur_ibs0 - (vecw_srli(cur_ibs0, 1) & m1);
          cur_hethet = cur_hethet - (vecw_srli(cur_hethet, 1) & m1);
          cur_het2hom1 = cur_het2hom1 - (vecw_srli(cur_het2hom1, 1) & m1);
          cur_het1hom2 = cur_het1hom2 - (vecw_srli(cur_het1hom2, 1) & m1);
          agg_homhom += (cur_homhom & m2) + (vecw_srli(cur_homhom, 2) & m2);
          agg_ibs0 += (cur_ibs0 & m2) + (vecw_srli(cur_ibs0, 2) & m2);
          agg_hethet += (cur_hethet & m2) + (vecw_srli(cur_hethet, 2) & m2);
          agg_het2hom1 += (cur_het2hom1 & m2) + (vecw_srli(cur_het2hom1, 2) & m2);
          agg_het1hom2 += (cur_het1hom2 & m2) + (vecw_srli(cur_het1hom2, 2) & m2);
        }
        acc_homhom.vw = acc_homhom.vw + (agg_homhom & m4) + (vecw_srli(agg_homhom, 4) & m4);
        acc_ibs0.vw = acc_ibs0.vw + (agg_ibs0 & m4) + (vecw_srli(agg_ibs0, 4) & m4);
        acc_hethet.vw = acc_hethet.vw + (agg_hethet & m4) + (vecw_srli(agg_hethet, 4) & m4);
        acc_het2hom1.vw = acc_het2hom1.vw + (agg_het2hom1 & m4) + (vecw_srli(agg_het2hom1, 4) & m4);
        acc_het1hom2.vw = acc_het1hom2.vw + (agg_het1hom2 & m4) + (vecw_srli(agg_het1hom2, 4) & m4);
      }
      const VecW m8 = VCONST_W(kMask00FF);
      acc_homhom.vw = (acc_homhom.vw & m8) + (vecw_srli(acc_homhom.vw, 8) & m8);
      acc_ibs0.vw = (acc_ibs0.vw & m8) + (vecw_srli(acc_ibs0.vw, 8) & m8);
      acc_hethet.vw = (acc_hethet.vw & m8) + (vecw_srli(acc_hethet.vw, 8) & m8);
      acc_het2hom1.vw = (acc_het2hom1.vw & m8) + (vecw_srli(acc_het2hom1.vw, 8) & m8);
      acc_het1hom2.vw = (acc_het1hom2.vw & m8) + (vecw_srli(acc_het1hom2.vw, 8) & m8);
      king_counts_iter[kKingOffsetIbs0] += UniVecHsum16(acc_ibs0);
      king_counts_iter[kKingOffsetHethet] += UniVecHsum16(acc_hethet);
      king_counts_iter[kKingOffsetHet2Hom1] += UniVecHsum16(acc_het2hom1);
      king_counts_iter[kKingOffsetHet1Hom2] += UniVecHsum16(acc_het1hom2);
      king_counts_iter[kKingOffsetHomhom] += UniVecHsum16(acc_homhom);
      king_counts_iter = &(king_counts_iter[5]);

      first_hom_iter = &(first_hom_iter[kKingMultiplexVecs]);
      first_ref2het_iter = &(first_ref2het_iter[kKingMultiplexVecs]);
    }
  }
}
#endif
static_assert(!(kKingMultiplexWords % 2), "kKingMultiplexWords must be even for safe bit-transpose.");

typedef struct CalcKingDenseCtxStruct {
  uintptr_t* smaj_hom[2];
  uintptr_t* smaj_ref2het[2];
  uint32_t homhom_needed;

  uint32_t* thread_start;

  uint32_t* king_counts;
} CalcKingDenseCtx;

THREAD_FUNC_DECL CalcKingDenseThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcKingDenseCtx* ctx = S_CAST(CalcKingDenseCtx*, arg->sharedp->context);

  const uint64_t mem_start_idx = ctx->thread_start[0];
  const uint64_t start_idx = ctx->thread_start[tidx];
  const uint32_t end_idx = ctx->thread_start[tidx + 1];
  const uint32_t homhom_needed = ctx->homhom_needed;
  uint32_t parity = 0;
  do {
    if (homhom_needed) {
      IncrKingHomhom(ctx->smaj_hom[parity], ctx->smaj_ref2het[parity], start_idx, end_idx, &(ctx->king_counts[((start_idx * (start_idx - 1) - mem_start_idx * (mem_start_idx - 1)) / 2) * 5]));
    } else {
      IncrKing(ctx->smaj_hom[parity], ctx->smaj_ref2het[parity], start_idx, end_idx, &(ctx->king_counts[(start_idx * (start_idx - 1) - mem_start_idx * (mem_start_idx - 1)) * 2]));
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

/*
double ComputeKinship(const uint32_t* king_counts_entry) {
  const uint32_t ibs0_ct = king_counts_entry[kKingOffsetIbs0];
  const uint32_t hethet_ct = king_counts_entry[kKingOffsetHethet];
  const uint32_t het2hom1_ct = king_counts_entry[kKingOffsetHet2Hom1];
  const uint32_t het1hom2_ct = king_counts_entry[kKingOffsetHet1Hom2];
  const intptr_t smaller_het_ct = hethet_ct + MINV(het1hom2_ct, het2hom1_ct);
  return 0.5 - (S_CAST(double, 4 * S_CAST(intptr_t, ibs0_ct) + het1hom2_ct + het2hom1_ct) / S_CAST(double, 4 * smaller_het_ct));
}
*/

// '2' refers to the larger index here
double ComputeKinship(const uint32_t* king_counts_entry, uint32_t singleton_het1_ct, uint32_t singleton_hom1_ct, uint32_t singleton_het2_ct, uint32_t singleton_hom2_ct) {
  const uint32_t ibs0_ct = king_counts_entry[kKingOffsetIbs0] + singleton_hom1_ct + singleton_hom2_ct;
  const uint32_t hethet_ct = king_counts_entry[kKingOffsetHethet];
  const uint32_t het2hom1_ct = king_counts_entry[kKingOffsetHet2Hom1] + singleton_het2_ct;
  const uint32_t het1hom2_ct = king_counts_entry[kKingOffsetHet1Hom2] + singleton_het1_ct;
  const intptr_t smaller_het_ct = hethet_ct + MINV(het1hom2_ct, het2hom1_ct);
  return 0.5 - (S_CAST(double, 4 * S_CAST(intptr_t, ibs0_ct) + het1hom2_ct + het2hom1_ct) / S_CAST(double, 4 * smaller_het_ct));
}

// could also return pointer to end?
void SetKingMatrixFname(KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, char* outname_end) {
  if (!(king_flags & (kfKingMatrixBin | kfKingMatrixBin4))) {
    char* outname_end2 = strcpya_k(outname_end, ".king");
    const uint32_t output_zst = king_flags & kfKingMatrixZs;
    if (parallel_tot != 1) {
      *outname_end2++ = '.';
      outname_end2 = u32toa(parallel_idx + 1, outname_end2);
    }
    if (output_zst) {
      outname_end2 = strcpya_k(outname_end2, ".zst");
    }
    *outname_end2 = '\0';
    return;
  }
  char* outname_end2 = strcpya_k(outname_end, ".king.bin");
  if (parallel_tot != 1) {
    *outname_end2++ = '.';
    outname_end2 = u32toa(parallel_idx + 1, outname_end2);
  }
  *outname_end2 = '\0';
}

void SetKingTableFname(KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, char* outname_end) {
  char* outname_end2 = strcpya_k(outname_end, ".kin0");
  const uint32_t output_zst = king_flags & kfKingTableZs;
  if (parallel_tot != 1) {
    *outname_end2++ = '.';
    outname_end2 = u32toa(parallel_idx + 1, outname_end2);
  }
  if (output_zst) {
    outname_end2 = strcpya_k(outname_end2, ".zst");
  }
  *outname_end2 = '\0';
}

char* AppendKingTableHeader(KingFlags king_flags, uint32_t king_col_fid, uint32_t king_col_sid, char* cswritep) {
  *cswritep++ = '#';
  if (king_flags & kfKingColId) {
    if (king_col_fid) {
      cswritep = strcpya_k(cswritep, "FID1\t");
    }
    // Was 'ID1' before alpha 3, but that's inconsistent with other plink2
    // commands, and in the meantime the header line still doesn't perfectly
    // match KING due to e.g. capitalization.
    cswritep = strcpya_k(cswritep, "IID1\t");
    if (king_col_sid) {
      cswritep = strcpya_k(cswritep, "SID1\t");
    }
    if (king_col_fid) {
      cswritep = strcpya_k(cswritep, "FID2\t");
    }
    cswritep = strcpya_k(cswritep, "IID2\t");
    if (king_col_sid) {
      cswritep = strcpya_k(cswritep, "SID2\t");
    }
  }
  if (king_flags & kfKingColNsnp) {
    cswritep = strcpya_k(cswritep, "NSNP\t");
  }
  if (king_flags & kfKingColHethet) {
    cswritep = strcpya_k(cswritep, "HETHET\t");
  }
  if (king_flags & kfKingColIbs0) {
    cswritep = strcpya_k(cswritep, "IBS0\t");
  }
  if (king_flags & kfKingColIbs1) {
    cswritep = strcpya_k(cswritep, "HET1_HOM2\tHET2_HOM1\t");
  }
  if (king_flags & kfKingColHamming) {
    cswritep = strcpya_k(cswritep, "IBS\t");
  }
  if (king_flags & kfKingColKinship) {
    cswritep = strcpya_k(cswritep, "KINSHIP\t");
  }
  DecrAppendBinaryEoln(&cswritep);
  return cswritep;
}

uint32_t KingMaxSparseCt(uint32_t row_end_idx) {
#ifdef USE_AVX2
  return row_end_idx / 33;
#else
  return row_end_idx / 30;
#endif
}

PglErr CalcKing(const SampleIdInfo* siip, const uintptr_t* variant_include_orig, const ChrInfo* cip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, double king_cutoff, double king_table_filter, KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, PgenReader* simple_pgrp, uintptr_t* sample_include, uint32_t* sample_ct_ptr, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  char* cswritetp = nullptr;
  CompressStreamState css;
  CompressStreamState csst;
  ThreadGroup tg;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  PreinitCstream(&csst);
  PreinitThreads(&tg);
  {
    const KingFlags matrix_shape = king_flags & kfKingMatrixShapemask;
    const char* flagname = matrix_shape? "--make-king" : ((king_flags & kfKingColAll)? "--make-king-table" : "--king-cutoff");
    if (unlikely(IsSet(cip->haploid_mask, 0))) {
      logerrprintf("Error: %s cannot be used on haploid genomes.\n", flagname);
      goto CalcKing_ret_INCONSISTENT_INPUT;
    }
    uint32_t sample_ct = *sample_ct_ptr;
    if (unlikely(sample_ct < 2)) {
      logerrprintf("Error: %s requires at least 2 samples.\n", flagname);
      goto CalcKing_ret_DEGENERATE_DATA;
    }
#ifdef __LP64__
    // there's also a UINT32_MAX / kKingMultiplexWords limit, but that's not
    // relevant for now
    if (unlikely(sample_ct > 134000000)) {
      // for text output, 134m * 16 is just below kMaxLongLine
      logerrprintf("Error: %s does not support > 134000000 samples.\n", flagname);
      reterr = kPglRetNotYetSupported;
      goto CalcKing_ret_1;
    }
#endif
    const uintptr_t sample_ctl = BitCtToWordCt(sample_ct);
    uintptr_t* kinship_table = nullptr;
    if (king_cutoff != -1) {
      if (unlikely(bigstack_calloc_w(sample_ct * sample_ctl, &kinship_table))) {
        goto CalcKing_ret_NOMEM;
      }
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uint32_t non_autosomal_variant_ct = CountNonAutosomalVariants(variant_include_orig, cip, 1, 1);
    if (non_autosomal_variant_ct) {
      uintptr_t* variant_include_next;
      if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include_next))) {
        goto CalcKing_ret_NOMEM;
      }
      logprintf("Excluding %u variant%s on non-autosomes from KING-robust calculation.\n", non_autosomal_variant_ct, (non_autosomal_variant_ct == 1)? "" : "s");
      variant_ct -= non_autosomal_variant_ct;
      if (!variant_ct) {
        logerrprintf("Error: No variants remaining for KING-robust calculation.\n");
        goto CalcKing_ret_DEGENERATE_DATA;
      }
      memcpy(variant_include_next, variant_include_orig, raw_variant_ctl * sizeof(intptr_t));
      ExcludeNonAutosomalVariants(cip, variant_include_next);
      variant_include_orig = variant_include_next;
    }
    uintptr_t* variant_include;

    if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include))) {
      goto CalcKing_ret_NOMEM;
    }

    uint32_t grand_row_start_idx;
    uint32_t grand_row_end_idx;
    ParallelBounds(sample_ct, 1, parallel_idx, parallel_tot, R_CAST(int32_t*, &grand_row_start_idx), R_CAST(int32_t*, &grand_row_end_idx));

    // possible todo: allow this to change between passes
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > sample_ct / 32) {
      calc_thread_ct = sample_ct / 32;
    }
    if (!calc_thread_ct) {
      calc_thread_ct = 1;
    }
    const uint32_t homhom_needed = (king_flags & kfKingColNsnp) || ((!(king_flags & kfKingCounts)) && (king_flags & (kfKingColHethet | kfKingColIbs0 | kfKingColIbs1 | kfKingColHamming)));
    CalcKingSparseCtx sparse_ctx;
    uint32_t sparse_read_block_size = 0;
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    // These values are now permitted to underflow from sparse-optimization.
    // Might want to change them to int32_t*.
    uint32_t* singleton_het_cts;
    uint32_t* singleton_hom_cts;
    uint32_t* singleton_missing_cts;
    {
      sparse_ctx.variant_include_orig = variant_include_orig;
      sparse_ctx.homhom_needed = homhom_needed;
      const uint32_t max_sparse_ct = KingMaxSparseCt(grand_row_end_idx);
      // Ok for this to be a slight underestimate, since bigstack_left()/8 is
      // an arbitrary limit anyway.
      const uintptr_t thread_xalloc_cacheline_ct = DivUp((3 * k1LU) * (max_sparse_ct + grand_row_end_idx), kInt32PerCacheline) + ((kPglVblockSize * 2) / kBitsPerCacheline);
      if (unlikely(PgenMtLoadInit(variant_include_orig, grand_row_end_idx, variant_ct, bigstack_left() / 8, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &sparse_ctx.genovecs, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &sparse_read_block_size, nullptr, main_loadbufs, &sparse_ctx.pgr_ptrs, &sparse_ctx.read_variant_uidx_starts))) {
        goto CalcKing_ret_NOMEM;
      }
      sparse_ctx.read_block_size = sparse_read_block_size;
      sparse_ctx.err_info = (~0LLU) << 32;
      if (unlikely(bigstack_alloc_u32p(calc_thread_ct, &sparse_ctx.thread_idx_bufs) ||
                   bigstack_alloc_u32p(calc_thread_ct, &sparse_ctx.thread_singleton_het_cts) ||
                   bigstack_alloc_u32p(calc_thread_ct, &sparse_ctx.thread_singleton_hom_cts) ||
                   bigstack_alloc_u32p(calc_thread_ct, &sparse_ctx.thread_singleton_missing_cts) ||
                   bigstack_alloc_u32(calc_thread_ct, &sparse_ctx.thread_skip_cts) ||
                   bigstack_alloc_wp(calc_thread_ct, &sparse_ctx.thread_sparse_excludes[0]) ||
                   bigstack_alloc_wp(calc_thread_ct, &sparse_ctx.thread_sparse_excludes[1]))) {
        goto CalcKing_ret_NOMEM;
      }
      const uint32_t read_block_sizel = sparse_read_block_size / kBitsPerWord;
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        if (unlikely(bigstack_alloc_u32(3 * max_sparse_ct, &(sparse_ctx.thread_idx_bufs[tidx])) ||
                     bigstack_alloc_u32(grand_row_end_idx, &(sparse_ctx.thread_singleton_het_cts[tidx])) ||
                     bigstack_alloc_u32(grand_row_end_idx, &(sparse_ctx.thread_singleton_hom_cts[tidx])) ||
                     bigstack_alloc_u32(grand_row_end_idx, &(sparse_ctx.thread_singleton_missing_cts[tidx])) ||
                     bigstack_alloc_w(read_block_sizel, &(sparse_ctx.thread_sparse_excludes[0][tidx])) ||
                     bigstack_alloc_w(read_block_sizel, &(sparse_ctx.thread_sparse_excludes[1][tidx])))) {
          goto CalcKing_ret_NOMEM;
        }
      }
      singleton_het_cts = sparse_ctx.thread_singleton_het_cts[0];
      singleton_hom_cts = sparse_ctx.thread_singleton_hom_cts[0];
      singleton_missing_cts = sparse_ctx.thread_singleton_missing_cts[0];
    }

    CalcKingDenseCtx dense_ctx;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg) ||
                 bigstack_alloc_u32(calc_thread_ct + 1, &dense_ctx.thread_start))) {
      goto CalcKing_ret_NOMEM;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t grei_ctaw = BitCtToAlignedWordCt(grand_row_end_idx);
    const uint32_t grei_ctaw2 = NypCtToAlignedWordCt(grand_row_end_idx);
    dense_ctx.homhom_needed = homhom_needed;
    const uint32_t king_bufsizew = kKingMultiplexWords * grand_row_end_idx;
    const uint32_t homhom_needed_p4 = dense_ctx.homhom_needed + 4;
    uintptr_t* cur_sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* loadbuf;
    uintptr_t* splitbuf_hom;
    uintptr_t* splitbuf_ref2het;
    VecW* vecaligned_buf;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &cur_sample_include) ||
                 bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_w(grei_ctaw2, &loadbuf) ||
                 bigstack_alloc_w(kPglBitTransposeBatch * grei_ctaw, &splitbuf_hom) ||
                 bigstack_alloc_w(kPglBitTransposeBatch * grei_ctaw, &splitbuf_ref2het) ||
                 bigstack_alloc_w(king_bufsizew, &(dense_ctx.smaj_hom[0])) ||
                 bigstack_alloc_w(king_bufsizew, &(dense_ctx.smaj_ref2het[0])) ||
                 bigstack_alloc_w(king_bufsizew, &(dense_ctx.smaj_hom[1])) ||
                 bigstack_alloc_w(king_bufsizew, &(dense_ctx.smaj_ref2het[1])) ||
                 bigstack_alloc_v(kPglBitTransposeBufvecs, &vecaligned_buf))) {
      goto CalcKing_ret_NOMEM;
    }

    // Make this automatically multipass when there's insufficient memory.  So
    // we open the output file(s) here, and just append in the main loop.
    unsigned char* numbuf = nullptr;
    if (matrix_shape) {
      SetKingMatrixFname(king_flags, parallel_idx, parallel_tot, outname_end);
      if (!(king_flags & (kfKingMatrixBin | kfKingMatrixBin4))) {
        // text matrix
        // won't be >2gb since sample_ct <= 134m
        const uint32_t overflow_buf_size = kCompressStreamBlock + 16 * sample_ct;
        reterr = InitCstreamAlloc(outname, 0, king_flags & kfKingMatrixZs, max_thread_ct, overflow_buf_size, &css, &cswritep);
        if (unlikely(reterr)) {
          goto CalcKing_ret_1;
        }
      } else {
        if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
          goto CalcKing_ret_OPEN_FAIL;
        }
        if (unlikely(bigstack_alloc_uc(sample_ct * 4 * (2 - ((king_flags / kfKingMatrixBin4) & 1)), &numbuf))) {
          goto CalcKing_ret_OPEN_FAIL;
        }
      }
    }
    uint32_t king_col_fid = 0;
    uint32_t king_col_sid = 0;
    uintptr_t max_sample_fmtid_blen = 0;
    char* collapsed_sample_fmtids = nullptr;
    if (king_flags & kfKingColAll) {
      const uint32_t overflow_buf_size = kCompressStreamBlock + kMaxMediumLine;
      SetKingTableFname(king_flags, parallel_idx, parallel_tot, outname_end);
      reterr = InitCstreamAlloc(outname, 0, king_flags & kfKingTableZs, max_thread_ct, overflow_buf_size, &csst, &cswritetp);
      if (unlikely(reterr)) {
        goto CalcKing_ret_1;
      }

      king_col_fid = FidColIsRequired(siip, king_flags / kfKingColMaybefid);
      king_col_sid = SidColIsRequired(siip->sids, king_flags / kfKingColMaybesid);
      if (!parallel_idx) {
        cswritetp = AppendKingTableHeader(king_flags, king_col_fid, king_col_sid, cswritetp);
      }
      if (unlikely(CollapsedSampleFmtidInitAlloc(sample_include, siip, grand_row_end_idx, king_col_fid, king_col_sid, &collapsed_sample_fmtids, &max_sample_fmtid_blen))) {
        goto CalcKing_ret_NOMEM;
      }
    }
    uint64_t king_table_filter_ct = 0;
    const uintptr_t cells_avail = bigstack_left() / (sizeof(int32_t) * homhom_needed_p4);
    const uint32_t pass_ct = CountTrianglePasses(grand_row_start_idx, grand_row_end_idx, 1, cells_avail);
    if (unlikely(!pass_ct)) {
      goto CalcKing_ret_NOMEM;
    }
    if (unlikely((pass_ct > 1) && (king_flags & kfKingMatrixSq))) {
      logerrputs("Insufficient memory for --make-king square output.  Try square0 or triangle\nshape instead.\n");
      goto CalcKing_ret_NOMEM;
    }
    uint32_t row_end_idx = grand_row_start_idx;
    sparse_ctx.king_counts = R_CAST(uint32_t*, g_bigstack_base);
    dense_ctx.king_counts = sparse_ctx.king_counts;
    for (uint32_t pass_idx_p1 = 1; pass_idx_p1 <= pass_ct; ++pass_idx_p1) {
      const uint32_t row_start_idx = row_end_idx;
      row_end_idx = NextTrianglePass(row_start_idx, grand_row_end_idx, 1, cells_avail);
      TriangleLoadBalance(calc_thread_ct, row_start_idx, row_end_idx, 1, dense_ctx.thread_start);
      memcpy(cur_sample_include, sample_include, raw_sample_ctl * sizeof(intptr_t));
      // bugfix (20 Nov 2019): forgot that --parallel could cause the old
      // row_end_idx != grand_row_end_idx comparison not work
      if (row_end_idx != orig_sample_ct) {
        uint32_t sample_uidx_end = IdxToUidxBasic(sample_include, row_end_idx);
        ClearBitsNz(sample_uidx_end, raw_sample_ct, cur_sample_include);
      }
      FillCumulativePopcounts(cur_sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      pgfip->block_base = main_loadbufs[0];  // needed after first pass
      // Update (9 Nov 2019): The one-time singleton/monomorphic scan has been
      // replaced with a more effective sparse-variant scan which happens on
      // every pass.
      sparse_ctx.sample_include = cur_sample_include;
      sparse_ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
      sparse_ctx.row_start_idx = row_start_idx;
      sparse_ctx.row_end_idx = row_end_idx;
      sparse_ctx.max_sparse_ct = KingMaxSparseCt(row_end_idx);
      logprintf("%s pass %u/%u: Scanning for rare variants... ", flagname, pass_idx_p1, pass_ct);
      fputs("0%", stdout);
      fflush(stdout);
      SetThreadFuncAndData(CalcKingSparseThread, &sparse_ctx, &tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto CalcKing_ret_THREAD_CREATE_FAIL;
      }
      memcpy(variant_include, variant_include_orig, raw_variant_ctl * sizeof(intptr_t));

      {
        const uint32_t read_block_sizel = sparse_read_block_size / kBitsPerWord;
        uint32_t prev_read_block_idx = 0;
        uint32_t read_block_idx = 0;
        uint32_t pct = 0;
        uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
        uint32_t parity = 0;
        for (uint32_t variant_idx = 0; ; ) {
          const uint32_t cur_block_size = MultireadNonempty(variant_include_orig, &tg, raw_variant_ct, sparse_read_block_size, pgfip, &read_block_idx, &reterr);
          if (unlikely(reterr)) {
            goto CalcKing_ret_PGR_FAIL;
          }
          JoinThreads(&tg);
          reterr = S_CAST(PglErr, sparse_ctx.err_info);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, sparse_ctx.err_info >> 32);
            goto CalcKing_ret_1;
          }
          if (!IsLastBlock(&tg)) {
            sparse_ctx.cur_block_size = cur_block_size;
            ComputeUidxStartPartition(variant_include_orig, cur_block_size, calc_thread_ct, read_block_idx * sparse_read_block_size, sparse_ctx.read_variant_uidx_starts);
            PgrCopyBaseAndOffset(pgfip, calc_thread_ct, sparse_ctx.pgr_ptrs);
            if (variant_idx + cur_block_size == variant_ct) {
              DeclareLastThreadBlock(&tg);
            }
            SpawnThreads(&tg);
          }
          parity = 1 - parity;
          if (variant_idx) {
            uintptr_t* variant_include_update = &(variant_include[prev_read_block_idx * read_block_sizel]);
            for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
              BitvecInvmask(sparse_ctx.thread_sparse_excludes[parity][tidx], read_block_sizel, variant_include_update);
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
              next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
            }
          }
          prev_read_block_idx = read_block_idx;
          ++read_block_idx;
          variant_idx += cur_block_size;
          pgfip->block_base = main_loadbufs[parity];
        }
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
      }
      fputs("\b\b", stdout);
      logputs("done.\n");
      const uint32_t cur_variant_ct = PopcountWords(variant_include, raw_variant_ctl);
      uint32_t sparse_variant_ct = variant_ct - cur_variant_ct;
      logprintf("%u variant%s handled by initial scan (%u remaining).\n", sparse_variant_ct, (sparse_variant_ct == 1)? "" : "s", cur_variant_ct);
      uint32_t skip_ct = sparse_ctx.thread_skip_cts[0];
      const uint32_t vec_ct = DivUp(row_end_idx, kInt32PerVec);
      for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
        U32CastVecAdd(sparse_ctx.thread_singleton_het_cts[tidx], vec_ct, singleton_het_cts);
        U32CastVecAdd(sparse_ctx.thread_singleton_hom_cts[tidx], vec_ct, singleton_hom_cts);
        U32CastVecAdd(sparse_ctx.thread_singleton_missing_cts[tidx], vec_ct, singleton_missing_cts);
        skip_ct += sparse_ctx.thread_skip_cts[tidx];
      }
      sparse_variant_ct -= skip_ct;
      if (cur_variant_ct) {
        SetThreadFuncAndData(CalcKingDenseThread, &dense_ctx, &tg);
        const uint32_t row_end_idxaw = BitCtToAlignedWordCt(row_end_idx);
        const uint32_t row_end_idxaw2 = NypCtToAlignedWordCt(row_end_idx);
        if (row_end_idxaw % 2) {
          const uint32_t cur_king_bufsizew = kKingMultiplexWords * row_end_idx;
          uintptr_t* smaj_hom0_last = &(dense_ctx.smaj_hom[0][kKingMultiplexWords - 1]);
          uintptr_t* smaj_ref2het0_last = &(dense_ctx.smaj_ref2het[0][kKingMultiplexWords - 1]);
          uintptr_t* smaj_hom1_last = &(dense_ctx.smaj_hom[1][kKingMultiplexWords - 1]);
          uintptr_t* smaj_ref2het1_last = &(dense_ctx.smaj_ref2het[1][kKingMultiplexWords - 1]);
          for (uint32_t offset = 0; offset < cur_king_bufsizew; offset += kKingMultiplexWords) {
            smaj_hom0_last[offset] = 0;
            smaj_ref2het0_last[offset] = 0;
            smaj_hom1_last[offset] = 0;
            smaj_ref2het1_last[offset] = 0;
          }
        }
        uintptr_t variant_uidx_base = 0;
        uintptr_t cur_bits = variant_include[0];
        uint32_t variants_completed = 0;
        uint32_t parity = 0;
        const uint32_t sample_batch_ct_m1 = (row_end_idx - 1) / kPglBitTransposeBatch;
        // Similar to plink 1.9 --genome.  For each pair of samples S1-S2, we
        // need to determine counts of the following:
        //   * S1 hom-S2 opposite hom
        //   * het-het
        //   * S1 hom-S2 het
        //   * S2 hom-S1 het
        //   * sometimes S1 hom-S2 same hom
        //   * (nonmissing determined via subtraction)
        // We handle this as follows:
        //   1. set n=0, reader thread loads first kKingMultiplex variants and
        //      converts+transposes the data to a sample-major format suitable
        //      for multithreaded computation.
        //   2. spawn threads
        //
        //   3. increment n by 1
        //   4. load block n unless eof
        //   5. permit threads to continue to next block, unless eof
        //   6. goto step 3 unless eof
        //
        //   7. write results
        // Results are always reported in lower-triangular order, rather than
        // KING's upper-triangular order, since the former plays more nicely
        // with incremental addition of samples.
        PgrSampleSubsetIndex pssi;
        PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
        do {
          const uint32_t cur_block_size = MINV(cur_variant_ct - variants_completed, kKingMultiplex);
          uintptr_t* cur_smaj_hom = dense_ctx.smaj_hom[parity];
          uintptr_t* cur_smaj_ref2het = dense_ctx.smaj_ref2het[parity];
          // "block" = distance computation granularity, usually 1024 or 1536
          //           variants
          // "batch" = variant-major-to-sample-major transpose granularity,
          //           currently 512 variants
          uint32_t variant_batch_size = kPglBitTransposeBatch;
          uint32_t variant_batch_size_rounded_up = kPglBitTransposeBatch;
          const uint32_t write_batch_ct_m1 = (cur_block_size - 1) / kPglBitTransposeBatch;
          for (uint32_t write_batch_idx = 0; ; ++write_batch_idx) {
            if (write_batch_idx >= write_batch_ct_m1) {
              if (write_batch_idx > write_batch_ct_m1) {
                break;
              }
              variant_batch_size = ModNz(cur_block_size, kPglBitTransposeBatch);
              variant_batch_size_rounded_up = variant_batch_size;
              const uint32_t variant_batch_size_rem = variant_batch_size % kBitsPerWord;
              if (variant_batch_size_rem) {
                const uint32_t trailing_variant_ct = kBitsPerWord - variant_batch_size_rem;
                variant_batch_size_rounded_up += trailing_variant_ct;
                ZeroWArr(trailing_variant_ct * row_end_idxaw, &(splitbuf_hom[variant_batch_size * row_end_idxaw]));
                ZeroWArr(trailing_variant_ct * row_end_idxaw, &(splitbuf_ref2het[variant_batch_size * row_end_idxaw]));
              }
            }
            uintptr_t* hom_iter = splitbuf_hom;
            uintptr_t* ref2het_iter = splitbuf_ref2het;
            for (uint32_t uii = 0; uii != variant_batch_size; ++uii) {
              const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
              // Model does not cleanly generalize to multiallelic variants
              // (unless there's something I overlooked, which is quite
              // possible).
              // Thought about using major allele counts in that case, but that
              // sacrifices a really nice property of this method: estimated
              // relationship coefficient between each pair of samples is
              // independent of estimated allele frequencies.  And the accuracy
              // improvement we'd get in return is microscopic.  So we stick to
              // REF/ALT allele counts instead.
              reterr = PgrGet(cur_sample_include, pssi, row_end_idx, variant_uidx, simple_pgrp, loadbuf);
              if (unlikely(reterr)) {
                PgenErrPrintNV(reterr, variant_uidx);
                goto CalcKing_ret_1;
              }
              SetTrailingNyps(row_end_idx, loadbuf);
              SplitHomRef2hetUnsafeW(loadbuf, row_end_idxaw2, hom_iter, ref2het_iter);
              hom_iter = &(hom_iter[row_end_idxaw]);
              ref2het_iter = &(ref2het_iter[row_end_idxaw]);
            }
            // uintptr_t* read_iter = loadbuf;
            uintptr_t* write_hom_iter = &(cur_smaj_hom[write_batch_idx * kPglBitTransposeWords]);
            uintptr_t* write_ref2het_iter = &(cur_smaj_ref2het[write_batch_idx * kPglBitTransposeWords]);
            uint32_t write_batch_size = kPglBitTransposeBatch;
            for (uint32_t sample_batch_idx = 0; ; ++sample_batch_idx) {
              if (sample_batch_idx >= sample_batch_ct_m1) {
                if (sample_batch_idx > sample_batch_ct_m1) {
                  break;
                }
                write_batch_size = ModNz(row_end_idx, kPglBitTransposeBatch);
              }
              // bugfix: read_batch_size must be rounded up to word boundary,
              // since we want to one-out instead of zero-out the trailing bits
              //
              // bugfix: if we always use kPglBitTransposeBatch instead of
              // variant_batch_size_rounded_up, we read/write past the
              // kKingMultiplex limit and clobber the first variants of the
              // next sample with garbage.
              TransposeBitblock(&(splitbuf_hom[sample_batch_idx * kPglBitTransposeWords]), row_end_idxaw, kKingMultiplexWords, variant_batch_size_rounded_up, write_batch_size, write_hom_iter, vecaligned_buf);
              TransposeBitblock(&(splitbuf_ref2het[sample_batch_idx * kPglBitTransposeWords]), row_end_idxaw, kKingMultiplexWords, variant_batch_size_rounded_up, write_batch_size, write_ref2het_iter, vecaligned_buf);
              write_hom_iter = &(write_hom_iter[kKingMultiplex * kPglBitTransposeWords]);
              write_ref2het_iter = &(write_ref2het_iter[kKingMultiplex * kPglBitTransposeWords]);
            }
          }
          const uint32_t cur_block_sizew = BitCtToWordCt(cur_block_size);
          if (cur_block_sizew < kKingMultiplexWords) {
            uintptr_t* write_hom_iter = &(cur_smaj_hom[cur_block_sizew]);
            uintptr_t* write_ref2het_iter = &(cur_smaj_ref2het[cur_block_sizew]);
            const uint32_t write_word_ct = kKingMultiplexWords - cur_block_sizew;
            for (uint32_t sample_idx = 0; sample_idx != row_end_idx; ++sample_idx) {
              ZeroWArr(write_word_ct, write_hom_iter);
              ZeroWArr(write_word_ct, write_ref2het_iter);
              write_hom_iter = &(write_hom_iter[kKingMultiplexWords]);
              write_ref2het_iter = &(write_ref2het_iter[kKingMultiplexWords]);
            }
          }
          if (variants_completed) {
            JoinThreads(&tg);
            // CalcKingThread() never errors out
          }
          // this update must occur after JoinThreads() call
          if (variants_completed + cur_block_size == cur_variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto CalcKing_ret_THREAD_CREATE_FAIL;
          }
          printf("\r%s pass %u/%u: %u variants complete.", flagname, pass_idx_p1, pass_ct, variants_completed);
          fflush(stdout);
          variants_completed += cur_block_size;
          parity = 1 - parity;
        } while (!IsLastBlock(&tg));
        JoinThreads(&tg);
      }
      if (matrix_shape || (king_flags & kfKingColAll)) {
        printf("\r%s pass %u/%u: Writing...                   \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", flagname, pass_idx_p1, pass_ct);
        fflush(stdout);
        // allow simultaneous --make-king + --make-king-table
        if (matrix_shape) {
          if (!(king_flags & (kfKingMatrixBin | kfKingMatrixBin4))) {
            const uint32_t is_squarex = king_flags & (kfKingMatrixSq | kfKingMatrixSq0);
            const uint32_t is_square0 = king_flags & kfKingMatrixSq0;
            uint32_t* results_iter = dense_ctx.king_counts;
            uint32_t sample_idx1 = row_start_idx;
            if (is_squarex && (!parallel_idx) && (pass_idx_p1)) {
              // dump "empty" first row
              sample_idx1 = 0;
            }
            for (; sample_idx1 != row_end_idx; ++sample_idx1) {
              const uint32_t singleton_het1_ct = singleton_het_cts[sample_idx1];
              const uint32_t singleton_hom1_ct = singleton_hom_cts[sample_idx1];
              for (uint32_t sample_idx2 = 0; sample_idx2 < sample_idx1; ++sample_idx2) {
                const double kinship_coeff = ComputeKinship(results_iter, singleton_het_cts[sample_idx2], singleton_hom_cts[sample_idx2], singleton_het1_ct, singleton_hom1_ct);
                if (kinship_table && (kinship_coeff > king_cutoff)) {
                  SetBit(sample_idx2, &(kinship_table[sample_idx1 * sample_ctl]));
                  SetBit(sample_idx1, &(kinship_table[sample_idx2 * sample_ctl]));
                }
                cswritep = dtoa_g(kinship_coeff, cswritep);
                *cswritep++ = '\t';
                results_iter = &(results_iter[homhom_needed_p4]);
              }
              if (is_squarex) {
                cswritep = strcpya_k(cswritep, "0.5");
                if (is_square0) {
                  // (roughly same performance as creating a tab-zero constant
                  // buffer in advance)
                  const uint32_t zcount = sample_ct - sample_idx1 - 1;
                  const uint32_t wct = DivUp(zcount, kBytesPerWord / 2);
                  // assumes little-endian
                  const uintptr_t tabzero_word = 0x3009 * kMask0001;
                  unsigned char* tabzero_write_iter = R_CAST(unsigned char*, cswritep);
                  for (uintptr_t widx = 0; widx != wct; ++widx) {
                    AppendW(tabzero_word, &tabzero_write_iter);
                  }
                  cswritep = &(cswritep[zcount * 2]);
                } else {
                  const uint32_t* results_iter2 = &(results_iter[sample_idx1 * homhom_needed_p4]);
                  // 0
                  // 1  2
                  // 3  4  5
                  // 6  7  8  9
                  // 10 11 12 13 14

                  // sample_idx1 = 0: [0] 0 1 3 6 10...
                  // sample_idx1 = 1: [1] 2 4 7 11...
                  // sample_idx1 = 2: [3] 5 8 12...
                  // sample_idx1 = 3: [6] 9 13...
                  for (uint32_t sample_idx2 = sample_idx1 + 1; sample_idx2 != sample_ct; ++sample_idx2) {
                    *cswritep++ = '\t';
                    cswritep = dtoa_g(ComputeKinship(results_iter2, singleton_het1_ct, singleton_hom1_ct, singleton_het_cts[sample_idx2], singleton_hom_cts[sample_idx2]), cswritep);
                    results_iter2 = &(results_iter2[sample_idx2 * homhom_needed_p4]);
                  }
                }
                ++cswritep;
              }
              DecrAppendBinaryEoln(&cswritep);
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto CalcKing_ret_WRITE_FAIL;
              }
            }
          } else {
            // binary matrix output
            // er, probably want to revise this so there's less duplicated code
            // from text matrix output...
            const uint32_t is_squarex = king_flags & (kfKingMatrixSq | kfKingMatrixSq0);
            const uint32_t is_square0 = king_flags & kfKingMatrixSq0;
            uint32_t* results_iter = dense_ctx.king_counts;
            uint32_t sample_idx1 = row_start_idx;
            if (is_squarex && (!parallel_idx)) {
              sample_idx1 = 0;
            }
            if (king_flags & kfKingMatrixBin4) {
              float* write_row = R_CAST(float*, numbuf);
              uintptr_t row_byte_ct = sample_ct * sizeof(float);
              for (; sample_idx1 != row_end_idx; ++sample_idx1) {
                const uint32_t singleton_het1_ct = singleton_het_cts[sample_idx1];
                const uint32_t singleton_hom1_ct = singleton_hom_cts[sample_idx1];
                for (uint32_t sample_idx2 = 0; sample_idx2 != sample_idx1; ++sample_idx2) {
                  const double kinship_coeff = ComputeKinship(results_iter, singleton_het_cts[sample_idx2], singleton_hom_cts[sample_idx2], singleton_het1_ct, singleton_hom1_ct);
                  if (kinship_table && (kinship_coeff > king_cutoff)) {
                    SetBit(sample_idx2, &(kinship_table[sample_idx1 * sample_ctl]));
                    SetBit(sample_idx1, &(kinship_table[sample_idx2 * sample_ctl]));
                  }
                  write_row[sample_idx2] = S_CAST(float, kinship_coeff);
                  results_iter = &(results_iter[homhom_needed_p4]);
                }
                if (is_squarex) {
                  write_row[sample_idx1] = 0.5f;
                  if (is_square0) {
                    const uint32_t right_fill_idx = sample_idx1 + 1;
                    ZeroFArr(sample_ct - right_fill_idx, &(write_row[right_fill_idx]));
                  } else {
                    const uint32_t* results_iter2 = &(results_iter[sample_idx1 * homhom_needed_p4]);
                    for (uint32_t sample_idx2 = sample_idx1 + 1; sample_idx2 != sample_ct; ++sample_idx2) {
                      write_row[sample_idx2] = S_CAST(float, ComputeKinship(results_iter2, singleton_het1_ct, singleton_hom1_ct, singleton_het_cts[sample_idx2], singleton_hom_cts[sample_idx2]));
                      results_iter2 = &(results_iter2[sample_idx2 * homhom_needed_p4]);
                    }
                  }
                } else {
                  row_byte_ct = sample_idx1 * sizeof(float);
                }
                if (unlikely(fwrite_checked(write_row, row_byte_ct, outfile))) {
                  goto CalcKing_ret_WRITE_FAIL;
                }
              }
            } else {
              double* write_row = R_CAST(double*, numbuf);
              uintptr_t row_byte_ct = sample_ct * sizeof(double);
              for (; sample_idx1 != row_end_idx; ++sample_idx1) {
                const uint32_t singleton_het1_ct = singleton_het_cts[sample_idx1];
                const uint32_t singleton_hom1_ct = singleton_hom_cts[sample_idx1];
                for (uint32_t sample_idx2 = 0; sample_idx2 != sample_idx1; ++sample_idx2) {
                  const double kinship_coeff = ComputeKinship(results_iter, singleton_het_cts[sample_idx2], singleton_hom_cts[sample_idx2], singleton_het1_ct, singleton_hom1_ct);
                  if (kinship_table && (kinship_coeff > king_cutoff)) {
                    SetBit(sample_idx2, &(kinship_table[sample_idx1 * sample_ctl]));
                    SetBit(sample_idx1, &(kinship_table[sample_idx2 * sample_ctl]));
                  }
                  write_row[sample_idx2] = kinship_coeff;
                  results_iter = &(results_iter[homhom_needed_p4]);
                }
                if (is_squarex) {
                  write_row[sample_idx1] = 0.5;
                  if (is_square0) {
                    const uint32_t right_fill_idx = sample_idx1 + 1;
                    ZeroDArr(sample_ct - right_fill_idx, &(write_row[right_fill_idx]));
                  } else {
                    const uint32_t* results_iter2 = &(results_iter[sample_idx1 * homhom_needed_p4]);
                    for (uint32_t sample_idx2 = sample_idx1 + 1; sample_idx2 != sample_ct; ++sample_idx2) {
                      write_row[sample_idx2] = ComputeKinship(results_iter2, singleton_het1_ct, singleton_hom1_ct, singleton_het_cts[sample_idx2], singleton_hom_cts[sample_idx2]);
                      results_iter2 = &(results_iter2[sample_idx2 * homhom_needed_p4]);
                    }
                  }
                } else {
                  row_byte_ct = sample_idx1 * sizeof(double);
                }
                if (unlikely(fwrite_checked(write_row, row_byte_ct, outfile))) {
                  goto CalcKing_ret_WRITE_FAIL;
                }
              }
            }
          }
        }
        if (king_flags & kfKingColAll) {
          uintptr_t* kinship_table_backup = nullptr;
          if (matrix_shape) {
            // We already updated the table; don't do it again.
            kinship_table_backup = kinship_table;
            kinship_table = nullptr;
          }
          const uint32_t king_col_id = king_flags & kfKingColId;
          const uint32_t king_col_nsnp = king_flags & kfKingColNsnp;
          const uint32_t king_col_hethet = king_flags & kfKingColHethet;
          const uint32_t king_col_ibs0 = king_flags & kfKingColIbs0;
          const uint32_t king_col_ibs1 = king_flags & kfKingColIbs1;
          const uint32_t king_col_hamming = king_flags & kfKingColHamming;
          const uint32_t king_col_kinship = king_flags & kfKingColKinship;
          const uint32_t report_counts = king_flags & kfKingCounts;
          uint32_t* results_iter = dense_ctx.king_counts;
          double nonmiss_recip = 0.0;
          for (uint32_t sample_idx1 = row_start_idx; sample_idx1 != row_end_idx; ++sample_idx1) {
            const char* sample_fmtid1 = &(collapsed_sample_fmtids[max_sample_fmtid_blen * sample_idx1]);
            const uint32_t singleton_het1_ct = singleton_het_cts[sample_idx1];
            const uint32_t singleton_hom1_ct = singleton_hom_cts[sample_idx1];
            const uint32_t sample_fmtid1_slen = strlen(sample_fmtid1);
            for (uint32_t sample_idx2 = 0; sample_idx2 != sample_idx1; ++sample_idx2, results_iter = &(results_iter[homhom_needed_p4])) {
              const uint32_t singleton_het2_ct = singleton_het_cts[sample_idx2];
              const uint32_t singleton_hom2_ct = singleton_hom_cts[sample_idx2];
              const uint32_t ibs0_ct = results_iter[kKingOffsetIbs0] + singleton_hom2_ct + singleton_hom1_ct;
              const uint32_t hethet_ct = results_iter[kKingOffsetHethet];
              // '2' here refers to the larger index, so this is swapped
              const uint32_t het2hom1_ct = results_iter[kKingOffsetHet2Hom1] + singleton_het1_ct;
              const uint32_t het1hom2_ct = results_iter[kKingOffsetHet1Hom2] + singleton_het2_ct;
              const intptr_t smaller_het_ct = hethet_ct + MINV(het1hom2_ct, het2hom1_ct);
              const double kinship_coeff = 0.5 - (S_CAST(double, 4 * S_CAST(intptr_t, ibs0_ct) + het1hom2_ct + het2hom1_ct) / S_CAST(double, 4 * smaller_het_ct));
              if (kinship_table && (kinship_coeff > king_cutoff)) {
                SetBit(sample_idx2, &(kinship_table[sample_idx1 * sample_ctl]));
                SetBit(sample_idx1, &(kinship_table[sample_idx2 * sample_ctl]));
              }
              // edge case fix (18 Nov 2017): kinship_coeff can be -inf when
              // smaller_het_ct is zero.  Don't filter those lines out when
              // --king-table-filter wasn't specified.
              if ((king_table_filter != -DBL_MAX) && (kinship_coeff < king_table_filter)) {
                ++king_table_filter_ct;
                continue;
              }
              if (king_col_id) {
                cswritetp = memcpyax(cswritetp, sample_fmtid1, sample_fmtid1_slen, '\t');
                cswritetp = strcpyax(cswritetp, &(collapsed_sample_fmtids[max_sample_fmtid_blen * sample_idx2]), '\t');
              }
              if (homhom_needed_p4 == 5) {
                const uint32_t homhom_ct = results_iter[kKingOffsetHomhom] + sparse_variant_ct - singleton_het2_ct - singleton_missing_cts[sample_idx2] - singleton_het1_ct - singleton_missing_cts[sample_idx1];
                const uint32_t nonmiss_ct = het1hom2_ct + het2hom1_ct + homhom_ct + hethet_ct;
                if (king_col_nsnp) {
                  cswritetp = u32toa_x(nonmiss_ct, '\t', cswritetp);
                }
                if (!report_counts) {
                  nonmiss_recip = 1.0 / u31tod(nonmiss_ct);
                }
              }
              if (king_col_hethet) {
                if (report_counts) {
                  cswritetp = u32toa(hethet_ct, cswritetp);
                } else {
                  cswritetp = dtoa_g(nonmiss_recip * u31tod(hethet_ct), cswritetp);
                }
                *cswritetp++ = '\t';
              }
              if (king_col_ibs0) {
                if (report_counts) {
                  cswritetp = u32toa(ibs0_ct, cswritetp);
                } else {
                  cswritetp = dtoa_g(nonmiss_recip * u31tod(ibs0_ct), cswritetp);
                }
                *cswritetp++ = '\t';
              }
              if (king_col_ibs1) {
                if (report_counts) {
                  cswritetp = u32toa_x(het1hom2_ct, '\t', cswritetp);
                  cswritetp = u32toa(het2hom1_ct, cswritetp);
                } else {
                  cswritetp = dtoa_g(nonmiss_recip * u31tod(het1hom2_ct), cswritetp);
                  *cswritetp++ = '\t';
                  cswritetp = dtoa_g(nonmiss_recip * u31tod(het2hom1_ct), cswritetp);
                }
                *cswritetp++ = '\t';
              }
              if (king_col_hamming) {
                const uint32_t hamming_dist = 2 * ibs0_ct + het1hom2_ct + het2hom1_ct;
                if (report_counts) {
                  cswritetp = u32toa(hamming_dist, cswritetp);
                } else {
                  cswritetp = dtoa_g(nonmiss_recip * 0.5 * S_CAST(double, hamming_dist), cswritetp);
                }
                *cswritetp++ = '\t';
              }
              if (king_col_kinship) {
                cswritetp = dtoa_g(kinship_coeff, cswritetp);
                ++cswritetp;
              }
              DecrAppendBinaryEoln(&cswritetp);
              if (unlikely(Cswrite(&csst, &cswritetp))) {
                goto CalcKing_ret_WRITE_FAIL;
              }
            }
          }

          if (matrix_shape) {
            kinship_table = kinship_table_backup;
          }
        }
      } else {
        printf("\r%s pass %u/%u: Condensing...                \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", flagname, pass_idx_p1, pass_ct);
        fflush(stdout);
        uint32_t* results_iter = dense_ctx.king_counts;
        for (uint32_t sample_idx1 = row_start_idx; sample_idx1 != row_end_idx; ++sample_idx1) {
          const uint32_t singleton_het1_ct = singleton_het_cts[sample_idx1];
          const uint32_t singleton_hom1_ct = singleton_hom_cts[sample_idx1];
          for (uint32_t sample_idx2 = 0; sample_idx2 != sample_idx1; ++sample_idx2) {
            const double kinship_coeff = ComputeKinship(results_iter, singleton_het_cts[sample_idx2], singleton_hom_cts[sample_idx2], singleton_het1_ct, singleton_hom1_ct);
            if (kinship_coeff > king_cutoff) {
              SetBit(sample_idx2, &(kinship_table[sample_idx1 * sample_ctl]));
              SetBit(sample_idx1, &(kinship_table[sample_idx2 * sample_ctl]));
            }
            results_iter = &(results_iter[homhom_needed_p4]);
          }
        }
      }
      fputs(" done.\n", stdout);
    }
    logprintf("%s: %u variant%s processed.\n", flagname, variant_ct, (variant_ct == 1)? "" : "s");
    // end-of-loop operations
    if (matrix_shape) {
      if (!(king_flags & (kfKingMatrixBin | kfKingMatrixBin4))) {
        if (unlikely(CswriteCloseNull(&css, cswritep))) {
          goto CalcKing_ret_WRITE_FAIL;
        }
      } else {
        if (unlikely(fclose_null(&outfile))) {
          goto CalcKing_ret_WRITE_FAIL;
        }
      }
      // Necessary to regenerate filename since it may have been overwritten by
      // --make-king-table.
      SetKingMatrixFname(king_flags, parallel_idx, parallel_tot, outname_end);

      char* write_iter = strcpya_k(g_logbuf, "Results written to ");
      const uint32_t outname_base_slen = S_CAST(uintptr_t, outname_end - outname);
      write_iter = memcpya(write_iter, outname, outname_base_slen + strlen(outname_end));
      write_iter = strcpya_k(write_iter, " and ");
      strcpy_k(&(outname_end[5]), ".id");
      write_iter = memcpya(write_iter, outname, outname_base_slen + 8);
      strcpy_k(write_iter, " .\n");
      WordWrapB(0);
      logputsb();
      reterr = WriteSampleIds(sample_include, siip, outname, sample_ct);
      if (unlikely(reterr)) {
        goto CalcKing_ret_1;
      }
    }
    if (king_flags & kfKingColAll) {
      if (unlikely(CswriteCloseNull(&csst, cswritetp))) {
        goto CalcKing_ret_WRITE_FAIL;
      }
      SetKingTableFname(king_flags, parallel_idx, parallel_tot, outname_end);
      char* write_iter = strcpya_k(g_logbuf, "Results written to ");
      const uint32_t outname_base_slen = S_CAST(uintptr_t, outname_end - outname);
      write_iter = memcpya(write_iter, outname, outname_base_slen + strlen(outname_end));
      if ((!parallel_idx) && (!(king_flags & kfKingColId))) {
        write_iter = strcpya_k(write_iter, " and ");
        strcpy_k(&(outname_end[5]), ".id");
        write_iter = memcpya(write_iter, outname, outname_base_slen + 8);
        strcpy_k(write_iter, " .\n");
        WordWrapB(0);
        logputsb();
        reterr = WriteSampleIds(sample_include, siip, outname, sample_ct);
        if (unlikely(reterr)) {
          goto CalcKing_ret_1;
        }
      } else {
        strcpy_k(write_iter, " .\n");
        WordWrapB(0);
        logputsb();
      }
      if (king_table_filter != -DBL_MAX) {
        const uint64_t grand_tot_cells = (S_CAST(uint64_t, grand_row_end_idx) * (grand_row_end_idx - 1) - S_CAST(uint64_t, grand_row_start_idx) * (grand_row_start_idx - 1)) / 2;
        const uint64_t reported_ct = grand_tot_cells - king_table_filter_ct;
        logprintf("--king-table-filter: %" PRIu64 " relationship%s reported (%" PRIu64 " filtered out).\n", reported_ct, (reported_ct == 1)? "" : "s", king_table_filter_ct);
      }
    }
    if (kinship_table) {
      BigstackReset(sample_include_cumulative_popcounts);
      *sample_ct_ptr = sample_ct;
      if (unlikely(KinshipPruneDestructive(kinship_table, sample_include, sample_ct_ptr))) {
        goto CalcKing_ret_NOMEM;
      }
    }
  }
  while (0) {
  CalcKing_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CalcKing_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  CalcKing_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  CalcKing_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  CalcKing_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  CalcKing_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  CalcKing_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 CalcKing_ret_1:
  CleanupThreads(&tg);
  CswriteCloseCond(&csst, cswritetp);
  CswriteCloseCond(&css, cswritep);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

#ifdef USE_SSE42
void IncrKingSubset(const uint32_t* loaded_sample_idx_pairs, const uintptr_t* smaj_hom, const uintptr_t* smaj_ref2het, uint32_t start_idx, uint32_t end_idx, uint32_t* king_counts) {
  const uint32_t* sample_idx_pair_iter = &(loaded_sample_idx_pairs[(2 * k1LU) * start_idx]);
  const uint32_t* sample_idx_pair_stop = &(loaded_sample_idx_pairs[(2 * k1LU) * end_idx]);
  uint32_t* king_counts_iter = &(king_counts[(4 * k1LU) * start_idx]);
  while (sample_idx_pair_iter != sample_idx_pair_stop) {
    // technically overflows for huge sample_ct
    const uint32_t first_offset = (*sample_idx_pair_iter++) * kKingMultiplexWords;
    const uint32_t second_offset = (*sample_idx_pair_iter++) * kKingMultiplexWords;
    const uintptr_t* first_hom = &(smaj_hom[first_offset]);
    const uintptr_t* first_ref2het = &(smaj_ref2het[first_offset]);
    const uintptr_t* second_hom = &(smaj_hom[second_offset]);
    const uintptr_t* second_ref2het = &(smaj_ref2het[second_offset]);
    uint32_t acc_ibs0 = 0;
    uint32_t acc_hethet = 0;
    uint32_t acc_het2hom1 = 0;
    uint32_t acc_het1hom2 = 0;
    for (uint32_t widx = 0; widx != kKingMultiplexWords; ++widx) {
      const uintptr_t hom1 = first_hom[widx];
      const uintptr_t hom2 = second_hom[widx];
      const uintptr_t ref2het1 = first_ref2het[widx];
      const uintptr_t ref2het2 = second_ref2het[widx];
      const uintptr_t homhom = hom1 & hom2;
      const uintptr_t het1 = ref2het1 & (~hom1);
      const uintptr_t het2 = ref2het2 & (~hom2);
      acc_ibs0 += PopcountWord((ref2het1 ^ ref2het2) & homhom);
      acc_hethet += PopcountWord(het1 & het2);
      acc_het2hom1 += PopcountWord(hom1 & het2);
      acc_het1hom2 += PopcountWord(hom2 & het1);
    }
    *king_counts_iter++ += acc_ibs0;
    *king_counts_iter++ += acc_hethet;
    *king_counts_iter++ += acc_het2hom1;
    *king_counts_iter++ += acc_het1hom2;
  }
}

void IncrKingSubsetHomhom(const uint32_t* loaded_sample_idx_pairs, const uintptr_t* smaj_hom, const uintptr_t* smaj_ref2het, uint32_t start_idx, uint32_t end_idx, uint32_t* king_counts) {
  const uint32_t* sample_idx_pair_iter = &(loaded_sample_idx_pairs[(2 * k1LU) * start_idx]);
  const uint32_t* sample_idx_pair_stop = &(loaded_sample_idx_pairs[(2 * k1LU) * end_idx]);
  uint32_t* king_counts_iter = &(king_counts[(5 * k1LU) * start_idx]);
  while (sample_idx_pair_iter != sample_idx_pair_stop) {
    // technically overflows for huge sample_ct
    const uint32_t first_offset = (*sample_idx_pair_iter++) * kKingMultiplexWords;
    const uint32_t second_offset = (*sample_idx_pair_iter++) * kKingMultiplexWords;
    const uintptr_t* first_hom = &(smaj_hom[first_offset]);
    const uintptr_t* first_ref2het = &(smaj_ref2het[first_offset]);
    const uintptr_t* second_hom = &(smaj_hom[second_offset]);
    const uintptr_t* second_ref2het = &(smaj_ref2het[second_offset]);
    uint32_t acc_homhom = 0;
    uint32_t acc_ibs0 = 0;
    uint32_t acc_hethet = 0;
    uint32_t acc_het2hom1 = 0;
    uint32_t acc_het1hom2 = 0;
    for (uint32_t widx = 0; widx != kKingMultiplexWords; ++widx) {
      const uintptr_t hom1 = first_hom[widx];
      const uintptr_t hom2 = second_hom[widx];
      const uintptr_t ref2het1 = first_ref2het[widx];
      const uintptr_t ref2het2 = second_ref2het[widx];
      const uintptr_t homhom = hom1 & hom2;
      const uintptr_t het1 = ref2het1 & (~hom1);
      const uintptr_t het2 = ref2het2 & (~hom2);
      acc_homhom += PopcountWord(homhom);
      acc_ibs0 += PopcountWord((ref2het1 ^ ref2het2) & homhom);
      acc_hethet += PopcountWord(het1 & het2);
      acc_het2hom1 += PopcountWord(hom1 & het2);
      acc_het1hom2 += PopcountWord(hom2 & het1);
    }
    *king_counts_iter++ += acc_ibs0;
    *king_counts_iter++ += acc_hethet;
    *king_counts_iter++ += acc_het2hom1;
    *king_counts_iter++ += acc_het1hom2;
    *king_counts_iter++ += acc_homhom;
  }
}
#else
void IncrKingSubset(const uint32_t* loaded_sample_idx_pairs, const uintptr_t* smaj_hom, const uintptr_t* smaj_ref2het, uint32_t start_idx, uint32_t end_idx, uint32_t* king_counts) {
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const uint32_t* sample_idx_pair_iter = &(loaded_sample_idx_pairs[(2 * k1LU) * start_idx]);
  const uint32_t* sample_idx_pair_stop = &(loaded_sample_idx_pairs[(2 * k1LU) * end_idx]);
  uint32_t* king_counts_iter = &(king_counts[(4 * k1LU) * start_idx]);
  while (sample_idx_pair_iter != sample_idx_pair_stop) {
    // technically overflows for huge sample_ct
    const uint32_t first_offset = (*sample_idx_pair_iter++) * kKingMultiplexWords;
    const uint32_t second_offset = (*sample_idx_pair_iter++) * kKingMultiplexWords;
    const VecW* first_hom = R_CAST(const VecW*, &(smaj_hom[first_offset]));
    const VecW* first_ref2het = R_CAST(const VecW*, &(smaj_ref2het[first_offset]));
    const VecW* second_hom = R_CAST(const VecW*, &(smaj_hom[second_offset]));
    const VecW* second_ref2het = R_CAST(const VecW*, &(smaj_ref2het[second_offset]));
    UniVec acc_ibs0;
    UniVec acc_hethet;
    UniVec acc_het2hom1;
    UniVec acc_het1hom2;
    acc_ibs0.vw = vecw_setzero();
    acc_hethet.vw = vecw_setzero();
    acc_het2hom1.vw = vecw_setzero();
    acc_het1hom2.vw = vecw_setzero();
    for (uint32_t vec_idx = 0; vec_idx < kKingMultiplexVecs; vec_idx += 3) {
      VecW hom1 = first_hom[vec_idx];
      VecW hom2 = second_hom[vec_idx];
      VecW ref2het1 = first_ref2het[vec_idx];
      VecW ref2het2 = second_ref2het[vec_idx];
      VecW het1 = vecw_and_notfirst(hom1, ref2het1);
      VecW het2 = vecw_and_notfirst(hom2, ref2het2);
      VecW agg_ibs0 = (ref2het1 ^ ref2het2) & (hom1 & hom2);
      VecW agg_hethet = het1 & het2;
      VecW agg_het2hom1 = hom1 & het2;
      VecW agg_het1hom2 = hom2 & het1;
      agg_ibs0 = agg_ibs0 - (vecw_srli(agg_ibs0, 1) & m1);
      agg_hethet = agg_hethet - (vecw_srli(agg_hethet, 1) & m1);
      agg_het2hom1 = agg_het2hom1 - (vecw_srli(agg_het2hom1, 1) & m1);
      agg_het1hom2 = agg_het1hom2 - (vecw_srli(agg_het1hom2, 1) & m1);
      agg_ibs0 = (agg_ibs0 & m2) + (vecw_srli(agg_ibs0, 2) & m2);
      agg_hethet = (agg_hethet & m2) + (vecw_srli(agg_hethet, 2) & m2);
      agg_het2hom1 = (agg_het2hom1 & m2) + (vecw_srli(agg_het2hom1, 2) & m2);
      agg_het1hom2 = (agg_het1hom2 & m2) + (vecw_srli(agg_het1hom2, 2) & m2);

      for (uint32_t offset = 1; offset != 3; ++offset) {
        hom1 = first_hom[vec_idx + offset];
        hom2 = second_hom[vec_idx + offset];
        ref2het1 = first_ref2het[vec_idx + offset];
        ref2het2 = second_ref2het[vec_idx + offset];
        het1 = vecw_and_notfirst(hom1, ref2het1);
        het2 = vecw_and_notfirst(hom2, ref2het2);
        VecW cur_ibs0 = (ref2het1 ^ ref2het2) & (hom1 & hom2);
        VecW cur_hethet = het1 & het2;
        VecW cur_het2hom1 = hom1 & het2;
        VecW cur_het1hom2 = hom2 & het1;
        cur_ibs0 = cur_ibs0 - (vecw_srli(cur_ibs0, 1) & m1);
        cur_hethet = cur_hethet - (vecw_srli(cur_hethet, 1) & m1);
        cur_het2hom1 = cur_het2hom1 - (vecw_srli(cur_het2hom1, 1) & m1);
        cur_het1hom2 = cur_het1hom2 - (vecw_srli(cur_het1hom2, 1) & m1);
        agg_ibs0 += (cur_ibs0 & m2) + (vecw_srli(cur_ibs0, 2) & m2);
        agg_hethet += (cur_hethet & m2) + (vecw_srli(cur_hethet, 2) & m2);
        agg_het2hom1 += (cur_het2hom1 & m2) + (vecw_srli(cur_het2hom1, 2) & m2);
        agg_het1hom2 += (cur_het1hom2 & m2) + (vecw_srli(cur_het1hom2, 2) & m2);
      }
      acc_ibs0.vw = acc_ibs0.vw + (agg_ibs0 & m4) + (vecw_srli(agg_ibs0, 4) & m4);
      acc_hethet.vw = acc_hethet.vw + (agg_hethet & m4) + (vecw_srli(agg_hethet, 4) & m4);
      acc_het2hom1.vw = acc_het2hom1.vw + (agg_het2hom1 & m4) + (vecw_srli(agg_het2hom1, 4) & m4);
      acc_het1hom2.vw = acc_het1hom2.vw + (agg_het1hom2 & m4) + (vecw_srli(agg_het1hom2, 4) & m4);
    }
    const VecW m8 = VCONST_W(kMask00FF);
    acc_ibs0.vw = (acc_ibs0.vw & m8) + (vecw_srli(acc_ibs0.vw, 8) & m8);
    acc_hethet.vw = (acc_hethet.vw & m8) + (vecw_srli(acc_hethet.vw, 8) & m8);
    acc_het2hom1.vw = (acc_het2hom1.vw & m8) + (vecw_srli(acc_het2hom1.vw, 8) & m8);
    acc_het1hom2.vw = (acc_het1hom2.vw & m8) + (vecw_srli(acc_het1hom2.vw, 8) & m8);
    *king_counts_iter++ += UniVecHsum16(acc_ibs0);
    *king_counts_iter++ += UniVecHsum16(acc_hethet);
    *king_counts_iter++ += UniVecHsum16(acc_het2hom1);
    *king_counts_iter++ += UniVecHsum16(acc_het1hom2);
  }
}

void IncrKingSubsetHomhom(const uint32_t* loaded_sample_idx_pairs, const uintptr_t* smaj_hom, const uintptr_t* smaj_ref2het, uint32_t start_idx, uint32_t end_idx, uint32_t* king_counts) {
  const VecW m1 = VCONST_W(kMask5555);
  const VecW m2 = VCONST_W(kMask3333);
  const VecW m4 = VCONST_W(kMask0F0F);
  const uint32_t* sample_idx_pair_iter = &(loaded_sample_idx_pairs[(2 * k1LU) * start_idx]);
  const uint32_t* sample_idx_pair_stop = &(loaded_sample_idx_pairs[(2 * k1LU) * end_idx]);
  uint32_t* king_counts_iter = &(king_counts[(5 * k1LU) * start_idx]);
  while (sample_idx_pair_iter != sample_idx_pair_stop) {
    // technically overflows for huge sample_ct
    const uint32_t first_offset = (*sample_idx_pair_iter++) * kKingMultiplexWords;
    const uint32_t second_offset = (*sample_idx_pair_iter++) * kKingMultiplexWords;
    const VecW* first_hom = R_CAST(const VecW*, &(smaj_hom[first_offset]));
    const VecW* first_ref2het = R_CAST(const VecW*, &(smaj_ref2het[first_offset]));
    const VecW* second_hom = R_CAST(const VecW*, &(smaj_hom[second_offset]));
    const VecW* second_ref2het = R_CAST(const VecW*, &(smaj_ref2het[second_offset]));
    UniVec acc_homhom;
    UniVec acc_ibs0;
    UniVec acc_hethet;
    UniVec acc_het2hom1;
    UniVec acc_het1hom2;
    acc_homhom.vw = vecw_setzero();
    acc_ibs0.vw = vecw_setzero();
    acc_hethet.vw = vecw_setzero();
    acc_het2hom1.vw = vecw_setzero();
    acc_het1hom2.vw = vecw_setzero();
    for (uint32_t vec_idx = 0; vec_idx < kKingMultiplexVecs; vec_idx += 3) {
      VecW hom1 = first_hom[vec_idx];
      VecW hom2 = second_hom[vec_idx];
      VecW ref2het1 = first_ref2het[vec_idx];
      VecW ref2het2 = second_ref2het[vec_idx];
      VecW agg_homhom = hom1 & hom2;
      VecW het1 = vecw_and_notfirst(hom1, ref2het1);
      VecW het2 = vecw_and_notfirst(hom2, ref2het2);
      VecW agg_ibs0 = (ref2het1 ^ ref2het2) & agg_homhom;
      VecW agg_hethet = het1 & het2;
      VecW agg_het2hom1 = hom1 & het2;
      VecW agg_het1hom2 = hom2 & het1;
      agg_homhom = agg_homhom - (vecw_srli(agg_homhom, 1) & m1);
      agg_ibs0 = agg_ibs0 - (vecw_srli(agg_ibs0, 1) & m1);
      agg_hethet = agg_hethet - (vecw_srli(agg_hethet, 1) & m1);
      agg_het2hom1 = agg_het2hom1 - (vecw_srli(agg_het2hom1, 1) & m1);
      agg_het1hom2 = agg_het1hom2 - (vecw_srli(agg_het1hom2, 1) & m1);
      agg_homhom = (agg_homhom & m2) + (vecw_srli(agg_homhom, 2) & m2);
      agg_ibs0 = (agg_ibs0 & m2) + (vecw_srli(agg_ibs0, 2) & m2);
      agg_hethet = (agg_hethet & m2) + (vecw_srli(agg_hethet, 2) & m2);
      agg_het2hom1 = (agg_het2hom1 & m2) + (vecw_srli(agg_het2hom1, 2) & m2);
      agg_het1hom2 = (agg_het1hom2 & m2) + (vecw_srli(agg_het1hom2, 2) & m2);

      for (uint32_t offset = 1; offset != 3; ++offset) {
        hom1 = first_hom[vec_idx + offset];
        hom2 = second_hom[vec_idx + offset];
        ref2het1 = first_ref2het[vec_idx + offset];
        ref2het2 = second_ref2het[vec_idx + offset];
        VecW cur_homhom = hom1 & hom2;
        het1 = vecw_and_notfirst(hom1, ref2het1);
        het2 = vecw_and_notfirst(hom2, ref2het2);
        VecW cur_ibs0 = (ref2het1 ^ ref2het2) & cur_homhom;
        VecW cur_hethet = het1 & het2;
        VecW cur_het2hom1 = hom1 & het2;
        VecW cur_het1hom2 = hom2 & het1;
        cur_homhom = cur_homhom - (vecw_srli(cur_homhom, 1) & m1);
        cur_ibs0 = cur_ibs0 - (vecw_srli(cur_ibs0, 1) & m1);
        cur_hethet = cur_hethet - (vecw_srli(cur_hethet, 1) & m1);
        cur_het2hom1 = cur_het2hom1 - (vecw_srli(cur_het2hom1, 1) & m1);
        cur_het1hom2 = cur_het1hom2 - (vecw_srli(cur_het1hom2, 1) & m1);
        agg_homhom += (cur_homhom & m2) + (vecw_srli(cur_homhom, 2) & m2);
        agg_ibs0 += (cur_ibs0 & m2) + (vecw_srli(cur_ibs0, 2) & m2);
        agg_hethet += (cur_hethet & m2) + (vecw_srli(cur_hethet, 2) & m2);
        agg_het2hom1 += (cur_het2hom1 & m2) + (vecw_srli(cur_het2hom1, 2) & m2);
        agg_het1hom2 += (cur_het1hom2 & m2) + (vecw_srli(cur_het1hom2, 2) & m2);
      }
      acc_homhom.vw = acc_homhom.vw + (agg_homhom & m4) + (vecw_srli(agg_homhom, 4) & m4);
      acc_ibs0.vw = acc_ibs0.vw + (agg_ibs0 & m4) + (vecw_srli(agg_ibs0, 4) & m4);
      acc_hethet.vw = acc_hethet.vw + (agg_hethet & m4) + (vecw_srli(agg_hethet, 4) & m4);
      acc_het2hom1.vw = acc_het2hom1.vw + (agg_het2hom1 & m4) + (vecw_srli(agg_het2hom1, 4) & m4);
      acc_het1hom2.vw = acc_het1hom2.vw + (agg_het1hom2 & m4) + (vecw_srli(agg_het1hom2, 4) & m4);
    }
    const VecW m8 = VCONST_W(kMask00FF);
    acc_homhom.vw = (acc_homhom.vw & m8) + (vecw_srli(acc_homhom.vw, 8) & m8);
    acc_ibs0.vw = (acc_ibs0.vw & m8) + (vecw_srli(acc_ibs0.vw, 8) & m8);
    acc_hethet.vw = (acc_hethet.vw & m8) + (vecw_srli(acc_hethet.vw, 8) & m8);
    acc_het2hom1.vw = (acc_het2hom1.vw & m8) + (vecw_srli(acc_het2hom1.vw, 8) & m8);
    acc_het1hom2.vw = (acc_het1hom2.vw & m8) + (vecw_srli(acc_het1hom2.vw, 8) & m8);
    *king_counts_iter++ += UniVecHsum16(acc_ibs0);
    *king_counts_iter++ += UniVecHsum16(acc_hethet);
    *king_counts_iter++ += UniVecHsum16(acc_het2hom1);
    *king_counts_iter++ += UniVecHsum16(acc_het1hom2);
    *king_counts_iter++ += UniVecHsum16(acc_homhom);
  }
}
#endif

typedef struct CalcKingTableSubsetCtxStruct {
  uintptr_t* smaj_hom[2];
  uintptr_t* smaj_ref2het[2];
  uint32_t* loaded_sample_idx_pairs;
  uint32_t homhom_needed;

  uint32_t* thread_start;

  uint32_t* king_counts;
} CalcKingTableSubsetCtx;

THREAD_FUNC_DECL CalcKingTableSubsetThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcKingTableSubsetCtx* ctx = S_CAST(CalcKingTableSubsetCtx*, arg->sharedp->context);

  const uint32_t start_idx = ctx->thread_start[tidx];
  const uint32_t end_idx = ctx->thread_start[tidx + 1];
  const uint32_t homhom_needed = ctx->homhom_needed;
  uint32_t parity = 0;
  do {
    if (homhom_needed) {
      IncrKingSubsetHomhom(ctx->loaded_sample_idx_pairs, ctx->smaj_hom[parity], ctx->smaj_ref2het[parity], start_idx, end_idx, ctx->king_counts);
    } else {
      IncrKingSubset(ctx->loaded_sample_idx_pairs, ctx->smaj_hom[parity], ctx->smaj_ref2het[parity], start_idx, end_idx, ctx->king_counts);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr KingTableSubsetLoad(const char* sorted_xidbox, const uint32_t* xid_map, const uintptr_t* sample_require, uintptr_t max_xid_blen, uintptr_t orig_sample_ct, double king_table_subset_thresh, XidMode xid_mode, RelConcordanceCheckMode rel_or_concordance_check, uint32_t require_xor, uint32_t kinship_skip, uint32_t is_first_parallel_scan, uint64_t pair_idx_start, uint64_t pair_idx_stop, uintptr_t line_idx, TextStream* txsp, uint64_t* pair_idx_ptr, uint32_t* loaded_sample_idx_pairs, char* idbuf) {
  PglErr reterr = kPglRetSuccess;
  {
    // In --king-table-require case, sample1_in_require + sample2_in_require
    // can be either 1 or 2.  In --king-table-require-xor case, it must be 1.
    const uint32_t require_sum_mask = require_xor? 1 : 3;
    uint64_t pair_idx = *pair_idx_ptr;
    // Assumes header line already read if pair_idx == 0, and if pair_idx is
    // positive, we're that far into the file.
    uint32_t* loaded_sample_idx_pairs_iter = loaded_sample_idx_pairs;
    ++line_idx;
    for (char* line_iter = TextLineEnd(txsp); TextGetUnsafe2(txsp, &line_iter); line_iter = AdvPastDelim(line_iter, '\n'), ++line_idx) {
      const char* linebuf_iter = line_iter;
      uint32_t sample_uidx1;
      if (SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &linebuf_iter, &sample_uidx1, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto KingTableSubsetLoad_ret_MISSING_TOKENS;
        }
        line_iter = K_CAST(char*, linebuf_iter);
        continue;
      }
      linebuf_iter = FirstNonTspace(linebuf_iter);
      if (rel_or_concordance_check) {
        // linebuf_iter must point to the start of the second FID, while
        // line_iter points to the start of the first.
        const uint32_t first_fid_slen = CurTokenEnd(line_iter) - line_iter;
        const uint32_t second_fid_slen = CurTokenEnd(linebuf_iter) - linebuf_iter;
        if ((first_fid_slen != second_fid_slen) || (!memequal(line_iter, linebuf_iter, first_fid_slen))) {
          line_iter = K_CAST(char*, linebuf_iter);
          continue;
        }
        if (rel_or_concordance_check == kRcCheckConcordance) {
          const char* first_iid = FirstNonTspace(&(line_iter[first_fid_slen + 1]));
          const char* second_iid = FirstNonTspace(&(linebuf_iter[second_fid_slen]));
          const uint32_t first_iid_slen = CurTokenEnd(first_iid) - first_iid;
          const uint32_t second_iid_slen = FirstSpaceOrEoln(second_iid) - second_iid;
          if ((first_iid_slen != second_iid_slen) || (!memequal(first_iid, second_iid, first_iid_slen))) {
            line_iter = K_CAST(char*, second_iid);
            continue;
          }
        }
      }
      uint32_t sample_uidx2;
      if (SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, 0, xid_mode, &linebuf_iter, &sample_uidx2, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto KingTableSubsetLoad_ret_MISSING_TOKENS;
        }
        line_iter = K_CAST(char*, linebuf_iter);
        continue;
      }
      if (unlikely(sample_uidx1 == sample_uidx2)) {
        // could technically be due to unloaded SID, so use inconsistent-input
        // error code
        snprintf(g_logbuf, kLogbufSize, "Error: Identical sample IDs on line %" PRIuPTR " of --king-table-subset file.\n", line_idx);
        goto KingTableSubsetLoad_ret_INCONSISTENT_INPUT_WW;
      }
      if (sample_require) {
        const uint32_t sample1_in_require = IsSet(sample_require, sample_uidx1);
        const uint32_t sample2_in_require = IsSet(sample_require, sample_uidx2);
        if (!((sample1_in_require + sample2_in_require) & require_sum_mask)) {
          line_iter = K_CAST(char*, linebuf_iter);
          continue;
        }
      }
      if (king_table_subset_thresh != -DBL_MAX) {
        linebuf_iter = FirstNonTspace(linebuf_iter);
        linebuf_iter = NextTokenMult0(linebuf_iter, kinship_skip);
        if (unlikely(!linebuf_iter)) {
          goto KingTableSubsetLoad_ret_MISSING_TOKENS;
        }
        double cur_kinship;
        const char* kinship_end = ScanadvDouble(linebuf_iter, &cur_kinship);
        if (!kinship_end) {
          line_iter = K_CAST(char*, linebuf_iter);
          continue;
        }
        if (unlikely(!IsSpaceOrEoln(*kinship_end))) {
          kinship_end = CurTokenEnd(kinship_end);
          *K_CAST(char*, kinship_end) = '\0';
          logerrprintfww("Error: Invalid numeric token '%s' on line %" PRIuPTR " of --king-table-subset file.\n", linebuf_iter, line_idx);
          goto KingTableSubsetLoad_ret_MALFORMED_INPUT;
        }
        if (cur_kinship < king_table_subset_thresh) {
          line_iter = K_CAST(char*, kinship_end);
          continue;
        }
      }
      line_iter = K_CAST(char*, linebuf_iter);
      if (pair_idx < pair_idx_start) {
        ++pair_idx;
        continue;
      }
      *loaded_sample_idx_pairs_iter++ = sample_uidx1;
      *loaded_sample_idx_pairs_iter++ = sample_uidx2;
      ++pair_idx;
      if (pair_idx == pair_idx_stop) {
        if (!is_first_parallel_scan) {
          TextSetPos(AdvPastDelim(line_iter, '\n'), txsp);
          goto KingTableSubsetLoad_finish;
        }
        // large --parallel job, first pass: count number of valid pairs, don't
        // save the remainder
        pair_idx_start = ~0LLU;
      }
    }
    if (unlikely(TextStreamErrcode2(txsp, &reterr))) {
      goto KingTableSubsetLoad_ret_TSTREAM_FAIL;
    }
  KingTableSubsetLoad_finish:
    *pair_idx_ptr = pair_idx;
  }
  while (0) {
  KingTableSubsetLoad_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--king-table-subset file", txsp);
    break;
  KingTableSubsetLoad_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  KingTableSubsetLoad_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --king-table-subset file has fewer tokens than expected.\n", line_idx);
  KingTableSubsetLoad_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
  return reterr;
}

typedef struct FidPairIteratorStruct {
  uint32_t block_start_idx;
  uint32_t block_end_idx;
  uint32_t idx1;
  uint32_t idx2;
} FidPairIterator;

void InitFidPairIterator(FidPairIterator* fpip) {
  fpip->block_start_idx = UINT32_MAX;  // deliberate overflow
  fpip->block_end_idx = 0;
  fpip->idx1 = 0;
  fpip->idx2 = 0;  // defensive
}

uint64_t CountRcCheckPairs(const char* nsorted_xidbox, RelConcordanceCheckMode rel_or_concordance_check, uintptr_t max_xid_blen, uintptr_t orig_sample_ct, char* idbuf) {
  uint64_t total = 0;
  for (uintptr_t block_start_idx = 0; block_start_idx != orig_sample_ct; ) {
    const char* fid_start = &(nsorted_xidbox[block_start_idx * max_xid_blen]);
    const char* fid_end = AdvToDelim(fid_start, '\t');
    uint32_t id_slen;
    if (rel_or_concordance_check == kRcCheckRel) {
      id_slen = fid_end - fid_start;
    } else {
      // concordance check
      const char* iid_end = AdvToDelim(&(fid_end[1]), '\t');
      id_slen = iid_end - fid_start;
    }
    memcpy(idbuf, fid_start, id_slen);
    idbuf[id_slen] = ' ';
    // bugfix (14 Jan 2020): forgot that natural-sorting was used...
    idbuf[id_slen + 1] = '\0';
    const uintptr_t block_end_idx = ExpsearchNsortStrLb(idbuf, nsorted_xidbox, max_xid_blen, orig_sample_ct, block_start_idx + 1);
    const uint64_t cur_block_size = block_end_idx - block_start_idx;
    total += (cur_block_size * (cur_block_size - 1)) / 2;
    block_start_idx = block_end_idx;
  }
  return total;
}

uint64_t CountRcCheckPairsSampleRequire(const char* nsorted_xidbox, const uint32_t* xid_map, const uintptr_t* sample_require, RelConcordanceCheckMode rel_or_concordance_check, uintptr_t max_xid_blen, uintptr_t orig_sample_ct, uint32_t require_xor, char* idbuf) {
  uint64_t total = 0;
  for (uintptr_t block_start_idx = 0; block_start_idx != orig_sample_ct; ) {
    const char* fid_start = &(nsorted_xidbox[block_start_idx * max_xid_blen]);
    const char* fid_end = AdvToDelim(fid_start, '\t');
    uint32_t id_slen;
    if (rel_or_concordance_check == kRcCheckRel) {
      id_slen = fid_end - fid_start;
    } else {
      const char* iid_end = AdvToDelim(&(fid_end[1]), '\t');
      id_slen = iid_end - fid_start;
    }
    memcpy(idbuf, fid_start, id_slen);
    idbuf[id_slen] = ' ';
    idbuf[id_slen + 1] = '\0';
    const uintptr_t block_end_idx = ExpsearchNsortStrLb(idbuf, nsorted_xidbox, max_xid_blen, orig_sample_ct, block_start_idx + 1);
    const uint64_t cur_block_size = block_end_idx - block_start_idx;
    uintptr_t required_ct = 0;
    for (uintptr_t block_idx = block_start_idx; block_idx != block_end_idx; ++block_idx) {
      const uint32_t sample_uidx = xid_map[block_idx];
      required_ct += IsSet(sample_require, sample_uidx);
    }
    const uint64_t optional_ct = cur_block_size - required_ct;
    if (!require_xor) {
      total += ((cur_block_size * (cur_block_size - 1)) - (optional_ct * (optional_ct - 1))) / 2;
    } else {
      total += required_ct * optional_ct;
    }
    block_start_idx = block_end_idx;
  }
  return total;
}

void GetRelCheckOrKTRequirePairs(const char* nsorted_xidbox, const uint32_t* xid_map, const uintptr_t* sample_require, uintptr_t max_xid_blen, uintptr_t orig_sample_ct, RelConcordanceCheckMode rel_or_concordance_check, uint32_t require_xor, uint32_t is_first_parallel_scan, uint64_t pair_idx_start, uint64_t pair_idx_stop, FidPairIterator* fpip, uint64_t* pair_idx_ptr, uint32_t* loaded_sample_idx_pairs, char* idbuf) {
  // Support "--make-king-table rel-check" and "--make-king-table
  // concordance-check" without an actual subset-file.  These loops aren't
  // really optimized, since actual rel-check use cases should have small
  // families where it doesn't matter.
  //
  // In the is_first_parallel_scan case, *pair_idx_ptr is set to the total
  // number of eligible pairs, even if it's greater than pair_idx_stop.
  // Otherwise, *pair_idx_ptr will not be set to larger than pair_idx_stop.
  uint32_t block_start_idx = fpip->block_start_idx;
  uint32_t block_end_idx = fpip->block_end_idx;
  uint32_t idx1 = fpip->idx1;
  uint32_t idx2 = fpip->idx2;
  uint64_t pair_idx = *pair_idx_ptr;
  uint32_t* loaded_sample_idx_pairs_iter = loaded_sample_idx_pairs;
  if (!sample_require) {
    while (1) {
      for (; idx1 != block_end_idx; ++idx1) {
        // idx1 >= idx2.
        uint32_t cur_pair_ct = idx1 - idx2;
        uint32_t idx2_stop = idx1;
        if (pair_idx_stop - pair_idx < cur_pair_ct) {
          cur_pair_ct = pair_idx_stop - pair_idx;
          idx2_stop = idx2 + cur_pair_ct;
        }
        if (pair_idx < pair_idx_start) {
          // possible in later --parallel shards
          const uint64_t skip_ct = pair_idx_start - pair_idx;
          if (skip_ct >= cur_pair_ct) {
            idx2 = idx2_stop;
          } else {
            idx2 += skip_ct;
          }
          // pair_idx is updated correctly after the inner loop
        }
        const uint32_t sample_uidx1 = xid_map[idx1];
        for (; idx2 != idx2_stop; ++idx2) {
          const uint32_t sample_uidx2 = xid_map[idx2];
          *loaded_sample_idx_pairs_iter++ = sample_uidx1;
          *loaded_sample_idx_pairs_iter++ = sample_uidx2;
        }
        pair_idx += cur_pair_ct;
        if (pair_idx == pair_idx_stop) {
          if (is_first_parallel_scan) {
            pair_idx = CountRcCheckPairs(nsorted_xidbox, rel_or_concordance_check, max_xid_blen, orig_sample_ct, idbuf);
          }
          goto GetRelCheckOrKTRequirePairs_early_exit;
        }
        idx2 = block_start_idx;
      }
      block_start_idx = block_end_idx;
      if (block_start_idx == orig_sample_ct) {
        break;
      }
      idx2 = block_start_idx;
      const char* fid_start = &(nsorted_xidbox[block_start_idx * max_xid_blen]);
      const char* fid_end = AdvToDelim(fid_start, '\t');
      uint32_t id_slen;
      if (rel_or_concordance_check == kRcCheckRel) {
        id_slen = fid_end - fid_start;
      } else {
        const char* iid_end = AdvToDelim(&(fid_end[1]), '\t');
        id_slen = iid_end - fid_start;
      }
      memcpy(idbuf, fid_start, id_slen);
      idbuf[id_slen] = ' ';
      idbuf[id_slen + 1] = '\0';
      block_end_idx = ExpsearchNsortStrLb(idbuf, nsorted_xidbox, max_xid_blen, orig_sample_ct, block_start_idx + 1);
    }
  } else {
    // --king-table-require[-xor]
    while (1) {
      for (; idx1 != block_end_idx; ++idx1) {
        // idx1 >= idx2.
        uint32_t cur_pair_ct = 0;
        const uint32_t sample_uidx1 = xid_map[idx1];
        if (IsSet(sample_require, sample_uidx1)) {
          if (!require_xor) {
            cur_pair_ct = idx1 - idx2;
            uint32_t idx2_stop = idx1;
            if (pair_idx_stop - pair_idx < cur_pair_ct) {
              cur_pair_ct = pair_idx_stop - pair_idx;
              idx2_stop = idx2 + cur_pair_ct;
            }
            if (pair_idx < pair_idx_start) {
              const uint64_t skip_ct = pair_idx_start - pair_idx;
              if (skip_ct >= cur_pair_ct) {
                idx2 = idx2_stop;
              } else {
                idx2 += skip_ct;
              }
            }
            for (; idx2 != idx2_stop; ++idx2) {
              const uint32_t sample_uidx2 = xid_map[idx2];
              *loaded_sample_idx_pairs_iter++ = sample_uidx1;
              *loaded_sample_idx_pairs_iter++ = sample_uidx2;
            }
          } else {
            // --king-table-require-xor
            cur_pair_ct += idx1 - idx2;
            for (uint32_t idx = idx2; idx != idx1; ++idx) {
              const uint32_t sample_uidx = xid_map[idx];
              cur_pair_ct -= IsSet(sample_require, sample_uidx);
            }
            uint32_t idx2_stop = idx1;
            if (pair_idx_stop - pair_idx < cur_pair_ct) {
              cur_pair_ct = pair_idx_stop - pair_idx;
              idx2_stop = idx2;
              uint32_t remaining_required_pair_ct = cur_pair_ct;
              while (1) {
                const uint32_t sample_uidx = xid_map[idx2_stop];
                ++idx2_stop;
                if (!IsSet(sample_require, sample_uidx)) {
                  --remaining_required_pair_ct;
                  if (!remaining_required_pair_ct) {
                    break;
                  }
                }
              }
            }
            if (pair_idx < pair_idx_start) {
              const uint64_t skip_ct = pair_idx_start - pair_idx;
              if (skip_ct >= cur_pair_ct) {
                idx2 = idx2_stop;
              } else {
                uint32_t remaining_skip_ct = skip_ct;
                while (1) {
                  const uint32_t sample_uidx = xid_map[idx2];
                  ++idx2;
                  if (!IsSet(sample_require, sample_uidx)) {
                    --remaining_skip_ct;
                    if (!remaining_skip_ct) {
                      break;
                    }
                  }
                }
              }
            }
            for (; idx2 != idx2_stop; ++idx2) {
              const uint32_t sample_uidx2 = xid_map[idx2];
              if (!IsSet(sample_require, sample_uidx2)) {
                *loaded_sample_idx_pairs_iter++ = sample_uidx1;
                *loaded_sample_idx_pairs_iter++ = sample_uidx2;
              }
            }
          }
        } else {
          // !IsSet(sample_require, sample_uidx1)
          for (uint32_t idx = idx2; idx != idx1; ++idx) {
            const uint32_t sample_uidx = xid_map[idx];
            cur_pair_ct += IsSet(sample_require, sample_uidx);
          }
          uint32_t idx2_stop = idx1;
          if (pair_idx_stop - pair_idx < cur_pair_ct) {
            cur_pair_ct = pair_idx_stop - pair_idx;
            idx2_stop = idx2;
            uint32_t remaining_required_pair_ct = cur_pair_ct;
            while (1) {
              const uint32_t sample_uidx = xid_map[idx2_stop];
              ++idx2_stop;
              if (IsSet(sample_require, sample_uidx)) {
                --remaining_required_pair_ct;
                if (!remaining_required_pair_ct) {
                  break;
                }
              }
            }
          }
          if (pair_idx < pair_idx_start) {
            const uint64_t skip_ct = pair_idx_start - pair_idx;
            if (skip_ct >= cur_pair_ct) {
              idx2 = idx2_stop;
            } else {
              uint32_t remaining_skip_ct = skip_ct;
              while (1) {
                const uint32_t sample_uidx = xid_map[idx2];
                ++idx2;
                if (IsSet(sample_require, sample_uidx)) {
                  --remaining_skip_ct;
                  if (!remaining_skip_ct) {
                    break;
                  }
                }
              }
            }
          }
          for (; idx2 != idx2_stop; ++idx2) {
            const uint32_t sample_uidx2 = xid_map[idx2];
            if (IsSet(sample_require, sample_uidx2)) {
              *loaded_sample_idx_pairs_iter++ = sample_uidx1;
              *loaded_sample_idx_pairs_iter++ = sample_uidx2;
            }
          }
        }
        pair_idx += cur_pair_ct;
        // bugfix (21 Jun 2023): forgot to handle the usual
        // rel_or_concordance_check==0 case.
        // (didn't notice since I was testing on IID-only data.)
        if (pair_idx == pair_idx_stop) {
          if (is_first_parallel_scan) {
            if (!rel_or_concordance_check) {
              const uint64_t orig_sample_ct64 = orig_sample_ct;
              const uintptr_t orig_sample_ctl = BitCtToWordCt(orig_sample_ct);
              const uint64_t optional_ct = orig_sample_ct - PopcountWords(sample_require, orig_sample_ctl);
              if (!require_xor) {
                pair_idx = ((orig_sample_ct64 * (orig_sample_ct - 1)) - (optional_ct * (optional_ct - 1))) / 2;
              } else {
                pair_idx = orig_sample_ct64 * optional_ct;
              }
            } else {
              pair_idx = CountRcCheckPairsSampleRequire(nsorted_xidbox, xid_map, sample_require, rel_or_concordance_check, max_xid_blen, orig_sample_ct, require_xor, idbuf);
            }
          }
          goto GetRelCheckOrKTRequirePairs_early_exit;
        }
        idx2 = block_start_idx;
      }
      block_start_idx = block_end_idx;
      if (block_start_idx == orig_sample_ct) {
        break;
      }
      idx2 = block_start_idx;
      if (!rel_or_concordance_check) {
        block_end_idx = orig_sample_ct;
      } else {
        const char* fid_start = &(nsorted_xidbox[block_start_idx * max_xid_blen]);
        const char* fid_end = AdvToDelim(fid_start, '\t');
        uint32_t id_slen;
        if (rel_or_concordance_check == kRcCheckRel) {
          id_slen = fid_end - fid_start;
        } else {
          const char* iid_end = AdvToDelim(&(fid_end[1]), '\t');
          id_slen = iid_end - fid_start;
        }
        memcpy(idbuf, fid_start, id_slen);
        idbuf[id_slen] = ' ';
        idbuf[id_slen + 1] = '\0';
        block_end_idx = ExpsearchNsortStrLb(idbuf, nsorted_xidbox, max_xid_blen, orig_sample_ct, block_start_idx + 1);
      }
    }
  }
 GetRelCheckOrKTRequirePairs_early_exit:
  *pair_idx_ptr = pair_idx;
  fpip->block_start_idx = block_start_idx;
  fpip->block_end_idx = block_end_idx;
  fpip->idx1 = idx1;
  fpip->idx2 = idx2;
}

PglErr CalcKingTableSubset(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const char* subset_fname, const char* require_fnames, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, double king_table_filter, double king_table_subset_thresh, RelConcordanceCheckMode rel_or_concordance_check, KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  // subset_fname permitted to be nullptr when rel_or_concordance_check is
  // nonzero.
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  CompressStreamState css;
  ThreadGroup tg;
  PreinitTextStream(&txs);
  PreinitCstream(&css);
  PreinitThreads(&tg);
  {
    if (unlikely(IsSet(cip->haploid_mask, 0))) {
      logerrputs("Error: --make-king-table cannot be used on haploid genomes.\n");
      goto CalcKingTableSubset_ret_INCONSISTENT_INPUT;
    }
    reterr = ConditionalAllocateNonAutosomalVariants(cip, "--make-king-table", raw_variant_ct, &variant_include, &variant_ct);
    if (unlikely(reterr)) {
      goto CalcKingTableSubset_ret_1;
    }

    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uintptr_t* sample_require = nullptr;
    const uint32_t require_xor = (king_flags / kfKingTableRequireXor) & 1;
    if (require_fnames) {
      uintptr_t* sample_require_tmp;
      reterr = LoadSampleIds(require_fnames, orig_sample_include, siip, "make-king-table", raw_sample_ct, orig_sample_ct, kfLoadSampleIdsMultifile, &sample_require_tmp, nullptr);
      if (unlikely(reterr)) {
        goto CalcKingTableSubset_ret_1;
      }
      sample_require = sample_require_tmp;
      const uint32_t require_set_size = PopcountWords(sample_require, raw_sample_ctl);
      if (unlikely(!require_set_size)) {
        logerrputs("Error: No sample ID(s) in --king-table-require[-xor] file(s) are in the current\ndataset.\n");
        goto CalcKingTableSubset_ret_INCONSISTENT_INPUT;
      }
      if (unlikely(require_xor && (require_set_size == orig_sample_ct))) {
        logerrputs("Error: All sample ID(s) in --king-table-require-xor file(s) are in the current\ndataset.\n");
        goto CalcKingTableSubset_ret_INCONSISTENT_INPUT;
      }
      logprintf("--king-table-require%s: %u sample ID%s loaded.\n", require_xor? "-xor" : "", require_set_size, (require_set_size == 1)? "" : "s");
    }
    // 1. Write output header line if necessary.
    // 2. Count number of relevant sample pairs (higher uidx in high 32 bits),
    //    and load as much as may be useful during first pass (usually there
    //    will be only one pass).
    // 3. If list is empty, error out.
    // 4. If --parallel, discard part of the list, then exit if remainder
    //    empty.
    // 5. If remainder of list is too large to process in one pass, determine
    //    number of necessary passes.  If output filename refers to the same
    //    thing as input file, append ~ to input filename.
    // Loop:
    // * Determine which sample indexes appear in this part of the list.
    //   Compute current cumulative_popcounts, perform uidx -> idx conversion.
    //   (Don't bother sorting the pairs, since that prevents
    //   --parallel/multipass mode from delivering the same results.)
    // * Execute usual KING-robust computation, write .kin0 entries.
    // * If not last pass, reload input .kin0, etc.
    //
    // Could store the pairs in a more compact manner, but can live with 50%
    // space bloat for now.
    uint32_t sample_ctaw = BitCtToAlignedWordCt(orig_sample_ct);
    uint32_t sample_ctaw2 = NypCtToAlignedWordCt(orig_sample_ct);
    uint32_t king_bufsizew = kKingMultiplexWords * orig_sample_ct;
    uintptr_t* cur_sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* loadbuf;
    uintptr_t* splitbuf_hom;
    uintptr_t* splitbuf_ref2het;
    VecW* vecaligned_buf;
    // ok if allocations are a bit oversized
    CalcKingTableSubsetCtx ctx;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &cur_sample_include) ||
                 bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_w(sample_ctaw2, &loadbuf) ||
                 bigstack_alloc_w(kPglBitTransposeBatch * sample_ctaw, &splitbuf_hom) ||
                 bigstack_alloc_w(kPglBitTransposeBatch * sample_ctaw, &splitbuf_ref2het) ||
                 bigstack_alloc_w(king_bufsizew, &(ctx.smaj_hom[0])) ||
                 bigstack_alloc_w(king_bufsizew, &(ctx.smaj_ref2het[0])) ||
                 bigstack_alloc_w(king_bufsizew, &(ctx.smaj_hom[1])) ||
                 bigstack_alloc_w(king_bufsizew, &(ctx.smaj_ref2het[1])) ||
                 bigstack_alloc_v(kPglBitTransposeBufvecs, &vecaligned_buf))) {
      goto CalcKingTableSubset_ret_NOMEM;
    }
    SetKingTableFname(king_flags, parallel_idx, parallel_tot, outname_end);
    if (subset_fname) {
      uint32_t fname_slen;
#ifdef _WIN32
      fname_slen = GetFullPathName(subset_fname, kPglFnamesize, g_textbuf, nullptr);
      if (unlikely((!fname_slen) || (fname_slen > kPglFnamesize)))
#else
      if (unlikely(!realpath(subset_fname, g_textbuf)))
#endif
      {
        logerrprintfww(kErrprintfFopen, subset_fname, strerror(errno));
        goto CalcKingTableSubset_ret_OPEN_FAIL;
      }
      if (RealpathIdentical(outname, g_textbuf, &(g_textbuf[kPglFnamesize + 64]))) {
        logerrputs("Warning: --king-table-subset input filename matches --make-king-table output\nfilename.  Appending '~' to input filename.\n");
        fname_slen = strlen(subset_fname);
        memcpy(g_textbuf, subset_fname, fname_slen);
        strcpy_k(&(g_textbuf[fname_slen]), "~");
        if (unlikely(rename(subset_fname, g_textbuf))) {
          logerrprintf("Error: Failed to append '~' to --king-table-subset input filename: %s.\n", strerror(errno));
          goto CalcKingTableSubset_ret_OPEN_FAIL;
        }
        subset_fname = g_textbuf;
      }
    }

    // Safe to "write" the header line now, if necessary.
    reterr = InitCstreamAlloc(outname, 0, king_flags & kfKingTableZs, max_thread_ct, kMaxMediumLine + kCompressStreamBlock, &css, &cswritep);
    if (unlikely(reterr)) {
      goto CalcKingTableSubset_ret_1;
    }
    const uint32_t king_col_fid = FidColIsRequired(siip, king_flags / kfKingColMaybefid);
    const uint32_t king_col_sid = SidColIsRequired(siip->sids, king_flags / kfKingColMaybesid);
    if (!parallel_idx) {
      cswritep = AppendKingTableHeader(king_flags, king_col_fid, king_col_sid, cswritep);
    }
    const uintptr_t max_sample_fmtid_blen = GetMaxSampleFmtidBlen(siip, king_col_fid, king_col_sid);
    char* collapsed_sample_fmtids;
    if (unlikely(bigstack_alloc_c(max_sample_fmtid_blen * orig_sample_ct, &collapsed_sample_fmtids))) {
      goto CalcKingTableSubset_ret_NOMEM;
    }
    // possible todo: allow this to change between passes
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > orig_sample_ct / 32) {
      calc_thread_ct = orig_sample_ct / 32;
    }
    if (!calc_thread_ct) {
      calc_thread_ct = 1;
    }
    // could eventually have 64-bit ctx.thread_start?
    if (unlikely(SetThreadCt(calc_thread_ct, &tg) ||
                 bigstack_alloc_u32(calc_thread_ct + 1, &ctx.thread_start))) {
      goto CalcKingTableSubset_ret_NOMEM;
    }

    uintptr_t line_idx = 0;
    uint32_t kinship_skip = 0;
    // we overwrite this in subset_fname case, so no need to check for
    // --strict-sid0
    XidMode xid_mode = siip->sids? kfXidModeFidIidSid : kfXidModeIidSid;
    if (subset_fname) {
      reterr = InitTextStream(subset_fname, kTextStreamBlenFast, 1, &txs);
      if (unlikely(reterr)) {
        if (reterr == kPglRetEof) {
          logerrputs("Error: Empty --king-table-subset file.\n");
          goto CalcKingTableSubset_ret_MALFORMED_INPUT;
        }
        goto CalcKingTableSubset_ret_TSTREAM_FAIL;
      }
      ++line_idx;
      const char* linebuf_iter = TextGet(&txs);
      if (unlikely(!linebuf_iter)) {
        if (!TextStreamErrcode2(&txs, &reterr)) {
          logerrputs("Error: Empty --king-table-subset file.\n");
          goto CalcKingTableSubset_ret_MALFORMED_INPUT;
        }
        goto CalcKingTableSubset_ret_TSTREAM_FAIL;
      }
      const char* token_end = CurTokenEnd(linebuf_iter);
      uint32_t token_slen = token_end - linebuf_iter;
      // Make this work with both KING- and plink2-generated .kin0 files.
      uint32_t fid_present = strequal_k(linebuf_iter, "#FID1", token_slen) || strequal_k(linebuf_iter, "FID", token_slen);
      if (fid_present) {
        linebuf_iter = FirstNonTspace(token_end);
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
        xid_mode = kfXidModeFidIid;
      } else {
        if (unlikely(*linebuf_iter != '#')) {
          goto CalcKingTableSubset_ret_INVALID_HEADER;
        }
        ++linebuf_iter;
        --token_slen;
        xid_mode = kfXidModeIid;
      }
      if (unlikely((!strequal_k(linebuf_iter, "ID1", token_slen)) && (!strequal_k(linebuf_iter, "IID1", token_slen)))) {
        goto CalcKingTableSubset_ret_INVALID_HEADER;
      }
      linebuf_iter = FirstNonTspace(token_end);
      token_end = CurTokenEnd(linebuf_iter);
      token_slen = token_end - linebuf_iter;
      if (strequal_k(linebuf_iter, "SID1", token_slen)) {
        if (siip->sids) {
          xid_mode = fid_present? kfXidModeFidIidSid : kfXidModeIidSid;
        } else {
          xid_mode |= kfXidModeFlagSkipSid;
        }
        linebuf_iter = FirstNonTspace(token_end);
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
      }
      if (fid_present) {
        if (unlikely(!strequal_k(linebuf_iter, "FID2", token_slen))) {
          goto CalcKingTableSubset_ret_INVALID_HEADER;
        }
        linebuf_iter = FirstNonTspace(token_end);
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
      }
      if (unlikely((!strequal_k(linebuf_iter, "ID2", token_slen)) && (!strequal_k(linebuf_iter, "IID2", token_slen)))) {
        goto CalcKingTableSubset_ret_INVALID_HEADER;
      }
      if (xid_mode & (kfXidModeFlagSid | kfXidModeFlagSkipSid)) {
        linebuf_iter = FirstNonTspace(token_end);
        token_end = CurTokenEnd(linebuf_iter);
        token_slen = token_end - linebuf_iter;
        if (unlikely(!strequal_k(linebuf_iter, "SID2", token_slen))) {
          goto CalcKingTableSubset_ret_INVALID_HEADER;
        }
      }
      if (king_table_subset_thresh != -DBL_MAX) {
        king_table_subset_thresh *= 1.0 - kSmallEpsilon;
        while (1) {
          linebuf_iter = FirstNonTspace(token_end);
          token_end = CurTokenEnd(linebuf_iter);
          token_slen = token_end - linebuf_iter;
          if (unlikely(!token_slen)) {
            logerrputs("Error: No kinship-coefficient column in --king-table-subset file.\n");
            goto CalcKingTableSubset_ret_INCONSISTENT_INPUT;
          }
          if (strequal_k(linebuf_iter, "KINSHIP", token_slen) || strequal_k(linebuf_iter, "Kinship", token_slen)) {
            break;
          }
          ++kinship_skip;
        }
      }
    }

    uint32_t* xid_map;  // IDs not collapsed
    char* sorted_xidbox;
    uintptr_t max_xid_blen;
    // may as well use natural-sort order in rel-check-only case
    reterr = SortedXidboxInitAlloc(orig_sample_include, siip, orig_sample_ct, xid_mode, (!subset_fname), &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto CalcKingTableSubset_ret_1;
    }
    char* idbuf;
    if (unlikely(bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto CalcKingTableSubset_ret_NOMEM;
    }

    ctx.homhom_needed = (king_flags & kfKingColNsnp) || ((!(king_flags & kfKingCounts)) && (king_flags & (kfKingColHethet | kfKingColIbs0 | kfKingColIbs1 | kfKingColHamming)));
    const uint32_t homhom_needed_p4 = ctx.homhom_needed + 4;
    // if homhom_needed, 8 + 20 bytes per pair, otherwise 8 + 16
    uintptr_t pair_buf_capacity = bigstack_left();
    if (unlikely(pair_buf_capacity < 2 * kCacheline)) {
      goto CalcKingTableSubset_ret_NOMEM;
    }
    // adverse rounding
    pair_buf_capacity = (pair_buf_capacity - 2 * kCacheline) / (24 + 4 * ctx.homhom_needed);
    if (pair_buf_capacity > 0xffffffffU) {
      // 32-bit ctx.thread_start[] for now
      pair_buf_capacity = 0xffffffffU;
    }
    ctx.loaded_sample_idx_pairs = S_CAST(uint32_t*, bigstack_alloc_raw_rd(pair_buf_capacity * 2 * sizeof(int32_t)));
    ctx.king_counts = R_CAST(uint32_t*, g_bigstack_base);
    SetThreadFuncAndData(CalcKingTableSubsetThread, &ctx, &tg);

    FidPairIterator fpi;
    InitFidPairIterator(&fpi);

    uint64_t pair_idx = 0;
    if (!subset_fname) {
      GetRelCheckOrKTRequirePairs(sorted_xidbox, xid_map, sample_require, max_xid_blen, orig_sample_ct, rel_or_concordance_check, require_xor, (parallel_tot != 1), 0, pair_buf_capacity, &fpi, &pair_idx, ctx.loaded_sample_idx_pairs, idbuf);
    } else {
      fputs("Scanning --king-table-subset file...", stdout);
      fflush(stdout);
      reterr = KingTableSubsetLoad(sorted_xidbox, xid_map, sample_require, max_xid_blen, orig_sample_ct, king_table_subset_thresh, xid_mode, rel_or_concordance_check, require_xor, kinship_skip, (parallel_tot != 1), 0, pair_buf_capacity, line_idx, &txs, &pair_idx, ctx.loaded_sample_idx_pairs, idbuf);
      if (unlikely(reterr)) {
        goto CalcKingTableSubset_ret_1;
      }
    }
    uint64_t pair_idx_global_start = 0;
    uint64_t pair_idx_global_stop = ~0LLU;
    if (parallel_tot != 1) {
      const uint64_t parallel_pair_ct = pair_idx;
      pair_idx_global_start = (parallel_idx * parallel_pair_ct) / parallel_tot;
      pair_idx_global_stop = ((parallel_idx + 1) * parallel_pair_ct) / parallel_tot;
      if (pair_idx > pair_buf_capacity) {
        // may as well document possible overflow
        if (unlikely(parallel_pair_ct > ((~0LLU) / kParallelMax))) {
          if (!subset_fname) {
            // This is easy to support if there's ever a need, of course.
            logerrputs("Error: Too many \"--make-king-table rel-check\" sample pairs for this " PROG_NAME_STR "\nbuild.\n");
          } else {
            logerrputs("Error: Too many --king-table-subset sample pairs for this " PROG_NAME_STR " build.\n");
          }
          reterr = kPglRetNotYetSupported;
          goto CalcKingTableSubset_ret_1;
        }
        if (pair_idx_global_stop > pair_buf_capacity) {
          // large --parallel job
          pair_idx = 0;
          if (!subset_fname) {
            InitFidPairIterator(&fpi);
            GetRelCheckOrKTRequirePairs(sorted_xidbox, xid_map, sample_require, max_xid_blen, orig_sample_ct, rel_or_concordance_check, require_xor, 0, pair_idx_global_start, MINV(pair_idx_global_stop, pair_idx_global_start + pair_buf_capacity), &fpi, &pair_idx, ctx.loaded_sample_idx_pairs, idbuf);
          } else {
            reterr = TextRewind(&txs);
            if (unlikely(reterr)) {
              goto CalcKingTableSubset_ret_TSTREAM_FAIL;
            }
            // bugfix (4 Oct 2019): forgot a bunch of reinitialization here
            line_idx = 1;
            char* header_throwaway;
            reterr = TextNextLineLstrip(&txs, &header_throwaway);
            if (unlikely(reterr)) {
              goto CalcKingTableSubset_ret_TSTREAM_REWIND_FAIL;
            }
            reterr = KingTableSubsetLoad(sorted_xidbox, xid_map, sample_require, max_xid_blen, orig_sample_ct, king_table_subset_thresh, xid_mode, rel_or_concordance_check, require_xor, kinship_skip, 0, pair_idx_global_start, MINV(pair_idx_global_stop, pair_idx_global_start + pair_buf_capacity), line_idx, &txs, &pair_idx, ctx.loaded_sample_idx_pairs, idbuf);
            if (unlikely(reterr)) {
              goto CalcKingTableSubset_ret_1;
            }
          }
        } else {
          pair_idx = pair_idx_global_stop;
          if (pair_idx_global_start) {
            memmove(ctx.loaded_sample_idx_pairs, &(ctx.loaded_sample_idx_pairs[pair_idx_global_start * 2]), (pair_idx_global_stop - pair_idx_global_start) * 2 * sizeof(int32_t));
          }
        }
      } else {
        pair_idx = pair_idx_global_stop;
        if (pair_idx_global_start) {
          memmove(ctx.loaded_sample_idx_pairs, &(ctx.loaded_sample_idx_pairs[pair_idx_global_start * 2]), (pair_idx_global_stop - pair_idx_global_start) * 2 * sizeof(int32_t));
        }
      }
    }
    uint64_t pair_idx_cur_start = pair_idx_global_start;
    uint64_t king_table_filter_ct = 0;
    uintptr_t pass_idx = 1;
    while (pair_idx_cur_start < pair_idx) {
      ZeroWArr(raw_sample_ctl, cur_sample_include);
      const uintptr_t cur_pair_ct = pair_idx - pair_idx_cur_start;
      const uintptr_t cur_pair_ct_x2 = 2 * cur_pair_ct;
      for (uintptr_t ulii = 0; ulii != cur_pair_ct_x2; ++ulii) {
        SetBit(ctx.loaded_sample_idx_pairs[ulii], cur_sample_include);
      }
      FillCumulativePopcounts(cur_sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      const uint32_t cur_sample_ct = sample_include_cumulative_popcounts[raw_sample_ctl - 1] + PopcountWord(cur_sample_include[raw_sample_ctl - 1]);
      const uint32_t cur_sample_ctaw = BitCtToAlignedWordCt(cur_sample_ct);
      const uint32_t cur_sample_ctaw2 = NypCtToAlignedWordCt(cur_sample_ct);
      if (cur_sample_ct != raw_sample_ct) {
        for (uintptr_t ulii = 0; ulii != cur_pair_ct_x2; ++ulii) {
          ctx.loaded_sample_idx_pairs[ulii] = RawToSubsettedPos(cur_sample_include, sample_include_cumulative_popcounts, ctx.loaded_sample_idx_pairs[ulii]);
        }
      }
      ZeroU32Arr(cur_pair_ct * homhom_needed_p4, ctx.king_counts);
      CollapsedSampleFmtidInit(cur_sample_include, siip, cur_sample_ct, king_col_fid, king_col_sid, max_sample_fmtid_blen, collapsed_sample_fmtids);
      for (uint32_t tidx = 0; tidx <= calc_thread_ct; ++tidx) {
        ctx.thread_start[tidx] = (tidx * S_CAST(uint64_t, cur_pair_ct)) / calc_thread_ct;
      }
      if (pass_idx != 1) {
        ReinitThreads(&tg);
      }
      // possible todo: singleton/monomorphic optimization for sufficiently
      // large jobs
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      uint32_t variants_completed = 0;
      uint32_t parity = 0;
      const uint32_t sample_batch_ct_m1 = (cur_sample_ct - 1) / kPglBitTransposeBatch;
      PgrSampleSubsetIndex pssi;
      PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
      do {
        const uint32_t cur_block_size = MINV(variant_ct - variants_completed, kKingMultiplex);
        uintptr_t* cur_smaj_hom = ctx.smaj_hom[parity];
        uintptr_t* cur_smaj_ref2het = ctx.smaj_ref2het[parity];
        // "block" = distance computation granularity, usually 1024 or 1536
        //           variants
        // "batch" = variant-major-to-sample-major transpose granularity,
        //           currently 512 variants
        uint32_t variant_batch_size = kPglBitTransposeBatch;
        uint32_t variant_batch_size_rounded_up = kPglBitTransposeBatch;
        const uint32_t write_batch_ct_m1 = (cur_block_size - 1) / kPglBitTransposeBatch;
        for (uint32_t write_batch_idx = 0; ; ++write_batch_idx) {
          if (write_batch_idx >= write_batch_ct_m1) {
            if (write_batch_idx > write_batch_ct_m1) {
              break;
            }
            variant_batch_size = ModNz(cur_block_size, kPglBitTransposeBatch);
            variant_batch_size_rounded_up = variant_batch_size;
            const uint32_t variant_batch_size_rem = variant_batch_size % kBitsPerWord;
            if (variant_batch_size_rem) {
              const uint32_t trailing_variant_ct = kBitsPerWord - variant_batch_size_rem;
              variant_batch_size_rounded_up += trailing_variant_ct;
              ZeroWArr(trailing_variant_ct * cur_sample_ctaw, &(splitbuf_hom[variant_batch_size * cur_sample_ctaw]));
              ZeroWArr(trailing_variant_ct * cur_sample_ctaw, &(splitbuf_ref2het[variant_batch_size * cur_sample_ctaw]));
            }
          }
          uintptr_t* hom_iter = splitbuf_hom;
          uintptr_t* ref2het_iter = splitbuf_ref2het;
          for (uint32_t uii = 0; uii != variant_batch_size; ++uii) {
            const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
            reterr = PgrGet(cur_sample_include, pssi, cur_sample_ct, variant_uidx, simple_pgrp, loadbuf);
            if (unlikely(reterr)) {
              PgenErrPrintNV(reterr, variant_uidx);
              goto CalcKingTableSubset_ret_1;
            }
            // may want to support some sort of low-MAF optimization here
            SetTrailingNyps(cur_sample_ct, loadbuf);
            SplitHomRef2hetUnsafeW(loadbuf, cur_sample_ctaw2, hom_iter, ref2het_iter);
            hom_iter = &(hom_iter[cur_sample_ctaw]);
            ref2het_iter = &(ref2het_iter[cur_sample_ctaw]);
          }
          // uintptr_t* read_iter = loadbuf;
          uintptr_t* write_hom_iter = &(cur_smaj_hom[write_batch_idx * kPglBitTransposeWords]);
          uintptr_t* write_ref2het_iter = &(cur_smaj_ref2het[write_batch_idx * kPglBitTransposeWords]);
          uint32_t write_batch_size = kPglBitTransposeBatch;
          for (uint32_t sample_batch_idx = 0; ; ++sample_batch_idx) {
            if (sample_batch_idx >= sample_batch_ct_m1) {
              if (sample_batch_idx > sample_batch_ct_m1) {
                break;
              }
              write_batch_size = ModNz(cur_sample_ct, kPglBitTransposeBatch);
            }
            // bugfix: read_batch_size must be rounded up to word boundary,
            // since we want to one-out instead of zero-out the trailing bits
            //
            // bugfix: if we always use kPglBitTransposeBatch instead of
            // variant_batch_size_rounded_up, we read/write past the
            // kKingMultiplex limit and clobber the first variants of the next
            // sample with garbage.
            TransposeBitblock(&(splitbuf_hom[sample_batch_idx * kPglBitTransposeWords]), cur_sample_ctaw, kKingMultiplexWords, variant_batch_size_rounded_up, write_batch_size, write_hom_iter, vecaligned_buf);
            TransposeBitblock(&(splitbuf_ref2het[sample_batch_idx * kPglBitTransposeWords]), cur_sample_ctaw, kKingMultiplexWords, variant_batch_size_rounded_up, write_batch_size, write_ref2het_iter, vecaligned_buf);
            write_hom_iter = &(write_hom_iter[kKingMultiplex * kPglBitTransposeWords]);
            write_ref2het_iter = &(write_ref2het_iter[kKingMultiplex * kPglBitTransposeWords]);
          }
        }
        const uint32_t cur_block_sizew = BitCtToWordCt(cur_block_size);
        if (cur_block_sizew < kKingMultiplexWords) {
          uintptr_t* write_hom_iter = &(cur_smaj_hom[cur_block_sizew]);
          uintptr_t* write_ref2het_iter = &(cur_smaj_ref2het[cur_block_sizew]);
          const uint32_t write_word_ct = kKingMultiplexWords - cur_block_sizew;
          for (uint32_t sample_idx = 0; sample_idx != cur_sample_ct; ++sample_idx) {
            ZeroWArr(write_word_ct, write_hom_iter);
            ZeroWArr(write_word_ct, write_ref2het_iter);
            write_hom_iter = &(write_hom_iter[kKingMultiplexWords]);
            write_ref2het_iter = &(write_ref2het_iter[kKingMultiplexWords]);
          }
        }
        if (variants_completed) {
          JoinThreads(&tg);
          // CalcKingTableSubsetThread() never errors out
        }
        // this update must occur after JoinThreads() call
        if (variants_completed + cur_block_size == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto CalcKingTableSubset_ret_THREAD_CREATE_FAIL;
        }
        printf("\r--make-king-table pass %" PRIuPTR ": %u variants complete.", pass_idx, variants_completed);
        fflush(stdout);
        variants_completed += cur_block_size;
        parity = 1 - parity;
      } while (!IsLastBlock(&tg));
      JoinThreads(&tg);
      printf("\r--make-king-table pass %" PRIuPTR ": Writing...                   \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", pass_idx);
      fflush(stdout);

      const uint32_t king_col_id = king_flags & kfKingColId;
      const uint32_t king_col_nsnp = king_flags & kfKingColNsnp;
      const uint32_t king_col_hethet = king_flags & kfKingColHethet;
      const uint32_t king_col_ibs0 = king_flags & kfKingColIbs0;
      const uint32_t king_col_ibs1 = king_flags & kfKingColIbs1;
      const uint32_t king_col_hamming = king_flags & kfKingColHamming;
      const uint32_t king_col_kinship = king_flags & kfKingColKinship;
      const uint32_t report_counts = king_flags & kfKingCounts;
      uint32_t* results_iter = ctx.king_counts;
      double nonmiss_recip = 0.0;
      for (uintptr_t cur_pair_idx = 0; cur_pair_idx != cur_pair_ct; ++cur_pair_idx, results_iter = &(results_iter[homhom_needed_p4])) {
        const uint32_t ibs0_ct = results_iter[kKingOffsetIbs0];
        const uint32_t hethet_ct = results_iter[kKingOffsetHethet];
        const uint32_t het2hom1_ct = results_iter[kKingOffsetHet2Hom1];
        const uint32_t het1hom2_ct = results_iter[kKingOffsetHet1Hom2];
        const intptr_t smaller_het_ct = hethet_ct + MINV(het1hom2_ct, het2hom1_ct);
        const double kinship_coeff = 0.5 - (S_CAST(double, 4 * S_CAST(intptr_t, ibs0_ct) + het1hom2_ct + het2hom1_ct) / S_CAST(double, 4 * smaller_het_ct));
        if ((king_table_filter != -DBL_MAX) && (kinship_coeff < king_table_filter)) {
          ++king_table_filter_ct;
          continue;
        }
        const uint32_t sample_idx1 = ctx.loaded_sample_idx_pairs[2 * cur_pair_idx];
        const uint32_t sample_idx2 = ctx.loaded_sample_idx_pairs[2 * cur_pair_idx + 1];
        if (king_col_id) {
          cswritep = strcpyax(cswritep, &(collapsed_sample_fmtids[max_sample_fmtid_blen * sample_idx1]), '\t');
          cswritep = strcpyax(cswritep, &(collapsed_sample_fmtids[max_sample_fmtid_blen * sample_idx2]), '\t');
        }
        if (homhom_needed_p4 == 5) {
          const uint32_t homhom_ct = results_iter[kKingOffsetHomhom];
          const uint32_t nonmiss_ct = het1hom2_ct + het2hom1_ct + homhom_ct + hethet_ct;
          if (king_col_nsnp) {
            cswritep = u32toa_x(nonmiss_ct, '\t', cswritep);
          }
          if (!report_counts) {
            nonmiss_recip = 1.0 / u31tod(nonmiss_ct);
          }
        }
        if (king_col_hethet) {
          if (report_counts) {
            cswritep = u32toa(hethet_ct, cswritep);
          } else {
            cswritep = dtoa_g(nonmiss_recip * u31tod(hethet_ct), cswritep);
          }
          *cswritep++ = '\t';
        }
        if (king_col_ibs0) {
          if (report_counts) {
            cswritep = u32toa(ibs0_ct, cswritep);
          } else {
            cswritep = dtoa_g(nonmiss_recip * u31tod(ibs0_ct), cswritep);
          }
          *cswritep++ = '\t';
        }
        if (king_col_ibs1) {
          if (report_counts) {
            cswritep = u32toa_x(het1hom2_ct, '\t', cswritep);
            cswritep = u32toa(het2hom1_ct, cswritep);
          } else {
            cswritep = dtoa_g(nonmiss_recip * u31tod(het1hom2_ct), cswritep);
            *cswritep++ = '\t';
            cswritep = dtoa_g(nonmiss_recip * u31tod(het2hom1_ct), cswritep);
          }
          *cswritep++ = '\t';
        }
        if (king_col_hamming) {
          const uint32_t hamming_dist = 2 * ibs0_ct + het1hom2_ct + het2hom1_ct;
          if (report_counts) {
            cswritep = u32toa(hamming_dist, cswritep);
          } else {
            cswritep = dtoa_g(nonmiss_recip * 0.5 * S_CAST(double, hamming_dist), cswritep);
          }
          *cswritep++ = '\t';
        }
        if (king_col_kinship) {
          cswritep = dtoa_g(kinship_coeff, cswritep);
          ++cswritep;
        }
        DecrAppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto CalcKingTableSubset_ret_WRITE_FAIL;
        }
      }

      putc_unlocked('\r', stdout);
      const uint64_t pair_complete_ct = pair_idx - pair_idx_global_start;
      logprintf("Subsetted --make-king-table: %" PRIu64 " pair%s complete.\n", pair_complete_ct, (pair_complete_ct == 1)? "" : "s");
      if (TextEof(&txs) || (pair_idx == pair_idx_global_stop)) {
        break;
      }
      pair_idx_cur_start = pair_idx;
      if (!subset_fname) {
        GetRelCheckOrKTRequirePairs(sorted_xidbox, xid_map, sample_require, max_xid_blen, orig_sample_ct, rel_or_concordance_check, require_xor, 0, pair_idx_global_start, MINV(pair_idx_global_stop, pair_idx_global_start + pair_buf_capacity), &fpi, &pair_idx, ctx.loaded_sample_idx_pairs, idbuf);
      } else {
        fputs("Scanning --king-table-subset file...", stdout);
        fflush(stdout);
        reterr = KingTableSubsetLoad(sorted_xidbox, xid_map, sample_require, max_xid_blen, orig_sample_ct, king_table_subset_thresh, xid_mode, rel_or_concordance_check, require_xor, kinship_skip, 0, pair_idx_cur_start, MINV(pair_idx_global_stop, pair_idx_cur_start + pair_buf_capacity), line_idx, &txs, &pair_idx, ctx.loaded_sample_idx_pairs, idbuf);
        if (unlikely(reterr)) {
          goto CalcKingTableSubset_ret_1;
        }
      }
      ++pass_idx;
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto CalcKingTableSubset_ret_WRITE_FAIL;
    }
    logprintfww("Results written to %s .\n", outname);
    if (king_table_filter != -DBL_MAX) {
      const uint64_t reported_ct = pair_idx - pair_idx_global_start - king_table_filter_ct;
      logprintf("--king-table-filter: %" PRIu64 " relationship%s reported (%" PRIu64 " filtered out).\n", reported_ct, (reported_ct == 1)? "" : "s", king_table_filter_ct);
    }
  }
  while (0) {
  CalcKingTableSubset_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CalcKingTableSubset_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  CalcKingTableSubset_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind("--king-table-subset file", &txs, &reterr);
    break;
  CalcKingTableSubset_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--king-table-subset file", &txs);
    break;
  CalcKingTableSubset_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  CalcKingTableSubset_ret_INVALID_HEADER:
    logerrputs("Error: Invalid header line in --king-table-subset file.\n");
  CalcKingTableSubset_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  CalcKingTableSubset_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  CalcKingTableSubset_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 CalcKingTableSubset_ret_1:
  CleanupThreads(&tg);
  CleanupTextStream2("--king-table-subset file", &txs, &reterr);
  CswriteCloseCond(&css, cswritep);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

// assumes trailing bits of genovec are zeroed out
PglErr ExpandCenteredVarmaj(const uintptr_t* genovec, const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t variance_standardize, uint32_t is_haploid, uint32_t sample_ct, uint32_t dosage_ct, double ref_freq, double* normed_dosages) {
  const double alt_freq = 1.0 - ref_freq;
  double inv_stdev;
  if (variance_standardize) {
    const double variance = 2 * ref_freq * alt_freq;
    if (!(variance > kSmallEpsilon)) {
      // See LoadMultiallelicCenteredVarmaj().  This check was tightened up in
      // alpha 3 to reject all-het and monomorphic-wrong-allele variants.
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      GenoarrCountFreqsUnsafe(genovec, sample_ct, genocounts);
      if (unlikely(dosage_ct || genocounts[1])) {
        return kPglRetDegenerateData;
      }
      if (variance != variance) {
        if (unlikely(genocounts[0] || genocounts[2])) {
          return kPglRetDegenerateData;
        }
      } else {
        if (ref_freq > 0.5) {
          if (unlikely(genocounts[2])) {
            return kPglRetDegenerateData;
          }
        } else {
          if (unlikely(genocounts[0])) {
            return kPglRetDegenerateData;
          }
        }
      }
      ZeroDArr(sample_ct, normed_dosages);
      return kPglRetSuccess;
    }
    inv_stdev = 1.0 / sqrt(variance);
    if (is_haploid) {
      // For our purposes, variance is doubled in haploid case.
      inv_stdev *= (1.0 / kSqrt2);
    }
    // possible todo:
    // * Could use one inv_stdev for males and one for nonmales for chrX
    //   --score (while still leaving that out of GRM... or just leave males
    //   out there?).  This depends on dosage compensation model; discussed in
    //   e.g. GCTA paper.
  } else {
    // Extra factor of 2 removed from haploid 'cov' formula in alpha 3.
    inv_stdev = is_haploid? 0.5 : 1.0;
  }
  PopulateRescaledDosage(genovec, dosage_present, dosage_main, inv_stdev, -2 * alt_freq * inv_stdev, 0.0, sample_ct, dosage_ct, normed_dosages);
  return kPglRetSuccess;
}

// This breaks the "don't pass pssi between functions" rule since it's a thin
// wrapper around PgrGetInv1D().
PglErr LoadBiallelicCenteredVarmaj(const uintptr_t* sample_include, PgrSampleSubsetIndex pssi, uint32_t variance_standardize, uint32_t is_haploid, uint32_t sample_ct, uint32_t variant_uidx, double ref_freq, PgenReader* simple_pgrp, uint32_t* missing_presentp, double* normed_dosages, uintptr_t* genovec_buf, uintptr_t* dosage_present_buf, Dosage* dosage_main_buf) {
  uint32_t dosage_ct;
  PglErr reterr = PgrGetD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, genovec_buf, dosage_present_buf, dosage_main_buf, &dosage_ct);
  if (unlikely(reterr)) {
    // don't print malformed-.pgen error message here, since this is called
    // from multithreaded loops
    return reterr;
  }
  ZeroTrailingNyps(sample_ct, genovec_buf);
  if (missing_presentp) {
    // missing_present assumed to be initialized to 0
    // this should probably be a library function...
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    if (!dosage_ct) {
      for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
        const uintptr_t detect_11 = Word11(genovec_buf[widx]);
        if (detect_11) {
          *missing_presentp = 1;
          break;
        }
      }
    } else {
      Halfword* dosage_present_alias = R_CAST(Halfword*, dosage_present_buf);
      for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
        const uintptr_t detect_11 = Word11(genovec_buf[widx]);
        if (detect_11) {
          if (PackWordToHalfword(detect_11) & (~dosage_present_alias[widx])) {
            *missing_presentp = 1;
            break;
          }
        }
      }
    }
  }
  return ExpandCenteredVarmaj(genovec_buf, dosage_present_buf, dosage_main_buf, variance_standardize, is_haploid, sample_ct, dosage_ct, ref_freq, normed_dosages);
}

double ComputeDiploidMultiallelicVariance(const double* cur_allele_freqs, uint32_t cur_allele_ct) {
  const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
  double variance = 0.0;
  double freq_sum = 0.0;
  for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct_m1; ++allele_idx) {
    const double cur_allele_freq = cur_allele_freqs[allele_idx];
    variance += cur_allele_freq * (1.0 - cur_allele_freq);
    freq_sum += cur_allele_freq;
  }
  if (freq_sum < 1.0 - kSmallEpsilon) {
    const double last_allele_freq = 1.0 - freq_sum;
    variance += freq_sum * last_allele_freq;
  }
  return variance;
}

// Assumes trailing bits of pgvp->genovec have been zeroed out.
BoolErr CheckMultiallelicDegenVariance(const PgenVariant* pgvp, const double* cur_allele_freqs, uint32_t sample_ct, uint32_t cur_allele_ct, double variance) {
  // One allele has 100% frequency (or all frequencies are NaN).
  // If it's the REF allele, error out unless all nonmissing genotypes are
  // homozygous-ref, in which case this row can be filled with zeroes (or
  // omitted).
  // If it's ALT1, error out unless all nonmissing genotypes are hom-ALT1, etc.
  const uintptr_t* genovec_buf = pgvp->genovec;
  STD_ARRAY_DECL(uint32_t, 4, genocounts);
  GenoarrCountFreqsUnsafe(genovec_buf, sample_ct, genocounts);
  if (unlikely(pgvp->dosage_ct || genocounts[1])) {
    return 1;
  }
  const uint32_t nm_sample_ct = genocounts[2];
  if (variance != variance) {
    // NaN frequency is possible when all founder genotypes/dosages are
    // missing.  Error out in this case unless all other genotypes/dosages are
    // also missing.
    return (genocounts[0] || nm_sample_ct);
  }
  if (cur_allele_freqs[0] > 0.5) {
    return (nm_sample_ct != 0);
  }
  if (unlikely(genocounts[0])) {
    return 1;
  }
  if (cur_allele_freqs[1] > 0.5) {
    return (pgvp->patch_10_ct != 0);
  }
  if (pgvp->patch_10_ct != nm_sample_ct) {
    return 0;
  }
  const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
  uint32_t mono_allele_idx;
  for (mono_allele_idx = 2; mono_allele_idx != cur_allele_ct_m1; ++mono_allele_idx) {
    if (cur_allele_freqs[mono_allele_idx] > 0.5) {
      break;
    }
  }
  return !AllBytesAreX(pgvp->patch_10_vals, mono_allele_idx, 2 * nm_sample_ct);
}

PglErr LoadMultiallelicCenteredVarmaj(const uintptr_t* sample_include, PgrSampleSubsetIndex pssi, const double* cur_allele_freqs, uint32_t variance_standardize, uint32_t is_haploid, uint32_t sample_ct, uint32_t variant_uidx, uint32_t cur_allele_ct, uint32_t allele_idx_start, uint32_t allele_idx_end, PgenReader* simple_pgrp, uint32_t* missing_presentp, double* normed_dosages, PgenVariant* pgvp, double* allele_1copy_buf) {
  // This handles cur_allele_ct == 2 correctly.  But we typically don't use it
  // in that case since it does ~2x as much work as necessary: the two
  // normed_dosages[] rows are identical except for opposite sign, so it's best
  // to combine them into one row.
  PglErr reterr = PgrGetMD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgvp);
  if (unlikely(reterr)) {
    return reterr;
  }
  ZeroTrailingNyps(sample_ct, pgvp->genovec);
  const uintptr_t* genovec_buf = pgvp->genovec;
  if (missing_presentp) {
    // missing_present assumed to be initialized to 0
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    if (!pgvp->dosage_ct) {
      for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
        const uintptr_t detect_11 = Word11(genovec_buf[widx]);
        if (detect_11) {
          *missing_presentp = 1;
          break;
        }
      }
    } else {
      Halfword* dosage_present_alias = R_CAST(Halfword*, pgvp->dosage_present);
      for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
        const uintptr_t detect_11 = Word11(genovec_buf[widx]);
        if (detect_11) {
          if (PackWordToHalfword(detect_11) & (~dosage_present_alias[widx])) {
            *missing_presentp = 1;
            break;
          }
        }
      }
    }
  }
  const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
  double freq_sum = cur_allele_freqs[0];
  for (uint32_t uii = 1; uii != cur_allele_ct_m1; ++uii) {
    freq_sum += cur_allele_freqs[uii];
  }
  const double last_allele_freq = 1.0 - freq_sum;
  double inv_stdev;
  if (variance_standardize) {
    const double variance = ComputeDiploidMultiallelicVariance(cur_allele_freqs, cur_allele_ct);
    if (!(variance > kSmallEpsilon)) {
      if (unlikely(CheckMultiallelicDegenVariance(pgvp, cur_allele_freqs, sample_ct, cur_allele_ct, variance))) {
        return kPglRetDegenerateData;
      }
      ZeroDArr(S_CAST(uintptr_t, sample_ct) * (allele_idx_end - allele_idx_start), normed_dosages);
      return kPglRetSuccess;
    }
    inv_stdev = (1.0 / kSqrt2) / sqrt(variance);
    if (is_haploid) {
      inv_stdev *= (1.0 / kSqrt2);
    }
  } else {
    inv_stdev = is_haploid? (0.5 / kSqrt2) : (1.0 / kSqrt2);
  }
  if (!pgvp->dosage_ct) {
    // diploid:
    //   \sum_i x_i * (1 - x_i)
    double lookup_vals[32] ALIGNV16;
    double* normed_dosages0 = normed_dosages - (allele_idx_start * S_CAST(uintptr_t, sample_ct));
    double alt1_intercept = 0.0;
    for (uint32_t allele_idx = allele_idx_start; allele_idx != allele_idx_end; ++allele_idx) {
      double cur_allele_freq;
      if (allele_idx != cur_allele_ct_m1) {
        cur_allele_freq = cur_allele_freqs[allele_idx];
      } else {
        cur_allele_freq = last_allele_freq;
      }
      const double intercept = -2 * cur_allele_freq * inv_stdev;
      if (!allele_idx) {
        // genovec entry of 0 corresponds to 2 copies of REF allele, etc.
        lookup_vals[0] = intercept + 2 * inv_stdev;
        lookup_vals[2] = intercept + inv_stdev;
        lookup_vals[4] = intercept;
        lookup_vals[6] = 0.0;
        InitLookup16x8bx2(lookup_vals);
        GenoarrLookup16x8bx2(genovec_buf, lookup_vals, sample_ct, normed_dosages0);
        continue;
      }
      allele_1copy_buf[allele_idx] = intercept + inv_stdev;
      if (allele_idx == 1) {
        alt1_intercept = intercept;
        lookup_vals[0] = intercept;
        lookup_vals[2] = intercept + inv_stdev;
        lookup_vals[4] = intercept + 2 * inv_stdev;
        lookup_vals[6] = 0.0;
        InitLookup16x8bx2(lookup_vals);
        GenoarrLookup16x8bx2(genovec_buf, lookup_vals, sample_ct, &(normed_dosages0[sample_ct]));
      } else {
        // bugfix (28 May 2023): missing genotypes in multiallelic variants
        // were not handled correctly here
        lookup_vals[0] = intercept;
        lookup_vals[2] = intercept;
        lookup_vals[4] = intercept;
        lookup_vals[6] = 0.0;
        InitLookup16x8bx2(lookup_vals);
        double* normed_dosages_cur_allele = &(normed_dosages0[allele_idx * S_CAST(uintptr_t, sample_ct)]);
        GenoarrLookup16x8bx2(genovec_buf, lookup_vals, sample_ct, normed_dosages_cur_allele);
      }
    }
    const uintptr_t* patch_01_set = pgvp->patch_01_set;
    const AlleleCode* patch_01_vals = pgvp->patch_01_vals;
    const uintptr_t* patch_10_set = pgvp->patch_10_set;
    const AlleleCode* patch_10_vals = pgvp->patch_10_vals;
    const uint32_t patch_01_ct = pgvp->patch_01_ct;
    const uint32_t patch_10_ct = pgvp->patch_10_ct;
    if ((allele_idx_start < 2) && (allele_idx_end >= 2)) {
      if (patch_01_ct) {
        uintptr_t sample_idx_base = 0;
        uintptr_t cur_bits = patch_01_set[0];
        if (cur_allele_ct == allele_idx_end) {
          for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
            normed_dosages0[sample_ct + sample_idx] = alt1_intercept;
            const uintptr_t cur_allele_code = patch_01_vals[uii];
            normed_dosages0[cur_allele_code * sample_ct + sample_idx] = allele_1copy_buf[cur_allele_code];
          }
        } else {
          for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
            normed_dosages0[sample_ct + sample_idx] = alt1_intercept;
            const uintptr_t cur_allele_code = patch_01_vals[uii];
            if (cur_allele_code < allele_idx_end) {
              normed_dosages0[cur_allele_code * sample_ct + sample_idx] = allele_1copy_buf[cur_allele_code];
            }
          }
        }
      }
      if (patch_10_ct) {
        uintptr_t sample_idx_base = 0;
        uintptr_t cur_bits = patch_10_set[0];
        if (cur_allele_ct == allele_idx_end) {
          for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
            normed_dosages0[sample_ct + sample_idx] = alt1_intercept;
            const uintptr_t ac0 = patch_10_vals[2 * uii];
            const uintptr_t ac1 = patch_10_vals[2 * uii + 1];
            const double ac0_1copy_val = allele_1copy_buf[ac0];
            if (ac0 == ac1) {
              normed_dosages0[ac0 * sample_ct + sample_idx] = ac0_1copy_val + inv_stdev;
            } else {
              normed_dosages0[ac0 * sample_ct + sample_idx] = ac0_1copy_val;
              normed_dosages0[ac1 * sample_ct + sample_idx] = allele_1copy_buf[ac1];
            }
          }
        } else {
          for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
            normed_dosages0[sample_ct + sample_idx] = alt1_intercept;
            const uintptr_t ac0 = patch_10_vals[2 * uii];
            if (ac0 >= allele_idx_end) {
              continue;
            }
            const uintptr_t ac1 = patch_10_vals[2 * uii + 1];
            const double ac0_1copy_val = allele_1copy_buf[ac0];
            if (ac0 == ac1) {
              normed_dosages0[ac0 * sample_ct + sample_idx] = ac0_1copy_val + inv_stdev;
            } else {
              normed_dosages0[ac0 * sample_ct + sample_idx] = ac0_1copy_val;
              if (ac1 < allele_idx_end) {
                normed_dosages0[ac1 * sample_ct + sample_idx] = allele_1copy_buf[ac1];
              }
            }
          }
        }
      }
    } else {
      if (patch_01_ct) {
        uintptr_t sample_idx_base = 0;
        uintptr_t cur_bits = patch_01_set[0];
        for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
          const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &cur_bits);
          const uintptr_t cur_allele_code = patch_01_vals[uii];
          if ((cur_allele_code >= allele_idx_start) && (cur_allele_code < allele_idx_end)) {
            normed_dosages0[cur_allele_code * sample_ct + sample_idx] = allele_1copy_buf[cur_allele_code];
          }
        }
      }
      if (patch_10_ct) {
        uintptr_t sample_idx_base = 0;
        uintptr_t cur_bits = patch_10_set[0];
        for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
          const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &cur_bits);
          const uintptr_t ac0 = patch_10_vals[2 * uii];
          if (ac0 >= allele_idx_end) {
            continue;
          }
          const uintptr_t ac1 = patch_10_vals[2 * uii + 1];
          if (ac1 < allele_idx_start) {
            continue;
          }
          const double ac0_1copy_val = allele_1copy_buf[ac0];
          if (ac0 == ac1) {
            normed_dosages0[ac0 * sample_ct + sample_idx] = ac0_1copy_val + inv_stdev;
          } else {
            if (ac0 >= allele_idx_start) {
              normed_dosages0[ac0 * sample_ct + sample_idx] = ac0_1copy_val;
            }
            if (ac1 < allele_idx_end) {
              normed_dosages0[ac1 * sample_ct + sample_idx] = allele_1copy_buf[ac1];
            }
          }
        }
      }
    }
    return kPglRetSuccess;
  }
  fputs("true multiallelic dosages not yet supported by LoadMultiallelicCenteredVarmaj()\n", stderr);
  exit(S_CAST(int32_t, kPglRetNotYetSupported));
  return kPglRetSuccess;
}

// This function handles error-message-logging.
PglErr LoadCenteredVarmajBlock(const uintptr_t* sample_include, PgrSampleSubsetIndex pssi, const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const double* allele_freqs, uint32_t variance_standardize, uint32_t is_haploid, uint32_t sample_ct, uint32_t variant_ct, PgenReader* simple_pgrp, double* normed_vmaj_iter, uintptr_t* variant_include_has_missing, uint32_t* cur_batch_sizep, uint32_t* variant_idxp, uintptr_t* variant_uidxp, uintptr_t* allele_idx_basep, uint32_t* cur_allele_ctp, uint32_t* incomplete_allele_idxp, PgenVariant* pgvp, double* allele_1copy_buf) {
  const uint32_t std_batch_size = *cur_batch_sizep;
  uint32_t variant_idx = *variant_idxp;
  uintptr_t variant_uidx = *variant_uidxp;
  uintptr_t allele_idx_base = *allele_idx_basep;
  uint32_t cur_allele_ct = *cur_allele_ctp;
  uint32_t incomplete_allele_idx = *incomplete_allele_idxp;
  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, variant_uidx + (incomplete_allele_idx != 0), &variant_uidx_base, &cur_bits);
  for (uint32_t allele_bidx = 0; allele_bidx != std_batch_size; ) {
    uint32_t missing_present = 0;
    if (!incomplete_allele_idx) {
      variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (!allele_idx_offsets) {
        allele_idx_base = variant_uidx;
      } else {
        allele_idx_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_base;
        allele_idx_base -= variant_uidx;
      }
    }
    uint32_t allele_idx_stop;
    uint32_t allele_idx_end;
    PglErr reterr;
    if (cur_allele_ct == 2) {
      allele_idx_stop = 1;
      allele_idx_end = 1;
      reterr = LoadBiallelicCenteredVarmaj(sample_include, pssi, variance_standardize, is_haploid, sample_ct, variant_uidx, allele_freqs[allele_idx_base], simple_pgrp, variant_include_has_missing? (&missing_present) : nullptr, normed_vmaj_iter, pgvp->genovec, pgvp->dosage_present, pgvp->dosage_main);
    } else {
      allele_idx_end = cur_allele_ct;
      allele_idx_stop = std_batch_size + incomplete_allele_idx - allele_bidx;
      if (allele_idx_stop > allele_idx_end) {
        allele_idx_stop = allele_idx_end;
      }
      reterr = LoadMultiallelicCenteredVarmaj(sample_include, pssi, &(allele_freqs[allele_idx_base]), variance_standardize, is_haploid, sample_ct, variant_uidx, cur_allele_ct, incomplete_allele_idx, allele_idx_stop, simple_pgrp, variant_include_has_missing? (&missing_present) : nullptr, normed_vmaj_iter, pgvp, allele_1copy_buf);
    }
    if (unlikely(reterr)) {
      if (reterr == kPglRetDegenerateData) {
        logputs("\n");
        logerrputs("Error: Zero-MAF variant is not actually monomorphic.  (This is possible when\ne.g. MAF is estimated from founders, but the minor allele was only observed in\nnonfounders.  In any case, you should be using e.g. --maf to filter out all\nvery-low-MAF variants, since the relationship matrix distance formula does not\nhandle them well.)\n");
      } else {
        PgenErrPrintNV(reterr, variant_uidx);
      }
      return reterr;
    }
    if (missing_present) {
      SetBit(variant_uidx, variant_include_has_missing);
    }
    const uintptr_t incr = allele_idx_stop - incomplete_allele_idx;
    normed_vmaj_iter = &(normed_vmaj_iter[incr * sample_ct]);
    allele_bidx += incr;
    if (allele_idx_stop == allele_idx_end) {
      if (++variant_idx == variant_ct) {
        *cur_batch_sizep = allele_bidx;
        break;
      }
      incomplete_allele_idx = 0;
    } else {
      incomplete_allele_idx = allele_idx_stop;
    }
  }
  *variant_idxp = variant_idx;
  *variant_uidxp = variant_uidx + (incomplete_allele_idx == 0);
  *allele_idx_basep = allele_idx_base;
  *cur_allele_ctp = cur_allele_ct;
  *incomplete_allele_idxp = incomplete_allele_idx;
  return kPglRetSuccess;
}

CONSTI32(kGrmVariantBlockSize, 144);

typedef struct CalcGrmPartCtxStruct {
  uint32_t* thread_start;
  uint32_t sample_ct;

  uint32_t cur_batch_size;
  double* normed_dosage_vmaj_bufs[2];
  double* normed_dosage_smaj_bufs[2];

  double* grm;
} CalcGrmPartCtx;

// turns out dsyrk_ does exactly what we want here
THREAD_FUNC_DECL CalcGrmThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  assert(!arg->tidx);
  CalcGrmPartCtx* ctx = S_CAST(CalcGrmPartCtx*, arg->sharedp->context);
  const uint32_t sample_ct = ctx->sample_ct;
  double* grm = ctx->grm;
  uint32_t parity = 0;
  do {
    const uint32_t cur_batch_size = ctx->cur_batch_size;
    if (cur_batch_size) {
      TransposeMultiplySelfIncr(ctx->normed_dosage_vmaj_bufs[parity], sample_ct, cur_batch_size, grm);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

// can't use dsyrk_, so we manually partition the GRM piece we need to compute
// into an appropriate number of sub-pieces
THREAD_FUNC_DECL CalcGrmPartThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcGrmPartCtx* ctx = S_CAST(CalcGrmPartCtx*, arg->sharedp->context);

  const uintptr_t sample_ct = ctx->sample_ct;
  const uintptr_t first_thread_row_start_idx = ctx->thread_start[0];
  const uintptr_t row_start_idx = ctx->thread_start[tidx];
  const uintptr_t sample_idx_end = ctx->thread_start[tidx + 1];
  const uintptr_t row_ct = sample_idx_end - row_start_idx;
  double* grm_piece = &(ctx->grm[(row_start_idx - first_thread_row_start_idx) * sample_ct]);
  uint32_t parity = 0;
  do {
    const uintptr_t cur_batch_size = ctx->cur_batch_size;
    if (cur_batch_size) {
      double* normed_vmaj = ctx->normed_dosage_vmaj_bufs[parity];
      double* normed_smaj = ctx->normed_dosage_smaj_bufs[parity];
      // quasi-bugfix (18 Mar 2022): forgot to skip extra right-side columns
      RowMajorMatrixMultiplyStridedIncr(&(normed_smaj[row_start_idx * cur_batch_size]), normed_vmaj, row_ct, cur_batch_size, sample_idx_end, sample_ct, cur_batch_size, sample_ct, grm_piece);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

CONSTI32(kDblMissingBlockWordCt, 128 / kBitsPerWord);
CONSTI32(kDblMissingBlockSize, 128);

typedef struct CalcDblMissingCtxStruct {
  uint32_t* thread_start;
  // missing_nz bit is set iff that sample has at least one missing entry in
  // current block
  uintptr_t* missing_nz[2];
  uintptr_t* missing_smaj[2];
  uint32_t* missing_dbl_exclude_cts;
} CalcDblMissingCtx;

THREAD_FUNC_DECL CalcDblMissingThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcDblMissingCtx* ctx = S_CAST(CalcDblMissingCtx*, arg->sharedp->context);

  const uint64_t first_thread_row_start_idx = ctx->thread_start[0];
  const uint64_t dbl_exclude_offset = (first_thread_row_start_idx * (first_thread_row_start_idx - 1)) / 2;
  const uint32_t row_start_idx = ctx->thread_start[tidx];
  const uintptr_t row_end_idx = ctx->thread_start[tidx + 1];
  uint32_t* missing_dbl_exclude_cts = ctx->missing_dbl_exclude_cts;
  uint32_t parity = 0;
  do {
    const uintptr_t* missing_nz = ctx->missing_nz[parity];
    const uintptr_t* missing_smaj = ctx->missing_smaj[parity];
    const uint32_t first_idx = AdvBoundedTo1Bit(missing_nz, 0, row_end_idx);
    uint32_t sample_idx = first_idx;
    uint32_t prev_missing_nz_ct = 0;
    if (sample_idx < row_start_idx) {
      sample_idx = AdvBoundedTo1Bit(missing_nz, row_start_idx, row_end_idx);
      if (sample_idx != row_end_idx) {
        prev_missing_nz_ct = PopcountBitRange(missing_nz, 0, row_start_idx);
      }
    }
    while (sample_idx < row_end_idx) {
      // todo: compare this explicit unroll with ordinary iteration over a
      // cur_words[] array
      // todo: try 1 word at a time, and 30 words at a time
      const uintptr_t cur_word0 = missing_smaj[sample_idx * kDblMissingBlockWordCt];
      const uintptr_t cur_word1 = missing_smaj[sample_idx * kDblMissingBlockWordCt + 1];
#ifndef __LP64__
      const uintptr_t cur_word2 = missing_smaj[sample_idx * kDblMissingBlockWordCt + 2];
      const uintptr_t cur_word3 = missing_smaj[sample_idx * kDblMissingBlockWordCt + 3];
#endif
      uintptr_t sample_idx2_base;
      uintptr_t cur_bits;
      BitIter1Start(missing_nz, first_idx, &sample_idx2_base, &cur_bits);
      // (sample_idx - 1) underflow ok
      uint32_t* write_base = &(missing_dbl_exclude_cts[((S_CAST(uint64_t, sample_idx) * (sample_idx - 1)) / 2) - dbl_exclude_offset]);
      for (uint32_t uii = 0; uii != prev_missing_nz_ct; ++uii) {
        const uint32_t sample_idx2 = BitIter1(missing_nz, &sample_idx2_base, &cur_bits);
        const uintptr_t* cur_missing_smaj_base = &(missing_smaj[sample_idx2 * kDblMissingBlockWordCt]);
        const uintptr_t cur_and0 = cur_word0 & cur_missing_smaj_base[0];
        const uintptr_t cur_and1 = cur_word1 & cur_missing_smaj_base[1];
#ifdef __LP64__
        if (cur_and0 || cur_and1) {
          write_base[sample_idx2] += Popcount2Words(cur_and0, cur_and1);
        }
#else
        const uintptr_t cur_and2 = cur_word2 & cur_missing_smaj_base[2];
        const uintptr_t cur_and3 = cur_word3 & cur_missing_smaj_base[3];
        if (cur_and0 || cur_and1 || cur_and2 || cur_and3) {
          write_base[sample_idx2] += Popcount4Words(cur_and0, cur_and1, cur_and2, cur_and3);
        }
#endif
      }
      ++prev_missing_nz_ct;
      sample_idx = AdvBoundedTo1Bit(missing_nz, sample_idx + 1, row_end_idx);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr CalcMissingMatrix(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* variant_include, uint32_t sample_ct, uint32_t variant_ct, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t row_start_idx, uintptr_t row_end_idx, uint32_t max_thread_ct, PgenReader* simple_pgrp, uint32_t** missing_cts_ptr, uint32_t** missing_dbl_exclude_cts_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  ThreadGroup tg;
  PreinitThreads(&tg);
  PglErr reterr = kPglRetSuccess;
  {
    const uintptr_t row_end_idxl = BitCtToWordCt(row_end_idx);
    // bugfix (1 Oct 2017): missing_vmaj rows must be vector-aligned
    const uintptr_t row_end_idxaw = BitCtToAlignedWordCt(row_end_idx);
    uintptr_t* missing_vmaj = nullptr;
    uintptr_t* genovec_buf = nullptr;
    CalcDblMissingCtx ctx;
    if (unlikely(bigstack_calloc_u32(row_end_idx, missing_cts_ptr) ||
                 bigstack_calloc_u32((S_CAST(uint64_t, row_end_idx) * (row_end_idx - 1) - S_CAST(uint64_t, row_start_idx) * (row_start_idx - 1)) / 2, missing_dbl_exclude_cts_ptr) ||
                 bigstack_calloc_w(row_end_idxl, &ctx.missing_nz[0]) ||
                 bigstack_calloc_w(row_end_idxl, &ctx.missing_nz[1]) ||
                 bigstack_alloc_w(NypCtToWordCt(row_end_idx), &genovec_buf) ||
                 bigstack_alloc_w(row_end_idxaw * (k1LU * kDblMissingBlockSize), &missing_vmaj) ||
                 bigstack_alloc_w(RoundUpPow2(row_end_idx, 2) * kDblMissingBlockWordCt, &ctx.missing_smaj[0]) ||
                 bigstack_alloc_w(RoundUpPow2(row_end_idx, 2) * kDblMissingBlockWordCt, &ctx.missing_smaj[1]))) {
      goto CalcMissingMatrix_ret_NOMEM;
    }
    uint32_t* missing_cts = *missing_cts_ptr;
    uint32_t* missing_dbl_exclude_cts = *missing_dbl_exclude_cts_ptr;
    ctx.missing_dbl_exclude_cts = missing_dbl_exclude_cts;
    VecW* transpose_bitblock_wkspace = S_CAST(VecW*, bigstack_alloc_raw(kPglBitTransposeBufbytes));
    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg) ||
                 bigstack_alloc_u32(calc_thread_ct + 1, &ctx.thread_start))) {
      goto CalcMissingMatrix_ret_NOMEM;
    }
    // note that this ctx.thread_start[] may have different values than the one
    // computed by CalcGrm(), since calc_thread_ct changes in the MTBLAS and
    // OS X cases.
    TriangleFill(sample_ct, calc_thread_ct, parallel_idx, parallel_tot, 0, 1, ctx.thread_start);
    assert(ctx.thread_start[0] == row_start_idx);
    assert(ctx.thread_start[calc_thread_ct] == row_end_idx);
    SetThreadFuncAndData(CalcDblMissingThread, &ctx, &tg);
    const uint32_t sample_transpose_batch_ct_m1 = (row_end_idx - 1) / kPglBitTransposeBatch;

    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    // caller's responsibility to print this
    // logputs("Correcting for missingness: ");
    fputs("0%", stdout);
    fflush(stdout);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    for (uint32_t cur_variant_idx_start = 0; ; ) {
      uint32_t cur_batch_size = 0;
      if (!IsLastBlock(&tg)) {
        cur_batch_size = kDblMissingBlockSize;
        uint32_t cur_variant_idx_end = cur_variant_idx_start + cur_batch_size;
        if (cur_variant_idx_end > variant_ct) {
          cur_batch_size = variant_ct - cur_variant_idx_start;
          cur_variant_idx_end = variant_ct;
          ZeroWArr((kDblMissingBlockSize - cur_batch_size) * row_end_idxaw, &(missing_vmaj[cur_batch_size * row_end_idxaw]));
        }
        uintptr_t* missing_vmaj_iter = missing_vmaj;
        for (uint32_t variant_idx = cur_variant_idx_start; variant_idx != cur_variant_idx_end; ++variant_idx) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          reterr = PgrGetMissingnessD(sample_include, pssi, row_end_idx, variant_uidx, simple_pgrp, nullptr, missing_vmaj_iter, nullptr, genovec_buf);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto CalcMissingMatrix_ret_1;
          }
          missing_vmaj_iter = &(missing_vmaj_iter[row_end_idxaw]);
        }
        uintptr_t* cur_missing_smaj_iter = ctx.missing_smaj[parity];
        uint32_t sample_batch_size = kPglBitTransposeBatch;
        for (uint32_t sample_transpose_batch_idx = 0; ; ++sample_transpose_batch_idx) {
          if (sample_transpose_batch_idx >= sample_transpose_batch_ct_m1) {
            if (sample_transpose_batch_idx > sample_transpose_batch_ct_m1) {
              break;
            }
            sample_batch_size = ModNz(row_end_idx, kPglBitTransposeBatch);
          }
          // missing_smaj offset needs to be 64-bit if kDblMissingBlockWordCt
          // increases
          TransposeBitblock(&(missing_vmaj[sample_transpose_batch_idx * kPglBitTransposeWords]), row_end_idxaw, kDblMissingBlockWordCt, kDblMissingBlockSize, sample_batch_size, &(cur_missing_smaj_iter[sample_transpose_batch_idx * kPglBitTransposeBatch * kDblMissingBlockWordCt]), transpose_bitblock_wkspace);
        }
        uintptr_t* cur_missing_nz = ctx.missing_nz[parity];
        ZeroWArr(row_end_idxl, cur_missing_nz);
        for (uint32_t sample_idx = 0; sample_idx != row_end_idx; ++sample_idx) {
          const uintptr_t cur_word0 = *cur_missing_smaj_iter++;
          const uintptr_t cur_word1 = *cur_missing_smaj_iter++;
#ifdef __LP64__
          if (cur_word0 || cur_word1) {
            SetBit(sample_idx, cur_missing_nz);
            missing_cts[sample_idx] += Popcount2Words(cur_word0, cur_word1);
          }
#else
          const uintptr_t cur_word2 = *cur_missing_smaj_iter++;
          const uintptr_t cur_word3 = *cur_missing_smaj_iter++;
          if (cur_word0 || cur_word1 || cur_word2 || cur_word3) {
            SetBit(sample_idx, cur_missing_nz);
            missing_cts[sample_idx] += Popcount4Words(cur_word0, cur_word1, cur_word2, cur_word3);
          }
#endif
        }
      }
      if (cur_variant_idx_start) {
        JoinThreads(&tg);
        // CalcDblMissingThread() never errors out
        if (IsLastBlock(&tg)) {
          break;
        }
        if (cur_variant_idx_start >= next_print_variant_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (cur_variant_idx_start * 100LLU) / variant_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
        }
      }
      if (cur_variant_idx_start + cur_batch_size == variant_ct) {
        DeclareLastThreadBlock(&tg);
      }
      if (unlikely(SpawnThreads(&tg))) {
        goto CalcMissingMatrix_ret_THREAD_CREATE_FAIL;
      }
      cur_variant_idx_start += cur_batch_size;
      parity = 1 - parity;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    bigstack_mark = R_CAST(unsigned char*, ctx.missing_nz[0]);
  }
  while (0) {
  CalcMissingMatrix_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CalcMissingMatrix_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 CalcMissingMatrix_ret_1:
  CleanupThreads(&tg);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr CalcGrm(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const double* allele_freqs, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, GrmFlags grm_flags, double grm_sparse_cutoff, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end, double** grm_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  CompressStreamState css;
  ThreadGroup tg;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  PreinitThreads(&tg);
  {
    assert(variant_ct);
#if defined(__APPLE__) || defined(USE_MTBLAS)
    uint32_t calc_thread_ct = 1;
#else
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct * parallel_tot > sample_ct / 32) {
      calc_thread_ct = sample_ct / (32 * parallel_tot);
      if (!calc_thread_ct) {
        calc_thread_ct = 1;
      }
    }
#endif
    if (unlikely(sample_ct < 2)) {
      logerrputs("Error: GRM construction requires at least two samples.\n");
      goto CalcGrm_ret_DEGENERATE_DATA;
    }
    const uintptr_t* sample_include = orig_sample_include;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uint32_t row_start_idx = 0;
    uintptr_t row_end_idx = sample_ct;
    uint32_t* thread_start = nullptr;
    if ((calc_thread_ct != 1) || (parallel_tot != 1)) {
      // note that grm should be allocated on bottom if no --parallel, since it
      // may continue to be used after function exit.  So we allocate this on
      // top.
      if (unlikely(bigstack_end_alloc_u32(calc_thread_ct + 1, &thread_start))) {
        goto CalcGrm_ret_NOMEM;
      }
      // slightly different from plink 1.9 since we don't bother to treat the
      // diagonal as a special case any more.
      TriangleFill2(sample_ct, calc_thread_ct, parallel_idx, parallel_tot, 0, thread_start);
      row_start_idx = thread_start[0];
      row_end_idx = thread_start[calc_thread_ct];
      if (row_end_idx < sample_ct) {
        // 0
        // 0 0
        // 0 0 0
        // 0 0 0 0
        // 1 1 1 1 1
        // 1 1 1 1 1 1
        // 2 2 2 2 2 2 2
        // 2 2 2 2 2 2 2 2
        // If we're computing part 0, we never need to load the last 4 samples;
        // if part 1, we don't need the last two; etc.
        uintptr_t* new_sample_include;
        if (unlikely(bigstack_alloc_w(raw_sample_ctl, &new_sample_include))) {
          goto CalcGrm_ret_NOMEM;
        }
        const uint32_t sample_uidx_end = 1 + IdxToUidxBasic(orig_sample_include, row_end_idx - 1);
        memcpy(new_sample_include, orig_sample_include, RoundUpPow2(sample_uidx_end, kBitsPerWord) / CHAR_BIT);
        ClearBitsNz(sample_uidx_end, raw_sample_ctl * kBitsPerWord, new_sample_include);
        sample_include = new_sample_include;
      }
      if ((!parallel_idx) && (calc_thread_ct == 1)) {
        thread_start = nullptr;
      }
    }

    CalcGrmPartCtx ctx;
    ctx.thread_start = thread_start;
    double* grm;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto CalcGrm_ret_NOMEM;
    }
    if (unlikely(bigstack_calloc64_d(S_CAST(uint64_t, row_end_idx - row_start_idx) * row_end_idx, &grm))) {
      if (!grm_ptr) {
        logerrputs("Error: Out of memory.  If you are SURE you are performing the right matrix\ncomputation, you can split it into smaller pieces with --parallel, and then\nconcatenate the results.  But before you try this, make sure the program you're\nproviding the matrix to can actually handle such a large input file.\n");
      } else {
        // Need to edit this if there are ever non-PCA ways to get here.
        if (!(grm_flags & kfGrmOutputMask)) {
          logerrputs("Error: Out of memory.  Consider \"--pca approx\" instead.\n");
        } else {
          logerrputs("Error: Out of memory.  Consider \"--pca approx\" (and not writing the GRM to\ndisk) instead.\n");
        }
      }
      goto CalcGrm_ret_NOMEM_CUSTOM;
    }
    ctx.sample_ct = row_end_idx;
    ctx.grm = grm;
    uint32_t* sample_include_cumulative_popcounts;
    PgenVariant pgv;
    double* allele_1copy_buf;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 BigstackAllocPgv(row_end_idx, allele_idx_offsets != nullptr, PgrGetGflags(simple_pgrp), &pgv) ||
                 bigstack_alloc_d(max_allele_ct, &allele_1copy_buf))) {
      goto CalcGrm_ret_NOMEM;
    }
    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    reterr = ConditionalAllocateNonAutosomalVariants(cip, "GRM construction", raw_variant_ct, &variant_include, &variant_ct);
    if (unlikely(reterr)) {
      goto CalcGrm_ret_1;
    }
    if (unlikely(bigstack_alloc_d(row_end_idx * kGrmVariantBlockSize, &ctx.normed_dosage_vmaj_bufs[0]) ||
                 bigstack_alloc_d(row_end_idx * kGrmVariantBlockSize, &ctx.normed_dosage_vmaj_bufs[1]))) {
      goto CalcGrm_ret_NOMEM;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* variant_include_has_missing = nullptr;
    if (!(grm_flags & kfGrmMeanimpute)) {
      if (unlikely(bigstack_calloc_w(raw_variant_ctl, &variant_include_has_missing))) {
        goto CalcGrm_ret_NOMEM;
      }
    }
    if (thread_start) {
      if (unlikely(bigstack_alloc_d(row_end_idx * kGrmVariantBlockSize, &ctx.normed_dosage_smaj_bufs[0]) ||
                   bigstack_alloc_d(row_end_idx * kGrmVariantBlockSize, &ctx.normed_dosage_smaj_bufs[1]))) {
        goto CalcGrm_ret_NOMEM;
      }
      SetThreadFuncAndData(CalcGrmPartThread, &ctx, &tg);
    } else {
      // defensive
      ctx.normed_dosage_smaj_bufs[0] = nullptr;
      ctx.normed_dosage_smaj_bufs[1] = nullptr;
      SetThreadFuncAndData(CalcGrmThread, &ctx, &tg);
    }
#ifdef USE_MTBLAS
    const uint32_t blas_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    BLAS_SET_NUM_THREADS(blas_thread_ct);
#endif
    // Main workflow:
    // 1. Set n=0, load batch 0
    //
    // 2. Spawn threads processing batch n
    // 3. Increment n by 1
    // 4. Load batch n unless eof
    // 5. Join threads
    // 6. Goto step 2 unless eof
    const uint32_t variance_standardize = !(grm_flags & kfGrmCov);
    const uint32_t is_haploid = cip->haploid_mask[0] & 1;
    uint32_t cur_batch_size = kGrmVariantBlockSize;
    uint32_t variant_idx_start = 0;
    uint32_t variant_idx = 0;
    uintptr_t variant_uidx = 0;
    uintptr_t allele_idx_base = 0;
    uint32_t cur_allele_ct = 2;
    uint32_t incomplete_allele_idx = 0;
    uint32_t parity = 0;
    uint32_t is_not_first_block = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    logputs("Constructing GRM: ");
    fputs("0%", stdout);
    fflush(stdout);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    while (1) {
      if (!IsLastBlock(&tg)) {
        double* normed_vmaj = ctx.normed_dosage_vmaj_bufs[parity];
        reterr = LoadCenteredVarmajBlock(sample_include, pssi, variant_include, allele_idx_offsets, allele_freqs, variance_standardize, is_haploid, row_end_idx, variant_ct, simple_pgrp, normed_vmaj, variant_include_has_missing, &cur_batch_size, &variant_idx, &variant_uidx, &allele_idx_base, &cur_allele_ct, &incomplete_allele_idx, &pgv, allele_1copy_buf);
        if (unlikely(reterr)) {
          goto CalcGrm_ret_1;
        }
        if (thread_start) {
          MatrixTransposeCopy(normed_vmaj, cur_batch_size, row_end_idx, ctx.normed_dosage_smaj_bufs[parity]);
        }
      }
      if (is_not_first_block) {
        JoinThreads(&tg);
        // CalcGrmPartThread() and CalcGrmThread() never error out
        if (IsLastBlock(&tg)) {
          break;
        }
        if (variant_idx_start >= next_print_variant_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (variant_idx_start * 100LLU) / variant_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
        }
      }
      ctx.cur_batch_size = cur_batch_size;
      if (variant_idx == variant_ct) {
        DeclareLastThreadBlock(&tg);
        cur_batch_size = 0;
      }
      if (unlikely(SpawnThreads(&tg))) {
        goto CalcGrm_ret_THREAD_CREATE_FAIL;
      }
      is_not_first_block = 1;
      variant_idx_start = variant_idx;
      parity = 1 - parity;
    }
    BLAS_SET_NUM_THREADS(1);
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    uint32_t* missing_cts = nullptr;  // stays null iff meanimpute
    uint32_t* missing_dbl_exclude_cts = nullptr;
    if (variant_include_has_missing) {
      const uint32_t variant_ct_with_missing = PopcountWords(variant_include_has_missing, raw_variant_ctl);
      // if no missing calls at all, act as if meanimpute was on
      if (variant_ct_with_missing) {
        logputs("Correcting for missingness... ");
        reterr = CalcMissingMatrix(sample_include, sample_include_cumulative_popcounts, variant_include_has_missing, sample_ct, variant_ct_with_missing, parallel_idx, parallel_tot, row_start_idx, row_end_idx, max_thread_ct, simple_pgrp, &missing_cts, &missing_dbl_exclude_cts);
        if (unlikely(reterr)) {
          goto CalcGrm_ret_1;
        }
      }
    }
    if (missing_cts) {
      // could parallelize this loop if it ever matters
      const uint32_t* missing_dbl_exclude_iter = missing_dbl_exclude_cts;
      for (uintptr_t row_idx = row_start_idx; row_idx != row_end_idx; ++row_idx) {
        const uint32_t variant_ct_base = variant_ct - missing_cts[row_idx];
        double* grm_iter = &(grm[(row_idx - row_start_idx) * row_end_idx]);
        for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
          *grm_iter++ /= u31tod(variant_ct_base - missing_cts[col_idx] + (*missing_dbl_exclude_iter++));
        }
        *grm_iter++ /= u31tod(variant_ct_base);
      }
    } else {
      const double variant_ct_recip = 1.0 / u31tod(variant_ct);
      for (uintptr_t row_idx = row_start_idx; row_idx != row_end_idx; ++row_idx) {
        double* grm_iter = &(grm[(row_idx - row_start_idx) * row_end_idx]);
        for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
          *grm_iter++ *= variant_ct_recip;
        }
      }
    }
    // N.B. Only the lower right of grm[] is valid when parallel_tot == 1.

    // possible todo: allow simultaneous --make-rel and
    // --make-grm-list/--make-grm-bin/--make-grm-sparse
    // (note that this routine may also be called by --pca, which may not write
    // a matrix to disk at all.)
    if (grm_flags & kfGrmOutputMask) {
      const GrmFlags matrix_shape = grm_flags & kfGrmMatrixShapemask;
      char* log_write_iter;
      if (matrix_shape) {
        // --make-rel
        fputs("--make-rel: Writing...", stdout);
        fflush(stdout);
        if (grm_flags & kfGrmMatrixBin) {
          char* outname_end2 = strcpya_k(outname_end, ".rel.bin");
          if (parallel_tot != 1) {
            *outname_end2++ = '.';
            outname_end2 = u32toa(parallel_idx + 1, outname_end2);
          }
          *outname_end2 = '\0';
          if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
            goto CalcGrm_ret_OPEN_FAIL;
          }
          double* write_double_buf = nullptr;
          if (matrix_shape == kfGrmMatrixSq0) {
            write_double_buf = R_CAST(double*, g_textbuf);
            ZeroDArr(kTextbufMainSize / sizeof(double), write_double_buf);
          } else if (matrix_shape == kfGrmMatrixSq) {
            if (unlikely(bigstack_alloc_d(row_end_idx - row_start_idx - 1, &write_double_buf))) {
              goto CalcGrm_ret_NOMEM;
            }
          }
          for (uintptr_t row_idx = row_start_idx; ; ) {
            const double* grm_row = &(grm[(row_idx - row_start_idx) * row_end_idx]);
            ++row_idx;
            if (unlikely(fwrite_checked(grm_row, row_idx * sizeof(double), outfile))) {
              goto CalcGrm_ret_WRITE_FAIL;
            }
            if (row_idx == row_end_idx) {
              break;
            }
            if (matrix_shape == kfGrmMatrixSq0) {
              uintptr_t zbytes_to_dump = (sample_ct - row_idx) * sizeof(double);
              while (zbytes_to_dump >= kTextbufMainSize) {
                if (unlikely(fwrite_checked(write_double_buf, kTextbufMainSize, outfile))) {
                  goto CalcGrm_ret_WRITE_FAIL;
                }
                zbytes_to_dump -= kTextbufMainSize;
              }
              if (zbytes_to_dump) {
                if (unlikely(fwrite_checked(write_double_buf, zbytes_to_dump, outfile))) {
                  goto CalcGrm_ret_WRITE_FAIL;
                }
              }
            } else if (matrix_shape == kfGrmMatrixSq) {
              double* write_double_iter = write_double_buf;
              const double* grm_col = &(grm[row_idx - 1]);
              for (uintptr_t row_idx2 = row_idx; row_idx2 != sample_ct; ++row_idx2) {
                *write_double_iter++ = grm_col[(row_idx2 - row_start_idx) * sample_ct];
              }
              if (unlikely(fwrite_checked(write_double_buf, (sample_ct - row_idx) * sizeof(double), outfile))) {
                goto CalcGrm_ret_WRITE_FAIL;
              }
            }
          }
          if (unlikely(fclose_null(&outfile))) {
            goto CalcGrm_ret_WRITE_FAIL;
          }
        } else if (grm_flags & kfGrmMatrixBin4) {
          // downcode all entries to floats
          char* outname_end2 = strcpya_k(outname_end, ".rel.bin");
          if (parallel_tot != 1) {
            *outname_end2++ = '.';
            outname_end2 = u32toa(parallel_idx + 1, outname_end2);
          }
          *outname_end2 = '\0';
          if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
            goto CalcGrm_ret_OPEN_FAIL;
          }
          float* write_float_buf;
          if (unlikely(bigstack_alloc_f(row_end_idx, &write_float_buf))) {
            goto CalcGrm_ret_NOMEM;
          }
          uintptr_t row_idx = row_start_idx;
          do {
            const double* grm_iter = &(grm[(row_idx - row_start_idx) * row_end_idx]);
            float* write_float_iter = write_float_buf;
            for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
              *write_float_iter++ = S_CAST(float, *grm_iter++);
            }
            ++row_idx;
            if (matrix_shape == kfGrmMatrixSq0) {
              ZeroFArr(sample_ct - row_idx, write_float_iter);
              write_float_iter = &(write_float_buf[sample_ct]);
            } else if (matrix_shape == kfGrmMatrixSq) {
              const double* grm_col = &(grm[row_idx - 1]);
              for (uintptr_t row_idx2 = row_idx; row_idx2 != sample_ct; ++row_idx2) {
                *write_float_iter++ = S_CAST(float, grm_col[(row_idx2 - row_start_idx) * sample_ct]);
              }
            }
            if (unlikely(fwrite_checked(write_float_buf, sizeof(float) * S_CAST(uintptr_t, write_float_iter - write_float_buf), outfile))) {
              goto CalcGrm_ret_WRITE_FAIL;
            }
          } while (row_idx < row_end_idx);
          if (unlikely(fclose_null(&outfile))) {
            goto CalcGrm_ret_WRITE_FAIL;
          }
        } else {
          char* outname_end2 = strcpya_k(outname_end, ".rel");
          if (parallel_tot != 1) {
            *outname_end2++ = '.';
            outname_end2 = u32toa(parallel_idx + 1, outname_end2);
          }
          const uint32_t output_zst = (grm_flags / kfGrmMatrixZs) & 1;
          if (output_zst) {
            outname_end2 = strcpya_k(outname_end2, ".zst");
          }
          *outname_end2 = '\0';
          reterr = InitCstreamAlloc(outname, 0, output_zst, max_thread_ct, kCompressStreamBlock + 16 * row_end_idx, &css, &cswritep);
          if (unlikely(reterr)) {
            goto CalcGrm_ret_1;
          }
          uintptr_t row_idx = row_start_idx;
          do {
            const double* grm_iter = &(grm[(row_idx - row_start_idx) * row_end_idx]);
            ++row_idx;
            for (uint32_t col_idx = 0; col_idx != row_idx; ++col_idx) {
              cswritep = dtoa_g(*grm_iter++, cswritep);
              *cswritep++ = '\t';
            }
            if (matrix_shape == kfGrmMatrixSq0) {
              AppendZerotabsUnsafe(sample_ct - row_idx, &cswritep);
            } else if (matrix_shape == kfGrmMatrixSq) {
              const double* grm_col = &(grm[row_idx - 1]);
              for (uintptr_t row_idx2 = row_idx; row_idx2 != sample_ct; ++row_idx2) {
                cswritep = dtoa_g(grm_col[(row_idx2 - row_start_idx) * sample_ct], cswritep);
                *cswritep++ = '\t';
              }
            }
            DecrAppendBinaryEoln(&cswritep);
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto CalcGrm_ret_WRITE_FAIL;
            }
          } while (row_idx < row_end_idx);
          if (unlikely(CswriteCloseNull(&css, cswritep))) {
            goto CalcGrm_ret_WRITE_FAIL;
          }
        }
        putc_unlocked('\r', stdout);
        log_write_iter = strcpya_k(g_logbuf, "--make-rel: GRM ");
        if (parallel_tot != 1) {
          log_write_iter = strcpya_k(log_write_iter, "component ");
        }
        log_write_iter = strcpya_k(log_write_iter, "written to ");
        log_write_iter = strcpya(log_write_iter, outname);
      } else {
        const uint32_t* missing_dbl_exclude_iter = missing_dbl_exclude_cts;
        if (grm_flags & kfGrmBin) {
          // --make-grm-bin
          float* write_float_buf;
          if (unlikely(bigstack_alloc_f(row_end_idx, &write_float_buf))) {
            goto CalcGrm_ret_NOMEM;
          }
          char* outname_end2 = strcpya_k(outname_end, ".grm.bin");
          if (parallel_tot != 1) {
            *outname_end2++ = '.';
            outname_end2 = u32toa(parallel_idx + 1, outname_end2);
          }
          *outname_end2 = '\0';
          if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
            goto CalcGrm_ret_OPEN_FAIL;
          }
          fputs("--make-grm-bin: Writing...", stdout);
          fflush(stdout);
          for (uintptr_t row_idx = row_start_idx; row_idx != row_end_idx; ++row_idx) {
            const double* grm_iter = &(grm[(row_idx - row_start_idx) * row_end_idx]);
            for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
              write_float_buf[col_idx] = S_CAST(float, *grm_iter++);
            }
            if (unlikely(fwrite_checked(write_float_buf, (row_idx + 1) * sizeof(float), outfile))) {
              goto CalcGrm_ret_WRITE_FAIL;
            }
          }
          if (unlikely(fclose_null(&outfile))) {
            goto CalcGrm_ret_WRITE_FAIL;
          }

          outname_end2 = strcpya_k(outname_end, ".grm.N.bin");
          if (parallel_tot != 1) {
            *outname_end2++ = '.';
            outname_end2 = u32toa(parallel_idx + 1, outname_end2);
          }
          *outname_end2 = '\0';
          if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
            goto CalcGrm_ret_OPEN_FAIL;
          }
          if (!missing_cts) {
            // trivial case: write the same number repeatedly
            const uintptr_t tot_cells = (S_CAST(uint64_t, row_end_idx) * (row_end_idx + 1) - S_CAST(uint64_t, row_start_idx) * (row_start_idx + 1)) / 2;
            const float variant_ctf = u31tof(variant_ct);
            write_float_buf = R_CAST(float*, g_textbuf);
            for (uint32_t uii = 0; uii != (kTextbufMainSize / sizeof(float)); ++uii) {
              write_float_buf[uii] = variant_ctf;
            }
            const uintptr_t full_write_ct = tot_cells / (kTextbufMainSize / sizeof(float));
            for (uintptr_t ulii = 0; ulii != full_write_ct; ++ulii) {
              if (unlikely(fwrite_checked(write_float_buf, kTextbufMainSize, outfile))) {
                goto CalcGrm_ret_WRITE_FAIL;
              }
            }
            const uintptr_t remainder = tot_cells % (kTextbufMainSize / sizeof(float));
            if (remainder) {
              if (unlikely(fwrite_checked(write_float_buf, remainder * sizeof(float), outfile))) {
                goto CalcGrm_ret_WRITE_FAIL;
              }
            }
          } else {
            for (uintptr_t row_idx = row_start_idx; row_idx != row_end_idx; ++row_idx) {
              const uint32_t variant_ct_base = variant_ct - missing_cts[row_idx];
              for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
                uint32_t cur_obs_ct = variant_ct_base;
                if (col_idx != row_idx) {
                  cur_obs_ct = cur_obs_ct - missing_cts[col_idx] + (*missing_dbl_exclude_iter++);
                }
                write_float_buf[col_idx] = u31tof(cur_obs_ct);
              }
              if (unlikely(fwrite_checked(write_float_buf, (row_idx + 1) * sizeof(float), outfile))) {
                goto CalcGrm_ret_WRITE_FAIL;
              }
            }
          }
          if (unlikely(fclose_null(&outfile))) {
            goto CalcGrm_ret_WRITE_FAIL;
          }
          putc_unlocked('\r', stdout);
          const uint32_t outname_copy_byte_ct = 5 + S_CAST(uintptr_t, outname_end - outname);
          log_write_iter = strcpya_k(g_logbuf, "--make-grm-bin: GRM ");
          if (parallel_tot != 1) {
            log_write_iter = strcpya_k(log_write_iter, "component ");
          }
          log_write_iter = strcpya_k(log_write_iter, "written to ");
          log_write_iter = memcpya(log_write_iter, outname, outname_copy_byte_ct);
          log_write_iter = strcpya_k(log_write_iter, "bin");
          if (parallel_tot != 1) {
            *log_write_iter++ = '.';
            log_write_iter = u32toa(parallel_idx + 1, log_write_iter);
          }
          log_write_iter = strcpya_k(log_write_iter, " , ");
          if (parallel_idx) {
            log_write_iter = strcpya_k(log_write_iter, "and ");
          }
          log_write_iter = strcpya_k(log_write_iter, "observation counts to ");
          log_write_iter = memcpya(log_write_iter, outname, outname_end2 - outname);
        } else {
          // --make-grm-list/--make-grm-sparse
          char* outname_end2 = strcpya_k(outname_end, ".grm");
          const uint32_t is_sparse = (grm_flags / kfGrmSparse) & 1;
          if (is_sparse) {
            outname_end2 = strcpya_k(outname_end2, ".sp");
          }
          if (parallel_tot != 1) {
            *outname_end2++ = '.';
            outname_end2 = u32toa(parallel_idx + 1, outname_end2);
          }
          const uint32_t is_zst = (grm_flags / kfGrmListZs) & 1;
          if (is_zst) {
            outname_end2 = strcpya_k(outname_end2, ".zst");
          }
          *outname_end2 = '\0';
          reterr = InitCstreamAlloc(outname, 0, is_zst, max_thread_ct, kCompressStreamBlock + kMaxMediumLine, &css, &cswritep);
          if (unlikely(reterr)) {
            goto CalcGrm_ret_1;
          }
          printf("--make-grm-%s: Writing...", is_sparse? "sparse" : "list");
          fflush(stdout);
          if (is_sparse) {
            for (uintptr_t row_idx = row_start_idx; row_idx != row_end_idx; ++row_idx) {
              const double* grm_iter = &(grm[(row_idx - row_start_idx) * row_end_idx]);
              for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
                const double rel_coeff = *grm_iter++;
                if (rel_coeff < grm_sparse_cutoff) {
                  continue;
                }
                cswritep = u32toa_x(row_idx, '\t', cswritep);
                cswritep = u32toa_x(col_idx, '\t', cswritep);
                // GCTA uses 8-digit precision here.
                cswritep = dtoa_g_p8(rel_coeff, cswritep);
                AppendBinaryEoln(&cswritep);
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto CalcGrm_ret_WRITE_FAIL;
                }
              }
            }
          } else {
            for (uintptr_t row_idx = row_start_idx; row_idx != row_end_idx; ++row_idx) {
              uint32_t variant_ct_base = variant_ct;
              if (missing_cts) {
                variant_ct_base -= missing_cts[row_idx];
              }
              const double* grm_iter = &(grm[(row_idx - row_start_idx) * row_end_idx]);
              for (uint32_t col_idx = 0; col_idx <= row_idx; ++col_idx) {
                cswritep = u32toa_x(row_idx + 1, '\t', cswritep);
                cswritep = u32toa_x(col_idx + 1, '\t', cswritep);
                if (missing_cts) {
                  uint32_t cur_obs_ct = variant_ct_base;
                  if (col_idx != row_idx) {
                    cur_obs_ct = cur_obs_ct - missing_cts[col_idx] + (*missing_dbl_exclude_iter++);
                  }
                  cswritep = u32toa(cur_obs_ct, cswritep);
                } else {
                  cswritep = u32toa(variant_ct_base, cswritep);
                }
                *cswritep++ = '\t';
                cswritep = dtoa_g(*grm_iter++, cswritep);
                AppendBinaryEoln(&cswritep);
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto CalcGrm_ret_WRITE_FAIL;
                }
              }
            }
          }
          if (unlikely(CswriteCloseNull(&css, cswritep))) {
            goto CalcGrm_ret_WRITE_FAIL;
          }
          putc_unlocked('\r', stdout);
          log_write_iter = strcpya_k(g_logbuf, "--make-grm-");
          if (is_sparse) {
            log_write_iter = strcpya_k(log_write_iter, "sparse");
          } else {
            log_write_iter = strcpya_k(log_write_iter, "list");
          }
          log_write_iter = strcpya_k(log_write_iter, ": GRM ");
          if (parallel_tot != 1) {
            log_write_iter = strcpya_k(log_write_iter, "component ");
          }
          log_write_iter = strcpya_k(log_write_iter, "written to ");
          log_write_iter = strcpya(log_write_iter, outname);
        }
      }
      if (!parallel_idx) {
        SampleIdFlags id_print_flags = siip->flags & kfSampleIdFidPresent;
        if (grm_flags & kfGrmNoIdHeader) {
          id_print_flags |= kfSampleIdNoIdHeader;
          if (grm_flags & kfGrmNoIdHeaderIidOnly) {
            id_print_flags |= kfSampleIdNoIdHeaderIidOnly;
          }
        }
        snprintf(&(outname_end[4]), kMaxOutfnameExtBlen - 4, ".id");
        reterr = WriteSampleIdsOverride(orig_sample_include, siip, outname, sample_ct, id_print_flags);
        if (unlikely(reterr)) {
          goto CalcGrm_ret_1;
        }
        log_write_iter = strcpya_k(log_write_iter, " , and IDs to ");
        log_write_iter = strcpya(log_write_iter, outname);
      }
      snprintf(log_write_iter, kLogbufSize - 2 * kPglFnamesize - 256, " .\n");
      WordWrapB(0);
      logputsb();
    }

    if (grm_ptr) {
      *grm_ptr = grm;
      // allocation right on top of grm[]
      bigstack_mark = R_CAST(unsigned char*, sample_include_cumulative_popcounts);
    }
  }
  while (0) {
  CalcGrm_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CalcGrm_ret_NOMEM_CUSTOM:
    reterr = kPglRetNomemCustomMsg;
    break;
  CalcGrm_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  CalcGrm_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  CalcGrm_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  CalcGrm_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 CalcGrm_ret_1:
  CswriteCloseCond(&css, cswritep);
  fclose_cond(outfile);
  CleanupThreads(&tg);
  BLAS_SET_NUM_THREADS(1);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

// should be able to remove NOLAPACK later since we already have a non-LAPACK
// SVD implementation
#ifndef NOLAPACK
// This seems to be better than 256 (due to avoidance of cache critical
// stride?)
// (still want this to be a multiple of 8, for cleaner multithreading)
// TODO: try to allow run to proceed with a smaller value when 240 induces an
// out-of-memory error.  It's cutting it close for UK Biobank.  Though we don't
// want to reduce the default value, I've tested that and it degrades
// performance.
CONSTI32(kPcaVariantBlockSize, 240);

typedef struct CalcPcaCtxStruct {
  uint32_t sample_ct;
  uint32_t pc_ct;

  double* yy_bufs[2];

  uint32_t cur_batch_size;

  double* g1;
  double* qq;
  double** y_transpose_bufs;
  double** g2_bb_part_bufs;
} CalcPcaCtx;

THREAD_FUNC_DECL CalcPcaXtxaThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcPcaCtx* ctx = S_CAST(CalcPcaCtx*, arg->sharedp->context);

  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t pc_ct_x2 = ctx->pc_ct * 2;
  const uintptr_t qq_col_ct = (ctx->pc_ct + 1) * pc_ct_x2;
  const uint32_t vidx_offset = tidx * kPcaVariantBlockSize;
  const double* g1 = ctx->g1;
  double* qq_iter = ctx->qq;
  double* y_transpose_buf = ctx->y_transpose_bufs[tidx];
  double* g2_part_buf = ctx->g2_bb_part_bufs[tidx];
  uint32_t parity = 0;
  do {
    const uint32_t cur_batch_size = ctx->cur_batch_size;
    if (vidx_offset < cur_batch_size) {
      uint32_t cur_thread_batch_size = cur_batch_size - vidx_offset;
      if (cur_thread_batch_size > kPcaVariantBlockSize) {
        cur_thread_batch_size = kPcaVariantBlockSize;
      }
      const double* yy_buf = &(ctx->yy_bufs[parity][S_CAST(uintptr_t, vidx_offset) * sample_ct]);
      double* cur_qq = &(qq_iter[vidx_offset * qq_col_ct]);
      RowMajorMatrixMultiplyStrided(yy_buf, g1, cur_thread_batch_size, sample_ct, pc_ct_x2, pc_ct_x2, sample_ct, qq_col_ct, cur_qq);
      MatrixTransposeCopy(yy_buf, cur_thread_batch_size, sample_ct, y_transpose_buf);
      RowMajorMatrixMultiplyStridedIncr(y_transpose_buf, cur_qq, sample_ct, cur_thread_batch_size, pc_ct_x2, qq_col_ct, cur_thread_batch_size, pc_ct_x2, g2_part_buf);
      qq_iter = &(qq_iter[cur_batch_size * qq_col_ct]);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

THREAD_FUNC_DECL CalcPcaXaThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcPcaCtx* ctx = S_CAST(CalcPcaCtx*, arg->sharedp->context);

  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t pc_ct_x2 = ctx->pc_ct * 2;
  const uintptr_t qq_col_ct = (ctx->pc_ct + 1) * pc_ct_x2;
  const uint32_t vidx_offset = tidx * kPcaVariantBlockSize;
  const double* g1 = ctx->g1;
  double* qq_iter = ctx->qq;
  uint32_t parity = 0;
  do {
    const uint32_t cur_batch_size = ctx->cur_batch_size;
    if (vidx_offset < cur_batch_size) {
      uint32_t cur_thread_batch_size = cur_batch_size - vidx_offset;
      if (cur_thread_batch_size > kPcaVariantBlockSize) {
        cur_thread_batch_size = kPcaVariantBlockSize;
      }
      const double* yy_buf = &(ctx->yy_bufs[parity][S_CAST(uintptr_t, vidx_offset) * sample_ct]);
      double* cur_qq = &(qq_iter[vidx_offset * qq_col_ct]);
      RowMajorMatrixMultiplyStrided(yy_buf, g1, cur_thread_batch_size, sample_ct, pc_ct_x2, pc_ct_x2, sample_ct, qq_col_ct, cur_qq);
      qq_iter = &(qq_iter[cur_batch_size * qq_col_ct]);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

THREAD_FUNC_DECL CalcPcaXtbThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcPcaCtx* ctx = S_CAST(CalcPcaCtx*, arg->sharedp->context);

  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t pc_ct_x2 = ctx->pc_ct * 2;
  const uintptr_t qq_col_ct = (ctx->pc_ct + 1) * pc_ct_x2;
  const uint32_t vidx_offset = tidx * kPcaVariantBlockSize;
  const double* qq_iter = &(ctx->qq[vidx_offset * qq_col_ct]);
  double* y_transpose_buf = ctx->y_transpose_bufs[tidx];
  double* bb_part_buf = ctx->g2_bb_part_bufs[tidx];
  uint32_t parity = 0;
  do {
    const uint32_t cur_batch_size = ctx->cur_batch_size;
    if (vidx_offset < cur_batch_size) {
      uint32_t cur_thread_batch_size = cur_batch_size - vidx_offset;
      if (cur_thread_batch_size > kPcaVariantBlockSize) {
        cur_thread_batch_size = kPcaVariantBlockSize;
      }
      const double* yy_buf = &(ctx->yy_bufs[parity][S_CAST(uintptr_t, vidx_offset) * sample_ct]);
      MatrixTransposeCopy(yy_buf, cur_thread_batch_size, sample_ct, y_transpose_buf);
      RowMajorMatrixMultiplyIncr(y_transpose_buf, qq_iter, sample_ct, qq_col_ct, cur_thread_batch_size, bb_part_buf);
      qq_iter = &(qq_iter[cur_batch_size * qq_col_ct]);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

typedef struct CalcPcaVarWtsCtxStruct {
  uint32_t sample_ct;
  uint32_t pc_ct;

  double* sample_wts_smaj;

  double* yy_bufs[2];

  uint32_t cur_batch_size;

  double* var_wts;
} CalcPcaVarWtsCtx;

THREAD_FUNC_DECL CalcPcaVarWtsThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcPcaVarWtsCtx* ctx = S_CAST(CalcPcaVarWtsCtx*, arg->sharedp->context);

  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t pc_ct = ctx->pc_ct;
  const uint32_t vidx_offset = tidx * kPcaVariantBlockSize;

  // either first batch size is calc_thread_ct * kPcaVariantBlockSize, or there
  // is only one batch
  const uintptr_t var_wts_part_size = S_CAST(uintptr_t, pc_ct) * ctx->cur_batch_size;

  const double* sample_wts = ctx->sample_wts_smaj;  // sample-major, pc_ct columns
  uint32_t parity = 0;
  do {
    const uint32_t cur_batch_size = ctx->cur_batch_size;
    if (vidx_offset < cur_batch_size) {
      uint32_t cur_thread_batch_size = cur_batch_size - vidx_offset;
      if (cur_thread_batch_size > kPcaVariantBlockSize) {
        cur_thread_batch_size = kPcaVariantBlockSize;
      }
      const double* yy_buf = &(ctx->yy_bufs[parity][S_CAST(uintptr_t, vidx_offset) * sample_ct]);
      // Variant weight matrix = X^T * S * D^{-1/2}, where X^T is the
      // variance-standardized genotype matrix, S is the sample weight matrix,
      // and D is a diagonal eigenvalue matrix.
      // We postpone the D^{-1/2} part for now, but it's straightforward to
      // switch to using precomputed (S * D^{-1/2}).
      double* cur_var_wts_part = &(ctx->var_wts[parity * var_wts_part_size + vidx_offset * S_CAST(uintptr_t, pc_ct)]);
      RowMajorMatrixMultiply(yy_buf, sample_wts, cur_thread_batch_size, pc_ct, sample_ct, cur_var_wts_part);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr FlushBiallelicVarWts(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* nonref_flags, const AlleleCode* maj_alleles, const double* var_wts_iter, const double* eigval_inv_sqrts, uint32_t batch_size, uint32_t pc_ct, PcaFlags pca_flags, uint32_t provref_col, uint32_t all_nonref, CompressStreamState* cssp, char** cswritepp, char* chr_buf, uint32_t* variant_idxp, uintptr_t* variant_uidxp, uint32_t* chr_fo_idxp, uint32_t* chr_endp, uint32_t* chr_buf_blenp) {
  char* cswritep = *cswritepp;
  uint32_t variant_idx = *variant_idxp;
  uintptr_t variant_uidx = *variant_uidxp;
  uint32_t chr_fo_idx = *chr_fo_idxp;
  uint32_t chr_end = *chr_endp;
  uint32_t chr_buf_blen = *chr_buf_blenp;

  const uint32_t variant_idx_stop = variant_idx + batch_size;
  const uint32_t ref_col = pca_flags & kfPcaVcolRef;
  const uint32_t alt1_col = pca_flags & kfPcaVcolAlt1;
  const uint32_t alt_col = pca_flags & kfPcaVcolAlt;
  const uint32_t maj_col = pca_flags & kfPcaVcolMaj;
  const uint32_t nonmaj_col = pca_flags & kfPcaVcolNonmaj;

  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, variant_uidx, &variant_uidx_base, &cur_bits);
  for (; variant_idx != variant_idx_stop; ++variant_idx) {
    variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    if (chr_buf) {
      // ok to skip this logic if chr_col not printed
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
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
      // guaranteed biallelic
      cswritep = strcpya(cswritep, cur_alleles[1]);
    }
    if (provref_col) {
      *cswritep++ = '\t';
      *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
    }
    const uint32_t maj_allele_idx = maj_alleles[variant_uidx];
    if (maj_col) {
      if (unlikely(Cswrite(cssp, &cswritep))) {
        return kPglRetWriteFail;
      }
      *cswritep++ = '\t';
      cswritep = strcpya(cswritep, cur_alleles[maj_allele_idx]);
    }
    if (nonmaj_col) {
      *cswritep++ = '\t';
      cswritep = strcpya(cswritep, cur_alleles[1 - maj_allele_idx]);
    }
    if (!maj_allele_idx) {
      for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
        *cswritep++ = '\t';
        // could avoid these multiplications by premultiplying the
        // sample weight matrix
        cswritep = dtoa_g((*var_wts_iter++) * eigval_inv_sqrts[pc_idx], cswritep);
      }
    } else {
      for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
        *cswritep++ = '\t';
        cswritep = dtoa_g((*var_wts_iter++) * (-eigval_inv_sqrts[pc_idx]), cswritep);
      }
    }
    AppendBinaryEoln(&cswritep);
    if (unlikely(Cswrite(cssp, &cswritep))) {
      // bugfix (15 Dec 2017): prevent buffer overflow when ALT, MAJ,
      // and NONMAJ columns all missing.
      return kPglRetWriteFail;
    }
  }
  *cswritepp = cswritep;
  *variant_idxp = variant_idx_stop;
  *variant_uidxp = variant_uidx + 1;
  *chr_fo_idxp = chr_fo_idx;
  *chr_endp = chr_end;
  *chr_buf_blenp = chr_buf_blen;
  return kPglRetSuccess;
}

PglErr FlushAlleleWts(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* nonref_flags, const double* var_wts_iter, const double* eigval_inv_sqrts, uint32_t batch_size, uint32_t pc_ct, PcaFlags pca_flags, uint32_t provref_col, uint32_t all_nonref, CompressStreamState* cssp, char** cswritepp, char* chr_buf, uint32_t* variant_idxp, uintptr_t* variant_uidxp, uintptr_t* allele_idx_offset_basep, uint32_t* cur_allele_ctp, uint32_t* incomplete_allele_idxp, uint32_t* chr_fo_idxp, uint32_t* chr_endp, uint32_t* chr_buf_blenp) {
  char* cswritep = *cswritepp;
  uint32_t variant_idx = *variant_idxp;
  uintptr_t variant_uidx = *variant_uidxp;
  uintptr_t allele_idx_offset_base = *allele_idx_offset_basep;
  uint32_t cur_allele_ct = *cur_allele_ctp;
  uint32_t incomplete_allele_idx = *incomplete_allele_idxp;
  uint32_t chr_fo_idx = *chr_fo_idxp;
  uint32_t chr_end = *chr_endp;
  uint32_t chr_buf_blen = *chr_buf_blenp;

  const uint32_t ref_col = pca_flags & kfPcaVcolRef;
  const uint32_t alt1_col = pca_flags & kfPcaVcolAlt1;
  const uint32_t alt_col = pca_flags & kfPcaVcolAlt;
  const uint32_t ax_col = pca_flags & kfPcaVcolAx;

  uintptr_t variant_uidx_base;
  uintptr_t cur_bits;
  BitIter1Start(variant_include, variant_uidx + (incomplete_allele_idx != 0), &variant_uidx_base, &cur_bits);
  for (uint32_t allele_bidx = 0; allele_bidx != batch_size; ) {
    if (!incomplete_allele_idx) {
      variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (chr_buf && (variant_uidx >= chr_end)) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
        *chr_name_end = '\t';
        chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
      }
      if (!allele_idx_offsets) {
        allele_idx_offset_base = variant_uidx * 2;
      } else {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
    }
    uint32_t allele_idx_end = cur_allele_ct;
    uint32_t allele_idx_stop;
    uint32_t incr;
    if (cur_allele_ct == 2) {
      allele_idx_stop = 2;
      incr = 1;
    } else {
      allele_idx_stop = batch_size + incomplete_allele_idx - allele_bidx;
      if (allele_idx_stop > allele_idx_end) {
        allele_idx_stop = allele_idx_end;
      }
      incr = allele_idx_stop - incomplete_allele_idx;
    }
    const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
    for (uint32_t allele_idx = incomplete_allele_idx; allele_idx != allele_idx_stop; ++allele_idx) {
      if (chr_buf) {
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
        for (uint32_t allele_idx2 = 1; allele_idx2 != cur_allele_ct; ++allele_idx2) {
          if (unlikely(Cswrite(cssp, &cswritep))) {
            return kPglRetWriteFail;
          }
          cswritep = strcpyax(cswritep, cur_alleles[allele_idx2], ',');
        }
        --cswritep;
      }
      if (provref_col) {
        *cswritep++ = '\t';
        *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
      }
      // A1 col always present
      if (unlikely(Cswrite(cssp, &cswritep))) {
        return kPglRetWriteFail;
      }
      *cswritep++ = '\t';
      cswritep = strcpya(cswritep, cur_alleles[allele_idx]);
      if (ax_col) {
        *cswritep++ = '\t';
        for (uint32_t allele_idx2 = 0; allele_idx2 != cur_allele_ct; ++allele_idx2) {
          if (allele_idx2 == allele_idx) {
            continue;
          }
          if (unlikely(Cswrite(cssp, &cswritep))) {
            return kPglRetWriteFail;
          }
          cswritep = strcpyax(cswritep, cur_alleles[allele_idx2], ',');
        }
        --cswritep;
      }
      if (cur_allele_ct == 2) {
        const double mult = allele_idx? -0.5 : 0.5;
        for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
          *cswritep++ = '\t';
          cswritep = dtoa_g((*var_wts_iter++) * mult * eigval_inv_sqrts[pc_idx], cswritep);
        }
        if (!allele_idx) {
          var_wts_iter -= pc_ct;
        }
      } else {
        for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
          *cswritep++ = '\t';
          cswritep = dtoa_g((*var_wts_iter++) * eigval_inv_sqrts[pc_idx], cswritep);
        }
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(cssp, &cswritep))) {
        return kPglRetWriteFail;
      }
    }
    allele_bidx += incr;
    if (allele_idx_stop == allele_idx_end) {
      ++variant_idx;
      incomplete_allele_idx = 0;
    } else {
      incomplete_allele_idx = allele_idx_stop;
    }
  }
  *cswritepp = cswritep;
  *variant_idxp = variant_idx;
  *variant_uidxp = variant_uidx + (incomplete_allele_idx == 0);
  *allele_idx_offset_basep = allele_idx_offset_base;
  *cur_allele_ctp = cur_allele_ct;
  *incomplete_allele_idxp = incomplete_allele_idx;
  *chr_fo_idxp = chr_fo_idx;
  *chr_endp = chr_end;
  *chr_buf_blenp = chr_buf_blen;
  return kPglRetSuccess;
}

PglErr CalcPca(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const double* allele_freqs, uint32_t raw_sample_ct, uintptr_t pca_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_allele_slen, uint32_t pc_ct, PcaFlags pca_flags, uint32_t max_thread_ct, PgenReader* simple_pgrp, sfmt_t* sfmtp, double* grm, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  CompressStreamState css;
  ThreadGroup tg;
  PreinitThreads(&tg);
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  {
    const uint32_t write_fid = FidColIsRequired(siip, pca_flags / kfPcaScolMaybefid);
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    const uint32_t write_sid = SidColIsRequired(sids, pca_flags / kfPcaScolMaybesid);
    const uint32_t is_approx = (pca_flags / kfPcaApprox) & 1;
    reterr = ConditionalAllocateNonAutosomalVariants(cip, is_approx? "PCA approximation" : "PCA", raw_variant_ct, &variant_include, &variant_ct);
    if (unlikely(reterr)) {
      goto CalcPca_ret_1;
    }
#ifdef __APPLE__
    // min OS X version is 10.7, so we can take Grand Central Dispatch dgemm
    // for granted
    // (tried this with Linux MKL + OpenMP as well, but results were inferior)
    uint32_t calc_thread_ct = 1;
#else
    // I/O thread generally has <1/8 of workload
    // TODO: recheck this, now that I/O thread is also responsible for fully
    // expanding dosages.  Still shouldn't be a big deal, but we probably want
    // sample_ct to affect the decision boundary now.
    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if ((calc_thread_ct - 1) * kPcaVariantBlockSize >= variant_ct) {
      calc_thread_ct = 1 + (variant_ct - 1) / kPcaVariantBlockSize;
    }
#endif
    if (unlikely(pc_ct > pca_sample_ct)) {
      // minor update (alpha 3): just error out here instead of trying to
      // auto-adjust PC count, number of .eigenvec output columns should be
      // easily predictable
      logerrprintf("Error: Too few samples to compute %u PCs with \"--pca approx\".\n", pc_ct);
      goto CalcPca_ret_DEGENERATE_DATA;
    }
    const uint32_t wts_requested = ((pca_flags & (kfPcaAlleleWts | kfPcaBiallelicVarWts)) != 0);
    const uint32_t biallelic_variant_ct = CountBiallelicVariants(variant_include, allele_idx_offsets, variant_ct);
    double* cur_var_wts = nullptr;
    double* eigval_inv_sqrts = nullptr;
    char* chr_buf = nullptr;
    uintptr_t overflow_buf_size = 3 * kMaxMediumLine;
    if (wts_requested) {
      if (pca_flags & kfPcaBiallelicVarWts) {
        if (unlikely(biallelic_variant_ct != variant_ct)) {
          logerrputs("Error: Multiallelic variant present in \"--pca biallelic-var-wts\" run.\n");
          goto CalcPca_ret_INCONSISTENT_INPUT;
        }
      }
      if (unlikely(bigstack_alloc_d(pc_ct, &cur_var_wts) ||
                   bigstack_alloc_d(pc_ct, &eigval_inv_sqrts))) {
        goto CalcPca_ret_NOMEM;
      }
      uint32_t max_chr_blen = 0;
      if (pca_flags & kfPcaVcolChrom) {
        max_chr_blen = GetMaxChrSlen(cip) + 1;
        if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
          goto CalcPca_ret_NOMEM;
        }
      }
      const uintptr_t overflow_buf_size2 = RoundUpPow2(kCompressStreamBlock + max_chr_blen + kMaxIdSlen + 2 * max_allele_slen + 32 + 16 * pc_ct, kCacheline);
      if (overflow_buf_size2 > overflow_buf_size) {
        overflow_buf_size = overflow_buf_size2;
      }
    }
    uintptr_t writebuf_alloc = overflow_buf_size;
    if (pca_flags & kfPcaVarZs) {
      writebuf_alloc += CstreamWkspaceReq(overflow_buf_size);
    }
    // temporary
    // todo: additional --pca-clusters allocations
    const uintptr_t* pca_sample_include = sample_include;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uint32_t* pca_sample_include_cumulative_popcounts;
    PgenVariant pgv;
    double* allele_1copy_buf = nullptr; // spurious g++ 4.8 warning
    double* eigvals;
    CalcPcaCtx ctx;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &pca_sample_include_cumulative_popcounts) ||
                 BigstackAllocPgv(pca_sample_ct, allele_idx_offsets != nullptr, PgrGetGflags(simple_pgrp), &pgv) ||
                 bigstack_alloc_d(max_allele_ct, &allele_1copy_buf) ||
                 bigstack_alloc_d(pc_ct, &eigvals) ||
                 SetThreadCt(calc_thread_ct, &tg))) {
      goto CalcPca_ret_NOMEM;
    }
    FillCumulativePopcounts(pca_sample_include, raw_sample_ctl, pca_sample_include_cumulative_popcounts);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(pca_sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    ctx.sample_ct = pca_sample_ct;
    ctx.pc_ct = pc_ct;
    const uintptr_t pca_row_ct = CountAlleles(variant_include, allele_idx_offsets, raw_variant_ct, variant_ct) - biallelic_variant_ct;
    const uint32_t is_haploid = cip->haploid_mask[0] & 1;
    uint32_t cur_allele_ct = 2;
    double* qq = nullptr;
    double* eigvecs_smaj;
    char* writebuf;
    if (is_approx) {
      if (pca_sample_ct <= 5000) {
        logerrputs("Warning: \"--pca approx\" is only recommended for analysis of >5000 samples.\n");
      }
      if (pca_row_ct > 5000000) {
        logerrputs("Warning: Use of \"--pca approx\" on >5m rows is not advisable.  Apply a MAF\nfilter if you haven't done so yet, and consider LD-pruning your variant set as\nwell.\n");
      }
      // This is ported from EIGENSOFT 6 src/ksrc/kjg_fpca.c , which is in turn
      // primarily based on Halko N, Martinsson P, Shkolnisky Y, Tygert M
      // (2011) An Algorithm for the Principal Component Analysis of Large Data
      // Sets.
      const uintptr_t pc_ct_x2 = pc_ct * 2;
      const uintptr_t qq_col_ct = (pc_ct + 1) * pc_ct_x2;
      // bugfix (30 Jan 2019): First SvdRectFused() call returns
      // min(variant_ct, qq_col_ct) singular vectors; this was previously
      // assumed to always be qq_col_ct, and very inaccurate results were
      // produced when the assumption wasn't true.
      // Simplest solution is to force the user to request fewer PCs, since the
      // final PCs wouldn't be accurate anyway.
      if (qq_col_ct > variant_ct) {
        logerrprintfww("Error: Too few variants to compute %u PCs with \"--pca approx\" (%u required).\n", pc_ct, qq_col_ct);
        goto CalcPca_ret_DEGENERATE_DATA;
      }
#ifndef LAPACK_ILP64
      if (unlikely((pca_row_ct * S_CAST(uint64_t, qq_col_ct)) > 0x7effffff)) {
        logerrputs("Error: \"--pca approx\" problem instance too large for this " PROG_NAME_STR " build.  If\nthis is really the computation you want, use a " PROG_NAME_STR " build with large-matrix\nsupport.\n");
        goto CalcPca_ret_INCONSISTENT_INPUT;
      }
#endif
      const double variant_ct_recip = 1.0 / u31tod(variant_ct);

      const uintptr_t gg_size = pca_sample_ct * pc_ct_x2;
      lapack_int svd_rect_lwork;
#ifdef LAPACK_ILP64
      GetSvdRectLwork(MAXV(pca_sample_ct, pca_row_ct), qq_col_ct, &svd_rect_lwork);
#else
      if (unlikely(GetSvdRectLwork(MAXV(pca_sample_ct, pca_row_ct), qq_col_ct, &svd_rect_lwork))) {
        logerrputs("Error: \"--pca approx\" problem instance too large for this " PROG_NAME_STR " build.  If\nthis is really the computation you want, use a " PROG_NAME_STR " build with large-matrix\nsupport.\n");
        goto CalcPca_ret_INCONSISTENT_INPUT;
      }
#endif
      uintptr_t svd_rect_wkspace_size = (svd_rect_lwork + qq_col_ct * qq_col_ct) * sizeof(double);
      if (svd_rect_wkspace_size < writebuf_alloc) {
        // used as writebuf later
        svd_rect_wkspace_size = writebuf_alloc;
      }

      unsigned char* svd_rect_wkspace;
      double* ss;
      double* g1;
      if (unlikely(bigstack_alloc_d(qq_col_ct, &ss) ||
                   bigstack_alloc_d(pca_row_ct * qq_col_ct, &qq) ||
                   bigstack_alloc_dp(calc_thread_ct, &ctx.y_transpose_bufs) ||
                   bigstack_alloc_dp(calc_thread_ct, &ctx.g2_bb_part_bufs) ||
                   bigstack_alloc_uc(svd_rect_wkspace_size, &svd_rect_wkspace) ||
                   bigstack_alloc_d(gg_size, &g1))) {
        goto CalcPca_ret_NOMEM;
      }
      const uintptr_t yy_alloc_incr = RoundUpPow2(kPcaVariantBlockSize * pca_sample_ct * sizeof(double), kCacheline);
      const uintptr_t b_size = pca_sample_ct * qq_col_ct;
      const uintptr_t g2_bb_part_alloc = RoundUpPow2(b_size * sizeof(double), kCacheline);
      // bugfix (16 Jan 2020)
      const uintptr_t per_thread_alloc = 3 * yy_alloc_incr + g2_bb_part_alloc;

      const uintptr_t bigstack_avail = bigstack_left();
      if (per_thread_alloc * calc_thread_ct > bigstack_avail) {
        if (unlikely(bigstack_avail < per_thread_alloc)) {
          goto CalcPca_ret_NOMEM;
        }
        calc_thread_ct = bigstack_avail / per_thread_alloc;
      }
      const uintptr_t yy_main_alloc = RoundUpPow2(kPcaVariantBlockSize * calc_thread_ct * pca_sample_ct * sizeof(double), kCacheline);
      ctx.yy_bufs[0] = S_CAST(double*, bigstack_alloc_raw(yy_main_alloc));
      ctx.yy_bufs[1] = S_CAST(double*, bigstack_alloc_raw(yy_main_alloc));
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        ctx.y_transpose_bufs[tidx] = S_CAST(double*, bigstack_alloc_raw(yy_alloc_incr));
        ctx.g2_bb_part_bufs[tidx] = S_CAST(double*, bigstack_alloc_raw(g2_bb_part_alloc));
      }
      FillGaussianDArr(gg_size / 2, max_thread_ct, sfmtp, g1);
      ctx.g1 = g1;
#ifdef __APPLE__
      fputs("Projecting random vectors... ", stdout);
#else
      printf("Projecting random vectors (%u compute thread%s)... ", calc_thread_ct, (calc_thread_ct == 1)? "" : "s");
#endif
      fflush(stdout);
      for (uint32_t iter_idx = 0; iter_idx <= pc_ct; ++iter_idx) {
        // kjg_fpca_XTXA(), kjg_fpca_XA()
        if (iter_idx < pc_ct) {
          SetThreadFuncAndData(CalcPcaXtxaThread, &ctx, &tg);
        } else {
          SetThreadFuncAndData(CalcPcaXaThread, &ctx, &tg);
        }
        for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
          ZeroDArr(gg_size, ctx.g2_bb_part_bufs[tidx]);
        }
        double* qq_iter = &(qq[iter_idx * pc_ct_x2]);  // offset on first row
        ctx.qq = qq_iter;

        // Main workflow:
        // 1. Set n=0, load batch 0
        //
        // 2. Spawn threads processing batch n
        // 3. Increment n by 1
        // 4. Load batch n unless eof
        // 5. Join threads
        // 6. Goto step 2 unless eof
        //
        // 7. Assemble next g1 by summing g2_parts
        uint32_t cur_batch_size = calc_thread_ct * kPcaVariantBlockSize;
        uint32_t variant_idx = 0;
        uintptr_t variant_uidx = 0;
        uintptr_t allele_idx_base = 0;
        uint32_t incomplete_allele_idx = 0;
        uint32_t parity = 0;
        uint32_t is_not_first_block = 0;
        while (1) {
          if (!IsLastBlock(&tg)) {
            reterr = LoadCenteredVarmajBlock(pca_sample_include, pssi, variant_include, allele_idx_offsets, allele_freqs, 1, is_haploid, pca_sample_ct, variant_ct, simple_pgrp, ctx.yy_bufs[parity], nullptr, &cur_batch_size, &variant_idx, &variant_uidx, &allele_idx_base, &cur_allele_ct, &incomplete_allele_idx, &pgv, allele_1copy_buf);
            if (unlikely(reterr)) {
              goto CalcPca_ret_1;
            }
          }
          if (is_not_first_block) {
            JoinThreads(&tg);
            if (IsLastBlock(&tg)) {
              break;
            }
          }
          ctx.cur_batch_size = cur_batch_size;
          if (variant_idx == variant_ct) {
            DeclareLastThreadBlock(&tg);
            cur_batch_size = 0;
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto CalcPca_ret_THREAD_CREATE_FAIL;
          }
          is_not_first_block = 1;
          parity = 1 - parity;
        }
        if (iter_idx < pc_ct) {
          memcpy(g1, ctx.g2_bb_part_bufs[0], gg_size * sizeof(double));
          for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
            const double* cur_g2_part = ctx.g2_bb_part_bufs[tidx];
            for (uintptr_t ulii = 0; ulii != gg_size; ++ulii) {
              g1[ulii] += cur_g2_part[ulii];
            }
          }
          for (uintptr_t ulii = 0; ulii != gg_size; ++ulii) {
            g1[ulii] *= variant_ct_recip;
          }
        }
#ifdef __APPLE__
        printf("\rProjecting random vectors... %u/%u", iter_idx + 1, pc_ct + 1);
#else
        printf("\rProjecting random vectors (%u compute thread%s)... %u/%u", calc_thread_ct, (calc_thread_ct == 1)? "" : "s", iter_idx + 1, pc_ct + 1);
#endif
        fflush(stdout);
      }
      fputs(".\n", stdout);
      logputs("Computing SVD of Krylov matrix... ");
      fflush(stdout);
      BLAS_SET_NUM_THREADS(max_thread_ct);
      IntErr svd_rect_err = SvdRectFused(pca_row_ct, qq_col_ct, svd_rect_lwork, qq, ss, svd_rect_wkspace);
      if (unlikely(svd_rect_err)) {
        logputs("\n");
        snprintf(g_logbuf, kLogbufSize, "Error: Failed to compute SVD of Krylov matrix (DGESVD info=%d).\n", S_CAST(int32_t, svd_rect_err));
        goto CalcPca_ret_DEGENERATE_DATA_2;
      }
      BLAS_SET_NUM_THREADS(1);
      logputs("done.\nRecovering top PCs from range approximation... ");
      fflush(stdout);

      // kjg_fpca_XTB()
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        ZeroDArr(b_size, ctx.g2_bb_part_bufs[tidx]);
      }
      SetThreadFuncAndData(CalcPcaXtbThread, &ctx, &tg);
      ctx.qq = qq;
      uint32_t cur_batch_size = calc_thread_ct * kPcaVariantBlockSize;
      uint32_t variant_idx = 0;
      uintptr_t variant_uidx = 0;
      uintptr_t allele_idx_base = 0;
      uint32_t incomplete_allele_idx = 0;
      uint32_t parity = 0;
      uint32_t is_not_first_block = 0;
      while (1) {
        if (!IsLastBlock(&tg)) {
          reterr = LoadCenteredVarmajBlock(pca_sample_include, pssi, variant_include, allele_idx_offsets, allele_freqs, 1, is_haploid, pca_sample_ct, variant_ct, simple_pgrp, ctx.yy_bufs[parity], nullptr, &cur_batch_size, &variant_idx, &variant_uidx, &allele_idx_base, &cur_allele_ct, &incomplete_allele_idx, &pgv, allele_1copy_buf);
          if (unlikely(reterr)) {
            // This error *didn't* happen on an earlier pass, so assign blame
            // to I/O.  (This may be additive with an error message printed by
            // LoadCenteredVarmajBlock().)
            goto CalcPca_ret_REWIND_FAIL;
          }
        }
        if (is_not_first_block) {
          JoinThreads(&tg);
          if (IsLastBlock(&tg)) {
            break;
          }
        }
        ctx.cur_batch_size = cur_batch_size;
        if (variant_idx == variant_ct) {
          DeclareLastThreadBlock(&tg);
          cur_batch_size = 0;
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto CalcPca_ret_THREAD_CREATE_FAIL;
        }
        is_not_first_block = 1;
        parity = 1 - parity;
      }
      double* bb = ctx.g2_bb_part_bufs[0];
      for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
        const double* cur_bb_part = ctx.g2_bb_part_bufs[tidx];
        for (uintptr_t ulii = 0; ulii != b_size; ++ulii) {
          bb[ulii] += cur_bb_part[ulii];
        }
      }
      BLAS_SET_NUM_THREADS(max_thread_ct);
      svd_rect_err = SvdRectFused(pca_sample_ct, qq_col_ct, svd_rect_lwork, bb, ss, svd_rect_wkspace);
      if (unlikely(svd_rect_err)) {
        logputs("\n");
        snprintf(g_logbuf, kLogbufSize, "Error: Failed to compute SVD of final matrix (DGESVD info=%d).\n", S_CAST(int32_t, svd_rect_err));
        goto CalcPca_ret_DEGENERATE_DATA_2;
      }
      BLAS_SET_NUM_THREADS(1);
      logputs("done.\n");
      eigvecs_smaj = g1;
      for (uint32_t sample_idx = 0; sample_idx != pca_sample_ct; ++sample_idx) {
        memcpy(&(eigvecs_smaj[sample_idx * S_CAST(uintptr_t, pc_ct)]), &(bb[sample_idx * qq_col_ct]), pc_ct * sizeof(double));
      }
      for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
        eigvals[pc_idx] = ss[pc_idx] * ss[pc_idx] * variant_ct_recip;
      }
      writebuf = R_CAST(char*, svd_rect_wkspace);
      // bugfix (25 Jun 2018): eigvals[] computation was missing a divide-by-2
      // somewhere, in both diploid and haploid cases.
      // update (30 Jan 2019): er, actually, no.
      if (is_haploid) {
        for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
          eigvals[pc_idx] *= 0.5;
        }
      }
    } else {
#ifndef LAPACK_ILP64
      if (unlikely(pca_sample_ct > 46340)) {
        logerrputs("Error: --pca non-approximate problem instance too large for this " PROG_NAME_STR " build.\nTry \"--pca approx\" instead, or a " PROG_NAME_STR " build with large-matrix support.\n");
        goto CalcPca_ret_INCONSISTENT_INPUT;
      }
#endif
      {
        // If there are NaNs in grm, ExtractEigvecs should fail; but we want to
        // print a customized error message in that case, and grm is destroyed
        // by ExtractEigvecs().  Okay, checking it for NaNs after the fact
        // probably works in practice, but I'd rather not risk it; checking for
        // NaNs in advance is cheap enough.
        uintptr_t row_idx;
        uintptr_t col_idx;
        if (LowerTriangularFirstInfOrNan(grm, pca_sample_ct, &row_idx, &col_idx)) {
          const double invalid_val = grm[row_idx * pca_sample_ct + col_idx];
          char* id_pair_str = g_textbuf;
          {
            char* write_iter = id_pair_str;
            if (row_idx == col_idx) {
              write_iter = strcpya_k(write_iter, " '");
            } else {
              write_iter = strcpya_k(write_iter, "s '");
            }
            const uintptr_t sample_uidx1 = IdxToUidxBasic(sample_include, col_idx);
            write_iter = AppendSpacedXid(sample_ids, sids, write_fid, write_sid, max_sample_id_blen, max_sid_blen, sample_uidx1, write_iter);
            if (row_idx == col_idx) {
              strcpy_k(write_iter, "' and itself");
            } else {
              write_iter = strcpya_k(write_iter, "' and '");
              const uintptr_t sample_uidx2 = IdxToUidxBasic(sample_include, row_idx);
              write_iter = AppendSpacedXid(sample_ids, sids, write_fid, write_sid, max_sample_id_blen, max_sid_blen, sample_uidx2, write_iter);
              strcpy_k(write_iter, "'");
            }
          }
          if (invalid_val != invalid_val) {
            logerrprintfww("Error: NaN value present in GRM, between sample%s. This usually occurs when you have samples with an extremely high amount of missing data; use e.g. --mind to filter them out before retrying.\n", id_pair_str);
            goto CalcPca_ret_DEGENERATE_DATA;
          }
          // I'll call this an internal error, even though in principle this
          // could result from malicious --read-freq input.
          logerrprintfww("Error: Infinite value present in GRM, between sample%s. This should never happen; you have probably encountered a plink2 bug. If you report the error on GitHub or the plink2-users Google group (make sure to include the full .log file in your report), we'll try to address it.\n", id_pair_str);
          reterr = kPglRetInternalError;
          goto CalcPca_ret_1;
        }
      }
      lapack_int lwork;
      lapack_int liwork;
      uintptr_t wkspace_byte_ct;
      if (unlikely(GetExtractEigvecsLworks(pca_sample_ct, pc_ct, &lwork, &liwork, &wkspace_byte_ct))) {
        goto CalcPca_ret_NOMEM;
      }
      const uintptr_t eigvecs_smaj_alloc = pc_ct * pca_sample_ct * sizeof(double);
      if (wkspace_byte_ct < eigvecs_smaj_alloc) {
        wkspace_byte_ct = eigvecs_smaj_alloc;
      }
      double* reverse_eigvecs_pcmaj;
      unsigned char* extract_eigvecs_wkspace;
      if (unlikely(bigstack_alloc_d(pc_ct * pca_sample_ct, &reverse_eigvecs_pcmaj) ||
                   bigstack_alloc_uc(wkspace_byte_ct, &extract_eigvecs_wkspace))) {
        goto CalcPca_ret_NOMEM;
      }
      logprintf("Extracting eigenvalue%s and eigenvector%s... ", (pc_ct == 1)? "" : "s", (pc_ct == 1)? "" : "s");
      fflush(stdout);
      BLAS_SET_NUM_THREADS(max_thread_ct);
      // not putting unlikely() here for now.
      if (ExtractEigvecs(pca_sample_ct, pc_ct, lwork, liwork, grm, eigvals, reverse_eigvecs_pcmaj, extract_eigvecs_wkspace)) {
        logputs("\n");
        logerrputs("Error: Failed to extract eigenvector(s) from GRM.\n");
        goto CalcPca_ret_DEGENERATE_DATA;
      }
      BLAS_SET_NUM_THREADS(1);
      logputs("done.\n");
      eigvecs_smaj = R_CAST(double*, extract_eigvecs_wkspace);
      BigstackShrinkTop(eigvecs_smaj, eigvecs_smaj_alloc);
      if (unlikely(bigstack_alloc_c(writebuf_alloc, &writebuf))) {
        goto CalcPca_ret_NOMEM;
      }

      // ExtractEigvecs() results are in reverse order, and we also need to
      // transpose eigenvectors to sample-major
      const uint32_t pc_ct_m1 = pc_ct - 1;
      const uint32_t pc_ct_div2 = pc_ct / 2;
      for (uint32_t pc_idx = 0; pc_idx != pc_ct_div2; ++pc_idx) {
        double tmp_eigval = eigvals[pc_idx];
        eigvals[pc_idx] = eigvals[pc_ct_m1 - pc_idx];
        eigvals[pc_ct_m1 - pc_idx] = tmp_eigval;
      }
      double* eigvecs_smaj_iter = eigvecs_smaj;
      for (uint32_t sample_idx = 0; sample_idx != pca_sample_ct; ++sample_idx) {
        uintptr_t pc_inv_idx = pc_ct;
        const double* reverse_eigvecs_col = &(reverse_eigvecs_pcmaj[sample_idx]);
        do {
          --pc_inv_idx;
          *eigvecs_smaj_iter++ = reverse_eigvecs_col[pc_inv_idx * pca_sample_ct];
        } while (pc_inv_idx);
      }
    }
    // (later: --pca-cluster-names, --pca-clusters)
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);

    if (wts_requested) {
      CalcPcaVarWtsCtx vwctx;
      vwctx.sample_ct = pca_sample_ct;
      vwctx.pc_ct = pc_ct;
      vwctx.sample_wts_smaj = eigvecs_smaj;
      for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
        eigval_inv_sqrts[pc_idx] = 1.0 / sqrt(eigvals[pc_idx]);
      }

      const uint32_t allele_wts = (pca_flags / kfPcaAlleleWts) & 1;
      const uint32_t output_zst = (pca_flags / kfPcaVarZs) & 1;
      if (allele_wts) {
        OutnameZstSet(".eigenvec.allele", output_zst, outname_end);
      } else {
        OutnameZstSet(".eigenvec.var", output_zst, outname_end);
      }
      reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, writebuf, R_CAST(unsigned char*, &(writebuf[overflow_buf_size])), &css);
      if (unlikely(reterr)) {
        goto CalcPca_ret_1;
      }
      cswritep = writebuf;
      *cswritep++ = '#';
      if (chr_buf) {
        cswritep = strcpya_k(cswritep, "CHROM\t");
      }
      if (pca_flags & kfPcaVcolPos) {
        cswritep = strcpya_k(cswritep, "POS\t");
      } else {
        variant_bps = nullptr;
      }
      cswritep = strcpya_k(cswritep, "ID");
      if (pca_flags & kfPcaVcolRef) {
        cswritep = strcpya_k(cswritep, "\tREF");
      }
      if (pca_flags & kfPcaVcolAlt1) {
        cswritep = strcpya_k(cswritep, "\tALT1");
      }
      if (pca_flags & kfPcaVcolAlt) {
        cswritep = strcpya_k(cswritep, "\tALT");
      }
      const uintptr_t* nonref_flags = PgrGetNonrefFlags(simple_pgrp);
      const uint32_t all_nonref = (PgrGetGflags(simple_pgrp) & kfPgenGlobalAllNonref) && (!nonref_flags);
      const uint32_t provref_col = (pca_flags & kfPcaVcolRef) && ProvrefCol(variant_include, nonref_flags, pca_flags / kfPcaVcolMaybeprovref, raw_variant_ct, all_nonref);
      if (provref_col) {
        cswritep = strcpya_k(cswritep, "\tPROVISIONAL_REF?");
      }
      if (allele_wts) {
        cswritep = strcpya_k(cswritep, "\tA1");
      }
      if (pca_flags & kfPcaVcolAx) {
        cswritep = strcpya_k(cswritep, "\tAX");
      }
      if (pca_flags & kfPcaVcolMaj) {
        cswritep = strcpya_k(cswritep, "\tMAJ");
      }
      if (pca_flags & kfPcaVcolNonmaj) {
        cswritep = strcpya_k(cswritep, "\tNONMAJ");
      }
      for (uint32_t pc_idx = 1; pc_idx <= pc_ct; ++pc_idx) {
        cswritep = strcpya_k(cswritep, "\tPC");
        cswritep = u32toa(pc_idx, cswritep);
      }
      AppendBinaryEoln(&cswritep);

      // Main workflow:
      // 1. Set n=0, load batch 0
      //
      // 2. Spawn threads processing batch n
      // 3. If n>0, write results and update projection for block (n-1)
      // 4. Increment n by 1
      // 5. Load batch n unless eof
      // 6. Join threads
      // 7. Goto step 2 unless eof
      //
      // 8. Write results and update projection for last block
#ifndef __APPLE__
      if (output_zst) {
        // compression is relatively expensive?
        calc_thread_ct = 1;
      }
#endif
      uintptr_t var_wts_part_size;
      double* var_wts = qq;
      if (var_wts) {
        var_wts_part_size = (MINV(pca_row_ct, calc_thread_ct * kPcaVariantBlockSize)) * S_CAST(uintptr_t, pc_ct);
        vwctx.yy_bufs[0] = ctx.yy_bufs[0];
        vwctx.yy_bufs[1] = ctx.yy_bufs[1];
        vwctx.var_wts = ctx.qq;
      } else {
        // non-approximate PCA, some buffers have not been allocated yet

        // if grm[] (which we no longer need) has at least as much remaining
        // space as bigstack, allocate from grm
        unsigned char* arena_bottom = R_CAST(unsigned char*, grm);
        unsigned char* arena_top = bigstack_mark;
        uintptr_t arena_avail = arena_top - arena_bottom;
        if (arena_avail < bigstack_left()) {
          arena_bottom = g_bigstack_base;
          arena_top = g_bigstack_end;
          arena_avail = bigstack_left();
        }
        const uintptr_t var_wts_part_alloc = RoundUpPow2(2 * kPcaVariantBlockSize * sizeof(double) * pc_ct, kCacheline);
        const uintptr_t yy_alloc_incr = RoundUpPow2(kPcaVariantBlockSize * pca_sample_ct * sizeof(double), kCacheline);
        const uintptr_t per_thread_alloc = 2 * yy_alloc_incr + var_wts_part_alloc;
        if (per_thread_alloc * calc_thread_ct > arena_avail) {
          if (unlikely(arena_avail < per_thread_alloc)) {
            goto CalcPca_ret_NOMEM;
          }
          calc_thread_ct = arena_avail / per_thread_alloc;
        }
        const uintptr_t yy_main_alloc = RoundUpPow2(kPcaVariantBlockSize * calc_thread_ct * pca_sample_ct * sizeof(double), kCacheline);
        vwctx.yy_bufs[0] = S_CAST(double*, arena_alloc_raw(yy_main_alloc, &arena_bottom));
        vwctx.yy_bufs[1] = S_CAST(double*, arena_alloc_raw(yy_main_alloc, &arena_bottom));
        var_wts_part_size = (MINV(pca_row_ct, calc_thread_ct * kPcaVariantBlockSize)) * S_CAST(uintptr_t, pc_ct);
        var_wts = S_CAST(double*, arena_alloc_raw_rd(2 * var_wts_part_size * sizeof(double), &arena_bottom));
        vwctx.var_wts = var_wts;
#ifndef NDEBUG
        if (arena_top == g_bigstack_end) {
          // we shouldn't make any more allocations, but just in case...
          g_bigstack_base = arena_bottom;
        }
#endif
      }
      if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
        goto CalcPca_ret_NOMEM;
      }
      SetThreadFuncAndData(CalcPcaVarWtsThread, &vwctx, &tg);
      uint32_t prev_batch_size = 0;
      uint32_t cur_batch_size = calc_thread_ct * kPcaVariantBlockSize;

      uint32_t variant_idx_load = 0;
      uintptr_t variant_uidx_load = 0;
      uintptr_t allele_idx_base_load = 0;
      uint32_t cur_allele_ct_load = 2;
      uint32_t incomplete_allele_idx_load = 0;

      uint32_t variant_idx_write = 0;
      uintptr_t variant_uidx_write = 0;
      uintptr_t allele_idx_offset_write = 0;
      // cur_allele_ct = 2;
      uint32_t incomplete_allele_idx_write = 0;
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;

      uint32_t parity = 0;
      uint32_t is_not_first_block = 0;
      while (1) {
        if (!IsLastBlock(&tg)) {
          reterr = LoadCenteredVarmajBlock(pca_sample_include, pssi, variant_include, allele_idx_offsets, allele_freqs, 1, is_haploid, pca_sample_ct, variant_ct, simple_pgrp, vwctx.yy_bufs[parity], nullptr, &cur_batch_size, &variant_idx_load, &variant_uidx_load, &allele_idx_base_load, &cur_allele_ct_load, &incomplete_allele_idx_load, &pgv, allele_1copy_buf);
          if (unlikely(reterr)) {
            goto CalcPca_ret_1;
          }
        }
        if (is_not_first_block) {
          JoinThreads(&tg);
        }
        if (!IsLastBlock(&tg)) {
          vwctx.cur_batch_size = cur_batch_size;
          if (variant_idx_load == variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto CalcPca_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (is_not_first_block) {
          // write *previous* block results
          const double* var_wts_iter = &(var_wts[parity * var_wts_part_size]);
          // (todo: update projection here)
          if (allele_wts) {
            reterr = FlushAlleleWts(variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, nonref_flags, var_wts_iter, eigval_inv_sqrts, prev_batch_size, pc_ct, pca_flags, provref_col, all_nonref, &css, &cswritep, chr_buf, &variant_idx_write, &variant_uidx_write, &allele_idx_offset_write, &cur_allele_ct, &incomplete_allele_idx_write, &chr_fo_idx, &chr_end, &chr_buf_blen);
          } else {
            reterr = FlushBiallelicVarWts(variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, nonref_flags, maj_alleles, var_wts_iter, eigval_inv_sqrts, prev_batch_size, pc_ct, pca_flags, provref_col, all_nonref, &css, &cswritep, chr_buf, &variant_idx_write, &variant_uidx_write, &chr_fo_idx, &chr_end, &chr_buf_blen);
          }
          if (unlikely(reterr)) {
            // only write_fail possible in practice
            goto CalcPca_ret_1;
          }
          if (variant_idx_write == variant_ct) {
            break;
          }
        }
        is_not_first_block = 1;
        prev_batch_size = cur_batch_size;
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto CalcPca_ret_WRITE_FAIL;
      }
      logprintfww("--pca%s: %s weights written to %s .\n", is_approx? " approx" : "", allele_wts? "Allele" : "Variant", outname);
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".eigenvec");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto CalcPca_ret_OPEN_FAIL;
    }
    char* write_iter = writebuf;
    *write_iter++ = '#';
    if (write_fid) {
      write_iter = strcpya_k(write_iter, "FID\t");
    }
    write_iter = strcpya_k(write_iter, "IID");
    if (write_sid) {
      write_iter = strcpya_k(write_iter, "\tSID");
    }
    for (uint32_t pc_idx = 1; pc_idx <= pc_ct; ++pc_idx) {
      write_iter = strcpya_k(write_iter, "\tPC");
      write_iter = u32toa(pc_idx, write_iter);
    }
    AppendBinaryEoln(&write_iter);
    const uint32_t sample_ct = pca_sample_ct;
    uintptr_t sample_uidx_base = 0;
    uintptr_t sample_include_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
      write_iter = AppendXid(sample_ids, sids, write_fid, write_sid, max_sample_id_blen, max_sid_blen, sample_uidx, write_iter);
      double* sample_wts_iter = &(eigvecs_smaj[sample_idx * pc_ct]);
      // todo: read from proj_sample_wts instead when pca_sample_include bit
      // not set
      for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
        *write_iter++ = '\t';
        write_iter = dtoa_g(*sample_wts_iter++, write_iter);
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto CalcPca_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto CalcPca_ret_WRITE_FAIL;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".eigenval");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto CalcPca_ret_OPEN_FAIL;
    }
    write_iter = writebuf;
    for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
      write_iter = dtoa_g(eigvals[pc_idx], write_iter);
      AppendBinaryEoln(&write_iter);
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto CalcPca_ret_WRITE_FAIL;
    }
    *outname_end = '\0';
    logprintfww("--pca%s: Eigenvector%s written to %s.eigenvec , and eigenvalue%s written to %s.eigenval .\n", is_approx? " approx" : "", (pc_ct == 1)? "" : "s", outname, (pc_ct == 1)? "" : "s", outname);
  }
  while (0) {
  CalcPca_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CalcPca_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  CalcPca_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, ".pgen file");
    break;
  CalcPca_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  CalcPca_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  CalcPca_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  CalcPca_ret_DEGENERATE_DATA_2:
    logerrputsb();
  CalcPca_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 CalcPca_ret_1:
  CleanupThreads(&tg);
  BLAS_SET_NUM_THREADS(1);
  CswriteCloseCond(&css, cswritep);
  fclose_cond(outfile);
  if (grm) {
    // nothing after --pca in the plink2 order of operations uses grm[]
    BigstackReset(grm);
  } else {
    BigstackReset(bigstack_mark);
  }
  return reterr;
}
#endif

PglErr ScoreListLoad(const char* list_fname, TextStream* score_txsp, LlStr** fname_iterp, uint32_t* infile_ctp) {
  PglErr reterr = kPglRetSuccess;
  {
    LlStr** llstr_endp = fname_iterp;
    uintptr_t infile_ct = 0;
    while (1) {
      const char* fname = TextGet(score_txsp);
      if (!fname) {
        break;
      }
      ++infile_ct;
      const char* fname_end = CurTokenEnd(fname);
      const uint32_t fname_slen = fname_end - fname;
      LlStr* new_entry;
      if (unlikely(bigstack_end_alloc_llstr(fname_slen + 1, &new_entry))) {
        goto ScoreListLoad_ret_NOMEM;
      }
      new_entry->next = nullptr;
      memcpyx(new_entry->str, fname, fname_slen, '\0');
      *llstr_endp = new_entry;
      llstr_endp = &(new_entry->next);
    }
    if (unlikely(!infile_ct)) {
      logerrprintfww("Error: --score-list: %s is empty.\n", list_fname);
      goto ScoreListLoad_ret_MALFORMED_INPUT;
    }
#ifdef __LP64__
    if (unlikely(infile_ct > 0xffffffffU)) {
      logerrprintfww("Error: --score-list: Too many entries in %s.\n", list_fname);
      reterr = kPglRetNotYetSupported;
      goto ScoreListLoad_ret_1;
    }
#endif
    reterr = TextRetarget((*fname_iterp)->str, score_txsp);
    if (unlikely(reterr)) {
      goto ScoreListLoad_ret_TSTREAM_FAIL;
    }
    *infile_ctp = infile_ct;
  }
  while (0) {
  ScoreListLoad_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ScoreListLoad_ret_TSTREAM_FAIL:
    TextStreamErrPrint(list_fname, score_txsp);
    break;
  ScoreListLoad_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
#ifdef __LP64__
 ScoreListLoad_ret_1:
#endif
  return reterr;
}

// to test: do we actually want cur_dosage_ints to be uint64_t* instead of
// uint32_t*?
void FillDdosageInts(const uintptr_t* genovec_buf, const uintptr_t* dosage_present, const Dosage* dosage_main_buf, uint32_t sample_ct, uint32_t dosage_ct, uint32_t is_diploid, uint64_t* cur_ddosage_ints) {
  uint64_t lookup_table[32] ALIGNV16;
  lookup_table[0] = 0;
  lookup_table[2] = (kDosageMid * k1LU) << is_diploid;
  lookup_table[4] = (kDosageMax * k1LU) << is_diploid;
  lookup_table[6] = 0;
  InitLookup16x8bx2(lookup_table);
  GenoarrLookup16x8bx2(genovec_buf, lookup_table, sample_ct, cur_ddosage_ints);
  if (!dosage_ct) {
    return;
  }
  uintptr_t sample_idx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_idx = BitIter1(dosage_present, &sample_idx_base, &cur_bits);
    cur_ddosage_ints[sample_idx] = (dosage_main_buf[dosage_idx] * k1LU) << is_diploid;
  }
}

CONSTI32(kScoreVariantBlockSize, 240);

// must be a power of 2, >= kBitsPerWord
#ifdef USE_AVX2
CONSTI32(kScoreSampleShardSizeAlign, 256);
#else
CONSTI32(kScoreSampleShardSizeAlign, 128);
#endif

// The previous implementation had poor parallelism; the main thread was
// responsible for too many per-sample uint64_t and floating-point operations.
// We now shard those operations across all worker threads, and leave only
// much faster bitvector and similar operations in the main thread.
typedef struct CalcScoreCtxStruct {
  const uintptr_t* variant_include;
  const uint32_t* variant_include_cumulative_popcounts;
  const uintptr_t* qsr_include;
  const uintptr_t* sex_nonmale_collapsed;
  const uintptr_t* sex_female_collapsed;
  uint32_t score_final_col_ct;
  uint32_t sample_shard_size;
  uint32_t sample_ct;
  uint32_t max_difflist_len;
  ScoreFlags flags;
  uint32_t qsr_ct;

  uint32_t* variant_uidxs[2];
  uintptr_t* genovecs[2];
  uintptr_t* raregenos[2];
  uint32_t* difflist_lens[2];
  uint32_t* difflist_sample_ids[2];
  int8_t* difflist_common_genos[2];
  uintptr_t* dosage_presents[2];
  Dosage* dosage_mains[2];
  uint32_t* dosage_cts[2];
  uintptr_t* missing_bitvecs[2];
  uintptr_t* missing_male_bitvecs[2];
  uintptr_t* missing_nonfemale_bitvecs[2];
  double* allele_freqs[2];
  double* geno_slopes[2];
  double* geno_intercepts[2];
  double* score_dense_coefs_cmaj[2];
  double* score_sparse_coefs_vmaj[2];
  unsigned char is_nonx_haploids[2][kScoreVariantBlockSize];
  unsigned char is_relevant_xs[2][kScoreVariantBlockSize];
  unsigned char is_ys[2][kScoreVariantBlockSize];
  unsigned char ploidy_m1s[2][kScoreVariantBlockSize];
  uint32_t cur_variant_batch_size;

  // per-thread buffers
  uint64_t** ddosage_incrs;
  double** dosages_vmaj;

  // final results
  uint64_t* ddosage_sums;
  double** sharded_final_scores_cmaj;
} CalcScoreCtx;

THREAD_FUNC_DECL CalcScoreThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  CalcScoreCtx* ctx = S_CAST(CalcScoreCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const uint32_t* variant_include_cumulative_popcounts = ctx->variant_include_cumulative_popcounts;
  const uintptr_t* qsr_include = ctx->qsr_include;

  const uint32_t score_final_col_ct = ctx->score_final_col_ct;
  uint32_t sample_shard_size = ctx->sample_shard_size;
  const uint32_t sample_idx_start = tidx * sample_shard_size;
  double* shard_final_scores_cmaj = ctx->sharded_final_scores_cmaj[tidx];
  const uint32_t sample_ct = ctx->sample_ct;
  if (sample_shard_size + sample_idx_start > sample_ct) {
    sample_shard_size = sample_ct - sample_idx_start;
  }
  const uint32_t sample_idx_startl2 = NypCtToWordCt(sample_idx_start);
  const uint32_t sample_idx_startl = BitCtToWordCt(sample_idx_start);
  const uint32_t shard_sizel = BitCtToWordCt(sample_shard_size);
  const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  const uint32_t max_difflist_len = ctx->max_difflist_len;
  const uintptr_t raregeno_stride = NypCtToAlignedWordCt(max_difflist_len);
  const uintptr_t difflist_sample_ids_stride = RoundUpPow2(max_difflist_len, kInt32PerVec);
  const uint32_t dosage_main_stride = RoundUpPow2(sample_ct, kDosagePerVec);
  const uintptr_t* shard_sex_nonmale_collapsed = &(ctx->sex_nonmale_collapsed[sample_idx_startl]);
  const uintptr_t* shard_sex_female_collapsed = &(ctx->sex_female_collapsed[sample_idx_startl]);
  const uint32_t shard_nonmale_ct = PopcountWords(shard_sex_nonmale_collapsed, shard_sizel);
  const uint32_t shard_male_ct = sample_shard_size - shard_nonmale_ct;
  const uint32_t shard_female_ct = PopcountWords(shard_sex_female_collapsed, shard_sizel);
  const uint32_t shard_nonfemale_ct = sample_shard_size - shard_female_ct;
  const ScoreFlags flags = ctx->flags;
  const uint32_t model_dominant = (flags / kfScoreDominant) & 1;
  const uint32_t domrec = model_dominant || (flags & kfScoreRecessive);
  const uint32_t no_meanimpute = (flags / kfScoreNoMeanimpute) & 1;
  const uint32_t se_mode = (flags / kfScoreSe) & 1;
  const uint32_t qsr_ct = ctx->qsr_ct;

  uint64_t* ddosage_incrs = ctx->ddosage_incrs[tidx];
  double* dosages_vmaj = ctx->dosages_vmaj[tidx];

  uint64_t* shard_ddosage_sums = &(ctx->ddosage_sums[sample_idx_start]);

  double cur_allele_freq = 0.0;
  double geno_slope = kRecipDosageMax;
  double geno_intercept = 0.0;
  uint32_t parity = 0;
  do {
    const uint32_t cur_variant_batch_size = ctx->cur_variant_batch_size;
    if (cur_variant_batch_size) {
      const uint32_t* variant_uidxs = ctx->variant_uidxs[parity];
      const uintptr_t* genovec_iter = &(ctx->genovecs[parity][sample_idx_startl2]);
      const uintptr_t* raregeno_iter = ctx->raregenos[parity];
      const uint32_t* difflist_lens = ctx->difflist_lens[parity];
      const uint32_t* difflist_sample_ids_iter = ctx->difflist_sample_ids[parity];
      const int8_t* difflist_common_genos = ctx->difflist_common_genos[parity];
      const uintptr_t* dosage_present_iter = ctx->dosage_presents[parity];
      const Dosage* dosage_main_iter = ctx->dosage_mains[parity];
      const uint32_t* dosage_cts = ctx->dosage_cts[parity];
      const uintptr_t* missing_bitvec_iter = &(ctx->missing_bitvecs[parity][sample_idx_startl]);
      const uintptr_t* missing_male_bitvec_iter = &(ctx->missing_male_bitvecs[parity][sample_idx_startl]);
      const uintptr_t* missing_nonfemale_bitvec_iter = &(ctx->missing_nonfemale_bitvecs[parity][sample_idx_startl]);
      const double* allele_freqs = ctx->allele_freqs[parity];
      const double* geno_slopes = ctx->geno_slopes[parity];
      const double* geno_intercepts = ctx->geno_intercepts[parity];
      const double* score_sparse_coefs_vmaj_iter = ctx->score_sparse_coefs_vmaj[parity];
      const unsigned char* is_nonx_haploids = &(ctx->is_nonx_haploids[parity][0]);
      const unsigned char* is_relevant_xs = &(ctx->is_relevant_xs[parity][0]);
      const unsigned char* is_ys = &(ctx->is_ys[parity][0]);
      const unsigned char* ploidy_m1s = &(ctx->ploidy_m1s[parity][0]);
      const uintptr_t* shard_dosage_present = nullptr;
      const Dosage* shard_dosage_main = nullptr;
      double* dosages_vmaj_iter = dosages_vmaj;

      uint32_t sparse_vidx = 0;
      uint32_t shard_dosage_ct = 0;
      for (uint32_t vidx = 0; vidx != cur_variant_batch_size; ++vidx) {
        const uint32_t is_nonx_haploid = is_nonx_haploids[vidx];
        const uint32_t is_relevant_x = is_relevant_xs[vidx];
        const uint32_t is_y = is_ys[vidx];
        if (dosage_cts) {
          // defensive
          shard_dosage_present = nullptr;
          shard_dosage_main = nullptr;
          // bugfix (3 Jun 2022)
          shard_dosage_ct = 0;
          const uint32_t dosage_ct = dosage_cts[vidx];
          if (dosage_ct) {
            const uint32_t dosage_main_offset = PopcountWords(dosage_present_iter, sample_idx_startl);
            if (dosage_main_offset != dosage_ct) {
              shard_dosage_ct = PopcountWords(&(dosage_present_iter[sample_idx_startl]), shard_sizel);
              if (shard_dosage_ct) {
                shard_dosage_present = &(dosage_present_iter[sample_idx_startl]);
                shard_dosage_main = &(dosage_main_iter[dosage_main_offset]);
              }
            }
          }
        } else if (!domrec) {
          // For now, the sparse-genotype optimization is only considered if
          // there are no dosages in the entire .pgen, and
          // 'dominant'/'recessive' wasn't specified.
          // It's actually straightforward to support domrec
          // (modifying difflist_common_geno and raregeno upfront is almost
          // enough), but let's get this working first.
          const uint32_t difflist_common_geno = S_CAST(int32_t, difflist_common_genos[vidx]);
          if (difflist_common_geno != UINT32_MAX) {
            const uint32_t difflist_len = difflist_lens[sparse_vidx];
            if (difflist_len) {
              uint32_t difflist_start_idx = 0;
              if (sample_idx_start) {
                difflist_start_idx = LowerBoundNonemptyU32(difflist_sample_ids_iter, difflist_len, sample_idx_start);
              }
              if (difflist_len > difflist_start_idx) {
                const uint32_t shard_difflist_len = LowerBoundNonemptyU32(&(difflist_sample_ids_iter[difflist_start_idx]), difflist_len - difflist_start_idx, sample_idx_start + sample_shard_size);
                if (shard_difflist_len) {
                  // Shard contains at least one uncommon genotype, so we have
                  // some work to do.
                  if (allele_freqs) {
                    cur_allele_freq = allele_freqs[vidx];
                    geno_slope = geno_slopes[vidx];
                    geno_intercept = geno_intercepts[vidx];
                  }
                  uint64_t ddosage_incr_lookup_table[4];
                  ddosage_incr_lookup_table[0] = 0;
                  ddosage_incr_lookup_table[1] = kDosageMax;
                  ddosage_incr_lookup_table[2] = kDosageMax * 2LL;
                  ddosage_incr_lookup_table[3] = 0;
                  const uint64_t common_ddosage = ddosage_incr_lookup_table[difflist_common_geno];
                  if (common_ddosage) {
                    for (uint32_t uii = 0; uii != 4; ++uii) {
                      // deliberate underflow
                      ddosage_incr_lookup_table[uii] -= common_ddosage;
                    }
                  }
                  double dosage_incr_lookup_table[4];
                  dosage_incr_lookup_table[0] = geno_intercept;
                  const double geno_incr = kDosageMax * geno_slope;
                  dosage_incr_lookup_table[1] = geno_incr + geno_intercept;
                  dosage_incr_lookup_table[2] = dosage_incr_lookup_table[1] + geno_incr;
                  dosage_incr_lookup_table[3] = 0.0;
                  if (!no_meanimpute) {
                    dosage_incr_lookup_table[3] = kDosageMax * 2LL * cur_allele_freq * geno_slope;
                  }
                  if (se_mode) {
                    for (uint32_t uii = 0; uii != 4; ++uii) {
                      dosage_incr_lookup_table[uii] *= dosage_incr_lookup_table[uii];
                    }
                  }
                  const double common_dosage = dosage_incr_lookup_table[difflist_common_geno];
                  for (uint32_t uii = 0; uii != 4; ++uii) {
                    dosage_incr_lookup_table[uii] -= common_dosage;
                  }

                  const uint32_t difflist_end_idx = difflist_start_idx + shard_difflist_len;
                  const uint32_t widx_last = (difflist_end_idx - 1) / kBitsPerWordD2;
                  uint32_t widx = difflist_start_idx / kBitsPerWordD2;
                  const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids_iter[widx * kBitsPerWordD2]);
                  uint32_t difflist_idx_lowbits = difflist_start_idx % kBitsPerWordD2;
                  uint32_t difflist_idx_lowbits_end = kBitsPerWordD2;
                  if (widx == widx_last) {
                    difflist_idx_lowbits_end = ModNz(difflist_end_idx, kBitsPerWordD2);
                  }
                  uintptr_t raregeno_word = raregeno_iter[widx] >> (2 * difflist_idx_lowbits);
                  while (1) {
                    for (; difflist_idx_lowbits != difflist_idx_lowbits_end; ++difflist_idx_lowbits) {
                      const uintptr_t cur_geno = raregeno_word & 3;
                      const uint32_t shard_sample_idx = cur_difflist_sample_ids[difflist_idx_lowbits] - sample_idx_start;
                      uint64_t* ddosage_sum_col_base = &(shard_ddosage_sums[shard_sample_idx]);
                      const uint64_t cur_ddosage_incr = ddosage_incr_lookup_table[cur_geno];
                      if (!qsr_ct) {
                        ddosage_sum_col_base[0] += cur_ddosage_incr;
                      } else {
                        const uint32_t variant_uidx = variant_uidxs[vidx];
                        const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                        for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                          if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                            ddosage_sum_col_base[qsr_idx * sample_ct] += cur_ddosage_incr;
                          }
                        }
                      }
                      double* final_score_col_iter = &(shard_final_scores_cmaj[shard_sample_idx]);
                      const double cur_dosage_incr = dosage_incr_lookup_table[cur_geno];
                      for (uint32_t uii = 0; uii != score_final_col_ct; ++uii) {
                        *final_score_col_iter += cur_dosage_incr * score_sparse_coefs_vmaj_iter[uii];
                        final_score_col_iter = &(final_score_col_iter[sample_shard_size]);
                      }
                      raregeno_word = raregeno_word >> 2;
                    }
                    ++widx;
                    if (widx >= widx_last) {
                      if (widx > widx_last) {
                        break;
                      }
                      difflist_idx_lowbits_end = ModNz(difflist_end_idx, kBitsPerWordD2);
                    }
                    cur_difflist_sample_ids = &(cur_difflist_sample_ids[kBitsPerWordD2]);
                    difflist_idx_lowbits = 0;
                    raregeno_word = raregeno_iter[widx];
                  }
                }
              }
            }

            genovec_iter = &(genovec_iter[sample_ctaw2]);
            missing_bitvec_iter = &(missing_bitvec_iter[sample_ctaw]);
            missing_male_bitvec_iter = &(missing_male_bitvec_iter[sample_ctaw]);
            missing_nonfemale_bitvec_iter = &(missing_nonfemale_bitvec_iter[sample_ctaw]);
            score_sparse_coefs_vmaj_iter = &(score_sparse_coefs_vmaj_iter[score_final_col_ct]);
            ++sparse_vidx;
            raregeno_iter = &(raregeno_iter[raregeno_stride]);
            difflist_sample_ids_iter = &(difflist_sample_ids_iter[difflist_sample_ids_stride]);
            continue;
          }
        }
        FillDdosageInts(genovec_iter, shard_dosage_present, shard_dosage_main, sample_shard_size, shard_dosage_ct, 1 - is_nonx_haploid, ddosage_incrs);
        if (is_y) {
          if (shard_female_ct) {
            uintptr_t shard_sample_idx_base = 0;
            uintptr_t shard_sex_female_collapsed_bits = shard_sex_female_collapsed[0];
            for (uint32_t shard_female_idx = 0; shard_female_idx != shard_female_ct; ++shard_female_idx) {
              const uintptr_t shard_sample_idx = BitIter1(shard_sex_female_collapsed, &shard_sample_idx_base, &shard_sex_female_collapsed_bits);
              ddosage_incrs[shard_sample_idx] = 0;
            }
          }
        } else if (is_relevant_x) {
          uintptr_t shard_sample_idx_base = 0;
          uintptr_t shard_sex_nonmale_collapsed_inv_bits = ~shard_sex_nonmale_collapsed[0];
          for (uint32_t shard_male_idx = 0; shard_male_idx != shard_male_ct; ++shard_male_idx) {
            const uintptr_t shard_sample_idx = BitIter0(shard_sex_nonmale_collapsed, &shard_sample_idx_base, &shard_sex_nonmale_collapsed_inv_bits);
            ddosage_incrs[shard_sample_idx] /= 2;
          }
        }
        if (domrec) {
          if (model_dominant) {
            for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
              if (ddosage_incrs[shard_sample_idx] > kDosageMax) {
                ddosage_incrs[shard_sample_idx] = kDosageMax;
              }
            }
          } else {
            // recessive
            for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
              uint64_t cur_ddosage_incr = ddosage_incrs[shard_sample_idx];
              if (cur_ddosage_incr <= kDosageMax) {
                cur_ddosage_incr = 0;
              } else {
                cur_ddosage_incr -= kDosageMax;
              }
              ddosage_incrs[shard_sample_idx] = cur_ddosage_incr;
            }
          }
        }
        if (!qsr_ct) {
          for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
            shard_ddosage_sums[shard_sample_idx] += ddosage_incrs[shard_sample_idx];
          }
        } else {
          const uint32_t variant_uidx = variant_uidxs[vidx];
          const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
          for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
            if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
              uint64_t* cur_ddosage_sums = &(shard_ddosage_sums[qsr_idx * sample_ct]);
              for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
                cur_ddosage_sums[shard_sample_idx] += ddosage_incrs[shard_sample_idx];
              }
            }
          }
        }
        if (allele_freqs) {
          cur_allele_freq = allele_freqs[vidx];
          geno_slope = geno_slopes[vidx];
          geno_intercept = geno_intercepts[vidx];
        }
        if ((!shard_dosage_ct) && (!is_nonx_haploid) && (!is_relevant_x)) {
          // Fast path for common no-dosage case.
          // (We wait till this point to branch since the shard_ddosage_sums
          // update isn't that expensive.)
          double lookup_table[32] ALIGNV16;
          lookup_table[0] = geno_intercept;
          const double geno_incr = kDosageMax * geno_slope;
          const double geno_one = geno_incr + geno_intercept;
          if (!domrec) {
            lookup_table[2] = geno_one;
            lookup_table[4] = geno_incr + geno_one;
          } else {
            lookup_table[2] = model_dominant? geno_one : geno_intercept;
            lookup_table[4] = geno_one;
          }
          lookup_table[6] = 0.0;
          if (!no_meanimpute) {
            double missing_effect = kDosageMax * cur_allele_freq * geno_slope;
            if (!domrec) {
              missing_effect *= 2;
            }
            lookup_table[6] = missing_effect;
          }
          if (se_mode) {
            lookup_table[0] *= lookup_table[0];
            lookup_table[2] *= lookup_table[2];
            lookup_table[4] *= lookup_table[4];
            lookup_table[6] *= lookup_table[6];
          }
          InitLookup16x8bx2(lookup_table);
          GenoarrLookup16x8bx2(genovec_iter, lookup_table, sample_shard_size, dosages_vmaj_iter);
        } else {
          const uint32_t shard_missing_ct = PopcountWords(missing_bitvec_iter, shard_sizel);
          if (shard_missing_ct) {
            double missing_effect = 0.0;
            if (!no_meanimpute) {
              missing_effect = kDosageMax * cur_allele_freq * geno_slope;
            }
            if (is_y) {
              ZeroDArr(sample_shard_size, dosages_vmaj_iter);
              if ((!no_meanimpute) && shard_nonfemale_ct) {
                for (uint32_t shard_widx = 0; shard_widx != shard_sizel; ++shard_widx) {
                  uintptr_t cur_missing_nonfemale_bits = missing_nonfemale_bitvec_iter[shard_widx];
                  if (!cur_missing_nonfemale_bits) {
                    continue;
                  }
                  double* cur_dosages_vmaj_iter = &(dosages_vmaj_iter[shard_widx * kBitsPerWord]);
                  do {
                    const uint32_t sample_idx_lowbits = ctzw(cur_missing_nonfemale_bits);
                    cur_dosages_vmaj_iter[sample_idx_lowbits] = missing_effect;
                    cur_missing_nonfemale_bits &= cur_missing_nonfemale_bits - 1;
                  } while (cur_missing_nonfemale_bits);
                }
              }
            } else if (is_relevant_x) {
              ZeroDArr(sample_shard_size, dosages_vmaj_iter);
              if (!no_meanimpute) {
                if (shard_male_ct) {
                  for (uint32_t shard_widx = 0; shard_widx != shard_sizel; ++shard_widx) {
                    uintptr_t cur_missing_male_bits = missing_male_bitvec_iter[shard_widx];
                    if (!cur_missing_male_bits) {
                      continue;
                    }
                    double* cur_dosages_vmaj_iter = &(dosages_vmaj_iter[shard_widx * kBitsPerWord]);
                    do {
                      const uint32_t sample_idx_lowbits = ctzw(cur_missing_male_bits);
                      cur_dosages_vmaj_iter[sample_idx_lowbits] = missing_effect;
                      cur_missing_male_bits &= cur_missing_male_bits - 1;
                    } while (cur_missing_male_bits);
                  }
                }
                if (shard_nonmale_ct) {
                  missing_effect *= 2;
                  for (uint32_t shard_widx = 0; shard_widx != shard_sizel; ++shard_widx) {
                    uintptr_t cur_missing_nonmale_bits = missing_bitvec_iter[shard_widx] & shard_sex_nonmale_collapsed[shard_widx];
                    if (!cur_missing_nonmale_bits) {
                      continue;
                    }
                    double* cur_dosages_vmaj_iter = &(dosages_vmaj_iter[shard_widx * kBitsPerWord]);
                    do {
                      const uint32_t sample_idx_lowbits = ctzw(cur_missing_nonmale_bits);
                      cur_dosages_vmaj_iter[sample_idx_lowbits] = missing_effect;
                      cur_missing_nonmale_bits &= cur_missing_nonmale_bits - 1;
                    } while (cur_missing_nonmale_bits);
                  }
                }
              }
            } else {
              const double ploidy_d = ploidy_m1s[vidx]? 2.0 : 1.0;
              missing_effect *= ploidy_d;
              uintptr_t shard_sample_idx_base = 0;
              uintptr_t missing_bits = missing_bitvec_iter[0];
              for (uint32_t missing_idx = 0; missing_idx != shard_missing_ct; ++missing_idx) {
                const uintptr_t shard_sample_idx = BitIter1(missing_bitvec_iter, &shard_sample_idx_base, &missing_bits);
                dosages_vmaj_iter[shard_sample_idx] = missing_effect;
              }
            }
          }
          const uint32_t shard_nm_sample_ct = sample_shard_size - shard_missing_ct;
          uintptr_t shard_sample_idx_base = 0;
          uintptr_t missing_inv_bits = ~missing_bitvec_iter[0];
          for (uint32_t shard_nm_sample_idx = 0; shard_nm_sample_idx != shard_nm_sample_ct; ++shard_nm_sample_idx) {
            const uintptr_t shard_sample_idx = BitIter0(missing_bitvec_iter, &shard_sample_idx_base, &missing_inv_bits);
            dosages_vmaj_iter[shard_sample_idx] = u63tod(ddosage_incrs[shard_sample_idx]) * geno_slope + geno_intercept;
          }
          if (se_mode) {
            // Suppose our score coefficients are drawn from independent
            // Gaussians.  Then the variance of the final score average is the
            // sum of the variances of the individual terms, divided by (T^2)
            // where T is the number of terms.  These individual variances are
            // of the form (<genotype value> * <stdev>)^2.
            //
            // Thus, we can use the same inner loop to compute standard errors,
            // as long as
            //   1. we square the genotypes and the standard errors before
            //      matrix multiplication, and
            //   2. we take the square root of the sums at the end.
            for (uint32_t shard_sample_idx = 0; shard_sample_idx != sample_shard_size; ++shard_sample_idx) {
              dosages_vmaj_iter[shard_sample_idx] *= dosages_vmaj_iter[shard_sample_idx];
            }
          }
        }
        dosages_vmaj_iter = &(dosages_vmaj_iter[sample_shard_size]);

        genovec_iter = &(genovec_iter[sample_ctaw2]);
        missing_bitvec_iter = &(missing_bitvec_iter[sample_ctaw]);
        missing_male_bitvec_iter = &(missing_male_bitvec_iter[sample_ctaw]);
        missing_nonfemale_bitvec_iter = &(missing_nonfemale_bitvec_iter[sample_ctaw]);
        if (dosage_present_iter) {
          dosage_present_iter = &(dosage_present_iter[sample_ctaw]);
          dosage_main_iter = &(dosage_main_iter[dosage_main_stride]);
        }
      }

      const uint32_t dense_variant_ct = cur_variant_batch_size - sparse_vidx;
      if (dense_variant_ct) {
        const double* score_dense_coefs_cmaj = ctx->score_dense_coefs_cmaj[parity];
        RowMajorMatrixMultiplyStridedIncr(score_dense_coefs_cmaj, dosages_vmaj, score_final_col_ct, kScoreVariantBlockSize, sample_shard_size, sample_shard_size, dense_variant_ct, sample_shard_size, shard_final_scores_cmaj);
      }
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

typedef struct ParsedQscoreRangeStruct {
  char* range_name;
  double lbound;
  double ubound;
} ParsedQscoreRange;

PglErr ScoreReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const double* allele_freqs, const ScoreInfo* score_info_ptr, const char* output_missing_pheno, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t nosex_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t xchr_model, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  const char* cur_input_fname = nullptr;
  char* cswritep = nullptr;
  FILE* score_tmpfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream score_txs;
  ThreadGroup tg;
  CompressStreamState css;
  PreinitTextStream(&score_txs);
  PreinitThreads(&tg);
  PreinitCstream(&css);
  {
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    if (!xchr_model) {
      uint32_t x_code;
      if (XymtExists(cip, kChrOffsetX, &x_code)) {
        uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
        uint32_t x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
        uint32_t x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
        if (!AllBitsAreZero(variant_include, x_start, x_end)) {
          uintptr_t* variant_include_no_x;
          if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include_no_x))) {
            goto ScoreReport_ret_NOMEM;
          }
          memcpy(variant_include_no_x, variant_include, raw_variant_ctl * sizeof(intptr_t));
          ClearBitsNz(x_start, x_end, variant_include_no_x);
          variant_include = variant_include_no_x;
        }
      }
    } else if (xchr_model == 2) {
      xchr_model = 0;
    }
    // now xchr_model is set iff it's 1

    const ScoreFlags flags = score_info_ptr->flags;
    const uint32_t multi_input = (flags / kfScoreMultiInput) & 1;
    if (multi_input) {
      assert(!(flags & (kfScoreListVariants | kfScoreColNallele | kfScoreColDenom | kfScoreColDosageSum)));
    }
    const uint32_t output_zst = (flags / kfScoreZs) & 1;
    uint32_t* variant_id_htable = nullptr;
    uint32_t variant_id_htable_size = 0;
    uint32_t* variant_include_cumulative_popcounts = nullptr;
    uintptr_t* qsr_include = nullptr;
    char** range_names = nullptr;
    uintptr_t qsr_ct = 0;
    if (multi_input || score_info_ptr->qsr_range_fname) {
      // Limit this to ~1/8 of available memory, since memory may be tight.
      variant_id_htable_size = GetHtableFastSize(variant_ct);
      const uintptr_t htable_size_limit = bigstack_left() / (8 * sizeof(int32_t));
      if (variant_id_htable_size > htable_size_limit) {
        variant_id_htable_size = htable_size_limit;
        const uint32_t htable_size_min = GetHtableMinSize(variant_ct);
        if (htable_size_min > variant_id_htable_size) {
          variant_id_htable_size = htable_size_min;
        }
      }
      if (unlikely(bigstack_alloc_u32(variant_id_htable_size, &variant_id_htable))) {
        goto ScoreReport_ret_NOMEM;
      }
      reterr = PopulateIdHtableMt(nullptr, variant_include, variant_ids, variant_ct, 0, variant_id_htable_size, max_thread_ct, nullptr, variant_id_htable, nullptr);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
    }
    if (score_info_ptr->qsr_range_fname) {
      // Strictly speaking, textFILE would be more appropriate for the range
      // file since it should be tiny, but it doesn't really matter.
      // We still reserve bigstack_left() / 8 for the line-buffer since the
      // data file usually contains allele codes, and we use TextRetarget()
      // below to use the buffer allocated here for the data file too (and for
      // the score file later).
      reterr = SizeAndInitTextStream(score_info_ptr->qsr_range_fname, bigstack_left() / 8, 1, &score_txs);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_QSR_RANGE_TSTREAM_FAIL;
      }
      unsigned char* bigstack_mark2 = g_bigstack_base;
      // strlen("<prefix>.<range name>.sscore[.zst]") < kPglFnamesize
      const uint32_t max_name_slen = kPglFnamesize - S_CAST(uintptr_t, outname_end - outname) - 9 - output_zst * 4;
      ParsedQscoreRange* parsed_qscore_ranges = R_CAST(ParsedQscoreRange*, g_bigstack_base);
      unsigned char* tmp_alloc_end = g_bigstack_end;
      uintptr_t miss_ct = 0;
      while (1) {
        ++line_idx;
        const char* line_start = TextGet(&score_txs);
        if (!line_start) {
          if (likely(!TextStreamErrcode2(&score_txs, &reterr))) {
            break;
          }
          goto ScoreReport_ret_QSR_RANGE_TSTREAM_FAIL;
        }
        // range name, p-value lower bound, p-value upper bound
        const char* range_name_end = CurTokenEnd(line_start);
        const char* lbound_start = FirstNonTspace(range_name_end);
        double lbound;
        // PLINK 1.9 documentation promises that lines with too few entries or
        // nonnumeric values in the second and third column are ignored.
        const char* lbound_end = ScantokDouble(lbound_start, &lbound);
        if (!lbound_end) {
          continue;
        }
        const char* ubound_start = FirstNonTspace(lbound_end);
        double ubound;
        const char* ubound_end = ScantokDouble(ubound_start, &ubound);
        if (!ubound_end) {
          continue;
        }
        if (unlikely(lbound > ubound)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Upper bound < lower bound on line %" PRIuPTR " of --q-score-range range file.\n", line_idx);
          goto ScoreReport_ret_MALFORMED_INPUT_WW;
        }
        const uint32_t name_slen = range_name_end - line_start;
        if (name_slen > max_name_slen) {
          snprintf(g_logbuf, kLogbufSize, "Error: Name too long on line %" PRIuPTR " of --q-score-range range file.\n", line_idx);
        }
        unsigned char* tmp_alloc_base = R_CAST(unsigned char*, &(parsed_qscore_ranges[qsr_ct]));
        if (S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) <= name_slen + sizeof(ParsedQscoreRange)) {
          goto ScoreReport_ret_NOMEM;
        }
        tmp_alloc_end -= name_slen + 1;
        char* stored_name = R_CAST(char*, tmp_alloc_end);
        memcpyx(stored_name, line_start, name_slen, '\0');
        parsed_qscore_ranges[qsr_ct].range_name = stored_name;
        parsed_qscore_ranges[qsr_ct].lbound = lbound;
        parsed_qscore_ranges[qsr_ct].ubound = ubound;
        ++qsr_ct;
      }
      if (unlikely(!qsr_ct)) {
        logerrputs("Error: Empty --q-score-range range file.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
      BigstackEndSet(tmp_alloc_end);
      if (unlikely(BigstackBaseSetChecked(&(parsed_qscore_ranges[qsr_ct])))) {
        goto ScoreReport_ret_NOMEM;
      }
#ifndef LAPACK_ILP64
      if (unlikely(qsr_ct > (0x7fffffff / kScoreVariantBlockSize))) {
        logerrputs("Error: --q-score-range range count too large for this " PROG_NAME_STR " build.  If this is\nreally the computation you want, use a " PROG_NAME_STR " build with large-matrix support.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
#  ifndef __LP64__
      const uint64_t bit_ct = S_CAST(uint64_t, qsr_ct) * variant_ct;
      if (unlikely(bit_ct > 0xffffffffU)) {
        goto ScoreReport_ret_NOMEM;
      }
#  endif
#endif
      if (unlikely(bigstack_end_alloc_u32(raw_variant_ctl, &variant_include_cumulative_popcounts) ||
                   bigstack_end_calloc_w(BitCtToWordCt(S_CAST(uint64_t, qsr_ct) * variant_ct), &qsr_include) ||
                   bigstack_end_alloc_cp(qsr_ct, &range_names))) {
        goto ScoreReport_ret_NOMEM;
      }
      for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
        range_names[qsr_idx] = parsed_qscore_ranges[qsr_idx].range_name;
      }
      const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
      uintptr_t* already_seen;
      if (unlikely(bigstack_calloc_w(variant_ctl, &already_seen))) {
        goto ScoreReport_ret_NOMEM;
      }
      FillCumulativePopcounts(variant_include, raw_variant_ctl, variant_include_cumulative_popcounts);
      reterr = TextRetarget(score_info_ptr->qsr_data_fname, &score_txs);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_QSR_RANGE_TSTREAM_FAIL;
      }
      const uint32_t colid_first = (score_info_ptr->qsr_varid_col_p1 < score_info_ptr->qsr_val_col_p1);
      uint32_t colmin;
      uint32_t coldiff;
      if (colid_first) {
        colmin = score_info_ptr->qsr_varid_col_p1 - 1;
        coldiff = score_info_ptr->qsr_val_col_p1 - score_info_ptr->qsr_varid_col_p1;
      } else {
        colmin = score_info_ptr->qsr_val_col_p1 - 1;
        coldiff = score_info_ptr->qsr_varid_col_p1 - score_info_ptr->qsr_val_col_p1;
      }
      line_idx = 0;
      miss_ct = 0;
      if (flags & kfScoreQsrHeader) {
        ++line_idx;
        if (unlikely(!TextGet(&score_txs))) {
          if (!TextStreamErrcode2(&score_txs, &reterr)) {
            logerrputs("Error: Empty --q-score-range data file.\n");
            goto ScoreReport_ret_MALFORMED_INPUT;
          }
          goto ScoreReport_ret_QSR_DATA_TSTREAM_FAIL;
        }
      }
      double* min_vals = nullptr;
      if (flags & kfScoreQsrMin) {
        // something like this is needed to handle --glm output for
        // multiallelic variants.
        // (possible todo: --glm modifier which requests all-allele joint tests
        // for multiallelic variants)
        if (unlikely(bigstack_alloc_d(variant_ct, &min_vals))) {
          goto ScoreReport_ret_NOMEM;
        }
      }
      while (1) {
        ++line_idx;
        const char* line_start = TextGet(&score_txs);
        if (!line_start) {
          if (likely(!TextStreamErrcode2(&score_txs, &reterr))) {
            break;
          }
          goto ScoreReport_ret_QSR_DATA_TSTREAM_FAIL;
        }
        const char* colid_ptr;
        const char* colval_ptr;
        if (colid_first) {
          colid_ptr = NextTokenMult0(line_start, colmin);
          colval_ptr = NextTokenMult(colid_ptr, coldiff);
          if (unlikely(!colval_ptr)) {
            goto ScoreReport_ret_QSR_DATA_MISSING_TOKENS;
          }
        } else {
          colval_ptr = NextTokenMult0(line_start, colmin);
          colid_ptr = NextTokenMult(colval_ptr, coldiff);
          if (unlikely(!colid_ptr)) {
            goto ScoreReport_ret_QSR_DATA_MISSING_TOKENS;
          }
        }
        const uint32_t varid_slen = strlen_se(colid_ptr);
        const uint32_t variant_uidx = VariantIdDupflagHtableFind(colid_ptr, variant_ids, variant_id_htable, varid_slen, variant_id_htable_size, max_variant_id_slen);
        if ((variant_uidx >> 31) || (!IsSet(variant_include, variant_uidx))) {
          ++miss_ct;
          continue;
        }
        double cur_val;
        if (!ScantokDouble(colval_ptr, &cur_val)) {
          // Tolerate NA without erroring out.  (Could count this as seen, but
          // that would force the min_vals logic to be more complicated.)
          const char* colval_end = CurTokenEnd(colval_ptr);
          if (likely(IsNanStr(colval_ptr, colval_end - colval_ptr))) {
            continue;
          }
          *K_CAST(char*, colval_end) = '\0';
          logerrprintfww("Error: Invalid value '%s' on line %" PRIuPTR " of --q-score-range data file.\n", colval_ptr, line_idx);
          goto ScoreReport_ret_MALFORMED_INPUT;
        }
        const uint32_t variant_idx = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx);
        const uintptr_t bit_idx_base = variant_idx * qsr_ct;
        if (min_vals) {
          if (IsSet(already_seen, variant_idx)) {
            if (min_vals[variant_idx] <= cur_val) {
              continue;
            }
            ClearBitsNz(bit_idx_base, bit_idx_base + qsr_ct, qsr_include);
          }
          min_vals[variant_idx] = cur_val;
        } else {
          if (IsSet(already_seen, variant_idx)) {
            logerrprintfww("Error: Duplicate variant ID '%s' in --q-score-range data file. (Add the 'min' modifier if this is a multiallelic variant that you want to use the minimum p-value for.)\n", variant_ids[variant_uidx]);
            goto ScoreReport_ret_MALFORMED_INPUT;
          }
        }
        SetBit(variant_idx, already_seen);
        for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
          if ((cur_val < parsed_qscore_ranges[qsr_idx].lbound) || (cur_val > parsed_qscore_ranges[qsr_idx].ubound)) {
            continue;
          }
          SetBit(bit_idx_base + qsr_idx, qsr_include);
        }
      }
      const uint32_t qsr_variant_ct = PopcountWords(already_seen, variant_ctl);
      if (unlikely(!qsr_variant_ct)) {
        logerrputs("Error: No valid entries in --q-score-range data file.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
      logprintf("--q-score-range: %" PRIuPTR " range%s and %u variant%s loaded.\n", qsr_ct, (qsr_ct == 1)? "" : "s", qsr_variant_ct, (qsr_variant_ct == 1)? "" : "s");
      if (miss_ct) {
        logerrprintf("Warning: %" PRIuPTR " line%s skipped in --q-score-range data file.\n", miss_ct, (miss_ct == 1)? "" : "s");
      }
      // possible todo: replace variant_include with already_seen, and compact
      // qsr_include.
      // but for now, we just free already_seen, and in the common use cases
      // this should be fine.
      reterr = TextRetarget(score_info_ptr->input_fname, &score_txs);
      if (unlikely(reterr)) {
        cur_input_fname = score_info_ptr->input_fname;
        goto ScoreReport_ret_TSTREAM_FAIL;
      }
      BigstackReset(bigstack_mark2);
      line_idx = 0;
    } else {
      reterr = SizeAndInitTextStream(score_info_ptr->input_fname, bigstack_left() / 8, 1, &score_txs);
      if (unlikely(reterr)) {
        cur_input_fname = score_info_ptr->input_fname;
        goto ScoreReport_ret_TSTREAM_FAIL;
      }
    }

    uint32_t varid_col_idx = score_info_ptr->varid_col_p1 - 1;
    uint32_t allele_col_idx = score_info_ptr->allele_col_p1 - 1;
    uint32_t relevant_col_idx_end = MAXV(varid_col_idx, allele_col_idx);
    uint32_t* score_col_idx_deltas = nullptr;
    uintptr_t score_col_ct = 1;
    if (!score_info_ptr->input_col_idx_range_list.name_ct) {
      // catch edge case
      const uint32_t score_col_idx = allele_col_idx + 1;
      if (unlikely(score_col_idx == varid_col_idx)) {
        logerrprintf("Error: --score%s variant ID column index matches a coefficient column index.\n", multi_input? "-list" : "");
        goto ScoreReport_ret_INVALID_CMDLINE;
      }
      relevant_col_idx_end += (score_col_idx > varid_col_idx);
      if (unlikely(bigstack_alloc_u32(1, &score_col_idx_deltas))) {
        goto ScoreReport_ret_NOMEM;
      }
      score_col_idx_deltas[0] = score_col_idx;
    } else {
      const uint32_t range_list_max = NumericRangeListMax(&(score_info_ptr->input_col_idx_range_list));
      if (range_list_max > relevant_col_idx_end) {
        relevant_col_idx_end = range_list_max;
      }
      unsigned char* bigstack_end_mark2 = g_bigstack_end;
      const uint32_t last_col_idxl = BitCtToWordCt(relevant_col_idx_end);
      uintptr_t* score_col_bitarr;
      if (unlikely(bigstack_end_calloc_w(last_col_idxl, &score_col_bitarr))) {
        goto ScoreReport_ret_NOMEM;
      }
      if (unlikely(NumericRangeListToBitarr(&(score_info_ptr->input_col_idx_range_list), relevant_col_idx_end, 1, 0, score_col_bitarr))) {
        goto ScoreReport_ret_MISSING_TOKENS;
      }
      if (unlikely(IsSet(score_col_bitarr, varid_col_idx))) {
        logerrprintf("Error: --score%s variant ID column index matches a coefficient column index.\n", multi_input? "-list" : "");
        goto ScoreReport_ret_INVALID_CMDLINE;
      }
      if (unlikely(IsSet(score_col_bitarr, allele_col_idx))) {
        logerrprintf("Error: --score%s allele column index matches a coefficient column index.\n", multi_input? "-list" : "");
        goto ScoreReport_ret_INVALID_CMDLINE;
      }
      score_col_ct = PopcountWords(score_col_bitarr, last_col_idxl);
      if (unlikely(bigstack_alloc_u32(score_col_ct, &score_col_idx_deltas))) {
        goto ScoreReport_ret_NOMEM;
      }
      uintptr_t col_uidx_base = 0;
      uintptr_t score_col_bitarr_bits = score_col_bitarr[0];
      for (uintptr_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
        const uint32_t col_uidx = BitIter1(score_col_bitarr, &col_uidx_base, &score_col_bitarr_bits);
        score_col_idx_deltas[score_col_idx] = col_uidx;
      }
      // now convert to deltas
      for (uintptr_t score_col_idx = score_col_ct - 1; score_col_idx; --score_col_idx) {
        score_col_idx_deltas[score_col_idx] -= score_col_idx_deltas[score_col_idx - 1];
      }
      BigstackEndReset(bigstack_end_mark2);
    }

    LlStr* fname_iter = nullptr;
    uint32_t infile_ct = 1;
    if (multi_input) {
      reterr = ScoreListLoad(score_info_ptr->input_fname, &score_txs, &fname_iter, &infile_ct);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
    } else {
      const uint32_t fname_blen = strlen(score_info_ptr->input_fname) + 1;
      if (unlikely(bigstack_alloc_llstr(fname_blen, &fname_iter))) {
        goto ScoreReport_ret_NOMEM;
      }
      fname_iter->next = nullptr;
      memcpy(fname_iter->str, score_info_ptr->input_fname, fname_blen);
    }

    uintptr_t score_final_col_ct = score_col_ct;
    if (qsr_ct) {
      const uint64_t prod = S_CAST(uint64_t, qsr_ct) * score_col_ct;
      if (prod > 0x7fffffff) {
        // little point in supporting this even in large-matrix build
        logerrputs("Error: <--score[-list] column count> * <--q-score-range range count> too large.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
#ifndef LAPACK_ILP64
      if (unlikely(prod > (0x7fffffff / kScoreVariantBlockSize))) {
        logerrputs("Error: <--score[-list] column count> * <--q-score-range range count> too large\nfor this " PROG_NAME_STR " build.  If this is really the computation you want, use a\n" PROG_NAME_STR "\nbuild with large-matrix support.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
#endif
      score_final_col_ct = qsr_ct * score_col_ct;
#ifndef LAPACK_ILP64
    } else {
      if (unlikely(score_final_col_ct > (0x7fffffff / kScoreVariantBlockSize))) {
        logerrputs("Error: --score[-list] column count too large for this " PROG_NAME_STR " build.  If this is\nreally the computation you want, use a " PROG_NAME_STR " build with large-matrix support.\n");
        goto ScoreReport_ret_INCONSISTENT_INPUT;
      }
#endif
    }
    const uint32_t calc_thread_ct = MAXV(1, max_thread_ct - 1);
    uint32_t sample_shard_size;
    uint32_t sample_shard_ct;
    if (sample_ct <= calc_thread_ct * kScoreSampleShardSizeAlign) {
      sample_shard_size = kScoreSampleShardSizeAlign;
      sample_shard_ct = DivUp(sample_ct, kScoreSampleShardSizeAlign);
    } else {
      sample_shard_size = kScoreSampleShardSizeAlign * DivUp(sample_ct, calc_thread_ct * kScoreSampleShardSizeAlign);
      sample_shard_ct = DivUp(sample_ct, sample_shard_size);
    }
#ifndef LAPACK_ILP64
    if (unlikely(score_final_col_ct * sample_shard_size > 0x7fffffff)) {
      logerrputs("Error: <--score[-list] column count> * <sample count> too large for this " PROG_NAME_STR "\nbuild.  If this is really the computation you want, use a " PROG_NAME_STR " build with\nlarge-matrix support.\n");
      goto ScoreReport_ret_INCONSISTENT_INPUT;
    }
#endif
    const uint32_t sample_shard_ct_m1 = sample_shard_ct - 1;
    CalcScoreCtx ctx;
    ctx.variant_include = qsr_ct? variant_include : nullptr;
    ctx.variant_include_cumulative_popcounts = variant_include_cumulative_popcounts;
    ctx.qsr_include = qsr_include;
    ctx.score_final_col_ct = score_final_col_ct;
    ctx.sample_shard_size = sample_shard_size;
    ctx.sample_ct = sample_ct;
    ctx.flags = flags;
    ctx.qsr_ct = qsr_ct;
    ctx.cur_variant_batch_size = kScoreVariantBlockSize;
    if (unlikely(SetThreadCt(sample_shard_ct, &tg))) {
      goto ScoreReport_ret_NOMEM;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
    const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
    const uintptr_t dosage_main_stride = RoundUpPow2(sample_ct, kDosagePerVec);
    const uint32_t acc1_vec_ct = BitCtToVecCt(sample_ct);
    const uint32_t acc4_vec_ct = acc1_vec_ct * 4;
    const uint32_t acc8_vec_ct = acc1_vec_ct * 8;
    const uint32_t write_score_avgs = (flags / kfScoreColScoreAvgs) & 1;
    const uint32_t write_score_sums = (flags / kfScoreColScoreSums) & 1;
    if (unlikely(S_CAST(uint64_t, score_col_ct) * infile_ct > 0x7fffff)) {
      goto ScoreReport_ret_TOO_MANY_COLUMNS;
    }
    const uintptr_t global_score_col_ct = score_col_ct * infile_ct;
    const uintptr_t overflow_buf_size = RoundUpPow2((global_score_col_ct * (write_score_avgs + write_score_sums) + pheno_ct) * 16 + 3 * kMaxIdSlen + kCompressStreamBlock + 64, kCacheline);
    if (unlikely(overflow_buf_size > S_CAST(uintptr_t, kMaxLongLine) + S_CAST(uintptr_t, kCompressStreamBlock))) {
      goto ScoreReport_ret_TOO_MANY_COLUMNS;
    }
    uintptr_t overflow_buf_alloc = overflow_buf_size;
    if (flags & (kfScoreZs | kfScoreListVariantsZs)) {
      overflow_buf_alloc += CstreamWkspaceReq(overflow_buf_size);
    }
    uintptr_t raw_allele_ct = 2 * raw_variant_ct;
    if (allele_idx_offsets) {
      raw_allele_ct = allele_idx_offsets[raw_variant_ct];
    }
    const uintptr_t raw_allele_ctl = BitCtToWordCt(raw_allele_ct);
    const uintptr_t qsr_ct_nz = qsr_ct + (qsr_ct == 0);
    uint32_t* sample_include_cumulative_popcounts = nullptr;
    uintptr_t* sex_nonmale_collapsed = nullptr;
    VecW* missing_diploid_accx = nullptr;
    VecW* missing_haploid_accx = nullptr;
    uint32_t* allele_ct_bases;
    uint32_t* male_allele_ct_decrs;
    uint32_t* nonfemale_allele_ct_incrs;
    uint32_t* variant_ct_rems;
    uint32_t* variant_hap_ct_rems;
    uintptr_t* already_seen_variants;
    uintptr_t* already_seen_alleles;
    double* tmpfile_buf;
    char* overflow_buf = nullptr;
    if (unlikely(bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw2, &ctx.genovecs[0]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw2, &ctx.genovecs[1]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_bitvecs[0]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_bitvecs[1]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_male_bitvecs[0]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_male_bitvecs[1]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_nonfemale_bitvecs[0]) ||
                 bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.missing_nonfemale_bitvecs[1]) ||
                 bigstack_alloc_d(kScoreVariantBlockSize * score_final_col_ct, &(ctx.score_dense_coefs_cmaj[0])) ||
                 bigstack_alloc_d(kScoreVariantBlockSize * score_final_col_ct, &(ctx.score_dense_coefs_cmaj[1])) ||
                 bigstack_alloc_dp(sample_shard_ct, &ctx.sharded_final_scores_cmaj) ||
                 bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_w(sample_ctl, &sex_nonmale_collapsed) ||
                 bigstack_end_clalloc_v(11 * acc4_vec_ct * qsr_ct_nz, &missing_diploid_accx) ||
                 bigstack_end_clalloc_v(11 * acc4_vec_ct * qsr_ct_nz, &missing_haploid_accx) ||
                 bigstack_alloc_u32(qsr_ct_nz, &allele_ct_bases) ||
                 bigstack_alloc_u32(qsr_ct_nz, &male_allele_ct_decrs) ||
                 bigstack_alloc_u32(qsr_ct_nz, &nonfemale_allele_ct_incrs) ||
                 bigstack_alloc_u32(qsr_ct_nz * 2, &variant_ct_rems) ||
                 bigstack_alloc_u32(qsr_ct_nz * 2, &variant_hap_ct_rems) ||
                 bigstack_alloc_u64p(sample_shard_ct, &ctx.ddosage_incrs) ||
                 bigstack_alloc_dp(sample_shard_ct, &ctx.dosages_vmaj) ||
                 bigstack_alloc_u64(qsr_ct_nz * sample_ct, &ctx.ddosage_sums) ||
                 bigstack_alloc_w(raw_variant_ctl, &already_seen_variants) ||
                 bigstack_alloc_w(raw_allele_ctl, &already_seen_alleles) ||
                 bigstack_end_clalloc_d(score_col_ct, &tmpfile_buf) ||
                 bigstack_end_clalloc_c(overflow_buf_alloc, &overflow_buf))) {
      goto ScoreReport_ret_NOMEM;
    }
    {
      uintptr_t cur_sample_shard_size = sample_shard_size;
      for (uint32_t shard_idx = 0; ; ++shard_idx) {
        if (shard_idx >= sample_shard_ct_m1) {
          if (shard_idx > sample_shard_ct_m1) {
            break;
          }
          cur_sample_shard_size = sample_ct - sample_shard_ct_m1 * sample_shard_size;
        }
        if (unlikely(bigstack_alloc_d(score_final_col_ct * cur_sample_shard_size, &ctx.sharded_final_scores_cmaj[shard_idx]) ||
                     bigstack_alloc_u64(cur_sample_shard_size, &ctx.ddosage_incrs[shard_idx]) ||
                     bigstack_alloc_d(kScoreVariantBlockSize * cur_sample_shard_size, &ctx.dosages_vmaj[shard_idx]))) {
          goto ScoreReport_ret_NOMEM;
        }
      }
    }
    if (qsr_ct) {
      if (unlikely(bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.variant_uidxs[0]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.variant_uidxs[1]))) {
        goto ScoreReport_ret_NOMEM;
      }
    } else {
      ctx.variant_uidxs[0] = nullptr;
      ctx.variant_uidxs[1] = nullptr;
    }
    if (PgrGetGflags(simple_pgrp) & kfPgenGlobalDosagePresent) {
      // We may want to buffer less than 2x240 variants at a time in this case
      // if there are >500k samples...
      if (unlikely(bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.dosage_presents[0]) ||
                   bigstack_alloc_w(kScoreVariantBlockSize * sample_ctaw, &ctx.dosage_presents[1]) ||
                   bigstack_alloc_dosage(kScoreVariantBlockSize * dosage_main_stride, &ctx.dosage_mains[0]) ||
                   bigstack_alloc_dosage(kScoreVariantBlockSize * dosage_main_stride, &ctx.dosage_mains[1]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.dosage_cts[0]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.dosage_cts[1]))) {
        goto ScoreReport_ret_NOMEM;
      }
    } else {
      ctx.dosage_presents[0] = nullptr;
      ctx.dosage_presents[1] = nullptr;
      ctx.dosage_mains[0] = nullptr;
      ctx.dosage_mains[1] = nullptr;
      ctx.dosage_cts[0] = nullptr;
      ctx.dosage_cts[1] = nullptr;
    }
    const uint32_t domrec = !!(flags & (kfScoreDominant | kfScoreRecessive));
    uint32_t* common_geno_sum_incrs = nullptr;
    double* common_score_incrs = nullptr;
    uint32_t max_difflist_len = 0;
    uintptr_t raregeno_stride = 0;
    uintptr_t difflist_sample_ids_stride = 0;
    if ((!ctx.dosage_presents[0]) && (!domrec)) {
      max_difflist_len = (score_final_col_ct == 1)? (sample_ct / 16) : (sample_ct / 32);
      raregeno_stride = NypCtToAlignedWordCt(max_difflist_len);
      difflist_sample_ids_stride = RoundUpPow2(max_difflist_len, kInt32PerVec);
      if (unlikely(bigstack_alloc_u32(qsr_ct_nz, &common_geno_sum_incrs) ||
                   bigstack_alloc_d(score_final_col_ct, &common_score_incrs) ||
                   bigstack_alloc_w(kScoreVariantBlockSize * raregeno_stride, &ctx.raregenos[0]) ||
                   bigstack_alloc_w(kScoreVariantBlockSize * raregeno_stride, &ctx.raregenos[1]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.difflist_lens[0]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize, &ctx.difflist_lens[1]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize * difflist_sample_ids_stride, &ctx.difflist_sample_ids[0]) ||
                   bigstack_alloc_u32(kScoreVariantBlockSize * difflist_sample_ids_stride, &ctx.difflist_sample_ids[1]) ||
                   bigstack_alloc_i8(kScoreVariantBlockSize, &ctx.difflist_common_genos[0]) ||
                   bigstack_alloc_i8(kScoreVariantBlockSize, &ctx.difflist_common_genos[1]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize * score_final_col_ct, &ctx.score_sparse_coefs_vmaj[0]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize * score_final_col_ct, &ctx.score_sparse_coefs_vmaj[1]))) {
        goto ScoreReport_ret_NOMEM;
      }
    } else {
      ctx.raregenos[0] = nullptr;
      ctx.raregenos[1] = nullptr;
      ctx.difflist_lens[0] = nullptr;
      ctx.difflist_lens[1] = nullptr;
      ctx.difflist_sample_ids[0] = nullptr;
      ctx.difflist_sample_ids[1] = nullptr;
      ctx.difflist_common_genos[0] = nullptr;
      ctx.difflist_common_genos[1] = nullptr;
      ctx.score_sparse_coefs_vmaj[0] = nullptr;
      ctx.score_sparse_coefs_vmaj[1] = nullptr;
    }
    ctx.max_difflist_len = max_difflist_len;
    if (allele_freqs) {
      if (unlikely(bigstack_alloc_d(kScoreVariantBlockSize, &ctx.allele_freqs[0]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.allele_freqs[1]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.geno_slopes[0]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.geno_slopes[1]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.geno_intercepts[0]) ||
                   bigstack_alloc_d(kScoreVariantBlockSize, &ctx.geno_intercepts[1]))) {
        goto ScoreReport_ret_NOMEM;
      }
    } else {
      ctx.allele_freqs[0] = nullptr;
      ctx.allele_freqs[1] = nullptr;
      ctx.geno_slopes[0] = nullptr;
      ctx.geno_slopes[1] = nullptr;
      ctx.geno_intercepts[0] = nullptr;
      ctx.geno_intercepts[1] = nullptr;
    }
    SetThreadFuncAndData(CalcScoreThread, &ctx, &tg);

    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_nonmale_collapsed);
    AlignedBitarrInvert(sample_ct, sex_nonmale_collapsed);
    ctx.sex_nonmale_collapsed = sex_nonmale_collapsed;
    const uintptr_t* sex_nonfemale = sex_male;
    uintptr_t* sex_female_collapsed = sex_nonmale_collapsed;
    if (nosex_ct) {
      uintptr_t* sex_female_tmp;
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &sex_female_tmp) ||
                   bigstack_alloc_w(sample_ctl, &sex_female_collapsed))) {
        goto ScoreReport_ret_NOMEM;
      }
      BitvecInvmaskCopy(sex_nm, sex_male, raw_sample_ctl, sex_female_tmp);
      CopyBitarrSubset(sex_female_tmp, sample_include, sample_ct, sex_female_collapsed);
      // now invert to nonfemale
      BitvecInvert(raw_sample_ctl, sex_female_tmp);
      // (not actually necessary here, but let's keep sex_nonfemale definition
      // consistent since failure to do so has resulted in a bug)
      BitvecAnd(sample_include, raw_sample_ctl, sex_female_tmp);
      sex_nonfemale = sex_female_tmp;
    }
    ctx.sex_female_collapsed = sex_nonmale_collapsed;

    if (!variant_id_htable) {
      reterr = AllocAndPopulateIdHtableMt(variant_include, variant_ids, variant_ct, 0, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".sscore.tmp");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &score_tmpfile))) {
      goto ScoreReport_ret_OPEN_FAIL;
    }

    const uint32_t ignore_dup_ids = (flags / kfScoreIgnoreDupIds) & 1;
    const uint32_t list_variants = (flags / kfScoreListVariants) & 1;
    if (list_variants) {
      const uint32_t list_variants_zst = (flags / kfScoreListVariantsZs) & 1;
      OutnameZstSet(".sscore.vars", list_variants_zst, outname_end);
      reterr = InitCstream(outname, 0, list_variants_zst, 1, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
      cswritep = overflow_buf;
    }

    char** score_col_names;
    if (unlikely(bigstack_end_clalloc_cp(global_score_col_ct, &score_col_names))) {
      goto ScoreReport_ret_NOMEM;
    }
    if (!(flags & kfScoreHeaderRead)) {
      // just perform this check once, rather than once per file.
      // errored out with "Too many score columns" earlier if LHS > 2^31.
      if (global_score_col_ct * 16 > bigstack_left()) {
        goto ScoreReport_ret_NOMEM;
      }
    }

#ifdef USE_MTBLAS
    if (sample_shard_ct < calc_thread_ct) {
      const uint32_t matrix_multiply_thread_ct = 1 + ((calc_thread_ct - 1) / sample_shard_ct);
      BLAS_SET_NUM_THREADS(matrix_multiply_thread_ct);
    }
#endif

    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    const uint32_t lines_to_skip_p1 = 1 + ((flags / kfScoreHeaderIgnore) & 1);
    // Stricter threshold for multi-score since standard matrix-multiply
    // workflow is more likely to pay off
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t variance_standardize = (flags / kfScoreVarianceStandardize) & 1;
    const uint32_t center = variance_standardize || (flags & kfScoreCenter);
    const uint32_t no_meanimpute = (flags / kfScoreNoMeanimpute) & 1;
    const uint32_t y_prohibited = nosex_ct && (!no_meanimpute);
    const uint32_t se_mode = (flags / kfScoreSe) & 1;
    const uintptr_t tmpfile_buf_byte_ct = score_col_ct * sizeof(double);
    // If we're not attempting the sparse-genotype optimization, this is always
    // UINT32_MAX.
    uint32_t difflist_common_geno = UINT32_MAX;
    uint32_t cur_allele_ct = 2;
    uint32_t dosage_ct = 0;
    double geno_slope = kRecipDosageMax;
    double geno_intercept = 0.0;
    char* score_name_write_iter = R_CAST(char*, g_bigstack_end);
    for (uint32_t file_idx1 = 1; file_idx1 <= infile_ct; ++file_idx1, fname_iter = fname_iter->next) {
      cur_input_fname = fname_iter->str;
      if (file_idx1 > 1) {
        reterr = TextRetarget(cur_input_fname, &score_txs);
        if (unlikely(reterr)) {
          goto ScoreReport_ret_TSTREAM_FAIL;
        }
      }
      ZeroVecArr(11 * acc4_vec_ct * qsr_ct_nz, missing_diploid_accx);
      ZeroVecArr(11 * acc4_vec_ct * qsr_ct_nz, missing_haploid_accx);
      ZeroU32Arr(qsr_ct_nz, allele_ct_bases);
      ZeroU32Arr(qsr_ct_nz, male_allele_ct_decrs);
      ZeroU32Arr(qsr_ct_nz, nonfemale_allele_ct_incrs);
      ZeroU64Arr(qsr_ct_nz * sample_ct, ctx.ddosage_sums);
      ZeroWArr(raw_variant_ctl, already_seen_variants);
      ZeroWArr(raw_allele_ctl, already_seen_alleles);
      {
        uintptr_t cur_sample_shard_size = sample_shard_size;
        for (uint32_t shard_idx = 0; ; ++shard_idx) {
          if (shard_idx >= sample_shard_ct_m1) {
            if (shard_idx > sample_shard_ct_m1) {
              break;
            }
            cur_sample_shard_size = sample_ct - sample_shard_ct_m1 * sample_shard_size;
          }
          ZeroDArr(score_final_col_ct * cur_sample_shard_size, ctx.sharded_final_scores_cmaj[shard_idx]);
        }
      }
      if (common_geno_sum_incrs) {
        ZeroU32Arr(qsr_ct_nz, common_geno_sum_incrs);
        ZeroDArr(score_final_col_ct, common_score_incrs);
      }
      for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct_nz; ++qsr_idx) {
        variant_ct_rems[2 * qsr_idx] = 15;
        variant_ct_rems[2 * qsr_idx + 1] = 17;
        variant_hap_ct_rems[2 * qsr_idx] = 15;
        variant_hap_ct_rems[2 * qsr_idx + 1] = 17;
      }
      ctx.cur_variant_batch_size = kScoreVariantBlockSize;
      line_idx = 0;
      char* line_start;
      for (uint32_t uii = 0; uii != lines_to_skip_p1; ++uii) {
        ++line_idx;
        line_start = TextGet(&score_txs);
        if (unlikely(!line_start)) {
          if (!TextStreamErrcode2(&score_txs, &reterr)) {
            logerrprintf("Error: --score%s: %s is empty.\n", multi_input? "-list" : "", cur_input_fname);
            goto ScoreReport_ret_MALFORMED_INPUT;
          }
          goto ScoreReport_ret_TSTREAM_FAIL;
        }
      }
      {
        const uintptr_t global_score_col_idx_start = (file_idx1 - 1) * score_col_ct;
        if (flags & kfScoreHeaderRead) {
          char* line_end = TextLineEnd(&score_txs);
          if (S_CAST(uintptr_t, line_end - line_start) > bigstack_left()) {
            goto ScoreReport_ret_NOMEM;
          }
          char** cur_score_col_names = &(score_col_names[global_score_col_idx_start]);
          char* read_iter = line_start;
          for (uintptr_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
            read_iter = NextTokenMult0(read_iter, score_col_idx_deltas[score_col_idx]);
            if (unlikely(!read_iter)) {
              goto ScoreReport_ret_MISSING_TOKENS;
            }
            char* token_end = CurTokenEnd(read_iter);
            const uint32_t slen = token_end - read_iter;
            score_name_write_iter -= 1 + slen;
            cur_score_col_names[score_col_idx] = score_name_write_iter;
            memcpyx(score_name_write_iter, read_iter, slen, '\0');
          }
          if (file_idx1 == infile_ct) {
            BigstackEndSet(score_name_write_iter);
            bigstack_end_clalign();
            unsigned char* bigstack_mark2 = g_bigstack_base;
            uintptr_t* col_subset_mask;
            if (unlikely(bigstack_alloc_w(BitCtToWordCt(score_col_ct), &col_subset_mask))) {
              goto ScoreReport_ret_NOMEM;
            }
            SetAllBits(score_col_ct, col_subset_mask);
            uint32_t dup_found;
            reterr = CheckIdUniqueness(g_bigstack_base, g_bigstack_end, col_subset_mask, TO_CONSTCPCONSTP(score_col_names), score_col_ct, max_thread_ct, &dup_found);
            if (unlikely(reterr)) {
              goto ScoreReport_ret_1;
            }
            if (unlikely(dup_found)) {
              logerrprintf("Error: --score%s: Score IDs are not unique.\n", multi_input? "-list" : "");
              goto ScoreReport_ret_MALFORMED_INPUT;
            }
            BigstackReset(bigstack_mark2);
          }
        } else {
          const uintptr_t col_idx_stop = global_score_col_idx_start + score_col_ct;
          for (uintptr_t global_score_col_idx = global_score_col_idx_start; global_score_col_idx != col_idx_stop; ++global_score_col_idx) {
            const uint32_t str_blen = 6 + UintSlen(global_score_col_idx + 1);
            score_name_write_iter -= str_blen;
            score_col_names[global_score_col_idx] = score_name_write_iter;
            char* tmp_write_iter = strcpya_k(score_name_write_iter, "SCORE");
            u32toa_x(global_score_col_idx + 1, '\0', tmp_write_iter);
          }
          if (file_idx1 == infile_ct) {
            BigstackEndSet(score_name_write_iter);
            bigstack_end_clalign();
          }
        }
      }

      ReinitThreads(&tg);
      uint32_t block_vidx = 0;
      uint32_t sparse_vidx = 0;
      uint32_t parity = 0;
      uint32_t prev_variant_uidx = UINT32_MAX;
      uintptr_t* genovec_iter = ctx.genovecs[0];
      uintptr_t* raregeno_iter = ctx.raregenos[0];
      uint32_t* difflist_lens = ctx.difflist_lens[0];
      uint32_t* difflist_sample_ids_iter = ctx.difflist_sample_ids[0];
      int8_t* difflist_common_genos = ctx.difflist_common_genos[0];
      uintptr_t* dosage_present_iter = ctx.dosage_presents[0];
      Dosage* dosage_main_iter = ctx.dosage_mains[0];
      uint32_t* dosage_cts = ctx.dosage_cts[0];
      uintptr_t* missing_bitvec_iter = ctx.missing_bitvecs[0];
      uintptr_t* missing_male_bitvec_iter = ctx.missing_male_bitvecs[0];
      uintptr_t* missing_nonfemale_bitvec_iter = ctx.missing_nonfemale_bitvecs[0];
      double* cur_allele_freqs = ctx.allele_freqs[0];
      double* geno_slopes = ctx.geno_slopes[0];
      double* geno_intercepts = ctx.geno_intercepts[0];
      double* score_dense_coefs_cmaj = ctx.score_dense_coefs_cmaj[0];
      double* score_sparse_coefs_vmaj = ctx.score_sparse_coefs_vmaj[0];
      uint32_t* variant_uidxs = ctx.variant_uidxs[0];
      unsigned char* is_nonx_haploids = &(ctx.is_nonx_haploids[0][0]);
      unsigned char* is_relevant_xs = &(ctx.is_relevant_xs[0][0]);
      unsigned char* is_ys = &(ctx.is_ys[0][0]);
      unsigned char* ploidy_m1s = &(ctx.ploidy_m1s[0][0]);
      uint32_t valid_variant_ct = 0;
      uint32_t difflist_len = 0;
      uintptr_t missing_var_id_ct = 0;
      uintptr_t missing_allele_code_ct = 0;
      uintptr_t duplicated_var_id_ct = 0;
      if (flags & kfScoreHeaderRead) {
        ++line_idx;
        line_start = TextGet(&score_txs);
      }
      for (; line_start; ++line_idx, line_start = TextGet(&score_txs)) {
        // varid_col_idx and allele_col_idx will almost always be very small
        char* variant_id_start = NextTokenMult0(line_start, varid_col_idx);
        if (unlikely(!variant_id_start)) {
          goto ScoreReport_ret_MISSING_TOKENS;
        }
        char* variant_id_token_end = CurTokenEnd(variant_id_start);
        const uint32_t variant_id_slen = variant_id_token_end - variant_id_start;
        uint32_t variant_uidx = VariantIdDupflagHtableFind(variant_id_start, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
        if (variant_uidx >> 31) {
          ++missing_var_id_ct;
          if (variant_uidx != UINT32_MAX) {
            if (unlikely(!ignore_dup_ids)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --score%s variant ID '%s' appears multiple times in main dataset.\n", multi_input? "-list" : "", variant_ids[variant_uidx & 0x7fffffff]);
              goto ScoreReport_ret_INCONSISTENT_INPUT_WW;
            }
            ++duplicated_var_id_ct;
            // subtract this from missing_var_id_ct later
          }
          continue;
        }
        char* allele_start = NextTokenMult0(line_start, allele_col_idx);
        if (unlikely(!allele_start)) {
          goto ScoreReport_ret_MISSING_TOKENS;
        }
        uintptr_t allele_idx_offset_base;
        if (!allele_idx_offsets) {
          allele_idx_offset_base = variant_uidx * 2;
        } else {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        }
        char* allele_end = CurTokenEnd(allele_start);
        const char allele_end_char = *allele_end;
        *allele_end = '\0';
        const uint32_t allele_blen = 1 + S_CAST(uintptr_t, allele_end - allele_start);
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);

        uint32_t cur_aidx = 0;
        for (; cur_aidx != cur_allele_ct; ++cur_aidx) {
          if (memequal(allele_start, cur_alleles[cur_aidx], allele_blen)) {
            break;
          }
        }
        // compiler is smart enough to avoid repeating this test
        if (cur_aidx == cur_allele_ct) {
          ++missing_allele_code_ct;
          *allele_end = allele_end_char;
          continue;
        }
        const uintptr_t allele_idx = allele_idx_offset_base + cur_aidx;
        if (unlikely(IsSet(already_seen_alleles, allele_idx))) {
          char* errwrite_iter = strcpya_k(g_logbuf, "Error: --score");
          if (multi_input) {
            errwrite_iter = strcpya_k(errwrite_iter, "-list");
          }
          errwrite_iter = strcpya_k(errwrite_iter, ": ");
          // Don't write allele code, since it might be too long for the
          // buffer.
          if (!cur_aidx) {
            errwrite_iter = strcpya_k(errwrite_iter, "REF");
          } else {
            errwrite_iter = strcpya_k(errwrite_iter, "ALT");
            errwrite_iter = u32toa(cur_aidx, errwrite_iter);
          }
          errwrite_iter = strcpya_k(errwrite_iter, " allele for variant '");
          errwrite_iter = strcpya(errwrite_iter, variant_ids[variant_uidx]);
          errwrite_iter = strcpya_k(errwrite_iter, "' appears multiple times in ");
          errwrite_iter = strcpya(errwrite_iter, cur_input_fname);
          if (multi_input) {
            errwrite_iter = strcpya_k(errwrite_iter, "-list");
          }
          strcpy_k(errwrite_iter, " file.\n");
          goto ScoreReport_ret_MALFORMED_INPUT_WW;
        }
        SetBit(allele_idx, already_seen_alleles);
        const uint32_t is_new_variant = 1 - IsSet(already_seen_variants, variant_uidx);
        SetBit(variant_uidx, already_seen_variants);

        // okay, the variant and allele are in our dataset.
        const uint32_t chr_idx = GetVariantChr(cip, variant_uidx);
        uint32_t is_nonx_haploid = IsSet(cip->haploid_mask, chr_idx);
        if (unlikely(domrec && is_nonx_haploid)) {
          fputs("\n", stdout);
          logerrputs("Error: --score[-list] 'dominant' and 'recessive' modifiers cannot be used with\nhaploid chromosomes.\n");
          goto ScoreReport_ret_INCONSISTENT_INPUT;
        }
        uint32_t is_relevant_x = (chr_idx == x_code);
        if (unlikely(variance_standardize && (is_relevant_x || (chr_idx == mt_code)))) {
          fputs("\n", stdout);
          logerrputs("Error: --score[-list] 'variance-standardize' cannot be used with chrX or MT.\n");
          goto ScoreReport_ret_INCONSISTENT_INPUT;
        }
        is_nonx_haploid = (!is_relevant_x) && is_nonx_haploid;

        // only if --xchr-model 1 (which is no longer the default)
        is_relevant_x = is_relevant_x && xchr_model;

        const uint32_t is_y = (chr_idx == y_code);
        if (unlikely(is_y && y_prohibited)) {
          fputs("\n", stdout);
          logerrputs("Error: When both chrY variants and unknown-sex samples are present,\n--score[-list] can only be run with the 'no-mean-imputation' modifier.\n");
          goto ScoreReport_ret_INCONSISTENT_INPUT;
        }
        uint32_t ploidy_m1;
        if ((variant_uidx == prev_variant_uidx) && (cur_allele_ct == 2)) {
          // Previous record is the exact inverse of what we now want, so just
          // copy and invert it.
          // This is worth the trouble since PgrGet1D() is usually the main
          // bottleneck, and the repeated-variant case definitely comes up from
          // at least the recommended PCA-projection workflow.
          if (difflist_common_geno != UINT32_MAX) {
            memcpy(raregeno_iter, &(raregeno_iter[-S_CAST(intptr_t, raregeno_stride)]), NypCtToWordCt(difflist_len) * sizeof(intptr_t));
            memcpy(difflist_sample_ids_iter, &(difflist_sample_ids_iter[-S_CAST(intptr_t, difflist_sample_ids_stride)]), difflist_len * sizeof(int32_t));
            GenovecInvertUnsafe(difflist_len, raregeno_iter);
            difflist_common_geno = (6 - difflist_common_geno) & 3;
          } else {
            memcpy(genovec_iter, &(genovec_iter[-S_CAST(intptr_t, sample_ctaw2)]), sample_ctaw2 * sizeof(intptr_t));
            GenovecInvertUnsafe(sample_ct, genovec_iter);
            if (dosage_present_iter) {
              memcpy(dosage_present_iter, &(dosage_present_iter[-S_CAST(intptr_t, sample_ctaw)]), sample_ctaw * sizeof(intptr_t));
              dosage_ct = dosage_cts[block_vidx - 1];
              memcpy(dosage_main_iter, &(dosage_main_iter[-S_CAST(intptr_t, dosage_main_stride)]), dosage_ct * sizeof(Dosage));
              BiallelicDosage16Invert(dosage_ct, dosage_main_iter);
            }
          }
        } else if (dosage_present_iter || (cur_allele_ct > 2) || domrec || is_nonx_haploid || is_relevant_x) {
          difflist_common_geno = UINT32_MAX;
          reterr = PgrGet1D(sample_include, pssi, sample_ct, variant_uidx, cur_aidx, simple_pgrp, genovec_iter, dosage_present_iter, dosage_main_iter, &dosage_ct);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ScoreReport_ret_1;
          }
        } else {
          // Sparse-genotype optimization.  We limit this to the no-dosage,
          // all-diploid, biallelic, no 'dominant'/'recessive' case to keep
          // code size under control.
          reterr = PgrGetDifflistOrGenovec(sample_include, pssi, sample_ct, max_difflist_len, variant_uidx, simple_pgrp, genovec_iter, &difflist_common_geno, raregeno_iter, difflist_sample_ids_iter, &difflist_len);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ScoreReport_ret_1;
          }
          if (!cur_aidx) {
            if (difflist_common_geno == UINT32_MAX) {
              GenovecInvertUnsafe(sample_ct, genovec_iter);
            } else {
              GenovecInvertUnsafe(difflist_len, raregeno_iter);
              difflist_common_geno = (6 - difflist_common_geno) & 3;
            }
          }
        }
        if (difflist_common_geno != UINT32_MAX) {
          ZeroTrailingNyps(difflist_len, raregeno_iter);
          if (is_new_variant) {
            uint32_t missing_exists = (difflist_common_geno == 3);
            if ((!missing_exists) && difflist_len) {
              STD_ARRAY_DECL(uint32_t, 4, genocounts);
              GenoarrCountFreqsUnsafe(raregeno_iter, difflist_len, genocounts);
              missing_exists = (genocounts[3] != 0);
            }
            // This must be kept in sync with how the main branch processes
            // missingness.
            if (missing_exists) {
              SparseToMissingness(raregeno_iter, difflist_sample_ids_iter, sample_ct, difflist_common_geno, difflist_len, missing_bitvec_iter);
            }
            if (!qsr_ct) {
              allele_ct_bases[0] += 2;
              if (missing_exists) {
                VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, variant_ct_rems, missing_diploid_accx);
              }
            } else {
              const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
              for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                  allele_ct_bases[qsr_idx] += 2;
                  if (missing_exists) {
                    VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, &(variant_ct_rems[2 * qsr_idx]), &(missing_diploid_accx[acc4_vec_ct * 11 * qsr_idx]));
                  }
                }
              }
            }
          }
          // Enable sparse ddosage_sums update.
          if ((difflist_common_geno == 1) || (difflist_common_geno == 2)) {
            if (!qsr_ct) {
              common_geno_sum_incrs[0] += difflist_common_geno;
            } else {
              const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
              for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                  common_geno_sum_incrs[qsr_idx] += difflist_common_geno;
                }
              }
            }
          }
          ploidy_m1 = 1;
        } else {
          ZeroTrailingNyps(sample_ct, genovec_iter);
          GenoarrToMissingnessUnsafe(genovec_iter, sample_ct, missing_bitvec_iter);
          if (dosage_ct) {
            BitvecInvmask(dosage_present_iter, sample_ctl, missing_bitvec_iter);
          }
          if (is_nonx_haploid) {
            if (is_y) {
              if (is_new_variant) {
                if (!qsr_ct) {
                  nonfemale_allele_ct_incrs[0] += 1;
                } else {
                  const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                  for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                    if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                      nonfemale_allele_ct_incrs[qsr_idx] += 1;
                    }
                  }
                }
              }
              BitvecInvmask(sex_female_collapsed, sample_ctl, missing_bitvec_iter);
            } else if (is_new_variant) {
              if (!qsr_ct) {
                allele_ct_bases[0] += 1;
              } else {
                const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                  if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                    allele_ct_bases[qsr_idx] += 1;
                  }
                }
              }
            }
            if (is_new_variant) {
              if (!qsr_ct) {
                VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, variant_hap_ct_rems, missing_haploid_accx);
              } else {
                const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                  if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                    VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, &(variant_ct_rems[2 * qsr_idx]), &(missing_haploid_accx[acc4_vec_ct * 11 * qsr_idx]));
                  }
                }
              }
            }
            if (is_y) {
              memcpy(missing_nonfemale_bitvec_iter, missing_bitvec_iter, sample_ctl * sizeof(intptr_t));
              BitvecOr(sex_female_collapsed, sample_ctl, missing_bitvec_iter);
            }
            ploidy_m1 = 0;
          } else {
            // !is_nonx_haploid
            if (is_relevant_x) {
              // Handle nonmales first.  Remove males from missing_bitvec_iter,
              // to be restored at the end of this block.
              BitvecInvmaskCopy(missing_bitvec_iter, sex_nonmale_collapsed, sample_ctl, missing_male_bitvec_iter);
              BitvecAnd(sex_nonmale_collapsed, sample_ctl, missing_bitvec_iter);
            }
            if (is_new_variant) {
              if (!qsr_ct) {
                allele_ct_bases[0] += 2;
                VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, variant_ct_rems, missing_diploid_accx);
              } else {
                const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                  if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                    allele_ct_bases[qsr_idx] += 2;
                    VerticalCounterUpdate(missing_bitvec_iter, acc1_vec_ct, &(variant_ct_rems[2 * qsr_idx]), &(missing_diploid_accx[acc4_vec_ct * 11 * qsr_idx]));
                  }
                }
              }
            }
            if (is_relevant_x) {
              // Now count males.
              if (is_new_variant) {
                if (!qsr_ct) {
                  male_allele_ct_decrs[0] += 1;
                  VerticalCounterUpdate(missing_male_bitvec_iter, acc1_vec_ct, variant_hap_ct_rems, missing_haploid_accx);
                } else {
                  const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
                  for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                    if (IsSet(qsr_include, qsr_idx + bit_idx_base)) {
                      male_allele_ct_decrs[qsr_idx] += 1;
                      VerticalCounterUpdate(missing_male_bitvec_iter, acc1_vec_ct, &(variant_hap_ct_rems[2 * qsr_idx]), &(missing_haploid_accx[acc4_vec_ct * 11 * qsr_idx]));
                    }
                  }
                }
              }
              // Restore missing males to missing_bitvec_iter.
              BitvecOr(missing_male_bitvec_iter, sample_ctl, missing_bitvec_iter);
            }
            ploidy_m1 = domrec? 0 : 1;
          }
        }
        const double ploidy_d = ploidy_m1? 2.0 : 1.0;
        if (allele_freqs) {
          const double cur_allele_freq = GetAlleleFreq(&(allele_freqs[allele_idx_offset_base - variant_uidx]), cur_aidx, cur_allele_ct);
          if (center) {
            // Note that we prohibit 'dominant'/'recessive' from being used
            // with 'center'/'variance-standardize'.
            if (variance_standardize) {
              const double variance = ploidy_d * 0.5 * ComputeDiploidMultiallelicVariance(&(allele_freqs[allele_idx_offset_base - variant_uidx]), cur_allele_ct);
              if (!(variance > kSmallEpsilon)) {
                // ZeroTrailingNyps(sample_ct, genovec_buf);
                STD_ARRAY_DECL(uint32_t, 4, genocounts);
                if (difflist_common_geno != UINT32_MAX) {
                  GenoarrCountFreqsUnsafe(raregeno_iter, difflist_len, genocounts);
                  genocounts[difflist_common_geno] = sample_ct - difflist_len;
                } else {
                  GenoarrCountFreqsUnsafe(genovec_iter, sample_ct, genocounts);
                }
                if (unlikely(dosage_ct || genocounts[1] || genocounts[2])) {
                  snprintf(g_logbuf, kLogbufSize, "Error: --score[-list] variance-standardize failure for variant '%s': estimated allele frequency is zero or NaN, but not all dosages are zero. (This is possible when e.g. allele frequencies are estimated from founders, but the allele is only observed in nonfounders.)\n", variant_ids[variant_uidx]);
                  goto ScoreReport_ret_DEGENERATE_DATA_WW;
                }
                geno_slope = 0.0;
              } else {
                geno_slope = kRecipDosageMax / sqrt(variance);
              }
            }
            // (ploidy * cur_allele_freq * kDosageMax) * geno_slope +
            //   geno_intercept == 0
            // bugfix: must use "-1.0 *" instead of - to avoid unsigned int
            //   wraparound
            geno_intercept = (-1.0 * kDosageMax) * ploidy_d * cur_allele_freq * geno_slope;
          }
          cur_allele_freqs[block_vidx] = cur_allele_freq;
          geno_slopes[block_vidx] = geno_slope;
          geno_intercepts[block_vidx] = geno_intercept;
        }

        *allele_end = allele_end_char;
        const uint32_t dense_vidx = block_vidx - sparse_vidx;
        double* score_coefs_iter = (difflist_common_geno == UINT32_MAX)? (&(score_dense_coefs_cmaj[dense_vidx])) : (&(score_sparse_coefs_vmaj[sparse_vidx * score_final_col_ct]));
        double* common_score_incrs_iter = common_score_incrs;
        const char* read_iter = line_start;
        for (uint32_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
          read_iter = NextTokenMult0(read_iter, score_col_idx_deltas[score_col_idx]);
          if (unlikely(!read_iter)) {
            goto ScoreReport_ret_MISSING_TOKENS;
          }
          double raw_coef;
          const char* token_end = ScantokDouble(read_iter, &raw_coef);
          if (unlikely(!token_end)) {
            char* errwrite_iter = strcpya_k(g_logbuf, "Error: --score");
            if (multi_input) {
              errwrite_iter = strcpya_k(errwrite_iter, "-list");
            }
            errwrite_iter = strcpya_k(errwrite_iter, ": Invalid coefficient ");
            // don't want to write this unconditionally, since it can
            // theoretically overflow the buffer.
            token_end = CurTokenEnd(read_iter);
            const uint32_t token_slen = token_end - read_iter;
            if (token_slen <= 77) {
              *errwrite_iter++ = '\'';
              errwrite_iter = memcpya(errwrite_iter, read_iter, token_slen);
              errwrite_iter = strcpya_k(errwrite_iter, "' ");
            }
            errwrite_iter = strcpya_k(errwrite_iter, "on line ");
            errwrite_iter = wtoa(line_idx, errwrite_iter);
            errwrite_iter = strcpya_k(errwrite_iter, " of ");
            errwrite_iter = strcpya(errwrite_iter, cur_input_fname);
            strcpy_k(errwrite_iter, " .\n");
            goto ScoreReport_ret_MALFORMED_INPUT_WW;
          }
          if (difflist_common_geno == UINT32_MAX) {
            // dense case: fill column for later matrix-multiply
            if (!qsr_ct) {
              *score_coefs_iter = raw_coef;
              score_coefs_iter = &(score_coefs_iter[kScoreVariantBlockSize]);
            } else {
              const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
              for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                double cur_coef = raw_coef * u31tod(IsSet(qsr_include, qsr_idx + bit_idx_base));
                *score_coefs_iter = cur_coef;
                score_coefs_iter = &(score_coefs_iter[kScoreVariantBlockSize]);
              }
            }
          } else {
            // sparse case: use variant-major representation instead
            if (!qsr_ct) {
              score_coefs_iter[0] = raw_coef;
            } else {
              const uintptr_t bit_idx_base = RawToSubsettedPos(variant_include, variant_include_cumulative_popcounts, variant_uidx) * qsr_ct;
              for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct; ++qsr_idx) {
                const double cur_coef = raw_coef * u31tod(IsSet(qsr_include, qsr_idx + bit_idx_base));
                score_coefs_iter[qsr_idx] = cur_coef;
              }
            }
            double common_dosage = 0.0;
            if (difflist_common_geno != 3) {
              common_dosage = geno_intercept + kDosageMax * geno_slope * difflist_common_geno;
            } else if (!no_meanimpute) {
              common_dosage = kDosageMax * 2LL * cur_allele_freqs[block_vidx] * geno_slope;
            }
            if (!se_mode) {
              for (uint32_t uii = 0; uii != qsr_ct_nz; ++uii) {
                common_score_incrs_iter[uii] += common_dosage * score_coefs_iter[uii];
              }
            } else {
              for (uint32_t uii = 0; uii != qsr_ct_nz; ++uii) {
                score_coefs_iter[uii] *= score_coefs_iter[uii];
              }
              const double common_dosage_sq = common_dosage * common_dosage;
              for (uint32_t uii = 0; uii != qsr_ct_nz; ++uii) {
                common_score_incrs_iter[uii] += common_dosage_sq * score_coefs_iter[uii];
              }
            }
            score_coefs_iter = &(score_coefs_iter[qsr_ct_nz]);
            common_score_incrs_iter = &(common_score_incrs_iter[qsr_ct_nz]);
          }
          read_iter = token_end;
        }
        if (difflist_common_genos) {
          difflist_common_genos[block_vidx] = S_CAST(int8_t, difflist_common_geno);
          if (difflist_common_geno != UINT32_MAX) {
            difflist_lens[sparse_vidx] = difflist_len;
            raregeno_iter = &(raregeno_iter[raregeno_stride]);
            difflist_sample_ids_iter = &(difflist_sample_ids_iter[difflist_sample_ids_stride]);
            ++sparse_vidx;
          }
        }
        if (is_new_variant) {
          if (list_variants) {
            cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
            AppendBinaryEoln(&cswritep);
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto ScoreReport_ret_WRITE_FAIL;
            }
          }
          ++valid_variant_ct;
          if (!(valid_variant_ct % 10000)) {
            // Technically still processing a few of these variants, but this
            // is less misleading than the previous 'loaded' message.
            const uint32_t k_ct = valid_variant_ct / 1000;
            if (multi_input) {
              printf("\r--score-list file %u/%u: %uk variants processed.", file_idx1, infile_ct, k_ct);
            } else {
              printf("\r--score: %uk variants processed.", k_ct);
            }
            fflush(stdout);
          }
        }
        if (variant_uidxs) {
          variant_uidxs[block_vidx] = variant_uidx;
        }
        is_nonx_haploids[block_vidx] = is_nonx_haploid;
        is_relevant_xs[block_vidx] = is_relevant_x;
        is_ys[block_vidx] = is_y;
        ploidy_m1s[block_vidx] = ploidy_m1;
        prev_variant_uidx = variant_uidx;
        genovec_iter = &(genovec_iter[sample_ctaw2]);
        missing_bitvec_iter = &(missing_bitvec_iter[sample_ctaw]);
        missing_male_bitvec_iter = &(missing_male_bitvec_iter[sample_ctaw]);
        missing_nonfemale_bitvec_iter = &(missing_nonfemale_bitvec_iter[sample_ctaw]);
        if (dosage_present_iter) {
          dosage_present_iter = &(dosage_present_iter[sample_ctaw]);
          dosage_main_iter = &(dosage_main_iter[dosage_main_stride]);
          dosage_cts[block_vidx] = dosage_ct;
        }
        ++block_vidx;
        if (block_vidx == kScoreVariantBlockSize) {
          if (se_mode) {
            const uint32_t dense_ct = block_vidx - sparse_vidx;
            if (dense_ct) {
              for (uintptr_t score_final_col_idx = 0; score_final_col_idx != score_final_col_ct; ++score_final_col_idx) {
                double* score_dense_coefs_row = &(score_dense_coefs_cmaj[score_final_col_idx * kScoreVariantBlockSize]);
                for (uint32_t uii = 0; uii != dense_ct; ++uii) {
                  score_dense_coefs_row[uii] *= score_dense_coefs_row[uii];
                }
              }
            }
            for (uintptr_t ulii = 0; ulii != sparse_vidx * score_final_col_ct; ++ulii) {
              score_sparse_coefs_vmaj[ulii] *= score_sparse_coefs_vmaj[ulii];
            }
          }
          parity = 1 - parity;
          const uint32_t is_not_first_block = ThreadsAreActive(&tg);
          if (is_not_first_block) {
            JoinThreads(&tg);
            // CalcScoreThread() never errors out
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto ScoreReport_ret_THREAD_CREATE_FAIL;
          }
          // we could instead have prev_genovec, prev_dosage_present, etc.
          // pointers, but this should be good enough
          prev_variant_uidx = UINT32_MAX;

          variant_uidxs = ctx.variant_uidxs[parity];
          genovec_iter = ctx.genovecs[parity];
          raregeno_iter = ctx.raregenos[parity];
          difflist_lens = ctx.difflist_lens[parity];
          difflist_sample_ids_iter = ctx.difflist_sample_ids[parity];
          difflist_common_genos = ctx.difflist_common_genos[parity];
          dosage_present_iter = ctx.dosage_presents[parity];
          dosage_main_iter = ctx.dosage_mains[parity];
          dosage_cts = ctx.dosage_cts[parity];
          missing_bitvec_iter = ctx.missing_bitvecs[parity];
          missing_male_bitvec_iter = ctx.missing_male_bitvecs[parity];
          missing_nonfemale_bitvec_iter = ctx.missing_nonfemale_bitvecs[parity];
          cur_allele_freqs = ctx.allele_freqs[parity];
          geno_slopes = ctx.geno_slopes[parity];
          geno_intercepts = ctx.geno_intercepts[parity];
          score_dense_coefs_cmaj = ctx.score_dense_coefs_cmaj[parity];
          score_sparse_coefs_vmaj = ctx.score_sparse_coefs_vmaj[parity];
          is_nonx_haploids = &(ctx.is_nonx_haploids[parity][0]);
          is_relevant_xs = &(ctx.is_relevant_xs[parity][0]);
          is_ys = &(ctx.is_ys[parity][0]);
          ploidy_m1s = &(ctx.ploidy_m1s[parity][0]);
          block_vidx = 0;
          sparse_vidx = 0;
        }
      }
      if (unlikely(TextStreamErrcode2(&score_txs, &reterr))) {
        goto ScoreReport_ret_TSTREAM_FAIL;
      }
      for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct_nz; ++qsr_idx) {
        VecW* missing_diploid_acc4 = &(missing_diploid_accx[acc4_vec_ct * 11 * qsr_idx]);
        VecW* missing_diploid_acc8 = &(missing_diploid_acc4[acc4_vec_ct]);
        VecW* missing_diploid_acc32 = &(missing_diploid_acc8[acc8_vec_ct]);
        VecW* missing_haploid_acc4 = &(missing_haploid_accx[acc4_vec_ct * 11 * qsr_idx]);
        VecW* missing_haploid_acc8 = &(missing_haploid_acc4[acc4_vec_ct]);
        VecW* missing_haploid_acc32 = &(missing_haploid_acc8[acc8_vec_ct]);
        VcountIncr4To8(missing_diploid_acc4, acc4_vec_ct, missing_diploid_acc8);
        VcountIncr8To32(missing_diploid_acc8, acc8_vec_ct, missing_diploid_acc32);
        VcountIncr4To8(missing_haploid_acc4, acc4_vec_ct, missing_haploid_acc8);
        VcountIncr8To32(missing_haploid_acc8, acc8_vec_ct, missing_haploid_acc32);
      }
      const uint32_t is_not_first_block = ThreadsAreActive(&tg);
      putc_unlocked('\r', stdout);
      if (missing_var_id_ct || missing_allele_code_ct || duplicated_var_id_ct) {
        missing_var_id_ct -= duplicated_var_id_ct;
        if (!missing_var_id_ct) {
          if (missing_allele_code_ct) {
            if (missing_allele_code_ct == 1) {
              snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: 1 entry in %s was skipped due to a mismatching allele code.\n", multi_input? "-list" : "", cur_input_fname);
            } else {
              snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: %" PRIuPTR " entries in %s were skipped due to mismatching allele codes.\n", multi_input? "-list" : "", missing_allele_code_ct, cur_input_fname);
            }
          }
        } else if (!missing_allele_code_ct) {
          if (missing_var_id_ct == 1) {
            snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: 1 entry in %s was skipped due to a missing variant ID.\n", multi_input? "-list" : "", cur_input_fname);
          } else {
            snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: %" PRIuPTR " entries in %s were skipped due to missing variant IDs.\n", multi_input? "-list" : "", missing_var_id_ct, cur_input_fname);
          }
        } else {
          char* warning_write_iter = strcpya_k(g_logbuf, "Warning: --score");
          if (multi_input) {
            warning_write_iter = strcpya_k(warning_write_iter, "-list");
          }
          warning_write_iter = strcpya_k(warning_write_iter, ": ");
          if (missing_var_id_ct == 1) {
            warning_write_iter = strcpya_k(warning_write_iter, "1 entry in ");
            warning_write_iter = strcpya(warning_write_iter, cur_input_fname);
            warning_write_iter = strcpya_k(warning_write_iter," was skipped due to a missing variant ID");
          } else {
            warning_write_iter = i64toa(missing_var_id_ct, warning_write_iter);
            warning_write_iter = strcpya_k(warning_write_iter, " entries in ");
            warning_write_iter = strcpya(warning_write_iter, cur_input_fname);
            warning_write_iter = strcpya_k(warning_write_iter, " were skipped due to missing variant IDs");
          }
          warning_write_iter = strcpya_k(warning_write_iter, ", and ");
          if (missing_allele_code_ct == 1) {
            warning_write_iter = strcpya_k(warning_write_iter, "1 was skipped due to a mismatching allele code");
          } else {
            warning_write_iter = i64toa(missing_allele_code_ct, warning_write_iter);
            warning_write_iter = strcpya_k(warning_write_iter, " were skipped due to mismatching allele codes");
          }
          strcpy_k(warning_write_iter, ".\n");
        }
        WordWrapB(0);
        logerrputsb();
        if (duplicated_var_id_ct) {
          if (duplicated_var_id_ct == 1) {
            snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: 1 entry in %s was skipped since its variant ID appears multiple times in the main dataset.\n", multi_input? "-list" : "", cur_input_fname);
          } else {
            snprintf(g_logbuf, kLogbufSize, "Warning: --score%s: %" PRIuPTR " entries in %s were skipped since their variant IDs appear multiple times in the main dataset.\n", multi_input? "-list" : "", duplicated_var_id_ct, cur_input_fname);
          }
          WordWrapB(0);
          logerrputsb();
        }
        if (!list_variants) {
          logerrputs("(Add the 'list-variants' modifier to see which variants were actually used for\nscoring.)\n");
        }
      }
      if (block_vidx) {
        if (is_not_first_block) {
          JoinThreads(&tg);
        }
      } else if (unlikely(!valid_variant_ct)) {
        logerrprintf("Error: --score%s: No valid variants in %s.\n", multi_input? "-list" : "", cur_input_fname);
        goto ScoreReport_ret_DEGENERATE_DATA;
      } else {
        JoinThreads(&tg);
      }
      DeclareLastThreadBlock(&tg);
      ctx.cur_variant_batch_size = block_vidx;
      if (se_mode) {
        const uint32_t dense_ct = block_vidx - sparse_vidx;
        if (dense_ct) {
          for (uintptr_t score_final_col_idx = 0; score_final_col_idx != score_final_col_ct; ++score_final_col_idx) {
            double* score_dense_coefs_row = &(score_dense_coefs_cmaj[score_final_col_idx * kScoreVariantBlockSize]);
            for (uint32_t uii = 0; uii != dense_ct; ++uii) {
              score_dense_coefs_row[uii] *= score_dense_coefs_row[uii];
            }
          }
        }
        for (uintptr_t ulii = 0; ulii != sparse_vidx * score_final_col_ct; ++ulii) {
          score_sparse_coefs_vmaj[ulii] *= score_sparse_coefs_vmaj[ulii];
        }
      }
      if (unlikely(SpawnThreads(&tg))) {
        goto ScoreReport_ret_THREAD_CREATE_FAIL;
      }
      JoinThreads(&tg);
      // Sparse-optimization postprocessing.
      if (common_geno_sum_incrs) {
        uint64_t* ddosage_sums_iter = ctx.ddosage_sums;
        for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
          for (uint32_t qsr_idx = 0; qsr_idx != qsr_ct_nz; ++qsr_idx) {
            *ddosage_sums_iter += common_geno_sum_incrs[qsr_idx] * S_CAST(uint64_t, kDosageMax);
            ++ddosage_sums_iter;
          }
        }
        uintptr_t cur_sample_shard_size = sample_shard_size;
        for (uint32_t shard_idx = 0; ; ++shard_idx) {
          if (shard_idx >= sample_shard_ct_m1) {
            if (shard_idx > sample_shard_ct_m1) {
              break;
            }
            cur_sample_shard_size = sample_ct - sample_shard_ct_m1 * sample_shard_size;
          }
          double* shard_final_scores_cmaj_iter = ctx.sharded_final_scores_cmaj[shard_idx];
          for (uintptr_t ulii = 0; ulii != score_final_col_ct; ++ulii) {
            const double cur_score_incr = common_score_incrs[ulii];
            for (uint32_t shard_sample_idx = 0; shard_sample_idx != cur_sample_shard_size; ++shard_sample_idx) {
              *shard_final_scores_cmaj_iter += cur_score_incr;
              ++shard_final_scores_cmaj_iter;
            }
          }
        }
      }
      if (se_mode) {
        // sample_ct * score_final_col_ct
        uintptr_t cur_sample_shard_size = sample_shard_size;
        for (uint32_t shard_idx = 0; ; ++shard_idx) {
          if (shard_idx >= sample_shard_ct_m1) {
            if (shard_idx > sample_shard_ct_m1) {
              break;
            }
            cur_sample_shard_size = sample_ct - sample_shard_ct_m1 * sample_shard_size;
          }
          double* shard_final_scores_cmaj = ctx.sharded_final_scores_cmaj[shard_idx];
          for (uintptr_t ulii = 0; ulii != cur_sample_shard_size * score_final_col_ct; ++ulii) {
            shard_final_scores_cmaj[ulii] = sqrt(shard_final_scores_cmaj[ulii]);
          }
        }
      }
      // Now compute and spill scores to temporary file.
      for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct_nz; ++qsr_idx) {
        const uint32_t* scrambled_missing_diploid_cts = &(R_CAST(uint32_t*, missing_diploid_accx)[acc4_vec_ct * (3 + 11 * qsr_idx) * kInt32PerVec]);
        const uint32_t* scrambled_missing_haploid_cts = &(R_CAST(uint32_t*, missing_haploid_accx)[acc4_vec_ct * (3 + 11 * qsr_idx) * kInt32PerVec]);

        uint32_t cur_shard_start = 0;
        uintptr_t cur_sample_shard_size = sample_shard_size;
        uintptr_t sample_uidx_base = 0;
        uintptr_t sample_include_bits = sample_include[0];
        for (uint32_t shard_idx = 0; ; ++shard_idx) {
          if (shard_idx >= sample_shard_ct_m1) {
            if (shard_idx > sample_shard_ct_m1) {
              break;
            }
            cur_sample_shard_size = sample_ct - cur_shard_start;
          }
          const double* shard_final_scores_cmaj = ctx.sharded_final_scores_cmaj[shard_idx];
          for (uint32_t shard_sample_idx = 0; shard_sample_idx != cur_sample_shard_size; ++shard_sample_idx) {
            const uint32_t sample_idx = shard_sample_idx + cur_shard_start;
            const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
            const uint32_t scrambled_idx = VcountScramble1(sample_idx);
            uint32_t denom = allele_ct_bases[qsr_idx] + IsSet(sex_nonfemale, sample_uidx) * nonfemale_allele_ct_incrs[qsr_idx] - IsSet(sex_male, sample_uidx) * male_allele_ct_decrs[qsr_idx];
            if (no_meanimpute) {
              const uint32_t nallele = denom - 2 * scrambled_missing_diploid_cts[scrambled_idx] - scrambled_missing_haploid_cts[scrambled_idx];
              denom = nallele;
            }
            const double* final_score_col = &(shard_final_scores_cmaj[shard_sample_idx]);
            if (write_score_avgs) {
              const double denom_recip = 1.0 / S_CAST(double, denom);
              double* tmpfile_write_iter = tmpfile_buf;
              for (uintptr_t score_final_col_idx = qsr_idx; score_final_col_idx < score_final_col_ct; score_final_col_idx += qsr_ct_nz) {
                *tmpfile_write_iter++ = final_score_col[score_final_col_idx * cur_sample_shard_size] * denom_recip;
              }
              if (unlikely(fwrite_checked(tmpfile_buf, tmpfile_buf_byte_ct, score_tmpfile))) {
                goto ScoreReport_ret_WRITE_FAIL;
              }
            }
            if (write_score_sums) {
              double* tmpfile_write_iter = tmpfile_buf;
              for (uint32_t score_final_col_idx = qsr_idx; score_final_col_idx < score_final_col_ct; score_final_col_idx += qsr_ct_nz) {
                *tmpfile_write_iter++ = final_score_col[score_final_col_idx * cur_sample_shard_size];
              }
              if (unlikely(fwrite_checked(tmpfile_buf, tmpfile_buf_byte_ct, score_tmpfile))) {
                goto ScoreReport_ret_WRITE_FAIL;
              }
            }
          }
          cur_shard_start += cur_sample_shard_size;
        }
      }
      if (multi_input) {
        logprintf("--score-list file %u/%u: %u variant%s processed.\n", file_idx1, infile_ct, valid_variant_ct, (valid_variant_ct == 1)? "" : "s");
      } else {
        logprintf("--score: %u variant%s processed.\n", valid_variant_ct, (valid_variant_ct == 1)? "" : "s");
      }
    }
    if (unlikely(fclose_null(&score_tmpfile))) {
      goto ScoreReport_ret_WRITE_FAIL;
    }
    if (list_variants) {
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto ScoreReport_ret_WRITE_FAIL;
      }
      cswritep = nullptr;
      logprintf("--score: Variant list written to %s .\n", outname);
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".sscore.tmp");
    if (unlikely(fopen_checked(outname, FOPEN_RB, &score_tmpfile))) {
      goto ScoreReport_ret_OPEN_FAIL;
    }

    // major: file_idx
    // then qsr_idx (out of qsr_ct_nz)
    // then sample_idx (out of sample_ct)
    // then score_avgs vs. score_sums (out of write_type_ct)
    // then score_col_idx within a single input (out of score_col_ct doubles)
    const uint32_t write_type_ct = write_score_avgs + write_score_sums;
    const uint64_t sample_idx_stride = write_type_ct * score_col_ct;
    const uint64_t file_idx_stride = sample_idx_stride * sample_ct * qsr_ct_nz;
    const uint32_t* scrambled_missing_diploid_cts = nullptr;
    const uint32_t* scrambled_missing_haploid_cts = nullptr;
    double* saved_scores = nullptr;
    double* cur_saved_scores = tmpfile_buf;
    for (uintptr_t qsr_idx = 0; qsr_idx != qsr_ct_nz; ++qsr_idx) {
      BigstackReset(bigstack_mark);
      char* outname_end2 = outname_end;
      if (range_names) {
        *outname_end2++ = '.';
        outname_end2 = strcpya(outname_end2, range_names[qsr_idx]);
      }
      OutnameZstSet(".sscore", output_zst, outname_end2);
      reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
      cswritep = overflow_buf;

      if (!qsr_idx) {
        const uint64_t saved_scores_byte_ct = S_CAST(uint64_t, score_final_col_ct) * infile_ct * sample_ct * write_type_ct * sizeof(double);
        if (saved_scores_byte_ct <= bigstack_left()) {
          // Guaranteed to not exhaust memory, since there are no allocations
          // after this point, and if qsr_ct > 1 the subsequent InitCstream()
          // calls will allocate the same memory as the first one.
          saved_scores = S_CAST(double*, bigstack_end_alloc_raw_rd(saved_scores_byte_ct));
          if (unlikely(fread_checked(saved_scores, saved_scores_byte_ct, score_tmpfile))) {
            goto ScoreReport_ret_READ_TMP_FAIL;
          }
        }
      }
      const uint32_t write_fid = FidColIsRequired(siip, flags / kfScoreColMaybefid);
      const char* sample_ids = siip->sample_ids;
      const char* sids = siip->sids;
      const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
      const uintptr_t max_sid_blen = siip->max_sid_blen;
      const uint32_t write_sid = SidColIsRequired(sids, flags / kfScoreColMaybesid);
      const uint32_t write_empty_pheno = (flags & kfScoreColPheno1) && (!pheno_ct);
      const uint32_t write_phenos = (flags & (kfScoreColPheno1 | kfScoreColPhenos)) && pheno_ct;
      if (write_phenos && (!(flags & kfScoreColPhenos))) {
        pheno_ct = 1;
      }
      *cswritep++ = '#';
      if (write_fid) {
        cswritep = strcpya_k(cswritep, "FID\t");
      }
      cswritep = strcpya_k(cswritep, "IID");
      if (write_sid) {
        cswritep = strcpya_k(cswritep, "\tSID");
      }
      if (write_phenos) {
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, &(pheno_names[pheno_idx * max_pheno_name_blen]));
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto ScoreReport_ret_WRITE_FAIL;
          }
        }
      } else if (write_empty_pheno) {
        cswritep = strcpya_k(cswritep, "\tPHENO1");
      }
      const uint32_t write_nallele = (flags / kfScoreColNallele) & 1;
      if (write_nallele) {
        cswritep = strcpya_k(cswritep, "\tALLELE_CT");
      }
      const uint32_t write_denom = (flags / kfScoreColDenom) & 1;
      if (write_denom) {
        cswritep = strcpya_k(cswritep, "\tDENOM");
      }
      const uint32_t write_dosage_sum = (flags / kfScoreColDosageSum) & 1;
      if (write_dosage_sum) {
        cswritep = strcpya_k(cswritep, "\tNAMED_ALLELE_DOSAGE_SUM");
      }
      if (write_score_avgs) {
        for (uint32_t global_score_col_idx = 0; global_score_col_idx != global_score_col_ct; ++global_score_col_idx) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, score_col_names[global_score_col_idx]);
          cswritep = strcpya_k(cswritep, "_AVG");
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto ScoreReport_ret_WRITE_FAIL;
          }
        }
      }
      if (write_score_sums) {
        for (uint32_t global_score_col_idx = 0; global_score_col_idx != global_score_col_ct; ++global_score_col_idx) {
          *cswritep++ = '\t';
          cswritep = strcpya(cswritep, score_col_names[global_score_col_idx]);
          cswritep = strcpya_k(cswritep, "_SUM");
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto ScoreReport_ret_WRITE_FAIL;
          }
        }
      }
      AppendBinaryEoln(&cswritep);
      if (!multi_input) {
        scrambled_missing_diploid_cts = &(R_CAST(uint32_t*, missing_diploid_accx)[acc4_vec_ct * (3 + 11 * qsr_idx) * kInt32PerVec]);
        scrambled_missing_haploid_cts = &(R_CAST(uint32_t*, missing_haploid_accx)[acc4_vec_ct * (3 + 11 * qsr_idx) * kInt32PerVec]);
      }
      const uint32_t omp_slen = strlen(output_missing_pheno);

      uintptr_t sample_uidx_base = 0;
      uintptr_t sample_include_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
        cswritep = AppendXid(sample_ids, sids, write_fid, write_sid, max_sample_id_blen, max_sid_blen, sample_uidx, cswritep);
        if (write_phenos) {
          // er, this probably belongs in its own function
          for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
            const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
            const PhenoDtype type_code = cur_pheno_col->type_code;
            *cswritep++ = '\t';
            if (type_code <= kPhenoDtypeQt) {
              if (!IsSet(cur_pheno_col->nonmiss, sample_uidx)) {
                cswritep = memcpya(cswritep, output_missing_pheno, omp_slen);
              } else if (type_code == kPhenoDtypeCc) {
                *cswritep++ = '1' + IsSet(cur_pheno_col->data.cc, sample_uidx);
              } else {
                cswritep = dtoa_g(cur_pheno_col->data.qt[sample_uidx], cswritep);
              }
            } else {
              // category index guaranteed to be zero for missing values
              cswritep = strcpya(cswritep, cur_pheno_col->category_names[cur_pheno_col->data.cat[sample_uidx]]);
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto ScoreReport_ret_WRITE_FAIL;
              }
            }
          }
        } else if (write_empty_pheno) {
          *cswritep++ = '\t';
          cswritep = memcpya(cswritep, output_missing_pheno, omp_slen);
        }

        if (!multi_input) {
          const uint32_t scrambled_idx = VcountScramble1(sample_idx);
          uint32_t denom = allele_ct_bases[qsr_idx] + IsSet(sex_nonfemale, sample_uidx) * nonfemale_allele_ct_incrs[qsr_idx] - IsSet(sex_male, sample_uidx) * male_allele_ct_decrs[qsr_idx];
          const uint32_t nallele = denom - 2 * scrambled_missing_diploid_cts[scrambled_idx] - scrambled_missing_haploid_cts[scrambled_idx];
          if (write_nallele) {
            *cswritep++ = '\t';
            cswritep = u32toa(nallele, cswritep);
          }
          if (no_meanimpute) {
            denom = nallele;
          }
          if (write_denom) {
            *cswritep++ = '\t';
            cswritep = u32toa(denom, cswritep);
          }
          if (write_dosage_sum) {
            *cswritep++ = '\t';
            cswritep = ddosagetoa(ctx.ddosage_sums[sample_idx + qsr_idx * sample_ct], cswritep);
          }
        }
        uint64_t qsr_sample_idx_offset = (qsr_idx * sample_ct + sample_idx) * sample_idx_stride;
        for (uint32_t write_type_idx = 0; write_type_idx != write_type_ct; ++write_type_idx) {
          if (write_type_idx == 1) {
            qsr_sample_idx_offset += score_col_ct;
          }
          for (uintptr_t file_idx = 0; file_idx != infile_ct; ++file_idx) {
            const uint64_t tmp_idx = file_idx * file_idx_stride + qsr_sample_idx_offset;
            if (saved_scores) {
              cur_saved_scores = &(saved_scores[tmp_idx]);
            } else {
              if (unlikely(fseeko(score_tmpfile, tmp_idx * sizeof(double), SEEK_SET) ||
                           fread_checked(cur_saved_scores, tmpfile_buf_byte_ct, score_tmpfile))) {
                goto ScoreReport_ret_READ_TMP_FAIL;
              }
            }
            for (uintptr_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
              *cswritep++ = '\t';
              cswritep = dtoa_g(cur_saved_scores[score_col_idx], cswritep);
            }
          }
        }
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto ScoreReport_ret_WRITE_FAIL;
        }
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto ScoreReport_ret_WRITE_FAIL;
      }
    }
    if (!qsr_ct) {
      logprintfww("--score%s: Results written to %s .\n", multi_input? "-list" : "", outname);
    } else {
      *outname_end = '\0';
      logprintfww("--score%s + --q-score-range: Results written to %s.<range name>.sscore%s .\n", multi_input? "-list" : "", outname, output_zst? ".zst" : "");
    }
    if (unlikely(fclose_null(&score_tmpfile))) {
      goto ScoreReport_ret_READ_TMP_FAIL;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".sscore.tmp");
    if (unlink(outname)) {
      logerrprintfww("Error: Failed to delete %s : %s.\n", outname, strerror(errno));
      goto ScoreReport_ret_WRITE_FAIL;
    }
  }
  while (0) {
  ScoreReport_ret_TSTREAM_FAIL:
    TextStreamErrPrint(cur_input_fname, &score_txs);
    break;
  ScoreReport_ret_QSR_RANGE_TSTREAM_FAIL:
    TextStreamErrPrint("--q-score-range range file", &score_txs);
    break;
  ScoreReport_ret_QSR_DATA_TSTREAM_FAIL:
    TextStreamErrPrint("--q-score-range data file", &score_txs);
    break;
  ScoreReport_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ScoreReport_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ScoreReport_ret_READ_TMP_FAIL:
    logerrprintfww(kErrprintfFread, "--score[-list] temporary file", rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  ScoreReport_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ScoreReport_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  ScoreReport_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    fputs("\n", stdout);
    logerrputsb();
  ScoreReport_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ScoreReport_ret_QSR_DATA_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of --q-score-range data file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetInconsistentInput;
    break;
  ScoreReport_ret_MISSING_TOKENS:
    fputs("\n", stdout);
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, cur_input_fname);
    reterr = kPglRetInconsistentInput;
    break;
  ScoreReport_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    fputs("\n", stdout);
    logerrputsb();
  ScoreReport_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ScoreReport_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  ScoreReport_ret_DEGENERATE_DATA_WW:
    WordWrapB(0);
    fputs("\n", stdout);
    logerrputsb();
  ScoreReport_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  ScoreReport_ret_TOO_MANY_COLUMNS:
    logerrprintf("Error: --score%s: Too many score columns.\n", (score_info_ptr->flags & kfScoreMultiInput)? "-list" : "");
    reterr = kPglRetNotYetSupported;
    break;
  }
 ScoreReport_ret_1:
  fclose_cond(score_tmpfile);
  CswriteCloseCond(&css, cswritep);
  CleanupThreads(&tg);
  BLAS_SET_NUM_THREADS(1);
  CleanupTextStream2("--score[-list] file", &score_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

typedef struct VscoreCtxStruct {
  const uintptr_t* variant_include;
  const ChrInfo* cip;
  const uintptr_t* allele_idx_offsets;
  const double* allele_freqs;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed;
  const uintptr_t* sex_male_interleaved_vec;
  const float* wts_f_smaj;
  const double* wts_d_smaj;
  uint32_t vscore_ct;
  uint32_t sample_ct;
  uint32_t male_ct;
  uint32_t is_xchr_model_1;
  uint32_t max_difflist_len;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** raregenos;
  uint32_t** difflist_sample_id_bufs;
  uintptr_t** dosage_presents;
  Dosage** dosage_mains;
  uint32_t* read_variant_uidx_starts;

  uint32_t cur_block_size;

  float** dosage_f_vmaj_bufs;
  float** tmp_f_result_bufs;
  double** dosage_d_vmaj_bufs;
  double** tmp_d_result_bufs;
#ifdef USE_CUDA
  CublasFmultiplier* cfms;
#endif

  // variant-major
  float* results_f[2];
  double* results_d[2];

  uint32_t* missing_cts[2];

  // high 32 bits = unused for now
  // low 32 bits = uint32_t(PglErr)
  uint64_t err_info;
} VscoreCtx;

// This setting seems optimal on my Mac (smaller doesn't take full advantage of
// AVX, larger creates cache problems?).
CONSTI32(kVscoreBlockSize, 32);

THREAD_FUNC_DECL VscoreThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  VscoreCtx* ctx = S_CAST(VscoreCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const ChrInfo* cip = ctx->cip;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const double* allele_freqs = ctx->allele_freqs;
  const uintptr_t* sample_include = ctx->sample_include;
  const uintptr_t* sex_male = ctx->sex_male_collapsed;
  const uintptr_t* sex_male_interleaved_vec = ctx->sex_male_interleaved_vec;
  const float* wts_f_smaj = ctx->wts_f_smaj;
  const double* wts_d_smaj = ctx->wts_d_smaj;

  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  uintptr_t* genovec = ctx->genovecs[tidx];
  uintptr_t* raregeno = ctx->raregenos[tidx];
  uint32_t* difflist_sample_ids = ctx->difflist_sample_id_bufs[tidx];
  uintptr_t* dosage_present = nullptr;
  Dosage* dosage_main = nullptr;
  if (ctx->dosage_presents) {
    dosage_present = ctx->dosage_presents[tidx];
    dosage_main = ctx->dosage_mains[tidx];
  }

  const uintptr_t vscore_ct = ctx->vscore_ct;
  const uintptr_t sample_ct = ctx->sample_ct;
  const uint32_t male_ct = ctx->male_ct;
  const uint32_t nonmale_ct = sample_ct - male_ct;
  const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
  const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
  const uint32_t is_xchr_model_1 = ctx->is_xchr_model_1;
  const uint32_t max_difflist_len = ctx->max_difflist_len;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);

  const uint32_t single_prec = (wts_f_smaj != nullptr);
#ifdef USE_CUDA
  CublasFmultiplier* cfmp = ctx->cfms? &(ctx->cfms[tidx]) : nullptr;
#endif

  float* tmp_f_result_buf = nullptr;
  double* tmp_d_result_buf = nullptr;
  uint16_t cur_bidxs[kVscoreBlockSize];

  float* dosage_f_vmaj = nullptr;
  double* dosage_d_vmaj = nullptr;
  if (single_prec) {
    tmp_f_result_buf = ctx->tmp_f_result_bufs[tidx];
    dosage_f_vmaj = ctx->dosage_f_vmaj_bufs[tidx];
#ifdef USE_CUDA
    if (CudaSetDevice(cfmp->device_idx)) {
      fputs("unexpected CudaSetDevice() failure in VscoreThread()\n", stderr);
      exit(S_CAST(int, kPglRetGpuFail));
    }
#endif
  } else {
    tmp_d_result_buf = ctx->tmp_d_result_bufs[tidx];
    dosage_d_vmaj = ctx->dosage_d_vmaj_bufs[tidx];
  }

  uint32_t is_y = 0;
  uint32_t is_x_or_y = 0;
  uint32_t is_nonxy_haploid = 0;
  uint32_t chr_end = 0;
  double slope = 0.0;

  uint32_t dosage_ct = 0;

  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_size = ctx->cur_block_size;
    const uint32_t bidx_end = ((tidx + 1) * cur_block_size) / calc_thread_ct;
    float* cur_results_f = ctx->results_f[parity];
    double* cur_results_d = ctx->results_d[parity];
    uint32_t* missing_cts = ctx->missing_cts[parity];
    uintptr_t row_idx = 0;
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);
    for (uint32_t variant_bidx = (tidx * cur_block_size) / calc_thread_ct; variant_bidx != bidx_end; ++variant_bidx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      if (variant_uidx >= chr_end) {
        const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        is_y = 0;
        is_nonxy_haploid = 0;
        if (chr_idx == x_code) {
          is_x_or_y = is_xchr_model_1;
        } else if (chr_idx == y_code) {
          is_x_or_y = 1;
          is_y = 1;
        } else {
          is_x_or_y = 0;
          is_nonxy_haploid = IsSet(cip->haploid_mask, chr_idx);
        }
        slope = (is_nonxy_haploid || is_y)? 0.5 : 1.0;
      }
      double ref_freq;
      if (!allele_idx_offsets) {
        ref_freq = allele_freqs[variant_uidx];
      } else {
        ref_freq = allele_freqs[allele_idx_offsets[variant_uidx] - variant_uidx];
      }
      const double missing_val = slope * 2 * (1.0 - ref_freq);
      if (!dosage_present) {
        uint32_t difflist_common_geno;
        uint32_t difflist_len;
        PglErr reterr = PgrGetDifflistOrGenovec(sample_include, pssi, sample_ct, max_difflist_len, variant_uidx, pgrp, genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto VscoreThread_err;
        }
        if (difflist_common_geno != UINT32_MAX) {
          if ((!is_x_or_y) && (!difflist_common_geno)) {

            float* target_f = nullptr;
            double* target_d = nullptr;
            if (single_prec) {
              target_f = &(cur_results_f[variant_bidx * vscore_ct]);
            } else {
              target_d = &(cur_results_d[variant_bidx * vscore_ct]);
            }
            uint32_t missing_ct = 0;
            if (!difflist_len) {
              if (single_prec) {
                ZeroFArr(vscore_ct, target_f);
              } else {
                ZeroDArr(vscore_ct, target_d);
              }
            } else {
              ZeroTrailingNyps(difflist_len, raregeno);
              if (single_prec) {
                ZeroFArr(vscore_ct * 3, tmp_f_result_buf);
              } else {
                ZeroDArr(vscore_ct * 3, tmp_d_result_buf);
              }
              const uint32_t word_ct_m1 = (difflist_len - 1) / kBitsPerWordD2;
              uint32_t loop_len = kBitsPerWordD2;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= word_ct_m1) {
                  if (widx > word_ct_m1) {
                    break;
                  }
                  loop_len = ModNz(difflist_len, kBitsPerWordD2);
                }
                // slightly nicer to work with 2..0 than 1..3 row-indexes
                uintptr_t raregeno_word = raregeno[widx];
                uintptr_t raregeno_invword = ~raregeno_word;
                missing_ct += Popcount01Word(raregeno_word & (raregeno_word >> 1) & kMask5555);
                const uint32_t* cur_difflist_sample_ids = &(difflist_sample_ids[widx * kBitsPerWordD2]);
                if (single_prec) {
                  for (uint32_t uii = 0; uii != loop_len; ++uii) {
                    const uintptr_t sample_idx = cur_difflist_sample_ids[uii];
                    const uint32_t cur_invgeno = raregeno_invword & 3;
                    const float* incr_src = &(wts_f_smaj[sample_idx * vscore_ct]);
                    float* incr_dst = &(tmp_f_result_buf[cur_invgeno * vscore_ct]);
                    for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                      incr_dst[ulii] += incr_src[ulii];
                    }
                    raregeno_invword = raregeno_invword >> 2;
                  }
                } else {
                  for (uint32_t uii = 0; uii != loop_len; ++uii) {
                    const uintptr_t sample_idx = cur_difflist_sample_ids[uii];
                    const uint32_t cur_invgeno = raregeno_invword & 3;
                    const double* incr_src = &(wts_d_smaj[sample_idx * vscore_ct]);
                    double* incr_dst = &(tmp_d_result_buf[cur_invgeno * vscore_ct]);
                    for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                      incr_dst[ulii] += incr_src[ulii];
                    }
                    raregeno_invword = raregeno_invword >> 2;
                  }
                }
              }
              if (single_prec) {
                if (!is_nonxy_haploid) {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_f[ulii] = 2 * tmp_f_result_buf[ulii + vscore_ct] + tmp_f_result_buf[ulii + 2 * vscore_ct];
                  }
                } else {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_f[ulii] = tmp_f_result_buf[ulii + vscore_ct] + S_CAST(float, 0.5) * tmp_f_result_buf[ulii + 2 * vscore_ct];
                  }
                }
                if (missing_ct) {
                  const float missing_val_f = S_CAST(float, missing_val);
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_f[ulii] += missing_val_f * tmp_f_result_buf[ulii];
                  }
                }
              } else {
                if (!is_nonxy_haploid) {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_d[ulii] = 2 * tmp_d_result_buf[ulii + vscore_ct] + tmp_d_result_buf[ulii + 2 * vscore_ct];
                  }
                } else {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_d[ulii] = tmp_d_result_buf[ulii + vscore_ct] + 0.5 * tmp_d_result_buf[ulii + 2 * vscore_ct];
                  }
                }
                if (missing_ct) {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_d[ulii] += missing_val * tmp_d_result_buf[ulii];
                  }
                }
              }
            }
            if (missing_cts) {
              missing_cts[variant_bidx] = missing_ct;
            }
            continue;
          }
          PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, sample_ct, difflist_len, genovec);
        }
      } else {
        PglErr reterr = PgrGetD(sample_include, pssi, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
        if (unlikely(reterr)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
          goto VscoreThread_err;
        }
        if ((!is_x_or_y) && (dosage_ct <= max_difflist_len)) {
          STD_ARRAY_DECL(uint32_t, 4, genocounts);
          ZeroTrailingNyps(sample_ct, genovec);
          if (!dosage_ct) {
            // dosage_present contains garbage if dosage_ct == 0; might want to
            // append 'Unsafe' to PgrGetD and similar function names...
            ZeroWArr(BitCtToWordCt(sample_ct), dosage_present);
          }
          GenoarrCountInvsubsetFreqs2(genovec, dosage_present, sample_ct, sample_ct - dosage_ct, genocounts);
          if (genocounts[0] >= sample_ct - max_difflist_len) {
            float* target_f = nullptr;
            double* target_d = nullptr;
            if (single_prec) {
              target_f = &(cur_results_f[variant_bidx * vscore_ct]);
            } else {
              target_d = &(cur_results_d[variant_bidx * vscore_ct]);
            }
            if (genocounts[0] == sample_ct) {
              if (single_prec) {
                ZeroFArr(vscore_ct, target_f);
              } else {
                ZeroDArr(vscore_ct, target_d);
              }
            } else {
              if (single_prec) {
                ZeroFArr(vscore_ct * 3, tmp_f_result_buf);
              } else {
                ZeroDArr(vscore_ct * 3, tmp_d_result_buf);
              }
              const Halfword* dosage_present_alias = R_CAST(Halfword*, dosage_present);
              const uint32_t sample_ctl2 = DivUp(sample_ct, kBitsPerWordD2);
              for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
                uintptr_t geno_word = genovec[widx];
                if (!geno_word) {
                  continue;
                }
                const uintptr_t dosage_mask = UnpackHalfwordToWord(dosage_present_alias[widx]);
                geno_word = geno_word & (~(dosage_mask * 3));
                if (!geno_word) {
                  continue;
                }
                if (single_prec) {
                  const float* cur_wts_smaj = &(wts_f_smaj[widx * kBitsPerWordD2 * vscore_ct]);
                  do {
                    const uint32_t shift_ct = ctzw(geno_word) & (~1);
                    const uintptr_t cur_invgeno = 3 & (~(geno_word >> shift_ct));
                    const float* incr_src = &(cur_wts_smaj[(shift_ct / 2) * vscore_ct]);
                    float* incr_dst = &(tmp_f_result_buf[cur_invgeno * vscore_ct]);
                    for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                      incr_dst[ulii] += incr_src[ulii];
                    }
                    geno_word &= ~((3 * k1LU) << shift_ct);
                  } while (geno_word);
                } else {
                  const double* cur_wts_smaj = &(wts_d_smaj[widx * kBitsPerWordD2 * vscore_ct]);
                  do {
                    const uint32_t shift_ct = ctzw(geno_word) & (~1);
                    const uintptr_t cur_invgeno = 3 & (~(geno_word >> shift_ct));
                    const double* incr_src = &(cur_wts_smaj[(shift_ct / 2) * vscore_ct]);
                    double* incr_dst = &(tmp_d_result_buf[cur_invgeno * vscore_ct]);
                    for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                      incr_dst[ulii] += incr_src[ulii];
                    }
                    geno_word &= ~((3 * k1LU) << shift_ct);
                  } while (geno_word);
                }
              }
              if (single_prec) {
                if (!is_nonxy_haploid) {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_f[ulii] = 2 * tmp_f_result_buf[ulii + vscore_ct] + tmp_f_result_buf[ulii + 2 * vscore_ct];
                  }
                } else {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_f[ulii] = tmp_f_result_buf[ulii + vscore_ct] + S_CAST(float, 0.5) * tmp_f_result_buf[ulii + 2 * vscore_ct];
                  }
                }
                if (genocounts[3]) {
                  const float missing_val_f = S_CAST(float, missing_val);
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_f[ulii] += missing_val_f * tmp_f_result_buf[ulii];
                  }
                }
              } else {
                if (!is_nonxy_haploid) {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_d[ulii] = 2 * tmp_d_result_buf[ulii + vscore_ct] + tmp_d_result_buf[ulii + 2 * vscore_ct];
                  }
                } else {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_d[ulii] = tmp_d_result_buf[ulii + vscore_ct] + 0.5 * tmp_d_result_buf[ulii + 2 * vscore_ct];
                  }
                }
                if (genocounts[3]) {
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_d[ulii] += missing_val * tmp_d_result_buf[ulii];
                  }
                }
              }
              uintptr_t sample_idx_base = 0;
              uintptr_t dosage_present_bits = dosage_present[0];
              if (single_prec) {
                const float slope_f = S_CAST(float, slope);
                for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
                  const uintptr_t sample_idx = BitIter1(dosage_present, &sample_idx_base, &dosage_present_bits);
                  const float* incr_src = &(wts_f_smaj[sample_idx * vscore_ct]);
                  const float cur_dosage = slope_f * S_CAST(float, kRecipDosageMid) * u31tof(dosage_main[dosage_idx]);
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_f[ulii] += cur_dosage * incr_src[ulii];
                  }
                }
              } else {
                for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
                  const uintptr_t sample_idx = BitIter1(dosage_present, &sample_idx_base, &dosage_present_bits);
                  const double* incr_src = &(wts_d_smaj[sample_idx * vscore_ct]);
                  const double cur_dosage = slope * kRecipDosageMid * u31tod(dosage_main[dosage_idx]);
                  for (uintptr_t ulii = 0; ulii != vscore_ct; ++ulii) {
                    target_d[ulii] += cur_dosage * incr_src[ulii];
                  }
                }
              }
            }
            if (missing_cts) {
              missing_cts[variant_bidx] = genocounts[3];
            }
            continue;
          }
        }
      }

      if (row_idx == kVscoreBlockSize) {
        if (single_prec) {
#ifdef USE_CUDA
          if (cfmp) {
            if (unlikely(CublasFmultiplyRowMajor1(dosage_f_vmaj, cfmp, tmp_f_result_buf))) {
              new_err_info = S_CAST(uint32_t, kPglRetGpuFail);
              goto VscoreThread_err;
            }
          } else {
            RowMajorFmatrixMultiply(dosage_f_vmaj, wts_f_smaj, kVscoreBlockSize, vscore_ct, sample_ct, tmp_f_result_buf);
          }
#else
          RowMajorFmatrixMultiply(dosage_f_vmaj, wts_f_smaj, kVscoreBlockSize, vscore_ct, sample_ct, tmp_f_result_buf);
#endif
          const float* tmp_f_result_iter = tmp_f_result_buf;
          for (uintptr_t ulii = 0; ulii != kVscoreBlockSize; ++ulii) {
            const uintptr_t cur_bidx = cur_bidxs[ulii];
            memcpy(&(cur_results_f[cur_bidx * vscore_ct]), tmp_f_result_iter, vscore_ct * sizeof(float));
            tmp_f_result_iter = &(tmp_f_result_iter[vscore_ct]);
          }
        } else {
          RowMajorMatrixMultiply(dosage_d_vmaj, wts_d_smaj, kVscoreBlockSize, vscore_ct, sample_ct, tmp_d_result_buf);
          const double* tmp_d_result_iter = tmp_d_result_buf;
          for (uintptr_t ulii = 0; ulii != kVscoreBlockSize; ++ulii) {
            const uintptr_t cur_bidx = cur_bidxs[ulii];
            memcpy(&(cur_results_d[cur_bidx * vscore_ct]), tmp_d_result_iter, vscore_ct * sizeof(double));
            tmp_d_result_iter = &(tmp_d_result_iter[vscore_ct]);
          }
        }
        row_idx = 0;
      }
      cur_bidxs[row_idx] = variant_bidx;
      if (single_prec) {
        float* cur_row = &(dosage_f_vmaj[row_idx * sample_ct]);
        PopulateRescaledDosageF(genovec, dosage_present, dosage_main, S_CAST(float, slope), S_CAST(float, 0.0), S_CAST(float, missing_val), sample_ct, dosage_ct, cur_row);
        if (is_x_or_y) {
          // Instead of doing this for every variant, we could precompute
          // chrX/chrY weight matrices with male weights halved/nonmale weights
          // zeroed out.  But the number of chrY variants is typically small
          // enough (and how often will --xchr-model 1 be used, anyway?) that I
          // don't think it's worth it.
          uintptr_t sample_uidx_base = 0;
          if (is_y) {
            // zero out nonmale values
            uintptr_t sex_male_invbits = ~sex_male[0];
            for (uint32_t nonmale_idx = 0; nonmale_idx != nonmale_ct; ++nonmale_idx) {
              const uintptr_t sample_uidx = BitIter0(sex_male, &sample_uidx_base, &sex_male_invbits);
              cur_row[sample_uidx] = S_CAST(float, 0.0);
            }
          } else {
            // xchr_model 1: halve male values
            uintptr_t sex_male_bits = sex_male[0];
            for (uint32_t male_idx = 0; male_idx != male_ct; ++male_idx) {
              const uintptr_t sample_uidx = BitIter1(sex_male, &sample_uidx_base, &sex_male_bits);
              cur_row[sample_uidx] *= S_CAST(float, 0.5);
            }
          }
        }
      } else {
        double* cur_row = &(dosage_d_vmaj[row_idx * sample_ct]);
        PopulateRescaledDosage(genovec, dosage_present, dosage_main, slope, 0.0, missing_val, sample_ct, dosage_ct, cur_row);
        if (is_x_or_y) {
          // Instead of doing this for every variant, we could precompute
          // chrX/chrY weight matrices with male weights halved/nonmale weights
          // zeroed out.  But the number of chrY variants is typically small
          // enough (and how often will --xchr-model 1 be used, anyway?) that I
          // don't think it's worth it.
          uintptr_t sample_uidx_base = 0;
          if (is_y) {
            // zero out nonmale values
            uintptr_t sex_male_invbits = ~sex_male[0];
            for (uint32_t nonmale_idx = 0; nonmale_idx != nonmale_ct; ++nonmale_idx) {
              const uintptr_t sample_uidx = BitIter0(sex_male, &sample_uidx_base, &sex_male_invbits);
              cur_row[sample_uidx] = 0.0;
            }
          } else {
            // xchr_model 1: halve male values
            uintptr_t sex_male_bits = sex_male[0];
            for (uint32_t male_idx = 0; male_idx != male_ct; ++male_idx) {
              const uintptr_t sample_uidx = BitIter1(sex_male, &sample_uidx_base, &sex_male_bits);
              cur_row[sample_uidx] *= 0.5;
            }
          }
        }
      }
      ++row_idx;
      if (missing_cts) {
        ZeroTrailingNyps(sample_ct, genovec);
        uint32_t missing_ct;
        if (!dosage_ct) {
          if (!is_y) {
            missing_ct = GenoarrCountMissingUnsafe(genovec, sample_ct);
          } else {
            missing_ct = GenoarrCountMissingSubset(genovec, sex_male_interleaved_vec, sample_ct);
          }
        } else {
          if (!is_y) {
            missing_ct = GenoarrCountMissingInvsubsetUnsafe(genovec, dosage_present, sample_ct);
          } else {
            // include males, exclude dosages
            const uint32_t fullword_ct = (sample_ct + kBitsPerWordD2 - 1) / kBitsPerWord;
            missing_ct = 0;
            for (uint32_t widx = 0; widx != fullword_ct; ++widx) {
              uintptr_t w1 = genovec[2 * widx];
              uintptr_t w2 = genovec[2 * widx + 1];
              w1 = w1 & (w1 >> 1);
              w2 = w2 & (w2 >> 1);
              w1 = PackWordToHalfwordMask5555(w1);
              w2 = PackWordToHalfwordMask5555(w2);
              const uintptr_t ww = w1 | (w2 << kBitsPerWordD2);
              missing_ct += PopcountWord(ww & sex_male[widx] & (~dosage_present[widx]));
            }
            if (sample_ct > fullword_ct * kBitsPerWord) {
              uintptr_t w1 = genovec[2 * fullword_ct];
              w1 = w1 & (w1 >> 1);
              w1 = PackWordToHalfwordMask5555(w1);
              missing_ct += PopcountWord(w1 & sex_male[fullword_ct] & (~dosage_present[fullword_ct]));
            }
          }
        }
        missing_cts[variant_bidx] = missing_ct;
      }
    }
    if (row_idx) {
      if (single_prec) {
        RowMajorFmatrixMultiply(dosage_f_vmaj, wts_f_smaj, row_idx, vscore_ct, sample_ct, tmp_f_result_buf);
        const float* tmp_f_result_iter = tmp_f_result_buf;
        for (uintptr_t ulii = 0; ulii != row_idx; ++ulii) {
          uintptr_t cur_bidx = cur_bidxs[ulii];
          memcpy(&(cur_results_f[cur_bidx * vscore_ct]), tmp_f_result_iter, vscore_ct * sizeof(float));
          tmp_f_result_iter = &(tmp_f_result_iter[vscore_ct]);
        }
      } else {
        RowMajorMatrixMultiply(dosage_d_vmaj, wts_d_smaj, row_idx, vscore_ct, sample_ct, tmp_d_result_buf);
        const double* tmp_d_result_iter = tmp_d_result_buf;
        for (uintptr_t ulii = 0; ulii != row_idx; ++ulii) {
          uintptr_t cur_bidx = cur_bidxs[ulii];
          memcpy(&(cur_results_d[cur_bidx * vscore_ct]), tmp_d_result_iter, vscore_ct * sizeof(double));
          tmp_d_result_iter = &(tmp_d_result_iter[vscore_ct]);
        }
      }
    }
    parity = 1 - parity;
    while (0) {
    VscoreThread_err:
      UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr Vscore(const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_male, const double* allele_freqs, const char* in_fname, const RangeList* col_idx_range_listp, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t nosex_ct, uint32_t max_allele_slen, VscoreFlags flags, uint32_t xchr_model, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  uint32_t calc_thread_ct = max_thread_ct;
  char* cswritep = nullptr;
  FILE* binfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  ThreadGroup tg;
  CompressStreamState css;
  VscoreCtx ctx;
  PreinitTextStream(&txs);
  PreinitThreads(&tg);
  PreinitCstream(&css);
#ifdef USE_CUDA
  ctx.cfms = nullptr;
#endif
  {
    // unsurprisingly, lots of overlap with --score
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    if (nosex_ct) {
      // forced mean-imputation doesn't mix well with inclusion of chrY
      // unknown-sex samples.
      if (unlikely(XymtIsNonempty(variant_include, cip, kChrOffsetY))) {
        logerrputs("Error: When chrY is present, --variant-score cannot be used with unknown-sex\nsamples.\n");
        goto Vscore_ret_INCONSISTENT_INPUT;
      }
    }
    if (!xchr_model) {
      uint32_t x_code;
      if (XymtExists(cip, kChrOffsetX, &x_code)) {
        uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
        uint32_t x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
        uint32_t x_end = cip->chr_fo_vidx_start[x_chr_fo_idx + 1];
        if (!AllBitsAreZero(variant_include, x_start, x_end)) {
          uintptr_t* variant_include_no_x;
          if (unlikely(bigstack_alloc_w(raw_variant_ctl, &variant_include_no_x))) {
            goto Vscore_ret_NOMEM;
          }
          memcpy(variant_include_no_x, variant_include, raw_variant_ctl * sizeof(intptr_t));
          variant_ct -= PopcountBitRange(variant_include, x_start, x_end);
          if (!variant_ct) {
            logerrputs("Error: No --variant-score variants remaining after --xchr-model 0.\n");
            goto Vscore_ret_INCONSISTENT_INPUT;
          }
          ClearBitsNz(x_start, x_end, variant_include_no_x);
          variant_include = variant_include_no_x;
        }
      }
    } else if (xchr_model == 2) {
      xchr_model = 0;
    }
    // now xchr_model is set iff it's 1

    // see KeepColMatch() and SampleSortFileMap()
    char* line_start;
    XidMode xid_mode;
    reterr = OpenAndLoadXidHeader(in_fname, "variant-score", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeaderFixedWidth : kfXidHeaderFixedWidthIgnoreSid, kTextStreamBlenFast, &txs, &xid_mode, &line_idx, &line_start);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: Empty --variant-score file.\n");
        reterr = kPglRetMalformedInput;
      }
      goto Vscore_ret_1;
    }
    const uint32_t id_col_ct = GetXidColCt(xid_mode);
    const uint32_t col_ct = CountTokens(line_start);
    if (unlikely(id_col_ct == col_ct)) {
      logerrputs("Error: No score columns in --variant-score file.\n");
      goto Vscore_ret_MALFORMED_INPUT;
    }
    uintptr_t vscore_ct;
    uint32_t* col_idx_deltas;
    if (!col_idx_range_listp->name_ct) {
      vscore_ct = col_ct - id_col_ct;
      if (unlikely(bigstack_alloc_u32(vscore_ct, &col_idx_deltas))) {
        goto Vscore_ret_NOMEM;
      }
      for (uint32_t uii = 0; uii != vscore_ct; ++uii) {
        col_idx_deltas[uii] = 1;
      }
    } else {
      const uint32_t col_ctl = BitCtToWordCt(col_ct);
      uintptr_t* vscore_col_bitarr;
      if (unlikely(bigstack_calloc_w(col_ctl, &vscore_col_bitarr))) {
        goto Vscore_ret_NOMEM;
      }
      if (unlikely(NumericRangeListToBitarr(col_idx_range_listp, col_ct, 1, 0, vscore_col_bitarr))) {
        goto Vscore_ret_MISSING_TOKENS;
      }
      if (vscore_col_bitarr[0] & ((1 << id_col_ct) - 1)) {
        logerrputs("Error: --vscore-col-nums argument overlaps with ID columns.\n");
        goto Vscore_ret_INCONSISTENT_INPUT;
      }
      vscore_ct = PopcountWords(vscore_col_bitarr, col_ctl);
      // since we don't allow overflow, this should be guaranteed to be
      // positive
      assert(vscore_ct);
      if (unlikely(bigstack_alloc_u32(vscore_ct, &col_idx_deltas))) {
        goto Vscore_ret_NOMEM;
      }
      uintptr_t col_uidx_base = 0;
      uintptr_t vscore_col_bitarr_bits = vscore_col_bitarr[0];
      for (uintptr_t vscore_idx = 0; vscore_idx != vscore_ct; ++vscore_idx) {
        const uint32_t col_uidx = BitIter1(vscore_col_bitarr, &col_uidx_base, &vscore_col_bitarr_bits);
        col_idx_deltas[vscore_idx] = col_uidx;
      }
      // now convert to deltas
      for (uintptr_t vscore_idx = vscore_ct - 1; vscore_idx; --vscore_idx) {
        col_idx_deltas[vscore_idx] -= col_idx_deltas[vscore_idx - 1];
      }
      col_idx_deltas[0] -= id_col_ct - 1;
    }
    char** vscore_names;
    if (unlikely(bigstack_end_alloc_cp(vscore_ct, &vscore_names))) {
      goto Vscore_ret_NOMEM;
    }
    const uint32_t is_header_line = (line_start[0] == '#');
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = g_bigstack_end;
    if (is_header_line) {
      const char* name_iter = NextTokenMult0(line_start, id_col_ct - 1);
      for (uintptr_t vscore_idx = 0; vscore_idx != vscore_ct; ++vscore_idx) {
        name_iter = NextTokenMult(name_iter, col_idx_deltas[vscore_idx]);
        const char* name_end = CurTokenEnd(name_iter);
        // don't actually need to enforce unique names, though we could print a
        // warning later
        const uint32_t cur_slen = name_end - name_iter;
        if (cur_slen > kMaxIdSlen) {
          snprintf(g_logbuf, kLogbufSize, "Error: Variant-score name in column %" PRIuPTR " of %s is too long.\n", vscore_idx + id_col_ct + 1, in_fname);
          goto Vscore_ret_MALFORMED_INPUT_WW;
        }
        if (StoreStringAtEnd(tmp_alloc_base, name_iter, cur_slen, &tmp_alloc_end, &(vscore_names[vscore_idx]))) {
          goto Vscore_ret_NOMEM;
        }
        name_iter = name_end;
      }
      ++line_idx;
      line_start = TextGet(&txs);
    } else {
      for (uintptr_t vscore_num = 1; vscore_num <= vscore_ct; ++vscore_num) {
        const uint32_t cur_blen = 7 + UintSlen(vscore_num);
        if (PtrWSubCk(tmp_alloc_base, cur_blen, &tmp_alloc_end)) {
          goto Vscore_ret_NOMEM;
        }
        char* cur_name_iter = R_CAST(char*, tmp_alloc_end);
        vscore_names[vscore_num - 1] = cur_name_iter;
        cur_name_iter = strcpya_k(cur_name_iter, "VSCORE");
        cur_name_iter = u32toa(vscore_num, cur_name_iter);
        *cur_name_iter = '\0';
      }
    }
    BigstackEndSet(tmp_alloc_end);
    uint32_t* xid_map;
    char* sorted_xidbox;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto Vscore_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    char* idbuf;
    uint32_t* sample_uidx_order;
    if (unlikely(bigstack_alloc_c(siip->max_sample_id_blen, &idbuf) ||
                 bigstack_alloc_u32(sample_ct, &sample_uidx_order))) {
      goto Vscore_ret_NOMEM;
    }
    const uint32_t single_prec = (flags / kfVscoreSinglePrec) & 1;
    double* raw_wts_d = nullptr;
    float* raw_wts_f = nullptr;
    if (flags & kfVscoreSinglePrec) {
      if (unlikely(bigstack_alloc64_f(sample_ct * S_CAST(uint64_t, vscore_ct), &raw_wts_f))) {
        goto Vscore_ret_NOMEM;
      }
    } else {
      if (unlikely(bigstack_alloc64_d(sample_ct * S_CAST(uint64_t, vscore_ct), &raw_wts_d))) {
        goto Vscore_ret_NOMEM;
      }
    }
    uintptr_t* already_seen;
    if (unlikely(bigstack_end_calloc_w(raw_sample_ctl, &already_seen))) {
      goto Vscore_ret_NOMEM;
    }
    uintptr_t miss_ct = 0;
    uint32_t hit_ct = 0;

    double* raw_wts_d_iter = raw_wts_d;
    float* raw_wts_f_iter = raw_wts_f;
    for (; line_start; ++line_idx, line_start = TextGet(&txs)) {
      if (unlikely(line_start[0] == '#')) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --variant-score file starts with a '#'. (This is only permitted before the first nonheader line, and if a #FID/IID header line is present it must denote the end of the header block.)\n", line_idx);
        goto Vscore_ret_MALFORMED_INPUT_WW;
      }
      const char* linebuf_iter = line_start;
      uint32_t sample_uidx;
      if (SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, sample_ct, 0, xid_mode, &linebuf_iter, &sample_uidx, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto Vscore_ret_MISSING_TOKENS;
        }
        ++miss_ct;
        continue;
      }
      if (unlikely(IsSet(already_seen, sample_uidx))) {
        TabsToSpaces(idbuf);
        snprintf(g_logbuf, kLogbufSize, "Error: Duplicate sample ID \"%s\" in --variant-score file.\n", idbuf);
        goto Vscore_ret_MALFORMED_INPUT_WW;
      }
      SetBit(sample_uidx, already_seen);
      sample_uidx_order[hit_ct] = sample_uidx;
      for (uintptr_t vscore_idx = 0; vscore_idx != vscore_ct; ++vscore_idx) {
        linebuf_iter = NextTokenMult(linebuf_iter, col_idx_deltas[vscore_idx]);
        if (unlikely(!linebuf_iter)) {
          goto Vscore_ret_MISSING_TOKENS;
        }
        double dxx;
        const char* token_end = ScantokDouble(linebuf_iter, &dxx);
        if (unlikely((!token_end) || (single_prec && (fabs(dxx) > 3.4028235677973362e38)))) {
          token_end = CurTokenEnd(linebuf_iter);
          *K_CAST(char*, token_end) = '\0';
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid coefficient '%s' on line %" PRIuPTR " of --variant-score file.\n", linebuf_iter, line_idx);
          goto Vscore_ret_MALFORMED_INPUT_WW;
        }
        if (single_prec) {
          *raw_wts_f_iter++ = S_CAST(float, dxx);
        } else {
          *raw_wts_d_iter++ = dxx;
        }
        linebuf_iter = token_end;
      }
      ++hit_ct;
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto Vscore_ret_TSTREAM_FAIL;
    }
    if (unlikely(CleanupTextStream2(in_fname, &txs, &reterr))) {
      goto Vscore_ret_1;
    }
    if (!hit_ct) {
      logerrputs("Error: No valid entries in --variant-score file.\n");
      goto Vscore_ret_INCONSISTENT_INPUT;
    }
    sample_include = already_seen;
    sample_ct = hit_ct;
#if defined(__LP64__) && !defined(LAPACK_ILP64)
    if (sample_ct * vscore_ct > 0x7fffffff) {
      logerrputs("Error: --variant-score input matrix too large for this " PROG_NAME_STR " build.  If this\nis really the computation you want, use a " PROG_NAME_STR " build with large-matrix\nsupport.\n");
      goto Vscore_ret_INCONSISTENT_INPUT;
    }
#endif
    ctx.variant_include = variant_include;
    ctx.cip = cip;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.allele_freqs = allele_freqs;
    ctx.sample_include = sample_include;
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    uint32_t compress_thread_ct = 1;
    const uint32_t output_zst = (flags / kfVscoreZs) & 1;
    snprintf(outname_end, kMaxOutfnameExtBlen, ".vscore");
    if (flags & (kfVscoreBin | kfVscoreBin4)) {
      snprintf(&(outname_end[7]), kMaxOutfnameExtBlen - 7, ".cols");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &binfile))) {
        goto Vscore_ret_OPEN_FAIL;
      }
      for (uintptr_t vscore_idx = 0; vscore_idx != vscore_ct; ++vscore_idx) {
        fputs(vscore_names[vscore_idx], binfile);
#ifdef _WIN32
        putc_unlocked('\r', binfile);
#endif
        putc_unlocked('\n', binfile);
      }
      if (unlikely(fclose_null(&binfile))) {
        goto Vscore_ret_WRITE_FAIL;
      }
      snprintf(&(outname_end[7]), kMaxOutfnameExtBlen - 7, ".bin");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &binfile))) {
        goto Vscore_ret_OPEN_FAIL;
      }
      snprintf(&(outname_end[7]), kMaxOutfnameExtBlen - 7, ".vars");
      if (output_zst) {
        snprintf(&(outname_end[12]), kMaxOutfnameExtBlen - 12, ".zst");
      }
    } else if (output_zst) {
      snprintf(&(outname_end[7]), kMaxOutfnameExtBlen - 7, ".zst");
      if (calc_thread_ct > 1) {
        // The more samples there are, the higher the compute:compress ratio we
        // want, though this is not a linear relationship due to the sparse
        // optimization.
        // 1:1 split seems to work well for a few thousand samples; I'm
        // guessing that ~7:1 is better for hundreds of thousands.
        if (sample_ct < 8192) {
          compress_thread_ct = calc_thread_ct / 2;
        } else {
          const uint32_t log2_sample_ct_m10 = bsru32(sample_ct) - 10;
          // 3/8, 4/16, 5/24, ...
          compress_thread_ct = (calc_thread_ct * log2_sample_ct_m10) / (8 * (log2_sample_ct_m10 - 2));
          if (!compress_thread_ct) {
            compress_thread_ct = 1;
          }
        }
        calc_thread_ct -= compress_thread_ct;
      }
    }
    {
      uint32_t* sample_include_cumulative_popcounts;
      if (unlikely(bigstack_end_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts))) {
        goto Vscore_ret_NOMEM;
      }
      float* wts_f_smaj = nullptr;
      double* wts_d_smaj = nullptr;
      if (single_prec) {
        if (unlikely(bigstack_end_alloc_f(sample_ct * vscore_ct, &wts_f_smaj))) {
          goto Vscore_ret_NOMEM;
        }
      } else {
        if (unlikely(bigstack_end_alloc_d(sample_ct * vscore_ct, &wts_d_smaj))) {
          goto Vscore_ret_NOMEM;
        }
      }
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
      logprintfww("--variant-score: %" PRIuPTR " score-vector%s loaded for %u sample%s.\n", vscore_ct, (vscore_ct == 1)? "" : "s", sample_ct, (sample_ct == 1)? "" : "s");
      if (miss_ct) {
        logerrprintf("Warning: %" PRIuPTR " line%s skipped in --variant-score file.\n", miss_ct, (miss_ct == 1)? "" : "s");
      }
      if (single_prec) {
        const float* wts_read_iter = raw_wts_f;
        for (uint32_t uii = 0; uii != sample_ct; ++uii) {
          const uint32_t sample_uidx = sample_uidx_order[uii];
          const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_uidx);
          memcpy(&(wts_f_smaj[sample_idx * vscore_ct]), wts_read_iter, vscore_ct * sizeof(float));
          wts_read_iter = &(wts_read_iter[vscore_ct]);
        }
        ctx.wts_f_smaj = wts_f_smaj;
        ctx.wts_d_smaj = nullptr;
      } else {
        const double* wts_read_iter = raw_wts_d;
        for (uint32_t uii = 0; uii != sample_ct; ++uii) {
          const uint32_t sample_uidx = sample_uidx_order[uii];
          const uint32_t sample_idx = RawToSubsettedPos(sample_include, sample_include_cumulative_popcounts, sample_uidx);
          memcpy(&(wts_d_smaj[sample_idx * vscore_ct]), wts_read_iter, vscore_ct * sizeof(double));
          wts_read_iter = &(wts_read_iter[vscore_ct]);
        }
        ctx.wts_f_smaj = nullptr;
        ctx.wts_d_smaj = wts_d_smaj;
      }
      BigstackReset(bigstack_mark);
#ifdef USE_CUDA
      if (single_prec && (vscore_ct >= 80)) {
        const uint32_t device_count = CudaGetDeviceCount();
        if (device_count) {
          if (unlikely(BIGSTACK_ALLOC_X(CublasFmultiplier, calc_thread_ct, &ctx.cfms))) {
            goto Vscore_ret_NOMEM;
          }
          uint32_t tidx = 0;
          CublasFmultiplierPreinit(&ctx.cfms[0]);
          for (uint32_t device_idx = 0; device_idx != device_count; ++device_idx) {
            if (unlikely(CudaSetDevice(device_idx))) {
              reterr = kPglRetGpuFail;
              logerrputs("Error: GPU operation failure.\n");
              goto Vscore_ret_1;
            }
            if (CublasFmultiplierRowMajorInit(kVscoreBlockSize, vscore_ct, sample_ct, &ctx.cfms[tidx])) {
              // Insufficient memory on this device.  Try the next one.
              continue;
            }
            CublasFmultiplierPreloadRowMajor2(ctx.wts_f_smaj, &ctx.cfms[tidx]);
            ++tidx;
            if (tidx == calc_thread_ct) {
              break;
            }
            CublasFmultiplierPreinit(&ctx.cfms[tidx]);
            CublasFmultiplierBorrowRowMajor2(&ctx.cfms[tidx - 1], &ctx.cfms[tidx]);
            if (CublasFmultiplierRowMajorInit(kVscoreBlockSize, vscore_ct, sample_ct, &ctx.cfms[tidx])) {
              continue;
            }
            ++tidx;
            if (tidx == calc_thread_ct) {
              break;
            }
            CublasFmultiplierPreinit(&ctx.cfms[tidx]);
          }
          if (!tidx) {
            logputs("Note: Not using GPU for --variant-score computation, due to insufficient\ndevice memory.  If this is a problem, try providing fewer scores at a time.\n");
            ctx.cfms = nullptr;
          } else {
            calc_thread_ct = tidx;
            logprintf("--variant-score: Using %u GPU handle%s.\n", calc_thread_ct, (calc_thread_ct == 1)? "" : "s");
          }
        }
      }
#endif
      const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
      uintptr_t* sex_male_collapsed;
      uintptr_t* sex_male_interleaved_vec;
      if (unlikely(bigstack_alloc_w(sample_ctl, &sex_male_collapsed) ||
                   bigstack_alloc_w(sample_ctv * kWordsPerVec, &sex_male_interleaved_vec) ||
                   bigstack_alloc_wp(calc_thread_ct, &ctx.raregenos) ||
                   bigstack_alloc_u32p(calc_thread_ct, &ctx.difflist_sample_id_bufs))) {
        goto Vscore_ret_NOMEM;
      }
      if (single_prec) {
        if (unlikely(bigstack_alloc_fp(calc_thread_ct, &ctx.dosage_f_vmaj_bufs) ||
                     bigstack_alloc_fp(calc_thread_ct, &ctx.tmp_f_result_bufs))) {
          goto Vscore_ret_NOMEM;
        }
        ctx.dosage_d_vmaj_bufs = nullptr;
        ctx.tmp_d_result_bufs = nullptr;
      } else {
        if (unlikely(bigstack_alloc_dp(calc_thread_ct, &ctx.dosage_d_vmaj_bufs) ||
                     bigstack_alloc_dp(calc_thread_ct, &ctx.tmp_d_result_bufs))) {
          goto Vscore_ret_NOMEM;
        }
        ctx.dosage_f_vmaj_bufs = nullptr;
        ctx.tmp_f_result_bufs = nullptr;
      }
      CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_male_collapsed);
      // bugfix (18 Oct 2023): forgot to clear trailing words
      ZeroTrailingWords(sample_ctl, sex_male_collapsed);
      FillInterleavedMaskVec(sex_male_collapsed, sample_ctv, sex_male_interleaved_vec);
      ctx.sex_male_collapsed = sex_male_collapsed;
      ctx.sex_male_interleaved_vec = sex_male_interleaved_vec;
    }
    ctx.vscore_ct = vscore_ct;
    ctx.sample_ct = sample_ct;
    const uint32_t male_ct = PopcountWords(ctx.sex_male_collapsed, sample_ctl);
    ctx.male_ct = male_ct;
    ctx.is_xchr_model_1 = xchr_model;

    const uint32_t chr_col = (flags / kfVscoreColChrom) & 1;
    char* chr_buf = nullptr;
    uint32_t max_chr_blen = 0;
    if (chr_col) {
      max_chr_blen = GetMaxChrSlen(cip) + 1;
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto Vscore_ret_NOMEM;
      }
    }
    const uint32_t ref_col = (flags / kfVscoreColRef) & 1;
    const uint32_t alt1_col = (flags / kfVscoreColAlt1) & 1;
    const uint32_t alt_col = (flags / kfVscoreColAlt) & 1;
    uintptr_t overflow_buf_size;
    if (binfile) {
      overflow_buf_size = kCompressStreamBlock + kMaxIdSlen + 16;
    } else {
      overflow_buf_size = kCompressStreamBlock + max_chr_blen * chr_col + kMaxIdSlen + 128 + (24 * k1LU) * vscore_ct + MAXV(ref_col + alt1_col, alt_col) * max_allele_slen;
    }
    reterr = InitCstreamAlloc(outname, 0, output_zst, compress_thread_ct, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto Vscore_ret_1;
    }
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t all_nonref = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t provref_col = ref_col && ProvrefCol(variant_include, nonref_flags, flags / kfVscoreColMaybeprovref, raw_variant_ct, all_nonref);
    const uint32_t nmiss_col = (flags / kfVscoreColNmiss) & 1;
    const uint32_t nobs_col = (flags / kfVscoreColNobs) & 1;
    if (!binfile) {
      *cswritep++ = '#';
      if (chr_col) {
        cswritep = strcpya_k(cswritep, "CHROM\t");
      }
      if (flags & kfVscoreColPos) {
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
      if (flags & kfVscoreColAltfreq) {
        cswritep = strcpya_k(cswritep, "\tALT_FREQ");
      } else {
        allele_freqs = nullptr;
      }
      if (nmiss_col) {
        cswritep = strcpya_k(cswritep, "\tMISSING_CT");
      }
      if (nobs_col) {
        cswritep = strcpya_k(cswritep, "\tOBS_CT");
      }
      for (uintptr_t vscore_idx = 0; vscore_idx != vscore_ct; ++vscore_idx) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, vscore_names[vscore_idx]);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto Vscore_ret_WRITE_FAIL;
        }
      }
      AppendBinaryEoln(&cswritep);
    }

    if (nmiss_col || nobs_col) {
      if (unlikely(bigstack_alloc_u32(kPglVblockSize, &ctx.missing_cts[0]) ||
                   bigstack_alloc_u32(kPglVblockSize, &ctx.missing_cts[1]))) {
        goto Vscore_ret_NOMEM;
      }
    } else {
      ctx.missing_cts[0] = nullptr;
      ctx.missing_cts[1] = nullptr;
    }

#ifdef USE_CUDA
    const uint32_t max_difflist_len = sample_ct / (single_prec? (ctx.cfms? 64 : 32) : 16);
#else
    const uint32_t max_difflist_len = sample_ct / (single_prec? 32 : 16);
#endif
    ctx.max_difflist_len = max_difflist_len;
    // * Per-thread raregeno buffers must have space for
    //   max_difflist_len nyps, and difflist_sample_ids buffers need space for
    //   that many uint32s.
    // * Per-thread dosage_vmaj buffers must have space for
    //   kVscoreBlockSize * sample_ct elements.
    // * Per-thread result buffers must have space for kVscoreBlockSize *
    //   vscore_ct elements.
    const uintptr_t thread_xalloc_cacheline_ct = DivUp(max_difflist_len, kNypsPerCacheline) + DivUp(max_difflist_len, kInt32PerCacheline) + DivUp(kVscoreBlockSize * S_CAST(uintptr_t, sample_ct) * sizeof(double), kCacheline) + DivUp(kVscoreBlockSize * vscore_ct * sizeof(double), kCacheline);

    // ctx.results must have space for 2 * vscore_ct * read_block_size values.
    const uintptr_t per_variant_xalloc_byte_ct = 2 * vscore_ct * ((8 * k1LU) >> single_prec);
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    // defensive
    ctx.dosage_presents = nullptr;
    ctx.dosage_mains = nullptr;
    uint32_t read_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, per_variant_xalloc_byte_ct, 0, pgfip, &calc_thread_ct, &ctx.genovecs, nullptr, nullptr, nullptr, dosage_is_present? (&ctx.dosage_presents) : nullptr, dosage_is_present? (&ctx.dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto Vscore_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto Vscore_ret_NOMEM;
    }
    {
      // could vector-align individual allocations and only cacheline-align at
      // thread boundaries, but the savings are microscopic
      const uintptr_t raregeno_alloc = kCacheline * DivUp(max_difflist_len, kNypsPerCacheline);
      const uintptr_t difflist_sample_ids_alloc = RoundUpPow2(max_difflist_len * sizeof(int32_t), kCacheline);
      const uintptr_t dosage_vmaj_alloc = RoundUpPow2(kVscoreBlockSize * S_CAST(uintptr_t, sample_ct) * ((8 * k1LU) >> single_prec), kCacheline);
      const uintptr_t tmp_result_alloc = RoundUpPow2(kVscoreBlockSize * vscore_ct * ((8 * k1LU) >> single_prec), kCacheline);
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        ctx.raregenos[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(raregeno_alloc));
        ctx.difflist_sample_id_bufs[tidx] = S_CAST(uint32_t*, bigstack_alloc_raw(difflist_sample_ids_alloc));
        if (single_prec) {
          ctx.dosage_f_vmaj_bufs[tidx] = S_CAST(float*, bigstack_alloc_raw(dosage_vmaj_alloc));
          ctx.tmp_f_result_bufs[tidx] = S_CAST(float*, bigstack_alloc_raw(tmp_result_alloc));
        } else {
          ctx.dosage_d_vmaj_bufs[tidx] = S_CAST(double*, bigstack_alloc_raw(dosage_vmaj_alloc));
          ctx.tmp_d_result_bufs[tidx] = S_CAST(double*, bigstack_alloc_raw(tmp_result_alloc));
        }
      }
    }
    const uintptr_t results_byte_ct = RoundUpPow2(per_variant_xalloc_byte_ct * read_block_size, kCacheline);
    if (single_prec) {
      ctx.results_f[0] = S_CAST(float*, bigstack_alloc_raw(results_byte_ct));
      ctx.results_f[1] = S_CAST(float*, bigstack_alloc_raw(results_byte_ct));
      ctx.results_d[0] = nullptr;
      ctx.results_d[1] = nullptr;
    } else {
      ctx.results_f[0] = nullptr;
      ctx.results_f[1] = nullptr;
      ctx.results_d[0] = S_CAST(double*, bigstack_alloc_raw(results_byte_ct));
      ctx.results_d[1] = S_CAST(double*, bigstack_alloc_raw(results_byte_ct));
    }
    assert(g_bigstack_base <= g_bigstack_end);
    ctx.err_info = (~0LLU) << 32;
    SetThreadFuncAndData(VscoreThread, &ctx, &tg);

    fputs("--variant-score: 0%", stdout);
    fflush(stdout);

    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
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
    // 8. Write results for last block
    uintptr_t write_variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t prev_block_size = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t cur_sample_ct = 0;
    uint32_t cur_allele_ct = 2;
    for (uint32_t variant_idx = 0; ; ++read_block_idx) {
      const uint32_t cur_block_size = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto Vscore_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, ctx.err_info);
        if (unlikely(reterr)) {
#ifdef USE_CUDA
          if (reterr == kPglRetGpuFail) {
            logputs("\n");
            logerrputs("Error: GPU operation failure.\n");
            goto Vscore_ret_1;
          }
#endif
          PgenErrPrintNV(reterr, ctx.err_info >> 32);
          goto Vscore_ret_1;
        }
      }
      if (!IsLastBlock(&tg)) {
        // it may make sense to put this boilerplate into its own function,
        // too...
        ctx.cur_block_size = cur_block_size;
        ComputeUidxStartPartition(variant_include, cur_block_size, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_size == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto Vscore_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const float* cur_results_f_iter = ctx.results_f[parity];
        const double* cur_results_d_iter = ctx.results_d[parity];
        const uint32_t* cur_missing_cts = ctx.missing_cts[parity];
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_size; ++variant_bidx) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            cur_sample_ct = (chr_idx == y_code)? male_ct : sample_ct;
            if (chr_buf) {
              char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
              *chr_name_end = '\t';
              chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
            }
          }
          if (binfile) {
            // may as well write variant-ID file in this loop
            cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
            AppendBinaryEoln(&cswritep);
            if (unlikely(Cswrite(&css, &cswritep))) {
              goto Vscore_ret_WRITE_FAIL;
            }
            continue;
          }
          if (chr_col) {
            cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
          }
          if (variant_bps) {
            cswritep = u32toa_x(variant_bps[write_variant_uidx], '\t', cswritep);
          }
          cswritep = strcpya(cswritep, variant_ids[write_variant_uidx]);
          uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            cur_allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offset_base;
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
                goto Vscore_ret_WRITE_FAIL;
              }
              cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
            }
            --cswritep;
          }
          if (provref_col) {
            *cswritep++ = '\t';
            *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, write_variant_uidx)))? 'Y' : 'N';
          }
          if (allele_freqs) {
            *cswritep++ = '\t';
            cswritep = dtoa_g(1.0 - allele_freqs[allele_idx_offset_base - write_variant_uidx], cswritep);
          }
          if (nmiss_col) {
            *cswritep++ = '\t';
            cswritep = u32toa(cur_missing_cts[variant_bidx], cswritep);
          }
          if (nobs_col) {
            *cswritep++ = '\t';
            cswritep = u32toa(cur_sample_ct - cur_missing_cts[variant_bidx], cswritep);
          }
          if (single_prec) {
            for (uintptr_t vscore_idx = 0; vscore_idx != vscore_ct; ++vscore_idx) {
              *cswritep++ = '\t';
              cswritep = ftoa_g(*cur_results_f_iter++, cswritep);
            }
          } else {
            for (uintptr_t vscore_idx = 0; vscore_idx != vscore_ct; ++vscore_idx) {
              *cswritep++ = '\t';
              cswritep = dtoa_g(*cur_results_d_iter++, cswritep);
            }
          }
          AppendBinaryEoln(&cswritep);
          if (unlikely(Cswrite(&css, &cswritep))) {
            goto Vscore_ret_WRITE_FAIL;
          }
        }
        if (binfile) {
          const uintptr_t entry_ct = vscore_ct * prev_block_size;
          if (single_prec) {
            // must be bin4
            if (unlikely(fwrite_checked(cur_results_f_iter, entry_ct * sizeof(float), binfile))) {
              goto Vscore_ret_WRITE_FAIL;
            }
          } else if (flags & kfVscoreBin) {
            if (unlikely(fwrite_checked(cur_results_d_iter, entry_ct * sizeof(double), binfile))) {
              goto Vscore_ret_WRITE_FAIL;
            }
          } else {
            const uintptr_t fullgroup_ct = entry_ct / 32768;
            float buf[32768];
            for (uintptr_t group_idx = 0; group_idx != fullgroup_ct; ++group_idx) {
              for (uint32_t uii = 0; uii != 32768; ++uii) {
                buf[uii] = S_CAST(float, *cur_results_d_iter++);
              }
              if (unlikely(fwrite_checked(buf, 32768 * sizeof(float), binfile))) {
                goto Vscore_ret_WRITE_FAIL;
              }
            }
            const uint32_t remainder = entry_ct % 32768;
            if (remainder) {
              for (uint32_t uii = 0; uii != remainder; ++uii) {
                buf[uii] = S_CAST(float, *cur_results_d_iter++);
              }
              if (unlikely(fwrite_checked(buf, remainder * sizeof(float), binfile))) {
                goto Vscore_ret_WRITE_FAIL;
              }
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
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
        }
      }
      prev_block_size = cur_block_size;
      variant_idx += cur_block_size;
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(CswriteCloseNull(&css, cswritep))) {
      goto Vscore_ret_WRITE_FAIL;
    }
    putc_unlocked('\r', stdout);
    if (!binfile) {
      logprintfww("--variant-score: Results written to %s .\n", outname);
    } else {
      if (unlikely(fclose_null(&binfile))) {
        goto Vscore_ret_WRITE_FAIL;
      }
      outname_end[8] = '\0';
      logprintfww("--variant-score: Score matrix written to %sbin , and associated column and variant ID labels written to %scols and %svars%s .\n", outname, outname, outname, output_zst? ".zst" : "");
    }
  }
  while (0) {
  Vscore_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Vscore_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  Vscore_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--variant-score file", &txs);
    break;
  Vscore_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  Vscore_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  Vscore_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  Vscore_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  Vscore_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of --variant-score file has fewer tokens than expected.\n", line_idx);
  Vscore_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  Vscore_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 Vscore_ret_1:
#ifdef USE_CUDA
  if (ctx.cfms) {
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      if (unlikely(CublasFmultiplierCleanup(&ctx.cfms[tidx]))) {
        if (!reterr) {
          reterr = kPglRetGpuFail;
        }
      }
    }
  }
#endif
  fclose_cond(binfile);
  CswriteCloseCond(&css, cswritep);
  CleanupThreads(&tg);
  CleanupTextStream2("--variant-score file", &txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

#ifndef NOLAPACK
// may need to increase textbuf size if this stops being true
static_assert(kMaxPc < 100000, "PhenoSvd() needs to be updated.");
PglErr PhenoSvd(const PhenoSvdInfo* psip, const uintptr_t* sample_include, const SampleIdInfo* siip, uint32_t raw_sample_ct, uint32_t orig_sample_ct, __attribute__((unused)) uint32_t max_thread_ct, char** pheno_names_ptr, uint32_t* pheno_ct_ptr, uintptr_t* max_pheno_name_blen_ptr, PhenoCol* pheno_cols, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  const uintptr_t orig_pheno_ct = *pheno_ct_ptr;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely(orig_pheno_ct < 2)) {
      logerrprintf("Error: --pheno-svd invoked %s.\n", orig_pheno_ct? "with only 1 phenotype" : "without any phenotypes");
      goto PhenoSvd_ret_INCONSISTENT_INPUT;
    }
    if (unlikely(orig_sample_ct <= 1)) {
      assert(orig_sample_ct == 1);
      logerrputs("Error: --pheno-svd invoked with only one sample.\n");
      goto PhenoSvd_ret_INCONSISTENT_INPUT;
    }
    const uint32_t raw_sample_ctaw = BitCtToAlignedWordCt(raw_sample_ct);
    uintptr_t* sample_intersect;
    if (unlikely(bigstack_alloc_w(raw_sample_ctaw, &sample_intersect))) {
      goto PhenoSvd_ret_NOMEM;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    memcpy(sample_intersect, sample_include, raw_sample_ctl * sizeof(intptr_t));
    ZeroTrailingWords(raw_sample_ctl, sample_intersect);
    const uintptr_t orig_max_pheno_name_blen = *max_pheno_name_blen_ptr;
    char* orig_pheno_names = *pheno_names_ptr;
    for (uintptr_t pheno_idx = 0; pheno_idx != orig_pheno_ct; ++pheno_idx) {
      const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
      if (unlikely(cur_pheno_col->type_code == kPhenoDtypeCat)) {
        logerrprintfww("Error: --pheno-svd: phenotype '%s' is categorical; it must be preprocessed with e.g. --split-cat-pheno.\n", &(orig_pheno_names[pheno_idx * orig_max_pheno_name_blen]));
        goto PhenoSvd_ret_INCONSISTENT_INPUT;
      }
      BitvecAnd(cur_pheno_col->nonmiss, raw_sample_ctl, sample_intersect);
    }
    const uint32_t new_sample_ct = PopcountWords(sample_intersect, raw_sample_ctl);
    const PhenoSvdFlags flags = psip->flags;
    if (new_sample_ct * 2 < orig_sample_ct) {
      if (unlikely(!(flags & kfPhenoSvdForce))) {
        logerrprintfww("Error: --pheno-svd: Only %u/%u sample%s have no missing phenotype values. Consider imputing some missing phenotype values, and/or excluding phenotypes with many missing values.\n", new_sample_ct, orig_sample_ct, (new_sample_ct == 1)? "" : "s");
        goto PhenoSvd_ret_INCONSISTENT_INPUT;
      }
    }
    const uintptr_t svd_dim = MINV(orig_pheno_ct, new_sample_ct);
    uint32_t new_pheno_ct = psip->ct;
    if (unlikely(svd_dim < new_pheno_ct)) {
      if (svd_dim == orig_pheno_ct) {
        logerrprintf("Error: --pheno-svd %u invoked with only %" PRIuPTR " phenotypes.\n", new_pheno_ct, orig_pheno_ct);
      } else {
        logerrprintf("Error: --pheno-svd %u invoked with only %u samples.\n", new_pheno_ct, new_sample_ct);
      }
      goto PhenoSvd_ret_INCONSISTENT_INPUT;
    }
    const uint64_t svd_rect_size = S_CAST(uint64_t, orig_pheno_ct) * new_sample_ct;
    lapack_int svd_rect_lwork;
#ifdef LAPACK_ILP64
    GetSvdRectLwork(orig_pheno_ct, new_sample_ct, &svd_rect_lwork);
#else
    if (unlikely((svd_rect_size > 0x7effffff) ||
                 GetSvdRectLwork(orig_pheno_ct, new_sample_ct, &svd_rect_lwork))) {
      logerrputs("Error: --pheno-svd problem instance too large for this " PROG_NAME_STR " build.  If this\nis really the computation you want, use a " PROG_NAME_STR " build with large-matrix support.\n");
      goto PhenoSvd_ret_INCONSISTENT_INPUT;
    }
#endif
    double* matrix; // orig_pheno x new_sample on input, pheno-pheno afterward
    double* singular_values;
    double* svd_sample_result;
    unsigned char* svd_rect_wkspace;
    if (unlikely(bigstack_alloc_d(svd_rect_size, &matrix) ||
                 bigstack_alloc_d(svd_dim, &singular_values) ||
                 bigstack_alloc_d(svd_rect_size, &svd_sample_result) ||
                 bigstack_alloc_uc(svd_rect_lwork * sizeof(double), &svd_rect_wkspace))) {
      goto PhenoSvd_ret_NOMEM;
    }
    {
      double* dwrite_iter = matrix;
      for (uintptr_t pheno_idx = 0; pheno_idx != orig_pheno_ct; ++pheno_idx) {
        const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);

        uintptr_t sample_uidx_base = 0;
        uintptr_t cur_bits = sample_intersect[0];
        if (cur_pheno_col->type_code == kPhenoDtypeQt) {
          const double* pheno_qt = cur_pheno_col->data.qt;
          for (uint32_t sample_idx = 0; sample_idx != new_sample_ct; ++sample_idx) {
            const uint32_t sample_uidx = BitIter1(sample_intersect, &sample_uidx_base, &cur_bits);
            dwrite_iter[sample_idx] = pheno_qt[sample_uidx];
          }
        } else {
          const uintptr_t* pheno_cc = cur_pheno_col->data.cc;
          for (uint32_t sample_idx = 0; sample_idx != new_sample_ct; ++sample_idx) {
            const uint32_t sample_uidx = BitIter1(sample_intersect, &sample_uidx_base, &cur_bits);
            dwrite_iter[sample_idx] = kSmallDoubles[IsSet(pheno_cc, sample_uidx)];
          }
        }
        dwrite_iter = &(dwrite_iter[new_sample_ct]);
      }
    }

    BLAS_SET_NUM_THREADS(max_thread_ct);
    IntErr svd_rect_err = SvdRect(orig_pheno_ct, new_sample_ct, svd_rect_lwork, matrix, singular_values, svd_sample_result, svd_rect_wkspace);
    if (unlikely(svd_rect_err)) {
      logerrprintf("Error: --pheno-svd: DGESVD failure, info=%d.\n", S_CAST(int32_t, svd_rect_err));
      goto PhenoSvd_ret_DEGENERATE_DATA;
    }

    if (!new_pheno_ct) {
      const double min_variance_explained = psip->min_variance_explained;
      if (min_variance_explained == 1.0) {
        new_pheno_ct = svd_dim;
      } else {
        double ssq = 0.0;
        for (uintptr_t sv_idx = 0; sv_idx != svd_dim; ++sv_idx) {
          const double cur_sv = singular_values[sv_idx];
          ssq += cur_sv * cur_sv;
        }
        const double target_ssq = min_variance_explained * ssq;
        ssq = 0.0;
        do {
          const double cur_sv = singular_values[new_pheno_ct];
          ssq += cur_sv * cur_sv;
          ++new_pheno_ct;
        } while (ssq < target_ssq);
        assert(new_pheno_ct <= svd_dim);
      }
      if (unlikely(new_pheno_ct > kMaxPc)) {
        logerrprintfww("Error: --pheno-svd variance=%g: Too many new phenotypes (%u; max %u).\n", min_variance_explained, new_pheno_ct, kMaxPc);
        goto PhenoSvd_ret_INCONSISTENT_INPUT;
      }
      logprintf("--pheno-svd variance=%g: %u/%" PRIuPTR " phenotype%s kept.\n", min_variance_explained, new_pheno_ct, svd_dim, (new_pheno_ct == 1)? "" : "s");
    }

    BigstackReset(svd_rect_wkspace);
    const uint32_t textbuf_req = kMaxMediumLine + MAXV(kMaxMediumLine, 3 * kMaxIdBlen + (kMaxDoubleGSlen + 1) * orig_pheno_ct);
    char* textbuf;
    if (unlikely(bigstack_alloc_c(textbuf_req, &textbuf))) {
      goto PhenoSvd_ret_NOMEM;
    }
    char* textbuf_flush = &(textbuf[kMaxMediumLine]);

    snprintf(outname_end, kMaxOutfnameExtBlen, ".svd.pheno");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto PhenoSvd_ret_OPEN_FAIL;
    }
    char* write_iter = textbuf;
    *write_iter++ = '#';
    const uint32_t col_fid = FidColIsRequired(siip, flags / kfPhenoSvdScolMaybefid);
    if (col_fid) {
      write_iter = strcpya_k(write_iter, "FID\t");
    }
    write_iter = strcpya_k(write_iter, "IID");
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    const uint32_t col_sid = SidColIsRequired(sids, flags / kfPhenoSvdScolMaybesid);
    if (col_sid) {
      write_iter = strcpya_k(write_iter, "\tSID");
    }
    for (uint32_t new_pheno_idx = 0; new_pheno_idx != new_pheno_ct; ++new_pheno_idx) {
      write_iter = strcpya_k(write_iter, "\tSVDPHENO");
      write_iter = u32toa(1 + new_pheno_idx, write_iter);
    }
    AppendBinaryEoln(&write_iter);
    if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
      goto PhenoSvd_ret_WRITE_FAIL;
    }
    {
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = sample_intersect[0];
      for (uint32_t sample_idx = 0; sample_idx != new_sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_intersect, &sample_uidx_base, &cur_bits);
        write_iter = AppendXid(sample_ids, sids, col_fid, col_sid, max_sample_id_blen, max_sid_blen, sample_uidx, write_iter);
        const double* svd_sample_result_row = &(svd_sample_result[sample_idx * svd_dim]);
        for (uint32_t new_pheno_idx = 0; new_pheno_idx != new_pheno_ct; ++new_pheno_idx) {
          *write_iter++ = '\t';
          write_iter = dtoa_g(svd_sample_result_row[new_pheno_idx], write_iter);
        }
        AppendBinaryEoln(&write_iter);
        if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
          goto PhenoSvd_ret_WRITE_FAIL;
        }
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto PhenoSvd_ret_WRITE_FAIL;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".svd.pheno_wts");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto PhenoSvd_ret_OPEN_FAIL;
    }
    write_iter = textbuf;
    *write_iter++ = '#';
    const uint32_t col_id = (flags / kfPhenoSvdPcolId) & 1;
    if (col_id) {
      write_iter = strcpya_k(write_iter, "NEW_PHENO_ID\t");
    }
    const uint32_t col_sv = (flags / kfPhenoSvdPcolSv) & 1;
    if (col_sv) {
      write_iter = strcpya_k(write_iter, "SINGULAR_VALUE\t");
    }
    for (uint32_t orig_pheno_idx = 0; orig_pheno_idx != orig_pheno_ct; ++orig_pheno_idx) {
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto PhenoSvd_ret_WRITE_FAIL;
      }
      write_iter = strcpyax(write_iter, &(orig_pheno_names[orig_pheno_idx * orig_max_pheno_name_blen]), '\t');
    }
    --write_iter;
    AppendBinaryEoln(&write_iter);
    if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
      goto PhenoSvd_ret_WRITE_FAIL;
    }
    for (uint32_t new_pheno_idx = 0; new_pheno_idx != new_pheno_ct; ++new_pheno_idx) {
      if (col_id) {
        write_iter = strcpya_k(write_iter, "SVDPHENO");
        write_iter = u32toa_x(new_pheno_idx + 1, '\t', write_iter);
      }
      if (col_sv) {
        write_iter = dtoa_g(singular_values[new_pheno_idx], write_iter);
        *write_iter++ = '\t';
      }
      const double* matrix_col_iter = &(matrix[new_pheno_idx]);
      for (uint32_t orig_pheno_idx = 0; orig_pheno_idx != orig_pheno_ct; ++orig_pheno_idx) {
        write_iter = dtoa_g(*matrix_col_iter, write_iter);
        *write_iter++ = '\t';
        matrix_col_iter = &(matrix_col_iter[svd_dim]);
      }
      --write_iter;
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto PhenoSvd_ret_WRITE_FAIL;
      }
    }

    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto PhenoSvd_ret_WRITE_FAIL;
    }
    outname_end[strlen(".svd.pheno")] = '\0';
    logprintfww("--pheno-svd: Results written to %s + %s_wts .\n", outname, outname);

    // ensure cleanup works if new phenotype initialization fails in the middle
    for (uint32_t pheno_idx = new_pheno_ct; pheno_idx != orig_pheno_ct; ++pheno_idx) {
      pheno_cols[pheno_idx].nonmiss = nullptr;
    }

    const uintptr_t old_pheno_names_size = orig_pheno_ct * orig_max_pheno_name_blen;
    const uintptr_t new_max_pheno_name_blen = 9 + UintSlen(new_pheno_ct + 1);
    const uintptr_t new_pheno_names_size = new_pheno_ct * new_max_pheno_name_blen;
    char* new_pheno_names = orig_pheno_names;
    if (new_pheno_names_size > old_pheno_names_size) {
      if (unlikely(pgl_malloc(new_pheno_names_size, &new_pheno_names))) {
        goto PhenoSvd_ret_NOMEM;
      }
      *pheno_names_ptr = new_pheno_names;
    }
    for (uint32_t pheno_idx = 0; pheno_idx != new_pheno_ct; ) {
      char* pheno_write_iter = strcpya_k(&(new_pheno_names[pheno_idx * new_max_pheno_name_blen]), "SVDPHENO");
      ++pheno_idx;
      write_iter = u32toa(pheno_idx, pheno_write_iter);
      *write_iter = '\0';
    }
    const uintptr_t nonmiss_vec_ct = BitCtToVecCt(raw_sample_ct);
    const uintptr_t data_vec_ct = DblCtToVecCt(raw_sample_ct);
    for (uint32_t pheno_idx = 0; pheno_idx != new_pheno_ct; ++pheno_idx) {
      // could try to reuse, but let's keep this a bit simpler for now
      vecaligned_free_cond(pheno_cols[pheno_idx].nonmiss);
      pheno_cols[pheno_idx].nonmiss = nullptr;
      uintptr_t* new_pheno_data_iter;
      if (unlikely(vecaligned_malloc((nonmiss_vec_ct + data_vec_ct) * kBytesPerVec, &new_pheno_data_iter))) {
        goto PhenoSvd_ret_NOMEM;
      }
      pheno_cols[pheno_idx].nonmiss = new_pheno_data_iter;
      memcpy(new_pheno_data_iter, sample_intersect, raw_sample_ctaw * sizeof(intptr_t));
      new_pheno_data_iter = &(new_pheno_data_iter[raw_sample_ctaw]);
      double* qt = R_CAST(double*, new_pheno_data_iter);
      pheno_cols[pheno_idx].data.qt = qt;

      const double* svd_sample_result_col = &(svd_sample_result[pheno_idx]);
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = sample_intersect[0];
      for (uint32_t sample_idx = 0; sample_idx != new_sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_intersect, &sample_uidx_base, &cur_bits);
        qt[sample_uidx] = svd_sample_result_col[sample_idx * svd_dim];
      }
    }
    *pheno_ct_ptr = new_pheno_ct;
    *max_pheno_name_blen_ptr = new_max_pheno_name_blen;
  }
  while (0) {
  PhenoSvd_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  PhenoSvd_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  PhenoSvd_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  PhenoSvd_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  PhenoSvd_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
  BLAS_SET_NUM_THREADS(1);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}
#endif

#ifdef __cplusplus
}  // namespace plink2
#endif
