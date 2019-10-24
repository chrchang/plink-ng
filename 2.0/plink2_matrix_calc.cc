// This file is part of PLINK 2.00, copyright (C) 2005-2019 Shaun Purcell,
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
#include "plink2_matrix.h"
#include "plink2_matrix_calc.h"
#include "plink2_random.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitScore(ScoreInfo* score_info_ptr) {
  score_info_ptr->flags = kfScore0;
  score_info_ptr->varid_col_p1 = 1;
  score_info_ptr->allele_col_p1 = 0;
  score_info_ptr->input_fname = nullptr;
  InitRangeList(&(score_info_ptr->input_col_idx_range_list));
}

void CleanupScore(ScoreInfo* score_info_ptr) {
  free_cond(score_info_ptr->input_fname);
  CleanupRangeList(&(score_info_ptr->input_col_idx_range_list));
}


uint32_t TriangleDivide(int64_t cur_prod_x2, int32_t modif) {
  // return smallest integer vv for which (vv * (vv + modif)) is no smaller
  // than cur_prod_x2, and neither term in the product is negative.
  int64_t vv;
  if (cur_prod_x2 == 0) {
    if (modif < 0) {
      return -modif;
    }
    return 0;
  }
  vv = S_CAST(int64_t, sqrt(S_CAST(double, cur_prod_x2)));
  while ((vv - 1) * (vv + modif - 1) >= cur_prod_x2) {
    vv--;
  }
  while (vv * (vv + modif) < cur_prod_x2) {
    vv++;
  }
  return vv;
}

void ParallelBounds(uint32_t ct, int32_t start, uint32_t parallel_idx, uint32_t parallel_tot, int32_t* __restrict bound_start_ptr, int32_t* __restrict bound_end_ptr) {
  int32_t modif = 1 - start * 2;
  int64_t ct_tot = S_CAST(int64_t, ct) * (ct + modif);
  *bound_start_ptr = TriangleDivide((ct_tot * parallel_idx) / parallel_tot, modif);
  *bound_end_ptr = TriangleDivide((ct_tot * (parallel_idx + 1)) / parallel_tot, modif);
}

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
    if (unlikely(
            bigstack_calloc_w(orig_sample_ctl, &sample_include_collapsed_nz) ||
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

PglErr KingCutoffBatch(const SampleIdInfo* siip, uint32_t raw_sample_ct, double king_cutoff, uintptr_t* sample_include, char* king_cutoff_fprefix, uint32_t* sample_ct_ptr) {
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
    if (unlikely(
            bigstack_calloc_w(sample_ct * orig_sample_ctl, &kinship_table) ||
            bigstack_alloc_u32(raw_sample_ct, &sample_uidx_to_king_uidx))) {
      goto KingCutoffBatch_ret_NOMEM;
    }

    snprintf(fprefix_end, 9, ".king.id");
    reterr = InitTextStream(king_cutoff_fprefix, kTextStreamBlenFast, 1, &txs);
    if (unlikely(reterr)) {
      goto KingCutoffBatch_ret_TSTREAM_FAIL;
    }
    // bugfix (18 Aug 2018): this missed some xid_mode possibilities
    // todo: try to simplify this interface, it's bordering on incomprehensible
    char* line_start;
    XidMode xid_mode;
    reterr = LoadXidHeader("king-cutoff", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeader0 : kfXidHeaderIgnoreSid, &line_idx, &txs, &xid_mode, &line_start, nullptr);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: Empty --king-cutoff ID file.\n");
        goto KingCutoffBatch_ret_MALFORMED_INPUT;
      }
      goto KingCutoffBatch_ret_TSTREAM_XID_FAIL;
    }

    uint32_t* xid_map;  // IDs not collapsed
    char* sorted_xidbox;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, sample_ct, 0, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto KingCutoffBatch_ret_1;
    }
    char* idbuf;
    if (unlikely(bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto KingCutoffBatch_ret_NOMEM;
    }
    SetAllU32Arr(raw_sample_ct, sample_uidx_to_king_uidx);
    if (*line_start == '#') {
      *line_start = '\0';
    }
    uintptr_t king_id_ct = 0;
    while (1) {
      if (!IsEolnKns(*line_start)) {
        const char* linebuf_iter = line_start;
        uint32_t sample_uidx;
        if (!SortedXidboxReadFind(sorted_xidbox, xid_map, max_xid_blen, sample_ct, 0, xid_mode, &linebuf_iter, &sample_uidx, idbuf)) {
          if (unlikely(sample_uidx_to_king_uidx[sample_uidx] != UINT32_MAX)) {
            char* first_tab = AdvToDelim(idbuf, '\t');
            char* second_tab = strchr(&(first_tab[1]), '\t');
            *first_tab = ' ';
            if (second_tab) {
              *second_tab = ' ';
            }
            snprintf(g_logbuf, kLogbufSize, "Error: Duplicate ID '%s' in %s .\n", idbuf, king_cutoff_fprefix);
            goto KingCutoffBatch_ret_MALFORMED_INPUT_WW;
          }
          sample_uidx_to_king_uidx[sample_uidx] = king_id_ct;
          ++king_id_ct;
        } else {
          if (unlikely(!linebuf_iter)) {
            goto KingCutoffBatch_ret_MISSING_TOKENS;
          }
        }
      }
      ++line_idx;
      reterr = TextNextLineLstrip(&txs, &line_start);
      if (reterr) {
        if (likely(reterr == kPglRetEof)) {
          reterr = kPglRetSuccess;
          break;
        }
        goto KingCutoffBatch_ret_TSTREAM_FAIL;
      }
    }

    BigstackReset(TextStreamMemStart(&txs));
    if (CleanupTextStream2(king_cutoff_fprefix, &txs, &reterr)) {
      goto KingCutoffBatch_ret_1;
    }
    uintptr_t* king_include;
    uint32_t* king_uidx_to_sample_idx;
    if (unlikely(
            bigstack_calloc_w(BitCtToWordCt(king_id_ct), &king_include) ||
            bigstack_alloc_u32(king_id_ct, &king_uidx_to_sample_idx))) {
      goto KingCutoffBatch_ret_NOMEM;
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
      goto KingCutoffBatch_ret_OPEN_FAIL;
    }
    if (unlikely(fseeko(binfile, 0, SEEK_END))) {
      goto KingCutoffBatch_ret_READ_FAIL;
    }
    const uint64_t fsize = ftello(binfile);
    const uint64_t fsize_double_expected = (king_id_ct * (S_CAST(uint64_t, king_id_ct) - 1) * (sizeof(double) / 2));
    const uint32_t is_double = (fsize == fsize_double_expected);
    rewind(binfile);
    const uint32_t first_king_uidx = AdvBoundedTo1Bit(king_include, 0, king_id_ct);
    uintptr_t king_uidx = AdvBoundedTo1Bit(king_include, first_king_uidx + 1, king_id_ct);
    if (king_uidx > 1) {
      if (fseeko(binfile, king_uidx * (S_CAST(uint64_t, king_uidx) - 1) * (2 + (2 * is_double)), SEEK_SET)) {
        goto KingCutoffBatch_ret_READ_FAIL;
      }
    }
    uintptr_t constraint_ct = 0;
    if (is_double) {
      // fread limit
      assert(king_id_ct <= ((kMaxBytesPerIO / sizeof(double)) + 1));
      double* king_drow;
      if (unlikely(bigstack_alloc_d(king_id_ct - 1, &king_drow))) {
        goto KingCutoffBatch_ret_NOMEM;
      }
      for (uint32_t king_idx = 1; king_uidx != king_id_ct; ++king_idx, ++king_uidx) {
        if (!IsSet(king_include, king_uidx)) {
          king_uidx = AdvBoundedTo1Bit(king_include, king_uidx + 1, king_id_ct);
          if (king_uidx == king_id_ct) {
            break;
          }
          if (unlikely(fseeko(binfile, S_CAST(uint64_t, king_uidx) * (king_uidx - 1) * (sizeof(double) / 2), SEEK_SET))) {
            goto KingCutoffBatch_ret_READ_FAIL;
          }
        }
        if (unlikely(!fread_unlocked(king_drow, king_uidx * sizeof(double), 1, binfile))) {
          goto KingCutoffBatch_ret_READ_FAIL;
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
        logerrprintfww("Error: Invalid --king-cutoff .bin file size (expected %" PRIu64 " or %" PRIu64 " bytes).\n", fsize_double_expected / 2, fsize_double_expected);
        goto KingCutoffBatch_ret_MALFORMED_INPUT;
      }
      assert(king_id_ct <= ((0x7ffff000 / sizeof(float)) + 1));
      const float king_cutoff_f = S_CAST(float, king_cutoff);
      float* king_frow;
      if (unlikely(bigstack_alloc_f(king_id_ct - 1, &king_frow))) {
        goto KingCutoffBatch_ret_NOMEM;
      }
      for (uint32_t king_idx = 1; king_uidx != king_id_ct; ++king_idx, ++king_uidx) {
        if (!IsSet(king_include, king_uidx)) {
          king_uidx = AdvBoundedTo1Bit(king_include, king_uidx + 1, king_id_ct);
          if (king_uidx == king_id_ct) {
            break;
          }
          if (unlikely(fseeko(binfile, S_CAST(uint64_t, king_uidx) * (king_uidx - 1) * (sizeof(float) / 2), SEEK_SET))) {
            goto KingCutoffBatch_ret_READ_FAIL;
          }
        }
        if (unlikely(!fread_unlocked(king_frow, king_uidx * sizeof(float), 1, binfile))) {
          goto KingCutoffBatch_ret_READ_FAIL;
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
      goto KingCutoffBatch_ret_NOMEM;
    }
  }
  while (0) {
  KingCutoffBatch_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  KingCutoffBatch_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  KingCutoffBatch_ret_READ_FAIL:
    logerrprintfww(kErrprintfFread, king_cutoff_fprefix, strerror(errno));
    reterr = kPglRetReadFail;
    break;
  KingCutoffBatch_ret_MISSING_TOKENS:
    logerrprintfww("Error: Fewer tokens than expected on line %" PRIuPTR " of %s .\n", line_idx, king_cutoff_fprefix);
    reterr = kPglRetMalformedInput;
    break;
  KingCutoffBatch_ret_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&txs)) {
      break;
    }
  KingCutoffBatch_ret_TSTREAM_FAIL:
    TextStreamErrPrint(king_cutoff_fprefix, &txs);
    break;
  KingCutoffBatch_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  KingCutoffBatch_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 KingCutoffBatch_ret_1:
  fclose_cond(binfile);
  if (CleanupTextStream(&txs, &reterr)) {
    snprintf(fprefix_end, 9, ".king.id");
    logerrprintfww(kErrprintfFread, king_cutoff_fprefix, strerror(errno));
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

// multithread globals
static uintptr_t* g_smaj_hom[2] = {nullptr, nullptr};
static uintptr_t* g_smaj_ref2het[2] = {nullptr, nullptr};
static uint32_t* g_thread_start = nullptr;
static uint32_t* g_king_counts = nullptr;
static uint32_t* g_loaded_sample_idx_pairs = nullptr;
static uint32_t g_homhom_needed = 0;

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
      *king_counts_iter++ += acc_ibs0;
      *king_counts_iter++ += acc_hethet;
      *king_counts_iter++ += acc_het2hom1;
      *king_counts_iter++ += acc_het1hom2;

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
      *king_counts_iter++ += acc_ibs0;
      *king_counts_iter++ += acc_hethet;
      *king_counts_iter++ += acc_het2hom1;
      *king_counts_iter++ += acc_het1hom2;
      *king_counts_iter++ += acc_homhom;

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
      *king_counts_iter++ += UniVecHsum16(acc_ibs0);
      *king_counts_iter++ += UniVecHsum16(acc_hethet);
      *king_counts_iter++ += UniVecHsum16(acc_het2hom1);
      *king_counts_iter++ += UniVecHsum16(acc_het1hom2);

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
      *king_counts_iter++ += UniVecHsum16(acc_ibs0);
      *king_counts_iter++ += UniVecHsum16(acc_hethet);
      *king_counts_iter++ += UniVecHsum16(acc_het2hom1);
      *king_counts_iter++ += UniVecHsum16(acc_het1hom2);
      *king_counts_iter++ += UniVecHsum16(acc_homhom);

      first_hom_iter = &(first_hom_iter[kKingMultiplexVecs]);
      first_ref2het_iter = &(first_ref2het_iter[kKingMultiplexVecs]);
    }
  }
}
#endif
static_assert(!(kKingMultiplexWords % 2), "kKingMultiplexWords must be even for safe bit-transpose.");

THREAD_FUNC_DECL CalcKingThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uint64_t mem_start_idx = g_thread_start[0];
  const uint64_t start_idx = g_thread_start[tidx];
  const uint32_t end_idx = g_thread_start[tidx + 1];
  const uint32_t homhom_needed = g_homhom_needed;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    if (homhom_needed) {
      IncrKingHomhom(g_smaj_hom[parity], g_smaj_ref2het[parity], start_idx, end_idx, &(g_king_counts[((start_idx * (start_idx - 1) - mem_start_idx * (mem_start_idx - 1)) / 2) * 5]));
    } else {
      IncrKing(g_smaj_hom[parity], g_smaj_ref2het[parity], start_idx, end_idx, &(g_king_counts[(start_idx * (start_idx - 1) - mem_start_idx * (mem_start_idx - 1)) * 2]));
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

/*
double ComputeKinship(const uint32_t* king_counts_entry) {
  const uint32_t ibs0_ct = king_counts_entry[0];
  const uint32_t hethet_ct = king_counts_entry[1];
  const uint32_t het2hom1_ct = king_counts_entry[2];
  const uint32_t het1hom2_ct = king_counts_entry[3];
  // const uint32_t homhom_ct = king_counts_entry[4];
  const intptr_t smaller_het_ct = hethet_ct + MINV(het1hom2_ct, het2hom1_ct);
  return 0.5 - (S_CAST(double, 4 * S_CAST(intptr_t, ibs0_ct) + het1hom2_ct + het2hom1_ct) / S_CAST(double, 4 * smaller_het_ct));
}
*/

// '2' refers to the larger index here
double ComputeKinship(const uint32_t* king_counts_entry, uint32_t singleton_het1_ct, uint32_t singleton_hom1_ct, uint32_t singleton_het2_ct, uint32_t singleton_hom2_ct) {
  const uint32_t ibs0_ct = king_counts_entry[0] + singleton_hom1_ct + singleton_hom2_ct;
  const uint32_t hethet_ct = king_counts_entry[1];
  const uint32_t het2hom1_ct = king_counts_entry[2] + singleton_het2_ct;
  const uint32_t het1hom2_ct = king_counts_entry[3] + singleton_het1_ct;
  // const uint32_t homhom_ct = king_counts_entry[4];
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
    cswritep = strcpya_k(cswritep, "ID1\t");
    if (king_col_sid) {
      cswritep = strcpya_k(cswritep, "SID1\t");
    }
    if (king_col_fid) {
      // bugfix (15 Feb 2018): yesterday's build had \n instead of \t here
      cswritep = strcpya_k(cswritep, "FID2\t");
    }
    cswritep = strcpya_k(cswritep, "ID2\t");
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
  if (king_flags & kfKingColKinship) {
    cswritep = strcpya_k(cswritep, "KINSHIP\t");
  }
  DecrAppendBinaryEoln(&cswritep);
  return cswritep;
}

PglErr CalcKing(const SampleIdInfo* siip, const uintptr_t* variant_include_orig, const ChrInfo* cip, uint32_t raw_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, double king_cutoff, double king_table_filter, KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, uintptr_t* sample_include, uint32_t* sample_ct_ptr, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  char* cswritetp = nullptr;
  CompressStreamState css;
  CompressStreamState csst;
  ThreadsState ts;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  PreinitCstream(&csst);
  InitThreads3z(&ts);
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
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* variant_include;
    uint32_t* singleton_het_cts;
    uint32_t* singleton_hom_cts;
    uint32_t* singleton_missing_cts;
    if (unlikely(
            bigstack_alloc_w(raw_variant_ctl, &variant_include) ||
            bigstack_calloc_u32(sample_ct, &singleton_het_cts) ||
            bigstack_calloc_u32(sample_ct, &singleton_hom_cts) ||
            bigstack_calloc_u32(sample_ct, &singleton_missing_cts))) {
      goto CalcKing_ret_NOMEM;
    }
    memcpy(variant_include, variant_include_orig, raw_variant_ctl * sizeof(intptr_t));
    const uint32_t non_autosomal_variant_ct = CountNonAutosomalVariants(variant_include, cip, 1, 1);
    if (non_autosomal_variant_ct) {
      logprintf("Excluding %u variant%s on non-autosomes from KING-robust calculation.\n", non_autosomal_variant_ct, (non_autosomal_variant_ct == 1)? "" : "s");
      variant_ct -= non_autosomal_variant_ct;
      if (!variant_ct) {
        logerrprintf("Error: No variants remaining for KING-robust calculation.\n");
        goto CalcKing_ret_DEGENERATE_DATA;
      }
      ExcludeNonAutosomalVariants(cip, variant_include);
    }

    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > sample_ct / 32) {
      calc_thread_ct = sample_ct / 32;
    }
    if (!calc_thread_ct) {
      calc_thread_ct = 1;
    }
    // possible todo: allow this to change between passes
    ts.calc_thread_ct = calc_thread_ct;
    if (unlikely(
            bigstack_alloc_u32(calc_thread_ct + 1, &g_thread_start) ||
            bigstack_alloc_thread(calc_thread_ct, &ts.threads))) {
      goto CalcKing_ret_NOMEM;
    }
    const uintptr_t sample_ctl = BitCtToWordCt(sample_ct);
    uintptr_t* kinship_table = nullptr;
    if (king_cutoff != -1) {
      if (unlikely(bigstack_calloc_w(sample_ct * sample_ctl, &kinship_table))) {
        goto CalcKing_ret_NOMEM;
      }
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
    const uint32_t sample_ctaw2 = QuaterCtToAlignedWordCt(sample_ct);
    g_homhom_needed = (king_flags & kfKingColNsnp) || ((!(king_flags & kfKingCounts)) && (king_flags & (kfKingColHethet | kfKingColIbs0 | kfKingColIbs1)));
    uint32_t grand_row_start_idx;
    uint32_t grand_row_end_idx;
    ParallelBounds(sample_ct, 1, parallel_idx, parallel_tot, R_CAST(int32_t*, &grand_row_start_idx), R_CAST(int32_t*, &grand_row_end_idx));
    const uint32_t king_bufsizew = kKingMultiplexWords * grand_row_end_idx;
    const uint32_t homhom_needed_p4 = g_homhom_needed + 4;
    uintptr_t* cur_sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* loadbuf;
    uintptr_t* splitbuf_hom;
    uintptr_t* splitbuf_ref2het;
    VecW* vecaligned_buf;
    if (unlikely(
            bigstack_alloc_w(raw_sample_ctl, &cur_sample_include) ||
            bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
            bigstack_alloc_w(sample_ctaw2, &loadbuf) ||
            bigstack_alloc_w(kPglBitTransposeBatch * sample_ctaw, &splitbuf_hom) ||
            bigstack_alloc_w(kPglBitTransposeBatch * sample_ctaw, &splitbuf_ref2het) ||
            bigstack_alloc_w(king_bufsizew, &(g_smaj_hom[0])) ||
            bigstack_alloc_w(king_bufsizew, &(g_smaj_ref2het[0])) ||
            bigstack_alloc_w(king_bufsizew, &(g_smaj_hom[1])) ||
            bigstack_alloc_w(king_bufsizew, &(g_smaj_ref2het[1])) ||
            bigstack_alloc_v(kPglBitTransposeBufvecs, &vecaligned_buf))) {
      goto CalcKing_ret_NOMEM;
    }
    uint32_t sparse_variant_ct = 0;
    // Count singleton/monomorphic variants, record all relevant statistics,
    // and remove them from the main loop.
    {
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      PgrClearLdCache(simple_pgrp);
      logprintf("%s: Scanning for singletons and monomorphic variants...", flagname);
      fflush(stdout);
      const uint32_t sample_ct_m1 = sample_ct - 1;
      uint32_t skip_ct = 0;
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      for (uint32_t uii = 0; uii != variant_ct; ++uii) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        reterr = PgrGet(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, loadbuf);
        if (unlikely(reterr)) {
          goto CalcKing_ret_PGR_FAIL;
        }
        ZeroTrailingQuaters(sample_ct, loadbuf);
        STD_ARRAY_DECL(uint32_t, 4, genocounts);
        GenovecCountFreqsUnsafe(loadbuf, sample_ct, genocounts);
        if (genocounts[0] >= sample_ct_m1) {
          if (genocounts[0] == sample_ct_m1) {
            for (uint32_t widx = 0; ; ++widx) {
              uintptr_t geno_word = loadbuf[widx];
              if (geno_word) {
                const uint32_t sample_idx_lowbits = ctzw(geno_word) / 2;
                // no need to mask out high bits, they're guaranteed to be zero
                geno_word = geno_word >> (sample_idx_lowbits * 2);
                const uint32_t sample_idx = widx * kBitsPerWordD2 + sample_idx_lowbits;
                if (geno_word == 1) {
                  singleton_het_cts[sample_idx] += 1;
                } else if (geno_word == 2) {
                  singleton_hom_cts[sample_idx] += 1;
                } else {
                  singleton_missing_cts[sample_idx] += 1;
                }
                break;
              }
            }
          }
        } else if (genocounts[2] >= sample_ct_m1) {
          if (genocounts[2] == sample_ct_m1) {
            for (uint32_t widx = 0; ; ++widx) {
              uintptr_t xor_word = loadbuf[widx] ^ kMaskAAAA;
              if (xor_word) {
                const uint32_t sample_idx_lowbits = ctzw(xor_word) / 2;
                xor_word = (xor_word >> (sample_idx_lowbits * 2)) & 3;
                const uint32_t sample_idx = widx * kBitsPerWordD2 + sample_idx_lowbits;
                if (xor_word == 3) {
                  singleton_het_cts[sample_idx] += 1;
                } else if (xor_word == 2) {
                  singleton_hom_cts[sample_idx] += 1;
                } else {
                  singleton_missing_cts[sample_idx] += 1;
                }
                break;
              }
            }
          }
        } else if ((genocounts[3] >= sample_ct_m1) || ((!g_homhom_needed) && ((genocounts[0] + genocounts[3] == sample_ct) || (genocounts[2] + genocounts[3] == sample_ct)))) {
          ++skip_ct;
        } else {
          continue;
        }
        ++sparse_variant_ct;
        ClearBit(variant_uidx, variant_include);
      }
      logputs(" done.\n");
      if (!sparse_variant_ct) {
        logputs("No singletons/monomorphics found.\n");
      } else {
        logprintf("%u variant%s handled by initial scan.\n", sparse_variant_ct, (sparse_variant_ct == 1)? "" : "s");
      }
      variant_ct -= sparse_variant_ct;
      sparse_variant_ct -= skip_ct;
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
    g_king_counts = R_CAST(uint32_t*, g_bigstack_base);
    for (uint32_t pass_idx_p1 = 1; pass_idx_p1 <= pass_ct; ++pass_idx_p1) {
      const uint32_t row_start_idx = row_end_idx;
      row_end_idx = NextTrianglePass(row_start_idx, grand_row_end_idx, 1, cells_avail);
      TriangleLoadBalance(calc_thread_ct, row_start_idx, row_end_idx, 1, g_thread_start);
      const uintptr_t tot_cells = (S_CAST(uint64_t, row_end_idx) * (row_end_idx - 1) - S_CAST(uint64_t, row_start_idx) * (row_start_idx - 1)) / 2;
      ZeroU32Arr(tot_cells * homhom_needed_p4, g_king_counts);

      if (variant_ct) {
        // possible todo: doubleton optimization
        const uint32_t row_end_idxaw = BitCtToAlignedWordCt(row_end_idx);
        const uint32_t row_end_idxaw2 = QuaterCtToAlignedWordCt(row_end_idx);
        if (row_end_idxaw % 2) {
          const uint32_t cur_king_bufsizew = kKingMultiplexWords * row_end_idx;
          uintptr_t* smaj_hom0_last = &(g_smaj_hom[0][kKingMultiplexWords - 1]);
          uintptr_t* smaj_ref2het0_last = &(g_smaj_ref2het[0][kKingMultiplexWords - 1]);
          uintptr_t* smaj_hom1_last = &(g_smaj_hom[1][kKingMultiplexWords - 1]);
          uintptr_t* smaj_ref2het1_last = &(g_smaj_ref2het[1][kKingMultiplexWords - 1]);
          for (uint32_t offset = 0; offset < cur_king_bufsizew; offset += kKingMultiplexWords) {
            smaj_hom0_last[offset] = 0;
            smaj_ref2het0_last[offset] = 0;
            smaj_hom1_last[offset] = 0;
            smaj_ref2het1_last[offset] = 0;
          }
        }
        memcpy(cur_sample_include, sample_include, raw_sample_ctl * sizeof(intptr_t));
        if (row_end_idx != grand_row_end_idx) {
          uint32_t sample_uidx_end = IdxToUidxBasic(sample_include, row_end_idx);
          ClearBitsNz(sample_uidx_end, raw_sample_ct, cur_sample_include);
        }
        FillCumulativePopcounts(cur_sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
        if (pass_idx_p1 != 1) {
          ReinitThreads3z(&ts);
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
        PgrClearLdCache(simple_pgrp);
        do {
          const uint32_t cur_block_size = MINV(variant_ct - variants_completed, kKingMultiplex);
          uintptr_t* cur_smaj_hom = g_smaj_hom[parity];
          uintptr_t* cur_smaj_ref2het = g_smaj_ref2het[parity];
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
              reterr = PgrGet(cur_sample_include, sample_include_cumulative_popcounts, row_end_idx, variant_uidx, simple_pgrp, loadbuf);
              if (unlikely(reterr)) {
                goto CalcKing_ret_PGR_FAIL;
              }
              SetTrailingQuaters(row_end_idx, loadbuf);
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
            JoinThreads3z(&ts);
            // CalcKingThread() never errors out
          } else {
            ts.thread_func_ptr = CalcKingThread;
          }
          // this update must occur after JoinThreads3z() call
          ts.is_last_block = (variants_completed + cur_block_size == variant_ct);
          if (unlikely(SpawnThreads3z(variants_completed, &ts))) {
            goto CalcKing_ret_THREAD_CREATE_FAIL;
          }
          printf("\r%s pass %u/%u: %u variants complete.", flagname, pass_idx_p1, pass_ct, variants_completed);
          fflush(stdout);
          variants_completed += cur_block_size;
          parity = 1 - parity;
        } while (!ts.is_last_block);
        JoinThreads3z(&ts);
      }
      if (matrix_shape || (king_flags & kfKingColAll)) {
        printf("\r%s pass %u/%u: Writing...                   \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", flagname, pass_idx_p1, pass_ct);
        fflush(stdout);
        // allow simultaneous --make-king + --make-king-table
        if (matrix_shape) {
          if (!(king_flags & (kfKingMatrixBin | kfKingMatrixBin4))) {
            const uint32_t is_squarex = king_flags & (kfKingMatrixSq | kfKingMatrixSq0);
            const uint32_t is_square0 = king_flags & kfKingMatrixSq0;
            uint32_t* results_iter = g_king_counts;
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
#ifdef __arm__
#  error "Unaligned accesses in CalcKing()."
#endif
                  uintptr_t* writep_alias = R_CAST(uintptr_t*, cswritep);
                  for (uintptr_t widx = 0; widx != wct; ++widx) {
                    *writep_alias++ = tabzero_word;
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
            uint32_t* results_iter = g_king_counts;
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
          const uint32_t king_col_kinship = king_flags & kfKingColKinship;
          const uint32_t report_counts = king_flags & kfKingCounts;
          uint32_t* results_iter = g_king_counts;
          double nonmiss_recip = 0.0;
          for (uint32_t sample_idx1 = row_start_idx; sample_idx1 != row_end_idx; ++sample_idx1) {
            const char* sample_fmtid1 = &(collapsed_sample_fmtids[max_sample_fmtid_blen * sample_idx1]);
            const uint32_t singleton_het1_ct = singleton_het_cts[sample_idx1];
            const uint32_t singleton_hom1_ct = singleton_hom_cts[sample_idx1];
            const uint32_t sample_fmtid1_slen = strlen(sample_fmtid1);
            for (uint32_t sample_idx2 = 0; sample_idx2 != sample_idx1; ++sample_idx2, results_iter = &(results_iter[homhom_needed_p4])) {
              const uint32_t singleton_het2_ct = singleton_het_cts[sample_idx2];
              const uint32_t singleton_hom2_ct = singleton_hom_cts[sample_idx2];
              const uint32_t ibs0_ct = results_iter[0] + singleton_hom2_ct + singleton_hom1_ct;
              const uint32_t hethet_ct = results_iter[1];
              // '2' here refers to the larger index, so this is swapped
              const uint32_t het2hom1_ct = results_iter[2] + singleton_het1_ct;
              const uint32_t het1hom2_ct = results_iter[3] + singleton_het2_ct;
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
                const uint32_t homhom_ct = results_iter[4] + sparse_variant_ct - singleton_het2_ct - singleton_missing_cts[sample_idx2] - singleton_het1_ct - singleton_missing_cts[sample_idx1];
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
        uint32_t* results_iter = g_king_counts;
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
    }
    putc_unlocked('\n', stdout);
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
  CleanupThreads3z(&ts, nullptr);
  CswriteCloseCond(&csst, cswritetp);
  CswriteCloseCond(&css, cswritep);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
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

THREAD_FUNC_DECL CalcKingTableSubsetThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uint32_t start_idx = g_thread_start[tidx];
  const uint32_t end_idx = g_thread_start[tidx + 1];
  const uint32_t homhom_needed = g_homhom_needed;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    if (homhom_needed) {
      IncrKingSubsetHomhom(g_loaded_sample_idx_pairs, g_smaj_hom[parity], g_smaj_ref2het[parity], start_idx, end_idx, g_king_counts);
    } else {
      IncrKingSubset(g_loaded_sample_idx_pairs, g_smaj_hom[parity], g_smaj_ref2het[parity], start_idx, end_idx, g_king_counts);
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

PglErr KingTableSubsetLoad(const char* sorted_xidbox, const uint32_t* xid_map, uintptr_t max_xid_blen, uintptr_t orig_sample_ct, double king_table_subset_thresh, XidMode xid_mode, uint32_t skip_sid, uint32_t kinship_skip, uint32_t is_first_parallel_scan, uint64_t pair_idx_start, uint64_t pair_idx_stop, uintptr_t line_idx, TextStream* txsp, uint64_t* pair_idx_ptr, uint32_t* loaded_sample_idx_pairs, char* idbuf) {
  PglErr reterr = kPglRetSuccess;
  {
    uint64_t pair_idx = *pair_idx_ptr;
    // Assumes header line already read if pair_idx == 0, and if pair_idx is
    // positive, we're that far into the file.
    uint32_t* loaded_sample_idx_pairs_iter = loaded_sample_idx_pairs;
    for (char* line_iter = TextLineEnd(txsp); ; line_iter = AdvPastDelim(line_iter, '\n')) {
      reterr = TextNextNonemptyLineLstripUnsafe(txsp, &line_idx, &line_iter);
      if (reterr) {
        if (likely(reterr == kPglRetEof)) {
          reterr = kPglRetSuccess;
          break;
        }
        goto KingTableSubsetLoad_ret_TSTREAM_FAIL;
      }
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
      if (skip_sid) {
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto KingTableSubsetLoad_ret_MISSING_TOKENS;
        }
        linebuf_iter = FirstNonTspace(CurTokenEnd(linebuf_iter));
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
      if (king_table_subset_thresh != -DBL_MAX) {
        linebuf_iter = FirstNonTspace(linebuf_iter);
        linebuf_iter = NextTokenMult0(linebuf_iter, kinship_skip);
        if (unlikely(!linebuf_iter)) {
          goto KingTableSubsetLoad_ret_MISSING_TOKENS;
        }
        double cur_kinship;
        if ((!ScanadvDouble(linebuf_iter, &cur_kinship)) || (cur_kinship < king_table_subset_thresh)) {
          line_iter = K_CAST(char*, linebuf_iter);
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
          break;
        }
        // large --parallel job, first pass: count number of valid pairs, don't
        // save the remainder
        pair_idx_start = ~0LLU;
      }
    }
    *pair_idx_ptr = pair_idx;
  }
  while (0) {
  KingTableSubsetLoad_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--king-table-subset file", txsp);
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

PglErr CalcKingTableSubset(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const char* subset_fname, uint32_t raw_sample_ct, uint32_t orig_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, double king_table_filter, double king_table_subset_thresh, KingFlags king_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  CompressStreamState css;
  ThreadsState ts;
  PreinitTextStream(&txs);
  PreinitCstream(&css);
  InitThreads3z(&ts);
  {
    if (unlikely(IsSet(cip->haploid_mask, 0))) {
      logerrputs("Error: --make-king-table cannot be used on haploid genomes.\n");
      goto CalcKingTableSubset_ret_INCONSISTENT_INPUT;
    }
    reterr = ConditionalAllocateNonAutosomalVariants(cip, "--make-king-table", raw_variant_ct, &variant_include, &variant_ct);
    if (unlikely(reterr)) {
      goto CalcKingTableSubset_ret_1;
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
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uint32_t sample_ctaw = BitCtToAlignedWordCt(orig_sample_ct);
    uint32_t sample_ctaw2 = QuaterCtToAlignedWordCt(orig_sample_ct);
    uint32_t king_bufsizew = kKingMultiplexWords * orig_sample_ct;
    uintptr_t* cur_sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* loadbuf;
    uintptr_t* splitbuf_hom;
    uintptr_t* splitbuf_ref2het;
    VecW* vecaligned_buf;
    // ok if allocations are a bit oversized
    if (unlikely(
            bigstack_alloc_w(raw_sample_ctl, &cur_sample_include) ||
            bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
            bigstack_alloc_w(sample_ctaw2, &loadbuf) ||
            bigstack_alloc_w(kPglBitTransposeBatch * sample_ctaw, &splitbuf_hom) ||
            bigstack_alloc_w(kPglBitTransposeBatch * sample_ctaw, &splitbuf_ref2het) ||
            bigstack_alloc_w(king_bufsizew, &(g_smaj_hom[0])) ||
            bigstack_alloc_w(king_bufsizew, &(g_smaj_ref2het[0])) ||
            bigstack_alloc_w(king_bufsizew, &(g_smaj_hom[1])) ||
            bigstack_alloc_w(king_bufsizew, &(g_smaj_ref2het[1])) ||
            bigstack_alloc_v(kPglBitTransposeBufvecs, &vecaligned_buf))) {
      goto CalcKingTableSubset_ret_NOMEM;
    }
    SetKingTableFname(king_flags, parallel_idx, parallel_tot, outname_end);
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
        logerrputs("Error: Failed to append '~' to --king-table-subset input filename.\n");
        goto CalcKingTableSubset_ret_OPEN_FAIL;
      }
      subset_fname = g_textbuf;
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
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > orig_sample_ct / 32) {
      calc_thread_ct = orig_sample_ct / 32;
    }
    if (!calc_thread_ct) {
      calc_thread_ct = 1;
    }
    // possible todo: allow this to change between passes
    ts.calc_thread_ct = calc_thread_ct;
    // could eventually have 64-bit g_thread_start?
    if (unlikely(
            bigstack_alloc_u32(calc_thread_ct + 1, &g_thread_start) ||
            bigstack_alloc_thread(calc_thread_ct, &ts.threads))) {
      goto CalcKingTableSubset_ret_NOMEM;
    }

    reterr = InitTextStream(subset_fname, kTextStreamBlenFast, 1, &txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: Empty --king-table-subset file.\n");
        goto CalcKingTableSubset_ret_MALFORMED_INPUT;
      }
      goto CalcKingTableSubset_ret_TSTREAM_FAIL;
    }
    const char* linebuf_iter;
    uintptr_t line_idx = 0;
    reterr = TextNextNonemptyLineLstripK(&txs, &line_idx, &linebuf_iter);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
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
    } else {
      if (unlikely(*linebuf_iter != '#')) {
        goto CalcKingTableSubset_ret_INVALID_HEADER;
      }
      ++linebuf_iter;
      --token_slen;
    }
    if (unlikely(!strequal_k(linebuf_iter, "ID1", token_slen))) {
      goto CalcKingTableSubset_ret_INVALID_HEADER;
    }
    linebuf_iter = FirstNonTspace(token_end);
    token_end = CurTokenEnd(linebuf_iter);
    token_slen = token_end - linebuf_iter;
    uint32_t skip_sid = 0;
    XidMode xid_mode = fid_present? kfXidModeFidIid : kfXidModeIid;
    if (strequal_k(linebuf_iter, "SID1", token_slen)) {
      if (siip->sids) {
        xid_mode = fid_present? kfXidModeFidIidSid : kfXidModeIidSid;
      } else {
        skip_sid = 1;
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
    if (unlikely(!strequal_k(linebuf_iter, "ID2", token_slen))) {
      goto CalcKingTableSubset_ret_INVALID_HEADER;
    }
    if (xid_mode == kfXidModeFidIidSid) {
      // technically don't need to check this in skip_sid case
      linebuf_iter = FirstNonTspace(token_end);
      token_end = CurTokenEnd(linebuf_iter);
      token_slen = token_end - linebuf_iter;
      if (unlikely(!strequal_k(linebuf_iter, "SID2", token_slen))) {
        goto CalcKingTableSubset_ret_INVALID_HEADER;
      }
    }
    uint32_t kinship_skip = 0;
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

    uint32_t* xid_map;  // IDs not collapsed
    char* sorted_xidbox;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(orig_sample_include, siip, orig_sample_ct, 0, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto CalcKingTableSubset_ret_1;
    }
    char* idbuf;
    if (unlikely(bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto CalcKingTableSubset_ret_NOMEM;
    }

    g_homhom_needed = (king_flags & kfKingColNsnp) || ((!(king_flags & kfKingCounts)) && (king_flags & (kfKingColHethet | kfKingColIbs0 | kfKingColIbs1)));
    const uint32_t homhom_needed_p4 = g_homhom_needed + 4;
    // if homhom_needed, 8 + 20 bytes per pair, otherwise 8 + 16
    uintptr_t pair_buf_capacity = bigstack_left();
    if (unlikely(pair_buf_capacity < 2 * kCacheline)) {
      goto CalcKingTableSubset_ret_NOMEM;
    }
    // adverse rounding
    pair_buf_capacity = (pair_buf_capacity - 2 * kCacheline) / (24 + 4 * g_homhom_needed);
    if (pair_buf_capacity > 0xffffffffU) {
      // 32-bit g_thread_start[] for now
      pair_buf_capacity = 0xffffffffU;
    }
    g_loaded_sample_idx_pairs = S_CAST(uint32_t*, bigstack_alloc_raw_rd(pair_buf_capacity * 2 * sizeof(int32_t)));
    g_king_counts = R_CAST(uint32_t*, g_bigstack_base);
    uint64_t pair_idx = 0;
    fputs("Scanning --king-table-subset file...", stdout);
    fflush(stdout);
    reterr = KingTableSubsetLoad(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, king_table_subset_thresh, xid_mode, skip_sid, kinship_skip, (parallel_tot != 1), 0, pair_buf_capacity, line_idx, &txs, &pair_idx, g_loaded_sample_idx_pairs, idbuf);
    if (unlikely(reterr)) {
      goto CalcKingTableSubset_ret_1;
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
          logerrputs("Error: Too many --king-table-subset sample pairs for this " PROG_NAME_STR " build.\n");
          reterr = kPglRetNotYetSupported;
          goto CalcKingTableSubset_ret_1;
        }
        if (pair_idx_global_stop > pair_buf_capacity) {
          // large --parallel job
          reterr = TextRewind(&txs);
          if (unlikely(reterr)) {
            goto CalcKingTableSubset_ret_TSTREAM_FAIL;
          }
          // bugfix (4 Oct 2019): forgot a bunch of reinitialization here
          line_idx = 0;
          char* header_throwaway;
          reterr = TextNextNonemptyLineLstrip(&txs, &line_idx, &header_throwaway);
          if (unlikely(reterr)) {
            goto CalcKingTableSubset_ret_TSTREAM_REWIND_FAIL;
          }
          pair_idx = 0;
          reterr = KingTableSubsetLoad(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, king_table_subset_thresh, xid_mode, skip_sid, kinship_skip, 0, pair_idx_global_start, MINV(pair_idx_global_stop, pair_idx_global_start + pair_buf_capacity), line_idx, &txs, &pair_idx, g_loaded_sample_idx_pairs, idbuf);
          if (unlikely(reterr)) {
            goto CalcKingTableSubset_ret_1;
          }
        } else {
          pair_idx = pair_idx_global_stop;
          if (pair_idx_global_start) {
            memmove(g_loaded_sample_idx_pairs, &(g_loaded_sample_idx_pairs[pair_idx_global_start * 2]), (pair_idx_global_stop - pair_idx_global_start) * 2 * sizeof(int32_t));
          }
        }
      } else {
        pair_idx = pair_idx_global_stop;
        if (pair_idx_global_start) {
          memmove(g_loaded_sample_idx_pairs, &(g_loaded_sample_idx_pairs[pair_idx_global_start * 2]), (pair_idx_global_stop - pair_idx_global_start) * 2 * sizeof(int32_t));
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
        SetBit(g_loaded_sample_idx_pairs[ulii], cur_sample_include);
      }
      FillCumulativePopcounts(cur_sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      const uint32_t cur_sample_ct = sample_include_cumulative_popcounts[raw_sample_ctl - 1] + PopcountWord(cur_sample_include[raw_sample_ctl - 1]);
      const uint32_t cur_sample_ctaw = BitCtToAlignedWordCt(cur_sample_ct);
      const uint32_t cur_sample_ctaw2 = QuaterCtToAlignedWordCt(cur_sample_ct);
      if (cur_sample_ct != raw_sample_ct) {
        for (uintptr_t ulii = 0; ulii != cur_pair_ct_x2; ++ulii) {
          g_loaded_sample_idx_pairs[ulii] = RawToSubsettedPos(cur_sample_include, sample_include_cumulative_popcounts, g_loaded_sample_idx_pairs[ulii]);
        }
      }
      ZeroU32Arr(cur_pair_ct * homhom_needed_p4, g_king_counts);
      CollapsedSampleFmtidInit(cur_sample_include, siip, cur_sample_ct, king_col_fid, king_col_sid, max_sample_fmtid_blen, collapsed_sample_fmtids);
      for (uint32_t tidx = 0; tidx <= calc_thread_ct; ++tidx) {
        g_thread_start[tidx] = (tidx * S_CAST(uint64_t, cur_pair_ct)) / calc_thread_ct;
      }
      if (pass_idx != 1) {
        ReinitThreads3z(&ts);
      }
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      uint32_t variants_completed = 0;
      uint32_t parity = 0;
      const uint32_t sample_batch_ct_m1 = (cur_sample_ct - 1) / kPglBitTransposeBatch;
      PgrClearLdCache(simple_pgrp);
      do {
        const uint32_t cur_block_size = MINV(variant_ct - variants_completed, kKingMultiplex);
        uintptr_t* cur_smaj_hom = g_smaj_hom[parity];
        uintptr_t* cur_smaj_ref2het = g_smaj_ref2het[parity];
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
            reterr = PgrGet(cur_sample_include, sample_include_cumulative_popcounts, cur_sample_ct, variant_uidx, simple_pgrp, loadbuf);
            if (unlikely(reterr)) {
              goto CalcKingTableSubset_ret_PGR_FAIL;
            }
            // may want to support some sort of low-MAF optimization here
            SetTrailingQuaters(cur_sample_ct, loadbuf);
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
          JoinThreads3z(&ts);
          // CalcKingTableSubsetThread() never errors out
        } else {
          ts.thread_func_ptr = CalcKingTableSubsetThread;
        }
        // this update must occur after JoinThreads3z() call
        ts.is_last_block = (variants_completed + cur_block_size == variant_ct);
        if (unlikely(SpawnThreads3z(variants_completed, &ts))) {
          goto CalcKingTableSubset_ret_THREAD_CREATE_FAIL;
        }
        printf("\r--make-king-table pass %" PRIuPTR ": %u variants complete.", pass_idx, variants_completed);
        fflush(stdout);
        variants_completed += cur_block_size;
        parity = 1 - parity;
      } while (!ts.is_last_block);
      JoinThreads3z(&ts);
      printf("\r--make-king-table pass %" PRIuPTR ": Writing...                   \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", pass_idx);
      fflush(stdout);

      const uint32_t king_col_id = king_flags & kfKingColId;
      const uint32_t king_col_nsnp = king_flags & kfKingColNsnp;
      const uint32_t king_col_hethet = king_flags & kfKingColHethet;
      const uint32_t king_col_ibs0 = king_flags & kfKingColIbs0;
      const uint32_t king_col_ibs1 = king_flags & kfKingColIbs1;
      const uint32_t king_col_kinship = king_flags & kfKingColKinship;
      const uint32_t report_counts = king_flags & kfKingCounts;
      uint32_t* results_iter = g_king_counts;
      double nonmiss_recip = 0.0;
      for (uintptr_t cur_pair_idx = 0; cur_pair_idx != cur_pair_ct; ++cur_pair_idx, results_iter = &(results_iter[homhom_needed_p4])) {
        const uint32_t ibs0_ct = results_iter[0];
        const uint32_t hethet_ct = results_iter[1];
        const uint32_t het2hom1_ct = results_iter[2];
        const uint32_t het1hom2_ct = results_iter[3];
        const intptr_t smaller_het_ct = hethet_ct + MINV(het1hom2_ct, het2hom1_ct);
        const double kinship_coeff = 0.5 - (S_CAST(double, 4 * S_CAST(intptr_t, ibs0_ct) + het1hom2_ct + het2hom1_ct) / S_CAST(double, 4 * smaller_het_ct));
        if ((king_table_filter != -DBL_MAX) && (kinship_coeff < king_table_filter)) {
          ++king_table_filter_ct;
          continue;
        }
        const uint32_t sample_idx1 = g_loaded_sample_idx_pairs[2 * cur_pair_idx];
        const uint32_t sample_idx2 = g_loaded_sample_idx_pairs[2 * cur_pair_idx + 1];
        if (king_col_id) {
          cswritep = strcpyax(cswritep, &(collapsed_sample_fmtids[max_sample_fmtid_blen * sample_idx1]), '\t');
          cswritep = strcpyax(cswritep, &(collapsed_sample_fmtids[max_sample_fmtid_blen * sample_idx2]), '\t');
        }
        if (homhom_needed_p4 == 5) {
          const uint32_t homhom_ct = results_iter[4];
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
      fputs("Scanning --king-table-subset file...", stdout);
      fflush(stdout);
      reterr = KingTableSubsetLoad(sorted_xidbox, xid_map, max_xid_blen, orig_sample_ct, king_table_subset_thresh, xid_mode, skip_sid, kinship_skip, 0, pair_idx_cur_start, MINV(pair_idx_global_stop, pair_idx_cur_start + pair_buf_capacity), line_idx, &txs, &pair_idx, g_loaded_sample_idx_pairs, idbuf);
      if (unlikely(reterr)) {
        goto CalcKingTableSubset_ret_1;
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
  CalcKingTableSubset_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
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
  CleanupThreads3z(&ts, nullptr);
  CleanupTextStream2("--king-table-subset file", &txs, &reterr);
  CswriteCloseCond(&css, cswritep);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

// assumes trailing bits of genovec are zeroed out
PglErr ExpandCenteredVarmaj(const uintptr_t* genovec, const uintptr_t* dosage_present, const Dosage* dosage_main, uint32_t variance_standardize, uint32_t is_haploid, uint32_t sample_ct, uint32_t dosage_ct, double maj_freq, double* normed_dosages) {
  const double nonmaj_freq = 1.0 - maj_freq;
  double inv_stdev;
  if (variance_standardize) {
    const double variance = 2 * maj_freq * nonmaj_freq;
    if (variance < kSmallEpsilon) {
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      GenovecCountFreqsUnsafe(genovec, sample_ct, genocounts);
      // remove unlikely() if any caller ever handles this case gracefully
      if (unlikely(dosage_ct || (genocounts[0] && (genocounts[1] || genocounts[2])) || (genocounts[1] && genocounts[2]))) {
        return kPglRetDegenerateData;
      }
      ZeroDArr(sample_ct, normed_dosages);
      return kPglRetSuccess;
    }
    inv_stdev = 1.0 / sqrt(variance);
    if (is_haploid) {
      // Variance is doubled in haploid case.
      inv_stdev *= (1.0 / kSqrt2);
    }
    // possible todo:
    // * Could use one inv_stdev for males and one for nonmales for chrX
    //   --score (while still leaving that out of GRM... or just leave males
    //   out there?).  This depends on dosage compensation model; discussed in
    //   e.g. GCTA paper.
  } else {
    inv_stdev = 1.0;
  }
  PopulateRescaledDosage(genovec, dosage_present, dosage_main, inv_stdev, -2 * nonmaj_freq * inv_stdev, 0.0, sample_ct, dosage_ct, normed_dosages);
  return kPglRetSuccess;
}

PglErr LoadCenteredVarmaj(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, uint32_t variance_standardize, uint32_t is_haploid, uint32_t sample_ct, uint32_t variant_uidx, AlleleCode maj_allele_idx, double maj_freq, PgenReader* simple_pgrp, uint32_t* missing_presentp, double* normed_dosages, uintptr_t* genovec_buf, uintptr_t* dosage_present_buf, Dosage* dosage_main_buf) {
  uint32_t dosage_ct;
  PglErr reterr = PgrGetInv1D(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, maj_allele_idx, simple_pgrp, genovec_buf, dosage_present_buf, dosage_main_buf, &dosage_ct);
  if (unlikely(reterr)) {
    // don't print malformed-.pgen error message here, since this is called
    // from multithreaded loops
    return reterr;
  }
  ZeroTrailingQuaters(sample_ct, genovec_buf);
  if (missing_presentp) {
    // missing_present assumed to be initialized to 0
    // this should probably be a library function...
    const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
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
  return ExpandCenteredVarmaj(genovec_buf, dosage_present_buf, dosage_main_buf, variance_standardize, is_haploid, sample_ct, dosage_ct, maj_freq, normed_dosages);
}

// multithread globals
static double* g_normed_dosage_vmaj_bufs[2] = {nullptr, nullptr};
static double* g_normed_dosage_smaj_bufs[2] = {nullptr, nullptr};

static double* g_grm = nullptr;

static uint32_t g_pca_sample_ct = 0;
static uint32_t g_cur_batch_size = 0;

CONSTI32(kGrmVariantBlockSize, 144);

// turns out dsyrk_ does exactly what we want here
THREAD_FUNC_DECL CalcGrmThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  assert(!tidx);
  const uint32_t sample_ct = g_pca_sample_ct;
  double* grm = g_grm;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_batch = g_is_last_thread_block;
    const uint32_t cur_batch_size = g_cur_batch_size;
    if (cur_batch_size) {
      TransposeMultiplySelfIncr(g_normed_dosage_vmaj_bufs[parity], sample_ct, cur_batch_size, grm);
    }
    if (is_last_batch) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

// can't use dsyrk_, so we manually partition the GRM piece we need to compute
// into an appropriate number of sub-pieces
THREAD_FUNC_DECL CalcGrmPartThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uintptr_t sample_ct = g_pca_sample_ct;
  const uintptr_t first_thread_row_start_idx = g_thread_start[0];
  const uintptr_t row_start_idx = g_thread_start[tidx];
  const uintptr_t row_ct = g_thread_start[tidx + 1] - row_start_idx;
  double* grm_piece = &(g_grm[(row_start_idx - first_thread_row_start_idx) * sample_ct]);
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_batch = g_is_last_thread_block;
    const uintptr_t cur_batch_size = g_cur_batch_size;
    if (cur_batch_size) {
      double* normed_vmaj = g_normed_dosage_vmaj_bufs[parity];
      double* normed_smaj = g_normed_dosage_smaj_bufs[parity];
      RowMajorMatrixMultiplyIncr(&(normed_smaj[row_start_idx * cur_batch_size]), normed_vmaj, row_ct, sample_ct, cur_batch_size, grm_piece);
    }
    if (is_last_batch) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

// missing_nz bit is set iff that sample has at least one missing entry in
// current block
static uintptr_t* g_missing_nz[2] = {nullptr, nullptr};
static uintptr_t* g_missing_smaj[2] = {nullptr, nullptr};
static uint32_t* g_missing_dbl_exclude_cts = nullptr;

CONSTI32(kDblMissingBlockWordCt, 2);
CONSTI32(kDblMissingBlockSize, kDblMissingBlockWordCt * kBitsPerWord);

THREAD_FUNC_DECL CalcDblMissingThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uint64_t first_thread_row_start_idx = g_thread_start[0];
  const uint64_t dbl_exclude_offset = (first_thread_row_start_idx * (first_thread_row_start_idx - 1)) / 2;
  const uint32_t row_start_idx = g_thread_start[tidx];
  const uintptr_t row_end_idx = g_thread_start[tidx + 1];
  uint32_t* missing_dbl_exclude_cts = g_missing_dbl_exclude_cts;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_batch = g_is_last_thread_block;

    // currently only care about zero vs. nonzero (I/O error)
    const uint32_t cur_batch_size = g_cur_batch_size;
    if (cur_batch_size) {
      const uintptr_t* missing_nz = g_missing_nz[parity];
      const uintptr_t* missing_smaj = g_missing_smaj[parity];
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
    }
    if (is_last_batch) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

PglErr CalcMissingMatrix(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* variant_include, uint32_t sample_ct, uint32_t variant_ct, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t row_start_idx, uintptr_t row_end_idx, uint32_t max_thread_ct, PgenReader* simple_pgrp, uint32_t** missing_cts_ptr, uint32_t** missing_dbl_exclude_cts_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  ThreadsState ts;
  InitThreads3z(&ts);
  PglErr reterr = kPglRetSuccess;
  {
    const uintptr_t row_end_idxl = BitCtToWordCt(row_end_idx);
    // bugfix (1 Oct 2017): missing_vmaj rows must be vector-aligned
    const uintptr_t row_end_idxaw = BitCtToAlignedWordCt(row_end_idx);
    uintptr_t* missing_vmaj = nullptr;
    uintptr_t* genovec_buf = nullptr;
    if (bigstack_calloc_u32(row_end_idx, missing_cts_ptr) ||
        bigstack_calloc_u32((S_CAST(uint64_t, row_end_idx) * (row_end_idx - 1) - S_CAST(uint64_t, row_start_idx) * (row_start_idx - 1)) / 2, missing_dbl_exclude_cts_ptr) ||
        bigstack_calloc_w(row_end_idxl, &g_missing_nz[0]) ||
        bigstack_calloc_w(row_end_idxl, &g_missing_nz[1]) ||
        bigstack_alloc_w(QuaterCtToWordCt(row_end_idx), &genovec_buf) ||
        bigstack_alloc_w(row_end_idxaw * (k1LU * kDblMissingBlockSize), &missing_vmaj) ||
        bigstack_alloc_w(RoundUpPow2(row_end_idx, 2) * kDblMissingBlockWordCt, &g_missing_smaj[0]) ||
        bigstack_alloc_w(RoundUpPow2(row_end_idx, 2) * kDblMissingBlockWordCt, &g_missing_smaj[1])) {
      goto CalcMissingMatrix_ret_NOMEM;
    }
    uint32_t* missing_cts = *missing_cts_ptr;
    uint32_t* missing_dbl_exclude_cts = *missing_dbl_exclude_cts_ptr;
    g_missing_dbl_exclude_cts = missing_dbl_exclude_cts;
    VecW* transpose_bitblock_wkspace = S_CAST(VecW*, bigstack_alloc_raw(kPglBitTransposeBufbytes));
    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    ts.calc_thread_ct = calc_thread_ct;
    if (bigstack_alloc_u32(calc_thread_ct + 1, &g_thread_start) ||
        bigstack_alloc_thread(calc_thread_ct, &ts.threads)) {
      goto CalcMissingMatrix_ret_NOMEM;
    }
    // note that this g_thread_start[] may have different values than the one
    // computed by CalcGrm(), since calc_thread_ct changes in the MTBLAS and
    // OS X cases.
    TriangleFill(sample_ct, calc_thread_ct, parallel_idx, parallel_tot, 0, 1, g_thread_start);
    assert(g_thread_start[0] == row_start_idx);
    assert(g_thread_start[calc_thread_ct] == row_end_idx);
    const uint32_t sample_transpose_batch_ct_m1 = (row_end_idx - 1) / kPglBitTransposeBatch;

    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    // caller's responsibility to print this
    // logputs("Correcting for missingness: ");
    fputs("0%", stdout);
    fflush(stdout);
    PgrClearLdCache(simple_pgrp);
    for (uint32_t cur_variant_idx_start = 0; ; ) {
      uint32_t cur_batch_size = 0;
      if (!ts.is_last_block) {
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
          reterr = PgrGetMissingnessD(sample_include, sample_include_cumulative_popcounts, row_end_idx, variant_uidx, simple_pgrp, nullptr, missing_vmaj_iter, nullptr, genovec_buf);
          if (unlikely(reterr)) {
            goto CalcMissingMatrix_ret_PGR_FAIL;
          }
          missing_vmaj_iter = &(missing_vmaj_iter[row_end_idxaw]);
        }
        uintptr_t* cur_missing_smaj_iter = g_missing_smaj[parity];
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
        uintptr_t* cur_missing_nz = g_missing_nz[parity];
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
        JoinThreads3z(&ts);
        // CalcDblMissingThread() never errors out
        if (ts.is_last_block) {
          break;
        }
        if (cur_variant_idx_start >= next_print_variant_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (cur_variant_idx_start * 100LLU) / variant_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
        }
      }
      ts.is_last_block = (cur_variant_idx_start + cur_batch_size == variant_ct);
      g_cur_batch_size = cur_batch_size;
      ts.thread_func_ptr = CalcDblMissingThread;
      if (unlikely(SpawnThreads3z(cur_variant_idx_start, &ts))) {
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
    bigstack_mark = R_CAST(unsigned char*, g_missing_nz[0]);
  }
  while (0) {
  CalcMissingMatrix_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  CalcMissingMatrix_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  CalcMissingMatrix_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads3z(&ts, &g_cur_batch_size);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr CalcGrm(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const AlleleCode* maj_alleles, const double* allele_freqs, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, GrmFlags grm_flags, uint32_t parallel_idx, uint32_t parallel_tot, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end, double** grm_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  CompressStreamState css;
  ThreadsState ts;
  PglErr reterr = kPglRetSuccess;
  PreinitCstream(&css);
  InitThreads3z(&ts);
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
    ts.calc_thread_ct = calc_thread_ct;
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
      TriangleFill(sample_ct, calc_thread_ct, parallel_idx, parallel_tot, 0, 1, thread_start);
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
    }
    g_thread_start = thread_start;
    double* grm;
    if (unlikely(bigstack_calloc_d((row_end_idx - row_start_idx) * row_end_idx, &grm))) {
      goto CalcGrm_ret_NOMEM;
    }
    g_pca_sample_ct = row_end_idx;
    g_grm = grm;
    const uint32_t row_end_idxl2 = QuaterCtToWordCt(row_end_idx);
    const uint32_t row_end_idxl = BitCtToWordCt(row_end_idx);
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* genovec_buf;
    uintptr_t* dosage_present_buf;
    Dosage* dosage_main_buf;
    if (unlikely(
            bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
            bigstack_alloc_thread(calc_thread_ct, &ts.threads) ||
            bigstack_alloc_w(row_end_idxl2, &genovec_buf) ||
            bigstack_alloc_w(row_end_idxl, &dosage_present_buf) ||
            bigstack_alloc_dosage(row_end_idx, &dosage_main_buf))) {
      goto CalcGrm_ret_NOMEM;
    }
    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    reterr = ConditionalAllocateNonAutosomalVariants(cip, "GRM construction", raw_variant_ct, &variant_include, &variant_ct);
    if (unlikely(reterr)) {
      goto CalcGrm_ret_1;
    }
    if (unlikely(
            bigstack_alloc_d(row_end_idx * kGrmVariantBlockSize, &g_normed_dosage_vmaj_bufs[0]) ||
            bigstack_alloc_d(row_end_idx * kGrmVariantBlockSize, &g_normed_dosage_vmaj_bufs[1]))) {
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
      if (unlikely(
              bigstack_alloc_d(row_end_idx * kGrmVariantBlockSize, &g_normed_dosage_smaj_bufs[0]) ||
              bigstack_alloc_d(row_end_idx * kGrmVariantBlockSize, &g_normed_dosage_smaj_bufs[1]))) {
        goto CalcGrm_ret_NOMEM;
      }
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
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t parity = 0;
    uint32_t cur_allele_ct = 2;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    logputs("Constructing GRM: ");
    fputs("0%", stdout);
    fflush(stdout);
    PgrClearLdCache(simple_pgrp);
    for (uint32_t cur_variant_idx_start = 0; ; ) {
      uint32_t cur_batch_size = 0;
      if (!ts.is_last_block) {
        cur_batch_size = kGrmVariantBlockSize;
        uint32_t cur_variant_idx_end = cur_variant_idx_start + cur_batch_size;
        if (cur_variant_idx_end > variant_ct) {
          cur_batch_size = variant_ct - cur_variant_idx_start;
          cur_variant_idx_end = variant_ct;
        }
        double* normed_vmaj_iter = g_normed_dosage_vmaj_bufs[parity];
        for (uint32_t variant_idx = cur_variant_idx_start; variant_idx != cur_variant_idx_end; ++variant_idx) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          const uint32_t maj_allele_idx = maj_alleles[variant_uidx];
          uint32_t missing_present = 0;
          uintptr_t allele_idx_base;
          if (!allele_idx_offsets) {
            allele_idx_base = variant_uidx;
          } else {
            allele_idx_base = allele_idx_offsets[variant_uidx];
            cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_base;
            allele_idx_base -= variant_uidx;
          }
          reterr = LoadCenteredVarmaj(sample_include, sample_include_cumulative_popcounts, variance_standardize, is_haploid, row_end_idx, variant_uidx, maj_allele_idx, GetAlleleFreq(&(allele_freqs[allele_idx_base]), maj_allele_idx, cur_allele_ct), simple_pgrp, variant_include_has_missing? (&missing_present) : nullptr, normed_vmaj_iter, genovec_buf, dosage_present_buf, dosage_main_buf);
          if (unlikely(reterr)) {
            if (reterr == kPglRetDegenerateData) {
              logputs("\n");
              logerrputs("Error: Zero-MAF variant is not actually monomorphic.  (This is possible when\ne.g. MAF is estimated from founders, but the minor allele was only observed in\nnonfounders.  In any case, you should be using e.g. --maf to filter out all\nvery-low-MAF variants, since the relationship matrix distance formula does not\nhandle them well.)\n");
            }
            goto CalcGrm_ret_PGR_FAIL;
          }
          if (missing_present) {
            SetBit(variant_uidx, variant_include_has_missing);
          }
          normed_vmaj_iter = &(normed_vmaj_iter[row_end_idx]);
        }
        if (thread_start) {
          MatrixTransposeCopy(g_normed_dosage_vmaj_bufs[parity], cur_batch_size, row_end_idx, g_normed_dosage_smaj_bufs[parity]);
        }
      }
      if (cur_variant_idx_start) {
        JoinThreads3z(&ts);
        // CalcGrmPartThread() and CalcGrmThread() never error out
        if (ts.is_last_block) {
          break;
        }
        if (cur_variant_idx_start >= next_print_variant_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (cur_variant_idx_start * 100LLU) / variant_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
        }
      }
      ts.is_last_block = (cur_variant_idx_start + cur_batch_size == variant_ct);
      g_cur_batch_size = cur_batch_size;
      if (!ts.thread_func_ptr) {
        if (thread_start) {
          ts.thread_func_ptr = CalcGrmPartThread;
        } else {
          ts.thread_func_ptr = CalcGrmThread;
        }
      }
      if (unlikely(SpawnThreads3z(cur_variant_idx_start, &ts))) {
        goto CalcGrm_ret_THREAD_CREATE_FAIL;
      }
      cur_variant_idx_start += cur_batch_size;
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
    // --make-grm-list/--make-grm-bin
    // (note that this routine may also be called by --pca, which may not write
    // a matrix to disk at all.)
    if (grm_flags & (kfGrmMatrixShapemask | kfGrmListmask | kfGrmBin)) {
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
              // (roughly same performance as creating a zero-tab constant
              // buffer in advance)
              const uint32_t zcount = sample_ct - row_idx;
              const uint32_t wct = DivUp(zcount, kBytesPerWord / 2);
              // assumes little-endian
              const uintptr_t zerotab_word = 0x930 * kMask0001;
#ifdef __arm__
#  error "Unaligned accesses in CalcGrm()."
#endif
              uintptr_t* writep_alias = R_CAST(uintptr_t*, cswritep);
              for (uintptr_t widx = 0; widx != wct; ++widx) {
                *writep_alias++ = zerotab_word;
              }
              cswritep = &(cswritep[zcount * 2]);
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
            const uintptr_t tot_cells = (S_CAST(uint64_t, row_end_idx) * (row_end_idx - 1) - S_CAST(uint64_t, row_start_idx) * (row_start_idx - 1)) / 2;
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
          // --make-grm-list
          char* outname_end2 = strcpya_k(outname_end, ".grm");
          if (parallel_tot != 1) {
            *outname_end2++ = '.';
            outname_end2 = u32toa(parallel_idx + 1, outname_end2);
          }
          if (grm_flags & kfGrmListZs) {
            outname_end2 = strcpya_k(outname_end2, ".zst");
          }
          *outname_end2 = '\0';
          reterr = InitCstreamAlloc(outname, 0, !(grm_flags & kfGrmListNoGz), max_thread_ct, kCompressStreamBlock + kMaxMediumLine, &css, &cswritep);
          if (unlikely(reterr)) {
            goto CalcGrm_ret_1;
          }
          fputs("--make-grm-list: Writing...", stdout);
          fflush(stdout);
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
          if (unlikely(CswriteCloseNull(&css, cswritep))) {
            goto CalcGrm_ret_WRITE_FAIL;
          }
          putc_unlocked('\r', stdout);
          log_write_iter = strcpya_k(g_logbuf, "--make-grm-list: GRM ");
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
  CalcGrm_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  CalcGrm_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
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
  CleanupThreads3z(&ts, &g_cur_batch_size);
  BLAS_SET_NUM_THREADS(1);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

// should be able to remove NOLAPACK later since we already have a non-LAPACK
// SVD implementation
#ifndef NOLAPACK
// this seems to be better than 256 (due to avoidance of cache critical
// stride?)
// (still want this to be a multiple of 8, for cleaner multithreading)
CONSTI32(kPcaVariantBlockSize, 240);

// multithread globals
static uintptr_t* g_genovecs[2] = {nullptr, nullptr};
static uint32_t* g_dosage_cts[2] = {nullptr, nullptr};
static uintptr_t* g_dosage_presents[2] = {nullptr, nullptr};
static Dosage* g_dosage_mains[2] = {nullptr, nullptr};
static double* g_cur_maj_freqs[2] = {nullptr, nullptr};
static double** g_yy_bufs = nullptr;
static double** g_y_transpose_bufs = nullptr;
static double** g_g2_bb_part_bufs = nullptr;
static double* g_g1 = nullptr;
static double* g_qq = nullptr;

static uint32_t g_pc_ct = 0;
static uint32_t g_is_haploid = 0;
static PglErr g_error_ret = kPglRetSuccess;

THREAD_FUNC_DECL CalcPcaXtxaThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uint32_t pca_sample_ct = g_pca_sample_ct;
  const uintptr_t pca_sample_ctaw2 = QuaterCtToAlignedWordCt(pca_sample_ct);
  const uintptr_t pca_sample_ctaw = BitCtToAlignedWordCt(pca_sample_ct);
  const uint32_t pc_ct_x2 = g_pc_ct * 2;
  const uintptr_t qq_col_ct = (g_pc_ct + 1) * pc_ct_x2;
  const uint32_t vidx_offset = tidx * kPcaVariantBlockSize;
  const double* g1 = g_g1;
  double* qq_iter = g_qq;
  double* yy_buf = g_yy_bufs[tidx];
  double* y_transpose_buf = g_y_transpose_bufs[tidx];
  double* g2_part_buf = g_g2_bb_part_bufs[tidx];
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_batch = g_is_last_thread_block;
    const uint32_t cur_batch_size = g_cur_batch_size;
    if (vidx_offset < cur_batch_size) {
      uint32_t cur_thread_batch_size = cur_batch_size - vidx_offset;
      if (cur_thread_batch_size > kPcaVariantBlockSize) {
        cur_thread_batch_size = kPcaVariantBlockSize;
      }
      const uintptr_t* genovec_iter = &(g_genovecs[parity][vidx_offset * pca_sample_ctaw2]);
      const uint32_t* cur_dosage_cts = &(g_dosage_cts[parity][vidx_offset]);
      const uintptr_t* dosage_present_iter = &(g_dosage_presents[parity][vidx_offset * pca_sample_ctaw]);
      const Dosage* dosage_main_iter = &(g_dosage_mains[parity][vidx_offset * pca_sample_ct]);
      const double* cur_maj_freqs_iter = &(g_cur_maj_freqs[parity][vidx_offset]);
      double* yy_iter = yy_buf;
      for (uint32_t uii = 0; uii != cur_thread_batch_size; ++uii) {
        // instead of setting is_haploid here, we just divide eigenvalues by 2
        // at the end if necessary
        PglErr reterr = ExpandCenteredVarmaj(genovec_iter, dosage_present_iter, dosage_main_iter, 1, 0, pca_sample_ct, cur_dosage_cts[uii], cur_maj_freqs_iter[uii], yy_iter);
        if (unlikely(reterr)) {
          g_error_ret = reterr;
          break;
        }
        yy_iter = &(yy_iter[pca_sample_ct]);
        genovec_iter = &(genovec_iter[pca_sample_ctaw2]);
        dosage_present_iter = &(dosage_present_iter[pca_sample_ctaw]);
        dosage_main_iter = &(dosage_main_iter[pca_sample_ct]);
      }
      double* cur_qq = &(qq_iter[vidx_offset * qq_col_ct]);
      RowMajorMatrixMultiplyStrided(yy_buf, g1, cur_thread_batch_size, pca_sample_ct, pc_ct_x2, pc_ct_x2, pca_sample_ct, qq_col_ct, cur_qq);
      MatrixTransposeCopy(yy_buf, cur_thread_batch_size, pca_sample_ct, y_transpose_buf);
      RowMajorMatrixMultiplyStridedIncr(y_transpose_buf, cur_qq, pca_sample_ct, cur_thread_batch_size, pc_ct_x2, qq_col_ct, cur_thread_batch_size, pc_ct_x2, g2_part_buf);
      qq_iter = &(qq_iter[cur_batch_size * qq_col_ct]);
    }
    if (is_last_batch) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

THREAD_FUNC_DECL CalcPcaXaThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uint32_t pca_sample_ct = g_pca_sample_ct;
  const uintptr_t pca_sample_ctaw2 = QuaterCtToAlignedWordCt(pca_sample_ct);
  const uintptr_t pca_sample_ctaw = BitCtToAlignedWordCt(pca_sample_ct);
  const uint32_t pc_ct_x2 = g_pc_ct * 2;
  const uintptr_t qq_col_ct = (g_pc_ct + 1) * pc_ct_x2;
  const uint32_t vidx_offset = tidx * kPcaVariantBlockSize;
  const double* g1 = g_g1;
  double* qq_iter = g_qq;
  double* yy_buf = g_yy_bufs[tidx];
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_batch = g_is_last_thread_block;
    const uint32_t cur_batch_size = g_cur_batch_size;
    if (vidx_offset < cur_batch_size) {
      uint32_t cur_thread_batch_size = cur_batch_size - vidx_offset;
      if (cur_thread_batch_size > kPcaVariantBlockSize) {
        cur_thread_batch_size = kPcaVariantBlockSize;
      }
      const uintptr_t* genovec_iter = &(g_genovecs[parity][vidx_offset * pca_sample_ctaw2]);
      const uint32_t* cur_dosage_cts = &(g_dosage_cts[parity][vidx_offset]);
      const uintptr_t* dosage_present_iter = &(g_dosage_presents[parity][vidx_offset * pca_sample_ctaw]);
      const Dosage* dosage_main_iter = &(g_dosage_mains[parity][vidx_offset * pca_sample_ct]);
      const double* cur_maj_freqs_iter = &(g_cur_maj_freqs[parity][vidx_offset]);
      double* yy_iter = yy_buf;
      for (uint32_t uii = 0; uii != cur_thread_batch_size; ++uii) {
        PglErr reterr = ExpandCenteredVarmaj(genovec_iter, dosage_present_iter, dosage_main_iter, 1, 0, pca_sample_ct, cur_dosage_cts[uii], cur_maj_freqs_iter[uii], yy_iter);
        if (unlikely(reterr)) {
          g_error_ret = reterr;
          break;
        }
        yy_iter = &(yy_iter[pca_sample_ct]);
        genovec_iter = &(genovec_iter[pca_sample_ctaw2]);
        dosage_present_iter = &(dosage_present_iter[pca_sample_ctaw]);
        dosage_main_iter = &(dosage_main_iter[pca_sample_ct]);
      }
      double* cur_qq = &(qq_iter[vidx_offset * qq_col_ct]);
      RowMajorMatrixMultiplyStrided(yy_buf, g1, cur_thread_batch_size, pca_sample_ct, pc_ct_x2, pc_ct_x2, pca_sample_ct, qq_col_ct, cur_qq);
      qq_iter = &(qq_iter[cur_batch_size * qq_col_ct]);
    }
    if (is_last_batch) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

THREAD_FUNC_DECL CalcPcaXtbThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uint32_t pca_sample_ct = g_pca_sample_ct;
  const uintptr_t pca_sample_ctaw2 = QuaterCtToAlignedWordCt(pca_sample_ct);
  const uintptr_t pca_sample_ctaw = BitCtToAlignedWordCt(pca_sample_ct);
  const uint32_t pc_ct_x2 = g_pc_ct * 2;
  const uintptr_t qq_col_ct = (g_pc_ct + 1) * pc_ct_x2;
  const uint32_t vidx_offset = tidx * kPcaVariantBlockSize;
  const double* qq_iter = &(g_qq[vidx_offset * qq_col_ct]);
  double* yy_buf = g_yy_bufs[tidx];
  double* y_transpose_buf = g_y_transpose_bufs[tidx];
  double* bb_part_buf = g_g2_bb_part_bufs[tidx];
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_batch = g_is_last_thread_block;
    const uint32_t cur_batch_size = g_cur_batch_size;
    if (vidx_offset < cur_batch_size) {
      uint32_t cur_thread_batch_size = cur_batch_size - vidx_offset;
      if (cur_thread_batch_size > kPcaVariantBlockSize) {
        cur_thread_batch_size = kPcaVariantBlockSize;
      }
      const uintptr_t* genovec_iter = &(g_genovecs[parity][vidx_offset * pca_sample_ctaw2]);
      const uint32_t* cur_dosage_cts = &(g_dosage_cts[parity][vidx_offset]);
      const uintptr_t* dosage_present_iter = &(g_dosage_presents[parity][vidx_offset * pca_sample_ctaw]);
      const Dosage* dosage_main_iter = &(g_dosage_mains[parity][vidx_offset * pca_sample_ct]);
      const double* cur_maj_freqs_iter = &(g_cur_maj_freqs[parity][vidx_offset]);
      double* yy_iter = yy_buf;
      for (uint32_t uii = 0; uii != cur_thread_batch_size; ++uii) {
        PglErr reterr = ExpandCenteredVarmaj(genovec_iter, dosage_present_iter, dosage_main_iter, 1, 0, pca_sample_ct, cur_dosage_cts[uii], cur_maj_freqs_iter[uii], yy_iter);
        if (unlikely(reterr)) {
          g_error_ret = reterr;
          break;
        }
        yy_iter = &(yy_iter[pca_sample_ct]);
        genovec_iter = &(genovec_iter[pca_sample_ctaw2]);
        dosage_present_iter = &(dosage_present_iter[pca_sample_ctaw]);
        dosage_main_iter = &(dosage_main_iter[pca_sample_ct]);
      }
      MatrixTransposeCopy(yy_buf, cur_thread_batch_size, pca_sample_ct, y_transpose_buf);
      RowMajorMatrixMultiplyIncr(y_transpose_buf, qq_iter, pca_sample_ct, qq_col_ct, cur_thread_batch_size, bb_part_buf);
      qq_iter = &(qq_iter[cur_batch_size * qq_col_ct]);
    }
    if (is_last_batch) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

THREAD_FUNC_DECL CalcPcaVarWtsThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  const uint32_t pca_sample_ct = g_pca_sample_ct;
  const uintptr_t pca_sample_ctaw2 = QuaterCtToAlignedWordCt(pca_sample_ct);
  const uintptr_t pca_sample_ctaw = BitCtToAlignedWordCt(pca_sample_ct);
  const uint32_t pc_ct = g_pc_ct;
  const uint32_t is_haploid = g_is_haploid;
  const uint32_t vidx_offset = tidx * kPcaVariantBlockSize;

  // either first batch size is calc_thread_ct * kPcaVariantBlockSize, or there
  // is only one batch
  const uintptr_t var_wts_part_size = S_CAST(uintptr_t, pc_ct) * g_cur_batch_size;

  const double* sample_wts = g_g1;  // sample-major, pc_ct columns
  double* yy_buf = g_yy_bufs[tidx];
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_batch = g_is_last_thread_block;
    const uint32_t cur_batch_size = g_cur_batch_size;
    if (vidx_offset < cur_batch_size) {
      uint32_t cur_thread_batch_size = cur_batch_size - vidx_offset;
      if (cur_thread_batch_size > kPcaVariantBlockSize) {
        cur_thread_batch_size = kPcaVariantBlockSize;
      }
      const uintptr_t* genovec_iter = &(g_genovecs[parity][vidx_offset * pca_sample_ctaw2]);
      const uint32_t* cur_dosage_cts = &(g_dosage_cts[parity][vidx_offset]);
      const uintptr_t* dosage_present_iter = &(g_dosage_presents[parity][vidx_offset * pca_sample_ctaw]);
      const Dosage* dosage_main_iter = &(g_dosage_mains[parity][vidx_offset * pca_sample_ct]);
      const double* cur_maj_freqs_iter = &(g_cur_maj_freqs[parity][vidx_offset]);
      double* yy_iter = yy_buf;
      for (uint32_t uii = 0; uii != cur_thread_batch_size; ++uii) {
        PglErr reterr = ExpandCenteredVarmaj(genovec_iter, dosage_present_iter, dosage_main_iter, 1, is_haploid, pca_sample_ct, cur_dosage_cts[uii], cur_maj_freqs_iter[uii], yy_iter);
        if (unlikely(reterr)) {
          g_error_ret = reterr;
          break;
        }
        yy_iter = &(yy_iter[pca_sample_ct]);
        genovec_iter = &(genovec_iter[pca_sample_ctaw2]);
        dosage_present_iter = &(dosage_present_iter[pca_sample_ctaw]);
        dosage_main_iter = &(dosage_main_iter[pca_sample_ct]);
      }
      // Variant weight matrix = X^T * S * D^{-1/2}, where X^T is the
      // variance-standardized genotype matrix, S is the sample weight matrix,
      // and D is a diagonal eigenvalue matrix.
      // We postpone the D^{-1/2} part for now, but it's straightforward to
      // switch to using precomputed (S * D^{-1/2}).
      double* cur_var_wts_part = &(g_qq[parity * var_wts_part_size + vidx_offset * S_CAST(uintptr_t, pc_ct)]);
      RowMajorMatrixMultiply(yy_buf, sample_wts, cur_thread_batch_size, pc_ct, pca_sample_ct, cur_var_wts_part);
    }
    if (is_last_batch) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

PglErr CalcPca(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* maj_alleles, const double* allele_freqs, uint32_t raw_sample_ct, uintptr_t pca_sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t pc_ct, PcaFlags pca_flags, uint32_t max_thread_ct, PgenReader* simple_pgrp, sfmt_t* sfmtp, double* grm, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  char* cswritep = nullptr;
  CompressStreamState css;
  ThreadsState ts;
  InitThreads3z(&ts);
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
    uint32_t calc_thread_ct = (max_thread_ct > 8)? (max_thread_ct - 1) : max_thread_ct;
    if ((calc_thread_ct - 1) * kPcaVariantBlockSize >= variant_ct) {
      calc_thread_ct = 1 + (variant_ct - 1) / kPcaVariantBlockSize;
    }
#endif
    ts.calc_thread_ct = calc_thread_ct;
    if (pc_ct > pca_sample_ct) {
      if (pca_sample_ct <= variant_ct) {
        pc_ct = pca_sample_ct;
        snprintf(g_logbuf, kLogbufSize, "Warning: calculating %u PCs, since there are only %u samples.\n", pc_ct, pc_ct);
      } else {
        pc_ct = variant_ct;
        snprintf(g_logbuf, kLogbufSize, "Warning: calculating %u PCs, since there are only %u autosomal variants.\n", pc_ct, pc_ct);
      }
      if (unlikely(pc_ct < 2)) {
        logerrputs("Error: Too few samples or autosomal variants for PCA.\n");
        goto CalcPca_ret_DEGENERATE_DATA;
      }
      logerrputsb();
    }
    const uint32_t var_wts = (pca_flags / kfPcaBiallelicVarWts) & 1;
    const uint32_t require_biallelic = var_wts && (!(pca_flags & kfPcaIgnoreBiallelicVarWtsRestriction));
    const uint32_t chr_col = pca_flags & kfPcaVcolChrom;
    const uint32_t ref_col = pca_flags & kfPcaVcolRef;
    const uint32_t alt1_col = pca_flags & kfPcaVcolAlt1;
    const uint32_t alt_col = pca_flags & kfPcaVcolAlt;
    const uint32_t maj_col = pca_flags & kfPcaVcolMaj;
    const uint32_t nonmaj_col = pca_flags & kfPcaVcolNonmaj;
    double* cur_var_wts = nullptr;
    double* eigval_inv_sqrts = nullptr;
    char* chr_buf = nullptr;
    uintptr_t overflow_buf_size = 3 * kMaxMediumLine;
    if (var_wts) {
      if (unlikely(
              bigstack_alloc_d(pc_ct, &cur_var_wts) ||
              bigstack_alloc_d(pc_ct, &eigval_inv_sqrts))) {
        goto CalcPca_ret_NOMEM;
      }
      uint32_t max_chr_blen = 0;
      if (chr_col) {
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
    const uint32_t pca_sample_ctaw2 = QuaterCtToAlignedWordCt(pca_sample_ct);
    const uint32_t pca_sample_ctaw = BitCtToAlignedWordCt(pca_sample_ct);
    uint32_t* pca_sample_include_cumulative_popcounts;
    double* eigvals;
    if (unlikely(
            bigstack_alloc_u32(raw_sample_ctl, &pca_sample_include_cumulative_popcounts) ||
            bigstack_alloc_d(pc_ct, &eigvals) ||
            bigstack_alloc_thread(calc_thread_ct, &ts.threads) ||
            bigstack_alloc_dp(calc_thread_ct, &g_yy_bufs))) {
      goto CalcPca_ret_NOMEM;
    }
    FillCumulativePopcounts(pca_sample_include, raw_sample_ctl, pca_sample_include_cumulative_popcounts);
    g_pca_sample_ct = pca_sample_ct;
    g_pc_ct = pc_ct;
    g_error_ret = kPglRetSuccess;
    g_is_haploid = cip->haploid_mask[0] & 1;
    uint32_t cur_allele_ct = 2;
    double* qq = nullptr;
    double* eigvecs_smaj;
    char* writebuf;
    if (is_approx) {
      if (pca_sample_ct <= 5000) {
        logerrputs("Warning: \"--pca approx\" is only recommended for analysis of >5000 samples.\n");
      }
      if (variant_ct > 5000000) {
        logerrputs("Warning: Use of \"--pca approx\" on >5m variants is not advisable.  Apply a MAF\nfilter if you haven't done so yet, and consider LD-pruning your variant set as\nwell.\n");
      }
      // This is ported from EIGENSOFT 6 src/ksrc/kjg_fpca.c , which is in turn
      // primarily based on Halko N, Martinsson P, Shkolnisky Y, Tygert M
      // (2011) An Algorithm for the Principal Component Analysis of Large Data
      // Sets.
      const uintptr_t pc_ct_x2 = pc_ct * 2;
      const uintptr_t qq_col_ct = (pc_ct + 1) * pc_ct_x2;
      // bugfix (30 Jan 2019): First SvdRect() call returns min(variant_ct,
      // qq_col_ct) singular vectors; this was previously assumed to always be
      // qq_col_ct, and very inaccurate results were produced when the
      // assumption wasn't true.
      // Simplest solution is to force the user to request fewer PCs, since the
      // final PCs wouldn't be accurate anyway.
      if (qq_col_ct > variant_ct) {
        logerrprintfww("Error: Too few variants to compute %u PCs with \"--pca approx\" (%u required).\n", pc_ct, qq_col_ct);
        goto CalcPca_ret_DEGENERATE_DATA;
      }
#ifndef LAPACK_ILP64
      if (unlikely((variant_ct * S_CAST(uint64_t, qq_col_ct)) > 0x7effffff)) {
        logerrputs("Error: \"--pca approx\" problem instance too large for this " PROG_NAME_STR " build.  If\nthis is really the computation you want, use a " PROG_NAME_STR " build with large-matrix\nsupport.\n");
        goto CalcPca_ret_INCONSISTENT_INPUT;
      }
#endif
      const double variant_ct_recip = 1.0 / u31tod(variant_ct);

      const uintptr_t g_size = pca_sample_ct * pc_ct_x2;
      __CLPK_integer svd_rect_lwork;
#ifdef LAPACK_ILP64
      GetSvdRectLwork(MAXV(pca_sample_ct, variant_ct), qq_col_ct, &svd_rect_lwork);
#else
      if (unlikely(GetSvdRectLwork(MAXV(pca_sample_ct, variant_ct), qq_col_ct, &svd_rect_lwork))) {
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
      if (unlikely(
              bigstack_alloc_d(qq_col_ct, &ss) ||
              bigstack_alloc_d(variant_ct * qq_col_ct, &qq) ||
              bigstack_alloc_dp(calc_thread_ct, &g_y_transpose_bufs) ||
              bigstack_alloc_dp(calc_thread_ct, &g_g2_bb_part_bufs) ||
              bigstack_alloc_uc(svd_rect_wkspace_size, &svd_rect_wkspace) ||
              bigstack_alloc_d(g_size, &g1))) {
        goto CalcPca_ret_NOMEM;
      }
      const uintptr_t genovecs_alloc = RoundUpPow2(pca_sample_ctaw2 * kPcaVariantBlockSize * sizeof(intptr_t), kCacheline);
      const uintptr_t dosage_cts_alloc = RoundUpPow2(kPcaVariantBlockSize * sizeof(int32_t), kCacheline);
      const uintptr_t dosage_presents_alloc = RoundUpPow2(pca_sample_ctaw * kPcaVariantBlockSize * sizeof(intptr_t), kCacheline);
      const uintptr_t dosage_main_alloc = RoundUpPow2(pca_sample_ct * kPcaVariantBlockSize * sizeof(Dosage), kCacheline);
      const uintptr_t cur_maj_freqs_alloc = RoundUpPow2(kPcaVariantBlockSize * sizeof(double), kCacheline);
      const uintptr_t yy_alloc = RoundUpPow2(kPcaVariantBlockSize * pca_sample_ct * sizeof(double), kCacheline);
      const uintptr_t b_size = pca_sample_ct * qq_col_ct;
      const uintptr_t g2_bb_part_alloc = RoundUpPow2(b_size * sizeof(double), kCacheline);
      const uintptr_t per_thread_alloc = 2 * (genovecs_alloc + dosage_cts_alloc + dosage_presents_alloc + dosage_main_alloc + cur_maj_freqs_alloc + yy_alloc) + g2_bb_part_alloc;

      const uintptr_t bigstack_avail = bigstack_left();
      if (per_thread_alloc * calc_thread_ct > bigstack_avail) {
        if (unlikely(bigstack_avail < per_thread_alloc)) {
          goto CalcPca_ret_NOMEM;
        }
        calc_thread_ct = bigstack_avail / per_thread_alloc;
      }
      for (uint32_t parity = 0; parity != 2; ++parity) {
        g_genovecs[parity] = S_CAST(uintptr_t*, bigstack_alloc_raw(genovecs_alloc * calc_thread_ct));
        g_dosage_cts[parity] = S_CAST(uint32_t*, bigstack_alloc_raw(dosage_cts_alloc * calc_thread_ct));
        g_dosage_presents[parity] = S_CAST(uintptr_t*, bigstack_alloc_raw(dosage_presents_alloc * calc_thread_ct));
        g_dosage_mains[parity] = S_CAST(Dosage*, bigstack_alloc_raw(dosage_main_alloc * calc_thread_ct));
        g_cur_maj_freqs[parity] = S_CAST(double*, bigstack_alloc_raw(cur_maj_freqs_alloc * calc_thread_ct));
      }
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        g_yy_bufs[tidx] = S_CAST(double*, bigstack_alloc_raw(yy_alloc));
        g_y_transpose_bufs[tidx] = S_CAST(double*, bigstack_alloc_raw(yy_alloc));
        g_g2_bb_part_bufs[tidx] = S_CAST(double*, bigstack_alloc_raw(g2_bb_part_alloc));
      }
      FillGaussianDArr(g_size / 2, max_thread_ct, sfmtp, g1);
      g_g1 = g1;
#ifdef __APPLE__
      fputs("Projecting random vectors... ", stdout);
#else
      printf("Projecting random vectors (%u compute thread%s)... ", calc_thread_ct, (calc_thread_ct == 1)? "" : "s");
#endif
      fflush(stdout);
      PgrClearLdCache(simple_pgrp);
      for (uint32_t iter_idx = 0; iter_idx <= pc_ct; ++iter_idx) {
        // kjg_fpca_XTXA(), kjg_fpca_XA()
        for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
          ZeroDArr(g_size, g_g2_bb_part_bufs[tidx]);
        }
        double* qq_iter = &(qq[iter_idx * pc_ct_x2]);  // offset on first row
        g_qq = qq_iter;

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
        uintptr_t variant_uidx_base = 0;
        uintptr_t cur_bits = variant_include[0];
        uint32_t parity = 0;
        for (uint32_t cur_variant_idx_start = 0; ; ) {
          uint32_t cur_batch_size = 0;
          if (!ts.is_last_block) {
            cur_batch_size = calc_thread_ct * kPcaVariantBlockSize;
            uint32_t cur_variant_idx_end = cur_variant_idx_start + cur_batch_size;
            if (cur_variant_idx_end > variant_ct) {
              cur_batch_size = variant_ct - cur_variant_idx_start;
              cur_variant_idx_end = variant_ct;
            }
            uintptr_t* genovec_iter = g_genovecs[parity];
            uint32_t* dosage_ct_iter = g_dosage_cts[parity];
            uintptr_t* dosage_present_iter = g_dosage_presents[parity];
            Dosage* dosage_main_iter = g_dosage_mains[parity];
            double* maj_freqs_write_iter = g_cur_maj_freqs[parity];
            for (uint32_t variant_idx = cur_variant_idx_start; variant_idx != cur_variant_idx_end; ++variant_idx) {
              const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
              const uint32_t maj_allele_idx = maj_alleles[variant_uidx];
              uint32_t dosage_ct;
              reterr = PgrGetInv1D(pca_sample_include, pca_sample_include_cumulative_popcounts, pca_sample_ct, variant_uidx, maj_allele_idx, simple_pgrp, genovec_iter, dosage_present_iter, dosage_main_iter, &dosage_ct);
              if (unlikely(reterr)) {
                goto CalcPca_ret_PGR_FAIL;
              }
              ZeroTrailingQuaters(pca_sample_ct, genovec_iter);
              genovec_iter = &(genovec_iter[pca_sample_ctaw2]);
              *dosage_ct_iter++ = dosage_ct;
              dosage_present_iter = &(dosage_present_iter[pca_sample_ctaw]);
              dosage_main_iter = &(dosage_main_iter[pca_sample_ct]);
              uintptr_t allele_idx_base;
              if (!allele_idx_offsets) {
                allele_idx_base = variant_uidx;
              } else {
                allele_idx_base = allele_idx_offsets[variant_uidx];
                cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_base;
                if (require_biallelic && (cur_allele_ct != 2)) {
                  logputs("\n");
                  logerrputs("Error: Multiallelic variant present in \"--pca biallelic-var-wts\" run.\n");
                  goto CalcPca_ret_INCONSISTENT_INPUT;
                }
                allele_idx_base -= variant_uidx;
              }
              // bugfix (23 Jul 2017): we already subtracted variant_uidx
              *maj_freqs_write_iter++ = GetAlleleFreq(&(allele_freqs[allele_idx_base]), maj_allele_idx, cur_allele_ct);
            }
          }
          if (cur_variant_idx_start) {
            JoinThreads3z(&ts);
            reterr = g_error_ret;
            if (unlikely(reterr)) {
              logputs("\n");
              logerrputs("Error: Zero-MAF variant is not actually monomorphic.  (This is possible when\ne.g. MAF is estimated from founders, but the minor allele was only observed in\nnonfounders.  In any case, you should be using e.g. --maf to filter out all\nvery-low-MAF variants, since the relationship matrix distance formula does not\nhandle them well.)\n");
              goto CalcPca_ret_1;
            }
            if (ts.is_last_block) {
              break;
            }
          }
          if (!cur_variant_idx_start) {
            if (iter_idx < pc_ct) {
              ts.thread_func_ptr = CalcPcaXtxaThread;
            } else {
              ts.thread_func_ptr = CalcPcaXaThread;
            }
          }
          ts.is_last_block = (cur_variant_idx_start + cur_batch_size == variant_ct);
          g_cur_batch_size = cur_batch_size;
          if (unlikely(SpawnThreads3z(cur_variant_idx_start, &ts))) {
            goto CalcPca_ret_THREAD_CREATE_FAIL;
          }
          cur_variant_idx_start += cur_batch_size;
          parity = 1 - parity;
        }
        if (iter_idx < pc_ct) {
          memcpy(g1, g_g2_bb_part_bufs[0], g_size * sizeof(double));
          for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
            const double* cur_g2_part = g_g2_bb_part_bufs[tidx];
            for (uintptr_t ulii = 0; ulii != g_size; ++ulii) {
              g1[ulii] += cur_g2_part[ulii];
            }
          }
          for (uintptr_t ulii = 0; ulii != g_size; ++ulii) {
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
      IntErr svd_rect_err = SvdRect(variant_ct, qq_col_ct, svd_rect_lwork, qq, ss, svd_rect_wkspace);
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
        ZeroDArr(b_size, g_g2_bb_part_bufs[tidx]);
      }
      uintptr_t variant_uidx_base = 0;
      uintptr_t cur_bits = variant_include[0];
      uint32_t parity = 0;
      ReinitThreads3z(&ts);
      g_qq = qq;
      for (uint32_t cur_variant_idx_start = 0; ; ) {
        uint32_t cur_batch_size = 0;
        if (!ts.is_last_block) {
          // probable todo: move this boilerplate in its own function
          cur_batch_size = calc_thread_ct * kPcaVariantBlockSize;
          uint32_t cur_variant_idx_end = cur_variant_idx_start + cur_batch_size;
          if (cur_variant_idx_end > variant_ct) {
            cur_batch_size = variant_ct - cur_variant_idx_start;
            cur_variant_idx_end = variant_ct;
          }
          uintptr_t* genovec_iter = g_genovecs[parity];
          uint32_t* dosage_ct_iter = g_dosage_cts[parity];
          uintptr_t* dosage_present_iter = g_dosage_presents[parity];
          Dosage* dosage_main_iter = g_dosage_mains[parity];
          double* maj_freqs_write_iter = g_cur_maj_freqs[parity];
          for (uint32_t variant_idx = cur_variant_idx_start; variant_idx != cur_variant_idx_end; ++variant_idx) {
            const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
            const uint32_t maj_allele_idx = maj_alleles[variant_uidx];
            uint32_t dosage_ct;
            reterr = PgrGetInv1D(pca_sample_include, pca_sample_include_cumulative_popcounts, pca_sample_ct, variant_uidx, maj_allele_idx, simple_pgrp, genovec_iter, dosage_present_iter, dosage_main_iter, &dosage_ct);
            if (unlikely(reterr)) {
              goto CalcPca_ret_PGR_FAIL;
            }
            ZeroTrailingQuaters(pca_sample_ct, genovec_iter);
            genovec_iter = &(genovec_iter[pca_sample_ctaw2]);
            *dosage_ct_iter++ = dosage_ct;
            dosage_present_iter = &(dosage_present_iter[pca_sample_ctaw]);
            dosage_main_iter = &(dosage_main_iter[pca_sample_ct]);
            uintptr_t allele_idx_base;
            if (!allele_idx_offsets) {
              allele_idx_base = variant_uidx;
            } else {
              allele_idx_base = allele_idx_offsets[variant_uidx];
              cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_base;
              allele_idx_base -= variant_uidx;
            }
            *maj_freqs_write_iter++ = GetAlleleFreq(&(allele_freqs[allele_idx_base]), maj_allele_idx, cur_allele_ct);
          }
        }
        if (cur_variant_idx_start) {
          JoinThreads3z(&ts);
          if (unlikely(g_error_ret)) {
            // this error *didn't* happen on an earlier pass, so assign blame
            // to I/O instead
            goto CalcPca_ret_REWIND_FAIL;
          }
          if (ts.is_last_block) {
            break;
          }
        }
        ts.is_last_block = (cur_variant_idx_start + cur_batch_size == variant_ct);
        g_cur_batch_size = cur_batch_size;
        ts.thread_func_ptr = CalcPcaXtbThread;
        if (unlikely(SpawnThreads3z(cur_variant_idx_start, &ts))) {
          goto CalcPca_ret_THREAD_CREATE_FAIL;
        }
        cur_variant_idx_start += cur_batch_size;
        parity = 1 - parity;
      }
      double* bb = g_g2_bb_part_bufs[0];
      for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
        const double* cur_bb_part = g_g2_bb_part_bufs[tidx];
        for (uintptr_t ulii = 0; ulii != b_size; ++ulii) {
          bb[ulii] += cur_bb_part[ulii];
        }
      }
      BLAS_SET_NUM_THREADS(max_thread_ct);
      svd_rect_err = SvdRect(pca_sample_ct, qq_col_ct, svd_rect_lwork, bb, ss, svd_rect_wkspace);
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
      if (g_is_haploid) {
        for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
          eigvals[pc_idx] *= 0.5;
        }
      }
    } else {
      if (require_biallelic) {
        if (unlikely(MultiallelicVariantPresent(variant_include, allele_idx_offsets, variant_ct))) {
          logerrputs("Error: Multiallelic variant present in \"--pca biallelic-var-wts\" run.\n");
          goto CalcPca_ret_INCONSISTENT_INPUT;
        }
      }
      __CLPK_integer lwork;
      __CLPK_integer liwork;
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
      if (unlikely(
              bigstack_alloc_d(pc_ct * pca_sample_ct, &reverse_eigvecs_pcmaj) ||
              bigstack_alloc_uc(wkspace_byte_ct, &extract_eigvecs_wkspace))) {
        goto CalcPca_ret_NOMEM;
      }
      logprintf("Extracting eigenvalue%s and eigenvector%s... ", (pc_ct == 1)? "" : "s", (pc_ct == 1)? "" : "s");
      fflush(stdout);
      BLAS_SET_NUM_THREADS(max_thread_ct);
      // not putting unlikely() here for now.
      if (ExtractEigvecs(pca_sample_ct, pc_ct, lwork, liwork, grm, eigvals, reverse_eigvecs_pcmaj, extract_eigvecs_wkspace)) {
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

    if (var_wts) {
      g_g1 = eigvecs_smaj;
      for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
        eigval_inv_sqrts[pc_idx] = 1.0 / sqrt(eigvals[pc_idx]);
      }

      const uint32_t output_zst = (pca_flags / kfPcaVarZs) & 1;
      OutnameZstSet(".eigenvec.var", output_zst, outname_end);
      reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, writebuf, R_CAST(unsigned char*, &(writebuf[overflow_buf_size])), &css);
      if (unlikely(reterr)) {
        goto CalcPca_ret_1;
      }
      cswritep = writebuf;
      *cswritep++ = '#';
      if (chr_col) {
        cswritep = strcpya_k(cswritep, "CHROM\t");
      }
      if (pca_flags & kfPcaVcolPos) {
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
      if (maj_col) {
        cswritep = strcpya_k(cswritep, "\tMAJ");
      }
      if (nonmaj_col) {
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
        ts.calc_thread_ct = 1;
      }
#endif
      uintptr_t var_wts_part_size;
      if (qq) {
        var_wts_part_size = (MINV(variant_ct, calc_thread_ct * kPcaVariantBlockSize)) * S_CAST(uintptr_t, pc_ct);
      } else {
        // non-approximate PCA, bunch of buffers have not been allocated yet

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
        const uintptr_t genovecs_alloc = RoundUpPow2(pca_sample_ctaw2 * kPcaVariantBlockSize * sizeof(intptr_t), kCacheline);
        const uintptr_t dosage_cts_alloc = RoundUpPow2(kPcaVariantBlockSize * sizeof(int32_t), kCacheline);
        const uintptr_t dosage_presents_alloc = RoundUpPow2(pca_sample_ctaw * kPcaVariantBlockSize * sizeof(intptr_t), kCacheline);
        const uintptr_t dosage_main_alloc = RoundUpPow2(pca_sample_ct * kPcaVariantBlockSize * sizeof(Dosage), kCacheline);
        const uintptr_t cur_maj_freqs_alloc = RoundUpPow2(kPcaVariantBlockSize * sizeof(double), kCacheline);
        const uintptr_t yy_alloc = RoundUpPow2(kPcaVariantBlockSize * pca_sample_ct * sizeof(double), kCacheline);
        const uintptr_t per_thread_alloc = 2 * (genovecs_alloc + dosage_cts_alloc + dosage_presents_alloc + dosage_main_alloc + cur_maj_freqs_alloc) + yy_alloc + var_wts_part_alloc;
        if (per_thread_alloc * calc_thread_ct > arena_avail) {
          if (unlikely(arena_avail < per_thread_alloc)) {
            goto CalcPca_ret_NOMEM;
          }
          calc_thread_ct = arena_avail / per_thread_alloc;
        }
        ts.calc_thread_ct = calc_thread_ct;
        for (uint32_t parity = 0; parity != 2; ++parity) {
          g_genovecs[parity] = S_CAST(uintptr_t*, arena_alloc_raw(genovecs_alloc * calc_thread_ct, &arena_bottom));
          g_dosage_cts[parity] = S_CAST(uint32_t*, arena_alloc_raw(dosage_cts_alloc * calc_thread_ct, &arena_bottom));
          g_dosage_presents[parity] = S_CAST(uintptr_t*, arena_alloc_raw(dosage_presents_alloc * calc_thread_ct, &arena_bottom));
          g_dosage_mains[parity] = S_CAST(Dosage*, arena_alloc_raw(dosage_main_alloc * calc_thread_ct, &arena_bottom));
          g_cur_maj_freqs[parity] = S_CAST(double*, arena_alloc_raw(cur_maj_freqs_alloc * calc_thread_ct, &arena_bottom));
        }
        for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
          g_yy_bufs[tidx] = S_CAST(double*, arena_alloc_raw(yy_alloc, &arena_bottom));
        }
        var_wts_part_size = (MINV(variant_ct, calc_thread_ct * kPcaVariantBlockSize)) * S_CAST(uintptr_t, pc_ct);
        qq = S_CAST(double*, arena_alloc_raw_rd(2 * var_wts_part_size * sizeof(double), &arena_bottom));
        g_qq = qq;
#ifndef NDEBUG
        if (arena_top == g_bigstack_end) {
          // we shouldn't make any more allocations, but just in case...
          g_bigstack_base = arena_bottom;
        }
#endif
      }
      uint32_t prev_batch_size = 0;
      uintptr_t variant_uidx_load_base = 0;
      uintptr_t load_bits = variant_include[0];
      uintptr_t variant_uidx_write_base = 0;
      uintptr_t write_bits = variant_include[0];
      uint32_t parity = 0;
      ReinitThreads3z(&ts);
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_buf_blen = 0;
      for (uint32_t cur_variant_idx_start = 0; ; ) {
        uint32_t cur_batch_size = 0;
        if (!ts.is_last_block) {
          cur_batch_size = calc_thread_ct * kPcaVariantBlockSize;
          uint32_t cur_variant_idx_end = cur_variant_idx_start + cur_batch_size;
          if (cur_variant_idx_end > variant_ct) {
            cur_batch_size = variant_ct - cur_variant_idx_start;
            cur_variant_idx_end = variant_ct;
          }
          uintptr_t* genovec_iter = g_genovecs[parity];
          uint32_t* dosage_ct_iter = g_dosage_cts[parity];
          uintptr_t* dosage_present_iter = g_dosage_presents[parity];
          Dosage* dosage_main_iter = g_dosage_mains[parity];
          double* maj_freqs_write_iter = g_cur_maj_freqs[parity];
          for (uint32_t variant_idx = cur_variant_idx_start; variant_idx != cur_variant_idx_end; ++variant_idx) {
            const uintptr_t variant_uidx_load = BitIter1(variant_include, &variant_uidx_load_base, &load_bits);
            const uint32_t maj_allele_idx = maj_alleles[variant_uidx_load];
            uint32_t dosage_ct;
            reterr = PgrGetInv1D(pca_sample_include, pca_sample_include_cumulative_popcounts, pca_sample_ct, variant_uidx_load, maj_allele_idx, simple_pgrp, genovec_iter, dosage_present_iter, dosage_main_iter, &dosage_ct);
            if (unlikely(reterr)) {
              goto CalcPca_ret_PGR_FAIL;
            }
            ZeroTrailingQuaters(pca_sample_ct, genovec_iter);
            genovec_iter = &(genovec_iter[pca_sample_ctaw2]);
            *dosage_ct_iter++ = dosage_ct;
            dosage_present_iter = &(dosage_present_iter[pca_sample_ctaw]);
            dosage_main_iter = &(dosage_main_iter[pca_sample_ct]);
            uintptr_t allele_idx_base;
            if (!allele_idx_offsets) {
              allele_idx_base = variant_uidx_load;
            } else {
              allele_idx_base = allele_idx_offsets[variant_uidx_load];
              cur_allele_ct = allele_idx_offsets[variant_uidx_load + 1] - allele_idx_base;
              allele_idx_base -= variant_uidx_load;
            }
            *maj_freqs_write_iter++ = GetAlleleFreq(&(allele_freqs[allele_idx_base]), maj_allele_idx, cur_allele_ct);
          }
        }
        if (cur_variant_idx_start) {
          JoinThreads3z(&ts);
          if (unlikely(g_error_ret)) {
            goto CalcPca_ret_REWIND_FAIL;
          }
        }
        if (!ts.is_last_block) {
          g_cur_batch_size = cur_batch_size;
          ts.is_last_block = (cur_variant_idx_start + cur_batch_size == variant_ct);
          ts.thread_func_ptr = CalcPcaVarWtsThread;
          if (unlikely(SpawnThreads3z(cur_variant_idx_start, &ts))) {
            goto CalcPca_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (cur_variant_idx_start) {
          // write *previous* block results
          const double* var_wts_iter = &(qq[parity * var_wts_part_size]);
          // (todo: update projection here)
          for (uint32_t vidx = cur_variant_idx_start - prev_batch_size; vidx != cur_variant_idx_start; ++vidx) {
            const uint32_t variant_uidx_write = BitIter1(variant_include, &variant_uidx_write_base, &write_bits);
            if (chr_col) {
              // ok to skip this logic if chr_col not printed
              if (variant_uidx_write >= chr_end) {
                do {
                  ++chr_fo_idx;
                  chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
                } while (variant_uidx_write >= chr_end);
                const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
                char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
                *chr_name_end = '\t';
                chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
              }
              cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
            }
            if (variant_bps) {
              cswritep = u32toa_x(variant_bps[variant_uidx_write], '\t', cswritep);
            }
            cswritep = strcpya(cswritep, variant_ids[variant_uidx_write]);
            uintptr_t allele_idx_offset_base = variant_uidx_write * 2;
            if (allele_idx_offsets) {
              allele_idx_offset_base = allele_idx_offsets[variant_uidx_write];
              cur_allele_ct = allele_idx_offsets[variant_uidx_write + 1] - allele_idx_offset_base;
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
                  goto CalcPca_ret_WRITE_FAIL;
                }
                cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
              }
              --cswritep;
            }
            const uint32_t maj_allele_idx = maj_alleles[variant_uidx_write];
            if (maj_col) {
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto CalcPca_ret_WRITE_FAIL;
              }
              *cswritep++ = '\t';
              cswritep = strcpya(cswritep, cur_alleles[maj_allele_idx]);
            }
            if (nonmaj_col) {
              *cswritep++ = '\t';
              for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
                if (allele_idx == maj_allele_idx) {
                  continue;
                }
                if (unlikely(Cswrite(&css, &cswritep))) {
                  goto CalcPca_ret_WRITE_FAIL;
                }
                cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
              }
              --cswritep;
            }
            for (uint32_t pc_idx = 0; pc_idx != pc_ct; ++pc_idx) {
              *cswritep++ = '\t';
              // could avoid these multiplications by premultiplying the
              // sample weight matrix
              cswritep = dtoa_g((*var_wts_iter++) * eigval_inv_sqrts[pc_idx], cswritep);
            }
            AppendBinaryEoln(&cswritep);
            if (unlikely(Cswrite(&css, &cswritep))) {
              // bugfix (15 Dec 2017): prevent buffer overflow when ALT, MAJ,
              // and NONMAJ columns all missing.
              goto CalcPca_ret_WRITE_FAIL;
            }
          }
        }
        if (cur_variant_idx_start == variant_ct) {
          break;
        }
        cur_variant_idx_start += cur_batch_size;
        prev_batch_size = cur_batch_size;
      }
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto CalcPca_ret_WRITE_FAIL;
      }
      logprintfww("--pca%s: Variant weights written to %s .\n", is_approx? " approx" : "", outname);
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
  CalcPca_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
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
  CleanupThreads3z(&ts, &g_cur_batch_size);
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

// to test: do we actually want cur_dosage_ints to be uint64_t* instead of
// uint32_t*?
// also, should this be moved to plink2_common?
void FillCurDosageInts(const uintptr_t* genovec_buf, const uintptr_t* dosage_present, const Dosage* dosage_main_buf, uint32_t sample_ct, uint32_t dosage_ct, uint32_t is_diploid_p1, uint64_t* cur_dosage_ints) {
  uint64_t lookup_table[32] ALIGNV16;
  lookup_table[0] = 0;
  lookup_table[2] = is_diploid_p1 * kDosageMid;
  lookup_table[4] = is_diploid_p1 * kDosageMax;
  lookup_table[6] = 0;
  InitLookup16x8bx2(lookup_table);
  GenoarrLookup16x8bx2(genovec_buf, lookup_table, sample_ct, cur_dosage_ints);
  if (!dosage_ct) {
    return;
  }
  uintptr_t sample_idx_base = 0;
  uintptr_t cur_bits = dosage_present[0];
  for (uint32_t dosage_idx = 0; dosage_idx != dosage_ct; ++dosage_idx) {
    const uintptr_t sample_idx = BitIter1(dosage_present, &sample_idx_base, &cur_bits);
    cur_dosage_ints[sample_idx] = dosage_main_buf[dosage_idx] * is_diploid_p1;
  }
}

CONSTI32(kScoreVariantBlockSize, 240);
static double* g_dosages_vmaj[2] = {nullptr, nullptr};
static double* g_score_coefs_cmaj[2] = {nullptr, nullptr};
// don't bother to explicitly multithread for now
static double* g_final_scores_cmaj = nullptr;
static uint32_t g_score_col_ct = 0;
static uint32_t g_sample_ct = 0;

THREAD_FUNC_DECL CalcScoreThread(void* arg) {
  const uintptr_t tidx = R_CAST(uintptr_t, arg);
  assert(!tidx);
  double* final_scores_cmaj = g_final_scores_cmaj;
  const uint32_t score_col_ct = g_score_col_ct;
  const uint32_t sample_ct = g_sample_ct;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_batch = g_is_last_thread_block;
    const uint32_t cur_batch_size = g_cur_batch_size;
    if (cur_batch_size) {
      RowMajorMatrixMultiplyStridedIncr(g_score_coefs_cmaj[parity], g_dosages_vmaj[parity], score_col_ct, kScoreVariantBlockSize, sample_ct, sample_ct, cur_batch_size, sample_ct, final_scores_cmaj);
    }
    if (is_last_batch) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH_OLD(tidx);
    parity = 1 - parity;
  }
}

PglErr ScoreReport(const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const ChrInfo* cip, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const double* allele_freqs, const ScoreInfo* score_info_ptr, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t xchr_model, uint32_t max_thread_ct, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  char* cswritep = nullptr;
  PglErr reterr = kPglRetSuccess;
  TextStream score_txs;
  ThreadsState ts;
  CompressStreamState css;
  PreinitTextStream(&score_txs);
  InitThreads3z(&ts);
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

    const ScoreFlags score_flags = score_info_ptr->flags;
    reterr = SizeAndInitTextStream(score_info_ptr->input_fname, bigstack_left() / 8, 1, &score_txs);
    if (unlikely(reterr)) {
      goto ScoreReport_ret_TSTREAM_FAIL;
    }
    uint32_t nonempty_lines_to_skip_p1 = 1 + ((score_flags / kfScoreHeaderIgnore) & 1);
    char* line_start;
    for (uint32_t uii = 0; uii != nonempty_lines_to_skip_p1; ++uii) {
      reterr = TextNextNonemptyLineLstrip(&score_txs, &line_idx, &line_start);
      if (unlikely(reterr)) {
        if (reterr == kPglRetEof) {
          logerrputs("Error: Empty --score file.\n");
          goto ScoreReport_ret_MALFORMED_INPUT;
        }
        goto ScoreReport_ret_TSTREAM_FAIL;
      }
    }
    uint32_t last_col_idx = CountTokens(line_start);
    const uint32_t varid_col_idx = score_info_ptr->varid_col_p1 - 1;
    const uint32_t allele_col_idx = score_info_ptr->allele_col_p1 - 1;
    if (unlikely(MAXV(varid_col_idx, allele_col_idx) >= last_col_idx)) {
      goto ScoreReport_ret_MISSING_TOKENS;
    }
    uint32_t* score_col_idx_deltas = nullptr;
    uintptr_t score_col_ct = 1;
    if (!score_info_ptr->input_col_idx_range_list.name_ct) {
      if (unlikely(allele_col_idx == last_col_idx)) {
        goto ScoreReport_ret_MISSING_TOKENS;
      }
      if (unlikely(bigstack_alloc_u32(1, &score_col_idx_deltas))) {
        goto ScoreReport_ret_NOMEM;
      }
      // catch edge case
      if (unlikely(allele_col_idx + 1 == varid_col_idx)) {
        logerrputs("Error: --score variant ID column index matches a coefficient column index.\n");
        goto ScoreReport_ret_INVALID_CMDLINE;
      }
      score_col_idx_deltas[0] = allele_col_idx + 1;
    } else {
      const uint32_t last_col_idxl = BitCtToWordCt(last_col_idx);
      uintptr_t* score_col_bitarr;
      if (unlikely(bigstack_end_calloc_w(last_col_idxl, &score_col_bitarr))) {
        goto ScoreReport_ret_NOMEM;
      }
      if (unlikely(NumericRangeListToBitarr(&(score_info_ptr->input_col_idx_range_list), last_col_idx, 1, 0, score_col_bitarr))) {
        goto ScoreReport_ret_MISSING_TOKENS;
      }
      if (unlikely(IsSet(score_col_bitarr, varid_col_idx))) {
        logerrputs("Error: --score variant ID column index matches a coefficient column index.\n");
        goto ScoreReport_ret_INVALID_CMDLINE;
      }
      if (unlikely(IsSet(score_col_bitarr, allele_col_idx))) {
        logerrputs("Error: --score allele column index matches a coefficient column index.\n");
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
      BigstackEndReset(bigstack_end_mark);
    }
    char** score_col_names;
    if (unlikely(bigstack_alloc_cp(score_col_ct, &score_col_names))) {
      goto ScoreReport_ret_NOMEM;
    }
    char* write_iter = R_CAST(char*, g_bigstack_base);
    // don't have to worry about overflow, since linebuf was limited to 1/8
    // of available workspace.
    if (score_flags & kfScoreHeaderRead) {
      char* read_iter = line_start;
      for (uintptr_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
        read_iter = NextTokenMult0(read_iter, score_col_idx_deltas[score_col_idx]);
        score_col_names[score_col_idx] = write_iter;
        char* token_end = CurTokenEnd(read_iter);
        const uint32_t slen = token_end - read_iter;
        write_iter = memcpyax(write_iter, read_iter, slen, '\0');
      }

      // don't reparse this line
      line_start = K_CAST(char*, &(g_one_char_strs[0]));
    } else {
      for (uintptr_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
        score_col_names[score_col_idx] = write_iter;
        write_iter = strcpya_k(write_iter, "SCORE");
        write_iter = u32toa_x(score_col_idx + 1, '\0', write_iter);
      }
    }
    BigstackBaseSet(write_iter);

    g_score_col_ct = score_col_ct;
    g_sample_ct = sample_ct;
    g_cur_batch_size = kScoreVariantBlockSize;
    ts.calc_thread_ct = 1;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t sample_ctl2 = QuaterCtToWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t acc1_vec_ct = BitCtToVecCt(sample_ct);
    const uint32_t acc4_vec_ct = acc1_vec_ct * 4;
    const uint32_t acc8_vec_ct = acc1_vec_ct * 8;
    const uint32_t write_score_avgs = (score_flags / kfScoreColScoreAvgs) & 1;
    const uint32_t write_score_sums = (score_flags / kfScoreColScoreSums) & 1;
    const uintptr_t overflow_buf_size = RoundUpPow2((score_col_ct * (write_score_avgs + write_score_sums) + pheno_ct) * 16 + 3 * kMaxIdSlen + kCompressStreamBlock + 64, kCacheline);
    uintptr_t overflow_buf_alloc = overflow_buf_size;
    if (score_flags & (kfScoreZs | kfScoreListVariantsZs)) {
      overflow_buf_alloc += CstreamWkspaceReq(overflow_buf_size);
    }
    uint32_t* sample_include_cumulative_popcounts = nullptr;
    uintptr_t* sex_nonmale_collapsed = nullptr;
    uintptr_t* genovec_buf = nullptr;
    uintptr_t* dosage_present_buf = nullptr;
    Dosage* dosage_main_buf = nullptr;
    uintptr_t* missing_acc1 = nullptr;
    uintptr_t* missing_male_acc1 = nullptr;
    uint64_t* dosage_sums;
    uint64_t* dosage_incrs;
    uintptr_t* already_seen;
    char* overflow_buf = nullptr;
    if (unlikely(
            bigstack_alloc_thread(1, &ts.threads) ||
            bigstack_alloc_d((kScoreVariantBlockSize * k1LU) * sample_ct, &(g_dosages_vmaj[0])) ||
            bigstack_alloc_d((kScoreVariantBlockSize * k1LU) * sample_ct, &(g_dosages_vmaj[1])) ||
            bigstack_alloc_d(kScoreVariantBlockSize * score_col_ct, &(g_score_coefs_cmaj[0])) ||
            bigstack_alloc_d(kScoreVariantBlockSize * score_col_ct, &(g_score_coefs_cmaj[1])) ||
            bigstack_calloc_d(score_col_ct * sample_ct, &g_final_scores_cmaj) ||
            // bugfix (4 Nov 2017): need raw_sample_ctl here, not sample_ctl
            bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
            bigstack_alloc_w(sample_ctl, &sex_nonmale_collapsed) ||
            bigstack_alloc_w(sample_ctl2, &genovec_buf) ||
            bigstack_alloc_w(sample_ctl, &dosage_present_buf) ||
            bigstack_alloc_dosage(sample_ct, &dosage_main_buf) ||
            bigstack_alloc_w(45 * acc1_vec_ct * kWordsPerVec, &missing_acc1) ||
            bigstack_alloc_w(45 * acc1_vec_ct * kWordsPerVec, &missing_male_acc1) ||
            bigstack_calloc_u64(sample_ct, &dosage_sums) ||
            bigstack_calloc_u64(sample_ct, &dosage_incrs) ||
            bigstack_calloc_w(raw_variant_ctl, &already_seen) ||
            bigstack_alloc_c(overflow_buf_alloc, &overflow_buf))) {
      goto ScoreReport_ret_NOMEM;
    }
    VecW* missing_diploid_acc4 = &(R_CAST(VecW*, missing_acc1)[acc1_vec_ct]);
    VecW* missing_diploid_acc8 = &(missing_diploid_acc4[acc4_vec_ct]);
    VecW* missing_diploid_acc32 = &(missing_diploid_acc8[acc8_vec_ct]);
    VecW* missing_haploid_acc4 = &(R_CAST(VecW*, missing_male_acc1)[acc1_vec_ct]);
    VecW* missing_haploid_acc8 = &(missing_haploid_acc4[acc4_vec_ct]);
    VecW* missing_haploid_acc32 = &(missing_haploid_acc8[acc8_vec_ct]);
    ZeroVecArr(acc4_vec_ct, missing_diploid_acc4);
    ZeroVecArr(acc8_vec_ct, missing_diploid_acc8);
    ZeroVecArr(acc8_vec_ct * 4, missing_diploid_acc32);
    ZeroVecArr(acc4_vec_ct, missing_haploid_acc4);
    ZeroVecArr(acc8_vec_ct, missing_haploid_acc8);
    ZeroVecArr(acc8_vec_ct * 4, missing_haploid_acc32);
    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_nonmale_collapsed);
    AlignedBitarrInvert(sample_ct, sex_nonmale_collapsed);
    const uint32_t nonmale_ct = PopcountWords(sex_nonmale_collapsed, sample_ctl);
    const uint32_t male_ct = sample_ct - nonmale_ct;
    uint32_t* variant_id_htable = nullptr;
    uint32_t variant_id_htable_size;
    reterr = AllocAndPopulateIdHtableMt(variant_include, variant_ids, variant_ct, 0, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
    if (unlikely(reterr)) {
      goto ScoreReport_ret_1;
    }

    const uint32_t ignore_dup_ids = (score_flags / kfScoreIgnoreDupIds) & 1;
    const uint32_t list_variants = (score_flags / kfScoreListVariants) & 1;
    if (list_variants) {
      const uint32_t output_zst = (score_flags / kfScoreListVariantsZs) & 1;
      OutnameZstSet(".sscore.vars", output_zst, outname_end);
      reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
      if (unlikely(reterr)) {
        goto ScoreReport_ret_1;
      }
      cswritep = overflow_buf;
    }

    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    const uint32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    const uint32_t model_dominant = (score_flags / kfScoreDominant) & 1;
    const uint32_t domrec = model_dominant || (score_flags & kfScoreRecessive);
    const uint32_t variance_standardize = (score_flags / kfScoreVarianceStandardize) & 1;
    const uint32_t center = variance_standardize || (score_flags & kfScoreCenter);
    const uint32_t no_meanimpute = (score_flags / kfScoreNoMeanimpute) & 1;
    const uint32_t se_mode = (score_flags / kfScoreSe) & 1;
    uint32_t block_vidx = 0;
    uint32_t parity = 0;
    uint32_t cur_allele_ct = 2;
    double* cur_dosages_vmaj_iter = g_dosages_vmaj[0];
    double* cur_score_coefs_cmaj = g_score_coefs_cmaj[0];
    double geno_slope = kRecipDosageMax;
    double geno_intercept = 0.0;
    double cur_allele_freq = 0.0;
    uint32_t variant_ct_rem15 = 15;
    uint32_t variant_ct_rem255d15 = 17;
    uint32_t variant_hap_ct_rem15 = 15;
    uint32_t variant_hap_ct_rem255d15 = 17;
    uint32_t allele_ct_base = 0;
    int32_t male_allele_ct_delta = 0;
    uint32_t valid_variant_ct = 0;
    uintptr_t missing_var_id_ct = 0;
    uintptr_t duplicated_var_id_ct = 0;
    uintptr_t missing_allele_code_ct = 0;
#ifdef USE_MTBLAS
    const uint32_t matrix_multiply_thread_ct = (max_thread_ct > 1)? (max_thread_ct - 1) : 1;
    BLAS_SET_NUM_THREADS(matrix_multiply_thread_ct);
#endif
    PgrClearLdCache(simple_pgrp);
    while (1) {
      if (!IsEolnKns(*line_start)) {
        // varid_col_idx and allele_col_idx will almost always be very small
        char* variant_id_start = NextTokenMult0(line_start, varid_col_idx);
        if (unlikely(!variant_id_start)) {
          goto ScoreReport_ret_MISSING_TOKENS;
        }
        char* variant_id_token_end = CurTokenEnd(variant_id_start);
        const uint32_t variant_id_slen = variant_id_token_end - variant_id_start;
        uint32_t variant_uidx = VariantIdDupflagHtableFind(variant_id_start, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
        if (!(variant_uidx >> 31)) {
          if (unlikely(IsSet(already_seen, variant_uidx))) {
            // todo for first alpha 3 build: Support --pca multiallelic
            // .var.wts output.  Scores for a single variant can be split
            // across multiple lines, as long as they're all consecutive.
            // Variant isn't loaded from .pgen and scored until we've seen the
            // next variant ID, or have reached EOF.
            snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --score file.\n", variant_ids[variant_uidx]);
            goto ScoreReport_ret_MALFORMED_INPUT_WW;
          }
          SetBit(variant_uidx, already_seen);
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
          char allele_end_char = *allele_end;
          *allele_end = '\0';
          const uint32_t allele_blen = 1 + S_CAST(uintptr_t, allele_end - allele_start);
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          uint32_t cur_allele_idx = 0;
          for (; cur_allele_idx != cur_allele_ct; ++cur_allele_idx) {
            if (memequal(allele_start, cur_alleles[cur_allele_idx], allele_blen)) {
              break;
            }
          }
          if (cur_allele_idx != cur_allele_ct) {
            // okay, the variant and allele are in our dataset.  Load it.
            uint32_t dosage_ct;
            reterr = PgrGet1D(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, cur_allele_idx, simple_pgrp, genovec_buf, dosage_present_buf, dosage_main_buf, &dosage_ct);
            if (unlikely(reterr)) {
              goto ScoreReport_ret_PGR_FAIL;
            }
            const uint32_t chr_idx = GetVariantChr(cip, variant_uidx);
            uint32_t is_nonx_haploid = IsSet(cip->haploid_mask, chr_idx);
            if (unlikely(domrec && is_nonx_haploid)) {
              logerrputs("Error: --score 'dominant' and 'recessive' modifiers cannot be used with haploid\nchromosomes.\n");
              goto ScoreReport_ret_INCONSISTENT_INPUT;
            }
            uint32_t is_relevant_x = (chr_idx == x_code);
            if (unlikely(variance_standardize && (is_relevant_x || (chr_idx == mt_code)))) {
              logerrputs("Error: --score 'variance-standardize' cannot be used with chrX or MT.\n");
              goto ScoreReport_ret_INCONSISTENT_INPUT;
            }
            is_nonx_haploid = (!is_relevant_x) && is_nonx_haploid;

            // only if --xchr-model 1 (which is no longer the default)
            is_relevant_x = is_relevant_x && xchr_model;

            const uint32_t is_y = (chr_idx == y_code);
            ZeroTrailingQuaters(sample_ct, genovec_buf);
            GenovecToMissingnessUnsafe(genovec_buf, sample_ct, missing_acc1);
            if (dosage_ct) {
              BitvecInvmask(dosage_present_buf, sample_ctl, missing_acc1);
            }
            FillCurDosageInts(genovec_buf, dosage_present_buf, dosage_main_buf, sample_ct, dosage_ct, 2 - is_nonx_haploid, dosage_incrs);
            double ploidy_d;
            if (is_nonx_haploid) {
              if (is_y) {
                uintptr_t sample_idx_base = 0;
                uintptr_t sex_nonmale_collapsed_bits = sex_nonmale_collapsed[0];
                for (uint32_t nonmale_idx = 0; nonmale_idx != nonmale_ct; ++nonmale_idx) {
                  const uintptr_t sample_idx = BitIter1(sex_nonmale_collapsed, &sample_idx_base, &sex_nonmale_collapsed_bits);
                  dosage_incrs[sample_idx] = 0;
                }
                ++male_allele_ct_delta;
                BitvecInvmask(sex_nonmale_collapsed, sample_ctl, missing_acc1);
              } else {
                ++allele_ct_base;
              }
              VcountIncr1To4(missing_acc1, acc1_vec_ct, missing_haploid_acc4);
              if (!(--variant_hap_ct_rem15)) {
                Vcount0Incr4To8(acc4_vec_ct, missing_haploid_acc4, missing_haploid_acc8);
                variant_hap_ct_rem15 = 15;
                if (!(--variant_hap_ct_rem255d15)) {
                  Vcount0Incr8To32(acc8_vec_ct, missing_haploid_acc8, missing_haploid_acc32);
                  variant_hap_ct_rem255d15 = 17;
                }
              }
              if (is_y) {
                memcpy(missing_male_acc1, missing_acc1, sample_ctl * sizeof(intptr_t));
                BitvecOr(sex_nonmale_collapsed, sample_ctl, missing_acc1);
              }
              ploidy_d = 1.0;
            } else {
              if (is_relevant_x) {
                uintptr_t sample_idx_base = 0;
                uintptr_t sex_nonmale_collapsed_inv_bits = ~sex_nonmale_collapsed[0];
                for (uint32_t male_idx = 0; male_idx != male_ct; ++male_idx) {
                  const uintptr_t sample_idx = BitIter0(sex_nonmale_collapsed, &sample_idx_base, &sex_nonmale_collapsed_inv_bits);
                  dosage_incrs[sample_idx] /= 2;
                }
                BitvecInvmaskCopy(missing_acc1, sex_nonmale_collapsed, sample_ctl, missing_male_acc1);
                BitvecAnd(sex_nonmale_collapsed, sample_ctl, missing_acc1);
              }
              VcountIncr1To4(missing_acc1, acc1_vec_ct, missing_diploid_acc4);
              if (!(--variant_ct_rem15)) {
                Vcount0Incr4To8(acc4_vec_ct, missing_diploid_acc4, missing_diploid_acc8);
                variant_ct_rem15 = 15;
                if (!(--variant_ct_rem255d15)) {
                  Vcount0Incr8To32(acc8_vec_ct, missing_diploid_acc8, missing_diploid_acc32);
                  variant_ct_rem255d15 = 17;
                }
              }
              allele_ct_base += 2;
              if (is_relevant_x) {
                --male_allele_ct_delta;
                VcountIncr1To4(missing_male_acc1, acc1_vec_ct, missing_haploid_acc4);
                if (!(--variant_hap_ct_rem15)) {
                  Vcount0Incr4To8(acc4_vec_ct, missing_haploid_acc4, missing_haploid_acc8);
                  variant_hap_ct_rem15 = 15;
                  if (!(--variant_hap_ct_rem255d15)) {
                    Vcount0Incr8To32(acc8_vec_ct, missing_haploid_acc8, missing_haploid_acc32);
                    variant_hap_ct_rem255d15 = 17;
                  }
                }
                BitvecOr(missing_male_acc1, sample_ctl, missing_acc1);
              }
              if (!domrec) {
                ploidy_d = 2.0;
              } else {
                if (model_dominant) {
                  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
                    if (dosage_incrs[sample_idx] > kDosageMax) {
                      dosage_incrs[sample_idx] = kDosageMax;
                    }
                  }
                } else {
                  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
                    uint64_t cur_dosage_incr = dosage_incrs[sample_idx];
                    if (cur_dosage_incr <= kDosageMax) {
                      cur_dosage_incr = 0;
                    } else {
                      cur_dosage_incr -= kDosageMax;
                    }
                    dosage_incrs[sample_idx] = cur_dosage_incr;
                  }
                }
                ploidy_d = 1.0;
              }
            }
            for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
              dosage_sums[sample_idx] += dosage_incrs[sample_idx];
            }
            if (allele_freqs) {
              cur_allele_freq = GetAlleleFreq(&(allele_freqs[allele_idx_offset_base - variant_uidx]), cur_allele_idx, cur_allele_ct);
            }
            if (center) {
              if (variance_standardize) {
                const double variance = ploidy_d * cur_allele_freq * (1.0 - cur_allele_freq);
                if (variance < kSmallEpsilon) {
                  // ZeroTrailingQuaters(sample_ct, genovec_buf);
                  STD_ARRAY_DECL(uint32_t, 4, genocounts);
                  GenovecCountFreqsUnsafe(genovec_buf, sample_ct, genocounts);
                  if (unlikely(dosage_ct || genocounts[1] || genocounts[2])) {
                    snprintf(g_logbuf, kLogbufSize, "Error: --score variance-standardize failure for ID '%s': estimated allele frequency is zero, but not all dosages are zero. (This is possible when e.g. allele frequencies are estimated from founders, but the allele is only observed in nonfounders.)\n", variant_ids[variant_uidx]);
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
            const uint32_t missing_ct = PopcountWords(missing_acc1, sample_ctl);
            const uint32_t nm_sample_ct = sample_ct - missing_ct;
            if (missing_ct) {
              double missing_effect = 0.0;
              if (!no_meanimpute) {
                missing_effect = kDosageMax * cur_allele_freq * geno_slope;
              }
              uintptr_t sample_idx_base = 0;
              if (is_y || is_relevant_x) {
                ZeroDArr(sample_ct, cur_dosages_vmaj_iter);
                if (!no_meanimpute) {
                  const uint32_t male_missing_ct = PopcountWords(missing_male_acc1, sample_ctl);
                  uintptr_t missing_male_acc1_bits = missing_male_acc1[0];
                  for (uint32_t male_missing_idx = 0; male_missing_idx != male_missing_ct; ++male_missing_idx) {
                    const uintptr_t sample_idx = BitIter1(missing_male_acc1, &sample_idx_base, &missing_male_acc1_bits);
                    cur_dosages_vmaj_iter[sample_idx] = missing_effect;
                  }
                  if (is_relevant_x) {
                    // missing_male_acc1 not used after this point, so okay to
                    // use buffer for nonmales
                    BitvecAndCopy(missing_acc1, sex_nonmale_collapsed, sample_ctl, missing_male_acc1);
                    missing_effect *= 2;
                    // bugfix (8 Jul 2018): need to reset sample_idx
                    sample_idx_base = 0;
                    missing_male_acc1_bits = missing_male_acc1[0];
                    const uint32_t nonmale_missing_ct = PopcountWords(missing_male_acc1, sample_ctl);
                    for (uint32_t nonmale_missing_idx = 0; nonmale_missing_idx != nonmale_missing_ct; ++nonmale_missing_idx) {
                      const uintptr_t sample_idx = BitIter1(missing_male_acc1, &sample_idx_base, &missing_male_acc1_bits);
                      cur_dosages_vmaj_iter[sample_idx] = missing_effect;
                    }
                  }
                }
              } else {
                missing_effect *= ploidy_d;
                uintptr_t missing_acc1_bits = missing_acc1[0];
                for (uint32_t missing_idx = 0; missing_idx != missing_ct; ++missing_idx) {
                  const uintptr_t sample_idx = BitIter1(missing_acc1, &sample_idx_base, &missing_acc1_bits);
                  cur_dosages_vmaj_iter[sample_idx] = missing_effect;
                }
              }
            }
            uintptr_t sample_idx_base = 0;
            uintptr_t missing_acc1_inv_bits = ~missing_acc1[0];
            for (uint32_t nm_sample_idx = 0; nm_sample_idx != nm_sample_ct; ++nm_sample_idx) {
              const uintptr_t sample_idx = BitIter0(missing_acc1, &sample_idx_base, &missing_acc1_inv_bits);
              cur_dosages_vmaj_iter[sample_idx] = u63tod(dosage_incrs[sample_idx]) * geno_slope + geno_intercept;
            }
            if (se_mode) {
              // Suppose our score coefficients are drawn from independent
              // Gaussians.  Then the variance of the final score average is
              // the sum of the variances of the individual terms, divided by
              // (T^2) where T is the number of terms.  These individual
              // variances are of the form (<genotype value> * <stdev>)^2.
              //
              // Thus, we can use the same inner loop to compute standard
              // errors, as long as
              //   1. we square the genotypes and the standard errors before
              //      matrix multiplication, and
              //   2. we take the square root of the sums at the end.
              for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
                cur_dosages_vmaj_iter[sample_idx] *= cur_dosages_vmaj_iter[sample_idx];
              }
            }
            cur_dosages_vmaj_iter = &(cur_dosages_vmaj_iter[sample_ct]);

            *allele_end = allele_end_char;
            double* cur_score_coefs_iter = &(cur_score_coefs_cmaj[block_vidx]);
            const char* read_iter = line_start;
            for (uint32_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
              read_iter = NextTokenMult0(read_iter, score_col_idx_deltas[score_col_idx]);
              if (unlikely(!read_iter)) {
                goto ScoreReport_ret_MISSING_TOKENS;
              }
              double raw_coef;
              const char* token_end = ScanadvDouble(read_iter, &raw_coef);
              if (unlikely(!token_end)) {
                snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --score file has an invalid coefficient.\n", line_idx);
                goto ScoreReport_ret_MALFORMED_INPUT_2;
              }
              *cur_score_coefs_iter = raw_coef;
              cur_score_coefs_iter = &(cur_score_coefs_iter[kScoreVariantBlockSize]);
              read_iter = token_end;
            }
            if (list_variants) {
              cswritep = strcpya(cswritep, variant_ids[variant_uidx]);
              AppendBinaryEoln(&cswritep);
              if (unlikely(Cswrite(&css, &cswritep))) {
                goto ScoreReport_ret_WRITE_FAIL;
              }
            }
            ++valid_variant_ct;
            if (!(valid_variant_ct % 10000)) {
              printf("\r--score: %uk variants loaded.", valid_variant_ct / 1000);
              fflush(stdout);
            }
            ++block_vidx;
            if (block_vidx == kScoreVariantBlockSize) {
              if (se_mode) {
                for (uintptr_t ulii = 0; ulii != kScoreVariantBlockSize * score_col_ct; ++ulii) {
                  cur_score_coefs_cmaj[ulii] *= cur_score_coefs_cmaj[ulii];
                }
              }
              parity = 1 - parity;
              const uint32_t is_not_first_block = (ts.thread_func_ptr != nullptr);
              if (is_not_first_block) {
                JoinThreads3z(&ts);
                // CalcScoreThread() never errors out
              } else {
                ts.thread_func_ptr = CalcScoreThread;
              }
              if (unlikely(SpawnThreads3z(is_not_first_block, &ts))) {
                goto ScoreReport_ret_THREAD_CREATE_FAIL;
              }
              cur_dosages_vmaj_iter = g_dosages_vmaj[parity];
              cur_score_coefs_cmaj = g_score_coefs_cmaj[parity];
              block_vidx = 0;
            }
          } else {
            ++missing_allele_code_ct;
          }
        } else {
          ++missing_var_id_ct;
          if (variant_uidx != UINT32_MAX) {
            if (unlikely(!ignore_dup_ids)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --score variant ID '%s' appears multiple times in main dataset.\n", variant_ids[variant_uidx & 0x7fffffff]);
              goto ScoreReport_ret_INCONSISTENT_INPUT_WW;
            }
            ++duplicated_var_id_ct;
          }
        }
      }
      ++line_idx;
      reterr = TextNextLineLstrip(&score_txs, &line_start);
      if (reterr) {
        if (likely(reterr == kPglRetEof)) {
          reterr = kPglRetSuccess;
          break;
        }
        goto ScoreReport_ret_TSTREAM_FAIL;
      }
    }
    VcountIncr4To8(missing_diploid_acc4, acc4_vec_ct, missing_diploid_acc8);
    VcountIncr8To32(missing_diploid_acc8, acc8_vec_ct, missing_diploid_acc32);
    VcountIncr4To8(missing_haploid_acc4, acc4_vec_ct, missing_haploid_acc8);
    VcountIncr8To32(missing_haploid_acc8, acc8_vec_ct, missing_haploid_acc32);
    const uint32_t is_not_first_block = (ts.thread_func_ptr != nullptr);
    putc_unlocked('\r', stdout);
    if (missing_var_id_ct || missing_allele_code_ct || duplicated_var_id_ct) {
      if (!missing_var_id_ct) {
        snprintf(g_logbuf, kLogbufSize, "Warning: %" PRIuPTR " --score file entr%s.\n", missing_allele_code_ct, (missing_allele_code_ct == 1)? "y was skipped due to a mismatching allele code" : "ies were skipped due to mismatching allele codes");
      } else if (!missing_allele_code_ct) {
        snprintf(g_logbuf, kLogbufSize, "Warning: %" PRIuPTR " --score file entr%s.\n", missing_var_id_ct, (missing_var_id_ct == 1)? "y was skipped due to a missing variant ID" : "ies were skipped due to missing variant IDs");
      } else {
        snprintf(g_logbuf, kLogbufSize, "Warning: %" PRIuPTR " --score file entr%s, and %" PRIuPTR " %s.\n", missing_var_id_ct, (missing_var_id_ct == 1)? "y was skipped due to a missing variant ID" : "ies were skipped due to missing variant IDs", missing_allele_code_ct, (missing_allele_code_ct == 1)? "was skipped due to a mismatching allele code" : "were skipped due to mismatching allele codes");
      }
      WordWrapB(0);
      logerrputsb();
      if (duplicated_var_id_ct) {
        logerrprintfww("Warning: %" PRIuPTR " --score file entr%s appear multiple times in the main dataset.\n", duplicated_var_id_ct, (duplicated_var_id_ct == 1)? "y was skipped since its variant ID" : "ies were skipped since their variant IDs");
      }
      if (!list_variants) {
        logerrputs("(Add the 'list-variants' modifier to see which variants were actually used for\nscoring.)\n");
      }
    }
    if (block_vidx) {
      if (is_not_first_block) {
        JoinThreads3z(&ts);
      } else {
        ts.thread_func_ptr = CalcScoreThread;
      }
    } else if (unlikely(!valid_variant_ct)) {
      logerrputs("Error: No valid variants in --score file.\n");
      goto ScoreReport_ret_DEGENERATE_DATA;
    } else {
      JoinThreads3z(&ts);
    }
    ts.is_last_block = 1;
    g_cur_batch_size = block_vidx;
    if (se_mode) {
      for (uintptr_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
        double* cur_score_coefs_row = &(cur_score_coefs_cmaj[score_col_idx * kScoreVariantBlockSize]);
        for (uint32_t uii = 0; uii != block_vidx; ++uii) {
          cur_score_coefs_row[uii] *= cur_score_coefs_row[uii];
        }
      }
    }
    if (unlikely(SpawnThreads3z(is_not_first_block, &ts))) {
      goto ScoreReport_ret_THREAD_CREATE_FAIL;
    }
    JoinThreads3z(&ts);
    if (se_mode) {
      // sample_ct * score_col_ct
      for (uintptr_t ulii = 0; ulii != sample_ct * score_col_ct; ++ulii) {
        g_final_scores_cmaj[ulii] = sqrt(g_final_scores_cmaj[ulii]);
      }
    }
    logprintf("--score: %u variant%s processed.\n", valid_variant_ct, (valid_variant_ct == 1)? "" : "s");
    if (list_variants) {
      if (unlikely(CswriteCloseNull(&css, cswritep))) {
        goto ScoreReport_ret_WRITE_FAIL;
      }
      cswritep = nullptr;
      logprintf("Variant list written to %s .\n", outname);
    }

    const uint32_t output_zst = (score_flags / kfScoreZs) & 1;
    OutnameZstSet(".sscore", output_zst, outname_end);
    reterr = InitCstream(outname, 0, output_zst, max_thread_ct, overflow_buf_size, overflow_buf, R_CAST(unsigned char*, &(overflow_buf[overflow_buf_size])), &css);
    if (unlikely(reterr)) {
      goto ScoreReport_ret_1;
    }
    cswritep = overflow_buf;
    const uint32_t write_fid = FidColIsRequired(siip, score_flags / kfScoreColMaybefid);
    const char* sample_ids = siip->sample_ids;
    const char* sids = siip->sids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    const uint32_t write_sid = SidColIsRequired(sids, score_flags / kfScoreColMaybesid);
    const uint32_t write_empty_pheno = (score_flags & kfScoreColPheno1) && (!pheno_ct);
    const uint32_t write_phenos = (score_flags & (kfScoreColPheno1 | kfScoreColPhenos)) && pheno_ct;
    if (write_phenos && (!(score_flags & kfScoreColPhenos))) {
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
    const uint32_t write_nmiss_allele = (score_flags / kfScoreColNmissAllele) & 1;
    if (write_nmiss_allele) {
      cswritep = strcpya_k(cswritep, "\tNMISS_ALLELE_CT");
    }
    const uint32_t write_denom = (score_flags / kfScoreColDenom) & 1;
    if (write_denom) {
      cswritep = strcpya_k(cswritep, "\tDENOM");
    }
    const uint32_t write_dosage_sum = (score_flags / kfScoreColDosageSum) & 1;
    if (write_dosage_sum) {
      cswritep = strcpya_k(cswritep, "\tNAMED_ALLELE_DOSAGE_SUM");
    }
    if (write_score_avgs) {
      for (uint32_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, score_col_names[score_col_idx]);
        cswritep = strcpya_k(cswritep, "_AVG");
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto ScoreReport_ret_WRITE_FAIL;
        }
      }
    }
    if (write_score_sums) {
      for (uint32_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
        *cswritep++ = '\t';
        cswritep = strcpya(cswritep, score_col_names[score_col_idx]);
        cswritep = strcpya_k(cswritep, "_SUM");
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto ScoreReport_ret_WRITE_FAIL;
        }
      }
    }
    AppendBinaryEoln(&cswritep);
    const uint32_t* scrambled_missing_diploid_cts = R_CAST(uint32_t*, missing_diploid_acc32);
    const uint32_t* scrambled_missing_haploid_cts = R_CAST(uint32_t*, missing_haploid_acc32);
    const char* output_missing_pheno = g_output_missing_pheno;
    const uint32_t omp_slen = strlen(output_missing_pheno);

    uintptr_t sample_uidx_base = 0;
    uintptr_t sample_include_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
      const char* cur_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
      if (!write_fid) {
        cur_sample_id = AdvPastDelim(cur_sample_id, '\t');
      }
      cswritep = strcpya(cswritep, cur_sample_id);
      if (write_sid) {
        *cswritep++ = '\t';
        if (sids) {
          cswritep = strcpya(cswritep, &(sids[max_sid_blen * sample_uidx]));
        } else {
          *cswritep++ = '0';
        }
      }
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
      const uint32_t scrambled_idx = VcountScramble1(sample_idx);
      uint32_t denom = allele_ct_base + IsSet(sex_male, sample_uidx) * male_allele_ct_delta;
      const uint32_t nmiss_allele_ct = denom - 2 * scrambled_missing_diploid_cts[scrambled_idx] - scrambled_missing_haploid_cts[scrambled_idx];
      if (write_nmiss_allele) {
        *cswritep++ = '\t';
        cswritep = u32toa(nmiss_allele_ct, cswritep);
      }
      if (no_meanimpute) {
        denom = nmiss_allele_ct;
      }
      if (write_denom) {
        *cswritep++ = '\t';
        cswritep = u32toa(denom, cswritep);
      }
      if (write_dosage_sum) {
        *cswritep++ = '\t';
        cswritep = dosagetoa(dosage_sums[sample_idx], cswritep);
      }
      const double* final_score_col = &(g_final_scores_cmaj[sample_idx]);
      if (write_score_avgs) {
        const double denom_recip = 1.0 / S_CAST(double, denom);
        for (uint32_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(final_score_col[score_col_idx * sample_ct] * denom_recip, cswritep);
        }
      }
      if (write_score_sums) {
        for (uint32_t score_col_idx = 0; score_col_idx != score_col_ct; ++score_col_idx) {
          *cswritep++ = '\t';
          cswritep = dtoa_g(final_score_col[score_col_idx * sample_ct], cswritep);
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
    logprintfww("--score: Results written to %s .\n", outname);
  }
  while (0) {
  ScoreReport_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--score file", &score_txs);
    break;
  ScoreReport_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ScoreReport_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ScoreReport_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ScoreReport_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  ScoreReport_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
  ScoreReport_ret_MALFORMED_INPUT_2:
    logputs("\n");
    logerrputsb();
  ScoreReport_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ScoreReport_ret_MISSING_TOKENS:
    logputs("\n");
    logerrprintfww("Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, score_info_ptr->input_fname);
    reterr = kPglRetInconsistentInput;
    break;
  ScoreReport_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logputs("\n");
    logerrputsb();
  ScoreReport_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ScoreReport_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  ScoreReport_ret_DEGENERATE_DATA_WW:
    WordWrapB(0);
    logputs("\n");
    logerrputsb();
  ScoreReport_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 ScoreReport_ret_1:
  CswriteCloseCond(&css, cswritep);
  CleanupThreads3z(&ts, &g_cur_batch_size);
  BLAS_SET_NUM_THREADS(1);
  CleanupTextStream2("--score file", &score_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
