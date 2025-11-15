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

#include "plink2_filter.h"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <string.h>

#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_stats.h"  // HweThresh(), etc.
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "plink2_decompress.h"
#include "plink2_random.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitExtractColCond(ExtractColCondInfo* eccip) {
  eccip->params = nullptr;
  eccip->match_flattened = nullptr;
  eccip->mismatch_flattened = nullptr;
  eccip->match_substr = 0;
  eccip->min = 0.0;
  eccip->max = DBL_MAX;
}

void CleanupExtractColCond(ExtractColCondInfo* eccip) {
  free_cond(eccip->params);
  free_cond(eccip->match_flattened);
  free_cond(eccip->mismatch_flattened);
}

PglErr FromToFlag(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* varid_from, const char* varid_to, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uintptr_t* variant_include, ChrInfo* cip, uint32_t* variant_ct_ptr) {
  PglErr reterr = kPglRetSuccess;
  {
    if (!(*variant_ct_ptr)) {
      goto FromToFlag_ret_1;
    }
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t variant_uidx_start = UINT32_MAX;
    if (varid_from) {
      uint32_t cur_llidx;
      variant_uidx_start = VariantIdDupHtableFind(varid_from, variant_ids, variant_id_htable, htable_dup_base, strlen(varid_from), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
      if (unlikely(variant_uidx_start == UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --from variant '%s' not found.\n", varid_from);
        goto FromToFlag_ret_INCONSISTENT_INPUT_WW;
      }
      // do *not* check variant_include here.  variant ID uniqueness should not
      // be dependent on the order in which filters are applied.
      if (unlikely(cur_llidx != UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --from variant ID '%s' appears multiple times.\n", varid_from);
        goto FromToFlag_ret_INCONSISTENT_INPUT_WW;
      }
      chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx_start);
    }
    uint32_t variant_uidx_end = 0;
    if (varid_to) {
      uint32_t cur_llidx;
      variant_uidx_end = VariantIdDupHtableFind(varid_to, variant_ids, variant_id_htable, htable_dup_base, strlen(varid_to), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
      if (unlikely(variant_uidx_end == UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --to variant '%s' not found.\n", varid_to);
        goto FromToFlag_ret_INCONSISTENT_INPUT_WW;
      }
      if (unlikely(cur_llidx != UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --to variant ID '%s' appears multiple times.\n", varid_to);
        goto FromToFlag_ret_INCONSISTENT_INPUT_WW;
      }
      uint32_t chr_fo_idx2 = GetVariantChrFoIdx(cip, variant_uidx_end);
      if (variant_uidx_start == UINT32_MAX) {
        chr_fo_idx = chr_fo_idx2;
        variant_uidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
      } else {
        if (unlikely(chr_fo_idx != chr_fo_idx2)) {
          logerrputs("Error: --from and --to variants are not on the same chromosome.\n");
          goto FromToFlag_ret_INCONSISTENT_INPUT;
        }
        if (variant_uidx_start > variant_uidx_end) {
          // permit order to be reversed
          uint32_t uii = variant_uidx_start;
          variant_uidx_start = variant_uidx_end;
          variant_uidx_end = uii;
        }
      }
      ++variant_uidx_end;  // convert to half-open interval
    } else {
      variant_uidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
    }
    if (variant_uidx_start) {
      ClearBitsNz(0, variant_uidx_start, variant_include);
    }
    if (variant_uidx_end < raw_variant_ct) {
      ClearBitsNz(variant_uidx_end, raw_variant_ct, variant_include);
    }
    ZeroWArr(kChrMaskWords, cip->chr_mask);
    SetBit(cip->chr_file_order[chr_fo_idx], cip->chr_mask);
    const uint32_t new_variant_ct = PopcountBitRange(variant_include, variant_uidx_start, variant_uidx_end);
    logprintf("--from/--to: %u variant%s remaining.\n", new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  FromToFlag_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  FromToFlag_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 FromToFlag_ret_1:
  return reterr;
}

PglErr SnpFlag(const uint32_t* variant_bps, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* varid_snp, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, int32_t window_bp, uintptr_t* variant_include, ChrInfo* cip, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    if (!(*variant_ct_ptr)) {
      goto SnpFlag_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uint32_t cur_llidx;
    uint32_t variant_uidx = VariantIdDupHtableFind(varid_snp, variant_ids, variant_id_htable, htable_dup_base, strlen(varid_snp), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
    if (unlikely(variant_uidx == UINT32_MAX)) {
      snprintf(g_logbuf, kLogbufSize, "Error: --%ssnp variant '%s' not found.\n", do_exclude? "exclude-" : "", varid_snp);
      goto SnpFlag_ret_INCONSISTENT_INPUT_WW;
    }
    if (window_bp == -1) {
      // duplicates ok

      uintptr_t* seen_uidxs;
      // not actually necessary in --exclude-snp case, but this is still fast
      // enough relative to hash table construction that there's no point in
      // complicating the code further to conditionally optimize this out
      if (unlikely(bigstack_calloc_w(raw_variant_ctl, &seen_uidxs))) {
        goto SnpFlag_ret_NOMEM;
      }
      for (; ; cur_llidx = htable_dup_base[cur_llidx + 1]) {
        SetBit(variant_uidx, seen_uidxs);
        if (cur_llidx == UINT32_MAX) {
          break;
        }
        variant_uidx = htable_dup_base[cur_llidx];
      }
      if (do_exclude) {
        BitvecInvmask(seen_uidxs, raw_variant_ctl, variant_include);
      } else {
        BitvecAnd(seen_uidxs, raw_variant_ctl, variant_include);
      }
    } else {
      if (unlikely(cur_llidx != UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --%ssnp + --window central variant ID '%s' appears multiple times.\n", do_exclude? "exclude-" : "", varid_snp);
        goto SnpFlag_ret_INCONSISTENT_INPUT_WW;
      }
      const uint32_t chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
      const uint32_t chr_vidx_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      const uint32_t center_bp = variant_bps[variant_uidx];
      const uint32_t window_bp_u = window_bp;
      uint32_t vidx_start = cip->chr_fo_vidx_start[chr_fo_idx];
      if (center_bp > window_bp_u) {
        vidx_start = LowerBoundConstrainedNonemptyU32(variant_bps, vidx_start, chr_vidx_end, center_bp - window_bp_u);
      }
      const uint32_t bp_end = 1 + center_bp + window_bp_u;
      const uint32_t vidx_end = LowerBoundConstrainedNonemptyU32(variant_bps, vidx_start, chr_vidx_end, bp_end);
      if (do_exclude) {
        ClearBitsNz(vidx_start, vidx_end, variant_include);
      } else {
        if (vidx_start) {
          ClearBitsNz(0, vidx_start, variant_include);
        }
        if (vidx_end < raw_variant_ct) {
          ClearBitsNz(vidx_end, raw_variant_ct, variant_include);
        }
        ZeroWArr(kChrMaskWords, cip->chr_mask);
        SetBit(cip->chr_file_order[chr_fo_idx], cip->chr_mask);
      }
    }
    const uint32_t new_variant_ct = PopcountWords(variant_include, raw_variant_ctl);
    logprintf("--%ssnp%s: %u variant%s remaining.\n", do_exclude? "exclude-" : "", (window_bp == -1)? "" : " + --window", new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  SnpFlag_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SnpFlag_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 SnpFlag_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr InterpretVariantRangeList(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const RangeList* snps_range_list_ptr, const char* flagname, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uintptr_t* seen_uidxs) {
  PglErr reterr = kPglRetSuccess;
  {
    const char* varid_strbox = snps_range_list_ptr->names;
    const unsigned char* starts_range = snps_range_list_ptr->starts_range;
    const uint32_t varid_ct = snps_range_list_ptr->name_ct;
    const uintptr_t varid_max_blen = snps_range_list_ptr->name_max_blen;
    uint32_t range_start_vidx = UINT32_MAX;
    for (uint32_t varid_idx = 0; varid_idx != varid_ct; ++varid_idx) {
      const char* cur_varid = &(varid_strbox[varid_idx * varid_max_blen]);
      uint32_t cur_llidx;
      uint32_t variant_uidx = VariantIdDupHtableFind(cur_varid, variant_ids, variant_id_htable, htable_dup_base, strlen(cur_varid), variant_id_htable_size, max_variant_id_slen, &cur_llidx);
      if (unlikely(variant_uidx == UINT32_MAX)) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s variant '%s' not found.\n", flagname, cur_varid);
        goto InterpretVariantRangeList_ret_INCONSISTENT_INPUT_WW;
      }
      if (starts_range[varid_idx]) {
        if (unlikely(cur_llidx != UINT32_MAX)) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s range-starting variant ID '%s' appears multiple times.\n", flagname, cur_varid);
          goto InterpretVariantRangeList_ret_INCONSISTENT_INPUT_WW;
        }
        range_start_vidx = variant_uidx;
      } else {
        if (range_start_vidx != UINT32_MAX) {
          if (unlikely(cur_llidx != UINT32_MAX)) {
            snprintf(g_logbuf, kLogbufSize, "Error: %s range-ending variant ID '%s' appears multiple times.\n", flagname, cur_varid);
            goto InterpretVariantRangeList_ret_INCONSISTENT_INPUT_WW;
          }
          if (variant_uidx < range_start_vidx) {
            const uint32_t uii = variant_uidx;
            variant_uidx = range_start_vidx;
            range_start_vidx = uii;
          }
          FillBitsNz(range_start_vidx, variant_uidx + 1, seen_uidxs);
        } else {
          for (; ; cur_llidx = htable_dup_base[cur_llidx + 1]) {
            SetBit(variant_uidx, seen_uidxs);
            if (cur_llidx == UINT32_MAX) {
              break;
            }
            variant_uidx = htable_dup_base[cur_llidx];
          }
        }
        range_start_vidx = UINT32_MAX;
      }
    }
  }
  while (0) {
  InterpretVariantRangeList_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
  return reterr;
}

PglErr SnpsFlag(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const RangeList* snps_range_list_ptr, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t do_exclude, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    if (!(*variant_ct_ptr)) {
      goto SnpsFlag_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* seen_uidxs;
    if (unlikely(bigstack_calloc_w(raw_variant_ctl, &seen_uidxs))) {
      goto SnpsFlag_ret_NOMEM;
    }
    reterr = InterpretVariantRangeList(variant_ids, variant_id_htable, htable_dup_base, snps_range_list_ptr, do_exclude? "--exclude-snps" : "--snps", max_variant_id_slen, variant_id_htable_size, seen_uidxs);
    if (unlikely(reterr)) {
      goto SnpsFlag_ret_1;
    }
    if (do_exclude) {
      BitvecInvmask(seen_uidxs, raw_variant_ctl, variant_include);
    } else {
      BitvecAnd(seen_uidxs, raw_variant_ctl, variant_include);
    }
    const uint32_t new_variant_ct = PopcountWords(variant_include, raw_variant_ctl);
    logprintf("--%ssnps: %u variant%s remaining.\n", do_exclude? "exclude-" : "", new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  SnpsFlag_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 SnpsFlag_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

void ExtractExcludeProcessTokens(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* shard_start, const char* shard_end, uint32_t variant_id_htable_size, uint32_t max_variant_id_slen, uintptr_t* already_seen) {
  const char* shard_iter = shard_start;
  while (1) {
    shard_iter = FirstPostspaceBounded(shard_iter, shard_end);
    if (shard_iter == shard_end) {
      return;
    }
    const char* token_end = CurTokenEnd(shard_iter);
    uint32_t cur_llidx;
    uint32_t variant_uidx = VariantIdDupHtableFind(shard_iter, variant_ids, variant_id_htable, htable_dup_base, token_end - shard_iter, variant_id_htable_size, max_variant_id_slen, &cur_llidx);
    shard_iter = token_end;
    if (variant_uidx == UINT32_MAX) {
      continue;
    }
    if (IsSet(already_seen, variant_uidx)) {
      continue;
    }
    for (; ; cur_llidx = htable_dup_base[cur_llidx + 1]) {
      SetBit(variant_uidx, already_seen);
      if (cur_llidx == UINT32_MAX) {
        break;
      }
      variant_uidx = htable_dup_base[cur_llidx];
    }
  }
}

CONSTI32(kMaxExtractExcludeThreads, 8);

typedef struct ExtractExcludeCtxStruct {
  const char* const* variant_ids;
  const uint32_t* variant_id_htable;
  const uint32_t* htable_dup_base;
  uintptr_t variant_id_htable_size;
  uint32_t max_variant_id_slen;

  char* shard_boundaries[kMaxExtractExcludeThreads + 1];
  uintptr_t* already_seens[kMaxExtractExcludeThreads];
} ExtractExcludeCtx;

// probable todo: rename and move into plink2_cmdline (or plink2_htable if
// VariantIdDupHtableFind is also moved).
THREAD_FUNC_DECL ExtractExcludeThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx_p1 = arg->tidx + 1;
  ExtractExcludeCtx* ctx = S_CAST(ExtractExcludeCtx*, arg->sharedp->context);

  const char* const* variant_ids = ctx->variant_ids;
  const uint32_t* variant_id_htable = ctx->variant_id_htable;
  const uint32_t* htable_dup_base = ctx->htable_dup_base;
  const uintptr_t variant_id_htable_size = ctx->variant_id_htable_size;
  const uint32_t max_variant_id_slen = ctx->max_variant_id_slen;
  uintptr_t* already_seen = ctx->already_seens[tidx_p1];
  do {
    ExtractExcludeProcessTokens(variant_ids, variant_id_htable, htable_dup_base, ctx->shard_boundaries[tidx_p1], ctx->shard_boundaries[tidx_p1 + 1], variant_id_htable_size, max_variant_id_slen, already_seen);
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr TokenExtractExclude(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const char* fnames, const char* flagname, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, VfilterType vft, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  const char* fname_tks = nullptr;
  PglErr reterr = kPglRetSuccess;
  TokenStream tks;
  ThreadGroup tg;
  PreinitTokenStream(&tks);
  PreinitThreads(&tg);
  ExtractExcludeCtx ctx;
  {
    const uint32_t calc_thread_ct_m1 = MINV(max_thread_ct, kMaxExtractExcludeThreads) - 1;
    if (unlikely(SetThreadCt0(calc_thread_ct_m1, &tg))) {
      goto TokenExtractExclude_ret_NOMEM;
    }
    if (!(*variant_ct_ptr)) {
      goto TokenExtractExclude_ret_1;
    }
    uint32_t decompress_thread_ct = 1;
    if (max_thread_ct > calc_thread_ct_m1 + 2) {
      decompress_thread_ct = max_thread_ct - calc_thread_ct_m1 - 1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    for (uint32_t tidx = 0; tidx <= calc_thread_ct_m1; ++tidx) {
      if (unlikely(bigstack_calloc_w(raw_variant_ctl, &(ctx.already_seens[tidx])))) {
        goto TokenExtractExclude_ret_NOMEM;
      }
    }
    if (calc_thread_ct_m1) {
      ctx.variant_ids = variant_ids;
      ctx.variant_id_htable = variant_id_htable;
      ctx.htable_dup_base = htable_dup_base;
      ctx.variant_id_htable_size = variant_id_htable_size;
      ctx.max_variant_id_slen = max_variant_id_slen;
      SetThreadFuncAndData(ExtractExcludeThread, &ctx, &tg);
    }
    const char* fnames_iter = fnames;
    do {
      if (fnames_iter == fnames) {
        fname_tks = fnames_iter;
        reterr = InitTokenStream(fnames_iter, decompress_thread_ct, &tks);
        if (unlikely(reterr)) {
          goto TokenExtractExclude_ret_TKSTREAM_FAIL;
        }
      } else {
        reterr = TokenStreamRetarget(fnames_iter, &tks);
        if (unlikely(reterr)) {
          goto TokenExtractExclude_ret_TKSTREAM_FAIL;
        }
        fname_tks = fnames_iter;
      }
      while (1) {
        reterr = TksNext(&tks, calc_thread_ct_m1 + 1, ctx.shard_boundaries);
        if (reterr) {
          break;
        }
        if (calc_thread_ct_m1) {
          if (unlikely(SpawnThreads(&tg))) {
            goto TokenExtractExclude_ret_THREAD_CREATE_FAIL;
          }
        }
        ExtractExcludeProcessTokens(variant_ids, variant_id_htable, htable_dup_base, ctx.shard_boundaries[0], ctx.shard_boundaries[1], variant_id_htable_size, max_variant_id_slen, ctx.already_seens[0]);
        JoinThreads0(&tg);
      }
      if (unlikely(reterr != kPglRetEof)) {
        goto TokenExtractExclude_ret_TKSTREAM_FAIL;
      }
      if (vft == kVfilterExtractIntersect) {
        for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
          BitvecOr(ctx.already_seens[tidx], raw_variant_ctl, ctx.already_seens[0]);
          ZeroWArr(raw_variant_ctl, ctx.already_seens[tidx]);
        }
        BitvecAnd(ctx.already_seens[0], raw_variant_ctl, variant_include);
        ZeroWArr(raw_variant_ctl, ctx.already_seens[0]);
      }
      fnames_iter = strnul(fnames_iter);
      ++fnames_iter;
    } while (*fnames_iter);
    reterr = kPglRetSuccess;
    if (vft == kVfilterExclude) {
      for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
        BitvecOr(ctx.already_seens[tidx], raw_variant_ctl, ctx.already_seens[0]);
      }
      BitvecInvmask(ctx.already_seens[0], raw_variant_ctl, variant_include);
    } else if (vft == kVfilterExtract) {
      for (uint32_t tidx = 1; tidx <= calc_thread_ct_m1; ++tidx) {
        BitvecOr(ctx.already_seens[tidx], raw_variant_ctl, ctx.already_seens[0]);
      }
      BitvecAnd(ctx.already_seens[0], raw_variant_ctl, variant_include);
    }
    const uint32_t new_variant_ct = PopcountWords(variant_include, raw_variant_ctl);
    logprintf("--%s: %u variant%s remaining.\n", flagname, new_variant_ct, (new_variant_ct == 1)? "" : "s");
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  TokenExtractExclude_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TokenExtractExclude_ret_TKSTREAM_FAIL:
    TokenStreamErrPrint(fname_tks, &tks);
    break;
  TokenExtractExclude_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 TokenExtractExclude_ret_1:
  CleanupThreads(&tg);
  if (fname_tks) {
    CleanupTokenStream2(fname_tks, &tks, &reterr);
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr ExtractColCond(const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const ExtractColCondInfo* eccip, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t htable_size, uint32_t max_thread_ct, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    // Very similar to e.g. UpdateVarBps().
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* variant_include_new;
    uintptr_t* already_seen;
    if (unlikely(bigstack_calloc_w(raw_variant_ctl, &variant_include_new) ||
                 bigstack_calloc_w(raw_variant_ctl, &already_seen))) {
      goto ExtractColCond_ret_NOMEM;
    }
    char* match_strbox = nullptr;
    uint32_t match_str_ct = 0;
    uintptr_t max_match_blen = 0;
    // We don't actually need to sort in the substring-match case, but
    // additional cost is negligible.
    if (eccip->match_flattened) {
      if (unlikely(MultistrToStrboxDedupAlloc(eccip->match_flattened, &match_strbox, &match_str_ct, &max_match_blen))) {
        goto ExtractColCond_ret_NOMEM;
      }
    }
    char* mismatch_strbox = nullptr;
    uint32_t mismatch_str_ct = 0;
    uintptr_t max_mismatch_blen = 0;
    if (eccip->mismatch_flattened) {
      if (unlikely(MultistrToStrboxDedupAlloc(eccip->mismatch_flattened, &mismatch_strbox, &mismatch_str_ct, &max_mismatch_blen))) {
        goto ExtractColCond_ret_NOMEM;
      }
    }
    reterr = SizeAndInitTextStream(eccip->params->fname, bigstack_left(), MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto ExtractColCond_ret_TSTREAM_FAIL;
    }
    reterr = TextSkip(eccip->params->skip_ct, &txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: Fewer lines than expected in --extract-col-cond file.\n");
        goto ExtractColCond_ret_INCONSISTENT_INPUT;
      }
      goto ExtractColCond_ret_TSTREAM_FAIL;
    }
    line_idx = eccip->params->skip_ct;
    const uint32_t colid_first = (eccip->params->colid < eccip->params->colx);
    uint32_t colmin;
    uint32_t coldiff;
    if (colid_first) {
      colmin = eccip->params->colid - 1;
      coldiff = eccip->params->colx - eccip->params->colid;
    } else {
      colmin = eccip->params->colx - 1;
      coldiff = eccip->params->colid - eccip->params->colx;
    }
    const char skipchar = eccip->params->skipchar;
    const uint32_t match_substr = eccip->match_substr;
    double val_min = eccip->min;
    double val_max = eccip->max;
    uintptr_t miss_ct = 0;
    while (1) {
      ++line_idx;
      char* line_start = TextGet(&txs);
      if (!line_start) {
        if (likely(!TextStreamErrcode2(&txs, &reterr))) {
          break;
        }
        goto ExtractColCond_ret_TSTREAM_FAIL;
      }
      char cc = *line_start;
      if (cc == skipchar) {
        continue;
      }
      char* colid_ptr;
      char* colval_ptr;
      if (colid_first) {
        colid_ptr = NextTokenMult0(line_start, colmin);
        colval_ptr = NextTokenMult(colid_ptr, coldiff);
        if (unlikely(!colval_ptr)) {
          goto ExtractColCond_ret_MISSING_TOKENS;
        }
      } else {
        colval_ptr = NextTokenMult0(line_start, colmin);
        colid_ptr = NextTokenMult(colval_ptr, coldiff);
        if (unlikely(!colid_ptr)) {
          goto ExtractColCond_ret_MISSING_TOKENS;
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
      if (unlikely(IsSet(already_seen, variant_uidx))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --extract-col-cond file.\n", cur_var_id);
        goto ExtractColCond_ret_INCONSISTENT_INPUT_WW;
      }
      SetBit(variant_uidx, already_seen);
      if (mismatch_str_ct) {
        char* token_end = CurTokenEnd(colval_ptr);
        if (!match_substr) {
          const int32_t ii = bsearch_strbox(colval_ptr, mismatch_strbox, token_end - colval_ptr, max_mismatch_blen, mismatch_str_ct);
          if (ii != -1) {
            continue;
          }
        } else {
          // todo: benchmark memmem, use it when it's available and better
          *token_end = '\0';
          uint32_t uii = 0;
          for (; uii != mismatch_str_ct; ++uii) {
            if (strstr(colval_ptr, &(mismatch_strbox[uii * max_mismatch_blen]))) {
              break;
            }
          }
          if (uii != mismatch_str_ct) {
            continue;
          }
          if (match_str_ct) {
            for (uii = 0; uii != match_str_ct; ++uii) {
              if (strstr(colval_ptr, &(match_strbox[uii * max_match_blen]))) {
                break;
              }
            }
            if (uii == match_str_ct) {
              continue;
            }
          }
        }
      } else if (match_str_ct) {
        // --extract-col-cond-match
        char* token_end = CurTokenEnd(colval_ptr);
        if (!match_substr) {
          const int32_t ii = bsearch_strbox(colval_ptr, match_strbox, token_end - colval_ptr, max_match_blen, match_str_ct);
          if (ii == -1) {
            continue;
          }
        } else {
          *token_end = '\0';
          uint32_t uii = 0;
          for (; uii != match_str_ct; ++uii) {
            if (strstr(colval_ptr, &(match_strbox[uii * max_match_blen]))) {
              break;
            }
          }
          if (uii == match_str_ct) {
            continue;
          }
        }
      } else {
        // min-max
        double val;
        if ((!ScantokDouble(colval_ptr, &val)) || (val < val_min) || (val > val_max)) {
          continue;
        }
      }
      for (; ; cur_llidx = htable_dup_base[cur_llidx + 1]) {
        SetBit(variant_uidx, variant_include_new);
        if (cur_llidx == UINT32_MAX) {
          break;
        }
        variant_uidx = htable_dup_base[cur_llidx];
      }
    }
    BitvecAnd(variant_include_new, raw_variant_ctl, variant_include);
    const uint32_t new_variant_ct = PopcountWords(variant_include, raw_variant_ctl);
    if (miss_ct) {
      logprintfww("--extract-col-cond: %u variant%s remaining, %" PRIuPTR " ID%s missing.\n", new_variant_ct, (new_variant_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
    } else {
      logprintf("--extract-col-cond: %u variant%s remaining.\n", new_variant_ct, (new_variant_ct == 1)? "" : "s");
    }
    *variant_ct_ptr = new_variant_ct;
  }
  while (0) {
  ExtractColCond_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExtractColCond_ret_TSTREAM_FAIL:
    TextStreamErrPrint(eccip->params->fname, &txs);
    break;
  ExtractColCond_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of --extract-col-cond file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  ExtractColCond_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  ExtractColCond_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  CleanupTextStream2(eccip->params->fname, &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

// could permit split-chromosome here
PglErr RmDup(const uintptr_t* sample_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uint32_t* variant_id_htable, const uint32_t* htable_dup_base, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const char* missing_varid_match, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t max_variant_id_slen, uintptr_t variant_id_htable_size, uint32_t orig_dup_ct, RmDupMode rmdup_mode, uint32_t save_list, uint32_t max_thread_ct, PgenReader* simple_pgrp, uintptr_t* variant_include, uint32_t* variant_ct_ptr, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  TextStream pvar_txs;
  PreinitTextStream(&pvar_txs);
  FILE* list_file = nullptr;
  FILE* mismatch_file = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (!orig_dup_ct) {
      logputs("Note: Skipping --rm-dup since no duplicate IDs are present.\n");
      goto RmDup_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uintptr_t* orig_dups;
    uintptr_t* already_seen;
    if (unlikely(bigstack_calloc_w(raw_variant_ctl, &orig_dups) ||
                 bigstack_calloc_w(raw_variant_ctl, &already_seen))) {
      goto RmDup_ret_NOMEM;
    }
    // variant_id_htable + htable_dup_base is structured as follows:
    // - Empty cells in variant_id_htable are represented by UINT32_MAX.
    // - Non-duplicate variant IDs are usually stored at
    //   variant_id_htable[hash(ID)].  Collisions are resolved with linear
    //   probing.
    // - Duplicate IDs are represented by a (2^31 + HALF the first
    //   htable_dup_base cell index) entry in variant_id_htable[hash(ID)].
    // - A htable_dup_base cell contains two parts:
    //   - First element (array index [k]) is a variant_uidx.
    //   - Second element (array index [k+1]) is another htable_dup_base cell
    //     index (not halved), or UINT32_MAX if we've reached the last ID.
    // Thus, we can initialize a bitarray with just the positions correspond to
    // duplicate-ID variant_uidxs set by looking at just the even entries of
    // htable_dup_base[].  (The bitarray may have some false positives from
    // recently-filtered variants and/or the missing-ID code.)
    for (uint32_t uii = 0; uii != orig_dup_ct; ++uii) {
      SetBit(htable_dup_base[uii * 2], orig_dups);
    }
    const uint32_t orig_variant_ct = *variant_ct_ptr;
    uint32_t dup_recheck_needed = 0;
    if (orig_variant_ct != raw_variant_ct) {
      BitvecAnd(variant_include, raw_variant_ctl, orig_dups);
      const uint32_t subsetted_dup_ct = PopcountWords(orig_dups, raw_variant_ctl);
      if (subsetted_dup_ct != orig_dup_ct) {
        if (!subsetted_dup_ct) {
          logputs("Note: Skipping --rm-dup since no duplicate IDs remain.\n");
          goto RmDup_ret_1;
        }
        orig_dup_ct = subsetted_dup_ct;
        dup_recheck_needed = 1;
      }
    }
    if (!missing_varid_match) {
      missing_varid_match = &(g_one_char_strs[92]);
    }
    const uint32_t missing_varid_blen = strlen(missing_varid_match) + 1;
    uint32_t* orig_dups_cumulative_popcounts = nullptr;
    // Load all relevant INFO lines before main loop, since main loop can
    // perform out-of-order lookups.
    const char** dup_info_strs = nullptr;
    if (pvar_info_reload && (rmdup_mode < kRmDupExcludeAll)) {
      if (unlikely(bigstack_alloc_u32(raw_variant_ctl, &orig_dups_cumulative_popcounts) ||
                   bigstack_alloc_kcp(orig_dup_ct, &dup_info_strs))) {
        goto RmDup_ret_NOMEM;
      }
      FillCumulativePopcounts(orig_dups, raw_variant_ctl, orig_dups_cumulative_popcounts);
      unsigned char* bigstack_mark2 = g_bigstack_base;
      const uint32_t decompress_thread_ct = ClipU32(max_thread_ct - 1, 1, 4);
      reterr = SizeAndInitTextStream(pvar_info_reload, bigstack_left() / 4, decompress_thread_ct, &pvar_txs);
      if (unlikely(reterr)) {
        goto RmDup_ret_TSTREAM_FAIL;
      }
      logputs("--rm-dup: Loading INFO field... ");
      fflush(stdout);
      unsigned char* tmp_alloc_end = g_bigstack_end;
      char* line_iter = nullptr;  // gcc 14 warning
      // if INFO column exists, #CHROM header line guaranteed
      do {
        reterr = TextNextLineLstrip(&pvar_txs, &line_iter);
        if (unlikely(reterr)) {
          goto RmDup_ret_TSTREAM_FAIL;
        }
      } while (!tokequal_k(line_iter, "#CHROM"));
      uint32_t info_col_idx = 1;
      {
        line_iter = &(line_iter[6]);
        for (; ; ++info_col_idx) {
          char* token_start = FirstNonTspace(line_iter);
          if (IsEolnKns(*token_start)) {
            reterr = kPglRetRewindFail;
            logerrprintfww(kErrprintfRewind, pvar_info_reload);
            goto RmDup_ret_1;
          }
          line_iter = CurTokenEnd(token_start);
          if (strequal_k(token_start, "INFO", line_iter - token_start)) {
            break;
          }
        }
      }
      unsigned char* tmp_alloc_base = g_bigstack_base;
      uint32_t dup_variant_idx = 0;
      for (uint32_t variant_uidx = 0; ; ++variant_uidx) {
        reterr = TextNextLineLstrip(&pvar_txs, &line_iter);
        if (unlikely(reterr)) {
          goto RmDup_ret_TSTREAM_FAIL;
        }
        if (IsSet(orig_dups, variant_uidx)) {
          if (!memequal(variant_ids[variant_uidx], missing_varid_match, missing_varid_blen)) {
            char* info_start = NextTokenMult(line_iter, info_col_idx);
            line_iter = CurTokenEnd(info_start);
            const uint32_t info_slen = line_iter - info_start;
            if (StoreStringAtEndK(tmp_alloc_base, info_start, info_slen, &tmp_alloc_end, &(dup_info_strs[dup_variant_idx]))) {
              goto RmDup_ret_NOMEM;
            }
          }
          ++dup_variant_idx;
          if (dup_variant_idx == orig_dup_ct) {
            break;
          }
        }
      }
      if (CleanupTextStream2(pvar_info_reload, &pvar_txs, &reterr)) {
        goto RmDup_ret_1;
      }
      logputs("done.\n");
      BigstackEndSet(tmp_alloc_end);
      BigstackReset(bigstack_mark2);
    }
    char* list_write_iter = g_textbuf;
    char* list_flush = &(list_write_iter[kMaxMediumLine]);
    char* mismatch_write_iter = nullptr;
    char* mismatch_flush = nullptr;
    if (rmdup_mode < kRmDupExcludeMismatch) {
      if (unlikely(bigstack_alloc_c(kMaxMediumLine + kMaxIdBlen, &mismatch_write_iter))) {
        goto RmDup_ret_NOMEM;
      }
      mismatch_flush = &(mismatch_write_iter[kMaxMediumLine]);
    }
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = orig_dups[0];
    uint32_t duplicate_ct = 0;
    uint32_t mismatch_ct = 0;
    uint32_t first_allele_ct = 2;
    uint32_t first_qual_is_present = 0;
    float first_qual = 0.0;
    uint32_t first_filter_is_present = 0;
    uint32_t first_filter_npass = 0;
    const char* first_filter_str = nullptr;
    const char* first_info_str = nullptr;
    double first_cm = 0.0;

    const uint32_t sample_ctb2 = NypCtToWordCt(sample_ct) * sizeof(intptr_t);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t sample_ctb = sample_ctl * sizeof(intptr_t);
    PgenVariant first_pgv;
    PreinitPgv(&first_pgv);
    PgenVariant cur_pgv;
    PreinitPgv(&cur_pgv);
    PgrSampleSubsetIndex pssi;
    if (simple_pgrp && (rmdup_mode < kRmDupExcludeAll)) {
      const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
      uint32_t* sample_include_cumulative_popcounts;
      if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                   BigstackAllocPgv(sample_ct, allele_idx_offsets != nullptr, PgrGetGflags(simple_pgrp), &first_pgv) ||
                   BigstackAllocPgv(sample_ct, allele_idx_offsets != nullptr, PgrGetGflags(simple_pgrp), &cur_pgv))) {
        goto RmDup_ret_NOMEM;
      }
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    } else {
      PgrClearSampleSubsetIndex(simple_pgrp, &pssi);
    }
    uint32_t missing_ct = 0;
    for (uint32_t variant_idx = 0; variant_idx != orig_dup_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(orig_dups, &variant_uidx_base, &cur_bits);
      if (IsSet(already_seen, variant_uidx)) {
        continue;
      }
      const char* cur_id = variant_ids[variant_uidx];
      if (memequal(cur_id, missing_varid_match, missing_varid_blen)) {
        ++missing_ct;
        continue;
      }
      uint32_t first_llidx;
      const uint32_t variant_uidx_ll_first = VariantIdDupHtableFind(cur_id, variant_ids, variant_id_htable, htable_dup_base, strlen(cur_id), variant_id_htable_size, max_variant_id_slen, &first_llidx);
      assert(first_llidx != UINT32_MAX);
      // 1. Verify this is still a duplicate in the current filtering state.
      // 2. If exclude-all or force-first mode, we're done; otherwise:
      //   3. Load variant information and genotype data for variant_uidx.
      //   4. Check for inequality.  If any inequality found, write variant ID
      //      to {output prefix}.rmdup.mismatch (lazy-opening the file before
      //      first write) if not in exclude-mismatch mode.
      if (dup_recheck_needed) {
        uint32_t is_still_dup = 0;
        uint32_t dupcheck_llidx = first_llidx;
        for (uint32_t dupcheck_vidx = variant_uidx_ll_first; ; dupcheck_llidx = htable_dup_base[dupcheck_llidx + 1]) {
          if ((dupcheck_vidx != variant_uidx) && IsSet(orig_dups, dupcheck_vidx)) {
            is_still_dup = 1;
            break;
          }
          if (dupcheck_llidx == UINT32_MAX) {
            break;
          }
          dupcheck_vidx = htable_dup_base[dupcheck_llidx];
        }
        if (!is_still_dup) {
          continue;
        }
      }
      if (save_list) {
        if (!list_file) {
          snprintf(outname_end, kMaxOutfnameExtBlen, ".rmdup.list");
          if (unlikely(fopen_checked(outname, FOPEN_WB, &list_file))) {
            goto RmDup_ret_OPEN_FAIL;
          }
        }
        list_write_iter = strcpya(list_write_iter, variant_ids[variant_uidx]);
        AppendBinaryEoln(&list_write_iter);
        if (unlikely(fwrite_ck(list_flush, list_file, &list_write_iter))) {
          goto RmDup_ret_WRITE_FAIL;
        }
      }
      ++duplicate_ct;
      if (rmdup_mode >= kRmDupExcludeAll) {
        uint32_t dupset_vidx = variant_uidx_ll_first;
        for (uint32_t cur_llidx = first_llidx; ; cur_llidx = htable_dup_base[cur_llidx + 1]) {
          SetBit(dupset_vidx, already_seen);
          ClearBit(dupset_vidx, variant_include);
          if (cur_llidx == UINT32_MAX) {
            break;
          }
          dupset_vidx = htable_dup_base[cur_llidx];
        }
        if (rmdup_mode == kRmDupForceFirst) {
          SetBit(variant_uidx, variant_include);
        }
        continue;
      }
      const uint32_t first_chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx);
      const uint32_t first_bp = variant_bps[variant_uidx];
      uintptr_t allele_idx_offset_base;
      if (!allele_idx_offsets) {
        allele_idx_offset_base = 2 * variant_uidx;
      } else {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        first_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* first_alleles = &(allele_storage[allele_idx_offset_base]);
      if (pvar_qual_present) {
        first_qual_is_present = IsSet(pvar_qual_present, variant_uidx);
        if (first_qual_is_present) {
          first_qual = pvar_quals[variant_uidx];
        }
      }
      if (pvar_filter_present) {
        first_filter_is_present = IsSet(pvar_filter_present, variant_uidx);
        if (first_filter_is_present) {
          first_filter_npass = IsSet(pvar_filter_npass, variant_uidx);
          if (first_filter_npass) {
            first_filter_str = pvar_filter_storage[variant_uidx];
          }
        }
      }
      if (dup_info_strs) {
        const uint32_t subsetted_idx = RawToSubsettedPos(orig_dups, orig_dups_cumulative_popcounts, variant_uidx);
        first_info_str = dup_info_strs[subsetted_idx];
      }
      if (variant_cms) {
        first_cm = variant_cms[variant_uidx];
      }

      uint32_t cur_llidx = first_llidx;
      uint32_t is_mismatch = 0;
      for (uint32_t ll_variant_uidx = variant_uidx_ll_first; ; ll_variant_uidx = htable_dup_base[cur_llidx], cur_llidx = htable_dup_base[cur_llidx + 1]) {
        if ((variant_uidx != ll_variant_uidx) && IsSet(orig_dups, ll_variant_uidx)) {
          SetBit(ll_variant_uidx, already_seen);
          ClearBit(ll_variant_uidx, variant_include);
          if (is_mismatch) {
            continue;
          }
          // Check .pvar fields for equality.
          is_mismatch = 1;
          if ((GetVariantChrFoIdx(cip, ll_variant_uidx) != first_chr_fo_idx) ||
              (variant_bps[ll_variant_uidx] != first_bp)) {
            continue;
          }
          if (!allele_idx_offsets) {
            allele_idx_offset_base = 2 * ll_variant_uidx;
          } else {
            allele_idx_offset_base = allele_idx_offsets[ll_variant_uidx];
            if (first_allele_ct != allele_idx_offsets[ll_variant_uidx + 1] - allele_idx_offset_base) {
              continue;
            }
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          uint32_t aidx = 0;
          for (; aidx != first_allele_ct; ++aidx) {
            if (!strequal_overread(first_alleles[aidx], cur_alleles[aidx])) {
              break;
            }
          }
          if (aidx != first_allele_ct) {
            continue;
          }
          if (pvar_qual_present) {
            if ((IsSet(pvar_qual_present, ll_variant_uidx) != first_qual_is_present) || (first_qual_is_present && (pvar_quals[ll_variant_uidx] != first_qual))) {
              continue;
            }
          }
          if (pvar_filter_present) {
            if (IsSet(pvar_filter_present, ll_variant_uidx) != first_filter_is_present) {
              continue;
            }
            if (first_filter_is_present) {
              if (IsSet(pvar_filter_npass, ll_variant_uidx) != first_filter_npass) {
                continue;
              }
              if (first_filter_npass && (!strequal_overread(first_filter_str, pvar_filter_storage[ll_variant_uidx]))) {
                continue;
              }
            }
          }
          if (first_info_str) {
            const uint32_t subsetted_idx = RawToSubsettedPos(orig_dups, orig_dups_cumulative_popcounts, ll_variant_uidx);
            if (!strequal_overread(first_info_str, dup_info_strs[subsetted_idx])) {
              continue;
            }
          }
          if (variant_cms) {
            if (variant_cms[ll_variant_uidx] != first_cm) {
              continue;
            }
          }
          is_mismatch = 0;
        }
        if (cur_llidx == UINT32_MAX) {
          break;
        }
      }
      if ((!is_mismatch) && first_pgv.genovec) {
        // Avoid loading genotypes when possible.
        reterr = PgrGetMDp(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &first_pgv);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, variant_uidx);
          goto RmDup_ret_1;
        }
        ZeroTrailingNyps(sample_ct, first_pgv.genovec);
        if (first_pgv.phasepresent_ct) {
          BitvecAnd(first_pgv.phasepresent, sample_ctl, first_pgv.phaseinfo);
        }
        cur_llidx = first_llidx;
        is_mismatch = 1;
        for (uint32_t ll_variant_uidx = variant_uidx_ll_first; ; ll_variant_uidx = htable_dup_base[cur_llidx], cur_llidx = htable_dup_base[cur_llidx + 1]) {
          if ((variant_uidx != ll_variant_uidx) && IsSet(orig_dups, ll_variant_uidx)) {
            reterr = PgrGetMDp(sample_include, pssi, sample_ct, ll_variant_uidx, simple_pgrp, &cur_pgv);
            if (unlikely(reterr)) {
              PgenErrPrintNV(reterr, ll_variant_uidx);
              goto RmDup_ret_1;
            }
            // todo: multidosage, multidphase
            if ((first_pgv.patch_01_ct != cur_pgv.patch_01_ct) ||
                (first_pgv.patch_10_ct != cur_pgv.patch_10_ct) ||
                (first_pgv.phasepresent_ct != cur_pgv.phasepresent_ct) ||
                (first_pgv.dosage_ct != cur_pgv.dosage_ct) ||
                (first_pgv.dphase_ct != cur_pgv.dphase_ct)) {
              break;
            }
            ZeroTrailingNyps(sample_ct, cur_pgv.genovec);
            if (!memequal(first_pgv.genovec, cur_pgv.genovec, sample_ctb2)) {
              break;
            }
            if (first_pgv.patch_01_ct) {
              if ((!memequal(first_pgv.patch_01_set, cur_pgv.patch_01_set, sample_ctb)) ||
                  (!memequal(first_pgv.patch_01_vals, cur_pgv.patch_01_vals, first_pgv.patch_01_ct * sizeof(AlleleCode)))) {
                break;
              }
            }
            if (first_pgv.patch_10_ct) {
              if ((!memequal(first_pgv.patch_10_set, cur_pgv.patch_10_set, sample_ctb)) ||
                  (!memequal(first_pgv.patch_10_vals, cur_pgv.patch_10_vals, first_pgv.patch_10_ct * sizeof(AlleleCode) * 2))) {
                break;
              }
            }
            if (first_pgv.phasepresent_ct) {
              BitvecAnd(cur_pgv.phasepresent, sample_ctl, cur_pgv.phaseinfo);
              if ((!memequal(first_pgv.phasepresent, cur_pgv.phasepresent, sample_ctb)) ||
                  (!memequal(first_pgv.phaseinfo, cur_pgv.phaseinfo, sample_ctb))) {
                break;
              }
            }
            if (first_pgv.dosage_ct) {
              if ((!memequal(first_pgv.dosage_present, cur_pgv.dosage_present, sample_ctb)) ||
                  (!memequal(first_pgv.dosage_main, cur_pgv.dosage_main, first_pgv.dosage_ct * sizeof(Dosage)))) {
                break;
              }
              if (first_pgv.dphase_ct) {
                if ((!memequal(first_pgv.dphase_present, cur_pgv.dphase_present, sample_ctb)) ||
                    (!memequal(first_pgv.dphase_delta, cur_pgv.dphase_delta, first_pgv.dphase_ct * sizeof(SDosage)))) {
                  break;
                }
              }
            }
          }
          if (cur_llidx == UINT32_MAX) {
            is_mismatch = 0;
            break;
          }
        }
      }
      if (is_mismatch) {
        ++mismatch_ct;
        if (rmdup_mode == kRmDupExcludeMismatch) {
          ClearBit(variant_uidx, variant_include);
        } else {
          if (rmdup_mode == kRmDupRetainMismatch) {
            cur_llidx = first_llidx;
            for (uint32_t ll_variant_uidx = variant_uidx_ll_first; ; ll_variant_uidx = htable_dup_base[cur_llidx], cur_llidx = htable_dup_base[cur_llidx + 1]) {
              if ((variant_uidx != ll_variant_uidx) && IsSet(orig_dups, ll_variant_uidx)) {
                SetBit(ll_variant_uidx, variant_include);
              }
              if (cur_llidx == UINT32_MAX) {
                break;
              }
            }
          }
          if (mismatch_file == nullptr) {
            snprintf(outname_end, kMaxOutfnameExtBlen, ".rmdup.mismatch");
            if (unlikely(fopen_checked(outname, FOPEN_WB, &mismatch_file))) {
              goto RmDup_ret_OPEN_FAIL;
            }
          }
          mismatch_write_iter = strcpya(mismatch_write_iter, variant_ids[variant_uidx]);
          AppendBinaryEoln(&mismatch_write_iter);
          if (unlikely(fwrite_ck(mismatch_flush, mismatch_file, &mismatch_write_iter))) {
            goto RmDup_ret_WRITE_FAIL;
          }
        }
      }
    }
    if (missing_ct) {
      logprintf("--rm-dup: %u missing-ID variant%s skipped.\n", missing_ct, (missing_ct == 1)? "" : "s");
    }
    if (mismatch_file != nullptr) {
      if (unlikely(fclose_flush_null(mismatch_flush, mismatch_write_iter, &mismatch_file))) {
        goto RmDup_ret_WRITE_FAIL;
      }
      logerrprintfww("%s: %u duplicate ID%s with inconsistent %svariant information detected by --rm-dup; see %s .\n", (rmdup_mode == kRmDupError)? "Error" : "Warning", mismatch_ct, (mismatch_ct == 1)? "" : "s", simple_pgrp? "genotype data or " : "", outname);
      if (rmdup_mode == kRmDupError) {
        reterr = kPglRetInconsistentInput;
        goto RmDup_ret_1;
      }
    } else if (mismatch_ct) {
      logprintfww("Note: %u duplicate ID%s with inconsistent %svariant information detected by --rm-dup exclude-mismatch; all copies removed.\n", mismatch_ct, (mismatch_ct == 1)? "" : "s", simple_pgrp? "genotype data or " : "");
    }
    *variant_ct_ptr = PopcountWords(variant_include, raw_variant_ctl);
    const uint32_t removed_variant_ct = orig_variant_ct - (*variant_ct_ptr);
    logprintfww("--rm-dup: %u duplicated ID%s, %u variant%s removed.\n", duplicate_ct, (duplicate_ct == 1)? "" : "s", removed_variant_ct, (removed_variant_ct == 1)? "" : "s");
    if (list_file != nullptr) {
      if (unlikely(fclose_flush_null(list_flush, list_write_iter, &list_file))) {
        goto RmDup_ret_WRITE_FAIL;
      }
      *outname_end = '\0';
      logprintfww("Full duplicate ID list written to %s.rmdup.list .\n", outname);
    }
  }
  while (0) {
  RmDup_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  RmDup_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  RmDup_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pvar_info_reload, &pvar_txs);
    break;
  RmDup_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 RmDup_ret_1:
  fclose_cond(mismatch_file);
  fclose_cond(list_file);
  CleanupTextStream2(pvar_info_reload, &pvar_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

void RandomThinProb(const char* flagname_p, const char* unitname, double thin_keep_prob, uint32_t raw_item_ct, sfmt_t* sfmtp, uintptr_t* item_include, uint32_t* item_ct_ptr) {
  // possible todo: try using truncated geometric distribution, like --dummy
  // can also parallelize this
  const uint32_t orig_item_ct = *item_ct_ptr;
  if (!orig_item_ct) {
    return;
  }
  const uint32_t uint32_thresh = S_CAST(uint32_t, thin_keep_prob * 4294967296.0 + 0.5);
  uintptr_t item_widx = 0;
  uintptr_t cur_bits = item_include[0];
  for (uint32_t item_idx = 0; item_idx != orig_item_ct; ++item_idx) {
    const uintptr_t lowbit = BitIter1y(item_include, &item_widx, &cur_bits);
    if (sfmt_genrand_uint32(sfmtp) >= uint32_thresh) {
      item_include[item_widx] ^= lowbit;
    }
  }
  const uint32_t new_item_ct = PopcountWords(item_include, BitCtToWordCt(raw_item_ct));
  *item_ct_ptr = new_item_ct;
  const uint32_t removed_ct = orig_item_ct - new_item_ct;
  logprintf("--%s: %u %s%s removed (%u remaining).\n", flagname_p, removed_ct, unitname, (removed_ct == 1)? "" : "s", new_item_ct);
  return;
}

PglErr RandomThinCt(const char* flagname_p, const char* unitname, uint32_t thin_keep_ct, uint32_t raw_item_ct, sfmt_t* sfmtp, uintptr_t* item_include, uint32_t* item_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t orig_item_ct = *item_ct_ptr;
    if (thin_keep_ct >= orig_item_ct) {
      logerrprintf("Warning: --%s parameter exceeds # of remaining %ss; skipping.\n", flagname_p, unitname);
      goto RandomThinCt_ret_1;
    }
    const uint32_t removed_ct = orig_item_ct - thin_keep_ct;
    const uint32_t raw_item_ctl = BitCtToWordCt(raw_item_ct);
    uintptr_t* perm_buf;
    uintptr_t* new_item_include;
    if (unlikely(bigstack_alloc_w(BitCtToWordCt(orig_item_ct), &perm_buf) ||
                 bigstack_alloc_w(raw_item_ctl, &new_item_include))) {
      goto RandomThinCt_ret_NOMEM;
    }
    // no actual interleaving here, but may as well use this function
    // note that this requires marker_ct >= 2
    GeneratePerm1Interleaved(orig_item_ct, thin_keep_ct, 0, 1, perm_buf, sfmtp);
    ExpandBytearr(perm_buf, item_include, raw_item_ctl, orig_item_ct, 0, new_item_include);
    memcpy(item_include, new_item_include, raw_item_ctl * sizeof(intptr_t));
    *item_ct_ptr = thin_keep_ct;
    logprintf("--%s: %u %s%s removed (%u remaining).\n", flagname_p, removed_ct, unitname, (removed_ct == 1)? "" : "s", thin_keep_ct);
  }
  while (0) {
  RandomThinCt_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 RandomThinCt_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

static const char kKeepRemoveFlagStrs[4][11] = {"keep", "remove", "keep-fam", "remove-fam"};

PglErr KeepOrRemove(const char* fnames, const SampleIdInfo* siip, uint32_t raw_sample_ct, KeepFlags flags, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto KeepOrRemove_ret_1;
    }
    const char* flag_name = kKeepRemoveFlagStrs[flags % 4];
    const LoadSampleIdsFlags load_sample_ids_flags = (flags & kfKeepFam)? (kfLoadSampleIdsMultifile | kfLoadSampleIdsFamOnly) : kfLoadSampleIdsMultifile;
    uintptr_t* seen_uidxs;
    uint32_t duplicate_ct;
    reterr = LoadSampleIds(fnames, sample_include, siip, flag_name, raw_sample_ct, orig_sample_ct, load_sample_ids_flags, &seen_uidxs, &duplicate_ct);
    if (unlikely(reterr)) {
      goto KeepOrRemove_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    if (flags & kfKeepRemove) {
      BitvecInvmask(seen_uidxs, raw_sample_ctl, sample_include);
    } else {
      memcpy(sample_include, seen_uidxs, raw_sample_ctl * sizeof(intptr_t));
    }
    const uint32_t sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    *sample_ct_ptr = sample_ct;
    logprintf("--%s: %u sample%s remaining.\n", flag_name, sample_ct, (sample_ct == 1)? "" : "s");
    if (duplicate_ct) {
      // "At least" since this does not count duplicate IDs absent from the
      // .fam.
      logerrprintf("Warning: At least %" PRIuPTR " duplicate ID%s in --%s file(s).\n", duplicate_ct, (duplicate_ct == 1)? "" : "s", flag_name);
    }
  }
 KeepOrRemove_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

static_assert(kTextbufSize >= 2 * kMaxIdBlen, "KeepOneId() needs to be updated.");
void KeepOneId(const char* sample_id_flattened, const SampleIdInfo* siip, uint32_t raw_sample_ct, uint32_t iid_sid, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  // Assumes sample_id_flattened is a multistr (sequence of nonempty
  // null-terminated strings, end-of-sequence is denoted by an empty string)
  // with 1-3 parts.
  // If it has 1 part, interpret as an IID, setting FID to 0.
  // If it has 2 parts, interpret as FID-IID unless iid_sid is true.
  const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
  {
    const char* cur_sid = nullptr;
    uintptr_t cur_sid_blen = 0;
    const char* cur_fid_iid;
    uintptr_t cur_fid_iid_blen;
    {
      const char* part1_start = sample_id_flattened;
      const uint32_t part1_blen = strlen(part1_start) + 1;
      const char* part2_start = &(part1_start[part1_blen]);
      char* fid_iid_write_iter = g_textbuf;
      if (!(*part2_start)) {
        fid_iid_write_iter = strcpya_k(fid_iid_write_iter, "0\t");
        fid_iid_write_iter = memcpya(fid_iid_write_iter, part1_start, part1_blen);
      } else {
        const uint32_t part2_blen = strlen(part2_start) + 1;
        const char* sid_start = &(part2_start[part2_blen]);
        if (!(*sid_start)) {
          if (!iid_sid) {
            fid_iid_write_iter = memcpyax(fid_iid_write_iter, part1_start, part1_blen - 1, '\t');
            fid_iid_write_iter = memcpya(fid_iid_write_iter, part2_start, part2_blen);

          } else {
            fid_iid_write_iter = strcpya_k(fid_iid_write_iter, "0\t");
            fid_iid_write_iter = memcpya(fid_iid_write_iter, part1_start, part1_blen);
            cur_sid = part2_start;
            cur_sid_blen = part2_blen;
          }
        } else {
          fid_iid_write_iter = memcpyax(fid_iid_write_iter, part1_start, part1_blen - 1, '\t');
          fid_iid_write_iter = memcpya(fid_iid_write_iter, part2_start, part2_blen);
          cur_sid = sid_start;
          cur_sid_blen = strlen(sid_start) + 1;
        }
      }
      cur_fid_iid_blen = fid_iid_write_iter - g_textbuf;
      cur_fid_iid = g_textbuf;
    }
    const char* sids = siip->sids;
    if (cur_sid) {
      if (!sids) {
        // No per-sample SID comparisons.  Either SID is always ok, or it never
        // is and we can immediately return sample_ct = 0.
        if ((siip->flags & kfSampleIdStrictSid0) && (!strequal_k_unsafe(cur_sid, "0"))) {
          goto KeepOneId_match_none;
        }
        cur_sid = nullptr;
        cur_sid_blen = 0;
      }
    } else {
      if (sids) {
        cur_sid = &(g_one_char_strs[96]);  // "0"
        cur_sid_blen = 2;
      }
    }
    const char* sample_ids = siip->sample_ids;
    const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
    const uintptr_t max_sid_blen = siip->max_sid_blen;
    if ((cur_fid_iid_blen > max_sample_id_blen) || (cur_sid_blen > max_sid_blen)) {
      goto KeepOneId_match_none;
    }
    uint32_t new_sample_ct = 0;
    for (uint32_t widx = 0; widx != raw_sample_ctl; ++widx) {
      uintptr_t cur_word = sample_include[widx];
      if (cur_word) {
        const uintptr_t sample_uidx_base = widx * kBitsPerWord;
        uintptr_t new_word = 0;
        do {
          const uint32_t sample_uidx_lowbits = ctzw(cur_word);
          const uintptr_t sample_uidx = sample_uidx_base + sample_uidx_lowbits;
          if (memequal(&(sample_ids[sample_uidx * max_sample_id_blen]), cur_fid_iid, cur_fid_iid_blen) &&
              ((!cur_sid) || memequal(&(sids[sample_uidx * max_sid_blen]), cur_sid, cur_sid_blen))) {
            new_word |= k1LU << sample_uidx_lowbits;
            ++new_sample_ct;
          }
          cur_word &= cur_word - 1;
        } while (cur_word);
        sample_include[widx] = new_word;
      }
    }
    logprintf("--indv: %u sample%s remaining.\n", new_sample_ct, (new_sample_ct == 1)? "" : "s");
    *sample_ct_ptr = new_sample_ct;
    return;
  }
 KeepOneId_match_none:
  ZeroWArr(raw_sample_ctl, sample_include);
  *sample_ct_ptr = 0;
  logputs("--indv: 0 samples remaining.\n");
}

// Minor extension of PLINK 1.x --filter.  (Renamed since --filter is not
// sufficiently self-describing; PLINK has lots of other filters on both
// samples and variants.  --filter is automatically converted to
// --keep-col-match for backward compatibility, though.)
PglErr KeepColMatch(const char* fname, const SampleIdInfo* siip, const char* strs_flattened, const char* col_name, uint32_t raw_sample_ct, uint32_t col_num, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto KeepColMatch_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* seen_xid_idxs;
    uintptr_t* keep_uidxs;
    if (unlikely(bigstack_calloc_w(BitCtToWordCt(orig_sample_ct), &seen_xid_idxs) ||
                 bigstack_calloc_w(raw_sample_ctl, &keep_uidxs))) {
      goto KeepColMatch_ret_NOMEM;
    }

    char* sorted_strbox;
    uintptr_t max_str_blen;
    uint32_t str_ct;
    if (unlikely(MultistrToStrboxDedupAlloc(strs_flattened, &sorted_strbox, &str_ct, &max_str_blen))) {
      goto KeepColMatch_ret_NOMEM;
    }

    reterr = SizeAndInitTextStream(fname, bigstack_left() - (bigstack_left() / 4), 1, &txs);
    if (unlikely(reterr)) {
      goto KeepColMatch_ret_TSTREAM_FAIL;
    }
    char* line_start;
    XidMode xid_mode;
    reterr = LoadXidHeader("keep-col-match", (siip->sids || (siip->flags & kfSampleIdStrictSid0))? kfXidHeaderFixedWidth : kfXidHeaderFixedWidthIgnoreSid, &line_idx, &txs, &xid_mode, &line_start);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        logerrputs("Error: Empty --keep-col-match file.\n");
        goto KeepColMatch_ret_MALFORMED_INPUT;
      }
      goto KeepColMatch_ret_TSTREAM_XID_FAIL;
    }
    const uint32_t id_col_ct = GetXidColCt(xid_mode);
    uint32_t postid_col_idx = 0;
    if (!col_name) {
      if (!col_num) {
        if (unlikely(id_col_ct == 3)) {
          logerrputs("Error: You must specify a --keep-col-match column with --keep-col-match-name or\n--keep-col-match-num.\n");
          goto KeepColMatch_ret_INCONSISTENT_INPUT;
        }
        col_num = 3;
      }
      if (unlikely(id_col_ct >= col_num)) {
        logerrputs("Error: --keep-col-match-num parameter too small (it refers to a sample ID\ncolumn in the --keep-col-match file).\n");
        goto KeepColMatch_ret_INCONSISTENT_INPUT;
      }
      postid_col_idx = col_num - id_col_ct;
    } else {
      if (unlikely(*line_start != '#')) {
        logerrputs("Error: --keep-col-match-name requires the --keep-col-match file to have a\nheader line starting with #FID or #IID.\n");
        goto KeepColMatch_ret_INCONSISTENT_INPUT;
      }
      const char* linebuf_iter = NextTokenMult(line_start, id_col_ct);
      if (unlikely(!linebuf_iter)) {
        logerrputs("Error: --keep-col-match-name column not found in --keep-col-match file.\n");
        goto KeepColMatch_ret_INCONSISTENT_INPUT;
      }
      const uint32_t col_name_slen = strlen(col_name);
      uint32_t cur_col_idx = 0;
      do {
        ++cur_col_idx;
        const char* token_end = CurTokenEnd(linebuf_iter);
        if ((S_CAST(uintptr_t, token_end - linebuf_iter) == col_name_slen) && memequal(linebuf_iter, col_name, col_name_slen)) {
          if (unlikely(postid_col_idx)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Multiple columns in --keep-col-match file are named '%s'.\n", col_name);
            goto KeepColMatch_ret_INCONSISTENT_INPUT_WW;
          }
          postid_col_idx = cur_col_idx;
        }
        linebuf_iter = FirstNonTspace(token_end);
      } while (!IsEolnKns(*linebuf_iter));
      if (unlikely(!postid_col_idx)) {
        logerrputs("Error: --keep-col-match-name column not found in --keep-col-match file.\n");
        goto KeepColMatch_ret_INCONSISTENT_INPUT;
      }
    }
    uint32_t* xid_map = nullptr;
    char* sorted_xidbox = nullptr;
    uintptr_t max_xid_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, orig_sample_ct, xid_mode, 0, &sorted_xidbox, &xid_map, &max_xid_blen);
    if (unlikely(reterr)) {
      goto KeepColMatch_ret_1;
    }
    char* idbuf = nullptr;
    if (unlikely(bigstack_alloc_c(max_xid_blen, &idbuf))) {
      goto KeepColMatch_ret_NOMEM;
    }
    if (*line_start == '#') {
      ++line_idx;
      line_start = TextGet(&txs);
    }
    for (; line_start; ++line_idx, line_start = TextGet(&txs)) {
      const char* linebuf_iter = line_start;
      uint32_t xid_idx_start;
      uint32_t xid_idx_end;
      if (SortedXidboxReadMultifind(sorted_xidbox, max_xid_blen, orig_sample_ct, 0, xid_mode, &linebuf_iter, &xid_idx_start, &xid_idx_end, idbuf)) {
        if (unlikely(!linebuf_iter)) {
          goto KeepColMatch_ret_MISSING_TOKENS;
        }
        continue;
      }
      if (unlikely(IsSet(seen_xid_idxs, xid_idx_start))) {
        logerrprintfww("Error: Sample ID on line %" PRIuPTR " of --keep-col-match file duplicates one earlier in the file.\n", line_idx);
        goto KeepColMatch_ret_MALFORMED_INPUT;
      }
      SetBit(xid_idx_start, seen_xid_idxs);
      linebuf_iter = NextTokenMult(linebuf_iter, postid_col_idx);
      if (unlikely(!linebuf_iter)) {
        goto KeepColMatch_ret_MISSING_TOKENS;
      }
      const char* token_end = CurTokenEnd(linebuf_iter);
      const int32_t ii = bsearch_strbox(linebuf_iter, sorted_strbox, token_end - linebuf_iter, max_str_blen, str_ct);
      if (ii != -1) {
        for (; xid_idx_start != xid_idx_end; ++xid_idx_start) {
          const uint32_t sample_uidx = xid_map[xid_idx_start];
          SetBit(sample_uidx, keep_uidxs);
        }
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto KeepColMatch_ret_TSTREAM_FAIL;
    }
    memcpy(sample_include, keep_uidxs, raw_sample_ctl * sizeof(intptr_t));
    const uint32_t sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    *sample_ct_ptr = sample_ct;
    logprintf("--keep-col-match: %u sample%s remaining.\n", sample_ct, (sample_ct == 1)? "" : "s");
  }
  while (0) {
  KeepColMatch_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  KeepColMatch_ret_TSTREAM_XID_FAIL:
    if (!TextStreamErrcode(&txs)) {
      break;
    }
  KeepColMatch_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--keep-col-match file", &txs);
    break;
  KeepColMatch_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  KeepColMatch_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --keep-col-match file has fewer tokens than expected.\n", line_idx);
  KeepColMatch_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  KeepColMatch_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 KeepColMatch_ret_1:
  BigstackReset(bigstack_mark);
  CleanupTextStream2("--keep-col-match file", &txs, &reterr);
  return reterr;
}

PglErr RequirePheno(const PhenoCol* pheno_cols, const char* pheno_names, const char* require_pheno_flattened, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t is_covar, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto RequirePheno_ret_1;
    }
    uint32_t required_pheno_ct = 0;
    uintptr_t max_required_pheno_blen = 2;
    uintptr_t* matched_phenos = nullptr;
    char* sorted_required_pheno_names = nullptr;
    if (require_pheno_flattened) {
      if (unlikely(MultistrToStrboxDedupAlloc(require_pheno_flattened, &sorted_required_pheno_names, &required_pheno_ct, &max_required_pheno_blen))) {
        goto RequirePheno_ret_NOMEM;
      }
      if (unlikely(bigstack_calloc_w(1 + (required_pheno_ct / kBitsPerWord), &matched_phenos))) {
        goto RequirePheno_ret_NOMEM;
      }
    } else {
      if (!pheno_ct) {
        logerrputs(is_covar? "Warning: No covariates loaded; ignoring --require-covar.\n" : "Warning: No phenotypes loaded; ignoring --require-pheno.\n");
        goto RequirePheno_ret_1;
      }
      required_pheno_ct = pheno_ct;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      if (sorted_required_pheno_names) {
        const char* cur_pheno_name = &(pheno_names[pheno_idx * max_pheno_name_blen]);
        const int32_t ii = bsearch_strbox(cur_pheno_name, sorted_required_pheno_names, strlen(cur_pheno_name), max_required_pheno_blen, required_pheno_ct);
        if (ii == -1) {
          continue;
        }
        SetBitI(ii, matched_phenos);
      }
      BitvecAnd(pheno_cols[pheno_idx].nonmiss, raw_sample_ctl, sample_include);
    }
    if (matched_phenos) {
      const uint32_t first_unmatched_idx = AdvTo0Bit(matched_phenos, 0);
      if (unlikely(first_unmatched_idx < required_pheno_ct)) {
        logerrprintfww("Error: --require-%s '%s' not loaded.\n", is_covar? "covar covariate" : "pheno phenotype", &(sorted_required_pheno_names[first_unmatched_idx * max_required_pheno_blen]));
        goto RequirePheno_ret_INCONSISTENT_INPUT;
      }
    }
    const uint32_t new_sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    const uint32_t removed_sample_ct = orig_sample_ct - new_sample_ct;
    logprintf("--require-%s: %u sample%s removed.\n", is_covar? "covar" : "pheno", removed_sample_ct, (removed_sample_ct == 1)? "" : "s");
    *sample_ct_ptr = new_sample_ct;
  }
  while (0) {
  RequirePheno_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  RequirePheno_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 RequirePheno_ret_1:
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

typedef struct KeepIfInternalTStruct {
  const PhenoCol* pheno_cols;
  const PhenoCol* covar_cols;
  const char* pheno_names;
  const char* covar_names;
  const char* errmsg_start;  // "Error: --keep-if " or "Error: --remove-if "
  uintptr_t max_pheno_name_blen;
  uintptr_t max_covar_name_blen;
  uint32_t raw_sample_ctl;
  uint32_t pheno_ct;
  uint32_t covar_ct;
  uint32_t affection_01;
} KeepIfInternalT;

const PhenoCol* GetPhenoCovarCol(const char* cur_name, const KeepIfInternalT* ctx) {
  const uint32_t name_blen = strlen(cur_name) + 1;
  const uintptr_t max_pheno_name_blen = ctx->max_pheno_name_blen;
  if (name_blen <= max_pheno_name_blen) {
    const char* pheno_names = ctx->pheno_names;
    const uint32_t pheno_ct = ctx->pheno_ct;
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      if (memequal(cur_name, &(pheno_names[pheno_idx * max_pheno_name_blen]), name_blen)) {
        return &(ctx->pheno_cols[pheno_idx]);
      }
    }
  }
  const uintptr_t max_covar_name_blen = ctx->max_covar_name_blen;
  if (name_blen <= max_covar_name_blen) {
    const char* covar_names = ctx->covar_names;
    const uint32_t covar_ct = ctx->pheno_ct;
    for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
      if (memequal(cur_name, &(covar_names[covar_idx * max_covar_name_blen]), name_blen)) {
        return &(ctx->covar_cols[covar_idx]);
      }
    }
  }
  return nullptr;
}

PglErr KeepIfInternal(const CmpExpr* cmp_expr, const KeepIfInternalT* ctx, uintptr_t* cur_include) {
  unsigned char* bigstack_mark = g_bigstack_base;
  const char* key = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = ctx->raw_sample_ctl;
    const CmpExprType etype = cmp_expr->etype;
    if (etype == kCmpExprTypeExists) {
      key = cmp_expr->args.k.key;
      const PhenoCol* cur_pheno_col = GetPhenoCovarCol(key, ctx);
      if (unlikely(!cur_pheno_col)) {
        goto KeepIfInternal_ret_PHENO_NAME_NOT_FOUND;
      }
      const uintptr_t* pheno_nm = cur_pheno_col->nonmiss;
      BitvecAnd(pheno_nm, raw_sample_ctl, cur_include);
    } else if (etype <= kCmpExprTypeStrEq) {
      uint32_t is_neq;
      if (etype <= kCmpExprTypeGeq) {
        key = cmp_expr->args.kn.key;
        is_neq = (etype == kCmpExprTypeNoteq);
      } else {
        key = cmp_expr->args.ks.key;
        is_neq = (etype == kCmpExprTypeStrNoteq);
      }
      const PhenoCol* cur_pheno_col = GetPhenoCovarCol(key, ctx);
      if (unlikely(!cur_pheno_col)) {
        goto KeepIfInternal_ret_PHENO_NAME_NOT_FOUND;
      }
      const uintptr_t* pheno_nm = cur_pheno_col->nonmiss;
      if (!is_neq) {
        BitvecAnd(pheno_nm, raw_sample_ctl, cur_include);
      }
      uintptr_t* cur_include_intersect;
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &cur_include_intersect))) {
        goto KeepIfInternal_ret_NOMEM;
      }
      memcpy(cur_include_intersect, cur_include, raw_sample_ctl * sizeof(intptr_t));
      if (is_neq) {
        BitvecAnd(pheno_nm, raw_sample_ctl, cur_include_intersect);
      }
      const uint32_t intersect_ct = PopcountWords(cur_include_intersect, raw_sample_ctl);
      if (cur_pheno_col->type_code == kPhenoDtypeQt) {
        if (unlikely(etype > kCmpExprTypeGeq)) {
          snprintf(g_logbuf, kLogbufSize, "%squantitative phenotype/covariate '%s' must be compared to a number, not '%s'.", ctx->errmsg_start, key, cmp_expr->args.ks.str_value);
          goto KeepIfInternal_ret_INCONSISTENT_INPUT_WW;
        }
        const double* pheno_vals = cur_pheno_col->data.qt;
        const double val = cmp_expr->args.kn.value;
        uintptr_t sample_uidx_base = 0;
        uintptr_t cur_bits = cur_include_intersect[0];
        if (etype == kCmpExprTypeLe) {
          for (uint32_t sample_idx = 0; sample_idx != intersect_ct; ++sample_idx) {
            const uintptr_t sample_uidx = BitIter1(cur_include_intersect, &sample_uidx_base, &cur_bits);
            if (pheno_vals[sample_uidx] >= val) {
              ClearBit(sample_uidx, cur_include);
            }
          }
        } else if (etype == kCmpExprTypeLeq) {
          for (uint32_t sample_idx = 0; sample_idx != intersect_ct; ++sample_idx) {
            const uintptr_t sample_uidx = BitIter1(cur_include_intersect, &sample_uidx_base, &cur_bits);
            if (pheno_vals[sample_uidx] > val) {
              ClearBit(sample_uidx, cur_include);
            }
          }
        } else if (etype == kCmpExprTypeNoteq) {
          for (uint32_t sample_idx = 0; sample_idx != intersect_ct; ++sample_idx) {
            const uintptr_t sample_uidx = BitIter1(cur_include_intersect, &sample_uidx_base, &cur_bits);
            if (pheno_vals[sample_uidx] == val) {
              ClearBit(sample_uidx, cur_include);
            }
          }
        } else if (etype == kCmpExprTypeEq) {
          for (uint32_t sample_idx = 0; sample_idx != intersect_ct; ++sample_idx) {
            const uintptr_t sample_uidx = BitIter1(cur_include_intersect, &sample_uidx_base, &cur_bits);
            if (pheno_vals[sample_uidx] != val) {
              ClearBit(sample_uidx, cur_include);
            }
          }
        } else if (etype == kCmpExprTypeGe) {
          for (uint32_t sample_idx = 0; sample_idx != intersect_ct; ++sample_idx) {
            const uintptr_t sample_uidx = BitIter1(cur_include_intersect, &sample_uidx_base, &cur_bits);
            if (pheno_vals[sample_uidx] <= val) {
              ClearBit(sample_uidx, cur_include);
            }
          }
        } else {
          assert(etype == kCmpExprTypeGeq);
          for (uint32_t sample_idx = 0; sample_idx != intersect_ct; ++sample_idx) {
            const uintptr_t sample_uidx = BitIter1(cur_include_intersect, &sample_uidx_base, &cur_bits);
            if (pheno_vals[sample_uidx] < val) {
              ClearBit(sample_uidx, cur_include);
            }
          }
        }
      } else if (cur_pheno_col->type_code == kPhenoDtypeCc) {
        const uint32_t affection_01 = ctx->affection_01;
        uint32_t val_12 = 0;  // 1 = control, 2 = case
        if (etype > kCmpExprTypeGeq) {
          const char* str_value = cmp_expr->args.ks.str_value;
          const uint32_t val_slen = strlen(str_value);
          if (val_slen == 4) {
            if (MatchUpperK(str_value, "CASE")) {
              val_12 = 2;
            } else if (MatchUpperK(str_value, "CTRL")) {
              val_12 = 1;
            }
          } else if (MatchUpperKLen(str_value, "CONTROL", val_slen)) {
            val_12 = 1;
          }
        } else {
          if (unlikely((!is_neq) && (etype != kCmpExprTypeEq))) {
            snprintf(g_logbuf, kLogbufSize, "%sbinary phenotype/covariate '%s' must be compared with == or !=; <, <=, >=, and > are not allowed.\n", ctx->errmsg_start, key);
            goto KeepIfInternal_ret_INCONSISTENT_INPUT_WW;
          }
          const double val = cmp_expr->args.kn.value;
          if (val == u31tod(1 - affection_01)) {
            val_12 = 1;
          } else if (val == u31tod(2 - affection_01)) {
            val_12 = 2;
          }
        }
        if (unlikely(!val_12)) {
          snprintf(g_logbuf, kLogbufSize, "%sbinary phenotype/covariate must be compared to case/%c or control/ctrl/%c.\n", ctx->errmsg_start, '2' - affection_01, '1' - affection_01);
          goto KeepIfInternal_ret_INCONSISTENT_INPUT_WW;
        }
        if (val_12 == 2) {
          BitvecAnd(cur_pheno_col->data.cc, raw_sample_ctl, cur_include);
        } else {
          BitvecInvmask(cur_pheno_col->data.cc, raw_sample_ctl, cur_include);
        }
      } else {
        assert(cur_pheno_col->type_code == kPhenoDtypeCat);
        if (unlikely(etype <= kCmpExprTypeGeq)) {
          snprintf(g_logbuf, kLogbufSize, "%scategorical phenotype/covariate '%s' must be compared to a category name, not a number.", ctx->errmsg_start, key);
          goto KeepIfInternal_ret_INCONSISTENT_INPUT_WW;
        }
        const char* str_value = cmp_expr->args.ks.str_value;
        const uint32_t nonnull_cat_ct = cur_pheno_col->nonnull_category_ct;
        const uint32_t cat_idx = 1 + bsearch_strptr_natural(str_value, &(cur_pheno_col->category_names[1]), nonnull_cat_ct);
        if (!cat_idx) {
          logerrprintfww("Warning%scategorical phenotype/covariate '%s' does not have a category named '%s'.\n", &(ctx->errmsg_start[5]), key, str_value);
          if (!is_neq) {
            ZeroWArr(raw_sample_ctl, cur_include);
          }
        } else {
          const uint32_t* cur_cats = cur_pheno_col->data.cat;
          uintptr_t sample_uidx_base = 0;
          uintptr_t cur_bits = cur_include_intersect[0];
          if (!is_neq) {
            for (uint32_t sample_idx = 0; sample_idx != intersect_ct; ++sample_idx) {
              const uintptr_t sample_uidx = BitIter1(cur_include_intersect, &sample_uidx_base, &cur_bits);
              if (cur_cats[sample_uidx] != cat_idx) {
                ClearBit(sample_uidx, cur_include);
              }
            }
          } else {
            for (uint32_t sample_idx = 0; sample_idx != intersect_ct; ++sample_idx) {
              const uintptr_t sample_uidx = BitIter1(cur_include_intersect, &sample_uidx_base, &cur_bits);
              if (cur_cats[sample_uidx] == cat_idx) {
                ClearBit(sample_uidx, cur_include);
              }
            }
          }
        }
      }
    } else {
      uintptr_t* other_include;
      if (unlikely(bigstack_alloc_w(raw_sample_ctl, &other_include))) {
        goto KeepIfInternal_ret_NOMEM;
      }
      memcpy(other_include, cur_include, raw_sample_ctl * sizeof(intptr_t));
      reterr = KeepIfInternal(cmp_expr->args.jct.children[0], ctx, other_include);
      if (unlikely(reterr)) {
        goto KeepIfInternal_ret_1;
      }
      if (etype == kCmpExprTypeNot) {
        BitvecInvmask(other_include, raw_sample_ctl, cur_include);
      } else {
        reterr = KeepIfInternal(cmp_expr->args.jct.children[1], ctx, cur_include);
        if (unlikely(reterr)) {
          goto KeepIfInternal_ret_1;
        }
        if (etype == kCmpExprTypeAnd) {
          BitvecAnd(other_include, raw_sample_ctl, cur_include);
        } else {
          assert(etype == kCmpExprTypeOr);
          BitvecOr(other_include, raw_sample_ctl, cur_include);
        }
      }
    }
  }
  while (0) {
  KeepIfInternal_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  KeepIfInternal_ret_PHENO_NAME_NOT_FOUND:
    snprintf(g_logbuf, kLogbufSize, "%sphenotype/covariate '%s' not loaded.\n", ctx->errmsg_start, key);
  KeepIfInternal_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 KeepIfInternal_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr KeepRemoveIf(const CmpExpr* cmp_expr, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t affection_01, uint32_t is_remove, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  KeepIfInternalT ctx;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto KeepRemoveIf_ret_1;
    }
    ctx.pheno_cols = pheno_cols;
    ctx.covar_cols = covar_cols;
    ctx.pheno_names = pheno_names;
    ctx.covar_names = covar_names;
    ctx.errmsg_start = is_remove? "Error: --remove-if " : "Error: --keep-if ";
    ctx.max_pheno_name_blen = max_pheno_name_blen;
    ctx.max_covar_name_blen = max_covar_name_blen;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    ctx.raw_sample_ctl = raw_sample_ctl;
    ctx.pheno_ct = pheno_ct;
    ctx.covar_ct = covar_ct;
    ctx.affection_01 = affection_01;
    CmpExpr cmp_expr_slot;
    if (is_remove) {
      cmp_expr_slot.etype = kCmpExprTypeNot;
      // not thrilled by this cast, but the alternatives all have comparable
      // drawbacks
      cmp_expr_slot.args.jct.children[0] = K_CAST(CmpExpr*, cmp_expr);
      cmp_expr = &cmp_expr_slot;
    }
    reterr = KeepIfInternal(cmp_expr, &ctx, sample_include);
    if (unlikely(reterr)) {
      goto KeepRemoveIf_ret_1;
    }
    const uint32_t new_sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    const uint32_t removed_sample_ct = orig_sample_ct - new_sample_ct;
    logprintf("--%s-if: %u sample%s removed.\n", is_remove? "remove" : "keep", removed_sample_ct, (removed_sample_ct == 1)? "" : "s");
    *sample_ct_ptr = new_sample_ct;
  }
 KeepRemoveIf_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr KeepRemoveCatsInternal(const PhenoCol* cur_pheno_col, const char* cats_fname, const char* cat_names_flattened, uint32_t raw_sample_ct, uint32_t is_remove, uint32_t max_thread_ct, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  char file_descrip[32];
  if (is_remove) {
    strcpy_k(file_descrip, "--remove-cats file");
  } else {
    strcpy_k(file_descrip, "--keep-cats file");
  }
  TokenStream tks;
  PreinitTokenStream(&tks);
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto KeepRemoveCatsInternal_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t cat_ct = cur_pheno_col->nonnull_category_ct + 1;
    const uint32_t cat_ctl = BitCtToWordCt(cat_ct);
    uintptr_t* affected_samples;
    uintptr_t* cat_include;
    if (unlikely(bigstack_calloc_w(raw_sample_ctl, &affected_samples) ||
                 bigstack_alloc_w(cat_ctl, &cat_include))) {
      goto KeepRemoveCatsInternal_ret_NOMEM;
    }
    SetAllBits(cat_ct, cat_include);
    const char* const* category_names = cur_pheno_col->category_names;
    uint32_t* cat_id_htable;
    uint32_t id_htable_size;
    reterr = AllocAndPopulateIdHtableMt(cat_include, category_names, cat_ct, 0, max_thread_ct, &cat_id_htable, nullptr, &id_htable_size, nullptr);
    if (unlikely(reterr)) {
      goto KeepRemoveCatsInternal_ret_1;
    }
    ZeroWArr(cat_ctl, cat_include);
    if (cats_fname) {
      reterr = InitTokenStream(cats_fname, MAXV(max_thread_ct - 1, 1), &tks);
      if (unlikely(reterr)) {
        goto KeepRemoveCatsInternal_ret_TKSTREAM_FAIL;
      }
      uintptr_t skip_ct = 0;
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
          *token_end = '\0';
          // can't overread, category_names not in main workspace
          const uint32_t cur_cat_idx = IdHtableFind(shard_iter, category_names, cat_id_htable, token_end - shard_iter, id_htable_size);
          if (cur_cat_idx == UINT32_MAX) {
            ++skip_ct;
          } else {
            SetBit(cur_cat_idx, cat_include);
          }
          shard_iter = token_end;
        }
      }
      if (unlikely(reterr != kPglRetEof)) {
        goto KeepRemoveCatsInternal_ret_TKSTREAM_FAIL;
      }
      if (CleanupTokenStream3(file_descrip, &tks, &reterr)) {
        goto KeepRemoveCatsInternal_ret_1;
      }
      if (skip_ct) {
        logerrprintf("Warning: %" PRIuPTR " --%s-cats categor%s not present.\n", skip_ct, is_remove? "remove" : "keep", (skip_ct == 1)? "y" : "ies");
      }
    }
    if (cat_names_flattened) {
      uint32_t skip_ct = 0;
      const char* cat_names_iter = cat_names_flattened;
      do {
        const uint32_t cat_name_slen = strlen(cat_names_iter);
        const uint32_t cur_cat_idx = IdHtableFind(cat_names_iter, category_names, cat_id_htable, cat_name_slen, id_htable_size);
        if (cur_cat_idx == UINT32_MAX) {
          ++skip_ct;
        } else {
          SetBit(cur_cat_idx, cat_include);
        }
        cat_names_iter = &(cat_names_iter[cat_name_slen + 1]);
      } while (*cat_names_iter);
      if (skip_ct) {
        logerrprintf("Warning: %u --%s-cat-names categor%s not present.\n", skip_ct, is_remove? "remove" : "keep", (skip_ct == 1)? "y" : "ies");
      }
    }
    const uint32_t selected_cat_ct = PopcountWords(cat_include, cat_ctl);
    if (!selected_cat_ct) {
      logerrprintf("Warning: No matching --%s-cat-names category names.\n", is_remove? "remove-cats/--remove" : "keep-cats/--keep");
    } else {
      const uint32_t* cur_cats = cur_pheno_col->data.cat;
      uintptr_t sample_uidx_base = 0;
      uintptr_t cur_bits = sample_include[0];
      for (uint32_t sample_idx = 0; sample_idx != orig_sample_ct; ++sample_idx) {
        const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
        const uint32_t cur_cat_idx = cur_cats[sample_uidx];
        if (IsSet(cat_include, cur_cat_idx)) {
          SetBit(sample_uidx, affected_samples);
        }
      }
      if (is_remove) {
        BitvecInvmask(affected_samples, raw_sample_ctl, sample_include);
      } else {
        BitvecAnd(affected_samples, raw_sample_ctl, sample_include);
      }
      const uint32_t new_sample_ct = PopcountWords(sample_include, raw_sample_ctl);
      const uint32_t removed_sample_ct = orig_sample_ct - new_sample_ct;
      logprintfww("--%s-cat-names: %u categor%s selected, %u sample%s removed.\n", is_remove? "remove-cats/--remove" : "keep-cats/--keep", selected_cat_ct, (selected_cat_ct == 1)? "y" : "ies", removed_sample_ct, (removed_sample_ct == 1)? "" : "s");
      *sample_ct_ptr = new_sample_ct;
    }
  }
  while (0) {
  KeepRemoveCatsInternal_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  KeepRemoveCatsInternal_ret_TKSTREAM_FAIL:
    TokenStreamErrPrint(file_descrip, &tks);
    break;
  }
 KeepRemoveCatsInternal_ret_1:
  CleanupTokenStream2(file_descrip, &tks, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr KeepRemoveCats(const char* cats_fname, const char* cat_names_flattened, const char* cat_phenoname, const PhenoCol* pheno_cols, const char* pheno_names, const PhenoCol* covar_cols, const char* covar_names, uint32_t raw_sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t covar_ct, uintptr_t max_covar_name_blen, uint32_t is_remove, uint32_t max_thread_ct, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  PglErr reterr = kPglRetSuccess;
  {
    if (!(*sample_ct_ptr)) {
      goto KeepRemoveCats_ret_1;
    }
    if (!cat_phenoname) {
      // Default behavior:
      // 1. If at least one categorical phenotype exists, fail on >= 2, select
      //    it if one.
      // 2. Otherwise, fail if 0 or >= 2 categorical covariates, select the
      //    categorical covariate if there's exactly one.
      uint32_t cat_pheno_idx = UINT32_MAX;
      const PhenoCol* cur_pheno_col = nullptr;
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        if (pheno_cols[pheno_idx].type_code == kPhenoDtypeCat) {
          if (unlikely(cat_pheno_idx != UINT32_MAX)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Multiple categorical phenotypes present. Use --%s-cat-pheno to specify which phenotype/covariate you want to filter on.\n", is_remove? "remove" : "keep");
            goto KeepRemoveCats_ret_INCONSISTENT_INPUT_WW;
          }
          cat_pheno_idx = pheno_idx;
        }
      }
      if (cat_pheno_idx != UINT32_MAX) {
        cur_pheno_col = &(pheno_cols[cat_pheno_idx]);
      } else {
        for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
          if (covar_cols[covar_idx].type_code == kPhenoDtypeCat) {
            if (unlikely(cat_pheno_idx != UINT32_MAX)) {
              snprintf(g_logbuf, kLogbufSize, "Error: Multiple categorical covariates and no categorical phenotype present. Use --%s-cat-pheno to specify which phenotype/covariate you want to filter on.\n", is_remove? "remove" : "keep");
              goto KeepRemoveCats_ret_INCONSISTENT_INPUT_WW;
            }
            cat_pheno_idx = covar_idx;
          }
        }
        if (unlikely(cat_pheno_idx == UINT32_MAX)) {
          snprintf(g_logbuf, kLogbufSize, "Error: --%s-cat-names requires a categorical phenotype or covariate.\n", is_remove? "remove-cats/--remove" : "keep-cats/--keep");
          goto KeepRemoveCats_ret_INCONSISTENT_INPUT_WW;
        }
        cur_pheno_col = &(covar_cols[cat_pheno_idx]);
      }
      reterr = KeepRemoveCatsInternal(cur_pheno_col, cats_fname, cat_names_flattened, raw_sample_ct, is_remove, max_thread_ct, sample_include, sample_ct_ptr);
      if (unlikely(reterr)) {
        goto KeepRemoveCats_ret_1;
      }
    } else {
      const uintptr_t name_blen = 1 + strlen(cat_phenoname);
      uint32_t success = 0;
      if (name_blen <= max_pheno_name_blen) {
        for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
          if (memequal(cat_phenoname, &(pheno_names[pheno_idx * max_pheno_name_blen]), name_blen)) {
            const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
            if (unlikely(cur_pheno_col->type_code != kPhenoDtypeCat)) {
              snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a categorical phenotype.\n", cat_phenoname);
              goto KeepRemoveCats_ret_INCONSISTENT_INPUT_WW;
            }
            reterr = KeepRemoveCatsInternal(cur_pheno_col, cats_fname, cat_names_flattened, raw_sample_ct, is_remove, max_thread_ct, sample_include, sample_ct_ptr);
            if (unlikely(reterr)) {
              goto KeepRemoveCats_ret_1;
            }
            success = 1;
            break;
          }
        }
      }
      if (name_blen <= max_covar_name_blen) {
        for (uint32_t covar_idx = 0; covar_idx != covar_ct; ++covar_idx) {
          if (memequal(cat_phenoname, &(covar_names[covar_idx * max_covar_name_blen]), name_blen)) {
            const PhenoCol* cur_pheno_col = &(covar_cols[covar_idx]);
            if (unlikely(cur_pheno_col->type_code != kPhenoDtypeCat)) {
              snprintf(g_logbuf, kLogbufSize, "Error: '%s' is not a categorical covariate.\n", cat_phenoname);
              goto KeepRemoveCats_ret_INCONSISTENT_INPUT_WW;
            }
            reterr = KeepRemoveCatsInternal(cur_pheno_col, cats_fname, cat_names_flattened, raw_sample_ct, is_remove, max_thread_ct, sample_include, sample_ct_ptr);
            if (unlikely(reterr)) {
              goto KeepRemoveCats_ret_1;
            }
            success = 1;
            break;
          }
        }
      }
      if (unlikely(!success)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --%s-cat-pheno phenotype/covariate not loaded.\n", is_remove? "remove" : "keep");
        goto KeepRemoveCats_ret_INCONSISTENT_INPUT_2;
      }
    }
  }
  while (0) {
  KeepRemoveCats_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
  KeepRemoveCats_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 KeepRemoveCats_ret_1:
  return reterr;
}

void ComputeAlleleFreqs(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const uint64_t* founder_allele_ddosages, uint32_t variant_ct, double af_pseudocount, double* allele_freqs) {
  // ok for maj_alleles or allele_freqs to be nullptr
  // note that founder_allele_ddosages is in 32768ths
  // could multithread this, but not a high priority since it would only tend
  // to reduce wall-clock time by a fraction of a second
  const double af_pseudocount_ddosage = af_pseudocount * u31tod(kDosageMax);
  uint32_t cur_allele_ct = 2;
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
    const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    uintptr_t allele_idx_offset_base;
    if (!allele_idx_offsets) {
      allele_idx_offset_base = 2 * variant_uidx;
    } else {
      allele_idx_offset_base = allele_idx_offsets[variant_uidx];
      cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
    }
    const uint64_t* cur_founder_allele_ddosages = &(founder_allele_ddosages[allele_idx_offset_base]);
    uint64_t tot_ddosage = 0;
    for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
      tot_ddosage += cur_founder_allele_ddosages[allele_idx];
    }
    double* cur_allele_freqs_base = &(allele_freqs[allele_idx_offset_base - variant_uidx]);
    const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
    if (!tot_ddosage) {
      const double cur_allele_ct_m1_recip = 1.0 / u31tod(cur_allele_ct);
      for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct_m1; ++allele_idx) {
        cur_allele_freqs_base[allele_idx] = cur_allele_ct_m1_recip;
      }
    } else {
      const double adj_tot_ddosage_recip = 1.0 / (u63tod(tot_ddosage) + af_pseudocount_ddosage * u31tod(cur_allele_ct));
      for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct_m1; ++allele_idx) {
        const double cur_ddosage = u63tod(cur_founder_allele_ddosages[allele_idx]) + af_pseudocount_ddosage;
        cur_allele_freqs_base[allele_idx] = cur_ddosage * adj_tot_ddosage_recip;
      }
    }
  }
}

CONSTI32(kMaxReadFreqAlleles, 255);

// relevant column types:
// 0: variant ID
// 1: ref allele code
// 2: all alt allele codes (potentially just alt1)
//
// (3-4 are --freq only)
// 3: ref freq/count
// 4: either all freqs/counts, or all-but-ref
//
// 5: obs ct (only relevant for --freq, but can be in --geno-counts)
//
// (6-11 are --geno-counts/--freqx only)
// 6: hom-ref count
// 7: het ref-alt counts (worst case, just ref-alt1)
// 8: altx-alty counts (worst case, just hom-alt1), or all pairs
// 9: hap-ref count
// 10: hap-alt counts (worst case, just hap-alt1), or all hap counts
// 11: --geno-counts numeq (if present, ignore 6..10)
//
// overrideable:
// 12->2: ALT1
// 13->4: ALT1_FREQ/ALT1_CT
// 14->7: HET_REF_ALT1_CT
// 15->8: HOM_ALT1_CT
// 16->10: HAP_ALT1_CT
ENUM_U31_DEF_START()
  kReadFreqColVarId = 0,
  kReadFreqColRefAllele,
  kReadFreqColAltAlleles,

  kReadFreqColRefFreq,
  kReadFreqColAltFreqs,

  kReadFreqColObsCt,

  kReadFreqColHomRefCt,
  kReadFreqColHetRefAltCts,
  kReadFreqColNonrefDiploidCts,
  kReadFreqColHapRefCt,
  kReadFreqColHapAltCts,
  kReadFreqColGenoCtNumeq,

  kReadFreqColAlt1Allele,
  kReadFreqColAlt1Freq,
  kReadFreqColHetRefAlt1Ct,
  kReadFreqColHomAlt1Ct,
  kReadFreqColHapAlt1Ct,

  kReadFreqColNull
ENUM_U31_DEF_END(ReadFreqColidx);

FLAGSET_DEF_START()
  kfReadFreqColset0,
  kfReadFreqColsetVarId = (1 << kReadFreqColVarId),
  kfReadFreqColsetRefAllele = (1 << kReadFreqColRefAllele),
  kfReadFreqColsetAltAlleles = (1 << kReadFreqColAltAlleles),
  kfReadFreqColsetBase = (kfReadFreqColsetVarId | kfReadFreqColsetRefAllele | kfReadFreqColsetAltAlleles),

  kfReadFreqColsetRefFreq = (1 << kReadFreqColRefFreq),
  kfReadFreqColsetAltFreqs = (1 << kReadFreqColAltFreqs),
  kfReadFreqColsetAfreqOnly = (kfReadFreqColsetRefFreq | kfReadFreqColsetAltFreqs),

  kfReadFreqColsetObsCt = (1 << kReadFreqColObsCt),

  kfReadFreqColsetHomRefCt = (1 << kReadFreqColHomRefCt),
  kfReadFreqColsetHetRefAltCts = (1 << kReadFreqColHetRefAltCts),
  kfReadFreqColsetNonrefDiploidCts = (1 << kReadFreqColNonrefDiploidCts),
  kfReadFreqColsetHapRefCt = (1 << kReadFreqColHapRefCt),
  kfReadFreqColsetHapAltCts = (1 << kReadFreqColHapAltCts),
  kfReadFreqColsetGcountDefault = ((kfReadFreqColsetHapAltCts * 2) - kfReadFreqColsetHomRefCt),

  kfReadFreqColsetGenoCtNumeq = (1 << kReadFreqColGenoCtNumeq),
  kfReadFreqColsetGcountOnly = (kfReadFreqColsetGcountDefault | kfReadFreqColsetGenoCtNumeq),

  kfReadFreqColsetAlt1Allele = (1 << kReadFreqColAlt1Allele),
  kfReadFreqColsetAlt1Freq = (1 << kReadFreqColAlt1Freq),
  kfReadFreqColsetHetRefAlt1Ct = (1 << kReadFreqColHetRefAlt1Ct),
  kfReadFreqColsetHomAlt1Ct = (1 << kReadFreqColHomAlt1Ct),
  kfReadFreqColsetHapAlt1Ct = (1 << kReadFreqColHapAlt1Ct)
FLAGSET_DEF_END(ReadFreqColFlags);

// Support exact reconstruction of original allele frequencies in --freq counts
// case.
static inline double ForceCountToDosage(double raw_count) {
  return u63tod(S_CAST(int64_t, raw_count * kDosageMax + 0.5)) * kRecipDosageMax;
}

PglErr ReadAlleleFreqs(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const char* read_freq_fname, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, double af_pseudocount, uint32_t max_thread_ct, double* allele_freqs, uintptr_t** variant_afreqcalcp) {
  // support PLINK 1.9 --freq/--freqx, and 2.0 --freq/--geno-counts.
  // GCTA-format no longer supported since it inhibits the allele consistency
  // check.
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream read_freq_txs;
  PreinitTextStream(&read_freq_txs);
  {
    if (!variant_ct) {
      logerrputs("Warning: Skipping --read-freq since no variants remain.\n");
      goto ReadAlleleFreqs_ret_1;
    }
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    double* cur_allele_freqs;
    uintptr_t* matched_loaded_alleles = nullptr; // spurious g++ 4.8 warning
    uintptr_t* matched_internal_alleles = nullptr;
    uint32_t* loaded_to_internal_allele_idx;
    uintptr_t* already_seen;
    if (unlikely(bigstack_calloc_w(raw_variant_ctl, variant_afreqcalcp) ||
                 bigstack_calloc_d(kMaxReadFreqAlleles, &cur_allele_freqs) ||
                 bigstack_alloc_w(BitCtToWordCt(kMaxReadFreqAlleles), &matched_loaded_alleles) ||
                 bigstack_alloc_w(BitCtToWordCt(max_allele_ct), &matched_internal_alleles) ||
                 bigstack_alloc_u32(kMaxReadFreqAlleles, &loaded_to_internal_allele_idx) ||
                 bigstack_calloc_w(raw_variant_ctl, &already_seen))) {
      goto ReadAlleleFreqs_ret_NOMEM;
    }
    bigstack_mark = R_CAST(unsigned char*, cur_allele_freqs);

    reterr = SizeAndInitTextStream(read_freq_fname, bigstack_left() / 8, MAXV(max_thread_ct - 1, 1), &read_freq_txs);
    if (unlikely(reterr)) {
      goto ReadAlleleFreqs_ret_TSTREAM_FAIL;
    }
    uint32_t* variant_id_htable = nullptr;
    uint32_t variant_id_htable_size;
    reterr = AllocAndPopulateIdHtableMt(variant_include, variant_ids, variant_ct, 0, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
    if (unlikely(reterr)) {
      goto ReadAlleleFreqs_ret_1;
    }
    char* line_start;
    do {
      ++line_idx;
      line_start = TextGet(&read_freq_txs);
      if (unlikely(!line_start)) {
        if (!TextStreamErrcode2(&read_freq_txs, &reterr)) {
          logerrputs("Error: Empty --read-freq file.\n");
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
        }
        goto ReadAlleleFreqs_ret_TSTREAM_FAIL;
      }
      // automatically skip header lines that start with '##' or '# '
    } while ((*line_start == '#') && (ctou32(line_start[1]) <= '#'));
    uint32_t col_skips[kReadFreqColNull];
    ReadFreqColidx col_types[kReadFreqColNull];
    uint32_t overrideable_pos[kReadFreqColNull - kReadFreqColAlt1Allele];
    uint32_t geno_counts = 0;
    uint32_t main_eq = 0;
    uint32_t is_numeq = 0;
    uint32_t use_obs_ct = 0;
    uint32_t infer_one_freq = 0;
    uint32_t infer_freq_loaded_idx = 0;
    uint32_t relevant_col_ct = 0;

    // interpretation of ColAltAlleles
    uint32_t allele_list_just_alt1 = 1;

    uint32_t is_frac = 0;  // if true, one frequency can be missing
    // could add consistency check (can't mix FREQ and CT)

    // interpretation of ColAltFreqs and ColNonrefDiploidCts
    uint32_t main_allele_idx_start = 1;
    uint32_t main_list_just_alt1 = 1;

    // interpretation of ColHapAltCts
    uint32_t hap_allele_idx_start = 1;
    uint32_t hap_list_just_alt1 = 1;

    uint32_t het_list_just_alt1 = 1;  // ColHetRefAltCts

    uint32_t biallelic_only = 0;

    ReadFreqColFlags header_cols = kfReadFreqColset0;
    uint32_t skip_header = 1;
    if (*line_start == '#') {
      // PLINK 2.0
      // guaranteed nonspace
      const char* linebuf_iter = &(line_start[1]);

      uint32_t col_idx = 0;
      while (1) {
        const char* token_end = CurTokenEnd(linebuf_iter);
        const uint32_t token_slen = token_end - linebuf_iter;
        ReadFreqColidx cur_colidx = kReadFreqColNull;
        if (token_slen <= 4) {
          if (strequal_k(linebuf_iter, "ID", token_slen)) {
            cur_colidx = kReadFreqColVarId;
          } else if (token_slen == 3) {
            if (memequal_sk(linebuf_iter, "REF")) {
              cur_colidx = kReadFreqColRefAllele;
            } else if (memequal_sk(linebuf_iter, "ALT")) {
              cur_colidx = kReadFreqColAltAlleles;
              if (allele_list_just_alt1) {
                header_cols &= ~kfReadFreqColsetAlt1Allele;
                allele_list_just_alt1 = 0;
              }
            } else if (memequal_sk(linebuf_iter, "CTS")) {
              goto ReadAlleleFreqs_freqmain_found1;
            }
          } else if (strequal_k(linebuf_iter, "ALT1", token_slen) && allele_list_just_alt1) {
            cur_colidx = kReadFreqColAlt1Allele;
          }
        } else if (strequal_k(linebuf_iter, "REF_FREQ", token_slen) ||
                   strequal_k(linebuf_iter, "REF_CT", token_slen)) {
          cur_colidx = kReadFreqColRefFreq;
          if (linebuf_iter[4] == 'F') {
            is_frac = 1;
          }
        } else if ((strequal_k(linebuf_iter, "ALT1_FREQ", token_slen) ||
                    strequal_k(linebuf_iter, "ALT1_CT", token_slen)) &&
                   main_list_just_alt1) {
          cur_colidx = kReadFreqColAlt1Freq;
          if (linebuf_iter[5] == 'F') {
            is_frac = 1;
          }
        } else if (strequal_k(linebuf_iter, "ALT_FREQS", token_slen) ||
                   strequal_k(linebuf_iter, "ALT_CTS", token_slen)) {
          if (linebuf_iter[4] == 'F') {
            is_frac = 1;
          }
          goto ReadAlleleFreqs_freqmain_found2;
        } else if (strequal_k(linebuf_iter, "FREQS", token_slen)) {
          is_frac = 1;
          goto ReadAlleleFreqs_freqmain_found1;
        } else if (strequal_k(linebuf_iter, "ALT_NUM_FREQS", token_slen) ||
                   strequal_k(linebuf_iter, "ALT_NUM_CTS", token_slen)) {
          is_numeq = 1;
          goto ReadAlleleFreqs_freqmain_found2;
        } else if (strequal_k(linebuf_iter, "NUM_FREQS", token_slen) ||
                   strequal_k(linebuf_iter, "NUM_CTS", token_slen)) {
          is_numeq = 1;
        ReadAlleleFreqs_freqmain_found1:
          main_allele_idx_start = 0;
        ReadAlleleFreqs_freqmain_found2:
          cur_colidx = kReadFreqColAltFreqs;
          if (main_list_just_alt1) {
            header_cols &= ~kfReadFreqColsetAlt1Freq;
            main_list_just_alt1 = 0;
          }
        } else if (strequal_k(linebuf_iter, "OBS_CT", token_slen)) {
          cur_colidx = kReadFreqColObsCt;
        } else if (strequal_k(linebuf_iter, "HOM_REF_CT", token_slen)) {
          cur_colidx = kReadFreqColHomRefCt;
        } else if (strequal_k(linebuf_iter, "HET_REF_ALT1_CT", token_slen) && het_list_just_alt1) {
          cur_colidx = kReadFreqColHetRefAlt1Ct;
        } else if (strequal_k(linebuf_iter, "HET_REF_ALT_CTS", token_slen)) {
          cur_colidx = kReadFreqColHetRefAltCts;
          if (het_list_just_alt1) {
            header_cols &= ~kfReadFreqColsetHetRefAlt1Ct;
            het_list_just_alt1 = 0;
          }
        } else if (strequal_k(linebuf_iter, "HOM_ALT1_CT", token_slen) && main_list_just_alt1) {
          cur_colidx = kReadFreqColHomAlt1Ct;
        } else if (strequal_k(linebuf_iter, "TWO_ALT_GENO_CTS", token_slen) ||
                   strequal_k(linebuf_iter, "NONREF_DIPLOID_GENO_CTS", token_slen)) {
          goto ReadAlleleFreqs_countmain_found;
        } else if (strequal_k(linebuf_iter, "DIPLOID_GENO_CTS", token_slen)) {
          main_allele_idx_start = 0;
        ReadAlleleFreqs_countmain_found:
          cur_colidx = kReadFreqColNonrefDiploidCts;
          if (main_list_just_alt1) {
            header_cols &= ~kfReadFreqColsetHomAlt1Ct;
            // could make this use a different variable than FREQS does
            main_list_just_alt1 = 0;
          }
        } else if (strequal_k(linebuf_iter, "HAP_REF_CT", token_slen)) {
          cur_colidx = kReadFreqColHapRefCt;
        } else if (strequal_k(linebuf_iter, "HAP_ALT1_CT", token_slen) && hap_list_just_alt1) {
          cur_colidx = kReadFreqColHapAlt1Ct;
        } else if (strequal_k(linebuf_iter, "HAP_ALT_CTS", token_slen)) {
          goto ReadAlleleFreqs_hapmain_found;
        } else if (strequal_k(linebuf_iter, "HAP_CTS", token_slen)) {
          hap_allele_idx_start = 0;
        ReadAlleleFreqs_hapmain_found:
          cur_colidx = kReadFreqColHapAltCts;
          if (hap_list_just_alt1) {
            header_cols &= ~kfReadFreqColsetHapAlt1Ct;
            hap_list_just_alt1 = 0;
          }
        } else if (strequal_k(linebuf_iter, "GENO_NUM_CTS", token_slen)) {
          cur_colidx = kReadFreqColGenoCtNumeq;
          is_numeq = 1;
        }
        if (cur_colidx != kReadFreqColNull) {
          const ReadFreqColFlags cur_colset = S_CAST(ReadFreqColFlags, 1U << cur_colidx);
          if (unlikely(header_cols & cur_colset)) {
            logerrputs("Error: Conflicting columns in header line of --read-freq file.\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
          if (cur_colidx >= kReadFreqColAlt1Allele) {
            overrideable_pos[cur_colidx - kReadFreqColAlt1Allele] = relevant_col_ct;
          }
          header_cols |= cur_colset;
          col_skips[relevant_col_ct] = col_idx;
          col_types[relevant_col_ct++] = cur_colidx;
        }
        linebuf_iter = FirstNonTspace(token_end);
        if (IsEolnKns(*linebuf_iter)) {
          break;
        }
        ++col_idx;
      }
      ReadFreqColFlags semifinal_header_cols = header_cols;
      if (header_cols & kfReadFreqColsetAlt1Allele) {
        header_cols ^= kfReadFreqColsetAltAlleles | kfReadFreqColsetAlt1Allele;
        col_types[overrideable_pos[0]] = kReadFreqColAltAlleles;
      }
      if (header_cols & kfReadFreqColsetAlt1Freq) {
        header_cols ^= kfReadFreqColsetAltFreqs | kfReadFreqColsetAlt1Freq;
        col_types[overrideable_pos[kReadFreqColAlt1Freq - kReadFreqColAlt1Allele]] = kReadFreqColAltFreqs;
      }
      if (header_cols & kfReadFreqColsetHetRefAlt1Ct) {
        header_cols ^= kfReadFreqColsetHetRefAltCts | kfReadFreqColsetHetRefAlt1Ct;
        col_types[overrideable_pos[kReadFreqColHetRefAlt1Ct - kReadFreqColAlt1Allele]] = kReadFreqColHetRefAltCts;
      }
      if (header_cols & kfReadFreqColsetHomAlt1Ct) {
        header_cols ^= kfReadFreqColsetNonrefDiploidCts | kfReadFreqColsetHomAlt1Ct;
        col_types[overrideable_pos[kReadFreqColHomAlt1Ct - kReadFreqColAlt1Allele]] = kReadFreqColNonrefDiploidCts;
      }
      if (header_cols & kfReadFreqColsetHapAlt1Ct) {
        header_cols ^= kfReadFreqColsetHapAltCts | kfReadFreqColsetHapAlt1Ct;
        col_types[overrideable_pos[kReadFreqColHapAlt1Ct - kReadFreqColAlt1Allele]] = kReadFreqColHapAltCts;
      }
      if ((semifinal_header_cols != header_cols) && (!(header_cols & kfReadFreqColsetGenoCtNumeq))) {
        // we're treating at least one ALT1 column as if it spoke for all ALT
        // alleles
        biallelic_only = 1;
      }

      main_eq = is_numeq;
      semifinal_header_cols = header_cols;
      if (header_cols & kfReadFreqColsetAfreqOnly) {
        if (unlikely(header_cols & kfReadFreqColsetGcountOnly)) {
          logerrputs("Error: Conflicting columns in header line of --read-freq file (--freq and\n--geno-counts values mixed together).\n");
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
        }
        ReadFreqColFlags header_cols_exempt = kfReadFreqColset0;
        if ((header_cols & kfReadFreqColsetAltFreqs) && (!is_numeq)) {
          // [ALT_]FREQS can be formatted as either
          //   0.5,0,0.2
          // or
          //   A=0.5,G=0.2
          // Look at the first nonheader line to distinguish between these two.
          ++line_idx;
          line_start = TextGet(&read_freq_txs);
          if (unlikely(!line_start)) {
            if (!TextStreamErrcode2(&read_freq_txs, &reterr)) {
              logerrputs("Error: Empty --read-freq file.\n");
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
            }
            goto ReadAlleleFreqs_ret_TSTREAM_FAIL;
          }
          linebuf_iter = line_start;
          const char* alt_freq_str = nullptr;
          for (uint32_t relevant_col_idx = 0; relevant_col_idx != relevant_col_ct; ++relevant_col_idx) {
            if (col_types[relevant_col_idx] == kReadFreqColAltFreqs) {
              alt_freq_str = NextTokenMult0(linebuf_iter, col_skips[relevant_col_idx]);
              break;
            }
          }
          if (unlikely(!alt_freq_str)) {
            goto ReadAlleleFreqs_ret_MISSING_TOKENS;
          }
          const uint32_t alt_freq_slen = CurTokenEnd(alt_freq_str) - alt_freq_str;
          // bare '.' can only appear in eq formats
          main_eq = ((alt_freq_slen == 1) && (*alt_freq_str == '.')) || (memchr(alt_freq_str, '=', alt_freq_slen) != nullptr);
          if (main_eq) {
            header_cols_exempt = kfReadFreqColsetAltAlleles;
            if (!main_allele_idx_start) {
              header_cols_exempt |= kfReadFreqColsetRefAllele;
            }
            header_cols &= ~header_cols_exempt;
          }
          skip_header = 0;
        }
        if (unlikely(((header_cols & kfReadFreqColsetBase) | header_cols_exempt) != kfReadFreqColsetBase)) {
          logerrputs("Error: Missing column(s) in --read-freq file (ID, REF, ALT[1] usually\nrequired).\n");
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
        }
        if (!main_allele_idx_start) {
          header_cols &= ~kfReadFreqColsetRefFreq;
        } else {
          if ((header_cols & (kfReadFreqColsetRefFreq | kfReadFreqColsetAltFreqs)) != (kfReadFreqColsetRefFreq | kfReadFreqColsetAltFreqs)) {
            if (main_list_just_alt1) {
              biallelic_only = 1;
            }
            infer_one_freq = 1;
            infer_freq_loaded_idx = (header_cols / kfReadFreqColsetRefFreq) & 1;
            if (!is_frac) {
              if (unlikely(!(header_cols & kfReadFreqColsetObsCt))) {
                logerrputs("Error: Missing column(s) in --read-freq file (at least two of {REF_CT, ALT1_CT,\nALT_CTS, OBS_CT} must be present).\n");
                goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
              }
              use_obs_ct = 1;
            }
          }
        }
        logputs("--read-freq: PLINK 2 --freq file detected.\n");
      } else if (likely(header_cols & kfReadFreqColsetGcountOnly)) {
        if (unlikely((header_cols & kfReadFreqColsetBase) != kfReadFreqColsetBase)) {
          logerrputs("Error: Missing column(s) in --read-freq file (ID, REF, ALT[1] required).\n");
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
        }
        // possible todo: allow one frequency/count to be missing.  (not really
        // necessary since PLINK 1.9 --freqx does not leave anything out,
        // unlike PLINK 1.x --freq)
        if (header_cols & kfReadFreqColsetGenoCtNumeq) {
          // don't need anything but GENO_NUM_CTS
          header_cols &= ~kfReadFreqColsetGcountDefault;
        } else {
          // require both diploid and haploid columns for now.  (could
          // conditionally drop one of these requirements later.)
          if (unlikely(!(header_cols & kfReadFreqColsetNonrefDiploidCts))) {
            logerrputs("Error: Missing column(s) in --read-freq file (HOM_ALT1_CT,\nTWO_ALT_GENO_CTS, or DIPLOID_GENO_CTS required).\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
          if (!main_allele_idx_start) {
            header_cols &= ~(kfReadFreqColsetHomRefCt | kfReadFreqColsetHetRefAltCts);
          } else if (unlikely((header_cols & (kfReadFreqColsetHomRefCt | kfReadFreqColsetHetRefAltCts)) != (kfReadFreqColsetHomRefCt | kfReadFreqColsetHetRefAltCts))) {
            logerrputs("Error: Missing column(s) in --read-freq file (HOM_REF_CT, HET_REF_ALT1_CT, or\nHET_REF_ALT_CTS required unless {DIPLOID_}GENO_CTS present).\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
          if (unlikely(!(header_cols & kfReadFreqColsetHapAltCts))) {
            logerrputs("Error: Missing column(s) in --read-freq file (HAP_ALT1_CT, HAP_ALT_CTS, or\nHAP_CTS required).\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
          if (!hap_allele_idx_start) {
            header_cols &= ~kfReadFreqColsetHapRefCt;
          } else if (unlikely(!(header_cols & kfReadFreqColsetHapRefCt))) {
            logerrputs("Error: Missing column(s) in --read-freq file (HAP_REF_CT required unless\nHAP_CTS or GENO_CTS present).\n");
            goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
          }
        }
        geno_counts = 1;
        logputs("--read-freq: PLINK 2 --geno-counts file detected.\n");
      } else {
        logerrputs("Error: Missing column(s) in --read-freq file (no frequencies/counts).\n");
        goto ReadAlleleFreqs_ret_MALFORMED_INPUT;
      }
      if (!use_obs_ct) {
        header_cols &= ~kfReadFreqColsetObsCt;
      }
      if (semifinal_header_cols != header_cols) {
        // remove redundant columns
        uint32_t relevant_col_idx_read = 0;
        while ((S_CAST(uint32_t, header_cols) >> col_types[relevant_col_idx_read]) & 1) {
          ++relevant_col_idx_read;
        }
        uint32_t relevant_col_idx_write = relevant_col_idx_read++;
        for (; relevant_col_idx_read != relevant_col_ct; ++relevant_col_idx_read) {
          const ReadFreqColidx cur_colidx = col_types[relevant_col_idx_read];
          if ((S_CAST(uint32_t, header_cols) >> cur_colidx) & 1) {
            col_types[relevant_col_idx_write] = cur_colidx;
            col_skips[relevant_col_idx_write] = col_skips[relevant_col_idx_read];
            ++relevant_col_idx_write;
          }
        }
        relevant_col_ct = relevant_col_idx_write;
      }
      for (uint32_t uii = relevant_col_ct - 1; uii; --uii) {
        col_skips[uii] -= col_skips[uii - 1];
      }
    } else {
      // PLINK 1.x
      // .frq:       CHR  SNP  A1  A2  MAF        NCHROBS
      // .frq.count: CHR  SNP  A1  A2  C1         C2       G0
      // .frqx:      CHR  SNP  A1  A2  C(HOM A1)  C(HET)   C(HOM A2)  C(HAP A1)
      //   C(HAP A2)  C(MISSING)
      // (yeah, the spaces in the .frqx header were a mistake, should have used
      // underscores.  oh well, live and learn.)
      col_skips[0] = 1;
      col_skips[1] = 1;
      col_skips[2] = 1;
      col_skips[3] = 1;

      col_types[0] = kReadFreqColVarId;
      // doesn't matter if we treat A1 or A2 as ref
      col_types[1] = kReadFreqColRefAllele;
      col_types[2] = kReadFreqColAltAlleles;
      biallelic_only = 1;
      if (StrStartsWithUnsafe(line_start, "CHR\tSNP\tA1\tA2\tC(HOM A1)\tC(HET)\tC(HOM A2)\tC(HAP A1)\tC(HAP A2)\tC(MISSING)")) {
        col_skips[4] = 1;
        col_skips[5] = 1;
        col_skips[6] = 1;
        col_skips[7] = 1;

        col_types[3] = kReadFreqColHomRefCt;
        col_types[4] = kReadFreqColHetRefAltCts;
        col_types[5] = kReadFreqColNonrefDiploidCts;
        col_types[6] = kReadFreqColHapRefCt;
        col_types[7] = kReadFreqColHapAltCts;
        header_cols = kfReadFreqColsetBase | kfReadFreqColsetGcountOnly;
        geno_counts = 1;
        relevant_col_ct = 8;
        logputs("--read-freq: PLINK 1.9 --freqx file detected.\n");
      } else {
        if (unlikely(!tokequal_k(line_start, "CHR"))) {
          goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
        }
        const char* linebuf_iter = FirstNonTspace(&(line_start[3]));
        if (unlikely(!tokequal_k(linebuf_iter, "SNP"))) {
          goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
        }
        linebuf_iter = FirstNonTspace(&(linebuf_iter[3]));
        if (unlikely(!tokequal_k(linebuf_iter, "A1"))) {
          goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
        }
        linebuf_iter = FirstNonTspace(&(linebuf_iter[2]));
        if (unlikely(!tokequal_k(linebuf_iter, "A2"))) {
          goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
        }
        linebuf_iter = FirstNonTspace(&(linebuf_iter[2]));
        col_types[3] = kReadFreqColRefFreq;
        if (tokequal_k(linebuf_iter, "MAF")) {
          is_frac = 1;
          infer_one_freq = 1;
          infer_freq_loaded_idx = 1;
          header_cols = kfReadFreqColsetBase | kfReadFreqColsetRefFreq;
          relevant_col_ct = 4;
          logputs("--read-freq: PLINK 1.x --freq file detected.\n");
        } else {
          if (unlikely(!tokequal_k(linebuf_iter, "C1"))) {
            goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
          }
          linebuf_iter = FirstNonTspace(&(linebuf_iter[2]));
          if (unlikely(!tokequal_k(linebuf_iter, "C2"))) {
            goto ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER;
          }
          col_skips[4] = 1;
          col_types[4] = kReadFreqColAltFreqs;
          header_cols = kfReadFreqColsetBase | kfReadFreqColsetAfreqOnly;
          relevant_col_ct = 5;
          logputs("--read-freq: PLINK 1.x '--freq counts' file detected.\n");
        }
      }
    }
    assert(relevant_col_ct <= 8);

    double freq_max = 4294967295.0;
    if (is_frac) {
      af_pseudocount = 0.0;
      freq_max = 1.0;
    }
    uintptr_t skipped_variant_ct = 0;
    uint32_t loaded_variant_ct = 0;
    uint32_t cur_allele_ct = 2;
    uint32_t variant_uidx = 0;
    char* line_iter = line_start;
    if (skip_header) {
      line_iter = TextLineEnd(&read_freq_txs);
      ++line_idx;
    }
    for (; TextGetUnsafe2(&read_freq_txs, &line_iter); ++line_idx) {
      {
        // not const since tokens may be null-terminated or comma-terminated
        // later
        char* token_ptrs[12];
        uint32_t token_slens[12];
        line_iter = TokenLex0(line_iter, R_CAST(uint32_t*, col_types), col_skips, relevant_col_ct, token_ptrs, token_slens);
        if (unlikely(!line_iter)) {
          goto ReadAlleleFreqs_ret_MISSING_TOKENS;
        }
        // characters may be modified, so better find the \n now just in case
        // it gets clobbered
        line_iter = AdvPastDelim(line_iter, '\n');
        const char* variant_id_start = token_ptrs[kReadFreqColVarId];
        const uint32_t variant_id_slen = token_slens[kReadFreqColVarId];
        variant_uidx = VariantIdDupflagHtableFind(variant_id_start, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
        if (variant_uidx >> 31) {
          if (likely(variant_uidx == UINT32_MAX)) {
            ++skipped_variant_ct;
            continue;
          }
          snprintf(g_logbuf, kLogbufSize, "Error: --read-freq variant ID '%s' appears multiple times in main dataset.\n", variant_ids[variant_uidx & 0x7fffffff]);
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
        }
        if (unlikely(IsSet(already_seen, variant_uidx))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Variant ID '%s' appears multiple times in --read-freq file.\n", variant_ids[variant_uidx]);
          goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
        }
        SetBit(variant_uidx, already_seen);

        uintptr_t allele_idx_offset_base;
        if (!allele_idx_offsets) {
          allele_idx_offset_base = variant_uidx * 2;
        } else {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
          if (biallelic_only) {
            goto ReadAlleleFreqs_skip_variant;
          }
        }
        ZeroWArr(BitCtToWordCt(cur_allele_ct), matched_internal_alleles);
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        uint32_t loaded_allele_idx_end = 0;
        {
          uint32_t unmatched_allele_ct = cur_allele_ct;
          uint32_t cur_loaded_allele_code_slen;
          char* cur_loaded_allele_code;
          char* loaded_allele_code_iter;
          char* loaded_allele_code_end;
          if (header_cols & kfReadFreqColsetRefAllele) {
            cur_loaded_allele_code_slen = token_slens[kReadFreqColRefAllele];
            cur_loaded_allele_code = token_ptrs[kReadFreqColRefAllele];
            // No special handling of missing allele code needed (if the
            // fileset is valid).
            cur_loaded_allele_code[cur_loaded_allele_code_slen] = '\0';
            if (header_cols & kfReadFreqColsetAltAlleles) {
              loaded_allele_code_iter = token_ptrs[kReadFreqColAltAlleles];
              loaded_allele_code_end = &(loaded_allele_code_iter[token_slens[kReadFreqColAltAlleles]]);
              *loaded_allele_code_end++ = ',';
            } else {
              // special case: with --freq alteq or alteqz column, we only need
              // to scrape REF here
              loaded_allele_code_iter = &(cur_loaded_allele_code[cur_loaded_allele_code_slen + 1]);
              loaded_allele_code_end = loaded_allele_code_iter;
            }
          } else {
            // no REF, only ALT.  skip the first iteration of the loop.
            loaded_allele_code_iter = token_ptrs[kReadFreqColAltAlleles];
            loaded_allele_code_end = &(loaded_allele_code_iter[token_slens[kReadFreqColAltAlleles]]);
            *loaded_allele_code_end++ = ',';

            cur_loaded_allele_code = loaded_allele_code_iter;
            loaded_allele_code_iter = S_CAST(char*, memchr(loaded_allele_code_iter, ',', loaded_allele_code_end - cur_loaded_allele_code));
            cur_loaded_allele_code_slen = loaded_allele_code_iter - cur_loaded_allele_code;
            *loaded_allele_code_iter++ = '\0';
            loaded_allele_idx_end = 1;
            matched_loaded_alleles[0] = 0;
          }
          uint32_t widx = 0;
          while (1) {
            if (!(loaded_allele_idx_end % kBitsPerWord)) {
              widx = loaded_allele_idx_end / kBitsPerWord;
              matched_loaded_alleles[widx] = 0;
            }
            if (cur_loaded_allele_code_slen <= max_allele_slen) {
              const uint32_t cur_blen = cur_loaded_allele_code_slen + 1;
              uintptr_t internal_allele_idx_base = 0;
              uintptr_t cur_inv_bits = ~matched_internal_alleles[0];
              for (uint32_t unmatched_allele_idx = 0; unmatched_allele_idx != unmatched_allele_ct; ++unmatched_allele_idx) {
                const uintptr_t internal_allele_idx = BitIter0(matched_internal_alleles, &internal_allele_idx_base, &cur_inv_bits);
                if (memequal(cur_loaded_allele_code, cur_alleles[internal_allele_idx], cur_blen)) {
                  if (unlikely(IsSet(matched_internal_alleles, internal_allele_idx))) {
                    snprintf(g_logbuf, kLogbufSize, "Error: Duplicate allele code on line %" PRIuPTR " of --read-freq file.\n", line_idx);
                    goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
                  }
                  SetBit(internal_allele_idx, matched_internal_alleles);
                  SetBit(loaded_allele_idx_end, matched_loaded_alleles);
                  loaded_to_internal_allele_idx[loaded_allele_idx_end] = internal_allele_idx;
                  break;
                }
              }
            }
            ++loaded_allele_idx_end;
            if (loaded_allele_code_iter == loaded_allele_code_end) {
              break;
            }
            if (unlikely(loaded_allele_idx_end == kMaxReadFreqAlleles)) {
              snprintf(g_logbuf, kLogbufSize, "Error: --read-freq file entry for variant ID '%s' has more than %u ALT alleles.\n", variant_ids[variant_uidx], kMaxReadFreqAlleles - 1);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
            }
            cur_loaded_allele_code = loaded_allele_code_iter;
            loaded_allele_code_iter = S_CAST(char*, memchr(loaded_allele_code_iter, ',', loaded_allele_code_end - cur_loaded_allele_code));
            cur_loaded_allele_code_slen = loaded_allele_code_iter - cur_loaded_allele_code;
            *loaded_allele_code_iter++ = '\0';
          }
        }
        // bugfix (31 May 2023): forgot to skip variant when not enough alleles
        // were loaded.
        {
          const uint32_t matched_allele_ct = PopcountWords(matched_loaded_alleles, BitCtToWordCt(loaded_allele_idx_end));
          if (matched_allele_ct != cur_allele_ct) {
            assert(matched_allele_ct < cur_allele_ct);
            goto ReadAlleleFreqs_skip_variant;
          }
        }

        double* allele_freqs_write = &(allele_freqs[allele_idx_offset_base - variant_uidx]);
        if (geno_counts) {
          ZeroDArr(cur_allele_ct, cur_allele_freqs);
          if (is_numeq) {
            const uint32_t full_slen = token_slens[kReadFreqColGenoCtNumeq];
            char* geno_num_cts = token_ptrs[kReadFreqColGenoCtNumeq];
            if (full_slen > 1) {
              geno_num_cts[full_slen] = ',';
              const char* geno_num_cts_iter = geno_num_cts;
              const char* geno_num_cts_end = &(geno_num_cts[full_slen]);
#ifndef __LP64__
              const uint32_t cap_div_10 = (loaded_allele_idx_end - 1) / 10;
              const uint32_t cap_mod_10 = (loaded_allele_idx_end - 1) % 10;
#endif
              while (1) {
                uint32_t second_loaded_allele_idx = UINT32_MAX;
                uint32_t first_loaded_allele_idx;
#ifdef __LP64__
                if (unlikely(ScanmovUintCapped(loaded_allele_idx_end - 1, &geno_num_cts_iter, &first_loaded_allele_idx))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (*geno_num_cts_iter == '/') {
                  ++geno_num_cts_iter;
                  if (unlikely(ScanmovUintCapped(loaded_allele_idx_end - 1, &geno_num_cts_iter, &second_loaded_allele_idx))) {
                    goto ReadAlleleFreqs_ret_INVALID_FREQS;
                  }
                }
#else
                if (unlikely(ScanmovUintCapped32(cap_div_10, cap_mod_10, &geno_num_cts_iter, &first_loaded_allele_idx))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (*geno_num_cts_iter == '/') {
                  ++geno_num_cts_iter;
                  if (unlikely(ScanmovUintCapped32(cap_div_10, cap_mod_10, &geno_num_cts_iter, &second_loaded_allele_idx))) {
                    goto ReadAlleleFreqs_ret_INVALID_FREQS;
                  }
                }
#endif
                if (unlikely(*geno_num_cts_iter != '=')) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                ++geno_num_cts_iter;
                double dxx;
                const char* cur_ct_end = ScanadvDouble(geno_num_cts_iter, &dxx);
                if (unlikely((!cur_ct_end) || (*cur_ct_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                dxx = ForceCountToDosage(dxx);
                if (IsSet(matched_loaded_alleles, first_loaded_allele_idx)) {
                  cur_allele_freqs[loaded_to_internal_allele_idx[first_loaded_allele_idx]] += dxx;
                }
                if ((second_loaded_allele_idx != UINT32_MAX) && IsSet(matched_loaded_alleles, second_loaded_allele_idx)) {
                  cur_allele_freqs[loaded_to_internal_allele_idx[second_loaded_allele_idx]] += dxx;
                }
                geno_num_cts_iter = cur_ct_end;
                if (geno_num_cts_iter == geno_num_cts_end) {
                  break;
                }
                ++geno_num_cts_iter;
              }
            } else if (unlikely(*geno_num_cts != '.')) {
              goto ReadAlleleFreqs_ret_INVALID_FREQS;
            }
          } else {
            const uint32_t internal0 = IsSet(matched_loaded_alleles, 0)? loaded_to_internal_allele_idx[0] : UINT32_MAX;
            if (header_cols & kfReadFreqColsetHomRefCt) {
              if (internal0 != UINT32_MAX) {
                const char* hom_ref_str = token_ptrs[kReadFreqColHomRefCt];
                double dxx;
                const char* hom_ref_end = ScantokDouble(hom_ref_str, &dxx);
                if (unlikely((!hom_ref_end) || (dxx < 0.0) || (dxx > 4294967295.0))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                dxx = ForceCountToDosage(dxx);
                cur_allele_freqs[internal0] += 2 * dxx;
              }

              char* het_refalt = token_ptrs[kReadFreqColHetRefAltCts];
              const uint32_t het_refalt_slen = token_slens[kReadFreqColHetRefAltCts];
              het_refalt[het_refalt_slen] = ',';
              const char* het_refalt_iter = het_refalt;
              const char* het_refalt_end = &(het_refalt[het_refalt_slen]);
              for (uint32_t alt_allele_idx = 1; alt_allele_idx != cur_allele_ct; ++alt_allele_idx) {
                if (unlikely(het_refalt_iter >= het_refalt_end)) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                double dxx;
                const char* cur_entry_end = ScanadvDouble(het_refalt_iter, &dxx);
                if (unlikely((!cur_entry_end) || (*cur_entry_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                dxx = ForceCountToDosage(dxx);
                if (internal0 != UINT32_MAX) {
                  cur_allele_freqs[internal0] += dxx;
                }
                if (IsSet(matched_loaded_alleles, alt_allele_idx)) {
                  cur_allele_freqs[loaded_to_internal_allele_idx[alt_allele_idx]] += dxx;
                }
                het_refalt_iter = &(cur_entry_end[1]);
              }
            }
            // ColNonrefDiploidCts required
            char* diploid_cts = token_ptrs[kReadFreqColNonrefDiploidCts];
            const uint32_t diploid_cts_slen = token_slens[kReadFreqColNonrefDiploidCts];
            diploid_cts[diploid_cts_slen] = ',';
            const char* diploid_cts_iter = diploid_cts;
            const char* diploid_cts_end = &(diploid_cts[diploid_cts_slen]);
            for (uint32_t second_allele_idx = main_allele_idx_start; second_allele_idx != cur_allele_ct; ++second_allele_idx) {
              uint32_t internalx = UINT32_MAX;
              if (IsSet(matched_loaded_alleles, second_allele_idx)) {
                internalx = loaded_to_internal_allele_idx[second_allele_idx];
              }
              // 1/1, 1/2, 2/2, 1/3, ...
              for (uint32_t first_allele_idx = main_allele_idx_start; first_allele_idx <= second_allele_idx; ++first_allele_idx) {
                if (unlikely(diploid_cts_iter >= diploid_cts_end)) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                double dxx;
                const char* cur_entry_end = ScanadvDouble(diploid_cts_iter, &dxx);
                if (unlikely((!cur_entry_end) || (*cur_entry_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                dxx = ForceCountToDosage(dxx);
                if (IsSet(matched_loaded_alleles, first_allele_idx)) {
                  cur_allele_freqs[loaded_to_internal_allele_idx[first_allele_idx]] += dxx;
                }
                if (internalx != UINT32_MAX) {
                  cur_allele_freqs[internalx] += dxx;
                }
                diploid_cts_iter = &(cur_entry_end[1]);
              }
            }

            if ((header_cols & kfReadFreqColsetHapRefCt) && (internal0 != UINT32_MAX)) {
              const char* hap_ref_str = token_ptrs[kReadFreqColHapRefCt];
              double dxx;
              const char* hap_ref_end = ScantokDouble(hap_ref_str, &dxx);
              if (unlikely((!hap_ref_end) || (dxx < 0.0) || (dxx > 4294967295.0))) {
                goto ReadAlleleFreqs_ret_INVALID_FREQS;
              }
              dxx = ForceCountToDosage(dxx);
              cur_allele_freqs[internal0] += dxx;
            }
            // ColHapAltCts required
            char* hap_alt = token_ptrs[kReadFreqColHapAltCts];
            const uint32_t hap_alt_slen = token_slens[kReadFreqColHapAltCts];
            hap_alt[hap_alt_slen] = ',';
            const char* hap_alt_iter = hap_alt;
            const char* hap_alt_end = &(hap_alt[hap_alt_slen]);
            for (uint32_t alt_allele_idx = 1; alt_allele_idx != cur_allele_ct; ++alt_allele_idx) {
              if (unlikely(hap_alt_iter >= hap_alt_end)) {
                goto ReadAlleleFreqs_ret_INVALID_FREQS;
              }
              double dxx;
              const char* cur_entry_end = ScanadvDouble(hap_alt_iter, &dxx);
              if (unlikely((!cur_entry_end) || (*cur_entry_end != ',') || (dxx < 0.0) || (dxx > 4294967295.0))) {
                goto ReadAlleleFreqs_ret_INVALID_FREQS;
              }
              if (IsSet(matched_loaded_alleles, alt_allele_idx)) {
                dxx = ForceCountToDosage(dxx);
                cur_allele_freqs[loaded_to_internal_allele_idx[alt_allele_idx]] += dxx;
              }
              hap_alt_iter = &(cur_entry_end[1]);
            }
          }
        } else {
          if ((header_cols & kfReadFreqColsetRefFreq) && IsSet(matched_loaded_alleles, 0)) {
            const char* ref_freq_str = token_ptrs[kReadFreqColRefFreq];
            double dxx;
            if (!ScantokDouble(ref_freq_str, &dxx)) {
              if (likely(IsNanStr(ref_freq_str, token_slens[kReadFreqColRefFreq]))) {
                goto ReadAlleleFreqs_skip_variant;
              }
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid REF frequency/count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
            }
            if (unlikely((dxx < 0.0) || (dxx > freq_max))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid REF frequency/count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_WW;
            }
            if (!is_frac) {
              dxx = ForceCountToDosage(dxx);
            }
            cur_allele_freqs[loaded_to_internal_allele_idx[0]] = dxx;
          }
          if (header_cols & kfReadFreqColsetAltFreqs) {
            const uint32_t full_slen = token_slens[kReadFreqColAltFreqs];
            char* alt_freq_start = token_ptrs[kReadFreqColAltFreqs];
            alt_freq_start[full_slen] = ',';
            if (!main_eq) {
              const char* alt_freq_iter = alt_freq_start;
              const char* alt_freq_end = &(alt_freq_start[full_slen]);
              for (uint32_t allele_idx = main_allele_idx_start; allele_idx != loaded_allele_idx_end; ++allele_idx, ++alt_freq_iter) {
                if (unlikely(alt_freq_iter >= alt_freq_end)) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (!IsSet(matched_loaded_alleles, allele_idx)) {
                  alt_freq_iter = AdvToDelim(alt_freq_iter, ',');
                  continue;
                }
                double dxx;
                const char* cur_freq_end = ScanadvDouble(alt_freq_iter, &dxx);
                if (!cur_freq_end) {
                  cur_freq_end = AdvToDelim(alt_freq_iter, ',');
                  if (likely(IsNanStr(alt_freq_iter, cur_freq_end - alt_freq_iter))) {
                    goto ReadAlleleFreqs_skip_variant;
                  }
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (unlikely((*cur_freq_end != ',') || (dxx < 0.0) || (dxx > freq_max))) {
                  goto ReadAlleleFreqs_ret_INVALID_FREQS;
                }
                if (!is_frac) {
                  dxx = ForceCountToDosage(dxx);
                }
                alt_freq_iter = cur_freq_end;
                cur_allele_freqs[loaded_to_internal_allele_idx[allele_idx]] = dxx;
              }
            } else {
              ZeroDArr(cur_allele_ct, cur_allele_freqs);
              if ((full_slen > 1) || (*alt_freq_start != '.')) {
                if (is_numeq) {
                  const char* alt_freq_iter = alt_freq_start;
                  const char* alt_freq_end = &(alt_freq_start[full_slen]);
#ifndef __LP64__
                  const uint32_t cap_div_10 = (loaded_allele_idx_end - 1) / 10;
                  const uint32_t cap_mod_10 = (loaded_allele_idx_end - 1) % 10;
#endif
                  while (1) {
                    const char* cur_entry_end = AdvToDelim(alt_freq_iter, ',');
                    uint32_t loaded_allele_idx;
#ifdef __LP64__
                    if (unlikely(ScanmovUintCapped(loaded_allele_idx_end - 1, &alt_freq_iter, &loaded_allele_idx))) {
                      goto ReadAlleleFreqs_ret_INVALID_FREQS;
                    }
#else
                    if (unlikely(ScanmovUintCapped32(cap_div_10, cap_mod_10, &alt_freq_iter, &loaded_allele_idx))) {
                      goto ReadAlleleFreqs_ret_INVALID_FREQS;
                    }
#endif
                    if (unlikely(*alt_freq_iter != '=')) {
                      goto ReadAlleleFreqs_ret_INVALID_FREQS;
                    }
                    if (IsSet(matched_loaded_alleles, loaded_allele_idx)) {
                      const uint32_t internal_allele_idx = loaded_to_internal_allele_idx[loaded_allele_idx];
                      if (unlikely(cur_allele_freqs[internal_allele_idx] != 0.0)) {
                        snprintf(g_logbuf, kLogbufSize, "Error: Duplicate entry on line %" PRIuPTR " of --read-freq file.\n", line_idx);
                        goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
                      }
                      ++alt_freq_iter;
                      double dxx;
                      const char* cur_freq_end = ScantokDouble(alt_freq_iter, &dxx);
                      if (!cur_freq_end) {
                        if (likely(IsNanStr(alt_freq_iter, cur_entry_end - alt_freq_iter))) {
                          goto ReadAlleleFreqs_skip_variant;
                        }
                        goto ReadAlleleFreqs_ret_INVALID_FREQS;
                      }
                      if (unlikely((dxx < 0.0) || (dxx > freq_max))) {
                        goto ReadAlleleFreqs_ret_INVALID_FREQS;
                      }
                      if (!is_frac) {
                        dxx = ForceCountToDosage(dxx);
                      }
                      cur_allele_freqs[internal_allele_idx] = dxx;
                    }
                    alt_freq_iter = cur_entry_end;
                    if (alt_freq_iter == alt_freq_end) {
                      break;
                    }
                    ++alt_freq_iter;
                  }
                } else {
                  char* alt_freq_iter = alt_freq_start;
                  char* alt_freq_end = &(alt_freq_start[full_slen]);
                  while (1) {
                    char* cur_entry_end = AdvToDelim(alt_freq_iter, ',');
                    const uint32_t cur_entry_slen = cur_entry_end - alt_freq_iter;
                    char* eq_ptr = S_CAST(char*, memchr(alt_freq_iter, '=', cur_entry_slen));
                    if (unlikely(!eq_ptr)) {
                      goto ReadAlleleFreqs_ret_INVALID_FREQS;
                    }

                    // necessary for string comparison
                    *eq_ptr++ = '\0';

                    const uint32_t cur_blen = eq_ptr - alt_freq_iter;
                    // O(n^2), may want to replace with O(n log n)
                    for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
                      if (memequal(alt_freq_iter, cur_alleles[internal_allele_idx], cur_blen)) {
                        if (unlikely(cur_allele_freqs[internal_allele_idx] != 0.0)) {
                          snprintf(g_logbuf, kLogbufSize, "Error: Duplicate entry on line %" PRIuPTR " of --read-freq file.\n", line_idx);
                          goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
                        }
                        alt_freq_iter = eq_ptr;
                        double dxx;
                        const char* cur_freq_end = ScantokDouble(alt_freq_iter, &dxx);
                        if (!cur_freq_end) {
                          if (likely(IsNanStr(alt_freq_iter, cur_entry_end - alt_freq_iter))) {
                            goto ReadAlleleFreqs_skip_variant;
                          }
                          goto ReadAlleleFreqs_ret_INVALID_FREQS;
                        }
                        if (unlikely((dxx < 0.0) || (dxx > freq_max))) {
                          goto ReadAlleleFreqs_ret_INVALID_FREQS;
                        }
                        if (!is_frac) {
                          dxx = ForceCountToDosage(dxx);
                        }
                        cur_allele_freqs[internal_allele_idx] = dxx;
                        break;
                      }
                    }
                    alt_freq_iter = cur_entry_end;
                    if (alt_freq_iter == alt_freq_end) {
                      break;
                    }
                    ++alt_freq_iter;
                  }
                }
              }
            }
          }
        }
        // The IsSet condition is false when the following happens:
        // - --read-freq file has multiple ALT alleles.
        // - These ALT alleles collectively match *all* alleles in the main
        //   dataset's version of the variant, including the main dataset's
        //   REF.  The --read-freq file's REF allele is not in the main
        //   dataset at all.
        // - So, while we normally infer the REF frequency by subtracting the
        //   sum of the other frequencies from 1, in this subcase there is
        //   nothing to infer.
        if (infer_one_freq && IsSet(matched_loaded_alleles, infer_freq_loaded_idx)) {
          double adj_obs_ct_recip = 1.0;
          if (header_cols & kfReadFreqColsetObsCt) {
            uint32_t obs_ct_raw;
            if (unlikely(ScanUintCapped(token_ptrs[kReadFreqColObsCt], UINT32_MAX, &obs_ct_raw))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid allele count on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
            }
            if ((!obs_ct_raw) && (af_pseudocount == 0.0)) {
              goto ReadAlleleFreqs_skip_variant;
            }
            adj_obs_ct_recip = 1.0 / (u63tod(obs_ct_raw) + af_pseudocount * cur_allele_ct);
          }
          const uint32_t infer_freq_internal_idx = loaded_to_internal_allele_idx[infer_freq_loaded_idx];
          if (cur_allele_ct == 2) {
            // optimize common case
            double known_freq_d = cur_allele_freqs[1 - infer_freq_internal_idx] + af_pseudocount;
            double known_scaled_freq = known_freq_d * adj_obs_ct_recip;
            if (known_scaled_freq <= 1.0) {
              if (infer_freq_internal_idx) {
                allele_freqs_write[0] = known_scaled_freq;
              } else {
                allele_freqs_write[0] = 1.0 - known_scaled_freq;
              }
            } else if (likely(known_scaled_freq < ((1.0 + kSmallEpsilon) / 0.99))) {
              if (infer_freq_internal_idx) {
                allele_freqs_write[0] = 1.0;
              } else {
                allele_freqs_write[0] = 0.0;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Frequency/count too large on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
            }
          } else {
            if (af_pseudocount != 0.0) {
              for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
                cur_allele_freqs[internal_allele_idx] += af_pseudocount;
              }
            }
            cur_allele_freqs[infer_freq_internal_idx] = 0.0;
            double known_freq_sum_d = 0.0;
            for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
              known_freq_sum_d += cur_allele_freqs[internal_allele_idx];
            }
            double known_scaled_freq_sum = known_freq_sum_d * adj_obs_ct_recip;
            if (likely(known_scaled_freq_sum < ((1.0 + kSmallEpsilon) / 0.99))) {
              if (known_scaled_freq_sum > 1.0) {
                // possible rounding error, rescale
                adj_obs_ct_recip = 1.0 / known_scaled_freq_sum;
                known_scaled_freq_sum = 1.0;
              }
              const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
              for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct_m1; ++internal_allele_idx) {
                double dxx;
                if (internal_allele_idx == infer_freq_internal_idx) {
                  dxx = 1.0 - known_scaled_freq_sum;
                } else {
                  dxx = adj_obs_ct_recip * cur_allele_freqs[internal_allele_idx];
                }
                allele_freqs_write[internal_allele_idx] = dxx;
              }
            } else {
              snprintf(g_logbuf, kLogbufSize, "Error: Frequency/count too large on line %" PRIuPTR " of --read-freq file.\n", line_idx);
              goto ReadAlleleFreqs_ret_MALFORMED_INPUT_2;
            }
          }
        } else {
          // complete frequency or count data
          if (af_pseudocount != 0.0) {
            for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
              cur_allele_freqs[internal_allele_idx] += af_pseudocount;
            }
          }
          double tot_freq = 0.0;
          for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct; ++internal_allele_idx) {
            tot_freq += cur_allele_freqs[internal_allele_idx];
          }
          if (tot_freq == 0.0) {
            goto ReadAlleleFreqs_skip_variant;
          }
          const double tot_freq_recip = 1.0 / tot_freq;
          const uint32_t cur_allele_ct_m1 = cur_allele_ct - 1;
          for (uint32_t internal_allele_idx = 0; internal_allele_idx != cur_allele_ct_m1; ++internal_allele_idx) {
            allele_freqs_write[internal_allele_idx] = tot_freq_recip * cur_allele_freqs[internal_allele_idx];
          }
        }

        ++loaded_variant_ct;
        if (!(loaded_variant_ct % 10000)) {
          printf("\r--read-freq: Frequencies for %uk variants loaded.", loaded_variant_ct / 1000);
          fflush(stdout);
        }
      }
      while (0) {
      ReadAlleleFreqs_skip_variant:
        SetBit(variant_uidx, *variant_afreqcalcp);
        ++skipped_variant_ct;
      }
    }
    if (unlikely(TextStreamErrcode2(&read_freq_txs, &reterr))) {
      goto ReadAlleleFreqs_ret_TSTREAM_FAIL;
    }
    // Set variant_afreqcalc to the set of variants in variant_include, but
    // without loaded frequencies.
    BitvecInvertAndMask(variant_include, raw_variant_ctl, already_seen);
    // already_seen is now the set of remaining variants that *weren't* seen,
    // and variant_afreqcalc is the set of variants which were seen but
    // skipped due to incomplete information.
    BitvecOr(already_seen, raw_variant_ctl, *variant_afreqcalcp);
    putc_unlocked('\r', stdout);
    logprintf("--read-freq: Frequencies for %u variant%s loaded.\n", loaded_variant_ct, (loaded_variant_ct == 1)? "" : "s");
    if (skipped_variant_ct) {
      logerrprintfww("Warning: %" PRIuPTR " entr%s skipped due to missing variant IDs, mismatching allele codes, and/or zero observations.\n", skipped_variant_ct, (skipped_variant_ct == 1)? "y" : "ies");
    }
  }
  while (0) {
  ReadAlleleFreqs_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--read-freq file", &read_freq_txs);
    break;
  ReadAlleleFreqs_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ReadAlleleFreqs_ret_UNRECOGNIZED_HEADER:
    logerrputs("Error: Unrecognized header line in --read-freq file.\n");
    reterr = kPglRetMalformedInput;
    break;
  ReadAlleleFreqs_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of --read-freq file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  ReadAlleleFreqs_ret_INVALID_FREQS:
    snprintf(g_logbuf, kLogbufSize, "Error: Invalid frequencies/counts on line %" PRIuPTR " of --read-freq file.\n", line_idx);
  ReadAlleleFreqs_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
  ReadAlleleFreqs_ret_MALFORMED_INPUT_2:
    logputs("\n");
    logerrputsb();
  ReadAlleleFreqs_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 ReadAlleleFreqs_ret_1:
  CleanupTextStream2("--read-freq file", &read_freq_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

void ComputeMajAlleles(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const double* allele_freqs, uint32_t variant_ct, AlleleCode* maj_alleles) {
  uint32_t cur_allele_ct = 2;
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
    const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    uintptr_t allele_idx_base;
    if (!allele_idx_offsets) {
      allele_idx_base = variant_uidx;
    } else {
      allele_idx_base = allele_idx_offsets[variant_uidx];
      cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_base;
      allele_idx_base -= variant_uidx;
    }
    maj_alleles[variant_uidx] = GetMajIdx(&(allele_freqs[allele_idx_base]), cur_allele_ct);
  }
}

PglErr MindFilter(const uint32_t* sample_missing_cts, const uint32_t* sample_hethap_cts, const SampleIdInfo* siip, uint32_t raw_sample_ct, uint32_t variant_ct, uint32_t variant_ct_y, double mind_thresh, uintptr_t* sample_include, uintptr_t* sex_male, uint32_t* sample_ct_ptr, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto MindFilter_ret_1;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);

    uint32_t max_missing_cts[2];
    mind_thresh *= 1 + kSmallEpsilon;
    max_missing_cts[0] = S_CAST(int32_t, u31tod(variant_ct - variant_ct_y) * mind_thresh);
    max_missing_cts[1] = S_CAST(int32_t, u31tod(variant_ct) * mind_thresh);
    uintptr_t* newly_excluded;
    if (unlikely(bigstack_calloc_w(raw_sample_ctl, &newly_excluded))) {
      goto MindFilter_ret_NOMEM;
    }
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != orig_sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      uint32_t cur_missing_geno_ct = sample_missing_cts[sample_uidx];
      if (sample_hethap_cts) {
        cur_missing_geno_ct += sample_hethap_cts[sample_uidx];
      }
      if (cur_missing_geno_ct > max_missing_cts[IsSet(sex_male, sample_uidx)]) {
        SetBit(sample_uidx, newly_excluded);
      }
    }
    const uint32_t removed_ct = PopcountWords(newly_excluded, raw_sample_ctl);
    // don't bother with allow_no_samples check here, better to have that in
    // just one place
    logprintf("%u sample%s removed due to missing genotype data (--mind).\n", removed_ct, (removed_ct == 1)? "" : "s");
    if (removed_ct) {
      BitvecInvmask(newly_excluded, raw_sample_ctl, sample_include);
      BitvecInvmask(newly_excluded, raw_sample_ctl, sex_male);
      snprintf(outname_end, kMaxOutfnameExtBlen, ".mindrem.id");
      reterr = WriteSampleIds(newly_excluded, siip, outname, removed_ct);
      if (unlikely(reterr)) {
        goto MindFilter_ret_1;
      }
      logprintfww("ID%s written to %s .\n", (removed_ct == 1)? "" : "s", outname);
      *sample_ct_ptr -= removed_ct;
    }
  }
  while (0) {
  MindFilter_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 MindFilter_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr SelectSidRepresentatives(const uintptr_t* sex_nm, const uintptr_t* sex_male, const uint32_t* sample_missing_cts, const uint32_t* sample_hethap_cts, const PedigreeIdInfo* piip, SelectSidMissingnessMode missingness_mode, SelectSidTiebreakMode tiebreak_mode, uint32_t parents_only, uint32_t raw_sample_ct, uintptr_t* sample_include, uint32_t* sample_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t orig_sample_ct = *sample_ct_ptr;
    if (!orig_sample_ct) {
      goto SelectSidRepresentatives_ret_1;
    }
    const SampleIdInfo* siip = &(piip->sii);
    if (!siip->sids) {
      logputs("Skipping --select-sid-representatives: all FID+IID pairs are unique.\n");
      goto SelectSidRepresentatives_ret_1;
    }
    const uint32_t use_nsort = (tiebreak_mode == kSelectSidTiebreakFirst) || (tiebreak_mode == kSelectSidTiebreakLast);
    char* sorted_3idbox;
    uint32_t* xid_map;
    uintptr_t max_3id_blen;
    reterr = SortedXidboxInitAlloc(sample_include, siip, orig_sample_ct, kfXidModeFidIidSid, use_nsort, &sorted_3idbox, &xid_map, &max_3id_blen);
    if (unlikely(reterr)) {
      goto SelectSidRepresentatives_ret_1;
    }
    char* idbuf;
    if (unlikely(bigstack_alloc_c(siip->max_sample_id_blen + 1, &idbuf))) {
      goto SelectSidRepresentatives_ret_NOMEM;
    }
    uintptr_t* all_parents = nullptr;
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    if (parents_only) {
      if (unlikely(bigstack_calloc_w(raw_sample_ctl, &all_parents))) {
        goto SelectSidRepresentatives_ret_NOMEM;
      }
      MarkParents(&(piip->parental_id_info), sorted_3idbox, xid_map, orig_sample_ct, max_3id_blen, use_nsort, all_parents, idbuf);
      if (AllWordsAreZero(all_parents, raw_sample_ctl)) {
        logputs("Skipping --select-sid-representatives parents-only: no parents remaining.\n");
        goto SelectSidRepresentatives_ret_1;
      }
    }
    if (missingness_mode == kSelectSidMissingness0) {
      sample_missing_cts = nullptr;
    }
    const int32_t tiebreak_dir = (tiebreak_mode >= kSelectSidTiebreakLast)? -1 : 1;
    uint32_t xid_idx_end = 0;
    do {
      const uint32_t xid_idx_start = xid_idx_end;
      const char* fid_start = &(sorted_3idbox[xid_idx_start * max_3id_blen]);
      const char* iid_end = AdvToDelim(AdvPastDelim(fid_start, '\t'), '\t');
      const uint32_t sample_2id_slen = iid_end - fid_start;
      memcpy(idbuf, fid_start, sample_2id_slen);
      idbuf[sample_2id_slen] = ' ';
      const uint32_t sample_2id_blen = sample_2id_slen + 1;
      idbuf[sample_2id_blen] = '\0';
      if (use_nsort) {
        xid_idx_end = ExpsearchNsortStrLb(idbuf, sorted_3idbox, max_3id_blen, orig_sample_ct, xid_idx_start + 1);
      } else {
        xid_idx_end = ExpsearchStrLb(idbuf, sorted_3idbox, sample_2id_blen, max_3id_blen, orig_sample_ct, xid_idx_start + 1);
      }
      if (xid_idx_end == xid_idx_start + 1) {
        continue;
      }
      uint32_t sample_uidx = xid_map[xid_idx_start];
      if (all_parents && (!IsSet(all_parents, sample_uidx))) {
        continue;
      }
      // intentional underflow on xid_idx_start == 0
      uint32_t xid_idx = (tiebreak_dir == 1)? (xid_idx_start - 1) : xid_idx_end;
      const uint32_t xid_idx_stop = (tiebreak_dir == 1)? (xid_idx_end - 1) : xid_idx_start;
      const uint32_t sex_nm_plus_male = IsSet(sex_nm, sample_uidx) + IsSet(sex_male, sample_uidx);
      uint32_t best_missing_ct = UINT32_MAX;
      // keep correct sample in sid-only case
      uint32_t cur_missing_ct = UINT32_MAXM1;
      uint32_t best_sample_uidx = 0;
      do {
        xid_idx += tiebreak_dir;
        sample_uidx = xid_map[xid_idx];
        if (unlikely(sex_nm_plus_male != IsSet(sex_nm, sample_uidx) + IsSet(sex_male, sample_uidx))) {
          idbuf[sample_2id_slen] = '\0';
          TabsToSpaces(idbuf);
          logerrprintf("Error: Inconsistent sex information for sample \"%s\".\n", idbuf);
          goto SelectSidRepresentatives_ret_INCONSISTENT_INPUT;
        }
        ClearBit(sample_uidx, sample_include);
        if (sample_missing_cts) {
          cur_missing_ct = sample_missing_cts[sample_uidx];
          if (sample_hethap_cts) {
            cur_missing_ct += sample_hethap_cts[sample_uidx];
          }
        }
        if (cur_missing_ct >= best_missing_ct) {
          continue;
        }
        best_missing_ct = cur_missing_ct;
        best_sample_uidx = sample_uidx;
      } while (xid_idx != xid_idx_stop);
      SetBit(best_sample_uidx, sample_include);
    } while (xid_idx_end != orig_sample_ct);
    const uint32_t sample_ct = PopcountWords(sample_include, raw_sample_ctl);
    const uint32_t removed_ct = orig_sample_ct - sample_ct;
    logprintf("--select-sid-representatives: %u sample%s removed.\n", removed_ct, (removed_ct == 1)? "" : "s");
    *sample_ct_ptr = sample_ct;
  }
  while (0) {
  SelectSidRepresentatives_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SelectSidRepresentatives_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 SelectSidRepresentatives_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

void EnforceGenoThresh(const ChrInfo* cip, const uint32_t* variant_missing_cts, const uint32_t* variant_hethap_cts, uint32_t sample_ct, uint32_t male_ct, uint32_t first_hap_uidx, double geno_thresh, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t prefilter_variant_ct = *variant_ct_ptr;
  geno_thresh *= 1 + kSmallEpsilon;
  const uint32_t missing_max_ct_nony = S_CAST(int32_t, geno_thresh * u31tod(sample_ct));
  const uint32_t missing_max_ct_y = S_CAST(int32_t, geno_thresh * u31tod(male_ct));
  uint32_t cur_missing_max_ct = missing_max_ct_nony;
  uint32_t removed_ct = 0;
  uint32_t y_thresh = UINT32_MAX;
  uint32_t y_end = UINT32_MAX;
  uint32_t y_code;
  if (XymtExists(cip, kChrOffsetY, &y_code)) {
    const uint32_t y_chr_fo_idx = cip->chr_idx_to_foidx[y_code];
    y_thresh = cip->chr_fo_vidx_start[y_chr_fo_idx];
    y_end = cip->chr_fo_vidx_start[y_chr_fo_idx + 1];
  }
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  for (uint32_t variant_idx = 0; variant_idx != prefilter_variant_ct; ++variant_idx) {
    const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    if (variant_uidx >= y_thresh) {
      if (variant_uidx < y_end) {
        y_thresh = y_end;
        cur_missing_max_ct = missing_max_ct_y;
      } else {
        y_thresh = UINT32_MAX;
        cur_missing_max_ct = missing_max_ct_nony;
      }
    }
    uint32_t cur_missing_ct = variant_missing_cts[variant_uidx];
    if (variant_uidx >= first_hap_uidx) {
      cur_missing_ct += variant_hethap_cts[variant_uidx - first_hap_uidx];
    }
    if (cur_missing_ct > cur_missing_max_ct) {
      ClearBit(variant_uidx, variant_include);
      ++removed_ct;
    }
  }
  logprintf("--geno: %u variant%s removed due to missing genotype data.\n", removed_ct, (removed_ct == 1)? "" : "s");
  *variant_ct_ptr -= removed_ct;
}

PglErr EnforceHweThresh(const ChrInfo* cip, const uintptr_t* allele_idx_offsets, const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_raw_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, autosomal_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_male_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 3, founder_x_nosex_geno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_knownsex_xgeno_cts), const STD_ARRAY_PTR_DECL(uint32_t, 2, x_male_xgeno_cts), const double* hwe_x_ln_pvals, MiscFlags misc_flags, double hwe_ln_thresh_base, double hwe_sample_size_term, uint32_t nonfounders, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  uint32_t prefilter_variant_ct = *variant_ct_ptr;
  const uint32_t midp = (misc_flags / kfMiscHweMidp) & 1;
  const uint32_t keep_fewhet = (misc_flags / kfMiscHweKeepFewhet) & 1;
  hwe_ln_thresh_base *= 1 + kSmallEpsilon;
  const double hwe_sample_size_coeff = (hwe_sample_size_term <= 0)? 0 : (hwe_sample_size_term * (-kLn10));
  uint32_t removed_ct = 0;
  uint32_t min_obs = UINT32_MAX;
  uint32_t max_obs = 0;
  uint32_t male_a1_ct = 0;
  uint32_t male_ax_ct = 0;
  const double* hwe_x_ln_pvals_iter = hwe_x_ln_pvals;
  uintptr_t autosomal_xgeno_idx = 0;
  uintptr_t x_xgeno_idx = 0;
  uint32_t x_skip_code = UINT32_MAX;
  uint32_t x_start = 0;
  uint32_t x_code;
  if (XymtExists(cip, kChrOffsetX, &x_code)) {
    const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
    x_start = cip->chr_fo_vidx_start[x_chr_fo_idx];
    if (!hwe_x_ln_pvals) {
      x_skip_code = x_code;  // only set this if we're skipping chrX
      prefilter_variant_ct -= PopcountBitRange(variant_include, x_start, cip->chr_fo_vidx_start[x_chr_fo_idx + 1]);
    }
  }
  // bugfix (2 Apr 2019): wasn't skipping chrY/chrM properly
  uint32_t y_code;
  if (XymtExists(cip, kChrOffsetY, &y_code)) {
    prefilter_variant_ct -= CountChrVariantsUnsafe(variant_include, cip, y_code);
  }
  uint32_t mt_code;
  if (XymtExists(cip, kChrOffsetMT, &mt_code)) {
    prefilter_variant_ct -= CountChrVariantsUnsafe(variant_include, cip, mt_code);
  }
  uint32_t chr_fo_idx = UINT32_MAX;
  uint32_t chr_end = 0;
  uint32_t is_x = 0;
  uintptr_t variant_uidx_base = 0;
  uintptr_t variant_include_bits = variant_include[0];
  uint32_t xallele_ct = 0;
  for (uint32_t variant_idx = 0; variant_idx != prefilter_variant_ct; ++variant_idx) {
    uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
    if (variant_uidx >= chr_end) {
      uint32_t chr_idx;
      do {
        ++chr_fo_idx;
        chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        chr_idx = cip->chr_file_order[chr_fo_idx];
      } while ((variant_uidx >= chr_end) || (chr_idx == y_code) || (chr_idx == mt_code) || (chr_idx == x_skip_code));
      BitIter1Start(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], &variant_uidx_base, &variant_include_bits);
      variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      is_x = (chr_idx == x_code);
    }
    STD_ARRAY_KREF(uint32_t, 3) cur_geno_cts = founder_raw_geno_cts[variant_uidx];
    if (allele_idx_offsets) {
      const uint32_t allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx];
      xallele_ct = (allele_ct == 2)? 0 : (allele_ct - 1);
    }
    uint32_t hom_a1_ct = cur_geno_cts[0];
    uint32_t het_a1_ct = cur_geno_cts[1];
    uint32_t two_ax_ct = cur_geno_cts[2];
    uint32_t test_failed = 0;
    uint32_t pval_computed = 0;
    uint32_t cur_obs_ct;
    if (!is_x) {
      cur_obs_ct = hom_a1_ct + het_a1_ct + two_ax_ct;
      if (!cur_obs_ct) {
        if (xallele_ct) {
          autosomal_xgeno_idx += xallele_ct;
        }
        continue;
      }
      const double hwe_ln_thresh = hwe_ln_thresh_base + u31tod(cur_obs_ct) * hwe_sample_size_coeff;
      const double hwe_thresh = exp(hwe_ln_thresh);
      for (uint32_t xallele_idx = 0; ; ++xallele_idx) {
        if (keep_fewhet) {
          if (het_a1_ct * S_CAST(uint64_t, het_a1_ct) <= (4LLU * hom_a1_ct) * two_ax_ct) {
            goto EnforceHweThresh_skip_autosomal;
          }
        }
        pval_computed = 1;
        test_failed = HweThreshLn(het_a1_ct, hom_a1_ct, two_ax_ct, midp, hwe_thresh, hwe_ln_thresh);
        if (test_failed) {
          break;
        }
      EnforceHweThresh_skip_autosomal:
        if (xallele_idx == xallele_ct) {
          break;
        }
        STD_ARRAY_KREF(uint32_t, 2) cur_xgeno_cts = autosomal_xgeno_cts[autosomal_xgeno_idx + xallele_idx];
        hom_a1_ct = cur_xgeno_cts[0];
        het_a1_ct = cur_xgeno_cts[1];
        two_ax_ct = cur_obs_ct - hom_a1_ct - het_a1_ct;
      }
      autosomal_xgeno_idx += xallele_ct;
    } else {
      if (founder_x_male_geno_cts) {
        STD_ARRAY_KREF(uint32_t, 3) cur_male_geno_cts = founder_x_male_geno_cts[variant_uidx - x_start];
        male_a1_ct = cur_male_geno_cts[0];
        hom_a1_ct -= male_a1_ct;
        het_a1_ct -= cur_male_geno_cts[1];
        male_ax_ct = cur_male_geno_cts[2];
        two_ax_ct -= male_ax_ct;
      }
      if (founder_x_nosex_geno_cts) {
        STD_ARRAY_KREF(uint32_t, 3) cur_nosex_geno_cts = founder_x_nosex_geno_cts[variant_uidx - x_start];
        hom_a1_ct -= cur_nosex_geno_cts[0];
        het_a1_ct -= cur_nosex_geno_cts[1];
        two_ax_ct -= cur_nosex_geno_cts[2];
      }
      const uint32_t female_obs_ct = hom_a1_ct + het_a1_ct + two_ax_ct;
      cur_obs_ct = female_obs_ct + male_a1_ct + male_ax_ct;
      const double hwe_ln_thresh = hwe_ln_thresh_base + u31tod(cur_obs_ct) * hwe_sample_size_coeff;
      pval_computed = 1;
      double joint_ln_pval = *hwe_x_ln_pvals_iter++;
      for (uint32_t xallele_idx = 0; ; ++xallele_idx) {
        test_failed = (joint_ln_pval < hwe_ln_thresh);
        if (test_failed && keep_fewhet && (het_a1_ct * S_CAST(uint64_t, het_a1_ct) < (4LLU * hom_a1_ct) * two_ax_ct)) {
          // In chrX keep_fewhet case, we only keep a
          // Graffelman/Weir-test-failing variant if there are few female hets,
          // *and* the male/female allele-frequency imbalance isn't severe
          // enough to make the G/W test fail on its own.
          joint_ln_pval = joint_ln_pval - hwe_ln_thresh;
          test_failed = !HweThreshLn(het_a1_ct, hom_a1_ct, two_ax_ct, midp, exp(joint_ln_pval), joint_ln_pval);
        }
        // bugfix (27 Jun 2020): don't clobber previous allele-test failure if
        // variant is multiallelic
        if (test_failed || (xallele_idx == xallele_ct)) {
          break;
        }
        STD_ARRAY_KREF(uint32_t, 2) cur_xgeno_cts = x_knownsex_xgeno_cts[x_xgeno_idx + xallele_idx];
        hom_a1_ct = cur_xgeno_cts[0];
        het_a1_ct = cur_xgeno_cts[1];
        if (x_male_xgeno_cts) {
          STD_ARRAY_KREF(uint32_t, 2) cur_male_xgeno_cts = x_male_xgeno_cts[x_xgeno_idx + xallele_idx];
          hom_a1_ct -= cur_male_xgeno_cts[0];
          het_a1_ct -= cur_male_xgeno_cts[1];
        }
        two_ax_ct = female_obs_ct - hom_a1_ct - het_a1_ct;
        joint_ln_pval = hwe_x_ln_pvals_iter[xallele_idx];
      }
      x_xgeno_idx += xallele_ct;
      hwe_x_ln_pvals_iter = &(hwe_x_ln_pvals_iter[xallele_ct]);
    }
    if (test_failed) {
      ClearBit(variant_uidx, variant_include);
      ++removed_ct;
    }
    if (pval_computed) {
      if (cur_obs_ct < min_obs) {
        min_obs = cur_obs_ct;
      }
      if (cur_obs_ct > max_obs) {
        max_obs = cur_obs_ct;
      }
    }
  }
  // Greer et al. "A reassessment of Hardy-Weinberg equilibrium filtering in
  // large sample Genomic studies" preprint recommends
  // hwe_sample_size_term=0.001.  If neither hwe_sample_size_term nor
  // 'keep-fewhet' were specified, warn if the filter is effectively more than
  // twice as stringent.
  if ((hwe_sample_size_term == -1) && (!keep_fewhet) && (u31tod(max_obs) * 0.0005 * (-kLn10) < hwe_ln_thresh_base)) {
    logerrputs("Error: --hwe filter is suspiciously strict for the sample size; you may be\nfiltering out many variants with association signals.\n");
    logerrprintfww("If you are performing routine QC, consider using --hwe with a sample-size term (e.g. \"--hwe 1e-5 0.001%s\" results in a threshold of 1e-6 with 1000 observations, 1e-7 with 2000 observations, etc.), and/or specifying 'keep-fewhet' to remove only variants with heterozygosity excess.\n", midp? " midp" : "");
    logerrputs("Alternatively, if the strictness is intentional (e.g. you are preparing data\nfor a method that requires variants to be near Hardy-Weinberg equilibrium),\nexplicitly specify a sample-size term of 0.\n");
    return kPglRetInconsistentInput;
  } else if ((hwe_sample_size_term <= 0) && (S_CAST(uint64_t, max_obs) * 9 > S_CAST(uint64_t, min_obs) * 10)) {
    logerrprintfww("Warning: --hwe observation counts vary by more than 10%. Consider using --hwe with a sample-size term (e.g. \"--hwe 1e-5 0.001%s%s\" results in a threshold of 1e-6 with 1000 observations, 1e-7 with 2000 observations, etc.), and/or filtering out high-missingness variants with --geno.\n", midp? " midp" : "", keep_fewhet? " keep-fewhet" : "");
  }
  logprintfww("--hwe%s%s: %u variant%s removed due to Hardy-Weinberg exact test (%s).\n", midp? " midp" : "", keep_fewhet? " keep-fewhet" : "", removed_ct, (removed_ct == 1)? "" : "s", nonfounders? "all samples" : "founders only");
  *variant_ct_ptr -= removed_ct;
  return kPglRetSuccess;
}

double GetTypedFreq(const double* cur_allele_freqs, uint32_t allele_ct, FreqFilterMode mode) {
  if ((allele_ct == 2) || (mode == kFreqFilterNref)) {
    const double ref_freq = cur_allele_freqs[0];
    const double nonref_freq = 1.0 - ref_freq;
    if (mode & (kFreqFilterNref | kFreqFilterAlt1)) {
      return nonref_freq;
    }
    return MINV(nonref_freq, ref_freq);
  }
  if (mode == kFreqFilterAlt1) {
    return cur_allele_freqs[1];
  }
  double tot_nonlast_freq = cur_allele_freqs[0];
  const uint32_t allele_ct_m1 = allele_ct - 1;
  if (mode == kFreqFilterNonmajor) {
    double max_freq = tot_nonlast_freq;
    for (uint32_t allele_idx = 1; allele_idx != allele_ct_m1; ++allele_idx) {
      const double cur_alt_freq = cur_allele_freqs[allele_idx];
      tot_nonlast_freq += cur_alt_freq;
      if (cur_alt_freq > max_freq) {
        max_freq = cur_alt_freq;
      }
    }
    const double nonmajor_freq = 1.0 - max_freq;
    return MINV(nonmajor_freq, tot_nonlast_freq);
  }
  // minor
  double min_freq = tot_nonlast_freq;
  for (uint32_t allele_idx = 1; allele_idx != allele_ct_m1; ++allele_idx) {
    const double cur_alt_freq = cur_allele_freqs[allele_idx];
    tot_nonlast_freq += cur_alt_freq;
    if (cur_alt_freq < min_freq) {
      min_freq = cur_alt_freq;
    }
  }
  const double last_freq = MAXV(0.0, 1.0 - tot_nonlast_freq);
  return MINV(min_freq, last_freq);
}

uint64_t GetTypedDdosage(const uint64_t* cur_allele_ddosages, uint32_t allele_ct, FreqFilterMode mode) {
  if (mode == kFreqFilterAlt1) {
    return cur_allele_ddosages[1];
  }
  if (mode == kFreqFilterNref) {
    uint64_t nref_ddosage = cur_allele_ddosages[1];
    for (uint32_t allele_idx = 2; allele_idx != allele_ct; ++allele_idx) {
      nref_ddosage += cur_allele_ddosages[allele_idx];
    }
    return nref_ddosage;
  }
  if (allele_ct == 2) {
    return MINV(cur_allele_ddosages[0], cur_allele_ddosages[1]);
  }
  if (mode == kFreqFilterNonmajor) {
    uint64_t max_ddosage = cur_allele_ddosages[0];
    uint64_t tot_ddosage = max_ddosage;
    for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
      const uint64_t cur_ddosage = cur_allele_ddosages[allele_idx];
      tot_ddosage += cur_ddosage;
      if (cur_ddosage > max_ddosage) {
        max_ddosage = cur_ddosage;
      }
    }
    return (tot_ddosage - max_ddosage);
  }
  // minor
  uint64_t min_ddosage = cur_allele_ddosages[0];
  for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
    const uint64_t cur_ddosage = cur_allele_ddosages[allele_idx];
    if (cur_ddosage < min_ddosage) {
      min_ddosage = cur_ddosage;
    }
  }
  return min_ddosage;
}

void EnforceFreqConstraints(const uintptr_t* allele_idx_offsets, const uint64_t* founder_allele_ddosages, const double* allele_freqs, STD_ARRAY_KREF(FreqFilterMode, 4) filter_modes, double min_maf, double max_maf, uint64_t min_allele_ddosage, uint64_t max_allele_ddosage, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t prefilter_variant_ct = *variant_ct_ptr;
  uint32_t removed_ct = 0;
  if ((min_maf != 0.0) || (max_maf != 1.0)) {
    // defend against floating point error
    min_maf *= 1.0 - kSmallEpsilon;
    max_maf *= 1.0 + kSmallEpsilon;
  } else {
    allele_freqs = nullptr;
  }
  const uint32_t dosage_filter = min_allele_ddosage || (max_allele_ddosage != (~0LLU));

  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  uint32_t allele_ct = 2;
  for (uint32_t variant_idx = 0; variant_idx != prefilter_variant_ct; ++variant_idx) {
    const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    uintptr_t allele_idx_offset_base;
    if (!allele_idx_offsets) {
      allele_idx_offset_base = 2 * variant_uidx;
    } else {
      allele_idx_offset_base = allele_idx_offsets[variant_uidx];
      allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
    }
    if (allele_freqs) {
      const double* cur_allele_freqs = &(allele_freqs[allele_idx_offset_base - variant_uidx]);
      if (min_maf != 0.0) {
        const double cur_typed_freq = GetTypedFreq(cur_allele_freqs, allele_ct, filter_modes[0]);
        if (cur_typed_freq < min_maf) {
          goto EnforceFreqConstraints_remove;
        }
      }
      if (max_maf < 1.0) {
        // technically a bit inefficient
        const double cur_typed_freq = GetTypedFreq(cur_allele_freqs, allele_ct, filter_modes[1]);
        if (cur_typed_freq > max_maf) {
          goto EnforceFreqConstraints_remove;
        }
      }
    }
    if (dosage_filter) {
      const uint64_t* cur_founder_allele_ddosages = &(founder_allele_ddosages[allele_idx_offset_base]);
      if (min_allele_ddosage) {
        const uint64_t cur_typed_ddosage = GetTypedDdosage(cur_founder_allele_ddosages, allele_ct, filter_modes[2]);
        if (cur_typed_ddosage < min_allele_ddosage) {
          goto EnforceFreqConstraints_remove;
        }
      }
      if (max_allele_ddosage != (~0LLU)) {
        const uint64_t cur_typed_ddosage = GetTypedDdosage(cur_founder_allele_ddosages, allele_ct, filter_modes[3]);
        if (cur_typed_ddosage > max_allele_ddosage) {
          goto EnforceFreqConstraints_remove;
        }
      }
    }
    continue;
  EnforceFreqConstraints_remove:
    ClearBit(variant_uidx, variant_include);
    ++removed_ct;
  }
  logprintfww("%u variant%s removed due to allele frequency threshold(s) (--maf/--max-maf/--mac/--max-mac).\n", removed_ct, (removed_ct == 1)? "" : "s");
  *variant_ct_ptr -= removed_ct;
}

void EnforceImpR2Thresh(const ChrInfo* cip, const double* imp_r2_vals, double imp_r2_min, double imp_r2_max, uint32_t is_minimac3_r2, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t prefilter_variant_ct = *variant_ct_ptr;
  imp_r2_min *= 1 - k2m35;
  imp_r2_max *= 1 + k2m35;
  uint32_t removed_ct = 0;
  uint32_t relevant_variant_ct = prefilter_variant_ct;
  // skip X, MT
  uint32_t x_code;
  if (XymtExists(cip, kChrOffsetX, &x_code)) {
    const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[x_code];
    relevant_variant_ct -= PopcountBitRange(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1]);
  }
  uint32_t mt_code;
  if (XymtExists(cip, kChrOffsetMT, &mt_code)) {
    const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[mt_code];
    relevant_variant_ct -= PopcountBitRange(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1]);
  }
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  uint32_t chr_fo_idx = UINT32_MAX;
  uint32_t chr_end = 0;
  for (uint32_t variant_idx = 0; variant_idx != relevant_variant_ct; ++variant_idx) {
    uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    if (variant_uidx >= chr_end) {
      do {
        uint32_t chr_idx;
        do {
          chr_idx = cip->chr_file_order[++chr_fo_idx];
          // bugfix (9 Jul 2018): this boolean condition was reversed
        } while ((chr_idx == x_code) || (chr_idx == mt_code));
        chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        variant_uidx = AdvBoundedTo1Bit(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], chr_end);
      } while (variant_uidx >= chr_end);
      const uint32_t variant_widx = variant_uidx / kBitsPerWord;
      variant_uidx_base = variant_widx * kBitsPerWord;
      cur_bits = variant_include[variant_widx] & ((~k1LU) << (variant_uidx % kBitsPerWord));
    }
    const double cur_imp_r2 = imp_r2_vals[variant_uidx];
    // er, we *don't* want to filter out NaN here, those variants may be worth
    // filtering but not on imputation quality grounds
    if ((cur_imp_r2 < imp_r2_min) || (cur_imp_r2 > imp_r2_max)) {
      ClearBit(variant_uidx, variant_include);
      ++removed_ct;
    }
  }
  logprintf("--%s-r2-filter: %u variant%s removed.\n", is_minimac3_r2? "minimac3" : "mach", removed_ct, (removed_ct == 1)? "" : "s");
  *variant_ct_ptr -= removed_ct;
}

void EnforceMinBpSpace(const ChrInfo* cip, const uint32_t* variant_bps, uint32_t min_bp_space, uintptr_t* variant_include, uint32_t* variant_ct_ptr) {
  const uint32_t orig_variant_ct = *variant_ct_ptr;
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = variant_include[0];
  uint32_t chr_fo_idx_p1 = 0;
  uint32_t chr_end = 0;
  uint32_t last_bp = 0;
  uint32_t removed_ct = 0;
  for (uint32_t variant_idx = 0; variant_idx != orig_variant_ct; ++variant_idx) {
    const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
    const uint32_t cur_bp = variant_bps[variant_uidx];
    if (variant_uidx >= chr_end) {
      do {
        chr_end = cip->chr_fo_vidx_start[++chr_fo_idx_p1];
      } while (variant_uidx >= chr_end);
      last_bp = cur_bp;
    } else {
      if (cur_bp < last_bp + min_bp_space) {
        ClearBit(variant_uidx, variant_include);
        ++removed_ct;
      } else {
        last_bp = cur_bp;
      }
    }
  }
  const uint32_t new_variant_ct = orig_variant_ct - removed_ct;
  logprintf("--bp-space: %u variant%s removed (%u remaining).\n", removed_ct, (removed_ct == 1)? "" : "s", new_variant_ct);
  *variant_ct_ptr = new_variant_ct;
}

// Generalization of plink 1.9 load_ax_alleles().
//
// Note that, when a variant has 2 (or more) nonmissing allele codes, we print
// a warning if the loaded allele doesn't match one of them (and the variant is
// unchanged).  However, if the variant is biallelic with one or two missing
// allele codes, a missing allele code is filled in.  This maintains
// compatibility with plink 1.9 --a1-allele/--a2-allele's behavior, and is good
// enough for most real-world use cases; however, it's a bit annoying to be
// unable to e.g. add the ref allele when dealing with a single-sample file
// with a het alt1/alt2 call.
// (Workaround for that case when merge is implemented: generate a
// single-sample file with all the right reference alleles, and merge with
// that.)
PglErr SetRefalt1FromFile(const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const TwoColParams* allele_flag_info, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_ct, Refalt1Mode refalt1_mode, uint32_t force, char input_missing_geno_char, uint32_t max_thread_ct, const char** allele_storage, uint32_t* max_allele_slen_ptr, AlleleCode* allele_permute, uintptr_t* nonref_flags, uintptr_t* previously_seen) {
  // temporary allocations on bottom, "permanent" allocations on top (so we
  // don't reset g_bigstack_end).
  // previously_seen[] should be preallocated iff both --ref-allele and
  // --alt[1]-allele are present in the same run.  when it is, this errors out
  // when the flags produce conflicting results.
  unsigned char* bigstack_mark = g_bigstack_base;

  // unaligned in middle of loop
  unsigned char* bigstack_end = g_bigstack_end;

  const uint32_t is_alt_or_alt1 = (refalt1_mode != kRefalt1ModeRef);
  const uint32_t is_alt1 = (refalt1_mode == kRefalt1ModeAlt1);
  const char* flagstr;
  if (refalt1_mode == kRefalt1ModeRef) {
    flagstr = "--ref-allele";
  } else if (refalt1_mode == kRefalt1ModeAlt) {
    flagstr = "--alt-allele";
  } else {
    flagstr = "--alt1-allele";
  }
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    const uint32_t max_allele_ctl = BitCtToWordCt(max_allele_ct);
    uintptr_t* already_seen;
    AlleleCode* loaded_ac_order;
    uintptr_t* alleles_seen;
    if (unlikely(bigstack_calloc_w(raw_variant_ctl, &already_seen) ||
                 bigstack_alloc_ac(max_allele_ct, &loaded_ac_order) ||
                 bigstack_alloc_w(max_allele_ctl, &alleles_seen))) {
      goto SetRefalt1FromFile_ret_NOMEM;
    }
    reterr = SizeAndInitTextStream(allele_flag_info->fname, bigstack_left() / 4, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto SetRefalt1FromFile_ret_TSTREAM_FAIL;
    }
    const uint32_t skip_ct = allele_flag_info->skip_ct;
    reterr = TextSkip(skip_ct, &txs);
    if (unlikely(reterr)) {
      if (reterr == kPglRetEof) {
        snprintf(g_logbuf, kLogbufSize, "Error: Fewer lines than expected in %s.\n", allele_flag_info->fname);
        goto SetRefalt1FromFile_ret_INCONSISTENT_INPUT_WW;
      }
      goto SetRefalt1FromFile_ret_TSTREAM_FAIL;
    }
    uint32_t* variant_id_htable = nullptr;
    uint32_t variant_id_htable_size;
    reterr = AllocAndPopulateIdHtableMt(variant_include, variant_ids, variant_ct, bigstack_left() / 8, max_thread_ct, &variant_id_htable, nullptr, &variant_id_htable_size, nullptr);
    if (unlikely(reterr)) {
      goto SetRefalt1FromFile_ret_1;
    }
    unsigned char* main_bigstack_base = g_bigstack_base;

    const uint32_t colid_first = (allele_flag_info->colid < allele_flag_info->colx);
    const char skipchar = allele_flag_info->skipchar;
    uint32_t colmin;
    uint32_t coldiff;
    if (colid_first) {
      colmin = allele_flag_info->colid - 1;
      coldiff = allele_flag_info->colx - allele_flag_info->colid;
    } else {
      colmin = allele_flag_info->colx - 1;
      coldiff = allele_flag_info->colid - allele_flag_info->colx;
    }
    line_idx = allele_flag_info->skip_ct;
    uintptr_t skipped_variant_ct = 0;
    uintptr_t missing_allele_ct = 0;
    uint32_t allele_mismatch_warning_ct = 0;
    uint32_t permuted_variant_ct = 0;
    uint32_t fillin_variant_ct = 0;
    uint32_t max_allele_blen = 1 + (*max_allele_slen_ptr);
    uint32_t cur_allele_ct = 2;
    while (1) {
      ++line_idx;
      char* line_start = TextGet(&txs);
      if (!line_start) {
        if (likely(!TextStreamErrcode2(&txs, &reterr))) {
          break;
        }
        goto SetRefalt1FromFile_ret_TSTREAM_FAIL;
      }
      char cc = *line_start;
      if (cc == skipchar) {
        continue;
      }
      // er, replace this with TokenLex0()...
      char* variant_id_start;
      char* allele_start;
      if (colid_first) {
        variant_id_start = NextTokenMult0(line_start, colmin);
        allele_start = NextTokenMult(variant_id_start, coldiff);
        if (unlikely(!allele_start)) {
          goto SetRefalt1FromFile_ret_MISSING_TOKENS;
        }
      } else {
        allele_start = NextTokenMult0(line_start, colmin);
        variant_id_start = NextTokenMult(allele_start, coldiff);
        if (unlikely(!variant_id_start)) {
          goto SetRefalt1FromFile_ret_MISSING_TOKENS;
        }
      }
      char* token_end = CurTokenEnd(variant_id_start);
      const uint32_t variant_id_slen = token_end - variant_id_start;
      const uint32_t variant_uidx = VariantIdDupflagHtableFind(variant_id_start, variant_ids, variant_id_htable, variant_id_slen, variant_id_htable_size, max_variant_id_slen);
      if (variant_uidx >> 31) {
        if (likely(variant_uidx == UINT32_MAX)) {
          ++skipped_variant_ct;
          continue;
        }
        snprintf(g_logbuf, kLogbufSize, "Error: %s variant ID '%s' appears multiple times in main dataset.\n", flagstr, variant_ids[variant_uidx & 0x7fffffff]);
        goto SetRefalt1FromFile_ret_MALFORMED_INPUT_WW;
      }
      token_end = CurTokenEnd(allele_start);
      const uint32_t allele_col_slen = token_end - allele_start;
      if (allele_col_slen == 1) {
        const char allele_char = *allele_start;
        // don't overwrite anything with missing code
        if ((allele_char == '.') || (allele_char == input_missing_geno_char)) {
          ++missing_allele_ct;
          continue;
        }
      }
      if (unlikely(IsSet(already_seen, variant_uidx))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Duplicate variant ID '%s' in %s file.\n", variant_ids[variant_uidx], flagstr);
        goto SetRefalt1FromFile_ret_MALFORMED_INPUT_WW;
      }
      SetBit(variant_uidx, already_seen);
      if (unlikely((refalt1_mode != kRefalt1ModeAlt) && memchr(allele_start, ',', allele_col_slen))) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s: Multiple alleles on line %" PRIuPTR ".%s\n", flagstr, line_idx, is_alt1? " (Did you mean --alt-allele?)" : "");
        goto SetRefalt1FromFile_ret_MALFORMED_INPUT_WW;
      }
      uintptr_t allele_idx_offset_base;
      if (!allele_idx_offsets) {
        allele_idx_offset_base = variant_uidx * 2;
      } else {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char** cur_alleles = &(allele_storage[allele_idx_offset_base]);
      ZeroWArr(max_allele_ctl, alleles_seen);
      // --ref-allele and --alt1-allele: There must be exactly one allele here.
      // --alt-allele: There can be anywhere from 1 to (cur_allele_ct - 1)
      //               alleles here.  If less than (cur_allele_ct - 1),
      //               remaining alleles retain original order.
      //               TODO: For cur_allele_ct > 20 or so, construct hash table
      //               instead of performing O(n^2) comparisons.
      uint32_t loaded_allele_ct = 0;
      uint32_t last_allele_blen;
      for (char* allele_iter = allele_start; ; ) {
        char* allele_end = AdvToDelimOrEnd(allele_iter, token_end, ',');
        *allele_end = '\0';
        last_allele_blen = 1 + S_CAST(uintptr_t, allele_end - allele_iter);
        for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
          if (memequal(allele_iter, cur_alleles[allele_idx], last_allele_blen)) {
            if (unlikely(IsSet(alleles_seen, allele_idx))) {
              snprintf(g_logbuf, kLogbufSize, "Error: --alt-allele: Duplicate allele on line %" PRIuPTR ".\n", line_idx);
              goto SetRefalt1FromFile_ret_MALFORMED_INPUT_2;
            }
            SetBit(allele_idx, alleles_seen);
            loaded_ac_order[loaded_allele_ct] = allele_idx;
            ++loaded_allele_ct;
            break;
          }
        }
        if (allele_end == token_end) {
          break;
        }
        allele_iter = &(allele_end[1]);
      }
      if (unlikely(loaded_allele_ct == cur_allele_ct)) {
        snprintf(g_logbuf, kLogbufSize, "Error: --alt-allele: All alleles appear on line %" PRIuPTR "; nothing left over to assign to REF.\n", line_idx);
        goto SetRefalt1FromFile_ret_MALFORMED_INPUT_2;
      }
      // Note that when both --ref-allele and --alt[1]-allele are present in
      // the same run, --alt[1]-allele must deal with the possibility of a
      // pre-altered allele_permute[].
      AlleleCode* cur_allele_permute = &(allele_permute[allele_idx_offset_base]);
      // This is usually 0 for --ref-allele, and 1 for --alt[1]-allele.
      // The exception is when --alt[1]-allele is being resolved after
      // --ref-allele swapped alleles: then it's 0.
      const uint32_t orig_main_allele_idx = cur_allele_permute[is_alt_or_alt1];

      if (loaded_allele_ct == 0) {
        if (cur_allele_ct > 2) {
          // could happen millions of times, so micromanage this instead of
          // using logpreprintfww()
          char* write_iter = strcpya_k(g_logbuf, "Warning: ");
          // strlen("--ref-allele") == strlen("--alt-allele")== 12,
          // strlen("--alt1-allele") == 13
          write_iter = memcpya(write_iter, flagstr, 12 + is_alt1);
          write_iter = strcpya_k(write_iter, " mismatch for multiallelic variant");
          // If we put this all on one line, its length would be
          // 60 + is_alt_or_alt1 + variant_id_slen.  Split into two if this is
          // >79.
          *write_iter++ = (variant_id_slen + is_alt1 < 20)? ' ' : '\n';
          *write_iter++ = '\'';
          write_iter = memcpya(write_iter, variant_ids[variant_uidx], variant_id_slen);
          strcpy_k(write_iter, "'.\n");
          if (allele_mismatch_warning_ct < 3) {
            logerrputsb();
          } else {
            logputs_silent(g_logbuf);
          }
          ++allele_mismatch_warning_ct;
          continue;
        }
        const char** new_allele_ptr = &(cur_alleles[orig_main_allele_idx]);
        const char* main_allele = *new_allele_ptr;
        const uint32_t multi_alt = (last_allele_blen <= allele_col_slen);
        uint32_t is_ref_changing = 1 - orig_main_allele_idx;
        // If main allele is missing, we fill it in (unless there were multiple
        // comma-separated ALT values).
        if ((!strequal_k_unsafe(main_allele, ".")) || multi_alt) {
          new_allele_ptr = &(cur_alleles[1 - orig_main_allele_idx]);
          const char* other_allele = *new_allele_ptr;
          if ((!strequal_k_unsafe(other_allele, ".")) || multi_alt) {
            char* write_iter = strcpya_k(g_logbuf, "Warning: ");
            write_iter = memcpya(write_iter, flagstr, 12 + is_alt1);
            write_iter = strcpya_k(write_iter, " mismatch for biallelic variant '");
            write_iter = strcpya(write_iter, variant_ids[variant_uidx]);
            strcpy_k(write_iter, "'.\n");
            if (allele_mismatch_warning_ct < 3) {
              logerrputsb();
            } else {
              logputs_silent(g_logbuf);
            }
            ++allele_mismatch_warning_ct;
            continue;
          }
          // If we get here during --ref-allele processing, that means REF is
          // nonmissing and unequal, and ALT is missing.  Fill ALT and then
          // swap REF/ALT.
          // If we get here during --alt[1]-allele processing:
          // * If --ref-allele previously swapped the alleles, that would mean
          //   REF is nonmissing and unequal, and ALT is simultaneously missing
          //   and equal to a nonmissing --ref-allele setting... impossible.
          // * So that didn't happen; instead, ALT is nonmissing and unequal,
          //   and REF is missing.  Fill REF and then swap REF/ALT.  Also
          //   (bugfix, 3 Jul 2024) set nonref_flags bit unless --ref-allele
          //   explicitly agreed with the ALT allele.
          cur_allele_permute[0] = 1;
          cur_allele_permute[1] = 0;
          ++permuted_variant_ct;
          is_ref_changing = 1;
        }
        if (is_ref_changing) {
          if (!is_alt_or_alt1) {
            if (!IsSet(nonref_flags, variant_uidx)) {
              if (unlikely(!force)) {
                goto SetRefalt1FromFile_ret_NOFORCE;
              }
            } else {
              ClearBit(variant_uidx, nonref_flags);
            }
          } else {
            // If --ref-allele had something to say about this variant,
            // --alt[1]-allele 'force' should not apply, and we should not set
            // a nonref_flags bit.
            // Otherwise, both are in play.
            if ((!previously_seen) || (!IsSet(previously_seen, variant_uidx))) {
              if (!IsSet(nonref_flags, variant_uidx)) {
                if (unlikely(!force)) {
                  goto SetRefalt1FromFile_ret_NOFORCE;
                }
                SetBit(variant_uidx, nonref_flags);
              }
            }
          }
        }
        if (last_allele_blen == 2) {
          *new_allele_ptr = &(g_one_char_strs[2 * ctou32(allele_start[0])]);
        } else {
          // No in-place-overwrite case here since
          // 1. *new_allele_ptr was always '.'
          // 2. More importantly, when --loop-cats is used,
          //    allele_storage_backup[n] must point to the unaltered original
          //    string.
          if (unlikely(S_CAST(uintptr_t, bigstack_end - main_bigstack_base) < last_allele_blen)) {
            goto SetRefalt1FromFile_ret_NOMEM;
          }
          if (last_allele_blen > max_allele_blen) {
            max_allele_blen = last_allele_blen;
          }
          bigstack_end -= last_allele_blen;
          memcpy(bigstack_end, allele_start, last_allele_blen);
          *new_allele_ptr = R_CAST(const char*, bigstack_end);
        }
        ++fillin_variant_ct;
        continue;
      }
      const uint32_t first_allele_idx = loaded_ac_order[0];
      if (first_allele_idx == orig_main_allele_idx) {
        uint32_t mismatch_found = 0;
        for (uint32_t alt_idx = 1; alt_idx != loaded_allele_ct; ++alt_idx) {
          if (loaded_ac_order[alt_idx] != cur_allele_permute[alt_idx + 1]) {
            mismatch_found = 1;
          }
        }
        if (!mismatch_found) {
          continue;
        }
      }
      // We're permuting some alleles.
      if (is_alt_or_alt1 && previously_seen && IsSet(previously_seen, variant_uidx)) {
        if (unlikely(IsSet(alleles_seen, cur_allele_permute[0]))) {
          // both --ref-allele and --alt[1]-allele in current run, and they
          // contradict each other.  error out instead of producing a
          // order-of-operations dependent result.
          snprintf(g_logbuf, kLogbufSize, "Error: --ref-allele and --alt[1]-allele assignments conflict for variant '%s'.\n", variant_ids[variant_uidx]);
          goto SetRefalt1FromFile_ret_INCONSISTENT_INPUT_WW;
        }
        // minor bugfix (3 Jul 2024): do not perform further nonref_flags
        // check; only --ref-allele 'force' setting should apply.
      } else {
        // * If !is_alt_or_alt1, we are directly changing the REF allele.
        // * If IsSet(alleles_seen, 0), we're processing --alt[1]-allele
        //   without a --ref-allele observation at this variant, and the REF
        //   allele also has to change.  We set it to the first free ALT
        //   allele, and mark it as provisional-reference.
        if ((!is_alt_or_alt1) || IsSet(alleles_seen, 0)) {
          if (!IsSet(nonref_flags, variant_uidx)) {
            if (unlikely(!force)) {
              goto SetRefalt1FromFile_ret_NOFORCE;
            }
          } else if (!is_alt_or_alt1) {
            ClearBit(variant_uidx, nonref_flags);
          }
          if (is_alt_or_alt1) {
            const uint32_t first_free_ac = AdvTo0Bit(alleles_seen, 0);
            cur_allele_permute[0] = first_free_ac;
            SetBit(variant_uidx, nonref_flags);
          }
        }
      }
      if (!is_alt_or_alt1) {
        cur_allele_permute[0] = first_allele_idx;
        for (uint32_t write_ac = 0; write_ac != first_allele_idx; ++write_ac) {
          cur_allele_permute[1 + write_ac] = write_ac;
        }
      } else {
        memcpy(&(cur_allele_permute[1]), loaded_ac_order, loaded_allele_ct * sizeof(AlleleCode));
        if (loaded_allele_ct + 1 != cur_allele_ct) {
          SetBit(cur_allele_permute[0], alleles_seen);

          // No need to update final unobserved alleles, if any.
          const uint32_t write_aidx_stop = FindLast1BitBefore(alleles_seen, cur_allele_ct) + 1;

          uintptr_t allele_idx_base = 0;
          uintptr_t cur_inv_bits = ~alleles_seen[0];
          for (uint32_t write_aidx = loaded_allele_ct + 1; write_aidx != write_aidx_stop; ++write_aidx) {
            const uintptr_t ac = BitIter0(alleles_seen, &allele_idx_base, &cur_inv_bits);
            cur_allele_permute[write_aidx] = ac;
          }
        }
      }
      ++permuted_variant_ct;
    }
    if (allele_mismatch_warning_ct > 3) {
      fprintf(stderr, "%u more allele-mismatch warning%s: see log file.\n", allele_mismatch_warning_ct - 3, (allele_mismatch_warning_ct == 4)? "" : "s");
    }
    if (fillin_variant_ct) {
      if (permuted_variant_ct) {
        logprintfww("%s: %u set%s of allele codes permuted, %u allele code%s filled in.\n", flagstr, permuted_variant_ct, (permuted_variant_ct == 1)? "" : "s", fillin_variant_ct, (fillin_variant_ct == 1)? "" : "s");
      } else {
        logprintf("%s: %u allele code%s filled in.\n", flagstr, fillin_variant_ct, (fillin_variant_ct == 1)? "" : "s");
      }
    } else if (permuted_variant_ct) {
      logprintf("%s: %u set%s of allele codes permuted.\n", flagstr, permuted_variant_ct, (permuted_variant_ct == 1)? "" : "s");
    } else {
      logprintf("%s: No variants changed.\n", flagstr);
    }
    if (skipped_variant_ct) {
      logerrprintfww("Warning: %" PRIuPTR " variant ID%s in %s file missing from main dataset.\n", skipped_variant_ct, (skipped_variant_ct == 1)? "" : "s", flagstr);
    }
    if (missing_allele_ct) {
      logerrprintfww("Warning: %" PRIuPTR " allele code%s in %s file were missing (%s skipped).\n", missing_allele_ct, (missing_allele_ct == 1)? "" : "s", flagstr, (missing_allele_ct == 1)? "this entry was" : "these entries were");
    }
    if (previously_seen && (!is_alt_or_alt1)) {
      memcpy(previously_seen, already_seen, raw_variant_ctl * sizeof(intptr_t));
    }
    // bugfix (19 Jun 2018): forgot to update max_allele_slen.
    *max_allele_slen_ptr = max_allele_blen - 1;
  }
  while (0) {
  SetRefalt1FromFile_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SetRefalt1FromFile_ret_TSTREAM_FAIL:
    {
      char file_descrip_buf[32];
      snprintf(file_descrip_buf, 32, "%s file", flagstr);
      TextStreamErrPrint(file_descrip_buf, &txs);
    }
    break;
  SetRefalt1FromFile_ret_MISSING_TOKENS:
    logerrprintfww("Error: Line %" PRIuPTR " of %s file has fewer tokens than expected.\n", line_idx, flagstr);
    reterr = kPglRetMalformedInput;
    break;
  SetRefalt1FromFile_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
  SetRefalt1FromFile_ret_MALFORMED_INPUT_2:
    logerrputsb();
    reterr = kPglRetMalformedInput;
    break;
  SetRefalt1FromFile_ret_NOFORCE:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s file contradicts 'known' reference allele. Add the 'force' modifier to %s to force an allele swap anyway.\n", line_idx, flagstr, flagstr);
  SetRefalt1FromFile_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 SetRefalt1FromFile_ret_1:
  if (CleanupTextStream(&txs, &reterr)) {
    logerrprintfww("Error: %s file read failure: %s.\n", flagstr, rstrerror(errno));
  }
  BigstackReset(bigstack_mark);
  BigstackEndSet(bigstack_end);
  return reterr;
}

PglErr MakeFounders(const uintptr_t* sample_include, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t require_two, PedigreeIdInfo* piip, uintptr_t* founder_info) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    char* idbuf;
    uintptr_t* nf_bitvec;
    if (bigstack_alloc_c(max_sample_id_blen, &idbuf) ||
        bigstack_alloc_w(raw_sample_ctl, &nf_bitvec)) {
      goto MakeFounders_ret_NOMEM;
    }
    BitvecInvmaskCopy(sample_include, founder_info, raw_sample_ctl, nf_bitvec);
    uint32_t sample_uidx = AdvBoundedTo1Bit(nf_bitvec, 0, raw_sample_ct);
    if (sample_uidx == raw_sample_ct) {
      logputs("Note: Skipping --make-founders since there are no nonfounders.\n");
      goto MakeFounders_ret_1;
    }
    const char* sample_ids = piip->sii.sample_ids;
    uint32_t* id_htable;
    uint32_t id_htable_size;
    if (unlikely(HtableGoodSizeAlloc(sample_ct, bigstack_left(), &id_htable, &id_htable_size))) {
      goto MakeFounders_ret_NOMEM;
    }
    PopulateStrboxSubsetHtable(sample_ids, sample_include, sample_ct, max_sample_id_blen, 1, id_htable_size, id_htable);
    const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
    char* paternal_ids = piip->parental_id_info.paternal_ids;
    char* maternal_ids = piip->parental_id_info.maternal_ids;
    uint32_t new_founder_ct = 0;
    do {
      const char* fid_start = &(sample_ids[sample_uidx * max_sample_id_blen]);
      const char* fid_postend = AdvPastDelim(fid_start, '\t');
      const uint32_t fid_blen = fid_postend - fid_start;
      char* iid_write_start = memcpya(idbuf, fid_start, fid_blen);

      uint32_t missing_parent_ct = 0;
      char* pat_ptr = &(paternal_ids[sample_uidx * max_paternal_id_blen]);
      const uint32_t pat_slen = strlen(pat_ptr);
      if (fid_blen + pat_slen >= max_sample_id_blen) {
        ++missing_parent_ct;
      } else {
        memcpy(iid_write_start, pat_ptr, pat_slen + 1);
        const uint32_t uii = StrboxHtableFind(idbuf, sample_ids, id_htable, max_sample_id_blen, fid_blen + pat_slen, id_htable_size);
        if (uii == UINT32_MAX) {
          ++missing_parent_ct;
        }
      }
      char* mat_ptr = &(maternal_ids[sample_uidx * max_maternal_id_blen]);
      const uint32_t mat_slen = strlen(mat_ptr);
      if (fid_blen + mat_slen >= max_sample_id_blen) {
        ++missing_parent_ct;
      } else {
        memcpy(iid_write_start, mat_ptr, mat_slen + 1);
        const uint32_t uii = StrboxHtableFind(idbuf, sample_ids, id_htable, max_sample_id_blen, fid_blen + mat_slen, id_htable_size);
        if (uii == UINT32_MAX) {
          ++missing_parent_ct;
        }
      }
      if (missing_parent_ct > require_two) {
        SetBit(sample_uidx, founder_info);
        strcpy_k(pat_ptr, "0");
        strcpy_k(mat_ptr, "0");
        ++new_founder_ct;
      }
      sample_uidx = AdvBoundedTo1Bit(nf_bitvec, sample_uidx + 1, raw_sample_ct);
    } while (sample_uidx != raw_sample_ct);
    logprintf("--make-founders: %u sample%s affected.\n", new_founder_ct, (new_founder_ct == 1)? "" : "s");
  }
  while (0) {
  MakeFounders_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 MakeFounders_ret_1:
  BigstackReset(bigstack_mark);
  return reterr;
}

static_assert(sizeof(Dosage) + sizeof(AlleleCode) <= 4, "MajRef() needs to be updated.");
PglErr MajRef(const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const uint64_t* main_allele_ddosages, uint32_t variant_ct, uint32_t max_allele_ct, uint32_t skip_real_ref, AlleleCode* allele_permute, uintptr_t* nonref_flags) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  {
    uint64_t* sortbuf = nullptr;
    if (allele_idx_offsets) {
      if (unlikely(bigstack_alloc_u64(max_allele_ct, &sortbuf))) {
        goto MajRef_ret_NOMEM;
      }
    }
    // todo: multithread this
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t allele_ct = 2;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      // bugfix (25 Jun 2024): nonref_flags check was backwards
      if (skip_real_ref && (!IsSet(nonref_flags, variant_uidx))) {
        continue;
      }
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const uint64_t* cur_allele_ddosages = &(main_allele_ddosages[allele_idx_offset_base]);
      if (allele_ct == 2) {
        if (cur_allele_ddosages[1] > cur_allele_ddosages[0]) {
          allele_permute[allele_idx_offset_base] = 1;
          allele_permute[allele_idx_offset_base + 1] = 0;
          if (nonref_flags) {
            SetBit(variant_uidx, nonref_flags);
          }
        }
        continue;
      }
      // Sort by allele frequency, breaking ties in favor of lower aidx.
      // sortbuf entries have a (shared) multiple of the allele frequency in
      // high bits, and (kPglMaxAlleleCt - aidx) in low bits.
      for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
        sortbuf[aidx] = cur_allele_ddosages[aidx] * (kPglMaxAlleleCt + 1) + kPglMaxAlleleCt - aidx;
      }
      STD_SORT(allele_ct, u64cmp, sortbuf);
      AlleleCode* cur_allele_permute = &(allele_permute[allele_idx_offset_base]);
      const uint32_t allele_ct_m1 = allele_ct - 1;
      for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
        cur_allele_permute[aidx] = kPglMaxAlleleCt - S_CAST(AlleleCode, sortbuf[allele_ct_m1 - aidx]);
      }
      if (nonref_flags && cur_allele_permute[0]) {
        SetBit(variant_uidx, nonref_flags);
      }
    }
  }
  while (0) {
  MajRef_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
