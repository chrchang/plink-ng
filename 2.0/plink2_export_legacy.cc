// This file is part of PLINK 2.00, copyright (C) 2005-2022 Shaun Purcell,
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


#include "plink2_common.h"
#include "plink2_compress_stream.h"

// This covers formats that are fully supported by PLINK 1.x (no multiallelic
// variants, dosages, or phase information).

#ifdef __cplusplus
namespace plink2 {
#endif

typedef struct TransposeToSmajReadCtxStruct {
  const uintptr_t* variant_include;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select);
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  uint32_t sample_ct;

  PgenReader** pgr_ptrs;

  uint32_t* variant_uidx_starts;
  uint32_t cur_block_write_ct;

  uintptr_t* vmaj_readbuf;

  uint64_t err_info;
} TransposeToSmajReadCtx;

THREAD_FUNC_DECL TransposeToSmajReadThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  TransposeToSmajReadCtx* ctx = S_CAST(TransposeToSmajReadCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select) = ctx->refalt1_select;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uintptr_t* sample_include = ctx->sample_include;
  const uint32_t read_sample_ct = ctx->sample_ct;
  const uintptr_t read_sample_ctaw2 = NypCtToAlignedWordCt(read_sample_ct);
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  uintptr_t prev_copy_ct = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_copy_ct = ctx->cur_block_write_ct;
    const uint32_t cur_idx_end = ((tidx + 1) * cur_block_copy_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, ctx->variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
    const uint32_t cur_idx_start = (tidx * cur_block_copy_ct) / calc_thread_ct;
    uintptr_t* vmaj_readbuf_iter = &(ctx->vmaj_readbuf[(prev_copy_ct + cur_idx_start) * read_sample_ctaw2]);
    for (uint32_t cur_idx = cur_idx_start; cur_idx != cur_idx_end; ++cur_idx) {
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      // todo: multiallelic case
      const PglErr reterr = PgrGet(sample_include, pssi, read_sample_ct, variant_uidx, pgrp, vmaj_readbuf_iter);
      if (unlikely(reterr)) {
        new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
        goto TransposeToSmajReadThread_err;
      }
      if (refalt1_select && (refalt1_select[variant_uidx][0] == 1)) {
        GenovecInvertUnsafe(read_sample_ct, vmaj_readbuf_iter);
        // don't need ZeroTrailingNyps()
      }
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[read_sample_ctaw2]);
    }
    prev_copy_ct += cur_block_copy_ct;
    while (0) {
    TransposeToSmajReadThread_err:
      UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

typedef struct TransposeToPlink1SmajWriteCtxStruct {
  const uintptr_t* variant_include;
  uint32_t variant_ct;
  uint32_t sample_ct;

  const uintptr_t* vmaj_readbuf;

  VecW** thread_vecaligned_bufs;

  uint32_t sample_batch_size;

  uintptr_t* smaj_writebufs[2];
} TransposeToPlink1SmajWriteCtx;

THREAD_FUNC_DECL TransposeToPlink1SmajWriteThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  TransposeToPlink1SmajWriteCtx* ctx = S_CAST(TransposeToPlink1SmajWriteCtx*, arg->sharedp->context);

  const uint32_t variant_ct = ctx->variant_ct;
  const uintptr_t variant_batch_ct = DivUp(variant_ct, kPglNypTransposeBatch);
  const uintptr_t variant_batch_word_ct = variant_batch_ct * kPglNypTransposeWords;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uintptr_t variant_batch_idx_start = (S_CAST(uint64_t, tidx) * variant_batch_ct) / calc_thread_ct;
  VecW* vecaligned_buf = ctx->thread_vecaligned_bufs[tidx];
  uintptr_t variant_batch_idx_full_end = ((S_CAST(uint64_t, tidx) + 1) * variant_batch_ct) / calc_thread_ct;
  uint32_t variant_idx_end;
  if (tidx + 1 < calc_thread_ct) {
    variant_idx_end = variant_batch_idx_full_end * kPglNypTransposeBatch;
  } else {
    variant_idx_end = variant_ct;
    if (variant_ct % kPglNypTransposeBatch) {
      --variant_batch_idx_full_end;
    }
  }
  const uint32_t thread_variant_ct = variant_idx_end - variant_batch_idx_start * kPglNypTransposeBatch;
  const uint32_t read_sample_ct = ctx->sample_ct;
  const uintptr_t read_sample_ctaw2 = NypCtToAlignedWordCt(read_sample_ct);
  const uintptr_t* vmaj_readbuf = ctx->vmaj_readbuf;
  uint32_t sample_widx = 0;
  uint32_t parity = 0;
  do {
    const uint32_t sample_batch_size = ctx->sample_batch_size;
    const uintptr_t* vmaj_readbuf_iter = &(vmaj_readbuf[variant_batch_idx_start * kPglNypTransposeBatch * read_sample_ctaw2 + sample_widx]);
    uintptr_t* smaj_writebuf_start = &(ctx->smaj_writebufs[parity][variant_batch_idx_start * kPglNypTransposeWords]);
    uintptr_t* smaj_writebuf_iter = smaj_writebuf_start;
    uint32_t variant_batch_size = kPglNypTransposeBatch;
    for (uintptr_t variant_batch_idx = variant_batch_idx_start; ; ++variant_batch_idx) {
      if (variant_batch_idx >= variant_batch_idx_full_end) {
        if (variant_batch_idx * kPglNypTransposeBatch >= variant_idx_end) {
          break;
        }
        variant_batch_size = variant_idx_end - variant_batch_idx * kPglNypTransposeBatch;
      }
      TransposeNypblock(vmaj_readbuf_iter, read_sample_ctaw2, variant_batch_word_ct, variant_batch_size, sample_batch_size, smaj_writebuf_iter, vecaligned_buf);
      smaj_writebuf_iter = &(smaj_writebuf_iter[kPglNypTransposeWords]);
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[variant_batch_size * read_sample_ctaw2]);
    }
    smaj_writebuf_iter = smaj_writebuf_start;
    for (uint32_t sample_idx = 0; sample_idx != sample_batch_size; ++sample_idx) {
      // could fold this into TransposeNypblock(), but I won't bother,
      // we're already saturating at ~3 threads
      PgrPlink2ToPlink1InplaceUnsafe(thread_variant_ct, smaj_writebuf_iter);
      ZeroTrailingNyps(thread_variant_ct, smaj_writebuf_iter);
      smaj_writebuf_iter = &(smaj_writebuf_iter[variant_batch_word_ct]);
    }
    parity = 1 - parity;
    sample_widx += sample_batch_size / kBitsPerWordD2;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

PglErr ExportIndMajorBed(const uintptr_t* orig_sample_include, const uintptr_t* variant_include, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup read_tg;
  ThreadGroup write_tg;
  PreinitThreads(&read_tg);
  PreinitThreads(&write_tg);
  TransposeToSmajReadCtx read_ctx;
  TransposeToPlink1SmajWriteCtx write_ctx;
  {
    // Possible special case: if the input file is a variant-major .bed, we do
    // not have enough memory to just load the whole file at once, and there
    // are more than ~20k samples, there can be a performance advantage to not
    // loading an entire variant at a time; we can use smaller fread calls and
    // reduce the number of (typically 4096 byte) disk blocks which need to be
    // read on each pass.  But let's get .pgen -> sample-major humming first.
    snprintf(outname_end, kMaxOutfnameExtBlen, ".bed");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportIndMajorBed_ret_OPEN_FAIL;
    }
    if (unlikely(!fwrite_unlocked("l\x1b", 3, 1, outfile))) {
      goto ExportIndMajorBed_ret_WRITE_FAIL;
    }
    if (variant_ct && sample_ct) {
      const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
      uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      // todo: if only 1 pass is needed, and no subsetting is happening, this
      // saturates at ~4 threads?
      STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
      uint32_t read_block_size;
      // note that this is restricted to half of available workspace
      if (unlikely(PgenMtLoadInit(variant_include, sample_ct, variant_ct, bigstack_left() / 2, pgr_alloc_cacheline_ct, 0, 0, 0, pgfip, &calc_thread_ct, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &read_ctx.pgr_ptrs, &read_ctx.variant_uidx_starts))) {
        goto ExportIndMajorBed_ret_NOMEM;
      }
      if (unlikely(SetThreadCt(calc_thread_ct, &read_tg))) {
        goto ExportIndMajorBed_ret_NOMEM;
      }
      read_ctx.variant_include = variant_include;
      read_ctx.refalt1_select = refalt1_select;
      read_ctx.err_info = (~0LLU) << 32;
      SetThreadFuncAndData(TransposeToSmajReadThread, &read_ctx, &read_tg);

      const uintptr_t variant_cacheline_ct = NypCtToCachelineCt(variant_ct);
      uint32_t output_calc_thread_ct = MINV(calc_thread_ct, variant_cacheline_ct);
      // 4 still seems to be best in AVX2 case
      if (output_calc_thread_ct > 4) {
        output_calc_thread_ct = 4;
      }
      uintptr_t* sample_include;
      uint32_t* sample_include_cumulative_popcounts;
      if (unlikely(SetThreadCt(output_calc_thread_ct, &write_tg) ||
                   bigstack_alloc_w(raw_sample_ctl, &sample_include) ||
                   bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                   bigstack_alloc_vp(output_calc_thread_ct, &write_ctx.thread_vecaligned_bufs))) {
        goto ExportIndMajorBed_ret_NOMEM;
      }
      for (uint32_t tidx = 0; tidx != output_calc_thread_ct; ++tidx) {
        write_ctx.thread_vecaligned_bufs[tidx] = S_CAST(VecW*, bigstack_alloc_raw(kPglNypTransposeBufbytes));
      }
      // each of the two write buffers should use <= 1/8 of the remaining
      // workspace
      const uintptr_t writebuf_cachelines_avail = bigstack_left() / (kCacheline * 8);
      uint32_t sample_batch_size = kPglNypTransposeBatch;
      if (variant_cacheline_ct * kPglNypTransposeBatch > writebuf_cachelines_avail) {
        sample_batch_size = RoundDownPow2(writebuf_cachelines_avail / variant_cacheline_ct, kBitsPerWordD2);
        if (unlikely(!sample_batch_size)) {
          goto ExportIndMajorBed_ret_NOMEM;
        }
      }
      write_ctx.smaj_writebufs[0] = S_CAST(uintptr_t*, bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size));
      write_ctx.smaj_writebufs[1] = S_CAST(uintptr_t*, bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size));
      const uintptr_t readbuf_vecs_avail = (bigstack_left() / kCacheline) * kVecsPerCacheline;
      if (unlikely(readbuf_vecs_avail < variant_ct)) {
        goto ExportIndMajorBed_ret_NOMEM;
      }
      uintptr_t read_sample_ctv2 = readbuf_vecs_avail / variant_ct;
      uint32_t read_sample_ct;
      if (read_sample_ctv2 >= NypCtToVecCt(sample_ct)) {
        read_sample_ct = sample_ct;
      } else {
        read_sample_ct = read_sample_ctv2 * kNypsPerVec;
      }
      uintptr_t read_sample_ctaw2 = NypCtToAlignedWordCt(read_sample_ct);
      uintptr_t* vmaj_readbuf = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(variant_ct * read_sample_ctaw2 * kBytesPerWord));
      read_ctx.vmaj_readbuf = vmaj_readbuf;
      write_ctx.variant_include = variant_include;
      write_ctx.variant_ct = variant_ct;
      write_ctx.vmaj_readbuf = vmaj_readbuf;
      SetThreadFuncAndData(TransposeToPlink1SmajWriteThread, &write_ctx, &write_tg);
      uint32_t sample_uidx_start = AdvTo1Bit(orig_sample_include, 0);
      const uintptr_t variant_ct4 = NypCtToByteCt(variant_ct);
      const uintptr_t variant_ctaclw2 = variant_cacheline_ct * kWordsPerCacheline;
      const uint32_t pass_ct = 1 + (sample_ct - 1) / read_sample_ct;
      for (uint32_t pass_idx = 0; pass_idx != pass_ct; ++pass_idx) {
        memcpy(sample_include, orig_sample_include, raw_sample_ctl * sizeof(intptr_t));
        if (sample_uidx_start) {
          ClearBitsNz(0, sample_uidx_start, sample_include);
        }
        uint32_t sample_uidx_end;
        if (pass_idx + 1 == pass_ct) {
          read_sample_ct = sample_ct - pass_idx * read_sample_ct;
          read_sample_ctaw2 = NypCtToAlignedWordCt(read_sample_ct);
          sample_uidx_end = raw_sample_ct;
        } else {
          sample_uidx_end = FindNth1BitFrom(orig_sample_include, sample_uidx_start + 1, read_sample_ct);
          ClearBitsNz(sample_uidx_end, raw_sample_ct, sample_include);
        }
        FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
        read_ctx.sample_include = sample_include;
        read_ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
        read_ctx.sample_ct = read_sample_ct;
        write_ctx.sample_ct = read_sample_ct;
        if (pass_idx) {
          pgfip->block_base = main_loadbufs[0];
          // er, don't need SetBaseAndOffset0?
          PgrSetBaseAndOffset0(main_loadbufs[0], calc_thread_ct, read_ctx.pgr_ptrs);
        }
        uint32_t parity = 0;
        uint32_t read_block_idx = 0;
        ReinitThreads(&read_tg);
        uint32_t pct = 0;
        uint32_t next_print_idx = variant_ct / 100;
        putc_unlocked('\r', stdout);
        printf("--export ind-major-bed pass %u/%u: loading... 0%%", pass_idx + 1, pass_ct);
        fflush(stdout);
        for (uint32_t variant_idx = 0; ; ) {
          const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &read_tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
          if (unlikely(reterr)) {
            goto ExportIndMajorBed_ret_PGR_FAIL;
          }
          if (variant_idx) {
            JoinThreads(&read_tg);
            reterr = S_CAST(PglErr, read_ctx.err_info);
            if (unlikely(reterr)) {
              PgenErrPrintNV(reterr, read_ctx.err_info >> 32);
              goto ExportIndMajorBed_ret_1;
            }
          }
          if (!IsLastBlock(&read_tg)) {
            read_ctx.cur_block_write_ct = cur_block_write_ct;
            ComputeUidxStartPartition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, read_ctx.variant_uidx_starts);
            PgrCopyBaseAndOffset(pgfip, calc_thread_ct, read_ctx.pgr_ptrs);
            if (variant_idx + cur_block_write_ct == variant_ct) {
              DeclareLastThreadBlock(&read_tg);
            }
            if (unlikely(SpawnThreads(&read_tg))) {
              goto ExportIndMajorBed_ret_THREAD_CREATE_FAIL;
            }
          }
          parity = 1 - parity;
          if (variant_idx == variant_ct) {
            break;
          }
          if (variant_idx >= next_print_idx) {
            if (pct > 10) {
              putc_unlocked('\b', stdout);
            }
            pct = (variant_idx * 100LLU) / variant_ct;
            printf("\b\b%u%%", pct++);
            fflush(stdout);
            next_print_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
          }

          ++read_block_idx;
          variant_idx += cur_block_write_ct;
          pgfip->block_base = main_loadbufs[parity];
        }
        // 2. Transpose and write.  (Could parallelize some of the transposing
        //    with the read loop, but since we can't write a single row until
        //    the read loop is done, and both write speed and write buffer
        //    space are bottlenecks, that can't be expected to help much.)
        ReinitThreads(&write_tg);
        write_ctx.sample_batch_size = sample_batch_size;
        parity = 0;
        if (pct > 10) {
          fputs("\b \b", stdout);
        }
        fputs("\b\b\b\b\b\b\b\b\b\b\b\b\bwriting... 0%", stdout);
        fflush(stdout);
        pct = 0;
        uint32_t flush_sample_idx = 0;
        next_print_idx = read_sample_ct / 100;
        for (uint32_t flush_sample_idx_end = 0; ; ) {
          if (!IsLastBlock(&write_tg)) {
            if (flush_sample_idx_end + sample_batch_size >= read_sample_ct) {
              DeclareLastThreadBlock(&write_tg);
              write_ctx.sample_batch_size = read_sample_ct - flush_sample_idx_end;
            }
            if (unlikely(SpawnThreads(&write_tg))) {
              goto ExportIndMajorBed_ret_THREAD_CREATE_FAIL;
            }
          }
          if (flush_sample_idx_end) {
            uintptr_t* smaj_writebuf_iter = write_ctx.smaj_writebufs[1 - parity];
            for (; flush_sample_idx != flush_sample_idx_end; ++flush_sample_idx) {
              fwrite_unlocked(smaj_writebuf_iter, variant_ct4, 1, outfile);
              smaj_writebuf_iter = &(smaj_writebuf_iter[variant_ctaclw2]);
            }
            if (unlikely(ferror_unlocked(outfile))) {
              goto ExportIndMajorBed_ret_WRITE_FAIL;
            }
            if (flush_sample_idx_end == read_sample_ct) {
              break;
            }
            if (flush_sample_idx_end >= next_print_idx) {
              if (pct > 10) {
                putc_unlocked('\b', stdout);
              }
              pct = (flush_sample_idx_end * 100LLU) / read_sample_ct;
              printf("\b\b%u%%", pct++);
              fflush(stdout);
              next_print_idx = (pct * S_CAST(uint64_t, read_sample_ct)) / 100;
            }
          }
          JoinThreads(&write_tg);
          parity = 1 - parity;
          flush_sample_idx_end += sample_batch_size;
          if (flush_sample_idx_end > read_sample_ct) {
            flush_sample_idx_end = read_sample_ct;
          }
        }
        if (pct > 10) {
          fputs("\b \b", stdout);
        }
        sample_uidx_start = sample_uidx_end;
      }
      fputs("\b\bdone.\n", stdout);
    }
    if (unlikely(fclose_null(&outfile))) {
      goto ExportIndMajorBed_ret_WRITE_FAIL;
    }
    logprintfww("--export ind-major-bed: %s written.\n", outname);
  }
  while (0) {
  ExportIndMajorBed_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportIndMajorBed_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportIndMajorBed_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ExportIndMajorBed_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportIndMajorBed_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 ExportIndMajorBed_ret_1:
  CleanupThreads(&write_tg);
  CleanupThreads(&read_tg);
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr ExportTped(const char* outname, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const double* variant_cms, uint32_t sample_ct, uint32_t variant_ct, uint32_t max_allele_slen, char exportf_delim, char lomg_char, PgenReader* simple_pgrp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportTped_ret_OPEN_FAIL;
    }
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const uint32_t max_chr_blen = 1 + GetMaxChrSlen(cip);
    const uintptr_t writebuf_size = RoundUpPow2(kMaxMediumLine + MAXV(2 * kMaxIdSlen + 32, MAXV(2 * max_allele_slen + 2, 4 * sample_ct)), kCacheline);
    char* chr_buf;
    char* writebuf;
    uintptr_t* genovec;
    char* genotext_buf;
    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf) ||
                 bigstack_alloc_c(writebuf_size, &writebuf) ||
                 bigstack_alloc_w(sample_ctl2, &genovec) ||
                 bigstack_alloc_c((6 * k1LU) * max_allele_slen + 10, &genotext_buf))) {
      goto ExportTped_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    char* write_iter = writebuf;

    logprintfww5("--export tped to %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t geno_to_text32[4];
    // delim, missing code, delim, missing code
    geno_to_text32[3] = (ctou32(exportf_delim) + 0x100 * ctou32(lomg_char)) * 0x10001;

    char* geno_to_str[4];
    uint32_t geno_slen[4];
    geno_to_str[3] = &(genotext_buf[(6 * k1LU) * max_allele_slen + 6]);
    geno_to_str[3][0] = exportf_delim;
    geno_to_str[3][1] = lomg_char;
    geno_to_str[3][2] = exportf_delim;
    geno_to_str[3][3] = lomg_char;
    geno_slen[3] = 4;

    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_blen = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uintptr_t variant_uidx_base = 0;
    uintptr_t variant_include_bits = variant_include[0];
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
        *chr_name_end++ = exportf_delim;
        chr_blen = chr_name_end - chr_buf;
      }
      char* line_start = write_iter;
      write_iter = memcpya(write_iter, chr_buf, chr_blen);
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], exportf_delim);
      if (variant_cms) {
        write_iter = dtoa_g_p8(variant_cms[variant_uidx], write_iter);
      } else {
        *write_iter++ = '0';
      }
      *write_iter++ = exportf_delim;
      write_iter = u32toa(variant_bps[variant_uidx], write_iter);
      uint64_t line_blen = write_iter - line_start;
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto ExportTped_ret_WRITE_FAIL;
      }

      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        if (unlikely(allele_idx_offsets[variant_uidx + 1] != allele_idx_offset_base + 2)) {
          logputs("\n");
          logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
          goto ExportTped_ret_INCONSISTENT_INPUT;
        }
      }
      reterr = PgrGet(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, genovec);
      if (unlikely(reterr)) {
        PgenErrPrintNV(reterr, variant_uidx);
        goto ExportTped_ret_1;
      }
      const char* ref_allele = allele_storage[allele_idx_offset_base];
      const char* alt_allele = allele_storage[allele_idx_offset_base + 1];
      const uint32_t ref_slen = strlen(ref_allele);
      const uint32_t alt_slen = strlen(alt_allele);
      uint32_t loop_len = kBitsPerWordD2;
      if ((ref_slen == 1) && (alt_slen == 1)) {
        // easy case: each genotype corresponds to 4 bytes
        // could set up a geno_pair_to_text64 array, but I think we're already
        // close enough to I/O-bound as-is
        geno_to_text32[0] = (ctou32(exportf_delim) + 0x100 * ctou32(ref_allele[0])) * 0x10001;
        // plink 1.x put the A1 allele first in the heterozygous case, so put
        // ALT first to approximate that
        geno_to_text32[1] = ctou32(exportf_delim) * 0x10001 + ctou32(alt_allele[0]) * 0x100 + ctou32(ref_allele[0]) * 0x1000000;

        geno_to_text32[2] = (ctou32(exportf_delim) + 0x100 * ctou32(alt_allele[0])) * 0x10001;
#ifdef NO_UNALIGNED
#  error "Unaligned accesses in ExportTped()."
#endif
        uint32_t* write_iter_alias = R_CAST(uint32_t*, write_iter);
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            loop_len = ModNz(sample_ct, kBitsPerWordD2);
          }
          uintptr_t geno_word = genovec[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
            const uintptr_t cur_geno = geno_word & 3;
            *write_iter_alias++ = geno_to_text32[cur_geno];
            geno_word >>= 2;
          }
        }
        write_iter = R_CAST(char*, write_iter_alias);
        line_blen += sample_ct * (4 * k1LU);
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto ExportTped_ret_WRITE_FAIL;
        }
      } else {
        char* genotext_write_iter = genotext_buf;
        geno_to_str[0] = genotext_write_iter;
        *genotext_write_iter++ = exportf_delim;
        genotext_write_iter = memcpyax(genotext_write_iter, ref_allele, ref_slen, exportf_delim);
        genotext_write_iter = memcpya(genotext_write_iter, ref_allele, ref_slen);
        geno_slen[0] = 2 * ref_slen + 2;

        geno_to_str[1] = genotext_write_iter;
        *genotext_write_iter++ = exportf_delim;
        genotext_write_iter = memcpyax(genotext_write_iter, alt_allele, alt_slen, exportf_delim);
        genotext_write_iter = memcpya(genotext_write_iter, ref_allele, ref_slen);
        geno_slen[1] = ref_slen + alt_slen + 2;

        geno_to_str[2] = genotext_write_iter;
        *genotext_write_iter++ = exportf_delim;
        genotext_write_iter = memcpyax(genotext_write_iter, alt_allele, alt_slen, exportf_delim);
        memcpy(genotext_write_iter, alt_allele, alt_slen);
        geno_slen[2] = 2 * alt_slen + 2;

        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            loop_len = ModNz(sample_ct, kBitsPerWordD2);
          }
          uintptr_t geno_word = genovec[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
            const uintptr_t cur_geno = geno_word & 3;
            const uintptr_t genotext_slen = geno_slen[cur_geno];
            write_iter = memcpya(write_iter, geno_to_str[cur_geno], genotext_slen);
            line_blen += genotext_slen;
            if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
              goto ExportTped_ret_WRITE_FAIL;
            }
            geno_word >>= 2;
          }
        }
      }

      AppendBinaryEoln(&write_iter);
      line_blen += strlen(EOLN_STR);
      if (unlikely(line_blen > kMaxLongLine)) {
        logputs("\n");
        logerrputs("Error: --export tped would create an excessively long line.\n");
        goto ExportTped_ret_INCONSISTENT_INPUT;
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
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto ExportTped_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
  }
  while (0) {
  ExportTped_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportTped_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportTped_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportTped_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 ExportTped_ret_1:
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr ExportPed(const char* outname, const uintptr_t* orig_sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const char* legacy_output_missing_pheno, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t compound_genotypes, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, char exportf_delim, char lomg_char, PgenFileInfo* pgfip) {
  // Similar to ExportIndMajorBed() and Export012Smaj().
  //
  // Could be sped up by parallelizing the final rendering step (using
  // __atomic_fetch_add to distribute work in a manner similar to the bgzf
  // compressor), but that isn't realistically worth the effort since nobody
  // should be using this format when speed matters.
  //
  // ...more precisely, the general case isn't worth the effort.  The
  // compound-genotypes and all-SNP ped subcases have better-behaved line
  // lengths, so they're substantially easier to accelerate.
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup read_tg;
  ThreadGroup transpose_tg;
  PreinitThreads(&read_tg);
  PreinitThreads(&transpose_tg);
  TransposeToSmajReadCtx read_ctx;
  TransposeToPlink1SmajWriteCtx transpose_ctx;
  {
    if (unlikely(((4 * k1LU) - compound_genotypes) * variant_ct > kMaxLongLine - 4 * kMaxIdSlen - 32)) {
      logerrprintf("Error: --export %s would create an excessively long line.\n", compound_genotypes? "compound-genotypes" : "ped");
      goto ExportPed_ret_INCONSISTENT_INPUT;
    }
    // Write header line; then fully load-and-transpose the first X samples,
    // flush them, load-and-transpose the next X, etc.
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportPed_ret_OPEN_FAIL;
    }

    uint32_t multichar_allele_present = 0;
    if (max_allele_slen > 1) {
      // Still possible for multichar_allele_present to be false, since
      // max_allele_slen is not updated after variant filtering.
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          // Unnecessary to verify the variant is biallelic, we'll catch
          // violations later.
        }
        const char* ref_allele = allele_storage[allele_idx_offset_base];
        const char* alt_allele = allele_storage[allele_idx_offset_base + 1];
        if (ref_allele[1] || alt_allele[1]) {
          multichar_allele_present = 1;
          break;
        }
      }
    }
    // Preprocessing: precompute all possible genotype strings.  (Yes, this can
    // eat a lot of memory, but ~100% of the time that's actually a problem,
    // the user has no business exporting a .ped anyway.)
    //
    // all-SNPs case: index into this array is of the form (variant_idx * 8) +
    //                cur_geno_pair
    const char* geno_pair_matrix = nullptr;
    // variable-length case: index into this array is of the form
    //                       (variant_idx * 4) + cur_geno.
    const char** genotype_strs = nullptr;

    char* writebuf = g_textbuf;
    if (!multichar_allele_present) {
      // Inner loop writes kBitsPerWordD4 genotypes at a time, and we just move
      // the write pointer backward at the end if the last block is short.
      // So, to avoid possibly copying uninitialized memory, we initialize
      // geno_pair_matrix up to variant_idx_end/2 instead of just
      // (variant_ct+1)/2.
      const uint32_t variant_idx_end = RoundUpPow2(variant_ct, kBitsPerWordD4);
      const uint32_t geno_pair_alloc_ct = variant_idx_end / 2;
      char* geno_pair_matrix_start;
      if (unlikely(bigstack_alloc_c(geno_pair_alloc_ct * ((128 - 32 * compound_genotypes) * k1LU), &geno_pair_matrix_start))) {
        goto ExportPed_ret_NOMEM;
      }
      char geno_strs[2][4][4];
      const uint32_t second_offset = 3 - compound_genotypes;
      for (uint32_t pair_part = 0; pair_part != 2; ++pair_part) {
        for (uint32_t geno_code = 0; geno_code != 4; ++geno_code) {
          geno_strs[pair_part][geno_code][0] = exportf_delim;
          geno_strs[pair_part][geno_code][1] = lomg_char;
          geno_strs[pair_part][geno_code][2] = exportf_delim;
          geno_strs[pair_part][geno_code][second_offset] = lomg_char;
        }
      }

      const uint32_t variant_ct_m1 = variant_ct - 1;
      char* geno_pair_matrix_iter = geno_pair_matrix_start;
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      uint32_t pair_part_ct = 2;
      for (uint32_t variant_idx = 0; ; variant_idx += 2) {
        if (variant_idx >= variant_ct_m1) {
          if (variant_idx > variant_ct_m1) {
            if (variant_idx == variant_idx_end) {
              break;
            }
            pair_part_ct = 0;
          } else {
            pair_part_ct = 1; // last variant
          }
        }
        for (uint32_t pair_part = 0; pair_part != pair_part_ct; ++pair_part) {
          const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
          uintptr_t allele_idx_offset_base = variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[variant_uidx];
            if (unlikely((allele_idx_offsets[variant_uidx + 1] != allele_idx_offset_base + 2))) {
              logputs("\n");
              logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
              goto ExportPed_ret_INCONSISTENT_INPUT;
            }
          }
          const char ref_char = allele_storage[allele_idx_offset_base][0];
          const char alt_char = allele_storage[allele_idx_offset_base + 1][0];
          geno_strs[pair_part][0][1] = alt_char;
          geno_strs[pair_part][0][second_offset] = alt_char;
          geno_strs[pair_part][2][1] = alt_char;
          geno_strs[pair_part][2][second_offset] = ref_char;
          geno_strs[pair_part][3][1] = ref_char;
          geno_strs[pair_part][3][second_offset] = ref_char;
        }
        if (!compound_genotypes) {
          for (uint32_t geno1_code = 0; geno1_code != 4; ++geno1_code) {
            for (uint32_t geno0_code = 0; geno0_code != 4; ++geno0_code) {
              geno_pair_matrix_iter = memcpya_k(geno_pair_matrix_iter, geno_strs[0][geno0_code], 4);
              geno_pair_matrix_iter = memcpya_k(geno_pair_matrix_iter, geno_strs[1][geno1_code], 4);
            }
          }
        } else {
          for (uint32_t geno1_code = 0; geno1_code != 4; ++geno1_code) {
            for (uint32_t geno0_code = 0; geno0_code != 4; ++geno0_code) {
              geno_pair_matrix_iter = memcpya_k(geno_pair_matrix_iter, geno_strs[0][geno0_code], 3);
              geno_pair_matrix_iter = memcpya_k(geno_pair_matrix_iter, geno_strs[1][geno1_code], 3);
            }
          }
        }
      }
      geno_pair_matrix = geno_pair_matrix_start;
    } else {
      if (unlikely(compound_genotypes)) {
        logputs("\n");
        logerrprintfww("Error: %s cannot contain multi-character allele codes.\n", outname);
        goto ExportPed_ret_INCONSISTENT_INPUT;
      }
      if (unlikely(bigstack_alloc_kcp(variant_ct * (4 * k1LU), &genotype_strs) ||
                   // (3n+1) strings, each string is null-terminated and has
                   // length >= 4
                   (bigstack_left() < 5 + (15 * k1LU) * variant_ct))) {
        goto ExportPed_ret_NOMEM;
      }
      unsigned char* tmp_alloc_base = g_bigstack_base;
      unsigned char* tmp_alloc_end = &(bigstack_end_mark[-5]);
      tmp_alloc_end[0] = exportf_delim;
      tmp_alloc_end[1] = lomg_char;
      tmp_alloc_end[2] = exportf_delim;
      tmp_alloc_end[3] = lomg_char;
      tmp_alloc_end[4] = '\0';
      const char* missing_geno_str = R_CAST(const char*, tmp_alloc_end);

      const char** genotype_strs_iter = genotype_strs;
      uintptr_t variant_uidx_base = 0;
      uintptr_t variant_include_bits = variant_include[0];
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          // Unlike .bim export, .map export does NOT verify that all remaining
          // variants are biallelic.
          if (unlikely(allele_idx_offsets[variant_uidx + 1] != allele_idx_offset_base + 2)) {
            logputs("\n");
            logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
            goto ExportPed_ret_INCONSISTENT_INPUT;
          }
        }
        const char* ref_allele = allele_storage[allele_idx_offset_base];
        const char* alt_allele = allele_storage[allele_idx_offset_base + 1];
        const uint32_t ref_slen = strlen(ref_allele);
        const uint32_t alt_slen = strlen(alt_allele);
        // Note that this is based on PLINK 1 .bed, not PLINK 2, 2-bit
        // encoding.  i.e. 0 = alt-alt, 2 = alt-ref (rendered in that order
        // in the .ped).
        if (PtrWSubCk(tmp_alloc_base, 3 * (ref_slen + alt_slen) + 9, &tmp_alloc_end)) {
          goto ExportPed_ret_NOMEM;
        }
        char* geno0_start = R_CAST(char*, tmp_alloc_end);
        char* geno0_iter = geno0_start;
        *geno0_iter++ = exportf_delim;
        geno0_iter = memcpyax(geno0_iter, alt_allele, alt_slen, exportf_delim);
        char* geno2_start = memcpya(geno0_iter, alt_allele, alt_slen + 1);
        *genotype_strs_iter++ = geno0_start;

        *genotype_strs_iter++ = missing_geno_str;

        char* geno2_iter = geno2_start;
        *geno2_iter++ = exportf_delim;
        geno2_iter = memcpyax(geno2_iter, alt_allele, alt_slen, exportf_delim);
        char* geno3_start = memcpya(geno2_iter, ref_allele, ref_slen + 1);
        *genotype_strs_iter++ = geno2_start;

        char* geno3_iter = geno3_start;
        *geno3_iter++ = exportf_delim;
        geno3_iter = memcpyax(geno3_iter, ref_allele, ref_slen, exportf_delim);
        memcpy(geno3_iter, ref_allele, ref_slen + 1);
        *genotype_strs_iter++ = geno3_start;
      }
      BigstackEndSet(tmp_alloc_end);
      if (max_allele_slen >= kMaxMediumLine / kBitsPerWord) {
        // We only check writebuf_flush after each block of kBitsPerWordD2
        // genotypes.
        if (unlikely(bigstack_alloc_c(kMaxMediumLine + (1 + max_allele_slen) * kBitsPerWord, &writebuf))) {
          goto ExportPed_ret_NOMEM;
        }
      }
    }

    char* write_iter = writebuf;
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);

    uintptr_t* pheno_nm = nullptr;
    uintptr_t* pheno_cc = nullptr;
    double* pheno_qt = nullptr;
    // .ped files don't support categorical phenotypes
    const uint32_t pheno_idx = FirstCcOrQtPhenoIdx(pheno_cols, pheno_ct);
    if (pheno_idx != UINT32_MAX) {
      const PhenoDtype type_code = pheno_cols[pheno_idx].type_code;
      pheno_nm = pheno_cols[pheno_idx].nonmiss;
      if (type_code == kPhenoDtypeCc) {
        pheno_cc = pheno_cols[pheno_idx].data.cc;
      } else {
        pheno_qt = pheno_cols[pheno_idx].data.qt;
      }
    }
    const uint32_t lomp_slen = strlen(legacy_output_missing_pheno);

    // * Read phase: main thread reads raw bytes with MultireadNonempty(),
    //   while other thread(s) decode to standard 2-bit genovecs
    // * Write phase: one thread transposes+rotates to sample-major .bed, and
    //   one thread renders and writes the final .ped text.  (Strictly
    //   speaking, the rotation is unnecessary, but its cost is negligible and
    //   we want to reuse TransposeToPlink1SmajWriteThread.)

    // todo: check when this saturates
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    uint32_t read_block_size;

    // note that this is restricted to half of available workspace
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct, variant_ct, bigstack_left() / 2, pgr_alloc_cacheline_ct, 0, 0, 0, pgfip, &calc_thread_ct, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &read_ctx.pgr_ptrs, &read_ctx.variant_uidx_starts))) {
      goto ExportPed_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &read_tg))) {
      goto ExportPed_ret_NOMEM;
    }
    read_ctx.variant_include = variant_include;
    read_ctx.refalt1_select = nullptr;
    read_ctx.err_info = (~0LLU) << 32;
    SetThreadFuncAndData(TransposeToSmajReadThread, &read_ctx, &read_tg);

    const uintptr_t variant_cacheline_ct = NypCtToCachelineCt(variant_ct);
    // transpose_calc_thread_ct == 1, since transposition is much cheaper than
    // the final render+write step
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    if (unlikely(SetThreadCt(1, &transpose_tg) ||
                 bigstack_alloc_w(raw_sample_ctl, &sample_include) ||
                 bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_vp(1, &transpose_ctx.thread_vecaligned_bufs))) {
      goto ExportPed_ret_NOMEM;
    }
    transpose_ctx.thread_vecaligned_bufs[0] = S_CAST(VecW*, bigstack_alloc_raw(kPglNypTransposeBufbytes));
    // each of the two transpose-result buffers should use <= 1/8 of the
    // remaining workspace
    const uintptr_t indmaj_bed_cachelines_avail = bigstack_left() / (kCacheline * 8);
    uint32_t sample_batch_size = kPglNypTransposeBatch;
    if (variant_cacheline_ct * kPglNypTransposeBatch > indmaj_bed_cachelines_avail) {
      sample_batch_size = RoundDownPow2(indmaj_bed_cachelines_avail / variant_cacheline_ct, kBitsPerWordD2);
      if (unlikely(!sample_batch_size)) {
        goto ExportPed_ret_NOMEM;
      }
    }
    transpose_ctx.smaj_writebufs[0] = S_CAST(uintptr_t*, bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size));
    transpose_ctx.smaj_writebufs[1] = S_CAST(uintptr_t*, bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size));

    const uintptr_t readbuf_vecs_avail = (bigstack_left() / kCacheline) * kVecsPerCacheline;
    if (unlikely(readbuf_vecs_avail < variant_ct)) {
      goto ExportPed_ret_NOMEM;
    }
    uintptr_t read_sample_ctv2 = readbuf_vecs_avail / variant_ct;
    uint32_t read_sample_ct;
    if (read_sample_ctv2 >= NypCtToVecCt(sample_ct)) {
      read_sample_ct = sample_ct;
    } else {
      read_sample_ct = read_sample_ctv2 * kNypsPerVec;
    }
    uintptr_t read_sample_ctaw2 = NypCtToAlignedWordCt(read_sample_ct);
    uintptr_t* vmaj_readbuf = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(variant_ct * read_sample_ctaw2 * kBytesPerWord));
    read_ctx.vmaj_readbuf = vmaj_readbuf;
    transpose_ctx.variant_include = variant_include;
    transpose_ctx.variant_ct = variant_ct;
    transpose_ctx.vmaj_readbuf = vmaj_readbuf;
    SetThreadFuncAndData(TransposeToPlink1SmajWriteThread, &transpose_ctx, &transpose_tg);

    const char* sample_ids = piip->sii.sample_ids;
    const char* paternal_ids = piip->parental_id_info.paternal_ids;
    const char* maternal_ids = piip->parental_id_info.maternal_ids;
    const uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;

    uint32_t sample_uidx_start = AdvTo1Bit(orig_sample_include, 0);
    const uint32_t variant_ctl2 = NypCtToWordCt(variant_ct);
    const uintptr_t variant_ctaclw2 = variant_cacheline_ct * kWordsPerCacheline;
    const uint32_t final_backtrack_byte_ct = (4 - compound_genotypes) * ((-variant_ct) & (kBitsPerWordD2 - 1));
    const uint32_t pass_ct = 1 + (sample_ct - 1) / read_sample_ct;
    for (uint32_t pass_idx = 0; pass_idx != pass_ct; ++pass_idx) {
      memcpy(sample_include, orig_sample_include, raw_sample_ctl * sizeof(intptr_t));
      if (sample_uidx_start) {
        ClearBitsNz(0, sample_uidx_start, sample_include);
      }
      uint32_t sample_uidx_end;
      if (pass_idx + 1 == pass_ct) {
        read_sample_ct = sample_ct - pass_idx * read_sample_ct;
        read_sample_ctaw2 = NypCtToAlignedWordCt(read_sample_ct);
        sample_uidx_end = raw_sample_ct;
      } else {
        sample_uidx_end = FindNth1BitFrom(orig_sample_include, sample_uidx_start + 1, read_sample_ct);
        ClearBitsNz(sample_uidx_end, raw_sample_ct, sample_include);
      }
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      read_ctx.sample_include = sample_include;
      read_ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
      read_ctx.sample_ct = read_sample_ct;
      transpose_ctx.sample_ct = read_sample_ct;
      if (pass_idx) {
        pgfip->block_base = main_loadbufs[0];
        // er, don't need SetBaseAndOffset0?
        PgrSetBaseAndOffset0(main_loadbufs[0], calc_thread_ct, read_ctx.pgr_ptrs);
      }
      uint32_t parity = 0;
      uint32_t read_block_idx = 0;
      ReinitThreads(&read_tg);
      uint32_t pct = 0;
      uint32_t next_print_idx = variant_ct / 100;
      putc_unlocked('\r', stdout);
      printf("--export %s pass %u/%u: loading... 0%%", compound_genotypes? "compound-genotypes" : "ped", pass_idx + 1, pass_ct);
      fflush(stdout);
      for (uint32_t variant_idx = 0; ; ) {
        const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &read_tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
        if (unlikely(reterr)) {
          goto ExportPed_ret_PGR_FAIL;
        }
        if (variant_idx) {
          JoinThreads(&read_tg);
          reterr = S_CAST(PglErr, read_ctx.err_info);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, read_ctx.err_info >> 32);
            goto ExportPed_ret_1;
          }
        }
        if (!IsLastBlock(&read_tg)) {
          read_ctx.cur_block_write_ct = cur_block_write_ct;
          ComputeUidxStartPartition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, read_ctx.variant_uidx_starts);
          PgrCopyBaseAndOffset(pgfip, calc_thread_ct, read_ctx.pgr_ptrs);
          if (variant_idx + cur_block_write_ct == variant_ct) {
            DeclareLastThreadBlock(&read_tg);
          }
          if (unlikely(SpawnThreads(&read_tg))) {
            goto ExportPed_ret_THREAD_CREATE_FAIL;
          }
        }
        parity = 1 - parity;
        if (variant_idx == variant_ct) {
          break;
        }
        if (variant_idx >= next_print_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (variant_idx * 100LLU) / variant_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
        }

        ++read_block_idx;
        variant_idx += cur_block_write_ct;
        pgfip->block_base = main_loadbufs[parity];
      }
      // 2. Transpose and write.
      ReinitThreads(&transpose_tg);
      transpose_ctx.sample_batch_size = sample_batch_size;
      parity = 0;
      if (pct > 10) {
        fputs("\b \b", stdout);
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\bwriting... 0%", stdout);
      fflush(stdout);
      pct = 0;
      uintptr_t sample_uidx_base;
      uintptr_t sample_include_bits;
      BitIter1Start(sample_include, sample_uidx_start, &sample_uidx_base, &sample_include_bits);
      uint32_t flush_sample_idx = 0;
      next_print_idx = read_sample_ct / 100;
      for (uint32_t flush_sample_idx_end = 0; ; ) {
        if (!IsLastBlock(&transpose_tg)) {
          if (flush_sample_idx_end + sample_batch_size >= read_sample_ct) {
            DeclareLastThreadBlock(&transpose_tg);
            transpose_ctx.sample_batch_size = read_sample_ct - flush_sample_idx_end;
          }
          if (unlikely(SpawnThreads(&transpose_tg))) {
            goto ExportPed_ret_THREAD_CREATE_FAIL;
          }
        }
        if (flush_sample_idx_end) {
          uintptr_t* indmaj_bed_iter = transpose_ctx.smaj_writebufs[1 - parity];
          for (; flush_sample_idx != flush_sample_idx_end; ++flush_sample_idx) {
            char* line_uncounted_start = write_iter;
            const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
            const char* cur_sample_fid = &(sample_ids[sample_uidx * max_sample_id_blen]);
            const char* fid_end = AdvToDelim(cur_sample_fid, '\t');
            write_iter = memcpyax(write_iter, cur_sample_fid, fid_end - cur_sample_fid, exportf_delim);
            write_iter = strcpyax(write_iter, &(fid_end[1]), exportf_delim);
            write_iter = strcpyax(write_iter, &(paternal_ids[sample_uidx * max_paternal_id_blen]), exportf_delim);
            write_iter = strcpyax(write_iter, &(maternal_ids[sample_uidx * max_maternal_id_blen]), exportf_delim);
            *write_iter++ = Sexchar(sex_nm, sex_male, sample_uidx);
            *write_iter++ = exportf_delim;
            if ((!pheno_nm) || (!IsSet(pheno_nm, sample_uidx))) {
              write_iter = memcpya(write_iter, legacy_output_missing_pheno, lomp_slen);
            } else if (pheno_cc) {
              *write_iter++ = '1' + IsSet(pheno_cc, sample_uidx);
            } else {
              write_iter = dtoa_g(pheno_qt[sample_uidx], write_iter);
            }
            if (geno_pair_matrix) {
              const uint32_t render_block_wsize = kMaxMediumLine / (4 * kBitsPerWordD2);
              const char* geno_pair_matrix_iter = geno_pair_matrix;
              for (uint32_t widx = 0; ; ) {
                uint32_t widx_stop = widx + render_block_wsize;
                if (widx_stop > variant_ctl2) {
                  if (widx == variant_ctl2) {
                    break;
                  }
                  widx_stop = variant_ctl2;
                }
                if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
                  goto ExportPed_ret_WRITE_FAIL;
                }
                if (!compound_genotypes) {
                  for (; widx != widx_stop; ++widx) {
                    uintptr_t geno_word = indmaj_bed_iter[widx];
                    for (uint32_t pair_idx = 0; pair_idx != kBitsPerWordD4; ++pair_idx) {
                      const uintptr_t geno_pair = geno_word & 15;
                      write_iter = memcpya_k(write_iter, &(geno_pair_matrix_iter[geno_pair * 8]), 8);
                      geno_pair_matrix_iter = &(geno_pair_matrix_iter[128]);
                      geno_word >>= 4;
                    }
                  }
                } else {
                  for (; widx != widx_stop; ++widx) {
                    uintptr_t geno_word = indmaj_bed_iter[widx];
                    for (uintptr_t pair_idx = 0; pair_idx != kBitsPerWordD4; ++pair_idx) {
                      const uintptr_t geno_pair = geno_word & 15;
                      write_iter = memcpya_k(write_iter, &(geno_pair_matrix_iter[geno_pair * 6]), 6);
                      geno_pair_matrix_iter = &(geno_pair_matrix_iter[96]);
                      geno_word >>= 4;
                    }
                  }
                }
              }
              write_iter -= final_backtrack_byte_ct;
              if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
                goto ExportPed_ret_WRITE_FAIL;
              }
              AppendBinaryEoln(&write_iter);
            } else {
              uint64_t line_blen = write_iter - line_uncounted_start;
              uint32_t variant_idx = 0;
              for (uint32_t widx = 0; ; ++widx) {
                uint32_t variant_idx_stop = variant_idx + kBitsPerWordD2;
                if (variant_idx_stop > variant_ct) {
                  if (variant_idx == variant_ct) {
                    break;
                  }
                  variant_idx_stop = variant_ct;
                }

                line_blen += write_iter - line_uncounted_start;
                if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
                  goto ExportPed_ret_WRITE_FAIL;
                }
                line_uncounted_start = write_iter;

                uintptr_t geno_word = indmaj_bed_iter[widx];
                for (; variant_idx != variant_idx_stop; ++variant_idx) {
                  const uintptr_t cur_geno = geno_word & 3;
                  // no overflow risk, since we errored out earlier if
                  // variant_ct * 4 > kMaxLongLine - small constant
                  write_iter = strcpya(write_iter, genotype_strs[variant_idx * 4 + cur_geno]);
                  geno_word >>= 2;
                }
              }
              line_blen += write_iter - line_uncounted_start;
              if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
                goto ExportPed_ret_WRITE_FAIL;
              }
              AppendBinaryEoln(&write_iter);
              line_blen += strlen(EOLN_STR);
              if (unlikely(line_blen > kMaxLongLine)) {
                logputs("\n");
                logerrputs("Error: --export ped would create an excessively long line.\n");
                goto ExportPed_ret_INCONSISTENT_INPUT;
              }
            }
            indmaj_bed_iter = &(indmaj_bed_iter[variant_ctaclw2]);
          }
          if (flush_sample_idx_end == read_sample_ct) {
            break;
          }
          if (flush_sample_idx_end >= next_print_idx) {
            if (pct > 10) {
              putc_unlocked('\b', stdout);
            }
            pct = (flush_sample_idx_end * 100LLU) / read_sample_ct;
            printf("\b\b%u%%", pct++);
            fflush(stdout);
            next_print_idx = (pct * S_CAST(uint64_t, read_sample_ct)) / 100;
          }
        }
        JoinThreads(&transpose_tg);
        parity = 1 - parity;
        flush_sample_idx_end += sample_batch_size;
        if (flush_sample_idx_end > read_sample_ct) {
          flush_sample_idx_end = read_sample_ct;
        }
      }
      sample_uidx_start = sample_uidx_end;
      if (pct > 10) {
        fputs("\b \b", stdout);
      }
      sample_uidx_start = sample_uidx_end;
    }

    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto ExportPed_ret_WRITE_FAIL;
    }
    fputs("\b\bdone.\n", stdout);
    logprintfww("--export %s: %s written.\n", compound_genotypes? "compound-genotypes" : "ped", outname);
  }
  while (0) {
  ExportPed_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportPed_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportPed_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ExportPed_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportPed_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ExportPed_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 ExportPed_ret_1:
  CleanupThreads(&transpose_tg);
  CleanupThreads(&read_tg);
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
