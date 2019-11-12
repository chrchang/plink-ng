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
#include "plink2_export.h"

#include <time.h>

#ifdef __cplusplus
namespace plink2 {
#endif

void InitExportf(ExportfInfo* exportf_info_ptr) {
  exportf_info_ptr->flags = kfExportf0;
  exportf_info_ptr->idpaste_flags = kfIdpaste0;
  exportf_info_ptr->id_delim = '\0';
  exportf_info_ptr->bgen_bits = 0;
  exportf_info_ptr->vcf_mode = kVcfExport0;
  exportf_info_ptr->export_allele_fname = nullptr;
}

void CleanupExportf(ExportfInfo* exportf_info_ptr) {
  free_cond(exportf_info_ptr->export_allele_fname);
}

PglErr ExportAlleleLoad(const char* fname, const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t max_thread_ct, uint32_t* max_export_allele_slen_ptr, STD_ARRAY_PTR_DECL(AlleleCode, 2, export_allele), const char*** allele_missing_ptr) {
  // "Permanent" allocations occur on high end of bigstack.
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  TextStream txs;
  PreinitTextStream(&txs);
  {
    uintptr_t* already_seen;
    if (unlikely(bigstack_calloc_w(BitCtToWordCt(raw_variant_ct), &already_seen))) {
      goto ExportAlleleLoad_ret_NOMEM;
    }
    const uint32_t max_allele_slen = *max_export_allele_slen_ptr;
    uint32_t max_export_allele_slen = max_allele_slen;
    const uintptr_t linebuf_min = 64 + kMaxIdSlen + max_export_allele_slen;
    uint32_t max_line_blen;
    if (StandardizeMaxLineBlenEx(MAXV(bigstack_left() / 8, linebuf_min), linebuf_min, &max_line_blen)) {
      goto ExportAlleleLoad_ret_NOMEM;
    }
    reterr = InitTextStream(fname, max_line_blen, MAXV(max_thread_ct - 1, 1), &txs);
    if (unlikely(reterr)) {
      goto ExportAlleleLoad_ret_TSTREAM_FAIL;
    }

    // use minimum instead of fast hash table size, since allele_missing could
    // be very large
    // (if this path was hotter, we could check if the input file was
    // uncompressed and small enough)
    const uint32_t variant_id_htable_size = GetHtableMinSize(variant_ct);
    uint32_t* variant_id_htable;
    if (unlikely(bigstack_alloc_u32(variant_id_htable_size * sizeof(int32_t), &variant_id_htable))) {
      goto ExportAlleleLoad_ret_NOMEM;
    }
    reterr = PopulateIdHtableMt(g_bigstack_end, variant_include, variant_ids, variant_ct, 0, variant_id_htable_size, max_thread_ct, &g_bigstack_base, variant_id_htable, nullptr);
    if (unlikely(reterr)) {
      goto ExportAlleleLoad_ret_1;
    }

    const char input_missing_geno_char = *g_input_missing_geno_ptr;
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = g_bigstack_end;
    const char** allele_missing = nullptr;
    uintptr_t duplicate_ct = 0;
    uintptr_t miss_ct = 0;
    uintptr_t line_idx = 0;
    uint32_t hit_ct = 0;
    uint32_t cur_allele_ct = 2;
    while (1) {
      ++line_idx;
      char* line_start = TextGet(&txs);
      if (!line_start) {
        if (likely(!TextStreamErrcode2(&txs, &reterr))) {
          break;
        }
        goto ExportAlleleLoad_ret_TSTREAM_FAIL;
      }
      char* token_end = CurTokenEnd(line_start);
      const uint32_t variant_uidx = VariantIdDupflagHtableFind(line_start, variant_ids, variant_id_htable, token_end - line_start, variant_id_htable_size, max_variant_id_slen);
      if (variant_uidx & 0x80000000U) {
        if (likely(variant_uidx == UINT32_MAX)) {
          ++miss_ct;
          continue;
        }
        *token_end = '\0';
        snprintf(g_logbuf, kLogbufSize, "Error: --export-allele variant ID '%s' appears multiple times in dataset.\n", line_start);
        goto ExportAlleleLoad_ret_INCONSISTENT_INPUT_WW;
      }
      char* allele_start = FirstNonTspace(token_end);
      char* allele_end = FirstSpaceOrEoln(allele_start);
      const uint32_t allele_slen = allele_end - allele_start;
      if (!allele_slen) {
        logerrprintfww("Error: Line %" PRIuPTR " of --export-allele file has fewer tokens than expected.\n", line_idx);
        goto ExportAlleleLoad_ret_MALFORMED_INPUT;
      }
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      int32_t cur_allele_idx = -1;
      if (allele_slen <= max_allele_slen) {
        for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
          if (memequal(allele_start, cur_alleles[allele_idx], allele_slen) && (!cur_alleles[allele_idx][allele_slen])) {
            cur_allele_idx = allele_idx;
            break;
          }
        }
      }
      if (IsSet(already_seen, variant_uidx)) {
        // Only ok if allele is the same as before.
        ++duplicate_ct;
        if (cur_allele_idx == -1) {
          if (likely(allele_missing && allele_missing[variant_uidx] && memequal(allele_start, allele_missing[variant_uidx], allele_slen) && (!allele_missing[variant_uidx][allele_slen]))) {
            continue;
          }
        } else {
          STD_ARRAY_KREF(AlleleCode, 2) cur_export_allele = export_allele[variant_uidx];
          if (likely(((!allele_missing) || (!allele_missing[variant_uidx])) && (S_CAST(uint32_t, cur_allele_idx) == cur_export_allele[0]))) {
            continue;
          }
        }
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --export-allele file has a variant ID duplicating an earlier line, and with a different allele code.\n", line_idx);
        goto ExportAlleleLoad_ret_INCONSISTENT_INPUT_WW;
      }
      SetBit(variant_uidx, already_seen);
      ++hit_ct;
      if (cur_allele_idx != -1) {
        STD_ARRAY_REF(AlleleCode, 2) cur_export_allele = export_allele[variant_uidx];
        cur_export_allele[0] = S_CAST(uint32_t, cur_allele_idx);
        // cur_export_allele[1] doesn't matter
        continue;
      }
      if (!allele_missing) {
        if (unlikely(bigstack_end_calloc_kcp(raw_variant_ct, &allele_missing))) {
          goto ExportAlleleLoad_ret_NOMEM;
        }
        tmp_alloc_end = g_bigstack_end;
      }
      if (allele_slen == 1) {
        char allele_char = *allele_start;
        if (allele_char == input_missing_geno_char) {
          allele_char = '.';
        }
        allele_missing[variant_uidx] = &(g_one_char_strs[2 * ctou32(allele_char)]);
        continue;
      }
      // no underflow danger since allele must fit in line-buffer, which is
      // positioned earlier in bigstack
      tmp_alloc_end -= allele_slen + 1;
      if (unlikely(tmp_alloc_end < tmp_alloc_base)) {
        goto ExportAlleleLoad_ret_NOMEM;
      }
      memcpyx(tmp_alloc_end, allele_start, allele_slen, '\0');
      allele_missing[variant_uidx] = R_CAST(char*, tmp_alloc_end);
    }
    if (duplicate_ct) {
      logerrprintfww("Warning: %" PRIuPTR " duplicate variant ID%s in --export-allele file (not an error since allele%s consistent).\n", duplicate_ct, (duplicate_ct == 1)? "" : "s", (duplicate_ct == 1)? " was" : "s were");
    }
    if (miss_ct) {
      snprintf(g_logbuf, kLogbufSize, "--export-allele: %u variant%s updated, %" PRIuPTR " variant ID%s not present.\n", hit_ct, (hit_ct == 1)? "" : "s", miss_ct, (miss_ct == 1)? "" : "s");
      WordWrapB(0);
    } else {
      snprintf(g_logbuf, kLogbufSize, "--export-allele: %u variant%s updated.\n", hit_ct, (hit_ct == 1)? "" : "s");
    }
    logputsb();
    *max_export_allele_slen_ptr = max_export_allele_slen;
    *allele_missing_ptr = allele_missing;
    BigstackEndSet(tmp_alloc_end);
  }
  while (0) {
 ExportAlleleLoad_ret_TSTREAM_FAIL:
    TextStreamErrPrint("--export-allele file", &txs);
    break;
 ExportAlleleLoad_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
 ExportAlleleLoad_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
 ExportAlleleLoad_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  }
 ExportAlleleLoad_ret_1:
  CleanupTextStream2("--export-allele file", &txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __arm__
#  error "Unaligned accesses in Export012Vmaj()."
#endif
PglErr Export012Vmaj(const char* outname, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const char* sample_ids, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(AlleleCode, 2, export_allele), const char* const* export_allele_missing, const double* variant_cms, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t variant_ct, uint32_t max_allele_slen, char exportf_delim, PgenReader* simple_pgrp) {
  // todo: --recode-allele
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t max_chr_blen = 1 + GetMaxChrSlen(cip);
    char* chr_buf;  // includes trailing tab
    char* writebuf;
    uintptr_t* genovec;
    // dosages are limited to 7 characters (x.yyyyy)
    if (unlikely(
            bigstack_alloc_c(max_chr_blen, &chr_buf) ||
            bigstack_alloc_c(kMaxMediumLine + max_chr_blen + 2 * kMaxIdSlen + 48 + 2 * max_allele_slen + (8 * k1LU) * sample_ct, &writebuf) ||
            bigstack_alloc_w(sample_ctl2, &genovec))) {
      goto Export012Vmaj_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    const uint32_t dosage_is_present = PgrGetGflags(simple_pgrp) & kfPgenGlobalDosagePresent;
    uintptr_t* dosage_present = nullptr;
    Dosage* dosage_main = nullptr;
    if (dosage_is_present) {
      if (unlikely(
              bigstack_alloc_w(sample_ctl, &dosage_present) ||
              bigstack_alloc_dosage(sample_ct, &dosage_main))) {
        goto Export012Vmaj_ret_NOMEM;
      }
    }
    uintptr_t* missingness_dosage = nullptr;
    if (export_allele_missing) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &missingness_dosage))) {
        goto Export012Vmaj_ret_NOMEM;
      }
    }
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto Export012Vmaj_ret_OPEN_FAIL;
    }
    char* write_iter;
    if (exportf_delim == '\t') {
      write_iter = strcpya_k(writebuf, "CHR\tSNP\t(C)M\tPOS\tCOUNTED\tALT");
    } else {
      write_iter = strcpya_k(writebuf, "CHR SNP (C)M POS COUNTED ALT");
    }
    uintptr_t sample_uidx_base = 0;
    uintptr_t sample_include_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &sample_include_bits);
      *write_iter++ = exportf_delim;
      const char* fid_start = &(sample_ids[sample_uidx * max_sample_id_blen]);
      const char* fid_end = AdvToDelim(fid_start, exportf_delim);
      write_iter = memcpyax(write_iter, fid_start, fid_end - fid_start, '_');
      write_iter = strcpya(write_iter, &(fid_end[1]));
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto Export012Vmaj_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);
    logprintfww5("--export A-transpose to %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_blen = 0;
    uint32_t exported_allele_idx = 0;
    uint32_t cur_allele_ct = 2;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;

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
      write_iter = memcpya(write_iter, chr_buf, chr_blen);
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], exportf_delim);
      if (variant_cms) {
        write_iter = dtoa_g(variant_cms[variant_uidx], write_iter);
      } else {
        *write_iter++ = '0';
      }
      *write_iter++ = exportf_delim;
      write_iter = u32toa_x(variant_bps[variant_uidx], exportf_delim, write_iter);

      if (export_allele) {
        exported_allele_idx = export_allele[variant_uidx][0];
        // alt1_allele_idx doesn't matter here
      }
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (export_allele_missing && export_allele_missing[variant_uidx]) {
        reterr = PgrGetMissingnessD(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, nullptr, missingness_dosage, nullptr, genovec);
        if (unlikely(reterr)) {
          goto Export012Vmaj_ret_PGR_FAIL;
        }
        write_iter = strcpyax(write_iter, export_allele_missing[variant_uidx], exportf_delim);
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto Export012Vmaj_ret_WRITE_FAIL;
        }
        // in this case, exported_allele_idx corresponds to the REF allele.
        // May as well put it first.
        write_iter = strcpya(write_iter, cur_alleles[exported_allele_idx]);
        for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
          if (allele_idx == exported_allele_idx) {
            continue;
          }
          *write_iter++ = ',';
          write_iter = strcpya(write_iter, cur_alleles[allele_idx]);
          if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
            goto Export012Vmaj_ret_WRITE_FAIL;
          }
        }
        // We only care about missingness here.
        const uint32_t sample_ctl_m1 = sample_ctl - 1;
        uint32_t loop_len = kBitsPerWord;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl_m1) {
            if (widx > sample_ctl_m1) {
              break;
            }
            loop_len = ModNz(sample_ct, kBitsPerWord);
          }
          uintptr_t missingness_word = missingness_dosage[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
            *write_iter++ = exportf_delim;
            if (missingness_word & 1) {
              write_iter = strcpya_k(write_iter, "NA");
            } else {
              *write_iter++ = '0';
            }
            missingness_word >>= 1;
          }
        }
      } else {
        uint32_t dosage_ct;
        reterr = PgrGet1D(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, exported_allele_idx, simple_pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
        if (unlikely(reterr)) {
          goto Export012Vmaj_ret_PGR_FAIL;
        }
        write_iter = strcpyax(write_iter, cur_alleles[exported_allele_idx], exportf_delim);
        const uint32_t first_alt_idx = (exported_allele_idx == 0);
        write_iter = strcpya(write_iter, cur_alleles[first_alt_idx]);
        if (cur_allele_ct > 2) {
          for (uint32_t allele_idx = first_alt_idx + 1; allele_idx != cur_allele_ct; ++allele_idx) {
            if (allele_idx == exported_allele_idx) {
              continue;
            }
            if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
              goto Export012Vmaj_ret_WRITE_FAIL;
            }
            *write_iter++ = ',';
            write_iter = strcpya(write_iter, cur_alleles[allele_idx]);
          }
        }
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto Export012Vmaj_ret_WRITE_FAIL;
        }
        uint32_t loop_len = kBitsPerWordD2;
        if (!dosage_ct) {
          for (uint32_t widx = 0; ; ++widx) {
            if (widx >= sample_ctl2_m1) {
              if (widx > sample_ctl2_m1) {
                break;
              }
              loop_len = ModNz(sample_ct, kBitsPerWordD2);
            }
            uintptr_t geno_word = genovec[widx];
            for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
              *write_iter++ = exportf_delim;
              uintptr_t cur_geno = geno_word & 3;
              if (cur_geno != 3) {
                *write_iter++ = '0' + cur_geno;
              } else {
                write_iter = strcpya_k(write_iter, "NA");
              }
              geno_word >>= 2;
            }
          }
        } else {
          Dosage* dosage_main_iter = dosage_main;
          for (uint32_t widx = 0; ; ++widx) {
            if (widx >= sample_ctl2_m1) {
              if (widx > sample_ctl2_m1) {
                break;
              }
              loop_len = ModNz(sample_ct, kBitsPerWordD2);
            }
            uintptr_t geno_word = genovec[widx];
            uint32_t dosage_present_hw = R_CAST(Halfword*, dosage_present)[widx];
            for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
              *write_iter++ = exportf_delim;
              if (dosage_present_hw & 1) {
                write_iter = PrintSmallDosage(*dosage_main_iter++, write_iter);
              } else {
                uintptr_t cur_geno = geno_word & 3;
                if (cur_geno != 3) {
                  *write_iter++ = '0' + cur_geno;
                } else {
                  write_iter = strcpya_k(write_iter, "NA");
                }
              }
              geno_word >>= 2;
              dosage_present_hw >>= 1;
            }
          }
        }
      }
      AppendBinaryEoln(&write_iter);
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
      goto Export012Vmaj_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
  }
  while (0) {
  Export012Vmaj_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Export012Vmaj_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  Export012Vmaj_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  Export012Vmaj_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

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

  PglErr reterr;
} TransposeToSmajReadCtx;

THREAD_FUNC_DECL TransposeToSmajReadThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  TransposeToSmajReadCtx* ctx = S_CAST(TransposeToSmajReadCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select) = ctx->refalt1_select;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uintptr_t* sample_include = ctx->sample_include;
  const uint32_t* sample_include_cumulative_popcounts = ctx->sample_include_cumulative_popcounts;
  const uint32_t read_sample_ct = ctx->sample_ct;
  const uintptr_t read_sample_ctaw2 = NypCtToAlignedWordCt(read_sample_ct);
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  uintptr_t prev_copy_ct = 0;
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
      const PglErr reterr = PgrGet(sample_include, sample_include_cumulative_popcounts, read_sample_ct, variant_uidx, pgrp, vmaj_readbuf_iter);
      if (unlikely(reterr)) {
        // no synchronization needed for now since only kPglRetMalformedInput
        // possible?
        ctx->reterr = reterr;
        break;
      }
      if (refalt1_select && (refalt1_select[variant_uidx][0] == 1)) {
        GenovecInvertUnsafe(read_sample_ct, vmaj_readbuf_iter);
        // don't need ZeroTrailingNyps()
      }
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[read_sample_ctaw2]);
    }
    prev_copy_ct += cur_block_copy_ct;
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
    if (unlikely(fwrite_checked("l\x1b\0", 3, outfile))) {
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
      // read_ctx.allele_idx_offsets = allele_idx_offsets;
      read_ctx.refalt1_select = refalt1_select;
      read_ctx.reterr = kPglRetSuccess;
      SetThreadFuncAndData(TransposeToSmajReadThread, &read_ctx, &read_tg);

      const uintptr_t variant_cacheline_ct = NypCtToCachelineCt(variant_ct);
      uint32_t output_calc_thread_ct = MINV(calc_thread_ct, variant_cacheline_ct);
      // 4 still seems to be best in AVX2 case
      if (output_calc_thread_ct > 4) {
        output_calc_thread_ct = 4;
      }
      uintptr_t* sample_include;
      uint32_t* sample_include_cumulative_popcounts;
      if (unlikely(
              SetThreadCt(output_calc_thread_ct, &write_tg) ||
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
            reterr = read_ctx.reterr;
            if (unlikely(reterr)) {
              goto ExportIndMajorBed_ret_PGR_FAIL;
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
  CleanupThreads(&write_tg);
  CleanupThreads(&read_tg);
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  BigstackReset(bigstack_mark);
  return reterr;
}

static_assert(kDosageMid == 16384, "PrintGenDosage() needs to be updated.");
char* PrintGenDosage(uint32_t rawval, char* start) {
  // Similar to PrintSmallDosage(), but it's complicated by .gen import's
  // quantization step (instead of rounding the numbers directly, they're first
  // converted to bgen-1.1 equivalents).  We check
  //   ((n - 0.75)/16384, (n + 0.75)/16384) for even n
  //   ((n - 0.25)/16384, (n + 0.25)/16384) for odd n
  // due to banker's rounding.
  *start++ = '0' + (rawval / 16384);
  rawval = rawval % 16384;
  if (!rawval) {
    return start;
  }
  *start++ = '.';
  const uint32_t halfwidth_65536ths = 3 - (2 * (rawval % 2));
  // (rawval * 4) is in 65536ths
  // 65536 * 625 = 40960k

  const uint32_t range_top_40960k = (rawval * 4 + halfwidth_65536ths) * 625;
  // this is technically checking a half-open rather than a fully-open
  // interval, but that's fine since we never hit the boundary points
  if ((range_top_40960k % 4096) < 1250 * halfwidth_65536ths) {
    // when this is true, the four-decimal-place approximation is in the range
    // which round-trips back to our original number.
    const uint32_t four_decimal_places = range_top_40960k / 4096;
    return u32toa_trunc4(four_decimal_places, start);
  }

  // we wish to print (100000 * remainder + 8192) / 16384, left-0-padded.  and
  // may as well banker's round too.
  //
  // banker's rounding yields a different result than regular rounding for n/64
  // when n is congruent to 1 mod 4:
  //   1/64 = .015625 -> print 0.01562
  //   3/64 = .046875 -> print 0.04688
  //   5/64 = .078125 -> print 0.07812
  const uint32_t five_decimal_places = ((3125 * rawval + 256) / 512) - ((rawval % 1024) == 256);
  const uint32_t first_decimal_place = five_decimal_places / 10000;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_four_digits = five_decimal_places - first_decimal_place * 10000;
  if (last_four_digits) {
    return u32toa_trunc4(last_four_digits, start);
  }
  return start;
}

static inline BoolErr bgzfwrite_ck(char* buf_flush, BgzfCompressStream* bgzfp, char** write_iter_ptr) {
  if ((*write_iter_ptr) < buf_flush) {
    return 0;
  }
  char* buf = &(buf_flush[-kMaxMediumLine]);
  char* buf_end = *write_iter_ptr;
  *write_iter_ptr = buf;
  return BgzfWrite(buf, buf_end - buf, bgzfp);
}

BoolErr bgzfclose_flush(char* buf_flush, char* write_iter, BgzfCompressStream* bgzfp, PglErr* reterrp) {
  char* buf = &(buf_flush[-kMaxMediumLine]);

  // safe to ignore this error-return since CleanupBgzfCompressStream will also
  // error-return
  BgzfWrite(buf, write_iter - buf, bgzfp);

  return CleanupBgzfCompressStream(bgzfp, reterrp);
}

PglErr ExportOxGen(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uint32_t sample_ct, uint32_t variant_ct, uint32_t max_allele_slen, __maybe_unused uint32_t max_thread_ct, ExportfFlags exportf_flags, PgenReader* simple_pgrp, char* outname, char* outname_end, uint32_t* sample_missing_geno_cts) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  BgzfCompressStream bgzf;
  PreinitBgzfCompressStream(&bgzf);
  {
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    uintptr_t* genovec;
    uintptr_t* sex_male_collapsed;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (unlikely(
            bigstack_alloc_w(sample_ctl2, &genovec) ||
            bigstack_alloc_w(sample_ctl, &sex_male_collapsed))) {
      goto ExportOxGen_ret_NOMEM;
    }
    CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_male_collapsed);

    // See LoadSampleMissingCtsThread() in plink2_filter.cc.
    // Yes, this is overkill, but the obvious alternative of incrementing
    // sample_missing_geno_cts[] when writing a missing call requires a bit of
    // custom chrY logic anyway.
    const uint32_t acc1_vec_ct = BitCtToVecCt(sample_ct);
    uintptr_t* missing_acc1;
    if (unlikely(bigstack_calloc_w(acc1_vec_ct * kWordsPerVec * 45, &missing_acc1))) {
      goto ExportOxGen_ret_NOMEM;
    }
    const uint32_t acc4_vec_ct = acc1_vec_ct * 4;
    const uint32_t acc8_vec_ct = acc1_vec_ct * 8;
    VecW* missing_acc4 = &(R_CAST(VecW*, missing_acc1)[acc1_vec_ct]);
    VecW* missing_acc8 = &(missing_acc4[acc4_vec_ct]);
    VecW* missing_acc32 = &(missing_acc8[acc8_vec_ct]);

    uintptr_t* dosage_present = nullptr;
    Dosage* dosage_main = nullptr;
    if (PgrGetGflags(simple_pgrp) & kfPgenGlobalDosagePresent) {
      const uint32_t multiallelic_present = (allele_idx_offsets != nullptr);
      if (unlikely(
              bigstack_alloc_dosage(sample_ct * (1 + multiallelic_present), &dosage_main) ||
              bigstack_alloc_w(sample_ctl, &dosage_present))) {
        goto ExportOxGen_ret_NOMEM;
      }
    }
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    // if no dosages, all genotypes are 6 bytes (missing = " 0 0 0")
    // with dosages, we print up to 5 digits past the decimal point, so 7 bytes
    //   + space for each number, 24 bytes max
    const uintptr_t max_geno_slen = 6 + (dosage_present != nullptr) * 18;
    char* chr_buf;  // includes trailing space
    char* writebuf;
    if (unlikely(
            bigstack_alloc_c(max_chr_blen, &chr_buf) ||
            bigstack_alloc_c(kMaxMediumLine + max_chr_blen + kMaxIdSlen + 16 + 2 * max_allele_slen + max_geno_slen * sample_ct, &writebuf))) {
      goto ExportOxGen_ret_NOMEM;
    }
    {
      uint32_t clvl = 0;
      if (!(exportf_flags & kfExportfBgz)) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".gen");
      } else {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".gen.gz");
        clvl = kBgzfDefaultClvl;
      }
      reterr = InitBgzfCompressStreamEx(outname, 0, clvl, max_thread_ct, &bgzf);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto ExportOxGen_ret_1;
      }
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = writebuf;
    uint32_t chr_blen = 0;

    // although we don't support --set-hh-missing, etc. here, we do still want
    // to be aware of chrY so we can exclude nonmales from the
    // sample_missing_geno_cts update there.
    uint32_t is_y = 0;

    uint32_t chr_fo_idx = UINT32_MAX;
    const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
    uint32_t chr_end = 0;
    uint32_t vidx_rem15 = 15;
    uint32_t vidx_rem255d15 = 17;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    // " 1 0 0", " 0 1 0", " 0 0 1", " 0 0 0"
    const uint64_t hardcall_strs[4] = {0x302030203120LLU, 0x302031203020LLU, 0x312030203020LLU, 0x302030203020LLU};
    const uint32_t ref_allele_last = !(exportf_flags & kfExportfRefFirst);
    logprintfww5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t ref_allele_idx = 0;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
        // Oxford spec doesn't seem to require spaces for .gen (only .sample),
        // but in practice spaces always seem to be used, and plink 1.9 doesn't
        // let you toggle this, so let's not worry about supporting tabs here
        *chr_name_end++ = ' ';
        chr_blen = chr_name_end - chr_buf;
        is_y = (chr_idx == y_code);
      }
      write_iter = memcpya(write_iter, chr_buf, chr_blen);
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], ' ');
      write_iter = u32toa_x(variant_bps[variant_uidx], ' ', write_iter);
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        if (allele_idx_offsets[variant_uidx + 1] != allele_idx_offset_base + 2) {
          logputs("\n");
          logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
          goto ExportOxGen_ret_INCONSISTENT_INPUT;
        }
      }
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx][0];
      }
      // No dosage rescaling here, too messy to put that in more than one
      // place.
      uint32_t dosage_ct;
      reterr = PgrGetD(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
      if (unlikely(reterr)) {
        goto ExportOxGen_ret_PGR_FAIL;
      }
      if (ref_allele_idx != ref_allele_last) {
        GenovecInvertUnsafe(sample_ct, genovec);
        BiallelicDosage16Invert(dosage_ct, dosage_main);
      }

      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (ref_allele_last) {
        write_iter = strcpyax(write_iter, cur_alleles[1 - ref_allele_idx], ' ');
        write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
      } else {
        write_iter = strcpyax(write_iter, cur_alleles[ref_allele_idx], ' ');
        write_iter = strcpya(write_iter, cur_alleles[1 - ref_allele_idx]);
      }
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      if (!dosage_ct) {
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t geno_word = genovec[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            // could process 2 genotypes at a time and do 8 byte + 4 byte copy
            // for higher speed (may need to back up write_iter by 6 bytes at
            // end).  but this is not a bottleneck.
            memcpy(write_iter, &(hardcall_strs[geno_word & 3]), 8);
            write_iter = &(write_iter[6]);
            geno_word >>= 2;
          }
        }
      } else {
        const Halfword* dosage_present_alias = R_CAST(Halfword*, dosage_present);
        const Dosage* dosage_main_iter = dosage_main;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t geno_word = genovec[widx];
          uint32_t dosage_present_hw = dosage_present_alias[widx];
          if (!dosage_present_hw) {
            for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
              memcpy(write_iter, &(hardcall_strs[geno_word & 3]), 8);
              write_iter = &(write_iter[6]);
              geno_word >>= 2;
            }
          } else {
            for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
              if (dosage_present_hw & 1) {
                const uint32_t dosage_int = *dosage_main_iter++;
                if (dosage_int <= kDosageMid) {
                  *write_iter++ = ' ';
                  write_iter = PrintGenDosage(kDosageMid - dosage_int, write_iter);
                  *write_iter++ = ' ';
                  write_iter = PrintGenDosage(dosage_int, write_iter);
                  write_iter = strcpya_k(write_iter, " 0");
                } else {
                  assert(dosage_int <= kDosageMax);
                  write_iter = strcpya_k(write_iter, " 0 ");
                  write_iter = PrintGenDosage(kDosageMax - dosage_int, write_iter);
                  *write_iter++ = ' ';
                  write_iter = PrintGenDosage(dosage_int - kDosageMid, write_iter);
                }
              } else {
                memcpy(write_iter, &(hardcall_strs[geno_word & 3]), 8);
                write_iter = &(write_iter[6]);
              }
              geno_word >>= 2;
              dosage_present_hw >>= 1;
            }
          }
        }
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
        goto ExportOxGen_ret_WRITE_FAIL;
      }
      // bugfix (13 Apr 2018): this missingness calculation was only taking
      // hardcalls into account, which is inappropriate for .gen/.bgen
      GenoarrToMissingnessUnsafe(genovec, sample_ct, missing_acc1);
      if (dosage_ct) {
        BitvecInvmask(dosage_present, sample_ctl, missing_acc1);
      }
      if (is_y) {
        BitvecAnd(sex_male_collapsed, sample_ctl, missing_acc1);
      }
      VcountIncr1To4(missing_acc1, acc1_vec_ct, missing_acc4);
      if (!(--vidx_rem15)) {
        Vcount0Incr4To8(acc4_vec_ct, missing_acc4, missing_acc8);
        vidx_rem15 = 15;
        if (!(--vidx_rem255d15)) {
          Vcount0Incr8To32(acc8_vec_ct, missing_acc8, missing_acc32);
          vidx_rem255d15 = 17;
        }
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
    if (unlikely(bgzfclose_flush(writebuf_flush, write_iter, &bgzf, &reterr))) {
      goto ExportOxGen_ret_1;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    VcountIncr4To8(missing_acc4, acc4_vec_ct, missing_acc8);
    VcountIncr8To32(missing_acc8, acc8_vec_ct, missing_acc32);
    uint32_t* scrambled_missing_cts = R_CAST(uint32_t*, missing_acc32);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = VcountScramble1(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  ExportOxGen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportOxGen_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ExportOxGen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportOxGen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 ExportOxGen_ret_1:
  CleanupBgzfCompressStream(&bgzf, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __arm__
#  error "Unaligned accesses in ExportOxHapslegend()."
#endif
PglErr ExportOxHapslegend(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male_collapsed, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, __maybe_unused uint32_t max_thread_ct, ExportfFlags exportf_flags, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  assert(sample_ct);
  assert(variant_ct);
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  BgzfCompressStream bgzf;
  PreinitBgzfCompressStream(&bgzf);
  {
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uint32_t just_haps = (exportf_flags / kfExportfHaps) & 1;
    const uint32_t male_ct = PopcountWords(sex_male_collapsed, sample_ctl);
    if (unlikely(XymtIsNonempty(variant_include, cip, kChrOffsetY) && (male_ct != sample_ct))) {
      logerrprintf("Error: '--export haps%s' must exclude chrY unless the dataset is all-male.\n", just_haps? "" : "legend");
      goto ExportOxHapslegend_ret_INCONSISTENT_INPUT;
    }
    const uint32_t ref_allele_last = !(exportf_flags & kfExportfRefFirst);
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    char* chr_buf = nullptr;
    uint32_t is_x = 0;
    uint32_t is_haploid = 0;
    const uint32_t variant_uidx_start = AdvTo1Bit(variant_include, 0);
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t allele_idx0 = ref_allele_last;
    uint32_t allele_idx1 = 1 - ref_allele_last;
    uintptr_t writebuf_alloc = 0;
    if (!just_haps) {
      // .legend doesn't have a chromosome column, so verify we only need to
      // export a single chromosome
      chr_fo_idx = GetVariantChrFoIdx(cip, variant_uidx_start);
      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      if (unlikely((chr_end != raw_variant_ct) && (PopcountBitRange(variant_include, variant_uidx_start, chr_end) != variant_ct))) {
        logerrputs("Error: '--export hapslegend' does not support multiple chromosomes.\n");
        goto ExportOxHapslegend_ret_INCONSISTENT_INPUT;
      }
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      is_x = (chr_idx == x_code);
      is_haploid = IsSet(cip->haploid_mask, chr_idx);
      snprintf(outname_end, kMaxOutfnameExtBlen, ".legend");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
        goto ExportOxHapslegend_ret_OPEN_FAIL;
      }
      char* writebuf;
      if (unlikely(bigstack_alloc_c(kMaxMediumLine + kMaxIdSlen + 32 + 2 * max_allele_slen, &writebuf))) {
        goto ExportOxHapslegend_ret_NOMEM;
      }
      char* writebuf_flush = &(writebuf[kMaxMediumLine]);
      char* write_iter = strcpya_k(writebuf, "id position a0 a1" EOLN_STR);
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      uintptr_t variant_uidx_base;
      uintptr_t cur_bits;
      BitIter1Start(variant_include, variant_uidx_start, &variant_uidx_base, &cur_bits);
      for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
        const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
        write_iter = strcpyax(write_iter, variant_ids[variant_uidx], ' ');
        write_iter = u32toa_x(variant_bps[variant_uidx], ' ', write_iter);
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          if (unlikely((!refalt1_select) && (allele_idx_offsets[variant_uidx + 1] != allele_idx_offset_base + 2))) {
            logputs("\n");
            logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
            goto ExportOxHapslegend_ret_INCONSISTENT_INPUT;
          }
        }
        if (refalt1_select) {
          allele_idx0 = refalt1_select[variant_uidx][ref_allele_last];
          allele_idx1 = refalt1_select[variant_uidx][1 - ref_allele_last];
        }
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        write_iter = strcpyax(write_iter, cur_alleles[allele_idx0], ' ');
        write_iter = strcpya(write_iter, cur_alleles[allele_idx1]);
        AppendBinaryEoln(&write_iter);
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto ExportOxHapslegend_ret_WRITE_FAIL;
        }
      }
      if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
        goto ExportOxHapslegend_ret_WRITE_FAIL;
      }
      logputs("done.\n");
      BigstackReset(writebuf);
    } else {
      const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
      if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf))) {
        goto ExportOxHapslegend_ret_NOMEM;
      }
      writebuf_alloc = max_chr_blen + kMaxIdSlen + 32 + 2 * max_allele_slen;
    }
    writebuf_alloc += kMaxMediumLine + (4 * k1LU) * sample_ct + kCacheline;
    const uint32_t sample_ctv = BitCtToVecCt(sample_ct);
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    char* writebuf;
    uintptr_t* sex_male_collapsed_interleaved;
    uintptr_t* genovec;
    uintptr_t* phasepresent;
    uintptr_t* phaseinfo;
    if (unlikely(
            bigstack_alloc_w(sample_ctv * kWordsPerVec, &sex_male_collapsed_interleaved) ||
            bigstack_alloc_c(writebuf_alloc, &writebuf) ||
            bigstack_alloc_w(sample_ctl2, &genovec) ||
            bigstack_alloc_w(sample_ctl, &phasepresent) ||
            bigstack_alloc_w(sample_ctl, &phaseinfo))) {
      goto ExportOxHapslegend_ret_NOMEM;
    }
    // sex_male_collapsed had trailing bits zeroed out
    FillInterleavedMaskVec(sex_male_collapsed, sample_ctv, sex_male_collapsed_interleaved);
    // these could also be compile-time constants
    uint32_t genotext_haploid[32];
    // this can be more efficient, but don't worry about it for now
    genotext_haploid[0] = 0x202d2030;  // "0 - "
    genotext_haploid[2] = 0x21475542;  // "BUG!"
    genotext_haploid[4] = 0x202d2031;
    genotext_haploid[6] = 0x21475542;
    InitLookup16x4bx2(genotext_haploid);
    uint32_t genotext_diploid[112];
    genotext_diploid[0] = 0x20302030;
    genotext_diploid[2] = 0x21475542;
    genotext_diploid[4] = 0x20312031;
    genotext_diploid[6] = 0x21475542;
    genotext_diploid[34] = 0x20312030;
    genotext_diploid[38] = 0x20302031;
    InitPhaseLookup4b(genotext_diploid);
    uint32_t genotext_x[128];
    genotext_x[0] = 0x20302030;
    genotext_x[2] = 0x21475542;
    genotext_x[4] = 0x20312031;
    genotext_x[6] = 0x21475542;
    genotext_x[32] = 0x202d2030;
    genotext_x[34] = 0x20312030;
    genotext_x[36] = 0x202d2031;
    genotext_x[38] = 0x20302031;
    InitPhaseXNohhLookup4b(genotext_x);

    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = writebuf;
    {
      uint32_t clvl = 0;
      if (!(exportf_flags & kfExportfBgz)) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".haps");
      } else {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".haps.gz");
        clvl = kBgzfDefaultClvl;
      }
      reterr = InitBgzfCompressStreamEx(outname, 0, clvl, max_thread_ct, &bgzf);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto ExportOxHapslegend_ret_1;
      }
    }
    logprintfww5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, variant_uidx_start, &variant_uidx_base, &cur_bits);
    uint32_t chr_blen = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
        *chr_name_end++ = ' ';
        chr_blen = chr_name_end - chr_buf;
        is_x = (chr_idx == x_code);
        is_haploid = IsSet(cip->haploid_mask, chr_idx);
      }
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        if (unlikely((!refalt1_select) && (allele_idx_offsets[variant_uidx + 1] != allele_idx_offset_base + 2))) {
          logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
          goto ExportOxHapslegend_ret_INCONSISTENT_INPUT;
        }
      }
      if (refalt1_select) {
        allele_idx0 = refalt1_select[variant_uidx][ref_allele_last];
        allele_idx1 = refalt1_select[variant_uidx][1 - ref_allele_last];
      }
      uint32_t phasepresent_ct;
      reterr = PgrGet2P(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, allele_idx0, allele_idx1, simple_pgrp, genovec, phasepresent, phaseinfo, &phasepresent_ct);
      if (unlikely(reterr)) {
        goto ExportOxHapslegend_ret_PGR_FAIL;
      }
      ZeroTrailingNyps(sample_ct, genovec);
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      GenoarrCountFreqsUnsafe(genovec, sample_ct, genocounts);
      if (unlikely(phasepresent_ct != genocounts[1])) {
        logputs("\n");
        logerrprintf("Error: '--export haps%s' must be used with a fully phased dataset.\n", just_haps? "" : "legend");
        goto ExportOxHapslegend_ret_INCONSISTENT_INPUT;
      } else if (unlikely(genocounts[3])) {
        logputs("\n");
        logerrprintf("Error: '--export haps%s' cannot be used with missing genotype calls.\n", just_haps? "" : "legend");
        goto ExportOxHapslegend_ret_INCONSISTENT_INPUT;
      }
      if (is_haploid) {
        // verify that there are no het haploids/mixed MTs
        if (is_x) {
          GenoarrCountSubsetFreqs(genovec, sex_male_collapsed_interleaved, sample_ct, male_ct, genocounts);
        }
        if (unlikely(genocounts[1])) {
          logputs("\n");
          logerrprintfww("Error: '--export haps%s' cannot be used when heterozygous haploid or mixed MT calls are present.%s\n", just_haps? "" : "legend", (is_x && (variant_bps[variant_uidx] <= 2781479))? " (Did you forget --split-par?)" : "");
          goto ExportOxHapslegend_ret_INCONSISTENT_INPUT;
        }
      }
      if (just_haps) {
        write_iter = memcpya(write_iter, chr_buf, chr_blen);
        write_iter = strcpyax(write_iter, variant_ids[variant_uidx], ' ');
        write_iter = u32toa_x(variant_bps[variant_uidx], ' ', write_iter);
        const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
        write_iter = strcpyax(write_iter, cur_alleles[allele_idx0], ' ');
        write_iter = strcpyax(write_iter, cur_alleles[allele_idx1], ' ');
      }
      if (!is_x) {
        if (!phasepresent_ct) {
          GenoarrLookup16x4bx2(genovec, is_haploid? genotext_haploid : genotext_diploid, sample_ct, write_iter);
        } else {
          BitvecAnd(phasepresent, sample_ctl, phaseinfo);
          PhaseLookup4b(genovec, phasepresent, phaseinfo, genotext_diploid, sample_ct, write_iter);
        }
      } else {
        if (!phasepresent_ct) {
          GenoarrSexLookup4b(genovec, sex_male_collapsed, genotext_x, sample_ct, write_iter);
        } else {
          BitvecAnd(phasepresent, sample_ctl, phaseinfo);
          PhaseXNohhLookup4b(genovec, phasepresent, phaseinfo, sex_male_collapsed, genotext_x, sample_ct, write_iter);
        }
      }
      write_iter = &(write_iter[sample_ct * 4]);
      DecrAppendBinaryEoln(&write_iter);
      if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
        goto ExportOxHapslegend_ret_WRITE_FAIL;
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
    if (unlikely(bgzfclose_flush(writebuf_flush, write_iter, &bgzf, &reterr))) {
      goto ExportOxHapslegend_ret_1;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
  }
  while (0) {
  ExportOxHapslegend_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportOxHapslegend_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportOxHapslegend_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportOxHapslegend_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ExportOxHapslegend_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  }
 ExportOxHapslegend_ret_1:
  fclose_cond(outfile);
  CleanupBgzfCompressStream(&bgzf, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

static const uint16_t kBgen11HardcallUsis[16] = {
  32768, 0, 0, 0,
  0, 32768, 0, 0,
  0, 0, 32768, 0,
  0, 0, 0, 0};

typedef struct ExportBgen11CtxStruct {
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select);
  uintptr_t* sex_male_collapsed;
  uint32_t sample_ct;
  uint32_t y_start;
  uint32_t y_end;
  uint32_t ref_allele_last;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** dosage_presents;
  Dosage** dosage_mains;
  uint32_t* read_variant_uidx_starts;

  uint32_t cur_block_write_ct;

  struct libdeflate_compressor** libdeflate_compressors;
  uintptr_t** missing_acc1;
  uint16_t** bgen_geno_bufs;
  uint32_t bgen_compressed_buf_max;

  unsigned char* writebufs[2];
  uint32_t* variant_bytects[2];

  // high 32 bits = variant_uidx, earlier one takes precedence
  // low 32 bits = uint32_t(PglErr)
  uint64_t err_info;
} ExportBgen11Ctx;

THREAD_FUNC_DECL ExportBgen11Thread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  ExportBgen11Ctx* ctx = S_CAST(ExportBgen11Ctx*, arg->sharedp->context);

  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  uintptr_t* genovec = ctx->genovecs[tidx];
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t acc1_vec_ct = BitCtToVecCt(sample_ct);
  const uint32_t acc4_vec_ct = acc1_vec_ct * 4;
  const uint32_t acc8_vec_ct = acc1_vec_ct * 8;
  uintptr_t* missing_acc1 = ctx->missing_acc1[tidx];
  VecW* missing_acc4 = &(R_CAST(VecW*, missing_acc1)[acc1_vec_ct]);
  VecW* missing_acc8 = &(missing_acc4[acc4_vec_ct]);
  VecW* missing_acc32 = &(missing_acc8[acc8_vec_ct]);
  uintptr_t* dosage_present = ctx->dosage_presents? ctx->dosage_presents[tidx] : nullptr;
  Dosage* dosage_main = dosage_present? ctx->dosage_mains[tidx] : nullptr;
  uint16_t* bgen_geno_buf = ctx->bgen_geno_bufs[tidx];
  struct libdeflate_compressor* compressor = ctx->libdeflate_compressors[tidx];
  const uintptr_t* variant_include = ctx->variant_include;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uintptr_t* sample_include = ctx->sample_include;
  const uint32_t* sample_include_cumulative_popcounts = ctx->sample_include_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed = ctx->sex_male_collapsed;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uint32_t sample_ctl2_m1 = NypCtToWordCt(sample_ct) - 1;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t bgen_geno_buf_blen = 6 * sample_ct;
  const uint32_t bgen_compressed_buf_max = ctx->bgen_compressed_buf_max;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select) = ctx->refalt1_select;
  uint32_t is_y = 0;
  uint32_t y_thresh = ctx->y_start;
  const uint32_t y_end = ctx->y_end;
  const uint32_t ref_allele_last = ctx->ref_allele_last;
  uint32_t vidx_rem15 = 15;
  uint32_t vidx_rem255d15 = 17;
  uint32_t ref_allele_idx = 0;
  uint32_t parity = 0;
  ZeroWArr(acc1_vec_ct * kWordsPerVec * 45, missing_acc1);
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
    uint32_t write_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t write_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* writebuf_iter = &(ctx->writebufs[parity][write_idx * S_CAST(uintptr_t, bgen_compressed_buf_max)]);
    uint32_t* variant_bytect_iter = &(ctx->variant_bytects[parity][write_idx]);
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
    for (; write_idx != write_idx_end; ++write_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= y_thresh) {
        if (variant_uidx < y_end) {
          y_thresh = y_end;
          is_y = 1;
        } else {
          y_thresh = UINT32_MAX;
          is_y = 0;
        }
      }
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx][0];
      }
      if (allele_idx_offsets) {
        if (allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx] != 2) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, kPglRetInconsistentInput);
          goto ExportBgen11Thread_err;
        }
      }
      uint32_t dosage_ct;
      PglErr reterr = PgrGetD(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
      if (unlikely(reterr)) {
        new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
        goto ExportBgen11Thread_err;
      }
      if (ref_allele_idx + ref_allele_last == 1) {
        GenovecInvertUnsafe(sample_ct, genovec);
        BiallelicDosage16Invert(dosage_ct, dosage_main);
      }
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint16_t* bgen_geno_buf_iter = bgen_geno_buf;
      if (!dosage_ct) {
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t geno_word = genovec[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            // possible todo: change this to 16x6bx2 lookup table
            memcpy_k(bgen_geno_buf_iter, &(kBgen11HardcallUsis[(geno_word & 3) * 4]), 6);
            bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
            geno_word >>= 2;
          }
        }
      } else {
        const Halfword* dosage_present_alias = R_CAST(Halfword*, dosage_present);
        const Dosage* dosage_main_iter = dosage_main;
        for (uint32_t widx = 0; ; ++widx) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t geno_word = genovec[widx];
          uint32_t dosage_present_hw = dosage_present_alias[widx];
          if (!dosage_present_hw) {
            for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
              memcpy_k(bgen_geno_buf_iter, &(kBgen11HardcallUsis[(geno_word & 3) * 4]), 6);
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
              geno_word >>= 2;
            }
          } else {
            for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
              if (dosage_present_hw & 1) {
                uint32_t dosage_int = *dosage_main_iter++;
                dosage_int *= 2;
                if (dosage_int <= kDosageMax) {
                  *bgen_geno_buf_iter++ = kDosageMax - dosage_int;
                  *bgen_geno_buf_iter++ = dosage_int;
                  *bgen_geno_buf_iter++ = 0;
                } else {
                  dosage_int -= kDosageMax;
                  *bgen_geno_buf_iter++ = 0;
                  *bgen_geno_buf_iter++ = kDosageMax - dosage_int;
                  *bgen_geno_buf_iter++ = dosage_int;
                }
              } else {
                memcpy_k(bgen_geno_buf_iter, &(kBgen11HardcallUsis[(geno_word & 3) * 4]), 6);
                bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
              }
              geno_word >>= 2;
              dosage_present_hw >>= 1;
            }
          }
        }
      }
      uintptr_t compressed_blen = libdeflate_zlib_compress(compressor, bgen_geno_buf, bgen_geno_buf_blen, writebuf_iter, bgen_compressed_buf_max);
      assert(compressed_blen);
      *variant_bytect_iter++ = compressed_blen;
      writebuf_iter = &(writebuf_iter[bgen_compressed_buf_max]);
      // bugfix (13 Apr 2018): this missingness calculation was only taking
      // hardcalls into account, which is inappropriate for .gen/.bgen
      GenoarrToMissingnessUnsafe(genovec, sample_ct, missing_acc1);
      if (dosage_ct) {
        BitvecInvmask(dosage_present, sample_ctl, missing_acc1);
      }
      if (is_y) {
        BitvecAnd(sex_male_collapsed, sample_ctl, missing_acc1);
      }
      VcountIncr1To4(missing_acc1, acc1_vec_ct, missing_acc4);
      if (!(--vidx_rem15)) {
        Vcount0Incr4To8(acc4_vec_ct, missing_acc4, missing_acc8);
        vidx_rem15 = 15;
        if (!(--vidx_rem255d15)) {
          Vcount0Incr8To32(acc8_vec_ct, missing_acc8, missing_acc32);
          vidx_rem255d15 = 17;
        }
      }
    }
    while (0) {
    ExportBgen11Thread_err:
      UpdateU64IfSmaller(new_err_info, &ctx->err_info);
      break;
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  VcountIncr4To8(missing_acc4, acc4_vec_ct, missing_acc8);
  VcountIncr8To32(missing_acc8, acc8_vec_ct, missing_acc32);
  THREAD_RETURN;
}

PglErr ExportBgen11(const char* outname, const uintptr_t* sample_include, uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, ExportfFlags exportf_flags, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uint32_t* sample_missing_geno_cts) {
  // isomorphic to ExportOxGen().
  assert(sample_ct);
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  ExportBgen11Ctx ctx;
  {
    const uint32_t max_chr_slen = GetMaxChrSlen(cip);
    if (unlikely(BIGSTACK_ALLOC_X(struct libdeflate_compressor*, max_thread_ct, &ctx.libdeflate_compressors))) {
      goto ExportBgen11_ret_NOMEM;
    }
    ZeroPtrArr(max_thread_ct, ctx.libdeflate_compressors);
    // allocate the first compressor so we can call
    // libdeflate_zlib_compress_bound().
    // could add a --bgz-level flag analogous to --zst-level, of course.
    ctx.libdeflate_compressors[0] = libdeflate_alloc_compressor(6);
    if (unlikely(!ctx.libdeflate_compressors[0])) {
      goto ExportBgen11_ret_NOMEM;
    }
    const uintptr_t bgen_compressed_buf_max = libdeflate_zlib_compress_bound(ctx.libdeflate_compressors[0], 6LU * sample_ct);
#ifdef __LP64__
    if (unlikely(bgen_compressed_buf_max > UINT32_MAX)) {
      logerrputs("Error: Too many samples for .bgen format.\n");
      goto ExportBgen11_ret_INCONSISTENT_INPUT;
    }
#endif
    ctx.bgen_compressed_buf_max = bgen_compressed_buf_max;
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    const uintptr_t writebuf_len = bgen_compressed_buf_max + 2 * max_allele_slen + 2 * kMaxIdSlen + 32;
    char* chr_buf;
    unsigned char* writebuf;
    if (unlikely(
            bigstack_alloc_c(max_chr_slen, &chr_buf) ||
            bigstack_alloc_uc(writebuf_len, &writebuf) ||
            bigstack_alloc_w(sample_ctl, &ctx.sex_male_collapsed))) {
      goto ExportBgen11_ret_NOMEM;
    }
    CopyBitarrSubset(sex_male, sample_include, sample_ct, ctx.sex_male_collapsed);

    const uintptr_t max_write_block_byte_ct = bigstack_left() / 4;
    uint32_t max_write_block_size = kPglVblockSize;
    for (; ; max_write_block_size /= 2) {
      // limit each write buffer to 1/4 of remaining workspace
      if (S_CAST(uint64_t, bgen_compressed_buf_max + sizeof(int32_t)) * max_write_block_size <= max_write_block_byte_ct) {
        break;
      }
      if (unlikely(max_write_block_size <= kBitsPerVec)) {
        goto ExportBgen11_ret_NOMEM;
      }
    }
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // seems to saturate around this point
    // todo: retest with libdeflate
    if (calc_thread_ct > 15) {
      calc_thread_ct = 15;
    }
    if (unlikely(
            bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(ctx.writebufs[0])) ||
            bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(ctx.writebufs[1])) ||
            bigstack_alloc_u32(max_write_block_size, &(ctx.variant_bytects[0])) ||
            bigstack_alloc_u32(max_write_block_size, &(ctx.variant_bytects[1])) ||
            bigstack_alloc_wp(calc_thread_ct, &ctx.missing_acc1) ||
            bigstack_alloc_u16p(calc_thread_ct, &ctx.bgen_geno_bufs))) {
      goto ExportBgen11_ret_NOMEM;
    }
    // we allocated [0] earlier
    for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
      ctx.libdeflate_compressors[tidx] = libdeflate_alloc_compressor(6);
      if (unlikely(!ctx.libdeflate_compressors[tidx])) {
        goto ExportBgen11_ret_NOMEM;
      }
    }

    const uint32_t acc1_vec_ct = BitCtToVecCt(sample_ct);
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    const uintptr_t track_missing_cacheline_ct = VecCtToCachelineCt(acc1_vec_ct * 45);
    // no overflow risk here, thanks to compress_bound() check above
    const uintptr_t bgen_geno_cacheline_ct = DivUp(6 * sample_ct, kCacheline);
    const uintptr_t thread_xalloc_cacheline_ct = track_missing_cacheline_ct + bgen_geno_cacheline_ct;
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    // defensive
    ctx.dosage_presents = nullptr;
    ctx.dosage_mains = nullptr;
    uint32_t read_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &ctx.genovecs, nullptr, nullptr, nullptr, dosage_is_present? (&ctx.dosage_presents) : nullptr, dosage_is_present? (&ctx.dosage_mains) : nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto ExportBgen11_ret_NOMEM;
    }
    if (read_block_size > max_write_block_size) {
      read_block_size = max_write_block_size;
    }

    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportBgen11_ret_OPEN_FAIL;
    }
    // bgen 1.1 header
    // note that \xxx character constants are interpreted in octal, so \24 is
    // decimal 20, etc.
    memcpy(writebuf, "\24\0\0\0\24\0\0", 8);
    memcpy(&(writebuf[8]), &variant_ct, 4);
    memcpy(&(writebuf[12]), &sample_ct, 4);
    memcpy(&(writebuf[16]), "bgen\5\0\0", 8);
    if (unlikely(fwrite_checked(writebuf, 24, outfile))) {
      goto ExportBgen11_ret_WRITE_FAIL;
    }

    const uint32_t ref_allele_last = !(exportf_flags & kfExportfRefFirst);
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      ctx.missing_acc1[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(track_missing_cacheline_ct * kCacheline));
      ctx.bgen_geno_bufs[tidx] = S_CAST(uint16_t*, bigstack_alloc_raw(bgen_geno_cacheline_ct * kCacheline));
    }
    ctx.sample_ct = sample_ct;
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.sample_include = sample_include;
    ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto ExportBgen11_ret_NOMEM;
    }
    ctx.refalt1_select = refalt1_select;
    GetXymtStartAndEnd(cip, kChrOffsetY, &ctx.y_start, &ctx.y_end);
    ctx.ref_allele_last = ref_allele_last;
    ctx.err_info = (~0LLU) << 32;
    SetThreadFuncAndData(ExportBgen11Thread, &ctx, &tg);

    // 6 bytes present at start of every bgen-1.1 variant record
    memcpy(writebuf, &sample_ct, 4);
    memcpy(&(writebuf[4]), "\0", 2);

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
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_slen = 0;

    uint32_t prev_block_write_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    logprintfww5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto ExportBgen11_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, ctx.err_info);
        if (unlikely(reterr)) {
          if (reterr == kPglRetInconsistentInput) {
            logputs("\n");
            logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
          }
          goto ExportBgen11_ret_PGR_FAIL;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_write_ct = cur_block_write_ct;
        ComputeUidxStartPartition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_write_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto ExportBgen11_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const unsigned char* compressed_data_iter = ctx.writebufs[parity];
        const uint32_t* variant_bytect_iter = ctx.variant_bytects[parity];
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_write_ct; ++variant_bidx) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            const char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
            chr_slen = chr_name_end - chr_buf;
          }
          const char* cur_variant_id = variant_ids[write_variant_uidx];
          const uint32_t id_slen = strlen(cur_variant_id);
          unsigned char* writebuf_iter = memcpyua(&(writebuf[6]), &id_slen, 2);
          writebuf_iter = memcpyua(writebuf_iter, cur_variant_id, id_slen);
          AppendU16(chr_slen, &writebuf_iter);
          writebuf_iter = memcpyua(writebuf_iter, chr_buf, chr_slen);
          AppendU32(variant_bps[write_variant_uidx], &writebuf_iter);
          uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            if (unlikely((!refalt1_select) && (allele_idx_offsets[write_variant_uidx + 1] != allele_idx_offset_base + 2))) {
              logputs("\n");
              logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
              goto ExportBgen11_ret_INCONSISTENT_INPUT;
            }
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          if (refalt1_select) {
            ref_allele_idx = refalt1_select[write_variant_uidx][0];
            alt1_allele_idx = refalt1_select[write_variant_uidx][1];
          }
          const char* first_allele;
          const char* second_allele;
          if (ref_allele_last) {
            first_allele = cur_alleles[alt1_allele_idx];
            second_allele = cur_alleles[ref_allele_idx];
          } else {
            first_allele = cur_alleles[ref_allele_idx];
            second_allele = cur_alleles[alt1_allele_idx];
          }
          uint32_t allele_slen = strlen(first_allele);
          AppendU32(allele_slen, &writebuf_iter);
          writebuf_iter = memcpyua(writebuf_iter, first_allele, allele_slen);
          allele_slen = strlen(second_allele);
          AppendU32(allele_slen, &writebuf_iter);
          writebuf_iter = memcpyua(writebuf_iter, second_allele, allele_slen);
          const uint32_t cur_variant_bytect = *variant_bytect_iter++;
          AppendU32(cur_variant_bytect, &writebuf_iter);
          writebuf_iter = memcpyua(writebuf_iter, compressed_data_iter, cur_variant_bytect);
          compressed_data_iter = &(compressed_data_iter[bgen_compressed_buf_max]);
          if (unlikely(fwrite_checked(writebuf, writebuf_iter - writebuf, outfile))) {
            goto ExportBgen11_ret_WRITE_FAIL;
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
      prev_block_write_ct = cur_block_write_ct;
      variant_idx += cur_block_write_ct;
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(fclose_null(&outfile))) {
      goto ExportBgen11_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    const uint32_t sample_ctav = acc1_vec_ct * kBitsPerVec;
    const uintptr_t acc32_offset = acc1_vec_ct * (13 * k1LU * kWordsPerVec);
    uint32_t* scrambled_missing_cts = R_CAST(uint32_t*, &(ctx.missing_acc1[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
      const uint32_t* thread_scrambled_missing_cts = R_CAST(uint32_t*, &(ctx.missing_acc1[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii != sample_ctav; ++uii) {
        scrambled_missing_cts[uii] += thread_scrambled_missing_cts[uii];
      }
    }
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = VcountScramble1(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  ExportBgen11_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportBgen11_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportBgen11_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ExportBgen11_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportBgen11_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ExportBgen11_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  if (ctx.libdeflate_compressors) {
    for (uint32_t tidx = 0; tidx != max_thread_ct; ++tidx) {
      if (!ctx.libdeflate_compressors[tidx]) {
        break;
      }
      libdeflate_free_compressor(ctx.libdeflate_compressors[tidx]);
    }
  }
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

// For constant-ploidy cases.  Writes up to 7 bytes past the end.
// For what it's worth, bgen_geno_buf_iter is always word-aligned on entry (8
// bytes after beginning of buffer, which is at the start of a cacheline), and
// it's easy to force vector-alignment if it would be useful.
unsigned char* FillBgen13PloidyAndMissingness(const uintptr_t* genovec, const uintptr_t* dosage_present, uintptr_t ploidy, uint32_t sample_ct, uint32_t dosage_ct, unsigned char* bgen_geno_buf_iter) {
  uint64_t* sample_ploidy_and_missingness_alias = R_CAST(uint64_t*, bgen_geno_buf_iter);
  // bugfix (30 Aug 2018): need to spell out mask in 32-bit case
  const uint64_t ploidy_u64 = ploidy * 0x101010101010101LLU;
  const uint32_t sample_ct8 = DivUp(sample_ct, 8);
  const uint16_t* genovec_alias = R_CAST(const uint16_t*, genovec);
  if (!dosage_ct) {
    for (uint32_t write_widx = 0; write_widx != sample_ct8; ++write_widx) {
      const uint64_t cur_geno8 = genovec_alias[write_widx];
      uint64_t cur_missing8 = cur_geno8 & (cur_geno8 >> 1);
#ifdef USE_AVX2
      // 0,2,4,6...14 -> 7,15,23,...,63
      // todo: try inverse-movespreadmask
      cur_missing8 = _pdep_u64(cur_missing8 & 0x5555, 0x80 * kMask1111);
#else
      // 0,2,4,6,8,10,12,14 -> 0,2,4,6,32,34,36,38
      cur_missing8 = (cur_missing8 | (cur_missing8 << 24)) & 0x5500000055LLU;
      // Multiply by number with bits 7,13,19,25 set, then mask out all but the
      // bits we want.  There are potential carries at bits 13, 19, 25, 45, 51,
      // and 57, but none of them hurt us.
      cur_missing8 = (cur_missing8 * 0x2082080) & (0x8080808080808080LLU);
#endif
      sample_ploidy_and_missingness_alias[write_widx] = cur_missing8 + ploidy_u64;
    }
  } else {
    // don't bother with word-based 32-bit code because this case becomes
    // more annoying
    const unsigned char* dosage_present_alias = R_CAST(const unsigned char*, dosage_present);
    for (uint32_t write_widx = 0; write_widx != sample_ct8; ++write_widx) {
      const uint64_t cur_geno8 = genovec_alias[write_widx];
      uint64_t cur_dosage_missing8 = ~dosage_present_alias[write_widx];
      uint64_t cur_hardcall_missing8 = cur_geno8 & (cur_geno8 >> 1);
#ifdef USE_AVX2
      // 0,2,4,6...14 -> 7,15,23,...,63
      cur_hardcall_missing8 = _pdep_u64(cur_hardcall_missing8 & 0x5555, 0x80 * kMask1111);
      cur_dosage_missing8 = _pdep_u64(cur_dosage_missing8, 0x80 * kMask0101);
      const uint64_t cur_missing8 = cur_hardcall_missing8 & cur_dosage_missing8;
#else
      // 0,2,4,6,8,10,12,14 -> 0,2,4,6,32,34,36,38
      cur_hardcall_missing8 = (cur_hardcall_missing8 | (cur_hardcall_missing8 << 24)) & 0x5500000055LLU;
      // -> 0,8,16,24,32,40,48,56, with extraneous bits set
      cur_hardcall_missing8 = cur_hardcall_missing8 * 0x41041;
      // 0,1,2,3,4,5,6,7 -> 0,8,16,24,32,40,48,56, with extraneous bits set
      // this operation also appears in GflagsVfilter(), may want to wrap it
      cur_dosage_missing8 = ((cur_dosage_missing8 & 0xfe) * 0x2040810204080LLU) | (cur_dosage_missing8 & 1);
      const uint64_t cur_missing8 = ((cur_hardcall_missing8 & cur_dosage_missing8) & 0x101010101010101LLU) << 7;
#endif
      sample_ploidy_and_missingness_alias[write_widx] = cur_missing8 + ploidy_u64;
    }
  }
  return &(bgen_geno_buf_iter[sample_ct]);
}

uint32_t NoFemaleMissing(const uintptr_t* genovec, const uintptr_t* dosage_present, const uintptr_t* sex_female, uint32_t sample_ctl2, uint32_t dosage_ct) {
  if (dosage_ct) {
    const uint32_t sample_ctl = DivUp(sample_ctl2, 2);
    if (!IntersectionIsEmpty(sex_female, dosage_present, sample_ctl)) {
      return 0;
    }
  }
  const Halfword* sex_female_alias = R_CAST(const Halfword*, sex_female);
  // todo: try vectorizing this loop
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uintptr_t geno_word = genovec[widx];
    const uintptr_t cur_female_01 = UnpackHalfwordToWord(sex_female_alias[widx]);
    if (geno_word & (geno_word >> 1) & cur_female_01) {
      return 0;
    }
  }
  return 1;
}

// These may belong in plink2_base or plink2_common.
// cur_write_bit_idx stays in [0, kBitsPerWord - 1]
static inline void AppendBits(uint32_t bit_ct, uintptr_t payload, uintptr_t* cur_write_bits_ptr, uint32_t* cur_write_bit_idx_ptr, uintptr_t** probs_write_iter_ptr) {
  uint32_t cur_write_bit_idx = *cur_write_bit_idx_ptr;
  *cur_write_bits_ptr |= payload << cur_write_bit_idx;
  cur_write_bit_idx += bit_ct;
  if (cur_write_bit_idx >= kBitsPerWord) {
    **probs_write_iter_ptr = *cur_write_bits_ptr;
    *probs_write_iter_ptr += 1;
    cur_write_bit_idx -= kBitsPerWord;
    *cur_write_bits_ptr = payload >> (bit_ct - cur_write_bit_idx);
  }
  *cur_write_bit_idx_ptr = cur_write_bit_idx;
}

static inline void Append0Bits(uint32_t bit_ct, uintptr_t* cur_write_bits_ptr, uint32_t* cur_write_bit_idx_ptr, uintptr_t** probs_write_iter_ptr) {
  uint32_t cur_write_bit_idx = *cur_write_bit_idx_ptr;
  cur_write_bit_idx += bit_ct;
  if (cur_write_bit_idx >= kBitsPerWord) {
    **probs_write_iter_ptr = *cur_write_bits_ptr;
    *probs_write_iter_ptr += 1;
    cur_write_bit_idx -= kBitsPerWord;
    *cur_write_bits_ptr = 0;
  }
  *cur_write_bit_idx_ptr = cur_write_bit_idx;
}

typedef struct Bgen13TablesStruct {
  // Precomputed 1-genotype-at-a-time tables (4 entries)
  uint32_t* diploid_basic_table;
  uint16_t* haploid_basic_table;

  // Precomputed 4-genotype-at-a-time tables (256 entries; multiply index by
  // two and copy two values at a time for table16).
  uint64_t* diploid_hardcall_table8;
  uint64_t* diploid_hardcall_table16;

  // Precomputed 1-genotype-at-a-time table (16 entries: phasepresent = bit 2,
  // phaseinfo = bit 3)
  uint32_t* diploid_phased_hardcall_table;

  // 4-genotype-at-a-time tables.
  uint32_t* haploid_hardcall_table8;
  uint64_t* haploid_hardcall_table16;

  uint32_t bit_precision;
} Bgen13Tables;

BoolErr ConstructBgen13LookupTables(uint32_t exportf_bits, Bgen13Tables* tablesp) {
  const uint32_t max_val = (1U << exportf_bits) - 1;
  const uint32_t half_val_roundeven = ((max_val + 1) / 2) & (~1);
  if (bigstack_alloc_u32(4, &tablesp->diploid_basic_table) ||
      bigstack_alloc_u32(16, &tablesp->diploid_phased_hardcall_table) ||
      bigstack_alloc_u16(4, &tablesp->haploid_basic_table)) {
    return 1;
  }
  // reference allele is exported as second allele so default import and
  // export settings work together
  // (yeah, this will make multiallelic case annoying)
  // unphased diploid:
  //   [0] = 0, 0
  //   [1] = 0, max
  //   [2] = max, 0
  //   [3] = 0, 0
  tablesp->diploid_basic_table[0] = 0;
  tablesp->diploid_basic_table[1] = max_val << exportf_bits;
  tablesp->diploid_basic_table[2] = max_val;
  tablesp->diploid_basic_table[3] = 0;

  // phased diploid:
  //   [0] = 0, 0
  //   [1] unphased = half, half
  //   [1] phased, phaseinfo clear = 0, max
  //   [1] phased, phaseinfo set = max, 0
  //   [2] = max, max
  //   [3] = 0, 0
  tablesp->diploid_phased_hardcall_table[0] = 0;
  tablesp->diploid_phased_hardcall_table[1] = half_val_roundeven * ((1U << exportf_bits) + 1);
  tablesp->diploid_phased_hardcall_table[2] = max_val * ((1U << exportf_bits) + 1);
  tablesp->diploid_phased_hardcall_table[3] = 0;
  // might need to defend against garbage phaseinfo
  // cheap to defend against garbage phasepresent
  memcpy(&(tablesp->diploid_phased_hardcall_table[4]), tablesp->diploid_phased_hardcall_table, 4 * sizeof(int32_t));
  memcpy(&(tablesp->diploid_phased_hardcall_table[8]), tablesp->diploid_phased_hardcall_table, 8 * sizeof(int32_t));
  tablesp->diploid_phased_hardcall_table[5] = max_val << exportf_bits;
  tablesp->diploid_phased_hardcall_table[13] = max_val;

  // haploid:
  //   [0] = 0
  //   [1] = half, may as well round-to-even
  //   [2] = max
  //   [3] = 0
  tablesp->haploid_basic_table[0] = 0;
  tablesp->haploid_basic_table[1] = half_val_roundeven;
  tablesp->haploid_basic_table[2] = max_val;
  tablesp->haploid_basic_table[3] = 0;

  if (exportf_bits <= 8) {
    // can conditionally skip some of these tables, but whatever
    if (unlikely(
            bigstack_alloc_u64(256, &tablesp->diploid_hardcall_table8) ||
            bigstack_alloc_u32(256, &tablesp->haploid_hardcall_table8))) {
      return 1;
    }
    tablesp->diploid_hardcall_table16 = nullptr;
    tablesp->haploid_hardcall_table16 = nullptr;
    uint64_t* write64_iter = tablesp->diploid_hardcall_table8;
    for (uint32_t uii = 0; uii != 4; ++uii) {
      const uint64_t entry3 = S_CAST(uint64_t, tablesp->diploid_basic_table[uii]) << (exportf_bits * 6);
      for (uint32_t ujj = 0; ujj != 4; ++ujj) {
        const uint64_t entry23 = entry3 | (S_CAST(uint64_t, tablesp->diploid_basic_table[ujj]) << (exportf_bits * 4));
        for (uint32_t ukk = 0; ukk != 4; ++ukk) {
          const uint64_t entry123 = entry23 | (S_CAST(uint64_t, tablesp->diploid_basic_table[ukk]) << (exportf_bits * 2));
          for (uint32_t umm = 0; umm != 4; ++umm) {
            *write64_iter++ = entry123 | tablesp->diploid_basic_table[umm];
          }
        }
      }
    }

    uint32_t* write32_iter = tablesp->haploid_hardcall_table8;
    for (uint32_t uii = 0; uii != 4; ++uii) {
      const uint32_t entry3 = S_CAST(uint32_t, tablesp->haploid_basic_table[uii]) << (exportf_bits * 3);
      for (uint32_t ujj = 0; ujj != 4; ++ujj) {
        const uint32_t entry23 = entry3 | (S_CAST(uint32_t, tablesp->haploid_basic_table[ujj]) << (exportf_bits * 2));
        for (uint32_t ukk = 0; ukk != 4; ++ukk) {
          const uint32_t entry123 = entry23 | (S_CAST(uint32_t, tablesp->haploid_basic_table[ukk]) << exportf_bits);
          for (uint32_t umm = 0; umm != 4; ++umm) {
            *write32_iter++ = entry123 | tablesp->haploid_basic_table[umm];
          }
        }
      }
    }
  } else {
    if (unlikely(
            bigstack_alloc_u64(512, &tablesp->diploid_hardcall_table16) ||
            bigstack_alloc_u64(256, &tablesp->haploid_hardcall_table16))) {
      return 1;
    }
    tablesp->diploid_hardcall_table8 = nullptr;
    tablesp->haploid_hardcall_table8 = nullptr;
    uint64_t* write64_iter = tablesp->diploid_hardcall_table16;
    uint64_t entry23_low = 0;
    for (uint32_t uii = 0; uii != 4; ++uii) {
      const uint64_t entry3 = S_CAST(uint64_t, tablesp->diploid_basic_table[uii]) << (exportf_bits * 2);
      for (uint32_t ujj = 0; ujj != 4; ++ujj) {
        uint64_t entry23_high = entry3 | tablesp->diploid_basic_table[ujj];
        if (exportf_bits < 16) {
          entry23_low = entry23_high << (4 * exportf_bits);
          entry23_high = entry23_high >> (64 - 4 * exportf_bits);
        }
        for (uint32_t ukk = 0; ukk != 4; ++ukk) {
          const uint64_t entry123_low = entry23_low | (S_CAST(uint64_t, tablesp->diploid_basic_table[ukk]) << (exportf_bits * 2));
          for (uint32_t umm = 0; umm != 4; ++umm) {
            *write64_iter++ = entry123_low | tablesp->diploid_basic_table[umm];
            *write64_iter++ = entry23_high;
          }
        }
      }
    }

    write64_iter = tablesp->haploid_hardcall_table16;
    for (uint32_t uii = 0; uii != 4; ++uii) {
      const uint64_t entry3 = S_CAST(uint64_t, tablesp->haploid_basic_table[uii]) << (exportf_bits * 3);
      for (uint32_t ujj = 0; ujj != 4; ++ujj) {
        const uint64_t entry23 = entry3 | (S_CAST(uint64_t, tablesp->haploid_basic_table[ujj]) << (exportf_bits * 2));
        for (uint32_t ukk = 0; ukk != 4; ++ukk) {
          const uint64_t entry123 = entry23 | (S_CAST(uint64_t, tablesp->haploid_basic_table[ukk]) << exportf_bits);
          for (uint32_t umm = 0; umm != 4; ++umm) {
            *write64_iter++ = entry123 | tablesp->haploid_basic_table[umm];
          }
        }
      }
    }
  }
  tablesp->bit_precision = exportf_bits;
  return 0;
}

typedef struct ExportBgen13CtxStruct {
  const uintptr_t* variant_include;
  const ChrInfo* cip;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select);
  uintptr_t* sex_male_collapsed;
  uintptr_t* sex_female_collapsed;
  Bgen13Tables tables;
  uint32_t sample_ct;
  uint32_t ref_allele_last;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** phasepresents;
  uintptr_t** phaseinfos;
  uintptr_t** dosage_presents;
  Dosage** dosage_mains;
  uintptr_t** dphase_presents;
  SDosage** dphase_deltas;
  uint32_t* read_variant_uidx_starts;

  uint32_t cur_block_write_ct;

  struct libdeflate_compressor** libdeflate_compressors;
  uintptr_t** missing_acc1;
  unsigned char** uncompressed_bgen_geno_bufs;
  uint32_t bgen_compressed_buf_max;

  unsigned char* writebufs[2];
  uint32_t* variant_bytects[2];

  // high 32 bits = variant_uidx, earlier one takes precedence
  // low 32 bits = uint32_t(PglErr)
  uint64_t err_info;
} ExportBgen13Ctx;

THREAD_FUNC_DECL ExportBgen13Thread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  ExportBgen13Ctx* ctx = S_CAST(ExportBgen13Ctx*, arg->sharedp->context);

  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgenVariant pgv;
  pgv.genovec = ctx->genovecs[tidx];
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t acc1_vec_ct = BitCtToVecCt(sample_ct);
  const uint32_t acc4_vec_ct = acc1_vec_ct * 4;
  const uint32_t acc8_vec_ct = acc1_vec_ct * 8;
  uintptr_t* missing_acc1 = ctx->missing_acc1[tidx];
  VecW* missing_acc4 = &(R_CAST(VecW*, missing_acc1)[acc1_vec_ct]);
  VecW* missing_acc8 = &(missing_acc4[acc4_vec_ct]);
  VecW* missing_acc32 = &(missing_acc8[acc8_vec_ct]);
  // todo: multiallelic support
  pgv.phasepresent = nullptr;
  pgv.phaseinfo = nullptr;
  if (ctx->phasepresents) {
    pgv.phasepresent = ctx->phasepresents[tidx];
    pgv.phaseinfo = ctx->phaseinfos[tidx];
  }
  pgv.dosage_present = ctx->dosage_presents? ctx->dosage_presents[tidx] : nullptr;
  pgv.dosage_main = pgv.dosage_present? ctx->dosage_mains[tidx] : nullptr;
  pgv.dphase_present = ctx->dphase_presents? ctx->dphase_presents[tidx] : nullptr;
  pgv.dphase_delta = pgv.dphase_present? ctx->dphase_deltas[tidx] : nullptr;
  // Note that we may write up to 12 bytes past the end
  unsigned char* uncompressed_bgen_geno_buf = ctx->uncompressed_bgen_geno_bufs[tidx];
  struct libdeflate_compressor* compressor = ctx->libdeflate_compressors? ctx->libdeflate_compressors[tidx] : nullptr;
  const uint32_t zst_level = g_zst_level;
  const uintptr_t* variant_include = ctx->variant_include;
  const ChrInfo* cip = ctx->cip;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uintptr_t* sample_include = ctx->sample_include;
  const uint32_t* sample_include_cumulative_popcounts = ctx->sample_include_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed = ctx->sex_male_collapsed;
  const uintptr_t* sex_female_collapsed = ctx->sex_female_collapsed;
  const uint16_t* bgen_haploid_basic_table = ctx->tables.haploid_basic_table;
  const uint32_t* bgen_diploid_basic_table = ctx->tables.diploid_basic_table;
  const uint64_t* bgen_diploid_hardcall_table8 = ctx->tables.diploid_hardcall_table8;
  const uint64_t* bgen_diploid_hardcall_table16 = ctx->tables.diploid_hardcall_table16;
  const uint32_t* bgen_diploid_phased_hardcall_table = ctx->tables.diploid_phased_hardcall_table;
  const uint32_t* bgen_haploid_hardcall_table8 = ctx->tables.haploid_hardcall_table8;
  const uint64_t* bgen_haploid_hardcall_table16 = ctx->tables.haploid_hardcall_table16;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
  const uint32_t sample_ct4 = DivUp(sample_ct, 4);
  const uint32_t bit_precision = ctx->tables.bit_precision;
  // note that this assumes bit_precision <= 16
  const uint32_t two_byte_probs = (bit_precision > 8);
  const uint32_t max_output_val = (1U << bit_precision) - 1;
  const uint32_t bgen_compressed_buf_max = ctx->bgen_compressed_buf_max;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select) = ctx->refalt1_select;
  uint32_t chr_fo_idx = UINT32_MAX;  // deliberate overflow
  uint32_t chr_end = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t cur_y = 0;

  uint32_t is_haploid = 0;  // includes chrX
  // Unlike the VCF/BCF export case,
  // * we do not save heterozygous haploid calls as ploidy=2 unless there is
  //   also phase information.
  // * chrY female (but not unknown-sex) ploidy is 0 when genotype is missing.
  const uint32_t male_ct = PopcountWords(sex_male_collapsed, sample_ctl);
  const uint32_t female_ct = PopcountWords(sex_female_collapsed, sample_ctl);
  const uint32_t x_code = (male_ct != sample_ct)? cip->xymt_codes[kChrOffsetX] : UINT32_MAXM1;
  const uint32_t y_code = female_ct? cip->xymt_codes[kChrOffsetY] : UINT32_MAXM1;
  const uint32_t ref_allele_last = ctx->ref_allele_last;
  uint32_t vidx_rem15 = 15;
  uint32_t vidx_rem255d15 = 17;
  uint32_t ref_allele_idx = 0;
  uint32_t parity = 0;
  ZeroWArr(acc1_vec_ct * kWordsPerVec * 45, missing_acc1);
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
    uint32_t write_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t write_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* writebuf_iter = &(ctx->writebufs[parity][write_idx * S_CAST(uintptr_t, bgen_compressed_buf_max)]);
    // [2n] = 4 + compressed_len, [2n + 1] = uncompressed_len
    uint32_t* variant_bytect_iter = &(ctx->variant_bytects[parity][2 * write_idx]);
    uintptr_t variant_uidx_base;
    uintptr_t variant_include_bits;
    BitIter1Start(variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);
    for (; write_idx != write_idx_end; ++write_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        is_y = (chr_idx == y_code);
        is_haploid = IsSet(cip->haploid_mask, chr_idx);
        is_x = 0;
        if (chr_idx == x_code) {
          if (male_ct) {
            is_x = 1;
          } else {
            is_haploid = 0;
          }
        }
        cur_y = 0;
      }
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx][0];
      }
      if (allele_idx_offsets) {
        // todo: multiallelic cases
        if (allele_idx_offsets[variant_uidx + 1] - allele_idx_offsets[variant_uidx] != 2) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, kPglRetInconsistentInput);
          goto ExportBgen13Thread_err;
        }
      }
      // (note that multiallelic variants should always be stored using the
      // phased format; just set both sides equal when the call is actually
      // unphased.  this way there's no O([allele count]^2) bloat problem.)
      PglErr reterr = PgrGetDp(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, &pgv);
      if (unlikely(reterr)) {
        new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
        goto ExportBgen13Thread_err;
      }
      // note that this is inverted from bgen-1.1
      // todo: switch to allele-rotation function call
      if (ref_allele_idx + ref_allele_last != 1) {
        GenovecInvertUnsafe(sample_ct, pgv.genovec);
        if (pgv.phasepresent_ct) {
          BitvecInvert(sample_ctl, pgv.phaseinfo);
        }
        if (pgv.dosage_ct) {
          BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
          if (pgv.dphase_ct) {
            BiallelicDphase16Invert(pgv.dphase_ct, pgv.dphase_delta);
          }
        }
      }
      unsigned char* bgen_geno_buf_iter = uncompressed_bgen_geno_buf;
      // 4 bytes: # of samples
      // 2 bytes: # of alleles
      // 1 byte: minimum ploidy
      // 1 byte: maximum ploidy
      // sample_ct bytes: high bit = missing, low bits = ploidy
      // 1 byte: is_phased
      // 1 byte: bit_precision

      AppendU32(sample_ct, &bgen_geno_buf_iter);
      if (is_y) {
        // All-nonfemales case treated as haploid general case.
        // The difference here is that missing female calls are saved as
        // ploidy 0; this basically forces us to write a single value at
        // a time.
        cur_y = !NoFemaleMissing(pgv.genovec, pgv.dosage_present, sex_female_collapsed, sample_ctl2, pgv.dosage_ct);
      }
      if (pgv.dphase_ct && (bit_precision < 15)) {
        // Theoretically possible for all dphase_delta values to be too small
        // to ever make left dosage != right dosage.  If so, allow unphased
        // output.
        uintptr_t sample_widx = 0;
        uintptr_t dosage_present_bits = pgv.dosage_present[0];
        uint32_t dphase_idx = 0;
        uint32_t dosage_idx = 0;
        for (; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
          const uintptr_t lowbit = BitIter1y(pgv.dosage_present, &sample_widx, &dosage_present_bits);
          if (pgv.dphase_present[sample_widx] & lowbit) {
            const uint32_t cur_dosage = pgv.dosage_main[dosage_idx];
            const int32_t cur_dphase_delta_val = pgv.dphase_delta[dphase_idx++];
            const uint32_t left_raw = (cur_dosage + cur_dphase_delta_val) >> 1;
            const uint32_t right_raw = (cur_dosage - cur_dphase_delta_val) >> 1;
            const uint32_t left_output_preshift = left_raw * max_output_val + kDosage4th;
            const uint32_t right_output_preshift = right_raw * max_output_val + kDosage4th;
            if ((left_output_preshift ^ right_output_preshift) / kDosageMid) {
              break;
            }
          }
        }
        if (dosage_idx >= pgv.dosage_ct) {
          pgv.dphase_ct = 0;
        }
      }
      const uint32_t use_phased_format = pgv.phasepresent_ct || pgv.dphase_ct;
      if (use_phased_format && pgv.dosage_ct) {
        if (!pgv.phasepresent_ct) {
          ZeroWArr(sample_ctl, pgv.phasepresent);
        } else if (!pgv.dphase_ct) {
          ZeroWArr(sample_ctl, pgv.dphase_present);
        }
      }
      if ((!is_haploid) || (pgv.phasepresent_ct == sample_ct) || (pgv.dphase_ct && UnionIsFull(pgv.phasepresent, pgv.dphase_present, sample_ct))) {
        // 2 alleles, min ploidy == max ploidy == 2
        // (this includes the chrX no-sex-info case)
        // dosages can be patched in after the fact
        AppendU32(0x2020002, &bgen_geno_buf_iter);
        bgen_geno_buf_iter = FillBgen13PloidyAndMissingness(pgv.genovec, pgv.dosage_present, 2, sample_ct, pgv.dosage_ct, bgen_geno_buf_iter);
        const unsigned char* genovec_alias = R_CAST(unsigned char*, pgv.genovec);
        if (!use_phased_format) {
          *bgen_geno_buf_iter++ = 0;
          *bgen_geno_buf_iter++ = bit_precision;
          unsigned char* probs_write_citer = bgen_geno_buf_iter;
          if (!two_byte_probs) {
            for (uint32_t geno_byte_idx = 0; geno_byte_idx != sample_ct4; ++geno_byte_idx) {
              const unsigned char cur_geno4 = genovec_alias[geno_byte_idx];
              memcpy(probs_write_citer, &(bgen_diploid_hardcall_table8[cur_geno4]), sizeof(int64_t));
              probs_write_citer = &(probs_write_citer[bit_precision]);
            }
          } else {
            // 9..16
            for (uint32_t geno_byte_idx = 0; geno_byte_idx != sample_ct4; ++geno_byte_idx) {
              const unsigned char cur_geno4 = genovec_alias[geno_byte_idx];
              memcpy(probs_write_citer, &(bgen_diploid_hardcall_table16[2 * cur_geno4]), 2 * sizeof(int64_t));
              probs_write_citer = &(probs_write_citer[bit_precision]);
            }
          }
        } else {
          // use_phased_format
          *bgen_geno_buf_iter++ = 1;
          *bgen_geno_buf_iter++ = bit_precision;
          const unsigned char* phasepresent_alias = R_CAST(const unsigned char*, pgv.phasepresent);
          const unsigned char* phaseinfo_alias = R_CAST(const unsigned char*, pgv.phaseinfo);
          const uint32_t bit_precision_x2 = bit_precision * 2;

          uintptr_t* probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
          uintptr_t cur_write_bits = 0;
          uint32_t cur_write_bit_idx = 0;
          // possible todo: iterate over genovec halfwords instead, and convert
          // diploid_phased_hardcall_table to support two-at-a-time lookup.
          for (uint32_t geno_byte_idx = 0; geno_byte_idx != sample_ct4; ++geno_byte_idx) {
            const uint32_t geno_byte_idx_d2 = geno_byte_idx / 2;
            uint32_t cur_geno4 = genovec_alias[geno_byte_idx];
            uint32_t cur_phasepresent4 = phasepresent_alias[geno_byte_idx_d2];
            uint32_t cur_phaseinfo4 = phaseinfo_alias[geno_byte_idx_d2];
            if (geno_byte_idx % 2) {
              cur_phasepresent4 >>= 4;
              cur_phaseinfo4 >>= 4;
            }
            for (uint32_t uii = 0; uii != 4; ++uii) {
              const uint32_t cur_index = (cur_geno4 & 3) | ((cur_phasepresent4 & 1) * 4) | ((cur_phaseinfo4 & 1) * 8);
              const uintptr_t payload = bgen_diploid_phased_hardcall_table[cur_index];
              AppendBits(bit_precision_x2, payload, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
              cur_geno4 >>= 2;
              cur_phasepresent4 >>= 1;
              cur_phaseinfo4 >>= 1;
            }
          }
          *probs_write_witer = cur_write_bits;
        }
        if (pgv.dosage_ct) {
          const uintptr_t bit_precision_x2 = bit_precision * 2;
          uintptr_t sample_uidx_base = 0;
          uintptr_t dosage_present_bits = pgv.dosage_present[0];
          uintptr_t* probs_write = R_CAST(uintptr_t*, bgen_geno_buf_iter);
          if (!use_phased_format) {
            // first value is P(geno=2), second value is P(geno=1)
            for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
              // multiply by e.g. 255/16384 and round
              // there are multiple ways to interpret "round to even" here,
              // so I give up and just round 0.5 up (especially since it
              // should be very rare).
              const uintptr_t sample_uidx = BitIter1(pgv.dosage_present, &sample_uidx_base, &dosage_present_bits);
              // note that this needs to be changed to uint64_t if internal
              // dosage representation is no longer 16-bit
              const uint32_t cur_dosage = pgv.dosage_main[dosage_idx];
              uint32_t output_prob2 = 0;
              uint32_t output_prob1;
              if (cur_dosage > kDosageMid) {
                output_prob2 = ((cur_dosage - kDosageMid) * max_output_val + kDosage4th) / kDosageMid;
                output_prob1 = max_output_val - output_prob2;
              } else {
                output_prob1 = (cur_dosage * max_output_val + kDosage4th) / kDosageMid;
              }
              CopyBits(output_prob2 | (output_prob1 << bit_precision), sample_uidx * S_CAST(uint64_t, bit_precision_x2), bit_precision_x2, probs_write);
            }
          } else {
            // phased format, need to take hphase and/or dphase into account
            uint32_t dphase_idx = 0;
            for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
              const uintptr_t sample_uidx = BitIter1(pgv.dosage_present, &sample_uidx_base, &dosage_present_bits);
              const uint32_t cur_dosage = pgv.dosage_main[dosage_idx];
              uint32_t output_prob1;
              uint32_t output_prob2;
              if (IsSet(pgv.dphase_present, sample_uidx)) {
                const int32_t cur_dphase_delta_val = pgv.dphase_delta[dphase_idx++];
                const uint32_t left_raw = (cur_dosage + cur_dphase_delta_val) >> 1;
                const uint32_t right_raw = (cur_dosage - cur_dphase_delta_val) >> 1;
                output_prob1 = (left_raw * max_output_val + kDosage4th) / kDosageMid;
                output_prob2 = (right_raw * max_output_val + kDosage4th) / kDosageMid;
              } else if (IsSet(pgv.phasepresent, sample_uidx)) {
                if (cur_dosage > kDosageMid) {
                  output_prob1 = ((cur_dosage - kDosageMid) * max_output_val + kDosage4th) / kDosageMid;
                  output_prob2 = max_output_val;
                } else {
                  output_prob1 = 0;
                  output_prob2 = (cur_dosage * max_output_val + kDosage4th) / kDosageMid;
                }
                if (IsSet(pgv.phaseinfo, sample_uidx)) {
                  const uint32_t tmpval = output_prob1;
                  output_prob1 = output_prob2;
                  output_prob2 = tmpval;
                }
              } else {
                // unphased
                output_prob1 = (cur_dosage * max_output_val + kDosageMid) / kDosageMax;
                output_prob2 = output_prob1;
              }
              CopyBits(output_prob1 | (output_prob2 << bit_precision), sample_uidx * S_CAST(uint64_t, bit_precision_x2), bit_precision_x2, probs_write);
            }
          }
        }
        const uint64_t tot_prob_bit_ct = sample_ct * 2 * S_CAST(uint64_t, bit_precision);
        bgen_geno_buf_iter = &(bgen_geno_buf_iter[tot_prob_bit_ct / CHAR_BIT]);
        const uint32_t remainder = tot_prob_bit_ct % 8;
        if (remainder) {
          *bgen_geno_buf_iter &= (1 << remainder) - 1;
          ++bgen_geno_buf_iter;
        }
      } else if ((!use_phased_format) && (!is_x) && (!cur_y)) {
        // 2 alleles, min ploidy == max ploidy == 1
        // dosages can be patched in after the fact
        AppendU32(0x1010002, &bgen_geno_buf_iter);
        bgen_geno_buf_iter = FillBgen13PloidyAndMissingness(pgv.genovec, pgv.dosage_present, 1, sample_ct, pgv.dosage_ct, bgen_geno_buf_iter);
        *bgen_geno_buf_iter++ = 0;
        *bgen_geno_buf_iter++ = bit_precision;
        unsigned char* probs_write_citer = bgen_geno_buf_iter;
        const uint16_t* genovec_alias = R_CAST(uint16_t*, pgv.genovec);
        const uint32_t sample_ct8 = DivUp(sample_ct, 8);
        const uint32_t bit_precision_x4 = bit_precision * 4;
        if (!two_byte_probs) {
          for (uint32_t geno_u16_idx = 0; geno_u16_idx != sample_ct8; ++geno_u16_idx) {
            const uint32_t cur_geno8 = genovec_alias[geno_u16_idx];
            const uint64_t new_bytes = bgen_haploid_hardcall_table8[cur_geno8 & 255] | (S_CAST(uint64_t, bgen_haploid_hardcall_table8[cur_geno8 >> 8]) << bit_precision_x4);
            memcpy(probs_write_citer, &new_bytes, sizeof(int64_t));
            probs_write_citer = &(probs_write_citer[bit_precision]);
          }
        } else {
          const uint32_t backfill = 64 - bit_precision_x4;
          for (uint32_t geno_u16_idx = 0; geno_u16_idx != sample_ct8; ++geno_u16_idx) {
            // this may write 14 bytes past the end
            const uint32_t cur_geno8 = genovec_alias[geno_u16_idx];
            uint64_t low_bits = bgen_haploid_hardcall_table16[cur_geno8 & 255];
            uint64_t high_bits = bgen_haploid_hardcall_table16[cur_geno8 >> 8];
            if (backfill) {
              low_bits |= high_bits << bit_precision_x4;
              high_bits = high_bits >> backfill;
            }
            memcpy(probs_write_citer, &low_bits, sizeof(int64_t));
            memcpy(&(probs_write_citer[8]), &high_bits, sizeof(int64_t));
            probs_write_citer = &(probs_write_citer[bit_precision]);
          }
        }
        if (pgv.dosage_ct) {
          uintptr_t sample_uidx_base = 0;
          uintptr_t dosage_present_bits = pgv.dosage_present[0];
          const uint32_t half_int = kDosageMid - (bit_precision == 1);
          uintptr_t* probs_write = R_CAST(uintptr_t*, bgen_geno_buf_iter);
          for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
            // multiply by e.g. 255/32768 and round-to-even
            const uintptr_t sample_uidx = BitIter1(pgv.dosage_present, &sample_uidx_base, &dosage_present_bits);
            // note that this needs to be changed to uint64_t if internal
            // dosage representation is no longer 16-bit
            const uint32_t cur_dosage = pgv.dosage_main[dosage_idx];
            const uint32_t output_prob = (cur_dosage * max_output_val + half_int) / kDosageMax;
            CopyBits(output_prob, sample_uidx * S_CAST(uint64_t, bit_precision), bit_precision, probs_write);
          }
        }
        const uint64_t tot_prob_bit_ct = sample_ct * S_CAST(uint64_t, bit_precision);
        bgen_geno_buf_iter = &(bgen_geno_buf_iter[tot_prob_bit_ct / CHAR_BIT]);
        const uint32_t remainder = tot_prob_bit_ct % 8;
        if (remainder) {
          *bgen_geno_buf_iter &= (1 << remainder) - 1;
          ++bgen_geno_buf_iter;
        }
      } else {
        uintptr_t cur_write_bits = 0;
        uint32_t cur_write_bit_idx = 0;
        uint32_t loop_len = kBitsPerWordD2;
        uintptr_t* probs_write_witer;
        if (!pgv.dosage_ct) {
          if (!pgv.phasepresent_ct) {
            if (is_x) {
              // if male_ct == 0, we used the all-diploid code path instead
              // if male_ct == sample_ct, we treat as haploid general case
              AppendU32(0x2010002, &bgen_geno_buf_iter);
              unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
              *bgen_geno_buf_iter++ = 0;
              *bgen_geno_buf_iter++ = bit_precision;
              probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                uint32_t male_hw = R_CAST(const Halfword*, sex_male_collapsed)[widx];
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
                  const uintptr_t cur_geno = geno_word & 3;
                  const uint32_t cur_male = male_hw & 1;
                  *ploidy_and_missingness_iter++ = (cur_geno == 3) * 128 + 2 - cur_male;
                  if (cur_male) {
                    AppendBits(bit_precision, bgen_haploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                  } else {
                    AppendBits(2 * bit_precision, bgen_diploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                  }
                  geno_word >>= 2;
                  male_hw >>= 1;
                }
              }
            } else {
              assert(cur_y);
              if ((female_ct == sample_ct) && AllBitsAreOne(pgv.genovec, sample_ct * 2)) {
                // chrY, all females, all missing = max ploidy 0
                AppendU32(0x0000002, &bgen_geno_buf_iter);
                bgen_geno_buf_iter = memsetua(bgen_geno_buf_iter, 128, sample_ct);
                *bgen_geno_buf_iter++ = 0;
                *bgen_geno_buf_iter++ = bit_precision;
                probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
              } else {
                AppendU32(0x1000002, &bgen_geno_buf_iter);
                unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
                bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
                *bgen_geno_buf_iter++ = 0;
                *bgen_geno_buf_iter++ = bit_precision;
                probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
                for (uint32_t widx = 0; ; ++widx) {
                  if (widx >= sample_ctl2_m1) {
                    if (widx > sample_ctl2_m1) {
                      break;
                    }
                    loop_len = ModNz(sample_ct, kBitsPerWordD2);
                  }
                  uintptr_t geno_word = pgv.genovec[widx];
                  uint32_t female_hw = R_CAST(const Halfword*, sex_female_collapsed)[widx];
                  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
                    const uintptr_t cur_geno = geno_word & 3;
                    const uint32_t cur_female = female_hw & 1;
                    const uint32_t ploidy_and_missingness = (cur_geno == 3) * (128 - cur_female) + 1;
                    *ploidy_and_missingness_iter++ = ploidy_and_missingness;
                    if (ploidy_and_missingness != 128) {
                      AppendBits(bit_precision, bgen_haploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                    }
                    geno_word >>= 2;
                    female_hw >>= 1;
                  }
                }
              }
            }
          } else {
            // phase present, no dosages present, variable ploidy
            // lookup tables still usable
            if (is_x) {
              AppendU32(0x2010002, &bgen_geno_buf_iter);
              unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
              *bgen_geno_buf_iter++ = 1;
              *bgen_geno_buf_iter++ = bit_precision;
              probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                uint32_t nonmale_hw = ~R_CAST(const Halfword*, sex_male_collapsed)[widx];
                uint32_t phasepresent_hw = R_CAST(const Halfword*, pgv.phasepresent)[widx];
                uint32_t phaseinfo_hw = R_CAST(const Halfword*, pgv.phaseinfo)[widx];
                // todo: pretty sure this can be sped up by pre-shifting
                // phaseinfo_hw by 3.  not sure what else can be done,
                // though...
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits, geno_word >>= 2, nonmale_hw >>= 1, phasepresent_hw >>= 1, phaseinfo_hw >>= 1) {
                  const uintptr_t cur_geno = geno_word & 3;
                  const uint32_t cur_nonmale = nonmale_hw & 1;
                  if (cur_geno == 3) {
                    *ploidy_and_missingness_iter++ = 129 + cur_nonmale;
                    Append0Bits(bit_precision << cur_nonmale, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                    continue;
                  }
                  const uint32_t cur_phasepresent = phasepresent_hw & 1;
                  if (cur_phasepresent | cur_nonmale) {
                    *ploidy_and_missingness_iter++ = 2;
                    const uint32_t cur_index = cur_geno | (cur_phasepresent * 4) | ((phaseinfo_hw & 1) * 8);
                    AppendBits(bit_precision * 2, bgen_diploid_phased_hardcall_table[cur_index], &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                    continue;
                  }
                  *ploidy_and_missingness_iter++ = 1;
                  AppendBits(bit_precision, bgen_haploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                }
              }
            } else {
              // can't be all-female all-missing since phasepresent_ct > 0
              AppendU32(0x2010002 - cur_y * 0x10000, &bgen_geno_buf_iter);
              unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
              *bgen_geno_buf_iter++ = 1;
              *bgen_geno_buf_iter++ = bit_precision;
              probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
              uint32_t female_hw = 0;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                if (cur_y) {
                  female_hw = R_CAST(const Halfword*, sex_female_collapsed)[widx];
                }
                uint32_t phasepresent_hw = R_CAST(const Halfword*, pgv.phasepresent)[widx];
                uint32_t phaseinfo_hw = R_CAST(const Halfword*, pgv.phaseinfo)[widx];
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits, geno_word >>= 2, phasepresent_hw >>= 1, phaseinfo_hw >>= 1) {
                  const uintptr_t cur_geno = geno_word & 3;
                  if (cur_geno == 3) {
                    // female_hw not incrementally shifted for now
                    if (!(female_hw & (1U << sample_idx_lowbits))) {
                      *ploidy_and_missingness_iter++ = 129;
                      Append0Bits(bit_precision, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                    } else {
                      *ploidy_and_missingness_iter++ = 128;
                    }
                    continue;
                  }
                  const uint32_t cur_phasepresent = phasepresent_hw & 1;
                  if (cur_phasepresent) {
                    *ploidy_and_missingness_iter++ = 2;
                    const uint32_t cur_index = cur_geno | (cur_phasepresent * 4) | ((phaseinfo_hw & 1) * 8);
                    AppendBits(bit_precision * 2, bgen_diploid_phased_hardcall_table[cur_index], &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                    continue;
                  }
                  *ploidy_and_missingness_iter++ = 1;
                  AppendBits(bit_precision, bgen_haploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                }
              }
            }
          }
        } else {
          // dosages present, variable ploidy
          // don't bother with lookup tables here
          const Dosage* dosage_main_iter = pgv.dosage_main;
          if (!use_phased_format) {
            if (is_x) {
              AppendU32(0x2010002, &bgen_geno_buf_iter);
              unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
              *bgen_geno_buf_iter++ = 0;
              *bgen_geno_buf_iter++ = bit_precision;
              probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                const uint32_t male_hw = R_CAST(const Halfword*, sex_male_collapsed)[widx];
                const uint32_t dosage_present_hw = R_CAST(const Halfword*, pgv.dosage_present)[widx];
                uint32_t shifted_bit = 1;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits, shifted_bit *= 2, geno_word >>= 2) {
                  const uint32_t cur_male = male_hw & shifted_bit;
                  uint32_t cur_dosage;
                  if (dosage_present_hw & shifted_bit) {
                    cur_dosage = *dosage_main_iter++;
                  } else {
                    const uintptr_t cur_geno = geno_word & 3;
                    if (cur_geno == 3) {
                      const uint32_t cur_nonmale = !cur_male;
                      *ploidy_and_missingness_iter++ = 129 + cur_nonmale;
                      Append0Bits(bit_precision << cur_nonmale, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                      continue;
                    }
                    cur_dosage = cur_geno * kDosageMid;
                  }
                  if (cur_male) {
                    *ploidy_and_missingness_iter++ = 1;
                    uint32_t output_prob = (cur_dosage * max_output_val + kDosageMid) / kDosageMax;
                    AppendBits(bit_precision, output_prob, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                    continue;
                  }
                  *ploidy_and_missingness_iter++ = 2;
                  uint32_t output_prob2 = 0;
                  uint32_t output_prob1;
                  if (cur_dosage > kDosageMid) {
                    output_prob2 = ((cur_dosage - kDosageMid) * max_output_val + kDosage4th) / kDosageMid;
                    output_prob1 = max_output_val - output_prob2;
                  } else {
                    output_prob1 = (cur_dosage * max_output_val + kDosage4th) / kDosageMid;
                  }
                  AppendBits(bit_precision * 2, output_prob2 | (output_prob1 << bit_precision), &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                }
              }
            } else {
              assert(cur_y);
              // can't be all-female all-missing since dosage_ct > 0
              AppendU32(0x1000002, &bgen_geno_buf_iter);
              unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
              *bgen_geno_buf_iter++ = 0;
              *bgen_geno_buf_iter++ = bit_precision;
              probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                const uint32_t female_hw = R_CAST(const Halfword*, sex_female_collapsed)[widx];
                const uint32_t dosage_present_hw = R_CAST(const Halfword*, pgv.dosage_present)[widx];
                uint32_t shifted_bit = 1;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits, shifted_bit *= 2, geno_word >>= 2) {
                  uint32_t cur_dosage;
                  if (dosage_present_hw & shifted_bit) {
                    cur_dosage = *dosage_main_iter++;
                  } else {
                    const uintptr_t cur_geno = geno_word & 3;
                    if (cur_geno == 3) {
                      if (!(female_hw & shifted_bit)) {
                        *ploidy_and_missingness_iter++ = 129;
                        Append0Bits(bit_precision, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                      } else {
                        *ploidy_and_missingness_iter++ = 128;
                      }
                      continue;
                    }
                    cur_dosage = cur_geno * kDosageMid;
                  }
                  *ploidy_and_missingness_iter++ = 1;
                  uint32_t output_prob = (cur_dosage * max_output_val + kDosageMid) / kDosageMax;
                  AppendBits(bit_precision, output_prob, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                }
              }
            }
          } else {
            const SDosage* dphase_delta_iter = pgv.dphase_delta;
            if (is_x) {
              AppendU32(0x2010002, &bgen_geno_buf_iter);
              unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
              *bgen_geno_buf_iter++ = 1;
              *bgen_geno_buf_iter++ = bit_precision;
              probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                const uint32_t male_hw = R_CAST(const Halfword*, sex_male_collapsed)[widx];
                const uint32_t phaseinfo_hw = R_CAST(const Halfword*, pgv.phaseinfo)[widx];
                const uint32_t dosage_present_hw = R_CAST(const Halfword*, pgv.dosage_present)[widx];
                const uint32_t dphase_present_hw = R_CAST(const Halfword*, pgv.dphase_present)[widx];
                const uint32_t either_phase_present_hw = dphase_present_hw | (R_CAST(const Halfword*, pgv.phasepresent)[widx]);
                uint32_t shifted_bit = 1;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits, shifted_bit *= 2, geno_word >>= 2) {
                  const uint32_t cur_male = male_hw & shifted_bit;
                  uint32_t cur_dosage;
                  if (dosage_present_hw & shifted_bit) {
                    cur_dosage = *dosage_main_iter++;
                  } else {
                    const uintptr_t cur_geno = geno_word & 3;
                    if (cur_geno == 3) {
                      const uint32_t cur_nonmale = !cur_male;
                      *ploidy_and_missingness_iter++ = 129 + cur_nonmale;
                      Append0Bits(bit_precision << cur_nonmale, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                      continue;
                    }
                    cur_dosage = cur_geno * kDosageMid;
                  }
                  uint32_t output_prob1;
                  uint32_t output_prob2;
                  if (either_phase_present_hw & shifted_bit) {
                    if (dphase_present_hw & shifted_bit) {
                      const int32_t cur_dphase_delta_val = *dphase_delta_iter++;
                      const uint32_t left_raw = (cur_dosage + cur_dphase_delta_val) >> 1;
                      const uint32_t right_raw = (cur_dosage - cur_dphase_delta_val) >> 1;
                      output_prob1 = (left_raw * max_output_val + kDosage4th) / kDosageMid;
                      output_prob2 = (right_raw * max_output_val + kDosage4th) / kDosageMid;
                    } else {
                      if (cur_dosage > kDosageMid) {
                        output_prob1 = ((cur_dosage - kDosageMid) * max_output_val + kDosage4th) / kDosageMid;
                        output_prob2 = max_output_val;
                      } else {
                        output_prob1 = 0;
                        output_prob2 = (cur_dosage * max_output_val + kDosage4th) / kDosageMid;
                      }
                      if (phaseinfo_hw & shifted_bit) {
                        const uint32_t tmpval = output_prob1;
                        output_prob1 = output_prob2;
                        output_prob2 = tmpval;
                      }
                    }
                  } else {
                    // unphased, nonmissing
                    output_prob1 = (cur_dosage * max_output_val + kDosageMid) / kDosageMax;
                    if (cur_male) {
                      *ploidy_and_missingness_iter++ = 1;
                      AppendBits(bit_precision, output_prob1, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                      continue;
                    }
                    output_prob2 = output_prob1;
                  }
                  *ploidy_and_missingness_iter++ = 2;
                  AppendBits(bit_precision * 2, output_prob1 | (output_prob2 << bit_precision), &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                }
              }
            } else {
              // can't be all-female all-missing since dosage_ct > 0
              AppendU32(0x2010002 - cur_y * 0x10000, &bgen_geno_buf_iter);
              unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
              *bgen_geno_buf_iter++ = 1;
              *bgen_geno_buf_iter++ = bit_precision;
              probs_write_witer = R_CAST(uintptr_t*, bgen_geno_buf_iter);
              uint32_t female_hw = 0;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                if (cur_y) {
                  female_hw = R_CAST(const Halfword*, sex_female_collapsed)[widx];
                }
                const uint32_t phaseinfo_hw = R_CAST(const Halfword*, pgv.phaseinfo)[widx];
                const uint32_t dosage_present_hw = R_CAST(const Halfword*, pgv.dosage_present)[widx];
                const uint32_t dphase_present_hw = R_CAST(const Halfword*, pgv.dphase_present)[widx];
                const uint32_t either_phase_present_hw = dphase_present_hw | (R_CAST(const Halfword*, pgv.phasepresent)[widx]);
                uint32_t shifted_bit = 1;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits, shifted_bit *= 2, geno_word >>= 2) {
                  uint32_t cur_dosage;
                  if (dosage_present_hw & shifted_bit) {
                    cur_dosage = *dosage_main_iter++;
                  } else {
                    const uintptr_t cur_geno = geno_word & 3;
                    if (cur_geno == 3) {
                      if (!(female_hw & shifted_bit)) {
                        *ploidy_and_missingness_iter++ = 129;
                        Append0Bits(bit_precision, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                      } else {
                        *ploidy_and_missingness_iter++ = 128;
                      }
                      continue;
                    }
                    cur_dosage = cur_geno * kDosageMid;
                  }
                  if (either_phase_present_hw & shifted_bit) {
                    *ploidy_and_missingness_iter++ = 2;
                    uint32_t output_prob1;
                    uint32_t output_prob2;
                    if (dphase_present_hw & shifted_bit) {
                      const int32_t cur_dphase_delta_val = *dphase_delta_iter++;
                      const uint32_t left_raw = (cur_dosage + cur_dphase_delta_val) >> 1;
                      const uint32_t right_raw = (cur_dosage - cur_dphase_delta_val) >> 1;
                      output_prob1 = (left_raw * max_output_val + kDosage4th) / kDosageMid;
                      output_prob2 = (right_raw * max_output_val + kDosage4th) / kDosageMid;
                    } else {
                      if (cur_dosage > kDosageMid) {
                        output_prob1 = ((cur_dosage - kDosageMid) * max_output_val + kDosage4th) / kDosageMid;
                        output_prob2 = max_output_val;
                      } else {
                        output_prob1 = 0;
                        output_prob2 = (cur_dosage * max_output_val + kDosage4th) / kDosageMid;
                      }
                      if (phaseinfo_hw & shifted_bit) {
                        const uint32_t tmpval = output_prob1;
                        output_prob1 = output_prob2;
                        output_prob2 = tmpval;
                      }
                    }
                    AppendBits(bit_precision * 2, output_prob1 | (output_prob2 << bit_precision), &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                    continue;
                  }
                  *ploidy_and_missingness_iter++ = 1;
                  uint32_t output_prob = (cur_dosage * max_output_val + kDosageMid) / kDosageMax;
                  AppendBits(bit_precision, output_prob, &cur_write_bits, &cur_write_bit_idx, &probs_write_witer);
                }
              }
            }
          }
        }
        *probs_write_witer = cur_write_bits;
        bgen_geno_buf_iter = R_CAST(unsigned char*, probs_write_witer);
        bgen_geno_buf_iter = &(bgen_geno_buf_iter[DivUp(cur_write_bit_idx, CHAR_BIT)]);
      }
      const uint32_t uncompressed_bytect = bgen_geno_buf_iter - uncompressed_bgen_geno_buf;
      /*
      for (uint32_t uii = 0; uii != uncompressed_bytect; ++uii) {
        printf("%u ", uncompressed_bgen_geno_buf[uii]);
      }
      printf("\n");
      exit(1);
      */
      uintptr_t compressed_bytect;
      if (compressor) {
        compressed_bytect = libdeflate_zlib_compress(compressor, uncompressed_bgen_geno_buf, uncompressed_bytect, writebuf_iter, bgen_compressed_buf_max);
      } else {
        compressed_bytect = ZSTD_compress(writebuf_iter, bgen_compressed_buf_max, uncompressed_bgen_geno_buf, uncompressed_bytect, zst_level);
      }
      assert(compressed_bytect);
      *variant_bytect_iter++ = 4 + compressed_bytect;
      *variant_bytect_iter++ = uncompressed_bytect;
      writebuf_iter = &(writebuf_iter[bgen_compressed_buf_max]);
      GenoarrToMissingnessUnsafe(pgv.genovec, sample_ct, missing_acc1);
      if (pgv.dosage_ct) {
        BitvecInvmask(pgv.dosage_present, sample_ctl, missing_acc1);
      }
      if (is_y) {
        BitvecAnd(sex_male_collapsed, sample_ctl, missing_acc1);
        // ignore missingness state of unknown-sex chrY calls for now
      }
      VcountIncr1To4(missing_acc1, acc1_vec_ct, missing_acc4);
      if (!(--vidx_rem15)) {
        Vcount0Incr4To8(acc4_vec_ct, missing_acc4, missing_acc8);
        vidx_rem15 = 15;
        if (!(--vidx_rem255d15)) {
          Vcount0Incr8To32(acc8_vec_ct, missing_acc8, missing_acc32);
          vidx_rem255d15 = 17;
        }
      }
    }
    while (0) {
    ExportBgen13Thread_err:
      UpdateU64IfSmaller(new_err_info, &ctx->err_info);
      break;
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  VcountIncr4To8(missing_acc4, acc4_vec_ct, missing_acc8);
  VcountIncr8To32(missing_acc8, acc8_vec_ct, missing_acc32);
  THREAD_RETURN;
}

// This allocates exported_sample_ids on top.
BoolErr ExportIdpaste(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* ftypename, uint32_t sample_ct, IdpasteFlags exportf_id_paste, char exportf_id_delim, uintptr_t* max_exported_sample_id_blen_ptr, char** exported_sample_ids_ptr, uint32_t** exported_id_htable_ptr) {
  const uint32_t write_fid = DataFidColIsRequired(sample_include, siip, sample_ct, exportf_id_paste / kfIdpasteMaybefid);
  const char* sample_ids = siip->sample_ids;
  const char* sids = siip->sids;
  const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
  uintptr_t max_sid_blen = siip->max_sid_blen;
  const uint32_t write_sid = DataSidColIsRequired(sample_include, sids, sample_ct, max_sid_blen, exportf_id_paste / kfIdpasteMaybesid);
  if (write_sid && (!sids)) {
    max_sid_blen = 2;
  }
  uint32_t id_delim_warning = 0;
  char id_delim = exportf_id_delim? exportf_id_delim : '_';
  const uintptr_t max_exported_sample_id_blen = max_sample_id_blen + write_sid * max_sid_blen;
  const uint32_t exported_id_htable_size = GetHtableMinSize(sample_ct);
  // check for duplicates
  if (unlikely(
          bigstack_alloc_c(sample_ct * max_exported_sample_id_blen, exported_sample_ids_ptr) ||
          bigstack_alloc_u32(exported_id_htable_size, exported_id_htable_ptr))) {
    return 1;
  }
  *max_exported_sample_id_blen_ptr = max_exported_sample_id_blen;
  char* exported_sample_ids = *exported_sample_ids_ptr;
  uint32_t* exported_id_htable = *exported_id_htable_ptr;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    const char* orig_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
    const char* orig_fid_end = AdvToDelim(orig_sample_id, '\t');
    char* exported_sample_ids_iter = &(exported_sample_ids[sample_idx * max_exported_sample_id_blen]);
    if (write_fid) {
      const uint32_t fid_slen = orig_fid_end - orig_sample_id;
      if ((!id_delim_warning) && memchr(orig_sample_id, id_delim, fid_slen)) {
        id_delim_warning = 1;
      }
      exported_sample_ids_iter = memcpyax(exported_sample_ids_iter, orig_sample_id, fid_slen, id_delim);
    }
    if (exportf_id_paste & kfIdpasteIid) {
      const char* orig_iid = &(orig_fid_end[1]);
      const uint32_t iid_slen = strlen(orig_iid);
      if ((!id_delim_warning) && memchr(orig_iid, id_delim, iid_slen)) {
        id_delim_warning = 1;
      }
      exported_sample_ids_iter = memcpyax(exported_sample_ids_iter, orig_iid, iid_slen, id_delim);
    }
    if (write_sid) {
      if (sids) {
        const char* orig_sid = &(sids[sample_uidx * max_sid_blen]);
        const uint32_t sid_slen = strlen(orig_sid);
        if ((!id_delim_warning) && memchr(orig_sid, id_delim, sid_slen)) {
          id_delim_warning = 1;
        }
        exported_sample_ids_iter = memcpya(exported_sample_ids_iter, orig_sid, sid_slen);
      } else {
        *exported_sample_ids_iter++ = '0';
      }
      ++exported_sample_ids_iter;
    }
    exported_sample_ids_iter[-1] = '\0';
  }
  // todo: revise this warning condition?
  if (id_delim_warning) {
    if (exportf_id_delim) {
      logerrprintfww("Warning: '%c' present in original sample IDs; --export %s will not be able to reconstruct them. Consider rerunning with a different --export id-delim= value.\n", exportf_id_delim, ftypename);
    } else {
      logerrprintfww("Warning: '_' present in original sample IDs; --export %s will not be able to reconstruct them. Consider rerunning with a suitable --export id-delim= value.\n", ftypename);
    }
  }
  if (PopulateStrboxHtable(exported_sample_ids, sample_ct, max_exported_sample_id_blen, exported_id_htable_size, exported_id_htable)) {
    logerrprintfww("Warning: Duplicate sample ID(s) are being written to --export %s file.\n", ftypename);
  }
  return 0;
}

PglErr ExportBgen13(const char* outname, const uintptr_t* sample_include, uint32_t* sample_include_cumulative_popcounts, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, ExportfFlags exportf_flags, uint32_t exportf_bits, IdpasteFlags exportf_id_paste, char exportf_id_delim, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uint32_t* sample_missing_geno_cts) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  ExportBgen13Ctx ctx;
  {
    const uint32_t use_zstd_compression = !(exportf_flags & kfExportfBgen12);
    if (use_zstd_compression) {
      ctx.libdeflate_compressors = nullptr;
    } else {
      if (unlikely(BIGSTACK_ALLOC_X(struct libdeflate_compressor*, max_thread_ct, &ctx.libdeflate_compressors))) {
        goto ExportBgen13_ret_NOMEM;
      }
      ZeroPtrArr(max_thread_ct, ctx.libdeflate_compressors);
    }
    const uint32_t phase_is_present = pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent);
    if (!exportf_bits) {
      // default
      exportf_bits = 16;
    } else if (phase_is_present) {
      if (exportf_bits < 15) {
        if (exportf_bits == 1) {
          logerrputs("Warning: Unphased heterozygous calls in partially-phased variants cannot be\nexported with bits=1.\n");
        } else {
          logerrprintf("Warning: Unphased heterozygous hardcalls in partially-phased variants are\npoorly represented with bits=%u.\n", exportf_bits);
          const double suggested_thresh = 0.00009375 * (16384 >> exportf_bits);
          if (exportf_bits < 4) {
            logerrprintfww("It is necessary to use e.g. --hard-call-threshold %g + --dosage-erase-threshold %g to re-import them cleanly.\n", suggested_thresh, suggested_thresh);
          } else {
            logerrprintfww("It is necessary to use e.g. --dosage-erase-threshold %g to re-import them cleanly.\n", suggested_thresh);
          }
        }
      }
    }
    if (unlikely(ConstructBgen13LookupTables(exportf_bits, &ctx.tables))) {
      goto ExportBgen13_ret_NOMEM;
    }
    const uint32_t max_chr_slen = GetMaxChrSlen(cip);
    uintptr_t bgen_geno_cacheline_ct;
    uintptr_t bgen_compressed_buf_max;
    {
      const uint64_t tot_prob_bit_ct = S_CAST(uint64_t, sample_ct) * 2 * exportf_bits;
      uint64_t bgen_geno_buf_size = 10 + sample_ct + DivUp(tot_prob_bit_ct, 8);
      if (unlikely(bgen_geno_buf_size > 0xffffffffU - 4)) {
        // could return VarRecordTooLarge error instead
        logerrputs("Error: Too many samples for .bgen format.\n");
        goto ExportBgen13_ret_INCONSISTENT_INPUT;
      }
      if (!use_zstd_compression) {
        bgen_compressed_buf_max = libdeflate_deflate_compress_bound(nullptr, bgen_geno_buf_size);
      } else {
        bgen_compressed_buf_max = ZSTD_compressBound(bgen_geno_buf_size);
      }
      // +14 since we may write that far past the end
      bgen_geno_buf_size += 14;
      bgen_geno_cacheline_ct = DivUp(bgen_geno_buf_size, kCacheline);
    }
    ctx.bgen_compressed_buf_max = bgen_compressed_buf_max;
    // When writing sample ID block, we flush-check after each sample ID.
    // When writing variant data blocks, we flush-check at the beginning of
    // each allele, and at the beginning and end of each genotype data block.
    uintptr_t writebuf_len = 16 + kMaxIdSlen + max_chr_slen;
    if (writebuf_len < max_allele_slen + 4) {
      writebuf_len = max_allele_slen + 4;
    }
    if (writebuf_len < bgen_compressed_buf_max) {
      writebuf_len = bgen_compressed_buf_max;
    }
    writebuf_len += kMaxMediumLine;
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    char* chr_buf;
    unsigned char* writebuf;
    if (unlikely(
            bigstack_alloc_c(max_chr_slen, &chr_buf) ||
            bigstack_alloc_uc(writebuf_len, &writebuf) ||
            bigstack_alloc_w(sample_ctl, &ctx.sex_male_collapsed) ||
            bigstack_alloc_w(sample_ctl, &ctx.sex_female_collapsed))) {
      goto ExportBgen13_ret_NOMEM;
    }
    CopyBitarrSubset(sex_male, sample_include, sample_ct, ctx.sex_male_collapsed);
    CopyBitarrSubset(sex_nm, sample_include, sample_ct, ctx.sex_female_collapsed);
    BitvecInvmask(ctx.sex_male_collapsed, sample_ctl, ctx.sex_female_collapsed);

    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportBgen13_ret_OPEN_FAIL;
    }
    // bgen 1.2-1.3 header
    // note that \xxx character constants are interpreted in octal, so \200 is
    // decimal 128, etc.
    unsigned char* write_iter = &(writebuf[4]);
    AppendU32(20, &write_iter);
    AppendU32(variant_ct, &write_iter);
    AppendU32(sample_ct, &write_iter);
    write_iter = memcpyua(write_iter, "bgen\0\0\0\200", 8);
    // compression mode (1 + use_zstd_compression), layout=2
    writebuf[20] = 9 + use_zstd_compression;
    unsigned char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    // always save sample IDs...
    char* exported_sample_ids;
    uint32_t* exported_id_htable;
    uintptr_t max_exported_sample_id_blen;
    if (unlikely(ExportIdpaste(sample_include, siip, use_zstd_compression? "bgen-1.3" : "bgen-1.2", sample_ct, exportf_id_paste, exportf_id_delim, &max_exported_sample_id_blen, &exported_sample_ids, &exported_id_htable))) {
      goto ExportBgen13_ret_NOMEM;
    }
    // Compute total length of sample ID block now; this is necessary to fill
    // in bytes 0-3 and 24-27 correctly; and we must do this now since we may
    // flush them before we're done rendering the block.
    uintptr_t sample_id_block_len = 2 * sample_ct + 8;
    const char* exported_sample_ids_iter = exported_sample_ids;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      sample_id_block_len += strlen(exported_sample_ids_iter);
      exported_sample_ids_iter = &(exported_sample_ids_iter[max_exported_sample_id_blen]);
    }
#ifdef __LP64__
    if (sample_id_block_len > 0xffffffffU - 20) {
      // ...unless combined sample ID length is actually greater than 4 GiB, in
      // which case bgen-1.2/1.3 can't actually represent it.
      logerrputs("Warning: Omitting sample ID block from .bgen file since it would overflow (more\nthan 4 GiB).  Consider using shorter IDs.\n");
      memcpy(writebuf, "\24\0\0", 4);
    } else {
#endif
      uint32_t initial_bgen_offset = sample_id_block_len + 20;
      memcpy(writebuf, &initial_bgen_offset, 4);
      AppendU32(sample_id_block_len, &write_iter);
      AppendU32(sample_ct, &write_iter);
      exported_sample_ids_iter = exported_sample_ids;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uint16_t cur_slen = strlen(exported_sample_ids_iter);
        AppendU16(cur_slen, &write_iter);
        write_iter = memcpyua(write_iter, exported_sample_ids_iter, cur_slen);
        exported_sample_ids_iter = &(exported_sample_ids_iter[max_exported_sample_id_blen]);
        if (unlikely(fwrite_uflush2(writebuf_flush, outfile, &write_iter))) {
          goto ExportBgen13_ret_WRITE_FAIL;
        }
      }
#ifdef __LP64__
    }
#endif
    BigstackReset(exported_sample_ids);

    const uintptr_t max_write_block_byte_ct = bigstack_left() / 4;
    uint32_t max_write_block_size = kPglVblockSize;
    for (; ; max_write_block_size /= 2) {
      // limit each write buffer to 1/4 of remaining workspace
      if ((S_CAST(uint64_t, bgen_compressed_buf_max + 2 * sizeof(int32_t))) * max_write_block_size <= max_write_block_byte_ct) {
        break;
      }
      // 5 bytes per sample * 500k samples = ~2.5 MB per variant; with
      // max_write_block_size lower limit of 256, minimum workspace is ~2.5 GB.
      // That's perfectly reasonable (and it's okay if PgenMtLoadInit raises
      // the requirement to, say, 4 times that).
      // Memory will get tighter once multiallelic variants are supported,
      // though.
      if (unlikely(max_write_block_size <= kBitsPerVec)) {
        goto ExportBgen13_ret_NOMEM;
      }
    }

    // todo: test when this saturates
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct > max_write_block_size) {
      calc_thread_ct = max_write_block_size;
    }

    if (unlikely(
            bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(ctx.writebufs[0])) ||
            bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(ctx.writebufs[1])) ||
            bigstack_alloc_u32(max_write_block_size * 2, &(ctx.variant_bytects[0])) ||
            bigstack_alloc_u32(max_write_block_size * 2, &(ctx.variant_bytects[1])) ||
            bigstack_alloc_wp(calc_thread_ct, &ctx.missing_acc1) ||
            bigstack_alloc_ucp(calc_thread_ct, &ctx.uncompressed_bgen_geno_bufs))) {
      goto ExportBgen13_ret_NOMEM;
    }

    const uint32_t acc1_vec_ct = BitCtToVecCt(sample_ct);
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    const uintptr_t track_missing_cacheline_ct = VecCtToCachelineCt(acc1_vec_ct * 45);
    const uintptr_t thread_xalloc_cacheline_ct = track_missing_cacheline_ct + bgen_geno_cacheline_ct;
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    ctx.phasepresents = nullptr;
    ctx.phaseinfos = nullptr;
    ctx.dosage_presents = nullptr;
    ctx.dosage_mains = nullptr;
    ctx.dphase_presents = nullptr;
    ctx.dphase_deltas = nullptr;
    uint32_t read_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct, raw_variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &ctx.genovecs, nullptr, phase_is_present? (&ctx.phasepresents) : nullptr, phase_is_present? (&ctx.phaseinfos) : nullptr, dosage_is_present? (&ctx.dosage_presents) : nullptr, dosage_is_present? (&ctx.dosage_mains) : nullptr, phase_is_present? (&ctx.dphase_presents) : nullptr, phase_is_present? (&ctx.dphase_deltas) : nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto ExportBgen13_ret_NOMEM;
    }
    // Ran this with and without --memory 640 on a 1000G phase 1 dataset
    // (chosen since it has dosages), and there was zero performance penalty in
    // the former case.
    //
    // So it probably makes sense to default to a much smaller chunk size, and
    // use the same type of thread-pool logic as the new bgzf writer.
    if (read_block_size > max_write_block_size) {
      read_block_size = max_write_block_size;
    }

    const uint32_t ref_allele_last = !(exportf_flags & kfExportfRefFirst);
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      ctx.missing_acc1[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(track_missing_cacheline_ct * kCacheline));
      ctx.uncompressed_bgen_geno_bufs[tidx] = S_CAST(unsigned char*, bigstack_alloc_raw(bgen_geno_cacheline_ct * kCacheline));
    }
    ctx.sample_ct = sample_ct;
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.sample_include = sample_include;
    ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto ExportBgen13_ret_NOMEM;
    }
    ctx.refalt1_select = refalt1_select;
    ctx.ref_allele_last = ref_allele_last;
    ctx.cip = cip;
    ctx.err_info = (~0LLU) << 32;
    if (!use_zstd_compression) {
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        ctx.libdeflate_compressors[tidx] = libdeflate_alloc_compressor(6);
        if (unlikely(!ctx.libdeflate_compressors[tidx])) {
          goto ExportBgen13_ret_NOMEM;
        }
      }
    }
    SetThreadFuncAndData(ExportBgen13Thread, &ctx, &tg);

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
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_slen = 0;

    uint32_t prev_block_write_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    logprintfww5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uint32_t allele_ct = 2;
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto ExportBgen13_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, ctx.err_info);
        if (unlikely(reterr)) {
          if (reterr == kPglRetInconsistentInput) {
            // temporary
            logputs("\n");
            logerrputs("Error: --export bgen-1.2/1.3 multiallelic variant support is under development.\n");
          }
          goto ExportBgen13_ret_PGR_FAIL;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_write_ct = cur_block_write_ct;
        ComputeUidxStartPartition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, ctx.read_variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_write_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto ExportBgen13_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const unsigned char* compressed_data_iter = ctx.writebufs[parity];
        const uint32_t* variant_bytect_iter = ctx.variant_bytects[parity];
        for (uint32_t variant_bidx = 0; variant_bidx != prev_block_write_ct; ++variant_bidx) {
          const uint32_t write_variant_uidx = BitIter1(variant_include, &write_variant_uidx_base, &cur_bits);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            const char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
            chr_slen = chr_name_end - chr_buf;
          }
          const char* cur_variant_id = variant_ids[write_variant_uidx];
          const uint32_t id_slen = strlen(cur_variant_id);
          // low 16 bits = null "SNP ID"
          AppendU32(id_slen << 16, &write_iter);
          write_iter = memcpyua(write_iter, cur_variant_id, id_slen);
          AppendU16(chr_slen, &write_iter);
          write_iter = memcpyua(write_iter, chr_buf, chr_slen);
          AppendU32(variant_bps[write_variant_uidx], &write_iter);
          uintptr_t allele_idx_offset_base = write_variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[write_variant_uidx];
            allele_ct = allele_idx_offsets[write_variant_uidx + 1] - allele_idx_offset_base;
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          if (refalt1_select) {
            ref_allele_idx = refalt1_select[write_variant_uidx][0];
            alt1_allele_idx = refalt1_select[write_variant_uidx][1];
          }
          AppendU16(allele_ct, &write_iter);
          const char* ref_allele = cur_alleles[ref_allele_idx];
          const uint32_t ref_allele_slen = strlen(ref_allele);
          if (!ref_allele_last) {
            if (unlikely(fwrite_uflush2(writebuf_flush, outfile, &write_iter))) {
              goto ExportBgen13_ret_WRITE_FAIL;
            }
            AppendU32(ref_allele_slen, &write_iter);
            write_iter = memcpyua(write_iter, ref_allele, ref_allele_slen);
          }
          if (unlikely(fwrite_uflush2(writebuf_flush, outfile, &write_iter))) {
            goto ExportBgen13_ret_WRITE_FAIL;
          }
          const char* cur_alt_allele = cur_alleles[alt1_allele_idx];
          uint32_t alt_allele_slen = strlen(cur_alt_allele);
          AppendU32(alt_allele_slen, &write_iter);
          write_iter = memcpyua(write_iter, cur_alt_allele, alt_allele_slen);
          if (allele_ct > 2) {
            for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
              if ((allele_idx == ref_allele_idx) || (allele_idx == alt1_allele_idx)) {
                continue;
              }
              if (unlikely(fwrite_uflush2(writebuf_flush, outfile, &write_iter))) {
                goto ExportBgen13_ret_WRITE_FAIL;
              }
              cur_alt_allele = cur_alleles[allele_idx];
              alt_allele_slen = strlen(cur_alt_allele);
              AppendU32(alt_allele_slen, &write_iter);
              write_iter = memcpyua(write_iter, cur_alt_allele, alt_allele_slen);
            }
          }
          if (ref_allele_last) {
            if (unlikely(fwrite_uflush2(writebuf_flush, outfile, &write_iter))) {
              goto ExportBgen13_ret_WRITE_FAIL;
            }
            AppendU32(ref_allele_slen, &write_iter);
            write_iter = memcpyua(write_iter, ref_allele, ref_allele_slen);
          }
          // compressed data length, then uncompressed data length
          write_iter = memcpyua(write_iter, variant_bytect_iter, 8);
          const uint32_t cur_compressed_bytect = *variant_bytect_iter - 4;
          variant_bytect_iter = &(variant_bytect_iter[2]);
          if (unlikely(fwrite_uflush2(writebuf_flush, outfile, &write_iter))) {
            goto ExportBgen13_ret_WRITE_FAIL;
          }
          // may want to elide this memcpy when possible
          write_iter = memcpyua(write_iter, compressed_data_iter, cur_compressed_bytect);
          compressed_data_iter = &(compressed_data_iter[bgen_compressed_buf_max]);
          if (unlikely(fwrite_uflush2(writebuf_flush, outfile, &write_iter))) {
            goto ExportBgen13_ret_WRITE_FAIL;
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
      prev_block_write_ct = cur_block_write_ct;
      variant_idx += cur_block_write_ct;
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(fclose_uflush_null(writebuf_flush, write_iter, &outfile))) {
      goto ExportBgen13_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    const uint32_t sample_ctav = acc1_vec_ct * kBitsPerVec;
    const uintptr_t acc32_offset = acc1_vec_ct * (13 * k1LU * kWordsPerVec);
    uint32_t* scrambled_missing_cts = R_CAST(uint32_t*, &(ctx.missing_acc1[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
      const uint32_t* thread_scrambled_missing_cts = R_CAST(uint32_t*, &(ctx.missing_acc1[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii != sample_ctav; ++uii) {
        scrambled_missing_cts[uii] += thread_scrambled_missing_cts[uii];
      }
    }
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = VcountScramble1(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  ExportBgen13_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportBgen13_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportBgen13_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ExportBgen13_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportBgen13_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ExportBgen13_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  if (ctx.libdeflate_compressors) {
    for (uint32_t tidx = 0; tidx != max_thread_ct; ++tidx) {
      if (!ctx.libdeflate_compressors[tidx]) {
        break;
      }
      libdeflate_free_compressor(ctx.libdeflate_compressors[tidx]);
    }
  }
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr ExportOxSample(const char* outname, const uintptr_t* sample_include, const char* sample_ids, const uint32_t* sample_missing_geno_cts, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t variant_ct, uint32_t y_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t pheno_ctl = BitCtToWordCt(pheno_ct);
    char* writebuf;
    uintptr_t* is_basic_categorical;
    // if phenotype is categorical, and all (non-null) category names are of
    // the form P[positive integer], then it's best to emit the positive
    // integer in the name string instead of the internal index.
    if (unlikely(
            bigstack_calloc_w(pheno_ctl, &is_basic_categorical) ||
            bigstack_alloc_c(kMaxMediumLine + max_sample_id_blen + 32 + pheno_ct * MAXV(kMaxMissingPhenostrBlen, 16), &writebuf))) {
      goto ExportOxSample_ret_NOMEM;
    }

    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportOxSample_ret_OPEN_FAIL;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = strcpya_k(writebuf, "ID_1 ID_2 missing sex");
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      *write_iter++ = ' ';
      write_iter = strcpya(write_iter, &(pheno_names[pheno_idx * max_pheno_name_blen]));
      const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
      if (cur_pheno_col->type_code == kPhenoDtypeCat) {
        const uint32_t nn_cat_ct = cur_pheno_col->nonnull_category_ct;
        const char* const* cur_cat_names = cur_pheno_col->category_names;
        uint32_t cat_idx;
        for (cat_idx = 1; cat_idx <= nn_cat_ct; ++cat_idx) {
          const char* cat_name_iter = cur_cat_names[cat_idx];
          if (*cat_name_iter == 'C') {
            uint32_t char_code = *(++cat_name_iter);
            if ((char_code - 49) < 9) {
              uint32_t uii;
              if (!ScanPosintCapped(cat_name_iter, 0x7fffffff, &uii)) {
                continue;
              }
            }
          }
          break;
        }
        if (cat_idx == nn_cat_ct + 1) {
          SetBit(pheno_idx, is_basic_categorical);
        }
      }
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto ExportOxSample_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);

    write_iter = strcpya_k(write_iter, "0 0 0 D");
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      *write_iter++ = ' ';
      const PhenoDtype cur_type_code = pheno_cols[pheno_idx].type_code;
      if (cur_type_code == kPhenoDtypeCc) {
        *write_iter++ = 'B';
      } else if (cur_type_code == kPhenoDtypeQt) {
        // .psam file does not distinguish between "continuous covariate" and
        // "continuous phenotype", that's lost on round-trip
        *write_iter++ = 'P';
      } else {
        *write_iter++ = 'D';
      }
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto ExportOxSample_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);

    const double nonmale_geno_ct_recip = 1.0 / u31tod(variant_ct - y_ct);
    const double male_geno_ct_recip = 1.0 / u31tod(variant_ct);
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      const char* cur_sample_id = &(sample_ids[max_sample_id_blen * sample_uidx]);
      const char* fid_end = AdvToDelim(cur_sample_id, '\t');
      write_iter = memcpyax(write_iter, cur_sample_id, fid_end - cur_sample_id, ' ');
      write_iter = strcpya(write_iter, &(fid_end[1]));
      *write_iter++ = ' ';
      const int32_t cur_missing_geno_ct = sample_missing_geno_cts[sample_idx];
      if (IsSet(sex_male, sample_uidx)) {
        write_iter = dtoa_g(cur_missing_geno_ct * male_geno_ct_recip, write_iter);
        write_iter = strcpya_k(write_iter, " 1");
      } else {
        write_iter = dtoa_g(cur_missing_geno_ct * nonmale_geno_ct_recip, write_iter);
        *write_iter++ = ' ';
        if (IsSet(sex_nm, sample_uidx)) {
          *write_iter++ = '2';
        } else {
          write_iter = strcpya_k(write_iter, "NA");
        }
      }
      for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
        *write_iter++ = ' ';
        const PhenoCol* cur_pheno_col = &(pheno_cols[pheno_idx]);
        if (!IsSet(cur_pheno_col->nonmiss, sample_uidx)) {
          write_iter = strcpya_k(write_iter, "NA");
        } else {
          const PhenoDtype cur_type_code = cur_pheno_col->type_code;
          if (cur_type_code == kPhenoDtypeCc) {
            *write_iter++ = '0' + IsSet(cur_pheno_col->data.cc, sample_uidx);
          } else if (cur_type_code == kPhenoDtypeQt) {
            write_iter = dtoa_g(cur_pheno_col->data.qt[sample_uidx], write_iter);
          } else {
            const uint32_t cur_cat_idx = cur_pheno_col->data.cat[sample_uidx];
            if (IsSet(is_basic_categorical, pheno_idx)) {
              write_iter = strcpya(write_iter, &(cur_pheno_col->category_names[cur_cat_idx][1]));
            } else {
              write_iter = u32toa(cur_cat_idx, write_iter);
            }
          }
        }
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto ExportOxSample_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto ExportOxSample_ret_WRITE_FAIL;
    }
  }
  while (0) {
  ExportOxSample_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportOxSample_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportOxSample_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

uint32_t ValidVcfAlleleCode(const char* allele_code_iter) {
  // returns 1 if probably valid (angle-bracket case is not exhaustively
  // checked), 0 if definitely not
  uint32_t uii = ctou32(*allele_code_iter);
  if ((uii == '<') || ((uii == '*') && (!allele_code_iter[1]))) {
    return 1;
  }
  do {
    uii -= 64;
    // A = 1, C = 3, G = 7, N = 14, T = 20, so (0x10408a >> ucc) & 1 works as a
    // set membership test
#ifdef __LP64__
    if ((uii > 63) || (!((0x10408a0010408aLLU >> uii) & 1))) {
      // if '[', ']', or '.', assume breakend
      return ((uii == 27) || (uii == 29) || (uii == 0xffffffeeU))? 1 : 0;
    }
#else
    if ((uii > 63) || (!((0x10408a >> (uii % 32)) & 1))) {
      return ((uii == 27) || (uii == 29) || (uii == 0xffffffeeU))? 1 : 0;
    }
#endif
    uii = ctou32(*(++allele_code_iter));
  } while (uii);
  return 1;
}

char* PrintDiploidVcfDosage(uint32_t dosage_int, uint32_t write_ds, char* write_iter) {
  if (write_ds) {
    return PrintSmallDosage(dosage_int, write_iter);
  }
  // GP
  if (dosage_int <= kDosageMid) {
    write_iter = PrintSmallDosage(kDosageMid - dosage_int, write_iter);
    *write_iter++ = ',';
    write_iter = PrintSmallDosage(dosage_int, write_iter);
    return strcpya_k(write_iter, ",0");
  }
  write_iter = strcpya_k(write_iter, "0,");
  write_iter = PrintSmallDosage(kDosageMax - dosage_int, write_iter);
  *write_iter++ = ',';
  return PrintSmallDosage(dosage_int - kDosageMid, write_iter);
}

static_assert(kDosageMax == 32768, "PrintHdsPair() needs to be updated.");
char* PrintHdsPair(uint32_t dosage_int, int32_t dphase_delta, char* start) {
  // We print left and right dosages in [0, 1] separated by a comma, with the
  // minimum number of decimal places such that the sum
  // ([0, 32768] <-> [0.0, 2.0]) and difference
  // ([-16384, 16384] <-> [-1.0, 1.0]) round-trip.  Since it's sum/difference
  // that's relevant, it makes sense to force the same number of digits on both
  // sides (though trailing zeroes are still omitted), rather than perform a
  // multiparameter search (especially since this needs to be extended to
  // multiallelic variants soon).
  //
  // The 4-decimal-place representation is constructed as follows, when it
  // exists:
  // 1. Check if multiples of 10^{-4} are in
  //    ((dosage_int - 0.5)/16384, (dosage_int + 0.5)/16384)
  //    and
  //    ((dphase_delta - 0.5)/16384, (dphase_delta + 0.5)/16384)
  //    If not, we can give up.
  // 2. Check if these have the same parity.  If they do, (x+y)/2 and (x-y)/2
  //    get the job done; otherwise we abort.
  // When this doesn't work, we render the left and right sides separately to
  // 5 decimal places.  The difference and the sum can each be off by up to
  // 1/100000.
  // More generally, with k nonzero-dosage alleles and n decimal places, the
  // *omitted* allele difference and sum can each be off by up to (k-1)/(10^n),
  // and once that's >= 1/32768, we lose reliable round-tripping.  So we will
  // go up to a 6-decimal-place fallback at k=5, 7 decimal places at k=32, etc.
  //
  // For diploid unphased dosages, we want (k-1)/(2 * 10^n) to be less than
  // 1/32768, so we go to six digits at k=8, seven digits at k=63, etc.
  // For haploid dosages, we want (k-1)/(10^n) to be less than 1/65536, so we
  // go to six digits at k=3, seven digits at k=17, etc.

  // (dosage_int * 2) and (dphase_delta * 2) are in 32768ths
  // 32768 * 625 = 20480k, smallest common denominator with 10^4
  const uint32_t sum_rangetop_20480k = (dosage_int * 2 + 1) * 625;

  // +32768 to force this to be nonnegative while keeping decimal part the same
  const uint32_t diffp1_rangetop_20480k = S_CAST(uint32_t, dphase_delta * 2 + 32769) * 625;
  if (((sum_rangetop_20480k % 2048) < 1250) && ((diffp1_rangetop_20480k % 20480) < 1250) && ((sum_rangetop_20480k & 2048) == (diffp1_rangetop_20480k & 2048))) {
    const uint32_t sum_x10k = sum_rangetop_20480k / 2048;
    const int32_t diff_x10k = (diffp1_rangetop_20480k / 2048) - 10000;
    const uint32_t left_x10k = (sum_x10k + diff_x10k) / 2;
    const uint32_t right_x10k = (sum_x10k - diff_x10k) / 2;
    if (left_x10k == 10000) {
      *start++ = '1';
    } else {
      *start++ = '0';
      if (left_x10k) {
        *start++ = '.';
        start = u32toa_trunc4(left_x10k, start);
      }
    }
    *start++ = ',';
    if (right_x10k == 10000) {
      *start++ = '1';
    } else {
      *start++ = '0';
      if (right_x10k) {
        *start++ = '.';
        start = u32toa_trunc4(right_x10k, start);
      }
    }
    return start;
  }
  // see PrintHaploidNonintDosage() banker's rounding
  uint32_t val_x32768 = dosage_int + dphase_delta;
  // bugfix (11 Dec 2018): Forgot to print leading digit and dot here.
  if (!(val_x32768 & 32767)) {
    *start++ = '0' + (val_x32768 == 32768);
  } else {
    const uint32_t five_decimal_places = ((3125 * val_x32768 + 512) / 1024) - ((val_x32768 % 2048) == 512);
    const uint32_t first_decimal_place = five_decimal_places / 10000;
    start = strcpya_k(start, "0.");
    *start++ = '0' + first_decimal_place;
    const uint32_t last_four_digits = five_decimal_places - first_decimal_place * 10000;
    if (last_four_digits) {
      start = u32toa_trunc4(last_four_digits, start);
    }
  }
  *start++ = ',';
  val_x32768 = dosage_int - dphase_delta;
  if (!(val_x32768 & 32767)) {
    *start++ = '0' + (val_x32768 == 32768);
    return start;
  }
  const uint32_t five_decimal_places = ((3125 * val_x32768 + 512) / 1024) - ((val_x32768 % 2048) == 512);
  const uint32_t first_decimal_place = five_decimal_places / 10000;
  start = strcpya_k(start, "0.");
  *start++ = '0' + first_decimal_place;
  const uint32_t last_four_digits = five_decimal_places - first_decimal_place * 10000;
  if (last_four_digits) {
    return u32toa_trunc4(last_four_digits, start);
  }
  return start;
}


char* AppendVcfMultiallelicDsForce01(uint32_t allele_ct_m2, uint32_t hds_force, AlleleCode ac, uint32_t is_hethap, char* write_iter) {
  write_iter = strcpya_k(write_iter, "\t0/");
  write_iter = u32toa(ac, write_iter);
  // DS
  const uint32_t trailing_ac = allele_ct_m2 + 1 - ac;
  write_iter = strcpya_k(write_iter, ":0");
  write_iter = u16setsa(write_iter, 0x302c, ac - 2);
  write_iter = strcpya_k(write_iter, ",1");
  write_iter = u16setsa(write_iter, 0x302c, trailing_ac);
  if (hds_force) {
    write_iter = strcpya_k(write_iter, ":0");
    write_iter = u16setsa(write_iter, 0x302c, ac - 2);
    write_iter = strcpya_k(write_iter, ",0.5");
    if (!is_hethap) {
      // trailing_ac after, 1 + (ac - 2) before
      write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
      write_iter = strcpya_k(write_iter, ",0.5");
    }
    write_iter = u16setsa(write_iter, 0x302c, trailing_ac);
  }
  return write_iter;
}

char* AppendVcfMultiallelicDsForce10Het(uint32_t allele_ct_m2, uint32_t hds_force, AlleleCode ac0, AlleleCode ac1, uint32_t is_hethap, char* write_iter) {
  *write_iter++ = '\t';
  write_iter = u32toa_x(ac0, '/', write_iter);
  write_iter = u32toa_x(ac1, ':', write_iter);
  const uint32_t leading_ac = ac0 - 1;
  const uint32_t middle_ac = ac1 - ac0 - 1;
  const uint32_t trailing_ac = allele_ct_m2 + 1 - ac1;
  write_iter = u16setsa(write_iter, 0x2c30, leading_ac);
  *write_iter++ = '1';
  write_iter = u16setsa(write_iter, 0x302c, middle_ac);
  write_iter = strcpya_k(write_iter, ",1");
  write_iter = u16setsa(write_iter, 0x302c, trailing_ac);
  if (hds_force) {
    *write_iter++ = ':';
    write_iter = u16setsa(write_iter, 0x2c30, leading_ac);
    write_iter = strcpya_k(write_iter, "0.5");
    write_iter = u16setsa(write_iter, 0x302c, middle_ac);
    write_iter = strcpya_k(write_iter, ",0.5");
    if (!is_hethap) {
      write_iter = u16setsa(write_iter, 0x302c, trailing_ac + leading_ac);
      write_iter = strcpya_k(write_iter, ",0.5");
      write_iter = u16setsa(write_iter, 0x302c, middle_ac);
      write_iter = strcpya_k(write_iter, ",0.5");
    }
    write_iter = u16setsa(write_iter, 0x302c, trailing_ac);
  }
  return write_iter;
}

char* AppendVcfMultiallelicDsForce10HomDiploid(uint32_t allele_ct_m2, uint32_t hds_force, AlleleCode ac, char phase_char, char* write_iter) {
  *write_iter++ = '\t';
  write_iter = u32toa_x(ac, phase_char, write_iter);
  write_iter = u32toa(ac, write_iter);
  // ac >= 2 guaranteed
  const uint32_t trailing_ac = allele_ct_m2 + 1 - ac;
  write_iter = strcpya_k(write_iter, ":0");
  write_iter = u16setsa(write_iter, 0x302c, ac - 2);
  write_iter = strcpya_k(write_iter, ",2");
  write_iter = u16setsa(write_iter, 0x302c, trailing_ac);
  if (hds_force) {
    write_iter = strcpya_k(write_iter, ":0");
    write_iter = u16setsa(write_iter, 0x302c, ac - 2);
    write_iter = strcpya_k(write_iter, ",1");
    write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
    write_iter = strcpya_k(write_iter, ",1");
    write_iter = u16setsa(write_iter, 0x302c, trailing_ac);
  }
  return write_iter;
}

char* AppendVcfMultiallelicDsForce10Haploid(uint32_t allele_ct_m2, uint32_t hds_force, AlleleCode ac, char* write_iter) {
  *write_iter++ = '\t';
  write_iter = u32toa(ac, write_iter);
  // ac >= 2 guaranteed
  const uint32_t trailing_ac = allele_ct_m2 + 1 - ac;
  write_iter = strcpya_k(write_iter, ":0");
  write_iter = u16setsa(write_iter, 0x302c, ac - 2);
  write_iter = strcpya_k(write_iter, ",1");
  write_iter = u16setsa(write_iter, 0x302c, trailing_ac);
  if (hds_force) {
    write_iter = strcpya_k(write_iter, ":0");
    write_iter = u16setsa(write_iter, 0x302c, ac - 2);
    write_iter = strcpya_k(write_iter, ",1");
    write_iter = u16setsa(write_iter, 0x302c, trailing_ac);
  }
  return write_iter;
}

char* AppendVcfMultiallelicDsForcePhased(uint32_t allele_ct_m2, uint32_t hds_force, AlleleCode ac0, AlleleCode ac1, uint32_t phaseinfo_hw_masked, char* write_iter) {
  AlleleCode ac_left = ac0;
  AlleleCode ac_right = ac1;
  if (phaseinfo_hw_masked) {
    ac_left = ac1;
    ac_right = ac0;
  }
  *write_iter++ = '\t';
  write_iter = u32toa_x(ac_left, '|', write_iter);
  write_iter = u32toa_x(ac_right, ':', write_iter);
  // DS
  // Don't define leading_ac since ac0 == 0 is possible.
  const uint32_t middle_ac = ac1 - ac0 - 1;
  const uint32_t trailing_ac = allele_ct_m2 + 1 - ac1;
  if (ac0) {
    write_iter = u16setsa(write_iter, 0x2c30, ac0 - 1);
    write_iter = strcpya_k(write_iter, "1,");
  }
  write_iter = u16setsa(write_iter, 0x2c30, middle_ac);
  *write_iter++ = '1';
  write_iter = u16setsa(write_iter, 0x302c, trailing_ac);
  if (hds_force) {
    *write_iter++ = ':';
    if (ac_left) {
      write_iter = u16setsa(write_iter, 0x2c30, ac_left - 1);
      write_iter = strcpya_k(write_iter, "1,");
    }
    write_iter = u16setsa(write_iter, 0x2c30, allele_ct_m2 + 1 - ac_left);
    if (ac_right) {
      write_iter = u16setsa(write_iter, 0x2c30, ac_right - 1);
      write_iter = strcpya_k(write_iter, "1,");
    }
    write_iter = u16setsa(write_iter, 0x2c30, allele_ct_m2 + 1 - ac_right);
    // remove trailing comma
    --write_iter;
  }
  return write_iter;
}

#ifdef __arm__
#  error "Unaligned accesses in ExportVcf()."
#endif
PglErr ExportVcf(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const SampleIdInfo* siip, const uintptr_t* sex_male_collapsed, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, ExportfFlags exportf_flags, VcfExportMode vcf_mode, IdpasteFlags exportf_id_paste, char exportf_id_delim, char* xheader, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  TextStream pvar_reload_txs;
  BgzfCompressStream bgzf;
  PreinitTextStream(&pvar_reload_txs);
  PreinitBgzfCompressStream(&bgzf);
  {
    {
      uint32_t clvl = 0;
      if (!(exportf_flags & kfExportfBgz)) {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".vcf");
      } else {
        snprintf(outname_end, kMaxOutfnameExtBlen, ".vcf.gz");
        clvl = kBgzfDefaultClvl;
      }
      reterr = InitBgzfCompressStreamEx(outname, 0, clvl, max_thread_ct, &bgzf);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto ExportVcf_ret_1;
      }
    }
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    // CHROM, POS, ID, REF, one ALT, eoln
    uintptr_t writebuf_blen = kMaxIdSlen + 32 + max_chr_blen + 2 * max_allele_slen;
    // QUAL, FILTER, INFO, FORMAT, genotypes, eoln
    // needs to be larger for >9 alt alleles
    uint32_t write_some_dosage = 0;
    uint32_t write_ds = 0;
    uint32_t write_hds = 0;
    uint32_t ds_force = 0;
    uint32_t hds_force = 0;
    if (vcf_mode != kVcfExport0) {
      write_some_dosage = 1;
      if (vcf_mode != kVcfExportGp) {
        write_ds = 1;
        if (vcf_mode == kVcfExportDsForce) {
          ds_force = 1;
        } else if (vcf_mode != kVcfExportDs) {
          write_hds = 1;
          if (vcf_mode == kVcfExportHdsForce) {
            ds_force = 1;
            hds_force = 1;
          }
        }
      }
    }
    if ((!ds_force) && write_some_dosage && (!(pgfip->gflags & kfPgenGlobalDosagePresent))) {
      write_some_dosage = 0;
      logerrprintf("Warning: No dosage data present.  %s will not be exported.\n", write_hds? "DS and HDS fields" : (write_ds? "DS field" : "GP field"));
      write_ds = 0;
      write_hds = 0;
    }
    // max_allele_ct == 2:
    //   GT: 4 bytes
    //   GP: 3 limited-precision numbers, up to (7 chars + delim) * 3
    //   DS + HDS: also 3 limited precision numbers
    //   DS only: 1 limited-precision number
    // max_allele_ct > 2:
    //   GT: 2 + 2 * UintSlen(max_allele_ct + 1)
    //   GP: if requested, check max_gp_allele_ct among multiallelic-dosage
    //         variants.
    //       if none found (always true for now), 24 chars.
    //       if present, (max_gp_allele_ct * (max_gp_allele_ct + 1) / 2) * 8,
    //         which sucks ass
    //   DS + HDS: 3 * (max_allele_ct - 1) limited-precision numbers
    //   DS only: (max_allele_ct - 1) limited-precision numbers
    uintptr_t output_bytes_per_sample;
    if (!allele_idx_offsets) {
      if (write_some_dosage) {
        output_bytes_per_sample = (write_ds && (!write_hds))? 12 : 28;
      } else {
        output_bytes_per_sample = 4;
      }
    } else {
      const uint32_t max_allele_ct = PgrGetMaxAlleleCt(simple_pgrp);
      output_bytes_per_sample = 2 + 2 * UintSlen(max_allele_ct + 1);
      if (write_some_dosage) {
        if (write_ds) {
          if (write_hds) {
            output_bytes_per_sample += 24 * (max_allele_ct - 1);
          } else {
            output_bytes_per_sample += 8 * (max_allele_ct - 1);
          }
        } else {
          // GP should only be written for biallelic variants, for now
          // would need to take max of this and GT length, but this is always
          // larger
          output_bytes_per_sample = 28;
        }
      }
    }
    const uintptr_t writebuf_blen_lbound = sample_ct * output_bytes_per_sample + 32 + max_filter_slen + info_reload_slen;
    if (writebuf_blen < writebuf_blen_lbound) {
      writebuf_blen = writebuf_blen_lbound;
    }
    writebuf_blen += kMaxMediumLine;
    char* writebuf;
    if (unlikely(bigstack_alloc_c(writebuf_blen, &writebuf))) {
      goto ExportVcf_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = strcpya_k(writebuf, "##fileformat=VCFv4.");
    *write_iter++ = (exportf_flags & kfExportfVcf43)? '3' : '2';
    write_iter = strcpya_k(write_iter, EOLN_STR "##fileDate=");
    time_t rawtime;
    time(&rawtime);
    struct tm* loctime;
    loctime = localtime(&rawtime);
    write_iter += strftime(write_iter, kMaxMediumLine, "%Y%m%d", loctime);
    write_iter = strcpya_k(write_iter, EOLN_STR "##source=PLINKv2.00" EOLN_STR);
    if (cip->chrset_source) {
      AppendChrsetLine(cip, &write_iter);
    }
    if (unlikely(BgzfWrite(writebuf, write_iter - writebuf, &bgzf))) {
      goto ExportVcf_ret_WRITE_FAIL;
    }
    const uint32_t chr_ctl = BitCtToWordCt(cip->chr_ct);
    uintptr_t* written_contig_header_lines;
    if (unlikely(bigstack_calloc_w(chr_ctl, &written_contig_header_lines))) {
      goto ExportVcf_ret_NOMEM;
    }
    // bugfix (30 Aug 2018): remove extraneous PAR1/PAR2 ##contig header lines,
    // count them as part of chrX
    const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
    const uint32_t par1_code = cip->xymt_codes[kChrOffsetPAR1];
    const uint32_t par2_code = cip->xymt_codes[kChrOffsetPAR2];
    uint32_t x_contig_line_written = 0;
    if (xheader) {
      memcpyao_k(writebuf, "##contig=<ID=", 13);
      char* xheader_iter = xheader;
      char* xheader_end = &(xheader[xheader_blen]);
      char* line_end = xheader;
      while (line_end != xheader_end) {
        xheader_iter = line_end;
        line_end = AdvPastDelim(xheader_iter, '\n');
        const uint32_t slen = line_end - xheader_iter;
        if ((slen > 14) && StrStartsWithUnsafe(xheader_iter, "##contig=<ID=")) {
          char* contig_name_start = &(xheader_iter[13]);
          char* contig_name_end = S_CAST(char*, memchr(contig_name_start, ',', slen - 14));
          if (!contig_name_end) {
            // if this line is technically well-formed (ends in '>'), it's
            // useless anyway, throw it out
            continue;
          }
          // if GetChrCodeCounted() is modified to not mutate
          // contig_name_start[], xheader can be changed to const char*
          const uint32_t chr_idx = GetChrCodeCounted(cip, contig_name_end - contig_name_start, contig_name_start);
          // bugfix (8 Sep 2018): must exclude ##contig lines not present in
          // input, otherwise chr_fo_idx == 0xffffffffU, etc.
          if (IsI32Neg(chr_idx) || (!IsSet(cip->chr_mask, chr_idx)) || (chr_idx == par1_code) || (chr_idx == par2_code)) {
            continue;
          }
          if (chr_idx == x_code) {
            x_contig_line_written = 1;
          }
          const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
          if (unlikely(IsSet(written_contig_header_lines, chr_fo_idx))) {
            logerrputs("Error: Duplicate ##contig line in .pvar file.\n");
            goto ExportVcf_ret_MALFORMED_INPUT;
          }
          SetBit(chr_fo_idx, written_contig_header_lines);
          // if --output-chr was used at some point, we need to sync the
          // ##contig chromosome code with the code in the VCF body.
          write_iter = chrtoa(cip, chr_idx, &(writebuf[13]));
          if (unlikely(BgzfWrite(writebuf, write_iter - writebuf, &bgzf))) {
            goto ExportVcf_ret_WRITE_FAIL;
          }
          if (unlikely(BgzfWrite(contig_name_end, line_end - contig_name_end, &bgzf))) {
            goto ExportVcf_ret_WRITE_FAIL;
          }
        } else {
          if (unlikely(BgzfWrite(xheader_iter, slen, &bgzf))) {
            goto ExportVcf_ret_WRITE_FAIL;
          }
        }
      }
    }
    write_iter = writebuf;
    // fill in the missing ##contig lines
    uint32_t contig_zero_written = 0;
    uint32_t chrx_end = 0;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != cip->chr_ct; ++chr_fo_idx) {
      if (IsSet(written_contig_header_lines, chr_fo_idx)) {
        continue;
      }
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if ((!IsSet(cip->chr_mask, chr_idx)) || AllBitsAreZero(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1])) {
        continue;
      }
      if ((chr_idx == x_code) || (chr_idx == par1_code) || (chr_idx == par2_code)) {
        // if variant included, should add max allele length instead of 1
        const uint32_t pos_end = variant_bps[cip->chr_fo_vidx_start[chr_fo_idx + 1] - 1] + 1;
        if (pos_end > chrx_end) {
          chrx_end = pos_end;
        }
        continue;
      }
      char* chr_name_write_start = strcpya_k(write_iter, "##contig=<ID=");
      char* chr_name_write_end = chrtoa(cip, chr_idx, chr_name_write_start);
      if ((*chr_name_write_start == '0') && (chr_name_write_end == &(chr_name_write_start[1]))) {
        // --allow-extra-chr 0 special case
        if (contig_zero_written) {
          continue;
        }
        contig_zero_written = 1;
        write_iter = strcpya_k(chr_name_write_end, ",length=2147483645");
      } else {
        if (unlikely(memchr(chr_name_write_start, ':', chr_name_write_end - chr_name_write_start))) {
          logerrputs("Error: VCF chromosome codes may not include the ':' character.\n");
          goto ExportVcf_ret_MALFORMED_INPUT;
        }
        write_iter = strcpya_k(chr_name_write_end, ",length=");
        if (1) {
          // er, if variant included, should look at lengths of last few allele
          // codes, and compute max(pos + <max allele len in variant>)...
          write_iter = u32toa(variant_bps[cip->chr_fo_vidx_start[chr_fo_idx + 1] - 1] + 1, write_iter);
        } else {
          // todo: scan for maximum pos in unsorted map case
        }
      }
      *write_iter++ = '>';
      AppendBinaryEoln(&write_iter);
      if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
        goto ExportVcf_ret_WRITE_FAIL;
      }
    }
    if (chrx_end && (!x_contig_line_written)) {
      char* chr_name_write_start = strcpya_k(write_iter, "##contig=<ID=");
      char* chr_name_write_end = chrtoa(cip, x_code, chr_name_write_start);
      write_iter = strcpya_k(chr_name_write_end, ",length=");
      if (1) {
        write_iter = u32toa(chrx_end, write_iter);
      } else {
      }
      *write_iter++ = '>';
      AppendBinaryEoln(&write_iter);
      if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
        goto ExportVcf_ret_WRITE_FAIL;
      }
    }
    BigstackReset(written_contig_header_lines);
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t all_nonref = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uint32_t write_pr = all_nonref;
    if (nonref_flags) {
      write_pr = !IntersectionIsEmpty(variant_include, nonref_flags, raw_variant_ctl);
    }
    const uint32_t info_pr_flag_present = (info_flags / kfInfoPrFlagPresent) & 1;
    if (write_pr) {
      if (unlikely(info_flags & kfInfoPrNonflagPresent)) {
        logputs("\n");
        logerrputs("Error: Conflicting INFO:PR definitions.  Either fix all REF alleles so that the\n'provisional reference' flag is no longer needed, or remove/rename the other\nuse of the INFO:PR key.\n");
        goto ExportVcf_ret_INCONSISTENT_INPUT;
      }
      if (!info_pr_flag_present) {
        write_iter = strcpya_k(write_iter, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">" EOLN_STR);
      }
    }
    if (write_ds) {
      write_iter = strcpya_k(write_iter, "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">" EOLN_STR);
      if (write_hds) {
        // bugfix (3 May 2018): 'Number=2' was inaccurate for haploid calls.

        // Note that HDS ploidy intentionally does NOT match GT ploidy in the
        // unphased het haploid case.
        write_iter = strcpya_k(write_iter, "##FORMAT=<ID=HDS,Number=.,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage \">" EOLN_STR);
      }
    } else if (write_some_dosage) {
      write_iter = strcpya_k(write_iter, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \">" EOLN_STR);
    }
    // possible todo: optionally export .psam information as
    // PEDIGREE/META/SAMPLE lines in header, and make --vcf be able to read it
    write_iter = strcpya_k(write_iter, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" EOLN_STR "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    char* exported_sample_ids;
    uint32_t* exported_id_htable;
    uintptr_t max_exported_sample_id_blen;
    if (unlikely(ExportIdpaste(sample_include, siip, "vcf", sample_ct, exportf_id_paste, exportf_id_delim, &max_exported_sample_id_blen, &exported_sample_ids, &exported_id_htable))) {
      goto ExportVcf_ret_NOMEM;
    }
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, &(exported_sample_ids[sample_idx * max_exported_sample_id_blen]));
      if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
        goto ExportVcf_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);
    BigstackReset(exported_sample_ids);

    logprintfww5("--export vcf%s to %s ... ", (exportf_flags & kfExportfBgz)? " bgz" : "", outname);
    fputs("0%", stdout);
    fflush(stdout);

    // includes trailing tab
    char* chr_buf;

    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    PgenVariant pgv;
    if (unlikely(
            bigstack_alloc_c(max_chr_blen, &chr_buf) ||
            // if we weren't using bigstack_alloc, this would need to be
            // sample_ctaw2
            bigstack_alloc_w(sample_ctl2, &pgv.genovec))) {
      goto ExportVcf_ret_NOMEM;
    }
    pgv.patch_01_set = nullptr;
    pgv.patch_01_vals = nullptr;
    pgv.patch_10_set = nullptr;
    pgv.patch_10_vals = nullptr;
    if (allele_idx_offsets) {
      if (unlikely(
              bigstack_alloc_w(sample_ctl, &(pgv.patch_01_set)) ||
              bigstack_alloc_ac(sample_ct, &(pgv.patch_01_vals)) ||
              bigstack_alloc_w(sample_ctl, &(pgv.patch_10_set)) ||
              bigstack_alloc_ac(sample_ct * 2, &(pgv.patch_10_vals)))) {
        goto ExportVcf_ret_NOMEM;
      }
    }
    // For now, if phased data is present, each homozygous call is represented
    // as phased iff the previous heterozygous call was phased.  (If no
    // previous heterozygous call exists, it's treated as phased.)  This does
    // the right thing when the entire genome is phased, and it induces about
    // as good a phase set approximation as you can get without explicitly
    // saving that info.  But that approximation is still pretty inaccurate; as
    // soon as we have any use for them, explicit phase set support should be
    // added to pgenlib.
    const uint32_t load_dphase = write_hds && (pgfip->gflags & kfPgenGlobalDosagePhasePresent);
    const uint32_t some_phased = (pgfip->gflags & kfPgenGlobalHardcallPhasePresent) || load_dphase;
    uintptr_t* prev_phased = nullptr;
    pgv.phasepresent = nullptr;
    pgv.phaseinfo = nullptr;
    if (some_phased) {
      if (unlikely(
              bigstack_alloc_w(sample_ctl, &prev_phased) ||
              bigstack_alloc_w(sample_ctl, &(pgv.phasepresent)) ||
              bigstack_alloc_w(sample_ctl, &(pgv.phaseinfo)))) {
        goto ExportVcf_ret_NOMEM;
      }
      SetAllBits(sample_ct, prev_phased);
    }

    pgv.dosage_present = nullptr;
    pgv.dosage_main = nullptr;
    pgv.dphase_present = nullptr;
    pgv.dphase_delta = nullptr;
    if (write_some_dosage && (pgfip->gflags & kfPgenGlobalDosagePresent)) {
      if (unlikely(
              bigstack_alloc_w(sample_ctl, &(pgv.dosage_present)) ||
              bigstack_alloc_dosage(sample_ct, &(pgv.dosage_main)))) {
        goto ExportVcf_ret_NOMEM;
      }
      if (load_dphase) {
        if (unlikely(
                bigstack_alloc_w(sample_ctl, &(pgv.dphase_present)) ||
                bigstack_alloc_dphase(sample_ct, &(pgv.dphase_delta)))) {
          goto ExportVcf_ret_NOMEM;
        }
      }
    }

    char* pvar_reload_line_iter = nullptr;
    uint32_t info_col_idx = 0;
    if (pvar_info_reload) {
      reterr = PvarInfoOpenAndReloadHeader(pvar_info_reload, 1 + (max_thread_ct > 1), &pvar_reload_txs, &pvar_reload_line_iter, &info_col_idx);
      if (unlikely(reterr)) {
        goto ExportVcf_ret_TSTREAM_REWIND_FAIL;
      }
    }

    // assumes little-endian
    // maybe want to pass references to these arrays later, but leave as
    // C-style for now
    // \t0/0  \t0/1  \t1/1  \t./.
    const uint32_t basic_genotext[4] = {0x302f3009, 0x312f3009, 0x312f3109, 0x2e2f2e09};

    // 4-genotype-at-a-time lookup in basic case
    uint32_t basic_genotext4[1024] ALIGNV16;
    basic_genotext4[0] = basic_genotext[0];
    basic_genotext4[4] = basic_genotext[1];
    basic_genotext4[8] = basic_genotext[2];
    basic_genotext4[12] = basic_genotext[3];
    InitLookup256x4bx4(basic_genotext4);

    uint32_t haploid_genotext_blen[8];
    // lengths [0], [2], and [3] reset at beginning of chromosome
    // chrX: nonmales are diploid and use indices 0..3, males are 4..7
    // other: all haploid, indices 0..3
    haploid_genotext_blen[1] = 4;
    haploid_genotext_blen[4] = 2;
    haploid_genotext_blen[5] = 4;
    haploid_genotext_blen[6] = 2;
    haploid_genotext_blen[7] = 2;
    // don't bother exporting GP for hardcalls
    // usually don't bother for DS, but DS-force is an exception

    // 4..7 = haploid, 5 should never be looked up
    // :0 :1 :2 :. :0 !! :1 :.
    const uint16_t ds_inttext[8] = {0x303a, 0x313a, 0x323a, 0x2e3a, 0x303a, 0x2121, 0x313a, 0x2e3a};

    // 0..3 = diploid unphased
    // 4..5 = phased, phaseinfo=0 first
    // update (1 Aug 2018): HDS ploidy should be more trustworthy than GT
    // ploidy when there's a conflict, since GT must render unphased het
    // haploids as 0/1.  So under HDS-force, diploid missing value is '.,.'
    // instead of '.'.  (Since we're only interested in storing a trustworthy
    // ploidy, we don't add more missing entries in the multiallelic case.)

    // could make this include DS text in front
    // could expand this to e.g. 10 cases and make [4], [5] less fiddly
    // :0,0  :0.5,0.5  :1,1  :.,.  :0,1  :1,0
    const uint64_t hds_inttext[6] = {0x302c303a, 0x352e302c352e303aLLU, 0x312c313a, 0x2e2c2e3a, 0x312c303a, 0x302c313a};

    // [4]..[7] = haploid unphased; == 4 for phased
    uint32_t hds_inttext_blen[8];

    // lengths [0]-[3] reset at beginning of chromosome
    hds_inttext_blen[4] = 2;
    hds_inttext_blen[5] = 4;
    hds_inttext_blen[6] = 2;
    hds_inttext_blen[7] = 2;
    const char* dot_ptr = &(g_one_char_strs[92]);
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t is_x = 0;
    uint32_t is_haploid = 0;  // includes chrX and chrY
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t rls_variant_uidx = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uint32_t allele_ct = 2;
    uint32_t invalid_allele_code_seen = 0;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      // a lot of this is redundant with write_pvar(), may want to factor the
      // commonalities out
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        is_x = (chr_idx == cip->xymt_codes[kChrOffsetX]);
        is_haploid = IsSet(cip->haploid_mask, chr_idx);
        // forced --merge-par, with diploid male output (is_x NOT set, but
        // chromosome code is X/chrX)
        if ((chr_idx == cip->xymt_codes[kChrOffsetPAR1]) || (chr_idx == cip->xymt_codes[kChrOffsetPAR2])) {
          chr_idx = cip->xymt_codes[kChrOffsetX];
        }
        char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
        *chr_name_end = '\t';
        chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
        // bugfix (3 May 2018): forgot to update hds_inttext_blen[]
        hds_inttext_blen[0] = 4;
        hds_inttext_blen[1] = 8;
        hds_inttext_blen[2] = 4;
        hds_inttext_blen[3] = 4;
        if (is_haploid) {
          if (is_x) {
            haploid_genotext_blen[0] = 4;
            haploid_genotext_blen[2] = 4;
            haploid_genotext_blen[3] = 4;
          } else {
            haploid_genotext_blen[0] = 2;
            haploid_genotext_blen[2] = 2;
            haploid_genotext_blen[3] = 2;
            hds_inttext_blen[0] = 2;
            hds_inttext_blen[1] = 4;
            hds_inttext_blen[2] = 2;
            hds_inttext_blen[3] = 2;
          }
        }
      }
      // #CHROM
      write_iter = memcpya(write_iter, chr_buf, chr_buf_blen);

      // POS
      write_iter = u32toa_x(variant_bps[variant_uidx], '\t', write_iter);

      // ID
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], '\t');

      // REF, ALT
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx][0];
        alt1_allele_idx = refalt1_select[variant_uidx][1];
      }
      if (cur_alleles[ref_allele_idx] != dot_ptr) {
        write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
        if (!invalid_allele_code_seen) {
          invalid_allele_code_seen = !ValidVcfAlleleCode(cur_alleles[ref_allele_idx]);
        }
      } else {
        *write_iter++ = 'N';
      }
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, cur_alleles[alt1_allele_idx]);
      if (!invalid_allele_code_seen) {
        invalid_allele_code_seen = !ValidVcfAlleleCode(cur_alleles[alt1_allele_idx]);
      }
      if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
        goto ExportVcf_ret_WRITE_FAIL;
      }
      if (allele_ct > 2) {
        for (uint32_t cur_allele_uidx = 0; cur_allele_uidx != allele_ct; ++cur_allele_uidx) {
          if ((cur_allele_uidx == ref_allele_idx) || (cur_allele_uidx == alt1_allele_idx)) {
            // if this is noticeably suboptimal, have two loops, with inner
            // loop going up to cur_allele_stop.
            // (also wrap this in a function, this comes up a bunch of times)
            continue;
          }
          *write_iter++ = ',';
          write_iter = strcpya(write_iter, cur_alleles[cur_allele_uidx]);
          if (!invalid_allele_code_seen) {
            invalid_allele_code_seen = !ValidVcfAlleleCode(cur_alleles[cur_allele_uidx]);
          }
          if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
            goto ExportVcf_ret_WRITE_FAIL;
          }
        }
      }

      // QUAL
      *write_iter++ = '\t';
      if ((!pvar_qual_present) || (!IsSet(pvar_qual_present, variant_uidx))) {
        *write_iter++ = '.';
      } else {
        write_iter = ftoa_g(pvar_quals[variant_uidx], write_iter);
      }

      // FILTER
      *write_iter++ = '\t';
      if ((!pvar_filter_present) || (!IsSet(pvar_filter_present, variant_uidx))) {
        *write_iter++ = '.';
      } else if (!IsSet(pvar_filter_npass, variant_uidx)) {
        write_iter = strcpya_k(write_iter, "PASS");
      } else {
        write_iter = strcpya(write_iter, pvar_filter_storage[variant_uidx]);
      }

      // INFO
      *write_iter++ = '\t';
      const uint32_t is_pr = all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx));
      if (pvar_reload_line_iter) {
        reterr = PvarInfoReloadAndWrite(info_pr_flag_present, info_col_idx, variant_uidx, is_pr, &pvar_reload_txs, &pvar_reload_line_iter, &write_iter, &rls_variant_uidx);
        if (unlikely(reterr)) {
          goto ExportVcf_ret_TSTREAM_REWIND_FAIL;
        }
      } else {
        if (is_pr) {
          write_iter = strcpya_k(write_iter, "PR");
        } else {
          *write_iter++ = '.';
        }
      }

      // FORMAT
      write_iter = strcpya_k(write_iter, "\tGT");

      // could defensively zero out more counts
      pgv.dosage_ct = 0;
      pgv.dphase_ct = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      if (allele_ct == 2) {
        if (!some_phased) {
          // biallelic, nothing phased in entire file
          // (technically possible for dosage-phase to be present, if no
          // hardcalls are phased and HDS output not requested)
          if (!write_some_dosage) {
            reterr = PgrGet(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, pgv.genovec);
          } else {
            reterr = PgrGetD(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
          }
          if (unlikely(reterr)) {
            goto ExportVcf_ret_PGR_FAIL;
          }
          if (!alt1_allele_idx) {
            // assumes biallelic
            GenovecInvertUnsafe(sample_ct, pgv.genovec);
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
            }
          }
          if ((!pgv.dosage_ct) && (!ds_force)) {
            if (!is_haploid) {
              // always 4 bytes wide, exploit that
              GenoarrLookup256x4bx4(pgv.genovec, basic_genotext4, sample_ct, write_iter);
              write_iter = &(write_iter[sample_ct * 4]);
            } else {
              // chrX: male homozygous/missing calls use only one character +
              //       tab
              // other haploid/MT: this is true for nonmales too
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                uint32_t sex_male_hw = is_x * (R_CAST(const Halfword*, sex_male_collapsed)[widx]);
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  const uint32_t cur_is_male = sex_male_hw & 1;
                  memcpy(write_iter, &(basic_genotext[cur_geno]), 4);
                  write_iter = &(write_iter[haploid_genotext_blen[cur_geno + cur_is_male * 4]]);
                  genovec_word >>= 2;
                  sex_male_hw >>= 1;
                }
              }
            }
          } else {
            // some dosages present, or {H}DS-force; unphased
            if (write_ds) {
              write_iter = strcpya_k(write_iter, ":DS");
              if (hds_force) {
                write_iter = strcpya_k(write_iter, ":HDS");
              }
            } else {
              write_iter = strcpya_k(write_iter, ":GP");
            }
            Dosage* dosage_main_iter = pgv.dosage_main;
            uint32_t dosage_present_hw = 0;
            if (!is_haploid) {
              // autosomal diploid, unphased
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                if (pgv.dosage_ct) {
                  dosage_present_hw = R_CAST(Halfword*, pgv.dosage_present)[widx];
                }
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
                  if (dosage_present_hw & 1) {
                    *write_iter++ = ':';
                    const uint32_t dosage_int = *dosage_main_iter++;
                    write_iter = PrintDiploidVcfDosage(dosage_int, write_ds, write_iter);
                    if (hds_force) {
                      *write_iter++ = ':';
                      char* write_iter2 = PrintHaploidNonintDosage(dosage_int, write_iter);
                      write_iter2[0] = ',';
                      write_iter = memcpya(&(write_iter2[1]), write_iter, write_iter2 - write_iter);
                    }
                  } else if (ds_force) {
                    write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno]), 2);
                    if (hds_force) {
                      memcpy(write_iter, &(hds_inttext[cur_geno]), 8);
                      write_iter = &(write_iter[hds_inttext_blen[cur_geno]]);
                    }
                  }
                  genovec_word >>= 2;
                  dosage_present_hw >>= 1;
                }
              }
            } else {
              // at least partly haploid, unphased
              uint32_t sex_male_hw = 0;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                if (is_x) {
                  sex_male_hw = R_CAST(const Halfword*, sex_male_collapsed)[widx];
                }
                if (pgv.dosage_ct) {
                  dosage_present_hw = R_CAST(Halfword*, pgv.dosage_present)[widx];
                }
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  const uint32_t cur_is_male = sex_male_hw & 1;
                  const uint32_t cur_genotext_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                  memcpy(write_iter, &(basic_genotext[cur_geno]), 4);
                  write_iter = &(write_iter[cur_genotext_blen]);
                  if (dosage_present_hw & 1) {
                    *write_iter++ = ':';
                    uint32_t dosage_int = *dosage_main_iter++;
                    if (cur_genotext_blen == 2) {
                      // render current hardcall as haploid
                      if (write_ds) {
                        char* write_iter2 = PrintHaploidNonintDosage(dosage_int, write_iter);
                        if (hds_force) {
                          write_iter2[0] = ':';
                          write_iter2 = memcpya(&(write_iter2[1]), write_iter, write_iter2 - write_iter);
                        }
                        write_iter = write_iter2;
                      } else {
                        // GP
                        write_iter = PrintHaploidNonintDosage(kDosageMax - dosage_int, write_iter);
                        *write_iter++ = ',';
                        write_iter = PrintHaploidNonintDosage(dosage_int, write_iter);
                      }
                    } else {
                      // render current hardcall as diploid (female X, or het
                      // haploid)
                      write_iter = PrintDiploidVcfDosage(dosage_int, write_ds, write_iter);
                      if (hds_force) {
                        *write_iter++ = ':';
                        char* write_iter2 = PrintHaploidNonintDosage(dosage_int, write_iter);
                        if (is_x && (!cur_is_male)) {
                          // but do not render phased-dosage as diploid in het
                          // haploid case
                          write_iter2[0] = ',';
                          write_iter2 = memcpya(&(write_iter2[1]), write_iter, write_iter2 - write_iter);
                        }
                        write_iter = write_iter2;
                      }
                    }
                  } else if (ds_force) {
                    write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno + 8 - 2 * cur_genotext_blen]), 2);
                    if (hds_force) {
                      memcpy(write_iter, &(hds_inttext[cur_geno]), 8);
                      write_iter = &(write_iter[hds_inttext_blen[cur_is_male * 4 + cur_geno]]);
                    }
                  }
                  genovec_word >>= 2;
                  sex_male_hw >>= 1;
                  dosage_present_hw >>= 1;
                }
              }
            }
          }
        } else {
          // biallelic, phased
          if (!write_some_dosage) {
            reterr = PgrGetP(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &(pgv.phasepresent_ct));
          } else {
            reterr = PgrGetDp(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, &pgv);
          }
          if (unlikely(reterr)) {
            goto ExportVcf_ret_PGR_FAIL;
          }
          if (!alt1_allele_idx) {
            // assumes biallelic
            GenovecInvertUnsafe(sample_ct, pgv.genovec);
            if (pgv.phasepresent_ct) {
              BitvecInvert(sample_ctl, pgv.phaseinfo);
            }
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
              if (pgv.dphase_ct) {
                BiallelicDphase16Invert(pgv.dphase_ct, pgv.dphase_delta);
              }
            }
          }
          uint32_t phasepresent_hw = 0;
          if ((!pgv.dosage_ct) && (!ds_force)) {
            if (!is_haploid) {
              // can't just use PhaseLookup4b(), thanks to prev_phased_halfword
              uint32_t* write_iter_u32_alias = R_CAST(uint32_t*, write_iter);
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                uint32_t prev_phased_halfword = R_CAST(Halfword*, prev_phased)[widx];

                if (pgv.phasepresent_ct) {
                  phasepresent_hw = R_CAST(Halfword*, pgv.phasepresent)[widx];
                }
                const uint32_t phaseinfo_hw = R_CAST(Halfword*, pgv.phaseinfo)[widx];
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uintptr_t cur_geno = genovec_word & 3;

                  // usually "\t0/0", etc.
                  uint32_t cur_basic_genotext = basic_genotext[cur_geno];
                  if (cur_geno == 1) {
                    const uint32_t cur_shift = (1U << sample_idx_lowbits);
                    if (phasepresent_hw & cur_shift) {
                      prev_phased_halfword |= cur_shift;
                      if (phaseinfo_hw & cur_shift) {
                        cur_basic_genotext ^= 0x1000100;  // 0|1 -> 1|0
                      }
                    } else {
                      prev_phased_halfword &= ~cur_shift;
                    }
                  }
                  // '/' = ascii 47, '|' = ascii 124
                  *write_iter_u32_alias++ = cur_basic_genotext + 0x4d0000 * ((prev_phased_halfword >> sample_idx_lowbits) & 1);
                  genovec_word >>= 2;
                }
                R_CAST(Halfword*, prev_phased)[widx] = prev_phased_halfword;
              }
              write_iter = R_CAST(char*, write_iter_u32_alias);
            } else {
              uint32_t is_male_hw = 0;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                if (is_x) {
                  is_male_hw = R_CAST(const Halfword*, sex_male_collapsed)[widx];
                }
                uint32_t prev_phased_halfword = R_CAST(Halfword*, prev_phased)[widx];

                if (pgv.phasepresent_ct) {
                  phasepresent_hw = R_CAST(Halfword*, pgv.phasepresent)[widx];
                }

                const uint32_t phaseinfo_hw = R_CAST(Halfword*, pgv.phaseinfo)[widx];
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  const uint32_t cur_is_male = is_male_hw & 1;
                  const uint32_t cur_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                  memcpy(write_iter, &(basic_genotext[cur_geno]), 4);
                  write_iter = &(write_iter[cur_blen]);
                  if (cur_blen == 4) {
                    if (cur_geno == 1) {
                      // a bit redundant with how is_male_hw is handled, but
                      // updating this on every loop iteration doesn't seem
                      // better
                      const uint32_t cur_shift = (1U << sample_idx_lowbits);
                      if (phasepresent_hw & cur_shift) {
                        prev_phased_halfword |= cur_shift;
                        if (phaseinfo_hw & cur_shift) {
                          memcpy(&(write_iter[-4]), "\t1|0", 4);
                        } else {
                          write_iter[-2] = '|';
                        }
                      } else {
                        prev_phased_halfword &= ~cur_shift;
                      }
                    } else if ((prev_phased_halfword >> sample_idx_lowbits) & 1) {
                      write_iter[-2] = '|';
                    }
                  }
                  genovec_word >>= 2;
                  is_male_hw >>= 1;
                }
                R_CAST(Halfword*, prev_phased)[widx] = prev_phased_halfword;
              }
            }
          } else {
            // both dosage (or {H}DS-force) and phase present
            if (write_ds) {
              write_iter = strcpya_k(write_iter, ":DS");
              if (hds_force || pgv.dphase_ct ||
                  (write_hds && pgv.phasepresent_ct && pgv.dosage_ct && (!IntersectionIsEmpty(pgv.phasepresent, pgv.dosage_present, sample_ctl)))) {
                write_iter = strcpya_k(write_iter, ":HDS");
                // no need to clear dphase_present, since we zero-initialize
                // dphase_present_hw and never refresh it when dphase_ct == 0
              }
            } else {
              write_iter = strcpya_k(write_iter, ":GP");
            }
            Dosage* dosage_main_iter = pgv.dosage_main;
            SDosage* dphase_delta_iter = pgv.dphase_delta;
            uint32_t dosage_present_hw = 0;
            uint32_t dphase_present_hw = 0;
            if (!is_haploid) {
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                uint32_t prev_phased_halfword = R_CAST(Halfword*, prev_phased)[widx];

                if (pgv.phasepresent_ct) {
                  phasepresent_hw = R_CAST(Halfword*, pgv.phasepresent)[widx];
                }

                const uint32_t phaseinfo_hw = R_CAST(Halfword*, pgv.phaseinfo)[widx];
                if (pgv.dosage_ct) {
                  dosage_present_hw = R_CAST(Halfword*, pgv.dosage_present)[widx];
                  if (pgv.dphase_ct) {
                    dphase_present_hw = R_CAST(Halfword*, pgv.dphase_present)[widx];
                  }
                }
                uint32_t cur_shift = 1;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
                  if (cur_geno == 1) {
                    if (phasepresent_hw & cur_shift) {
                      prev_phased_halfword |= cur_shift;
                      if (phaseinfo_hw & cur_shift) {
                        memcpy(&(write_iter[-4]), "\t1|0", 4);
                      }
                    } else {
                      prev_phased_halfword &= ~cur_shift;
                    }
                  }
                  if (prev_phased_halfword & cur_shift) {
                    write_iter[-2] = '|';
                  }
                  if (dosage_present_hw & cur_shift) {
                    *write_iter++ = ':';
                    const uint32_t dosage_int = *dosage_main_iter++;
                    write_iter = PrintDiploidVcfDosage(dosage_int, write_ds, write_iter);
                    // bugfix (29 May 2018): don't print HDS field if not
                    // requested
                    if (write_hds) {
                      if ((phasepresent_hw | dphase_present_hw) & cur_shift) {
                        int32_t cur_dphase_delta;
                        if (dphase_present_hw & cur_shift) {
                          cur_dphase_delta = *dphase_delta_iter++;
                        } else {
                          cur_dphase_delta = DosageHomdist(dosage_int);
                          if (!(phaseinfo_hw & cur_shift)) {
                            cur_dphase_delta = -cur_dphase_delta;
                          }
                        }
                        *write_iter++ = ':';
                        write_iter = PrintHdsPair(dosage_int, cur_dphase_delta, write_iter);
                      } else if (hds_force) {
                        *write_iter++ = ':';
                        char* write_iter2 = PrintHaploidNonintDosage(dosage_int, write_iter);
                        write_iter2[0] = ',';
                        write_iter = memcpya(&(write_iter2[1]), write_iter, write_iter2 - write_iter);
                      }
                    }
                  } else if (ds_force) {
                    write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno]), 2);
                    if (hds_force) {
                      uint32_t hds_inttext_index = cur_geno;
                      uint32_t tmp_blen = hds_inttext_blen[cur_geno];
                      if (phasepresent_hw & cur_shift) {
                        // do we want to remove this branch?  doubt it's
                        // worthwhile since variable-length memcpy will branch
                        // anyway...
                        hds_inttext_index = 4 + ((phaseinfo_hw >> sample_idx_lowbits) & 1);
                        tmp_blen = 4;
                      }
                      memcpy(write_iter, &(hds_inttext[hds_inttext_index]), 8);
                      write_iter = &(write_iter[tmp_blen]);
                    }
                  }
                  genovec_word >>= 2;
                  cur_shift <<= 1;
                }
                R_CAST(Halfword*, prev_phased)[widx] = prev_phased_halfword;
              }
            } else {
              // dosage (or {H}DS-force) and phase present, partly/fully
              // haploid
              uint32_t is_male_hw = 0;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                if (is_x) {
                  is_male_hw = R_CAST(const Halfword*, sex_male_collapsed)[widx];
                }
                uint32_t prev_phased_halfword = R_CAST(Halfword*, prev_phased)[widx];

                if (pgv.phasepresent_ct) {
                  phasepresent_hw = R_CAST(Halfword*, pgv.phasepresent)[widx];
                }

                const uint32_t phaseinfo_hw = R_CAST(Halfword*, pgv.phaseinfo)[widx];
                if (pgv.dosage_ct) {
                  dosage_present_hw = R_CAST(Halfword*, pgv.dosage_present)[widx];
                  if (pgv.dphase_ct) {
                    dphase_present_hw = R_CAST(Halfword*, pgv.dphase_present)[widx];
                  }
                }
                uint32_t cur_shift = 1;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  const uint32_t cur_is_male = is_male_hw & 1;
                  const uint32_t cur_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                  memcpy(write_iter, &(basic_genotext[cur_geno]), 4);
                  write_iter = &(write_iter[cur_blen]);
                  if (cur_blen == 4) {
                    // render current hardcall as diploid (chrX nonmale, or het
                    // haploid)
                    if (cur_geno == 1) {
                      if (phasepresent_hw & cur_shift) {
                        prev_phased_halfword |= cur_shift;
                        if (phaseinfo_hw & cur_shift) {
                          memcpy(&(write_iter[-4]), "\t1|0", 4);
                        }
                      } else {
                        prev_phased_halfword &= ~cur_shift;
                      }
                    }
                    if (prev_phased_halfword & cur_shift) {
                      write_iter[-2] = '|';
                    }
                    if (dosage_present_hw & cur_shift) {
                      *write_iter++ = ':';
                      const uint32_t dosage_int = *dosage_main_iter++;
                      write_iter = PrintDiploidVcfDosage(dosage_int, write_ds, write_iter);
                      if (write_hds) {
                        if ((phasepresent_hw | dphase_present_hw) & cur_shift) {
                          int32_t cur_dphase_delta;
                          if (dphase_present_hw & cur_shift) {
                            cur_dphase_delta = *dphase_delta_iter++;
                          } else {
                            cur_dphase_delta = DosageHomdist(dosage_int);
                            if (!(phaseinfo_hw & cur_shift)) {
                              cur_dphase_delta = -cur_dphase_delta;
                            }
                          }
                          *write_iter++ = ':';
                          write_iter = PrintHdsPair(dosage_int, cur_dphase_delta, write_iter);
                        } else if (hds_force) {
                          *write_iter++ = ':';
                          char* write_iter2 = PrintHaploidNonintDosage(dosage_int, write_iter);
                          // do not render phased-dosage as diploid in unphased
                          // het haploid case
                          if (is_x && (!cur_is_male)) {
                            write_iter2[0] = ',';
                            write_iter2 = memcpya(&(write_iter2[1]), write_iter, write_iter2 - write_iter);
                          }
                          write_iter = write_iter2;
                        }
                      }
                    } else if (ds_force) {
                      write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno]), 2);
                      if (hds_force) {
                        uint32_t hds_inttext_index = cur_geno;
                        uint32_t tmp_blen;
                        if (phasepresent_hw & cur_shift) {
                          // do we want to remove this branch?  doubt it's
                          // worthwhile since variable-length memcpy will
                          // branch anyway...
                          hds_inttext_index = 4 + ((phaseinfo_hw >> sample_idx_lowbits) & 1);
                          tmp_blen = 4;
                        } else {
                          // do not render phased-dosage as diploid in het
                          // haploid case
                          tmp_blen = hds_inttext_blen[cur_is_male * 4 + cur_geno];
                        }
                        memcpy(write_iter, &(hds_inttext[hds_inttext_index]), 8);
                        write_iter = &(write_iter[tmp_blen]);
                      }
                    }
                  } else {
                    // render current hardcall as haploid
                    // (can't get here for hardcall-phased)
                    if (dosage_present_hw & cur_shift) {
                      *write_iter++ = ':';
                      const uint32_t dosage_int = *dosage_main_iter++;
                      if (write_ds) {
                        write_iter = PrintHaploidNonintDosage(dosage_int, write_iter);
                        if (dphase_present_hw & cur_shift) {
                          // explicit dosage-phase, so render HDS as diploid
                          const int32_t cur_dphase_delta = *dphase_delta_iter++;
                          *write_iter++ = ':';
                          write_iter = PrintHdsPair(dosage_int, cur_dphase_delta, write_iter);
                        } else if (hds_force) {
                          *write_iter++ = ':';
                          write_iter = PrintHaploidNonintDosage(dosage_int, write_iter);
                        }
                      } else {
                        // GP
                        write_iter = PrintHaploidNonintDosage(kDosageMax - dosage_int, write_iter);
                        *write_iter++ = ',';
                        write_iter = PrintHaploidNonintDosage(dosage_int, write_iter);
                      }
                    } else if (ds_force) {
                      write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno + 4]), 2);
                      if (hds_force) {
                        memcpy(write_iter, &(hds_inttext[cur_geno]), 8);
                        write_iter = &(write_iter[hds_inttext_blen[cur_geno + 4]]);
                      }
                    }
                  }
                  genovec_word >>= 2;
                  is_male_hw >>= 1;
                  cur_shift <<= 1;
                }
                R_CAST(Halfword*, prev_phased)[widx] = prev_phased_halfword;
              }
            }
          }
        }
      } else {
        // multiallelic cases
        // multiallelic dosage not supported yet
        if (ref_allele_idx || (alt1_allele_idx != 1)) {
          // todo: rotation function that can also be used by --make-pgen
          // maybe add a PgrGetM2() function too
          logputs("\n");
          logerrputs("Error: VCF-export multiallelic rotation is under development.\n");
          reterr = kPglRetNotYetSupported;
          goto ExportVcf_ret_1;
        }
        if (!some_phased) {
          reterr = PgrGetM(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, &pgv);
          if (unlikely(reterr)) {
            goto ExportVcf_ret_PGR_FAIL;
          }
          if (!ds_force) {
            if (!is_haploid) {
              if (allele_ct <= 10) {
                GenoarrLookup256x4bx4(pgv.genovec, basic_genotext4, sample_ct, write_iter);
                if (pgv.patch_01_ct) {
                  // Patch some 0/1 entries to 0/x.
                  char* genotext_offset3 = &(write_iter[3]);
                  uintptr_t sample_idx_base = 0;
                  uintptr_t patch_01_bits = pgv.patch_01_set[0];
                  for (uint32_t uii = 0; uii != pgv.patch_01_ct; ++uii) {
                    const uintptr_t sample_idx = BitIter1(pgv.patch_01_set, &sample_idx_base, &patch_01_bits);
                    genotext_offset3[4 * sample_idx] = '0' + pgv.patch_01_vals[uii];
                  }
                }
                if (pgv.patch_10_ct) {
                  // Patch some 1/1 entries to x/y.
                  char* genotext_offset1 = &(write_iter[1]);
                  uintptr_t sample_idx_base = 0;
                  uintptr_t patch_10_bits = pgv.patch_10_set[0];
                  for (uint32_t uii = 0; uii != pgv.patch_10_ct; ++uii) {
                    const uintptr_t sample_idx = BitIter1(pgv.patch_10_set, &sample_idx_base, &patch_10_bits);
                    genotext_offset1[4 * sample_idx] = '0' + pgv.patch_10_vals[2 * uii];
                    genotext_offset1[4 * sample_idx + 2] = '0' + pgv.patch_10_vals[2 * uii + 1];
                  }
                }
                write_iter = &(write_iter[sample_ct * 4]);
              } else {
                if (!pgv.patch_01_ct) {
                  ZeroWArr(sample_ctl, pgv.patch_01_set);
                }
                if (!pgv.patch_10_ct) {
                  ZeroWArr(sample_ctl, pgv.patch_10_set);
                }
                const AlleleCode* patch_01_vals_iter = pgv.patch_01_vals;
                const AlleleCode* patch_10_vals_iter = pgv.patch_10_vals;
                for (uint32_t widx = 0; ; ++widx) {
                  if (widx >= sample_ctl2_m1) {
                    if (widx > sample_ctl2_m1) {
                      break;
                    }
                    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                  }
                  uintptr_t genovec_word = pgv.genovec[widx];
                  uint32_t multiallelic_hw = (R_CAST(const Halfword*, pgv.patch_01_set)[widx]) | (R_CAST(const Halfword*, pgv.patch_10_set)[widx]);
                  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                    const uint32_t cur_geno = genovec_word & 3;
                    if (!(multiallelic_hw & 1)) {
                      write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
                    } else if (cur_geno == 1) {
                      write_iter = strcpya_k(write_iter, "\t0/");
                      const AlleleCode ac = *patch_01_vals_iter++;
                      write_iter = u32toa(ac, write_iter);
                    } else {
                      AlleleCode ac = *patch_10_vals_iter++;
                      *write_iter++ = '\t';
                      write_iter = u32toa_x(ac, '/', write_iter);
                      ac = *patch_10_vals_iter++;
                      write_iter = u32toa(ac, write_iter);
                    }
                    genovec_word >>= 2;
                    multiallelic_hw >>= 1;
                  }
                }
              }
            } else {
              if (!pgv.patch_01_ct) {
                ZeroWArr(sample_ctl, pgv.patch_01_set);
              }
              if (!pgv.patch_10_ct) {
                ZeroWArr(sample_ctl, pgv.patch_10_set);
              }
              // at least partially haploid, !ds_force
              // We don't separately the allele_ct <= 10 vs. > 10 cases when
              // entries are already variable-width in the <= 10 case.
              const AlleleCode* patch_01_vals_iter = pgv.patch_01_vals;
              const AlleleCode* patch_10_vals_iter = pgv.patch_10_vals;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                uint32_t sex_male_hw = is_x * (R_CAST(const Halfword*, sex_male_collapsed)[widx]);
                // no need to separate patch_01 and patch_10 since this
                // information is redundant with genovec_word contents
                uint32_t multiallelic_hw = (R_CAST(const Halfword*, pgv.patch_01_set)[widx]) | (R_CAST(const Halfword*, pgv.patch_10_set)[widx]);
                // probable todo: multiallelic_hw == 0 fast path, check
                // whether <= inner_loop_last vs. < inner_loop_end makes a
                // difference, etc.
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  const uint32_t cur_is_male = sex_male_hw & 1;
                  if (!(multiallelic_hw & 1)) {
                    memcpy(write_iter, &(basic_genotext[cur_geno]), 4);
                    write_iter = &(write_iter[haploid_genotext_blen[cur_geno + cur_is_male * 4]]);
                  } else if (cur_geno == 1) {
                    // always heterozygous, 4 characters
                    write_iter = strcpya_k(write_iter, "\t0/");
                    const AlleleCode ac = *patch_01_vals_iter++;
                    write_iter = u32toa(ac, write_iter);
                  } else {
                    const AlleleCode ac0 = *patch_10_vals_iter++;
                    const AlleleCode ac1 = *patch_10_vals_iter++;
                    *write_iter++ = '\t';
                    write_iter = u32toa(ac0, write_iter);
                    if ((ac0 != ac1) || (haploid_genotext_blen[2 + cur_is_male * 4] == 4)) {
                      *write_iter++ = '/';
                      write_iter = u32toa(ac1, write_iter);
                    }
                  }
                  genovec_word >>= 2;
                  sex_male_hw >>= 1;
                  multiallelic_hw >>= 1;
                }
              }
            }
          } else {
            if (!pgv.patch_01_ct) {
              ZeroWArr(sample_ctl, pgv.patch_01_set);
            }
            if (!pgv.patch_10_ct) {
              ZeroWArr(sample_ctl, pgv.patch_10_set);
            }
            // ds_force, !some_phased
            write_iter = strcpya_k(write_iter, ":DS");
            if (hds_force) {
              write_iter = strcpya_k(write_iter, ":HDS");
            }
            const uint32_t allele_ct_m2 = allele_ct - 2;
            const AlleleCode* patch_01_vals_iter = pgv.patch_01_vals;
            const AlleleCode* patch_10_vals_iter = pgv.patch_10_vals;
            if (!is_haploid) {
              // DS-force, autosomal diploid, unphased, dosage_ct == 0
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                uint32_t multiallelic_hw = (R_CAST(const Halfword*, pgv.patch_01_set)[widx]) | (R_CAST(const Halfword*, pgv.patch_10_set)[widx]);
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  if (!(multiallelic_hw & 1)) {
                    write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
                    // DS
                    write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno]), 2);
                    if (cur_geno != 3) {
                      // repeat ",0" if nonmissing
                      write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                      if (hds_force) {
                        const uint64_t cur_inttext = hds_inttext[cur_geno];
                        const uint32_t hap_blen = hds_inttext_blen[cur_geno + 4];
                        *R_CAST(uint32_t*, write_iter) = cur_inttext;
                        write_iter = &(write_iter[hap_blen]);
                        write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                        *R_CAST(uint32_t*, write_iter) = cur_inttext >> (8 * hap_blen);
                        write_iter = &(write_iter[hap_blen]);
                        write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                      }
                    } else {
                      if (hds_force) {
                        write_iter = strcpya_k(write_iter, ":.,.");
                      }
                    }
                  } else if (cur_geno == 1) {
                    const AlleleCode ac = *patch_01_vals_iter++;
                    write_iter = AppendVcfMultiallelicDsForce01(allele_ct_m2, hds_force, ac, 0, write_iter);
                  } else {
                    const AlleleCode ac0 = *patch_10_vals_iter++;
                    const AlleleCode ac1 = *patch_10_vals_iter++;
                    if (ac0 != ac1) {
                      write_iter = AppendVcfMultiallelicDsForce10Het(allele_ct_m2, hds_force, ac0, ac1, 0, write_iter);
                    } else {
                      write_iter = AppendVcfMultiallelicDsForce10HomDiploid(allele_ct_m2, hds_force, ac0, '/', write_iter);
                    }
                  }
                  genovec_word >>= 2;
                  multiallelic_hw >>= 1;
                }
              }
            } else {
              // DS-force, at least partly haploid, unphased
              uint32_t sex_male_hw = 0;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                if (is_x) {
                  sex_male_hw = R_CAST(const Halfword*, sex_male_collapsed)[widx];
                }
                uint32_t multiallelic_hw = (R_CAST(const Halfword*, pgv.patch_01_set)[widx]) | (R_CAST(const Halfword*, pgv.patch_10_set)[widx]);
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  const uint32_t cur_is_male = sex_male_hw & 1;
                  const uint32_t should_be_diploid = is_x && (!cur_is_male);
                  if (!(multiallelic_hw & 1)) {
                    const uint32_t cur_genotext_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                    memcpy(write_iter, &(basic_genotext[cur_geno]), 4);
                    write_iter = &(write_iter[cur_genotext_blen]);
                    // DS
                    write_iter = memcpya(write_iter, &(ds_inttext[cur_geno + 8 - 2 * cur_genotext_blen]), 2);
                    if (cur_geno != 3) {
                      // repeat ",0"
                      write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                      if (hds_force) {
                        const uint64_t cur_inttext = hds_inttext[cur_geno];
                        const uint32_t hap_blen = hds_inttext_blen[4 + cur_geno];
                        *R_CAST(uint32_t*, write_iter) = cur_inttext;
                        write_iter = &(write_iter[hap_blen]);
                        if (should_be_diploid) {
                          write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                          *R_CAST(uint32_t*, write_iter) = cur_inttext >> (8 * hap_blen);
                          write_iter = &(write_iter[hap_blen]);
                        }
                        write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                      }
                    } else {
                      if (hds_force) {
                        strcpy_k(write_iter, ":.,.");
                        write_iter = &(write_iter[2 + 2 * (is_x & (~cur_is_male))]);
                      }
                    }
                  } else if (cur_geno == 1) {
                    const AlleleCode ac = *patch_01_vals_iter++;
                    write_iter = AppendVcfMultiallelicDsForce01(allele_ct_m2, hds_force, ac, !should_be_diploid, write_iter);
                  } else {
                    const AlleleCode ac0 = *patch_10_vals_iter++;
                    const AlleleCode ac1 = *patch_10_vals_iter++;
                    if (ac0 != ac1) {
                      write_iter = AppendVcfMultiallelicDsForce10Het(allele_ct_m2, hds_force, ac0, ac1, !should_be_diploid, write_iter);
                    } else if (haploid_genotext_blen[cur_geno + cur_is_male * 4] == 2) {
                      write_iter = AppendVcfMultiallelicDsForce10Haploid(allele_ct_m2, hds_force, ac0, write_iter);
                    } else {
                      write_iter = AppendVcfMultiallelicDsForce10HomDiploid(allele_ct_m2, hds_force, ac0, '/', write_iter);
                    }
                  }
                  genovec_word >>= 2;
                  multiallelic_hw >>= 1;
                }
              }
            }
          }
        } else {
          // multiallelic, phased
          reterr = PgrGetMP(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, &pgv);
          if (unlikely(reterr)) {
            goto ExportVcf_ret_PGR_FAIL;
          }
          if (!pgv.patch_01_ct) {
            ZeroWArr(sample_ctl, pgv.patch_01_set);
          }
          if (!pgv.patch_10_ct) {
            ZeroWArr(sample_ctl, pgv.patch_10_set);
          }
          const AlleleCode* patch_01_vals_iter = pgv.patch_01_vals;
          const AlleleCode* patch_10_vals_iter = pgv.patch_10_vals;
          uint32_t phasepresent_hw = 0;
          if (!ds_force) {
            if (!is_haploid) {
              if (allele_ct <= 10) {
                uint32_t* write_iter_u32_alias = R_CAST(uint32_t*, write_iter);
                for (uint32_t widx = 0; ; ++widx) {
                  if (widx >= sample_ctl2_m1) {
                    if (widx > sample_ctl2_m1) {
                      break;
                    }
                    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                  }
                  uintptr_t genovec_word = pgv.genovec[widx];
                  const uint32_t multiallelic_hw = (R_CAST(const Halfword*, pgv.patch_01_set)[widx]) | (R_CAST(const Halfword*, pgv.patch_10_set)[widx]);
                  uint32_t prev_phased_halfword = R_CAST(Halfword*, prev_phased)[widx];

                  if (pgv.phasepresent_ct) {
                    phasepresent_hw = R_CAST(Halfword*, pgv.phasepresent)[widx];
                  }
                  const uint32_t phaseinfo_hw = R_CAST(Halfword*, pgv.phaseinfo)[widx];
                  uint32_t cur_shift = 1;
                  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                    const uintptr_t cur_geno = genovec_word & 3;
                    uint32_t cur_basic_genotext;
                    if (!(multiallelic_hw & cur_shift)) {
                      // usually "\t0/0", etc.
                      cur_basic_genotext = basic_genotext[cur_geno];
                      if (cur_geno == 1) {
                        if (phasepresent_hw & cur_shift) {
                          prev_phased_halfword |= cur_shift;
                          if (phaseinfo_hw & cur_shift) {
                            cur_basic_genotext ^= 0x1000100;  // 0|1 -> 1|0
                          }
                        } else {
                          prev_phased_halfword &= ~cur_shift;
                        }
                      }
                    } else {
                      AlleleCode ac0;
                      AlleleCode ac1;
                      if (cur_geno == 1) {
                        ac0 = 0;
                        ac1 = *patch_01_vals_iter++;
                      } else {
                        ac0 = *patch_10_vals_iter++;
                        ac1 = *patch_10_vals_iter++;
                      }
                      if (ac0 != ac1) {
                        if (phasepresent_hw & cur_shift) {
                          prev_phased_halfword |= cur_shift;
                          if (phaseinfo_hw & cur_shift) {
                            const AlleleCode ac_swap = ac0;
                            ac0 = ac1;
                            ac1 = ac_swap;
                          }
                        } else {
                          prev_phased_halfword &= ~cur_shift;
                        }
                      }
                      cur_basic_genotext = 0x302f3009 + (ac0 * 256) + (ac1 * 0x1000000);
                    }
                    // '/' = ascii 47, '|' = ascii 124
                    *write_iter_u32_alias++ = cur_basic_genotext + 0x4d0000 * ((prev_phased_halfword >> sample_idx_lowbits) & 1);
                    genovec_word >>= 2;
                    cur_shift = cur_shift * 2;
                  }
                  R_CAST(Halfword*, prev_phased)[widx] = prev_phased_halfword;
                }
                write_iter = R_CAST(char*, write_iter_u32_alias);
              } else {
                // phased, allele_ct > 10, !is_haploid, !ds_force
                for (uint32_t widx = 0; ; ++widx) {
                  if (widx >= sample_ctl2_m1) {
                    if (widx > sample_ctl2_m1) {
                      break;
                    }
                    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                  }
                  uintptr_t genovec_word = pgv.genovec[widx];
                  const uint32_t multiallelic_hw = (R_CAST(const Halfword*, pgv.patch_01_set)[widx]) | (R_CAST(const Halfword*, pgv.patch_10_set)[widx]);
                  uint32_t prev_phased_halfword = R_CAST(Halfword*, prev_phased)[widx];

                  if (pgv.phasepresent_ct) {
                    phasepresent_hw = R_CAST(Halfword*, pgv.phasepresent)[widx];
                  }

                  const uint32_t phaseinfo_hw = R_CAST(Halfword*, pgv.phaseinfo)[widx];
                  uint32_t cur_shift = 1;
                  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                    const uint32_t cur_geno = genovec_word & 3;
                    if (!(multiallelic_hw & cur_shift)) {
                      write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
                      if (cur_geno == 1) {
                        if (phasepresent_hw & cur_shift) {
                          prev_phased_halfword |= cur_shift;
                          if (phaseinfo_hw & cur_shift) {
                            memcpy(&(write_iter[-4]), "\t1|0", 4);
                          }
                        } else {
                          prev_phased_halfword &= ~cur_shift;
                        }
                      }
                      if (prev_phased_halfword & cur_shift) {
                        write_iter[-2] = '|';
                      }
                    } else {
                      AlleleCode ac0;
                      AlleleCode ac1;
                      if (cur_geno == 1) {
                        ac0 = 0;
                        ac1 = *patch_01_vals_iter++;
                      } else {
                        ac0 = *patch_10_vals_iter++;
                        ac1 = *patch_10_vals_iter++;
                      }
                      if (ac0 != ac1) {
                        if (phasepresent_hw & cur_shift) {
                          prev_phased_halfword |= cur_shift;
                          if (phaseinfo_hw & cur_shift) {
                            const AlleleCode ac_swap = ac0;
                            ac0 = ac1;
                            ac1 = ac_swap;
                          }
                        } else {
                          prev_phased_halfword &= ~cur_shift;
                        }
                      }
                      *write_iter++ = '\t';
                      write_iter = u32toa(ac0, write_iter);
                      if (prev_phased_halfword & cur_shift) {
                        *write_iter++ = '|';
                      } else {
                        *write_iter++ = '/';
                      }
                      write_iter = u32toa(ac1, write_iter);
                    }
                    genovec_word >>= 2;
                    cur_shift = cur_shift * 2;
                  }
                  R_CAST(Halfword*, prev_phased)[widx] = prev_phased_halfword;
                }
              }
            } else {
              // phased, at least partially haploid, !ds_force
              // probable todo: merge this with biallelic code, do the same for
              // other less-common biallelic-variable-width cases
              uint32_t is_male_hw = 0;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                if (is_x) {
                  is_male_hw = R_CAST(const Halfword*, sex_male_collapsed)[widx];
                }
                const uint32_t multiallelic_hw = (R_CAST(const Halfword*, pgv.patch_01_set)[widx]) | (R_CAST(const Halfword*, pgv.patch_10_set)[widx]);
                uint32_t prev_phased_halfword = R_CAST(Halfword*, prev_phased)[widx];

                if (pgv.phasepresent_ct) {
                  phasepresent_hw = R_CAST(Halfword*, pgv.phasepresent)[widx];
                }

                const uint32_t phaseinfo_hw = R_CAST(Halfword*, pgv.phaseinfo)[widx];
                uint32_t cur_shift = 1;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  const uint32_t cur_is_male = is_male_hw & 1;
                  if (!(multiallelic_hw & cur_shift)) {
                    const uint32_t cur_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                    memcpy(write_iter, &(basic_genotext[cur_geno]), 4);
                    write_iter = &(write_iter[cur_blen]);
                    if (cur_blen == 4) {
                      if (cur_geno == 1) {
                        if (phasepresent_hw & cur_shift) {
                          prev_phased_halfword |= cur_shift;
                          if (phaseinfo_hw & cur_shift) {
                            memcpy(&(write_iter[-4]), "\t1|0", 4);
                          } else {
                            write_iter[-2] = '|';
                          }
                        } else {
                          prev_phased_halfword &= ~cur_shift;
                        }
                      } else if ((prev_phased_halfword >> sample_idx_lowbits) & 1) {
                        write_iter[-2] = '|';
                      }
                    }
                  } else {
                    AlleleCode ac0;
                    AlleleCode ac1;
                    if (cur_geno == 1) {
                      ac0 = 0;
                      ac1 = *patch_01_vals_iter++;
                    } else {
                      ac0 = *patch_10_vals_iter++;
                      ac1 = *patch_10_vals_iter++;
                    }
                    *write_iter++ = '\t';
                    if ((ac0 != ac1) || (haploid_genotext_blen[2 + cur_is_male * 4] == 4)) {
                      if (ac0 != ac1) {
                        if (phasepresent_hw & cur_shift) {
                          prev_phased_halfword |= cur_shift;
                          if (phaseinfo_hw & cur_shift) {
                            const AlleleCode ac_swap = ac0;
                            ac0 = ac1;
                            ac1 = ac_swap;
                          }
                        } else {
                          prev_phased_halfword &= ~cur_shift;
                        }
                      }
                      write_iter = u32toa(ac0, write_iter);
                      if (prev_phased_halfword & cur_shift) {
                        *write_iter++ = '|';
                      } else {
                        *write_iter++ = '/';
                      }
                    }
                    write_iter = u32toa(ac1, write_iter);
                  }
                  genovec_word >>= 2;
                  is_male_hw >>= 1;
                  cur_shift = cur_shift * 2;
                }
                R_CAST(Halfword*, prev_phased)[widx] = prev_phased_halfword;
              }
            }
          } else {
            // phased, ds_force
            write_iter = strcpya_k(write_iter, ":DS");
            if (hds_force) {
              write_iter = strcpya_k(write_iter, ":HDS");
            }
            const uint32_t allele_ct_m2 = allele_ct - 2;
            if (!is_haploid) {
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                const uint32_t multiallelic_hw = (R_CAST(const Halfword*, pgv.patch_01_set)[widx]) | (R_CAST(const Halfword*, pgv.patch_10_set)[widx]);
                uint32_t prev_phased_halfword = R_CAST(Halfword*, prev_phased)[widx];

                if (pgv.phasepresent_ct) {
                  phasepresent_hw = R_CAST(Halfword*, pgv.phasepresent)[widx];
                }

                const uint32_t phaseinfo_hw = R_CAST(Halfword*, pgv.phaseinfo)[widx];
                uint32_t cur_shift = 1;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  if (!(multiallelic_hw & cur_shift)) {
                    write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
                    if (cur_geno == 1) {
                      if (phasepresent_hw & cur_shift) {
                        prev_phased_halfword |= cur_shift;
                        if (phaseinfo_hw & cur_shift) {
                          memcpy(&(write_iter[-4]), "\t1|0", 4);
                        }
                      } else {
                        prev_phased_halfword &= ~cur_shift;
                      }
                    }
                    if (prev_phased_halfword & cur_shift) {
                      write_iter[-2] = '|';
                    }
                    // DS
                    write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno]), 2);
                    if (cur_geno != 3) {
                      // repeat ",0" if nonmissing
                      write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                      if (hds_force) {
                        uint32_t hds_inttext_index = cur_geno;
                        if (phasepresent_hw & cur_shift) {
                          hds_inttext_index = 4 + ((phaseinfo_hw >> sample_idx_lowbits) & 1);
                        }
                        const uint64_t cur_inttext = hds_inttext[hds_inttext_index];
                        const uint32_t hap_blen = hds_inttext_blen[hds_inttext_index] / 2;
                        *R_CAST(uint32_t*, write_iter) = cur_inttext;
                        write_iter = &(write_iter[hap_blen]);
                        write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                        *R_CAST(uint32_t*, write_iter) = cur_inttext >> (8 * hap_blen);
                        write_iter = &(write_iter[hap_blen]);
                        write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                      }
                    } else {
                      if (hds_force) {
                        write_iter = strcpya_k(write_iter, ":.,.");
                      }
                    }
                  } else {
                    AlleleCode ac0;
                    AlleleCode ac1;
                    if (cur_geno == 1) {
                      ac0 = 0;
                      ac1 = *patch_01_vals_iter++;
                    } else {
                      ac0 = *patch_10_vals_iter++;
                      ac1 = *patch_10_vals_iter++;
                    }
                    if (ac0 != ac1) {
                      if (phasepresent_hw & cur_shift) {
                        prev_phased_halfword |= cur_shift;
                        write_iter = AppendVcfMultiallelicDsForcePhased(allele_ct_m2, hds_force, ac0, ac1, phaseinfo_hw & cur_shift, write_iter);
                      } else {
                        prev_phased_halfword &= ~cur_shift;
                        if (!ac0) {
                          write_iter = AppendVcfMultiallelicDsForce01(allele_ct_m2, hds_force, ac1, 0, write_iter);
                        } else {
                          write_iter = AppendVcfMultiallelicDsForce10Het(allele_ct_m2, hds_force, ac0, ac1, 0, write_iter);
                        }
                      }
                    } else {
                      write_iter = AppendVcfMultiallelicDsForce10HomDiploid(allele_ct_m2, hds_force, ac0, (prev_phased_halfword & cur_shift)? '|' : '/', write_iter);
                    }
                  }
                  genovec_word >>= 2;
                  cur_shift <<= 1;
                }
                R_CAST(Halfword*, prev_phased)[widx] = prev_phased_halfword;
              }
            } else {
              // dosage (or {H}DS-force) and phase present, partly/fully
              // haploid
              uint32_t is_male_hw = 0;
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                }
                uintptr_t genovec_word = pgv.genovec[widx];
                if (is_x) {
                  is_male_hw = R_CAST(const Halfword*, sex_male_collapsed)[widx];
                }
                const uint32_t multiallelic_hw = (R_CAST(const Halfword*, pgv.patch_01_set)[widx]) | (R_CAST(const Halfword*, pgv.patch_10_set)[widx]);
                uint32_t prev_phased_halfword = R_CAST(Halfword*, prev_phased)[widx];

                if (pgv.phasepresent_ct) {
                  phasepresent_hw = R_CAST(Halfword*, pgv.phasepresent)[widx];
                }

                const uint32_t phaseinfo_hw = R_CAST(Halfword*, pgv.phaseinfo)[widx];
                uint32_t cur_shift = 1;
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  const uint32_t cur_is_male = is_male_hw & 1;
                  const uint32_t should_be_diploid = is_x && (!cur_is_male);
                  // bugfix (29 Dec 2018): need to check correct bit here
                  if (!(multiallelic_hw & cur_shift)) {
                    const uint32_t cur_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                    memcpy(write_iter, &(basic_genotext[cur_geno]), 4);
                    write_iter = &(write_iter[cur_blen]);
                    if (cur_blen == 4) {
                      // render current hardcall as diploid (chrX nonmale, or
                      // het haploid)
                      if (cur_geno == 1) {
                        if (phasepresent_hw & cur_shift) {
                          prev_phased_halfword |= cur_shift;
                          if (phaseinfo_hw & cur_shift) {
                            memcpy(&(write_iter[-4]), "\t1|0", 4);
                          }
                        } else {
                          prev_phased_halfword &= ~cur_shift;
                        }
                      }
                      if (prev_phased_halfword & cur_shift) {
                        write_iter[-2] = '|';
                      }
                      // DS
                      write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno]), 2);
                      if (cur_geno != 3) {
                        // repeat ",0"
                        write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                        if (hds_force) {
                          uint32_t hds_inttext_index = cur_geno;
                          if (phasepresent_hw & cur_shift) {
                            hds_inttext_index = 4 + ((phaseinfo_hw >> sample_idx_lowbits) & 1);
                          }
                          const uint64_t cur_inttext = hds_inttext[hds_inttext_index];
                          const uint32_t hap_blen = hds_inttext_blen[hds_inttext_index] / 2;
                          *R_CAST(uint32_t*, write_iter) = cur_inttext;
                          write_iter = &(write_iter[hap_blen]);
                          write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                          // Don't render unphased het haploid as diploid here
                          if (should_be_diploid || (phasepresent_hw & cur_shift)) {
                            *R_CAST(uint32_t*, write_iter) = cur_inttext >> (8 * hap_blen);
                            write_iter = &(write_iter[hap_blen]);
                            write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                          }
                        }
                      } else {
                        if (hds_force) {
                          strcpy_k(write_iter, ":.,.");
                          write_iter = &(write_iter[2 + 2 * should_be_diploid]);
                        }
                      }
                    } else {
                      // render current hardcall as haploid
                      // (can't get here for hardcall-phased)
                      write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno + 4]), 2);
                      if (cur_geno != 3) {
                        // repeat ",0"
                        write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                        if (hds_force) {
                          // don't need to perform generic copy from
                          // hds_inttext, since only possibilities are :0 and
                          // :1
                          *R_CAST(uint16_t*, write_iter) = 0x303a + 128 * cur_geno;
                          write_iter = &(write_iter[2]);
                          write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                        }
                      } else {
                        if (hds_force) {
                          write_iter = strcpya_k(write_iter, ":.");
                        }
                      }
                    }
                  } else {
                    AlleleCode ac0;
                    AlleleCode ac1;
                    if (cur_geno == 1) {
                      ac0 = 0;
                      ac1 = *patch_01_vals_iter++;
                    } else {
                      ac0 = *patch_10_vals_iter++;
                      ac1 = *patch_10_vals_iter++;
                    }
                    if (ac0 != ac1) {
                      if (phasepresent_hw & cur_shift) {
                        prev_phased_halfword |= cur_shift;
                        write_iter = AppendVcfMultiallelicDsForcePhased(allele_ct_m2, hds_force, ac0, ac1, phaseinfo_hw & cur_shift, write_iter);
                      } else {
                        prev_phased_halfword &= ~cur_shift;
                        if (!ac0) {
                          write_iter = AppendVcfMultiallelicDsForce01(allele_ct_m2, hds_force, ac1, !should_be_diploid, write_iter);
                        } else {
                          write_iter = AppendVcfMultiallelicDsForce10Het(allele_ct_m2, hds_force, ac0, ac1, !should_be_diploid, write_iter);
                        }
                      }
                    } else {
                      if (should_be_diploid) {
                        write_iter = AppendVcfMultiallelicDsForce10HomDiploid(allele_ct_m2, hds_force, ac0, (prev_phased_halfword & cur_shift)? '|' : '/', write_iter);
                      } else {
                        write_iter = AppendVcfMultiallelicDsForce10Haploid(allele_ct_m2, hds_force, ac0, write_iter);
                      }
                    }
                  }
                  genovec_word >>= 2;
                  is_male_hw >>= 1;
                  cur_shift <<= 1;
                }
                R_CAST(Halfword*, prev_phased)[widx] = prev_phased_halfword;
              }
            }
          }
        }
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
        goto ExportVcf_ret_WRITE_FAIL;
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
    if (unlikely(bgzfclose_flush(writebuf_flush, write_iter, &bgzf, &reterr))) {
      goto ExportVcf_ret_1;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
    if (invalid_allele_code_seen) {
      logerrputs("Warning: At least one VCF allele code violates the official specification;\nother tools may not accept the file.  (Valid codes must either start with a\n'<', only contain characters in {A,C,G,T,N,a,c,g,t,n}, be an isolated '*', or\nrepresent a breakend.)\n");
    }
  }
  while (0) {
  ExportVcf_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportVcf_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind(pvar_info_reload, &pvar_reload_txs, &reterr);
    break;
  ExportVcf_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportVcf_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ExportVcf_ret_INCONSISTENT_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ExportVcf_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  }
 ExportVcf_ret_1:
  CleanupTextStream2(pvar_info_reload, &pvar_reload_txs, &reterr);
  CleanupBgzfCompressStream(&bgzf, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

static const Dosage kGenoToDosage[4] = {0, 16384, 32768, 65535};

typedef struct DosageTransposeCtxStruct {
  const uintptr_t* variant_include;
  const char* const* export_allele_missing;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select);
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  uint32_t sample_ct;
  uintptr_t stride;

  PgenReader** pgr_ptrs;
  uint32_t* read_variant_uidx_starts;
  uint32_t* write_vidx_starts;
  Dosage* smaj_dosagebuf;

  uint32_t cur_block_write_ct;

  uintptr_t** thread_write_genovecs;
  uintptr_t** thread_write_dosagepresents;
  Dosage** thread_write_dosagevals;

  uint64_t err_info;
} DosageTransposeCtx;

THREAD_FUNC_DECL DosageTransposeThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  DosageTransposeCtx* ctx = S_CAST(DosageTransposeCtx*, arg->sharedp->context);

  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t sample_ctd4 = sample_ct / 4;
  const uint32_t sample_rem = sample_ct % 4;
  const uintptr_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
  const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  const uint32_t sample_ctab2 = kBytesPerWord * sample_ctaw2;
  const uintptr_t stride = ctx->stride;
  const uintptr_t* variant_include = ctx->variant_include;
  const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select) = ctx->refalt1_select;
  const char* const* export_allele_missing = ctx->export_allele_missing;
  const uintptr_t* sample_include = ctx->sample_include;
  const uint32_t* sample_include_cumulative_popcounts = ctx->sample_include_cumulative_popcounts;

  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  uintptr_t* genovec_buf = ctx->thread_write_genovecs[tidx];
  uintptr_t* dosagepresent_buf = ctx->thread_write_dosagepresents[tidx];
  Dosage* dosagevals_buf = ctx->thread_write_dosagevals[tidx];
  uint32_t ref_allele_idx = 0;
  uint64_t new_err_info;
  do {
    const uint32_t cur_block_write_ct = ctx->cur_block_write_ct;
    const uint32_t vidx_end = ctx->write_vidx_starts[tidx + 1];
    uint32_t vidx_start = ctx->write_vidx_starts[tidx];
    if (cur_block_write_ct && (vidx_end != vidx_start)) {
      uintptr_t variant_uidx_base;
      uintptr_t variant_include_bits;
      BitIter1Start(variant_include, ctx->read_variant_uidx_starts[tidx], &variant_uidx_base, &variant_include_bits);
      Dosage* smaj_dosagebuf_iter = &(ctx->smaj_dosagebuf[vidx_start]);
      uint32_t dosage_cts[kDosagePerCacheline];
      do {
        uint32_t vidx_block_end = RoundDownPow2(vidx_start, kDosagePerCacheline) + kDosagePerCacheline;
        if (vidx_block_end > vidx_end) {
          vidx_block_end = vidx_end;
        }
        const uint32_t vidx_block_size = vidx_block_end - vidx_start;

        // part 1: decompress data
        uintptr_t* genovec_iter = genovec_buf;
        uintptr_t* dosage_present_iter = dosagepresent_buf;
        Dosage* dosage_main_iter = dosagevals_buf;
        for (uint32_t vidx_offset = 0; vidx_offset != vidx_block_size; ++vidx_offset) {
          const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
          if (refalt1_select) {
            ref_allele_idx = refalt1_select[variant_uidx][0];
          }
          uint32_t dosage_ct;
          const PglErr reterr = PgrGet1D(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, ref_allele_idx, pgrp, genovec_iter, dosage_present_iter, R_CAST(uint16_t*, dosage_main_iter), &dosage_ct);
          if (unlikely(reterr)) {
            new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
            goto DosageTransposeThread_err;
          }
          if (export_allele_missing && export_allele_missing[variant_uidx]) {
            // only keep missing vs. nonmissing distinction, zero everything
            // else out
            const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
            if (dosage_ct) {
              const Halfword* dosage_present_hw = R_CAST(const Halfword*, dosage_present_iter);
              for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
                genovec_iter[widx] = (Word11(genovec_iter[widx]) & (~UnpackHalfwordToWord(dosage_present_hw[widx]))) * 3;
              }
              dosage_ct = 0;
            } else {
              for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
                genovec_iter[widx] = Word11(genovec_iter[widx]) * 3;
              }
            }
          }
          genovec_iter = &(genovec_iter[sample_ctaw2]);
          dosage_present_iter = &(dosage_present_iter[sample_ctaw]);
          dosage_main_iter = &(dosage_main_iter[sample_ct]);
          dosage_cts[vidx_offset] = dosage_ct;
        }

        // part 2: process hardcalls for 4 samples at a time
        Dosage* dosagebuf_write_iter0 = smaj_dosagebuf_iter;
        for (uint32_t sample4_idx = 0; sample4_idx != sample_ctd4; ++sample4_idx) {
          Dosage* dosagebuf_write_iter1 = &(dosagebuf_write_iter0[stride]);
          Dosage* dosagebuf_write_iter2 = &(dosagebuf_write_iter1[stride]);
          Dosage* dosagebuf_write_iter3 = &(dosagebuf_write_iter2[stride]);
          const unsigned char* geno_read_iter = &(R_CAST(const unsigned char*, genovec_buf)[sample4_idx]);
          for (uint32_t vidx_offset = 0; vidx_offset != vidx_block_size; ++vidx_offset) {
            uint32_t cur_geno = *geno_read_iter;
            dosagebuf_write_iter0[vidx_offset] = kGenoToDosage[cur_geno & 3];
            dosagebuf_write_iter1[vidx_offset] = kGenoToDosage[(cur_geno >> 2) & 3];
            dosagebuf_write_iter2[vidx_offset] = kGenoToDosage[(cur_geno >> 4) & 3];
            dosagebuf_write_iter3[vidx_offset] = kGenoToDosage[(cur_geno >> 6) & 3];
            geno_read_iter = &(geno_read_iter[sample_ctab2]);
          }
          dosagebuf_write_iter0 = &(dosagebuf_write_iter3[stride]);
        }
        if (sample_rem) {
          const unsigned char* geno_read_iter = &(R_CAST(const unsigned char*, genovec_buf)[sample_ctd4]);
          for (uint32_t vidx_offset = 0; vidx_offset != vidx_block_size; ++vidx_offset) {
            uint32_t cur_geno = *geno_read_iter;
            Dosage* dosagebuf_write_iterx = &(dosagebuf_write_iter0[vidx_offset]);
            for (uint32_t sample_idx_lowbits = 0; ; ) {
              *dosagebuf_write_iterx = kGenoToDosage[cur_geno & 3];
              if (++sample_idx_lowbits == sample_rem) {
                break;
              }
              cur_geno >>= 2;
              dosagebuf_write_iterx = &(dosagebuf_write_iterx[stride]);
            }
            geno_read_iter = &(geno_read_iter[sample_ctab2]);
          }
        }
        // part 3: patch in dosages
        for (uint32_t vidx_offset = 0; vidx_offset != vidx_block_size; ++vidx_offset) {
          const uint32_t cur_dosage_ct = dosage_cts[vidx_offset];
          if (cur_dosage_ct) {
            const uintptr_t* dosage_present = &(dosagepresent_buf[vidx_offset * sample_ctaw]);
            const Dosage* dosage_main = &(dosagevals_buf[vidx_offset * sample_ct]);
            Dosage* cur_dosage_write = &(smaj_dosagebuf_iter[vidx_offset]);
            uintptr_t sample_idx_base = 0;
            uintptr_t dosage_present_bits = dosage_present[0];
            for (uint32_t dosage_idx = 0; dosage_idx != cur_dosage_ct; ++dosage_idx) {
              const uintptr_t sample_idx = BitIter1(dosage_present, &sample_idx_base, &dosage_present_bits);
              cur_dosage_write[sample_idx * stride] = dosage_main[dosage_idx];
            }
          }
        }
        vidx_start = vidx_block_end;
        smaj_dosagebuf_iter = &(smaj_dosagebuf_iter[vidx_block_size]);
      } while (vidx_start != vidx_end);
    }
    while (0) {
    DosageTransposeThread_err:
      UpdateU64IfSmaller(new_err_info, &ctx->err_info);
      break;
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

static_assert(sizeof(Dosage) == 2, "Export012Smaj() needs to be updated.");
PglErr Export012Smaj(const char* outname, const uintptr_t* orig_sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(AlleleCode, 2, export_allele), const char* const* export_allele_missing, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t include_dom, uint32_t include_uncounted, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, char exportf_delim, PgenFileInfo* pgfip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  DosageTransposeCtx ctx;
  {
    // Write header line; then fully load-and-transpose the first X samples,
    // flush them, load-and-transpose the next X, etc.
    // Similar to ExportIndMajorBed() and Export012Vmaj().
    // Priority is making the with-dosage case work well, since plink 1.9
    // already handles the no-dosage case.
    // (possible todo: have a separate no-dosage fast path)
    if (unlikely(variant_ct * (1 + include_dom) > (kMaxLongLine - 4 * kMaxIdSlen - 64) / 8)) {
      snprintf(g_logbuf, kLogbufSize, "Error: Too many variants for --export A%s.  (Try to work with A-transpose\ninstead.)\n", include_dom? "D" : "");
      goto Export012Smaj_ret_INCONSISTENT_INPUT_2;
    }
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto Export012Smaj_ret_OPEN_FAIL;
    }
    char* writebuf = g_textbuf;
    if (max_allele_slen > kMaxMediumLine - 5) {
      if (unlikely(bigstack_alloc_c(kMaxMediumLine + 5 + max_allele_slen, &writebuf))) {
        goto Export012Smaj_ret_NOMEM;
      }
    }
    char* write_iter = writebuf;
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    write_iter = strcpya_k(write_iter, "FID");
    *write_iter++ = exportf_delim;
    write_iter = strcpya_k(write_iter, "IID");
    *write_iter++ = exportf_delim;
    write_iter = strcpya_k(write_iter, "PAT");
    *write_iter++ = exportf_delim;
    write_iter = strcpya_k(write_iter, "MAT");
    *write_iter++ = exportf_delim;
    write_iter = strcpya_k(write_iter, "SEX");
    *write_iter++ = exportf_delim;
    write_iter = strcpya_k(write_iter, "PHENOTYPE");
    uintptr_t variant_uidx_base = 0;
    uintptr_t variant_include_bits = variant_include[0];
    uint32_t exported_allele_idx = 0;
    uint32_t cur_allele_ct = 2;
    uint64_t bytes_written = 0;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &variant_include_bits);
      *write_iter++ = exportf_delim;
      const char* cur_var_id = variant_ids[variant_uidx];
      const uint32_t cur_slen = strlen(cur_var_id);
      write_iter = memcpyax(write_iter, cur_var_id, cur_slen, '_');
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      if (export_allele) {
        exported_allele_idx = export_allele[variant_uidx][0];
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (export_allele_missing && export_allele_missing[variant_uidx]) {
        write_iter = strcpya(write_iter, export_allele_missing[variant_uidx]);
        if (write_iter >= writebuf_flush) {
          bytes_written += write_iter - writebuf;
          if (unlikely(fwrite_flush2(writebuf_flush, outfile, &write_iter))) {
            goto Export012Smaj_ret_WRITE_FAIL;
          }
        }
        if (include_uncounted) {
          write_iter = strcpya_k(write_iter, "(/");
          write_iter = strcpya(write_iter, cur_alleles[exported_allele_idx]);
          if (write_iter >= writebuf_flush) {
            bytes_written += write_iter - writebuf;
            if (unlikely(fwrite_flush2(writebuf_flush, outfile, &write_iter))) {
              goto Export012Smaj_ret_WRITE_FAIL;
            }
          }
        }
      } else {
        write_iter = strcpya(write_iter, cur_alleles[exported_allele_idx]);
        if (write_iter >= writebuf_flush) {
          bytes_written += write_iter - writebuf;
          if (unlikely(fwrite_flush2(writebuf_flush, outfile, &write_iter))) {
            goto Export012Smaj_ret_WRITE_FAIL;
          }
        }
        if (include_uncounted) {
          write_iter = strcpya_k(write_iter, "(/");
        }
      }
      if (include_uncounted) {
        for (uint32_t allele_idx = 0; allele_idx != cur_allele_ct; ++allele_idx) {
          if (allele_idx == exported_allele_idx) {
            continue;
          }
          write_iter = strcpya(write_iter, cur_alleles[allele_idx]);
          if (write_iter >= writebuf_flush) {
            bytes_written += write_iter - writebuf;
            if (unlikely(fwrite_flush2(writebuf_flush, outfile, &write_iter))) {
              goto Export012Smaj_ret_WRITE_FAIL;
            }
          }
          *write_iter++ = ',';
        }
        write_iter[-1] = ')';
      }
      if (include_dom) {
        *write_iter++ = exportf_delim;
        write_iter = memcpya(write_iter, cur_var_id, cur_slen);
        write_iter = strcpya_k(write_iter, "_HET");
        if (write_iter >= writebuf_flush) {
          bytes_written += write_iter - writebuf;
          if (unlikely(fwrite_flush2(writebuf_flush, outfile, &write_iter))) {
            goto Export012Smaj_ret_WRITE_FAIL;
          }
        }
      }
    }
    AppendBinaryEoln(&write_iter);
    bytes_written += write_iter - writebuf;
    if (unlikely(bytes_written > kMaxLongLine)) {
      snprintf(g_logbuf, kLogbufSize, "Error: --export A%s header line too long (>2GiB).\n", include_dom? "D" : "");
      goto Export012Smaj_ret_INCONSISTENT_INPUT_2;
    }

    uintptr_t* pheno_nm = nullptr;
    uintptr_t* pheno_cc = nullptr;
    double* pheno_qt = nullptr;
    // .raw files don't support categorical phenotypes
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
    const char* legacy_output_missing_pheno = g_legacy_output_missing_pheno;
    const uint32_t lomp_slen = strlen(legacy_output_missing_pheno);

    // Initially had the main read thread also perform decompression and
    // inversion, and the worker threads only unpacked/transposed the data to a
    // sample-major dosage matrix; but this made it pointless to have more than
    // 1 worker thread.
    //
    // So, moving closer to ExportIndMajorBed() strategy:
    // * Main read thread only loads raw bytes with PgfiMultiread().
    // * Worker thread(s) select subsets of the loaded variants to process in a
    //   cacheline-aware manner: first worker's interval ends with variant_idx
    //   divisible by (kCacheline / sizeof(Dosage)) == 32 (unless it's the
    //   only active worker), and all other interval(s) start with variant_idx
    //   divisible by 32.

    // todo: check when this saturates
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct * kDosagePerCacheline > variant_ct) {
      calc_thread_ct = DivUp(variant_ct, kDosagePerCacheline);
    }
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    uint32_t read_block_size;

    // note that we only allow this to use 1/4 of remaining memory
    if (unlikely(PgenMtLoadInit(variant_include, raw_sample_ct, variant_ct, bigstack_left() / 4, pgr_alloc_cacheline_ct, 0, 0, 0, pgfip, &calc_thread_ct, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto Export012Smaj_ret_NOMEM;
    }

    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    if (unlikely(
            bigstack_alloc_w(raw_sample_ctl, &sample_include) ||
            bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
            bigstack_alloc_u32(calc_thread_ct + 1, &ctx.write_vidx_starts) ||
            bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_genovecs) ||
            bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_dosagepresents) ||
            bigstack_alloc_dosagep(calc_thread_ct, &ctx.thread_write_dosagevals))) {
      goto Export012Smaj_ret_NOMEM;
    }

    // Remaining memory byte requirements:
    //   calc_thread_ct * kDosagePerCacheline * kBytesPerWord *
    //     read_sample_ctaw2 for per-thread genovecs buffers
    //   calc_thread_ct * kDosagePerCacheline * kBytesPerWord *
    //     read_sample_ctaw for per-thread dosage_presents buffers
    //   calc_thread_ct * kDosagePerCacheline * read_sample_ct *
    //     sizeof(Dosage) for per-thread dosage_main buffers
    //   (dosage_ct buffers just go on the thread stacks)
    //   read_sample_ct * variant_ct * sizeof(Dosage) for ctx.smaj_dosagebuf
    // This is about
    //   read_sample_ct *
    //     (calc_thread_ct * kDosagePerCacheline * 2.375 + variant_ct * 2)
    uintptr_t bytes_avail = bigstack_left();
    // account for rounding
    if (unlikely(bytes_avail < kCacheline + calc_thread_ct * kDosagePerCacheline * (2 * kCacheline))) {
      goto Export012Smaj_ret_NOMEM;
    }
    bytes_avail -= kCacheline + calc_thread_ct * kDosagePerCacheline * (2 * kCacheline);
    uint32_t read_sample_ct = sample_ct;
    uint32_t pass_ct = 1;
    const uintptr_t bytes_per_sample = calc_thread_ct * (kDosagePerCacheline / 8) * (3LLU + 8 * sizeof(Dosage)) + variant_ct * sizeof(Dosage);
    if ((sample_ct * S_CAST(uint64_t, bytes_per_sample)) > bytes_avail) {
      read_sample_ct = bytes_avail / bytes_per_sample;
      if (unlikely(!read_sample_ct)) {
        goto Export012Smaj_ret_NOMEM;
      }
      if (read_sample_ct > 4) {
        read_sample_ct = RoundDownPow2(read_sample_ct, 4);
      }
      pass_ct = 1 + (sample_ct - 1) / read_sample_ct;
    }
    uintptr_t read_sample_ctaw = BitCtToAlignedWordCt(read_sample_ct);
    uintptr_t read_sample_ctaw2 = NypCtToAlignedWordCt(read_sample_ct);
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      ctx.thread_write_genovecs[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(kDosagePerCacheline * sizeof(intptr_t) * read_sample_ctaw2));
      ctx.thread_write_dosagepresents[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(kDosagePerCacheline * sizeof(intptr_t) * read_sample_ctaw));
      ctx.thread_write_dosagevals[tidx] = S_CAST(Dosage*, bigstack_alloc_raw(kDosagePerCacheline * sizeof(Dosage) * read_sample_ct));
    }
    ctx.variant_include = variant_include;
    ctx.refalt1_select = export_allele;
    ctx.export_allele_missing = export_allele_missing;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto Export012Smaj_ret_NOMEM;
    }
    ctx.sample_ct = read_sample_ct;
    ctx.stride = RoundUpPow2(variant_ct, kDosagePerCacheline);
    ctx.smaj_dosagebuf = S_CAST(Dosage*, bigstack_alloc_raw_rd(read_sample_ct * S_CAST(uintptr_t, ctx.stride) * sizeof(Dosage)));
    ctx.err_info = (~0LLU) << 32;
    SetThreadFuncAndData(DosageTransposeThread, &ctx, &tg);

    const char* sample_ids = piip->sii.sample_ids;
    const char* paternal_ids = piip->parental_id_info.paternal_ids;
    const char* maternal_ids = piip->parental_id_info.maternal_ids;
    const uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
    uint32_t sample_uidx_start = AdvTo1Bit(orig_sample_include, 0);
    for (uint32_t pass_idx = 0; pass_idx != pass_ct; ++pass_idx) {
      memcpy(sample_include, orig_sample_include, raw_sample_ctl * sizeof(intptr_t));
      if (sample_uidx_start) {
        ClearBitsNz(0, sample_uidx_start, sample_include);
      }
      uint32_t sample_uidx_end;
      if (pass_idx + 1 == pass_ct) {
        read_sample_ct = sample_ct - pass_idx * read_sample_ct;
        ctx.sample_ct = read_sample_ct;
        sample_uidx_end = raw_sample_ct;
        read_sample_ctaw = BitCtToAlignedWordCt(read_sample_ct);
        read_sample_ctaw2 = NypCtToAlignedWordCt(read_sample_ct);
      } else {
        sample_uidx_end = FindNth1BitFrom(orig_sample_include, sample_uidx_start + 1, read_sample_ct);
        ClearBitsNz(sample_uidx_end, raw_sample_ct, sample_include);
      }
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      ctx.sample_include = sample_include;
      ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
      if (pass_idx) {
        ReinitThreads(&tg);
        pgfip->block_base = main_loadbufs[0];
        PgrSetBaseAndOffset0(main_loadbufs[0], calc_thread_ct, ctx.pgr_ptrs);
      }
      putc_unlocked('\r', stdout);
      printf("--export A%s pass %u/%u: loading... 0%%", include_dom? "D" : "", pass_idx + 1, pass_ct);
      fflush(stdout);
      // Main workflow:
      // 1. Set n=0, load first calc_thread_ct * kDosagePerCacheline
      //    post-filtering variants
      //
      // 2. Spawn threads processing batch n
      // 3. Load batch (n+1) unless eof
      // 4. Join threads
      // 5. Increment n by 1
      // 6. Goto step 2 unless eof
      uint32_t parity = 0;
      uint32_t read_block_idx = 0;
      uint32_t pct = 0;
      uint32_t next_print_idx = variant_ct / 100;
      for (uint32_t variant_idx = 0; ; ) {
        const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
        if (unlikely(reterr)) {
          goto Export012Smaj_ret_PGR_FAIL;
        }
        if (variant_idx) {
          JoinThreads(&tg);
          reterr = S_CAST(PglErr, ctx.err_info);
          if (unlikely(reterr)) {
            goto Export012Smaj_ret_PGR_FAIL;
          }
        }
        if (!IsLastBlock(&tg)) {
          ctx.cur_block_write_ct = cur_block_write_ct;
          ComputePartitionAligned(variant_include, calc_thread_ct, read_block_idx * read_block_size, variant_idx, cur_block_write_ct, kDosagePerCacheline, ctx.read_variant_uidx_starts, ctx.write_vidx_starts);
          PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
          if (variant_idx + cur_block_write_ct == variant_ct) {
            DeclareLastThreadBlock(&tg);
          }
          if (unlikely(SpawnThreads(&tg))) {
            goto Export012Smaj_ret_THREAD_CREATE_FAIL;
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
      if (pct > 10) {
        fputs("\b \b", stdout);
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\bwriting... 0%", stdout);
      fflush(stdout);
      pct = 0;
      next_print_idx = read_sample_ct / 100;
      uintptr_t sample_uidx_base;
      uintptr_t sample_include_bits;
      BitIter1Start(sample_include, sample_uidx_start, &sample_uidx_base, &sample_include_bits);
      const Dosage* cur_dosage_row = ctx.smaj_dosagebuf;
      for (uint32_t sample_idx = 0; sample_idx != read_sample_ct; ++sample_idx) {
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
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto Export012Smaj_ret_WRITE_FAIL;
        }
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          *write_iter++ = exportf_delim;
          uint32_t cur_dosage_val = cur_dosage_row[variant_idx];
          if (cur_dosage_val != 65535) {
            write_iter = PrintSmallDosage(cur_dosage_val, write_iter);
            if (include_dom) {
              *write_iter++ = exportf_delim;
              write_iter = PrintSmallDosage(16384 - abs_i32(cur_dosage_val - 16384), write_iter);
            }
          } else {
            write_iter = strcpya_k(write_iter, "NA");
            if (include_dom) {
              *write_iter++ = exportf_delim;
              write_iter = strcpya_k(write_iter, "NA");
            }
          }
          // todo: try making this check less frequently
          if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
            goto Export012Smaj_ret_WRITE_FAIL;
          }
        }
        AppendBinaryEoln(&write_iter);
        cur_dosage_row = &(cur_dosage_row[ctx.stride]);
        if (sample_idx >= next_print_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (sample_idx * 100LLU) / read_sample_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_idx = (pct * S_CAST(uint64_t, read_sample_ct)) / 100;
        }
      }
      sample_uidx_start = sample_uidx_end;
      if (pct > 10) {
        fputs("\b \b", stdout);
      }
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto Export012Smaj_ret_WRITE_FAIL;
    }
    fputs("\b\bdone.\n", stdout);
    logprintfww("--export A%s: %s written.\n", include_dom? "D" : "", outname);
  }
  while (0) {
  Export012Smaj_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Export012Smaj_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  Export012Smaj_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  Export012Smaj_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  Export012Smaj_ret_INCONSISTENT_INPUT_2:
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  Export012Smaj_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr Exportf(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const STD_ARRAY_PTR_DECL(AlleleCode, 2, refalt1_select), const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const ExportfInfo* eip, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, MakePlink2Flags make_plink2_flags, uintptr_t pgr_alloc_cacheline_ct, char* xheader, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* sex_male_collapsed;
    if (unlikely(
            bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
            bigstack_alloc_w(sample_ctaw, &sex_male_collapsed))) {
      goto Exportf_ret_NOMEM;
    }
    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_male_collapsed);
    ZeroTrailingWords(sample_ctl, sex_male_collapsed);
    uint32_t* sample_missing_geno_cts = nullptr;
    ExportfFlags flags = eip->flags;
    if (flags & (kfExportfOxGen | kfExportfHaps | kfExportfHapsLegend | kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13)) {
      if (unlikely(bigstack_alloc_u32(sample_ct, &sample_missing_geno_cts))) {
        goto Exportf_ret_NOMEM;
      }
    }
    if (flags & (kfExportf01 | kfExportf12)) {
      // todo
    }
    if (flags & (kfExportfTypemask - kfExportfIndMajorBed - kfExportfVcf - kfExportfOxGen - kfExportfBgen11 - kfExportfBgen12 - kfExportfBgen13 - kfExportfHaps - kfExportfHapsLegend - kfExportfATranspose - kfExportfA - kfExportfAD)) {
      logerrputs("Error: Only VCF, oxford, bgen-1.x, haps, hapslegend, A, AD, A-transpose, and\nind-major-bed output have been implemented so far.\n");
      reterr = kPglRetNotYetSupported;
      goto Exportf_ret_1;
    }
    const char exportf_delim = (flags & kfExportfSpaces)? ' ' : '\t';
    // Move this first, so we error out more quickly when any variants are
    // still multiallelic.
    if ((!(make_plink2_flags & kfMakeBim)) && (flags & kfExportfIndMajorBed)) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".bim");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteMapOrBim(outname, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, nullptr, refalt1_select, variant_cms, variant_ct, max_allele_slen, exportf_delim, 0, max_thread_ct);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
      logputs("done.\n");
    }

    const STD_ARRAY_PTR_DECL(AlleleCode, 2, export_allele) = refalt1_select;
    const char** export_allele_missing = nullptr;
    uint32_t max_export_allele_slen = max_allele_slen;
    if (eip->export_allele_fname) {
      STD_ARRAY_PTR_DECL(AlleleCode, 2, new_export_allele);
      if (unlikely(BIGSTACK_ALLOC_STD_ARRAY(AlleleCode, 2, raw_variant_ct, &new_export_allele))) {
        goto Exportf_ret_NOMEM;
      }
      if (refalt1_select) {
        memcpy(new_export_allele, refalt1_select, raw_variant_ct * 2 * sizeof(AlleleCode));
      } else {
        memset(new_export_allele, 0, raw_variant_ct * 2 * sizeof(AlleleCode));
      }
      reterr = ExportAlleleLoad(eip->export_allele_fname, variant_include, variant_ids, allele_idx_offsets, allele_storage, raw_variant_ct, variant_ct, max_variant_id_slen, max_thread_ct, &max_export_allele_slen, new_export_allele, &export_allele_missing);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
      export_allele = new_export_allele;
    }
    if (flags & kfExportfATranspose) {
      // multiallelic ok
      snprintf(outname_end, kMaxOutfnameExtBlen, ".traw");
      PgrClearLdCache(simple_pgrp);
      reterr = Export012Vmaj(outname, sample_include, sample_include_cumulative_popcounts, piip->sii.sample_ids, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, export_allele, export_allele_missing, variant_cms, sample_ct, piip->sii.max_sample_id_blen, variant_ct, max_export_allele_slen, exportf_delim, simple_pgrp);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & kfExportfIndMajorBed) {
      // multiallelic not ok, but already checked
      reterr = ExportIndMajorBed(sample_include, variant_include, refalt1_select, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, max_thread_ct, pgr_alloc_cacheline_ct, pgfip, outname, outname_end);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & kfExportfOxGen) {
      // multiallelic really not ok
      PgrClearLdCache(simple_pgrp);
      reterr = ExportOxGen(sample_include, sample_include_cumulative_popcounts, sex_male, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, refalt1_select, sample_ct, variant_ct, max_allele_slen, max_thread_ct, flags, simple_pgrp, outname, outname_end, sample_missing_geno_cts);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & (kfExportfHaps | kfExportfHapsLegend)) {
      // multiallelic not ok
      PgrClearLdCache(simple_pgrp);
      reterr = ExportOxHapslegend(sample_include, sample_include_cumulative_popcounts, sex_male_collapsed, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, refalt1_select, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, flags, simple_pgrp, outname, outname_end);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
      ZeroU32Arr(sample_ct, sample_missing_geno_cts);
    }
    const IdpasteFlags idpaste_flags = eip->idpaste_flags;
    const char id_delim = eip->id_delim;
    if (flags & kfExportfBgen11) {
      // multiallelic not ok
      assert(PopcountWords(sample_include, raw_sample_ctl) == sample_ct);
      snprintf(outname_end, kMaxOutfnameExtBlen, ".bgen");
      reterr = ExportBgen11(outname, sample_include, sample_include_cumulative_popcounts, sex_male, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, refalt1_select, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, flags, pgr_alloc_cacheline_ct, pgfip, sample_missing_geno_cts);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    } else if (flags & (kfExportfBgen12 | kfExportfBgen13)) {
      // multiallelic ok
      snprintf(outname_end, kMaxOutfnameExtBlen, ".bgen");
      reterr = ExportBgen13(outname, sample_include, sample_include_cumulative_popcounts, &(piip->sii), sex_nm, sex_male, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, refalt1_select, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, flags, eip->bgen_bits, idpaste_flags, id_delim, pgr_alloc_cacheline_ct, pgfip, sample_missing_geno_cts);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & (kfExportfOxGen | kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13 | kfExportfHaps | kfExportfHapsLegend)) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".sample");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      uint32_t y_ct = 0;
      const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
      if ((!IsI32Neg(y_code)) && IsSet(cip->chr_mask, y_code)) {
        y_ct = CountChrVariantsUnsafe(variant_include, cip, y_code);
      }
      assert(PopcountWords(sample_include, raw_sample_ctl) == sample_ct);
      reterr = ExportOxSample(outname, sample_include, piip->sii.sample_ids, sample_missing_geno_cts, sex_nm, sex_male, pheno_cols, pheno_names, sample_ct, piip->sii.max_sample_id_blen, pheno_ct, max_pheno_name_blen, variant_ct, y_ct);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
      logputs("done.\n");
    }
    if (flags & kfExportfVcf) {
      // multiallelic ok
      PgrClearLdCache(simple_pgrp);
      reterr = ExportVcf(sample_include, sample_include_cumulative_popcounts, &(piip->sii), sex_male_collapsed, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pvar_info_reload, xheader_blen, info_flags, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, max_thread_ct, flags, eip->vcf_mode, idpaste_flags, id_delim, xheader, pgfip, simple_pgrp, outname, outname_end);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    // todo: everything else
    // sample-major output should share a (probably multithreaded) transpose
    // routine
    if (flags & (kfExportfA | kfExportfAD)) {
      // multiallelic ok
      snprintf(outname_end, kMaxOutfnameExtBlen, ".raw");
      reterr = Export012Smaj(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, variant_include, variant_ids, allele_idx_offsets, allele_storage, export_allele, export_allele_missing, raw_sample_ct, sample_ct, pheno_ct, raw_variant_ct, variant_ct, max_export_allele_slen, (flags / kfExportfAD) & 1, (flags / kfExportfIncludeAlt) & 1, max_thread_ct, pgr_alloc_cacheline_ct, exportf_delim, pgfip);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }

    if ((!(make_plink2_flags & kfMakeFam)) && (flags & kfExportfIndMajorBed)) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".fam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteFam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, nullptr, sample_ct, pheno_ct, exportf_delim);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
      logputs("done.\n");
    }
  }
  while (0) {
  Exportf_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 Exportf_ret_1:
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
