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

#include "plink2_export.h"

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include "include/plink2_bgzf.h"
#include "include/plink2_bits.h"
#include "include/plink2_htable.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "plink2_cmdline.h"
#include "plink2_compress_stream.h"
#include "plink2_decompress.h"
#include "plink2_export_legacy.h"
#include "zstd/lib/zstd.h"

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

PglErr ExportAlleleLoad(const char* fname, const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, char input_missing_geno_char, uint32_t max_thread_ct, uint32_t* max_export_allele_slen_ptr, AlleleCode* export_allele, const char*** allele_missing_ptr) {
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
          if (strequal_unsafe(cur_alleles[allele_idx], allele_start, allele_slen)) {
            cur_allele_idx = allele_idx;
            break;
          }
        }
      }
      if (IsSet(already_seen, variant_uidx)) {
        // Only ok if allele is the same as before.
        ++duplicate_ct;
        if (cur_allele_idx == -1) {
          if (likely(allele_missing && allele_missing[variant_uidx] && strequal_unsafe(allele_missing[variant_uidx], allele_start, allele_slen))) {
            continue;
          }
        } else {
          if (likely(((!allele_missing) || (!allele_missing[variant_uidx])) && (S_CAST(uint32_t, cur_allele_idx) == export_allele[variant_uidx]))) {
            continue;
          }
        }
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --export-allele file has a variant ID duplicating an earlier line, and with a different allele code.\n", line_idx);
        goto ExportAlleleLoad_ret_INCONSISTENT_INPUT_WW;
      }
      SetBit(variant_uidx, already_seen);
      ++hit_ct;
      if (cur_allele_idx != -1) {
        export_allele[variant_uidx] = S_CAST(uint32_t, cur_allele_idx);
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
      if (StoreStringAtEndK(tmp_alloc_base, allele_start, allele_slen, &tmp_alloc_end, &(allele_missing[variant_uidx]))) {
        goto ExportAlleleLoad_ret_NOMEM;
      }
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

PglErr Export012Vmaj(const char* outname, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const char* sample_ids, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* export_allele, const char* const* export_allele_missing, const double* variant_cms, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t variant_ct, uint32_t max_allele_slen, char exportf_delim, PgenReader* simple_pgrp) {
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
    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf) ||
                 bigstack_alloc_c(kMaxMediumLine + max_chr_blen + 2 * kMaxIdSlen + 48 + 2 * max_allele_slen + (8 * k1LU) * sample_ct, &writebuf) ||
                 bigstack_alloc_w(sample_ctl2, &genovec))) {
      goto Export012Vmaj_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    const uint32_t dosage_is_present = PgrGetGflags(simple_pgrp) & kfPgenGlobalDosagePresent;
    uintptr_t* dosage_present = nullptr;
    Dosage* dosage_main = nullptr;
    if (dosage_is_present) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &dosage_present) ||
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
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    logprintfww5("--export Av to %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_blen = 0;
    uint32_t exported_allele_idx = 0;
    uint32_t cur_allele_ct = 2;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;

    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
        write_iter = dtoa_g_p8(variant_cms[variant_uidx], write_iter);
      } else {
        *write_iter++ = '0';
      }
      *write_iter++ = exportf_delim;
      write_iter = u32toa_x(variant_bps[variant_uidx], exportf_delim, write_iter);

      if (export_allele) {
        exported_allele_idx = export_allele[variant_uidx];
      }
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (export_allele_missing && export_allele_missing[variant_uidx]) {
        reterr = PgrGetMissingnessD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, nullptr, missingness_dosage, nullptr, genovec);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, variant_uidx);
          goto Export012Vmaj_ret_1;
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
        reterr = PgrGet1D(sample_include, pssi, sample_ct, variant_uidx, exported_allele_idx, simple_pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, variant_uidx);
          goto Export012Vmaj_ret_1;
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
            uint32_t dosage_present_hw = DowncastWToHW(dosage_present)[widx];
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
  Export012Vmaj_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 Export012Vmaj_ret_1:
  fclose_cond(outfile);
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

PglErr ExportOxGen(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_nonfemale_collapsed, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, uint32_t sample_ct, uint32_t variant_ct, uint32_t max_allele_slen, __maybe_unused uint32_t max_thread_ct, ExportfFlags exportf_flags, PgenReader* simple_pgrp, char* outname, char* outname_end, uint32_t* sample_missing_geno_cts) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  BgzfCompressStream bgzf;
  PreinitBgzfCompressStream(&bgzf);
  {
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    uintptr_t* genovec;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (unlikely(bigstack_alloc_w(sample_ctl2, &genovec))) {
      goto ExportOxGen_ret_NOMEM;
    }

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
      if (unlikely(bigstack_alloc_dosage(sample_ct * (1 + multiallelic_present), &dosage_main) ||
                   bigstack_alloc_w(sample_ctl, &dosage_present))) {
        goto ExportOxGen_ret_NOMEM;
      }
    }
    const uint32_t max_chr_blen = GetMaxChrSlen(cip) + 1;
    const uint32_t is_v2 = (exportf_flags / kfExportfOxGenV2) & 1;
    // if no dosages, all genotypes are 6 bytes (missing = " 0 0 0")
    // with dosages, we print up to 5 digits past the decimal point, so 7 bytes
    //   + space for each number, 24 bytes max
    const uintptr_t max_geno_slen = 6 + (dosage_present != nullptr) * 18;
    char* chr_buf;  // includes trailing space
    char* writebuf;
    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf) ||
                 bigstack_alloc_c(kMaxMediumLine + max_chr_blen + (kMaxIdSlen << is_v2) + 16 + 2 * max_allele_slen + max_geno_slen * sample_ct, &writebuf))) {
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

    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
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
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
      const char* variant_id = variant_ids[variant_uidx];
      const uint32_t variant_id_slen = strlen(variant_id);
      write_iter = memcpyax(write_iter, variant_id, variant_id_slen, ' ');
      if (is_v2) {
        write_iter = memcpyax(write_iter, variant_id, variant_id_slen, ' ');
      }
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
      if (allele_permute) {
        ref_allele_idx = allele_permute[allele_idx_offset_base];
      }
      // No dosage rescaling here, too messy to put that in more than one
      // place.
      uint32_t dosage_ct;
      reterr = PgrGetD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
      if (unlikely(reterr)) {
        PgenErrPrintNV(reterr, variant_uidx);
        goto ExportOxGen_ret_1;
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
        const Halfword* dosage_present_alias = DowncastWToHW(dosage_present);
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
        BitvecAnd(sex_nonfemale_collapsed, sample_ctl, missing_acc1);
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
    uint32_t* scrambled_missing_cts = DowncastVecWToU32(missing_acc32);
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = VcountScramble1(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  ExportOxGen_ret_NOMEM:
    reterr = kPglRetNomem;
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

PglErr ExportOxHapslegend(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male_collapsed, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, __maybe_unused uint32_t max_thread_ct, ExportfFlags exportf_flags, PgenReader* simple_pgrp, char* outname, char* outname_end) {
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
          // Multiallelic variants currently tolerated if allele_permute
          // specified, and only its top two alleles are ever present.
          if (unlikely((!allele_permute) && (allele_idx_offsets[variant_uidx + 1] != allele_idx_offset_base + 2))) {
            logputs("\n");
            logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
            goto ExportOxHapslegend_ret_INCONSISTENT_INPUT;
          }
        }
        if (allele_permute) {
          allele_idx0 = allele_permute[allele_idx_offset_base + ref_allele_last];
          allele_idx1 = allele_permute[allele_idx_offset_base + 1 - ref_allele_last];
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
    if (unlikely(bigstack_alloc_w(sample_ctv * kWordsPerVec, &sex_male_collapsed_interleaved) ||
                 bigstack_alloc_c(writebuf_alloc, &writebuf) ||
                 bigstack_alloc_w(sample_ctl2, &genovec) ||
                 bigstack_alloc_w(sample_ctl, &phasepresent) ||
                 bigstack_alloc_w(sample_ctl, &phaseinfo))) {
      goto ExportOxHapslegend_ret_NOMEM;
    }
    // sex_male_collapsed had trailing words zeroed out
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
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    uint32_t chr_blen = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
        if (unlikely((!allele_permute) && (allele_idx_offsets[variant_uidx + 1] != allele_idx_offset_base + 2))) {
          logputs("\n");
          logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
          goto ExportOxHapslegend_ret_INCONSISTENT_INPUT;
        }
      }
      if (allele_permute) {
        allele_idx0 = allele_permute[allele_idx_offset_base + ref_allele_last];
        allele_idx1 = allele_permute[allele_idx_offset_base + 1 - ref_allele_last];
      }
      uint32_t phasepresent_ct;
      reterr = PgrGet2P(sample_include, pssi, sample_ct, variant_uidx, allele_idx0, allele_idx1, simple_pgrp, genovec, phasepresent, phaseinfo, &phasepresent_ct);
      if (unlikely(reterr)) {
        PgenErrPrintNV(reterr, variant_uidx);
        goto ExportOxHapslegend_ret_1;
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
  const AlleleCode* allele_permute;
  const uintptr_t* sex_nonfemale_collapsed;
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
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  const uintptr_t* sex_nonfemale_collapsed = ctx->sex_nonfemale_collapsed;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uint32_t sample_ctl2_m1 = NypCtToWordCt(sample_ct) - 1;
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t bgen_geno_buf_blen = 6 * sample_ct;
  const uint32_t bgen_compressed_buf_max = ctx->bgen_compressed_buf_max;
  const AlleleCode* allele_permute = ctx->allele_permute;
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
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        if (allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base != 2) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, kPglRetInconsistentInput);
          goto ExportBgen11Thread_err;
        }
      }
      if (allele_permute) {
        ref_allele_idx = allele_permute[allele_idx_offset_base];
      }
      uint32_t dosage_ct;
      PglErr reterr = PgrGetD(sample_include, pssi, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_main, &dosage_ct);
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
        const Halfword* dosage_present_alias = DowncastWToHW(dosage_present);
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
        BitvecAnd(sex_nonfemale_collapsed, sample_ctl, missing_acc1);
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
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  {
    VcountIncr4To8(missing_acc4, acc4_vec_ct, missing_acc8);
    VcountIncr8To32(missing_acc8, acc8_vec_ct, missing_acc32);
  }
  while (0) {
  ExportBgen11Thread_err:
    UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    THREAD_BLOCK_FINISH(arg);
    break;
  }
  THREAD_RETURN;
}

PglErr ExportBgen11(const char* outname, const uintptr_t* sample_include, uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_nonfemale_collapsed, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, ExportfFlags exportf_flags, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uint32_t* sample_missing_geno_cts) {
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
    const uintptr_t writebuf_len = bgen_compressed_buf_max + 2 * max_allele_slen + 2 * kMaxIdSlen + 32;
    char* chr_buf;
    unsigned char* writebuf;
    if (unlikely(bigstack_alloc_c(max_chr_slen, &chr_buf) ||
                 bigstack_alloc_uc(writebuf_len, &writebuf))) {
      goto ExportBgen11_ret_NOMEM;
    }
    ctx.sex_nonfemale_collapsed = sex_nonfemale_collapsed;

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
    if (unlikely(bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(ctx.writebufs[0])) ||
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
    ctx.allele_permute = allele_permute;
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
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
          } else {
            PgenErrPrintNV(reterr, ctx.err_info >> 32);
          }
          goto ExportBgen11_ret_1;
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
            if (unlikely((!allele_permute) && (allele_idx_offsets[write_variant_uidx + 1] != allele_idx_offset_base + 2))) {
              logputs("\n");
              logerrprintfww("Error: %s cannot contain multiallelic variants.\n", outname);
              goto ExportBgen11_ret_INCONSISTENT_INPUT;
            }
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          if (allele_permute) {
            ref_allele_idx = allele_permute[allele_idx_offset_base];
            alt1_allele_idx = allele_permute[allele_idx_offset_base + 1];
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
    uint32_t* scrambled_missing_cts = DowncastWToU32(&(ctx.missing_acc1[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
      const uint32_t* thread_scrambled_missing_cts = DowncastWToU32(&(ctx.missing_acc1[tidx][acc32_offset]));
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
 ExportBgen11_ret_1:
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
  pgfip->block_base = nullptr;
  return reterr;
}

// For constant-ploidy cases.  Writes up to 7 bytes past the end.
// For what it's worth, bgen_geno_buf_iter is always word-aligned on entry (8
// bytes after beginning of buffer, which is at the start of a cacheline), and
// it's easy to force vector-alignment if it would be useful.
unsigned char* FillBgen13PloidyAndMissingness(const uintptr_t* genovec, const uintptr_t* dosage_present, uintptr_t ploidy, uint32_t sample_ct, uint32_t dosage_ct, unsigned char* bgen_geno_buf_iter) {
  // bugfix (30 Aug 2018): need to spell out mask in 32-bit case
  const uint64_t ploidy_u64 = ploidy * 0x101010101010101LLU;
  const uint32_t sample_ct8 = DivUp(sample_ct, 8);
  const uint16_t* genovec_alias = DowncastKWToU16(genovec);
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
      const uint64_t write_u64 = cur_missing8 + ploidy_u64;
      CopyToUnalignedOffsetU64(bgen_geno_buf_iter, &write_u64, write_widx);
    }
  } else {
    // don't bother with word-based 32-bit code because this case becomes
    // more annoying
    const unsigned char* dosage_present_alias = DowncastKToUc(dosage_present);
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
      const uint64_t write_u64 = cur_missing8 + ploidy_u64;
      CopyToUnalignedOffsetU64(bgen_geno_buf_iter, &write_u64, write_widx);
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
  const Halfword* sex_female_alias = DowncastKWToHW(sex_female);
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
static inline void AppendBits(uint32_t bit_ct, uintptr_t payload, uintptr_t* cur_write_bits_ptr, uint32_t* cur_write_bit_idx_ptr, unsigned char** probs_write_iter_ptr) {
  uint32_t cur_write_bit_idx = *cur_write_bit_idx_ptr;
  *cur_write_bits_ptr |= payload << cur_write_bit_idx;
  cur_write_bit_idx += bit_ct;
  if (cur_write_bit_idx >= kBitsPerWord) {
    AppendW(*cur_write_bits_ptr, probs_write_iter_ptr);
    cur_write_bit_idx -= kBitsPerWord;
    *cur_write_bits_ptr = payload >> (bit_ct - cur_write_bit_idx);
  }
  *cur_write_bit_idx_ptr = cur_write_bit_idx;
}

static inline void Append0Bits(uint32_t bit_ct, uintptr_t* cur_write_bits_ptr, uint32_t* cur_write_bit_idx_ptr, unsigned char** probs_write_iter_ptr) {
  uint32_t cur_write_bit_idx = *cur_write_bit_idx_ptr;
  cur_write_bit_idx += bit_ct;
  if (cur_write_bit_idx >= kBitsPerWord) {
    AppendW(*cur_write_bits_ptr, probs_write_iter_ptr);
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
    if (unlikely(bigstack_alloc_u64(256, &tablesp->diploid_hardcall_table8) ||
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
    if (unlikely(bigstack_alloc_u64(512, &tablesp->diploid_hardcall_table16) ||
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
  const AlleleCode* allele_permute;
  const uintptr_t* sex_male_collapsed;
  const uintptr_t* sex_female_collapsed;
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
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
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
  const AlleleCode* allele_permute = ctx->allele_permute;
  uint32_t chr_fo_idx = UINT32_MAX;  // deliberate overflow
  uint32_t chr_end = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;
  uint32_t cur_y = 0;

  uint32_t is_haploid = 0;  // includes chrX
  // Unlike the VCF/BCF export case, we do not save heterozygous haploid calls
  // as ploidy=2 unless there is also phase information.
  // Unlike VCF, chrY female (but not unknown-sex) ploidy is 0 when genotype is
  // missing.
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
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        // todo: multiallelic cases
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        if (allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base != 2) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, kPglRetInconsistentInput);
          goto ExportBgen13Thread_err;
        }
      }
      if (allele_permute) {
        ref_allele_idx = allele_permute[allele_idx_offset_base];
      }
      // (note that multiallelic variants should always be stored using the
      // phased format; just set both sides equal when the call is actually
      // unphased.  this way there's no O([allele count]^2) bloat problem.)
      PglErr reterr = PgrGetDp(sample_include, pssi, sample_ct, variant_uidx, pgrp, &pgv);
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
        *bgen_geno_buf_iter++ = use_phased_format;
        *bgen_geno_buf_iter++ = bit_precision;
        const unsigned char* genovec_alias = DowncastKToUc(pgv.genovec);
        unsigned char* probs_write_citer = bgen_geno_buf_iter;
        if (!use_phased_format) {
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
          const unsigned char* phasepresent_alias = DowncastKToUc(pgv.phasepresent);
          const unsigned char* phaseinfo_alias = DowncastKToUc(pgv.phaseinfo);
          const uint32_t bit_precision_x2 = bit_precision * 2;

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
              AppendBits(bit_precision_x2, payload, &cur_write_bits, &cur_write_bit_idx, &probs_write_citer);
              cur_geno4 >>= 2;
              cur_phasepresent4 >>= 1;
              cur_phaseinfo4 >>= 1;
            }
          }
          CopyToUnalignedW(probs_write_citer, &cur_write_bits);
        }
        if (pgv.dosage_ct) {
          const uintptr_t bit_precision_x2 = bit_precision * 2;
          uintptr_t sample_uidx_base = 0;
          uintptr_t dosage_present_bits = pgv.dosage_present[0];
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
              CopyBits(output_prob2 | (output_prob1 << bit_precision), sample_uidx * S_CAST(uint64_t, bit_precision_x2), bit_precision_x2, bgen_geno_buf_iter);
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
              CopyBits(output_prob1 | (output_prob2 << bit_precision), sample_uidx * S_CAST(uint64_t, bit_precision_x2), bit_precision_x2, bgen_geno_buf_iter);
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
        const uint16_t* genovec_u16 = DowncastKWToU16(pgv.genovec);
        const uint32_t sample_ct8 = DivUp(sample_ct, 8);
        const uint32_t bit_precision_x4 = bit_precision * 4;
        if (!two_byte_probs) {
          for (uint32_t geno_u16_idx = 0; geno_u16_idx != sample_ct8; ++geno_u16_idx) {
            const uint32_t cur_geno8 = genovec_u16[geno_u16_idx];
            const uint64_t new_bytes = bgen_haploid_hardcall_table8[cur_geno8 & 255] | (S_CAST(uint64_t, bgen_haploid_hardcall_table8[cur_geno8 >> 8]) << bit_precision_x4);
            memcpy(probs_write_citer, &new_bytes, sizeof(int64_t));
            probs_write_citer = &(probs_write_citer[bit_precision]);
          }
        } else {
          const uint32_t backfill = 64 - bit_precision_x4;
          for (uint32_t geno_u16_idx = 0; geno_u16_idx != sample_ct8; ++geno_u16_idx) {
            // this may write 14 bytes past the end
            const uint32_t cur_geno8 = genovec_u16[geno_u16_idx];
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
          for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ++dosage_idx) {
            // multiply by e.g. 255/32768 and round-to-even
            const uintptr_t sample_uidx = BitIter1(pgv.dosage_present, &sample_uidx_base, &dosage_present_bits);
            // note that this needs to be changed to uint64_t if internal
            // dosage representation is no longer 16-bit
            const uint32_t cur_dosage = pgv.dosage_main[dosage_idx];
            const uint32_t output_prob = (cur_dosage * max_output_val + half_int) / kDosageMax;
            CopyBits(output_prob, sample_uidx * S_CAST(uint64_t, bit_precision), bit_precision, bgen_geno_buf_iter);
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
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                uint32_t male_hw = DowncastKWToHW(sex_male_collapsed)[widx];
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
                  const uintptr_t cur_geno = geno_word & 3;
                  const uint32_t cur_male = male_hw & 1;
                  *ploidy_and_missingness_iter++ = (cur_geno == 3) * 128 + 2 - cur_male;
                  if (cur_male) {
                    AppendBits(bit_precision, bgen_haploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                  } else {
                    AppendBits(2 * bit_precision, bgen_diploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
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
              } else {
                AppendU32(0x1000002, &bgen_geno_buf_iter);
                unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
                bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
                *bgen_geno_buf_iter++ = 0;
                *bgen_geno_buf_iter++ = bit_precision;
                for (uint32_t widx = 0; ; ++widx) {
                  if (widx >= sample_ctl2_m1) {
                    if (widx > sample_ctl2_m1) {
                      break;
                    }
                    loop_len = ModNz(sample_ct, kBitsPerWordD2);
                  }
                  uintptr_t geno_word = pgv.genovec[widx];
                  uint32_t female_hw = DowncastKWToHW(sex_female_collapsed)[widx];
                  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
                    const uintptr_t cur_geno = geno_word & 3;
                    const uint32_t cur_female = female_hw & 1;
                    const uint32_t ploidy_and_missingness = (cur_geno == 3) * (128 - cur_female) + 1;
                    *ploidy_and_missingness_iter++ = ploidy_and_missingness;
                    if (ploidy_and_missingness != 128) {
                      AppendBits(bit_precision, bgen_haploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
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
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                uint32_t nonmale_hw = ~DowncastKWToHW(sex_male_collapsed)[widx];
                uint32_t phasepresent_hw = DowncastKWToHW(pgv.phasepresent)[widx];
                uint32_t phaseinfo_hw = DowncastKWToHW(pgv.phaseinfo)[widx];
                // todo: pretty sure this can be sped up by pre-shifting
                // phaseinfo_hw by 3.  not sure what else can be done,
                // though...
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits, geno_word >>= 2, nonmale_hw >>= 1, phasepresent_hw >>= 1, phaseinfo_hw >>= 1) {
                  const uintptr_t cur_geno = geno_word & 3;
                  const uint32_t cur_nonmale = nonmale_hw & 1;
                  if (cur_geno == 3) {
                    *ploidy_and_missingness_iter++ = 129 + cur_nonmale;
                    Append0Bits(bit_precision << cur_nonmale, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                    continue;
                  }
                  const uint32_t cur_phasepresent = phasepresent_hw & 1;
                  if (cur_phasepresent | cur_nonmale) {
                    *ploidy_and_missingness_iter++ = 2;
                    const uint32_t cur_index = cur_geno | (cur_phasepresent * 4) | ((phaseinfo_hw & 1) * 8);
                    AppendBits(bit_precision * 2, bgen_diploid_phased_hardcall_table[cur_index], &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                    continue;
                  }
                  *ploidy_and_missingness_iter++ = 1;
                  AppendBits(bit_precision, bgen_haploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                }
              }
            } else {
              // can't be all-female all-missing since phasepresent_ct > 0
              AppendU32(0x2010002 - cur_y * 0x10000, &bgen_geno_buf_iter);
              unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
              *bgen_geno_buf_iter++ = 1;
              *bgen_geno_buf_iter++ = bit_precision;
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
                  female_hw = DowncastKWToHW(sex_female_collapsed)[widx];
                }
                uint32_t phasepresent_hw = DowncastKWToHW(pgv.phasepresent)[widx];
                uint32_t phaseinfo_hw = DowncastKWToHW(pgv.phaseinfo)[widx];
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits, geno_word >>= 2, phasepresent_hw >>= 1, phaseinfo_hw >>= 1) {
                  const uintptr_t cur_geno = geno_word & 3;
                  if (cur_geno == 3) {
                    // female_hw not incrementally shifted for now
                    if (!(female_hw & (1U << sample_idx_lowbits))) {
                      *ploidy_and_missingness_iter++ = 129;
                      Append0Bits(bit_precision, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                    } else {
                      *ploidy_and_missingness_iter++ = 128;
                    }
                    continue;
                  }
                  const uint32_t cur_phasepresent = phasepresent_hw & 1;
                  if (cur_phasepresent) {
                    *ploidy_and_missingness_iter++ = 2;
                    const uint32_t cur_index = cur_geno | (cur_phasepresent * 4) | ((phaseinfo_hw & 1) * 8);
                    AppendBits(bit_precision * 2, bgen_diploid_phased_hardcall_table[cur_index], &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                    continue;
                  }
                  *ploidy_and_missingness_iter++ = 1;
                  AppendBits(bit_precision, bgen_haploid_basic_table[cur_geno], &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
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
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                const uint32_t male_hw = DowncastKWToHW(sex_male_collapsed)[widx];
                const uint32_t dosage_present_hw = DowncastKWToHW(pgv.dosage_present)[widx];
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
                      Append0Bits(bit_precision << cur_nonmale, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                      continue;
                    }
                    cur_dosage = cur_geno * kDosageMid;
                  }
                  if (cur_male) {
                    *ploidy_and_missingness_iter++ = 1;
                    uint32_t output_prob = (cur_dosage * max_output_val + kDosageMid) / kDosageMax;
                    AppendBits(bit_precision, output_prob, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
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
                  AppendBits(bit_precision * 2, output_prob2 | (output_prob1 << bit_precision), &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
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
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                const uint32_t female_hw = DowncastKWToHW(sex_female_collapsed)[widx];
                const uint32_t dosage_present_hw = DowncastKWToHW(pgv.dosage_present)[widx];
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
                        Append0Bits(bit_precision, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                      } else {
                        *ploidy_and_missingness_iter++ = 128;
                      }
                      continue;
                    }
                    cur_dosage = cur_geno * kDosageMid;
                  }
                  *ploidy_and_missingness_iter++ = 1;
                  uint32_t output_prob = (cur_dosage * max_output_val + kDosageMid) / kDosageMax;
                  AppendBits(bit_precision, output_prob, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
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
              for (uint32_t widx = 0; ; ++widx) {
                if (widx >= sample_ctl2_m1) {
                  if (widx > sample_ctl2_m1) {
                    break;
                  }
                  loop_len = ModNz(sample_ct, kBitsPerWordD2);
                }
                uintptr_t geno_word = pgv.genovec[widx];
                const uint32_t male_hw = DowncastKWToHW(sex_male_collapsed)[widx];
                const uint32_t phaseinfo_hw = DowncastKWToHW(pgv.phaseinfo)[widx];
                const uint32_t dosage_present_hw = DowncastKWToHW(pgv.dosage_present)[widx];
                const uint32_t dphase_present_hw = DowncastKWToHW(pgv.dphase_present)[widx];
                const uint32_t either_phase_present_hw = dphase_present_hw | (DowncastKWToHW(pgv.phasepresent)[widx]);
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
                      Append0Bits(bit_precision << cur_nonmale, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
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
                      AppendBits(bit_precision, output_prob1, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                      continue;
                    }
                    output_prob2 = output_prob1;
                  }
                  *ploidy_and_missingness_iter++ = 2;
                  AppendBits(bit_precision * 2, output_prob1 | (output_prob2 << bit_precision), &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                }
              }
            } else {
              // can't be all-female all-missing since dosage_ct > 0
              AppendU32(0x2010002 - cur_y * 0x10000, &bgen_geno_buf_iter);
              unsigned char* ploidy_and_missingness_iter = bgen_geno_buf_iter;
              bgen_geno_buf_iter = &(bgen_geno_buf_iter[sample_ct]);
              *bgen_geno_buf_iter++ = 1;
              *bgen_geno_buf_iter++ = bit_precision;
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
                  female_hw = DowncastKWToHW(sex_female_collapsed)[widx];
                }
                const uint32_t phaseinfo_hw = DowncastKWToHW(pgv.phaseinfo)[widx];
                const uint32_t dosage_present_hw = DowncastKWToHW(pgv.dosage_present)[widx];
                const uint32_t dphase_present_hw = DowncastKWToHW(pgv.dphase_present)[widx];
                const uint32_t either_phase_present_hw = dphase_present_hw | (DowncastKWToHW(pgv.phasepresent)[widx]);
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
                        Append0Bits(bit_precision, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
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
                    AppendBits(bit_precision * 2, output_prob1 | (output_prob2 << bit_precision), &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                    continue;
                  }
                  *ploidy_and_missingness_iter++ = 1;
                  uint32_t output_prob = (cur_dosage * max_output_val + kDosageMid) / kDosageMax;
                  AppendBits(bit_precision, output_prob, &cur_write_bits, &cur_write_bit_idx, &bgen_geno_buf_iter);
                }
              }
            }
          }
        }
        // bugfix (7 Oct 2024)
        CopyToUnalignedW(bgen_geno_buf_iter, &cur_write_bits);
        bgen_geno_buf_iter = &(bgen_geno_buf_iter[DivUp(cur_write_bit_idx, CHAR_BIT)]);
      }
      const uint32_t uncompressed_bytect = bgen_geno_buf_iter - uncompressed_bgen_geno_buf;
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
        BitvecInvmask(sex_female_collapsed, sample_ctl, missing_acc1);
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
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  {
    VcountIncr4To8(missing_acc4, acc4_vec_ct, missing_acc8);
    VcountIncr8To32(missing_acc8, acc8_vec_ct, missing_acc32);
  }
  while (0) {
  ExportBgen13Thread_err:
    UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    THREAD_BLOCK_FINISH(arg);
    break;
  }
  THREAD_RETURN;
}

// This allocates exported_sample_ids on top.  It can easily be modified to
// export an ID hash table as well (since such a hash table is always
// constructed for the sake of detecting and warning about duplicate
// post-id-paste IDs).
BoolErr ExportIdpaste(const uintptr_t* sample_include, const SampleIdInfo* siip, const char* ftypename, const char* reimport_flagname, uint32_t sample_ct, IdpasteFlags exportf_id_paste, char exportf_id_delim, uintptr_t* max_exported_sample_id_blen_ptr, char** exported_sample_ids_ptr) {
  const uint32_t write_fid = DataFidColIsRequired(sample_include, siip, sample_ct, exportf_id_paste / kfIdpasteMaybefid);
  const char* sample_ids = siip->sample_ids;
  const char* sids = siip->sids;
  const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
  uintptr_t max_sid_blen = siip->max_sid_blen;
  const uint32_t write_iid = (exportf_id_paste / kfIdpasteIid) & 1;
  const uint32_t write_sid = DataSidColIsRequired(sample_include, sids, sample_ct, max_sid_blen, exportf_id_paste / kfIdpasteMaybesid);
  if (write_sid && (!sids)) {
    max_sid_blen = 2;
  }
  uint32_t id_delim_warning = 0;
  char id_delim = exportf_id_delim? exportf_id_delim : '_';
  const uintptr_t max_exported_sample_id_blen = max_sample_id_blen + write_sid * max_sid_blen;
  const uint32_t exported_id_htable_size = GetHtableMinSize(sample_ct);
  uint32_t* new_id_htable;
  // check for duplicates
  if (unlikely(bigstack_alloc_c(sample_ct * max_exported_sample_id_blen, exported_sample_ids_ptr) ||
               bigstack_alloc_u32(exported_id_htable_size, &new_id_htable))) {
    return 1;
  }
  *max_exported_sample_id_blen_ptr = max_exported_sample_id_blen;
  char* exported_sample_ids = *exported_sample_ids_ptr;
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
    if (write_iid) {
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
  if (id_delim_warning && (write_fid + write_iid + write_sid >= 2)) {
    if (reimport_flagname) {
      if (exportf_id_delim) {
        logerrprintfww("Warning: '%c' present in original sample IDs; --%s will not be able to reconstruct them. Consider rerunning with a different --export id-delim= value.\n", exportf_id_delim, reimport_flagname);
      } else {
        logerrprintfww("Warning: '_' present in original sample IDs; --%s will not be able to reconstruct them. Consider rerunning with a suitable --export id-delim= value.\n", reimport_flagname);
      }
    } else {
      if (exportf_id_delim) {
        logerrprintfww("Warning: '%c' present in original sample IDs; this may prevent them from being reconstructed. Consider rerunning with a different --export id-delim= value.\n", exportf_id_delim);
      } else {
        logerrputs("Warning: '_' present in original sample IDs; this may prevent them from being\nreconstructed. Consider rerunning with a suitable --export id-delim= value.\n");
      }
    }
  }
  if (PopulateStrboxHtable(exported_sample_ids, sample_ct, max_exported_sample_id_blen, exported_id_htable_size, new_id_htable)) {
    logerrprintfww("Warning: Duplicate sample ID(s) are being written to --export %s file.\n", ftypename);
  }
  BigstackReset(new_id_htable);
  return 0;
}

PglErr ExportBgen13(const char* outname, const uintptr_t* sample_include, uint32_t* sample_include_cumulative_popcounts, const SampleIdInfo* siip, const uintptr_t* sex_male_collapsed, const uintptr_t* sex_female_collapsed, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, ExportfFlags exportf_flags, uint32_t exportf_bits, IdpasteFlags exportf_id_paste, char exportf_id_delim, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, uint32_t* sample_missing_geno_cts) {
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
    char* chr_buf;
    unsigned char* writebuf;
    if (unlikely(bigstack_alloc_c(max_chr_slen, &chr_buf) ||
                 bigstack_alloc_uc(writebuf_len, &writebuf))) {
      goto ExportBgen13_ret_NOMEM;
    }
    ctx.sex_male_collapsed = sex_male_collapsed;
    ctx.sex_female_collapsed = sex_female_collapsed;

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
    // Save sample IDs unless bgen-omit-sample-id-block specified...
    if (exportf_flags & kfExportfBgenOmitSampleIdBlock) {
      memcpy(writebuf, "\24\0\0", 4);
      writebuf[23] = 0;
    } else {
      char* exported_sample_ids;
      uintptr_t max_exported_sample_id_blen;
      if (unlikely(ExportIdpaste(sample_include, siip, use_zstd_compression? "bgen-1.3" : "bgen-1.2", "bgen", sample_ct, exportf_id_paste, exportf_id_delim, &max_exported_sample_id_blen, &exported_sample_ids))) {
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
        // ...or combined sample ID length is actually greater than 4 GiB, in
        // which case bgen-1.2/1.3 can't actually represent it.
        logerrputs("Warning: Omitting sample ID block from .bgen file since it would overflow (more\nthan 4 GiB).  Consider using shorter IDs.\n");
        memcpy(writebuf, "\24\0\0", 4);
        writebuf[23] = 0; // bugfix (6 Aug 2024): forgot this
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
    }

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

    if (unlikely(bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(ctx.writebufs[0])) ||
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
    ctx.allele_permute = allele_permute;
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
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
          } else {
            PgenErrPrintNV(reterr, ctx.err_info >> 32);
          }
          goto ExportBgen13_ret_1;
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
          if (allele_permute) {
            ref_allele_idx = allele_permute[allele_idx_offset_base];
            alt1_allele_idx = allele_permute[allele_idx_offset_base + 1];
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
    uint32_t* scrambled_missing_cts = DowncastWToU32(&(ctx.missing_acc1[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx != calc_thread_ct; ++tidx) {
      const uint32_t* thread_scrambled_missing_cts = DowncastWToU32(&(ctx.missing_acc1[tidx][acc32_offset]));
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
 ExportBgen13_ret_1:
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
  pgfip->block_base = nullptr;
  return reterr;
}

PglErr ExportOxSample(const char* outname, const uintptr_t* sample_include, const char* sample_ids, const uint32_t* sample_missing_geno_cts, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t variant_ct, uint32_t y_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    logprintfww5("Writing %s ... ", outname);
    fflush(stdout);
    const uint32_t pheno_ctl = BitCtToWordCt(pheno_ct);
    char* writebuf;
    uintptr_t* is_basic_categorical;
    // if phenotype is categorical, and all (non-null) category names are of
    // the form C[positive integer], then it's best to emit the positive
    // integer in the name string instead of the internal index.
    if (unlikely(bigstack_calloc_w(pheno_ctl, &is_basic_categorical) ||
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
    logputs("done.\n");
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

PglErr ExportOxSampleV2(const char* outname, const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uint32_t* sample_missing_geno_cts, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t variant_ct, uint32_t y_ct, IdpasteFlags exportf_id_paste, char exportf_id_delim) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    const char* paternal_ids = piip->parental_id_info.paternal_ids;
    const char* maternal_ids = piip->parental_id_info.maternal_ids;
    const uintptr_t max_sample_id_blen = piip->sii.max_sample_id_blen;
    const uintptr_t max_paternal_id_blen = piip->parental_id_info.max_paternal_id_blen;
    const uintptr_t max_maternal_id_blen = piip->parental_id_info.max_maternal_id_blen;
    uintptr_t writebuf_size = kMaxMediumLine + max_sample_id_blen + 32;
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      const PhenoDtype type_code = pheno_cols[pheno_idx].type_code;
      if (type_code == kPhenoDtypeCc) {
        writebuf_size += 2;
      } else if (type_code == kPhenoDtypeQt) {
        writebuf_size += kMaxDoubleGSlen + 1;
      } else {
        const uint32_t cat_ct = pheno_cols[pheno_idx].nonnull_category_ct + 1;
        const char** cat_names = pheno_cols[pheno_idx].category_names;
        uint32_t max_cat_slen = 2;  // "NA"
        for (uint32_t cat_idx = 1; cat_idx != cat_ct; ++cat_idx) {
          const uint32_t cur_slen = strlen(cat_names[cat_idx]);
          if (cur_slen > max_cat_slen) {
            max_cat_slen = cur_slen;
          }
        }
        writebuf_size += max_cat_slen + 1;
      }
    }
    // possible todo: support column sets here
    const uint32_t write_parents = DataParentalColsAreRequired(sample_include, &(piip->sii), &(piip->parental_id_info), sample_ct, 1);
    if (write_parents) {
      writebuf_size += max_paternal_id_blen + max_maternal_id_blen;
    }

    char* writebuf;
    if (unlikely(bigstack_alloc_c(writebuf_size, &writebuf))) {
      goto ExportOxSampleV2_ret_NOMEM;
    }

    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportOxSampleV2_ret_OPEN_FAIL;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = strcpya_k(writebuf, "ID missing");
    if (write_parents) {
      write_iter = strcpya_k(write_iter, " father mother");
    }
    write_iter = strcpya_k(write_iter, " sex");
    for (uint32_t pheno_idx = 0; pheno_idx != pheno_ct; ++pheno_idx) {
      *write_iter++ = ' ';
      write_iter = strcpya(write_iter, &(pheno_names[pheno_idx * max_pheno_name_blen]));
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto ExportOxSampleV2_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);

    write_iter = strcpya_k(write_iter, "0 0");
    if (write_parents) {
      write_iter = strcpya_k(write_iter, " D D");
    }
    write_iter = strcpya_k(write_iter, " D");
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
        goto ExportOxSampleV2_ret_WRITE_FAIL;
      }
    }
    AppendBinaryEoln(&write_iter);

    char* exported_sample_ids;
    uintptr_t max_exported_sample_id_blen;
    if (unlikely(ExportIdpaste(sample_include, &piip->sii, "sample-v2", "sample", sample_ct, exportf_id_paste, exportf_id_delim, &max_exported_sample_id_blen, &exported_sample_ids))) {
      goto ExportOxSampleV2_ret_NOMEM;
    }

    logprintfww5("Writing %s ... ", outname);
    fflush(stdout);
    const double nonmale_geno_ct_recip = 1.0 / u31tod(variant_ct - y_ct);
    const double male_geno_ct_recip = 1.0 / u31tod(variant_ct);
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      write_iter = strcpyax(write_iter, &(exported_sample_ids[sample_idx * max_exported_sample_id_blen]), ' ');
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      const int32_t cur_missing_geno_ct = sample_missing_geno_cts[sample_idx];
      const uint32_t is_male = IsSet(sex_male, sample_uidx);
      if (is_male) {
        write_iter = dtoa_g(cur_missing_geno_ct * male_geno_ct_recip, write_iter);
      } else {
        write_iter = dtoa_g(cur_missing_geno_ct * nonmale_geno_ct_recip, write_iter);
      }
      *write_iter++ = ' ';
      if (write_parents) {
        write_iter = strcpyax(write_iter, &(paternal_ids[max_paternal_id_blen * sample_uidx]), ' ');
        write_iter = strcpyax(write_iter, &(maternal_ids[max_maternal_id_blen * sample_uidx]), ' ');
      }
      if (IsSet(sex_nm, sample_uidx)) {
        *write_iter++ = '2' - is_male;
      } else {
        write_iter = strcpya_k(write_iter, "NA");
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
            write_iter = strcpya(write_iter, cur_pheno_col->category_names[cur_cat_idx]);
          }
        }
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
        goto ExportOxSampleV2_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
      goto ExportOxSampleV2_ret_WRITE_FAIL;
    }
    logputs("done.\n");
  }
  while (0) {
  ExportOxSampleV2_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportOxSampleV2_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportOxSampleV2_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

const unsigned char kValidVcf43ContigChars[256] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

uint32_t ValidVcfContigName(const char* name_start, const char* name_end, uint32_t v43) {
  // Name already guaranteed to not start with '#'.
  // With either VCF version, we don't allow the contig name to contain a ',',
  // even though the current specification doesn't prohibit it: it messes up
  // parsing of the ##contig header line too often in practice.
  // With v4.3 only, we also prohibit '<', '>', '[', ']', and '*', all
  // characters with ASCII code >126, and first character '='.
  // Colons are *not* prohibited any more, as of 1 Feb 2020.
  const uint32_t slen = name_end - name_start;
  if (unlikely(memchr(name_start, ',', slen))) {
    logerrputs("Error: " PROG_NAME_STR " does not permit contig names in exported VCF files to contain\ncommas, since they are too likely to create parsing problems.\n");
    return 0;
  }
  if (!v43) {
    return 1;
  }
  if (unlikely((name_start[0] == '*') || (name_start[0] == '='))) {
    logerrputs("Error: VCFv4.3 contig names cannot start with '*' or '='.\n");
    return 0;
  }
  for (uint32_t uii = 0; uii != slen; ++uii) {
    unsigned char ucc = name_start[uii];
    if (unlikely(!kValidVcf43ContigChars[ucc])) {
      char* write_iter = strcpya_k(g_logbuf, "Error: Contig name '");
      write_iter = memcpya(write_iter, name_start, slen);
      strcpy_k(write_iter, "' is prohibited by the VCFv4.3 specification ('<', '>', '[', ']', and ASCII codes outside 33..126 not allowed).\n");
      WordWrapB(0);
      logerrputsb();
      return 0;
    }
  }
  return 1;
}

// Similar to WriteRestOfHeaderLine() in plink2_data.
BoolErr BgzfWriteRestOfHeaderLine(const char* hkvline_iter, const char* idval, const char* line_end, uint32_t id_slen, BgzfCompressStream* bgzfp) {
  if (idval > hkvline_iter) {
    // Copy entries before ",ID=".
    const char* idkey_prestart = &(idval[-4]);
    // 1-byte BgzfWrite call is annoying, but this is not the common case.
    if (unlikely(BgzfWrite(",", 1, bgzfp) ||
                 BgzfWrite(hkvline_iter, idkey_prestart - hkvline_iter, bgzfp))) {
      return 1;
    }
  }
  const char* idval_end = &(idval[id_slen]);
  return BgzfWrite(idval_end, line_end - idval_end, bgzfp);
}

// Similar to FixAndWriteContigHeaderLine() in plink2_data.  However, there are
// two changes to the input assumptions:
// * Trailing comma after ID field is assumed to be written.
// * line_end (pointing after the input '\n') replaces closing_gt, since we
//   never need to append an IDX= entry.
// These reduce the number of BgzfWrite() calls.
PglErr BgzfFixAndWriteContigHeaderLine(const char* hkvline_iter, const char* idval, const char* line_end, uint32_t id_slen, uint32_t contig_len, BgzfCompressStream* bgzfp) {
  if (hkvline_iter[-1] == '>') {
    // Only ID field present in input.
    // 20 = strlen("length=") + max 10 bytes for contig_len + '>' + Windows
    //      EOLN
    // need to increase this buffer size if contig_len is changed to uint64_t
    char writebuf[20];
    char* write_iter = strcpya_k(writebuf, "length=");
    write_iter = u32toa_x(contig_len, '>', write_iter);
    AppendBinaryEoln(&write_iter);
    return BgzfWrite(writebuf, write_iter - writebuf, bgzfp)? kPglRetWriteFail : kPglRetSuccess;
  }
  const char* idval_end = &(idval[id_slen]);
  const char* lenvalstr_start;
  uint32_t lenval_slen;
  IntErr interr = HkvlineFindK(hkvline_iter, "length", &lenvalstr_start, &lenval_slen);
  if (interr) {
    if (unlikely(S_CAST(int32_t, interr) == 2)) {
      logerrputs("Error: Malformed ##contig line in .pvar file.\n");
      return kPglRetMalformedInput;
    }
    // Existing header line doesn't contain a length.  Place the length
    // immediately after the ID, and then append the rest of the existing line.
    char writebuf[18];
    char* write_iter = strcpya_k(writebuf, "length=");
    write_iter = u32toa_x(contig_len, ',', write_iter);
    if (unlikely(BgzfWrite(writebuf, write_iter - writebuf, bgzfp))) {
      return kPglRetWriteFail;
    }
    if (idval < hkvline_iter) {
      return BgzfWrite(hkvline_iter, line_end - hkvline_iter, bgzfp)? kPglRetWriteFail : kPglRetSuccess;
    }
    // -3 to exclude "ID="
    const char* idkey_start = &(idval[-3]);
    if (unlikely(BgzfWrite(hkvline_iter, idkey_start - hkvline_iter, bgzfp) ||
                 BgzfWrite(idval_end, line_end - idval_end, bgzfp))) {
      return kPglRetWriteFail;
    }
    return kPglRetSuccess;
  }
  const char* lenvalstr_end = &(lenvalstr_start[lenval_slen]);
  char writebuf[10];
  char* write_iter = u32toa(contig_len, writebuf);
  if (idval < lenvalstr_start) {
    // Usual case: ID was before length in input.
    if (idval < hkvline_iter) {
      if (unlikely(BgzfWrite(hkvline_iter, lenvalstr_start - hkvline_iter, bgzfp))) {
        return kPglRetWriteFail;
      }
    } else {
      const char* idkey_start = &(idval[-3]);
      if (unlikely(BgzfWrite(hkvline_iter, idkey_start - hkvline_iter, bgzfp) ||
                   BgzfWrite(idval_end, lenvalstr_start - idval_end, bgzfp))) {
        return kPglRetWriteFail;
      }
    }
    if (unlikely(BgzfWrite(writebuf, write_iter - writebuf, bgzfp) ||
                 BgzfWrite(lenvalstr_end, line_end - lenvalstr_end, bgzfp))) {
      return kPglRetWriteFail;
    }
    return kPglRetSuccess;
  }
  // ID is later.
  const char* idkey_prestart = &(idval[-4]);
  if (unlikely(BgzfWrite(hkvline_iter, lenvalstr_start - hkvline_iter, bgzfp) ||
               BgzfWrite(writebuf, write_iter - writebuf, bgzfp) ||
               BgzfWrite(lenvalstr_end, idkey_prestart - lenvalstr_end, bgzfp) ||
               BgzfWrite(idval_end, line_end - idval_end, bgzfp))) {
    return kPglRetWriteFail;
  }
  return kPglRetSuccess;
}

uint32_t ValidVcfMetaLine(const char* line_start, uint32_t line_slen, uint32_t v43) {
  const char* first_eq = S_CAST(const char*, memchr(line_start, '=', line_slen));
  if (!first_eq) {
    logerrprintf("Error: Header line in .pvar file does not conform to the VCFv4.%c specification\n(no '=').\n", '2' + v43);
    return 0;
  }
  if (v43 && (first_eq[1] == '<')) {
    // If header line starts with "##key=<", ID is required to conform to the
    // VCFv4.3 (but not 4.2) specification.
    const char* hkvline_iter = &(first_eq[2]);
    const char* idval;
    uint32_t id_slen;
    if (HkvlineIdK(&hkvline_iter, &idval, &id_slen)) {
      logerrprintf("Error: Header line in .pvar file does not conform to the VCFv4.3 specification\n(value starts with '<', but that is not followed by an ID).  This is permitted\nin VCFv4.2, so you may want to export that format instead.\n");
      return 0;
    }
  }
  return 1;
}

uint32_t ValidVcfAlleleCode(const char* allele_code_iter) {
  // returns 1 if probably valid (angle-bracket case is not exhaustively
  // checked), 0 if definitely not
  // TODO: try vectorized loop for long strings (see
  // https://github.com/grailbio/bio/blob/master/biosimd/biosimd_amd64.s ),
  // combined with lookup table for short strings (always use lookup table for
  // first char, then continue using it up to vector boundary).
  uint32_t uii = ctou32(*allele_code_iter);
  if ((uii == '<') || ((uii == '*') && (!allele_code_iter[1]))) {
    return 1;
  }
  do {
    uii -= 64;
    // A = 1, C = 3, G = 7, N = 14, T = 20, so (0x10408a >> ucc) & 1 works as a
    // set membership test
    // (maybe should replace this with another 256-element lookup table...)
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

// Assumes trailing bits of genovec are clear up to word-PAIR boundary.
void UpdateVcfPrevPhased(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, uint32_t sample_ct, uintptr_t* prev_phased) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    const uintptr_t geno_lo = genovec[2 * widx];
    const uintptr_t geno_hi = genovec[2 * widx + 1];
    if (!(geno_lo || geno_hi)) {
      continue;
    }
    const uintptr_t het_word = PackWordToHalfwordMask5555(geno_lo & (~(geno_lo >> 1))) | (S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_hi & (~(geno_hi >> 1)))) << kBitsPerWordD2);
    const uintptr_t phasepresent_word = phasepresent[widx];
    prev_phased[widx] = (prev_phased[widx] & (~het_word)) | phasepresent_word;
  }
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
  // see PrintDdosageDecimal() banker's rounding
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


char* AppendVcfMultiallelicDsForce01(uint32_t allele_ct_m2, uint32_t ds_only, uint32_t hds_force, AlleleCode ac, uint32_t is_hethap, char* write_iter) {
  *write_iter++ = '\t';
  if (!ds_only) {
    write_iter = strcpya_k(write_iter, "0/");
    write_iter = u32toa(ac, write_iter);
    *write_iter++ = ':';
  }
  // DS
  const uint32_t trailing_ac = allele_ct_m2 + 1 - ac;
  *write_iter++ = '0';
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

char* AppendVcfMultiallelicDsForce10Het(uint32_t allele_ct_m2, uint32_t ds_only, uint32_t hds_force, AlleleCode ac0, AlleleCode ac1, uint32_t is_hethap, char* write_iter) {
  *write_iter++ = '\t';
  if (!ds_only) {
    write_iter = u32toa_x(ac0, '/', write_iter);
    write_iter = u32toa_x(ac1, ':', write_iter);
  }
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

char* AppendVcfMultiallelicDsForce10HomDiploid(uint32_t allele_ct_m2, uint32_t ds_only, uint32_t hds_force, AlleleCode ac, char phase_char, char* write_iter) {
  *write_iter++ = '\t';
  if (!ds_only) {
    write_iter = u32toa_x(ac, phase_char, write_iter);
    write_iter = u32toa(ac, write_iter);
    *write_iter++ = ':';
  }
  // ac >= 2 guaranteed
  const uint32_t trailing_ac = allele_ct_m2 + 1 - ac;
  *write_iter++ = '0';
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

char* AppendVcfMultiallelicDsForce10Haploid(uint32_t allele_ct_m2, uint32_t ds_only, uint32_t hds_force, AlleleCode ac, char* write_iter) {
  *write_iter++ = '\t';
  if (!ds_only) {
    write_iter = u32toa(ac, write_iter);
    *write_iter++ = ':';
  }
  // ac >= 2 guaranteed
  const uint32_t trailing_ac = allele_ct_m2 + 1 - ac;
  *write_iter++ = '0';
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

PglErr ExportVcf(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const SampleIdInfo* siip, const uintptr_t* sex_male_collapsed, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const uint32_t* contig_lens, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, ExportfFlags exportf_flags, VcfExportMode vcf_mode, IdpasteFlags exportf_id_paste, char exportf_id_delim, char* xheader, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
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
    uint32_t ds_only = 0;
    uint32_t hds_force = 0;
    if (vcf_mode != kVcfExport0) {
      write_some_dosage = 1;
      if (vcf_mode != kVcfExportGp) {
        write_ds = 1;
        if (vcf_mode == kVcfExportDsForce) {
          ds_force = 1;
        } else if (vcf_mode == kVcfExportDsOnly) {
          // quasi-bugfix (16 Mar 2022): DS-only is now prohibited when partly-
          // haploid or all-haploid variants are present, since without GT
          // there to indicate ploidy, 0..1 vs. 0..2 scaling is not adequately
          // defined.
          if (unlikely((cip->haploid_mask[0] & 1) ||
                       XymtIsNonempty(variant_include, cip, kChrOffsetX) ||
                       XymtIsNonempty(variant_include, cip, kChrOffsetY) ||
                       XymtIsNonempty(variant_include, cip, kChrOffsetMT))) {
            logerrputs("Error: --export vcf-dosage=DS-only mode can only be used with all-diploid data.\n");
            goto ExportVcf_ret_INCONSISTENT_INPUT;
          }
          ds_force = 1;
          ds_only = 1;
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
    //   GT: 4 bytes (unless ds_only)
    //   GP: 3 limited-precision numbers, up to (7 chars + delim) * 3
    //   DS + HDS: also 3 limited precision numbers
    //   DS only: 1 limited-precision number
    // max_allele_ct > 2:
    //   GT: 2 + 2 * UintSlen(max_allele_ct + 1) (unless ds_only)
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
        if (ds_only) {
          output_bytes_per_sample = 8;
        } else {
          output_bytes_per_sample = (write_ds && (!write_hds))? 12 : 28;
        }
      } else {
        output_bytes_per_sample = 4;
      }
    } else {
      const uint32_t max_allele_ct = PgrGetMaxAlleleCt(simple_pgrp);
      if (ds_only) {
        output_bytes_per_sample = 8 * (max_allele_ct - 1);
      } else {
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
    char* write_iter = writebuf;
    const uint32_t v43 = (exportf_flags / kfExportfVcf43) & 1;
    AppendVcfHeaderStart(v43, &write_iter);
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
    uint32_t contig_zero_written = 0;
    uint32_t x_contig_line_written = 0;
    if (xheader) {
      memcpyo_k(writebuf, "##contig=<ID=", 13);
      char* xheader_end = &(xheader[xheader_blen]);
      for (char* line_end = xheader; line_end != xheader_end; ) {
        char* xheader_iter = line_end;
        line_end = AdvPastDelim(xheader_iter, '\n');
        if (StrStartsWithUnsafe(xheader_iter, "##contig=<")) {
          xheader_iter = &(xheader_iter[strlen("##contig=<")]);
          char* idval;
          uint32_t id_slen;
          if (unlikely(HkvlineId(&xheader_iter, &idval, &id_slen))) {
            logerrputs("Error: Malformed ##contig line in .pvar file.\n");
            goto ExportVcf_ret_MALFORMED_INPUT;
          }
          if (xheader_iter[-1] == '>') {
            // regenerate instead of keeping
            continue;
          }
          // if GetChrCodeCounted() is modified to not mutate idval[], xheader
          // can be changed to const char*
          const uint32_t chr_idx = GetChrCodeCounted(cip, id_slen, idval);
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
          char* contig_write_start = &(writebuf[13]);
          write_iter = chrtoa(cip, chr_idx, contig_write_start);
          if ((*contig_write_start == '0') && (write_iter == &(contig_write_start[1]))) {
            // --allow-extra-chr 0 special case
            contig_zero_written = 1;  // technically we write this a bit later
            continue;
          }
          if (unlikely(!ValidVcfContigName(contig_write_start, write_iter, v43))) {
            goto ExportVcf_ret_MALFORMED_INPUT;
          }
          if (contig_lens && contig_lens[chr_idx]) {
            *write_iter++ = ',';
            if (unlikely(BgzfWrite(writebuf, write_iter - writebuf, &bgzf))) {
              goto ExportVcf_ret_WRITE_FAIL;
            }
            reterr = BgzfFixAndWriteContigHeaderLine(xheader_iter, idval, line_end, id_slen, contig_lens[chr_idx], &bgzf);
            if (unlikely(reterr)) {
              goto ExportVcf_ret_1;
            }
          } else {
            if (unlikely(BgzfWrite(writebuf, write_iter - writebuf, &bgzf) ||
                         BgzfWriteRestOfHeaderLine(xheader_iter, idval, line_end, id_slen, &bgzf))) {
              goto ExportVcf_ret_WRITE_FAIL;
            }
          }
        } else if (xheader_iter[1] == '#') {
          // quasi-bugfix (23 Jan 2021): don't export pre-#CHROM header lines
          // starting with a single '#'
          const uint32_t slen = line_end - xheader_iter;
          if (unlikely(!ValidVcfMetaLine(xheader_iter, slen, v43))) {
            goto ExportVcf_ret_INCONSISTENT_INPUT;
          }
          if (unlikely(BgzfWrite(xheader_iter, slen, &bgzf))) {
            goto ExportVcf_ret_WRITE_FAIL;
          }
        }
      }
    }
    write_iter = writebuf;
    // fill in the missing ##contig lines
    if (contig_zero_written) {
      write_iter = strcpya_k(write_iter, "##contig=<ID=0>" EOLN_STR);
    }
    uint32_t chrx_exists = 0;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx != cip->chr_ct; ++chr_fo_idx) {
      if (IsSet(written_contig_header_lines, chr_fo_idx)) {
        continue;
      }
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if ((!IsSet(cip->chr_mask, chr_idx)) || AllBitsAreZero(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1])) {
        continue;
      }
      if ((chr_idx == x_code) || (chr_idx == par1_code) || (chr_idx == par2_code)) {
        chrx_exists = 1;
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
        write_iter = chr_name_write_end;
      } else {
        if (unlikely(!ValidVcfContigName(chr_name_write_start, chr_name_write_end, v43))) {
          goto ExportVcf_ret_MALFORMED_INPUT;
        }
        write_iter = chr_name_write_end;
        if (contig_lens && contig_lens[chr_idx]) {
          write_iter = strcpya_k(write_iter, ",length=");
          write_iter = u32toa(contig_lens[chr_idx], write_iter);
        }
      }
      *write_iter++ = '>';
      AppendBinaryEoln(&write_iter);
      if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
        goto ExportVcf_ret_WRITE_FAIL;
      }
    }
    if (chrx_exists && (!x_contig_line_written)) {
      char* chr_name_write_start = strcpya_k(write_iter, "##contig=<ID=");
      char* chr_name_write_end = chrtoa(cip, x_code, chr_name_write_start);
      write_iter = strcpya_k(chr_name_write_end, ">" EOLN_STR);
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
        logerrputs("Error: Conflicting INFO/PR definitions.  Either fix all REF alleles so that the\n'provisional reference' flag is no longer needed, or remove/rename the other\nuse of the INFO/PR key.\n");
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
    // possible todo: if there is ever a real standard for encoding pedigree
    // information in the VCF header, do it here and make --vcf be able to read
    // it
    if (!ds_only) {
      write_iter = strcpya_k(write_iter, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" EOLN_STR);
    }
    write_iter = strcpya_k(write_iter, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    char* exported_sample_ids;
    uintptr_t max_exported_sample_id_blen;
    if (unlikely(ExportIdpaste(sample_include, siip, "vcf", "vcf", sample_ct, exportf_id_paste, exportf_id_delim, &max_exported_sample_id_blen, &exported_sample_ids))) {
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

    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    PgenVariant pgv;
    if (unlikely(bigstack_alloc_c(max_chr_blen, &chr_buf) ||
                 bigstack_alloc_w(sample_ctl * 2, &pgv.genovec))) {
      goto ExportVcf_ret_NOMEM;
    }
    pgv.genovec[sample_ctl * 2 - 1] = 0;
    pgv.patch_01_set = nullptr;
    pgv.patch_01_vals = nullptr;
    pgv.patch_10_set = nullptr;
    pgv.patch_10_vals = nullptr;
    if (allele_idx_offsets) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &(pgv.patch_01_set)) ||
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
    // Note that prev_phased is NOT reinitialized at the beginning of each
    // chromosome.
    const uint32_t load_dphase = write_hds && (pgfip->gflags & kfPgenGlobalDosagePhasePresent);
    const uint32_t some_phased = (!ds_only) && ((pgfip->gflags & kfPgenGlobalHardcallPhasePresent) || load_dphase);
    uintptr_t* prev_phased = nullptr;
    pgv.phasepresent = nullptr;
    pgv.phaseinfo = nullptr;
    if (some_phased) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &prev_phased) ||
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
      if (unlikely(bigstack_alloc_w(sample_ctl, &(pgv.dosage_present)) ||
                   bigstack_alloc_dosage(sample_ct, &(pgv.dosage_main)))) {
        goto ExportVcf_ret_NOMEM;
      }
      if (load_dphase) {
        if (unlikely(bigstack_alloc_w(sample_ctl, &(pgv.dphase_present)) ||
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
        goto ExportVcf_ret_TSTREAM_FAIL;
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

    // 2-genotype-at-a-time lookup in basic phased case
    uint32_t phased_genotext2[492];
    phased_genotext2[0] = basic_genotext[0];
    phased_genotext2[2] = basic_genotext[1];
    phased_genotext2[4] = basic_genotext[2];
    phased_genotext2[6] = basic_genotext[3];
    phased_genotext2[32] = 0x307c3009;
    phased_genotext2[34] = 0x317c3009;
    phased_genotext2[36] = 0x317c3109;
    phased_genotext2[38] = 0x2e7c2e09;
    phased_genotext2[162] = 0x307c3109;  // 1|0
    InitVcfPhaseLookup4b(phased_genotext2);

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
    // ':' replaced with '\t' if ds_only
    uint16_t ds_inttext[8] = {0x303a, 0x313a, 0x323a, 0x2e3a, 0x303a, 0x2121, 0x313a, 0x2e3a};
    uint16_t ds_inttext4[1024];
    if (ds_only) {
      ds_inttext[0] = 0x3009;
      ds_inttext[1] = 0x3109;
      ds_inttext[2] = 0x3209;
      ds_inttext[3] = 0x2e09;
      ds_inttext[4] = 0x3009;
      ds_inttext[6] = 0x3109;
      ds_inttext[7] = 0x2e09;

      ds_inttext4[0] = 0x3009;
      ds_inttext4[4] = 0x3109;
      ds_inttext4[8] = 0x3209;
      ds_inttext4[12] = 0x2e09;
      InitLookup256x2bx4(ds_inttext4);
    }

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
    const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t is_x = 0;
    uint32_t is_haploid = 0;  // includes chrX and chrY
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
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
      if (allele_permute) {
        ref_allele_idx = allele_permute[allele_idx_offset_base];
        alt1_allele_idx = allele_permute[allele_idx_offset_base + 1];
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
          goto ExportVcf_ret_TSTREAM_FAIL;
        }
      } else {
        if (is_pr) {
          write_iter = strcpya_k(write_iter, "PR");
        } else {
          *write_iter++ = '.';
        }
      }

      // FORMAT
      if (!ds_only) {
        write_iter = strcpya_k(write_iter, "\tGT");
      }

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
            reterr = PgrGet(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec);
          } else {
            reterr = PgrGetD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
          }
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ExportVcf_ret_1;
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
                uint32_t sex_male_hw = is_x * (DowncastKWToHW(sex_male_collapsed)[widx]);
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
          } else if ((!pgv.dosage_ct) && ds_only && (!is_haploid)) {
            write_iter = strcpya_k(write_iter, "\tDS");
            // always 2 bytes wide, exploit that
            GenoarrLookup256x2bx4(pgv.genovec, ds_inttext4, sample_ct, write_iter);
            write_iter = &(write_iter[sample_ct * 2]);
          } else {
            // some dosages present, or {H}DS-force; unphased
            if (write_ds) {
              if (ds_only) {
                write_iter = strcpya_k(write_iter, "\tDS");
              } else {
                write_iter = strcpya_k(write_iter, ":DS");
                if (hds_force) {
                  write_iter = strcpya_k(write_iter, ":HDS");
                }
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
                  dosage_present_hw = DowncastWToHW(pgv.dosage_present)[widx];
                }
                if (ds_only) {
                  for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                    const uint32_t cur_geno = genovec_word & 3;
                    if (dosage_present_hw & 1) {
                      *write_iter++ = '\t';
                      const uint32_t dosage_int = *dosage_main_iter++;
                      write_iter = PrintDiploidVcfDosage(dosage_int, write_ds, write_iter);
                    } else {
                      write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno]), 2);
                    }
                    genovec_word >>= 2;
                    dosage_present_hw >>= 1;
                  }
                } else {
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
              }
            } else {
              // at least partly haploid, unphased
              const char ds_prechar = ds_only? '\t' : ':';
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
                  sex_male_hw = DowncastKWToHW(sex_male_collapsed)[widx];
                }
                if (pgv.dosage_ct) {
                  dosage_present_hw = DowncastWToHW(pgv.dosage_present)[widx];
                }
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  const uint32_t cur_is_male = sex_male_hw & 1;
                  const uint32_t cur_genotext_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                  memcpy(write_iter, &(basic_genotext[cur_geno]), 4);
                  write_iter = &(write_iter[cur_genotext_blen]);
                  if (dosage_present_hw & 1) {
                    *write_iter++ = ds_prechar;
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
            reterr = PgrGetP(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &(pgv.phasepresent_ct));
          } else {
            reterr = PgrGetDp(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
          }
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ExportVcf_ret_1;
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
          // little point in trying to optimize this case, thanks to
          // prev_phased.  Instead, we might want to have a fast path for the
          // all-phased case.
          if (!pgv.phasepresent_ct) {
            ZeroWArr(sample_ctl, pgv.phasepresent);
          }
          if ((!pgv.dosage_ct) && (!ds_force)) {
            if (!is_haploid) {
              ZeroTrailingNyps(sample_ct, pgv.genovec);
              UpdateVcfPrevPhased(pgv.genovec, pgv.phasepresent, sample_ct, prev_phased);
              BitvecAnd(pgv.phasepresent, sample_ctl, pgv.phaseinfo);
              VcfPhaseLookup4b(pgv.genovec, prev_phased, pgv.phaseinfo, phased_genotext2, sample_ct, write_iter);
              write_iter = &(write_iter[sample_ct * 4]);
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
                  is_male_hw = DowncastKWToHW(sex_male_collapsed)[widx];
                }
                uint32_t prev_phased_halfword = DowncastWToHW(prev_phased)[widx];

                const uint32_t phasepresent_hw = DowncastWToHW(pgv.phasepresent)[widx];
                const uint32_t phaseinfo_hw = DowncastWToHW(pgv.phaseinfo)[widx];
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
                DowncastWToHW(prev_phased)[widx] = prev_phased_halfword;
              }
            }
          } else {
            // both dosage (or {H}DS-force) and phase present
            if (write_ds) {
              write_iter = strcpya_k(write_iter, ":DS");
              if (hds_force || pgv.dphase_ct ||
                  (write_hds && pgv.phasepresent_ct && pgv.dosage_ct && (!IntersectionIsEmpty(pgv.phasepresent, pgv.dosage_present, sample_ctl)))) {
                write_iter = strcpya_k(write_iter, ":HDS");
                // dphase_present can be nullptr, so we zero-initialize
                // dphase_present_hw and never refresh it when dphase_ct == 0.
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
                uint32_t prev_phased_halfword = DowncastWToHW(prev_phased)[widx];
                const uint32_t phasepresent_hw = DowncastWToHW(pgv.phasepresent)[widx];
                const uint32_t phaseinfo_hw = DowncastWToHW(pgv.phaseinfo)[widx];
                if (pgv.dosage_ct) {
                  dosage_present_hw = DowncastWToHW(pgv.dosage_present)[widx];
                  if (pgv.dphase_ct) {
                    dphase_present_hw = DowncastWToHW(pgv.dphase_present)[widx];
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
                DowncastWToHW(prev_phased)[widx] = prev_phased_halfword;
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
                  is_male_hw = DowncastKWToHW(sex_male_collapsed)[widx];
                }
                uint32_t prev_phased_halfword = DowncastWToHW(prev_phased)[widx];
                const uint32_t phasepresent_hw = DowncastWToHW(pgv.phasepresent)[widx];
                const uint32_t phaseinfo_hw = DowncastWToHW(pgv.phaseinfo)[widx];
                if (pgv.dosage_ct) {
                  dosage_present_hw = DowncastWToHW(pgv.dosage_present)[widx];
                  if (pgv.dphase_ct) {
                    dphase_present_hw = DowncastWToHW(pgv.dphase_present)[widx];
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
                DowncastWToHW(prev_phased)[widx] = prev_phased_halfword;
              }
            }
          }
        }
      } else {
        // multiallelic cases
        // multiallelic dosage not supported yet
        // don't bother supporting multiallelic GP export when dosages present:
        // no clean way to handle e.g. dosage(A) = dosage(C) = dosage(G) =
        // dosage(T)=0.5
        if (ref_allele_idx || (alt1_allele_idx != 1)) {
          // TODO: imitate --make-pgen allele rotation
          // maybe add a PgrGetM2() function too
          logputs("\n");
          logerrputs("Error: VCF-export multiallelic rotation is under development.\n");
          reterr = kPglRetNotYetSupported;
          goto ExportVcf_ret_1;
        }
        if (!some_phased) {
          reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ExportVcf_ret_1;
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
                  uint32_t multiallelic_hw = (DowncastKWToHW(pgv.patch_01_set)[widx]) | (DowncastKWToHW(pgv.patch_10_set)[widx]);
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
                uint32_t sex_male_hw = is_x * (DowncastKWToHW(sex_male_collapsed)[widx]);
                // no need to separate patch_01 and patch_10 since this
                // information is redundant with genovec_word contents
                uint32_t multiallelic_hw = (DowncastKWToHW(pgv.patch_01_set)[widx]) | (DowncastKWToHW(pgv.patch_10_set)[widx]);
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
            if (ds_only) {
              write_iter = strcpya_k(write_iter, "\tDS");
            } else {
              write_iter = strcpya_k(write_iter, ":DS");
              if (hds_force) {
                write_iter = strcpya_k(write_iter, ":HDS");
              }
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
                uint32_t multiallelic_hw = (DowncastKWToHW(pgv.patch_01_set)[widx]) | (DowncastKWToHW(pgv.patch_10_set)[widx]);
                for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                  const uint32_t cur_geno = genovec_word & 3;
                  if (!(multiallelic_hw & 1)) {
                    if (!ds_only) {
                      write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
                    }
                    // DS
                    write_iter = memcpya_k(write_iter, &(ds_inttext[cur_geno]), 2);
                    if (cur_geno != 3) {
                      // repeat ",0" if nonmissing
                      write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                      if (hds_force) {
                        uint64_t cur_inttext = hds_inttext[cur_geno];
                        const uint32_t hap_blen = hds_inttext_blen[cur_geno + 4];
                        memcpy(write_iter, &cur_inttext, 4);
                        write_iter = &(write_iter[hap_blen]);
                        write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                        cur_inttext = cur_inttext >> (8 * hap_blen);
                        memcpy(write_iter, &cur_inttext, 4);
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
                    write_iter = AppendVcfMultiallelicDsForce01(allele_ct_m2, ds_only, ds_force, ac, 0, write_iter);
                  } else {
                    const AlleleCode ac0 = *patch_10_vals_iter++;
                    const AlleleCode ac1 = *patch_10_vals_iter++;
                    if (ac0 != ac1) {
                      write_iter = AppendVcfMultiallelicDsForce10Het(allele_ct_m2, ds_only, hds_force, ac0, ac1, 0, write_iter);
                    } else {
                      write_iter = AppendVcfMultiallelicDsForce10HomDiploid(allele_ct_m2, ds_only, hds_force, ac0, '/', write_iter);
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
                  sex_male_hw = DowncastKWToHW(sex_male_collapsed)[widx];
                }
                uint32_t multiallelic_hw = (DowncastKWToHW(pgv.patch_01_set)[widx]) | (DowncastKWToHW(pgv.patch_10_set)[widx]);
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
                        uint64_t cur_inttext = hds_inttext[cur_geno];
                        const uint32_t hap_blen = hds_inttext_blen[4 + cur_geno];
                        memcpy(write_iter, &cur_inttext, 4);
                        write_iter = &(write_iter[hap_blen]);
                        if (should_be_diploid) {
                          write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                          cur_inttext = cur_inttext >> (8 * hap_blen);
                          memcpy(write_iter, &cur_inttext, 4);
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
                    write_iter = AppendVcfMultiallelicDsForce01(allele_ct_m2, ds_only, hds_force, ac, !should_be_diploid, write_iter);
                  } else {
                    const AlleleCode ac0 = *patch_10_vals_iter++;
                    const AlleleCode ac1 = *patch_10_vals_iter++;
                    if (ac0 != ac1) {
                      write_iter = AppendVcfMultiallelicDsForce10Het(allele_ct_m2, ds_only, hds_force, ac0, ac1, !should_be_diploid, write_iter);
                    } else if (haploid_genotext_blen[cur_geno + cur_is_male * 4] == 2) {
                      write_iter = AppendVcfMultiallelicDsForce10Haploid(allele_ct_m2, ds_only, hds_force, ac0, write_iter);
                    } else {
                      write_iter = AppendVcfMultiallelicDsForce10HomDiploid(allele_ct_m2, ds_only, hds_force, ac0, '/', write_iter);
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
          reterr = PgrGetMP(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ExportVcf_ret_1;
          }
          if (!pgv.patch_01_ct) {
            ZeroWArr(sample_ctl, pgv.patch_01_set);
          }
          if (!pgv.patch_10_ct) {
            ZeroWArr(sample_ctl, pgv.patch_10_set);
          }
          if (!pgv.phasepresent_ct) {
            ZeroWArr(sample_ctl, pgv.phasepresent);
          }
          const AlleleCode* patch_01_vals_iter = pgv.patch_01_vals;
          const AlleleCode* patch_10_vals_iter = pgv.patch_10_vals;
          if (!ds_force) {
            if (!is_haploid) {
              if (allele_ct <= 10) {
                for (uint32_t widx = 0; ; ++widx) {
                  if (widx >= sample_ctl2_m1) {
                    if (widx > sample_ctl2_m1) {
                      break;
                    }
                    inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
                  }
                  uintptr_t genovec_word = pgv.genovec[widx];
                  const uint32_t multiallelic_hw = (DowncastKWToHW(pgv.patch_01_set)[widx]) | (DowncastKWToHW(pgv.patch_10_set)[widx]);
                  uint32_t prev_phased_halfword = DowncastWToHW(prev_phased)[widx];
                  const uint32_t phasepresent_hw = DowncastWToHW(pgv.phasepresent)[widx];
                  const uint32_t phaseinfo_hw = DowncastWToHW(pgv.phaseinfo)[widx];
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
                    CAppendU32(cur_basic_genotext + 0x4d0000 * ((prev_phased_halfword >> sample_idx_lowbits) & 1), &write_iter);
                    genovec_word >>= 2;
                    cur_shift = cur_shift * 2;
                  }
                  DowncastWToHW(prev_phased)[widx] = prev_phased_halfword;
                }
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
                  const uint32_t multiallelic_hw = (DowncastKWToHW(pgv.patch_01_set)[widx]) | (DowncastKWToHW(pgv.patch_10_set)[widx]);
                  uint32_t prev_phased_halfword = DowncastWToHW(prev_phased)[widx];
                  const uint32_t phasepresent_hw = DowncastWToHW(pgv.phasepresent)[widx];
                  const uint32_t phaseinfo_hw = DowncastWToHW(pgv.phaseinfo)[widx];
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
                  DowncastWToHW(prev_phased)[widx] = prev_phased_halfword;
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
                  is_male_hw = DowncastKWToHW(sex_male_collapsed)[widx];
                }
                const uint32_t multiallelic_hw = (DowncastKWToHW(pgv.patch_01_set)[widx]) | (DowncastKWToHW(pgv.patch_10_set)[widx]);
                uint32_t prev_phased_halfword = DowncastWToHW(prev_phased)[widx];
                const uint32_t phasepresent_hw = DowncastWToHW(pgv.phasepresent)[widx];
                const uint32_t phaseinfo_hw = DowncastWToHW(pgv.phaseinfo)[widx];
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
                DowncastWToHW(prev_phased)[widx] = prev_phased_halfword;
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
                const uint32_t multiallelic_hw = (DowncastKWToHW(pgv.patch_01_set)[widx]) | (DowncastKWToHW(pgv.patch_10_set)[widx]);
                uint32_t prev_phased_halfword = DowncastWToHW(prev_phased)[widx];
                const uint32_t phasepresent_hw = DowncastWToHW(pgv.phasepresent)[widx];
                const uint32_t phaseinfo_hw = DowncastWToHW(pgv.phaseinfo)[widx];
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
                        // bugfix (8 Mar 2020): hap_blen was set incorrectly in
                        // phased case
                        uint32_t hap_blen = 2;
                        if (phasepresent_hw & cur_shift) {
                          hds_inttext_index = 4 + ((phaseinfo_hw >> sample_idx_lowbits) & 1);
                        } else {
                          hap_blen = hds_inttext_blen[hds_inttext_index] / 2;
                        }
                        uint64_t cur_inttext = hds_inttext[hds_inttext_index];
                        memcpy(write_iter, &cur_inttext, 4);
                        write_iter = &(write_iter[hap_blen]);
                        write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                        cur_inttext = cur_inttext >> (8 * hap_blen);
                        memcpy(write_iter, &cur_inttext, 4);
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
                          write_iter = AppendVcfMultiallelicDsForce01(allele_ct_m2, 0, hds_force, ac1, 0, write_iter);
                        } else {
                          write_iter = AppendVcfMultiallelicDsForce10Het(allele_ct_m2, 0, hds_force, ac0, ac1, 0, write_iter);
                        }
                      }
                    } else {
                      write_iter = AppendVcfMultiallelicDsForce10HomDiploid(allele_ct_m2, 0, hds_force, ac0, (prev_phased_halfword & cur_shift)? '|' : '/', write_iter);
                    }
                  }
                  genovec_word >>= 2;
                  cur_shift <<= 1;
                }
                DowncastWToHW(prev_phased)[widx] = prev_phased_halfword;
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
                  is_male_hw = DowncastKWToHW(sex_male_collapsed)[widx];
                }
                const uint32_t multiallelic_hw = (DowncastKWToHW(pgv.patch_01_set)[widx]) | (DowncastKWToHW(pgv.patch_10_set)[widx]);
                uint32_t prev_phased_halfword = DowncastWToHW(prev_phased)[widx];
                const uint32_t phasepresent_hw = DowncastWToHW(pgv.phasepresent)[widx];
                const uint32_t phaseinfo_hw = DowncastWToHW(pgv.phaseinfo)[widx];
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
                          uint32_t hap_blen = 2;
                          if (phasepresent_hw & cur_shift) {
                            hds_inttext_index = 4 + ((phaseinfo_hw >> sample_idx_lowbits) & 1);
                          } else {
                            hap_blen = hds_inttext_blen[hds_inttext_index] / 2;
                          }
                          uint64_t cur_inttext = hds_inttext[hds_inttext_index];
                          memcpy(write_iter, &cur_inttext, 4);
                          write_iter = &(write_iter[hap_blen]);
                          write_iter = u16setsa(write_iter, 0x302c, allele_ct_m2);
                          // Don't render unphased het haploid as diploid here
                          if (should_be_diploid || (phasepresent_hw & cur_shift)) {
                            cur_inttext = cur_inttext >> (8 * hap_blen);
                            memcpy(write_iter, &cur_inttext, 4);
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
                          CAppendU16(0x303a + 128 * cur_geno, &write_iter);
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
                          write_iter = AppendVcfMultiallelicDsForce01(allele_ct_m2, 0, hds_force, ac1, !should_be_diploid, write_iter);
                        } else {
                          write_iter = AppendVcfMultiallelicDsForce10Het(allele_ct_m2, 0, hds_force, ac0, ac1, !should_be_diploid, write_iter);
                        }
                      }
                    } else {
                      if (should_be_diploid) {
                        write_iter = AppendVcfMultiallelicDsForce10HomDiploid(allele_ct_m2, 0, hds_force, ac0, (prev_phased_halfword & cur_shift)? '|' : '/', write_iter);
                      } else {
                        write_iter = AppendVcfMultiallelicDsForce10Haploid(allele_ct_m2, 0, hds_force, ac0, write_iter);
                      }
                    }
                  }
                  genovec_word >>= 2;
                  is_male_hw >>= 1;
                  cur_shift <<= 1;
                }
                DowncastWToHW(prev_phased)[widx] = prev_phased_halfword;
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
        next_print_variant_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
  ExportVcf_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pvar_info_reload, &pvar_reload_txs);
    break;
  ExportVcf_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportVcf_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ExportVcf_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 ExportVcf_ret_1:
  CleanupTextStream2(pvar_info_reload, &pvar_reload_txs, &reterr);
  CleanupBgzfCompressStream(&bgzf, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr AddToFifHtable(unsigned char* arena_bottom, const char* key, uint32_t htable_size, uint32_t key_slen, unsigned char prechar, unsigned char** arena_top_ptr, char** keys, uint32_t* htable, uint32_t* key_ct_ptr, uint32_t* cur_idx_ptr) {
  const uint32_t key_ct = *key_ct_ptr;
  const uint32_t cur_idx = IdHtableAddNnt(key, TO_CONSTCPCONSTP(keys), key_slen, htable_size, key_ct, htable);
  if (cur_idx == UINT32_MAX) {
    if (StoreStringAndPrecharAtEnd(arena_bottom, key, prechar, key_slen, arena_top_ptr, &(keys[key_ct]))) {
      return kPglRetNomem;
    }
    *cur_idx_ptr = key_ct;
    *key_ct_ptr = key_ct + 1;
    return kPglRetSuccess;
  }
  // Key already present.  Only permit this if previous instance was FILTER and
  // this is INFO, or vice versa, or this is a FORMAT key.
  // INFO entries always have prechar bit 0 set and 3 unset, while FILTER
  // entries are the reverse.  FORMAT keys have prechar == 0.
  char* existing_key = keys[cur_idx];
  const unsigned char new_prechar = S_CAST(unsigned char, existing_key[-1]) ^ prechar;
  if (unlikely(prechar && ((new_prechar & 9) != 9))) {
    char* write_iter = strcpya_k(g_logbuf, "Error: Multiple ");
    if (prechar & 8) {
      write_iter = strcpya_k(write_iter, "FILTER");
    } else {
      write_iter = strcpya_k(write_iter, "INFO");
    }
    *write_iter++ = ':';
    write_iter = memcpya(write_iter, key, key_slen);
    strcpy_k(write_iter, " header lines.\n");
    WordWrapB(0);
    logerrputsb();
    return kPglRetMalformedInput;
  }
  existing_key[-1] = new_prechar;
  // bugfix (3 Jan 2021): forgot to return this
  *cur_idx_ptr = cur_idx;
  return kPglRetSuccess;
}

// Assumes ii is not in 0x80000000..0x80000007.
char* AppendBcfTypedInt(int32_t ii, char* write_iter) {
  // 0x80..0x87 reserved
  if ((ii >= -120) && (ii < 128)) {
    *write_iter++ = 0x11;
    *write_iter++ = ii;
    return write_iter;
  }
  if ((ii >= -32760) && (ii < 32768)) {
    *write_iter++ = 0x12;
    int16_t sii = ii;
    memcpy(write_iter, &sii, sizeof(int16_t));
    return &(write_iter[2]);
  }
  *write_iter++ = 0x13;
  memcpy(write_iter, &ii, sizeof(int32_t));
  return &(write_iter[4]);
}

char* AppendBcfVecType(uint32_t type_code, uint32_t len, char* write_iter) {
  if (len < 15) {
    *write_iter++ = (len << 4) | type_code;
    return write_iter;
  }
  *write_iter++ = 0xf0 | type_code;
  return AppendBcfTypedInt(len, write_iter);
}

char* AppendBcfCountedString(const char* ss, uint32_t slen, char* write_iter) {
  write_iter = AppendBcfVecType(7, slen, write_iter);
  write_iter = memcpya(write_iter, ss, slen);
  return write_iter;
}

static inline char* AppendBcfString(const char* ss, char* write_iter) {
  const uint32_t slen = strlen(ss);
  return AppendBcfCountedString(ss, slen, write_iter);
}


// TODO: SSE4-optimized versions of many of these functions.  See
// InverseMovemaskFF().
// For each male non-het, replace the second GT byte with 0x81 (END_OF_VECTOR).
void FixBcfMaleXGtPloidy(const uintptr_t* genovec, const uintptr_t* sex_male, uint32_t sample_ct, char* __restrict gt_start) {
  const Halfword* sex_male_alias = DowncastKWToHW(sex_male);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  char* gt_offset1 = &(gt_start[1]);
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uint32_t cur_male_hw = sex_male_alias[widx];
    const uintptr_t geno_word = genovec[widx];
    const uint32_t geno_nonhet_hw = PackWordToHalfwordMask5555((~geno_word) | (geno_word >> 1));
    uint32_t male_nonhet_hw = cur_male_hw & geno_nonhet_hw;
    if (male_nonhet_hw) {
      char* cur_gt_offset1 = &(gt_offset1[widx * kBitsPerWord]);
      do {
        const uint32_t sample_idx_lowbits = ctzu32(male_nonhet_hw);
        cur_gt_offset1[sample_idx_lowbits * 2] = 0x81;
        male_nonhet_hw &= male_nonhet_hw - 1;
      } while (male_nonhet_hw);
    }
  }
}

// For each female missing call, replace all byte(s) with 0x81 END_OF_VECTOR.
void FixBcfFemaleYGtPloidy2(const uintptr_t* genovec, const uintptr_t* sex_female, uint32_t sample_ct, char* __restrict gt_start_orig) {
  const Halfword* sex_female_alias = DowncastKWToHW(sex_female);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint16_t end_of_vector = 0x8181;
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uint32_t cur_female_hw = sex_female_alias[widx];
    const uint32_t geno_missing_hw = Pack11ToHalfword(genovec[widx]);
    uint32_t female_missing_hw = cur_female_hw & geno_missing_hw;
    if (female_missing_hw) {
      unsigned char* cur_gt_start = CToUc(&(gt_start_orig[widx * kBitsPerWord]));
      do {
        const uint32_t sample_idx_lowbits = ctzu32(female_missing_hw);
        CopyToUnalignedOffsetU16(cur_gt_start, &end_of_vector, sample_idx_lowbits);
        female_missing_hw &= female_missing_hw - 1;
      } while (female_missing_hw);
    }
  }
}

void FixBcfFemaleYGtPloidy1(const uintptr_t* genovec, const uintptr_t* sex_female, uint32_t sample_ct, char* __restrict gt_start) {
  const Halfword* sex_female_alias = DowncastKWToHW(sex_female);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uint32_t cur_female_hw = sex_female_alias[widx];
    const uint32_t geno_missing_hw = Pack11ToHalfword(genovec[widx]);
    uint32_t female_missing_hw = cur_female_hw & geno_missing_hw;
    if (female_missing_hw) {
      char* cur_gt_start = &(gt_start[widx * kBitsPerWordD2]);
      do {
        const uint32_t sample_idx_lowbits = ctzu32(female_missing_hw);
        cur_gt_start[sample_idx_lowbits] = 0x81;
        female_missing_hw &= female_missing_hw - 1;
      } while (female_missing_hw);
    }
  }
}

void FillBcfDsPloidy1(const PgenVariant* pgvp, const uintptr_t* __restrict sex_female, const uint32_t* __restrict ds_y_genobytes2, const uint32_t* __restrict ds_haploid_genobytes4, uint32_t sample_ct, uint32_t forced, uint32_t is_y, void* __restrict ds_start_orig) {
  const uintptr_t* genovec = pgvp->genovec;
  // ds_start not required to be 4-byte aligned
  if (forced) {
    if (is_y) {
      // most dosages 0..1, het-hap = 0..2, female missing -> END_OF_VECTOR
      GenoarrSexLookup4b(genovec, sex_female, ds_y_genobytes2, sample_ct, ds_start_orig);
    } else {
      // most dosages 0..1, het-hap = 0..2
      GenoarrLookup256x4bx4(genovec, ds_haploid_genobytes4, sample_ct, ds_start_orig);
    }
  }
  const uint32_t dosage_ct = pgvp->dosage_ct;
  if (!dosage_ct) {
    return;
  }
  const Halfword* dosage_present_alias = DowncastKWToHW(pgvp->dosage_present);
  const Dosage* dosage_main_stop = &(pgvp->dosage_main[dosage_ct]);
  unsigned char* ds_start_uc = S_CAST(unsigned char*, ds_start_orig);
  uint32_t widx = UINT32_MAX;  // deliberate overflow
  for (const Dosage* dosage_main_iter = pgvp->dosage_main; dosage_main_iter != dosage_main_stop; ) {
    uint32_t dosage_present_hw;
    do {
      dosage_present_hw = dosage_present_alias[++widx];
    } while (!dosage_present_hw);
    const uintptr_t geno_word = genovec[widx];
    const uint32_t cur_hets = PackWordToHalfwordMask5555(geno_word & (~(geno_word >> 1)));
    unsigned char* cur_ds = &(ds_start_uc[widx * (kBitsPerWord * 2)]);
    do {
      const uint32_t sample_idx_lowbits = ctzu32(dosage_present_hw);
      const uint32_t dosage_int = *dosage_main_iter++;
      const uint32_t scaled_dosage_int = dosage_int << ((cur_hets >> sample_idx_lowbits) & 1);  // het-haps on 0..2 scale
      const float dosage_f = S_CAST(float, scaled_dosage_int) * S_CAST(float, kRecipDosageMax);
      CopyToUnalignedOffsetF(cur_ds, &dosage_f, sample_idx_lowbits);
      dosage_present_hw &= dosage_present_hw - 1;
    } while (dosage_present_hw);
  }
}

// Biallelic.
void FillBcfDs(const PgenVariant* pgvp, const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_female, const uint32_t* __restrict ds_genobytes4, const uint32_t* __restrict ds_x_genobytes2, const uint32_t* __restrict ds_y_genobytes2, const uint32_t* __restrict ds_haploid_genobytes4, uint32_t sample_ct, uint32_t ds_force, uint32_t is_x, uint32_t is_y, uint32_t is_haploid, void* __restrict ds_start_orig) {
  unsigned char* ds_start_uc = S_CAST(unsigned char*, ds_start_orig);
  if (!ds_force) {
    // ds_start not required to be 4-byte aligned
    const uint32_t end_of_vector = 0x7f800002;
    for (uint32_t uii = 0; uii != sample_ct; ++uii) {
      // END_OF_VECTOR is more appropriate than MISSING here
      // see "Vectors of mixed length" section in the spec
      CopyToUnalignedOffsetU32(ds_start_uc, &end_of_vector, uii);
    }
  }
  const uintptr_t* genovec = pgvp->genovec;
  const uintptr_t* dosage_present = pgvp->dosage_present;
  const Dosage* dosage_main = pgvp->dosage_main;
  if (!is_haploid) {
    if (ds_force) {
      GenoarrLookup256x4bx4(genovec, ds_genobytes4, sample_ct, ds_start_uc);
    }
    // these loops are optimized for sparse dosages.  possible todo: branch on
    // dosage_ct/sample_ct ratio.
    uint32_t widx = UINT32_MAX;  // deliberate overflow
    for (uint32_t dosage_idx = 0; dosage_idx != pgvp->dosage_ct; ) {
      uintptr_t dosage_present_word;
      do {
        dosage_present_word = dosage_present[++widx];
      } while (!dosage_present_word);
      unsigned char* cur_ds = &(ds_start_uc[widx * (kBitsPerWord * k1LU * 4)]);
      do {
        const uint32_t sample_idx_lowbits = ctzw(dosage_present_word);
        const uint32_t dosage_int = dosage_main[dosage_idx++];
        const float dosage_f = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMid);
        CopyToUnalignedOffsetF(cur_ds, &dosage_f, sample_idx_lowbits);
        dosage_present_word &= dosage_present_word - 1;
      } while (dosage_present_word);
    }
    return;
  }
  if (!is_x) {
    FillBcfDsPloidy1(pgvp, sex_female, ds_y_genobytes2, ds_haploid_genobytes4, sample_ct, ds_force, is_y, ds_start_uc);
    return;
  }
  if (ds_force) {
    // most male dosages 0..1, het-hap = 0..2
    GenoarrSexLookup4b(genovec, sex_male, ds_x_genobytes2, sample_ct, ds_start_uc);
  }
  const uint32_t dosage_ct = pgvp->dosage_ct;
  if (!dosage_ct) {
    return;
  }
  const Halfword* sex_male_alias = DowncastKWToHW(sex_male);
  const Halfword* dosage_present_alias = DowncastKWToHW(dosage_present);
  const Dosage* dosage_main_stop = &(dosage_main[dosage_ct]);
  uint32_t widx = UINT32_MAX;  // deliberate overflow
  for (const Dosage* dosage_main_iter = dosage_main; dosage_main_iter != dosage_main_stop; ) {
    uint32_t dosage_present_hw;
    do {
      dosage_present_hw = dosage_present_alias[++widx];
    } while (!dosage_present_hw);
    const uintptr_t geno_word = genovec[widx];
    const uint32_t cur_hets = PackWordToHalfwordMask5555(geno_word & (~(geno_word >> 1)));
    const uint32_t male_hw = sex_male_alias[widx];
    const uint32_t ploidy_2_hw = cur_hets | (~male_hw);
    unsigned char* cur_ds = &(ds_start_uc[widx * (kBitsPerWord * 2)]);
    do {
      const uint32_t sample_idx_lowbits = ctzu32(dosage_present_hw);
      const uint32_t is_ploidy_2 = 1 & (ploidy_2_hw >> sample_idx_lowbits);
      const uint32_t scaled_dosage_int = (*dosage_main_iter++) << is_ploidy_2;
      const float dosage_f = S_CAST(float, scaled_dosage_int) * S_CAST(float, kRecipDosageMax);
      CopyToUnalignedOffsetF(cur_ds, &dosage_f, sample_idx_lowbits);
      dosage_present_hw &= dosage_present_hw - 1;
    } while (dosage_present_hw);
  }
}

void FillBcfPhasedGt(const uintptr_t* __restrict genovec, const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict phaseinfo, const uint16_t* __restrict phased_genobytes2, uint32_t sample_ct, uintptr_t* __restrict prev_phased, char* __restrict gt_start) {
  // Assumes trailing bits of genovec zeroed up to word-PAIR, and phaseinfo has
  // been anded with phasepresent.
  // No need to handle chrX separately since FixBcfMaleXGtPloidy() converts all
  // male phased-hom calls generated by this function into ploidy-1.
  // Similarly, no need to handle chrY separately.
  UpdateVcfPrevPhased(genovec, phasepresent, sample_ct, prev_phased);
  VcfPhaseLookup2b(genovec, prev_phased, phaseinfo, phased_genobytes2, sample_ct, gt_start);
}

void FixBcfMaleXHdsPloidy(const uintptr_t* __restrict phasepresent, const uintptr_t* __restrict dphase_present, const uintptr_t* __restrict sex_male, uint32_t sample_ct, void* __restrict hds_start) {
  // Keep male HDS encoded as ploidy 2 when either phasepresent or
  // dphase_present bit set.
  // Otherwise replace the second value with 0x7f800002.
  // This does some unnecessary work when the missingness rate is high, but we
  // can live with that.
  unsigned char* hds_offset1_uc = &(S_CAST(unsigned char*, hds_start)[4]);
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t end_of_vector = 0x7f800002;
  if (!dphase_present) {
    for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
      const uintptr_t male_word = sex_male[widx];
      const uintptr_t phasepresent_word = phasepresent[widx];
      uintptr_t male_unphased_word = male_word & (~phasepresent_word);
      if (male_unphased_word) {
        unsigned char* cur_hds_offset1 = &(hds_offset1_uc[widx * (kBitsPerWord * k1LU * 8)]);
        do {
          const uint32_t sample_idx_lowbits = ctzw(male_unphased_word);
          CopyToUnalignedOffsetU32(cur_hds_offset1, &end_of_vector, sample_idx_lowbits * 2);
          male_unphased_word &= male_unphased_word - 1;
        } while (male_unphased_word);
      }
    }
    return;
  }
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    const uintptr_t male_word = sex_male[widx];
    const uintptr_t phasepresent_word = phasepresent[widx] | dphase_present[widx];
    uintptr_t male_unphased_word = male_word & (~phasepresent_word);
    if (male_unphased_word) {
      unsigned char* cur_hds_offset1 = &(hds_offset1_uc[widx * (kBitsPerWord * k1LU * 8)]);
      do {
        const uint32_t sample_idx_lowbits = ctzw(male_unphased_word);
        CopyToUnalignedOffsetU32(cur_hds_offset1, &end_of_vector, sample_idx_lowbits * 2);
        male_unphased_word &= male_unphased_word - 1;
      } while (male_unphased_word);
    }
  }
}

// Replace first HDS entry for each female missing call with 0x7f800002
// END_OF_VECTOR.
void FixBcfFemaleYHdsPloidy2(const uintptr_t* __restrict genovec, const uintptr_t* __restrict dosage_present, const uintptr_t* __restrict sex_female, uint32_t sample_ct, void* __restrict hds_start) {
  unsigned char* hds_uc = S_CAST(unsigned char*, hds_start);
  const Halfword* sex_female_alias = DowncastKWToHW(sex_female);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t end_of_vector = 0x7f800002;
  if (!dosage_present) {
    for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
      const uint32_t cur_female_hw = sex_female_alias[widx];
      const uint32_t geno_missing_hw = Pack11ToHalfword(genovec[widx]);
      uint32_t female_missing_hw = cur_female_hw & geno_missing_hw;
      if (female_missing_hw) {
        unsigned char* cur_hds = &(hds_uc[widx * (kBitsPerWord * k1LU * 4)]);
        do {
          const uint32_t sample_idx_lowbits = ctzu32(female_missing_hw);
          CopyToUnalignedOffsetU32(cur_hds, &end_of_vector, sample_idx_lowbits * 2);
          female_missing_hw &= female_missing_hw - 1;
        } while (female_missing_hw);
      }
    }
    return;
  }
  const Halfword* dosage_present_alias = DowncastKWToHW(dosage_present);
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uint32_t cur_female_hw = sex_female_alias[widx];
    const uint32_t geno_missing_hw = Pack11ToHalfword(genovec[widx]);
    const uint32_t dosage_nm_hw = dosage_present_alias[widx];
    uint32_t female_missing_hw = cur_female_hw & geno_missing_hw & (~dosage_nm_hw);
    if (female_missing_hw) {
      unsigned char* cur_hds = &(hds_uc[widx * (kBitsPerWord * k1LU * 4)]);
      do {
        const uint32_t sample_idx_lowbits = ctzu32(female_missing_hw);
        CopyToUnalignedOffsetU32(cur_hds, &end_of_vector, sample_idx_lowbits * 2);
        female_missing_hw &= female_missing_hw - 1;
      } while (female_missing_hw);
    }
  }
}

void FixBcfFemaleY64allelicGtPloidy2(const uintptr_t* genovec, const uintptr_t* sex_female, uint32_t sample_ct, char* __restrict gt_start_orig) {
  const Halfword* sex_female_alias = DowncastKWToHW(sex_female);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint32_t end_of_vector = 0x80018001U;
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uint32_t cur_female_hw = sex_female_alias[widx];
    const uint32_t geno_missing_hw = Pack11ToHalfword(genovec[widx]);
    uint32_t female_missing_hw = cur_female_hw & geno_missing_hw;
    if (female_missing_hw) {
      unsigned char* cur_gt_start = CToUc(&(gt_start_orig[widx * (kBitsPerWord * 2)]));
      do {
        const uint32_t sample_idx_lowbits = ctzu32(female_missing_hw);
        CopyToUnalignedOffsetU32(cur_gt_start, &end_of_vector, sample_idx_lowbits);
        female_missing_hw &= female_missing_hw - 1;
      } while (female_missing_hw);
    }
  }
}

void FixBcfFemaleY64allelicGtPloidy1(const uintptr_t* genovec, const uintptr_t* sex_female, uint32_t sample_ct, char* __restrict gt_start_orig) {
  const Halfword* sex_female_alias = DowncastKWToHW(sex_female);
  const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
  const uint16_t end_of_vector = 0x8001;
  for (uint32_t widx = 0; widx != sample_ctl2; ++widx) {
    const uint32_t cur_female_hw = sex_female_alias[widx];
    const uint32_t geno_missing_hw = Pack11ToHalfword(genovec[widx]);
    uint32_t female_missing_hw = cur_female_hw & geno_missing_hw;
    if (female_missing_hw) {
      unsigned char* cur_gt_start = CToUc(&(gt_start_orig[widx * kBitsPerWord]));
      do {
        const uint32_t sample_idx_lowbits = ctzu32(female_missing_hw);
        CopyToUnalignedOffsetU16(cur_gt_start, &end_of_vector, sample_idx_lowbits);
        female_missing_hw &= female_missing_hw - 1;
      } while (female_missing_hw);
    }
  }
}

// Just like DS, except with no het-hap special case.
void FillBcfHdsPloidy1(const PgenVariant* pgvp, const uintptr_t* __restrict sex_female, const uint32_t* __restrict hds_y_genobytes2, const uint32_t* __restrict hds_haploid_genobytes4, uint32_t sample_ct, uint32_t forced, uint32_t is_y, void* __restrict hds_start_orig) {
  const uintptr_t* genovec = pgvp->genovec;
  if (forced) {
    if (is_y) {
      // most dosages 0..1, female missing -> END_OF_VECTOR
      GenoarrSexLookup4b(genovec, sex_female, hds_y_genobytes2, sample_ct, hds_start_orig);
    } else {
      // all dosages 0..1
      GenoarrLookup256x4bx4(genovec, hds_haploid_genobytes4, sample_ct, hds_start_orig);
    }
  }
  const uint32_t dosage_ct = pgvp->dosage_ct;
  if (!dosage_ct) {
    return;
  }
  const uintptr_t* dosage_present = pgvp->dosage_present;
  const Dosage* dosage_main_stop = &(pgvp->dosage_main[dosage_ct]);
  unsigned char* hds_start_uc = S_CAST(unsigned char*, hds_start_orig);
  uint32_t widx = UINT32_MAX;  // deliberate overflow
  for (const Dosage* dosage_main_iter = pgvp->dosage_main; dosage_main_iter != dosage_main_stop; ) {
    uintptr_t dosage_present_word;
    do {
      dosage_present_word = dosage_present[++widx];
    } while (!dosage_present_word);
    unsigned char* cur_hds_start = &(hds_start_uc[widx * (kBitsPerWord * k1LU * 4)]);
    do {
      const uint32_t sample_idx_lowbits = ctzw(dosage_present_word);
      const uint32_t dosage_int = *dosage_main_iter++;
      const float dosage_f = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMax);
      CopyToUnalignedOffsetF(cur_hds_start, &dosage_f, sample_idx_lowbits);
      dosage_present_word &= dosage_present_word - 1;
    } while (dosage_present_word);
  }
}

void FillBcfGpPloidy2(const PgenVariant* pgvp, const uintptr_t* __restrict sex_male, uint32_t sample_ct, uint32_t is_x, void* __restrict gp_start_orig) {
  unsigned char* gp_start_uc = S_CAST(unsigned char*, gp_start_orig);
  const uint32_t eov_bits = 0x7f800002;
  for (uint32_t uii = 0; uii != 3 * sample_ct; ++uii) {
    CopyToUnalignedOffsetU32(gp_start_uc, &eov_bits, uii);
  }
  const uint32_t dosage_ct = pgvp->dosage_ct;
  if (!dosage_ct) {
    return;
  }
  // dosage_main, etc. now guaranteed to be non-null
  const uintptr_t* dosage_present = pgvp->dosage_present;
  const Dosage* dosage_main_stop = &(pgvp->dosage_main[dosage_ct]);
  uint32_t widx = UINT32_MAX;  // deliberate overflow
  const float zero = 0.0;
  if (!is_x) {
    for (const Dosage* dosage_main_iter = pgvp->dosage_main; dosage_main_iter != dosage_main_stop; ) {
      uintptr_t dosage_present_word;
      do {
        dosage_present_word = dosage_present[++widx];
      } while (!dosage_present_word);
      unsigned char* gp_word_uc = &(gp_start_uc[widx * (kBitsPerWord * k1LU * 12)]);
      do {
        const uint32_t sample_idx_lowbits = ctzw(dosage_present_word);
        unsigned char* cur_gp = &(gp_word_uc[12 * sample_idx_lowbits]);
        const uint32_t dosage_int = *dosage_main_iter++;
        if (dosage_int < kDosageMid) {
          const float dosage0 = S_CAST(float, kDosageMid - dosage_int) * S_CAST(float, kRecipDosageMid);
          const float dosage1 = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMid);
          CopyToUnalignedF(cur_gp, &dosage0);
          CopyToUnalignedOffsetF(cur_gp, &dosage1, 1);
          CopyToUnalignedOffsetF(cur_gp, &zero, 2);
        } else {
          const float dosage1 = S_CAST(float, kDosageMax - dosage_int) * S_CAST(float, kRecipDosageMid);
          const float dosage2 = S_CAST(float, dosage_int - kDosageMid) * S_CAST(float, kRecipDosageMid);
          CopyToUnalignedF(cur_gp, &zero);
          CopyToUnalignedOffsetF(cur_gp, &dosage1, 1);
          CopyToUnalignedOffsetF(cur_gp, &dosage2, 2);
        }
        dosage_present_word &= dosage_present_word - 1;
      } while (dosage_present_word);
    }
    return;
  }
  for (const Dosage* dosage_main_iter = pgvp->dosage_main; dosage_main_iter != dosage_main_stop; ) {
    uintptr_t dosage_present_word;
    do {
      dosage_present_word = dosage_present[++widx];
    } while (!dosage_present_word);
    const uintptr_t male_word = sex_male[widx];
    unsigned char* gp_word_uc = &(gp_start_uc[widx * (kBitsPerWord * k1LU * 12)]);
    do {
      const uint32_t sample_idx_lowbits = ctzw(dosage_present_word);
      unsigned char* cur_gp = &(gp_word_uc[12 * sample_idx_lowbits]);
      const uint32_t dosage_int = *dosage_main_iter++;
      if ((male_word >> sample_idx_lowbits) & 1) {
        const float ref_prob = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMax);
        const float dosage1 = S_CAST(float, 1.0) - ref_prob;
        CopyToUnalignedF(cur_gp, &ref_prob);
        CopyToUnalignedOffsetF(cur_gp, &dosage1, 1);
        CopyToUnalignedOffsetU32(cur_gp, &eov_bits, 2);
      } else {
        if (dosage_int < kDosageMid) {
          const float het_prob = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMid);
          const float dosage0 = S_CAST(float, 1.0) - het_prob;
          CopyToUnalignedF(cur_gp, &dosage0);
          CopyToUnalignedOffsetF(cur_gp, &het_prob, 1);
          CopyToUnalignedOffsetF(cur_gp, &zero, 2);
        } else {
          const float het_prob = S_CAST(float, kDosageMax - dosage_int) * S_CAST(float, kRecipDosageMid);
          const float dosage2 = S_CAST(float, 1.0) - het_prob;
          CopyToUnalignedF(cur_gp, &zero);
          CopyToUnalignedOffsetF(cur_gp, &het_prob, 1);
          CopyToUnalignedOffsetF(cur_gp, &dosage2, 2);
        }
      }
      dosage_present_word &= dosage_present_word - 1;
    } while (dosage_present_word);
  }
}

void FillBcfGpPloidy1(const PgenVariant* pgvp, uint32_t sample_ct, void* __restrict gp_start_orig) {
  unsigned char* gp_start_uc = S_CAST(unsigned char*, gp_start_orig);
  const uint64_t end_of_vector = 0x7f8000027f800002LLU;
  for (uint32_t uii = 0; uii != sample_ct; ++uii) {
    CopyToUnalignedOffsetU64(gp_start_uc, &end_of_vector, uii);
  }
  const uint32_t dosage_ct = pgvp->dosage_ct;
  if (!dosage_ct) {
    return;
  }
  const uintptr_t* dosage_present = pgvp->dosage_present;
  const Dosage* dosage_main_stop = &(pgvp->dosage_main[dosage_ct]);
  uint32_t widx = UINT32_MAX;  // deliberate overflow
  for (const Dosage* dosage_main_iter = pgvp->dosage_main; dosage_main_iter != dosage_main_stop; ) {
    uintptr_t dosage_present_word;
    do {
      dosage_present_word = dosage_present[++widx];
    } while (!dosage_present_word);
    unsigned char* gp_word_uc = &(gp_start_uc[widx * (kBitsPerWord * k1LU * 8)]);
    do {
      const uint32_t sample_idx_lowbits = ctzw(dosage_present_word);
      unsigned char* cur_gp = &(gp_word_uc[8 * sample_idx_lowbits]);
      const uint32_t dosage_int = *dosage_main_iter++;
      const float alt_prob = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMax);
      const float dosage0 = S_CAST(float, 1.0) - alt_prob;
      CopyToUnalignedF(cur_gp, &dosage0);
      CopyToUnalignedOffsetF(cur_gp, &alt_prob, 1);
      dosage_present_word &= dosage_present_word - 1;
    } while (dosage_present_word);
  }
}


// todo: benchmark this against the biallelic-specialized function, and delete
// the latter if it isn't significantly faster.
void FixBcfMaleX3allelicGtPloidy(const uintptr_t* __restrict sex_male, uint32_t sample_ct, char* __restrict gt_start) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    uintptr_t male_word = sex_male[widx];
    if (male_word) {
      char* cur_gt = &(gt_start[widx * kBitsPerWord * 2]);
      do {
        const uint32_t sample_idx_lowbits = ctzw(male_word);
        const uint32_t gt_low = ctou32(cur_gt[2 * sample_idx_lowbits]);
        // clear phase bit
        const uint32_t gt_high = ctou32(cur_gt[2 * sample_idx_lowbits + 1]) & 0xfe;
        if (gt_low == gt_high) {
          cur_gt[2 * sample_idx_lowbits + 1] = 0x81;
        }
        male_word &= male_word - 1;
      } while (male_word);
    }
  }
}

void FixBcfMaleX64allelicGtPloidy(const uintptr_t* __restrict sex_male, uint32_t sample_ct, char* __restrict gt_start) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint16_t end_of_vector = 0x8001;
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    uintptr_t male_word = sex_male[widx];
    if (male_word) {
      unsigned char* cur_gt = CToUc(&(gt_start[widx * (kBitsPerWord * k1LU * 4)]));
      do {
        const uint32_t sample_idx_lowbits = ctzw(male_word);
        const uint32_t gt_low = CopyFromUnalignedOffsetU16ZX(cur_gt, 2 * sample_idx_lowbits);
        const uint32_t gt_high = CopyFromUnalignedOffsetU16ZX(cur_gt, 2 * sample_idx_lowbits + 1) & 0xfffe;
        if (gt_low == gt_high) {
          CopyToUnalignedOffsetU16(cur_gt, &end_of_vector, 2 * sample_idx_lowbits + 1);
        }
        male_word &= male_word - 1;
      } while (male_word);
    }
  }
}

void FillBcf3allelicGtPloidy2(const PgenVariant* pgvp, const uint16_t* __restrict basic_genobytes4, uint32_t sample_ct, char* __restrict gt_start) {
  GenoarrLookup256x2bx4(pgvp->genovec, basic_genobytes4, sample_ct, gt_start);
  const uint32_t patch_01_ct = pgvp->patch_01_ct;
  if (patch_01_ct) {
    const uintptr_t* patch_01_set = pgvp->patch_01_set;
    const AlleleCode* patch_01_vals = pgvp->patch_01_vals;
    char* gt_offset1 = &(gt_start[1]);
    uintptr_t sample_idx_base = 0;
    uintptr_t patch_01_bits = patch_01_set[0];
    for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &patch_01_bits);
      gt_offset1[2 * sample_idx] = patch_01_vals[uii] * 2 + 2;
    }
  }
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  if (patch_10_ct) {
    const uintptr_t* patch_10_set = pgvp->patch_10_set;
    const DoubleAlleleCode* patch_10_vals_alias = R_CAST(const DoubleAlleleCode*, pgvp->patch_10_vals);
    uintptr_t sample_idx_base = 0;
    uintptr_t patch_10_bits = patch_10_set[0];
    for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &patch_10_bits);
      const uint16_t cur_u16 = patch_10_vals_alias[uii] * 2 + 0x202;
      CopyToUnalignedOffsetU16(CToUc(gt_start), &cur_u16, sample_idx);
    }
  }
}

// Convert homozygous and missing GT to ploidy 1.
// Homozygous GT may have been phased.
void FixBcf3allelicGtHh(uint32_t sample_ct, char* __restrict gt_start) {
#ifdef USE_SSE2
  const uint32_t fullvec_ct = sample_ct / (kBytesPerVec / 2);
  const VecU16 inv_m8 = vecu16_set1(0xff00);
  const VecU16 high_nophase_mask = vecu16_set1(0xfe00);
  const VecU16 eov = vecu16_set1(0x8100);
  char* gt_stop = &(gt_start[fullvec_ct * kBytesPerVec]);
  for (char* gt_iter = gt_start; gt_iter != gt_stop; gt_iter += kBytesPerVec) {
    const VecU16 vv_orig = vecu16_loadu(gt_iter);
    const VecU16 vv_lshift8 = vecu16_slli(vv_orig, 8);
    const VecU16 vv_high_nophase = vv_orig & high_nophase_mask;
    const VecU16 vv_replace_mask = (vv_lshift8 == vv_high_nophase) & inv_m8;
    const VecU16 vv_final = vecu16_blendv(vv_orig, eov, vv_replace_mask);
    vecu16_storeu(gt_iter, vv_final);
  }
  uint32_t sample_idx = fullvec_ct * (kBytesPerVec / 2);
#else  // !USE_SSE2
  uint32_t sample_idx = 0;
#endif
  for (; sample_idx != sample_ct; ++sample_idx) {
    const uint32_t gt_low = ctou32(gt_start[2 * sample_idx]);
    const uint32_t gt_high = ctou32(gt_start[2 * sample_idx + 1]) & 0xfe;
    if (gt_low == gt_high) {
      gt_start[2 * sample_idx + 1] = 0x81;
    }
  }
}

void FillBcf3allelicGtPloidy1(const PgenVariant* pgvp, const unsigned char* __restrict haploid_genobytes4, uint32_t sample_ct, char* __restrict gt_start) {
  GenoarrLookup256x1bx4(pgvp->genovec, haploid_genobytes4, sample_ct, gt_start);
  // patch_01_ct guaranteed to be zero; all patch_10 entries guaranteed to be
  // homozygous
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  const uintptr_t* patch_10_set = pgvp->patch_10_set;
  const AlleleCode* patch_10_vals = pgvp->patch_10_vals;
  uintptr_t sample_idx_base = 0;
  uintptr_t patch_10_bits = patch_10_set[0];
  for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
    const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &patch_10_bits);
    gt_start[sample_idx] = patch_10_vals[2 * uii] * 2 + 2;
  }
}

void FillBcf64allelicGtPloidy2(const PgenVariant* pgvp, const uint32_t* __restrict wide_genobytes4, uint32_t sample_ct, char* __restrict gt_start) {
  GenoarrLookup256x4bx4(pgvp->genovec, wide_genobytes4, sample_ct, gt_start);
  const uint32_t patch_01_ct = pgvp->patch_01_ct;
  if (patch_01_ct) {
    const uintptr_t* patch_01_set = pgvp->patch_01_set;
    const AlleleCode* patch_01_vals = pgvp->patch_01_vals;
    unsigned char* gt_offset1 = CToUc(&(gt_start[2]));
    uintptr_t sample_idx_base = 0;
    uintptr_t patch_01_bits = patch_01_set[0];
    for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_idx_base, &patch_01_bits);
      const uint16_t cur_u16 = patch_01_vals[uii] * 2 + 2;
      CopyToUnalignedOffsetU16(gt_offset1, &cur_u16, 2 * sample_idx);
    }
  }
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  if (patch_10_ct) {
    const uintptr_t* patch_10_set = pgvp->patch_10_set;
    const AlleleCode* patch_10_vals_iter = pgvp->patch_10_vals;
    uintptr_t sample_idx_base = 0;
    uintptr_t patch_10_bits = patch_10_set[0];
    for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
      const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &patch_10_bits);
      const uint32_t ac0 = *patch_10_vals_iter++;
      const uint32_t ac1 = *patch_10_vals_iter++;
      uint16_t cur_u16 = ac0 * 2 + 2;
      CopyToUnalignedOffsetU16(CToUc(gt_start), &cur_u16, 2 * sample_idx);
      cur_u16 = ac1 * 2 + 2;
      CopyToUnalignedOffsetU16(CToUc(gt_start), &cur_u16, 2 * sample_idx + 1);
    }
  }
}

void FillBcf64allelicGtPloidy1(const PgenVariant* pgvp, const uint16_t* __restrict wide_haploid_genobytes4, uint32_t sample_ct, char* __restrict gt_start) {
  GenoarrLookup256x2bx4(pgvp->genovec, wide_haploid_genobytes4, sample_ct, gt_start);
  // patch_01_ct guaranteed to be zero; all patch_10 entries guaranteed to be
  // homozygous
  const uint32_t patch_10_ct = pgvp->patch_10_ct;
  const uintptr_t* patch_10_set = pgvp->patch_10_set;
  const AlleleCode* patch_10_vals = pgvp->patch_10_vals;
  uintptr_t sample_idx_base = 0;
  uintptr_t patch_10_bits = patch_10_set[0];
  for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
    const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_idx_base, &patch_10_bits);
    const uint32_t ac = patch_10_vals[2 * uii];
    const uint16_t cur_u16 = ac * 2 + 2;
    CopyToUnalignedOffsetU16(CToUc(gt_start), &cur_u16, sample_idx);
  }
}

// Convert homozygous and missing GT to ploidy 1.
void FixBcf64allelicGtHh(uint32_t sample_ct, char* __restrict gt_start) {
#ifdef USE_SSE2
  const uint32_t fullvec_ct = sample_ct / (kBytesPerVec / 4);
  const VecU32 inv_m16 = vecu32_set1(0xffff0000U);
  const VecU32 high_nophase_mask = vecu32_set1(0xfffe0000U);
  const VecU32 eov = vecu32_set1(0x80010000U);
  char* gt_stop = &(gt_start[fullvec_ct * (k1LU * kBytesPerVec)]);
  for (char* gt_iter = gt_start; gt_iter != gt_stop; gt_iter += kBytesPerVec) {
    const VecU32 vv_orig = vecu32_loadu(gt_iter);
    const VecU32 vv_lshift16 = vecu32_slli(vv_orig, 16);
    const VecU32 vv_high_nophase = vv_orig & high_nophase_mask;
    const VecU32 vv_replace_mask = (vv_lshift16 == vv_high_nophase) & inv_m16;
    const VecU32 vv_final = vecu32_blendv(vv_orig, eov, vv_replace_mask);
    vecu32_storeu(gt_iter, vv_final);
  }
  uint32_t sample_idx = fullvec_ct * (kBytesPerVec / 4);
#else  // !USE_SSE2
  uint32_t sample_idx = 0;
#endif
  unsigned char* gt_start_uc = CToUc(gt_start);
  const uint16_t end_of_vector = 0x8001;
  for (; sample_idx != sample_ct; ++sample_idx) {
    const uint32_t gt_low = CopyFromUnalignedOffsetU16ZX(gt_start_uc, 2 * sample_idx);
    const uint32_t gt_high = CopyFromUnalignedOffsetU16ZX(gt_start_uc, 2 * sample_idx + 1) & 0xfffe;
    if (gt_low == gt_high) {
      CopyToUnalignedOffsetU16(gt_start_uc, &end_of_vector, 2 * sample_idx + 1);
    }
  }
}

// Assumes trailing bits of genovec are zeroed out to next word-PAIR.
// Assumes patch_01_set is zeroed out if patch_01_ct is zero; ditto for
// patch_10.
// See AppendVcfMultiallelicDsForce...().
void FillBcfMultiallelicDsForce(const PgenVariant* pgvp, const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_female, uint32_t sample_ct, uint32_t alt_allele_ct, uint32_t is_x, uint32_t is_y, uint32_t is_haploid, void* __restrict ds_start_orig) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uint32_t ds_wstride = alt_allele_ct * (kBitsPerWord * 4);
  const uintptr_t* genovec = pgvp->genovec;
  const uintptr_t* patch_01_set = pgvp->patch_01_set;
  const AlleleCode* patch_01_vals_iter = pgvp->patch_01_vals;
  const uintptr_t* patch_10_set = pgvp->patch_10_set;
  const AlleleCode* patch_10_vals_iter = pgvp->patch_10_vals;
  // float 1.0 = 0x3f800000  float 2.0 = 0x40000000
  const uint32_t one_fbits = 0x3f800000;
  const uint32_t ds_hom_val = ((!is_haploid) || is_x)? 0x40000000 : one_fbits;
  memset(ds_start_orig, 0, sample_ct * sizeof(int32_t) * alt_allele_ct);
  unsigned char* ds_start_uc = S_CAST(unsigned char*, ds_start_orig);
  const uint32_t missing_fbits = 0x7f800001;
  const uint32_t end_of_vector = 0x7f800002;
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    const uintptr_t geno_lo = genovec[2 * widx];
    const uintptr_t geno_hi = genovec[2 * widx + 1];
    if (!(geno_lo || geno_hi)) {
      continue;
    }
    // overflow of uint32 doesn't matter since we'd error out before writing
    // this variant
    unsigned char* ds_word_uc = &(ds_start_uc[widx * ds_wstride]);

    uintptr_t geno_word1 = PackWordToHalfwordMask5555(geno_lo) | (S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_hi)) << kBitsPerWordD2);
    uintptr_t geno_word2 = PackWordToHalfwordMaskAAAA(geno_lo) | (S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(geno_hi)) << kBitsPerWordD2);
    uintptr_t geno_missing_word = geno_word1 & geno_word2;
    geno_word1 ^= geno_missing_word;
    if (geno_word1) {
      uintptr_t patch_01_word = patch_01_set[widx];
      // 0/1
      geno_word1 ^= patch_01_word;
      while (geno_word1) {
        const uint32_t sample_idx_lowbits = ctzw(geno_word1);
        CopyToUnalignedOffsetU32(ds_word_uc, &one_fbits, sample_idx_lowbits * alt_allele_ct);
        geno_word1 &= geno_word1 - 1;
      }
      // 0/x
      if (patch_01_word) {
        unsigned char* ds_word_uc_sub1 = &(ds_word_uc[-4]);
        do {
          const uint32_t sample_idx_lowbits = ctzw(patch_01_word);
          const uint32_t ac = *patch_01_vals_iter++;
          CopyToUnalignedOffsetU32(ds_word_uc_sub1, &one_fbits, sample_idx_lowbits * alt_allele_ct + ac);
          patch_01_word &= patch_01_word - 1;
        } while (patch_01_word);
      }
    }
    geno_word2 ^= geno_missing_word;
    if (geno_word2) {
      uintptr_t patch_10_word = patch_10_set[widx];
      geno_word2 ^= patch_10_word;
      if (!is_x) {
        // 1/1
        while (geno_word2) {
          const uint32_t sample_idx_lowbits = ctzw(geno_word2);
          CopyToUnalignedOffsetU32(ds_word_uc, &ds_hom_val, sample_idx_lowbits * alt_allele_ct);
          geno_word2 &= geno_word2 - 1;
        }
        // x/y
        if (patch_10_word) {
          unsigned char* ds_word_uc_sub1 = &(ds_word_uc[-4]);
          do {
            const uint32_t sample_idx_lowbits = ctzw(patch_10_word);
            const uint32_t ac0 = *patch_10_vals_iter++;
            const uint32_t ac1 = *patch_10_vals_iter++;
            unsigned char* cur_ds_sub1 = &(ds_word_uc_sub1[sample_idx_lowbits * alt_allele_ct * 4]);
            if (ac0 == ac1) {
              CopyToUnalignedOffsetU32(cur_ds_sub1, &ds_hom_val, ac0);
            } else {
              CopyToUnalignedOffsetU32(cur_ds_sub1, &one_fbits, ac0);
              CopyToUnalignedOffsetU32(cur_ds_sub1, &one_fbits, ac1);
            }
            patch_10_word &= patch_10_word - 1;
          } while (patch_10_word);
        }
      } else {
        const uintptr_t male_word = sex_male[widx];
        while (geno_word2) {
          const uint32_t sample_idx_lowbits = ctzw(geno_word2);
          const uint32_t is_male = (male_word >> sample_idx_lowbits) & 1;
          const uint32_t cur_u32 = 0x40000000 - (is_male << 23);
          CopyToUnalignedOffsetU32(ds_word_uc, &cur_u32, sample_idx_lowbits * alt_allele_ct);
          geno_word2 &= geno_word2 - 1;
        }
        if (patch_10_word) {
          unsigned char* ds_word_uc_sub1 = &(ds_word_uc[-4]);
          do {
            const uint32_t sample_idx_lowbits = ctzw(patch_10_word);
            const uint32_t ac0 = *patch_10_vals_iter++;
            const uint32_t ac1 = *patch_10_vals_iter++;
            unsigned char* cur_ds_sub1 = &(ds_word_uc_sub1[sample_idx_lowbits * alt_allele_ct * 4]);
            if (ac0 == ac1) {
              const uint32_t is_male = (male_word >> sample_idx_lowbits) & 1;
              const uint32_t cur_u32 = 0x40000000 - (is_male << 23);
              CopyToUnalignedOffsetU32(cur_ds_sub1, &cur_u32, ac0);
            } else {
              CopyToUnalignedOffsetU32(cur_ds_sub1, &one_fbits, ac0);
              CopyToUnalignedOffsetU32(cur_ds_sub1, &one_fbits, ac1);
            }
            patch_10_word &= patch_10_word - 1;
          } while (patch_10_word);
        }
      }
    }
    if (geno_missing_word) {
      if (!is_y) {
        do {
          const uint32_t sample_idx_lowbits = ctzw(geno_missing_word);
          unsigned char* cur_ds_uc = &(ds_word_uc[sample_idx_lowbits * alt_allele_ct * 4]);
          CopyToUnalignedU32(cur_ds_uc, &missing_fbits);
          // To match VCF-encoding behavior, we store just one missing value
          // followed by END_OF_VECTOR, instead of a whole bunch of missing
          // values.
          for (uint32_t uii = 1; uii != alt_allele_ct; ++uii) {
            CopyToUnalignedOffsetU32(cur_ds_uc, &end_of_vector, uii);
          }
          geno_missing_word &= geno_missing_word - 1;
        } while (geno_missing_word);
      } else {
        const uintptr_t female_word = sex_female[widx];
        do {
          const uint32_t sample_idx_lowbits = ctzw(geno_missing_word);
          const uint32_t is_female = (female_word >> sample_idx_lowbits) & 1;
          unsigned char* cur_ds_uc = &(ds_word_uc[sample_idx_lowbits * alt_allele_ct * 4]);
          // end-of-vector instead of missing if female
          const uint32_t cur_u32 = 0x7f800001 + is_female;
          CopyToUnalignedU32(cur_ds_uc, &cur_u32);
          for (uint32_t uii = 1; uii != alt_allele_ct; ++uii) {
            CopyToUnalignedOffsetU32(cur_ds_uc, &end_of_vector, uii);
          }
          geno_missing_word &= geno_missing_word - 1;
        } while (geno_missing_word);
      }
    }
  }
}

void FillBcf3allelicPhasedGt(const PgenVariant* pgvp, uint32_t sample_ct, uint32_t is_haploid, uintptr_t* __restrict prev_phased, char* __restrict gt_start) {
  // * For every homozygous call, use the existing prev_phased bit.
  // * For every heterozygous call, use the phasepresent bit, update
  //   prev_phased, and if phased, swap GT entry from x|y to y|x when phaseinfo
  //   bit is set.
  // * For every missing call, do nothing.
  // * No need to handle chrX separately since FixBcfMaleXGtPloidy() converts
  //   all male phased-hom calls generated by this function into ploidy-1.
  //   Similarly, no need to handle chrY separately.
  // If we want to change this to a FillBcfPhasedGt() call followed by a
  // sparse-patch, phased[_hh]_genobytes2 should be extended to at least 251x2
  // entries.
  const uintptr_t* genovec = pgvp->genovec;
  const Halfword* patch_01_set_alias = DowncastKWToHW(pgvp->patch_01_set);
  const AlleleCode* patch_01_iter = pgvp->patch_01_vals;
  const Halfword* patch_10_set_alias = DowncastKWToHW(pgvp->patch_10_set);
  const AlleleCode* patch_10_iter = pgvp->patch_10_vals;
  const Halfword* phasepresent_alias = DowncastKWToHW(pgvp->phasepresent);
  const Halfword* phaseinfo_alias = DowncastKWToHW(pgvp->phaseinfo);
  uint16_t gt_base[4] = {0x202, 0x402, 0x404, 0};
  if (is_haploid) {
    gt_base[0] = 0x8102;
    gt_base[2] = 0x8104;
    gt_base[3] = 0x8100;
  }
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  Halfword* prev_phased_alias = DowncastWToHW(prev_phased);
  uint32_t loop_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genovec[widx];
    const uint32_t patch_01_set_hw = patch_01_set_alias[widx];
    const uint32_t patch_10_set_hw = patch_10_set_alias[widx];
    const uint32_t phasepresent_hw = phasepresent_alias[widx];
    const uint32_t phaseinfo_hw = phaseinfo_alias[widx];
    uint32_t prev_phased_hw = prev_phased_alias[widx];
    unsigned char* cur_gt = CToUc(&(gt_start[widx * kBitsPerWord]));
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
      const uintptr_t cur_geno = geno_word & 3;
      uint32_t cur_gt_base = gt_base[cur_geno];
      if ((cur_geno == 1) || (cur_geno == 2)) {
        const uint32_t cur_shift = 1U << sample_idx_lowbits;
        uint32_t is_het;
        if (cur_geno == 1) {
          is_het = 1;
          if (patch_01_set_hw & cur_shift) {
            const uint32_t ac = *patch_01_iter++;
            cur_gt_base = ((ac + 1) << 9) | 2;
          }
        } else {
          is_het = 0;
          if (patch_10_set_hw & cur_shift) {
            const uint32_t ac0 = *patch_10_iter++;
            const uint32_t ac1 = *patch_10_iter++;
            const uint32_t ac0_p1 = ac0 + 1;
            if (ac0 == ac1) {
              if (!is_haploid) {
                cur_gt_base = ac0_p1 * 0x202;
              } else {
                cur_gt_base = 0x8100 | (ac0_p1 << 1);
              }
            } else {
              is_het = 1;
              cur_gt_base = (ac0_p1 << 1) | ((ac1 + 1) << 9);
            }
          }
        }
        if (is_het) {
          if (phasepresent_hw & cur_shift) {
            prev_phased_hw |= cur_shift;
            if (phaseinfo_hw & cur_shift) {
              cur_gt_base = (cur_gt_base >> 8) | ((cur_gt_base & 0xff) << 8);
            }
          } else {
            prev_phased_hw &= ~cur_shift;
          }
        }
      }
      const uint16_t cur_u16 = cur_gt_base | (((prev_phased_hw >> sample_idx_lowbits) & 1) << 8);
      CopyToUnalignedOffsetU16(cur_gt, &cur_u16, sample_idx_lowbits);
      geno_word >>= 2;
    }
    prev_phased_alias[widx] = prev_phased_hw;
  }
}

void FillBcf64allelicPhasedGt(const PgenVariant* pgvp, uint32_t sample_ct, uint32_t is_haploid, uintptr_t* __restrict prev_phased, char* __restrict gt_start) {
  const uintptr_t* genovec = pgvp->genovec;
  const Halfword* patch_01_set_alias = DowncastKWToHW(pgvp->patch_01_set);
  const AlleleCode* patch_01_iter = pgvp->patch_01_vals;
  const Halfword* patch_10_set_alias = DowncastKWToHW(pgvp->patch_10_set);
  const AlleleCode* patch_10_iter = pgvp->patch_10_vals;
  const Halfword* phasepresent_alias = DowncastKWToHW(pgvp->phasepresent);
  const Halfword* phaseinfo_alias = DowncastKWToHW(pgvp->phaseinfo);
  uint32_t gt_base[4] = {0x20002, 0x40002, 0x40004, 0};
  if (is_haploid) {
    gt_base[0] = 0x80010002U;
    gt_base[2] = 0x80010004U;
    gt_base[3] = 0x80010000U;
  }
  const uint32_t sample_ctl2_m1 = (sample_ct - 1) / kBitsPerWordD2;
  Halfword* prev_phased_alias = DowncastWToHW(prev_phased);
  uint32_t loop_len = kBitsPerWordD2;
  for (uint32_t widx = 0; ; ++widx) {
    if (widx >= sample_ctl2_m1) {
      if (widx > sample_ctl2_m1) {
        return;
      }
      loop_len = ModNz(sample_ct, kBitsPerWordD2);
    }
    uintptr_t geno_word = genovec[widx];
    const uint32_t patch_01_set_hw = patch_01_set_alias[widx];
    const uint32_t patch_10_set_hw = patch_10_set_alias[widx];
    const uint32_t phasepresent_hw = phasepresent_alias[widx];
    const uint32_t phaseinfo_hw = phaseinfo_alias[widx];
    uint32_t prev_phased_hw = prev_phased_alias[widx];
    unsigned char* cur_gt = CToUc(&(gt_start[widx * (kBitsPerWord * 2)]));
    for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits != loop_len; ++sample_idx_lowbits) {
      const uintptr_t cur_geno = geno_word & 3;
      uint32_t cur_gt_base = gt_base[cur_geno];
      if ((cur_geno == 1) || (cur_geno == 2)) {
        const uint32_t cur_shift = 1U << sample_idx_lowbits;
        uint32_t is_het;
        if (cur_geno == 1) {
          is_het = 1;
          if (patch_01_set_hw & cur_shift) {
            const uint32_t ac = *patch_01_iter++;
            cur_gt_base = ((ac + 1) << 17) | 2;
          }
        } else {
          is_het = 0;
          if (patch_10_set_hw & cur_shift) {
            const uint32_t ac0 = *patch_10_iter++;
            const uint32_t ac1 = *patch_10_iter++;
            const uint32_t ac0_p1 = ac0 + 1;
            if (ac0 == ac1) {
              if (!is_haploid) {
                cur_gt_base = ac0_p1 * 0x20002;
              } else {
                cur_gt_base = 0x80010000U | (ac0_p1 << 1);
              }
            } else {
              is_het = 1;
              cur_gt_base = (ac0_p1 << 1) | ((ac1 + 1) << 17);
            }
          }
        }
        if (is_het) {
          if (phasepresent_hw & cur_shift) {
            prev_phased_hw |= cur_shift;
            if (phaseinfo_hw & cur_shift) {
              cur_gt_base = (cur_gt_base >> 16) | ((cur_gt_base & 0xffff) << 16);
            }
          } else {
            prev_phased_hw &= ~cur_shift;
          }
        }
      }
      const uint32_t cur_u32 = cur_gt_base | (((prev_phased_hw >> sample_idx_lowbits) & 1) << 16);
      CopyToUnalignedOffsetU32(cur_gt, &cur_u32, sample_idx_lowbits);
      geno_word >>= 2;
    }
    prev_phased_alias[widx] = prev_phased_hw;
  }
}

void FillBcfMultiallelicHdsForce(const PgenVariant* pgvp, const uintptr_t* __restrict sex_male, const uintptr_t* __restrict sex_female, uint32_t sample_ct, uint32_t alt_allele_ct, uint32_t is_x, uint32_t is_y, uint32_t is_haploid, void* __restrict hds_start_orig) {
  const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
  const uintptr_t* genovec = pgvp->genovec;
  const uintptr_t* patch_01_set = pgvp->patch_01_set;
  const AlleleCode* patch_01_vals_iter = pgvp->patch_01_vals;
  const uintptr_t* patch_10_set = pgvp->patch_10_set;
  const AlleleCode* patch_10_vals_iter = pgvp->patch_10_vals;
  uint32_t* hds_u32 = S_CAST(uint32_t*, hds_start_orig);
  if ((!is_haploid) || is_x || pgvp->phasepresent_ct) {
    const uintptr_t* phasepresent = pgvp->phasepresent;
    const uintptr_t* phaseinfo = pgvp->phaseinfo;
    memset(hds_start_orig, 0, sample_ct * sizeof(int32_t) * 2 * alt_allele_ct);
    const uint32_t hds_wstride = alt_allele_ct * 2 * kBitsPerWord;
    for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
      const uintptr_t geno_lo = genovec[2 * widx];
      const uintptr_t geno_hi = genovec[2 * widx + 1];
      if (!(geno_lo || geno_hi)) {
        continue;
      }
      // overflow of uint32 doesn't matter since we'd error out before writing
      // this variant
      uint32_t* hds_word_u32 = &(hds_u32[widx * hds_wstride]);
      const uintptr_t phasepresent_word = phasepresent[widx];
      const uintptr_t phaseinfo_word = phaseinfo[widx];

      uintptr_t geno_word1 = PackWordToHalfwordMask5555(geno_lo) | (S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_hi)) << kBitsPerWordD2);
      uintptr_t geno_word2 = PackWordToHalfwordMaskAAAA(geno_lo) | (S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(geno_hi)) << kBitsPerWordD2);
      uintptr_t geno_missing_word = geno_word1 & geno_word2;
      geno_word1 ^= geno_missing_word;
      if (geno_word1) {
        uintptr_t patch_01_word = patch_01_set[widx];
        // 0/1
        geno_word1 ^= patch_01_word;
        while (geno_word1) {
          const uint32_t sample_idx_lowbits = ctzw(geno_word1);
          const uintptr_t lowbit = geno_word1 & (-geno_word1);
          uint32_t* cur_hds = &(hds_word_u32[sample_idx_lowbits * alt_allele_ct * 2]);
          if (!(phasepresent_word & lowbit)) {
            // float 0.5
            cur_hds[0] = 0x3f000000;
            cur_hds[alt_allele_ct] = 0x3f000000;
          } else {
            // if phaseinfo bit unset, second entry is 1.0; otherwise, first
            // entry.
            cur_hds[alt_allele_ct * (1 & (~(phaseinfo_word >> sample_idx_lowbits)))] = 0x3f800000;
          }
          geno_word1 ^= lowbit;
        }
        // 0/x
        if (patch_01_word) {
          uint32_t* hds_word_u32_sub1 = &(hds_word_u32[-1]);
          do {
            const uint32_t sample_idx_lowbits = ctzw(patch_01_word);
            const uintptr_t lowbit = patch_01_word & (-patch_01_word);
            const uint32_t ac = *patch_01_vals_iter++;
            uint32_t* cur_hds_ac = &(hds_word_u32_sub1[sample_idx_lowbits * alt_allele_ct * 2 + ac]);
            if (!(phasepresent_word & lowbit)) {
              cur_hds_ac[0] = 0x3f000000;
              cur_hds_ac[alt_allele_ct] = 0x3f000000;
            } else {
              cur_hds_ac[alt_allele_ct * (1 & (~(phaseinfo_word >> sample_idx_lowbits)))] = 0x3f800000;
            }
            patch_01_word ^= lowbit;
          } while (patch_01_word);
        }
      }
      geno_word2 ^= geno_missing_word;
      if (geno_word2) {
        uintptr_t patch_10_word = patch_10_set[widx];
        geno_word2 ^= patch_10_word;
        // 1/1
        while (geno_word2) {
          const uint32_t sample_idx_lowbits = ctzw(geno_word2);
          uint32_t* cur_hds = &(hds_word_u32[sample_idx_lowbits * 2 * alt_allele_ct]);
          cur_hds[0] = 0x3f800000;
          cur_hds[alt_allele_ct] = 0x3f800000;
          geno_word2 &= geno_word2 - 1;
        }
        // x/y
        if (patch_10_word) {
          do {
            const uint32_t sample_idx_lowbits = ctzw(patch_10_word);
            const uint32_t ac0 = *patch_10_vals_iter++;
            const uint32_t ac1 = *patch_10_vals_iter++;
            const uint32_t ac0_m1 = ac0 - 1;
            uint32_t* cur_hds = &(hds_word_u32[sample_idx_lowbits * 2 * alt_allele_ct]);
            if (ac0 == ac1) {
              cur_hds[ac0_m1] = 0x3f800000;
              cur_hds[ac0_m1 + alt_allele_ct] = 0x3f800000;
            } else {
              const uintptr_t lowbit = patch_10_word & (-patch_10_word);
              const uint32_t ac1_m1 = ac1 - 1;
              if (!(phasepresent_word & lowbit)) {
                cur_hds[ac0_m1] = 0x3f000000;
                cur_hds[ac1_m1] = 0x3f000000;
                cur_hds[alt_allele_ct + ac0_m1] = 0x3f000000;
                cur_hds[alt_allele_ct + ac1_m1] = 0x3f000000;
              } else {
                const uintptr_t phaseinfo_bit = (phaseinfo_word >> sample_idx_lowbits) & 1;
                cur_hds[ac0_m1 + alt_allele_ct * phaseinfo_bit] = 0x3f800000;
                cur_hds[ac1_m1 + alt_allele_ct * (1 ^ phaseinfo_bit)] = 0x3f800000;
              }
            }
            patch_10_word &= patch_10_word - 1;
          } while (patch_10_word);
        }
      }
      if (geno_missing_word) {
        if (!is_y) {
          do {
            const uint32_t sample_idx_lowbits = ctzw(geno_missing_word);
            uint32_t* cur_hds_u32 = &(hds_word_u32[sample_idx_lowbits * 2 * alt_allele_ct]);
            // To match VCF-encoding behavior, we only store 1 or 2 missing
            // values (depending on ploidy), instead of a whole bunch of
            // missing values.
            cur_hds_u32[0] = 0x7f800001;
            cur_hds_u32[1] = 0x7f800001;
            for (uint32_t uii = 2; uii != 2 * alt_allele_ct; ++uii) {
              cur_hds_u32[uii] = 0x7f800002;
            }
            geno_missing_word &= geno_missing_word - 1;
          } while (geno_missing_word);
        } else {
          const uintptr_t female_word = sex_female[widx];
          do {
            const uint32_t sample_idx_lowbits = ctzw(geno_missing_word);
            const uint32_t is_female = (female_word >> sample_idx_lowbits) & 1;
            uint32_t* cur_hds_u32 = &(hds_word_u32[sample_idx_lowbits * 2 * alt_allele_ct]);
            // end-of-vector instead of missing if female
            cur_hds_u32[0] = 0x7f800001 + is_female;
            // Only need to fill first half; second half will be filled with
            // END_OF_VECTOR values later, since this is chrY
            for (uint32_t uii = 1; uii != alt_allele_ct; ++uii) {
              cur_hds_u32[uii] = 0x7f800002;
            }
            geno_missing_word &= geno_missing_word - 1;
          } while (geno_missing_word);
        }
      }
    }
    if (!is_haploid) {
      return;
    }
    uint32_t* hds_second = &(hds_u32[alt_allele_ct]);
    if (!is_x) {
      // Set all ploidies to 1 except when phase present.
      for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
        uintptr_t unphased_word = sex_male[widx] & (~(phasepresent[widx]));
        if (!unphased_word) {
          continue;
        }
        uint32_t* hds_word_second_u32 = &(hds_second[widx * hds_wstride]);
        do {
          const uint32_t sample_idx_lowbits = ctzw(unphased_word);
          uint32_t* cur_hds_u32 = &(hds_word_second_u32[sample_idx_lowbits * 2 * alt_allele_ct]);
          for (uint32_t uii = 0; uii != alt_allele_ct; ++uii) {
            cur_hds_u32[uii] = 0x7f800002;  // END_OF_VECTOR
          }
          unphased_word &= unphased_word - 1;
        } while (unphased_word);
      }
    } else {
      // Set male ploidy to 1 except when phase present.
      for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
        uintptr_t male_unphased_word = sex_male[widx] & (~(phasepresent[widx]));
        if (!male_unphased_word) {
          continue;
        }
        uint32_t* hds_word_second_u32 = &(hds_second[widx * hds_wstride]);
        do {
          const uint32_t sample_idx_lowbits = ctzw(male_unphased_word);
          uint32_t* cur_hds_u32 = &(hds_word_second_u32[sample_idx_lowbits * 2 * alt_allele_ct]);
          for (uint32_t uii = 0; uii != alt_allele_ct; ++uii) {
            cur_hds_u32[uii] = 0x7f800002;
          }
          male_unphased_word &= male_unphased_word - 1;
        } while (male_unphased_word);
      }
    }
    return;
  }
  memset(hds_start_orig, 0, sample_ct * sizeof(int32_t) * alt_allele_ct);
  const uint32_t hds_wstride = alt_allele_ct * kBitsPerWord;
  for (uint32_t widx = 0; widx != sample_ctl; ++widx) {
    const uintptr_t geno_lo = genovec[2 * widx];
    const uintptr_t geno_hi = genovec[2 * widx + 1];
    if (!(geno_lo || geno_hi)) {
      continue;
    }
    uint32_t* hds_word_u32 = &(hds_u32[widx * hds_wstride]);

    uintptr_t geno_word1 = PackWordToHalfwordMask5555(geno_lo) | (S_CAST(uintptr_t, PackWordToHalfwordMask5555(geno_hi)) << kBitsPerWordD2);
    uintptr_t geno_word2 = PackWordToHalfwordMaskAAAA(geno_lo) | (S_CAST(uintptr_t, PackWordToHalfwordMaskAAAA(geno_hi)) << kBitsPerWordD2);
    uintptr_t geno_missing_word = geno_word1 & geno_word2;
    geno_word1 ^= geno_missing_word;
    if (geno_word1) {
      uintptr_t patch_01_word = patch_01_set[widx];
      // 0/1
      geno_word1 ^= patch_01_word;
      while (geno_word1) {
        const uint32_t sample_idx_lowbits = ctzw(geno_word1);
        hds_word_u32[sample_idx_lowbits * alt_allele_ct] = 0x3f000000;
        geno_word1 &= geno_word1 - 1;
      }
      // 0/x
      if (patch_01_word) {
        uint32_t* hds_word_u32_sub1 = &(hds_word_u32[-1]);
        do {
          const uint32_t sample_idx_lowbits = ctzw(patch_01_word);
          const uint32_t ac = *patch_01_vals_iter++;
          hds_word_u32_sub1[sample_idx_lowbits * alt_allele_ct + ac] = 0x3f000000;
          patch_01_word &= patch_01_word - 1;
        } while (patch_01_word);
      }
    }
    geno_word2 ^= geno_missing_word;
    if (geno_word2) {
      uintptr_t patch_10_word = patch_10_set[widx];
      geno_word2 ^= patch_10_word;
      // 1/1
      while (geno_word2) {
        const uint32_t sample_idx_lowbits = ctzw(geno_word2);
        hds_word_u32[sample_idx_lowbits * alt_allele_ct] = 0x3f800000;
        geno_word2 &= geno_word2 - 1;
      }
      // x/y
      if (patch_10_word) {
        do {
          const uint32_t sample_idx_lowbits = ctzw(patch_10_word);
          const uint32_t ac0 = *patch_10_vals_iter++;
          const uint32_t ac1 = *patch_10_vals_iter++;
          uint32_t* cur_hds = &(hds_word_u32[sample_idx_lowbits * alt_allele_ct]);
          if (ac0 == ac1) {
            cur_hds[ac0 - 1] = 0x3f800000;
          } else {
            cur_hds[ac0 - 1] = 0x3f000000;
            cur_hds[ac1 - 1] = 0x3f000000;
          }
          patch_10_word &= patch_10_word - 1;
        } while (patch_10_word);
      }
    }
    if (geno_missing_word) {
      if (!is_y) {
        do {
          const uint32_t sample_idx_lowbits = ctzw(geno_missing_word);
          uint32_t* cur_hds_u32 = &(hds_word_u32[sample_idx_lowbits * alt_allele_ct]);
          cur_hds_u32[0] = 0x7f800001;
          for (uint32_t uii = 1; uii != alt_allele_ct; ++uii) {
            cur_hds_u32[uii] = 0x7f800002;
          }
          geno_missing_word &= geno_missing_word - 1;
        } while (geno_missing_word);
      } else {
        const uintptr_t female_word = sex_female[widx];
        do {
          const uint32_t sample_idx_lowbits = ctzw(geno_missing_word);
          const uint32_t is_female = (female_word >> sample_idx_lowbits) & 1;
          uint32_t* cur_hds_u32 = &(hds_word_u32[sample_idx_lowbits * alt_allele_ct]);
          cur_hds_u32[0] = 0x7f800001 + is_female;
          for (uint32_t uii = 1; uii != alt_allele_ct; ++uii) {
            cur_hds_u32[uii] = 0x7f800002;
          }
          geno_missing_word &= geno_missing_word - 1;
        } while (geno_missing_word);
      }
    }
  }
}

void MakeRlenWarningStr(const char* variant_id) {
  char* err_write_iter = strcpya_k(g_logbuf, "Warning: Unexpected INFO/END value for variant");
  // 50 chars + strlen(variant_id), split into two lines if >79
  const uint32_t variant_id_slen = strlen(variant_id);
  *err_write_iter++ = (variant_id_slen < 30)? ' ' : '\n';
  *err_write_iter++ = '\'';
  err_write_iter = memcpya(err_write_iter, variant_id, variant_id_slen);
  strcpy_k(err_write_iter, "'.\n");
}

static_assert(sizeof(AlleleCode) == 1, "ExportBcf() needs to be updated.");
PglErr ExportBcf(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const SampleIdInfo* siip, const uintptr_t* sex_male_collapsed, const uintptr_t* sex_female_collapsed, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const uint32_t* contig_lens, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, ExportfFlags exportf_flags, VcfExportMode vcf_mode, IdpasteFlags exportf_id_paste, char exportf_id_delim, char* xheader, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  TextStream pvar_reload_txs;
  BgzfCompressStream bgzf;
  PreinitTextStream(&pvar_reload_txs);
  PreinitBgzfCompressStream(&bgzf);
  {
    {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".bcf");
      reterr = InitBgzfCompressStreamEx(outname, 0, kBgzfDefaultClvl, max_thread_ct, &bgzf);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, outname, strerror(errno));
        }
        goto ExportBcf_ret_1;
      }
    }
    const uint32_t max_chr_slen = GetMaxChrSlen(cip);
    uint32_t write_some_dosage = 0;
    uint32_t write_ds = 0;
    uint32_t write_hds = 0;
    uint32_t ds_force = 0;
    uint32_t ds_only = 0;
    uint32_t hds_force = 0;
    if (vcf_mode != kVcfExport0) {
      write_some_dosage = 1;
      if (vcf_mode != kVcfExportGp) {
        write_ds = 1;
        if (vcf_mode == kVcfExportDsForce) {
          ds_force = 1;
        } else if (vcf_mode == kVcfExportDsOnly) {
          if (unlikely((cip->haploid_mask[0] & 1) ||
                       XymtIsNonempty(variant_include, cip, kChrOffsetX) ||
                       XymtIsNonempty(variant_include, cip, kChrOffsetY) ||
                       XymtIsNonempty(variant_include, cip, kChrOffsetMT))) {
            logerrputs("Error: --export vcf-dosage=DS-only mode can only be used with all-diploid data.\n");
            goto ExportBcf_ret_INCONSISTENT_INPUT;
          }
          ds_force = 1;
          ds_only = 1;
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
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t all_nonref = (pgfip->gflags & kfPgenGlobalAllNonref) && (!nonref_flags);
    const uint32_t raw_variant_ctl = BitCtToWordCt(raw_variant_ct);
    uint32_t write_pr = all_nonref;
    if (nonref_flags) {
      write_pr = !IntersectionIsEmpty(variant_include, nonref_flags, raw_variant_ctl);
    }
    const uint32_t info_pr_flag_present = (info_flags / kfInfoPrFlagPresent) & 1;
    // Count FILTER/INFO/FORMAT keys for hash table sizing purposes.  Slight
    // overestimate ok, so we don't bother to avoid double-counting
    // FILTER/PASS, or counting FORMAT/HDS when it won't actually be used, etc.
    uint32_t fif_key_ct_ubound = 5;
    if (xheader) {
      const char* xheader_end = &(xheader[xheader_blen]);
      for (const char* xheader_iter = xheader; xheader_iter != xheader_end; xheader_iter = AdvPastDelim(xheader_iter, '\n')) {
        fif_key_ct_ubound += StrStartsWithUnsafe(xheader_iter, "##FILTER=<") || StrStartsWithUnsafe(xheader_iter, "##INFO=<");
      }
    }
    const uint32_t fif_keys_htable_size = GetHtableFastSize(fif_key_ct_ubound);
    // The byte size of the header must be stored in bytes [5,9).  So we
    // generate the entire header, fill in its size, and only then do we flush.
    // (In principle, we could compute its exact size and then generate the
    // contents incrementally to save memory in extreme cases, but it's limited
    // to 4 GiB anyway so I won't bother.)
    // Ways for the actual header to exceed xheader_blen bytes:
    // * New FILTER/PASS, INFO/PR, and FORMAT header lines, and beginning of
    //   #CHROM line.  These add up to less than 1 KB.
    // * ",IDX=<#>" appended to existing FILTER/INFO header lines.  The number
    //   of added bytes is less than (5 + UintSlen(fif_key_ct_ubound - 1)) *
    //   fif_key_ct_ubound.
    // * ",IDX=<#>" appended to existing contig lines, or entirely new contig
    //   lines.  This is less than chr_ct * <small constant> + sum of extra-chr
    //   string string lengths.
    // * ",length=<#>" inserted into existing contig lines.  This is also
    //   accounted for in the previous term.
    // * Sample IDs.
    uintptr_t header_ubound = xheader_blen + 1024 + (5 + UintSlen(fif_key_ct_ubound - 1)) * S_CAST(uintptr_t, fif_key_ct_ubound);
    const uint32_t chr_ct = cip->chr_ct;
    {
      const uint32_t max_code = cip->max_code;
      for (uint32_t chr_fo_idx = 0; chr_fo_idx != chr_ct; ++chr_fo_idx) {
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        if (!IsSet(cip->chr_mask, chr_idx)) {
          continue;
        }
        // ##contig=<ID=chr99,length=1000000000,IDX=12345>
        // 37 + kMaxChrCodeDigits + strlen(EOLN_STR) + length of contig ID
        header_ubound += 37 + kMaxChrCodeDigits + strlen(EOLN_STR);
        if (chr_idx <= max_code) {
          // regular code
          header_ubound += kMaxChrTextnumSlen + 3;
        } else {
          header_ubound += strlen(cip->nonstd_names[chr_idx]);
        }
      }
      const uintptr_t max_sample_id_blen = siip->max_sample_id_blen;
      uintptr_t max_sid_blen = siip->max_sid_blen;
      const char* sids = siip->sids;
      const uint32_t write_sid = DataSidColIsRequired(sample_include, sids, sample_ct, max_sid_blen, exportf_id_paste / kfIdpasteMaybesid);
      if (write_sid && (!sids)) {
        max_sid_blen = 2;
      }
      const uintptr_t max_exported_sample_id_blen = max_sample_id_blen + write_sid * max_sid_blen;
      header_ubound += max_exported_sample_id_blen * sample_ct;
    }
    // We need these three data structures in the main BCF-writing loop.
    char** fif_keys_mutable;
    uint32_t* fif_keys_htable;
    ChrIdx* chr_fo_to_bcf_idx;
    // We'll also need writebuf, but it'll be resized.
    // All bigstack-base allocations after writebuf in this section are
    // temporary.
    char* writebuf;
    if (unlikely(bigstack_alloc_u32(fif_keys_htable_size, &fif_keys_htable) ||
                 bigstack_calloc_cp(fif_key_ct_ubound, &fif_keys_mutable) ||
                 bigstack_alloc_chridx(chr_ct, &chr_fo_to_bcf_idx) ||
                 bigstack_alloc_c(header_ubound + 9, &writebuf))) {
      goto ExportBcf_ret_NOMEM;
    }
    SetAllU32Arr(fif_keys_htable_size, fif_keys_htable);
    char* write_iter = memcpya_k(writebuf, "BCF\2\2\0\0\0", 9);
    const uint32_t v43 = (exportf_flags / kfExportfBcf43) & 1;
    uint32_t pr_key_idx = UINT32_MAX;
    uint32_t end_key_idx = UINT32_MAX;
    uint32_t dosage_key_idx = UINT32_MAX;
    uint32_t hds_key_idx = UINT32_MAX;
    uint32_t gt_key_idx = UINT32_MAX;
    uint32_t info_end_exists = 0;
    uint32_t fif_key_ct = 0;
    {
      // Now generate the header body and build the string dictionaries.
      // We store additional information about each FILTER/INFO/FORMAT key in
      // fif_keys[index][-1]:
      //   bits 0-2: INFO type (1 = Flag, 3 = Integer, 5 = Float,
      //             7 = String/Character)
      //   bit 3: FILTER?
      //   bits 4-7 unused for now
      char* header_start = write_iter;
      AppendVcfHeaderStart(v43, &write_iter);
      if (cip->chrset_source) {
        AppendChrsetLine(cip, &write_iter);
      }
      const uint32_t chr_ctl = BitCtToWordCt(chr_ct);
      uintptr_t* written_contig_header_lines;
      if (unlikely(bigstack_calloc_w(chr_ctl, &written_contig_header_lines))) {
        goto ExportBcf_ret_NOMEM;
      }
      unsigned char* tmp_alloc_end = g_bigstack_end;
      // Match htslib/bcftools behavior of always putting an explicit
      // FILTER/PASS header line first.
      write_iter = strcpya_k(write_iter, "##FILTER=<ID=PASS,Description=\"All filters passed\",IDX=0>" EOLN_STR);
      {
        uint32_t cur_idx;
        reterr = AddToFifHtable(g_bigstack_base, "PASS", fif_keys_htable_size, strlen("PASS"), 8, &tmp_alloc_end, fif_keys_mutable, fif_keys_htable, &fif_key_ct, &cur_idx);
        if (unlikely(reterr)) {
          goto ExportBcf_ret_1;
        }
      }

      const uint32_t x_code = cip->xymt_codes[kChrOffsetX];
      const uint32_t par1_code = cip->xymt_codes[kChrOffsetPAR1];
      const uint32_t par2_code = cip->xymt_codes[kChrOffsetPAR2];
      uint32_t contig_zero_written = 0;
      uint32_t x_contig_line_written = 0;
      uint32_t contig_string_idx_end = 0;
      if (xheader) {
        char* xheader_end = &(xheader[xheader_blen]);
        for (char* line_end = xheader; line_end != xheader_end; ) {
          char* line_start = line_end;
          line_end = AdvPastDelim(&(line_start[1]), '\n');
          if (line_start[1] != '#') {
            continue;
          }
          char* line_main = &(line_start[2]);
          const uint32_t is_filter_line = StrStartsWithUnsafe(line_main, "FILTER=<");
          const uint32_t is_info_line = StrStartsWithUnsafe(line_main, "INFO=<");
          if ((!is_filter_line) && (!is_info_line)) {
            if (!StrStartsWithUnsafe(line_main, "contig=<")) {
              const uint32_t slen = line_end - line_start;
              if (unlikely(!ValidVcfMetaLine(line_start, slen, v43))) {
                goto ExportBcf_ret_INCONSISTENT_INPUT;
              }
              write_iter = memcpya(write_iter, line_start, slen);
              continue;
            }
            char* hkvline_iter = &(line_main[strlen("contig=<")]);
            char* idval;
            uint32_t id_slen;
            if (unlikely(HkvlineId(&hkvline_iter, &idval, &id_slen))) {
              logerrputs("Error: Malformed ##contig header line in .pvar file.\n");
              goto ExportBcf_ret_MALFORMED_INPUT;
            }
            if (hkvline_iter[-1] == '>') {
              // Otherwise-empty line.  Regenerate instead of copying it.
              continue;
            }
            // if GetChrCodeCounted() is modified to not mutate idval, xheader
            // can be changed to const char*
            const uint32_t chr_idx = GetChrCodeCounted(cip, id_slen, idval);
            if (IsI32Neg(chr_idx) || (!IsSet(cip->chr_mask, chr_idx)) || (chr_idx == par1_code) || (chr_idx == par2_code)) {
              continue;
            }
            if (chr_idx == x_code) {
              x_contig_line_written = 1;
            }
            const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[chr_idx];
            if (unlikely(IsSet(written_contig_header_lines, chr_fo_idx))) {
              logerrputs("Error: Duplicate ##contig line in .pvar file.\n");
              goto ExportBcf_ret_MALFORMED_INPUT;
            }
            SetBit(chr_fo_idx, written_contig_header_lines);
            char* contig_write_start = strcpya_k(write_iter, "##contig=<ID=");
            char* contig_write_end = chrtoa(cip, chr_idx, contig_write_start);
            if ((*contig_write_start == '0') && (contig_write_end == &(contig_write_start[1]))) {
              // --allow-extra-chr 0 special case
              // note that write_iter has *not* been advanced
              contig_zero_written = 1;  // technically we write this later
              continue;
            }
            write_iter = contig_write_end;
            if (unlikely(!ValidVcfContigName(contig_write_start, write_iter, v43))) {
              goto ExportBcf_ret_MALFORMED_INPUT;
            }
            char* closing_gt = &(line_end[-2]);
            while (ctou32(*closing_gt) <= 32) {
              --closing_gt;
            }
            if (unlikely(*closing_gt != '>')) {
              logerrputs("Error: ##contig header line in .pvar file doesn't end with '>'.\n");
              goto ExportBcf_ret_MALFORMED_INPUT;
            }
            if (contig_lens && contig_lens[chr_idx]) {
              if (unlikely(FixAndWriteContigHeaderLine(hkvline_iter, idval, closing_gt, id_slen, contig_lens[chr_idx], &write_iter))) {
                goto ExportBcf_ret_MALFORMED_INPUT;
              }
            } else {
              write_iter = WriteRestOfHeaderLine(hkvline_iter, idval, closing_gt, id_slen, write_iter);
            }
            write_iter = strcpya_k(write_iter, ",IDX=");
            write_iter = u32toa(contig_string_idx_end, write_iter);
            write_iter = strcpya_k(write_iter, ">" EOLN_STR);
            chr_fo_to_bcf_idx[chr_fo_idx] = contig_string_idx_end++;
            continue;
          }
          char* hkvline_iter = &(line_main[8 - 2 * is_info_line]);
          char* idval;
          uint32_t id_slen;
          if (unlikely(HkvlineId(&hkvline_iter, &idval, &id_slen))) {
            // This part of INFO line should have already been lexed
            // successfully during .pvar load.
            assert(is_filter_line);
            logerrputs("Error: Malformed FILTER header line in .pvar file.\n");
            goto ExportBcf_ret_MALFORMED_INPUT;
          }
          char* idval_end = &(idval[id_slen]);
          if (unlikely(id_slen > kMaxIdSlen)) {
            logerrputs("Error: " PROG_NAME_STR " does not support FILTER/INFO keys longer than " MAX_ID_SLEN_STR " characters.\n");
            goto ExportBcf_ret_MALFORMED_INPUT;
          }
          unsigned char prechar = 8;
          if (is_info_line) {
            char* typestr;
            uint32_t type_slen;
            if (unlikely(HkvlineFind(hkvline_iter, "Type", &typestr, &type_slen))) {
              *idval_end = '\0';
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid INFO/%s header line in .pvar file.\n", idval);
              goto ExportBcf_ret_MALFORMED_INPUT_WW;
            }
            if (strequal_k(typestr, "Flag", type_slen)) {
              prechar = 1;
            } else if (strequal_k(typestr, "Integer", type_slen)) {
              prechar = 3;
              if (!info_end_exists) {
                // ignore when this isn't of integer type
                info_end_exists = strequal_k(idval, "END", id_slen);
              }
            } else if (strequal_k(typestr, "Float", type_slen)) {
              prechar = 5;
            } else if (likely(strequal_k(typestr, "String", type_slen) ||
                              strequal_k(typestr, "Character", type_slen))) {
              prechar = 7;
            } else {
              *idval_end = '\0';
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid INFO/%s header line in .pvar file.\n", idval);
              goto ExportBcf_ret_MALFORMED_INPUT_WW;
            }
          } else if (strequal_k(idval, "PASS", id_slen)) {
            // Don't duplicate this.
            continue;
          }
          char* closing_gt = &(line_end[-2]);
          while (ctou32(*closing_gt) <= 32) {
            --closing_gt;
          }
          if (unlikely(*closing_gt != '>')) {
            *idval_end = '\0';
            snprintf(g_logbuf, kLogbufSize, "Error: %s:%s header line in .pvar file doesn't end with '>'.\n", is_info_line? "INFO" : "FILTER", idval);
            goto ExportBcf_ret_MALFORMED_INPUT;
          }
          uint32_t cur_idx;
          reterr = AddToFifHtable(g_bigstack_base, idval, fif_keys_htable_size, id_slen, prechar, &tmp_alloc_end, fif_keys_mutable, fif_keys_htable, &fif_key_ct, &cur_idx);
          if (unlikely(reterr)) {
            goto ExportBcf_ret_1;
          }
          write_iter = memcpya(write_iter, line_start, closing_gt - line_start);
          write_iter = strcpya_k(write_iter, ",IDX=");
          write_iter = u32toa(cur_idx, write_iter);
          write_iter = strcpya_k(write_iter, ">" EOLN_STR);
        }
      }
      // fill in the missing ##contig lines
      if (contig_zero_written) {
        write_iter = strcpya_k(write_iter, "##contig=<ID=0>" EOLN_STR);
      }
      uint32_t chrx_exists = 0;
      for (uint32_t chr_fo_idx = 0; chr_fo_idx != cip->chr_ct; ++chr_fo_idx) {
        if (IsSet(written_contig_header_lines, chr_fo_idx)) {
          continue;
        }
        const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        if ((!IsSet(cip->chr_mask, chr_idx)) || AllBitsAreZero(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1])) {
          continue;
        }
        if ((chr_idx == x_code) || (chr_idx == par1_code) || (chr_idx == par2_code)) {
          chrx_exists = 1;
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
          write_iter = chr_name_write_end;
        } else {
          write_iter = chr_name_write_end;
          if (contig_lens && contig_lens[chr_idx]) {
            write_iter = strcpya_k(write_iter, ",length=");
            write_iter = u32toa(contig_lens[chr_idx], write_iter);
          }
        }
        write_iter = strcpya_k(write_iter, ",IDX=");
        write_iter = u32toa(contig_string_idx_end, write_iter);
        write_iter = strcpya_k(write_iter, ">" EOLN_STR);
        chr_fo_to_bcf_idx[chr_fo_idx] = contig_string_idx_end++;
      }
      if (chrx_exists) {
        uint32_t x_contig_string_idx;
        if (!x_contig_line_written) {
          char* chr_name_write_start = strcpya_k(write_iter, "##contig=<ID=");
          char* chr_name_write_end = chrtoa(cip, x_code, chr_name_write_start);
          write_iter = strcpya_k(chr_name_write_end, ",IDX=");
          write_iter = u32toa(contig_string_idx_end, write_iter);
          write_iter = strcpya_k(write_iter, ">" EOLN_STR);
          x_contig_string_idx = contig_string_idx_end;
          if ((!IsI32Neg(x_code)) && IsSet(cip->chr_mask, x_code)) {
            const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[x_code];
            chr_fo_to_bcf_idx[chr_fo_idx] = contig_string_idx_end;
          }
          ++contig_string_idx_end;
        } else {
          const uint32_t x_chr_fo_idx = cip->chr_idx_to_foidx[x_code];
          x_contig_string_idx = chr_fo_to_bcf_idx[x_chr_fo_idx];
        }
        if ((!IsI32Neg(par1_code)) && IsSet(cip->chr_mask, par1_code)) {
          const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[par1_code];
          chr_fo_to_bcf_idx[chr_fo_idx] = x_contig_string_idx;
        }
        if ((!IsI32Neg(par2_code)) && IsSet(cip->chr_mask, par2_code)) {
          const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[par2_code];
          chr_fo_to_bcf_idx[chr_fo_idx] = x_contig_string_idx;
        }
      }
      BigstackReset(written_contig_header_lines);
      if (write_pr) {
        if (unlikely(info_flags & kfInfoPrNonflagPresent)) {
          logerrputs("Error: Conflicting INFO/PR definitions.  Either fix all REF alleles so that the\n'provisional reference' flag is no longer needed, or remove/rename the other\nuse of the INFO/PR key.\n");
          goto ExportBcf_ret_INCONSISTENT_INPUT;
        }
        if (!info_pr_flag_present) {
          write_iter = strcpya_k(write_iter, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\",IDX=");
          reterr = AddToFifHtable(g_bigstack_base, "PR", fif_keys_htable_size, strlen("PR"), 1, &tmp_alloc_end, fif_keys_mutable, fif_keys_htable, &fif_key_ct, &pr_key_idx);
          if (unlikely(reterr)) {
            goto ExportBcf_ret_1;
          }
          write_iter = u32toa(pr_key_idx, write_iter);
          write_iter = strcpya_k(write_iter, ">" EOLN_STR);
        } else {
          pr_key_idx = IdHtableFind("PR", TO_CONSTCPCONSTP(fif_keys_mutable), fif_keys_htable, strlen("PR"), fif_keys_htable_size);
          assert(pr_key_idx != UINT32_MAX);
        }
      }
      end_key_idx = IdHtableFind("END", TO_CONSTCPCONSTP(fif_keys_mutable), fif_keys_htable, strlen("END"), fif_keys_htable_size);
      if (end_key_idx != UINT32_MAX) {
        if ((fif_keys_mutable[end_key_idx][-1] & 7) != 3) {
          // Not an INFO key of type integer; ignore.
          end_key_idx = UINT32_MAX;
        }
      }
      if (write_ds) {
        write_iter = strcpya_k(write_iter, "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\",IDX=");
        reterr = AddToFifHtable(g_bigstack_base, "DS", fif_keys_htable_size, strlen("DS"), 0, &tmp_alloc_end, fif_keys_mutable, fif_keys_htable, &fif_key_ct, &dosage_key_idx);
        if (unlikely(reterr)) {
          goto ExportBcf_ret_1;
        }
        write_iter = u32toa(dosage_key_idx, write_iter);
        write_iter = strcpya_k(write_iter, ">" EOLN_STR);
        if (write_hds) {
          // Note that HDS ploidy intentionally does NOT match GT ploidy in the
          // unphased het haploid case.
          write_iter = strcpya_k(write_iter, "##FORMAT=<ID=HDS,Number=.,Type=Float,Description=\"Estimated Haploid Alternate Allele Dosage \",IDX=");
          reterr = AddToFifHtable(g_bigstack_base, "HDS", fif_keys_htable_size, strlen("HDS"), 0, &tmp_alloc_end, fif_keys_mutable, fif_keys_htable, &fif_key_ct, &hds_key_idx);
          if (unlikely(reterr)) {
            goto ExportBcf_ret_1;
          }
          write_iter = u32toa(hds_key_idx, write_iter);
          write_iter = strcpya_k(write_iter, ">" EOLN_STR);
        }
      } else if (write_some_dosage) {
        write_iter = strcpya_k(write_iter, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 \",IDX=");
        reterr = AddToFifHtable(g_bigstack_base, "GP", fif_keys_htable_size, strlen("GP"), 0, &tmp_alloc_end, fif_keys_mutable, fif_keys_htable, &fif_key_ct, &dosage_key_idx);
        if (unlikely(reterr)) {
          goto ExportBcf_ret_1;
        }
        write_iter = u32toa(dosage_key_idx, write_iter);
        write_iter = strcpya_k(write_iter, ">" EOLN_STR);
      }
      // possible todo: optionally export .psam information as
      // PEDIGREE/META/SAMPLE lines in header, and make --vcf/--bcf be able to
      // read it
      if (!ds_only) {
        write_iter = strcpya_k(write_iter, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\",IDX=");
        reterr = AddToFifHtable(g_bigstack_base, "GT", fif_keys_htable_size, strlen("GT"), 0, &tmp_alloc_end, fif_keys_mutable, fif_keys_htable, &fif_key_ct, &gt_key_idx);
        if (unlikely(reterr)) {
          goto ExportBcf_ret_1;
        }
        write_iter = u32toa(gt_key_idx, write_iter);
        write_iter = strcpya_k(write_iter, ">" EOLN_STR);
      }
      write_iter = strcpya_k(write_iter, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
      BigstackEndSet(tmp_alloc_end);

      char* exported_sample_ids;
      uintptr_t max_exported_sample_id_blen;
      if (unlikely(ExportIdpaste(sample_include, siip, "bcf", "bcf", sample_ct, exportf_id_paste, exportf_id_delim, &max_exported_sample_id_blen, &exported_sample_ids))) {
        goto ExportBcf_ret_NOMEM;
      }
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        *write_iter++ = '\t';
        write_iter = strcpya(write_iter, &(exported_sample_ids[sample_idx * max_exported_sample_id_blen]));
      }
      AppendBinaryEoln(&write_iter);
      *write_iter++ = '\0';
      assert(S_CAST(uintptr_t, write_iter - header_start) <= header_ubound);
#ifdef __LP64__
      if (unlikely(S_CAST(uintptr_t, write_iter - header_start) > 0xffffffffU)) {
        logerrputs("Error: Cannot export BCF, since header would be too large.\n");
        goto ExportBcf_ret_INCONSISTENT_INPUT;
      }
#endif
      const uint32_t l_header = write_iter - header_start;
      memcpy(&(writebuf[5]), &l_header, sizeof(int32_t));
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
      goto ExportBcf_ret_WRITE_FAIL;
    }
    // max_allele_ct == 2:
    //   GT: 2 bytes
    //   GP: 3 floats = 12 bytes
    //   DS + HDS: 3 floats = 12 bytes
    //   DS only: 1 float = 4 bytes
    // max_allele_ct > 2:
    //   GT: 2 bytes if max_allele_ct < 64, 4 bytes otherwise
    //   GP: if requested, check max_gp_allele_ct among multiallelic-dosage
    //         variants.
    //       if none found (always true for now), 12 bytes.
    //       if present, (max_gp_allele_ct * (max_gp_allele_ct + 1) / 2) * 8,
    //         which sucks ass
    //   DS + HDS: 3 * (max_allele_ct - 1) floats
    //   DS only: (max_allele_ct - 1) floats
    const uint32_t max_allele_ct = PgrGetMaxAlleleCt(simple_pgrp);
    uintptr_t output_bytes_per_sample;
    if (!allele_idx_offsets) {
      if (write_some_dosage) {
        if (ds_only) {
          output_bytes_per_sample = 4;
        } else {
          output_bytes_per_sample = (write_ds && (!write_hds))? 6 : 14;
        }
      } else {
        output_bytes_per_sample = 2;
      }
    } else {
      if (ds_only) {
        output_bytes_per_sample = 4 * (max_allele_ct - 1);
      } else {
        output_bytes_per_sample = 2 + 2 * (max_allele_ct > 63);
        if (write_some_dosage) {
          if (write_ds) {
            if (write_hds) {
              output_bytes_per_sample += 12 * (max_allele_ct - 1);
            } else {
              output_bytes_per_sample += 4 * (max_allele_ct - 1);
            }
          } else {
            // GP should only be written for biallelic variants, for now
            // would need to take max of this and GT length, but this is always
            // larger
            output_bytes_per_sample = 14;
          }
        }
      }
    }
    // FILTER entries past the first always require 2+ .pvar bytes, and never
    // require more than 4 output bytes.
    // INFO worst case is also a 2:4 ratio.
    uint64_t writebuf_blen = max_chr_slen + S_CAST(uint64_t, kMaxIdSlen) + (max_filter_slen + info_reload_slen) * 2 + S_CAST(uint64_t, sample_ct) * output_bytes_per_sample + 128;
    {
      // Also need to fit all alleles simultaneously; we can't follow our usual
      // strategy of flushing after every allele.
      uintptr_t max_allele_length_bytes = 1;
      if (max_allele_slen > 14) {
        // Need 2-5 additional bytes to encode overflow-size.
        if (max_allele_slen < 128) {
          max_allele_length_bytes = 3;
        } else if (max_allele_slen < 32768) {
          max_allele_length_bytes = 4;
        } else {
          max_allele_length_bytes = 6;
        }
      }
      uint64_t allele_ubound = (S_CAST(uint64_t, max_allele_slen) + max_allele_length_bytes) * max_allele_ct;
      if (allele_ubound > 1000000000) {
        // Compute a tighter bound.  Could multithread this if it ever comes up
        // regularly.
        allele_ubound = 0;
        uintptr_t variant_uidx_base = 0;
        uintptr_t cur_bits = variant_include[0];
        uint32_t allele_ct = 2;
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          uintptr_t allele_idx_offset_base = variant_uidx * 2;
          if (allele_idx_offsets) {
            allele_idx_offset_base = allele_idx_offsets[variant_uidx];
            allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
          }
          const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
          uint64_t cur_ubound = allele_ct * max_allele_length_bytes;
          for (uint32_t allele_idx = 0; allele_idx != allele_ct; ++allele_idx) {
            cur_ubound += strlen(cur_alleles[allele_idx]);
          }
          if (cur_ubound > allele_ubound) {
            allele_ubound = cur_ubound;
          }
        }
      }
      writebuf_blen += allele_ubound;
    }
    writebuf_blen += kMaxMediumLine;
    BigstackReset(writebuf);
    if (unlikely(bigstack_alloc64_c(writebuf_blen, &writebuf))) {
      goto ExportBcf_ret_NOMEM;
    }

    logprintfww5("--export bcf to %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    PgenVariant pgv;
    // a few branches assume genovec is allocated up to word-pair boundary, and
    // that extra word is zeroed out if there is one
    if (unlikely(bigstack_alloc_w(sample_ctl * 2, &pgv.genovec))) {
      goto ExportBcf_ret_NOMEM;
    }
    pgv.genovec[sample_ctl * 2 - 1] = 0;
    pgv.patch_01_set = nullptr;
    pgv.patch_01_vals = nullptr;
    pgv.patch_10_set = nullptr;
    pgv.patch_10_vals = nullptr;
    if (allele_idx_offsets) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &(pgv.patch_01_set)) ||
                   bigstack_alloc_ac(sample_ct, &(pgv.patch_01_vals)) ||
                   bigstack_alloc_w(sample_ctl, &(pgv.patch_10_set)) ||
                   bigstack_alloc_ac(sample_ct * 2, &(pgv.patch_10_vals)))) {
        goto ExportBcf_ret_NOMEM;
      }
    }
    // See ExportVcf() common on how prev_phased works.
    const uint32_t load_dphase = write_hds && (pgfip->gflags & kfPgenGlobalDosagePhasePresent);
    const uint32_t some_phased = (!ds_only) && ((pgfip->gflags & kfPgenGlobalHardcallPhasePresent) || load_dphase);
    uintptr_t* prev_phased = nullptr;
    pgv.phasepresent = nullptr;
    pgv.phaseinfo = nullptr;
    if (some_phased) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &prev_phased) ||
                   bigstack_alloc_w(sample_ctl, &(pgv.phasepresent)) ||
                   bigstack_alloc_w(sample_ctl, &(pgv.phaseinfo)))) {
        goto ExportBcf_ret_NOMEM;
      }
      SetAllBits(sample_ct, prev_phased);
    }

    pgv.dosage_present = nullptr;
    pgv.dosage_main = nullptr;
    pgv.dphase_present = nullptr;
    pgv.dphase_delta = nullptr;
    if (write_some_dosage && (pgfip->gflags & kfPgenGlobalDosagePresent)) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &(pgv.dosage_present)) ||
                   bigstack_alloc_dosage(sample_ct, &(pgv.dosage_main)))) {
        goto ExportBcf_ret_NOMEM;
      }
      if (load_dphase) {
        if (unlikely(bigstack_alloc_w(sample_ctl, &(pgv.dphase_present)) ||
                     bigstack_alloc_dphase(sample_ct, &(pgv.dphase_delta)))) {
          goto ExportBcf_ret_NOMEM;
        }
      }
    }
    const uint32_t fif_key_ctl = BitCtToWordCt(fif_key_ct);
    uintptr_t* fif_seen;
    if (unlikely(bigstack_alloc_w(fif_key_ctl, &fif_seen))) {
      goto ExportBcf_ret_NOMEM;
    }

    // assumes little-endian
    // maybe want to pass references to these arrays later, but leave as
    // C-style for now
    // 4-genotype-at-a-time lookup in basic case
    // these lookup tables are getting large enough to belong in workspace
    // rather than stack memory.
    uint16_t* basic_genobytes4;
    unsigned char* haploid_genobytes4;
    uint16_t* haploid_hh_genobytes4;
    uint16_t* phased_genobytes2;
    uint16_t* phased_hh_genobytes2;
    uint32_t* ds_genobytes4;
    uint32_t* ds_haploid_genobytes4;
    uint32_t* hds_haploid_genobytes4;
    if (unlikely(bigstack_alloc_u16(1024, &basic_genobytes4) ||
                 bigstack_alloc_uc(kLookup256x1bx4Size, &haploid_genobytes4) ||
                 bigstack_alloc_u16(1024, &haploid_hh_genobytes4) ||
                 bigstack_alloc_u16(492, &phased_genobytes2) ||
                 bigstack_alloc_u16(492, &phased_hh_genobytes2) ||
                 bigstack_alloc_u32(1024, &ds_genobytes4) ||
                 bigstack_alloc_u32(1024, &ds_haploid_genobytes4) ||
                 bigstack_alloc_u32(1024, &hds_haploid_genobytes4))) {
      goto ExportBcf_ret_NOMEM;
    }
    // 0/0  0/1  1/1  ./.
    basic_genobytes4[0] = 0x0202;
    basic_genobytes4[4] = 0x0402;
    basic_genobytes4[8] = 0x0404;
    basic_genobytes4[12] = 0;
    InitLookup256x2bx4(basic_genobytes4);

    haploid_genobytes4[0] = 2;
    haploid_genobytes4[4] = 0;
    haploid_genobytes4[8] = 4;
    haploid_genobytes4[12] = 0;
    InitLookup256x1bx4(haploid_genobytes4);

    haploid_hh_genobytes4[0] = 0x8102;
    haploid_hh_genobytes4[4] = 0x0402;
    haploid_hh_genobytes4[8] = 0x8104;
    haploid_hh_genobytes4[12] = 0x8100;
    InitLookup256x2bx4(haploid_hh_genobytes4);

    phased_genobytes2[0] = 0x0202;
    phased_genobytes2[2] = 0x0402;
    phased_genobytes2[4] = 0x0404;
    phased_genobytes2[6] = 0;
    phased_genobytes2[32] = 0x0302;
    phased_genobytes2[34] = 0x0502;
    phased_genobytes2[36] = 0x0504;
    phased_genobytes2[38] = 0x100;
    phased_genobytes2[162] = 0x0304;
    InitVcfPhaseLookup2b(phased_genobytes2);

    phased_hh_genobytes2[0] = 0x8102;
    phased_hh_genobytes2[2] = 0x0402;
    phased_hh_genobytes2[4] = 0x8104;
    phased_hh_genobytes2[6] = 0x8100;
    phased_hh_genobytes2[32] = 0x8102;
    phased_hh_genobytes2[34] = 0x0502;
    phased_hh_genobytes2[36] = 0x8104;
    phased_hh_genobytes2[38] = 0x8100;
    phased_hh_genobytes2[162] = 0x0304;
    InitVcfPhaseLookup2b(phased_hh_genobytes2);

    // don't bother exporting GP for hardcalls
    // usually don't bother for DS, but DS-force is an exception
    ds_genobytes4[0] = 0;
    ds_genobytes4[4] = 0x3f800000;  // float 1.0
    ds_genobytes4[8] = 0x40000000;  // float 2.0
    // in ds_force case, use MISSING instead of END_OF_VECTOR
    ds_genobytes4[12] = 0x7f800001;
    InitLookup256x4bx4(ds_genobytes4);

    // N.B. het-haps are interpreted as ploidy-2, so DS is 1 instead of 0.5.
    ds_haploid_genobytes4[0] = 0;
    ds_haploid_genobytes4[4] = 0x3f800000;
    ds_haploid_genobytes4[8] = 0x3f800000;
    ds_haploid_genobytes4[12] = 0x7f800001;
    InitLookup256x4bx4(ds_haploid_genobytes4);

    // HDS has no special treatment of het-haploids.
    hds_haploid_genobytes4[0] = 0;
    hds_haploid_genobytes4[4] = 0x3f000000;  // float 0.5
    hds_haploid_genobytes4[8] = 0x3f800000;
    hds_haploid_genobytes4[12] = 0x7f800001;
    InitLookup256x4bx4(hds_haploid_genobytes4);

    uint32_t* wide_genobytes4 = nullptr;
    uint16_t* wide_haploid_genobytes4 = nullptr;
    if (max_allele_ct > 63) {
      if (unlikely(bigstack_alloc_u32(1024, &wide_genobytes4) ||
                   bigstack_alloc_u16(1024, &wide_haploid_genobytes4))) {
        goto ExportBcf_ret_NOMEM;
      }
      wide_genobytes4[0] = 0x20002;
      wide_genobytes4[4] = 0x40002;
      wide_genobytes4[8] = 0x40004;
      wide_genobytes4[12] = 0;
      InitLookup256x4bx4(wide_genobytes4);

      wide_haploid_genobytes4[0] = 2;
      wide_haploid_genobytes4[4] = 0;
      wide_haploid_genobytes4[8] = 4;
      wide_haploid_genobytes4[12] = 0;
      InitLookup256x2bx4(wide_haploid_genobytes4);
    }

    int32_t* info_int_buf = nullptr;
    char* pvar_reload_line_iter = nullptr;
    uint32_t info_col_idx = 0;
    if (pvar_info_reload) {
      // Defend against the worst case for now.  Could parse INFO-header-line
      // Number fields and enforce associated limitations in the future, but
      // we'd still use this integer-buffer size if there are any
      // Type=Integer,Number=. declarations.
      if (unlikely(bigstack_alloc_i32((info_reload_slen + 1) / 2, &info_int_buf))) {
        goto ExportBcf_ret_NOMEM;
      }
      reterr = PvarInfoOpenAndReloadHeader(pvar_info_reload, 1 + (max_thread_ct > 1), &pvar_reload_txs, &pvar_reload_line_iter, &info_col_idx);
      if (unlikely(reterr)) {
        goto ExportBcf_ret_TSTREAM_FAIL;
      }
    }

    uint32_t ds_x_genobytes2[128];
    ds_x_genobytes2[0] = 0;
    ds_x_genobytes2[2] = 0x3f800000;
    ds_x_genobytes2[4] = 0x40000000;
    ds_x_genobytes2[6] = 0x7f800001;
    ds_x_genobytes2[32] = 0;
    ds_x_genobytes2[34] = 0x3f800000;  // male het-hap
    ds_x_genobytes2[36] = 0x3f800000;
    ds_x_genobytes2[38] = 0x7f800001;
    // this performs the correct initialization for GenoarrSexLookup4b(), even
    // though phase is not relevant
    InitPhaseXNohhLookup4b(ds_x_genobytes2);

    uint32_t ds_y_genobytes2[128];
    ds_y_genobytes2[0] = 0;
    ds_y_genobytes2[2] = 0x3f800000;
    ds_y_genobytes2[4] = 0x3f800000;
    ds_y_genobytes2[6] = 0x7f800001;
    ds_y_genobytes2[32] = 0;
    ds_y_genobytes2[34] = 0x3f800000;
    ds_y_genobytes2[36] = 0x3f800000;
    ds_y_genobytes2[38] = 0x7f800002;
    InitPhaseXNohhLookup4b(ds_y_genobytes2);

    // HDS ploidy should be more trustworthy than GT ploidy when there's a
    // conflict, since GT must render unphased het haploids as 0/1.  So under
    // HDS-force, diploid missing value is '.,.' instead of '.'.  (Since we're
    // only interested in storing a trustworthy ploidy, we don't add more
    // missing entries in the multiallelic case.)

    uint64_t hds_genobytes2[112] ALIGNV16;
    hds_genobytes2[0] = 0;
    hds_genobytes2[2] = 0x3f0000003f000000LLU;   // 0.5,0.5
    hds_genobytes2[4] = 0x3f8000003f800000LLU;   // 1,1
    hds_genobytes2[6] = 0x7f8000017f800001LLU;   // .,.
    hds_genobytes2[34] = 0x3f80000000000000LLU;  // 0,1
    hds_genobytes2[38] = 0x3f800000;             // 1,0
    InitPhaseLookup8b(hds_genobytes2);

    uint64_t hds_unphased_x_genobytes2[128] ALIGNV16;
    hds_unphased_x_genobytes2[0] = 0;
    hds_unphased_x_genobytes2[2] = 0x3f0000003f000000LLU;  // 0.5,0.5
    hds_unphased_x_genobytes2[4] = 0x3f8000003f800000LLU;  // 1,1
    hds_unphased_x_genobytes2[6] = 0x7f8000017f800001LLU;  // .,.
    hds_unphased_x_genobytes2[32] = 0x7f80000200000000LLU;  // 0
    hds_unphased_x_genobytes2[34] = 0x7f8000023f000000LLU;  // 0.5
    hds_unphased_x_genobytes2[36] = 0x7f8000023f800000LLU;  // 1
    hds_unphased_x_genobytes2[38] = 0x7f8000027f800001LLU;  // .
    InitPhaseXNohhLookup8b(hds_unphased_x_genobytes2);

    uint64_t hds_hh_genobytes2[112] ALIGNV16;
    hds_hh_genobytes2[0] = 0x7f80000200000000LLU;
    hds_hh_genobytes2[2] = 0x7f8000023f000000LLU;   // 0.5
    hds_hh_genobytes2[4] = 0x7f8000023f800000LLU;   // 1
    hds_hh_genobytes2[6] = 0x7f8000027f800001LLU;   // .
    hds_hh_genobytes2[34] = 0x3f80000000000000LLU;  // 0,1
    hds_hh_genobytes2[38] = 0x3f800000;             // 1,0
    InitPhaseLookup8b(hds_hh_genobytes2);

    uint32_t hds_y_genobytes2[128];
    hds_y_genobytes2[0] = 0;
    hds_y_genobytes2[2] = 0x3f000000;
    hds_y_genobytes2[4] = 0x3f800000;
    hds_y_genobytes2[6] = 0x7f800001;
    hds_y_genobytes2[32] = 0;
    hds_y_genobytes2[34] = 0x3f000000;
    hds_y_genobytes2[36] = 0x3f800000;
    hds_y_genobytes2[38] = 0x7f800002;
    InitPhaseXNohhLookup4b(hds_y_genobytes2);

    const char* const* fif_keys = TO_CONSTCPCONSTP(fif_keys_mutable);
    const uint32_t inf_bits[2] = {0x7f800000, 0xff800000U};
    const char* dot_ptr = &(g_one_char_strs[92]);
    const uint32_t male_ct = PopcountWords(sex_male_collapsed, sample_ctl);
    const uint32_t female_ct = PopcountWords(sex_female_collapsed, sample_ctl);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_bcf_idx = 0;
    uint32_t is_x = 0;
    uint32_t is_y = 0;
    uint32_t is_haploid = 0;  // includes chrX and chrY
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    uint32_t rls_variant_uidx = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uint32_t allele_ct = 2;
    uint32_t invalid_allele_code_seen = 0;
    uint32_t first_rlen_warning_idxs[3];
    uint32_t rlen_warning_ct = 0;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        is_x = 0;
        is_y = 0;
        is_haploid = IsSet(cip->haploid_mask, chr_idx);
        if (chr_idx == cip->xymt_codes[kChrOffsetX]) {
          if (!male_ct) {
            // No male chrX ploidy fix needed if there are no males.
            is_haploid = 0;
          } else if (male_ct != sample_ct) {
            // Can use generic haploid encoder if all males.
            is_x = 1;
          }
        } else if ((chr_idx == cip->xymt_codes[kChrOffsetY]) && female_ct) {
          // No female chrY ploidy fix needed if there are no females.
          is_y = 1;
        }
        chr_bcf_idx = chr_fo_to_bcf_idx[chr_fo_idx];
      }
      char* rec_start = write_iter;
      // first 8 bytes are l_shared and l_indiv, fill them later
      memcpy(&(rec_start[8]), &chr_bcf_idx, sizeof(int32_t));  // CHROM
      const int32_t bp0 = S_CAST(int32_t, variant_bps[variant_uidx] - 1);
      memcpy(&(rec_start[12]), &bp0, sizeof(int32_t));  // POS
      // [16,20) is rlen, which may depend on INFO/END... er, VCF 4.5 spec says
      // INFO/SVLEN and FORMAT/LEN should be consulted instead

      if ((!pvar_quals) || (!IsSet(pvar_qual_present, variant_uidx))) {
        const uint32_t missing_val = 0x7f800001;
        memcpy(&(rec_start[20]), &missing_val, 4);
      } else {
        memcpy(&(rec_start[20]), &(pvar_quals[variant_uidx]), 4);
      }

      // [24,26) is n_info
      // [26,28) is n_allele
      memcpy(&(rec_start[28]), &sample_ct, sizeof(int32_t));
      // [31] is n_fmt
      write_iter = &(rec_start[32]);

      // ID
      const char* variant_id = variant_ids[variant_uidx];
      if (strequal_k_unsafe(variant_id, ".")) {
        *write_iter++ = '\7';
      } else {
        write_iter = AppendBcfString(variant_id, write_iter);
      }

      // REF, ALT
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (allele_permute) {
        ref_allele_idx = allele_permute[allele_idx_offset_base];
        alt1_allele_idx = allele_permute[allele_idx_offset_base + 1];
      }
      const char* ref_allele = cur_alleles[ref_allele_idx];
      uint32_t rlen = strlen(ref_allele);
      if (ref_allele != dot_ptr) {
        write_iter = AppendBcfCountedString(ref_allele, rlen, write_iter);
        if (!invalid_allele_code_seen) {
          invalid_allele_code_seen = !ValidVcfAlleleCode(ref_allele);
        }
      } else {
        write_iter = AppendBcfString("N", write_iter);
      }
      const char* alt1_allele = cur_alleles[alt1_allele_idx];
      uint16_t n_allele = allele_ct;
      if (alt1_allele != dot_ptr) {
        write_iter = AppendBcfString(alt1_allele, write_iter);
        if (!invalid_allele_code_seen) {
          invalid_allele_code_seen = !ValidVcfAlleleCode(alt1_allele);
        }
        if (allele_ct > 2) {
          for (uint32_t cur_allele_uidx = 0; cur_allele_uidx != allele_ct; ++cur_allele_uidx) {
            if ((cur_allele_uidx == ref_allele_idx) || (cur_allele_uidx == alt1_allele_idx)) {
              // if this is noticeably suboptimal, have two loops, with inner
              // loop going up to cur_allele_stop.
              // (also wrap this in a function, this comes up a bunch of times)
              continue;
            }
            const char* cur_allele = cur_alleles[cur_allele_uidx];
            write_iter = AppendBcfString(cur_allele, write_iter);
            if (!invalid_allele_code_seen) {
              invalid_allele_code_seen = !ValidVcfAlleleCode(cur_allele);
            }
          }
        }
      } else {
        if (unlikely(allele_ct > 2)) {
          snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: Multiallelic variant '%s' has a missing ALT allele.\n", variant_id);
          goto ExportBcf_ret_MALFORMED_INPUT_WW;
        }
        n_allele = 1;
      }

      // QUAL already filled

      // FILTER
      if ((!pvar_filter_present) || (!IsSet(pvar_filter_present, variant_uidx))) {
        *write_iter++ = '\0';
      } else if (!IsSet(pvar_filter_npass, variant_uidx)) {
        *write_iter++ = 0x11;
        *write_iter++ = '\0';
      } else {
        const char* cur_filter_iter = pvar_filter_storage[variant_uidx];
        // two passes; determine type descriptor during first pass (and may as
        // well deduplicate and sort while we're at it)
        ZeroWArr(fif_key_ctl, fif_seen);
        while (1) {
          // bugfix (3 Jan 2021): this is semicolon, not comma, delimited
          const char* tok_end = Strchrnul(cur_filter_iter, ';');
          const uint32_t fif_idx = IdHtableFindNnt(cur_filter_iter, fif_keys, fif_keys_htable, tok_end - cur_filter_iter, fif_keys_htable_size);
          // Second predicate verifies this is actually a FILTER key.
          if (unlikely((fif_idx == UINT32_MAX) || (!(fif_keys[fif_idx][-1] & 8)))) {
            snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: Variant '%s' has a FILTER value with no corresponding header line.\n", variant_id);
            goto ExportBcf_ret_MALFORMED_INPUT_WW;
          }
          SetBit(fif_idx, fif_seen);
          if ((*tok_end) == '\0') {
            break;
          }
          cur_filter_iter = &(tok_end[1]);
        }
        const uint32_t nfilter = PopcountWords(fif_seen, fif_key_ctl);
        const uint32_t max_idx = FindLast1BitBefore(fif_seen, fif_key_ct);
        uintptr_t fif_idx_base = 0;
        uintptr_t fif_seen_bits = fif_seen[0];
        if (max_idx <= 127) {
          write_iter = AppendBcfVecType(1, nfilter, write_iter);
          for (uint32_t uii = 0; uii != nfilter; ++uii) {
            *write_iter++ = BitIter1(fif_seen, &fif_idx_base, &fif_seen_bits);
          }
        } else if (max_idx <= 32767) {
          write_iter = AppendBcfVecType(2, nfilter, write_iter);
          for (uint32_t uii = 0; uii != nfilter; ++uii) {
            const uint16_t fif_idx = BitIter1(fif_seen, &fif_idx_base, &fif_seen_bits);
            memcpy(write_iter, &fif_idx, sizeof(int16_t));
            write_iter += sizeof(int16_t);
          }
        } else {
          write_iter = AppendBcfVecType(3, nfilter, write_iter);
          for (uint32_t uii = 0; uii != nfilter; ++uii) {
            const uint32_t fif_idx = BitIter1(fif_seen, &fif_idx_base, &fif_seen_bits);
            memcpy(write_iter, &fif_idx, sizeof(int32_t));
            write_iter += sizeof(int32_t);
          }
        }
      }

      // INFO
      const uint32_t is_pr = all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx));
      uint16_t n_info = 0;
      if (pvar_reload_line_iter) {
        reterr = PvarInfoReload(info_col_idx, variant_uidx, &pvar_reload_txs, &pvar_reload_line_iter, &rls_variant_uidx);
        if (unlikely(reterr)) {
          goto ExportBcf_ret_TSTREAM_FAIL;
        }
        // BCF translation of PvarInfoWrite().
        char* info_iter = pvar_reload_line_iter;
        char* info_end = CurTokenEnd(info_iter);
        uint32_t n_info_raw = 0;
        ZeroWArr(fif_key_ctl, fif_seen);
        if ((*info_iter != '.') || (&(info_iter[1]) != info_end)) {
          const char orig_info_end_char = *info_end;
          *info_end = ';';
          char* info_stop = &(info_end[1]);
          do {
            char* key_end = S_CAST(char*, rawmemchr2(info_iter, '=', ';'));
            const uint32_t fif_idx = IdHtableFindNnt(info_iter, fif_keys, fif_keys_htable, key_end - info_iter, fif_keys_htable_size);
            // Second predicate verifies this is actually an INFO key.
            uint32_t info_type_code = 0;
            if (fif_idx != UINT32_MAX) {
              info_type_code = ctou32(fif_keys[fif_idx][-1]) & 7;
            }
            if (unlikely(!(info_type_code & 1))) {
              snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: Variant '%s' has an INFO key with no corresponding header line.\n", variant_id);
              goto ExportBcf_ret_MALFORMED_INPUT_WW;
            }
            if (unlikely(IsSet(fif_seen, fif_idx))) {
              snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: Variant '%s' has multiple INFO/%s entries.\n", variant_id, fif_keys[fif_idx]);
              goto ExportBcf_ret_MALFORMED_INPUT_WW;
            }
            SetBit(fif_idx, fif_seen);
            ++n_info_raw;
            write_iter = AppendBcfTypedInt(fif_idx, write_iter);
            if (*key_end == ';') {
              // bugfix (8 Jan 2023): this previously errored out on non-flag
              // keys, but VCFv4.3 (though not v4.2) allows non-flag keys to
              // have no value.

              // Anything goes here.  As of this writing, the spec suggests
              // 0x11 0x01 (single int8 equal to 1), but bcftools emits the
              // more compact 0x00 (typeless missing value); we imitate
              // bcftools.
              *write_iter++ = '\0';
              info_iter = &(key_end[1]);
              continue;
            }
            char* value_start = &(key_end[1]);
            char* value_end = S_CAST(char*, rawmemchr(value_start, ';'));
            if (info_type_code == 7) {
              // No need to count commas in this case.
              write_iter = AppendBcfCountedString(value_start, value_end - value_start, write_iter);
            } else {
              const uint32_t value_ct = 1 + CountByte(value_start, ',', value_end - value_start);
              *value_end = ',';
              if (info_type_code == 3) {
                // Integer
                // We don't know upfront whether int8, int16, or int32 is best.
                // We could use two passes, but ScanmovIntBounded32() is slower
                // than the other steps here.  So we instead decode to an
                // int32[] buffer while tracking whether any numbers are
                // actually out of int8/int16 range, and then convert to the
                // final encoding afterward.
                const char* value_iter = value_start;
                uint32_t encode_as_small_ints = 2;  // 2 = int8, 1 = int16
                for (uint32_t value_idx = 0; value_idx != value_ct; ++value_idx) {
                  int32_t cur_val;
                  if (*value_iter == '.') {
                    cur_val = -2147483647 - 1;  // missing
                    ++value_iter;
                  } else {
                    if (unlikely(ScanmovIntBounded(0x7ffffff8, 0x7fffffff, &value_iter, &cur_val))) {
                      snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: Variant '%s' has an invalid INFO/%s value.\n", variant_id, fif_keys[fif_idx]);
                      goto ExportBcf_ret_MALFORMED_INPUT_WW;
                    }
                    if (encode_as_small_ints) {
                      if ((cur_val < -32760) || (cur_val > 32767)) {
                        encode_as_small_ints = 0;
                      } else if ((encode_as_small_ints == 2) && ((cur_val < -120) || (cur_val > 127))) {
                        encode_as_small_ints = 1;
                      }
                    }
                  }
                  if (unlikely(*value_iter != ',')) {
                    snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: Variant '%s' has an invalid INFO/%s value.\n", variant_id, fif_keys[fif_idx]);
                    goto ExportBcf_ret_MALFORMED_INPUT_WW;
                  }
                  info_int_buf[value_idx] = cur_val;
                  ++value_iter;
                }
                write_iter = AppendBcfVecType(3 - encode_as_small_ints, value_ct, write_iter);
                if (!encode_as_small_ints) {
                  // int32
                  write_iter = memcpya(write_iter, info_int_buf, value_ct * sizeof(int32_t));
                } else if (encode_as_small_ints == 1) {
                  // int16
                  for (uint32_t value_idx = 0; value_idx != value_ct; ++value_idx) {
                    int32_t cur_val = info_int_buf[value_idx];
                    if (cur_val == (-2147483647 - 1)) {
                      cur_val = -32768;
                    }
                    const int16_t cur_i16 = cur_val;
                    CopyToUnalignedOffsetI16(CToUc(write_iter), &cur_i16, value_idx);
                  }
                  write_iter = &(write_iter[value_ct * sizeof(int16_t)]);
                } else {
                  // int8
                  for (uint32_t value_idx = 0; value_idx != value_ct; ++value_idx) {
                    int32_t cur_val = info_int_buf[value_idx];
                    if (cur_val == (-2147483647 - 1)) {
                      cur_val = -128;
                    }
                    write_iter[value_idx] = cur_val;
                  }
                  write_iter = &(write_iter[value_ct]);
                }
                if ((fif_idx == end_key_idx) && (value_ct == 1) && (info_int_buf[0] != (-2147483647 - 1))) {
                  // Only save INFO/END-derived rlen when REF is a single base,
                  // and info_end_rlen is positive.
                  // Otherwise, print a warning if the two rlens aren't equal.
                  // TODO: process INFO/SVLEN and FORMAT/LEN as directed by the
                  // VCF 4.5 specification instead.
                  const uint32_t info_end_rlen = 1 + info_int_buf[0] - variant_bps[variant_uidx];
                  if ((rlen == 1) && (S_CAST(int32_t, info_end_rlen) > 0)) {
                    rlen = info_end_rlen;
                  } else if (rlen != info_end_rlen) {
                    MakeRlenWarningStr(variant_id);
                    if (!rlen_warning_ct) {
                      logputs_silent("\n");
                    }
                    logputs_silent(g_logbuf);
                    if (rlen_warning_ct < 3) {
                      // Don't print to console now, since that has an annoying
                      // interaction with the progress indicator.
                      first_rlen_warning_idxs[rlen_warning_ct] = variant_uidx;
                    }
                    ++rlen_warning_ct;
                  }
                }
              } else if (likely(info_type_code == 5)) {
                // Float
                // Permit nan and +/-inf, but prohibit floating-point overflow.
                write_iter = AppendBcfVecType(5, value_ct, write_iter);
                const char* value_iter = value_start;
                const uint32_t missing_bits = 0x7f800001;
                const uint32_t nan_bits = 0x7fc00000;
                unsigned char* write_uc = CToUc(write_iter);
                for (uint32_t value_idx = 0; value_idx != value_ct; ++value_idx) {
                  if (*value_iter == '.') {
                    CopyToUnalignedOffsetU32(write_uc, &missing_bits, value_idx);
                    ++value_iter;
                  } else {
                    double dxx;
                    const char* floatstr_end = ScanadvDouble(value_iter, &dxx);
                    if (floatstr_end) {
                      if (unlikely(fabs(dxx) > 3.4028235677973362e38)) {
                        snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: Variant '%s' has an out-of-range INFO/%s value.\n", variant_id, fif_keys[fif_idx]);
                        goto ExportBcf_ret_MALFORMED_INPUT_WW;
                      }
                      const float fxx = S_CAST(float, dxx);
                      CopyToUnalignedOffsetF(write_uc, &fxx, value_idx);
                      value_iter = floatstr_end;
                    } else {
                      const char* next_comma = S_CAST(const char*, rawmemchr(value_iter, ','));
                      const uint32_t slen = next_comma - value_iter;
                      if (IsNanStr(value_iter, slen)) {
                        CopyToUnalignedOffsetU32(write_uc, &nan_bits, value_idx);
                      } else {
                        uint32_t is_neg = 0;
                        if (unlikely(!IsInfStr(value_iter, slen, &is_neg))) {
                          snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: Variant '%s' has an invalid INFO/%s value.\n", variant_id, fif_keys[fif_idx]);
                          goto ExportBcf_ret_MALFORMED_INPUT_WW;
                        }
                        CopyToUnalignedOffsetU32(write_uc, &(inf_bits[is_neg]), value_idx);
                      }
                      value_iter = &(next_comma[1]);
                      continue;
                    }
                  }
                  if (unlikely(*value_iter != ',')) {
                    snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: Variant '%s' has an invalid INFO/%s value.\n", variant_id, fif_keys[fif_idx]);
                    goto ExportBcf_ret_MALFORMED_INPUT_WW;
                  }
                  ++value_iter;
                }
                write_iter = &(write_iter[value_ct * sizeof(float)]);
              } else {
                snprintf(g_logbuf, kLogbufSize, "Error: --export bcf: INFO/%s for variant '%s' has a value, but is of type Flag.\n", fif_keys[fif_idx], variant_id);
                goto ExportBcf_ret_MALFORMED_INPUT_WW;
              }
              // *value_end = ';';  // don't actually need this
            }
            info_iter = &(value_end[1]);
          } while (info_iter != info_stop);
          *info_end = orig_info_end_char;
        }
        if (is_pr && (!IsSet(fif_seen, pr_key_idx))) {
          ++n_info_raw;
          write_iter = AppendBcfTypedInt(pr_key_idx, write_iter);
          *write_iter++ = '\0';
        }
        if (unlikely(n_info_raw > 0xffff)) {
          logputs("\n");
          logerrprintfww("Error: INFO component of variant '%s' has too many keys to represent as BCF.\n", variant_id);
          goto ExportBcf_ret_INCONSISTENT_INPUT;
        }
        n_info = n_info_raw;
      } else if (is_pr) {
        n_info = 1;
        write_iter = AppendBcfTypedInt(pr_key_idx, write_iter);
        *write_iter++ = '\0';
      }

#ifdef __LP64__
      if (unlikely(S_CAST(uintptr_t, write_iter - &(rec_start[8])) > UINT32_MAX)) {
        logputs("\n");
        logerrprintfww("Error: CHROM..INFO component of variant '%s' is too large to represent as BCF.\n", variant_id);
        goto ExportBcf_ret_INCONSISTENT_INPUT;
      }
#endif
      const uint32_t l_shared = write_iter - &(rec_start[8]);
      memcpy(rec_start, &l_shared, sizeof(int32_t));
      memcpy(&(rec_start[16]), &rlen, sizeof(int32_t));
      memcpy(&(rec_start[24]), &n_info, sizeof(int16_t));
      memcpy(&(rec_start[26]), &n_allele, sizeof(int16_t));
      // Simplest for l_indiv [4,8) and n_fmt [31] calculation to wait until
      // the entire variant record is complete.

      // could defensively zero out more counts
      pgv.dosage_ct = 0;
      pgv.dphase_ct = 0;
      char* indiv_start = write_iter;
      uint32_t n_fmt = 0;
      if (!ds_only) {
        write_iter = AppendBcfTypedInt(gt_key_idx, write_iter);
        n_fmt = 1;
      }
      if (allele_ct == 2) {
        if (!some_phased) {
          // biallelic, nothing phased in entire file (or ds_only)
          // (technically possible for dosage-phase to be present, if no
          // hardcalls are phased and HDS output not requested)
          if (!write_some_dosage) {
            reterr = PgrGet(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec);
          } else {
            reterr = PgrGetD(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec, pgv.dosage_present, pgv.dosage_main, &(pgv.dosage_ct));
          }
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ExportBcf_ret_1;
          }
          if (!alt1_allele_idx) {
            GenovecInvertUnsafe(sample_ct, pgv.genovec);
            if (pgv.dosage_ct) {
              BiallelicDosage16Invert(pgv.dosage_ct, pgv.dosage_main);
            }
          }
          if (!ds_only) {
            if ((!is_haploid) || is_x) {
              *write_iter++ = 0x21;
              GenoarrLookup256x2bx4(pgv.genovec, basic_genobytes4, sample_ct, write_iter);
              if (is_x) {
                FixBcfMaleXGtPloidy(pgv.genovec, sex_male_collapsed, sample_ct, write_iter);
              }
              write_iter = &(write_iter[sample_ct * 2]);
            } else {
              // Ploidy 1 unless there's a het-haploid.
              ZeroTrailingNyps(sample_ct, pgv.genovec);
              if (AtLeastOneHetUnsafe(pgv.genovec, sample_ct)) {
                *write_iter++ = 0x21;
                GenoarrLookup256x2bx4(pgv.genovec, haploid_hh_genobytes4, sample_ct, write_iter);
                if (is_y) {
                  FixBcfFemaleYGtPloidy2(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                }
                write_iter = &(write_iter[sample_ct * 2]);
              } else {
                *write_iter++ = 0x11;
                GenoarrLookup256x1bx4(pgv.genovec, haploid_genobytes4, sample_ct, write_iter);
                if (is_y) {
                  FixBcfFemaleYGtPloidy1(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                }
                write_iter = &(write_iter[sample_ct]);
              }
            }
          }
          if (pgv.dosage_ct || ds_force) {
            ++n_fmt;
            write_iter = AppendBcfTypedInt(dosage_key_idx, write_iter);
            if (write_ds) {
              *write_iter++ = 0x15;
              FillBcfDs(&pgv, sex_male_collapsed, sex_female_collapsed, ds_genobytes4, ds_x_genobytes2, ds_y_genobytes2, ds_haploid_genobytes4, sample_ct, ds_force, is_x, is_y, is_haploid, write_iter);
              write_iter = &(write_iter[sample_ct * sizeof(int32_t)]);
              if (hds_force) {
                ++n_fmt;
                write_iter = AppendBcfTypedInt(hds_key_idx, write_iter);
                if (!is_haploid) {
                  *write_iter++ = 0x25;
                  GenoarrLookup16x8bx2(pgv.genovec, hds_genobytes2, sample_ct, write_iter);
                  uint32_t widx = 0;
                  for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ) {
                    uintptr_t dosage_present_word;
                    do {
                      dosage_present_word = pgv.dosage_present[widx++];
                    } while (!dosage_present_word);
                    unsigned char* cur_hds = CToUc(&(write_iter[widx * (kBitsPerWord * 8 * k1LU)]));
                    do {
                      const uint32_t sample_idx_lowbits = ctzw(dosage_present_word);
                      const uint32_t dosage_int = pgv.dosage_main[dosage_idx++];
                      const float dosage_f = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMid);
                      CopyToUnalignedOffsetF(cur_hds, &dosage_f, 2 * sample_idx_lowbits);
                      CopyToUnalignedOffsetF(cur_hds, &dosage_f, 2 * sample_idx_lowbits + 1);
                      dosage_present_word &= dosage_present_word - 1;
                    } while (dosage_present_word);
                  }
                  write_iter = &(write_iter[sample_ct * sizeof(int32_t) * 2]);
                } else if (is_x) {
                  *write_iter++ = 0x25;
                  // male dosages 0..1
                  GenoarrSexLookup8b(pgv.genovec, sex_male_collapsed, hds_unphased_x_genobytes2, sample_ct, write_iter);
                  uint32_t widx = 0;
                  for (uint32_t dosage_idx = 0; dosage_idx != pgv.dosage_ct; ) {
                    uintptr_t dosage_present_word;
                    do {
                      dosage_present_word = pgv.dosage_present[widx++];
                    } while (!dosage_present_word);
                    unsigned char* cur_hds = CToUc(&(write_iter[widx * (kBitsPerWord * 8 * k1LU)]));
                    do {
                      const uint32_t sample_idx_lowbits = ctzw(dosage_present_word);
                      const uint32_t dosage_int = pgv.dosage_main[dosage_idx++];
                      const float dosage_f = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMid);
                      CopyToUnalignedOffsetF(cur_hds, &dosage_f, 2 * sample_idx_lowbits);
                      CopyToUnalignedOffsetF(cur_hds, &dosage_f, 2 * sample_idx_lowbits + 1);
                      dosage_present_word &= dosage_present_word - 1;
                    } while (dosage_present_word);
                  }
                  write_iter = &(write_iter[sample_ct * sizeof(int32_t) * 2]);
                } else {
                  // max ploidy 1.
                  *write_iter++ = 0x15;
                  FillBcfHdsPloidy1(&pgv, sex_female_collapsed, hds_y_genobytes2, hds_haploid_genobytes4, sample_ct, 1, is_y, write_iter);
                  write_iter = &(write_iter[sample_ct * sizeof(int32_t)]);
                }
              }
            } else {
              // GP
              if ((!is_haploid) || is_x) {
                *write_iter++ = 0x35;
                FillBcfGpPloidy2(&pgv, sex_male_collapsed, sample_ct, is_x, write_iter);
                write_iter = &(write_iter[sample_ct * sizeof(int32_t) * 3]);
              } else {
                *write_iter++ = 0x25;
                FillBcfGpPloidy1(&pgv, sample_ct, write_iter);
                write_iter = &(write_iter[sample_ct * sizeof(int32_t) * 2]);
              }
            }
          }
        } else {
          // biallelic, phased
          if (!write_some_dosage) {
            reterr = PgrGetP(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, pgv.genovec, pgv.phasepresent, pgv.phaseinfo, &(pgv.phasepresent_ct));
          } else {
            reterr = PgrGetDp(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
          }
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ExportBcf_ret_1;
          }
          if (!alt1_allele_idx) {
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
          // little point in trying to optimize this case, thanks to
          // prev_phased.  Instead, we might want to have a fast path for the
          // all-phased case.
          if (!pgv.phasepresent_ct) {
            ZeroWArr(sample_ctl, pgv.phasepresent);
          }
          ZeroTrailingNyps(sample_ct, pgv.genovec);
          BitvecAnd(pgv.phasepresent, sample_ctl, pgv.phaseinfo);
          if ((!is_haploid) || is_x) {
            *write_iter++ = 0x21;
            FillBcfPhasedGt(pgv.genovec, pgv.phasepresent, pgv.phaseinfo, phased_genobytes2, sample_ct, prev_phased, write_iter);
            if (is_x) {
              FixBcfMaleXGtPloidy(pgv.genovec, sex_male_collapsed, sample_ct, write_iter);
            }
            write_iter = &(write_iter[sample_ct * 2]);
          } else {
            // Ploidy 1 unless there's a het-haploid.
            if (AtLeastOneHetUnsafe(pgv.genovec, sample_ct)) {
              *write_iter++ = 0x21;
              FillBcfPhasedGt(pgv.genovec, pgv.phasepresent, pgv.phaseinfo, phased_hh_genobytes2, sample_ct, prev_phased, write_iter);
              if (is_y) {
                FixBcfFemaleYGtPloidy2(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
              }
              write_iter = &(write_iter[sample_ct * 2]);
            } else {
              *write_iter++ = 0x11;
              GenoarrLookup256x1bx4(pgv.genovec, haploid_genobytes4, sample_ct, write_iter);
              if (is_y) {
                FixBcfFemaleYGtPloidy1(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
              }
              write_iter = &(write_iter[sample_ct]);
            }
          }
          if (pgv.dosage_ct || ds_force) {
            // both dosage (or {H}DS-force) and phase present
            ++n_fmt;
            write_iter = AppendBcfTypedInt(dosage_key_idx, write_iter);
            if (write_ds) {
              *write_iter++ = 0x15;
              FillBcfDs(&pgv, sex_male_collapsed, sex_female_collapsed, ds_genobytes4, ds_x_genobytes2, ds_y_genobytes2, ds_haploid_genobytes4, sample_ct, ds_force, is_x, is_y, is_haploid, write_iter);
              write_iter = &(write_iter[sample_ct * sizeof(int32_t)]);
              const uint32_t dphase_ct = pgv.dphase_ct;
              if (hds_force || dphase_ct ||
                  (write_hds && pgv.phasepresent_ct && pgv.dosage_ct && (!IntersectionIsEmpty(pgv.phasepresent, pgv.dosage_present, sample_ctl)))) {
                ++n_fmt;
                write_iter = AppendBcfTypedInt(hds_key_idx, write_iter);
                // dphase_present can be nullptr, so we zero-initialize
                // dphase_present_word and never refresh it when dphase_ct==0.
                uintptr_t dphase_present_word = 0;
                if ((!is_haploid) || is_x) {
                  *write_iter++ = 0x25;
                  if (!hds_force) {
                    const uint64_t end_of_vector = 0x7f8000027f800002LLU;
                    for (uint32_t uii = 0; uii != sample_ct; ++uii) {
                      CopyToUnalignedOffsetU64(CToUc(write_iter), &end_of_vector, uii);
                    }
                  } else {
                    PhaseLookup8b(pgv.genovec, pgv.phasepresent, pgv.phaseinfo, hds_genobytes2, sample_ct, write_iter);
                  }
                  // Note that dosage_main may be nullptr in HDS-force case.
                  // Guard against undefined behavior (pointer arithmetic on
                  // nullptr).
                  if (pgv.dosage_ct) {
                    const uintptr_t* phasepresent = pgv.phasepresent;
                    const uintptr_t* phaseinfo = pgv.phaseinfo;
                    const uintptr_t* dosage_present = pgv.dosage_present;
                    const Dosage* dosage_main_stop = &(pgv.dosage_main[pgv.dosage_ct]);
                    const uintptr_t* dphase_present = pgv.dphase_present;
                    const SDosage* dphase_delta_iter = pgv.dphase_delta;
                    uint32_t widx = UINT32_MAX;  // deliberate overflow
                    for (const Dosage* dosage_main_iter = pgv.dosage_main; dosage_main_iter != dosage_main_stop; ) {
                      uintptr_t dosage_present_word;
                      do {
                        dosage_present_word = dosage_present[++widx];
                      } while (!dosage_present_word);
                      const uintptr_t phasepresent_word = phasepresent[widx];
                      const uintptr_t phaseinfo_word = phaseinfo[widx];
                      if (dphase_ct) {
                        dphase_present_word = dphase_present[widx];
                      }
                      unsigned char* cur_hds = CToUc(&(write_iter[widx * (kBitsPerWord * 8 * k1LU)]));
                      do {
                        const uint32_t sample_idx_lowbits = ctzw(dosage_present_word);
                        const uintptr_t lowbit = dosage_present_word & (-dosage_present_word);
                        const uint32_t dosage_int = *dosage_main_iter++;
                        const uintptr_t is_phased = (phasepresent_word | dphase_present_word) & lowbit;
                        if ((!hds_force) && (!is_phased)) {
                          dosage_present_word ^= lowbit;
                          continue;
                        }
                        float dosage0 = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMax);
                        float dosage1 = dosage0;
                        if (is_phased) {
                          int32_t cur_dphase_delta;
                          if (dphase_present_word & lowbit) {
                            cur_dphase_delta = *dphase_delta_iter++;
                          } else {
                            cur_dphase_delta = DosageHomdist(dosage_int);
                            if (!(phaseinfo_word & lowbit)) {
                              cur_dphase_delta = -cur_dphase_delta;
                            }
                          }
                          const float dosage0_incr_f = S_CAST(float, cur_dphase_delta) * S_CAST(float, kRecipDosageMax);
                          dosage0 += dosage0_incr_f;
                          dosage1 -= dosage0_incr_f;
                        }
                        CopyToUnalignedOffsetF(cur_hds, &dosage0, 2 * sample_idx_lowbits);
                        CopyToUnalignedOffsetF(cur_hds, &dosage1, 2 * sample_idx_lowbits + 1);
                        dosage_present_word ^= lowbit;
                      } while (dosage_present_word);
                    }
                  }
                  if (is_x) {
                    FixBcfMaleXHdsPloidy(pgv.phasepresent, dphase_ct? pgv.dphase_present : nullptr, sex_male_collapsed, sample_ct, write_iter);
                  }
                  write_iter = &(write_iter[sample_ct * sizeof(int32_t) * 2]);
                } else if (pgv.phasepresent_ct || pgv.dphase_ct) {
                  // Supposed to be all haploid, but we must render as ploidy-2
                  // anyway.
                  // probable todo: merge this with the regular diploid case.
                  *write_iter++ = 0x25;
                  if (!hds_force) {
                    const uint64_t end_of_vector = 0x7f8000027f800002LLU;
                    for (uint32_t uii = 0; uii != sample_ct; ++uii) {
                      CopyToUnalignedOffsetU64(CToUc(write_iter), &end_of_vector, uii);
                    }
                  } else {
                    PhaseLookup8b(pgv.genovec, pgv.phasepresent, pgv.phaseinfo, hds_hh_genobytes2, sample_ct, write_iter);
                  }
                  if (pgv.dosage_ct) {
                    const uintptr_t* phasepresent = pgv.phasepresent;
                    const uintptr_t* phaseinfo = pgv.phaseinfo;
                    const uintptr_t* dosage_present = pgv.dosage_present;
                    const Dosage* dosage_main_stop = &(pgv.dosage_main[pgv.dosage_ct]);
                    const uintptr_t* dphase_present = pgv.dphase_present;
                    const SDosage* dphase_delta_iter = pgv.dphase_delta;
                    const uint32_t eov_bits = 0x7f800002;
                    uint32_t widx = UINT32_MAX;  // deliberate overflow
                    for (const Dosage* dosage_main_iter = pgv.dosage_main; dosage_main_iter != dosage_main_stop; ) {
                      uintptr_t dosage_present_word;
                      do {
                        dosage_present_word = dosage_present[++widx];
                      } while (!dosage_present_word);
                      const uintptr_t phasepresent_word = phasepresent[widx];
                      const uintptr_t phaseinfo_word = phaseinfo[widx];
                      if (dphase_ct) {
                        dphase_present_word = dphase_present[widx];
                      }
                      unsigned char* cur_hds = CToUc(&(write_iter[widx * (kBitsPerWord * 8 * k1LU)]));
                      do {
                        const uint32_t sample_idx_lowbits = ctzw(dosage_present_word);
                        const uintptr_t lowbit = dosage_present_word & (-dosage_present_word);
                        const uint32_t dosage_int = *dosage_main_iter++;
                        const uintptr_t is_phased = (phasepresent_word | dphase_present_word) & lowbit;
                        if ((!hds_force) && (!is_phased)) {
                          dosage_present_word ^= lowbit;
                          continue;
                        }
                        float dosage0 = S_CAST(float, dosage_int) * S_CAST(float, kRecipDosageMax);
                        if (is_phased) {
                          int32_t cur_dphase_delta;
                          if (dphase_present_word & lowbit) {
                            cur_dphase_delta = *dphase_delta_iter++;
                          } else {
                            cur_dphase_delta = DosageHomdist(dosage_int);
                            if (!(phaseinfo_word & lowbit)) {
                              cur_dphase_delta = -cur_dphase_delta;
                            }
                          }
                          const float dosage0_incr_f = S_CAST(float, cur_dphase_delta) * S_CAST(float, kRecipDosageMax);
                          const float fxx = dosage0 + dosage0_incr_f;
                          const float fyy = dosage0 - dosage0_incr_f;
                          CopyToUnalignedOffsetF(cur_hds, &fxx, 2 * sample_idx_lowbits);
                          CopyToUnalignedOffsetF(cur_hds, &fyy, 2 * sample_idx_lowbits + 1);
                        } else {
                          CopyToUnalignedOffsetF(cur_hds, &dosage0, 2 * sample_idx_lowbits);
                          CopyToUnalignedOffsetU32(cur_hds, &eov_bits, 2 * sample_idx_lowbits + 1);
                        }
                        dosage_present_word ^= lowbit;
                      } while (dosage_present_word);
                    }
                  }
                  if (is_y) {
                    FixBcfFemaleYHdsPloidy2(pgv.genovec, pgv.dosage_ct? pgv.dosage_present : nullptr, sex_female_collapsed, sample_ct, write_iter);
                  }
                  write_iter = &(write_iter[sample_ct * sizeof(int32_t) * 2]);
                } else {
                  *write_iter++ = 0x15;
                  FillBcfHdsPloidy1(&pgv, sex_female_collapsed, hds_y_genobytes2, hds_haploid_genobytes4, sample_ct, 1, is_y, write_iter);
                  write_iter = &(write_iter[sample_ct * sizeof(int32_t)]);
                }
              }
            } else {
              // GP ignores phase.
              if ((!is_haploid) || is_x) {
                *write_iter++ = 0x35;
                FillBcfGpPloidy2(&pgv, sex_male_collapsed, sample_ct, is_x, write_iter);
                write_iter = &(write_iter[sample_ct * sizeof(int32_t) * 3]);
              } else {
                *write_iter++ = 0x25;
                FillBcfGpPloidy1(&pgv, sample_ct, write_iter);
                write_iter = &(write_iter[sample_ct * sizeof(int32_t) * 2]);
              }
            }
          }
        }
      } else {
        // multiallelic cases
        // multiallelic dosage not supported yet
        if (ref_allele_idx || (alt1_allele_idx != 1)) {
          logputs("\n");
          logerrputs("Error: BCF-export multiallelic rotation is under development.\n");
          reterr = kPglRetNotYetSupported;
          goto ExportBcf_ret_1;
        }
        if (!some_phased) {
          reterr = PgrGetM(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ExportBcf_ret_1;
          }
          if (!ds_only) {
            if (allele_ct <= 63) {
              // Still 1 byte per entry.  Almost identical to the biallelic
              // case, just need to patch in allele codes from patch_{01,10}
              // at the end.
              if ((!is_haploid) || is_x) {
                *write_iter++ = 0x21;
                FillBcf3allelicGtPloidy2(&pgv, basic_genobytes4, sample_ct, write_iter);
                if (is_x) {
                  FixBcfMaleX3allelicGtPloidy(sex_male_collapsed, sample_ct, write_iter);
                }
                write_iter = &(write_iter[sample_ct * 2]);
              } else {
                // Ploidy 1 unless there's a het-haploid.
                if (AtLeastOneMultiallelicHet(&pgv, sample_ct)) {
                  *write_iter++ = 0x21;
                  FillBcf3allelicGtPloidy2(&pgv, basic_genobytes4, sample_ct, write_iter);
                  FixBcf3allelicGtHh(sample_ct, write_iter);
                  if (is_y) {
                    FixBcfFemaleYGtPloidy2(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                  }
                  write_iter = &(write_iter[sample_ct * 2]);
                } else {
                  *write_iter++ = 0x11;
                  FillBcf3allelicGtPloidy1(&pgv, haploid_genobytes4, sample_ct, write_iter);
                  if (is_y) {
                    FixBcfFemaleYGtPloidy1(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                  }
                  write_iter = &(write_iter[sample_ct]);
                }
              }
            } else {
              // 2 bytes per entry, since allele_ct <= 255 for now.
              if ((!is_haploid) || is_x) {
                *write_iter++ = 0x22;
                FillBcf64allelicGtPloidy2(&pgv, wide_genobytes4, sample_ct, write_iter);
                if (is_x) {
                  FixBcfMaleX64allelicGtPloidy(sex_male_collapsed, sample_ct, write_iter);
                }
                write_iter = &(write_iter[sample_ct * 4]);
              } else {
                // Ploidy 1 unless there's a het-haploid.
                if (AtLeastOneMultiallelicHet(&pgv, sample_ct)) {
                  *write_iter++ = 0x22;
                  FillBcf64allelicGtPloidy2(&pgv, wide_genobytes4, sample_ct, write_iter);
                  FixBcf64allelicGtHh(sample_ct, write_iter);
                  if (is_y) {
                    FixBcfFemaleY64allelicGtPloidy2(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                  }
                  write_iter = &(write_iter[sample_ct * 4]);
                } else {
                  *write_iter++ = 0x12;
                  FillBcf64allelicGtPloidy1(&pgv, wide_haploid_genobytes4, sample_ct, write_iter);
                  if (is_y) {
                    FixBcfFemaleY64allelicGtPloidy1(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                  }
                  write_iter = &(write_iter[sample_ct * 2]);
                }
              }
            }
          }
          // yes, this will need to be reworked when proper multiallelic dosage
          // support is added
          if (ds_force) {
            ++n_fmt;
            write_iter = AppendBcfTypedInt(dosage_key_idx, write_iter);
            const uint32_t alt_allele_ct = allele_ct - 1;
            write_iter = AppendBcfVecType(5, alt_allele_ct, write_iter);
            ZeroTrailingNyps(sample_ct, pgv.genovec);
            if (!pgv.patch_01_ct) {
              ZeroWArr(sample_ctl, pgv.patch_01_set);
            }
            if (!pgv.patch_10_ct) {
              ZeroWArr(sample_ctl, pgv.patch_10_set);
            }
            FillBcfMultiallelicDsForce(&pgv, sex_male_collapsed, sex_female_collapsed, sample_ct, alt_allele_ct, is_x, is_y, is_haploid, write_iter);
            write_iter = &(write_iter[sample_ct * sizeof(int32_t) * alt_allele_ct]);
            if (hds_force) {
              ++n_fmt;
              write_iter = AppendBcfTypedInt(hds_key_idx, write_iter);
              const uint32_t max_ploidy = ((!is_haploid) || is_x)? 2 : 1;
              write_iter = AppendBcfVecType(5, max_ploidy * alt_allele_ct, write_iter);
              FillBcfMultiallelicHdsForce(&pgv, sex_male_collapsed, sex_female_collapsed, sample_ct, alt_allele_ct, is_x, is_y, is_haploid, write_iter);
              write_iter = &(write_iter[sample_ct * sizeof(int32_t) * alt_allele_ct * max_ploidy]);
            }
          }
        } else {
          // multiallelic phased
          reterr = PgrGetMP(sample_include, pssi, sample_ct, variant_uidx, simple_pgrp, &pgv);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, variant_uidx);
            goto ExportBcf_ret_1;
          }
          if (!pgv.patch_01_ct) {
            ZeroWArr(sample_ctl, pgv.patch_01_set);
          }
          if (!pgv.patch_10_ct) {
            ZeroWArr(sample_ctl, pgv.patch_10_set);
          }
          if (!pgv.phasepresent_ct) {
            ZeroWArr(sample_ctl, pgv.phasepresent);
          }
          if (allele_ct <= 63) {
            if ((!is_haploid) || is_x) {
              *write_iter++ = 0x21;
              FillBcf3allelicPhasedGt(&pgv, sample_ct, 0, prev_phased, write_iter);
              if (is_x) {
                FixBcfMaleX3allelicGtPloidy(sex_male_collapsed, sample_ct, write_iter);
              }
              write_iter = &(write_iter[sample_ct * 2]);
            } else {
              // Ploidy 1 unless there's a het-haploid.
              if (AtLeastOneMultiallelicHet(&pgv, sample_ct)) {
                *write_iter++ = 0x21;
                FillBcf3allelicPhasedGt(&pgv, sample_ct, 1, prev_phased, write_iter);
                FixBcf3allelicGtHh(sample_ct, write_iter);
                if (is_y) {
                  FixBcfFemaleYGtPloidy2(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                }
                write_iter = &(write_iter[sample_ct * 2]);
              } else {
                *write_iter++ = 0x11;
                FillBcf3allelicGtPloidy1(&pgv, haploid_genobytes4, sample_ct, write_iter);
                if (is_y) {
                  FixBcfFemaleYGtPloidy1(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                }
                write_iter = &(write_iter[sample_ct]);
              }
            }
          } else {
            // 2 bytes per entry, since allele_ct <= 255 for now.
            if ((!is_haploid) || is_x) {
              *write_iter++ = 0x22;
              FillBcf64allelicPhasedGt(&pgv, sample_ct, 0, prev_phased, write_iter);
              if (is_x) {
                FixBcfMaleX64allelicGtPloidy(sex_male_collapsed, sample_ct, write_iter);
              }
              write_iter = &(write_iter[sample_ct * 4]);
            } else {
              // Ploidy 1 unless there's a het-haploid.
              if (AtLeastOneMultiallelicHet(&pgv, sample_ct)) {
                *write_iter++ = 0x22;
                FillBcf64allelicPhasedGt(&pgv, sample_ct, 1, prev_phased, write_iter);
                FixBcf64allelicGtHh(sample_ct, write_iter);
                if (is_y) {
                  FixBcfFemaleY64allelicGtPloidy2(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                }
                write_iter = &(write_iter[sample_ct * 4]);
              } else {
                *write_iter++ = 0x12;
                FillBcf64allelicGtPloidy1(&pgv, wide_haploid_genobytes4, sample_ct, write_iter);
                if (is_y) {
                  FixBcfFemaleY64allelicGtPloidy1(pgv.genovec, sex_female_collapsed, sample_ct, write_iter);
                }
                write_iter = &(write_iter[sample_ct * 2]);
              }
            }
          }
          // yes, this will need to be reworked when proper multiallelic dosage
          // support is added
          if (ds_force) {
            ++n_fmt;
            write_iter = AppendBcfTypedInt(dosage_key_idx, write_iter);
            const uint32_t alt_allele_ct = allele_ct - 1;
            write_iter = AppendBcfVecType(5, alt_allele_ct, write_iter);
            ZeroTrailingNyps(sample_ct, pgv.genovec);
            FillBcfMultiallelicDsForce(&pgv, sex_male_collapsed, sex_female_collapsed, sample_ct, alt_allele_ct, is_x, is_y, is_haploid, write_iter);
            write_iter = &(write_iter[sample_ct * sizeof(int32_t) * alt_allele_ct]);
            if (hds_force) {
              ++n_fmt;
              write_iter = AppendBcfTypedInt(hds_key_idx, write_iter);
              const uint32_t max_ploidy = ((!is_haploid) || is_x || pgv.phasepresent_ct)? 2 : 1;
              write_iter = AppendBcfVecType(5, max_ploidy * alt_allele_ct, write_iter);
              FillBcfMultiallelicHdsForce(&pgv, sex_male_collapsed, sex_female_collapsed, sample_ct, alt_allele_ct, is_x, is_y, is_haploid, write_iter);
              write_iter = &(write_iter[sample_ct * sizeof(int32_t) * alt_allele_ct * max_ploidy]);
            }
          }
        }
      }
#ifdef __LP64__
      if (S_CAST(uintptr_t, write_iter - indiv_start) > UINT32_MAX) {
        logerrprintfww("Error: Genotype/dosage component of variant '%s' is too large to represent as BCF.\n", variant_id);
        goto ExportBcf_ret_INCONSISTENT_INPUT;
      }
#endif
      const uint32_t l_indiv = write_iter - indiv_start;
      memcpy(&(rec_start[4]), &l_indiv, sizeof(int32_t));
      rec_start[31] = n_fmt;
      if (unlikely(bgzfwrite_ck(writebuf_flush, &bgzf, &write_iter))) {
        goto ExportBcf_ret_WRITE_FAIL;
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
    if (unlikely(bgzfclose_flush(writebuf_flush, write_iter, &bgzf, &reterr))) {
      goto ExportBcf_ret_1;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    if (rlen_warning_ct) {
      fputs("done.\n", stdout);
      const uint32_t rlen_warning_print_ct = MINV(3, rlen_warning_ct);
      for (uint32_t uii = 0; uii != rlen_warning_print_ct; ++uii) {
        const char* variant_id = variant_ids[first_rlen_warning_idxs[uii]];
        MakeRlenWarningStr(variant_id);
        logerrputsb();
      }
      if (rlen_warning_ct > 3) {
        fprintf(stderr, "%u more INFO/END warning%s; see log file.\n", rlen_warning_ct - 3, (rlen_warning_ct == 4)? "" : "s");
      }
      logputs_silent("BCF export done.\n");
    } else {
      logputs("done.\n");
    }
    if (invalid_allele_code_seen) {
      logerrputs("Warning: At least one BCF allele code violates the official specification;\nother tools may not accept the file.  (Valid codes must either start with a\n'<', only contain characters in {A,C,G,T,N,a,c,g,t,n}, be an isolated '*', or\nrepresent a breakend.)\n");
    }
  }
  while (0) {
  ExportBcf_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportBcf_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pvar_info_reload, &pvar_reload_txs);
    break;
  ExportBcf_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportBcf_ret_MALFORMED_INPUT_WW:
    logputs("\n");
    WordWrapB(0);
    logerrputsb();
  ExportBcf_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ExportBcf_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 ExportBcf_ret_1:
  CleanupTextStream2(pvar_info_reload, &pvar_reload_txs, &reterr);
  CleanupBgzfCompressStream(&bgzf, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

static const Dosage kGenoToDosage[4] = {0, 16384, 32768, 65535};

typedef struct DosageTransposeCtxStruct {
  const uintptr_t* variant_include;
  const char* const* export_allele_missing;
  const AlleleCode* export_allele;
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
  const AlleleCode* export_allele = ctx->export_allele;
  const char* const* export_allele_missing = ctx->export_allele_missing;
  const uintptr_t* sample_include = ctx->sample_include;

  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
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
          if (export_allele) {
            ref_allele_idx = export_allele[variant_uidx];
          }
          uint32_t dosage_ct;
          const PglErr reterr = PgrGet1D(sample_include, pssi, sample_ct, variant_uidx, ref_allele_idx, pgrp, genovec_iter, dosage_present_iter, R_CAST(uint16_t*, dosage_main_iter), &dosage_ct);
          if (unlikely(reterr)) {
            new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
            goto DosageTransposeThread_err;
          }
          if (export_allele_missing && export_allele_missing[variant_uidx]) {
            // only keep missing vs. nonmissing distinction, zero everything
            // else out
            const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
            if (dosage_ct) {
              const Halfword* dosage_present_hw = DowncastKWToHW(dosage_present_iter);
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
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

static_assert(sizeof(Dosage) == 2, "Export012Smaj() needs to be updated.");
PglErr Export012Smaj(const char* outname, const uintptr_t* orig_sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const uintptr_t* variant_include, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* export_allele, const char* const* export_allele_missing, const char* legacy_output_missing_pheno, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t include_dom, uint32_t include_uncounted, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, char exportf_delim, PgenFileInfo* pgfip) {
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
      snprintf(g_logbuf, kLogbufSize, "Error: Too many variants for --export A%s.  (Try to work with Av instead.)\n", include_dom? "D" : "");
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
        exported_allele_idx = export_allele[variant_uidx];
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
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct, variant_ct, bigstack_left() / 4, pgr_alloc_cacheline_ct, 0, 0, 0, pgfip, &calc_thread_ct, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.read_variant_uidx_starts))) {
      goto Export012Smaj_ret_NOMEM;
    }

    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    if (unlikely(bigstack_alloc_w(raw_sample_ctl, &sample_include) ||
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
    //   read_sample_ct * stride * sizeof(Dosage) for ctx.smaj_dosagebuf
    // This is roughly
    //   read_sample_ct *
    //     (calc_thread_ct * kDosagePerCacheline * 2.375 + variant_ct * 2)
    uintptr_t bytes_avail = bigstack_left();
    // account for rounding
    const uintptr_t round_ceil = kCacheline + calc_thread_ct * (3 * kCacheline + kDosagePerCacheline * (2 * sizeof(intptr_t)));
    if (unlikely(bytes_avail < round_ceil)) {
      goto Export012Smaj_ret_NOMEM;
    }
    bytes_avail -= round_ceil;
    const uintptr_t stride = RoundUpPow2(variant_ct, kDosagePerCacheline);
    uint32_t read_sample_ct = sample_ct;
    uint32_t pass_ct = 1;
    const uintptr_t bytes_per_sample = calc_thread_ct * (kDosagePerCacheline / 8) * (3LLU + 8 * sizeof(Dosage)) + stride * sizeof(Dosage);
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
    ctx.export_allele = export_allele;
    ctx.export_allele_missing = export_allele_missing;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto Export012Smaj_ret_NOMEM;
    }
    ctx.sample_ct = read_sample_ct;
    ctx.stride = stride;
    ctx.smaj_dosagebuf = S_CAST(Dosage*, bigstack_alloc_raw_rd(read_sample_ct * S_CAST(uintptr_t, ctx.stride) * sizeof(Dosage)));
    assert(g_bigstack_base <= g_bigstack_end);
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
      uint32_t next_print_idx = (variant_ct + 99) / 100;
      for (uint32_t variant_idx = 0; ; ) {
        const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
        if (unlikely(reterr)) {
          goto Export012Smaj_ret_PGR_FAIL;
        }
        if (variant_idx) {
          JoinThreads(&tg);
          reterr = S_CAST(PglErr, ctx.err_info);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, ctx.err_info >> 32);
            goto Export012Smaj_ret_1;
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
          next_print_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
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
      next_print_idx = (read_sample_ct + 99) / 100;
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
          next_print_idx = (pct * S_CAST(uint64_t, read_sample_ct) + 99) / 100;
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
 Export012Smaj_ret_1:
  CleanupThreads(&tg);
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef USE_SHUFFLE8
// Assumes outvec_ct is positive.
void Subst4bitTo8(const void* src_vec, const VecW lookup, uint32_t outvec_ct, void* dst_vec) {
  // See Expand4bitTo8(); this just applies a substitution operation on top of
  // that.
  const VecW m4 = VCONST_W(kMask0F0F);
  const VecW* src_iter = S_CAST(const VecW*, src_vec);
  VecW* dst_iter = S_CAST(VecW*, dst_vec);
  VecW* dst_last = &(dst_iter[outvec_ct - 1]);
  while (1) {
    VecW cur_vec = *src_iter++;
    VecW vec_lo;
    VecW vec_hi;
    vecw_lo_and_hi_nybbles(cur_vec, m4, &vec_lo, &vec_hi);
    const VecW result_lo = vecw_shuffle8(lookup, vec_lo);
    const VecW result_hi = vecw_shuffle8(lookup, vec_hi);
    *dst_iter++ = result_lo;
    if (dst_iter >= dst_last) {
      if (dst_iter == dst_last) {
        *dst_iter = result_hi;
      }
      return;
    }
    *dst_iter++ = result_hi;
  }
}

// Maps 0 -> nybble0, 1 -> (nybble0 | nybble2), 2 -> nybble2, 3 -> 15.
void GenovecToNybblesSSE4(const uintptr_t* genovec, uint32_t sample_ct, uint32_t nybble0, uint32_t nybble2, uintptr_t* result_vec) {
  const uint32_t lo_u32 = (nybble0 * 0x101) | (nybble2 * 0x10100) | 0xf000000;
  const uint32_t first_u32 = lo_u32 | (nybble0 * 0x10101010);
  const uint32_t third_u32 = lo_u32 | (nybble2 * 0x10101010);
  const VecW lookup = vecw_setr32(first_u32, first_u32 | third_u32, third_u32, lo_u32 | 0xf0f0f0f0U);
  const uint32_t outvec_ct = NybbleCtToVecCt(sample_ct);
  Subst4bitTo8(genovec, lookup, outvec_ct, result_vec);
}

// writes up to vector boundary
static inline void NybblesToIupac(const uintptr_t* nybblevec, uint32_t entry_ct, uintptr_t* result_vec) {
  const VecW lookup = vecw_setr8('.', 'A', 'C', 'M', 'G', 'R', 'S', '.', 'T', 'W', 'Y', '.', 'K', '.', '.', 'N');
  const uint32_t outvec_ct = DivUp(entry_ct, kBytesPerVec);
  Subst4bitTo8(nybblevec, lookup, outvec_ct, result_vec);
}
#else
// some of this may be moved to pgenlib_misc

void GenovecToNybblesNoSSE4(const uint16_t* const* genovec_to_nybble_tables, const uintptr_t* genovec, uint32_t sample_ct, uint32_t nybble0, uint32_t nybble2, uintptr_t* result_vec) {
  const uint32_t nybble0_table_idx = ctzu32(nybble0) + 4 * (nybble0 == 15);
  const uint32_t nybble2_table_idx = ctzu32(nybble2) + 4 * (nybble2 == 15);
  const uint32_t table_idx = nybble0_table_idx + 5 * nybble2_table_idx;
  const uint16_t* cur_genovec_to_nybble_table = genovec_to_nybble_tables[table_idx];

  const unsigned char* genovec_alias = DowncastKToUc(genovec);
  uint16_t* result_alias = DowncastWToU16(result_vec);
  const uint32_t inbyte_ct = NypCtToByteCt(sample_ct);
  for (uint32_t inbyte_idx = 0; inbyte_idx != inbyte_ct; ++inbyte_idx) {
    result_alias[inbyte_idx] = cur_genovec_to_nybble_table[genovec_alias[inbyte_idx]];
  }
}

#  define PAIR_TABLE256(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p) \
  {(a), (a), (b), (a), (c), (a), (d), (a), (e), (a), (f), (a), (g), (a), (h), (a), (i), (a), (j), (a), (k), (a), (l), (a), (m), (a), (n), (a), (o), (a), (p), (a), \
  (a), (b), (b), (b), (c), (b), (d), (b), (e), (b), (f), (b), (g), (b), (h), (b), (i), (b), (j), (b), (k), (b), (l), (b), (m), (b), (n), (b), (o), (b), (p), (b), \
  (a), (c), (b), (c), (c), (c), (d), (c), (e), (c), (f), (c), (g), (c), (h), (c), (i), (c), (j), (c), (k), (c), (l), (c), (m), (c), (n), (c), (o), (c), (p), (c), \
  (a), (d), (b), (d), (c), (d), (d), (d), (e), (d), (f), (d), (g), (d), (h), (d), (i), (d), (j), (d), (k), (d), (l), (d), (m), (d), (n), (d), (o), (d), (p), (d), \
  (a), (e), (b), (e), (c), (e), (d), (e), (e), (e), (f), (e), (g), (e), (h), (e), (i), (e), (j), (e), (k), (e), (l), (e), (m), (e), (n), (e), (o), (e), (p), (e), \
  (a), (f), (b), (f), (c), (f), (d), (f), (e), (f), (f), (f), (g), (f), (h), (f), (i), (f), (j), (f), (k), (f), (l), (f), (m), (f), (n), (f), (o), (f), (p), (f), \
  (a), (g), (b), (g), (c), (g), (d), (g), (e), (g), (f), (g), (g), (g), (h), (g), (i), (g), (j), (g), (k), (g), (l), (g), (m), (g), (n), (g), (o), (g), (p), (g), \
  (a), (h), (b), (h), (c), (h), (d), (h), (e), (h), (f), (h), (g), (h), (h), (h), (i), (h), (j), (h), (k), (h), (l), (h), (m), (h), (n), (h), (o), (h), (p), (h), \
  (a), (i), (b), (i), (c), (i), (d), (i), (e), (i), (f), (i), (g), (i), (h), (i), (i), (i), (j), (i), (k), (i), (l), (i), (m), (i), (n), (i), (o), (i), (p), (i), \
  (a), (j), (b), (j), (c), (j), (d), (j), (e), (j), (f), (j), (g), (j), (h), (j), (i), (j), (j), (j), (k), (j), (l), (j), (m), (j), (n), (j), (o), (j), (p), (j), \
  (a), (k), (b), (k), (c), (k), (d), (k), (e), (k), (f), (k), (g), (k), (h), (k), (i), (k), (j), (k), (k), (k), (l), (k), (m), (k), (n), (k), (o), (k), (p), (k), \
  (a), (l), (b), (l), (c), (l), (d), (l), (e), (l), (f), (l), (g), (l), (h), (l), (i), (l), (j), (l), (k), (l), (l), (l), (m), (l), (n), (l), (o), (l), (p), (l), \
  (a), (m), (b), (m), (c), (m), (d), (m), (e), (m), (f), (m), (g), (m), (h), (m), (i), (m), (j), (m), (k), (m), (l), (m), (m), (m), (n), (m), (o), (m), (p), (m), \
  (a), (n), (b), (n), (c), (n), (d), (n), (e), (n), (f), (n), (g), (n), (h), (n), (i), (n), (j), (n), (k), (n), (l), (n), (m), (n), (n), (n), (o), (n), (p), (n), \
  (a), (o), (b), (o), (c), (o), (d), (o), (e), (o), (f), (o), (g), (o), (h), (o), (i), (o), (j), (o), (k), (o), (l), (o), (m), (o), (n), (o), (o), (o), (p), (o), \
  (a), (p), (b), (p), (c), (p), (d), (p), (e), (p), (f), (p), (g), (p), (h), (p), (i), (p), (j), (p), (k), (p), (l), (p), (m), (p), (n), (p), (o), (p), (p), (p)}

void NybblearrLookup256x1bx2(const uintptr_t* nybblearr, const void* table256x1bx2, uint32_t entry_ct, void* __restrict result) {
  const uint16_t* table_alias = S_CAST(const uint16_t*, table256x1bx2);
  const unsigned char* nybblearr_alias = DowncastKToUc(nybblearr);
  unsigned char* resultb = S_CAST(unsigned char*, result);
  const uint32_t full_byte_ct = entry_ct / 2;
  for (uint32_t byte_idx = 0; byte_idx != full_byte_ct; ++byte_idx) {
    CopyToUnalignedOffsetU16(resultb, &(table_alias[nybblearr_alias[byte_idx]]), byte_idx);
  }
  if (entry_ct % 2) {
    resultb[full_byte_ct * sizeof(int16_t)] = table_alias[nybblearr_alias[full_byte_ct]];
  }
}

void InitLookup256x05bx4(uint32_t nybble0, uint32_t nybble1, uint32_t nybble2, uint32_t nybble3, void* table256x05bx4) {
  unsigned char* table_iter = S_CAST(unsigned char*, table256x05bx4);
  const uint32_t nybble0_shifted = nybble0 << 4;
  const uint32_t nybble1_shifted = nybble1 << 4;
  const uint32_t nybble2_shifted = nybble2 << 4;
  const uint32_t nybble3_shifted = nybble3 << 4;
  unsigned char vals[16];
  vals[0] = nybble0 | nybble0_shifted;
  vals[1] = nybble1 | nybble0_shifted;
  vals[2] = nybble2 | nybble0_shifted;
  vals[3] = nybble3 | nybble0_shifted;
  vals[4] = nybble0 | nybble1_shifted;
  vals[5] = nybble1 | nybble1_shifted;
  vals[6] = nybble2 | nybble1_shifted;
  vals[7] = nybble3 | nybble1_shifted;
  vals[8] = nybble0 | nybble2_shifted;
  vals[9] = nybble1 | nybble2_shifted;
  vals[10] = nybble2 | nybble2_shifted;
  vals[11] = nybble3 | nybble2_shifted;
  vals[12] = nybble0 | nybble3_shifted;
  vals[13] = nybble1 | nybble3_shifted;
  vals[14] = nybble2 | nybble3_shifted;
  vals[15] = nybble3 | nybble3_shifted;
  for (uint32_t high_idx = 0; high_idx != 16; ++high_idx) {
    const uint32_t cur_high = vals[high_idx];
    for (uint32_t low_idx = 0; low_idx != 16; ++low_idx) {
      *table_iter++ = vals[low_idx];
      *table_iter++ = cur_high;
    }
  }
}

static const unsigned char kNybblePairToIupac[512] = PAIR_TABLE256('.', 'A', 'C', 'M', 'G', 'R', 'S', '.', 'T', 'W', 'Y', '.', 'K', '.', '.', 'N');

static inline void NybblesToIupac(const uintptr_t* nybblevec, uint32_t entry_ct, uintptr_t* result_vec) {
  NybblearrLookup256x1bx2(nybblevec, kNybblePairToIupac, entry_ct, result_vec);
}
#endif

typedef struct PhylipReadCtxStruct {
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const char* const* allele_storage;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  unsigned char allele_nybble_table[256];
  // Precomputed 256-entry lookup tables.
  // Format depends on phased vs. unphased.  No lookup table needed in unphased
  // SSE4+ case.
  const void* genovec_to_nybble_tables[25];
  uint32_t sample_ct;

  PgenReader** pgr_ptrs;
  uintptr_t** genovecs;
  uintptr_t** thread_read_mhc;
  uintptr_t** phasepresents;
  uintptr_t** phaseinfos;

  uint32_t* variant_uidx_starts;
  uint32_t cur_block_write_ct;

  // nybble-encoded IUPAC codes
  uintptr_t* vmaj_readbuf;

  uint32_t* write_sample_nm_cts;

  uint64_t err_info;
} PhylipReadCtx;

THREAD_FUNC_DECL PhylipReadThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  PhylipReadCtx* ctx = S_CAST(PhylipReadCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const char* const* allele_storage = ctx->allele_storage;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uintptr_t* sample_include = ctx->sample_include;
  const unsigned char* allele_nybble_table = ctx->allele_nybble_table;
#ifndef USE_SHUFFLE8
  const uint16_t* const* genovec_to_nybble_tables = R_CAST(const uint16_t* const*, ctx->genovec_to_nybble_tables);
#endif
  const uint32_t read_sample_ct = ctx->sample_ct;
  const uintptr_t read_sample_ctaw4 = NybbleCtToAlignedWordCt(read_sample_ct);
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  uintptr_t* genovec = ctx->genovecs[tidx];
  PgenVariant pgv;
  pgv.genovec = genovec;
  SetPgvThreadMhcNull(read_sample_ct, tidx, ctx->thread_read_mhc, &pgv);
  uint32_t* write_sample_nm_cts_iter = nullptr;
  uintptr_t prev_copy_ct = 0;
  uint64_t new_err_info = 0;
  uint32_t cur_allele_ct = 2;
  uint32_t allele_nybbles[4];
  do {
    const uintptr_t cur_block_copy_ct = ctx->cur_block_write_ct;
    const uint32_t cur_idx_end = ((tidx + 1) * cur_block_copy_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, ctx->variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
    const uint32_t cur_idx_start = (tidx * cur_block_copy_ct) / calc_thread_ct;
    uintptr_t* vmaj_readbuf_iter = &(ctx->vmaj_readbuf[(prev_copy_ct + cur_idx_start) * read_sample_ctaw4]);
    if (ctx->write_sample_nm_cts) {
      write_sample_nm_cts_iter = &(ctx->write_sample_nm_cts[prev_copy_ct + cur_idx_start]);
    }
    for (uint32_t cur_idx = cur_idx_start; cur_idx != cur_idx_end; ++cur_idx) {
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        if (unlikely(cur_allele_ct > 4)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, kPglRetInconsistentInput);
          goto PhylipReadThread_err;
        }
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      for (uint32_t aidx = 0; aidx != cur_allele_ct; ++aidx) {
        const char* cur_allele = cur_alleles[aidx];
        allele_nybbles[aidx] = allele_nybble_table[S_CAST(unsigned char, cur_allele[0])];
        if (unlikely((!allele_nybbles[aidx]) || cur_allele[1])) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, kPglRetInconsistentInput);
          goto PhylipReadThread_err;
        }
      }
      PglErr reterr;
      if (cur_allele_ct == 2) {
        reterr = PgrGet(sample_include, pssi, read_sample_ct, variant_uidx, pgrp, genovec);
      } else {
        reterr = PgrGetM(sample_include, pssi, read_sample_ct, variant_uidx, pgrp, &pgv);
      }
      if (unlikely(reterr)) {
        new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
        goto PhylipReadThread_err;
      }
      ZeroTrailingNyps(read_sample_ct, genovec);
      const uint32_t nybble0 = allele_nybbles[0];
      const uint32_t nybble1 = allele_nybbles[1];
#ifdef USE_SHUFFLE8
      GenovecToNybblesSSE4(genovec, read_sample_ct, nybble0, nybble1, vmaj_readbuf_iter);
#else
      GenovecToNybblesNoSSE4(genovec_to_nybble_tables, genovec, read_sample_ct, nybble0, nybble1, vmaj_readbuf_iter);
#endif
      uint32_t missing_allele_present = (nybble0 == 15) || (nybble1 == 15);
      if (cur_allele_ct > 2) {
        const uint32_t patch_01_ct = pgv.patch_01_ct;
        if (patch_01_ct) {
          uint32_t patch_01_nybbles[4];
          patch_01_nybbles[2] = nybble0 | allele_nybbles[2];
          missing_allele_present |= (allele_nybbles[2] == 15);
          if (cur_allele_ct == 4) {
            patch_01_nybbles[3] = nybble0 | allele_nybbles[3];
            missing_allele_present |= (allele_nybbles[3] == 15);
          }
          const uintptr_t* patch_01_set = pgv.patch_01_set;
          const AlleleCode* patch_01_vals = pgv.patch_01_vals;
          uintptr_t sample_widx = 0;
          uintptr_t patch_01_bits = patch_01_set[0];
          for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_widx, &patch_01_bits);
            const AlleleCode ac = patch_01_vals[uii];
            AssignNybblearrEntry(sample_idx, patch_01_nybbles[ac], vmaj_readbuf_iter);
          }
        }
        const uint32_t patch_10_ct = pgv.patch_10_ct;
        if (patch_10_ct) {
          const uintptr_t* patch_10_set = pgv.patch_10_set;
          const AlleleCode* patch_10_vals = pgv.patch_10_vals;
          uintptr_t sample_widx = 0;
          uintptr_t patch_10_bits = patch_10_set[0];
          for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_widx, &patch_10_bits);
            const AlleleCode ac0 = patch_10_vals[uii * 2];
            const AlleleCode ac1 = patch_10_vals[uii * 2 + 1];
            AssignNybblearrEntry(sample_idx, allele_nybbles[ac0] | allele_nybbles[ac1], vmaj_readbuf_iter);
          }
        }
      }
      if (write_sample_nm_cts_iter) {
        uint32_t missing_ct;
        if (missing_allele_present) {
          // REF=N or ALT=N is possible (see e.g. chrM POS=3107 in the original
          // 1000 Genomes phase 3 variant calls).  They should also be counted
          // as missing.
          // (This is intentionally different from vcf2phylip.py's behavior.)
          missing_ct = CountNybble(vmaj_readbuf_iter, ~k0LU, read_sample_ct);
        } else {
          missing_ct = GenoarrCountMissingUnsafe(genovec, read_sample_ct);
        }
        *write_sample_nm_cts_iter += read_sample_ct - missing_ct;
        ++write_sample_nm_cts_iter;
      }
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[read_sample_ctaw4]);
    }
    prev_copy_ct += cur_block_copy_ct;
    while (0) {
    PhylipReadThread_err:
      UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

THREAD_FUNC_DECL PhylipPhasedReadThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  PhylipReadCtx* ctx = S_CAST(PhylipReadCtx*, arg->sharedp->context);

  const uintptr_t* variant_include = ctx->variant_include;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const char* const* allele_storage = ctx->allele_storage;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uintptr_t* sample_include = ctx->sample_include;
  const unsigned char* allele_nybble_table = ctx->allele_nybble_table;
  const uint32_t* const* genovec_to_nybblepair_tables = R_CAST(const uint32_t* const*, ctx->genovec_to_nybble_tables);
  const uint32_t read_sample_ct = ctx->sample_ct;
  const uint32_t read_sample_ctl = BitCtToWordCt(read_sample_ct);
  const uintptr_t read_sample_ctaw8 = DivUp(read_sample_ct, kBytesPerVec) * kWordsPerVec;
  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  uintptr_t* genovec = ctx->genovecs[tidx];
  PgenVariant pgv;
  pgv.genovec = genovec;
  SetPgvThreadMhcNull(read_sample_ct, tidx, ctx->thread_read_mhc, &pgv);
  uintptr_t* phasepresent = nullptr;
  uintptr_t* phaseinfo = nullptr;
  if (ctx->phasepresents) {
    phasepresent = ctx->phasepresents[tidx];
    phaseinfo = ctx->phaseinfos[tidx];
  }
  pgv.phasepresent = phasepresent;
  pgv.phaseinfo = phaseinfo;
  uint32_t* write_sample_nm_cts_iter = nullptr;
  uintptr_t prev_copy_ct = 0;
  uint64_t new_err_info = 0;
  uint32_t cur_allele_ct = 2;
  uint32_t allele_nybbles[4];
  do {
    const uintptr_t cur_block_copy_ct = ctx->cur_block_write_ct;
    const uint32_t cur_idx_end = ((tidx + 1) * cur_block_copy_ct) / calc_thread_ct;
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, ctx->variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
    const uint32_t cur_idx_start = (tidx * cur_block_copy_ct) / calc_thread_ct;
    uintptr_t* vmaj_readbuf_iter = &(ctx->vmaj_readbuf[(prev_copy_ct + cur_idx_start) * read_sample_ctaw8]);
    if (ctx->write_sample_nm_cts) {
      write_sample_nm_cts_iter = &(ctx->write_sample_nm_cts[prev_copy_ct + cur_idx_start]);
    }
    for (uint32_t cur_idx = cur_idx_start; cur_idx != cur_idx_end; ++cur_idx) {
      const uintptr_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        cur_allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
        if (unlikely(cur_allele_ct > 4)) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, kPglRetInconsistentInput);
          goto PhylipPhasedReadThread_err;
        }
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      for (uint32_t aidx = 0; aidx != cur_allele_ct; ++aidx) {
        const char* cur_allele = cur_alleles[aidx];
        const uint32_t cur_nybble = allele_nybble_table[S_CAST(unsigned char, cur_allele[0])];
        if (unlikely((!cur_nybble) || cur_allele[1])) {
          new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, kPglRetInconsistentInput);
          goto PhylipPhasedReadThread_err;
        }
        allele_nybbles[aidx] = cur_nybble;
      }
      PglErr reterr;
      if (cur_allele_ct == 2) {
        reterr = PgrGetP(sample_include, pssi, read_sample_ct, variant_uidx, pgrp, genovec, phasepresent, phaseinfo, &pgv.phasepresent_ct);
      } else {
        reterr = PgrGetMP(sample_include, pssi, read_sample_ct, variant_uidx, pgrp, &pgv);
      }
      if (unlikely(reterr)) {
        new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
        goto PhylipPhasedReadThread_err;
      }
      ZeroTrailingNyps(read_sample_ct, genovec);
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      GenoarrCountFreqsUnsafe(genovec, read_sample_ct, genocounts);
      const uint32_t nybble0 = allele_nybbles[0];
      const uint32_t nybble1 = allele_nybbles[1];
      {
        const uint32_t allele0_table_idx = ctzu32(nybble0) + 4 * (nybble0 == 15);
        const uint32_t allele1_table_idx = ctzu32(nybble1) + 4 * (nybble1 == 15);
        const uint32_t table_idx = allele0_table_idx + 5 * allele1_table_idx;
        GenoarrLookup256x1bx4(genovec, genovec_to_nybblepair_tables[table_idx], read_sample_ct, vmaj_readbuf_iter);
      }
      uint32_t het_ct = genocounts[1];
      uint32_t missing_allele_present = (nybble0 == 15) || (nybble1 == 15);
      if (cur_allele_ct > 2) {
        unsigned char* nybblepairs = DowncastToUc(vmaj_readbuf_iter);
        const uint32_t patch_01_ct = pgv.patch_01_ct;
        uint32_t shifted_nybbles[4];
        shifted_nybbles[2] = allele_nybbles[2] << 4;
        missing_allele_present |= (allele_nybbles[2] == 15);
        if (cur_allele_ct == 4) {
          shifted_nybbles[3] = allele_nybbles[3] << 4;
          missing_allele_present |= (allele_nybbles[3] == 15);
        }
        if (patch_01_ct) {
          const uintptr_t* patch_01_set = pgv.patch_01_set;
          const AlleleCode* patch_01_vals = pgv.patch_01_vals;
          uintptr_t sample_widx = 0;
          uintptr_t patch_01_bits = patch_01_set[0];
          for (uint32_t uii = 0; uii != patch_01_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(patch_01_set, &sample_widx, &patch_01_bits);
            const AlleleCode ac = patch_01_vals[uii];
            nybblepairs[sample_idx] = nybble0 | shifted_nybbles[ac];
          }
        }
        const uint32_t patch_10_ct = pgv.patch_10_ct;
        if (patch_10_ct) {
          const uintptr_t* patch_10_set = pgv.patch_10_set;
          const AlleleCode* patch_10_vals = pgv.patch_10_vals;
          uintptr_t sample_widx = 0;
          uintptr_t patch_10_bits = patch_10_set[0];
          for (uint32_t uii = 0; uii != patch_10_ct; ++uii) {
            const uintptr_t sample_idx = BitIter1(patch_10_set, &sample_widx, &patch_10_bits);
            const AlleleCode ac0 = patch_10_vals[uii * 2];
            const AlleleCode ac1 = patch_10_vals[uii * 2 + 1];
            het_ct += (ac0 != ac1);
            nybblepairs[sample_idx] = allele_nybbles[ac0] | shifted_nybbles[ac1];
          }
        }
      }
      if (write_sample_nm_cts_iter) {
        uint32_t missing_ct;
        if (missing_allele_present) {
          missing_ct = CountNybble(vmaj_readbuf_iter, ~k0LU, read_sample_ct * 2);
        } else {
          missing_ct = genocounts[3] * 2;
        }
        *write_sample_nm_cts_iter += read_sample_ct * 2 - missing_ct;
        ++write_sample_nm_cts_iter;
      }
      const uint32_t phasepresent_ct = pgv.phasepresent_ct;
      if (unlikely(phasepresent_ct != het_ct)) {
        // another type of InconsistentInput
        new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, kPglRetInternalUse1);
        goto PhylipPhasedReadThread_err;
      }
      if (phasepresent_ct) {
        for (uint32_t widx = 0; widx != read_sample_ctl; ++widx) {
          const uintptr_t phasepresent_word = phasepresent[widx];
          uintptr_t phaseinfo_word = phasepresent_word & phaseinfo[widx];
          if (phaseinfo_word) {
            unsigned char* cur_nybblepairs = DowncastToUc(vmaj_readbuf_iter);
            cur_nybblepairs = &(cur_nybblepairs[widx * kBitsPerWord]);
            do {
              const uint32_t sample_idx_lowbits = ctzw(phaseinfo_word);
              const uint32_t nybblepair = cur_nybblepairs[sample_idx_lowbits];
              cur_nybblepairs[sample_idx_lowbits] = (nybblepair << 4) | (nybblepair >> 4);
              phaseinfo_word &= phaseinfo_word - 1;
            } while (phaseinfo_word);
          }
        }
      }
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[read_sample_ctaw8]);
    }
    prev_copy_ct += cur_block_copy_ct;
    while (0) {
    PhylipPhasedReadThread_err:
      UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    }
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

typedef struct PhylipTransposeCtxStruct {
  const uintptr_t* variant_include;
  uint32_t variant_ct;
  uint32_t write_sample_ct;

  const uintptr_t* vmaj_readbuf;

  VecW** thread_vecaligned_bufs;
  uintptr_t** thread_smaj_nybble_bufs;

  uint32_t sample_batch_size;

  uintptr_t* smaj_writebufs[2];
} PhylipTransposeCtx;

THREAD_FUNC_DECL PhylipTransposeWriteThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  PhylipTransposeCtx* ctx = S_CAST(PhylipTransposeCtx*, arg->sharedp->context);

  const uint32_t variant_ct = ctx->variant_ct;
  const uintptr_t variant_batch_ct = DivUp(variant_ct, kPglNybbleTransposeBatch);
  const uintptr_t write_word_stride = RoundUpPow2(variant_ct, kCacheline) / kBytesPerWord;
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uintptr_t variant_batch_idx_start = (S_CAST(uint64_t, tidx) * variant_batch_ct) / calc_thread_ct;
  VecW* vecaligned_buf = ctx->thread_vecaligned_bufs[tidx];
  uintptr_t* smaj_nybble_buf = ctx->thread_smaj_nybble_bufs[tidx];
  uintptr_t variant_batch_idx_full_end = ((S_CAST(uint64_t, tidx) + 1) * variant_batch_ct) / calc_thread_ct;
  uint32_t variant_idx_end;
  if (tidx + 1 < calc_thread_ct) {
    variant_idx_end = variant_batch_idx_full_end * kPglNybbleTransposeBatch;
  } else {
    variant_idx_end = variant_ct;
    if (variant_ct % kPglNybbleTransposeBatch) {
      --variant_batch_idx_full_end;
    }
  }
  const uint32_t sample_ct = ctx->write_sample_ct;
  const uintptr_t sample_ctaw4 = NybbleCtToAlignedWordCt(sample_ct);
  const uintptr_t* vmaj_readbuf = ctx->vmaj_readbuf;
  uint32_t sample_widx = 0;
  uint32_t parity = 0;
  do {
    const uint32_t sample_batch_size = ctx->sample_batch_size;
    const uintptr_t* vmaj_readbuf_iter = &(vmaj_readbuf[variant_batch_idx_start * kPglNybbleTransposeBatch * sample_ctaw4 + sample_widx]);
    uintptr_t* smaj_writebuf_start = &(ctx->smaj_writebufs[parity][variant_batch_idx_start * kPglNybbleTransposeWords * 2]);
    uintptr_t* smaj_writebuf_iter = smaj_writebuf_start;
    uint32_t variant_batch_size = kPglNybbleTransposeBatch;
    for (uintptr_t variant_batch_idx = variant_batch_idx_start; ; ++variant_batch_idx) {
      if (variant_batch_idx >= variant_batch_idx_full_end) {
        if (variant_batch_idx * kPglNybbleTransposeBatch >= variant_idx_end) {
          break;
        }
        variant_batch_size = variant_idx_end - variant_batch_idx * kPglNybbleTransposeBatch;
      }
      TransposeNybbleblock(vmaj_readbuf_iter, sample_ctaw4, kPglNybbleTransposeWords, variant_batch_size, sample_batch_size, smaj_nybble_buf, vecaligned_buf);
      const uintptr_t* nybble_read_iter = smaj_nybble_buf;
      uintptr_t* ascii_write_iter = smaj_writebuf_iter;
      for (uint32_t uii = 0; uii != sample_batch_size; ++uii) {
        NybblesToIupac(nybble_read_iter, variant_batch_size, ascii_write_iter);
        nybble_read_iter = &(nybble_read_iter[kPglNybbleTransposeWords]);
        ascii_write_iter = &(ascii_write_iter[write_word_stride]);
      }
      smaj_writebuf_iter = &(smaj_writebuf_iter[kPglNybbleTransposeWords * 2]);
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[variant_batch_size * sample_ctaw4]);
    }
    parity = 1 - parity;
    sample_widx += sample_batch_size / kBitsPerWordD4;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

static_assert(kTextbufSize >= kMaxMediumLine + 2 * kMaxIdSlen + 385, "ExportPhylip() used-sites writer must be updated.");
PglErr ExportPhylip(const uintptr_t* orig_sample_include, const SampleIdInfo* siip, const uintptr_t* sex_male_collapsed, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t export_phased, uint32_t export_used_sites, uint32_t max_thread_ct, IdpasteFlags exportf_id_paste, char exportf_id_delim,  uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip, char* outname, char* outname_end) {
  // Most similar to ExportPed().
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup read_tg;
  ThreadGroup transpose_tg;
  PreinitThreads(&read_tg);
  PreinitThreads(&transpose_tg);
  PhylipReadCtx read_ctx;
  PhylipTransposeCtx transpose_ctx;
  {
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    // 1. Fill genovec_to_nybble_tables, and perform initial validation if
    //    applicable.
    const unsigned char acgtn_to_nybble[5] = {1, 2, 4, 8, 15};
    if (export_phased) {
      // Verify all data is diploid.  In main loops, we also verify that all
      // data is phased.
      const uint32_t males_exist = !AllWordsAreZero(sex_male_collapsed, sample_ctl);
      if (unlikely((cip->haploid_mask[0] & 1) ||
                   (males_exist && XymtIsNonempty(variant_include, cip, kChrOffsetX)) ||
                   XymtIsNonempty(variant_include, cip, kChrOffsetY) ||
                   XymtIsNonempty(variant_include, cip, kChrOffsetMT))) {
        logerrputs("Error: --export phylip-phased can only be used with all-diploid data.\n");
        goto ExportPhylip_ret_INCONSISTENT_INPUT;
      }
      uint32_t* genovec_to_nybblepair_iter;
      if (unlikely(bigstack_alloc_u32(21 * (kLookup256x1bx4Size / sizeof(int32_t)), &genovec_to_nybblepair_iter))) {
        goto ExportPhylip_ret_NOMEM;
      }
      for (uint32_t high_acgtn = 0; high_acgtn != 5; ++high_acgtn) {
        const uint32_t high_nybble = acgtn_to_nybble[high_acgtn];
        for (uint32_t low_acgtn = 0; low_acgtn != 5; ++low_acgtn) {
          const uint32_t table_idx = high_acgtn * 5 + low_acgtn;
          if ((low_acgtn == high_acgtn) && (low_acgtn != 4)) {
            read_ctx.genovec_to_nybble_tables[table_idx] = nullptr;
            continue;
          }
          const uint32_t low_nybble = acgtn_to_nybble[low_acgtn];
          genovec_to_nybblepair_iter[0] = low_nybble * 17;
          genovec_to_nybblepair_iter[1] = low_nybble | (high_nybble << 4);
          genovec_to_nybblepair_iter[2] = high_nybble * 17;
          genovec_to_nybblepair_iter[3] = 255;
          InitLookup256x1bx4(genovec_to_nybblepair_iter);
          read_ctx.genovec_to_nybble_tables[table_idx] = genovec_to_nybblepair_iter;
          genovec_to_nybblepair_iter = &(genovec_to_nybblepair_iter[kLookup256x1bx4Size / sizeof(int32_t)]);
        }
      }
    } else {
#ifdef USE_SHUFFLE8
      // defensive
      for (uint32_t table_idx = 0; table_idx != 25; ++table_idx) {
        read_ctx.genovec_to_nybble_tables[table_idx] = nullptr;
      }
#else
      uint16_t* genovec_to_nybble_iter;
      if (unlikely(bigstack_alloc_u16(21 * 256, &genovec_to_nybble_iter))) {
        goto ExportPhylip_ret_NOMEM;
      }
      for (uint32_t low_acgtn = 0; low_acgtn != 5; ++low_acgtn) {
        const uint32_t low_nybble = acgtn_to_nybble[low_acgtn];
        for (uint32_t high_acgtn = 0; high_acgtn != 5; ++high_acgtn) {
          const uint32_t table_idx = high_acgtn * 5 + low_acgtn;
          if ((low_acgtn == high_acgtn) && (low_acgtn != 4)) {
            read_ctx.genovec_to_nybble_tables[table_idx] = nullptr;
            continue;
          }
          const uint32_t high_nybble = acgtn_to_nybble[high_acgtn];
          InitLookup256x05bx4(low_nybble, low_nybble | high_nybble, high_nybble, 15, genovec_to_nybble_iter);
          read_ctx.genovec_to_nybble_tables[table_idx] = genovec_to_nybble_iter;
          genovec_to_nybble_iter = &(genovec_to_nybble_iter[256]);
        }
      }
#endif
    }
    memset(read_ctx.allele_nybble_table, 0, 256);
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 'A')] = 1;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 'a')] = 1;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 'C')] = 2;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 'c')] = 2;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 'G')] = 4;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 'g')] = 4;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 'T')] = 8;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 't')] = 8;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 'N')] = 15;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, 'n')] = 15;
    read_ctx.allele_nybble_table[S_CAST(unsigned char, '.')] = 15;
    read_ctx.write_sample_nm_cts = nullptr;
    if (export_used_sites) {
      if (unlikely(bigstack_calloc_u32(variant_ct, &read_ctx.write_sample_nm_cts))) {
        goto ExportPhylip_ret_NOMEM;
      }
    }

    char* exported_sample_ids;
    uintptr_t max_exported_sample_id_blen;
    if (unlikely(ExportIdpaste(orig_sample_include, siip, export_phased? "phylip-phased" : "phylip", nullptr, sample_ct, exportf_id_paste, exportf_id_delim, &max_exported_sample_id_blen, &exported_sample_ids))) {
      goto ExportPhylip_ret_NOMEM;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".phy");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportPhylip_ret_OPEN_FAIL;
    }
    // yes, it's trivial to allocate this on bigstack instead
    char* textbuf = g_textbuf;
    char* textbuf_line_start;
    char* write_id_stop;
    {
      char* write_iter = u32toa_x(sample_ct << export_phased, ' ', textbuf);
      write_iter = u32toa(variant_ct, write_iter);
      if (unlikely(fwrite_checked(textbuf, write_iter - textbuf, outfile))) {
        goto ExportPhylip_ret_WRITE_FAIL;
      }
      textbuf_line_start = textbuf;
      AppendBinaryEoln(&textbuf_line_start);
      // vcf2phylip.py pads second column to start after <max ID length>+3
      // characters.  max_exported_sample_id_blen can be greater than
      // <max ID length>+1, so we need to determine true max ID length here.
      uint32_t max_id_slen = 0;
      const char* exported_sample_ids_iter = exported_sample_ids;
      for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
        const uint32_t cur_id_slen = strlen(exported_sample_ids_iter);
        if (cur_id_slen > max_id_slen) {
          max_id_slen = cur_id_slen;
        }
        exported_sample_ids_iter = &(exported_sample_ids_iter[max_exported_sample_id_blen]);
      }
      // add 2 for export_phased, since we append _A and _B to the IDs
      write_id_stop = &(textbuf_line_start[max_id_slen + 2 * export_phased + 3]);
    }

    const char id_delim = exportf_id_delim? exportf_id_delim : '_';

    // * Read phase: main thread reads raw bytes with MultireadNonempty(),
    //   while other thread(s) decode to 4-bit IUPAC
    // * Write phase: most threads transpose and then expand to ASCII, main
    //   thread flushes ASCII IUPAC codes (and adds the tiny amount of other
    //   text we need).
    const uint32_t max_allele_ct = pgfip->max_allele_ct;
    const uint32_t mhc_needed = (max_allele_ct > 2);

    // todo: check when this saturates for each phase.
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    read_ctx.phasepresents = nullptr;
    read_ctx.phaseinfos = nullptr;
    read_ctx.thread_read_mhc = nullptr;
    uint32_t read_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct, variant_ct, bigstack_left() / 4, pgr_alloc_cacheline_ct, 0, 0, 0, pgfip, &calc_thread_ct, &read_ctx.genovecs, mhc_needed? (&read_ctx.thread_read_mhc) : nullptr, export_phased? (&read_ctx.phasepresents) : nullptr, export_phased? (&read_ctx.phaseinfos) : nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &read_ctx.pgr_ptrs, &read_ctx.variant_uidx_starts))) {
      goto ExportPhylip_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &read_tg))) {
      goto ExportPhylip_ret_NOMEM;
    }
    read_ctx.variant_include = variant_include;
    read_ctx.allele_idx_offsets = allele_idx_offsets;
    read_ctx.allele_storage = allele_storage;
    read_ctx.err_info = (~0LLU) << 32;
    if (export_phased) {
      SetThreadFuncAndData(PhylipPhasedReadThread, &read_ctx, &read_tg);
    } else {
      SetThreadFuncAndData(PhylipReadThread, &read_ctx, &read_tg);
    }

    const uintptr_t variant_write_cacheline_ct = DivUp(variant_ct, kCacheline);
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    if (unlikely(SetThreadCt(calc_thread_ct, &transpose_tg) ||
                 bigstack_alloc_w(raw_sample_ctl, &sample_include) ||
                 bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_vp(calc_thread_ct, &transpose_ctx.thread_vecaligned_bufs) ||
                 bigstack_alloc_wp(calc_thread_ct, &transpose_ctx.thread_smaj_nybble_bufs))) {
      goto ExportPhylip_ret_NOMEM;
    }
    for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
      // PgenLoadMtInit() call succeeded with just 1/4 of workspace, so these
      // can't fail
      transpose_ctx.thread_vecaligned_bufs[tidx] = S_CAST(VecW*, bigstack_alloc_raw(kPglNybbleTransposeBufbytes));
      transpose_ctx.thread_smaj_nybble_bufs[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(kPglNybbleTransposeBufbytes));
    }
    // each of the two transpose-result buffers should use <= 1/6 of the
    // remaining workspace.  (this is more than for ExportPed since
    // transpose-result buffers are 4x as large, while vmaj_readbuf is only
    // twice as large.)
    const uintptr_t smaj_writebuf_cachelines_avail = bigstack_left() / (kCacheline * 6);
    uint32_t write_sample_batch_size = kPglNybbleTransposeBatch;
    if (variant_write_cacheline_ct * kPglNybbleTransposeBatch > smaj_writebuf_cachelines_avail) {
      // Until the final iteration, TransposeNybbleblock must transpose at
      // least a word's worth of nybbles per variant at a time.
      write_sample_batch_size = RoundDownPow2(smaj_writebuf_cachelines_avail / variant_write_cacheline_ct, kBitsPerWordD4);
      if (unlikely(!write_sample_batch_size)) {
        goto ExportPhylip_ret_NOMEM;
      }
    }
    transpose_ctx.smaj_writebufs[0] = S_CAST(uintptr_t*, bigstack_alloc_raw(variant_write_cacheline_ct * kCacheline * write_sample_batch_size));
    transpose_ctx.smaj_writebufs[1] = S_CAST(uintptr_t*, bigstack_alloc_raw(variant_write_cacheline_ct * kCacheline * write_sample_batch_size));

    const uintptr_t readbuf_vecs_avail = (bigstack_left() / kCacheline) * kVecsPerCacheline;
    if (unlikely(readbuf_vecs_avail < variant_ct)) {
      goto ExportPhylip_ret_NOMEM;
    }
    uintptr_t read_sample_ctv4 = readbuf_vecs_avail / variant_ct;
    uint32_t read_sample_ct;
    if (read_sample_ctv4 >= NybbleCtToVecCt(sample_ct)) {
      read_sample_ct = sample_ct;
    } else {
      read_sample_ct = read_sample_ctv4 * kNybblesPerVec;
    }
    uintptr_t read_sample_ctaw4 = NybbleCtToAlignedWordCt(read_sample_ct);
    uintptr_t* vmaj_readbuf = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(variant_ct * read_sample_ctaw4 * kBytesPerWord));
    read_ctx.vmaj_readbuf = vmaj_readbuf;
    transpose_ctx.variant_include = variant_include;
    transpose_ctx.variant_ct = variant_ct;
    transpose_ctx.vmaj_readbuf = vmaj_readbuf;
    SetThreadFuncAndData(PhylipTransposeWriteThread, &transpose_ctx, &transpose_tg);

    uint32_t sample_uidx_start = AdvTo1Bit(orig_sample_include, 0);
    const uintptr_t variant_ctaclw8 = variant_write_cacheline_ct * kWordsPerCacheline;
    const uint32_t pass_ct = 1 + (sample_ct - 1) / read_sample_ct;
    for (uint32_t pass_idx = 0; pass_idx != pass_ct; ++pass_idx) {
      memcpy(sample_include, orig_sample_include, raw_sample_ctl * sizeof(intptr_t));
      if (sample_uidx_start) {
        ClearBitsNz(0, sample_uidx_start, sample_include);
      }
      const uint32_t read_sample_idx_offset = pass_idx * read_sample_ct;
      uint32_t sample_uidx_end;
      if (pass_idx + 1 == pass_ct) {
        read_sample_ct = sample_ct - read_sample_idx_offset;
        read_sample_ctaw4 = NybbleCtToAlignedWordCt(read_sample_ct);
        sample_uidx_end = raw_sample_ct;
      } else {
        sample_uidx_end = FindNth1BitFrom(orig_sample_include, sample_uidx_start + 1, read_sample_ct);
        ClearBitsNz(sample_uidx_end, raw_sample_ct, sample_include);
      }
      FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      read_ctx.sample_include = sample_include;
      read_ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
      read_ctx.sample_ct = read_sample_ct;
      const uint32_t write_sample_ct = read_sample_ct << export_phased;
      transpose_ctx.write_sample_ct = write_sample_ct;
      if (pass_idx) {
        pgfip->block_base = main_loadbufs[0];
        // er, don't need SetBaseAndOffset0?
        PgrSetBaseAndOffset0(main_loadbufs[0], calc_thread_ct, read_ctx.pgr_ptrs);
      }
      uint32_t parity = 0;
      uint32_t read_block_idx = 0;
      ReinitThreads(&read_tg);
      uint32_t pct = 0;
      uint32_t next_print_idx = (variant_ct + 99) / 100;
      putc_unlocked('\r', stdout);
      printf("--export phylip%s pass %u/%u: loading... 0%%", export_phased? "-phased" : "", pass_idx + 1, pass_ct);
      fflush(stdout);
      for (uint32_t variant_idx = 0; ; ) {
        const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &read_tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
        if (unlikely(reterr)) {
          goto ExportPhylip_ret_PGR_FAIL;
        }
        if (variant_idx) {
          JoinThreads(&read_tg);
          reterr = S_CAST(PglErr, read_ctx.err_info);
          if (unlikely(reterr)) {
            const uint32_t variant_uidx = read_ctx.err_info >> 32;
            if (reterr == kPglRetInternalUse1) {
              logprintfww("\nError: --export phylip-phased: 0-based variant #%u is not fully phased.\n", variant_uidx);
              reterr = kPglRetInconsistentInput;
            } else if (reterr == kPglRetInconsistentInput) {
              logprintfww("\nError: --export phylip%s: 0-based variant #%u has allele code(s) outside {A,C,G,T,missing}. (Did you forget --snps-only?)\n", export_phased? "-phased" : "", variant_uidx);
            } else {
              PgenErrPrintNV(reterr, variant_uidx);
            }
            goto ExportPhylip_ret_1;
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
            goto ExportPhylip_ret_THREAD_CREATE_FAIL;
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
          next_print_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
        }
        ++read_block_idx;
        variant_idx += cur_block_write_ct;
        pgfip->block_base = main_loadbufs[parity];
      }
      // 2. Transpose and write.
      ReinitThreads(&transpose_tg);
      transpose_ctx.sample_batch_size = write_sample_batch_size;
      parity = 0;
      if (pct > 10) {
        fputs("\b \b", stdout);
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\bwriting... 0%", stdout);
      fflush(stdout);
      pct = 0;
      uint32_t flush_write_sample_idx = 0;
      next_print_idx = (write_sample_ct + 99) / 100;
      for (uint32_t flush_write_sample_idx_end = 0; ; ) {
        if (!IsLastBlock(&transpose_tg)) {
          if (flush_write_sample_idx_end + write_sample_batch_size >= write_sample_ct) {
            DeclareLastThreadBlock(&transpose_tg);
            transpose_ctx.sample_batch_size = write_sample_ct - flush_write_sample_idx_end;
          }
          if (unlikely(SpawnThreads(&transpose_tg))) {
            goto ExportPhylip_ret_THREAD_CREATE_FAIL;
          }
        }
        if (flush_write_sample_idx_end) {
          const uintptr_t* ascii_iupac_iter = transpose_ctx.smaj_writebufs[1 - parity];
          for (; flush_write_sample_idx != flush_write_sample_idx_end; ++flush_write_sample_idx) {
            const uint32_t orig_sample_idx = read_sample_idx_offset + (flush_write_sample_idx >> export_phased);
            char* write_iter = strcpya(textbuf_line_start, &(exported_sample_ids[orig_sample_idx * max_exported_sample_id_blen]));
            if (export_phased) {
              *write_iter++ = id_delim;
              *write_iter++ = 'A' + (flush_write_sample_idx % 2);
            }
            memset(write_iter, ' ', write_id_stop - write_iter);
            if (unlikely(fwrite_checked(textbuf, write_id_stop - textbuf, outfile) ||
                         fwrite_checked(ascii_iupac_iter, variant_ct, outfile))) {
              goto ExportPhylip_ret_WRITE_FAIL;
            }
            ascii_iupac_iter = &(ascii_iupac_iter[variant_ctaclw8]);
          }
          if (flush_write_sample_idx_end == write_sample_ct) {
            break;
          }
          if (flush_write_sample_idx_end >= next_print_idx) {
            if (pct > 10) {
              putc_unlocked('\b', stdout);
            }
            pct = (flush_write_sample_idx_end * 100LLU) / write_sample_ct;
            printf("\b\b%u%%", pct++);
            fflush(stdout);
            next_print_idx = (pct * S_CAST(uint64_t, write_sample_ct) + 99) / 100;
          }
        }
        JoinThreads(&transpose_tg);
        parity = 1 - parity;
        flush_write_sample_idx_end += write_sample_batch_size;
        if (flush_write_sample_idx_end > write_sample_ct) {
          flush_write_sample_idx_end = write_sample_ct;
        }
      }
      sample_uidx_start = sample_uidx_end;
      if (pct > 10) {
        fputs("\b \b", stdout);
      }
    }

#ifdef _WIN32
    if (unlikely((putc_unlocked('\r', outfile) == EOF))) {
      goto ExportPhylip_ret_WRITE_FAIL;
    }
#endif
    if (unlikely((putc_unlocked('\n', outfile) == EOF) ||
                 fclose_null(&outfile))) {
      goto ExportPhylip_ret_WRITE_FAIL;
    }
    fputs("\b\bdone.\n", stdout);
    logprintfww("--export phylip%s: %s written.\n", export_phased? "-phased" : "", outname);

    if (export_used_sites) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".used_sites.tsv");
      if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
        goto ExportPhylip_ret_OPEN_FAIL;
      }
      const uint32_t* write_sample_nm_cts = read_ctx.write_sample_nm_cts;
      char* writebuf_flush = &(textbuf[kMaxMediumLine]);
      char* write_iter = strcpya_k(textbuf, "#CHROM\tPOS\tNUM_SAMPLES" EOLN_STR);
      uint32_t chr_fo_idx = UINT32_MAX;
      uint32_t chr_end = 0;
      uint32_t chr_blen = 0;
      uint32_t min_ct = UINT32_MAX;
      char* chr_buf = &(writebuf_flush[kMaxIdSlen + 384]);
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
          *chr_name_end++ = '\t';
          chr_blen = chr_name_end - chr_buf;
        }
        write_iter = memcpya(write_iter, chr_buf, chr_blen);
        write_iter = u32toa_x(variant_bps[variant_uidx], '\t', write_iter);
        const uint32_t cur_ct = write_sample_nm_cts[variant_idx];
        if (cur_ct < min_ct) {
          min_ct = cur_ct;
        }
        write_iter = u32toa(cur_ct, write_iter);
        AppendBinaryEoln(&write_iter);
        if (unlikely(fwrite_ck(writebuf_flush, outfile, &write_iter))) {
          goto ExportPhylip_ret_WRITE_FAIL;
        }
      }
      if (unlikely(fclose_flush_null(writebuf_flush, write_iter, &outfile))) {
        goto ExportPhylip_ret_WRITE_FAIL;
      }
      logprintfww("--export phylip%s: %s written.\n", export_phased? "-phased" : "", outname);
      logprintf("Minimum non-N sample count across loci: %u/%u\n", min_ct, sample_ct << export_phased);
    }
  }
  while (0) {
  ExportPhylip_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportPhylip_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportPhylip_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ExportPhylip_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportPhylip_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  ExportPhylip_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 ExportPhylip_ret_1:
  CleanupThreads(&transpose_tg);
  CleanupThreads(&read_tg);
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr ExportEigSnp(const char* outname, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, const double* variant_cms, uint32_t variant_ct, char exportf_delim, uint32_t* hash_ptr) {
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t autosome_ct = cip->autosome_ct;
    if (unlikely(autosome_ct > 87)) {
      logerrputs("Error: EIGENSOFT .snp format does not support autosome_ct > 87.\n");
      goto ExportEigSnp_ret_INCONSISTENT_INPUT;
    }
    char chr_buf[8];
    ChrOutput output_encoding = cip->output_encoding;
    // TODO: only respect prefix setting, since ADMIXTOOLS 2 does not support
    // e.g. "X".
    if (output_encoding & (kfChrOutputM | kfChrOutputMT | kfChrOutput0M)) {
      output_encoding &= kfChrOutputPrefix;
      logprintf("Note: --output-chr setting incompatible with ADMIXTOOLS 2 interpretation of\n.snp format; treating as '%s26'.\n", (output_encoding == kfChrOutputPrefix)? "chr" : "");
    }
    logprintfww5("Writing %s ... ", outname);
    fflush(stdout);
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportEigSnp_ret_OPEN_FAIL;
    }
    // multicharacter allele codes prohibited, so lines can't be long
    char* write_iter = g_textbuf;
    char* textbuf_flush = &(write_iter[kMaxMediumLine]);
    uint32_t ref_allele_idx = 0;
    uintptr_t variant_uidx_base = 0;
    uintptr_t cur_bits = variant_include[0];
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t eighash = 0;
    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        if (chr_idx > cip->max_numeric_code) {
          if (chr_idx <= cip->max_code) {
            // PAR1/PAR2 -> XY
            chr_idx = autosome_ct + (kChrOffsetXY + 1);
          } else {
            if (unlikely(!cip->zero_extra_chrs)) {
              snprintf(g_logbuf, kLogbufSize, "Error: %s cannot contain nonstandard chromosome codes.\n", outname);
              goto ExportEigSnp_ret_INCONSISTENT_INPUT_WW_N;
            }
            chr_idx = 0;
          }
        }
        char* chr_name_end = chr_buf;
        if (!chr_idx) {
          *chr_name_end++ = '0';
        } else if ((chr_idx < autosome_ct + (kChrOffsetXY + 1)) || (chr_idx > autosome_ct + (kChrOffsetMT + 1))) {
          chr_name_end = ChrNameStdEx(cip, chr_idx, output_encoding, chr_buf);
        } else {
          // MT -> 90, XY -> 91
          chr_name_end = chr_buf;
          if (output_encoding & kfChrOutputPrefix) {
            chr_name_end = strcpya_k(chr_name_end, "chr");
          }
          *chr_name_end++ = '9';
          *chr_name_end++ = '0' + (chr_idx == autosome_ct + (kChrOffsetXY + 1));
        }
        *chr_name_end = exportf_delim;
        chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
      }
      // ID, CHROM, CM, POS, REF, ALT
      // convertf documentation states that sample IDs are limited to 39
      // characters, but says nothing about variant IDs even though it's unable
      // to handle longer IDs.  So we don't enforce a variant ID length limit
      // for now.
      const char* variant_id = variant_ids[variant_uidx];
      const uint32_t variant_id_slen = strlen(variant_id);
      UpdateEighash(variant_id, variant_id_slen, &eighash);
      write_iter = memcpyax(write_iter, variant_id, variant_id_slen, exportf_delim);
      write_iter = memcpya(write_iter, chr_buf, chr_buf_blen);
      if (variant_cms) {
        write_iter = dtoa_g_p8(variant_cms[variant_uidx], write_iter);
      } else {
        *write_iter++ = '0';
      }
      *write_iter++ = exportf_delim;
      write_iter = u32toa_x(variant_bps[variant_uidx], exportf_delim, write_iter);

      uintptr_t allele_idx_offset_base = variant_uidx * 2;
      if (allele_idx_offsets) {
        allele_idx_offset_base = allele_idx_offsets[variant_uidx];
        if (unlikely(allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base != 2)) {
          snprintf(g_logbuf, kLogbufSize, "Error: %s cannot contain multiallelic variants.\n", outname);
          goto ExportEigSnp_ret_INCONSISTENT_INPUT_WW_N;
        }
      }
      const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
      if (allele_permute) {
        ref_allele_idx = allele_permute[allele_idx_offset_base];
      }
      const char* ref_allele = cur_alleles[ref_allele_idx];
      const char* alt_allele = cur_alleles[1 - ref_allele_idx];
      if (unlikely(ref_allele[1] || alt_allele[1])) {
        snprintf(g_logbuf, kLogbufSize, "Error: %s cannot contain multi-character allele codes.\n", outname);
        goto ExportEigSnp_ret_INCONSISTENT_INPUT_WW_N;
      }
      *write_iter++ = ref_allele[0];
      *write_iter++ = exportf_delim;
      char alt_char = alt_allele[0];
      if (alt_char == '.') {
        alt_char = 'X';
      }
      *write_iter++ = alt_char;
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto ExportEigSnp_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto ExportEigSnp_ret_WRITE_FAIL;
    }
    logputs("done.\n");
    *hash_ptr = eighash;
  }
  while (0) {
  ExportEigSnp_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportEigSnp_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportEigSnp_ret_INCONSISTENT_INPUT_WW_N:
    logputs("\n");
    WordWrapB(0);
    logerrputsb();
  ExportEigSnp_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  fclose_cond(outfile);
  return reterr;
}

PglErr ExportEigInd(const char* outname, const uintptr_t* sample_include, const SampleIdInfo* siip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, uint32_t sample_ct, uint32_t pheno_ct, char exportf_delim, IdpasteFlags exportf_id_paste, char exportf_id_delim, uint32_t* hash_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportEigInd_ret_OPEN_FAIL;
    }
    const uintptr_t* pheno_nm = nullptr;
    const uintptr_t* pheno_cc = nullptr;
    const double* pheno_qt = nullptr;
    const uint32_t* pheno_cat = nullptr;
    const char** cat_names = nullptr;
    if (pheno_ct) {
      const PhenoDtype type_code = pheno_cols[0].type_code;
      pheno_nm = pheno_cols[0].nonmiss;
      if (type_code == kPhenoDtypeCc) {
        pheno_cc = pheno_cols[0].data.cc;
      } else if (type_code == kPhenoDtypeQt) {
        pheno_qt = pheno_cols[0].data.qt;
      } else {
        pheno_cat = pheno_cols[0].data.cat;
        cat_names = pheno_cols[0].category_names;
      }
    }
    char* exported_sample_ids;
    uintptr_t max_exported_sample_id_blen;
    if (unlikely(ExportIdpaste(sample_include, siip, ".ind", "eigfile", sample_ct, exportf_id_paste, exportf_id_delim, &max_exported_sample_id_blen, &exported_sample_ids))) {
      goto ExportEigInd_ret_NOMEM;
    }
    if (unlikely(max_exported_sample_id_blen > 40)) {
      logerrputs("Error: Sample IDs in .ind files are limited to 39 characters.\n");
      goto ExportEigInd_ret_INCONSISTENT_INPUT;
    }
    logprintfww5("Writing %s ... ", outname);
    fflush(stdout);
    char* write_iter = g_textbuf;
    char* textbuf_flush = &(write_iter[kMaxMediumLine]);
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    uint32_t eighash = 0;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uint32_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      // ID, SEX, PHENO
      const char* exported_sample_id = &(exported_sample_ids[sample_idx * max_exported_sample_id_blen]);
      const uint32_t exported_sample_id_slen = strlen(exported_sample_id);
      UpdateEighash(exported_sample_id, exported_sample_id_slen, &eighash);
      write_iter = memcpyax(write_iter, exported_sample_id, exported_sample_id_slen, exportf_delim);
      char sex_char = 'U';
      if (IsSet(sex_nm, sample_uidx)) {
        sex_char = IsSet(sex_male, sample_uidx)? 'M' : 'F';
      }
      *write_iter++ = sex_char;
      *write_iter++ = exportf_delim;
      if (pheno_nm && IsSet(pheno_nm, sample_uidx)) {
        if (pheno_cc) {
          if (IsSet(pheno_cc, sample_uidx)) {
            write_iter = strcpya_k(write_iter, "Case");
          } else {
            write_iter = strcpya_k(write_iter, "Control");
          }
        } else if (pheno_qt) {
          write_iter = dtoa_g(pheno_qt[sample_uidx], write_iter);
        } else {
          write_iter = strcpya(write_iter, cat_names[pheno_cat[sample_uidx]]);
        }
      } else {
        // May want to make this configurable.
        write_iter = strcpya_k(write_iter, "Ignore");
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_ck(textbuf_flush, outfile, &write_iter))) {
        goto ExportEigInd_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_flush_null(textbuf_flush, write_iter, &outfile))) {
      goto ExportEigInd_ret_WRITE_FAIL;
    }
    logputs("done.\n");
    *hash_ptr = eighash;
  }
  while (0) {
  ExportEigInd_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportEigInd_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportEigInd_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportEigInd_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct ExportEigGenoCtxStruct {
  const uintptr_t* variant_include;
  const uintptr_t* allele_idx_offsets;
  const uintptr_t* sample_include;
  const uint32_t* sample_include_cumulative_popcounts;
  const AlleleCode* allele_permute;
  uint32_t sample_ct;

  PgenReader** pgr_ptrs;
  uintptr_t** refcountvecs;
  uint32_t* variant_uidx_starts;

  uint32_t cur_block_write_ct;

  unsigned char* writebufs[2];

  // high 32 bits = variant_uidx, earlier one takes precedence
  // low 32 bits = uint32_t(PglErr)
  uint64_t err_info;
} ExportEigGenoCtx;

static_assert(kBytesPerVec <= 48, "ExportEigGenoThread() needs to be updated.");
THREAD_FUNC_DECL ExportEigGenoThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  ExportEigGenoCtx* ctx = S_CAST(ExportEigGenoCtx*, arg->sharedp->context);

  PgenReader* pgrp = ctx->pgr_ptrs[tidx];
  uintptr_t* refcountvec = ctx->refcountvecs[tidx];
  const uint32_t sample_ct = ctx->sample_ct;
  const uintptr_t* variant_include = ctx->variant_include;
  const uintptr_t* allele_idx_offsets = ctx->allele_idx_offsets;
  const uintptr_t* sample_include = ctx->sample_include;
  PgrSampleSubsetIndex pssi;
  PgrSetSampleSubsetIndex(ctx->sample_include_cumulative_popcounts, pgrp, &pssi);
  const uint32_t calc_thread_ct = GetThreadCt(arg->sharedp);
  const uint32_t sample_ct4 = DivUp(sample_ct, 4);
  const uintptr_t output_rec_blen = MAXV(48, sample_ct4);
#ifdef USE_SSE2
  const uint32_t vec_ct_m1 = (sample_ct - 1) / kNypsPerVec;
  const uint32_t final_byte_offset = (sample_ct4 < kBytesPerVec)? 0 : (sample_ct4 - kBytesPerVec);
  const VecW mask_0f0f = vecw_set1(kMask0F0F);
  // tried EigGenoToPgenThread()'s lookup-table approach, no noticeable
  // efficiency improvement.
  const VecW mask_3333 = vecw_set1(kMask3333);
#else
  const uint32_t word_ct_m1 = (sample_ct - 1) / kBitsPerWordD2;
  const uint32_t final_byte_offset = (sample_ct4 < kBytesPerWord)? 0 : (sample_ct4 - kBytesPerWord);
#endif
  const AlleleCode* allele_permute = ctx->allele_permute;
  uint32_t ref_allele_idx = 0;
  uint32_t parity = 0;
  uint64_t new_err_info = 0;
  do {
    const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
    uint32_t write_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t write_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* writebuf_iter = &(ctx->writebufs[parity][write_idx * output_rec_blen]);
    uintptr_t variant_uidx_base;
    uintptr_t cur_bits;
    BitIter1Start(variant_include, ctx->variant_uidx_starts[tidx], &variant_uidx_base, &cur_bits);
    for (; write_idx != write_idx_end; ++write_idx) {
      const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
      if (allele_permute) {
        uintptr_t allele_idx_offset_base = variant_uidx * 2;
        if (allele_idx_offsets) {
          allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          // variant should already be validated to be biallelic by
          // ExportEigSnp()
        }
        ref_allele_idx = allele_permute[allele_idx_offset_base];
      }
      PglErr reterr = PgrGet1(sample_include, pssi, sample_ct, variant_uidx, ref_allele_idx, pgrp, refcountvec);
      if (unlikely(reterr)) {
        new_err_info = (S_CAST(uint64_t, variant_uidx) << 32) | S_CAST(uint32_t, reterr);
        goto ExportEigGenoThread_err;
      }
      // We already got REF allele counts, so only need to reverse nyp order.
      ZeroTrailingNyps(sample_ct, refcountvec);
      unsigned char* refcountvec_uc = DowncastToUc(refcountvec);
      unsigned char* write_uc = writebuf_iter;
      uint32_t byte_offset = 0;
#ifdef USE_SSE2
      for (uint32_t vidx = 0; ; ++vidx) {
        if (vidx >= vec_ct_m1) {
          if (vidx > vec_ct_m1) {
            break;
          }
          byte_offset = final_byte_offset;
        }
        const VecW input_vec = vecw_loadu(&(refcountvec_uc[byte_offset]));
        const VecW high_nybbles_shifted = vecw_srli(input_vec, 4) & mask_0f0f;
        const VecW low_nybbles = input_vec & mask_0f0f;
        const VecW nybble_swapped = high_nybbles_shifted | vecw_slli(low_nybbles, 4);
        const VecW high_nyps = vecw_and_notfirst(mask_3333, nybble_swapped);
        const VecW low_nyps = mask_3333 & nybble_swapped;
        const VecW result = vecw_srli(high_nyps, 2) | vecw_slli(low_nyps, 2);
        vecw_storeu(&(write_uc[byte_offset]), result);
        byte_offset += kBytesPerVec;
      }
#else
      for (uint32_t widx = 0; ; ++widx) {
        if (widx >= word_ct_m1) {
          if (widx > word_ct_m1) {
            break;
          }
          byte_offset = final_byte_offset;
        }
        uintptr_t input_word;
        CopyFromUnalignedW(&input_word, &(refcountvec_uc[byte_offset]));
        const uintptr_t high_nybbles = input_word & (~kMask0F0F);
        const uintptr_t low_nybbles = input_word & kMask0F0F;
        const uintptr_t nybble_swapped = (high_nybbles >> 4) | (low_nybbles << 4);
        const uintptr_t high_nyps = nybble_swapped & (~kMask3333);
        const uintptr_t low_nyps = nybble_swapped & kMask3333;
        const uintptr_t result = (high_nyps >> 2) | (low_nyps << 2);
        CopyToUnalignedW(&(write_uc[byte_offset]), &result);
        byte_offset += kBytesPerWord;
      }
#endif
      if (sample_ct4 < 48) {
        // in non-vectorized case, could skip this after first two iterations
        memset(&(writebuf_iter[sample_ct4]), 0, 48 - sample_ct4);
      }
      writebuf_iter = &(writebuf_iter[output_rec_blen]);
    }
    parity = 1 - parity;
  } while (!THREAD_BLOCK_FINISH(arg));
  while (0) {
  ExportEigGenoThread_err:
    UpdateU64IfSmaller(new_err_info, &ctx->err_info);
    THREAD_BLOCK_FINISH(arg);
    break;
  }
  THREAD_RETURN;
}

PglErr ExportEigGeno(const char* outname, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const AlleleCode* allele_permute, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t s_hash, uint32_t v_hash, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  ExportEigGenoCtx ctx;
  {
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportEigGeno_ret_OPEN_FAIL;
    }
    const uintptr_t output_rec_blen = MAXV(48, DivUp(sample_ct, 4));
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // todo: wouldn't be surprised if this saturated at 1 thread without
    // sample-subsetting, and less than 4 even with subsetting...
    if (calc_thread_ct > 4) {
      calc_thread_ct = 4;
    }

    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    uint32_t read_block_size;
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct, variant_ct, bigstack_left(), pgr_alloc_cacheline_ct, 0, output_rec_blen * 2, 0, pgfip, &calc_thread_ct, &ctx.refcountvecs, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &ctx.pgr_ptrs, &ctx.variant_uidx_starts))) {
      goto ExportEigGeno_ret_NOMEM;
    }
    if (unlikely(bigstack_alloc_uc(output_rec_blen * read_block_size, &(ctx.writebufs[0])) ||
                 bigstack_alloc_uc(output_rec_blen * read_block_size, &(ctx.writebufs[1])))) {
      // this shouldn't be possible
      assert(0);
      goto ExportEigGeno_ret_NOMEM;
    }
    // write header
    memset(ctx.writebufs[0], 0, output_rec_blen);
    snprintf(R_CAST(char*, ctx.writebufs[0]), output_rec_blen, "GENO %7d %7d %x %x", sample_ct, variant_ct, s_hash, v_hash);
    if (unlikely(fwrite_checked(ctx.writebufs[0], output_rec_blen, outfile))) {
      goto ExportEigGeno_ret_WRITE_FAIL;
    }

    ctx.sample_ct = sample_ct;
    ctx.variant_include = variant_include;
    ctx.allele_idx_offsets = allele_idx_offsets;
    ctx.sample_include = sample_include;
    ctx.sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
    if (unlikely(SetThreadCt(calc_thread_ct, &tg))) {
      goto ExportEigGeno_ret_NOMEM;
    }
    ctx.allele_permute = allele_permute;
    ctx.err_info = (~0LLU) << 32;
    SetThreadFuncAndData(ExportEigGenoThread, &ctx, &tg);

    // Main workflow:
    // 1. Set n=0, load/skip block 0
    //
    // 2. Spawn threads processing block n
    // 3. If n>0, write results for block (n-1) (this could be done by another
    //    thread?)
    // 4. Increment n by 1
    // 5. Load/skip block n unless eof
    // 6. Join threads
    // 7. Goto step 2 unless eof
    //
    // 8. Write results for last block
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;

    uint32_t prev_block_write_ct = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = (variant_ct + 99) / 100;
    logprintfww5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    for (uint32_t variant_idx = 0; ; ) {
      const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
      if (unlikely(reterr)) {
        goto ExportEigGeno_ret_PGR_FAIL;
      }
      if (variant_idx) {
        JoinThreads(&tg);
        reterr = S_CAST(PglErr, ctx.err_info);
        if (unlikely(reterr)) {
          PgenErrPrintNV(reterr, ctx.err_info >> 32);
          goto ExportEigGeno_ret_1;
        }
      }
      if (!IsLastBlock(&tg)) {
        ctx.cur_block_write_ct = cur_block_write_ct;
        ComputeUidxStartPartition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, ctx.variant_uidx_starts);
        PgrCopyBaseAndOffset(pgfip, calc_thread_ct, ctx.pgr_ptrs);
        if (variant_idx + cur_block_write_ct == variant_ct) {
          DeclareLastThreadBlock(&tg);
        }
        if (unlikely(SpawnThreads(&tg))) {
          goto ExportEigGeno_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const unsigned char* outdata = ctx.writebufs[parity];
        if (unlikely(fwrite_checked(outdata, output_rec_blen * prev_block_write_ct, outfile))) {
          goto ExportEigGeno_ret_WRITE_FAIL;
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
      ++read_block_idx;
      prev_block_write_ct = cur_block_write_ct;
      variant_idx += cur_block_write_ct;
      pgfip->block_base = main_loadbufs[parity];
    }
    if (unlikely(fclose_null(&outfile))) {
      goto ExportEigGeno_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    logputs("done.\n");
  }
  while (0) {
  ExportEigGeno_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportEigGeno_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportEigGeno_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ExportEigGeno_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportEigGeno_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 ExportEigGeno_ret_1:
  CleanupThreads(&tg);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

typedef struct TransposeToEigTgenoWriteCtxStruct {
  uint32_t variant_ct;
  uint32_t sample_ct;

  const uintptr_t* vmaj_readbuf;

  VecW** thread_vecaligned_bufs;

  uint32_t sample_batch_size;

  uintptr_t* smaj_writebufs[2];
} TransposeToEigTgenoWriteCtx;

THREAD_FUNC_DECL TransposeToEigTgenoWriteThread(void* raw_arg) {
  // nearly identical to TransposeToPlink1SmajWriteThread().
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uintptr_t tidx = arg->tidx;
  TransposeToEigTgenoWriteCtx* ctx = S_CAST(TransposeToEigTgenoWriteCtx*, arg->sharedp->context);

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
#ifdef USE_SSE2
  const uint32_t last_word_idx = (thread_variant_ct - 1) / kBitsPerWordD2;
#endif
  uintptr_t last_word_mask = ~k0LU;
  uint32_t trailing_zero_word_stop = 0;
  if (variant_idx_end == variant_ct) {
    // mask out trailing bits... except that we mask out *low* bits in the last
    // mixed byte
    const uint32_t mask_out_nyp_ct = (-thread_variant_ct) & (kBitsPerWordD2 - 1);
    const uint32_t mask_out_nyp_remainder_ct = mask_out_nyp_ct % 4;
    last_word_mask = ((k1LU - (k1LU << (2 * mask_out_nyp_remainder_ct))) << (kBitsPerWord - 8)) - 1;
    const uint32_t mask_out_fullbyte_ct = mask_out_nyp_ct / 4;
    last_word_mask >>= mask_out_fullbyte_ct * 8;
    if (variant_ct <= 192 - 4 * kBytesPerWord) {
      trailing_zero_word_stop = 48 / kBytesPerWord;
    }
  }
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
#ifdef USE_SSE2
    const uint32_t vec_ct = NypCtToVecCt(thread_variant_ct);
    const VecW mask_0f0f = vecw_set1(kMask0F0F);
#  ifdef USE_SHUFFLE8
    const VecW eig_lookup_high = vecw_setr8(10, 6, 2, 14, 9, 5, 1, 13,
                                            8, 4, 0, 12, 11, 7, 3, 15);
    const VecW eig_lookup_low = vecw_slli(eig_lookup_high, 4);
#  else
    const VecW mask_3333 = vecw_set1(kMask3333);
    const VecW mask_5555 = vecw_set1(kMask5555);
#  endif
#else
    const uint32_t word_ct = NypCtToWordCt(thread_variant_ct);
#endif
    for (uint32_t sample_idx = 0; sample_idx != sample_batch_size; ++sample_idx) {
#ifdef USE_SSE2
      VecW* __attribute__((may_alias)) smaj_writebuf_viter = R_CAST(VecW*, smaj_writebuf_iter);
      for (uint32_t vidx = 0; vidx != vec_ct; ++vidx) {
        const VecW input_vec = *smaj_writebuf_viter;
        const VecW high_nybbles_shifted = vecw_srli(input_vec, 4) & mask_0f0f;
        const VecW low_nybbles = input_vec & mask_0f0f;
#  ifdef USE_SHUFFLE8
        const VecW low_transformed = vecw_shuffle8(eig_lookup_low, low_nybbles);
        const VecW high_transformed = vecw_shuffle8(eig_lookup_high, high_nybbles_shifted);
        const VecW result = low_transformed | high_transformed;
#  else
        const VecW nybble_swapped = high_nybbles_shifted | vecw_slli(low_nybbles, 4);
        const VecW high_nyps = vecw_and_notfirst(mask_3333, nybble_swapped);
        const VecW low_nyps = mask_3333 & nybble_swapped;
        const VecW nyp_reversed = vecw_srli(high_nyps, 2) | vecw_slli(low_nyps, 2);
        const VecW is_even = vecw_and_notfirst(nyp_reversed, mask_5555);
        const VecW result = nyp_reversed ^ vecw_slli(is_even, 1);
#  endif
        *smaj_writebuf_viter++ = result;
      }
      smaj_writebuf_iter[last_word_idx] &= last_word_mask;
      if (trailing_zero_word_stop) {
        for (uint32_t widx = last_word_idx + 1; widx != trailing_zero_word_stop; ++widx) {
          smaj_writebuf_iter[widx] = 0;
        }
      }
#else
      for (uint32_t widx = 0; widx != word_ct; ++widx) {
        uintptr_t input_word = smaj_writebuf_iter[widx];
        const uintptr_t high_nybbles = input_word & (~kMask0F0F);
        const uintptr_t low_nybbles = input_word & kMask0F0F;
        const uintptr_t nybble_swapped = (high_nybbles >> 4) | (low_nybbles << 4);
        const uintptr_t high_nyps = nybble_swapped & (~kMask3333);
        const uintptr_t low_nyps = nybble_swapped & kMask3333;
        const uintptr_t nyp_reversed = (high_nyps >> 2) | (low_nyps << 2);
        const uintptr_t is_even = (~nyp_reversed) & kMask5555;
        smaj_writebuf_iter[widx] = nyp_reversed ^ (is_even << 1);
      }
      smaj_writebuf_iter[word_ct - 1] &= last_word_mask;
      if (trailing_zero_word_stop) {
        for (uint32_t widx = word_ct; widx != trailing_zero_word_stop; ++widx) {
          smaj_writebuf_iter[widx] = 0;
        }
      }
#endif
      // todo: 48-byte minimum
      smaj_writebuf_iter = &(smaj_writebuf_iter[variant_batch_word_ct]);
    }
    parity = 1 - parity;
    sample_widx += sample_batch_size / kBitsPerWordD2;
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

static_assert(kCacheline >= 48, "ExportEigTgeno() needs to be updated.");
PglErr ExportEigTgeno(const char* outname, const uintptr_t* orig_sample_include, const uintptr_t* variant_include, const uintptr_t* allele_idx_offsets, const AlleleCode* allele_permute, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t s_hash, uint32_t v_hash, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, PgenFileInfo* pgfip) {
  // Nearly identical to ExportIndMajorBed().
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup read_tg;
  ThreadGroup write_tg;
  PreinitThreads(&read_tg);
  PreinitThreads(&write_tg);
  MTPgenReadCtx read_ctx;
  TransposeToEigTgenoWriteCtx write_ctx;
  {
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto ExportEigTgeno_ret_OPEN_FAIL;
    }
    {
      char header[48];
      memset(header, 0, 48);
      snprintf(header, 48, "TGENO %7d %7d %x %x", sample_ct, variant_ct, s_hash, v_hash);
      if (unlikely(!fwrite_unlocked(header, 48, 1, outfile))) {
        goto ExportEigTgeno_ret_WRITE_FAIL;
      }
    }
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // todo: if only 1 pass is needed, and no subsetting is happening, read
    // saturates at ~4 threads?
    STD_ARRAY_DECL(unsigned char*, 2, main_loadbufs);
    uint32_t read_block_size;
    // note that this is restricted to half of available workspace
    if (unlikely(PgenMtLoadInit(variant_include, sample_ct, variant_ct, bigstack_left() / 2, pgr_alloc_cacheline_ct, 0, 0, 0, pgfip, &calc_thread_ct, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, &read_block_size, nullptr, main_loadbufs, &read_ctx.pgr_ptrs, &read_ctx.variant_uidx_starts))) {
      goto ExportEigTgeno_ret_NOMEM;
    }
    if (unlikely(SetThreadCt(calc_thread_ct, &read_tg))) {
      goto ExportEigTgeno_ret_NOMEM;
    }
    read_ctx.variant_include = variant_include;
    read_ctx.allele_idx_offsets = allele_idx_offsets;
    read_ctx.allele_permute = allele_permute;
    read_ctx.err_info = (~0LLU) << 32;
    SetThreadFuncAndData(MTPgenReadThread, &read_ctx, &read_tg);

    const uintptr_t variant_cacheline_ct = NypCtToCachelineCt(variant_ct);
    uint32_t output_calc_thread_ct = MINV(calc_thread_ct, variant_cacheline_ct);
    // todo: recheck this threshold
    if (output_calc_thread_ct > 4) {
      output_calc_thread_ct = 4;
    }
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    uintptr_t* sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    if (unlikely(SetThreadCt(output_calc_thread_ct, &write_tg) ||
                 bigstack_alloc_w(raw_sample_ctl, &sample_include) ||
                 bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_vp(output_calc_thread_ct, &write_ctx.thread_vecaligned_bufs))) {
      goto ExportEigTgeno_ret_NOMEM;
    }
    for (uint32_t tidx = 0; tidx != output_calc_thread_ct; ++tidx) {
      write_ctx.thread_vecaligned_bufs[tidx] = S_CAST(VecW*, bigstack_alloc_raw(kPglNypTransposeBufbytes));
    }
    if (unlikely(g_bigstack_base > g_bigstack_end)) {
      goto ExportEigTgeno_ret_NOMEM;
    }
    // each of the two write buffers should use <= 1/8 of the remaining
    // workspace
    const uintptr_t writebuf_cachelines_avail = bigstack_left() / (kCacheline * 8);
    uint32_t sample_batch_size = kPglNypTransposeBatch;
    if (variant_cacheline_ct * kPglNypTransposeBatch > writebuf_cachelines_avail) {
      sample_batch_size = RoundDownPow2(writebuf_cachelines_avail / variant_cacheline_ct, kBitsPerWordD2);
      if (unlikely(!sample_batch_size)) {
        goto ExportEigTgeno_ret_NOMEM;
      }
    }
    write_ctx.smaj_writebufs[0] = S_CAST(uintptr_t*, bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size));
    write_ctx.smaj_writebufs[1] = S_CAST(uintptr_t*, bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size));
    const uintptr_t readbuf_vecs_avail = (bigstack_left() / kCacheline) * kVecsPerCacheline;
    if (unlikely(readbuf_vecs_avail < variant_ct)) {
      goto ExportEigTgeno_ret_NOMEM;
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
    write_ctx.variant_ct = variant_ct;
    write_ctx.vmaj_readbuf = vmaj_readbuf;
    SetThreadFuncAndData(TransposeToEigTgenoWriteThread, &write_ctx, &write_tg);
    uint32_t sample_uidx_start = AdvTo1Bit(orig_sample_include, 0);
    const uintptr_t output_rec_blen = MAXV(48, NypCtToByteCt(variant_ct));
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
      uint32_t next_print_idx = (variant_ct + 99) / 100;
      putc_unlocked('\r', stdout);
      printf("--export eigt pass %u/%u: loading... 0%%", pass_idx + 1, pass_ct);
      fflush(stdout);
      for (uint32_t variant_idx = 0; ; ) {
        const uint32_t cur_block_write_ct = MultireadNonempty(variant_include, &read_tg, raw_variant_ct, read_block_size, pgfip, &read_block_idx, &reterr);
        if (unlikely(reterr)) {
          goto ExportEigTgeno_ret_PGR_FAIL;
        }
        if (variant_idx) {
          JoinThreads(&read_tg);
          reterr = S_CAST(PglErr, read_ctx.err_info);
          if (unlikely(reterr)) {
            PgenErrPrintNV(reterr, read_ctx.err_info >> 32);
            goto ExportEigTgeno_ret_1;
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
            goto ExportEigTgeno_ret_THREAD_CREATE_FAIL;
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
          next_print_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
        }

        ++read_block_idx;
        variant_idx += cur_block_write_ct;
        pgfip->block_base = main_loadbufs[parity];
      }
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
      next_print_idx = (read_sample_ct + 99) / 100;
      for (uint32_t flush_sample_idx_end = 0; ; ) {
        if (!IsLastBlock(&write_tg)) {
          if (flush_sample_idx_end + sample_batch_size >= read_sample_ct) {
            DeclareLastThreadBlock(&write_tg);
            write_ctx.sample_batch_size = read_sample_ct - flush_sample_idx_end;
          }
          if (unlikely(SpawnThreads(&write_tg))) {
            goto ExportEigTgeno_ret_THREAD_CREATE_FAIL;
          }
        }
        if (flush_sample_idx_end) {
          uintptr_t* smaj_writebuf_iter = write_ctx.smaj_writebufs[1 - parity];
          for (; flush_sample_idx != flush_sample_idx_end; ++flush_sample_idx) {
            fwrite_unlocked(smaj_writebuf_iter, output_rec_blen, 1, outfile);
            smaj_writebuf_iter = &(smaj_writebuf_iter[variant_ctaclw2]);
          }
          if (unlikely(ferror_unlocked(outfile))) {
            goto ExportEigTgeno_ret_WRITE_FAIL;
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
            next_print_idx = (pct * S_CAST(uint64_t, read_sample_ct) + 99) / 100;
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
    if (unlikely(fclose_null(&outfile))) {
      goto ExportEigTgeno_ret_WRITE_FAIL;
    }
    logprintfww("--export eigt: %s written.\n", outname);
  }
  while (0) {
  ExportEigTgeno_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ExportEigTgeno_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  ExportEigTgeno_ret_PGR_FAIL:
    PgenErrPrintN(reterr);
    break;
  ExportEigTgeno_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  ExportEigTgeno_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 ExportEigTgeno_ret_1:
  CleanupThreads(&write_tg);
  CleanupThreads(&read_tg);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  pgfip->block_base = nullptr;
  return reterr;
}

PglErr Exportf(const uintptr_t* sample_include, const PedigreeIdInfo* piip, const uintptr_t* sex_nm, const uintptr_t* sex_male, const PhenoCol* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* allele_permute, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, const char* const* pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, const ExportfInfo* eip, const char* legacy_output_missing_pheno, const uint32_t* contig_lens, uintptr_t xheader_blen, InfoFlags info_flags, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_variant_id_slen, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, char input_missing_geno_char, char output_missing_geno_char, char legacy_output_missing_geno_char, uint32_t max_thread_ct, MakePlink2Flags make_plink2_flags, uintptr_t pgr_alloc_cacheline_ct, char* xheader, PgenFileInfo* pgfip, PgenReader* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
    const uint32_t sample_ctaw = BitCtToAlignedWordCt(sample_ct);
    const uint32_t sample_ctl = BitCtToWordCt(sample_ct);
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* sex_male_collapsed;
    if (unlikely(bigstack_alloc_u32(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
                 bigstack_alloc_w(sample_ctaw, &sex_male_collapsed))) {
      goto Exportf_ret_NOMEM;
    }
    FillCumulativePopcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, simple_pgrp, &pssi);
    CopyBitarrSubset(sex_male, sample_include, sample_ct, sex_male_collapsed);
    ZeroTrailingWords(sample_ctl, sex_male_collapsed);
    // possible todo: remove one of {sex_female_collapsed,
    // sex_nonfemale_collapsed}
    // unfortunately, the awkward sex_nonfemale is a slightly more convenient
    // representation (sex_nonfemale_cumulative_popcounts for chrY).
    uintptr_t* sex_female_collapsed = nullptr;
    ExportfFlags flags = eip->flags;
    if (flags & (kfExportfBcf | kfExportfBgen12 | kfExportfBgen13)) {
      if (unlikely(bigstack_alloc_w(sample_ctl, &sex_female_collapsed))) {
        goto Exportf_ret_NOMEM;
      }
      CopyBitarrSubset(sex_nm, sample_include, sample_ct, sex_female_collapsed);
      BitvecInvmask(sex_male_collapsed, sample_ctl, sex_female_collapsed);
    }
    uintptr_t* sex_nonfemale_collapsed = nullptr;
    if (flags & (kfExportfOxGen | kfExportfBgen11)) {
      uintptr_t* sex_nonfemale_tmp;
      if (unlikely(bigstack_alloc_w(sample_ctl, &sex_nonfemale_collapsed) ||
                   bigstack_alloc_w(raw_sample_ctl, &sex_nonfemale_tmp))) {
        goto Exportf_ret_NOMEM;
      }
      AlignedBitarrOrnotCopy(sex_male, sex_nm, raw_sample_ct, sex_nonfemale_tmp);
      CopyBitarrSubset(sex_nonfemale_tmp, sample_include, sample_ct, sex_nonfemale_collapsed);
      BigstackReset(sex_nonfemale_tmp);
    }
    uint32_t* sample_missing_geno_cts = nullptr;
    if (flags & (kfExportfOxGen | kfExportfHaps | kfExportfHapsLegend | kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13)) {
      if (unlikely(bigstack_alloc_u32(sample_ct, &sample_missing_geno_cts))) {
        goto Exportf_ret_NOMEM;
      }
    }
    if (flags & (kfExportf01 | kfExportf12)) {
      const uint32_t is_12 = (flags / kfExportf12) & 1;
      if (unlikely((flags & (kfExportfTped | kfExportfPed | kfExportfCompound)) && ((legacy_output_missing_geno_char == '1') || (ctou32(legacy_output_missing_geno_char) == '0' + (2 * is_12))))) {
        // Could move enforcement of this to initial command-line parsing if it
        // turns out to be a common mistake.
        logerrputs("Error: --export '01'/'12' conflicts with the --output-missing-genotype setting.\n");
        goto Exportf_ret_INVALID_CMDLINE;
      }
      // yes, this could be time- and memory-optimized...
      uintptr_t global_allele_ct = raw_variant_ct * 2;
      if (allele_idx_offsets) {
        global_allele_ct = allele_idx_offsets[raw_variant_ct];
      }
      const char** subst_allele_storage;
      if (unlikely(bigstack_alloc_kcp(global_allele_ct, &subst_allele_storage))) {
        goto Exportf_ret_NOMEM;
      }
      // '0' or '1'
      const char* new_ref_allele = &(g_one_char_strs[96 + 2 * is_12]);

      const char* new_alt_allele = &(new_ref_allele[2]);
      if (!allele_idx_offsets) {
        const char** subst_allele_storage_iter = subst_allele_storage;
        for (uint32_t variant_uidx = 0; variant_uidx != raw_variant_ct; ++variant_uidx) {
          *subst_allele_storage_iter++ = new_ref_allele;
          *subst_allele_storage_iter++ = new_alt_allele;
        }
      } else {
        uintptr_t variant_uidx_base = 0;
        uintptr_t cur_bits = variant_include[0];
        for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
          const uint32_t variant_uidx = BitIter1(variant_include, &variant_uidx_base, &cur_bits);
          const uintptr_t allele_idx_offset_base = allele_idx_offsets[variant_uidx];
          const uintptr_t allele_idx_offset_end = allele_idx_offsets[variant_uidx + 1];
          if (allele_idx_offset_end != allele_idx_offset_base + 2) {
            // Could permit this with ALT2 being encoded as '2'/'3', ALT3 being
            // encoded as '3'/'4', etc.  But let's not implement it
            // until/unless someone wants it.
            logerrprintf("Error: --export '01'/'12' cannot be used with multiallelic variants.  (Contact\nus if you want ALT2 to be encoded as '%c', ALT3 to be encoded as '%c', etc.)\n", '2' + is_12, '3' + is_12);
            goto Exportf_ret_INCONSISTENT_INPUT;
          }
          subst_allele_storage[allele_idx_offset_base] = new_ref_allele;
          subst_allele_storage[allele_idx_offset_base + 1] = new_alt_allele;
        }
      }
      allele_storage = subst_allele_storage;
    }
    if (flags & (kfExportfTypemask - kfExportfIndMajorBed - kfExportfVcf - kfExportfBcf - kfExportfOxGen - kfExportfBgen11 - kfExportfBgen12 - kfExportfBgen13 - kfExportfHaps - kfExportfHapsLegend - kfExportfAv - kfExportfA - kfExportfAD - kfExportfTped - kfExportfPed - kfExportfPhylip - kfExportfPhylipPhased - kfExportfEig - kfExportfEigt - kfExportfCompound)) {
      logerrputs("Error: Only VCF, BCF, oxford, bgen-1.x, haps, hapslegend, A, AD, Av, ped, tped,\ncompound-genotypes, phylip, phylip-phased, eig, eigt, and ind-major-bed output\nhave been implemented so far.\n");
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
      reterr = WriteMapOrBim(outname, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, variant_cms, variant_ct, max_allele_slen, exportf_delim, output_missing_geno_char, 0, max_thread_ct);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
      logputs("done.\n");
    }

    const AlleleCode* export_allele = nullptr;
    const char** export_allele_missing = nullptr;
    uint32_t max_export_allele_slen = max_allele_slen;
    if ((flags & (kfExportfAv | kfExportfA | kfExportfAD)) && (allele_permute || eip->export_allele_fname)) {
      AlleleCode* new_export_allele;
      if (unlikely(bigstack_alloc_ac(raw_variant_ct, &new_export_allele))) {
        goto Exportf_ret_NOMEM;
      }
      if (allele_permute) {
        if (!allele_idx_offsets) {
          // possible todo: SIMD
          for (uint32_t variant_uidx = 0; variant_uidx != raw_variant_ct; ++variant_uidx) {
            new_export_allele[variant_uidx] = allele_permute[variant_uidx * 2];
          }
        } else {
          for (uint32_t variant_uidx = 0; variant_uidx != raw_variant_ct; ++variant_uidx) {
            new_export_allele[variant_uidx] = allele_permute[allele_idx_offsets[variant_uidx]];
          }
        }
      } else {
        memset(new_export_allele, 0, raw_variant_ct * sizeof(AlleleCode));
      }
      if (eip->export_allele_fname) {
        reterr = ExportAlleleLoad(eip->export_allele_fname, variant_include, variant_ids, allele_idx_offsets, allele_storage, raw_variant_ct, variant_ct, max_variant_id_slen, input_missing_geno_char, max_thread_ct, &max_export_allele_slen, new_export_allele, &export_allele_missing);
        if (unlikely(reterr)) {
          goto Exportf_ret_1;
        }
      }
      export_allele = new_export_allele;
    }
    if (flags & kfExportfAv) {
      // multiallelic ok
      snprintf(outname_end, kMaxOutfnameExtBlen, ".traw");
      reterr = Export012Vmaj(outname, sample_include, sample_include_cumulative_popcounts, piip->sii.sample_ids, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, export_allele, export_allele_missing, variant_cms, sample_ct, piip->sii.max_sample_id_blen, variant_ct, max_export_allele_slen, exportf_delim, simple_pgrp);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & kfExportfIndMajorBed) {
      // multiallelic not ok, but already checked
      reterr = ExportIndMajorBed(sample_include, variant_include, allele_idx_offsets, allele_permute, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, max_thread_ct, pgr_alloc_cacheline_ct, pgfip, outname, outname_end);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & kfExportfOxGen) {
      // multiallelic really not ok
      reterr = ExportOxGen(sample_include, sample_include_cumulative_popcounts, sex_nonfemale_collapsed, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, sample_ct, variant_ct, max_allele_slen, max_thread_ct, flags, simple_pgrp, outname, outname_end, sample_missing_geno_cts);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & (kfExportfHaps | kfExportfHapsLegend)) {
      // multiallelic not ok
      reterr = ExportOxHapslegend(sample_include, sample_include_cumulative_popcounts, sex_male_collapsed, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, flags, simple_pgrp, outname, outname_end);
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
      reterr = ExportBgen11(outname, sample_include, sample_include_cumulative_popcounts, sex_nonfemale_collapsed, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, flags, pgr_alloc_cacheline_ct, pgfip, sample_missing_geno_cts);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    } else if (flags & (kfExportfBgen12 | kfExportfBgen13)) {
      // multiallelic ok... in the future, anyway.
      snprintf(outname_end, kMaxOutfnameExtBlen, ".bgen");
      reterr = ExportBgen13(outname, sample_include, sample_include_cumulative_popcounts, &(piip->sii), sex_male_collapsed, sex_female_collapsed, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, flags, eip->bgen_bits, idpaste_flags, id_delim, pgr_alloc_cacheline_ct, pgfip, sample_missing_geno_cts);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & (kfExportfOxGen | kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13 | kfExportfHaps | kfExportfHapsLegend)) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".sample");
      uint32_t y_ct = 0;
      const uint32_t y_code = cip->xymt_codes[kChrOffsetY];
      if ((!IsI32Neg(y_code)) && IsSet(cip->chr_mask, y_code)) {
        y_ct = CountChrVariantsUnsafe(variant_include, cip, y_code);
      }
      assert(PopcountWords(sample_include, raw_sample_ctl) == sample_ct);
      if (flags & kfExportfSampleV2) {
        reterr = ExportOxSampleV2(outname, sample_include, piip, sample_missing_geno_cts, sex_nm, sex_male, pheno_cols, pheno_names, sample_ct, pheno_ct, max_pheno_name_blen, variant_ct, y_ct, idpaste_flags, id_delim);
      } else {
        reterr = ExportOxSample(outname, sample_include, piip->sii.sample_ids, sample_missing_geno_cts, sex_nm, sex_male, pheno_cols, pheno_names, sample_ct, piip->sii.max_sample_id_blen, pheno_ct, max_pheno_name_blen, variant_ct, y_ct);
      }
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & kfExportfVcf) {
      // multiallelic ok
      reterr = ExportVcf(sample_include, sample_include_cumulative_popcounts, &(piip->sii), sex_male_collapsed, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pvar_info_reload, contig_lens, xheader_blen, info_flags, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, max_thread_ct, flags, eip->vcf_mode, idpaste_flags, id_delim, xheader, pgfip, simple_pgrp, outname, outname_end);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }
    if (flags & kfExportfBcf) {
      reterr = ExportBcf(sample_include, sample_include_cumulative_popcounts, &(piip->sii), sex_male_collapsed, sex_female_collapsed, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pvar_info_reload, contig_lens, xheader_blen, info_flags, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, max_thread_ct, flags, eip->vcf_mode, idpaste_flags, id_delim, xheader, pgfip, simple_pgrp, outname, outname_end);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }

    if (flags & (kfExportfA | kfExportfAD)) {
      // multiallelic ok
      snprintf(outname_end, kMaxOutfnameExtBlen, ".raw");
      reterr = Export012Smaj(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, variant_include, variant_ids, allele_idx_offsets, allele_storage, export_allele, export_allele_missing, legacy_output_missing_pheno, raw_sample_ct, sample_ct, pheno_ct, raw_variant_ct, variant_ct, max_export_allele_slen, (flags / kfExportfAD) & 1, (flags / kfExportfIncludeAlt) & 1, max_thread_ct, pgr_alloc_cacheline_ct, exportf_delim, pgfip);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }

    if (flags & kfExportfTped) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".tfam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteFam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, nullptr, legacy_output_missing_pheno, sample_ct, pheno_ct, exportf_delim);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
      logputs("done.\n");

      snprintf(outname_end, kMaxOutfnameExtBlen, ".tped");
      reterr = ExportTped(outname, sample_include, sample_include_cumulative_popcounts, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, variant_cms, sample_ct, variant_ct, max_allele_slen, exportf_delim, legacy_output_missing_geno_char, simple_pgrp);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }

    if (flags & (kfExportfPed | kfExportfCompound)) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".map");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteMapOrBim(outname, variant_include, cip, variant_bps, variant_ids, nullptr, nullptr, nullptr, variant_cms, variant_ct, 0, exportf_delim, '\0', 0, max_thread_ct);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
      logputs("done.\n");

      snprintf(outname_end, kMaxOutfnameExtBlen, ".ped");
      reterr = ExportPed(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, variant_include, allele_idx_offsets, allele_storage, legacy_output_missing_pheno, raw_sample_ct, sample_ct, pheno_ct, raw_variant_ct, variant_ct, max_allele_slen, (flags / kfExportfCompound) & 1, max_thread_ct, pgr_alloc_cacheline_ct, exportf_delim, legacy_output_missing_geno_char, pgfip);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }

    if (flags & (kfExportfPhylip | kfExportfPhylipPhased)) {
      reterr = ExportPhylip(sample_include, &(piip->sii), sex_male_collapsed, variant_include, cip, variant_bps, allele_idx_offsets, allele_storage, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, (flags / kfExportfPhylipPhased) & 1, (flags / kfExportfPhylipUsedSites) & 1, max_thread_ct, idpaste_flags, id_delim, pgr_alloc_cacheline_ct, pgfip, outname, outname_end);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }

    if (flags & (kfExportfEig | kfExportfEigt)) {
      // .snp first, since it'll error out on multiallelic variants and
      // multi-character allele codes.
      snprintf(outname_end, kMaxOutfnameExtBlen, ".snp");
      uint32_t v_hash;
      reterr = ExportEigSnp(outname, variant_include, cip, variant_bps, variant_ids, allele_idx_offsets, allele_storage, allele_permute, variant_cms, variant_ct, exportf_delim, &v_hash);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }

      snprintf(outname_end, kMaxOutfnameExtBlen, ".ind");
      uint32_t s_hash;
      reterr = ExportEigInd(outname, sample_include, &(piip->sii), sex_nm, sex_male, pheno_cols, sample_ct, pheno_ct, exportf_delim, idpaste_flags, id_delim, &s_hash);
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }

      snprintf(outname_end, kMaxOutfnameExtBlen, ".geno");
      if (flags & kfExportfEig) {
        reterr = ExportEigGeno(outname, sample_include, sample_include_cumulative_popcounts, variant_include, allele_idx_offsets, allele_permute, sample_ct, raw_variant_ct, variant_ct, s_hash, v_hash, max_thread_ct, pgr_alloc_cacheline_ct, pgfip);
      } else {
        reterr = ExportEigTgeno(outname, sample_include, variant_include, allele_idx_offsets, allele_permute, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, s_hash, v_hash, max_thread_ct, pgr_alloc_cacheline_ct, pgfip);
      }
      if (unlikely(reterr)) {
        goto Exportf_ret_1;
      }
    }

    if ((!(make_plink2_flags & kfMakeFam)) && (flags & kfExportfIndMajorBed)) {
      snprintf(outname_end, kMaxOutfnameExtBlen, ".fam");
      logprintfww5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = WriteFam(outname, sample_include, piip, sex_nm, sex_male, pheno_cols, nullptr, legacy_output_missing_pheno, sample_ct, pheno_ct, exportf_delim);
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
  Exportf_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  Exportf_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 Exportf_ret_1:
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
