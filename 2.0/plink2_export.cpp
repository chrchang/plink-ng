// This file is part of PLINK 2.00, copyright (C) 2005-2017 Shaun Purcell,
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

#include "bgzf.h"
#include "zstd/lib/zstd.h"

#include <time.h>

#ifdef __cplusplus
namespace plink2 {
#endif

// a few multithread globals
static uint16_t** g_bgen_geno_bufs = nullptr; // per-thread

static uint32_t g_sample_ct = 0;
static uint32_t g_calc_thread_ct = 0;
static uint32_t g_cur_block_write_ct = 0;
static pglerr_t g_error_ret = kPglRetSuccess;

// more multithread globals
static pgen_reader_t** g_pgr_ptrs = nullptr;
static uintptr_t** g_genovecs = nullptr;
static uintptr_t** g_dosage_presents = nullptr;
static dosage_t** g_dosage_val_bufs = nullptr;
static uint32_t* g_read_variant_uidx_starts = nullptr; // size calc_thread_ct

static unsigned char* g_writebufs[2] = {nullptr, nullptr};

static const uintptr_t* g_variant_include = nullptr;
static const chr_info_t* g_cip = nullptr;
static const uintptr_t* g_sample_include = nullptr;
static uint32_t* g_sample_include_cumulative_popcounts = nullptr;
static uintptr_t* g_sex_male_collapsed_interleaved = nullptr;
static const uintptr_t* g_variant_allele_idxs = nullptr;
static const alt_allele_ct_t* g_refalt1_select = nullptr;

static vul_t** g_thread_vecaligned_bufs = nullptr;
static uintptr_t** g_thread_write_genovecs = nullptr;
static uintptr_t** g_thread_write_dosagepresents = nullptr;
static dosage_t** g_thread_write_dosagevals = nullptr;

static uint32_t g_stride = 0;


FLAGSET_DEF_START()
  kfPlink2Write0,
  kfPlink2WriteSetHhMissing = (1 << 0),
  kfPlink2WriteSetMixedMtMissing = (1 << 1),
  kfPlink2WriteMeMissing = (1 << 2),
  kfPlink2WriteZeroCluster = (1 << 3),
  kfPlink2WriteFillRef = (1 << 4),
  kfPlink2WriteLateDosageErase = (1 << 5),
  // no need for sample_sort, determined by g_collapsed_sort_map != nullptr?
  kfPlink2WritePlink1 = (1 << 6)
FLAGSET_DEF_END(plink2_write_flags_t);
// todo: add .pgen-specific stuff


// assumes rawval is in [0, 163839]
static_assert(kDosageMid == 16384, "print_small_dosage() needs to be updated.");
char* print_small_dosage(uint32_t rawval, char* start) {
  // Instead of constant 5-digit precision, we print fewer digits whenever that
  // doesn't interfere with proper round-tripping.  I.e. we search for the
  // shortest string in
  //   ((n - 0.5)/16384, (n + 0.5)/16384).
  // E.g. 3277/16384 is 0.20001 when printed with 5-digit precision, but we'd
  // print that as 0.2 since that's still in (3276.5/16384, 3277.5/16384).
  *start++ = '0' + (rawval / 16384);
  rawval = rawval % 16384;
  if (!rawval) {
    return start;
  }
  *start++ = '.';
  // (rawval * 2) is in 32768ths
  // 32768 * 625 = 20480k

  const uint32_t range_top_20480k = (rawval * 2 + 1) * 625;
  // this is technically checking a half-open rather than a fully-open
  // interval, but that's fine since we never hit the boundary points
  if ((range_top_20480k % 2048) < 1250) {
    // when this is true, the four-decimal-place approximation is in the range
    // which round-trips back to our original number.
    const uint32_t four_decimal_places = range_top_20480k / 2048;
    return uitoa_trunc4(four_decimal_places, start);
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
    return uitoa_trunc4(last_four_digits, start);
  }
  return start;
}

#ifdef __arm__
  #error "Unaligned accesses in export_012_vmaj()."
#endif
pglerr_t export_012_vmaj(const char* outname, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const char* sample_ids, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, const double* variant_cms, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t variant_ct, uint32_t max_allele_slen, char exportf_delim, pgen_reader_t* simple_pgrp) {
  // todo: --recode-allele?
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    const uint32_t max_chr_blen = 1 + get_max_chr_slen(cip);
    char* chr_buf; // includes trailing tab
    char* writebuf;
    uintptr_t* genovec;
    // dosages are limited to 7 characters (x.yyyyy)
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
        bigstack_alloc_c(kMaxMediumLine + max_chr_blen + 2 * kMaxIdSlen + 48 + 2 * max_allele_slen + (8 * k1LU) * sample_ct, &writebuf) ||
        bigstack_alloc_ul(sample_ctl2, &genovec)) {
      goto export_012_vmaj_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    const uint32_t dosage_is_present = simple_pgrp->fi.gflags & kfPgenGlobalDosagePresent;
    uintptr_t* dosage_present = nullptr;
    dosage_t* dosage_vals = nullptr;
    if (dosage_is_present) {
      if (bigstack_alloc_ul(sample_ctl, &dosage_present) ||
          bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
        goto export_012_vmaj_ret_NOMEM;
      }
    }
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_012_vmaj_ret_OPEN_FAIL;
    }
    char* write_iter = strcpya(writebuf, (exportf_delim == '\t')? "CHR\tSNP\t(C)M\tPOS\tCOUNTED\tALT" : "CHR SNP (C)M POS COUNTED ALT");
    uint32_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_unsafe_ck(sample_include, &sample_uidx);
      *write_iter++ = exportf_delim;
      const char* fid_start = &(sample_ids[sample_uidx * max_sample_id_blen]);
      const char* fid_end = (const char*)rawmemchr(fid_start, exportf_delim);
      write_iter = memcpyax(write_iter, fid_start, fid_end - fid_start, '_');
      write_iter = strcpya(write_iter, &(fid_end[1]));
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto export_012_vmaj_ret_WRITE_FAIL;
      }
    }
    append_binary_eoln(&write_iter);
    LOGPRINTFWW5("--export A-transpose to %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_blen = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t cur_allele_ct = 2;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;

    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t variant_uidx = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
        *chr_name_end++ = exportf_delim;
        chr_blen = (uintptr_t)(chr_name_end - chr_buf);
      }
      write_iter = memcpya(write_iter, chr_buf, chr_blen);
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], exportf_delim);
      if (variant_cms) {
        write_iter = dtoa_g(variant_cms[variant_uidx], write_iter);
      } else {
        *write_iter++ = '0';
      }
      *write_iter++ = exportf_delim;
      write_iter = uint32toa_x(variant_bps[variant_uidx], exportf_delim, write_iter);
      // todo: multiallelic case
      uint32_t dosage_ct;
      uint32_t is_explicit_alt1;
      reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
      if (reterr) {
        if (reterr != kPglRetReadFail) {
          logprint("\n");
          logerrprint("Error: Malformed .pgen file.\n");
        }
        goto export_012_vmaj_ret_1;
      }
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[2 * variant_uidx];
      }
      if (!ref_allele_idx) {
        // we *usually* invert, since COUNTED = REF.
        genovec_invert_unsafe(sample_ct, genovec);
        biallelic_dosage16_invert(dosage_ct, dosage_vals);
      }
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
        variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      }
      char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
      write_iter = strcpyax(write_iter, cur_alleles[ref_allele_idx], exportf_delim);
      const uint32_t first_alt_idx = (ref_allele_idx == 0);
      write_iter = strcpya(write_iter, cur_alleles[first_alt_idx]);
      if (cur_allele_ct > 2) {
        for (uint32_t allele_idx = first_alt_idx + 1; allele_idx < cur_allele_ct; ++allele_idx) {
          if (allele_idx == ref_allele_idx) {
            continue;
          }
          if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
            goto export_012_vmaj_ret_WRITE_FAIL;
          }
          *write_iter++ = ',';
          write_iter = strcpya(write_iter, cur_alleles[allele_idx]);
        }
      }
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto export_012_vmaj_ret_WRITE_FAIL;
      }
      uint32_t widx = 0;
      uint32_t loop_len = kBitsPerWordD2;
      if (!dosage_ct) {
        while (1) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            loop_len = MOD_NZ(sample_ct, kBitsPerWordD2);
          }
          uintptr_t geno_word = genovec[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits < loop_len; ++sample_idx_lowbits) {
            *write_iter++ = exportf_delim;
            uintptr_t cur_geno = geno_word & 3;
            if (cur_geno != 3) {
              *write_iter++ = '0' + cur_geno;
            } else {
              write_iter = strcpya(write_iter, "NA");
            }
            geno_word >>= 2;
          }
          ++widx;
        }
      } else {
        dosage_t* dosage_vals_iter = dosage_vals;
        while (1) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            loop_len = MOD_NZ(sample_ct, kBitsPerWordD2);
          }
          uintptr_t geno_word = genovec[widx];
          uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits < loop_len; ++sample_idx_lowbits) {
            *write_iter++ = exportf_delim;
            if (dosage_present_hw & 1) {
              write_iter = print_small_dosage(*dosage_vals_iter++, write_iter);
            } else {
              uintptr_t cur_geno = geno_word & 3;
              if (cur_geno != 3) {
                *write_iter++ = '0' + cur_geno;
              } else {
                write_iter = strcpya(write_iter, "NA");
              }
            }
            geno_word >>= 2;
            dosage_present_hw >>= 1;
          }
          ++widx;
        }
      }
      append_binary_eoln(&write_iter);
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * ((uint64_t)variant_ct)) / 100;
      }
    }
    if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
      goto export_012_vmaj_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
  }
  while (0) {
  export_012_vmaj_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_012_vmaj_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_012_vmaj_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 export_012_vmaj_ret_1:
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

// more multithread globals
static uintptr_t* g_vmaj_readbuf = nullptr;

THREAD_FUNC_DECL transpose_to_smaj_read_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  const uintptr_t* variant_include = g_variant_include;
  // const uintptr_t* variant_allele_idxs = g_variant_allele_idxs;
  const alt_allele_ct_t* refalt1_select = g_refalt1_select;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uintptr_t* sample_include = g_sample_include;
  const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
  const uint32_t read_sample_ct = g_sample_ct;
  const uintptr_t read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
  uintptr_t prev_copy_ct = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_copy_ct = g_cur_block_write_ct;
    const uint32_t cur_idx_end = ((tidx + 1) * cur_block_copy_ct) / calc_thread_ct;
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    uint32_t cur_idx = (tidx * cur_block_copy_ct) / calc_thread_ct;
    uintptr_t* vmaj_readbuf_iter = &(g_vmaj_readbuf[(prev_copy_ct + cur_idx) * read_sample_ctaw2]);
    for (; cur_idx < cur_idx_end; ++cur_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      // todo: multiallelic case
      const pglerr_t reterr = pgr_read_refalt1_genovec_subset_unsafe(sample_include, sample_include_cumulative_popcounts, read_sample_ct, variant_uidx, pgrp, vmaj_readbuf_iter);
      if (reterr) {
        g_error_ret = reterr;
        break;
      }
      if (refalt1_select && (refalt1_select[2 * variant_uidx] == 1)) {
        genovec_invert_unsafe(read_sample_ct, vmaj_readbuf_iter);
        // don't need zero_trailing_quaters()
      }
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[read_sample_ctaw2]);
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    prev_copy_ct += cur_block_copy_ct;
    THREAD_BLOCK_FINISH(tidx);
  }
}

// more multithread globals
static uintptr_t* g_smaj_writebufs[2] = {nullptr, nullptr};
static uint32_t g_variant_ct = 0;
static uint32_t g_sample_batch_size = 0;
static uint32_t g_output_calc_thread_ct = 0;

THREAD_FUNC_DECL transpose_to_plink1_smaj_write_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t variant_ct = g_variant_ct;
  const uintptr_t variant_batch_ct = DIV_UP(variant_ct, kPglQuaterTransposeBatch);
  const uintptr_t variant_batch_word_ct = variant_batch_ct * kPglQuaterTransposeWords;
  const uint32_t calc_thread_ct = g_output_calc_thread_ct;
  const uint32_t variant_batch_idx_start = (((uint64_t)tidx) * variant_batch_ct) / calc_thread_ct;
  vul_t* vecaligned_buf = g_thread_vecaligned_bufs[tidx];
  uintptr_t variant_batch_idx_full_end = ((((uint64_t)tidx) + 1) * variant_batch_ct) / calc_thread_ct;
  uint32_t variant_idx_end;
  if (tidx + 1 < calc_thread_ct) {
    variant_idx_end = variant_batch_idx_full_end * kPglQuaterTransposeBatch;
  } else {
    variant_idx_end = variant_ct;
    if (variant_ct % kPglQuaterTransposeBatch) {
      --variant_batch_idx_full_end;
    }
  }
  const uint32_t thread_variant_ct = variant_idx_end - variant_batch_idx_start * kPglQuaterTransposeBatch;
  const uint32_t read_sample_ct = g_sample_ct;
  const uintptr_t read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
  const uintptr_t* vmaj_readbuf = g_vmaj_readbuf;
  uint32_t sample_widx = 0;
  uint32_t parity = 0;
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    uintptr_t variant_batch_idx = variant_batch_idx_start;
    uint32_t variant_batch_size = kPglQuaterTransposeBatch;
    const uintptr_t* vmaj_readbuf_iter = &(vmaj_readbuf[variant_batch_idx * kPglQuaterTransposeBatch * read_sample_ctaw2 + sample_widx]);
    const uint32_t sample_batch_size = g_sample_batch_size;
    uintptr_t* smaj_writebuf_start = &(g_smaj_writebufs[parity][variant_batch_idx * kPglQuaterTransposeWords]);
    uintptr_t* smaj_writebuf_iter = smaj_writebuf_start;
    while (1) {
      if (variant_batch_idx >= variant_batch_idx_full_end) {
        if (variant_batch_idx * kPglQuaterTransposeBatch >= variant_idx_end) {
          break;
        }
        variant_batch_size = variant_idx_end - variant_batch_idx * kPglQuaterTransposeBatch;
      }
      transpose_quaterblock(vmaj_readbuf_iter, read_sample_ctaw2, variant_batch_word_ct, variant_batch_size, sample_batch_size, smaj_writebuf_iter, vecaligned_buf);
      smaj_writebuf_iter = &(smaj_writebuf_iter[kPglQuaterTransposeWords]);
      vmaj_readbuf_iter = &(vmaj_readbuf_iter[variant_batch_size * read_sample_ctaw2]);
      ++variant_batch_idx;
    }
    smaj_writebuf_iter = smaj_writebuf_start;
    for (uint32_t sample_idx = 0; sample_idx < sample_batch_size; ++sample_idx) {
      // could fold this into transpose_quaterblock(), but I won't bother,
      // we're already saturating at ~3 threads
      pgr_plink2_to_plink1_inplace_unsafe(thread_variant_ct, smaj_writebuf_iter);
      zero_trailing_quaters(thread_variant_ct, smaj_writebuf_iter);
      smaj_writebuf_iter = &(smaj_writebuf_iter[variant_batch_word_ct]);
    }
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
    sample_widx += sample_batch_size / kBitsPerWordD2;
  }
}

pglerr_t export_ind_major_bed(const uintptr_t* orig_sample_include, const uintptr_t* variant_include, const uintptr_t* variant_allele_idxs, const alt_allele_ct_t* refalt1_select, uint32_t raw_sample_ct, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    // Possible special case: if the input file is a variant-major .bed, we do
    // not have enough memory to just load the whole file at once, and there
    // are more than ~20k samples, there can be a performance advantage to not
    // loading an entire variant at a time; we can use smaller fread calls and
    // reduce the number of (typically 4096 byte) disk blocks which need to be
    // read on each pass.  But let's get .pgen -> sample-major humming first.
    strcpy(outname_end, ".bed");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_ind_major_bed_ret_OPEN_FAIL;
    }
    if (fwrite_checked("l\x1b\0", 3, outfile)) {
      goto export_ind_major_bed_ret_WRITE_FAIL;
    }
    if (variant_ct && sample_ct) {
      const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
      uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      // todo: if only 1 pass is needed, and no subsetting is happening, this
      // saturates at ~4 threads?
      unsigned char* bigstack_end_mark = g_bigstack_end;
      // restrict multithread_load_init() to half of available workspace
      g_bigstack_end = &(g_bigstack_base[round_up_pow2(bigstack_left() / 2, kCacheline)]);
      unsigned char* main_loadbufs[2];
      pthread_t* threads;
      uint32_t read_block_size;
      if (multithread_load_init(variant_include, sample_ct, variant_ct, pgr_alloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, &g_genovecs, nullptr, nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
        goto export_ind_major_bed_ret_NOMEM;
      }
      g_bigstack_end = bigstack_end_mark;
      g_variant_include = variant_include;
      g_variant_allele_idxs = variant_allele_idxs;
      g_refalt1_select = refalt1_select;
      g_calc_thread_ct = calc_thread_ct;

      const uintptr_t variant_cacheline_ct = QUATERCT_TO_CLCT(variant_ct);
      uint32_t output_calc_thread_ct = MINV(calc_thread_ct, variant_cacheline_ct);
      // 4 still seems to be best in AVX2 case
      if (output_calc_thread_ct > 4) {
        output_calc_thread_ct = 4;
      }
      uintptr_t* sample_include;
      uint32_t* sample_include_cumulative_popcounts;
      if (bigstack_alloc_ul(raw_sample_ctl, &sample_include) ||
          bigstack_alloc_ui(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
          bigstack_alloc_vp(output_calc_thread_ct, &g_thread_vecaligned_bufs)) {
        goto export_ind_major_bed_ret_NOMEM;
      }
      for (uint32_t tidx = 0; tidx < output_calc_thread_ct; ++tidx) {
        g_thread_vecaligned_bufs[tidx] = (vul_t*)bigstack_alloc_raw(kPglQuaterTransposeBufbytes);
      }
      // each of the two write buffers should use <= 1/8 of the remaining
      // workspace
      const uintptr_t writebuf_cachelines_avail = bigstack_left() / (kCacheline * 8);
      uint32_t sample_batch_size = kPglQuaterTransposeBatch;
      if (variant_cacheline_ct * kPglQuaterTransposeBatch > writebuf_cachelines_avail) {
        sample_batch_size = round_down_pow2(writebuf_cachelines_avail / variant_cacheline_ct, kBitsPerWordD2);
        if (!sample_batch_size) {
          goto export_ind_major_bed_ret_NOMEM;
        }
      }
      g_smaj_writebufs[0] = (uintptr_t*)bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size);
      g_smaj_writebufs[1] = (uintptr_t*)bigstack_alloc_raw(variant_cacheline_ct * kCacheline * sample_batch_size);
      const uintptr_t readbuf_vecs_avail = (bigstack_left() / kCacheline) * kVecsPerCacheline;
      if (readbuf_vecs_avail < variant_ct) {
        goto export_ind_major_bed_ret_NOMEM;
      }
      uintptr_t read_sample_ctv2 = readbuf_vecs_avail / variant_ct;
      uint32_t read_sample_ct;
      if (read_sample_ctv2 >= QUATERCT_TO_VECCT(sample_ct)) {
        read_sample_ct = sample_ct;
      } else {
        read_sample_ct = read_sample_ctv2 * kQuatersPerVec;
      }
      uintptr_t read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
      uintptr_t* vmaj_readbuf = (uintptr_t*)bigstack_alloc_raw_rd(variant_ct * read_sample_ctaw2 * kBytesPerWord);
      g_variant_ct = variant_ct;
      g_output_calc_thread_ct = output_calc_thread_ct;
      g_error_ret = kPglRetSuccess;
      uint32_t sample_uidx_start = next_set_unsafe(orig_sample_include, 0);
      const uintptr_t variant_ct4 = QUATERCT_TO_BYTECT(variant_ct);
      const uintptr_t variant_ctaclw2 = variant_cacheline_ct * kWordsPerCacheline;
      const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
      const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
      const uint32_t pass_ct = 1 + (sample_ct - 1) / read_sample_ct;
      for (uint32_t pass_idx = 0; pass_idx < pass_ct; ++pass_idx) {
        memcpy(sample_include, orig_sample_include, raw_sample_ctl * sizeof(intptr_t));
        if (sample_uidx_start) {
          clear_bits_nz(0, sample_uidx_start, sample_include);
        }
        uint32_t sample_uidx_end;
        if (pass_idx + 1 == pass_ct) {
          read_sample_ct = sample_ct - pass_idx * read_sample_ct;
          read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
          sample_uidx_end = raw_sample_ct;
        } else {
          sample_uidx_end = jump_forward_set_unsafe(orig_sample_include, sample_uidx_start + 1, read_sample_ct);
          clear_bits_nz(sample_uidx_end, raw_sample_ct, sample_include);
        }
        fill_cumulative_popcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
        g_sample_include = sample_include;
        g_sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
        g_vmaj_readbuf = vmaj_readbuf;
        g_sample_ct = read_sample_ct;
        if (pass_idx) {
          pgfip->block_base = main_loadbufs[0];
          for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
            pgr_clear_ld_cache(g_pgr_ptrs[tidx]);
            g_pgr_ptrs[tidx]->fi.block_base = main_loadbufs[0];
            g_pgr_ptrs[tidx]->fi.block_offset = 0;
          }
        }
        uint32_t parity = 0;
        uint32_t read_block_idx = 0;
        uint32_t variant_idx = 0;
        uint32_t is_last_block = 0;
        uint32_t cur_read_block_size = read_block_size;
        uint32_t pct = 0;
        uint32_t next_print_idx = variant_ct / 100;
        putc_unlocked('\r', stdout);
        printf("--export ind-major-bed pass %u/%u: loading... 0%%", pass_idx + 1, pass_ct);
        fflush(stdout);
        while (1) {
          uintptr_t cur_block_write_ct = 0;
          if (!is_last_block) {
            while (read_block_idx < read_block_ct_m1) {
              cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
              if (cur_block_write_ct) {
                break;
              }
              ++read_block_idx;
            }
            if (read_block_idx == read_block_ct_m1) {
              cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
              cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
            }
            if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
              if (variant_idx) {
                join_threads2z(calc_thread_ct, 0, threads);
                g_cur_block_write_ct = 0;
                error_cleanup_threads2z(transpose_to_smaj_read_thread, calc_thread_ct, threads);
              }
              goto export_ind_major_bed_ret_THREAD_CREATE_FAIL;
            }
          }
          if (variant_idx) {
            join_threads2z(calc_thread_ct, is_last_block, threads);
            reterr = g_error_ret;
            if (reterr) {
              if (!is_last_block) {
                g_cur_block_write_ct = 0;
                error_cleanup_threads2z(transpose_to_smaj_read_thread, calc_thread_ct, threads);
              }
              if (reterr == kPglRetMalformedInput) {
                logprint("\n");
                logerrprint("Error: Malformed .pgen file.\n");
              }
              goto export_ind_major_bed_ret_1;
            }
          }
          if (!is_last_block) {
            g_cur_block_write_ct = cur_block_write_ct;
            compute_uidx_start_partition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
            for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
              g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
              g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
            }
            is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
            if (spawn_threads2z(transpose_to_smaj_read_thread, calc_thread_ct, is_last_block, threads)) {
              goto export_ind_major_bed_ret_THREAD_CREATE_FAIL;
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
            next_print_idx = (pct * ((uint64_t)variant_ct)) / 100;
          }

          ++read_block_idx;
          variant_idx += cur_block_write_ct;
          pgfip->block_base = main_loadbufs[parity];
        }
        // 2. Transpose and write.  (Could parallelize some of the transposing
        //    with the read loop, but since we can't write a single row until
        //    the read loop is done, and both write speed and write buffer
        //    space are bottlenecks, that can't be expected to help much.)
        g_sample_batch_size = sample_batch_size;
        parity = 0;
        is_last_block = 0;
        if (pct > 10) {
          fputs("\b \b", stdout);
        }
        fputs("\b\b\b\b\b\b\b\b\b\b\b\b\bwriting... 0%", stdout);
        fflush(stdout);
        pct = 0;
        uint32_t flush_sample_idx = 0;
        uint32_t flush_sample_idx_end = 0;
        next_print_idx = read_sample_ct / 100;
        while (1) {
          if (!is_last_block) {
            is_last_block = (flush_sample_idx_end + sample_batch_size >= read_sample_ct);
            if (is_last_block) {
              g_sample_batch_size = read_sample_ct - flush_sample_idx_end;
            }
            if (spawn_threads2z(transpose_to_plink1_smaj_write_thread, output_calc_thread_ct, is_last_block, threads)) {
              goto export_ind_major_bed_ret_THREAD_CREATE_FAIL;
            }
          }
          if (flush_sample_idx_end) {
            uintptr_t* smaj_writebuf_iter = g_smaj_writebufs[1 - parity];
            for (; flush_sample_idx < flush_sample_idx_end; ++flush_sample_idx) {
              fwrite_unlocked(smaj_writebuf_iter, variant_ct4, 1, outfile);
              smaj_writebuf_iter = &(smaj_writebuf_iter[variant_ctaclw2]);
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
              next_print_idx = (pct * ((uint64_t)read_sample_ct)) / 100;
            }
          }
          join_threads2z(output_calc_thread_ct, is_last_block, threads);
          if (ferror_unlocked(outfile)) {
            // may as well put this check when there are no threads to clean up
            goto export_ind_major_bed_ret_WRITE_FAIL;
          }
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
    if (fclose_null(&outfile)) {
      goto export_ind_major_bed_ret_WRITE_FAIL;
    }
    LOGPRINTFWW("--export ind-major-bed: %s written.\n", outname);
  }
  while (0) {
  export_ind_major_bed_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_ind_major_bed_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_ind_major_bed_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  export_ind_major_bed_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 export_ind_major_bed_ret_1:
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  bigstack_reset(bigstack_mark);
  return reterr;
}

static_assert(kDosageMid == 16384, "print_gen_dosage() needs to be updated.");
char* print_gen_dosage(uint32_t rawval, char* start) {
  // Similar to small_dosage_print(), but it's complicated by .gen import's
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
    return uitoa_trunc4(four_decimal_places, start);
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
    return uitoa_trunc4(last_four_digits, start);
  }
  return start;
}

// note that this is also in plink2_filter; may belong in plink2_common
static inline void incr_missing_row(const uintptr_t* genovec, uint32_t acc2_vec_ct, uintptr_t* missing_acc2) {
  const vul_t* genovvec = (const vul_t*)genovec;
  vul_t* missing_acc2v = (vul_t*)missing_acc2;
  const vul_t m1 = VCONST_UL(kMask5555);
  for (uint32_t vidx = 0; vidx < acc2_vec_ct; ++vidx) {
    const vul_t geno_vword = genovvec[vidx];
    const vul_t geno_vword_shifted_masked = vul_rshift(geno_vword, 1) & m1;
    missing_acc2v[vidx] = missing_acc2v[vidx] + (geno_vword & geno_vword_shifted_masked);
  }
}

pglerr_t export_ox_gen(const char* outname, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, uint32_t sample_ct, uint32_t variant_ct, uint32_t max_allele_slen, exportf_flags_t exportf_modifier, pgen_reader_t* simple_pgrp, uint32_t* sample_missing_geno_cts) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
    uintptr_t* genovec;
    uintptr_t* sex_male_collapsed_interleaved;
    uintptr_t* sex_male_collapsed_tmp;
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
    if (bigstack_alloc_ul(sample_ctl2, &genovec) ||

        bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_interleaved) ||
        bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_tmp)) {
      goto export_ox_gen_ret_NOMEM;
    }
    copy_bitarr_subset(sex_male, sample_include, sample_ct, sex_male_collapsed_tmp);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    zero_trailing_words(sample_ctl, sex_male_collapsed_tmp);
    fill_interleaved_mask_vec(sex_male_collapsed_tmp, sample_ctv, sex_male_collapsed_interleaved);
    bigstack_reset(sex_male_collapsed_tmp);

    // See load_sample_missing_cts in plink2_filter.cpp.
    // Yes, this is overkill, but the obvious alternative of incrementing
    // sample_missing_geno_cts[] when writing a missing call requires a bit of
    // custom chrY logic anyway.
    const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
    uintptr_t* missing_acc2;
    if (bigstack_calloc_ul(acc2_vec_ct * kWordsPerVec * 23, &missing_acc2)) {
      goto export_ox_gen_ret_NOMEM;
    }
    const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
    const uint32_t acc8_vec_ct = acc2_vec_ct * 4;
    uintptr_t* missing_acc4 = &(missing_acc2[acc2_vec_ct * kWordsPerVec]);
    uintptr_t* missing_acc8 = &(missing_acc4[acc4_vec_ct * kWordsPerVec]);
    uintptr_t* missing_acc32 = &(missing_acc8[acc8_vec_ct * kWordsPerVec]);

    uintptr_t* dosage_present = nullptr;
    dosage_t* dosage_vals = nullptr;
    if (simple_pgrp->fi.gflags & kfPgenGlobalDosagePresent) {
      const uint32_t multiallelic_present = (variant_allele_idxs != nullptr);
      if (bigstack_alloc_dosage(sample_ct * (1 + multiallelic_present), &dosage_vals) ||
          bigstack_alloc_ul(sample_ctl, &dosage_present)) {
        goto export_ox_gen_ret_NOMEM;
      }
    }
    const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
    // if no dosages, all genotypes are 6 bytes (missing = " 0 0 0")
    // with dosages, we print up to 5 digits past the decimal point, so 7 bytes
    //   + space for each number, 24 bytes max
    const uintptr_t max_geno_slen = 6 + (dosage_present != nullptr) * 18;
    char* chr_buf; // includes trailing space
    char* writebuf;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
        bigstack_alloc_c(kMaxMediumLine + max_chr_blen + kMaxIdSlen + 16 + 2 * max_allele_slen + max_geno_slen * sample_ct, &writebuf)) {
      goto export_ox_gen_ret_NOMEM;
    }
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_ox_gen_ret_OPEN_FAIL;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = writebuf;
    uint32_t variant_uidx = 0;
    uint32_t chr_blen = 0;

    // although we don't support --set-hh-missing, etc. here, we do still want
    // to be aware of chrY so we can exclude nonmales from the
    // sample_missing_geno_cts update there.
    uint32_t is_y = 0;

    uint32_t chr_fo_idx = UINT32_MAX;
    const int32_t y_code = cip->xymt_codes[kChrOffsetY];
    uint32_t chr_end = 0;
    uint32_t vidx_rem3 = 3;
    uint32_t vidx_rem15d3 = 5;
    uint32_t vidx_rem255d15 = 17;
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    const char hardcall_strs[] = " 1 0 0   0 1 0   0 0 1   0 0 0";
    const uint32_t ref_allele_second = !(exportf_modifier & kfExportfRefFirst);
    LOGPRINTFWW5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
        // Oxford spec doesn't seem to require spaces for .gen (only .sample),
        // but in practice spaces always seem to be used, and plink 1.9 doesn't
        // let you toggle this, so let's not worry about supporting tabs here
        *chr_name_end++ = ' ';
        chr_blen = (uintptr_t)(chr_name_end - chr_buf);
        is_y = (chr_idx == y_code);
      }
      write_iter = memcpya(write_iter, chr_buf, chr_blen);
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], ' ');
      write_iter = uint32toa_x(variant_bps[variant_uidx], ' ', write_iter);
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
        variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      }
      uint32_t is_explicit_alt1;
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx * 2];
        alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
      }
      // todo: multiallelic case
      uint32_t dosage_ct;
      reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
      if (reterr) {
        if (reterr != kPglRetReadFail) {
          logprint("\n");
          logerrprint("Error: Malformed .pgen file.\n");
        }
        goto export_ox_gen_ret_1;
      }
      if (ref_allele_idx + ref_allele_second == 1) {
        assert((!dosage_ct) || (!is_explicit_alt1));
        genovec_invert_unsafe(sample_ct, genovec);
        biallelic_dosage16_invert(dosage_ct, dosage_vals);
      }

      char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
      if (ref_allele_second) {
        write_iter = strcpyax(write_iter, cur_alleles[alt1_allele_idx], ' ');
        write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
      } else {
        write_iter = strcpyax(write_iter, cur_alleles[ref_allele_idx], ' ');
        write_iter = strcpya(write_iter, cur_alleles[alt1_allele_idx]);
      }
      uint32_t widx = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      if (!dosage_ct) {
        while (1) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t geno_word = genovec[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            write_iter = memcpya(write_iter, &(hardcall_strs[(geno_word & 3) * 8]), 6);
            geno_word >>= 2;
          }
          ++widx;
        }
      } else {
        const halfword_t* dosage_present_alias = (halfword_t*)dosage_present;
        const dosage_t* dosage_vals_iter = dosage_vals;
        if (!is_explicit_alt1) {
          while (1) {
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
                write_iter = memcpya(write_iter, &(hardcall_strs[(geno_word & 3) * 8]), 6);
                geno_word >>= 2;
              }
            } else {
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                if (dosage_present_hw & 1) {
                  const uint32_t dosage_int = *dosage_vals_iter++;
                  if (dosage_int <= kDosageMid) {
                    *write_iter++ = ' ';
                    write_iter = print_gen_dosage(kDosageMid - dosage_int, write_iter);
                    *write_iter++ = ' ';
                    write_iter = print_gen_dosage(dosage_int, write_iter);
                    write_iter = strcpya(write_iter, " 0");
                  } else {
                    assert(dosage_int <= kDosageMax);
                    write_iter = memcpyl3a(write_iter, " 0 ");
                    write_iter = print_gen_dosage(kDosageMax - dosage_int, write_iter);
                    *write_iter++ = ' ';
                    write_iter = print_gen_dosage(dosage_int - kDosageMid, write_iter);
                  }
                } else {
                  write_iter = memcpya(write_iter, &(hardcall_strs[(geno_word & 3) * 8]), 6);
                }
                geno_word >>= 2;
                dosage_present_hw >>= 1;
              }
            }
            ++widx;
          }
        } else {
          // todo
          // In multiallelic case, if ref/alt1 dosages sum to less than 2 (but
          // more than 0), we first internally rescale them to sum to 2, to
          // make .gen and bgen-1.1 export isomorphic, and bgen-1.2 export as
          // similar as possible.
          assert(0);
        }
      }
      append_binary_eoln(&write_iter);
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto export_ox_gen_ret_WRITE_FAIL;
      }
      if (is_y) {
        interleaved_mask_zero(sex_male_collapsed_interleaved, acc2_vec_ct, genovec);
      }
      incr_missing_row(genovec, acc2_vec_ct, missing_acc2);
      if (!(--vidx_rem3)) {
        unroll_zero_incr_2_4(acc2_vec_ct, missing_acc2, missing_acc4);
        vidx_rem3 = 3;
        if (!(--vidx_rem15d3)) {
          unroll_zero_incr_4_8(acc4_vec_ct, missing_acc4, missing_acc8);
          vidx_rem15d3 = 5;
          if (!(--vidx_rem255d15)) {
            unroll_zero_incr_8_32(acc8_vec_ct, missing_acc8, missing_acc32);
            vidx_rem255d15 = 17;
          }
        }
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * ((uint64_t)variant_ct)) / 100;
      }
    }
    if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
      goto export_ox_gen_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
    unroll_incr_2_4(missing_acc2, acc2_vec_ct, missing_acc4);
    unroll_incr_4_8(missing_acc4, acc4_vec_ct, missing_acc8);
    unroll_incr_8_32(missing_acc8, acc8_vec_ct, missing_acc32);
    uint32_t* scrambled_missing_cts = (uint32_t*)missing_acc32;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = scramble_2_4_8_32(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  export_ox_gen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_ox_gen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_ox_gen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 export_ox_gen_ret_1:
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

#ifdef __arm__
  #error "Unaligned accesses in export_ox_hapslegend()."
#endif
pglerr_t export_ox_hapslegend(const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male_collapsed, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, exportf_flags_t exportf_modifier, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  assert(sample_ct);
  assert(variant_ct);
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    const uint32_t just_haps = (exportf_modifier / kfExportfHaps) & 1;
    const uint32_t male_ct = popcount_longs(sex_male_collapsed, sample_ctl);
    if (xymt_is_nonempty(variant_include, cip, kChrOffsetY) && (male_ct != sample_ct)) {
      LOGERRPRINTF("Error: '--export haps%s' must exclude chrY unless the dataset is all-male.\n", just_haps? "" : "legend");
      goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
    }
    const uint32_t ref_allele_second = !(exportf_modifier & kfExportfRefFirst);
    const int32_t x_code = cip->xymt_codes[kChrOffsetX];
    const int32_t mt_code = cip->xymt_codes[kChrOffsetMT];
    char* chr_buf = nullptr;
    uint32_t is_x = 0;
    uint32_t is_haploid_or_mt = 0;
    uint32_t variant_uidx = next_set_unsafe(variant_include, 0);
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uintptr_t writebuf_alloc = 0;
    if (!just_haps) {
      // .legend doesn't have a chromosome column, so verify we only need to
      // export a single chromosome
      const uint32_t variant_uidx_start = variant_uidx;
      chr_fo_idx = get_variant_chr_fo_idx(cip, variant_uidx_start);
      chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      if ((chr_end != raw_variant_ct) && (popcount_bit_idx(variant_include, variant_uidx_start, chr_end) != variant_ct)) {
        logerrprint("Error: '--export hapslegend' does not support multiple chromosomes.\n");
        goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
      }
      const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      is_x = (chr_idx == x_code);
      is_haploid_or_mt = is_set(cip->haploid_mask, chr_idx) || (chr_idx == mt_code);
      strcpy(outname_end, ".legend");
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
        goto export_ox_hapslegend_ret_OPEN_FAIL;
      }
      char* writebuf;
      if (bigstack_alloc_c(kMaxMediumLine + kMaxIdSlen + 32 + 2 * max_allele_slen, &writebuf)) {
        goto export_ox_hapslegend_ret_NOMEM;
      }
      char* writebuf_flush = &(writebuf[kMaxMediumLine]);
      char* write_iter = strcpya(writebuf, "id position a0 a1" EOLN_STR);
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
        next_set_unsafe_ck(variant_include, &variant_uidx);
        write_iter = strcpyax(write_iter, variant_ids[variant_uidx], ' ');
        write_iter = uint32toa_x(variant_bps[variant_uidx], ' ', write_iter);
        if (refalt1_select) {
          ref_allele_idx = refalt1_select[variant_uidx * 2];
          alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
        }
        uintptr_t variant_allele_idx_base = variant_uidx * 2;
        if (variant_allele_idxs) {
          variant_allele_idx_base = variant_allele_idxs[variant_uidx];
        }
        char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
        if (ref_allele_second) {
          write_iter = strcpyax(write_iter, cur_alleles[alt1_allele_idx], ' ');
          write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
        } else {
          write_iter = strcpyax(write_iter, cur_alleles[ref_allele_idx], ' ');
          write_iter = strcpya(write_iter, cur_alleles[alt1_allele_idx]);
        }
        append_binary_eoln(&write_iter);
        if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
          goto export_ox_hapslegend_ret_WRITE_FAIL;
        }
      }
      if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
        goto export_ox_hapslegend_ret_WRITE_FAIL;
      }
      logprint("done.\n");
      variant_uidx = variant_uidx_start;
      bigstack_reset(writebuf);
    } else {
      const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
      if (bigstack_alloc_c(max_chr_blen, &chr_buf)) {
        goto export_ox_hapslegend_ret_NOMEM;
      }
      writebuf_alloc = max_chr_blen + kMaxIdSlen + 32 + 2 * max_allele_slen;
    }
    writebuf_alloc += kMaxMediumLine + (4 * k1LU) * sample_ct + kCacheline;
    const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    char* writebuf;
    uintptr_t* sex_male_collapsed_interleaved;
    uintptr_t* genovec;
    uintptr_t* phasepresent;
    uintptr_t* phaseinfo;
    if (bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_interleaved) ||
        bigstack_alloc_c(writebuf_alloc, &writebuf) ||
        bigstack_alloc_ul(sample_ctl2, &genovec) ||
        bigstack_alloc_ul(sample_ctl, &phasepresent) ||
        bigstack_alloc_ul(sample_ctl, &phaseinfo)) {
      goto export_ox_hapslegend_ret_NOMEM;
    }
    // sex_male_collapsed had trailing bits zeroed out
    fill_interleaved_mask_vec(sex_male_collapsed, sample_ctv, sex_male_collapsed_interleaved);
    // assumes little-endian
    // 3 = 1|0, not missing
    // 4..7 = male chrX
    // user's responsibility to split off PARs
    uint32_t genotext[7];
    genotext[0] = 0x20302030;
    genotext[2] = 0x20312031;
    genotext[4] = 0x202d2030;
    genotext[6] = 0x202d2031;
    if (ref_allele_second) {
      genotext[1] = 0x20302031;
      genotext[3] = 0x20312030;
    } else {
      genotext[1] = 0x20312030;
      genotext[3] = 0x20302031;
    }
#ifndef NDEBUG
    genotext[5] = 0x21475542; // "BUG!"
#endif
    uint32_t* cur_genotext = genotext;
    if (is_haploid_or_mt && (!is_x)) {
      cur_genotext = &(genotext[4]);
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = writebuf;
    strcpy(outname_end, ".haps");
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_ox_hapslegend_ret_OPEN_FAIL;
    }
    LOGPRINTFWW5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t chr_blen = 0;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
        *chr_name_end++ = ' ';
        chr_blen = (uintptr_t)(chr_name_end - chr_buf);
        is_x = (chr_idx == x_code);
        is_haploid_or_mt = is_set(cip->haploid_mask, chr_idx) || (chr_idx == mt_code);
        if ((!is_haploid_or_mt) || is_x) {
          cur_genotext = genotext;
        } else {
          cur_genotext = &(genotext[4]);
        }
      }
      uint32_t phasepresent_ct;
      reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, phasepresent, phaseinfo, &phasepresent_ct);
      if (reterr) {
        goto export_ox_hapslegend_ret_PGR_FAIL;
      }
      zero_trailing_quaters(sample_ct, genovec);
      if (!phasepresent_ct) {
        // phaseinfo is NOT cleared in this case
        fill_ulong_zero(sample_ctl, phaseinfo);
      }
      uint32_t genocounts[4];
      genovec_count_freqs_unsafe(genovec, sample_ct, genocounts);
      if (phasepresent_ct != genocounts[1]) {
        logprint("\n");
        LOGERRPRINTF("Error: '--export haps%s' must be used with a fully phased dataset.\n", just_haps? "" : "legend");
        goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
      } else if (genocounts[3]) {
        logprint("\n");
        LOGERRPRINTF("Error: '--export haps%s' cannot be used with missing genotype calls.\n", just_haps? "" : "legend");
        goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
      }
      if (is_haploid_or_mt) {
        // verify that there are no het haploids (treating MT as haploid here)
        if (is_x) {
          genovec_count_subset_freqs(genovec, sex_male_collapsed_interleaved, sample_ct, male_ct, genocounts);
        }
        if (genocounts[1]) {
          logprint("\n");
          LOGERRPRINTFWW("Error: '--export haps%s' cannot be used when heterozygous haploid/MT calls are present.%s\n", just_haps? "" : "legend", (is_x && (variant_bps[variant_uidx] <= 2781479))? " (Did you forget --split-par?)" : "");
          goto export_ox_hapslegend_ret_INCONSISTENT_INPUT;
        }
      }
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
        variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      }
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx * 2];
        alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
      }
      // this logic only works in the biallelic case
      if (ref_allele_second + ref_allele_idx == 1) {
        genovec_invert_unsafe(sample_ct, genovec);
        zero_trailing_quaters(sample_ct, genovec);
      }
      if (just_haps) {
        write_iter = memcpya(write_iter, chr_buf, chr_blen);
        write_iter = strcpyax(write_iter, variant_ids[variant_uidx], ' ');
        write_iter = uint32toa_x(variant_bps[variant_uidx], ' ', write_iter);
        char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
        if (ref_allele_second) {
          write_iter = strcpyax(write_iter, cur_alleles[alt1_allele_idx], ' ');
          write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
        } else {
          write_iter = strcpyax(write_iter, cur_alleles[ref_allele_idx], ' ');
          write_iter = strcpya(write_iter, cur_alleles[alt1_allele_idx]);
        }
        *write_iter++ = ' ';
      }
      uint32_t* write_iter_ui_alias = (uint32_t*)write_iter;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      if (!is_x) {
        while (1) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t genovec_word = genovec[widx];
          const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            const uintptr_t cur_geno = genovec_word & 3;
            *write_iter_ui_alias++ = cur_genotext[cur_geno + 2 * ((phaseinfo_halfword >> sample_idx_lowbits) & 1)];
            genovec_word >>= 2;
          }
          ++widx;
        }
      } else {
        while (1) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t genovec_word = genovec[widx];
          const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
          const uint32_t male_halfword = ((const halfword_t*)sex_male_collapsed)[widx];

          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            const uintptr_t cur_geno = genovec_word & 3;
            if (cur_geno == 2) {
              assert(!((phaseinfo_halfword >> sample_idx_lowbits) & 1));
            }
            *write_iter_ui_alias++ = cur_genotext[cur_geno + 2 * ((phaseinfo_halfword >> sample_idx_lowbits) & 1) + 4 * ((male_halfword >> sample_idx_lowbits) & 1)];
            genovec_word >>= 2;
          }
          ++widx;
        }
      }
      write_iter = (char*)write_iter_ui_alias;
      decr_append_binary_eoln(&write_iter);
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto export_ox_hapslegend_ret_WRITE_FAIL;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * ((uint64_t)variant_ct)) / 100;
      }
    }
    if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
      goto export_ox_hapslegend_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
  }
  while (0) {
  export_ox_hapslegend_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_ox_hapslegend_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_ox_hapslegend_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  export_ox_hapslegend_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  export_ox_hapslegend_ret_PGR_FAIL:
    if (reterr != kPglRetReadFail) {
      logprint("\n");
      logerrprint("Error: Malformed .pgen file.\n");
    }
  }
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

// more multithread globals
static uintptr_t** g_missing_acc2 = nullptr;
static uint32_t* g_variant_bytects[2] = {nullptr, nullptr};
static uint32_t g_ref_allele_second = 0;
static uint32_t g_bgen_compressed_buf_max = 0;
static uint32_t g_y_start = 0;
static uint32_t g_y_end = 0;

static const uint16_t bgen11_hardcall_usis[] = {32768, 0, 0, 0,
                                                0, 32768, 0, 0,
                                                0, 0, 32768, 0,
                                                0, 0, 0, 0};

THREAD_FUNC_DECL export_bgen11_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec = g_genovecs[tidx];
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
  const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
  const uint32_t acc8_vec_ct = acc2_vec_ct * 4;
  uintptr_t* missing_acc2 = g_missing_acc2[tidx];
  uintptr_t* missing_acc4 = &(missing_acc2[acc2_vec_ct * kWordsPerVec]);
  uintptr_t* missing_acc8 = &(missing_acc4[acc4_vec_ct * kWordsPerVec]);
  uintptr_t* missing_acc32 = &(missing_acc8[acc8_vec_ct * kWordsPerVec]);
  uintptr_t* dosage_present = g_dosage_presents? g_dosage_presents[tidx] : nullptr;
  dosage_t* dosage_vals = dosage_present? g_dosage_val_bufs[tidx] : nullptr;
  uint16_t* bgen_geno_buf = g_bgen_geno_bufs[tidx];
  const uintptr_t* variant_include = g_variant_include;
  const uintptr_t* sample_include = g_sample_include;
  const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed_interleaved = g_sex_male_collapsed_interleaved;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t sample_ctl2_m1 = QUATERCT_TO_WORDCT(sample_ct) - 1;
  const uint32_t bgen_geno_buf_blen = 6 * sample_ct;
  const uint32_t bgen_compressed_buf_max = g_bgen_compressed_buf_max;
  const alt_allele_ct_t* refalt1_select = g_refalt1_select;
  uint32_t is_y = 0;
  uint32_t y_thresh = g_y_start;
  const uint32_t y_end = g_y_end;
  const uint32_t ref_allele_second = g_ref_allele_second;
  uint32_t vidx_rem3 = 3;
  uint32_t vidx_rem15d3 = 5;
  uint32_t vidx_rem255d15 = 17;
  uint32_t ref_allele_idx = 0;
  uint32_t parity = 0;
  fill_ulong_zero(acc2_vec_ct * kWordsPerVec * 23, missing_acc2);
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t write_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t write_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* writebuf_iter = &(g_writebufs[parity][write_idx * ((uintptr_t)bgen_compressed_buf_max)]);
    uint32_t* variant_bytect_iter = &(g_variant_bytects[parity][write_idx]);
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    for (; write_idx < write_idx_end; ++write_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
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
        ref_allele_idx = refalt1_select[variant_uidx * 2];
      }
      uint32_t dosage_ct;
      uint32_t is_explicit_alt1;
      pglerr_t reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
      if (reterr) {
        g_error_ret = reterr;
        break;
      }
      if (ref_allele_idx + ref_allele_second == 1) {
        genovec_invert_unsafe(sample_ct, genovec);
        biallelic_dosage16_invert(dosage_ct, dosage_vals);
      }
      uint32_t widx = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint16_t* bgen_geno_buf_iter = bgen_geno_buf;
      if (!dosage_ct) {
        while (1) {
          if (widx >= sample_ctl2_m1) {
            if (widx > sample_ctl2_m1) {
              break;
            }
            inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
          }
          uintptr_t geno_word = genovec[widx];
          for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
            memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
            bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
            geno_word >>= 2;
          }
          ++widx;
        }
      } else {
        const halfword_t* dosage_present_alias = (halfword_t*)dosage_present;
        const dosage_t* dosage_vals_iter = dosage_vals;
        if (!is_explicit_alt1) {
          while (1) {
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
                memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
                bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
                geno_word >>= 2;
              }
            } else {
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                if (dosage_present_hw & 1) {
                  uint32_t dosage_int = *dosage_vals_iter++;
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
                  memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
                  bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
                }
                geno_word >>= 2;
                dosage_present_hw >>= 1;
              }
            }
            ++widx;
          }
        } else {
          // todo
          assert(0);
        }
      }
      uLongf compressed_blen = bgen_compressed_buf_max;
      if (compress(writebuf_iter, &compressed_blen, (const unsigned char*)bgen_geno_buf, bgen_geno_buf_blen)) {
        // is this actually possible?
        g_error_ret = kPglRetNomem;
        break;
      }
      *variant_bytect_iter++ = compressed_blen;
      writebuf_iter = &(writebuf_iter[bgen_compressed_buf_max]);
      if (is_y) {
        interleaved_mask_zero(sex_male_collapsed_interleaved, acc2_vec_ct, genovec);
      }
      incr_missing_row(genovec, acc2_vec_ct, missing_acc2);
      if (!(--vidx_rem3)) {
        unroll_zero_incr_2_4(acc2_vec_ct, missing_acc2, missing_acc4);
        vidx_rem3 = 3;
        if (!(--vidx_rem15d3)) {
          unroll_zero_incr_4_8(acc4_vec_ct, missing_acc4, missing_acc8);
          vidx_rem15d3 = 5;
          if (!(--vidx_rem255d15)) {
            unroll_zero_incr_8_32(acc8_vec_ct, missing_acc8, missing_acc32);
            vidx_rem255d15 = 17;
          }
        }
      }
    }
    if (is_last_block) {
      unroll_incr_2_4(missing_acc2, acc2_vec_ct, missing_acc4);
      unroll_incr_4_8(missing_acc4, acc4_vec_ct, missing_acc8);
      unroll_incr_8_32(missing_acc8, acc8_vec_ct, missing_acc32);
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

pglerr_t export_bgen11(const char* outname, const uintptr_t* sample_include, uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, exportf_flags_t exportf_modifier, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, uint32_t* sample_missing_geno_cts) {
  // isomorphic to export_ox_gen().
  assert(sample_ct);
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  // use gzip instead of zstd here.
  ZWRAP_useZSTDcompression(0);
  {
    const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
    const uint32_t max_chr_slen = get_max_chr_slen(cip);
    const uintptr_t bgen_compressed_buf_max = compressBound(6LU * sample_ct);
#ifdef __LP64__
    if (bgen_compressed_buf_max > UINT32_MAX) {
      logerrprint("Error: Too many samples for .bgen format.\n");
      goto export_bgen11_ret_INCONSISTENT_INPUT;
    }
#endif
    g_bgen_compressed_buf_max = bgen_compressed_buf_max;
    const uintptr_t writebuf_len = bgen_compressed_buf_max + 2 * max_allele_slen + 2 * kMaxIdSlen + 32;
    char* chr_buf;
    unsigned char* writebuf;
    uintptr_t* sex_male_collapsed_tmp;
    if (bigstack_alloc_c(max_chr_slen, &chr_buf) ||
        bigstack_alloc_uc(writebuf_len, &writebuf) ||
        bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male_collapsed_interleaved) ||
        bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_tmp)) {
      goto export_bgen11_ret_NOMEM;
    }
    copy_bitarr_subset(sex_male, sample_include, sample_ct, sex_male_collapsed_tmp);
    zero_trailing_words(BITCT_TO_WORDCT(sample_ct), sex_male_collapsed_tmp);
    fill_interleaved_mask_vec(sex_male_collapsed_tmp, sample_ctv, g_sex_male_collapsed_interleaved);
    bigstack_reset(sex_male_collapsed_tmp);

    const uintptr_t max_write_block_byte_ct = bigstack_left() / 4;
    uint32_t max_write_block_size = kPglVblockSize;
    while (1) {
      // limit each write buffer to 1/4 of remaining workspace
      if (((uint64_t)(bgen_compressed_buf_max + sizeof(int32_t))) * max_write_block_size <= max_write_block_byte_ct) {
        break;
      }
      if (max_write_block_size <= kBitsPerVec) {
        goto export_bgen11_ret_NOMEM;
      }
      max_write_block_size /= 2;
    }
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // seems to saturate around this point
    if (calc_thread_ct > 15) {
      calc_thread_ct = 15;
    }
    if (bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(g_writebufs[0])) ||
        bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(g_writebufs[1])) ||
        bigstack_alloc_ui(max_write_block_size, &(g_variant_bytects[0])) ||
        bigstack_alloc_ui(max_write_block_size, &(g_variant_bytects[1])) ||
        bigstack_alloc_ulp(calc_thread_ct, &g_missing_acc2) ||
        bigstack_alloc_usip(calc_thread_ct, &g_bgen_geno_bufs)) {
      goto export_bgen11_ret_NOMEM;
    }

    const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    const uintptr_t track_missing_cacheline_ct = VECCT_TO_CLCT(acc2_vec_ct * 23);
    const uintptr_t bgen_geno_cacheline_ct = DIV_UP(6 * sample_ct, (kCacheline * k1LU));
    const uintptr_t thread_xalloc_cacheline_ct = track_missing_cacheline_ct + bgen_geno_cacheline_ct;
    unsigned char* main_loadbufs[2];
    pthread_t* threads;
    uint32_t read_block_size;
    if (multithread_load_init(variant_include, sample_ct, variant_ct, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, pgfip, &calc_thread_ct, &g_genovecs, dosage_is_present? (&g_dosage_presents) : nullptr, dosage_is_present? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto export_bgen11_ret_NOMEM;
    }
    if (read_block_size > max_write_block_size) {
      read_block_size = max_write_block_size;
    }

    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_bgen11_ret_OPEN_FAIL;
    }
    // bgen 1.1 header
    // note that \xxx character constants are interpreted in octal, so \24 is
    // decimal 20, etc.
    memcpy(writebuf, "\24\0\0\0\24\0\0\0", 8);
    memcpy(&(writebuf[8]), &variant_ct, 4);
    memcpy(&(writebuf[12]), &sample_ct, 4);
    memcpy(&(writebuf[16]), "bgen\5\0\0\0", 8);
    if (fwrite_checked(writebuf, 24, outfile)) {
      goto export_bgen11_ret_WRITE_FAIL;
    }

    const uint32_t ref_allele_second = !(exportf_modifier & kfExportfRefFirst);
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_missing_acc2[tidx] = (uintptr_t*)bigstack_alloc_raw(track_missing_cacheline_ct * kCacheline);
      g_bgen_geno_bufs[tidx] = (uint16_t*)bigstack_alloc_raw(bgen_geno_cacheline_ct * kCacheline);
    }
    g_sample_ct = sample_ct;
    g_variant_include = variant_include;
    g_sample_include = sample_include;
    g_sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
    g_calc_thread_ct = calc_thread_ct;
    g_refalt1_select = refalt1_select;
    get_xymt_start_and_end(cip, kChrOffsetY, &g_y_start, &g_y_end);
    g_ref_allele_second = ref_allele_second;
    g_cip = cip;

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
    const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
    const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t write_variant_uidx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_slen = 0;

    uint32_t prev_block_write_ct = 0;
    uint32_t variant_idx = 0;
    uint32_t is_last_block = 0;
    uint32_t cur_read_block_size = read_block_size;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    LOGPRINTFWW5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    while (1) {
      uintptr_t cur_block_write_ct = 0;
      if (!is_last_block) {
        while (read_block_idx < read_block_ct_m1) {
          cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
          if (cur_block_write_ct) {
            break;
          }
          ++read_block_idx;
        }
        if (read_block_idx == read_block_ct_m1) {
          cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
          cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
        }
        if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
          if (variant_idx) {
            join_threads2z(calc_thread_ct, 0, threads);
            g_cur_block_write_ct = 0;
            error_cleanup_threads2z(export_bgen11_thread, calc_thread_ct, threads);
          }
          goto export_bgen11_ret_READ_FAIL;
        }
      }
      if (variant_idx) {
        join_threads2z(calc_thread_ct, is_last_block, threads);
        reterr = g_error_ret;
        if (reterr) {
          if (!is_last_block) {
            g_cur_block_write_ct = 0;
            error_cleanup_threads2z(export_bgen11_thread, calc_thread_ct, threads);
          }
          if (reterr == kPglRetMalformedInput) {
            logprint("\n");
            logerrprint("Error: Malformed .pgen file.\n");
          }
          goto export_bgen11_ret_1;
        }
      }
      if (!is_last_block) {
        g_cur_block_write_ct = cur_block_write_ct;
        compute_uidx_start_partition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
        for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
          g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
          g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
        }
        is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
        if (spawn_threads2z(export_bgen11_thread, calc_thread_ct, is_last_block, threads)) {
          goto export_bgen11_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const unsigned char* compressed_data_iter = g_writebufs[parity];
        const uint32_t* variant_bytect_iter = g_variant_bytects[parity];
        for (uint32_t variant_bidx = 0; variant_bidx < prev_block_write_ct; ++variant_bidx, ++write_variant_uidx) {
          next_set_unsafe_ck(variant_include, &write_variant_uidx);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
            chr_slen = (uintptr_t)(chr_name_end - chr_buf);
          }
          const char* cur_variant_id = variant_ids[write_variant_uidx];
          const uint32_t id_slen = strlen(cur_variant_id);
          memcpy(&(writebuf[6]), &id_slen, 4);
          // deliberately clobber top two bytes
          unsigned char* writebuf_iter = (unsigned char*)memcpya(&(writebuf[8]), cur_variant_id, id_slen);
          memcpy(writebuf_iter, &chr_slen, 4);
          writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[2]), chr_buf, chr_slen);
          memcpy(writebuf_iter, &(variant_bps[write_variant_uidx]), 4);
          writebuf_iter = &(writebuf_iter[4]);
          uintptr_t variant_allele_idx_base = write_variant_uidx * 2;
          if (variant_allele_idxs) {
            variant_allele_idx_base = variant_allele_idxs[write_variant_uidx];
          }
          char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
          if (refalt1_select) {
            ref_allele_idx = refalt1_select[write_variant_uidx * 2];
            alt1_allele_idx = refalt1_select[write_variant_uidx * 2 + 1];
          }
          const char* first_allele;
          const char* second_allele;
          if (ref_allele_second) {
            first_allele = cur_alleles[alt1_allele_idx];
            second_allele = cur_alleles[ref_allele_idx];
          } else {
            first_allele = cur_alleles[ref_allele_idx];
            second_allele = cur_alleles[alt1_allele_idx];
          }
          uint32_t allele_slen = strlen(first_allele);
          memcpy(writebuf_iter, &allele_slen, 4);
          writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[4]), first_allele, allele_slen);
          allele_slen = strlen(second_allele);
          memcpy(writebuf_iter, &allele_slen, 4);
          writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[4]), second_allele, allele_slen);
          const uint32_t cur_variant_bytect = *variant_bytect_iter++;
          memcpy(writebuf_iter, &cur_variant_bytect, 4);
          writebuf_iter = &(writebuf_iter[4]);
          memcpy(writebuf_iter, compressed_data_iter, cur_variant_bytect);
          writebuf_iter = &(writebuf_iter[cur_variant_bytect]);
          compressed_data_iter = &(compressed_data_iter[bgen_compressed_buf_max]);
          if (fwrite_checked(writebuf, writebuf_iter - writebuf, outfile)) {
            if (variant_idx < variant_ct) {
              join_threads2z(calc_thread_ct, is_last_block, threads);
              if (!is_last_block) {
                g_cur_block_write_ct = 0;
                error_cleanup_threads2z(export_bgen11_thread, calc_thread_ct, threads);
              }
            }
            goto export_bgen11_ret_WRITE_FAIL;
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
        next_print_variant_idx = (pct * ((uint64_t)variant_ct)) / 100;
      }
      ++read_block_idx;
      prev_block_write_ct = cur_block_write_ct;
      variant_idx += cur_block_write_ct;
      pgfip->block_base = main_loadbufs[parity];
    }
    if (fclose_null(&outfile)) {
      goto export_bgen11_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
    const uint32_t sample_ctav2 = acc2_vec_ct * kQuatersPerVec;
    const uintptr_t acc32_offset = acc2_vec_ct * (7 * k1LU * kWordsPerVec);
    uint32_t* scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx < calc_thread_ct; ++tidx) {
      const uint32_t* thread_scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii < sample_ctav2; ++uii) {
        scrambled_missing_cts[uii] += thread_scrambled_missing_cts[uii];
      }
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = scramble_2_4_8_32(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  export_bgen11_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_bgen11_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_bgen11_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  export_bgen11_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
#ifdef __LP64__
  export_bgen11_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
#endif
  export_bgen11_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 export_bgen11_ret_1:
  ZWRAP_useZSTDcompression(1);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

/*
static uintptr_t** g_phasepresents = nullptr;
static uintptr_t** g_phaseinfos = nullptr;
static uintptr_t** g_dphase_presents = nullptr;

static uint32_t g_bgen_bit_precision = 0;
static uint32_t g_bgen_uncompressed_buf_max = 0;

// memcpy(target,
//   &(g_bgen_hardcall_write[cur_geno * biallelic_diploid_byte_ct]),
//   biallelic_diploid_byte_ct) should do the right thing
static unsigned char* g_bgen_hardcall_write = nullptr;

THREAD_FUNC_DECL export_bgen13_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec = g_genovecs[tidx];
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
  const uint32_t acc4_vec_ct = acc2_vec_ct * 2;
  const uint32_t acc8_vec_ct = acc2_vec_ct * 4;
  uintptr_t* missing_acc2 = g_missing_acc2[tidx];
  uintptr_t* missing_acc4 = &(missing_acc2[acc2_vec_ct * kWordsPerVec]);
  uintptr_t* missing_acc8 = &(missing_acc4[acc4_vec_ct * kWordsPerVec]);
  uintptr_t* missing_acc32 = &(missing_acc8[acc8_vec_ct * kWordsPerVec]);
  uintptr_t* phasepresent = nullptr;
  uintptr_t* phaseinfo = nullptr;
  if (g_phasepresents) {
    phasepresent = g_phasepresents[tidx];
    phaseinfo = g_phaseinfos[tidx];
  }
  uintptr_t* dosage_present = g_dosage_presents? g_dosage_presents[tidx] : nullptr;
  uintptr_t* dphase_present = g_dphase_presents? g_dphase_presents[tidx] : nullptr;
  dosage_t* dosage_vals = dosage_present? g_dosage_val_bufs[tidx] : nullptr;
  unsigned char* uncompressed_bgen_geno_buf = g_thread_wkspaces[tidx];
  const uintptr_t* variant_include = g_variant_include;
  const chr_info_t* cip = g_cip;
  const uintptr_t* sample_include = g_sample_include;
  const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
  const uintptr_t* sex_male_collapsed_interleaved = g_sex_male_collapsed_interleaved;
  const unsigned char* bgen_diploid_hardcall_write = g_bgen_hardcall_write;
  const uint32_t calc_thread_ct = g_calc_thread_ct;
  const uint32_t sample_ctl2_m1 = QUATERCT_TO_WORDCT(sample_ct) - 1;
  const uint32_t bit_precision = g_bgen_bit_precision;
  const uint32_t bytes_per_prob = DIV_UP(bit_precision, CHAR_BIT);

  // note that this applies to both unphased and phased output, for different
  // reasons
  const uint32_t biallelic_diploid_byte_ct = 2 * bytes_per_prob;
  const unsigned char* bgen_haploid_hardcall_write = &(bgen_diploid_hardcall_write[4 * biallelic_diploid_byte_ct]);
  ;;;

  const uint32_t bgen_uncompressed_buf_max = g_bgen_uncompressed_buf_max;
  const uint32_t bgen_compressed_buf_max = g_bgen_compressed_buf_max;
  const alt_allele_ct_t* refalt1_select = g_refalt1_select;
  uint32_t chr_fo_idx = UINT32_MAX; // deliberate overflow
  uint32_t chr_end = 0;
  uint32_t is_x = 0;
  uint32_t is_y = 0;

  uint32_t is_haploid_or_mt = 0; // includes chrX and chrY
  // for bgen-1.2/1.3 and VCF/BCF export, MT ploidy is 1 unless the call is
  //   heterozygous (i.e. it's treated the same way as an ordinary haploid
  //   chromosome); similarly for chrX male ploidy
  // for bgen-1.2/1.3, chrY female (but not unknown-sex) ploidy is 0 when
  //   genotype is missing

  const uint32_t ref_allele_second = g_ref_allele_second;
  uint32_t vidx_rem3 = 3;
  uint32_t vidx_rem15d3 = 5;
  uint32_t vidx_rem255d15 = 17;
  uint32_t ref_allele_idx = 0;
  uint32_t parity = 0;
  fill_ulong_zero(acc2_vec_ct * kWordsPerVec * 23, missing_acc2);
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uintptr_t cur_block_write_ct = g_cur_block_write_ct;
    uint32_t write_idx = (tidx * cur_block_write_ct) / calc_thread_ct;
    const uint32_t write_idx_end = ((tidx + 1) * cur_block_write_ct) / calc_thread_ct;
    unsigned char* writebuf_iter = &(g_writebufs[parity][write_idx * ((uintptr_t)bgen_compressed_buf_max)]);
    uint32_t* variant_bytect_iter = &(g_variant_bytects[parity][write_idx]);
    uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
    for (; write_idx < write_idx_end; ++write_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        is_x = (chr_idx == x_code);
        is_y = (chr_idx == y_code);
        is_haploid_or_mt = is_set(cip->haploid_mask, chr_idx) || (chr_idx == mt_code);
      }
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx * 2];
      }
      // todo: export phase info
      // todo: multiallelic cases
      uint32_t dosage_ct;
      uint32_t is_explicit_alt1;
      pglerr_t reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
      if (reterr) {
        g_error_ret = reterr;
        break;
      }
      if (ref_allele_idx + ref_allele_second == 1) {
        genovec_invert_unsafe(sample_ct, genovec);
        biallelic_dosage16_invert(dosage_ct, dosage_vals);
      }
      unsigned char* bgen_geno_buf_iter = uncompressed_bgen_geno_buf;
      // 4 bytes: # of samples
      // 2 bytes: # of alleles
      // 1 byte: minimum ploidy
      // 1 byte: maximum ploidy
      // sample_ct bytes: high bit = missing, low bits = ploidy
      // 1 byte: is_phased
      // 1 byte: bit_precision
      bgen_geno_buf_iter = memcpya(bgen_geno_buf_iter, &sample_ct, 4);
      uint32_t widx = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      if (!dosage_ct) {
        if (!is_haploid_or_mt) {
          // 2 alleles, min ploidy == max ploidy == 2
          *((uint32_t*)bgen_geno_buf_iter)++ = 0x2020002;
          unsigned char* sample_ploidy_and_missingness = bgen_geno_buf_iter;
          bgen_geno_buf_iter = memseta(bgen_geno_buf_iter, 2, sample_ct);
          *bgen_geno_buf_iter++ = 0; // not phased
          *bgen_geno_buf_iter++ = bit_precision;
          while (1) {
            if (widx >= sample_ctl2_m1) {
              if (widx > sample_ctl2_m1) {
                break;
              }
              inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
            }
            uintptr_t geno_word = genovec[widx];
            for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
              const uintptr_t cur_geno = geno_word & 3;
              bgen_geno_buf_iter = memcpya(bgen_geno_buf_iter, &(bgen_hardcall_write[cur_geno * biallelic_diploid_byte_ct]), biallelic_diploid_byte_ct);
              if (cur_geno == 3) {
                // maybe handle this in a different loop?
                sample_ploidy_and_missingness_iter[sample_idx_lowbits] = 130;
              }
              geno_word >>= 2;
            }
            ++widx;
            sample_ploidy_and_missingness_iter = &(sample_ploidy_and_missingness_iter[kBitsPerWordD2]);
          }
        } else if (is_x) {
        } else {
          // ...
        }
      } else {
        const halfword_t* dosage_present_alias = (halfword_t*)dosage_present;
        const dosage_t* dosage_vals_iter = dosage_vals;
        if (!is_explicit_alt1) {
          while (1) {
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
                memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
                bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
                geno_word >>= 2;
              }
            } else {
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                if (dosage_present_hw & 1) {
                  uint32_t dosage_int = *dosage_vals_iter++;
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
                  memcpy(bgen_geno_buf_iter, &(bgen11_hardcall_usis[(geno_word & 3) * 4]), 6);
                  bgen_geno_buf_iter = &(bgen_geno_buf_iter[3]);
                }
                geno_word >>= 2;
                dosage_present_hw >>= 1;
              }
            }
            ++widx;
          }
        } else {
          // todo
          assert(0);
        }
      }
      uLongf compressed_blen = bgen_compressed_buf_max;
      if (compress(writebuf_iter, &compressed_blen, (const unsigned char*)bgen_geno_buf, bgen_geno_buf_blen)) {
        // is this actually possible?
        g_error_ret = kPglRetNomem;
        break;
      }
      *variant_bytect_iter++ = compressed_blen;
      writebuf_iter = &(writebuf_iter[bgen_compressed_buf_max]);
      if (is_y) {
        interleaved_mask_zero(sex_male_collapsed_interleaved, acc2_vec_ct, genovec);
      }
      incr_missing_row(genovec, acc2_vec_ct, missing_acc2);
      if (!(--vidx_rem3)) {
        unroll_zero_incr_2_4(acc2_vec_ct, missing_acc2, missing_acc4);
        vidx_rem3 = 3;
        if (!(--vidx_rem15d3)) {
          unroll_zero_incr_4_8(acc4_vec_ct, missing_acc4, missing_acc8);
          vidx_rem15d3 = 5;
          if (!(--vidx_rem255d15)) {
            unroll_zero_incr_8_32(acc8_vec_ct, missing_acc8, missing_acc32);
            vidx_rem255d15 = 17;
          }
        }
      }
    }
    if (is_last_block) {
      unroll_incr_2_4(missing_acc2, acc2_vec_ct, missing_acc4);
      unroll_incr_4_8(missing_acc4, acc4_vec_ct, missing_acc8);
      unroll_incr_8_32(missing_acc8, acc8_vec_ct, missing_acc32);
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
    parity = 1 - parity;
  }
}

pglerr_t export_bgen13(const char* outname, const uintptr_t* sample_include, uint32_t* sample_include_cumulative_popcounts, const uintptr_t* sex_male, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, uint32_t sample_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_thread_ct, exportf_flags_t exportf_modifier, uint32_t exportf_bits, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, uint32_t* sample_missing_geno_cts) {
  assert(sample_ct);
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  // use gzip iff v1.2
  if (exportf_modifier & kfExportfBgen12) {
    ZWRAP_useZSTDcompression(0);
  }
  {
    if (exportf_bits > 16) {
      logerrprint("Error: bits= parameter is currently limited to 16.  (This is sufficient to\ncapture all information in a .pgen file.)\n");
      reterr = kPglRetNotYetSupported;
      goto export_bgen13_ret_1;
    }
    if (pgfip->gflags & (kfPgenGlobalHardcallPhasePresent | kfPgenGlobalDosagePhasePresent)) {
      logerrprint("Error: Export of phase information to .bgen files is currently under\ndevelopment.\n");
      reterr = kPglRetNotYetSupported;
      goto export_bgen13_ret_1;
    }
    const uint32_t sample_ctv = BITCT_TO_VECCT(sample_ct);
    const uint32_t max_chr_slen = get_max_chr_slen(cip);
    const uintptr_t bgen_compressed_buf_max = compressBound(6LU * sample_ct);
#ifdef __LP64__
    if (bgen_compressed_buf_max > UINT32_MAX) {
      logerrprint("Error: Too many samples for .bgen format.\n");
      goto export_bgen13_ret_INCONSISTENT_INPUT;
    }
#endif
    g_bgen_compressed_buf_max = bgen_compressed_buf_max;
    const uintptr_t writebuf_len = bgen_compressed_buf_max + 2 * max_allele_slen + 2 * kMaxIdSlen + 32;
    char* chr_buf;
    unsigned char* writebuf;
    uintptr_t* sex_male_collapsed_tmp;
    if (bigstack_alloc_c(max_chr_slen, &chr_buf) ||
        bigstack_alloc_uc(writebuf_len, &writebuf) ||
        bigstack_alloc_ul(sample_ctv * kWordsPerVec, &g_sex_male_collapsed_interleaved) ||
        bigstack_alloc_ul(sample_ctv * kWordsPerVec, &sex_male_collapsed_tmp)) {
      goto export_bgen13_ret_NOMEM;
    }
    copy_bitarr_subset(sex_male, sample_include, sample_ct, sex_male_collapsed_tmp);
    zero_trailing_words(BITCT_TO_WORDCT(sample_ct), sex_male_collapsed_tmp);
    fill_interleaved_mask_vec(sex_male_collapsed_tmp, sample_ctv, g_sex_male_collapsed_interleaved);
    bigstack_reset(sex_male_collapsed_tmp);

    const uintptr_t max_write_block_byte_ct = bigstack_left() / 4;
    uint32_t max_write_block_size = kPglVblockSize;
    while (1) {
      // limit each write buffer to 1/4 of remaining workspace
      if (((uint64_t)(bgen_compressed_buf_max + sizeof(int32_t))) * max_write_block_size <= max_write_block_byte_ct) {
        break;
      }
      if (max_write_block_size <= kBitsPerVec) {
        goto export_bgen13_ret_NOMEM;
      }
      max_write_block_size /= 2;
    }
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    // seems to saturate around this point
    if (calc_thread_ct > 15) {
      calc_thread_ct = 15;
    }
    if (bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(g_writebufs[0])) ||
        bigstack_alloc_uc(bgen_compressed_buf_max * max_write_block_size, &(g_writebufs[1])) ||
        bigstack_alloc_ui(max_write_block_size, &(g_variant_bytects[0])) ||
        bigstack_alloc_ui(max_write_block_size, &(g_variant_bytects[1])) ||
        bigstack_alloc_ulp(calc_thread_ct, &g_missing_acc2) ||
        bigstack_alloc_usip(calc_thread_ct, &g_bgen_geno_bufs)) {
      goto export_bgen13_ret_NOMEM;
    }

    const uint32_t acc2_vec_ct = QUATERCT_TO_VECCT(sample_ct);
    const uint32_t dosage_is_present = pgfip->gflags & kfPgenGlobalDosagePresent;
    const uintptr_t track_missing_cacheline_ct = VECCT_TO_CLCT(acc2_vec_ct * 23);
    const uintptr_t bgen_geno_cacheline_ct = DIV_UP(6 * sample_ct, (kCacheline * k1LU));
    const uintptr_t thread_xalloc_cacheline_ct = track_missing_cacheline_ct + bgen_geno_cacheline_ct;
    unsigned char* main_loadbufs[2];
    pthread_t* threads;
    uint32_t read_block_size;
    if (multithread_load_init(variant_include, sample_ct, raw_variant_ct, pgr_alloc_cacheline_ct, thread_xalloc_cacheline_ct, 0, pgfip, &calc_thread_ct, &g_genovecs, dosage_is_present? (&g_dosage_presents) : nullptr, dosage_is_present? (&g_dosage_val_bufs) : nullptr, &read_block_size, main_loadbufs, &threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto export_bgen13_ret_NOMEM;
    }
    if (read_block_size > max_write_block_size) {
      read_block_size = max_write_block_size;
    }

    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_bgen13_ret_OPEN_FAIL;
    }
    // bgen 1.1 header
    // note that \xxx character constants are interpreted in octal, so \24 is
    // decimal 20, etc.
    memcpy(writebuf, "\24\0\0\0\24\0\0\0", 8);
    memcpy(&(writebuf[8]), &variant_ct, 4);
    memcpy(&(writebuf[12]), &sample_ct, 4);
    memcpy(&(writebuf[16]), "bgen\5\0\0\0", 8);
    if (fwrite_checked(writebuf, 24, outfile)) {
      goto export_bgen13_ret_WRITE_FAIL;
    }

    const uint32_t ref_allele_second = !(exportf_modifier & kfExportfRefFirst);
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_missing_acc2[tidx] = (uintptr_t*)bigstack_alloc_raw(track_missing_cacheline_ct * kCacheline);
      g_bgen_geno_bufs[tidx] = (uint16_t*)bigstack_alloc_raw(bgen_geno_cacheline_ct * kCacheline);
    }
    g_sample_ct = sample_ct;
    g_variant_include = variant_include;
    g_sample_include = sample_include;
    g_sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
    g_calc_thread_ct = calc_thread_ct;
    g_refalt1_select = refalt1_select;
    get_xymt_start_and_end(cip, kChrOffsetY, &g_y_start, &g_y_end);
    g_ref_allele_second = ref_allele_second;
    g_cip = cip;

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
    const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
    const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
    uint32_t parity = 0;
    uint32_t read_block_idx = 0;
    uint32_t write_variant_uidx = 0;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_slen = 0;

    uint32_t prev_block_write_ct = 0;
    uint32_t variant_idx = 0;
    uint32_t is_last_block = 0;
    uint32_t cur_read_block_size = read_block_size;
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    LOGPRINTFWW5("Writing %s ... ", outname);
    fputs("0%", stdout);
    fflush(stdout);
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    while (1) {
      uintptr_t cur_block_write_ct = 0;
      if (!is_last_block) {
        while (read_block_idx < read_block_ct_m1) {
          cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
          if (cur_block_write_ct) {
            break;
          }
          ++read_block_idx;
        }
        if (read_block_idx == read_block_ct_m1) {
          cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
          cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
        }
        if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
          if (variant_idx) {
            join_threads2z(calc_thread_ct, 0, threads);
            g_cur_block_write_ct = 0;
            error_cleanup_threads2z(export_bgen13_thread, calc_thread_ct, threads);
          }
          goto export_bgen13_ret_READ_FAIL;
        }
      }
      if (variant_idx) {
        join_threads2z(calc_thread_ct, is_last_block, threads);
        reterr = g_error_ret;
        if (reterr) {
          if (!is_last_block) {
            g_cur_block_write_ct = 0;
            error_cleanup_threads2z(export_bgen13_thread, calc_thread_ct, threads);
          }
          if (reterr == kPglRetMalformedInput) {
            logprint("\n");
            logerrprint("Error: Malformed .pgen file.\n");
          }
          goto export_bgen13_ret_1;
        }
      }
      if (!is_last_block) {
        g_cur_block_write_ct = cur_block_write_ct;
        compute_uidx_start_partition(variant_include, cur_block_write_ct, calc_thread_ct, read_block_idx * read_block_size, g_read_variant_uidx_starts);
        for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
          g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
          g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
        }
        is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
        if (spawn_threads2z(export_bgen13_thread, calc_thread_ct, is_last_block, threads)) {
          goto export_bgen13_ret_THREAD_CREATE_FAIL;
        }
      }
      parity = 1 - parity;
      if (variant_idx) {
        // write *previous* block results
        const unsigned char* compressed_data_iter = g_writebufs[parity];
        const uint32_t* variant_bytect_iter = g_variant_bytects[parity];
        for (uint32_t variant_bidx = 0; variant_bidx < prev_block_write_ct; ++variant_bidx, ++write_variant_uidx) {
          next_set_unsafe_ck(variant_include, &write_variant_uidx);
          if (write_variant_uidx >= chr_end) {
            do {
              ++chr_fo_idx;
              chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
            } while (write_variant_uidx >= chr_end);
            const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
            char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
            chr_slen = (uintptr_t)(chr_name_end - chr_buf);
          }
          const char* cur_variant_id = variant_ids[write_variant_uidx];
          const uint32_t id_slen = strlen(cur_variant_id);
          memcpy(&(writebuf[6]), &id_slen, 4);
          // deliberately clobber top two bytes
          unsigned char* writebuf_iter = (unsigned char*)memcpya(&(writebuf[8]), cur_variant_id, id_slen);
          memcpy(writebuf_iter, &chr_slen, 4);
          writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[2]), chr_buf, chr_slen);
          memcpy(writebuf_iter, &(variant_bps[write_variant_uidx]), 4);
          writebuf_iter = &(writebuf_iter[4]);
          uintptr_t variant_allele_idx_base = write_variant_uidx * 2;
          if (variant_allele_idxs) {
            variant_allele_idx_base = variant_allele_idxs[write_variant_uidx];
          }
          char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
          if (refalt1_select) {
            ref_allele_idx = refalt1_select[write_variant_uidx * 2];
            alt1_allele_idx = refalt1_select[write_variant_uidx * 2 + 1];
          }
          const char* first_allele;
          const char* second_allele;
          if (ref_allele_second) {
            first_allele = cur_alleles[alt1_allele_idx];
            second_allele = cur_alleles[ref_allele_idx];
          } else {
            first_allele = cur_alleles[ref_allele_idx];
            second_allele = cur_alleles[alt1_allele_idx];
          }
          uint32_t allele_slen = strlen(first_allele);
          memcpy(writebuf_iter, &allele_slen, 4);
          writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[4]), first_allele, allele_slen);
          allele_slen = strlen(second_allele);
          memcpy(writebuf_iter, &allele_slen, 4);
          writebuf_iter = (unsigned char*)memcpya(&(writebuf_iter[4]), second_allele, allele_slen);
          const uint32_t cur_variant_bytect = *variant_bytect_iter++;
          memcpy(writebuf_iter, &cur_variant_bytect, 4);
          writebuf_iter = &(writebuf_iter[4]);
          memcpy(writebuf_iter, compressed_data_iter, cur_variant_bytect);
          writebuf_iter = &(writebuf_iter[cur_variant_bytect]);
          compressed_data_iter = &(compressed_data_iter[bgen_compressed_buf_max]);
          if (fwrite_checked(writebuf, writebuf_iter - writebuf, outfile)) {
            if (variant_idx < variant_ct) {
              join_threads2z(calc_thread_ct, is_last_block, threads);
              if (!is_last_block) {
                g_cur_block_write_ct = 0;
                error_cleanup_threads2z(export_bgen13_thread, calc_thread_ct, threads);
              }
            }
            goto export_bgen13_ret_WRITE_FAIL;
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
        next_print_variant_idx = (pct * ((uint64_t)variant_ct)) / 100;
      }
      ++read_block_idx;
      prev_block_write_ct = cur_block_write_ct;
      variant_idx += cur_block_write_ct;
      pgfip->block_base = main_loadbufs[parity];
    }
    if (fclose_null(&outfile)) {
      goto export_bgen13_ret_WRITE_FAIL;
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
    const uint32_t sample_ctav2 = acc2_vec_ct * kQuatersPerVec;
    const uintptr_t acc32_offset = acc2_vec_ct * (7 * k1LU * kWordsPerVec);
    uint32_t* scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[0][acc32_offset]));
    for (uint32_t tidx = 1; tidx < calc_thread_ct; ++tidx) {
      const uint32_t* thread_scrambled_missing_cts = (uint32_t*)(&(g_missing_acc2[tidx][acc32_offset]));
      for (uint32_t uii = 0; uii < sample_ctav2; ++uii) {
        scrambled_missing_cts[uii] += thread_scrambled_missing_cts[uii];
      }
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      const uint32_t scrambled_idx = scramble_2_4_8_32(sample_idx);
      sample_missing_geno_cts[sample_idx] = scrambled_missing_cts[scrambled_idx];
    }
  }
  while (0) {
  export_bgen13_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_bgen13_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_bgen13_ret_READ_FAIL:
    reterr = kPglRetReadFail;
    break;
  export_bgen13_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
#ifdef __LP64__
  export_bgen13_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
#endif
  export_bgen13_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 export_bgen13_ret_1:
  ZWRAP_useZSTDcompression(1);
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}
*/

pglerr_t export_ox_sample(const char* outname, const uintptr_t* sample_include, const char* sample_ids, const uint32_t* sample_missing_geno_cts, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, uint32_t sample_ct, uintptr_t max_sample_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t variant_ct, uint32_t y_ct) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t pheno_ctl = BITCT_TO_WORDCT(pheno_ct);
    char* writebuf;
    uintptr_t* is_basic_categorical;
    // if phenotype is categorical, and all (non-null) category names are of
    // the form P[positive integer], then it's best to emit the positive
    // integer in the name string instead of the internal index.
    if (bigstack_calloc_ul(pheno_ctl, &is_basic_categorical) ||
        bigstack_alloc_c(kMaxMediumLine + max_sample_id_blen + 32 + pheno_ct * MAXV(kMaxMissingPhenostrBlen, 16), &writebuf)) {
      goto export_ox_sample_ret_NOMEM;
    }

    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_ox_sample_ret_OPEN_FAIL;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = strcpya(writebuf, "ID_1 ID_2 missing sex");
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      *write_iter++ = ' ';
      write_iter = strcpya(write_iter, &(pheno_names[pheno_idx * max_pheno_name_blen]));
      const pheno_col_t* cur_pheno_col = &(pheno_cols[pheno_idx]);
      if (cur_pheno_col->type_code == kPhenoDtypeCat) {
        const uint32_t nn_cat_ct = cur_pheno_col->nonnull_category_ct;
        char** cur_cat_names = cur_pheno_col->category_names;
        uint32_t cat_idx;
        for (cat_idx = 1; cat_idx <= nn_cat_ct; ++cat_idx) {
          const char* cat_name_iter = cur_cat_names[cat_idx];
          if (*cat_name_iter == 'C') {
            uint32_t char_code = *(++cat_name_iter);
            if ((char_code - 49) < 9) {
              uint32_t uii;
              if (!scan_posint_capped(cat_name_iter, 0x7fffffff, &uii)) {
                continue;
              }
            }
          }
          break;
        }
        if (cat_idx == nn_cat_ct + 1) {
          set_bit(pheno_idx, is_basic_categorical);
        }
      }
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto export_ox_sample_ret_WRITE_FAIL;
      }
    }
    append_binary_eoln(&write_iter);

    write_iter = strcpya(write_iter, "0 0 0 D");
    for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
      *write_iter++ = ' ';
      const pheno_dtype_t cur_type_code = pheno_cols[pheno_idx].type_code;
      if (cur_type_code == kPhenoDtypeCc) {
        *write_iter++ = 'B';
      } else if (cur_type_code == kPhenoDtypeQt) {
        // .psam file does not distinguish between "continuous covariate" and
        // "continuous phenotype", that's lost on round-trip
        *write_iter++ = 'P';
      } else {
        *write_iter++ = 'D';
      }
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto export_ox_sample_ret_WRITE_FAIL;
      }
    }
    append_binary_eoln(&write_iter);

    const double nonmale_geno_ct_recip = 1.0 / ((double)((int32_t)(variant_ct - y_ct)));
    const double male_geno_ct_recip = 1.0 / ((double)((int32_t)variant_ct));
    uintptr_t sample_uidx = 0;
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_ul_unsafe_ck(sample_include, &sample_uidx);
      const char* cur_sample_id = &(sample_ids[max_sample_id_blen * sample_uidx]);
      const char* fid_end = (const char*)rawmemchr(cur_sample_id, '\t');
      write_iter = memcpyax(write_iter, cur_sample_id, (uintptr_t)(fid_end - cur_sample_id), ' ');
      write_iter = strcpya(write_iter, &(fid_end[1]));
      *write_iter++ = ' ';
      const int32_t cur_missing_geno_ct = sample_missing_geno_cts[sample_idx];
      if (is_set(sex_male, sample_uidx)) {
        write_iter = dtoa_g(cur_missing_geno_ct * male_geno_ct_recip, write_iter);
        write_iter = strcpya(write_iter, " 1");
      } else {
        write_iter = dtoa_g(cur_missing_geno_ct * nonmale_geno_ct_recip, write_iter);
        *write_iter++ = ' ';
        if (is_set(sex_nm, sample_uidx)) {
          *write_iter++ = '2';
        } else {
          write_iter = strcpya(write_iter, "NA");
        }
      }
      for (uint32_t pheno_idx = 0; pheno_idx < pheno_ct; ++pheno_idx) {
        *write_iter++ = ' ';
        const pheno_col_t* cur_pheno_col = &(pheno_cols[pheno_idx]);
        if (!is_set(cur_pheno_col->nonmiss, sample_uidx)) {
          write_iter = strcpya(write_iter, "NA");
        } else {
          const pheno_dtype_t cur_type_code = cur_pheno_col->type_code;
          if (cur_type_code == kPhenoDtypeCc) {
            *write_iter++ = '0' + is_set(cur_pheno_col->data.cc, sample_uidx);
          } else if (cur_type_code == kPhenoDtypeQt) {
            write_iter = dtoa_g(cur_pheno_col->data.qt[sample_uidx], write_iter);
          } else {
            const uint32_t cur_cat_idx = cur_pheno_col->data.cat[sample_uidx];
            if (is_set(is_basic_categorical, pheno_idx)) {
              write_iter = strcpya(write_iter, &(cur_pheno_col->category_names[cur_cat_idx][1]));
            } else {
              write_iter = uint32toa(cur_cat_idx, write_iter);
            }
          }
        }
      }
      append_binary_eoln(&write_iter);
      if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
        goto export_ox_sample_ret_WRITE_FAIL;
      }
    }
    if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
      goto export_ox_sample_ret_WRITE_FAIL;
    }
  }
  while (0) {
  export_ox_sample_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_ox_sample_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_ox_sample_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
  fclose_cond(outfile);
  bigstack_reset(bigstack_mark);
  return reterr;
}

uint32_t valid_vcf_allele_code(const char* allele_code) {
  // returns 1 if probably valid (angle-bracket case is not exhaustively
  // checked), 0 if definitely not
  uint32_t uii = (unsigned char)(*allele_code);
  if ((uii == '<') || ((uii == '*') && (!allele_code[1]))) {
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
    uii = (unsigned char)(*(++allele_code));
  } while (uii);
  return 1;
}

char* diploid_vcf_dosage_print(uint32_t dosage_int, uint32_t write_ds, char* write_iter) {
  if (write_ds) {
    return print_small_dosage(dosage_int, write_iter);
  }
  if (dosage_int <= kDosageMid) {
    write_iter = print_small_dosage(kDosageMid - dosage_int, write_iter);
    *write_iter++ = ',';
    write_iter = print_small_dosage(dosage_int, write_iter);
    return strcpya(write_iter, ",0");
  }
  write_iter = strcpya(write_iter, "0,");
  write_iter = print_small_dosage(kDosageMax - dosage_int, write_iter);
  *write_iter++ = ',';
  return print_small_dosage(dosage_int - kDosageMid, write_iter);
}

// assumes rawval is in [0, 327679]
static_assert(kDosageMax == 32768, "haploid_dosage_print() needs to be updated.");
char* haploid_dosage_print(uint32_t rawval, char* start) {
  // Instead of constant 5-digit precision, we print fewer digits whenever that
  // doesn't interfere with proper round-tripping.  I.e. we search for the
  // shortest string in
  //   ((n - 0.5)/32768, (n + 0.5)/32768).
  *start++ = '0' + (rawval / 32768);
  rawval = rawval % 32768;
  if (!rawval) {
    // this shouldn't come up for now?
    return start;
  }
  *start++ = '.';

  // (rawval * 2) is in 65536ths
  // 65536 * 625 = 40960k
  const uint32_t range_top_40960k = rawval * 1250 + 625;
  // ok to check half-open interval since we never hit boundary
  if ((range_top_40960k % 4096) < 1250) {
    // when this is true, the four-decimal-place approximation is in the range
    // which round-trips back to our original number.
    const uint32_t four_decimal_places = range_top_40960k / 4096;
    return uitoa_trunc4(four_decimal_places, start);
  }

  // we wish to print (100000 * remainder + 16384) / 32768, left-0-padded.  and
  // may as well banker's round too.
  //
  // banker's rounding yields a different result than regular rounding for n/64
  // when n is congruent to 1 mod 4.  32768/64 = 512.
  const uint32_t five_decimal_places = ((3125 * rawval + 512) / 1024) - ((rawval % 2048) == 512);
  const uint32_t first_decimal_place = five_decimal_places / 10000;
  *start++ = '0' + first_decimal_place;
  const uint32_t last_four_digits = five_decimal_places - first_decimal_place * 10000;
  if (last_four_digits) {
    return uitoa_trunc4(last_four_digits, start);
  }
  return start;
}

boolerr_t flexbwrite_flush(char* buf, size_t len, FILE* outfile, BGZF* bgz_outfile) {
  if (outfile) {
    return fwrite_checked(buf, len, outfile);
  }
  return (bgzf_write(bgz_outfile, buf, len) < 0);
}


// these assume buf_flush - buf = kMaxMediumLine
// outfile should be nullptr iff we're doing bgzf compression
boolerr_t flexbwrite_flush2(char* buf_flush, FILE* outfile, BGZF* bgz_outfile, char** write_iter_ptr) {
  char* buf = &(buf_flush[-((int32_t)kMaxMediumLine)]);
  char* buf_end = *write_iter_ptr;
  *write_iter_ptr = buf;
  return flexbwrite_flush(buf, (uintptr_t)(buf_end - buf), outfile, bgz_outfile);
}

static inline boolerr_t flexbwrite_ck(char* buf_flush, FILE* outfile, BGZF* bgz_outfile, char** write_iter_ptr) {
  if ((*write_iter_ptr) < buf_flush) {
    return 0;
  }
  return flexbwrite_flush2(buf_flush, outfile, bgz_outfile, write_iter_ptr);
}


#ifdef __arm__
  #error "Unaligned accesses in export_vcf()."
#endif
pglerr_t export_vcf(char* xheader, const uintptr_t* sample_include, const uint32_t* sample_include_cumulative_popcounts, const char* sample_ids, const char* sids, const uintptr_t* sex_male_collapsed, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, __maybe_unused uint32_t max_thread_ct, exportf_flags_t exportf_modifier, idpaste_t exportf_id_paste, char exportf_id_delim, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  gzFile gz_pvar_reload = nullptr;
  BGZF* bgz_outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    if (!(exportf_modifier & kfExportfBgz)) {
      strcpy(outname_end, ".vcf");
      if (fopen_checked(outname, FOPEN_WB, &outfile)) {
        goto export_vcf_ret_OPEN_FAIL;
      }
    } else {
      strcpy(outname_end, ".vcf.gz");
      bgz_outfile = bgzf_open(outname, "w");
      if (!bgz_outfile) {
        goto export_vcf_ret_OPEN_FAIL;
      }
#ifndef _WIN32
      if (max_thread_ct > 1) {
        // 128 doesn't seem any worse than 256 (and is clearly better than 64)
        // also tried reducing thread count by 1, that seems worse
        if (bgzf_mt2(g_bigstack_end, MINV(128, max_thread_ct), 128, &g_bigstack_base, bgz_outfile)) {
          goto export_vcf_ret_NOMEM;
        }
      }
#endif
    }
    const uint32_t max_chr_blen = get_max_chr_slen(cip) + 1;
    // CHROM, POS, ID, REF, one ALT, eoln
    uintptr_t writebuf_blen = kMaxIdSlen + 32 + max_chr_blen + 2 * max_allele_slen;
    // QUAL, FILTER, INFO, FORMAT, genotypes, eoln
    // needs to be larger for >9 alt alleles
    const uint32_t dosage_force = (exportf_modifier / kfExportfVcfDosageForce) & 1;
    uint32_t write_ds = (exportf_modifier / kfExportfVcfDosageDs) & 1;
    uint32_t write_gp_or_ds = write_ds || (exportf_modifier & kfExportfVcfDosageGp);
    if ((!dosage_force) && write_gp_or_ds && (!(pgfip->gflags & kfPgenGlobalDosagePresent))) {
      write_gp_or_ds = 0;
      LOGERRPRINTF("Warning: No dosage data present.  %s field will not be exported.\n", write_ds? "DS" : "GP");
      write_ds = 0;
    }
    if (writebuf_blen < ((4 * k1LU) + write_gp_or_ds * 24 - write_ds * 16) * sample_ct + 32 + max_filter_slen + info_reload_slen) {
      writebuf_blen = ((4 * k1LU) + write_gp_or_ds * 24 - write_ds * 16) * sample_ct + 32 + max_filter_slen + info_reload_slen;
    }
    writebuf_blen += kMaxMediumLine;
    char* writebuf;
    if (bigstack_alloc_c(writebuf_blen, &writebuf)) {
      goto export_vcf_ret_NOMEM;
    }
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    char* write_iter = strcpya(writebuf, "##fileformat=VCFv4.3" EOLN_STR "##fileDate=");
    time_t rawtime;
    time(&rawtime);
    struct tm* loctime;
    loctime = localtime(&rawtime);
    write_iter += strftime(write_iter, kMaxMediumLine, "%Y%m%d", loctime);
    write_iter = strcpya(write_iter, EOLN_STR "##source=PLINKv2.00" EOLN_STR);
    if (cip->chrset_source) {
      append_chrset_line(cip, &write_iter);
    }
    if (flexbwrite_flush(writebuf, write_iter - writebuf, outfile, bgz_outfile)) {
      goto export_vcf_ret_WRITE_FAIL;
    }
    const uint32_t chr_ctl = BITCT_TO_WORDCT(cip->chr_ct);
    uintptr_t* written_contig_header_lines;
    if (bigstack_calloc_ul(chr_ctl, &written_contig_header_lines)) {
      goto export_vcf_ret_NOMEM;
    }
    if (xheader) {
      memcpy(writebuf, "##contig=<ID=", 13);
      char* xheader_iter = xheader;
      char* xheader_end = &(xheader[xheader_blen]);
      char* line_end = xheader;
      while (line_end != xheader_end) {
        xheader_iter = line_end;
        line_end = (char*)rawmemchr(xheader_iter, '\n');
        ++line_end;
        const uint32_t slen = (uintptr_t)(line_end - xheader_iter);
        if ((slen > 14) && str_startswith2(xheader_iter, "##contig=<ID=")) {
          char* contig_name_start = &(xheader_iter[13]);
          char* contig_name_end = (char*)memchr(contig_name_start, ',', slen - 14);
          if (!contig_name_end) {
            // if this line is technically well-formed (ends in '>'), it's
            // useless anyway, throw it out
            continue;
          }
          const int32_t chr_idx = get_chr_code_counted(cip, contig_name_end - contig_name_start, contig_name_start);
          if (chr_idx < 0) {
            continue;
          }
          const uint32_t chr_fo_idx = cip->chr_idx_to_foidx[(uint32_t)chr_idx];
          if (IS_SET(written_contig_header_lines, chr_fo_idx)) {
            logerrprint("Error: Duplicate ##contig line in .pvar file.\n");
            goto export_vcf_ret_MALFORMED_INPUT;
          }
          SET_BIT(chr_fo_idx, written_contig_header_lines);
          // if --output-chr was used at some point, we need to sync the
          // ##contig chromosome code with the code in the VCF body.
          write_iter = chr_name_write(cip, chr_idx, &(writebuf[13]));
          if (flexbwrite_flush(writebuf, write_iter - writebuf, outfile, bgz_outfile)) {
            goto export_vcf_ret_WRITE_FAIL;
          }
          if (flexbwrite_flush(contig_name_end, (uintptr_t)(line_end - contig_name_end), outfile, bgz_outfile)) {
            goto export_vcf_ret_WRITE_FAIL;
          }
        } else {
          if (flexbwrite_flush(xheader_iter, slen, outfile, bgz_outfile)) {
            goto export_vcf_ret_WRITE_FAIL;
          }
        }
      }
    }
    write_iter = writebuf;
    // fill in the missing ##contig lines
    uint32_t contig_zero_written = 0;
    for (uint32_t chr_fo_idx = 0; chr_fo_idx < cip->chr_ct; ++chr_fo_idx) {
      if (IS_SET(written_contig_header_lines, chr_fo_idx)) {
        continue;
      }
      const int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if ((!IS_SET(cip->chr_mask, chr_idx)) || are_all_bits_zero(variant_include, cip->chr_fo_vidx_start[chr_fo_idx], cip->chr_fo_vidx_start[chr_fo_idx + 1])) {
        continue;
      }
      char* chr_name_write_start = strcpya(write_iter, "##contig=<ID=");
      char* chr_name_write_end = chr_name_write(cip, chr_idx, chr_name_write_start);
      if ((*chr_name_write_start == '0') && (chr_name_write_end == &(chr_name_write_start[1]))) {
        // --allow-extra-chr 0 special case
        if (contig_zero_written) {
          continue;
        }
        contig_zero_written = 1;
        write_iter = strcpya(chr_name_write_end, ",length=2147483645");
      } else {
        if (memchr(chr_name_write_start, ':', chr_name_write_end - chr_name_write_start)) {
          logerrprint("Error: VCF chromosome codes may not include the ':' character.\n");
          goto export_vcf_ret_MALFORMED_INPUT;
        }
        write_iter = strcpya(chr_name_write_end, ",length=");
        if (1) {
          write_iter = uint32toa(variant_bps[cip->chr_fo_vidx_start[chr_fo_idx + 1] - 1] + 1, write_iter);
        } else {
          // todo: unsorted map case
        }
      }
      *write_iter++ = '>';
      append_binary_eoln(&write_iter);
      if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
        goto export_vcf_ret_WRITE_FAIL;
      }
    }
    bigstack_reset(written_contig_header_lines);
    const uint32_t all_nonref = pgfip->gflags & kfPgenGlobalAllNonref;
    const uintptr_t* nonref_flags = pgfip->nonref_flags;
    const uint32_t raw_variant_ctl = BITCT_TO_WORDCT(raw_variant_ct);
    uint32_t write_pr = all_nonref;
    if (nonref_flags) {
      for (uint32_t widx = 0; widx < raw_variant_ctl; ++widx) {
        if (variant_include[widx] & nonref_flags[widx]) {
          write_pr = 1;
          break;
        }
      }
    }
    if (write_pr) {
      if (xheader_info_pr_nonflag) {
        logprint("\n");
        logerrprint("Error: Conflicting INFO:PR fields.  Either fix all REF alleles so that the\n'provisional reference' field is no longer needed, or remove/rename the other\nINFO:PR field.\n");
        goto export_vcf_ret_INCONSISTENT_INPUT;
      }
      if (!xheader_info_pr) {
        write_iter = strcpya(write_iter, "##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\">" EOLN_STR);
      }
    }
    if (write_ds) {
      write_iter = strcpya(write_iter, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]\">" EOLN_STR);
    } else if (write_gp_or_ds) {
      write_iter = strcpya(write_iter, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Phred-scaled Genotype Likelihoods\">" EOLN_STR);
    }
    // possible todo: optionally export .psam information as
    // PEDIGREE/META/SAMPLE lines in header, and make --vcf be able to read it
    write_iter = strcpya(write_iter, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" EOLN_STR "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    uint32_t write_sid = 0;
    // possible for both MAYBESID and SID to be set
    if (exportf_id_paste & kfIdpasteSid) {
      write_sid = 1;
      if (!sids) {
        max_sid_blen = 2;
      }
    } else if ((exportf_id_paste & kfIdpasteMaybesid) && sids) {
      // no nonzero check in load_psam(), so we have to do it here
      uint32_t sample_uidx = 0;
      for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
        next_set_unsafe_ck(sample_include, &sample_uidx);
        if (!strequal_k2(&(sids[sample_uidx * max_sid_blen]), "0")) {
          write_sid = 1;
          break;
        }
      }
    }
    uint32_t sample_uidx = 0;
    uint32_t id_delim_warning = 0;
    char id_delim = exportf_id_delim? exportf_id_delim : '_';
    const uintptr_t max_exported_sample_id_blen = max_sample_id_blen + write_sid * max_sid_blen;
    char* exported_sample_ids;
    const uint32_t exported_id_htable_size = get_htable_min_size(sample_ct);
    uint32_t* exported_id_htable;
    // check for duplicates
    if (bigstack_alloc_c(sample_ct * max_exported_sample_id_blen, &exported_sample_ids) ||
        bigstack_alloc_ui(exported_id_htable_size, &exported_id_htable)) {
      goto export_vcf_ret_NOMEM;
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx, ++sample_uidx) {
      next_set_unsafe_ck(sample_include, &sample_uidx);
      const char* orig_sample_id = &(sample_ids[sample_uidx * max_sample_id_blen]);
      const char* orig_fid_end = (const char*)rawmemchr(orig_sample_id, '\t');
      char* exported_sample_ids_iter = &(exported_sample_ids[sample_idx * max_exported_sample_id_blen]);
      if (exportf_id_paste & kfIdpasteFid) {
        const uint32_t fid_slen = (uintptr_t)(orig_fid_end - orig_sample_id);
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
    if (id_delim_warning) {
      if (exportf_id_delim) {
        LOGERRPRINTF("Warning: '%c' present in original sample IDs; --vcf will not be able to\nreconstruct them.  Consider rerunning with a different --export id-delim=\nvalue.\n", exportf_id_delim);
      } else {
        logerrprint("Warning: '_' present in original sample IDs; --vcf will not be able to\nreconstruct them.  Consider rerunning with a suitable --export id-delim= value.\n");
      }
    }
    if (populate_strbox_htable(exported_sample_ids, sample_ct, max_exported_sample_id_blen, exported_id_htable_size, exported_id_htable)) {
      logerrprint("Warning: Duplicate sample ID(s) present in exported VCF file.\n");
    }
    for (uint32_t sample_idx = 0; sample_idx < sample_ct; ++sample_idx) {
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, &(exported_sample_ids[sample_idx * max_exported_sample_id_blen]));
      if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
        goto export_vcf_ret_WRITE_FAIL;
      }
    }
    append_binary_eoln(&write_iter);
    bigstack_reset(exported_sample_ids);

    LOGPRINTFWW5("--export vcf%s to %s ... ", bgz_outfile? " bgz" : "", outname);
    fputs("0%", stdout);
    fflush(stdout);

    // includes trailing tab
    char* chr_buf;

    const uint32_t sample_ctl2 = QUATERCT_TO_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uintptr_t* genovec;
    uintptr_t* allele_include;
    if (bigstack_alloc_c(max_chr_blen, &chr_buf) ||
    // if we weren't using bigstack_alloc, this would need to be sample_ctaw2
        bigstack_alloc_ul(sample_ctl2, &genovec) ||
        bigstack_alloc_ul(BITCT_TO_WORDCT(kPglMaxAltAlleleCt), &allele_include)) {
      goto export_vcf_ret_NOMEM;
    }
    // For now, if phased data is present, each homozygous call is represented
    // as phased iff the previous heterozygous call was phased.  (If no
    // previous heterozygous call exists, it's treated as phased.)  This does
    // the right thing when the entire genome is phased, and it induces about
    // as good a phase set approximation as you can get without explicitly
    // saving that info.  But that approximation is still pretty inaccurate; as
    // soon as we have any use for them, explicit phase set support should be
    // added to pgenlib.
    const uint32_t some_phased = (pgfip->gflags / kfPgenGlobalHardcallPhasePresent) & 1;
    uintptr_t* prev_phased = nullptr;
    uintptr_t* phasepresent = nullptr;
    uintptr_t* phaseinfo = nullptr;
    if (some_phased) {
      if (bigstack_alloc_ul(sample_ctl, &prev_phased) ||
          bigstack_alloc_ul(sample_ctl, &phasepresent) ||
          bigstack_alloc_ul(sample_ctl, &phaseinfo)) {
        goto export_vcf_ret_NOMEM;
      }
      fill_all_bits(sample_ct, prev_phased);
    }

    uintptr_t* dosage_present = nullptr;
    dosage_t* dosage_vals = nullptr;
    if (write_gp_or_ds) {
      if (bigstack_alloc_ul(sample_ctl, &dosage_present) ||
          bigstack_alloc_dosage(sample_ct, &dosage_vals)) {
        goto export_vcf_ret_NOMEM;
      }
    }

    char* loadbuf = nullptr;
    uintptr_t loadbuf_size = 0;
    uint32_t info_col_idx = 0;
    if (pvar_info_reload) {
      reterr = pvar_info_open_and_reload_header(pvar_info_reload, &gz_pvar_reload, &loadbuf, &loadbuf_size, &info_col_idx);
      if (reterr) {
        goto export_vcf_ret_1;
      }
    }

    // assumes little-endian
    uint32_t basic_genotext[4];
    basic_genotext[0] = 0x302f3009; // \t0/0
    basic_genotext[1] = 0x312f3009; // \t0/1
    basic_genotext[2] = 0x312f3109; // \t1/1
    basic_genotext[3] = 0x2e2f2e09; // \t./.
    char haploid_genotext[4][4];
    uint32_t haploid_genotext_blen[8]; // 4..7 = male chrX
    memcpy(haploid_genotext[0], "\t0/0", 4);
    memcpy(haploid_genotext[1], "\t0/1", 4);
    memcpy(haploid_genotext[2], "\t1/1", 4);
    memcpy(haploid_genotext[3], "\t./.", 4);
    haploid_genotext_blen[1] = 4;
    haploid_genotext_blen[4] = 2;
    haploid_genotext_blen[5] = 4;
    haploid_genotext_blen[6] = 2;
    haploid_genotext_blen[7] = 2;
    // don't bother exporting GP for hardcalls
    // usually don't bother for DS, but DS-force is an exception
    char dosage_inttext[16]; // 4..7 = haploid, 5 should never be looked up
    memcpy(dosage_inttext, ":0:1:2:.:0:.:1:.", 16);
    const char* dot_ptr = &(g_one_char_strs[92]);
    const uint32_t sample_ctl2_m1 = sample_ctl2 - 1;
    uint32_t chr_fo_idx = UINT32_MAX;
    uint32_t chr_end = 0;
    uint32_t chr_buf_blen = 0;
    uint32_t variant_uidx = 0;
    uint32_t is_x = 0;
    uint32_t is_haploid_or_mt = 0; // includes chrX and chrY
    uint32_t pct = 0;
    uint32_t next_print_variant_idx = variant_ct / 100;
    uint32_t gz_variant_uidx = 0;
    uint32_t ref_allele_idx = 0;
    uint32_t alt1_allele_idx = 1;
    uint32_t cur_allele_ct = 2;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      // a lot of this is redundant with write_pvar(), may want to factor the
      // commonalities out
      next_set_unsafe_ck(variant_include, &variant_uidx);
      if (variant_uidx >= chr_end) {
        do {
          ++chr_fo_idx;
          chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
        } while (variant_uidx >= chr_end);
        int32_t chr_idx = cip->chr_file_order[chr_fo_idx];
        is_x = (chr_idx == cip->xymt_codes[kChrOffsetX]);
        is_haploid_or_mt = is_set(cip->haploid_mask, chr_idx) || (chr_idx == cip->xymt_codes[kChrOffsetMT]);
        // forced --merge-par, with diploid male output (is_x NOT set, but
        // chromosome code is X/chrX)
        if ((chr_idx == cip->xymt_codes[kChrOffsetPAR1]) || (chr_idx == cip->xymt_codes[kChrOffsetPAR2])) {
          chr_idx = cip->xymt_codes[kChrOffsetX];
        }
        char* chr_name_end = chr_name_write(cip, chr_idx, chr_buf);
        *chr_name_end = '\t';
        chr_buf_blen = 1 + (uintptr_t)(chr_name_end - chr_buf);
        if (is_haploid_or_mt) {
          if (is_x) {
            haploid_genotext_blen[0] = 4;
            haploid_genotext_blen[2] = 4;
            haploid_genotext_blen[3] = 4;
          } else {
            haploid_genotext_blen[0] = 2;
            haploid_genotext_blen[2] = 2;
            haploid_genotext_blen[3] = 2;
          }
        }
      }
      // #CHROM
      write_iter = memcpya(write_iter, chr_buf, chr_buf_blen);

      // POS
      write_iter = uint32toa_x(variant_bps[variant_uidx], '\t', write_iter);

      // ID
      write_iter = strcpyax(write_iter, variant_ids[variant_uidx], '\t');

      // REF, ALT
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
        variant_allele_idx_base = variant_allele_idxs[variant_uidx];
        cur_allele_ct = variant_allele_idxs[variant_uidx + 1] - variant_allele_idx_base;
      }
      char** cur_alleles = &(allele_storage[variant_allele_idx_base]);
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[variant_uidx * 2];
        alt1_allele_idx = refalt1_select[variant_uidx * 2 + 1];
        // this logic only works in the biallelic case
        assert(cur_allele_ct == 2);
        if (!is_haploid_or_mt) {
          if (alt1_allele_idx) {
            basic_genotext[0] = 0x302f3009;
            basic_genotext[2] = 0x312f3109;
          } else {
            basic_genotext[0] = 0x312f3109;
            basic_genotext[2] = 0x302f3009;
          }
        } else {
          if (alt1_allele_idx) {
            memcpy(haploid_genotext[0], "\t0/0", 4);
            memcpy(haploid_genotext[2], "\t1/1", 4);
          } else {
            memcpy(haploid_genotext[0], "\t1/1", 4);
            memcpy(haploid_genotext[2], "\t0/0", 4);
          }
        }
        if (alt1_allele_idx) {
          memcpy(dosage_inttext, ":0:1:2:.:0:.:1:.", 16);
        } else {
          memcpy(dosage_inttext, ":2:1:0:.:1:.:0:.", 16);
        }
      }
      if (cur_alleles[ref_allele_idx] != dot_ptr) {
        write_iter = strcpya(write_iter, cur_alleles[ref_allele_idx]);
      } else {
        *write_iter++ = 'N';
      }
      *write_iter++ = '\t';
      write_iter = strcpya(write_iter, cur_alleles[alt1_allele_idx]);
      if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
        goto export_vcf_ret_WRITE_FAIL;
      }
      if (cur_allele_ct > 2) {
        fill_all_bits(cur_allele_ct, allele_include);
        CLEAR_BIT(ref_allele_idx, allele_include);
        CLEAR_BIT(alt1_allele_idx, allele_include);
        uint32_t cur_allele_uidx = 0;
        uint32_t alt_allele_idx = 2;
        do {
          *write_iter++ = ',';
          next_set_unsafe_ck(allele_include, &cur_allele_uidx);
          write_iter = strcpya(write_iter, cur_alleles[cur_allele_uidx++]);
          if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
            goto export_vcf_ret_WRITE_FAIL;
          }
        } while (++alt_allele_idx < cur_allele_ct);
      }

      // QUAL
      *write_iter++ = '\t';
      if ((!pvar_qual_present) || (!IS_SET(pvar_qual_present, variant_uidx))) {
        *write_iter++ = '.';
      } else {
        write_iter = ftoa_g(pvar_quals[variant_uidx], write_iter);
      }

      // FILTER
      *write_iter++ = '\t';
      if ((!pvar_filter_present) || (!IS_SET(pvar_filter_present, variant_uidx))) {
        *write_iter++ = '.';
      } else if (!IS_SET(pvar_filter_npass, variant_uidx)) {
        write_iter = strcpya(write_iter, "PASS");
      } else {
        write_iter = strcpya(write_iter, pvar_filter_storage[variant_uidx]);
      }

      // INFO
      *write_iter++ = '\t';
      const uint32_t is_pr = all_nonref || (nonref_flags && IS_SET(nonref_flags, variant_uidx));
      if (gz_pvar_reload) {
        reterr = pvar_info_reload_and_write(loadbuf_size, xheader_info_pr, info_col_idx, variant_uidx, is_pr, gz_pvar_reload, &write_iter, &gz_variant_uidx, loadbuf);
        if (reterr) {
          goto export_vcf_ret_1;
        }
      } else {
        if (is_pr) {
          write_iter = strcpya(write_iter, "PR");
        } else {
          *write_iter++ = '.';
        }
      }

      // FORMAT
      write_iter = memcpyl3a(write_iter, "\tGT");

      uint32_t dosage_ct = 0;
      uint32_t is_explicit_alt1 = 0;
      uint32_t inner_loop_last = kBitsPerWordD2 - 1;
      uint32_t widx = 0;
      if (!some_phased) {
        // biallelic, nothing phased in entire file
        if (!write_gp_or_ds) {
          reterr = pgr_read_refalt1_genovec_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec);
        } else {
          reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
        }
        if (reterr) {
          goto export_vcf_ret_PGR_FAIL;
        }
        if ((!dosage_ct) && (!dosage_force)) {
          if (!is_haploid_or_mt) {
            // always 4 bytes wide, exploit that
            uint32_t* write_iter_ui_alias = (uint32_t*)write_iter;
            while (1) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = genovec[widx];
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                *write_iter_ui_alias++ = basic_genotext[genovec_word & 3];
                genovec_word >>= 2;
              }
              ++widx;
            }
            write_iter = (char*)write_iter_ui_alias;
          } else {
            // chrX: male homozygous/missing calls use only one character + tab
            // other haploid/MT: this is true for nonmales too
            while (1) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = genovec[widx];
              uint32_t sex_male_hw = is_x * (((const halfword_t*)sex_male_collapsed)[widx]);
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                const uint32_t cur_geno = genovec_word & 3;
                const uint32_t cur_is_male = sex_male_hw & 1;
                write_iter = memcpya(write_iter, haploid_genotext[cur_geno], haploid_genotext_blen[cur_geno + cur_is_male * 4]);
                genovec_word >>= 2;
                sex_male_hw >>= 1;
              }
              ++widx;
            }
          }
        } else {
          // some dosages present, or DS-force
          if (write_ds) {
            write_iter = memcpyl3a(write_iter, ":DS");
            if (!dosage_ct) {
              // DS-force, need to clear this
              fill_ulong_zero(sample_ctl, dosage_present);
            }
          } else {
            write_iter = memcpyl3a(write_iter, ":GP");
          }
          if (!alt1_allele_idx) {
            biallelic_dosage16_invert(dosage_ct, dosage_vals);
          }
          dosage_t* dosage_vals_iter = dosage_vals;
          if (!is_haploid_or_mt) {
            while (1) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = genovec[widx];
              uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                const uint32_t cur_geno = genovec_word & 3;
                write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
                if (dosage_present_hw & 1) {
                  *write_iter++ = ':';
                  const uint32_t dosage_int = *dosage_vals_iter++;
                  write_iter = diploid_vcf_dosage_print(dosage_int, write_ds, write_iter);
                } else if (dosage_force) {
                  write_iter = memcpya(write_iter, &(dosage_inttext[cur_geno * 2]), 2);
                }
                genovec_word >>= 2;
                dosage_present_hw >>= 1;
              }
              ++widx;
            }
          } else {
            while (1) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = genovec[widx];
              uint32_t sex_male_hw = is_x * (((const halfword_t*)sex_male_collapsed)[widx]);
              uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                const uint32_t cur_geno = genovec_word & 3;
                const uint32_t cur_is_male = sex_male_hw & 1;
                const uint32_t cur_genotext_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                write_iter = memcpya(write_iter, haploid_genotext[cur_geno], cur_genotext_blen);
                if (dosage_present_hw & 1) {
                  *write_iter++ = ':';
                  uint32_t dosage_int = *dosage_vals_iter++;
                  if (cur_genotext_blen == 2) {
                    if (write_ds) {
                      write_iter = haploid_dosage_print(dosage_int, write_iter);
                    } else {
                      write_iter = haploid_dosage_print(kDosageMax - dosage_int, write_iter);
                      *write_iter++ = ',';
                      write_iter = haploid_dosage_print(dosage_int, write_iter);
                    }
                  } else {
                    // het haploid, or female X
                    write_iter = diploid_vcf_dosage_print(dosage_int, write_ds, write_iter);
                  }
                } else if (dosage_force) {
                  write_iter = memcpya(write_iter, &(dosage_inttext[2 * cur_geno + 16 - 4 * cur_genotext_blen]), 2);
                }
                genovec_word >>= 2;
                sex_male_hw >>= 1;
                dosage_present_hw >>= 1;
              }
              ++widx;
            }
          }
        }
      } else {
        // biallelic, phased
        uint32_t at_least_one_phase_present;
        if (!write_gp_or_ds) {
          reterr = pgr_read_refalt1_genovec_hphase_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, phasepresent, phaseinfo, &at_least_one_phase_present);
        } else {
          reterr = pgr_read_refalt1_genovec_hphase_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, simple_pgrp, genovec, phasepresent, phaseinfo, &at_least_one_phase_present, dosage_present, dosage_vals, &dosage_ct, &is_explicit_alt1);
        }
        if (reterr) {
          goto export_vcf_ret_PGR_FAIL;
        }
        at_least_one_phase_present = (at_least_one_phase_present != 0);
        if ((!dosage_ct) && (!dosage_force)) {
          if (!is_haploid_or_mt) {
            uint32_t* write_iter_ui_alias = (uint32_t*)write_iter;
            while (1) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = genovec[widx];
              uint32_t prev_phased_halfword = ((halfword_t*)prev_phased)[widx];

              // zero this out if phasepresent_ct == 0
              const uint32_t phasepresent_halfword = at_least_one_phase_present * (((halfword_t*)phasepresent)[widx]);

              const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                const uintptr_t cur_geno = genovec_word & 3;

                // usually "\t0/0", etc.
                uint32_t cur_basic_genotext = basic_genotext[cur_geno];
                if (cur_geno == 1) {
                  const uint32_t cur_shift = (1U << sample_idx_lowbits);
                  if (phasepresent_halfword & cur_shift) {
                    prev_phased_halfword |= cur_shift;
                    if (phaseinfo_halfword & cur_shift) {
                      cur_basic_genotext ^= 0x1000100; // 0|1 -> 1|0
                    }
                  } else {
                    prev_phased_halfword &= ~cur_shift;
                  }
                }
                // '/' = ascii 47, '|' = ascii 124
                *write_iter_ui_alias++ = cur_basic_genotext + 0x4d0000 * ((prev_phased_halfword >> sample_idx_lowbits) & 1);
                genovec_word >>= 2;
              }
              ((halfword_t*)prev_phased)[widx] = prev_phased_halfword;
              ++widx;
            }
            write_iter = (char*)write_iter_ui_alias;
          } else {
            while (1) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = genovec[widx];
              uint32_t is_male_hw = is_x * (((const halfword_t*)sex_male_collapsed)[widx]);
              uint32_t prev_phased_halfword = ((halfword_t*)prev_phased)[widx];

              // zero this out if phasepresent_ct == 0
              const uint32_t phasepresent_halfword = at_least_one_phase_present * (((halfword_t*)phasepresent)[widx]);

              const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                const uint32_t cur_geno = genovec_word & 3;
                const uint32_t cur_is_male = is_male_hw & 1;
                const uint32_t cur_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                write_iter = memcpya(write_iter, haploid_genotext[cur_geno], cur_blen);
                if (cur_blen == 4) {
                  if (cur_geno == 1) {
                    // a bit redundant with how is_male_hw is handled, but
                    // updating this on every loop iteration doesn't seem better
                    const uint32_t cur_shift = (1U << sample_idx_lowbits);
                    if (phasepresent_halfword & cur_shift) {
                      prev_phased_halfword |= cur_shift;
                      if (phaseinfo_halfword & cur_shift) {
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
              ((halfword_t*)prev_phased)[widx] = prev_phased_halfword;
              ++widx;
            }
          }
        } else {
          // both dosage (or DS-force) and phase present
          if (write_ds) {
            write_iter = memcpyl3a(write_iter, ":DS");
            if (!dosage_ct) {
              fill_ulong_zero(sample_ctl, dosage_present);
            }
          } else {
            write_iter = memcpyl3a(write_iter, ":GP");
          }
          if (!alt1_allele_idx) {
            biallelic_dosage16_invert(dosage_ct, dosage_vals);
          }
          dosage_t* dosage_vals_iter = dosage_vals;
          if (!is_haploid_or_mt) {
            while (1) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = genovec[widx];
              uint32_t prev_phased_halfword = ((halfword_t*)prev_phased)[widx];

              // zero this out if phasepresent_ct == 0
              const uint32_t phasepresent_halfword = at_least_one_phase_present * (((halfword_t*)phasepresent)[widx]);

              const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
              const uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
              uint32_t cur_shift = 1;
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                const uint32_t cur_geno = genovec_word & 3;
                write_iter = memcpya(write_iter, &(basic_genotext[cur_geno]), 4);
                if (cur_geno == 1) {
                  if (phasepresent_halfword & cur_shift) {
                    prev_phased_halfword |= cur_shift;
                    if (phaseinfo_halfword & cur_shift) {
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
                  const uint32_t dosage_int = *dosage_vals_iter++;
                  write_iter = diploid_vcf_dosage_print(dosage_int, write_ds, write_iter);
                } else if (dosage_force) {
                  write_iter = memcpya(write_iter, &(dosage_inttext[cur_geno * 2]), 2);
                }
                genovec_word >>= 2;
                cur_shift <<= 1;
              }
              ((halfword_t*)prev_phased)[widx] = prev_phased_halfword;
              ++widx;
            }
          } else {
            while (1) {
              if (widx >= sample_ctl2_m1) {
                if (widx > sample_ctl2_m1) {
                  break;
                }
                inner_loop_last = (sample_ct - 1) % kBitsPerWordD2;
              }
              uintptr_t genovec_word = genovec[widx];
              uint32_t is_male_hw = is_x * (((const halfword_t*)sex_male_collapsed)[widx]);
              uint32_t prev_phased_halfword = ((halfword_t*)prev_phased)[widx];

              // zero this out if phasepresent_ct == 0
              const uint32_t phasepresent_halfword = at_least_one_phase_present * (((halfword_t*)phasepresent)[widx]);

              const uint32_t phaseinfo_halfword = ((halfword_t*)phaseinfo)[widx];
              const uint32_t dosage_present_hw = ((halfword_t*)dosage_present)[widx];
              uint32_t cur_shift = 1;
              for (uint32_t sample_idx_lowbits = 0; sample_idx_lowbits <= inner_loop_last; ++sample_idx_lowbits) {
                const uint32_t cur_geno = genovec_word & 3;
                const uint32_t cur_is_male = is_male_hw & 1;
                const uint32_t cur_blen = haploid_genotext_blen[cur_geno + cur_is_male * 4];
                write_iter = memcpya(write_iter, haploid_genotext[cur_geno], cur_blen);
                if (cur_blen == 4) {
                  if (cur_geno == 1) {
                    if (phasepresent_halfword & cur_shift) {
                      prev_phased_halfword |= cur_shift;
                      if (phaseinfo_halfword & cur_shift) {
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
                    const uint32_t dosage_int = *dosage_vals_iter++;
                    write_iter = diploid_vcf_dosage_print(dosage_int, write_ds, write_iter);
                  } else if (dosage_force) {
                    write_iter = memcpya(write_iter, &(dosage_inttext[cur_geno * 2]), 2);
                  }
                } else {
                  if (dosage_present_hw & cur_shift) {
                    *write_iter++ = ':';
                    const uint32_t dosage_int = *dosage_vals_iter++;
                    if (write_ds) {
                      write_iter = haploid_dosage_print(dosage_int, write_iter);
                    } else {
                      write_iter = haploid_dosage_print(kDosageMax - dosage_int, write_iter);
                      *write_iter++ = ',';
                      write_iter = haploid_dosage_print(dosage_int, write_iter);
                    }
                  } else if (dosage_force) {
                    write_iter = memcpya(write_iter, &(dosage_inttext[cur_geno * 2 + 8]), 2);
                  }
                }
                genovec_word >>= 2;
                is_male_hw >>= 1;
                cur_shift <<= 1;
              }
              ((halfword_t*)prev_phased)[widx] = prev_phased_halfword;
              ++widx;
            }
          }
        }
      }
      // todo: multiallelic cases (separate out cur_allele_ct <= 10)
      append_binary_eoln(&write_iter);
      if (flexbwrite_ck(writebuf_flush, outfile, bgz_outfile, &write_iter)) {
        goto export_vcf_ret_WRITE_FAIL;
      }
      if (variant_idx >= next_print_variant_idx) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_variant_idx = (pct * ((uint64_t)variant_ct)) / 100;
      }
    }
    if (write_iter != writebuf) {
      if (flexbwrite_flush(writebuf, write_iter - writebuf, outfile, bgz_outfile)) {
        goto export_vcf_ret_WRITE_FAIL;
      }
    }
    if (bgz_outfile) {
      if (bgzf_close(bgz_outfile)) {
        bgz_outfile = nullptr;
        goto export_vcf_ret_WRITE_FAIL;
      }
      bgz_outfile = nullptr;
    } else {
      if (fclose_null(&outfile)) {
        goto export_vcf_ret_WRITE_FAIL;
      }
    }
    if (pct > 10) {
      putc_unlocked('\b', stdout);
    }
    fputs("\b\b", stdout);
    LOGPRINTF("done.\n");
  }
  while (0) {
  export_vcf_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_vcf_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_vcf_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  export_vcf_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  export_vcf_ret_INCONSISTENT_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  export_vcf_ret_PGR_FAIL:
    if (reterr != kPglRetReadFail) {
      logprint("\n");
      logerrprint("Error: Malformed .pgen file.\n");
    }
  }
 export_vcf_ret_1:
  fclose_cond(outfile);
  gzclose_cond(gz_pvar_reload);
  if (bgz_outfile) {
    bgzf_close(bgz_outfile);
  }
  bigstack_reset(bigstack_mark);
  return reterr;
}

// more multithread globals
static dosage_t* g_smaj_dosagebuf = nullptr;
static uint32_t* g_write_vidx_starts = nullptr;
static const dosage_t kGenoToDosage[4] = {0, 16384, 32768, 65535};

THREAD_FUNC_DECL dosage_transpose_thread(void* arg) {
  const uintptr_t tidx = (uintptr_t)arg;
  const uint32_t sample_ct = g_sample_ct;
  const uint32_t sample_ctd4 = sample_ct / 4;
  const uint32_t sample_rem = sample_ct % 4;
  const uintptr_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(sample_ct);
  const uint32_t sample_ctab2 = kBytesPerWord * sample_ctaw2;
  const uintptr_t stride = g_stride;
  const uintptr_t* variant_include = g_variant_include;
  const alt_allele_ct_t* refalt1_select = g_refalt1_select;
  const uintptr_t* sample_include = g_sample_include;
  const uint32_t* sample_include_cumulative_popcounts = g_sample_include_cumulative_popcounts;
  pgen_reader_t* pgrp = g_pgr_ptrs[tidx];
  uintptr_t* genovec_buf = g_thread_write_genovecs[tidx];
  uintptr_t* dosagepresent_buf = g_thread_write_dosagepresents[tidx];
  dosage_t* dosagevals_buf = g_thread_write_dosagevals[tidx];
  while (1) {
    const uint32_t is_last_block = g_is_last_thread_block;
    const uint32_t cur_block_write_ct = g_cur_block_write_ct;
    const uint32_t vidx_end = g_write_vidx_starts[tidx + 1];
    uint32_t vidx_start = g_write_vidx_starts[tidx];
    if (cur_block_write_ct && (vidx_end != vidx_start)) {
      uint32_t variant_uidx = g_read_variant_uidx_starts[tidx];
      dosage_t* smaj_dosagebuf_iter = &(g_smaj_dosagebuf[vidx_start]);
      uint32_t dosage_cts[kDosagePerCacheline];
      do {
        uint32_t vidx_block_end = round_down_pow2(vidx_start, kDosagePerCacheline) + kDosagePerCacheline;
        if (vidx_block_end > vidx_end) {
          vidx_block_end = vidx_end;
        }
        const uint32_t vidx_block_size = vidx_block_end - vidx_start;

        // part 1: decompress data
        uintptr_t* genovec_iter = genovec_buf;
        uintptr_t* dosage_present_iter = dosagepresent_buf;
        dosage_t* dosage_vals_iter = dosagevals_buf;
        for (uint32_t vidx_offset = 0; vidx_offset < vidx_block_size; ++vidx_offset, ++variant_uidx) {
          next_set_unsafe_ck(variant_include, &variant_uidx);
          // todo: multiallelic case
          uint32_t dosage_ct;
          uint32_t is_explicit_alt1;
          const pglerr_t reterr = pgr_read_refalt1_genovec_dosage16_subset_unsafe(sample_include, sample_include_cumulative_popcounts, sample_ct, variant_uidx, pgrp, genovec_iter, dosage_present_iter, (uint16_t*)dosage_vals_iter, &dosage_ct, &is_explicit_alt1);
          if (reterr) {
            g_error_ret = reterr;
            goto dosage_transpose_thread_fail;
          }
          if ((!refalt1_select) || (!refalt1_select[variant_uidx * 2])) {
            genovec_invert_unsafe(sample_ct, genovec_iter);
            biallelic_dosage16_invert(dosage_ct, dosage_vals_iter);
          }
          genovec_iter = &(genovec_iter[sample_ctaw2]);
          dosage_present_iter = &(dosage_present_iter[sample_ctaw]);
          dosage_vals_iter = &(dosage_vals_iter[sample_ct]);
          dosage_cts[vidx_offset] = dosage_ct;
        }

        // part 2: process hardcalls for 4 samples at a time
        dosage_t* dosagebuf_write_iter0 = smaj_dosagebuf_iter;
        for (uint32_t sample4_idx = 0; sample4_idx < sample_ctd4; ++sample4_idx) {
          dosage_t* dosagebuf_write_iter1 = &(dosagebuf_write_iter0[stride]);
          dosage_t* dosagebuf_write_iter2 = &(dosagebuf_write_iter1[stride]);
          dosage_t* dosagebuf_write_iter3 = &(dosagebuf_write_iter2[stride]);
          const unsigned char* geno_read_iter = &(((const unsigned char*)genovec_buf)[sample4_idx]);
          for (uint32_t vidx_offset = 0; vidx_offset < vidx_block_size; ++vidx_offset) {
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
          const unsigned char* geno_read_iter = &(((const unsigned char*)genovec_buf)[sample_ctd4]);
          for (uint32_t vidx_offset = 0; vidx_offset < vidx_block_size; ++vidx_offset) {
            uint32_t cur_geno = *geno_read_iter;
            dosage_t* dosagebuf_write_iterx = &(dosagebuf_write_iter0[vidx_offset]);
            uint32_t sample_idx_lowbits = 0;
            while (1) {
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
        for (uint32_t vidx_offset = 0; vidx_offset < vidx_block_size; ++vidx_offset) {
          const uint32_t cur_dosage_ct = dosage_cts[vidx_offset];
          if (cur_dosage_ct) {
            const uintptr_t* dosage_present = &(dosagepresent_buf[vidx_offset * sample_ctaw]);
            const dosage_t* dosage_vals = &(dosagevals_buf[vidx_offset * sample_ct]);
            dosage_t* cur_dosage_write = &(smaj_dosagebuf_iter[vidx_offset]);
            uint32_t sample_idx = 0;
            for (uint32_t dosage_idx = 0; dosage_idx < cur_dosage_ct; ++dosage_idx, ++sample_idx) {
              next_set_unsafe_ck(dosage_present, &sample_idx);
              cur_dosage_write[sample_idx * stride] = dosage_vals[dosage_idx];
            }
          }
        }
        vidx_start = vidx_block_end;
        smaj_dosagebuf_iter = &(smaj_dosagebuf_iter[vidx_block_size]);
      } while (vidx_start != vidx_end);
    }
  dosage_transpose_thread_fail:
    if (is_last_block) {
      THREAD_RETURN;
    }
    THREAD_BLOCK_FINISH(tidx);
  }
}

static_assert(sizeof(dosage_t) == 2, "export_012_smaj() needs to be updated.");
pglerr_t export_012_smaj(const char* outname, const uintptr_t* orig_sample_include, const char* sample_ids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const uintptr_t* variant_include, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t include_dom, uint32_t include_uncounted, uint32_t max_thread_ct, uintptr_t pgr_alloc_cacheline_ct, char exportf_delim, pgen_file_info_t* pgfip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  threads_state_t ts;
  init_threads3z(&ts);
  FILE* outfile = nullptr;
  pglerr_t reterr = kPglRetSuccess;
  {
    // Write header line; then fully load-and-transpose the first X samples,
    // flush them, load-and-transpose the next X, etc.
    // Similar to export_ind_major_bed() and export_012_vmaj().
    // Priority is making the with-dosage case work well, since plink 1.9
    // already handles the no-dosage case.
    // (possible todo: have a separate no-dosage fast path)
    if (variant_ct * (1 + include_dom) > (kMaxLongLine - 4 * kMaxIdSlen - 64) / 8) {
      sprintf(g_logbuf, "Error: Too many variants for --export A%s.  (Try to work with A-transpose\ninstead.)\n", include_dom? "D" : "");
      goto export_012_smaj_ret_INCONSISTENT_INPUT_2;
    }
    if (fopen_checked(outname, FOPEN_WB, &outfile)) {
      goto export_012_smaj_ret_OPEN_FAIL;
    }
    char* writebuf = g_textbuf;
    if (max_allele_slen > kMaxMediumLine - 5) {
      if (bigstack_alloc_c(kMaxMediumLine + 5 + max_allele_slen, &writebuf)) {
        goto export_012_smaj_ret_NOMEM;
      }
    }
    char* write_iter = writebuf;
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    write_iter = memcpyl3a(write_iter, "FID");
    *write_iter++ = exportf_delim;
    write_iter = memcpyl3a(write_iter, "IID");
    *write_iter++ = exportf_delim;
    write_iter = memcpyl3a(write_iter, "PAT");
    *write_iter++ = exportf_delim;
    write_iter = memcpyl3a(write_iter, "MAT");
    *write_iter++ = exportf_delim;
    write_iter = memcpyl3a(write_iter, "SEX");
    *write_iter++ = exportf_delim;
    write_iter = strcpya(write_iter, "PHENOTYPE");
    uint32_t ref_allele_idx = 0;
    uint32_t variant_uidx = 0;
    uint64_t bytes_written = 0;
    for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx, ++variant_uidx) {
      next_set_unsafe_ck(variant_include, &variant_uidx);
      *write_iter++ = exportf_delim;
      const char* cur_var_id = variant_ids[variant_uidx];
      const uint32_t cur_slen = strlen(cur_var_id);
      write_iter = memcpyax(write_iter, cur_var_id, cur_slen, '_');
      uintptr_t variant_allele_idx_base = variant_uidx * 2;
      if (variant_allele_idxs) {
        variant_allele_idx_base = variant_allele_idxs[variant_uidx];
      }
      if (refalt1_select) {
        ref_allele_idx = refalt1_select[2 * variant_uidx];
      }
      write_iter = strcpya(write_iter, allele_storage[variant_allele_idx_base + ref_allele_idx]);
      if (write_iter >= writebuf_flush) {
        bytes_written += (uintptr_t)(write_iter - writebuf);
        if (fwrite_flush2(writebuf_flush, outfile, &write_iter)) {
          goto export_012_smaj_ret_WRITE_FAIL;
        }
      }
      if (include_uncounted) {
        write_iter = strcpya(write_iter, "(/");
        // todo: multiallelic case
        write_iter = strcpya(write_iter, allele_storage[variant_allele_idx_base + 1 - ref_allele_idx]);
        *write_iter++ = ')';
        if (write_iter >= writebuf_flush) {
          bytes_written += (uintptr_t)(write_iter - writebuf);
          if (fwrite_flush2(writebuf_flush, outfile, &write_iter)) {
            goto export_012_smaj_ret_WRITE_FAIL;
          }
        }
      }
      if (include_dom) {
        *write_iter++ = exportf_delim;
        write_iter = memcpya(write_iter, cur_var_id, cur_slen);
        write_iter = strcpya(write_iter, "_HET");
        if (write_iter >= writebuf_flush) {
          bytes_written += (uintptr_t)(write_iter - writebuf);
          if (fwrite_flush2(writebuf_flush, outfile, &write_iter)) {
            goto export_012_smaj_ret_WRITE_FAIL;
          }
        }
      }
    }
    append_binary_eoln(&write_iter);
    bytes_written += (uintptr_t)(write_iter - writebuf);
    if (bytes_written > kMaxLongLine) {
      sprintf(g_logbuf, "Error: --export A%s header line too long (>2GB).\n", include_dom? "D" : "");
      goto export_012_smaj_ret_INCONSISTENT_INPUT_2;
    }

    uintptr_t* pheno_nm = nullptr;
    uintptr_t* pheno_cc = nullptr;
    double* pheno_qt = nullptr;
    // .raw files don't support categorical phenotypes
    const uint32_t pheno_idx = first_cc_or_qt_pheno_idx(pheno_cols, pheno_ct);
    if (pheno_idx != UINT32_MAX) {
      const pheno_dtype_t type_code = pheno_cols[pheno_idx].type_code;
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
    // So, moving closer to export_ind_major_bed() strategy:
    // * Main read thread only loads raw bytes with pgfi_multiread().
    // * Worker thread(s) select subsets of the loaded variants to process in a
    //   cacheline-aware manner: first worker's interval ends with variant_idx
    //   divisible by (kCacheline / sizeof(dosage_t)) == 32 (unless it's the
    //   only active worker), and all other interval(s) start with variant_idx
    //   divisible by 32.

    // todo: check when this saturates
    uint32_t calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
    if (calc_thread_ct * kDosagePerCacheline > variant_ct) {
      calc_thread_ct = DIV_UP(variant_ct, kDosagePerCacheline);
    }
    unsigned char* main_loadbufs[2];
    uint32_t read_block_size;

    // temporarily hide 3/4 of memory
    unsigned char* bigstack_end_mark = g_bigstack_end;
    g_bigstack_end = &(g_bigstack_base[round_down_pow2(((uintptr_t)(g_bigstack_end - g_bigstack_base)) / 4, kCacheline)]);
    if (multithread_load_init(variant_include, raw_sample_ct, variant_ct, pgr_alloc_cacheline_ct, 0, 0, pgfip, &calc_thread_ct, nullptr, nullptr, nullptr, &read_block_size, main_loadbufs, &ts.threads, &g_pgr_ptrs, &g_read_variant_uidx_starts)) {
      goto export_012_smaj_ret_NOMEM;
    }
    g_bigstack_end = bigstack_end_mark;

    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    uintptr_t* sample_include;
    uint32_t* sample_include_cumulative_popcounts;
    if (bigstack_alloc_ul(raw_sample_ctl, &sample_include) ||
        bigstack_alloc_ui(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
        bigstack_alloc_ui(calc_thread_ct + 1, &g_write_vidx_starts) ||
        bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_genovecs) ||
        bigstack_alloc_ulp(calc_thread_ct, &g_thread_write_dosagepresents) ||
        bigstack_alloc_dosagep(calc_thread_ct, &g_thread_write_dosagevals)) {
      goto export_012_smaj_ret_NOMEM;
    }

    // Remaining memory byte requirements:
    //   calc_thread_ct * kDosagePerCacheline * kBytesPerWord *
    //     read_sample_ctaw2 for per-thread genovecs buffers
    //   calc_thread_ct * kDosagePerCacheline * kBytesPerWord *
    //     read_sample_ctaw for per-thread dosage_presents buffers
    //   calc_thread_ct * kDosagePerCacheline * read_sample_ct *
    //     sizeof(dosage_t) for per-thread dosage_vals buffers
    //   (dosage_ct buffers just go on the thread stacks)
    //   read_sample_ct * variant_ct * sizeof(dosage_t) for g_smaj_dosagebuf
    // This is about
    //   read_sample_ct *
    //     (calc_thread_ct * kDosagePerCacheline * 2.375 + variant_ct * 2)
    ts.calc_thread_ct = calc_thread_ct;
    uintptr_t bytes_avail = bigstack_left();
    // account for rounding
    if (bytes_avail < kCacheline + calc_thread_ct * kDosagePerCacheline * (2 * kCacheline)) {
      goto export_012_smaj_ret_NOMEM;
    }
    bytes_avail -= kCacheline + calc_thread_ct * kDosagePerCacheline * (2 * kCacheline);
    uint32_t read_sample_ct = sample_ct;
    uint32_t pass_ct = 1;
    const uintptr_t bytes_per_sample = calc_thread_ct * (kDosagePerCacheline / 8) * (3LLU + 8 * sizeof(dosage_t)) + variant_ct * sizeof(dosage_t);
    if ((sample_ct * ((uint64_t)bytes_per_sample)) > bytes_avail) {
      read_sample_ct = bytes_avail / bytes_per_sample;
      if (!read_sample_ct) {
        goto export_012_smaj_ret_NOMEM;
      }
      if (read_sample_ct > 4) {
        read_sample_ct = round_down_pow2(read_sample_ct, 4);
      }
      pass_ct = 1 + (sample_ct - 1) / read_sample_ct;
    }
    uintptr_t read_sample_ctaw = BITCT_TO_ALIGNED_WORDCT(read_sample_ct);
    uintptr_t read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
    for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
      g_thread_write_genovecs[tidx] = (uintptr_t*)bigstack_alloc_raw(kDosagePerCacheline * sizeof(intptr_t) * read_sample_ctaw2);
      g_thread_write_dosagepresents[tidx] = (uintptr_t*)bigstack_alloc_raw(kDosagePerCacheline * sizeof(intptr_t) * read_sample_ctaw);
      g_thread_write_dosagevals[tidx] = (dosage_t*)bigstack_alloc_raw(kDosagePerCacheline * sizeof(dosage_t) * read_sample_ct);
    }
    g_variant_include = variant_include;
    g_refalt1_select = refalt1_select;
    g_calc_thread_ct = calc_thread_ct;
    g_sample_ct = read_sample_ct;
    g_stride = round_up_pow2(variant_ct, kDosagePerCacheline);
    g_smaj_dosagebuf = (dosage_t*)bigstack_alloc_raw_rd(read_sample_ct * ((uintptr_t)g_stride) * sizeof(dosage_t));
    g_error_ret = kPglRetSuccess;

    const uint32_t read_block_sizel = BITCT_TO_WORDCT(read_block_size);
    const uint32_t read_block_ct_m1 = (raw_variant_ct - 1) / read_block_size;
    uint32_t sample_uidx_start = next_set_unsafe(orig_sample_include, 0);
    for (uint32_t pass_idx = 0; pass_idx < pass_ct; ++pass_idx) {
      memcpy(sample_include, orig_sample_include, raw_sample_ctl * sizeof(intptr_t));
      if (sample_uidx_start) {
        clear_bits_nz(0, sample_uidx_start, sample_include);
      }
      uint32_t sample_uidx_end;
      if (pass_idx + 1 == pass_ct) {
        read_sample_ct = sample_ct - pass_idx * read_sample_ct;
        g_sample_ct = read_sample_ct;
        sample_uidx_end = raw_sample_ct;
        read_sample_ctaw = BITCT_TO_ALIGNED_WORDCT(read_sample_ct);
        read_sample_ctaw2 = QUATERCT_TO_ALIGNED_WORDCT(read_sample_ct);
      } else {
        sample_uidx_end = jump_forward_set_unsafe(orig_sample_include, sample_uidx_start + 1, read_sample_ct);
        clear_bits_nz(sample_uidx_end, raw_sample_ct, sample_include);
      }
      fill_cumulative_popcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
      g_sample_include = sample_include;
      g_sample_include_cumulative_popcounts = sample_include_cumulative_popcounts;
      if (pass_idx) {
        reinit_threads3z(&ts);
        pgfip->block_base = main_loadbufs[0];
        for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
          pgr_clear_ld_cache(g_pgr_ptrs[tidx]);
          g_pgr_ptrs[tidx]->fi.block_base = main_loadbufs[0];
          g_pgr_ptrs[tidx]->fi.block_offset = 0;
        }
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
      uint32_t variant_idx = 0;
      uint32_t cur_read_block_size = read_block_size;
      uint32_t pct = 0;
      uint32_t next_print_idx = variant_ct / 100;
      while (1) {
        uintptr_t cur_block_write_ct = 0;
        if (!ts.is_last_block) {
          while (read_block_idx < read_block_ct_m1) {
            cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), read_block_sizel);
            if (cur_block_write_ct) {
              break;
            }
            ++read_block_idx;
          }
          if (read_block_idx == read_block_ct_m1) {
            cur_read_block_size = raw_variant_ct - (read_block_idx * read_block_size);
            cur_block_write_ct = popcount_longs(&(variant_include[read_block_idx * read_block_sizel]), BITCT_TO_WORDCT(cur_read_block_size));
          }
          if (pgfi_multiread(variant_include, read_block_idx * read_block_size, read_block_idx * read_block_size + cur_read_block_size, cur_block_write_ct, pgfip)) {
            goto export_012_smaj_ret_THREAD_CREATE_FAIL;
          }
        }
        if (variant_idx) {
          join_threads3z(&ts);
          reterr = g_error_ret;
          if (reterr) {
            if (reterr == kPglRetMalformedInput) {
              logprint("\n");
              logerrprint("Error: Malformed .pgen file.\n");
            }
            goto export_012_smaj_ret_1;
          }
        }
        if (!ts.is_last_block) {
          g_cur_block_write_ct = cur_block_write_ct;
          compute_partition_aligned(variant_include, calc_thread_ct, read_block_idx * read_block_size, variant_idx, cur_block_write_ct, kDosagePerCacheline, g_read_variant_uidx_starts, g_write_vidx_starts);
          for (uint32_t tidx = 0; tidx < calc_thread_ct; ++tidx) {
            g_pgr_ptrs[tidx]->fi.block_base = pgfip->block_base;
            g_pgr_ptrs[tidx]->fi.block_offset = pgfip->block_offset;
          }
          ts.is_last_block = (variant_idx + cur_block_write_ct == variant_ct);
          ts.thread_func_ptr = dosage_transpose_thread;
          if (spawn_threads3z(read_block_idx, &ts)) {
            goto export_012_smaj_ret_THREAD_CREATE_FAIL;
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
          next_print_idx = (pct * ((uint64_t)variant_ct)) / 100;
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
      uint32_t sample_uidx = sample_uidx_start;
      const dosage_t* cur_dosage_row = g_smaj_dosagebuf;
      for (uint32_t sample_idx = 0; sample_idx < read_sample_ct; ++sample_idx, ++sample_uidx) {
        next_set_unsafe_ck(sample_include, &sample_uidx);
        const char* cur_sample_fid = &(sample_ids[sample_uidx * max_sample_id_blen]);
        const char* fid_end = (const char*)rawmemchr(cur_sample_fid, '\t');
        write_iter = memcpyax(write_iter, cur_sample_fid, (uintptr_t)(fid_end - cur_sample_fid), exportf_delim);
        write_iter = strcpyax(write_iter, &(fid_end[1]), exportf_delim);
        write_iter = strcpyax(write_iter, &(paternal_ids[sample_uidx * max_paternal_id_blen]), exportf_delim);
        write_iter = strcpyax(write_iter, &(maternal_ids[sample_uidx * max_maternal_id_blen]), exportf_delim);
        *write_iter++ = sexchar(sex_nm, sex_male, sample_uidx);
        *write_iter++ = exportf_delim;
        if ((!pheno_nm) || (!IS_SET(pheno_nm, sample_uidx))) {
          write_iter = memcpya(write_iter, legacy_output_missing_pheno, lomp_slen);
        } else if (pheno_cc) {
          *write_iter++ = '1' + IS_SET(pheno_cc, sample_uidx);
        } else {
          write_iter = dtoa_g(pheno_qt[sample_uidx], write_iter);
        }
        if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
          goto export_012_smaj_ret_WRITE_FAIL;
        }
        for (uint32_t variant_idx = 0; variant_idx < variant_ct; ++variant_idx) {
          *write_iter++ = exportf_delim;
          uint32_t cur_dosage_val = cur_dosage_row[variant_idx];
          if (cur_dosage_val != 65535) {
            write_iter = print_small_dosage(cur_dosage_val, write_iter);
            if (include_dom) {
              *write_iter++ = exportf_delim;
              write_iter = print_small_dosage(16384 - abs_int32(cur_dosage_val - 16384), write_iter);
            }
          } else {
            write_iter = strcpya(write_iter, "NA");
            if (include_dom) {
              *write_iter++ = exportf_delim;
              write_iter = strcpya(write_iter, "NA");
            }
          }
          // todo: try making this check less frequently
          if (fwrite_ck(writebuf_flush, outfile, &write_iter)) {
            goto export_012_smaj_ret_WRITE_FAIL;
          }
        }
        append_binary_eoln(&write_iter);
        cur_dosage_row = &(cur_dosage_row[g_stride]);
        if (sample_idx >= next_print_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (sample_idx * 100LLU) / read_sample_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_idx = (pct * ((uint64_t)read_sample_ct)) / 100;
        }
      }
      sample_uidx_start = sample_uidx_end;
      if (pct > 10) {
        fputs("\b \b", stdout);
      }
    }
    if (fclose_flush_null(writebuf_flush, write_iter, &outfile)) {
      goto export_012_smaj_ret_WRITE_FAIL;
    }
    fputs("\b\bdone.\n", stdout);
    LOGPRINTFWW("--export A%s: %s written.\n", include_dom? "D" : "", outname);
  }
  while (0) {
  export_012_smaj_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  export_012_smaj_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  export_012_smaj_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  export_012_smaj_ret_INCONSISTENT_INPUT_2:
    logerrprintb();
    reterr = kPglRetInconsistentInput;
    break;
  export_012_smaj_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 export_012_smaj_ret_1:
  threads3z_cleanup(&ts, &g_cur_block_write_ct);
  fclose_cond(outfile);
  pgfip->block_base = nullptr;
  bigstack_reset(bigstack_mark);
  return reterr;
}

pglerr_t exportf(char* xheader, const uintptr_t* sample_include, const char* sample_ids, const char* sids, const char* paternal_ids, const char* maternal_ids, const uintptr_t* sex_nm, const uintptr_t* sex_male, const pheno_col_t* pheno_cols, const char* pheno_names, const uintptr_t* variant_include, const chr_info_t* cip, const uint32_t* variant_bps, char** variant_ids, const uintptr_t* variant_allele_idxs, char** allele_storage, const alt_allele_ct_t* refalt1_select, const uintptr_t* pvar_qual_present, const float* pvar_quals, const uintptr_t* pvar_filter_present, const uintptr_t* pvar_filter_npass, char** pvar_filter_storage, const char* pvar_info_reload, const double* variant_cms, uintptr_t xheader_blen, uint32_t xheader_info_pr, uint32_t xheader_info_pr_nonflag, uint32_t raw_sample_ct, uint32_t sample_ct, uintptr_t max_sample_id_blen, uintptr_t max_sid_blen, uintptr_t max_paternal_id_blen, uintptr_t max_maternal_id_blen, uint32_t pheno_ct, uintptr_t max_pheno_name_blen, uint32_t raw_variant_ct, uint32_t variant_ct, uint32_t max_allele_slen, uint32_t max_filter_slen, uint32_t info_reload_slen, uint32_t max_thread_ct, make_plink2_t make_plink2_modifier, exportf_flags_t exportf_modifier, idpaste_t exportf_id_paste, char exportf_id_delim, __maybe_unused uint32_t exportf_bits, uintptr_t pgr_alloc_cacheline_ct, pgen_file_info_t* pgfip, pgen_reader_t* simple_pgrp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  pglerr_t reterr = kPglRetSuccess;
  {
    const uint32_t raw_sample_ctl = BITCT_TO_WORDCT(raw_sample_ct);
    const uint32_t sample_ctaw = BITCT_TO_ALIGNED_WORDCT(sample_ct);
    const uint32_t sample_ctl = BITCT_TO_WORDCT(sample_ct);
    uint32_t* sample_include_cumulative_popcounts;
    uintptr_t* sex_male_collapsed;
    if (bigstack_alloc_ui(raw_sample_ctl, &sample_include_cumulative_popcounts) ||
        bigstack_alloc_ul(sample_ctaw, &sex_male_collapsed)) {
      goto exportf_ret_NOMEM;
    }
    fill_cumulative_popcounts(sample_include, raw_sample_ctl, sample_include_cumulative_popcounts);
    copy_bitarr_subset(sex_male, sample_include, sample_ct, sex_male_collapsed);
    zero_trailing_words(sample_ctl, sex_male_collapsed);
    uint32_t* sample_missing_geno_cts = nullptr;
    if (exportf_modifier & (kfExportfOxGen | kfExportfHaps | kfExportfHapsLegend | kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13)) {
      if (bigstack_alloc_ui(sample_ct, &sample_missing_geno_cts)) {
        goto exportf_ret_NOMEM;
      }
    }
    if (exportf_modifier & (kfExportf01 | kfExportf12)) {
      // todo
    }
    if (exportf_modifier & (kfExportfTypemask - kfExportfIndMajorBed - kfExportfVcf - kfExportfOxGen - kfExportfBgen11 - kfExportfHaps - kfExportfHapsLegend - kfExportfATranspose - kfExportfA - kfExportfAD)) {
      logerrprint("Error: Only VCF, oxford, bgen-1.1, haps, hapslegend, A, AD, A-transpose, and\nind-major-bed output have been implemented so far.\n");
      reterr = kPglRetNotYetSupported;
      goto exportf_ret_1;
    }
    const char exportf_delim = (exportf_modifier & kfExportfSpaces)? ' ' : '\t';
    if (exportf_modifier & kfExportfATranspose) {
      strcpy(outname_end, ".traw");
      pgr_clear_ld_cache(simple_pgrp);
      reterr = export_012_vmaj(outname, sample_include, sample_include_cumulative_popcounts, sample_ids, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, variant_cms, sample_ct, max_sample_id_blen, variant_ct, max_allele_slen, exportf_delim, simple_pgrp);
      if (reterr) {
        goto exportf_ret_1;
      }
    }
    if (exportf_modifier & kfExportfIndMajorBed) {
      reterr = export_ind_major_bed(sample_include, variant_include, variant_allele_idxs, refalt1_select, raw_sample_ct, sample_ct, raw_variant_ct, variant_ct, max_thread_ct, pgr_alloc_cacheline_ct, pgfip, outname, outname_end);
      if (reterr) {
        goto exportf_ret_1;
      }
    }
    if (exportf_modifier & kfExportfOxGen) {
      strcpy(outname_end, ".gen");
      pgr_clear_ld_cache(simple_pgrp);
      reterr = export_ox_gen(outname, sample_include, sample_include_cumulative_popcounts, sex_male, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, sample_ct, variant_ct, max_allele_slen, exportf_modifier, simple_pgrp, sample_missing_geno_cts);
      if (reterr) {
        goto exportf_ret_1;
      }
    }
    if (exportf_modifier & (kfExportfHaps | kfExportfHapsLegend)) {
      pgr_clear_ld_cache(simple_pgrp);
      reterr = export_ox_hapslegend(sample_include, sample_include_cumulative_popcounts, sex_male_collapsed, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, exportf_modifier, simple_pgrp, outname, outname_end);
      if (reterr) {
        goto exportf_ret_1;
      }
      fill_uint_zero(sample_ct, sample_missing_geno_cts);
    }
    if (exportf_modifier & kfExportfBgen11) {
      assert(popcount_longs(sample_include, raw_sample_ctl) == sample_ct);
      strcpy(outname_end, ".bgen");
      reterr = export_bgen11(outname, sample_include, sample_include_cumulative_popcounts, sex_male, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, exportf_modifier, pgr_alloc_cacheline_ct, pgfip, sample_missing_geno_cts);
      if (reterr) {
        goto exportf_ret_1;
      }
      /*
    } else if (exportf_modifier & (kfExportfBgen12 | kfExportfBgen13)) {
      strcpy(outname_end, ".bgen");
      reterr = export_bgen13(outname, sample_include, sample_include_cumulative_popcounts, sex_male, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, sample_ct, raw_variant_ct, variant_ct, max_allele_slen, max_thread_ct, exportf_modifier, exportf_bits, pgr_alloc_cacheline_ct, pgfip, sample_missing_geno_cts);
      if (reterr) {
        goto exportf_ret_1;
      }
      */
    }
    if (exportf_modifier & (kfExportfOxGen | kfExportfBgen11 | kfExportfBgen12 | kfExportfBgen13 | kfExportfHaps | kfExportfHapsLegend)) {
      strcpy(outname_end, ".sample");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      uint32_t y_ct = 0;
      int32_t y_code = cip->xymt_codes[kChrOffsetY];
      if ((y_code >= 0) && is_set(cip->chr_mask, y_code)) {
        y_ct = count_chr_variants_unsafe(variant_include, cip, y_code);
      }
      assert(popcount_longs(sample_include, raw_sample_ctl) == sample_ct);
      reterr = export_ox_sample(outname, sample_include, sample_ids, sample_missing_geno_cts, sex_nm, sex_male, pheno_cols, pheno_names, sample_ct, max_sample_id_blen, pheno_ct, max_pheno_name_blen, variant_ct, y_ct);
      if (reterr) {
        goto exportf_ret_1;
      }
      logprint("done.\n");
    }
    if (exportf_modifier & kfExportfVcf) {
      pgr_clear_ld_cache(simple_pgrp);
      reterr = export_vcf(xheader, sample_include, sample_include_cumulative_popcounts, sample_ids, sids, sex_male_collapsed, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, pvar_qual_present, pvar_quals, pvar_filter_present, pvar_filter_npass, pvar_filter_storage, pvar_info_reload, xheader_blen, xheader_info_pr, xheader_info_pr_nonflag, sample_ct, max_sample_id_blen, max_sid_blen, raw_variant_ct, variant_ct, max_allele_slen, max_filter_slen, info_reload_slen, max_thread_ct, exportf_modifier, exportf_id_paste, exportf_id_delim, pgfip, simple_pgrp, outname, outname_end);
      if (reterr) {
        goto exportf_ret_1;
      }
    }
    // todo: everything else
    // sample-major output should share a (probably multithreaded) transpose
    // routine
    if (exportf_modifier & (kfExportfA | kfExportfAD)) {
      strcpy(outname_end, ".raw");
      reterr = export_012_smaj(outname, sample_include, sample_ids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, variant_include, variant_ids, variant_allele_idxs, allele_storage, refalt1_select, raw_sample_ct, sample_ct, max_sample_id_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, raw_variant_ct, variant_ct, max_allele_slen, (exportf_modifier / kfExportfAD) & 1, (exportf_modifier / kfExportfIncludeAlt) & 1, max_thread_ct, pgr_alloc_cacheline_ct, exportf_delim, pgfip);
      if (reterr) {
        goto exportf_ret_1;
      }
    }

    if ((!(make_plink2_modifier & kfMakeFam)) && (exportf_modifier & kfExportfIndMajorBed)) {
      strcpy(outname_end, ".fam");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = write_fam(outname, sample_include, sample_ids, paternal_ids, maternal_ids, sex_nm, sex_male, pheno_cols, nullptr, sample_ct, max_sample_id_blen, max_paternal_id_blen, max_maternal_id_blen, pheno_ct, exportf_delim);
      if (reterr) {
        goto exportf_ret_1;
      }
      logprint("done.\n");
    }
    if ((!(make_plink2_modifier & kfMakeBim)) && (exportf_modifier & kfExportfIndMajorBed)) {
      strcpy(outname_end, ".bim");
      LOGPRINTFWW5("Writing %s ... ", outname);
      fflush(stdout);
      reterr = write_map_or_bim(outname, variant_include, cip, variant_bps, variant_ids, variant_allele_idxs, allele_storage, nullptr, refalt1_select, variant_cms, variant_ct, max_allele_slen, exportf_delim, 0, max_thread_ct);
      if (reterr) {
        goto exportf_ret_1;
      }
      logprint("done.\n");
    }
  }
  while (0) {
  exportf_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  }
 exportf_ret_1:
  bigstack_reset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
