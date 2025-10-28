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

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>  // unlink()

#include "include/pgenlib_misc.h"
#include "include/pgenlib_read.h"
#include "include/pgenlib_write.h"
#include "include/plink2_base.h"
#include "include/plink2_bits.h"
#include "include/plink2_string.h"
#include "include/plink2_text.h"
#include "include/plink2_thread.h"
#include "plink2_cmdline.h"
#include "plink2_common.h"
#include "plink2_compress_stream.h"
#include "plink2_data.h"
#include "plink2_decompress.h"
#include "plink2_psam.h"

// This covers formats that are fully supported by PLINK 1.x (no multiallelic
// variants, dosages, or phase information).

#ifdef __cplusplus
namespace plink2 {
#endif

// Counts number of variants (both raw and with chromosome filter applied),
// checks whether a nonzero CM column exists, and allocates+returns the
// chromosome-and-negative-bp-filter bitvector.
// Assumes FinalizeChrset() has already been called.
// Errors out on .bim file.
// variant_include assumed to be initialized to nullptr.
PglErr ScanMap(const char* mapname, MiscFlags misc_flags, LoadFilterLogFlags load_filter_log_import_flags, ChrInfo* cip, uint32_t* raw_variant_ct_ptr, uint32_t* variant_ct_ptr, uint32_t* neg_bp_seen_ptr, uint32_t* at_least_one_nzero_cm_ptr, uintptr_t** variant_include_ptr) {
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream map_txs;
  PreinitTextStream(&map_txs);
  {
    // Even though no .map line should exceed kMaxMediumLine bytes (128 KiB) in
    // practice, TextStream actually doesn't permit a max_line_blen value
    // smaller than kDecompressMinBlen (1 MiB), and there is little reason to
    // use a value below kTextStreamBlenFast (11 MiB) in practice unless memory
    // is very tight.
    reterr = InitTextStreamEx(mapname, 1, kTextStreamBlenFast, kTextStreamBlenFast, 1, &map_txs);
    if (unlikely(reterr)) {
      goto ScanMap_ret_TSTREAM_FAIL;
    }
    uintptr_t raw_variant_ct_limit = CHAR_BIT * S_CAST(uintptr_t, g_bigstack_end - g_bigstack_base);
    if (raw_variant_ct_limit > kPglMaxVariantCt) {
      raw_variant_ct_limit = kPglMaxVariantCt;
    }
    uintptr_t* variant_include = R_CAST(uintptr_t*, g_bigstack_base);

    char* line_start;
    do {
      ++line_idx;
      line_start = TextGet(&map_txs);
      if (unlikely(!line_start)) {
        if (!TextStreamErrcode2(&map_txs, &reterr)) {
          logerrputs("Error: Empty .map file.\n");
          goto ScanMap_ret_DEGENERATE_DATA;
        }
        goto ScanMap_ret_TSTREAM_FAIL;
      }
    } while (*line_start == '#');
    uint32_t map_cols = 3;
    {
      const char* linebuf_iter = NextTokenMult(line_start, 2);
      if (unlikely(!linebuf_iter)) {
        goto ScanMap_ret_MISSING_TOKENS;
      }
      linebuf_iter = NextToken(linebuf_iter);
      if (linebuf_iter) {
        linebuf_iter = NextToken(linebuf_iter);
        if (likely(!linebuf_iter)) {
          map_cols = 4;
        } else {
          // do NOT permit >4 columns here.
          snprintf(g_logbuf, kLogbufSize, "Error: %s is not a .map file (too many columns).\n", mapname);
          goto ScanMap_ret_MALFORMED_INPUT_WW;
        }
      }
    }

    const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
    uint32_t at_least_one_nzero_cm = 0;
    uint32_t raw_variant_ct = 0;
    uint32_t variant_ct = 0;
    uint32_t neg_bp_seen = 0;
    uintptr_t variant_include_word = 0;
    uint32_t variant_uidx_lowbits = 0;
    char* line_iter = line_start;
    while (1) {
      if (unlikely(raw_variant_ct == raw_variant_ct_limit)) {
        if (raw_variant_ct_limit == kPglMaxVariantCt) {
          logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
          goto ScanMap_ret_MALFORMED_INPUT;
        }
        goto ScanMap_ret_NOMEM;
      }
      {
        // chrom, id, (cm?), pos
        char* chr_code_end = CurTokenEnd(line_iter);
        char* variant_id_start = FirstNonTspace(chr_code_end);
        uint32_t cur_chr_code;
        reterr = GetOrAddChrCodeDestructive(".map file", line_idx, prohibit_extra_chr, line_iter, chr_code_end, cip, &cur_chr_code);
        if (unlikely(reterr)) {
          goto ScanMap_ret_1;
        }
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          line_iter = variant_id_start;
          goto ScanMap_skip_variant;
        }
        char* third_token_start = FirstNonTspace(FirstSpaceOrEoln(variant_id_start));
        char* bp_start = third_token_start;
        if (map_cols == 4) {
          bp_start = FirstNonTspace(FirstSpaceOrEoln(third_token_start));
        }
        if (IsEolnKns(*bp_start)) {
          goto ScanMap_ret_MISSING_TOKENS;
        }
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(bp_start, &cur_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, mapname);
          goto ScanMap_ret_MALFORMED_INPUT_WW;
        }
        line_iter = bp_start;
        if (cur_bp < 0) {
          neg_bp_seen = 1;
          goto ScanMap_skip_variant;
        }
        if ((map_cols == 4) && (!at_least_one_nzero_cm)) {
          double cur_cm;
          char* cm_end = ScantokDouble(third_token_start, &cur_cm);
          if (unlikely(!cm_end)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, mapname);
            goto ScanMap_ret_MALFORMED_INPUT_WW;
          }
          at_least_one_nzero_cm = (cur_cm != 0.0);
        }
        variant_include_word |= k1LU << variant_uidx_lowbits;
        ++variant_ct;
      }
    ScanMap_skip_variant:
      if (++variant_uidx_lowbits == kBitsPerWord) {
        variant_include[raw_variant_ct / kBitsPerWord] = variant_include_word;
        variant_uidx_lowbits = 0;
        variant_include_word = 0;
      }
      ++raw_variant_ct;
      do {
        ++line_idx;
        line_iter = AdvPastDelim(line_iter, '\n');
        if (!TextGetUnsafe2(&map_txs, &line_iter)) {
          if (!TextStreamErrcode2(&map_txs, &reterr)) {
            goto ScanMap_eof;
          }
          goto ScanMap_ret_TSTREAM_FAIL;
        }
        // bugfix (4 Aug 2024): original .map specification permits interior
        // comment lines
      } while (line_iter[0] == '#');
    }
  ScanMap_eof:
    if (unlikely(!variant_ct)) {
      char* write_iter = strcpya_k(g_logbuf, "Error: All ");
      write_iter = u32toa(raw_variant_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (raw_variant_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " in .map file excluded by ");
      if (load_filter_log_import_flags) {
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
      }
      if (neg_bp_seen) {
        if (load_filter_log_import_flags) {
          write_iter = strcpya_k(write_iter, " and/or ");
        }
        write_iter = strcpya_k(write_iter, "negative bp coordinates");
      }
      strcpy_k(write_iter, ".\n");
      goto ScanMap_ret_INCONSISTENT_INPUT_WW;
    }
    *raw_variant_ct_ptr = raw_variant_ct;
    *variant_ct_ptr = variant_ct;
    *neg_bp_seen_ptr = neg_bp_seen;
    *at_least_one_nzero_cm_ptr = at_least_one_nzero_cm;
    if (raw_variant_ct != variant_ct) {
      if (variant_uidx_lowbits) {
        variant_include[raw_variant_ct / kBitsPerWord] = variant_include_word;
      }
      BigstackFinalizeW(variant_include, DivUp(raw_variant_ct, kBitsPerWord));
      *variant_include_ptr = variant_include;
    }
  }
  while (0) {
  ScanMap_ret_TSTREAM_FAIL:
    TextStreamErrPrint(mapname, &map_txs);
    break;
  ScanMap_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  ScanMap_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, mapname);
  ScanMap_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  ScanMap_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  ScanMap_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  ScanMap_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 ScanMap_ret_1:
  CleanupTextStream2(mapname, &map_txs, &reterr);
  BigstackEndReset(bigstack_end_mark);
  return reterr;
}

// Assumes ScanMap() was previously called on the .map.
PglErr MapToPvar(const char* mapname, const ChrInfo* cip, const char* const* allele_storage, uint32_t variant_ct, uint32_t max_allele_slen, ImportFlags import_flags, uint32_t at_least_one_nzero_cm, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream map_txs;
  PreinitTextStream(&map_txs);

  char* cswritep = nullptr;
  CompressStreamState css;
  PreinitCstream(&css);
  {
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    const uint32_t output_zst = (import_flags / kfImportKeepAutoconvVzs) & 1;
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + 2 * kMaxIdSlen + 32 + 2 * max_allele_slen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &css, &cswritep);
    if (unlikely(reterr)) {
      goto MapToPvar_ret_1;
    }

    reterr = InitTextStreamEx(mapname, 0, kTextStreamBlenFast, kTextStreamBlenFast, 1, &map_txs);
    if (unlikely(reterr)) {
      goto MapToPvar_ret_TSTREAM_REWIND_FAIL;
    }

    char* line_start;
    do {
      ++line_idx;
      line_start = TextGet(&map_txs);
      if (unlikely(!line_start)) {
        goto MapToPvar_ret_TSTREAM_REWIND_FAIL;
      }
    } while (*line_start == '#');
    uint32_t map_cols = 3;
    {
      const char* linebuf_iter = NextTokenMult(line_start, 2);
      if (unlikely(!linebuf_iter)) {
        goto MapToPvar_ret_REWIND_FAIL;
      }
      linebuf_iter = NextToken(linebuf_iter);
      if (linebuf_iter) {
        linebuf_iter = NextToken(linebuf_iter);
        if (likely(!linebuf_iter)) {
          map_cols = 4;
        } else {
          goto MapToPvar_ret_REWIND_FAIL;
        }
      }
    }
    cswritep = strcpya_k(cswritep, "#CHROM\tPOS\tID\tREF\tALT");
    if (at_least_one_nzero_cm) {
      cswritep = strcpya_k(cswritep, "\tCM");
    }
    AppendBinaryEoln(&cswritep);

    char* line_iter = line_start;
    for (uint32_t variant_idx = 0; ; ) {
      {
        // chrom, id, (cm?), pos
        char* chr_code_end = CurTokenEnd(line_iter);
        const uint32_t cur_chr_code = GetChrCodeCounted(cip, chr_code_end - line_iter, line_iter);
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          line_iter = chr_code_end;
          goto MapToPvar_skip_variant;
        }
        char* variant_id_start = FirstNonTspace(chr_code_end);
        char* variant_id_end = CurTokenEnd(variant_id_start);
        char* third_token_start = FirstNonTspace(variant_id_end);
        char* bp_start = third_token_start;
        if (map_cols == 4) {
          bp_start = FirstNonTspace(FirstSpaceOrEoln(third_token_start));
        }
        if (IsEolnKns(*bp_start)) {
          goto MapToPvar_ret_REWIND_FAIL;
        }
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(bp_start, &cur_bp))) {
          goto MapToPvar_ret_REWIND_FAIL;
        }
        line_iter = bp_start;
        if (cur_bp < 0) {
          goto MapToPvar_skip_variant;
        }
        cswritep = chrtoa(cip, cur_chr_code, cswritep);
        *cswritep++ = '\t';
        cswritep = u32toa_x(cur_bp, '\t', cswritep);
        const uint32_t variant_id_slen = variant_id_end - variant_id_start;
        if (unlikely(variant_id_slen > kMaxIdSlen)) {
          logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
          goto MapToPvar_ret_MALFORMED_INPUT;
        }
        cswritep = memcpyax(cswritep, variant_id_start, variant_id_slen, '\t');
        const char* ref_allele = allele_storage[variant_idx * 2];
        if (ref_allele) {
          cswritep = strcpya(cswritep, ref_allele);
        } else {
          *cswritep++ = '.';
        }
        *cswritep++ = '\t';
        const char* alt_allele = allele_storage[variant_idx * 2 + 1];
        if (alt_allele) {
          cswritep = strcpya(cswritep, alt_allele);
        } else {
          *cswritep++ = '.';
        }
        if (at_least_one_nzero_cm) {
          double cur_cm;
          char* cm_end = ScantokDouble(third_token_start, &cur_cm);
          if (unlikely(!cm_end)) {
            logerrprintfww("Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, mapname);
            goto MapToPvar_ret_MALFORMED_INPUT;
          }
          *cswritep++ = '\t';
          cswritep = dtoa_g_p8(cur_cm, cswritep);
        }
        AppendBinaryEoln(&cswritep);
        if (unlikely(Cswrite(&css, &cswritep))) {
          goto MapToPvar_ret_WRITE_FAIL;
        }
        ++variant_idx;
        if (variant_idx == variant_ct) {
          break;
        }
      }
    MapToPvar_skip_variant:
      do {
        ++line_idx;
        line_iter = AdvPastDelim(line_iter, '\n');
        if (unlikely(!TextGetUnsafe2(&map_txs, &line_iter))) {
          TextStreamErrcode2(&map_txs, &reterr);
          goto MapToPvar_ret_TSTREAM_REWIND_FAIL;
        }
      } while (line_iter[0] == '#');
    }
  }
  while (0) {
  MapToPvar_ret_TSTREAM_REWIND_FAIL:
    TextStreamErrPrintRewind(".map file", &map_txs, &reterr);
    break;
  MapToPvar_ret_REWIND_FAIL:
    logerrprintfww(kErrprintfRewind, ".map file");
    reterr = kPglRetRewindFail;
    break;
  MapToPvar_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  MapToPvar_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  }
 MapToPvar_ret_1:
  CswriteCloseCond(&css, cswritep);
  CleanupTextStream2(mapname, &map_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

// could add an option to LoadPvar() to not require allele columns, but .map
// is easy enough to write a separate loader for...
CONSTI32(kLoadMapBlockSize, 65536);

// assumes FinalizeChrset() has already been called.
// .bim ok
PglErr LoadMap(const char* mapname, MiscFlags misc_flags, LoadFilterLogFlags load_filter_log_import_flags, ChrInfo* cip, uint32_t* max_variant_id_slen_ptr, ChrIdx** variant_chr_codes_ptr, uint32_t** variant_bps_ptr, char*** variant_ids_ptr, double** variant_cms_ptr, uint32_t* variant_ct_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  uintptr_t line_idx = 0;
  PglErr reterr = kPglRetSuccess;
  TextStream map_txs;
  PreinitTextStream(&map_txs);
  {
    // Workspace used as follows:
    // |--linebuf--|--temp-->----|----<- variant IDs --|
    //            1/4                                 end
    // linebuf is overwritten with the main return arrays at the end.
    reterr = SizeAndInitTextStream(mapname, bigstack_left() / 4, 1, &map_txs);
    if (unlikely(reterr)) {
      goto LoadMap_ret_TSTREAM_FAIL;
    }
    char* line_start;
    do {
      ++line_idx;
      line_start = TextGet(&map_txs);
      if (unlikely(!line_start)) {
        if (!TextStreamErrcode2(&map_txs, &reterr)) {
          logerrputs("Error: Empty .map file.\n");
          goto LoadMap_ret_DEGENERATE_DATA;
        }
        goto LoadMap_ret_TSTREAM_FAIL;
      }
    } while (*line_start == '#');
    uint32_t map_cols = 3;
    {
      const char* linebuf_iter = NextTokenMult(line_start, 2);
      if (unlikely(!linebuf_iter)) {
        goto LoadMap_ret_MISSING_TOKENS;
      }
      linebuf_iter = NextToken(linebuf_iter);
      if (linebuf_iter) {
        linebuf_iter = NextToken(linebuf_iter);
        if (!linebuf_iter) {
          map_cols = 4;
        } else {
          linebuf_iter = NextToken(linebuf_iter);
          if (linebuf_iter) {
            if (unlikely(NextToken(linebuf_iter))) {
              // do NOT permit >6 columns, .bim is ok but .pvar is not
              // (pointless to support .pvar for legacy formats)
              snprintf(g_logbuf, kLogbufSize, "Error: %s is not a .map/.bim file (too many columns).\n", mapname);
              goto LoadMap_ret_MALFORMED_INPUT_WW;
            }
            map_cols = 4;
          }
        }
      }
    }

    const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
    uint32_t max_variant_id_slen = *max_variant_id_slen_ptr;
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = bigstack_end_mark;
    uint16_t* cur_chr_codes = nullptr;
    uint32_t* cur_bps = nullptr;
    char** cur_ids = nullptr;
    double* cur_cms = nullptr;
    double cur_cm = 0.0;
    uint32_t at_least_one_nzero_cm = 0;
    uint32_t variant_ct = 0;
    uint32_t neg_bp_seen = 0;
    uintptr_t variant_skip_ct = 0;
    char* line_iter = line_start;
    while (1) {
      {
        // chrom, id, (cm?), pos
        char* chr_code_end = CurTokenEnd(line_iter);
        char* linebuf_iter = FirstNonTspace(chr_code_end);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto LoadMap_ret_MISSING_TOKENS;
        }
        uint32_t cur_chr_code;
        reterr = GetOrAddChrCodeDestructive(".map file", line_idx, prohibit_extra_chr, line_iter, chr_code_end, cip, &cur_chr_code);
        if (unlikely(reterr)) {
          goto LoadMap_ret_1;
        }
        if (!IsSet(cip->chr_mask, cur_chr_code)) {
          line_iter = linebuf_iter;
          goto LoadMap_skip_variant;
        }
        char* token_end = CurTokenEnd(linebuf_iter);
        uint32_t id_slen = token_end - linebuf_iter;
        if (id_slen > max_variant_id_slen) {
          max_variant_id_slen = id_slen;
        }
        tmp_alloc_end -= id_slen + 1;
        if (unlikely(tmp_alloc_end < tmp_alloc_base)) {
          goto LoadMap_ret_NOMEM;
        }
        memcpyx(tmp_alloc_end, linebuf_iter, id_slen, '\0');
        linebuf_iter = FirstNonTspace(token_end);
        if (unlikely(IsEolnKns(*linebuf_iter))) {
          goto LoadMap_ret_MISSING_TOKENS;
        }

        if (map_cols == 4) {
          char* cm_end = ScantokDouble(linebuf_iter, &cur_cm);
          if (unlikely(!cm_end)) {
            snprintf(g_logbuf, kLogbufSize, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, mapname);
            goto LoadMap_ret_MALFORMED_INPUT_WW;
          }
          at_least_one_nzero_cm = (cur_cm != 0.0);
          linebuf_iter = NextToken(cm_end);
          if (unlikely(!linebuf_iter)) {
            goto LoadMap_ret_MISSING_TOKENS;
          }
        }
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(linebuf_iter, &cur_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, mapname);
          goto LoadMap_ret_MALFORMED_INPUT_WW;
        }
        line_iter = linebuf_iter;
        if (cur_bp < 0) {
          neg_bp_seen = 1;
          goto LoadMap_skip_variant;
        }

        const uint32_t variant_idx_lowbits = variant_ct % kLoadMapBlockSize;
        if (!variant_idx_lowbits) {
          if (unlikely(S_CAST(uintptr_t, tmp_alloc_end - tmp_alloc_base) <= kLoadMapBlockSize * (sizeof(int16_t) + sizeof(int32_t) + sizeof(intptr_t) + sizeof(double)))) {
            goto LoadMap_ret_NOMEM;
          }
          cur_chr_codes = R_CAST(uint16_t*, tmp_alloc_base);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_chr_codes[kLoadMapBlockSize]));
          cur_bps = R_CAST(uint32_t*, tmp_alloc_base);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_bps[kLoadMapBlockSize]));
          cur_ids = R_CAST(char**, tmp_alloc_base);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_ids[kLoadMapBlockSize]));
          cur_cms = R_CAST(double*, tmp_alloc_base);
          tmp_alloc_base = R_CAST(unsigned char*, &(cur_cms[kLoadMapBlockSize]));
        }
        cur_chr_codes[variant_idx_lowbits] = cur_chr_code;
        cur_ids[variant_idx_lowbits] = R_CAST(char*, tmp_alloc_end);
        cur_cms[variant_idx_lowbits] = cur_cm;
        cur_bps[variant_idx_lowbits] = cur_bp;
        ++variant_ct;
      }
    LoadMap_skip_variant:
      ++variant_skip_ct;
      do {
        ++line_idx;
        line_iter = AdvPastDelim(line_iter, '\n');
        if (!TextGetUnsafe2(&map_txs, &line_iter)) {
          if (!TextStreamErrcode2(&map_txs, &reterr)) {
            goto LoadMap_eof;
          }
          goto LoadMap_ret_TSTREAM_FAIL;
        }
      } while (line_iter[0] == '#');
    }
  LoadMap_eof:
    if (unlikely(max_variant_id_slen > kMaxIdSlen)) {
      logerrputs("Error: Variant names are limited to " MAX_ID_SLEN_STR " characters.\n");
      goto LoadMap_ret_MALFORMED_INPUT;
    }

    if (unlikely(!variant_ct)) {
      char* write_iter = strcpya_k(g_logbuf, "Error: All ");
      write_iter = wtoa(variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (variant_skip_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " in .map file excluded by ");
      if (load_filter_log_import_flags) {
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
      }
      if (neg_bp_seen) {
        if (load_filter_log_import_flags) {
          write_iter = strcpya_k(write_iter, " and/or ");
        }
        write_iter = strcpya_k(write_iter, "negative bp coordinates");
      }
      strcpy_k(write_iter, ".\n");
      goto LoadMap_ret_INCONSISTENT_INPUT_WW;
    }
    // true requirement is weaker, but whatever
    g_bigstack_end = g_bigstack_base;
    g_bigstack_base = TextStreamMemStart(&map_txs);
    if (unlikely(CleanupTextStream2(mapname, &map_txs, &reterr))) {
      goto LoadMap_ret_1;
    }

    if (unlikely(bigstack_alloc_chridx(variant_ct, variant_chr_codes_ptr) ||
                 bigstack_alloc_u32(variant_ct, variant_bps_ptr) ||
                 bigstack_alloc_cp(variant_ct, variant_ids_ptr))) {
      goto LoadMap_ret_NOMEM;
    }
    ChrIdx* variant_chr_codes = *variant_chr_codes_ptr;
    uint32_t* variant_bps = *variant_bps_ptr;
    char** variant_ids = *variant_ids_ptr;
    double* variant_cms = nullptr;
    if (at_least_one_nzero_cm) {
      if (unlikely(bigstack_alloc_d(variant_ct, variant_cms_ptr))) {
        goto LoadMap_ret_NOMEM;
      }
      variant_cms = *variant_cms_ptr;
    } else {
      *variant_cms_ptr = nullptr;
    }
    *max_variant_id_slen_ptr = max_variant_id_slen;
    *variant_ct_ptr = variant_ct;
    const uint32_t full_block_ct = variant_ct / kLoadMapBlockSize;
    bigstack_mark = g_bigstack_base;
    unsigned char* read_iter = g_bigstack_end;  // bugfix (30 May 2018)
    BigstackEndSet(tmp_alloc_end);
    bigstack_end_mark = g_bigstack_end;

    for (uint32_t block_idx = 0; block_idx != full_block_ct; ++block_idx) {
      memcpy(&(variant_chr_codes[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(ChrIdx));
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(ChrIdx)]);
      memcpy(&(variant_bps[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(int32_t));
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int32_t)]);
      memcpy(&(variant_ids[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(intptr_t));
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(intptr_t)]);
      if (at_least_one_nzero_cm) {
        memcpy(&(variant_cms[block_idx * kLoadMapBlockSize]), read_iter, kLoadMapBlockSize * sizeof(double));
      }
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(double)]);
    }
    const uint32_t variant_ct_lowbits = variant_ct % kLoadMapBlockSize;
    memcpy(&(variant_chr_codes[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(ChrIdx));
    read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int16_t)]);
    memcpy(&(variant_bps[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(int32_t));
    read_iter = &(read_iter[kLoadMapBlockSize * sizeof(int32_t)]);
    memcpy(&(variant_ids[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(intptr_t));
    if (at_least_one_nzero_cm) {
      read_iter = &(read_iter[kLoadMapBlockSize * sizeof(intptr_t)]);
      memcpy(&(variant_cms[full_block_ct * kLoadMapBlockSize]), read_iter, variant_ct_lowbits * sizeof(double));
    }
  }
  while (0) {
  LoadMap_ret_TSTREAM_FAIL:
    TextStreamErrPrint(mapname, &map_txs);
    break;
  LoadMap_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  LoadMap_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, mapname);
  LoadMap_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  LoadMap_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  LoadMap_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
    reterr = kPglRetInconsistentInput;
    break;
  LoadMap_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 LoadMap_ret_1:
  CleanupTextStream2(mapname, &map_txs, &reterr);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

// Ok for in_psamname to alias outname.
PglErr RewritePsam(const char* in_psamname, const char* missing_catname, MiscFlags misc_flags, FamCol fam_cols, int32_t missing_pheno, uint32_t psam_01, uint32_t max_thread_ct, char* outname, char* outname_end, uint32_t* raw_sample_ctp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PhenoCol* pheno_cols = nullptr;
  char* pheno_names = nullptr;
  uint32_t pheno_ct = 0;
  PglErr reterr = kPglRetSuccess;
  {
    PedigreeIdInfo pii;
    InitPedigreeIdInfo(misc_flags, &pii);
    uint32_t raw_sample_ct = 0;
    uintptr_t* sample_include = nullptr;
    uintptr_t* sex_nm = nullptr;
    uintptr_t* sex_male = nullptr;
    uintptr_t* founder_info = nullptr;
    uintptr_t max_pheno_name_blen = 0;
    reterr = LoadPsam(in_psamname, nullptr, missing_catname, fam_cols, 0x7fffffff, missing_pheno, (misc_flags / kfMiscAffection01) & 1, (misc_flags / kfMiscNoCategorical) & 1, (misc_flags / kfMiscNeg9PhenoReallyMissing) & 1, max_thread_ct, &pii, &sample_include, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
    if (unlikely(reterr)) {
      goto RewritePsam_ret_1;
    }
    if (raw_sample_ctp) {
      *raw_sample_ctp = raw_sample_ct;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    // Note that --output-missing-phenotype doesn't apply to autoconversion.
    reterr = WritePsam(outname, sample_include, &(pii.sii), (&pii.parental_id_info), sex_nm, sex_male, pheno_cols, pheno_names, nullptr, "NA", raw_sample_ct, pheno_ct, max_pheno_name_blen, kfPsamColDefault, psam_01);
  }
 RewritePsam_ret_1:
  CleanupPhenoCols(pheno_ct, pheno_cols);
  free_cond(pheno_names);
  BigstackReset(bigstack_mark);
  return reterr;
}

BoolErr TpedToPgenSnp(uint32_t sample_idx_start, uint32_t sample_idx_stop, char allele1_char, char allele2_char, char input_missing_geno_char, char** text_iter_ptr, uintptr_t* genovec) {
  char* text_iter = *text_iter_ptr;
  for (uint32_t sample_idx = sample_idx_start; sample_idx != sample_idx_stop; ++sample_idx) {
    const char first_delim = text_iter[0];
    const char first_allele_char = text_iter[1];
    const char second_delim = text_iter[2];
    const char second_allele_char = text_iter[3];
    text_iter = &(text_iter[4]);
    if (unlikely(((first_delim != '\t') && (first_delim != ' ')) || ((second_delim != '\t') && (second_delim != ' ')))) {
      return 1;
    }
    uintptr_t cur_geno;
    if (first_allele_char == allele1_char) {
      if (second_allele_char == allele1_char) {
        continue;
      }
      if (unlikely(second_allele_char != allele2_char)) {
        return 1;
      }
      cur_geno = 1;
    } else if (first_allele_char == allele2_char) {
      if (second_allele_char == allele1_char) {
        cur_geno = 1;
      } else if (likely(second_allele_char == allele2_char)) {
        cur_geno = 2;
      } else {
        return 1;
      }
    } else {
      if (unlikely(((first_allele_char != '.') && (first_allele_char != input_missing_geno_char)) || ((second_allele_char != '.') && (second_allele_char != input_missing_geno_char)))) {
        return 1;
      }
      cur_geno = 3;
    }
    genovec[sample_idx / kBitsPerWordD2] |= cur_geno << ((sample_idx % kBitsPerWordD2) * 2);
  }
  *text_iter_ptr = text_iter;
  return 0;
}

// psam_generated assumed to be initialized to 1.
// Unlike plink 1.9, this does not support lines longer than 2 GiB.
// It's possible to parallelize this more, but it isn't realistically worth the
// effort, since there's so little reason to use this format over VCF for
// larger datasets.
PglErr TpedToPgen(const char* tpedname, const char* tfamname, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, FamCol fam_cols, int32_t missing_pheno, char input_missing_geno_char, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* psam_generated_ptr) {
  unsigned char* bigstack_mark = g_bigstack_base;
  TextStream tped_txs;
  PreinitTextStream(&tped_txs);
  uintptr_t line_idx = 0;

  char* pvar_cswritep = nullptr;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  STPgenWriter spgw;
  PglErr reterr = kPglRetSuccess;
  PreinitSpgw(&spgw);
  {
    uint32_t tfam_sample_ct = 0;
    if (import_flags & kfImportKeepAutoconv) {
      // Only need to generate a .psam if this is a conversion-only run, or
      // --keep-autoconv was specified.  Otherwise Plink2Core() can simply
      // interpret the .tfam as a psam file.
      reterr = RewritePsam(tfamname, missing_catname, misc_flags, fam_cols, missing_pheno, 0, max_thread_ct, outname, outname_end, &tfam_sample_ct);
      if (unlikely(reterr)) {
        goto TpedToPgen_ret_1;
      }
    } else {
      // Could scan the .tfam for the sake of performing a consistency check,
      // but that check will happen later anyway.
      *psam_generated_ptr = 0;
    }

    // First pass: determine variant_ct (applying chromosome filter),
    // sample_ct, and maximum line length; check whether CM column needs to be
    // in .pvar file.
    // (Tried making this single-pass, did not improve performance.)
    const uint32_t decompress_thread_ct = MAXV(1, max_thread_ct - 1);
    reterr = SizeAndInitTextStream(tpedname, bigstack_left(), decompress_thread_ct, &tped_txs);
    if (unlikely(reterr)) {
      goto TpedToPgen_ret_TSTREAM_FAIL;
    }
    ++line_idx;
    char* tped_line_start = TextGet(&tped_txs);
    if (unlikely(!tped_line_start)) {
      if (TextStreamErrcode2(&tped_txs, &reterr)) {
        goto TpedToPgen_ret_TSTREAM_FAIL;
      }
      snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", tpedname);
      goto TpedToPgen_ret_DEGENERATE_DATA;
    }
    uint32_t sample_ct;
    {
      const uint32_t token_ct = CountTokens(tped_line_start);
      if (unlikely(token_ct < 6)) {
        logerrputs("Error: Too few columns in .tped file.\n");
        goto TpedToPgen_ret_MALFORMED_INPUT;
      }
      sample_ct = (token_ct - 4) / 2;
      if (tfam_sample_ct) {
        if (unlikely(tfam_sample_ct != sample_ct)) {
          logerrprintfww("Error: .tped file has %u sample%s, while .tfam file has %u.\n", sample_ct, (sample_ct == 1)? "" : "s", tfam_sample_ct);
          goto TpedToPgen_ret_INCONSISTENT_INPUT;
        }
      } else {
        if (unlikely((token_ct % 2) == 1)) {
          snprintf(g_logbuf, kLogbufSize, "Error: Unexpected number of columns in .tped file (%u; even number expected).\n", token_ct);
          goto TpedToPgen_ret_MALFORMED_INPUT_WW;
        }
      }
      logprintf("--tped: %u sample%s present.\n", sample_ct, (sample_ct == 1)? "" : "s");
    }
    FinalizeChrset(load_filter_log_import_flags, cip);
    const uint32_t prohibit_extra_chr = (misc_flags / kfMiscProhibitExtraChr) & 1;
    uint32_t max_line_blen = TextLineEnd(&tped_txs) - tped_line_start;
    uint32_t variant_ct = 0;
    uint32_t neg_bp_seen = 0;
    uint32_t at_least_one_nzero_cm = 0;
    while (1) {
      char* chr_code_end = CurTokenEnd(tped_line_start);

      // must do this before chromosome-code null-termination
      char* cm_start = FirstNonTspace(FirstSpaceOrEoln(FirstNonTspace(chr_code_end)));
      char* bp_start = FirstNonTspace(FirstSpaceOrEoln(cm_start));
      if (unlikely(IsEolnKns(*bp_start))) {
        goto TpedToPgen_ret_MISSING_TOKENS;
      }

      uint32_t cur_chr_code;
      reterr = GetOrAddChrCodeDestructive(".tped file", line_idx, prohibit_extra_chr, tped_line_start, chr_code_end, cip, &cur_chr_code);
      if (unlikely(reterr)) {
        goto TpedToPgen_ret_1;
      }
      if (IsSet(cip->chr_mask, cur_chr_code)) {
        int32_t cur_bp;
        if (unlikely(ScanIntAbsDefcap(bp_start, &cur_bp))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, tpedname);
          goto TpedToPgen_ret_MALFORMED_INPUT_WW;
        }
        if (cur_bp >= 0) {
          if (variant_ct == kPglMaxVariantCt) {
            logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 3 variants.  We recommend using\nother software for very deep studies of small numbers of genomes.\n");
            goto TpedToPgen_ret_MALFORMED_INPUT;
          }
          ++variant_ct;

          if (!at_least_one_nzero_cm) {
            double cur_cm;
            if (unlikely(!ScantokDouble(cm_start, &cur_cm))) {
              snprintf(g_logbuf, kLogbufSize, "Error: Invalid centimorgan position on line %" PRIuPTR " of %s.\n", line_idx, tpedname);
              goto TpedToPgen_ret_MALFORMED_INPUT_WW;
            }
            at_least_one_nzero_cm = (cur_cm != 0.0);
          }
        } else {
          neg_bp_seen = 1;
        }
      }

      ++line_idx;
      tped_line_start = TextGet(&tped_txs);
      if (!tped_line_start) {
        break;
      }
      const uint32_t cur_line_blen = TextLineEnd(&tped_txs) - tped_line_start;
      if (max_line_blen < cur_line_blen) {
        max_line_blen = cur_line_blen;
      }
    }
    if (unlikely(TextStreamErrcode2(&tped_txs, &reterr))) {
      goto TpedToPgen_ret_TSTREAM_FAIL;
    }
    const uintptr_t variant_skip_ct = line_idx - 1 - variant_ct;
    if (unlikely(!variant_ct)) {
      if (!variant_skip_ct) {
        logerrputs("Error: No variants in --tped file.\n");
        goto TpedToPgen_ret_INCONSISTENT_INPUT;
      }
      char* write_iter = strcpya_k(g_logbuf, "Error: All ");
      write_iter = wtoa(variant_skip_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (variant_skip_ct != 1) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " in --tped file excluded by ");
      if (load_filter_log_import_flags) {
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
      }
      if (neg_bp_seen) {
        if (load_filter_log_import_flags) {
          write_iter = strcpya_k(write_iter, " and/or ");
        }
        write_iter = strcpya_k(write_iter, "negative bp coordinates");
      }
      strcpy_k(write_iter, ".\n");
      goto TpedToPgen_ret_INCONSISTENT_INPUT_WW;
    }
    if (unlikely(CleanupTextStream2(tpedname, &tped_txs, &reterr))) {
      goto TpedToPgen_ret_1;
    }
    BigstackReset(bigstack_mark);
    {
      char* write_iter = strcpya_k(g_logbuf, "--tped: ");
      write_iter = wtoa(line_idx - 1, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (line_idx != 2) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " scanned");
      if (variant_skip_ct) {
        write_iter = strcpya_k(write_iter, "; ");
        write_iter = wtoa(variant_skip_ct, write_iter);
        write_iter = strcpya_k(write_iter, " excluded by ");
        if (load_filter_log_import_flags) {
          AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        }
        if (neg_bp_seen) {
          if (load_filter_log_import_flags) {
            write_iter = strcpya_k(write_iter, " and/or ");
          }
          write_iter = strcpya_k(write_iter, "negative bp coordinates");
        }
        write_iter = strcpya_k(write_iter, ", ");
        write_iter = u32toa(variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " remaining");
      } else if (load_filter_log_import_flags) {
        write_iter = strcpya_k(write_iter, " (");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, " had no effect)");
      }
      strcpy_k(write_iter, ".\n");
      WordWrapB(0);
      logputsb();
    }
    reterr = InitTextStream(tpedname, MAXV(max_line_blen, kTextStreamBlenFast), decompress_thread_ct, &tped_txs);
    if (unlikely(reterr)) {
      goto TpedToPgen_ret_TSTREAM_REWIND_FAIL;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    const uint32_t output_zst = (import_flags / kfImportKeepAutoconvVzs) & 1;
    if (output_zst) {
      snprintf(&(outname_end[5]), kMaxOutfnameExtBlen - 5, ".zst");
    }
    const uintptr_t overflow_buf_size = kCompressStreamBlock + 64 + max_line_blen;
    reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, &pvar_css, &pvar_cswritep);
    if (unlikely(reterr)) {
      goto TpedToPgen_ret_1;
    }
    pvar_cswritep = strcpya_k(pvar_cswritep, "#CHROM\tPOS\tID\tREF\tALT");
    if (at_least_one_nzero_cm) {
      pvar_cswritep = strcpya_k(pvar_cswritep, "\tCM");
    }
    AppendBinaryEoln(&pvar_cswritep);

    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, kfPgenGlobal0, 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, outname, strerror(errno));
      }
      goto TpedToPgen_ret_1;
    }
    unsigned char* spgw_alloc;
    if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
      goto TpedToPgen_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    uintptr_t* genovec;
    if (unlikely(bigstack_alloc_w(sample_ctl2, &genovec))) {
      goto TpedToPgen_ret_NOMEM;
    }

    fputs("--tped: 0%", stdout);
    char* tped_line_iter = TextLineEnd(&tped_txs);
    uint32_t variant_idx = 0;
    uint32_t pct = 0;
    uint32_t next_print_idx = (variant_ct + 99) / 100;
    for (line_idx = 1; ; ++line_idx) {
      reterr = TextGetUnsafe(&tped_txs, &tped_line_iter);
      if (unlikely(reterr)) {
        goto TpedToPgen_ret_TSTREAM_REWIND_FAIL;
      }
      char* chr_code_end = CurTokenEnd(tped_line_iter);
      const uint32_t cur_chr_code = GetChrCodeCounted(cip, chr_code_end - tped_line_iter, tped_line_iter);
      if (!IsSet(cip->chr_mask, cur_chr_code)) {
        tped_line_iter = AdvPastDelim(chr_code_end, '\n');
        continue;
      }
      char* variant_id_start = FirstNonTspace(chr_code_end);
      char* variant_id_end = CurTokenEnd(variant_id_start);
      char* cm_start = FirstNonTspace(variant_id_end);
      char* cm_end = CurTokenEnd(cm_start);
      char* bp_start = FirstNonTspace(cm_end);
      // previously verified that bp token exists
      int32_t cur_bp;
      if (unlikely(ScanIntAbsDefcap(bp_start, &cur_bp))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, tpedname);
        goto TpedToPgen_ret_MALFORMED_INPUT_WW;
      }
      if (cur_bp < 0) {
        tped_line_iter = AdvPastDelim(bp_start, '\n');
        continue;
      }
      pvar_cswritep = chrtoa(cip, cur_chr_code, pvar_cswritep);
      *pvar_cswritep++ = '\t';
      pvar_cswritep = u32toa_x(cur_bp, '\t', pvar_cswritep);
      pvar_cswritep = memcpyax(pvar_cswritep, variant_id_start, variant_id_end - variant_id_start, '\t');
      tped_line_iter = FirstNonTspace(CurTokenEnd(bp_start));

      ZeroWArr(sample_ctl2, genovec);
      const char* allele1 = nullptr;
      const char* allele2 = nullptr;
      uint32_t allele1_slen = 0;
      uint32_t allele2_slen = 0;
      uint32_t sample_idx = 0;
      for (; sample_idx != sample_ct; ++sample_idx) {
        if (unlikely(IsEolnKns(*tped_line_iter))) {
          goto TpedToPgen_ret_MISSING_TOKENS;
        }
        char* first_allele_end = CurTokenEnd(tped_line_iter);
        const uint32_t first_allele_slen = first_allele_end - tped_line_iter;
        uint32_t missing_ct = 0;
        if ((first_allele_slen == allele1_slen) && memequal(tped_line_iter, allele1, allele1_slen)) {
          // do nothing
        } else if ((first_allele_slen == 1) && ((tped_line_iter[0] == '.') || (tped_line_iter[0] == input_missing_geno_char))) {
          missing_ct = 1;
        } else if (!allele1_slen) {
          allele1 = tped_line_iter;
          allele1_slen = first_allele_slen;
        } else {
          allele2 = tped_line_iter;
          allele2_slen = first_allele_slen;
          break;
        }
        char* second_allele_start = FirstNonTspace(first_allele_end);
        if (unlikely(IsEolnKns(*second_allele_start))) {
          goto TpedToPgen_ret_MISSING_TOKENS;
        }
        char* second_allele_end = CurTokenEnd(second_allele_start);
        const uint32_t second_allele_slen = second_allele_end - second_allele_start;
        if ((second_allele_slen == allele1_slen) && memequal(second_allele_start, allele1, allele1_slen)) {
        } else if ((second_allele_slen == 1) && ((second_allele_start[0] == '.') || (second_allele_start[0] == input_missing_geno_char))) {
          ++missing_ct;
        } else if (allele1_slen) {
          // only way allele1_slen can be zero here is if first allele was
          // missing, in which case we'll immediately exit on half-missing
          allele2 = second_allele_start;
          allele2_slen = second_allele_slen;
          break;
        }
        if (missing_ct) {
          if (unlikely(missing_ct == 1)) {
            goto TpedToPgen_ret_HALF_MISSING;
          }
          genovec[sample_idx / kBitsPerWordD2] |= (3 * k1LU) << ((sample_idx % kBitsPerWordD2) * 2);
        }
        tped_line_iter = FirstNonTspace(second_allele_end);
      }
      // At this point, both allele codes are locked in.
      if ((allele1_slen == 1) && (allele2_slen == 1)) {
        char* eoln_ptr = AdvToDelim(tped_line_iter, '\n');
        if (eoln_ptr[-1] == '\r') {
          // don't punish Windows eoln
          --eoln_ptr;
        }
        const uint32_t remaining_sample_ct = sample_ct - sample_idx;
        char* text_iter = &(tped_line_iter[-1]);
        if (eoln_ptr != &(text_iter[remaining_sample_ct * 4])) {
          goto TpedToPgen_general_case;
        }
        // Common-case optimization: if both allele codes are length-1, and the
        // line ends at the earliest possible byte, valid lines must strictly
        // alternate between delimiter-bytes and allele-code bytes.
        // We don't increment sample_idx or tped_line_iter until the loop
        // finishes, since if anything unexpected happens, we back up and
        // reparse the line the slow way to generate the correct error message.
        const char allele1_char = allele1[0];
        const char allele2_char = allele2[0];
        uint32_t sample_idx2 = sample_idx;
#ifdef __LP64__
        // Try to use vector instructions to process most of the rest of the
        // line, falling back on the scalar loop for the beginning and the end.
        //
        // We use movemask to efficiently convert between text-space and
        // genovec-space.  Its return type fits kInt16PerVec (16 for AVX2, 8
        // otherwise) 2-bit genotype values, so that's the natural chunk size.
        // (Yes, this is overkill for the immediate problem, but it's a
        // reasonable testbed for parsing ideas.)
        const uint32_t sample_idx_interior_start = RoundUpPow2(sample_idx, kInt16PerVec);
        const uint32_t sample_idx_interior_stop = RoundDownPow2(sample_ct, kInt16PerVec);
        if (sample_idx_interior_stop > sample_idx_interior_start) {
          if (TpedToPgenSnp(sample_idx2, sample_idx_interior_start, allele1_char, allele2_char, input_missing_geno_char, &text_iter, genovec)) {
            goto TpedToPgen_general_case;
          }
          const uint32_t vidx_stop = sample_idx_interior_stop / kInt16PerVec;

          const VecUc vvec_tab = vecuc_set1('\t');
          const VecUc vvec_space = vecuc_set1(' ');
          const VecUc vvec_ref = vecuc_set1(allele1_char);
          const VecUc vvec_alt = vecuc_set1(allele2_char);
          const VecUc vvec_dot = vecuc_set1('.');
          const VecUc vvec_missing = vecuc_set1(input_missing_geno_char);
          const VecUc m8 = vecuc_set1_epi16(0x00FF);
          // We bitwise-and all validity checks with ok_acc, and then perform
          // a single check of ok_acc's contents at loop-end.
          VecUc ok_acc = vecuc_set1(-1);

          Vec8thUint* genovec_alias = R_CAST(Vec8thUint*, genovec);
          for (uint32_t vidx = sample_idx_interior_start / kInt16PerVec; vidx != vidx_stop; ++vidx) {
            const VecUc text0 = vecuc_loadu(text_iter);
            const VecUc text1 = vecuc_loadu(&(text_iter[kBytesPerVec]));
            text_iter = &(text_iter[2 * kBytesPerVec]);

            const VecUc even_all = vecuc_gather_even(text0, text1, m8);
            ok_acc = ok_acc & ((even_all == vvec_tab) | (even_all == vvec_space));
            const VecUc odd_all = vecuc_gather_odd(text0, text1);

            const VecUc cur_ref = (odd_all == vvec_ref);
            const VecUc cur_alt = (odd_all == vvec_alt);
            const VecUc cur_missing = (odd_all == vvec_dot) | (odd_all == vvec_missing);
            // All allele codes must match REF, ALT, or missing.
            ok_acc = ok_acc & (cur_ref | cur_alt | cur_missing);

            // Half-missing check could be done in vector-space instead (so we
            // only have one movemask operation instead of two per loop, and we
            // remove an if-statement), but that code is more complicated, and
            // doesn't seem to be any faster.

            // todo: better ARM implementation
            const Vec8thUint alt_bits = vecuc_movemask(cur_alt);
            const Vec8thUint missing_bits = vecuc_movemask(cur_missing);
            // Even missing bytes must match odd missing bits.
            if ((missing_bits ^ (missing_bits >> 1)) & (kVec8thUintMax / 3)) {
              // Don't jump directly to HALF_MISSING, since there may be
              // another earlier error.
              goto TpedToPgen_general_case;
            }
            genovec_alias[vidx] = ((alt_bits & (kVec8thUintMax / 3)) + ((alt_bits >> 1) & (kVec8thUintMax / 3))) | missing_bits;
          }
          if (!vec0255_is_all_set(ok_acc)) {
            goto TpedToPgen_general_case;
          }
          sample_idx2 = sample_idx_interior_stop;
        }
#endif
        if (unlikely(TpedToPgenSnp(sample_idx2, sample_ct, allele1_char, allele2_char, input_missing_geno_char, &text_iter, genovec))) {
          goto TpedToPgen_general_case;
        }
        tped_line_iter = eoln_ptr;
      } else {
      TpedToPgen_general_case:
        for (; sample_idx != sample_ct; ++sample_idx) {
          if (unlikely(IsEolnKns(*tped_line_iter))) {
            goto TpedToPgen_ret_MISSING_TOKENS;
          }
          char* first_allele_end = CurTokenEnd(tped_line_iter);
          const uint32_t first_allele_slen = first_allele_end - tped_line_iter;
          uintptr_t allele2_ct = 0;
          if ((first_allele_slen == allele1_slen) && memequal(tped_line_iter, allele1, allele1_slen)) {
            // do nothing
          } else if ((first_allele_slen == allele2_slen) && memequal(tped_line_iter, allele2, allele2_slen)) {
            allele2_ct = 1;
          } else if (likely((first_allele_slen == 1) && ((tped_line_iter[0] == '.') || (tped_line_iter[0] == input_missing_geno_char)))) {
            allele2_ct = 4;
          } else {
            goto TpedToPgen_ret_MULTIALLELIC;
          }
          char* second_allele_start = FirstNonTspace(first_allele_end);
          if (unlikely(IsEolnKns(*second_allele_start))) {
            goto TpedToPgen_ret_MISSING_TOKENS;
          }
          char* second_allele_end = CurTokenEnd(second_allele_start);
          const uint32_t second_allele_slen = second_allele_end - second_allele_start;
          if ((second_allele_slen == allele1_slen) && memequal(second_allele_start, allele1, allele1_slen)) {
            // do nothing
          } else if ((second_allele_slen == allele2_slen) && memequal(second_allele_start, allele2, allele2_slen)) {
            ++allele2_ct;
          } else if (likely((second_allele_slen == 1) && ((second_allele_start[0] == '.') || (second_allele_start[0] == input_missing_geno_char)))) {
            allele2_ct += 4;
          } else {
            goto TpedToPgen_ret_MULTIALLELIC;
          }
          if (allele2_ct) {
            if (allele2_ct >= 4) {
              if (unlikely(allele2_ct != 8)) {
                // Value of 4 or 5 corresponds to 0/. or 1/.
                goto TpedToPgen_ret_HALF_MISSING;
              }
              // Turn this into a genotype code.
              allele2_ct = 3;
            }
            genovec[sample_idx / kBitsPerWordD2] |= allele2_ct << ((sample_idx % kBitsPerWordD2) * 2);
          }

          tped_line_iter = FirstNonTspace(second_allele_end);
        }
      }
      // Count alleles, and swap allele2 to REF if it's more common than
      // allele1.
      // (This could be done slightly more efficiently with e.g. a function
      // that only tracked 0b00 and 0b10 counts, but it isn't a bottleneck.)
      STD_ARRAY_DECL(uint32_t, 4, genocounts);
      GenoarrCountFreqsUnsafe(genovec, sample_ct, genocounts);
      if (genocounts[2] > genocounts[0]) {
        const char* tmp_allele = allele1;
        allele1 = allele2;
        allele2 = tmp_allele;

        const uint32_t tmp_slen = allele1_slen;
        allele1_slen = allele2_slen;
        allele2_slen = tmp_slen;

        GenovecInvertUnsafe(sample_ct, genovec);
        ZeroTrailingNyps(sample_ct, genovec);
      }

      if (allele1) {
        pvar_cswritep = memcpya(pvar_cswritep, allele1, allele1_slen);
      } else {
        *pvar_cswritep++ = '.';
      }
      *pvar_cswritep++ = '\t';
      if (allele2) {
        pvar_cswritep = memcpya(pvar_cswritep, allele2, allele2_slen);
      } else {
        *pvar_cswritep++ = '.';
      }

      if (at_least_one_nzero_cm) {
        double cur_cm;
        ScanadvDouble(cm_start, &cur_cm);
        *pvar_cswritep++ = '\t';
        pvar_cswritep = dtoa_g_p8(cur_cm, pvar_cswritep);
      }
      AppendBinaryEoln(&pvar_cswritep);
      if (unlikely(Cswrite(&pvar_css, &pvar_cswritep))) {
        goto TpedToPgen_ret_WRITE_FAIL;
      }

      if (unlikely(SpgwAppendBiallelicGenovec(genovec, &spgw))) {
        goto TpedToPgen_ret_WRITE_FAIL;
      }

      ++variant_idx;
      if (variant_idx >= next_print_idx) {
        if (variant_idx == variant_ct) {
          break;
        }
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        pct = (variant_idx * 100LLU) / variant_ct;
        printf("\b\b%u%%", pct++);
        fflush(stdout);
        next_print_idx = (pct * S_CAST(uint64_t, variant_ct) + 99) / 100;
      }
      tped_line_iter = AdvPastDelim(tped_line_iter, '\n');
    }
    reterr = SpgwFinish(&spgw);
    if (unlikely(reterr)) {
      goto TpedToPgen_ret_1;
    }
    putc_unlocked('\r', stdout);
    char* write_iter = strcpya_k(g_logbuf, "--tped: ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen + 5);
    write_iter = strcpya_k(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    if (output_zst) {
      write_iter = strcpya_k(write_iter, ".zst");
    }
    if (import_flags & kfImportKeepAutoconv) {
      write_iter = strcpya_k(write_iter, " + ");
      write_iter = memcpya(write_iter, outname, outname_base_slen);
      write_iter = strcpya_k(write_iter, ".psam");
    }
    strcpy_k(write_iter, " written.\n");
    WordWrapB(0);
    logputsb();
  }
  while (0) {
  TpedToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  TpedToPgen_ret_TSTREAM_FAIL:
    TextStreamErrPrint(tpedname, &tped_txs);
    break;
  TpedToPgen_ret_TSTREAM_REWIND_FAIL:
    if ((reterr == kPglRetOpenFail) || (reterr == kPglRetEof)) {
      reterr = kPglRetRewindFail;
    }
    TextStreamErrPrintRewind(tpedname, &tped_txs, &reterr);
    break;
  TpedToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  TpedToPgen_ret_HALF_MISSING:
    putc_unlocked('\n', stdout);
    logerrprintfww("Error: Half-missing genotype on line %" PRIuPTR " of %s.\n", line_idx, tpedname);
    reterr = kPglRetMalformedInput;
    break;
  TpedToPgen_ret_MULTIALLELIC:
    putc_unlocked('\n', stdout);
    logerrprintfww("Error: Multiallelic variant on line %" PRIuPTR " of %s. This violates the .tped specification; please reformat this as e.g. VCF.\n", line_idx, tpedname);
    reterr = kPglRetMalformedInput;
    break;
  TpedToPgen_ret_MISSING_TOKENS:
    putc_unlocked('\n', stdout);
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, tpedname);
  TpedToPgen_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  TpedToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  TpedToPgen_ret_INCONSISTENT_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  TpedToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  TpedToPgen_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 TpedToPgen_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CswriteCloseCond(&pvar_css, pvar_cswritep);
  CleanupTextStream2(tpedname, &tped_txs, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct Plink1SmajTransposeCtxStruct {
  uint32_t sample_ct;
  uint32_t loadbuf_ul_stride;

  uintptr_t* plink1_smaj_loadbuf_iter;

  const uintptr_t* allele_flips_iter;

  VecW** thread_vecaligned_bufs;
  uintptr_t** thread_write_genovecs;
  PgenWriterCommon** pwcs;

  uint32_t cur_block_write_ct;
} Plink1SmajTransposeCtx;

void Plink1SmajTransposeMain(uint32_t tidx, Plink1SmajTransposeCtx* ctx) {
  const uint32_t sample_ct = ctx->sample_ct;
  const uint32_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
  const uintptr_t transpose_block_ct_m1 = (sample_ct - 1) / kPglNypTransposeBatch;
  PgenWriterCommon* pwcp = ctx->pwcs[tidx];
  VecW* vecaligned_buf = ctx->thread_vecaligned_bufs[tidx];
  uintptr_t* write_genovec = ctx->thread_write_genovecs[tidx];
  const uintptr_t cur_block_write_ct = ctx->cur_block_write_ct;
  const uintptr_t loadbuf_ul_stride = ctx->loadbuf_ul_stride;
  uint32_t write_idx = tidx * kPglVblockSize;
  uintptr_t* read_iter = &(ctx->plink1_smaj_loadbuf_iter[write_idx / kBitsPerWordD2]);
  const uintptr_t* allele_flips_iter = ctx->allele_flips_iter;
  if (allele_flips_iter) {
    allele_flips_iter = &(allele_flips_iter[write_idx / kBitsPerWord]);
  }
  const uint32_t write_idx_end = MINV(write_idx + kPglVblockSize, cur_block_write_ct);
  while (write_idx < write_idx_end) {
    const uintptr_t* read_iter2 = read_iter;
    // uintptr_t* write_iter = write_genovec;
    const uint32_t vblock_size = MINV(kPglNypTransposeBatch, write_idx_end - write_idx);
    uint32_t read_batch_size = kPglNypTransposeBatch;
    for (uintptr_t transpose_block_idx = 0; ; ++transpose_block_idx) {
      if (transpose_block_idx >= transpose_block_ct_m1) {
        if (transpose_block_idx > transpose_block_ct_m1) {
          break;
        }
        read_batch_size = ModNz(sample_ct, kPglNypTransposeBatch);
      }
      TransposeNypblock(read_iter2, loadbuf_ul_stride, sample_ctaw2, read_batch_size, vblock_size, &(write_genovec[transpose_block_idx * kPglNypTransposeWords]), vecaligned_buf);
      read_iter2 = &(read_iter2[kPglNypTransposeBatch * loadbuf_ul_stride]);
    }
    uintptr_t* cur_write_genovec = write_genovec;
    if (!allele_flips_iter) {
      for (uint32_t uii = 0; uii != vblock_size; ++uii, cur_write_genovec = &(cur_write_genovec[sample_ctaw2])) {
        PgrPlink1ToPlink2InplaceUnsafe(sample_ct, cur_write_genovec);
        ZeroTrailingNyps(sample_ct, cur_write_genovec);
        PwcAppendBiallelicGenovec(cur_write_genovec, pwcp);
      }
    } else {
      for (uint32_t uii = 0; uii != vblock_size; ++uii, cur_write_genovec = &(cur_write_genovec[sample_ctaw2])) {
        PgrPlink1ToPlink2InplaceUnsafe(sample_ct, cur_write_genovec);
        if (IsSet(allele_flips_iter, uii)) {
          // could have a dedicated flip-and-invert function
          GenovecInvertUnsafe(sample_ct, cur_write_genovec);
        }
        ZeroTrailingNyps(sample_ct, cur_write_genovec);
        PwcAppendBiallelicGenovec(cur_write_genovec, pwcp);
      }
    }
    write_idx += vblock_size;
    read_iter = &(read_iter[kPglNypTransposeWords]);
    if (allele_flips_iter) {
      allele_flips_iter = &(allele_flips_iter[kPglNypTransposeWords / 2]);
    }
  }
}

THREAD_FUNC_DECL Plink1SmajTransposeThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  const uint32_t tidx = arg->tidx;
  Plink1SmajTransposeCtx* ctx = S_CAST(Plink1SmajTransposeCtx*, arg->sharedp->context);
  do {
    Plink1SmajTransposeMain(tidx, ctx);
  } while (!THREAD_BLOCK_FINISH(arg));
  THREAD_RETURN;
}

static_assert((kPglNypTransposeBatch % kNypsPerVec == 0) && (kPglNypTransposeBatch % kVecsPerCacheline == 0), "Plink1SampleMajorToPgenLowmem() needs to be updated.");
PglErr Plink1SampleMajorToPgenLowmem(const char* pgenname, const uintptr_t* allele_flips, uintptr_t variant_ct, uintptr_t sample_ct, uint32_t real_ref_alleles, uint32_t raw_load_batch_size, uint32_t raw_load_batch_ct, FILE* infile, unsigned char* raw_loadbuf) {
  // caller expected to free memory
  PglErr reterr = kPglRetSuccess;
  STPgenWriter spgw;
  PreinitSpgw(&spgw);
  {
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(pgenname, nullptr, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, kfPgenGlobal0, 2 - real_ref_alleles, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (unlikely(reterr)) {
      if (reterr == kPglRetOpenFail) {
        logerrprintfww(kErrprintfFopen, pgenname, strerror(errno));
      }
      goto Plink1SampleMajorToPgenLowmem_ret_1;
    }
    const uintptr_t sample_ctaw2 = NypCtToAlignedWordCt(sample_ct);
    unsigned char* spgw_alloc;
    VecW* vecaligned_buf;
    uintptr_t* write_genovec;
    if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc) ||
                 bigstack_alloc_v(kPglNypTransposeBufbytes / kBytesPerVec, &vecaligned_buf) ||
                 bigstack_alloc_w(sample_ctaw2 * kPglNypTransposeBatch, &write_genovec))) {
      goto Plink1SampleMajorToPgenLowmem_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

    uintptr_t cachelines_avail = bigstack_left() / kCacheline;
    const uint64_t full_load_vecs_req = sample_ct * S_CAST(uint64_t, NypCtToVecCt(variant_ct));
    // Tried making this a double-buffer so that load and transpose-compress
    // could occur simultaneously; that doesn't seem any better than loading as
    // much as much as possible per pass.
    uintptr_t* plink1_smaj_loadbuf = R_CAST(uintptr_t*, g_bigstack_base);
    uint32_t cur_vidx_ct;
    if (full_load_vecs_req > cachelines_avail * kVecsPerCacheline) {
      // Load the largest multiple of kPglNypTransposeBatch variants at a time
      // that fits into the remaining workspace.
      const uint64_t min_load_cl_req = DivUpU64(sample_ct * NypCtToVecCt(kPglNypTransposeBatch), kVecsPerCacheline);
      const uint32_t vbatches_per_load = cachelines_avail / min_load_cl_req;
      if (unlikely(!vbatches_per_load)) {
        g_failed_alloc_attempt_size = min_load_cl_req * kCacheline;
        goto Plink1SampleMajorToPgenLowmem_ret_NOMEM;
      }
      cur_vidx_ct = vbatches_per_load * kPglNypTransposeBatch;
    } else {
      cur_vidx_ct = variant_ct;
    }
    // can't have any allocations past this point unless we update
    // g_bigstack_base

    uint32_t cur_vidx_ct4 = NypCtToByteCt(cur_vidx_ct);
    uintptr_t cur_vidx_ctaw2 = NypCtToAlignedWordCt(cur_vidx_ct);
    const uint32_t variant_ct4 = NypCtToByteCt(variant_ct);
    const uint32_t raw_load_batch_ct_m1 = raw_load_batch_ct - 1;
    const uintptr_t transpose_block_ct_m1 = (sample_ct - 1) / kPglNypTransposeBatch;
    const uint32_t pass_ct = 1 + (variant_ct - 1) / cur_vidx_ct;
    uint32_t pass_idx1 = 0;
    for (uint32_t cur_vidx_base = 0; ; ) {
      uint32_t cur_raw_load_batch_size = raw_load_batch_size;
      uintptr_t* smaj_loadbuf_iter = plink1_smaj_loadbuf;
      putc_unlocked('\r', stdout);
      ++pass_idx1;
      printf("Pass %u/%u: loading... 0%%", pass_idx1, pass_ct);
      fflush(stdout);
      uint32_t pct = 0;
      uint32_t next_print_idx = (raw_load_batch_ct + 99) / 100;
      const uint64_t seek_addl_offset = 3 + cur_vidx_base / 4;
      for (uint32_t raw_load_batch_idx = 0; ; ) {
        // possible todo: check if multithreaded reading is faster if we know
        // we're reading from a SSD
        if (raw_load_batch_size == 1) {
          if (unlikely(fseeko(infile, seek_addl_offset + raw_load_batch_idx * S_CAST(uint64_t, variant_ct4), SEEK_SET) ||
                       (!fread_unlocked(smaj_loadbuf_iter, cur_vidx_ct4, 1, infile)))) {
            goto Plink1SampleMajorToPgenLowmem_ret_READ_FAIL;
          }
          smaj_loadbuf_iter = &(smaj_loadbuf_iter[cur_vidx_ctaw2]);
        } else {
          if (unlikely(!fread_unlocked(raw_loadbuf, cur_raw_load_batch_size * variant_ct4, 1, infile))) {
            goto Plink1SampleMajorToPgenLowmem_ret_READ_FAIL;
          }
          unsigned char* raw_loadbuf_iter = &(raw_loadbuf[cur_vidx_base / 4]);
          for (uint32_t uii = 0; uii != cur_raw_load_batch_size; ++uii) {
            memcpy(smaj_loadbuf_iter, raw_loadbuf_iter, cur_vidx_ct4);
            raw_loadbuf_iter = &(raw_loadbuf_iter[variant_ct4]);
            smaj_loadbuf_iter = &(smaj_loadbuf_iter[cur_vidx_ctaw2]);
          }
        }
        ++raw_load_batch_idx;
        if (raw_load_batch_idx >= raw_load_batch_ct_m1) {
          if (raw_load_batch_idx > raw_load_batch_ct_m1) {
            break;
          }
          cur_raw_load_batch_size = sample_ct - raw_load_batch_idx * raw_load_batch_size;
        }
        if (raw_load_batch_idx >= next_print_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (raw_load_batch_idx * 100LLU) / raw_load_batch_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_idx = (pct * S_CAST(uint64_t, raw_load_batch_ct) + 99) / 100;
        }
      }
      putc_unlocked('\r', stdout);
      printf("Pass %u/%u: transposing and compressing... 0%%", pass_idx1, pass_ct);
      fflush(stdout);
      const uintptr_t* allele_flips_iter = allele_flips? &(allele_flips[cur_vidx_base / kBitsPerWord]) : nullptr;
      uintptr_t* read_iter = plink1_smaj_loadbuf;
      pct = 0;
      next_print_idx = (cur_vidx_ct + 99) / 100;
      for (uint32_t write_idx = 0; write_idx < cur_vidx_ct; ) {
        // loadbuf_ul_stride = cur_vidx_ctaw2
        const uintptr_t* read_iter2 = read_iter;
        const uint32_t vblock_size = MINV(kPglNypTransposeBatch, cur_vidx_ct - write_idx);
        uint32_t read_batch_size = kPglNypTransposeBatch;
        for (uintptr_t transpose_block_idx = 0; ; ++transpose_block_idx) {
          if (transpose_block_idx >= transpose_block_ct_m1) {
            if (transpose_block_idx > transpose_block_ct_m1) {
              break;
            }
            read_batch_size = ModNz(sample_ct, kPglNypTransposeBatch);
          }
          TransposeNypblock(read_iter2, cur_vidx_ctaw2, sample_ctaw2, read_batch_size, vblock_size, &(write_genovec[transpose_block_idx * kPglNypTransposeWords]), vecaligned_buf);
          read_iter2 = &(read_iter2[kPglNypTransposeBatch * cur_vidx_ctaw2]);
        }
        uintptr_t* cur_write_genovec = write_genovec;
        for (uint32_t uii = 0; uii != vblock_size; ++uii, cur_write_genovec = &(cur_write_genovec[sample_ctaw2])) {
          PgrPlink1ToPlink2InplaceUnsafe(sample_ct, cur_write_genovec);
          if (allele_flips_iter && IsSet(allele_flips_iter, uii)) {
            GenovecInvertUnsafe(sample_ct, cur_write_genovec);
          }
          ZeroTrailingNyps(sample_ct, cur_write_genovec);
          reterr = SpgwAppendBiallelicGenovec(cur_write_genovec, &spgw);
          if (unlikely(reterr)) {
            goto Plink1SampleMajorToPgenLowmem_ret_WRITE_FAIL;
          }
        }
        if (write_idx >= next_print_idx) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          pct = (write_idx * 100LLU) / cur_vidx_ct;
          printf("\b\b%u%%", pct++);
          fflush(stdout);
          next_print_idx = (pct * S_CAST(uint64_t, cur_vidx_ct) + 99) / 100;
        }
        write_idx += vblock_size;
        read_iter = &(read_iter[kPglNypTransposeWords]);
        if (allele_flips_iter) {
          allele_flips_iter = &(allele_flips_iter[kPglNypTransposeWords / 2]);
        }
      }
      cur_vidx_base += cur_vidx_ct;
      if (cur_vidx_base == variant_ct) {
        if (pct > 10) {
          putc_unlocked('\b', stdout);
        }
        fputs("\b\bdone.\n", stdout);
        break;
      }
      fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                     ", stdout);
      if (unlikely(fseeko(infile, 3, SEEK_SET))) {
        goto Plink1SampleMajorToPgenLowmem_ret_READ_FAIL;
      }
      if (variant_ct - cur_vidx_base <= cur_vidx_ct) {
        cur_vidx_ct = variant_ct - cur_vidx_base;
        cur_vidx_ct4 = NypCtToByteCt(cur_vidx_ct);
        cur_vidx_ctaw2 = NypCtToAlignedWordCt(cur_vidx_ct);
      }
    }
    reterr = SpgwFinish(&spgw);
    if (unlikely(reterr)) {
      goto Plink1SampleMajorToPgenLowmem_ret_1;
    }
    logprintf("Transpose complete.\n");
  }
  while (0) {
  Plink1SampleMajorToPgenLowmem_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Plink1SampleMajorToPgenLowmem_ret_READ_FAIL:
    if (feof_unlocked(infile)) {
      errno = 0;
    }
    logputs("\n");
    logerrprintfww(kErrprintfFread, ".bed file", rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  Plink1SampleMajorToPgenLowmem_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  }
 Plink1SampleMajorToPgenLowmem_ret_1:
  CleanupSpgw(&spgw, &reterr);
  return reterr;
}

PglErr Plink1SampleMajorToPgen(const char* pgenname, const uintptr_t* allele_flips, uintptr_t variant_ct, uintptr_t sample_ct, uint32_t real_ref_alleles, uint32_t max_thread_ct, FILE* infile) {
  unsigned char* bigstack_mark = g_bigstack_base;
  MTPgenWriter* mpgwp = nullptr;
  ThreadGroup tg;
  PreinitThreads(&tg);
  Plink1SmajTransposeCtx ctx;
  PglErr reterr = kPglRetSuccess;
  {
    // If allele_flips is nullptr, we are being called by Plink2Core().
    // Otherwise, we are being called by PedmapToPgen().
    if (!allele_flips) {
      // file size already validated by PgfiInitPhase1()
      logprintfww("Sample-major .bed file detected.  Transposing to %s .\n", pgenname);
    } else {
      logprintfww("Transposing sample-major .bed to %s , and setting major alleles to provisional-REF.\n", pgenname);
    }
    if (unlikely((!variant_ct) || (!sample_ct))) {
      logputs("\n");
      logerrputs("Error: Zero-variant/zero-sample .pgen writing is not supported.\n");
      goto Plink1SampleMajorToPgen_ret_INCONSISTENT_INPUT;
    }
    const uint32_t variant_ct4 = NypCtToByteCt(variant_ct);
    unsigned char* raw_loadbuf = nullptr;
    uint32_t raw_load_batch_size = 1;
#ifndef __APPLE__
    if (variant_ct4 < 5120) {
      // assuming 4K block size, fseek won't let us avoid reading many
      // unnecessary disk blocks
      raw_load_batch_size += 131071 / variant_ct4;
      if (unlikely(bigstack_alloc_uc(raw_load_batch_size * variant_ct4, &raw_loadbuf))) {
        goto Plink1SampleMajorToPgen_ret_NOMEM;
      }
      // bugfix (2 Jun 2022): forgot to skip header bytes
      if (unlikely(fseeko(infile, 3, SEEK_SET))) {
        goto Plink1SampleMajorToPgen_ret_READ_FAIL;
      }
    }
#else
    // macOS seems to suck at seek/read interleaving
    // haven't carefully tuned these constants, but I know they shouldn't be
    // much smaller on my test Mac.
    if (variant_ct4 < 1048576) {
      raw_load_batch_size += 2097151 / variant_ct4;
      if (unlikely(bigstack_alloc_uc(raw_load_batch_size * variant_ct4, &raw_loadbuf))) {
        goto Plink1SampleMajorToPgen_ret_NOMEM;
      }
      if (unlikely(fseeko(infile, 3, SEEK_SET))) {
        goto Plink1SampleMajorToPgen_ret_READ_FAIL;
      }
    }
#endif
    const uint32_t raw_load_batch_ct_m1 = (sample_ct - 1) / raw_load_batch_size;
    if (!raw_load_batch_ct_m1) {
      raw_load_batch_size = sample_ct;
    }
    const uint32_t raw_load_batch_ct = raw_load_batch_ct_m1 + 1;
    unsigned char* bigstack_mark2 = g_bigstack_base;
    if ((max_thread_ct > 1) && (variant_ct > kPglVblockSize * 2)) {
      uintptr_t alloc_base_cacheline_ct;
      uint64_t mpgw_per_thread_cacheline_ct;
      uint32_t vrec_len_byte_ct;
      uint64_t vblock_cacheline_ct;
      MpgwInitPhase1(nullptr, variant_ct, sample_ct, kfPgenGlobal0, &alloc_base_cacheline_ct, &mpgw_per_thread_cacheline_ct, &vrec_len_byte_ct, &vblock_cacheline_ct);
#ifndef __LP64__
      if ((mpgw_per_thread_cacheline_ct > (0x7fffffff / kCacheline)) || (vblock_cacheline_ct > (0x7fffffff / kCacheline))) {
        goto Plink1SampleMajorToPgen_fallback;
      }
#endif

      uint32_t calc_thread_ct = DivUp(variant_ct, kPglVblockSize);
      if (calc_thread_ct >= max_thread_ct) {
        calc_thread_ct = (max_thread_ct > 2)? (max_thread_ct - 1) : max_thread_ct;
      }
      // note that BIGSTACK_ALLOC_X() doesn't work here due to variable-size
      // array at end
      mpgwp = S_CAST(MTPgenWriter*, bigstack_alloc((calc_thread_ct + DivUp(sizeof(MTPgenWriter), kBytesPerWord)) * sizeof(intptr_t)));
      if (!mpgwp) {
        goto Plink1SampleMajorToPgen_fallback;
      }
      PreinitMpgw(mpgwp);
      if (bigstack_alloc_vp(calc_thread_ct, &ctx.thread_vecaligned_bufs) ||
          bigstack_alloc_wp(calc_thread_ct, &ctx.thread_write_genovecs)) {
        goto Plink1SampleMajorToPgen_fallback;
      }
      ctx.pwcs = &(mpgwp->pwcs[0]);
      uintptr_t cachelines_avail = bigstack_left() / kCacheline;
      // inner loop transposes kPglNypTransposeBatch variants at a time
      const uintptr_t transpose_thread_cacheline_ct = kPglNypTransposeBufbytes / kCacheline + NypCtToVecCt(sample_ct) * (kPglNypTransposeBatch / kVecsPerCacheline);
      if (cachelines_avail < calc_thread_ct * S_CAST(uint64_t, transpose_thread_cacheline_ct)) {
        goto Plink1SampleMajorToPgen_fallback;
      }
      for (uint32_t tidx = 0; tidx != calc_thread_ct; ++tidx) {
        ctx.thread_vecaligned_bufs[tidx] = S_CAST(VecW*, bigstack_alloc_raw(kPglNypTransposeBufbytes));
        ctx.thread_write_genovecs[tidx] = S_CAST(uintptr_t*, bigstack_alloc_raw(NypCtToVecCt(sample_ct) * kBytesPerVec * kPglNypTransposeBatch));
      }
      cachelines_avail = bigstack_left() / kCacheline;
      // Main workflow:
      // 1. Load next calc_thread_ct * vblock_group_ct * kPglVblockSize
      //    variants.
      //    calc_thread_ct is reduced as necessary to ensure the compression
      //    write buffers use <= 1/8 of total workspace.
      //    with calc_thread_ct determined, vblock_group_ct is then chosen to
      //    use as much of the remaining workspace as possible.
      // 2. Repeat vblock_group_ct times:
      //    a. Spawn threads processing calc_thread_ct vblocks
      //    b. Join threads
      //    c. Flush results
      // 3. Goto step 1 unless eof.  (vblock_group_ct may be smaller on last
      //    iteration.)
      // No double-buffering here since main bottleneck is how many variants we
      // can load at once.
      if ((cachelines_avail / 8) < alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) {
        if ((cachelines_avail / 8) >= alloc_base_cacheline_ct + 2 * mpgw_per_thread_cacheline_ct) {
          calc_thread_ct = ((cachelines_avail / 8) - alloc_base_cacheline_ct) / mpgw_per_thread_cacheline_ct;
        } else {
          goto Plink1SampleMajorToPgen_fallback;
        }
      }
      // todo: determine appropriate calc_thread_ct limit.  (should not be less
      // than 7-8.)
      unsigned char* mpgw_alloc = S_CAST(unsigned char*, bigstack_alloc_raw((alloc_base_cacheline_ct + mpgw_per_thread_cacheline_ct * calc_thread_ct) * kCacheline));
      reterr = MpgwInitPhase2(pgenname, nullptr, variant_ct, sample_ct, kPgenWriteBackwardSeek, kfPgenGlobal0, 2 - real_ref_alleles, vrec_len_byte_ct, vblock_cacheline_ct, calc_thread_ct, mpgw_alloc, mpgwp);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          logerrprintfww(kErrprintfFopen, pgenname, strerror(errno));
        }
        goto Plink1SampleMajorToPgen_ret_1;
      }

      cachelines_avail = bigstack_left() / kCacheline;
      const uint64_t full_load_vecs_req = sample_ct * S_CAST(uint64_t, NypCtToVecCt(variant_ct));
      uintptr_t* plink1_smaj_loadbuf;
      uint32_t vblock_group_ct;
      uint32_t cur_vidx_ct;
      if (full_load_vecs_req > cachelines_avail * kVecsPerCacheline) {
        // each iteration requires ((kPglVblockSize / 4) * calc_thread_ct *
        //   sample_ct) bytes to be loaded
        vblock_group_ct = cachelines_avail / ((kPglVblockSize / (4 * kCacheline)) * calc_thread_ct * S_CAST(uintptr_t, sample_ct));
        assert(vblock_group_ct);
        cur_vidx_ct = vblock_group_ct * calc_thread_ct * kPglVblockSize;
        plink1_smaj_loadbuf = S_CAST(uintptr_t*, bigstack_alloc_raw_rd((cur_vidx_ct / 4) * S_CAST(uintptr_t, sample_ct)));
        // bugfix (18 Nov 2017): this may be larger than variant_ct
        if (cur_vidx_ct > variant_ct) {
          cur_vidx_ct = variant_ct;
          vblock_group_ct = 1 + (cur_vidx_ct - 1) / (kPglVblockSize * calc_thread_ct);
        }
      } else {
        vblock_group_ct = 1 + ((variant_ct - 1) / (calc_thread_ct * kPglVblockSize));
        cur_vidx_ct = variant_ct;
        plink1_smaj_loadbuf = S_CAST(uintptr_t*, bigstack_alloc_raw_rd(full_load_vecs_req * kBytesPerVec));
      }
      uint32_t cur_vidx_ct4 = NypCtToByteCt(cur_vidx_ct);
      uint32_t cur_vidx_ctaw2 = NypCtToAlignedWordCt(cur_vidx_ct);
      const uint32_t pass_ct = 1 + (variant_ct - 1) / cur_vidx_ct;
      ctx.sample_ct = sample_ct;
      ctx.loadbuf_ul_stride = NypCtToVecCt(cur_vidx_ct) * kWordsPerVec;
      ctx.allele_flips_iter = nullptr;
      if (unlikely(SetThreadCt0(calc_thread_ct - 1, &tg))) {
        goto Plink1SampleMajorToPgen_ret_NOMEM;
      }
      SetThreadFuncAndData(Plink1SmajTransposeThread, &ctx, &tg);
      uint32_t pass_idx1 = 0;
      for (uint32_t cur_vidx_base = 0; ; ) {
        uint32_t cur_raw_load_batch_size = raw_load_batch_size;
        uintptr_t* smaj_loadbuf_iter = plink1_smaj_loadbuf;
        putc_unlocked('\r', stdout);
        ++pass_idx1;
        printf("Pass %u/%u: loading... 0%%", pass_idx1, pass_ct);
        fflush(stdout);
        uint32_t pct = 0;
        uint32_t next_print_idx = (raw_load_batch_ct + 99) / 100;
        const uint64_t seek_addl_offset = 3 + cur_vidx_base / 4;
        for (uint32_t raw_load_batch_idx = 0; ; ) {
          if (raw_load_batch_size == 1) {
            if (unlikely(fseeko(infile, seek_addl_offset + raw_load_batch_idx * S_CAST(uint64_t, variant_ct4), SEEK_SET))) {
              goto Plink1SampleMajorToPgen_ret_READ_FAIL;
            }
            if (unlikely(!fread_unlocked(smaj_loadbuf_iter, cur_vidx_ct4, 1, infile))) {
              goto Plink1SampleMajorToPgen_ret_READ_FAIL;
            }
            smaj_loadbuf_iter = &(smaj_loadbuf_iter[cur_vidx_ctaw2]);
          } else {
            if (unlikely(!fread_unlocked(raw_loadbuf, cur_raw_load_batch_size * variant_ct4, 1, infile))) {
              goto Plink1SampleMajorToPgen_ret_READ_FAIL;
            }
            unsigned char* raw_loadbuf_iter = &(raw_loadbuf[cur_vidx_base / 4]);
            for (uint32_t uii = 0; uii != cur_raw_load_batch_size; ++uii) {
              memcpy(smaj_loadbuf_iter, raw_loadbuf_iter, cur_vidx_ct4);
              raw_loadbuf_iter = &(raw_loadbuf_iter[variant_ct4]);
              smaj_loadbuf_iter = &(smaj_loadbuf_iter[cur_vidx_ctaw2]);
            }
          }
          ++raw_load_batch_idx;
          if (raw_load_batch_idx >= raw_load_batch_ct_m1) {
            if (raw_load_batch_idx > raw_load_batch_ct_m1) {
              break;
            }
            cur_raw_load_batch_size = sample_ct - raw_load_batch_idx * raw_load_batch_size;
          }
          if (raw_load_batch_idx >= next_print_idx) {
            if (pct > 10) {
              putc_unlocked('\b', stdout);
            }
            pct = (raw_load_batch_idx * 100LLU) / raw_load_batch_ct;
            printf("\b\b%u%%", pct++);
            fflush(stdout);
            next_print_idx = (pct * S_CAST(uint64_t, raw_load_batch_ct) + 99) / 100;
          }
        }
        const uintptr_t last_tidx = calc_thread_ct - 1;
        uint32_t load_idx = 0;
        ctx.cur_block_write_ct = calc_thread_ct * kPglVblockSize;
        putc_unlocked('\r', stdout);
        printf("Pass %u/%u: transposing and compressing... 0%%", pass_idx1, pass_ct);
        ReinitThreads(&tg);
        pct = 0;
        next_print_idx = (vblock_group_ct + 99) / 100;
        do {
          if (load_idx >= next_print_idx) {
            if (pct > 10) {
              putc_unlocked('\b', stdout);
            }
            pct = (load_idx * 100LLU) / vblock_group_ct;
            printf("\b\b%u%%", pct++);
            fflush(stdout);
            next_print_idx = (pct * S_CAST(uint64_t, vblock_group_ct) + 99) / 100;
          }
          ctx.plink1_smaj_loadbuf_iter = &(plink1_smaj_loadbuf[load_idx * calc_thread_ct * (kPglVblockSize / kBitsPerWordD2)]);
          if (allele_flips) {
            ctx.allele_flips_iter = &(allele_flips[load_idx * calc_thread_ct * (kPglVblockSize / kBitsPerWord)]);
          }
          if (++load_idx == vblock_group_ct) {
            DeclareLastThreadBlock(&tg);
            ctx.cur_block_write_ct = cur_vidx_ct - (load_idx - 1) * calc_thread_ct * kPglVblockSize;
          }
          if (last_tidx) {
            if (unlikely(SpawnThreads(&tg))) {
              goto Plink1SampleMajorToPgen_ret_THREAD_CREATE_FAIL;
            }
          }
          Plink1SmajTransposeMain(last_tidx, &ctx);
          if (last_tidx) {
            JoinThreads(&tg);
            // Plink1SmajTransposeThread() never errors out
          }
          reterr = MpgwFlush(mpgwp);
          if (unlikely(reterr)) {
            goto Plink1SampleMajorToPgen_ret_WRITE_FAIL;
          }
        } while (!IsLastBlock(&tg));
        cur_vidx_base += cur_vidx_ct;
        if (cur_vidx_base == variant_ct) {
          if (pct > 10) {
            putc_unlocked('\b', stdout);
          }
          break;
        }
        fputs("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b                     ", stdout);
        // assumes PgfiInitPhase1() leaves file pointer at byte 3; otherwise,
        // necessary to put this at top of main loop
        if (unlikely(fseeko(infile, 3, SEEK_SET))) {
          goto Plink1SampleMajorToPgen_ret_READ_FAIL;
        }
        if (variant_ct - cur_vidx_base <= cur_vidx_ct) {
          cur_vidx_ct = variant_ct - cur_vidx_base;
          cur_vidx_ct4 = NypCtToByteCt(cur_vidx_ct);
          cur_vidx_ctaw2 = NypCtToAlignedWordCt(cur_vidx_ct);
          ctx.loadbuf_ul_stride = NypCtToVecCt(cur_vidx_ct) * kWordsPerVec;
          vblock_group_ct = 1 + (cur_vidx_ct - 1) / (kPglVblockSize * calc_thread_ct);
        }
      }
      mpgwp = nullptr;
      fputs("\b\bdone.\n", stdout);
      logprintf("Transpose complete.\n");
    } else {
    Plink1SampleMajorToPgen_fallback:
      BigstackReset(bigstack_mark2);
      mpgwp = nullptr;
      reterr = Plink1SampleMajorToPgenLowmem(pgenname, allele_flips, variant_ct, sample_ct, real_ref_alleles, raw_load_batch_size, raw_load_batch_ct, infile, raw_loadbuf);
    }
  }
  while (0) {
  Plink1SampleMajorToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  Plink1SampleMajorToPgen_ret_READ_FAIL:
    if (feof_unlocked(infile)) {
      errno = 0;
    }
    logputs("\n");
    logerrprintfww(kErrprintfFread, ".bed file", rstrerror(errno));
    reterr = kPglRetReadFail;
    break;
  Plink1SampleMajorToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  Plink1SampleMajorToPgen_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  Plink1SampleMajorToPgen_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
 Plink1SampleMajorToPgen_ret_1:
  CleanupThreads(&tg);
  CleanupMpgw(mpgwp, &reterr);
  BigstackReset(bigstack_mark);
  return reterr;
}

PglErr PedmapToPgen(const char* pedname, const char* mapname, const char* missing_catname, MiscFlags misc_flags, ImportFlags import_flags, LoadFilterLogFlags load_filter_log_import_flags, uint32_t psam_01, FamCol fam_cols, int32_t missing_pheno, char input_missing_geno_char, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip) {
  unsigned char* bigstack_mark = g_bigstack_base;
  unsigned char* bigstack_end_mark = g_bigstack_end;
  FILE* indmaj_bed_file = nullptr;
  FILE* tmp_fam_file = nullptr;
  TextStream ped_txs;
  PreinitTextStream(&ped_txs);
  uintptr_t line_idx = 0;

  char* pvar_cswritep = nullptr;
  CompressStreamState pvar_css;
  PreinitCstream(&pvar_css);
  PglErr reterr = kPglRetSuccess;
  {
    // 1. Scan .map file, count number of variants and check whether nonzero CM
    //    values exist.
    // 2. Convert .ped to temporary sample-major .bed file, using the variant
    //    count to disambiguate the compound-genotypes case from the regular
    //    case.  Write temporary .fam file as well, to be reread with same
    //    fam_cols setting.
    // 3. RewritePsam() call.
    // 4. Now that we know the alleles and their counts, populate allele_flips.
    // 5. Rescan .map file, write .pvar; then free allele storage space.
    // 6. Transpose sample-major .bed to .pgen, while flipping the necessary
    //    variants.
    // 7. Delete temporary files.
    FinalizeChrset(load_filter_log_import_flags, cip);
    uintptr_t* variant_include = nullptr;
    uint32_t raw_variant_ct;
    uint32_t variant_ct;
    uint32_t neg_bp_seen;
    uint32_t at_least_one_nzero_cm;
    reterr = ScanMap(mapname, misc_flags, load_filter_log_import_flags, cip, &raw_variant_ct, &variant_ct, &neg_bp_seen, &at_least_one_nzero_cm, &variant_include);
    if (unlikely(reterr)) {
      goto PedmapToPgen_ret_1;
    }
    if (unlikely(raw_variant_ct > (kMaxLongLine / 4))) {
      logerrputs("Error: Too many variants for .ped file converter.\n");
      goto PedmapToPgen_ret_MALFORMED_INPUT;
    }
    {
      char* write_iter = strcpya_k(g_logbuf, "--pedmap: ");
      write_iter = u32toa(raw_variant_ct, write_iter);
      write_iter = strcpya_k(write_iter, " variant");
      if (line_idx != 2) {
        *write_iter++ = 's';
      }
      write_iter = strcpya_k(write_iter, " in .map file");
      if (raw_variant_ct != variant_ct) {
        write_iter = strcpya_k(write_iter, "; ");
        write_iter = u32toa(raw_variant_ct - variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " excluded by ");
        if (load_filter_log_import_flags) {
          AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        }
        if (neg_bp_seen) {
          if (load_filter_log_import_flags) {
            write_iter = strcpya_k(write_iter, " and/or ");
          }
          write_iter = strcpya_k(write_iter, "negative bp coordinates");
        }
        write_iter = strcpya_k(write_iter, ", ");
        write_iter = u32toa(variant_ct, write_iter);
        write_iter = strcpya_k(write_iter, " remaining");
      } else if (load_filter_log_import_flags) {
        write_iter = strcpya_k(write_iter, " (");
        AppendLoadFilterFlagnames(load_filter_log_import_flags, &write_iter);
        write_iter = strcpya_k(write_iter, " had no effect)");
      }
      strcpy_k(write_iter, ".\n");
      WordWrapB(0);
      logputsb();
    }

    const uint32_t variant_ctl = BitCtToWordCt(variant_ct);
    const uint32_t variant_ctl2 = NypCtToWordCt(variant_ct);
    const uint32_t variant_ct4 = NypCtToByteCt(variant_ct);
    uintptr_t* allele_flips;
    const char** allele_codes;
    uint32_t* second_allele_plus_missing_cts;
    uintptr_t* genovec;
    if (unlikely(bigstack_calloc_w(variant_ctl, &allele_flips) ||
                 bigstack_end_alloc_kcp(2 * variant_ct, &allele_codes) ||
                 bigstack_end_calloc_u32(variant_ct, &second_allele_plus_missing_cts) ||
                 bigstack_end_alloc_w(variant_ctl2, &genovec))) {
      goto PedmapToPgen_ret_NOMEM;
    }
    const char* null_str = &(g_one_char_strs[0]);
    for (uint32_t uii = 0; uii != 2 * variant_ct; ++uii) {
      allele_codes[uii] = null_str;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".bed.smaj");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &indmaj_bed_file))) {
      goto PedmapToPgen_ret_OPEN_FAIL;
    }
    if (unlikely(!fwrite_unlocked("l\x1b", 3, 1, indmaj_bed_file))) {
      goto PedmapToPgen_ret_WRITE_FAIL;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".fam.tmp");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &tmp_fam_file))) {
      goto PedmapToPgen_ret_OPEN_FAIL;
    }
    char* fam_writebuf = g_textbuf;
    char* fam_writebuf_flush = &(fam_writebuf[kMaxMediumLine]);
    char* fam_write_iter = fam_writebuf;
    const uint32_t decompress_thread_ct = MAXV(1, max_thread_ct - 1);
    reterr = SizeAndInitTextStream(pedname, bigstack_left() / 4, decompress_thread_ct, &ped_txs);
    if (unlikely(reterr)) {
      goto PedmapToPgen_ret_TSTREAM_FAIL;
    }
    unsigned char* tmp_alloc_base = g_bigstack_base;
    unsigned char* tmp_alloc_end = g_bigstack_end;
    const uint32_t fam_col_ct_m1 = ((fam_cols / kfFamCol1) & 1) + ((fam_cols / (kfFamCol34 / 2)) & 2) + ((fam_cols / kfFamCol5) & 1) + ((fam_cols / kfFamCol6) & 1);
    char* ped_line_start;
    do {
      ++line_idx;
      ped_line_start = TextGet(&ped_txs);
      if (unlikely(!ped_line_start)) {
        if (TextStreamErrcode2(&ped_txs, &reterr)) {
          goto PedmapToPgen_ret_TSTREAM_FAIL;
        }
        snprintf(g_logbuf, kLogbufSize, "Error: %s is empty.\n", pedname);
        goto PedmapToPgen_ret_DEGENERATE_DATA;
      }
    } while (ped_line_start[0] == '#');
    uint32_t compound_genotypes = 0;
    {
      const uint32_t token_ct = CountTokens(ped_line_start);
      const uint32_t compound_genotypes_token_ct = fam_col_ct_m1 + 1 + raw_variant_ct;
      if (token_ct == compound_genotypes_token_ct) {
        compound_genotypes = 1;
      } else {
        const uint32_t regular_token_ct = compound_genotypes_token_ct + raw_variant_ct;
        if (unlikely(token_ct != regular_token_ct)) {
          logerrprintfww("Error: Unexpected number of columns in .ped file (%u or %u expected).\n", compound_genotypes_token_ct, regular_token_ct);
          goto PedmapToPgen_ret_INCONSISTENT_INPUT;
        }
      }
    }
    uint32_t max_allele_slen = 1;
    uint32_t sample_idx = 0;
    while (1) {
      char* fam_last = NextTokenMult0(ped_line_start, fam_col_ct_m1);
      if (unlikely(!fam_last)) {
        goto PedmapToPgen_ret_MISSING_TOKENS;
      }
      char* fam_end = CurTokenEnd(fam_last);
      const uintptr_t fam_copy_byte_ct = fam_end - ped_line_start;
      if (unlikely(fam_copy_byte_ct > kMaxLongLine - strlen(EOLN_STR))) {
        logerrprintfww("Error: A leading column of line %" PRIuPTR " of %s is too long.\n", line_idx, pedname);
        goto PedmapToPgen_ret_MALFORMED_INPUT;
      }
      fam_write_iter = memcpya(fam_write_iter, ped_line_start, fam_copy_byte_ct);
      AppendBinaryEoln(&fam_write_iter);
      if (unlikely(fwrite_ck(fam_writebuf_flush, tmp_fam_file, &fam_write_iter))) {
        goto PedmapToPgen_ret_WRITE_FAIL;
      }
      char* ped_iter = fam_end;
      uint32_t variant_idx = 0;
      uint32_t variant_idx_lowbits_x2 = 0;
      uintptr_t geno_word = 0;
      if (!compound_genotypes) {
        for (uint32_t raw_variant_idx = 0; ; ++raw_variant_idx) {
          char* first_allele_start = FirstNonTspace(ped_iter);
          char* first_allele_end = FirstSpaceOrEoln(first_allele_start);
          char* second_allele_start = FirstNonTspace(first_allele_end);
          if (unlikely(IsSpaceOrEoln(*second_allele_start))) {
            goto PedmapToPgen_ret_MISSING_TOKENS;
          }
          ped_iter = CurTokenEnd(second_allele_start);
          if (variant_include && (!IsSet(variant_include, raw_variant_idx))) {
            continue;
          }
          const uint32_t first_allele_slen = first_allele_end - first_allele_start;
          const uint32_t second_allele_slen = ped_iter - second_allele_start;
          const uint32_t is_het = (first_allele_slen != second_allele_slen) || (!memequal(first_allele_start, second_allele_start, first_allele_slen));
          const uint32_t variant_idx_x2 = variant_idx * 2;
          const uint32_t variant_idx_x2_p1 = variant_idx_x2 + 1;
          const char* prov_ref_allele = allele_codes[variant_idx_x2];
          const char* prov_alt_allele = allele_codes[variant_idx_x2_p1];
          // Note that this uses PLINK 1 encoding.
          uintptr_t cur_geno;
          if (strequal_unsafe(prov_ref_allele, first_allele_start, first_allele_slen)) {
          PedmapToPgen_ref_x:
            cur_geno = 3 - is_het;
            if (is_het) {
              if (!strequal_unsafe(prov_alt_allele, second_allele_start, second_allele_slen)) {
                if (unlikely(prov_alt_allele != null_str)) {
                  if ((second_allele_slen == 1) && ((second_allele_start[0] == '.') || (second_allele_start[0] == input_missing_geno_char))) {
                    goto PedmapToPgen_ret_HALF_MISSING;
                  }
                  goto PedmapToPgen_ret_MULTIALLELIC;
                }
                if (second_allele_slen == 1) {
                  allele_codes[variant_idx_x2_p1] = &(g_one_char_strs[ctou32(second_allele_start[0]) * 2]);
                } else {
                  if (unlikely(StoreStringAtEndK(tmp_alloc_base, second_allele_start, second_allele_slen, &tmp_alloc_end, &(allele_codes[variant_idx_x2_p1])))) {
                    goto PedmapToPgen_ret_NOMEM;
                  }
                  if (second_allele_slen > max_allele_slen) {
                    max_allele_slen = second_allele_slen;
                  }
                }
              }
            }
          } else if (strequal_unsafe(prov_alt_allele, first_allele_start, first_allele_slen)) {
          PedmapToPgen_alt_x:
            cur_geno = is_het * 2;
            if (is_het) {
              if (unlikely(!strequal_unsafe(prov_ref_allele, second_allele_start, second_allele_slen))) {
                if ((second_allele_slen == 1) && ((second_allele_start[0] == '.') || (second_allele_start[0] == input_missing_geno_char))) {
                  goto PedmapToPgen_ret_HALF_MISSING;
                }
                goto PedmapToPgen_ret_MULTIALLELIC;
              }
            }
          } else if ((first_allele_slen == 1) && ((first_allele_start[0] == '.') || (first_allele_start[0] == input_missing_geno_char))) {
            if (unlikely(is_het)) {
              goto PedmapToPgen_ret_HALF_MISSING;
            }
            cur_geno = 1;
          } else if (prov_ref_allele == null_str) {
            if (first_allele_slen == 1) {
              allele_codes[variant_idx_x2] = &(g_one_char_strs[ctou32(first_allele_start[0]) * 2]);
            } else {
              if (unlikely(StoreStringAtEndK(tmp_alloc_base, first_allele_start, first_allele_slen, &tmp_alloc_end, &(allele_codes[variant_idx_x2])))) {
                goto PedmapToPgen_ret_NOMEM;
              }
              if (first_allele_slen > max_allele_slen) {
                max_allele_slen = first_allele_slen;
              }
            }
            goto PedmapToPgen_ref_x;
          } else if (likely(prov_alt_allele == null_str)) {
            if (first_allele_slen == 1) {
              allele_codes[variant_idx_x2_p1] = &(g_one_char_strs[ctou32(first_allele_start[0]) * 2]);
            } else {
              if (unlikely(StoreStringAtEndK(tmp_alloc_base, first_allele_start, first_allele_slen, &tmp_alloc_end, &(allele_codes[variant_idx_x2_p1])))) {
                goto PedmapToPgen_ret_NOMEM;
              }
              if (first_allele_slen > max_allele_slen) {
                max_allele_slen = first_allele_slen;
              }
            }
            goto PedmapToPgen_alt_x;
          } else {
            goto PedmapToPgen_ret_MULTIALLELIC;
          }
          second_allele_plus_missing_cts[variant_idx] += (4 - cur_geno) >> 1;
          geno_word |= cur_geno << variant_idx_lowbits_x2;
          variant_idx_lowbits_x2 += 2;
          if (variant_idx_lowbits_x2 == kBitsPerWord) {
            genovec[variant_idx / kBitsPerWordD2] = geno_word;
            geno_word = 0;
            variant_idx_lowbits_x2 = 0;
          }
          if (++variant_idx == variant_ct) {
            break;
          }
        }
      } else {
        for (uint32_t raw_variant_idx = 0; ; ++raw_variant_idx) {
          char* alleles_start = FirstNonTspace(ped_iter);
          if (unlikely(IsSpaceOrEoln(*alleles_start))) {
            goto PedmapToPgen_ret_MISSING_TOKENS;
          }
          ped_iter = CurTokenEnd(alleles_start);
          if (unlikely(ped_iter != &(alleles_start[2]))) {
            putc_unlocked('\n', stdout);
            logerrprintfww("Error: --pedmap: .map file and number of tokens in first line of .ped file imply that the latter is in the compound-genotypes format, but line %" PRIuPTR " has a genotype that isn't length-2.\n", line_idx);
            goto PedmapToPgen_ret_INCONSISTENT_INPUT;
          }
          if (variant_include && (!IsSet(variant_include, raw_variant_idx))) {
            continue;
          }
          const char first_allele_char = alleles_start[0];
          const char second_allele_char = alleles_start[1];
          const uint32_t is_het = (first_allele_char != second_allele_char);
          const uint32_t variant_idx_x2 = variant_idx * 2;
          const uint32_t variant_idx_x2_p1 = variant_idx_x2 + 1;
          const char prov_ref_char = allele_codes[variant_idx_x2][0];
          const char prov_alt_char = allele_codes[variant_idx_x2_p1][0];
          // Note that this uses PLINK 1 encoding.
          uintptr_t cur_geno;
          if (prov_ref_char == first_allele_char) {
          PedmapToPgen_refchar_x:
            cur_geno = 3 - is_het;
            if (is_het) {
              if (prov_alt_char != second_allele_char) {
                if (unlikely(prov_alt_char)) {
                  if ((second_allele_char == '.') || (second_allele_char == input_missing_geno_char)) {
                    goto PedmapToPgen_ret_HALF_MISSING;
                  }
                  goto PedmapToPgen_ret_MULTIALLELIC;
                }
                allele_codes[variant_idx_x2_p1] = &(g_one_char_strs[ctou32(second_allele_char) * 2]);
              }
            }
          } else if (prov_alt_char == first_allele_char) {
          PedmapToPgen_altchar_x:
            cur_geno = is_het * 2;
            if (is_het) {
              if (unlikely(prov_ref_char != second_allele_char)) {
                if ((second_allele_char == '.') || (second_allele_char == input_missing_geno_char)) {
                  goto PedmapToPgen_ret_HALF_MISSING;
                }
                goto PedmapToPgen_ret_MULTIALLELIC;
              }
            }
          } else if ((first_allele_char == '.') || (first_allele_char == input_missing_geno_char)) {
            if (unlikely(is_het)) {
              goto PedmapToPgen_ret_HALF_MISSING;
            }
            cur_geno = 1;
          } else if (!prov_ref_char) {
            allele_codes[variant_idx_x2] = &(g_one_char_strs[ctou32(first_allele_char) * 2]);
            goto PedmapToPgen_refchar_x;
          } else if (likely(!prov_alt_char)) {
            allele_codes[variant_idx_x2_p1] = &(g_one_char_strs[ctou32(first_allele_char) * 2]);
            goto PedmapToPgen_altchar_x;
          } else {
            goto PedmapToPgen_ret_MULTIALLELIC;
          }
          second_allele_plus_missing_cts[variant_idx] += (4 - cur_geno) >> 1;
          geno_word |= cur_geno << variant_idx_lowbits_x2;
          variant_idx_lowbits_x2 += 2;
          if (variant_idx_lowbits_x2 == kBitsPerWord) {
            genovec[variant_idx / kBitsPerWordD2] = geno_word;
            geno_word = 0;
            variant_idx_lowbits_x2 = 0;
          }
          if (++variant_idx == variant_ct) {
            break;
          }
        }
      }
      if (variant_idx_lowbits_x2) {
        genovec[variant_ct / kBitsPerWordD2] = geno_word;
      }
      if (unlikely(!fwrite_unlocked(genovec, variant_ct4, 1, indmaj_bed_file))) {
        goto PedmapToPgen_ret_WRITE_FAIL;
      }
      ++sample_idx;
      if ((sample_idx % 100) == 0) {
        printf("\r--pedmap: %u samples scanned.", sample_idx);
        fflush(stdout);
      }
      do {
        ++line_idx;
        ped_line_start = TextGet(&ped_txs);
        if (!ped_line_start) {
          goto PedmapToPgen_ped_eof;
        }
      } while (ped_line_start[0] == '#');
      if (unlikely(sample_idx == 0x7fffffff)) {
        logerrputs("Error: " PROG_NAME_STR " does not support more than 2^31 - 2 samples.\n");
        goto PedmapToPgen_ret_MALFORMED_INPUT;
      }
    }
  PedmapToPgen_ped_eof:
    if (unlikely(TextStreamErrcode2(&ped_txs, &reterr))) {
      goto PedmapToPgen_ret_TSTREAM_FAIL;
    }
    if (unlikely(fclose_flush_null(fam_writebuf_flush, fam_write_iter, &tmp_fam_file) ||
                 fclose_null(&indmaj_bed_file))) {
      goto PedmapToPgen_ret_WRITE_FAIL;
    }
    const uint32_t sample_ct = sample_idx;
    *outname_end = '\0';
    putc_unlocked('\r', stdout);
    logprintfww("--pedmap: %u sample%s present, genotypes extracted to %s.bed.smaj .\n", sample_ct, (sample_ct == 1)? "" : "s", outname);
    *outname_end = '.';
    BigstackEndSet(tmp_alloc_end);

    reterr = RewritePsam(outname, missing_catname, misc_flags, fam_cols, missing_pheno, psam_01, max_thread_ct, outname, outname_end, nullptr);
    if (unlikely(reterr)) {
      goto PedmapToPgen_ret_1;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".fam.tmp");
    if (unlikely(unlink(outname))) {
      goto PedmapToPgen_ret_WRITE_FAIL;
    }

    for (uint32_t variant_idx = 0; variant_idx != variant_ct; ++variant_idx) {
      if (second_allele_plus_missing_cts[variant_idx] > sample_ct) {
        SetBit(variant_idx, allele_flips);
        const char* tmp_allele_code = allele_codes[variant_idx * 2];
        allele_codes[variant_idx * 2] = allele_codes[variant_idx * 2 + 1];
        allele_codes[variant_idx * 2 + 1] = tmp_allele_code;
      } else if (allele_codes[variant_idx * 2 + 1] == null_str) {
        // Restore standard representation.
        allele_codes[variant_idx * 2 + 1] = nullptr;
        if (allele_codes[variant_idx * 2] == null_str) {
          allele_codes[variant_idx * 2] = nullptr;
        }
      }
    }

    reterr = MapToPvar(mapname, cip, allele_codes, variant_ct, max_allele_slen, import_flags, at_least_one_nzero_cm, outname, outname_end);
    if (unlikely(reterr)) {
      goto PedmapToPgen_ret_1;
    }
    BigstackEndReset(bigstack_end_mark);

    snprintf(outname_end, kMaxOutfnameExtBlen, ".bed.smaj");
    if (unlikely(fopen_checked(outname, FOPEN_RB, &indmaj_bed_file))) {
      goto PedmapToPgen_ret_OPEN_FAIL;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    reterr = Plink1SampleMajorToPgen(outname, allele_flips, variant_ct, sample_ct, 0, max_thread_ct, indmaj_bed_file);
    if (unlikely(reterr)) {
      goto PedmapToPgen_ret_1;
    }
    // can't do anything meaningful with an fclose error here
    fclose(indmaj_bed_file);
    indmaj_bed_file = nullptr;
    snprintf(outname_end, kMaxOutfnameExtBlen, ".bed.smaj");
    if (unlikely(unlink(outname))) {
      goto PedmapToPgen_ret_WRITE_FAIL;
    }

    char* write_iter = strcpya_k(g_logbuf, "--pedmap: ");
    const uint32_t outname_base_slen = outname_end - outname;
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pgen + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    write_iter = strcpya_k(write_iter, ".pvar");
    if (import_flags & kfImportKeepAutoconvVzs) {
      write_iter = strcpya_k(write_iter, ".zst");
    }
    write_iter = strcpya_k(write_iter, " + ");
    write_iter = memcpya(write_iter, outname, outname_base_slen);
    strcpy_k(write_iter, ".psam written. .bed.smaj and .fam.tmp temporary files deleted.\n");
    WordWrapB(0);
    logputsb();
  }
  while (0) {
  PedmapToPgen_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  PedmapToPgen_ret_TSTREAM_FAIL:
    TextStreamErrPrint(pedname, &ped_txs);
    break;
  PedmapToPgen_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  PedmapToPgen_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  PedmapToPgen_ret_MISSING_TOKENS:
    putc_unlocked('\n', stdout);
    logerrprintf("Error: Line %" PRIuPTR " of .ped file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  PedmapToPgen_ret_HALF_MISSING:
    putc_unlocked('\n', stdout);
    logerrprintfww("Error: Half-missing genotype on line %" PRIuPTR " of %s.\n", line_idx, pedname);
    reterr = kPglRetMalformedInput;
    break;
  PedmapToPgen_ret_MULTIALLELIC:
    putc_unlocked('\n', stdout);
    logerrprintfww("Error: Multiallelic variant in %s. This violates the .ped specification; please reformat the file as e.g. VCF.\n", pedname);
    reterr = kPglRetMalformedInput;
    break;
  PedmapToPgen_ret_MALFORMED_INPUT:
    putc_unlocked('\n', stdout);
    reterr = kPglRetMalformedInput;
    break;
  PedmapToPgen_ret_INCONSISTENT_INPUT:
    putc_unlocked('\n', stdout);
    reterr = kPglRetInconsistentInput;
    break;
  PedmapToPgen_ret_DEGENERATE_DATA:
    reterr = kPglRetDegenerateData;
    break;
  }
 PedmapToPgen_ret_1:
  CswriteCloseCond(&pvar_css, pvar_cswritep);
  CleanupTextStream2(pedname, &ped_txs, &reterr);
  fclose_cond(tmp_fam_file);
  fclose_cond(indmaj_bed_file);
  BigstackDoubleReset(bigstack_mark, bigstack_end_mark);
  return reterr;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
