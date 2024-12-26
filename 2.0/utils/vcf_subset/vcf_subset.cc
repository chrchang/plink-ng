// Copyright (C) 2024 Christopher Chang.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "../../include/plink2_bits.h"
#include "../../include/plink2_htable.h"
#include "../../include/plink2_memory.h"
#include "../../include/plink2_text.h"

#ifdef __cplusplus
namespace plink2 {
#endif

static const char ver_str[] = "vcf_subset v1.0.0"

#ifdef __LP64__
#  ifdef USE_AVX2
    " AVX2"
#  elif defined(__APPLE__) && !defined(__x86_64__)
    " M1"
#  elif defined(USE_SSE42)
    " SSE4.2"
#  else
    " 64-bit"
#  endif
#else
  " 32-bit"
#endif

  " (26 Dec 2024)";
static const char ver_str2[] =
  // include leading space if day < 10, so character length stays the same
  ""

  // length of architecture string + this string should be 7 characterrs
#ifdef __LP64__
#  ifdef USE_AVX2
    "  "
#  elif defined(__APPLE__) && !defined(__x86_64__)
    "    "
#  elif defined(USE_SSE42)
    ""
#  else
    ""
#  endif
#else
  ""
#endif
  "              cog-genomics.org/plink/2.0/\n"
  "(C) 2024 Christopher Chang                                          GNU LGPL v3\n";

void PrintVer() {
  fputs(ver_str, stdout);
  fputs(ver_str2, stdout);
}

void DispUsage(FILE* dst_stream) {
  PrintVer();
  fputs(
"Usage: vcf_subset [options...] --in <input VCF path> --out <output VCF path>\n"
"Optional flags:\n"
"  --allow-any-paths  : By default, the program errors out if the input and\n"
"                       output paths don't end in '.vcf', '.vcf.gz', or\n"
"                       '.vcf.bgz'.  Use this flag to remove the restriction.\n"
"  --allow-no-samples : By default, the program errors out if --indv and/or\n"
"                       --keep causes no samples to remain.  Use this flag to\n"
"                       remove the restriction.\n"
"  --bgz-level <#>    : Set bgzip compression level (1-12, default 6).  Note\n"
"                       that compression occurs iff the --out path ends in '.gz'\n"
"                       or '.bgz'.\n"
"  --erase-info       : Remove INFO header lines, set all INFO entries to '.'.\n"
"  --help             : Just display this help text.\n"
"  --indv <ID>        : Keep only the sample with the given ID.\n"
"  --keep <filename>  : Keep only samples with IDs in the given file.\n"
"  --memory <#>       : Set size, in MiB, of initial workspace malloc attempt.\n"
"  --no-samples       : Remove all samples.\n"
"  --threads <#>      : Set maximum number of compute threads.\n"
"  --version          : Just display version string.\n"
, dst_stream);
}

// Mostly used as a buffer for error messages.
CONSTI32(kTextbufSize, 2 * kMaxMediumLine + 256);
char g_textbuf[kTextbufSize];

void TextErrPrint(const char* file_descrip, const char* errmsg, PglErr reterr) {
  assert(reterr != kPglRetSuccess);
  if (reterr == kPglRetOpenFail) {
    fprintf(stderr, kErrprintfFopen, file_descrip, errmsg);
  } else if (reterr == kPglRetReadFail) {
    fprintf(stderr, kErrprintfFread, file_descrip, errmsg);
  } else if (reterr == kPglRetDecompressFail) {
    fprintf(stderr, kErrprintfDecompress, file_descrip, errmsg);
  } else if (reterr == kPglRetMalformedInput) {
    if (errmsg == kShortErrInteriorEmptyLine) {
      fprintf(stderr, "Error: Unexpected interior empty line in %s.\n", file_descrip);
    } else {
      assert(errmsg == kShortErrLongLine);
      fprintf(stderr, "Error: Pathologically long line in %s.\n", file_descrip);
    }
  }
  // kPglRetRewindFail not used here
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

PglErr VcfSubset(unsigned char* bigstack_base, unsigned char* bigstack_end, const char* inpath, const char* indv_str, const char* keeppath, const char* outpath, uint32_t outpath_slen, uint32_t allow_no_samples, uint32_t bgz_level, uint32_t erase_info, uint32_t no_samples, uint32_t max_thread_ct) {
  uintptr_t line_idx = 1;
  PglErr reterr = kPglRetSuccess;
  TextStream in_txs;
  BgzfCompressStream out_bgzf;
  PreinitTextStream(&in_txs);
  PreinitBgzfCompressStream(&out_bgzf);
  {
    const uint32_t max_line_blen = MINV((bigstack_end - bigstack_base) / 4, kMaxLongLine);
    uintptr_t max_header_blen = S_CAST(uintptr_t, bigstack_end - bigstack_base) / 2;
#ifdef __LP64__
    if (max_header_blen > (1LLU << 32)) {
      max_header_blen = 1LLU << 32;
    }
#endif
    char* readbuf;
    char* writebuf;
    if (unlikely(arena_alloc_c(bigstack_end, max_line_blen + kDecompressChunkSize, &bigstack_base, &readbuf) ||
                 arena_alloc_c(bigstack_end, max_line_blen + kMaxMediumLine, &bigstack_base, &writebuf))) {
      goto VcfSubset_ret_NOMEM;
    }
    char* write_iter = writebuf;
    char* writebuf_flush = &(writebuf[kMaxMediumLine]);
    // If output is uncompressed, allow reader to use MAX(1, max_thread_ct - 2)
    // threads, leaving one thread for copy-subsetting and one thread for
    // flushing.
    // Otherwise, start by reserving ~half of threads to reader.  Then, scan
    // final header line, estimate what fraction of the file we're keeping
    // (currently affected by --indv, --keep, --no-samples), and set the number
    // of compression threads under the assumption that compression is 5x as
    // expensive as decompression (and consequently some of the decompression
    // threads may be expected to be fairly idle).
    uint32_t clvl = 0;  // 0 = no compression
    uint32_t decompress_thread_ct = 1;
    if (StrEndsWith(outpath, ".gz", outpath_slen) ||
        StrEndsWith(outpath, ".bgz", outpath_slen)) {
      clvl = bgz_level;
      if (max_thread_ct >= 4) {
        decompress_thread_ct = max_thread_ct / 2;
      }
    } else {
      if (max_thread_ct >= 4) {
        decompress_thread_ct = max_thread_ct - 2;
      }
    }
    reterr = TextStreamOpenEx(inpath, kMaxLongLine, max_line_blen + kDecompressChunkSize, decompress_thread_ct, nullptr, readbuf, &in_txs);
    if (unlikely(reterr)) {
      goto VcfSubset_ret_IN_TSTREAM_FAIL;
    }
    // Can only be >1 if bgzf.  (But that's the expected case.)
    decompress_thread_ct = MAXV(1, TextDecompressThreadCt(&in_txs));

    uintptr_t* sample_include = nullptr;
    char* line_start = TextLineEnd(&in_txs);
    uint32_t raw_sample_ct = 0;
    uint32_t raw_sample_ctl = 0;
    uint32_t sample_ct = 0;
    {
      unsigned char* bigstack_mark2 = bigstack_base;
      char* header = R_CAST(char*, bigstack_base);
      char* header_limit = &(header[max_header_blen]);
      char* header_write_iter = header;
      for (; ; ++line_idx) {
        reterr = TextNextLineUnsafe(&in_txs, &line_start);
        if (unlikely(reterr)) {
          if (reterr == kPglRetEof) {
            fputs("Error: No #CHROM header line or variant records in --in file.\n", stderr);
            goto VcfSubset_ret_MALFORMED_INPUT;
          }
          goto VcfSubset_ret_IN_TSTREAM_FAIL;
        }
        if (unlikely(*line_start != '#')) {
          if ((line_idx == 1) && memequal_sk(line_start, "BCF")) {
            if (line_start[3] == 2) {
              snprintf(g_textbuf, kTextbufSize, "Error: %s appears to be a BCF2 file; this program only supports VCF for now. You can use 'bcftools view' to convert it to a VCF.\n", inpath);
              goto VcfSubset_ret_MALFORMED_INPUT_WW;
            }
            if (line_start[3] == 4) {
              snprintf(g_textbuf, kTextbufSize, "Error: %s appears to be a BCF1 file. You can use 'bcftools view' to convert it to a VCF.\n", inpath);
              goto VcfSubset_ret_MALFORMED_INPUT_WW;
            }
          }
          fputs("Error: No #CHROM header line in --in file.\n", stderr);
          goto VcfSubset_ret_MALFORMED_INPUT;
        }
        if (line_start[1] != '#') {
          break;
        }
        // convert LF / CR-LF into OS-appropriate newline
        char* line_lf = AdvToDelim(line_start, '\n');
        char* line_end = &(line_lf[1]);
        if (erase_info && StrStartsWithUnsafe(&(line_start[2]), "INFO=<")) {
          line_start = line_end;
          continue;
        }
        uintptr_t line_slen = line_lf - line_start;
        if (line_lf[-1] == '\r') {
          --line_slen;
        }
        if (unlikely(S_CAST(uintptr_t, header_limit - header_write_iter) < line_slen + strlen(EOLN_STR))) {
#ifdef __LP64__
          if (max_header_blen == (1LLU << 32)) {
            fputs("Error: VCF header too large (>= 4 GiB, note that this isn't allowed by the BCF\nspecification).\n", stderr);
            goto VcfSubset_ret_MALFORMED_INPUT;
          }
#endif
          fputs("Error: VCF header too large to fit in memory.\n", stderr);
          goto VcfSubset_ret_NOMEM;
        }
        header_write_iter = memcpya(header_write_iter, line_start, line_slen);
        AppendBinaryEoln(&header_write_iter);
        line_start = line_end;
      }

      if (unlikely(!StrStartsWithUnsafe(line_start, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"))) {
        snprintf(g_textbuf, kTextbufSize, "Error: Header line %" PRIuPTR " of --in file does not have expected field sequence after #CHROM.\n", line_idx);
        goto VcfSubset_ret_MALFORMED_INPUT_WW;
      }
      char* chrom_line_lf = AdvToDelim(line_start, '\n');
      uintptr_t chrom_line_slen = chrom_line_lf - line_start;
      if (chrom_line_lf[-1] == '\r') {
        --chrom_line_slen;
      }
      if (unlikely(S_CAST(uintptr_t, header_limit - header_write_iter) < chrom_line_slen + strlen(EOLN_STR))) {
#ifdef __LP64__
        if (max_header_blen == (1LLU << 32)) {
          fputs("Error: VCF header too large (>= 4 GiB, note that this isn't allowed by the BCF\nspecification).\n", stderr);
          goto VcfSubset_ret_MALFORMED_INPUT;
        }
#endif
        fputs("Error: VCF header too large to fit in memory.\n", stderr);
        goto VcfSubset_ret_NOMEM;
      }
      // Note that some sample columns in this last header line may not be
      // copied, so this may overallocate.  Not a big deal, we soon free this
      // memory.
      ArenaBaseSet(&(header_write_iter[chrom_line_slen + strlen(EOLN_STR)]), &bigstack_base);

      char* info_end = &(line_start[strlen("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")]);
      header_write_iter = memcpya(header_write_iter, line_start, info_end - line_start);
      if (StrStartsWithUnsafe(info_end, "\tFORMAT\t")) {
        char* first_sample_id_start = &(info_end[strlen("\tFORMAT\t")]);
        // 1. count number of samples
        // 2. allocate sample_include off bigstack_end
        // 3. construct sample ID hash table, error out if any ID duplicated,
        //    apply --indv/--keep if present
        // 4. generate output header line
        raw_sample_ct = 1 + CountByte(first_sample_id_start, '\t', chrom_line_lf - first_sample_id_start);
        raw_sample_ctl = BitCtToWordCt(raw_sample_ct);
        const char** sample_ids;
        if (unlikely(arena_alloc_kcp(bigstack_end, raw_sample_ct, &bigstack_base, &sample_ids) ||
                     arena_end_alloc_w(bigstack_base, raw_sample_ctl, &bigstack_end, &sample_include))) {
          goto VcfSubset_ret_NOMEM;
        }
        SetAllBits(raw_sample_ct, sample_include);
        const uint32_t raw_sample_ctm1 = raw_sample_ct - 1;
        // Construct array of C-strings, pointing to (now-null-terminated)
        // IDs in the input line buffer.
        char* sample_id_iter = first_sample_id_start;
        for (uint32_t sample_uidx = 0; sample_uidx != raw_sample_ctm1; ++sample_uidx) {
          char* sample_id_end = AdvToDelim(sample_id_iter, '\t');
          *sample_id_end = '\0';
          sample_ids[sample_uidx] = sample_id_iter;
          sample_id_iter = &(sample_id_end[1]);
        }
        char* last_sample_id_end = AdvToDelim(sample_id_iter, '\n');
        if (last_sample_id_end[-1] == '\r') {
          --last_sample_id_end;
        }
        *last_sample_id_end = '\0';
        sample_ids[raw_sample_ctm1] = sample_id_iter;

        uint32_t dup_found;
        if (keeppath) {
          reterr = NondupIdLoad(bigstack_base, bigstack_end, sample_ids, keeppath, raw_sample_ct, raw_sample_ct, MAXV(1, max_thread_ct - 1), sample_include, &dup_found, g_textbuf);
          if (unlikely(reterr)) {
            if (g_textbuf[0]) {
              fputs(g_textbuf, stderr);
            }
            goto VcfSubset_ret_1;
          }
          if (unlikely(dup_found)) {
            fputs("Error: #CHROM header line in --in file contains duplicate sample IDs.\n", stderr);
            goto VcfSubset_ret_MALFORMED_INPUT;
          }
          sample_ct = PopcountWords(sample_include, raw_sample_ctl);
          if (!sample_ct) {
            if (unlikely(!allow_no_samples)) {
              snprintf(g_textbuf, kTextbufSize, "Error: #CHROM header line does not any sample IDs in %s. (Use --allow-no-samples to permit this.)\n", keeppath);
              goto VcfSubset_ret_INCONSISTENT_INPUT_WW;
            }
          }
        } else {
          uint32_t* sample_id_htable;
          uint32_t sample_id_htable_size;
          reterr = AllocAndPopulateNondupHtableMt(bigstack_end, sample_include, sample_ids, raw_sample_ct, MAXV(1, max_thread_ct - 1), &bigstack_base, &sample_id_htable, &sample_id_htable_size, &dup_found);
          if (unlikely(reterr)) {
            goto VcfSubset_ret_1;
          }
          if (unlikely(dup_found)) {
            fputs("Error: #CHROM header line in --in file contains duplicate sample IDs.\n", stderr);
            goto VcfSubset_ret_MALFORMED_INPUT;
          }
          if (no_samples || indv_str) {
            ZeroWArr(raw_sample_ctl, sample_include);
            sample_ct = 0;
            if (indv_str) {
              const uint32_t sample_uidx = IdHtableFind(indv_str, sample_ids, sample_id_htable, strlen(indv_str), sample_id_htable_size);
              if (sample_uidx == UINT32_MAX) {
                if (unlikely(!allow_no_samples)) {
                  snprintf(g_textbuf, kTextbufSize, "Error: #CHROM header line does not contain sample ID '%s'.\n", indv_str);
                  goto VcfSubset_ret_INCONSISTENT_INPUT_WW;
                }
              } else {
                SetBit(sample_uidx, sample_include);
                sample_ct = 1;
              }
            }
          } else {
            sample_ct = raw_sample_ct;
          }
        }
        // Only copy FORMAT column if there's at least one sample.
        if (sample_ct != 0) {
          header_write_iter = memcpya(header_write_iter, info_end, first_sample_id_start - info_end);
          uintptr_t sample_uidx_base = 0;
          uintptr_t cur_bits = sample_include[0];
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            const uint32_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
            header_write_iter = strcpyax(header_write_iter, sample_ids[sample_uidx], '\t');
          }
          --header_write_iter;
        }
      }
      // Now we can provide a good thread count to the compressor (in the
      // likely case where we're compressing).  Open --out file for writing,
      // flush header.
      uint32_t compress_thread_ct = 1;
      if ((max_thread_ct > 2) && clvl) {
        // In this estimate, we treat the INFO column as if it were as
        // expensive as 40 sample columns, and the other header columns as if
        // they were like 10 samples.  Could tune this more and/or make this
        // configurable, but I expect this is good enough the vast majority of
        // the time.
        const uint64_t decompress_weight = 50 + raw_sample_ct;
        const uint64_t compress_weight = 5LLU * (50 - 40 * erase_info + sample_ct);
        const uint64_t denom = compress_weight + decompress_weight;
        compress_thread_ct = MAXV(1, (compress_weight * max_thread_ct + denom / 2) / denom);
      }
      reterr = InitBgzfCompressStreamEx(outpath, 0, clvl, compress_thread_ct, &out_bgzf);
      if (unlikely(reterr)) {
        if (reterr == kPglRetOpenFail) {
          snprintf(g_textbuf, kTextbufSize, kErrprintfFopen, outpath, strerror(errno));
          WordWrap(0, g_textbuf);
          fputs(g_textbuf, stderr);
        }
        goto VcfSubset_ret_1;
      }
      AppendBinaryEoln(&header_write_iter);
      if ((raw_sample_ct == 0) || (sample_ct > 0)) {
        if (unlikely(BgzfWrite(header, header_write_iter - header, &out_bgzf))) {
          goto VcfSubset_ret_WRITE_FAIL;
        }
      } else {
        // Also remove FORMAT header lines.
        for (char* header_reread_iter = header; header_reread_iter < header_write_iter; ) {
          char* read_line_end = AdvPastDelim(header_reread_iter, '\n');
          if (StrStartsWithUnsafe(header_reread_iter, "##FORMAT=<")) {
            header_reread_iter = read_line_end;
            continue;
          }
          write_iter = memcpya(write_iter, header_reread_iter, read_line_end - header_reread_iter);
          if (unlikely(bgzfwrite_ck(writebuf_flush, &out_bgzf, &write_iter))) {
            goto VcfSubset_ret_WRITE_FAIL;
          }
          header_reread_iter = read_line_end;
        }
      }
      bigstack_base = bigstack_mark2;
      line_start = &(chrom_line_lf[1]);
    }
    ++line_idx;
    const uintptr_t first_variant_line_idx = line_idx;
    for (; TextGetUnsafe2(&in_txs, &line_start); ++line_idx) {
      char* line_end = AdvPastDelim(line_start, '\n');
      char* filter_end = AdvToNthDelimChecked(line_start, line_end, 7, '\t');
      if (unlikely(!filter_end)) {
        goto VcfSubset_ret_MISSING_TOKENS;
      }
      char* info_start = &(filter_end[1]);
      char* info_end = strchrnul_n(info_start, '\t');
      if (!erase_info) {
        write_iter = memcpya(write_iter, line_start, info_end - line_start);
      } else {
        write_iter = memcpyax(write_iter, line_start, info_start - line_start, '.');
      }
      if (sample_ct != 0) {
        if (unlikely(*info_end == '\n')) {
          goto VcfSubset_ret_MISSING_TOKENS;
        }
        char* format_start = &(info_end[1]);
        // this count is efficient enough that we may as well do it up front,
        // and remove the corresponding check from the sample-filtering inner
        // loop.
        const uint32_t line_raw_sample_ct = CountByte(format_start, '\t', line_end - format_start);
        if (unlikely(line_raw_sample_ct != raw_sample_ct)) {
          if (line_raw_sample_ct < raw_sample_ct) {
            goto VcfSubset_ret_MISSING_TOKENS;
          } else {
            fprintf(stderr, "Error: Line %" PRIuPTR " of --in file has more tokens than expected.\n", line_idx);
            goto VcfSubset_ret_MALFORMED_INPUT;
          }
        }
        char* line_crlf = &(line_end[-1]);
        if (*line_crlf == '\r') {
          --line_crlf;
        }
        if (sample_ct == raw_sample_ct) {
          // trivial copy-all case
          write_iter = memcpya(write_iter, info_end, line_crlf - info_end);
        } else {
          char* sample_read_iter = AdvPastDelim(format_start, '\t');
          write_iter = memcpya(write_iter, info_end, sample_read_iter - info_end);
          // This lets us avoid special-casing the last sample.
          *line_crlf = '\t';
          if (sample_ct > raw_sample_ct / 8) {
            // dense case
            // todo: tune threshold
            // todo: try accelerating with SIMD
            for (uint32_t sample_uidx = 0; sample_uidx != raw_sample_ct; ++sample_uidx) {
              char* next_sample_start = AdvPastDelim(sample_read_iter, '\t');
              if (IsSet(sample_include, sample_uidx)) {
                write_iter = memcpya(write_iter, sample_read_iter, next_sample_start - sample_read_iter);
              }
              sample_read_iter = next_sample_start;
            }
          } else {
            // sparse case
            uint32_t prev_sample_uidx_p1 = 0;
            uintptr_t sample_uidx_base = 0;
            uintptr_t cur_bits = sample_include[0];
            for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
              const uint32_t cur_sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
              const uint32_t skip_ct = cur_sample_uidx - prev_sample_uidx_p1;
              if (skip_ct) {
                sample_read_iter = AdvToNthDelim(sample_read_iter, skip_ct, '\t');
                ++sample_read_iter;
              }
              char* next_sample_start = AdvPastDelim(sample_read_iter, '\t');
              write_iter = memcpya(write_iter, sample_read_iter, next_sample_start - sample_read_iter);
              prev_sample_uidx_p1 = cur_sample_uidx + 1;
            }
          }
          --write_iter;
        }
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(bgzfwrite_ck(writebuf_flush, &out_bgzf, &write_iter))) {
        goto VcfSubset_ret_WRITE_FAIL;
      }
      line_start = line_end;
    }
    if (unlikely(TextStreamErrcode2(&in_txs, &reterr))) {
      goto VcfSubset_ret_IN_TSTREAM_FAIL;
    }
    if (unlikely(bgzfclose_flush(writebuf_flush, write_iter, &out_bgzf, &reterr))) {
      goto VcfSubset_ret_1;
    }
    const uintptr_t variant_ct = line_idx - first_variant_line_idx;
    if (raw_sample_ct == 0) {
      snprintf(g_textbuf, kTextbufSize, "vcf_subset: %" PRIuPTR " variant%s written to %s .\n", variant_ct, (variant_ct == 1)? "" : "s", outpath);
    } else {
      snprintf(g_textbuf, kTextbufSize, "vcf_subset: %u/%u sample%s and %" PRIuPTR " variant%s written to %s .\n", sample_ct, raw_sample_ct, (raw_sample_ct == 1)? "" : "s", variant_ct, (variant_ct == 1)? "" : "s", outpath);
    }
    WordWrap(0, g_textbuf);
    fputs(g_textbuf, stdout);
  }
  while (0) {
  VcfSubset_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  VcfSubset_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  VcfSubset_ret_MISSING_TOKENS:
    fprintf(stderr, "Error: Line %" PRIuPTR " of --in file has fewer tokens than expected.\n", line_idx);
    reterr = kPglRetMalformedInput;
    break;
  VcfSubset_ret_MALFORMED_INPUT_WW:
    WordWrap(0, g_textbuf);
    fputs(g_textbuf, stderr);
  VcfSubset_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  VcfSubset_ret_INCONSISTENT_INPUT_WW:
    WordWrap(0, g_textbuf);
    fputs(g_textbuf, stderr);
    reterr = kPglRetInconsistentInput;
    break;
  VcfSubset_ret_IN_TSTREAM_FAIL:
    putc_unlocked('\n', stdout);
    TextErrPrint("--in file", TextStreamError(&in_txs), TextStreamErrcode(&in_txs));
    break;
  }
 VcfSubset_ret_1:
  CleanupBgzfCompressStream(&out_bgzf, &reterr);
  if (unlikely(CleanupTextStream(&in_txs, &reterr))) {
    snprintf(g_textbuf, kTextbufSize, kErrprintfFread, inpath, strerror(errno));
    WordWrap(0, g_textbuf);
    fputs(g_textbuf, stderr);
  }
  return reterr;
}

void PrintExitMsg(PglErr reterr) {
  if (!reterr) {
    return;
  }
  if (reterr == kPglRetNomem) {
    putc_unlocked('\n', stdout);
    fputs(kErrstrNomem, stderr);
    if (g_failed_alloc_attempt_size) {
      fprintf(stderr, "Failed allocation size: %" PRIu64 "\n", g_failed_alloc_attempt_size);
    }
  } else if (reterr == kPglRetWriteFail) {
    putc_unlocked('\n', stdout);
    if (errno) {
      fprintf(stderr, kErrstrWrite, strerror(errno));
    } else {
      // Defensive.
      fputs("Error: File write failure: Untracked cause.", stderr);
    }
  } else if (reterr == kPglRetThreadCreateFail) {
    putc_unlocked('\n', stdout);
    fputs(kErrstrThreadCreate, stderr);
  } else if (reterr == kPglRetLongLine) {
    putc_unlocked('\n', stdout);
    fputs("Error: Unhandled internal line-too-long message.\n", stderr);
  } else if (reterr == kPglRetEof) {
    putc_unlocked('\n', stdout);
    fputs("Error: Unhandled internal EOF message.\n", stderr);
  }
}

#ifndef __LP64__
  // 2047 seems to consistently fail on both OS X and Windows
#  ifdef _WIN32
CONSTI32(kMalloc32bitMibMax, 1792);
#  else
#    ifdef __APPLE__
CONSTI32(kMalloc32bitMibMax, 1952);
#    else
CONSTI32(kMalloc32bitMibMax, 2047);
#    endif
#  endif
#endif

#ifdef __cplusplus
}  // namespace plink2
#endif

int main(int argc, char** argv) {
#ifdef __cplusplus
  using namespace plink2;
#endif
  unsigned char* bigstack_ua = nullptr;
  PglErr reterr = kPglRetSuccess;
  {
    // not actually unlikely, but may as well be consistent re: unlikely() and
    // returning a nonzero error code
    if (unlikely(argc <= 1)) {
      DispUsage(stderr);
      reterr = kPglRetSkipped;
      goto main_ret_1;
    }
    // doesn't use plink2_cmdline; flags may not be processed in alphabetical
    // order, etc.
    const uintptr_t bigstack_min_mib = 640;
    const uintptr_t nonstack_min_mib = 64;
    // 2 GiB for input line buffer, 2 GiB for output line buffer, other stuff
    // should take less than 4 GiB combined (including temporary header load:
    // BCF does not support l_text >= 2^32).
    const uintptr_t max_default_mib = 8192;

    const uint32_t argc_u = S_CAST(uint32_t, argc);
    // These just point to argv when initialized, so they must not be freed on
    // exit.
    char* inpath = nullptr;
    char* indv_str = nullptr;
    char* keeppath = nullptr;
    char* outpath = nullptr;

    uint32_t allow_any_paths = 0;
    uint32_t allow_no_samples = 0;
    uint32_t bgz_level = 0;
    uint32_t erase_info = 0;
    uintptr_t malloc_size_mib = 0;
    uint32_t no_samples = 0;
    uint32_t thread_ct = 0;
    for (uint32_t arg_idx = 1; arg_idx < argc_u; ++arg_idx) {
      const char* cur_arg = argv[arg_idx];
      if (unlikely((cur_arg[0] != '-') || ((cur_arg[1] == '-') && (cur_arg[2] == '-')))) {
        snprintf(g_textbuf, kTextbufSize, "Error: '%s' is not a flag (should start with one or two dashes).\n", cur_arg);
        goto main_ret_INVALID_CMDLINE_WW;
      }
      const char* flagname_p = &(cur_arg[1 + (cur_arg[1] == '-')]);
      const uint32_t flagname_slen = strlen(flagname_p);
      if (strequal_k(flagname_p, "allow-any-paths", flagname_slen)) {
        allow_any_paths = 1;
        continue;
      }
      if (strequal_k(flagname_p, "allow-no-samples", flagname_slen)) {
        allow_no_samples = 1;
        continue;
      }
      if (strequal_k(flagname_p, "bgz-level", flagname_slen)) {
        if (unlikely(bgz_level)) {
          fputs("Error: Multiple instances of --bgz-level.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        ++arg_idx;
        if (unlikely(arg_idx == argc_u)) {
          fputs("Error: Missing --bgz-level parameter.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        if (unlikely(ScanPosintDefcapx(argv[arg_idx], &bgz_level) || (bgz_level > 12))) {
          snprintf(g_textbuf, kTextbufSize, "Error: Invalid --bgz-level argument '%s' (must be an integer 1-12).\n", argv[arg_idx]);
          goto main_ret_INVALID_CMDLINE_WW;
        }
        continue;
      }
      if (strequal_k(flagname_p, "erase-info", flagname_slen)) {
        // doesn't matter if there are multiple instances of this flag
        erase_info = 1;
        continue;
      }
      if (strequal_k(flagname_p, "help", flagname_slen)) {
        if (argc_u > 2) {
          // No need for plink-style keyword-based subsetting.
          fputs("Note: --help present; ignoring rest of command line.\n", stdout);
        }
        DispUsage(stdout);
        goto main_ret_1;
      }
      if (strequal_k(flagname_p, "in", flagname_slen)) {
        if (unlikely(inpath)) {
          fputs("Error: Multiple instances of --in.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        ++arg_idx;
        if (unlikely(arg_idx == argc_u)) {
          fputs("Error: Missing --in parameter.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        inpath = argv[arg_idx];
        continue;
      }
      if (strequal_k(flagname_p, "indv", flagname_slen)) {
        if (unlikely(indv_str)) {
          fputs("Error: Multiple instances of --indv.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        ++arg_idx;
        if (unlikely(arg_idx == argc_u)) {
          fputs("Error: Missing --indv parameter.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        indv_str = argv[arg_idx];
        continue;
      }
      if (strequal_k(flagname_p, "keep", flagname_slen)) {
        if (unlikely(keeppath)) {
          fputs("Error: Multiple instances of --keep.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        ++arg_idx;
        if (unlikely(arg_idx == argc_u)) {
          fputs("Error: Missing --keep parameter.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        keeppath = argv[arg_idx];
        continue;
      }
      if (strequal_k(flagname_p, "memory", flagname_slen)) {
        if (unlikely(malloc_size_mib)) {
          fputs("Error: Multiple instances of --memory.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        ++arg_idx;
        if (unlikely(arg_idx == argc_u)) {
          fputs("Error: Missing --memory parameter.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        if (unlikely(ScanPosintptrx(argv[arg_idx], &malloc_size_mib))) {
          snprintf(g_textbuf, kTextbufSize, "Error: Invalid --memory argument '%s'.\n", argv[arg_idx]);
          goto main_ret_INVALID_CMDLINE_WW;
        }
        if (unlikely(malloc_size_mib > bigstack_min_mib)) {
          snprintf(g_textbuf, kTextbufSize, "Error: Invalid --memory argument '%s' (minimum %" PRIuPTR ").\n", argv[arg_idx], bigstack_min_mib);
          goto main_ret_INVALID_CMDLINE_WW;
        }
#ifndef __LP64__
        if (unlikely(malloc_size_mib > kMalloc32bitMibMax)) {
          fprintf(fprintf, "Error: --memory argument too large for 32-bit version (max %u).\n", kMalloc32bitMibMax);
          goto main_ret_INVALID_CMDLINE;
        }
#endif
        continue;
      }
      if (strequal_k(flagname_p, "no-samples", flagname_slen)) {
        no_samples = 1;
        continue;
      }
      if (strequal_k(flagname_p, "out", flagname_slen)) {
        if (unlikely(outpath)) {
          fputs("Error: Multiple instances of --out.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        ++arg_idx;
        if (unlikely(arg_idx == argc_u)) {
          fputs("Error: Missing --out parameter.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        outpath = argv[arg_idx];
        continue;
      }
      if (strequal_k(flagname_p, "threads", flagname_slen)) {
        if (unlikely(thread_ct)) {
          fputs("Error: Multiple instances of --threads.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        ++arg_idx;
        if (unlikely(arg_idx == argc_u)) {
          fputs("Error: Missing --threads parameter.\n", stderr);
          goto main_ret_INVALID_CMDLINE;
        }
        if (unlikely(ScanPosintDefcapx(argv[arg_idx], &thread_ct))) {
          snprintf(g_textbuf, kTextbufSize, "Error: Invalid --threads argument '%s'.\n", argv[arg_idx]);
          goto main_ret_INVALID_CMDLINE_WW;
        }
        continue;
      }
      if (likely(strequal_k(flagname_p, "version", flagname_slen))) {
        if (argc_u > 2) {
          fputs("Warning: --version present; ignoring rest of command line.\n", stderr);
        }
        fputs(ver_str, stdout);
        putc_unlocked('\n', stdout);
        goto main_ret_1;
      }
      snprintf(g_textbuf, kTextbufSize, "Error: Unrecognized flag '%s'.\n", argv[arg_idx]);
      goto main_ret_INVALID_CMDLINE_WW;
    }  // end of initial command-line parsing loop
    PrintVer();

    // 1. require --in and --out
    // 2. enforce path restrictions
    // 3. count cores, allocate workspace
    // 4. open --in file, load all ## header lines
    // 5. count number of samples, apply --indv or --keep
    // 6. open --out file, process remaining lines
    if ((!inpath) || (!outpath)) {
      fputs("Error: --in and --out are required.\n", stderr);
      goto main_ret_INVALID_CMDLINE;
    }
    if (indv_str && keeppath) {
      fputs("Error: --indv and --keep cannot be used simultaneously.\n", stderr);
      goto main_ret_INVALID_CMDLINE;
    }

    const uint32_t inpath_slen = strlen(inpath);
    const uint32_t outpath_slen = strlen(outpath);
    if (!allow_any_paths) {
      if (unlikely((!StrEndsWith(inpath, ".vcf", inpath_slen)) &&
                   (!StrEndsWith(inpath, ".vcf.gz", inpath_slen)) &&
                   (!StrEndsWith(inpath, ".vcf.bgz", inpath_slen)))) {
        fputs("Error: --in argument does not end with an expected file extension (.vcf,\n.vcf.gz, .vcf.bgz).  You can override this sanity check with --allow-any-paths.\n", stderr);
        goto main_ret_INVALID_CMDLINE;
      }
      if (unlikely((!StrEndsWith(outpath, ".vcf", outpath_slen)) &&
                   (!StrEndsWith(outpath, ".vcf.gz", outpath_slen)) &&
                   (!StrEndsWith(outpath, ".vcf.bgz", outpath_slen)))) {
        fputs("Error: --out argument does not end with an expected file extension (.vcf,\n.vcf.gz, .vcf.bgz).  You can override this sanity check with --allow-any-paths.\n", stderr);
        goto main_ret_INVALID_CMDLINE;
      }
    }

    if (!thread_ct) {
      thread_ct = NumCpu(nullptr);
      if (thread_ct > kMaxThreads) {
        thread_ct = kMaxThreads;
      }
    }
    if (thread_ct > 8) {
      printf("Using up to %u threads (change this with --threads).\n", thread_ct);
    } else {
      printf("Using %s%u compute thread%s.\n", (thread_ct > 1)? "up to " : "", thread_ct, (thread_ct == 1)? "" : "s");
    }

    // todo: move most/all of this logic from CmdlineParsePhase3() to helper(s)
    // in plink2_memory
    const uint64_t total_mib = DetectMib();
    if (!malloc_size_mib) {
      if (!total_mib) {
        malloc_size_mib = max_default_mib;
      } else if (total_mib < (bigstack_min_mib * 2)) {
        malloc_size_mib = bigstack_min_mib;
      } else {
#ifdef __LP64__
        malloc_size_mib = total_mib / 2;
#else
        malloc_size_mib = MINV(2047, total_mib / 2);
#endif
        if (malloc_size_mib > max_default_mib) {
          malloc_size_mib = max_default_mib;
        }
      }
    }
    assert(malloc_size_mib > bigstack_min_mib);
#ifndef __LP64__
    if (malloc_size_mib > kMalloc32bitMibMax) {
      malloc_size_mib = kMalloc32bitMibMax;
    }
#endif
    if (total_mib) {
      const uint64_t mem_available_kib = GetMemAvailableKib(kMaxMediumLine, g_textbuf);
      if (mem_available_kib == (~0LLU)) {
        printf("%" PRIu64 " MiB RAM detected; reserving %" PRIuPTR " MiB for main workspace.\n", total_mib, malloc_size_mib);
      } else {
        const uint64_t mem_available_mib = mem_available_kib / 1024;
        if (mem_available_mib < malloc_size_mib + nonstack_min_mib) {
          if (mem_available_mib < bigstack_min_mib + nonstack_min_mib) {
            malloc_size_mib = bigstack_min_mib;
          } else {
            malloc_size_mib = mem_available_mib - nonstack_min_mib;
          }
        }
        printf("%" PRIu64 " MiB RAM detected, ~%" PRIu64 " available; reserving %" PRIuPTR " MiB for main workspace.\n", total_mib, mem_available_mib, malloc_size_mib);
      }
    } else {
      printf("Failed to determine total system memory.  Attempting to reserve %" PRIuPTR " MiB.\n", malloc_size_mib);
    }

    // guarantee contiguous address space outside of main workspace
    // ...though this is pointless with overcommit on, may want to
    // conditionally compile this out, and/or conditionally skip this
    unsigned char* bigstack_base;
    unsigned char* bigstack_end;
    uintptr_t malloc_mib_final = malloc_size_mib;
    {
      unsigned char* bubble;
      if (unlikely(pgl_malloc(nonstack_min_mib << 20, &bubble))) {
        goto main_ret_NOMEM;
      }

      // don't use pgl_malloc here since we don't automatically want to set
      // g_failed_alloc_attempt_size on failure
      bigstack_ua = S_CAST(unsigned char*, malloc(malloc_mib_final << 20));
      while (!bigstack_ua) {
        malloc_mib_final = (malloc_mib_final * 3) / 4;
        if (malloc_mib_final < bigstack_min_mib) {
          malloc_mib_final = bigstack_min_mib;
        }
        bigstack_ua = S_CAST(unsigned char*, malloc(malloc_mib_final << 20));
        if (unlikely((!bigstack_ua) && (malloc_mib_final == bigstack_min_mib))) {
          g_failed_alloc_attempt_size = bigstack_min_mib << 20;
          free(bubble);
          goto main_ret_NOMEM;
        }
      }
      // force 64-byte align to make cache line sensitivity work
      unsigned char* bigstack_initial_base = R_CAST(unsigned char*, RoundUpPow2(R_CAST(uintptr_t, bigstack_ua), kCacheline));
      bigstack_base = bigstack_initial_base;
      // no g_one_char_strs here, just have a 64 byte overread buffer
      bigstack_end = &(bigstack_initial_base[RoundDownPow2((malloc_mib_final << 20) - 64 - S_CAST(uintptr_t, bigstack_initial_base - bigstack_ua), kCacheline)]);
      free(bubble);
    }
    if (malloc_size_mib != malloc_mib_final) {
      printf("Allocated %" PRIuPTR " MiB successfully, after larger attempt(s) failed.\n", malloc_mib_final);
    }

    if (!bgz_level) {
      bgz_level = kBgzfDefaultClvl;
    }
    reterr = VcfSubset(bigstack_base, bigstack_end, inpath, indv_str, keeppath, outpath, outpath_slen, allow_no_samples, bgz_level, erase_info, no_samples, thread_ct);
    if (unlikely(reterr)) {
      goto main_ret_1;
    }
  }
  while (0) {
  main_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  main_ret_INVALID_CMDLINE_WW:
    WordWrap(0, g_textbuf);
    fputs(g_textbuf, stderr);
  main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
 main_ret_1:
  PrintExitMsg(reterr);
  free_cond(bigstack_ua);
  return S_CAST(int32_t, reterr);
}
