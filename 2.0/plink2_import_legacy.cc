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


#include "include/pgenlib_write.h"
#include "plink2_compress_stream.h"
#include "plink2_import.h"
#include "plink2_psam.h"

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr PedmapToPgen(__attribute__((unused)) const char* pedname, __attribute__((unused)) const char* mapname, __attribute__((unused)) MiscFlags misc_flags, __attribute__((unused)) ImportFlags import_flags, __attribute__((unused)) uint32_t psam_01, __attribute__((unused)) FamCol fam_cols, __attribute__((unused)) int32_t missing_pheno, __attribute__((unused)) uint32_t max_thread_ct, __attribute__((unused)) char* outname, __attribute__((unused)) char* outname_end, __attribute__((unused)) ChrInfo* cip) {
  logerrputs("Error: .ped import is under development.\n");
  return kPglRetNotYetSupported;
}

PglErr RewritePsam(const char* in_psamname, MiscFlags misc_flags, FamCol fam_cols, int32_t missing_pheno, uint32_t max_thread_ct, char* outname, char* outname_end, uint32_t* raw_sample_ctp) {
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
    reterr = LoadPsam(in_psamname, nullptr, fam_cols, 0x7fffffff, missing_pheno, (misc_flags / kfMiscAffection01) & 1, max_thread_ct, &pii, &sample_include, &founder_info, &sex_nm, &sex_male, &pheno_cols, &pheno_names, &raw_sample_ct, &pheno_ct, &max_pheno_name_blen);
    if (unlikely(reterr)) {
      goto RewritePsam_ret_1;
    }
    *raw_sample_ctp = raw_sample_ct;

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    // Note that --output-missing-phenotype doesn't apply to autoconversion.
    reterr = WritePsam(outname, sample_include, &(pii.sii), (&pii.parental_id_info), sex_nm, sex_male, pheno_cols, pheno_names, nullptr, "NA", raw_sample_ct, pheno_ct, max_pheno_name_blen, kfPsamColDefault);
  }
 RewritePsam_ret_1:
  CleanupPhenoCols(pheno_ct, pheno_cols);
  free_cond(pheno_names);
  BigstackReset(bigstack_mark);
  return reterr;
}

// psam_generated assumed to be initialized to 1.
// Unlike plink 1.9, this does not support lines longer than 2 GiB.
PglErr TpedToPgen(const char* tpedname, const char* tfamname, MiscFlags misc_flags, ImportFlags import_flags, FamCol fam_cols, int32_t missing_pheno, uint32_t max_thread_ct, char* outname, char* outname_end, ChrInfo* cip, uint32_t* psam_generated_ptr) {
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
      reterr = RewritePsam(tfamname, misc_flags, fam_cols, missing_pheno, max_thread_ct, outname, outname_end, &tfam_sample_ct);
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
    // possible todo: benchmark against generating a .bed temporary file.
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
    FinalizeChrset(misc_flags, cip);
    const uint32_t allow_extra_chrs = (misc_flags / kfMiscAllowExtraChrs) & 1;
    uint32_t max_line_blen = TextLineEnd(&tped_txs) - tped_line_start;
    uint32_t variant_ct = 0;
    uint32_t at_least_one_nzero_cm = 0;
    while (1) {
      char* chr_code_end = CurTokenEnd(tped_line_start);

      // must do this before chromosome-code null-termination
      char* cm_start = FirstNonTspace(FirstSpaceOrEoln(FirstNonTspace(chr_code_end)));
      if (unlikely(IsEolnKns(*cm_start))) {
        goto TpedToPgen_ret_MISSING_TOKENS;
      }

      uint32_t cur_chr_code;
      reterr = GetOrAddChrCodeDestructive(".tped file", line_idx, allow_extra_chrs, tped_line_start, chr_code_end, cip, &cur_chr_code);
      if (unlikely(reterr)) {
        goto TpedToPgen_ret_1;
      }
      if (IsSet(cip->chr_mask, cur_chr_code)) {
        if (variant_ct == 0x7ffffffd) {
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
    if (unlikely(!variant_ct)) {
      logerrputs("Error: All .tped variants excluded by chromosome filter.\n");
      goto TpedToPgen_ret_INVALID_CMDLINE;
    }
    if (unlikely(CleanupTextStream2(tpedname, &tped_txs, &reterr))) {
      goto TpedToPgen_ret_1;
    }
    BigstackReset(bigstack_mark);
    const uint32_t variant_skip_ct = line_idx - 1 - variant_ct;
    if (variant_skip_ct) {
      logprintfww("--tped: %" PRIuPTR " variants scanned; %u excluded by chromosome filter, %u remaining.\n", line_idx - 1, variant_skip_ct, variant_ct);
    } else {
      logprintf("--tped: %u variant%s scanned.\n", variant_ct, (variant_ct == 1)? "" : "s");
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
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, 0, kfPgenGlobal0, 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
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
    uint32_t next_print_idx = variant_ct / 100;
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
      pvar_cswritep = chrtoa(cip, cur_chr_code, pvar_cswritep);

      char* variant_id_start = FirstNonTspace(chr_code_end);
      char* variant_id_end = FirstSpaceOrEoln(variant_id_start);
      char* cm_start = FirstNonTspace(variant_id_end);
      char* cm_end = FirstSpaceOrEoln(cm_start);
      char* bp_start = FirstNonTspace(cm_end);
      if (unlikely(IsEolnKns(*bp_start))) {
        goto TpedToPgen_ret_MISSING_TOKENS;
      }
      uint32_t cur_bp;
      if (unlikely(ScanUintDefcap(bp_start, &cur_bp))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid bp coordinate on line %" PRIuPTR " of %s.\n", line_idx, tpedname);
        goto TpedToPgen_ret_MALFORMED_INPUT_WW;
      }
      *pvar_cswritep++ = '\t';
      pvar_cswritep = u32toa_x(cur_bp, '\t', pvar_cswritep);
      pvar_cswritep = memcpyax(pvar_cswritep, variant_id_start, variant_id_end - variant_id_start, '\t');
      tped_line_iter = FirstNonTspace(CurTokenEnd(bp_start));

      ZeroWArr(sample_ctl2, genovec);
      const char input_missing_geno_char = g_input_missing_geno_ptr[0];
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
        dtoa_g_p8(cur_cm, pvar_cswritep);
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
        next_print_idx = (pct * S_CAST(uint64_t, variant_ct)) / 100;
      }
      tped_line_iter = AdvPastDelim(tped_line_iter, '\n');
    }
    SpgwFinish(&spgw);
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
    ForgetExtraChrNames(1, cip);
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
  TpedToPgen_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  TpedToPgen_ret_HALF_MISSING:
    logerrprintfww("Error: Half-missing genotype on line %" PRIuPTR " of %s.\n", line_idx, tpedname);
    reterr = kPglRetMalformedInput;
    break;
  TpedToPgen_ret_MULTIALLELIC:
    logerrprintfww("Error: Multiallelic variant on line %" PRIuPTR " of %s. This violates the .tped specification; please reformat this as e.g. VCF.\n", line_idx, tpedname);
    reterr = kPglRetMalformedInput;
    break;
  TpedToPgen_ret_MISSING_TOKENS:
    snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of %s has fewer tokens than expected.\n", line_idx, tpedname);
  TpedToPgen_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  TpedToPgen_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
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

#ifdef __cplusplus
}  // namespace plink2
#endif
