#include "include/pgenlib_misc.h"
#include "include/pgenlib_read.h"
#include "include/pgenlib_write.h"
#include "include/plink2_base.h"
#include "include/plink2_bits.h"

// #define SUBSET_TEST

int32_t main(int32_t argc, char** argv) {
#ifdef __cplusplus
  using namespace plink2;
#endif
  PglErr reterr = kPglRetSuccess;
  unsigned char* pgfi_alloc = nullptr;
  unsigned char* pgr_alloc = nullptr;
  unsigned char* spgw_alloc = nullptr;
  uintptr_t* genovec = nullptr;
  uintptr_t* raregeno = nullptr;
  uintptr_t* sample_include = nullptr;
  uint32_t* sample_include_cumulative_popcounts = nullptr;
  uint32_t* difflist_sample_ids = nullptr;
  FILE* outfile = nullptr;
  PgenHeaderCtrl header_ctrl;
  STPgenWriter spgw;
  uint32_t write_sample_ct;
  PgenFileInfo pgfi;
  PgenReader pgr;
  PreinitPgfi(&pgfi);
  PreinitPgr(&pgr);
  PreinitSpgw(&spgw);
  {
    if ((argc < 3) || (argc > 6)) {
      fputs(
"Usage:\n"
"pgen_compress [-i] <input .bed or .pgen> <output filename> [sample_ct]\n"
"  * sample_ct is required when loading a .bed file\n"
"  * -i causes the index to be saved to a separate .pgen.pgi file (usually it's\n"
"    embedded at the front of the .pgen); this is compatible with fully\n"
"    sequential .pgen writing\n"
"pgen_compress -u <input .pgen> <output .bed>\n"
            , stdout);
      goto main_ret_INVALID_CMDLINE;
    }
    const uint32_t write_separate_index = (argv[1][0] == '-') && (argv[1][1] == 'i') && (argv[1][2] == '\0');
    const uint32_t decompress = (argv[1][0] == '-') && (argv[1][1] == 'u') && (argv[1][2] == '\0');
    const uint32_t input_idx = 1 + write_separate_index + decompress;
    uint32_t sample_ct = 0xffffffffU;
    if (S_CAST(uint32_t, argc) == input_idx + 3) {
      if (ScanPosintDefcap(argv[input_idx + 2], &sample_ct)) {
        fprintf(stderr, "error: invalid sample_ct\n");
        goto main_ret_INVALID_CMDLINE;
      }
    }
    char errstr_buf[kPglErrstrBufBlen];
    uintptr_t cur_alloc_cacheline_ct;
    reterr = PgfiInitPhase1(argv[input_idx], nullptr, 0xffffffffU, sample_ct, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, errstr_buf);
    if (reterr) {
      fputs(errstr_buf, stderr);
      goto main_ret_1;
    }
    sample_ct = pgfi.raw_sample_ct;
    if (!sample_ct) {
      fprintf(stderr, "error: sample_ct == 0\n");
      goto main_ret_INVALID_CMDLINE;
    }
    const uint32_t variant_ct = pgfi.raw_variant_ct;
    if (!variant_ct) {
      fprintf(stderr, "error: variant_ct == 0\n");
      goto main_ret_INVALID_CMDLINE;
    }
    if (cachealigned_malloc(cur_alloc_cacheline_ct * kCacheline, &pgfi_alloc)) {
      goto main_ret_NOMEM;
    }
    uint32_t max_vrec_width;
    // todo: test block-fread
    reterr = PgfiInitPhase2(header_ctrl, 0, 0, 0, 0, variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, &cur_alloc_cacheline_ct, errstr_buf);
    if (reterr) {
      fputs(errstr_buf, stderr);
      goto main_ret_1;
    }
    if (cachealigned_malloc(cur_alloc_cacheline_ct * kCacheline, &pgr_alloc)) {
      goto main_ret_NOMEM;
    }

    // modify this when trying block-fread
    reterr = PgrInit(argv[input_idx], max_vrec_width, &pgfi, &pgr, pgr_alloc);
    if (reterr) {
      fprintf(stderr, "PgrInit error %u\n", S_CAST(uint32_t, reterr));
      goto main_ret_1;
    }

    if (S_CAST(uint32_t, argc) == input_idx + 3) {
      printf("%u variant%s detected.\n", variant_ct, (variant_ct == 1)? "" : "s");
    } else {
      printf("%u variant%s and %u sample%s detected.\n", variant_ct, (variant_ct == 1)? "" : "s", sample_ct, (sample_ct == 1)? "" : "s");
    }
    if (cachealigned_malloc(NypCtToVecCt(sample_ct) * kBytesPerVec, &genovec)) {
      goto main_ret_NOMEM;
    }
    if (decompress) {
      outfile = fopen(argv[3], FOPEN_WB);
      if (!outfile) {
        goto main_ret_OPEN_FAIL;
      }
      PgrSampleSubsetIndex pssi;
      PgrClearSampleSubsetIndex(&pgr, &pssi);
      const uintptr_t final_mask = (k1LU << ((sample_ct % kBitsPerWordD2) * 2)) - k1LU;
      const uint32_t final_widx = NypCtToWordCt(sample_ct) - 1;
      const uint32_t variant_byte_ct = (sample_ct + 3) / 4;
      fwrite("l\x1b\x01", 3, 1, outfile);
      for (uint32_t vidx = 0; vidx < variant_ct; ) {
        reterr = PgrGet(nullptr, pssi, sample_ct, vidx, &pgr, genovec);
        if (reterr) {
          fprintf(stderr, "\nread error %u, vidx=%u\n", S_CAST(uint32_t, reterr), vidx);
          goto main_ret_1;
        }
        PgrPlink2ToPlink1InplaceUnsafe(sample_ct, genovec);
        if (final_mask) {
          genovec[final_widx] &= final_mask;
        }
        fwrite(genovec, variant_byte_ct, 1, outfile);
        ++vidx;
        if (!(vidx % 100000)) {
          printf("\r%u.%um variants decompressed.", vidx / 1000000, (vidx / 100000) % 10);
          fflush(stdout);
        }
      }
      if (fclose_null(&outfile)) {
        goto main_ret_WRITE_FAIL;
      }
      printf("\n");
      goto main_ret_1;
    }
#ifdef SUBSET_TEST
    // write_sample_ct = sample_ct - 3;
    write_sample_ct = 3;
#else
    write_sample_ct = sample_ct;
#endif
    // Demonstrate that, when write_separate_index is true, variant_ct_limit
    // can be an overestimate.
    // Also demonstrate automatic 8-bit -> 4-bit index compaction when
    // write_separate_index is true, we declare that hardcall-phase is present,
    // but we never write any hardcall-phase data.
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(argv[input_idx + 1], nullptr, nullptr, write_separate_index? (variant_ct * 2) : variant_ct, write_sample_ct, 0, write_separate_index? kPgenWriteSeparateIndex : kPgenWriteBackwardSeek, write_separate_index? kfPgenGlobalHardcallPhasePresent : kfPgenGlobal0, 2, &spgw, &cur_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      fprintf(stderr, "compression phase 1 error %u\n", S_CAST(uint32_t, reterr));
      goto main_ret_1;
    }
    if (cachealigned_malloc(cur_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto main_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t max_read_difflist_len = 2 * (sample_ct / kPglMaxDifflistLenDivisor);
    const uint32_t max_write_difflist_len = 2 * (write_sample_ct / kPglMaxDifflistLenDivisor);
    if (cachealigned_malloc(RoundUpPow2((max_read_difflist_len + 3) / 4, kCacheline), &raregeno) ||
        cachealigned_malloc(RoundUpPow2((sample_ct + 7) / 8, kCacheline), &sample_include) ||
        cachealigned_malloc(RoundUpPow2((1 + (sample_ct / kBitsPerWord)) * sizeof(int32_t), kCacheline), &sample_include_cumulative_popcounts) ||
        cachealigned_malloc(RoundUpPow2((max_read_difflist_len + 1) * sizeof(int32_t), kCacheline), &difflist_sample_ids)) {
      goto main_ret_NOMEM;
    }
#ifdef SUBSET_TEST
    fill_ulong_zero(BitCtToWordCt(sample_ct), sample_include);
    SetBit(123, sample_include);
    SetBit(127, sample_include);
    SetBit(320, sample_include);
    // ClearBit(123, sample_include);
    // ClearBit(127, sample_include);
    // ClearBit(320, sample_include);
#else
    SetAllBits(sample_ct, sample_include);
#endif
    FillCumulativePopcounts(sample_include, 1 + (sample_ct / kBitsPerWord), sample_include_cumulative_popcounts);
    uint32_t difflist_common_geno = 0;  // spurious uninit-variable warning
    PgrSampleSubsetIndex pssi;
    PgrSetSampleSubsetIndex(sample_include_cumulative_popcounts, &pgr, &pssi);
    for (uint32_t vidx = 0; vidx < variant_ct; ) {
      uint32_t difflist_len;
      reterr = PgrGetDifflistOrGenovec(sample_include, pssi, write_sample_ct, max_read_difflist_len, vidx, &pgr, genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
      if (reterr) {
        fprintf(stderr, "\nread error %u, vidx=%u\n", S_CAST(uint32_t, reterr), vidx);
        goto main_ret_1;
      }
      if (difflist_common_geno == 0xffffffffU) {
        ZeroTrailingBits(write_sample_ct * 2, genovec);
        reterr = SpgwAppendBiallelicGenovec(genovec, &spgw);
      } else if (difflist_len <= max_write_difflist_len) {
        ZeroTrailingBits(2 * difflist_len, raregeno);
        difflist_sample_ids[difflist_len] = write_sample_ct;
        reterr = SpgwAppendBiallelicDifflistLimited(raregeno, difflist_sample_ids, difflist_common_geno, difflist_len, &spgw);
      } else {
        PgrDifflistToGenovecUnsafe(raregeno, difflist_sample_ids, difflist_common_geno, write_sample_ct, difflist_len, genovec);
        ZeroTrailingBits(write_sample_ct * 2, genovec);
        reterr = SpgwAppendBiallelicGenovec(genovec, &spgw);
      }
      if (reterr) {
        fprintf(stderr, "\ncompress/write error %u, vidx=%u\n", S_CAST(uint32_t, reterr), vidx);
        goto main_ret_1;
      }
      ++vidx;
      if (!(vidx % 100000)) {
        printf("\r%u.%um variants compressed.", vidx / 1000000, (vidx / 100000) % 10);
        fflush(stdout);
      }
    }
  }
  printf("\n");

  reterr = SpgwFinish(&spgw);
  while (0) {
  main_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  main_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  main_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  main_ret_INVALID_CMDLINE:
    reterr = kPglRetInvalidCmdline;
    break;
  }
 main_ret_1:
  CleanupPgr(&pgr, &reterr);
  CleanupPgfi(&pgfi, &reterr);
  CleanupSpgw(&spgw, &reterr);
  if (pgfi_alloc) {
    aligned_free(pgfi_alloc);
  }
  if (pgr_alloc) {
    aligned_free(pgr_alloc);
  }
  if (spgw_alloc) {
    aligned_free(spgw_alloc);
  }
  if (genovec) {
    aligned_free(genovec);
  }
  if (raregeno) {
    aligned_free(raregeno);
  }
  if (sample_include) {
    aligned_free(sample_include);
  }
  if (sample_include_cumulative_popcounts) {
    aligned_free(sample_include_cumulative_popcounts);
  }
  if (difflist_sample_ids) {
    aligned_free(difflist_sample_ids);
  }
  if (outfile) {
    fclose(outfile);
  }
  return S_CAST(int32_t, reterr);
}
