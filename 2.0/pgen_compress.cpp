#include "pgenlib_internal.h"

// #define SUBSET_TEST

int32_t main(int32_t argc, char** argv) {
#ifdef __cplusplus
  using namespace plink2;
#endif
  pglerr_t reterr = kPglRetSuccess;
  unsigned char* pgfi_alloc = nullptr;
  unsigned char* pgr_alloc = nullptr;
  unsigned char* spgw_alloc = nullptr;
  uintptr_t* genovec = nullptr;
  uintptr_t* raregeno = nullptr;
  uintptr_t* sample_include = nullptr;
  uint32_t* sample_include_cumulative_popcounts = nullptr;
  uint32_t* difflist_sample_ids = nullptr;
  FILE* outfile = nullptr;
  pgen_header_ctrl_t header_ctrl;
  st_pgen_writer_t spgw;
  uint32_t write_sample_ct;
  pgen_file_info_t pgfi;
  pgen_reader_t pgr;
  pgfi_preinit(&pgfi);
  pgr_preinit(&pgr);
  spgw_preinit(&spgw);
  {
    const uint32_t use_mmap = 0;
    if ((argc < 3) || (argc > 5)) {
      fputs(
"Usage:\n"
"pgen_compress [input .bed or .pgen] [output filename] {sample_ct}\n"
"  (sample_ct is required when loading a .bed file)\n"
"pgen_compress -u [input .pgen] [output .bed]\n"
	    , stdout);
      goto main_ret_INVALID_CMDLINE;
    }
    const uint32_t decompress = (argv[1][0] == '-') && (argv[1][1] == 'u') && (argv[1][2] == '\0');
    uint32_t sample_ct = 0xffffffffU;
    if (((uint32_t)argc) == 4 + decompress) {
      if (scan_posint_defcap(argv[3 + decompress], &sample_ct)) {
	goto main_ret_INVALID_CMDLINE;
      }
    }
    char errstr_buf[kPglErrstrBufBlen];
    uintptr_t cur_alloc_cacheline_ct;
    reterr = pgfi_init_phase1(argv[1 + decompress], 0xffffffffU, sample_ct, use_mmap, &header_ctrl, &pgfi, &cur_alloc_cacheline_ct, errstr_buf);
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
    reterr = pgfi_init_phase2(header_ctrl, 0, 0, 0, 0, variant_ct, &max_vrec_width, &pgfi, pgfi_alloc, &cur_alloc_cacheline_ct, errstr_buf);
    if (reterr) {
      fputs(errstr_buf, stderr);
      goto main_ret_1;
    }
    if (cachealigned_malloc(cur_alloc_cacheline_ct * kCacheline, &pgr_alloc)) {
      goto main_ret_NOMEM;
    }

    // modify this when trying block-fread
    reterr = pgr_init(use_mmap? nullptr : argv[1 + decompress], max_vrec_width, &pgfi, &pgr, pgr_alloc);
    if (reterr) {
      fprintf(stderr, "pgr_init error %u\n", (uint32_t)reterr);
      goto main_ret_1;
    }

    if (((uint32_t)argc) == 4 + decompress) {
      printf("%u variant%s detected.\n", variant_ct, (variant_ct == 1)? "" : "s");
    } else {
      printf("%u variant%s and %u sample%s detected.\n", variant_ct, (variant_ct == 1)? "" : "s", sample_ct, (sample_ct == 1)? "" : "s");
    }
    if (cachealigned_malloc(QUATERCT_TO_VECCT(sample_ct) * kBytesPerVec, &genovec)) {
      goto main_ret_NOMEM;
    }
    if (decompress) {
      outfile = fopen(argv[3], FOPEN_WB);
      if (!outfile) {
	goto main_ret_OPEN_FAIL;
      }
      const uintptr_t final_mask = (k1LU << ((sample_ct % kBitsPerWordD2) * 2)) - k1LU;
      const uint32_t final_widx = QUATERCT_TO_WORDCT(sample_ct) - 1;
      const uint32_t variant_byte_ct = (sample_ct + 3) / 4;
      fwrite("l\x1b\x01", 3, 1, outfile);
      for (uint32_t vidx = 0; vidx < variant_ct;) {
	reterr = pgr_read_refalt1_genovec_subset_unsafe(nullptr, nullptr, sample_ct, vidx, &pgr, genovec);
	if (reterr) {
	  fprintf(stderr, "\nread error %u, vidx=%u\n", (uint32_t)reterr, vidx);
	  goto main_ret_1;
	}
	pgr_plink2_to_plink1_inplace_unsafe(sample_ct, genovec);
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
    uint32_t max_vrec_len;
    reterr = spgw_init_phase1(argv[2], nullptr, nullptr, variant_ct, write_sample_ct, kfPgenGlobal0, 2, &spgw, &cur_alloc_cacheline_ct, &max_vrec_len);
    if (reterr) {
      fprintf(stderr, "compression phase 1 error %u\n", (uint32_t)reterr);
      goto main_ret_1;
    }
    if (cachealigned_malloc(cur_alloc_cacheline_ct * kCacheline, &spgw_alloc)) {
      goto main_ret_NOMEM;
    }
    spgw_init_phase2(max_vrec_len, &spgw, spgw_alloc);

    const uint32_t max_simple_difflist_len = sample_ct / kBitsPerWordD2;
    const uint32_t max_returned_difflist_len = max_simple_difflist_len + (sample_ct / kPglMaxDifflistLenDivisor);
    const uint32_t max_difflist_len = 2 * (write_sample_ct / kPglMaxDifflistLenDivisor);
    if (cachealigned_malloc(round_up_pow2((max_returned_difflist_len + 3) / 4, kCacheline), &raregeno) ||
        cachealigned_malloc(round_up_pow2((sample_ct + 7) / 8, kCacheline), &sample_include) ||
	cachealigned_malloc(round_up_pow2((1 + (sample_ct / kBitsPerWord)) * sizeof(int32_t), kCacheline), &sample_include_cumulative_popcounts) ||
	cachealigned_malloc(round_up_pow2((max_returned_difflist_len + 1) * sizeof(int32_t), kCacheline), &difflist_sample_ids)) {
      goto main_ret_NOMEM;
    }
    fill_all_bits(sample_ct, sample_include);
#ifdef SUBSET_TEST
    fill_ulong_zero(BITCT_TO_WORDCT(sample_ct), sample_include);
    set_bit(123, sample_include);
    set_bit(127, sample_include);
    set_bit(320, sample_include);
    // clear_bit(123, sample_include);
    // clear_bit(127, sample_include);
    // clear_bit(320, sample_include);
    fill_cumulative_popcounts(sample_include, 1 + (sample_ct / kBitsPerWord), sample_include_cumulative_popcounts);
#endif
    for (uint32_t vidx = 0; vidx < variant_ct;) {
      uint32_t difflist_common_geno;
      uint32_t difflist_len;
      reterr = pgr_read_refalt1_difflist_or_genovec_subset_unsafe(sample_include, sample_include_cumulative_popcounts, write_sample_ct, max_simple_difflist_len, vidx, &pgr, genovec, &difflist_common_geno, raregeno, difflist_sample_ids, &difflist_len);
      if (reterr) {
	fprintf(stderr, "\nread error %u, vidx=%u\n", (uint32_t)reterr, vidx);
	goto main_ret_1;
      }
      if (difflist_common_geno == 0xffffffffU) {
        zero_trailing_bits(write_sample_ct * 2, genovec);
	reterr = spgw_append_biallelic_genovec(genovec, &spgw);
      } else if (difflist_len <= max_difflist_len) {
	zero_trailing_bits(2 * difflist_len, raregeno);
	difflist_sample_ids[difflist_len] = write_sample_ct;
	reterr = spgw_append_biallelic_difflist_limited(raregeno, difflist_sample_ids, difflist_common_geno, difflist_len, &spgw);
      } else {
	pgr_difflist_to_genovec_unsafe(raregeno, difflist_sample_ids, difflist_common_geno, write_sample_ct, difflist_len, genovec);
        zero_trailing_bits(write_sample_ct * 2, genovec);
	reterr = spgw_append_biallelic_genovec(genovec, &spgw);
      }
      if (reterr) {
	fprintf(stderr, "\ncompress/write error %u, vidx=%u\n", (uint32_t)reterr, vidx);
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

  spgw_finish(&spgw);
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
  pgr_cleanup(&pgr);
#ifndef NO_MMAP
  pgfi_cleanup(&pgfi);
#endif
  spgw_cleanup(&spgw);
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
  return (uint32_t)reterr;
}
