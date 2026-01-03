// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
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

#include "plink2_perm.h"

#include "include/plink2_bits.h"
#include "include/plink2_string.h"
#include "plink2_cmdline.h"
#include "plink2_random.h"

#ifdef __cplusplus
namespace plink2 {
#endif

void InitPermConfig(PermConfig* perm_config_ptr) {
  perm_config_ptr->within_phenoname = nullptr;
  perm_config_ptr->aperm_min = 6;
  perm_config_ptr->aperm_max = 1000000;
  perm_config_ptr->aperm_alpha = 0.0;
  perm_config_ptr->aperm_beta = 0.0001;
  perm_config_ptr->aperm_init_interval = 1.0;
  perm_config_ptr->aperm_interval_slope = 0.001 * (1 + kSmallEpsilon);
  perm_config_ptr->flags = kfPermColDefault;
}

void CleanupPermConfig(PermConfig* perm_config_ptr) {
  free_cond(perm_config_ptr->within_phenoname);
}

BoolErr PermuteWithinInit(const uintptr_t* sample_include, const PhenoCol* permute_within_phenocol, const uintptr_t* orig_pheno_cc, uint32_t raw_sample_ctl, uint32_t sample_ct, uint32_t* cat_ct_ptr, uint32_t** perm_idx_to_orig_uidx_ptr, uint32_t** cat_cumulative_sizes_ptr, uint32_t** cat_case_cts_ptr) {
  if (!permute_within_phenocol) {
    if (unlikely(bigstack_alloc_u32(sample_ct, perm_idx_to_orig_uidx_ptr) ||
                 bigstack_alloc_u32(2, cat_cumulative_sizes_ptr))) {
      return 1;
    }
    uint32_t* perm_idx_to_orig_uidx = *perm_idx_to_orig_uidx_ptr;
    uint32_t* cat_cumulative_sizes = *cat_cumulative_sizes_ptr;
    cat_cumulative_sizes[0] = 0;
    cat_cumulative_sizes[1] = sample_ct;
    uintptr_t sample_uidx_base = 0;
    uintptr_t cur_bits = sample_include[0];
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
      perm_idx_to_orig_uidx[sample_idx] = sample_uidx;
    }
    if (orig_pheno_cc) {
      if (unlikely(bigstack_alloc_u32(1, cat_case_cts_ptr))) {
        return 1;
      }
      **cat_case_cts_ptr = PopcountWordsIntersect(sample_include, orig_pheno_cc, raw_sample_ctl);
    }
    *cat_ct_ptr = 1;
    return 0;
  }
  uint32_t cat_ct = permute_within_phenocol->nonnull_category_ct + 1;
  if (unlikely(bigstack_alloc_u32(sample_ct, perm_idx_to_orig_uidx_ptr) ||
               bigstack_alloc_u32(cat_ct + 1, cat_cumulative_sizes_ptr))) {
    return 1;
  }
  uint32_t* cat_case_cts = nullptr;
  if (orig_pheno_cc) {
    if (unlikely(bigstack_calloc_u32(cat_ct, cat_case_cts_ptr))) {
      return 1;
    }
    cat_case_cts = *cat_case_cts_ptr;
  }
  uint32_t* perm_idx_to_orig_uidx = *perm_idx_to_orig_uidx_ptr;
  uint32_t* cat_cumulative_sizes = *cat_cumulative_sizes_ptr;
  // 1. compute cat_sizes
  // 2. transform to working write-offsets, fill perm_idx_to_orig_uidx
  // 3. reinterpret final array as cat_cumulative_sizes, set [0] to 0
  // 4. remove empty categories
  uint32_t* cat_sizes = &(cat_cumulative_sizes[1]);
  IdentifyRemainingCatsAndMostCommon(sample_include, permute_within_phenocol, sample_ct, nullptr, cat_sizes);

  uint32_t cat_offset = 0;
  for (uint32_t cat_idx = 0; cat_idx != cat_ct; ++cat_idx) {
    const uint32_t cur_cat_size = cat_sizes[cat_idx];
    cat_sizes[cat_idx] = cat_offset;
    cat_offset += cur_cat_size;
  }
  uint32_t* cat_offsets = cat_sizes;
  const uint32_t* phenocats = permute_within_phenocol->data.cat;
  uintptr_t sample_uidx_base = 0;
  uintptr_t cur_bits = sample_include[0];
  for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
    const uintptr_t sample_uidx = BitIter1(sample_include, &sample_uidx_base, &cur_bits);
    const uint32_t cur_cat_idx = phenocats[sample_uidx];
    perm_idx_to_orig_uidx[cat_offsets[cur_cat_idx]] = sample_uidx;
    cat_offsets[cur_cat_idx] += 1;
    if (cat_case_cts) {
      cat_case_cts[cur_cat_idx] += IsSet(orig_pheno_cc, sample_uidx);
    }
  }

  cat_cumulative_sizes[0] = 0;

  cat_offset = 0;
  uint32_t read_idx = 1;
  for (; read_idx <= cat_ct; ++read_idx) {
    if (cat_offsets[read_idx] == cat_offset) {
      break;
    }
    cat_offset = cat_offsets[read_idx];
  }
  if (read_idx <= cat_ct) {
    uint32_t write_case_ct = cat_case_cts? cat_case_cts[read_idx] : 0;
    uint32_t write_idx = read_idx - 1;
    while (++read_idx <= cat_ct) {
      if (cat_offsets[read_idx] == cat_offset) {
        continue;
      }
      if (cat_case_cts) {
        cat_case_cts[write_idx] = write_case_ct;
        write_case_ct = cat_case_cts[read_idx];
      }
      cat_offsets[write_idx++] = cat_offset;
      cat_offset = cat_offsets[read_idx];
    }
    if (cat_case_cts) {
      cat_case_cts[write_idx] = write_case_ct;
    }
    cat_offsets[write_idx] = cat_offset;
    cat_ct = write_idx;
  }
  *cat_ct_ptr = cat_ct;
  return 0;
}

void PermuteWithinB(const uint32_t* perm_idx_to_orig_uidx, const uint32_t* cat_cumulative_sizes, const uint32_t* cat_case_cts, uint32_t raw_sample_ctl, uint32_t cat_ct, sfmt_t* sfmtp, uintptr_t* permuted_pheno_cc) {
  // not currently a bottleneck, so we don't replicate plink 1.9's
  // more-than-half-cases branches
  ZeroWArr(raw_sample_ctl, permuted_pheno_cc);
  uint32_t cur_cat_offset = 0;
  for (uint32_t cat_idx = 0; cat_idx != cat_ct; ++cat_idx) {
    uint32_t remaining_case_ct = cat_case_cts[cat_idx];
    if (!remaining_case_ct) {
      continue;
    }
    const uint32_t* cat_sample_idx_to_uidx = &(perm_idx_to_orig_uidx[cur_cat_offset]);
    const uint32_t next_cat_offset = cat_cumulative_sizes[cat_idx + 1];
    const uint32_t cur_cat_size = next_cat_offset - cur_cat_offset;
    for (uint32_t cat_sample_idx = 0; cat_sample_idx != cur_cat_size; ++cat_sample_idx) {
      const uint32_t rand_val = RandU32(cur_cat_size - cat_sample_idx, sfmtp);
      if (rand_val < remaining_case_ct) {
        SetBit(cat_sample_idx_to_uidx[cat_sample_idx], permuted_pheno_cc);
        if (--remaining_case_ct == 0) {
          break;
        }
      }
    }
    cur_cat_offset = next_cat_offset;
  }
}

void PermuteWithinD(const uint32_t* perm_idx_to_orig_uidx, const uint32_t* cat_cumulative_sizes, const double* orig_pheno_qt, uint32_t cat_ct, sfmt_t* sfmtp, double* permuted_pheno_qt, uint32_t* perm_buf) {
  const uint32_t sample_ct = cat_cumulative_sizes[cat_ct];
  memcpy(perm_buf, perm_idx_to_orig_uidx, sample_ct * sizeof(int32_t));
  uint32_t cur_cat_offset = 0;
  for (uint32_t cat_idx = 0; cat_idx != cat_ct; ++cat_idx) {
    const uint32_t* cat_sample_idx_to_uidx = &(perm_idx_to_orig_uidx[cur_cat_offset]);
    const uint32_t next_cat_offset = cat_cumulative_sizes[cat_idx + 1];
    const uint32_t cur_cat_size = next_cat_offset - cur_cat_offset;
    uint32_t* cur_perm_buf = &(perm_buf[cur_cat_offset]);
    PermuteU32(cur_cat_size, sfmtp, cur_perm_buf);
    for (uint32_t cat_sample_idx = 0; cat_sample_idx != cur_cat_size; ++cat_sample_idx) {
      permuted_pheno_qt[cat_sample_idx_to_uidx[cat_sample_idx]] = orig_pheno_qt[cur_perm_buf[cat_sample_idx]];
    }
    cur_cat_offset = next_cat_offset;
  }
}

// On entry, *perms_actual_ptr must be the maximum allowed number of
// permutations.
// Returns remaining_variant_ct.
uint32_t UpdateAdaptiveRemainingVariants(const uintptr_t* valid_alleles, const uintptr_t* valid_alleles_cumulative_popcounts_w, const uintptr_t* allele_idx_offsets, const AlleleCode* omitted_alleles, uint32_t raw_variant_ctl, uintptr_t valid_allele_ct, uint32_t perm_idx_start, uint32_t subbatch_size, double adaptive_intercept, double adaptive_slope, uint32_t* valid_allele_emp1_denoms, uintptr_t* remaining_variants, uint32_t* first_adapt_check_pidx_ptr, uint32_t* perms_actual_ptr) {
  const uint32_t perm_idx_start_p1 = perm_idx_start + 1;
  uint32_t first_adapt_check_pidx = *first_adapt_check_pidx_ptr;
  while (first_adapt_check_pidx < subbatch_size) {
    const uint32_t perm_ct = first_adapt_check_pidx + perm_idx_start_p1;
    first_adapt_check_pidx += S_CAST(int32_t, adaptive_intercept + u31tod(perm_ct) * adaptive_slope);
  }
  *first_adapt_check_pidx_ptr = first_adapt_check_pidx - subbatch_size;
  // Update remaining_variants from valid_allele_emp1_denoms.
  uint32_t omitted_allele_idx = 0;
  for (uint32_t variant_widx = 0; variant_widx != raw_variant_ctl; ++variant_widx) {
    uintptr_t prev_variants_word = remaining_variants[variant_widx];
    if (!prev_variants_word) {
      continue;
    }
    uintptr_t next_variants_word = 0;
    if (!allele_idx_offsets) {
      const uint32_t variant_uidx_offset = variant_widx * kBitsPerWord;
      do {
        const uint32_t variant_uidx_lowbits = ctzw(prev_variants_word);
        const uint32_t variant_uidx = variant_uidx_offset + variant_uidx_lowbits;
        // will need to add another case here for non-glm
        if (omitted_alleles) {
          omitted_allele_idx = omitted_alleles[variant_uidx];
        }
        const uint32_t allele_uidx = 2 * variant_uidx + 1 - omitted_allele_idx;
        const uint32_t valid_allele_idx = RawToSubsettedPosW(valid_alleles, valid_alleles_cumulative_popcounts_w, allele_uidx);
        if (!valid_allele_emp1_denoms[valid_allele_idx]) {
          next_variants_word |= k1LU << variant_uidx_lowbits;
        }
        prev_variants_word &= prev_variants_word - 1;
      } while (prev_variants_word);
    } else {
      const uintptr_t* cur_allele_idx_offsets = &(allele_idx_offsets[variant_widx * kBitsPerWord]);
      do {
        const uint32_t variant_uidx_lowbits = ctzw(prev_variants_word);
        const uintptr_t allele_uidx_start = cur_allele_idx_offsets[variant_uidx_lowbits];
        const uintptr_t allele_uidx_stop = cur_allele_idx_offsets[variant_uidx_lowbits + 1];
        const uintptr_t valid_allele_idx_start = RawToSubsettedPosW(valid_alleles, valid_alleles_cumulative_popcounts_w, allele_uidx_start);
        const uintptr_t valid_allele_idx_stop = RawToSubsettedPosW(valid_alleles, valid_alleles_cumulative_popcounts_w, allele_uidx_stop);
        for (uint32_t valid_allele_idx = valid_allele_idx_start; valid_allele_idx != valid_allele_idx_stop; ++valid_allele_idx) {
          if (!valid_allele_emp1_denoms[valid_allele_idx]) {
            next_variants_word |= k1LU << variant_uidx_lowbits;
            break;
          }
        }
        prev_variants_word &= prev_variants_word - 1;
      } while (prev_variants_word);
    }
    remaining_variants[variant_widx] = next_variants_word;
  }
  const uint32_t remaining_variant_ct = PopcountWords(remaining_variants, raw_variant_ctl);
  if (!remaining_variant_ct) {
    *perms_actual_ptr = MaxElementU32(valid_allele_emp1_denoms, valid_allele_ct) - 1;
    return 0;
  }
  const uint32_t perms_total = *perms_actual_ptr;
  if (perm_idx_start + subbatch_size < perms_total) {
    return remaining_variant_ct;
  }
  const uint32_t perms_total_p1 = perms_total + 1;
  for (uintptr_t valid_allele_idx = 0; valid_allele_idx != valid_allele_ct; ++valid_allele_idx) {
    if (!valid_allele_emp1_denoms[valid_allele_idx]) {
      valid_allele_emp1_denoms[valid_allele_idx] = perms_total_p1;
    }
  }
  return 0;
}

PglErr InitPermReportWriter(PermFlags perm_flags, uint32_t perm_adapt, uint32_t provref_col, uintptr_t overflow_buf_size, char* outname, char** outname_end2_ptr, CompressStreamState* css_ptr, char** cswritep_ptr) {
  const uint32_t output_zst = (perm_flags / kfPermZs) & 1;
  const uint32_t chr_col = perm_flags & kfPermColChrom;
  const uint32_t pos_col = perm_flags & kfPermColPos;
  const uint32_t ref_col = perm_flags & kfPermColRef;
  const uint32_t alt1_col = perm_flags & kfPermColAlt1;
  const uint32_t alt_col = perm_flags & kfPermColAlt;
  const uint32_t omitted_col = perm_flags & kfPermColOmitted;
  const uint32_t perm_count = perm_flags & kfPermCounts;

  char* outname_end2 = *outname_end2_ptr;
  *outname_end2++ = '.';
  *outname_end2++ = perm_adapt? 'a' : 'm';
  outname_end2 = strcpya_k(outname_end2, "perm");
  *outname_end2_ptr = outname_end2;
  *outname_end2 = '\0';
  if (output_zst) {
    strcpy_k(outname_end2, ".zst");
  }
  // forced-singlethreaded
  PglErr reterr = InitCstreamAlloc(outname, 0, output_zst, 1, overflow_buf_size, css_ptr, cswritep_ptr);
  if (unlikely(reterr)) {
    return reterr;
  }
  char* cswritep = *cswritep_ptr;
  *cswritep++ = '#';
  if (chr_col) {
    cswritep = strcpya_k(cswritep, "CHROM\t");
  }
  if (pos_col) {
    cswritep = strcpya_k(cswritep, "POS\t");
  }
  cswritep = strcpya_k(cswritep, "ID\t");
  if (ref_col) {
    cswritep = strcpya_k(cswritep, "REF\t");
  }
  if (alt1_col) {
    cswritep = strcpya_k(cswritep, "ALT1\t");
  }
  if (alt_col) {
    cswritep = strcpya_k(cswritep, "ALT\t");
  }
  if (provref_col) {
    cswritep = strcpya_k(cswritep, "PROVISIONAL_REF?\t");
  }
  cswritep = strcpya_k(cswritep, "A1\t");
  if (omitted_col) {
    cswritep = strcpya_k(cswritep, "OMITTED\t");
  }
  if (perm_adapt) {
    if (perm_count) {
      cswritep = strcpya_k(cswritep, "EMP1_CT\tPERM_CT");
    } else {
      cswritep = strcpya_k(cswritep, "EMP1\tPERM_CT");
    }
  } else {
    if (perm_count) {
      cswritep = strcpya_k(cswritep, "EMP1_CT\tEMP2_CT");
    } else {
      cswritep = strcpya_k(cswritep, "EMP1\tEMP2");
    }
  }
  AppendBinaryEoln(&cswritep);
  *cswritep_ptr = cswritep;
  return kPglRetSuccess;
}

BoolErr WritePermReportBody(const uintptr_t* orig_variant_include, const ChrInfo* cip, const uint32_t* variant_bps, const char* const* variant_ids, const uintptr_t* valid_alleles, const uintptr_t* allele_idx_offsets, const char* const* allele_storage, const AlleleCode* omitted_alleles, const uintptr_t* nonref_flags, const double* orig_permstats, const uint32_t* valid_allele_emp1_ctx2s, const uint32_t* valid_allele_emp1_denoms, uint32_t orig_variant_ct, uint32_t all_nonref, uint32_t provref_col, double ln_pfilter, PermFlags perm_flags, uint32_t perms_total, uint32_t lower_stat_is_more_extreme, double* mperm_best_stats, CompressStreamState* css_ptr, char** cswritep_ptr, char* chr_buf) {
  char* cswritep = *cswritep_ptr;
  const uint32_t ref_col = perm_flags & kfPermColRef;
  const uint32_t alt1_col = perm_flags & kfPermColAlt1;
  const uint32_t alt_col = perm_flags & kfPermColAlt;
  const uint32_t omitted_col = perm_flags & kfPermColOmitted;
  const uint32_t perm_count = perm_flags & kfPermCounts;
  uint32_t allele_ct = 2;
  uint32_t omitted_allele_idx = 0;

  const uint32_t perm_adapt = !mperm_best_stats;
  if (!perm_adapt) {
    STD_SORT(perms_total, double_cmp, mperm_best_stats);
  }

  const double pfilter = (ln_pfilter <= 0.0)? exp(ln_pfilter) : 2.0;
  const double emp2_denomx2_recip = 0.5 / u31tod(perms_total + 1);
  uint32_t chr_fo_idx = UINT32_MAX;
  uint32_t chr_end = 0;
  uint32_t chr_buf_blen = 0;
  uint32_t emp1_ct_x2 = 0;
  uintptr_t valid_allele_idx = 0;
  uintptr_t variant_uidx_base = 0;
  uintptr_t cur_bits = orig_variant_include[0];
  for (uint32_t variant_idx = 0; variant_idx != orig_variant_ct; ++variant_idx) {
    const uint32_t variant_uidx = BitIter1(orig_variant_include, &variant_uidx_base, &cur_bits);
    uintptr_t allele_idx_offset_base = variant_uidx * 2;
    if (allele_idx_offsets) {
      allele_idx_offset_base = allele_idx_offsets[variant_uidx];
      allele_ct = allele_idx_offsets[variant_uidx + 1] - allele_idx_offset_base;
    }
    if (omitted_alleles) {
      omitted_allele_idx = omitted_alleles[variant_uidx];
    }
    if (variant_uidx >= chr_end) {
      do {
        ++chr_fo_idx;
        chr_end = cip->chr_fo_vidx_start[chr_fo_idx + 1];
      } while (variant_uidx >= chr_end);
      const uint32_t chr_idx = cip->chr_file_order[chr_fo_idx];
      if (chr_buf) {
        char* chr_name_end = chrtoa(cip, chr_idx, chr_buf);
        *chr_name_end = '\t';
        chr_buf_blen = 1 + S_CAST(uintptr_t, chr_name_end - chr_buf);
      }
    }
    const char* const* cur_alleles = &(allele_storage[allele_idx_offset_base]);
    for (uint32_t aidx = 0; aidx != allele_ct; ++aidx) {
      if (aidx == omitted_allele_idx) {
        continue;
      }
      const uint32_t is_valid = IsSet(valid_alleles, aidx + allele_idx_offset_base);
      double emp1 = 1.5;
      if (is_valid) {
        emp1_ct_x2 = valid_allele_emp1_ctx2s[valid_allele_idx];
        if (perm_adapt) {
          emp1 = u31tod(emp1_ct_x2 + 2) / u31tod(2 * valid_allele_emp1_denoms[valid_allele_idx]);
        } else {
          emp1 = u31tod(emp1_ct_x2 + 2) * emp2_denomx2_recip;
        }
      }
      if (emp1 > pfilter) {
        valid_allele_idx += is_valid;
        continue;
      }
      if (chr_buf) {
        cswritep = memcpya(cswritep, chr_buf, chr_buf_blen);
      }
      if (variant_bps) {
        cswritep = u32toa_x(variant_bps[variant_uidx], '\t', cswritep);
      }
      cswritep = strcpyax(cswritep, variant_ids[variant_uidx], '\t');
      if (ref_col) {
        cswritep = strcpyax(cswritep, cur_alleles[0], '\t');
      }
      if (alt1_col) {
        cswritep = strcpyax(cswritep, cur_alleles[1], '\t');
      }
      if (alt_col) {
        for (uint32_t allele_idx = 1; allele_idx != allele_ct; ++allele_idx) {
          if (unlikely(Cswrite(css_ptr, &cswritep))) {
            *cswritep_ptr = cswritep;
            return 1;
          }
          cswritep = strcpyax(cswritep, cur_alleles[allele_idx], ',');
        }
        cswritep[-1] = '\t';
      }
      if (provref_col) {
        *cswritep++ = (all_nonref || (nonref_flags && IsSet(nonref_flags, variant_uidx)))? 'Y' : 'N';
        *cswritep++ = '\t';
      }
      cswritep = strcpyax(cswritep, cur_alleles[aidx], '\t');
      if (omitted_col) {
        cswritep = strcpyax(cswritep, cur_alleles[omitted_allele_idx], '\t');
      }
      if (perm_adapt) {
        if (!is_valid) {
          cswritep = strcpya_k(cswritep, "NA\tNA");
        } else {
          if (perm_count) {
            cswritep = u32toa(emp1_ct_x2 / 2, cswritep);
            if (emp1_ct_x2 % 2) {
              cswritep = strcpya_k(cswritep, ".5");
            }
          } else {
            cswritep = dtoa_g(emp1, cswritep);
          }
          *cswritep++ = '\t';
          cswritep = u32toa(valid_allele_emp1_denoms[valid_allele_idx] - 1, cswritep);
        }
      } else {
        if (!is_valid) {
          cswritep = strcpya_k(cswritep, "NA\tNA");
        } else {
          const double orig_permstat = orig_permstats[valid_allele_idx];
          const uint32_t emp2_lower_bound = LowerBoundNonemptyD(mperm_best_stats, perms_total, orig_permstat);
          uint32_t emp2_upper_bound = emp2_lower_bound;
          while ((emp2_upper_bound != perms_total) && (mperm_best_stats[emp2_upper_bound] == orig_permstat)) {
            ++emp2_upper_bound;
          }
          uint32_t emp2_ct_x2 = emp2_lower_bound + emp2_upper_bound;
          if (!lower_stat_is_more_extreme) {
            emp2_ct_x2 = perms_total * 2 - emp2_ct_x2;
          }
          if (perm_count) {
            cswritep = u32toa(emp1_ct_x2 / 2, cswritep);
            if (emp1_ct_x2 % 2) {
              cswritep = strcpya_k(cswritep, ".5");
            }
            *cswritep++ = '\t';
            cswritep = u32toa(emp2_ct_x2 / 2, cswritep);
            if (emp2_ct_x2 % 2) {
              cswritep = strcpya_k(cswritep, ".5");
            }
          } else {
            cswritep = dtoa_g(emp1, cswritep);
            *cswritep++ = '\t';
            cswritep = dtoa_g(u31tod(emp2_ct_x2 + 2) * emp2_denomx2_recip, cswritep);
          }
        }
      }
      AppendBinaryEoln(&cswritep);
      if (unlikely(Cswrite(css_ptr, &cswritep))) {
        *cswritep_ptr = cswritep;
        return 1;
      }
      valid_allele_idx += is_valid;
    }
  }
  *cswritep_ptr = cswritep;
  return CswriteCloseNull(css_ptr, cswritep);
}

#ifdef __cplusplus
}  // namespace plink2
#endif
