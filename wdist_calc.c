#include "wdist_common.h"

// The key ideas behind this calculator's design are:
//
// 1. Incremental processing of SNPs.  Each element A_{jk} of a distance or
// relationship matrix is a sum of M terms, one for each SNP, multiplied by a
// missingness correction at the end.  We can walk through the SNPs
// sequentially without keeping much in memory beyond partial sums;
// conveniently, this plays well with the decision made by the PLINK team a few
// years ago to switch to SNP-major binary files.
//
// 2. Multiplexing of markers using bitwise operations.  For instance, there
// are only seven possible ways SNP i can affect the relationship matrix entry
// between individuals j and k:
//    a. j and k are both homozygous for the rare allele
//    b. one is homozygous rare, one is heterozygous
//    c. one is homozygous rare, one is homozygous common
//    d. both are heterozygous
//    e. one is heterozygous, one is homozygous common
//    f. both are homozygous common
//    g. one or both has missing genotype data
// Seven cases can be distinguished by three bits, so one can compose a
// sequence of bitwise operations that maps a pair of padded 2-bit PLINK
// genotypes to seven different 3-bit values according to this breakdown.
// On 64-bit machines, this allows integer operations to act on 20 markers
// simultaneously.  (There's space for a 21st, but we currently choose not to
// use it.)
//
// 3. Lookup tables describing the effect of 5-7 markers at a time on a
// distance or relationship, rather than just one.  For relationship matrices,
// idea #2 allows computation of a single 64-bit integer where bits 0-2
// describe the relationship on marker #0, bits 3-5 describe the relationship
// on marker #1, ..., all the way up to bits 57-59 describing the relationship
// on marker #19.  We then want to perform the update
//    A_{jk} := A_{jk} + f_0(x_0) + f_1(x_1) + ... + f_19(x_19)
// where the x_i's are bit trios, and the f_i's map them to floating point
// terms.  We could do this with 20 table lookups and floating point additions.
// Or, we could structure this as
//    A_{jk} := A_{jk} + f_{0-4}(x_{0-4}) + ... + f_{15-19}(x_{15-19})
// where x_{0-4} denotes the lowest order *15* bits, and f_{0-4} maps them
// directly to f_0(x_0) + f_1(x_1) + f_2(x_2) + f_3(x_3) + f_4(x_4); similarly
// for f_{5-9}, f_{10-14}, and f_{15-19}.  This requires precomputation of four
// lookup tables of size 2^15, total size 1 MB (which fits comfortably in a
// typical L2 cache these days), and licenses the use of four table lookups and
// adds instead of twenty.
//
// 4. Bitslicing algorithms for especially fast calculation of unweighted
// distances and SNP covariances.  Zero-exponent distance matrices and IBS
// matrices are special cases, reducing to Hamming weight computations plus a
// bit of dressing to handle missing markers.  Hamming weight computation on
// x86 CPUs has been studied extensively by others; a good reference as of this
// writing is
//    http://www.dalkescientific.com/writings/diary/archive/2011/11/02/faster_popcount_update.html .
// We use a variant of the Kim Walisch/Cedric Lauradoux bitslicing algorithm
// discussed there (with most 64-bit integer operations replaced by SSE2
// instructions), which runs several times faster than our corresponding
// nonzero exponent distance matrix computation.
//
// We also reduce SNP covariance calculation (used in LD-based pruning) to a
// few Hamming weights.  (This can also be done for covariances between
// individuals, but only when there are no missing markers, so WDIST does not
// include an implementation of that.)
//
// 5. Splitting the distance/relationship matrix into pieces of roughly equal
// size and assigning one thread to each piece.  This is an "embarrassingly
// parallel" problem with no need for careful synchronization.  Cluster
// computing is supported in essentially the same manner.
//
//
//
// In the end, we can't get around the O(MN^2) nature of these calculations,
// but we believe we have beaten down the leading constant by a large enough
// factor to meaningfully help researchers.

void update_rel_ibc(double* rel_ibc, uintptr_t* geno, double* set_allele_freqs, int32_t ibc_type, uint32_t indiv_ct) {
  // first calculate weight array, then loop
  int32_t ii;
  int32_t jj;
  int32_t kk;
  uint32_t uii;
  double twt;
  double* wtptr;
  double mult = 1.0;
  uintptr_t ulii;
  double weights[BITCT * 10];
  double* weights1 = &(weights[64]);
  double* weights2 = &(weights[128]);
  double* weights3 = &(weights[192]);
  double* weights4 = &(weights[256]);
#if __LP64__
  double* weights5 = &(weights[320]);
  double* weights6 = &(weights[384]);
  double* weights7 = &(weights[448]);
  double* weights8 = &(weights[512]);
  double* weights9 = &(weights[576]);
#endif
  double wtarr[BITCT2 * 5];
  double *wptr = weights;
  fill_double_zero(wtarr, BITCT2 * 5);
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if ((set_allele_freqs[ii] != 0.0) && (set_allele_freqs[ii] < (1.0 - EPSILON))) {
      if (ibc_type) {
        if (ibc_type == 2) {
          wtarr[ii * 8] = 2;
          wtarr[ii * 8 + 2] = 2.0 - 1.0 / (2 * set_allele_freqs[ii] * (1.0 - set_allele_freqs[ii]));
          wtarr[ii * 8 + 3] = 2;
        } else {
          twt = 2 * set_allele_freqs[ii];
          if (ibc_type == 1) {
            mult = 1 / (twt * (1.0 - set_allele_freqs[ii]));
          }
          wtarr[ii * 8] = twt * twt * mult;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt) * mult;
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt) * mult;
        }
      } else {
        twt = 1.0 - set_allele_freqs[ii];
        mult = 1 / (set_allele_freqs[ii] * twt);
        wtarr[ii * 8] = 1.0 + set_allele_freqs[ii] * set_allele_freqs[ii] * mult;
        wtarr[ii * 8 + 3] = 1.0 + twt * twt * mult;
      }
    } else {
      if (ibc_type) {
        if (ibc_type == -1) {
          twt = 2 * set_allele_freqs[ii];
          wtarr[ii * 8] = twt * twt;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt);
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt);
        } else if (ibc_type == 1) {
	  wtarr[ii * 8 + 2] = INFINITY;
          if (set_allele_freqs[ii] == 0.0) {
            wtarr[ii * 8] = 0;
            wtarr[ii * 8 + 3] = INFINITY;
          } else {
            wtarr[ii * 8] = INFINITY;
            wtarr[ii * 8 + 3] = 0;
          }
        } else {
          // need to set to 1 instead of 2 for agreement with GCTA
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 2] = -INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      } else {
        if (set_allele_freqs[ii] == 0.0) {
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 3] = INFINITY;
        } else {
          wtarr[ii * 8] = INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      }
    }
  }
  for (kk = 0; kk < (BITCT * 5) / 32; kk++) {
    wtptr = &(wtarr[16 * kk]);
    for (ii = 0; ii < 8; ii++) {
      twt = wtptr[ii + 8];
      for (jj = 0; jj < 8; jj++) {
        *wptr++ = twt + wtptr[jj];
      }
    }
  }
  for (uii = 0; uii < indiv_ct; uii++) {
    ulii = *geno++;
#if __LP64__
    *rel_ibc += weights9[ulii >> 54] + weights8[(ulii >> 48) & 63] + weights7[(ulii >> 42) & 63] + weights6[(ulii >> 36) & 63] + weights5[(ulii >> 30) & 63] + weights4[(ulii >> 24) & 63] + weights3[(ulii >> 18) & 63] + weights2[(ulii >> 12) & 63] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#else
    *rel_ibc += weights4[ulii >> 24] + weights3[(ulii >> 18) & 63] + weights2[(ulii >> 12) & 63] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#endif
    rel_ibc++;
  }
}

void update_rel_f_ibc(float* rel_ibc, uintptr_t* geno, float* set_allele_freqs, int32_t ibc_type, uint32_t indiv_ct) {
  // first calculate weight array, then loop
  int32_t ii;
  int32_t jj;
  int32_t kk;
  uint32_t uii;
  float twt;
  float* wtptr;
  float mult = 1.0;
  uintptr_t ulii;
  float weights[BITCT * 10];
  float* weights1 = &(weights[64]);
  float* weights2 = &(weights[128]);
  float* weights3 = &(weights[192]);
  float* weights4 = &(weights[256]);
#if __LP64__
  float* weights5 = &(weights[320]);
  float* weights6 = &(weights[384]);
  float* weights7 = &(weights[448]);
  float* weights8 = &(weights[512]);
  float* weights9 = &(weights[576]);
#endif
  float wtarr[BITCT2 * 5];
  float *wptr = weights;
  fill_float_zero(wtarr, BITCT2 * 5);
  for (ii = 0; ii < MULTIPLEX_REL / 3; ii += 1) {
    if ((set_allele_freqs[ii] != 0.0) && (set_allele_freqs[ii] < (1.0 - EPSILON))) {
      if (ibc_type) {
        if (ibc_type == 2) {
          wtarr[ii * 8] = 2;
          wtarr[ii * 8 + 2] = 2.0 - 1.0 / (2 * set_allele_freqs[ii] * (1.0 - set_allele_freqs[ii]));
          wtarr[ii * 8 + 3] = 2;
        } else {
          twt = 2 * set_allele_freqs[ii];
          if (ibc_type == 1) {
            mult = 1 / (twt * (1.0 - set_allele_freqs[ii]));
          }
          wtarr[ii * 8] = twt * twt * mult;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt) * mult;
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt) * mult;
        }
      } else {
        twt = 1.0 - set_allele_freqs[ii];
        mult = 1 / (set_allele_freqs[ii] * twt);
        wtarr[ii * 8] = 1.0 + set_allele_freqs[ii] * set_allele_freqs[ii] * mult;
        wtarr[ii * 8 + 3] = 1.0 + twt * twt * mult;
      }
    } else {
      if (ibc_type) {
        if (ibc_type == -1) {
          twt = 2 * set_allele_freqs[ii];
          wtarr[ii * 8] = twt * twt;
          wtarr[ii * 8 + 2] = (1.0 - twt) * (1.0 - twt);
          wtarr[ii * 8 + 3] = (2.0 - twt) * (2.0 - twt);
        } else if (ibc_type == 1) {
	  wtarr[ii * 8 + 2] = INFINITY;
          if (set_allele_freqs[ii] == 0.0) {
            wtarr[ii * 8] = 0;
            wtarr[ii * 8 + 3] = INFINITY;
          } else {
            wtarr[ii * 8] = INFINITY;
            wtarr[ii * 8 + 3] = 0;
          }
        } else {
          // need to set to 1 instead of 2 for agreement with GCTA
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 2] = -INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      } else {
        if (set_allele_freqs[ii] == 0.0) {
          wtarr[ii * 8] = 1;
          wtarr[ii * 8 + 3] = INFINITY;
        } else {
          wtarr[ii * 8] = INFINITY;
          wtarr[ii * 8 + 3] = 1;
        }
      }
    }
  }
  for (kk = 0; kk < (BITCT * 5) / 32; kk++) {
    wtptr = &(wtarr[16 * kk]);
    for (ii = 0; ii < 8; ii++) {
      twt = wtptr[ii + 8];
      for (jj = 0; jj < 8; jj++) {
        *wptr++ = twt + wtptr[jj];
      }
    }
  }
  for (uii = 0; uii < indiv_ct; uii++) {
    ulii = *geno++;
#if __LP64__
    *rel_ibc += weights9[ulii >> 54] + weights8[(ulii >> 48) & 63] + weights7[(ulii >> 42) & 63] + weights6[(ulii >> 36) & 63] + weights5[(ulii >> 30) & 63] + weights4[(ulii >> 24) & 63] + weights3[(ulii >> 18) & 63] + weights2[(ulii >> 12) & 63] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#else
    *rel_ibc += weights4[ulii >> 24] + weights3[(ulii >> 18) & 63] + weights2[(ulii >> 12) & 63] + weights1[(ulii >> 6) & 63] + weights[ulii & 63];
#endif
    rel_ibc++;
  }
}

/*
void collapse_chrom_marker_idxs(Chrom_info* chrom_info_ptr, uintptr_t* marker_exclude, int32_t unfiltered_marker_ct) {
  uint32_t* chrom_fo = chrom_info_ptr->chrom_file_order;
  uint32_t* chrom_fo_midxs = chrom_info_ptr->chrom_file_order_marker_idx;
  uint32_t chrom_ct = chrom_info_ptr->chrom_ct;
  int32_t midx = 0;
  int32_t new_midx = 0;
  int32_t chrom_end_midx;
  uint32_t chrom_fo_idx;
  for (chrom_fo_idx = 0; chrom_fo_idx < chrom_ct; chrom_fo_idx++) {
    chrom_fo_midxs[chrom_fo_idx] = new_midx;
    chrom_info_ptr->chrom_start[chrom_fo[chrom_fo_idx]] = new_midx;
    chrom_end_midx = chrom_fo_midxs[chrom_fo_idx + 1];
    for (; midx < chrom_end_midx; midx++) {
      if (!is_set(marker_exclude, midx)) {
	new_midx++;
      }
    }
    // todo: collapse when chromosome completely eliminated
    chrom_info_ptr->chrom_end[chrom_fo[chrom_fo_idx]] = new_midx;
  }
  chrom_fo_midxs[chrom_fo_idx] = new_midx;
}
*/

  /*
// (using this as a dumping ground for old code that is likely enough to be
// useful later that I don't want to be forced to dredge it from git)

    collapse_arr(marker_alleles, 2, marker_exclude, unfiltered_marker_ct);
    sscanf(output_missing_pheno, "%lg", &missing_phenod);
    // if this becomes much more of a maintenance nightmare, consider exiting
    // function and reloading from .bed from scratch
    if (g_pheno_c) {
      collapse_arr(g_pheno_c, sizeof(char), indiv_exclude, unfiltered_indiv_ct);
    } else if (g_pheno_d) {
      collapse_arr((char*)g_pheno_d, sizeof(double), indiv_exclude, unfiltered_indiv_ct);
    }
    if (sex_info) {
      collapse_arr((char*)sex_info, sizeof(char), indiv_exclude, unfiltered_indiv_ct);
    }
    collapse_bitarr(marker_reverse, marker_exclude, unfiltered_marker_ct);
    if ((calculation_type & CALC_WRITE_SNPLIST) || ((calculation_type & CALC_RECODE) && (recode_modifier & RECODE_LGEN))) {
      collapse_arr(marker_ids, max_marker_id_len, marker_exclude, unfiltered_marker_ct);
    }
    if (calculation_type & (CALC_WRITE_SNPLIST | CALC_GENOME | CALC_LD_PRUNE)) {
      collapse_chrom_marker_idxs(chrom_info_ptr, marker_exclude, unfiltered_marker_ct);
      if (calculation_type & (CALC_GENOME | CALC_LD_PRUNE)) {
	collapse_arr((char*)marker_pos, sizeof(int32_t), marker_exclude, unfiltered_marker_ct);
      }
    }
    collapse_arr(person_ids, max_person_id_len, indiv_exclude, unfiltered_indiv_ct);
    if (fam_col_34) {
      collapse_arr(paternal_ids, max_paternal_id_len, indiv_exclude, unfiltered_indiv_ct);
      collapse_arr(maternal_ids, max_maternal_id_len, indiv_exclude, unfiltered_indiv_ct);
      collapse_bitarr(founder_info, indiv_exclude, unfiltered_indiv_ct);
    }
    if (wt_needed) {
      collapse_arr((char*)g_marker_weights, sizeof(double), marker_exclude, unfiltered_marker_ct);
    }
    collapse_arr((char*)set_allele_freqs, sizeof(double), marker_exclude, unfiltered_marker_ct);
    unfiltered_marker_ct -= marker_exclude_ct;
    marker_exclude_ct = 0;
    fill_ulong_zero(marker_exclude, (unfiltered_marker_ct + (BITCT - 1)) / BITCT);
    unfiltered_indiv_ct -= indiv_exclude_ct;
    unfiltered_indiv_ct4 = (unfiltered_indiv_ct + 3) / 4;
    indiv_exclude_ct = 0;
    fill_ulong_zero(indiv_exclude, (unfiltered_indiv_ct + (BITCT - 1)) / BITCT);
  }
  if (!allow_no_sex) {
    ii = indiv_exclude_ct;
    exclude_ambiguous_sex(unfiltered_indiv_ct, indiv_exclude, &indiv_exclude_ct, sex_info);
    ii = indiv_exclude_ct - ii;
    if (ii) {
      sprintf(logbuf, "%d individual%s with unknown sex removed (prevent with --allow-no-sex).\n", ii, (ii == 1)? "" : "s");
      logprintb();
    }
    g_indiv_ct = unfiltered_indiv_ct - indiv_exclude_ct;
    if (!g_indiv_ct) {
      logprint("Error: No people remaining.\n");
      goto wdist_ret_INVALID_FORMAT;
    }
  */
