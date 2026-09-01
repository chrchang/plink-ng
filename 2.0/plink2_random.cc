// This library is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang, Benjamin Demaille.
//
// This library is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library.  If not, see <http://www.gnu.org/licenses/>.

#include "plink2_random.h"

#include <assert.h>
#include <math.h>

#include "include/plink2_thread.h"
#include "plink2_cmdline.h"
#include "include/pgenlib_write.h"
#include "plink2_decompress.h"

#ifdef __cplusplus
namespace plink2 {
#endif

double RandNormal(sfmt_t* sfmtp, double* secondval_ptr) {
  // Box-Muller.  try changing this to e.g. ziggurat if it's ever a serious
  // bottleneck.
  const double dxx = sqrt(-2 * log(RandUnif(sfmtp)));
  const double dyy = (2 * kPi) * RandUnif(sfmtp);
  *secondval_ptr = dxx * cos(dyy);
  return dxx * sin(dyy);
}

BoolErr InitAllocSfmtpArr(uint32_t thread_ct, uint32_t use_main_sfmt_as_element_zero, sfmt_t* sfmtp, sfmt_t*** sfmtp_arrp) {
  if (unlikely(BIGSTACK_ALLOC_X(sfmt_t*, thread_ct, sfmtp_arrp))) {
    return 1;
  }
  sfmt_t** sfmtp_arr = *sfmtp_arrp;
  if (use_main_sfmt_as_element_zero) {
    sfmtp_arr[0] = sfmtp;
  }
  if (thread_ct > use_main_sfmt_as_element_zero) {
    uint32_t uibuf[4];
    for (uint32_t tidx = use_main_sfmt_as_element_zero; tidx != thread_ct; ++tidx) {
      if (unlikely(BIGSTACK_ALLOC_X(sfmt_t, 1, &(sfmtp_arr[tidx])))) {
        return 1;
      }
      for (uint32_t uii = 0; uii != 4; ++uii) {
        uibuf[uii] = sfmt_genrand_uint32(sfmtp);
      }
      sfmt_init_by_array(sfmtp_arr[tidx], uibuf, 4);
    }
  }
  return 0;
}

typedef struct FillGaussianDArrCtxStruct {
  sfmt_t** sfmtp_arr;
  uintptr_t entry_pair_ct;

  double* dst;
} FillGaussianDArrCtx;

void FillGaussianDArrMain(uintptr_t tidx, uintptr_t thread_ct, FillGaussianDArrCtx* ctx) {
  const uintptr_t entry_pair_ct = ctx->entry_pair_ct;
  sfmt_t* sfmtp = ctx->sfmtp_arr[tidx];
  // 32-bit overflow fix (12 Oct 2019): forgot to cast to uint64_t
  const uintptr_t idx_start = (S_CAST(uint64_t, tidx) * entry_pair_ct) / thread_ct;
  const uintptr_t idx_ct = ((S_CAST(uint64_t, tidx + 1) * entry_pair_ct) / thread_ct) - idx_start;
  double* dst_iter = &(ctx->dst[idx_start * 2]);
  for (uintptr_t ulii = 0; ulii != idx_ct; ++ulii) {
    double dxx;
    *dst_iter++ = RandNormal(sfmtp, &dxx);
    *dst_iter++ = dxx;
  }
}

THREAD_FUNC_DECL FillGaussianDArrThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  FillGaussianDArrCtx* ctx = S_CAST(FillGaussianDArrCtx*, arg->sharedp->context);
  FillGaussianDArrMain(arg->tidx, GetThreadCt(arg->sharedp) + 1, ctx);
  THREAD_RETURN;
}

PglErr FillGaussianDArr(uintptr_t entry_pair_ct, uint32_t thread_ct, sfmt_t* sfmtp, double* darray) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  FillGaussianDArrCtx ctx;
  {
    const uintptr_t max_useful_thread_ct = DivUp(entry_pair_ct, 262144);
    if (thread_ct > max_useful_thread_ct) {
      thread_ct = max_useful_thread_ct;
    }
    if (unlikely(SetThreadCt0(thread_ct - 1, &tg) ||
                 InitAllocSfmtpArr(thread_ct, 1, sfmtp, &ctx.sfmtp_arr))) {
      goto FillGaussianDArr_ret_NOMEM;
    }
    ctx.dst = darray;
    ctx.entry_pair_ct = entry_pair_ct;
    if (thread_ct > 1) {
      SetThreadFuncAndData(FillGaussianDArrThread, &ctx, &tg);
      DeclareLastThreadBlock(&tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto FillGaussianDArr_ret_THREAD_CREATE_FAIL;
      }
    }
    FillGaussianDArrMain(thread_ct - 1, thread_ct, &ctx);
    JoinThreads0(&tg);
  }
  while (0) {
  FillGaussianDArr_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  FillGaussianDArr_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  BigstackReset(bigstack_mark);
  return reterr;
}

typedef struct RandomizeArenaCtxStruct {
  sfmt_t** sfmtp_arr;

  // assumed to be at least 8-byte-aligned
  unsigned char* arena_bottom;
  unsigned char* arena_top;
} RandomizeArenaCtx;

void RandomizeArenaMain(uintptr_t tidx, uintptr_t thread_ct, RandomizeArenaCtx* ctx) {
  unsigned char* arena_bottom = ctx->arena_bottom;
  unsigned char* arena_top = ctx->arena_top;
  const uint64_t arena_int64_ct = S_CAST(uintptr_t, arena_top - arena_bottom) / sizeof(int64_t);
  uint64_t* arena_bottom_alias = R_CAST(uint64_t*, arena_bottom);
  assert(arena_int64_ct >= thread_ct);
  const uint64_t start_idx = RoundDownPow2((tidx * arena_int64_ct) / thread_ct, kInt64PerCacheline);
  uint64_t end_idx = ((tidx + 1) * arena_int64_ct) / thread_ct;
  if (tidx + 1 != thread_ct) {
    end_idx = RoundDownPow2(end_idx, kInt64PerCacheline);
  }
  sfmt_t* sfmtp = ctx->sfmtp_arr[tidx];
  for (uintptr_t ulii = start_idx; ulii != end_idx; ++ulii) {
    arena_bottom_alias[ulii] = sfmt_genrand_uint64(sfmtp);
  }
}

THREAD_FUNC_DECL RandomizeArenaThread(void* raw_arg) {
  ThreadGroupFuncArg* arg = S_CAST(ThreadGroupFuncArg*, raw_arg);
  RandomizeArenaCtx* ctx = S_CAST(RandomizeArenaCtx*, arg->sharedp->context);
  RandomizeArenaMain(arg->tidx, GetThreadCt(arg->sharedp) + 1, ctx);
  THREAD_RETURN;
}

PglErr RandomizeBigstack(uint32_t thread_ct, sfmt_t* sfmtp) {
  unsigned char* bigstack_mark = g_bigstack_base;
  PglErr reterr = kPglRetSuccess;
  ThreadGroup tg;
  PreinitThreads(&tg);
  RandomizeArenaCtx ctx;
  {
    if (thread_ct > 16) {
      thread_ct = 16;
    }
    if (unlikely(SetThreadCt0(thread_ct - 1, &tg) ||
                 InitAllocSfmtpArr(thread_ct, 1, sfmtp, &ctx.sfmtp_arr))) {
      goto RandomizeBigstack_ret_NOMEM;
    }
    ctx.arena_bottom = g_bigstack_base;
    ctx.arena_top = g_bigstack_end;
    if (thread_ct > 1) {
      SetThreadFuncAndData(RandomizeArenaThread, &ctx, &tg);
      DeclareLastThreadBlock(&tg);
      if (unlikely(SpawnThreads(&tg))) {
        goto RandomizeBigstack_ret_THREAD_CREATE_FAIL;
      }
    }
    RandomizeArenaMain(thread_ct - 1, thread_ct, &ctx);
    JoinThreads0(&tg);
    // now ensure the bytes reserved by InitAllocSfmtpArr() are also properly
    // randomized (some of them already are, but there are gaps)
    uint64_t* initial_segment_end = R_CAST(uint64_t*, g_bigstack_base);
    for (uint64_t* initial_segment_iter = R_CAST(uint64_t*, bigstack_mark); initial_segment_iter != initial_segment_end; ++initial_segment_iter) {
      *initial_segment_iter = sfmt_genrand_uint64(sfmtp);
    }
  }
  while (0) {
  RandomizeBigstack_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  RandomizeBigstack_ret_THREAD_CREATE_FAIL:
    reterr = kPglRetThreadCreateFail;
    break;
  }
  CleanupThreads(&tg);
  BigstackReset(bigstack_mark);
  return reterr;
}

void PermuteU32(uint32_t entry_ct, sfmt_t* sfmtp, uint32_t* u32arr) {
  uint32_t* u32arr_iter = u32arr;
  for (uintptr_t remaining_entry_ct = entry_ct; remaining_entry_ct > 1; --remaining_entry_ct) {
    const uint32_t swap_offset = RandU32(remaining_entry_ct, sfmtp);
    // note that swap_offset can be 0, so we can't call swap_u32() if we add
    // __restrict to its definition
    swap_u32(&(u32arr_iter[0]), &(u32arr_iter[swap_offset]));
    ++u32arr_iter;
  }
}

void GeneratePerm1Interleaved(uint32_t tot_bit_ct, uint32_t set_bit_ct, uintptr_t perm_start_idx, uintptr_t perm_end_idx, uintptr_t* perm_buf, sfmt_t* sfmtp) {
  assert(tot_bit_ct > 1);
  const uintptr_t tot_bit_ctl = BitCtToWordCt(tot_bit_ct);
  const uintptr_t perm_ct = perm_end_idx - perm_start_idx;
  const uint32_t tot_quotient = 0x100000000LLU / tot_bit_ct;
  const uint32_t upper_bound = tot_bit_ct * tot_quotient - 1;
  uint32_t totq_preshift;
  uint64_t totq_magic;
  uint32_t totq_postshift;
  uint32_t totq_incr;
  // seeing as how we're gonna divide by the same number a billion times or so,
  // it just might be worth optimizing that division...
  DivisionMagicNums(tot_quotient, &totq_magic, &totq_preshift, &totq_postshift, &totq_incr);
  if (set_bit_ct * 2 < tot_bit_ct) {
    for (uintptr_t widx = 0; widx != tot_bit_ctl; ++widx) {
      ZeroWArr(perm_ct, &(perm_buf[perm_start_idx + (widx * perm_end_idx)]));
    }
    for (uintptr_t perm_idx = perm_start_idx; perm_idx != perm_end_idx; ++perm_idx) {
      uintptr_t* pbptr = &(perm_buf[perm_idx]);
      for (uint32_t num_set = 0; num_set != set_bit_ct; ++num_set) {
        uintptr_t widx;
        uintptr_t shifted_bit;
	do {
          uint32_t urand;
	  do {
	    urand = sfmt_genrand_uint32(sfmtp);
	  } while (urand > upper_bound);
	  // this is identical to lowbits = urand / tot_quotient
	  uintptr_t lowbits = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	  widx = lowbits / kBitsPerWord;
	  lowbits &= (kBitsPerWord - 1);
          shifted_bit = k1LU << lowbits;
	} while (pbptr[widx * perm_end_idx] & shifted_bit);
	pbptr[widx * perm_end_idx] |= shifted_bit;
      }
    }
  } else {
    for (uintptr_t widx = 0; widx != tot_bit_ctl; ++widx) {
      SetAllWArr(perm_ct, &(perm_buf[perm_start_idx + (widx * perm_end_idx)]));
    }
    // "set" has reversed meaning here
    set_bit_ct = tot_bit_ct - set_bit_ct;
    for (uintptr_t perm_idx = perm_start_idx; perm_idx != perm_end_idx; ++perm_idx) {
      uintptr_t* pbptr = &(perm_buf[perm_idx]);
      for (uint32_t num_set = 0; num_set != set_bit_ct; ++num_set) {
        uintptr_t widx;
        uintptr_t shifted_bit;
	do {
          uint32_t urand;
	  do {
	    urand = sfmt_genrand_uint32(sfmtp);
	  } while (urand > upper_bound);
	  uintptr_t lowbits = (totq_magic * ((urand >> totq_preshift) + totq_incr)) >> totq_postshift;
	  widx = lowbits / kBitsPerWord;
	  lowbits &= (kBitsPerWord - 1);
          shifted_bit = k1LU << lowbits;
	} while (!(pbptr[widx * perm_end_idx] & shifted_bit));
	pbptr[widx * perm_end_idx] ^= shifted_bit;
      }
    }
    const uint32_t remaining_bit_ct = tot_bit_ct % kBitsPerWord;
    if (remaining_bit_ct) {
      const uintptr_t final_mask = (~k0LU) >> (kBitsPerWord - remaining_bit_ct);
      uintptr_t* pbptr = &(perm_buf[(tot_bit_ctl - 1) * perm_end_idx + perm_start_idx]);
      for (uintptr_t perm_idx = perm_start_idx; perm_idx != perm_end_idx; ++perm_idx) {
	*pbptr &= final_mask;
	pbptr++;
      }
    }
  }
}

// --simulate / --simulate-qt.  The frequency-table helpers below are ported
// from PLINK 1.9's simulate_init_freqs_{cc,qt}() and
// simulate_cc_get_conditional_probs(), and the generator keeps the same
// random-number call sequence, so the same --seed produces the same dataset in
// both programs.
static const double kTwo63 = 9223372036854775808.0;
CONSTI32(kSimulateMaxAlleleCt, 13);

static void SimulateCcGetConditionalProbs(double prevalence, double g0, double g1, double g2, double het_odds, double hom0_odds, double* f0p, double* f1p, double* f2p) {
  // PLINK 1.07 interprets het_odds and hom0_odds as probability ratios instead
  // of odds ratios.  The two are nearly identical if prevalence is small, but
  // it's still better to avoid that approximation.  This requires solving a
  // cubic equation in X := odds(hom2).
  //
  // prevalence = g0 * (hom0_odds * X) / (1 + hom0_odds * X)
  //            + g1 * (het_odds * X) / (1 + het_odds * X)
  //            + g2 * X / (1 + X)
  //
  // prevalence * (1 + X) * (1 + het_odds * X) * (1 + hom0_odds * X) =
  //   X * (g0 * hom0_odds * (1 + het_odds * X) * (1 + X) +
  //        g1 * het_odds * (1 + hom0_odds * X) * (1 + X) +
  //        g2 * (1 + hom0_odds * X) * (1 + het_odds * X))
  //
  //   X^3 * prevalence * het_odds * hom0_odds
  // + X^2 * prevalence * (het_odds * hom0_odds + het_odds + hom0_odds)
  // + X * prevalence * (1 + het_odds + hom0_odds)
  // + prevalence =
  //   X^3 * het_odds * hom0_odds
  // + X^2 * (g0 * hom0_odds * (1 + het_odds) +
  //          g1 * het_odds * (1 + hom0_odds) +
  //          g2 * (hom0_odds + het_odds))
  // + X * (g0 * hom0_odds + g1 * het_odds + g2)
  //
  // 0 = X^3 * het_odds * hom0_odds * (1 - prevalence)
  //   + X^2 * (g0 * hom0_odds * (1 + het_odds) +
  //            g1 * het_odds * (1 + hom0_odds) +
  //            g2 * (hom0_odds + het_odds) -
  //            prevalence * (het_odds * hom0_odds + het_odds + hom0_odds))
  //   + X * (g0 * hom0_odds + g1 * het_odds + g2 -
  //          prevalence * (1 + het_odds + hom0_odds)
  //   - prevalence
  STD_ARRAY_DECL(double, 3, solutions);
  double coef_recip;
  double cur_f2_odds;
  double cur_f0_odds;
  double cur_f1_odds;
  uint32_t root_ct;
  uint32_t root_idx;
  if ((prevalence == 0) || (prevalence == 1)) {
    *f0p = prevalence;
    *f1p = prevalence;
    *f2p = prevalence;
    return;
  }
  coef_recip = 1.0 / (het_odds * hom0_odds * (1.0 - prevalence));
  // this always has a positive solution since f(0) is negative
  root_ct = CubicRealRoots(coef_recip * (g0 * hom0_odds * (1 + het_odds) + g1 * het_odds + (1 + hom0_odds) + g2 * (hom0_odds + het_odds) - prevalence * (het_odds * hom0_odds + het_odds + hom0_odds)), coef_recip * (g0 * hom0_odds + g1 * het_odds + g2 - prevalence * (1 + het_odds + hom0_odds)), coef_recip * (-prevalence), solutions);
  cur_f2_odds = solutions[0];
  root_idx = 0;
  while ((cur_f2_odds <= 0) && (root_idx + 1 < root_ct)) {
    cur_f2_odds = solutions[++root_idx];
  }
  // odds = p / (1 - p)
  // -> p = odds / (1 + odds)
  cur_f0_odds = cur_f2_odds * hom0_odds;
  cur_f1_odds = cur_f2_odds * het_odds;
  *f0p = cur_f0_odds / (1 + cur_f0_odds);
  *f1p = cur_f1_odds / (1 + cur_f1_odds);
  *f2p = cur_f2_odds / (1 + cur_f2_odds);
}


static void SimulateInitFreqsQt(uint32_t tags_or_haps, double dprime, double qt_var, double qt_dom, double missing_freq, double* freqs, uint64_t* thresholds, double* qt_adj) {
  // Initialize frequency table for current SNP.  Similar to instanceSNP_QT()
  // in PLINK 1.07 simul.cpp.  (The difference is that heterozygote frequencies
  // are not stored in a redundant manner.)
  // thresholds[0]: P(causal 11, marker 11) * 2^63
  // thresholds[1]: P(causal 11, marker 11 or missing) * 2^63
  // thresholds[2]: P(causal 11, marker 11 or missing or 12) * 2^63
  // ...
  // thresholds[14]: [1 - P(causal 22, marker 22)] * 2^63
  // We use 2^63 instead of 2^64 to avoid integer overflow headaches.
  double freq = freqs[0];
  double mfreq = freqs[1];
  double nonmissing_freq = 1 - missing_freq;
  double scaled_miss_freq = missing_freq * kTwo63;
  double scaled_nm_freq = nonmissing_freq * kTwo63;
  double ld = freq * (1 - mfreq);

  double qq = 1 - freq;
  double aa = sqrt(qt_var / (2 * freq * qq * ((1 + qt_dom * (qq - freq)) * (1 + qt_dom * (qq - freq)) + qt_dom * 2 * freq * qq * qt_dom)));
  double dd = qt_dom * aa;
  double dxx = aa * (1 - 2 * freq * (1 + qq * qt_dom));

  // Joint causal variant/marker probs, considering one allele at a time.
  // First digit indicates causal allele, second digit indicates marker.
  double h21 = qq * mfreq;
  double h11;
  double h12;
  double h22;

  // first 2 digits indicate causal variant genotype, last 2 indicate marker.
  // this is DIFFERENT from instanceSNP()'s usage.
  double h_11_11;
  double h_11_12;
  double h_11_22;
  double h_12_11;
  double h_12_12;
  double h_12_22;
  double h_22_11;
  double h_22_12;
  double h_22_22;

  // Constraints:
  // 1. mean = (qq^2) * qt_adj[3] + 2 * freq * qq * qt_adj[2] +
  //           (freq^2) * qt_adj[0]
  //         = 0
  // 2. variance = (qq^2) * (qt_adj[3]^2) + 2 * freq * qq * (qt_adj[2]^2) +
  //               (freq^2) * (qt_adj[0]^2)
  //             = qt_var
  // 3. qt_adj[2] = ((qt_adj[0] + qt_adj[3]) / 2) +
  //                qt_dom * ((qt_adj[0] - qt_adj[3]) / 2)
  // PLINK 1.07 computes the correct values, but incorrectly reverses their
  // positions, and also incorrectly uses sp[s].gAB in some places where .gBB
  // should be used and vice versa.
  qt_adj[0] = dxx + aa;
  qt_adj[2] = dxx + dd;
  qt_adj[3] = dxx - aa;

  if (h21 < ld) {
    ld = h21;
  }
  ld *= dprime;
  h11 = freq * mfreq + ld;
  h12 = freq * (1 - mfreq) - ld;
  h21 -= ld;
  h22 = qq * (1 - mfreq) + ld;
  h_11_11 = h11 * h11;
  h_11_12 = h11 * h12 * 2;
  h_11_22 = h12 * h12;
  h_12_11 = h21 * h11 * 2;
  h_12_12 = (h22 * h11 + h21 * h12) * 2;
  h_12_22 = h22 * h12 * 2;
  h_22_11 = h21 * h21;
  h_22_12 = h22 * h21 * 2;
  h_22_22 = h22 * h22;
  // determination of causal variant missing/nonmissing status is deferred, to
  // simplify phenotype updating
  if (!tags_or_haps) {
    thresholds[0] = S_CAST(uint64_t, (h_11_11 + h_12_11 + h_22_11) * kTwo63);
    thresholds[1] = thresholds[0] + S_CAST(uint64_t, (h_11_12 + h_12_12 + h_22_12) * kTwo63);
  } else {
    thresholds[0] = S_CAST(uint64_t, h_11_11 * scaled_nm_freq);
    thresholds[1] = thresholds[0] + S_CAST(uint64_t, (h_11_11 + h_11_12 + h_11_22) * scaled_miss_freq);
    thresholds[2] = thresholds[1] + S_CAST(uint64_t, h_11_12 * scaled_nm_freq);
    thresholds[3] = thresholds[2] + S_CAST(uint64_t, h_11_22 * scaled_nm_freq);
    thresholds[4] = thresholds[3] + S_CAST(uint64_t, h_12_11 * scaled_nm_freq);
    thresholds[5] = thresholds[4] + S_CAST(uint64_t, (h_12_11 + h_12_12 + h_12_22) * scaled_miss_freq);
    thresholds[6] = thresholds[5] + S_CAST(uint64_t, h_12_12 * scaled_nm_freq);
    thresholds[7] = thresholds[6] + S_CAST(uint64_t, h_12_22 * scaled_nm_freq);
    thresholds[8] = thresholds[7] + S_CAST(uint64_t, h_22_11 * scaled_nm_freq);
    thresholds[9] = thresholds[8] + S_CAST(uint64_t, (h_22_11 + h_22_12 + h_22_22) * scaled_miss_freq);
    thresholds[10] = thresholds[9] + S_CAST(uint64_t, h_22_12 * scaled_nm_freq);
  }
}


static void SimulateInitFreqsCc(uint32_t do_haps, double dprime, double* freqs, double prevalence, double het_odds, double hom0_odds, double missing_freq, uint64_t* ctrl_thresholds, uint64_t* case_thresholds) {
  // Similar to instanceSNP().
  double freq = freqs[0];
  double mfreq = freqs[1];
  double nonmissing_freq = 1 - missing_freq;
  double scaled_missing_freq = missing_freq * kTwo63;
  double scaled_nonmissing_freq = nonmissing_freq * kTwo63;

  // P(causal variant genotype)
  double g0 = freq * freq;
  double g1 = 2 * freq * (1 - freq);
  double g2 = 1 - g0 - g1;
  // P(marker genotype)
  double mg0 = mfreq * mfreq;
  double mg1 = 2 * mfreq * (1 - mfreq);
  double mg2 = 1 - mg0 - mg1;

  double ld = freq * (1 - mfreq);
  double h21 = (1 - freq) * mfreq;
  double h11;
  double h12;
  double h22;

  double f0;
  double f1;
  double f2;
  double mf0;
  double mf1;
  double mf2;
  // first 2 digits indicate causal variant genotype, last 2 indicate marker.
  // this is DIFFERENT from instanceSNP()'s usage.
  double h_11_11;
  double h_11_12;
  double h_11_22;
  double h_12_11;
  double h_12_12;
  double h_12_22;
  double h_22_11;
  double h_22_12;
  double h_22_22;

  double a0;
  double a1;
  double a2;
  double u0;
  double u1;
  double u2;

  double xh_11_11;
  double xh_11_12;
  double xh_11_22;
  double xh_12_11;
  double xh_12_12;
  double xh_12_22;
  double xh_22_11;
  double xh_22_12;
  double xh_22_22;

  double tot_recip;
  double tot_recip_nmsq;
  double tot_recip_1miss;

  if (h21 < ld) {
    ld = h21;
  }
  ld *= dprime;

  // Joint causal variant/marker probs, considering one allele at a time.
  h11 = freq * mfreq + ld;
  h12 = freq * (1 - mfreq) - ld;
  h21 -= ld;
  h22 = (1 - freq) * (1 - mfreq) + ld;


  // Now considering both alleles simultaneously.
  h_11_11 = h11 * h11;
  h_11_12 = h11 * h12 * 2;
  h_11_22 = h12 * h12;
  h_12_11 = h21 * h11 * 2;
  h_12_12 = (h22 * h11 + h21 * h12) * 2;
  h_12_22 = h22 * h12 * 2;
  h_22_11 = h21 * h21;
  h_22_12 = h22 * h21 * 2;
  h_22_22 = h22 * h22;

  // P(case | causal variant genotype)
  SimulateCcGetConditionalProbs(prevalence, g0, g1, g2, het_odds, hom0_odds, &f0, &f1, &f2);
  // P(case | marker genotype)
  mf0 = (f0 * h_11_11 + f1 * h_12_11 + f2 * h_22_11) / mg0;
  mf1 = (f0 * h_11_12 + f1 * h_12_12 + f2 * h_22_12) / mg1;
  mf2 = (f0 * h_11_22 + f1 * h_12_22 + f2 * h_22_22) / mg2;
  // P(marker genotype | affection status)
  a0 = mg0 * mf0;
  a1 = mg1 * mf1;
  a2 = mg2 * mf2;
  tot_recip = 1.0 / (a0 + a1 + a2);
  a0 *= tot_recip;
  a1 *= tot_recip;
  a2 *= tot_recip;
  u0 = mg0 * (1 - mf0);
  u1 = mg1 * (1 - mf1);
  u2 = mg2 * (1 - mf2);
  tot_recip = 1.0 / (u0 + u1 + u2);
  u0 *= tot_recip;
  u1 *= tot_recip;
  u2 *= tot_recip;

  if (!do_haps) {
    case_thresholds[0] = S_CAST(uint64_t, a0 * scaled_nonmissing_freq);
    case_thresholds[1] = case_thresholds[0] + S_CAST(uint64_t, scaled_missing_freq);
    case_thresholds[2] = case_thresholds[1] + S_CAST(uint64_t, a1 * scaled_nonmissing_freq);
    ctrl_thresholds[0] = S_CAST(uint64_t, u0 * scaled_nonmissing_freq);
    ctrl_thresholds[1] = ctrl_thresholds[0] + S_CAST(uint64_t, scaled_missing_freq);
    ctrl_thresholds[2] = ctrl_thresholds[1] + S_CAST(uint64_t, u1 * scaled_nonmissing_freq);
  } else {
    xh_11_11 = h_11_11 * f0;
    xh_11_12 = h_11_12 * f0;
    xh_11_22 = h_11_22 * f0;
    xh_12_11 = h_12_11 * f1;
    xh_12_12 = h_12_12 * f1;
    xh_12_22 = h_12_22 * f1;
    xh_22_11 = h_22_11 * f2;
    xh_22_12 = h_22_12 * f2;
    xh_22_22 = h_22_22 * f2;
    tot_recip = 1.0 / (xh_11_11 + xh_11_12 + xh_11_22 + xh_12_11 + xh_12_12 + xh_12_22 + xh_22_11 + xh_22_12 + xh_22_22);
    tot_recip_nmsq = tot_recip * nonmissing_freq * scaled_nonmissing_freq;
    tot_recip_1miss = tot_recip * nonmissing_freq * scaled_missing_freq;

    case_thresholds[0] = S_CAST(uint64_t, xh_11_11 * tot_recip_nmsq);
    case_thresholds[1] = case_thresholds[0] + S_CAST(uint64_t, (xh_11_11 + xh_11_12 + xh_11_22) * tot_recip_1miss);
    case_thresholds[2] = case_thresholds[1] + S_CAST(uint64_t, xh_11_12 * tot_recip_nmsq);
    case_thresholds[3] = case_thresholds[2] + S_CAST(uint64_t, xh_11_22 * tot_recip_nmsq);
    case_thresholds[4] = case_thresholds[3] + S_CAST(uint64_t, (xh_11_11 + xh_12_11 + xh_22_11) * tot_recip_1miss);
    case_thresholds[5] = case_thresholds[4] + S_CAST(uint64_t, missing_freq * scaled_missing_freq);
    case_thresholds[6] = case_thresholds[5] + S_CAST(uint64_t, (xh_11_12 + xh_12_12 + xh_22_12) * tot_recip_1miss);
    case_thresholds[7] = case_thresholds[6] + S_CAST(uint64_t, (xh_11_22 + xh_12_22 + xh_22_22) * tot_recip_1miss);
    case_thresholds[8] = case_thresholds[7] + S_CAST(uint64_t, xh_12_11 * tot_recip_nmsq);
    case_thresholds[9] = case_thresholds[8] + S_CAST(uint64_t, (xh_12_11 + xh_12_12 + xh_12_22) * tot_recip_1miss);
    case_thresholds[10] = case_thresholds[9] + S_CAST(uint64_t, xh_12_12 * tot_recip_nmsq);
    case_thresholds[11] = case_thresholds[10] + S_CAST(uint64_t, xh_12_22 * tot_recip_nmsq);
    case_thresholds[12] = case_thresholds[11] + S_CAST(uint64_t, xh_22_11 * tot_recip_nmsq);
    case_thresholds[13] = case_thresholds[12] + S_CAST(uint64_t, (xh_22_11 + xh_22_12 + xh_22_22) * tot_recip_1miss);
    case_thresholds[14] = case_thresholds[13] + S_CAST(uint64_t, xh_22_12 * tot_recip_nmsq);

    xh_11_11 = h_11_11 * (1 - f0);
    xh_11_12 = h_11_12 * (1 - f0);
    xh_11_22 = h_11_22 * (1 - f0);
    xh_12_11 = h_12_11 * (1 - f1);
    xh_12_12 = h_12_12 * (1 - f1);
    xh_12_22 = h_12_22 * (1 - f1);
    xh_22_11 = h_22_11 * (1 - f2);
    xh_22_12 = h_22_12 * (1 - f2);
    xh_22_22 = h_22_22 * (1 - f2);
    tot_recip = 1.0 / (xh_11_11 + xh_11_12 + xh_11_22 + xh_12_11 + xh_12_12 + xh_12_22 + xh_22_11 + xh_22_12 + xh_22_22);
    tot_recip_nmsq = tot_recip * nonmissing_freq * scaled_nonmissing_freq;
    tot_recip_1miss = tot_recip * nonmissing_freq * scaled_missing_freq;

    ctrl_thresholds[0] = S_CAST(uint64_t, xh_11_11 * tot_recip_nmsq);
    ctrl_thresholds[1] = ctrl_thresholds[0] + S_CAST(uint64_t, (xh_11_11 + xh_11_12 + xh_11_22) * tot_recip_1miss);
    ctrl_thresholds[2] = ctrl_thresholds[1] + S_CAST(uint64_t, xh_11_12 * tot_recip_nmsq);
    ctrl_thresholds[3] = ctrl_thresholds[2] + S_CAST(uint64_t, xh_11_22 * tot_recip_nmsq);
    ctrl_thresholds[4] = ctrl_thresholds[3] + S_CAST(uint64_t, (xh_11_11 + xh_12_11 + xh_22_11) * tot_recip_1miss);
    ctrl_thresholds[5] = ctrl_thresholds[4] + S_CAST(uint64_t, missing_freq * scaled_missing_freq);
    ctrl_thresholds[6] = ctrl_thresholds[5] + S_CAST(uint64_t, (xh_11_12 + xh_12_12 + xh_22_12) * tot_recip_1miss);
    ctrl_thresholds[7] = ctrl_thresholds[6] + S_CAST(uint64_t, (xh_11_22 + xh_12_22 + xh_22_22) * tot_recip_1miss);
    ctrl_thresholds[8] = ctrl_thresholds[7] + S_CAST(uint64_t, xh_12_11 * tot_recip_nmsq);
    ctrl_thresholds[9] = ctrl_thresholds[8] + S_CAST(uint64_t, (xh_12_11 + xh_12_12 + xh_12_22) * tot_recip_1miss);
    ctrl_thresholds[10] = ctrl_thresholds[9] + S_CAST(uint64_t, xh_12_12 * tot_recip_nmsq);
    ctrl_thresholds[11] = ctrl_thresholds[10] + S_CAST(uint64_t, xh_12_22 * tot_recip_nmsq);
    ctrl_thresholds[12] = ctrl_thresholds[11] + S_CAST(uint64_t, xh_22_11 * tot_recip_nmsq);
    ctrl_thresholds[13] = ctrl_thresholds[12] + S_CAST(uint64_t, (xh_22_11 + xh_22_12 + xh_22_22) * tot_recip_1miss);
    ctrl_thresholds[14] = ctrl_thresholds[13] + S_CAST(uint64_t, xh_22_12 * tot_recip_nmsq);
  }
}


// Number of entries in sorted_u64_arr which ullii is greater than.
static uintptr_t SimulateU64ArrGreaterThan(const uint64_t* sorted_u64_arr, uintptr_t arr_length, uint64_t ullii) {
  intptr_t min_idx = 0;
  intptr_t max_idx = arr_length - 1;
  while (min_idx < max_idx) {
    const uintptr_t mid_idx = (S_CAST(uintptr_t, min_idx) + S_CAST(uintptr_t, max_idx)) / 2;
    if (ullii > sorted_u64_arr[mid_idx]) {
      min_idx = mid_idx + 1;
    } else {
      max_idx = mid_idx - 1;
    }
  }
  if (ullii > sorted_u64_arr[S_CAST(uintptr_t, min_idx)]) {
    return min_idx + 1;
  }
  return min_idx;
}

// Exchanges 00 and 11 in PLINK 1 coding, leaving 01 and 10 untouched.
static void SimulateReverseGenobuf(uint32_t sample_ct, uintptr_t* genobuf) {
  const uint32_t word_ct = NypCtToWordCt(sample_ct);
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t cur_word = genobuf[widx];
    const uintptr_t match_word = (~(cur_word ^ (cur_word >> 1))) & kMask5555;
    genobuf[widx] = cur_word ^ (match_word | (match_word << 1));
  }
}

// PLINK 1 coding (0 = hom A1, 1 = missing, 2 = het, 3 = hom A2) to plink2's
// (0 = hom REF, 1 = het, 2 = hom ALT, 3 = missing), with A1 = ALT.
static void SimulatePlink1ToGenovec(uint32_t sample_ct, const uintptr_t* src, uintptr_t* dst) {
  const uint32_t word_ct = NypCtToWordCt(sample_ct);
  for (uint32_t widx = 0; widx != word_ct; ++widx) {
    const uintptr_t cur_word = src[widx];
    const uintptr_t lo = cur_word & kMask5555;
    const uintptr_t hi = (cur_word >> 1) & kMask5555;
    dst[widx] = (((~hi) & kMask5555) << 1) | (lo ^ hi);
  }
}

PglErr SimulateDataset(const char* simulate_fname, const char* name_prefix, SimulateFlags flags, uint32_t case_ct, uint32_t ctrl_ct, uint32_t qt_sample_ct, double prevalence, double missing_freq, sfmt_t* sfmtp, char* outname, char* outname_end) {
  unsigned char* bigstack_mark = g_bigstack_base;
  FILE* outfile = nullptr;
  FILE* simfreq_file = nullptr;
  uintptr_t line_idx = 0;
  TextStream txs;
  PreinitTextStream(&txs);
  STPgenWriter spgw;
  PreinitSpgw(&spgw);
  PglErr reterr = kPglRetSuccess;
  {
    const uint32_t is_qt = (flags / kfSimulateQt) & 1;
    const uint32_t do_tags = (flags / kfSimulateTags) & 1;
    const uint32_t do_haps = (flags / kfSimulateHaps) & 1;
    const uint32_t tags_or_haps = do_tags | do_haps;
    const uint32_t randomize_alleles = (flags & (kfSimulateAcgt | kfSimulate1234 | kfSimulate12))? 1 : 0;
    const uint32_t simulate_12 = (flags / kfSimulate12) & 1;
    const uint32_t missing_thresh = S_CAST(uint32_t, missing_freq * 4294967296.0);
    const uint32_t sample_ct = is_qt? qt_sample_ct : (case_ct + ctrl_ct);
    if (unlikely(!sample_ct)) {
      logerrputs("Error: --simulate requires at least one sample.\n");
      goto SimulateDataset_ret_INCONSISTENT_INPUT;
    }
    const uint32_t sample_ctl2 = NypCtToWordCt(sample_ct);
    char alleles[kSimulateMaxAlleleCt];
    if (randomize_alleles) {
      if (flags & kfSimulateAcgt) {
        memcpy(alleles, "ACAGATCGCTGTA", 13);
      } else if (flags & kfSimulate1234) {
        memcpy(alleles, "1213142324341", 13);
      } else {
        memcpy(alleles, "121", 3);
      }
    } else {
      if (is_qt) {
        memcpy(alleles, "HLAB", 4);
      } else {
        memcpy(alleles, "DdAB", 4);
      }
    }
    double* qt_vals = nullptr;
    uintptr_t* writebuf;
    uintptr_t* writebuf2 = nullptr;
    uintptr_t* genovec;
    if (unlikely(bigstack_alloc_w(sample_ctl2, &writebuf) ||
                 bigstack_alloc_w(sample_ctl2, &genovec))) {
      goto SimulateDataset_ret_NOMEM;
    }
    if (do_haps) {
      if (unlikely(bigstack_alloc_w(sample_ctl2, &writebuf2))) {
        goto SimulateDataset_ret_NOMEM;
      }
    }
    if (is_qt) {
      if (unlikely(bigstack_calloc_d(sample_ct, &qt_vals))) {
        goto SimulateDataset_ret_NOMEM;
      }
    }

    // sfmt_genrand_uint64() must not be mixed with sfmt_genrand_uint32() on
    // the same generator, so a second one is seeded from the main one.
    sfmt_t* sfmt64p;
    if (unlikely(BIGSTACK_ALLOC_X(sfmt_t, 1, &sfmt64p))) {
      goto SimulateDataset_ret_NOMEM;
    }
    {
      uint32_t init_arr[4];
      for (uint32_t uii = 0; uii != 4; ++uii) {
        init_arr[uii] = sfmt_genrand_uint32(sfmtp);
      }
      sfmt_init_by_array(sfmt64p, init_arr, 4);
    }

    // First pass: total variant count, for the .pgen header.
    uint64_t variant_ct_u64 = 0;
    reterr = SizeAndInitTextStream(simulate_fname, bigstack_left() / 4, 1, &txs);
    if (unlikely(reterr)) {
      goto SimulateDataset_ret_TSTREAM_FAIL;
    }
    while (1) {
      ++line_idx;
      const char* line_start = TextGet(&txs);
      if (!line_start) {
        break;
      }
      if (IsEolnKns(*line_start)) {
        continue;
      }
      uint32_t cur_ct;
      if (unlikely(ScanUintCappedx(line_start, 0x7ffffffd, &cur_ct))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid SNP count on line %" PRIuPTR " of --simulate%s input file.\n", line_idx, is_qt? "-qt" : "");
        goto SimulateDataset_ret_MALFORMED_INPUT_WW;
      }
      variant_ct_u64 += cur_ct;
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto SimulateDataset_ret_TSTREAM_FAIL;
    }
    if (unlikely(!variant_ct_u64)) {
      snprintf(g_logbuf, kLogbufSize, "Error: --simulate%s input file specifies zero SNPs.\n", is_qt? "-qt" : "");
      goto SimulateDataset_ret_MALFORMED_INPUT_WW;
    }
    if (unlikely(variant_ct_u64 > (do_haps? 0x3ffffffeLLU : 0x7ffffffdLLU))) {
      snprintf(g_logbuf, kLogbufSize, "Error: --simulate%s input file specifies too many SNPs.\n", is_qt? "-qt" : "");
      goto SimulateDataset_ret_MALFORMED_INPUT_WW;
    }
    const uint32_t variant_ct = variant_ct_u64 * (1 + do_haps);
    if (unlikely(CleanupTextStream2(simulate_fname, &txs, &reterr))) {
      goto SimulateDataset_ret_1;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".pvar");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto SimulateDataset_ret_OPEN_FAIL;
    }
    fputs("#CHROM\tPOS\tID\tREF\tALT" EOLN_STR, outfile);
    snprintf(outname_end, kMaxOutfnameExtBlen, ".simfreq");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &simfreq_file))) {
      goto SimulateDataset_ret_OPEN_FAIL;
    }
    snprintf(outname_end, kMaxOutfnameExtBlen, ".pgen");
    uintptr_t spgw_alloc_cacheline_ct;
    uint32_t max_vrec_len;
    reterr = SpgwInitPhase1(outname, nullptr, nullptr, variant_ct, sample_ct, 0, kPgenWriteBackwardSeek, kfPgenGlobal0, 2, &spgw, &spgw_alloc_cacheline_ct, &max_vrec_len);
    if (unlikely(reterr)) {
      goto SimulateDataset_ret_1;
    }
    unsigned char* spgw_alloc;
    if (unlikely(bigstack_alloc_uc(spgw_alloc_cacheline_ct * kCacheline, &spgw_alloc))) {
      goto SimulateDataset_ret_NOMEM;
    }
    SpgwInitPhase2(max_vrec_len, &spgw, spgw_alloc);
    *outname_end = '\0';
    logprintfww5("Writing --simulate%s dataset to %s.pgen + %s.pvar + %s.psam ... ", is_qt? "-qt" : "", outname, outname, outname);
    fflush(stdout);

    reterr = SizeAndInitTextStream(simulate_fname, bigstack_left() / 4, 1, &txs);
    if (unlikely(reterr)) {
      goto SimulateDataset_ret_TSTREAM_FAIL;
    }
    char* textbuf = g_textbuf;
    char snp_label[kMaxIdSlen + 16];
    uint64_t thresholds[15];
    uint64_t case_thresholds[15];
    double freqs[6];
    double qt_adj[4];
    char cur_alleles[4];
    double qt_totvar = 0.0;
    uint32_t bp = 1;
    uint32_t zero_odds_ratio_warning_given = 0;
    line_idx = 0;
    while (1) {
      ++line_idx;
      char* line_start = TextGet(&txs);
      if (!line_start) {
        break;
      }
      if (IsEolnKns(*line_start)) {
        continue;
      }
      char* snp_label_ptr = NextToken(line_start);
      char* freq_lb_ptr = NextToken(snp_label_ptr);
      char* marker_freq_lb_ptr = nullptr;
      char* marker_ld_ptr = nullptr;
      char* penult_ptr;
      if (tags_or_haps) {
        marker_freq_lb_ptr = NextTokenMult(freq_lb_ptr, 2);
        marker_ld_ptr = NextTokenMult(marker_freq_lb_ptr, 2);
        penult_ptr = NextToken(marker_ld_ptr);
      } else {
        penult_ptr = NextTokenMult(freq_lb_ptr, 2);
      }
      char* last_ptr = penult_ptr? NextToken(penult_ptr) : nullptr;
      if (unlikely((!last_ptr) || IsEolnKns(*last_ptr))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --simulate%s file has fewer tokens than expected.\n", line_idx, is_qt? "-qt" : "");
        goto SimulateDataset_ret_MALFORMED_INPUT_WW;
      }
      {
        const char* extra_ptr = NextToken(last_ptr);
        if (unlikely(extra_ptr && (!IsEolnKns(*extra_ptr)))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Line %" PRIuPTR " of --simulate%s file has more tokens than expected.\n", line_idx, is_qt? "-qt" : "");
          goto SimulateDataset_ret_MALFORMED_INPUT_WW;
        }
      }
      uint32_t cur_variant_ct;
      ScanUintCappedx(line_start, 0x7ffffffd, &cur_variant_ct);
      if (!cur_variant_ct) {
        continue;
      }
      uint32_t snp_label_slen = CurTokenEnd(snp_label_ptr) - snp_label_ptr;
      if (unlikely(snp_label_slen > kMaxIdSlen - 16)) {
        snprintf(g_logbuf, kLogbufSize, "Error: Overlong SNP label on line %" PRIuPTR " of --simulate%s file.\n", line_idx, is_qt? "-qt" : "");
        goto SimulateDataset_ret_MALFORMED_INPUT_WW;
      }
      memcpy(snp_label, snp_label_ptr, snp_label_slen);
      snp_label[snp_label_slen++] = '_';
      double freq_lb;
      double freq_delta;
      if (unlikely((!ScantokDouble(freq_lb_ptr, &freq_lb)) || (!ScantokDouble(NextToken(freq_lb_ptr), &freq_delta)) || (freq_lb < 0) || (freq_delta < freq_lb) || (freq_delta > 1))) {
        snprintf(g_logbuf, kLogbufSize, "Error: Invalid allele frequency bound on line %" PRIuPTR " of --simulate%s file.\n", line_idx, is_qt? "-qt" : "");
        goto SimulateDataset_ret_MALFORMED_INPUT_WW;
      }
      freq_delta -= freq_lb;
      double marker_freq_lb = 0.0;
      double marker_freq_ub = 0.0;
      double dprime = 1.0;
      if (tags_or_haps) {
        if (unlikely((!ScantokDouble(marker_freq_lb_ptr, &marker_freq_lb)) || (!ScantokDouble(NextToken(marker_freq_lb_ptr), &marker_freq_ub)) || (marker_freq_lb < 0) || (marker_freq_ub < marker_freq_lb) || (marker_freq_ub > 1))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid marker allele frequency bound on line %" PRIuPTR " of --simulate%s file.\n", line_idx, is_qt? "-qt" : "");
          goto SimulateDataset_ret_MALFORMED_INPUT_WW;
        }
        if (unlikely((!ScantokDouble(marker_ld_ptr, &dprime)) || (dprime < 0) || (dprime > 1))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid d-prime on line %" PRIuPTR " of --simulate%s input file.\n", line_idx, is_qt? "-qt" : "");
          goto SimulateDataset_ret_MALFORMED_INPUT_WW;
        }
      }
      double qt_var = 0.0;
      double qt_dom = 0.0;
      double het_odds = 0.0;
      double hom0_odds = 0.0;
      if (is_qt) {
        if (unlikely((!ScantokDouble(penult_ptr, &qt_var)) || (qt_var < 0) || (qt_var > 1))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid variance value on line %" PRIuPTR " of --simulate-qt file.\n", line_idx);
          goto SimulateDataset_ret_MALFORMED_INPUT_WW;
        }
        if (unlikely((qt_var > 0) && (((freq_delta == 0) && ((freq_lb == 0) || (freq_lb == 1))) || (tags_or_haps && (marker_freq_lb == marker_freq_ub) && ((marker_freq_lb == 0) || (marker_freq_lb == 1)))))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Nonzero variance with fixed 0/1 allele frequency on line %" PRIuPTR " of --simulate-qt file.\n", line_idx);
          goto SimulateDataset_ret_MALFORMED_INPUT_WW;
        }
        qt_totvar += S_CAST(double, S_CAST(intptr_t, cur_variant_ct)) * qt_var;
        if (unlikely(qt_totvar > 1 + kSmallEpsilon)) {
          logputs("\n");
          logerrputs("Error: --simulate-qt input file specifies QTL variance greater than 1.\n");
          goto SimulateDataset_ret_MALFORMED_INPUT;
        }
        if (unlikely(!ScantokDouble(last_ptr, &qt_dom))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid dominance deviation value on line %" PRIuPTR " of --simulate-qt file.\n", line_idx);
          goto SimulateDataset_ret_MALFORMED_INPUT_WW;
        }
      } else {
        if (unlikely((!ScantokDouble(penult_ptr, &het_odds)) || (het_odds < 0))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid heterozygote disease odds ratio on line %" PRIuPTR " of --simulate file.\n", line_idx);
          goto SimulateDataset_ret_MALFORMED_INPUT_WW;
        }
        const uint32_t last_slen = CurTokenEnd(last_ptr) - last_ptr;
        if ((last_slen == 4) && MatchUpperCounted(last_ptr, "MULT", 4)) {
          hom0_odds = het_odds * het_odds;
        } else if (unlikely((!ScantokDouble(last_ptr, &hom0_odds)) || (hom0_odds < 0))) {
          snprintf(g_logbuf, kLogbufSize, "Error: Invalid homozygote disease odds ratio on line %" PRIuPTR " of --simulate file.\n", line_idx);
          goto SimulateDataset_ret_MALFORMED_INPUT_WW;
        }
        if ((!zero_odds_ratio_warning_given) && ((het_odds == 0) || (hom0_odds == 0))) {
          logputs("\n");
          logerrputs("Warning: Zero odds ratio present in --simulate input file.  Did you mean\n--simulate-qt instead?\n");
          zero_odds_ratio_warning_given = 1;
        }
      }
      for (uint32_t cur_variant_idx = 0; cur_variant_idx != cur_variant_ct; ++cur_variant_idx) {
        freqs[0] = freq_lb + RandUnif(sfmtp) * freq_delta;
        if (tags_or_haps) {
          // Keep marker freq / causal freq and their complements' ratio inside
          // [dprime, 1 / dprime], unless that puts the marker frequency out of
          // range.
          double dxx;
          double dyy;
          if (dprime > 0) {
            dxx = 1 - (1 - freqs[0]) / dprime;
            double dzz = freqs[0] * dprime;
            if (dxx < dzz) {
              dxx = dzz;
            }
            if (dxx < marker_freq_lb) {
              dxx = marker_freq_lb;
            }
            dyy = 1 - (1 - freqs[0]) * dprime;
            dzz = freqs[0] / dprime;
            if (dyy > dzz) {
              dyy = dzz;
            }
            if (dyy > marker_freq_ub) {
              dyy = marker_freq_ub;
            }
            if (dyy < dxx) {
              if (dyy < marker_freq_lb) {
                dyy = marker_freq_lb;
              }
              dxx = dyy;
            }
          } else {
            dxx = marker_freq_lb;
            dyy = marker_freq_ub;
          }
          freqs[1] = dxx + RandUnif(sfmtp) * (dyy - dxx);
        } else {
          freqs[1] = freqs[0];
        }
        if (is_qt) {
          SimulateInitFreqsQt(do_haps, dprime, qt_var, qt_dom, missing_freq, freqs, thresholds, qt_adj);
        } else {
          SimulateInitFreqsCc(do_haps, dprime, freqs, prevalence, het_odds, hom0_odds, missing_freq, thresholds, case_thresholds);
        }
        // .simfreq: the realized parameters for this variant.
        {
          char* write_iter = textbuf;
          *write_iter++ = '1';
          *write_iter++ = ' ';
          if (cur_variant_ct > 1) {
            write_iter = memcpya(write_iter, snp_label, snp_label_slen);
            write_iter = u32toa(cur_variant_idx, write_iter);
          } else {
            write_iter = memcpya(write_iter, snp_label, snp_label_slen - 1);
          }
          *write_iter++ = '\t';
          write_iter = dtoa_g(freqs[0], write_iter); *write_iter++ = ' ';
          write_iter = dtoa_g(freqs[0], write_iter); *write_iter++ = '\t';
          if (tags_or_haps) {
            write_iter = dtoa_g(freqs[1], write_iter); *write_iter++ = ' ';
            write_iter = dtoa_g(freqs[1], write_iter); *write_iter++ = '\t';
            write_iter = dtoa_g(dprime, write_iter); *write_iter++ = '\t';
          }
          if (is_qt) {
            write_iter = dtoa_g(qt_var, write_iter); *write_iter++ = '\t';
            write_iter = dtoa_g(qt_dom, write_iter); *write_iter++ = '\n';
          } else {
            write_iter = dtoa_g(het_odds, write_iter); *write_iter++ = '\t';
            write_iter = dtoa_g(hom0_odds, write_iter); *write_iter++ = '\n';
          }
          if (unlikely(fwrite_checked(textbuf, write_iter - textbuf, simfreq_file))) {
            goto SimulateDataset_ret_WRITE_FAIL;
          }
        }
        if (randomize_alleles) {
          uint32_t uii;
          uint32_t ujj;
          if (!simulate_12) {
            do {
              uii = sfmt_genrand_uint32(sfmtp);
            } while (uii >= 4294967184U);
            uii = uii % 144U;
            ujj = uii / 12;
            uii -= ujj * 12;
          } else {
            uii = sfmt_genrand_uint32(sfmtp) & 3;
            ujj = uii >> 1;
            uii &= 1;
          }
          memcpy(cur_alleles, &(alleles[uii]), 2);
          memcpy(&(cur_alleles[2]), &(alleles[ujj]), 2);
        } else {
          memcpy(cur_alleles, alleles, 4);
        }
        uintptr_t ulii = 0;
        uintptr_t uljj = 0;
        uint32_t ukk = 0;
        uintptr_t* wbptr = writebuf;
        uintptr_t* wbptr2 = writebuf2;
        if (!do_haps) {
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            const uint64_t ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
            uintptr_t ulkk;
            if (is_qt) {
              if (ullii > thresholds[1]) {
                ulkk = 3;
              } else if (ullii > thresholds[0]) {
                ulkk = 2;
              } else {
                ulkk = 0;
              }
              qt_vals[sample_idx] += qt_adj[ulkk];
              if (sfmt_genrand_uint32(sfmtp) < missing_thresh) {
                ulkk = 1;
              }
            } else {
              const uint64_t* cur_thresholds = (sample_idx < case_ct)? case_thresholds : thresholds;
              if (ullii > cur_thresholds[1]) {
                if (ullii > cur_thresholds[2]) {
                  ulkk = 3;
                } else {
                  ulkk = 2;
                }
              } else if (ullii > cur_thresholds[0]) {
                ulkk = 1;
              } else {
                ulkk = 0;
              }
            }
            ulii |= ulkk << ukk;
            ukk += 2;
            if (ukk == kBitsPerWord) {
              *wbptr++ = ulii;
              ulii = 0;
              ukk = 0;
            }
          }
          if (ukk) {
            *wbptr = ulii;
          }
        } else {
          for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
            const uint64_t ullii = sfmt_genrand_uint64(sfmt64p) >> 1;
            uintptr_t ulkk;
            uintptr_t ulmm;
            if (is_qt) {
              ulkk = SimulateU64ArrGreaterThan(thresholds, 11, ullii);
              ulmm = ulkk & 3;
              ulkk /= 4;
              ulkk += (ulkk + 1) >> 1;
              qt_vals[sample_idx] += qt_adj[ulkk];
              if (sfmt_genrand_uint32(sfmtp) < missing_thresh) {
                ulkk = 1;
              }
            } else {
              const uint64_t* cur_thresholds = (sample_idx < case_ct)? case_thresholds : thresholds;
              ulkk = SimulateU64ArrGreaterThan(cur_thresholds, 15, ullii);
              ulmm = ulkk & 3;
              ulkk /= 4;
            }
            ulii |= ulkk << ukk;
            uljj |= ulmm << ukk;
            ukk += 2;
            if (ukk == kBitsPerWord) {
              *wbptr++ = ulii;
              *wbptr2++ = uljj;
              ulii = 0;
              uljj = 0;
              ukk = 0;
            }
          }
          if (ukk) {
            *wbptr = ulii;
            *wbptr2 = uljj;
          }
        }
        // PLINK 1.x makes the more common allele A2.
        if (PopcountWords(writebuf, sample_ctl2) < sample_ct) {
          SimulateReverseGenobuf(sample_ct, writebuf);
          const char cc = cur_alleles[0];
          cur_alleles[0] = cur_alleles[1];
          cur_alleles[1] = cc;
        }
        {
          char* write_iter = textbuf;
          write_iter = strcpya_k(write_iter, "1\t");
          write_iter = u32toa_x(bp++, '\t', write_iter);
          if (cur_variant_ct > 1) {
            write_iter = memcpya(write_iter, snp_label, snp_label_slen);
            write_iter = u32toa(cur_variant_idx, write_iter);
          } else {
            write_iter = memcpya(write_iter, snp_label, snp_label_slen - 1);
          }
          if (do_tags) {
            write_iter = strcpya_k(write_iter, "_M");
          }
          *write_iter++ = '\t';
          *write_iter++ = cur_alleles[1];
          *write_iter++ = '\t';
          *write_iter++ = cur_alleles[0];
          AppendBinaryEoln(&write_iter);
          if (unlikely(fwrite_checked(textbuf, write_iter - textbuf, outfile))) {
            goto SimulateDataset_ret_WRITE_FAIL;
          }
        }
        SimulatePlink1ToGenovec(sample_ct, writebuf, genovec);
        ZeroTrailingNyps(sample_ct, genovec);
        if (unlikely(SpgwAppendBiallelicGenovec(genovec, &spgw))) {
          goto SimulateDataset_ret_WRITE_FAIL;
        }
        if (do_haps) {
          if (PopcountWords(writebuf2, sample_ctl2) < sample_ct) {
            SimulateReverseGenobuf(sample_ct, writebuf2);
            const char cc = cur_alleles[2];
            cur_alleles[2] = cur_alleles[3];
            cur_alleles[3] = cc;
          }
          char* write_iter = textbuf;
          write_iter = strcpya_k(write_iter, "1\t");
          write_iter = u32toa_x(bp++, '\t', write_iter);
          write_iter = memcpya(write_iter, snp_label, snp_label_slen);
          write_iter = u32toa(cur_variant_idx, write_iter);
          write_iter = strcpya_k(write_iter, "_M\t");
          *write_iter++ = cur_alleles[3];
          *write_iter++ = '\t';
          *write_iter++ = cur_alleles[2];
          AppendBinaryEoln(&write_iter);
          if (unlikely(fwrite_checked(textbuf, write_iter - textbuf, outfile))) {
            goto SimulateDataset_ret_WRITE_FAIL;
          }
          SimulatePlink1ToGenovec(sample_ct, writebuf2, genovec);
          ZeroTrailingNyps(sample_ct, genovec);
          if (unlikely(SpgwAppendBiallelicGenovec(genovec, &spgw))) {
            goto SimulateDataset_ret_WRITE_FAIL;
          }
        }
      }
    }
    if (unlikely(TextStreamErrcode2(&txs, &reterr))) {
      goto SimulateDataset_ret_TSTREAM_FAIL;
    }
    reterr = SpgwFinish(&spgw);
    if (unlikely(reterr)) {
      goto SimulateDataset_ret_1;
    }
    if (unlikely(fclose_null(&outfile) || fclose_null(&simfreq_file))) {
      goto SimulateDataset_ret_WRITE_FAIL;
    }

    snprintf(outname_end, kMaxOutfnameExtBlen, ".psam");
    if (unlikely(fopen_checked(outname, FOPEN_WB, &outfile))) {
      goto SimulateDataset_ret_OPEN_FAIL;
    }
    fputs(is_qt? ("#FID\tIID\tSEX\tPHENO1" EOLN_STR) : ("#FID\tIID\tSEX\tPHENO1" EOLN_STR), outfile);
    const uint32_t name_prefix_slen = name_prefix? strlen(name_prefix) : 0;
    double qt_resid_sd = 0.0;
    if (is_qt && (qt_totvar < 1 - kSmallEpsilon)) {
      qt_resid_sd = sqrt(1 - qt_totvar);
    }
    double secondval = 0.0;
    for (uint32_t sample_idx = 0; sample_idx != sample_ct; ++sample_idx) {
      char* write_iter = textbuf;
      for (uint32_t uii = 0; uii != 2; ++uii) {
        if (name_prefix_slen) {
          write_iter = memcpyax(write_iter, name_prefix, name_prefix_slen, '-');
        }
        write_iter = strcpya_k(write_iter, "per");
        write_iter = u32toa_x(sample_idx, '\t', write_iter);
      }
      write_iter = strcpya_k(write_iter, "2\t");
      if (is_qt) {
        double cur_pheno;
        if (sample_idx & 1) {
          cur_pheno = qt_vals[sample_idx] + qt_resid_sd * secondval;
        } else {
          cur_pheno = qt_vals[sample_idx] + qt_resid_sd * RandNormal(sfmtp, &secondval);
        }
        write_iter = dtoa_g(cur_pheno, write_iter);
      } else {
        *write_iter++ = (sample_idx < case_ct)? '2' : '1';
      }
      AppendBinaryEoln(&write_iter);
      if (unlikely(fwrite_checked(textbuf, write_iter - textbuf, outfile))) {
        goto SimulateDataset_ret_WRITE_FAIL;
      }
    }
    if (unlikely(fclose_null(&outfile))) {
      goto SimulateDataset_ret_WRITE_FAIL;
    }
    *outname_end = '\0';
    logputs("done.\n");
    logprintfww("Realized simulation parameters saved to %s.simfreq .\n", outname);
  }
  while (0) {
  SimulateDataset_ret_NOMEM:
    reterr = kPglRetNomem;
    break;
  SimulateDataset_ret_OPEN_FAIL:
    reterr = kPglRetOpenFail;
    break;
  SimulateDataset_ret_TSTREAM_FAIL:
    TextStreamErrPrint(simulate_fname, &txs);
    break;
  SimulateDataset_ret_WRITE_FAIL:
    reterr = kPglRetWriteFail;
    break;
  SimulateDataset_ret_MALFORMED_INPUT_WW:
    WordWrapB(0);
    logerrputsb();
  SimulateDataset_ret_MALFORMED_INPUT:
    reterr = kPglRetMalformedInput;
    break;
  SimulateDataset_ret_INCONSISTENT_INPUT:
    reterr = kPglRetInconsistentInput;
    break;
  }
 SimulateDataset_ret_1:
  CleanupSpgw(&spgw, &reterr);
  CleanupTextStream2(simulate_fname, &txs, &reterr);
  fclose_cond(simfreq_file);
  fclose_cond(outfile);
  BigstackReset(bigstack_mark);
  return reterr;
}

#ifdef __cplusplus
}
#endif
