// This library is part of PLINK 2.00, copyright (C) 2005-2020 Shaun Purcell,
// Christopher Chang.
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
    if (unlikely(
            SetThreadCt0(thread_ct - 1, &tg) ||
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
    if (unlikely(
            SetThreadCt0(thread_ct - 1, &tg) ||
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

#ifdef __cplusplus
}
#endif
